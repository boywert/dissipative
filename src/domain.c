/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/domain.c
 * \date        MM/YYYY
 * \author     
 * \brief        
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"
#include "domain.h"
#include "voronoi.h"

/*! \file domain.c
 *  \brief code for domain decomposition
 *
 *  This file contains the code for the domain decomposition of the
 *  simulation volume.  The domains are constructed from disjoint subsets
 *  of the leaves of a fiducial top-level tree that covers the full
 *  simulation volume. Domain boundaries hence run along tree-node
 *  divisions of a fiducial global BH tree. As a result of this method, the
 *  tree force are in principle strictly independent of the way the domains
 *  are cut. The domain decomposition can be carried out for an arbitrary
 *  number of CPUs. Individual domains are not cubical, but spatially
 *  coherent since the leaves are traversed in a Peano-Hilbert order and
 *  individual domains form segments along this order.  This also ensures
 *  that each domain has a small surface to volume ratio, which minimizes
 *  communication.
 */



/*! This is the main routine for the domain decomposition.  It acts as a
 *  driver routine that allocates various temporary buffers, maps the
 *  particles back onto the periodic box if needed, and then does the
 *  domain decomposition, and a final Peano-Hilbert order of all particles
 *  as a tuning measure.
 */
void domain_Decomposition(void)
{
  TIMER_START(CPU_DOMAIN);

  double t0 = second();

  mpi_printf("DOMAIN: Begin domain decomposition (sync-point %d).\n", All.NumCurrentTiStep);

  domain_prechecks();
  domain_prepare_voronoi_dynamic_update();
  /* map the particles back onto the box */
  do_box_wrapping();
  domain_init_sum_cost();
  domain_allocate();
  domain_allocate_lists();
  topNodes = (struct local_topnode_data *) mymalloc_movable(&topNodes, "topNodes", (MaxTopNodes * sizeof(struct local_topnode_data)));
  /* find total cost factors */
  domain_find_total_cost();
  /* determine global dimensions of domain grid */
  domain_findExtent();

  /* determine top-level tree */
  domain_determineTopTree();

  /* find the split of the top-level tree */
  domain_combine_topleaves_to_domains(All.MultipleDomains * NTask, NTopleaves);

  /* combine on each MPI task several of the domains (namely the number All.MultipleDomains) */
  domain_combine_multipledomains();

  /* permutate the task assignment such that the smallest number of particles needs to be moved */
  domain_optimize_domain_to_task_mapping();

  double ta = second();
  /* in case we retain the neighbor connectivity, do some preparatory flagging */
  domain_voronoi_dynamic_flag_particles();
  /* eliminate cells that might have been eliminated or were turned into stars */
  domain_rearrange_particle_sequence();
  /* determine for each cpu how many particles have to be shifted to other cpus */
  domain_countToGo();
  double tb = second();
  mpi_printf("DOMAIN: particle rearrangement work took %g sec\n", timediff(ta, tb));

  /* finally, carry out the actual particle exchange */
  domain_exchange();

  /* copy what we need for the topnodes */
  domain_preserve_relevant_topnode_data();
  myfree(topNodes);
  domain_free_lists();
  TimeOfLastDomainConstruction = All.Time;

  double t1 = second();
  mpi_printf("DOMAIN: domain decomposition done. (took in total %g sec)\n", timediff(t0, t1));

  TIMER_STOP(CPU_DOMAIN);
  TIMER_START(CPU_PEANO);

  peano_hilbert_order();
  myfree(Key);

  TIMER_STOPSTART(CPU_PEANO, CPU_DOMAIN);

  myfree(DomainListOfLocalTopleaves);

#ifndef AMR
#ifdef ONEDIMS
  voronoi_1D_order();
#endif
#endif

  TopNodes = (struct topnode_data *) myrealloc_movable(TopNodes, NTopnodes * sizeof(struct topnode_data));
  DomainTask = (int *) myrealloc_movable(DomainTask, NTopleaves * sizeof(int));

  domain_voronoi_dynamic_update_execute();

  DomainListOfLocalTopleaves = (int *) mymalloc_movable(&DomainListOfLocalTopleaves, "DomainListOfLocalTopleaves", (NTopleaves * sizeof(int)));

  memset(DomainNLocalTopleave, 0, NTask * sizeof(int));

  for(int i = 0; i < NTopleaves; i++)
    DomainNLocalTopleave[DomainTask[i]]++;

  DomainFirstLocTopleave[0] = 0;
  for(int i = 1; i < NTask; i++)
    DomainFirstLocTopleave[i] = DomainFirstLocTopleave[i - 1] + DomainNLocalTopleave[i - 1];

  memset(DomainNLocalTopleave, 0, NTask * sizeof(int));

  for(int i = 0; i < NTopleaves; i++)
    {
      int task = DomainTask[i];
      int off = DomainFirstLocTopleave[task] + DomainNLocalTopleave[task]++;
      DomainListOfLocalTopleaves[off] = i;
    }

  reconstruct_timebins();

  for(int i = 0; i < GRAVCOSTLEVELS; i++)
    All.LevelHasBeenMeasured[i] = 0;

  domain_post_checks();

  domain_report_balance();

  TIMER_STOP(CPU_DOMAIN);
}



void domain_prepare_voronoi_dynamic_update(void)
{
#if defined(VORONOI_DYNAMIC_UPDATE) || defined(AMR_CONNECTIONS)
  /* prepare storage for translation table */
  N_trans = NumGas;             /* length of translation table */
  trans_table = mymalloc_movable(&trans_table, "trans_table", N_trans * sizeof(struct trans_data));
  MPI_Allreduce(&Nvc, &Largest_Nvc, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif
}

void domain_voronoi_dynamic_flag_particles(void)
{
#if defined(VORONOI_DYNAMIC_UPDATE) || defined(AMR_CONNECTIONS)
  /* flag the particles that need to be exported */
  for(int i = 0; i < NumPart; i++)
    {
      int no = 0;

      while(topNodes[no].Daughter >= 0)
        no = topNodes[no].Daughter + (Key[i] - topNodes[no].StartKey) / (topNodes[no].Size >> 3);

      no = topNodes[no].Leaf;

      int task = DomainTask[no];
      domain_mark_in_trans_table(i, task);
    }
#endif
}

void domain_voronoi_dynamic_update_execute(void)
{
#if defined(VORONOI_DYNAMIC_UPDATE) || defined(AMR_CONNECTIONS)
  CPU_Step[CPU_DOMAIN] += measure_time();
  if(Largest_Nvc > 0)
    domain_exchange_and_update_DC();

  myfree_movable(trans_table);

  CPU_Step[CPU_MESH_DYNAMIC] += measure_time();
#endif
}

void domain_preserve_relevant_topnode_data(void)
{
#pragma omp parallel for
  for(int i = 0; i < NTopnodes; i++)
    {
      TopNodes[i].StartKey = topNodes[i].StartKey;
      TopNodes[i].Size = topNodes[i].Size;
      TopNodes[i].Daughter = topNodes[i].Daughter;
      TopNodes[i].Leaf = topNodes[i].Leaf;

      int bits = my_ffsll(TopNodes[i].Size);
      int blocks = (bits - 1) / 3 - 1;

      for(int j = 0; j < 8; j++)
        {
          peano1D xb, yb, zb;
          peano_hilbert_key_inverse(TopNodes[i].StartKey + j * (TopNodes[i].Size >> 3), BITS_PER_DIMENSION, &xb, &yb, &zb);
          xb >>= blocks;
          yb >>= blocks;
          zb >>= blocks;
          int idx = (xb & 1) | ((yb & 1) << 1) | ((zb & 1) << 2);
          if(idx < 0 || idx > 7)
            terminate("j=%d  idx=%d", j, idx);

          TopNodes[i].MortonToPeanoSubnode[idx] = j;
        }
    }
}


void domain_find_total_cost(void)
{
  if(All.MultipleDomains < 1 || All.MultipleDomains > 512)
    terminate("All.MultipleDomains < 1 || All.MultipleDomains > 512");

  gravcost = sphcost = 0;
  double partcount = 0;
  double sphpartcount = 0;

  for(int i = 0; i < NumPart; i++)
    {
#ifdef ADDBACKGROUNDGRID
      if(P[i].Type != 0)
        continue;
#endif
      partcount += 1.0;

      gravcost += domain_grav_tot_costfactor(i);

      double hydrocost = domain_hydro_tot_costfactor(i);
      sphcost += hydrocost;

      if(hydrocost > 0)
        sphpartcount += 1.0;
    }

  double loc[4] = {gravcost, sphcost, partcount, sphpartcount}, sum[4];

  MPI_Allreduce(loc, sum, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  totgravcost = sum[0];
  totsphcost = sum[1];
  totpartcount = sum[2];
  double totsphpartcount = sum[3];

  if(totsphcost > 0 && totgravcost > 0 && totsphpartcount > (All.TopNodeFactor * All.MultipleDomains * NTask))
    {
      /* in this case we give equal weight to gravitational work-load, hydro work load, and particle load.
       */
      normsum_work = 0.333333;
      normsum_load = 0.333333;
      normsum_worksph = 0.333333;
      fac_work = normsum_work / totgravcost;
      fac_load = normsum_load / totpartcount;
      fac_worksph = normsum_worksph / totsphcost;
    }
  else if(totgravcost > 0)
    {
      /* in this case we give equal weight to gravitational work-load and particle load.
       * The final pieces should have at most imbalance 2.0 in either of the two
       */
      normsum_work = 0.5;
      normsum_load = 0.5;
      normsum_worksph = 0;
      fac_work = normsum_work / totgravcost;
      fac_load = normsum_load / totpartcount;
      fac_worksph = 0.0;
    }
  else if(totsphcost > 0)
    {
      /* here we only appear to do hydrodynamics. We hence give equal weight to SPH cost and
       * particle load.
       */
      normsum_work = 0;
      normsum_load = 0.5;
      normsum_worksph = 0.5;
      fac_work = 0.0;
      fac_load = normsum_load / totpartcount;
      fac_worksph = normsum_worksph / totsphcost;
    }
  else
    terminate("strange: totsphcost=%g  totgravcost=%g\n", totsphcost, totgravcost);
}




peano1D domain_double_to_int(double d)
{
  union
  {
    double d;
    unsigned long long ull;
  } u;
  u.d = d;
  return (peano1D) ((u.ull & 0xFFFFFFFFFFFFFllu) >> (52 - BITS_PER_DIMENSION));
}



/*! This function allocates all the stuff that will be required for the tree-construction/walk later on */
void domain_allocate(void)
{
  MaxTopNodes = (int) (All.TopNodeAllocFactor * All.MaxPart + 1);

  if(DomainStartList)
    terminate("domain storage already allocated");

  DomainStartList = (int *) mymalloc_movable(&DomainStartList, "DomainStartList", (NTask * All.MultipleDomains * sizeof(int)));
  DomainEndList = (int *) mymalloc_movable(&DomainEndList, "DomainEndList", (NTask * All.MultipleDomains * sizeof(int)));
  DomainFirstLocTopleave = (int *) mymalloc_movable(&DomainFirstLocTopleave, "DomainFirstLocTopleave", NTask * sizeof(int));
  DomainNLocalTopleave = (int *) mymalloc_movable(&DomainNLocalTopleave, "DomainNLocalTopleave", NTask * sizeof(int));
  TopNodes = (struct topnode_data *) mymalloc_movable(&TopNodes, "TopNodes", (MaxTopNodes * sizeof(struct topnode_data)));
  DomainTask = (int *) mymalloc_movable(&DomainTask, "DomainTask", (MaxTopNodes * sizeof(int)));
  DomainListOfLocalTopleaves = (int *) mymalloc_movable(&DomainListOfLocalTopleaves, "DomainListOfLocalTopleaves", (MaxTopNodes * sizeof(int)));
}



void domain_free(void)
{
  if(!DomainStartList)
    terminate("domain storage not allocated");

  myfree_movable(DomainListOfLocalTopleaves);
  myfree_movable(DomainTask);
  myfree_movable(TopNodes);
  myfree_movable(DomainNLocalTopleave);
  myfree_movable(DomainFirstLocTopleave);
  myfree_movable(DomainEndList);
  myfree_movable(DomainStartList);

  DomainTask = NULL;
  TopNodes = NULL;
  DomainNLocalTopleave = NULL;
  DomainFirstLocTopleave = NULL;
  DomainEndList = NULL;
  DomainStartList = NULL;
}

void domain_printf(char *buf)
{
  if(RestartFlag <= 2)
    {
      fprintf(FdDomain, "%s", buf);
    }
}


void domain_report_balance(void)
{
  /* get total particle counts */
  long long loc_count[2 * TIMEBINS], glob_count[2 * TIMEBINS];

  for(int i = 0; i < TIMEBINS; i++)
    {
      loc_count[i] = TimeBinsGravity.TimeBinCount[i];
      loc_count[TIMEBINS + i] = TimeBinsHydro.TimeBinCount[i];
    }

  MPI_Reduce(loc_count, glob_count, 2 * TIMEBINS, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  double loc_max_data[2 * TIMEBINS + 3], glob_max_data[2 * TIMEBINS + 3];
  loc_max_data[2 * TIMEBINS + 0] = NumPart;
  loc_max_data[2 * TIMEBINS + 1] = NumGas;
  loc_max_data[2 * TIMEBINS + 2] = NumPart - NumGas;

  double glob_sum_data[2 * TIMEBINS];

  double *loc_HydroCost = &loc_max_data[0];
  double *loc_GravCost = &loc_max_data[TIMEBINS];
  double *max_HydroCost = &glob_max_data[0];
  double *max_GravCost = &glob_max_data[TIMEBINS];
  double *glob_HydroCost = &glob_sum_data[0];
  double *glob_GravCost = &glob_sum_data[TIMEBINS];


  for(int i = 0; i < TIMEBINS; i++)
    {
      loc_GravCost[i] = 0;
      loc_HydroCost[i] = 0;
    }

#ifdef SELFGRAVITY
  for(int i = 0; i < NumPart; i++)
    {
      for(int bin = All.LowestOccupiedTimeBin; bin <= All.HighestOccupiedTimeBin; bin++)
        {
 #ifdef HIERARCHICAL_GRAVITY
          if(bin >= P[i].TimeBinGrav)
#endif
            {
              if(domain_bintolevel[bin] >= 0)
                loc_GravCost[bin] += MIN_FLOAT_NUMBER + domain_grav_weight[bin] * P[i].GravCost[domain_bintolevel[bin]];
              else
                {
                  if(domain_refbin[bin] >= 0)
                    loc_GravCost[bin] += MIN_FLOAT_NUMBER + domain_grav_weight[bin] * P[i].GravCost[domain_bintolevel[domain_refbin[bin]]];
                  else
                    loc_GravCost[bin] += 1.0;
                }
            }
        }
    }
#endif

  for(int i = 0; i < NumPart; i++)
    if(P[i].Type == 0)
      loc_HydroCost[P[i].TimeBinHydro] += 1.0;

  /* now determine the cumulative cost for the hydrodynamics */
  for(int i = 1; i <= All.HighestOccupiedTimeBin; i++)
    loc_HydroCost[i] += loc_HydroCost[i - 1];


  MPI_Reduce(loc_max_data, glob_sum_data, 2 * TIMEBINS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(loc_max_data, glob_max_data, 2 * TIMEBINS + 3, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      double max_tot = glob_max_data[2*TIMEBINS + 0];
      double max_sph = glob_max_data[2*TIMEBINS + 1];
      double max_dm  = glob_max_data[2*TIMEBINS + 2];

      long long *tot_count = &glob_count[0];
      long long *tot_count_sph = &glob_count[TIMEBINS];

      long long tot_cumulative[TIMEBINS];
      tot_cumulative[0] = tot_count[0];

      for(int i = 1; i < TIMEBINS; i++)
        tot_cumulative[i] = tot_count[i] + tot_cumulative[i - 1];

      double tot_gravcost = 0, max_gravcost = 0, tot_hydrocost = 0, max_hydrocost = 0;

      All.TotGravCost = 0;

      for(int i = 0; i < TIMEBINS; i++)
        {
          All.TotGravCost += domain_to_be_balanced[i] * glob_GravCost[i] / NTask;

          tot_gravcost += domain_to_be_balanced[i] * glob_GravCost[i] / NTask;
          max_gravcost += domain_to_be_balanced[i] * max_GravCost[i];

          tot_hydrocost += domain_to_be_balanced[i] * glob_HydroCost[i] / NTask;
          max_hydrocost += domain_to_be_balanced[i] * max_HydroCost[i];
        }

      double bal_grav_bin[TIMEBINS], bal_grav_bin_rel[TIMEBINS];
      double bal_hydro_bin[TIMEBINS], bal_hydro_bin_rel[TIMEBINS];

      for(int i = 0; i < TIMEBINS; i++)
        {
          if(tot_count[i] > 0)
            {
              bal_grav_bin[i] = max_GravCost[i] / (glob_GravCost[i] / NTask + 1.0e-60);
              bal_grav_bin_rel[i] = (tot_gravcost + domain_to_be_balanced[i] * (max_GravCost[i] - glob_GravCost[i]
                  / NTask)) / (tot_gravcost + 1.0e-60);
            }
          else
            {
              bal_grav_bin[i] = 0.0;
              bal_grav_bin_rel[i] = 0.0;
            }

          if(tot_count_sph[i] > 0)
            {
              bal_hydro_bin[i] = max_HydroCost[i] / (glob_HydroCost[i] / NTask + 1.0e-60);
              bal_hydro_bin_rel[i] = (tot_hydrocost + domain_to_be_balanced[i] * (max_HydroCost[i] - glob_HydroCost[i]
                  / NTask)) / (tot_hydrocost + 1.0e-60);
            }
          else
            {
              bal_hydro_bin[i] = 0.0;
              bal_hydro_bin_rel[i] = 0.0;
            }
        }

      char buf[1000];
      sprintf(buf, "\nDOMAIN BALANCE, Sync-Point %d, Time: %g\n", All.NumCurrentTiStep, All.Time);
      domain_printf(buf);
      sprintf(buf, "Timebins:       Gravity       Hydro  cumulative      grav-balance       hydro-balance\n");
      domain_printf(buf);

      long long tot = 0, tot_sph = 0;

      for(int i = TIMEBINS - 1; i >= 0; i--)
        {
#if (defined(SELFGRAVITY) || defined(EXTERNALGRAVITY) || defined(EXACT_GRAVITY_FOR_PARTICLE_TYPE)) && !defined(MESHRELAX)
          if(tot_count_sph[i] > 0 || tot_count[i] > 0)
#else
          if(tot_count[i] > 0)
            tot += tot_count[i];

          if(tot_count_sph[i] > 0)
#endif
            {
              char buf[1000];
              sprintf(buf, "%c%cbin=%2d     %10llu  %10llu  %10llu  %c %6.3f |%6.3f  %c   %6.3f |%6.3f\n",
                      i == All.HighestActiveTimeBin ? '>' : ' ',
                      i >= All.SmallestTimeBinWithDomainDecomposition ? '|' : ' ',
                      i, tot_count[i], tot_count_sph[i], tot_cumulative[i],
                      domain_bintolevel[i] >= 0 ? 'm' : ' ', bal_grav_bin[i], bal_grav_bin_rel[i], domain_to_be_balanced[i] > 0 ? '*' : ' ', bal_hydro_bin[i], bal_hydro_bin_rel[i]);
              domain_printf(buf);

              tot += tot_count[i];
              tot_sph += tot_count_sph[i];
            }
        }

      sprintf(buf, "-------------------------------------------------------------------------------------\n");
      domain_printf(buf);
      sprintf(buf, "BALANCE,  LOAD:  %6.3f      %6.3f      %6.3f  WORK:     %6.3f              %6.3f\n",
              max_dm / (tot - tot_sph + 1.0e-60) * NTask,
              max_sph / (tot_sph + 1.0e-60) * NTask, max_tot / (tot + 1.0e-60) * NTask, max_gravcost / (tot_gravcost + 1.0e-60), max_hydrocost / (tot_hydrocost + 1.0e-60));
      domain_printf(buf);
      sprintf(buf, "-------------------------------------------------------------------------------------\n");
      domain_printf(buf);
      sprintf(buf, "\n");
      domain_printf(buf);
      myflush(FdDomain);
    }

  /* the following needs to be known by all the tasks */
  MPI_Bcast(&All.TotGravCost, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}
