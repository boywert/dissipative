/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/blackhole/blackhole_neighbors.c
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
#include <math.h>
#include <gsl/gsl_math.h>

#include "../allvars.h"
#include "../proto.h"

#ifdef BLACK_HOLES



static int blackhole_tree_ngb_evaluate(int i, int mode, int threadid);



#if defined(MEASURE_POTMIN_AROUND_BH) || defined(BH_BONDI_CAPTURE)

static struct bhresultsimported_data
{
  int index;
#ifdef MEASURE_POTMIN_AROUND_BH
  MyDouble BH_MinPotPos[3];
  MyFloat BH_MinPot;
  MyFloat BH_MinPot_ActiveM;
  MyFloat BH_MinPot_TotalM;
#endif
#ifdef BH_FRICTION
  MyFloat BH_RhoTot;
  MyDouble BH_MinPotPos_Extended[3];
  MyFloat BH_MinPot_Extended;
#endif
#ifdef BH_BONDI_CAPTURE
  MyFloat BH_CaptureMass;
#endif
}
 *Tree_BhngbResultsImported;

#endif



/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble Pos[3];
  MyDouble Mass;
  MyFloat Vel[3];
  MyFloat BH_Hsml;
  MyFloat BH_U;
  MyIDType ID;
#ifdef MASSIVE_SEEDS_MERGER
  MyFloat HostHaloMass;
#endif

  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;



/* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in * in, int i, int firstnode)
{
  int k;

  if(i < NumPart)
    {
      for(k = 0; k < 3; k++)
        {
          in->Pos[k] = P[i].Pos[k];
          in->Vel[k] = P[i].Vel[k];
        }
      in->Mass = P[i].Mass;
      in->BH_Hsml = BPP(i).BH_Hsml;
      in->BH_U = BPP(i).BH_U;
      in->ID = P[i].ID;
#ifdef MASSIVE_SEEDS_MERGER
      in->HostHaloMass = BPP(i).HostHaloMass;
#endif
    }
  else
    {
      i -= Tree_ImportedNodeOffset;

      for(k = 0; k < 3; k++)
        {
          in->Pos[k] = Tree_Points[i].Pos[k];
          in->Vel[k] = TBPP(i).Vel[k];
        }
      in->Mass = Tree_Points[i].Mass;
      in->BH_Hsml = TBPP(i).BH_Hsml;
      in->BH_U = TBPP(i).BH_U;
      in->ID = TBPP(i).ID;
#ifdef MASSIVE_SEEDS_MERGER
      in->HostHaloMass = TBPP(i).HostHaloMass;
#endif
    }

  in->Firstnode = firstnode;
}




/* local data structure that holds results acquired on remote processors */
typedef struct
{
#ifdef MEASURE_POTMIN_AROUND_BH
  MyDouble BH_MinPotPos[3];
  MyFloat BH_MinPot;
  MyFloat BH_MinPot_ActiveM;
  MyFloat BH_MinPot_TotalM;
#endif
#ifdef BH_FRICTION
  MyFloat BH_RhoTot;
  MyDouble BH_MinPotPos_Extended[3];
  MyFloat BH_MinPot_Extended;
#endif
#ifdef BH_BONDI_CAPTURE
  MyFloat BH_CaptureMass;
#endif
} data_out;

static data_out *DataResult, *DataOut;




 /* routine to store or combine result data */
static void out2particle(data_out * out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES)      /* initial store */
    {
      if(i < NumPart)
        {
#ifdef MEASURE_POTMIN_AROUND_BH
          int k;
          BPP(i).BH_MinPot = out->BH_MinPot;
          for(k = 0; k < 3; k++)
            BPP(i).BH_MinPotPos[k] = out->BH_MinPotPos[k];
          BPP(i).BH_MinPot_ActiveM = out->BH_MinPot_ActiveM;
          BPP(i).BH_MinPot_TotalM = out->BH_MinPot_TotalM;
#endif
#ifdef BH_FRICTION
          BPP(i).BH_RhoTot = out->BH_RhoTot;
          BPP(i).BH_MinPot_Extended = out->BH_MinPot_Extended;
          for(k = 0; k < 3; k++)
            BPP(i).BH_MinPotPos_Extended[k] = out->BH_MinPotPos_Extended[k];
#endif
#ifdef BH_BONDI_CAPTURE
          BPP(i).BH_CaptureMass = out->BH_CaptureMass;
#endif
        }
      else
        {
#ifdef MEASURE_POTMIN_AROUND_BH
          int k, idx = Tree_ResultIndexList[i - Tree_ImportedNodeOffset];

          Tree_BhngbResultsImported[idx].BH_MinPot = out->BH_MinPot;
          for(k = 0; k < 3; k++)
            Tree_BhngbResultsImported[idx].BH_MinPotPos[k] = out->BH_MinPotPos[k];
          Tree_BhngbResultsImported[idx].BH_MinPot_ActiveM = out->BH_MinPot_ActiveM;
          Tree_BhngbResultsImported[idx].BH_MinPot_TotalM = out->BH_MinPot_TotalM;
#endif
#ifdef BH_FRICTION
          Tree_BhngbResultsImported[idx].BH_RhoTot = out->BH_RhoTot;
          Tree_BhngbResultsImported[idx].BH_MinPot_Extended = out->BH_MinPot_Extended;
          for(k = 0; k < 3; k++)
            Tree_BhngbResultsImported[idx].BH_MinPotPos_Extended[k] = out->BH_MinPotPos_Extended[k];
#endif
#ifdef BH_BONDI_CAPTURE
          Tree_BhngbResultsImported[idx].BH_CaptureMass = out->BH_CaptureMass;
#endif
        }
    }
  else                          /* combine */
    {
      if(i < NumPart)
        {
#ifdef MEASURE_POTMIN_AROUND_BH
          int k;

          if(BPP(i).BH_MinPot > out->BH_MinPot)
            {
              BPP(i).BH_MinPot = out->BH_MinPot;
              for(k = 0; k < 3; k++)
                BPP(i).BH_MinPotPos[k] = out->BH_MinPotPos[k];
            }
          BPP(i).BH_MinPot_ActiveM += out->BH_MinPot_ActiveM;
          BPP(i).BH_MinPot_TotalM += out->BH_MinPot_TotalM;
#endif
#ifdef BH_FRICTION
          BPP(i).BH_RhoTot += out->BH_RhoTot;
          BPP(i).BH_MinPot_Extended = out->BH_MinPot_Extended;
          for(k = 0; k < 3; k++)
            BPP(i).BH_MinPotPos_Extended[k] = out->BH_MinPotPos_Extended[k];
#endif
#ifdef BH_BONDI_CAPTURE
          BPP(i).BH_CaptureMass += out->BH_CaptureMass;
#endif
        }
      else
        {
#ifdef MEASURE_POTMIN_AROUND_BH
          int k, idx = Tree_ResultIndexList[i - Tree_ImportedNodeOffset];

          if(Tree_BhngbResultsImported[idx].BH_MinPot > out->BH_MinPot)
            {
              Tree_BhngbResultsImported[idx].BH_MinPot = out->BH_MinPot;
              for(k = 0; k < 3; k++)
                Tree_BhngbResultsImported[idx].BH_MinPotPos[k] = out->BH_MinPotPos[k];
            }
          Tree_BhngbResultsImported[idx].BH_MinPot_ActiveM += out->BH_MinPot_ActiveM;
          Tree_BhngbResultsImported[idx].BH_MinPot_TotalM += out->BH_MinPot_TotalM;
#endif
#ifdef BH_FRICTION
          Tree_BhngbResultsImported[idx].BH_RhoTot += out->BH_RhoTot;
          Tree_BhngbResultsImported[idx].BH_MinPot_Extended = out->BH_MinPot_Extended;
          for(k = 0; k < 3; k++)
            Tree_BhngbResultsImported[idx].BH_MinPotPos_Extended[k] = out->BH_MinPotPos_Extended[k];
#endif
#ifdef BH_BONDI_CAPTURE
          Tree_BhngbResultsImported[idx].BH_CaptureMass += out->BH_CaptureMass;
#endif
        }
    }
}


#include "../generic_comm_helpers2.h"

static void kernel_local(void)
{
  int i;
#ifdef GENERIC_ASYNC
  int flag = 0;
#endif
#pragma omp parallel private(i)
  {
    int j, threadid = get_thread_num();
#ifdef GENERIC_ASYNC
    int count = 0;
#endif

    for(j = 0; j < NTask; j++)
      Thread[threadid].Exportflag[j] = -1;

    while(1)
      {
        if(Thread[threadid].ExportSpace < MinSpace)
          break;

#ifdef GENERIC_ASYNC
        if(threadid == 0)
          {
            if((count & POLLINGINTERVAL) == 0)
              if(generic_polling_primary(count, Nforces))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        i = NextParticle++;

        if(i >= Nforces)
          break;

        int idx = TargetList[i];

        blackhole_tree_ngb_evaluate(idx, MODE_LOCAL_PARTICLES, threadid);
      }
  }
}


static void kernel_imported(void)
{
  /* now do the particles that were sent to us */
  int i, cnt = 0;
#pragma omp parallel private(i)
  {
    int threadid = get_thread_num();
#ifdef GENERIC_ASYNC
    int count = 0;
#endif

    while(1)
      {
#pragma omp atomic capture
        i = cnt++;

        if(i >= Nimport)
          break;

#ifdef GENERIC_ASYNC
        if(threadid == 0)
          {
            if((count & POLLINGINTERVAL) == 0)
              generic_polling_secondary();
          }

        count++;
#endif

        blackhole_tree_ngb_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

/* This routines uses the gravitational tree to search in the BH_Hsml neighborhood around each
 * active BH for other BHs (which become candidates for merging), and if desired, it also determines
 * the location of the minimum gravitational potential in this environment.
 */
void blackhole_find_neighboring_holes_and_potmin(void)
{
  int idx, i, ncount;

  TIMER_START(CPU_BH_NGB);

  mpi_printf("BLACK_HOLES: Begin finding BH neighbours\n");

#ifdef MEASURE_POTMIN_AROUND_BH
  forcetree_update_exported_potential_values();
#ifdef BH_FRICTION
  blackhole_friction_store_previous_minimum();
#endif
#endif

  generic_set_MaxNexport();


  /* Create list of targets. We do this here to simplify the treatment of the two possible locations of source points */
  TargetList = mymalloc("TargetList", (NumPart + Tree_NumPartImported) * sizeof(int));
  Tree_ResultIndexList = mymalloc("Tree_ResultIndexList", Tree_NumPartImported * sizeof(int));

  Nforces = 0;
  for(idx = 0; idx < TimeBinsBHAccretion.NActiveParticles; idx++)
    {
      i = TimeBinsBHAccretion.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Type == 5 && Tree_Task_list[i] == ThisTask)
        TargetList[Nforces++] = i;
#ifdef BH_BONDI_CAPTURE
      if(P[i].Type == 5 && Tree_Task_list[i] == ThisTask)
        BPP(i).BH_CaptureMass = 0;
#endif
    }

  for(i = 0, ncount = 0; i < Tree_NumPartImported; i++)
#ifndef HIERARCHICAL_GRAVITY
    if(Tree_Points[i].ActiveFlag)
#endif
    if(Tree_Points[i].Type == 5)
      {
        Tree_ResultIndexList[i] = ncount++;
        TargetList[Nforces++] = i + Tree_ImportedNodeOffset;
      }

#if defined(MEASURE_POTMIN_AROUND_BH) || defined(BH_BONDI_CAPTURE)
  Tree_BhngbResultsImported = mymalloc("Tree_BhngbResultsImported", ncount * sizeof(struct bhresultsimported_data));
#endif


  generic_comm_pattern(Nforces, kernel_local, kernel_imported);


#if defined(MEASURE_POTMIN_AROUND_BH) || defined(BH_BONDI_CAPTURE)
  int ngrp, recvTask, nexport, nimport, k, j, n;

  /* now communicate the results in Tree_BhngbResultsImported */
  for(j = 0; j < NTask; j++)
    Recv_count[j] = 0;

  for(i = 0, n = 0, k = 0; i < NTask; i++)
    for(j = 0; j < Force_Recv_count[i]; j++, n++)
      {
#ifndef HIERARCHICAL_GRAVITY
        if(Tree_Points[n].ActiveFlag)
#endif
        if(Tree_Points[n].Type == 5)
          {
            Tree_BhngbResultsImported[k].index = Tree_Points[n].index;
            Recv_count[i]++;
            k++;
          }
      }

  MPI_Alltoall(Recv_count, 1, MPI_INT, Send_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nexport = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nexport += Send_count[j];
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  struct bhresultsimported_data *tmp_results = mymalloc("tmp_results", nexport * sizeof(struct bhresultsimported_data));
  memset(tmp_results, -1, nexport * sizeof(struct bhresultsimported_data));

  /* exchange  data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          MPI_Sendrecv(&Tree_BhngbResultsImported[Recv_offset[recvTask]],
                       Recv_count[recvTask] * sizeof(struct bhresultsimported_data), MPI_BYTE, recvTask, TAG_FOF_A,
                       &tmp_results[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct bhresultsimported_data), MPI_BYTE, recvTask, TAG_FOF_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

  for(i = 0; i < nexport; i++)
    {
      int target = tmp_results[i].index;
#ifdef MEASURE_POTMIN_AROUND_BH
      BPP(target).BH_MinPot = tmp_results[i].BH_MinPot;
      for(k = 0; k < 3; k++)
        BPP(target).BH_MinPotPos[k] = tmp_results[i].BH_MinPotPos[k];

      BPP(target).BH_MinPot_ActiveM = tmp_results[i].BH_MinPot_ActiveM;
      BPP(target).BH_MinPot_TotalM = tmp_results[i].BH_MinPot_TotalM;
#endif
#ifdef BH_FRICTION
      BPP(target).BH_MinPot_Extended = tmp_results[i].BH_MinPot_Extended;
      for(k = 0; k < 3; k++)
        BPP(target).BH_MinPotPos_Extended[k] = tmp_results[i].BH_MinPotPos_Extended[k];

      BPP(target).BH_RhoTot += tmp_results[i].BH_RhoTot;
#endif
#ifdef BH_BONDI_CAPTURE
      BPP(target).BH_CaptureMass = tmp_results[i].BH_CaptureMass;
#endif
    }

  myfree(tmp_results);

  myfree(Tree_BhngbResultsImported);
#endif
  myfree(Tree_ResultIndexList);
  myfree(TargetList);

#ifdef BH_FRICTION
  blackhole_friction_update_vel_pot_minimum();
#endif

  mpi_printf("BLACK_HOLES: done with BH neighbours\n");
#ifdef BLACKHOLE_POTMIN_DIAGNOSTIC
  blackhole_potmin_diagnostic();
#endif

  TIMER_STOP(CPU_BH_NGB);
}

#ifdef REPOSITION_ON_POTMIN
void blackhole_reposition(void)
{
  TIMER_START(CPU_BH_NGB);

  int idx, i, k;

  /* do the repositionening */
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Type == 5)
        {
          if(BPP(i).BH_MinPot_ActiveM >= BH_REPOSITION_POTMIN_TRUST_THRESHOLD * BPP(i).BH_MinPot_TotalM)
            {
              if((BPP(i).BH_MinPot < 0.5 * BHPOTVALUEINIT) && (BPP(i).BH_MinPot != 0))
                {
                  for(k = 0; k < 3; k++)
                    P[i].Pos[k] = BPP(i).BH_MinPotPos[k];
                }
            }
        }
    }

  TIMER_STOP(CPU_BH_NGB);
}
#endif


static int blackhole_tree_ngb_evaluate(int target, int mode, int threadid)
{
  int k, numnodes, *firstnode;
  int no;
  double h, h2;
#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif
  double dx, dy, dz, r2;
  MyIDType id;
  MyDouble *pos;
#ifdef MEASURE_POTMIN_AROUND_BH
  MyFloat minpotpos[3] = { 0, 0, 0 }, minpot = BHPOTVALUEINIT;
#endif
#ifdef BH_FRICTION
  MyFloat minpotpos_extended[3] = { 0, 0, 0 }, minpot_extended = BHPOTVALUEINIT;
  double wk, rhotot = 0;
#endif
#ifdef BH_BONDI_CAPTURE
  MyFloat mass_capture = 0;
#endif
#ifdef MASSIVE_SEEDS_MERGER
  double mass, mass_max, Rsearch;
#endif

  data_in local, *target_data;
  data_out out;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      target_data = &local;

      numnodes = 1;
      firstnode = NULL;
    }
  else
    {
      target_data = &DataGet[target];

      generic_get_numnodes(target, &numnodes, &firstnode);
    }

  pos = target_data->Pos;
  h = target_data->BH_Hsml;
  id = target_data->ID;

#if !defined(BH_DO_NOT_PREVENT_MERGERS) && !defined(MASSIVE_SEEDS_MERGER)
  double csnd = sqrt(GAMMA * GAMMA_MINUS1 * target_data->BH_U);
#endif

#ifdef MASSIVE_SEEDS_MERGER
  mass = target_data->Mass;
  Rsearch = target_data->HostHaloMass;
  Rsearch = pow(All.G * Rsearch / (100.0 * All.cf_H * All.cf_H), 1.0 / 3) / All.cf_atime;
  h = 0.05 * Rsearch;
#endif

  h2 = h * h;

#ifdef BH_FRICTION
  double hinv = 1.0 / h;
#ifndef  TWODIMS
  double hinv3 = hinv * hinv * hinv;
#else
  double hinv3 = hinv * hinv / boxSize_Z;
#endif
#endif

#ifdef BH_FRICTION
  double h2_extended = FAC_TWO_TO_TWO_THIRDS * h2;
  h *= sqrt(FAC_TWO_TO_TWO_THIRDS);
#endif

#ifdef MEASURE_POTMIN_AROUND_BH
  out.BH_MinPot_ActiveM = 0;
  out.BH_MinPot_TotalM = 0;
#endif

  for(k = 0; k < numnodes; k++)
    {
      if(mode == MODE_LOCAL_PARTICLES)
        {
          no = Tree_MaxPart;    /* root node */
        }
      else
        {
          no = firstnode[k];
          no = Nodes[no].u.d.nextnode;  /* open it */
        }

      while(no >= 0)
        {
          if(no < Tree_MaxPart) /* single particle */
            {
              dx = GRAVITY_NEAREST_X(Tree_Pos_list[3 * no + 0] - pos[0]);
              dy = GRAVITY_NEAREST_Y(Tree_Pos_list[3 * no + 1] - pos[1]);
              dz = GRAVITY_NEAREST_Z(Tree_Pos_list[3 * no + 2] - pos[2]);

              r2 = dx * dx + dy * dy + dz * dz;

              if(r2 < h2)
                {
                  /* minpot */
#ifdef MEASURE_POTMIN_AROUND_BH
                  out.BH_MinPot_TotalM += P[no].Mass;
                  if(TimeBinSynchronized[P[no].TimeBinGrav])    /* only use active particles who just had their potential updated */
                    {
                      out.BH_MinPot_ActiveM += P[no].Mass;
                      if(minpot > P[no].Potential)
                        {
                          minpot = P[no].Potential;

                          for(k = 0; k < 3; k++)
                            minpotpos[k] = P[no].Pos[k];
                        }
                    }
#endif
#ifdef BH_FRICTION
                  double r = sqrt(r2);
                  double u = r * hinv;

                  if(u < 0.5)
                    wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
                  else
                    wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

                  rhotot += P[no].Mass * wk;
#endif
#ifdef BH_BONDI_CAPTURE
                  if(P[no].Type == 0 || P[no].Type == 1 || P[no].Type == 4)     /* sum the total mass in gas, stars and dark matter */
                    mass_capture += P[no].Mass;
#endif
                  if(P[no].Type == 5)   /* we have a potential black hole merger */
                    {
                      if(id != P[no].ID)
                        {
#ifndef BH_DO_NOT_PREVENT_MERGERS
                          double vrel;
                          /* compute relative velocity of BHs */
#error "this doesn't work at present since we are not guaranteed to have access to velocities of particles in the tree if they are an imported TreePoint"
                          for(k = 0, vrel = 0; k < 3; k++)
                            vrel += (P[no].Vel[k] - vel[k]) * (P[no].Vel[k] - vel[k]);

                          vrel = sqrt(vrel) / All.cf_atime;

#ifdef MASSIVE_SEEDS_MERGER
                          mass_max = mass + P[no].Mass;

                          if(vrel < sqrt(All.G * mass_max / h))
#else
                          if(vrel < CSND_FRAC_BH_MERGE * csnd)
#endif
#endif
                            {
                              if(BPP(no).SwallowID < id && P[no].ID < id)
                                BPP(no).SwallowID = id;
                            }
                        }
                    }
                }
#ifdef BH_FRICTION
              else if(r2 < h2_extended)
                {
                  out.BH_MinPot_TotalM += P[no].Mass;
                  if(TimeBinSynchronized[P[no].TimeBinGrav])    /* only use active particles who just had their potential updated */
                    {
                      out.BH_MinPot_ActiveM += P[no].Mass;
                      if(minpot_extended > P[no].Potential)
                        {
                          minpot_extended = P[no].Potential;

                          for(k = 0; k < 3; k++)
                            minpotpos_extended[k] = P[no].Pos[k];
                        }
                    }
                }
#endif
              no = Nextnode[no];
            }
          else if(no < Tree_MaxPart + Tree_MaxNodes)    /* internal node */
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                {
                  if(no < Tree_FirstNonTopLevelNode)    /* we reached a top-level node again, which means that we are done with the branch */
                    break;
                }

              struct NODE *current = &Nodes[no];

              no = current->u.d.sibling;        /* in case the node can be discarded */

              double dist = h + 0.5 * current->len;
              dx = NGB_PERIODIC_LONG_X(current->center[0] - pos[0]);
              if(dx > dist)
                continue;
              dy = NGB_PERIODIC_LONG_Y(current->center[1] - pos[1]);
              if(dy > dist)
                continue;
              dz = NGB_PERIODIC_LONG_Z(current->center[2] - pos[2]);
              if(dz > dist)
                continue;
              /* now test against the minimal sphere enclosing everything */
              dist += FACT1 * current->len;
              if(dx * dx + dy * dy + dz * dz > dist * dist)
                continue;

              no = current->u.d.nextnode;       /* ok, we need to open the node */
            }
          else if(no >= Tree_ImportedNodeOffset)        /* point from imported nodelist */
            {
              int n = no - Tree_ImportedNodeOffset;

              dx = GRAVITY_NEAREST_X(Tree_Points[n].Pos[0] - pos[0]);
              dy = GRAVITY_NEAREST_Y(Tree_Points[n].Pos[1] - pos[1]);
              dz = GRAVITY_NEAREST_Z(Tree_Points[n].Pos[2] - pos[2]);

              r2 = dx * dx + dy * dy + dz * dz;

              if(r2 < h2)
                {
#ifdef MEASURE_POTMIN_AROUND_BH
                  out.BH_MinPot_TotalM += Tree_Points[n].Mass;

                  if(Tree_Points[n].Type & 16)  /* only use active particles who had their potential updated */
                    {
                      out.BH_MinPot_ActiveM += Tree_Points[n].Mass;

                      if(minpot > Tree_Points[n].Potential)
                        {
                          minpot = Tree_Points[n].Potential;

                          for(k = 0; k < 3; k++)
                            minpotpos[k] = Tree_Points[n].Pos[k];
                        }
                    }
#endif
#ifdef BH_FRICTION
                  double r = sqrt(r2);
                  double u = r * hinv;

                  if(u < 0.5)
                    wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
                  else
                    wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

                  rhotot += Tree_Points[n].Mass * wk;
#endif

#ifdef BH_BONDI_CAPTURE
                  /* sum the total mass in gas, stars and dark matter */
                  if(Tree_Points[n].Type == 0 || Tree_Points[n].Type == 1 || Tree_Points[n].Type == 4)
                    mass_capture += Tree_Points[n].Mass;
#endif
                  if(Tree_Points[n].Type == 5)   /* we have a potential black hole merger */
                    {
                      if(id != TBPP(n).ID)
                        {
#ifndef BH_DO_NOT_PREVENT_MERGERS
                          double vrel;
                          /* compute relative velocity of BHs */
                          for(k = 0, vrel = 0; k < 3; k++)
                            vrel += (TBPP(n).Vel[k] - vel[k]) * (TBPP(n).Vel[k] - vel[k]);

                          vrel = sqrt(vrel) / All.cf_atime;

#ifdef MASSIVE_SEEDS_MERGER
                          mass_max = mass + Tree_Points[n].Mass;

                          if(vrel < sqrt(All.G * mass_max / h))
#else
                          if(vrel < CSND_FRAC_BH_MERGE * csnd)
#endif
#endif
                            {
                              if(TBPP(n).SwallowID < id && TBPP(n).ID < id)
                                TBPP(n).SwallowID = id;
                            }
                        }
                    }
                }
#ifdef BH_FRICTION
              else if(r2 < h2_extended)
                {
                  out.BH_MinPot_TotalM += Tree_Points[n].Mass;

                  if(Tree_Points[n].Type & 16)  /* only use active particles who had their potential updated */
                    {
                      out.BH_MinPot_ActiveM += Tree_Points[n].Mass;

                      if(minpot_extended > Tree_Points[n].Potential)
                        {
                          minpot_extended = Tree_Points[n].Potential;

                          for(k = 0; k < 3; k++)
                            minpotpos_extended[k] = Tree_Points[n].Pos[k];
                        }
                    }
                }
#endif
              no = Nextnode[no - Tree_MaxNodes];
            }
          else                  /* pseudo particle */
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                terminate("mode == MODE_IMPORTED_PARTICLES");

              if(target >= 0)
                tree_treefind_export_node_threads(no, target, threadid);

              no = Nextnode[no - Tree_MaxNodes];
              continue;
            }
        }

    }


#ifdef MEASURE_POTMIN_AROUND_BH
  out.BH_MinPot = minpot;
  for(k = 0; k < 3; k++)
    out.BH_MinPotPos[k] = minpotpos[k];
#endif
#ifdef BH_FRICTION
  out.BH_RhoTot = rhotot;
  out.BH_MinPot_Extended = minpot_extended;
  for(k = 0; k < 3; k++)
    out.BH_MinPotPos_Extended[k] = minpotpos_extended[k];
#endif
#ifdef BH_BONDI_CAPTURE
  out.BH_CaptureMass = mass_capture;
#endif

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#ifdef BLACKHOLE_POTMIN_DIAGNOSTIC
void blackhole_potmin_diagnostic(void)
{
  int i;

  struct minpot_data
  {
    MyDouble pos[3];
    double minpot;
    double bh_pos[3];
    double bh_hsml;
    double bh_minpot;
    double bh_rhotot;
    double bh_minpotvel[3];
    double Pos[3], Vel[3];
  } my, *mylist;

  my.minpot = 1.0e30;
  my.bh_hsml = 0;
  my.bh_minpot = 0;

  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type == 5)
        {
          my.bh_pos[0] = BPP(i).BH_MinPotPos[0];
          my.bh_pos[1] = BPP(i).BH_MinPotPos[1];
          my.bh_pos[2] = BPP(i).BH_MinPotPos[2];
          my.bh_hsml = BPP(i).BH_Hsml;
          my.bh_minpot = BPP(i).BH_MinPot;
          my.bh_rhotot = BPP(i).BH_RhoTot;
          my.bh_minpotvel[0] = BPP(i).BH_MinPotVel[0];
          my.bh_minpotvel[1] = BPP(i).BH_MinPotVel[1];
          my.bh_minpotvel[2] = BPP(i).BH_MinPotVel[2];

          my.Pos[0] = P[i].Pos[0];
          my.Pos[1] = P[i].Pos[1];
          my.Pos[2] = P[i].Pos[2];
          my.Vel[0] = P[i].Vel[0];
          my.Vel[1] = P[i].Vel[1];
          my.Vel[2] = P[i].Vel[2];
        }

      if(my.minpot > P[i].Potential)
        {
          my.minpot = P[i].Potential;
#ifdef CELL_CENTER_GRAVITY
          if(P[i].Type == 0)
            {
              my.pos[0] = SphP[i].Center[0];
              my.pos[1] = SphP[i].Center[1];
              my.pos[2] = SphP[i].Center[2];
            }
          else
#endif
            {
              my.pos[0] = P[i].Pos[0];
              my.pos[1] = P[i].Pos[1];
              my.pos[2] = P[i].Pos[2];
            }
        }
    }


  mylist = mymalloc("mylist", NTask * sizeof(struct minpot_data));

  MPI_Allgather(&my, sizeof(struct minpot_data), MPI_BYTE, mylist, sizeof(struct minpot_data), MPI_BYTE, MPI_COMM_WORLD);

  int task = 0;

  for(i = 0; i < NTask; i++)
    {
      if(mylist[i].minpot < my.minpot)
        {
          my.minpot = mylist[i].minpot;
          my.pos[0] = mylist[i].pos[0];
          my.pos[1] = mylist[i].pos[1];
          my.pos[2] = mylist[i].pos[2];
          task = i;
        }
      if(mylist[i].bh_hsml > 0)
        {
          my.bh_hsml = mylist[i].bh_hsml;
          my.bh_pos[0] = mylist[i].bh_pos[0];
          my.bh_pos[1] = mylist[i].bh_pos[1];
          my.bh_pos[2] = mylist[i].bh_pos[2];
          my.bh_minpot = mylist[i].bh_minpot;
          my.bh_rhotot = mylist[i].bh_rhotot;
          my.bh_minpotvel[0] = mylist[i].bh_minpotvel[0];
          my.bh_minpotvel[1] = mylist[i].bh_minpotvel[1];
          my.bh_minpotvel[2] = mylist[i].bh_minpotvel[2];

          my.Pos[0] = mylist[i].Pos[0];
          my.Pos[1] = mylist[i].Pos[1];
          my.Pos[2] = mylist[i].Pos[2];
          my.Vel[0] = mylist[i].Vel[0];
          my.Vel[1] = mylist[i].Vel[1];
          my.Vel[2] = mylist[i].Vel[2];
        }
    }

  if(ThisTask == 0)
    {
      fprintf(FdBHDiag, "%g %d  %g %g %g   %12.7f %g  %g %g %g  %12.7f  %g %g %g %g  %g %g %g  %g %g %g \n", All.Time, task, my.pos[0], my.pos[1], my.pos[2], my.minpot, my.bh_hsml, my.bh_pos[0],
              my.bh_pos[1], my.bh_pos[2], my.bh_minpot, my.bh_rhotot, my.bh_minpotvel[0], my.bh_minpotvel[1], my.bh_minpotvel[2], my.Pos[0], my.Pos[1], my.Pos[2], my.Vel[0], my.Vel[1], my.Vel[2]);
      fflush(FdBHDiag);
    }

  myfree(mylist);
}
#endif


#endif
