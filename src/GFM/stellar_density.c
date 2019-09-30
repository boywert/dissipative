/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/GFM/stellar_density.c
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


#ifdef GFM_STELLAR_EVOLUTION

static int find_cells_to_enrich_evaluate(int target, int mode, int thread_id);


#ifdef DETAILEDTIMINGS
static double tstart;
static int current_timebin;
#endif

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble Pos[3];
  MyFloat Hsml;

  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

 /* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in * in, int i, int firstnode)
{
  in->Pos[0] = P[StarParticle[i].index].Pos[0];
  in->Pos[1] = P[StarParticle[i].index].Pos[1];
  in->Pos[2] = P[StarParticle[i].index].Pos[2];

  in->Hsml = STP(StarParticle[i].index).Hsml;

  in->Firstnode = firstnode;
}


 /* local data structure that holds results acquired on remote processors */
typedef struct
{
  MyFloat NumNgb;
  MyFloat NormSph;
  MyFloat Dhsmlrho;
  MyFloat ClosestNeighbourDistance;
#ifdef FM_VAR_SN_EFF
  MyFloat AvgMetalNgb;
#endif
#ifdef FM_MASS_WEIGHT_SN
  MyFloat TotNgbMass;
#endif
#ifdef DELAYED_COOLING
  MyDouble AvgPress;
  MyDouble AvgHDens;
  MyFloat MinBlastRadius;
#endif
#ifdef FM_RADIATION_FEEDBACK
  MyFloat RadFeed_MinGasDist;
#endif
#if defined(NON_STOCHASTIC_MOMENTUM_FEEDBACK) && defined(INJECT_INTO_SINGLE_CELL)
  MyIDType ClosestNeighbourID;
#endif
#if defined(FM_SN_COOLING_RADIUS_BOOST) || defined(FM_RADIATION_FEEDBACK)
  MyFloat LocISMdens; 
  MyFloat LocISMZdens;
#endif
} data_out;

static data_out *DataResult, *DataOut;


 /* routine to store or combine result data */
static void out2particle(data_out * out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES)      /* initial store */
    {
      StarParticle[i].NumNgb = out->NumNgb;
      StarParticle[i].NormSph = out->NormSph;
      StarParticle[i].Dhsmlrho = out->Dhsmlrho;
      StarParticle[i].ClosestNeighbourDistance = out->ClosestNeighbourDistance;
#if defined(NON_STOCHASTIC_MOMENTUM_FEEDBACK) && defined(INJECT_INTO_SINGLE_CELL)
      StarParticle[i].ClosestNeighbourID = out->ClosestNeighbourID;
#endif

#ifdef DELAYED_COOLING
      StarParticle[i].AvgPress = out->AvgPress;
      StarParticle[i].AvgHDens = out->AvgHDens;
      StarParticle[i].MinBlastRadius = out->MinBlastRadius;
#endif
#ifdef FM_VAR_SN_EFF
      StarParticle[i].AvgMetalNgb = out->AvgMetalNgb;
#endif
#ifdef FM_MASS_WEIGHT_SN
      StarParticle[i].TotNgbMass = out->TotNgbMass;
#endif
#ifdef  FM_RADIATION_FEEDBACK
      StarParticle[i].RadFeed_MinGasDist = out->RadFeed_MinGasDist;
#endif
#if defined(FM_STAR_FEEDBACK) && (defined(FM_SN_COOLING_RADIUS_BOOST) || defined(FM_RADIATION_FEEDBACK))
      StarParticle[i].LocISMdens  = out->LocISMdens;
      StarParticle[i].LocISMZdens = out->LocISMZdens;
#endif
    }
  else                          /* combine */
    {
      StarParticle[i].NumNgb += out->NumNgb;
      StarParticle[i].NormSph += out->NormSph;
      StarParticle[i].Dhsmlrho += out->Dhsmlrho;

      if(out->ClosestNeighbourDistance < StarParticle[i].ClosestNeighbourDistance)
        {
          StarParticle[i].ClosestNeighbourDistance = out->ClosestNeighbourDistance;
#if defined(NON_STOCHASTIC_MOMENTUM_FEEDBACK) && defined(INJECT_INTO_SINGLE_CELL)
          StarParticle[i].ClosestNeighbourID = out->ClosestNeighbourID;
#endif
        }

#ifdef DELAYED_COOLING
      StarParticle[i].AvgPress += out->AvgPress;
      StarParticle[i].AvgHDens += out->AvgHDens;
      StarParticle[i].MinBlastRadius = dmin(out->MinBlastRadius, StarParticle[i].MinBlastRadius);
#endif
#ifdef FM_VAR_SN_EFF
      StarParticle[i].AvgMetalNgb += out->AvgMetalNgb;
#endif
#ifdef FM_MASS_WEIGHT_SN
      StarParticle[i].TotNgbMass += out->TotNgbMass;
#endif
#ifdef FM_RADIATION_FEEDBACK
      StarParticle[i].RadFeed_MinGasDist = dmin(out->RadFeed_MinGasDist, StarParticle[i].RadFeed_MinGasDist);
#endif
#if defined(FM_STAR_FEEDBACK) && (defined(FM_SN_COOLING_RADIUS_BOOST) || defined(FM_RADIATION_FEEDBACK))
      StarParticle[i].LocISMdens  += out->LocISMdens;
      StarParticle[i].LocISMZdens += out->LocISMZdens;
#endif
    }
}

#include "../generic_comm_helpers2.h"


static int Npart;
static MyFloat *Left, *Right;
static unsigned char *Todo;

static void kernel_local(void)
{
#ifdef DETAILEDTIMINGS
  double t0 = second();
#endif

  /* do local particles */
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
              if(generic_polling_primary(count, Npart))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        i = NextParticle++;

        if(i >= Npart)
          break;

        if(Todo[i] > 0)         /* do we already have hsml for this star? */
          {
            int p = StarParticle[i].index;

            if(P[p].Ti_Current != All.Ti_Current)
              {
                terminate("we should not get here");
#if (NUM_THREADS > 1)
                omp_set_lock(&ParticleLocks[p]);

                if(P[p].Ti_Current != All.Ti_Current)
                  {
#endif
                    drift_particle(p, All.Ti_Current);
#if (NUM_THREADS > 1)
                  }
                omp_unset_lock(&ParticleLocks[p]);
#endif
              }

            find_cells_to_enrich_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
          }
      }
  }


#ifdef DETAILEDTIMINGS
  double t1 = second();
  fprintf(FdDetailed, "%d %d %d %d %g %g\n", All.NumCurrentTiStep, current_timebin, DETAILED_TIMING_STELLARDENSITY,
          MODE_LOCAL_PARTICLES, timediff(tstart, t0), timediff(tstart, t1));
#endif
}

static void kernel_imported(void)
{
#ifdef DETAILEDTIMINGS
  double t0 = second();
#endif

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

        find_cells_to_enrich_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }

#ifdef DETAILEDTIMINGS
  double t1 = second();
  fprintf(FdDetailed, "%d %d %d %d %g %g\n", All.NumCurrentTiStep, current_timebin, DETAILED_TIMING_STELLARDENSITY,
          MODE_IMPORTED_PARTICLES, timediff(tstart, t0), timediff(tstart, t1));
#endif
}

void find_cells_dump(int npart)
{
  int i, npleft, iter = 0;
  long long ntot, npartall;

  Npart = npart;
  sumup_large_ints(1, &npart, &npartall);
  if(npartall == 0)
    return;

  CPU_Step[CPU_MISC] += measure_time();
  double t0 = second();

  Left = (MyFloat *) mymalloc("Left", npart * sizeof(MyFloat));
  Right = (MyFloat *) mymalloc("Right", npart * sizeof(MyFloat));
  Todo = (unsigned char *) mymalloc("Todo", npart * sizeof(unsigned char));

  for(i = 0; i < npart; i++)
    {
      Left[i] = Right[i] = 0;
      Todo[i] = 1;
    }

  generic_set_MaxNexport();

#ifdef DETAILEDTIMINGS
  tstart = second();
  current_timebin = All.HighestActiveTimeBin;
#endif

  /* we will repeat the whole thing for those stars where we didn't find enough neighbours */
  do
    {
      double tA = second();

      generic_comm_pattern(Npart, kernel_local, kernel_imported);

      /* do final operations on results */
      for(i = 0, npleft = 0; i < npart; i++)
        {
          if(Todo[i] > 0)
            {
              /* now check whether we had enough neighbours */
              if(StarParticle[i].NumNgb < (All.DesNumNgbEnrichment - All.MaxNumNgbDeviationEnrichment) || (StarParticle[i].NumNgb > (All.DesNumNgbEnrichment + All.MaxNumNgbDeviationEnrichment)))
                {
                  /* need to redo this particle */
                  npleft++;

                  if(StarParticle[i].NumNgb > 0)
                    {
                      StarParticle[i].Dhsmlrho *= STP(StarParticle[i].index).Hsml / (NUMDIMS * StarParticle[i].NumNgb / (NORM_COEFF * pow(STP(StarParticle[i].index).Hsml, 3)));

                      if(StarParticle[i].Dhsmlrho > -0.9)       /* note: this would be -1 if only a single particle at zero lag is found */
                        StarParticle[i].Dhsmlrho = 1 / (1 + StarParticle[i].Dhsmlrho);
                      else
                        StarParticle[i].Dhsmlrho = 1;
                    }
                  else
                    StarParticle[i].Dhsmlrho = 1;

#ifndef GFM_EXACT_NUMNGB
                  if(Left[i] > 0 && Right[i] > 0)
                    if((Right[i] - Left[i]) < 1.0e-3 * Left[i])
                      {
                        /* this one should be ok */
                        npleft--;
                        Todo[i] = 0;    /* done */
                        continue;
                      }
#endif

                  if(StarParticle[i].NumNgb < (All.DesNumNgbEnrichment - All.MaxNumNgbDeviationEnrichment))
                    Left[i] = dmax(STP(StarParticle[i].index).Hsml, Left[i]);
                  else
                    {
                      if(Right[i] != 0)
                        {
                          if(STP(StarParticle[i].index).Hsml < Right[i])
                            Right[i] = STP(StarParticle[i].index).Hsml;
                        }
                      else
                        Right[i] = STP(StarParticle[i].index).Hsml;
                    }

                  if(iter >= MAXITER - 10)
                    {
                      printf("i=%d task=%d ID=%d Hsml=%g Left=%g Right=%g Ngbs=%g  Right-Left=%g\n   pos=(%g|%g|%g)\n",
                             i, ThisTask, (int) P[StarParticle[i].index].ID, STP(StarParticle[i].index).Hsml, Left[i],
                             Right[i], (float) StarParticle[i].NumNgb, Right[i] - Left[i], P[StarParticle[i].index].Pos[0], P[StarParticle[i].index].Pos[1], P[StarParticle[i].index].Pos[2]);
                      myflush(stdout);
                    }

                  if(Right[i] > 0 && Left[i] > 0)
                    STP(StarParticle[i].index).Hsml = pow(0.5 * (pow(Left[i], 3) + pow(Right[i], 3)), 1.0 / 3);
                  else
                    {
                      if(Right[i] == 0 && Left[i] == 0)
                        terminate("should not occur");  /* can't occur */

                      if(Right[i] == 0 && Left[i] > 0)
                        {
                          double fac = 1.26;

                          if(fabs(StarParticle[i].NumNgb - All.DesNumNgbEnrichment) < 0.5 * All.DesNumNgbEnrichment)
                            {
                              fac = 1 - (StarParticle[i].NumNgb - All.DesNumNgbEnrichment) / (NUMDIMS * StarParticle[i].NumNgb) * StarParticle[i].Dhsmlrho;

                              if(fac > 1.26)
                                fac = 1.26;
                            }

                          STP(StarParticle[i].index).Hsml *= fac;
                        }

                      if(Right[i] > 0 && Left[i] == 0)
                        {
                          double fac = 1 / 1.26;

                          if(fabs(StarParticle[i].NumNgb - All.DesNumNgbEnrichment) < 0.5 * All.DesNumNgbEnrichment)
                            {
                              fac = 1 - (StarParticle[i].NumNgb - All.DesNumNgbEnrichment) / (NUMDIMS * StarParticle[i].NumNgb) * StarParticle[i].Dhsmlrho;

                              if(fac < 1 / 1.26)
                                fac = 1 / 1.26;
                            }
                          STP(StarParticle[i].index).Hsml *= fac;
                        }
                    }
                }
              else
                Todo[i] = 0;    /* done */
            }
        }

      sumup_large_ints(1, &npleft, &ntot);

      double tB = second();

      if(ntot > 0)
        {
          iter++;

          if(iter > 0)
            mpi_printf("GFM_STELLAR_EVOLUTION: star ngb iteration %3d: need to repeat for %12lld particles. (previous iteration took %g sec)\n", iter, ntot, timediff(tA, tB));

          if(iter > MAXITER)
            terminate("failed to converge in neighbour iteration in find_cells_dump()\n");
        }
    }
  while(ntot > 0);

  myfree(Todo);
  myfree(Right);
  myfree(Left);

#if defined(FM_SN_COOLING_RADIUS_BOOST) || defined(FM_VAR_SN_EFF) || defined(FM_MASS_WEIGHT_SN)
  for(i = 0; i < npart; i++)
    {
#ifdef FM_VAR_SN_EFF
#ifdef FM_MASS_WEIGHT_SN
         StarParticle[i].AvgMetalNgb /= StarParticle[i].TotNgbMass;
#else
         StarParticle[i].AvgMetalNgb /= StarParticle[i].NumNgb;
#endif
         STP(StarParticle[i].index).AvgMetalNgb = StarParticle[i].AvgMetalNgb;
#endif

#ifdef FM_MASS_WEIGHT_SN
         STP(StarParticle[i].index).TotNgbMass = StarParticle[i].TotNgbMass;
#endif
#if defined(FM_STAR_FEEDBACK) && (defined(FM_SN_COOLING_RADIUS_BOOST) || defined(FM_RADIATION_FEEDBACK))
         STP(StarParticle[i].index).LocISMdens = StarParticle[i].LocISMdens;
         STP(StarParticle[i].index).LocISMZdens= StarParticle[i].LocISMZdens;
#endif
    }
#endif

#ifdef DETAILEDTIMINGS
  double tend = second();
  fprintf(FdDetailed, "%d %d %d %d %g %g\n", All.NumCurrentTiStep, current_timebin, DETAILED_TIMING_STELLARDENSITY,
           MODE_FINISHED, timediff(tstart, tend), timediff(tstart, tend));
  fflush(FdDetailed);
#endif

  CPU_Step[CPU_GFM_ENRICH] += measure_time();
  double t1 = second();

  mpi_printf("GFM_STELLAR_EVOLUTION: active particles %lld, stellar density iterations took = %g sec\n", npartall, timediff(t0, t1));
}


/*! This function represents the core of the star density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
int find_cells_to_enrich_evaluate(int target, int mode, int thread_id)
{
  int numnodes, *firstnode;
  double wk, dwk;
#ifdef DELAYED_COOLING
  double avg_pressure = 0;
  double avg_H_density = 0;
  double XH = HYDROGEN_MASSFRAC;
  double min_blast_radius = MAX_REAL_NUMBER;
#endif
#ifdef FM_VAR_SN_EFF
  double sum_ngb_metallicity = 0;
#endif
#ifdef FM_MASS_WEIGHT_SN
  double tot_ngb_mass = 0;
#endif
#if defined(FM_SN_COOLING_RADIUS_BOOST) || defined(FM_RADIATION_FEEDBACK)
  double loc_ism_dens = 0;
  double loc_ism_z_dens = 0;
#endif

#if defined(NON_STOCHASTIC_MOMENTUM_FEEDBACK) && defined(INJECT_INTO_SINGLE_CELL)
  MyIDType minDistGasID = 0;
#endif

  double minDistGas = MAX_DOUBLE_NUMBER;
  double tol = 1.0e-2;

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

  MyDouble *pos = target_data->Pos;
  double h = target_data->Hsml;

  double hinv = 1.0 / h;
#ifndef  TWODIMS
  double hinv3 = hinv * hinv * hinv;
#else
  double hinv3 = hinv * hinv / boxSize_Z;
#endif

  double h3 = 1.0 / hinv3;

#ifdef DELAYED_COOLING
  wk = hinv3 / NORM_COEFF;
#endif

  double hinv4 = hinv3 * hinv;

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, thread_id, numnodes, firstnode);

  double normsph = 0;
  double weighted_numngb = 0;
  double dhsmlrho = 0;
#ifdef GFM_EXACT_NUMNGB
  int numngb = 0;
#endif

  for(int n = 0; n < nfound; n++)
    {
      int j = Thread[thread_id].Ngblist[n];

      if(P[j].Mass > 0 && P[j].ID != 0)
        {
          double r2 = Thread[thread_id].R2list[n];

#ifdef GFM_EXACT_NUMNGB
          numngb++;
#endif
          double r = sqrt(r2);

          double u = r * hinv;

          if(u < 0.5)
            {
              wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
              dwk = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
            }
          else
            {
              wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
              dwk = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
            }

          weighted_numngb += NORM_COEFF * wk * h3;      /* 4.0/3 * PI = 4.188790204786, swallowed cells are not counted here, because they have mass_j=0 */
          dhsmlrho += (-(NUMDIMS * hinv * wk + u * dwk));

#if defined(FM_SN_COOLING_RADIUS_BOOST) || defined(FM_RADIATION_FEEDBACK)
          loc_ism_dens   += NORM_COEFF * P[j].Mass * wk;
          loc_ism_z_dens += NORM_COEFF * P[j].Mass * wk * SphP[j].Metallicity;
#endif

#ifndef GFM_TOPHAT_KERNEL
#ifdef FM_MASS_WEIGHT_SN
          tot_ngb_mass += NORM_COEFF * P[j].Mass * wk * h3;
#endif
          normsph += SphP[j].Volume * wk;
#else
#ifdef FM_MASS_WEIGHT_SN
          tot_ngb_mass += NORM_COEFF * P[j].Mass * wk * h3; /* this ensures that the feedback is spread the same way the neighbourse are counted */
#endif
          normsph += SphP[j].Volume;
#endif

#ifdef DELAYED_COOLING
          avg_pressure += GAMMA_MINUS1 * P[j].Mass * SphP[j].Utherm * wk;
#ifdef GFM_COOLING_METAL
          XH = SphP[j].MetalsFraction[element_index_Hydrogen];
#endif
          avg_H_density += P[j].Mass * XH / PROTONMASS * wk;
          min_blast_radius = dmin(min_blast_radius, r);
#endif
#ifdef FM_VAR_SN_EFF
#ifdef FM_MASS_WEIGHT_SN
          sum_ngb_metallicity += SphP[j].Metallicity * NORM_COEFF * P[j].Mass * wk * h3;
#else
          sum_ngb_metallicity += SphP[j].Metallicity * NORM_COEFF * wk * h3;
#endif
#endif

          if(r < minDistGas)
            {
              minDistGas = r;
#if defined(NON_STOCHASTIC_MOMENTUM_FEEDBACK) && defined(INJECT_INTO_SINGLE_CELL)
              minDistGasID = P[j].ID;
#endif
            }
        }
    }

#ifdef GFM_EXACT_NUMNGB
  weighted_numngb = numngb;
#endif

  out.NumNgb = weighted_numngb;
  out.NormSph = normsph;
  out.Dhsmlrho = dhsmlrho;
#ifdef DELAYED_COOLING
  out.AvgPress = avg_pressure;
  out.AvgHDens = avg_H_density;
  out.MinBlastRadius = min_blast_radius;
#endif
#ifdef FM_VAR_SN_EFF
  out.AvgMetalNgb = sum_ngb_metallicity;
#endif
#ifdef FM_MASS_WEIGHT_SN
  out.TotNgbMass = tot_ngb_mass;
#endif
#ifdef FM_RADIATION_FEEDBACK
  out.RadFeed_MinGasDist = minDistGas;
#endif
#if defined(FM_SN_COOLING_RADIUS_BOOST) || defined(FM_RADIATION_FEEDBACK)
  out.LocISMdens  = loc_ism_dens; 
  out.LocISMZdens = loc_ism_z_dens; 
#endif
#if defined(NON_STOCHASTIC_MOMENTUM_FEEDBACK) && defined(INJECT_INTO_SINGLE_CELL)
  out.ClosestNeighbourID = minDistGasID;
#endif
  out.ClosestNeighbourDistance = (1.0 + tol) * minDistGas;

  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}


void find_cells_to_enrich(void)
{
  find_cells_dump(Nstar);
}

#ifdef GFM_WINDS_STRIPPING
void find_cells_to_strip(void)
{
  find_cells_dump(Nwinds_to_strip);
}

#endif

#endif
