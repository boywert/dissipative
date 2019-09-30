/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dust_live/dust_density.c
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

#ifdef DUST_LIVE

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble Pos[3];
  MyFloat Hsml;
  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

/* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in *in, int i, int firstnode)
{
  in->Pos[0] = P[DustParticle[i].index].Pos[0];
  in->Pos[1] = P[DustParticle[i].index].Pos[1];
  in->Pos[2] = P[DustParticle[i].index].Pos[2];
  in->Hsml = DTP(DustParticle[i].index).Hsml;
  in->Firstnode = firstnode;
}

/* local data structure that holds results acquired on remote processors */
typedef struct
{
  MyFloat NumNgb;
  MyFloat NormSph;
  MyFloat TotNgbMass;
  MyFloat Dhsmlrho;
  int MinGasTimeBin;
} data_out;

static data_out *DataResult, *DataOut;

 /* routine to store or combine result data */
static void out2particle(data_out * out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES)      /* initial store */
    {
      DustParticle[i].NumNgb = out->NumNgb;
      DustParticle[i].NormSph = out->NormSph;
      DustParticle[i].TotNgbMass = out->TotNgbMass;
      DustParticle[i].Dhsmlrho = out->Dhsmlrho;
      DTP(DustParticle[i].index).MinGasTimeBin = out->MinGasTimeBin;
    }
  else                          /* combine */
    {
      DustParticle[i].NumNgb += out->NumNgb;
      DustParticle[i].NormSph += out->NormSph;
      DustParticle[i].TotNgbMass += out->TotNgbMass;
      DustParticle[i].Dhsmlrho += out->Dhsmlrho;
      DTP(DustParticle[i].index).MinGasTimeBin = imin(DTP(DustParticle[i].index).MinGasTimeBin, out->MinGasTimeBin);
    }
}

#include "../generic_comm_helpers2.h"

static int Npart;
static MyFloat *Left, *Right;
static unsigned char *Todo;

static void kernel_local(void)
{
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

        if(Todo[i] > 0) /* do we already have hsml for this dust particle? */
          {
            int p = DustParticle[i].index;

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

            find_drag_cells_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
          }
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

        find_drag_cells_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

void find_drag_cells(int npart)
{
  int i, npleft, iter = 0;
  long long ntot, npartall;

  Npart = npart;
  sumup_large_ints(1, &npart, &npartall);
  if(npartall == 0)
    return;

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

  /* we will repeat the whole thing for those dust particles where we didn't find enough neighbours */
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
              if(DustParticle[i].NumNgb < (All.DesNumNgbDust - All.MaxNumNgbDeviationDust) || (DustParticle[i].NumNgb > (All.DesNumNgbDust + All.MaxNumNgbDeviationDust)))
                {
                  /* need to redo this particle */
                  npleft++;

                  if(DustParticle[i].NumNgb > 0)
                    {
                      DustParticle[i].Dhsmlrho *= DTP(DustParticle[i].index).Hsml / (NUMDIMS * DustParticle[i].NumNgb / (NORM_COEFF * pow(DTP(DustParticle[i].index).Hsml, 3)));

                      if(DustParticle[i].Dhsmlrho > -0.9) /* note: this would be -1 if only a single particle at zero lag is found */
                        DustParticle[i].Dhsmlrho = 1 / (1 + DustParticle[i].Dhsmlrho);
                      else
                        DustParticle[i].Dhsmlrho = 1;
                    }
                  else
                    DustParticle[i].Dhsmlrho = 1;

                  if(Left[i] > 0 && Right[i] > 0)
                    if((Right[i] - Left[i]) < 1.0e-3 * Left[i])
                      {
                        /* this one should be ok */
                        npleft--;
                        Todo[i] = 0; /* done */
                        continue;
                      }

                  if(DustParticle[i].NumNgb < (All.DesNumNgbDust - All.MaxNumNgbDeviationDust))
                    Left[i] = dmax(DTP(DustParticle[i].index).Hsml, Left[i]);
                  else
                    {
                      if(Right[i] != 0)
                        {
                          if(DTP(DustParticle[i].index).Hsml < Right[i])
                            Right[i] = DTP(DustParticle[i].index).Hsml;
                        }
                      else
                        Right[i] = DTP(DustParticle[i].index).Hsml;
                    }

                  if(iter >= MAXITER - 10)
                    {
                      printf("i=%d task=%d ID=%d Hsml=%g Left=%g Right=%g Ngbs=%g  Right-Left=%g\n   pos=(%g|%g|%g)\n",
                             i, ThisTask, (int) P[DustParticle[i].index].ID, DTP(DustParticle[i].index).Hsml, Left[i],
                             Right[i], (float) DustParticle[i].NumNgb, Right[i] - Left[i], P[DustParticle[i].index].Pos[0], P[DustParticle[i].index].Pos[1], P[DustParticle[i].index].Pos[2]);
                      myflush(stdout);
                    }

                  if(Right[i] > 0 && Left[i] > 0)
                    DTP(DustParticle[i].index).Hsml = pow(0.5 * (pow(Left[i], 3) + pow(Right[i], 3)), 1.0 / 3);
                  else
                    {
                      if(Right[i] == 0 && Left[i] == 0)
                        terminate("should not occur"); /* can't occur */

                      if(Right[i] == 0 && Left[i] > 0)
                        {
                          double fac = 1.26;

                          if(fabs(DustParticle[i].NumNgb - All.DesNumNgbDust) < 0.5 * All.DesNumNgbDust)
                            {
                              fac = 1 - (DustParticle[i].NumNgb - All.DesNumNgbDust) / (NUMDIMS * DustParticle[i].NumNgb) * DustParticle[i].Dhsmlrho;

                              if(fac > 1.26)
                                fac = 1.26;
                            }

                          DTP(DustParticle[i].index).Hsml *= fac;
                        }

                      if(Right[i] > 0 && Left[i] == 0)
                        {
                          double fac = 1 / 1.26;

                          if(fabs(DustParticle[i].NumNgb - All.DesNumNgbDust) < 0.5 * All.DesNumNgbDust)
                            {
                              fac = 1 - (DustParticle[i].NumNgb - All.DesNumNgbDust) / (NUMDIMS * DustParticle[i].NumNgb) * DustParticle[i].Dhsmlrho;

                              if(fac < 1 / 1.26)
                                fac = 1 / 1.26;
                            }
                          DTP(DustParticle[i].index).Hsml *= fac;
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
            mpi_printf("DUST_LIVE: dust ngb iteration %3d: need to repeat for %12lld particles. (previous iteration took %g sec)\n", iter, ntot, timediff(tA, tB));

          if(iter > MAXITER)
            terminate("failed to converge in neighbour iteration in find_drag_cells()\n");
        }
    }
  while(ntot > 0);

  myfree(Todo);
  myfree(Right);
  myfree(Left);

  double t1 = second();

  mpi_printf("DUST_LIVE: active particles %lld, dust density iterations took = %g sec\n", npartall, timediff(t0, t1));
}


/*! This function represents the core of the dust neighbor density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
int find_drag_cells_evaluate(int target, int mode, int thread_id)
{
  int numnodes, *firstnode;
  double wk, dwk;

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
  double hinv4 = hinv3 * hinv;

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, thread_id, numnodes, firstnode);

  double normsph = 0;
  double tot_ngb_mass = 0.0;
  double weighted_numngb = 0;
  double dhsmlrho = 0;
  int min_gas_time_bin = TIMEBINS;

  for(int n = 0; n < nfound; n++)
    {
      int j = Thread[thread_id].Ngblist[n];

      if(P[j].Mass > 0 && P[j].ID != 0)
        {
          double r2 = Thread[thread_id].R2list[n];

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

#ifndef GFM_TOPHAT_KERNEL
          tot_ngb_mass += NORM_COEFF * P[j].Mass * wk * h3;
          normsph += SphP[j].Volume * wk;
#else
          tot_ngb_mass += P[j].Mass;
          normsph += SphP[j].Volume;
#endif

          min_gas_time_bin = imin(min_gas_time_bin, P[j].TimeBinHydro);
        }
    }

  out.NumNgb = weighted_numngb;
  out.NormSph = normsph;
  out.TotNgbMass = tot_ngb_mass;
  out.Dhsmlrho = dhsmlrho;
  out.MinGasTimeBin = min_gas_time_bin;

  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#endif
