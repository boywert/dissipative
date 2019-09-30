/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/cosmic_rays/cosmic_rays_find_ngbs.c
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


#ifdef COSMIC_RAYS_SN_INJECTION

int *NgbCount;
MyFloat *NumNgb;
MyFloat *NormSph;
unsigned char *Todo;

static void find_cells_evaluate(int target, int mode, int thread_id);

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
  in->Pos[0] = P[i].Pos[0];
  in->Pos[1] = P[i].Pos[1];
  in->Pos[2] = P[i].Pos[2];

  in->Hsml = P[i].Hsml;

  in->Firstnode = firstnode;
}


 /* local data structure that holds results acquired on remote processors */
typedef struct
{
  int NgbCount;
  MyFloat NumNgb;
  MyFloat NormSph;
} data_out;

static data_out *DataResult, *DataOut;


 /* routine to store or combine result data */
static void out2particle(data_out * out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES)      /* initial store */
    {
      NgbCount[i] = out->NgbCount;
      NumNgb[i] = out->NumNgb;
      NormSph[i] = out->NormSph;
    }
  else                          /* combine */
    {
      NgbCount[i] += out->NgbCount;
      NumNgb[i] += out->NumNgb;
      NormSph[i] += out->NormSph;
    }
}

#include "../generic_comm_helpers2.h"


static void kernel_local(void)
{
  int idx;
#ifdef GENERIC_ASYNC
  int flag = 0;
#endif
  /* do local particles */
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
              if(generic_polling_primary(count, TimeBinsGravity.NActiveParticles))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        idx = NextParticle++;

        if(idx >= TimeBinsGravity.NActiveParticles)
          break;

        int i = TimeBinsGravity.ActiveParticleList[idx];
        if(i < 0)
          continue;

        if(P[i].Type != 4)
          continue;

        if(Todo[i] > 0)         /* do we already have hsml for this star? */
          find_cells_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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

        find_cells_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}


void cosmic_rays_find_ngbs(int *ngbCount, MyFloat * numNgb, MyFloat * normSph)
{
  NgbCount = ngbCount;
  NumNgb = numNgb;
  NormSph = normSph;

  MyFloat *Left, *Right;
  int idx, i, npleft, iter = 0;
  long long ntot;

  Left = (MyFloat *) mymalloc("Left", NumPart * sizeof(MyFloat));
  Right = (MyFloat *) mymalloc("Right", NumPart * sizeof(MyFloat));
  Todo = (unsigned char *) mymalloc("Todo", NumPart * sizeof(unsigned char));

  memset(Todo, 0, NumPart * sizeof(unsigned char));

  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0 || P[i].Type != 4)
        continue;

      if(P[i].CRInjection == 1)
        {
          Left[i] = Right[i] = 0;
          Todo[i] = 1;
        }
    }

  generic_set_MaxNexport();



  /* we will repeat the whole thing for those stars where we didn't find enough neighbours */
  do
    {
      generic_comm_pattern(TimeBinsGravity.NActiveParticles, kernel_local, kernel_imported);

      for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
        {
          i = TimeBinsGravity.ActiveParticleList[idx];
          if(i < 0 || P[i].Type != 4)
            continue;
          if(Todo[i] > 0)
            NumNgb[i] /= P[i].Mass;
        }

      /* do final operations on results */
      npleft = 0;
      for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
        {
          i = TimeBinsGravity.ActiveParticleList[idx];
          if(i < 0 || P[i].Type != 4)
            continue;

          if(Todo[i] > 0)
            {
              /* now check that we have at least one neighbour */
              if(NgbCount[i] == 0)
                {
                  if(Todo[i] == 2)
                    terminate("In principle this should not happen. Iter=%d, particle=%d, Right=%g, Left=%g, Hsml=%g", iter, i, Right[i], Left[i], P[i].Hsml);

                  npleft++;
                  P[i].Hsml *= 2.0;
                  continue;
                }

              /* now check whether we had enough neighbours */
              if(NumNgb[i] < (All.DesNumNgbCRInjection - All.MaxNumNgbDeviationCRInjection) || (NumNgb[i] > (All.DesNumNgbCRInjection + All.MaxNumNgbDeviationCRInjection)))
                {
                  /* need to redo this particle */
                  npleft++;

                  /* This means that the minimum Hsml has already been computed but it is still not converged (from now on the particle has at least 1 neighbour) */
                  Todo[i] = 2;

                  if(Left[i] > 0 && Right[i] > 0)
                    if((Right[i] - Left[i]) < 1.0e-3 * Left[i])
                      {
                        /* this one should be ok */
                        npleft--;
                        Todo[i] = 0;    /* done */
                        continue;
                      }

                  if(NumNgb[i] < (All.DesNumNgbCRInjection - All.MaxNumNgbDeviationCRInjection))
                    Left[i] = dmax(P[i].Hsml, Left[i]);
                  else
                    {
                      if(Right[i] != 0)
                        {
                          if(P[i].Hsml < Right[i])
                            Right[i] = P[i].Hsml;
                        }
                      else
                        Right[i] = P[i].Hsml;

                      if(NgbCount[i] == 1)
                        {
                          /* this one should be ok. Even with only one neighbour we are above the target mass (low-res region) */
                          npleft--;
                          Todo[i] = 0;  /* done */
                          continue;
                        }
                    }

                  if(iter >= MAXITER - 10)
                    {
                      printf("i=%d task=%d ID=%d Hsml=%g Left=%g Right=%g Ngbs=%g, NgbsCount=%d Right-Left=%g\n   pos=(%g|%g|%g)\n",
                             i, ThisTask, (int) P[i].ID, P[i].Hsml, Left[i], Right[i], (float) NumNgb[i], NgbCount[i], Right[i] - Left[i], P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
                      myflush(stdout);
                    }

                  if(Right[i] > 0 && Left[i] > 0)
                    P[i].Hsml = pow(0.5 * (pow(Left[i], 3) + pow(Right[i], 3)), 1.0 / 3);
                  else
                    {
                      if(Right[i] == 0 && Left[i] == 0)
                        terminate("should not occur");  /* can't occur */

                      if(Right[i] == 0 && Left[i] > 0)
                        P[i].Hsml *= 1.26;

                      if(Right[i] > 0 && Left[i] == 0)
                        P[i].Hsml /= 1.26;
                    }
                }
              else
                Todo[i] = 0;    /* done */
            }
        }

      sumup_large_ints(1, &npleft, &ntot);

      if(ntot > 0)
        {
          iter++;

          if(iter > 0)
            mpi_printf("COSMIC_RAYS: star ngb iteration %3d: need to repeat for %12lld particles.\n", iter, ntot);

          if(iter > MAXITER)
            terminate("failed to converge in neighbour iteration in cosmic_rays_find_ngbs()\n");
        }
    }
  while(ntot > 0);

  myfree(Todo);
  myfree(Right);
  myfree(Left);
}

void find_cells_evaluate(int target, int mode, int thread_id)
{
  int j, n;
  int numngb, numnodes, *firstnode;
  double h, h2;
  double wk, hinv, hinv3;
  double u;
  double dx, dy, dz, r, r2;
  MyFloat weighted_numngb, normsph;
  MyDouble *pos;
  double mass_j;
#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
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
  h = target_data->Hsml;


  h2 = h * h;

  hinv = 1.0 / h;
#ifndef TWODIMS
  hinv3 = hinv * hinv * hinv;
#else
  hinv3 = hinv * hinv / boxSize_Z;
#endif

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, thread_id, numnodes, firstnode);

  numngb = 0;
  normsph = 0;
  weighted_numngb = 0;

  for(n = 0; n < nfound; n++)
    {
      j = Thread[thread_id].Ngblist[n];

      if(P[j].Mass > 0 && P[j].ID != 0)
        {
          dx = NGB_PERIODIC_LONG_X(pos[0] - P[j].Pos[0]);
          dy = NGB_PERIODIC_LONG_Y(pos[1] - P[j].Pos[1]);
          dz = NGB_PERIODIC_LONG_Z(pos[2] - P[j].Pos[2]);

          r2 = dx * dx + dy * dy + dz * dz;

          if(r2 < h2)
            {
              numngb++;

              r = sqrt(r2);
              u = r * hinv;

              if(u < 0.5)
                wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
              else
                wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

              mass_j = P[j].Mass;
              weighted_numngb += FLT(NORM_COEFF * mass_j * wk / hinv3); /* 4.0/3 * PI = 4.188790204786, swallowed cells are not counted here, because they have mass_j=0 */
              normsph += SphP[j].Volume * wk;
            }
        }
    }

  out.NgbCount = numngb;
  out.NumNgb = weighted_numngb;
  out.NormSph = normsph;

  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;
}

#endif
