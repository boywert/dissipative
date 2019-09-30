/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/GFM/stellar_feedback_find_cells.c
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


#if defined(FM_STAR_FEEDBACK) && defined(DELAYED_COOLING)

/* communication structures */
typedef struct
{
  MyDouble Pos[3];
  MyFloat BlastRadius;
  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

typedef struct
{
  MyDouble NormSphFeedback;
} data_out;

static data_out *DataResult, *DataOut;

static void particle2in(data_in * in, int i, int firstnode)
{
  in->Pos[0] = P[StarParticle[i].index].Pos[0];
  in->Pos[1] = P[StarParticle[i].index].Pos[1];
  in->Pos[2] = P[StarParticle[i].index].Pos[2];

  in->BlastRadius = STP(StarParticle[i].index).BlastRadius;

  in->Firstnode = firstnode;
}

static void out2particle(data_out * out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES)
    {
      StarParticle[i].NormSphFeedback = out->NormSphFeedback;
    }
  else
    {
      StarParticle[i].NormSphFeedback += out->NormSphFeedback;
    }
}


#include "../generic_comm_helpers2.h"

static void kernel_local(void)
{
  int i;
#ifdef GENERIC_ASYNC
  int flag = 0;
#endif

#pragma omp parallel private(idx)
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
              if(generic_polling_primary(count, Nstar))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        i = NextParticle++;

        if(i >= Nstar)
          break;

        find_feedback_cells_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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

        find_feedback_cells_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}


void find_feedback_cells(void)
{
  long long ntot;

  sumup_large_ints(1, &Nstar, &ntot);
  if(ntot == 0)
    return;

  mpi_printf("GFM_FEEDBACK: Finding cells for delayed cooling.\n");

  generic_set_MaxNexport();

  double t0 = second();

  generic_comm_pattern(Nstar, kernel_local, kernel_imported);

  double t1 = second();

  mpi_printf("GFM_FEEDBACK: Done! Calculation took %g sec\n", timediff(t0, t1));
}


/*! This function represents the core of the star density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
int find_feedback_cells_evaluate(int target, int mode, int thread_id)
{
  int numnodes, *firstnode;
  data_in local, *in;
  data_out out;

  double h, h2;
  double dx, dy, dz, r2;

#ifndef GFM_TOPHAT_KERNEL
  double wk, u, hinv, hinv3;
#endif

#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif

  MyFloat normsph;
  MyDouble *pos;


  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      in = &local;

      numnodes = 1;
      firstnode = NULL;

      pos = in->Pos;
      h = in->BlastRadius;
    }
  else
    {
      in = &DataGet[target];
      generic_get_numnodes(target, &numnodes, &firstnode);

      pos = in->Pos;
      h = in->BlastRadius;
    }

  normsph = 0;
  h2 = h * h;

#ifndef GFM_TOPHAT_KERNEL
  hinv = 1.0 / h;
#ifndef  TWODIMS
  hinv3 = hinv * hinv * hinv;
#else
  hinv3 = hinv * hinv / boxSize_Z;
#endif
#endif

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, thread_id, numnodes, firstnode);

  for(int n = 0; n < nfound; n++)
    {
      int j = Thread[thread_id].Ngblist[n];

      if(P[j].Mass > 0 && P[j].ID != 0) /* skip cells that have been swallowed or dissolved */
        {
          dx = NGB_PERIODIC_LONG_X(pos[0] - P[j].Pos[0]);
          dy = NGB_PERIODIC_LONG_Y(pos[1] - P[j].Pos[1]);
          dz = NGB_PERIODIC_LONG_Z(pos[2] - P[j].Pos[2]);

          r2 = dx * dx + dy * dy + dz * dz;

          if(r2 < h2)
            {
#ifndef GFM_TOPHAT_KERNEL
              u = sqrt(r2) * hinv;
              if(u < 0.5)
                wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
              else
                wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

              normsph += SphP[j].Volume * wk;
#else
              normsph += SphP[j].Volume;
#endif
            }
        }
    }

  out.NormSphFeedback = normsph;

  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#endif
