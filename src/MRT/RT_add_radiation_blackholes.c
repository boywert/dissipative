/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/MRT/RT_add_radiation.c
 * \date        10/2017
 * \author      Federico Marinacci, Rahul Kannan, David Barnes
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

#include "../allvars.h"
#include "../proto.h"


#ifdef MRT_SOURCES

#ifdef MRT_BH

static int add_photons_evaluate(int target, int mode, int threadid);

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble Pos[3];
  MyFloat Hsml;
  MyFloat NormSph;
  MyFloat PhotReleased[MRT_BINS];

  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

/* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in * in, int i, int firstnode)
{
  in->Pos[0]  = P[BHParticle[i].index].Pos[0];
  in->Pos[1]  = P[BHParticle[i].index].Pos[1];
  in->Pos[2]  = P[BHParticle[i].index].Pos[2];
  in->Hsml    = BPP(BHParticle[i].index).BH_Hsml;
  in->NormSph = BHParticle[i].NormSph;

  for(int bin = 0; bin < MRT_BINS; bin++)
    in->PhotReleased[bin] = BHParticle[i].TotalPhotReleased[bin];

  in->Firstnode = firstnode;
}

typedef struct
{
  char dummy;
} data_out;

static data_out *DataResult, *DataOut;

/* routine to store or combine result data */
static void out2particle(data_out * out, int i, int mode)
{
}

#include "../generic_comm_helpers2.h"


static int Ncount;

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
              if(generic_polling_primary(count, Ncount))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        i = NextParticle++;

        if(i >= Ncount)
          break;

        add_photons_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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

        add_photons_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

void add_ionizing_blackhole_radiation(void)
{
  Ncount = Nsource;

  long long ntot;
  sumup_large_ints(1, &Nsource, &ntot);
  if(ntot == 0)
    return;

  generic_set_MaxNexport();

  double t0 = second();

  generic_comm_pattern(Nsource, kernel_local, kernel_imported);

  double t1 = second();
  mpi_printf("RT: adding ionizing radiation took %g sec\n", timediff(t0, t1));
}


static int add_photons_evaluate(int target, int mode, int threadid)
{
  int j, n;
  int numnodes, *firstnode;
  double h;
#ifndef GFM_TOPHAT_KERNEL
  double wk, u, r, hinv, hinv3;
#endif
  double weight_fac;
  MyDouble *pos;
  MyFloat *TotalPhotReleased;
  MyFloat normsph;

  data_in local, *in;
  data_out out;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      in = &local;

      numnodes = 1;
      firstnode = NULL;
    }
  else
    {
      in = &DataGet[target];

      generic_get_numnodes(target, &numnodes, &firstnode);
    }

  pos = in->Pos;
  h = in->Hsml;
  normsph = in->NormSph;
  TotalPhotReleased = in->PhotReleased;

#ifndef GFM_TOPHAT_KERNEL
  hinv = 1.0 / h;
#ifndef  TWODIMS
  hinv3 = hinv * hinv * hinv;
#else
  hinv3 = hinv * hinv / boxSize_Z;
#endif
#endif

  out.dummy = 0;

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, threadid, numnodes, firstnode);

  for(n = 0; n < nfound; n++)
    {
      j = Thread[threadid].Ngblist[n];

      if(P[j].Mass > 0 && P[j].ID != 0) /* skip cells that have been swallowed or dissolved */
        {
#ifndef GFM_TOPHAT_KERNEL
          double r2 = Thread[threadid].R2list[n];

          r = sqrt(r2);
          u = r * hinv;

          if(u < 0.5)
            wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
          else
            wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

          weight_fac = SphP[j].Volume * wk / normsph;
#else
          weight_fac = SphP[j].Volume / normsph;
#endif

          for(int bin = 0; bin < MRT_BINS; bin++)
            SphP[j].Cons_DensPhot[bin] += weight_fac * TotalPhotReleased[bin];
        }
    }


  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#endif /* MRT_BH */

#endif /* MRT_SOURCES */
