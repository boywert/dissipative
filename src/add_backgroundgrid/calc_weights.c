/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/add_backgroundgrid/calc_weights.c
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
#include "../allvars.h"
#include "../proto.h"
#include "../domain.h"
#include "add_bggrid.h"

#ifdef ADDBACKGROUNDGRID

static int find_cells_evaluate(int target, int mode, int thread_id);

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

  in->Hsml = SphP[i].Hsml;

  in->Firstnode = firstnode;
}


 /* local data structure that holds results acquired on remote processors */
typedef struct
{
  MyFloat Weight;
} data_out;

static data_out *DataResult, *DataOut;


 /* routine to store or combine result data */
static void out2particle(data_out * out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES)      /* initial store */
    {
      SphP[i].Weight = out->Weight;
    }
  else                          /* combine */
    {
      SphP[i].Weight += out->Weight;
    }
}

#include "../generic_comm_helpers2.h"

static void kernel_local(void)
{
  int idx;
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


void calculate_weights()
{
  domain_free();
  domain_Decomposition();       /* do new domain decomposition, will also make a new chained-list of synchronized particles */

  ngb_treeallocate();
  ngb_treebuild(NumGas);

  mpi_printf("ADD BACKGROUND GRID: distribution of fluid quantities in a SPH-like fashion\n");
  mpi_printf("ADD BACKGROUND GRID: finding the normalization factors\n");

  TimeBinsGravity.NActiveParticles = 0;

  int i;
  for(i = 0; i < NumGas; i++)
    {
      if(P[i].Mass > 0)
        {
          TimeBinsGravity.ActiveParticleList[TimeBinsGravity.NActiveParticles] = i;
          TimeBinsGravity.NActiveParticles++;
        }
    }

  generic_set_MaxNexport();

  generic_comm_pattern(TimeBinsGravity.NActiveParticles, kernel_local, kernel_imported);

  mpi_printf("ADD BACKGROUND GRID: done\n");
}

int find_cells_evaluate(int target, int mode, int thread_id)
{
  int j, n, numnodes, *firstnode;
  double h, h2, hinv, hinv3;
  MyDouble dx, dy, dz, r;
  MyDouble *pos;
#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif

  double weight = 0;

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

  for(n = 0; n < nfound; n++)
    {
      j = Thread[thread_id].Ngblist[n];

      if(P[j].ID >= IDNew)
        {
          dx = NGB_PERIODIC_LONG_X(pos[0] - P[j].Pos[0]);
          dy = NGB_PERIODIC_LONG_Y(pos[1] - P[j].Pos[1]);
          dz = NGB_PERIODIC_LONG_Z(pos[2] - P[j].Pos[2]);

          double r2 = dx * dx + dy * dy + dz * dz;

          if(r2 < h2)
            {
              r = sqrt(r2);

              double u = r * hinv;
              double wk;
              if(u < 0.5)
                wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
              else
                wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

              weight += wk * SphP[j].Volume;
            }
        }
    }

  out.Weight = weight;

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#endif
