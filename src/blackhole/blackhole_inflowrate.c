/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/blackhole/blackhole_inflowrate.c
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


#if defined(BLACK_HOLES) && defined(BH_INFLOW_RATE)

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble Pos[3];
  MyFloat Vel[3];

  int Firstnode;
} data_in static data_in *DataIn, *DataGet;


/* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in * in, int i, int firstnode)
{
  int k;

  for(k = 0; k < 3; k++)
    {
      in->Pos[k] = P[i].Pos[k];
      in->Vel[k] = P[i].Vel[k];
    }

  in->Firstnode = firstnode;
}


 /* local data structure that holds results acquired on remote processors */
typedef struct
{
  MyFloat VolSum;
  MyFloat VrSum;
} data_out;

static data_out *DataResult, *DataOut;


 /* routine to store or combine result data */
static void out2particle(data_out * out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES)      /* initial store */
    {
      out->VolSum = vol_sum;
      out->VrSum = vr_sum;
    }
  else                          /* combine */
    {
      out->VolSum += vol_sum;
      out->VrSum += vr_sum;
    }
}


#include "../generic_comm_helpers2.h"


static MyFloat *BH_VolSum, *BH_VrSum;
static double VelToPhys, HubbleOfa, Atime, Atime2, Atime3;

static int blackhole_evaluate_inflowrate(int target, int mode, int threadid);


static void kernel_local(void)
{
  int idx;
#ifdef GENERIC_ASYNC
  int flag = 0;
#endif
#pragma omp parallel private(i, idx)
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
              if(generic_polling_primary(count, TimeBinsBHAccretion.NActiveParticles))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        idx = NextParticle++;

        if(idx >= TimeBinsBHAccretion.NActiveParticles)
          break;

        int i = TimeBinsBHAccretion.ActiveParticleList[idx];
        if(i < 0)
          continue;

        blackhole_evaluate_inflowrate(i, MODE_LOCAL_PARTICLES, threadid);
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

        blackhole_evaluate_inflowrate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}



void blackhole_inflowrate(void)
{
  int idx, i, j, k, n, nexport, nimport, ndone, ndone_flag;
  int ngrp, recvTask, dummy;
  double local_inflowrate_max = 0.0, local_inflowrate_min = MAX_REAL_NUMBER;
  double global_inflowrate_max, global_inflowrate_min;
  double local_outflowrate_max = 0.0, local_outflowrate_min = MAX_REAL_NUMBER;
  double global_outflowrate_max, global_outflowrate_min;
  double unit_fac = (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

  mpi_printf("BH_INFLOW_RATE: Measure gas inflow rate around BHs\n");

  for(idx = 0; idx < TimeBinsBHAccretion.NActiveParticles; idx++)
    {
      i = TimeBinsBHAccretion.ActiveParticleList[idx];
      if(i < 0)
        continue;

      BPP(n).BH_InflowRate = 0.0;
    }

  VelToPhys = 1.0 / All.cf_atime;
  Atime = All.cf_atime;
  Atime2 = Atime * Atime;
  Atime3 = Atime * Atime2;
  if(All.ComovingIntegrationOn)
    HubbleOfa = All.cf_hubble_a;
  else
    HubbleOfa = 0;

  BH_VolSum = (MyFloat *) mymalloc("BH_VolSum", NumBHs * sizeof(MyFloat));
  memset(BH_VolSum, 0, NumBHs * sizeof(MyDouble));
  BH_VrSum = (MyFloat *) mymalloc("BH_VrSum", NumBHs * sizeof(MyFloat));
  memset(BH_VrSum, 0, NumBHs * sizeof(MyFloat));

  generic_set_MaxNexport();

  generic_comm_pattern(TimeBinsBHAccretion.NActiveParticles, kernel_local, kernel_imported);


  for(idx = 0; idx < TimeBinsBHAccretion.NActiveParticles; idx++)
    {
      n = TimeBinsBHAccretion.ActiveParticleList[idx];
      if(n < 0)
        continue;

      if(BH_VolSum[P[n].AuxDataID] > 0)
        {
          BPP(n).BH_InflowRate = 4.0 * M_PI * BH_VrSum[P[n].AuxDataID] / BH_VolSum[P[n].AuxDataID];

          if(BPP(n).BH_InflowRate > 0)
            {
              if(BPP(n).BH_InflowRate > local_inflowrate_max)
                local_inflowrate_max = BPP(n).BH_InflowRate;
              if(BPP(n).BH_InflowRate < local_inflowrate_min)
                local_inflowrate_min = BPP(n).BH_InflowRate;
            }
          if(BPP(n).BH_InflowRate < 0)
            {
              if(fabs(BPP(n).BH_InflowRate) > local_outflowrate_max)
                local_outflowrate_max = fabs(BPP(n).BH_InflowRate);
              if(fabs(BPP(n).BH_InflowRate) < local_outflowrate_min)
                local_outflowrate_min = fabs(BPP(n).BH_InflowRate);
            }
        }
      else
        BPP(n).BH_InflowRate = 0;
    }

  myfree(BH_VrSum);
  myfree(BH_VolSum);

  MPI_Allreduce(&local_inflowrate_min, &global_inflowrate_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&local_inflowrate_max, &global_inflowrate_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  MPI_Allreduce(&local_outflowrate_min, &global_outflowrate_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&local_outflowrate_max, &global_outflowrate_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  if(global_inflowrate_max > 0.0)
    mpi_printf("BH_INFLOW_RATE: min(inflow rate) = %g M_sun/yr max(inflow rate) = %g M_sun/yr\n", global_inflowrate_min * unit_fac, global_inflowrate_max * unit_fac);
  else
    mpi_printf("BH_INFLOW_RATE: no inflows\n");

  if(global_outflowrate_max > 0.0)
    mpi_printf("BH_INFLOW_RATE: min(outflow rate) = %g M_sun/yr max(outflow rate) = %g M_sun/yr\n", global_outflowrate_min * unit_fac, global_outflowrate_max * unit_fac);
  else
    mpi_printf("BH_INFLOW_RATE: no outflows\n");

}


static int blackhole_evaluate_inflowrate(int target, int mode, int *nexport, int *nSend_local)
{
  int startnode, numngb, j, n, listindex = 0;
  MyDouble *pos;
  MyFloat *vel;
  double dx, dy, dz, dvx, dvy, dvz, r, r2;
  double mass, vol_sum, vr, vr_sum;
#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif

  vol_sum = 0.0;
  vr_sum = 0.0;

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
  vel = target_data->Vel;

  int nfound = ngb_treefind_variable_threads(pos, All.BlackHoleInflowRadius, target, mode, threadid, numnodes, firstnode);

  for(n = 0; n < nfound; n++)
    {
      j = Thread[threadid].Ngblist[n];

      /* physical */
      dx = Atime * NEAREST_X(P[j].Pos[0] - pos[0]);
      dy = Atime * NEAREST_Y(P[j].Pos[1] - pos[1]);
      dz = Atime * NEAREST_Z(P[j].Pos[2] - pos[2]);

      dvx = VelToPhys * (P[j].Vel[0] - vel[0]) + HubbleOfa * dx;
      dvy = VelToPhys * (P[j].Vel[1] - vel[1]) + HubbleOfa * dy;
      dvz = VelToPhys * (P[j].Vel[2] - vel[2]) + HubbleOfa * dz;

      r2 = dx * dx + dy * dy + dz * dz;

      r = sqrt(r2);

      if(P[j].Type == 0 && r < Atime * All.BlackHoleInflowRadius && P[j].Mass != 0 && P[j].ID != 0)
        {
          mass = P[j].Mass;
          if(r > 0)
            vr = -(dx * dvx + dy * dvy + dz * dvz) / r;
          else
            vr = 0;
          vol_sum += Atime3 * SphP[j].Volume;
          vr_sum += P[j].Mass * vr * r2;
        }
    }

  out.VolSum = vol_sum;
  out.VrSum = vr_sum;

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#endif
