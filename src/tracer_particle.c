/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/tracer_particle.c
 * \date        08/2010
 * \author      Mark Vogelsberger, Shy Genel, Dylan Nelson
 * \brief       Velocity Field Tracer Particles
 * \details     
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>

#include "allvars.h"
#include "proto.h"
#include "voronoi.h"

#ifdef TRACER_PARTICLE

static void find_primary_cell_evaluate(int target, int mode, int threadid);

/* temporary particle arrays */
static MyDouble *tracer_nearest_dist;
static MyDouble MinTracerHsml;
static int *tracer_newTimeBin;

typedef struct
{
  MyDouble pos[3];              /* tracer particle position */
  MyFloat hsml;                 /* current search radius */

  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

static void particle2in(data_in * in, int i, int firstnode)
{
  in->pos[0] = P[i].Pos[0];
  in->pos[1] = P[i].Pos[1];
  in->pos[2] = P[i].Pos[2];

  in->hsml = P[i].TracerHsml;

  in->Firstnode = firstnode;
};

typedef struct
{
  MyDouble nearestDist;         /* distance to closest cell on task */
  MyFloat vel[3];
  int newTimeBin;;
#ifdef TRACER_PART_NUM_FLUID_QUANTITIES
  MyFloat fluid_quantities[TRACER_PART_NUM_FLUID_QUANTITIES];
#endif
#ifdef TRACER_TRAJECTORY
  MyFloat rho, temp, utherm;
#ifdef TRACER_TRAJECTORY_EXTENDED_OUTPUT
  MyFloat composition[EOS_NSPECIES], dedt;
#endif
#endif
} data_out;

static data_out *DataResult, *DataOut;

static void out2particle(data_out * out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES)
    {
      tracer_nearest_dist[i] = out->nearestDist;
      P[i].Vel[0] = out->vel[0];
      P[i].Vel[1] = out->vel[1];
      P[i].Vel[2] = out->vel[2];
      tracer_newTimeBin[i] = out->newTimeBin;

#ifdef TRACER_PART_NUM_FLUID_QUANTITIES
      for(int k = 0; k < TRACER_PART_NUM_FLUID_QUANTITIES; k++)
        P[i].fluid_quantities[k] = out->fluid_quantities[k];
#endif
#ifdef TRACER_TRAJECTORY
      P[i].tRho = out->rho;
      P[i].tTemp = out->temp;
      P[i].tUtherm = out->utherm;
#ifdef TRACER_TRAJECTORY_EXTENDED_OUTPUT
      for(int k = 0; k < EOS_NSPECIES; k++)
        P[i].tComposition[k] = out->composition[k];
      P[i].tDedt = out->dedt;
#endif
#endif
    }
  else
    {
      if(out->nearestDist < tracer_nearest_dist[i])     /* closer cell on other task? */
        {
          tracer_nearest_dist[i] = out->nearestDist;
          P[i].Vel[0] = out->vel[0];
          P[i].Vel[1] = out->vel[1];
          P[i].Vel[2] = out->vel[2];
          tracer_newTimeBin[i] = out->newTimeBin;

#ifdef TRACER_TRAJECTORY
          P[i].tRho = out->rho;
          P[i].tTemp = out->temp;
          P[i].tUtherm = out->utherm;
#ifdef TRACER_TRAJECTORY_EXTENDED_OUTPUT
          for(int k = 0; k < EOS_NSPECIES; k++)
            P[i].tComposition[k] = out->composition[k];
          P[i].tDedt = out->dedt;
#endif
#endif


#ifdef TRACER_PART_NUM_FLUID_QUANTITIES
#if (TRACER_PART_TMAX)
          if(out->fluid_quantities[TracerPartTmaxIndex] > P[i].fluid_quantities[TracerPartTmaxIndex])
            {
              P[i].fluid_quantities[TracerPartTmaxIndex] = out->fluid_quantities[TracerPartTmaxIndex];
#if (TRACER_PART_TMAX_TIME)
              P[i].fluid_quantities[TracerPartTmaxTimeIndex] = All.Time;
#endif
#if (TRACER_PART_TMAX_RHO)
              P[i].fluid_quantities[TracerPartTmaxRhoIndex] = out->fluid_quantities[TracerPartTmaxRhoIndex];
#endif
            }
#endif
#if (TRACER_PART_RHOMAX)
          if(out->fluid_quantities[TracerPartRhomaxIndex] > P[i].fluid_quantities[TracerPartRhomaxIndex])
            {
              P[i].fluid_quantities[TracerPartRhomaxIndex] = out->fluid_quantities[TracerPartRhomaxIndex];
#if (TRACER_PART_RHOMAX_TIME)
              P[i].fluid_quantities[TracerPartRhomaxTimeIndex] = All.Time;
#endif
            }
#endif
#if (TRACER_PART_MACHMAX)
          if(out->fluid_quantities[TracerPartMachmaxIndex] > P[i].fluid_quantities[TracerPartMachmaxIndex])
            P[i].fluid_quantities[TracerPartMachmaxIndex] = out->fluid_quantities[TracerPartMachmaxIndex];
#endif
#if (TRACER_PART_ENTMAX)
          if(out->fluid_quantities[TracerPartEntmaxIndex] > P[i].fluid_quantities[TracerPartEntmaxIndex])
            {
              P[i].fluid_quantities[TracerPartEntmaxIndex] = out->fluid_quantities[TracerPartEntmaxIndex];
#if (TRACER_PART_ENTMAX_TIME)
              P[i].fluid_quantities[TracerPartEntmaxTimeIndex] = All.Time;
#endif
            }
#endif
#endif
        }
    }
};

#include "generic_comm_helpers2.h"


static void kernel_local(void)
{
  int idx;
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
              if(generic_polling_primary(count, TimeBinsTracer.NActiveParticles))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        idx = NextParticle++;

        if(idx >= TimeBinsTracer.NActiveParticles)
          break;

        int i = TimeBinsTracer.ActiveParticleList[idx];
        if(i < 0)
          continue;

        if(P[i].Type != TRACER_PARTICLE)
          continue;

        find_primary_cell_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
      }
  }
}

static void kernel_imported(void)
{
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

        find_primary_cell_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

void find_primary_cell(void)
{
  int idx, i;
  int ntot, iter = 0;

  generic_set_MaxNexport();

  /* we will repeat the whole thing for those points where we did not find a nearest neighbor */
  do
    {
      generic_comm_pattern(TimeBinsTracer.NActiveParticles, kernel_local, kernel_imported);

      int npleft = 0;

      /* do final operations on results */
      for(idx = 0; idx < TimeBinsTracer.NActiveParticles; idx++)
        {
          i = TimeBinsTracer.ActiveParticleList[idx];
          if(i < 0)
            continue;

          if(P[i].Type == TRACER_PARTICLE)
            {
              if(tracer_nearest_dist[i] == MAX_REAL_NUMBER)
                {
                  npleft++;
                  P[i].TracerHsml *= 2.0;
                  if(iter >= MAXITER - 10)
                    {
                      printf("i=%d task=%d hsml=%g nearest dist=%g pos=(%g|%g|%g)\n", i, ThisTask, P[i].TracerHsml, tracer_nearest_dist[i], P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
                      myflush(stdout);
                    }
                  if(iter > MAXITER)
                    terminate("TRACER_PARTICLE: iter > MAXITER");
                }
              else
                {
                  if(tracer_nearest_dist[i] > 0)
                    P[i].TracerHsml = 1.5 * dmax(tracer_nearest_dist[i], MinTracerHsml);
                  else
                    P[i].TracerHsml = 1.5 * MinTracerHsml;
                }
            }
        }

      /* sum up the left overs */
      MPI_Allreduce(&npleft, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      if(ntot > 0)              /* ok, we need to repeat for a few tracer particles */
        {
          iter++;
          if(iter > 0 && ThisTask == 0)
            {
              printf("TRACER_PARTICLE: iteration %d: need to repeat for %d points.\n", iter, ntot);
              myflush(stdout);
            }

          if(iter > MAXITER)
            terminate("TRACER_PARTICLE: failed to converge in tracer particles\n");
        }
    }
  while(ntot > 0);
}


void find_primary_cell_evaluate(int target, int mode, int threadid)
{
  int j, n, k, newTimeBin;
  double dt, dx, dy, dz, r, h, mu, nearestDist;
  MyDouble *pos;
  double vel[3], grad_vx[3], grad_vy[3], grad_vz[3];
#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif
#ifdef TRACER_PART_NUM_FLUID_QUANTITIES
  MyFloat fluid_quantities[TRACER_PART_NUM_FLUID_QUANTITIES];
#endif
#ifdef TRACER_TRAJECTORY
  MyFloat rho, temp, utherm;
#ifdef TRACER_TRAJECTORY_EXTENDED_OUTPUT
  MyFloat composition[EOS_NSPECIES], dedt;
#endif
#endif

  /* reset */
  nearestDist = MAX_REAL_NUMBER;
  newTimeBin = 0;
  vel[0] = vel[1] = vel[2] = 0;
#ifdef TRACER_PART_NUM_FLUID_QUANTITIES
  memset(fluid_quantities, 0, TRACER_PART_NUM_FLUID_QUANTITIES * sizeof(MyFloat));
#endif

#ifdef TRACER_TRAJECTORY
  rho = temp = utherm = 0.0;
#ifdef TRACER_TRAJECTORY_EXTENDED_OUTPUT
  for(k = 0; k < EOS_NSPECIES; k++)
    composition[k] = 0;
  dedt = 0;
#endif
#endif

  data_in local, *in;
  data_out out;
  int numnodes, *firstnode;

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

  pos = in->pos;
  h = in->hsml;

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, threadid, numnodes, firstnode);

  for(n = 0; n < nfound; n++)
    {
      j = Thread[threadid].Ngblist[n];

      dx = NEAREST_X(pos[0] - P[j].Pos[0]);
      dy = NEAREST_Y(pos[1] - P[j].Pos[1]);
      dz = NEAREST_Z(pos[2] - P[j].Pos[2]);

      r = sqrt(dx * dx + dy * dy + dz * dz);

      if(r < nearestDist && r < h && P[j].ID != 0 && P[j].Mass > 0)
        {
          nearestDist = r;
          /* updated stored fluid quantities from parent gas cell */
#ifdef COOLING
          mu = 4.0 / (1.0 + 3.0 * HYDROGEN_MASSFRAC + 4.0 * HYDROGEN_MASSFRAC * SphP[j].Ne);
#else
          mu = 1.0;
#endif

#ifdef TRACER_PART_NUM_FLUID_QUANTITIES
#if (TRACER_PART_TMAX)
#ifdef USE_SFR
          if(SphP[j].Sfr == 0)
            {
#endif
              fluid_quantities[TracerPartTmaxIndex] = GAMMA_MINUS1 / BOLTZMANN * SphP[j].Utherm * PROTONMASS * mu;
#if (TRACER_PART_TMAX_TIME)
              fluid_quantities[TracerPartTmaxTimeIndex] = All.Time;
#endif
#if (TRACER_PART_TMAX_RHO)
              fluid_quantities[TracerPartTmaxRhoIndex] = SphP[j].Density;
#endif
#ifdef USE_SFR
            }
#endif
#endif
#if (TRACER_PART_RHOMAX)
          fluid_quantities[TracerPartRhomaxIndex] = SphP[j].Density;
#if (TRACER_PART_RHOMAX_TIME)
          fluid_quantities[TracerPartRhomaxTimeIndex] = All.Time;
#endif
#endif
#if (TRACER_PART_MACHMAX)
          fluid_quantities[TracerPartMachmaxIndex] = SphP[j].MaxMach;
#endif
#if (TRACER_PART_ENTMAX)
          fluid_quantities[TracerPartEntmaxIndex] = (SphP[j].Pressure * All.cf_a3inv) / pow(SphP[j].Density * All.cf_a3inv, GAMMA);
#if (TRACER_PART_ENTMAX_TIME)
          fluid_quantities[TracerPartEntmaxTimeIndex] = All.Time;
#endif
#endif
#endif
#ifdef TRACER_TRAJECTORY
          rho = SphP[j].Density;
#if defined(EOS_DEGENERATE) || defined(EOS_OPAL)
          temp = SphP[j].EOSTemperature;
#endif
          utherm = SphP[j].Utherm;
#ifdef TRACER_TRAJECTORY_EXTENDED_OUTPUT
          for(k = 0; k < EOS_NSPECIES; k++)
            composition[k] = SphP[j].Composition[k];
          dedt = SphP[j].dedt;
#endif
#endif
          /* begin velocity calculation for tracers */
          for(k = 0; k < 3; k++)
            {
              grad_vx[k] = All.cf_atime * SphP[j].Grad.dvel[0][k];
              grad_vy[k] = All.cf_atime * SphP[j].Grad.dvel[1][k];
              grad_vz[k] = All.cf_atime * SphP[j].Grad.dvel[2][k];
            }

          if(All.ComovingIntegrationOn)
            {
              grad_vx[0] -= All.cf_atime * All.cf_atime * All.cf_hubble_a;
              grad_vy[1] -= All.cf_atime * All.cf_atime * All.cf_hubble_a;
              grad_vz[2] -= All.cf_atime * All.cf_atime * All.cf_hubble_a;
            }

          dt = (P[j].TimeBinHydro ? (((integertime) 1) << P[j].TimeBinHydro) : 0) * All.Timebase_interval;
          dt /= All.cf_hubble_a;

          /* spatial reconstruction */
          vel[0] = SphP[j].Momentum[0] / P[j].Mass + grad_vx[0] * dx + grad_vx[1] * dy + grad_vx[2] * dz;
          vel[1] = SphP[j].Momentum[1] / P[j].Mass + grad_vy[0] * dx + grad_vy[1] * dy + grad_vy[2] * dz;
          vel[2] = SphP[j].Momentum[2] / P[j].Mass + grad_vz[0] * dx + grad_vz[1] * dy + grad_vz[2] * dz;

          /* half step prediction */
          vel[0] -= 0.5 * dt * (SphP[j].Grad.dpress[0] / SphP[j].Density);
          vel[1] -= 0.5 * dt * (SphP[j].Grad.dpress[1] / SphP[j].Density);
          vel[2] -= 0.5 * dt * (SphP[j].Grad.dpress[2] / SphP[j].Density);

          newTimeBin = P[j].TimeBinHydro;
        }
    }

  out.nearestDist = nearestDist;
  for(k = 0; k < 3; k++)
    out.vel[k] = vel[k];

  out.newTimeBin = newTimeBin;

#ifdef TRACER_PART_NUM_FLUID_QUANTITIES
  for(k = 0; k < TRACER_PART_NUM_FLUID_QUANTITIES; k++)
    out.fluid_quantities[k] = fluid_quantities[k];
#endif

#ifdef TRACER_TRAJECTORY
  out.rho = rho;
  out.temp = temp;
  out.utherm = utherm;
#ifdef TRACER_TRAJECTORY_EXTENDED_OUTPUT
  for(k = 0; k < EOS_NSPECIES; k++)
    out.composition[k] = composition[k];
  out.dedt = dedt;
#endif
#endif

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;
}


/* put tracers into time integration */
void tracer_find_timesteps(void)
{
  int idx, i, bin, binold;

  for(idx = 0; idx < TimeBinsTracer.NActiveParticles; idx++)
    {
      i = TimeBinsTracer.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Type != TRACER_PARTICLE)
        continue;

      bin = tracer_newTimeBin[i];
      binold = P[i].TimeBinHydro;


      if(bin > binold)          /* timestep wants to increase */
        {
          while(TimeBinSynchronized[bin] == 0 && bin > binold)  /* make sure the new step is synchronized */
            bin--;
        }

      P[i].TimeBinHydro = bin;
      timebin_move_particle(&TimeBinsTracer, i, binold, bin);
    }
}


void tracer_particle_assign_cell_properties_and_timestep(void)
{
  int idx, i;

  MinTracerHsml = All.MinimumTracerHsml;

  CPU_Step[CPU_MISC] += measure_time();
  mpi_printf("TRACER_PARTICLE: updating tracers particles...\n");

  /* allocate */
  tracer_nearest_dist = (MyDouble *) mymalloc("tracer_nearest_dist", NumPart * sizeof(MyDouble));
  tracer_newTimeBin = (int *) mymalloc("tracer_newTimeBin", NumPart * sizeof(int));

  for(idx = 0; idx < TimeBinsTracer.NActiveParticles; idx++)
    {
      i = TimeBinsTracer.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Type == TRACER_PARTICLE)
        {
          P[i].Vel[0] = 0;
          P[i].Vel[1] = 0;
          P[i].Vel[2] = 0;
          tracer_nearest_dist[i] = MAX_REAL_NUMBER;
        }
    }

  /* find primary cell of tracer particles */
  find_primary_cell();

  /* assign time step to tracer particles */
  tracer_find_timesteps();

  /* free */
  myfree(tracer_newTimeBin);
  myfree(tracer_nearest_dist);

  mpi_printf("TRACER_PARTICLE: done.\n");
  CPU_Step[CPU_TRACERS] += measure_time();
}

#endif /* TRACER_PARTICLE */

#if defined(TRACER_PART_NUM_FLUID_QUANTITIES) || defined(TRACER_MC_NUM_FLUID_QUANTITIES)
/*! \brief Parse configuration options related to tracers.
 */
void set_tracer_part_indices(void)
{
  int j;

#ifdef TRACER_PART_NUM_FLUID_QUANTITIES
  j = 0;

  TracerPartTmaxIndex = TracerPartTmaxTimeIndex = TracerPartTmaxRhoIndex = TracerPartRhomaxIndex = TracerPartRhomaxTimeIndex =
    TracerPartMachmaxIndex = TracerPartEntmaxIndex = TracerPartEntmaxTimeIndex = -1;

  if(TRACER_PART_TMAX)
    TracerPartTmaxIndex = j++;

  if(TRACER_PART_TMAX_TIME)
    TracerPartTmaxTimeIndex = j++;

  if(TRACER_PART_TMAX_RHO)
    TracerPartTmaxRhoIndex = j++;

  if(TRACER_PART_RHOMAX)
    TracerPartRhomaxIndex = j++;

  if(TRACER_PART_RHOMAX_TIME)
    TracerPartRhomaxTimeIndex = j++;

  if(TRACER_PART_MACHMAX)
    TracerPartMachmaxIndex = j++;

  if(TRACER_PART_ENTMAX)
    TracerPartEntmaxIndex = j++;

  if(TRACER_PART_ENTMAX_TIME)
    TracerPartEntmaxTimeIndex = j++;

  if(j != TRACER_PART_NUM_FLUID_QUANTITIES)
    terminate("TRACER_PARTICLE: fluid quantities do not match -> expected %d got %d", TRACER_PART_NUM_FLUID_QUANTITIES, j);
  mpi_printf
    ("\nTRACER_PARTICLE: tracked quantity: \nTRACER_PARTICLE: TRACER_PART_TMAX=%d TRACER_PART_TMAX_TIME=%d TRACER_PART_TMAX_RHO=%d TRACER_PART_RHOMAX=%d TRACER_PART_RHOMAX_TIME=%d TRACER_PART_MACHMAX=%d TRACER_PART_ENTMAX=%d TRACER_PART_ENTMAX_TIME=%d\n",
     TRACER_PART_TMAX, TRACER_PART_TMAX_TIME, TRACER_PART_TMAX_RHO, TRACER_PART_RHOMAX, TRACER_PART_RHOMAX_TIME, TRACER_PART_MACHMAX, TRACER_PART_ENTMAX, TRACER_PART_ENTMAX_TIME);
  mpi_printf
    ("TRACER_PARTICLE: tracked quantity indices: \nTRACER_PARTICLE: TracerPartTmaxIndex=%d TracerPartTmaxTimeIndex=%d TracerPartTmaxRhoIndex=%d TracerPartRhomaxIndex=%d TracerPartRhomaxTimeIndex=%d TracerPartMachmaxIndex=%d TracerPartEntmaxIndex=%d TracerPartEntmaxTimeIndex=%d\n\n",
     TracerPartTmaxIndex, TracerPartTmaxTimeIndex, TracerPartTmaxRhoIndex, TracerPartRhomaxIndex, TracerPartRhomaxTimeIndex, TracerPartMachmaxIndex, TracerPartEntmaxIndex, TracerPartEntmaxTimeIndex);
#endif

#ifdef TRACER_MC_NUM_FLUID_QUANTITIES
  j = 0;

  TracerMCTmaxIndex = TracerMCTmaxTimeIndex = TracerMCTmaxRhoIndex = TracerMCRhomaxIndex = TracerMCRhomaxTimeIndex = TracerMCMachmaxIndex =
    TracerMCEntmaxIndex = TracerMCEntmaxTimeIndex = TracerMCLastStarTimeIndex = TracerMCWindCounterIndex = TracerMCExchangeCounterIndex =
    TracerMCExchangeDistanceIndex = TracerMCExchangeDistanceErrorIndex = TracerMCShockMachMaxIndex = -1;

  if(TRACER_MC_TMAX)
    TracerMCTmaxIndex = j++;

  if(TRACER_MC_TMAX_TIME)
    TracerMCTmaxTimeIndex = j++;

  if(TRACER_MC_TMAX_RHO)
    TracerMCTmaxRhoIndex = j++;

  if(TRACER_MC_RHOMAX)
    TracerMCRhomaxIndex = j++;

  if(TRACER_MC_RHOMAX_TIME)
    TracerMCRhomaxTimeIndex = j++;

  if(TRACER_MC_MACHMAX)
    TracerMCMachmaxIndex = j++;

  if(TRACER_MC_ENTMAX)
    TracerMCEntmaxIndex = j++;

  if(TRACER_MC_ENTMAX_TIME)
    TracerMCEntmaxTimeIndex = j++;

  if(TRACER_MC_LAST_STAR_TIME)
    TracerMCLastStarTimeIndex = j++;

  if(TRACER_MC_WIND_COUNTER)
    TracerMCWindCounterIndex = j++;

  if(TRACER_MC_EXCHANGE_COUNTER)
    TracerMCExchangeCounterIndex = j++;

  if(TRACER_MC_EXCHANGE_DISTANCE)
    TracerMCExchangeDistanceIndex = j++;

  if(TRACER_MC_EXCHANGE_DISTANCE_ERROR)
    TracerMCExchangeDistanceErrorIndex = j++;

  if(TRACER_MC_SHOCKMACHNUM_MAX)
    TracerMCShockMachMaxIndex = j++;

  if(j != TRACER_MC_NUM_FLUID_QUANTITIES)
    terminate("TRACER_MC: fluid quantities do not match -> expected %d got %d", TRACER_MC_NUM_FLUID_QUANTITIES, j);
  mpi_printf
    ("\nTRACER_MC: tracked quantity: \nTRACER_MC: TRACER_MC_TMAX=%d TRACER_MC_TMAX_TIME=%d TRACER_MC_TMAX_RHO=%d TRACER_MC_RHOMAX=%d TRACER_MC_RHOMAX_TIME=%d TRACER_MC_MACHMAX=%d TRACER_MC_ENTMAX=%d TRACER_MC_ENTMAX_TIME=%d TRACER_MC_LAST_STAR_TIME=%d TRACER_MC_WIND_COUNTER=%d TRACER_MC_EXCHANGE_COUNTER=%d TRACER_MC_EXCHANGE_DISTANCE=%d TRACER_MC_EXCHANGE_DISTANCE_ERROR=%d TRACER_MC_SHOCKMACHNUM_MAX=%d\n",
     TRACER_MC_TMAX, TRACER_MC_TMAX_TIME, TRACER_MC_TMAX_RHO, TRACER_MC_RHOMAX, TRACER_MC_RHOMAX_TIME, TRACER_MC_MACHMAX, TRACER_MC_ENTMAX,
     TRACER_MC_ENTMAX_TIME, TRACER_MC_LAST_STAR_TIME, TRACER_MC_WIND_COUNTER, TRACER_MC_EXCHANGE_COUNTER, TRACER_MC_EXCHANGE_DISTANCE, TRACER_MC_EXCHANGE_DISTANCE_ERROR, TRACER_MC_SHOCKMACHNUM_MAX);
  mpi_printf
    ("TRACER_MC: tracked quantity indices: \nTRACER_MC: TracerMCTmaxIndex=%d TracerMCTmaxTimeIndex=%d TracerMCTmaxRhoIndex=%d TracerMCRhomaxIndex=%d TracerMCRhomaxTimeIndex=%d TracerMCMachmaxIndex=%d TracerMCEntmaxIndex=%d TracerMCEntmaxTimeIndex=%d TracerMCLastStarTimeIndex=%d TracerMCWindCounterIndex=%d TracerMCExchangeCounterIndex=%d TracerMCExchangeDistanceIndex=%d TracerMCExchangeDistanceErrorIndex=%d TracerMCShockMachMaxIndex=%d\n\n",
     TracerMCTmaxIndex, TracerMCTmaxTimeIndex, TracerMCTmaxRhoIndex, TracerMCRhomaxIndex, TracerMCRhomaxTimeIndex, TracerMCMachmaxIndex,
     TracerMCEntmaxIndex, TracerMCEntmaxTimeIndex, TracerMCLastStarTimeIndex, TracerMCWindCounterIndex, TracerMCExchangeCounterIndex,
     TracerMCExchangeDistanceIndex, TracerMCExchangeDistanceErrorIndex, TracerMCShockMachMaxIndex);
#endif

}
#endif /* TRACER_PART_NUM_FLUID_QUANTITIES || TRACER_MC_NUM_FLUID_QUANTITIES */

#if defined(TRACER_MC) || defined(TRACER_PARTICLE)

/*! zero all stored fluid quantities after they are written, such that the "maximum" values are 
 *  actually local maxima between each successive snapshot and the last
 */
void reset_tracer_parent_fluid_properties(void)
{
#if defined(TRACER_MC_NUM_FLUID_QUANTITIES) || defined(TRACER_PART_NUM_FLUID_QUANTITIES)
  int i, k;
#endif

#if defined(TRACER_MC) && defined(TRACER_MC_NUM_FLUID_QUANTITIES)
  for(i = 0; i < N_tracer; i++)
    for(k = 0; k < TRACER_MC_NUM_FLUID_QUANTITIES; k++)
      {
#if (TRACER_MC_LAST_STAR_TIME)
        if(k == TracerMCLastStarTimeIndex)
          continue;
#endif
#if (TRACER_MC_WIND_COUNTER)
        if(k == TracerMCWindCounterIndex)
          continue;
#endif
#if (TRACER_MC_EXCHANGE_COUNTER)
        if(k == TracerMCExchangeCounterIndex)
          continue;
#endif
#if (TRACER_MC_EXCHANGE_DISTANCE)
        if(k == TracerMCExchangeDistanceIndex)
          continue;
#endif
#if (TRACER_MC_EXCHANGE_DISTANCE_ERROR)
        if(k == TracerMCExchangeDistanceErrorIndex)
          continue;
#endif

        TracerLinkedList[i].fluid_quantities[k] = 0.0;
      }
#endif

#if defined(TRACER_PARTICLE) && defined(TRACER_PART_NUM_FLUID_QUANTITIES)
  for(i = 0; i < NumPart; i++)
    for(k = 0; k < TRACER_PART_NUM_FLUID_QUANTITIES; k++)
      P[i].fluid_quantities[k] = 0.0;
#endif
}

#endif /* TRACER_MC || TRACER_PARTICLE */
