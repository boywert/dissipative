/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/predict.c
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
#include "allvars.h"
#include "proto.h"

/*! \file predict.c
 *  \brief Contains the routines to find the next sync point,
 * manage the list of active timebins/active particles and to drift particles.
 */

/*! \brief This function (re)builds the time bin lists.
 *
 * It counts the number of particles in each timebin and updates the
 * linked lists containing the particles of each time bin. Afterwards the
 * linked list of active particles is updated by make_list_of_active_particles().
 *
 * The linked lists for each timebin are stored in #FirstInTimeBin[], #LastInTimeBin[],
 * #PrevInTimeBin[] and #NextInTimeBin[]. The counters of particles per timebin are
 * #TimeBinCount and #TimeBinCountSph.
 */
void reconstruct_timebins(void)
{
  TIMER_START(CPU_TIMELINE);

  int i, bin;

  for(bin = 0; bin < TIMEBINS; bin++)
    {
      TimeBinsHydro.TimeBinCount[bin] = 0;
      TimeBinsHydro.FirstInTimeBin[bin] = -1;
      TimeBinsHydro.LastInTimeBin[bin] = -1;

      TimeBinsGravity.TimeBinCount[bin] = 0;
      TimeBinsGravity.FirstInTimeBin[bin] = -1;
      TimeBinsGravity.LastInTimeBin[bin] = -1;
#ifdef USE_SFR
      TimeBinSfr[bin] = 0;
#endif
#ifdef BLACK_HOLES
      TimeBinsBHAccretion.TimeBinCount[bin] = 0;
      TimeBinsBHAccretion.FirstInTimeBin[bin] = -1;
      TimeBinsBHAccretion.LastInTimeBin[bin] = -1;

      TimeBin_BH_mass[bin] = 0;
      TimeBin_BH_dynamicalmass[bin] = 0;
      TimeBin_BH_Mdot[bin] = 0;
      TimeBin_BH_Medd[bin] = 0;
#endif
#ifdef TRACER_PARTICLE
      TimeBinsTracer.TimeBinCount[bin] = 0;
      TimeBinsTracer.FirstInTimeBin[bin] = -1;
      TimeBinsTracer.LastInTimeBin[bin] = -1;
#endif
#ifdef SINKS
      TimeBinsSinksAccretion.TimeBinCount[bin] = 0;
      TimeBinsSinksAccretion.FirstInTimeBin[bin] = -1;
      TimeBinsSinksAccretion.LastInTimeBin[bin] = -1;
#endif
#ifdef DUST_LIVE
      TimeBinsDust.TimeBinCount[bin] = 0;
      TimeBinsDust.FirstInTimeBin[bin] = -1;
      TimeBinsDust.LastInTimeBin[bin] = -1;
#endif
    }

  for(i = 0; i < NumGas; i++)
    {
      if(P[i].ID == 0 && P[i].Mass == 0)
        continue;

      if(P[i].Type != 0)
        continue;

      bin = P[i].TimeBinHydro;

      if(TimeBinsHydro.TimeBinCount[bin] > 0)
        {
          TimeBinsHydro.PrevInTimeBin[i] = TimeBinsHydro.LastInTimeBin[bin];
          TimeBinsHydro.NextInTimeBin[i] = -1;
          TimeBinsHydro.NextInTimeBin[TimeBinsHydro.LastInTimeBin[bin]] = i;
          TimeBinsHydro.LastInTimeBin[bin] = i;
        }
      else
        {
          TimeBinsHydro.FirstInTimeBin[bin] = TimeBinsHydro.LastInTimeBin[bin] = i;
          TimeBinsHydro.PrevInTimeBin[i] = TimeBinsHydro.NextInTimeBin[i] = -1;
        }
      TimeBinsHydro.TimeBinCount[bin]++;

#ifdef USE_SFR
      TimeBinSfr[bin] += SphP[i].Sfr;
#endif
    }

  for(i = 0; i < NumPart; i++)
    {
      if(P[i].ID == 0 && P[i].Mass == 0)
        continue;

      bin = P[i].TimeBinGrav;

#ifdef TRACER_PARTICLE
      if(P[i].Type == TRACER_PARTICLE)
        continue;
#endif

      if(TimeBinsGravity.TimeBinCount[bin] > 0)
        {
          TimeBinsGravity.PrevInTimeBin[i] = TimeBinsGravity.LastInTimeBin[bin];
          TimeBinsGravity.NextInTimeBin[i] = -1;
          TimeBinsGravity.NextInTimeBin[TimeBinsGravity.LastInTimeBin[bin]] = i;
          TimeBinsGravity.LastInTimeBin[bin] = i;
        }
      else
        {
          TimeBinsGravity.FirstInTimeBin[bin] = TimeBinsGravity.LastInTimeBin[bin] = i;
          TimeBinsGravity.PrevInTimeBin[i] = TimeBinsGravity.NextInTimeBin[i] = -1;
        }
      TimeBinsGravity.TimeBinCount[bin]++;

#ifdef BLACK_HOLES
      if(P[i].Type == 5)
        {
          bin = P[i].TimeBinHydro;

          if(TimeBinsBHAccretion.TimeBinCount[bin] > 0)
            {
              TimeBinsBHAccretion.PrevInTimeBin[i] = TimeBinsBHAccretion.LastInTimeBin[bin];
              TimeBinsBHAccretion.NextInTimeBin[i] = -1;
              TimeBinsBHAccretion.NextInTimeBin[TimeBinsBHAccretion.LastInTimeBin[bin]] = i;
              TimeBinsBHAccretion.LastInTimeBin[bin] = i;
            }
          else
            {
              TimeBinsBHAccretion.FirstInTimeBin[bin] = TimeBinsBHAccretion.LastInTimeBin[bin] = i;
              TimeBinsBHAccretion.PrevInTimeBin[i] = TimeBinsBHAccretion.NextInTimeBin[i] = -1;
            }
          TimeBinsBHAccretion.TimeBinCount[bin]++;


          TimeBin_BH_mass[bin] += BPP(i).BH_Mass;
          TimeBin_BH_dynamicalmass[bin] += P[i].Mass;
          TimeBin_BH_Mdot[bin] += BPP(i).BH_Mdot;
          if(BPP(i).BH_Mass > 0)
            TimeBin_BH_Medd[bin] += BPP(i).BH_Mdot / BPP(i).BH_Mass;
        }
#endif

#ifdef SINKS
      if(P[i].Type == 5)
        {
          bin = P[i].TimeBinHydro;

          if(TimeBinsSinksAccretion.TimeBinCount[bin] > 0)
            {
              TimeBinsSinksAccretion.PrevInTimeBin[i] = TimeBinsSinksAccretion.LastInTimeBin[bin];
              TimeBinsSinksAccretion.NextInTimeBin[i] = -1;
              TimeBinsSinksAccretion.NextInTimeBin[TimeBinsSinksAccretion.LastInTimeBin[bin]] = i;
              TimeBinsSinksAccretion.LastInTimeBin[bin] = i;
            }
          else
            {
              TimeBinsSinksAccretion.FirstInTimeBin[bin] = TimeBinsSinksAccretion.LastInTimeBin[bin] = i;
              TimeBinsSinksAccretion.PrevInTimeBin[i] = TimeBinsSinksAccretion.NextInTimeBin[i] = -1;
            }
          TimeBinsSinksAccretion.TimeBinCount[bin]++;
        }
#endif
    }

#ifdef TRACER_PARTICLE
  for(i = 0; i < NumPart; i++)
    {
      if(P[i].ID == 0 && P[i].Mass == 0)
        continue;

      if(P[i].Type != TRACER_PARTICLE)
        continue;

      bin = P[i].TimeBinHydro;

      if(TimeBinsTracer.TimeBinCount[bin] > 0)
        {
          TimeBinsTracer.PrevInTimeBin[i] = TimeBinsTracer.LastInTimeBin[bin];
          TimeBinsTracer.NextInTimeBin[i] = -1;
          TimeBinsTracer.NextInTimeBin[TimeBinsTracer.LastInTimeBin[bin]] = i;
          TimeBinsTracer.LastInTimeBin[bin] = i;
        }
      else
        {
          TimeBinsTracer.FirstInTimeBin[bin] = TimeBinsTracer.LastInTimeBin[bin] = i;
          TimeBinsTracer.PrevInTimeBin[i] = TimeBinsTracer.NextInTimeBin[i] = -1;
        }
      TimeBinsTracer.TimeBinCount[bin]++;
    }
#endif

#ifdef DUST_LIVE
  for(i = 0; i < NumPart; i++)
    {
      if(P[i].ID == 0 && P[i].Mass == 0)
        continue;

      if(P[i].Type != DUST_LIVE)
        continue;

      bin = P[i].TimeBinHydro;

      if(TimeBinsDust.TimeBinCount[bin] > 0)
        {
          TimeBinsDust.PrevInTimeBin[i] = TimeBinsDust.LastInTimeBin[bin];
          TimeBinsDust.NextInTimeBin[i] = -1;
          TimeBinsDust.NextInTimeBin[TimeBinsDust.LastInTimeBin[bin]] = i;
          TimeBinsDust.LastInTimeBin[bin] = i;
        }
      else
        {
          TimeBinsDust.FirstInTimeBin[bin] = TimeBinsDust.LastInTimeBin[bin] = i;
          TimeBinsDust.PrevInTimeBin[i] = TimeBinsDust.NextInTimeBin[i] = -1;
        }
      TimeBinsDust.TimeBinCount[bin]++;
    }
#endif

  make_list_of_active_particles();

  TIMER_STOP(CPU_TIMELINE);
}


/*! \brief This function finds the next synchronization point of the system.
 * (i.e. the earliest point of time any of the particles needs a force
 * computation), and drifts the system to this point of time.
 *
 * This function drifts all particles, including inactive particles to the
 * next sync point. This is done by drift_particles(). Afterwards the linked
 * list of active particles is updated to the new sync point by
 * make_list_of_active_particles(). Particles become active/inactive here.
 */
void find_next_sync_point(void)
{
  int n;
  integertime ti_next_kick, ti_next_kick_global, ti_next_for_bin, dt_bin;
  double timeold;

  TIMER_START(CPU_DRIFTS);

  timeold = All.Time;

  All.NumCurrentTiStep++;

  /* find the next kick time */
  ti_next_kick = TIMEBASE;

  for(n = 0; n < TIMEBINS; n++)
    {
      int active = TimeBinsHydro.TimeBinCount[n];

#if (defined(SELFGRAVITY) || defined(EXTERNALGRAVITY) || defined(EXACT_GRAVITY_FOR_PARTICLE_TYPE)) && !defined(MESHRELAX)
      active += TimeBinsGravity.TimeBinCount[n];
#endif
#ifdef TRACER_PARTICLE
      active += TimeBinsTracer.TimeBinCount[n];
#endif
#ifdef BLACK_HOLES
      active += TimeBinsBHAccretion.TimeBinCount[n];
#endif
#ifdef SINKS
      active += TimeBinsSinksAccretion.TimeBinCount[n];
#endif
#ifdef DUST_LIVE
      active += TimeBinsDust.TimeBinCount[n];
#endif
      if(active)
        {
          if(n > 0)
            {
              dt_bin = (((integertime) 1) << n);
              ti_next_for_bin = (All.Ti_Current / dt_bin) * dt_bin + dt_bin;    /* next kick time for this timebin */
            }
          else
            {
              dt_bin = 0;
              ti_next_for_bin = All.Ti_Current;
            }

          if(ti_next_for_bin < ti_next_kick)
            ti_next_kick = ti_next_for_bin;
        }

    }

#ifdef ENLARGE_DYNAMIC_RANGE_IN_TIME
  minimum_large_ints(1, &ti_next_kick, &ti_next_kick_global);
#else
  MPI_Allreduce(&ti_next_kick, &ti_next_kick_global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
#endif

  All.Previous_Ti_Current = All.Ti_Current;
  All.Ti_Current = ti_next_kick_global;

  if(All.ComovingIntegrationOn)
    All.Time = All.TimeBegin * exp(All.Ti_Current * All.Timebase_interval);
  else
    All.Time = All.TimeBegin + All.Ti_Current * All.Timebase_interval;

  set_cosmo_factors_for_current_time();
  All.TimeStep = All.Time - timeold;

  mark_active_timebins();

  TIMER_STOP(CPU_DRIFTS);
}


void mark_active_timebins(void)
{
  int n;
  int lowest_active_bin = TIMEBINS, highest_active_bin = 0;
  int lowest_occupied_bin = TIMEBINS, highest_occupied_bin = 0;
  int lowest_occupied_gravity_bin = TIMEBINS, highest_occupied_gravity_bin = 0;
  int highest_synchronized_bin = 0;
  int nsynchronized_gravity = 0, nsynchronized_hydro = 0;
  integertime dt_bin;

  /* mark the bins that will be synchronized/active */

  for(n = 0; n < TIMEBINS; n++)
    {
      if(TimeBinsGravity.TimeBinCount[n])
        {
          if(highest_occupied_gravity_bin < n)
            highest_occupied_gravity_bin = n;

          if(lowest_occupied_gravity_bin > n)
            lowest_occupied_gravity_bin = n;
        }

      int active = TimeBinsHydro.TimeBinCount[n] + TimeBinsGravity.TimeBinCount[n];
#ifdef TRACER_PARTICLE
      active += TimeBinsTracer.TimeBinCount[n];
#endif
#ifdef BLACK_HOLES
      active += TimeBinsBHAccretion.TimeBinCount[n];
#endif
#ifdef SINKS
      active += TimeBinsSinksAccretion.TimeBinCount[n];
#endif
#ifdef DUST_LIVE
      active += TimeBinsDust.TimeBinCount[n];
#endif
      if(active)
        {
          if(highest_occupied_bin < n)
            highest_occupied_bin = n;

          if(lowest_occupied_bin > n)
            lowest_occupied_bin = n;
        }

      dt_bin = (((integertime) 1) << n);

      if((All.Ti_Current % dt_bin) == 0)
        {
          TimeBinSynchronized[n] = 1;
          All.Ti_begstep[n] = All.Ti_Current;

          nsynchronized_gravity += TimeBinsGravity.TimeBinCount[n];
          nsynchronized_hydro += TimeBinsHydro.TimeBinCount[n];

          if(highest_synchronized_bin < n)
            highest_synchronized_bin = n;

          if(active)
            {
              if(highest_active_bin < n)
                highest_active_bin = n;

              if(lowest_active_bin > n)
                lowest_active_bin = n;
            }
        }
      else
        TimeBinSynchronized[n] = 0;
    }



  int lowest_in[3], lowest_out[3];
  lowest_in[0] = lowest_occupied_bin;
  lowest_in[1] = lowest_occupied_gravity_bin;
  lowest_in[2] = lowest_active_bin;
  MPI_Allreduce(lowest_in, lowest_out, 3, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  All.LowestOccupiedTimeBin = lowest_out[0];
  All.LowestOccupiedGravTimeBin = lowest_out[1];
  All.LowestActiveTimeBin = lowest_out[2];



  int highest_in[4], highest_out[4];
  highest_in[0] = highest_occupied_bin;
  highest_in[1] = highest_occupied_gravity_bin;
  highest_in[2] = highest_active_bin;
  highest_in[3] = highest_synchronized_bin;
  MPI_Allreduce(highest_in, highest_out, 4, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  All.HighestOccupiedTimeBin = highest_out[0];
  All.HighestOccupiedGravTimeBin = highest_out[1];
  All.HighestActiveTimeBin = highest_out[2];
  All.HighestSynchronizedTimeBin = highest_out[3];


  /* note: the lowest synchronized bin is always 1 */


  int input_ints[2 + 2 * TIMEBINS];
  long long output_longs[2 + 2 * TIMEBINS];

  input_ints[0] = nsynchronized_hydro;
  input_ints[1] = nsynchronized_gravity;
  memcpy(input_ints + 2, TimeBinsGravity.TimeBinCount, TIMEBINS * sizeof(int));
  memcpy(input_ints + 2 + TIMEBINS, TimeBinsHydro.TimeBinCount, TIMEBINS * sizeof(int));

  sumup_large_ints(2 + 2 * TIMEBINS, input_ints, output_longs);

  All.GlobalNSynchronizedHydro = output_longs[0];
  All.GlobalNSynchronizedGravity = output_longs[1];
  long long *tot_count_grav = output_longs + 2;
  long long *tot_count_sph = output_longs + 2 + TIMEBINS;


  long long tot_grav = 0, tot_sph = 0;

  for(n = 0; n < TIMEBINS; n++)
    {
      tot_grav += tot_count_grav[n];
      tot_sph += tot_count_sph[n];

      if(n > 0)
        {
          tot_count_grav[n] += tot_count_grav[n - 1];
          tot_count_sph[n] += tot_count_sph[n - 1];
        }
    }

  All.SmallestTimeBinWithDomainDecomposition = All.HighestOccupiedTimeBin;

  for(n = All.HighestOccupiedTimeBin; n >= All.LowestOccupiedTimeBin; n--)
    {
      if(tot_count_grav[n] > All.ActivePartFracForNewDomainDecomp * tot_grav || tot_count_sph[n] > All.ActivePartFracForNewDomainDecomp * tot_sph)
        All.SmallestTimeBinWithDomainDecomposition = n;
    }



#ifdef GFM_AGN_RADIATION
  All.SmallestTimeBinWithAGNRad = All.LowestActiveTimeBin + GFM_MAX_TIMEBINS_WITHOUT_AGN_RAD;
  if(All.SmallestTimeBinWithAGNRad > All.HighestOccupiedTimeBin)
    All.SmallestTimeBinWithAGNRad = All.HighestOccupiedTimeBin;
#endif

#ifdef DMPIC
  integertime ti_next_pic = dmpic_get_next_pictime();
  while(All.PicCurrentNr < All.PicCount && ti_next_pic >= All.Previous_Ti_Current && ti_next_pic < All.Ti_Current)
    {
      drift_particles(All.Previous_Ti_Current, ti_next_pic);
      dmpic_make_image();
      All.Previous_Ti_Current = ti_next_pic;
      ti_next_pic = dmpic_get_next_pictime();
    }
#endif
}



void drift_all_particles(void)
{
  int i;

  TIMER_START(CPU_DRIFTS);

#pragma omp parallel for private(i)
  for(i = 0; i < NumPart; i++)
    drift_particle(i, All.Ti_Current);

  TIMER_STOP(CPU_DRIFTS);
}




/*! \brief This function drifts drifts a particle i to time1
 *
 * @param i     particle/cell index
 * @param time1 time to which particles get drifted
 */
void drift_particle(int i, integertime time1)
{
  int j;

  if(i < 0)
    terminate("i=%d  NumPart=%d", i, NumPart);

  integertime time0 = P[i].Ti_Current;

  if(time1 == time0)
    return;

  if(time1 < time0)
    terminate("no prediction into past allowed: time0=%lld time1=%lld\n", (long long) time0, (long long) time1);

  double dt_drift;

  if(All.ComovingIntegrationOn)
    dt_drift = get_drift_factor(time0, time1);
  else
    dt_drift = (time1 - time0) * All.Timebase_interval;

  if(P[i].Type == 0)
    {
      for(j = 0; j < 3; j++)
        {
          P[i].Pos[j] += SphP[i].VelVertex[j] * dt_drift;
        }
    }
  else
    {
#ifndef MESHRELAX
      for(j = 0; j < 3; j++)
        P[i].Pos[j] += P[i].Vel[j] * dt_drift;

#if defined(REFLECTIVE_X)
      if(P[i].Pos[0] < 0 || P[i].Pos[0] > boxSize_X)
        {
          P[i].Pos[0] = 2 * (P[i].Pos[0] > boxSize_X ? 1 : 0) * boxSize_X - P[i].Pos[0];
          P[i].Vel[0] *= -1;
        }
#endif
#if defined(REFLECTIVE_Y)
      if(P[i].Pos[1] < 0 || P[i].Pos[1] > boxSize_Y)
        {
          P[i].Pos[1] = 2 * (P[i].Pos[1] > boxSize_Y ? 1 : 0) * boxSize_Y - P[i].Pos[1];
          P[i].Vel[1] *= -1;
        }
#endif
#if defined(REFLECTIVE_Z)
      if(P[i].Pos[2] < 0 || P[i].Pos[2] > boxSize_Z)
        {
          P[i].Pos[2] = 2 * (P[i].Pos[2] > boxSize_Z ? 1 : 0) * boxSize_Z - P[i].Pos[2];
          P[i].Vel[2] *= -1;
        }
#endif

#endif
    }

#ifdef RELAXOBJECT_BINARY
  double vel;
  if(P[i].Pos[0] < 0.5 * All.BoxSize)
    vel = 1e7 * 0.9 / 1.1;
  else
    vel = -1e7;

  P[i].Pos[0] += vel * dt_drift;
#endif

  P[i].Ti_Current = time1;
}


static int int_compare(const void *a, const void *b)
{
  if(*((int *) a) < *((int *) b))
    return -1;

  if(*((int *) a) > *((int *) b))
    return +1;

  return 0;
}

/*! \brief This function builds the linear list of active particles.
 *
 * The list is stored in the array ActiveParticleList of the
 * TimeBinData structs
 */
void make_list_of_active_particles(void)
{
  TIMER_START(CPU_DRIFTS);

  int i, n;
  /* make a link list with the particles in the active time bins */
  TimeBinsHydro.NActiveParticles = 0;

  for(n = 0; n < TIMEBINS; n++)
    {
      if(TimeBinSynchronized[n])
        {
          for(i = TimeBinsHydro.FirstInTimeBin[n]; i >= 0; i = TimeBinsHydro.NextInTimeBin[i])
            if((P[i].Type == 0) && !((P[i].ID == 0) && (P[i].Mass == 0)))
              {
                if(P[i].Ti_Current != All.Ti_Current)
                  drift_particle(i, All.Ti_Current);

                TimeBinsHydro.ActiveParticleList[TimeBinsHydro.NActiveParticles] = i;
                TimeBinsHydro.NActiveParticles++;
              }
        }
    }

  TimeBinsGravity.NActiveParticles = 0;

  for(n = 0; n < TIMEBINS; n++)
    {
      if(TimeBinSynchronized[n])
        {
          for(i = TimeBinsGravity.FirstInTimeBin[n]; i >= 0; i = TimeBinsGravity.NextInTimeBin[i])
            {
              if(!((P[i].ID == 0) && (P[i].Mass == 0)))
                {
                  if(P[i].Ti_Current != All.Ti_Current)
                    drift_particle(i, All.Ti_Current);

                  TimeBinsGravity.ActiveParticleList[TimeBinsGravity.NActiveParticles] = i;
                  TimeBinsGravity.NActiveParticles++;
                }
            }
        }
    }

  /* sort both lists for better memory efficiency */
  mysort(TimeBinsHydro.ActiveParticleList, TimeBinsHydro.NActiveParticles, sizeof(int), int_compare);
  mysort(TimeBinsGravity.ActiveParticleList, TimeBinsGravity.NActiveParticles, sizeof(int), int_compare);

#ifdef TRACER_PARTICLE
  TimeBinsTracer.NActiveParticles = 0;

  for(n = 0; n < TIMEBINS; n++)
    {
      if(TimeBinSynchronized[n])
        {
          for(i = TimeBinsTracer.FirstInTimeBin[n]; i >= 0; i = TimeBinsTracer.NextInTimeBin[i])
            if(P[i].Type == TRACER_PARTICLE)
              {
                if(P[i].Ti_Current != All.Ti_Current)
                  drift_particle(i, All.Ti_Current);

                TimeBinsTracer.ActiveParticleList[TimeBinsTracer.NActiveParticles] = i;
                TimeBinsTracer.NActiveParticles++;
              }
        }
    }

  mysort(TimeBinsTracer.ActiveParticleList, TimeBinsTracer.NActiveParticles, sizeof(int), int_compare);
#endif

#ifdef BLACK_HOLES
  TimeBinsBHAccretion.NActiveParticles = 0;

  for(n = 0; n < TIMEBINS; n++)
    {
      if(TimeBinSynchronized[n])
        {
          for(i = TimeBinsBHAccretion.FirstInTimeBin[n]; i >= 0; i = TimeBinsBHAccretion.NextInTimeBin[i])
            if((P[i].Type == 5) && (P[i].ID != 0) && (P[i].Mass != 0))
              {
                if(P[i].Ti_Current != All.Ti_Current)
                  drift_particle(i, All.Ti_Current);

                TimeBinsBHAccretion.ActiveParticleList[TimeBinsBHAccretion.NActiveParticles] = i;
                TimeBinsBHAccretion.NActiveParticles++;
              }
        }
    }

  mysort(TimeBinsBHAccretion.ActiveParticleList, TimeBinsBHAccretion.NActiveParticles, sizeof(int), int_compare);
#endif

#ifdef SINKS
  TimeBinsSinksAccretion.NActiveParticles = 0;

  for(n = 0; n < TIMEBINS; n++)
    {
      if(TimeBinSynchronized[n])
        {
          for(i = TimeBinsSinksAccretion.FirstInTimeBin[n]; i >= 0; i = TimeBinsSinksAccretion.NextInTimeBin[i])
            if((P[i].Type == 5) && !((P[i].ID == 0) && (P[i].Mass == 0)))
              {
                if(P[i].Ti_Current != All.Ti_Current)
                  drift_particle(i, All.Ti_Current);

                TimeBinsSinksAccretion.ActiveParticleList[TimeBinsSinksAccretion.NActiveParticles] = i;
                TimeBinsSinksAccretion.NActiveParticles++;
              }
        }
    }

  mysort(TimeBinsSinksAccretion.ActiveParticleList, TimeBinsSinksAccretion.NActiveParticles, sizeof(int), int_compare);
#endif

#ifdef DUST_LIVE
  TimeBinsDust.NActiveParticles = 0;

  for(n = 0; n < TIMEBINS; n++)
    {
      if(TimeBinSynchronized[n])
        {
          for(i = TimeBinsDust.FirstInTimeBin[n]; i >= 0; i = TimeBinsDust.NextInTimeBin[i])
            if((P[i].Type == DUST_LIVE) && (P[i].ID != 0) && (P[i].Mass != 0))
              {
                if(P[i].Ti_Current != All.Ti_Current)
                  drift_particle(i, All.Ti_Current);

                TimeBinsDust.ActiveParticleList[TimeBinsDust.NActiveParticles] = i;
                TimeBinsDust.NActiveParticles++;
              }
        }
    }

  mysort(TimeBinsDust.ActiveParticleList, TimeBinsDust.NActiveParticles, sizeof(int), int_compare);
#endif

  int in[6];
  long long out[6];

  n = 2;
  in[0] = TimeBinsGravity.NActiveParticles;
  in[1] = TimeBinsHydro.NActiveParticles;

#ifdef BLACK_HOLES
  in[2] = TimeBinsBHAccretion.NActiveParticles;
  n = 3;
#endif
#ifdef TRACER_PARTICLE
  in[3] = TimeBinsTracer.NActiveParticles;
  n = 4;
#endif
#ifdef DUST_LIVE
  in[4] = TimeBinsDust.NActiveParticles;
  n = 5;
#endif
#ifdef SINKS
  in[5] = TimeBinsSinksAccretion.NActiveParticles;
  n = 6;
#endif

  sumup_large_ints(n, in, out);

  TimeBinsGravity.GlobalNActiveParticles = out[0];
  TimeBinsHydro.GlobalNActiveParticles = out[1];

#ifdef BLACK_HOLES
  TimeBinsBHAccretion.GlobalNActiveParticles = out[2];
#endif
#ifdef TRACER_PARTICLE
  TimeBinsTracer.GlobalNActiveParticles = out[3];
#endif
#ifdef DUST_LIVE
  TimeBinsDust.GlobalNActiveParticles = out[4];
#endif
#ifdef SINKS
  TimeBinsSinksAccretion.GlobalNActiveParticles = out[5];
#endif

  TIMER_STOP(CPU_DRIFTS);
}
