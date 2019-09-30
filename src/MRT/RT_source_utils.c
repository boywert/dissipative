/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/MRT/RT_source_utils.c
 * \date        09/2017
 * \author      Federico Marinacci, Rahul Kannan, David Barnes
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - 10.10.2017 Black hole routines added
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../proto.h"
#include "../allvars.h"

#ifdef MRT_SOURCES

#define PHOT_NORMALIZATION 1e63

#ifdef MRT_STARS

#define LIFETIME 0.1

double get_photon_released_star(int p, double dt)
{
  if(P[p].Type != 4)
    terminate("Particle should be a star! Type=%d, index %d", P[p].Type, p);

  double phot = 5e48 * dt * 1e3 * SEC_PER_MEGAYEAR / PHOT_NORMALIZATION;

  return phot;
}

void start_stellar_sources(void)
{
  int idx, i;
  double time_begstep, dt, dtime, dtime_in_Gyr, age_of_star_in_Gyr_endstep, age_of_star_in_Gyr;

  //TIMER_START(CPU_GFM_ENRICH);

  //  mpi_printf("GFM_STELLAR_EVOLUTION: GFM_STELLAR_EVOLUTION...\n");

  StarParticle = mymalloc("StarParticle", N_star * sizeof(struct star_particle));

  Nsource = 0;

  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Ti_Current != All.Ti_Current)
        {
          terminate("how can this be?");
          drift_particle(i, All.Ti_Current);
        }

      if(P[i].Type == 4 && P[i].Mass > 0 && STP(i).BirthTime > 0)
        {
          if(All.ComovingIntegrationOn)
            time_begstep = All.TimeBegin * exp(All.Ti_begstep[P[i].TimeBinGrav] * All.Timebase_interval);
          else
            time_begstep = All.TimeBegin + All.Ti_begstep[P[i].TimeBinGrav] * All.Timebase_interval;

          dt = (P[i].TimeBinGrav ? (((integertime) 1) << P[i].TimeBinGrav) : 0) * All.Timebase_interval;

          if(All.ComovingIntegrationOn)
            dtime = All.Time * dt / All.cf_time_hubble_a;
          else
            dtime = dt;

          dtime *= All.UnitTime_in_s / All.HubbleParam;
          dtime_in_Gyr = 1e-3 * dtime / SEC_PER_MEGAYEAR;

          age_of_star_in_Gyr_endstep = get_time_difference_in_Gyr(STP(i).BirthTime, time_begstep);      /* (Note: All.Ti_begstep[] has already been advanced for the next step at this point) */
          age_of_star_in_Gyr = age_of_star_in_Gyr_endstep - dtime_in_Gyr;

          if(age_of_star_in_Gyr > LIFETIME)
            continue;       
          else
            dt = ((age_of_star_in_Gyr_endstep < LIFETIME) ? dtime_in_Gyr : (age_of_star_in_Gyr_endstep - LIFETIME));

          StarParticle[Nsource].index = i;
          StarParticle[Nsource].NumNgb = 0;
          StarParticle[Nsource].NormSph = 0;

          for(int bin = 0; bin < MRT_BINS; bin++)
            StarParticle[Nsource].TotalPhotReleased[bin] = get_photon_released_star(i, dt);

          Nsource++;
        }
    }
}

void end_stellar_sources(void)
{
  myfree(StarParticle);

  mpi_printf("RT: stellar sources done.\n");

  //TIMER_STOP(CPU_GFM_ENRICH);
}

#endif /* MRT_STARS */

#ifdef MRT_BH

struct bh_particle *BHParticle;

double get_photon_released_blackholes(int p, double dt)
{
  if(P[p].Type != 5)
    terminate("Particle should be a black hole! Type=%d, index %d", P[p].Type, p);

  double phot = pow(10.0, All.LogAGNLuminosity) * All.UVLuminosityFraction *  dt * 1e3 * SEC_PER_MEGAYEAR / (PHOT_NORMALIZATION * 2.17896e-11);

  return phot;
}

void start_blackhole_sources(void)
{
  int idx, i;
  double dt, dtime, dtime_in_Gyr;

  BHParticle = mymalloc("BHParticle", NumBHs * sizeof(struct bh_particle));

  Nsource = 0;

  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Ti_Current != All.Ti_Current)
        {
          terminate("how can this be?");
          drift_particle(i, All.Ti_Current);
        }

      if(P[i].Type == 5 && P[i].Mass > 0)
        {

          dt = (P[i].TimeBinGrav ? (((integertime) 1) << P[i].TimeBinGrav) : 0) * All.Timebase_interval;

          if(All.ComovingIntegrationOn)
            dtime = All.Time * dt / All.cf_time_hubble_a;
          else
            dtime = dt;

          dtime *= All.UnitTime_in_s / All.HubbleParam;
          dtime_in_Gyr = 1e-3 * dtime / SEC_PER_MEGAYEAR;

          BHParticle[Nsource].index    = i;
          BHParticle[Nsource].NumNgb   = 0;
          BHParticle[Nsource].NormSph  = 0;
	  BHParticle[Nsource].Dhsmlrho = BPP(i).BH_Hsml;

          for(int bin = 0; bin < MRT_BINS; bin++)
            BHParticle[Nsource].TotalPhotReleased[bin] = get_photon_released_blackholes(i, dt);

          Nsource++;
        }
    }
}

void end_blackhole_sources(void)
{
  myfree(BHParticle);

  mpi_printf("RT: BH sources done.\n");
}

#endif /* MRT_BH */

#endif /* MRT_SOURCES */
