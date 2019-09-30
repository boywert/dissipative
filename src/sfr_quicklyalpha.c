/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sfr_eEOS.c
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
#include "allvars.h"
#include "proto.h"
#include "forcetree.h"


/*

#ifdef QUICK_LYALPHA
	  temp = u_to_temp_fac * (SphP[i].Entropy) /
	    GAMMA_MINUS1 * pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1);

	  if(SphP[i].d.Density > All.OverDensThresh && temp < 1.0e5)
	    flag = 0;
	  else
	    flag = 1;
00#endif

*/



#ifdef USE_SFR
#ifdef QUICK_LYALPHA

/** \brief Main driver for star formation and gas cooling.
 *
 *  This function loops over all the active gas cells. If a given cell
 *  meets the criteria for star formation to be active the multi-phase
 *  model is activated, the properties of the cell are updated according to
 *  the latter and the star formation rate computed. In the other case, the
 *  standard isochoric cooling is applied to the gas cell by calling the function
 *  cool_cell() and the star formation rate is set to 0.
 */
void cooling_and_starformation(void)
{
  TIMER_START(CPU_COOLINGSFR);

  int idx, i, bin;
  double dt, dtime;
  double unew, du;
  double dens;
  double tsfr;

  /* note: assuming FULL ionization */
  double u_to_temp_fac = (4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC))) * PROTONMASS / BOLTZMANN * GAMMA_MINUS1 * All.UnitEnergy_in_cgs / All.UnitMass_in_g;

  /* clear the SFR stored in the active timebins */
  for(bin = 0; bin < TIMEBINS; bin++)
    if(TimeBinSynchronized[bin])
      TimeBinSfr[bin] = 0;


  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Mass == 0 && P[i].ID == 0)
        continue;               /* skip cells that have been swallowed or eliminated */

      dens = SphP[i].Density;

      dt = (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;
      dtime = All.cf_atime * dt / All.cf_time_hubble_a;


      /* apply the temperature floor (and, if relevant, BH thermal feedback) */
      unew = dmax(All.MinEgySpec, SphP[i].Utherm);
      if(unew < 0)
        terminate("Invalid Temperature: Task=%d i=%d unew=%g\n", ThisTask, i, unew);
      du = unew - SphP[i].Utherm;
      SphP[i].Utherm += du;
      SphP[i].Energy += All.cf_atime * All.cf_atime * du * P[i].Mass;

      SphP[i].Sfr = 0;


#ifdef QUICK_LYALPHA_LATETIMEONLY
      if(All.Time < All.TimeOfCoolingStart)
        continue;
#endif

      /* do cooling */
      cool_cell(i);


      /* enable star formation if gas is above SF density threshold */
      if(All.ComovingIntegrationOn)
        if(dens > All.OverDensThresh && u_to_temp_fac * SphP[i].Utherm <= All.TemperatureThresh)
          {
            /* active star formation */

            /* form the stars on the local timestep timescale */
            tsfr = dtime;

            if(dt > 0)
              {
                if(P[i].TimeBinHydro)   /* upon start-up, we need to protect against dt==0 */
                  {
                    SphP[i].Sfr = P[i].Mass / tsfr * (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);
                    TimeBinSfr[P[i].TimeBinHydro] += SphP[i].Sfr;
                  }
              }
          }
    }                           /* end of main loop over active particles */

  TIMER_STOP(CPU_COOLINGSFR);
}

/** \brief Return the star formation rate associated with the gas cell i.
 *
 *  \param i the index of the gas cell
 *  \return star formation rate in solar masses / yr
 */
double get_starformation_rate(int i)
{
  if(RestartFlag == 3)
    return SphP[i].Sfr;

  double rateOfSF = 0;
  double tsfr;
  /* note: assuming FULL ionization */
  double u_to_temp_fac = (4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC))) * PROTONMASS / BOLTZMANN * GAMMA_MINUS1 * All.UnitEnergy_in_cgs / All.UnitMass_in_g;

  double dens = SphP[i].Density;


  if(All.ComovingIntegrationOn)
    if(dens > All.OverDensThresh && u_to_temp_fac * SphP[i].Utherm <= All.TemperatureThresh)
      {
        /* active star formation */

        /* form the stars on the local timestep timescale */
        double dt = (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;
        tsfr = All.cf_atime * dt / All.cf_time_hubble_a;

        if(dt > 0)
          {
            if(P[i].TimeBinHydro)       /* upon start-up, we need to protect against dt==0 */
              {
                rateOfSF = P[i].Mass / tsfr;
                TimeBinSfr[P[i].TimeBinHydro] += SphP[i].Sfr;
              }
          }
      }

  /* convert to solar masses per yr */
  rateOfSF *= (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

  return rateOfSF;
}



/** \brief Set the appropriate units for the parameters of the multi-phase model.

void set_units_sfr(void)
{
  All.OverDensThresh = All.CritOverDensity * All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
}

*/


#endif
#endif
