/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/fm_star_formation/radpressthick.c
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


#if defined(RADPRESS_OPT_THICK) && defined(OTVET)


double PI = 3.141592653589;
/* Use OTVET to give number of photons per cell and compute radiation pressure */
void radpressthick(void)
{
  int idx, target, ncount, n, i, k;
  double dt, dtime, dvel, lum_per_mass;
  double local_ekin_min = MAX_REAL_NUMBER, local_ekin_max = MIN_REAL_NUMBER;
  double global_ekin_min, global_ekin_max;
  double ekin_before, ekin_after, ekin_ratio;
  double local_dvel_max = MIN_REAL_NUMBER, global_dvel_max;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      for(k = 0; k < 3; k++)
        SphP[i].RadPress[k] = 0.;
    }

  /* [c_light] = UnitVelocity */
  double c_light = CLIGHT / All.UnitVelocity_in_cm_per_s;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      dt = (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;
/*
          dt = (All.otvet_Radiation_Ti_endstep - All.otvet_Radiation_Ti_begstep) * All.Timebase_interval;
*/
      if(All.ComovingIntegrationOn)
        dt /= All.cf_hubble_a;
      dtime = dt;

      /* subtract kin. energy */
      SphP[i].Energy -= 0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);
      ekin_before = 0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);

      /* do the kick for gas cells */
      for(k = 0; k < 3; k++)
        {
#ifndef OTVET_MULTI_FREQUENCY

          SphP[i].RadPress[k] = (SphP[i].n_gamma_abs[0] / SphP[i].Density) * SphP[i].vector_n[k];

          SphP[i].RadPress[k] *= 13.6 * ELECTRONVOLT_IN_ERGS * (All.HubbleParam / All.UnitEnergy_in_cgs) / c_light;
#endif
          SphP[i].RadPress[k] /= P[i].Mass;

          dvel = SphP[i].RadPress[k];

          /* add velocity */
          P[i].Vel[k] += dvel;

          /* add momentum */
          SphP[i].Momentum[k] += P[i].Mass * dvel;


          if(!gsl_finite(dvel))
            terminate("RADPRESS_OPT_THICK: bad radiation pressure dvel=%g P[i].Mass=%g SphP[i].RadPress[k]=%g dtime=%g ID=%d", dvel, P[i].Mass, SphP[i].RadPress[k], dtime, P[i].ID);
        }

      /* add kin. energy */
      SphP[i].Energy += 0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);
      ekin_after = 0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);

      ekin_ratio = (ekin_before > 0) ? (ekin_after / ekin_before) : 1.0;

      if(ekin_ratio < local_ekin_min)
        local_ekin_min = ekin_ratio;

      if(ekin_ratio > local_ekin_max)
        local_ekin_max = ekin_ratio;

    }

  MPI_Reduce(&local_ekin_min, &global_ekin_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_ekin_max, &global_ekin_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  mpi_printf("RADPRESS_OPT_THICK: kin. energy ratio = (%g - %g) \n", global_ekin_min, global_ekin_max);
}

#endif
