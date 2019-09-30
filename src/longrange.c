/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/longrange.c
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


/*! \file longrange.c
 *  \brief driver routines for computation of long-range gravitational PM force
 */

#ifdef PMGRID

/*! \brief Driver routine to call initialization of periodic or/and non-periodic FFT
 *  routines.
 *
 *  Either pm_init_periodic() or pm_init_nonperiodic() is called.
 */
void long_range_init(void)
{
#ifndef GRAVITY_NOT_PERIODIC
  pm_init_periodic();
#ifdef TWODIMS
  pm2d_init_periodic();
#endif
#ifdef PLACEHIGHRESREGION
  pm_init_nonperiodic();
#endif
#else
  pm_init_nonperiodic();
#endif
}

/*! \brief Driver routine to determine the extend of the non-
 * periodic or high resolution region
 *
 * The initialization is done by pm_init_regionsize(). Afterwards
 * the convolution kernels are computed by pm_setup_nonperiodic_kernel()
 */
void long_range_init_regionsize(void)
{
#ifndef GRAVITY_NOT_PERIODIC
#ifdef PLACEHIGHRESREGION
  if(RestartFlag != 1)
    pm_init_regionsize();
  pm_setup_nonperiodic_kernel();
#endif
#else
  if(RestartFlag != 1)
    pm_init_regionsize();
  pm_setup_nonperiodic_kernel();
#endif
}


/*! \brief This function computes the long-range PM force for all particles.
 *
 * In case of a periodic grid the force is calculated by pmforce_periodic()
 * otherwise by pmforce_nonperiodic(). If a high resolution region is specified
 * for the PM force, pmforce_nonperiodic() calculates that force in both cases.
 */
void long_range_force(void)
{
  int i;

  TIMER_START(CPU_PM_GRAVITY);


#ifdef GRAVITY_NOT_PERIODIC
  int j;
  double fac;
#endif


  for(i = 0; i < NumPart; i++)
    {
      P[i].GravPM[0] = P[i].GravPM[1] = P[i].GravPM[2] = 0;
#ifdef EVALPOTENTIAL
      P[i].PM_Potential = 0;
#endif
    }

#ifndef SELFGRAVITY
  return;
#endif


#ifndef GRAVITY_NOT_PERIODIC

#ifndef SIDM_MAXWELLIAN
#ifdef TWODIMS
  pm2d_force_periodic(0);
#else
  pmforce_periodic(0, NULL);
#endif
#endif

#ifdef PLACEHIGHRESREGION
  i = pmforce_nonperiodic(1);

  if(i == 1)                    /* this is returned if a particle lied outside allowed range */
    {
      pm_init_regionsize();
      pm_setup_nonperiodic_kernel();
      i = pmforce_nonperiodic(1);       /* try again */
    }
  if(i == 1)
    terminate("despite we tried to increase the region, we still don't fit all particles in it");
#endif


#else /* non periodic PM mesh */
  i = pmforce_nonperiodic(0);

  if(i == 1)                    /* this is returned if a particle lied outside allowed range */
    {
      pm_init_regionsize();
      pm_setup_nonperiodic_kernel();
      i = pmforce_nonperiodic(0);       /* try again */
    }
  if(i == 1)
    terminate("despite we tried to increase the region, somehow we still don't fit all particles in it");
#ifdef PLACEHIGHRESREGION
  i = pmforce_nonperiodic(1);

  if(i == 1)                    /* this is returned if a particle lied outside allowed range */
    {
      pm_init_regionsize();
      pm_setup_nonperiodic_kernel();

      /* try again */

      for(i = 0; i < NumPart; i++)
        P[i].GravPM[0] = P[i].GravPM[1] = P[i].GravPM[2] = 0;

      i = pmforce_nonperiodic(0) + pmforce_nonperiodic(1);

    }
  if(i != 0)
    terminate("despite we tried to increase the region, somehow we still don't fit all particles in it");
#endif
#endif


#ifdef GRAVITY_NOT_PERIODIC
  if(All.ComovingIntegrationOn)
    {
      fac = 0.5 * All.Hubble * All.Hubble * All.Omega0;

      for(i = 0; i < NumPart; i++)
        for(j = 0; j < 3; j++)
          P[i].GravPM[j] += fac * P[i].Pos[j];
    }


  /* Finally, the following factor allows a computation of cosmological simulation
     with vacuum energy in physical coordinates */

  if(All.ComovingIntegrationOn == 0)
    {
      fac = All.OmegaLambda * All.Hubble * All.Hubble;

      for(i = 0; i < NumPart; i++)
        for(j = 0; j < 3; j++)
          P[i].GravPM[j] += fac * P[i].Pos[j];
    }
#endif


  TIMER_STOP(CPU_PM_GRAVITY);

#ifndef LEGACY_DISPLACEMENT_CONSTRAINT
  find_long_range_step_constraint();
#endif
}


#endif
