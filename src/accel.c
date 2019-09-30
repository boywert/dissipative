/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/accel.c
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


/*! \file accel.c
 *  \brief Driver routines to carry out gravity force computation
 */


/*! \brief This routine computes the gravity accelerations for all active particles.
 *
 * If the particle mesh is used and the current time step
 * requires a PM force computation, new long range forces are
 * computed by long_range_force(). Then the shortrange tree forces
 * are computed by gravity(). The force tree is rebuild every time step.
 */
void compute_grav_accelerations(int timebin, int fullflag)
{
  if(TimeBinsGravity.GlobalNActiveParticles > 0)
    {
#ifdef FMM
      fmm_gravity(timebin);
#else
      if(All.TypeOfOpeningCriterion == 1 && All.Ti_Current == 0 && All.ErrTolTheta > 0)
        {
          /* For the first timestep, we do one gravity calculation up front
           * with the Barnes & Hut Criterion to allow usage of relative opening
           * criterion with consistent accuracy.
           */
#ifdef PMGRID
          long_range_force();
#endif
          gravity(timebin, fullflag);
        }

      gravity(timebin, fullflag);       /* computes (short-range) gravity accel. */
#endif

#ifdef FORCETEST
      gravity_forcetest();
#endif
    }

#ifdef SECOND_ORDER_ICS
  if(All.Ti_Current == 0 && RestartFlag == 0)
    second_order_ics();         /* produces the actual ICs from the special second order IC file */
#endif
}


/*! \brief main driver routine of tree force calculation
 *
 *  This routine handles the whole tree force calculation. First it
 *  build a new force tree force_treebuild() every timestep. This tree is then
 *  used to calculate a new tree force for every active particle ( gravity_tree() ).
 */
void gravity(int timebin, int fullflag)
{
  double tstart = second();

#if defined(SELFGRAVITY) || defined(TREE_RAD) || defined(RADCOOL)
  /* set new softening lengths on global steps to take into account possible cosmological time variation */
  if(timebin == All.HighestOccupiedGravTimeBin)
    set_softenings();


#ifdef ALLOW_DIRECT_SUMMATION
  if(TimeBinsGravity.GlobalNActiveParticles < DIRECT_SUMMATION_THRESHOLD)
    {
      gravity_direct(timebin);

#ifndef ONEDIMS_SPHERICAL
      gravity_force_finalize(timebin);
#endif

#ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE
  calc_exact_gravity_for_particle_type();
#endif


#ifdef EXTERNALGRAVITY
  gravity_external();
#endif

    }
  else
#endif
    {
#ifdef ONEDIMS_SPHERICAL
      gravity_monopole_1d_spherical();
#else
      force_treeallocate(NumPart, All.MaxPart);
      if(TimeBinsGravity.GlobalNActiveParticles >= 10 * NTask)
        force_treebuild(NumPart, 1, 0, timebin);
      else
        force_treebuild(NumPart, 0, 0, timebin);

      gravity_tree(timebin);

#ifndef ONEDIMS_SPHERICAL
      gravity_force_finalize(timebin);
#endif

#ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE
      calc_exact_gravity_for_particle_type();
#endif

#ifdef EXTERNALGRAVITY
      gravity_external();
#endif

      /* note: we here moved 'gravity_force_finalize' in front of the non-standard physics call so that the BH positioning
       * search sees the total gravitational potential
       */
      if(fullflag == FLAG_FULL_TREE && RestartFlag != 18)
        calculate_non_standard_physics_with_valid_gravity_tree();

      /* this is for runs which have the full tree at each time step; no HIERARCHICAL_GRAVITY */
      calculate_non_standard_physics_with_valid_gravity_tree_always();

      myfree(Father);
      myfree(Nextnode);
#ifdef BLACK_HOLES
      myfree(Tree_AuxBH_Points);
#endif
      myfree(Tree_Points);
      force_treefree();
#endif
    }

#else

  /* self-gravity is switched off */
  int idx, i, j;
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;
#ifdef EVALPOTENTIAL
      P[i].Potential = 0;
#endif
      for(j = 0; j < 3; j++)
        P[i].GravAccel[j] = 0;
    }

#ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE
    calc_exact_gravity_for_particle_type();
#endif


#ifdef EXTERNALGRAVITY
    gravity_external();
#endif

#endif

  double tend = second();
  mpi_printf("GRAVITY: done for timebin %d,  %lld particles  (took %g sec)\n", timebin, TimeBinsGravity.GlobalNActiveParticles, timediff(tstart, tend));
}



void gravity_force_finalize(int timebin)
{
  int i, j, idx;
  double ax, ay, az;

  TIMER_START(CPU_TREE);

  /* now add things for comoving integration */

#ifdef GRAVITY_NOT_PERIODIC
#ifndef PMGRID
  if(All.ComovingIntegrationOn)
    {
      double fac = 0.5 * All.Hubble * All.Hubble * All.Omega0 / All.G;

      for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
        {
          i = TimeBinsGravity.ActiveParticleList[idx];
          if(i < 0)
            continue;

          for(j = 0; j < 3; j++)
            P[i].GravAccel[j] += fac * P[i].Pos[j];
        }
    }
#endif
#endif

#ifdef HIERARCHICAL_GRAVITY
  if(timebin == All.HighestOccupiedGravTimeBin)
#endif
    {
      mpi_printf("GRAVTREE: Setting OldAcc!\n");

      for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
        {
          i = TimeBinsGravity.ActiveParticleList[idx];
          if(i < 0)
            continue;

#ifdef PMGRID
          ax = P[i].GravAccel[0] + P[i].GravPM[0] / All.G;
          ay = P[i].GravAccel[1] + P[i].GravPM[1] / All.G;
          az = P[i].GravAccel[2] + P[i].GravPM[2] / All.G;
#else
          ax = P[i].GravAccel[0];
          ay = P[i].GravAccel[1];
          az = P[i].GravAccel[2];
#endif

#ifdef SECOND_ORDER_ICS
          if(All.Ti_Current == 0 && RestartFlag == 0)
            continue;           /* to prevent that we overwrite OldAcc in the first evaluation */
#endif
          P[i].OldAcc = sqrt(ax * ax + ay * ay + az * az);
        }
    }


  /*  muliply by G */
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      for(j = 0; j < 3; j++)
        P[i].GravAccel[j] *= All.G;

#ifdef EVALPOTENTIAL

#if defined(PMGRID) && defined(PERIODIC) && !defined(GRAVITY_NOT_PERIODIC)
      P[i].Potential += All.MassPMregions[0] * M_PI / (All.Asmth[0] * All.Asmth[0] * boxSize_X * boxSize_Y * boxSize_Z);  
#ifdef PLACEHIGHRESREGION
      P[i].Potential += All.MassPMregions[1] * M_PI / (All.Asmth[1] * All.Asmth[1] * boxSize_X * boxSize_Y * boxSize_Z);
#endif
#endif

      /* It's better to not remove the self-potential here to get a smooth potential field for co-spatial particles with varying mass or softening.
       * For calculating the binding energy of a particle, the self-energy should then be removed as
       *
       *  P[i].Potential += P[i].Mass / (All.ForceSoftening[P[i].SofteningType] / 2.8);
       */

      P[i].Potential *= All.G;

#ifdef PMGRID
#ifndef FORCETEST_TESTFORCELAW
      P[i].Potential += P[i].PM_Potential;      /* add in long-range potential */
#endif
#endif
#endif
      if(All.ComovingIntegrationOn)
        {
#ifdef GRAVITY_NOT_PERIODIC
          double fac, r2;
          int k;

          fac = -0.5 * All.Omega0 * All.Hubble * All.Hubble;

          for(k = 0, r2 = 0; k < 3; k++)
            r2 += P[i].Pos[k] * P[i].Pos[k];
#ifdef EVALPOTENTIAL
          P[i].Potential += fac * r2;
#endif
#endif
        }
      else
        {
          double fac, r2;
          int k;

          fac = -0.5 * All.OmegaLambda * All.Hubble * All.Hubble;

          if(fac != 0)
            {
              for(k = 0, r2 = 0; k < 3; k++)
                r2 += P[i].Pos[k] * P[i].Pos[k];
#ifdef EVALPOTENTIAL
              P[i].Potential += fac * r2;
#endif
            }
        }
    }

  /* Finally, the following factor allows a computation of a cosmological simulation
     with vacuum energy in physical coordinates */
#ifdef GRAVITY_NOT_PERIODIC
#ifndef PMGRID
  if(All.ComovingIntegrationOn == 0)
    {
      double fac = All.OmegaLambda * All.Hubble * All.Hubble;

      for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
        {
          i = TimeBinsGravity.ActiveParticleList[idx];
          if(i < 0)
            continue;

          for(j = 0; j < 3; j++)
            P[i].GravAccel[j] += fac * P[i].Pos[j];
        }
    }
#endif
#endif

  TIMER_STOP(CPU_TREE);
}
