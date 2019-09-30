/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/subfind/subfind_density.c
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
#include <sys/stat.h>
#include <sys/types.h>

#include "../allvars.h"
#include "../proto.h"

#if defined(BLACK_HOLES) && defined(BH_FRICTION)

#ifdef BH_DRAG
#error "BH_DRAG and BH_FRICTION should not be used together"
#endif

static void blackhole_get_friction_acceleration(int i, MyFloat * Acc, double *dtmax);

void blackhole_friction_store_previous_minimum(void)
{
  int idx, i, k;

  for(idx = 0; idx < TimeBinsBHAccretion.NActiveParticles; idx++)
    {
      i = TimeBinsBHAccretion.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(BPP(i).BH_MinPot_ActiveM >= BH_REPOSITION_POTMIN_TRUST_THRESHOLD * BPP(i).BH_MinPot_TotalM)
        {
          BPP(i).BH_MinPotTime_Previous = BPP(i).BH_MinPotTime;
          BPP(i).BH_MinPotTime = All.Ti_Current;

          for(k = 0; k < 3; k++)
            BPP(i).BH_MinPotPos_Previous[k] = BPP(i).BH_MinPotPos[k];
        }
    }
}

void blackhole_friction_update_vel_pot_minimum(void)
{
#if defined(PERIODIC) && !defined(GRAVITY_NOT_PERIODIC)
  double xtmp, ytmp, ztmp;
#endif
  int idx, i;
  double dtavg;

  for(idx = 0; idx < TimeBinsBHAccretion.NActiveParticles; idx++)
    {
      i = TimeBinsBHAccretion.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(BPP(i).BH_MinPot_ActiveM >= BH_REPOSITION_POTMIN_TRUST_THRESHOLD * BPP(i).BH_MinPot_TotalM)
        {
          if(BPP(i).BH_MinPotTime - BPP(i).BH_MinPotTime_Previous > 0 && BPP(i).BH_MinPotTime >= 0 && BPP(i).BH_MinPotTime_Previous >= 0)
            {
              double dt, dloga = (BPP(i).BH_MinPotTime - BPP(i).BH_MinPotTime_Previous) * All.Timebase_interval;

              if(All.ComovingIntegrationOn)
                {
                  dt = dloga / All.cf_hubble_a;
                  dtavg = All.BHFrictionAvgTime / All.cf_hubble_a;
                }
              else
                {
                  dt = dloga;
                  dtavg = All.BHFrictionAvgTime / All.Hubble;
                }

              if(dtavg > BPP(i).BH_MinPotCumAvgTime)
                dtavg = BPP(i).BH_MinPotCumAvgTime;

              BPP(i).BH_MinPotCumAvgTime += dt;

              double dx = GRAVITY_NEAREST_X(BPP(i).BH_MinPotPos[0] - BPP(i).BH_MinPotPos_Previous[0]);
              double dy = GRAVITY_NEAREST_Y(BPP(i).BH_MinPotPos[1] - BPP(i).BH_MinPotPos_Previous[1]);
              double dz = GRAVITY_NEAREST_Z(BPP(i).BH_MinPotPos[2] - BPP(i).BH_MinPotPos_Previous[2]);

              if((dtavg + dt) > 0)
                {
                  BPP(i).BH_MinPotVel[0] = (BPP(i).BH_MinPotVel[0] * dtavg + All.cf_atime * All.cf_atime * dx) / (dtavg + dt);
                  BPP(i).BH_MinPotVel[1] = (BPP(i).BH_MinPotVel[1] * dtavg + All.cf_atime * All.cf_atime * dy) / (dtavg + dt);
                  BPP(i).BH_MinPotVel[2] = (BPP(i).BH_MinPotVel[2] * dtavg + All.cf_atime * All.cf_atime * dz) / (dtavg + dt);
                }
            }
        }
    }
}

void blackhole_friction_apply(void)
{
  int idx, i, j, k;
  double dt_gravkick, dtmax;

  TIMER_START(CPU_BH);

  for(idx = 0; idx < TimeBinsBHAccretion.NActiveParticles; idx++)
    {
      i = TimeBinsBHAccretion.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(BPP(i).BH_MinPot_ActiveM >= BH_REPOSITION_POTMIN_TRUST_THRESHOLD * BPP(i).BH_MinPot_TotalM)
        {
#ifdef BH_FRICTION_AGGRESSIVE
          if((BPP(i).BH_MinPot < 0.5 * BHPOTVALUEINIT) && (BPP(i).BH_MinPot != 0) && (BPP(i).BH_MinPot_Extended < 0.5 * BHPOTVALUEINIT) && (BPP(i).BH_MinPot_Extended != 0))
            {
              if(BPP(i).BH_MinPot_Extended < BPP(i).BH_MinPot)
                {
                  for(k = 0; k < 3; k++)
                    P[i].Pos[k] = BPP(i).BH_MinPotPos_Extended[k];
                }
              else
                {
                  for(k = 0; k < 3; k++)
                    P[i].Pos[k] = BPP(i).BH_MinPotPos[k];
                }
            }
#else
          if((BPP(i).BH_MinPot < 0.5 * BHPOTVALUEINIT) && (BPP(i).BH_MinPot != 0) &&
             (BPP(i).BH_MinPot_Extended < 0.5 * BHPOTVALUEINIT) && (BPP(i).BH_MinPot_Extended != 0) && (BPP(i).BH_MinPot_Extended < BPP(i).BH_MinPot))
            {
              for(k = 0; k < 3; k++)
                P[i].Pos[k] = BPP(i).BH_MinPotPos_Extended[k];
            }
#endif


          if(BPP(i).BH_MinPotTime >= 0 && BPP(i).BH_MinPotTime_Previous >= 0)
            {
              if(All.ComovingIntegrationOn)
                dt_gravkick = get_gravkick_factor(BPP(i).BH_MinPotTime_Previous, BPP(i).BH_MinPotTime);
              else
                dt_gravkick = (BPP(i).BH_MinPotTime - BPP(i).BH_MinPotTime_Previous) * All.Timebase_interval;

              MyFloat Accel[3];
              blackhole_get_friction_acceleration(i, Accel, &dtmax);

              if(dt_gravkick > dtmax)
                dt_gravkick = dtmax;

              /* apply the kick */
              for(j = 0; j < 3; j++)
                P[i].Vel[j] += Accel[j] * dt_gravkick;
            }
        }
    }

  TIMER_STOP(CPU_BH);
}

/* calculate friction acceleration */
static void blackhole_get_friction_acceleration(int i, MyFloat * Acc, double *dtmax)
{
  int k;

  if(BPP(i).BH_RhoTot > 0)
    {
      double tdyn = 1.0 / sqrt(4 * M_PI * All.G * BPP(i).BH_RhoTot);
      double vrel[3];

      double cosmofac = 1.0 / sqrt(All.cf_atime);

      for(k = 0; k < 3; k++)
        vrel[k] = P[i].Vel[k] - BPP(i).BH_MinPotVel[k];

      for(k = 0; k < 3; k++)
        Acc[k] = -cosmofac * All.BHFrictionCoefficient * vrel[k] / tdyn;

      *dtmax = tdyn / (cosmofac * All.BHFrictionCoefficient);
    }
  else
    {
      for(k = 0; k < 3; k++)
        Acc[k] = 0;

      *dtmax = 0;
    }
}




#ifdef BH_HARMONIC_OSCILLATOR_FORCE
void blackhole_harmonic_force(void)
{
  for(int idx = 0; idx < TimeBinsBHAccretion.NActiveParticles; idx++)
     {
       int i = TimeBinsBHAccretion.ActiveParticleList[idx];
       if(i < 0)
         continue;

       if((BPP(i).BH_MinPot < 0.5 * BHPOTVALUEINIT) && (BPP(i).BH_MinPot != 0) &&
                    (BPP(i).BH_MinPot_Extended < 0.5 * BHPOTVALUEINIT) && (BPP(i).BH_MinPot_Extended != 0) && (BPP(i).BH_MinPot_Extended > BPP(i).BH_MinPot))
	 {

	   for(int k = 0; k < 3; k++)
	     P[i].GravAccel[k] = - 4.0 * M_PI / 3.0 * All.G * BPP(i).BH_RhoTot * ( P[i].Pos[k] - BPP(i).BH_MinPotPos[k]);

	 }
     }
}
#endif







#endif
