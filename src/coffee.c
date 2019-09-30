/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/coffee.c
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"
#include "voronoi.h"



#ifdef COFFEE_PROBLEM

void coffee_overide_velocities(void)
{
  int idx, i;
  double vx, vy;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].ID >= 10000000)
        {
          coffee_get_velocity(P[i].Pos[0], P[i].Pos[1], &vx, &vy);

          SphP[i].VelVertex[0] = vx;
          SphP[i].VelVertex[1] = vy;
        }
    }
}

void coffee_get_velocity(double x, double y, double *vx, double *vy)
{
  double dt, r, phi, omega;

  x -= 0.5;
  y -= 0.5;

  phi = atan2(y, x);
  r = sqrt(x * x + y * y);

  omega = 2 * M_PI / 5;         /* angular velocity */

  /*  the actual time-step of the particle */
  dt = (P[0].TimeBinHydro ? (((integertime) 1) << P[0].TimeBinHydro) : 0) * All.Timebase_interval;

  if(P[0].TimeBinHydro == 0)
    {
      dt = 1.0e-5;
    }

  phi += omega * dt;

  *vx = (r * cos(phi) - x) / dt;
  *vy = (r * sin(phi) - y) / dt;
}






/*
void check_for_two_close_approaches(void)
{
  int i;
  MyIDType id1, id2;
  int q1, q2;
  double n[2], nn, vp1, vp2;

  for(i = 0; i < Nvf; i++)
    {
      id1 = VF[i].p1->ID;
      id2 = VF[i].p2->ID;

      if((id1 >= 10000000 && id1 < 20000000 && id2 < 10000000) ||
         (id2 >= 10000000 && id2 < 20000000 && id1 < 10000000)) 
        {
          if((id1 >= 10000000 && id1 < 20000000 && id2 < 10000000))
            {
              q1 = VF[i].p1->index;
              q2 = VF[i].p2->index;
            }
          else
            {
              q1 = VF[i].p2->index;
              q2 = VF[i].p1->index;
            }

          n[0] = P[q2].Pos[0] - P[q1].Pos[0];
          n[1] = P[q2].Pos[1] - P[q1].Pos[1];
          nn = sqrt(n[0]*n[0] + n[1]*n[1]);

          n[0]/=nn;
          n[1]/=nn;

          if(nn < 0.015)
            {
              vp1 = SphP[q1].VelVertex[0]*n[0] + SphP[q1].VelVertex[1]*n[1];
              vp2 = SphP[q2].VelVertex[0]*n[0] + SphP[q2].VelVertex[1]*n[1];
              
              if(vp1 > vp2) 
                {
                  SphP[q2].VelVertex[0] += -vp2*n[0] + vp1*n[0];
                  SphP[q2].VelVertex[1] += -vp2*n[1] + vp1*n[1];
                }
            }
        }
    }
}
*/



#endif
