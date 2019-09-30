/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/circumstellar/circumstellar_gravity.c
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

#include "../allvars.h"
#include "../proto.h"

#include "./circumstellar_proto.h"

#ifdef CIRCUMSTELLAR
void circumstellar_calc_gravity_from_stars_planets_only(void)
{
#if defined(GRAVITY_FROM_STARS_PLANETS_ONLY) && defined(EXTERNALGRAVITY)
  double wp, fac;
  double dx, dy, dz, r, r2;
  double h, h_inv, h3_inv, u;
  int idx, k, i;

  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      for(k = 0; k < All.TotPartSources; k++)
        {
          if(SourcePartListGlobal[k].SourceID == P[i].ID)
            continue;

          dx = P[i].Pos[0] - SourcePartListGlobal[k].Pos[0];
          dy = P[i].Pos[1] - SourcePartListGlobal[k].Pos[1];
          dz = P[i].Pos[2] - SourcePartListGlobal[k].Pos[2];

          r2 = dx * dx + dy * dy + dz * dz;

          r = sqrt(r2);
          h = All.ForceSoftening[SourcePartListGlobal[k].SofteningType];

          //using spline softening
          if(r >= h)
            {
              fac = 1 / (r2 * r);
              wp = -1 / r;
            }
          else
            {
              h_inv = 1.0 / h;
              h3_inv = h_inv * h_inv * h_inv;
              u = r * h_inv;

              if(u < 0.5)
                {
                  fac = h3_inv * (10.666666666667 + u * u * (32.0 * u - 38.4));
                  wp = h_inv * (-2.8 + u * u * (5.333333333333 + u * u * (6.4 * u - 9.6)));
                }
              else
                {
                  fac = h3_inv * (21.333333333333 - 48.0 * u + 38.4 * u * u - 10.666666666667 * u * u * u - 0.066666666667 / (u * u * u));
                  wp = h_inv * (-3.2 + 0.066666666667 / u + u * u * (10.666666666667 + u * (-16.0 + u * (9.6 - 2.133333333333 * u))));
                }
            }
          P[i].GravAccel[0] -= All.G * SourcePartListGlobal[k].Mass * fac * dx;
          P[i].GravAccel[1] -= All.G * SourcePartListGlobal[k].Mass * fac * dy;
          P[i].GravAccel[2] -= All.G * SourcePartListGlobal[k].Mass * fac * dz;

#ifdef EVALPOTENTIAL
          P[i].Potential += All.G * SourcePartListGlobal[k].Mass * wp;
#endif
          P[i].ExtPotential += All.G * SourcePartListGlobal[k].Mass * wp;


        }
    }

#endif
}
#endif
