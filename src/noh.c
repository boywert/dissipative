/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/noh.c
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


#ifdef NOH_PROBLEM


#ifdef TWODIMS
void set_special_noh_boundary_conditions(void)  /* 2D version */
{
  int i;
  double x, y, r, rho;

  for(i = 0; i < NumGas; i++)
    {
      x = SphP[i].Center[0] - 1.0;
      y = SphP[i].Center[1] - 1.0;

      r = sqrt(x * x + y * y);

      if(r > 0.8)
        {
          rho = 1.0 * pow(1 + All.Time / r, 1);
          P[i].Mass = SphP[i].Volume * rho;

          P[i].Vel[0] = -x / r;
          P[i].Vel[1] = -y / r;
          P[i].Vel[2] = 0;

          SphP[i].Momentum[0] = P[i].Mass * P[i].Vel[0];
          SphP[i].Momentum[1] = P[i].Mass * P[i].Vel[1];

          SphP[i].Utherm = 1.0e-6 * pow(rho, GAMMA_MINUS1) / GAMMA_MINUS1;

          SphP[i].Energy = P[i].Mass * SphP[i].Utherm + 0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);

          SphP[i].Pressure = GAMMA_MINUS1 * rho * SphP[i].Utherm;

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
          SphP[i].Entropy = P[i].Mass * log(SphP[i].Pressure / pow(SphP[i].Density, GAMMA));
#endif
        }
    }
}

#else

void set_special_noh_boundary_conditions(void)  /* 3D version */
{
  int i;
  double x, y, z, r, rho;

  for(i = 0; i < NumGas; i++)
    {
      x = SphP[i].Center[0] - 1.5;
      y = SphP[i].Center[1] - 1.5;
      z = SphP[i].Center[2] - 1.5;

      r = sqrt(x * x + y * y + z * z);

      if(r > 0.8)
        {
          rho = 1.0 * pow(1 + All.Time / r, 2);
          P[i].Mass = SphP[i].Volume * rho;

          P[i].Vel[0] = -x / r;
          P[i].Vel[1] = -y / r;
          P[i].Vel[2] = -z / r;

          SphP[i].Momentum[0] = P[i].Mass * P[i].Vel[0];
          SphP[i].Momentum[1] = P[i].Mass * P[i].Vel[1];
          SphP[i].Momentum[2] = P[i].Mass * P[i].Vel[2];

          SphP[i].Utherm = 1.0e-6 * pow(rho, GAMMA_MINUS1) / GAMMA_MINUS1;

          SphP[i].Energy = P[i].Mass * SphP[i].Utherm + 0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);

          SphP[i].Pressure = GAMMA_MINUS1 * rho * SphP[i].Utherm;
        }
    }
}
#endif

#endif
