/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/circumstellar/circumstellar_kepler.c
 * \date        11/2013
 * \author      Diego J. Munoz
 * \brief       Routines for solving Kepler's equation
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - 13.11.2013 Standard routine to iteratively solve Kepler's equation 
 *              for general orbital eccentricities using the Danby method.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../allvars.h"
#include "../proto.h"

#include "./circumstellar_proto.h"


void circumstellar_solve_kepler(double mean_anomaly, double ecc, double *x, double *y, double *vx, double *vy)
{
  double ecc_anomaly;

  /*Danby scheme to solve Kepler's equation */
  if(sin(mean_anomaly) >= 0)
    ecc_anomaly = mean_anomaly + 0.85 * ecc;
  else
    ecc_anomaly = mean_anomaly - 0.85 * ecc;

  double delta1, delta2, delta3;
  while(1)
    {
      delta1 = -(ecc_anomaly - ecc * sin(ecc_anomaly) - mean_anomaly) / (1 - ecc * cos(ecc_anomaly));
      delta2 = -(ecc_anomaly - ecc * sin(ecc_anomaly) - mean_anomaly) / (1 - ecc * cos(ecc_anomaly) + 0.5 * delta1 * ecc * sin(ecc_anomaly));
      delta3 = -(ecc_anomaly - ecc * sin(ecc_anomaly) - mean_anomaly) / (1 - ecc * cos(ecc_anomaly) + 0.5 * delta2 * ecc * sin(ecc_anomaly) + delta2 * delta2 * ecc * cos(ecc_anomaly) / 6);

      ecc_anomaly += delta3;

      if(delta3 < 1.0e-6)
        break;

    }

  *x = (cos(ecc_anomaly) - ecc);
  *y = sqrt(1 - ecc * ecc) * sin(ecc_anomaly);

  *vx = -sin(ecc_anomaly) / (1 - ecc * cos(ecc_anomaly));
  *vy = sqrt(1 - ecc * ecc) * cos(ecc_anomaly) / (1 - ecc * cos(ecc_anomaly));

  return;
}
