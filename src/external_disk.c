/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/external_disk.c
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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>

#include "allvars.h"
#include "proto.h"

#ifdef EXTERNALDISKPOTENTIAL

#define DISK_SCALE_R0      1.0
#define DISK_MASS_M0       1.0


#define EXTERNALDISK_TABLE_LENGTH 10000
#define WORKSIZE 100000

static double TabPotential[EXTERNALDISK_TABLE_LENGTH];

static int first_call = 0;


double externaldisk_dphidR(double r, void *param)
{
  double Sigma0 = DISK_MASS_M0 / (2 * M_PI * DISK_SCALE_R0 * DISK_SCALE_R0);

  double y = r / (2 * DISK_SCALE_R0);

  if(y < 1.0e-4)
    return 0;

  double dphidR = 2 * M_PI * All.G * Sigma0 * y * (gsl_sf_bessel_I0(y) * gsl_sf_bessel_K0(y) - gsl_sf_bessel_I1(y) * gsl_sf_bessel_K1(y));

  return dphidR;
}




double externaldisk_potential(double r)
{
  if(first_call == 0)
    {
      first_call = 1;

      int i;
      double r, result, abserr;
      gsl_function F;
      gsl_integration_workspace *workspace;

      workspace = gsl_integration_workspace_alloc(WORKSIZE);

      for(i = 0; i < EXTERNALDISK_TABLE_LENGTH; i++)
        {
          F.function = &externaldisk_dphidR;

          r = i * All.BoxSize / EXTERNALDISK_TABLE_LENGTH;

          gsl_integration_qag(&F, 0, r, 0.0001, 1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
          TabPotential[i] = result;

          /* printf("i=%d r=%g   pot=%g\n", i, r, TabPotential[i]); */
        }

      gsl_integration_workspace_free(workspace);
    }

  int bin;
  double x, dx;

  x = r / (All.BoxSize / EXTERNALDISK_TABLE_LENGTH);
  bin = (int) x;

  if(bin >= EXTERNALDISK_TABLE_LENGTH - 1)
    terminate("bin number exceeds maximum table length");

  dx = x - bin;

  return (1 - dx) * TabPotential[bin] + dx * TabPotential[bin + 1];
}

#endif
