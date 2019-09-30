/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/gradients.c
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

int N_Grad = 0;

struct grad_elements grad_elements[MAXGRADIENTS], *GDensity, *GVelx, *GVely, *GVelz, *GPressure, *GUtherm;

void init_gradients()
{
#if defined(MAXSCALARS) || defined(DVR_RENDER)
  int k;
#endif

  gradient_init(&SphP[0].Density, &PrimExch[0].Density, SphP[0].Grad.drho, GRADIENT_TYPE_DENSITY);

#if defined(SPECIAL_RELATIVITY) || defined(GENERAL_RELATIVITY) || defined(DEREFINE_GENTLY) || defined(CONDUCTION_SATURATION) || defined(NON_LINEAR_SLOPE_LIMITERS) || defined(CALCULATE_QUANTITIES_IN_POSTPROCESS)
  gradient_init(&SphP[0].Utherm, &PrimExch[0].Utherm, SphP[0].Grad.dutherm, GRADIENT_TYPE_UTHERM);
#endif

#ifdef DVR_RENDER
  for(k = 0; k < DVR_NUM_FIELDS; k++)
    gradient_init(&SphP[0].DvrFields[k], &PrimExch[0].DvrFields[k], SphP[0].Grad.dDvrFields[k], GRADIENT_TYPE_NORMAL);
#endif

  gradient_init(&P[0].Vel[0], &PrimExch[0].VelGas[0], SphP[0].Grad.dvel[0], GRADIENT_TYPE_VELX);
  gradient_init(&P[0].Vel[1], &PrimExch[0].VelGas[1], SphP[0].Grad.dvel[1], GRADIENT_TYPE_VELY);
  gradient_init(&P[0].Vel[2], &PrimExch[0].VelGas[2], SphP[0].Grad.dvel[2], GRADIENT_TYPE_VELZ);

  gradient_init(&SphP[0].Pressure, &PrimExch[0].Pressure, SphP[0].Grad.dpress, GRADIENT_TYPE_PRESSURE);

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
  gradient_init(&SphP[0].A, &PrimExch[0].A, SphP[0].Grad.dA, GRADIENT_TYPE_NORMAL);
#endif

#ifdef TGCHEM
  gradient_init(&SphP[0].Gamma, &PrimExch[0].Gamma, SphP[0].Grad.dgamma, GRADIENT_TYPE_NORMAL);
#endif

#ifdef VARIABLE_GAMMA
  gradient_init(&SphP[0].GammaE, &PrimExch[0].GammaE, SphP[0].Grad.dgammaE, GRADIENT_TYPE_NORMAL);
  gradient_init(&SphP[0].GammaC, &PrimExch[0].GammaC, SphP[0].Grad.dgammaC, GRADIENT_TYPE_NORMAL);
#endif

#ifdef MHD
  gradient_init(&SphP[0].B[0], &PrimExch[0].B[0], SphP[0].Grad.dB[0], GRADIENT_TYPE_NORMAL);
  gradient_init(&SphP[0].B[1], &PrimExch[0].B[1], SphP[0].Grad.dB[1], GRADIENT_TYPE_NORMAL);
  gradient_init(&SphP[0].B[2], &PrimExch[0].B[2], SphP[0].Grad.dB[2], GRADIENT_TYPE_NORMAL);
#endif

#ifdef MHD_CT
  gradient_init(&SphP[0].A[0], &PrimExch[0].A[0], SphP[0].Grad.dA[0], GRADIENT_TYPE_AX);
  gradient_init(&SphP[0].A[1], &PrimExch[0].A[1], SphP[0].Grad.dA[1], GRADIENT_TYPE_AY);
  gradient_init(&SphP[0].A[2], &PrimExch[0].A[2], SphP[0].Grad.dA[2], GRADIENT_TYPE_AZ);
#endif

#ifdef COSMIC_RAYS
  gradient_init(&SphP[0].CR_Pressure, &PrimExch[0].CR_Pressure, SphP[0].Grad.dcrPressure, GRADIENT_TYPE_NORMAL);
#endif

#if defined(VORONOI_PROJ_TAU)
  gradient_init(&SphP[0].Temperature, &PrimExch[0].Temperature, SphP[0].Grad.dtemp, GRADIENT_TYPE_NORMAL);
#endif

#ifdef MAXSCALARS
  MyFloat *addr;

  for(k = 0; k < N_Scalar; k++)
    {
      addr = (MyFloat *) (((char *) (&SphP[0])) + scalar_elements[k].offset);
      gradient_init(addr, &PrimExch[0].Scalars[k], SphP[0].Grad.dscalars[k], GRADIENT_TYPE_NORMAL);
    }
#endif

#ifdef TRACER_FIELD
  gradient_init(&SphP[0].Tracer, &PrimExch[0].Tracer, SphP[0].Grad.dtracer, GRADIENT_TYPE_NORMAL);
#endif

#ifdef FLD
  gradient_init(&SphP[0].n_gamma, &PrimExch[0].n_gamma, SphP[0].Grad.dngamma, GRADIENT_TYPE_FLD);
#endif

  mpi_printf("INIT: %d/%d Gradients used.\n", N_Grad, MAXGRADIENTS);
}

void gradient_init(MyFloat * addr, MyFloat * addr_exch, MySingle * addr_grad, int type)
{
  if(N_Grad == MAXGRADIENTS)
    {
      mpi_printf("Failed to register gradient, maximum of %d already reached\n", MAXGRADIENTS);
      terminate("MAXGRADIENTS reached");
    }

  grad_elements[N_Grad].type = type;

  if((type == GRADIENT_TYPE_VELX) || (type == GRADIENT_TYPE_VELY) || (type == GRADIENT_TYPE_VELZ))
    {
      /* basic structure is P */
      grad_elements[N_Grad].offset = ((char *) addr) - ((char *) &P[0]);
    }
  else
    {
      /* basic structure is SphP */
      grad_elements[N_Grad].offset = ((char *) addr) - ((char *) &SphP[0]);
    }

  grad_elements[N_Grad].offset_exch = ((char *) addr_exch) - ((char *) &PrimExch[0]);
  grad_elements[N_Grad].offset_grad = ((char *) addr_grad) - ((char *) &(SphP[0].Grad));

  switch (type)
    {
    case GRADIENT_TYPE_VELX:
      GVelx = &grad_elements[N_Grad];
      break;
    case GRADIENT_TYPE_VELY:
      GVely = &grad_elements[N_Grad];
      break;
    case GRADIENT_TYPE_VELZ:
      GVelz = &grad_elements[N_Grad];
      break;
    case GRADIENT_TYPE_DENSITY:
      GDensity = &grad_elements[N_Grad];
      break;
    case GRADIENT_TYPE_PRESSURE:
      GPressure = &grad_elements[N_Grad];
      break;
    case GRADIENT_TYPE_UTHERM:
      GUtherm = &grad_elements[N_Grad];
      break;
    default:
      break;
    }

  N_Grad++;
}
