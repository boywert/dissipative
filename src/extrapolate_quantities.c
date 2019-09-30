/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/extrapolate_quantities.c
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

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_sf_gamma.h>
#include "allvars.h"
#include "proto.h"
#include "voronoi.h"


void face_time_advance_gradients(struct grad_data *delta_grad, struct state *st)
{
#if defined(SECOND_DERIVATIVES) && defined(RECONSTRUCT_GRADIENTS)

#if defined (MESHRELAX) || defined (DISABLE_TIME_EXTRAPOLATION)
  /* do not do time extrapolation */
  return;
#endif

  int l;
  for(l = 0; l < 3; l++)
    {
#ifdef VISCOSITY
      delta_grad->dvel[0][l] =
        -st->dt_half * (st->grad->dvel[0][0] * st->grad->dvel[0][l] +
                        st->grad->dvel[1][l] * st->grad->dvel[0][1] +
                        st->grad->dvel[2][l] * st->grad->dvel[0][2] -
                        st->grad->drho[l] * st->grad->dpress[0] * 1.0 / st->rho * 1.0 / st->rho +
                        st->velx * st->hessian->ddvelx[0][l] + st->vely * st->hessian->ddvelx[1][l] + st->velz * st->hessian->ddvelx[2][l] + st->hessian->ddpress[0][l] * 1.0 / st->rho);

      delta_grad->dvel[1][l] =
        -st->dt_half * (st->grad->dvel[1][0] * st->grad->dvel[0][l] +
                        st->grad->dvel[1][l] * st->grad->dvel[1][1] +
                        st->grad->dvel[2][l] * st->grad->dvel[1][2] -
                        st->grad->drho[l] * st->grad->dpress[1] * 1.0 / st->rho * 1.0 / st->rho +
                        st->velx * st->hessian->ddvely[0][l] + st->vely * st->hessian->ddvely[1][l] + st->velz * st->hessian->ddvely[2][l] + st->hessian->ddpress[1][l] * 1.0 / st->rho);

      delta_grad->dvel[2][l] =
        -st->dt_half * (st->grad->dvel[2][0] * st->grad->dvel[0][l] +
                        st->grad->dvel[1][l] * st->grad->dvel[2][1] +
                        st->grad->dvel[2][l] * st->grad->dvel[2][2] -
                        st->grad->drho[l] * st->grad->dpress[2] * 1.0 / st->rho * 1.0 / st->rho +
                        st->velx * st->hessian->ddvelz[0][l] + st->vely * st->hessian->ddvelz[1][l] + st->velz * st->hessian->ddvelz[2][l] + st->hessian->ddpress[2][l] * 1.0 / st->rho);
#endif

#ifdef THERMAL_CONDUCTION
      double gamma = GAMMA;
      delta_grad->drho[l] =
        -st->dt_half * (st->grad->dvel[0][l] * st->grad->drho[0] + st->grad->drho[l] * st->grad->dvel[0][0] +
                        st->grad->dvel[1][l] * st->grad->drho[1] + st->grad->drho[l] * st->grad->dvel[1][1] +
                        st->grad->dvel[2][l] * st->grad->drho[2] + st->grad->drho[l] * st->grad->dvel[2][2] +
                        st->velx * st->hessian->ddrho[0][l] + st->rho * st->hessian->ddvelx[0][l] +
                        st->vely * st->hessian->ddrho[1][l] + st->rho * st->hessian->ddvely[1][l] + st->velz * st->hessian->ddrho[2][l] + st->rho * st->hessian->ddvelz[2][l]);

      delta_grad->dpress[l] =
        -st->dt_half * (gamma * st->grad->dpress[l] * st->grad->dvel[0][0] +
                        st->grad->dvel[0][l] * st->grad->dpress[0] +
                        gamma * st->grad->dpress[l] * st->grad->dvel[1][1] +
                        st->grad->dvel[1][l] * st->grad->dpress[1] +
                        gamma * st->grad->dpress[l] * st->grad->dvel[2][2] +
                        st->grad->dvel[2][l] * st->grad->dpress[2] +
                        gamma * st->press * st->hessian->ddvelx[0][l] +
                        st->velx * st->hessian->ddpress[0][l] +
                        gamma * st->press * st->hessian->ddvely[1][l] + st->vely * st->hessian->ddpress[1][l] + gamma * st->press * st->hessian->ddvelz[2][l] + st->velz * st->hessian->ddpress[2][l]);
#endif

#if defined(TRACER_FIELD) && defined(TRACER_DIFFUSION)
      if(l == 0)
        delta_grad->dtracer[l] =
          -st->dt_half * ((st->grad->dvel[0][0] + st->grad->dvel[1][1] + st->grad->dvel[2][2]) *
                          st->grad->dtracer[l] + st->velx * (st->hessian->ddtracer[0][0] + st->hessian->ddtracer[1][1] + st->hessian->ddtracer[2][2]));
      if(l == 1)
        delta_grad->dtracer[l] =
          -st->dt_half * ((st->grad->dvel[0][0] + st->grad->dvel[1][1] + st->grad->dvel[2][2]) *
                          st->grad->dtracer[l] + st->vely * (st->hessian->ddtracer[0][0] + st->hessian->ddtracer[1][1] + st->hessian->ddtracer[2][2]));
      if(l == 2)
        delta_grad->dtracer[l] =
          -st->dt_half * ((st->grad->dvel[0][0] + st->grad->dvel[1][1] + st->grad->dvel[2][2]) *
                          st->grad->dtracer[l] + st->velz * (st->hessian->ddtracer[0][0] + st->hessian->ddtracer[1][1] + st->hessian->ddtracer[2][2]));

#endif

    }
#endif
}


void face_space_extrapolate_gradients(struct grad_data *delta_grad, struct state *st)
{
#if defined(SECOND_DERIVATIVES) && defined(RECONSTRUCT_GRADIENTS)

#ifdef DISABLE_SPATIAL_RECONSTRUCTION
  return;
#endif


#ifdef VISCOSITY
  face_extrapolate_gradient(delta_grad->dvel[0], st->hessian->ddvelx[0], st->dx, st->dy, st->dz);
  face_extrapolate_gradient(delta_grad->dvel[1], st->hessian->ddvely[0], st->dx, st->dy, st->dz);
  face_extrapolate_gradient(delta_grad->dvel[2], st->hessian->ddvelz[0], st->dx, st->dy, st->dz);
#endif

#ifdef THERMAL_CONDUCTION
  face_extrapolate_gradient(delta_grad->drho, st->hessian->ddrho[0], st->dx, st->dy, st->dz);
  face_extrapolate_gradient(delta_grad->dpress, st->hessian->ddpress[0], st->dx, st->dy, st->dz);
#endif

#if defined(TRACER_FIELD) && defined(TRACER_DIFFUSION)
  face_extrapolate_gradient(delta_grad->dtracer, st->hessian->ddtracer[0], st->dx, st->dy, st->dz);
#endif

#endif
}


void face_extrapolate_gradient(double *delta, double *hessian, double dx, double dy, double dz)
{
#if defined(SECOND_DERIVATIVES) && defined(RECONSTRUCT_GRADIENTS)
  int l;
  for(l = 0; l < 3; l++)
    delta[l] = hessian[l * 3] * dx + hessian[l * 3 + 1] * dy + hessian[l * 3 + 2] * dz;
#endif
}


void face_add_gradient_extrapolations(struct state *st, struct grad_data *delta_grad_time, struct grad_data *delta_grad_space)
{
#if defined(SECOND_DERIVATIVES) && defined(RECONSTRUCT_GRADIENTS)
  int l;
  for(l = 0; l < 3; l++)
    {

#if !defined(MESHRELAX) && !defined(DISABLE_TIME_EXTRAPOLATION)
#ifdef VISCOSITY
      st->grad->dvel[0][l] += delta_grad_time->dvel[0][l];
      st->grad->dvel[1][l] += delta_grad_time->dvel[1][l];
      st->grad->dvel[2][l] += delta_grad_time->dvel[2][l];
#endif

#ifdef THERMAL_CONDUCTION
      st->grad->drho[l] += delta_grad_time->drho[l];
      st->grad->dpress[l] += delta_grad_time->dpress[l];
#endif

#if defined(TRACER_FIELD) && defined(TRACER_DIFFUSION)
      st->grad->dtracer[l] += delta_grad_time->dtracer[l];
#endif
#endif

#if !defined(DISABLE_SPATIAL_EXTRAPOLATION)
#ifdef VISCOSITY
      st->grad->dvel[0][l] += delta_grad_space->dvel[0][l];
      st->grad->dvel[1][l] += delta_grad_space->dvel[1][l];
      st->grad->dvel[2][l] += delta_grad_space->dvel[2][l];
#endif

#ifdef THERMAL_CONDUCTION

      st->grad->drho[l] += delta_grad_space->drho[l];
      st->grad->dpress[l] += delta_grad_space->dpress[l];
#endif

#if defined(TRACER_FIELD) && defined(TRACER_DIFFUSION)
      st->grad->dtracer[l] += delta_grad_space->dtracer[l];
#endif
#endif


    }

#endif
}
