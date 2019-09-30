/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/diffusion_fluxes.c
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

#define TCOND_MEAN_MOL_WT 18.01528      /* water */

void face_get_gradients(struct state *st_L, struct state *st_R, struct state_face *st_face, struct fluxes *flux)
{
#if defined(VISCOSITY) || defined(THERMAL_CONDUCTION) || defined(TRACER_DIFFUSION)

  int l;
#if defined (SECOND_DERIVATIVES) && defined(RECONSTRUCT_GRADIENTS)

#ifdef VISCOSITY
  int k;
  for(l = 0; l < 3; l++)
    for(k = 0; k < 3; k++)
      st_face->vel_grad[k][l] = 0.5 * (st_R->grad->dvel[k][l] + st_L->grad->dvel[k][l]);
#endif

#ifdef THERMAL_CONDUCTION
  double prefactor = TCOND_MEAN_MOL_WT * PROTONMASS / BOLTZMANN;
  for(l = 0; l < 3; l++)
    st_face->dTemp[l] =
      0.5 * prefactor * (st_R->grad->dpress[l] / st_R->rho -
                         st_R->grad->drho[l] * st_R->press / (st_R->rho * st_R->rho) + st_L->grad->dpress[l] / st_L->rho - st_L->grad->drho[l] * st_L->press / (st_L->rho * st_L->rho));
#endif

#ifdef TRACER_DIFFUSION
  for(l = 0; l < 3; l++)
    st_face->dConservedTracer[l] = 0.5 * (st_R->grad->dtracer[l] * st_R->rho + st_R->grad->drho[l] * st_R->tracer + st_L->grad->dtracer[l] * st_L->rho + st_L->grad->drho[l] * st_L->tracer);
#endif

#else
  /*If gradients are not being reconstructed and extrapolated, choose the upwind gradients. */

  if(flux->mass > 0)
    {

#ifdef VISCOSITY
      int k;
      for(l = 0; l < 3; l++)
        for(k = 0; k < 3; k++)
          st_face->vel_grad[k][l] = st_L->grad->dvel[k][l];
#endif
#ifdef THERMAL_CONDUCTION
      double prefactor = TCOND_MEAN_MOL_WT * PROTONMASS / BOLTZMANN;
      for(l = 0; l < 3; l++)
        st_face->dTemp[l] = prefactor * (st_L->grad->dpress[l] / st_L->rho - st_L->grad->drho[l] * st_L->press / (st_L->rho * st_L->rho));
#endif
#ifdef TRACER_DIFFUSION
      for(l = 0; l < 3; l++)
        st_face->dConservedTracer[l] = st_L->grad->dtracer[l] * st_L->rho + st_L->grad->drho[l] * st_L->tracer;
#endif

    }
  else
    {

#ifdef VISCOSITY
      int k;
      for(l = 0; l < 3; l++)
        for(k = 0; k < 3; k++)
          st_face->vel_grad[k][l] = st_R->grad->dvel[k][l];
#endif
#ifdef THERMAL_CONDUCTION
      double prefactor = TCOND_MEAN_MOL_WT * PROTONMASS / BOLTZMANN;
      for(l = 0; l < 3; l++)
        st_face->dTemp[l] = prefactor * (st_R->grad->dpress[l] / st_R->rho - st_R->grad->drho[l] * st_R->press / (st_R->rho * st_R->rho));
#endif
#ifdef TRACER_DIFFUSION
      for(l = 0; l < 3; l++)
        st_face->dConservedTracer[l] = st_R->grad->dtracer[l] * st_R->rho + st_R->grad->drho[l] * st_R->tracer;
#endif

    }
#endif


#endif
}

double local_get_dynvisc_coefficient(double x, double y, double z, double rho, double press)
{
#ifdef VISCOSITY
#ifdef GLOBAL_VISCOSITY
  return All.dyn_visc;
#endif
#ifdef USE_KINEMATIC_VISCOSITY
  return All.KinematicViscosity * rho;
#endif
#ifdef ALPHA_VISCOSITY
  return get_alpha_viscosity(x, y, z, rho, press) * rho;
#endif
#endif
  return 0;
}


double local_get_bulkvisc_coefficient(double x, double y, double z, double rho, double press)
{
#ifdef VISCOSITY
#ifdef GLOBAL_VISCOSITY
  return All.bulk_visc;
#endif
#ifdef USE_KINEMATIC_VISCOSITY
  return 0;
#endif
#ifdef ALPHA_VISCOSITY
  return 0;
#endif
#endif
  return 0;
}

void face_get_viscous_fluxes(struct state_face *st_face, struct fluxes *flux, struct geometry *geom, double dyn_visc, double bulk_visc)
/*Once we have the velocity and velocity gradients (regardless of how we
obtained them, we need to project them onto the normals of each face)*/
/*The effect of moving boundaries is already taken care of in the computation
of the advectve fluxes*/
{
#ifdef VISCOSITY
  double div_vel = st_face->vel_grad[0][0] + st_face->vel_grad[1][1] + st_face->vel_grad[2][2];

  double fac1 = 4.0 / 3.0 * dyn_visc;
  double fac2 = 2.0 / 3.0 * dyn_visc;
  double fac3 = div_vel * bulk_visc;

  flux->momentum[0] +=
    (geom->nx *
     (fac1 * st_face->vel_grad[0][0] - fac2 * (st_face->vel_grad[1][1] + st_face->vel_grad[2][2]) + fac3) +
     geom->ny * dyn_visc * (st_face->vel_grad[0][1] + st_face->vel_grad[1][0]) + geom->nz * dyn_visc * (st_face->vel_grad[0][2] + st_face->vel_grad[2][0]));


  flux->momentum[1] +=
    (geom->ny *
     (fac1 * st_face->vel_grad[1][1] - fac2 * (st_face->vel_grad[0][0] + st_face->vel_grad[2][2]) + fac3) +
     geom->nx * dyn_visc * (st_face->vel_grad[0][1] + st_face->vel_grad[1][0]) + geom->nz * dyn_visc * (st_face->vel_grad[1][2] + st_face->vel_grad[2][1]));


  flux->momentum[2] +=
    (geom->nz *
     (fac1 * st_face->vel_grad[2][2] - fac2 * (st_face->vel_grad[0][0] + st_face->vel_grad[1][1]) + fac3) +
     geom->nx * dyn_visc * (st_face->vel_grad[0][2] + st_face->vel_grad[2][0]) + geom->ny * dyn_visc * (st_face->vel_grad[1][2] + st_face->vel_grad[2][1]));


  flux->energy +=
    (geom->nx *
     (st_face->velx *
      (fac1 * st_face->vel_grad[0][0] - fac2 * (st_face->vel_grad[1][1] + st_face->vel_grad[2][2]) + fac3) +
      st_face->vely * dyn_visc * (st_face->vel_grad[1][0] + st_face->vel_grad[0][1]) +
      st_face->velz * dyn_visc * (st_face->vel_grad[2][0] + st_face->vel_grad[0][2])) +
     geom->ny * (st_face->vely *
                 (fac1 * st_face->vel_grad[1][1] -
                  fac2 * (st_face->vel_grad[0][0] + st_face->vel_grad[2][2]) + fac3) +
                 st_face->velx * dyn_visc * (st_face->vel_grad[1][0] + st_face->vel_grad[0][1]) +
                 st_face->velz * dyn_visc * (st_face->vel_grad[1][2] + st_face->vel_grad[2][1])) +
     geom->nz * (st_face->velz *
                 (fac1 * st_face->vel_grad[2][2] -
                  fac2 * (st_face->vel_grad[0][0] + st_face->vel_grad[1][1]) + fac3) +
                 st_face->velx * dyn_visc * (st_face->vel_grad[2][0] + st_face->vel_grad[0][2]) + st_face->vely * dyn_visc * (st_face->vel_grad[2][1] + st_face->vel_grad[1][2])));


  CPU_Step[CPU_VISCOUS_FLUXES] += measure_time();
#endif
}



void face_get_conduction_fluxes(struct state_face *st_face, struct fluxes *flux, struct geometry *geom)
{
#ifdef THERMAL_CONDUCTION
  flux->energy += -All.ThermalConductivity * (st_face->dTemp[0] * geom->nx + st_face->dTemp[1] * geom->ny + st_face->dTemp[2] * geom->nz);
#endif
}

void face_get_scalar_diffusion_fluxes(struct state_face *st_face, struct fluxes *flux, struct geometry *geom)
{
#ifdef TRACER_DIFFUSION

  flux->tracer = All.TracerDiffusivity * (st_face->dConservedTracer[0] * geom->nx + st_face->dConservedTracer[1] * geom->ny + st_face->dConservedTracer[2] * geom->nz);

#endif
}

void face_extrapolate_viscous_kick(struct state *st, double dyn_visc, double bulk_visc)
{
#if defined(VISCOSITY) && defined(SECOND_DERIVATIVES)
/* here we kick the state velocity of the cell by applying the
viscous forces as they appear in the Navier-Stokes equations i.e.
as a source term */

  double laplacian_vel[3], grad_div_vel[3];
  struct hessian_data *hessian;


  /*  dyn_visc = st->dyn_visc;
     bulk_visc = st->bulk_visc; */

  hessian = st->hessian;

  laplacian_vel[0] = hessian->ddvelx[0][0] + hessian->ddvelx[1][1] + hessian->ddvelx[2][2];
  laplacian_vel[1] = hessian->ddvely[0][0] + hessian->ddvely[1][1] + hessian->ddvely[2][2];
  laplacian_vel[2] = hessian->ddvelz[0][0] + hessian->ddvelz[1][1] + hessian->ddvelz[2][2];

  grad_div_vel[0] = hessian->ddvelx[0][0] + hessian->ddvely[1][0] + hessian->ddvelz[2][0];
  grad_div_vel[1] = hessian->ddvelx[0][1] + hessian->ddvely[1][1] + hessian->ddvelz[2][1];
  grad_div_vel[2] = hessian->ddvelx[0][2] + hessian->ddvely[1][2] + hessian->ddvelz[2][2];


  st->velx += st->dt_half * (dyn_visc * laplacian_vel[0] + (bulk_visc + dyn_visc / 3.0 * grad_div_vel[0])) / st->rho;

  st->vely += st->dt_half * (dyn_visc * laplacian_vel[1] + (bulk_visc + dyn_visc / 3.0 * grad_div_vel[1])) / st->rho;

  st->velz += st->dt_half * (dyn_visc * laplacian_vel[2] + (bulk_visc + dyn_visc / 3.0 * grad_div_vel[2])) / st->rho;
#endif
}


double get_alpha_viscosity(double x, double y, double z, double rho, double press)
{
#if defined(VISCOSITY) && defined(ALPHA_VISCOSITY)
  if(ALPHA_VISCOSITY == 1)      //alpha-viscosity based on pressure
    return All.AlphaCoefficient * press * rho;
  else if(ALPHA_VISCOSITY == 2) // alpha-viscosity based on position
    {
#ifdef CIRCUMSTELLAR
      return get_circumstellar_alpha_viscosity(x, y, z, rho, press);
#else
      return 0;
#endif
    }
#endif
  return 0;
}
