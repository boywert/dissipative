/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/riemann_hlle_special_relativity.c
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

#include "allvars.h"
#include "proto.h"
#include <math.h>

#ifdef SPECIAL_RELATIVITY

static void get_lambda(struct state *st, double *lambda); // get signal speed
static void get_conserved_quantities( struct state *st, double *d, double *sx, double *sy, double *sz, double *t);
static void get_fluxes_from_state(struct state *st, double d, double sx, double sy, double sz, double t, struct fluxes *flux);

#ifndef VORONOI_STATIC_MESH

double godunov_flux_3d_hlle_special_relativity(struct state *st_L, struct state *st_R, double *vel_face, struct state_face *st_face, struct fluxes *flux)
{
  // This is HLLE

  // we need the coordinate velocities in the lab frame. However, Arepo transformed to the face frame with a Galeian trafo
  // thus we just revert it here

  // When full GR we need to reconsider this issue !!!

  st_L->velx += vel_face[0];   // this is what is handed over as vel_face_turned
  st_L->vely += vel_face[1];
  st_L->velz += vel_face[2];
  st_R->velx += vel_face[0];
  st_R->vely += vel_face[1];
  st_R->velz += vel_face[2];

  double face_velx =  vel_face[0];

  double lambda_L[3], lambda_R[3];

  get_lambda( st_L, lambda_L );
  get_lambda( st_R, lambda_R );
  
  double a_L = dmin( 0., dmin( lambda_L[2], lambda_R[2] ) );
  double a_R = dmax( 0., dmax( lambda_L[1], lambda_R[1] ) );
  
  // below we need the real signal velocities - that's important

  double sa_L = dmin( lambda_L[2], lambda_R[2] );
  double sa_R = dmax( lambda_L[1], lambda_R[1] );

  double d_L, sx_L, sy_L, sz_L, t_L;
  get_conserved_quantities( st_L, &d_L, &sx_L, &sy_L, &sz_L, &t_L );
  double d_R, sx_R, sy_R, sz_R, t_R;
  get_conserved_quantities( st_R, &d_R, &sx_R, &sy_R, &sz_R, &t_R );
  
  struct fluxes flux_L, flux_R;
  get_fluxes_from_state( st_L, d_L, sx_L, sy_L, sz_L, t_L, &flux_L );
  get_fluxes_from_state( st_R, d_R, sx_R, sy_R, sz_R, t_R, &flux_R );
  
  double a = 1. / (a_R - a_L);
  flux->mass = a * (a_R * flux_L.mass - a_L * flux_R.mass + a_R * a_L * (d_R - d_L));
  flux->momentum[0] = a * (a_R * flux_L.momentum[0] - a_L * flux_R.momentum[0] + a_R * a_L * (sx_R - sx_L));
  flux->momentum[1] = a * (a_R * flux_L.momentum[1] - a_L * flux_R.momentum[1] + a_R * a_L * (sy_R - sy_L));
  flux->momentum[2] = a * (a_R * flux_L.momentum[2] - a_L * flux_R.momentum[2] + a_R * a_L * (sz_R - sz_L));
  flux->energy = a * (a_R * flux_L.energy - a_L * flux_R.energy + a_R * a_L * (t_R - t_L));
  
  // compare face velocity to signal velocities

#ifdef SR_HLLC_ZERO_COMPVEL

  double comp_vel =  0.0;  // but this one is apparently numerically more stable

#else

  double comp_vel =  face_velx;   // one should actually choose this velocity

#endif


  if ( comp_vel >= sa_R ) {

    flux->mass        =  flux->mass         - face_velx * d_R;// could be done in a numerically better way
    flux->momentum[0] =  flux->momentum[0]  - face_velx * sx_R;
    flux->momentum[1] =  flux->momentum[1]  - face_velx * sy_R;
    flux->momentum[2] =  flux->momentum[2]  - face_velx * sz_R;
    flux->energy      =  flux->energy       - face_velx * t_R;

  }
  else if ( comp_vel <= sa_L ) {

    flux->mass        =  flux->mass         - face_velx * d_L;
    flux->momentum[0] =  flux->momentum[0]  - face_velx * sx_L;
    flux->momentum[1] =  flux->momentum[1]  - face_velx * sy_L;
    flux->momentum[2] =  flux->momentum[2]  - face_velx * sz_L;
    flux->energy      =  flux->energy       - face_velx * t_L;

  }
  else if ( comp_vel > sa_L && comp_vel < sa_R ) { // that means we are in the starred region

    double sa =  1. / (sa_R - sa_L);

    double U_d_hll  = sa * (sa_R * d_R - sa_L * d_L + flux_L.mass - flux_R.mass);
    double U_sx_hll = sa * (sa_R * sx_R - sa_L * sx_L + flux_L.momentum[0] - flux_R.momentum[0]);
    double U_sy_hll = sa * (sa_R * sy_R - sa_L * sy_L + flux_L.momentum[1] - flux_R.momentum[1]);
    double U_sz_hll = sa * (sa_R * sz_R - sa_L * sz_L + flux_L.momentum[2] - flux_R.momentum[2]);
    double U_t_hll  = sa * (sa_R * t_R - sa_L * t_L + flux_L.energy - flux_R.energy);

    flux->mass        =  flux->mass         - face_velx * U_d_hll;
    flux->momentum[0] =  flux->momentum[0]  - face_velx * U_sx_hll;
    flux->momentum[1] =  flux->momentum[1]  - face_velx * U_sy_hll;
    flux->momentum[2] =  flux->momentum[2]  - face_velx * U_sz_hll;
    flux->energy      =  flux->energy       - face_velx * U_t_hll;

  }
  else {

     printf( "Problems in relativistic moving HLLE with velocities: %g, %g, %g, stopping.", a_L, face_velx, a_R );
     terminate( "Mmh" );

  }

  st_face->press = 1;  
  return st_face->press;
}

#else

double godunov_flux_3d_hlle_special_relativity(struct state *st_L, struct state *st_R, double *vel_face, struct state_face *st_face, struct fluxes *flux)
{
  // This is HLLE
  double lambda_L[3], lambda_R[3];
  
  get_lambda( st_L, lambda_L );
  get_lambda( st_R, lambda_R );
  
  double a_L = dmin( 0., dmin( lambda_L[2], lambda_R[2] ) );
  double a_R = dmax( 0., dmax( lambda_L[1], lambda_R[1] ) );
  
  double d_L, sx_L, sy_L, sz_L, t_L;
  get_conserved_quantities( st_L, &d_L, &sx_L, &sy_L, &sz_L, &t_L );
  double d_R, sx_R, sy_R, sz_R, t_R;
  get_conserved_quantities( st_R, &d_R, &sx_R, &sy_R, &sz_R, &t_R );
  
  struct fluxes flux_L, flux_R;
  get_fluxes_from_state( st_L, d_L, sx_L, sy_L, sz_L, t_L, &flux_L );
  get_fluxes_from_state( st_R, d_R, sx_R, sy_R, sz_R, t_R, &flux_R );
  
  double a = 1. / (a_R - a_L);
  flux->mass = a * (a_R * flux_L.mass - a_L * flux_R.mass + a_R * a_L * (d_R - d_L));
  flux->momentum[0] = a * (a_R * flux_L.momentum[0] - a_L * flux_R.momentum[0] + a_R * a_L * (sx_R - sx_L));
  flux->momentum[1] = a * (a_R * flux_L.momentum[1] - a_L * flux_R.momentum[1] + a_R * a_L * (sy_R - sy_L));
  flux->momentum[2] = a * (a_R * flux_L.momentum[2] - a_L * flux_R.momentum[2] + a_R * a_L * (sz_R - sz_L));
  flux->energy = a * (a_R * flux_L.energy - a_L * flux_R.energy + a_R * a_L * (t_R - t_L));
  
  st_face->press = 1;  
  return st_face->press;
}

#endif

void get_lambda(struct state *st, double *lambda)
{
  double v2 = st->velx*st->velx + st->vely*st->vely + st->velz*st->velz;
  double enthalpy = 1. + st->utherm + st->press / st->rho;
  double csnd = sqrt( GAMMA * st->press / (st->rho * enthalpy) );
  
  double t1 = st->velx * (1. - csnd*csnd) / (1. - v2 * csnd*csnd);
  double t2 = csnd / (1. - v2 * csnd*csnd);
  double t3 = (1.-v2) * (1. - v2 * csnd*csnd - st->velx*st->velx * (1.-csnd*csnd));
  
  lambda[0] = st->velx;
  lambda[1] = t1 + t2 * sqrt(t3);
  lambda[2] = t1 - t2 * sqrt(t3);
}

void get_conserved_quantities( struct state *st, double *d, double *sx, double *sy, double *sz, double *t )
{
  double v2 = st->velx*st->velx + st->vely*st->vely + st->velz*st->velz;
  double enthalpy = 1. + st->utherm + st->press / st->rho;
  double w = 1. / sqrt( 1. - v2 ); // lorentz factor
  
  *d = w * st->rho;
  *sx = *d * enthalpy * w * st->velx;
  *sy = *d * enthalpy * w * st->vely;
  *sz = *d * enthalpy * w * st->velz;
  
  *t = *d * enthalpy * w - st->press - *d;
}

void get_fluxes_from_state(struct state *st, double d, double sx, double sy, double sz, double t, struct fluxes *flux)
{
  flux->mass = d * st->velx;
  flux->momentum[0] = sx * st->velx + st->press;
  flux->momentum[1] = sy * st->velx;
  flux->momentum[2] = sz * st->velx;
  flux->energy = sx - d * st->velx;
}

#endif
