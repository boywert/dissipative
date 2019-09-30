/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/riemann_hllc_special_relativity.c
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
  // This is HLLC like Mignone & Bodo 2005, but with their bugs corrected, and expressed in our conserved variables

  // we need the coordinate velocities in the lab frame. However, Arepo transformed to the face frame with a Galeian trafo
  // thus we just revert it here

  // When full GR we need to reconsider this issue !!!

  st_L->velx += vel_face[0];  // this is what is handed over as vel_face_turned
  st_L->vely += vel_face[1];
  st_L->velz += vel_face[2];
  st_R->velx += vel_face[0];
  st_R->vely += vel_face[1];
  st_R->velz += vel_face[2];

  double face_velx =  vel_face[0];

  double lambda_L[3], lambda_R[3];
  
  get_lambda( st_L, lambda_L );
  get_lambda( st_R, lambda_R );
  
  double a_L = dmin( lambda_L[2], lambda_R[2] );  // check that this is exactly what we need for HLLC - yes wo checking zero (that's important)
  double a_R = dmax( lambda_L[1], lambda_R[1] );
  
  double d_L, sx_L, sy_L, sz_L, t_L;
  get_conserved_quantities( st_L, &d_L, &sx_L, &sy_L, &sz_L, &t_L );
  double d_R, sx_R, sy_R, sz_R, t_R;
  get_conserved_quantities( st_R, &d_R, &sx_R, &sy_R, &sz_R, &t_R );
  
  
  struct fluxes flux_L, flux_R;
  get_fluxes_from_state( st_L, d_L, sx_L, sy_L, sz_L, t_L, &flux_L );
  get_fluxes_from_state( st_R, d_R, sx_R, sy_R, sz_R, t_R, &flux_R );

  // set up quadratic eq for velocity of contact disc.
  
  // to this end we use the above HLLE fluxes and states  - here we omit dividing by 1/(a_R-a_L)
  
  double F_t_hll  = (a_R * flux_L.energy - a_L * flux_R.energy + a_R * a_L * (t_R - t_L));
  double F_d_hll  = (a_R * flux_L.mass - a_L * flux_R.mass + a_R * a_L * (d_R - d_L));
  double F_sx_hll = (a_R * flux_L.momentum[0] - a_L * flux_R.momentum[0] + a_R * a_L * (sx_R - sx_L));
  
  double U_t_hll  = (a_R * t_R - a_L * t_L + flux_L.energy - flux_R.energy);
  double U_d_hll  = (a_R * d_R - a_L * d_L + flux_L.mass - flux_R.mass);
  double U_sx_hll = (a_R * sx_R - a_L * sx_L + flux_L.momentum[0] - flux_R.momentum[0]);
  
  double qa = F_t_hll + F_d_hll;
  double qb = - U_t_hll - U_d_hll - F_sx_hll;
  double qc = U_sx_hll;
  
  double lamstar = 0.5 * ( - qb - sqrt(qb*qb - 4.*qa*qc) ) / qa;  // only the minus sign yields physically acceptable solution
  
  // maybe need to check for < eps instead of exact zero...
  if(qa == 0)
    lamstar = 0;
  
  // now we have the velocity of the contact disc. and compute the starred pressure (Mignone's error corrected, my energy formulation, for left state, which should equal right state)
  
  double pstar = ( lamstar * ( a_L * (t_L + d_L) - sx_L) - (sx_L * (a_L - st_L->velx) - st_L->press ) ) / ( 1. - a_L * lamstar);
  
  double dstar, sxstar, systar, szstar, tstar;
  
  // check in which region the face lies
  
#ifdef SR_HLLC_ZERO_COMPVEL

  double comp_vel =  0.0;  // but this one is apparently numerically more stable

#else

  double comp_vel =  face_velx;   // one should actually choose this velocity

#endif

  if ( a_L >= comp_vel ) { // F_L,  one could use get_fluxes here, but below we need lamstar

    /*    flux->mass = d_L * st_L->velx;
    flux->momentum[0] = sx_L * st_L->velx + st_L->press;
    flux->momentum[1] = sy_L * st_L->velx;
    flux->momentum[2] = sz_L * st_L->velx;
    flux->energy = sx_L - d_L * st_L->velx;

    flux->mass        =  flux->mass         - face_velx * d_L;
    flux->momentum[0] =  flux->momentum[0]  - face_velx * sx_L;
    flux->momentum[1] =  flux->momentum[1]  - face_velx * sy_L;
    flux->momentum[2] =  flux->momentum[2]  - face_velx * sz_L;
    flux->energy      =  flux->energy       - face_velx * t_L; */

    // numerically better ?? implementation

    double vdiff  = st_L->velx - face_velx;

    flux->mass        = d_L * vdiff;
    flux->momentum[0] = sx_L* vdiff + st_L->press;
    flux->momentum[1] = sy_L* vdiff; 
    flux->momentum[2] = sz_L* vdiff; 
    flux->energy      = sx_L - d_L * st_L->velx - face_velx * t_L;

} else if ( a_R <= comp_vel ) { // F_R

    /*   flux->mass = d_R * st_R->velx;
    flux->momentum[0] = sx_R * st_R->velx + st_R->press;
    flux->momentum[1] = sy_R * st_R->velx;
    flux->momentum[2] = sz_R * st_R->velx;
    flux->energy = sx_R - d_R * st_R->velx;

    flux->mass        =  flux->mass         - face_velx * d_R;
    flux->momentum[0] =  flux->momentum[0]  - face_velx * sx_R;
    flux->momentum[1] =  flux->momentum[1]  - face_velx * sy_R;
    flux->momentum[2] =  flux->momentum[2]  - face_velx * sz_R;
    flux->energy      =  flux->energy       - face_velx * t_R;*/

    // numerically better ?? implementation

    double vdiff  = st_R->velx - face_velx;

    flux->mass        = d_R * vdiff;
    flux->momentum[0] = sx_R* vdiff + st_R->press;
    flux->momentum[1] = sy_R* vdiff; 
    flux->momentum[2] = sz_R* vdiff; 
    flux->energy      = sx_R - d_R * st_R->velx - face_velx * t_R;
 
} else if (lamstar >= comp_vel ) {

    dstar  = d_L * (a_L- st_L->velx) / (a_L-lamstar);
    sxstar = (sx_L * (a_L - st_L->velx) - st_L->press + pstar ) / (a_L-lamstar);
    systar = sy_L * (a_L- st_L->velx) / (a_L-lamstar);
    szstar = sz_L * (a_L- st_L->velx) / (a_L-lamstar);
    tstar  = (t_L * (a_L- st_L->velx) + pstar * lamstar - st_L->press * st_L->velx ) / (a_L-lamstar);
    // could be also expressed as tstar = sxstar/lamstar - dstar - pstar
    //double tstar = sxstar / lamstar - dstar - pstar
    
    /*flux->mass = dstar * lamstar;
    flux->momentum[0] = sxstar * lamstar + pstar;
    flux->momentum[1] = systar * lamstar;
    flux->momentum[2] = szstar * lamstar;
    flux->energy = sxstar - dstar * lamstar;

    flux->mass        =  flux->mass         - face_velx * dstar;
    flux->momentum[0] =  flux->momentum[0]  - face_velx * sxstar;
    flux->momentum[1] =  flux->momentum[1]  - face_velx * systar;
    flux->momentum[2] =  flux->momentum[2]  - face_velx * szstar;
    flux->energy      =  flux->energy       - face_velx * tstar;*/

    // numerically better ?? implementation

    double vdiff  = lamstar - face_velx;

    flux->mass        = dstar * vdiff;
    flux->momentum[0] = sxstar * vdiff + pstar;
    flux->momentum[1] = systar * vdiff; 
    flux->momentum[2] = szstar * vdiff; 
    flux->energy      = sxstar - dstar * lamstar - face_velx * tstar;

  } else {
    dstar  = d_R * (a_R- st_R->velx) / (a_R-lamstar);
    sxstar = (sx_R * (a_R - st_R->velx) - st_R->press + pstar ) / (a_R-lamstar);
    systar = sy_R * (a_R- st_R->velx) / (a_R-lamstar);
    szstar = sz_R * (a_R- st_R->velx) / (a_R-lamstar);
    tstar  = (t_R * (a_R- st_R->velx) + pstar * lamstar - st_R->press * st_R->velx ) / (a_R-lamstar);
    // could be also expressed as tstar = sxstar/lamstar - dstar - pstar
    //double tstar = sxstar / lamstar - dstar - pstar
    
    /*flux->mass = dstar * lamstar;
    flux->momentum[0] = sxstar * lamstar + pstar;
    flux->momentum[1] = systar * lamstar;
    flux->momentum[2] = szstar * lamstar;
    flux->energy = sxstar - dstar * lamstar;

    flux->mass        =  flux->mass         - face_velx * dstar;
    flux->momentum[0] =  flux->momentum[0]  - face_velx * sxstar;
    flux->momentum[1] =  flux->momentum[1]  - face_velx * systar;
    flux->momentum[2] =  flux->momentum[2]  - face_velx * szstar;
    flux->energy      =  flux->energy       - face_velx * tstar;*/

    // numerically better ?? implementation

    double vdiff  = lamstar - face_velx;

    flux->mass        = dstar * vdiff;
    flux->momentum[0] = sxstar * vdiff + pstar;
    flux->momentum[1] = systar * vdiff; 
    flux->momentum[2] = szstar * vdiff; 
    flux->energy      = sxstar - dstar * lamstar - face_velx * tstar;

  }
  // if (check that no other conditions can be fulfilled) actually, one of the conditions must be true because a_L < lamstar < a_R
  
  if ( (lamstar < a_L || lamstar > a_R) && a_L < 0. && a_R > 0. ) { // probably this is not necessary
    printf( "Problems with weird eigenvalues/velocities in relativistic HLLC: %g, %g, %g, stopping.", a_L, lamstar, a_R );
    terminate( "YEAHYEAHYEAH" );
  }
    
  /* For the moving mesh we probably need the states in the starred region if the fae is moving within the starred region
     So we should compute the starred states explicitely as above and save them somewhere
     Or we immeidatedly apply the coorection of the flux due to the moving face but then we need it's velocity here
   
   */
    
  st_face->press = 1;
  //st_face->press = pstar;  // p at contact disc.; if this is returned it may cause terminate although the code runs fine actually
  
  return st_face->press;
}

#else

double godunov_flux_3d_hlle_special_relativity(struct state *st_L, struct state *st_R, double *vel_face, struct state_face *st_face, struct fluxes *flux)
{
  // This is HLLC like Mignone & Bodo 2005, but with their bugs corrected, and expressed in our conserved variables
  double lambda_L[3], lambda_R[3];
  
  get_lambda( st_L, lambda_L );
  get_lambda( st_R, lambda_R );
  
  double a_L = dmin( lambda_L[2], lambda_R[2] );  // check that this is exactly what we need for HLLC - yes wo checking zero (that's important)
  double a_R = dmax( lambda_L[1], lambda_R[1] );
  
  double d_L, sx_L, sy_L, sz_L, t_L;
  get_conserved_quantities( st_L, &d_L, &sx_L, &sy_L, &sz_L, &t_L );
  double d_R, sx_R, sy_R, sz_R, t_R;
  get_conserved_quantities( st_R, &d_R, &sx_R, &sy_R, &sz_R, &t_R );
  
  
  struct fluxes flux_L, flux_R;
  get_fluxes_from_state( st_L, d_L, sx_L, sy_L, sz_L, t_L, &flux_L );
  get_fluxes_from_state( st_R, d_R, sx_R, sy_R, sz_R, t_R, &flux_R );

  // set up quadratic eq for velocity of contact disc.
  
  // to this end we use the above HLLE fluxes and states  - here we omit dividing by 1/(a_R-a_L)
  
  double F_t_hll  = (a_R * flux_L.energy - a_L * flux_R.energy + a_R * a_L * (t_R - t_L));
  double F_d_hll  = (a_R * flux_L.mass - a_L * flux_R.mass + a_R * a_L * (d_R - d_L));
  double F_sx_hll = (a_R * flux_L.momentum[0] - a_L * flux_R.momentum[0] + a_R * a_L * (sx_R - sx_L));
  
  double U_t_hll  = (a_R * t_R - a_L * t_L + flux_L.energy - flux_R.energy);
  double U_d_hll  = (a_R * d_R - a_L * d_L + flux_L.mass - flux_R.mass);
  double U_sx_hll = (a_R * sx_R - a_L * sx_L + flux_L.momentum[0] - flux_R.momentum[0]);
  
  double qa = F_t_hll + F_d_hll;
  double qb = - U_t_hll - U_d_hll - F_sx_hll;
  double qc = U_sx_hll;
  
  double lamstar = 0.5 * ( - qb - sqrt(qb*qb - 4.*qa*qc) ) / qa;  // only the minus sign yields physically acceptable solution
  
  // maybe need to check for < eps instead of exact zero...
  if(qa == 0)
    lamstar = 0;
  
  // now we have the velocity of the contact disc. and compute the starred pressure (Mignone's error corrected, my energy formulation, for left state, which should equal right state)
  
  double pstar = ( lamstar * ( a_L * (t_L + d_L) - sx_L) - (sx_L * (a_L - st_L->velx) - st_L->press ) ) / ( 1. - a_L * lamstar);
  
  double dstar, sxstar, systar, szstar, tstar;
  
  // check in which region the face lies
  
  if ( a_L >= 0. ) { // F_L,  one could use get_fluxes here, but below we need lamstar
    flux->mass = d_L * st_L->velx;
    flux->momentum[0] = sx_L * st_L->velx + st_L->press;
    flux->momentum[1] = sy_L * st_L->velx;
    flux->momentum[2] = sz_L * st_L->velx;
    flux->energy = sx_L - d_L * st_L->velx;
  } else if ( a_R <= 0.) { // F_R
    flux->mass = d_R * st_R->velx;
    flux->momentum[0] = sx_R * st_R->velx + st_R->press;
    flux->momentum[1] = sy_R * st_R->velx;
    flux->momentum[2] = sz_R * st_R->velx;
    flux->energy = sx_R - d_R * st_R->velx;
  } else if (lamstar >= 0. ) {
    dstar  = d_L * (a_L- st_L->velx) / (a_L-lamstar);
    sxstar = (sx_L * (a_L - st_L->velx) - st_L->press + pstar ) / (a_L-lamstar);
    systar = sy_L * (a_L- st_L->velx) / (a_L-lamstar);
    szstar = sz_L * (a_L- st_L->velx) / (a_L-lamstar);
    tstar  = (t_L * (a_L- st_L->velx) + pstar * lamstar - st_L->press * st_L->velx ) / (a_L-lamstar);
    // could be also expressed as tstar = sxstar/lamstar - dstar - pstar
    //double tstar = sxstar / lamstar - dstar - pstar
    
    flux->mass = dstar * lamstar;
    flux->momentum[0] = sxstar * lamstar + pstar;
    flux->momentum[1] = systar * lamstar;
    flux->momentum[2] = szstar * lamstar;
    flux->energy = sxstar - dstar * lamstar;
  } else {
    dstar  = d_R * (a_R- st_R->velx) / (a_R-lamstar);
    sxstar = (sx_R * (a_R - st_R->velx) - st_R->press + pstar ) / (a_R-lamstar);
    systar = sy_R * (a_R- st_R->velx) / (a_R-lamstar);
    szstar = sz_R * (a_R- st_R->velx) / (a_R-lamstar);
    tstar  = (t_R * (a_R- st_R->velx) + pstar * lamstar - st_R->press * st_R->velx ) / (a_R-lamstar);
    // could be also expressed as tstar = sxstar/lamstar - dstar - pstar
    //double tstar = sxstar / lamstar - dstar - pstar
    
    flux->mass = dstar * lamstar;
    flux->momentum[0] = sxstar * lamstar + pstar;
    flux->momentum[1] = systar * lamstar;
    flux->momentum[2] = szstar * lamstar;
    flux->energy = sxstar - dstar * lamstar;
  }
  // if (check that no other conditions can be fulfilled) actually, one of the conditions must be true because a_L < lamstar < a_R
  
  if ( (lamstar < a_L || lamstar > a_R) && a_L < 0. && a_R > 0. ) { // probably this is not necessary
    printf( "Problems with weird eigenvalues/velocities in relativistic HLLC: %g, %g, %g, stopping.", a_L, lamstar, a_R );
    terminate( "YEAHYEAHYEAH" );
  }
    
  /* For the moving mesh we probably need the states in the starred region if the fae is moving within the starred region
     So we should compute the starred states explicitely as above and save them somewhere
     Or we immeidatedly apply the coorection of the flux due to the moving face but then we need it's velocity here
   
   */
    
  //st_face->press = 1;
  st_face->press = pstar;  // p at contact disc. 
  
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
