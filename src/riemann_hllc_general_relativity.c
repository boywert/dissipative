/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/riemann_hllc_general_relativity.c
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

#ifdef GENERAL_RELATIVITY

static void get_lambda(struct state *st, double *lambda, double alp, double psi, double bx_turned_to_xface); // get signal speed
static void get_conserved_quantities( struct state *st, double *d, double *sx, double *sy, double *sz, double *t, double psi);
static void get_fluxes_from_state(struct state *st, double d, double sx, double sy, double sz, double t, struct fluxes *flux, double alp, double psi, double bx_turned_to_xface);

#ifndef VORONOI_STATIC_MESH

double godunov_flux_3d_hlle_general_relativity(struct state *st_L, struct state *st_R, double *vel_face, struct state_face *st_face, struct fluxes *flux, struct geometry *geom, double facex, double facey, double facez)
{
    
  // This is HLLC here i need to think for GR: the fluxes are ok, but the term from the faces needs to be schecked
    
  double x = facex;
  double y = facey;
  double z = facez;
    
  // metric quantities in cfc - just declaration
  
  double alp, psi, bx, by, bz;

  get_metric_general_relativity(x,y,z,&alp,&psi,&bx,&by,&bz); 
  
  double bx_turned_to_xface = bx * geom->nx + by * geom->ny + bz * geom->nz;
      
  double det = psi*psi*psi*psi*psi*psi;  
  
  double ps4 = psi*psi*psi*psi;   
    
  // we need the coordinate velocities in the lab frame. However, Arepo transformed to the face frame with a Galeian trafo
  // thus we just revert it here

  // When full GR we need to reconsider this issue !!! -> we should use the cooridnate velocity here, not the primitive!!!
    

  st_L->velx += vel_face[0];   // this is what is handed over as vel_face_turned
  st_L->vely += vel_face[1];
  st_L->velz += vel_face[2];
  st_R->velx += vel_face[0];
  st_R->vely += vel_face[1];
  st_R->velz += vel_face[2];

  double face_velx =  vel_face[0] - bx_turned_to_xface / alp;   //  <- here we should use coordinate velocity

  // but careful if we refine the velocity faces as coordinate velocities, for the moment this is correct since the vels contain the primitive velocity and NOT the coordinate velocity
  
  // This is HLLC for GR
  double lambda_L[3], lambda_R[3];
  
  get_lambda( st_L, lambda_L, alp, psi, bx_turned_to_xface  );
  get_lambda( st_R, lambda_R, alp, psi, bx_turned_to_xface  );
  
  double a_L =  dmin( lambda_L[2], lambda_R[2] );  // important not to include zero here as in HLLE
  double a_R =  dmax( lambda_L[1], lambda_R[1] );
  
  double d_L, sx_L, sy_L, sz_L, t_L;
  get_conserved_quantities( st_L, &d_L, &sx_L, &sy_L, &sz_L, &t_L, psi );
  double d_R, sx_R, sy_R, sz_R, t_R;
  get_conserved_quantities( st_R, &d_R, &sx_R, &sy_R, &sz_R, &t_R, psi );

  // coefficients for quadratic eqs.
  
  double vbar  =  st_L->velx - bx_turned_to_xface / alp;

  double cl    =  t_L * ( a_L - alp * vbar ) - alp * det * st_L->velx * st_L->press;

  double el    =  d_L * ( a_L - alp * vbar );

  double bl    =  sx_L * ( a_L - alp * vbar ) - alp * det * st_L->press;

  vbar     =  st_R->velx - bx_turned_to_xface / alp;
  
  double cr    =   t_R * ( a_R - alp * vbar ) - alp * det * st_R->velx * st_R->press;

  double er    =   d_R * ( a_R - alp * vbar );

  double br    =  sx_R * ( a_R - alp * vbar ) - alp * det * st_R->press;

  
  double qa  = ps4 * ( -(cl+el)*(a_R+bx_turned_to_xface) + (cr+er)*(a_L+bx_turned_to_xface) )  / alp;
  
  double qb  = ( cl+el -cr-er + bl*(a_R+bx_turned_to_xface)/alp -br*(a_L+bx_turned_to_xface)/alp  ) ;

  double qc  = (br - bl)/ps4;
  
  // solve quadratic eq. and catch certain cases
  
  double chi;
  
  if ( qb*qb - 4.*qa*qc >= 0.0 ) {
      chi      =  0.5 * ( - qb - sqrt(qb*qb - 4.*qa*qc) ) / qa;
  }
  else {
      chi      =  0.5 * - qb  / qa;
  }
  
  if ( fabs(qa) < 1.e-10 ) {  // hand-waving choice 
      chi = -qc/qb;
  }
  
  // coordinate velocity of contact disc.
  
  double lamstar = chi - bx_turned_to_xface/alp;
  
  // check:
  
  if (!gsl_finite(lamstar) || !gsl_finite(chi)) {
      
      printf("lamstar chi %g %g \n",lamstar,chi);
      printf("qs %g %g %g \n",qa,qb,qc);
      terminate(" stupid lamstar");
  }
  
  // compute pressure at contact disc. - L and R should be the same by construction but may differ numerically
      
  double pstarl = ( bl - ps4*chi*(cl+el) ) / ( ps4*chi*(a_L+bx_turned_to_xface)/alp - 1.0 ) / alp / det;

  double pstarr = ( br - ps4*chi*(cr+er) ) / ( ps4*chi*(a_R+bx_turned_to_xface)/alp - 1.0 ) / alp / det;
  
  // locate region of face and compute fluxes
  
  
#ifdef SR_HLLC_ZERO_COMPVEL

  double comp_vel =  0.0;  // but this one is apparently numerically more stable

#else

  double comp_vel =  face_velx;   // one should actually choose this velocity

#endif

  if ( a_L >= comp_vel ) {

      struct fluxes flux_L;
      get_fluxes_from_state( st_L, d_L, sx_L, sy_L, sz_L, t_L, &flux_L, alp, psi, bx_turned_to_xface  );
  
      flux->mass        =  flux_L.mass;
      flux->momentum[0] =  flux_L.momentum[0];
      flux->momentum[1] =  flux_L.momentum[1];
      flux->momentum[2] =  flux_L.momentum[2];
      flux->energy      =  flux_L.energy; 
      
  }
  else if ( a_R <= comp_vel) {
      
      struct fluxes flux_R;
      get_fluxes_from_state( st_R, d_R, sx_R, sy_R, sz_R, t_R, &flux_R, alp, psi, bx_turned_to_xface  );
  
      flux->mass        =  flux_R.mass;
      flux->momentum[0] =  flux_R.momentum[0];
      flux->momentum[1] =  flux_R.momentum[1];
      flux->momentum[2] =  flux_R.momentum[2];
      flux->energy      =  flux_R.energy;    
  }
  else if ( lamstar >= comp_vel ) {
      
      //computed starred quantities
      
      vbar    = st_L->velx - bx_turned_to_xface / alp;

      double pstar   = pstarl;
      
      double dstar   = el / (a_L-alp*lamstar);

      double sxstar  = (bl + alp*det*pstar ) / (a_L - alp*lamstar);

      double systar  = sy_L * (a_L-alp*vbar) / (a_L-alp*lamstar);

      double szstar  = sz_L * (a_L-alp*vbar) / (a_L-alp*lamstar);

      double taustar = (cl + alp*det*chi*pstar) / (a_L-alp*lamstar);
      
      flux->mass        = alp * dstar * lamstar;
      flux->momentum[0] = alp * sxstar * lamstar + alp * det*pstar;
      flux->momentum[1] = alp * systar * lamstar;
      flux->momentum[2] = alp * szstar * lamstar;
      flux->energy      = alp * taustar * lamstar + alp * det * chi * pstar;
      
  }
  else {
      //computed starred quantities
      
      vbar    = st_R->velx - bx_turned_to_xface / alp;

      double pstar   = pstarr;
      
      double dstar   = er / (a_R-alp*lamstar);

      double sxstar  = (br + alp*det*pstar ) / (a_R - alp*lamstar);

      double systar  = sy_R * (a_R-alp*vbar) / (a_R-alp*lamstar);

      double szstar  = sz_R * (a_R-alp*vbar) / (a_R-alp*lamstar);

      double taustar = (cr + alp*det*chi*pstar) / (a_R-alp*lamstar);
      
      flux->mass        = alp * dstar * lamstar;
      flux->momentum[0] = alp * sxstar * lamstar + alp * det*pstar;
      flux->momentum[1] = alp * systar * lamstar;
      flux->momentum[2] = alp * szstar * lamstar;
      flux->energy      = alp * taustar * lamstar + alp * det * chi * pstar   ;
      
  }
      
      
   if(!gsl_finite(flux->energy))
   {
       printf("GUde, here lareday bulshit HLLC moving\n");
       printf("\n");
       printf("as %g %g %g \n",a_R,a_L,lamstar);      
       //printf("fluxes %g %g %g %g\n",flux_L.energy,flux_R.energy,t_R , t_L);
       printf("\n");
       printf("lams %g %g %g %g\n",lambda_L[2], lambda_R[2],lambda_L[1], lambda_R[1]);
       printf("metric %g %g %g %g\n",alp,psi,bx_turned_to_xface,lambda_L[0]);
       printf("\n");       
       printf("state L %g %g %g %g %g %g \n",st_L->rho,st_L->velx,st_L->vely,st_L->velz,st_L->press,st_L->utherm);
       
       printf("\n");
       printf("state R %g %g %g %g %g %g \n",st_R->rho,st_R->velx,st_R->vely,st_R->velz,st_R->press,st_R->utherm);       
       printf("\n");
       printf("position %g %g %g\n",x,y,z);
       printf("\n");
       printf("cons L %g %g %g %g %g \n",d_L, sx_L, sy_L, sz_L, t_L);       
       printf("\n");      
       printf("cons R %g %g %g %g %g \n",d_R, sx_R, sy_R, sz_R, t_R); 
       printf("\n");      
   }
  
  
  st_face->press = 1;  
  return st_face->press;
}

#else

double godunov_flux_3d_hlle_general_relativity(struct state *st_L, struct state *st_R, double *vel_face, struct state_face *st_face, struct fluxes *flux, struct geometry *geom, double facex, double facey, double facez)
{

    
  double x = facex;
  double y = facey;
  double z = facez;
    
  // metric quantities in cfc - just declaration
  
  double alp, psi, bx, by, bz;

  get_metric_general_relativity(x,y,z,&alp,&psi,&bx,&by,&bz); 
  
  double bx_turned_to_xface = bx * geom->nx + by * geom->ny + bz * geom->nz;   
  
  double det = psi*psi*psi*psi*psi*psi;  
  
  double ps4 = psi*psi*psi*psi;  
    
  // This is HLLC for GR
  double lambda_L[3], lambda_R[3];
  
  get_lambda( st_L, lambda_L, alp, psi, bx_turned_to_xface  );
  get_lambda( st_R, lambda_R, alp, psi, bx_turned_to_xface  );
  
  double a_L =  dmin( lambda_L[2], lambda_R[2] );  // important not to include zero here as in HLLE
  double a_R =  dmax( lambda_L[1], lambda_R[1] );
  
  double d_L, sx_L, sy_L, sz_L, t_L;
  get_conserved_quantities( st_L, &d_L, &sx_L, &sy_L, &sz_L, &t_L, psi );
  double d_R, sx_R, sy_R, sz_R, t_R;
  get_conserved_quantities( st_R, &d_R, &sx_R, &sy_R, &sz_R, &t_R, psi );

  // coefficients for quadratic eqs.
  
  double vbar  =  st_L->velx - bx_turned_to_xface / alp;

  double cl    =  t_L * ( a_L - alp * vbar ) - alp * det * st_L->velx * st_L->press;

  double el    =  d_L * ( a_L - alp * vbar );

  double bl    =  sx_L * ( a_L - alp * vbar ) - alp * det * st_L->press;

  vbar     =  st_R->velx - bx_turned_to_xface / alp;
  
  double cr    =   t_R * ( a_R - alp * vbar ) - alp * det * st_R->velx * st_R->press;

  double er    =   d_R * ( a_R - alp * vbar );

  double br    =  sx_R * ( a_R - alp * vbar ) - alp * det * st_R->press;

  
  double qa  = ps4 * ( -(cl+el)*(a_R+bx_turned_to_xface) + (cr+er)*(a_L+bx_turned_to_xface) )  / alp;
  
  double qb  = ( cl+el -cr-er + bl*(a_R+bx_turned_to_xface)/alp -br*(a_L+bx_turned_to_xface)/alp  ) ;

  double qc  = (br - bl)/ps4;
  
  // solve quadratic eq. and catch certain cases
  
  double chi;
  
  if ( qb*qb - 4.*qa*qc >= 0.0 ) {
      chi      =  0.5 * ( - qb - sqrt(qb*qb - 4.*qa*qc) ) / qa;
  }
  else {
      chi      =  0.5 * - qb  / qa;
  }
  
  if ( fabs(qa) < 1.e-10 ) {  // hand-waving choice 
      chi = -qc/qb;
  }
  
  // coordinate velocity of contact disc.
  
  double lamstar = chi - bx_turned_to_xface/alp;
  
  // check:
  
  if (!gsl_finite(lamstar) || !gsl_finite(chi)) {
      
      printf("lamstar chi %g %g \n",lamstar,chi);
      printf("qs %g %g %g \n",qa,qb,qc);
      terminate(" stupid lamstar");
  }
  
  // compute pressure at contact disc. - L and R should be the same by construction but may differ numerically
      
  double pstarl = ( bl - ps4*chi*(cl+el) ) / ( ps4*chi*(a_L+bx_turned_to_xface)/alp - 1.0 ) / alp / det;

  double pstarr = ( br - ps4*chi*(cr+er) ) / ( ps4*chi*(a_R+bx_turned_to_xface)/alp - 1.0 ) / alp / det;
  
  // locate region of face and compute fluxes
  
  if ( a_L >= 0.0 ) {

      struct fluxes flux_L;
      get_fluxes_from_state( st_L, d_L, sx_L, sy_L, sz_L, t_L, &flux_L, alp, psi, bx_turned_to_xface  );
  
      flux->mass        =  flux_L.mass;
      flux->momentum[0] =  flux_L.momentum[0];
      flux->momentum[1] =  flux_L.momentum[1];
      flux->momentum[2] =  flux_L.momentum[2];
      flux->energy      =  flux_L.energy; 
      
  }
  else if ( a_R <= 0.0) {
      
      struct fluxes flux_R;
      get_fluxes_from_state( st_R, d_R, sx_R, sy_R, sz_R, t_R, &flux_R, alp, psi, bx_turned_to_xface  );
  
      flux->mass        =  flux_R.mass;
      flux->momentum[0] =  flux_R.momentum[0];
      flux->momentum[1] =  flux_R.momentum[1];
      flux->momentum[2] =  flux_R.momentum[2];
      flux->energy      =  flux_R.energy;    
  }
  else if ( lamstar >= 0.0 ) {
      
      //computed starred quantities
      
      vbar    = st_L->velx - bx_turned_to_xface / alp;

      double pstar   = pstarl;
      
      double dstar   = el / (a_L-alp*lamstar);

      double sxstar  = (bl + alp*det*pstar ) / (a_L - alp*lamstar);

      double systar  = sy_L * (a_L-alp*vbar) / (a_L-alp*lamstar);

      double szstar  = sz_L * (a_L-alp*vbar) / (a_L-alp*lamstar);

      double taustar = (cl + alp*det*chi*pstar) / (a_L-alp*lamstar);
      
      flux->mass        = alp * dstar * lamstar;
      flux->momentum[0] = alp * sxstar * lamstar + alp * det*pstar;
      flux->momentum[1] = alp * systar * lamstar;
      flux->momentum[2] = alp * szstar * lamstar;
      flux->energy      = alp * taustar * lamstar + alp * det * chi * pstar;
      
  }
  else {
      //computed starred quantities
      
      vbar    = st_R->velx - bx_turned_to_xface / alp;

      double pstar   = pstarr;
      
      double dstar   = er / (a_R-alp*lamstar);

      double sxstar  = (br + alp*det*pstar ) / (a_R - alp*lamstar);

      double systar  = sy_R * (a_R-alp*vbar) / (a_R-alp*lamstar);

      double szstar  = sz_R * (a_R-alp*vbar) / (a_R-alp*lamstar);

      double taustar = (cr + alp*det*chi*pstar) / (a_R-alp*lamstar);
      
      flux->mass        = alp * dstar * lamstar;
      flux->momentum[0] = alp * sxstar * lamstar + alp * det*pstar;
      flux->momentum[1] = alp * systar * lamstar;
      flux->momentum[2] = alp * szstar * lamstar;
      flux->energy      = alp * taustar * lamstar + alp * det * chi * pstar   ;
      
  }
      
      
   if(!gsl_finite(flux->energy))
   {
       printf("GUde, here lareday bulshit HLLC static\n");
       printf("\n");
       printf("as %g %g %g \n",a_R,a_L,lamstar);      
       //printf("fluxes %g %g %g %g\n",flux_L.energy,flux_R.energy,t_R , t_L);
       printf("\n");
       printf("lams %g %g %g %g\n",lambda_L[2], lambda_R[2],lambda_L[1], lambda_R[1]);
       printf("metric %g %g %g %g\n",alp,psi,bx_turned_to_xface,lambda_L[0]);
       printf("\n");       
       printf("state L %g %g %g %g %g %g \n",st_L->rho,st_L->velx,st_L->vely,st_L->velz,st_L->press,st_L->utherm);
       
       printf("\n");
       printf("state R %g %g %g %g %g %g \n",st_R->rho,st_R->velx,st_R->vely,st_R->velz,st_R->press,st_R->utherm);       
       printf("\n");
       printf("position %g %g %g\n",x,y,z);
       printf("\n");
       printf("cons L %g %g %g %g %g \n",d_L, sx_L, sy_L, sz_L, t_L);       
       printf("\n");      
       printf("cons R %g %g %g %g %g \n",d_R, sx_R, sy_R, sz_R, t_R); 
       printf("\n");      
   }
  
  
  st_face->press = 1;  
  return st_face->press;
}

#endif

void get_lambda(struct state *st, double *lambda, double alp, double psi, double bx_turned_to_xface )
{
  
  double ps4 = pow(psi,4.);
  
  double v2 = st->velx*st->velx + st->vely*st->vely + st->velz*st->velz;
  
  v2  = v2 * ps4; // here we need to be careful when the spatial metric is non-diagnonal
  
  double enthalpy = 1. + st->utherm + st->press / st->rho;
  double csnd = sqrt( GAMMA * st->press / (st->rho * enthalpy) );
  
  double t1 = alp * st->velx * (1.-csnd*csnd) / (1.-v2*csnd*csnd);
  
  double t2 = alp * csnd / (1.-v2*csnd*csnd);  // this is just the prefactor
  
  double gammaxx = 1./ps4;  // this is gamma^XX; for cfc quite simple
  
  double t3 = sqrt( (1.-v2) * ( gammaxx * (1.-v2*csnd*csnd) - st->velx * st->velx * (1.-csnd*csnd) ) );
  
  
  lambda[0] = alp * st->velx - bx_turned_to_xface; 
  
  lambda[1] = t1 + t2 * t3 - bx_turned_to_xface;

  lambda[2] = t1 - t2 * t3 - bx_turned_to_xface;

    
}

void get_conserved_quantities( struct state *st, double *d, double *sx, double *sy, double *sz, double *t, double psi )
{
  
  double det = pow(psi,6.);
    
  double ps4 = pow(psi,4.);
  
  
  double v2 = st->velx*st->velx + st->vely*st->vely + st->velz*st->velz;
  double enthalpy = 1. + st->utherm + st->press / st->rho;
  double w = 1. / sqrt( 1. - v2 * ps4 ); // lorentz factor
  
  *d = det * w * st->rho;
  *sx = *d * enthalpy * w * st->velx * ps4;
  *sy = *d * enthalpy * w * st->vely * ps4;
  *sz = *d * enthalpy * w * st->velz * ps4; // the last ps4 because v_j=ps4*v^j
    
    // here we need to be careful when the spatial metric is not diagonal any longer
  
  *t = *d * enthalpy * w - det * st->press - *d;
}

void get_fluxes_from_state(struct state *st, double d, double sx, double sy, double sz, double t, struct fluxes *flux, double alp, double psi, double bx_turned_to_xface)
{
  
  double det = pow(psi,6.);
  double ps4 = pow(psi,4.);

  double coord_vx = st->velx - bx_turned_to_xface / alp;
  
  flux->mass =   alp * d * coord_vx;
  
  flux->momentum[0] =  alp * ( sx * coord_vx + st->press * det ); 

  flux->momentum[1] =  alp * sy * coord_vx;

  flux->momentum[2] =  alp * sz * coord_vx; 
  
  //flux->energy = alp * ( t * coord_vx + st->press * st->velx * det );  // note that there are different formulations possible

  flux->energy = alp * (sx/ps4 - d*st->velx) - bx_turned_to_xface * t; // maybe better because tau less effective
}

#endif
