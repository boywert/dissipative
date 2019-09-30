/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/riemann_hlle_general_relativity.c
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
  // This is HLLE here i need to think for GR: the fluxes are ok, but the term from the faces needs to be schecked
    
  double x = facex;
  double y = facey;
  double z = facez;
    
  // metric quantities in cfc - just declaration
  
  double alp, psi, bx, by, bz;

  get_metric_general_relativity(x,y,z,&alp,&psi,&bx,&by,&bz); 
  
  double bx_turned_to_xface = bx * geom->nx + by * geom->ny + bz * geom->nz;
    
    
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

  printf("do an eos call before\n");
  terminate("eos call required");
  
  double lambda_L[3], lambda_R[3];

  get_lambda( st_L, lambda_L, alp, psi, bx_turned_to_xface );
  get_lambda( st_R, lambda_R, alp, psi, bx_turned_to_xface );
  
  double a_L = dmin( 0., dmin( lambda_L[2], lambda_R[2] ) );
  double a_R = dmax( 0., dmax( lambda_L[1], lambda_R[1] ) );
  
  // below we need the real signal velocities - that's important

  double sa_L = dmin( lambda_L[2], lambda_R[2] );
  double sa_R = dmax( lambda_L[1], lambda_R[1] );

  double d_L, sx_L, sy_L, sz_L, t_L;
  get_conserved_quantities( st_L, &d_L, &sx_L, &sy_L, &sz_L, &t_L, psi );
  double d_R, sx_R, sy_R, sz_R, t_R;
  get_conserved_quantities( st_R, &d_R, &sx_R, &sy_R, &sz_R, &t_R, psi );
  
  struct fluxes flux_L, flux_R;
  get_fluxes_from_state( st_L, d_L, sx_L, sy_L, sz_L, t_L, &flux_L, alp, psi, bx_turned_to_xface );
  get_fluxes_from_state( st_R, d_R, sx_R, sy_R, sz_R, t_R, &flux_R, alp, psi, bx_turned_to_xface );
  
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

  //todo: be careful which velocity is used for faces

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

double godunov_flux_3d_hlle_general_relativity(struct state *st_L, struct state *st_R, double *vel_face, struct state_face *st_face, struct fluxes *flux, struct geometry *geom, double facex, double facey, double facez)
{

    
  double x = facex;
  double y = facey;
  double z = facez;
    
  // metric quantities in cfc - just declaration
  
  double alp, psi, bx, by, bz;

  get_metric_general_relativity(x,y,z,&alp,&psi,&bx,&by,&bz); 
  
  double bx_turned_to_xface = bx * geom->nx + by * geom->ny + bz * geom->nz;   

  // TODO: implemnet here general eos!!!!

  st_L->press=(GAMMA-1.)*st_L->utherm*st_L->rho;  // this is just to check if there is a difference
  st_R->press=(GAMMA-1.)*st_R->utherm*st_R->rho;

  if (st_L->press < 100.*pow(st_L->rho,GAMMA)) 
    {
    st_L->press   = 100.*pow(st_L->rho,GAMMA);
    st_L->utherm  = st_L->press / st_L->rho / (GAMMA-1.);
    }
  if (st_R->press < 100.*pow(st_R->rho,GAMMA))
    {
    st_R->press=100.*pow(st_R->rho,GAMMA);
    st_R->utherm  = st_R->press/ st_R->rho / (GAMMA-1.);
    }

  //reconstruction via Pressure
  //  st_L->utherm=st_L->press/st_L->rho/(GAMMA-1.);
  //  st_R->utherm=st_R->press/st_R->rho/(GAMMA-1.);


  // use extracted values from my grid code

  // values for the orginal flux 
  /* st_R->rho=1.100601261920919E-003;
   st_L->rho=1.084848585078480E-003;

   st_R->velz=1.201064353895558E-003;
    st_L->velz=1.174269068496179E-003;

    st_R->vely=1.201064353895576E-003;
    st_L->vely=1.174269068496182E-003;

    st_R->velx=9.013679755931743E-004;
    st_L->velx=1.174269068496197E-003;

    st_R->press=1.288410674224794E-004;
    st_L->press=1.256082036076200E-004;

    st_R->utherm=0.117064255584814;
    st_L->utherm=0.115784087600145;
    
    alp=0.680721267064839;
    psi=1.18671726588024;
    bx_turned_to_xface=0.0;*/

  /*st_R->rho=1.192464092887163E-003;
  st_L->rho=1.165723042357372E-003;
  st_R->velz=-7.197496753939645E-004;
  st_L->velz=-8.049052282711178E-004;
  st_R->vely=-7.197496753939315E-004;
  st_L->vely=-8.049052282711173E-004;
  st_R->velx=-5.827880229619364E-004;
  st_L->velx=-8.049052282711034E-004;
  st_R->press=1.297157302847943E-004;
  st_L->press=1.260612824961022E-004;
  st_R->utherm=0.108779569178247;
  st_L->utherm=0.108139993733997;

  alp=0.680721267064839;
  psi=1.18671726588024;
  bx_turned_to_xface=0.0;*/


  // This is HLLE for GR done
  double lambda_L[3], lambda_R[3];
  
  get_lambda( st_L, lambda_L, alp, psi, bx_turned_to_xface  );
  get_lambda( st_R, lambda_R, alp, psi, bx_turned_to_xface  );
  
  double a_L = dmin( 0., dmin( lambda_L[2], lambda_R[2] ) );
  double a_R = dmax( 0., dmax( lambda_L[1], lambda_R[1] ) );
  
  double d_L, sx_L, sy_L, sz_L, t_L;
  get_conserved_quantities( st_L, &d_L, &sx_L, &sy_L, &sz_L, &t_L, psi );
  double d_R, sx_R, sy_R, sz_R, t_R;
  get_conserved_quantities( st_R, &d_R, &sx_R, &sy_R, &sz_R, &t_R, psi );
  
  struct fluxes flux_L, flux_R;
  get_fluxes_from_state( st_L, d_L, sx_L, sy_L, sz_L, t_L, &flux_L, alp, psi, bx_turned_to_xface  );
  get_fluxes_from_state( st_R, d_R, sx_R, sy_R, sz_R, t_R, &flux_R, alp, psi, bx_turned_to_xface  );
  
  double a = 1. / (a_R - a_L);
  flux->mass = a * (a_R * flux_L.mass - a_L * flux_R.mass + a_R * a_L * (d_R - d_L));
  flux->momentum[0] = a * (a_R * flux_L.momentum[0] - a_L * flux_R.momentum[0] + a_R * a_L * (sx_R - sx_L));
  flux->momentum[1] = a * (a_R * flux_L.momentum[1] - a_L * flux_R.momentum[1] + a_R * a_L * (sy_R - sy_L));
  flux->momentum[2] = a * (a_R * flux_L.momentum[2] - a_L * flux_R.momentum[2] + a_R * a_L * (sz_R - sz_L));
  flux->energy = a * (a_R * flux_L.energy - a_L * flux_R.energy + a_R * a_L * (t_R - t_L));

  //printf("sx %g %g\n",sx_R,sx_L);
  //printf("taus %g %g\n",t_R,t_L);
  //printf("fluex %g %g %g %g %g\n",flux->mass,flux->momentum[0],flux->momentum[1],flux->momentum[2],flux->energy);
  //terminate("extracted values\n");
  
   if(!gsl_finite(flux->energy))
   {
       printf("GUde, here lareday bulshit staic \n");
       printf("\n");
       printf("as %g %g %g\n",a,a_R,a_L);      
       printf("fluxes %g %g %g %g\n",flux_L.energy,flux_R.energy,t_R , t_L);
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
  
  //    double det = pow(psi,6.);
//    
//    double coord_vx = st->velx - bx_turned_to_xface / alp;
//    
//    flux->mass =  det * alp * d * coord_vx;
//    
//    flux->momentum[0] = det * alp * ( sx * coord_vx + st->press ); 
//  
//    flux->momentum[1] = det * alp * sy * coord_vx;
//  
//    flux->momentum[2] = det * alp * sz * coord_vx; 
//    
//    flux->energy = det * alp * ( t * coord_vx + st->press * st->velx );  // note that there are different formulations possible

 double det = pow(psi,6.);
 
 double coord_vx = st->velx - bx_turned_to_xface / alp;
 
 flux->mass =   alp * d * coord_vx;
 
 flux->momentum[0] =  alp * ( sx * coord_vx + det * st->press ); 

 flux->momentum[1] =  alp * sy * coord_vx;

 flux->momentum[2] =  alp * sz * coord_vx; 
 
 flux->energy =  alp * ( t * coord_vx + det * st->press * st->velx );  // note that there are different formulations possible -> probably not a big difference if cons are consistent with prims

  //printf("in HLLE ad %g %g %g %g %g %g\n",alp , t  ,bx_turned_to_xface, det , st->press , st->velx);

  double ps4 = psi*psi*psi*psi;
 
 flux->energy =  alp * ( sx / ps4 - d * st->velx ) - bx_turned_to_xface * t; 
  //  printf("in HLLE %g %g %g %g %g %g %g\n",alp , sx , ps4 , d , st->velx , bx_turned_to_xface , t);    
}

#endif
