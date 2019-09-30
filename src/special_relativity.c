/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/special_relativity.c
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

void update_primitive_variables_special_relativity(struct particle_data *P, struct sph_particle_data *SphP, int i, struct pv_update_data *pvd)
{
  /* here we assume that the EoS is available in this routine either as table or analytically
   for the momemt this routine works only for analytic ideal gas EoS

   special relativstic with polytrope/ideal fluid

   -> no distinction between co and contravariant
   -> simple expressions and method for root finding*/
  
  double d = P[i].Mass / SphP[i].Volume;
  double sx = SphP[i].Momentum[0] / SphP[i].Volume;
  double sy = SphP[i].Momentum[1] / SphP[i].Volume;
  double sz = SphP[i].Momentum[2] / SphP[i].Volume;
  double t = SphP[i].Energy / SphP[i].Volume;

  double p0 = SphP[i].Pressure;

  double     r,p,e;     // primitives: rho, pressure, vxyz, internal specific energy
  double     kpolyy, gam;         // polytropic constant and index for analytic ideal gas eos

 // internal variables

  int                    iter;               // iterations   
  int                    nmax;            // maximum number of iterations
  double               s2, tpd, root;     // inetrnal variables for con2primSR
  double               dr, de;            // derivatives of pressure wrt r, e
  double               f, df;             // Newton method's function and derivative
  double               rel, prec;         // relative change in pressure and required precisio

  // EoS settings for analytic Eos only - maybe better taken from some module

  kpolyy  =  100.0;

  gam     =  GAMMA;

  // Newton method settings

  nmax = 25;      // probably 10 is sufficient

  prec = 1e-10;   // may need to be reduce for certain problems, e.g. 1e-8 for sieglers second problem; org 1e-10

  p  = p0;       // start value of pressure is its previous value

  s2 =  sx * sx + sy * sy + sz * sz;  

  // do some checks: -- maybe one should monitor how often this is imposed

  /*if ( t < 0.0) { // reset to the previous value by recalculating cons from prim

    double w = 1. / sqrt(1. - (P[i].Vel[0]*P[i].Vel[0] + P[i].Vel[1]*P[i].Vel[1] + P[i].Vel[2]*P[i].Vel[2] ) ); 

    double oldd = w * SphP[i].Density;

    double oldp = SphP[i].Pressure;

    double oldu = SphP[i].Utherm;

    double oldh = 1. + oldu + oldp/SphP[i].Density;

    t =  oldd * ( oldh * w -1. ) - oldp; //thats a bit more

    //printf("tau reset");

    //terminate("tau was reset");

    }*/


  for (iter = 0; iter < nmax; iter++) {

    tpd  =  t + p + d;

    root =  sqrt(tpd*tpd-s2);

    r    =  d / tpd * root ;

    e    =  (root-p*tpd/root-d) / d;

    dr   =  d * s2 / (root*tpd*tpd);

    de   =  p * s2 / (r*root*root*tpd);

    f    =  p - (gam-1.0)*r*e;

    df   =  1.0 - (gam-1.0)*e*dr - (gam-1.0)*r*de;
    
    // newton step
    rel = fabs(f/df/p);

    p   = p - f/df ;

    if ( !(p>1.e-8) ) {
      p = 1.e-8;
	}

  /*if ( !(p > 0.0) ) {
      printf("  ");
      printf("p<0;p=%g,r=%g,e=%g,f=%g,p0=%g,d=%g,s2=%g  ",p,r,e,f,p0,d,s2);
      printf("t=%g,tpd=%g,root=%g,rel=%g",t,tpd,root,rel);
      printf("  ");
      terminate("Scheisse");
      }*/


    if ( rel < prec ) break;  //  ! required accuracy reached
  }

  /*if ( iter == nmax ) {
    printf( "SPECIAL_RELATIVITY: Max iter reached in in con2prim Newton method, rel=%g, prec=%g, f=%g, df=%g, p=%g, nmax=%i, p0=%g\n, i=%i", rel, prec, f, df, p, nmax, p0, i );
    printf( "SPECIAL_RELATIVITY: cons: d=%g, s2=%g, t=%g, p0=%g",d,s2,t,p0 );
    print_particle_info(i);
    terminate( "special relativity update prim fail" );
  }*/

  // compute all variables based on current pressure a last time

  tpd  =  t + p + d;

  root =  sqrt(tpd*tpd-s2);

  SphP[i].Density = d / tpd * root ;

  SphP[i].Utherm = (root-p*tpd/root-d) / d;
  SphP[i].Pressure = p;

  P[i].Vel[0] = sx / tpd;  // also other numerical expression are possible here
  P[i].Vel[1] = sy / tpd;
  P[i].Vel[2] = sz / tpd;
  
  if(p < 0)
    {
      print_particle_info(i);
      terminate( "shit" );
    }
  
  // TODO: add scalars, ...
}

double get_cfl_sound_speed_special_relativity(int p)
{
  double vx = P[p].Vel[0];
  double vy = P[p].Vel[1];
  double vz = P[p].Vel[2];
  
  double v2 = vx*vx + vy*vy + vz*vz;
  
  double Enthalpy = 1. + SphP[p].Utherm + SphP[p].Pressure / SphP[p].Density;
  double csnd = sqrt( GAMMA * SphP[p].Pressure / (SphP[p].Density * Enthalpy) );
  double csnd2 = csnd * csnd;
  
  double lambdamax = 0;
  
  int k;
  for(k = 0; k < 3; k++)
    {
      double lambda0 = fabs( P[p].Vel[k] );
      
      if(lambda0 > lambdamax)
        lambdamax = lambda0;
  
      double t1 = P[p].Vel[k] * (1. - csnd2) / (1. - v2 * csnd2);
      double t2 = csnd / (1. - v2 * csnd2);
      double t3 = (1.-v2) * (1.-v2*csnd2 - P[p].Vel[k]*P[p].Vel[k] * (1.-csnd2) );
  
      double lambdap = fabs( t1 + t2 * sqrt(t3) );
      double lambdam = fabs( t1 - t2 * sqrt(t3) );
      
      if(lambdap > lambdamax)
        lambdamax = lambdap;
      
      if(lambdam > lambdamax)
        lambdamax = lambdam;
    }
  
  return lambdamax;
}

void compute_conserved_quantities_from_ICs_special_relativity(int i)
{
  double v2 = P[i].Vel[0]*P[i].Vel[0] + P[i].Vel[1]*P[i].Vel[1] + P[i].Vel[2]*P[i].Vel[2];
  double w = 1. / sqrt( 1. - v2 );
  
#ifdef MESHRELAX_DENSITY_IN_INPUT
  P[i].Mass = SphP[i].Density * w * SphP[i].Volume;
#else
  SphP[i].Density = P[i].Mass / (w * SphP[i].Volume);
#endif
  
  SphP[i].Pressure = (GAMMA-1.) * SphP[i].Utherm * SphP[i].Density;
  
  double Enthalpy = 1. + SphP[i].Utherm + SphP[i].Pressure / SphP[i].Density;
  
  int k;
  for(k = 0; k < 3; k++)
    SphP[i].Momentum[k] = P[i].Mass * w * Enthalpy * P[i].Vel[k];
  
  SphP[i].Energy = P[i].Mass * w * Enthalpy - SphP[i].Pressure * SphP[i].Volume - P[i].Mass;
}

#endif
