/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/MRT/RT_cooling.c
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../allvars.h"
#include "../proto.h"
#include "./RT.h"

#include "../cooling/cooling_proto.h"

#if defined(MRT_COOLING_HEATING)

/* rate1 : photoheating for a blackbody spectrum */
/* rate2 : recombination cooling rate */
/* rate3 : collisional ionization cooling rate */
/* rate4 : collisional excitation cooling rate */
/* rate5 : Bremsstrahlung cooling rate */


static double XH = HYDROGEN_MASSFRAC;   /* hydrogen abundance by mass */
static double yhelium;


static double ne, necgs, nHcgs;
static double nH0, nHp, nHep, nHe0, nHepp;

static double DoCool_u_old_input, DoCool_rho_input, DoCool_dt_input, DoCool_ne_guess_input;


/* now do the heating (note: we know how many photons we absorbed) */
double mrt_DoHeating(int i, double dt_internal)
{

  //  mpi_printf("Entered Here\n") ;
  double sigma, a3inv, nH;
  double hubble_a, rate, du, c_light;
  double n_gamma, nHI;
  double E;
  double Lambda_Heat, utherm;

  if(All.ComovingIntegrationOn)
    {
      a3inv = 1.0 / All.Time / All.Time / All.Time;
      hubble_a = hubble_function(All.Time);
    }
  else
    {
      a3inv = hubble_a = 1.0;
    }

  Lambda_Heat = mrt_get_heating_rate(i);

  du = Lambda_Heat * dt_internal / hubble_a / (SphP[i].Density * a3inv);        /*LVS: rt_DoHeating.c in arepo has no 1/hubble_a factor here?? */

  //  if(P[i].Pos[0]<All.BoxSize/30.0)
  //printf("Heating=%g\t%g\t%g\t%g\n", du, Lambda_Heat, dt_internal, SphP[i].Density) ;
  return du;
}


double mrt_DoCooling(int i, double dt_internal)
{
  double dtime, a3inv, iter;
  double lambda, du, de;
  double u_old, u_lower, u_upper, ratefact, u;
  double utherm;

  if(All.ComovingIntegrationOn)
    {
      dtime = dt_internal / hubble_function(All.Time);
      a3inv = 1.0 / All.Time / All.Time / All.Time;
    }
  else
    {
      dtime = dt_internal;
      a3inv = 1.0;
    }

  utherm = SphP[i].Utherm;

  /* do the cooling */
  lambda = mrt_get_cooling_rate(i, utherm);
  du = lambda * dtime / (SphP[i].Density * a3inv);

  if(fabs(du) < 0.2 * utherm)
    {
      /* cooling is slow, we can do it explicitly */
      return du;
    }
  else
    {
      /* rapid cooling. Better calculate an implicit solution, which is determined by bisection */
      u_old = utherm;
      u_lower = u_old / sqrt(1.1);
      u_upper = u_old * sqrt(1.1);
      ratefact = dtime / (SphP[i].Density * a3inv);
      iter = 0;

      /* bracketing */
      while(u_lower - u_old - ratefact * mrt_get_cooling_rate(i, u_lower) > 0)
        {
          u_upper /= 1.1;
          u_lower /= 1.1;


          if(iter++ >= 1000)    //MAXITER)
            terminate("bracketing failure");
        }

      /* bisection */
      iter = 0;
      do
        {
          u = 0.5 * (u_lower + u_upper);

          if(u - u_old - ratefact * mrt_get_cooling_rate(i, u) > 0)
            u_upper = u;
          else
            u_lower = u;

          du = u_upper - u_lower;

          iter++;

          if(iter >= (MAXITER - 10))
            printf("u= %g\n", u);

          if(iter >= MAXITER)
            terminate("convergence failure");
        }
      while(fabs(du / u) > 1.0e-6);

      du = u - u_old;

      return du;
    }

}

/* returns cooling rate */
double mrt_get_cooling_rate(int i, double utherm)
{
  double Lambda = 0.0 ;

  double temp, molecular_weight;
  double a3inv;
  double nH;
  double rate2, rate3, rate4, rate5;
  double de2, de3, de4, de5;
#ifdef MRT_INCLUDE_HE
  double rateHe2, rateHe3, rateHe4, rateHe5;
  double deHe2, deHe3, deHe4, deHe5;
#endif

  double fac = All.UnitTime_in_s / pow(All.UnitLength_in_cm, 3) / All.UnitEnergy_in_cgs * pow(All.HubbleParam, 3);

  if(All.ComovingIntegrationOn)
    a3inv = 1 / All.Time / All.Time / All.Time;
  else
    a3inv = 1;

  nH = HYDROGEN_MASSFRAC * SphP[i].Density * a3inv / PROTONMASS * All.UnitMass_in_g / All.HubbleParam;  //physical
  molecular_weight = 4 / (1 + 3 * HYDROGEN_MASSFRAC + 4 * HYDROGEN_MASSFRAC * SphP[i].Ne);

  temp = utherm * GAMMA_MINUS1 * (molecular_weight * PROTONMASS / All.UnitMass_in_g * All.HubbleParam) / (BOLTZMANN / All.UnitEnergy_in_cgs * All.HubbleParam);

#ifndef MRT_NO_UV
  //printf("temperature = %g", temp) ;

  /* all rates in erg cm^3 s^-1 in code units */
  /* recombination cooling rate */
  rate2 = 8.7e-27 * sqrt(temp) * pow(temp / 1e3, -0.2) / (1.0 + pow(temp / 1e6, 0.7)) * fac;
  de2 = SphP[i].HII * nH * SphP[i].Ne * nH * rate2;

  /* collisional ionization cooling rate */
  rate3 = 1.27e-21 * sqrt(temp) * exp(-157809.1 / temp) / (1.0 + sqrt(temp / 1e5)) * fac;
  de3 = SphP[i].HI * nH * SphP[i].Ne * nH * rate3;

  /* collisional excitation cooling rate */
  rate4 = 7.5e-19 / (1.0 + sqrt(temp / 1e5)) * exp(-118348 / temp) * fac;
  de4 = SphP[i].HI * nH * SphP[i].Ne * nH * rate4;

  /* Bremsstrahlung cooling rate */
  rate5 = 1.42e-27 * sqrt(temp) * fac;
  de5 = SphP[i].HII * nH * SphP[i].Ne * nH * rate5;

  Lambda = de2 + de3 + de4 + de5;

  /* inverse Compton cooling rate */
  if(All.ComovingIntegrationOn)
    Lambda += 5.406e-36 * SphP[i].Ne * (temp - (2.73 / All.Time)) / pow(All.Time, 4) * fac;

#ifdef MRT_INCLUDE_HE
  /* recombination cooling rate */
  rateHe2 = 1.55e-26 * pow(temp, 0.3647) * fac;
  deHe2 = SphP[i].HeII * nH * SphP[i].Ne * nH * rateHe2;

  rateHe2 = 3.48e-26 * sqrt(temp) * pow(temp / 1e3, -0.2) / (1.0 + pow(temp / 1e6, 0.7)) * fac;
  deHe2 += SphP[i].HeIII * nH * SphP[i].Ne * nH * rateHe2;

  /* collisional ionization cooling rate */
  rateHe3 = 9.38e-22 * sqrt(temp) * exp(-285335.4 / temp) / (1.0 + sqrt(temp / 1e5)) * fac;
  deHe3 = SphP[i].HeI * nH * SphP[i].Ne * nH * rateHe3;

  rateHe3 = 4.95e-22 * sqrt(temp) * exp(-631515 / temp) / (1.0 + sqrt(temp / 1e5)) * fac;
  deHe3 += SphP[i].HeII * nH * SphP[i].Ne * nH * rateHe3;

  rateHe3 = 5.01e-27 * pow(temp, -0.1687) / (1.0 + sqrt(temp / 1e5)) * exp(-55338 / temp) * fac;
  rateHe3 *= pow(All.HubbleParam / All.UnitLength_in_cm, 3);
  deHe3 += SphP[i].HeII * nH * SphP[i].Ne * nH * SphP[i].Ne * nH * rateHe3;

  /* collisional excitation cooling rate */
  rateHe4 = 5.54e-17 * pow(temp, -0.397) / (1.0 + sqrt(temp / 1e5)) * exp(-473638 / temp) * fac;
  deHe4 = SphP[i].HeII * nH * SphP[i].Ne * nH * rateHe4;

  rateHe4 = 9.10e-27 * pow(temp, -0.1687) / (1.0 + sqrt(temp / 1e5)) * exp(-13179 / temp) * fac;
  rateHe4 *= pow(All.HubbleParam / All.UnitLength_in_cm, 3);
  deHe4 += SphP[i].HeII * nH * SphP[i].Ne * nH * SphP[i].Ne * nH * rateHe4;

  /* Bremsstrahlung cooling rate */
  rateHe5 = 1.42e-27 * sqrt(temp) * fac;
  deHe5 = (SphP[i].HeII + 4.0 * SphP[i].HeIII + SphP[i].Ne) * nH * rateHe5;

  Lambda += deHe2 + deHe3 + deHe4 + deHe5;
#endif
#endif

  /*#ifdef MRT_IR_LTE
  double cspeed = 2.99792458e10 / All.UnitVelocity_in_cm_per_s ;
  for(int num1=0; num1<IR_BINS;num1++)
    Lambda += SphP[i].KappaIR[num1]*cspeed*radiation_constant*pow(temp,4) ;


  //  if(P[i].Pos[0]<All.BoxSize/30.0)
  // printf("cool=%g\t%g\t%g\t%g\t%g\n", SphP[i].KappaIR[0], cspeed, radiation_constant, pow(temp,4), SphP[i].KappaIR[0]*cspeed*radiation_constant*pow(temp,4) ) ;

  #endif*/

  return -Lambda;
}

/*#ifndef MRT_MULTI_FREQUENCY
double mrt_get_heating_rate(int i)
{
  double sigma, a3inv, nH;
  double hubble_a, rate, c_light;
  double n_gamma, nHI;
  double E;

  if(All.ComovingIntegrationOn)
    {
      a3inv = 1.0 / All.Time / All.Time / All.Time;
      hubble_a = hubble_function(All.Time);
    }
  else
    {
      a3inv = hubble_a = 1.0;
    }

  c_light = CLIGHT / All.UnitVelocity_in_cm_per_s;

  nH = HYDROGEN_MASSFRAC * SphP[i].Density * a3inv / PROTONMASS * All.UnitMass_in_g / All.HubbleParam;
  nHI = SphP[i].HI * nH;

  sigma = 1.63e-18 / All.UnitLength_in_cm / All.UnitLength_in_cm * All.HubbleParam * All.HubbleParam;
  //n_gamma = SphP[i].n_gamma[0] / P[i].Mass * a3inv;
  n_gamma = SphP[i].DensPhot * 1e63 ; /// SphP[i].Volume ;
  //E = 30.0 * ELECTRONVOLT_IN_ERGS / All.UnitEnergy_in_cgs * All.HubbleParam;
  E = 0.0 ;
  rate = nHI * c_light * E * sigma * n_gamma;

  return rate;
}

#else
*/
double mrt_get_heating_rate(int i)
{
  double rate = 0 ;
#ifndef MRT_NO_UV
  int j;
  double a3inv, nH;
  double hubble_a ;
  double nHI, n_gamma;
  double c_light;

#ifdef MRT_INCLUDE_HE
  double nHeI, nHeII, nHeIII;
#endif

  c_light = CLIGHT / All.UnitVelocity_in_cm_per_s;

  if(All.ComovingIntegrationOn)
    {
      a3inv = All.Time / All.Time / All.Time;
      hubble_a = hubble_function(All.Time);
    }
  else
    {
      a3inv = hubble_a = 1.0;
    }

  nH = HYDROGEN_MASSFRAC * SphP[i].Density * a3inv / PROTONMASS * All.UnitMass_in_g / All.HubbleParam;
  nHI = SphP[i].HI * nH;

#ifdef MRT_INCLUDE_HE
  nHeI = SphP[i].HeI * nH;
  nHeII = SphP[i].HeII * nH;
  nHeIII = SphP[i].HeIII * nH;
#endif

  for(j = 0; j < UV_BINS; j++)
    {
      //      n_gamma = SphP[i].n_gamma[j] / P[i].Mass * a3inv;
      n_gamma = SphP[i].DensPhot[j] * 1e63  ;

      //      if(nu[j] >= 13.6)
      rate += nHI * c_light * mrt_sigma_HI[j] * G_HI[j] * n_gamma ;
      
#ifdef MRT_INCLUDE_HE
	// if(nu[j] >= 24.6)
      rate += nHeI * c_light * mrt_sigma_HeI[j] * G_HeI[j] * n_gamma;
      
	//      if(nu[j] >= 54.4)
      rate += nHeII * c_light * mrt_sigma_HeII[j] * G_HeII[j] * n_gamma;
#endif
    }
#endif

  /*#ifdef MRT_IR_LTE
  for(int num1=UV_BINS;num1<(UV_BINS+IR_BINS);num1++)
      rate += SphP[i].KappaIR[num1-UV_BINS]*c_internal_units*SphP[i].DensPhot[num1] ;

  //  if(P[i].Pos[0]<All.BoxSize/30.0)
  // printf("heat=%g\t%g\t%g\t%g\n", SphP[i].KappaIR[0], c_internal_units, SphP[i].DensPhot[0], SphP[i].KappaIR[0]*c_internal_units*SphP[i].DensPhot[0] ) ;
  #endif*/
  //  return 0.0 ;
  return rate;
}


//#endif

/* -------------------------*/


/* returns cooling time.
 * NOTE: If we actually have heating, a cooling time of 0 is returned.
 */
double mrt_GetCoolingTime(int i, double u_old, double rho, double *ne_guess)
{
  double u;
  double ratefact;
  double LambdaNet, coolingtime;

  DoCool_u_old_input = u_old;
  DoCool_rho_input = rho;
  DoCool_ne_guess_input = *ne_guess;

  rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;    /* convert to physical cgs units */
  u_old *= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;


  nHcgs = XH * rho / PROTONMASS;        /* hydrogen number dens in cgs units */
  ratefact = nHcgs * nHcgs / rho;

  u = u_old;

  LambdaNet = mrt_get_heating_rate(i) + mrt_get_cooling_rate(i, u);

  /* bracketing */

  if(LambdaNet >= 0)            /* ups, we have actually heating due to UV background */
    return 0;

  coolingtime = u_old / (-ratefact * LambdaNet);

  coolingtime *= All.HubbleParam / All.UnitTime_in_s;

  return coolingtime;
}

#endif
