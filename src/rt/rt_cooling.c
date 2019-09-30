/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/rt/rt_cooling.c
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

#if defined(RT_COOLING_PHOTOHEATING)

static double c_light;
static double Lambda;

/* rate1 : photoheating for a blackbody spectrum */
/* rate2 : recombination cooling rate */
/* rate3 : collisional ionization cooling rate */
/* rate4 : collisional excitation cooling rate */
/* rate5 : Bremsstrahlung cooling rate */


double rt_DoHeating(int i, double dtime)
{
  int j;
  double sigma, nH;
  double rate, du, de, c_light;
  double n_gamma, nHI, E;

  set_cosmo_factors_for_current_time();

  c_light = CLIGHT / All.UnitVelocity_in_cm_per_s;
  nH = (HYDROGEN_MASSFRAC * SphP[i].Density * All.cf_a3inv) / (PROTONMASS / All.UnitMass_in_g * All.HubbleParam);
  sigma = 1.49e-18 / All.UnitLength_in_cm / All.UnitLength_in_cm * All.HubbleParam * All.HubbleParam;

  for(j = 0, n_gamma = 0; j < RT_N_DIR; j++)
    n_gamma += SphP[i].DensPhot[j] * All.cf_a3inv;

  nHI = SphP[i].nHI * nH;
  E = 6.4 * ELECTRONVOLT_IN_ERGS / All.UnitEnergy_in_cgs * All.HubbleParam;

  rate = nHI * c_light * E * sigma * n_gamma;

  du = rate * dtime / (SphP[i].Density * All.cf_a3inv);

  return du;
}

double rt_DoCooling(int i, double dtime)
{
  set_cosmo_factors_for_current_time();

  /* now do the cooling */
  double lambda = rt_get_cooling_rate(i, SphP[i].Utherm);
  double du = lambda * dtime / (SphP[i].Density * All.cf_a3inv);

  if(fabs(du) < 0.2 * SphP[i].Utherm)   /* cooling is slow, we can do it explicitly */
    {
      return du;
    }
  else
    {
      /* rapid cooling. Better calculate an implicit solution, which is determined by bisection */

      double u_old = SphP[i].Utherm;
      double u_lower = u_old / sqrt(1.1);
      double u_upper = u_old * sqrt(1.1);
      double ratefact = dtime / (SphP[i].Density * All.cf_a3inv);
      int iter = 0;

      /* bracketing */
      while(u_lower - u_old - ratefact * rt_get_cooling_rate(i, u_lower) > 0)
        {
          u_upper /= 1.1;
          u_lower /= 1.1;

          if(iter++ >= MAXITER)
            terminate("bracketing failure");
        }


      /* bisection */
      double u;
      iter = 0;
      do
        {
          u = 0.5 * (u_lower + u_upper);

          if(u - u_old - ratefact * rt_get_cooling_rate(i, u) > 0)
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
double rt_get_cooling_rate(int i, double utherm)
{
  double temp, molecular_weight;
  double nH;
  double rate2, rate3, rate4, rate5;
  double de2, de3, de4, de5;
#ifdef RT_INCLUDE_HE
  double rateHe2, rateHe3, rateHe4, rateHe5;
  double nHe, nHeIII, deHe2, deHe3, deHe4, deHe5;
#endif

  set_cosmo_factors_for_current_time();

  c_light = CLIGHT / All.UnitVelocity_in_cm_per_s;

  nH = (HYDROGEN_MASSFRAC * SphP[i].Density * All.cf_a3inv) / (PROTONMASS / All.UnitMass_in_g * All.HubbleParam);       //physical
  molecular_weight = 4 / (1 + 3 * HYDROGEN_MASSFRAC + 4 * HYDROGEN_MASSFRAC * SphP[i].n_elec);

  temp = utherm * GAMMA_MINUS1 * (molecular_weight * PROTONMASS / All.UnitMass_in_g * All.HubbleParam) / (BOLTZMANN / All.UnitEnergy_in_cgs * All.HubbleParam);

  /* all rates in erg cm^3 s^-1 in code units */

  /* recombination cooling rate */
  rate2 = 8.7e-27 * pow(temp, 0.5) * pow(temp / 1e3, -0.2) / (1.0 + pow(temp / 1e6, 0.7));
  rate2 *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
  rate2 *= All.HubbleParam * All.HubbleParam;
  rate2 /= All.UnitEnergy_in_cgs / All.HubbleParam;
  de2 = SphP[i].nHII * nH * SphP[i].n_elec * nH * rate2;

  /* collisional ionization cooling rate */
  rate3 = 1.27e-21 * pow(temp, 0.5) * exp(-157809.1 / temp) / (1.0 + pow(temp / 1e5, 0.5));
  rate3 *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
  rate3 *= All.HubbleParam * All.HubbleParam;
  rate3 /= All.UnitEnergy_in_cgs / All.HubbleParam;
  de3 = SphP[i].nHI * nH * SphP[i].n_elec * nH * rate3;

  /* collisional excitation cooling rate */
  rate4 = 7.5e-19 / (1.0 + pow(temp / 1e5, 0.5)) * exp(-118348 / temp);
  rate4 *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
  rate4 *= All.HubbleParam * All.HubbleParam;
  rate4 /= All.UnitEnergy_in_cgs / All.HubbleParam;
  de4 = SphP[i].nHI * nH * SphP[i].n_elec * nH * rate4;

  /* Bremsstrahlung cooling rate */
  rate5 = 1.42e-27 * pow(temp, 0.5);
  rate5 *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
  rate5 *= All.HubbleParam * All.HubbleParam;
  rate5 /= All.UnitEnergy_in_cgs / All.HubbleParam;
  de5 = SphP[i].nHII * nH * SphP[i].n_elec * nH * rate5;

  Lambda = de2 + de3 + de4 + de5;

#ifdef RT_INCLUDE_HE
  nHe = ((1.0 - HYDROGEN_MASSFRAC) * SphP[i].Density * All.cf_a3inv) / (4.0 * PROTONMASS / All.UnitMass_in_g * All.HubbleParam);
  nHeIII = 1 - (SphP[i].nHeI + SphP[i].nHeII);

  /* recombination cooling rate */
  rateHe2 = 1.55e-26 * pow(temp, 0.3647);
  rateHe2 *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
  rateHe2 *= All.HubbleParam * All.HubbleParam;
  rateHe2 /= All.UnitEnergy_in_cgs / All.HubbleParam;
  deHe2 = SphP[i].nHeII * nHe * SphP[i].n_elec * nH * rateHe2;

  rateHe2 = 3.48e-26 * pow(temp, 0.5) * pow(temp / 1e3, -0.2) / (1.0 + pow(temp / 1e6, 0.7));
  rateHe2 *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
  rateHe2 *= All.HubbleParam * All.HubbleParam;
  rateHe2 /= All.UnitEnergy_in_cgs / All.HubbleParam;
  deHe2 += nHeIII * nHe * SphP[i].n_elec * nH * rateHe2;

  /* collisional ionization cooling rate */
  rateHe3 = 9.38e-22 * pow(temp, 0.5) * exp(-285335.4 / temp) / (1.0 + pow(temp / 1e5, 0.5));
  rateHe3 *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
  rateHe3 *= All.HubbleParam * All.HubbleParam;
  rateHe3 /= All.UnitEnergy_in_cgs / All.HubbleParam;
  deHe3 = SphP[i].nHeI * nHe * SphP[i].n_elec * nH * rateHe3;

  rateHe3 = 4.95e-22 * pow(temp, 0.5) * exp(-631515 / temp) / (1.0 + pow(temp / 1e5, 0.5));
  rateHe3 *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
  rateHe3 *= All.HubbleParam * All.HubbleParam;
  rateHe3 /= All.UnitEnergy_in_cgs / All.HubbleParam;
  deHe3 += SphP[i].nHeII * nHe * SphP[i].n_elec * nH * rateHe3;

  /* collisional excitation cooling rate */
  rateHe4 = 5.54e-17 * pow(temp, -0.397) / (1.0 + pow(temp / 1e5, 0.5)) * exp(-473638 / temp);
  rateHe4 *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
  rateHe4 *= All.HubbleParam * All.HubbleParam;
  rateHe4 /= All.UnitEnergy_in_cgs / All.HubbleParam;
  deHe4 = SphP[i].nHeII * nHe * SphP[i].n_elec * nH * rateHe4;

  rateHe4 = 9.10e-27 * pow(temp, -0.1687) / (1.0 + pow(temp / 1e5, 0.5)) * exp(-13179 / temp);
  rateHe4 *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
  rateHe4 *= 1.0 / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
  rateHe4 *= All.HubbleParam * All.HubbleParam * All.HubbleParam * All.HubbleParam * All.HubbleParam;
  rateHe4 /= All.UnitEnergy_in_cgs / All.HubbleParam;
  deHe4 += SphP[i].nHeII * nHe * SphP[i].n_elec * nH * SphP[i].n_elec * nH * rateHe4;

  /* Bremsstrahlung cooling rate */
  rateHe5 = 1.42e-27 * pow(temp, 0.5);
  rateHe5 *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
  rateHe5 *= All.HubbleParam * All.HubbleParam;
  rateHe5 /= All.UnitEnergy_in_cgs / All.HubbleParam;
  deHe5 = (SphP[i].nHeII * nHe * SphP[i].n_elec * nH + 4.0 * nHeIII * nHe * SphP[i].n_elec * nH) * rateHe5;

  Lambda += deHe2 + deHe3 + deHe4 + deHe5;
#endif


  return -Lambda;
}



#endif
