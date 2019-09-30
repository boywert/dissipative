/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/MRT/RT_chem_ps2011.c
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

/* Solves the chemical network (H and He) computing the number of photon absorptions
 * and expected recombinations. Uses the method described in Petkova & Springel 2011 rather
 * than the original 2009.
 * You must use this function if using the RADPRESS_OPT_THICK option.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

#include "../allvars.h"
#include "../proto.h"
#include "../voronoi.h"


#if defined(MRT) && defined(MRT_CHEMISTRY_PS2011)
#ifndef MRT_MULTI_FREQUENCY

/* 
 * this function calculates the ODEs that describe the evolution of
 * neutral hydrogen fraction and photon abundance. The input vector is
 * y[ nHI, nGamma ], where nHI and nGamma are the normalized (relative to
 * the total hydrogen abundance) abundances. The ouput vector contains
 * the rates, f = dy/dt. The coefficient vector, coeff[], describes the
 * prefactors for recombination, collisional ionization, and photo ionization.
 */
int mrt_rate_ODEs(double t, const double y[], double f[], void *params)
{
  double *coeff = (double *) params;
  f[0] = coeff[0] * (1 - y[0]) * (1 - y[0]) - coeff[2] * (1 - y[0]) * y[0] - coeff[1] * y[0] * y[1];
  f[1] = -coeff[1] * y[0] * y[1];
  return GSL_SUCCESS;
}


/* 
 * This function advances the chemical abundances, and does the photo heating
 * and cooling, as necessary.
 */
void mrt_update_chemistry_ps2011(void)
{
  int idx, i, j;
  double nH, temp, alpha, gamma, molecular_weight, nHI, nHII;
  double sigma, c_light;
  double dt;
  double count_part_ode = 0, count_part_ode_iter = 0;
  double tot_count_part_ode, tot_count_part_ode_iter;
  double fac, nGamma, nGamma_new, nHI_new, nHII_new;
  double mean_nHI, mean_nHI_all;


#ifdef MRT_INCLUDE_HE
  double alpha_HeII, alpha_HeIII, gamma_HeI, gamma_HeII;
  double nHe, nHeII, nHeIII;
  double y, D, E, F, J;
  double mean_nHeI, mean_nHeI_all;
  mean_nHeI = 0;
#endif

  set_cosmo_factors_for_current_time();

  sigma = 1.49e-18 / All.UnitLength_in_cm / All.UnitLength_in_cm * All.HubbleParam * All.HubbleParam; 
  c_light = CLIGHT / All.UnitVelocity_in_cm_per_s; 

     mean_nHI = 0.; 
   for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++) 
     { 
       i = TimeBinsHydro.ActiveParticleList[idx]; 
       if(i < 0) 
         continue; 

       dt = (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval ; 
      

       if(All.ComovingIntegrationOn) 
         dt = dt / All.cf_hubble_a; 

       if(dt <= 0.) 
         continue; 

#ifdef RADPRESS_OPT_THICK 
       double n_gamma_old = SphP[i].n_gamma[0]; 
#endif 

       nH = (HYDROGEN_MASSFRAC * SphP[i].Density * All.cf_a3inv) / (PROTONMASS / All.UnitMass_in_g * All.HubbleParam); 

       molecular_weight = 4 / (1 + 3 * HYDROGEN_MASSFRAC + 4 * HYDROGEN_MASSFRAC * SphP[i].Ne); 

       //       molecular_weight = 1.0 ;

       temp = SphP[i].Utherm * GAMMA_MINUS1 * (molecular_weight * PROTONMASS / All.UnitMass_in_g * All.HubbleParam) / (BOLTZMANN / All.UnitEnergy_in_cgs * All.HubbleParam); 

       //\* collisional ionization rate *\/ 
       gamma = 5.85e-11 * pow(temp, 0.5) * exp(-157809.1 / temp) / (1.0 + pow(temp / 1e5, 0.5)); 
       gamma *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm; 
       gamma *= All.HubbleParam * All.HubbleParam; 
#ifdef MRT_NOCOLLISION_IONIZATION 
       gamma = 0.; 
#endif 

       ///\* alpha_B recombination coefficient *\/ */
       alpha = 2.59e-13 * pow(temp / 1e4, -0.7); 
       alpha *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm; 
       alpha *= All.HubbleParam * All.HubbleParam; 


#ifdef MRT_INCLUDE_HE 
       nHe = ((1.0 - HYDROGEN_MASSFRAC) * SphP[i].Density * All.cf_a3inv) / (4.0 * PROTONMASS / All.UnitMass_in_g * All.HubbleParam); 

	 /*       /\* collisional ionization rate *\/ */
       gamma_HeI = 2.38e-11 * pow(temp, 0.5) * exp(-285335.4 / temp) / (1.0 + pow(temp / 1e5, 0.5)); 
       gamma_HeI *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm; 
       gamma_HeI *= All.HubbleParam * All.HubbleParam; 

       gamma_HeII = 5.68e-12 * pow(temp, 0.5) * exp(-631515 / temp) / (1.0 + pow(temp / 1e5, 0.5)); 
       gamma_HeII *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm; 
       gamma_HeII *= All.HubbleParam * All.HubbleParam; 

/*       /\* alpha_B recombination coefficient *\/ */
       alpha_HeII = 1.5e-10 * pow(temp, -0.6353); 
       alpha_HeII *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm; 
       alpha_HeII *= All.HubbleParam * All.HubbleParam; 

       alpha_HeIII = 3.36e-10 * pow(temp, -0.5) * pow(temp / 1e3, -0.2) / (1.0 + pow(temp / 1e6, 0.7)); 
       alpha_HeIII *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm; 
       alpha_HeIII *= All.HubbleParam * All.HubbleParam; 
#ifdef MRT_NOCOLLISION_IONIZATION 
       gamma_HeI = 0.; 
       gamma_HeII = 0.; 
 #endif 
 #endif

      nHI = SphP[i].HI;
      nHII = SphP[i].HII;

      if(nHI + nHII < 0.95)
	terminate("nHI = %g, nHII = %g, nHI+nHII = %g\n", nHI, nHII) ;

      nGamma = SphP[i].Cons_DensPhot[0] * 1e63 / nH / SphP[i].Volume;
      sigma = mrt_sigma_HI[0];

      int calculated_flag = 0;

      double d_nGamma = -c_light * sigma * nH * dt * nHI * nGamma;      /* linear estimate of the photo abundance change */

      if(fabs(d_nGamma) <= 0.05 * nGamma)       /* we are nearly in photoionization equilibrium, and/or the photon density is so large that it should change little */
        {
          nGamma_new = nGamma + d_nGamma;

          if(isnan(nGamma_new))
	    {
	      printf("%g\t%g\t%g\n", nGamma, d_nGamma, nGamma_new) ;
	      terminate("LVS0: nGamma_new is nan\n");
	    }
          /* let's solve implicitely for new hydrogen abundance */

          double A = alpha * nH * dt;
          double B = c_light * sigma * nH * nGamma * dt;
          double D = gamma * nH * dt;

          //if(A != 0)
          if(A >= 0)
            {
              double b = 2 * A + B + D + 1;
              double det = b * b - 4 * (A + D) * (A + nHI);
              if(det >= 0)
                {
                  nHI_new = (b - sqrt(det)) / (2 * (A + D));

                  if(fabs(nHI_new - nHI) < 0.25 * nHI)  /* only accept if change is sufficiently small */
                    calculated_flag = 1;
                }
            }
          else
            {
              if(CLIGHT != 0)
                terminate("not implemented");

            }
        }
      if(calculated_flag == 0 && nGamma < 1.0e-2 * nHI) /* here photoionization cannot change the neutral fraction much */
        {
          /* first, calculate an implicit solution for the new gamma, neglecting a change in nHI */
          fac = c_light * sigma * nH * dt;      /* photo ionization factor */
          if(fac > 0)
            {
	      //nGamma_new = (-(1 + fac * (1 - nHII - nGamma)) + sqrt(pow(1 + fac * (1 - nHII - nGamma), 2) + 4 * fac * nGamma)) / (2 * fac); 
              /*LVS: changed this for eq. 28 in Petkova & Springel 2011b */
              nGamma_new = nGamma / (1 + c_light * sigma * nH * (1 - nHII) * dt);

              if(isnan(nGamma_new))
                terminate("LVS3: nGamma_new is nan\n");
            }
          else
            {
              nGamma_new = nGamma;
              if(isnan(nGamma_new))
                terminate("LVS4: nGamma_new is nan\n");
            }
          /* now calculate the new neutral and ionized hydrogen fraction */
          double A = alpha * nH * dt;   /* recombination factor */
          double B = gamma * nH * dt;   /* collisional ionization factor */

          if((A + B) != 0)
            {
              double det = (1 - B) * (1 - B) + 4 * (A + B) * (nHII + nGamma - nGamma_new);      /* Eq. 27 */
              if(det > 0)
                {
                  nHII_new = (-1 + B + sqrt(det)) / (2 * (A + B));
                  nHI_new = 1 - nHII_new;

                  if(fabs(nHI_new - nHI) < 0.25 * nHI)  /* only accept if change is sufficiently small */
                    calculated_flag = 1;
                }
            }
        }

      if(calculated_flag == 0)  /* if no accurate solution has been found yet, integrate the problem with an ODE solver */
        {
          count_part_ode++;

          /* let's here explicitely integrate the (stiff) differential equation */
          const gsl_odeiv_step_type *T = gsl_odeiv_step_rkf45;
          gsl_odeiv_step *s = gsl_odeiv_step_alloc(T, 2);
          gsl_odeiv_control *c = gsl_odeiv_control_y_new(1e-5, 0.0);
          gsl_odeiv_evolve *e = gsl_odeiv_evolve_alloc(2);

          double coeff[3] = { alpha * nH, c_light * sigma * nH, gamma * nH };
          gsl_odeiv_system sys = { mrt_rate_ODEs, NULL, 2, coeff };

          double t = 0.0, t1 = dt;
          //double h = 0.25 * dt;
          //double h = 0.025 * dt;
          double h = 0.0025 * dt;       /* LVS: this can be adjusted according to resolution needed. Smaller is better but slower */
          double y[2] = { nHI, nGamma };

          while(t < t1)
            {
              int status = gsl_odeiv_evolve_apply(e, c, s, &sys, &t, t1, &h, y);

              if(status != GSL_SUCCESS)
                terminate("we seem to have failed in the integration\n");

              count_part_ode_iter++;
            }

          gsl_odeiv_evolve_free(e);
          gsl_odeiv_control_free(c);
          gsl_odeiv_step_free(s);

          nHI_new = y[0];
          nGamma_new = y[1];
          if(isnan(nGamma_new))
            terminate("LVS5: nGamma_new is nan | HI = %g, HII = %g, nH = %g\n", SphP[i].HI, SphP[i].HII, nH);
        }


      /* make sure that we stay within physical value bounds */
      if(nGamma_new < 0)
        nGamma_new = 0;

      if(isnan(nGamma_new))
        terminate("nGamma_new is nan\n");

      if(isnan(nHI_new))
        terminate("nHI_new is nan\n");

      if(nHI_new < 0)
        nHI_new = 0;

      if(nHI_new > 1)
        nHI_new = 1;

      nHII_new = 1 - nHI_new;

      /* store the new abundance results */

      SphP[i].HI = nHI_new;
      SphP[i].HII = nHII_new;
      SphP[i].Ne = SphP[i].HII;


      //      double nH_times_volume = HYDROGEN_MASSFRAC * P[i].Mass / PROTONMASS * All.UnitMass_in_g / All.HubbleParam;
      double nH_times_volume = P[i].Mass ;
      SphP[i].nHI = SphP[i].HI * nH_times_volume ;
      SphP[i].nHII = SphP[i].HII * nH_times_volume ;
      SphP[i].ne = SphP[i].Ne * nH_times_volume ;



      double ratio = (nGamma_new / 1e63 * nH * SphP[i].Volume ) / SphP[i].Cons_DensPhot[0] ;
            SphP[i].Cons_DensPhot_absorbed[0] = SphP[i].Cons_DensPhot[0] * (1.0 - ratio) ;
      SphP[i].Cons_DensPhot[0] *= ratio ;

      int num1 ;
      for(num1=0;num1<3;num1++)
	  SphP[i].Cons_RT_F[0][num1] *= ratio ;	  


#ifdef RADPRESS_OPT_THICK
      SphP[i].n_gamma_abs[0] = n_gamma_old - SphP[i].n_gamma[0];        /* number density of photons absorbed in this timestep */
      if(SphP[i].n_gamma_abs[0] < 0)
        terminate("number of absorbed photons is negative=%g n_gamma_old=%g n_gamma_new=%g \n", SphP[i].n_gamma_abs[0], n_gamma_old, SphP[i].n_gamma[0]);
#endif

      mean_nHI += SphP[i].HI * SphP[i].Volume;


#ifdef MRT_INCLUDE_HE
      nHeIII = 1 - (SphP[i].HeI + SphP[i].HeII);

      SphP[i].Ne += SphP[i].HeII * (1.0 - HYDROGEN_MASSFRAC) / (4.0 * HYDROGEN_MASSFRAC);
      SphP[i].Ne += 2.0 * nHeIII * (1.0 - HYDROGEN_MASSFRAC) / (4.0 * HYDROGEN_MASSFRAC);

      D = dt * gamma_HeII * nH;
      E = dt * alpha_HeIII * nH;
      F = dt * gamma_HeI * nH;
      J = dt * alpha_HeII * nH;

      nHeII = SphP[i].HeII + F * SphP[i].Ne - ((F * SphP[i].Ne - E * SphP[i].Ne) / (1.0 + E * SphP[i].Ne) * nHeIII);

      nHeII /= 1.0 + F * SphP[i].Ne + D * SphP[i].Ne + J * SphP[i].Ne + ((F * SphP[i].Ne - E * SphP[i].Ne) / (1.0 + E * SphP[i].Ne) * D * SphP[i].Ne);

      if(nHeII < 0)
        nHeII = 0.0;

      if(nHeII > 1)
        nHeII = 1.0;

      nHeIII = (nHeIII + D * nHeII * SphP[i].Ne) / (1.0 + E * SphP[i].Ne);

      if(nHeIII < 0)
        nHeIII = 0.0;

      if(nHeIII > 1)
        nHeIII = 1.0;


      SphP[i].Ne = SphP[i].HII;
      SphP[i].Ne += nHeII * (1.0 - HYDROGEN_MASSFRAC) / (4.0 * HYDROGEN_MASSFRAC);
      SphP[i].Ne += 2.0 * nHeIII * (1.0 - HYDROGEN_MASSFRAC) / (4.0 * HYDROGEN_MASSFRAC);

      SphP[i].HeII = nHeII;

      SphP[i].HeI = 1.0 - (SphP[i].HeII + nHeIII);

      if(SphP[i].HeI < 0 || SphP[i].HeI > 1 || isnan(SphP[i].HeI))
        {
          terminate("ERROR, wrong value of nHeI %g\n", SphP[i].HeI);
        }

      mean_nHeI += SphP[i].HeI * SphP[i].Volume;
#endif

    }

  MPI_Reduce(&count_part_ode, &tot_count_part_ode, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&count_part_ode_iter, &tot_count_part_ode_iter, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  MPI_Reduce(&mean_nHI, &mean_nHI_all, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#ifdef MRT_INCLUDE_HE
  MPI_Reduce(&mean_nHeI, &mean_nHeI_all, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif

  if(ThisTask == 0)
    {
      fprintf(FdMRT, "%g %g ", All.Time, mean_nHI_all / All.BoxSize / All.BoxSize / All.BoxSize);
#ifndef MRT_INCLUDE_HE
      fprintf(FdMRT, "\n");
#else
      fprintf(FdMRT, "%g\n", mean_nHeI_all / All.BoxSize / All.BoxSize / All.BoxSize);
#endif
      myflush(FdMRT);
    }

  if(ThisTask == 0)
    {
      printf
        ("Ionization balance solver: fraction of particle requiring ODE: %g, average number of iterations=%g\n",
         tot_count_part_ode / All.TotNumGas, tot_count_part_ode_iter / (tot_count_part_ode + 1.0e-8));
    }

}

#endif /* MULTI-FREQUENCY */

#endif
