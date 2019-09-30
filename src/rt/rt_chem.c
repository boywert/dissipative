/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/rt/rt_chem.c
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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

#include "../allvars.h"
#include "../proto.h"
#include "../voronoi.h"


#ifdef RT_ADVECT

/* 
 * this function calculates the ODEs that describe the evolution of
 * neutral hydrogen fraction and photon abundance. The input vector is
 * y[ nHI, nGamma ], where nHI and nGamma are the normalized (relative to
 * the total hydrogen abundance) abundances. The ouput vector contains
 * the rates, f = dy/dt. The coefficient vector, coeff[], describes the
 * prefactors for recombination, collisional ionization, and photo ionization.
 */
int rt_rate_ODEs(double t, const double y[], double f[], void *params)
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
void rt_update_chemistry(double dt_internal)
{
  int i, j;
  double nH, temp, alpha, gamma, molecular_weight, nHI, nHII;
  double sigma, c_light;
  double dt;
  double count_part_ode = 0, count_part_ode_iter = 0;
  double tot_count_part_ode, tot_count_part_ode_iter;
  double fac, nGamma, nGamma_new, nHI_new, nHII_new;
  double mean_nHI, mean_nHI_all;


#ifdef RT_INCLUDE_HE
  double alpha_HeII, alpha_HeIII, gamma_HeI, gamma_HeII;
  double nHe, nHeII, nHeIII;
  double y, D, E, F, J;
  double mean_nHeI, mean_nHeI_all;
  mean_nHeI = 0;
#endif

  set_cosmo_factors_for_current_time();

  sigma = 1.49e-18 / All.UnitLength_in_cm / All.UnitLength_in_cm * All.HubbleParam * All.HubbleParam;
  c_light = CLIGHT / All.UnitVelocity_in_cm_per_s;

  for(i = 0, mean_nHI = 0; i < NumGas; i++)
    if(P[i].Type == 0)
      {
        dt = dt_internal / All.cf_hubble_a;

        nH = (HYDROGEN_MASSFRAC * SphP[i].Density * All.cf_a3inv) / (PROTONMASS / All.UnitMass_in_g * All.HubbleParam);

        molecular_weight = 4 / (1 + 3 * HYDROGEN_MASSFRAC + 4 * HYDROGEN_MASSFRAC * SphP[i].n_elec);

        temp = SphP[i].Utherm * GAMMA_MINUS1 * (molecular_weight * PROTONMASS / All.UnitMass_in_g * All.HubbleParam) / (BOLTZMANN / All.UnitEnergy_in_cgs * All.HubbleParam);

        /* collisional ionization rate */
        gamma = 5.85e-11 * pow(temp, 0.5) * exp(-157809.1 / temp) / (1.0 + pow(temp / 1e5, 0.5));
        gamma *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
        gamma *= All.HubbleParam * All.HubbleParam;

        /* alpha_B recombination coefficient */
        alpha = 2.59e-13 * pow(temp / 1e4, -0.7);
        alpha *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
        alpha *= All.HubbleParam * All.HubbleParam;

#ifdef RT_INCLUDE_HE
        nHe = ((1.0 - HYDROGEN_MASSFRAC) * SphP[i].Density * All.cf_a3inv) / (4.0 * PROTONMASS / All.UnitMass_in_g * All.HubbleParam);

        /* collisional ionization rate */
        gamma_HeI = 2.38e-11 * pow(temp, 0.5) * exp(-285335.4 / temp) / (1.0 + pow(temp / 1e5, 0.5));
        gamma_HeI *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
        gamma_HeI *= All.HubbleParam * All.HubbleParam;

        gamma_HeII = 5.68e-12 * pow(temp, 0.5) * exp(-631515 / temp) / (1.0 + pow(temp / 1e5, 0.5));
        gamma_HeII *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
        gamma_HeII *= All.HubbleParam * All.HubbleParam;

        /* alpha_B recombination coefficient */
        alpha_HeII = 1.5e-10 * pow(temp, -0.6353);
        alpha_HeII *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
        alpha_HeII *= All.HubbleParam * All.HubbleParam;

        alpha_HeIII = 3.36e-10 * pow(temp, -0.5) * pow(temp / 1e3, -0.2) / (1.0 + pow(temp / 1e6, 0.7));
        alpha_HeIII *= All.UnitTime_in_s / All.UnitLength_in_cm / All.UnitLength_in_cm / All.UnitLength_in_cm;
        alpha_HeIII *= All.HubbleParam * All.HubbleParam;
#endif

        double densphot = 0.0;
        for(j = 0; j < RT_N_DIR; j++)
          {
            SphP[i].DensPhot[j] = SphP[i].Photons[j] / SphP[i].Volume;
            densphot += SphP[i].Photons[j] / SphP[i].Volume;
          }

        nHI = SphP[i].nHI;
        nHII = SphP[i].nHII;
        nGamma = densphot * All.cf_a3inv / nH;

        int calculated_flag = 0;

        double d_nGamma = -c_light * sigma * nH * dt * nHI * nGamma;    /* linear estimate of the photo abundance change */

        if(fabs(d_nGamma) <= 0.05 * nGamma)     /* we are nearly in photoionization equilibrium, and/or the photon density is so large that it should change little */
          {
            nGamma_new = nGamma + d_nGamma;

            /* let's solve implicitely for new hydrogen abundance */

            double A = alpha * nH * dt;
            double B = c_light * sigma * nH * nGamma * dt;
            double D = gamma * nH * dt;

            if(A != 0)
              {
                double b = 2 * A + B + D + 1;
                double det = b * b - 4 * (A + D) * (A + nHI);
                if(det >= 0)
                  {
                    nHI_new = (b - sqrt(det)) / (2 * (A + D));

                    if(fabs(nHI_new - nHI) < 0.25 * nHI)        /* only accept if change is sufficiently small */
                      calculated_flag = 1;
                  }
              }
            else
              {
                if(CLIGHT != 0)
                  terminate("not implemented");

                /* if there are no recombinations (A=0), we can directly do an implicit solve */
                double B = c_light * sigma * nH * dt;
                double b = 1 + B * (nGamma - nHI);
                double det = b * b + 4 * B * nHI;
                if(B > 0 && det > 0)
                  {
                    nHI_new = (-b + sqrt(det)) / (2 * B);
                    nGamma_new = nGamma + nHI_new - nHI;
                  }
                else
                  {
                    nHI_new = nHI;
                    nGamma_new = nGamma;
                  }
                calculated_flag = 1;
              }
          }
        if(calculated_flag == 0 && nGamma < 1.0e-2 * nHI)       /* here photoionization cannot change the neutral fraction much */
          {
            /* first, calculate an implicit solution for the new gamma, neglecting a change in nHI */
            fac = c_light * sigma * nH * dt;    /* photo ionization factor */
            if(fac > 0)
              nGamma_new = (-(1 + fac * (1 - nHII - nGamma)) + sqrt(pow(1 + fac * (1 - nHII - nGamma), 2) + 4 * fac * nGamma)) / (2 * fac);
            else
              nGamma_new = nGamma;

            /* now calculate the new neutral and ionized hydrogen fraction */
            double A = alpha * nH * dt; /* recombination factor */
            double B = gamma * nH * dt; /* collisional ionization factor */

            if((A + B) != 0)
              {
                double det = (1 - B) * (1 - B) + 4 * (A + B) * (nHII + nGamma - nGamma_new);
                if(det > 0)
                  {
                    nHII_new = (-1 + B + sqrt(det)) / (2 * (A + B));
                    nHI_new = 1 - nHII_new;

                    if(fabs(nHI_new - nHI) < 0.25 * nHI)        /* only accept if change is sufficiently small */
                      calculated_flag = 1;
                  }
              }
          }

        if(calculated_flag == 0)        /* if no accurate solution has been found yet, integrate the problem with an ODE solver */
          {
            count_part_ode++;

            /* let's here explicitely integrate the (stiff) differential equation */
            const gsl_odeiv_step_type *T = gsl_odeiv_step_rkf45;
            gsl_odeiv_step *s = gsl_odeiv_step_alloc(T, 2);
            gsl_odeiv_control *c = gsl_odeiv_control_y_new(1e-5, 0.0);
            gsl_odeiv_evolve *e = gsl_odeiv_evolve_alloc(2);

            double coeff[3] = { alpha * nH, c_light * sigma * nH, gamma * nH };
            gsl_odeiv_system sys = { rt_rate_ODEs, NULL, 2, coeff };

            double t = 0.0, t1 = dt;
            double h = 0.25 * dt;
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

        for(j = 0; j < RT_N_DIR; j++)
          {
            if(densphot > 0)
              SphP[i].Photons[j] = nH / All.cf_a3inv * nGamma_new * SphP[i].Volume * SphP[i].DensPhot[j] / densphot;

            if(SphP[i].Photons[j] < 0)
              {
                printf("oops %g %g %g  \n", SphP[i].Photons[j], nGamma_new, nGamma);
                SphP[i].Photons[j] = 0.0;
              }

            SphP[i].DensPhot[j] = SphP[i].Photons[j] / SphP[i].Volume;
          }

        /* store the new abundance results */

        SphP[i].nHI = nHI_new;
        SphP[i].nHII = nHII_new;
        SphP[i].n_elec = SphP[i].nHII;

        mean_nHI += SphP[i].nHI * SphP[i].Volume;


#ifdef RT_INCLUDE_HE
        nHeIII = 1 - (SphP[i].nHeI + SphP[i].nHeII);

        SphP[i].n_elec += SphP[i].nHeII * (1.0 - HYDROGEN_MASSFRAC) / (4.0 * HYDROGEN_MASSFRAC);
        SphP[i].n_elec += 2.0 * nHeIII * (1.0 - HYDROGEN_MASSFRAC) / (4.0 * HYDROGEN_MASSFRAC);

        D = dt * gamma_HeII * nH;
        E = dt * alpha_HeIII * nH;
        F = dt * gamma_HeI * nH;
        J = dt * alpha_HeII * nH;

        nHeII = SphP[i].nHeII + F * SphP[i].n_elec - ((F * SphP[i].n_elec - E * SphP[i].n_elec) / (1.0 + E * SphP[i].n_elec) * nHeIII);

        nHeII /= 1.0 + F * SphP[i].n_elec + D * SphP[i].n_elec + J * SphP[i].n_elec + ((F * SphP[i].n_elec - E * SphP[i].n_elec) / (1.0 + E * SphP[i].n_elec) * D * SphP[i].n_elec);

        if(nHeII < 0)
          nHeII = 0.0;

        if(nHeII > 1)
          nHeII = 1.0;

        nHeIII = (nHeIII + D * nHeII * SphP[i].n_elec) / (1.0 + E * SphP[i].n_elec);

        if(nHeIII < 0)
          nHeIII = 0.0;

        if(nHeIII > 1)
          nHeIII = 1.0;


        SphP[i].n_elec = SphP[i].nHII;
        SphP[i].n_elec += nHeII * (1.0 - HYDROGEN_MASSFRAC) / (4.0 * HYDROGEN_MASSFRAC);
        SphP[i].n_elec += 2.0 * nHeIII * (1.0 - HYDROGEN_MASSFRAC) / (4.0 * HYDROGEN_MASSFRAC);

        SphP[i].nHeII = nHeII;

        SphP[i].nHeI = 1.0 - (SphP[i].nHeII + nHeIII);

        if(SphP[i].nHeI < 0 || SphP[i].nHeI > 1 || isnan(SphP[i].nHeI))
          {
            terminate("ERROR, wrong value of nHeI %g\n", SphP[i].nHeI);
          }

        mean_nHeI += SphP[i].nHeI * SphP[i].Volume;
#endif

      }

  MPI_Reduce(&count_part_ode, &tot_count_part_ode, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&count_part_ode_iter, &tot_count_part_ode_iter, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  MPI_Reduce(&mean_nHI, &mean_nHI_all, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#ifdef RT_INCLUDE_HE
  MPI_Reduce(&mean_nHeI, &mean_nHeI_all, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif

  if(ThisTask == 0)
    {
      fprintf(FdRad, "%g %g ", All.Time, mean_nHI_all / All.BoxSize / All.BoxSize / All.BoxSize);
#ifndef RT_INCLUDE_HE
      fprintf(FdRad, "\n");
#else
      fprintf(FdRad, "%g\n", mean_nHeI_all / All.BoxSize / All.BoxSize / All.BoxSize);
#endif
      myflush(FdRad);
    }

  if(ThisTask == 0)
    {
      printf
        ("Ionization balance solver: fraction of particle requiring ODE: %g, average number of iterations=%g\n",
         tot_count_part_ode / All.TotNumGas, tot_count_part_ode_iter / (tot_count_part_ode + 1.0e-8));
    }

}

/* this function sets up simple initial conditions for a single source in a uniform field of gas with constant density*/
void rt_set_simple_inits(void)
{
  int i, j;
  double nHeIII;

  mpi_printf("initialising RT variables...\n");

  for(i = 0; i < NumGas; i++)
    if(P[i].Type == 0)
      {
        SphP[i].nHII = 1e-8;
        SphP[i].nHI = 1 - SphP[i].nHII;
        SphP[i].n_elec = SphP[i].nHII;

#ifdef RT_INCLUDE_HE
        nHeIII = 1e-8;
        SphP[i].nHeII = 1e-8;
        SphP[i].nHeI = 1.0 - (SphP[i].nHeII + nHeIII);
        SphP[i].n_elec += (SphP[i].nHeII + 2.0 * nHeIII) * (1.0 - HYDROGEN_MASSFRAC) / 4.0 / HYDROGEN_MASSFRAC;
#endif

        for(j = 0; j < RT_N_DIR; j++)
          {
            SphP[i].Photons[j] = 0.0;
            SphP[i].DensPhot[j] = 0.0;
#ifndef RT_HEALPIX_NSIDE
            SphP[i].SourceID[j] = -1;
            SphP[i].SourcePos[j][0] = 0.0;
            SphP[i].SourcePos[j][1] = 0.0;
            SphP[i].SourcePos[j][2] = 0.0;
#endif
          }

      }
  mpi_printf("done initialising RT variables...\n");
}


#endif
