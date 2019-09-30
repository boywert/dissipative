/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/anisotropic_RT/RT_chem_ps2009.c
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

/* Solves the chemical network (H and He) following Petkova & Springel 2009. This
 * means that absorptions are computed during CG-method, this only updates the H and He
 * ratios, not the photon density per cell SphP.n_gamma. 
*/
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "../allvars.h"
#include "../proto.h"
#include "RT.h"

#if defined(MRT) && defined(MRT_CHEMISTRY_PS2009)



#ifndef MRT_MULTI_FREQUENCY
void mrt_update_chemistry_ps2009(void)
{
  int idx, i;
  double nH, temp, molecular_weight, rho;
  double nHII, nHI, nHI_new, nHII_new;
  double dt, dtime, a3inv, hubble_a, c_light;
  double A, B, CC;
  double n_gamma;
  double alpha_HII, gamma_HI;
  double fac;

#ifdef MRT_INCLUDE_HE
  double alpha_HeII, alpha_HeIII, gamma_HeI, gamma_HeII;
  double nHeII, nHeIII;
  double D, E, F, G, J, L;
  double y_fac;
#endif

#ifndef MRT_MULTI_FREQUENCY
  //  fac = 1.0 / All.UnitLength_in_cm / All.UnitLength_in_cm * All.HubbleParam * All.HubbleParam;
  // mrt_sigma_HI[0] = 6.3e-18 * fac;
  
#endif

  fac = All.UnitTime_in_s / pow(All.UnitLength_in_cm, 3) * All.HubbleParam * All.HubbleParam;

  if(All.ComovingIntegrationOn)
    {
      hubble_a = hubble_function(All.Time);
      a3inv = 1.0 / All.Time / All.Time / All.Time;
    }
  else
    hubble_a = a3inv = 1.0;

  c_light = CLIGHT / All.UnitVelocity_in_cm_per_s;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      dt = (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval ;

      if(All.ComovingIntegrationOn)
        dtime = dt / hubble_a;
      else
        dtime = dt;

      rho = SphP[i].Density * a3inv;

      nH = HYDROGEN_MASSFRAC * rho / PROTONMASS * All.UnitMass_in_g / All.HubbleParam;

      //molecular_weight = 4 / (1 + 3 * HYDROGEN_MASSFRAC + 4 * HYDROGEN_MASSFRAC * SphP[i].Ne);

      molecular_weight = 1.0 ;

      temp = SphP[i].Utherm * GAMMA_MINUS1 * (molecular_weight * PROTONMASS / All.UnitMass_in_g * All.HubbleParam) / (BOLTZMANN / All.UnitEnergy_in_cgs * All.HubbleParam);

      /* collisional ionization rate */
      gamma_HI = 5.85e-11 * sqrt(temp) * exp(-157809.1 / temp) / (1.0 + sqrt(temp / 1e5)) * fac;
    
      /* alpha_B recombination coefficient */
      alpha_HII = 2.59e-13 * pow(temp / 1e4, -0.7) * fac;

      n_gamma = SphP[i].Cons_DensPhot[0] * 1e63 / nH / SphP[i].Volume ;


      double sigma =  mrt_sigma_HI[0] ;

      /* number of photons should be positive */
      if(n_gamma < 0 || isnan(n_gamma))
        {
          printf("NEGATIVE n_gamma: %g %d %d \n", n_gamma, i, ThisTask);
          printf("n_gamma %g mass %g a3inv %g \n", SphP[i].Cons_DensPhot[0], P[i].Mass, a3inv);
          terminate("111");
        }

      A = dtime * gamma_HI * nH * SphP[i].Ne ;
      B = dtime * c_light * nH * n_gamma * sigma ;
      CC = dtime * alpha_HII * nH * SphP[i].Ne ;

      /* semi-implicit scheme for ionization */
      nHII = SphP[i].HII + B + A;

      nHII /= 1.0 + B + CC + A;

      if(nHII < 0 || nHII > 1 || isnan(nHII))
        {
          print_particle_info(i);
          printf("ERROR nHII %g \n", nHII);
          terminate("333");
        }
      
      
      //SphP[i].Ne = nHII;        // --need to check what happens with cooling's Ne??

      //SphP[i].HII = nHII;

      //SphP[i].HI = 1.0 - nHII;


      double nH_times_volume = P[i].Mass ;
      SphP[i].nHI = (1.0 - nHII) * nH_times_volume ;
      SphP[i].nHII = nHII * nH_times_volume ;
      SphP[i].ne = nHII * nH_times_volume ;      

      double nGamma_new = n_gamma/(1.0 + dtime * c_light * sigma * nH * (1.0 - nHII)) ;

      double ratio = (nGamma_new / 1e63 * nH * SphP[i].Volume) / SphP[i].Cons_DensPhot[0] ;
      //SphP[i].Cons_DensPhot_absorbed[0] *= (1.0 - ratio) ;
      SphP[i].Cons_DensPhot[0] *= ratio ;

      int num1 ;
      for(num1=0;num1<3;num1++)
          SphP[i].Cons_RT_F[0][num1] *= ratio ;

#ifdef MRT_INCLUDE_HE
      /* collisional ionization rate */
      gamma_HeI = 2.38e-11 * sqrt(temp) * exp(-285335.4 / temp) / (1.0 + sqrt(temp / 1e5)) * fac;
      gamma_HeII = 5.68e-12 * sqrt(temp) * exp(-631515 / temp) / (1.0 + sqrt(temp / 1e5)) * fac;

      /* alpha_B recombination coefficient */
      alpha_HeII = 1.5e-10 * pow(temp, -0.6353) * fac;
      alpha_HeIII = 3.36e-10 / sqrt(temp) * pow(temp / 1e3, -0.2) / (1.0 + pow(temp / 1e6, 0.7)) * fac;

      SphP[i].Ne += SphP[i].HeII + 2.0 * SphP[i].HeIII;

      D = dtime * gamma_HeII * nH * SphP[i].Ne;
      E = dtime * alpha_HeIII * nH * SphP[i].Ne;
      F = dtime * gamma_HeI * nH * SphP[i].Ne;
      J = dtime * alpha_HeII * nH * SphP[i].Ne;
      G = 0.0;
      L = 0.0;

      y_fac = (1.0 - HYDROGEN_MASSFRAC) / 4.0 / HYDROGEN_MASSFRAC;

      nHeII = SphP[i].HeII / y_fac;
      nHeIII = SphP[i].HeIII / y_fac;

      nHeII = nHeII + G + F - ((G + F - E) / (1.0 + E)) * nHeIII;

      nHeII /= 1.0 + G + F + D + J + ((G + F - E) / (1.0 + E)) * (D + L);

      if(nHeII < 0 || nHeII > 1 || isnan(nHeII))
        {
          printf("ERROR nHeII %g \n", nHeII);
          terminate("333");
        }

      nHeIII = nHeIII + (D + L) * nHeII;

      nHeIII /= 1.0 + E;

      if(nHeIII < 0 || nHeIII > 1 || isnan(nHeIII))
        {
          printf("ERROR nHeIII %g \n", nHeIII);
          terminate("333");
        }

      //SphP[i].Ne = SphP[i].HII + nHeII + 2.0 * nHeIII;

      nHeII *= y_fac;
      nHeIII *= y_fac;

      //      double nH_times_volume = P[i].Mass ;
      SphP[i].nHeII = nHeII * nH_times_volume ;
      SphP[i].nHeIII = nHeIII * nH_times_volume ;

      double nHeI = y_fac - nHeII - nHeII ;
      double ne = nHII + nHeII + 2.0 * nHeIII

       
	//SphP[i].Ne = SphP[i].HII + nHeII + 2.0 * nHeIII;

	//SphP[i].HeII = nHeII;
	//SphP[i].HeIII = nHeIII;

      //SphP[i].HeI = y_fac - SphP[i].HeII - SphP[i].HeIII;

      if(nHeI < 0)
        nHeI = 0.0;

      if(nHeI > y_fac)
        nHeI = y_fac;

      SphP[i].nHeI = nHeI * nH_times_volume ;
      SphP[i].ne = ne * nH_times_volume ;
#endif
    }

}

#else

/*---------------------------------------------------------------------*/
/* if the multi-frequency scheme is used*/
/*---------------------------------------------------------------------*/
void mrt_update_chemistry_ps2009(void)
{
  int idx, i, j;
  double nH, temp, molecular_weight, rho;
  double nHII, c_light, n_gamma;
  double dt, dtime, a3inv, hubble_a;
  double A, B, CC;
  double alpha_HII, gamma_HI;
  double fac;
  double k_HI;

#ifdef MRT_INCLUDE_HE
  double alpha_HeII, alpha_HeIII, gamma_HeI, gamma_HeII;
  double nHeII, nHeIII;
  double D, E, F, G, J, L;
  double k_HeI, k_HeII;
  double y_fac;
#endif

  fac = All.UnitTime_in_s / pow(All.UnitLength_in_cm, 3) * All.HubbleParam * All.HubbleParam;

  c_light = CLIGHT / All.UnitVelocity_in_cm_per_s;

  if(All.ComovingIntegrationOn)
    {
      hubble_a = hubble_function(All.Time);
      a3inv = All.Time / All.Time / All.Time;
    }
  else
    hubble_a = a3inv = 1.0;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      /* get the photo-ionization rates */
      k_HI = 0.0;
#ifdef MRT_INCLUDE_HE
      k_HeI = k_HeII = 0.0;
#endif


      rho = SphP[i].Density * a3inv;
      
      nH = HYDROGEN_MASSFRAC * rho / PROTONMASS * All.UnitMass_in_g / All.HubbleParam;
      
      molecular_weight = 4 / (1 + 3 * HYDROGEN_MASSFRAC + 4 * HYDROGEN_MASSFRAC * SphP[i].Ne);
      
      temp = SphP[i].Utherm * GAMMA_MINUS1 * (molecular_weight * PROTONMASS / All.UnitMass_in_g * All.HubbleParam) / (BOLTZMANN / All.UnitEnergy_in_cgs * All.HubbleParam);

      //temp = 1e4;
      //molecular_weight = 1.0 ;

      for(j = 0; j < UV_BINS; j++)
        {
	  //          n_gamma = SphP[i].n_gamma[j] / P[i].Mass * a3inv;

	  n_gamma = SphP[i].Cons_DensPhot[j] * 1e63 / SphP[i].Volume ;

	  //          if(nu[j] >= 13.6)
	  k_HI += c_light * mrt_sigma_HI[j] * n_gamma ;

#ifdef MRT_INCLUDE_HE
          //if(nu[j] >= 24.6)
	  k_HeI += c_light * mrt_sigma_HeI[j] * n_gamma ;

          //if(nu[j] >= 54.4)
	  k_HeII += c_light * mrt_sigma_HeII[j] * n_gamma ;
#endif
        }

      dt = (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval ;

      if(All.ComovingIntegrationOn)
        dtime = dt / hubble_a;
      else
        dtime = dt;


      //#if defined(OTVET_TEST_SST) && !defined(COOLING)
      //temp = 1e4;
      //#endif

      /* collisional ionization rate */
      gamma_HI = 5.85e-11 * sqrt(temp) * exp(-157809.1 / temp) / (1.0 + sqrt(temp / 1e5)) * fac;

      /* alpha_B recombination coefficient */
      alpha_HII = 2.59e-13 * pow(temp / 1e4, -0.7) * fac;

      A = dtime * gamma_HI * nH * SphP[i].Ne;
      B = dtime * k_HI;
      CC = dtime * alpha_HII * nH * SphP[i].Ne;

      /* semi-implicit scheme for ionization */
      nHII = SphP[i].HII + B + A;

      nHII /= 1.0 + B + CC + A;


      //            if(nHII < 0)
      //nHII = 0.0 ;
      
      if(nHII < 1e-20 || nHII > (1+1e-5) || isnan(nHII))
        {
          print_particle_info(i);
          printf("ERROR nHII %g %g %g %g\n", nHII, SphP[i].HI, SphP[i].HII, SphP[i].Ne);
	  printf("A, B, CC = %g \t %g \t %g \n", A, B, CC) ;
	  printf("Cons dens phot %g \n", SphP[i].Cons_DensPhot[0]) ;
          terminate("333");
        }

      //      SphP[i].Ne = nHII;

      //SphP[i].HII = nHII;

      //SphP[i].HI = 1.0 - nHII;


      double nH_times_volume = P[i].Mass ;
      SphP[i].nHI = (1.0 - nHII) * nH_times_volume ;
      SphP[i].nHII = nHII * nH_times_volume ;
      SphP[i].ne = nHII * nH_times_volume ;      


#ifdef MRT_INCLUDE_HE
      /* collisional ionization rate */
      gamma_HeI = 2.38e-11 * sqrt(temp) * exp(-285335.4 / temp) / (1.0 + sqrt(temp / 1e5)) * fac;
      gamma_HeII = 5.68e-12 * sqrt(temp) * exp(-631515 / temp) / (1.0 + sqrt(temp / 1e5)) * fac;

      /* alpha_B recombination coefficient */
      alpha_HeII = 1.5e-10 * pow(temp, -0.6353) * fac;
      alpha_HeIII = 3.36e-10 / sqrt(temp) * pow(temp / 1e3, -0.2) / (1.0 + pow(temp / 1e6, 0.7)) * fac;

      //SphP[i].Ne += SphP[i].HeII + 2.0 * SphP[i].HeIII;
      SphP[i].Ne = nHII + SphP[i].HeII + 2.0 * SphP[i].HeIII;

      D = dtime * gamma_HeII * nH * SphP[i].Ne;
      E = dtime * alpha_HeIII * nH * SphP[i].Ne;
      F = dtime * gamma_HeI * nH * SphP[i].Ne;
      J = dtime * alpha_HeII * nH * SphP[i].Ne;
      G = dtime * k_HeI;
      L = dtime * k_HeII;

      y_fac = (1.0 - HYDROGEN_MASSFRAC) / 4.0 / HYDROGEN_MASSFRAC;

      nHeII = SphP[i].HeII / y_fac;
      nHeIII = SphP[i].HeIII / y_fac;

      nHeII = nHeII + G + F - ((G + F - E) / (1.0 + E)) * nHeIII;

      nHeII /= 1.0 + G + F + D + J + L + ((G + F - E) / (1.0 + E)) * (D + L);

      if(isnan(nHeII))
        {
          printf("ERROR nHeII %g %g %g\n", nHeII, temp, k_HeI);
          terminate("333");
	}
      if(nHeII<0)
	nHeII = 0.0 ;
      if(nHeII>1)
	nHeII = 1.0 ;


      nHeIII = nHeIII + (D + L) * nHeII;

      nHeIII /= 1.0 + E;

      if(isnan(nHeIII))
        {
          printf("ERROR nHeIII %g %g %g\n", nHeIII, temp, k_HeII);
          terminate("333");
        }

      if(nHeIII<0)
	nHeIII = 0.0 ;
      if(nHeIII>1)
	nHeIII = 1.0 ;

      double nHeI = 1.0 - nHeII - nHeIII ;
      nHeII *= y_fac ;
      nHeIII *= y_fac ;
      nHeI *= y_fac ;


      //double nHeI = y_fac - nHeII - nHeII ;
      double ne = nHII + nHeII + 2.0 * nHeIII ;

       
      if(nHeI < 0)
        nHeI = 0.0;

      if(nHeI > y_fac)
        nHeI = y_fac;

      SphP[i].nHeI = nHeI * nH_times_volume ;
      SphP[i].ne = ne * nH_times_volume ;
      SphP[i].nHeII = nHeII * nH_times_volume ;
      SphP[i].nHeIII = nHeIII * nH_times_volume ;


      //      SphP[i].Ne = SphP[i].HII + nHeII + 2.0 * nHeIII;

      //SphP[i].HeII = nHeII;
      //SphP[i].HeIII = nHeIII;

      //SphP[i].HeI = y_fac - SphP[i].HeII - SphP[i].HeIII;

      //if(SphP[i].HeI < 0)
      //SphP[i].HeI = 0.0;

      //if(SphP[i].HeI > y_fac)
      //SphP[i].HeI = y_fac;
#endif
      double KK = dtime * c_light * nH ;
      double ratio ;
      for(j = 0; j < UV_BINS; j++)
	{
	  double sum_KK = mrt_sigma_HI[j] * SphP[i].HI 
#ifdef MRT_INCLUDE_HE
	    + mrt_sigma_HeI[j] * SphP[i].HeI + mrt_sigma_HeII[j] * SphP[i].HeII
#endif
	    ;   
	  ratio = 1.0 / (1.0 + KK*sum_KK) ;
	  SphP[i].Cons_DensPhot_absorbed[j] *= (1.0 - ratio) ;
	  SphP[i].Cons_DensPhot[j] *= ratio ;

	  int num1 ;
	  for(num1=0;num1<3;num1++)
	    SphP[i].Cons_RT_F[j][num1] *= ratio ;
	}
    }
}
#endif
/*
void mrt_write_stats(void)
{
  int i;
  double rho, a3inv;
  double total_nHI, total_V, total_nHI_all, total_V_all;
  total_nHI = 0.0;
  total_V = 0.0;

#ifndef MRT_MULTI_FREQUENCY
  double total_ng, total_ng_all, n_gamma;
  total_ng = 0.0;
#endif

#ifdef MRT_INCLUDE_HE
  double total_nHeI, total_nHeI_all;
  double total_nHeII, total_nHeII_all;
  total_nHeI = total_nHeII = 0.0;
#endif

  if(All.ComovingIntegrationOn)
    a3inv = All.Time / All.Time / All.Time;
  else
    a3inv = 1.0;

  for(i = 0; i < NumGas; i++)
    if(P[i].Type == 0)
      {
        rho = SphP[i].Density * a3inv;

#ifndef MRT_MULTI_FREQUENCY
        n_gamma = SphP[i].n_gamma[0] / P[i].Mass * a3inv;
        total_ng += n_gamma / 1e53 * P[i].Mass / rho;
#endif

#ifdef MRT_INCLUDE_HE
        total_nHeI += SphP[i].HeI * P[i].Mass / rho;
        total_nHeII += SphP[i].HeII * P[i].Mass / rho;
#endif
        total_nHI += SphP[i].HI * P[i].Mass / rho;
        total_V += P[i].Mass / rho;
      }

  MPI_Allreduce(&total_nHI, &total_nHI_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&total_V, &total_V_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#ifndef MRT_MULTI_FREQUENCY
  MPI_Allreduce(&total_ng, &total_ng_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
#ifdef MRT_INCLUDE_HE
  MPI_Allreduce(&total_nHeI, &total_nHeI_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&total_nHeII, &total_nHeII_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  if(ThisTask == 0)
    {
      if(All.Time == All.TimeBegin)
        {
          fprintf(FdMRT, " MRT: time, nHI");
#ifndef MRT_MULTI_FREQUENCY
          fprintf(FdMRT, ", n_gamma");
#endif
#ifdef MRT_INCLUDE_HE
          fprintf(FdMRT, ", nHeI, nHeII \n");
#else
          fprintf(FdMRT, "\n");
#endif
        }

      fprintf(FdMRT, "%g %g ", All.Time, total_nHI_all / total_V_all);
#ifndef MRT_MULTI_FREQUENCY
      fprintf(FdMRT, "%g ", total_ng_all / total_V_all);
#endif
#ifdef MRT_INCLUDE_HE
      fprintf(FdMRT, "%g %g\n", total_nHeI_all / total_V_all, total_nHeII_all / total_V_all);
#else
      fprintf(FdMRT, "\n");
#endif
      fflush(FdMRT);
    }

}
*/
#endif
