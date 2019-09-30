/*! * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/MRT/RT.c
 * \date        MM/YYYY
 * \author      Rahul Kannan
 * \brief        
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 */


/*! \file RT_semi_HYPRE.c
 *  \brief main driver for an moment based RT with the VET formalism
 *
 *  
 */



#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

#include "../allvars.h"
#include "../proto.h"
#include "../domain.h"
#include "../voronoi.h"



#ifdef MRT

 
void mrt_update_chemistry()
{
#ifndef MRT_NO_UV
#ifdef MRT_CHEMISTRY_PS2011
  mrt_update_chemistry_ps2011() ;
  mpi_printf("RT: Updated Chemistry with PS2011\n") ;
#endif
#ifdef MRT_CHEMISTRY_PS2009
  mrt_update_chemistry_ps2009() ;
  mpi_printf("RT: Updated Chemistry with PS2009\n") ;
#endif
#endif

}



void init_RT()
{
  mpi_printf("|||||||||||%d\t%d\t%d\t%d|||||||||\n\n\n\n", UV_BINS, IR_BINS, MRT_BINS, COUNT_GRAD_MRT) ;
#ifdef MRT_CONSTANT_KAPPA
  c_internal_units = 1.0 ;
#else
  c_internal_units = CLIGHT / All.UnitVelocity_in_cm_per_s ;
#endif

#ifdef MRT_IR_LTE
  radiation_constant = RADIATION_CONSTANT * pow(All.UnitLength_in_cm,3) / All.UnitEnergy_in_cgs ;
  mpi_printf("RT: INIT : Radiation Constant = %le\n", radiation_constant) ;
#endif

#ifdef MRT_RIEMANN_HLLE
  readin_hlle_eingenvalues() ;
#endif

#ifdef MRT_IR_GRAIN_KAPPA
  mpi_printf("MRT: Reading grain kappa data.\n");
  read_grain_kappa_data();
#endif

  mpi_printf("MRT: CLIGHT = %le\n", c_internal_units);
}


void add_source_fluxes() /*Correct this*/
{
  mpi_printf("RT: Adding source fluxes....\n") ;

#ifdef MRT_SOURCES

#ifdef MRT_STARS
  do_ionizing_stellar_sources();
#endif

#ifdef MRT_BH
  do_ionizing_blackhole_sources();
#endif

#else
  //#ifndef MRT_NO_UV
  double dt = 0.0 ;
    
 
  int i, idx ;

  double X, Y, Z, radius ;

  int Nincells = 0 ;
  int nincells = 0 ;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i > 0)
        break;
    }

  dt = (All.HighestActiveTimeBin ? (((integertime) 1) << All.HighestActiveTimeBin) : 0) * All.Timebase_interval;

  double xcm = 0.0 ;
  double ycm = 0.0 ;
  double zcm = 0.0 ;

  double Xcm, Ycm, Zcm ;



  for(i=0;i<NumGas;i++)
    {
      
      
      //if(isnan(SphP[i].Cons_DensPhot))
      //	terminate("Before WHAT!!!!\n") ;

      X = (P[i].Pos[0] - All.BoxSize/2.0) ; //- 0.6*c_internal_units*All.Time ;

      Y = P[i].Pos[1] - All.BoxSize/2.0 ;
#ifndef TWODIMS
      Z = P[i].Pos[2] - All.BoxSize/2.0 ;
#endif

#ifndef TWODIMS
      radius = sqrt(X*X + Y*Y + Z*Z) ;
#else
      radius = sqrt(X*X + Y*Y) ;
#endif

      //      if(radius<0.2)
      if(P[i].Pos[0]>20 && P[i].Pos[0]<25 && P[i].Pos[1]>20 && P[i].Pos[1]<25)	
	{
	  //	  printf("Entered Here\n") ;
	  nincells += 1 ;
	  xcm += P[i].Pos[0] ;
	  ycm += P[i].Pos[1] ;
	  zcm += P[i].Pos[2] ;
	  // printf("X, Y, ID = %g %g %d\n", P[i].Pos[0], P[i].Pos[1], P[i].ID) ;
	}
      //      terminate("NNNN\n") ;
    }

  MPI_Allreduce(&nincells, &Nincells, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&xcm, &Xcm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;
  MPI_Allreduce(&ycm, &Ycm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;
  MPI_Allreduce(&zcm, &Zcm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;



  mpi_printf("RT: %d number of cells injected with photons\t", Nincells) ;
  mpi_printf("RT: Center of mass of cells x = %g, y = %g, z = %g\n", Xcm/Nincells, Ycm/Nincells, Zcm/Nincells) ;

for(i=0; i<NumGas; i++)
    {
      X = (P[i].Pos[0] - All.BoxSize/2.0) ; // - 0.6*c_internal_units*All.Time ;
      Y = P[i].Pos[1] - All.BoxSize/2.0 ;
#ifndef TWODIMS
      Z = P[i].Pos[2] - All.BoxSize/2.0 ;
#endif

#ifndef TWODIMS
      radius = sqrt(X*X + Y*Y + Z*Z) ;
#else
      radius = sqrt(X*X + Y*Y) ;
#endif
      //      if(radius<0.2)
      if(P[i].Pos[0]>20 && P[i].Pos[0]<25 && P[i].Pos[1]>20 && P[i].Pos[1]<25)	 
	{
	  double Nphotons = 10; //1.6e37 * All.UnitTime_in_s / All.UnitEnergy_in_cgs ;
	  for(int num1=0;num1<UV_BINS;num1++)
	    {
	      //	      if(All.Time==0.0)
	      //	      SphP[i].Cons_DensPhot[num1] = Nphotons ;
	      //	      //	      //	      //	      double Nphotons = 1e50  ;
	      //	      SphP[i].Cons_DensPhot[num1] += lum[num1] * dt * All.UnitTime_in_s * Nphotons / ((double) (Nincells)) ;   

	      //SphP[i].DensPhot[num1] = SphP[i].Cons_DensPhot[num1] / SphP[i].Volume ; 
	      /* double nx = P[i].Pos[0] - All.BoxSize/2.0 ;
	      double ny = P[i].Pos[1] - All.BoxSize/2.0 ;
	      double nz = P[i].Pos[2] - All.BoxSize/2.0 ;

	      double modn = sqrt(nx*nx + ny*ny + nz*nz) ;
	      if(modn == 0.0)
	      	modn = 1.0 ;
	      nx = nx/modn ;
	      ny = ny/modn ;
	      nz = nz/modn ;
	      
	      SphP[i].Cons_RT_F[num1][0] +=  0.99999999 * c_internal_units * SphP[i].Cons_DensPhot[num1] ;
	      SphP[i].Cons_RT_F[num1][1] += 0.0 ;
	      SphP[i].Cons_RT_F[num1][2] += 0.0 ;

	      /*SphP[i].RT_F[num1][0] = SphP[i].Cons_RT_F[num1][0] / SphP[i].Volume ; 
	      SphP[i].RT_F[num1][1] = SphP[i].Cons_RT_F[num1][1] / SphP[i].Volume ; 
	      SphP[i].RT_F[num1][2] = SphP[i].Cons_RT_F[num1][2] / SphP[i].Volume ; */
	    }
	}
	      //	      printf("RT: Added flux\n") ;
	      /*double total_energy = SphP[i].Cons_DensPhot[num1] + SphP[i].Trapped_Cons_DensPhot[num1-UV_BINS] + dt * All.UnitTime_in_s * Nphotons / ((double) (Nincells)) ;
	      // SphP[i].Cons_DensPhot[num1] += dt * All.UnitTime_in_s * Nphotons / ((double) (Nincells)) ;   
	      double rad = 2.0*get_cell_radius(i) ;
	      double tau = SphP[i].KappaIR[num1-UV_BINS] * rad ;
	      //if(dt>0.0)
	      //terminate("rad = %g \t kappa = %g \t tau = %g \n", rad, tau, SphP[i].KappaIR[num1-UV_BINS]) ;
	      double Es = (1.0 - exp(-2.0/(3.0*tau)))*total_energy ;
	      double Et = total_energy - Es ;
	      SphP[i].Trapped_Cons_DensPhot[num1-UV_BINS] = Et ;
	      SphP[i].Cons_DensPhot[num1] = Es ;*/
      //	    }
      //	}
	    /*            for(int num1=0;num1<MRT_BINS;num1++)
	{
	  SphP[i].Cons_RT_F[num1][0] = 0.99*c_internal_units*SphP[i].Cons_DensPhot[num1] ;
	  SphP[i].Cons_RT_F[num1][1] = 0.0 ;
	  SphP[i].Cons_RT_F[num1][2] = 0.0;
	  }*/
      // }
//#endif
    }
#endif
}

#ifdef MRT_TIME_EXTRAPOLATION
void calculate_div_F(struct state *dl, struct state *st)
{
  /*div F = Trace (grad  F)*/
  /*grad F = N grad(F/N) + 1/N F grad N*/

  struct rt_grad_data *rtgrad = st->rtgrad ;
  
  int num1, j ;
  for(num1=0;num1<MRT_BINS;num1++)
    {
      double trace = 0.0 ;
      double top = st->DensPhot[num1] ;
      double bottom = 1.0/st->DensPhot[num1] ;
      {
	for(j=0;j<3;j++)
	  trace += top * rtgrad->dFN[num1][j][j] + bottom * st->RT_F[num1][j] * rtgrad->dDensPhot[num1][j] ;
      }
      dl->divF[num1] = trace ;
      //      if(trace != 0.0)
      //	terminate("Working\n") ;
    }
  return ;
}
#endif


#ifdef MRT_IR
#ifdef MRT_IR_GRAIN_KAPPA
int IR_N_pts;
double *IR_logT, *IR_logkappaP, *IR_logkappaR;
gsl_interp_accel *accIR_kappaP;
gsl_interp_accel *accIR_kappaR;
gsl_spline *splineIR_kappaP;
gsl_spline *splineIR_kappaR;

/* Read in tabulated grain opacity (both Planck and Rosseland averaged) values
 * at various temperatures.  Data is in log cgs units. */
void read_grain_kappa_data(void)
{
  char setname[MAXLEN_PATH];

  hid_t file_id = my_H5Fopen(All.GrainKappaPath, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t dataset;

  sprintf(setname, "/NData");
  dataset = my_H5Dopen(file_id, setname);
  my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &IR_N_pts, setname);
  my_H5Dclose(dataset, setname);

  IR_logT = (double *) mymalloc("IR_logT", IR_N_pts * sizeof(double));
  sprintf(setname, "/LogTemperature");
  dataset = my_H5Dopen(file_id, setname);
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, IR_logT, setname);
  my_H5Dclose(dataset, setname);

  IR_logkappaP = (double *) mymalloc("IR_logkappaP", IR_N_pts * sizeof(double));
  sprintf(setname, "/LogKappaP");
  dataset = my_H5Dopen(file_id, setname);
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, IR_logkappaP, setname);
  my_H5Dclose(dataset, setname);

  IR_logkappaR = (double *) mymalloc("IR_logkappaR", IR_N_pts * sizeof(double));
  sprintf(setname, "/LogKappaR");
  dataset = my_H5Dopen(file_id, setname);
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, IR_logkappaR, setname);
  my_H5Dclose(dataset, setname);

  accIR_kappaP = gsl_interp_accel_alloc();
  accIR_kappaR = gsl_interp_accel_alloc();
  splineIR_kappaP = gsl_spline_alloc(gsl_interp_linear, IR_N_pts);
  splineIR_kappaR = gsl_spline_alloc(gsl_interp_linear, IR_N_pts);

  gsl_spline_init(splineIR_kappaP, IR_logT, IR_logkappaP, IR_N_pts);
  gsl_spline_init(splineIR_kappaR, IR_logT, IR_logkappaR, IR_N_pts);

  my_H5Fclose(file_id, All.GrainKappaPath);
}
#endif

void set_kappa_times_rho_IR(int i, struct sph_particle_data *kappaSphP) /*Expand this function later to get kappa*rho from the properties of the dust in the cell*/
{
#ifndef MRT_IR_GRAIN_KAPPA
  for(int num1=0; num1<IR_BINS; num1++)
    //  kappaSphP[i].KappaIR[num1] = 1e-12 * All.UnitLength_in_cm ;
       {
	 double meanweight = 4. / (3 * HYDROGEN_MASSFRAC + 1 + 4 * HYDROGEN_MASSFRAC * SphP[i].Ne);
	 //double meanweight = 2.33 ;
	 double temperature = kappaSphP[i].Utherm * All.UnitEnergy_in_cgs / All.UnitMass_in_g * GAMMA_MINUS1 * meanweight * PROTONMASS / BOLTZMANN;

	 double val = 0.0316 * (temperature/10.0) * (temperature/10.0) * All.UnitMass_in_g / pow(All.UnitLength_in_cm, 2.0) * kappaSphP[i].Density ;
	 kappaSphP[i].KappaIR_R[num1] = val ;
	 kappaSphP[i].KappaIR_P[num1] = 0.1*val/0.0316 ;

	 //    kappaSphP[i].KappaIR_R[num1] = kappaSphP[i].Density * 10.0 * All.UnitMass_in_g / All.UnitLength_in_cm / All.UnitLength_in_cm ;
	 //kappaSphP[i].KappaIR_R[num1] = 6.4815586e-17 * All.UnitLength_in_cm ;
	 //kappaSphP[i].KappaIR_P[num1] = kappaSphP[i].KappaIR_R[num1];
      //      terminate("kappa = %g \n", kappaSphP[i].KappaIR[num1]) ;
    }
    //kappaSphP[i].KappaIR[num1] = 6.4815586e-17 * All.UnitLength_in_cm ;

#else
  double DGR = 0.0;
  for(int chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
    {
      for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
        {
          DGR += SphP[i].MetalsDustFraction[chan][k];
        }
    }

  double meanweight = 4. / (3 * HYDROGEN_MASSFRAC + 1 + 4 * HYDROGEN_MASSFRAC * SphP[i].Ne);
  double temperature = SphP[i].Utherm * All.UnitEnergy_in_cgs / All.UnitMass_in_g * GAMMA_MINUS1 * meanweight * PROTONMASS / BOLTZMANN;

  double logT = log10(temperature);
  if(logT < IR_logT[0])
    {
      logT = IR_logT[0];
    }
  if(logT > IR_logT[IR_N_pts-1])
    {
      logT = IR_logT[IR_N_pts-1];
    }

  for(int num1 = 0; num1 < IR_BINS; num1++)
    {
      /* Tabulated grain opacities given in cm^2 per unit mass of dust.
       * Multiply by dust-to-gas ratio to get per unit mass of gas. */
      double logkappaP_cgs = gsl_spline_eval(splineIR_kappaP, logT, accIR_kappaP);
      double kappaP_cgs = pow(10.0, logkappaP_cgs);
      double kappaP_code = kappaP_cgs * All.UnitMass_in_g / (All.UnitLength_in_cm * All.UnitLength_in_cm);
      SphP[i].KappaIR_P[num1] = kappaP_code * DGR * SphP[i].Density;

      double logkappaR_cgs = gsl_spline_eval(splineIR_kappaR, logT, accIR_kappaR);
      double kappaR_cgs = pow(10.0, logkappaR_cgs);
      double kappaR_code = kappaR_cgs * All.UnitMass_in_g / (All.UnitLength_in_cm * All.UnitLength_in_cm);
      SphP[i].KappaIR_R[num1] = kappaR_code * DGR * SphP[i].Density;
    }
#endif
}
#endif




#ifdef MRT_RADIATION_PRESSURE
void do_radiation_pressure_source_terms(void)
{
  mpi_printf("RT: Adding the radiation pressure source term...\n") ;
  /*Single frequency d(rho v)/dt = nHI * sigma * F * E ({energy of absorbed photons}) / CLIGHT*/
  int idx, i, num1;
  double a3inv, hubble_a ;
  double c_speed =  2.99792458e10 / All.UnitVelocity_in_cm_per_s ;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
	continue;

      double dt = 0.5 * (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval ; /*factor of 0.5 included because of spilt for RK time integration*/

      if(All.ComovingIntegrationOn)
	{
	  a3inv = 1.0 / All.Time / All.Time / All.Time;
	  hubble_a = hubble_function(All.Time);
	}
      else
	{
	  a3inv = hubble_a = 1.0;
	}

      double KE_old = (SphP[i].Momentum[0]*SphP[i].Momentum[0] + SphP[i].Momentum[1]*SphP[i].Momentum[1] + SphP[i].Momentum[2]*SphP[i].Momentum[2])/P[i].Mass/2.0 ;


      int  j;

#ifndef MRT_NO_UV


      //       nH = (HYDROGEN_MASSFRAC * SphP[i].Density * All.cf_a3inv) / (PROTONMASS / All.UnitMass_in_g * All.HubbleParam); 

      double molecular_weight = 4 / (1 + 3 * HYDROGEN_MASSFRAC + 4 * HYDROGEN_MASSFRAC * SphP[i].Ne); 
       
      double nH = HYDROGEN_MASSFRAC * SphP[i].Density * a3inv / PROTONMASS * All.UnitMass_in_g / All.HubbleParam;
      double nHI = SphP[i].HI * nH;
      
#ifdef MRT_INCLUDE_HE
      double nHeI = SphP[i].HeI * nH;
      double nHeII = SphP[i].HeII * nH;
      double nHeIII = SphP[i].HeIII * nH;
#endif
      

      for(j=0;j<UV_BINS;j++)
	{
	  double mom_inj = 0.0 ;
	  double n_gamma = SphP[i].DensPhot[j] * 1e63 * SphP[i].Volume ;
	  double modF = sqrt(SphP[i].RT_F[j][0]*SphP[i].RT_F[j][0] + SphP[i].RT_F[j][1]*SphP[i].RT_F[j][1] + SphP[i].RT_F[j][2]*SphP[i].RT_F[j][2]) ;
	  if(modF==0.0)
	    modF = 1.0 ;

	  mom_inj += nHI * c_internal_units * mrt_sigma_HI[j] * P_HI[j] * n_gamma / c_speed  ; 
#ifdef MRT_INCLUDE_HE
	  mom_inj += nHeI * c_internal_units * mrt_sigma_HeI[j] * P_HeI[j] * n_gamma / c_speed ;
	  mom_inj += nHeII * c_internal_units * mrt_sigma_HeII[j] * P_HeII[j] * n_gamma / c_speed ;
#endif
	  if(isnan(mom_inj))
	    terminate("rad pressure gone wrong!!! nHI = %g HI = %g \n\n", nHI, SphP[i].HI) ;


	  for(num1=0;num1<3;num1++)                                                   
	    SphP[i].Momentum[num1] += dt * mom_inj * SphP[i].RT_F[j][num1]/modF ;
	  //		  SphP[i].Momentum[0] += dt * mom_inj * X ;
	  //	  SphP[i].Momentum[1] += dt * mom_inj * Y ;
	  //	  SphP[i].Momentum[2] += dt * mom_inj * Z ;
	}
#endif

#ifdef MRT_IR_LTE
#ifdef BOUNDARY_REFL_SOLIDSIDE_MINID
      if(P[i].ID < BOUNDARY_REFL_SOLIDSIDE_MINID || P[i].ID > BOUNDARY_REFL_SOLIDSIDE_MAXID)
	{
#endif
      for(j=UV_BINS;j<(UV_BINS+IR_BINS);j++)
	for(num1=0;num1<3;num1++)
	    SphP[i].Momentum[num1] += dt * SphP[i].KappaIR_R[j-UV_BINS] * SphP[i].RT_F[j][num1] * SphP[i].Volume / c_speed ; 
#ifdef BOUNDARY_REFL_SOLIDSIDE_MINID
	}
#endif
#endif


      double KE_new = (SphP[i].Momentum[0]*SphP[i].Momentum[0] + SphP[i].Momentum[1]*SphP[i].Momentum[1] + SphP[i].Momentum[2]*SphP[i].Momentum[2])/P[i].Mass/2.0 ;

      SphP[i].Energy += KE_new - KE_old ;
    }  
}
#endif


#if defined(MRT_CHEMISTRY_PS2011) || defined(MRT_CHEMISTRY_PS2009)

void mrt_get_sigma(void)
{
  double fac = 1.0 / All.UnitLength_in_cm / All.UnitLength_in_cm * All.HubbleParam * All.HubbleParam;

#ifndef MRT_MULTI_FREQUENCY
  mrt_sigma_HI[0] = 6.3e-18 * fac;
  P_HI[0] = 13.6 * ELECTRONVOLT_IN_ERGS / All.UnitEnergy_in_cgs * All.HubbleParam ;

#else
  int i, j, integral;
  double e, d_nu, e_start, e_end;
  double sum_HI_sigma, sum_HI_G;
  double hc, T_eff, I_nu;
  double sig, f, fac_two;
#ifdef MRT_INCLUDE_HE
  double sum_HeI_sigma, sum_HeII_sigma;
  double sum_HeI_G, sum_HeII_G;
#endif

  T_eff = 1e5 ; // K
  hc = CLIGHT * PLANCK;

  integral = 10000;

  fac_two = ELECTRONVOLT_IN_ERGS / All.UnitEnergy_in_cgs * All.HubbleParam;

  nu[0] = 13.6;
  nu[1] = 100.0 ;
  //  nu[1] = 24.6;
  //nu[2] = 54.4;
  //nu[3] = 100.0;
  /*
  sum_HI_sigma = 0.0;
  sum_HI_G = 0.0;
#ifdef MRT_INCLUDE_HE
  sum_HeI_G = sum_HeII_G = 0.0;
  sum_HeI_sigma = 0.0;
  sum_HeII_sigma = 0.0;
#endif
  */
  double lum_tot = 0.0;

  for(i = 0; i < UV_BINS; i++)
    {


      sum_HI_sigma = 0.0;
      sum_HI_G = 0.0;
#ifdef MRT_INCLUDE_HE
      sum_HeI_G = sum_HeII_G = 0.0;
      sum_HeI_sigma = 0.0;
      sum_HeII_sigma = 0.0;
#endif

      e_start = nu[i];

      //      if(i == MRT_BINS - 1)
      //  e_end = 500.0;
      // else
      e_end = nu[i + 1];

      d_nu = (e_end - e_start) / (float) (integral - 1);

      mrt_sigma_HI[i] = 0.0;
      G_HI[i] = 0.0;
      P_HI[i] = 0.0 ;
#ifdef MRT_INCLUDE_HE
      mrt_sigma_HeI[i] = 0.0;
      mrt_sigma_HeII[i] = 0.0;
      G_HeI[i] = G_HeII[i] = 0.0;
      P_HeI[i] = P_HeII[i] = 0.0;
#endif

      lum[i] = 0.0  ;
      for(j = 0; j < integral; j++)
        {
          e = e_start + j * d_nu;

          I_nu = 2.0 * pow(e * ELECTRONVOLT_IN_ERGS, 3) / (hc * hc) / (exp(e * ELECTRONVOLT_IN_ERGS / (BOLTZMANN * T_eff)) - 1.0); // Blackbody with T_eff

	  lum[i] += d_nu * I_nu / e;

          if(nu[i] >= 13.6)
            {
              f = sqrt((e / 13.6) - 1.0);

              if(j == 0)
                sig = 6.3e-18;
              else
                sig = 6.3e-18 * pow(13.6 / e, 4) * exp(4 - (4 * atan(f) / f)) / (1.0 - exp(-2 * M_PI / f));

              mrt_sigma_HI[i] += d_nu * sig * I_nu / e;

              sum_HI_sigma += d_nu * I_nu / e;

              G_HI[i] += d_nu * sig * (e - 13.6) * I_nu / e;

              P_HI[i] += d_nu * sig * e  * I_nu / e;

              sum_HI_G += d_nu * sig * I_nu / e;
            }

#ifdef MRT_INCLUDE_HE
          if(nu[i] >= 24.6)
            {
              f = sqrt((e / 24.6) - 1.0);

              if(j == 0)
                sig = 7.83e-18;
              else
                sig = 7.83e-18 * pow(24.6 / e, 4) * exp(4 - (4 * atan(f) / f)) / (1.0 - exp(-2 * M_PI / f));

              mrt_sigma_HeI[i] += d_nu * sig * I_nu / e;

              sum_HeI_sigma += d_nu * I_nu / e;

              G_HeI[i] += d_nu * sig * (e - 24.6) * I_nu / e;

              P_HeI[i] += d_nu * sig * e * I_nu / e;

              sum_HeI_G += d_nu * sig * I_nu / e;
            }

          if(nu[i] >= 54.4)
            {
              f = sqrt((e / 54.4) - 1.0);

              if(j == 0)
                sig = 1.58e-18;
              else
                sig = 1.58e-18 * pow(54.4 / e, 4) * exp(4 - (4 * atan(f) / f)) / (1.0 - exp(-2 * M_PI / f));

              mrt_sigma_HeII[i] += d_nu * sig * I_nu / e;

              sum_HeII_sigma += d_nu * I_nu / e;

              G_HeII[i] += d_nu * sig * (e - 54.4) * I_nu / e;

	      P_HeII[i] += d_nu * sig * e * I_nu / e;

              sum_HeII_G += d_nu * sig * I_nu / e;
            }
#endif
        }
     if(nu[i] >= 13.6)
       {
          mrt_sigma_HI[i] *= fac / sum_HI_sigma;
          G_HI[i] *= fac_two / sum_HI_G;
          P_HI[i] *= fac_two / sum_HI_G;
        }

#ifdef MRT_INCLUDE_HE
      if(nu[i] >= 24.6)
        {
          mrt_sigma_HeI[i] *= fac / sum_HeI_sigma;
          G_HeI[i] *= fac_two / sum_HeI_G;
          P_HeI[i] *= fac_two / sum_HeI_G;
        }

      if(nu[i] >= 54.4)
        {
          mrt_sigma_HeII[i] *= fac / sum_HeII_sigma;
          G_HeII[i] *= fac_two / sum_HeII_G;
          P_HeII[i] *= fac_two / sum_HeII_G;
        }
#endif

      lum_tot += lum[i] ;
    }

  for(i = 0; i < UV_BINS; i++)
    {
      lum[i] /= lum_tot ;
    }
      /*      if(nu[i] >= 13.6)
        {
          mrt_sigma_HI[i] *= fac / sum_HI_sigma;
          G_HI[i] *= fac_two / sum_HI_G;
        }

#ifdef MRT_INCLUDE_HE
      if(nu[i] >= 24.6)
        {
          mrt_sigma_HeI[i] *= fac / sum_HeI_sigma;
          G_HeI[i] *= fac_two / sum_HeI_G;
        }

      if(nu[i] >= 54.4)
        {
          mrt_sigma_HeII[i] *= fac / sum_HeII_sigma;
          G_HeII[i] *= fac_two / sum_HeII_G;
        }
#endif
}*/


  //  lum[2] *= 10.0 ;
  if(ThisTask == 0)
    for(i = 0; i < UV_BINS; i++)
      {
	printf("RT SIGMA: %g %g %g | %g %g %g | %g %g %g \n", mrt_sigma_HI[i] / fac, G_HI[i] / fac_two, P_HI[i] / fac_two, mrt_sigma_HeI[i] / fac, G_HeI[i] / fac_two, P_HeI[i] / fac_two, mrt_sigma_HeII[i] / fac, G_HeII[i] / fac_two, P_HeII[i] / fac_two);
	printf("Luminosity fractions i = %d \t %g \n", i, lum[i]) ;
      }

#endif
}

#endif

#endif /*MRT */
