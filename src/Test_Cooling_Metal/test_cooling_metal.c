/*
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/Test_Cooling_Metal/test_cooling_metal.c
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

#include <string.h>
#include <stdio.h>

#include "../allvars.h"
#include "../proto.h"


#ifdef TEST_COOLING_METAL

#if !defined (COOLING)
#error TEST_COOLING_COOLING requires GFM_COOLING_METAL and COOLING enabled
#endif

/* set parameters */
static double HydrogenNumberDensity = 0.001;    /* cm^-3 */
static double Tmin = 1e4;
static double Tmax = 1e9;
static double MetallicityInSolar = 1.0;
#ifdef GFM_AGN_RADIATION
double AgnRadFlux = 100;        /* erg s^-1 cm^-2 */
#endif
#ifdef RADCOOL
double Phios = 1e-20;
double Phins = 1e-20;
#ifdef RADCOOL_HOTHALO
//double PhiT6 = 3.16228e20 ;
//double PhiT7 = 3.46228e20 ;
//double PhiTho8 = 3.16228e19 ;

double PhiT6 = 3.1623e22;
double PhiT7 = 3.1623e22;
double PhiT8 = 3.1623e22;
#endif
#endif



void test_cooling_metal_init()
{
  double meanweight;

  mpi_printf("\nThis is Arepo, version %s.\n\nRunning on %d processors.\n\nCode was compiled with settings:\n\n", AREPO_VERSION, NTask);

  if(ThisTask == 0)
    output_compile_time_options();

  //All.UnitLength_in_cm = 3.085678e21;
  // All.UnitMass_in_g = 1.989e33;
  // All.UnitVelocity_in_cm_per_s = 1.0e5;
  All.UnitLength_in_cm = 1.;
  All.UnitEnergy_in_cgs = 1.;
  All.UnitMass_in_g = 1.;
  All.UnitVelocity_in_cm_per_s = 1.;
  All.GravityConstantInternal = 0;

  set_units();

  All.HubbleParam = 1.0;

  All.MaxMemSize = 700;

  All.TimeBegin = 0.99;
  All.Time = All.TimeBegin;
  All.ComovingIntegrationOn = 0;
#ifdef RADCOOL
  All.OldStarsOn = 0;
  All.NewStarsOn = 0;
  All.SelfShieldingOn = 0;
  All.SelfShieldingDensity = 0.01;
#ifdef RADCOOL_HOTHALO
  All.HotHaloOn = 1;
#endif
#endif
  All.CoolingOn = 1;
  All.StarformationOn = 1;

  All.MinMetalTemp = 1e2;

  All.MinEgySpec = 0;
  All.MinGasTemp = 100;

#ifdef GFM_AGN_RADIATION
  /* turn off self-shielding */
  All.SelfShieldingDensity = 1e20;
#endif

  meanweight = 4.0 / (1 + 3 * HYDROGEN_MASSFRAC);       /* note: assuming NEUTRAL GAS */

  All.MinEgySpec = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.MinGasTemp;
  All.MinEgySpec *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

#ifdef GFM_COOLING_METAL
  sprintf(All.CoolingTablePath, "%s", "./../Arepo_GFM_Tables/Cooling/cooling_metal_UVB.hdf5");
#ifdef GFM_AGN_RADIATION
  sprintf(All.CoolingTablePath, "%s", "./../Arepo_GFM_Tables/Cooling/cooling_metal_AGN_Compton.hdf5");
#endif
#if defined (RADCOOL) && !defined (RADCOOL_HOTHALO)
  sprintf(All.CoolingTablePath, "%s", "./../Arepo_GFM_Tables/Cooling/cooltable_rad.hdf5");
#endif
#if defined (RADCOOL) && defined (RADCOOL_HOTHALO)
  sprintf(All.CoolingTablePath, "%s", "./../grp3_lvl9/cooltable_rad_hothalo.hdf5");
#endif
  sprintf(All.YieldTablePath, "%s", "./../Arepo_GFM_Tables/Yields");
#endif
#if defined (RADCOOL) && !defined (RADCOOL_HOTHALO)
  sprintf(All.TreecoolFile, "%s", "./data/TREECOOL_HM");
#else
  sprintf(All.TreecoolFile, "%s", "./../grp3_lvl9/data/TREECOOL_fg_dec11");
#endif

#ifdef GFM_AGN_RADIATION
  sprintf(All.TreecoolFileAGN, "%s", "./data/TREECOOL_AGN");
#endif
#if defined (RADCOOL) && !defined (RADCOOL_HOTHALO)
  sprintf(All.TreecoolFileRAD, "%s", "./data/TREECOOL_RAD");
#else
#if defined (RADCOOL) && defined (RADCOOL_HOTHALO)
  sprintf(All.TreecoolFileRAD, "%s", "./../grp3_lvl9/data/TREECOOL_RAD_HH");
#endif
#endif

#ifdef GFM_DUST_COOLING
  sprintf(All.SelfShieldingFile, "%s", "./data/SelfShielding_Rahmati12");
#endif
  mymalloc_init();

#ifdef GFM_COOLING_METAL
  init_cooling_metal();
  mpi_printf("GFM_COOLING_METAL: Metal line cooling rates initialized.\n");
#endif

  //All.Time = All.TimeBegin;
  set_cosmo_factors_for_current_time();

#ifdef GFM_COOLING_METAL
  read_cooling_tables_current_time();
  mpi_printf("GFM_COOLING_METAL: Metal line cooling rates loaded.\n");
#endif

  InitCool();
}

void test_cooling_function()
{
  FILE *fp;
  char buf[200], msg[200];
  int i, j;
  int bin_num = 400;
  double ne_guess, u, u_min, u_max, delta_u, mu;
  double temp, LambdaNet, dnH, lnH;
  double nHmin = -8.0;
  double nHmax = 2.0;
  double bin_num_den = 20;
  double XHe_primordial, XHe_solar, XHe, XH, density, yhelium, metal_cooling;
#ifdef GFM_DUST_COOLING
  double dust_cooling;
#endif
  //dnH = (nHmax-nHmin)/bin_num_den ;
  //for(j=0; j<bin_num_den; j++){
  //lnH = ((float) j)*dnH + nHmin ;
  //HydrogenNumberDensity = pow(10.0, lnH) ;
  XHe_primordial = 1 - HYDROGEN_MASSFRAC;
  XHe_solar = 4.0 * GFM_SOLAR_HE_ABUNDANCE / (1.0 + 4.0 * GFM_SOLAR_HE_ABUNDANCE + GFM_SOLAR_METALLICITY);
  XHe = XHe_primordial + (XHe_solar - XHe_primordial) * MetallicityInSolar;     /* this is an approximation that XHe scales linearly with metallicity */
  XH = 1 - XHe - MetallicityInSolar * GFM_SOLAR_METALLICITY;
  density = PROTONMASS / XH * HydrogenNumberDensity;
  yhelium = (1.0 - XH - MetallicityInSolar * GFM_SOLAR_METALLICITY) / (4.0 * XH);

  test_cooling_metal_init();

  density /= (All.HubbleParam * All.HubbleParam);       //so that nh = 1 cm^{-3}
#if defined(GFM_AGN_RADIATION) || defined(GFM_UVB_CORRECTIONS) || defined(RADCOOL)
#ifdef GFM_AGN_RADIATION
  update_radiation_state(density, XH, AgnRadFlux);
  mpi_printf("CellsWithAGNRadiation = %d\n", CellsWithAGNRadiation);
#else
#ifdef RADCOOL
  update_radiation_state(density, XH, Phios, Phins
#ifdef RADCOOL_HOTHALO
                         , PhiT6, PhiT7, PhiT8
#endif
    );
#else
  update_radiation_state(density, XH, 0);
#endif
#endif
#endif

#ifdef GFM_COOLING_METAL
  update_gas_state(density, XH, MetallicityInSolar * GFM_SOLAR_METALLICITY);    //so that the gas is at solar metallicity
#endif
  density *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;

  mpi_printf("starting cooling function test\n");

  if(ThisTask == 0)
    {
      sprintf(buf, "%s%s", All.OutputDir, "cooling_curve.txt");

      if(!(fp = fopen(buf, "w")))
        {
          sprintf(msg, "error in opening file '%s'\n", buf);
          terminate(msg);
        }

      mu = (1 + 4 * yhelium) / (1 + yhelium);   //gas is neutral
      u_min = BOLTZMANN * Tmin / (GAMMA_MINUS1 * mu * PROTONMASS);      //gas at 10^1 K
      u_min *= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;

      mu = (1 + 4 * yhelium) / (2 + 3.0 * yhelium);     //total ionization
      u_max = BOLTZMANN * Tmax / (GAMMA_MINUS1 * mu * PROTONMASS);      //gas at 10^ K
      u_max *= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;

      delta_u = log10(u_max / u_min) / (bin_num - 1);

      for(i = 0; i < bin_num; i++)
        {
          ne_guess = 1.0;
          u = u_min * pow(10., i * delta_u);
          temp = convert_u_to_temp(u, density, &ne_guess);
          metal_cooling = 0.0;

#ifdef GFM_COOLING_METAL
#if !defined(GFM_AGN_RADIATION) && !defined(RADCOOL)
          metal_cooling = get_CoolingMetalRate(log10(MetallicityInSolar), log10(HydrogenNumberDensity), log10(temp));
#else
#ifdef RADCOOL
          metal_cooling = get_CoolingMetalRate(log10(Phios), log10(Phins)
#ifdef RADCOOL_HOTHALO
                                               , log10(PhiT6), log10(PhiT7), log10(PhiT8)
#endif
                                               , log10(MetallicityInSolar), log10(HydrogenNumberDensity), log10(temp));
#else
          metal_cooling = get_CoolingMetalRate(log10(AgnRadFlux), log10(MetallicityInSolar), log10(HydrogenNumberDensity), log10(temp));
#endif
#endif
#else
          metal_cooling = 0.0;
#endif
#ifdef GFM_DUST_COOLING
          dust_cooling = get_CoolingDustRate(log10(temp), 1.0e-4);
#endif

          LambdaNet = -CoolingRate(log10(temp), density, &ne_guess);
#ifdef GFM_DUST_COOLING
          fprintf(fp, "%e %e %e %e %e\n", log10(HydrogenNumberDensity), log10(temp), LambdaNet, metal_cooling, dust_cooling);
#else
          fprintf(fp, "%e %e %e %e\n", log10(HydrogenNumberDensity), log10(temp), LambdaNet, metal_cooling);
#endif
        }

      fclose(fp);
    }
  //  }
  mpi_printf("output written in file %s\n", buf);
  mpi_printf("cooling function test finished\n");
}

void test_cooling()
{
  FILE *fp;
  char buf[200], msg[200];
  int i;
  //int bin_num = 2;
  double ne_guess, u, delta_t, mu, time;
  double temp = 2e4;

  double XHe_primordial = 1 - HYDROGEN_MASSFRAC;
  double XHe_solar = 4.0 * GFM_SOLAR_HE_ABUNDANCE / (1.0 + 4.0 * GFM_SOLAR_HE_ABUNDANCE + GFM_SOLAR_METALLICITY);
  double XHe = XHe_primordial + (XHe_solar - XHe_primordial) * MetallicityInSolar;      /* this is an approximation that XHe scales linearly with metallicity */
  double XH = 1 - XHe - MetallicityInSolar * GFM_SOLAR_METALLICITY;
  double density = PROTONMASS / XH * HydrogenNumberDensity;
  double yhelium = (1.0 - XH - MetallicityInSolar * GFM_SOLAR_METALLICITY) / (4.0 * XH);
  double u_new, cool_time;
  double t_new;


  density /= (All.HubbleParam * All.HubbleParam);       //so that nh = 1 cm^{-3}
#if defined(GFM_AGN_RADIATION) || defined(GFM_UVB_CORRECTIONS) || defined(RADCOOL)
#ifdef GFM_AGN_RADIATION
  update_radiation_state(density, XH, AgnRadFlux);
  mpi_printf("CellsWithAGNRadiation = %d\n", CellsWithAGNRadiation);
#else
#ifdef RADCOOL
  update_radiation_state(density, XH, Phios, Phins
#ifdef RADCOOL_HOTHALO
                         , PhiT6, PhiT7, PhiT8
#endif
    );
#else
  update_radiation_state(density, XH, 0);
#endif
#endif
#endif
#ifdef GFM_COOLING_METAL
  update_gas_state(density, XH, MetallicityInSolar * GFM_SOLAR_METALLICITY);    //so that the gas is at solar metallicity
#endif
  density *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;

  mpi_printf("starting cooling test\n");

  if(ThisTask == 0)
    {
      sprintf(buf, "%s%s", All.OutputDir, "cooling_test.txt");

      if(!(fp = fopen(buf, "w")))
        {
          sprintf(msg, "error in opening file '%s'\n", buf);
          terminate(msg);
        }
      //temp = 1e7 ;
      find_abundances_and_rates(log10(temp), density, &ne_guess);
      mu = (1 + 4 * yhelium) / (1 + yhelium + ne_guess);
      u = BOLTZMANN * temp / (mu * PROTONMASS * GAMMA_MINUS1);
      cool_time = GetCoolingTime(u, density, &ne_guess);
      printf("Cooling Time = %le\n", cool_time / 3.1556926e+13);
      //delta_t = log10(5.0 / 1.0e-3) / bin_num;
      time = 0.0;
      int bin_num = 1000;
      //if(cool_time =0.0)
      //cool_time = 3.1556926e+13*100.0 ;
      cool_time = 1e16;
      delta_t = (2.0 / ((float) bin_num)) * cool_time;
      fprintf(fp, "%18.15e %18.15e\n", time / 3.1556926e+13, temp);
      double pressure = temp * density;
      printf("pressure = %g\t density = %g\n", pressure, density);
      for(i = 0; i <= bin_num; i++)
        {
          u_new = DoCooling(u, density, delta_t, &ne_guess);
          t_new = convert_u_to_temp(u_new, density, &ne_guess);
          //time *= pow(10., delta_t);
          time += delta_t;
          fprintf(fp, "%18.15e %18.15e\n", time / 3.1556926e+13, t_new);
          //printf("tnew = %g\n", t_new) ;
          //density = pressure/t_new ;
          u = u_new;
        }

      fclose(fp);
    }
  mpi_printf("output written in file %s\n", buf);
  mpi_printf("cooling test finished\n");
}

void test_cooling_metal(double u, double density, double dt, double *ne)
{
  double XHe_primordial = 1 - HYDROGEN_MASSFRAC;
  double XHe_solar = 4.0 * GFM_SOLAR_HE_ABUNDANCE / (1.0 + 4.0 * GFM_SOLAR_HE_ABUNDANCE + GFM_SOLAR_METALLICITY);
  double XHe = XHe_primordial + (XHe_solar - XHe_primordial) * MetallicityInSolar;      /* this is an approximation that XHe scales linearly with metallicity */
  double XH = 1 - XHe - MetallicityInSolar * GFM_SOLAR_METALLICITY;

  mpi_printf("starting metal cooling test\n");

  test_cooling_metal_init();
#if defined(GFM_AGN_RADIATION) || defined(GFM_UVB_CORRECTIONS) || defined(RADCOOL)
#ifdef GFM_AGN_RADIATION
  update_radiation_state(density, XH, AgnRadFlux);
  mpi_printf("CellsWithAGNRadiation = %d\n", CellsWithAGNRadiation);
#else
#ifdef RADCOOL
  update_radiation_state(density, XH, Phios, Phins
#ifdef RADCOOL_HOTHALO
                         , PhiT6, PhiT7, PhiT8
#endif
    );
#else
  update_radiation_state(density, XH, 0);
#endif
#endif
#endif
#ifdef GFM_COOLING_METAL
  update_gas_state(density, XH, MetallicityInSolar * GFM_SOLAR_METALLICITY);    /*so that the gas is at solar metallicity */
#endif

  if(ThisTask == 0)
    {
      DoCooling(u, density, dt, ne);
    }

  mpi_printf("metal cooling test finished\n");
}

#endif
