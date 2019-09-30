/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/fm_star_formation/test_sfr.c
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

#ifdef TEST_SFR


#ifndef FM_SFR

#if !defined (COOLING) && !defined (USE_SFR)
#error TEST_SFR requires COOLING and USE_SFR enabled
#endif

void test_sfr_eEOS(void)
{
  All.UnitLength_in_cm = 3.085678e24;
  All.UnitMass_in_g = 1.989e43;
  All.UnitVelocity_in_cm_per_s = 1.0e5;
  All.GravityConstantInternal = 0;

  All.HubbleParam = 0.73;
  All.OmegaBaryon = 0.04;

  All.MaxMemSize = 700;

  All.TimeBegin = 0.0;
  All.ComovingIntegrationOn = 0;

  All.CoolingOn = 1;
  All.StarformationOn = 1;

  All.MinEgySpec = 0;
  All.MinGasTemp = 100;

  All.CritPhysDensity = 0;
  All.TemperatureThresh = 0.0;
  All.MaxSfrTimescale = 0.00227;
  All.CritOverDensity = 57.7;
  All.TempSupernova = 5.73e7;
  All.TempClouds = 1000.0;
  All.FactorSN = 0.264089;
  All.FactorEVP = 573.0;

#ifdef MODIFIED_EOS
  All.FactorDensThresh = 0.5;
  All.FactorUthermAtThresh = 3.0;
  All.FactorUthermJoin = 6.0;

  check_modified_eos_parameters();
#endif

  set_units();

  double meanweight = 4.0 / (1 + 3 * HYDROGEN_MASSFRAC);        /* note: assuming NEUTRAL GAS */

  All.MinEgySpec = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.MinGasTemp;
  All.MinEgySpec *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

#ifdef GFM_COOLING_METAL
  sprintf(All.CoolingTablePath, "%s", "./data/Arepo_GFM_Tables/Cooling");
#ifdef GFM_AGN_RADIATION
  sprintf(All.CoolingTablePath, "%s", "./data/Arepo_GFM_Tables/Cooling/cooling_metal_AGN_Compton_self_shielding_Rahmati12.hdf5");
#endif
  sprintf(All.YieldTablePath, "%s", "./data/Arepo_GFM_Tables/Yields");
#endif
  sprintf(All.TreecoolFile, "%s", "./data/TREECOOL_fg_dec11");
#ifdef GFM_AGN_RADIATION
  sprintf(All.TreecoolFileAGN, "%s", "./data/TREECOOL_AGN");
#endif
  mymalloc_init();

#ifdef GFM_COOLING_METAL
  init_cooling_metal();
  mpi_printf("GFM_COOLING_METAL: Metal line cooling rates initialized.\n");
#endif

#ifdef COOLING
  All.Time = All.TimeBegin;
  set_cosmo_factors_for_current_time();
  InitCool();
#endif

  init_clouds();

}
#endif

#ifdef FM_SFR

#if !defined (COOLING) && !defined (FM_SFR)
#error TEST_SFR requires COOLING and FM_SFR enabled
#endif

void test_fm_sfr(void)
{
  All.UnitLength_in_cm = 3.085678e21;
  All.UnitMass_in_g = 1.989e43;
  All.UnitVelocity_in_cm_per_s = 1.0e5;
  All.GravityConstantInternal = 0;

  All.HubbleParam = 0.7;

  set_units();

  All.MaxMemSize = 700;

  All.TimeBegin = 0.0;
  All.ComovingIntegrationOn = 0;

  All.CoolingOn = 1;
  All.StarformationOn = 1;

  All.MinEgySpec = 0;
  All.MinGasTemp = 100;

  All.DensThreshold = 0.0;
  All.SfrEfficiency = 0.02;

#ifdef USE_POLYTROPIC_EQSTATE
  All.UthermThreshold = 10000.0;
#endif

  All.SofteningGas = 0.647;
  All.TargetGasMassFactor = 1.0;
  All.ReferenceGasPartMass = 7.49999e-06;

  mymalloc_init();

  All.Time = All.TimeBegin;
  set_cosmo_factors_for_current_time();
  init_star_formation();
}
#endif

void test_sfr(void)
{
  mpi_printf("\nThis is Arepo, version %s.\n\nRunning on %d processors.\n\nCode was compiled with settings:\n\n", AREPO_VERSION, NTask);

  if(ThisTask == 0)
    output_compile_time_options();

#ifndef FM_SFR
  test_sfr_eEOS();
#endif

#ifdef FM_SFR
  test_fm_sfr();
#endif
}

#endif
