/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/grackle/grackle.c
 * \date        01/2016
 * \author      Matthew Smith
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

#ifdef COOLING
#ifdef GRACKLE

#ifdef GRACKLE_TAB
#ifdef GRACKLE_D
#error "GRACKLE: Cannot use tabulated mode with deuterium chemistry"
#elif defined(GRACKLE_H2)
#error "GRACKLE: Cannot use tabulated mode with H2 chemistry"
#endif
#endif

void initialise_grackle()
{
#ifdef GRACKLE_VERBOSE
  grackle_verbose = 1;
#endif

  my_grackle_units.comoving_coordinates = All.ComovingIntegrationOn;
  my_grackle_units.density_units = All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam * All.cf_a3inv;
  my_grackle_units.length_units = All.UnitLength_in_cm * All.HubbleParam * All.cf_atime;
  my_grackle_units.time_units = All.UnitTime_in_s / All.HubbleParam;
  my_grackle_units.velocity_units = All.UnitVelocity_in_cm_per_s;
  my_grackle_units.a_units = 1.0;


  if(set_default_chemistry_parameters() == 0)
    terminate("Grackle: Error in set_default_chemistry_parameters.\n");

  /* Set parameter values for chemistry. */
  grackle_data.use_grackle = All.GrackleOn;     /* chemistry on  */
  grackle_data.with_radiative_cooling = All.GrackleRadiativeCooling;    /* cooling on */
#ifdef GRACKLE_TAB
  grackle_data.primordial_chemistry = 0;        /* Tabulated mode */
#elif defined(GRACKLE_D)
  grackle_data.primordial_chemistry = 3;        /* H, He, H2, D */
#elif defined(GRACKLE_H2)
  grackle_data.primordial_chemistry = 2;        /* H, He, H2 */
#else
  grackle_data.primordial_chemistry = 1;        /* H, He */
#endif
  grackle_data.metal_cooling = All.GrackleMetalCooling; /* metal cooling on */
  grackle_data.UVbackground = All.GrackleUVB;   /* UV background on */
  grackle_data.grackle_data_file = All.GrackleDataFile; /* data file */

  if(grackle_data.Gamma != GAMMA)
    terminate("The value of Gamma in AREPO and GRACKLE must be the same (almost certainly 5/3)");

  if(initialize_chemistry_data(&my_grackle_units, All.cf_atime) == 0)
    terminate("Grackle: Error in initialize_chemistry_data.\n");

}

void cooling_only()             /* normal cooling routine when star formation is disabled */
{
  int nactive = 0;

  CPU_Step[CPU_MISC] += measure_time();
  mpi_printf("GRACKLE: Cooling active cells\n");
  cool_active_cells();
  CPU_Step[CPU_COOLINGSFR] += measure_time();
}

void cool_active_cells()
{
  gr_float *density, *energy, *dummy_velocity, *metal_density;
#ifndef GRACKLE_TAB
  gr_float *HI_density, *HII_density, *HeI_density, *HeII_density, *HeIII_density;
  gr_float *HM_density, *H2I_density, *H2II_density;
  gr_float *DI_density, *DII_density, *HDI_density;
  gr_float *e_density;
  gr_float *pressure;
#endif
  double dt, dtime;
  int grid_rank = 3;
  int grid_dimension[3], grid_start[3], grid_end[3];
  int ngas_active = 0;

  set_cosmo_factors_for_current_time();

  my_grackle_units.density_units = All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam * All.cf_a3inv;
  my_grackle_units.length_units = All.UnitLength_in_cm * All.HubbleParam * All.cf_atime;

  density = mymalloc("grackle_density", TimeBinsHydro.NActiveParticles * sizeof(gr_float));
  energy = mymalloc("grackle_energy", TimeBinsHydro.NActiveParticles * sizeof(gr_float));
  metal_density = mymalloc("grackle_metal_density", TimeBinsHydro.NActiveParticles * sizeof(gr_float));
  dummy_velocity = NULL;
#ifndef GRACKLE_TAB
  HI_density = mymalloc("grackle_HI_density", TimeBinsHydro.NActiveParticles * sizeof(gr_float));
  HII_density = mymalloc("grackle_HII_density", TimeBinsHydro.NActiveParticles * sizeof(gr_float));
  HeI_density = mymalloc("grackle_HeI_density", TimeBinsHydro.NActiveParticles * sizeof(gr_float));
  HeII_density = mymalloc("grackle_HeII_density", TimeBinsHydro.NActiveParticles * sizeof(gr_float));
  HeIII_density = mymalloc("grackle_HeIII_density", TimeBinsHydro.NActiveParticles * sizeof(gr_float));
  e_density = mymalloc("grackle_e_density", TimeBinsHydro.NActiveParticles * sizeof(gr_float));
#ifdef GRACKLE_H2
  HM_density = mymalloc("grackle_HM_density", TimeBinsHydro.NActiveParticles * sizeof(gr_float));
  H2I_density = mymalloc("grackle_H2I_density", TimeBinsHydro.NActiveParticles * sizeof(gr_float));
  H2II_density = mymalloc("grackle_H2II_density", TimeBinsHydro.NActiveParticles * sizeof(gr_float));
#else
  HM_density = NULL;
  H2I_density = NULL;
  H2II_density = NULL;
#endif
#ifdef GRACKLE_D
  DI_density = mymalloc("grackle_DI_density", TimeBinsHydro.NActiveParticles * sizeof(gr_float));
  DII_density = mymalloc("grackle_DII_density", TimeBinsHydro.NActiveParticles * sizeof(gr_float));
  HDI_density = mymalloc("grackle_HDI_density", TimeBinsHydro.NActiveParticles * sizeof(gr_float));
#else
  DI_density = NULL;
  DII_density = NULL;
  HDI_density = NULL;
#endif
  pressure = mymalloc("grackle_pressure", TimeBinsHydro.NActiveParticles * sizeof(gr_float));
#endif
  for(int timebin = All.LowestActiveTimeBin; timebin <= All.HighestActiveTimeBin; timebin++)
  {
    dt = (((integertime) 1) << timebin) * All.Timebase_interval;
    dtime = All.cf_atime * dt / All.cf_time_hubble_a;
    mpi_printf("GRACKLE: dtime = %g\n",dtime);
    for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;
      if(P[i].Mass == 0 && P[i].ID == 0)
        continue;           /* skip cells that have been swallowed or eliminated */
      if(P[i].TimeBinHydro != timebin)
        continue;

      density[ngas_active] = SphP[i].Density;
      energy[ngas_active] = SphP[i].Utherm * All.cf_atime * All.cf_atime; /* Grackle expects comoving */

#ifdef METALS
      metal_density[ngas_active] = SphP[i].Metallicity * density[ngas_active];
#else
      metal_density[ngas_active] = All.GrackleInitialMetallicity * density[ngas_active];
#endif

#ifndef GRACKLE_TAB
      HI_density[ngas_active] = SphP[i].GrackleSpeciesFraction[0] * density[ngas_active];
      HII_density[ngas_active] = SphP[i].GrackleSpeciesFraction[1] * density[ngas_active];
      HeI_density[ngas_active] = SphP[i].GrackleSpeciesFraction[2] * density[ngas_active];
      HeII_density[ngas_active] = SphP[i].GrackleSpeciesFraction[3] * density[ngas_active];
      HeIII_density[ngas_active] = SphP[i].GrackleSpeciesFraction[4] * density[ngas_active];
      e_density[ngas_active] = SphP[i].e_frac * density[ngas_active] * (PROTONMASS / ELECTRONMASS);
#ifdef GRACKLE_H2
      HM_density[ngas_active] = SphP[i].GrackleSpeciesFraction[5] * density[ngas_active];
      H2I_density[ngas_active] = SphP[i].GrackleSpeciesFraction[6] * density[ngas_active];
      H2II_density[ngas_active] = SphP[i].GrackleSpeciesFraction[7] * density[ngas_active];
#endif
#ifdef GRACKLE_D
      DI_density[ngas_active] = SphP[i].GrackleSpeciesFraction[8] * density[ngas_active];
      DII_density[ngas_active] = SphP[i].GrackleSpeciesFraction[9] * density[ngas_active];
      HDI_density[ngas_active] = SphP[i].GrackleSpeciesFraction[10] * density[ngas_active];
#endif
#endif
      ngas_active++;
    }

  for(int j = 1; j < 3; j++)
    {
      grid_dimension[j] = 1;
      grid_start[j] = 0;
      grid_end[j] = 0;
    }
  grid_dimension[0] = ngas_active;
  grid_end[0] = ngas_active - 1;
  mpi_printf("GRACKLE: ngas_active = %d\n",ngas_active);
  for(int j = 0; j < ngas_active; j++)
    mpi_printf("rho: %g u: %g HI: %g HII %g HeI %g HeII %g HeIII %g e: %g\n",density[j],energy[j],HI_density[j],HII_density[j],HeI_density[j],HeII_density[j],HeIII_density[j],e_density[j]);

#ifdef GRACKLE_TAB
  if(solve_chemistry_table(&my_grackle_units, All.cf_atime, dtime, grid_rank, grid_dimension, 
                          grid_start, grid_end, density, energy, dummy_velocity, dummy_velocity, dummy_velocity, 
                          metal_density) == 0)
    terminate("GRACKLE: Error in solve_chemistry_table.\n");

  /* Pressure in tabular case calculated later */
#else
  if(solve_chemistry(&my_grackle_units, All.cf_atime, dtime, grid_rank, grid_dimension, grid_start, 
                    grid_end, density, energy, dummy_velocity, dummy_velocity, dummy_velocity, HI_density, 
                    HII_density, HM_density, HeI_density, HeII_density, HeIII_density, H2I_density, 
                    H2II_density, DI_density, DII_density, HDI_density, e_density, metal_density) == 0)
    terminate("GRACKLE: Error in solve_chemistry.\n");

  if(calculate_pressure(&my_grackle_units, All.cf_atime, grid_rank, grid_dimension, grid_start, grid_end,
                        density, energy, HI_density, HII_density, HM_density, HeI_density, HeII_density, 
                        HeIII_density, H2I_density, H2II_density, DI_density, DII_density, HDI_density, 
                        e_density, metal_density, pressure) == 0)
    terminate("GRACKLE: Error in calculate_pressure");

#endif

  ngas_active = 0;
  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;
      if(P[i].Mass == 0 && P[i].ID == 0)
        continue;           /* skip cells that have been swallowed or eliminated */
      if(P[i].TimeBinHydro != timebin)
        continue;

      SphP[i].Utherm = energy[ngas_active] / (All.cf_atime * All.cf_atime);
#ifndef GRACKLE_TAB
      SphP[i].GrackleSpeciesFraction[0] = HI_density[ngas_active] / density[ngas_active];
      SphP[i].GrackleSpeciesFraction[1] = HII_density[ngas_active] / density[ngas_active];
      SphP[i].GrackleSpeciesFraction[2] = HeI_density[ngas_active] / density[ngas_active];
      SphP[i].GrackleSpeciesFraction[3] = HeII_density[ngas_active] / density[ngas_active];
      SphP[i].GrackleSpeciesFraction[4] = HeIII_density[ngas_active] / density[ngas_active];
      SphP[i].e_frac = e_density[ngas_active] / density[ngas_active] / (PROTONMASS / ELECTRONMASS);
#ifdef GRACKLE_H2
      SphP[i].GrackleSpeciesFraction[5] = HM_density[ngas_active] / density[ngas_active];
      SphP[i].GrackleSpeciesFraction[6] = H2I_density[ngas_active] / density[ngas_active];
      SphP[i].GrackleSpeciesFraction[7] = H2II_density[ngas_active] / density[ngas_active];
#endif
#ifdef GRACKLE_D
      SphP[i].GrackleSpeciesFraction[8] = DI_density[ngas_active] / density[ngas_active];
      SphP[i].GrackleSpeciesFraction[9] = DII_density[ngas_active] / density[ngas_active];
      SphP[i].GrackleSpeciesFraction[10] = HDI_density[ngas_active] / density[ngas_active];
#endif
      SphP[i].Pressure = pressure[ngas_active];

      for(int j = 0; j < GRACKLE_SPECIES_NUMBER; j++)
        SphP[i].GrackleSpeciesMass[j] = P[i].Mass * SphP[i].GrackleSpeciesFraction[j];

      SphP[i].e_mass = P[i].Mass * SphP[i].e_frac;
#else
      SphP[i].Pressure = GAMMA_MINUS1 * SphP[i].Density * SphP[i].Utherm; /* For tabular */
#endif
      ngas_active++;
        
    }
  }


#ifndef GRACKLE_TAB
#ifdef GRACKLE_D
  myfree(pressure);
  myfree(HDI_density);
  myfree(DII_density);
  myfree(DI_density);
#endif
#ifdef GRACKLE_H2
  myfree(H2II_density);
  myfree(H2I_density);
  myfree(HM_density);
#endif
  myfree(e_density);
  myfree(HeIII_density);
  myfree(HeII_density);
  myfree(HeI_density);
  myfree(HII_density);
  myfree(HI_density);
#endif
  myfree(metal_density);
  myfree(energy);
  myfree(density);
}

double get_temp_individual_cell_grackle(int i)
{
  gr_float density, energy, metal_density;
  gr_float temperature;
#ifndef GRACKLE_TAB
  gr_float HI_density, HII_density, HeI_density, HeII_density, HeIII_density;
  gr_float HM_density, H2I_density, H2II_density;
  gr_float DI_density, DII_density, HDI_density;
  gr_float e_density;
#endif
  int grid_rank = 3;
  int grid_dimension[3], grid_start[3], grid_end[3];

  set_cosmo_factors_for_current_time();

  density = SphP[i].Density;
  energy = SphP[i].Utherm * All.cf_atime * All.cf_atime;
#ifdef METALS
  metal_density = SphP[i].Metallicity * density;
#else
  metal_density = All.GrackleInitialMetallicity * density;
#endif

#ifndef GRACKLE_TAB
  HI_density = SphP[i].GrackleSpeciesFraction[0] * density;
  HII_density = SphP[i].GrackleSpeciesFraction[1] * density;
  HeI_density = SphP[i].GrackleSpeciesFraction[2] * density;
  HeII_density = SphP[i].GrackleSpeciesFraction[3] * density;
  HeIII_density = SphP[i].GrackleSpeciesFraction[4] * density;
  e_density = SphP[i].e_frac * density * (PROTONMASS / ELECTRONMASS);
#ifdef GRACKLE_H2
  HM_density = SphP[i].GrackleSpeciesFraction[5] * density;
  H2I_density = SphP[i].GrackleSpeciesFraction[6] * density;
  H2II_density = SphP[i].GrackleSpeciesFraction[7] * density;
#else
  HM_density = 0;
  H2I_density = 0;
  H2II_density = 0;
#endif
#ifdef GRACKLE_D
  DI_density = SphP[i].GrackleSpeciesFraction[8] * density;
  DII_density = SphP[i].GrackleSpeciesFraction[9] * density;
  HDI_density = SphP[i].GrackleSpeciesFraction[10] * density;
#else
  DI_density = 0;
  DII_density = 0;
  HDI_density = 0;
#endif
#endif

  for(int j = 1; j < 3; j++)
    {
      grid_dimension[j] = 1;
      grid_start[j] = 0;
      grid_end[j] = 0;
    }
  grid_dimension[0] = 1;
  grid_start[0] = 0;
  grid_end[0] = 0;

#ifdef GRACKLE_TAB
  if(calculate_temperature_table(&my_grackle_units, All.cf_atime, grid_rank, grid_dimension, grid_start, 
                                grid_end, &density, &energy, &metal_density, &temperature) == 0)
    terminate("GRACKLE: Error in calculate_temperature_table.\n");
#else
  if(calculate_temperature(&my_grackle_units, All.cf_atime, grid_rank, grid_dimension, grid_start, grid_end,
                           &density, &energy, &HI_density, &HII_density, &HM_density, &HeI_density, 
                           &HeII_density, &HeIII_density, &H2I_density, &H2II_density, &DI_density, 
                           &DII_density, &HDI_density, &e_density, &metal_density, &temperature) == 0)
    terminate("GRACKLE: Error in calculate_temperature.\n");
#endif
return temperature;
}

double get_cooling_time_individual_cell_grackle(int i)
{
  gr_float density, energy, dummy_velocity, metal_density;
  gr_float t_cool;
#ifndef GRACKLE_TAB
  gr_float HI_density, HII_density, HeI_density, HeII_density, HeIII_density;
  gr_float HM_density, H2I_density, H2II_density;
  gr_float DI_density, DII_density, HDI_density;
  gr_float e_density;
#endif
  int grid_rank = 3;
  int grid_dimension[3], grid_start[3], grid_end[3];

  set_cosmo_factors_for_current_time();

  density = SphP[i].Density;
  energy = SphP[i].Utherm * All.cf_atime * All.cf_atime;
#ifdef METALS
  metal_density = SphP[i].Metallicity * density;
#else
  metal_density = All.GrackleInitialMetallicity * density;
#endif
  dummy_velocity = 0;

#ifndef GRACKLE_TAB
  HI_density = SphP[i].GrackleSpeciesFraction[0] * density;
  HII_density = SphP[i].GrackleSpeciesFraction[1] * density;
  HeI_density = SphP[i].GrackleSpeciesFraction[2] * density;
  HeII_density = SphP[i].GrackleSpeciesFraction[3] * density;
  HeIII_density = SphP[i].GrackleSpeciesFraction[4] * density;
  e_density = SphP[i].e_frac * density * (PROTONMASS / ELECTRONMASS);
#ifdef GRACKLE_H2
  HM_density = SphP[i].GrackleSpeciesFraction[5] * density;
  H2I_density = SphP[i].GrackleSpeciesFraction[6] * density;
  H2II_density = SphP[i].GrackleSpeciesFraction[7] * density;
#else
  HM_density = 0;
  H2I_density = 0;
  H2II_density = 0;
#endif
#ifdef GRACKLE_D
  DI_density = SphP[i].GrackleSpeciesFraction[8] * density;
  DII_density = SphP[i].GrackleSpeciesFraction[9] * density;
  HDI_density = SphP[i].GrackleSpeciesFraction[10] * density;
#else
  DI_density = 0;
  DII_density = 0;
  HDI_density = 0;
#endif
#endif

  for(int j = 1; j < 3; j++)
    {
      grid_dimension[j] = 1;
      grid_start[j] = 0;
      grid_end[j] = 0;
    }
  grid_dimension[0] = 1;
  grid_start[0] = 0;
  grid_end[0] = 0;

#ifdef GRACKLE_TAB
  if(calculate_cooling_time_table(&my_grackle_units, All.cf_atime, grid_rank, grid_dimension, grid_start, 
                                grid_end, &density, &energy, &dummy_velocity, &dummy_velocity, &dummy_velocity, &metal_density, &t_cool) == 0)
    terminate("GRACKLE: Error in calculate_cooling_time_table.\n");
#else
  if(calculate_cooling_time(&my_grackle_units, All.cf_atime, grid_rank, grid_dimension, grid_start, grid_end,
                           &density, &energy, &dummy_velocity, &dummy_velocity, &dummy_velocity, &HI_density, &HII_density, &HM_density, &HeI_density, 
                           &HeII_density, &HeIII_density, &H2I_density, &H2II_density, &DI_density, 
                           &DII_density, &HDI_density, &e_density, &metal_density, &t_cool) == 0)
    terminate("GRACKLE: Error in calculate_cooling_time.\n");
#endif
return t_cool;
}

#ifndef GRACKLE_TAB
#ifndef GRACKLE_ABUNDANCE_IN_ICS
void grackle_initialise_abundances()
{
  MyFloat tinynumber = 1.0e-20;
  MyFloat metallicity;

  mpi_printf("GRACKLE: Initialising abundances\n");
#ifndef METALS
  metallicity = All.GrackleInitialMetallicity;
#endif
  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

#ifdef METALS
      metallicity = SphP[i].Metallicity;
#endif
  		SphP[i].GrackleSpeciesFraction[0] = HYDROGEN_MASSFRAC * (1.0 - metallicity);    /* HI    */
      SphP[i].GrackleSpeciesFraction[1] = tinynumber;   /* HII   */
      SphP[i].GrackleSpeciesFraction[2] = (1.0 - HYDROGEN_MASSFRAC) * (1.0 - metallicity);        /* HeI   */
      SphP[i].GrackleSpeciesFraction[3] = tinynumber;   /* HeII  */
      SphP[i].GrackleSpeciesFraction[4] = tinynumber;   /* HeIII */
#ifdef GRACKLE_H2
      SphP[i].GrackleSpeciesFraction[5] = tinynumber;   /* HM    */
      SphP[i].GrackleSpeciesFraction[6] = tinynumber;   /* H2I   */
      SphP[i].GrackleSpeciesFraction[7] = tinynumber;   /* H2II  */
#endif
#ifdef GRACKLE_D
      SphP[i].GrackleSpeciesFraction[8] = tinynumber;   /* DI    */
      SphP[i].GrackleSpeciesFraction[9] = tinynumber;   /* DII   */
      SphP[i].GrackleSpeciesFraction[10] = tinynumber;  /* HDI   */
#endif
      SphP[i].e_frac = tinynumber;      /* electrons */
    }
}
#endif

/* Used on startup. Chemistry solved iteratively (without internal energy being modified) until abundances converge */
void grackle_converge_abundances()
{
  gr_float *density, *energy, *x_velocity, *y_velocity, *z_velocity, *metal_density;
  gr_float *HI_density, *HII_density, *HeI_density, *HeII_density, *HeIII_density;
  gr_float *HM_density, *H2I_density, *H2II_density;
  gr_float *DI_density, *DII_density, *HDI_density;
  gr_float *e_density;
  gr_float new_species_frac[GRACKLE_SPECIES_NUMBER];
  gr_float new_e_frac;
  int *con_list;
  double dtime;
  int iter, con_local;
  int grid_rank = 3;
  int grid_dimension[3], grid_start[3], grid_end[3];
  int ngas_con, ngas_uncon;

  mpi_printf("GRACKLE: Converging abundances\n");
  set_cosmo_factors_for_current_time();

  dtime = 0.01 * SEC_PER_MEGAYEAR / All.UnitTime_in_s;

  my_grackle_units.density_units = All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam * All.cf_a3inv;
  my_grackle_units.length_units = All.UnitLength_in_cm * All.HubbleParam * All.cf_atime;

  /* Allocate buffers for abundances etc. */
  density = mymalloc("grackle_density", NumGas * sizeof(gr_float));
  energy = mymalloc("grackle_energy", NumGas * sizeof(gr_float));
  metal_density = mymalloc("grackle_metal_density", NumGas * sizeof(gr_float));
  x_velocity = NULL;
  y_velocity = NULL;
  z_velocity = NULL;
  HI_density = mymalloc("grackle_HI_density", NumGas * sizeof(gr_float));
  HII_density = mymalloc("grackle_HII_density", NumGas * sizeof(gr_float));
  HeI_density = mymalloc("grackle_HeI_density", NumGas * sizeof(gr_float));
  HeII_density = mymalloc("grackle_HeII_density", NumGas * sizeof(gr_float));
  HeIII_density = mymalloc("grackle_HeIII_density", NumGas * sizeof(gr_float));
  e_density = mymalloc("grackle_e_density", NumGas * sizeof(gr_float));
#ifdef GRACKLE_H2
  HM_density = mymalloc("grackle_HM_density", NumGas * sizeof(gr_float));
  H2I_density = mymalloc("grackle_H2I_density", NumGas * sizeof(gr_float));
  H2II_density = mymalloc("grackle_H2II_density", NumGas * sizeof(gr_float));
#else
  HM_density = NULL;
  H2I_density = NULL;
  H2II_density = NULL;
#endif
#ifdef GRACKLE_D
  DI_density = mymalloc("grackle_DI_density", NumGas * sizeof(gr_float));
  DII_density = mymalloc("grackle_DII_density", NumGas * sizeof(gr_float));
  HDI_density = mymalloc("grackle_HDI_density", NumGas * sizeof(gr_float));
#else
  DI_density = NULL;
  DII_density = NULL;
  HDI_density = NULL;
#endif
  con_list = mymalloc("grackle_converged_list", NumGas * sizeof(int));

  ngas_con = 0;
  iter = 0;
  while(ngas_con < NumGas)
    {
      ngas_uncon = 0;
      for(int i = 0; i < NumGas; i++)
        {
          if(P[i].Type != 0)
            continue;
          if(iter = 0)
            con_list[i] = 0;    /* On first pass assume all abundances are unconverged */
          else if(con_list[i] == 1)
            continue;           /* This cell already has converged abundances. Don't bother passing to GRACKLE again */

          density[ngas_uncon] = SphP[i].Density;
          energy[ngas_uncon] = SphP[i].Utherm * All.cf_atime * All.cf_atime;

#ifdef METALS
          metal_density[ngas_uncon] = SphP[i].Metallicity * density[ngas_uncon];
#else
          metal_density[ngas_uncon] = All.GrackleInitialMetallicity * density[ngas_uncon];
#endif


          HI_density[ngas_uncon] = SphP[i].GrackleSpeciesFraction[0] * density[ngas_uncon];
          HII_density[ngas_uncon] = SphP[i].GrackleSpeciesFraction[1] * density[ngas_uncon];
          HeI_density[ngas_uncon] = SphP[i].GrackleSpeciesFraction[2] * density[ngas_uncon];
          HeII_density[ngas_uncon] = SphP[i].GrackleSpeciesFraction[3] * density[ngas_uncon];
          HeIII_density[ngas_uncon] = SphP[i].GrackleSpeciesFraction[4] * density[ngas_uncon];
          e_density[ngas_uncon] = SphP[i].e_frac * density[ngas_uncon] * (PROTONMASS / ELECTRONMASS);
#ifdef GRACKLE_H2
          HM_density[ngas_uncon] = SphP[i].GrackleSpeciesFraction[5] * density[ngas_uncon];
          H2I_density[ngas_uncon] = SphP[i].GrackleSpeciesFraction[6] * density[ngas_uncon];
          H2II_density[ngas_uncon] = SphP[i].GrackleSpeciesFraction[7] * density[ngas_uncon];
#endif
#ifdef GRACKLE_D
          DI_density[ngas_uncon] = SphP[i].GrackleSpeciesFraction[8] * density[ngas_uncon];
          DII_density[ngas_uncon] = SphP[i].GrackleSpeciesFraction[9] * density[ngas_uncon];
          HDI_density[ngas_uncon] = SphP[i].GrackleSpeciesFraction[10] * density[ngas_uncon];
#endif
          ngas_uncon++;
        }

      for(int j = 1; j < 3; j++)
        {
          grid_dimension[j] = 1;
          grid_start[j] = 0;
          grid_end[j] = 0;
        }
      grid_dimension[0] = ngas_uncon;
      grid_start[0] = 0;
      grid_end[0] = ngas_uncon - 1;

      mpi_printf("GRACKLE: ngas_active = %d\n",ngas_uncon);
      for(int j = 0; j < ngas_uncon; j++)
        mpi_printf("rho: %g u: %g HI: %g HII %g HeI %g HeII %g HeIII %g e: %g\n",density[j],energy[j],HI_density[j],HII_density[j],HeI_density[j],HeII_density[j],HeIII_density[j],e_density[j]);

      if(solve_chemistry(&my_grackle_units, All.cf_atime, dtime, grid_rank, grid_dimension, grid_start, 
                        grid_end, density, energy, x_velocity, y_velocity, z_velocity, HI_density, 
                        HII_density, HM_density, HeI_density, HeII_density, HeIII_density, H2I_density, 
                        H2II_density, DI_density, DII_density, HDI_density, e_density, metal_density) == 0)
        terminate("GRACKLE: Error in solve_chemistry.\n");

      ngas_uncon = 0;
      for(int i= 0; i < NumGas; i++)
        {
          if(i < 0)
            continue;
          if(con_list[i] == 1)
            continue;

          new_species_frac[0] = HI_density[ngas_uncon] / density[ngas_uncon];
          new_species_frac[1] = HII_density[ngas_uncon] / density[ngas_uncon];
          new_species_frac[2] = HeI_density[ngas_uncon] / density[ngas_uncon];
          new_species_frac[3] = HeII_density[ngas_uncon] / density[ngas_uncon];
          new_species_frac[4] = HeIII_density[ngas_uncon] / density[ngas_uncon];
          new_e_frac = e_density[ngas_uncon] / density[ngas_uncon] / (PROTONMASS / ELECTRONMASS);
#ifdef GRACKLE_H2
          new_species_frac[5] = HM_density[ngas_uncon] / density[ngas_uncon];
          new_species_frac[6] = H2I_density[ngas_uncon] / density[ngas_uncon];
          new_species_frac[7] = H2II_density[ngas_uncon] / density[ngas_uncon];
#endif
#ifdef GRACKLE_D
          new_species_frac[8] = DI_density[ngas_uncon] / density[ngas_uncon];
          new_species_frac[9] = DII_density[ngas_uncon] / density[ngas_uncon];
          new_species_frac[10] = HDI_density[ngas_uncon] / density[ngas_uncon];
#endif
          con_local = 1;
          if(abs((new_e_frac - SphP[i].e_frac) / SphP[i].e_frac) > 0.01)
            {
              con_local = 0;
            }
          else
            {
              for(int j = 0; j < GRACKLE_SPECIES_NUMBER; j++)
                {
                  if(abs((new_species_frac[j] - SphP[i].GrackleSpeciesFraction[j]) / SphP[i].GrackleSpeciesFraction[j]) > 0.01)
                    {
                      con_local = 0;
                      break;
                    }
                }
            }

          if(con_local = 1)
            {
              ngas_con++;
              con_list[i] = 1;

              for(int j = 0; j < GRACKLE_SPECIES_NUMBER; j++)
                {
                  SphP[i].GrackleSpeciesFraction[j] = new_species_frac[j];
                  SphP[i].GrackleSpeciesMass[j] = P[i].Mass * SphP[i].GrackleSpeciesFraction[j];
                }
              SphP[i].e_frac = new_e_frac;
              SphP[i].e_mass = P[i].Mass * SphP[i].e_frac;
            }

          ngas_uncon++;
        }

      if(iter++ > 10000)
        terminate("GRACKLE: unable to converge abundances");
    }                           /* End of while loop, exit when number of converged cells equals NumGas */

  myfree(con_list);
#ifdef GRACKLE_D
  myfree(HDI_density);
  myfree(DII_density);
  myfree(DI_density);
#endif
#ifdef GRACKLE_H2
  myfree(H2II_density);
  myfree(H2I_density);
  myfree(HM_density);
#endif
  myfree(e_density);
  myfree(HeIII_density);
  myfree(HeII_density);
  myfree(HeI_density);
  myfree(HII_density);
  myfree(HI_density);
  myfree(metal_density);
  myfree(energy);
  myfree(density);
}
#endif

#endif
#endif
