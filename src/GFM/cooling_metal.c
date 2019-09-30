/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/GFM/cooling_metal.c
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
#include <mpi.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include "../allvars.h"
#include "../proto.h"

#ifdef GFM_COOLING_METAL

static MyFloat *redshift_bins, *hydrogenNumberDensity_bins, *temperature_bins;
static unsigned int *n_redshift, *n_hydrogenNumberDensity, *n_temperature;
#ifndef RADCOOL
static MyFloat *metallicityInSolar_bins;
static unsigned int *n_metallicityInSolar;
#endif
#ifdef GFM_AGN_RADIATION
static MyFloat *bolFlux_bins;
static unsigned int *n_bolFlux;
static MyFloat *****netCoolingRate;
#else
#ifdef RADCOOL
static MyFloat *Phins_bins;
static MyFloat *Phios_bins;
static unsigned int *n_Phins;
static unsigned int *n_Phios;
static float *****netCoolingRate;
static float *****netHeatingRate;
#ifdef RADCOOL_HOTHALO
static MyFloat *redshift_bins_lz;
static unsigned int *n_redshift_lz;
static MyFloat *Phios_bins_lz;
static MyFloat *Phins_bins_lz;
static unsigned int *n_Phios_lz;
static unsigned int *n_Phins_lz;
static MyFloat *PhiT6_bins;
static unsigned int *n_PhiT6;
static MyFloat *PhiT7_bins;
static unsigned int *n_PhiT7;
static MyFloat *PhiT8_bins;
static unsigned int *n_PhiT8;
static float ********netCoolingRateLowZ;
static float ********netHeatingRateLowZ;
#endif
#else
static MyFloat ****netCoolingRate;
#endif
#endif

static int cur_redshift_low_index;
#ifdef RADCOOL_HOTHALO
static int cur_redshift_low_index_lz;
#endif

#ifndef HAVE_HDF5
#error "Need HAVE_HDF5 enabled to load metal cooling tables"
#endif

#if defined(RADCOOL) && defined(GFM_AGN_RADIATION)
#error "Only actitivate either one of RADCOOL or GFM_AGN_RADIATION"
#endif

#if defined(RADCOOL_HOTHALO) && !defined(RADCOOL)
#error "RADCOOL_HOTHALO needs RADCOOL option to be activated"
#endif

void set_pointer(void)
{
  redshift_bins = coolingMetalTable.Redshift_bins;
#ifdef GFM_AGN_RADIATION
  bolFlux_bins = coolingMetalTable.BolFlux_bins;
#endif
#ifdef RADCOOL
  Phios_bins = coolingMetalTable.Phios_bins;
  Phins_bins = coolingMetalTable.Phins_bins;
#ifdef RADCOOL_HOTHALO
  redshift_bins_lz = coolingMetalTable.Redshift_bins_lz;
  Phios_bins_lz = coolingMetalTable.Phios_bins_lz;
  Phins_bins_lz = coolingMetalTable.Phins_bins_lz;
  PhiT6_bins = coolingMetalTable.PhiT6_bins;
  PhiT7_bins = coolingMetalTable.PhiT7_bins;
  PhiT8_bins = coolingMetalTable.PhiT8_bins;
#endif
#endif
#ifndef RADCOOL
  metallicityInSolar_bins = coolingMetalTable.MetallicityInSolar_bins;
#endif
  hydrogenNumberDensity_bins = coolingMetalTable.HydrogenNumberDensity_bins;
  temperature_bins = coolingMetalTable.Temperature_bins;

  netCoolingRate = coolingMetalTable.NetCoolingRate;
#ifdef RADCOOL
  netHeatingRate = coolingMetalTable.NetHeatingRate;
#ifdef RADCOOL_HOTHALO
  netCoolingRateLowZ = coolingMetalTable.NetCoolingRateLowZ;
  netHeatingRateLowZ = coolingMetalTable.NetHeatingRateLowZ;
#endif
#endif
  n_redshift = &coolingMetalTable.N_Redshift;
#ifdef GFM_AGN_RADIATION
  n_bolFlux = &coolingMetalTable.N_BolFlux;
#endif
#ifdef RADCOOL
  n_Phios = &coolingMetalTable.N_Phios;
  n_Phins = &coolingMetalTable.N_Phins;
#ifdef RADCOOL_HOTHALO
  n_redshift_lz = &coolingMetalTable.N_Redshift_lz;
  n_Phios_lz = &coolingMetalTable.N_Phios_lz;
  n_Phins_lz = &coolingMetalTable.N_Phins_lz;
  n_PhiT6 = &coolingMetalTable.N_PhiT6;
  n_PhiT7 = &coolingMetalTable.N_PhiT7;
  n_PhiT8 = &coolingMetalTable.N_PhiT8;
#endif
#endif
#ifndef RADCOOL
  n_metallicityInSolar = &coolingMetalTable.N_MetallicityInSolar;
#endif
  n_hydrogenNumberDensity = &coolingMetalTable.N_HydrogenNumberDensity;
  n_temperature = &coolingMetalTable.N_Temperature;
}

void allocate_cooling_tables(void)
{
  int i, k;
#ifdef GFM_AGN_RADIATION
  int j;
#endif
#ifdef RADCOOL
  int j, n;
#ifdef RADCOOL_HOTHALO
  int idx6, idx7, idx8;
#endif
#else
  int m;
#endif


  /* allocate redshift bins */
  coolingMetalTable.Redshift_bins = (MyFloat *) mymalloc("redshift_bins", *n_redshift * sizeof(MyFloat));

#ifdef GFM_AGN_RADIATION
  /* allocate bolflux bins */
  coolingMetalTable.BolFlux_bins = (MyFloat *) mymalloc("bolFlux_bins", *n_bolFlux * sizeof(MyFloat));
#endif

#ifdef RADCOOL
  /*allocate the bins for new and old stars */
  coolingMetalTable.Phios_bins = (MyFloat *) mymalloc("Phios_bins", *n_Phios * sizeof(MyFloat));
  coolingMetalTable.Phins_bins = (MyFloat *) mymalloc("Phins_bins", *n_Phins * sizeof(MyFloat));
#ifdef RADCOOL_HOTHALO
  coolingMetalTable.Redshift_bins_lz = (MyFloat *) mymalloc("redshift_bins_lz", *n_redshift_lz * sizeof(MyFloat));
  coolingMetalTable.Phios_bins_lz = (MyFloat *) mymalloc("Phios_bins_lz", *n_Phios_lz * sizeof(MyFloat));
  coolingMetalTable.Phins_bins_lz = (MyFloat *) mymalloc("Phins_bins_lz", *n_Phins_lz * sizeof(MyFloat));
  coolingMetalTable.PhiT6_bins = (MyFloat *) mymalloc("PhiT6_bins", *n_PhiT6 * sizeof(MyFloat));
  coolingMetalTable.PhiT7_bins = (MyFloat *) mymalloc("PhiT7_bins", *n_PhiT7 * sizeof(MyFloat));
  coolingMetalTable.PhiT8_bins = (MyFloat *) mymalloc("PhiT8_bins", *n_PhiT8 * sizeof(MyFloat));
#endif
#endif
#ifndef RADCOOL
  /* allocate metallicity bins */
  coolingMetalTable.MetallicityInSolar_bins = (MyFloat *) mymalloc("metallicityInSolar_bins", *n_metallicityInSolar * sizeof(MyFloat));
#endif

  /* allocate hydrogen number density bins */
  coolingMetalTable.HydrogenNumberDensity_bins = (MyFloat *) mymalloc("hydrogenNumberDensity_bins", *n_hydrogenNumberDensity * sizeof(MyFloat));

  /* allocate temperature bins */
  coolingMetalTable.Temperature_bins = (MyFloat *) mymalloc("temperature_bins", *n_temperature * sizeof(MyFloat));

#ifdef GFM_AGN_RADIATION
  /* allocate cooling rate array */
  coolingMetalTable.NetCoolingRate = (MyFloat *****) mymalloc("netCoolingRate", 2 * sizeof(MyFloat ****));
  for(i = 0; i < 2; i++)
    {
      coolingMetalTable.NetCoolingRate[i] = (MyFloat ****) mymalloc("netCoolingRate", *n_bolFlux * sizeof(MyFloat ***));
      for(j = 0; j < *n_bolFlux; j++)
        {
          coolingMetalTable.NetCoolingRate[i][j] = (MyFloat ***) mymalloc("netCoolingRate", *n_metallicityInSolar * sizeof(MyFloat **));
          for(m = 0; m < *n_metallicityInSolar; m++)
            {
              coolingMetalTable.NetCoolingRate[i][j][m] = (MyFloat **) mymalloc("netCoolingRate", *n_hydrogenNumberDensity * sizeof(MyFloat *));
              for(k = 0; k < *n_hydrogenNumberDensity; k++)
                {
                  coolingMetalTable.NetCoolingRate[i][j][m][k] = (MyFloat *) mymalloc("netCoolingRate", *n_temperature * sizeof(MyFloat));
                }
            }
        }
    }
#else
#ifdef RADCOOL
  coolingMetalTable.NetCoolingRate = (float *****) mymalloc("netCoolingRate", 2 * sizeof(float ****));
  coolingMetalTable.NetHeatingRate = (float *****) mymalloc("netHeatingRate", 2 * sizeof(float ****));
  for(i = 0; i < 2; i++)
    {
      coolingMetalTable.NetCoolingRate[i] = (float ****) mymalloc("netCoolingRate", *n_Phios * sizeof(float ***));
      coolingMetalTable.NetHeatingRate[i] = (float ****) mymalloc("netHeatingRate", *n_Phios * sizeof(float ***));
      for(j = 0; j < *n_Phios; j++)
        {
          coolingMetalTable.NetCoolingRate[i][j] = (float ***) mymalloc("netCoolingRate", *n_Phins * sizeof(float **));
          coolingMetalTable.NetHeatingRate[i][j] = (float ***) mymalloc("netHeatingRate", *n_Phins * sizeof(float **));
          for(n = 0; n < *n_Phins; n++)
            {
              coolingMetalTable.NetCoolingRate[i][j][n] = (float **) mymalloc("netCoolingRate", *n_hydrogenNumberDensity * sizeof(float *));
              coolingMetalTable.NetHeatingRate[i][j][n] = (float **) mymalloc("netHeatingRate", *n_hydrogenNumberDensity * sizeof(float *));
              for(k = 0; k < *n_hydrogenNumberDensity; k++)
                {
                  coolingMetalTable.NetCoolingRate[i][j][n][k] = (float *) mymalloc("netCoolingRate", *n_temperature * sizeof(float));
                  coolingMetalTable.NetHeatingRate[i][j][n][k] = (float *) mymalloc("netHeatingRate", *n_temperature * sizeof(float));
                }
            }
        }
    }
  mpi_printf("GFM:METAL COOLING: RADCOOL : Allocated memory for the cooling table\n");

#ifdef RADCOOL_HOTHALO

  coolingMetalTable.NetCoolingRateLowZ = (float ********) mymalloc("netCoolingRateLowZ", 2 * sizeof(float *******));
  coolingMetalTable.NetHeatingRateLowZ = (float ********) mymalloc("netHeatingRateLowZ", 2 * sizeof(float *******));
  for(i = 0; i < 2; i++)
    {
      coolingMetalTable.NetCoolingRateLowZ[i] = (float *******) mymalloc("netCoolingRateLowZ", *n_Phios_lz * sizeof(float ******));
      coolingMetalTable.NetHeatingRateLowZ[i] = (float *******) mymalloc("netHeatingRateLowZ", *n_Phios_lz * sizeof(float ******));
      for(j = 0; j < *n_Phios_lz; j++)
        {
          coolingMetalTable.NetCoolingRateLowZ[i][j] = (float ******) mymalloc("netCoolingRateLowZ", *n_Phins_lz * sizeof(float *****));
          coolingMetalTable.NetHeatingRateLowZ[i][j] = (float ******) mymalloc("netHeatingRateLowZ", *n_Phins_lz * sizeof(float *****));
          for(n = 0; n < *n_Phins_lz; n++)
            {
              coolingMetalTable.NetCoolingRateLowZ[i][j][n] = (float *****) mymalloc("netCoolingRateLowZ", *n_PhiT6 * sizeof(float ****));
              coolingMetalTable.NetHeatingRateLowZ[i][j][n] = (float *****) mymalloc("netHeatingRateLowZ", *n_PhiT6 * sizeof(float ****));
              for(idx6 = 0; idx6 < *n_PhiT6; idx6++)
                {
                  coolingMetalTable.NetCoolingRateLowZ[i][j][n][idx6] = (float ****) mymalloc("netCoolingRateLowZ", *n_PhiT7 * sizeof(float ***));
                  coolingMetalTable.NetHeatingRateLowZ[i][j][n][idx6] = (float ****) mymalloc("netHeatingRateLowZ", *n_PhiT7 * sizeof(float ***));
                  for(idx7 = 0; idx7 < *n_PhiT7; idx7++)
                    {
                      coolingMetalTable.NetCoolingRateLowZ[i][j][n][idx6][idx7] = (float ***) mymalloc("netCoolingRateLowZ", *n_PhiT8 * sizeof(float **));
                      coolingMetalTable.NetHeatingRateLowZ[i][j][n][idx6][idx7] = (float ***) mymalloc("netHeatingRateLowZ", *n_PhiT8 * sizeof(float **));
                      for(idx8 = 0; idx8 < *n_PhiT8; idx8++)
                        {
                          coolingMetalTable.NetCoolingRateLowZ[i][j][n][idx6][idx7][idx8] = (float **) mymalloc("netCoolingRateLowZ", *n_hydrogenNumberDensity * sizeof(float *));
                          coolingMetalTable.NetHeatingRateLowZ[i][j][n][idx6][idx7][idx8] = (float **) mymalloc("netHeatingRateLowZ", *n_hydrogenNumberDensity * sizeof(float *));
                          for(k = 0; k < *n_hydrogenNumberDensity; k++)
                            {
                              coolingMetalTable.NetCoolingRateLowZ[i][j][n][idx6][idx7][idx8][k] = (float *) mymalloc("netCoolingRateLowZ", *n_temperature * sizeof(float));
                              coolingMetalTable.NetHeatingRateLowZ[i][j][n][idx6][idx7][idx8][k] = (float *) mymalloc("netHeatingRateLowZ", *n_temperature * sizeof(float));
                            }
                        }
                    }
                }
            }
        }
    }

  mpi_printf("GFM: METAL COOLING : RADCOOL_HOTHALO : Allocated memory for the low Z cooling file\n");
#endif

#else
  /* allocate cooling rate array */
  coolingMetalTable.NetCoolingRate = (MyFloat ****) mymalloc("netCoolingRate", 2 * sizeof(MyFloat ***));
  for(i = 0; i < 2; i++)
    {
      coolingMetalTable.NetCoolingRate[i] = (MyFloat ***) mymalloc("netCoolingRate", *n_metallicityInSolar * sizeof(MyFloat **));
      for(m = 0; m < *n_metallicityInSolar; m++)
        {
          coolingMetalTable.NetCoolingRate[i][m] = (MyFloat **) mymalloc("netCoolingRate", *n_hydrogenNumberDensity * sizeof(MyFloat *));

          for(k = 0; k < *n_hydrogenNumberDensity; k++)
            {
              coolingMetalTable.NetCoolingRate[i][m][k] = (MyFloat *) mymalloc("netCoolingRate", *n_temperature * sizeof(MyFloat));
            }
        }
    }
#endif
#endif
}



/* call this to allocate and read in cooling tables */
void read_cooling_tables_init(void)
{
  int i, k, l;
#ifdef GFM_AGN_RADIATION
  int j;
#endif
#ifdef RADCOOL
  int j, n;
#ifdef RADCOOL_HOTHALO
  int idx6, idx7, idx8;
#endif
#else
  int m;
#endif


  char fname[MAXLEN_PATH];
  hid_t file_id, dataset;
  /* note: need to free temp heap memory to avoid setting uninialized values in other structures */

#ifndef RADCOOL
  double *tmpdbl;
#else
  double *tmpdbl;
  float *tmpfl, *tmpflh;
#ifdef RADCOOL_HOTHALO
  float *tmpfl_lz, *tmpflh_lz;
#endif
#endif

#ifndef RADCOOL
  int tmpdbllen;
#else
  int tmpdbllen;
  int tmpfllen;
#ifdef RADCOOL_HOTHALO
  int tmpfllen_lz;
#endif
#endif

  sprintf(fname, "%s", All.CoolingTablePath);

  n_redshift = &coolingMetalTable.N_Redshift;
#ifdef GFM_AGN_RADIATION
  n_bolFlux = &coolingMetalTable.N_BolFlux;
#endif
#ifdef RADCOOL
  n_Phios = &coolingMetalTable.N_Phios;
  n_Phins = &coolingMetalTable.N_Phins;
#ifdef RADCOOL_HOTHALO
  n_redshift_lz = &coolingMetalTable.N_Redshift_lz;
  n_Phios_lz = &coolingMetalTable.N_Phios_lz;
  n_Phins_lz = &coolingMetalTable.N_Phins_lz;
  n_PhiT6 = &coolingMetalTable.N_PhiT6;
  n_PhiT7 = &coolingMetalTable.N_PhiT7;
  n_PhiT8 = &coolingMetalTable.N_PhiT8;
#endif
#endif
#ifndef RADCOOL
  n_metallicityInSolar = &coolingMetalTable.N_MetallicityInSolar;
#endif
  n_hydrogenNumberDensity = &coolingMetalTable.N_HydrogenNumberDensity;
  n_temperature = &coolingMetalTable.N_Temperature;


  mpi_printf("GFM_COOLING_METAL: initializing from file %s\n", fname);
  file_id = my_H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  dummyvar = 0;                 /* this redshift initialization is just added to fix a strange linking problem on Mac OSX - without it, the file cooling_metal_vars.c is somehow not linked in */

  dataset = my_H5Dopen(file_id, "N_Redshift");
  my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, n_redshift, "N_Redshift");
  my_H5Dclose(dataset, "N_Redshift");

  //cur_redshift_low_index = *n_redshift - 2;

  dataset = my_H5Dopen_if_existing(file_id, "N_BolFlux");
#ifdef GFM_AGN_RADIATION
  if(dataset < 0)
    terminate("GFM_AGN_RADIATION: unable to open dataset N_BolFlux\n");
  my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, n_bolFlux, "N_BolFlux");
  my_H5Dclose(dataset, "N_BolFlux");
#else
  if(dataset >= 0)
    terminate("GFM_COOLING_METAL: dataset N_BolFlux exists but GFM_AGN_RADIATION is off\n");
#endif

  dataset = my_H5Dopen_if_existing(file_id, "N_Phios");
#ifdef RADCOOL
  if(dataset < 0)
    terminate("RADCOOL: unable to open dataset N_Phios\n");
  my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, n_Phios, "N_Phios");
  my_H5Dclose(dataset, "N_Phios");
#else
  if(dataset >= 0)
    terminate("RADCOOL: dataset N_Phios exists but RADCOOL is off\n");
#endif

  dataset = my_H5Dopen_if_existing(file_id, "N_Phins");
#ifdef RADCOOL
  if(dataset < 0)
    terminate("RADCOOL: unable to open dataset N_Phins\n");
  my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, n_Phins, "N_Phins");
  my_H5Dclose(dataset, "N_Phins");
#else
  if(dataset >= 0)
    terminate("RADCOOL: dataset N_Phins exists but RADCOOL is off\n");
#endif


  dataset = my_H5Dopen_if_existing(file_id, "N_Redshift_lz");
#ifdef RADCOOL_HOTHALO
  if(dataset < 0)
    terminate("RADCOOL: unable to open dataset N_Redshift_lz\n");
  my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, n_redshift_lz, "N_Redshift_lz");
  my_H5Dclose(dataset, "N_Redshift_lz");
#else
  if(dataset >= 0)
    terminate("RADCOOL_HOTHALO: dataset N_Redshift_lz exists but RADCOOL_HOTHALO is off\n");
#endif

  dataset = my_H5Dopen_if_existing(file_id, "N_Phios_lz");
#ifdef RADCOOL_HOTHALO
  if(dataset < 0)
    terminate("RADCOOL: unable to open dataset N_Phios_lz\n");
  my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, n_Phios_lz, "N_Phios_lz");
  my_H5Dclose(dataset, "N_Phios_lz");
#else
  if(dataset >= 0)
    terminate("RADCOOL_HOTHALO: dataset N_Phios_lz exists but RADCOOL_HOTHALO is off\n");
#endif

  dataset = my_H5Dopen_if_existing(file_id, "N_Phins_lz");
#ifdef RADCOOL_HOTHALO
  if(dataset < 0)
    terminate("RADCOOL: unable to open dataset N_Phins_lz\n");
  my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, n_Phins_lz, "N_Phins_lz");
  my_H5Dclose(dataset, "N_Phins_lz");
#else
  if(dataset >= 0)
    terminate("RADCOOL_HOTHALO: dataset N_Phins_lz exists but RADCOOL_HOTHALO is off\n");
#endif

  dataset = my_H5Dopen_if_existing(file_id, "N_PhiT6");
#ifdef RADCOOL_HOTHALO
  if(dataset < 0)
    terminate("RADCOOL: unable to open dataset N_PhiT6\n");
  my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, n_PhiT6, "N_PhiT6");
  my_H5Dclose(dataset, "N_PhiT6");
#else
  if(dataset >= 0)
    terminate("RADCOOL_HOTHALO: dataset N_PhiT6 exists but RADCOOL_HOTHALO is off\n");
#endif

  dataset = my_H5Dopen_if_existing(file_id, "N_PhiT7");
#ifdef RADCOOL_HOTHALO
  if(dataset < 0)
    terminate("RADCOOL: unable to open dataset N_PhiT7\n");
  my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, n_PhiT7, "N_PhiT7");
  my_H5Dclose(dataset, "N_PhiT7");
#else
  if(dataset >= 0)
    terminate("RADCOOL_HOTHALO: dataset N_PhiT7 exists but RADCOOL_HOTHALO is off\n");
#endif

  dataset = my_H5Dopen_if_existing(file_id, "N_PhiT8");
#ifdef RADCOOL_HOTHALO
  if(dataset < 0)
    terminate("RADCOOL: unable to open dataset N_PhiT8\n");
  my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, n_PhiT8, "N_PhiT8");
  my_H5Dclose(dataset, "N_PhiT8");
#else
  if(dataset >= 0)
    terminate("RADCOOL_HOTHALO: dataset N_PhiT8 exists but RADCOOL_HOTHALO is off\n");
#endif


#ifndef RADCOOL
  dataset = my_H5Dopen(file_id, "N_MetallicityInSolar");
  my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, n_metallicityInSolar, "N_MetallicityInSolar");
  my_H5Dclose(dataset, "N_MetallicityInSolar");
#endif
  dataset = my_H5Dopen(file_id, "N_HydrogenNumberDensity");
  my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, n_hydrogenNumberDensity, "N_HydrogenNumberDensity");
  my_H5Dclose(dataset, "N_HydrogenNumberDensity");

  dataset = my_H5Dopen(file_id, "N_Temperature");
  my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, n_temperature, "N_Temperature");
  my_H5Dclose(dataset, "N_Temperature");

  cur_redshift_low_index = *n_redshift - 2;
#ifdef RADCOOL_HOTHALO
  cur_redshift_low_index_lz = *n_redshift_lz - 2;
#endif

  if(ThisTask == 0)
    {
      printf("GFM_COOLING_METAL: N_Redshift              = %d\n", *n_redshift);
#ifdef GFM_AGN_RADIATION
      printf("GFM_COOLING_METAL: N_BolFlux               = %d\n", *n_bolFlux);
#endif
#ifdef RADCOOL
      printf("RADCOOL:GFM_COOLING_METAL: N_Phios         = %d\n", *n_Phios);
      printf("RADCOOL:GFM_COOLING_METAL: N_Phins         = %d\n", *n_Phins);
#ifdef RADCOOL_HOTHALO
      printf("RADCOOL_HOTHALO:GFM_COOLING_METAL: N_Redshift_lz     = %d\n", *n_redshift_lz);
      printf("RADCOOL_HOTHALO:GFM_COOLING_METAL: N_Phios_lz        = %d\n", *n_Phios_lz);
      printf("RADCOOL_HOTHALO:GFM_COOLING_METAL: N_Phins_lz        = %d\n", *n_Phins_lz);
      printf("RADCOOL_HOTHALO:GFM_COOLING_METAL: N_PhiT6           = %d\n", *n_PhiT6);
      printf("RADCOOL_HOTHALO:GFM_COOLING_METAL: N_PhiT7           = %d\n", *n_PhiT7);
      printf("RADCOOL_HOTHALO:GFM_COOLING_METAL: N_PhiT8           = %d\n", *n_PhiT8);
#endif
#endif
#ifndef RADCOOL
      printf("GFM_COOLING_METAL: N_MetallicityInSolar    = %d\n", *n_metallicityInSolar);
#endif
      printf("GFM_COOLING_METAL: N_HydrogenNumberDensity = %d\n", *n_hydrogenNumberDensity);
      printf("GFM_COOLING_METAL: N_Temperature           = %d\n", *n_temperature);
    }
  my_H5Fclose(file_id, fname);

#ifdef GFM_AGN_RADIATION
  if(*n_metallicityInSolar == 1)
    terminate("GFM_COOLING_METAL: N_MetallicityInSolar == 1. We better stop since this is most likely not what you wanted.");
#endif

  allocate_cooling_tables();

  set_pointer();

  mpi_printf("GFM_COOLING_METAL: allocated cooling arrays\n");

  file_id = my_H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  tmpdbllen = *n_redshift * sizeof(double);
  tmpdbl = (double *) mymalloc("tmpdbl", tmpdbllen);
  dataset = my_H5Dopen(file_id, "Redshift_bins");
  mpi_printf("GFM_COOLING_METAL: reading %010d/%010d Bytes\n", H5Dget_storage_size(dataset), tmpdbllen);
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, "Redshift_bins");
  my_H5Dclose(dataset, "Redshift_bins");
  for(i = 0; i < *n_redshift; i++)
    redshift_bins[i] = tmpdbl[i];
  memset(tmpdbl, 0, tmpdbllen);
  myfree(tmpdbl);

#ifdef GFM_AGN_RADIATION
  tmpdbllen = *n_bolFlux * sizeof(double);
  tmpdbl = (double *) mymalloc("tmpdbl", tmpdbllen);
  dataset = my_H5Dopen(file_id, "BolFlux_bins");
  mpi_printf("GFM_COOLING_METAL: reading %010d/%010d Bytes\n", H5Dget_storage_size(dataset), tmpdbllen);
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, "BolFlux_bins");
  my_H5Dclose(dataset, "BolFlux_bins");
  for(i = 0; i < *n_bolFlux; i++)
    bolFlux_bins[i] = tmpdbl[i];
  memset(tmpdbl, 0, tmpdbllen);
  myfree(tmpdbl);
#endif


#ifdef RADCOOL
  tmpdbllen = *n_Phios * sizeof(double);
  tmpdbl = (double *) mymalloc("tmpdbl", tmpdbllen);
  dataset = my_H5Dopen(file_id, "Phios_bins");
  mpi_printf("GFM_COOLING_METAL: reading %010d/%010d Bytes\n", H5Dget_storage_size(dataset), tmpdbllen);
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, "Phios_bins");
  my_H5Dclose(dataset, "Phios_bins");
  for(i = 0; i < *n_Phios; i++)
    Phios_bins[i] = tmpdbl[i];
  memset(tmpdbl, 0, tmpdbllen);
  myfree(tmpdbl);


  tmpdbllen = *n_Phins * sizeof(double);
  tmpdbl = (double *) mymalloc("tmpdbl", tmpdbllen);
  dataset = my_H5Dopen(file_id, "Phins_bins");
  mpi_printf("GFM_COOLING_METAL: reading %010d/%010d Bytes\n", H5Dget_storage_size(dataset), tmpdbllen);
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, "Phins_bins");
  my_H5Dclose(dataset, "Phins_bins");
  for(i = 0; i < *n_Phins; i++)
    Phins_bins[i] = tmpdbl[i];
  memset(tmpdbl, 0, tmpdbllen);
  myfree(tmpdbl);

#ifdef RADCOOL_HOTHALO

  tmpdbllen = *n_redshift_lz * sizeof(double);
  tmpdbl = (double *) mymalloc("tmpdbl", tmpdbllen);
  dataset = my_H5Dopen(file_id, "Redshift_bins_lz");
  mpi_printf("GFM_COOLING_METAL: reading %010d/%010d Bytes\n", H5Dget_storage_size(dataset), tmpdbllen);
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, "Redshift_bins_lz");
  my_H5Dclose(dataset, "Redshift_bins_lz");
  for(i = 0; i < *n_redshift_lz; i++)
    redshift_bins_lz[i] = tmpdbl[i];
  memset(tmpdbl, 0, tmpdbllen);
  myfree(tmpdbl);


  tmpdbllen = *n_Phios_lz * sizeof(double);
  tmpdbl = (double *) mymalloc("tmpdbl", tmpdbllen);
  dataset = my_H5Dopen(file_id, "Phios_bins_lz");
  mpi_printf("GFM_COOLING_METAL: reading %010d/%010d Bytes\n", H5Dget_storage_size(dataset), tmpdbllen);
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, "Phios_bins_lz");
  my_H5Dclose(dataset, "Phios_bins_lz");
  for(i = 0; i < *n_Phios_lz; i++)
    Phios_bins_lz[i] = tmpdbl[i];
  memset(tmpdbl, 0, tmpdbllen);
  myfree(tmpdbl);


  tmpdbllen = *n_Phins_lz * sizeof(double);
  tmpdbl = (double *) mymalloc("tmpdbl", tmpdbllen);
  dataset = my_H5Dopen(file_id, "Phins_bins_lz");
  mpi_printf("GFM_COOLING_METAL: reading %010d/%010d Bytes\n", H5Dget_storage_size(dataset), tmpdbllen);
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, "Phins_bins_lz");
  my_H5Dclose(dataset, "Phins_bins_lz");
  for(i = 0; i < *n_Phins_lz; i++)
    Phins_bins_lz[i] = tmpdbl[i];
  memset(tmpdbl, 0, tmpdbllen);
  myfree(tmpdbl);

  tmpdbllen = *n_PhiT6 * sizeof(double);
  tmpdbl = (double *) mymalloc("tmpdbl", tmpdbllen);
  dataset = my_H5Dopen(file_id, "PhiT6_bins");
  mpi_printf("GFM_COOLING_METAL: reading %010d/%010d Bytes\n", H5Dget_storage_size(dataset), tmpdbllen);
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, "PhiT6_bins");
  my_H5Dclose(dataset, "PhiT6_bins");
  for(i = 0; i < *n_PhiT6; i++)
    PhiT6_bins[i] = tmpdbl[i];
  memset(tmpdbl, 0, tmpdbllen);
  myfree(tmpdbl);

  tmpdbllen = *n_PhiT7 * sizeof(double);
  tmpdbl = (double *) mymalloc("tmpdbl", tmpdbllen);
  dataset = my_H5Dopen(file_id, "PhiT7_bins");
  mpi_printf("GFM_COOLING_METAL: reading %010d/%010d Bytes\n", H5Dget_storage_size(dataset), tmpdbllen);
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, "PhiT7_bins");
  my_H5Dclose(dataset, "PhiT7_bins");
  for(i = 0; i < *n_PhiT7; i++)
    PhiT7_bins[i] = tmpdbl[i];
  memset(tmpdbl, 0, tmpdbllen);
  myfree(tmpdbl);

  tmpdbllen = *n_PhiT8 * sizeof(double);
  tmpdbl = (double *) mymalloc("tmpdbl", tmpdbllen);
  dataset = my_H5Dopen(file_id, "PhiT8_bins");
  mpi_printf("GFM_COOLING_METAL: reading %010d/%010d Bytes\n", H5Dget_storage_size(dataset), tmpdbllen);
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, "PhiT8_bins");
  my_H5Dclose(dataset, "PhiT8_bins");
  for(i = 0; i < *n_PhiT8; i++)
    PhiT8_bins[i] = tmpdbl[i];
  memset(tmpdbl, 0, tmpdbllen);
  myfree(tmpdbl);

#endif
#endif

#ifndef RADCOOL
  tmpdbllen = *n_metallicityInSolar * sizeof(double);
  tmpdbl = (double *) mymalloc("tmpdbl", tmpdbllen);
  dataset = my_H5Dopen(file_id, "MetallicityInSolar_bins");
  mpi_printf("GFM_COOLING_METAL: reading %010d/%010d Bytes\n", H5Dget_storage_size(dataset), tmpdbllen);
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, "MetallicityInSolar_bins");
  my_H5Dclose(dataset, "MetallicityInSolar_bins");
  for(i = 0; i < *n_metallicityInSolar; i++)
    metallicityInSolar_bins[i] = tmpdbl[i];
  memset(tmpdbl, 0, tmpdbllen);
  myfree(tmpdbl);
#endif

  tmpdbllen = *n_hydrogenNumberDensity * sizeof(double);
  tmpdbl = (double *) mymalloc("tmpdbl", tmpdbllen);
  dataset = my_H5Dopen(file_id, "HydrogenNumberDensity_bins");
  mpi_printf("GFM_COOLING_METAL: reading %010d/%010d Bytes\n", H5Dget_storage_size(dataset), tmpdbllen);
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, "HydrogenNumberDensity_bins");
  my_H5Dclose(dataset, "HydrogenNumberDensity_bins");
  for(i = 0; i < *n_hydrogenNumberDensity; i++)
    hydrogenNumberDensity_bins[i] = tmpdbl[i];
  memset(tmpdbl, 0, tmpdbllen);
  myfree(tmpdbl);

  tmpdbllen = *n_temperature * sizeof(double);
  tmpdbl = (double *) mymalloc("tmpdbl", tmpdbllen);
  dataset = my_H5Dopen(file_id, "Temperature_bins");
  mpi_printf("GFM_COOLING_METAL: reading %010d/%010d Bytes\n", H5Dget_storage_size(dataset), tmpdbllen);
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, "Temperature_bins");
  my_H5Dclose(dataset, "Temperature_bins");
  for(i = 0; i < *n_temperature; i++)
    temperature_bins[i] = tmpdbl[i];
  memset(tmpdbl, 0, tmpdbllen);
  myfree(tmpdbl);

#ifdef GFM_AGN_RADIATION
  tmpdbllen = 2 * *n_bolFlux * *n_metallicityInSolar * *n_hydrogenNumberDensity * *n_temperature * sizeof(double);
  hsize_t offset[5] = { cur_redshift_low_index, 0, 0, 0, 0 };
  hsize_t count[5] = { 2, *n_bolFlux, *n_metallicityInSolar, *n_hydrogenNumberDensity, *n_temperature };
#else
#ifdef RADCOOL
  tmpfllen = 2 * *n_Phios * *n_Phins * *n_hydrogenNumberDensity * *n_temperature * sizeof(float);
  hsize_t offset[5] = { cur_redshift_low_index, 0, 0, 0, 0 };
  hsize_t count[5] = { 2, *n_Phios, *n_Phins, *n_hydrogenNumberDensity, *n_temperature };
#ifdef RADCOOL_HOTHALO
  tmpfllen_lz = 2 * *n_Phios_lz * *n_Phins_lz * *n_PhiT6 * *n_PhiT7 * *n_PhiT8 * *n_hydrogenNumberDensity * *n_temperature * sizeof(float);
  hsize_t offset_lz[8] = { cur_redshift_low_index_lz, 0, 0, 0, 0, 0, 0, 0 };
  hsize_t count_lz[8] = { 2, *n_Phios_lz, *n_Phins_lz, *n_PhiT6, *n_PhiT7, *n_PhiT8, *n_hydrogenNumberDensity, *n_temperature };
#endif
#else
  tmpdbllen = 2 * *n_metallicityInSolar * *n_hydrogenNumberDensity * *n_temperature * sizeof(double);
  hsize_t offset[4] = { cur_redshift_low_index, 0, 0, 0 };
  hsize_t count[4] = { 2, *n_metallicityInSolar, *n_hydrogenNumberDensity, *n_temperature };
#endif
#endif

#ifndef RADCOOL
  tmpdbl = (double *) mymalloc("tmpdbl", tmpdbllen);
  dataset = my_H5Dopen(file_id, "NetCoolingRate");
  mpi_printf("GFM_COOLING_METAL: reading %010d/%010d Bytes\n", 2 * H5Dget_storage_size(dataset) / *n_redshift, tmpdbllen);

  hsize_t len = (hsize_t) (tmpdbllen / sizeof(double));

  hid_t filespace = my_H5Dget_space(dataset, "NetCoolingRate");
  hid_t memspace = my_H5Screate_simple(1, &len, NULL);
  my_H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, tmpdbl, "NetCoolingRate");
  my_H5Dclose(dataset, "NetCoolingRate");
#else
  tmpfl = (float *) mymalloc("tmpfl", tmpfllen);
  dataset = my_H5Dopen(file_id, "NetCoolingRate");
  mpi_printf("GFM_COOLING_METAL: reading %010d/%010d Bytes\n", 2 * H5Dget_storage_size(dataset) / *n_redshift, tmpfllen);

  hsize_t len = (hsize_t) (tmpfllen / sizeof(float));

  hid_t filespace = my_H5Dget_space(dataset, "NetCoolingRate");
  hid_t memspace = my_H5Screate_simple(1, &len, NULL);
  my_H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
  my_H5Dread(dataset, H5T_NATIVE_FLOAT, memspace, filespace, H5P_DEFAULT, tmpfl, "NetCoolingRate");
  my_H5Dclose(dataset, "NetCoolingRate");

  tmpflh = (float *) mymalloc("tmpflh", tmpfllen);
  dataset = my_H5Dopen(file_id, "NetHeatingRate");
  mpi_printf("GFM_COOLING_METAL: reading %010d/%010d Bytes\n", 2 * H5Dget_storage_size(dataset) / *n_redshift, tmpfllen);

  len = (hsize_t) (tmpfllen / sizeof(float));

  filespace = my_H5Dget_space(dataset, "NetHeatingRate");
  memspace = my_H5Screate_simple(1, &len, NULL);
  my_H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
  my_H5Dread(dataset, H5T_NATIVE_FLOAT, memspace, filespace, H5P_DEFAULT, tmpflh, "NetHeatingRate");
  my_H5Dclose(dataset, "NetHeatingRate");
#ifdef RADCOOL_HOTHALO
  tmpfl_lz = (float *) mymalloc("tmpfl_lz", tmpfllen_lz);
  dataset = my_H5Dopen(file_id, "NetCoolingRateLowZ");
  mpi_printf("GFM_COOLING_METAL: reading %010d/%010d Bytes\n", 2 * H5Dget_storage_size(dataset) / *n_redshift_lz, tmpfllen_lz);

  hsize_t len_lz = (hsize_t) (tmpfllen_lz / sizeof(float));

  hid_t filespace_lz = my_H5Dget_space(dataset, "NetCoolingRateLowZ");
  hid_t memspace_lz = my_H5Screate_simple(1, &len_lz, NULL);
  my_H5Sselect_hyperslab(filespace_lz, H5S_SELECT_SET, offset_lz, NULL, count_lz, NULL);
  my_H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_lz, filespace_lz, H5P_DEFAULT, tmpfl_lz, "NetCoolingRateLowZ");
  my_H5Dclose(dataset, "NetCoolingRateLowZ");

  tmpflh_lz = (float *) mymalloc("tmpflh_lz", tmpfllen_lz);
  dataset = my_H5Dopen(file_id, "NetHeatingRateLowZ");
  mpi_printf("GFM_COOLING_METAL: reading %010d/%010d Bytes\n", 2 * H5Dget_storage_size(dataset) / *n_redshift_lz, tmpfllen_lz);

  len_lz = (hsize_t) (tmpfllen_lz / sizeof(float));

  filespace_lz = my_H5Dget_space(dataset, "NetHeatingRateLowZ");
  memspace_lz = my_H5Screate_simple(1, &len_lz, NULL);
  my_H5Sselect_hyperslab(filespace_lz, H5S_SELECT_SET, offset_lz, NULL, count_lz, NULL);
  my_H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_lz, filespace_lz, H5P_DEFAULT, tmpflh_lz, "NetHeatingRateLowZ");
  my_H5Dclose(dataset, "NetHeatingRateLowZ");
#endif

#endif

  mpi_printf("GFM_COOLING_METAL: initializing cooling table with highest redshift bins, between z=%g and z=%g\n", redshift_bins[cur_redshift_low_index], redshift_bins[cur_redshift_low_index + 1]);
#ifdef RADCOOL_HOTHALO
  mpi_printf("GFM_COOLING_METAL: initializing cooling table with highest redshift bins - LOW Z, between z=%g and z=%g\n",
             redshift_bins_lz[cur_redshift_low_index_lz], redshift_bins_lz[cur_redshift_low_index_lz + 1]);
#endif

  int index = 0;

#ifdef GFM_AGN_RADIATION
  for(i = 0; i < 2; i++)
    for(j = 0; j < *n_bolFlux; j++)
      for(m = 0; m < *n_metallicityInSolar; m++)
        for(k = 0; k < *n_hydrogenNumberDensity; k++)
          for(l = 0; l < *n_temperature; l++)
            netCoolingRate[i][j][m][k][l] = tmpdbl[index++];
#else
#ifdef RADCOOL
  /* there is only one metallicity value in the table */
  for(i = 0; i < 2; i++)
    for(j = 0; j < *n_Phios; j++)
      for(n = 0; n < *n_Phins; n++)
        for(k = 0; k < *n_hydrogenNumberDensity; k++)
          for(l = 0; l < *n_temperature; l++)
            {
              netCoolingRate[i][j][n][k][l] = tmpfl[index];
              netHeatingRate[i][j][n][k][l] = tmpflh[index++];
            }
  mpi_printf("RADCOOL : Initialized the main tables \n");
#ifdef RADCOOL_HOTHALO
  for(i = 0; i < 2; i++)
    for(j = 0; j < *n_Phios_lz; j++)
      for(n = 0; n < *n_Phins_lz; n++)
        for(idx6 = 0; idx6 < *n_PhiT6; idx6++)
          for(idx7 = 0; idx7 < *n_PhiT7; idx7++)
            for(idx8 = 0; idx8 < *n_PhiT8; idx8++)
              for(k = 0; k < *n_hydrogenNumberDensity; k++)
                for(l = 0; l < *n_temperature; l++)
                  {
                    netCoolingRateLowZ[i][j][n][idx6][idx7][idx8][k][l] = tmpfl_lz[index];
                    netHeatingRateLowZ[i][j][n][idx6][idx7][idx8][k][l] = tmpflh_lz[index++];
                  }
  mpi_printf("RADCOOL HOTHALO: initialized LOW Z tables\n");
#endif
#else
  for(i = 0; i < 2; i++)
    for(m = 0; m < *n_metallicityInSolar; m++)
      for(k = 0; k < *n_hydrogenNumberDensity; k++)
        for(l = 0; l < *n_temperature; l++)
          netCoolingRate[i][m][k][l] = tmpdbl[index++];
#endif
#endif

#ifdef RADCOOL
#ifdef RADCOOL_HOTHALO
  memset(tmpflh_lz, 0, tmpfllen_lz);
  myfree(tmpflh_lz);
  memset(tmpfl_lz, 0, tmpfllen_lz);
  myfree(tmpfl_lz);
#endif
  memset(tmpflh, 0, tmpfllen);
  myfree(tmpflh);
  memset(tmpfl, 0, tmpfllen);
  myfree(tmpfl);
#else
  memset(tmpdbl, 0, tmpdbllen);
  myfree(tmpdbl);
#endif

  my_H5Fclose(file_id, fname);


  mpi_printf("GFM_COOLING_METAL: cooling arrays read in\n");
  mpi_printf("GFM_COOLING_METAL: initialized cooling table with highest redshift bins, between z=%g and z=%g\n", redshift_bins[cur_redshift_low_index], redshift_bins[cur_redshift_low_index + 1]);
#ifdef RADCOOL_HOTHALO
  mpi_printf("GFM_COOLING_METAL: initialized cooling table with highest redshift bins LOW Z, between z=%g and z=%g\n",
             redshift_bins_lz[cur_redshift_low_index_lz], redshift_bins_lz[cur_redshift_low_index_lz + 1]);
#endif


#if (DVR_RENDER==0)
  return;
#endif

  if(ThisTask == 0)
    {
      FILE *fd = NULL;
      char buf[1000];
#ifdef GFM_AGN_RADIATION
      sprintf(buf, "%s/cooling_metal_bins_AGN_Compton.txt", All.OutputDir);
#else
#ifdef RADCOOL
      sprintf(buf, "%s/cooling_metal_bins_RADCOOL.txt", All.OutputDir);
#ifdef RADCOOL_HOTHALO
      sprintf(buf, "%s/cooling_metal_bins_RADCOOL_HOTHALO.txt", All.OutputDir);
#endif
#else
      sprintf(buf, "%s/cooling_metal_bins_UVB.txt", All.OutputDir);
#endif
#endif

      if(!(fd = fopen(buf, "w")))
        terminate("can't open file '%s'\n", buf);

      fprintf(fd, "#Redshift bins:\n");
      for(i = 0; i < *n_redshift; i++)
        fprintf(fd, "%d %g\n", i, redshift_bins[i]);

#ifdef GFM_AGN_RADIATION
      fprintf(fd, "#BolFlux bins:\n");
      for(i = 0; i < *n_bolFlux; i++)
        fprintf(fd, "%d %g\n", i, bolFlux_bins[i]);
#endif

#ifdef RADCOOL
      fprintf(fd, "Phios bins:\n");
      for(i = 0; i < *n_Phios; i++)
        fprintf(fd, "%d %g\n", i, Phios_bins[i]);

      fprintf(fd, "Phins bins:\n");
      for(i = 0; i < *n_Phins; i++)
        fprintf(fd, "%d %g\n", i, Phins_bins[i]);
#ifdef RADCOOL_HOTHALO

      fprintf(fd, "#Low Redshift bins:\n");
      for(i = 0; i < *n_redshift_lz; i++)
        fprintf(fd, "%d %g\n", i, redshift_bins_lz[i]);

      fprintf(fd, "Phios Low Redshift bins:\n");
      for(i = 0; i < *n_Phios_lz; i++)
        fprintf(fd, "%d %g\n", i, Phios_bins_lz[i]);

      fprintf(fd, "Phins Low Redshift bins:\n");
      for(i = 0; i < *n_Phins_lz; i++)
        fprintf(fd, "%d %g\n", i, Phins_bins_lz[i]);

      fprintf(fd, "PhiT6 bins:\n");
      for(i = 0; i < *n_PhiT6; i++)
        fprintf(fd, "%d %g\n", i, PhiT6_bins[i]);

      fprintf(fd, "PhiT7 bins:\n");
      for(i = 0; i < *n_PhiT7; i++)
        fprintf(fd, "%d %g\n", i, PhiT7_bins[i]);

      fprintf(fd, "PhiT8 bins:\n");
      for(i = 0; i < *n_PhiT8; i++)
        fprintf(fd, "%d %g\n", i, PhiT8_bins[i]);
#endif
#endif
#ifndef RADCOOL
      fprintf(fd, "#MetallicityInSolar bins:\n");
      for(i = 0; i < *n_metallicityInSolar; i++)
        fprintf(fd, "%d %g\n", i, metallicityInSolar_bins[i]);
#endif
      fprintf(fd, "#HydrogenNumberDensity bins:\n");
      for(i = 0; i < *n_hydrogenNumberDensity; i++)
        fprintf(fd, "%d %g\n", i, hydrogenNumberDensity_bins[i]);

      fprintf(fd, "#Temperature_bins:\n");
      for(i = 0; i < *n_temperature; i++)
        fprintf(fd, "%d %g\n", i, temperature_bins[i]);

      fclose(fd);
    }
}


/* call this to read in portion of cooling tables for current time */
void read_cooling_tables_current_time(void)
{
  int i, k, l;
#ifdef GFM_AGN_RADIATION
  int j;
#endif
#ifdef RADCOOL
  int j, n;
#ifdef RADCOOL_HOTHALO
  int idx6, idx7, idx8;
#endif
#else
  int m;
#endif
  char fname[MAXLEN_PATH];
  hid_t file_id, dataset;
#ifndef RADCOOL
  double *tmpdbl;
  int tmpdbllen;
#endif

#ifdef RADCOOL
  float *tmpfl, *tmpflh;
  int tmpfllen;
#endif

#ifndef RADCOOL_HOTHALO
  if(cur_redshift_low_index > 0 && All.cf_redshift < redshift_bins[cur_redshift_low_index])
    {
      while(cur_redshift_low_index > 0 && All.cf_redshift < redshift_bins[cur_redshift_low_index])
        cur_redshift_low_index--;

      sprintf(fname, "%s", All.CoolingTablePath);
      mpi_printf("GFM_COOLING_METAL: reading cooling table for current time, between z=%g and z=%g\n", redshift_bins[cur_redshift_low_index], redshift_bins[cur_redshift_low_index + 1]);

      file_id = my_H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

      dummyvar = 0;             /* this redshift initialization is just added to fix a strange linking problem on Mac OSX - without it, the file cooling_metal_vars.c is somehow not linked in */
#ifdef GFM_AGN_RADIATION
      tmpdbllen = 2 * *n_bolFlux * *n_metallicityInSolar * *n_hydrogenNumberDensity * *n_temperature * sizeof(double);
      hsize_t offset[5] = { cur_redshift_low_index, 0, 0, 0, 0 };
      hsize_t count[5] = { 2, *n_bolFlux, *n_metallicityInSolar, *n_hydrogenNumberDensity, *n_temperature };
#else
#ifdef RADCOOL
      tmpfllen = 2 * *n_Phios * *n_Phins * *n_hydrogenNumberDensity * *n_temperature * sizeof(float);
      hsize_t offset[5] = { cur_redshift_low_index, 0, 0, 0, 0 };
      hsize_t count[5] = { 2, *n_Phios, *n_Phins, *n_hydrogenNumberDensity, *n_temperature };
#else
      tmpdbllen = 2 * *n_metallicityInSolar * *n_hydrogenNumberDensity * *n_temperature * sizeof(double);
      hsize_t offset[4] = { cur_redshift_low_index, 0, 0, 0 };
      hsize_t count[4] = { 2, *n_metallicityInSolar, *n_hydrogenNumberDensity, *n_temperature };
#endif
#endif

#ifndef RADCOOL
      tmpdbl = (double *) mymalloc("tmpdbl", tmpdbllen);
      dataset = my_H5Dopen(file_id, "NetCoolingRate");
      mpi_printf("GFM_COOLING_METAL: reading %010d/%010d Bytes\n", 2 * H5Dget_storage_size(dataset) / *n_redshift, tmpdbllen);

      hsize_t len = (hsize_t) (tmpdbllen / sizeof(double));

      hid_t filespace = my_H5Dget_space(dataset, "NetCoolingRate");
      hid_t memspace = my_H5Screate_simple(1, &len, NULL);
      my_H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
      my_H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, tmpdbl, "NetCoolingRate");
      my_H5Dclose(dataset, "NetCoolingRate");
#else
      tmpfl = (float *) mymalloc("tmpfl", tmpfllen);
      dataset = my_H5Dopen(file_id, "NetCoolingRate");
      mpi_printf("GFM_COOLING_METAL: reading %010d/%010d Bytes\n", 2 * H5Dget_storage_size(dataset) / *n_redshift, tmpfllen);

      hsize_t len = (hsize_t) (tmpfllen / sizeof(float));

      hid_t filespace = my_H5Dget_space(dataset, "NetCoolingRate");
      hid_t memspace = my_H5Screate_simple(1, &len, NULL);
      my_H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
      my_H5Dread(dataset, H5T_NATIVE_FLOAT, memspace, filespace, H5P_DEFAULT, tmpfl, "NetCoolingRate");
      my_H5Dclose(dataset, "NetCoolingRate");

      tmpflh = (float *) mymalloc("tmpflh", tmpfllen);
      dataset = my_H5Dopen(file_id, "NetHeatingRate");
      mpi_printf("GFM_COOLING_METAL: reading %010d/%010d Bytes\n", 2 * H5Dget_storage_size(dataset) / *n_redshift, tmpfllen);

      len = (hsize_t) (tmpfllen / sizeof(float));

      filespace = my_H5Dget_space(dataset, "NetHeatingRate");
      memspace = my_H5Screate_simple(1, &len, NULL);
      my_H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
      my_H5Dread(dataset, H5T_NATIVE_FLOAT, memspace, filespace, H5P_DEFAULT, tmpflh, "NetHeatingRate");
      my_H5Dclose(dataset, "NetHeatingRate");
#endif

      int index = 0;

#ifdef GFM_AGN_RADIATION
      for(i = 0; i < 2; i++)
        for(j = 0; j < *n_bolFlux; j++)
          for(m = 0; m < *n_metallicityInSolar; m++)
            for(k = 0; k < *n_hydrogenNumberDensity; k++)
              for(l = 0; l < *n_temperature; l++)
                netCoolingRate[i][j][m][k][l] = tmpdbl[index++];
#else
      /* there is only one metallicity value in the table */
#ifdef RADCOOL
      for(i = 0; i < 2; i++)
        for(j = 0; j < *n_Phios; j++)
          for(n = 0; n < *n_Phins; n++)
            for(k = 0; k < *n_hydrogenNumberDensity; k++)
              for(l = 0; l < *n_temperature; l++)
                {
                  netCoolingRate[i][j][n][k][l] = tmpfl[index];
                  netHeatingRate[i][j][n][k][l] = tmpflh[index++];
                }
#else
      for(i = 0; i < 2; i++)
        for(m = 0; m < *n_metallicityInSolar; m++)
          for(k = 0; k < *n_hydrogenNumberDensity; k++)
            for(l = 0; l < *n_temperature; l++)
              netCoolingRate[i][m][k][l] = tmpdbl[index++];
#endif
#endif

#ifndef RADCOOL
      memset(tmpdbl, 0, tmpdbllen);
      myfree(tmpdbl);
#else
      memset(tmpflh, 0, tmpfllen);
      myfree(tmpflh);
      memset(tmpfl, 0, tmpfllen);
      myfree(tmpfl);
#endif
      my_H5Fclose(file_id, fname);
    }
#else
  if(All.cf_redshift < HOTHALO_ON_zFACTOR)
    {
      if(cur_redshift_low_index_lz > 0 && All.cf_redshift < redshift_bins_lz[cur_redshift_low_index_lz])
        {
          while(cur_redshift_low_index_lz > 0 && All.cf_redshift < redshift_bins_lz[cur_redshift_low_index_lz])
            cur_redshift_low_index_lz--;

          sprintf(fname, "%s", All.CoolingTablePath);
          mpi_printf("GFM_COOLING_METAL: reading cooling table for current time, between z=%g and z=%g\n",
                     redshift_bins_lz[cur_redshift_low_index_lz], redshift_bins_lz[cur_redshift_low_index_lz + 1]);


          file_id = my_H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

          dummyvar = 0;         /* this redshift initialization is just added to fix a strange linking problem on Mac OSX - without it, the file cooling_metal_vars.c is somehow not linked in */

          tmpfllen = 2 * *n_Phios_lz * *n_Phins_lz * *n_PhiT6 * *n_PhiT7 * *n_PhiT8 * *n_hydrogenNumberDensity * *n_temperature * sizeof(float);
          hsize_t offset[8] = { cur_redshift_low_index_lz, 0, 0, 0, 0, 0, 0, 0 };
          hsize_t count[8] = { 2, *n_Phios_lz, *n_Phins_lz, *n_PhiT6, *n_PhiT7, *n_PhiT8, *n_hydrogenNumberDensity, *n_temperature };

          tmpfl = (float *) mymalloc("tmpfl", tmpfllen);
          dataset = my_H5Dopen(file_id, "NetCoolingRateLowZ");
          mpi_printf("GFM_COOLING_METAL: reading %010d/%010d Bytes\n", 2 * H5Dget_storage_size(dataset) / *n_redshift_lz, tmpfllen);

          hsize_t len = (hsize_t) (tmpfllen / sizeof(float));

          hid_t filespace = my_H5Dget_space(dataset, "NetCoolingRateLowZ");
          hid_t memspace = my_H5Screate_simple(1, &len, NULL);
          my_H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
          my_H5Dread(dataset, H5T_NATIVE_FLOAT, memspace, filespace, H5P_DEFAULT, tmpfl, "NetCoolingRateLowZ");
          my_H5Dclose(dataset, "NetCoolingRateLowZ");

          tmpflh = (float *) mymalloc("tmpflh", tmpfllen);
          dataset = my_H5Dopen(file_id, "NetHeatingRateLowZ");
          mpi_printf("GFM_COOLING_METAL: reading %010d/%010d Bytes\n", 2 * H5Dget_storage_size(dataset) / *n_redshift_lz, tmpfllen);

          len = (hsize_t) (tmpfllen / sizeof(float));

          filespace = my_H5Dget_space(dataset, "NetHeatingRateLowZ");
          memspace = my_H5Screate_simple(1, &len, NULL);
          my_H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
          my_H5Dread(dataset, H5T_NATIVE_FLOAT, memspace, filespace, H5P_DEFAULT, tmpflh, "NetHeatingRateLowZ");
          my_H5Dclose(dataset, "NetHeatingRateLowZ");

          int index = 0;

          for(i = 0; i < 2; i++)
            for(j = 0; j < *n_Phios_lz; j++)
              for(n = 0; n < *n_Phins_lz; n++)
                for(idx6 = 0; idx6 < *n_PhiT6; idx6++)
                  for(idx7 = 0; idx7 < *n_PhiT7; idx7++)
                    for(idx8 = 0; idx8 < *n_PhiT8; idx8++)
                      for(k = 0; k < *n_hydrogenNumberDensity; k++)
                        for(l = 0; l < *n_temperature; l++)
                          {
                            netCoolingRateLowZ[i][j][n][idx6][idx7][idx8][k][l] = tmpfl[index];
                            netHeatingRateLowZ[i][j][n][idx6][idx7][idx8][k][l] = tmpflh[index++];
                          }

          memset(tmpflh, 0, tmpfllen);
          myfree(tmpflh);
          memset(tmpfl, 0, tmpfllen);
          myfree(tmpfl);
          my_H5Fclose(file_id, fname);
        }
    }
  else
    {
      if(cur_redshift_low_index > 0 && All.cf_redshift < redshift_bins[cur_redshift_low_index])
        {
          while(cur_redshift_low_index > 0 && All.cf_redshift < redshift_bins[cur_redshift_low_index])
            cur_redshift_low_index--;

          sprintf(fname, "%s", All.CoolingTablePath);
          mpi_printf("GFM_COOLING_METAL: reading cooling table for current time, between z=%g and z=%g\n", redshift_bins[cur_redshift_low_index], redshift_bins[cur_redshift_low_index + 1]);



          file_id = my_H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

          dummyvar = 0;         /* this redshift initialization is just added to fix a strange linking problem on Mac OSX - without it, the file cooling_metal_vars.c is somehow not linked in */
          tmpfllen = 2 * *n_Phios * *n_Phins * *n_hydrogenNumberDensity * *n_temperature * sizeof(float);
          hsize_t offset[5] = { cur_redshift_low_index, 0, 0, 0, 0 };
          hsize_t count[5] = { 2, *n_Phios, *n_Phins, *n_hydrogenNumberDensity, *n_temperature };

          tmpfl = (float *) mymalloc("tmpfl", tmpfllen);
          dataset = my_H5Dopen(file_id, "NetCoolingRate");
          mpi_printf("GFM_COOLING_METAL: reading %010d/%010d Bytes\n", 2 * H5Dget_storage_size(dataset) / *n_redshift, tmpfllen);

          hsize_t len = (hsize_t) (tmpfllen / sizeof(float));

          hid_t filespace = my_H5Dget_space(dataset, "NetCoolingRate");
          hid_t memspace = my_H5Screate_simple(1, &len, NULL);
          my_H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
          my_H5Dread(dataset, H5T_NATIVE_FLOAT, memspace, filespace, H5P_DEFAULT, tmpfl, "NetCoolingRate");
          my_H5Dclose(dataset, "NetCoolingRate");

          tmpflh = (float *) mymalloc("tmpflh", tmpfllen);
          dataset = my_H5Dopen(file_id, "NetHeatingRate");
          mpi_printf("GFM_COOLING_METAL: reading %010d/%010d Bytes\n", 2 * H5Dget_storage_size(dataset) / *n_redshift, tmpfllen);

          len = (hsize_t) (tmpfllen / sizeof(float));

          filespace = my_H5Dget_space(dataset, "NetHeatingRate");
          memspace = my_H5Screate_simple(1, &len, NULL);
          my_H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
          my_H5Dread(dataset, H5T_NATIVE_FLOAT, memspace, filespace, H5P_DEFAULT, tmpflh, "NetHeatingRate");
          my_H5Dclose(dataset, "NetHeatingRate");

          int index = 0;
          /* there is only one metallicity value in the table */
          for(i = 0; i < 2; i++)
            for(j = 0; j < *n_Phios; j++)
              for(n = 0; n < *n_Phins; n++)
                for(k = 0; k < *n_hydrogenNumberDensity; k++)
                  for(l = 0; l < *n_temperature; l++)
                    {
                      netCoolingRate[i][j][n][k][l] = tmpfl[index];
                      netHeatingRate[i][j][n][k][l] = tmpflh[index++];
                    }

          memset(tmpflh, 0, tmpfllen);
          myfree(tmpflh);
          memset(tmpfl, 0, tmpfllen);
          myfree(tmpfl);
          my_H5Fclose(file_id, fname);
        }
    }
#endif
}


/* 
 returns metal net cooling rate as Rate/nH^2 
 log_BolFlux, = log10(F/(erg/s/cm^2))
 log_MetallicityInSolar = log10(Z/Z_solar)
 log_HydrogenNumberDensity = log10(nH/(cm^-3))
 T comes from primordial network 
*/
#ifdef GFM_AGN_RADIATION
MyFloat get_CoolingMetalRate(MyFloat log_BolFlux, MyFloat log_MetallicityInSolar, MyFloat log_HydrogenNumberDensity, MyFloat log_Temperature)
#else
#ifdef RADCOOL
MyFloat get_CoolingMetalRate(MyFloat log_Phios, MyFloat log_Phins
#ifdef RADCOOL_HOTHALO
                             , MyFloat log_PhiT6, MyFloat log_PhiT7, MyFloat log_PhiT8
#endif
                             , MyFloat log_MetallicityInSolar, MyFloat log_HydrogenNumberDensity, MyFloat log_Temperature)
#else
MyFloat get_CoolingMetalRate(MyFloat log_MetallicityInSolar, MyFloat log_HydrogenNumberDensity, MyFloat log_Temperature)
#endif
#endif
{
  MyFloat coolingMetalRate = 0.0;
  double delta_i, d_i_local, delta_k, d_k_local, delta_l, d_l_local;
  double d_i, d_k, d_l, d_m;
  int i1, i2, k1, k2, l1, l2, m1;

#ifdef GFM_AGN_RADIATION
  double delta_j, d_j_local;
  double d_j;
  int j1, j2;
  double delta_m, d_m_local;
  int m2;
#endif

#ifdef RADCOOL
  double delta_j, d_j_local;
  double d_j;
  int j1, j2;
  double delta_n, d_n_local;
  double d_n;
  int n1, n2;
#ifdef RADCOOL_HOTHALO
  double delta_i6, d_i6_local;
  double d_i6;
  int i61, i62;
  double delta_i7, d_i7_local;
  double d_i7;
  int i71, i72;
  double delta_i8, d_i8_local;
  double d_i8;
  int i81, i82;
#endif
#endif


#ifdef RADCOOL
  double cool, heat;
#endif


/* no metal line cooling beyond the largest redshift */
  if(All.cf_redshift > redshift_bins[*n_redshift - 1])
    return 0.0;

/* no metal line cooling outside of temperature range */
  if(log_Temperature < temperature_bins[0] || log_Temperature > temperature_bins[*n_temperature - 1])
    return 0.0;

  /* search redshift */
#ifndef RADCOOL_HOTHALO
  if(All.cf_redshift > redshift_bins[0])
    {
      for(i1 = *n_redshift - 1; i1 > 0 && All.cf_redshift < redshift_bins[i1]; i1--);
      i2 = i1 + 1;
      if(i2 >= *n_redshift)
        i2 = *n_redshift - 1;
      i1 -= cur_redshift_low_index;
      i2 -= cur_redshift_low_index;
      if(i1 != 0)
        terminate("GFM_COOLING_METAL: i1 != 0: i1=%d i2=%d cur_redshift_low_index=%d redshift_bins[i1+cur_redshift_low_index]=%g", i1, i2,
                  cur_redshift_low_index, redshift_bins[i1 + cur_redshift_low_index]);
      if(i2 != 1)
        terminate("GFM_COOLING_METAL: i2 != 1: i1=%d i2=%d cur_redshift_low_index=%d redshift_bins[i2+cur_redshift_low_index]=%g", i1, i2,
                  cur_redshift_low_index, redshift_bins[i2 + cur_redshift_low_index]);
      if(All.cf_redshift >= redshift_bins[0] && All.cf_redshift <= redshift_bins[*n_redshift - 1])
        d_i_local = All.cf_redshift - redshift_bins[i1 + cur_redshift_low_index];
      else
        d_i_local = 0;
      delta_i = redshift_bins[i2 + cur_redshift_low_index] - redshift_bins[i1 + cur_redshift_low_index];
      if(delta_i > 0)
        d_i_local = d_i_local / delta_i;
      else
        d_i_local = 0;
    }
  else
    {
      i1 = 0;
      i2 = 0;
      d_i_local = 0.0;
    }
  d_i = d_i_local;
#else
  if(All.cf_redshift < HOTHALO_ON_zFACTOR)
    {
      if(All.cf_redshift > redshift_bins_lz[0])
        {
          for(i1 = *n_redshift_lz - 1; i1 > 0 && All.cf_redshift < redshift_bins_lz[i1]; i1--);
          i2 = i1 + 1;

          if(i2 >= *n_redshift_lz)
            i2 = *n_redshift_lz - 1;

          i1 -= cur_redshift_low_index_lz;
          i2 -= cur_redshift_low_index_lz;

          if(i1 != 0)
            terminate("GFM_COOLING_METAL: i1 != 0: i1=%d i2=%d cur_redshift_low_index=%d redshift_bins_lz[i1+cur_redshift_low_index]=%g", i1, i2,
                      cur_redshift_low_index_lz, redshift_bins_lz[i1 + cur_redshift_low_index_lz]);
          if(i2 != 1)
            terminate("GFM_COOLING_METAL: i2 != 1: i1=%d i2=%d cur_redshift_low_index=%d redshift_bins_lz[i2+cur_redshift_low_index]=%g", i1, i2,
                      cur_redshift_low_index_lz, redshift_bins_lz[i2 + cur_redshift_low_index_lz]);
          if(All.cf_redshift >= redshift_bins_lz[0] && All.cf_redshift <= redshift_bins_lz[*n_redshift_lz - 1])
            d_i_local = All.cf_redshift - redshift_bins_lz[i1 + cur_redshift_low_index_lz];
          else
            d_i_local = 0;
          delta_i = redshift_bins_lz[i2 + cur_redshift_low_index_lz] - redshift_bins_lz[i1 + cur_redshift_low_index_lz];
          if(delta_i > 0)
            d_i_local = d_i_local / delta_i;
          else
            d_i_local = 0;
        }
      else
        {
          i1 = 0;
          i2 = 0;
          d_i_local = 0.0;
        }
      d_i = d_i_local;
    }
  else
    {
      if(All.cf_redshift > redshift_bins[0])
        {
          for(i1 = *n_redshift - 1; i1 > 0 && All.cf_redshift < redshift_bins[i1]; i1--);
          i2 = i1 + 1;
          if(i2 >= *n_redshift)
            i2 = *n_redshift - 1;
          i1 -= cur_redshift_low_index;
          i2 -= cur_redshift_low_index;
          if(i1 != 0)
            terminate("GFM_COOLING_METAL: i1 != 0: i1=%d i2=%d cur_redshift_low_index=%d redshift_bins[i1+cur_redshift_low_index]=%g", i1, i2,
                      cur_redshift_low_index, redshift_bins[i1 + cur_redshift_low_index]);
          if(i2 != 1)
            terminate("GFM_COOLING_METAL: i2 != 1: i1=%d i2=%d cur_redshift_low_index=%d redshift_bins[i2+cur_redshift_low_index]=%g", i1, i2,
                      cur_redshift_low_index, redshift_bins[i2 + cur_redshift_low_index]);
          if(All.cf_redshift >= redshift_bins[0] && All.cf_redshift <= redshift_bins[*n_redshift - 1])
            d_i_local = All.cf_redshift - redshift_bins[i1 + cur_redshift_low_index];
          else
            d_i_local = 0;
          delta_i = redshift_bins[i2 + cur_redshift_low_index] - redshift_bins[i1 + cur_redshift_low_index];
          if(delta_i > 0)
            d_i_local = d_i_local / delta_i;
          else
            d_i_local = 0;
        }
      else
        {
          i1 = 0;
          i2 = 0;
          d_i_local = 0.0;
        }
      d_i = d_i_local;
    }
#endif


#ifdef GFM_AGN_RADIATION
  /* search bol. flux */
  if(log_BolFlux > bolFlux_bins[0])
    {
      for(j1 = 0; j1 < *n_bolFlux - 2 && log_BolFlux > bolFlux_bins[j1 + 1]; j1++);
      j2 = j1 + 1;
      if(j2 >= *n_bolFlux)
        j2 = *n_bolFlux - 1;
      if(log_BolFlux >= bolFlux_bins[0] && log_BolFlux <= bolFlux_bins[*n_bolFlux - 1])
        d_j_local = log_BolFlux - bolFlux_bins[j1];
      else
        d_j_local = 0;
      delta_j = bolFlux_bins[j2] - bolFlux_bins[j1];
      if(delta_j > 0)
        d_j_local = d_j_local / delta_j;
      else
        d_j_local = 0;
    }
  else
    {
      j1 = 0;
      j2 = 0;
      d_j_local = 0.0;
    }
  d_j = d_j_local;
#endif

#if defined (RADCOOL) && !defined (RADCOOL_HOTHALO)
  /* search Phios */
  if(log_Phios > Phios_bins[0])
    {
      for(j1 = 0; j1 < *n_Phios - 2 && log_Phios > Phios_bins[j1 + 1]; j1++);
      j2 = j1 + 1;
      if(j2 >= *n_Phios)
        j2 = *n_Phios - 1;
      if(log_Phios >= Phios_bins[1] && log_Phios <= Phios_bins[*n_Phios - 1])
        d_j_local = log_Phios - Phios_bins[j1];
      else
        d_j_local = 0;
      delta_j = Phios_bins[j2] - Phios_bins[j1];
      if(delta_j > 0)
        d_j_local = d_j_local / delta_j;
      else
        d_j_local = 0;
    }
  else
    {
      j1 = 0;
      j2 = 0;
      d_j_local = 0.0;
    }
  d_j = d_j_local;

  /*search Phins */
  if(log_Phins > Phins_bins[0])
    {
      for(n1 = 0; n1 < *n_Phins - 2 && log_Phins > Phins_bins[n1 + 1]; n1++);
      n2 = n1 + 1;
      if(n2 >= *n_Phins)
        n2 = *n_Phins - 1;
      if(log_Phins >= Phins_bins[1] && log_Phins <= Phins_bins[*n_Phins - 1])
        d_n_local = log_Phins - Phins_bins[n1];
      else
        d_n_local = 0;
      delta_n = Phins_bins[n2] - Phins_bins[n1];
      if(delta_n > 0)
        d_n_local = d_n_local / delta_n;
      else
        d_n_local = 0;
    }
  else
    {
      n1 = 0;
      n2 = 0;
      d_n_local = 0.0;
    }
  d_n = d_n_local;
#else
#if defined (RADCOOL) && defined (RADCOOL_HOTHALO)
  if(All.cf_redshift < HOTHALO_ON_zFACTOR)
    {
      /* search Phios */
      if(log_Phios > Phios_bins_lz[0])
        {
          for(j1 = 0; j1 < *n_Phios_lz - 2 && log_Phios > Phios_bins_lz[j1 + 1]; j1++);
          j2 = j1 + 1;
          if(j2 >= *n_Phios_lz)
            j2 = *n_Phios_lz - 1;
          if(log_Phios >= Phios_bins_lz[0] && log_Phios <= Phios_bins_lz[*n_Phios_lz - 1])
            d_j_local = log_Phios - Phios_bins_lz[j1];
          else
            d_j_local = 0;
          delta_j = Phios_bins_lz[j2] - Phios_bins_lz[j1];
          if(delta_j > 0)
            d_j_local = d_j_local / delta_j;
          else
            d_j_local = 0;
        }
      else
        {
          j1 = 0;
          j2 = 0;
          d_j_local = 0.0;
        }
      d_j = d_j_local;

      /*search Phins */
      if(log_Phins > Phins_bins_lz[0])
        {
          for(n1 = 0; n1 < *n_Phins_lz - 2 && log_Phins > Phins_bins_lz[n1 + 1]; n1++);
          n2 = n1 + 1;
          if(n2 >= *n_Phins_lz)
            n2 = *n_Phins_lz - 1;
          if(log_Phins >= Phins_bins_lz[0] && log_Phins <= Phins_bins_lz[*n_Phins_lz - 1])
            d_n_local = log_Phins - Phins_bins_lz[n1];
          else
            d_n_local = 0;
          delta_n = Phins_bins_lz[n2] - Phins_bins_lz[n1];
          if(delta_n > 0)
            d_n_local = d_n_local / delta_n;
          else
            d_n_local = 0;
        }
      else
        {
          n1 = 0;
          n2 = 0;
          d_n_local = 0.0;
        }
      d_n = d_n_local;

      /*search PhiT6 */
      if(log_PhiT6 > PhiT6_bins[0])
        {
          for(i61 = 0; i61 < *n_PhiT6 - 2 && log_PhiT6 > PhiT6_bins[i61 + 1]; i61++);
          i62 = i61 + 1;
          if(i62 >= *n_PhiT6)
            i62 = *n_PhiT6 - 1;
          if(log_PhiT6 >= PhiT6_bins[0] && log_PhiT6 <= PhiT6_bins[*n_PhiT6 - 1])
            d_i6_local = log_PhiT6 - PhiT6_bins[i61];
          else
            d_i6_local = 0;
          delta_i6 = PhiT6_bins[i62] - PhiT6_bins[i61];
          if(delta_i6 > 0)
            d_i6_local = d_i6_local / delta_i6;
          else
            d_i6_local = 0;
        }
      else
        {
          i61 = 0;
          i62 = 0;
          d_i6_local = 0.0;
        }
      d_i6 = d_i6_local;

      /*search PhiT7 */
      if(log_PhiT7 > PhiT7_bins[0])
        {
          for(i71 = 0; i71 < *n_PhiT7 - 2 && log_PhiT7 > PhiT7_bins[i71 + 1]; i71++);
          i72 = i71 + 1;
          if(i72 >= *n_PhiT7)
            i72 = *n_PhiT7 - 1;
          if(log_PhiT7 >= PhiT7_bins[0] && log_PhiT7 <= PhiT7_bins[*n_PhiT7 - 1])
            d_i7_local = log_PhiT7 - PhiT7_bins[i71];
          else
            d_i7_local = 0;
          delta_i7 = PhiT7_bins[i72] - PhiT7_bins[i71];
          if(delta_i7 > 0)
            d_i7_local = d_i7_local / delta_i7;
          else
            d_i7_local = 0;
        }
      else
        {
          i71 = 0;
          i72 = 0;
          d_i7_local = 0.0;
        }
      d_i7 = d_i7_local;

      if(log_PhiT8 > PhiT8_bins[0])
        {
          for(i81 = 0; i81 < *n_PhiT8 - 2 && log_PhiT8 > PhiT8_bins[i81 + 1]; i81++);
          i82 = i81 + 1;
          if(i82 >= *n_PhiT8)
            i82 = *n_PhiT8 - 1;
          if(log_PhiT8 >= PhiT8_bins[0] && log_PhiT8 <= PhiT8_bins[*n_PhiT8 - 1])
            d_i8_local = log_PhiT8 - PhiT8_bins[i81];
          else
            d_i8_local = 0;
          delta_i8 = PhiT8_bins[i82] - PhiT8_bins[i81];
          if(delta_i8 > 0)
            d_i8_local = d_i8_local / delta_i8;
          else
            d_i8_local = 0;
        }
      else
        {
          i81 = 0;
          i82 = 0;
          d_i8_local = 0.0;
        }
      d_i8 = d_i8_local;
    }
  else
    {
      /* search Phios */
      if(log_Phios > Phios_bins[0])
        {
          for(j1 = 0; j1 < *n_Phios - 2 && log_Phios > Phios_bins[j1 + 1]; j1++);
          j2 = j1 + 1;
          if(j2 >= *n_Phios)
            j2 = *n_Phios - 1;
          if(log_Phios >= Phios_bins[1] && log_Phios <= Phios_bins[*n_Phios - 1])
            d_j_local = log_Phios - Phios_bins[j1];
          else
            d_j_local = 0;
          delta_j = Phios_bins[j2] - Phios_bins[j1];
          if(delta_j > 0)
            d_j_local = d_j_local / delta_j;
          else
            d_j_local = 0;
        }
      else
        {
          j1 = 0;
          j2 = 0;
          d_j_local = 0.0;
        }
      d_j = d_j_local;

      /*search Phins */
      if(log_Phins > Phins_bins[0])
        {
          for(n1 = 0; n1 < *n_Phins - 2 && log_Phins > Phins_bins[n1 + 1]; n1++);
          n2 = n1 + 1;
          if(n2 >= *n_Phins)
            n2 = *n_Phins - 1;
          if(log_Phins >= Phins_bins[1] && log_Phins <= Phins_bins[*n_Phins - 1])
            d_n_local = log_Phins - Phins_bins[n1];
          else
            d_n_local = 0;
          delta_n = Phins_bins[n2] - Phins_bins[n1];
          if(delta_n > 0)
            d_n_local = d_n_local / delta_n;
          else
            d_n_local = 0;
        }
      else
        {
          n1 = 0;
          n2 = 0;
          d_n_local = 0.0;
        }
      d_n = d_n_local;
    }
#endif
#endif



#ifdef GFM_AGN_RADIATION
  /* search metallicity, and interpolate in linear space, not log space */
  if(log_MetallicityInSolar > metallicityInSolar_bins[0])
    {
      for(m1 = 0; m1 < *n_metallicityInSolar - 2 && log_MetallicityInSolar > metallicityInSolar_bins[m1 + 1]; m1++);
      m2 = m1 + 1;
      if(m2 >= *n_metallicityInSolar)
        m2 = *n_metallicityInSolar - 1;
      if(log_MetallicityInSolar >= metallicityInSolar_bins[0] && log_MetallicityInSolar <= metallicityInSolar_bins[*n_metallicityInSolar - 1])
        d_m_local = pow(10, log_MetallicityInSolar) - pow(10, metallicityInSolar_bins[m1]);
      else
        d_m_local = 0;
      delta_m = pow(10, metallicityInSolar_bins[m2]) - pow(10, metallicityInSolar_bins[m1]);
      if(delta_m > 0)
        d_m_local = d_m_local / delta_m;
      else
        d_m_local = 0;
    }
  else
    {
      m1 = 0;
      m2 = 0;
      d_m_local = 0.0;
    }
  d_m = d_m_local;
#else
  /* set metallicity bins, do not interpolate */
  m1 = 0;
  d_m = 0;
#endif

  /* search hydrogen number density */
  if(log_HydrogenNumberDensity > hydrogenNumberDensity_bins[0])
    {
      for(k1 = 0; k1 < *n_hydrogenNumberDensity - 2 && log_HydrogenNumberDensity > hydrogenNumberDensity_bins[k1 + 1]; k1++);
      k2 = k1 + 1;
      if(k2 >= *n_hydrogenNumberDensity)
        k2 = *n_hydrogenNumberDensity - 1;
      if(log_HydrogenNumberDensity >= hydrogenNumberDensity_bins[0] && log_HydrogenNumberDensity <= hydrogenNumberDensity_bins[*n_hydrogenNumberDensity - 1])
        d_k_local = log_HydrogenNumberDensity - hydrogenNumberDensity_bins[k1];
      else
        d_k_local = 0;
      delta_k = hydrogenNumberDensity_bins[k2] - hydrogenNumberDensity_bins[k1];
      if(delta_k > 0)
        d_k_local = d_k_local / delta_k;
      else
        d_k_local = 0;
    }
  else
    {
      k1 = 0;
      k2 = 0;
      d_k_local = 0.0;
    }
  d_k = d_k_local;


  /* search temperature */
  if(log_Temperature > temperature_bins[0])
    {
      for(l1 = 0; l1 < *n_temperature - 2 && log_Temperature > temperature_bins[l1 + 1]; l1++);
      l2 = l1 + 1;
      if(l2 >= *n_temperature)
        l2 = *n_temperature - 1;
      if(log_Temperature >= temperature_bins[0] && log_Temperature <= temperature_bins[*n_temperature - 1])
        d_l_local = log_Temperature - temperature_bins[l1];
      else
        d_l_local = 0;
      delta_l = temperature_bins[l2] - temperature_bins[l1];
      if(delta_l > 0)
        d_l_local = d_l_local / delta_l;
      else
        d_l_local = 0;
    }
  else
    {
      l1 = 0;
      l2 = 0;
      d_l_local = 0.0;
    }
  d_l = d_l_local;

  //printf("VALS: %g %g %g %g\n", All.cf_redshift, log_MetallicityInSolar, log_HydrogenNumberDensity, log_Temperature);
  //printf("BINS: %d %d %d    DELTA: %g %g %g\n", i1, k1, l1, d_i, d_k, d_l);
  //coolingMetalRate = netCoolingRate[i1][k1][l1];

#ifdef GFM_AGN_RADIATION
  coolingMetalRate = interpol_5d(netCoolingRate, i1, j1, m1, k1, l1, d_i, d_j, d_m, d_k, d_l);
#else
#if defined (RADCOOL) && !defined(RADCOOL_HOTHALO)
  cool = interpol_5d(netCoolingRate, i1, j1, n1, k1, l1, d_i, d_j, d_n, d_k, d_l);
  heat = interpol_5d(netHeatingRate, i1, j1, n1, k1, l1, d_i, d_j, d_n, d_k, d_l);
#else
#if defined (RADCOOL) && defined(RADCOOL_HOTHALO)
  if(All.cf_redshift < HOTHALO_ON_zFACTOR)
    {
      cool = interpol_8d(netCoolingRateLowZ, i1, j1, n1, i61, i71, i81, k1, l1, d_i, d_j, d_n, d_i6, d_i7, d_i8, d_k, d_l);
      heat = interpol_8d(netHeatingRateLowZ, i1, j1, n1, i61, i71, i81, k1, l1, d_i, d_j, d_n, d_i6, d_i7, d_i8, d_k, d_l);
    }
  else
    {
      cool = interpol_5d(netCoolingRate, i1, j1, n1, k1, l1, d_i, d_j, d_n, d_k, d_l);
      heat = interpol_5d(netHeatingRate, i1, j1, n1, k1, l1, d_i, d_j, d_n, d_k, d_l);
    }
#else
  coolingMetalRate = interpol_4d(netCoolingRate, i1, m1, k1, l1, d_i, d_m, d_k, d_l);
#endif
#endif
#endif

#ifdef RADCOOL
  coolingMetalRate = exp(cool) - exp(heat);
#endif

#if !defined(GFM_AGN_RADIATION)
  /* apply metallicity scaling */
  coolingMetalRate *= pow(10.0, log_MetallicityInSolar);
#endif

  return coolingMetalRate;
}

/* call this at beginning to initialize metal line cooling */
void init_cooling_metal()
{
  read_cooling_tables_init();
}


void dump_cooling_table()
{
  int i, k, l;
  char buf[255];

  sprintf(buf, "%s/%s", All.OutputDir, "coolingtable.txt");

  FILE *fd = fopen("coolingtable.txt", "w");
#ifdef GFM_AGN_RADIATION
  for(i = 0; i < *n_redshift; i++)
    for(k = 0; k < *n_hydrogenNumberDensity; k++)
      for(l = 0; l < *n_temperature; l++)
        fprintf(fd, "%d %d %d %.15e\n", i, k, l, netCoolingRate[i][0][0][k][l]);
#else
#ifdef RADCOOL
  for(i = 0; i < *n_redshift; i++)
    for(k = 0; k < *n_hydrogenNumberDensity; k++)
      for(l = 0; l < *n_temperature; l++)
        fprintf(fd, "%d %d %d %.15e\n", i, k, l, netCoolingRate[i][0][0][k][l]);
#else
  for(i = 0; i < *n_redshift; i++)
    for(k = 0; k < *n_hydrogenNumberDensity; k++)
      for(l = 0; l < *n_temperature; l++)
        fprintf(fd, "%d %d %d %.15e\n", i, k, l, netCoolingRate[i][0][k][l]);
#endif
#endif
  fclose(fd);
}

#endif
