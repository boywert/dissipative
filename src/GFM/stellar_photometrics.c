/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/GFM/stellar_photometrics.c
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
#include <gsl/gsl_math.h>

#include "../allvars.h"
#include "../proto.h"

#ifdef GFM_STELLAR_PHOTOMETRICS

#ifndef HAVE_HDF5
#error "Need HAVE_HDF5 enabled to load stellar photometrics tables"
#endif

int allocate_stellar_photometrics_tables(void)
{
  int i;

  stellarLuminosityTable.LogMetallicity_bins = (MyFloat *) mymalloc("stellarLuminosityTable.LogMetallicity_bins", stellarLuminosityTable.N_LogMetallicity * sizeof(MyFloat));

  stellarLuminosityTable.LogAgeInGyr_bins = (MyFloat *) mymalloc("stellarLuminosityTable.LogAgeInGyr_bins", stellarLuminosityTable.N_LogAgeInGyr * sizeof(MyFloat));

  stellarLuminosityTable.Magnitude_U = (MyFloat **) mymalloc("stellarLuminosityTable.Magnitude_U", stellarLuminosityTable.N_LogMetallicity * sizeof(MyFloat **));
  stellarLuminosityTable.Magnitude_B = (MyFloat **) mymalloc("stellarLuminosityTable.Magnitude_B", stellarLuminosityTable.N_LogMetallicity * sizeof(MyFloat **));
  stellarLuminosityTable.Magnitude_V = (MyFloat **) mymalloc("stellarLuminosityTable.Magnitude_V", stellarLuminosityTable.N_LogMetallicity * sizeof(MyFloat **));
  stellarLuminosityTable.Magnitude_K = (MyFloat **) mymalloc("stellarLuminosityTable.Magnitude_K", stellarLuminosityTable.N_LogMetallicity * sizeof(MyFloat **));
  stellarLuminosityTable.Magnitude_g = (MyFloat **) mymalloc("stellarLuminosityTable.Magnitude_g", stellarLuminosityTable.N_LogMetallicity * sizeof(MyFloat **));
  stellarLuminosityTable.Magnitude_r = (MyFloat **) mymalloc("stellarLuminosityTable.Magnitude_r", stellarLuminosityTable.N_LogMetallicity * sizeof(MyFloat **));
  stellarLuminosityTable.Magnitude_i = (MyFloat **) mymalloc("stellarLuminosityTable.Magnitude_i", stellarLuminosityTable.N_LogMetallicity * sizeof(MyFloat **));
  stellarLuminosityTable.Magnitude_z = (MyFloat **) mymalloc("stellarLuminosityTable.Magnitude_z", stellarLuminosityTable.N_LogMetallicity * sizeof(MyFloat **));

  for(i = 0; i < stellarLuminosityTable.N_LogMetallicity; i++)
    {
      stellarLuminosityTable.Magnitude_U[i] = (MyFloat *) mymalloc("stellarLuminosityTable.Magnitude_U", stellarLuminosityTable.N_LogAgeInGyr * sizeof(MyFloat *));
      stellarLuminosityTable.Magnitude_B[i] = (MyFloat *) mymalloc("stellarLuminosityTable.Magnitude_B", stellarLuminosityTable.N_LogAgeInGyr * sizeof(MyFloat *));
      stellarLuminosityTable.Magnitude_V[i] = (MyFloat *) mymalloc("stellarLuminosityTable.Magnitude_V", stellarLuminosityTable.N_LogAgeInGyr * sizeof(MyFloat *));
      stellarLuminosityTable.Magnitude_K[i] = (MyFloat *) mymalloc("stellarLuminosityTable.Magnitude_K", stellarLuminosityTable.N_LogAgeInGyr * sizeof(MyFloat *));
      stellarLuminosityTable.Magnitude_g[i] = (MyFloat *) mymalloc("stellarLuminosityTable.Magnitude_g", stellarLuminosityTable.N_LogAgeInGyr * sizeof(MyFloat *));
      stellarLuminosityTable.Magnitude_r[i] = (MyFloat *) mymalloc("stellarLuminosityTable.Magnitude_r", stellarLuminosityTable.N_LogAgeInGyr * sizeof(MyFloat *));
      stellarLuminosityTable.Magnitude_i[i] = (MyFloat *) mymalloc("stellarLuminosityTable.Magnitude_i", stellarLuminosityTable.N_LogAgeInGyr * sizeof(MyFloat *));
      stellarLuminosityTable.Magnitude_z[i] = (MyFloat *) mymalloc("stellarLuminosityTable.Magnitude_z", stellarLuminosityTable.N_LogAgeInGyr * sizeof(MyFloat *));

    }

  return 0;
}


int read_stellar_photometrics_tables()
{
  int i, k, l;
  char fname[MAXLEN_PATH];
  hid_t file_id, dataset;
  double *tmpdbl;

  mpi_printf("GFM_STELLAR_PHOTOMETRICS: opening file...\n");
  sprintf(fname, "%s/%s", All.PhotometricsTablePath, "stellar_photometrics.hdf5");
  mpi_printf("GFM_STELLAR_PHOTOMETRICS: done.\n");

  file_id = my_H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  dummyvarstellarphotometrics = 0;

  dataset = my_H5Dopen(file_id, "N_LogMetallicity");
  my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &stellarLuminosityTable.N_LogMetallicity, "N_LogMetallicity");
  my_H5Dclose(dataset, "N_LogMetallicity");

  dataset = my_H5Dopen(file_id, "N_LogAgeInGyr");
  my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &stellarLuminosityTable.N_LogAgeInGyr, "N_LogAgeInGyr");
  my_H5Dclose(dataset, "N_LogAgeInGyr");

  if(ThisTask == 0)
    {
      printf("GFM_STELLAR_PHOTOMETRICS: found %d LogMetallicity entries\n", stellarLuminosityTable.N_LogMetallicity);
      printf("GFM_STELLAR_PHOTOMETRICS: found %d LogAgeInGyr entries\n", stellarLuminosityTable.N_LogAgeInGyr);
    }
  my_H5Fclose(file_id, fname);


  allocate_stellar_photometrics_tables();

  mpi_printf("GFM_STELLAR_PHOTOMETRICS: allocated stellar photometrics arrays...\n");

  file_id = my_H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  tmpdbl = (double *) mymalloc("tmpdbl", stellarLuminosityTable.N_LogMetallicity * sizeof(double));
  dataset = my_H5Dopen(file_id, "LogMetallicity_bins");
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, "LogMetallicity_bins");
  my_H5Dclose(dataset, "LogMetallicity_bins");
  for(i = 0; i < stellarLuminosityTable.N_LogMetallicity; i++)
    stellarLuminosityTable.LogMetallicity_bins[i] = tmpdbl[i];
  myfree(tmpdbl);

  tmpdbl = (double *) mymalloc("tmpdbl", stellarLuminosityTable.N_LogAgeInGyr * sizeof(double));
  dataset = my_H5Dopen(file_id, "LogAgeInGyr_bins");
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, "LogAgeInGyr_bins");
  my_H5Dclose(dataset, "LogAgeInGyr_bins");
  for(i = 0; i < stellarLuminosityTable.N_LogAgeInGyr; i++)
    stellarLuminosityTable.LogAgeInGyr_bins[i] = tmpdbl[i];
  myfree(tmpdbl);


  tmpdbl = (double *) mymalloc("tmpdbl", stellarLuminosityTable.N_LogMetallicity * stellarLuminosityTable.N_LogAgeInGyr * sizeof(double));

  dataset = my_H5Dopen(file_id, "Magnitude_U");
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, "Magnitude_U");
  my_H5Dclose(dataset, "Magnitude_U");

  for(k = 0; k < stellarLuminosityTable.N_LogMetallicity; k++)
    for(l = 0; l < stellarLuminosityTable.N_LogAgeInGyr; l++)
      {
        int index = l + stellarLuminosityTable.N_LogAgeInGyr * k;
        stellarLuminosityTable.Magnitude_U[k][l] = tmpdbl[index];
      }


  dataset = my_H5Dopen(file_id, "Magnitude_B");
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, "Magnitude_B");
  my_H5Dclose(dataset, "Magnitude_B");

  for(k = 0; k < stellarLuminosityTable.N_LogMetallicity; k++)
    for(l = 0; l < stellarLuminosityTable.N_LogAgeInGyr; l++)
      {
        int index = l + stellarLuminosityTable.N_LogAgeInGyr * k;
        stellarLuminosityTable.Magnitude_B[k][l] = tmpdbl[index];
      }


  dataset = my_H5Dopen(file_id, "Magnitude_V");
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, "Magnitude_V");
  my_H5Dclose(dataset, "Magnitude_V");

  for(k = 0; k < stellarLuminosityTable.N_LogMetallicity; k++)
    for(l = 0; l < stellarLuminosityTable.N_LogAgeInGyr; l++)
      {
        int index = l + stellarLuminosityTable.N_LogAgeInGyr * k;
        stellarLuminosityTable.Magnitude_V[k][l] = tmpdbl[index];
      }


  dataset = my_H5Dopen(file_id, "Magnitude_K");
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, "Magnitude_K");
  my_H5Dclose(dataset, "Magnitude_K");

  for(k = 0; k < stellarLuminosityTable.N_LogMetallicity; k++)
    for(l = 0; l < stellarLuminosityTable.N_LogAgeInGyr; l++)
      {
        int index = l + stellarLuminosityTable.N_LogAgeInGyr * k;
        stellarLuminosityTable.Magnitude_K[k][l] = tmpdbl[index];
      }


  dataset = my_H5Dopen(file_id, "Magnitude_g");
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, "Magnitude_g");
  my_H5Dclose(dataset, "Magnitude_g");

  for(k = 0; k < stellarLuminosityTable.N_LogMetallicity; k++)
    for(l = 0; l < stellarLuminosityTable.N_LogAgeInGyr; l++)
      {
        int index = l + stellarLuminosityTable.N_LogAgeInGyr * k;
        stellarLuminosityTable.Magnitude_g[k][l] = tmpdbl[index];
      }


  dataset = my_H5Dopen(file_id, "Magnitude_r");
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, "Magnitude_r");
  my_H5Dclose(dataset, "Magnitude_r");

  for(k = 0; k < stellarLuminosityTable.N_LogMetallicity; k++)
    for(l = 0; l < stellarLuminosityTable.N_LogAgeInGyr; l++)
      {
        int index = l + stellarLuminosityTable.N_LogAgeInGyr * k;
        stellarLuminosityTable.Magnitude_r[k][l] = tmpdbl[index];
      }


  dataset = my_H5Dopen(file_id, "Magnitude_i");
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, "Magnitude_i");
  my_H5Dclose(dataset, "Magnitude_i");

  for(k = 0; k < stellarLuminosityTable.N_LogMetallicity; k++)
    for(l = 0; l < stellarLuminosityTable.N_LogAgeInGyr; l++)
      {
        int index = l + stellarLuminosityTable.N_LogAgeInGyr * k;
        stellarLuminosityTable.Magnitude_i[k][l] = tmpdbl[index];
      }


  dataset = my_H5Dopen(file_id, "Magnitude_z");
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, "Magnitude_z");
  my_H5Dclose(dataset, "Magnitude_z");

  for(k = 0; k < stellarLuminosityTable.N_LogMetallicity; k++)
    for(l = 0; l < stellarLuminosityTable.N_LogAgeInGyr; l++)
      {
        int index = l + stellarLuminosityTable.N_LogAgeInGyr * k;
        stellarLuminosityTable.Magnitude_z[k][l] = tmpdbl[index];
      }


  myfree(tmpdbl);
  my_H5Fclose(file_id, fname);

  mpi_printf("GFM_STELLAR_PHOTOMETRICS: stellar photometrics arrays read in...\n");

#if (DVR_RENDER==0)
  return 0;
#endif

  if(ThisTask == 0)
    {
      char buf[1000];
      sprintf(buf, "%s/stellar_photometrics_bins.txt", All.OutputDir);
      FILE *fd = fopen(buf, "w");
      fprintf(fd, "#LogMetallicity bins:\n");
      for(i = 0; i < stellarLuminosityTable.N_LogMetallicity; i++)
        fprintf(fd, "#%d %g\n", i, stellarLuminosityTable.LogMetallicity_bins[i]);

      fprintf(fd, "#LogAgeInGyr bins:\n");
      for(i = 0; i < stellarLuminosityTable.N_LogAgeInGyr; i++)
        fprintf(fd, "#%d %g\n", i, stellarLuminosityTable.LogAgeInGyr_bins[i]);

      fclose(fd);

      //test_stellar_photometrics(1.0, 0.02, -3.5, 1.2, 100);
    }

  return 0;
}



void get_stellar_photometrics(MyFloat ageInGyr, MyFloat metallicity, MyFloat mass, stellar_photometrics * st_photo)
{
  double delta_i, d_i_local, delta_k, d_k_local;
  double d_i, d_k;
  int i1, i2, k1, k2;

  double log_Metallicity;
  double log_AgeInGyr = log10(ageInGyr);
  /* SSP mass in units of M_sun */
  double mass_in_solar = mass * All.UnitMass_in_g / All.HubbleParam / SOLAR_MASS;

  if(metallicity < 0)
    log_Metallicity = GFM_MIN_METAL;
  else
    log_Metallicity = log10(metallicity);

  /* search stellar age */
  if(log_AgeInGyr > stellarLuminosityTable.LogAgeInGyr_bins[0])
    {
      for(i1 = 0; i1 < stellarLuminosityTable.N_LogAgeInGyr - 2 && log_AgeInGyr > stellarLuminosityTable.LogAgeInGyr_bins[i1 + 1]; i1++);

      i2 = i1 + 1;

      if(i2 >= stellarLuminosityTable.N_LogAgeInGyr)
        i2 = stellarLuminosityTable.N_LogAgeInGyr - 1;

      if(log_AgeInGyr >= stellarLuminosityTable.LogAgeInGyr_bins[0] && log_AgeInGyr <= stellarLuminosityTable.LogAgeInGyr_bins[stellarLuminosityTable.N_LogAgeInGyr - 1])
        d_i_local = log_AgeInGyr - stellarLuminosityTable.LogAgeInGyr_bins[i1];
      else
        d_i_local = 0;

      delta_i = stellarLuminosityTable.LogAgeInGyr_bins[i2] - stellarLuminosityTable.LogAgeInGyr_bins[i1];

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


  /* search metallicity */
  if(log_Metallicity > stellarLuminosityTable.LogMetallicity_bins[0])
    {
      for(k1 = 0; k1 < stellarLuminosityTable.N_LogMetallicity - 2 && log_Metallicity > stellarLuminosityTable.LogMetallicity_bins[k1 + 1]; k1++);

      k2 = k1 + 1;

      if(k2 >= stellarLuminosityTable.N_LogMetallicity)
        k2 = stellarLuminosityTable.N_LogMetallicity - 1;

      if(log_Metallicity >= stellarLuminosityTable.LogMetallicity_bins[0] && log_Metallicity <= stellarLuminosityTable.LogMetallicity_bins[stellarLuminosityTable.N_LogMetallicity - 1])
        d_k_local = log_Metallicity - stellarLuminosityTable.LogMetallicity_bins[k1];
      else
        d_k_local = 0;

      delta_k = stellarLuminosityTable.LogMetallicity_bins[k2] - stellarLuminosityTable.LogMetallicity_bins[k1];

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

  /* absolute band magnitudes */
  st_photo->Magnitude_U = interpol_2d(stellarLuminosityTable.Magnitude_U, k1, i1, d_k, d_i) - 2.5 * log10(mass_in_solar);
  st_photo->Magnitude_B = interpol_2d(stellarLuminosityTable.Magnitude_B, k1, i1, d_k, d_i) - 2.5 * log10(mass_in_solar);
  st_photo->Magnitude_V = interpol_2d(stellarLuminosityTable.Magnitude_V, k1, i1, d_k, d_i) - 2.5 * log10(mass_in_solar);
  st_photo->Magnitude_K = interpol_2d(stellarLuminosityTable.Magnitude_K, k1, i1, d_k, d_i) - 2.5 * log10(mass_in_solar);
  st_photo->Magnitude_g = interpol_2d(stellarLuminosityTable.Magnitude_g, k1, i1, d_k, d_i) - 2.5 * log10(mass_in_solar);
  st_photo->Magnitude_r = interpol_2d(stellarLuminosityTable.Magnitude_r, k1, i1, d_k, d_i) - 2.5 * log10(mass_in_solar);
  st_photo->Magnitude_i = interpol_2d(stellarLuminosityTable.Magnitude_i, k1, i1, d_k, d_i) - 2.5 * log10(mass_in_solar);
  st_photo->Magnitude_z = interpol_2d(stellarLuminosityTable.Magnitude_z, k1, i1, d_k, d_i) - 2.5 * log10(mass_in_solar);

}

void init_stellar_photometrics(void)
{
  read_stellar_photometrics_tables();
#ifdef SUBFIND
  generate_random_directions();
#endif
}


void assign_stellar_photometrics(int i, stellar_photometrics * st_photo)
{
  MyFloat ageInGyr, metallicity, mass;

  if(P[i].Type == 4 && P[i].Mass > 0 && STP(i).BirthTime > 0)
    {
      /* BC03 tabulates in initial stellar mass, i.e. to account for mass loss we need to look up STP(i).InitialMass */
      mass = STP(i).InitialMass;
      metallicity = STP(i).Metallicity;
      ageInGyr = get_time_difference_in_Gyr(STP(i).BirthTime, All.Time);
      get_stellar_photometrics(ageInGyr, metallicity, mass, st_photo);
    }
  else
    memset(st_photo, 0, GFM_STELLAR_PHOTOMETRICS_BANDS * sizeof(MyFloat));
}

void test_stellar_photometrics(double mass_in_solar, double metallicity, double log_AgeInGyr_min, double log_AgeInGyr_max, int log_AgeInGyr_bins)
{
  MyFloat log_AgeInGyr, dlog_AgeInGyr = (log_AgeInGyr_max - log_AgeInGyr_min) / log_AgeInGyr_bins;
  stellar_photometrics st_photo;
  double mass_internal;
  int i;
  FILE *fd = fopen("stellar_photometrics_sample.txt", "w");

  mass_internal = mass_in_solar * SOLAR_MASS * All.HubbleParam / All.UnitMass_in_g;

  fprintf(fd, "#stellar photometrics for: Z=%g\n", metallicity);
  fprintf(fd, "#log_AgeInGyr   Magnitude_U   Magnitude_B   Magnitude_V   Magnitude_K   Magnitude_g   Magnitude_r   Magnitude_i   Magnitude_z\n");
  for(i = 0; i < log_AgeInGyr_bins; i++)
    {
      log_AgeInGyr = log_AgeInGyr_min + i * dlog_AgeInGyr;
      get_stellar_photometrics(pow(10.0, log_AgeInGyr), metallicity, mass_internal, &st_photo);
      fprintf(fd, "%g %g %g %g %g %g %g %g %g\n", log_AgeInGyr, st_photo.Magnitude_U, st_photo.Magnitude_B, st_photo.Magnitude_V,
              st_photo.Magnitude_K, st_photo.Magnitude_g, st_photo.Magnitude_r, st_photo.Magnitude_i, st_photo.Magnitude_z);
    }
  fclose(fd);
}

void generate_random_directions(void)
{
  if(ThisTask == 0)
    for(int i = 0; i < GFM_STELLAR_PHOTOMETRICS_DIRECTIONS; i++)
      {
        StellarPhotometricsRandomAngles[i][0] = acos(2 * get_random_number() - 1);
        StellarPhotometricsRandomAngles[i][1] = 2 * M_PI * get_random_number();
        // mpi_printf("GFM_STELLAR_PHOTOMETRICS: StellarPhotometricsRandomAngles[%d]=[%g %g]\n", i, StellarPhotometricsRandomAngles[i][0], StellarPhotometricsRandomAngles[i][1]);
      }

  /* to get the same directions on all tasks */
  MPI_Bcast(&StellarPhotometricsRandomAngles[0][0], GFM_STELLAR_PHOTOMETRICS_DIRECTIONS * 2, MPI_FLOAT, 0, MPI_COMM_WORLD);
}

#endif
