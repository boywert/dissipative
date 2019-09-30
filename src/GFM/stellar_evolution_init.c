/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/GFM/stellar_evolution_init.c
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

#ifdef GFM_STELLAR_EVOLUTION

#ifndef HAVE_HDF5
#error "Need HAVE_HDF5 enabled to load metal cooling tables"
#endif

/*!
 * This routine allocates the yield table arrays.
 */
int allocate_yield_tables()
{
  int i, j;

  /* allocate SNIa table */
  /* allocate element name array */
  yieldsSNIa.ElementName = (char **) mymalloc("yieldsSNIa.ElementName", yieldsSNIa.N_ELEMENTS * sizeof(char *));

  for(i = 0; i < yieldsSNIa.N_ELEMENTS; i++)
    yieldsSNIa.ElementName[i] = (char *) mymalloc("yieldsSNIa.ElementName", GFM_EL_NAME_LENGTH * sizeof(char));

  /* allocate yield array */
  yieldsSNIa.Yield = (MyFloat *) mymalloc("yieldsSNIa.Yield", yieldsSNIa.N_ELEMENTS * sizeof(MyFloat));

  /* allocate SNII table */
  /* allocate element name array */
  yieldsSNII.ElementName = (char **) mymalloc("yieldsSNII.ElementName", yieldsSNII.N_ELEMENTS * sizeof(char *));

  for(i = 0; i < yieldsSNII.N_ELEMENTS; i++)
    yieldsSNII.ElementName[i] = (char *) mymalloc("yieldsSNII.ElementName", GFM_EL_NAME_LENGTH * sizeof(char));

  /* allocate mass array */
  yieldsSNII.Mass = (MyFloat *) mymalloc("yieldsSNII.Mass", yieldsSNII.N_MASS * sizeof(MyFloat));

  /* allocate metallicity array */
  yieldsSNII.Metallicity = (MyFloat *) mymalloc("yieldsSNII.Metallicity", yieldsSNII.N_Z * sizeof(MyFloat));

  /* allocate yield array */
  yieldsSNII.Yield = (MyFloat ***) mymalloc("yieldsSNII.Yield", yieldsSNII.N_Z * sizeof(MyFloat **));
  yieldsSNII.Ejecta = (MyFloat **) mymalloc("yieldsSNII.Ejecta", yieldsSNII.N_Z * sizeof(MyFloat *));
  yieldsSNII.TotalMetals = (MyFloat **) mymalloc("yieldsSNII.TotalMetals", yieldsSNII.N_Z * sizeof(MyFloat *));

  for(i = 0; i < yieldsSNII.N_Z; i++)
    {
      yieldsSNII.Yield[i] = (MyFloat **) mymalloc("yieldsSNII.Yield", yieldsSNII.N_ELEMENTS * sizeof(MyFloat *));
      yieldsSNII.Ejecta[i] = (MyFloat *) mymalloc("yieldsSNII.Ejecta", yieldsSNII.N_MASS * sizeof(MyFloat));
      yieldsSNII.TotalMetals[i] = (MyFloat *) mymalloc("yieldsSNII.TotalMetals", yieldsSNII.N_MASS * sizeof(MyFloat));

      for(j = 0; j < yieldsSNII.N_ELEMENTS; j++)
        yieldsSNII.Yield[i][j] = (MyFloat *) mymalloc("yieldsSNII.Yield", yieldsSNII.N_MASS * sizeof(MyFloat));
    }

  /* allocate AGB table */
  /* allocate element name array */
  yieldsAGB.ElementName = (char **) mymalloc("yieldsAGB.ElementName", yieldsAGB.N_ELEMENTS * sizeof(char *));

  for(i = 0; i < yieldsAGB.N_ELEMENTS; i++)
    yieldsAGB.ElementName[i] = (char *) mymalloc("yieldsAGB.ElementName", GFM_EL_NAME_LENGTH * sizeof(char));

  /* allocate mass array */
  yieldsAGB.Mass = (MyFloat *) mymalloc("yieldsAGB.Mass", yieldsAGB.N_MASS * sizeof(MyFloat));

  /* allocate metallicity array */
  yieldsAGB.Metallicity = (MyFloat *) mymalloc("yieldsAGB.Metallicity", yieldsAGB.N_Z * sizeof(MyFloat));

  /* allocate yield array */
  yieldsAGB.Yield = (MyFloat ***) mymalloc("yieldsAGB.Yield", yieldsAGB.N_Z * sizeof(MyFloat **));
  yieldsAGB.Ejecta = (MyFloat **) mymalloc("yieldsAGB.Ejecta", yieldsAGB.N_Z * sizeof(MyFloat *));
  yieldsAGB.TotalMetals = (MyFloat **) mymalloc("yieldsAGB.TotalMetals", yieldsAGB.N_Z * sizeof(MyFloat *));

  for(i = 0; i < yieldsAGB.N_Z; i++)
    {
      yieldsAGB.Yield[i] = (MyFloat **) mymalloc("yieldsAGB.Yield", yieldsAGB.N_ELEMENTS * sizeof(MyFloat *));
      yieldsAGB.Ejecta[i] = (MyFloat *) mymalloc("yieldsAGB.Ejecta", yieldsAGB.N_MASS * sizeof(MyFloat));
      yieldsAGB.TotalMetals[i] = (MyFloat *) mymalloc("yieldsAGB.TotalMetals", yieldsAGB.N_MASS * sizeof(MyFloat));

      for(j = 0; j < yieldsAGB.N_ELEMENTS; j++)
        yieldsAGB.Yield[i][j] = (MyFloat *) mymalloc("yieldsAGB.Yield", yieldsAGB.N_MASS * sizeof(MyFloat));
    }

  /* allocate Lifetime table */
  /* allocate mass array */
  Lifetimes.Mass = (MyFloat *) mymalloc("Lifetimes.Mass", Lifetimes.N_MASS * sizeof(MyFloat));

  /* allocate metallicity array */
  Lifetimes.Metallicity = (MyFloat *) mymalloc("Lifetimes.Metallicity", Lifetimes.N_Z * sizeof(MyFloat));

  /* allocate lifetime array */
  Lifetimes.Dyingtime = (MyFloat **) mymalloc("Lifetimes.Dyingtime", Lifetimes.N_Z * sizeof(MyFloat *));

  for(i = 0; i < Lifetimes.N_Z; i++)
    Lifetimes.Dyingtime[i] = (MyFloat *) mymalloc("Lifetimes.Dyingtime", Lifetimes.N_MASS * sizeof(MyFloat));


  return 0;
}

/*!
 * This routine reads in the yield tables (HDF5 format).
 */
int read_yield_tables()
{
  int i, j, k;

  char fname[MAXLEN_PATH], setname[MAXLEN_PATH];

  hid_t file_id, dataset, datatype;

  double *tmpdbl;

  /* read in table dimensions for SNIa */
  sprintf(fname, "%s/%s", All.YieldTablePath, GFM_SNIa_YIELDNAME);

  file_id = my_H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  dataset = my_H5Dopen(file_id, "Number_of_species");

  my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &yieldsSNIa.N_ELEMENTS, "Number_of_species");
  my_H5Dclose(dataset, "Number_of_species");

  my_H5Fclose(file_id, fname);

  mpi_printf("GFM_STELLAR_EVOLUTION: SNIa Header Complete...\n");

  /* read in table dimensions for SNII */
  sprintf(fname, "%s/%s", All.YieldTablePath, GFM_SNII_YIELDNAME);

  file_id = my_H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  dataset = my_H5Dopen(file_id, "Number_of_species");
  my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &yieldsSNII.N_ELEMENTS, "Number_of_species");
  my_H5Dclose(dataset, "Number_of_species");

  dataset = my_H5Dopen(file_id, "Number_of_masses");
  my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &yieldsSNII.N_MASS, "Number_of_masses");
  my_H5Dclose(dataset, "Number_of_masses");

  dataset = my_H5Dopen(file_id, "Number_of_metallicities");
  my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &yieldsSNII.N_Z, "Number_of_metallicities");
  my_H5Dclose(dataset, "Number_of_metallicities");

  my_H5Fclose(file_id, fname);

  mpi_printf("GFM_STELLAR_EVOLUTION: SNII Header Complete...\n");

  /* read in table dimensions for AGB */
  sprintf(fname, "%s/%s", All.YieldTablePath, GFM_AGB_YIELDNAME);

  file_id = my_H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  dataset = my_H5Dopen(file_id, "Number_of_species");
  my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &yieldsAGB.N_ELEMENTS, "Number_of_species");
  my_H5Dclose(dataset, "Number_of_species");

  dataset = my_H5Dopen(file_id, "Number_of_masses");
  my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &yieldsAGB.N_MASS, "Number_of_masses");
  my_H5Dclose(dataset, "Number_of_masses");

  dataset = my_H5Dopen(file_id, "Number_of_metallicities");
  my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &yieldsAGB.N_Z, "Number_of_metallicities");
  my_H5Dclose(dataset, "Number_of_metallicities");

  my_H5Fclose(file_id, fname);

  mpi_printf("GFM_STELLAR_EVOLUTION: AGB Header Complete...\n");

  /* read in table dimensions for Lifetime table */
  sprintf(fname, "%s/%s", All.YieldTablePath, GFM_LIFETIMENAME);

  file_id = my_H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  dataset = my_H5Dopen(file_id, "Number_of_masses");
  my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Lifetimes.N_MASS, "Number_of_masses");
  my_H5Dclose(dataset, "Number_of_masses");

  dataset = my_H5Dopen(file_id, "Number_of_metallicities");
  my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Lifetimes.N_Z, "Number_of_metallicities");
  my_H5Dclose(dataset, "Number_of_metallicities");

  my_H5Fclose(file_id, fname);

  /* allocate memory for yield tables (based on dimensions just read in) */
  allocate_yield_tables();

  /* read in SNIa yield table */
  sprintf(fname, "%s/%s", All.YieldTablePath, GFM_SNIa_YIELDNAME);

  file_id = my_H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  datatype = my_H5Tcopy(H5T_C_S1);
  my_H5Tset_size(datatype, H5T_VARIABLE);
  dataset = my_H5Dopen(file_id, "Species_names");
  my_H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, yieldsSNIa.ElementName, "Species_names");
  my_H5Dclose(dataset, "Species_names");
  my_H5Tclose(datatype);

  tmpdbl = (double *) mymalloc("yield", yieldsSNIa.N_ELEMENTS * sizeof(double));
  dataset = my_H5Dopen(file_id, "Yield");
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, "Yield");
  my_H5Dclose(dataset, "Yield");
  for(i = 0; i < yieldsSNIa.N_ELEMENTS; i++)
    yieldsSNIa.Yield[i] = tmpdbl[i];
  myfree(tmpdbl);

  tmpdbl = (double *) mymalloc("yield", sizeof(double));
  dataset = my_H5Dopen(file_id, "Total_Metals");
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, "Total_Metals");
  my_H5Dclose(dataset, "Total_Metals");
  yieldsSNIa.TotalMetals_spline = *tmpdbl;
  myfree(tmpdbl);

  my_H5Fclose(file_id, fname);

  mpi_printf("GFM_STELLAR_EVOLUTION: SNIa yields and total metals complete...\n");

  /* read in SNII mass and metallicity arrays and yield table */
  sprintf(fname, "%s/%s", All.YieldTablePath, GFM_SNII_YIELDNAME);

  file_id = my_H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  datatype = my_H5Tcopy(H5T_C_S1);
  my_H5Tset_size(datatype, H5T_VARIABLE);
  dataset = my_H5Dopen(file_id, "Species_names");
  my_H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, yieldsSNII.ElementName, "Species_names");
  my_H5Dclose(dataset, "Species_names");
  my_H5Tclose(datatype);

  tmpdbl = (double *) mymalloc("mass", yieldsSNII.N_MASS * sizeof(double));
  dataset = my_H5Dopen(file_id, "Masses");
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, "Masses");
  my_H5Dclose(dataset, "Masses");
  for(i = 0; i < yieldsSNII.N_MASS; i++)
    yieldsSNII.Mass[i] = tmpdbl[i];
  myfree(tmpdbl);

  tmpdbl = (double *) mymalloc("metallicity", yieldsSNII.N_Z * sizeof(double));
  dataset = my_H5Dopen(file_id, "Metallicities");
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, "Metallicities");
  my_H5Dclose(dataset, "Metallicities");
  for(i = 0; i < yieldsSNII.N_Z; i++)
    yieldsSNII.Metallicity[i] = tmpdbl[i];
  myfree(tmpdbl);

  double tempyield1[yieldsSNII.N_ELEMENTS][yieldsSNII.N_MASS];
  double tempej1[yieldsSNII.N_MASS], tempmet1[yieldsSNII.N_MASS];
  char *tempname1[yieldsSNII.N_Z];

  datatype = my_H5Tcopy(H5T_C_S1);
  my_H5Tset_size(datatype, H5T_VARIABLE);
  dataset = my_H5Dopen(file_id, "Yield_names");
  my_H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, tempname1, "Yield_names");
  my_H5Dclose(dataset, "Yield_names");
  my_H5Tclose(datatype);

  for(i = 0; i < yieldsSNII.N_Z; i++)
    {
      sprintf(setname, "/Yields/%s/Yield", tempname1[i]);
      dataset = my_H5Dopen(file_id, setname);
      my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tempyield1, setname);
      my_H5Dclose(dataset, setname);
      sprintf(setname, "/Yields/%s/Ejected_mass", tempname1[i]);
      dataset = my_H5Dopen(file_id, setname);
      my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tempej1, setname);
      my_H5Dclose(dataset, setname);
      sprintf(setname, "/Yields/%s/Total_Metals", tempname1[i]);
      dataset = my_H5Dopen(file_id, setname);
      my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tempmet1, setname);
      my_H5Dclose(dataset, setname);


      for(k = 0; k < yieldsSNII.N_MASS; k++)
        {
          yieldsSNII.Ejecta[i][k] = tempej1[k];
          yieldsSNII.TotalMetals[i][k] = tempmet1[k];
          for(j = 0; j < yieldsSNII.N_ELEMENTS; j++)
            yieldsSNII.Yield[i][j][k] = tempyield1[j][k];
        }
    }

  my_H5Fclose(file_id, fname);

  mpi_printf("GFM_STELLAR_EVOLUTION: SNII mass, metallicity, ejected mass, yields, and total metals complete...\n");

  /* read in AGB mass and metallicity arrays and yield table */
  sprintf(fname, "%s/%s", All.YieldTablePath, GFM_AGB_YIELDNAME);

  file_id = my_H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  /* allocate element name array */
  datatype = my_H5Tcopy(H5T_C_S1);
  my_H5Tset_size(datatype, H5T_VARIABLE);
  dataset = my_H5Dopen(file_id, "Species_names");
  my_H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, yieldsAGB.ElementName, "Species_names");
  my_H5Dclose(dataset, "Species_names");
  my_H5Tclose(datatype);

  tmpdbl = (double *) mymalloc("mass", yieldsAGB.N_MASS * sizeof(double));
  dataset = my_H5Dopen(file_id, "Masses");
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, "Masses");
  my_H5Dclose(dataset, "Masses");
  for(i = 0; i < yieldsAGB.N_MASS; i++)
    yieldsAGB.Mass[i] = tmpdbl[i];
  myfree(tmpdbl);

  tmpdbl = (double *) mymalloc("metallicity", yieldsAGB.N_Z * sizeof(double));
  dataset = my_H5Dopen(file_id, "Metallicities");
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, "Metallicities");
  my_H5Dclose(dataset, "Metallicities");
  for(i = 0; i < yieldsAGB.N_Z; i++)
    yieldsAGB.Metallicity[i] = tmpdbl[i];
  myfree(tmpdbl);

  double tempyield2[yieldsAGB.N_ELEMENTS][yieldsAGB.N_MASS];
  double tempej2[yieldsAGB.N_MASS], tempmet2[yieldsAGB.N_MASS];
  char *tempname2[yieldsAGB.N_Z];

  datatype = my_H5Tcopy(H5T_C_S1);
  my_H5Tset_size(datatype, H5T_VARIABLE);
  dataset = my_H5Dopen(file_id, "Yield_names");
  my_H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, tempname2, "Yield_names");
  my_H5Dclose(dataset, "Yield_names");
  my_H5Tclose(datatype);

  for(i = 0; i < yieldsAGB.N_Z; i++)
    {
      sprintf(setname, "/Yields/%s/Yield", tempname2[i]);
      dataset = my_H5Dopen(file_id, setname);
      my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tempyield2, setname);
      my_H5Dclose(dataset, setname);
      sprintf(setname, "/Yields/%s/Ejected_mass", tempname2[i]);
      dataset = my_H5Dopen(file_id, setname);
      my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tempej2, setname);
      my_H5Dclose(dataset, setname);
      sprintf(setname, "/Yields/%s/Total_Metals", tempname2[i]);
      dataset = my_H5Dopen(file_id, setname);
      my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tempmet2, setname);
      my_H5Dclose(dataset, setname);

      for(k = 0; k < yieldsAGB.N_MASS; k++)
        {
          yieldsAGB.Ejecta[i][k] = tempej2[k];
          yieldsAGB.TotalMetals[i][k] = tempmet2[k];

          for(j = 0; j < yieldsAGB.N_ELEMENTS; j++)
            yieldsAGB.Yield[i][j][k] = tempyield2[j][k];
        }
    }

  my_H5Fclose(file_id, fname);

  mpi_printf("GFM_STELLAR_EVOLUTION: AGB mass, metallicity, ejected mass, yields, and total metals complete...\n");

  /* read in Lifetime table */
  sprintf(fname, "%s/%s", All.YieldTablePath, GFM_LIFETIMENAME);

  file_id = my_H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  tmpdbl = (double *) mymalloc("mass", Lifetimes.N_MASS * sizeof(double));
  dataset = my_H5Dopen(file_id, "Masses");
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, "Masses");
  my_H5Dclose(dataset, "Masses");
  for(i = 0; i < Lifetimes.N_MASS; i++)
    Lifetimes.Mass[i] = tmpdbl[i];
  myfree(tmpdbl);

  tmpdbl = (double *) mymalloc("metallicity", Lifetimes.N_Z * sizeof(double));
  dataset = my_H5Dopen(file_id, "Metallicities");
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, "Metallicities");
  my_H5Dclose(dataset, "Metallicities");
  for(i = 0; i < Lifetimes.N_Z; i++)
    Lifetimes.Metallicity[i] = tmpdbl[i];
  myfree(tmpdbl);

  double temptime[Lifetimes.N_Z][Lifetimes.N_MASS];

  dataset = my_H5Dopen(file_id, "Lifetimes");
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temptime, "Lifetimes");
  my_H5Dclose(dataset, "Lifetimes");

  for(i = 0; i < Lifetimes.N_Z; i++)
    for(j = 0; j < Lifetimes.N_MASS; j++)
      Lifetimes.Dyingtime[i][j] = log10(temptime[i][j]);

  my_H5Fclose(file_id, fname);



  /* get these indices and store them to avoid look-up overhead in inner loops */
  element_index_Hydrogen = element_index("Hydrogen");
  element_index_Helium = element_index("Helium");
  element_index_Carbon = element_index("Carbon");
  element_index_Nitrogen = element_index("Nitrogen");
  element_index_Oxygen = element_index("Oxygen");
  element_index_Neon = element_index("Neon");
  element_index_Magnesium = element_index("Magnesium");
  element_index_Silicon = element_index("Silicon");
  element_index_Iron = element_index("Iron");

  return 0;
}



void spline_interpolate_yields(int spline_element_index)
{
  char element_name[GFM_EL_NAME_LENGTH];

  int i, k;

  double *yield, *mass, result;

  int element_index;

  gsl_interp_accel *my_accel_ptr;

  gsl_spline *my_spline_ptr;

  sprintf(element_name, "%s", ElementNames[spline_element_index]);      /* metals, including Helium */


  if(ThisTask == 0)
    printf("GFM_STELLAR_EVOLUTION: Computing yield for %-10s \t index=%02d\n", element_name, spline_element_index);

  /*  SNIa yields */
  element_index = get_element_index(yieldsSNIa.ElementName, yieldsSNIa.N_ELEMENTS, element_name);

  if(element_index < 0)
    terminate("GFM_STELLAR_EVOLUTION: SNIa: element not found %s\n", element_name);

  yieldsSNIa.spline[spline_element_index] = yieldsSNIa.Yield[element_index];

  /* SNII yields */
  element_index = get_element_index(yieldsSNII.ElementName, yieldsSNII.N_ELEMENTS, element_name);

  if(element_index < 0)
    terminate("GFM_STELLAR_EVOLUTION: SNII: element not found %s\n", element_name);

  yield = (double *) mymalloc("yield", yieldsSNII.N_MASS * sizeof(double));
  mass = (double *) mymalloc("SNIImass", yieldsSNII.N_MASS * sizeof(double));

  my_accel_ptr = gsl_interp_accel_alloc();
  my_spline_ptr = gsl_spline_alloc(gsl_interp_linear, yieldsSNII.N_MASS);

  for(i = 0; i < yieldsSNII.N_Z; i++)
    {
      for(k = 0; k < yieldsSNII.N_MASS; k++)
        {
          yield[k] = yieldsSNII.Yield[i][element_index][k] / pow(10.0, yieldsSNII.Mass[k]);
          mass[k] = yieldsSNII.Mass[k];
        }

      gsl_spline_init(my_spline_ptr, mass, yield, yieldsSNII.N_MASS);

      for(k = 0; k < GFM_N_MASS_BINS; k++)
        {
          if(yield_mass_bin[k] < mass[0])
            result = yield[0];
          else if(yield_mass_bin[k] > mass[yieldsSNII.N_MASS - 1])
            result = yield[yieldsSNII.N_MASS - 1];
          else
            result = gsl_spline_eval(my_spline_ptr, yield_mass_bin[k], my_accel_ptr);

          yieldsSNII.spline[i][spline_element_index][k] = pow(10.0, yield_mass_bin[k]) * result;
        }

    }

  myfree(mass);
  myfree(yield);
  gsl_spline_free(my_spline_ptr);
  gsl_interp_accel_free(my_accel_ptr);

  /* AGB yields */
  element_index = get_element_index(yieldsAGB.ElementName, yieldsAGB.N_ELEMENTS, element_name);

  if(element_index < 0)
    {
      if(ThisTask == 0)
        printf("GFM_STELLAR_EVOLUTION: AGB: element not found: %s\n", element_name);

      for(i = 0; i < yieldsAGB.N_Z; i++)
        {
          for(k = 0; k < GFM_N_MASS_BINS; k++)
            yieldsAGB.spline[i][spline_element_index][k] = 0.0;
        }
    }
  else
    {
      yield = (double *) mymalloc("yield", yieldsAGB.N_MASS * sizeof(double));
      mass = (double *) mymalloc("AGBmass", yieldsAGB.N_MASS * sizeof(double));

      my_accel_ptr = gsl_interp_accel_alloc();
      my_spline_ptr = gsl_spline_alloc(gsl_interp_linear, yieldsAGB.N_MASS);

      for(i = 0; i < yieldsAGB.N_Z; i++)
        {
          for(k = 0; k < yieldsAGB.N_MASS; k++)
            {
              yield[k] = yieldsAGB.Yield[i][element_index][k] / pow(10.0, yieldsAGB.Mass[k]);
              mass[k] = yieldsAGB.Mass[k];
            }

          gsl_spline_init(my_spline_ptr, mass, yield, yieldsAGB.N_MASS);

          for(k = 0; k < GFM_N_MASS_BINS; k++)
            {
              if(yield_mass_bin[k] < mass[0])
                result = yield[0];
              else if(yield_mass_bin[k] > mass[yieldsAGB.N_MASS - 1])
                result = yield[yieldsAGB.N_MASS - 1];
              else
                result = gsl_spline_eval(my_spline_ptr, yield_mass_bin[k], my_accel_ptr);

              yieldsAGB.spline[i][spline_element_index][k] = pow(10.0, yield_mass_bin[k]) * result;
            }

        }

      myfree(mass);
      myfree(yield);
      gsl_spline_free(my_spline_ptr);
      gsl_interp_accel_free(my_accel_ptr);
    }
}


#ifdef GFM_NORMALIZED_METAL_ADVECTION
void spline_augment_yields(void)
{
  int spline_element_index = GFM_N_CHEM_ELEMENTS - 1;   /* these are the other metals */

  mpi_printf("GFM_STELLAR_EVOLUTION: Augmenting yield for %-10s \t index=%02d\n", ElementNames[spline_element_index], spline_element_index);

  /* SNIa yields */
  yieldsSNIa.spline[spline_element_index] = yieldsSNIa.TotalMetals_spline;
  for(int i = 2; i < GFM_N_CHEM_ELEMENTS - 1; i++)
    yieldsSNIa.spline[spline_element_index] -= yieldsSNIa.spline[i];


  /* SNII yields */
  for(int i = 0; i < yieldsSNII.N_Z; i++)
    for(int k = 0; k < GFM_N_MASS_BINS; k++)
      {
        double sum_yields = 0;

        for(int j = 0; j < GFM_N_CHEM_ELEMENTS - 1; j++)
          sum_yields += yieldsSNII.spline[i][j][k];

        yieldsSNII.spline[i][spline_element_index][k] = -sum_yields;
      }


  /* AGB yields */
  for(int i = 0; i < yieldsAGB.N_Z; i++)
    for(int k = 0; k < GFM_N_MASS_BINS; k++)
      {
        double sum_yields = 0;

        for(int j = 0; j < GFM_N_CHEM_ELEMENTS - 1; j++)
          sum_yields += yieldsAGB.spline[i][j][k];

        yieldsAGB.spline[i][spline_element_index][k] = -sum_yields;
      }

}



#endif



void spline_interpolate_ejecta()
{
  int i, k;

  double *yield, *mass, result;

  gsl_interp_accel *my_accel_ptr;

  gsl_spline *my_spline_ptr;


  mpi_printf("GFM_STELLAR_EVOLUTION: Computing ejecta\n");

  /* SNII yields */
  yield = (double *) mymalloc("yield", yieldsSNII.N_MASS * sizeof(double));
  mass = (double *) mymalloc("SNIImass", yieldsSNII.N_MASS * sizeof(double));

  my_accel_ptr = gsl_interp_accel_alloc();
  my_spline_ptr = gsl_spline_alloc(gsl_interp_linear, yieldsSNII.N_MASS);

  for(i = 0; i < yieldsSNII.N_Z; i++)
    {
      for(k = 0; k < yieldsSNII.N_MASS; k++)
        {
          yield[k] = yieldsSNII.Ejecta[i][k] / pow(10.0, yieldsSNII.Mass[k]);
          mass[k] = yieldsSNII.Mass[k];
        }
      gsl_spline_init(my_spline_ptr, mass, yield, yieldsSNII.N_MASS);

      for(k = 0; k < GFM_N_MASS_BINS; k++)
        {
          if(yield_mass_bin[k] < mass[0])
            result = yield[0];
          else if(yield_mass_bin[k] > mass[yieldsSNII.N_MASS - 1])
            result = yield[yieldsSNII.N_MASS - 1];
          else
            result = gsl_spline_eval(my_spline_ptr, yield_mass_bin[k], my_accel_ptr);

          yieldsSNII.Ejecta_spline[i][k] = pow(10.0, yield_mass_bin[k]) * result;
        }
    }

  for(i = 0; i < yieldsSNII.N_Z; i++)
    {
      for(k = 0; k < yieldsSNII.N_MASS; k++)
        {
          yield[k] = yieldsSNII.TotalMetals[i][k] / pow(10.0, yieldsSNII.Mass[k]);
          mass[k] = yieldsSNII.Mass[k];
        }

      gsl_spline_init(my_spline_ptr, mass, yield, yieldsSNII.N_MASS);

      for(k = 0; k < GFM_N_MASS_BINS; k++)
        {
          if(yield_mass_bin[k] < mass[0])
            result = yield[0];
          else if(yield_mass_bin[k] > mass[yieldsSNII.N_MASS - 1])
            result = yield[yieldsSNII.N_MASS - 1];
          else
            result = gsl_spline_eval(my_spline_ptr, yield_mass_bin[k], my_accel_ptr);

          yieldsSNII.TotalMetals_spline[i][k] = pow(10.0, yield_mass_bin[k]) * result;
        }
    }

  myfree(mass);
  myfree(yield);
  gsl_spline_free(my_spline_ptr);
  gsl_interp_accel_free(my_accel_ptr);

  /* AGB yields */
  yield = (double *) mymalloc("yield", yieldsAGB.N_MASS * sizeof(double));
  mass = (double *) mymalloc("massAGB", yieldsAGB.N_MASS * sizeof(double));

  my_accel_ptr = gsl_interp_accel_alloc();
  my_spline_ptr = gsl_spline_alloc(gsl_interp_linear, yieldsAGB.N_MASS);

  for(i = 0; i < yieldsAGB.N_Z; i++)
    {
      for(k = 0; k < yieldsAGB.N_MASS; k++)
        {
          yield[k] = yieldsAGB.Ejecta[i][k] / pow(10.0, yieldsAGB.Mass[k]);
          mass[k] = yieldsAGB.Mass[k];
        }

      gsl_spline_init(my_spline_ptr, mass, yield, yieldsAGB.N_MASS);

      for(k = 0; k < GFM_N_MASS_BINS; k++)
        {
          if(yield_mass_bin[k] < mass[0])
            result = yield[0];
          else if(yield_mass_bin[k] > mass[yieldsAGB.N_MASS - 1])
            result = yield[yieldsAGB.N_MASS - 1];
          else
            result = gsl_spline_eval(my_spline_ptr, yield_mass_bin[k], my_accel_ptr);

          yieldsAGB.Ejecta_spline[i][k] = pow(10.0, yield_mass_bin[k]) * result;
        }
    }

  for(i = 0; i < yieldsAGB.N_Z; i++)
    {
      for(k = 0; k < yieldsAGB.N_MASS; k++)
        {
          yield[k] = yieldsAGB.TotalMetals[i][k] / pow(10.0, yieldsAGB.Mass[k]);
          mass[k] = yieldsAGB.Mass[k];
        }

      gsl_spline_init(my_spline_ptr, mass, yield, yieldsAGB.N_MASS);

      for(k = 0; k < GFM_N_MASS_BINS; k++)
        {
          if(yield_mass_bin[k] < mass[0])
            result = yield[0];
          else if(yield_mass_bin[k] > mass[yieldsAGB.N_MASS - 1])
            result = yield[yieldsAGB.N_MASS - 1];
          else
            result = gsl_spline_eval(my_spline_ptr, yield_mass_bin[k], my_accel_ptr);

          yieldsAGB.TotalMetals_spline[i][k] = pow(10.0, yield_mass_bin[k]) * result;
        }
    }

  myfree(mass);
  myfree(yield);
  gsl_spline_free(my_spline_ptr);
  gsl_interp_accel_free(my_accel_ptr);
}



void init_yields(void)
{
  MyFloat lm_min, lm_max, dlm;

  int i, j, k;

  double total_mass;

  /* let all tasks read the yield tables */
  read_yield_tables();


  /* set up stellar mass bins */
  lm_min = log10(All.IMF_MinMass_Msun); /* min mass in solar masses */
  lm_max = log10(All.IMF_MaxMass_Msun); /* max mass in solar masses */

  dlm = (lm_max - lm_min) / (double) (GFM_N_MASS_BINS - 1);

  for(i = 0; i < GFM_N_MASS_BINS; i++)
    yield_mass_bin[i] = dlm * i + lm_min;


  /* convert tables to log10 */
  /* SNII tables */
  for(k = 0; k < yieldsSNII.N_MASS; k++)
    yieldsSNII.Mass[k] = log10(yieldsSNII.Mass[k]);

  for(i = 0; i < yieldsSNII.N_Z; i++)
    {
      if(yieldsSNII.Metallicity[i] > 0)
        yieldsSNII.Metallicity[i] = log10(yieldsSNII.Metallicity[i]);
      else
        yieldsSNII.Metallicity[i] = GFM_MIN_METAL;
    }

  /* AGB tables */
  for(i = 0; i < yieldsAGB.N_MASS; i++)
    yieldsAGB.Mass[i] = log10(yieldsAGB.Mass[i]);

  for(i = 0; i < yieldsAGB.N_Z; i++)
    {
      if(yieldsAGB.Metallicity[i] > 0)
        yieldsAGB.Metallicity[i] = log10(yieldsAGB.Metallicity[i]);
      else
        yieldsAGB.Metallicity[i] = GFM_MIN_METAL;
    }


  /* allocate spline interpolated tables */
  /* SNIa tables */
  yieldsSNIa.spline = (MyFloat *) mymalloc("yieldsSNIa.spline", GFM_N_CHEM_ELEMENTS * sizeof(MyFloat));

  /* SNII tables */
  yieldsSNII.spline = (MyFloat ***) mymalloc("yieldsSNII.spline", yieldsSNII.N_Z * sizeof(MyFloat **));
  yieldsSNII.Ejecta_spline = (MyFloat **) mymalloc("yieldsSNII.Ejecta_spline", yieldsSNII.N_Z * sizeof(MyFloat *));
  yieldsSNII.TotalMetals_spline = (MyFloat **) mymalloc("yieldsSNII.TotalMetals_spline", yieldsSNII.N_Z * sizeof(MyFloat *));
  for(i = 0; i < yieldsSNII.N_Z; i++)
    {
      yieldsSNII.spline[i] = (MyFloat **) mymalloc("yieldsSNII.spline", GFM_N_CHEM_ELEMENTS * sizeof(MyFloat *));
      yieldsSNII.Ejecta_spline[i] = (MyFloat *) mymalloc("yieldsSNII.Ejecta_spline", GFM_N_MASS_BINS * sizeof(MyFloat));
      yieldsSNII.TotalMetals_spline[i] = (MyFloat *) mymalloc("yieldsSNII.TotalMetals_spline", GFM_N_MASS_BINS * sizeof(MyFloat));
      for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
        yieldsSNII.spline[i][j] = (MyFloat *) mymalloc("yieldsSNII.spline", GFM_N_MASS_BINS * sizeof(MyFloat));
    }

  /* AGB tables */
  yieldsAGB.spline = (MyFloat ***) mymalloc("yieldsAGB.spline", yieldsAGB.N_Z * sizeof(MyFloat **));
  yieldsAGB.Ejecta_spline = (MyFloat **) mymalloc("yieldsAGB.Ejecta_spline", yieldsAGB.N_Z * sizeof(MyFloat *));
  yieldsAGB.TotalMetals_spline = (MyFloat **) mymalloc("yieldsAGB.TotalMetals_spline", yieldsAGB.N_Z * sizeof(MyFloat *));
  for(i = 0; i < yieldsAGB.N_Z; i++)
    {
      yieldsAGB.Ejecta_spline[i] = (MyFloat *) mymalloc("yieldsAGB.Ejecta_spline", GFM_N_MASS_BINS * sizeof(MyFloat));
      yieldsAGB.TotalMetals_spline[i] = (MyFloat *) mymalloc("yieldsAGB.TotalMetals_spline", GFM_N_MASS_BINS * sizeof(MyFloat));
      yieldsAGB.spline[i] = (MyFloat **) mymalloc("yieldsAGB.spline", GFM_N_CHEM_ELEMENTS * sizeof(MyFloat *));
      for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
        yieldsAGB.spline[i][j] = (MyFloat *) mymalloc("yieldsAGB.spline", GFM_N_MASS_BINS * sizeof(MyFloat));
    }

  /* for each element retable the yields */
  for(i = 0; i < GFM_N_CHEM_ELEMENTS; i++)
    {
#ifdef GFM_NORMALIZED_METAL_ADVECTION
      if(i == GFM_N_CHEM_ELEMENTS - 1)
        continue;
#endif
      spline_interpolate_yields(i);
    }

#ifdef GFM_NORMALIZED_METAL_ADVECTION
  spline_augment_yields();
#endif

  /* retable the ejecta and total metal mass released */
  spline_interpolate_ejecta();

  /* double check scaling */
  for(i = 0; i < yieldsSNII.N_Z; i++)
    {
      for(k = 0; k < GFM_N_MASS_BINS; k++)
        {

          total_mass = yieldsSNII.TotalMetals_spline[i][k] + yieldsSNII.Metallicity[i] * yieldsSNII.Ejecta_spline[i][k];
          if(element_index("Hydrogen") > 0)
            total_mass += yieldsSNII.spline[i][element_index("Hydrogen")][k] + (1.0 - yieldsSNII.Metallicity[i]) * 0.75 * yieldsSNII.Ejecta_spline[i][k];
          if(element_index("Helium") > 0)
            total_mass += yieldsSNII.spline[i][element_index("Helium")][k] + (1.0 - yieldsSNII.Metallicity[i]) * 0.25 * yieldsSNII.Ejecta_spline[i][k];


          if(total_mass > pow(10.0, yield_mass_bin[k]))
            terminate("GFM_STELLAR_EVOLUTION: Not conserving mass! Your type II scale factors are probably too high");
        }
    }
}

#endif
