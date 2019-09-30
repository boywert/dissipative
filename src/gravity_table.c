/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/gravity_table.c
 * \date        MM/YYYY
 * \author      Robin Tress
 * \brief       Gravitational forces from an external look-up table
 * \details
 *
 *This module allows to read in a look-up table storing the gravitational accelerations on a R-z grid;
 *by calling the function grav_table_find_grav_acceleration() at the right position (in the grav_external.c file)
 *it is so possible to compute the acceleration at a generic position for every particle.
 * 
 *For usage set the config flag GRAVITY_TABLE and define the location of the file containing the look-up table 
 *as a parameter ExternalGravForcesFile in the parameter file. Then just make sure to call the function 
 *grav_table_find_grav_acceleration() where needed. 
 *
 *Pay particular attention when building the gravitational acceleration look-up table and storing it to a file:
 *  - the grid has to be equi-spaced both in z and in R 
 *  - the number of z bins have to be equal to the number of R bins (square table)
 *  - the file has to be formatted as an ASCII table in the form: 
 *       R_0 z_0 aR_00 az_00
 *       R_0 z_1 aR_01 az_01
 *       ...
 *       R_1 z_0 aR_10 az_10
 *       R_1 z_1 aR_11 az_11
 *       ...
 *  - R and z have to be given in kpc
 *  - aR and az have to be given in cm/s^2
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

/*! \brief reads in the gravity table and initializes the module 
 *
 * TODO: check for the integrity of the table just read in 
 */
void grav_table_init(void)
{
  
  FILE *gravity_table_file = fopen(All.ExternalGravForcesFile,"r"); 
  int ch;
  int i;
  
  if(gravity_table_file == NULL)
    terminate("GRAVITY_TABLE: error opening the gravitational forces file\n");

  All.NGravityTableBins = 0;

  while ( !feof(gravity_table_file))
    {
      ch = fgetc(gravity_table_file);
      if(ch == '\n')
        All.NGravityTableBins++;
    }

  GravT = (struct grav_table_data *) mymalloc_movable(&GravT, "GravT", All.NGravityTableBins*sizeof(struct grav_table_data));
  mpi_printf("GRAVITY_TABLE: number of bins = %i \n", All.NGravityTableBins);
  rewind(gravity_table_file);

/* R and z are only needed for security checks */
  for(i = 0; i < All.NGravityTableBins; i++)
    {
      fscanf(gravity_table_file, "%lf%lf%lf%lf", &GravT[i].R, &GravT[i].z, &GravT[i].acc_R, &GravT[i].acc_z);
      GravT[i].R *= KILOPARSEC / All.UnitLength_in_cm; 
      GravT[i].z *= KILOPARSEC / All.UnitLength_in_cm;
      GravT[i].acc_R *= All.UnitTime_in_s * (All.UnitTime_in_s / All.UnitLength_in_cm);
      GravT[i].acc_z *= All.UnitTime_in_s * (All.UnitTime_in_s / All.UnitLength_in_cm);

      if(i == 0)
        {
          All.MaxR = GravT[i].R;
          All.MaxZ = GravT[i].z;
          All.MinR = GravT[i].R;
          All.MinZ = GravT[i].z;
        }
      else
        {   
          All.MaxR = ((GravT[i].R > All.MaxR) ? GravT[i].R : All.MaxR); 
          All.MaxZ = ((GravT[i].z > All.MaxZ) ? GravT[i].z : All.MaxZ);

          All.MinR = ((GravT[i].R < All.MinR) ? GravT[i].R : All.MinR);
          All.MinZ = ((GravT[i].z < All.MinZ) ? GravT[i].z : All.MinZ);
        }

    }

  fclose(gravity_table_file);
  
  All.DeltaR = (All.MaxR - All.MinR) / (sqrt(All.NGravityTableBins) - 1.);
  All.DeltaZ = (All.MaxZ - All.MinZ) / (sqrt(All.NGravityTableBins) - 1.);
}

/*! \brief The gravitational acceleration is defined on a R-z grid,
 *         this function calculates its value for a generic position
 *         by bilinear interpolation.
 * \param xi,yi,zi Position of the point for which we need the gravitational acceleration
 * \param sp_acc pointer to the memory where the gravitational acceleration has to be saved to 
 */
void grav_table_find_grav_acceleration(double xi, double yi, double zi, double *sp_acc)
{
  double sp_acc_R1, sp_acc_R2, sp_acc_R, sp_acc_z1, sp_acc_z2;

  int j, k; /*index of the radial and z bin within the gravity table*/

  double cyl_Ri = sqrt(xi * xi + yi * yi);
  double abzi   = fabs(zi);

  /*if out of gravity table boundary then assume the last valid point*/

  if(cyl_Ri < All.MinR)
    cyl_Ri = All.MinR;

  if(abzi < All.MinZ)
    abzi = All.MinZ;
 
  if(cyl_Ri > All.MaxR)
    {
      cyl_Ri = All. MaxR;
      j = floor((cyl_Ri - All.MinR) / All.DeltaR) - 1;
    }
  else
    j = floor((cyl_Ri - All.MinR) / All.DeltaR);

  if(abzi > All.MaxZ)
    {
      abzi = All.MaxZ;
      k = floor((abzi   - All.MinZ) / All.DeltaZ) - 1;
    }
  else
    k = floor((abzi   - All.MinZ) / All.DeltaZ);

  double z2 = (k+1) * All.DeltaZ + All.MinZ;
  double r2 = (j+1) * All.DeltaR + All.MinR;
  double z1 = z2 - All.DeltaZ;
  double r1 = r2 - All.DeltaR;

  int nz = sqrt(All.NGravityTableBins); 

  sp_acc_R1 = grav_table_interpolation(abzi,   z1, z2, GravT[j*nz + k].acc_R,     GravT[j*nz + k+1].acc_R);
  sp_acc_R2 = grav_table_interpolation(abzi,   z1, z2, GravT[(j+1)*nz + k].acc_R, GravT[(j+1)*nz + k+1].acc_R); 

  sp_acc_R  = grav_table_interpolation(cyl_Ri, r1, r2, sp_acc_R1, sp_acc_R2);
  
  double theta = atan2(yi,xi);
  
  sp_acc[0] = sp_acc_R * cos(theta);
  sp_acc[1] = sp_acc_R * sin(theta);
  
  sp_acc_z1 = grav_table_interpolation(abzi, z1, z2, GravT[j*nz + k].acc_z,     GravT[j*nz + k+1].acc_z); 
  sp_acc_z2 = grav_table_interpolation(abzi, z1, z2, GravT[(j+1)*nz + k].acc_z, GravT[(j+1)*nz + k+1].acc_z); 
  
  sp_acc[2] = ((zi < 0) ? -1. : 1.) * grav_table_interpolation(cyl_Ri, r1, r2, sp_acc_z1, sp_acc_z2);

  /*Exceptions*/
  if (GravT[j*nz + k].R     > cyl_Ri || GravT[(j+1)*nz + k].R   < cyl_Ri ||
      GravT[j*nz + k+1].R   > cyl_Ri || GravT[(j+1)*nz + k+1].R < cyl_Ri ||
      GravT[j*nz + k].z     > abzi   || GravT[j*nz + k+1].z     < abzi   ||
      GravT[(j+1)*nz + k].z > abzi   || GravT[(j+1)*nz + k+1].z < abzi   ||
      ((j+1)*nz + k+1)     > All.NGravityTableBins)
    terminate("Corrupt gravity table\n");    

}

double grav_table_interpolation(double x, double x1, double x2, double f_x1, double f_x2)
{
  return (f_x1 + (f_x2 - f_x1) * (x - x1) / (x2 - x1));
}

