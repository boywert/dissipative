/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dust_live/dust_init.c
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
#include "../proto.h"
#include "../allvars.h"

#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#ifdef DUST_LIVE

#ifdef GFM_DUST
#error "DUST_LIVE is not compatible with GFM_DUST!  Turn one off!"
#endif

#ifdef DL_GRAIN_BINS

void init_dust(void)
{
  init_grain_size_distribution();

  mpi_printf("DUST_LIVE: Grain size distribution initialized.\n");
}

#if defined(DL_SHATTERING) || defined(DL_COAGULATION)
void load_grain_velocity_data(char *label, double *gsd_velocities)
{
  int N_pts;
  double *radius, *velocity;
  hid_t file_id, dataset;

  char setname[MAXLEN_PATH];

  file_id = my_H5Fopen(All.GrainDataPath, H5F_ACC_RDONLY, H5P_DEFAULT);

  sprintf(setname, "/%s/NPoints", label);
  dataset = my_H5Dopen(file_id, setname);
  my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &N_pts, setname);
  my_H5Dclose(dataset, setname);

  radius = (double *) mymalloc("radius", N_pts * sizeof(double));
  sprintf(setname, "/%s/Radius", label);
  dataset = my_H5Dopen(file_id, setname);
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, radius, setname);
  my_H5Dclose(dataset, setname);

  velocity = (double *) mymalloc("velocity", N_pts * sizeof(double));
  sprintf(setname, "/%s/Velocity", label);
  dataset = my_H5Dopen(file_id, setname);
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, velocity, setname);
  my_H5Dclose(dataset, setname);

  if(radius[0] > GSD.Midpoints[0])
    {
      terminate("DUST_LIVE: Grain velocity table does not extend to small enough radii.\n");
    }
  if(radius[N_pts-1] < GSD.Midpoints[DL_GRAIN_BINS-1])
    {
      terminate("DUST_LIVE: Grain velocity table does not extend to large enough radii.\n");
    }
  for(int i = 0, j = 0; i < DL_GRAIN_BINS; i++)
    {
      while(radius[j+1] < GSD.Midpoints[i])
        {
          j++;
        }

      /* Interpolate in log-space, since that's how the data is presented. */
      double x = log10(GSD.Midpoints[i]) - log10(radius[j]);
      double y = log10(radius[j+1]) - log10(GSD.Midpoints[i]);

      gsd_velocities[i] = (y / (x+y)) * velocity[j] + (x / (x+y)) * velocity[j+1];
    }

  myfree(velocity);
  myfree(radius);

  my_H5Fclose(file_id, All.GrainDataPath);
}
#endif

#ifdef DL_SNE_DESTRUCTION
void load_sn_destruction_data(void)
{
  int N_pts;
  double *arra;

  hid_t file_id, dataset;
  char setname[MAXLEN_PATH];

  file_id = my_H5Fopen(All.GrainDataPath, H5F_ACC_RDONLY, H5P_DEFAULT);

  sprintf(setname, "/SNDestruction/NRadius");
  dataset = my_H5Dopen(file_id, setname);
  my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &N_pts, setname);
  my_H5Dclose(dataset, setname);

  arra = (double *) mymalloc("arra", N_pts * sizeof(double));
  sprintf(setname, "/SNDestruction/Radius");
  dataset = my_H5Dopen(file_id, setname);
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, arra, setname);
  my_H5Dclose(dataset, setname);

  double arrxi[N_pts][N_pts];
  sprintf(setname, "/SNDestruction/XiFrac");
  dataset = my_H5Dopen(file_id, setname);
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, arrxi, setname);
  my_H5Dclose(dataset, setname);

  /* Array whose j-th entry is the index into GSD.Midpoints giving the midpoint
   * closest to arra[j]. */
  int table_to_bin[N_pts];
  for(int j = 0; j < N_pts; j++)
    {
      double min_dist = MAX_REAL_NUMBER;
      int min_ind;
      for(int jj = 0; jj < DL_GRAIN_BINS; jj++)
        {
          double dist = fabs(log10(GSD.Midpoints[jj]) - log10(arra[j]));
          if(dist < min_dist)
            {
              min_dist = dist;
              min_ind = jj;
            }
        }

      table_to_bin[j] = min_ind;
    }

  for(int i = 0; i < DL_GRAIN_BINS; i++)
    for(int j = 0; j < DL_GRAIN_BINS; j++)
      GSD.XiFrac[i][j] = 0.0;

  /* Because linear interpolation of these xi(a, a') fractions could get noisy,
   * let's just interpolate by taking each DL_GRAIN_BINS midpoint and assigning
   * it xi(a, a') values using the closest a' (start grain size) value. */
  for(int j = 0; j < DL_GRAIN_BINS; j++)
    {
      double min_dist = MAX_REAL_NUMBER;
      int min_ind;
      for(int jj = 0; jj < N_pts; jj++)
        {
          double dist = fabs(log10(GSD.Midpoints[j]) - log10(arra[jj]));
          if(dist < min_dist)
            {
              min_dist = dist;
              min_ind = jj;
            }
        }

      /* Index min_ind gives the entry in arra that is the closest starting
       * grain size to bin j.  Since the number of entries in arra may differ
       * from DL_GRAIN_BINS, we need to allocate the xi fractions in column
       * min_ind accordingly.  We assign each xi fraction in this column to its
       * closest DL_GRAIN_BIN. */
      for(int ii = 0; ii < N_pts; ii++)
        {
          GSD.XiFrac[table_to_bin[ii]][j] += arrxi[ii][min_ind] / GSD.Widths[table_to_bin[ii]];
        }
    }

  myfree(arra);

  my_H5Fclose(file_id, All.GrainDataPath);
}
#endif

#ifdef DL_PRODUCTION
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
void normalize_dnda(double *num_grains, double *bin_slopes)
#else
void normalize_dnda(double *num_grains)
#endif
{
  double tot_mass = 0.0;
  for(int i = 0; i < DL_GRAIN_BINS; i++)
    {
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
      tot_mass += num_grains[i] * bin_avg_mass(i, num_grains[i], bin_slopes[i]);
#else
      tot_mass += num_grains[i] * GSD.AvgMasses[i];
#endif
    }

  for(int i = 0; i < DL_GRAIN_BINS; i++)
    {
      num_grains[i] /= tot_mass;
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
      bin_slopes[i] /= tot_mass;
#endif
    }
}

double agb_dnda(double a)
{
  double a_AGB = 0.1; /* um */
  double sigma_AGB = 0.47;
  double lnfac = log(a / a_AGB);
  return 1.0/pow(a, 5) * exp(-lnfac*lnfac / (2.0*sigma_AGB*sigma_AGB));
}

void load_snii_dnda(void)
{
  int N_pts;
  double *loga, *logdnda;

  hid_t file_id, dataset;
  char setname[MAXLEN_PATH];

  file_id = my_H5Fopen(All.GrainDataPath, H5F_ACC_RDONLY, H5P_DEFAULT);

  sprintf(setname, "/InitialSNII/NPoints");
  dataset = my_H5Dopen(file_id, setname);
  my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &N_pts, setname);
  my_H5Dclose(dataset, setname);

  loga = (double *) mymalloc("loga", N_pts * sizeof(double));
  sprintf(setname, "/InitialSNII/Radius");
  dataset = my_H5Dopen(file_id, setname);
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, loga, setname);
  my_H5Dclose(dataset, setname);

  logdnda = (double *) mymalloc("logdnda", N_pts * sizeof(double));
  sprintf(setname, "/InitialSNII/GSD");
  dataset = my_H5Dopen(file_id, setname);
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, logdnda, setname);
  my_H5Dclose(dataset, setname);

  for(int i = 0; i < N_pts; i++)
    {
      /* Take logs so we can interpolate in log space. */
      loga[i] = log10(loga[i]);
      logdnda[i] = log10(logdnda[i]);
    }

  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear, N_pts);

  gsl_spline_init(spline, loga, logdnda, N_pts);

  for(int i = 0; i < DL_GRAIN_BINS; i++)
    {
      double dnda_lhs = gsl_spline_eval(spline, log10(GSD.Edges[i]), acc);
      double dnda_rhs = gsl_spline_eval(spline, log10(GSD.Edges[i+1]), acc);
      dnda_lhs = pow(10.0, dnda_lhs);
      dnda_rhs = pow(10.0, dnda_rhs);

      GSD.SNII_NumGrains[i] = (dnda_lhs + dnda_rhs) / 2.0 * GSD.Widths[i];
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
      GSD.SNII_BinSlopes[i] = (dnda_rhs - dnda_lhs) / GSD.Widths[i];
#endif
    }

  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);

  myfree(logdnda);
  myfree(loga);

  my_H5Fclose(file_id, All.GrainDataPath);
}

void spline_interpolate_condeff(int elem)
{
  gsl_interp_accel *my_accel_ptr = gsl_interp_accel_alloc();
  gsl_spline *my_spline_ptr = gsl_spline_alloc(gsl_interp_linear, GSD.AGB_CondEff_NM);
  double *mass = (double *) mymalloc("mass", GSD.AGB_CondEff_NM * sizeof(double));
  double *condeff = (double *) mymalloc("condeff", GSD.AGB_CondEff_NM * sizeof(double));
  double result;

  /* If element can condense into dust, have nonzero condensation efficiencies.
   * Otherwise, all condensation efficiencies are zero. */
  if((elem == element_index_Carbon) || (elem == element_index_Oxygen) || (elem == element_index_Magnesium) || (elem == element_index_Silicon) || (elem == element_index_Iron))
    {
      for(int i = 0; i < GSD.AGB_CondEff_NZ; i++)
        {
          for(int k = 0; k < GSD.AGB_CondEff_NM; k++)
            {
              mass[k] = GSD.AGB_CondEff_Mass[k]; /* already in log units */
              if(GSD.AGB_CondEff[i][elem][k] > 0.0)
                condeff[k] = log10(GSD.AGB_CondEff[i][elem][k]);
              else
                condeff[k] = -10.0;
            }

          gsl_spline_init(my_spline_ptr, mass, condeff, GSD.AGB_CondEff_NM);

          for(int k = 0; k < GFM_N_MASS_BINS; k++)
            {
              if(imf_mass_bin_log10[k] < mass[0])
                result = condeff[0];
              else if(imf_mass_bin_log10[k] > mass[GSD.AGB_CondEff_NM-1])
                result = condeff[GSD.AGB_CondEff_NM-1];
              else
                result = gsl_spline_eval(my_spline_ptr, imf_mass_bin_log10[k], my_accel_ptr);

              /* Cannot form more dust than overall metals. */
              GSD.AGB_CondEff_spline[i][elem][k] = dmin(1.0, pow(10.0, result));
            }
        }
    }
  else
    {
      for(int i = 0; i < GSD.AGB_CondEff_NZ; i++)
        {
          for(int k = 0; k < GFM_N_MASS_BINS; k++)
            {
              GSD.AGB_CondEff_spline[i][elem][k] = 0.0;
            }
        }
    }

    myfree(condeff);
    myfree(mass);
    gsl_spline_free(my_spline_ptr);
    gsl_interp_accel_free(my_accel_ptr);
}

/* Load dust condensation efficiencies for SNe and AGB stars.  Call after IMF
 * has been initialized so we know the masses for interpolation.
 *
 * Currently only AGB condensation efficiencies are functions of metallicity.
 * This HDF5 condensation efficiency table requires these values at the same
 * metallicity values used for AGB yields. */
void load_condensation_efficiencies(void)
{
  hid_t file_id, dataset;

  char setname[MAXLEN_PATH];
  double frac;

  file_id = my_H5Fopen(All.GrainDataPath, H5F_ACC_RDONLY, H5P_DEFAULT);

  /* SNII */
  for(int i = 0; i < GFM_N_CHEM_ELEMENTS; i++)
    {
      sprintf(setname, "/CondEffSNII/%s/CondEff", ElementNames[i]);
      dataset = my_H5Dopen(file_id, setname);
      my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &frac, setname);
      my_H5Dclose(dataset, setname);

      if((frac > 1.0) || (frac < 0.0))
        terminate("Invalid SNII condensation efficiency %g for element %s!", frac, ElementNames[i]);
      GSD.SNII_CondEff[i] = frac;
    }

  /* AGB */
  sprintf(setname, "/CondEffAGB/NMetallicities");
  dataset = my_H5Dopen(file_id, setname);
  my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &GSD.AGB_CondEff_NZ, setname);
  my_H5Dclose(dataset, setname);

  if(GSD.AGB_CondEff_NZ != yieldsAGB.N_Z)
    {
      terminate("Dust condensation efficiencies for AGB stars need to be tabulated at the same metallicities used for AGB yields!");
    }

  sprintf(setname, "/CondEffAGB/NMasses");
  dataset = my_H5Dopen(file_id, setname);
  my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &GSD.AGB_CondEff_NM, setname);
  my_H5Dclose(dataset, setname);

  /* Allocate space based on dimensions of grain data tables. */
  GSD.AGB_CondEff_Mass = (double *) mymalloc("AGB_CondEff_Mass", GSD.AGB_CondEff_NM * sizeof(double));
  GSD.AGB_CondEff_Metallicity = (double *) mymalloc("AGB_CondEff_Metallicity", GSD.AGB_CondEff_NZ * sizeof(double));

  GSD.AGB_CondEff = (double ***) mymalloc("AGB_CondEff", GSD.AGB_CondEff_NZ * sizeof(double **));
  for(int i = 0; i < GSD.AGB_CondEff_NZ; i++)
    {
      GSD.AGB_CondEff[i] = (double **) mymalloc("AGB_CondEff", GFM_N_CHEM_ELEMENTS * sizeof(double *));
      for(int j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
        {
          GSD.AGB_CondEff[i][j] = (double *) mymalloc("AGB_CondEff", GSD.AGB_CondEff_NM * sizeof(double));
        }
    }

  /* Read in the tables. */
  sprintf(setname, "/CondEffAGB/Metallicities");
  dataset = my_H5Dopen(file_id, setname);
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, GSD.AGB_CondEff_Metallicity, setname);
  my_H5Dclose(dataset, setname);

  sprintf(setname, "/CondEffAGB/Masses");
  dataset = my_H5Dopen(file_id, setname);
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, GSD.AGB_CondEff_Mass, setname);
  my_H5Dclose(dataset, setname);

  for(int i = 0; i < GSD.AGB_CondEff_NZ; i++)
    {
      for(int j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
        {
          sprintf(setname, "/CondEffAGB/Z%d/%s/CondEff", i, ElementNames[i]);
          dataset = my_H5Dopen(file_id, setname);
          my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, GSD.AGB_CondEff[i][j], "AGB_CondEff");
          my_H5Dclose(dataset, "AGB_CondEff");
        }
    }

  /* Allocate space for interpolated condensation efficiencies using GFM mass bins. */
  GSD.AGB_CondEff_spline = (double ***) mymalloc("AGB_CondEff_spline", GSD.AGB_CondEff_NZ * sizeof(double **));
  for(int i = 0; i < GSD.AGB_CondEff_NZ; i++)
    {
      GSD.AGB_CondEff_spline[i] = (double **) mymalloc("AGB_CondEff_spline", GFM_N_CHEM_ELEMENTS * sizeof(double *));
      for(int j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
        {
          GSD.AGB_CondEff_spline[i][j] = (double *) mymalloc("AGB_CondEff_spline", GFM_N_MASS_BINS * sizeof(double));
        }
    }

  /* To interpolate, use the logarithms of stellar mass and metallicity. */
  for(int i = 0; i < GSD.AGB_CondEff_NZ; i++)
    {
      if(GSD.AGB_CondEff_Metallicity[i] > 0.0)
        GSD.AGB_CondEff_Metallicity[i] = log10(GSD.AGB_CondEff_Metallicity[i]);
      else
        GSD.AGB_CondEff_Metallicity[i] = GFM_MIN_METAL;
    }

  for(int i = 0; i < GSD.AGB_CondEff_NM; i++)
    {
      GSD.AGB_CondEff_Mass[i] = log10(GSD.AGB_CondEff_Mass[i]);
    }

  for(int i = 0; i < GFM_N_CHEM_ELEMENTS; i++)
    {
      spline_interpolate_condeff(i);
    }

  for(int i = 0; i < GSD.AGB_CondEff_NZ; i++)
    {
      for(int j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
        {
          for(int k = 0; k < GFM_N_MASS_BINS; k++)
            {
              if((GSD.AGB_CondEff_spline[i][j][k] > 1.0) || (GSD.AGB_CondEff_spline[i][j][k] < 0.0))
                terminate("Invalid AGB condensation efficiency %g for element %s i=%d j=%d k=%d!", GSD.AGB_CondEff_spline[i][j][k], ElementNames[j], i, j, k);
            }
        }
    }

  my_H5Fclose(file_id, All.GrainDataPath);
}
#endif

void init_grain_size_distribution(void)
{
  GSD.Edges = (double *) mymalloc("GSD.Edges", (DL_GRAIN_BINS+1) * sizeof(double));
  GSD.Midpoints = (double *) mymalloc("GSD.Midpoints", DL_GRAIN_BINS * sizeof(double));
  GSD.Widths = (double *) mymalloc("GSD.Widths", DL_GRAIN_BINS * sizeof(double));
  GSD.AvgMasses = (double *) mymalloc("GSD.AvgMasses", DL_GRAIN_BINS * sizeof(double));
  GSD.EdgeMasses = (double *) mymalloc("GSD.EdgeMasses", (DL_GRAIN_BINS+1) * sizeof(double));
#if defined(DL_SHATTERING) || defined(DL_COAGULATION)
  GSD.VelocitiesCNM = (double *) mymalloc("GSD.VelocitiesCNM", DL_GRAIN_BINS * sizeof(double));
  GSD.VelocitiesWIM = (double *) mymalloc("GSD.VelocitiesWIM", DL_GRAIN_BINS * sizeof(double));
#endif

  double cm_per_um = 1.0e-4;
  GSD.InternalDensity = (All.GrainDensity * pow(cm_per_um, 3.0)) / (All.UnitMass_in_g / All.HubbleParam);

  double logdelta = (log10(All.MaxGrainSize) - log10(All.MinGrainSize)) / DL_GRAIN_BINS;
  double delta = pow(10.0, logdelta);

  GSD.Edges[0] = All.MinGrainSize;
  for(int i = 1; i <= DL_GRAIN_BINS; i++)
    {
      GSD.Edges[i] = GSD.Edges[i-1] * delta;
    }

  for(int i = 0; i <= DL_GRAIN_BINS; i++)
    {
      GSD.EdgeMasses[i] = (4.0*M_PI/3.0) * GSD.InternalDensity * pow(GSD.Edges[i], 3);
    }

  for(int i = 0; i < DL_GRAIN_BINS; i++)
    {
      GSD.Midpoints[i] = (GSD.Edges[i] + GSD.Edges[i+1]) / 2.0;
      GSD.Widths[i] = (GSD.Edges[i+1] - GSD.Edges[i]);
      double num = M_PI * GSD.InternalDensity * (pow(GSD.Edges[i+1], 4.0) - pow(GSD.Edges[i], 4.0));
      double denom = 3.0 * (GSD.Edges[i+1] - GSD.Edges[i]);
      GSD.AvgMasses[i] = num/denom;
    }

#if defined(DL_SHATTERING) || defined(DL_COAGULATION)
#ifdef DL_SHATTERING
  GSD.VelShat = 2.0 * 1.0e9; /* 2 km/s in um/s */
#endif

#ifdef DL_COAGULATION
  /* Follows Hirashita+ (2009), Equation 8. */
  double F_stick = 10.0;
  double gamma_coag = (25.0 + 12.0) / 2.0; /* erg / cm^2 */
  double E_coag = (5.4e11 + 3.4e10) / 2.0; /* dyn / cm^2 */
  double const_fac = pow(gamma_coag, 5.0/6.0) / (pow(E_coag, 1.0/3.0) * sqrt(All.GrainDensity)); /* all in cgs */
  for(int k = 0; k < DL_GRAIN_BINS; k++)
    {
      for(int i = 0; i <= k; i++)
        {
          double sum_3 = pow(GSD.Midpoints[k], 3) + pow(GSD.Midpoints[i], 3); /* um^3 */
          double totsum_3 = pow(GSD.Midpoints[k] + GSD.Midpoints[i], 3); /* um^3 */
          double a_fac = sqrt(sum_3 / totsum_3); /* dimensionless */
          double R_fac = GSD.Midpoints[k] * GSD.Midpoints[i] / (GSD.Midpoints[k] + GSD.Midpoints[i]); /* um */
          R_fac *= cm_per_um; /* cm */
          double v_ki = 2.14 * F_stick * a_fac / pow(R_fac, 5.0/6.0) * const_fac; /* cm/s, since everything's cgs */
          GSD.VelCoag[k][i] = v_ki / cm_per_um; /* um/s */
          GSD.VelCoag[i][k] = GSD.VelCoag[k][i]; /* um/s */
        }
    }
#endif

  /* Compute prefactors used to evaluate grain collision cross section
   * integrals in I_2_kj(). */
  for(int k = 0; k < DL_GRAIN_BINS; k++)
    {
      for(int j = 0; j <= k; j++)
        {
          double ak = GSD.Edges[k];
          double akp1 = GSD.Edges[k+1];
          double aj = GSD.Edges[j];
          double ajp1 = GSD.Edges[j+1];
          GSD.I_2_kj_Prefac[k][j] = ((aj - ajp1)*(ak - akp1)*(2*pow(aj,2) + 2*aj*ajp1 + 2*pow(ajp1,2) + 3*aj*(ak + akp1) + 3*ajp1*(ak + akp1) + 2*(pow(ak,2) + ak*akp1 + pow(akp1,2))))/6.;
          GSD.I_2_kj_Prefac[j][k] = GSD.I_2_kj_Prefac[k][j];
        }
    }

  /* Attempt to load an HDF5 file of grain velocity profiles.
   * Must have the following datasets:
   *   /VelocitiesCNM/NPoints
   *   /VelocitiesCNM/Radius (in um, sorted ascending)
   *   /VelocitiesCNM/Velocity (in um/s)
   *   /VelocitiesWIM/NPoints
   *   /VelocitiesWIM/Radius (in um, sorted ascending)
   *   /VelocitiesWIM/Velocity (in um/s)
   */
  mpi_printf("DUST_LIVE: Attempting to load grain velocity tables.\n");

  load_grain_velocity_data("VelocitiesCNM", GSD.VelocitiesCNM);
  load_grain_velocity_data("VelocitiesWIM", GSD.VelocitiesWIM);

#endif

#ifdef DL_SNE_DESTRUCTION
  mpi_printf("DUST_LIVE: Attempt to load supernova grain destruction tables.\n");

  load_sn_destruction_data();
#endif

#ifdef DL_PRODUCTION
  mpi_printf("DUST_LIVE: Attempting to load stellar grain size distributions.\n");

  /* Compute initial piecewise constant or piecewise linear grain size
   * distributions using values at bin edges. */
  for(int i = 0; i < DL_GRAIN_BINS; i++)
    {
      double dnda_lhs = agb_dnda(GSD.Edges[i]);
      double dnda_rhs = agb_dnda(GSD.Edges[i+1]);

      GSD.AGB_NumGrains[i] = (dnda_lhs + dnda_rhs) / 2.0 * GSD.Widths[i];
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
      GSD.AGB_BinSlopes[i] = (dnda_rhs - dnda_lhs) / GSD.Widths[i];
#endif
    }

  /* Attempt to load an HDF5 file of relative grain size distribution for SNII.
   * Must have the following datasets:
   *   /InitialSNII/NPoints
   *   /InitialSNII/Radius (in um, sorted ascending)
   *   /InitialSNII/GSD (relative)
   */
  load_snii_dnda();

#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
  normalize_dnda(GSD.AGB_NumGrains, GSD.AGB_BinSlopes);
  normalize_dnda(GSD.SNII_NumGrains, GSD.SNII_BinSlopes);
#else
  normalize_dnda(GSD.AGB_NumGrains);
  normalize_dnda(GSD.SNII_NumGrains);
#endif

  load_condensation_efficiencies();
#endif
}

#endif
#endif
