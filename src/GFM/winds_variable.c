/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/GFM/winds_variable.c
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

#include "../allvars.h"
#include "../proto.h"
#include "../voronoi.h"

#ifdef GFM_WINDS_VARIABLE

static double gfm_calc_variable_wind_parameters_specific_mass_loading(double velocity, double v_vir, double vel_disp, double utherm, double metallicity)
{
  double wind_momentum = All.VariableWindSpecMomentum;
#if defined(GFM_STELLAR_EVOLUTION) && (GFM_STELLAR_EVOLUTION == 0)

  /* calculate appropriate factors based on IMF */
  double FactorSN, WindEgySpecSN;
#ifdef GFM_CONST_IMF
  FactorSN = All.FactorSN;
  WindEgySpecSN = All.WindEgySpecSN;
#elif defined(GFM_VARIABLE_IMF)
#if (GFM_VARIABLE_IMF == 0)
  define_imf(vel_disp);
#else
#error "GFM_VARIABLE_IMF mode is not ok"
#endif
  FactorSN = calc_FactorSN();
  WindEgySpecSN = calc_WindEgySpecSN();
#endif

  /* wind_energy = factor * (SN energy per stellar mass formed of all stars) */
  double wind_energy = All.WindEnergyFactor * FactorSN * WindEgySpecSN;


#if defined(GFM_WIND_ENERGY_METAL_DEPENDENCE_TANH) && !defined(GFM_WIND_ENERGY_METAL_DEPENDENCE)
#error "GFM_WIND_ENERGY_METAL_DEPENDENCE_TANH requires GFM_WIND_ENERGY_METAL_DEPENDENCE"
#endif

#ifdef GFM_WIND_ENERGY_METAL_DEPENDENCE
#ifdef GFM_WIND_ENERGY_METAL_DEPENDENCE_TANH
  if(metallicity > 0)
    wind_energy *= 1.0 - (1.0 - All.WindEnergyReductionFactor) * 0.5 * (tanh( log10(metallicity / All.WindEnergyReductionMetallicity) / All.WindEnergyReductionExponent) + 1.0);
#else
  wind_energy *= (All.WindEnergyReductionFactor + (1 - All.WindEnergyReductionFactor) / (1 + pow(metallicity / All.WindEnergyReductionMetallicity, All.WindEnergyReductionExponent)));
#endif
#endif


#ifdef GFM_WINDS_HUBBLESCALING
  if(All.Time > 0)
    {
      double eff = hubble_function(All.Time) / hubble_function(1.0 / (All.WindSuppressionRedshift + 1.0));
      if(eff > 1)
        eff = 1;

      wind_energy = eff * All.WindEnergyFactor * FactorSN * WindEgySpecSN;
    }
#endif

#ifdef GFM_WINDS_MASSSCALING
  if(All.Time > 0)
    {
      double m_vir = pow(v_vir, 3) / (10 * All.G * All.cf_H);

      double eff = pow(hubble_function(All.Time) / All.Hubble, 2.0 / 3) * pow(All.VariableWindMassScale / m_vir, 1.0 / 3);

      if(eff > All.WindEnergyFactor)
        eff = All.WindEnergyFactor;

      wind_energy = eff * FactorSN * WindEgySpecSN;
    }
#endif

#else // for defined(GFM_STELLAR_EVOLUTION) && (GFM_STELLAR_EVOLUTION == 0)
  /* wind_energy = factor * (SN energy per stellar mass formed of long-lived stars) */
  double wind_energy = All.WindEnergyFactor * All.FactorSN * All.WindEgySpecSN / (1 - All.FactorSN);

#endif
  return (wind_energy + sqrt(wind_energy * wind_energy + velocity * velocity * wind_momentum * wind_momentum)) / (velocity * velocity + 2 * utherm);
}



void gfm_calc_variable_wind_parameters(double stellar_mass, double halo_val, double vel_disp, double metallicity, wind_parameter * wp)
{
  double v, v_vir, v_esc, utherm = 0;
  double halo_sigma;

#if (GFM_WINDS_VARIABLE==0)     /* infer wind velocity from halo mass */
  /* we take here the FOF mass */ .
    v_vir = pow(10.0 * All.G * All.cf_H * halo_mass, 1.0 / 3);

  v_vir *= 0.95;                /* empirical correction factor to better line up with the local velocity dispersion estimates obtained with GFM_WINDS_VARIABLE for well resolved halos */

  halo_sigma = v_vir / sqrt(2);
#endif


#if (GFM_WINDS_VARIABLE==1)     /* infer wind velocity from local velocity dipsersion */
  /* 3D -> 1D velocity dispersion (note that subfind_density() calculates local 3D dispersion) */
  halo_sigma = halo_val * GFM_FAC_SIGMA;
  v_vir = halo_sigma * sqrt(2); /* estimate ( found empirically to somewhere between halo_sigma*sqrt(3/2) and halo_sigma*sqrt(2) ) */
#endif

  v_esc = 3.5 * v_vir;          /* estimate */
  wp->v_esc_halo = v_esc;
  wp->v_vir_halo = v_vir;

#ifdef GFM_WINDS_VARIABLE_HUBBLE
  halo_sigma *= pow(All.Hubble / All.cf_H, 1.0 / 3);
#endif

  v = (All.VariableWindVelFactor * halo_sigma) > All.MinWindVel ? (All.VariableWindVelFactor * halo_sigma) : All.MinWindVel;

#ifdef GFM_WINDS_THERMAL_NEWDEF
  utherm = All.ThermalWindFraction / (1 - All.ThermalWindFraction) * v * v / 2;
#endif
#ifdef GFM_WINDS_THERMAL
  utherm = All.ThermalWindFactor * (3.0 / 2) * halo_sigma * halo_sigma;
#endif

  /* wind mass for this particle, assuming the wind is first given the energy wind_energy and then the momentum wind_momentum */
  wp->wind_mass = stellar_mass * gfm_calc_variable_wind_parameters_specific_mass_loading(v, v_vir, vel_disp, utherm, metallicity);
  wp->wind_velocity = v;
  wp->wind_utherm = utherm;
  wp->v_esc_halo = v_esc;
  wp->v_vir_halo = v_vir;
}


void init_variable_winds(void)
{
  int i, bins = 100;
  wind_parameter wp;
  double time_old = All.Time;
  All.Time = 1.0;               /* for gfm_calc_variable_wind_parameters to give z=0 results */
  set_cosmo_factors_for_current_time();

#if defined(DVR_RENDER) && (DVR_RENDER==0)
  return;
#endif

  char buf[1000];
  sprintf(buf, "%s/variable_wind_scaling.txt", All.OutputDir);
  FILE *fd = fopen(buf, "w");

#if defined(GFM_VARIABLE_IMF) && (GFM_VARIABLE_IMF==0)
  fprintf(fd, "NOTE: the following assumes 249.28 for vel_disp, i.e. assuming Salpeter IMF\n");
#endif

#if (GFM_WINDS_VARIABLE==0)
  double halo_mass, min_halo_mass = 1e-5, max_halo_mass = 1e5;
  double dlog10_halo_mass = log10(max_halo_mass / min_halo_mass) / bins;
  fprintf(fd, "#halo mass variable winds: eta(v_wind=400 [internal units]) = %g\n", gfm_calc_variable_wind_parameters_specific_mass_loading(400, 400.0, 249.28, 0, 0));
  fprintf(fd, "#variable wind scaling: min_halo_mass=%g    max_halo_mass=%g    bins=%d    redshift=0\n", min_halo_mass, max_halo_mass, bins);
  fprintf(fd, "#M_200[internal units]\tV_200[internal units]\tconcentration\tv_esc[internal units]\tv_wind[internal units]\teta\n");
  for(i = 0; i < bins; i++)
    {
      halo_mass = pow(10.0, (i + 0.5) * dlog10_halo_mass + log10(min_halo_mass));
      gfm_calc_variable_wind_parameters(1.0, halo_mass, 249.28, 0, &wp);
      fprintf(fd, "%e\t%e\t%e\t%e\t%e\t%e\n", halo_mass, wp.v_vir_halo, wp.c_halo, wp.v_esc_halo, wp.wind_velocity, wp.wind_mass);
    }
#endif
#if (GFM_WINDS_VARIABLE==1)
  double halo_sigma, min_halo_sigma = 1, max_halo_sigma = 10000;
  double dlog10_halo_sigma = log10(max_halo_sigma / min_halo_sigma) / bins;
  fprintf(fd, "#sigma variable winds: eta(v_wind=400 [internal units]) = %g\n", gfm_calc_variable_wind_parameters_specific_mass_loading(400.0, 400.0, 249.28, 0, 0));
  fprintf(fd, "#variable wind scaling: min_halo_sigma=%g    max_halo_sigma=%g    bins=%d\n", min_halo_sigma, max_halo_sigma, bins);
  fprintf(fd, "#3D sigma[internal units]\tV_200[internal units]\tv_wind[internal units]\teta\n");
  for(i = 0; i < bins; i++)
    {
      halo_sigma = pow(10.0, (i + 0.5) * dlog10_halo_sigma + log10(min_halo_sigma));
      gfm_calc_variable_wind_parameters(1.0, halo_sigma, 249.28, 0, &wp);
      fprintf(fd, "%e\t%e\t%e\t%e\n", halo_sigma, wp.v_vir_halo, wp.wind_velocity, wp.wind_mass);
    }
#endif
  fclose(fd);
  All.Time = time_old;          /* return to previous value */
  set_cosmo_factors_for_current_time();
}
#endif
