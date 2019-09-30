/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dust_live/dust_transfer.c
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

#ifdef DUST_LIVE
#ifdef DL_GRAIN_BINS

#if defined(DL_GROWTH) || defined(DL_SPUTTERING)
void do_dust_transfer(void)
{
  int i;
  double dt_base, dt_step, hubble_a;

  start_dust();
  find_drag_cells(Ndust);
  drag_kernel();
  TIMER_STOPSTART(CPU_DUST, CPU_DUST_TRANSFER);

  if(All.ComovingIntegrationOn)
    {
      hubble_a = hubble_function(All.Time);
    }

  for(i = 0; i < Ndust; i++)
    {
      int p = DustParticle[i].index;

      dt_base = (P[p].TimeBinHydro ? (((integertime) 1) << P[p].TimeBinHydro) : 0) * All.Timebase_interval;

      if(All.ComovingIntegrationOn)
        dt_step = dt_base / hubble_a;
      else
        dt_step = dt_base;

      /* Calculate adot from various physical processes.  Keep adot in units
       * of um/Gyr during each update. */
      double adot = 0.0;
#ifdef DL_GROWTH
      /* Hirashita & Voshchinnikov (2014), Equations 5 to 8 */
      double rho_H_ref = (1.0e3 * PROTONMASS) / (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam);
      double T_ref = 10.0;
      double Z_ref = 0.0127;
      double S_ref = 0.3;
      double S_local = S_ref;
      double rho_H_local = DustParticle[i].LocalGasDensityH * All.cf_a3inv;
      double T_local = DustParticle[i].LocalGasTemp;
      double Z_local = DustParticle[i].LocalGasZ;

      double adot_growth = (Z_local/Z_ref) * (rho_H_local/rho_H_ref) * sqrt(T_local/T_ref) * (S_local/S_ref);
      adot += adot_growth;
#endif
#ifdef DL_SPUTTERING
      /* Tsai & Matthews (1995), Equation 14, converted to um/Gyr */
      double rho_local_cgs = DTP(p).LocalGasDensity * All.cf_a3inv * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
      double T_fac = 1.0 + pow(2.0e6/DustParticle[i].LocalGasTemp, 2.5);
      double h_fac = 1009.28; /* cm^3 um / Gyr */
      double adot_sputter = -h_fac * (rho_local_cgs / PROTONMASS) / T_fac;
      adot += adot_sputter;
#endif

      /* Since [a] = um and [adot] = um / Gyr, need timestep in Gyr. */
      double dt_Gyr = dt_step * All.UnitTime_in_s / All.HubbleParam / SEC_PER_GIGAYEAR;
      update_grain_sizes(i, adot, dt_Gyr, dt_step);
    }

  /* Attempt to remove or add metal mass from or to surrounding gas cells. */
  dust_transfer_kernel();

  /* Correct the grain size distribution update in case there weren't enough
   * metals in some gas cells. */
  for(i = 0; i < Ndust; i++)
    {
      update_dust_element_fractions(i);
      correct_grain_conserved_quantities(i);
      check_dust_for_removal(i);

      /* Rescale the grain size distribution by a small amount if necessary to
       * ensure integrated mass matches the particle mass. */
      double mass = 0.0;
      int p = DustParticle[i].index;
      for(int j = 0; j < DL_GRAIN_BINS; j++)
        {
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
          mass += DTP(p).NumGrains[j] * bin_avg_mass(j, DTP(p).NumGrains[j], DTP(p).BinSlopes[j]);
#else
          mass += DTP(p).NumGrains[j] * GSD.AvgMasses[j];
#endif
        }
      double mass_ratio = (mass > 0.0) ? P[p].Mass / mass : 0.0;
      /* Just scale the grain size distribution if possible. */
      if(mass_ratio > 0.0)
        {
          for(int j = 0; j < DL_GRAIN_BINS; j++)
            {
              DTP(p).NumGrains[j] *= mass_ratio;
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
              DTP(p).BinSlopes[j] *= mass_ratio;
#endif
            }
        }
      /* Because of floating point arithmetic, if there are not grains left
       * matching the particle mass, put them in the smallest bin. */
      else
        {
          DTP(p).NumGrains[0] = P[p].Mass / GSD.AvgMasses[0];
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
          DTP(p).BinSlopes[0] = 0.0;
#endif
        }
    }

  TIMER_STOPSTART(CPU_DUST_TRANSFER, CPU_DUST);
  end_dust();
}
#endif

#ifdef DL_SNE_DESTRUCTION
void do_dust_transfer_sne(void)
{
  int i;
  double dt_base, dt_step, hubble_a;

  start_dust();
  find_drag_cells(Ndust);
  drag_kernel();
  TIMER_STOPSTART(CPU_DUST, CPU_DUST_TRANSFER);

  if(All.ComovingIntegrationOn)
    {
      hubble_a = hubble_function(All.Time);
    }

  for(i = 0; i < Ndust; i++)
    {
      int p = DustParticle[i].index;

      dt_base = (P[p].TimeBinHydro ? (((integertime) 1) << P[p].TimeBinHydro) : 0) * All.Timebase_interval;

      if(All.ComovingIntegrationOn)
        dt_step = dt_base / hubble_a;
      else
        dt_step = dt_base;

      /* Since [a] = um and [adot] = um / Gyr, need timestep in Gyr. */
      double dt_Gyr = dt_step * All.UnitTime_in_s / All.HubbleParam / SEC_PER_GIGAYEAR;

      /* Supernova destruction happens in a subgrid fashion, not by specifying
       * a form of da/dt as for growth or sputtering. */
      update_grain_sizes_sne(i, dt_Gyr, dt_step);
    }

  /* Attempt to remove or add metal mass from or to surrounding gas cells. */
  dust_transfer_kernel();

  /* Correct the grain size distribution update in case there weren't enough
   * metals in some gas cells. */
  for(i = 0; i < Ndust; i++)
    {
      update_dust_element_fractions(i);
      correct_grain_conserved_quantities(i);
      check_dust_for_removal(i);

      /* Rescale the grain size distribution by a small amount if necessary to
       * ensure integrated mass matches the particle mass. */
      double mass = 0.0;
      int p = DustParticle[i].index;
      for(int j = 0; j < DL_GRAIN_BINS; j++)
        {
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
          mass += DTP(p).NumGrains[j] * bin_avg_mass(j, DTP(p).NumGrains[j], DTP(p).BinSlopes[j]);
#else
          mass += DTP(p).NumGrains[j] * GSD.AvgMasses[j];
#endif
        }
      double mass_ratio = (mass > 0.0) ? P[p].Mass / mass : 0.0;
      /* Just scale the grain size distribution if possible. */
      if(mass_ratio > 0.0)
        {
          for(int j = 0; j < DL_GRAIN_BINS; j++)
            {
              DTP(p).NumGrains[j] *= mass_ratio;
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
              DTP(p).BinSlopes[j] *= mass_ratio;
#endif
            }
        }
      /* Because of floating point arithmetic, if there are not grains left
       * matching the particle mass, put them in the smallest bin. */
      else
        {
          DTP(p).NumGrains[0] = P[p].Mass / GSD.AvgMasses[0];
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
          DTP(p).BinSlopes[0] = 0.0;
#endif
        }
    }

  TIMER_STOPSTART(CPU_DUST_TRANSFER, CPU_DUST);
  end_dust();
}
#endif

#if defined(DL_GROWTH) || defined(DL_SPUTTERING) || defined(DL_SNE_DESTRUCTION)
void update_dust_element_fractions(int i)
{
  int p = DustParticle[i].index;

  double elem_masses[GFM_N_CHEM_ELEMENTS];
  for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
    {
      elem_masses[k] = P[p].Mass * DTP(p).MetalFractions[k] + DustParticle[i].DeltaMetalMasses[k];
      /* Reset to zero in the event of a tiny negative mass remaining from all
       * of a species' mass being lost and floating point arithmetic.  This
       * does not affect the particle mass, just the normalized dust metal
       * fractions. */
      if(elem_masses[k] < 0.0)
        elem_masses[k] = 0.0;
    }

  double sum = 0.0;
  for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
    {
      sum += elem_masses[k];
    }
  for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
    {
      DTP(p).MetalFractions[k] = (sum > 0.0 ? elem_masses[k] / sum : 0.0);
    }
}

void check_dust_for_removal(int i)
{
  int p = DustParticle[i].index;

  if(P[p].Mass <= 0.0)
    {
      P[p].Mass = 0.0;
      /* If we want to remove a dust particle, its idx in TimeBinsDust (stored
       * when filling in DustParticle) is not the same as its idx in
       * TimeBinsGravity.  Manually search TimeBinsGravity for the idx of this
       * dust particle.  This should not be too inefficient because the number
       * of removed dust particles is not huge. */
      int active_grav_idx = -1;
      for(int idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
        {
          int i_grav = TimeBinsGravity.ActiveParticleList[idx];
          if((i_grav >= 0) && (i_grav == p))
            {
              active_grav_idx = idx;
              break;
            }
        }
      if(active_grav_idx == -1)
        terminate("Could not find active dust particle in the gravity time bins!");

      timebin_remove_particle(&TimeBinsDust, DustParticle[i].active_idx, P[p].TimeBinHydro);
      timebin_remove_particle(&TimeBinsGravity, active_grav_idx, P[p].TimeBinGrav);
    }
}
#endif

#if defined(DL_SNE_DESTRUCTION) || defined(DL_SHATTERING) || defined(DL_COAGULATION)
void do_shattering_coagulation()
{
  start_dust();
  find_drag_cells(Ndust);
  drag_kernel();
  TIMER_STOPSTART(CPU_DUST, CPU_DUST_SHATTER);

  /* Prepare to do dust-dust neighbor searches. */
  begin_shattering();
  /* Do the dust-dust neighbor searches to get dust densities. */
  dust_findHsml();
  /* May need to communicate some results depending on the tree structure. */
  exchange_shattering_results();

  /* The above calculations obtained dust densities that may be used even when
   * shattering and coagulation are not active (e.g. in supernova destruction).
   * We only need to proceed if shattering and coagulation take place. */

#if defined(DL_SHATTERING) || defined(DL_COAGULATION)
  /* Using dust densities, actually update grain size distributions. */
  /* Entirely local, no mass transfer to gas.  DustP also contains recent
   * CloudFrac estimates. */
  int i;
  double dt_base, dt_step, hubble_a;

  if(All.ComovingIntegrationOn)
    {
      hubble_a = hubble_function(All.Time);
    }

  for(i = 0; i < Ndust; i++)
    {
      int p = DustParticle[i].index;

      dt_base = (P[p].TimeBinHydro ? (((integertime) 1) << P[p].TimeBinHydro) : 0) * All.Timebase_interval;

      if(All.ComovingIntegrationOn)
        dt_step = dt_base / hubble_a;
      else
        dt_step = dt_base;

      /* Use timestep in Gyr, with mass loss rates calculated as mass per Gyr. */
      double dt_Gyr = dt_step * All.UnitTime_in_s / All.HubbleParam / SEC_PER_GIGAYEAR;
#ifdef DL_SHATTERING
      update_grain_sizes_shattering(i, dt_Gyr, dt_step, GSD_SHATTERING);
#endif
#ifdef DL_COAGULATION
      update_grain_sizes_shattering(i, dt_Gyr, dt_step, GSD_COAGULATION);
#endif
    }
#endif

  /* Free temporary memory. */
  end_shattering();

  TIMER_STOPSTART(CPU_DUST_SHATTER, CPU_DUST);
  end_dust();
}
#endif

#ifdef DL_SNE_DESTRUCTION
void update_sn_rates(void)
{
  double hubble_a;
  if(All.ComovingIntegrationOn)
    hubble_a = hubble_function(All.Time);

  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Type == 0)
        {
          if((P[i].Mass == 0) && (P[i].ID == 0))
            continue;

          double dt_base = (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;
          double dt_step;

          if(All.ComovingIntegrationOn)
            dt_step = dt_base / hubble_a;
          else
            dt_step = dt_base;

          /* Use timestep in Gyr, with mass loss rates calculated as mass per Gyr. */
          double dt_Gyr = dt_step * All.UnitTime_in_s / All.HubbleParam / SEC_PER_GIGAYEAR;
          if(dt_Gyr > 0.0)
            SphP[i].SNRate = SphP[i].NumSNII / dt_Gyr;
          else
            SphP[i].SNRate = 0.0;

          SphP[i].NumSNII = 0.0;
        }
    }
}
#endif

#endif
#endif
