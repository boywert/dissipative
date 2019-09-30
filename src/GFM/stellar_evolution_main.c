/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/GFM/stellar_evolution_main.c
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


#ifdef GFM_STELLAR_EVOLUTION

#define WARN_TOLERANCE 0.00001



/*!
 * driver routine for stellar evolution, set up arguments for stellar_evolution(...) call
 */
void do_stellar_evolution(MyFloat age_of_star_in_Gyr, MyFloat dtime_in_Gyr, int iPart, stellar_evolution_data * sed)
{
  int k_elem;
  MyFloat initial_metals[GFM_N_CHEM_ELEMENTS];  /* original element abundances */
  MyFloat metallicity;
  MyDouble initial_mass;
  MyDouble log_min_mass, log_max_mass, log_metallicity;

  if(P[iPart].Type != 4)
    terminate("GFM_STELLAR_EVOLUTION: Not a stellar particle.\n");

  if(P[iPart].Mass <= 0)
    terminate("GFM_STELLAR_EVOLUTION: Bad mass.\n");

  if(STP(iPart).BirthTime <= 0)
    terminate("GFM_STELLAR_EVOLUTION: Wind particle.\n");


  /* initial fields from star particle properties */
  initial_mass = STP(iPart).InitialMass;
  for(k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
    initial_metals[k_elem] = STP(iPart).MassMetals[k_elem] / P[iPart].Mass;

  metallicity = STP(iPart).Metallicity;
  if(metallicity > 0)
    log_metallicity = log10(metallicity);
  else
    log_metallicity = GFM_MIN_METAL;

  /* zero variables and arrays */
  sed->total_mass_released = 0;
  sed->total_metal_mass_released = 0;
  for(k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
    sed->metal_mass_released[k_elem] = 0;

#if defined(GFM_DUST) || (defined(DUST_LIVE) && defined(DL_PRODUCTION))
  sed->total_dust_mass_released = 0.0;
#endif
#ifdef GFM_DUST
  for(int chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
    {
      for(k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
        {
          sed->dust_mass_released[chan][k_elem] = 0.0;
        }
    }
#endif
#if defined(DUST_LIVE) && defined(DL_PRODUCTION)
  for(k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
    {
      sed->dust_mass_released[k_elem] = 0.0;
    }
#endif
#ifdef GFM_CHEMTAGS
  for(k_elem = 0; k_elem < GFM_N_CHEM_TAGS; k_elem++)
    {
      sed->metal_mass_released_chemtags[k_elem] = 0;
    }
#endif

  sed->number_of_SNIa = 0;
  sed->number_of_SNII = 0;
  sed->AGB_mass_released = 0.0;
  sed->SNIa_mass_released = 0.0;
  sed->SNII_mass_released = 0.0;

  /* minimum and maximum mass of stars that will die during this time-step */
  if(get_dying_mass_in_Msun(age_of_star_in_Gyr, metallicity) <= 0 || get_dying_mass_in_Msun(age_of_star_in_Gyr + dtime_in_Gyr, metallicity) <= 0)
    terminate("GFM_STELLAR_EVOLUTION: Negative mass for dying stars");

  log_max_mass = log10(get_dying_mass_in_Msun(age_of_star_in_Gyr, metallicity));
  log_min_mass = log10(get_dying_mass_in_Msun(age_of_star_in_Gyr + dtime_in_Gyr, metallicity));

  if(log_min_mass > log_max_mass)
    terminate("GFM_STELLAR_EVOLUTION: Min mass larger than max mass for stellar evolution: i=%d  log_min_mass=%g  log_max_mass=%g  diff=%g\n", iPart,
              log_min_mass, log_max_mass, log_max_mass - log_min_mass);

  if(log_min_mass == log_max_mass)
    return;

  /* initialize the IMF */
#ifdef GFM_VARIABLE_IMF
#if (GFM_VARIABLE_IMF == 0)
  define_imf(STP(iPart).DMVelDisp);
#else
#error "GFM_VARIABLE_IMF mode is not ok"
#endif
#endif

  /* perform stellar evolution, returns: total mass, metall mass, element mass, SNIa, SNII rates via sed */
  evolve_AGB(log_min_mass, log_max_mass, log_metallicity, initial_metals, sed);

  sed->AGB_mass_released = sed->total_mass_released;

  evolve_SNIa(log_min_mass, log_max_mass, age_of_star_in_Gyr, dtime_in_Gyr, metallicity, sed);

  sed->SNIa_mass_released = sed->total_mass_released - sed->AGB_mass_released;

  evolve_SNII(log_min_mass, log_max_mass, log_metallicity, initial_metals, sed);

  sed->SNII_mass_released = sed->total_mass_released - sed->AGB_mass_released - sed->SNIa_mass_released;

  // Note there is no:  NSNS_mass_released parameter ... assuming Eu contributes nothing to mass or energy dump
#ifdef GFM_RPROCESS
  evolve_NSNS(log_min_mass, log_max_mass, age_of_star_in_Gyr, dtime_in_Gyr, metallicity, initial_mass, sed);
#endif

  /* convert mass and metals_released to code units */
  sed->total_mass_released *= initial_mass;
  sed->total_metal_mass_released *= initial_mass;
  for(k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
    sed->metal_mass_released[k_elem] *= initial_mass;

#if defined(GFM_DUST) || (defined(DUST_LIVE) && defined(DL_PRODUCTION))
  sed->total_dust_mass_released *= initial_mass;
#endif
#ifdef GFM_DUST
  for(int chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
    {
      for(k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
        {
          sed->dust_mass_released[chan][k_elem] *= initial_mass;
        }
    }
#endif
#if defined(DUST_LIVE) && defined(DL_PRODUCTION)
  for(k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
    {
      sed->dust_mass_released[k_elem] *= initial_mass;
    }
#endif
  sed->AGB_mass_released *= initial_mass;
  sed->SNIa_mass_released *= initial_mass;
  sed->SNII_mass_released *= initial_mass;

#ifdef GFM_CHEMTAGS
  sed->metal_mass_released_chemtags[GFM_SNIA_CHEMTAG] *= initial_mass;
  sed->metal_mass_released_chemtags[GFM_SNII_CHEMTAG] *= initial_mass;
  sed->metal_mass_released_chemtags[GFM_AGB_CHEMTAG] *= initial_mass;
#ifdef GFM_RPROCESS
  sed->metal_mass_released_chemtags[GFM_NSNS_CHEMTAG] *= initial_mass;
#endif
#ifdef GFM_SPLITFE
  sed->metal_mass_released_chemtags[GFM_FESNIA_CHEMTAG] *= initial_mass;
  sed->metal_mass_released_chemtags[GFM_FESNII_CHEMTAG] *= initial_mass;
#endif
#endif
  
  /* do same for number of SN of each type */
  sed->number_of_SNIa *= initial_mass * All.UnitMass_in_g / SOLAR_MASS / All.HubbleParam;
  sed->number_of_SNII *= initial_mass * All.UnitMass_in_g / SOLAR_MASS / All.HubbleParam;

  if(sed->total_mass_released > P[iPart].Mass)
    {
      terminate("GFM_STELLAR_EVOLUTION: released mass larger than stellar mass: ID=%lld P[iPart].Mass=%g sed->total_mass_released=%g initial_mass=%g age_Gyr=%g dt_Gyr=%g\n",
                (long long) P[iPart].ID, P[iPart].Mass, sed->total_mass_released, initial_mass, age_of_star_in_Gyr, dtime_in_Gyr);
    }

  if(sed->total_mass_released < 0)
    {
      if(sed->total_mass_released < -WARN_TOLERANCE * initial_mass)
        warn("GFM_STELLAR_EVOLUTION: negative mass returned: ID=%lld P[iPart].Mass=%g sed->total_mass_released=%g initial_mass=%g age_Gyr=%g dt_Gyr=%g\n",
             (long long) P[iPart].ID, P[iPart].Mass, sed->total_mass_released, initial_mass, age_of_star_in_Gyr, dtime_in_Gyr);

      /* fix */
      sed->total_mass_released = 0;
      sed->total_metal_mass_released = 0;
      for(int k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
        sed->metal_mass_released[k_elem] = 0;
#ifdef GFM_CHEMTAGS
      for(int k_elem = 0; k_elem < GFM_N_CHEM_TAGS; k_elem++)
        sed->metal_mass_released_chemtags[k_elem] = 0;
#endif
    }

  if(sed->total_metal_mass_released > sed->total_mass_released)
    {
      if(sed->total_metal_mass_released  > (1.0 + WARN_TOLERANCE) * sed->total_mass_released)
        warn
          ("GFM_STELLAR_EVOLUTION: released metal mass larger than total released mass: ID=%lld total_metal_mass_released=%g total_mass_released=%g",
           (long long) P[iPart].ID, sed->total_metal_mass_released, sed->total_mass_released);
      /* fix */
      sed->total_mass_released = 0;
      sed->total_metal_mass_released = 0;
      for(int k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
        sed->metal_mass_released[k_elem] = 0;
#ifdef GFM_CHEMTAGS
      for(int k_elem = 0; k_elem < GFM_N_CHEM_TAGS; k_elem++)
        sed->metal_mass_released_chemtags[k_elem] = 0;
#endif
    }

  if(sed->total_metal_mass_released < 0)
    {
      if(sed->total_metal_mass_released  < -WARN_TOLERANCE * initial_mass)
        terminate("GFM_STELLAR_EVOLUTION: negative metal mass returned: ID=%lld P[iPart].Mass=%g sed->total_metal_mass_released=%g initial_mass=%g\n",
             (long long) P[iPart].ID, P[iPart].Mass, sed->total_metal_mass_released, initial_mass);
      /* fix */
      sed->total_mass_released = 0;
      sed->total_metal_mass_released = 0;
      for(int k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
        sed->metal_mass_released[k_elem] = 0;
#ifdef GFM_CHEMTAGS
      for(int k_elem = 0; k_elem < GFM_N_CHEM_TAGS; k_elem++)
        sed->metal_mass_released_chemtags[k_elem] = 0;
#endif
    }

  for(k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)       /* H, He and metals */
    {
      if(sed->metal_mass_released[k_elem] < 0)
        {
          if(sed->metal_mass_released[k_elem]  < -WARN_TOLERANCE * initial_mass)
            {
              char buf[2000];
              sprintf(buf, "GFM_STELLAR_EVOLUTION: negative metal mass element returned: ID=%lld P[iPart].Mass=%g sed->metal_mass_released[k_elem]=%g initial_mass=%g k_elem=%d  age_of_star_in_Gyr=%g  dtime_in_Gyr=%g metallicity=%g  initial_metals=",
               (long long) P[iPart].ID, P[iPart].Mass, sed->metal_mass_released[k_elem], initial_mass, k_elem, age_of_star_in_Gyr, dtime_in_Gyr, metallicity);

              for(int kk = 0; kk < GFM_N_CHEM_ELEMENTS; kk++)
        	{
        	  char buf1[100];
        	  sprintf(buf1, "%g ", initial_metals[kk]);
        	  strcat(buf, buf1);
        	}
              strcat(buf, " metal_mass_released=");
              for(int kk = 0; kk < GFM_N_CHEM_ELEMENTS; kk++)
        	{
        	  char buf1[100];
        	  sprintf(buf1, "%g ", sed->metal_mass_released[kk]);
        	  strcat(buf, buf1);
        	}
              strcat(buf, "\n");
              warn(buf);
            }

          /* fix */
          sed->total_mass_released = 0;
          sed->total_metal_mass_released = 0;
          for(int k_elem2 = 0; k_elem2 < GFM_N_CHEM_ELEMENTS; k_elem2++)
            sed->metal_mass_released[k_elem2] = 0;
#ifdef GFM_CHEMTAGS
          for(int k_elem3 = 0; k_elem3 < GFM_N_CHEM_TAGS; k_elem3++)
            sed->metal_mass_released_chemtags[k_elem3] = 0;
#endif
        }
    }

#ifdef GFM_DUST
  for(int chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
    {
      for(k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
        {
          if(sed->dust_mass_released[chan][k_elem] < 0.0)
            {
              sed->dust_mass_released[chan][k_elem] = 0.0;
            }
        }
    }
#endif

#if (GFM_STELLAR_EVOLUTION == 1)
  sed->total_mass_released = 0;
#endif

#ifdef GFM_DISCRETE_ENRICHMENT
    double delta_M = sed->total_mass_released / P[iPart].Mass;
    double discrete_return_frac = 0.0001; // All.GFM_DiscreteEnrichFrac
    double discrete_min_age = 0.1; // Gyr, below which all returns are allowed (for fast evolution of early stellar populations)

    if(delta_M < discrete_return_frac && age_of_star_in_Gyr > discrete_min_age)
    {
      // zero all enrichment returns and quit early
      sed->total_mass_released = 0;
      sed->total_metal_mass_released = 0;
      for(k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
        sed->metal_mass_released[k_elem] = 0;

#ifdef GFM_CHEMTAGS
      for(k_elem = 0; k_elem < GFM_N_CHEM_TAGS; k_elem++)
        sed->metal_mass_released_chemtags[k_elem] = 0; 
#endif

      sed->number_of_SNIa = 0;
      sed->number_of_SNII = 0;
      sed->AGB_mass_released = 0.0;
      sed->SNIa_mass_released = 0.0;
      sed->SNII_mass_released = 0.0;
      return;
    }
#endif /* GFM_DISCRETE_ENRICHMENT */

#if defined(DUST_LIVE) && defined(DL_PRODUCTION)
  if(sed->total_dust_mass_released < 0.0)
    sed->total_dust_mass_released = 0.0;

  for(k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
    {
      if(sed->dust_mass_released[k_elem] < 0.0)
        {
          sed->dust_mass_released[k_elem] = 0.0;
        }
    }
#endif

  /* fractional mass loss */
  double frac = 1.0 - sed->total_mass_released / P[iPart].Mass;

  /* scale mass and metal content */
  P[iPart].Mass *= frac;

  for(k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
    STP(iPart).MassMetals[k_elem] *= frac;

#ifdef GFM_CHEMTAGS
    // NOTE: we assume NSNS production mass is negligable - does not factor into the total mass released
  for(int k_elem = 0; k_elem < GFM_N_CHEM_TAGS; k_elem++)
    STP(iPart).MassMetalsChemTags[k_elem] *= frac;
#endif
}

void evolve_active_stars(void)
{
  int i, iel;
  double time_begstep, dt, dtime, age_of_star_in_Gyr, age_of_star_in_Gyr_endstep, dtime_in_Gyr;
  stellar_evolution_data sed;
  double AGBLocalMassReleased = 0.0;
  double SNIaLocalMassReleased = 0.0;
  double SNIILocalMassReleased = 0.0;

#ifdef FM_STAR_FEEDBACK
  MyFloat SN_Number;
#endif
#ifdef GFM_DUST
  int chan;
#endif

  /* do local star particles */
  for(i = 0; i < Nstar; i++)
    {
      if(All.ComovingIntegrationOn)
        time_begstep = All.TimeBegin * exp(All.Ti_begstep[P[StarParticle[i].index].TimeBinGrav] * All.Timebase_interval);
      else
        time_begstep = All.TimeBegin + All.Ti_begstep[P[StarParticle[i].index].TimeBinGrav] * All.Timebase_interval;

      dt = (P[StarParticle[i].index].TimeBinGrav ? (((integertime) 1) << P[StarParticle[i].index].TimeBinGrav) : 0) * All.Timebase_interval;

      if(All.ComovingIntegrationOn)
        dtime = All.Time * dt / All.cf_time_hubble_a;
      else
        dtime = dt;

      dtime *= All.UnitTime_in_s / All.HubbleParam;
      dtime_in_Gyr = dtime / SEC_PER_MEGAYEAR / 1000;

      age_of_star_in_Gyr_endstep = get_time_difference_in_Gyr(STP(StarParticle[i].index).BirthTime, time_begstep);      /* (Note: All.Ti_begstep[] has already been advanced for the next step at this point) */
      age_of_star_in_Gyr = age_of_star_in_Gyr_endstep - dtime_in_Gyr;

#ifdef GFM_DISCRETE_ENRICHMENT
      /* we take the age at the beginning of the interval, and dtime is since the last enrichment event */
      age_of_star_in_Gyr = get_time_difference_in_Gyr(STP(StarParticle[i].index).BirthTime, STP(StarParticle[i].index).lastEnrichTime);
      dtime_in_Gyr = get_time_difference_in_Gyr(STP(StarParticle[i].index).lastEnrichTime, time_begstep);
#endif

      do_stellar_evolution(age_of_star_in_Gyr, dtime_in_Gyr, StarParticle[i].index, &sed);

#ifdef GFM_DISCRETE_ENRICHMENT
      /* if we will do an enrichment event, update the lastEnrichTime to the end of the current timestep */
      if(sed.total_mass_released > 0.0)
        STP(StarParticle[i].index).lastEnrichTime = time_begstep;
#endif

      AGBLocalMassReleased += sed.AGB_mass_released;
      SNIaLocalMassReleased += sed.SNIa_mass_released;
      SNIILocalMassReleased += sed.SNII_mass_released;

      if(dtime_in_Gyr > 0)
        {
          STP(StarParticle[i].index).SNIaRate = sed.number_of_SNIa / (dtime_in_Gyr * 1.e9);
          STP(StarParticle[i].index).SNIIRate = sed.number_of_SNII / (dtime_in_Gyr * 1.e9);
        }

      StarParticle[i].TotalMassReleased = sed.total_mass_released;
      StarParticle[i].TotalMetalMassReleased = sed.total_metal_mass_released;
      for(iel = 0; iel < GFM_N_CHEM_ELEMENTS; iel++)
        {
          StarParticle[i].MetalMassReleased[iel] = sed.metal_mass_released[iel];
        }
#ifdef GFM_DUST
      for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
        {
          for(iel = 0; iel < GFM_N_CHEM_ELEMENTS; iel++)
            {
              StarParticle[i].DustMassReleased[chan][iel] = sed.dust_mass_released[chan][iel];
            }
        }
#endif
#if defined(GFM_DUST) || defined(DUST_LIVE)
      StarParticle[i].NumSNII = sed.number_of_SNII;
#endif
#if defined(DUST_LIVE) && defined(DL_PRODUCTION)
      STP(StarParticle[i].index).DeltaDustMassTot = 0.0;
      for(iel = 0; iel < GFM_N_CHEM_ELEMENTS; iel++)
        {
          STP(StarParticle[i].index).DeltaDustMass[iel] = sed.dust_mass_released[iel];
          STP(StarParticle[i].index).DeltaDustMassTot += sed.dust_mass_released[iel];
        }
      /* Use total mass returned by stars of different types as a proxy for
       * determining dominant contributor to grain size distribution. */
      if(sed.AGB_mass_released > dmax(sed.SNII_mass_released, sed.SNIa_mass_released))
        STP(StarParticle[i].index).DndaType = GSD_DNDA_AGB;
      else if(sed.SNII_mass_released > sed.SNIa_mass_released)
        STP(StarParticle[i].index).DndaType = GSD_DNDA_SNII;
      else
        STP(StarParticle[i].index).DndaType = GSD_DNDA_SNIA;
#endif

#ifdef GFM_CHEMTAGS
      for(iel = 0; iel < GFM_N_CHEM_TAGS; iel++)
        {
          StarParticle[i].MetalMassReleasedChemTags[iel] = sed.metal_mass_released_chemtags[iel];
        }
#endif

#ifdef GFM_STELLAR_FEEDBACK
      StarParticle[i].SNIaEnergyReleased = sed.number_of_SNIa * All.EnergyPerSNIa * (All.HubbleParam / All.UnitEnergy_in_cgs);
      StarParticle[i].AGBMomentumReleased = sed.AGB_mass_released * All.AGBWindVelocity;
#endif
#ifdef GFM_WINDS_LOCAL
      StarParticle[i].WindEnergyReleased = sed.number_of_SNII * All.WindEnergyIn1e51erg * GFM_SNII_ENERGY * (All.HubbleParam / All.UnitEnergy_in_cgs);
#endif

#ifdef FM_STAR_FEEDBACK
#ifdef FM_SN_COOLING_RADIUS_BOOST
      StarParticle[i].NumSNII = sed.number_of_SNII;
      StarParticle[i].NumSNIa = sed.number_of_SNIa;
#endif
#ifndef INSTANTANEOUS_DEPOSITION
      SN_Number = sed.number_of_SNII;
#ifdef FM_VAR_SN_EFF       /*PAM*/
      double metallicity_loc = StarParticle[i].AvgMetalNgb;
      double metal_floor = 0.001*GFM_SOLAR_METALLICITY;
      if(metallicity_loc < metal_floor)
                 metallicity_loc = metal_floor;

      double sn_eff = 0.3 + (3.0 - 0.3)/(1. + metallicity_loc / (0.1 * GFM_SOLAR_METALLICITY));
      SN_Number *= sn_eff;  /* WARNING: now is not the SN number but is weighted by SN efficiency */
#endif
#ifdef DIRECT_MOMENTUM_INJECTION_FEEDBACK
      StarParticle[i].TotalMomentumReleased = compute_SN_momentum(SN_Number);
#endif
      StarParticle[i].TotalEnergyReleased = compute_SN_energy(SN_Number);
#else
      StarParticle[i].TotalEnergyReleased = compute_total_SN_energy(StarParticle[i].index, age_of_star_in_Gyr + dtime_in_Gyr, &SN_Number);
#endif
#ifdef OUTPUT_STELLAR_FEEDBACK
      STP(StarParticle[i].index).SNII_Num = SN_Number;
      STP(StarParticle[i].index).Cum_SNII_Num += SN_Number;
      STP(StarParticle[i].index).FeedbackEnergy = StarParticle[i].TotalEnergyReleased;
      STP(StarParticle[i].index).FeedbackThEnergy = (1.0 - All.EtaKineticEnergy) * StarParticle[i].TotalEnergyReleased;
      STP(StarParticle[i].index).FeedbackKinEnergy = All.EtaKineticEnergy * StarParticle[i].TotalEnergyReleased;
      STP(StarParticle[i].index).TotalMassReleased = StarParticle[i].TotalMassReleased;
      STP(StarParticle[i].index).TotalMassToEnrich = StarParticle[i].NumNgb;
#endif
#ifdef DELAYED_COOLING
      double tol = 1.0e-5;
      StarParticle[i].BlastRadius = compute_blast_radius(StarParticle[i].TotalEnergyReleased, StarParticle[i].AvgHDens * All.cf_a3inv, StarParticle[i].AvgPress * All.cf_a3inv) / All.cf_atime;
      StarParticle[i].BlastRadius = dmax(StarParticle[i].BlastRadius, StarParticle[i].MinBlastRadius + tol);
      StarParticle[i].BlastRadius = dmin(StarParticle[i].BlastRadius, STP(StarParticle[i].index).Hsml);
      StarParticle[i].CoolShutoffTime = compute_cooling_shutoff_time(StarParticle[i].TotalEnergyReleased, StarParticle[i].AvgHDens * All.cf_a3inv, StarParticle[i].AvgPress * All.cf_a3inv);
      STP(StarParticle[i].index).BlastRadius = StarParticle[i].BlastRadius;
#endif
#endif

#ifdef FM_RADIATION_FEEDBACK

      if(age_of_star_in_Gyr >= All.InputTimeHeatRadiationFeedback)
        {
          StarParticle[i].StromgrenRadius = -999.;
          StarParticle[i].RadiationMomentumReleased = 0.;
          StarParticle[i].RadCoolShutoffTime = 0.;
          //StarParticle[i].RadVelocityKick = 0.;
#ifdef FM_STOCHASTIC_HII_PHOTOIONIZATION
          StarParticle[i].StromgrenMass   = -999.;
#endif
        }
      else
        {
          double Lum = All.LumToMassRatioRadiationFeedback * P[StarParticle[i].index].Mass;     /* luminosity in code units */
          Lum *= (All.UnitMass_in_g / SOLAR_MASS) * SOLAR_LUM / All.HubbleParam;        /* in CGS, HubbleParam because of mass scaling, LVS?? */
          double avgGasdens = StarP[i].LocISMdens;						/* comoving density */

          double minGasDist = StarParticle[i].RadFeed_MinGasDist * All.cf_atime;        /* physical */
          StarParticle[i].StromgrenRadius = compute_stromgren_radius(Lum, avgGasdens * All.cf_a3inv, minGasDist);
          StarParticle[i].StromgrenRadius /= All.cf_atime;      /*comoving code units */
#ifdef FM_STOCHASTIC_HII_PHOTOIONIZATION
          StarParticle[i].StromgrenMass   = compute_stromgren_mass(Lum, avgGasdens * All.cf_a3inv);
#endif

#if !defined(FM_STOCHASTIC_HII_PHOTOIONIZATION)
          if(age_of_star_in_Gyr >= All.InputTimeMomRadiationFeedback && STP(StarParticle[i].index).RadFeed_Flag == 0)
            {
              StarParticle[i].RadiationMomentumReleased = Lum / CLIGHT * age_of_star_in_Gyr * SEC_PER_GIGAYEAR; /* gr (cm/s) h^-1 ... */
              StarParticle[i].RadiationMomentumReleased /= (All.UnitVelocity_in_cm_per_s * All.UnitMass_in_g);  /* code units, although missing cf_atime, added at use in radiation_stellar_feedback.c */
              StarParticle[i].RadiationMomentumReleased *= All.HubbleParam;

              STP(StarParticle[i].index).RadFeed_Flag = -1;     /* inputs all momentum only once */
            }
          else
            StarParticle[i].RadiationMomentumReleased = 0.;

          StarParticle[i].RadCoolShutoffTime = (All.InputTimeHeatRadiationFeedback - age_of_star_in_Gyr) * SEC_PER_GIGAYEAR * All.HubbleParam / All.UnitTime_in_s;      /* in code units */
#else
            if(age_of_star_in_Gyr <= All.InputTimeMomRadiationFeedback)
            {
                StarParticle[i].RadiationMomentumReleased = Lum / CLIGHT * dtime_in_Gyr * SEC_PER_GIGAYEAR; /* gr (cm/s) h^-1 ... */
                StarParticle[i].RadiationMomentumReleased /= (All.UnitVelocity_in_cm_per_s * All.UnitMass_in_g);  /* code units, although missing cf_atime, added at use in radiation_stellar_feedback.c */
                StarParticle[i].RadiationMomentumReleased *= All.HubbleParam;
                
                STP(StarParticle[i].index).RadFeed_Flag = 0;     /* inputs all momentum only once */
            }
            else
                StarParticle[i].RadiationMomentumReleased = 0.;
            StarParticle[i].RadCoolShutoffTime = dtime_in_Gyr * SEC_PER_GIGAYEAR * All.HubbleParam / All.UnitTime_in_s;      /* in code units */
#endif
        }

      STP(StarParticle[i].index).StromgrenRadius = StarParticle[i].StromgrenRadius;
#ifdef FM_STOCHASTIC_HII_PHOTOIONIZATION
      STP(StarParticle[i].index).StromgrenMass   = StarParticle[i].StromgrenMass;
#endif
#ifdef FM_RADIATION_FEEDBACK_DEBUG
      STP(StarParticle[i].index).RadiationMomentumReleased = StarParticle[i].RadiationMomentumReleased;
      STP(StarParticle[i].index).RadCoolShutoffTime = StarParticle[i].RadCoolShutoffTime;
#endif
#endif

#ifdef FM_EARLY_STAR_FEEDBACK
      //if (age_of_star_in_Gyr > 0.015)    /*After 15 Myr no more input */
      if(age_of_star_in_Gyr > 0.004)    /*After 4 Myr no more input */
        {
          StarParticle[i].EarlyTotalEnergyReleased = 0;
        }
      else
        {
          StarParticle[i].EarlyTotalEnergyReleased = compute_EarlyFeedback_energy(age_of_star_in_Gyr, P[StarParticle[i].index].Mass, dtime);    /* dtime in sec, returns in code units  */
          //printf("FM_EARLY ENERGY=, age=%g, dt=%g, Energy=%g \n",age_of_star_in_Gyr, dtime, StarParticle[i].EarlyTotalEnergyReleased/(All.HubbleParam / All.UnitEnergy_in_cgs));  /* output energy in cgs, time in sec */
        }
      STP(StarParticle[i].index).EarlyCumulativeFeedback += StarParticle[i].EarlyTotalEnergyReleased;
#endif
    }

#ifdef VERBOSE
  double AGBGlobalMassReleased, SNIaGlobalMassReleased, SNIIGlobalMassReleased;

  MPI_Allreduce(&AGBLocalMassReleased, &AGBGlobalMassReleased, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&SNIaLocalMassReleased, &SNIaGlobalMassReleased, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&SNIILocalMassReleased, &SNIIGlobalMassReleased, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  All.AGBMassReleased += AGBGlobalMassReleased;
  All.SNIaMassReleased += SNIaGlobalMassReleased;
  All.SNIIMassReleased += SNIIGlobalMassReleased;
#endif
}


#endif
