/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/ism/ism_sfr.c
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

#include <math.h>
#include "../allvars.h"
#include "../proto.h"

#ifdef ISM

#define CRIT_DENS_TOLERANCE  (1.0 / 3.0)

void init_star_formation()
{
  double XH = HYDROGEN_MASSFRAC;
  double tff_at_threshold;
  double mol_weight;
#ifdef REFINEMENT
  double avg_gas_mass, jeans_density, jeans_mass, c_sound;
#endif

  mpi_printf("\nInitializing SF module\n");

#ifdef GFM_COOLING_METAL
  XH = GFM_INITIAL_ABUNDANCE_HYDROGEN;
#endif

  if(All.ComovingIntegrationOn)
    All.OverDensThresh = All.CritOverDensity * All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

/*
 Converting the density threshold from h^2 cm^-3 to code units
 neglecting metal contribution and assuming neutrality for mean 
 molecular weight. If the threshold is set to zero the code computes,
 automatically this quantity according to the definition of the 
 Jeans mass (only possible if refinement is enabled). 
 It is the physical density, NOT the comoving one!!! 
*/
  mol_weight = 4.0 / (1.0 + 3.0 * XH);

#ifdef REFINEMENT
  avg_gas_mass = All.TargetGasMassFactor * All.ReferenceGasPartMass;
  c_sound = sqrt(GAMMA * BOLTZMANN * 600. / (mol_weight * PROTONMASS)) / All.UnitVelocity_in_cm_per_s;
  jeans_density = pow(M_PI, 5.0) * pow(c_sound * c_sound / All.G, 3.0) / pow(6 * sqrt(8) * avg_gas_mass, 2.0);
  double jeans_lenght = sqrt(M_PI * c_sound * c_sound / All.G / jeans_density);
  jeans_mass = 4.0 / 3.0 * M_PI * jeans_density * pow(jeans_lenght / 2.0, 3.0);

  /* adjusting the critical density such that the jeans mass = 8xavg_gas_mass */
  /* decreasing then this value by a factor of 3 for safety (see also Teyssier+12) */
  if(All.DensThreshold == 0)
    All.DensThreshold = jeans_density * CRIT_DENS_TOLERANCE;
  else
#endif
    All.DensThreshold = mol_weight * PROTONMASS * All.DensThreshold / All.UnitDensity_in_cgs;

  tff_at_threshold = 0.25 * sqrt(1.5 * M_PI / (All.G * All.DensThreshold));

#ifdef USE_SFR
  All.MaxSfrTimescale = tff_at_threshold / All.SfrEfficiency;

  integrate_sfr();              /* compute KS law for adopted setting */

  mpi_printf("KS law computed\n");
  mpi_printf("SF module initialized\n\n");
#endif

}

#ifdef USE_SFR

/*
 * This routine does cooling and star formation for
 * the ISM model.
 */
void cooling_and_starformation(void)
{
  CPU_Step[CPU_MISC] += measure_time();

  int idx, i, bin;
  double tsfr, cloudmass, dens, fH2;
  short flag_sf;

  /* clear the SFR stored in the active timebins */
  for(bin = 0; bin < TIMEBINS; bin++)
    if(TimeBinSynchronized[bin])
      TimeBinSfr[bin] = 0;

  set_cosmo_factors_for_current_time();

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Type == 0)
        {
          if(P[i].Mass == 0 && P[i].ID == 0)
            continue;           /* skip cells that have been swallowed or eliminated */

          flag_sf = 0;

          dens = SphP[i].Density;

#ifdef ISM_H2_SFR
          fH2 = get_ism_h2_frac(i);
#else
          fH2 = 1.;
#endif

          tsfr = All.MaxSfrTimescale * sqrt(All.DensThreshold / (dens * All.cf_a3inv));
          cloudmass = fH2 * P[i].Mass;

          /* SF only if gas density above threshold */
          SphP[i].Sfr = 0.0;

          if(dens * All.cf_a3inv >= All.DensThreshold)
            flag_sf = 1;

          if(All.ComovingIntegrationOn) /* to protect against SF at too high redshift */
            if(dens < All.OverDensThresh)
              flag_sf = 0;

          cool_cell(i);

          if(flag_sf == 1)
            {
              SphP[i].Sfr = cloudmass / tsfr;
              SphP[i].Sfr *= (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);
              TimeBinSfr[P[i].TimeBinHydro] += SphP[i].Sfr;
            }
        }
    }                           /* end of main loop over active particles */

  CPU_Step[CPU_COOLINGSFR] += measure_time();
}


double get_starformation_rate(int i)
{
  double tsfr, dens;
  double rateOfSF = 0.0;
  double fH2 = 1.0;
  double cloudmass = P[i].Mass;

  set_cosmo_factors_for_current_time();

  dens = SphP[i].Density;
  tsfr = All.MaxSfrTimescale * sqrt(All.DensThreshold / (dens * All.cf_a3inv));

#ifdef ISM_H2_SFR
  fH2 = get_ism_h2_frac(i);
#else
  fH2 = 1.;
#endif

  /* star fomation only above threshold */
  if(dens * All.cf_a3inv >= All.DensThreshold)
    rateOfSF = fH2 * cloudmass / tsfr;

  if(All.ComovingIntegrationOn) /* to protect against SF at too high redshift */
    if(dens < All.OverDensThresh)
      rateOfSF = 0.0;


  /* convert to solar masses per yr */
  rateOfSF *= (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

  return rateOfSF;
}

void integrate_sfr(void)
{
  double rho0, rho, q, dz, gam, sigma = 0, sigma_u4, sigmasfr = 0;
  double P, drho, dq;
  double tsfr, z, meanweight, u4;
  FILE *fd = NULL;

  meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));   /* note: assuming FULL ionization */
  u4 = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * 1.0e4;
  u4 *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

  if(WriteMiscFiles && (ThisTask == 0))
    fd = fopen("sfrrate.txt", "w");

  for(rho0 = All.DensThreshold; rho0 <= 10000 * All.DensThreshold; rho0 *= 1.02)
    {
      z = 0;
      rho = rho0;
      q = 0;
      dz = 0.001;
      gam = 1.0;                /* isothermal sheet */

      sigma = sigmasfr = sigma_u4 = 0;

      while(rho > 0.0001 * rho0)
        {
          if(rho > All.DensThreshold)
            {
              tsfr = All.MaxSfrTimescale * sqrt(All.DensThreshold / rho);
            }
          else
            {
              tsfr = 0;
              sigma_u4 += rho * dz;
            }

          P = GAMMA_MINUS1 * rho * u4;

          drho = q;
          dq = -(gam - 2) / rho * q * q - 4 * M_PI * All.G / (gam * P) * rho * rho * rho;

          sigma += rho * dz;
          if(tsfr > 0)
            {
              sigmasfr += rho / tsfr * dz;
            }

          rho += drho * dz;
          q += dq * dz;
        }


      sigma *= 2;               /* to include the other side */
      sigmasfr *= 2;
      sigma_u4 *= 2;

      sigma *= All.HubbleParam * (All.UnitMass_in_g / SOLAR_MASS) * PARSEC * PARSEC / (All.UnitLength_in_cm * All.UnitLength_in_cm);
      sigmasfr *= All.HubbleParam * All.HubbleParam * (All.UnitMass_in_g / SOLAR_MASS) * (SEC_PER_YEAR / All.UnitTime_in_s) * KILOPARSEC * KILOPARSEC / (All.UnitLength_in_cm * All.UnitLength_in_cm);
      sigma_u4 *= All.HubbleParam * (All.UnitMass_in_g / SOLAR_MASS) * PARSEC * PARSEC / (All.UnitLength_in_cm * All.UnitLength_in_cm);


      if(WriteMiscFiles && (ThisTask == 0))
        {
          fprintf(fd, "%g %g %g %g\n", rho0, sigma, sigmasfr, sigma_u4);
        }
    }


  if(All.ComovingIntegrationOn)
    {
      All.Time = All.TimeBegin;
      set_cosmo_factors_for_current_time();
      IonizeParams();
    }

  if(WriteMiscFiles && (ThisTask == 0))
    fclose(fd);
}


double get_ism_h2_frac(int i)
{
  return 1.;
}

#endif


#endif /* closes ISM */
