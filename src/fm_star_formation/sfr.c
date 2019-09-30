/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/fm_star_formation/sfr.c
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

#ifdef COMPUTE_SFR_FROM_H2
#include "H2_frac_proto.h"
#endif

#ifdef FM_SFR

#ifndef USE_SFR
#error "FM_SFR requires USE_SFR"
#endif

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

  if(RestartFlag != 1)
    {
#ifdef REFINEMENT
      avg_gas_mass = All.TargetGasMassFactor * All.ReferenceGasPartMass;
      c_sound = sqrt(GAMMA * BOLTZMANN * 600. / (mol_weight * PROTONMASS)) / All.UnitVelocity_in_cm_per_s;
/*    c_sound = 8. * GAMMA * GAMMA_MINUS1 * All.G * avg_gas_mass / All.SofteningGas;
      jeans_density = pow(3 * avg_gas_mass / (4.0 * M_PI * All.SofteningGas * All.SofteningGas * All.SofteningGas), 1.0 / 3.0);
      jeans_mass = pow(All.G, -1.5) * pow(M_PI, 2.5) * pow(c_sound, 1.5) / sqrt(jeans_density) / 6.0;

      mpi_printf("Avg density at softening %g (int units), %g (cm^{-3} h^2)\n", jeans_density,
  	     jeans_density * All.UnitDensity_in_cgs / (mol_weight * PROTONMASS));
      mpi_printf("Resolved Jeans mass %g (int units), avg gas mass %g (int units)\n", jeans_mass, avg_gas_mass);
      mpi_printf("Avg c_sound %g (int units), Jeans \'temperature\' %g\n", sqrt(c_sound),
	     c_sound * pow(All.UnitVelocity_in_cm_per_s, 2.0) * (mol_weight * PROTONMASS) / BOLTZMANN / GAMMA);
*/
      jeans_density = pow(M_PI, 5.0) * pow(c_sound * c_sound / All.G, 3.0) / pow(6 * sqrt(8) * avg_gas_mass, 2.0);
      //jeans_density = M_PI * c_sound * c_sound / All.G / pow(4 * 0.018 * 0.7, 2.0);
      //jeans_mass = jeans_density * pow(0.018 * 0.7, 3.0);
      double jeans_lenght = sqrt(M_PI * c_sound * c_sound / All.G / jeans_density);
      jeans_mass = 4.0 / 3.0 * M_PI * jeans_density * pow(jeans_lenght / 2.0, 3.0);

      mpi_printf("Avg density at softening %g (int units), %g (cm^{-3} h^2)\n", jeans_density, jeans_density * All.UnitDensity_in_cgs / (mol_weight * PROTONMASS));
      mpi_printf("Avg c_sound %g (int units)\n", c_sound);
      mpi_printf("Jeans mass %g (int units), cell mass %g \n", jeans_mass, avg_gas_mass);
      mpi_printf("Jeans lenght %g (int units)\n", jeans_lenght, avg_gas_mass);
      mpi_printf("Cell size %g (int units)\n", pow(3.0 * avg_gas_mass / (4.0 * M_PI * jeans_density), 1.0 / 3.0));
      /* adjusting the critical density such that the jeans mass = 8xavg_gas_mass */
      /* decreasing then this value by a factor of 3 for safety (see also Teyssier+12) */
      if(All.DensThreshold == 0)
        //    All.DensThreshold = jeans_density * pow(8. * avg_gas_mass / jeans_mass, 2.0) * CRIT_DENS_TOLERANCE;
        All.DensThreshold = jeans_density * CRIT_DENS_TOLERANCE;
      else
#endif
        All.DensThreshold = mol_weight * PROTONMASS * All.DensThreshold / All.UnitDensity_in_cgs / All.HubbleParam / All.HubbleParam;
    }

  tff_at_threshold = 0.25 * sqrt(1.5 * M_PI / (All.G * All.DensThreshold));
  All.MaxSfrTimescale = tff_at_threshold / All.SfrEfficiency;

  mpi_printf("density threshold %g (int units), density threshold %g (cm^{-3} h^2)\n", All.DensThreshold, All.DensThreshold * All.UnitDensity_in_cgs / (mol_weight * PROTONMASS));

  mpi_printf("free-fall time at threshold %g (int units), SF efficiency %g\n", tff_at_threshold, All.SfrEfficiency);

  mpi_printf("Max SFR time-scale (int units) % g, SFR timescale %g (Gyr h^-1)\n", All.MaxSfrTimescale, All.MaxSfrTimescale * All.UnitTime_in_s / (1.0e3 * SEC_PER_MEGAYEAR));

#ifdef USE_POLYTROPIC_EQSTATE
  if(RestartFlag != 1)
    {
      mpi_printf("Temperature at threshold %g K, ", All.UthermThreshold);

      All.UthermThreshold = (BOLTZMANN * All.UthermThreshold) / (GAMMA_MINUS1 * mol_weight * PROTONMASS);
      All.UthermThreshold /= (All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s);
    }

  mpi_printf("Utherm at threshold %g (int units)\n", All.UthermThreshold);
#endif

#if defined(FM_STAR_FEEDBACK) && defined(FM_SN_COOLING_RADIUS_BOOST)
  All.one_SNe_energy = All.FeedbackEfficiency * 1.0e51 * (All.HubbleParam / All.UnitEnergy_in_cgs);  /* in code units */
  All.one_SNe_mass   = 10.5*SOLAR_MASS/All.UnitMass_in_g * All.HubbleParam;                          /* in code units */
  All.SNe_velocity   = sqrt(2.0*All.one_SNe_energy/All.one_SNe_mass);    					     /* in code units */
#endif

  integrate_sfr();              /* compute KS law for adopted setting */

  mpi_printf("KS law computed\n");
  mpi_printf("SF module initialized\n\n");
}

/*
 * This routine does cooling and star formation for
 * the FM model.
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


  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Mass == 0 && P[i].ID == 0)
        continue;               /* skip cells that have been swallowed or eliminated */

      flag_sf = 0;

      dens = SphP[i].Density;

#ifdef COMPUTE_SFR_FROM_H2
      fH2 = MolecularHFrac(i);
#else
      fH2 = 1.;
#endif

#ifdef OUTPUT_MOLECULAR_FRACTION
      SphP[i].MolecularFrac = fH2;
#endif

      tsfr = All.MaxSfrTimescale * sqrt(All.DensThreshold / (dens * All.cf_a3inv));
      cloudmass = fH2 * P[i].Mass;

#ifdef DELAYED_COOLING
      if(SphP[i].CoolShutoffTime > 0)
        cloudmass = 0.0;        /* SF is not permitted for a cell where the cooling is disabled */
#endif
#ifdef DELAYED_COOLING_TURB
      if(check_if_turbulent_energy_above_threshold(i))
        cloudmass = 0.0;        /* SF is not permitted for a cell where the cooling is disabled */
#endif

#ifdef FM_RADIATION_FEEDBACK
      if(SphP[i].GasRadCoolShutoffTime > 0)
        cloudmass = 0.0;
#endif

      /* SF only if gas density above threshold */
      SphP[i].Sfr = 0.0;

      if(dens * All.cf_a3inv >= All.DensThreshold)
        flag_sf = 1;

      if(All.ComovingIntegrationOn)     /* to protect against SF at too high redshift */
        if(dens < All.OverDensThresh)
          flag_sf = 0;

      cool_cell(i);

      if(flag_sf == 1)
        {
          SphP[i].Sfr = cloudmass / tsfr;

#ifdef USE_POLYTROPIC_EQSTATE
#ifdef DELAYED_COOLING_TURB
          double umin = All.UthermThreshold * pow((dens * All.cf_a3inv) / All.DensThreshold, 1.0 / 3.0);
#else
          double umin = All.UthermThreshold * pow((dens * All.cf_a3inv) / All.DensThreshold, GAMMA_MINUS1);
#endif
          double unew = dmax(SphP[i].Utherm, umin);
          double du = unew - SphP[i].Utherm;
          SphP[i].Utherm += du;
          SphP[i].Energy += All.cf_atime * All.cf_atime * du * P[i].Mass;
          set_pressure_of_cell(i);
#endif
          SphP[i].Sfr *= (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);
          TimeBinSfr[P[i].TimeBinHydro] += SphP[i].Sfr;
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

  dens = SphP[i].Density;
  tsfr = All.MaxSfrTimescale * sqrt(All.DensThreshold / (dens * All.cf_a3inv));

#ifdef COMPUTE_SFR_FROM_H2
  fH2 = MolecularHFrac(i);
#endif
#ifdef DELAYED_COOLING
  if(SphP[i].CoolShutoffTime > 0)
    cloudmass = 0.0;            /* SF is not permitted for a cell where the cooling is disabled */
#endif
#ifdef DELAYED_COOLING_TURB
  if(check_if_turbulent_energy_above_threshold(i))
    cloudmass = 0.0;            /* SF is not permitted for a cell where the cooling is disabled */
#endif

  /* star fomation only above threshold */
  if(dens * All.cf_a3inv >= All.DensThreshold)
    rateOfSF = fH2 * cloudmass / tsfr;

  if(All.ComovingIntegrationOn) /* to protect against SF at too high redshift */
    if(dens < All.OverDensThresh)
      rateOfSF = 0.0;

#ifdef FM_MOLEC_COOLING
  MyFloat ne    = SphP[i].Ne;
  MyFloat temp  = convert_u_to_temp(SphP[i].Utherm * All.UnitPressure_in_cgs / All.UnitDensity_in_cgs,
                                    SphP[i].Density * All.cf_a3inv * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam,
                                    &ne);
  if(temp > 1e3) rateOfSF *= 0.0;
#endif
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

#ifdef OUTPUT_SF_PROBABILITY
double compute_sf_probability(int index)
{
  double prob = 0, p, mass_of_star;
  double dt, dtime;

  if(P[index].Type == 0)
    {
      if(P[index].Mass == 0 && P[index].ID == 0)
        return 0.0;             /* skip cells that have been swallowed or eliminated */

      dt = (P[index].TimeBin ? (((integertime) 1) << P[index].TimeBin) : 0) * All.Timebase_interval;

      /* the actual time-step */
      dtime = All.cf_atime * dt / All.cf_time_hubble_a;

      mass_of_star = 0;
      p = 0;

      if(SphP[index].Sfr > 0)
        {
          p = SphP[index].Sfr / ((All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR)) * dtime / P[index].Mass;

#if defined(REFINEMENT_SPLIT_CELLS) && defined(REFINEMENT_MERGE_CELLS)
          if(P[index].Mass < 2.0 * All.TargetGasMass)
#ifdef SFR_KEEP_CELLS
            mass_of_star = 0.9 * P[index].Mass;
#else
            mass_of_star = P[index].Mass;
#endif
          else
            mass_of_star = All.TargetGasMass;

#ifdef REFINEMENT_HIGH_RES_GAS
          if(SphP[index].HighResMass < HIGHRESMASSFAC * P[index].Mass)
            {
              /* this cell does not appear to be in the high-res region.
                 If we form a star, then it is given the mass of the cell,
                 and later we give the star the SofteningType=3 particle to give it large softening */
#ifdef SFR_KEEP_CELLS
              mass_of_star = 0.9 * P[index].Mass;
#else
              mass_of_star = P[index].Mass;
#endif
            }

#endif /* REFINEMENT_SPLIT_CELLS && REFINEMENT_MERGE_CELLS */

#else
          mass_of_star = P[index].Mass;
#endif

#ifdef SFR_KEEP_CELLS
          if(P[index].Mass < 0.5 * All.TargetGasMass)
            continue;           /* do not make stars from cells that should be derefined */
#endif
          prob = P[index].Mass / mass_of_star * (1 - exp(-p));
        }

      if(prob < 0)
        terminate("prob < 0");

      if(prob > 1)
        {
          printf
            ("SFR: Warning, need to make a heavier star than desired. Task=%d prob=%g P[index].Mass=%g mass_of_star=%g mass_of_star_new=%g p=%g\n",
             ThisTask, prob, P[index].Mass, mass_of_star, P[index].Mass * (1 - exp(-p)), p);
          mass_of_star = P[index].Mass * (1 - exp(-p));
          prob = 1.0;
        }
    }

  return prob;
}
#endif /* closes OUTPUT_SF_PROBABILITY */

#endif /* closes FM_SFR */
