/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sfr_eEOS.c
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
#include "allvars.h"
#include "proto.h"
#include "forcetree.h"

/** \file sfr_eEOS.c
 *  \brief Star formation rate routines for the effective multi-phase model.
 */

#ifdef USE_SFR

#if !defined(FM_SFR) && !defined(ISM) && !defined(LOCAL_FEEDBACK)
#if !defined(QUICK_LYALPHA)

/** \brief Main driver for star formation and gas cooling.
 *
 *  This function loops over all the active gas cells. If a given cell
 *  meets the criteria for star formation to be active the multi-phase
 *  model is activated, the properties of the cell are updated according to
 *  the latter and the star formation rate computed. In the other case, the
 *  standard isochoric cooling is applied to the gas cell by calling the function
 *  cool_cell() and the star formation rate is set to 0.
 */
void cooling_and_starformation(void)
{
  TIMER_START(CPU_COOLINGSFR);

  int idx, i, bin, flag;
  double dt, dtime, ne = 1;
  double unew, du;
  double cloudmass;
  double factorEVP, dens;
  double tsfr;
  double egyeff, x;

#ifdef SLOW_RELAX_TO_EOS
  double trelax, tcool;
#endif

#ifdef BH_THERMALFEEDBACK
  double temp;
#endif

  double eos_dens_threshold = All.PhysDensThresh;
#ifdef MODIFIED_EOS
  eos_dens_threshold *= All.FactorDensThresh;
#endif

  /* note: assuming FULL ionization */
  double u_to_temp_fac = (4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC))) * PROTONMASS / BOLTZMANN * GAMMA_MINUS1 * All.UnitEnergy_in_cgs / All.UnitMass_in_g;

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


      dens = SphP[i].Density;

      dt = (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;
      dtime = All.cf_atime * dt / All.cf_time_hubble_a;

      /* apply the temperature floor (and, if relevant, BH thermal feedback) */

      unew = dmax(All.MinEgySpec, SphP[i].Utherm);

#ifdef BH_THERMALFEEDBACK
      if(SphP[i].Injected_BH_Energy)
        {
          if(SphP[i].Injected_BH_Energy < 0)
            terminate("strange feedback energy: Thistask=%d i=%d SphP[i].Injected_BH_Energy=%g\n", ThisTask, i, SphP[i].Injected_BH_Energy);

          unew += SphP[i].Injected_BH_Energy / P[i].Mass;

          temp = u_to_temp_fac * unew;

          if(temp > 5.0e9)
            unew = 5.0e9 / u_to_temp_fac;

          AGNEnergyT_Is += SphP[i].Injected_BH_Energy;
          SphP[i].Injected_BH_Energy = 0;
        }
#endif

      if(unew < 0)
        terminate("Invalid Temperature: Task=%d i=%d unew=%g\n", ThisTask, i, unew);

      du = unew - SphP[i].Utherm;
      SphP[i].Utherm += du;
      SphP[i].Energy += All.cf_atime * All.cf_atime * du * P[i].Mass;

      egyeff = 0.;
      /* calculate the effective equation of state for gas above the density threshold */
      if(dens * All.cf_a3inv >= eos_dens_threshold)
        {
          ne = SphP[i].Ne;
          egyeff = calc_egyeff(i, dens * All.cf_a3inv, &ne, &x, &tsfr, &factorEVP);
        }

      /* do cooling, except for gas above the EOS density threshold that is colder than the eEOS */
      if(dens * All.cf_a3inv < eos_dens_threshold || (dens * All.cf_a3inv >= eos_dens_threshold && SphP[i].Utherm > egyeff))
        {
#ifdef GFM_LAMBDA
         cooling_set_partindex(i);
#endif

#ifdef SLOW_RELAX_TO_EOS
          tcool = GetCoolingTime(SphP[i].Utherm, dens * All.cf_a3inv, &ne);
#else
          GetCoolingTime(SphP[i].Utherm, dens * All.cf_a3inv, &ne);
#endif
          cool_cell(i);
        }
#ifdef SLOW_RELAX_TO_EOS
      else
        tcool = 0;
#endif

      /* check whether conditions for star formation are fulfilled.
       * f=1  normal cooling
       * f=0  star formation
       */

      flag = 1;                 /* default is normal cooling */

      /* enable star formation if gas is above SF density threshold */
      if(dens * All.cf_a3inv >= eos_dens_threshold)
        if(SphP[i].Utherm <= egyeff || u_to_temp_fac * SphP[i].Utherm <= All.TemperatureThresh)
          flag = 0;

      if(All.ComovingIntegrationOn)
        if(dens < All.OverDensThresh)
          flag = 1;

      if(P[i].Mass == 0)        /* tracer particles don't form stars */
        flag = 1;

      if(flag == 1)
        SphP[i].Sfr = 0;

      /* active star formation */
      if(flag == 0)
        {
          SphP[i].Ne = (HYDROGEN_MASSFRAC + 1) / 2 / HYDROGEN_MASSFRAC; /* note: assuming FULL ionization */

          cloudmass = x * P[i].Mass;

          if(tsfr < dtime)
            tsfr = dtime;

          if(dt > 0)
            {
              if(P[i].TimeBinHydro)     /* upon start-up, we need to protect against dt==0 */
                {
                  unew = SphP[i].Utherm;

                  if(SphP[i].Utherm < egyeff)
                    {
#ifdef SLOW_RELAX_TO_EOS
                      trelax = tsfr * (1 - x) / x / (All.FactorSN * (1 + factorEVP));
                      if(tcool < trelax && tcool > 0)
                        trelax = tcool;
                      unew = (egyeff + (SphP[i].Utherm - egyeff) * exp(-dtime / trelax));
#else
                      unew = egyeff;
#endif
                    }

                  du = unew - SphP[i].Utherm;
                  if(unew < All.MinEgySpec)
                    du = All.MinEgySpec - SphP[i].Utherm;

                  SphP[i].Utherm += du;
                  SphP[i].Energy += All.cf_atime * All.cf_atime * du * P[i].Mass;

#ifdef OUTPUT_COOLHEAT
                  if(dtime > 0)
                    SphP[i].CoolHeat = du * P[i].Mass / dtime;
#endif

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
                  SphP[i].A = (GAMMA - 1.0) * SphP[i].Utherm / pow(dens * All.cf_a3inv, GAMMA - 1);
                  SphP[i].Entropy = log(SphP[i].A) * P[i].Mass;
#endif
                  set_pressure_of_cell(i);
                }
            }

#if defined(GFM_STELLAR_EVOLUTION) && (GFM_STELLAR_EVOLUTION == 0)
          SphP[i].Sfr = cloudmass / tsfr * (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);
#else
          SphP[i].Sfr = (1 - All.FactorSN) * cloudmass / tsfr * (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);
#endif
#ifdef MODIFIED_EOS
          /* prevent star formation if gas is not dense enough */
          if(dens * All.cf_a3inv < All.PhysDensThresh)
            SphP[i].Sfr = 0;
#endif

#if defined(REFINEMENT_AROUND_BH) && defined(SUPPRESS_SF_IN_REFINEMENT_REGION)
        if(SphP[i].RefBHFlag) {
          SphP[i].Sfr = 0;
        }
#endif

          TimeBinSfr[P[i].TimeBinHydro] += SphP[i].Sfr;
        }
    }                           /* end of main loop over active particles */

  TIMER_STOP(CPU_COOLINGSFR);
}

/** \brief Return the star formation rate associated with the gas cell i.
 *
 *  \param i the index of the gas cell
 *  \return star formation rate in solar masses / yr
 */
double get_starformation_rate(int i)
{
  if(RestartFlag == 3)
    return SphP[i].Sfr;

  double rateOfSF;
  int flag;
  double tsfr;
  double factorEVP, egyeff, ne, x, cloudmass;
  /* note: assuming FULL ionization */
  double u_to_temp_fac = (4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC))) * PROTONMASS / BOLTZMANN * GAMMA_MINUS1 * All.UnitEnergy_in_cgs / All.UnitMass_in_g;

  double eos_dens_threshold = All.PhysDensThresh;
#ifdef MODIFIED_EOS
  eos_dens_threshold *= All.FactorDensThresh;
#endif

  flag = 1;                     /* default is normal cooling */
  egyeff = 0.0;

  if(SphP[i].Density * All.cf_a3inv >= eos_dens_threshold)
    {
      ne = SphP[i].Ne;
      egyeff = calc_egyeff(i, SphP[i].Density * All.cf_a3inv, &ne, &x, &tsfr, &factorEVP);
    }

  if(SphP[i].Density * All.cf_a3inv >= All.PhysDensThresh)
    if(SphP[i].Utherm <= 1.01 * egyeff || u_to_temp_fac * SphP[i].Utherm <= All.TemperatureThresh)
      flag = 0;

  if(All.ComovingIntegrationOn)
    if(SphP[i].Density < All.OverDensThresh)
      flag = 1;

  if(flag == 1)
    return 0;

  cloudmass = x * P[i].Mass;

#if defined(DUST_LIVE) && defined(DL_GRAIN_BINS) && (defined(DL_SNE_DESTRUCTION) || defined(DL_SHATTERING) || defined(DL_COAGULATION))
  SphP[i].CloudFrac = x;
#endif

#if defined(GFM_STELLAR_EVOLUTION) && (GFM_STELLAR_EVOLUTION == 0)
  rateOfSF = cloudmass / tsfr;
#else
  rateOfSF = (1 - All.FactorSN) * cloudmass / tsfr;
#endif

  /* convert to solar masses per yr */
  rateOfSF *= (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

  return rateOfSF;
}

#endif

/** \brief Initialize the parameters of effective multi-phase model.
 *
 *   In particular this function computes the value of PhysDensThresh, that is
 *   the physical density threshold above which star formation is active, if its
 *   value was set to 0 in the parameter file.
 */
void init_clouds(void)
{
  double A0, dens, tcool, ne, coolrate, egyhot, x, u4, meanweight;
  double tsfr, peff, fac, neff, egyeff, factorEVP, sigma, thresholdStarburst;

#ifdef SOFTEREQS
  meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));   /* note: assuming FULL ionization */
  All.UForSofterEQS = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.TempForSofterEQS;
  All.UForSofterEQS *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;       /* unit conversion */
#endif

  if(All.PhysDensThresh == 0)
    {
      A0 = All.FactorEVP;

      egyhot = All.EgySpecSN / A0;

      meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));       /* note: assuming FULL ionization */
      u4 = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * 1.0e4;
      u4 *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

      /* choose a high reference density to avoid that we pick up a compton cooling contribution */
      if(All.ComovingIntegrationOn)
        dens = 1.0e10 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
      else
        dens = 1.0e10 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

      if(All.ComovingIntegrationOn)
        {
          All.Time = 1.0;       /* to be guaranteed to get z=0 rate */
          set_cosmo_factors_for_current_time();
          IonizeParams();
        }

      ne = 1.0;
      SetZeroIonization();

#ifdef GFM_COOLING_METAL
      update_gas_state(dens, GFM_INITIAL_ABUNDANCE_HYDROGEN, 0.0);
#endif
      tcool = GetCoolingTime(egyhot, dens, &ne);

      coolrate = egyhot / tcool / dens;

      x = (egyhot - u4) / (egyhot - All.EgySpecCold);

      All.PhysDensThresh = x / pow(1 - x, 2) * (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold) / (All.MaxSfrTimescale * coolrate);

      mpi_printf
        ("USE_SFR: A0=%g   PhysDensThresh=%g (int units) %g h^2 cm^-3   expected fraction of cold gas at threshold=%g   tcool=%g   dens=%g   egyhot=%g\n",
         A0, All.PhysDensThresh, All.PhysDensThresh / (PROTONMASS / HYDROGEN_MASSFRAC / All.UnitDensity_in_cgs), x, tcool, dens, egyhot);

#ifdef MODIFIED_EOS
      All.UthermAtThresh = effective_eos(All.PhysDensThresh);
      All.JoinDens = invert_effective_eos(All.UthermAtThresh * All.FactorUthermJoin);
#endif

      dens = All.PhysDensThresh;

      do
        {
          ne = 0.5;
          egyeff = calc_egyeff(-1, dens, &ne, &x, &tsfr, &factorEVP);
          peff = GAMMA_MINUS1 * dens * egyeff;

          fac = 1 / (log(dens * 1.025) - log(dens));
          dens *= 1.025;

          neff = -log(peff) * fac;

          ne = 0.5;
          egyeff = calc_egyeff(-1, dens, &ne, &x, &tsfr, &factorEVP);
          peff = GAMMA_MINUS1 * dens * egyeff;

          neff += log(peff) * fac;
        }
      while(neff > 4.0 / 3);

      thresholdStarburst = dens;

#ifdef STEEPER_SFR_FOR_STARBURST
      All.PhysDensThreshStarburst = thresholdStarburst;
#endif

      mpi_printf("USE_SFR: run-away sets in for dens=%g   dynamic range for quiescent star formation=%g\n", thresholdStarburst, thresholdStarburst / All.PhysDensThresh);

      integrate_sfr();

      if(ThisTask == 0)
        {
          sigma = 10.0 / All.Hubble * 1.0e-10 / pow(1.0e-3, 2);

          printf("USE_SFR: isotherm sheet central density=%g   z0=%g\n", M_PI * All.G * sigma * sigma / (2 * GAMMA_MINUS1) / u4, GAMMA_MINUS1 * u4 / (2 * M_PI * All.G * sigma));
          myflush(stdout);

        }

#if defined(GFM_STELLAR_EVOLUTION) && (GFM_STELLAR_EVOLUTION == 0)
      mpi_printf("USE_SFR: SNII energy=%g [internal units] = %g [erg/M_sun] = %g [1e51 erg/Msun]\n", All.FactorSN * All.EgySpecSN,
                 All.FactorSN * All.EgySpecSN / (All.UnitMass_in_g / All.UnitEnergy_in_cgs) * SOLAR_MASS,
                 All.FactorSN * All.EgySpecSN / (All.UnitMass_in_g / All.UnitEnergy_in_cgs) * SOLAR_MASS / 1e51);
#else
      mpi_printf("USE_SFR: SNII energy=%g [internal units] = %g [erg/M_sun] = %g [1e51 erg/Msun]\n", All.FactorSN * All.EgySpecSN,
                 All.FactorSN * All.EgySpecSN / (1 - All.FactorSN) / (All.UnitMass_in_g / All.UnitEnergy_in_cgs) * SOLAR_MASS,
                 All.FactorSN * All.EgySpecSN / (1 - All.FactorSN) / (All.UnitMass_in_g / All.UnitEnergy_in_cgs) * SOLAR_MASS / 1e51);
#endif

#if defined(BH_PRESSURE_CRITERION)
#ifdef BH_BONDI_DENSITY
#error "BH_PRESSURE_CRITERION and BH_BONDI_DENSITY do not work together"
#endif
      /* choose a high density to avoid that we pick up a compton cooling contribution */
      dens = 1.0e10 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

      meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));       /* note: assuming FULL ionization */

      FILE *fd;
      if(WriteMiscFiles && (ThisTask == 0))
        {
          char buf[1000];
          sprintf(buf, "%s/bh_pressure_threshold.txt", All.OutputDir);
          fd = fopen(buf, "w");
        }
      else
        fd = 0;

      double u = u4 * 1.0e5;

      if(All.ComovingIntegrationOn)
        {
          All.Time = 1.0;       /* to be guaranteed to get z=0 rate */
          set_cosmo_factors_for_current_time();
          IonizeParams();
        }

      All.Ref_BH_Pressure = 0;
      All.Ref_BH_Mass = 0;

      while(u >= u4)
        {
          double temp = u * meanweight * GAMMA_MINUS1 * PROTONMASS / BOLTZMANN * All.UnitEnergy_in_cgs / All.UnitMass_in_g;

#ifdef GFM_COOLING_METAL
          update_gas_state(dens, GFM_INITIAL_ABUNDANCE_HYDROGEN, 0.0);
#endif
          double ne = 1.0;
          double tcool = GetCoolingTime(u, dens, &ne);

          double mgas = All.DesNumNgbBlackHole * All.ReferenceGasPartMass;

          double mbh = sqrt(u * mgas * pow(GAMMA * GAMMA_MINUS1 * u, 1.5) /
                            (dens * tcool * All.BlackHoleFeedbackFactor * All.BlackHoleAccretionFactor * All.BlackHoleRadiativeEfficiency * 4 * M_PI *
                             pow(CLIGHT / All.UnitVelocity_in_cm_per_s, 2) * All.G * All.G));

          double ref_press = GAMMA_MINUS1 * All.PhysDensThresh * u;
#ifdef MODIFIED_EOS
          ref_press *= All.FactorUthermAtThresh;
#endif

          if(All.Ref_BH_Pressure == 0)
            {
              All.Ref_BH_Pressure = ref_press;
              All.Ref_BH_Mass = mbh;
	      mpi_printf("BLACKHOLES:   Ref_BH_Pressure=%g   Ref_BH_Mass=%g\n", All.Ref_BH_Pressure, All.Ref_BH_Mass);
            }

          double fit_press = All.Ref_BH_Pressure * pow(mbh / All.Ref_BH_Mass, 1.0);

          if(fd)
            fprintf(fd, "%g %g %g %g %g %g\n", u, temp, mbh, ref_press, fit_press, u / (tcool * dens));

          u *= 0.95;
        }

      if(fd)
        fclose(fd);
#endif


#ifdef BH_NF_RADIO
      /* choose a high density to avoid that we pick up a compton cooling contribution */
      FILE *fdrad;
      if(WriteMiscFiles && (ThisTask == 0))
        {
          char buf[1000];
          sprintf(buf, "%s/radio_mode_eff_vs_vvir.txt", All.OutputDir);
          fdrad = fopen(buf, "w");
        }
      else
        fdrad = 0;

      double vcirc;

      for(vcirc = 20.0; vcirc <= 2000.0; vcirc *= 1.02)
        {
          double R = blackhole_get_radio_efficiency(vcirc);

          if(fdrad)
            fprintf(fdrad, "%g %g\n", vcirc, R);
        }

      if(fdrad)
        fclose(fdrad);
#endif

      if(All.ComovingIntegrationOn)
        {
          All.Time = All.TimeBegin;
          set_cosmo_factors_for_current_time();
          IonizeParams();
        }

    }
}

/** \brief Compute the effective equation of state for the gas and
 *         the integrated SFR per unit area.
 *
 *  This function computes the effective equation of state for the gas and
 *  the integrated SFR per unit area. It saves the results into two files:
 *  eos.txt for the equation of state and sfrrate.txt for the integrated SFR.
 *  In the latter case, the SFR is determined by integrating along the vertical
 *  direction the gas density of an infinite self-gravitating isothermal sheet.
 *  The integrated gas density is saved as well, so effectively sfrrate.txt
 *  contains the Kennicutt-Schmidt law of the star formation model.
 */
void integrate_sfr(void)
{
  double rho0, rho, rho2, q, dz, gam, sigma = 0, sigma_u4, sigmasfr = 0, ne, P1;
  double x = 0, P, P2, x2, tsfr2, factorEVP2, drho, dq;
  double meanweight, u4, tsfr, factorEVP, egyeff, egyeff2;
  FILE *fd;

  double eos_dens_threshold = All.PhysDensThresh;
#ifdef MODIFIED_EOS
  eos_dens_threshold *= All.FactorDensThresh;
#endif

  meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));   /* note: assuming FULL ionization */
  u4 = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * 1.0e4;
  u4 *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

  if(All.ComovingIntegrationOn)
    {
      All.Time = 1.0;           /* to be guaranteed to get z=0 rate */
      set_cosmo_factors_for_current_time();
      IonizeParams();
    }

  if(WriteMiscFiles && (ThisTask == 0))
    fd = fopen("eos.txt", "w");
  else
    fd = 0;

  for(rho = eos_dens_threshold; rho <= 1000 * eos_dens_threshold; rho *= 1.1)
    {
      ne = 1.0;
      egyeff = calc_egyeff(-1, rho, &ne, &x, &tsfr, &factorEVP);

      P = GAMMA_MINUS1 * rho * egyeff;

      if(WriteMiscFiles && (ThisTask == 0))
        {
          fprintf(fd, "%g %g %g\n", rho, P, x);
        }
    }

  if(WriteMiscFiles && (ThisTask == 0))
    fclose(fd);


  if(WriteMiscFiles && (ThisTask == 0))
    fd = fopen("sfrrate.txt", "w");
  else
    fd = 0;

  for(rho0 = eos_dens_threshold; rho0 <= 10000 * eos_dens_threshold; rho0 *= 1.02)
    {
      rho = rho0;
      q = 0;
      dz = 0.001;

      sigma = sigmasfr = sigma_u4 = 0;

      while(rho > 0.0001 * rho0)
        {
          if(rho > All.PhysDensThresh)
            {
              ne = 1.0;
              egyeff = calc_egyeff(-1, rho, &ne, &x, &tsfr, &factorEVP);

              P = P1 = GAMMA_MINUS1 * rho * egyeff;

              rho2 = 1.1 * rho;

              egyeff2 = calc_egyeff(-1, rho2, &ne, &x2, &tsfr2, &factorEVP2);

              P2 = GAMMA_MINUS1 * rho2 * egyeff2;

              gam = log(P2 / P1) / log(rho2 / rho);
            }
          else
            {
              tsfr = 0;

              P = GAMMA_MINUS1 * rho * u4;
              gam = 1.0;

              sigma_u4 += rho * dz;
            }


          drho = q;
          dq = -(gam - 2) / rho * q * q - 4 * M_PI * All.G / (gam * P) * rho * rho * rho;

          sigma += rho * dz;
          if(tsfr > 0)
            {
#if defined(GFM_STELLAR_EVOLUTION) && (GFM_STELLAR_EVOLUTION == 0)
              sigmasfr += rho * x / tsfr * dz;
#else
              sigmasfr += (1 - All.FactorSN) * rho * x / tsfr * dz;
#endif
            }

          rho += drho * dz;
          q += dq * dz;
        }


      sigma *= 2;               /* to include the other side */
      sigmasfr *= 2;
      sigma_u4 *= 2;

      sigma *= All.HubbleParam * (All.UnitMass_in_g / SOLAR_MASS) * PARSEC * PARSEC / (All.UnitLength_in_cm * All.UnitLength_in_cm);
      sigmasfr *= All.HubbleParam * All.HubbleParam * (All.UnitMass_in_g / SOLAR_MASS) * (SEC_PER_YEAR / All.UnitTime_in_s) * 1.0e6 * PARSEC * PARSEC / (All.UnitLength_in_cm * All.UnitLength_in_cm);
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

/** \brief Set the appropriate units for the parameters of the multi-phase model.
 */
void set_units_sfr(void)
{
  double meanweight;

  All.OverDensThresh = All.CritOverDensity * All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

  All.PhysDensThresh = All.CritPhysDensity * PROTONMASS / HYDROGEN_MASSFRAC / All.UnitDensity_in_cgs;

  meanweight = 4 / (1 + 3 * HYDROGEN_MASSFRAC); /* note: assuming NEUTRAL GAS */

  All.EgySpecCold = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.TempClouds;
  All.EgySpecCold *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

  meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));   /* note: assuming FULL ionization */

  All.EgySpecSN = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.TempSupernova;
  All.EgySpecSN *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;
}

/** \brief Calculate the effective energy of the multi-phase model
 */
double calc_egyeff(int i, double gasdens, double *ne, double *x, double *tsfr, double *factorEVP)
{
  double egyhot, egyeff, tcool, y;
  double rho = gasdens;

#ifndef MODIFIED_EOS
  if(rho < All.PhysDensThresh)
    terminate("Error in gas cell %d, Effective equation of state can be computed only for gas above SF density threshold", i);
#else
  rho = dmax(rho, All.PhysDensThresh);
#endif

  *tsfr = sqrt(All.PhysDensThresh / rho) * All.MaxSfrTimescale;

  *factorEVP = pow(rho / All.PhysDensThresh, -0.8) * All.FactorEVP;

  egyhot = All.EgySpecSN / (1 + *factorEVP) + All.EgySpecCold;

#if defined(GFM_AGN_RADIATION) || defined(GFM_UVB_CORRECTIONS) || defined(UVB_SELF_SHIELDING)
  update_radiation_state(rho, HYDROGEN_MASSFRAC, 0);
#endif

#ifdef GFM_COOLING_METAL
  update_gas_state(rho, GFM_INITIAL_ABUNDANCE_HYDROGEN, 0.0);
#endif

  if(i >= 0)
    {
#ifndef OTVET_COOLING_HEATING
      tcool = GetCoolingTime(egyhot, rho, ne);
#else
      tcool = otvet_GetCoolingTime(i, egyhot, rho, ne);
#endif
    }
  else
    tcool = GetCoolingTime(egyhot, rho, ne);

  y = *tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);

  *x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));

  egyeff = egyhot * (1 - *x) + All.EgySpecCold * (*x);

#ifdef MODIFIED_EOS
  /* modify the value of the effective energy to make the EOS harder, recompute the mass fraction of cold gas *x accordingly */
  if(gasdens < All.PhysDensThresh)
    egyeff = egyeff * (1.0 + (All.FactorUthermAtThresh - 1.0) * (gasdens - All.PhysDensThresh * All.FactorDensThresh) / (All.PhysDensThresh * (1.0 - All.FactorDensThresh)));
  else if(gasdens < All.JoinDens)
    egyeff = All.FactorUthermAtThresh * All.UthermAtThresh * (1.0 + (All.FactorUthermJoin / All.FactorUthermAtThresh - 1.0) * (gasdens - All.PhysDensThresh) / (All.JoinDens - All.PhysDensThresh));

  *x = (egyeff - egyhot) / (All.EgySpecCold - egyhot);

  if(*x < 0.0)
    terminate("Not possible");
#endif

#ifdef SOFTEREQS
  /* use an intermediate EQS, between isothermal and the full multiphase model */
  egyeff = All.FactorForSofterEQS * egyeff + (1 - All.FactorForSofterEQS) * All.UForSofterEQS;
#endif

#ifdef STEEPER_SFR_FOR_STARBURST
  /* modify sfr at high densities, but don't change the energy calculated above */
  if(rho > All.PhysDensThreshStarburst)
    *tsfr = pow(All.PhysDensThreshStarburst / rho, All.StarburstPowerLawIndex ) * sqrt(All.PhysDensThresh / All.PhysDensThreshStarburst) * All.MaxSfrTimescale;
#endif

  return egyeff;
}

#ifdef MODIFIED_EOS
/** \brief Calculate the effective energy of the multi-phase model given the gas density
 */
double effective_eos(double rho)
{
  double egyhot, egyeff, tcool, y, x, ne, tsfr, factorEVP;

  ne = 1.0;

  tsfr = sqrt(All.PhysDensThresh / rho) * All.MaxSfrTimescale;

  factorEVP = pow(rho / All.PhysDensThresh, -0.8) * All.FactorEVP;

  egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

#if defined(GFM_AGN_RADIATION) || defined(GFM_UVB_CORRECTIONS) || defined(UVB_SELF_SHIELDING)
  update_radiation_state(rho, HYDROGEN_MASSFRAC, 0);
#endif

#ifdef GFM_COOLING_METAL
  update_gas_state(rho, GFM_INITIAL_ABUNDANCE_HYDROGEN, 0.0);
#endif

  tcool = GetCoolingTime(egyhot, rho, &ne);

  y = tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);

  x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));

  egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

  return egyeff;
}

/** \brief Calculate the gas density given the effective energy of the multi-phase model
 */
double invert_effective_eos(double egyeff_old)
{
  double drho, rho, rho_up, rho_down, egyeff;
  int iter = 0;

  if(egyeff_old < All.UthermAtThresh)
    terminate("Impossible to invert EOS. egyeff=%g is less than %g the effective energy at the star formation density threshold", egyeff_old, All.UthermAtThresh);

  rho_down = All.PhysDensThresh;
  rho_up = All.PhysDensThresh;

  /* bracketing */
  if(effective_eos(rho_down) - egyeff_old <= 0.0)
    {
      while(effective_eos(rho_up) - egyeff_old < 0.0)
        {
          rho_down = rho_up;
          rho_up *= 1.1;
        }
    }
  else if(effective_eos(rho_up) - egyeff_old > 0.0)
    {
      while(effective_eos(rho_down) - egyeff_old > 0.0)
        {
          rho_up = rho_down;
          rho_down /= 1.1;
        }
    }
  else
    terminate("This should not happen. Not good initial guesses for bisection");

  do
    {
      rho = 0.5 * (rho_down + rho_up);

      egyeff = effective_eos(rho);

      if(egyeff - egyeff_old > 0)
        {
          rho_up = rho;
        }
      else
        {
          rho_down = rho;
        }

      drho = rho_up - rho_down;

      iter++;

      if(iter >= (MAXITER - 10))
        printf("rho= %g\n", rho);
    }
  while(fabs(drho / rho) > 1.0e-6 && iter < MAXITER);

  if(iter >= MAXITER)
    terminate("failed to converge in invert_effective_eos(): egyeff_old=%g\n", egyeff_old);

  return rho;
}

void check_modified_eos_parameters()
{
  if(All.FactorDensThresh > 1.0)
    terminate("Maximum value for FactorDensThresh is 1.0");

  if(All.FactorUthermAtThresh < 1.0)
    terminate("Minimum value for FactorUthermAtThresh is 1.0");

  if(All.FactorUthermJoin < All.FactorUthermAtThresh)
    terminate("Minimum value for FactorUthermJoin is the value of FactorUthermAtThresh. In this case %g", All.FactorUthermAtThresh);
}
#endif

#endif
#endif
