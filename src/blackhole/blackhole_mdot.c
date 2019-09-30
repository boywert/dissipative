/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/blackhole/blackhole_mdot.c
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


#ifdef BLACK_HOLES

static int blackhole_get_mode(double bhmass, double mdot, double mdot_adios, double meddington);
static double blackhole_continous_mode_switch(double mdot, double mdot_adios, double meddington);

void blackhole_calculate_mdot(void)
{
  int idx, n, flag;
  double mdot, meddington, dt, accreted_mass, mdot_adios = 0;
  double fraction_kinetic, fraction_thermal;

  for(idx = 0; idx < TimeBinsBHAccretion.NActiveParticles; idx++)
    {
      n = TimeBinsBHAccretion.ActiveParticleList[idx];
      if(n < 0)
        continue;

      /* Eddington rate */
      BPP(n).BH_MdotEddington = meddington = blackhole_mdot_eddington(BPP(n).BH_Mass);

      flag = 0;
      mdot = 0.0;
      accreted_mass = 0.0;

#ifdef BH_BONDI_DEFAULT
      BPP(n).BH_MdotBondi = mdot = blackhole_bondi_rate(n, 0);
      flag++;
#endif

#ifdef BH_BONDI_DENSITY
      BPP(n).BH_MdotBondi = mdot = blackhole_bondi_rate(n, 1);
      flag++;
#endif

#ifdef BH_BONDI_CAPTURE
      BPP(n).BH_MdotBondi = mdot = blackhole_bondi_rate(n, 2);
      flag++;
#endif

#ifdef BH_BONDI_DISK_VORTICITY
      BPP(n).BH_MdotBondi = mdot = blackhole_bondi_rate(n, 3);
      flag++;
#endif

      if(flag != 1)
        terminate("BLACK_HOLES: accretion type error");

#ifdef BH_ADIOS_WIND
      mdot_adios = mdot / All.BlackHoleAccretionFactor;
#endif

      int agn_mode = blackhole_get_mode(BPP(n).BH_Mass, mdot, mdot_adios, meddington);

#ifdef BH_CONTINOUS_MODE_SWITCH
      fraction_kinetic = blackhole_continous_mode_switch(mdot, mdot_adios, meddington);
      fraction_thermal = 1.0 - fraction_kinetic;

#else
      if(agn_mode == BH_QUASAR_MODE)
        {
          fraction_kinetic = 0.0;
          fraction_thermal = 1.0;

#ifdef BH_SPIN_EVOLUTION
#if (BH_SPIN_MODEL == 2)
          BPP(n).BH_SpinModel = 1;
#endif
#endif
        }
      else
        {
          fraction_kinetic = 1.0;
          fraction_thermal = 0.0;
#ifdef BH_SPIN_EVOLUTION
#if (BH_SPIN_MODEL == 2)
          BPP(n).BH_SpinModel = 0;
#endif
#endif
        }
#endif


#ifdef BH_ADIOS_WIND
      if(agn_mode == BH_RADIO_MODE)
        mdot = mdot_adios;
#endif


#ifdef BH_PRESSURE_CRITERION
      if(agn_mode == BH_QUASAR_MODE)
        {
          double press_thresh = All.Ref_BH_Pressure * pow(BPP(n).BH_Mass / All.Ref_BH_Mass, 1.0);
          double tot_physical_press = BPP(n).BH_Pressure * All.cf_a3inv;

#ifdef BH_USE_ALFVEN_SPEED_IN_BONDI
          /* convert magnetic pressure to physical (a^{-4}) */
          tot_physical_press += BPP(n).BH_Bpress * All.cf_a3inv / All.Time;
#endif

          if(tot_physical_press < press_thresh)
            mdot *= pow(tot_physical_press / press_thresh, 2);
        }
#endif


#ifdef BH_EXACT_INTEGRATION
      double lambda = 0.0, eta = 0.0, bonditime = 0.0;

      /* Check whether black hole has been merged */
      if(P[n].ID != 0 && P[n].Mass > 0)
        {
          lambda = mdot / BPP(n).BH_Mass / BPP(n).BH_Mass;
          eta = All.BlackHoleEddingtonFactor * meddington / BPP(n).BH_Mass;
          bonditime = (eta - BPP(n).BH_Mass * lambda) / (BPP(n).BH_Mass * lambda * eta);
        }
#endif

      if(mdot > All.BlackHoleEddingtonFactor * meddington)
        mdot = All.BlackHoleEddingtonFactor * meddington;

      dt = (P[n].TimeBinHydro ? (((integertime) 1) << P[n].TimeBinHydro) : 0) * All.Timebase_interval / All.cf_hubble_a;

#ifdef BH_EXACT_INTEGRATION
      if(bonditime >= 0.0)
        {
          if(dt > bonditime)
            {
              /* first phase Bondi accretion, second phase Eddington-limited accretion */
              accreted_mass = BPP(n).BH_Mass * BPP(n).BH_Mass * lambda * bonditime / (1.0 - BPP(n).BH_Mass * lambda * bonditime);
              accreted_mass += ((BPP(n).BH_Mass + (1. - All.BlackHoleRadiativeEfficiency) * accreted_mass) * (exp(eta * (dt - bonditime)) - 1.0));
            }
          else
            accreted_mass = BPP(n).BH_Mass * BPP(n).BH_Mass * lambda * dt / (1.0 - BPP(n).BH_Mass * lambda * dt);
        }
      else
        accreted_mass = BPP(n).BH_Mass * (exp(eta * dt) - 1.0);

      /* reset mdot as average accretion rate over the time-step */
      if(dt > 0)
        mdot = accreted_mass / dt;
      else
        mdot = 0;
#else
      accreted_mass = mdot * dt;
#endif


      BPP(n).BH_Mdot = mdot;

#ifdef BH_NF_RADIO
      BPP(n).BH_Mdot_quasar = mdot;
#endif

      if(!gsl_finite(mdot))
        terminate("mdot=%g", mdot);

      double deltaM = (1. - All.BlackHoleRadiativeEfficiency) * accreted_mass;

#ifdef BH_BIPOLAR_FEEDBACK
      deltaM *= (1. - All.BHBipolarEfficiency);
#endif

      BPP(n).BH_Mass += deltaM;

      if(fraction_thermal > 0.0)
        {
          BPP(n).BH_CumMass_QM += deltaM;

          double egy = fraction_thermal * All.BlackHoleFeedbackFactor * All.BlackHoleRadiativeEfficiency * accreted_mass * pow(CLIGHT / All.UnitVelocity_in_cm_per_s, 2);
#ifdef BH_THERMALFEEDBACK
          BPP(n).BH_ThermEnergy += egy;
          BPP(n).BH_CumEgy_QM += egy;
#endif
          AGNEnergyT_Should += egy;
        }

      if(fraction_kinetic > 0.0)
        {
          BPP(n).BH_CumMass_RM += deltaM;

#ifdef BH_ADIOS_WIND

#ifdef BH_ADIOS_DENS_DEP_EFFICIANCY
          double rad_efficiency = All.BlackHoleRadiativeEfficiency *
                (BPP(n).BH_Density * All.cf_a3inv / (All.RadioFeedbackFactorRefDensityFactor * All.PhysDensThresh)) *
                pow(BPP(n).BH_Mass / All.RadioFeedbackFactorPivotMass, All.RadioFeedbackFactorSlope);

          if(rad_efficiency > All.RadioFeedbackFactorMaxEfficiency)
            rad_efficiency = All.RadioFeedbackFactorMaxEfficiency;
#else
          double rad_efficiency = All.BlackHoleRadiativeEfficiency;
#endif

#if defined(BH_ADIOS_ONLY_ABOVE_MINIMUM_DENSITY)
          if(agn_mode == BH_RADIO_MODE)
            if((BPP(n).BH_Density * All.cf_a3inv) < All.RadioFeedbackMinDensityFactor * All.PhysDensThresh)
              {
                rad_efficiency *= pow((BPP(n).BH_Density * All.cf_a3inv) / (All.RadioFeedbackMinDensityFactor * All.PhysDensThresh), 1);
              }
#endif

          double egy = fraction_kinetic * All.RadioFeedbackFactor * rad_efficiency * accreted_mass * pow(CLIGHT / All.UnitVelocity_in_cm_per_s, 2);

          BPP(n).BH_WindEnergy += egy;
          BPP(n).BH_CumEgy_RM += egy;
#endif

#ifdef BH_BUBBLES
          BPP(n).BH_CumMass_RM += deltaM;
          BPP(n).BH_Mass_bubbles += fraction_kinetic * (1. - All.BlackHoleRadiativeEfficiency) * accreted_mass;
#endif
        }

    }
}





double blackhole_bondi_rate(int n, int type)
{
  double rho, soundspeed, bhvel, mdot;
  int flag;

  rho = BPP(n).BH_Density;

#ifndef BH_USE_ALFVEN_SPEED_IN_BONDI
  soundspeed = sqrt(GAMMA * GAMMA_MINUS1 * BPP(n).BH_U);
#else
  soundspeed = sqrt(GAMMA * GAMMA_MINUS1 * BPP(n).BH_U + 2.0 * BPP(n).BH_Bpress / rho / All.cf_atime);
#endif

#ifdef DRAINGAS
#if (DRAINGAS == 2)
  rho = BPP(n).CellDensity;
#endif
#endif

  flag = 0;
  mdot = 0.0;

#ifdef BH_USE_GASVEL_IN_BONDI
  bhvel = sqrt(pow(P[n].Vel[0] - BPP(n).BH_SurroundingGasVel[0], 2) + pow(P[n].Vel[1] - BPP(n).BH_SurroundingGasVel[1], 2) + pow(P[n].Vel[2] - BPP(n).BH_SurroundingGasVel[2], 2));
#else
  bhvel = 0;
#endif

  /* change to physical variables in comoving integration */
  bhvel /= All.cf_atime;
  rho *= All.cf_a3inv;

  /* default Bondi */
#ifdef BH_BONDI_DEFAULT
  if(type == 0)
    {
      if(rho > 0)
        mdot = 4. * M_PI * All.BlackHoleAccretionFactor * All.G * All.G * BPP(n).BH_Mass * BPP(n).BH_Mass * rho / pow((pow(soundspeed, 2) + pow(bhvel, 2)), 1.5);
      flag++;
    }
#endif

  /* density Bondi */
#ifdef BH_BONDI_DENSITY
  if(type == 1)
    {
      double nH = HYDROGEN_MASSFRAC * rho / PROTONMASS * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
      double nH_thresh = All.PhysDensThresh / (PROTONMASS / HYDROGEN_MASSFRAC / All.UnitDensity_in_cgs);
      double alpha_accretion;

      if(nH < nH_thresh)
        alpha_accretion = 1.0;
      else
        alpha_accretion = pow(nH / nH_thresh, All.BlackHoleAccretionSlope);

      if(rho > 0)
        mdot = 4. * M_PI * alpha_accretion * All.G * All.G * BPP(n).BH_Mass * BPP(n).BH_Mass * rho / pow((pow(soundspeed, 2) + pow(bhvel, 2)), 1.5);

      flag++;
    }
#endif

  /* capture Bondi */
#ifdef BH_BONDI_CAPTURE
  if(type == 2)
    {
      double sigma;
      sigma = pow(All.G * (BPP(n).BH_CaptureMass + BPP(n).BH_Mass) / BPP(n).BH_Hsml, 0.5);

      if(rho > 0 && sigma > 0)
        mdot =
          4. * M_PI * All.BlackHoleAccretionFactor * All.G * All.G * (BPP(n).BH_CaptureMass + BPP(n).BH_Mass) * (BPP(n).BH_CaptureMass +
                                                                                                                 BPP(n).BH_Mass) * rho / pow((pow(soundspeed, 2) + pow(sigma, 2)), 1.5);
      flag++;
    }
#endif

  /* vorticity Bondi */
#ifdef BH_BONDI_DISK_VORTICITY
  if(type == 3)
    {
      double omega_star = 0, omega_crit = 0, r_bondi = 0;

      if(rho > 0)
        mdot = 4. * M_PI * All.BlackHoleAccretionFactor * All.G * All.G * BPP(n).BH_Mass * BPP(n).BH_Mass * rho / pow((pow(soundspeed, 2) + pow(bhvel, 2)), 1.5);

      omega_crit = pow(2, 0.5) * soundspeed / (CLIGHT / All.UnitVelocity_in_cm_per_s);

      omega_star = sqrt(pow(BPP(n).BH_GasVort[0], 2) + pow(BPP(n).BH_GasVort[1], 2) + pow(BPP(n).BH_GasVort[2], 2));

      r_bondi =
        50.0 * (PARSEC / All.UnitLength_in_cm) * ((BPP(n).BH_Mass * All.UnitMass_in_g / All.HubbleParam) / (1.e7 * SOLAR_MASS)) *
        pow(((soundspeed * All.UnitVelocity_in_cm_per_s) / (30.0 * 1.e5)), -2);

      omega_star *= r_bondi / soundspeed;       /* omega_star = |vorticity| * r_bondi / c_snd */

      if(omega_star > omega_crit)
        {
          double factor_omega;
          factor_omega = 0;

          if(omega_star < 0.1)
            factor_omega = 0.34;

          if(omega_star >= 0.1)
            factor_omega = 0.34 * 2.0 * log(16.0 * omega_star) / (3.0 * M_PI * omega_star);

          if(!gsl_finite(factor_omega) || factor_omega <= 0)
            terminate("BH_BONDI_DISK_VORTICITY: factor_omega wrong factor_omega=%g\n", factor_omega);

          if(rho > 0)
            mdot *= factor_omega;
        }
      flag++;
    }
#endif

  if(BPP(n).BH_Mass > 0)
    fprintf(FdBlackHolesDetails, "BH=%llu %g %g %g %g %g\n", (long long) P[n].ID, All.Time, BPP(n).BH_Mass, mdot, rho, soundspeed);


#ifdef BH_BIPOLAR_FEEDBACK
  if(BPP(n).BH_Mass > 0) {
    MyDouble *bj = BPP(n).BH_BipolarJ;
    fprintf(FdBlackHolesBipolar, "BH=%llu %g %g %g %g %g %d %d %g\n", (long long) P[n].ID, All.Time, bj[0],bj[1],bj[2],BPP(n).BH_BipolarSum,(BPP(n).BH_BipolarColdFraction>All.BHBipolarColdFraction),BPP(n).BH_BipolarColdDisk, BPP(n).BH_BipolarColdFraction);
  }
#endif



  if(flag != 1)
    terminate("BLACK_HOLES: accretion type error");

  return mdot;
}


double blackhole_mdot_eddington(double bh_mass)
{
  return (4 * M_PI * GRAVITY * CLIGHT * PROTONMASS / (All.BlackHoleRadiativeEfficiency * CLIGHT * CLIGHT * THOMPSON)) * bh_mass * All.UnitTime_in_s / All.HubbleParam;
}

double blackhole_luminosity_eddington(double bh_mass)
{
  return (4 * M_PI * GRAVITY * CLIGHT * PROTONMASS / (CLIGHT * CLIGHT * THOMPSON)) * bh_mass * All.UnitTime_in_s / All.HubbleParam;
}


static int blackhole_get_mode(double bh_mass, double mdot, double mdot_adios, double meddington)
{
#if defined(UNIFIED_FEEDBACK)
  if(mdot < All.RadioThreshold * meddington)
    return BH_RADIO_MODE;
#endif

#if defined(BH_ADIOS_WIND)
#ifdef BH_ADIOS_WIND_WITH_QUASARTHRESHOLD

  double thresh = All.QuasarThreshold;

#ifdef BH_ADIOS_WIND_WITH_VARIABLE_QUASARTHRESHOLD
  thresh *= pow(bh_mass / (1.0e8 * SOLAR_MASS / (All.UnitMass_in_g / All.HubbleParam)), 2);
  if(thresh > 0.1)
    thresh = 0.1;
#endif

  if(mdot / All.BlackHoleAccretionFactor < thresh * meddington)
    return BH_RADIO_MODE;
#else
  if(All.RadioFeedbackFactor * mdot_adios > All.BlackHoleFeedbackFactor * mdot)
    return BH_RADIO_MODE;
#endif
#endif

  return BH_QUASAR_MODE;
}

#ifdef BH_CONTINOUS_MODE_SWITCH
double blackhole_continous_mode_switch(double mdot, double mdot_adios, double meddington)
{
  double mdot_max, f_k;
  if(All.BlackHoleFeedbackFactor * mdot > All.RadioFeedbackFactor * mdot_adios)
    mdot_max = mdot;
  else
    mdot_max = mdot_adios;
  if(mdot_max > meddington)
    mdot_max = meddington;
  f_k = pow(All.RadioFeedbackFactor * mdot_adios, 4.0) / (pow(All.BlackHoleFeedbackFactor * mdot, 4.0) + pow(All.RadioFeedbackFactor * mdot_adios, 4.0));
  f_k *= (1.0 - pow(mdot_max / meddington, 4.0));
  return f_k;
}
#endif

#ifdef BH_NF_RADIO

double blackhole_get_mdot_radio_from_radiolum(int n)
{
  double mdot = BPP(n).BH_RadioLum / (All.BlackHoleRadiativeEfficiency * pow(CLIGHT / All.UnitVelocity_in_cm_per_s, 2));

  double meddington = blackhole_mdot_eddington(BPP(n).BH_Mass);

  if(mdot > All.BlackHoleEddingtonFactor * meddington)
    mdot = All.BlackHoleEddingtonFactor * meddington;

  return mdot;
}

void blackhole_calculate_mdot_radiomode(void)
{
  int idx, n;

  for(idx = 0; idx < TimeBinsBHAccretion.NActiveParticles; idx++)
    {
      n = TimeBinsBHAccretion.ActiveParticleList[idx];
      if(n < 0)
        continue;

      BPP(n).BH_Mdot_radio = blackhole_get_mdot_radio_from_radiolum(n);

      double dt = (P[n].TimeBinHydro ? (((integertime) 1) << P[n].TimeBinHydro) : 0) * All.Timebase_interval / All.cf_hubble_a;

      /* update the total accretion rate */
      BPP(n).BH_Mdot += BPP(n).BH_Mdot_radio;
      BPP(n).BH_RadioEgyFeedback += All.BlackHoleRadiativeEfficiency * BPP(n).BH_Mdot_radio * dt * pow(CLIGHT / All.UnitVelocity_in_cm_per_s, 2);

      BPP(n).BH_Mass += (1. - All.BlackHoleRadiativeEfficiency) * BPP(n).BH_Mdot_radio * dt;
      BPP(n).BH_CumMass_RM += (1. - All.BlackHoleRadiativeEfficiency) * BPP(n).BH_Mdot_radio * dt;
    }
}

#endif



#endif
