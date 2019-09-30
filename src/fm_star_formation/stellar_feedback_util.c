/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/fm_star_formation/stellar_feedback_util.c
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

#include "../allvars.h"
#include "../proto.h"


#ifdef FM_STAR_FEEDBACK

#ifndef FM_SFR
#warning FM_STAR_FEEDBACK requires FM_SFR
#endif

#define SN_ENERGY   1.0e51      /* erg */
#define CM_PER_KM   1.0e5
#define SN_MOMENTUM (1.0e10 * SOLAR_MASS)      /* gr cm/s */

/*! \file fm_stellar_feedback.c
 *  \brief routines for stellar feedback
 */
MyDouble compute_SN_energy(MyDouble Number_of_SNae)
{
  MyDouble Injected_Energy = All.FeedbackEfficiency * Number_of_SNae * SN_ENERGY;       /* a fraction FeedbackEfficiency is given to the gas */
  Injected_Energy *= (All.HubbleParam / All.UnitEnergy_in_cgs); /* to code units */

  return Injected_Energy;
}

#ifdef DIRECT_MOMENTUM_INJECTION_FEEDBACK
MyDouble compute_SN_momentum(MyDouble Number_of_SNae)
{
  MyDouble Injected_Momentum = All.MomentumFactor * Number_of_SNae * SN_MOMENTUM;       
  Injected_Momentum *= (All.HubbleParam / All.UnitVelocity_in_cm_per_s / All.UnitMass_in_g); /* to code units */

  return Injected_Momentum;
}
#endif

/*  \param b = <e_r, P_i>, 
 *  \param c = 2[M\DeltaE + \deltaM E_f] = 2J
 *  E_f is the final kinetic energy \DeltaE the feedback energy
 *  \return the magnitude of the injected momentum in radial direction
 *          such that the kinetic energy of the cell is increased by \DeltaE
 */
MyDouble quadratic_equation_solver(MyDouble b, MyDouble c)
{
  MyDouble res;

  if(b <= 0.0)
    res = -(b - sqrt(b * b + c));
  else
    res = c / (b + sqrt(b * b + c));

  return res;
}

#ifdef DELAYED_COOLING

void set_cooling_shutoff_time(int i, MyFloat CoolShutoffTime)
{
#if (SHUTOFFTIME_UPDATE == 0)
  SphP[i].CoolShutoffTime = dmax(SphP[i].CoolShutoffTime, CoolShutoffTime);
#endif
#if (SHUTOFFTIME_UPDATE == 1)
  SphP[i].CoolShutoffTime += CoolShutoffTime;
#endif
}

/* compute the SN blast radius, following Chevalier (1974) */
MyFloat compute_blast_radius(MyFloat energy, MyFloat avg_H_density, MyFloat avg_pressure)
{
  /* assumes energy in 10^{51} erg avg_H_density in cgs, avg_pressure in 10^{-4} K cm^{-3} */
  energy *= (All.UnitEnergy_in_cgs / All.HubbleParam / 1.0e51);
  avg_H_density *= (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam);
  avg_pressure *= (All.UnitPressure_in_cgs * All.HubbleParam * All.HubbleParam / BOLTZMANN * 1.0e-4);

  MyFloat blast_radius = pow(10., 1.74) * pow(energy, 0.32) * pow(avg_H_density, -0.16) * pow(avg_pressure, -0.2);
  blast_radius *= (PARSEC * All.HubbleParam / All.UnitLength_in_cm);

  return blast_radius;
}

MyFloat compute_cooling_shutoff_time(MyFloat energy, MyFloat avg_H_density, MyFloat avg_pressure)
{
  /* assumes energy in 10^{51} erg avg_H_density in cgs, avg_pressure in 10^{-4} K cm^{-3} */
  energy *= (All.UnitEnergy_in_cgs / All.HubbleParam / 1.0e51);
  avg_H_density *= (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam);
  avg_pressure *= (All.UnitPressure_in_cgs * All.HubbleParam * All.HubbleParam / BOLTZMANN * 1.0e-4);

  MyFloat shutoff_time = pow(10., 6.85) * pow(energy, 0.32) * pow(avg_H_density, 0.34) * pow(avg_pressure, -0.7);
  shutoff_time *= (SEC_PER_YEAR * All.HubbleParam / All.UnitTime_in_s);

  return shutoff_time;
}

/* This function returns 0 if the cell cannot cool, otherwise returns the time for which cooling is allowed */
MyFloat update_cooling_shutoff_time(int i, MyFloat dt)
{
  MyFloat dtcool = 0.0;

  SphP[i].CoolShutoffTime -= dt;        /* particle has been evolved for a time step */

  if(SphP[i].CoolShutoffTime < 0.0)     /* cooling can happen again */
    {
      dtcool = -SphP[i].CoolShutoffTime;        /* part of the time step for which cooling is allowed */
      SphP[i].CoolShutoffTime = 0;      /* reset the value so that cell can undergo normal cooling */
    }

  return dtcool;
}

#endif /* DELAYED_COOLING */

#ifdef DELAYED_COOLING_TURB

void init_turbulent_energy()
{
  mpi_printf("Initializing delayed cooling module\n");

  if(RestartFlag != 1)
    {
      /* converting the dissipation time scale from Megayears to code units */
      All.DissipationTime *= (SEC_PER_MEGAYEAR / All.UnitTime_in_s * All.HubbleParam);

      /* converting sigma threshold from km/s to code units */
      All.SigmaThreshold *= (CM_PER_KM / All.UnitVelocity_in_cm_per_s);

      mpi_printf("Dissipation timescale %g (int units), Sigma threshold %g (int units)\n", All.DissipationTime, All.SigmaThreshold);

      /* from sigma threshold to turbulent energy threshold */
      All.SigmaThreshold *= (0.5 * All.SigmaThreshold);
    }

  mpi_printf("Turbulent energy threshold %g (int units)\n", All.SigmaThreshold);

  mpi_printf("Delayed cooling with turbulent energy advection initialized\n\n");
}

void convert_specific_turbulent_energy_to_turbulent_energy()
{
  int i;
  mpi_printf(" Converting specific turbulent energy to turbulent energy... \n");

  for(i = 0; i < NumGas; i++)
    SphP[i].MassUturb = SphP[i].Uturb * P[i].Mass;
}

int check_if_turbulent_energy_above_threshold(int i)
{
  return (SphP[i].Uturb > All.SigmaThreshold ? 1 : 0);
}

MyFloat dissipate_turbulent_energy(MyFloat Uturb_old, MyFloat dt)
{
  return (Uturb_old * exp(-dt / All.DissipationTime));
}
#endif

#ifdef INSTANTANEOUS_DEPOSITION
void init_instantaneous_deposition()
{
  mpi_printf("Initializing instantaneous deposition module\n");

  if(RestartFlag != 1)
    {
      /* converting stellar age from Megayears to Gyr */
      All.MinFeedAge *= 1.0e-3;
      All.SNaePerUnitMass = integrate_imf(log10(All.SNII_MinMass_Msun), log10(All.SNII_MaxMass_Msun), INTEGRATE_IMF_NUMBER, NULL);
    }

  mpi_printf("Stellar age after feedback energy is injected %g (Gyr)\n", All.MinFeedAge);
  mpi_printf("Number of type II supernovae per unit mass %g\n", All.SNaePerUnitMass);

  mpi_printf("Instantaneous deposition module initialized\n\n");
}

MyDouble compute_total_SN_energy(int index, MyDouble Age, MyFloat * SN_Number)
{
  MyDouble Injected_Energy = 0.0;

  /* initialize the number of SN events for feedback */
  *SN_Number = 0.0;

  /* instantaneosuly deposit all the type II SN energy */
  if(Age >= All.MinFeedAge && STP(index).FeedbackFlag > 0)
    {
      *SN_Number = STP(index).InitialMass * All.UnitMass_in_g * All.SNaePerUnitMass / SOLAR_MASS / All.HubbleParam / All.FeedbackInjectionEvents;
      Injected_Energy = *SN_Number * All.FeedbackEfficiency * SN_ENERGY;        /* a fraction FeedbackEfficiency is given to the gas */
      Injected_Energy *= (All.HubbleParam / All.UnitEnergy_in_cgs);     /* to code units */
      STP(index).FeedbackFlag--;        /* star is allowed to inject feedback only for FeedbackInjectionEvents times */
    }

  return Injected_Energy;
}
#endif

#endif
