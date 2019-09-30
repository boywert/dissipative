#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "../allvars.h"
#include "../proto.h"

#ifdef GFM_DUST
#define WARN_TOLERANCE 1.0e-3

/* Assumes set_cosmo_factors_for_current_time() has been called. */
double get_cell_dtime_Gyr(int i)
{
  double time_begstep, dt, dtime, dtime_in_Gyr;

  set_cosmo_factors_for_current_time();

  if(All.ComovingIntegrationOn)
    time_begstep = All.TimeBegin * exp(All.Ti_begstep[P[i].TimeBinGrav] * All.Timebase_interval);
  else
    time_begstep = All.TimeBegin + All.Ti_begstep[P[i].TimeBinGrav] * All.Timebase_interval;

  dt = (P[i].TimeBinGrav ? (((integertime) 1) << P[i].TimeBinGrav) : 0) * All.Timebase_interval;

  if(All.ComovingIntegrationOn)
    dtime = All.Time * dt / All.cf_time_hubble_a;
  else
    dtime = dt;

  dtime *= All.UnitTime_in_s / All.HubbleParam;
  dtime_in_Gyr = dtime / SEC_PER_MEGAYEAR / 1000;
  return dtime_in_Gyr;
}


/* Update the mass of dust in each species for each active cell. */
void dust_growth_and_destruction(void)
{
  int idx, i, k, chan;
  double dm_dust, dtime, mass_metal, mass_dust, dust_chan_fracs[GFM_DUST_N_CHANNELS];
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if((P[i].Type != 0) || (P[i].Mass <= 0.0))
        {
          continue;
        }

      dtime = get_cell_dtime_Gyr(i);
      for(k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
        {
          if((k == element_index_Hydrogen) || (k == element_index_Helium) || (k == element_index_Nitrogen) || (k == element_index_Neon))
            {
              continue;
            }

          mass_metal = SphP[i].MassMetals[k];
          mass_dust = 0.0;
          for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
            {
              mass_dust += SphP[i].MassMetalsDust[chan][k];
            }

          for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
            {
              if(mass_dust > 0.0)
                {
                  dust_chan_fracs[chan] = SphP[i].MassMetalsDust[chan][k] / mass_dust;
                }
              else
                {
                  dust_chan_fracs[chan] = 0.0;
                }
            }

          double f_dust;
          if(mass_dust + mass_metal > 0.0)
            f_dust = mass_dust / (mass_dust + mass_metal);
          else
            f_dust = 1.0;

          double utherm = SphP[i].Utherm;
          double yHe = (1.0 - HYDROGEN_MASSFRAC) / (4.0 * HYDROGEN_MASSFRAC);
          double mu = (1.0 + 4.0 * yHe) / (1.0 + yHe + SphP[i].Ne);
          double temperature_i = (mu * PROTONMASS) / BOLTZMANN * GAMMA_MINUS1 * utherm * All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;

#ifdef USE_SFR
          /* assume that ISM gas is at 10^4 K */
	  if (SphP[i].Sfr > 0.0)
            temperature_i = 1e4;
#endif

          // Kelvin, see Dwek 1998 or Bekki 2015.
          double temperature_0 = 50.0;
          // Internal units corresponding to n_H = 100 cm^-3, assuming mu = 1.4.
          double rho_0 = (1.4 * 1.0e2 * PROTONMASS) / (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam);
          double rho_i = SphP[i].Density * All.cf_a3inv;

          // All.Dust_Growth_Tau in Gyr, since dtime will be in Gyr.
          double a_1 = All.DustSingleGrainSize / 0.1;
#ifndef GFM_DUST_MRN
          double growth_tau = All.Dust_Growth_Tau * a_1 * (rho_0 / rho_i) * pow(temperature_0 / temperature_i, 0.5);
#else
#error "This is not implemented!"
#endif
	  SphP[i].DustTauGrowth = growth_tau;
          dm_dust = dtime * (1.0 - f_dust) * mass_dust / growth_tau;

          if(dm_dust > mass_metal)
            {
              dm_dust = mass_metal;
            }

          for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
            {
              SphP[i].MassMetalsDust[chan][k] += dm_dust * dust_chan_fracs[chan];
            }

          SphP[i].MassMetals[k] -= dm_dust;
          if(SphP[i].MassMetals[k] < 0.0)
            {
              SphP[i].MassMetals[k] = 0.0;
            }
          SphP[i].MassMetallicity -= dm_dust;
          if(SphP[i].MassMetallicity < 0.0)
            {
              SphP[i].MassMetallicity = 0.0;
            }
        }

      for(k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
        {
          mass_metal = SphP[i].MassMetals[k];
          mass_dust = 0.0;
          for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
            {
              mass_dust += SphP[i].MassMetalsDust[chan][k];
            }

          for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
            {
              if(mass_dust > 0.0)
                {
                  dust_chan_fracs[chan] = SphP[i].MassMetalsDust[chan][k] / mass_dust;
                }
              else
                {
                  dust_chan_fracs[chan] = 0.0;
                }
            }

#if (GFM_DUST_DESTMODE==0)
          /* Consult Lisenfeld & Ferrara (1998), Equation 13. */
          double epsilon = 0.3;
          double SNR;
          if(dtime > 0.0)
            {
              SNR = SphP[i].NumSNII / dtime;
            }
          else
            {
              SNR = 0.0;
            }
          /* Calculate shocked gas mass, in internal units. */
          double M_s100 = 1.09 * 6800 * SOLAR_MASS / (All.UnitMass_in_g / All.HubbleParam);
          dm_dust = mass_dust * (epsilon * M_s100) / P[i].Mass * SNR * dtime;

#elif (GFM_DUST_DESTMODE==1)
          dm_dust = dtime * mass_dust / All.Dust_Destruction_Tau;
#endif /* GFM_DUST_DESTMODE */

          if(dm_dust > mass_dust)
            {
              dm_dust = mass_dust;
            }

          for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
            {
              SphP[i].MassMetalsDust[chan][k] -= dm_dust * dust_chan_fracs[chan];
              if(SphP[i].MassMetalsDust[chan][k] < 0.0)
                {
                  SphP[i].MassMetalsDust[chan][k] = 0.0;
                }
            }
          SphP[i].MassMetals[k] += dm_dust;
          SphP[i].MassMetallicity += dm_dust;
        }

#ifdef GFM_DUST_SPUTTERING
#ifdef USE_SFR
      if (SphP[i].Sfr==0.0)
#endif
      {
#if (GFM_DUST_SPUTTERING==0)
      double a_grain = 1.0e-5;  /* 0.1 um */
      double H_s = 4.3 * ELECTRONVOLT_IN_ERGS;  /* 4.3 eV */
      double S_i;
      double amu_masses[] = { 1.0079, 4.0026, 12.0107, 14.0067, 15.9994,
        20.1797, 24.305, 28.0855, 55.845
      };

      /* In these collisions, k is the type of grain, and l is the type */
      /* of gas-phase atom. */
      int l;
      for(k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
        {
          /* No dust for these elements. */
          if((k == element_index_Hydrogen) || (k == element_index_Helium) || (k == element_index_Nitrogen) || (k == element_index_Neon))
            {
              continue;
            }

          for(l = 0; l < GFM_N_CHEM_ELEMENTS; l++)
            {
              mass_dust = 0.0;
              for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
                {
                  mass_dust += SphP[i].MassMetalsDust[chan][k];
                }
              double N_grains = mass_dust / (amu_masses[k] * 1.66e-24 / (All.UnitMass_in_g / All.Hubble));

              for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
                {
                  if(mass_dust > 0.0)
                    {
                      dust_chan_fracs[chan] = SphP[i].MassMetalsDust[chan][k] / mass_dust;
                    }
                  else
                    {
                      dust_chan_fracs[chan] = 0.0;
                    }
                }

              if(l == element_index_Hydrogen)
                {
                  S_i = 4.0e-5;
                }
              else if(l == element_index_Helium)
                {
                  S_i = 8.0e-4;
                }
              else
                {
                  S_i = 8.0e-3;
                }
              double gas_mi = amu_masses[l] * 1.66e-24;
              double gas_ni = SphP[i].Density * All.cf_a3inv * SphP[i].MetalsFraction[l] * (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam) / gas_mi;
              double v_thr = pow(8.0 * H_s / gas_mi, 0.5);

              double numerator =
                (-8.0 * H_s * BOLTZMANN * T + 4.0 * BOLTZMANN * BOLTZMANN * T * T - 4.0 * H_s * gas_mi * v_thr * v_thr + 2.0 * BOLTZMANN * gas_mi * T * v_thr * v_thr +
                 0.5 * gas_mi * gas_mi * v_thr * v_thr * v_thr * v_thr);
              double denominator = exp(0.5 * gas_mi * v_thr * v_thr / (BOLTZMANN * T)) * gas_mi * gas_mi;
              /* Number of particles sputtered per unit time. */
              double rate_sp =
                gas_ni * pow(gas_mi / (2.0 * 3.1415926 * BOLTZMANN * T), 1.5) * 4.0 * (3.1415926 * a_grain * a_grain) * S_i / H_s * BOLTZMANN * T * numerator / denominator * (2.0 * 3.1415926 / 2.0);
              /* add a number of grains parameter */
              double N_sp = rate_sp * N_grains * (dtime * SEC_PER_GIGAYEAR);

              dm_dust = N_sp * (amu_masses[k] * 1.66e-24 / (All.UnitMass_in_g / All.Hubble));

              if(dm_dust > mass_dust)
                {
                  dm_dust = mass_dust;
                }

              for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
                {
                  SphP[i].MassMetalsDust[chan][k] -= dm_dust * dust_chan_fracs[chan];
                  if(SphP[i].MassMetalsDust[chan][k] < 0.0)
                    {
                      SphP[i].MassMetalsDust[chan][k] = 0.0;
                    }
                }
              SphP[i].MassMetals[k] += dm_dust;
              SphP[i].MassMetallicity += dm_dust;
            }
        }
#elif (GFM_DUST_SPUTTERING==1)
      /* See Tsai & Matthews (1995). */
      for(k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
        {
          /* No dust for these elements. */
          if((k == element_index_Hydrogen) || (k == element_index_Helium) || (k == element_index_Nitrogen) || (k == element_index_Neon))
            {
              continue;
            }

          mass_dust = 0.0;
          for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
            {
              mass_dust += SphP[i].MassMetalsDust[chan][k];
            }

          for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
            {
              if(mass_dust > 0.0)
                {
                  dust_chan_fracs[chan] = SphP[i].MassMetalsDust[chan][k] / mass_dust;
                }
              else
                {
                  dust_chan_fracs[chan] = 0.0;
                }
            }


          double tau_sput = gfm_get_dust_thermal_sput_tau(i);
          SphP[i].DustTauSputter = tau_sput;
          dm_dust = mass_dust / tau_sput * dtime;

          if(dm_dust > mass_dust)
            {
              dm_dust = mass_dust;
            }

          for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
            {
              SphP[i].MassMetalsDust[chan][k] -= dm_dust * dust_chan_fracs[chan];
              if(SphP[i].MassMetalsDust[chan][k] < 0.0)
                {
                  SphP[i].MassMetalsDust[chan][k] = 0.0;
                }
            }
          SphP[i].MassMetals[k] += dm_dust;
          SphP[i].MassMetallicity += dm_dust;
        }
#endif
      }
#endif

      /* Reset number of supernovae to zero so that we can determine how */
      /* many occur in the next timestep for this cell. */
      SphP[i].NumSNII = 0.0;
    }
}

double gfm_get_dust_thermal_sput_tau(int i)
{
  double rho_cgs = SphP[i].Density * All.cf_a3inv * (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam);
  double rho_27  = rho_cgs / 1e-27;
  double T_tilde = 2.0e6;       /* K */
  double utherm = SphP[i].Utherm;
  double yHe = (1.0 - HYDROGEN_MASSFRAC) / (4.0 * HYDROGEN_MASSFRAC);
  double mu = (1.0 + 4.0 * yHe) / (1.0 + yHe + SphP[i].Ne);
  double T = (mu * PROTONMASS) / BOLTZMANN * GAMMA_MINUS1 * utherm  * All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;
#ifndef GFM_DUST_MRN
  double a_1 = All.DustSingleGrainSize / 0.1; //All.DustSingleGrainSize in mu
  return All.Dust_Sputter_Tau_Fac * 0.17 * (a_1 / rho_27) * (1.0 + pow(T_tilde / T, 2.5)) / 3; //Gyr
#else
  return All.Dust_Sputter_Tau_Fac * (0.0029 / rho_27) * (1.0 + pow(T_tilde / T, 2.5)) / 3; //Gyr
#endif

}


#endif /* GFM_DUST */
