/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/GFM/stellar_evolution_evolve.c
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
#include <mpi.h>
#include <hdf5.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include "../allvars.h"
#include "../proto.h"

#include "cooling_metal_vars.h"

#ifdef GFM_STELLAR_EVOLUTION

static double SNIa_Rate_Norm;

/*!
 * This routine computes yields (per unit solar mass) of AGB stars
*/
void evolve_AGB(MyDouble log_min_mass, MyDouble log_max_mass, MyFloat log_metallicity, MyFloat * initial_metals, stellar_evolution_data * sed)
{
  int ilow, ihigh, imass;
  int iz_low, iz_high, i;
  MyFloat dz, metallicity;
  MyFloat temp_mass;
  double imf_integrand_mass[GFM_N_MASS_BINS];
#ifdef GFM_DUST
  int carbon_index, oxygen_index, other_index;
  MyFloat carbon_ejecta, oxygen_ejecta;
  MyFloat ejecta_CO_ratio[GFM_N_MASS_BINS];
#endif

  if(All.AGB_MassTransferOn == 0)
    return;

  metallicity = pow(10.0, log_metallicity);

  /* determine integration range, limiting to stars that become AGB stars */
  if(log_max_mass > log10(All.SNII_MinMass_Msun))
    log_max_mass = log10(All.SNII_MinMass_Msun);

  if(log_min_mass >= log_max_mass)
    return;

  /* determine which mass bins will contribute */
  get_imf_bins(log_min_mass, log_max_mass, &ilow, &ihigh);

  /* determine yield of these bins */
  get_z_indicies(log_metallicity, yieldsAGB.Metallicity, yieldsAGB.N_Z, &iz_low, &iz_high, &dz);

  /* total ejected metal mass from ejected metal mass table */
  for(imass = ilow; imass < ihigh + 1; imass++)
    imf_integrand_mass[imass] = (1 - dz) * (yieldsAGB.TotalMetals_spline[iz_low][imass] + metallicity * yieldsAGB.Ejecta_spline[iz_low][imass])
      + dz * (yieldsAGB.TotalMetals_spline[iz_high][imass] + metallicity * yieldsAGB.Ejecta_spline[iz_high][imass]);

  /* total mass released */
  temp_mass = integrate_imf(log_min_mass, log_max_mass, INTEGRATE_IMF_YIELD, imf_integrand_mass);

  if(temp_mass >= 0)
    sed->total_metal_mass_released += temp_mass;

  /* total ejected mass from ejected mass table */
  for(imass = ilow; imass < ihigh + 1; imass++)
    imf_integrand_mass[imass] = (1 - dz) * yieldsAGB.Ejecta_spline[iz_low][imass] + dz * yieldsAGB.Ejecta_spline[iz_high][imass];

  double ejected_mass = integrate_imf(log_min_mass, log_max_mass, INTEGRATE_IMF_YIELD, imf_integrand_mass);       /* total mass released */

  sed->total_mass_released += ejected_mass;

  for(i = 0; i < GFM_N_CHEM_ELEMENTS; i++)
    {
      double met_released = initial_metals[i] * ejected_mass;

      sed->metal_mass_released[i] += met_released;

#ifdef GFM_CHEMTAGS
      if (i > 1) // Only metals
	sed->metal_mass_released_chemtags[GFM_AGB_CHEMTAG] += met_released;
#if defined(GFM_SPLITFE) && defined(GFM_SPLITFE_ADDINAGB)// divide negative Fe contribution between 2 tags
      if (i == element_index_Iron)
        {
          sed->metal_mass_released_chemtags[GFM_FESNIA_CHEMTAG] += 0.5 * met_released;
          sed->metal_mass_released_chemtags[GFM_FESNII_CHEMTAG] += 0.5 * met_released;
        }
#endif
#endif

      /* ejected metal mass of newly created metals */
      for(imass = ilow; imass < ihigh + 1; imass++)
        imf_integrand_mass[imass] = (1 - dz) * (yieldsAGB.spline[iz_low][i][imass]) + dz * (yieldsAGB.spline[iz_high][i][imass]);

      met_released = integrate_imf(log_min_mass, log_max_mass, INTEGRATE_IMF_YIELD, imf_integrand_mass);
      sed->metal_mass_released[i] += met_released;

#ifdef GFM_CHEMTAGS
      if (i > 1) // only metals
	sed->metal_mass_released_chemtags[GFM_AGB_CHEMTAG] += met_released;
#if defined(GFM_SPLITFE) && defined(GFM_SPLITFE_ADDINAGB) // divide negative Fe contribution between 2 tags
      if (i == element_index_Iron)
        {
          sed->metal_mass_released_chemtags[GFM_FESNIA_CHEMTAG] += 0.5 * met_released;
          sed->metal_mass_released_chemtags[GFM_FESNII_CHEMTAG] += 0.5 * met_released;
        }
#endif
#endif
    }



#ifdef GFM_DUST
  carbon_index = element_index_Carbon;
  oxygen_index = element_index_Oxygen;
  if((carbon_index < 0) || (oxygen_index < 0))
    {
      terminate("GFM_DUST requires carbon and oxygen tracking");
    }
  /* Need to calculate C/O ratio in the ejecta (from initial metallicity */
  /* plus new metals) for each mass bin. */
  for(imass = ilow; imass < ihigh + 1; imass++)
    {
      carbon_ejecta = ejecta_mass_for_bin(yieldsAGB.Ejecta_spline, yieldsAGB.spline, carbon_index, initial_metals, iz_low, iz_high, imass, dz);
      oxygen_ejecta = ejecta_mass_for_bin(yieldsAGB.Ejecta_spline, yieldsAGB.spline, oxygen_index, initial_metals, iz_low, iz_high, imass, dz);
      if(oxygen_ejecta < 1.0e-10)
        ejecta_CO_ratio[imass] = 1.0e10;
      else
        ejecta_CO_ratio[imass] = carbon_ejecta / oxygen_ejecta;
    }

  for(i = 0; i < GFM_N_CHEM_ELEMENTS; i++)
    {
      for(imass = ilow; imass < ihigh + 1; imass++)
        {
          if(ejecta_CO_ratio[imass] > 1.0)
            {
              if(i == carbon_index)
                {
                  imf_integrand_mass[imass] =
                    All.AGB_Dust_Delta_C *
                    (ejecta_mass_for_bin(yieldsAGB.Ejecta_spline, yieldsAGB.spline, carbon_index, initial_metals, iz_low, iz_high, imass, dz) -
                     0.75 * ejecta_mass_for_bin(yieldsAGB.Ejecta_spline, yieldsAGB.spline, oxygen_index, initial_metals, iz_low, iz_high, imass, dz));
                }
              else
                {
                  imf_integrand_mass[imass] = 0.0;
                }
            }
          else                  /* if (ejecta_CO_ratio[imass] <= 1.0) */
            {
              if((i == element_index_Carbon) || (i == element_index_Hydrogen) || (i == element_index_Helium) || (i == element_index_Nitrogen) || (i == element_index_Neon))
                {
                  imf_integrand_mass[imass] = 0.0;
                }
              else if(i == element_index_Oxygen)
                {
                  imf_integrand_mass[imass] = 0.0;
                  other_index = element_index_Magnesium;
                  imf_integrand_mass[imass] += 10.0 * All.AGB_Dust_Delta_Metal / GFM_DUST_AMU_MG
                    * (ejecta_mass_for_bin(yieldsAGB.Ejecta_spline, yieldsAGB.spline, other_index, initial_metals, iz_low, iz_high, imass, dz));
                  other_index = element_index_Silicon;
                  imf_integrand_mass[imass] += 10.0 * All.AGB_Dust_Delta_Metal / GFM_DUST_AMU_SI
                    * (ejecta_mass_for_bin(yieldsAGB.Ejecta_spline, yieldsAGB.spline, other_index, initial_metals, iz_low, iz_high, imass, dz));
                  other_index = element_index_Iron;
                  imf_integrand_mass[imass] += 10.0 * All.AGB_Dust_Delta_Metal / GFM_DUST_AMU_FE
                    * (ejecta_mass_for_bin(yieldsAGB.Ejecta_spline, yieldsAGB.spline, other_index, initial_metals, iz_low, iz_high, imass, dz));
                }
              else
                {
                  imf_integrand_mass[imass] = All.AGB_Dust_Delta_Metal * ejecta_mass_for_bin(yieldsAGB.Ejecta_spline, yieldsAGB.spline, i, initial_metals, iz_low, iz_high, imass, dz);
                }
            }

          /* Ensure that in each mass bin, no more dust is going to be */
          /* produced than overall metals (i.e. dust and gas-phase). */
          double bin_metals = ejecta_mass_for_bin(yieldsAGB.Ejecta_spline, yieldsAGB.spline, i, initial_metals, iz_low, iz_high, imass, dz);
          if(imf_integrand_mass[imass] > bin_metals)
            {
              imf_integrand_mass[imass] = bin_metals;
            }
        }                       /* mass bin loop */

      double released_dust = integrate_imf(log_min_mass, log_max_mass, INTEGRATE_IMF_YIELD, imf_integrand_mass);
      sed->dust_mass_released[GFM_DUST_AGB][i] += released_dust;
      sed->total_dust_mass_released += released_dust;
      /* With dust enabled, we want the metal variables to only refer to */
      /* gas-phase metals. */
      sed->metal_mass_released[i] -= released_dust;
      sed->total_metal_mass_released -= released_dust;
    }                           /* element loop */
#endif
#if defined(DUST_LIVE) && defined(DL_PRODUCTION)
  for(i = 0; i < GFM_N_CHEM_ELEMENTS; i++)
    {
      for(imass = ilow; imass < ihigh + 1; imass++)
        {
          double dust_lhs = GSD.AGB_CondEff_spline[iz_low][i][imass] * (yieldsAGB.Ejecta_spline[iz_low][imass] * initial_metals[i] + yieldsAGB.spline[iz_low][i][imass]);
          double dust_rhs = GSD.AGB_CondEff_spline[iz_high][i][imass] * (yieldsAGB.Ejecta_spline[iz_high][imass] * initial_metals[i] + yieldsAGB.spline[iz_high][i][imass]);
          imf_integrand_mass[imass] = (1-dz)*dust_lhs + dz*dust_rhs;
        }
      double released_dust = integrate_imf(log_min_mass, log_max_mass, INTEGRATE_IMF_YIELD, imf_integrand_mass);
      sed->dust_mass_released[i] += released_dust;
      sed->total_dust_mass_released += released_dust;
      sed->metal_mass_released[i] -= released_dust;
      sed->total_metal_mass_released -= released_dust;
      sed->total_mass_released -= released_dust;
    }
#endif
}


/*!
 * This routine computes yields from and number of SNII stars
 */
void evolve_SNII(MyDouble log_min_mass, MyDouble log_max_mass, MyFloat log_metallicity, MyFloat * initial_metals, stellar_evolution_data * sed)
{
  int ilow, ihigh, imass;
  int iz_low, iz_high, i;
  MyFloat dz, metallicity;
  MyFloat temp_mass;
  double imf_integrand_mass[GFM_N_MASS_BINS];
#ifdef GFM_DUST
  int other_index;
#endif
  sed->number_of_SNII = 0;

  if(All.SNII_MassTransferOn == 0)
    return;

  metallicity = pow(10.0, log_metallicity);

  /* determine integration range: make sure all these stars actually become SN of type II */
  if(log_min_mass < log10(All.SNII_MinMass_Msun))
    log_min_mass = log10(All.SNII_MinMass_Msun);

  if(log_max_mass > log10(All.SNII_MaxMass_Msun))
    log_max_mass = log10(All.SNII_MaxMass_Msun);

  if(log_min_mass >= log_max_mass)
    return;

  /* determine which mass bins will contribute */
  get_imf_bins(log_min_mass, log_max_mass, &ilow, &ihigh);
  sed->number_of_SNII = integrate_imf(log_min_mass, log_max_mass, INTEGRATE_IMF_NUMBER, NULL);

  /* determine yield of these bins */
  get_z_indicies(log_metallicity, yieldsSNII.Metallicity, yieldsSNII.N_Z, &iz_low, &iz_high, &dz);

  /* total ejected metal mass from ejected metal mass table */
  for(imass = ilow; imass < ihigh + 1; imass++)
    imf_integrand_mass[imass] = (1 - dz) * (yieldsSNII.TotalMetals_spline[iz_low][imass] + metallicity * yieldsSNII.Ejecta_spline[iz_low][imass])
      + dz * (yieldsSNII.TotalMetals_spline[iz_high][imass] + metallicity * yieldsSNII.Ejecta_spline[iz_high][imass]);

  /* total mass released */
  temp_mass = integrate_imf(log_min_mass, log_max_mass, INTEGRATE_IMF_YIELD, imf_integrand_mass);

  if(temp_mass >= 0)
    sed->total_metal_mass_released += temp_mass;

  /* total ejected mass from ejected mass table */
  for(imass = ilow; imass < ihigh + 1; imass++)
    imf_integrand_mass[imass] = (1 - dz) * yieldsSNII.Ejecta_spline[iz_low][imass] + dz * yieldsSNII.Ejecta_spline[iz_high][imass];

  double ejected_mass = integrate_imf(log_min_mass, log_max_mass, INTEGRATE_IMF_YIELD, imf_integrand_mass);       /* total mass released */

  sed->total_mass_released += ejected_mass;

  for(i = 0; i < GFM_N_CHEM_ELEMENTS; i++)
    {
      /* ejected metal mass based on initial metallicity */
      double met_released = initial_metals[i] * ejected_mass;

      sed->metal_mass_released[i] += met_released;

#ifdef GFM_CHEMTAGS
      if (i>1)
	sed->metal_mass_released_chemtags[GFM_SNII_CHEMTAG] += met_released;
#ifdef GFM_SPLITFE
      if(i == element_index_Iron)
	sed->metal_mass_released_chemtags[GFM_FESNII_CHEMTAG] += met_released;
#endif
#endif

      /* ejected metal mass of newly created metals */
      for(imass = ilow; imass < ihigh + 1; imass++)
        imf_integrand_mass[imass] = (1 - dz) * (yieldsSNII.spline[iz_low][i][imass]) + dz * (yieldsSNII.spline[iz_high][i][imass]);

      met_released = integrate_imf(log_min_mass, log_max_mass, INTEGRATE_IMF_YIELD, imf_integrand_mass);

      sed->metal_mass_released[i] += met_released;
#ifdef GFM_CHEMTAGS
      if (i > 1) // only metals
	sed->metal_mass_released_chemtags[GFM_SNII_CHEMTAG] += met_released;
#ifdef GFM_SPLITFE
      if(i == element_index_Iron)
	sed->metal_mass_released_chemtags[GFM_FESNII_CHEMTAG] += met_released;
#endif
#endif
    }


#ifdef GFM_DUST
  for(i = 0; i < GFM_N_CHEM_ELEMENTS; i++)
    {
      for(imass = ilow; imass < ihigh + 1; imass++)
        {
          if(i == element_index_Oxygen)
            {
              imf_integrand_mass[imass] = 0.0;
              other_index = element_index_Magnesium;
              imf_integrand_mass[imass] += 10.0 * All.SN_Dust_Delta_Metal / GFM_DUST_AMU_MG
                * (ejecta_mass_for_bin(yieldsSNII.Ejecta_spline, yieldsSNII.spline, other_index, initial_metals, iz_low, iz_high, imass, dz));
              other_index = element_index_Silicon;
              imf_integrand_mass[imass] += 10.0 * All.SN_Dust_Delta_Metal / GFM_DUST_AMU_SI
                * (ejecta_mass_for_bin(yieldsSNII.Ejecta_spline, yieldsSNII.spline, other_index, initial_metals, iz_low, iz_high, imass, dz));
              other_index = element_index_Iron;
              imf_integrand_mass[imass] += 10.0 * All.SN_Dust_Delta_Metal / GFM_DUST_AMU_FE
                * (ejecta_mass_for_bin(yieldsSNII.Ejecta_spline, yieldsSNII.spline, other_index, initial_metals, iz_low, iz_high, imass, dz));
            }
          else if((i == element_index_Hydrogen) || (i == element_index_Helium) || (i == element_index_Nitrogen) || (i == element_index_Neon))
            imf_integrand_mass[imass] = 0.0;
          else if(i == element_index_Carbon)
            imf_integrand_mass[imass] = All.SN_Dust_Delta_C * ejecta_mass_for_bin(yieldsSNII.Ejecta_spline, yieldsSNII.spline, i, initial_metals, iz_low, iz_high, imass, dz);
          else
            imf_integrand_mass[imass] = All.SN_Dust_Delta_Metal * ejecta_mass_for_bin(yieldsSNII.Ejecta_spline, yieldsSNII.spline, i, initial_metals, iz_low, iz_high, imass, dz);

          /* Ensure that in each mass bin, no more dust is going to be */
          /* produced than overall metals (i.e. dust and gas-phase). */
          double bin_metals = ejecta_mass_for_bin(yieldsSNII.Ejecta_spline, yieldsSNII.spline, i, initial_metals, iz_low, iz_high, imass, dz);
          if(imf_integrand_mass[imass] > bin_metals)
            {
              imf_integrand_mass[imass] = bin_metals;
            }
        }                       /* mass bin loop */

      double released_dust = integrate_imf(log_min_mass, log_max_mass, INTEGRATE_IMF_YIELD, imf_integrand_mass);
      sed->dust_mass_released[GFM_DUST_SNII][i] += released_dust;
      sed->total_dust_mass_released += released_dust;
      /* With dust enabled, we want the metal variables to only refer to */
      /* gas-phase metals. */
      sed->metal_mass_released[i] -= released_dust;
      sed->total_metal_mass_released -= released_dust;
    }                           /* element loop */
#endif
#if defined(DUST_LIVE) && defined(DL_PRODUCTION)
  for(i = 0; i < GFM_N_CHEM_ELEMENTS; i++)
    {
      for(imass = ilow; imass < ihigh + 1; imass++)
        {
          double dust_lhs = GSD.SNII_CondEff[i] * (yieldsSNII.Ejecta_spline[iz_low][imass] * initial_metals[i] + yieldsSNII.spline[iz_low][i][imass]);
          double dust_rhs = GSD.SNII_CondEff[i] * (yieldsSNII.Ejecta_spline[iz_high][imass] * initial_metals[i] + yieldsSNII.spline[iz_high][i][imass]);
          imf_integrand_mass[imass] = (1-dz)*dust_lhs + dz*dust_rhs;
        }
      double released_dust = integrate_imf(log_min_mass, log_max_mass, INTEGRATE_IMF_YIELD, imf_integrand_mass);
      sed->dust_mass_released[i] += released_dust;
      sed->total_dust_mass_released += released_dust;
      sed->metal_mass_released[i] -= released_dust;
      sed->total_metal_mass_released -= released_dust;
      sed->total_mass_released -= released_dust;
    }
#endif
}


MyFloat get_SNIa_number_dtd_efolding(MyFloat age_of_star_in_Gyr, MyFloat dt_in_Gyr)
{
  return SNIa_Rate_Norm * (exp(-age_of_star_in_Gyr / All.SNIa_Rate_TAU) - exp(-(age_of_star_in_Gyr + dt_in_Gyr) / All.SNIa_Rate_TAU));
}

MyFloat get_SNIa_number_dtd_powerlaw(MyFloat age_of_star_in_Gyr, MyFloat dt_in_Gyr)
{
  /* if WDs have not yet formed */
  if(age_of_star_in_Gyr + dt_in_Gyr <= All.SNIa_Rate_TAU)
    return 0.0;

  /* Place lower integration limit at All.SNIa_Rate_TAU Gyrs, and adjust integration timestep (dt) size */
  if(age_of_star_in_Gyr < All.SNIa_Rate_TAU)
    {
      dt_in_Gyr -= (All.SNIa_Rate_TAU - age_of_star_in_Gyr);
      age_of_star_in_Gyr += (All.SNIa_Rate_TAU - age_of_star_in_Gyr);
    }

  return SNIa_Rate_Norm * (pow(age_of_star_in_Gyr / All.SNIa_Rate_TAU, 1. - GFM_SNIA_DTD_POWERLAW_INDEX) - pow((age_of_star_in_Gyr + dt_in_Gyr) / All.SNIa_Rate_TAU, 1. - GFM_SNIA_DTD_POWERLAW_INDEX));
}


/*!
 * This routine computes yields from and number of SNIa stars
 */
void evolve_SNIa(MyDouble log_min_mass, MyDouble log_max_mass, MyFloat age_of_star_in_Gyr, MyFloat dt_in_Gyr, MyFloat metallicity, stellar_evolution_data * sed)
{
  int i;
#ifdef GFM_DUST
  int other_index;
#endif

#ifdef GFM_SNIA_DTD_EFOLDING
  sed->number_of_SNIa = get_SNIa_number_dtd_efolding(age_of_star_in_Gyr, dt_in_Gyr);
#endif

#ifdef GFM_SNIA_DTD_POWERLAW
  sed->number_of_SNIa = get_SNIa_number_dtd_powerlaw(age_of_star_in_Gyr, dt_in_Gyr);
#endif

  if(All.SNIa_MassTransferOn == 1)
    {
      for(i = 0; i < GFM_N_CHEM_ELEMENTS; i++)
	{
          sed->metal_mass_released[i] += sed->number_of_SNIa * yieldsSNIa.spline[i];
#ifdef GFM_CHEMTAGS
	if (i > 1) // only metals
	  sed->metal_mass_released_chemtags[GFM_SNIA_CHEMTAG] += sed->number_of_SNIa * yieldsSNIa.spline[i];
#endif
	}

      sed->total_mass_released += sed->number_of_SNIa * (yieldsSNIa.TotalMetals_spline + yieldsSNIa.spline[element_index_Hydrogen] + yieldsSNIa.spline[element_index_Helium]);

      sed->total_metal_mass_released += sed->number_of_SNIa * yieldsSNIa.TotalMetals_spline;

#ifdef GFM_SPLITFE
      sed->metal_mass_released_chemtags[GFM_FESNIA_CHEMTAG] += sed->number_of_SNIa * yieldsSNIa.spline[element_index_Iron];
#endif

#ifdef GFM_DUST
      for(i = 0; i < GFM_N_CHEM_ELEMENTS; i++)
        {
          double released_dust = 0.0;

          if(i == element_index_Oxygen)
            {
              other_index = element_index_Magnesium;
              released_dust += 10.0 * All.SN_Dust_Delta_Metal / GFM_DUST_AMU_MG * (sed->number_of_SNIa) * yieldsSNIa.spline[other_index];
              other_index = element_index_Silicon;
              released_dust += 10.0 * All.SN_Dust_Delta_Metal / GFM_DUST_AMU_SI * (sed->number_of_SNIa) * yieldsSNIa.spline[other_index];
              other_index = element_index_Iron;
              released_dust += 10.0 * All.SN_Dust_Delta_Metal / GFM_DUST_AMU_FE * (sed->number_of_SNIa) * yieldsSNIa.spline[other_index];
            }
          else if((i == element_index_Hydrogen) || (i == element_index_Helium) || (i == element_index_Nitrogen) || (i == element_index_Neon))
            released_dust += 0.0;
          else if(i == element_index_Carbon)
            released_dust += All.SN_Dust_Delta_C * (sed->number_of_SNIa) * yieldsSNIa.spline[i];
          else
            released_dust += All.SN_Dust_Delta_Metal * (sed->number_of_SNIa) * yieldsSNIa.spline[i];

          if(released_dust > (sed->number_of_SNIa * yieldsSNIa.spline[i]))
            {
              released_dust = sed->number_of_SNIa * yieldsSNIa.spline[i];
            }
          sed->dust_mass_released[GFM_DUST_SNIa][i] += released_dust;
          sed->total_dust_mass_released += released_dust;
          /* With dust enabled, we want the metal variables to only refer to */
          /* gas-phase metals. */
          sed->metal_mass_released[i] -= released_dust;
          sed->total_metal_mass_released -= released_dust;
        }
#endif
#if defined(DUST_LIVE) && defined(DL_PRODUCTION)
      for(i = 0; i < GFM_N_CHEM_ELEMENTS; i++)
        {
          double released_dust = GSD.SNII_CondEff[i] * sed->number_of_SNIa * yieldsSNIa.spline[i];
          sed->dust_mass_released[i] += released_dust;
          sed->total_dust_mass_released += released_dust;
          sed->metal_mass_released[i] -= released_dust;
          sed->total_metal_mass_released -= released_dust;
          sed->total_mass_released -= released_dust;
        }
#endif
    }
}

#ifdef GFM_RPROCESS
// Get NSNS merger number
MyFloat get_NSNS_number_dtd_powerlaw(MyFloat age_of_star_in_Gyr, MyFloat dt_in_Gyr)
{

  /* if NSNS have not yet formed */
  if(age_of_star_in_Gyr + dt_in_Gyr <= All.NSNS_Rate_TAU)
    return 0.0;

  /* Place lower integration limit at All.NSNS_Rate_TAU Gyrs, and adjust integration timestep (dt) size */
  if(age_of_star_in_Gyr < All.NSNS_Rate_TAU)
    {
      dt_in_Gyr -= (All.NSNS_Rate_TAU - age_of_star_in_Gyr);
      age_of_star_in_Gyr += (All.NSNS_Rate_TAU - age_of_star_in_Gyr);
    }

  // calculate integrated DTD => total number of NSNS mergers that have gone off
  return SNIa_Rate_Norm * All.NSNS_per_SNIa * (pow(age_of_star_in_Gyr / All.NSNS_Rate_TAU, 1. - GFM_NSNS_DTD_POWERLAW_INDEX) - pow((age_of_star_in_Gyr + dt_in_Gyr) / All.NSNS_Rate_TAU, 1. - GFM_NSNS_DTD_POWERLAW_INDEX));
}


/*!
 * This routine computes yields from and number of NS-NS mergers
 */
void evolve_NSNS(MyDouble log_min_mass, MyDouble log_max_mass, MyFloat age_of_star_in_Gyr, MyFloat dt_in_Gyr, MyFloat metallicity, MyFloat initial_mass, stellar_evolution_data * sed)
{
  MyFloat NSNS_number_floor;
  MyFloat NSNS_random;

  sed->number_of_NSNS = get_NSNS_number_dtd_powerlaw(age_of_star_in_Gyr, dt_in_Gyr);
  sed->number_of_NSNS *= initial_mass * All.UnitMass_in_g / SOLAR_MASS / All.HubbleParam;
  // now, scale everything by the fact that we might have a fraction of a NSNS explosion going off
  NSNS_number_floor = floor(sed->number_of_NSNS);
  NSNS_random = get_random_number_aux();   /* drawing from this random number doesn't change the sequence of the main code's random numbers */
  if(NSNS_random <= (sed->number_of_NSNS - NSNS_number_floor))
    {
      NSNS_number_floor += 1;
    }
  sed->number_of_NSNS = NSNS_number_floor;

  if(All.NSNS_MassTransferOn == 1)
    {
      sed->metal_mass_released_chemtags[GFM_NSNS_CHEMTAG] += sed->number_of_NSNS * All.NSNS_MassPerEvent
          / All.UnitMass_in_g * SOLAR_MASS * All.HubbleParam / initial_mass;
    }
}
#endif




/*!
 * compute mass of stars dying at some age_of_star_in_Gyr based on inversion of lifetime tables
 */

MyFloat get_dying_mass_in_Msun(MyFloat age_of_star_in_Gyr, MyFloat metallicity)
{
  MyFloat mass = 0, mass1, mass2, d_metal, d_time1 = 0, d_time2 = 0, logage;

  int i_metal, i, i_time1 = -1, i_time2 = -1;

  if(age_of_star_in_Gyr <= 0)
    {
      mass = All.IMF_MaxMass_Msun;
      return mass;
    }

  logage = log10(age_of_star_in_Gyr * 1.E9);

  if(metallicity <= Lifetimes.Metallicity[0])
    {
      i_metal = 0;
      d_metal = 0.0;
    }
  else if(metallicity >= Lifetimes.Metallicity[Lifetimes.N_Z - 1])
    {
      i_metal = Lifetimes.N_Z - 2;
      d_metal = 1.0;
    }
  else
    {
      for(i_metal = 0; i_metal < Lifetimes.N_Z - 1; i_metal++)
        if(Lifetimes.Metallicity[i_metal + 1] > metallicity)
          break;

      d_metal = (metallicity - Lifetimes.Metallicity[i_metal]) / (Lifetimes.Metallicity[i_metal + 1] - Lifetimes.Metallicity[i_metal]);
    }

  if(logage >= Lifetimes.Dyingtime[i_metal][0])
    {
      i_time1 = 0;
      d_time1 = 0.0;
    }
  else if(logage <= Lifetimes.Dyingtime[i_metal][Lifetimes.N_MASS - 1])
    {
      i_time1 = Lifetimes.N_MASS - 2;
      d_time1 = 1.0;
    }

  if(logage >= Lifetimes.Dyingtime[i_metal + 1][0])
    {
      i_time2 = 0;
      d_time2 = 0.0;
    }
  else if(logage <= Lifetimes.Dyingtime[i_metal + 1][Lifetimes.N_MASS - 1])
    {
      i_time2 = Lifetimes.N_MASS - 2;
      d_time2 = 1.0;
    }

  i = Lifetimes.N_MASS;
  while(i >= 0 && (i_time1 == -1 || i_time2 == -1))
    {
      i--;
      if(Lifetimes.Dyingtime[i_metal][i] >= logage && i_time1 == -1)
        {
          i_time1 = i;
          d_time1 = (logage - Lifetimes.Dyingtime[i_metal][i_time1]) / (Lifetimes.Dyingtime[i_metal][i_time1 + 1] - Lifetimes.Dyingtime[i_metal][i_time1]);
        }
      if(Lifetimes.Dyingtime[i_metal + 1][i] >= logage && i_time2 == -1)
        {
          i_time2 = i;
          d_time2 = (logage - Lifetimes.Dyingtime[i_metal + 1][i_time2]) / (Lifetimes.Dyingtime[i_metal + 1][i_time2 + 1] - Lifetimes.Dyingtime[i_metal + 1][i_time2]);
        }
    }

  mass1 = interpol_1d(Lifetimes.Mass, i_time1, d_time1);
  mass2 = interpol_1d(Lifetimes.Mass, i_time2, d_time2);

  mass = (1.0 - d_metal) * mass1 + d_metal * mass2;

  if(mass > All.IMF_MaxMass_Msun)
    mass = All.IMF_MaxMass_Msun;

  return mass;
}


/*!
 * get lifetime from table
 */
MyFloat get_lifetime_in_Gyr(MyDouble mass, MyFloat metallicity)
{
  MyFloat time = 0, d_mass, d_metal;

  int i_mass, i_metal;


  if(mass <= Lifetimes.Mass[0])
    {
      i_mass = 0;
      d_mass = 0.0;
    }
  else if(mass >= Lifetimes.Mass[Lifetimes.N_MASS - 1])
    {
      i_mass = Lifetimes.N_MASS - 2;
      d_mass = 1.0;
    }
  else
    {
      for(i_mass = 0; i_mass < Lifetimes.N_MASS - 1; i_mass++)
        if(Lifetimes.Mass[i_mass + 1] > mass)
          break;

      d_mass = (mass - Lifetimes.Mass[i_mass]) / (Lifetimes.Mass[i_mass + 1] - Lifetimes.Mass[i_mass]);
    }

  if(metallicity <= Lifetimes.Metallicity[0])
    {
      i_metal = 0;
      d_metal = 0.0;
    }
  else if(metallicity >= Lifetimes.Metallicity[Lifetimes.N_Z - 1])
    {
      i_metal = Lifetimes.N_Z - 2;
      d_metal = 1.0;
    }
  else
    {
      for(i_metal = 0; i_metal < Lifetimes.N_Z - 1; i_metal++)
        if(Lifetimes.Metallicity[i_metal + 1] > metallicity)
          break;

      d_metal = (metallicity - Lifetimes.Metallicity[i_metal]) / (Lifetimes.Metallicity[i_metal + 1] - Lifetimes.Metallicity[i_metal]);
    }

  /* time in years */
  time = pow(10.0, interpol_2d(Lifetimes.Dyingtime, i_metal, i_mass, d_metal, d_mass));
  /* now in Gyr */
  time /= 1.0e9;

  return time;
}





void get_z_indicies(MyFloat log_metallicity, MyFloat * metal_values, int N_Z, int *iz_low, int *iz_high, MyFloat * dz)
{
  MyFloat deltaz, dz_local;
  int i1, i2;

  if(log_metallicity > GFM_MIN_METAL)
    {
      for(i1 = 0; i1 < N_Z - 1 && log_metallicity > metal_values[i1 + 1]; i1++);

      i2 = i1 + 1;

      if(i2 >= N_Z)
        i2 = N_Z - 1;

      if(log_metallicity >= metal_values[0] && log_metallicity <= metal_values[N_Z - 1])
        dz_local = log_metallicity - metal_values[i1];
      else
        dz_local = 0;

      deltaz = metal_values[i2] - metal_values[i1];

      if(deltaz > 0)
        dz_local = dz_local / deltaz;
      else
        dz_local = 0;
    }
  else
    {
      i1 = 0;
      i2 = 0;
      dz_local = 0.0;
    }

  *iz_low = i1;
  *iz_high = i2;
  *dz = dz_local;
}




/*!
 * integrates tabulated IMF with mode=
 * INTEGRATE_IMF_NUMBER integrate number
 * INTEGRATE_IMF_MASS integrate mass
 * INTEGRATE_IMF_YIELD integrate yields
 *
 * Note imf_integrand_mass is only required for INTEGRATE_IMF_YIELD and may be null for other modes.
 */

MyDouble integrate_imf(MyDouble log_min_mass, MyDouble log_max_mass, int mode, double *imf_integrand_mass)
{
  MyDouble result = 0.0;
  int ilow, ihigh, index;
  MyDouble dm;
  MyDouble integrand[GFM_N_MASS_BINS];

  get_imf_bins(log_min_mass, log_max_mass, &ilow, &ihigh);

  for(index = ilow; index < ihigh + 1; index++)
    {
      if(mode == INTEGRATE_IMF_NUMBER)
        integrand[index] = imf_by_number[index] * imf_mass_bin[index];  /* integrate number */
      else if(mode == INTEGRATE_IMF_MASS)
        integrand[index] = imf_by_number[index] * imf_mass_bin[index] * imf_mass_bin[index];    /* integrate mass */
      else if(mode == INTEGRATE_IMF_YIELD)
        integrand[index] = imf_integrand_mass[index] * imf_by_number[index] * imf_mass_bin[index];      /* integrate number * yield */
      else
        terminate("fail.");

      /* integrate using trapezoidal rule */
      result += integrand[index];
    }

  result = result - 0.5 * integrand[ilow] - 0.5 * integrand[ihigh];

  /* correct first bin */
  dm = (log_min_mass - imf_mass_bin_log10[ilow]) / imf_dlog10_Msun;

  if(dm < 0.5)
    result -= dm * integrand[ilow];
  else
    {
      result -= 0.5 * integrand[ilow];
      result -= (dm - 0.5) * integrand[ilow + 1];
    }

  /* correct last bin */
  dm = (log_max_mass - imf_mass_bin_log10[ihigh - 1]) / imf_dlog10_Msun;

  if(dm < 0.5)
    {
      result -= 0.5 * integrand[ihigh];
      result -= (0.5 - dm) * integrand[ihigh - 1];
    }
  else
    result -= (1 - dm) * integrand[ihigh];

  result *= imf_dlog10_Msun * log(10.0);        /* log(10) since mass function tabulated as function of log_10(mass) */

  return result;
}



void get_imf_bins(MyDouble log_min_mass, MyDouble log_max_mass, int *ilow, int *ihigh)
{
  int i1, i2;

  if(log_min_mass < imf_mass_bin_log10[0])
    log_min_mass = imf_mass_bin_log10[0];

  if(log_min_mass > imf_mass_bin_log10[GFM_N_MASS_BINS - 1])
    log_min_mass = imf_mass_bin_log10[GFM_N_MASS_BINS - 1];

  if(log_max_mass < imf_mass_bin_log10[0])
    log_max_mass = imf_mass_bin_log10[0];

  if(log_max_mass > imf_mass_bin_log10[GFM_N_MASS_BINS - 1])
    log_max_mass = imf_mass_bin_log10[GFM_N_MASS_BINS - 1];

  for(i1 = 0; i1 < GFM_N_MASS_BINS - 2 && imf_mass_bin_log10[i1 + 1] < log_min_mass; i1++);
  for(i2 = 1; i2 < GFM_N_MASS_BINS - 1 && imf_mass_bin_log10[i2] < log_max_mass; i2++);

  *ilow = i1;
  *ihigh = i2;
}

/*!
 * define IMF arrays based on given property, and normalize such that total mass = 1
 */
double define_imf(double imf_param)
{
  int i;
  double mass;
#if defined(GFM_VARIABLE_IMF) && (GFM_VARIABLE_IMF == 0)
  double slope;
#endif

  for(i = 0; i < GFM_N_MASS_BINS; i++)
    {
      mass = imf_mass_bin[i];

#if defined(GFM_CONST_IMF) && (GFM_CONST_IMF == 0)
      /* Chabrier 2003 */
      if(mass > 1.0)
        imf_by_number[i] = 0.237912 * pow(mass, -2.3);
      else
        imf_by_number[i] = 0.852464 * exp(pow((log10(mass) - log10(0.079)), 2.0) / (-2.0 * pow(0.69, 2))) / mass;
#elif defined(GFM_CONST_IMF) && (GFM_CONST_IMF == 1)
      /* pure power-law IMF */
      imf_by_number[i] = pow(mass, All.IMFslope);

#elif defined(GFM_VARIABLE_IMF) && (GFM_VARIABLE_IMF == 0)
      /* pure power-law IMF that depends on the velocity dispersion */
      slope = -(2.3 * log10(imf_param / 200) + 2.13);   /* Spiniello et al. 2014, MNRAS, 438, 1483 */
      imf_by_number[i] = pow(mass, slope);

#else
#error "IMF mode is not defined"
#endif
    }

  double norm = integrate_imf(log10(All.IMF_MinMass_Msun), log10(All.IMF_MaxMass_Msun), INTEGRATE_IMF_MASS, NULL);

  for(i = 0; i < GFM_N_MASS_BINS; i++)
    imf_by_number[i] /= norm;

  return norm;
}

/*!
 * initialize global IMF
 */
void init_imf(void)
{
  int i;
  imf_dlog10_Msun = (log10(All.IMF_MaxMass_Msun) - log10(All.IMF_MinMass_Msun)) / (double) (GFM_N_MASS_BINS - 1);

  for(i = 0; i < GFM_N_MASS_BINS; i++)
    {
      imf_mass_bin_log10[i] = log10(All.IMF_MinMass_Msun) + i * imf_dlog10_Msun;
      imf_mass_bin[i] = pow(10, imf_mass_bin_log10[i]);
    }

#ifdef GFM_CONST_IMF
  double norm = define_imf(-1);

  mpi_printf("GFM_STELLAR_EVOLUTION: IMF normalization -> norm_factor=%g   N_tot=%g   M_tot=%g M_sun   N_SNII=%g   M_SNII=%g M_sun\n",
             norm,
             integrate_imf(log10(All.IMF_MinMass_Msun), log10(All.IMF_MaxMass_Msun), INTEGRATE_IMF_NUMBER, NULL),
             integrate_imf(log10(All.IMF_MinMass_Msun), log10(All.IMF_MaxMass_Msun), INTEGRATE_IMF_MASS, NULL),
             integrate_imf(log10(All.SNII_MinMass_Msun), log10(All.SNII_MaxMass_Msun), INTEGRATE_IMF_NUMBER, NULL),
             integrate_imf(log10(All.SNII_MinMass_Msun), log10(All.SNII_MaxMass_Msun), INTEGRATE_IMF_MASS, NULL));

#ifdef USE_SFR
  All.FactorSN = calc_FactorSN();
  mpi_printf("GFM_STELLAR_EVOLUTION: Rescaling FactorSN(beta) to %g\n", All.FactorSN);
#ifdef GFM_WINDS
  All.WindEgySpecSN = calc_WindEgySpecSN();
  mpi_printf("GFM_STELLAR_EVOLUTION: Rescaling WindEgySpecSN to %g (internal units) = %g (erg/Msun)\n", All.WindEgySpecSN,
             All.WindEgySpecSN / (All.UnitMass_in_g / SOLAR_MASS / All.UnitEnergy_in_cgs));
#endif
#endif
#endif
}

#ifdef USE_SFR
double calc_FactorSN(void)
{
  return integrate_imf(log10(All.SNII_MinMass_Msun), log10(All.SNII_MaxMass_Msun), INTEGRATE_IMF_MASS, NULL);
}
#endif

#ifdef GFM_WINDS
double calc_WindEgySpecSN(void)
{
  return GFM_SNII_ENERGY * integrate_imf(log10(All.SNII_MinMass_Msun), log10(All.SNII_MaxMass_Msun), INTEGRATE_IMF_NUMBER,
                                         NULL) / calc_FactorSN() * All.UnitMass_in_g / SOLAR_MASS / All.UnitEnergy_in_cgs;
}
#endif

/*!
 * initialize SNIa rates, and normalize such that total number of SNIa is equal to All.SNIa_Rate_Norm between SNIa_Rate_TAU
 * and the age of the universe
 */
void init_SNIa_rates(void)
{
  double age_of_universe_in_Gyr = 0.0;

  if(All.ComovingIntegrationOn)
    age_of_universe_in_Gyr += get_time_difference_in_Gyr(All.TimeBegin, 1.0);
  else
    age_of_universe_in_Gyr = 13.7;      /* taken from Maoz et al. (2012) */

#ifdef GFM_SNIA_DTD_POWERLAW
  SNIa_Rate_Norm = All.SNIa_Rate_Norm / (1.0 - pow(age_of_universe_in_Gyr / All.SNIa_Rate_TAU, 1. - GFM_SNIA_DTD_POWERLAW_INDEX));
#elif GFM_SNIA_DTD_EFOLDING
  SNIa_Rate_Norm = All.SNIa_Rate_Norm / (1.0 - exp(-age_of_universe_in_Gyr / All.SNIa_Rate_TAU));
#else
  terminate("GFM_STELLAR_EVOLUTION: Unknown delay time distribution");
#endif

  mpi_printf("GFM_STELLAR_EVOLUTION: age of the universe %g\n", age_of_universe_in_Gyr);
  mpi_printf("GFM_STELLAR_EVOLUTION: Rescaling SNIa_Rate_Norm to %g\n", SNIa_Rate_Norm);
}



/*!
 * routine to test the stellar evolution model:
 * -set GFM_STELLAR_EVOLUTION=2 in Config.sh to run only this routine
 */
void test_stellar_evolution(void)
{
  int iPart = 0, k_elem, fc;
  MyFloat dtime = 5e-6, dtime_in_Gyr;
  MyFloat dtime_output = 1e-4, next_outputtime = 0.0;
  MyFloat age_of_star, age_of_star_in_Gyr;
  MyFloat CumMetalsReleased[GFM_N_CHEM_ELEMENTS];
  MyFloat CumMetalMassReleased;
#ifdef GFM_DUST
  int chan;
  MyFloat DustMassReleased[GFM_N_CHEM_ELEMENTS];
  FILE *dust_file;
#endif
  FILE *fd[6];
  char fname[255];
  stellar_evolution_data sed;

  mpi_printf("GFM_STELLAR_EVOLUTION: STELLAR EVOLUTION TEST (dtime=%g, dtime_output=%g)\n", dtime, dtime_output);

    /* perform the required init steps */
  read_parameter_file(ParameterFile);

  /* make sure that we are running physical */
  All.ComovingIntegrationOn = 0;


  mymalloc_init();

  set_units();

  /* init stellar evolution model -> select which elements should be put at which index */
  strcpy(ElementNames[0], "Hydrogen");
  strcpy(ElementNames[1], "Helium");
  strcpy(ElementNames[2], "Carbon");
  strcpy(ElementNames[3], "Nitrogen");
  strcpy(ElementNames[4], "Oxygen");
  strcpy(ElementNames[5], "Neon");
  strcpy(ElementNames[6], "Magnesium");
  strcpy(ElementNames[7], "Silicon");
  strcpy(ElementNames[8], "Iron");
#ifdef GFM_NORMALIZED_METAL_ADVECTION
  strcpy(ElementNames[9], "OtherMetals");
#endif // end advection


  init_imf();
  mpi_printf("GFM_STELLAR_EVOLUTION: IMF initialized.\n");
  init_yields();
  mpi_printf("GFM_STELLAR_EVOLUTION: Yields initialized.\n");
  init_SNIa_rates();
  mpi_printf("GFM_STELLAR_EVOLUTION: Type Ia rates initialized.\n");

  /*
     -allocate particle data and star particle data (100 particles only)
     -needed because stellar evolution code works on it
   */
  StarP = (struct star_particle_data *) mymalloc("StarP", 100 * sizeof(struct star_particle_data));
  P = (struct particle_data *) mymalloc("P", 100 * sizeof(struct particle_data));

  for(fc = 0; fc < 6; fc++)
    {
      sprintf(fname, "stellar_evolution_%d.dat", fc);
      printf("GFM_STELLAR_EVOLUTION: opening %s...\n", fname);
      fd[fc] = fopen(fname, "w");
    }
#ifdef GFM_DUST
  printf("GFM_DUST: opening stellar_evolution_dust.dat...\n");
  dust_file = fopen("stellar_evolution_dust.dat", "w");
#endif

  P[iPart].Mass = 1.0;
  P[iPart].AuxDataID = 0;
  P[iPart].Type = 4;

  P[iPart].SofteningType = All.SofteningTypeOfPartType[4];

#ifdef INDIVIDUAL_GRAVITY_SOFTENING
  if(((1 << P[iPart].Type) & (INDIVIDUAL_GRAVITY_SOFTENING)))
    P[iPart].SofteningType = get_softening_type_from_mass(P[iPart].Mass);
#endif

  STP(iPart).InitialMass = 1.0;
  STP(iPart).BirthTime = 1.0;   /* mark as stellar particle, and not wind. not important, but do_stellar_evolution() would stop otherwise */

  /*
     -initial metallicity of star, used for lookup metal dependend yields
     -also used for total amount of ejected metals (=metallicity * ejected mass)
     -metallicity = metal mass fraction
     -solar = 0.0127
   */
  STP(iPart).Metallicity = GFM_SOLAR_METALLICITY;       // 0.0127 for diff test

  /*
     -set the total amount of initial metals (and H, He) in internal mass units
     -these metals pass through the star, they do not enter the yield table readout
     i.e. changing one element here does not effect any other element
   */
  /* primordial: Y=0.25, Z=0 (BBN)  linear to  solar: Y=0.2806, Z=0.0127  --> Y=0.25 + 2.41*Z */
  MyDouble Y = 1.0 - GFM_INITIAL_ABUNDANCE_HYDROGEN;
  MyDouble X = GFM_INITIAL_ABUNDANCE_HYDROGEN;
  //MyDouble Y = 0.25 + 2.41 * STP(iPart).Metallicity;
  //MyDouble X = 1.0 - Y - STP(iPart).Metallicity;
  STP(iPart).MassMetals[0] = P[iPart].Mass * X;
  STP(iPart).MassMetals[1] = P[iPart].Mass * Y;
  STP(iPart).MassMetals[2] = STP(iPart).Metallicity / GFM_SOLAR_METALLICITY * P[iPart].Mass * 2.06e-3;
  STP(iPart).MassMetals[3] = STP(iPart).Metallicity / GFM_SOLAR_METALLICITY * P[iPart].Mass * 8.36e-4;
  STP(iPart).MassMetals[4] = STP(iPart).Metallicity / GFM_SOLAR_METALLICITY * P[iPart].Mass * 5.49e-3;
  STP(iPart).MassMetals[5] = STP(iPart).Metallicity / GFM_SOLAR_METALLICITY * P[iPart].Mass * 1.41e-3;
  STP(iPart).MassMetals[6] = STP(iPart).Metallicity / GFM_SOLAR_METALLICITY * P[iPart].Mass * 5.91e-4;
  STP(iPart).MassMetals[7] = STP(iPart).Metallicity / GFM_SOLAR_METALLICITY * P[iPart].Mass * 6.83e-4;
  STP(iPart).MassMetals[8] = STP(iPart).Metallicity / GFM_SOLAR_METALLICITY * P[iPart].Mass * 1.10e-3;

#ifdef GFM_CHEMTAGS
  for(k_elem = 0; k_elem < GFM_N_CHEM_TAGS; k_elem++)
    STP(iPart).MassMetalsChemTags[k_elem] = 0.0;
#endif

  /* array for cumulative ejecta */
  for(k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
    CumMetalsReleased[k_elem] = 0.0;

  CumMetalMassReleased = 0.0;

#ifdef GFM_DUST
  for(k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
    {
      DustMassReleased[k_elem] = 0.0;
    }
#endif

  /* integrate over time */
  double NtotSNIa = 0.0;
  double NtotSNII = 0.0;
  for(age_of_star = 0.0; age_of_star <= 10.0; age_of_star += dtime)
    {
      age_of_star_in_Gyr = get_time_difference_in_Gyr(0, age_of_star);
      dtime_in_Gyr = get_time_difference_in_Gyr(0, dtime);

      /* main stellar evolution routine */
      do_stellar_evolution(age_of_star_in_Gyr, dtime_in_Gyr, iPart, &sed);

      for(k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
        CumMetalsReleased[k_elem] += sed.metal_mass_released[k_elem];

      CumMetalMassReleased += sed.total_metal_mass_released;

      NtotSNIa += sed.number_of_SNIa;
      NtotSNII += sed.number_of_SNII;

#ifdef GFM_DUST
      for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
        {
          for(k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
            {
              DustMassReleased[k_elem] += sed.dust_mass_released[chan][k_elem];
            }
        }
#endif

      /* do an output? */
      if(age_of_star >= next_outputtime)
        {
          fprintf(fd[0], "%g %g %g\n", age_of_star_in_Gyr, sed.number_of_SNIa / All.HubbleParam / (dtime_in_Gyr * 1.e9), NtotSNIa);

          /* [ejected metals over ejected iron] */
          if(CumMetalsReleased[8] > 0)
            fprintf(fd[1], "%g %g %g %g %g %g %g %g %g %g %g\n", age_of_star_in_Gyr, age_of_star,
                    (CumMetalsReleased[0] / CumMetalsReleased[8]) / (0.7065 / 1.1e-3),
                    (CumMetalsReleased[1] / CumMetalsReleased[8]) / (0.2806 / 1.1e-3),
                    (CumMetalsReleased[2] / CumMetalsReleased[8]) / (2.06e-3 / 1.1e-3),
                    (CumMetalsReleased[3] / CumMetalsReleased[8]) / (8.36e-4 / 1.1e-3),
                    (CumMetalsReleased[4] / CumMetalsReleased[8]) / (5.49e-3 / 1.1e-3),
                    (CumMetalsReleased[5] / CumMetalsReleased[8]) / (1.41e-3 / 1.1e-3),
                    (CumMetalsReleased[6] / CumMetalsReleased[8]) / (5.91e-4 / 1.1e-3),
                    (CumMetalsReleased[7] / CumMetalsReleased[8]) / (6.83e-4 / 1.1e-3), (CumMetalsReleased[8] / CumMetalsReleased[8]) / (1.1e-3 / 1.1e-3));
          else
            fprintf(fd[1], "%g %g 0 0 0 0 0 0 0 0 0\n", age_of_star_in_Gyr, age_of_star);

          /* [ejected metals over ejected hydrogen] */
          if(CumMetalsReleased[0] > 0)
            fprintf(fd[2], "%g %g %g %g %g %g %g %g %g %g %g\n", age_of_star_in_Gyr, age_of_star,
                    (CumMetalsReleased[0] / CumMetalsReleased[0]) / (0.7065 / 0.7065),
                    (CumMetalsReleased[1] / CumMetalsReleased[0]) / (0.2806 / 0.7065),
                    (CumMetalsReleased[2] / CumMetalsReleased[0]) / (2.06e-3 / 0.7065),
                    (CumMetalsReleased[3] / CumMetalsReleased[0]) / (8.36e-4 / 0.7065),
                    (CumMetalsReleased[4] / CumMetalsReleased[0]) / (5.49e-3 / 0.7065),
                    (CumMetalsReleased[5] / CumMetalsReleased[0]) / (1.41e-3 / 0.7065),
                    (CumMetalsReleased[6] / CumMetalsReleased[0]) / (5.91e-4 / 0.7065),
                    (CumMetalsReleased[7] / CumMetalsReleased[0]) / (6.83e-4 / 0.7065), (CumMetalsReleased[8] / CumMetalsReleased[0]) / (1.1e-3 / 0.7065));
          else
            fprintf(fd[2], "%g %g 0 0 0 0 0 0 0 0 0\n", age_of_star_in_Gyr, age_of_star);


          /* remaining stellar mass fraction and cumulative released mass */
          fprintf(fd[3], "%g %g %g %g %g %g %g %g %g %g %g %g %g\n", age_of_star_in_Gyr, age_of_star,
                  P[iPart].Mass / STP(iPart).InitialMass, CumMetalMassReleased, CumMetalsReleased[0], CumMetalsReleased[1],
                  CumMetalsReleased[2], CumMetalsReleased[3], CumMetalsReleased[4], CumMetalsReleased[5], CumMetalsReleased[6], CumMetalsReleased[7], CumMetalsReleased[8]);

          /* [ejected metals] */
          fprintf(fd[4], "%g %g %g %g %g %g %g %g %g %g %g\n", age_of_star_in_Gyr, age_of_star,
                  (CumMetalsReleased[0]),
                  (CumMetalsReleased[1]),
                  (CumMetalsReleased[2]), (CumMetalsReleased[3]), (CumMetalsReleased[4]), (CumMetalsReleased[5]), (CumMetalsReleased[6]), (CumMetalsReleased[7]), (CumMetalsReleased[8]));

          /* [SN explosions] */
          fprintf(fd[5], "%g %g %g\n", age_of_star_in_Gyr, NtotSNIa, NtotSNII);

#ifdef GFM_DUST
          fprintf(dust_file, "%g %g %g %g %g %g %g %g %g %g %g\n",
                  age_of_star_in_Gyr, age_of_star,
                  (DustMassReleased[0]),
                  (DustMassReleased[1]),
                  (DustMassReleased[2]), (DustMassReleased[3]), (DustMassReleased[4]), (DustMassReleased[5]), (DustMassReleased[6]), (DustMassReleased[7]), (DustMassReleased[8]));
#endif

          next_outputtime += dtime_output;
        }
    }

  printf("NSNIatot=%g, expected=%g\n", NtotSNIa, All.SNIa_Rate_Norm);
  printf("NSNIItot=%g, SNII/SNIa=%g\n", NtotSNII, NtotSNII / NtotSNIa);
  printf("ejected mass fraction: %g\n", (STP(iPart).InitialMass - P[iPart].Mass) / STP(iPart).InitialMass);
  printf("stellar age in Gyr: %g\n", age_of_star_in_Gyr);

  for(fc = 0; fc < 6; fc++)
    fclose(fd[fc]);

#ifdef GFM_DUST
  fclose(dust_file);
#endif

  myfree(P);
  myfree(StarP);
}

#endif
