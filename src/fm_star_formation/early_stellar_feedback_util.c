/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/fm_star_formation/early_stellar_feedback_util.c
 * \date        05/2014
 * \author      Laura Sales 
 * \brief        
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 *   routines for early stellar feedback
 */

#include "../allvars.h"
#include "../proto.h"


#ifdef FM_EARLY_STAR_FEEDBACK

int is_doing_early_feedback(int i)
{
  if(StarParticle[i].EarlyTotalEnergyReleased > 0)
    return 1;

  return 0;
}

MyDouble compute_EarlyFeedback_energy(MyDouble age_star_in_gyr, MyFloat mass, MyFloat dt_in_sec)
{
  MyDouble lm_ssp, log_age, stellum;

/* below is taken from Paul's Hopkins-like ism/ism_HII_heating.c */

  if(age_star_in_gyr < 0.00259710)
    {
      lm_ssp = 328.04726;
    }
  else
    {
      log_age = log10(age_star_in_gyr) - (-2.6226865);
      lm_ssp = 468.1832 * pow(10., -2.2408834 * log_age - 2.9534785 * log_age * log_age);
    }

  stellum = lm_ssp * mass;
  stellum *= 1.95 * All.UnitMass_in_g / All.HubbleParam;        /* CGS units */
  /*LVS: 1.95 to go from Msun to Lsun */
  //

/* ====================================================================================================*/

  MyDouble Injected_Energy = All.EarlyFeedbackEfficiency * dt_in_sec * stellum; /* a fraction EarlyFeedbackEfficiency is given to the gas */
  Injected_Energy *= (All.HubbleParam / All.UnitEnergy_in_cgs); /* to code units */

  return Injected_Energy;
}

#endif
