/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/sidm/sidm.c
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
#include "sidm_vars.h"

/* generic driver routine for the scatter process */
void sidm_evaluate_scatter_process(scatter_process_data_in * sprdata_in, scatter_process_data_out * sprdata_out, MyIDType pID)
{

#ifndef SIDM_NO_NGB_SEL
  /* check that the reaction matches the particles states; i.e. make sure we have not selected the wrong scatter partner */
  if((sprdata_in->State1 != SMSIDM[sprdata_in->Reaction].In1) || (sprdata_in->State2 != SMSIDM[sprdata_in->Reaction].In2))
    terminate("SIDM: state scatter reaction mismatch: reaction=%d state1=%d state2=%d\n", sprdata_in->State1, sprdata_in->State2, sprdata_in->Reaction);
#endif

  sidm_scatter_in_to_out(sprdata_in, sprdata_out, pID);
}


void sidm_scatter_in_to_out(scatter_process_data_in * sprdata_in, scatter_process_data_out * sprdata_out, MyIDType pID)
{
  double epsilon = 1.0;
  double vcm_in[3], dv_in[3], xunit[3], dvabs_in;
  double delta_E_1, delta_E_2, delta_E;
  unsigned char reaction = sprdata_in->Reaction;
  double in_mass1 = sprdata_in->Mass1;
  double in_mass2 = sprdata_in->Mass2;
  double out_mass1, out_mass2;

  sidm_get_delta_mass(reaction, in_mass1, in_mass2, &out_mass1, &out_mass2);

  sprdata_out->OutState1 = SMSIDM[reaction].Out1;
  sprdata_out->OutState2 = SMSIDM[reaction].Out2;

  sprdata_out->Mass1 = out_mass1;
  sprdata_out->Mass2 = out_mass2;

  sidm_get_delta_energy(reaction, &delta_E_1, &delta_E_2);

  delta_E = delta_E_1 + delta_E_2;

  sidm_get_velocities(sprdata_in, sprdata_out, vcm_in, dv_in, &dvabs_in);

  if(delta_E != 0.0)
    {
      double mu_in = (in_mass1 * in_mass2) / (in_mass1 + in_mass2);
      double mu_out = (out_mass1 * out_mass2) / (out_mass1 + out_mass2);
      double fac = mu_in/mu_out * (1 + (2.0 * delta_E) / (mu_in * dvabs_in*dvabs_in));   
      if (fac < 0.0)  //maximally inelastic case, all relative kin. energy absorbed  FIXME 
        fac = 0.0;
      epsilon = sqrt(fac);
    }

#ifdef SIDM_NO_ENERGYCHANGE
  /* force scattering to be elastic; this is useful to compare elastic vs. inelastic scattering with otherwise the same cross sections */
  epsilon = 1.0;
#endif

  sidm_get_random_unit_sphere_vector(get_random_number(), get_random_number(), xunit);

  sidm_set_velocities(sprdata_in, sprdata_out, dvabs_in, epsilon, vcm_in, xunit);

}


/* set velocities after scattering process */
void sidm_set_velocities(scatter_process_data_in * sprdata_in, scatter_process_data_out * sprdata_out, double dvabs_in, double epsilon, double *vcm_in, double *xunit)
{
  int k;

  double fac_cm = (sprdata_in->Mass1 + sprdata_in->Mass2) / (sprdata_out->Mass1 + sprdata_out->Mass2);
  double fac_1  = sprdata_out->Mass2 / (sprdata_out->Mass1 + sprdata_out->Mass2);
  double fac_2  = sprdata_out->Mass1 / (sprdata_out->Mass1 + sprdata_out->Mass2);

  for(k = 0; k < 3; k++)
    {
      /* physical velocity to internal velocity for both particles */
      sprdata_out->Vel1[k] = 1. / SIDM_avel * (fac_cm * vcm_in[k] - fac_1 * dvabs_in * epsilon * xunit[k]);
      sprdata_out->Vel2[k] = 1. / SIDM_avel * (fac_cm * vcm_in[k] + fac_2 * dvabs_in * epsilon * xunit[k]);
        
      if ((isnan(sprdata_out->Vel1[k])) || (isinf(sprdata_out->Vel1[k])))
        terminate("SIDM: BAD SCATTER VELOCITIES: k=%d vel=%g fac_1=%g epsilon=%g fac_cm=%g vcm_in=%g\n", k, sprdata_out->Vel1[k], fac_1, epsilon, fac_cm, vcm_in[k]);
      if ((isnan(sprdata_out->Vel2[k])) || (isinf(sprdata_out->Vel2[k])))
        terminate("SIDM: BAD SCATTER VELOCITIES: k=%d vel=%g fac_2=%g epsilon=%g fac_cm=%g vcm_in=%g\n", k, sprdata_out->Vel2[k], fac_2, epsilon, fac_cm, vcm_in[k]);
    }
}


/* compute center of mass velocities */
void sidm_get_velocities(scatter_process_data_in * sprdata_in, scatter_process_data_out * sprdata_out, double *vcm_in, double *dv_in, double *dvabs_in)
{
  int k;
  for(k = 0; k < 3; k++)
    {
      /* physical center-of-mass velocity */
      vcm_in[k] = SIDM_avel * (sprdata_in->Mass1 * sprdata_in->Vel1[k] + sprdata_in->Mass2 * sprdata_in->Vel2[k]) / (sprdata_in->Mass1 + sprdata_in->Mass2);

      /* physical relative velocity */
      dv_in[k] = SIDM_avel * (sprdata_in->Vel1[k] - sprdata_in->Vel2[k]);
    }

  *dvabs_in = sqrt(dv_in[0] * dv_in[0] + dv_in[1] * dv_in[1] + dv_in[2] * dv_in[2]);
}

 /* unit sphere vector for isotropic scattering */
void sidm_get_random_unit_sphere_vector(double randnum1, double randnum2, double *xunit)
{
  double phi = 2 * M_PI * randnum1;
  double theta = acos(randnum2 * 2 - 1);

  xunit[0] = sin(theta) * cos(phi);
  xunit[1] = sin(theta) * sin(phi);
  xunit[2] = cos(theta);
}
