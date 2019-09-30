/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/runge_kutta_full.h
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

#include "arepoconfig.h"
#include "allvars.h"
#include "proto.h"
#include "runge_kutta_full.h"

#ifdef RUNGE_KUTTA_FULL_UPDATE

void rk_save_conservative_variables()
{
  int idx, i, k;
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Type != 0 || (P[i].ID == 0 && P[i].Mass == 0))
        continue;

      SphP[i].rk.Mass = P[i].Mass;
      for(k = 0; k < 3; k++) SphP[i].rk.Momentum[k] = SphP[i].Momentum[k];
      SphP[i].rk.Energy = SphP[i].Energy;
#ifdef MAXSCALARS
      for(k = 0; k < MAXSCALARS; k++) SphP[i].rk.Scalars[k] = *(MyFloat *) (((char *) (&SphP[i])) + scalar_elements[k].offset_mass);
#endif
#ifdef MHD
      for(k = 0; k < 3; k++) SphP[i].rk.BConserved[k] = SphP[i].BConserved[k];
#endif
    }
}

void rk_finish_step()
{
#ifndef GENERAL_RELATIVITY
  update_primitive_variables(); //because it's done before calling this routine
#endif
  exchange_primitive_variables();
  calculate_gradients();
  exchange_primitive_variables_and_gradients();
  compute_interface_fluxes(&Mesh);
  
  int idx, i, k;
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Type != 0 || (P[i].ID == 0 && P[i].Mass == 0))
        continue;
        
      P[i].Mass = 0.5 * (SphP[i].rk.Mass + P[i].Mass);
      for(k = 0; k < 3; k++) SphP[i].Momentum[k] = 0.5 * (SphP[i].rk.Momentum[k] + SphP[i].Momentum[k]);
      SphP[i].Energy = 0.5 * (SphP[i].rk.Energy + SphP[i].Energy);
#ifdef MAXSCALARS
      for(k = 0; k < MAXSCALARS; k++) *(MyFloat *) (((char *) (&SphP[i])) + scalar_elements[k].offset_mass) = 0.5 * (SphP[i].rk.Scalars[k] + *(MyFloat *) (((char *) (&SphP[i])) + scalar_elements[k].offset_mass));
#endif
#ifdef MHD
      for(k = 0; k < 3; k++) SphP[i].BConserved[k] = 0.5 * (SphP[i].rk.BConserved[k] + SphP[i].BConserved[k]);
#endif
    }
}

void rk_derefinement_add( struct conservative_variables *target, struct conservative_variables *source, double fac )
{
  int k;

  target->Mass += fac * source->Mass;
  for(k = 0; k < 3; k++) target->Momentum[k] += fac * source->Momentum[k];
  target->Energy += fac * source->Energy;
#ifdef MAXSCALARS
  for(k = 0; k < MAXSCALARS; k++) target->Scalars[k] += fac * source->Scalars[k];
#endif
#ifdef MHD
  for(k = 0; k < 3; k++) target->BConserved[k] += fac * source->BConserved[k];
#endif
}

void rk_derefinement_set( struct conservative_variables *target, struct conservative_variables *source, double fac )
{
  int k;

  target->Mass = fac * source->Mass;
  for(k = 0; k < 3; k++) target->Momentum[k] = fac * source->Momentum[k];
  target->Energy = fac * source->Energy;
#ifdef MAXSCALARS
  for(k = 0; k < MAXSCALARS; k++) target->Scalars[k] = fac * source->Scalars[k];
#endif
#ifdef MHD
  for(k = 0; k < 3; k++) target->BConserved[k] = fac * source->BConserved[k];
#endif
}

void rk_multiply( struct conservative_variables *target, double fac )
{
  int k;

  target->Mass *= fac;
  for(k = 0; k < 3; k++) target->Momentum[k] *= fac;
  target->Energy *= fac;
#ifdef MAXSCALARS
  for(k = 0; k < MAXSCALARS; k++) target->Scalars[k] *= fac;
#endif
#ifdef MHD
  for(k = 0; k < 3; k++) target->BConserved[k] *= fac;
#endif
}

#endif
