/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/tgchem/tgchem_photo.c
 * \date        01/2013
 * \author      Thomas Greif
 * \brief       Primordial chemistry and cooling network
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 */

#ifdef TGCHEM_TEST
#include "tgchem_test.h"
#else
#include "../allvars.h"
#include "../proto.h"
#endif


void tgchem_photo()
{
  double lsob, N_H2, x, b5, f_shield;

  // H2 line cooling using radiative transfer
  TGCD.HRHeatRate = 0.;

#ifdef HEALRAY
  if(HRD.SourceFlag >= 0)
    {
      TGCHEM_TEMP();
      TGCHEM_COOLRATE(1);

      TGCD.HRHeatRate = -TGCD.AbH2 * TGCD.NH * TGCD.CoolRate[1] + SphP[TGCD.Index].HeatRate;
    }
#endif

  // H2 dissociation by LW radiation (Wolcott-Green et al. 2011)
  if(!TGCD.ChemMode)
    {
      if(TGCD.DivVel)
        lsob = dmin(TGCD.VThermal / dabs(TGCD.DivVel), TGCD.JeansLength);
      else
        lsob = TGCD.JeansLength;

      N_H2 = TGCD.AbH2 * TGCD.NH * lsob;

      x = N_H2 / 5e14;

      b5 = TGCD.VThermal / 1e5;

      f_shield = 0.965 / pow(1. + x / b5, 1.1) + 0.035 / sqrt(1. + x) * exp(-8.5e-4 * sqrt(1. + x));

      TGCD.PhRate[0] = TGCD.LWDissRate * TGCD.ChemJ21 * f_shield;
    }
  else
    TGCD.PhRate[0] = 0;

#ifdef HEALRAY
  TGCD.HIPhotonEnergy = HRD.HIPhotonEnergy;
#else
  TGCD.HIPhotonEnergy = 0;
#endif
}
