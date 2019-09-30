/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/tgchem/tgchem_chem.c
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


void tgchem_chem(N_Vector dspecies)
{
  int i, j, k;
  double rate, cr, ds, k1, k2;

  for(i = 0; i < TGCHEM_NUM_ABUNDANCES; i++)
    {
      for(j = 0; j < TGCHEM_NUM_RATES; j++)
        TGCD.CrRate[i][j] = TGCD.DsRate[i][j] = TGCD.ChemicalRate[i][j] = 0;

      TGCD.TotCrRate[i] = TGCD.TotDsRate[i] = TGCD.TotChemicalRate[i] = 0;
    }

  TGCD.ChemicalRateSum = 0;

  // H-

  // Creation

  i = 0;
  j = 0;
  k = 0;

  // H + e -> H- + ph
  TGCD.CrRate[i][j++] = rate = TGCD.ChemRate[0] * TGCD.AbHI * TGCD.AbE * TGCD.NH;
  //TGCD.ChemicalRate[i][k++] = -TGCHEM_CHI_HM * TGCD.NH * rate;

  // Destruction

  j = 0;

  // H + H- -> H2 + e
  TGCD.DsRate[i][j++] = rate = -TGCD.ChemRate[1] * TGCD.AbHI * TGCD.NH;

  for(j = cr = ds = 0; j < TGCHEM_NUM_RATES; j++)
    {
      cr += TGCD.CrRate[i][j];
      ds += dabs(TGCD.DsRate[i][j]);
    }

  if(ds == 0)
    TGCD.AbHM = 0;
  else
    TGCD.AbHM = dmin(cr / ds, TGCD.AbMax[0]);

  for(j = 0; j < TGCHEM_NUM_RATES; j++)
    TGCD.DsRate[i][j] *= TGCD.AbHM;

  // H2

  // Creation

  i++;
  j = 0;
  k = 0;

  // H + H- -> H2 + e
  TGCD.CrRate[i][j++] = rate = TGCD.ChemRate[1] * TGCD.AbHI * TGCD.AbHM * TGCD.NH;
  TGCD.ChemicalRate[i][k++] = TGCHEM_CHI_H2 * TGCD.NH * rate;

  // 3H -> H + H2
  TGCD.CrRate[i][j++] = rate = TGCD.ChemRate[2] * TGCD.AbHI * TGCD.AbHI * TGCD.AbHI * TGCD.NH * TGCD.NH;
  TGCD.ChemicalRate[i][k++] = TGCHEM_CHI_H2 * TGCD.NH * rate;

  // 2H + H2 -> 2H2
  TGCD.CrRate[i][j++] = rate = TGCD.ChemRate[2] / 8. * TGCD.AbHI * TGCD.AbHI * TGCD.AbH2 * TGCD.NH * TGCD.NH;
  TGCD.ChemicalRate[i][k++] = TGCHEM_CHI_H2 * TGCD.NH * rate;

  // Destruction

  j = 0;

  // H + H2 -> 3H
  TGCD.DsRate[i][j++] = rate = -TGCD.ChemRate[3] * TGCD.AbHI * TGCD.AbH2 * TGCD.NH;
  TGCD.ChemicalRate[i][k++] = TGCHEM_CHI_H2 * TGCD.NH * rate;

  // H2 + H2 -> 2H + H2
  TGCD.DsRate[i][j++] = rate = -TGCD.ChemRate[3] / 8. * TGCD.AbH2 * TGCD.AbH2 * TGCD.NH;
  TGCD.ChemicalRate[i][k++] = TGCHEM_CHI_H2 * TGCD.NH * rate;

  // H2 + ph -> 2H
  TGCD.DsRate[i][j++] = rate = -TGCD.PhRate[0] * TGCD.AbH2;
  TGCD.ChemicalRate[i][k++] = TGCHEM_CHI_H2 * TGCD.NH * rate;

  // H+

  // Creation

  i++;
  j = 0;
  k = 0;

  // H + e -> H+ + 2e
  TGCD.CrRate[i][j++] = rate = TGCD.ChemRate[4] * TGCD.AbHI * TGCD.AbE * TGCD.NH;

  // Destruction

  j = 0;

  // H+ + 2e -> H + e + ph
  k1 = TGCD.ChemRate[5] * TGCD.AbHII * TGCD.AbE * TGCD.NH;

  // H+ + 2e -> H + e + ph
  k2 = (TGCD.ChemRate[6] * TGCD.AbE * TGCD.NH) * TGCD.AbE * TGCD.AbE * TGCD.NH;

  TGCD.DsRate[i][j++] = rate = -pow(k1, TGCD.TransPower) * pow(k2, 1. - TGCD.TransPower);

  for(i = 0; i < TGCHEM_NUM_ABUNDANCES; i++)
    for(j = 0; j < TGCHEM_NUM_RATES; j++)
      {
        TGCD.TotCrRate[i] += TGCD.CrRate[i][j];
        TGCD.TotDsRate[i] += TGCD.DsRate[i][j];
        TGCD.TotChemicalRate[i] += TGCD.ChemicalRate[i][j];
      }

  TGCD.TotCrRate[0] = TGCD.TotDsRate[0] = TGCD.TotChemicalRate[0] = 0.;

  if(TGCD.ChemMode || TGCD.EqH2Flag)
    TGCD.TotCrRate[1] = TGCD.TotDsRate[1] = TGCD.TotChemicalRate[1] = 0.;

  if(TGCD.EqHIIFlag)
    TGCD.TotCrRate[2] = TGCD.TotDsRate[2] = TGCD.TotChemicalRate[2] = 0.;

  for(i = 0; i < TGCHEM_NUM_ABUNDANCES; i++)
    {
      NV_Ith_S(dspecies, i) = TGCD.TotCrRate[i] + TGCD.TotDsRate[i];

      TGCD.ChemicalRateSum += TGCD.TotChemicalRate[i];
    }
}
