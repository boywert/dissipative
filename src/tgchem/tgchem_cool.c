/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/tgchem/tgchem_cool.c
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


void tgchem_cool(N_Vector dspecies)
{
  int i;
  double n0_rate, lte_rate, x, fesc, tau;

  TGCD.CoolingRateSum = 0.;

  for(i = 0; i < TGCHEM_NUM_COOLING; i++)
    TGCD.CoolingRate[i] = 0.;

#ifdef HEALRAY
  if(HRD.SourceFlag < 0)
#endif
    {
      // H2 ro-vibrational cooling
      n0_rate = TGCD.CoolRate[0] * TGCD.AbH2 * TGCD.AbHI * TGCD.NH * TGCD.NH;
      lte_rate = TGCD.CoolRate[1] * TGCD.AbH2 * TGCD.NH;

      if(n0_rate)
        TGCD.CoolingRate[0] = lte_rate / (1. + lte_rate / n0_rate);

      // Escape fraction: Fitting function (based on Ripamonti & Abel 2004)
      if(TGCD.ChemH2Mode == 0)
        fesc = tgchem_h2_line_fit_fesc();

      // Escape fraction: Soboloev method (Yoshida et al. 2006, Clark et al. 2011)
      if(TGCD.ChemH2Mode == 1)
        fesc = tgchem_h2_line_sob_fesc();

      // Escape fraction: Optically thin
      if(TGCD.ChemH2Mode == 2)
        fesc = 1.;

      TGCD.CoolingRate[0] *= fesc;

      TGCD.H2EscFrac = fesc;
    }

  // H2 collision-induced emission
  TGCD.CoolingRate[1] = TGCD.CoolRate[2] * TGCD.AbH2 * TGCD.NH * TGCD.NH;

  // Escape fraction: Fitting function (Clark et al. 2011)
  x = TGCD.NH / TGCD.CIEOptThickNHThresh;

  tau = dmax(pow(x, 2.8), 1e-10);

  fesc = (1. - exp(-tau)) / tau;

  TGCD.CoolingRate[1] *= fesc;

  // HI electronic excitation cooling
  TGCD.CoolingRate[2] = TGCD.CoolRate[3] * TGCD.AbHI * TGCD.AbE * TGCD.NH * TGCD.NH;

  // Escape fraction: we need something better than this!
  x = TGCD.NH / 1e7;

  fesc = exp(-x);

  TGCD.CoolingRate[2] *= fesc;

  // H- continuum cooling
  TGCD.CoolingRate[3] = TGCHEM_CHI_HM * TGCD.NH * TGCD.ChemRate[0] * TGCD.AbHI * TGCD.AbE * TGCD.NH;

  // Escape fraction: Ad-hoc exponential cut-off to reproduce Omukai 2001
  x = TGCD.NH / TGCD.CIEOptThickNHThresh;

  fesc = exp(-x);

  TGCD.CoolingRate[3] *= fesc;

  // Inverse Compton cooling with the CMB (Peebles 1971)
  TGCD.CoolingRate[4] = 5.65e-36 * pow(TGCD.TempCMB, 4) * (TGCD.Temp - TGCD.TempCMB) * TGCD.AbE * TGCD.NH;

  if(TGCD.ChemMode)
    TGCD.CoolingRate[0] = TGCD.CoolingRate[1] = 0.;

  for(i = 0; i < TGCHEM_NUM_COOLING; i++)
    TGCD.CoolingRateSum += TGCD.CoolingRate[i];

  TGCD.HeatRate = TGCD.ChemicalRateSum - TGCD.CoolingRateSum + TGCD.HRHeatRate;

  NV_Ith_S(dspecies, TGCHEM_NUM_ABUNDANCES) = TGCD.HeatRate;
}


double tgchem_h2_line_fit_fesc()
{
  double x, fesc;

  x = TGCD.NH / TGCD.H2OptThickNHThresh;

  if(x >= 1.)
    fesc = TGCD.H2OptThickConst * x / (pow(x, TGCD.H2OptThickConst) + TGCD.H2OptThickConst - 1.);
  else
    fesc = 1.;

  return fesc;
}


double tgchem_h2_line_sob_fesc()
{
  double fesc = dmin(TGCD.H2SobEmiss / TGCD.CoolRate[1], 1.);

  return fesc;
}
