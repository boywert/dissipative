/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/tgchem/tgchem_utils.c
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

void tgchem_eq_flag();
void tgchem_jeans();
void tgchem_h2_column();


int tgchem_rates(double time, N_Vector species, N_Vector dspecies, void *user_data)
{
  tgchem_debug_rates(0, species, dspecies);

  tgchem_check_for_nan(0, species);

  tgchem_compute_vars(0, species);

  tgchem_debug_rates(1, species, dspecies);

  tgchem_compute_rates();

  tgchem_chem(dspecies);

  tgchem_cool(dspecies);

  tgchem_debug_rates(2, species, dspecies);

  TGCD.NumRateCalls++;

  return 0;
}


void tgchem_check_for_nan(int mode, N_Vector species)
{
  int i;

  for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
    if(NV_Ith_S(species, i) != NV_Ith_S(species, i))
      terminate("Mode %d: Species %d is NaN! ID: %d", mode, i, TGCD.ID);
}


void tgchem_compute_vars(int mode, N_Vector species)
{
  int i;
  double val, x, abn, ab, gvar, gvar_H2;

  for(i = 0; i < TGCHEM_NUM_ABUNDANCES; i++)
    {
      val = dmin(dmax(NV_Ith_S(species, i), 0), TGCD.AbMax[i]);

      if(i == 0)
        TGCD.AbHM = val;

      if(i == 1)
        TGCD.AbH2 = val;

      if(i == 2)
        TGCD.AbHII = val;
    }

  if(!mode)
    TGCD.Energy = dmax(NV_Ith_S(species, TGCHEM_NUM_ABUNDANCES), 0.);

  TGCD.AbHI = dmax(1. - 2. * TGCD.AbH2 - TGCD.AbHII, 0.);

  TGCD.AbE = TGCD.AbHII;

  TGCD.NTot = dmax(1. + HE_ABUND - TGCD.AbH2 + TGCD.AbE, 0.) * TGCD.NH;
  //TGCD.NTot = (1. + HE_ABUND) * TGCD.NH;

  TGCD.Mu = (1. + 4. * HE_ABUND) / dmax(1. + HE_ABUND - TGCD.AbH2 + TGCD.AbHII, 0.);
  //TGCD.Mu = (1. + 4. * HE_ABUND) / (1. + HE_ABUND);

  if(!mode)
    TGCD.Temp = dmin(dmax(GAMMA_MINUS1 * TGCD.Energy / BOLTZMANN / TGCD.NTot, TGCHEM_TEMP_MIN), TGCHEM_TEMP_MAX);

  x = 6.1e3 / TGCD.Temp;

  abn = dmax(1. + HE_ABUND - TGCD.AbH2 + TGCD.AbE, 0.);
  ab = dmax(1. + HE_ABUND - 2 * TGCD.AbH2 + TGCD.AbE, 0.);

  gvar = 1. / GAMMA_MINUS1;
  gvar_H2 = 5. / 2. + pow(x, 2) * exp(x) / pow(exp(x) - 1., 2);

  TGCD.Gamma = 1. + (abn / (ab * gvar + TGCD.AbH2 * gvar_H2));
  //TGCD.Gamma = GAMMA;

  if(mode)
    TGCD.Energy = NV_Ith_S(TGCD.Species, TGCHEM_NUM_ABUNDANCES) = BOLTZMANN * TGCD.Temp * TGCD.NTot / TGCD.Gamma;
  else
    TGCD.Temp = dmin(dmax((TGCD.Gamma - 1.) * TGCD.Energy / BOLTZMANN / TGCD.NTot, TGCHEM_TEMP_MIN), TGCHEM_TEMP_MAX);
}


void tgchem_compute_rates()
{
  int i;

  TGCHEM_TEMP();

  for(i = 0; i < TGCHEM_NUM_CHEM; i++)
    TGCHEM_CHEMRATE(i);

  for(i = 0; i < TGCHEM_NUM_COOL; i++)
    TGCHEM_COOLRATE(i);

  if(TGCD.ChemH2Mode == 1)
    TGCHEM_H2_SOBEMISS();
}


void tgchem_compute_aux()
{
  tgchem_eq_flag();

  tgchem_jeans();

  tgchem_h2_column();

  tgchem_photo();
}


void tgchem_eq_flag()
{
  TGCD.EqH2Flag = 0;
  TGCD.EqHIIFlag = 0;

  if(!TGCD.ChemMode && TGCD.NH > TGCD.ChemH2Thresh)
    TGCD.EqH2Flag = 1;

  if(TGCD.NH > TGCD.ChemHIIThresh)
    TGCD.EqHIIFlag = 2;

  if(TGCD.EqH2Flag || TGCD.EqHIIFlag)
    TGCD.NumEq++;
  else
    TGCD.NumNeq++;
}


void tgchem_jeans()
{
  TGCD.Csnd = sqrt(TGCD.Gamma * BOLTZMANN * TGCD.Temp / TGCD.Mu / PROTONMASS);
  TGCD.VThermal = sqrt(BOLTZMANN * TGCD.Temp / PROTONMASS);
  TGCD.JeansLength = TGCD.Csnd * sqrt(3. * M_PI / 32. / GRAVITY / TGCD.Density);
}


void tgchem_h2_column()
{
  double lsob;

  if(TGCD.ChemH2Mode == 1)
    {
      if(TGCD.DivVel)
        lsob = dmin(TGCD.VThermal / dabs(TGCD.DivVel), TGCD.JeansLength);
      else
        lsob = TGCD.JeansLength;

      TGCD.H2Column = dmin(dmax(TGCD.NH * lsob, TGCHEM_H2_COLUMN_MIN), TGCHEM_H2_COLUMN_MAX);

      TGCHEM_H2_COLUMN();
    }
}


int tgchem_check_species(N_Vector species)
{
  int i;

  tgchem_check_for_nan(1, TGCD.Species);

  for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
    {
      if(NV_Ith_S(TGCD.Species, i) < 0.)
        {
          if(NV_Ith_S(TGCD.Species, i) < -NV_Ith_S(TGCD.SpeciesTol, i))
            return -1;

          NV_Ith_S(TGCD.Species, i) = 0.;
        }
    }

  for(i = 0; i < TGCHEM_NUM_ABUNDANCES; i++)
    if(NV_Ith_S(TGCD.Species, i) > TGCD.AbMax[i])
      {
        if((NV_Ith_S(TGCD.Species, i) - TGCD.AbMax[i]) / NV_Ith_S(TGCD.Species, i) > TGCHEM_TOL)
          return -1;

        NV_Ith_S(TGCD.Species, i) = TGCD.AbMax[i];
      }

  return 0;
}


void tgchem_debug_rates(int mode, N_Vector species, N_Vector dspecies)
{
  int i, j;

  //if(TGCD.DebugFlag > -1)
  if(0)
    //if(TGCD.NumRateCalls == 0 && (TGCD.Index == TGCD.DebugFlag - 1 || TGCD.Index == TGCD.DebugFlag))
#ifdef TGCHEM_TEST
    if(TGCD.Index == TGCD.DebugFlag)
#else
    if(TGCD.ID == TGCD.DebugFlag)
#endif
      if(TGCD.NumRateCalls % 100 == 0)
        {
          if(mode == 0)
            {
              printf("\n");

              printf("NumRateCalls: %d\n", TGCD.NumRateCalls);

              printf("Species: ");

              for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
                printf("%g ", NV_Ith_S(species, i));

              printf("\n");
            }
          else if(mode == 1)
            printf("Abundances: %g %g %g %g %g\n", TGCD.AbHM, TGCD.AbH2, TGCD.AbHII, TGCD.Energy, TGCD.Temp);
          else
            {
              printf("Rates:\n");

              for(i = 0; i < TGCHEM_NUM_ABUNDANCES; i++)
                {
                  printf("\t Abundance %d:\n", i);

                  for(j = 0; j < TGCHEM_NUM_RATES; j++)
                    if(TGCD.CrRate[i][j])
                      printf("\t\t Creation Rate %d: %g\n", j, TGCD.CrRate[i][j]);

                  for(j = 0; j < TGCHEM_NUM_RATES; j++)
                    if(TGCD.DsRate[i][j])
                      printf("\t\t Destruction Rate %d: %g\n", j, TGCD.DsRate[i][j]);

                  for(j = 0; j < TGCHEM_NUM_RATES; j++)
                    if(TGCD.ChemicalRate[i][j])
                      printf("\t\t Chemical Rate %d: %g\n", j, TGCD.ChemicalRate[i][j]);
                }

              for(i = 0; i < TGCHEM_NUM_COOLING; i++)
                if(TGCD.CoolingRate[i])
                  printf("CoolingRate %d: %g\n", i, TGCD.CoolingRate[i]);

              printf("Total Rates: %g %g %g\n", TGCD.ChemicalRateSum, TGCD.CoolingRateSum, TGCD.HRHeatRate);

              printf("DSpecies: ");

              for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
                printf("%g ", NV_Ith_S(dspecies, i));

              printf("\n");

              //terminate("");
            }
        }
}


void tgchem_get_opac()
{
  int rho_idx, temp_idx, idx1, idx2;
  double log_rho, log_temp, log_opac, drho, dtemp;
  double opac1, opac2, opac3, opac4;

  log_rho = log10(TGCD.Density);
  rho_idx = log_rho + 16;

  if(rho_idx < 0)
    {
      rho_idx = 0;
      drho = 0;
    }
  else if(rho_idx > TGCHEM_OPAC_NUM_RHO - 1)
    {
      rho_idx = TGCHEM_OPAC_NUM_RHO - 1;
      drho = 0;
    }
  else
    drho = log_rho - (rho_idx - 16);

  log_temp = log10(TGCD.Temp);
  temp_idx = (log_temp - 1.8) / 0.1;

  if(temp_idx < 0)
    {
      temp_idx = 0;
      dtemp = 0;
    }
  else if(temp_idx > TGCHEM_OPAC_NUM_TEMP - 1)
    {
      temp_idx = TGCHEM_OPAC_NUM_TEMP - 1;
      dtemp = 0;
    }
  else
    dtemp = (log_temp - (1.8 + 0.1 * temp_idx)) / 0.1;

  idx1 = TGCHEM_OPAC_NUM_RHO * temp_idx + rho_idx;

  if(rho_idx < TGCHEM_OPAC_NUM_RHO - 1)
    idx2 = idx1 + 1;
  else
    idx2 = idx1;

  opac1 = TGCD.OpacTable[idx1];
  opac2 = TGCD.OpacTable[idx2];

  opac3 = opac1 + drho * (opac2 - opac1);

  if(temp_idx < TGCHEM_OPAC_NUM_TEMP - 1)
    {
      idx1 = TGCHEM_OPAC_NUM_RHO * (temp_idx + 1) + rho_idx;

      if(rho_idx < TGCHEM_OPAC_NUM_RHO - 1)
        idx2 = idx1 + 1;
      else
        idx2 = idx1;

      opac1 = TGCD.OpacTable[idx1];
      opac2 = TGCD.OpacTable[idx2];

      opac4 = opac1 + drho * (opac2 - opac1);
    }
  else
    opac4 = opac3;

  log_opac = opac3 + dtemp * (opac4 - opac3);

  TGCD.Opac = pow(10., log_opac);
}
