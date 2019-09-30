/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/tgchem/tgchem_init.c
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

struct LineData_struct
{
  int j;
  int k;
  int l;
  int n;
  int line_idx;
  double emiss;
};

struct LineData_struct *LineData;

void tgchem_alloc();
void tgchem_h2_fesc_out();
void tgchem_h2_tables();
void tgchem_init_opac();
int compare_emiss(const void *a, const void *b);

void tgchem_begrun()
{
#ifndef TGCHEM_TEST
  CPU_Step[CPU_MISC] += measure_time();

  MPI_Bcast(&TGCD, sizeof(struct TGCD_struct), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif

  int i = 0;

  TGCD.AbMax[i++] = 1.;
  TGCD.AbMax[i++] = 0.5;
  TGCD.AbMax[i++] = 1.;

  TGCD.Species = N_VNew_Serial(TGCHEM_NUM_SPECIES);
  TGCD.SpeciesSave = N_VNew_Serial(TGCHEM_NUM_SPECIES);
  TGCD.SpeciesStepSave = N_VNew_Serial(TGCHEM_NUM_SPECIES);
  TGCD.DSpecies = N_VNew_Serial(TGCHEM_NUM_SPECIES);
  TGCD.SpeciesTol = N_VNew_Serial(TGCHEM_NUM_SPECIES);

  i = 0;

  NV_Ith_S(TGCD.SpeciesTol, i++) = TGCHEM_TOL_ABHM;
  NV_Ith_S(TGCD.SpeciesTol, i++) = TGCHEM_TOL_ABH2;
  NV_Ith_S(TGCD.SpeciesTol, i++) = TGCHEM_TOL_ABHII;
  NV_Ith_S(TGCD.SpeciesTol, i++) = TGCHEM_TOL_ABENERGY;

  //tgchem_init_cvode();
  /*
     TGCD.CVODEMem = CVodeCreate(CV_BDF, CV_NEWTON);
     CVodeInit(TGCD.CVODEMem, tgchem_rates, 0, TGCD.Species);
     CVodeSVtolerances(TGCD.CVODEMem, TGCHEM_TOL, TGCD.SpeciesTol);
     CVodeSetErrFile(TGCD.CVODEMem, NULL);
     CVDense(TGCD.CVODEMem, TGCHEM_NUM_SPECIES);
   */
  tgchem_alloc();

  tgchem_init_rates();

#ifndef TGCHEM_TEST
  CPU_Step[CPU_TGCHEM] += measure_time();
#endif
}


void tgchem_init_cvode()
{
  TGCD.CVODEMem = CVodeCreate(CV_BDF, CV_NEWTON);
  CVodeInit(TGCD.CVODEMem, tgchem_rates, 0, TGCD.Species);
  CVodeSVtolerances(TGCD.CVODEMem, TGCHEM_TOL, TGCD.SpeciesTol);
  CVodeSetErrFile(TGCD.CVODEMem, NULL);
  CVDense(TGCD.CVODEMem, TGCHEM_NUM_SPECIES);
}


void tgchem_finish_cvode()
{
  CVodeFree(&TGCD.CVODEMem);
}


void tgchem_init_rates()
{
  int i, j, k, l, m, n, p, count, target_i, line_idx;
  int num_temp, num_h2, num_eq, num_chem, num_cool;
  double temp, tgchem_temp_min, tgchem_temp_max, ln_temp_ev;
  double Z_H2, Z_HI, dE, f1, f2, g1, g2, nu_D, emiss, cross, tau, fesc, norm_fac;
  char buf[200];
  FILE *file;

  // Debugging

  target_i = 1650;

  // Constants

  // H2 ro-vibrational cooling optical depth correction (based on Ripamonti & Abel 2004)
  TGCD.H2OptThickConst = 1.45;
  TGCD.H2OptThickNHThresh = HYDROGEN_MASSFRAC * 8e9 / pow(TGCD.H2OptThickConst, 1. / (TGCD.H2OptThickConst - 1.));

  // H2 collision-induced emission optical depth correction (Clark et al. 2011)
  TGCD.CIEOptThickNHThresh = 1.5e16;

  // Transition to HII equilibrium recombination rate
  TGCD.ChemHIITrans = 1e17;

  // H2 dissociation by LW radiation (Abel et al 1997)
  TGCD.LWDissRate = 1.38e-12;

  // H2 energy levels and spontaneous emission coefficients (Borysow et al. 1989; Turner et al. 1977)
  tgchem_h2_tables();

  // Planck mean opacities (Mayer & Duschl 2005)
  tgchem_init_opac();

  // Temperature-dependent tables

  for(i = 0; i < TGCHEM_NUM_TEMP; i++)
    TGCD.TempTable[i] = pow(10., log10(TGCHEM_TEMP_MIN) + i * TGCHEM_LOG_DTEMP);

  for(i = 0; i < TGCHEM_H2_NUM_COLUMN; i++)
    TGCD.H2ColumnTable[i] = pow(10., log10(TGCHEM_H2_COLUMN_MIN) + i * TGCHEM_H2_LOG_DCOLUMN);

  LineData = mymalloc_movable(&LineData, "LineData", TGCHEM_H2_TOT_NUM_LINES * sizeof(struct LineData_struct));

  memset(TGCD.EqTable, 0, TGCHEM_NUM_EQ * TGCHEM_NUM_TEMP * sizeof(double));
  memset(TGCD.ChemTable, 0, TGCHEM_NUM_CHEM * TGCHEM_NUM_TEMP * sizeof(double));
  memset(TGCD.CoolTable, 0, TGCHEM_NUM_COOL * TGCHEM_NUM_TEMP * sizeof(double));
  memset(TGCD.H2EmissTable, 0, TGCHEM_H2_TOT_NUM_LINES * TGCHEM_NUM_TEMP * sizeof(double));
  memset(TGCD.H2CrossTable, 0, TGCHEM_H2_TOT_NUM_LINES * TGCHEM_NUM_TEMP * sizeof(double));

  memset(TGCD.H2TotEmissTable, 0, TGCHEM_NUM_TEMP * sizeof(double));
  memset(TGCD.H2TotCrossTable, 0, TGCHEM_NUM_TEMP * sizeof(double));
  memset(TGCD.H2SobEmissTable, 0, TGCHEM_NUM_TEMP * TGCHEM_H2_NUM_COLUMN * sizeof(double));

  for(i = 0; i < TGCHEM_NUM_TEMP; i++)
    {
      temp = TGCD.TempTable[i];

      ln_temp_ev = log(BOLTZMANN * temp / ELECTRON_VOLT);

      // Partition functions

      for(j = Z_H2 = 0; j < TGCHEM_H2_NUM_V; j++)
        for(k = 0; k < TGCHEM_H2_NUM_J; k++)
          Z_H2 += (2. * k + 1.) * exp(-TGCD.H2Energy[TGCHEM_H2_NUM_J * j + k] / BOLTZMANN / temp);

      for(j = Z_HI = 0; j <= 5; j++)
        Z_HI += (2. * j * j) * exp(-TGCHEM_CHI_HI / j / j / BOLTZMANN / temp);

      // H2 emission and absorption

      for(j = count = 0; j < TGCHEM_H2_NUM_V; j++)
        for(k = 0; k < TGCHEM_H2_NUM_J; k++)
          for(l = 0; l <= j; l++)
            for(m = 0; m < TGCHEM_H2_NUM_T; m++)
              {
                line_idx = (TGCHEM_H2_NUM_J * TGCHEM_H2_NUM_V * TGCHEM_H2_NUM_T) * j + (TGCHEM_H2_NUM_V * TGCHEM_H2_NUM_T) * k + TGCHEM_H2_NUM_T * l + m;

                if(TGCD.H2SpontCoeff[line_idx])
                  {
                    if(m == 0)
                      n = k - 2;
                    if(m == 1)
                      n = k;
                    if(m == 2)
                      n = k + 2;

                    dE = TGCD.H2Energy[TGCHEM_H2_NUM_J * j + k] - TGCD.H2Energy[TGCHEM_H2_NUM_J * l + n];

                    g1 = (2. * n + 1.);
                    g2 = (2. * k + 1.);

                    f1 = g1 * exp(-TGCD.H2Energy[TGCHEM_H2_NUM_J * l + n] / BOLTZMANN / temp) / Z_H2;
                    f2 = g2 * exp(-TGCD.H2Energy[TGCHEM_H2_NUM_J * j + k] / BOLTZMANN / temp) / Z_H2;

                    nu_D = dE / PLANCK / CLIGHT * sqrt(BOLTZMANN * temp / PROTONMASS);

                    emiss = dE * f2 * TGCD.H2SpontCoeff[line_idx];

                    cross = pow(PLANCK * CLIGHT / dE, 2) / 8. / M_PI * g2 / g1 * TGCD.H2SpontCoeff[line_idx] * (1. - exp(-dE / BOLTZMANN / temp)) / sqrt(M_PI) / nu_D;

                    TGCD.H2EmissTable[line_idx * TGCHEM_NUM_TEMP + i] = emiss;
                    TGCD.H2CrossTable[line_idx * TGCHEM_NUM_TEMP + i] = f1 * cross;

                    TGCD.H2TotEmissTable[i] += emiss;
                    TGCD.H2TotCrossTable[i] += emiss * f1 * cross;

                    if(TGCD.ChemH2Mode == 1)
                      for(p = 0; p < TGCHEM_H2_NUM_COLUMN; p++)
                        {
                          tau = dmax(f1 * cross * TGCD.H2ColumnTable[p], 1e-10);

                          fesc = (1. - exp(-tau)) / tau;

                          TGCD.H2SobEmissTable[i * TGCHEM_H2_NUM_COLUMN + p] += emiss * fesc;

                          //if(i == target_i && (p == 150 || p == 151))
                          //printf("fesc: %d %g %d %d %d %d %d %g %g %g %g %g %g %g %g\n", i, temp, j, k, l, n, p, TGCD.H2ColumnTable[p], f1, cross, tau, fesc, emiss, emiss * fesc, TGCD.H2SobEmissTable[i * TGCHEM_H2_NUM_COLUMN + p]);
                        }

                    LineData[count].j = j;
                    LineData[count].k = k;
                    LineData[count].l = l;
                    LineData[count].n = n;
                    LineData[count].line_idx = line_idx;
                    LineData[count].emiss = emiss;

                    //if(i == target_i)
                    //mpi_printf("%d %g %d %d %d %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g\n", i, temp, j, k, l, n, line_idx, dE, g1, g2, TGCD.H2Energy[TGCHEM_H2_NUM_J * l + n], TGCD.H2Energy[TGCHEM_H2_NUM_J * j + k], Z_H2, f1, f2, nu_D, TGCD.H2SpontCoeff[line_idx], emiss, cross, f1 * cross);

                    count++;
                  }
              }
#ifdef HEALRAY
      if(HRD.RayNumLines > count)
        terminate("HRD.RayNumLines too large!");

      TGCD.H2TotCrossTable[i] /= TGCD.H2TotEmissTable[i];

      qsort(LineData, count, sizeof(struct LineData_struct), compare_emiss);

      for(j = 0; j < HRD.RayNumLines; j++)
        TGCD.H2LineIdxTable[j * TGCHEM_NUM_TEMP + i] = LineData[j].line_idx;

      if(i == target_i)
        {
          double tot_emiss = 0;

          for(j = tot_emiss = 0; j < HRD.RayNumLines; j++)
            tot_emiss += TGCD.H2EmissTable[TGCD.H2LineIdxTable[j * TGCHEM_NUM_TEMP + i] * TGCHEM_NUM_TEMP + i];

          //mpi_printf("tot emiss: %d %g %d %g %g %g\n", i, temp, count, TGCD.H2TotEmissTable[i], tot_emiss, TGCD.H2TotCrossTable[i]);

          //for(j = 0; j < HRD.RayNumLines; j++)
          //printf("emiss: %d %d %d %d %d %d %g %g\n", j, LineData[j].j, LineData[j].k, LineData[j].l, LineData[j].n, TGCD.H2LineIdxTable[j * TGCHEM_NUM_TEMP + i], TGCD.H2EmissTable[TGCD.H2LineIdxTable[j * TGCHEM_NUM_TEMP + i] * TGCHEM_NUM_TEMP + i], TGCD.H2CrossTable[TGCD.H2LineIdxTable[j * TGCHEM_NUM_TEMP + i] * TGCHEM_NUM_TEMP + i]);
        }
#endif
      // Equilibrium abundances

      TGCD.EqTable[0 * TGCHEM_NUM_TEMP + i] = pow(Z_HI, 2) / Z_H2 * pow(M_PI * PROTONMASS * BOLTZMANN * temp / pow(PLANCK, 2), 3. / 2.) * exp(-TGCHEM_CHI_H2 / BOLTZMANN / temp);

      TGCD.EqTable[1 * TGCHEM_NUM_TEMP + i] = log(2. / Z_HI * pow(2. * M_PI * ELECTRONMASS * BOLTZMANN * temp / pow(PLANCK, 2), 3. / 2.)) - TGCHEM_CHI_HI / BOLTZMANN / temp;

      // Chemical rates

      // H + e -> H- + ph (Galli & Palla 1998)
      TGCD.ChemTable[0 * TGCHEM_NUM_TEMP + i] = 1.4e-18 * pow(temp, 0.928) * exp(-temp / 1.62e4);

      // H + H- -> H2 + e (Kreckel et al. 2010)
      TGCD.ChemTable[1 * TGCHEM_NUM_TEMP + i] =
        1.35e-9 * (pow(temp, 9.8493e-2) + 3.2852e-1 * pow(temp, 5.561e-1) + 2.771e-7 * pow(temp, 2.1826)) / (1 + 6.191e-3 * pow(temp, 1.0461) +
                                                                                                             8.9712e-11 * pow(temp, 3.0424) + 3.2576e-14 * pow(temp, 3.7741));

      // 3H -> H + H2 (Forrey 2013)
      TGCD.ChemTable[2 * TGCHEM_NUM_TEMP + i] = 6e-32 * pow(temp, -1. / 4.) + 2e-31 / sqrt(temp);

      // H + H2 -> 3H (Detailed balance with three-body formation rate)
      TGCD.ChemTable[3 * TGCHEM_NUM_TEMP + i] = TGCD.ChemTable[2 * TGCHEM_NUM_TEMP + i] * TGCD.EqTable[0 * TGCHEM_NUM_TEMP + i];

      // H + e -> H+ + 2e (Janev et al. 1987)
      TGCD.ChemTable[4 * TGCHEM_NUM_TEMP + i] =
        -32.71396786 + 13.536556 * ln_temp_ev - 5.73932875 * pow(ln_temp_ev, 2) + 1.56315498 * pow(ln_temp_ev, 3) - 0.2877056 * pow(ln_temp_ev,
                                                                                                                                    4) +
        3.48255977e-2 * pow(ln_temp_ev, 5) - 2.63197617e-3 * pow(ln_temp_ev, 6) + 1.11954395e-4 * pow(ln_temp_ev, 7) - 2.03914985e-6 * pow(ln_temp_ev, 8);

      // H+ + 2e -> H + e + ph (Case B recombination: Ferland et al. 1992)
      TGCD.ChemTable[5 * TGCHEM_NUM_TEMP + i] = 2.753e-14 * pow(3.15614e5 / temp, 3. / 2.) * pow(1 + pow(1.15188e5 / temp, 4.07e-1), -2.242);

      // H+ + 2e -> H + e + ph (Detailed balance with collisional ionization rate)
      TGCD.ChemTable[6 * TGCHEM_NUM_TEMP + i] = exp(TGCD.ChemTable[4 * TGCHEM_NUM_TEMP + i] - TGCD.EqTable[1 * TGCHEM_NUM_TEMP + i]);

      // Revert to exponential

      TGCD.EqTable[1 * TGCHEM_NUM_TEMP + i] = exp(TGCD.EqTable[1 * TGCHEM_NUM_TEMP + i]);
      TGCD.ChemTable[4 * TGCHEM_NUM_TEMP + i] = exp(TGCD.ChemTable[4 * TGCHEM_NUM_TEMP + i]);

      // Cooling rates

      // H2 n -> 0 cooling (Galli & Palla 1998)
      TGCD.CoolTable[0 * TGCHEM_NUM_TEMP + i] = pow(10, -103 + 97.59 * log10(temp) - 48.05 * pow(log10(temp), 2) + 10.8 * pow(log10(temp), 3) - 0.9032 * pow(log10(temp), 4));

      // H2 LTE cooling
      TGCD.CoolTable[1 * TGCHEM_NUM_TEMP + i] = TGCD.H2TotEmissTable[i];

      // H2 collision-induced emission (Ripamonti & Abel 2004)
      TGCD.CoolTable[2 * TGCHEM_NUM_TEMP + i] = 5.3e-49 * pow(temp, 4);

      // HI electronic excitation cooling (Cen 1992)
      TGCD.CoolTable[3 * TGCHEM_NUM_TEMP + i] = 7.5e-19 * exp(-1.18348e5 / temp) / (1. + sqrt(temp / 1e5));
    }

  myfree_movable(LineData);

  for(i = 0; i < TGCHEM_NUM_TEMP - 1; i++)
    {
      for(j = 0; j < TGCHEM_NUM_EQ; j++)
        TGCD.EqTable[j * TGCHEM_NUM_TEMP + i] = dmin(TGCD.EqTable[j * TGCHEM_NUM_TEMP + i], TGCHEM_MAX_RATE);

      for(j = 0; j < TGCHEM_NUM_CHEM; j++)
        TGCD.ChemTable[j * TGCHEM_NUM_TEMP + i] = dmin(TGCD.ChemTable[j * TGCHEM_NUM_TEMP + i], TGCHEM_MAX_RATE);

      for(j = 0; j < TGCHEM_NUM_COOL; j++)
        TGCD.CoolTable[j * TGCHEM_NUM_TEMP + i] = dmin(TGCD.CoolTable[j * TGCHEM_NUM_TEMP + i], TGCHEM_MAX_RATE);
    }

  for(i = 0; i < TGCHEM_NUM_TEMP - 1; i++)
    {
      for(j = 0; j < TGCHEM_NUM_EQ; j++)
        if(TGCD.EqTable[j * TGCHEM_NUM_TEMP + i] != TGCD.EqTable[j * TGCHEM_NUM_TEMP + i])
          terminate("Encountered Nan in rate initialization: %d %d", i, j);

      for(j = 0; j < TGCHEM_NUM_CHEM; j++)
        if(TGCD.ChemTable[j * TGCHEM_NUM_TEMP + i] != TGCD.ChemTable[j * TGCHEM_NUM_TEMP + i])
          terminate("Encountered Nan in rate initialization: %d %d", i, j);

      for(j = 0; j < TGCHEM_NUM_COOL; j++)
        if(TGCD.CoolTable[j * TGCHEM_NUM_TEMP + i] != TGCD.CoolTable[j * TGCHEM_NUM_TEMP + i])
          terminate("Encountered Nan in rate initialization: %d %d", i, j);
    }

  for(i = 0; i < TGCHEM_NUM_TEMP - 1; i++)
    {
      for(j = 0; j < TGCHEM_NUM_EQ; j++)
        TGCD.DEqTable[j * TGCHEM_NUM_TEMP + i] = (TGCD.EqTable[j * TGCHEM_NUM_TEMP + i + 1] - TGCD.EqTable[j * TGCHEM_NUM_TEMP + i]) / (TGCD.TempTable[i + 1] - TGCD.TempTable[i]);

      for(j = 0; j < TGCHEM_NUM_CHEM; j++)
        TGCD.DChemTable[j * TGCHEM_NUM_TEMP + i] = (TGCD.ChemTable[j * TGCHEM_NUM_TEMP + i + 1] - TGCD.ChemTable[j * TGCHEM_NUM_TEMP + i]) / (TGCD.TempTable[i + 1] - TGCD.TempTable[i]);

      for(j = 0; j < TGCHEM_NUM_COOL; j++)
        TGCD.DCoolTable[j * TGCHEM_NUM_TEMP + i] = (TGCD.CoolTable[j * TGCHEM_NUM_TEMP + i + 1] - TGCD.CoolTable[j * TGCHEM_NUM_TEMP + i]) / (TGCD.TempTable[i + 1] - TGCD.TempTable[i]);

      for(j = 0; j < TGCHEM_H2_TOT_NUM_LINES; j++)
        TGCD.H2DEmissTable[j * TGCHEM_NUM_TEMP + i] = (TGCD.H2EmissTable[j * TGCHEM_NUM_TEMP + i + 1] - TGCD.H2EmissTable[j * TGCHEM_NUM_TEMP + i]) / (TGCD.TempTable[i + 1] - TGCD.TempTable[i]);

      TGCD.H2DTotEmissTable[i] = (TGCD.H2TotEmissTable[i + 1] - TGCD.H2TotEmissTable[i]) / (TGCD.TempTable[i + 1] - TGCD.TempTable[i]);

      for(j = 0; j < TGCHEM_H2_TOT_NUM_LINES; j++)
        TGCD.H2DCrossTable[j * TGCHEM_NUM_TEMP + i] = (TGCD.H2CrossTable[j * TGCHEM_NUM_TEMP + i + 1] - TGCD.H2CrossTable[j * TGCHEM_NUM_TEMP + i]) / (TGCD.TempTable[i + 1] - TGCD.TempTable[i]);

      TGCD.H2DTotCrossTable[i] = (TGCD.H2TotCrossTable[i + 1] - TGCD.H2TotCrossTable[i]) / (TGCD.TempTable[i + 1] - TGCD.TempTable[i]);
    }

  for(j = 0; j < TGCHEM_NUM_EQ; j++)
    TGCD.DEqTable[j * TGCHEM_NUM_TEMP + TGCHEM_NUM_TEMP - 1] = TGCD.DEqTable[j * TGCHEM_NUM_TEMP + TGCHEM_NUM_TEMP - 2];

  for(j = 0; j < TGCHEM_NUM_CHEM; j++)
    TGCD.DChemTable[j * TGCHEM_NUM_TEMP + TGCHEM_NUM_TEMP - 1] = TGCD.DChemTable[j * TGCHEM_NUM_TEMP + TGCHEM_NUM_TEMP - 2];

  for(j = 0; j < TGCHEM_NUM_COOL; j++)
    TGCD.DCoolTable[j * TGCHEM_NUM_TEMP + TGCHEM_NUM_TEMP - 1] = TGCD.DCoolTable[j * TGCHEM_NUM_TEMP + TGCHEM_NUM_TEMP - 2];

  for(j = 0; j < TGCHEM_H2_TOT_NUM_LINES; j++)
    TGCD.H2DEmissTable[j * TGCHEM_NUM_TEMP + TGCHEM_NUM_TEMP - 1] = TGCD.H2DEmissTable[j * TGCHEM_NUM_TEMP + TGCHEM_NUM_TEMP - 2];

  TGCD.H2DTotEmissTable[TGCHEM_NUM_TEMP - 1] = TGCD.H2DTotEmissTable[TGCHEM_NUM_TEMP - 2];

  for(j = 0; j < TGCHEM_H2_TOT_NUM_LINES; j++)
    TGCD.H2DCrossTable[j * TGCHEM_NUM_TEMP + TGCHEM_NUM_TEMP - 1] = TGCD.H2DCrossTable[j * TGCHEM_NUM_TEMP + TGCHEM_NUM_TEMP - 2];

  TGCD.H2DTotCrossTable[TGCHEM_NUM_TEMP - 1] = TGCD.H2DTotCrossTable[TGCHEM_NUM_TEMP - 2];

  for(i = 0; i < TGCHEM_NUM_TEMP; i++)
    {
      for(j = 0; j < TGCHEM_H2_NUM_COLUMN - 1; j++)
        TGCD.H2DSobEmissTable[i * TGCHEM_H2_NUM_COLUMN + j] =
          (TGCD.H2SobEmissTable[i * TGCHEM_H2_NUM_COLUMN + j + 1] - TGCD.H2SobEmissTable[i * TGCHEM_H2_NUM_COLUMN + j]) / (TGCD.H2ColumnTable[j + 1] - TGCD.H2ColumnTable[j]);
    }

  for(i = 0; i < TGCHEM_NUM_TEMP; i++)
    TGCD.H2DSobEmissTable[i * TGCHEM_H2_NUM_COLUMN + TGCHEM_H2_NUM_COLUMN - 1] = TGCD.H2DSobEmissTable[i * TGCHEM_H2_NUM_COLUMN + TGCHEM_H2_NUM_COLUMN - 2];

  if(TGCD.ChemH2Mode == 1)
    tgchem_h2_fesc_out();

#ifdef TGCHEM_TEST
  sprintf(buf, "../data/tgchem_rates.dat");
#else
  sprintf(buf, "data/tgchem_rates.dat");
#endif

  if(!(file = fopen(buf, "w")))
    terminate("Could not open file!\n");

  tgchem_temp_min = TGCHEM_TEMP_MIN;
  tgchem_temp_max = TGCHEM_TEMP_MAX;

  num_temp = TGCHEM_NUM_TEMP;
  num_eq = TGCHEM_NUM_EQ;
  num_chem = TGCHEM_NUM_CHEM;
  num_cool = TGCHEM_NUM_COOL;

  fwrite(&TGCD.ChemH2Mode, sizeof(int), 1, file);
  fwrite(&tgchem_temp_min, sizeof(double), 1, file);
  fwrite(&tgchem_temp_max, sizeof(double), 1, file);

  fwrite(&num_temp, sizeof(int), 1, file);
  fwrite(&num_eq, sizeof(int), 1, file);
  fwrite(&num_chem, sizeof(int), 1, file);
  fwrite(&num_cool, sizeof(int), 1, file);

  fwrite(TGCD.TempTable, sizeof(double), num_temp, file);

  fwrite(TGCD.EqTable, sizeof(double), num_temp * num_eq, file);
  fwrite(TGCD.ChemTable, sizeof(double), num_temp * num_chem, file);
  fwrite(TGCD.CoolTable, sizeof(double), num_temp * num_cool, file);

  if(TGCD.ChemH2Mode == 1)
    {
      fwrite(TGCD.H2TestSobXVal, sizeof(double), num_temp, file);
      fwrite(TGCD.H2TestFitEscFrac, sizeof(double), num_temp, file);
      fwrite(TGCD.H2TestSobEscFrac, sizeof(double), num_temp, file);
    }

  fclose(file);
}


void tgchem_alloc()
{
  TGCD.TempTable = mymalloc_movable(TGCD.TempTable, "TGCD.TempTable", TGCHEM_NUM_TEMP * sizeof(double));
  TGCD.EqTable = mymalloc_movable(TGCD.EqTable, "TGCD.EqTable", TGCHEM_NUM_EQ * TGCHEM_NUM_TEMP * sizeof(double));
  TGCD.DEqTable = mymalloc_movable(TGCD.DEqTable, "TGCD.DEqTable", TGCHEM_NUM_EQ * TGCHEM_NUM_TEMP * sizeof(double));
  TGCD.ChemTable = mymalloc_movable(TGCD.ChemTable, "TGCD.ChemTable", TGCHEM_NUM_CHEM * TGCHEM_NUM_TEMP * sizeof(double));
  TGCD.DChemTable = mymalloc_movable(TGCD.DChemTable, "TGCD.DChemTable", TGCHEM_NUM_CHEM * TGCHEM_NUM_TEMP * sizeof(double));
  TGCD.CoolTable = mymalloc_movable(TGCD.CoolTable, "TGCD.CoolTable", TGCHEM_NUM_COOL * TGCHEM_NUM_TEMP * sizeof(double));
  TGCD.DCoolTable = mymalloc_movable(TGCD.DCoolTable, "TGCD.DCoolTable", TGCHEM_NUM_COOL * TGCHEM_NUM_TEMP * sizeof(double));

  TGCD.H2SpontCoeff = mymalloc_movable(TGCD.H2SpontCoeff, "TGCD.H2SpontCoeff", TGCHEM_H2_TOT_NUM_LINES * sizeof(double));
  TGCD.H2ColumnTable = mymalloc_movable(TGCD.H2ColumnTable, "TGCD.H2ColumnTable", TGCHEM_H2_NUM_COLUMN * sizeof(double));
  TGCD.H2SobEmissTable = mymalloc_movable(TGCD.H2SobEmissTable, "TGCD.H2SobEmissTable", TGCHEM_NUM_TEMP * TGCHEM_H2_NUM_COLUMN * sizeof(double));
  TGCD.H2DSobEmissTable = mymalloc_movable(TGCD.H2DSobEmissTable, "TGCD.H2DSobEmissTable", TGCHEM_NUM_TEMP * TGCHEM_H2_NUM_COLUMN * sizeof(double));
  TGCD.H2TestSobXVal = mymalloc_movable(TGCD.H2TestSobXVal, "TGCD.H2TestSobXVal", TGCHEM_NUM_TEMP * sizeof(double));
  TGCD.H2TestFitEscFrac = mymalloc_movable(TGCD.H2TestFitEscFrac, "TGCD.H2TestFitEscFrac", TGCHEM_NUM_TEMP * sizeof(double));
  TGCD.H2TestSobEscFrac = mymalloc_movable(TGCD.H2TestSobEscFrac, "TGCD.H2TestSobEscFrac", TGCHEM_NUM_TEMP * sizeof(double));
  TGCD.H2EmissTable = mymalloc_movable(TGCD.H2EmissTable, "TGCD.H2EmissTable", TGCHEM_H2_TOT_NUM_LINES * TGCHEM_NUM_TEMP * sizeof(double));
  TGCD.H2DEmissTable = mymalloc_movable(TGCD.H2DEmissTable, "TGCD.H2DEmissTable", TGCHEM_H2_TOT_NUM_LINES * TGCHEM_NUM_TEMP * sizeof(double));
  TGCD.H2TotEmissTable = mymalloc_movable(TGCD.H2TotEmissTable, "TGCD.H2TotEmissTable", TGCHEM_NUM_TEMP * sizeof(double));
  TGCD.H2DTotEmissTable = mymalloc_movable(TGCD.H2DTotEmissTable, "TGCD.H2DTotEmissTable", TGCHEM_NUM_TEMP * sizeof(double));
  TGCD.H2CrossTable = mymalloc_movable(TGCD.H2CrossTable, "TGCD.H2CrossTable", TGCHEM_H2_TOT_NUM_LINES * TGCHEM_NUM_TEMP * sizeof(double));
  TGCD.H2DCrossTable = mymalloc_movable(TGCD.H2DCrossTable, "TGCD.H2DCrossTable", TGCHEM_H2_TOT_NUM_LINES * TGCHEM_NUM_TEMP * sizeof(double));
  TGCD.H2TotCrossTable = mymalloc_movable(TGCD.H2TotCrossTable, "TGCD.H2TotCrossTable", TGCHEM_NUM_TEMP * sizeof(double));
  TGCD.H2DTotCrossTable = mymalloc_movable(TGCD.H2DTotCrossTable, "TGCD.H2DTotCrossTable", TGCHEM_NUM_TEMP * sizeof(double));
#ifdef HEALRAY
  TGCD.H2LineIdxTable = mymalloc_movable(TGCD.H2LineIdxTable, "TGCD.H2LineIdxTable", HRD.RayNumLines * TGCHEM_NUM_TEMP * sizeof(int));
#endif

  TGCD.OpacTable = mymalloc_movable(TGCD.OpacTable, "TGCD.OpacTable", TGCHEM_OPAC_NUM_RHO * TGCHEM_OPAC_NUM_TEMP * sizeof(double));
}


void tgchem_h2_fesc_out()
{
  int i;
  double temp, nh_test, nh_min, nh_max, jeans_frac, fesc;
  double H2TotEmiss;

  nh_test = 1e9;
  jeans_frac = 0.1;

  nh_min = 1e9;
  nh_max = 1e18;

  TGCD.NH = nh_test;

  for(i = 0; i < TGCHEM_NUM_TEMP; i++)
    {
      //TGCD.Temp = TGCD.TempTable[i];
      TGCD.Temp = 1e3;

      TGCD.NH = pow(10., log10(nh_min) + i * log10(nh_max / nh_min) / TGCHEM_NUM_TEMP);
      TGCD.Density = TGCD.NH * PROTONMASS / HYDROGEN_MASSFRAC;

      TGCD.Csnd = sqrt(BOLTZMANN * TGCD.Temp / PROTONMASS);
      TGCD.JeansLength = TGCD.Csnd * sqrt(3. * M_PI / 32. / GRAVITY / TGCD.Density);

      TGCD.H2Column = dmin(dmax(TGCD.NH * jeans_frac * TGCD.JeansLength, TGCHEM_H2_COLUMN_MIN), TGCHEM_H2_COLUMN_MAX);

      TGCHEM_TEMP();
      TGCHEM_COOLRATE(1);
      TGCHEM_H2_COLUMN();
      TGCHEM_H2_SOBEMISS();

      TGCD.H2TestSobXVal[i] = TGCD.NH;
      TGCD.H2TestFitEscFrac[i] = tgchem_h2_line_fit_fesc();
      TGCD.H2TestSobEscFrac[i] = tgchem_h2_line_sob_fesc();
    }
}


void tgchem_h2_tables()
{
  int i, j, k, l, idx, line_idx;

  // H2 energy levels

  double min_energy;
  double coeff_v[TGCHEM_H2_ENERGY_SA];

  double prefac[TGCHEM_H2_ENERGY_SA] = { -1e5, 1., -1e-2, 1e-5, -1e-8, 1e-11, -1e-14, 1e-17 };

  double coeff[TGCHEM_H2_ENERGY_SA * TGCHEM_H2_ENERGY_SB] = { 3.8496e-1, -4.609e-2, 1.78e-3, -7e-5, 2.951e-6,
    5.4438e1, 4.6063, -2.005, 1.926e-1, -6.2953e-3,
    -7.2593e-1, 5.999, -1.5187, 1.2721e-1, -3.3391e-3,
    -1.2662e1, 2.0047e1, -4.7873, 3.69e-1, -9.003e-3,
    -2.4006e1, 3.3989e1, -7.6841, 5.2413e-1, -1.1297e-2,
    -2.2384e1, 3.115e1, -6.6139, 3.9226e-1, -9.496e-3,
    -1.0541e1, 1.4746e1, -2.9476, 1.6016e-1, -6.5005e-3,
    -2.0021, 2.84, -5.4654e-1, 3.1636e-2, -2.4398e-3
  };

  for(i = 0; i < TGCHEM_H2_NUM_V; i++)
    {
      for(j = 0; j < TGCHEM_H2_ENERGY_SA; j++)
        {
          coeff_v[j] = 0;

          for(k = 0; k < TGCHEM_H2_ENERGY_SB; k++)
            coeff_v[j] += coeff[j * TGCHEM_H2_ENERGY_SB + k] * pow(i + 1. / 2., k);
        }

      for(j = 0; j < TGCHEM_H2_NUM_J; j++)
        {
          idx = TGCHEM_H2_NUM_J * i + j;

          TGCD.H2Energy[idx] = 0;

          for(k = 0; k < TGCHEM_H2_ENERGY_SA; k++)
            TGCD.H2Energy[idx] += prefac[k] * coeff_v[k] * pow(j * (j + 1.), k);

          if(i == 0 && j == 0)
            min_energy = dabs(TGCD.H2Energy[idx]);

          TGCD.H2Energy[idx] += min_energy;

          TGCD.H2Energy[idx] *= PLANCK * CLIGHT;

          //mpi_printf("%d %d %g\n", i, j, TGCD.H2Energy[idx] / ELECTRON_VOLT);
        }
    }

  // H2 spontaneous emission coefficients

  // v0 -> v0

  double tmp0[TGCHEM_H2_NUM_J * TGCHEM_H2_NUM_T] = { 0., 0., 0.,
    0., 0., 0.,
    2.94e-11, 0., 0.,
    4.76e-10, 0., 0.,
    2.76e-9, 0., 0.,
    9.84e-9, 0., 0.,
    2.64e-8, 0., 0.,
    5.88e-8, 0., 0.,
    1.14e-7, 0., 0.,
    2e-7, 0., 0.,
    3.24e-7, 0., 0.,
    4.9e-7, 0., 0.,
    7.03e-7, 0., 0.,
    9.64e-7, 0., 0.,
    1.27e-6, 0., 0.,
    1.62e-6, 0., 0.,
    2e-6, 0., 0.,
    2.41e-6, 0., 0.,
    2.83e-6, 0., 0.,
    3.26e-6, 0., 0.
  };

  // v1 -> v0

  double tmp1[TGCHEM_H2_NUM_J * TGCHEM_H2_NUM_T] = { 0., 0., 8.54e-7,
    0., 4.29e-7, 4.23e-7,
    2.53e-7, 3.03e-7, 2.9e-7,
    3.47e-7, 2.78e-7, 2.09e-7,
    3.98e-7, 2.65e-7, 1.5e-7,
    4.21e-7, 2.55e-7, 1.06e-7,
    4.19e-7, 2.45e-7, 7.38e-8,
    3.96e-7, 2.34e-7, 4.98e-8,
    3.54e-7, 2.23e-7, 3.27e-8,
    2.98e-7, 2.12e-7, 2.07e-8,
    2.34e-7, 1.99e-7, 1.27e-8,
    1.68e-7, 1.87e-7, 7.44e-9,
    1.05e-7, 1.74e-7, 4.16e-9,
    5.3e-8, 1.61e-7, 2.19e-9,
    1.65e-8, 1.47e-7, 1.08e-9,
    4.26e-10, 1.34e-7, 4.87e-10,
    8.38e-9, 1.21e-7, 1.96e-10,
    4.27e-8, 1.08e-7, 6.71e-11,
    1.04e-7, 9.61e-8, 1.82e-11,
    1.93e-7, 8.43e-8, 3.38e-12
  };

  // v1 -> v1

  double tmp2[TGCHEM_H2_NUM_J * TGCHEM_H2_NUM_T] = { 0., 0., 0.,
    0., 0., 0.,
    2.79e-11, 0., 0.,
    4.5e-10, 0., 0.,
    2.59e-9, 0., 0.,
    9.21e-9, 0., 0.,
    2.46e-8, 0., 0.,
    5.44e-8, 0., 0.,
    1.05e-7, 0., 0.,
    1.82e-7, 0., 0.,
    2.92e-7, 0., 0.,
    4.38e-7, 0., 0.,
    6.21e-7, 0., 0.,
    8.42e-7, 0., 0.,
    1.1e-6, 0., 0.,
    1.38e-6, 0., 0.,
    1.69e-6, 0., 0.,
    2e-6, 0., 0.,
    2.32e-6, 0., 0.,
    2.64e-6, 0., 0.
  };

  // v2 -> v0

  double tmp3[TGCHEM_H2_NUM_J * TGCHEM_H2_NUM_T] = { 0., 0., 3.47e-7,
    0., 1.94e-7, 1.61e-7,
    1.27e-7, 1.38e-7, 1.03e-7,
    1.9e-7, 1.29e-7, 6.98e-8,
    2.38e-7, 1.25e-7, 4.72e-8,
    2.77e-7, 1.23e-7, 3.15e-8,
    3.07e-7, 1.21e-7, 2.05e-8,
    3.28e-7, 1.2e-7, 1.31e-8,
    3.39e-7, 1.18e-7, 8.07e-9,
    3.4e-7, 1.17e-7, 4.83e-9,
    3.3e-7, 1.15e-7, 2.79e-9,
    3.12e-7, 1.13e-7, 1.55e-9,
    2.85e-7, 1.11e-7, 8.27e-10,
    2.53e-7, 1.08e-7, 4.21e-10,
    2.16e-7, 1.05e-7, 2.04e-10,
    1.78e-7, 1.02e-7, 9.45e-11,
    1.39e-7, 9.88e-8, 4.23e-11,
    1.01e-7, 9.5e-8, 1.91e-11,
    6.8e-8, 9.07e-8, 9.63e-12,
    3.98e-8, 8.62e-8, 6.26e-12
  };

  // v2 -> v1

  double tmp4[TGCHEM_H2_NUM_J * TGCHEM_H2_NUM_T] = { 0., 0., 1.29e-6,
    0., 6.37e-7, 6.4e-7,
    3.68e-7, 4.5e-7, 4.41e-7,
    4.98e-7, 4.12e-7, 3.18e-7,
    5.6e-7, 3.91e-7, 2.28e-7,
    5.77e-7, 3.74e-7, 1.62e-7,
    5.57e-7, 3.58e-7, 1.12e-7,
    5.05e-7, 3.4e-7, 7.5e-8,
    4.3e-7, 3.22e-7, 4.88e-8,
    3.38e-7, 3.03e-7, 3.06e-8,
    2.41e-7, 2.84e-7, 1.85e-8,
    1.49e-7, 2.63e-7, 1.07e-8,
    7.2e-8, 2.42e-7, 5.84e-9,
    1.96e-8, 2.21e-7, 3e-9,
    1.49e-11, 2.01e-7, 1.43e-9,
    1.91e-8, 1.8e-7, 6.17e-10,
    8.07e-8, 1.6e-7, 2.33e-10,
    1.86e-7, 1.41e-7, 7.32e-11,
    3.35e-7, 1.22e-7, 1.73e-11,
    5.24e-7, 1.05e-7, 2.46e-12
  };

  // v2 -> v2

  double tmp5[TGCHEM_H2_NUM_J * TGCHEM_H2_NUM_T] = { 0., 0., 0.,
    0., 0., 0.,
    2.56e-11, 0., 0.,
    4.12e-10, 0., 0.,
    2.37e-9, 0., 0.,
    8.37e-9, 0., 0.,
    2.22e-8, 0., 0.,
    4.88e-8, 0., 0.,
    9.33e-8, 0., 0.,
    1.61e-7, 0., 0.,
    2.55e-7, 0., 0.,
    3.79e-7, 0., 0.,
    5.32e-7, 0., 0.,
    7.13e-7, 0., 0.,
    9.18e-7, 0., 0.,
    1.14e-6, 0., 0.,
    1.37e-6, 0., 0.,
    1.61e-6, 0., 0.,
    1.84e-6, 0., 0.,
    2.05e-6, 0., 0.
  };

  for(i = 0; i < TGCHEM_H2_NUM_V; i++)
    for(j = 0; j < TGCHEM_H2_NUM_J; j++)
      for(k = 0; k < TGCHEM_H2_NUM_V; k++)
        for(l = 0; l < TGCHEM_H2_NUM_T; l++)
          {
            line_idx = (TGCHEM_H2_NUM_J * TGCHEM_H2_NUM_V * TGCHEM_H2_NUM_T) * i + (TGCHEM_H2_NUM_V * TGCHEM_H2_NUM_T) * j + TGCHEM_H2_NUM_T * k + l;

            TGCD.H2SpontCoeff[line_idx] = 0;

            idx = TGCHEM_H2_NUM_T * j + l;

            if(i == 0 && k == 0)
              TGCD.H2SpontCoeff[line_idx] = tmp0[idx];

            if(i == 1 && k == 0)
              TGCD.H2SpontCoeff[line_idx] = tmp1[idx];

            if(i == 1 && k == 1)
              TGCD.H2SpontCoeff[line_idx] = tmp2[idx];

            if(i == 2 && k == 0)
              TGCD.H2SpontCoeff[line_idx] = tmp3[idx];

            if(i == 2 && k == 1)
              TGCD.H2SpontCoeff[line_idx] = tmp4[idx];

            if(i == 2 && k == 2)
              TGCD.H2SpontCoeff[line_idx] = tmp5[idx];

            //mpi_printf("line idx: %d %d %d %d: %g\n", i, j, k, l, TGCD.H2SpontCoeff[line_idx]);
          }
}


void tgchem_init_opac()
{
  double opac[TGCHEM_OPAC_NUM_RHO * TGCHEM_OPAC_NUM_TEMP] = { -4.36, -4.36, -4.36, -4.36, -4.36, -4.36, -4.36, -4.36, -4.36, -4.32, -4.07, -3.34, -2.38, -1.39, -0.39,
    -4.54, -4.54, -4.54, -4.54, -4.54, -4.54, -4.54, -4.53, -4.53, -4.45, -4.02, -3.16, -2.18, -1.18, -0.18,
    -4.71, -4.71, -4.71, -4.71, -4.71, -4.71, -4.71, -4.71, -4.69, -4.54, -3.95, -3.02, -2.03, -1.03, -0.03,
    -4.89, -4.89, -4.89, -4.89, -4.89, -4.89, -4.89, -4.88, -4.85, -4.60, -3.88, -2.92, -1.92, -0.92, 0.08,
    -5.06, -5.06, -5.06, -5.06, -5.06, -5.06, -5.05, -5.05, -4.99, -4.63, -3.81, -2.83, -1.84, -0.84, 0.16,
    -5.39, -5.29, -5.24, -5.22, -5.22, -5.22, -5.22, -5.21, -5.11, -4.65, -3.77, -2.79, -1.79, -0.79, 0.21,
    -6.01, -5.98, -5.90, -5.75, -5.57, -5.45, -5.40, -5.36, -5.23, -4.68, -3.77, -2.78, -1.78, -0.78, 0.22,
    -6.08, -6.08, -6.08, -6.07, -6.05, -5.98, -5.86, -5.68, -5.39, -4.73, -3.79, -2.80, -1.80, -0.80, 0.2,
    -6.13, -6.13, -6.13, -6.13, -6.13, -6.12, -6.10, -6.00, -5.60, -4.80, -3.83, -2.84, -1.84, -0.84, 0.16,
    -6.16, -6.16, -6.16, -6.16, -6.16, -6.16, -6.15, -6.08, -5.68, -4.84, -3.87, -2.87, -1.87, -0.87, 0.13,
    -6.15, -6.15, -6.15, -6.15, -6.15, -6.15, -6.14, -6.07, -5.69, -4.86, -3.88, -2.88, -1.88, -0.88, 0.12,
    -6.05, -6.05, -6.05, -6.05, -6.05, -6.05, -6.05, -5.99, -5.64, -4.82, -3.85, -2.85, -1.85, -0.85, 0.15,
    -5.91, -5.91, -5.91, -5.91, -5.91, -5.91, -5.91, -5.85, -5.52, -4.72, -3.74, -2.75, -1.75, -0.75, 0.25,
    -5.73, -5.57, -5.51, -5.48, -5.48, -5.48, -5.47, -5.44, -5.23, -4.55, -3.60, -2.60, -1.60, -0.60, 0.4,
    -8.02, -7.06, -6.31, -5.82, -5.27, -4.83, -4.65, -4.58, -4.52, -4.23, -3.46, -2.50, -1.50, -0.50, 0.5,
    -11.02, -10.02, -9.02, -8.02, -7.07, -6.32, -5.79, -5.08, -4.33, -3.84, -3.29, -2.41, -1.43, -0.43, 0.57,
    -12.05, -11.51, -10.91, -10.20, -9.34, -8.39, -7.42, -6.53, -5.62, -4.51, -3.44, -2.42, -1.43, -0.43, 0.57,
    -9.69, -9.21, -8.72, -8.22, -7.72, -7.22, -6.71, -6.18, -5.50, -4.59, -3.69, -2.57, -1.50, -0.47, 0.54,
    -6.65, -6.57, -6.39, -6.08, -5.66, -5.19, -4.70, -4.21, -3.71, -3.21, -2.72, -2.23, -1.52, -0.53, 0.49,
    -3.71, -3.52, -3.46, -3.42, -3.36, -3.23, -2.96, -2.58, -2.12, -1.63, -1.14, -0.69, -0.33, 0.12, 0.71,
    -3.49, -2.52, -1.66, -1.18, -1.00, -0.94, -0.89, -0.81, -0.63, -0.31, 0.11, 0.57, 1.01, 1.37, 1.76,
    -3.84, -2.93, -1.94, -0.95, -0.02, 0.62, 0.89, 0.99, 1.04, 1.12, 1.30, 1.61, 2.02, 2.40, 2.74,
    -3.87, -2.91, -2.08, -1.32, -0.38, 0.61, 1.52, 2.13, 2.38, 2.47, 2.54, 2.66, 2.89, 3.23, 3.54,
    -4.53, -3.53, -2.53, -1.54, -0.58, 0.28, 1.18, 2.14, 2.95, 3.39, 3.55, 3.55, 3.55, 3.55, 3.55,
    -5.08, -4.08, -3.08, -2.08, -1.08, -0.08, 0.90, 1.81, 2.73, 3.60, 4.14, 4.14, 4.14, 4.14, 4.14,
    -5.39, -4.50, -3.55, -2.55, -1.55, -0.56, 0.44, 1.44, 2.41, 3.33, 4.18, 4.18, 4.18, 4.18, 4.18,
    -5.75, -4.75, -3.76, -2.81, -1.93, -0.97, 0.02, 1.02, 2.02, 3.01, 3.94, 3.94, 3.94, 3.94, 3.94,
    -6.13, -5.13, -4.13, -3.13, -2.13, -1.15, -0.24, 0.69, 1.68, 2.68, 3.67, 3.67, 3.67, 3.67, 3.67,
    -6.45, -5.45, -4.45, -3.45, -2.45, -1.45, -0.45, 0.53, 1.47, 2.43, 3.42, 3.42, 3.42, 3.42, 3.42
  };

  memcpy(TGCD.OpacTable, opac, sizeof(opac));
}


int compare_emiss(const void *a, const void *b)
{
  if(((struct LineData_struct *) a)->emiss > (((struct LineData_struct *) b)->emiss))
    return -1;

  if(((struct LineData_struct *) a)->emiss < (((struct LineData_struct *) b)->emiss))
    return +1;

  return 0;
}
