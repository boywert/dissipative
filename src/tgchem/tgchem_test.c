/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/tgchem/tgchem_test.c
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

#include "tgchem_test.h"

int ThisTask;
int NTask;

double WallClockTime;

struct TGCD_struct TGCD;


void tgchem_test_pars()
{
  TGCD.ChemMode = 0;
  TGCD.ChemH2Mode = 0;
  TGCD.ChemH2Thresh = 1e15;
  TGCD.ChemHIIThresh = 1e21;
  TGCD.ChemInitAbH2 = 6.6e-7;
  TGCD.ChemInitAbHII = 2.6e-4;
  TGCD.ChemJ21 = 1e5;
}


void tgchem_test_init()
{
  char buf[MAX_STRING_LEN];

  TGCD.CollapseFac = 1;
  TGCD.TestStepSize = 1e-1;
  TGCD.RedShift = 20;
  TGCD.NHStart = 1e-3;
  TGCD.NHEnd = 1e30;
  TGCD.Temp = 2e2;

  TGCD.NH = TGCD.NHStart;
  TGCD.TempCMB = TEMP_CMB * (1 + TGCD.RedShift);

  NV_Ith_S(TGCD.Species, 0) = 0;
  NV_Ith_S(TGCD.Species, 1) = TGCD.ChemInitAbH2;
  NV_Ith_S(TGCD.Species, 2) = TGCD.ChemInitAbHII;

  tgchem_compute_vars(1, TGCD.Species);

  TGCD.ID = 0;
  TGCD.Task = 0;
  TGCD.Index = 0;
  TGCD.DivVel = 0.;
  TGCD.HRHeatRate = 0.;
}


int main(int argc, char **argv)
{
  char buf[MAX_STRING_LEN];
  int i, j, flag_break, num_abundances, num_rates, num_cooling;
  double dt, t_ff, dt_rho, drho_dt, dt_pdv, rho_max, dummy;
  FILE *file;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);

  WallClockTime = second();

  tgchem_test_pars();

  tgchem_begrun();

  tgchem_init_cvode();

  tgchem_test_init();

  rho_max = PROTONMASS * TGCD.NHEnd / HYDROGEN_MASSFRAC;

  sprintf(buf, "../../sator/data/tgchem.dat");

  if(!(file = fopen(buf, "w")))
    terminate("Could not open file!\n");

  num_abundances = TGCHEM_NUM_ABUNDANCES;
  num_rates = TGCHEM_NUM_RATES;
  num_cooling = TGCHEM_NUM_COOLING;

  fseek(file, 4, SEEK_SET);
  fwrite(&num_abundances, sizeof(int), 1, file);
  fwrite(&num_rates, sizeof(int), 1, file);
  fwrite(&num_cooling, sizeof(int), 1, file);
  fwrite(&TGCD.NHStart, sizeof(double), 1, file);
  fwrite(&TGCD.NHEnd, sizeof(double), 1, file);

  for(i = 0; i < TGCHEM_NUM_ABUNDANCES; i++)
    {
      dummy = NV_Ith_S(TGCD.SpeciesTol, i);
      fwrite(&dummy, sizeof(double), 1, file);
    }

  flag_break = 0;

  TGCD.DebugFlag = -1;

  while(!flag_break)
    {
      TGCD.NumRateCalls = 0;

      TGCD.Density = PROTONMASS * TGCD.NH / HYDROGEN_MASSFRAC;

      dt = t_ff = sqrt(3. * M_PI / 32. / GRAVITY / TGCD.Density);

      drho_dt = TGCD.CollapseFac * TGCD.Density / t_ff;

      dt = dmin(TGCD.TestStepSize * TGCD.Density / drho_dt, dt);

      TGCD.HydroHeatRate = 4. / 3. * TGCD.CollapseFac * NV_Ith_S(TGCD.Species, TGCHEM_NUM_ABUNDANCES) / t_ff;

      dt = dmin(TGCD.TestStepSize * NV_Ith_S(TGCD.Species, TGCHEM_NUM_ABUNDANCES) / TGCD.HydroHeatRate, dt);

      if(TGCD.Density + drho_dt * dt > rho_max)
        {
          dt = (rho_max - TGCD.Density) / drho_dt;

          flag_break = 1;
        }

      TGCD.Density += drho_dt * dt;
      TGCD.NH = HYDROGEN_MASSFRAC * TGCD.Density / PROTONMASS;
      TGCD.TransPower = 1. / (1. + pow(TGCD.NH / TGCD.ChemHIITrans, 2));

      NV_Ith_S(TGCD.Species, TGCHEM_NUM_ABUNDANCES) += TGCD.HydroHeatRate * dt;

      tgchem_compute_vars(0, TGCD.Species);

      tgchem_compute_aux();

      tgchem_step(dt);

      tgchem_rates(0, TGCD.Species, TGCD.DSpecies, 0);

      NV_Ith_S(TGCD.Species, 0) = TGCD.AbHM;

      //if(TGCD.Index == TGCD.DebugFlag - 1 || TGCD.Index == TGCD.DebugFlag)
      //if(TGCD.NH > 5e19 && TGCD.NH < 2e20)
      if(0)
        printf("iter = %d, nh = %g, temp = %g, gamma = %g, abhm = %g, abh2 = %g, abhp = %g, NIter = %d\n", TGCD.Index, TGCD.NH, TGCD.Temp, TGCD.Gamma,
               TGCD.AbHM, TGCD.AbH2, TGCD.AbHII, TGCD.NumRateCalls);

      fwrite(&TGCD.NH, sizeof(double), 1, file);
      fwrite(&TGCD.Temp, sizeof(double), 1, file);
      fwrite(&TGCD.Gamma, sizeof(double), 1, file);

      for(i = 0; i < TGCHEM_NUM_ABUNDANCES; i++)
        {
          dummy = NV_Ith_S(TGCD.Species, i);
          fwrite(&dummy, sizeof(double), 1, file);
        }

      for(i = 0; i < TGCHEM_NUM_ABUNDANCES; i++)
        {
          for(j = 0; j < TGCHEM_NUM_RATES; j++)
            {
              dummy = TGCD.CrRate[i][j] / TGCD.NH;
              fwrite(&dummy, sizeof(double), 1, file);
            }

          for(j = 0; j < TGCHEM_NUM_RATES; j++)
            {
              dummy = TGCD.DsRate[i][j] / TGCD.NH;
              fwrite(&dummy, sizeof(double), 1, file);
            }
        }

      for(i = 0; i < TGCHEM_NUM_ABUNDANCES; i++)
        for(j = 0; j < TGCHEM_NUM_RATES; j++)
          {
            dummy = TGCD.ChemicalRate[i][j] / pow(TGCD.NH, 2);
            fwrite(&dummy, sizeof(double), 1, file);
          }

      dummy = TGCD.HydroHeatRate / pow(TGCD.NH, 2);
      fwrite(&dummy, sizeof(double), 1, file);

      for(i = 0; i < TGCHEM_NUM_COOLING; i++)
        {
          dummy = TGCD.CoolingRate[i] / pow(TGCD.NH, 2);
          fwrite(&dummy, sizeof(double), 1, file);
        }

      if(TGCD.Index == TGCD.DebugFlag)
        terminate("");

      TGCD.Index++;
    }

  fseek(file, 0, SEEK_SET);
  fwrite(&TGCD.Index, sizeof(int), 1, file);

  fclose(file);

  tgchem_finish_cvode();

  printf("%d iterations, took %g seconds. Done!\n", TGCD.Index, measure_time());

  MPI_Finalize();
}
