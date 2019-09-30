/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/tgset.c
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

#include "allvars.h"
#include "proto.h"

void tgset_init(void)
{
  All.NHMax = 0;

  if(All.SnapDensThresh > 0)
    {
      tgset_get_nh_max();

      while(All.NHMax > All.SnapDensThresh)
        All.SnapDensThresh *= pow(10.0, 0.25);
    }


}

void tgset_check_for_snap(void)
{
  if(All.SnapNumFac > 0)
    {
#ifndef ISOTHERM_EQS
      if(All.NHMax > 1e2)
        {
#endif
          if((All.NumCurrentTiStep + 2) % All.SnapNumFac == 0)
            produce_dump();
#ifndef ISOTHERM_EQS
        }
      else
        {
          if(All.Ti_Current >= All.Ti_nextoutput && All.Ti_nextoutput >= 0)
            {
              produce_dump();

              All.Ti_nextoutput = find_next_outputtime(All.Ti_Current + 1);
            }
        }
#endif
    }

  if(All.SnapDensThresh > 0 && All.NHMax > All.SnapDensThresh)
    {
      produce_dump();

      All.SnapDensThresh *= pow(10.0, 0.25);

      if(ThisTask == 0)
        printf("SnapDensThresh = %g\n", All.SnapDensThresh);
    }

  if(All.SnapNumFac == 0 && All.SnapDensThresh == 0)
    if(All.Ti_Current >= All.Ti_nextoutput && All.Ti_nextoutput >= 0)
      {
        produce_dump();

        All.Ti_nextoutput = find_next_outputtime(All.Ti_Current + 1);
      }
}

void tgset_get_nh_max(void)
{
  int idx, i;
  double a3inv, rho, nh;

  struct
  {
    double nh_max;
    int task;
  } local, global;

  if(All.ComovingIntegrationOn)
    a3inv = 1.0 / (All.Time * All.Time * All.Time);
  else
    a3inv = 1;

  local.task = ThisTask;
  local.nh_max = 0;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;


      rho = SphP[i].Density * a3inv * All.HubbleParam * All.HubbleParam * All.UnitDensity_in_cgs;

      nh = HYDROGEN_MASSFRAC * rho / PROTONMASS;

      if(nh > local.nh_max)
        {
          local.nh_max = nh;
          All.NHMaxIdx = i;
        }
    }

  MPI_Allreduce(&local, &global, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

  if(global.nh_max > All.NHMax)
    {
      All.NHMax = global.nh_max;
      All.NHMaxTask = global.task;

      MPI_Bcast(&All.NHMaxIdx, 1, MPI_INT, All.NHMaxTask, MPI_COMM_WORLD);
    }

  mpi_printf("\nMaximum Density = %g\n\n", All.NHMax);

  if(All.NHMax > 1e16)
    terminate("Maximum Density reached!\n");
}

#ifdef REFINEMENT
void tgset_jeans_init(void)
{
  int i;
  double nh_min, nh_max, gamma_eff_min, gamma_eff_max;

  if(JEANS_TABLE_SIZE < 2)
    terminate("JEANS_TABLE_SIZE must be at least 2!\n");

  nh_min = 1e5;
  nh_max = 1e16;

  for(i = 0; i < JEANS_TABLE_SIZE; i++)
    {
      All.JeansLogNHTable[i] = log10(nh_min) + i * (log10(nh_max) - log10(nh_min)) / (JEANS_TABLE_SIZE - 1);
      All.TargetJeansNumberTable[i] = All.JeansNumber;
    }

  //for(i = 0; i < JEANS_TABLE_SIZE; i++)
  //mpi_printf("%g %g\n", All.JeansLogNHTable[i], All.TargetJeansNumberTable[i]);
}

int tgset_jeans_ref(int i)
{
  int index;
  double a, a3inv;
  double rho, cellrad, t_ff, mu, c_s, nh, log_nh, fac;
  double jeans_length, jeans_number, target_jeans_number;

  if(All.ComovingIntegrationOn)
    {
      a = All.Time;
      a3inv = 1 / (a * a * a);
    }
  else
    a = a3inv = 1;

  rho = SphP[i].Density * a3inv * All.HubbleParam * All.HubbleParam * All.UnitDensity_in_cgs;

  t_ff = sqrt(3. * M_PI / 32. / GRAVITY / rho);

#ifdef PRIMCHEM
  mu = (1. + 4. * HE_ABUND) / (1. + HE_ABUND - SphP[i].TracAbund[0]);
#else
  mu = 4. / (1. + 3. * HYDROGEN_MASSFRAC);
#endif

  c_s = sqrt(GAMMA * BOLTZMANN * All.JeansTemp / mu / PROTONMASS);

  jeans_length = c_s * t_ff;

  cellrad = get_cell_radius(i) * a / All.HubbleParam * All.UnitLength_in_cm;

  jeans_number = jeans_length / cellrad;

  nh = HYDROGEN_MASSFRAC * rho / PROTONMASS;

  log_nh = log10(nh);

  index = (int) ((log_nh - All.JeansLogNHTable[0]) / (All.JeansLogNHTable[JEANS_TABLE_SIZE - 1] - All.JeansLogNHTable[0]) * (JEANS_TABLE_SIZE - 1));

  if(index < 0)
    index = 0;

  if(index > JEANS_TABLE_SIZE - 1)
    index = JEANS_TABLE_SIZE - 1;

  target_jeans_number = All.TargetJeansNumberTable[index];

  if(log_nh > All.JeansLogNHTable[0] && log_nh < All.JeansLogNHTable[JEANS_TABLE_SIZE - 1] && index >= 0 && index < JEANS_TABLE_SIZE - 1)
    {
      fac = (log_nh - All.JeansLogNHTable[index]) / (All.JeansLogNHTable[index + 1] - All.JeansLogNHTable[index]);

      target_jeans_number += fac * (All.TargetJeansNumberTable[index + 1] - All.TargetJeansNumberTable[index]);
    }

  if(jeans_number < target_jeans_number)
    return 1;

  /*
     #ifdef ISOTHERM_EQS
     lambda_jeans = All.IsoSoundSpeed * All.UnitVelocity_in_cm_per_s * sqrt(M_PI / GRAVITY / rho);
     #else
     #ifdef PRIMCHEM
     lambda_jeans = sqrt(SphP[i].Gamma * (SphP[i].Gamma - 1) * SphP[i].Utherm * All.UnitEnergy_in_cgs / All.UnitMass_in_g * M_PI / GRAVITY / rho);
     #else
     lambda_jeans = sqrt(GAMMA * GAMMA_MINUS1 * SphP[i].Utherm * All.UnitEnergy_in_cgs / All.UnitMass_in_g * M_PI / GRAVITY / rho);
     #endif
     #endif
   */

  return 0;
}

#endif
