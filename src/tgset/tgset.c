/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/tgset/tgset.c
 * \date        01/2013
 * \author      Thomas Greif
 * \brief       Special settings for primordial runs
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 */

#include "../allvars.h"
#include "../proto.h"


void tgset_begrun()
{
  int i;
  double nh_low, nh_high, gamma_eff_min, gamma_eff_max;
  double jeans_number_low, jeans_number_high;

  MPI_Bcast(&TGD, sizeof(struct TGD_struct), MPI_BYTE, 0, MPI_COMM_WORLD);

  TGD.NHMax = 0.;

  if(TGSET_JEANS_TABLE_SIZE < 2)
    terminate("TGSET_JEANS_TABLE_SIZE must be at least 2!\n");

  if(TGD.JeansDensLow <= 0.)
    TGD.JeansDensLow = 1.;

  if(TGD.JeansDensHigh <= 0.)
    TGD.JeansDensHigh = 1.;

  if(TGD.JeansNumberLow <= 0.)
    TGD.JeansNumberLow = 8.;

  if(TGD.JeansNumberHigh <= 0.)
    TGD.JeansNumberHigh = 8.;

  for(i = 0; i < TGSET_JEANS_TABLE_SIZE; i++)
    {
      TGD.JeansLogNHTable[i] = log10(TGD.JeansDensLow) + i * (log10(TGD.JeansDensHigh) - log10(TGD.JeansDensLow)) / (TGSET_JEANS_TABLE_SIZE - 1);

      TGD.TargetJeansNumberTable[i] = pow(10., log10(TGD.JeansNumberLow) + i * (log10(TGD.JeansNumberHigh) - log10(TGD.JeansNumberLow)) / (TGSET_JEANS_TABLE_SIZE - 1));
    }
}


void tgset_dens_init(int i)
{
  if(TGD.NHInit)
    {
      P[i].Mass = SphP[i].Volume * TGD.NHInit * PROTONMASS / HYDROGEN_MASSFRAC / All.UnitDensity_in_cgs;

      if(All.ComovingIntegrationOn)
        P[i].Mass *= pow(All.Time, 3) / pow(All.HubbleParam, 2);
    }
}


void tgset_check_for_snap()
{
  double dt;

  if(All.TimeBetSnapshot)
    {
      if(All.Ti_Current >= All.Ti_nextoutput && All.Ti_nextoutput >= 0)
        {
          produce_dump();

          All.Ti_nextoutput = find_next_outputtime(All.Ti_Current + 1);
        }
    }
  else
    {
      if(All.NumCurrentTiStep > 0 && (All.NumCurrentTiStep) % TGD.SnapNumFac == 0)
        produce_dump();
    }
}


void tgset_get_nh_max()
{
  int idx, i;
  double rho, nh;

  struct
  {
    double nh_max;
    int task;
  } local, global;

  set_cosmo_factors_for_current_time();

  local.task = ThisTask;
  local.nh_max = 0;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      rho = SphP[i].Density * All.cf_a3inv * All.HubbleParam * All.HubbleParam * All.UnitDensity_in_cgs;

      nh = HYDROGEN_MASSFRAC * rho / PROTONMASS;

      if(nh > local.nh_max)
        {
          local.nh_max = nh;
          TGD.NHMaxIdx = i;
        }
    }

  MPI_Allreduce(&local, &global, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

  if(global.nh_max > TGD.NHMax)
    {
      TGD.NHMax = global.nh_max;
      TGD.NHMaxTask = global.task;

      MPI_Bcast(&TGD.NHMaxIdx, 1, MPI_INT, TGD.NHMaxTask, MPI_COMM_WORLD);
    }

  mpi_printf("\nTGSET: Maximum Density = %g\n\n", TGD.NHMax);

  if(TGD.NHTerm && TGD.NHMax > TGD.NHTerm)
    {
      produce_dump();

      mpi_printf("TGSET: Termination Density reached!\n");

      endrun();
    }
}


int tgset_jeans_ref(int mode, int i)
{
  int index;
  double rho, cellrad, t_ff, gamma, temp, mu, csnd, nh, log_nh, fac;
  double jeans_length, jeans_number, target_jeans_number;

  set_cosmo_factors_for_current_time();

  rho = SphP[i].Density * All.cf_a3inv * All.HubbleParam * All.HubbleParam * All.UnitDensity_in_cgs;

  nh = HYDROGEN_MASSFRAC * rho / PROTONMASS;

  t_ff = sqrt(3. * M_PI / 32. / GRAVITY / rho);

#ifdef TGCHEM
  TGCD.AbH2 = SphP[i].Abund[1];
  TGCD.AbHII = SphP[i].Abund[2];

  mu = (1. + 4. * HE_ABUND) / (1. + HE_ABUND - TGCD.AbH2 + TGCD.AbHII);

  gamma = SphP[i].Gamma;
#else
  if(dabs(TGD.JeansTemp) > 1e4)
    mu = 4. / (8. - 5. * (1. - HYDROGEN_MASSFRAC));
  else
    mu = 4. / (1. + 3. * HYDROGEN_MASSFRAC);

  gamma = GAMMA;
#endif

#ifdef ISOTHERM_EQS
  csnd = All.IsoSoundSpeed * All.UnitVelocity_in_cm_per_s;
#else
  temp = mu * (gamma - 1.) * PROTONMASS / BOLTZMANN * SphP[i].Utherm * All.UnitEnergy_in_cgs / All.UnitMass_in_g;

  if(TGD.JeansTemp < 0 && temp < dabs(TGD.JeansTemp) || TGD.JeansTemp > 0)
    temp = dabs(TGD.JeansTemp);

  csnd = sqrt(gamma * BOLTZMANN * temp / mu / PROTONMASS);
#endif

  jeans_length = csnd * t_ff;

  cellrad = get_cell_radius(i) * All.cf_atime / All.HubbleParam * All.UnitLength_in_cm;

  jeans_number = jeans_length / cellrad;

  log_nh = log10(nh);

  index = (int) ((log_nh - TGD.JeansLogNHTable[0]) / (TGD.JeansLogNHTable[TGSET_JEANS_TABLE_SIZE - 1] - TGD.JeansLogNHTable[0]) * (TGSET_JEANS_TABLE_SIZE - 1));

  index = imin(imax(index, 0), TGSET_JEANS_TABLE_SIZE - 1);

  target_jeans_number = TGD.TargetJeansNumberTable[index];

  if(log_nh > TGD.JeansLogNHTable[0] && log_nh < TGD.JeansLogNHTable[TGSET_JEANS_TABLE_SIZE - 1] && index >= 0 && index < TGSET_JEANS_TABLE_SIZE - 1)
    {
      fac = (log_nh - TGD.JeansLogNHTable[index]) / (TGD.JeansLogNHTable[index + 1] - TGD.JeansLogNHTable[index]);

      target_jeans_number += fac * (TGD.TargetJeansNumberTable[index + 1] - TGD.TargetJeansNumberTable[index]);
    }

  if(mode == 1)
    {
      if(jeans_number > 1.5 * target_jeans_number)
        return 1;
    }
  else
    {
      if(jeans_number < target_jeans_number)
        return 1;
    }

  return 0;
}


void tgset_get_image_limits(char **argv, int *xaxis, int *yaxis, int *zaxis, double *xmin, double *xmax, double *ymin, double *ymax, double *zmin, double *zmax, int *weight_flag)
{
  int center_flag, width_flag;
  double width, center[3];

  center_flag = atoi(argv[9]);
  *weight_flag = atoi(argv[10]);
  width_flag = atoi(argv[11]);
  width = atof(argv[12]);

  set_cosmo_factors_for_current_time();

  if(center_flag == 0)
    {
      tgset_get_nh_max();

      if(ThisTask == TGD.NHMaxTask)
        {
          center[0] = P[TGD.NHMaxIdx].Pos[0];
          center[1] = P[TGD.NHMaxIdx].Pos[1];
          center[2] = P[TGD.NHMaxIdx].Pos[2];
        }

      MPI_Bcast(center, 3, MPI_DOUBLE, TGD.NHMaxTask, MPI_COMM_WORLD);
    }
  else if(center_flag == 1)
    {
      center[0] = atof(argv[13]);
      center[1] = atof(argv[14]);
      center[2] = atof(argv[15]);
    }
  else if(center_flag == 2)
    {
      center[0] = All.BoxSize / 2;
      center[1] = All.BoxSize / 2;
      center[2] = All.BoxSize / 2;
    }
  else if(center_flag == 3)
    {
      int i, j;
      double dens_sum, dens_sum_all, denspos_sum[3], denspos_sum_all[3];

      dens_sum = 0;

      for(j = 0; j < 3; j++)
        denspos_sum[j] = 0;

      for(i = 0; i < NumGas; i++)
        {
          for(j = 0; j < 3; j++)
            denspos_sum[j] += sqrt(SphP[i].Density) * P[i].Pos[j];

          dens_sum += sqrt(SphP[i].Density);
        }

      MPI_Allreduce(&denspos_sum, &denspos_sum_all, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&dens_sum, &dens_sum_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      for(j = 0; j < 3; j++)
        center[j] = denspos_sum_all[j] / dens_sum_all;
    }
#ifdef SINKS
  else if(center_flag == 4)
    {
      int i, j;
      double mass_sum = 0;

      center[0] = center[1] = center[2] = 0;

      if(SKD.TotNumSinks)
        {
          for(i = 0; i < SKD.TotNumSinks; i++)
            {
              for(j = 0; j < 3; j++)
                center[j] += SKD.SINK[i].Mass * SKD.SINK[i].Pos[j];

              mass_sum += SKD.SINK[i].Mass;
            }

          for(j = 0; j < 3; j++)
            center[j] /= mass_sum;
        }
      else
        terminate("Using sinks for centering but found no sinks!\n");
    }
#endif
  else
    terminate("Unknown center flag!\n");

  if(width_flag == 1)
    width *= All.HubbleParam / All.cf_atime * PARSEC / All.UnitLength_in_cm;
  else if(width_flag == 2)
    width *= ASTRONOMICAL_UNIT / All.UnitLength_in_cm;
  else
    width *= All.HubbleParam;

  *xmin = center[*xaxis] - width / 2;
  *xmax = center[*xaxis] + width / 2;
  *ymin = center[*yaxis] - width / 2;
  *ymax = center[*yaxis] + width / 2;

  if(RestartFlag == 4)
    *zmin = center[*zaxis];
  else
    {
      *zmin = center[*zaxis] - width / 2;
      *zmax = center[*zaxis] + width / 2;
    }

  mpi_printf("center = %g %g %g, width = %g\n", center[0], center[1], center[2], width);
}
