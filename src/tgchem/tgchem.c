/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/tgchem/tgchem.c
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

#include "../allvars.h"
#include "../proto.h"

void tgchem_prepare_step(int i);
void tgchem_finish_step(int i);
void tgchem_compute_step_constants();
void tgchem_compute_cell_constants(int i);


void tgchem()
{
  int idx, i, j;
  int tot_num_eq, max_num_rate_calls, max_num_sub_steps;
  long tot_num_cells_done, tot_num_rate_calls, tot_num_sub_steps;
  long avg_num_rate_calls, avg_num_sub_steps;
  double dt, u_old, frac_cells_eq, imba, cells_per_s;
  double t0, t1, tdiff, tdiff_eq, tdiff_avg, tdiff_max, frac_tdiff_eq;

  CPU_Step[CPU_MISC] += measure_time();

  if(TGCD.ChemMode >= 0)
    {
      t0 = second();

      mpi_printf("TGCHEM: Doing chemistry step...\n");

      TGCD.DebugFlag = -1;

      tgchem_compute_step_constants();

      tgchem_init_cvode();

      for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
        {
          i = TimeBinsHydro.ActiveParticleList[idx];
          if(i < 0)
            continue;

          dt = (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;
          dt *= All.UnitTime_in_s / All.cf_hubble_a / All.HubbleParam;
#ifndef FORCE_EQUAL_TIMESTEPS
          if(TGCD.ChemRndMode)
            dt *= pow(2., 2. * get_random_number() - 1.);
#endif
          tgchem_prepare_step(i);

          tgchem_step(dt);

          tgchem_finish_step(i);
        }

      tgchem_finish_cvode();

      t1 = second();

      tdiff = timediff(t0, t1);

      MPI_Allreduce(&TGCD.NumEq, &tot_num_eq, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&TGCD.TotNumRateCalls, &tot_num_rate_calls, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&TGCD.TotNumSubSteps, &tot_num_sub_steps, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&TGCD.NumCellsDone, &tot_num_cells_done, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&TGCD.MaxNumRateCalls, &max_num_rate_calls, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(&TGCD.MaxNumSubSteps, &max_num_sub_steps, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

      if(tot_num_cells_done)
        {
          frac_cells_eq = (double) tot_num_eq / tot_num_cells_done;
          avg_num_rate_calls = tot_num_rate_calls / tot_num_cells_done;
          avg_num_sub_steps = tot_num_sub_steps / tot_num_cells_done;
        }
      else
        frac_cells_eq = avg_num_rate_calls = avg_num_sub_steps = 0;

      MPI_Allreduce(&TGCD.DtEq, &tdiff_eq, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&tdiff, &tdiff_avg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&tdiff, &tdiff_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

      cells_per_s = tot_num_cells_done / tdiff_avg;

      frac_tdiff_eq = tdiff_eq / tdiff_avg;

      tdiff_avg /= NTask;

      imba = (tdiff_max - tdiff_avg) / tdiff_avg;

      mpi_printf
        ("TGCHEM: Done! Rate calls per cell: %d (%d), Eq fraction: %g (%g), Imbalance: %g, Cells per task per second: %g, Took %g seconds\n\n",
         avg_num_rate_calls, max_num_rate_calls, frac_cells_eq, frac_tdiff_eq, imba, cells_per_s, tdiff_max);
#ifdef HEALRAY
      if(HRD.SourceFlag >= 0 && !HRD.RayTimeFac)
        {
          produce_dump();

          endrun();
        }
#endif
    }

  CPU_Step[CPU_TGCHEM] += measure_time();
}


void tgchem_compute_step_constants()
{
  set_cosmo_factors_for_current_time();

  TGCD.NumEq = 0;
  TGCD.NumNeq = 0;
  TGCD.TotNumRateCalls = 0;
  TGCD.TotNumSubSteps = 0;
  TGCD.NumCellsDone = 0;

  TGCD.DtEq = 0;
  TGCD.DtNEq = 0;

  TGCD.RedShift = All.cf_redshift;
  TGCD.TempCMB = TEMP_CMB * (1. + TGCD.RedShift);
  TGCD.ConvertEnergy = pow(All.UnitLength_in_cm, 3) / All.UnitEnergy_in_cgs;
}


void tgchem_prepare_step(int i)
{
  int j;

  tgchem_compute_cell_constants(i);

  for(j = 0; j < TGCHEM_NUM_ABUNDANCES; j++)
    NV_Ith_S(TGCD.Species, j) = SphP[i].Abund[j];

  NV_Ith_S(TGCD.Species, TGCHEM_NUM_ABUNDANCES) = SphP[i].Utherm * TGCD.IntDensity / TGCD.ConvertEnergy;

  for(j = 0; j < TGCHEM_NUM_SPECIES; j++)
    NV_Ith_S(TGCD.SpeciesStepSave, j) = NV_Ith_S(TGCD.Species, j);

  tgchem_compute_vars(0, TGCD.Species);

  tgchem_compute_aux();
}


void tgchem_compute_cell_constants(int i)
{
  TGCD.NumRateCalls = 0;
  TGCD.NumSubSteps = 0;

  TGCD.Task = ThisTask;
  TGCD.Index = i;
  TGCD.ID = P[i].ID;

  TGCD.IntDensity = SphP[i].Density * All.cf_a3inv * All.HubbleParam * All.HubbleParam;
  TGCD.Density = All.UnitDensity_in_cgs * TGCD.IntDensity;
  TGCD.NH = HYDROGEN_MASSFRAC * TGCD.Density / PROTONMASS;
  TGCD.TransPower = 1. / (1. + pow(TGCD.NH / TGCD.ChemHIITrans, 2));
  TGCD.DivVel = SphP[i].DivVel * (1 + TGCD.RedShift) * All.HubbleParam / All.UnitTime_in_s;
  TGCD.HydroHeatRate = SphP[i].HydroHeatRate * TGCD.IntDensity / TGCD.ConvertEnergy;
  TGCD.HydroHeatRate *= All.HubbleParam / All.UnitTime_in_s;
}


void tgchem_finish_step(int i)
{
  int j;
  double Utherm_new, du;

  if(tgchem_check_species())
    {
      printf("Something wrong with species! Mode 3: %d %d %d %g %g %g\n", TGCD.Task, TGCD.Index, TGCD.ID, NV_Ith_S(TGCD.Species, 0), NV_Ith_S(TGCD.Species, 1), NV_Ith_S(TGCD.Species, 2));

      terminate("");
    }

  tgchem_rates(0, TGCD.Species, TGCD.DSpecies, 0);

  SphP[i].Gamma = TGCD.Gamma;

  for(j = 0; j < TGCHEM_NUM_ABUNDANCES; j++)
    SphP[i].Abund[j] = NV_Ith_S(TGCD.Species, j);

  SphP[i].Abund[0] = TGCD.AbHM;

  if(!TGCD.HRHeatRate)
    SphP[i].EscFrac = TGCD.H2EscFrac;

  Utherm_new = dmax(TGCD.ConvertEnergy * NV_Ith_S(TGCD.Species, TGCHEM_NUM_ABUNDANCES) / TGCD.IntDensity, All.MinEgySpec);

  du = Utherm_new - SphP[i].Utherm;

  SphP[i].Utherm += du;

  SphP[i].Energy += du * P[i].Mass * pow(All.cf_atime, 2);
#ifdef USE_ENTROPY_FOR_COLD_FLOWS
  SphP[i].A = (SphP[i].Gamma - 1) * SphP[i].Utherm / pow(SphP[i].Density * All.cf_a3inv, GAMMA_MINUS1);
  SphP[i].Entropy = log(SphP[i].A) * P[i].Mass;
#endif
  set_pressure_of_cell(i);

  TGCD.MaxNumRateCalls = imax(TGCD.NumRateCalls, TGCD.MaxNumRateCalls);
  TGCD.MaxNumSubSteps = imax(TGCD.NumSubSteps, TGCD.MaxNumSubSteps);

  TGCD.TotNumRateCalls += TGCD.NumRateCalls;
  TGCD.TotNumSubSteps += TGCD.NumSubSteps;

  TGCD.NumCellsDone++;
}


double tgchem_get_timestep(int i, double dt)
{
  int j;
  double u, u_old, u_new, du;
  double dt_chem;

  return dt;

  /*
     dt *= All.UnitTime_in_s / All.HubbleParam;

     tgchem_prepare_step(i);

     tgchem_rates(0, TGCD.Species, TGCD.DSpecies, 0);

     dt_chem = 0.05 * NV_Ith_S(TGCD.Species, TGCHEM_NUM_ABUNDANCES) / dabs(NV_Ith_S(TGCD.DSpecies, TGCHEM_NUM_ABUNDANCES));

     dt = dmin(dt, dt_chem);
   */

  /*
     for(j = 0; j < TGCHEM_NUM_ABUNDANCES; j++)
     NV_Ith_S(TGCD.SpeciesStepSave, j) = NV_Ith_S(TGCD.Species, j);

     u_old = NV_Ith_S(TGCD.Species, TGCHEM_NUM_ABUNDANCES);

     if(TGCD.HydroHeatRate)
     du = 2. * TGCD.HydroHeatRate * dt;
     else
     du = All.CourantFac * u_old;

     while(1)
     {
     u = dmax(u_old + du, u_old / 2.);

     for(j = 0; j < TGCHEM_NUM_ABUNDANCES; j++)
     NV_Ith_S(TGCD.Species, j) = NV_Ith_S(TGCD.SpeciesStepSave, j);

     NV_Ith_S(TGCD.Species, TGCHEM_NUM_ABUNDANCES) = u;

     tgchem_compute_vars(0, TGCD.Species);

     tgchem_compute_aux();

     tgchem_step(dt);

     u_new = NV_Ith_S(TGCD.Species, TGCHEM_NUM_ABUNDANCES);

     if(dabs(u_new - u) / u > 0.05)
     {
     dt /= 2.;
     du /= 2.;
     }
     else
     break;
     }
   */
  dt *= All.HubbleParam / All.UnitTime_in_s;

  return dt;
}
