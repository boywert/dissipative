/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/healray/healray_finish.c
 * \date        01/2013
 * \author      Thomas Greif
 * \brief       Adaptive ray-tracing 
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 */

#include "../allvars.h"
#include "../proto.h"

void healray_update_cells();
void healray_restore_mesh();
void healray_test_error();
void healray_sanity_and_logging();
void healray_free();


void healray_finish_step()
{
  healray_update_cells();

  healray_finish_rayout();

  healray_test_error();

  healray_sanity_and_logging();

  omp_set_num_threads(NUM_THREADS);

  healray_free();
#ifndef FORCE_EQUAL_TIMESTEPS
  healray_restore_mesh();
#endif
}


void healray_update_cells()
{
  int i, j;

  int tot_num = 0;
  int all_num, all_gas;

  for(i = 0; i < NumGas; i++)
    {
      if(HRD.RayMode)
        {
          if(HRD.EnergyAux[i])
            SphP[i].HeatRate = HRD.Energy[i] / HRD.EnergyAux[i] * HRD.EnergyUnit;

          if(HRD.EnergyAux[i])
            SphP[i].HeatRate = HRD.Energy[i] * HRD.EnergyUnit;

          tot_num += HRD.Energy[i];
        }
      else
        {
          if(HRD.Energy[i])
            SphP[i].HeatRate = HRD.Energy[i] * HRD.EnergyUnit / (HRD.VolToCgs * SphP[i].Volume);

          if(HRD.EnergyFrac[i] >= 0)
            SphP[i].EscFrac = HRD.EnergyFrac[i];

          if(HRD.RayDebugMode && HRD.EnergyFrac[i] >= 0)
            printf("escfrac: %d %d %g\n", ThisTask, i, SphP[i].EscFrac);
        }
    }

  MPI_Allreduce(&tot_num, &all_num, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&NumGas, &all_gas, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  all_num /= all_gas;

  mpi_printf("avg rays = %d %d %d %d\n", NumGas, tot_num, all_gas, all_num);
}


void healray_test_error()
{
  int i, j, k, count, tot_count;
  double dr, r, r_cgs, r2, h, heat_rate, error, tot_error;

  if(HRD.RayNumSources)
    {
      count = 0;
      tot_error = 0.;

      for(i = 0; i < NumGas; i++)
        {
          if(P[i].ID == 0 && P[i].Mass == 0)
            continue;

          heat_rate = 0;

          for(j = 0; j < HRD.RayNumSources; j++)
            {
              r2 = 0;

              for(k = 0; k < 3; k++)
                {
                  dr = P[i].Pos[k] - HRD.SourcePosList[3 * j + k];

                  r2 += dr * dr;
                }

              r = sqrt(r2);

              h = get_cell_radius(i);

              if(!r)
                r = h;

              r_cgs = HRD.LenToCgs * r;

              heat_rate += HRD.SourceEgyList[j] * HRD.EnergyUnit / 4. / M_PI / r_cgs / r_cgs * exp(-HRD.Alpha * r_cgs) * HRD.Alpha;
            }

          count++;

          error = log10(dabs(SphP[i].HeatRate - heat_rate) / heat_rate);
          tot_error += error;
        }

      MPI_Allreduce(&count, &tot_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&tot_error, &HRD.Error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      HRD.Error = pow(10., HRD.Error / tot_count);
    }
}


void healray_sanity_and_logging()
{
  int i;
  double energy_absorb, tot_energy_absorb, tot_energy_esc, tot_energy_err;

  if(!HRD.RayMode)
    {
      energy_absorb = 0.;

      for(i = 0; i < NumGas; i++)
        energy_absorb += HRD.Energy[i];

      MPI_Allreduce(&HRD.TotEnergyEsc, &tot_energy_esc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&energy_absorb, &tot_energy_absorb, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      if(HRD.TotEnergyInit)
        tot_energy_err = dabs(HRD.TotEnergyInit - tot_energy_absorb - tot_energy_esc) / HRD.TotEnergyInit;
      else
        tot_energy_err = 0.;

      if(tot_energy_err > 1e-10)
        {
          mpi_printf("HEALRAY: Unusually large error in total energy: %g. Something is wrong!\n", tot_energy_err);

          endrun();
        }
    }

  MPI_Allreduce(&HRD.LocNumOps, &HRD.NumOps, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&HRD.LocNumAdvances, &HRD.NumAdvances, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&HRD.LocNumComms, &HRD.NumComms, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

  HRD.AdvancesPerRay = (double) HRD.NumAdvances / HRD.TotInitNumRays;
  HRD.CommsPerRay = (double) HRD.NumComms / HRD.TotInitNumRays;
}


void healray_free()
{
  int i;

  if(HRD.RayMode)
    myfree_movable(HRD.PartRayCount);

  myfree_movable(HRD.RayEnergy);
  myfree_movable(HRD.RayEnergyLineInit);
  myfree_movable(HRD.RayEnergyLine);
  myfree_movable(HRD.HRRD);

  if(HRD.RayMode)
    myfree_movable(HRD.EnergyAux);
  else
    myfree_movable(HRD.EnergyFrac);

  myfree_movable(HRD.Energy);

  if(HRD.RayMode)
    myfree_movable(HRD.IdxList);
  else
    myfree_movable(HRD.EnergyInit);

  if(HRD.RayNumSources)
    {
      myfree_movable(HRD.SourceEgyList);
      myfree_movable(HRD.SourcePosList);
      myfree_movable(HRD.SourceIDList);
    }

  if(HRD.RayNumTrace)
    {
#ifdef TGCHEM
      if(!TGCD.ChemMode && HRD.RayMultiFreqMode)
        {
          for(i = HRD.RayNumThreads - 1; i >= 0; i--)
            myfree_movable(HRD.TraceEnergyFracFreq[i]);

          myfree_movable(HRD.TraceEnergyFracFreq);
        }
#endif
      for(i = HRD.RayNumThreads - 1; i >= 0; i--)
        myfree_movable(HRD.TraceEnergyFrac[i]);

      myfree_movable(HRD.TraceEnergyFrac);
      myfree_movable(HRD.TraceList);
    }

  myfree_movable(HRD.NuDFac);
  myfree_movable(HRD.NuDPos);

  myfree_movable(HRD.ThreadPrevNumRays);
  myfree_movable(HRD.ThreadNumRays);
}


void healray_restore_mesh()
{
  int i;

  free_mesh();

  for(i = 0; i < TIMEBINS; i++)
    TimeBinSynchronized[i] = HRD.SaveTimeBinActive[i];

  for(i = 0; i < NumPart; i++)
    P[i].TimeBin = HRD.SaveTimeBin[i];

  CPU_Step[CPU_HEALRAY] += measure_time();

  reconstruct_timebins();

  create_mesh();

  mesh_setup_exchange();

  myfree_movable(HRD.SaveTimeBin);
}
