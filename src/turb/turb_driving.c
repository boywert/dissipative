/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/turb/turb_driving.c
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

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "../allvars.h"
#include "../proto.h"

#if (defined (VS_TURB) || defined (AB_TURB)) && defined(STATICNFW)
#define POTBINS 10000
static double pot[POTBINS];
#endif

#ifdef  VS_TURB


#define NVEC 3

double kvec[3][NVEC];
double ampl[3][NVEC];
double phase[3][NVEC];

static int TurbLastSet;

static double TurbDt = 0.0115;  //0.08;

static double TurbAmpl = 0.6741573197587302;    // 15;



void init_turb(void)
{
  kvec[0][0] = 2 * M_PI / All.BoxSize;
  kvec[1][0] = 0;
  kvec[2][0] = 0;

  kvec[0][1] = 0;
  kvec[1][1] = 2 * M_PI / All.BoxSize;
  kvec[2][1] = 0;

  kvec[0][2] = 0;
  kvec[1][2] = 0;
  kvec[2][2] = 2 * M_PI / All.BoxSize;

  TurbLastSet = -1;

  set_turb_ampl();

  init_static_nfw();
}

#endif

#if defined (VS_TURB) || defined (AB_TURB)

void init_static_nfw()
{
#ifdef STATICNFW
  double r, dr, m;
  int i;

  dr = (All.BoxSize / POTBINS);

  for(i = 1, pot[0] = 0; i <= POTBINS; i++)
    {
      r = dr * (i - 0.5);

      m = enclosed_mass(r);

#ifdef NFW_DARKFRACTION
      m *= NFW_DARKFRACTION;
#endif
      pot[i] = pot[i - 1] + All.G * m / (r * r) * dr;
    }
#endif

}

double get_turb_pot(double x, double y, double z)
{
#ifdef STATICNFW
  double dr = (All.BoxSize / POTBINS);

  double dx = x - boxHalf_X;
  double dy = y - boxHalf_Y;
  double dz = z - boxHalf_Z;

  double r = sqrt(dx * dx + dy * dy + dz * dz);

  double u = r / dr;
  int i = u;
  u -= i;

  if(i >= POTBINS)
    terminate("a");

  return pot[i] * (1 - u) + pot[i + 1] * u;
#else
  return 0;
#endif
}

#endif

#ifdef VS_TURB

void set_turb_ampl(void)
{
  int k, i, dim, slot;

  slot = All.Time / TurbDt;

  if(slot > TurbLastSet)
    {
      double dt = All.Time - All.SetLastTime;

      if(dt > 0)
        for(i = 0; i < NumGas; i++)
          {
            SphP[i].DuDt_diss = SphP[i].EgyDiss / P[i].Mass / dt;
            SphP[i].EgyDiss = 0;

            SphP[i].DuDt_drive = SphP[i].EgyDrive / P[i].Mass / dt;
            SphP[i].EgyDrive = 0;
          }

      All.SetLastTime = All.Time;


      gsl_rng_set(random_generator, (long) (slot * 100 + 42));

      for(k = 0; k < NVEC; k++)
        {
          for(dim = 0; dim < 3; dim++)
            {
              phase[dim][k] = 2 * M_PI * gsl_rng_uniform(random_generator);
              ampl[dim][k] = gsl_ran_gaussian(random_generator, 1.0);
            }

          /* now project out component parallel to k-vector to get solenoidal forcing */
          double k2 = kvec[0][k] * kvec[0][k] + kvec[1][k] * kvec[1][k] + kvec[2][k] * kvec[2][k];
          double ak2 = kvec[0][k] * ampl[0][k] + kvec[1][k] * ampl[1][k] + kvec[2][k] * ampl[2][k];

          for(dim = 0; dim < 3; dim++)
            ampl[dim][k] -= ak2 / k2 * kvec[dim][k];
        }

      TurbLastSet = slot;

      MPI_Bcast(&phase[0][0], 3 * NVEC, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&ampl[0][0], 3 * NVEC, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
}


void add_turb_accel(void)
{
  int idx, i, k, dim;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Type == 0)
        {
          for(dim = 0; dim < 3; dim++)
            SphP[i].TurbAccel[dim] = 0;

          for(k = 0; k < NVEC; k++)
            {
              double arg = P[i].Pos[0] * kvec[0][k] + P[i].Pos[1] * kvec[1][k] + P[i].Pos[2] * kvec[2][k];

              for(dim = 0; dim < 3; dim++)
                {
                  double co = cos(arg + phase[dim][k]);
                  double ac = TurbAmpl * ampl[dim][k] * co;

                  SphP[i].TurbAccel[dim] += ac;
                }
            }
        }
    }
}

#endif

#if defined (VS_TURB) || defined (AB_TURB)

void reset_turb_temp(void)
{
  TIMER_START(CPU_TURB_RESET) int idx, i;
  double esum = 0, globsum;

#ifndef DG
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Type == 0)
        {
          double u_target = All.RefEntropy * pow(SphP[i].Density, GAMMA_MINUS1) / GAMMA_MINUS1;

          double du = SphP[i].Utherm - u_target;

          SphP[i].Utherm = u_target;

          SphP[i].Energy =
            P[i].Mass * (SphP[i].Utherm + 0.5 * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]) + get_turb_pot(P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]));

#ifdef MHD
          SphP[i].Energy += 0.5 * (SphP[i].B[0] * SphP[i].B[0] + SphP[i].B[1] * SphP[i].B[1] + SphP[i].B[2] * SphP[i].B[2]) * SphP[i].Volume;
#endif
          SphP[i].Pressure = GAMMA_MINUS1 * SphP[i].Density * SphP[i].Utherm;

          SphP[i].EgyDiss += P[i].Mass * du;

          esum += P[i].Mass * du;
        }
    }
#else
  double w_cell[5];
  int l, j;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      double du = SphP[i].Utherm;

      for(l = 0; l < Nof_base_functions; l++)
        {
          SphP[i].Weights[l][W_E] = 0.;
        }

      for(j = 0; j < Nof_inner_quad_points; j++)        //loop over quadrature points
        {
          calc_state_at_quad_point(SphP[i].Weights, j, w_cell);

          for(l = 0; l < Nof_base_functions; l++)       //loop over base functions
            {
              SphP[i].Weights[l][W_E] +=
                (0.5 * (w_cell[W_PX] * w_cell[W_PX] + w_cell[W_PY] * w_cell[W_PY] + w_cell[W_PZ] * w_cell[W_PZ]) / w_cell[W_RHO] +
                 All.RefEntropy * pow(w_cell[W_RHO], GAMMA) / GAMMA_MINUS1) * GET_Inner_base_values(j, l) * GET_Inner_quad_points_weights(j);


              SphP[i].WeightsDiss[l] +=
                (w_cell[W_E] - 0.5 * (w_cell[W_PX] * w_cell[W_PX] + w_cell[W_PY] * w_cell[W_PY] + w_cell[W_PZ] * w_cell[W_PZ]) / w_cell[W_RHO] -
                 All.RefEntropy * pow(w_cell[W_RHO], GAMMA) / GAMMA_MINUS1) * GET_Inner_base_values(j, l) * GET_Inner_quad_points_weights(j) / DG_PROJ_NORM;
            }
        }

      for(l = 0; l < Nof_base_functions; l++)
        {
          SphP[i].Weights[l][W_E] /= DG_PROJ_NORM;
        }

      SphP[i].Utherm =
        SphP[i].Weights[0][W_E] / SphP[i].Weights[0][W_RHO] - 0.5 * (SphP[i].Weights[0][W_PX] * SphP[i].Weights[0][W_PX] + SphP[i].Weights[0][W_PY] * SphP[i].Weights[0][W_PY] +
                                                                     SphP[i].Weights[0][W_PZ] * SphP[i].Weights[0][W_PZ]) / (SphP[i].Weights[0][W_RHO] * SphP[i].Weights[0][W_RHO]);

      du -= SphP[i].Utherm;

      SphP[i].Energy = SphP[i].Weights[0][W_E] * SphP[i].Volume;

      SphP[i].Pressure = GAMMA_MINUS1 * SphP[i].Density * SphP[i].Utherm;

      SphP[i].EgyDiss += P[i].Mass * du;

      esum += P[i].Mass * du;
    }

  double (*a)[NOF_BASE_FUNCTIONS][5] = (double (*)[NOF_BASE_FUNCTIONS][5]) mymalloc("a", NumGas * Nof_base_functions * 5 * sizeof(double));
  copy_cell_weights_to_array(a);
  minmod_limiter(a);
  copy_array_to_cell_weights(a);
  myfree(a);
#endif

  MPI_Allreduce(&esum, &globsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  All.TurbDissipatedEnergy += globsum;

TIMER_STOP(CPU_TURB_RESET)}


void do_turb_driving_step_first_half(void)
{
  set_turb_ampl();


#ifndef DG
  int idx, i, j;
  integertime ti_step, tstart, tend;
  double dvel[3], dt_gravkick;


  CPU_Step[CPU_MISC] += measure_time();

  add_turb_accel();

  double esum = 0, globsum;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      ti_step = P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0;

      tstart = All.Ti_begstep[P[i].TimeBinHydro];       /* beginning of step */
      tend = tstart + ti_step / 2;      /* midpoint of step */

      if(All.ComovingIntegrationOn)
        dt_gravkick = get_gravkick_factor(tstart, tend);
      else
        dt_gravkick = (tend - tstart) * All.Timebase_interval;

      if(P[i].Type == 0)
        {
          double ekin0 = 0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);

          for(j = 0; j < 3; j++)
            {
              dvel[j] = SphP[i].TurbAccel[j] * dt_gravkick;
              SphP[i].Momentum[j] += P[i].Mass * dvel[j];
              P[i].Vel[j] += dvel[j];
            }

          double ekin1 = 0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);

          SphP[i].Energy += ekin1 - ekin0;

          SphP[i].EgyDrive += ekin1 - ekin0;

          esum += ekin1 - ekin0;
        }
    }

  MPI_Allreduce(&esum, &globsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  All.TurbInjectedEnergy += globsum;

  CPU_Step[CPU_DRIFTS] += measure_time();
#endif
}


void do_turb_driving_step_second_half(void)
{
#ifndef DG
  CPU_Step[CPU_MISC] += measure_time();
  int idx, i, j;
  integertime ti_step, tstart, tend;
  double dvel[3], dt_gravkick;

  double esum = 0, globsum;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      ti_step = P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0;

      tend = All.Ti_begstep[P[i].TimeBinHydro]; /* end of step (Note: All.Ti_begstep[] has already been advanced for the next step at this point)   */
      tstart = tend - ti_step / 2;      /* midpoint of step */

      if(All.ComovingIntegrationOn)
        dt_gravkick = get_gravkick_factor(tstart, tend);
      else
        dt_gravkick = (tend - tstart) * All.Timebase_interval;

      if(P[i].Type == 0)
        {
          double ekin0 = 0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);

          for(j = 0; j < 3; j++)
            {
              dvel[j] = SphP[i].TurbAccel[j] * dt_gravkick;
              SphP[i].Momentum[j] += P[i].Mass * dvel[j];
              P[i].Vel[j] += dvel[j];
            }

          double ekin1 = 0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);

          SphP[i].Energy += ekin1 - ekin0;

          SphP[i].EgyDrive += ekin1 - ekin0;

          esum += ekin1 - ekin0;
        }
    }
  MPI_Allreduce(&esum, &globsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  All.TurbInjectedEnergy += globsum;

  CPU_Step[CPU_DRIFTS] += measure_time();
#endif
}




void log_turb_temp(void)
{
  int i;
  double dudt_drive = 0;
  double dudt_diss = 0;
  double mass = 0;
  double ekin = 0;
  double etot = 0;

  double vmean[] = { 0., 0., 0. };
  double vmean_global[] = { 0., 0., 0. };
  int count;
  long long sum;

  for(i = 0; i < NumGas; i++)
    {
      if(P[i].Mass == 0 && P[i].ID == 0)
        continue;

      vmean[0] += P[i].Vel[0];
      vmean[1] += P[i].Vel[1];
      vmean[2] += P[i].Vel[2];
      count++;
    }

  MPI_Allreduce(vmean, vmean_global, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  sumup_large_ints(1, &count, &sum);

  vmean_global[0] /= sum;
  vmean_global[1] /= sum;
  vmean_global[2] /= sum;

  for(i = 0; i < NumGas; i++)
    {
      if(P[i].Mass == 0 && P[i].ID == 0)
        continue;

      dudt_drive += P[i].Mass * SphP[i].DuDt_drive;
      dudt_diss += P[i].Mass * SphP[i].DuDt_diss;
      ekin +=
        0.5 * P[i].Mass * ((P[i].Vel[0] - vmean_global[0]) * (P[i].Vel[0] - vmean_global[0]) + (P[i].Vel[1] - vmean_global[1]) * (P[i].Vel[1] - vmean_global[1]) +
                           (P[i].Vel[2] - vmean_global[2]) * (P[i].Vel[2] - vmean_global[2]));

      etot += P[i].Mass * SphP[i].Utherm + P[i].Mass * get_turb_pot(P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);

      mass += P[i].Mass;
    }

  etot += ekin;

  double glob_mass, glob_dudt_drive, glob_dudt_diss, glob_ekin, glob_etot;

  MPI_Allreduce(&mass, &glob_mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&dudt_drive, &glob_dudt_drive, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&dudt_diss, &glob_dudt_diss, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&ekin, &glob_ekin, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&etot, &glob_etot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  double mach = sqrt(2 * (glob_ekin / glob_mass)) / All.IsoSoundSpeed;

  if(ThisTask == 0)
    fprintf(FdTurb, "%g %g %g %g %g %g %g\n",
            All.Time, mach, glob_etot / glob_mass, glob_dudt_drive / glob_mass, glob_dudt_diss / glob_mass, All.TurbInjectedEnergy / glob_mass, All.TurbDissipatedEnergy / glob_mass);
}



#endif
