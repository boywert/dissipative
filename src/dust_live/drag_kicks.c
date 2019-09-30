/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dust_live/drag_kicks.c
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
#include <gsl/gsl_math.h>

#include "../allvars.h"
#include "../proto.h"

#ifdef DUST_LIVE

#ifdef DL_GRAIN_BINS
/* For particle index p and bin index i, calculate a quantity (in units of
 * um^2) at grain size a used in determining the effective stopping time-scale.
 * */
double t_s_eff_helper(int p, int i, double a)
{
  double fac1 = DTP(p).NumGrains[i]*a*a*a/3.0 / GSD.Widths[i];
  double fac2 = 0.0;
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
  fac2 = DTP(p).BinSlopes[i] * (pow(a, 4)/4.0 - GSD.Midpoints[i]*a*a*a/3.0);
#endif
  return fac1 + fac2;
}
#endif

void compute_drag_acceleration(void)
{
  for(int i = 0; i < Ndust; i++)
    {
      int p = DustParticle[i].index;
      if(P[p].Mass == 0.0)
        continue;

      /* Calculate quantities and stopping time-scale in physical units.
       * Will convert to comoving when drag acceleration is applied. */
#ifdef DL_GRAIN_BINS
      /* If we're evolving grain size distributions, we need to calculate an
       * effective stopping time-scale. */
      double rho_grain = All.GrainDensity / (All.HubbleParam * All.HubbleParam * All.UnitDensity_in_cgs);
      double a2 = 0.0;
      for(int k = 0; k < DL_GRAIN_BINS; k++)
        {
          a2 += (t_s_eff_helper(p, k, GSD.Edges[k+1]) - t_s_eff_helper(p, k, GSD.Edges[k]));
        }
      if(a2 <= 0.0)
        {
          terminate("DUST_LIVE: Nonpositive effective stopping time-scale, a2=%g, p=%d, mass=%g!", a2, p, P[p].Mass);
        }

      double a_grain = 3.0 * P[p].Mass / (4.0 * M_PI * GSD.InternalDensity * a2); /* um */
      if(a_grain > All.MaxGrainSize)
        a_grain = All.MaxGrainSize;
      if(a_grain < All.MinGrainSize)
        a_grain = All.MinGrainSize;
      a_grain *= 1.0e-4 * All.HubbleParam / All.UnitLength_in_cm; /* internal length, since 1 um = 1.0e-4 cm */
#else
      /* If we're not evolving grain size distributions, just fix the grain
       * density and radius. */
      double rho_grain = 2.4 / (All.HubbleParam * All.HubbleParam * All.UnitDensity_in_cgs);
      double a_grain = 1.0e-5 * All.HubbleParam / All.UnitLength_in_cm;
#endif
      double rho_gas = DTP(p).LocalGasDensity * All.cf_a3inv;
      double c_s = DTP(p).LocalSoundSpeed;
      /* Hopkins+ 2016, Equation 2 */
      double t_s = sqrt(M_PI * GAMMA) / sqrt(8.0) * (rho_grain * a_grain) / (rho_gas * c_s);
#ifdef DL_STOPPING_TIME_CORRECTION
      double v2 = 0.0;
      for(int j = 0; j < 3; j++)
        {
          double dvel = (DTP(p).LocalGasVelocity[j] - P[p].Vel[j]) / All.cf_atime;
          v2 += (dvel * dvel);
        }
      double corr = pow(1.0 + 9.0*M_PI/128.0 * v2 / (c_s * c_s), -0.5);
      t_s *= corr;
#endif

      DTP(p).StoppingTime = t_s;
      for(int j = 0; j < 3; j++)
        {
          /* Compute physical velocity difference. */
          double dvel = -(P[p].Vel[j] - DTP(p).LocalGasVelocity[j]) / All.cf_atime;
          DTP(p).DragAccel[j] = dvel / t_s;
        }
    }
}

void do_drag_step(void)
{
  int i, j;
  double dvel, dt_base, dt_kick, hubble_a;

  start_dust();
  find_drag_cells(Ndust);
  drag_kernel();
  TIMER_STOPSTART(CPU_DUST, CPU_DUST_DRAG);

  compute_drag_acceleration();

  if(All.ComovingIntegrationOn)
    {
      hubble_a = hubble_function(All.Time);
    }

  for(i = 0; i < Ndust; i++)
    {
      int p = DustParticle[i].index;
      if(P[p].Mass == 0.0)
        continue;

      dt_base = 0.5 * (P[p].TimeBinHydro ? (((integertime) 1) << P[p].TimeBinHydro) : 0) * All.Timebase_interval;

      if(All.ComovingIntegrationOn)
        dt_kick = dt_base / hubble_a;
      else
        dt_kick = dt_base;

      /* Convert drag acceleration back to comoving to apply kicks. */
#ifdef DL_DRAG_SEMI_IMPLICIT
      /* Loren-Aguilar+ 2015, Equation 17 */
      double t_s = DTP(p).StoppingTime;
      double xi = 1.0 - exp(-dt_kick / t_s);
      double fac1 = (dt_kick + t_s) * xi - dt_kick;
      for(j = 0; j < 3; j++)
        {
          double a_dust = P[p].GravAccel[j];
#ifdef PMGRID
          a_dust += P[p].GravPM[j];
#endif
          double fac2 = (a_dust - DustParticle[i].LocalGasAccel[j] + DustParticle[i].LocalGradP[j] / DTP(p).LocalGasDensity);
          double fac3 = -xi * (P[p].Vel[j] - DTP(p).LocalGasVelocity[j]);
          dvel = fac1 * fac2 + fac3;
          /* Could compute an effective physical acceleration, but we would
           * convert right back to comoving below. */
          P[p].Vel[j] += dvel;
        }

#else
      /* Simple explicit update, since dust timesteps have been limited
       * to account for stopping time. */
      for(j = 0; j < 3; j++)
        {
          /* DragAccel was computed as a physical acceleration earlier. */
          dvel = (dt_kick * DTP(p).DragAccel[j]) * All.cf_atime;
          P[p].Vel[j] += dvel;
        }
#endif
    }

  TIMER_STOPSTART(CPU_DUST_DRAG, CPU_DUST);
  end_dust();
}

void do_drag_step_first_half(void)
{
  do_drag_step();
}

void do_drag_step_second_half(void)
{
  do_drag_step();
}

#endif
