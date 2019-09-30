/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dust_live/dust_timestep.c
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

/* Update dust time bins of active dust particles using minimum hydro time bin
 * of nearby gas cells.
 */
void update_timebins_dust(void)
{
  start_dust();
  find_drag_cells(Ndust);
  drag_kernel();
  for(int i = 0; i < Ndust; i++)
    {
      int p = DustParticle[i].index;
      int bin_old = P[p].TimeBinHydro;
      int bin_new;
      integertime ti_step = get_timestep_dust(p);

      /* Since we have a time bin directly, and not dt, we could just handle
       * synchronization manually, but we'll use the helper for its validity
       * checks and error handling. */
      timebins_get_bin_and_do_validity_checks(ti_step, &bin_new, bin_old);

      timebin_move_particle(&TimeBinsDust, p, bin_old, bin_new);
      P[p].TimeBinHydro = bin_new;
    }
  end_dust();
}

integertime get_timestep_dust(int p)
{
  double dt; //, dt_courant;
  integertime ti_step;

  /* Start with timestep corresponding to minimum of nearby gas cells. */
  // TODO: comoving check
  double dt_gas_min = (DTP(p).MinGasTimeBin ? (((integertime) 1) << DTP(p).MinGasTimeBin) : 0) * All.Timebase_interval;
  // TODO: are these line below necessary?  MinGasTimeBin seems to account for All.CourantFac
  // since the gas bins already have it
  //dt_courant = All.CourantFac * dt_gas_min;
  //if(All.ComovingIntegrationOn)
  //  dt_courant *= All.Time;
  //dt = dt_courant;

  dt = dt_gas_min;

  double rel_vel2 = 0.0;
  for(int j = 0; j < 3; j++)
    {
      double dvel = (DTP(p).LocalGasVelocity[j] - P[p].Vel[j]) / All.cf_atime;
      rel_vel2 += (dvel * dvel);
    }

  double csnd_eff = sqrt(pow(DTP(p).LocalSoundSpeed, 2) + rel_vel2);
  if(csnd_eff <= 0.0)
    csnd_eff = 1.0e-30;

  double dt_courant = All.CourantFac * DTP(p).Hsml / csnd_eff;
  if(All.ComovingIntegrationOn)
    dt_courant *= All.Time;

  if(dt_courant < dt)
    dt = dt_courant;

#ifndef DL_DRAG_SEMI_IMPLICIT
  /* If we require explicit drag timesteps, make sure to resolve the stopping timescale. */
  if((DTP(p).StoppingTime > 0.0) && (dt > DTP(p).StoppingTime * All.StoppingTimeFrac))
    {
      dt = DTP(p).StoppingTime * All.StoppingTimeFrac;
    }
#endif

#ifdef DL_GRAIN_BINS
#if defined(DL_GROWTH) || defined(DL_SPUTTERING) || defined(DL_SNE_DESTRUCTION) || defined(DL_SHATTERING) || defined(DL_COAGULATION)
  /* If we are evolving a grain size distribution, choose a timestep limit
   * based on the maximum allowable fractional change in mass in the grain size
   * bins. */
  if(dt > All.MaxBinFracMassChg * DTP(p).BinMassChgTau)
    {
      dt = All.MaxBinFracMassChg * DTP(p).BinMassChgTau;
    }
  /* Reset in preparation for next timestep. */
  DTP(p).BinMassChgTau = MAX_REAL_NUMBER;

  /* Store current mass for use in timescale calculations taking place during
   * grain size evolution. */
  DTP(p).OrigMass = P[p].Mass;
#endif
#endif

  dt *= All.cf_hubble_a;

  if(dt >= All.MaxSizeTimestep)
    dt = All.MaxSizeTimestep;

#ifdef PMGRID
  if(dt >= All.DtDisplacement)
    dt = All.DtDisplacement;
#endif

  ti_step = (integertime) (dt / All.Timebase_interval);
  validate_timestep(dt, ti_step, p);

  return ti_step;
}

#endif
