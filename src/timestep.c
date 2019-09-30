/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/timestep.c
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
#include "allvars.h"
#include "proto.h"

/*! \file timestep.c
 *  \brief routines for 'kicking' particles in
 *  momentum space and assigning new timesteps
 */


#ifdef SIDM
static int local_sidm_limit_timestep, local_sidm_reduce_timestep;
void reset_timestep_counter(void)
{
  local_sidm_limit_timestep = local_sidm_reduce_timestep = 0;
}
#endif

#if defined(COSMIC_RAYS_DIFFUSION) && defined(COSMIC_RAYS_DIFFUSION_GLOBAL_TIMESTEP)
static void compute_timestep_cr_diffusion();
#endif

#if defined(TURBULENT_METALDIFFUSION) 
void compute_timestep_mm_diffusion();
#endif

#ifdef MONOTONE_CONDUCTION
static void compute_timestep_tc();
#endif

#ifdef IMPLICIT_OHMIC_DIFFUSION
static void compute_timestep_ohm();
static void set_timestep_ohm(integertime globTimeStep);
#endif



/*! This function sets various cosmological factors based on the current simulation time
 */

void set_cosmo_factors_for_current_time(void)
{
  if(All.ComovingIntegrationOn)
    {
      All.cf_atime = All.Time;
      All.cf_a2inv = 1 / (All.Time * All.Time);
      All.cf_a3inv = 1 / (All.Time * All.Time * All.Time);
      All.cf_afac1 = pow(All.Time, 3 * GAMMA_MINUS1);
      All.cf_afac2 = 1 / pow(All.Time, 3 * GAMMA - 2);
      All.cf_afac3 = pow(All.Time, 3 * (1 - GAMMA) / 2.0);
      All.cf_hubble_a = All.cf_H = All.cf_Hrate = hubble_function(All.Time);
      All.cf_time_hubble_a = All.Time * All.cf_hubble_a;
      All.cf_redshift = 1 / All.Time - 1;
    }
  else
    {
      All.cf_atime = 1;
      All.cf_a2inv = 1;
      All.cf_a3inv = 1;
      All.cf_afac1 = 1;
      All.cf_afac2 = 1;
      All.cf_afac3 = 1;
      All.cf_hubble_a = 1;
      All.cf_H = All.Hubble;
      All.cf_time_hubble_a = 1;
      All.cf_Hrate = 0;
      All.cf_redshift = 0;
    }
}

/*! This function advances the system in momentum space, i.e. it does apply the 'kick' operation after the
 *  forces have been computed. Additionally, it assigns new timesteps to particles. At start-up, a
 *  half-timestep is carried out, as well as at the end of the simulation. In between, the half-step kick that
 *  ends the previous timestep and the half-step kick for the new timestep are combined into one operation.
 */

void find_timesteps_without_gravity(void)
{
#ifdef TREE_BASED_TIMESTEPS
  tree_based_timesteps();
#endif

  TIMER_START(CPU_TIMELINE);

  int idx, i, bin, binold;
  integertime ti_step;


#ifdef TGCHEM
  tgchem_compute_step_constants();
#endif

#if defined(COSMIC_RAYS_DIFFUSION) && defined(COSMIC_RAYS_DIFFUSION_GLOBAL_TIMESTEP)
  compute_timestep_cr_diffusion();
#endif

#ifdef TURBULENT_METALDIFFUSION
  compute_timestep_mm_diffusion();
#endif

#ifdef MONOTONE_CONDUCTION
  compute_timestep_tc();
#endif

#if defined(IMPLICIT_OHMIC_DIFFUSION) && !defined(FORCE_EQUAL_TIMESTEPS)
  compute_timestep_ohm();
#endif


#ifdef LEGACY_DISPLACEMENT_CONSTRAINT
  if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin || All.DtDisplacement == 0)
    find_dt_displacement_constraint(All.cf_hubble_a * All.cf_atime * All.cf_atime);
#endif

#ifdef FORCE_EQUAL_TIMESTEPS
  integertime globTimeStep = TIMEBASE;

#ifdef PMGRID
  globTimeStep = get_timestep_pm();
#endif

#if (defined(SELFGRAVITY) || defined(EXTERNALGRAVITY) || defined(EXACT_GRAVITY_FOR_PARTICLE_TYPE) || defined (DG_EXTERNAL_ACCELERATION)) && !defined(MESHRELAX)
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      ti_step = get_timestep_gravity(i);
      if(ti_step < globTimeStep)
        globTimeStep = ti_step;
    }
#endif

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;


      ti_step = get_timestep_hydro(i);
      if(ti_step < globTimeStep)
        globTimeStep = ti_step;
    }

#ifdef BLACK_HOLES
  for(idx = 0; idx < TimeBinsBHAccretion.NActiveParticles; idx++)
    {
      i = TimeBinsBHAccretion.ActiveParticleList[idx];
      if(i < 0)
        continue;

      ti_step = get_timestep_bhaccretion(i);
      if(ti_step < globTimeStep)
        globTimeStep = ti_step;
    }
#endif

#ifdef SINKS
  for(idx = 0; idx < TimeBinsSinksAccretion.NActiveParticles; idx++)
    {
      i = TimeBinsSinksAccretion.ActiveParticleList[idx];
      if(i < 0)
        continue;

      ti_step = get_timestep_sinksaccretion(i);
      if(ti_step < globTimeStep)
        globTimeStep = ti_step;
    }
#endif

#ifdef DUST_LIVE
  start_dust();
  find_drag_cells(Ndust);
  drag_kernel();
  for(idx = 0; idx < TimeBinsDust.NActiveParticles; idx++)
    {
      i = TimeBinsDust.ActiveParticleList[idx];
      if(i < 0)
        continue;

      ti_step = get_timestep_dust(i);
      if(ti_step < globTimeStep)
        globTimeStep = ti_step;
    }
  end_dust();
#endif

#ifdef ENLARGE_DYNAMIC_RANGE_IN_TIME
  minimum_large_ints(1, &globTimeStep, &All.GlobalTimeStep);
#else
  MPI_Allreduce(&globTimeStep, &All.GlobalTimeStep, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
#endif

#ifdef NUCLEAR_NETWORK_TIMESTEP_LIMITER
  if(TimeBinsHydro.NActiveParticles > 0)
    {
      int idx = 0;
      while(idx < TimeBinsHydro.NActiveParticles)
        {
          i = TimeBinsHydro.ActiveParticleList[idx];
          if(i >= 0)
            {
              timebins_get_bin_and_do_validity_checks(All.GlobalTimeStep, &bin, P[i].TimeBinHydro);
              break;
            }
          idx++;
        }
    }
  else
    {
      bin = -1;
    }
  MPI_Allreduce(&bin, &binold, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  
  bin = binold;
  network_main( &bin );
  
  int binnew;
  MPI_Allreduce(&bin, &binnew, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  
  if(binnew < binold)
    {
      if(binnew <= 1)
        terminate( "Timestep <= 1 not allowed." );
      
      All.GlobalTimeStep = binnew ? (((integertime) 1) << binnew) : 0;
      mpi_printf( "NUCLEAR NETWORK: Reducing global timebin from %d to %d, new global timestep=%g.\n", binold, binnew, All.GlobalTimeStep * All.Timebase_interval );
    }
#endif

  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      timebins_get_bin_and_do_validity_checks(All.GlobalTimeStep, &bin, P[i].TimeBinGrav);
      binold = P[i].TimeBinGrav;
      timebin_move_particle(&TimeBinsGravity, i, binold, bin);
      P[i].TimeBinGrav = bin;
    }

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      timebins_get_bin_and_do_validity_checks(All.GlobalTimeStep, &bin, P[i].TimeBinHydro);
      binold = P[i].TimeBinHydro;
      timebin_move_particle(&TimeBinsHydro, i, binold, bin);
      P[i].TimeBinHydro = bin;
    }

#ifdef BLACK_HOLES
  for(idx = 0; idx < TimeBinsBHAccretion.NActiveParticles; idx++)
    {
      i = TimeBinsBHAccretion.ActiveParticleList[idx];
      if(i < 0)
        continue;

      timebins_get_bin_and_do_validity_checks(All.GlobalTimeStep, &bin, P[i].TimeBinHydro);
      binold = P[i].TimeBinHydro;
      timebin_move_particle(&TimeBinsBHAccretion, i, binold, bin);
      P[i].TimeBinHydro = bin;
    }
#endif

#ifdef SINKS
  for(idx = 0; idx < TimeBinsSinksAccretion.NActiveParticles; idx++)
    {
      i = TimeBinsSinksAccretion.ActiveParticleList[idx];
      if(i < 0)
        continue;
      
      timebins_get_bin_and_do_validity_checks(All.GlobalTimeStep, &bin, P[i].TimeBinHydro);
      binold = P[i].TimeBinHydro;
      timebin_move_particle(&TimeBinsSinksAccretion, i, binold, bin);
      P[i].TimeBinHydro = bin;
    }
#endif

#ifdef DUST_LIVE
  for(idx = 0; idx < TimeBinsDust.NActiveParticles; idx++)
    {
      i = TimeBinsDust.ActiveParticleList[idx];
      if(i < 0)
        continue;

      timebins_get_bin_and_do_validity_checks(All.GlobalTimeStep, &bin, P[i].TimeBinHydro);
      binold = P[i].TimeBinHydro;
      timebin_move_particle(&TimeBinsDust, i, binold, bin);
      P[i].TimeBinHydro = bin;
    }
#endif

#ifdef IMPLICIT_OHMIC_DIFFUSION
  set_timestep_ohm(All.GlobalTimeStep);
#endif


#else /* not FORCE_EQUAL_TIMESTEPS */
  /* Calculate and assign hydro timesteps */
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      ti_step = get_timestep_hydro(i);
#ifdef REFINEMENT_MERGE_PAIRS
      /* this needs to be fixed, probably its better to compare timebins instead of timesteps */
      if(SphP[i].DerefPartnerId > 0)
        {
          voronoi_derefinement_update_partner_index(i);
          if(SphP[i].DerefPartnerIndex >= 0)
            {
              double ti_step_partner = get_timestep_hydro(SphP[i].DerefPartnerIndex);
              ti_step = dmin(ti_step, ti_step_partner);
            }
        }
#endif
      timebins_get_bin_and_do_validity_checks(ti_step, &bin, P[i].TimeBinHydro);

      binold = P[i].TimeBinHydro;
      timebin_move_particle(&TimeBinsHydro, i, binold, bin);
      P[i].TimeBinHydro = bin;
    }

#ifdef BLACK_HOLES
  /* Calculate and assign accretion timesteps for black holes */
  for(idx = 0; idx < TimeBinsBHAccretion.NActiveParticles; idx++)
    {
      i = TimeBinsBHAccretion.ActiveParticleList[idx];
      if(i < 0)
        continue;

      ti_step = get_timestep_bhaccretion(i);

      timebins_get_bin_and_do_validity_checks(ti_step, &bin, P[i].TimeBinHydro);

      binold = P[i].TimeBinHydro;
      timebin_move_particle(&TimeBinsBHAccretion, i, binold, bin);
      P[i].TimeBinHydro = bin;
    }
#endif

#ifdef SINKS
  /* Calculate and assign accretion timesteps for black holes */
  for(idx = 0; idx < TimeBinsSinksAccretion.NActiveParticles; idx++)
    {
      i = TimeBinsSinksAccretion.ActiveParticleList[idx];
      if(i < 0)
        continue;

      ti_step = get_timestep_sinksaccretion(i);

      timebins_get_bin_and_do_validity_checks(ti_step, &bin, P[i].TimeBinHydro);

      binold = P[i].TimeBinHydro;
      timebin_move_particle(&TimeBinsSinksAccretion, i, binold, bin);
      P[i].TimeBinHydro = bin;
    }
#endif

#ifdef DUST_LIVE
    update_timebins_dust();
#endif

#endif /* FORCE_EQUAL_TIMESTEPS */

#ifdef SIDM
  long long global_sidm_limit_timestep, global_sidm_reduce_timestep;
  sumup_large_ints(1, &local_sidm_limit_timestep, &global_sidm_limit_timestep);
  sumup_large_ints(1, &local_sidm_reduce_timestep, &global_sidm_reduce_timestep);
  mpi_printf("SIDM: reduced by dt_sidm=%llu   dt_sidm limited by dt=%llu\n", global_sidm_reduce_timestep, global_sidm_limit_timestep);
#endif
  TIMER_STOP(CPU_TIMELINE);
}

void update_timesteps_from_gravity(void)
{
#ifdef FORCE_EQUAL_TIMESTEPS
  return;                       /* don't need to do this */
#endif

#if !((defined(SELFGRAVITY) || defined(EXTERNALGRAVITY) || defined(EXACT_GRAVITY_FOR_PARTICLE_TYPE))) || defined(MESHRELAX)
  return;
#endif

  TIMER_START(CPU_TIMELINE);

  int idx, i, binold;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].TimeBinGrav < P[i].TimeBinHydro)
        {
          binold = P[i].TimeBinHydro;
          timebin_move_particle(&TimeBinsHydro, i, binold, P[i].TimeBinGrav);
          P[i].TimeBinHydro = P[i].TimeBinGrav;
        }
    }

#ifdef BLACK_HOLES
  for(idx = 0; idx < TimeBinsBHAccretion.NActiveParticles; idx++)
    {
      i = TimeBinsBHAccretion.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].TimeBinGrav < P[i].TimeBinHydro)
        {
          binold = P[i].TimeBinHydro;
          timebin_move_particle(&TimeBinsBHAccretion, i, binold, P[i].TimeBinGrav);
          P[i].TimeBinHydro = P[i].TimeBinGrav;
        }
    }
#endif

#ifdef SINKS
  for(idx = 0; idx < TimeBinsSinksAccretion.NActiveParticles; idx++)
    {
      i = TimeBinsSinksAccretion.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].TimeBinGrav < P[i].TimeBinHydro)
        {
          binold = P[i].TimeBinHydro;
          timebin_move_particle(&TimeBinsSinksAccretion, i, binold, P[i].TimeBinGrav);
          P[i].TimeBinHydro = P[i].TimeBinGrav;
        }
    }
#endif

#ifdef DUST_LIVE
  for(idx = 0; idx < TimeBinsDust.NActiveParticles; idx++)
    {
      i = TimeBinsDust.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].TimeBinGrav < P[i].TimeBinHydro)
        {
          binold = P[i].TimeBinHydro;
          timebin_move_particle(&TimeBinsDust, i, binold, P[i].TimeBinGrav);
          P[i].TimeBinHydro = P[i].TimeBinGrav;
        }
    }
#endif

  TIMER_STOP(CPU_TIMELINE);
}

#ifdef PMGRID
integertime get_timestep_pm(void)
{
  integertime ti_step = TIMEBASE;
  while(ti_step > (All.DtDisplacement / All.Timebase_interval))
    ti_step >>= 1;

  if(ti_step > (All.PM_Ti_endstep - All.PM_Ti_begstep)) /* PM-timestep wants to increase */
    {
      int bin = get_timestep_bin(ti_step);
      int binold = get_timestep_bin(All.PM_Ti_endstep - All.PM_Ti_begstep);

      while(TimeBinSynchronized[bin] == 0 && bin > binold)      /* make sure the new step is synchronized */
        bin--;

      ti_step = bin ? (((integertime) 1) << bin) : 0;
    }

  if(All.Ti_Current == TIMEBASE)        /* we here finish the last timestep. */
    ti_step = 0;

  return ti_step;
}
#endif



integertime get_timestep_gravity(int p)
{
  double dt;
  integertime ti_step;

  double ax, ay, az, ac;

#ifdef SIDM
    double dt_sidm, dt_sidm_tmp;
    int state, reactiong;
#endif

#ifdef DG
if(P[p].Type == 0)
  dt = dg_time_step_cell_gravity(p);
else
#endif
  {

    /* calculate total acceleration */
    ax = All.cf_a2inv * P[p].GravAccel[0];
    ay = All.cf_a2inv * P[p].GravAccel[1];
    az = All.cf_a2inv * P[p].GravAccel[2];

#if defined(PMGRID) && !defined(NO_PMFORCE_IN_SHORT_RANGE_TIMESTEP)
    ax += All.cf_a2inv * P[p].GravPM[0];
    ay += All.cf_a2inv * P[p].GravPM[1];
    az += All.cf_a2inv * P[p].GravPM[2];
#endif

#ifdef AB_TURB
    if(P[p].Type == 0)
    {
      ax += All.cf_a2inv * SphP[p].TurbAccel[0];
      ay += All.cf_a2inv * SphP[p].TurbAccel[1];
      az += All.cf_a2inv * SphP[p].TurbAccel[2];
    }
#endif

    ac = sqrt(ax * ax + ay * ay + az * az);       /* this is now the physical acceleration */
    if(ac == 0)
      ac = 1.0e-30;

    switch (All.TypeOfTimestepCriterion)
    {
    case 0:
      /* only type 0 implemented at the moment -> remove type ? */
      dt = sqrt(2 * All.ErrTolIntAccuracy * All.cf_atime * All.ForceSoftening[P[p].SofteningType] / 2.8 / ac);
      break;

    default:
      terminate("Undefined timestep criterion");
      break;
    }
  }

#ifdef SIDM
  /* constraint time step to avoid too many multiple scatter events (prefactor 0.1 can be adopted) */
  if((1 << P[p].Type) & (SIDM))
    {
      dt_sidm = MAX_DOUBLE_NUMBER;

      for(state = 0; state < SIDM_STATES; state++)
        for(reactiong = 0; reactiong < SIDM_REACTIONS; reactiong++)
          {
            if(P[p].sidm_Density[state] > 0.0 && P[p].sidm_VelDisp[state] > 0.0)
              {
                dt_sidm_tmp =
                  All.DtimeFac / ((All.HubbleParam * All.HubbleParam * All.cf_a3inv * P[p].sidm_Density[state]) * sidm_cross_sigma(P[p].sidm_VelDisp[state] / All.cf_atime, reactiong) *
                                  (P[p].sidm_VelDisp[state] / All.cf_atime));
                if(dt_sidm_tmp < dt_sidm)
                  dt_sidm = dt_sidm_tmp;
              }
          }

      /* how much smaller can dt_sidm be compared to dt_grav (to limit extreme cases where dt_sidm would get extremely small). This behaviour is controlled by All.DtimeFacLim --> limit timestep */
      if(dt_sidm < All.DtimeFacLim * dt)
        {
          dt_sidm = All.DtimeFacLim * dt;
          local_sidm_limit_timestep++;
        }
#ifndef SIDM_NO_TIMESTEP
      /* check whether particle requires a smaller timestep than dt due to self-interactions --> reduce timestep */
      if(dt_sidm < dt)
        {
          dt = dt_sidm;
          local_sidm_reduce_timestep++;
        }
#endif
    }
#endif

#ifdef OTVET_FIXTIMESTEP
  dt = All.MaxSizeTimestep;
#endif

#ifdef BECDM
  double cosmofac = 1;
  if(All.ComovingIntegrationOn)
    cosmofac = All.cf_atime * All.cf_atime;  // a^2
  if( dt >= cosmofac * All.mAxion/(6.0*All.hbar) * pow( All.BoxSize / ((double) PMGRID),2) )
    dt = cosmofac * All.mAxion/(6.0*All.hbar) * pow( All.BoxSize / ((double) PMGRID),2);
#endif

  dt *= All.cf_hubble_a;

  if(P[p].Mass == 0 && P[p].ID == 0)
    dt = All.MaxSizeTimestep;   /* this particle has been swallowed or eliminated */

  if(dt >= All.MaxSizeTimestep)
    dt = All.MaxSizeTimestep;

  if(dt < All.MinSizeTimestep)
    {
#ifdef NOSTOP_WHEN_BELOW_MINTIMESTEP
      dt = All.MinSizeTimestep;
#else
      print_particle_info(p);
      terminate("Timestep dt=%g below All.MinSizeTimestep=%g", dt, All.MinSizeTimestep);
#endif
    }

#ifdef PMGRID
  if(dt >= All.DtDisplacement)
    dt = All.DtDisplacement;
#endif

  ti_step = (integertime) (dt / All.Timebase_interval);
  validate_timestep(dt, ti_step, p);

  return ti_step;
}

#ifdef BLACK_HOLES
integertime get_timestep_bhaccretion(int p)
{
  double dt, dt_courant, dt_accr;
  integertime ti_step;

  dt_courant = All.CourantFac * BPP(p).BH_DtGasNeighbor;
  if(All.ComovingIntegrationOn)
    dt_courant *= All.Time;

  dt = dt_courant;

  if(BPP(p).BH_Mdot > 0 && BPP(p).BH_Mass > 0)
    {
      dt_accr = 0.25 * BPP(p).BH_Mass / BPP(p).BH_Mdot;
      if(dt_accr < dt)
        dt = dt_accr;
    }


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

#ifdef SINKS
integertime get_timestep_sinksaccretion(int p)
{
  integertime ti_step;
  double dt;

  dt = (P[p].TimeBinSink ? (((integertime) 1) << P[p].TimeBinSink) : 0) * All.Timebase_interval;

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

integertime get_timestep_hydro(int p)
{
  double dt = 0, dt_courant = 0;
  integertime ti_step;

  assert(P[p].Type == 0);

#ifndef DG
  double csnd = get_sound_speed(p);

#if defined(VORONOI_STATIC_MESH) || defined(AMR)
  csnd += sqrt(P[p].Vel[0] * P[p].Vel[0] + P[p].Vel[1] * P[p].Vel[1] + P[p].Vel[2] * P[p].Vel[2]);
#else
#ifdef SPECIAL_BOUNDARY
  csnd += sqrt((P[p].Vel[0] - SphP[p].VelVertex[0]) * (P[p].Vel[0] - SphP[p].VelVertex[0]) +
               (P[p].Vel[1] - SphP[p].VelVertex[1]) * (P[p].Vel[1] - SphP[p].VelVertex[1]) + (P[p].Vel[2] - SphP[p].VelVertex[2]) * (P[p].Vel[2] - SphP[p].VelVertex[2]));
#endif
#endif

#ifdef SPECIAL_RELATIVITY
  csnd = get_cfl_sound_speed_special_relativity(p);
#endif

#ifdef GENERAL_RELATIVITY
  csnd = get_cfl_sound_speed_general_relativity(p);
#endif
  
#if defined(COSMIC_RAYS_STREAMING) && !defined(COSMIC_RAYS_STREAMING_EXPLICIT)
  csnd += sqrt( (SphP[p].B[0]*SphP[p].B[0] + SphP[p].B[1]*SphP[p].B[1] + SphP[p].B[2]*SphP[p].B[2]) / SphP[p].Density );
#endif
  
  double rad = get_cell_radius(p);

  if(csnd <= 0)
    csnd = 1.0e-30;

  dt_courant = rad / csnd;

#ifdef TREE_BASED_TIMESTEPS
  if(dt_courant > SphP[p].CurrentMaxTiStep)
    dt_courant = SphP[p].CurrentMaxTiStep;
#endif

  dt_courant *= All.CourantFac;
#else
  dt_courant = dg_time_step_cell(p);
#endif

  if(All.ComovingIntegrationOn)
    dt_courant *= All.Time;

  dt = dt_courant;

#if defined(USE_SFR) && !defined(LOCAL_FEEDBACK) && !defined(QUICK_LYALPHA)

  if(P[p].Type == 0)            /* to protect using a particle that has been turned into a star */
    {
      double sfr = get_starformation_rate(p);
      double mass_loading = 0;
      if(sfr > 0)
        {
#ifdef GFM_WINDS
#ifdef GFM_WINDS_VARIABLE
#if (GFM_WINDS_VARIABLE==0)
          if(SphP[p].w.HostHaloMass > 0 && P[p].Mass > 0)
            {
              wind_parameter wp;
              gfm_calc_variable_wind_parameters(P[p].Mass, SphP[p].w.HostHaloMass, SphP[p].w.DMVelDisp, SphP[p].Metallicity, &wp);
#endif
#if (GFM_WINDS_VARIABLE==1)
              if(SphP[p].w.DMVelDisp > 0 && P[p].Mass > 0)
                {
                  wind_parameter wp;
                  gfm_calc_variable_wind_parameters(P[p].Mass, SphP[p].w.DMVelDisp, SphP[p].w.DMVelDisp, SphP[p].Metallicity, &wp);
#endif
                  mass_loading = wp.wind_mass / P[p].Mass;
                }
#else
          mass_loading = All.WindEfficiency;
#endif
#endif
        }

      double dt_sfr = 0.1 * P[p].Mass / ((1 + mass_loading) * sfr / ((All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR)));
      if(dt_sfr < dt)
        dt = dt_sfr;
    }
#endif

#ifdef TGCHEM
  dt = tgchem_get_timestep(p, dt);
#endif

#ifdef OTVET_FIXTIMESTEP
  dt = All.MaxSizeTimestep;
#endif

#ifdef MRT
  double chi = c_internal_units ;
#if defined(MRT_COMOVING) && defined(VORONOI_STATIC_MESH)
  chi += sqrt(P[p].Vel[0]*P[p].Vel[0] + P[p].Vel[1]*P[p].Vel[1] + P[p].Vel[2]*P[p].Vel[2]) ;
#endif
  double dt_rt = ((double) (All.RTNumSubCycles)) * 0.3 * rad / chi ;
  //mpi_printf("%d\t%g\t%g\t%g\n", p, chi, dt_rt, dt) ;
  if(dt_rt < dt)
    dt = dt_rt ;
#endif

#ifdef NON_IDEAL_MHD
#if defined(OHMIC_DIFFUSION) || defined(AMBIPOLAR_DIFFUSION)
  double dt_nonidealmhd = MAX_REAL_NUMBER;
#ifdef OHMIC_DIFFUSION
  double diff_coeff = dmax(1e-30, All.OhmicDiffusionCoefficient);
  double dt_tmp = 0.2 * rad * rad / diff_coeff;
  if(dt_tmp < dt_nonidealmhd)
    dt_nonidealmhd = dt_tmp;
#endif
#ifdef AMBIPOLAR_DIFFUSION
  double dt_tmp = 0.2 * rad * rad / All.AmbipolarDiffusionCoefficient;
  if(dt_tmp < dt_nonidealmhd)
    dt_nonidealmhd = dt_tmp;
#endif

  if(All.ComovingIntegrationOn)
    dt_nonidealmhd *= All.Time * All.Time;

#ifdef NON_IDEAL_MHD_EXPLICIT_LIMIT_TIMESTEP
  double Bpress = 0.5 * (SphP[p].B[0] * SphP[p].B[0] + SphP[p].B[1] * SphP[p].B[1] + SphP[p].B[2] * SphP[p].B[2]);
  double v_alfven = sqrt(2.0 * Bpress / SphP[p].Density / All.cf_atime);
  double dt_alfven = All.NonidealMHDTimelimit * rad / v_alfven;
  if (dt_nonidealmhd < dt_alfven)
    dt_nonidealmhd = dt_alfven;
#endif

  if(dt_nonidealmhd < dt)
    dt = dt_nonidealmhd;
#endif

#endif

#ifdef MHD_POWELL_LIMIT_TIMESTEP
  double b = sqrt(SphP[p].B[0] * SphP[p].B[0] + SphP[p].B[1] * SphP[p].B[1] + SphP[p].B[2] * SphP[p].B[2]);
  double bmin = sqrt(2 * 0.01 * SphP[p].Utherm * SphP[p].Density * All.cf_atime);
  double v = sqrt(P[p].Vel[0] * P[p].Vel[0] + P[p].Vel[1] * P[p].Vel[1] + P[p].Vel[2] * P[p].Vel[2]) / All.cf_atime;
  double dt_powell = 0.5 * (b + bmin) / (fabs(SphP[p].DivB / All.cf_atime * v));

  if(dt_powell < dt)
    dt = dt_powell;
#endif


  /* convert the physical timestep to dloga if needed. Note: If comoving integration has not been selected,
     All.cf_hubble_a=1.
   */

  dt *= All.cf_hubble_a;

  if(dt >= All.MaxSizeTimestep)
    dt = All.MaxSizeTimestep;

#ifdef TIMESTEP_OUTPUT_LIMIT
  if(dt >= All.TimestepOutputLimit)
    dt = All.TimestepOutputLimit;
#endif

  if(dt < All.MinSizeTimestep)
    {
#ifdef NOSTOP_WHEN_BELOW_MINTIMESTEP
      dt = All.MinSizeTimestep;
#else
      print_particle_info(p);
      terminate("Timestep dt=%g below All.MinSizeTimestep=%g", dt, All.MinSizeTimestep);
#endif
    }


#ifdef PMGRID
  if(dt >= All.DtDisplacement)
    dt = All.DtDisplacement;
#endif


#ifdef COSMIC_RAYS_DIFFUSION
#ifndef COSMIC_RAYS_DIFFUSION_GLOBAL_TIMESTEP
#ifdef COSMIC_RAYS_DIFFUSION_EXPLICIT
  double coeff = 1e28;
  if(All.CR_Diffusion_Coefficient > 0)
    coeff = All.CR_Diffusion_Coefficient;
  coeff = coeff / (All.UnitLength_in_cm * All.UnitLength_in_cm) * All.UnitTime_in_s;

  double dt_diff = rad * rad / coeff;
  if(dt_diff < dt)
    dt = dt_diff;

  /* if we run the implicit solver, but without a global cr timestep, the timesteps
     are only determined by the hydro */
#endif
#else
  if(dt >= All.dt_cr_diffusion)
    dt = All.dt_cr_diffusion;
#endif
#endif

#if defined(COSMIC_RAYS_STREAMING) && defined(COSMIC_RAYS_STREAMING_EXPLICIT)
  double dt_streaming = get_timestep_streaming(p);
  if(dt_streaming < dt)
    dt = dt_streaming;
#endif

#ifdef TURBULENT_METALDIFFUSION
  if(dt >= All.dt_metaldiff)
    dt = All.dt_metaldiff;
#endif

#ifdef MONOTONE_CONDUCTION
  if(dt >= All.dt_max_conduction)
    dt = All.dt_max_conduction;
#endif

#if defined(IMPLICIT_OHMIC_DIFFUSION) && !defined(FORCE_EQUAL_TIMESTEPS)
  if(dt >= All.dt_max_ohmdiffusion)
    dt = All.dt_max_ohmdiffusion;
#endif


#ifdef SINK_PARTICLES
  if(P[p].Mass == 0 && P[p].ID == 0)
    dt = All.MaxSizeTimestep;   /* this particle has been swallowed or eliminated */
#endif


  ti_step = (integertime) (dt / All.Timebase_interval);
  validate_timestep(dt, ti_step, p);

  return ti_step;
}

#if defined(TURBULENT_METALDIFFUSION) 
void compute_timestep_mm_diffusion()
{
  if(All.metaldiff_Ti_endstep == All.Ti_Current)
    {
      mpi_printf("TURBULENT_METALDIFFUSION: All.metaldiff_Ti_endstep,All.Ti_Current=%d,%d\n",All.metaldiff_Ti_endstep,All.Ti_Current);
      double dt = MAX_DOUBLE_NUMBER;

#ifdef TURBULENT_METALDIFFUSION_EXPLICIT
      coeff = All.metaldiff_kappamax / (All.UnitLength_in_cm * All.UnitLength_in_cm) * All.UnitTime_in_s;
#endif
      int igas;
      for(igas = 0; igas < NumGas; igas++)
        {
#ifdef TURBULENT_METALDIFFUSION_EXPLICIT
          double rad = get_cell_radius(igas);
          double dt_diff = rad * rad / coeff;

          if(dt_diff < dt)
            dt = dt_diff;
#else
          double csnd = get_sound_speed(igas);

#ifdef VORONOI_STATIC_MESH
          csnd += sqrt(P[igas].Vel[0] * P[igas].Vel[0] + P[igas].Vel[1] * P[igas].Vel[1] + P[igas].Vel[2] * P[igas].Vel[2]);
#endif

          double rad = get_cell_radius(igas);

          if(csnd <= 0)
            csnd = 1.0e-30;

          double dt_courant = All.Courantmetaldiff * rad / csnd;

          if(dt_courant < dt)
            dt = dt_courant;
#endif
        }

      integertime ti_step;

#ifdef PMGRID
      integertime ti_step = (integertime) get_timestep_pm();
      dt = ti_step * All.Timebase_interval;
#endif

      if(dt >= All.MaxSizeTimestep)
        dt = All.MaxSizeTimestep;

      MPI_Allreduce(&dt, &All.dt_metaldiff, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

      ti_step = (integertime) (All.dt_metaldiff / All.Timebase_interval);

      int bin, binold;
      binold = get_timestep_bin(All.metaldiff_Ti_endstep - All.metaldiff_Ti_begstep);
      timebins_get_bin_and_do_validity_checks(ti_step, &bin, binold);

      ti_step = bin ? (((integertime) 1) << bin) : 0;

      All.metaldiff_Ti_begstep = All.metaldiff_Ti_endstep;                                         
      All.metaldiff_Ti_endstep = All.metaldiff_Ti_begstep + ti_step;

      mpi_printf("TURBULENT_METALDIFFUSION: dt=%g, bin=%d, ti_step=%d\n", All.dt_metaldiff, bin, ti_step);
      mpi_printf("TURBULENT_METALDIFFUSION: begstep = %llu\tendstep = %llu\tti_step = %llu\n", (long long) All.metaldiff_Ti_begstep, (long long) All.metaldiff_Ti_endstep,(long long) ti_step) ;

    }
}
#endif
 


#if defined(COSMIC_RAYS_DIFFUSION) && defined(COSMIC_RAYS_DIFFUSION_GLOBAL_TIMESTEP)
void compute_timestep_cr_diffusion()
{
  if(All.cr_diffusion_Ti_endstep == All.Ti_Current)
    {
#ifdef COSMIC_RAYS_DIFFUSION_CONSTANT_TIMESTEP
      All.dt_cr_diffusion = All.CRDiffusionTimestep;
#else
      double dt = MAX_DOUBLE_NUMBER;

      int igas;
      for(igas = 0; igas < NumGas; igas++)
        {
#ifdef COSMIC_RAYS_DIFFUSION_EXPLICIT
          double coeff = 1e28;
          if(All.CR_Diffusion_Coefficient > 0)
            coeff = All.CR_Diffusion_Coefficient;
          coeff = coeff / (All.UnitLength_in_cm * All.UnitLength_in_cm) * All.UnitTime_in_s;

          double rad = get_cell_radius(igas);
          double dt_diff = rad * rad / coeff;

          if(dt_diff < dt)
            dt = dt_diff;
#else
          double csnd = get_sound_speed(igas);

#ifdef VORONOI_STATIC_MESH
          csnd += sqrt(P[igas].Vel[0] * P[igas].Vel[0] + P[igas].Vel[1] * P[igas].Vel[1] + P[igas].Vel[2] * P[igas].Vel[2]);
#endif

          double rad = get_cell_radius(igas);

          if(csnd <= 0)
            csnd = 1.0e-30;

          double dt_courant = All.CourantCRDiffusion * rad / csnd;

          if(dt_courant < dt)
            dt = dt_courant;
#endif
        }

      if(dt >= All.MaxSizeTimestep)
        dt = All.MaxSizeTimestep;

      MPI_Allreduce(&dt, &All.dt_cr_diffusion, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif
      integertime ti_step = (integertime) (All.dt_cr_diffusion / All.Timebase_interval);

      int bin, binold;
      binold = get_timestep_bin(All.cr_diffusion_Ti_endstep - All.cr_diffusion_Ti_begstep);
      timebins_get_bin_and_do_validity_checks(ti_step, &bin, binold);

      ti_step = bin ? (((integertime) 1) << bin) : 0;

      All.cr_diffusion_Ti_begstep = All.cr_diffusion_Ti_endstep;
      All.cr_diffusion_Ti_endstep = All.cr_diffusion_Ti_begstep + ti_step;

      mpi_printf("COSMIC_RAYS: DIFFUSION: dt=%g, bin=%d, ti_step=%d\n", All.dt_cr_diffusion, bin, ti_step);
    }
}
#endif


#ifdef MONOTONE_CONDUCTION
void compute_timestep_tc()
{
  if(All.conduction_Ti_endstep == All.Ti_Current)
    {
      
      double dt = MAX_DOUBLE_NUMBER;
      double ncfl = 4.0;

#ifdef RESTRICT_KAPPA
      double alpha_max = All.Chi_max_int/All.cf_atime/All.cf_atime ;
#endif


      int igas;
      for(igas = 0; igas < NumGas; igas++)
        {

          double rad = get_cell_radius(igas);


#ifdef CONDUCTION_CONSTANT
          double alpha = All.ConductionCoeff;
#else
          double alpha = All.ConductionCoeff * pow(SphP[igas].Utherm, 2.5) / SphP[igas].Density;
#ifdef CONDUCTION_SATURATION
          double electron_free_path = All.ElectronFreePathFactor * SphP[igas].Utherm * SphP[igas].Utherm / (SphP[igas].Density * All.cf_a3inv);

          double mod_dU =
            sqrt(SphP[igas].Grad.dutherm[0] * SphP[igas].Grad.dutherm[0] + SphP[igas].Grad.dutherm[1] * SphP[igas].Grad.dutherm[1] + SphP[igas].Grad.dutherm[2] * SphP[igas].Grad.dutherm[2]);

          if(mod_dU != 0.0)
            {
              double temp_scale_length = All.Time * fabs(SphP[igas].Utherm) / mod_dU;

              alpha /= (1 + 4.2 * electron_free_path / temp_scale_length);
            }
#endif
#endif


          if(All.ComovingIntegrationOn)
            alpha *= All.cf_atime;

#ifdef RESTRICT_KAPPA
	  if(alpha > alpha_max)
	    alpha = alpha_max ;
#endif

          double dt_courant = ncfl * rad * rad / alpha ;

#ifdef USE_SFR
	  if(SphP[igas].Sfr > 0.0)
	    dt_courant = MAX_DOUBLE_NUMBER ;
#endif

          if(dt_courant < dt)
            dt = dt_courant;
        }


      MPI_Allreduce(&dt, &All.dt_conduction, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

      All.dt_max_conduction = All.dt_conduction * All.cf_hubble_a * ((double) (All.MaxConductionSubCycles));
      //All.dt_conduction *= All.cf_hubble_a;

      /*if(All.dt_conduction >= All.MaxSizeTimestep)
         All.dt_conduction = All.MaxSizeTimestep; */

      if(All.dt_max_conduction >= All.MaxSizeTimestep)
        All.dt_max_conduction = All.MaxSizeTimestep;

      //integertime ti_step = (integertime) (All.dt_conduction / All.Timebase_interval);


      integertime ti_step = (integertime) (All.dt_max_conduction / All.Timebase_interval);

      int bin, binold;
      binold = get_timestep_bin(All.conduction_Ti_endstep - All.conduction_Ti_begstep);
      timebins_get_bin_and_do_validity_checks(ti_step, &bin, binold);

      ti_step = bin ? (((integertime) 1) << bin) : 0;

      All.conduction_Ti_begstep = All.conduction_Ti_endstep;
      All.conduction_Ti_endstep = All.conduction_Ti_begstep + ti_step;

      mpi_printf("CONDUCTION: begstep = %llu\tendstep = %llu\tti_step = %llu\n", (long long) All.conduction_Ti_begstep, (long long) All.conduction_Ti_endstep, (long long) ti_step) ;
      mpi_printf("CONDUCTION: dt=%g, dt_single=%g, bin=%d, ti_step=%llu\n", All.dt_max_conduction, All.dt_conduction * All.cf_hubble_a, bin, (long long) ti_step);
    }
}
#endif

#ifdef IMPLICIT_OHMIC_DIFFUSION
void compute_timestep_ohm()
{
  if(All.ohmdiffusion_Ti_endstep == All.Ti_Current)
    {
      double dt = MAX_DOUBLE_NUMBER;
#ifdef OHM_CRANK_NICHOLSON
      double ncfl = 0.5;
#else
      double ncfl = 4.0;
#endif

      for(int igas = 0; igas < NumGas; igas++)
        {

          double rad = get_cell_radius(igas);

          double alpha = All.OhmicDiffusionCoefficient;

          if(All.ComovingIntegrationOn)
            alpha /= All.cf_atime * All.cf_atime;

          double dt_courant = ncfl * rad * rad / alpha;

          if(dt_courant < dt)
            dt = dt_courant;
        }


      MPI_Allreduce(&dt, &All.dt_ohmdiffusion, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

      All.dt_max_ohmdiffusion = All.dt_ohmdiffusion * All.cf_hubble_a;

      if(All.dt_max_ohmdiffusion >= All.MaxSizeTimestep)
        All.dt_max_ohmdiffusion = All.MaxSizeTimestep;

      integertime ti_step = (integertime) (All.dt_max_ohmdiffusion / All.Timebase_interval);

      int bin, binold;
      binold = get_timestep_bin(All.ohmdiffusion_Ti_endstep - All.ohmdiffusion_Ti_begstep);
      timebins_get_bin_and_do_validity_checks(ti_step, &bin, binold);

      ti_step = bin ? (((integertime) 1) << bin) : 0;

      All.ohmdiffusion_Ti_begstep = All.ohmdiffusion_Ti_endstep;
      All.ohmdiffusion_Ti_endstep = All.ohmdiffusion_Ti_begstep + ti_step;

      mpi_printf("NON-IDEAL MHD: Ohm diffusion, begstep = %llu endstep = %llu ti_step = %llu\n", (long long)All.ohmdiffusion_Ti_begstep, (long long)All.ohmdiffusion_Ti_endstep, (long long)ti_step);
      mpi_printf("NON-IDEAL MHD: Ohm diffusion, dt=%g, dt_single=%g, bin=%d, ti_step=%llu\n", All.dt_max_ohmdiffusion, All.dt_ohmdiffusion * All.cf_hubble_a, bin, (long long)ti_step);
    }
}

void set_timestep_ohm(integertime globTimeStep)
{
  if(All.ohmdiffusion_Ti_endstep == All.Ti_Current)
    {
      int bin, binold;
      integertime ti_step = globTimeStep;
      binold = get_timestep_bin(All.ohmdiffusion_Ti_endstep - All.ohmdiffusion_Ti_begstep);
      timebins_get_bin_and_do_validity_checks(ti_step, &bin, binold);

      ti_step = bin ? (((integertime) 1) << bin) : 0;

      All.ohmdiffusion_Ti_begstep = All.ohmdiffusion_Ti_endstep;
      All.ohmdiffusion_Ti_endstep = All.ohmdiffusion_Ti_begstep + ti_step;

      mpi_printf("NON-IDEAL MHD: Ohm diffusion, begstep = %llu endstep = %llu ti_step = %llu\n", (long long)All.ohmdiffusion_Ti_begstep, (long long)All.ohmdiffusion_Ti_endstep, (long long)ti_step);
    }
}
#endif

void validate_timestep(double dt, integertime ti_step, int p)
{
  if(!(ti_step > 0 && ti_step < TIMEBASE))
    {
      printf("\nError: A timestep of size zero was assigned on the integer timeline!\n" "We better stop.\n" "Task=%d Part-ID=%lld type=%d", ThisTask, (long long) P[p].ID, P[p].Type);

      printf("tibase=%g dt=%g ti_step=%d, xyz=(%g|%g|%g) vel=(%g|%g|%g) tree=(%g|%g|%g) mass=%g\n\n",
             All.Timebase_interval, dt, (int) ti_step, P[p].Pos[0], P[p].Pos[1], P[p].Pos[2], P[p].Vel[0], P[p].Vel[1], P[p].Vel[2], P[p].GravAccel[0], P[p].GravAccel[1], P[p].GravAccel[2], P[p].Mass);

      print_particle_info(p);
      myflush(stdout);
      terminate("integer timestep reached zero");
    }

  if(ti_step == 1)
    {
      printf("Time-step of integer size 1 found for particle i=%d, pos=(%g|%g|%g), ID=%lld, dt=%g\n", p, P[p].Pos[0], P[p].Pos[1], P[p].Pos[2], (long long) P[p].ID, dt);
      print_particle_info(p);
    }
}

int test_if_grav_timestep_is_too_large(int p, int bin)
{
  integertime ti_step_bin = bin ? (((integertime) 1) << bin) : 0;

  integertime ti_step = get_timestep_gravity(p);

  if(P[p].Type == 0
#if defined(BLACK_HOLES) || defined(SINKS)
    || P[p].Type == 5
#endif
#ifdef DUST_LIVE
    || P[p].Type == DUST_LIVE
#endif
  )
    {
      if((P[p].ID != 0) && (P[p].Mass != 0))
        {
          int bin_hydro = P[p].TimeBinHydro + MAX_TIMEBIN_DIFFERENCE;
          integertime ti_step_hydro = bin_hydro ? (((integertime) 1) << bin_hydro) : 0;
          if(ti_step_hydro < ti_step)
            ti_step = ti_step_hydro;
        }
    }

  if(ti_step < ti_step_bin)
    return 1;
  else
    return 0;
}


#ifndef LEGACY_DISPLACEMENT_CONSTRAINT
#ifdef PMGRID
void find_long_range_step_constraint(void)
{
  int p;
  double ax, ay, az, ac;
  double dt, dtmin = MAX_DOUBLE_NUMBER;

  for(p = 0; p < NumPart; p++)
    {
      if(P[p].Type == 0)
        continue;

#ifdef PM_TIMESTEP_BASED_ON_TYPES
      if(((1 << P[p].Type) & (PM_TIMESTEP_BASED_ON_TYPES)))
#endif
	{
	  /* calculate acceleration */
	  ax = All.cf_a2inv * P[p].GravPM[0];
	  ay = All.cf_a2inv * P[p].GravPM[1];
	  az = All.cf_a2inv * P[p].GravPM[2];

	  ac = sqrt(ax * ax + ay * ay + az * az);   /* this is now the physical acceleration */

	  if(ac < MIN_FLOAT_NUMBER)
	    ac = MIN_FLOAT_NUMBER;

	  dt = sqrt(2.0 * All.ErrTolIntAccuracy * All.cf_atime * All.ForceSoftening[P[p].SofteningType] / (2.8 * ac));

	  dt *= All.cf_hubble_a;

	  if(dt < dtmin)
	    dtmin = dt;
	}
    }

  dtmin *= 2.0;                 /* move it one timebin higher to prevent being too conservative */

  MPI_Allreduce(&dtmin, &All.DtDisplacement, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  mpi_printf("TIMESTEPS: displacement time constraint: %g  (%g)\n", All.DtDisplacement, All.MaxSizeTimestep);

  if(All.DtDisplacement > All.MaxSizeTimestep)
    All.DtDisplacement = All.MaxSizeTimestep;

  if(All.DtDisplacement < All.MinSizeTimestep)
    All.DtDisplacement = All.MinSizeTimestep;
}
#endif
#endif







#ifdef LEGACY_DISPLACEMENT_CONSTRAINT
/*! This function computes an upper limit ('dt_displacement') to the global timestep of the system based on
 *  the rms velocities of particles. For cosmological simulations, the criterion used is that the rms
 *  displacement should be at most a fraction MaxRMSDisplacementFac of the mean particle separation. Note that
 *  the latter is estimated using the assigned particle masses, separately for each particle type. If comoving
 *  integration is not used, the function imposes no constraint on the timestep.
 */
void find_dt_displacement_constraint(double hfac /*!<  should be  a^2*H(a)  */ )
{
  int i, type;
  int count[NTYPES];
  long long count_sum[NTYPES];
  double v[NTYPES], v_sum[NTYPES], msum[NTYPES], sum_mass[NTYPES], mim[NTYPES], min_mass[NTYPES];
  double dt, dmean, asmth = 0;

  All.DtDisplacement = All.MaxSizeTimestep;

  if(All.ComovingIntegrationOn)
    {
      for(type = 0; type < NTYPES; type++)
        {
          count[type] = 0;
          v[type] = 0;
          msum[type] = 0;
          mim[type] = 1.0e30;
        }

      for(i = 0; i < NumPart; i++)
        {
          v[P[i].Type] += P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2];

          msum[P[i].Type] += P[i].Mass;
          if(P[i].Mass > 0)
            {
              if(mim[P[i].Type] > P[i].Mass)
                mim[P[i].Type] = P[i].Mass;
            }
          count[P[i].Type]++;
        }

      MPI_Allreduce(v, v_sum, NTYPES, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(msum, sum_mass, NTYPES, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(mim, min_mass, NTYPES, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

      sumup_large_ints(NTYPES, count, count_sum);

#ifdef USE_SFR
      /* add star and gas particles together to treat them on equal footing, using the original gas particle
         spacing. */
      v_sum[0] += v_sum[4];
      count_sum[0] += count_sum[4];
      sum_mass[0] += sum_mass[4];
      v_sum[4] = v_sum[0];
      count_sum[4] = count_sum[0];
      sum_mass[4] = sum_mass[0];
#ifdef BLACK_HOLES
      v_sum[0] += v_sum[5];
      count_sum[0] += count_sum[5];
      sum_mass[0] += sum_mass[5];
      v_sum[5] = v_sum[0];
      count_sum[5] = count_sum[0];
      sum_mass[5] = sum_mass[0];
// OLD CODE: does not work, see below
//      min_mass[5] = min_mass[0];     
#endif
#endif

      for(type = 0; type < NTYPES; type++)
        {
          if(count_sum[type] > 0)
            {
              if(type == 0 || (type == 4 && All.StarformationOn))
                dmean = pow(sum_mass[type] / count_sum[type] / (All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)), 1.0 / 3);
              else
                dmean = pow(sum_mass[type] / count_sum[type] / ((All.Omega0 - All.OmegaBaryon) * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)), 1.0 / 3);

#ifdef BLACK_HOLES
// OLDE CODE: does not work if cells get strongly refined such that min_mass[5]=min_mass[0] gets very small leading to small min_mass[0]
//              if(type == 5)
//                dmean = pow(min_mass[type] / (All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)), 1.0 / 3);
              if(type == 5)
                dmean = pow(sum_mass[type] / count_sum[type] / (All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)), 1.0 / 3);
#endif
              dt = All.MaxRMSDisplacementFac * hfac * dmean / sqrt(v_sum[type] / count_sum[type]);

#ifdef PMGRID
              asmth = All.Asmth[0];
#ifdef PLACEHIGHRESREGION
              if(((1 << type) & (PLACEHIGHRESREGION)))
                asmth = All.Asmth[1];
#endif
              if(asmth < dmean)
                dt = All.MaxRMSDisplacementFac * hfac * asmth / sqrt(v_sum[type] / count_sum[type]);
#endif

              if(ThisTask == 0)
                printf
                  ("TIMESTEPS: type=%d  dmean=%g asmth=%g avgmass=%g sum_mass=%g n=%d  a=%g  sqrt(<p^2>)=%g  dlogmax=%g\n",
                   type, dmean, asmth, sum_mass[type] / count_sum[type], sum_mass[type], (int) count_sum[type], All.Time, sqrt(v_sum[type] / count_sum[type]), dt);


#ifdef TRACER_PARTICLE
              if(type != TRACER_PARTICLE)
#endif
                if(dt < All.DtDisplacement)
                  All.DtDisplacement = dt;
            }
        }

      mpi_printf("TIMESTEPS: displacement time constraint: %g  (%g)\n", All.DtDisplacement, All.MaxSizeTimestep);
    }
}
#endif

int get_timestep_bin(integertime ti_step)
{
  int bin = -1;

  if(ti_step == 0)
    return 0;

  if(ti_step == 1)
    terminate("time-step of integer size 1 not allowed\n");

  while(ti_step)
    {
      bin++;
      ti_step >>= 1;
    }

  return bin;
}

double get_time_difference_in_Gyr(double a0, double a1)
{
  double result, time_diff = 0, t0, t1, factor1, factor2, term1, term2;

  if(All.ComovingIntegrationOn)
    {
      if(All.OmegaLambda + All.Omega0 != 1)
        printf("only implemented for flat cosmology so far.");

      factor1 = 2.0 / (3.0 * sqrt(All.OmegaLambda));

      term1 = sqrt(All.OmegaLambda / All.Omega0) * pow(a0, 1.5);
      term2 = sqrt(1 + All.OmegaLambda / All.Omega0 * pow(a0, 3));
      factor2 = log(term1 + term2);

      t0 = factor1 * factor2;

      term1 = sqrt(All.OmegaLambda / All.Omega0) * pow(a1, 1.5);
      term2 = sqrt(1 + All.OmegaLambda / All.Omega0 * pow(a1, 3));
      factor2 = log(term1 + term2);

      t1 = factor1 * factor2;

      result = t1 - t0;

      time_diff = result / (HUBBLE * All.HubbleParam);  /* now in seconds */
      time_diff /= SEC_PER_MEGAYEAR * 1000;     /* now in gigayears */
    }
  else
    {
      time_diff = (a1 - a0) * All.UnitTime_in_s / All.HubbleParam;      /* now in seconds */
      time_diff /= SEC_PER_MEGAYEAR * 1000;     /* now in gigayears */
    }

  return time_diff;
}

void timebins_init(struct TimeBinData *tbData, const char *name, int *MaxPart)
{
  int i;
  tbData->NActiveParticles = 0;
  tbData->ActiveParticleList = 0;

  for(i = 0; i < TIMEBINS; i++)
    {
      tbData->FirstInTimeBin[i] = -1;
      tbData->LastInTimeBin[i] = -1;
    }

  tbData->NextInTimeBin = 0;
  tbData->PrevInTimeBin = 0;

  strncpy(tbData->Name, name, 99);
  tbData->Name[99] = 0;
  tbData->MaxPart = MaxPart;
}

void timebins_allocate(struct TimeBinData *tbData)
{
  char Identifier[200];
  Identifier[199] = 0;

  snprintf(Identifier, 199, "NextActiveParticle%s", tbData->Name);
  tbData->ActiveParticleList = (int *) mymalloc_movable(&tbData->ActiveParticleList, Identifier, *(tbData->MaxPart) * sizeof(int));

  snprintf(Identifier, 199, "NextInTimeBin%s", tbData->Name);
  tbData->NextInTimeBin = (int *) mymalloc_movable(&tbData->NextInTimeBin, Identifier, *(tbData->MaxPart) * sizeof(int));

  snprintf(Identifier, 199, "PrevInTimeBin%s", tbData->Name);
  tbData->PrevInTimeBin = (int *) mymalloc_movable(&tbData->PrevInTimeBin, Identifier, *(tbData->MaxPart) * sizeof(int));
}

void timebins_reallocate(struct TimeBinData *tbData)
{
  tbData->ActiveParticleList = (int *) myrealloc_movable(tbData->ActiveParticleList, *(tbData->MaxPart) * sizeof(int));
  tbData->NextInTimeBin = (int *) myrealloc_movable(tbData->NextInTimeBin, *(tbData->MaxPart) * sizeof(int));
  tbData->PrevInTimeBin = (int *) myrealloc_movable(tbData->PrevInTimeBin, *(tbData->MaxPart) * sizeof(int));
}

void timebins_get_bin_and_do_validity_checks(integertime ti_step, int *bin_new, int bin_old)
{
  /* make it a power 2 subdivision */
  integertime ti_min = TIMEBASE;
  while(ti_min > ti_step)
    ti_min >>= 1;
  ti_step = ti_min;

  /* get timestep bin */
  int bin = -1;

  if(ti_step == 0)
    bin = 0;

  if(ti_step == 1)
    terminate("time-step of integer size 1 not allowed\n");

  while(ti_step)
    {
      bin++;
      ti_step >>= 1;
    }

  if(bin > bin_old)             /* timestep wants to increase */
    {
      while(TimeBinSynchronized[bin] == 0 && bin > bin_old)     /* make sure the new step is synchronized */
        bin--;

      ti_step = bin ? (((integertime) 1) << bin) : 0;
    }

  if(All.Ti_Current >= TIMEBASE)        /* we here finish the last timestep. */
    {
      ti_step = 0;
      bin = 0;
    }

  if((TIMEBASE - All.Ti_Current) < ti_step)     /* check that we don't run beyond the end */
    {
      terminate("we are beyond the end of the timeline");       /* should not happen */
    }

  *bin_new = bin;
}

void timebin_move_particle(struct TimeBinData *tbData, int p, int timeBin_old, int timeBin_new)
{
  if(timeBin_old == timeBin_new)
    return;

  tbData->TimeBinCount[timeBin_old]--;

  int prev = tbData->PrevInTimeBin[p];
  int next = tbData->NextInTimeBin[p];

  if(tbData->FirstInTimeBin[timeBin_old] == p)
    tbData->FirstInTimeBin[timeBin_old] = next;
  if(tbData->LastInTimeBin[timeBin_old] == p)
    tbData->LastInTimeBin[timeBin_old] = prev;
  if(prev >= 0)
    tbData->NextInTimeBin[prev] = next;
  if(next >= 0)
    tbData->PrevInTimeBin[next] = prev;

  if(tbData->TimeBinCount[timeBin_new] > 0)
    {
      tbData->PrevInTimeBin[p] = tbData->LastInTimeBin[timeBin_new];
      tbData->NextInTimeBin[tbData->LastInTimeBin[timeBin_new]] = p;
      tbData->NextInTimeBin[p] = -1;
      tbData->LastInTimeBin[timeBin_new] = p;
    }
  else
    {
      tbData->FirstInTimeBin[timeBin_new] = tbData->LastInTimeBin[timeBin_new] = p;
      tbData->PrevInTimeBin[p] = tbData->NextInTimeBin[p] = -1;
    }

  tbData->TimeBinCount[timeBin_new]++;

#ifdef USE_SFR
  if((P[p].Type == 0) && (tbData == &TimeBinsHydro))
    timebin_move_sfr(p, timeBin_old, timeBin_new);
#endif

#ifdef BLACK_HOLES
  if((P[p].Type == 5) && (tbData == &TimeBinsBHAccretion))
    timebin_move_bh(p, timeBin_old, timeBin_new);
#endif
}

void timebin_remove_particle(struct TimeBinData *tbData, int idx, int bin)
{
  int p = tbData->ActiveParticleList[idx];
  tbData->ActiveParticleList[idx] = -1;

  if(bin == -1)
    {
      if(tbData == &TimeBinsGravity)
        bin = P[p].TimeBinGrav;
      else
        bin = P[p].TimeBinHydro;
    }

  tbData->TimeBinCount[bin]--;

  if(p >= 0)
    {
      int prev = tbData->PrevInTimeBin[p];
      int next = tbData->NextInTimeBin[p];

      if(prev >= 0)
        tbData->NextInTimeBin[prev] = next;
      if(next >= 0)
        tbData->PrevInTimeBin[next] = prev;

      if(tbData->FirstInTimeBin[bin] == p)
        tbData->FirstInTimeBin[bin] = next;
      if(tbData->LastInTimeBin[bin] == p)
        tbData->LastInTimeBin[bin] = prev;
    }
}

/* insert a particle into the timebin struct behind another already existing particle */
void timebin_add_particle(struct TimeBinData *tbData, int i_new, int i_old, int timeBin, int addToListOfActiveParticles)
{
  tbData->TimeBinCount[timeBin]++;

  if(i_old < 0)
    {
      /* if we don't have an existing particle to add if after, let's take the last one in this timebin */
      i_old = tbData->LastInTimeBin[timeBin];

      if(i_old < 0)
        {
          /* the timebin is empty at the moment, so just add the new particle */
          tbData->FirstInTimeBin[timeBin] = i_new;
          tbData->LastInTimeBin[timeBin] = i_new;
          tbData->NextInTimeBin[i_new] = -1;
          tbData->PrevInTimeBin[i_new] = -1;
        }
    }

  if(i_old >= 0)
    {
      /* otherwise we added it already */
      tbData->PrevInTimeBin[i_new] = i_old;
      tbData->NextInTimeBin[i_new] = tbData->NextInTimeBin[i_old];
      if(tbData->NextInTimeBin[i_old] >= 0)
        tbData->PrevInTimeBin[tbData->NextInTimeBin[i_old]] = i_new;
      tbData->NextInTimeBin[i_old] = i_new;
      if(tbData->LastInTimeBin[timeBin] == i_old)
        tbData->LastInTimeBin[timeBin] = i_new;
    }

  if(addToListOfActiveParticles)
    {
      tbData->ActiveParticleList[tbData->NActiveParticles] = i_new;
      tbData->NActiveParticles++;
    }
}

void timebin_cleanup_list_of_active_particles(struct TimeBinData *tbData)
{
  int idx, i;
  for(idx = 0; idx < tbData->NActiveParticles; idx++)
    {
      i = tbData->ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].ID == 0 && P[i].Mass == 0)
        timebin_remove_particle(tbData, idx, -1);
    }
}

#ifdef USE_SFR
void timebin_move_sfr(int p, int timeBin_old, int timeBin_new)
{
  TimeBinSfr[timeBin_old] -= SphP[p].Sfr;
  TimeBinSfr[timeBin_new] += SphP[p].Sfr;
}
#endif

#ifdef BLACK_HOLES
void timebin_move_bh(int p, int timeBin_old, int timeBin_new)
{
  TimeBin_BH_mass[timeBin_old] -= BPP(p).BH_Mass;
  TimeBin_BH_dynamicalmass[timeBin_old] -= P[p].Mass;
  TimeBin_BH_Mdot[timeBin_old] -= BPP(p).BH_Mdot;
  if(BPP(p).BH_Mass > 0)
    TimeBin_BH_Medd[timeBin_old] -= BPP(p).BH_Mdot / BPP(p).BH_Mass;

  TimeBin_BH_mass[timeBin_new] += BPP(p).BH_Mass;
  TimeBin_BH_dynamicalmass[timeBin_new] += P[p].Mass;
  TimeBin_BH_Mdot[timeBin_new] += BPP(p).BH_Mdot;
  if(BPP(p).BH_Mass > 0)
    TimeBin_BH_Medd[timeBin_new] += BPP(p).BH_Mdot / BPP(p).BH_Mass;
}
#endif

void timebin_make_list_of_active_particles_up_to_timebin(struct TimeBinData *tbData, int timebin)
{
  int tbin;
  tbData->NActiveParticles = 0;
  for(tbin = timebin; tbin >= 0; tbin--)
    timebin_add_particles_of_timebin_to_list_of_active_particles(tbData, tbin);
}

void timebin_add_particles_of_timebin_to_list_of_active_particles(struct TimeBinData *tbData, int timebin)
{
  int i;
  for(i = tbData->FirstInTimeBin[timebin]; i >= 0; i = tbData->NextInTimeBin[i])
    if(!(P[i].ID == 0 && P[i].Mass == 0))
      {
        tbData->ActiveParticleList[tbData->NActiveParticles] = i;
        tbData->NActiveParticles++;
      }
}
