/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/GFM/stellar_feedback.c
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

#include "stellar_feedback_kernels.h"

#if defined(GFM_STELLAR_FEEDBACK) || defined(GFM_WINDS_LOCAL) || defined(FM_STAR_FEEDBACK)


static data_in *DataIn, *DataGet;
static data_out *DataResult, *DataOut;

#ifdef FM_STAR_FEEDBACK
#ifndef DELAYED_COOLING_TURB
#ifdef DELAYED_COOLING
static int sum_kicked_cells, kicked_cells;
static double max_blast_radius, min_blast_radius;
static float global_max_blast_radius, global_min_blast_radius;
static float max_shut_time, min_shut_time;
static float global_max_shut_time, global_min_shut_time;
#else
static int kicked_cells;
static double massToKick, massKicked;
static float maxKickVel, minKickVel;
#if !defined(NON_STOCHASTIC_MOMENTUM_FEEDBACK) && !defined(DIRECT_MOMENTUM_INJECTION_FEEDBACK)
static int sum_kicked_cells;
static double sumMassToKick, sumMassKicked;
static double totalMassToKick = 0.0, totalMassKicked = 0.0;
static float global_maxKickVel, global_minKickVel;
static double global_max_energy_diff, global_min_energy_diff;
static double global_max_momentum_diff, global_min_momentum_diff;
#endif
#endif
#endif
#endif


static void particle2in(data_in * in, int i, int firstnode)
{
  in->Pos[0] = P[StarParticle[i].index].Pos[0];
  in->Pos[1] = P[StarParticle[i].index].Pos[1];
  in->Pos[2] = P[StarParticle[i].index].Pos[2];

  in->Hsml = STP(StarParticle[i].index).Hsml;
  in->NormSph = StarParticle[i].NormSph;
  in->TotalMassReleased = StarParticle[i].TotalMassReleased;

#if defined(FM_MASS_WEIGHT_SN) || defined(NON_STOCHASTIC_MOMENTUM_FEEDBACK)
  in->TotNgbMass = StarParticle[i].TotNgbMass;
#endif
#ifdef GFM_WINDS_LOCAL
  in->WindEnergyReleased = StarParticle[i].WindEnergyReleased;
#endif
#ifdef GFM_STELLAR_FEEDBACK
  in->SNIaEnergyReleased = StarParticle[i].SNIaEnergyReleased;
  in->AGBMomentumReleased = StarParticle[i].AGBMomentumReleased;
#endif
#ifdef FM_STAR_FEEDBACK
  in->TotalEnergyReleased = StarParticle[i].TotalEnergyReleased;
#ifdef FM_SN_COOLING_RADIUS_BOOST
  in->n_SNII = StarParticle[i].NumSNII;
  in->n_SNIa = StarParticle[i].NumSNIa;
  in->LocISMdens   = StarParticle[i].LocISMdens;      /* local ISM density (code units) */
  in->LocISMZdens  = StarParticle[i].LocISMZdens;     /* local ISM metallicity          */
#endif
#ifdef DIRECT_MOMENTUM_INJECTION_FEEDBACK
  in->TotalMomentumReleased = StarParticle[i].TotalMomentumReleased;
#endif
#if defined(NON_STOCHASTIC_MOMENTUM_FEEDBACK) && defined(INJECT_INTO_SINGLE_CELL)
  in->ClosestNeighbourID = StarParticle[i].ClosestNeighbourID;
#endif
#ifdef DELAYED_COOLING
  in->NormSphFeedback = StarParticle[i].NormSphFeedback;
  in->BlastRadius = StarParticle[i].BlastRadius;
  in->CoolShutoffTime = StarParticle[i].CoolShutoffTime;
#else
#ifndef DELAYED_COOLING_TURB
  in->NumNgb = StarParticle[i].NumNgb;
#endif
#endif
#endif

  in->Firstnode = firstnode;
}

static void out2particle(data_out * out, int i, int mode)
{
#ifdef FM_STAR_FEEDBACK
  if(mode == MODE_LOCAL_PARTICLES)
    {
#if !defined DELAYED_COOLING_TURB && !defined(NON_STOCHASTIC_MOMENTUM_FEEDBACK) && !defined(DIRECT_MOMENTUM_INJECTION_FEEDBACK)
#ifdef DELAYED_COOLING
      kicked_cells += out->kicked_cells;
      max_blast_radius = dmax(max_blast_radius, out->BlastRadius);
      min_blast_radius = dmin(min_blast_radius, out->BlastRadius);
      max_shut_time = dmax(max_shut_time, out->ShutoffTime);
      min_shut_time = dmin(min_shut_time, out->ShutoffTime);
#else
      kicked_cells += out->kicked_cells;
      massToKick += out->mass_to_kick;
      massKicked += out->mass_kicked;
      minKickVel = dmin(minKickVel, out->kick_vel);
      maxKickVel = dmax(maxKickVel, out->kick_vel);
      StarParticle[i].deltaEKin = out->deltaEKin;
      StarParticle[i].deltaMomentum[0] = out->deltaMomentum[0];
      StarParticle[i].deltaMomentum[1] = out->deltaMomentum[1];
      StarParticle[i].deltaMomentum[2] = out->deltaMomentum[2];
#endif
#endif
#if !defined(DELAYED_COOLING_TURB) && !defined (DELAYED_COOLING)
      StarParticle[i].TotalMomentumInjected = out->TotalMomentumInjected;
#endif
    }
  else
    {
#if !defined DELAYED_COOLING_TURB && !defined(NON_STOCHASTIC_MOMENTUM_FEEDBACK) && !defined(DIRECT_MOMENTUM_INJECTION_FEEDBACK)
#ifdef DELAYED_COOLING
      kicked_cells += out->kicked_cells;
      max_blast_radius = dmax(max_blast_radius, out->BlastRadius);
      min_blast_radius = dmin(min_blast_radius, out->BlastRadius);
      max_shut_time = dmax(max_shut_time, out->ShutoffTime);
      min_shut_time = dmin(min_shut_time, out->ShutoffTime);
#else
      kicked_cells += out->kicked_cells;
      massKicked += out->mass_kicked;
      minKickVel = dmin(minKickVel, out->kick_vel);
      maxKickVel = dmax(maxKickVel, out->kick_vel);
      StarParticle[i].deltaEKin += out->deltaEKin;
      StarParticle[i].deltaMomentum[0] += out->deltaMomentum[0];
      StarParticle[i].deltaMomentum[1] += out->deltaMomentum[1];
      StarParticle[i].deltaMomentum[2] += out->deltaMomentum[2];
#endif
#endif
#if !defined(DELAYED_COOLING_TURB) && !defined (DELAYED_COOLING)
      StarParticle[i].TotalMomentumInjected += out->TotalMomentumInjected;
#endif
    }
#endif
}


#include "../generic_comm_helpers2.h"

static void kernel_local(void)
{
  int i;
#ifdef GENERIC_ASYNC
  int flag = 0;
#endif
#pragma omp parallel private(idx)
  {
    int j, threadid = get_thread_num();
#ifdef GENERIC_ASYNC
    int count = 0;
#endif

    for(j = 0; j < NTask; j++)
      Thread[threadid].Exportflag[j] = -1;

    while(1)
      {
        if(Thread[threadid].ExportSpace < MinSpace)
          break;

#ifdef GENERIC_ASYNC
        if(threadid == 0)
          {
            if((count & POLLINGINTERVAL) == 0)
              if(generic_polling_primary(count, Nstar))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif


#pragma omp atomic capture
        i = NextParticle++;

        if(i >= Nstar)
          break;

        if(is_doing_stellar_feedback(i))
          stellar_feedback_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
      }
  }
}

static void kernel_imported(void)
{
  /* now do the particles that were sent to us */
  int i, cnt = 0;
#pragma omp parallel private(i)
  {
    int threadid = get_thread_num();
#ifdef GENERIC_ASYNC
    int count = 0;
#endif

    while(1)
      {
#pragma omp atomic capture
        i = cnt++;

        if(i >= Nimport)
          break;

#ifdef GENERIC_ASYNC
        if(threadid == 0)
          {
            if((count & POLLINGINTERVAL) == 0)
              generic_polling_secondary();
          }

        count++;
#endif

        stellar_feedback_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

void do_stellar_feedback(void)
{
  long long ntot;

  sumup_large_ints(1, &Nstar, &ntot);
  if(ntot == 0)
    return;

  TIMER_STOPSTART(CPU_GFM_ENRICH, CPU_GFM_FEEDBACK);

#ifdef GFM_WINDS_LOCAL
  for(int i = 0; i < NumGas; i++)
    SphP[i].WindEnergyReceived = 0;
#endif

#ifdef FM_STAR_FEEDBACK
#ifndef DELAYED_COOLING_TURB
#ifdef DELAYED_COOLING
  kicked_cells = 0;

  if(Nstar > 0)
    {
      max_blast_radius = -MAX_REAL_NUMBER;
      min_blast_radius = MAX_REAL_NUMBER;
      max_shut_time = -MAX_REAL_NUMBER;
      min_shut_time = MAX_REAL_NUMBER;

      /* it can happen that a star is active but no SN exploded */
      for(int i = 0; i < Nstar; i++)
        if(StarParticle[i].TotalEnergyReleased <= 0)
          {
            max_blast_radius = min_blast_radius = 0.0;
            max_shut_time = min_shut_time = 0.0;
          }
    }
  else
    {
      max_blast_radius = min_blast_radius = 0.0;
      max_shut_time = min_shut_time = 0.0;
    }
#else
  kicked_cells = 0;
  massToKick = massKicked = 0.0;

  if(Nstar > 0)
    {
      maxKickVel = -MAX_REAL_NUMBER;
      minKickVel = MAX_REAL_NUMBER;

      /* it can happen that a star is active but no SN exploded */
      for(int i = 0; i < Nstar; i++)
        if(StarParticle[i].TotalEnergyReleased <= 0)
          minKickVel = maxKickVel = 0.0;
    }
  else
    maxKickVel = minKickVel = 0.0;
#endif
#endif
#endif

  mpi_printf("GFM_FEEDBACK: Begin stellar feedback calculation.\n");

  generic_set_MaxNexport();

  double t0 = second();

  generic_comm_pattern(Nstar, kernel_local, kernel_imported);

#ifdef FM_STAR_FEEDBACK
#ifndef DELAYED_COOLING_TURB
#if defined(NON_STOCHASTIC_MOMENTUM_FEEDBACK) || defined(DIRECT_MOMENTUM_INJECTION_FEEDBACK)
#ifdef OUTPUT_STELLAR_FEEDBACK
  for(int i = 0; i < Nstar; i++)
    {
      STP(StarParticle[i].index).FeedbackMomentum = StarParticle[i].TotalMomentumInjected;
      STP(StarParticle[i].index).Cum_FeedbackMomentum += StarParticle[i].TotalMomentumInjected;
    }
#endif
#else
#ifdef DELAYED_COOLING
  float kicked_cells_per_star = 0;
  float max_shutoff = 0.0, min_shutoff = 0.0;
  float global_max_shutoff, global_min_shutoff;
  int active_stars = 0;
  int tot_active_stars = 0;
  int no_cool_cells = 0;
  int tot_no_cool_cells = 0;

  /* count only active star particles that had SN explosions */
  for(int i = 0; i < Nstar; i++)
    if(StarParticle[i].TotalEnergyReleased > 0)
      active_stars++;

  if(active_stars > 0)
    {
      max_shutoff = -MAX_REAL_NUMBER;
      min_shutoff = MAX_REAL_NUMBER;
    }

  /* count gas particles with deactivated cooling */
  for(int i = 0; i < NumGas; i++)
    if(SphP[i].CoolShutoffTime > 0)
      {
        max_shutoff = dmax(max_shutoff, SphP[i].CoolShutoffTime);
        min_shutoff = dmin(min_shutoff, SphP[i].CoolShutoffTime);
        no_cool_cells++;
      }

  MPI_Reduce(&kicked_cells, &sum_kicked_cells, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&active_stars, &tot_active_stars, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&no_cool_cells, &tot_no_cool_cells, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&min_blast_radius, &global_min_blast_radius, 1, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&max_blast_radius, &global_max_blast_radius, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&min_shut_time, &global_min_shut_time, 1, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&max_shut_time, &global_max_shut_time, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&min_shutoff, &global_min_shutoff, 1, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&max_shutoff, &global_max_shutoff, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);

  if(tot_active_stars > 0)
    kicked_cells_per_star = ((float) sum_kicked_cells / tot_active_stars);

  mpi_printf
    ("GFM_FEEDBACK: delayed cooling feedback (this step)  ---> number of kicked cells = %d, active stars (with SN explosions) = %d, kicked cells / star = %g\n",
     sum_kicked_cells, tot_active_stars, kicked_cells_per_star);
  mpi_printf("GFM_FEEDBACK: delayed cooling feedback (this step)  ---> number of cells with deactivated cooling = %d\n", tot_no_cool_cells);
  mpi_printf("GFM_FEEDBACK: delayed cooling feedback (this step)  ---> min blast radius (int units)= %g, max blast radius (int units)= %g\n", global_min_blast_radius, global_max_blast_radius);
  mpi_printf
    ("GFM_FEEDBACK: delayed cooling feedback (this step)  ---> as computed by feedback min shutoff time (int units)= %g, max shutoff time (int units)= %g\n",
     global_min_shut_time, global_max_shut_time);
  mpi_printf
    ("GFM_FEEDBACK: delayed cooling feedback (this step)  ---> as computed in SphP struct  min shutoff time (int units)= %g, max shutoff time (int units)= %g\n",
     global_min_shutoff, global_max_shutoff);
#else
  double energy_diff, energy_released, max_energy_diff, min_energy_diff;
  double momentum_diff, max_momentum_diff, min_momentum_diff;
  float avg_mass_kicked = 0.0, kicked_cells_per_star = 0;
  int active_stars = 0;
  int tot_active_stars = 0;
  if(Nstar > 0)
    {
      max_energy_diff = -MAX_REAL_NUMBER;
      min_energy_diff = MAX_REAL_NUMBER;
      max_momentum_diff = -MAX_REAL_NUMBER;
      min_momentum_diff = MAX_REAL_NUMBER;

      /* count only active star particles that had SN explosions */
      for(int i = 0; i < Nstar; i++)
        {
#ifdef OUTPUT_STELLAR_FEEDBACK
          STP(StarParticle[i].index).FeedbackMomentum = StarParticle[i].TotalMomentumInjected;
          STP(StarParticle[i].index).Cum_FeedbackMomentum += StarParticle[i].TotalMomentumInjected;
#endif
          if(StarParticle[i].TotalEnergyReleased > 0)
            {
              energy_released = All.EtaKineticEnergy * StarParticle[i].TotalEnergyReleased * All.cf_atime * All.cf_atime;
              energy_diff = (StarParticle[i].deltaEKin - energy_released) / energy_released;
              momentum_diff = sqrt(StarParticle[i].deltaMomentum[0] * StarParticle[i].deltaMomentum[0] +
                                   StarParticle[i].deltaMomentum[1] * StarParticle[i].deltaMomentum[1] + StarParticle[i].deltaMomentum[2] * StarParticle[i].deltaMomentum[2]);
              active_stars++;
            }
          else
            {
              energy_diff = 0.0;        // because no SN had occurred
              momentum_diff = 0.0;
              //continue;
            }

          max_energy_diff = dmax(max_energy_diff, energy_diff);
          min_energy_diff = dmin(min_energy_diff, energy_diff);
          max_momentum_diff = dmax(max_momentum_diff, momentum_diff);
          min_momentum_diff = dmin(min_momentum_diff, momentum_diff);
        }
    }
  else
    {
      max_energy_diff = 0.0;
      min_energy_diff = 0.0;
      max_momentum_diff = 0.0;
      min_momentum_diff = 0.0;
    }

  MPI_Reduce(&kicked_cells, &sum_kicked_cells, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&massToKick, &sumMassToKick, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&massKicked, &sumMassKicked, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&active_stars, &tot_active_stars, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&minKickVel, &global_minKickVel, 1, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&maxKickVel, &global_maxKickVel, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&max_energy_diff, &global_max_energy_diff, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&min_energy_diff, &global_min_energy_diff, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&max_momentum_diff, &global_max_momentum_diff, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&min_momentum_diff, &global_min_momentum_diff, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

  totalMassToKick += sumMassToKick;
  totalMassKicked += sumMassKicked;

  if(sum_kicked_cells > 0)
    avg_mass_kicked = sumMassKicked / sum_kicked_cells;

  if(tot_active_stars > 0)
    kicked_cells_per_star = ((float) sum_kicked_cells / tot_active_stars);

  float cells_to_kick = sumMassToKick / All.TargetGasMass;

  mpi_printf
    ("GFM_FEEDBACK: kinetic feedback (this step)  ---> number of kicked cells = %d, active stars (with SN explosions) = %d, kicked cells / star = %g\n",
     sum_kicked_cells, tot_active_stars, kicked_cells_per_star);
  mpi_printf("GFM_FEEDBACK: kinetic feedback (this step)  ---> number of cells to kick = %g, cells to kick / star = %g\n", cells_to_kick, cells_to_kick / tot_active_stars);
  mpi_printf("GFM_FEEDBACK: kinetic feedback (this step)  ---> total mass to kick = %g, total kicked mass = %g, avg mass kicked = %g\n", sumMassToKick, sumMassKicked, avg_mass_kicked);
  mpi_printf("                                                 target gas mass = %g\n", All.TargetGasMass);
  mpi_printf("GFM_FEEDBACK: kinetic feedback (cumulative) ---> total mass to kick = %g, total kicked mass = %g\n", totalMassToKick, totalMassKicked);
  mpi_printf("GFM_FEEDBACK: kinetic feedback (this step)  ---> max kick velocity = %g, min kick velocity = %g\n", global_maxKickVel, global_minKickVel);
  mpi_printf("GFM_FEEDBACK: kinetic feedback (this step)  ---> max energy diff = %g, min energy diff = %g\n", global_max_energy_diff, global_min_energy_diff);
  mpi_printf("GFM_FEEDBACK: kinetic feedback (this step)  ---> max momentum diff = %g, min momentum diff = %g\n", global_max_momentum_diff, global_min_momentum_diff);
#endif
#endif
#endif
#endif

  double t1 = second();

  mpi_printf("GFM_FEEDBACK: stellar feedback calculation done took %g sec\n", timediff(t0, t1));

  TIMER_STOPSTART(CPU_GFM_FEEDBACK, CPU_GFM_ENRICH);
}


int stellar_feedback_evaluate(int target, int mode, int thread_id)
{
  int numnodes, *firstnode;
  data_in local, *in;
  data_out out;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      in = &local;

      numnodes = 1;
      firstnode = NULL;
    }
  else
    {
      in = &DataGet[target];
      generic_get_numnodes(target, &numnodes, &firstnode);
    }

#ifdef GFM_STELLAR_FEEDBACK
  GFM_stellar_feedback(target, mode, thread_id, numnodes, firstnode, in);
#endif
#ifdef GFM_WINDS_LOCAL
  GFM_winds_local(target, mode, thread_id, numnodes, firstnode, in);
#endif
#ifdef FM_STAR_FEEDBACK
#ifdef DELAYED_COOLING
  delayed_cooling(target, mode, thread_id, numnodes, firstnode, in, &out);
#endif
#ifdef DELAYED_COOLING_TURB
  delayed_cooling_turbulence(target, mode, thread_id, numnodes, firstnode, in);
#endif
#ifdef FM_SN_COOLING_RADIUS_BOOST
  cooling_radius_momentum_feedback( target, mode, thread_id, numnodes, firstnode, in, &out);
#endif
#ifdef NON_STOCHASTIC_MOMENTUM_FEEDBACK
  momentum_feedback(target, mode, thread_id, numnodes, firstnode, in, &out);
#endif
#ifdef DIRECT_MOMENTUM_INJECTION_FEEDBACK
  direct_momentum_feedback(target, mode, thread_id, numnodes, firstnode, in, &out);
#endif
#if !defined(DELAYED_COOLING) && !defined(DELAYED_COOLING_TURB) && !defined(NON_STOCHASTIC_MOMENTUM_FEEDBACK) && !defined(DIRECT_MOMENTUM_INJECTION_FEEDBACK)
  stochastic_momentum_feedback(target, mode, thread_id, numnodes, firstnode, in, &out);
#endif
#endif

  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#endif
