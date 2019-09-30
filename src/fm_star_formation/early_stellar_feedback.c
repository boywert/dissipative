/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/fm_star_formation/early_stellar_feedback.c
 * \date        05/2014
 * \author      Laura Sales
 * \brief        
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - 14.05.2014 Description
 *   Implements an stochastic early feedback sampling the gas particles within the stellar kernel
 *   Velocity of wind is computed based on local escape speed with mass and radius within kernel
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "../allvars.h"
#include "../proto.h"


#if defined(FM_EARLY_STAR_FEEDBACK)

typedef struct
{
  MyDouble Pos[3];
  MyFloat Hsml;
  MyFloat NormSph;
  MyDouble EarlyTotalEnergyReleased;
  MyFloat NumNgb;
  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

typedef struct
{
  int kicked_cells;
  MyDouble mass_kicked;
  MyDouble mass_to_kick;
  MyFloat kick_vel;
  MyDouble deltaEarlyEKin;
  MyDouble deltaEarlyMomentum[3];
  MyDouble EarlyMomentumInjected;
} data_out;

static data_out *DataResult, *DataOut;

static int sum_kicked_cells, kicked_cells;
static double massToKick, massKicked;
static double sumMassToKick, sumMassKicked;
static double totalMassToKick = 0.0, totalMassKicked = 0.0;
static double maxKickVel, minKickVel, global_maxKickVel, global_minKickVel;
static double global_max_energy_diff, global_min_energy_diff;
static double global_max_momentum_diff, global_min_momentum_diff;

static void particle2in(data_in * in, int i, int firstnode)
{
  in->Pos[0] = P[StarParticle[i].index].Pos[0];
  in->Pos[1] = P[StarParticle[i].index].Pos[1];
  in->Pos[2] = P[StarParticle[i].index].Pos[2];

  in->Hsml = STP(StarParticle[i].index).Hsml;
  in->NormSph = StarParticle[i].NormSph;
  in->EarlyTotalEnergyReleased = StarParticle[i].EarlyTotalEnergyReleased;
  in->NumNgb = StarParticle[i].NumNgb;

  in->Firstnode = firstnode;
}

static void out2particle(data_out * out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES)
    {
      kicked_cells += out->kicked_cells;
      massToKick += out->mass_to_kick;
      massKicked += out->mass_kicked;
      minKickVel = dmin(minKickVel, out->kick_vel);
      maxKickVel = dmax(maxKickVel, out->kick_vel);
      StarParticle[i].deltaEarlyEKin = out->deltaEarlyEKin;
      StarParticle[i].deltaEarlyMomentum[0] = out->deltaEarlyMomentum[0];
      StarParticle[i].deltaEarlyMomentum[1] = out->deltaEarlyMomentum[1];
      StarParticle[i].deltaEarlyMomentum[2] = out->deltaEarlyMomentum[2];
      StarParticle[i].EarlyMomentumInjected = out->EarlyMomentumInjected;
    }
  else
    {
      kicked_cells += out->kicked_cells;
      massKicked += out->mass_kicked;
      minKickVel = dmin(minKickVel, out->kick_vel);
      maxKickVel = dmax(maxKickVel, out->kick_vel);
      StarParticle[i].deltaEarlyEKin += out->deltaEarlyEKin;
      StarParticle[i].deltaEarlyMomentum[0] += out->deltaEarlyMomentum[0];
      StarParticle[i].deltaEarlyMomentum[1] += out->deltaEarlyMomentum[1];
      StarParticle[i].deltaEarlyMomentum[2] += out->deltaEarlyMomentum[2];
      StarParticle[i].EarlyMomentumInjected += out->EarlyMomentumInjected;
    }
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

        if(is_doing_early_feedback(i))
          early_stellar_feedback_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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

        early_stellar_feedback_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}


void do_early_stellar_feedback(void)
{
  long long ntot;

  sumup_large_ints(1, &Nstar, &ntot);
  if(ntot == 0)
    return;

  TIMER_STOPSTART(CPU_GFM_ENRICH, CPU_GFM_FEEDBACK);

  kicked_cells = 0;
  massToKick = massKicked = 0.0;

  if(Nstar > 0)
    {
      maxKickVel = -MAX_REAL_NUMBER;
      minKickVel = MAX_REAL_NUMBER;

      /* it can happen that a star is active but no SN exploded */
      for(int i = 0; i < Nstar; i++)
        if(StarParticle[i].EarlyTotalEnergyReleased <= 0)
          minKickVel = maxKickVel = 0.0;
    }
  else
    maxKickVel = minKickVel = 0.0;

  mpi_printf("FM_EARLY_FEEDBACK: Begin early stellar feedback calculation.\n");

  generic_set_MaxNexport();

  double t0 = second();

  generic_comm_pattern(Nstar, kernel_local, kernel_imported);

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
#ifdef OUTPUT_EARLY_STELLAR_FEEDBACK
           STP(StarParticle[i].index).EarlyFeedbackMomentum = StarParticle[i].EarlyMomentumInjected;
           STP(StarParticle[i].index).Cum_EarlyFeedbackMomentum += StarParticle[i].EarlyMomentumInjected;
#endif
          if(StarParticle[i].EarlyTotalEnergyReleased > 0)
            {
              energy_released = All.EarlyEtaKineticEnergy * StarParticle[i].EarlyTotalEnergyReleased * All.cf_atime * All.cf_atime;
              energy_diff = (StarParticle[i].deltaEarlyEKin - energy_released) / energy_released;
              momentum_diff = sqrt(StarParticle[i].deltaEarlyMomentum[0] * StarParticle[i].deltaEarlyMomentum[0] +
                                   StarParticle[i].deltaEarlyMomentum[1] * StarParticle[i].deltaEarlyMomentum[1] + StarParticle[i].deltaEarlyMomentum[2] * StarParticle[i].deltaEarlyMomentum[2]);
              active_stars++;
            }
          else
            {
              energy_diff = 0.0;        // because no radiation from young stars has occurred
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
  MPI_Reduce(&minKickVel, &global_minKickVel, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&maxKickVel, &global_maxKickVel, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
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
    ("FM_EARLY_FEEDBACK: kinetic feedback (this step)  ---> number of kicked cells = %d, active stars (with SN explosions) = %d, kicked cells / star = %g\n",
     sum_kicked_cells, tot_active_stars, kicked_cells_per_star);
  mpi_printf("FM_EARLY_FEEDBACK: kinetic feedback (this step)  ---> number of cells to kick = %g, cells to kick / star = %g\n", cells_to_kick, cells_to_kick / tot_active_stars);
  mpi_printf("FM_EARLY_FEEDBACK: kinetic feedback (this step)  ---> total mass to kick = %g, total kicked mass = %g, avg mass kicked = %g\n", sumMassToKick, sumMassKicked, avg_mass_kicked);
  mpi_printf("                                                 target gas mass = %g\n", All.TargetGasMass);
  mpi_printf("FM_EARLY_FEEDBACK: kinetic feedback (cumulative) ---> total mass to kick = %g, total kicked mass = %g\n", totalMassToKick, totalMassKicked);
  mpi_printf("FM_EARLY_FEEDBACK: kinetic feedback (this step)  ---> max kick velocity = %g, min kick velocity = %g\n", global_maxKickVel, global_minKickVel);
  mpi_printf("FM_EARLY_FEEDBACK: kinetic feedback (this step)  ---> max energy diff = %g, min energy diff = %g\n", global_max_energy_diff, global_min_energy_diff);
  mpi_printf("FM_EARLY_FEEDBACK: kinetic feedback (this step)  ---> max momentum diff = %g, min momentum diff = %g\n", global_max_momentum_diff, global_min_momentum_diff);

  double t1 = second();

  mpi_printf("FM_EARLY_FEEDBACK: stellar feedback calculation done. Took % g sec\n", timediff(t0, t1));

  TIMER_STOPSTART(CPU_GFM_FEEDBACK, CPU_GFM_ENRICH);
}


int early_stellar_feedback_evaluate(int target, int mode, int thread_id)
{
  int numnodes, *firstnode;
  data_in local, *in;
  data_out out;

  double h, h2;
#ifndef GFM_TOPHAT_KERNEL
  double wk, u, hinv, hinv3;
#endif
#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif
  double weight_fac;
  double dx, dy, dz, r2, r;
  MyDouble *pos;
  MyDouble de_feedback;
  MyFloat normsph;

  MyFloat EKin, prob, velocity, tot_mass;
  MyDouble EarlyTotalEnergyReleased, EarlyTotalKinEnergyReleased;

  int n_cells_kicked = 0;
  MyDouble mass_kicked = 0.0;
  MyDouble sum_deltaEarlyEKin = 0.0;
  MyDouble sum_deltaEarlyMomentum[3] = {0.0, 0.0, 0.0};

  MyDouble dp_feedback, inj_mom[3], dp_tot = 0;

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

  pos = in->Pos;
  h = in->Hsml;
  normsph = in->NormSph;
  EarlyTotalEnergyReleased = in->EarlyTotalEnergyReleased * All.cf_atime * All.cf_atime;
  EarlyTotalKinEnergyReleased = All.EarlyEtaKineticEnergy * EarlyTotalEnergyReleased;
  tot_mass = in->NumNgb * All.ReferenceGasPartMass;

  //star_sph_radius = pow(0.75 * normsph / 3.14159,0.3333); 
  //velocity = pow(2. * GRAVITY * tot_mass / star_sph_radius,0.5);
  velocity = sqrt(All.G * tot_mass / h);
  velocity *= 2. * All.cf_atime;    /* factor of 2 to be on the safe side. LVS: check cosmology factor?? */

  //printf("FM_EARLY_VELOCITY0 vel=%g, tot_mass=%g, h_radius=%g \n",velocity, tot_mass, h);

  /* couldn't find any neighbour (it can happen in zoom runs in the low-low res region) */
  if(tot_mass <= 0)
    prob = 0.0;
  else
    {
      prob = 2.0 * EarlyTotalKinEnergyReleased / (velocity * velocity * tot_mass);  /*LVS: here velocity CANNOT be larger than set in parameterfile */
      velocity *= dmax(1.0, sqrt(prob) + 1.0e-6);   /*LVS: here velocity CAN be larger than set in parameterfile $$ */
      prob = 2.0 * EarlyTotalKinEnergyReleased / (velocity * velocity * tot_mass);  /*LVS: here velocity CAN be larger than set in parameter file $$ */
    }

#ifndef GFM_TOPHAT_KERNEL
  hinv = 1.0 / h;
#ifndef  TWODIMS
  hinv3 = hinv * hinv * hinv;
#else
  hinv3 = hinv * hinv / boxSize_Z;
#endif
#endif

  h2 = h * h;

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, thread_id, numnodes, firstnode);

  for(int n = 0; n < nfound; n++)
    {
      int j = Thread[thread_id].Ngblist[n];

      if(P[j].Mass > 0 && P[j].ID != 0) /* skip cells that have been swallowed or dissolved */
        {
          dx = NEAREST_X(pos[0] - P[j].Pos[0]);
          dy = NEAREST_Y(pos[1] - P[j].Pos[1]);
          dz = NEAREST_Z(pos[2] - P[j].Pos[2]);

          r2 = dx * dx + dy * dy + dz * dz;

          if(r2 < h2)
            {
              r = sqrt(r2);
#ifndef GFM_TOPHAT_KERNEL
              u = r * hinv;

              if(u < 0.5)
                wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
              else
                wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

              weight_fac = SphP[j].Volume * wk / normsph;
#else
              weight_fac = SphP[j].Volume / normsph;
#endif
              if(EarlyTotalEnergyReleased > 0)
                {
                  de_feedback = 0.0;
                  EKin = 0.0;

                  /* --- activate this if want to use all energy available, activate together with " $$ " above ---
                     if(prob > 1.0)
                     terminate("Early Feedback, Total mass within kernel is less than mass to be kicked: tot_mass %g, kick_mass %g",
                     tot_mass, prob * tot_mass); */
                  /* -- end activate $$ ---- */

                  /* do feedback */
                  if(get_random_number() < prob)
                    {
                      /* injected feedback momentum */
                      dp_feedback = P[j].Mass * velocity;
                      dp_tot += dp_feedback;
#if (FM_EARLY_STAR_FEEDBACK_KICK_TYPE == 0)
                      double theta = acos(2 * get_random_number() - 1);
                      double phi = 2 * M_PI * get_random_number();
                      inj_mom[0] = -dp_feedback * sin(theta) * cos(phi);
                      inj_mom[1] = -dp_feedback * sin(theta) * sin(phi);
                      inj_mom[2] = -dp_feedback * cos(theta);
#endif
#if (FM_EARLY_STAR_FEEDBACK_KICK_TYPE == 1)
                      inj_mom[0] = -dp_feedback * dx / r;
                      inj_mom[1] = -dp_feedback * dy / r;
                      inj_mom[2] = -dp_feedback * dz / r;
#endif
#if (FM_EARLY_STAR_FEEDBACK_KICK_TYPE == 2)
                      double dir[3];
                      double norm;
                      dir[0] = P[j].GravAccel[1] * SphP[j].Momentum[2] - P[j].GravAccel[2] * SphP[j].Momentum[1];
                      dir[1] = P[j].GravAccel[2] * SphP[j].Momentum[0] - P[j].GravAccel[0] * SphP[j].Momentum[2];
                      dir[2] = P[j].GravAccel[0] * SphP[j].Momentum[1] - P[j].GravAccel[1] * SphP[j].Momentum[0];
                      norm = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);

                      if(get_random_number() < 0.5)
                        norm = -norm;

                      /* what happens if norm is zero ??? */
                      inj_mom[0] = -dp_feedback * dir[0] / norm;
                      inj_mom[1] = -dp_feedback * dir[1] / norm;
                      inj_mom[2] = -dp_feedback * dir[2] / norm;
#endif
                      /* Compute internal energy */
                      EKin = 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] + SphP[j].Momentum[2] * SphP[j].Momentum[2]) / P[j].Mass;
                      SphP[j].Energy -= EKin;

                      double deltaEarlyEKin = -EKin;

                      /* this accounts for the thermal energy from feedback */
                      de_feedback = 0.5 * dp_feedback * dp_feedback / (P[j].Mass * All.EarlyEtaKineticEnergy);
                      SphP[j].Energy += (1.0 - All.EarlyEtaKineticEnergy) * de_feedback;

                      /* momentum is injected in radial direction */
                      SphP[j].Momentum[0] += inj_mom[0];
                      SphP[j].Momentum[1] += inj_mom[1];
                      SphP[j].Momentum[2] += inj_mom[2];

                      /* Compute new kinetic energy and add to the internal energy */
                      EKin = 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] + SphP[j].Momentum[2] * SphP[j].Momentum[2]) / P[j].Mass;
                      SphP[j].Energy += EKin;

                      //SphP[j].IsKicked = 1;

                      deltaEarlyEKin += EKin;
                      sum_deltaEarlyEKin += deltaEarlyEKin;

                      sum_deltaEarlyMomentum[0] += inj_mom[0];
                      sum_deltaEarlyMomentum[1] += inj_mom[1];
                      sum_deltaEarlyMomentum[2] += inj_mom[2];

                      mass_kicked += P[j].Mass;
                      /* SphP[j].KickTimeEarlyFeedback = -All.Time; *//*LVS current time but negative as flag */

                      n_cells_kicked++;
                    }
                }
            }
        }
    }

  out.kicked_cells = n_cells_kicked;
  out.mass_to_kick = prob * tot_mass;
  out.mass_kicked = mass_kicked;
  out.kick_vel = velocity;
  out.deltaEarlyEKin = sum_deltaEarlyEKin;
  out.deltaEarlyMomentum[0] = sum_deltaEarlyMomentum[0];
  out.deltaEarlyMomentum[1] = sum_deltaEarlyMomentum[1];
  out.deltaEarlyMomentum[2] = sum_deltaEarlyMomentum[2];
  out.EarlyMomentumInjected = dp_tot;

  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#endif
