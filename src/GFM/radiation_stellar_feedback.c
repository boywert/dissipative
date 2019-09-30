/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/GFM/radiation_stellar_feedback.c
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

#if defined(FM_RADIATION_FEEDBACK)

typedef struct
{
  MyDouble Pos[3];
  MyFloat NormSphRadFeedback;
  MyDouble RadiationMomentumReleased;
  MyFloat StromgrenRadius;
  MyFloat RadCoolShutoffTime;
  MyFloat Velocity;
#ifdef FM_STOCHASTIC_HII_PHOTOIONIZATION
  MyFloat NormSphRadFeedback_cold;
  MyFloat StromgrenMass;
  MyFloat Hsml;
  MyFloat LowestDensityDirection_x;
  MyFloat LowestDensityDirection_y;
  MyFloat LowestDensityDirection_z;
#endif
  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

typedef struct
{
  int kicked_cells;
  MyFloat StromgrenRadius;
  MyFloat ShutoffTime;
  MyDouble mass_kicked;
  MyFloat kick_vel;
  MyDouble deltaEKin;
  MyDouble deltaMomentum[3];
} data_out;

static data_out *DataResult, *DataOut;

static int sum_kicked_cells, kicked_cells;
static float max_strom_radius, min_strom_radius;
static float global_max_strom_radius, global_min_strom_radius;
static float max_shutoff = 0.0, min_shutoff = 0.0;
static float global_max_shutoff, global_min_shutoff;
static double massKicked, sumMassKicked;
static float maxKickVel, minKickVel, global_maxKickVel, global_minKickVel;
static double global_max_energy_diff, global_min_energy_diff;
static double global_max_momentum_diff, global_min_momentum_diff;

static void particle2in(data_in * in, int i, int firstnode)
{
  in->Pos[0] = P[StarParticle[i].index].Pos[0];
  in->Pos[1] = P[StarParticle[i].index].Pos[1];
  in->Pos[2] = P[StarParticle[i].index].Pos[2];

  in->NormSphRadFeedback = StarParticle[i].NormSphRadFeedback;
#ifdef FM_STOCHASTIC_HII_PHOTOIONIZATION
  in->NormSphRadFeedback_cold = StarParticle[i].NormSphRadFeedback_cold;
#endif
  in->RadiationMomentumReleased = StarParticle[i].RadiationMomentumReleased;
  in->StromgrenRadius = StarParticle[i].StromgrenRadius;
  in->RadCoolShutoffTime = StarParticle[i].RadCoolShutoffTime;

  double RadiationMomentumReleased = StarParticle[i].RadiationMomentumReleased * All.cf_atime;
#ifdef FM_STOCHASTIC_HII_PHOTOIONIZATION
  in->Velocity = compute_radiation_feedback_velocity(i, StarParticle[i].NormSphRadFeedback, STP(StarParticle[i].index).Hsml * All.cf_atime, RadiationMomentumReleased);
  in->StromgrenMass = StarParticle[i].StromgrenMass;
  in->Hsml = STP(StarParticle[i].index).Hsml;
  in->LowestDensityDirection_x = StarParticle[i].LowestDensityDirection_x;
  in->LowestDensityDirection_y = StarParticle[i].LowestDensityDirection_y;
  in->LowestDensityDirection_z = StarParticle[i].LowestDensityDirection_z;
#else
  in->Velocity = compute_radiation_feedback_velocity(i, StarParticle[i].NormSphRadFeedback, StarParticle[i].StromgrenRadius * All.cf_atime, RadiationMomentumReleased);
#endif

  in->Firstnode = firstnode;
}

static void out2particle(data_out * out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES)
    {
      kicked_cells += out->kicked_cells;
      max_strom_radius = dmax(max_strom_radius, out->StromgrenRadius);
      min_strom_radius = dmin(min_strom_radius, out->StromgrenRadius);
      massKicked += out->mass_kicked;
      minKickVel = dmin(minKickVel, out->kick_vel);
      maxKickVel = dmax(maxKickVel, out->kick_vel);
      StarParticle[i].deltaEKin = out->deltaEKin;
      StarParticle[i].deltaMomentum[0] = out->deltaMomentum[0];
      StarParticle[i].deltaMomentum[1] = out->deltaMomentum[1];
      StarParticle[i].deltaMomentum[2] = out->deltaMomentum[2];
    }
  else
    {
      kicked_cells += out->kicked_cells;
      max_strom_radius = dmax(max_strom_radius, out->StromgrenRadius);
      min_strom_radius = dmin(min_strom_radius, out->StromgrenRadius);
      massKicked += out->mass_kicked;
      minKickVel = dmin(minKickVel, out->kick_vel);
      maxKickVel = dmax(maxKickVel, out->kick_vel);
      StarParticle[i].deltaEKin += out->deltaEKin;
      StarParticle[i].deltaMomentum[0] += out->deltaMomentum[0];
      StarParticle[i].deltaMomentum[1] += out->deltaMomentum[1];
      StarParticle[i].deltaMomentum[2] += out->deltaMomentum[2];
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

#ifdef FM_STOCHASTIC_HII_PHOTOIONIZATION
        if((StarParticle[i].StromgrenMass > 0.0) && (StarParticle[i].RadiationMomentumReleased > 0.0) && (StarParticle[i].NormSphRadFeedback > 0) )
          radiation_feedback_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
#else
        if(StarParticle[i].StromgrenRadius > 0.0)
          radiation_feedback_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
#endif
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

        radiation_feedback_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}


#ifdef FM_STOCHASTIC_HII_PHOTOIONIZATION

void do_radiation_stellar_feedback(void)
{
  long long ntot;
  sumup_large_ints(1, &Nstar, &ntot);

  if(ntot == 0)  return;

  mpi_printf("RADIATION_FEEDBACK: Begin radiation stellar feedback calculation.\n");

  generic_set_MaxNexport();
  double t0 = second();
  generic_comm_pattern(Nstar, kernel_local, kernel_imported);
  double t1 = second();

  mpi_printf("FM_RADIATION_FEEDBACK: stellar feedback calculation took %g sec\n", timediff(t0, t1));
}






int radiation_feedback_evaluate(int target, int mode, int thread_id)
{
  int numnodes, *firstnode;
  data_in local, *in;
  data_out out;

  double h, h2, hinv, hinv3;
  double dx, dy, dz, r2, r;
  MyDouble *pos;
  MyDouble de_feedback;

  MyDouble RadiationMomentumReleased;
  MyFloat RadCoolShutoffTime;
  int n_cells_kicked = 0;
  MyFloat velocity, tot_mass;
  MyDouble mass_kicked, deltaEKin, sum_deltaEKin, sum_deltaMomentum[3];
  MyDouble dp_feedback, inj_mom[3];

  MyFloat meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC)) * PROTONMASS;      /* note: assuming FULL ionization */
  MyFloat u_to_temp = meanweight / BOLTZMANN * (GAMMA_MINUS1) * All.UnitEnergy_in_cgs / All.UnitMass_in_g;

  MyFloat u, temp, unew, du;
  MyFloat uphys, rho, ne;

#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif


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
  if(h > 0.1) h = 0.1;                  // ToDo:  At large distances the ionizing photon input rate needs to account for the geometric incidence of radiation 
  MyFloat strom_mass = in->StromgrenMass;
  MyFloat rp_kernel_mass = in->NormSphRadFeedback;
  MyFloat pi_kernel_mass = in->NormSphRadFeedback_cold;
  RadiationMomentumReleased = in->RadiationMomentumReleased * All.cf_atime; /* copied from AGBMomentum in GFM/stellar_feedbac.c */
  RadCoolShutoffTime = in->RadCoolShutoffTime;
  velocity = in->Velocity;

  MyFloat rad_dx = in->LowestDensityDirection_x;
  MyFloat rad_dy = in->LowestDensityDirection_y;
  MyFloat rad_dz = in->LowestDensityDirection_z;

  MyFloat rad_dr = sqrt( rad_dx * rad_dx + rad_dy * rad_dy + rad_dz * rad_dz );
  rad_dx /= rad_dr;
  rad_dy /= rad_dr;
  rad_dz /= rad_dr; 

  if(in->NormSphRadFeedback <= 0 && h > 0)
        terminate("tot_mass<=0 should not happen in radiation_stellar_feedback mode 0: pos %g|%g|%g, tot_mass %g, NumNgb %g, h %g , RadiationMomentumReleased %g, velocity %g,  ID = %d\n",
                  pos[0], pos[1], pos[2], tot_mass, StarParticle[target].NumNgb, h, RadiationMomentumReleased, velocity, P[StarParticle[target].index].ID);

  h2 = h * h;
  hinv = 1.0 / h;
#ifndef  TWODIMS
  hinv3 = hinv * hinv * hinv;
#else
  hinv3 = hinv * hinv / boxSize_Z;
#endif

  mass_kicked = 0.0;
  sum_deltaEKin = 0.0;
  sum_deltaMomentum[0] = 0.0;
  sum_deltaMomentum[1] = 0.0;
  sum_deltaMomentum[2] = 0.0;

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
              de_feedback = 0.;
              dp_feedback = 0.;
              /* First input momentum if applicable */
              if(RadiationMomentumReleased > 0)
                {
                  double Ekin = 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] + SphP[j].Momentum[2] * SphP[j].Momentum[2]) / P[j].Mass;
                  dp_feedback = P[j].Mass * velocity;   /* all particles within rs same velocity */

                  /* injected feedback momentum radially away */
                  inj_mom[0] = -dp_feedback * rad_dx;	//dx / r;
                  inj_mom[1] = -dp_feedback * rad_dy;	//dy / r;
                  inj_mom[2] = -dp_feedback * rad_dz;	//dz / r;

                  /* momentum due to stellar mass return is injected in radial direction */
                  SphP[j].Momentum[0] += inj_mom[0] ;
                  SphP[j].Momentum[1] += inj_mom[1] ;
                  SphP[j].Momentum[2] += inj_mom[2] ;

                  de_feedback = 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] + SphP[j].Momentum[2] * SphP[j].Momentum[2]) / P[j].Mass - Ekin;

                  deltaEKin += de_feedback;
                  sum_deltaEKin += deltaEKin;

                  sum_deltaMomentum[0] += inj_mom[0];
                  sum_deltaMomentum[1] += inj_mom[1];
                  sum_deltaMomentum[2] += inj_mom[2];

                  mass_kicked += P[j].Mass;
                  n_cells_kicked++;
                }               /* if RadMom > 0 */

              double p_ionize = strom_mass / pi_kernel_mass;

              u = SphP[j].Utherm;
              uphys = u * All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;
              rho = SphP[j].Density * All.cf_a3inv * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;        /* convert to physical cgs units */
              ne = SphP[j].Ne;
              temp = convert_u_to_temp(uphys, rho, &ne);

              /* identify gas that could be photoionized, but make sure we can "reselect" gas that has just returned to non-photoionized state */
              if(SphP[j].Utherm < 2.0*All.PhotoionizationEgySpec && SphP[j].GasRadCoolShutoffTime==0.0)
              {
                if(get_random_number() < p_ionize)
                {

                  du = dmax(All.PhotoionizationEgySpec - u, 0.0);
                  de_feedback += du * P[j].Mass;        /* add change in temperature to energy */
                  SphP[j].Utherm += du;
                  SphP[j].Ne = 1.;
                  SphP[j].GasRadCoolShutoffTime = RadCoolShutoffTime;
                  
                  //printf("   PHOTOIONIZATION WORKS!  strom_mass = %f  pi_kernel_mass = %f \n", strom_mass, pi_kernel_mass);
               }
             }
             SphP[j].Energy += de_feedback;    /* de_feedback has a momentum and heat component */
            }                   /* if r2<rs2) */
        }                       /* if P[j].ID != 0 */
    }

  out.kicked_cells = n_cells_kicked;
  out.StromgrenRadius = h;
  out.ShutoffTime = RadCoolShutoffTime;
  out.mass_kicked = mass_kicked;
  out.kick_vel = velocity;
  out.deltaEKin = sum_deltaEKin;
  out.deltaMomentum[0] = sum_deltaMomentum[0];
  out.deltaMomentum[1] = sum_deltaMomentum[1];
  out.deltaMomentum[2] = sum_deltaMomentum[2];

  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}


#else

void do_radiation_stellar_feedback(void)
{
  long long ntot;
  sumup_large_ints(1, &Nstar, &ntot);

  if(ntot == 0)
    return;

  int i;

  kicked_cells = 0;

  if(Nstar > 0)
    {
      max_strom_radius = -MAX_REAL_NUMBER;
      min_strom_radius = MAX_REAL_NUMBER;

      /* it can happen that a star is active but no radiation is emitted */
      for(i = 0; i < Nstar; i++)
        if(StarParticle[i].RadiationMomentumReleased <= 0)
          max_strom_radius = min_strom_radius = 0.0;
    }
  else
    max_strom_radius = min_strom_radius = 0.0;

  kicked_cells = 0;
  massKicked = 0.0;

  if(Nstar > 0)
    {
      maxKickVel = -MAX_REAL_NUMBER;
      minKickVel = MAX_REAL_NUMBER;

      /* it can happen that a star is active but no radiation is emitted */
      for(i = 0; i < Nstar; i++)
        if(StarParticle[i].RadiationMomentumReleased <= 0)
          minKickVel = maxKickVel = 0.0;
    }
  else
    maxKickVel = minKickVel = 0.0;


  mpi_printf("RADIATION_FEEDBACK: Begin radiation stellar feedback calculation.\n");

  generic_set_MaxNexport();


  double t0 = second();

  generic_comm_pattern(Nstar, kernel_local, kernel_imported);


  int kicked_cells_per_star = 0;
  int active_stars = 0;
  int tot_active_stars = 0;
  int no_cool_cells = 0;
  int tot_no_cool_cells = 0;

  /* count only active star particles that deposited radiation */
  for(i = 0; i < Nstar; i++)
    if(StarParticle[i].RadiationMomentumReleased > 0)
      active_stars++;

  if(active_stars > 0)
    {
      max_shutoff = -MAX_REAL_NUMBER;
      min_shutoff = MAX_REAL_NUMBER;
    }

  /* count gas particles with deactivated cooling */
  for(i = 0; i < NumGas; i++)
    if(SphP[i].GasRadCoolShutoffTime > 0)
      {
        max_shutoff = dmax(max_shutoff, SphP[i].GasRadCoolShutoffTime);
        min_shutoff = dmin(min_shutoff, SphP[i].GasRadCoolShutoffTime);
        no_cool_cells++;
      }

  MPI_Reduce(&kicked_cells, &sum_kicked_cells, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&active_stars, &tot_active_stars, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&no_cool_cells, &tot_no_cool_cells, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&min_strom_radius, &global_min_strom_radius, 1, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&max_strom_radius, &global_max_strom_radius, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&min_shutoff, &global_min_shutoff, 1, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&max_shutoff, &global_max_shutoff, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);

  if(tot_active_stars > 0)
    kicked_cells_per_star = ((float) sum_kicked_cells / tot_active_stars);

  mpi_printf("FM_RADIATION_FEEDBACK:  ---> number of kicked cells = %d, active stars (with SN explosions) = %d, kicked cells / star = %g\n", sum_kicked_cells, tot_active_stars, kicked_cells_per_star);
  mpi_printf("FM_RADIATION_FEEDBACK:  ---> number of cells with deactivated cooling = %d\n", tot_no_cool_cells);
  mpi_printf("FM_RADIATION_FEEDBACK:  ---> min stromgren radius (int units)= %g, max stromgren radius (int units)= %g\n", global_min_strom_radius, global_max_strom_radius);
  mpi_printf("FM_RADIATION_FEEDBACK:  ---> as computed in SphP struct  min shutoff time (int units)= %g, max shutoff time (int units)= %g\n", global_min_shutoff, global_max_shutoff);



  double energy_diff, energy_released, max_energy_diff, min_energy_diff;
  double momentum_diff, max_momentum_diff, min_momentum_diff;
  float avg_mass_kicked = 0.0;
  if(Nstar > 0)
    {
      max_energy_diff = -MAX_REAL_NUMBER;
      min_energy_diff = MAX_REAL_NUMBER;
      max_momentum_diff = -MAX_REAL_NUMBER;
      min_momentum_diff = MAX_REAL_NUMBER;

      /* count only active star particles that emitted radiation */
      for(i = 0; i < Nstar; i++)
        {
          if(StarParticle[i].RadiationMomentumReleased > 0)
            {
              energy_released = StarParticle[i].RadiationMomentumReleased * All.cf_atime * All.cf_atime;
              energy_diff = (StarParticle[i].deltaEKin - energy_released) / energy_released;
              momentum_diff = sqrt(StarParticle[i].deltaMomentum[0] * StarParticle[i].deltaMomentum[0] +
                                   StarParticle[i].deltaMomentum[1] * StarParticle[i].deltaMomentum[1] + StarParticle[i].deltaMomentum[2] * StarParticle[i].deltaMomentum[2]);
            }
          else
            {
              energy_diff = 0.0;        // because no radiation feedback had occurred
              momentum_diff = 0.0;
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

  MPI_Reduce(&massKicked, &sumMassKicked, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&minKickVel, &global_minKickVel, 1, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&maxKickVel, &global_maxKickVel, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&max_energy_diff, &global_max_energy_diff, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&min_energy_diff, &global_min_energy_diff, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&max_momentum_diff, &global_max_momentum_diff, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&min_momentum_diff, &global_min_momentum_diff, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

  if(sum_kicked_cells > 0)
    avg_mass_kicked = sumMassKicked / sum_kicked_cells;

  if(tot_active_stars > 0)
    kicked_cells_per_star = ((float) sum_kicked_cells / tot_active_stars);

  mpi_printf("FM_RADIATION_FEEDBACK:  ---> total kicked mass = %g, avg mass kicked = %g\n", sumMassKicked, avg_mass_kicked);
  mpi_printf("FM_RADIATION_FEEDBACK: ---> max kick velocity = %g, min kick velocity = %g\n", global_maxKickVel, global_minKickVel);
  mpi_printf("FM_RADIATION_FEEDBACK: ---> max energy diff = %g, min energy diff = %g\n", global_max_energy_diff, global_min_energy_diff);
  mpi_printf("FM_RADIATION_FEEDBACK:  ---> max momentum diff = %g, min momentum diff = %g\n", global_max_momentum_diff, global_min_momentum_diff);

  double t1 = second();

  mpi_printf("FM_RADIATION_FEEDBACK: stellar feedback calculation took %g sec\n", timediff(t0, t1));
}


int radiation_feedback_evaluate(int target, int mode, int thread_id)
{
  int numnodes, *firstnode;
  data_in local, *in;
  data_out out;

  double h, h2, hinv, hinv3;
  double dx, dy, dz, r2, r;
  MyDouble *pos;
  MyDouble de_feedback;

  MyDouble RadiationMomentumReleased;
  MyFloat RadCoolShutoffTime;
  int n_cells_kicked = 0;
  MyFloat velocity, tot_mass;
  MyDouble mass_kicked, deltaEKin, sum_deltaEKin, sum_deltaMomentum[3];
  MyDouble dp_feedback, inj_mom[3];

  MyFloat meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC)) * PROTONMASS;      /* note: assuming FULL ionization */
  MyFloat u_to_temp = meanweight / BOLTZMANN * (GAMMA_MINUS1) * All.UnitEnergy_in_cgs / All.UnitMass_in_g;

  ///MyFloat yhelium, meanweight;
  //MyFloat u_to_temp =  (GAMMA_MINUS1) / BOLTZMANN * All.UnitEnergy_in_cgs / All.UnitMass_in_g;
  MyFloat u, temp, unew, du;
  MyFloat uphys, rho, ne;

#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif


  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      in = &local;

      numnodes = 1;
      firstnode = NULL;

      pos = in->Pos;
      h = in->StromgrenRadius;
      RadiationMomentumReleased = in->RadiationMomentumReleased * All.cf_atime; /* copied from AGBMomentum in GFM/stellar_feedbac.c */
      RadCoolShutoffTime = in->RadCoolShutoffTime;
      velocity = in->Velocity;

      if(in->NormSphRadFeedback <= 0 && h > 0)
        terminate("tot_mass<=0 should not happen in radiation_stellar_feedback mode 0: pos %g|%g|%g, tot_mass %g, NumNgb %g, h %g , RadiationMomentumReleased %g, velocity %g,  ID = %d\n",
                  pos[0], pos[1], pos[2], tot_mass, StarParticle[target].NumNgb, h, RadiationMomentumReleased, velocity, P[StarParticle[target].index].ID);
    }
  else
    {
      in = &DataGet[target];
      generic_get_numnodes(target, &numnodes, &firstnode);

      pos = in->Pos;
      h = in->StromgrenRadius;
      RadiationMomentumReleased = in->RadiationMomentumReleased * All.cf_atime;
      RadCoolShutoffTime = in->RadCoolShutoffTime;
      velocity = in->Velocity;
    }


  h2 = h * h;

  hinv = 1.0 / h;
#ifndef  TWODIMS
  hinv3 = hinv * hinv * hinv;
#else
  hinv3 = hinv * hinv / boxSize_Z;
#endif

  mass_kicked = 0.0;
  sum_deltaEKin = 0.0;
  sum_deltaMomentum[0] = 0.0;
  sum_deltaMomentum[1] = 0.0;
  sum_deltaMomentum[2] = 0.0;


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
              de_feedback = 0.;
              dp_feedback = 0.;
              /* First input momentum if applicable */
              if(RadiationMomentumReleased > 0)
                {
                  double Ekin = 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] + SphP[j].Momentum[2] * SphP[j].Momentum[2]) / P[j].Mass;

                  if(r<1.0)  // TODO:  Handle diffuse gas (optically thin) limit
                      dp_feedback = P[j].Mass * velocity;   /*all particles within rs same velocity */
                  else
                      dp_feedback = 0.0;

                  /* injected feedback momentum radially away */
                  inj_mom[0] = -dp_feedback * dx / r;
                  inj_mom[1] = -dp_feedback * dy / r;
                  inj_mom[2] = -dp_feedback * dz / r;

                  /* momentum due to stellar mass return is injected in radial direction */
                  SphP[j].Momentum[0] += inj_mom[0];
                  SphP[j].Momentum[1] += inj_mom[1];
                  SphP[j].Momentum[2] += inj_mom[2];

                  de_feedback = 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] + SphP[j].Momentum[2] * SphP[j].Momentum[2]) / P[j].Mass - Ekin;

                  deltaEKin += de_feedback;
                  sum_deltaEKin += deltaEKin;

                  sum_deltaMomentum[0] += inj_mom[0];
                  sum_deltaMomentum[1] += inj_mom[1];
                  sum_deltaMomentum[2] += inj_mom[2];

                  mass_kicked += P[j].Mass;
                  n_cells_kicked++;
                }               /* if RadMom > 0 */

              /* Second, keep cells at 10^4 if applicable */
              /* Actually the ionization state should be changed accordingly */
              ///yhelium = (1 - SphP[j].MetalsFraction[element_index_Hydrogen] - SphP[j].Metallicity) / (4. * SphP[j].MetalsFraction[element_index_Hydrogen]);
              ///meanweight = (1 + 4 * yhelium) / (1 + yhelium + SphP[j].Ne);
              u = SphP[j].Utherm;
              uphys = u * All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;
              rho = SphP[j].Density * All.cf_a3inv * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;        /* convert to physical cgs units */
              ne = SphP[j].Ne;
              temp = convert_u_to_temp(uphys, rho, &ne);
              if(temp < 1.e4)
                {
                  temp = 1.e4;
                  unew = temp / u_to_temp;
                  du = unew - u;
                  de_feedback += du * P[j].Mass;        /* add change in temperature to energy */
                  SphP[j].Utherm += du;
                  SphP[j].Ne = 1.;
                }
              SphP[j].GasRadCoolShutoffTime = RadCoolShutoffTime;

              SphP[j].Energy += de_feedback;    /* de_feedback has a momentum and heat component */
/*
#ifdef USE_ENTROPY_FOR_COLD_FLOWS
		      SphP[j].Utherm =
			(SphP[j].Energy -

				SphP[j].Momentum[2] * SphP[j].Momentum[2]) / P[j].Mass) / P[j].Mass / (All.cf_atime * All.cf_atime);
		      SphP[j].A = (GAMMA - 1.0) * SphP[j].Utherm / pow(SphP[j].Density * All.cf_a3inv, GAMMA - 1);
		      SphP[j].Entropy = log(SphP[j].A) * P[j].Mass;
#endif
*/
            }                   /* if r2<rs2) */

        }                       /* if P[j].ID != 0 */
    }

  out.kicked_cells = n_cells_kicked;
  out.StromgrenRadius = h;
  out.ShutoffTime = RadCoolShutoffTime;
  out.mass_kicked = mass_kicked;
  out.kick_vel = velocity;
  out.deltaEKin = sum_deltaEKin;
  out.deltaMomentum[0] = sum_deltaMomentum[0];
  out.deltaMomentum[1] = sum_deltaMomentum[1];
  out.deltaMomentum[2] = sum_deltaMomentum[2];

  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#endif

MyFloat compute_radiation_feedback_velocity(int n, MyFloat stromgren_mass, MyFloat stromgren_radius, MyFloat radiatmom)
{
  MyFloat mass_phy, rs_phy, sigma_gas, tau;

  mass_phy = stromgren_mass / All.HubbleParam * All.UnitMass_in_g;
  rs_phy = stromgren_radius / All.HubbleParam * All.UnitLength_in_cm;
  if(rs_phy > 0)  sigma_gas = mass_phy / (3.14159 * rs_phy * rs_phy);
  else  sigma_gas = 0;
  tau = All.DustOppacityRadiationFeedback * sigma_gas;  /*Dust oppacity = \kappa and in physical units */
  StarParticle[n].RadFeedTau = tau;
  //printf("LVS-TAU=%g \n");
  STP(StarParticle[n].index).RadFeedTau = StarParticle[n].RadFeedTau;

  MyFloat velocity = radiatmom / stromgren_mass * (1. + tau);   /*code units and comoving */

  return velocity;
}


#endif
