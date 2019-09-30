/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/GFM/winds.c
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

#include "../allvars.h"
#include "../proto.h"
#include "../voronoi.h"


#if defined(GFM_WINDS) || defined(GFM_WINDS_LOCAL)


/* local number of wind particles and maximum number that can be stored on task */
static int Nwind, MaxNwind;
/* counter for recoupled local wind particles */


/* initialize wind energy */
void init_winds(void)
{
#ifdef GFM_STELLAR_EVOLUTION
  All.WindEnergyFactor = All.WindEnergyIn1e51erg;
  /* All.WindEgySpecSN is calculated (if needed for GFM_CONST_IMF) in init_imf() */
#else
  All.WindEnergyFactor = All.WindEnergyFraction;
  All.WindEgySpecSN = All.EgySpecSN;
#endif
}

/* main driver routine to process wind particles */
void do_winds(void)
{
  TIMER_START(CPU_GFM_WINDS);

  double t0 = second();
  int idx, i;
  int local_recoupled;
  double local_mass_recoupled;

#ifdef TRACER_MC_CHECKS
  long long check_total_tracers = get_total_number_of_tracers(-1);
  long long check_total_bh_tracers_rear = get_total_number_of_tracers(5);
#endif


  /* initialize local recoupling wind particle counter */
  Nwind = 0;

  /* maximum number of local recoupling wind particles */
  MaxNwind = NumPart;

  /* allocate structure for local, recoupling, wind particles */
  WindParticle = mymalloc("WindParticle", MaxNwind * sizeof(struct wind_particle));

  /* set up wind particle data */
  setup_windparticles();

  /* find Voronoi cells wind particles */
  find_wind_cells(Nwind);

  /* cool wind and check which wind particles to recouple */
  cool_and_check_recouple();

  /* do the recoupling, i.e. put gas back in cell */
  recouple_wind_particles(Nwind, &local_recoupled, &local_mass_recoupled);

#ifdef VERBOSE
  int tot_recoupled;
  double tot_mass_recoupled;
  MPI_Allreduce(&local_recoupled, &tot_recoupled, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&local_mass_recoupled, &tot_mass_recoupled, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  mpi_printf("GFM_WINDS: wind particles recoupled=%d  wind mass recoupled=%g\n", tot_recoupled, tot_mass_recoupled);
#endif

  for(i = 0; i < Nwind; i++)
    if(STP(WindParticle[i].index).BirthTime == 0)
      {
        P[WindParticle[i].index].Mass = 0;
      }

  /* remove wind particle that recoupled now from timebins */
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Type != 4)
        continue;

      if(STP(i).BirthTime == 0)
        timebin_remove_particle(&TimeBinsGravity, idx, P[i].TimeBinGrav);
    }

  /* free wind particles */
  myfree(WindParticle);


  All.CumWindEnergy_Is += WindEnergy_Is;
  All.CumWindEnergy_Should += WindEnergy_Should;

#ifdef VERBOSE
  /* sum up total stellar mass for statistics */
  double local_stellar_mass = 0.0, global_stellar_mass;
  double global_SNII_energy;
  for(i = 0; i < NumPart; i++)
    if(P[i].Type == 4 && STP(i).BirthTime > 0)
      local_stellar_mass += P[i].Mass;

  MPI_Allreduce(&local_stellar_mass, &global_stellar_mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  /* total stellar mass in M_sun */
  global_stellar_mass *= (All.UnitMass_in_g / SOLAR_MASS) / All.HubbleParam;
  /* total SNII energy in erg */
  global_SNII_energy = WindEnergy_Is * All.UnitEnergy_in_cgs / All.HubbleParam;

  mpi_printf
    ("GFM_WINDS: energy statistics (step)       EnergyIs[internal units]=%e EnergyShould[internal units]=%e Ratio=%g EnergyIs/formed stellar mass[erg/M_sun]=%e\n",
     WindEnergy_Is, WindEnergy_Should, (WindEnergy_Should > 0) ? (WindEnergy_Is / WindEnergy_Should) : (0), (global_stellar_mass > 0) ? (global_SNII_energy / global_stellar_mass) : (0));

  mpi_printf("GFM_WINDS: energy statistics (cumulative) EnergyIs[internal units]=%e EnergyShould[internal units]=%e Ratio=%g\n",
             All.CumWindEnergy_Is, All.CumWindEnergy_Should, (All.CumWindEnergy_Should > 0) ? (All.CumWindEnergy_Is / All.CumWindEnergy_Should) : (0));
#endif

#ifdef TRACER_MC_CHECKS
  if(check_total_tracers != get_total_number_of_tracers(-1))
    terminate("TRACER_MC: unconserved tracers during GFM_WINDS: total number of tracers BEFORE = %lld, AFTER = %lld\n", check_total_tracers, get_total_number_of_tracers(-1));
  if(check_total_bh_tracers_rear != get_total_number_of_tracers(5))
    terminate("TRACER_MC: strange, global number of BH tracers changed in GFM_WINDS, BEFORE = %lld, AFTER = %lld\n", check_total_bh_tracers_rear, get_total_number_of_tracers(5));
#endif

  double t1 = second();
  mpi_printf("GFM_WINDS: done. (took %g sec)\n", timediff(t0, t1));

  TIMER_STOP(CPU_GFM_WINDS);
}


void setup_windparticles(void)
{
  int idx, i;
  double dt, dtime;

  /* find wind particle that recoupled now */
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Type != 4)
        continue;

      if(STP(i).BirthTime >= 0) /* real star particle or already recoupled wind particle */
        continue;

      dt = (P[i].TimeBinGrav ? (((integertime) 1) << P[i].TimeBinGrav) : 0) * All.Timebase_interval;

      dtime = dt / All.cf_hubble_a;

      /* increase time */
      if(STP(i).BirthTime < 0)
        STP(i).BirthTime += dtime;

      if(Nwind >= MaxNwind)
        terminate("Nwind >= MaxNwind");

      WindParticle[Nwind].CellID = 0;   /* no cell assigned yet, we require IDs>0 for this to work */
      WindParticle[Nwind].index = i;    /* index of wind particle */
      Nwind++;
    }

#ifdef VERBOSE
  int Nwind_tot = 0;
  MPI_Reduce(&Nwind, &Nwind_tot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  mpi_printf("GFM_WINDS: total number of active wind particles=%d\n", Nwind_tot);
#endif
}



void cool_and_check_recouple(void)
{
  int i;

  for(i = 0; i < Nwind; i++)
    {
      /* check for recoupling */
      if((WindParticle[i].Density * All.cf_a3inv < All.WindFreeTravelDensFac * All.PhysDensThresh) || (STP(WindParticle[i].index).BirthTime >= 0))
        STP(WindParticle[i].index).BirthTime = 0;

#if !defined(GFM_WINDS_THERMAL) && !defined(GFM_WINDS_THERMAL_NEWDEF)
      /* now do the cooling */
      double dt = (P[i].TimeBinGrav ? (((integertime) 1) << P[i].TimeBinGrav) : 0) * All.Timebase_interval;
      double dtime = dt / All.cf_hubble_a;
#ifdef GFM_COOLING_METAL
      update_gas_state(WindParticle[i].Density * All.cf_a3inv,
                       STP(WindParticle[i].index).MassMetals[element_index_Hydrogen] / P[WindParticle[i].index].Mass, STP(WindParticle[i].index).Metallicity);
#endif
      double ne = 1.0;
      double unew = DoCooling(dmax(All.MinEgySpec, STP(WindParticle[i].index).Utherm), WindParticle[i].Density * All.cf_a3inv, dtime, &ne);
      double du = unew - STP(WindParticle[i].index).Utherm;
      STP(WindParticle[i].index).Utherm += du;
#endif
    }
}





#ifdef GFM_WINDS
/* calculate the loading factor of the wind for a star-forming cell */
void gfm_calc_wind_parameters(int i, double p, double *p_wind, double *v_wind, double *u_wind)
{
  *v_wind = 0;
  *u_wind = 0;
  *p_wind = 0;

  if(P[i].Type != 0)            /* to protect using a particle that has been turned into a star */
    return;

  double stellar_mass = p * P[i].Mass;  /* amount of stars expect to form */

#ifdef GFM_WINDS_VARIABLE
#if (GFM_WINDS_VARIABLE==0)
  if(SphP[i].w.HostHaloMass > 0 && stellar_mass > 0)
    {
      wind_parameter wp;
      gfm_calc_variable_wind_parameters(stellar_mass, SphP[i].w.HostHaloMass, SphP[i].w.DMVelDisp, SphP[i].Metallicity, &wp);
      *p_wind = wp.wind_mass / P[i].Mass;
      *v_wind = wp.wind_velocity;
      *u_wind = wp.wind_utherm;
    }
#endif
#if (GFM_WINDS_VARIABLE==1)
  if(SphP[i].w.DMVelDisp > 0 && stellar_mass > 0)
    {
      wind_parameter wp;
      gfm_calc_variable_wind_parameters(stellar_mass, SphP[i].w.DMVelDisp, SphP[i].w.DMVelDisp, SphP[i].Metallicity, &wp);
      *p_wind = wp.wind_mass / P[i].Mass;
      *v_wind = wp.wind_velocity;
      *u_wind = wp.wind_utherm;
    }
#endif

#else

  double FactorSN, WindEgySpecSN;

#ifdef GFM_CONST_IMF
  FactorSN = All.FactorSN;
  WindEgySpecSN = All.WindEgySpecSN;
#elif defined(GFM_VARIABLE_IMF)

#if (GFM_VARIABLE_IMF == 0)
  define_imf(SphP[i].w.DMVelDisp);
#else
#error "GFM_VARIABLE_IMF mode is not ok"
#endif
  FactorSN = calc_FactorSN();
  WindEgySpecSN = calc_WindEgySpecSN();
#endif

  *p_wind = All.WindEfficiency * stellar_mass / P[i].Mass;

#if defined(GFM_STELLAR_EVOLUTION) && (GFM_STELLAR_EVOLUTION == 0)
  *v_wind = sqrt(2 * All.WindEnergyFactor * FactorSN * WindEgySpecSN / All.WindEfficiency);
#else
  *v_wind = sqrt(2 * All.WindEnergyFactor * FactorSN * WindEgySpecSN / (1 - FactorSN) / All.WindEfficiency);
#endif

#endif
}
#endif


/* convert a cell into a wind particle */
void gfm_add_wind(int i, double v, double utherm)
{
  double norm, dir[3];
  int j;
#ifndef GFM_BIPOLAR_WINDS
  MyFloat theta, phi;
#endif

#ifdef GFM_BIPOLAR_WINDS

#if !(GFM_BIPOLAR_WINDS == 3)
#ifdef HIERARCHICAL_GRAVITY
  MyFloat *GravAccel = SphP[i].FullGravAccel;
#else
  MyFloat *GravAccel = P[i].GravAccel;
#endif
#endif

#if (GFM_BIPOLAR_WINDS == 0)
  dir[0] = GravAccel[1] * (SphP[i].Momentum[2] / P[i].Mass) - GravAccel[2] * (SphP[i].Momentum[1] / P[i].Mass);
  dir[1] = GravAccel[2] * (SphP[i].Momentum[0] / P[i].Mass) - GravAccel[0] * (SphP[i].Momentum[2] / P[i].Mass);
  dir[2] = GravAccel[0] * (SphP[i].Momentum[1] / P[i].Mass) - GravAccel[1] * (SphP[i].Momentum[0] / P[i].Mass);
#else
#if (GFM_BIPOLAR_WINDS == 3)
  dir[0] = SphP[i].DensGasAngMomentum[0];
  dir[1] = SphP[i].DensGasAngMomentum[1];
  dir[2] = SphP[i].DensGasAngMomentum[2];
  if(dir[0] == 0 && dir[1] == 0 && dir[2] == 0)
    {
      double theta = acos(2 * get_random_number() - 1);
      double phi = 2 * M_PI * get_random_number();
      dir[0] = sin(theta) * cos(phi);
      dir[1] = sin(theta) * sin(phi);
      dir[2] = cos(theta);
    }
#else
  dir[0] =
    (GravAccel[1] - SphP[i].GroupGravAcc[1]) * ((SphP[i].Momentum[2] / P[i].Mass) - SphP[i].GroupVel[2]) - (GravAccel[2] -
                                                                                                            SphP[i].GroupGravAcc[2]) * ((SphP[i].Momentum[1] / P[i].Mass) - SphP[i].GroupVel[1]);
  dir[1] =
    (GravAccel[2] - SphP[i].GroupGravAcc[2]) * ((SphP[i].Momentum[0] / P[i].Mass) - SphP[i].GroupVel[0]) - (GravAccel[0] -
                                                                                                            SphP[i].GroupGravAcc[0]) * ((SphP[i].Momentum[2] / P[i].Mass) - SphP[i].GroupVel[2]);
  dir[2] =
    (GravAccel[0] - SphP[i].GroupGravAcc[0]) * ((SphP[i].Momentum[1] / P[i].Mass) - SphP[i].GroupVel[1]) - (GravAccel[1] -
                                                                                                            SphP[i].GroupGravAcc[1]) * ((SphP[i].Momentum[0] / P[i].Mass) - SphP[i].GroupVel[0]);
#endif
#endif
#else
  theta = acos(2 * get_random_number() - 1);
  phi = 2 * M_PI * get_random_number();
  dir[0] = sin(theta) * cos(phi);
  dir[1] = sin(theta) * sin(phi);
  dir[2] = cos(theta);
#endif

  for(j = 0, norm = 0; j < 3; j++)
    norm += dir[j] * dir[j];

  norm = sqrt(norm);
  if(get_random_number() < 0.5)
    norm = -norm;

  /* here we turn the gas particle itself into a wind particle, which is a special star particle */
  Stars_converted++;

  if(norm != 0)
    {
      for(j = 0; j < 3; j++)
        dir[j] /= norm;

      for(j = 0; j < 3; j++)
        {
          P[i].Vel[j] += (v * All.cf_atime * dir[j]);
        }
    }

  convert_cell_into_star(i, -All.WindFreeTravelMaxTimeFactor / All.cf_hubble_a);

  STP(i).Hsml = 1.5 * get_cell_radius(i);
  STP(i).Utherm = SphP[i].Utherm + utherm;

#ifdef REFINEMENT_HIGH_RES_GAS
  STP(i).HighResMass = SphP[i].HighResMass;
#endif
}

/* spawn a wind particle from a cell */
void gfm_spawn_wind_from_cell(int igas, double v, double utherm, int istar, MyFloat mass_of_wind)
{
  double norm, dir[3];
  int k;
#ifndef GFM_BIPOLAR_WINDS
  MyFloat theta, phi;
#endif

  /* calculate kick */
#ifdef GFM_BIPOLAR_WINDS

#if !(GFM_BIPOLAR_WINDS == 3)
#ifdef HIERARCHICAL_GRAVITY
  MyFloat *GravAccel = SphP[igas].FullGravAccel;
#else
  MyFloat *GravAccel = P[igas].GravAccel;
#endif
#endif

#if (GFM_BIPOLAR_WINDS == 0)
  dir[0] = GravAccel[1] * (SphP[igas].Momentum[2] / P[igas].Mass) - GravAccel[2] * (SphP[igas].Momentum[1] / P[igas].Mass);
  dir[1] = GravAccel[2] * (SphP[igas].Momentum[0] / P[igas].Mass) - GravAccel[0] * (SphP[igas].Momentum[2] / P[igas].Mass);
  dir[2] = GravAccel[0] * (SphP[igas].Momentum[1] / P[igas].Mass) - GravAccel[1] * (SphP[igas].Momentum[0] / P[igas].Mass);
#else
#if (GFM_BIPOLAR_WINDS == 3)
  dir[0] = SphP[igas].DensGasAngMomentum[0];
  dir[1] = SphP[igas].DensGasAngMomentum[1];
  dir[2] = SphP[igas].DensGasAngMomentum[2];
  if(dir[0] == 0 && dir[1] == 0 && dir[2] == 0)
    {
      double theta = acos(2 * get_random_number() - 1);
      double phi = 2 * M_PI * get_random_number();
      dir[0] = sin(theta) * cos(phi);
      dir[1] = sin(theta) * sin(phi);
      dir[2] = cos(theta);
    }
#else

  dir[0] =
    (GravAccel[1] - SphP[igas].GroupGravAcc[1]) * ((SphP[igas].Momentum[2] / P[igas].Mass) - SphP[igas].GroupVel[2]) - (GravAccel[2] -
                                                                                                                        SphP
                                                                                                                        [igas].GroupGravAcc
                                                                                                                        [2]) * ((SphP[igas].Momentum[1] / P[igas].Mass) - SphP[igas].GroupVel[1]);
  dir[1] =
    (GravAccel[2] - SphP[igas].GroupGravAcc[2]) * ((SphP[igas].Momentum[0] / P[igas].Mass) - SphP[igas].GroupVel[0]) - (GravAccel[0] -
                                                                                                                        SphP
                                                                                                                        [igas].GroupGravAcc
                                                                                                                        [0]) * ((SphP[igas].Momentum[2] / P[igas].Mass) - SphP[igas].GroupVel[2]);
  dir[2] =
    (GravAccel[0] - SphP[igas].GroupGravAcc[0]) * ((SphP[igas].Momentum[1] / P[igas].Mass) - SphP[igas].GroupVel[1]) - (GravAccel[1] -
                                                                                                                        SphP
                                                                                                                        [igas].GroupGravAcc
                                                                                                                        [1]) * ((SphP[igas].Momentum[0] / P[igas].Mass) - SphP[igas].GroupVel[0]);
#endif
#endif
#else
  theta = acos(2 * get_random_number() - 1);
  phi = 2 * M_PI * get_random_number();
  dir[0] = sin(theta) * cos(phi);
  dir[1] = sin(theta) * sin(phi);
  dir[2] = cos(theta);
#endif

  for(k = 0, norm = 0; k < 3; k++)
    norm += dir[k] * dir[k];

  norm = sqrt(norm);
  if(get_random_number() < 0.5)
    norm = -norm;

  /* create new wind particle */
  P[istar] = P[igas];
  P[istar].Type = 4;
  P[istar].SofteningType = All.SofteningTypeOfPartType[4];


  timebin_add_particle(&TimeBinsGravity, istar, igas, P[istar].TimeBinGrav, TimeBinSynchronized[P[istar].TimeBinGrav]);

  if(mass_of_wind >= P[igas].Mass)
    terminate("mass_of_wind > P[igas].Mass");

  P[istar].Mass = mass_of_wind;

#ifdef INDIVIDUAL_GRAVITY_SOFTENING
  if(((1 << P[istar].Type) & (INDIVIDUAL_GRAVITY_SOFTENING)))
    P[istar].SofteningType = get_softening_type_from_mass(P[istar].Mass);
#endif


#ifdef STELLARAGE
  P[istar].StellarAge = -All.WindFreeTravelMaxTimeFactor / All.cf_hubble_a;
#endif

#ifdef GFM
  gfm_add_star(istar, igas, P[istar].Mass, -All.WindFreeTravelMaxTimeFactor / All.cf_hubble_a, 1.5 * get_cell_radius(igas));
#endif

#ifdef MHD
  double Emag = 0.5 * (SphP[igas].B[0] * SphP[igas].B[0] + SphP[igas].B[1] * SphP[igas].B[1] + SphP[igas].B[2] * SphP[igas].B[2]) * SphP[igas].Volume * All.cf_atime;
  SphP[igas].Energy -= Emag;
#endif


  STP(istar).Hsml = 1.5 * get_cell_radius(igas);
  STP(istar).Utherm = SphP[igas].Utherm + utherm;

  /* now change the conserved quantities in the cell in proportion */
  double fac = (P[igas].Mass - P[istar].Mass) / P[igas].Mass;


#ifdef REFINEMENT_HIGH_RES_GAS
  STP(istar).HighResMass = SphP[igas].HighResMass * (1 - fac);
  /* SphP[igas].HighResMass *= fac; this will be done with the scalars later */
#endif


  P[igas].Mass *= fac;
  SphP[igas].Energy *= fac;
  SphP[igas].Momentum[0] *= fac;
  SphP[igas].Momentum[1] *= fac;
  SphP[igas].Momentum[2] *= fac;
#ifdef USE_ENTROPY_FOR_COLD_FLOWS
  SphP[igas].Entropy *= fac;
#endif

#ifdef MHD
  SphP[igas].Energy += Emag;
#endif

#ifdef MAXSCALARS
  for(int s = 0; s < N_Scalar; s++) /* Note, the changes in MATERIALS, HIGHRESGASMASS, etc., are treated as part of the Scalars */
    *(MyFloat *) (((char *) (&SphP[igas])) + scalar_elements[s].offset_mass) *= fac;
#endif
#ifdef GFM_CHEMTAGS
  for(int k = 0; k < GFM_N_CHEM_TAGS; k++)
    SphP[igas].MassMetalsChemTags[k] *= fac;
#endif

#if defined(TRACER_MC)
  P[istar].TracerHead = -1;     /* new wind particle starts with no tracers */
  P[istar].NumberOfTracers = 0; /* new wind particle starts with no tracers */
  consider_moving_tracers_local(igas, istar, 1 - fac);
#endif

  /* and finally, give the kick */
  if(norm != 0)
    {
      for(k = 0; k < 3; k++)
        dir[k] /= norm;

      for(k = 0; k < 3; k++)
        {
          P[istar].Vel[k] += (v * All.cf_atime * dir[k]);
        }
    }

  return;
}

#endif
