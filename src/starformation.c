/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/starformation.c
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
#include "forcetree.h"

/** \file starformation.c
 *  \brief Generic creation routines for star and wind particles.
 *
 *  Star formation rates are calculated in sfr_eEOS for the multiphase model 
 *  and in fm_star_formation/sfr.c for stellar feedback ISM.
 */

#ifdef USE_SFR

static int stars_spawned;                   /**< local number of star particles spawned in the time step */
static int tot_stars_spawned;               /**< global number of star paricles spawned in the time step */
static int stars_converted;                 /**< local number of gas cells converted into stars in the time step */
static int tot_stars_converted;             /**< global number of gas cells converted into stars in the time step */
static int altogether_spawned;              /**< local number of star+wind particles spawned in the time step */
static int tot_altogether_spawned;          /**< global number of star+wind particles spawned in the time step */
static double cum_mass_stars = 0.0;         /**< cumulative mass of stars created in the time step (global value) */

#ifdef GFM_WINDS
static int wind_spawned;                    /**< local number of wind particles spawned in the time step */
static int wind_converted;                  /**< local number of gas cells convered into wind particles in the time step */
static int tot_wind_converted;              /**< global number of gas cells convered into wind particles in the time step */
static int tot_wind_spawned;                /**< global number of wind particles spawned in the time step */
static double windenergy_local_sum_should;  /**< local expected value for the wind energy */
static double windenergy_local_sum_is;      /**< local actual value of the wind energy for the created wind particles */
#ifdef VERBOSE
static double windenergy_global_sum_should; /**< global expected value for the wind energy */
static double windenergy_global_sum_is;     /**< global actual value of the wind energy for the created wind particles */
#endif
#endif

/** \brief Initialization routine
 */
#if !defined(LOCAL_FEEDBACK)
static int sfr_init_called = 0;
void sfr_init()
{
  if(sfr_init_called)
    return;

  sfr_init_called = 1;

#if !defined(FM_SFR) && !defined(ISM)
  init_clouds();
#else
  init_star_formation();
#if defined(FM_STAR_FEEDBACK) && defined(DELAYED_COOLING_TURB)
  init_turbulent_energy();
#endif
#if defined(FM_STAR_FEEDBACK) && defined(INSTANTANEOUS_DEPOSITION)
  init_instantaneous_deposition();
#endif
#endif
}
#endif

/** \brief This routine creates star/wind particles according to their respective rates.
 *
 *  This function loops over all the active gas cells. If in a given cell the SFR is 
 *  greater than zero, the probability of forming a star or a wind particle is computed
 *  and the corresponding particle is created stichastically according to the model 
 *  in Springel & Hernquist (2003, MNRAS). It also saves information about the formed stellar
 *  mass and the star formation rate in the file FdSfr.
 */
void sfr_create_star_particles(void)
{
  TIMER_START(CPU_COOLINGSFR);

  int idx, i, bin;
  double dt, dtime;
  MyDouble mass_of_star;
  double sum_sm, total_sm, rate, sum_mass_stars, total_sum_mass_stars;
  double p = 0, pall = 0, prob, p_decide;
  double rate_in_msunperyear;
  double sfrrate, totsfrrate;

#ifdef DO_NOT_CREATE_STAR_PARTICLES
  return;
#endif

#ifdef COSMIC_RAYS_EXTRA_DIAGNOSTICS
  double InitialCREnergy = 0;
  for(int i=0; i<NumGas; i++)
    if(P[i].Mass != 0 && P[i].ID != 0 && P[i].Type == 0)
      InitialCREnergy += SphP[i].CR_Energy;
#endif

#ifdef METALS
  double w = 0;
#endif

#ifdef GFM_WINDS
  double v_wind = 0, u_wind = 0, p_wind;
#endif

#ifdef GFM_WINDS
  WindEnergy_Should = 0.0;
  WindEnergy_Is = 0.0;
#endif

  stars_spawned = stars_converted = 0;
#ifdef GFM_WINDS
  wind_spawned = wind_converted = 0;
  windenergy_local_sum_should = 0.0, windenergy_local_sum_is = 0.0;
#endif
  sum_sm = sum_mass_stars = 0;

#ifdef REFINEMENT_AROUND_BH
#ifdef SUPPRESS_SF_IN_REFINEMENT_REGION
  mpi_printf("BLACK_HOLES: Updating refinement flag for cells around BHs prior to star formation\n");
  blackhole_mark_cells_for_refinement();
#endif
#endif

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i >= 0)
        {
          if(P[i].Mass == 0 && P[i].ID == 0)
            continue;           /* skip cells that have been swallowed or eliminated */

#ifdef REFINEMENT_AROUND_BH
#ifdef SUPPRESS_SF_IN_REFINEMENT_REGION
        if(SphP[i].RefBHFlag) {
            continue; /* skip cells that are being refined.  */
        }
#endif
#endif

#ifdef SFR_KEEP_CELLS
          if(P[i].Mass < 0.3 * All.TargetGasMass)
            continue;
#endif

#if defined(ISM_HII_HEATING) || defined(FM_RADIATION_FEEDBACK)
          if(SphP[i].GasRadCoolShutoffTime > 0)
            continue;
#endif

          dt = (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;

          /*  the actual time-step */

          dtime = All.cf_atime * dt / All.cf_time_hubble_a;

          mass_of_star = 0;
          prob = 0;
          p = 0;
          pall = 0;

          if(SphP[i].Sfr > 0)
            {
              p = SphP[i].Sfr / ((All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR)) * dtime / P[i].Mass;
              pall = p;
              sum_sm += P[i].Mass * (1 - exp(-p));

#ifdef METALS
              SphP[i].Metallicity = SphP[i].MassMetallicity / P[i].Mass;        /* update this here because primitive variables have not been updated */
              w = get_random_number();
              assert(SphP[i].Metallicity >= 0);
              SphP[i].Metallicity += w * METAL_YIELD * (1 - exp(-p));
              SphP[i].MassMetallicity = SphP[i].Metallicity * P[i].Mass;
              P[i].Metallicity = SphP[i].Metallicity;
              assert(P[i].Metallicity >= 0);
#endif

#if defined(REFINEMENT_SPLIT_CELLS) && defined(REFINEMENT_MERGE_CELLS)

              if(P[i].Mass < 2.0 * All.TargetGasMass)
#ifdef SFR_KEEP_CELLS
                mass_of_star = 0.9 * P[i].Mass;
#else
                mass_of_star = P[i].Mass;
#endif
              else
                mass_of_star = All.TargetGasMass;

#ifdef REFINEMENT_HIGH_RES_GAS
              if(SphP[i].HighResMass < HIGHRESMASSFAC * P[i].Mass)
                {
                  /* this cell does not appear to be in the high-res region.
                     If we form a star, then it is given the mass of the cell,
                     and later we give the star the SofteningType=3 particle to give it large softening */
#ifdef SFR_KEEP_CELLS
                  mass_of_star = 0.9 * P[i].Mass;
#else
                  mass_of_star = P[i].Mass;
#endif
                }

#endif /* REFINEMENT_SPLIT_CELLS && REFINEMENT_MERGE_CELLS */

#else
              mass_of_star = P[i].Mass;
#endif

#ifdef SFR_KEEP_CELLS
              if(P[i].Mass < 0.5 * All.TargetGasMass)
                continue;       /* do not make stars from cells that should be derefined */
#endif

#ifdef GFM_WINDS
              gfm_calc_wind_parameters(i, p, &p_wind, &v_wind, &u_wind);
              pall += p_wind;

              windenergy_local_sum_should += 0.5 * (p_wind * P[i].Mass) * v_wind * v_wind + (p_wind * P[i].Mass) * u_wind;
#endif
              prob = P[i].Mass / mass_of_star * (1 - exp(-pall));
            }

          if(prob == 0)
            continue;

          if(prob < 0)
            terminate("prob < 0");

          if(prob > 1)
            {
              printf
                ("SFR: Warning, need to make a heavier star than desired. Task=%d prob=%g P[i].Mass=%g mass_of_star=%g mass_of_star_new=%g p=%g pall=%g\n",
                 ThisTask, prob, P[i].Mass, mass_of_star, P[i].Mass * (1 - exp(-pall)), p, pall);
              mass_of_star = P[i].Mass * (1 - exp(-pall));
              prob = 1.0;
            }

          /* decide what process to consider (currently available: make a star or kick to wind) */
          p_decide = get_random_number();

          if(p_decide < p / pall)       /* ok, it is decided to consider star formation */
            make_star(idx, i, prob, mass_of_star, &sum_mass_stars);

#ifdef GFM_WINDS
          if(p_decide >= p / pall)      /* ok, it is decided to consider winds */
            make_wind(idx, i, prob, mass_of_star, v_wind, u_wind);
#endif

#ifdef METALS
          if(SphP[i].Sfr > 0)
            {
              if(P[i].Type == 0)        /* to protect using a particle that has been turned into a star */
                {
                  SphP[i].Metallicity += (1 - w) * METAL_YIELD * (1 - exp(-p));
                  SphP[i].MassMetallicity = SphP[i].Metallicity * P[i].Mass;
                  assert(SphP[i].Metallicity >= 0);

                }
            }
          P[i].Metallicity = SphP[i].Metallicity;
          assert(P[i].Metallicity >= 0);
#endif
        }
    }                           /* end of main loop over active gas particles */


  int in[4], out[4], cnt = 2;
  in[0] = stars_spawned;
  in[1] = stars_converted;
#ifdef GFM_WINDS
  in[2] = wind_spawned;
  in[3] = wind_converted;
  cnt = 4;
#endif

  MPI_Allreduce(in, out, cnt, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  tot_stars_spawned = out[0];
  tot_stars_converted = out[1];

#ifdef GFM
  All.TotN_star += tot_stars_spawned + tot_stars_converted;
#ifdef GFM_WINDS
  tot_wind_spawned = out[2];
  tot_wind_converted = out[3];
  All.TotN_star += tot_wind_spawned + tot_wind_converted;
#endif
#endif

#ifdef GFM_WINDS
  if(tot_wind_spawned > 0 || tot_wind_converted > 0)
    mpi_printf("GFM_WINDS: spawned %d wind particles, converted %d gas particles into wind particles\n", tot_wind_spawned, tot_wind_converted);
#endif

  if(tot_stars_spawned > 0 || tot_stars_converted > 0)
    mpi_printf("SFR: spawned %d stars, converted %d gas particles into stars\n", tot_stars_spawned, tot_stars_converted);

  tot_altogether_spawned = tot_stars_spawned;
  altogether_spawned = stars_spawned;
#ifdef GFM_WINDS
  tot_altogether_spawned += tot_wind_spawned;
  altogether_spawned += wind_spawned;
#endif
  if(tot_altogether_spawned)
    {
      /* need to assign new unique IDs to the spawned stars */

      int *list;

      if(All.MaxID == 0)        /* MaxID not calculated yet */
        calculate_maxid();

      list = mymalloc("list", NTask * sizeof(int));

      MPI_Allgather(&altogether_spawned, 1, MPI_INT, list, 1, MPI_INT, MPI_COMM_WORLD);

      MyIDType newid = All.MaxID + 1;

      for(i = 0; i < ThisTask; i++)
        newid += list[i];

      myfree(list);

      for(i = 0; i < altogether_spawned; i++)
        {
          P[NumPart + i].ID = newid;

#ifdef TRACER_MC_CHECKS
          int next = P[NumPart + i].TracerHead;

          while(next >= 0)      /* traverse forwards to end */
            {
              TracerLinkedList[next].ParentID = newid;

              next = TracerLinkedList[next].Next;
            }
#endif
          newid++;
        }

      All.MaxID += tot_altogether_spawned;
    }

  /* Note: New tree construction can be avoided because of  `force_add_star_to_tree()' */
  if(tot_stars_spawned > 0 || tot_stars_converted > 0)
    {
      All.TotNumPart += tot_stars_spawned;
      All.TotNumGas -= tot_stars_converted;
      NumPart += stars_spawned;
    }
#ifdef GFM_WINDS
  if(tot_wind_spawned > 0 || tot_wind_converted > 0)
    {
      All.TotNumPart += tot_wind_spawned;
      All.TotNumGas -= tot_wind_converted;
      NumPart += wind_spawned;
    }
#ifdef VERBOSE
  MPI_Reduce(&windenergy_local_sum_should, &windenergy_global_sum_should, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&windenergy_local_sum_is, &windenergy_global_sum_is, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  WindEnergy_Should += windenergy_global_sum_should;
  WindEnergy_Is += windenergy_global_sum_is;
#endif
#endif

  for(bin = 0, sfrrate = 0; bin < TIMEBINS; bin++)
    if(TimeBinsHydro.TimeBinCount[bin])
      sfrrate += TimeBinSfr[bin];

  double din[3] = {sfrrate, sum_sm, sum_mass_stars}, dout[3];

  MPI_Reduce(din, dout, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      totsfrrate = dout[0];
      total_sm = dout[1];
      total_sum_mass_stars = dout[2];

      if(All.TimeStep > 0)
        rate = total_sm / (All.TimeStep / All.cf_time_hubble_a);
      else
        rate = 0;

      /* compute the cumulative mass of stars (->>> CHECK ME!!!) */
      cum_mass_stars += total_sum_mass_stars;

      /* convert to solar masses per yr */
      rate_in_msunperyear = rate * (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

      fprintf(FdSfr, "%14e %14e %14e %14e %14e %14e\n", All.Time, total_sm, totsfrrate, rate_in_msunperyear, total_sum_mass_stars, cum_mass_stars);
      myflush(FdSfr);
    }

#ifdef COSMIC_RAYS_EXTRA_DIAGNOSTICS
  double FinalCREnergy = 0;
  for(int i=0; i<NumGas; i++)
    if(P[i].Mass != 0 && P[i].ID != 0 && P[i].Type == 0)
      FinalCREnergy += SphP[i].CR_Energy;
  
  double dCREnergy = FinalCREnergy - InitialCREnergy;
  All.TotalCREnergyLossSfr += dCREnergy;
#endif

  TIMER_STOP(CPU_COOLINGSFR);
}



/** \brief Convert a cell into a star. 
 *
 *  This function convertss an active star-forming gas cell into a star. 
 *  The particle information of the gas cell is copied to the
 *  location istar and the fields necessary for the creation of the star
 *  particle are initialized. 
 *

 *  \param i index of the gas cell to be converted
 *  \param birthtime time of birth (in code units) of the stellar particle
 */
void convert_cell_into_star(int i, double birthtime)
{
  P[i].Type = 4;
  P[i].SofteningType = All.SofteningTypeOfPartType[P[i].Type];

#if defined(REFINEMENT_HIGH_RES_GAS) && !defined(TGSET)
  if(SphP[i].HighResMass < HIGHRESMASSFAC * P[i].Mass)
    {
      /* this cell does not appear to be in the high-res region.
         We give the star the SofteningType=3 particle to give it large softening */
      P[i].SofteningType = All.SofteningTypeOfPartType[3];
    }
#endif

#ifdef INDIVIDUAL_GRAVITY_SOFTENING
  if(((1 << P[i].Type) & (INDIVIDUAL_GRAVITY_SOFTENING)))
    P[i].SofteningType = get_softening_type_from_mass(P[i].Mass);
#endif

#ifdef GFM
  gfm_add_star(i, i, P[i].Mass, birthtime, 1.5 * get_cell_radius(i));

#ifdef TRACER_MC
  int i_tracer = P[i].TracerHead;
  while(i_tracer >= 0)
    {
#if (TRACER_MC_LAST_STAR_TIME)
      if(birthtime < 0)         /* changed to a wind particle */
        TracerLinkedList[i_tracer].fluid_quantities[TracerMCLastStarTimeIndex] += 3 * All.TimeMax;
      else                      /* changed to a normal star particle */
        TracerLinkedList[i_tracer].fluid_quantities[TracerMCLastStarTimeIndex] = 2 * All.TimeMax;
#endif
#if (TRACER_MC_WIND_COUNTER)
      if(birthtime < 0)         /* changed to a wind particle */
        TracerLinkedList[i_tracer].fluid_quantities[TracerMCWindCounterIndex] += 1;
#endif
      i_tracer = TracerLinkedList[i_tracer].Next;
    }
#endif /* TRACER_MC */
#endif /* GFM */


  TimeBinSfr[P[i].TimeBinHydro] -= SphP[i].Sfr;

#ifdef STELLARAGE
  P[i].StellarAge = birthtime;
#endif

#ifdef VORONOI_DYNAMIC_UPDATE
  voronoi_remove_connection(i);
#endif

  return;
}


/** \brief Spawn a star particle from a gas cell.
 *
 *  This function spawns a star particle from an active star-forming
 *  cell. The particle information of the gas cell is copied to the
 *  location istar and the fields necessary for the creation of the star
 *  particle are initialized. The conserved variables of the gas cell
 *  are then updated according to the mass ratio between the two components
 *  to ensure conservation.
 *
 *  \param igas index of the gas cell from which the star is spawned
 *  \param birthtime time of birth (in code units) of the stellar particle
 *  \param istar index of the spawned stellar particle
 *  \param mass_of_star the mass of the spawned stellar particle
 */
void spawn_star_from_cell(int igas, double birthtime, int istar, MyDouble mass_of_star)
{
  P[istar] = P[igas];
  P[istar].Type = 4;
  P[istar].SofteningType = All.SofteningTypeOfPartType[P[istar].Type];
  P[istar].Mass = mass_of_star;


#if defined(REFINEMENT_HIGH_RES_GAS) && !defined(TGSET)
  if(SphP[igas].HighResMass < HIGHRESMASSFAC * P[igas].Mass)
    {
      /* this cell does not appear to be in the high-res region.
         We give the star the SofteningType=3 particle to give it large softening */
      P[istar].SofteningType = All.SofteningTypeOfPartType[3];
    }
#endif

#ifdef INDIVIDUAL_GRAVITY_SOFTENING
  if(((1 << P[istar].Type) & (INDIVIDUAL_GRAVITY_SOFTENING)))
    P[istar].SofteningType = get_softening_type_from_mass(P[istar].Mass);
#endif



  timebin_add_particle(&TimeBinsGravity, istar, igas, P[istar].TimeBinGrav, TimeBinSynchronized[P[istar].TimeBinGrav]);

#ifdef STELLARAGE
  P[istar].StellarAge = birthtime;
#endif

  /* now change the conserved quantities in the cell in proportion */
  double fac = (P[igas].Mass - P[istar].Mass) / P[igas].Mass;

#ifdef GFM
  gfm_add_star(istar, igas, mass_of_star, birthtime, 1.5 * get_cell_radius(igas));
#endif

#ifdef MHD
  double Emag = 0.5 * (SphP[igas].B[0] * SphP[igas].B[0] + SphP[igas].B[1] * SphP[igas].B[1] + SphP[igas].B[2] * SphP[igas].B[2]) * SphP[igas].Volume * All.cf_atime;
  SphP[igas].Energy -= Emag;
#endif

  P[igas].Mass *= fac;
  SphP[igas].Energy *= fac;
  SphP[igas].Momentum[0] *= fac;
  SphP[igas].Momentum[1] *= fac;
  SphP[igas].Momentum[2] *= fac;

#ifdef MHD
  SphP[igas].Energy += Emag;
#endif

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
  SphP[igas].Entropy *= fac;
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
  P[istar].TracerHead = -1;     /* new star particle starts with no tracer children */
  P[istar].NumberOfTracers = 0;
  consider_moving_tracers_local(igas, istar, 1 - fac);
#endif

  return;
}

/** \brief Make a star particle from a gas cell.
 *
 *  Given a gas cell where star formation is active and the probability
 *  of forming a star, this function selectes either to convert the gas 
 *  cell into a star particle or to spawn a star depending on the 
 *  target mass for the star.
 *  
 *  \param idx index of the gas cell in the hydro list of active cells  
 *  \param i index of the gas cell
 *  \param prob probability of making a star
 *  \param mass_of_star desired mass of the star particle
 *  \param sum_mass_stars holds the mass of all the stars created at the current time-step (for the local task)
 */
void make_star(int idx, int i, double prob, MyDouble mass_of_star, double *sum_mass_stars)
{
  if(mass_of_star > P[i].Mass)
    terminate("mass_of_star > P[i].Mass");

  if(get_random_number() < prob)
    {
      if(mass_of_star == P[i].Mass)
        {
          /* here we turn the gas particle itself into a star particle */
          Stars_converted++;
          stars_converted++;

          *sum_mass_stars += P[i].Mass;

          convert_cell_into_star(i, All.Time);
          timebin_remove_particle(&TimeBinsHydro, idx, P[i].TimeBinHydro);
#if defined(FM_STAR_FEEDBACK) && defined(OUTPUT_STELLAR_FEEDBACK)
          update_cumulative_feedback_energy_of_converted_cells(i);
#endif
#ifdef COSMIC_RAYS
          P[i].CRInjection = 1;
          P[i].Hsml = 4. * get_cell_radius(i);
#endif
        }
      else
        {
          /* in this case we spawn a new star particle, only reducing the mass in the cell by mass_of_star */
          altogether_spawned = stars_spawned;
#ifdef GFM_WINDS
          altogether_spawned += wind_spawned;
#endif
          if(NumPart + altogether_spawned >= All.MaxPart)
            terminate("NumPart=%d spwawn %d particles no space left (All.MaxPart=%d)\n", NumPart, altogether_spawned, All.MaxPart);

          int j = NumPart + altogether_spawned; /* index of new star */

          spawn_star_from_cell(i, All.Time, j, mass_of_star);

          *sum_mass_stars += mass_of_star;
          stars_spawned++;
#ifdef COSMIC_RAYS
          P[j].CRInjection = 1;
          P[j].Hsml = 4. * get_cell_radius(i);
#endif
        }
    }
}


#ifdef GFM_WINDS
/** \brief Make a wind particle from a gas cell.
 *
 *  Given a gas cell and the probability of generating a wind particle, 
 *  this function selectes either to convert the gas cell into a wind particle 
 *  or to spawn a wind particle depending on the target mass.
 *
 *  \param idx index of the gas cell in the hydro list of active cells
 *  \param i index of the gas cell
 *  \param prob probability of making a wind particle
 *  \param mass_of_wind desired mass of the wind particle
 *  \param v_wind wind speed
 */
void make_wind(int idx, int i, double prob, MyDouble mass_of_wind, double v_wind, double u_wind)
{
  if(mass_of_wind > P[i].Mass)
    terminate("mass_of_wind > P[i].Mass");

  if(get_random_number() < prob)
    {
      if(mass_of_wind == P[i].Mass)
        {
          /* here we turn the gas particle itself into a wind particle */
          gfm_add_wind(i, v_wind, u_wind);
          timebin_remove_particle(&TimeBinsHydro, idx, P[i].TimeBinHydro);
          windenergy_local_sum_is += 0.5 * P[i].Mass * v_wind * v_wind + P[i].Mass * u_wind;
          Stars_converted++;
          wind_converted++;
        }
      else
        {
          /* in this case we spawn a new wind particle, only reducing the mass in the cell by mass_of_star */
          altogether_spawned = stars_spawned + wind_spawned;
          if(NumPart + altogether_spawned >= All.MaxPart)
            terminate("NumPart=%d spwawn %d particles no space left (All.MaxPart=%d)\n", NumPart, altogether_spawned, All.MaxPart);

          int j = NumPart + altogether_spawned; /* index of new wind particle */

          gfm_spawn_wind_from_cell(i, v_wind, u_wind, j, mass_of_wind);
          windenergy_local_sum_is += 0.5 * mass_of_wind * v_wind * v_wind + mass_of_wind * u_wind;
          wind_spawned++;
        }
    }
}

#endif
#endif /* closes SFR */
