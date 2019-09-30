/* \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/blackhole/blackhole.c
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

/*! \file blackhole.c
*  \brief routines for gas accretion onto black holes, and black hole mergers
*/


#ifdef BLACK_HOLES

#if defined(BH_BUBBLES) && defined(BH_ADIOS_WIND)
#error "BH_ADIOS_WIND does not work with BH_BUBBLES"
#endif

 /* local data structure that holds results acquired on remote processors */
void blackhole_accretion(void)
{
  if(TimeBinsBHAccretion.GlobalNActiveParticles == 0)
    {
      mpi_printf("BLACK_HOLES: No active BHs\n");
      return;
    }

  TIMER_START(CPU_BH_ACCRETION);

  int idx, n, bin;
  double mass_real, total_mass_real, medd, total_mdoteddington;
  double mass_holes, total_mass_holes, mdot, total_mdot;
#if defined(BH_BUBBLES) || defined(BH_NF_RADIO)
  int num_activebh = 0;
#endif
#ifdef TRACER_MC_CHECKS
  long long check_total_tracers = get_total_number_of_tracers(-1);
  long long check_total_star_tracers = get_total_number_of_tracers(4);
  long long check_total_bh_tracers = get_total_number_of_tracers(5);
#endif

  /* for feedback statistics */
  AGNEnergyM_Should = AGNEnergyM_Is = AGNEnergyT_Should = AGNEnergyT_Is = 0.0;

  mpi_printf("BLACK_HOLES: Beginning black-hole accretion\n");

  /* quasar accretion rates */
  blackhole_calculate_mdot();

#if defined(BH_NF_RADIO)
  /* radio mode accretion rate */
  blackhole_calculate_mdot_radiomode();
#endif

#ifdef BH_DRAG
  blackhole_drag_force();
#endif

#ifdef BH_THERMALFEEDBACK_ACC
  blackhole_accumulate_energy();
#endif

#if defined(BH_THERMALFEEDBACK) || defined(BH_THERMALFEEDBACK_ACC)
  blackhole_assign_feedback();
#endif

#ifdef BH_ADIOS_WIND
  blackhole_blow_wind();
#endif

#ifdef BH_THERMALFEEDBACK_ACC
  blackhole_reset_energy();
#endif

#ifdef BH_SPIN_EVOLUTION
  blackhole_spin_evolution();
  myflush(FdBlackHolesSpin);
#endif

  blackhole_swallow_gas();

  /* recalculate timebin-based BH statistics for active bins */
  for(n = 0; n < TIMEBINS; n++)
    {
      if(TimeBinSynchronized[n])
        {
          TimeBin_BH_mass[n] = 0;
          TimeBin_BH_dynamicalmass[n] = 0;
          TimeBin_BH_Mdot[n] = 0;
          TimeBin_BH_Medd[n] = 0;
        }
    }

  /* find active radio mode BHs */
  for(idx = 0; idx < TimeBinsBHAccretion.NActiveParticles; idx++)
    {
      n = TimeBinsBHAccretion.ActiveParticleList[idx];
      if(n < 0)
        continue;

      bin = P[n].TimeBinHydro;
      TimeBin_BH_mass[bin] += BPP(n).BH_Mass;
      TimeBin_BH_dynamicalmass[bin] += P[n].Mass;
      TimeBin_BH_Mdot[bin] += BPP(n).BH_Mdot;
      if(BPP(n).BH_Mass > 0)
        TimeBin_BH_Medd[bin] += BPP(n).BH_Mdot / BPP(n).BH_Mass;
#ifdef BH_BUBBLES
      if(BPP(n).BH_Mass_bubbles > 0 && BPP(n).BH_Mass_bubbles > All.BlackHoleRadioTriggeringFactor * BPP(n).BH_Mass_ini)
        num_activebh++;
#endif
#ifdef BH_NF_RADIO
      if(BPP(n).BH_HaloVvir > 0)
        {
          double egy_thresh = blackhole_get_bubble_energy_thresh(n);

          if(BPP(n).BH_RadioEgyFeedback > egy_thresh)
            num_activebh++;
        }
#endif
    }

#ifdef BH_BUBBLES
  blackhole_do_bubbles(num_activebh);
#endif

#ifdef BH_NF_RADIO
  blackhole_do_bubbles_nf(num_activebh);
#endif

  /* assemble statistics for log-files */
  mdot = 0;
  mass_holes = 0;
  mass_real = 0;
  medd = 0;

  for(bin = 0; bin < TIMEBINS; bin++)
    if(TimeBinsBHAccretion.TimeBinCount[bin])
      {
        mass_holes += TimeBin_BH_mass[bin];
        mass_real += TimeBin_BH_dynamicalmass[bin];
        mdot += TimeBin_BH_Mdot[bin];
        medd += TimeBin_BH_Medd[bin];
      }

  double mass[4], mass_tot[4];
  mass[0] = mass_holes;
  mass[1] = mass_real;
  mass[2] = mdot;
  mass[3] = medd;
  MPI_Reduce(mass, mass_tot, 4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  total_mass_holes = mass_tot[0];
  total_mass_real = mass_tot[1];
  total_mdot = mass_tot[2];
  total_mdoteddington = mass_tot[3];

  if(ThisTask == 0)
    {
      /* convert to solar masses per yr */
      double mdot_in_msun_per_year = total_mdot * (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

      total_mdoteddington *= 1.0 / ((4 * M_PI * GRAVITY * CLIGHT * PROTONMASS / (All.BlackHoleRadiativeEfficiency * CLIGHT * CLIGHT * THOMPSON)) * All.UnitTime_in_s / All.HubbleParam);

      fprintf(FdBlackHoles, "%14e %6d  %14e %14e %14e %14e %14e\n", All.Time, All.TotNumBHs, total_mass_holes, total_mdot, mdot_in_msun_per_year, total_mass_real, total_mdoteddington);
      myflush(FdBlackHoles);
    }

  myflush(FdBlackHolesDetails);

#ifdef TRACER_MC_CHECKS
  if(check_total_tracers != get_total_number_of_tracers(-1))
    terminate("TRACER_MC: unconserved tracers during BLACK_HOLES (accretion): total number of tracers BEFORE = %lld, AFTER = %lld\n", check_total_tracers, get_total_number_of_tracers(-1));
  if(check_total_bh_tracers > get_total_number_of_tracers(5))
    terminate("TRACER_MC: strange, global number of BH tracers decreased in BHACC, BEFORE = %lld, AFTER = %lld\n", check_total_bh_tracers, get_total_number_of_tracers(5));
  if(check_total_star_tracers != get_total_number_of_tracers(4))
    terminate("TRACER_MC: strange, global number of STAR tracers changed in BHACC, BEFORE = %lld, AFTER = %lld\n", check_total_star_tracers, get_total_number_of_tracers(4));
#endif

  TIMER_STOP(CPU_BH_ACCRETION);
}


void blackhole_energy_log_info(void)
{
  return;                       /* the info below can be enabled for debugging purposes */


  TIMER_START(CPU_LOGS);

  double AGNEnergyEMobsInStep_Is, AGNEnergyEMInStep_Is;
  double totAGNEnergyT_Should, totAGNEnergyT_Is, totAGNEnergyM_Should, totAGNEnergyM_Is, totAGNEnergyEMobs_Is, totAGNEnergyEM_Is;

  AGNEnergyEMobsInStep_Is = AGNEnergyEMobs_Is * get_time_difference_in_Gyr(All.Time - All.TimeStep, All.Time) * SEC_PER_GIGAYEAR;
  AGNEnergyEMInStep_Is = AGNEnergyEM_Is * get_time_difference_in_Gyr(All.Time - All.TimeStep, All.Time) * SEC_PER_GIGAYEAR;

  /* per step */
  /*
     MPI_Reduce(&AGNEnergyT_Should, &totAGNEnergyT_Should, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
     MPI_Reduce(&AGNEnergyT_Is, &totAGNEnergyT_Is, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
     MPI_Reduce(&AGNEnergyM_Should, &totAGNEnergyM_Should, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
     MPI_Reduce(&AGNEnergyM_Is, &totAGNEnergyM_Is, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
     MPI_Reduce(&AGNEnergyEMobsInStep_Is, &totAGNEnergyEMobs_Is, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
     MPI_Reduce(&AGNEnergyEMInStep_Is, &totAGNEnergyEM_Is, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   */

  double in[6], out[6];
  in[0] = AGNEnergyT_Should;
  in[1] = AGNEnergyT_Is;
  in[2] = AGNEnergyM_Should;
  in[3] = AGNEnergyM_Is;
  in[4] = AGNEnergyEMobsInStep_Is;
  in[5] = AGNEnergyEMInStep_Is;

  MPI_Reduce(in, out, 6, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  totAGNEnergyT_Should = out[0];
  totAGNEnergyT_Is = out[1];
  totAGNEnergyM_Should = out[2];
  totAGNEnergyM_Is = out[3];
  totAGNEnergyEMobs_Is = out[3];
  totAGNEnergyEM_Is = out[5];


  /* cumulative energies */
  All.CumAGNEnergyT_Should += totAGNEnergyT_Should;
  All.CumAGNEnergyT_Is += totAGNEnergyT_Is;
  All.CumAGNEnergyM_Should += totAGNEnergyM_Should;
  All.CumAGNEnergyM_Is += totAGNEnergyM_Is;
  All.CumAGNEnergyEMobs_Is += totAGNEnergyEMobs_Is;
  All.CumAGNEnergyEM_Is += totAGNEnergyEM_Is;

  /* NOTE: the thermal energy per step is not expected to be balanced, because only active cells actually receive it in the same time step. The cumulative should be exact. */
  mpi_printf("BLACK_HOLES: AGN feedback channels (step)        Is      --> thermal = %g erg    mechanical = %g erg    electro-magnetic = (%g/%g) erg\n",
             totAGNEnergyT_Is * All.UnitEnergy_in_cgs, totAGNEnergyM_Is * All.UnitEnergy_in_cgs, totAGNEnergyEMobs_Is, totAGNEnergyEM_Is);
  mpi_printf("BLACK_HOLES: AGN feedback channels (step)        Should  --> thermal = %g erg    mechanical = %g erg\n",
             totAGNEnergyT_Should * All.UnitEnergy_in_cgs, totAGNEnergyM_Should * All.UnitEnergy_in_cgs);
  mpi_printf("BLACK_HOLES: AGN feedback channels (step)        Ratio   --> thermal = %g        mechanical = %g \n",
             (totAGNEnergyT_Should > 0) ? (totAGNEnergyT_Is / totAGNEnergyT_Should) : (0), (totAGNEnergyM_Should > 0) ? (totAGNEnergyM_Is / totAGNEnergyM_Should) : (0));

  mpi_printf("BLACK_HOLES: AGN feedback channels (cumulative)  Is      --> thermal = %g erg    mechanical = %g erg    electro-magnetic = (%g/%g) erg\n",
             All.CumAGNEnergyT_Is * All.UnitEnergy_in_cgs, All.CumAGNEnergyM_Is * All.UnitEnergy_in_cgs, All.CumAGNEnergyEMobs_Is, All.CumAGNEnergyEM_Is);
  mpi_printf("BLACK_HOLES: AGN feedback channels (cumulative)  Should  --> thermal = %g erg    mechanical = %g erg\n",
             All.CumAGNEnergyT_Should * All.UnitEnergy_in_cgs, All.CumAGNEnergyM_Should * All.UnitEnergy_in_cgs);
  mpi_printf("BLACK_HOLES: AGN feedback channels (cumulative)  Ratio   --> thermal = %g       mechanical = %g \n",
             (All.CumAGNEnergyT_Should > 0) ? (All.CumAGNEnergyT_Is / All.CumAGNEnergyT_Should) : (0), (All.CumAGNEnergyM_Should > 0) ? (All.CumAGNEnergyM_Is / All.CumAGNEnergyM_Should) : (0));

  TIMER_STOP(CPU_LOGS);
}

#endif
