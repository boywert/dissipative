/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/logs.c
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
#include <unistd.h>
#include <ctype.h>
#include <sys/stat.h>

#ifdef VTUNE_INSTRUMENT
#include <ittnotify.h>
#endif


#include "allvars.h"
#include "proto.h"
#include "voronoi.h"

/*! \file logs.c
 *  \brief Log-files handling
 */

#define CPU_STRING_LEN     120

/*! \brief Contains informations about the used CPU timers like it's name, symbols,...
 *
 */
struct timer_d Timer_data[CPU_LAST + 1];

#ifdef VTUNE_INSTRUMENT
__itt_domain* vtune_domain;
#endif



enum timers TimerStack[TIMER_STACK_DEPTH];
int TimerStackPos = 0;


/*! \brief Open files for logging.
 *
 *   This function opens various log-files that report on the status and
 *   performance of the simulation. Upon restart, the code will append to 
 *   these files.
 */
void open_logfiles(void)
{
  char mode[2], buf[1000], msg[1000];

  if(RestartFlag == 0)
    strcpy(mode, "w");
  else
    strcpy(mode, "a");

  if(ThisTask == 0)
    mkdir(All.OutputDir, 02755);
  MPI_Barrier(MPI_COMM_WORLD);

#ifdef DETAILEDTIMINGS
  sprintf(buf, "%stimings_detailed_%d.txt", All.OutputDir, ThisTask);
  if(!(FdDetailed = fopen(buf, mode)))
    terminate("error in opening file '%s'\n", buf);
#endif

#ifdef BLACK_HOLES
  if(ThisTask == 0)
    {
      sprintf(buf, "%sblackhole_details", All.OutputDir);
      mkdir(buf, 02755);
    }
  MPI_Barrier(MPI_COMM_WORLD);
  /* Note: This is done by everyone */
  sprintf(buf, "%sblackhole_details/blackhole_details_%d.txt", All.OutputDir, ThisTask);
  if(!(FdBlackHolesDetails = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }
#ifdef BH_SPIN_EVOLUTION
  sprintf(buf, "%sblackhole_details/blackhole_spin_%d.txt", All.OutputDir, ThisTask);
  if(!(FdBlackHolesSpin = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }
#endif
#ifdef BH_NEW_CENTERING
  sprintf(buf, "%sblackhole_details/blackhole_repositioning_%d.txt", All.OutputDir, ThisTask);
  if(!(FdBlackHolesRepos = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }
#endif
  if(ThisTask == 0)
    {
      sprintf(buf, "%sblackhole_mergers", All.OutputDir);
      mkdir(buf, 02755);
    }
  MPI_Barrier(MPI_COMM_WORLD);
  /* Note: This is done by everyone */
  sprintf(buf, "%sblackhole_mergers/blackhole_mergers_%d.txt", All.OutputDir, ThisTask);
  if(!(FdBlackHolesMergers = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }

#ifdef BH_BIPOLAR_FEEDBACK
  if(ThisTask == 0)
    {
      sprintf(buf, "%sblackhole_bipolar", All.OutputDir);
      mkdir(buf, 02755);
    }
  MPI_Barrier(MPI_COMM_WORLD);
  sprintf(buf, "%sblackhole_bipolar/blackhole_bipolar_%d.txt", All.OutputDir, ThisTask);
  if(!(FdBlackHolesBipolar = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }
#endif

#endif

  if(ThisTask != 0)             /* only the root processors writes to the log files */
    return;

  sprintf(buf, "%s%s", All.OutputDir, "cpu.txt");
  if(!(FdCPU = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }

  sprintf(buf, "%s%s", All.OutputDir, "info.txt");
  if(!(FdInfo = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }

  sprintf(buf, "%s%s", All.OutputDir, "energy.txt");
  if(!(FdEnergy = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }

  sprintf(buf, "%s%s", All.OutputDir, "timings.txt");
  if(!(FdTimings = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }

  sprintf(buf, "%s%s", All.OutputDir, "balance.txt");
  if(!(FdBalance = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }

  sprintf(buf, "%s%s", All.OutputDir, "timebins.txt");
  if(!(FdTimebin = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }

  sprintf(buf, "%s%s", All.OutputDir, "domain.txt");
  if(!(FdDomain = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }

  sprintf(buf, "%s%s", All.OutputDir, "memory.txt");
  if(!(FdMemory = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }

#ifdef LOCAL_FEEDBACK
  sprintf(buf, "%slocal_feedback.txt", All.OutputDir);
  if(!(FdLocalFeedback = fopen(buf, mode)))
    terminate("error in opening file '%s'\n", buf);
#endif

#ifdef FORCETEST
  sprintf(buf, "%s%s", All.OutputDir, "forcetest.txt");
  if(!(FdForceTest = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }
  fclose(FdForceTest);
#endif

#ifdef RESTART_DEBUG
  sprintf(buf, "%s%s", All.OutputDir, "restartdebug.txt");
  if(!(FdRestartTest = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }
#endif

#ifdef OUTPUT_CPU_CSV
  sprintf(buf, "%s%s", All.OutputDir, "cpu.csv");
  if(!(FdCPUCSV = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }
#endif

#if defined (VS_TURB) || defined (AB_TURB)
  sprintf(buf, "%s%s", All.OutputDir, "turb.txt");
  if(!(FdTurb = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }
#endif

#ifdef RT_ADVECT
  sprintf(buf, "%s%s", All.OutputDir, "radtransfer.txt");
  if(!(FdRad = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }
#endif

#ifdef GFM_STELLAR_EVOLUTION
  sprintf(buf, "%s%s", All.OutputDir, "metals_tot.txt");
  if(!(FdMetalsTot = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      terminate(msg);
    }
  sprintf(buf, "%s%s", All.OutputDir, "metals_gas.txt");
  if(!(FdMetalsGas = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      terminate(msg);
    }
  sprintf(buf, "%s%s", All.OutputDir, "metals_stars.txt");
  if(!(FdMetalsStars = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      terminate(msg);
    }
  sprintf(buf, "%s%s", All.OutputDir, "SN.txt");
  if(!(FdSN = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      terminate(msg);
    }
#endif

#if defined(FM_STAR_FEEDBACK) && defined(OUTPUT_STELLAR_FEEDBACK)
  sprintf(buf, "%s%s", All.OutputDir, "stellar_feedback.txt");
  if(!(FdFeedback = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      terminate(msg);
    }

  sprintf(buf, "%s%s", All.OutputDir, "gas_feedback.txt");
  if(!(FdGasFeed = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      terminate(msg);
    }
#endif

  int i = 0;
  fprintf(FdBalance, "\n");

#ifdef OUTPUT_CPU_CSV
  fprintf(FdCPUCSV, "STEP, TIME, CPUS, MULTIPLEDOMAIN, HIGHESTTIMEBIN, ");
#endif
  for(; i < CPU_LAST; i++)
    {
      if(Timer_data[i].symb != 0 && Timer_data[i].symbImbal != 0)
        {
          fprintf(FdBalance, "%-20s = '%c' / '%c'\n", Timer_data[i].longname, Timer_data[i].symb, Timer_data[i].symbImbal);
        }
#ifdef OUTPUT_CPU_CSV
      fprintf(FdCPUCSV, "%s1, %s2, %s3, ", Timer_data[i].shortname, Timer_data[i].shortname, Timer_data[i].shortname);
#endif
    }
  fprintf(FdBalance, "\n");

#ifdef OUTPUT_CPU_CSV
  fprintf(FdCPUCSV, "\n");
#endif

#ifdef USE_SFR
  sprintf(buf, "%s%s", All.OutputDir, "sfr.txt");
  if(!(FdSfr = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }
#endif

#ifdef BLACKHOLE_POTMIN_DIAGNOSTIC
  sprintf(buf, "%s%s", All.OutputDir, "potmin_diag.txt");
  if(!(FdBHDiag = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }
#endif

#ifdef BLACK_HOLES
  sprintf(buf, "%s%s", All.OutputDir, "blackholes.txt");
  if(!(FdBlackHoles = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }
#endif



#ifdef DARKENERGY
  sprintf(buf, "%s%s", All.OutputDir, "darkenergy.txt");
  if(!(FdDE = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }
  else
    {
      if(RestartFlag == 0)
        {
          fprintf(FdDE, "nstep time H(a) ");
#ifndef TIMEDEPDE
          fprintf(FdDE, "w0 Omega_L ");
#else
          fprintf(FdDE, "w(a) Omega_L ");
#endif
          fprintf(FdDE, "\n");
          myflush(FdDE);
        }
    }
#endif
#ifdef OTVET
  sprintf(buf, "%s%s", All.OutputDir, "otvet.txt");
  if(!(FdOTVET = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }

#ifdef OTVET_MULTI_FREQUENCY
  sprintf(buf, "%s%s", All.OutputDir, "otvet_star_lum.txt");
  if(!(FdOTVETStar = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }
#endif
#endif

#ifdef SNE_FEEDBACK
  sprintf(buf, "%s%s", All.OutputDir, "supernovas.txt");
  if(!(FdSNe = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }
  if(RestartFlag == 0) {
#ifdef INJECT_TRACER_INTO_SN
      fprintf(FdSNe, "# time, sne_x, sne_y, sne_z, n, n_tracer, mass, volume, radius, energy\n");
#else
      fprintf(FdSNe, "# time,\t\t\t\tx, \t\t\t\t\t\ty, \t\t\t\t\t\tz, \t\t\t\t\t\tn, scheme, radius, mean_density, mean_temperature\n");
#endif
      myflush(FdSNe);
  }
#endif

#if defined(MRT) && defined(MRT_CHEMISTRY_PS2011)
  sprintf(buf, "%s%s", All.OutputDir, "mrt.txt");
  if(!(FdMRT = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }
#endif

#ifdef BINARYLOG
  sprintf(buf, "%s%s", All.OutputDir, "binary.txt");
  if(!(FdBinary = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }
#endif

#ifdef DG_CON_VARS_SUM_TO_FILE
  sprintf(buf, "%s%s", All.OutputDir, "angular_momentum_dg.txt");
  if(!(FdAngularMomentumDG = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }

  sprintf(buf, "%s%s", All.OutputDir, "mass_dg.txt");
  if(!(FdMassDG = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }

  sprintf(buf, "%s%s", All.OutputDir, "energy_dg.txt");
  if(!(FdEnergyDG = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }
#endif

#ifdef DG
  sprintf(buf, "%s%s", All.OutputDir, "info_dg.txt");
  if(!(FdInfoDG = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }
#endif

#ifdef NUCLEAR_NETWORK
  sprintf(buf, "%s%s", All.OutputDir, "composition.txt");
  if(!(FdNetwork = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }
#endif

#ifdef GENERAL_RELATIVITY
  sprintf(buf, "%s%s", All.OutputDir, "grlog.txt");
  if(!(FdGR = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }
#endif

#ifdef SINK_PARTICLES
#ifdef DUMP_SINK_PARTICLE_INFO
    if( open_sink_particle_evolution_file() == 0 )
      {
        sprintf(msg, "error in opening sink file... Too many already?");
        terminate(msg);
      }
#endif
#endif

#ifdef COSMIC_RAYS
  sprintf(buf, "%s%s", All.OutputDir, "crenergy.txt");
  if(!(FdCREnergy = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }
#endif

#if defined(DUST_LIVE) && defined(DL_PRODUCTION)
  sprintf(buf, "%s%s", All.OutputDir, "dust.txt");
  if(!(FdDust = fopen(buf, mode)))
    {
      sprintf(msg, "error in opening file '%s'\n", buf);
      terminate(msg);
    }
#endif
}


/*! \brief  Close the global log-files.
 */
void close_logfiles(void)
{
  if(ThisTask != 0)             /* only the root processors writes to the log files */
    return;

  fclose(FdCPU);
  fclose(FdInfo);
  fclose(FdEnergy);
  fclose(FdTimings);
  fclose(FdBalance);
  fclose(FdTimebin);

#ifdef OUTOUT_CPU_CSV
  fclose(FdCPUCSV);
#endif

#ifdef RT_ADVECT
  fclose(FdRad);
#endif

#ifdef USE_SFR
  fclose(FdSfr);
#endif

#ifdef GFM_STELLAR_EVOLUTION
  fclose(FdMetalsGas);
  fclose(FdMetalsStars);
  fclose(FdMetalsTot);
  fclose(FdSN);
#endif

#if defined(FM_STAR_FEEDBACK) && defined(OUTPUT_STELLAR_FEEDBACK)
  fclose(FdFeedback);
  fclose(FdGasFeed);
#endif

#ifdef DARKENERGY
  fclose(FdDE);
#endif

#ifdef OTVET
  fclose(FdOTVET);
#ifdef OTVET_MULTI_FREQUENCY
  fclose(FdOTVETStar);
#endif
#endif

#ifdef SNE_FEEDBACK
  fclose(FdSNe);
#endif

#ifdef BINARYLOG
  fclose(FdBinary);
#endif

#ifdef DG_CON_VARS_SUM_TO_FILE
  fclose(FdAngularMomentumDG);
  fclose(FdMassDG);
  fclose(FdEnergyDG);
#endif

#ifdef DG
  fclose(FdInfoDG);
#endif

#ifdef NUCLEAR_NETWORK
  fclose(FdNetwork);
#endif

#ifdef GENERAL_RELATIVITY
  fclose(FdGR);
#endif

#ifdef SINK_PARTICLES
#ifdef DUMP_SINK_PARTICLE_INFO
  fclose(FdSinkPart);
#endif
#endif

#ifdef COSMIC_RAYS
  fclose(FdCREnergy);
#endif

#if defined(DUST_LIVE) && defined(DL_PRODUCTION)
  fclose(FdDust);
#endif
}



/*! \brief Write the FdInfo and FdTimeBin files.
 *
 * At each time step this function writes on to two log-files.
 * In FdInfo, it just lists the timesteps that have been done, while in
 * FdTimeBin it outputs information about the active and occupied time-bins.
 */
void output_log_messages(void)
{
  double z;
  int i, j;
  double sum, avg_CPU_TimeBin[TIMEBINS], frac_CPU_TimeBin[TIMEBINS];
  int weight, corr_weight;
  long long tot_cumulative_grav[TIMEBINS], tot_cumulative_sph[TIMEBINS];
  long long tot_grav, tot_sph;

  TIMER_START(CPU_LOGS);

  report_detailed_memory_usage_of_largest_task();

  long long count[4 * TIMEBINS], tot_count[4 * TIMEBINS];
  long long *tot_count_grav = &tot_count[0], *tot_count_sph = &tot_count[TIMEBINS];
  int nelem = 2 * TIMEBINS;

  for(int i = 0; i < TIMEBINS; i++)
    count[i] = TimeBinsGravity.TimeBinCount[i];

  for(int i = 0; i < TIMEBINS; i++)
    count[i + TIMEBINS] = TimeBinsHydro.TimeBinCount[i];

#ifdef TRACER_PARTICLE
  for(int i=0; i < TIMEBINS; i++)
  count[i + 2 * TIMEBINS] = TimeBinsTracer.TimeBinCount[i];

  nelem += TIMEBINS;

  long long *tot_count_tracer= &tot_count[2 * TIMEBINS], tot_tracer;
#endif
#ifdef DUST_LIVE
  for(int i = 0; i < TIMEBINS; i++)
    count[i + nelem] = TimeBinsDust.TimeBinCount[i];

  long long *tot_count_dust = &tot_count[nelem], tot_dust;
  nelem += TIMEBINS;
#endif

  MPI_Reduce(count, tot_count, nelem, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      if(All.ComovingIntegrationOn)
        {
          z = 1.0 / (All.Time) - 1;
          fprintf(FdInfo, "\nSync-Point %d, TimeBin=%d, Time: %g, Redshift: %g, Systemstep: %g, Dloga: %g, Nsync-grv: %10llu, Nsync-hyd: %10llu\n",
                  All.NumCurrentTiStep, All.HighestActiveTimeBin, All.Time, z, All.TimeStep, log(All.Time) - log(All.Time - All.TimeStep), All.GlobalNSynchronizedGravity, All.GlobalNSynchronizedHydro);
          printf("\n\nSync-Point %d, Time: %g, Redshift: %g, Systemstep: %g, Dloga: %g, Nsync-grv: %10llu, Nsync-hyd: %10llu\n", All.NumCurrentTiStep,
                 All.Time, z, All.TimeStep, log(All.Time) - log(All.Time - All.TimeStep), All.GlobalNSynchronizedGravity, All.GlobalNSynchronizedHydro);
          fprintf(FdTimebin, "\nSync-Point %d, Time: %g, Redshift: %g, Systemstep: %g, Dloga: %g\n", All.NumCurrentTiStep, All.Time, z, All.TimeStep, log(All.Time) - log(All.Time - All.TimeStep));
          myflush(FdInfo);
        }
      else
        {
          fprintf(FdInfo, "\nSync-Point %d, TimeBin=%d, Time: %g, Systemstep: %g, Nsync-grv: %10llu, Nsync-hyd: %10llu\n",
                  All.NumCurrentTiStep, All.HighestActiveTimeBin, All.Time, All.TimeStep, All.GlobalNSynchronizedGravity, All.GlobalNSynchronizedHydro);
          printf("\n\nSync-Point %d, Time: %g, Systemstep: %g, Nsync-grv: %10llu, Nsync-hyd: %10llu\n", All.NumCurrentTiStep,
                 All.Time, All.TimeStep, All.GlobalNSynchronizedGravity, All.GlobalNSynchronizedHydro);
          fprintf(FdTimebin, "\nSync-Point %d, Time: %g, Systemstep: %g\n", All.NumCurrentTiStep, All.Time, All.TimeStep);
          myflush(FdInfo);
        }

      for(i = 1, tot_cumulative_grav[0] = tot_count_grav[0], tot_cumulative_sph[0] = tot_count_sph[0]; i < TIMEBINS; i++)
        {
          tot_cumulative_grav[i] = tot_count_grav[i] + tot_cumulative_grav[i - 1];
          tot_cumulative_sph[i] = tot_count_sph[i] + tot_cumulative_sph[i - 1];
        }

      for(i = 0; i < TIMEBINS; i++)
        {
          for(j = 0, sum = 0; j < All.CPU_TimeBinCountMeasurements[i]; j++)
            sum += All.CPU_TimeBinMeasurements[i][j];
          if(All.CPU_TimeBinCountMeasurements[i])
            avg_CPU_TimeBin[i] = sum / All.CPU_TimeBinCountMeasurements[i];
          else
            avg_CPU_TimeBin[i] = 0;
        }

      for(i = All.HighestOccupiedTimeBin, weight = 1, sum = 0; i >= 0 && tot_count_grav[i] > 0; i--, weight *= 2)
        {
          if(weight > 1)
            corr_weight = weight / 2;
          else
            corr_weight = weight;

          frac_CPU_TimeBin[i] = corr_weight * avg_CPU_TimeBin[i];
          sum += frac_CPU_TimeBin[i];
        }

      for(i = All.HighestOccupiedTimeBin; i >= 0 && tot_count_grav[i] > 0; i--)
        {
          if(sum)
            frac_CPU_TimeBin[i] /= sum;
        }

      char tracerString[13];
#ifdef TRACER_PARTICLE
      tot_tracer = 0;
      sprintf(tracerString, "tracer      ");
#else
      sprintf(tracerString, "%s", "");
#endif

      char dustString[13];
#ifdef DUST_LIVE
      tot_dust = 0;
      sprintf(dustString, "dust        ");
#else
      sprintf(dustString, "%s", "");
#endif

      fprintf(FdTimebin, "Occupied timebins: gravity      hydro     %s     %s     dt              cumul-grav   cumul-sph A D    avg-time  cpu-frac\n", tracerString, dustString);
      for(i = TIMEBINS - 1, tot_grav = tot_sph = 0; i >= 0; i--)
        {
          int binUsed = 0;

#if (defined(SELFGRAVITY) || defined(EXTERNALGRAVITY) || defined(EXACT_GRAVITY_FOR_PARTICLE_TYPE)) && !defined(MESHRELAX)
          if(tot_count_grav[i] > 0)
            binUsed = 1;
#endif

          if(tot_count_sph[i] > 0)
            binUsed = 1;

#ifdef TRACER_PARTICLE
          if(tot_count_tracer[i] > 0)
            binUsed = 1;

          sprintf(tracerString, "%10llu  ", tot_count_tracer[i]);
#else
          sprintf(tracerString, "%s", "");
#endif
#ifdef DUST_LIVE
          if(tot_count_dust[i] > 0)
            binUsed = 1;

          sprintf(dustString, "%10llu  ", tot_count_dust[i]);
#else
          sprintf(dustString, "%s", "");
#endif

          if(binUsed)
            {
              fprintf(FdTimebin,
                      " %c  bin=%2d      %10llu  %10llu  %s  %s  %16.12f       %10llu  %10llu %c %c  %10.2f    %5.1f%%\n",
                      TimeBinSynchronized[i] ? 'X' : ' ',
                      i, tot_count_grav[i], tot_count_sph[i], tracerString, dustString,
                      i > 0 ? (((integertime) 1) << i) * All.Timebase_interval : 0.0,
                      tot_cumulative_grav[i], tot_cumulative_sph[i],
                      (i == All.HighestActiveTimeBin) ? '<' : ' ',
                      (All.HighestActiveTimeBin >= All.SmallestTimeBinWithDomainDecomposition && i == All.HighestActiveTimeBin) ? '*' : ' ', avg_CPU_TimeBin[i], 100.0 * frac_CPU_TimeBin[i]);

              if(TimeBinSynchronized[i])
                {
                  tot_grav += tot_count_grav[i];
                  tot_sph += tot_count_sph[i];
#ifdef TRACER_PARTICLE
                  tot_tracer += tot_count_tracer[i];
#endif
#ifdef DUST_LIVE
                  tot_dust += tot_count_dust[i];
#endif
                }
            }
        }

#ifdef TRACER_PARTICLE
      fprintf(FdTimebin, "               ------------------------------------\n");
#else
      fprintf(FdTimebin, "               ------------------------\n");
#endif

#ifdef TRACER_PARTICLE
      sprintf(tracerString, "%10llu  ", tot_tracer);
#else
      sprintf(tracerString, "%s", "");
#endif

#ifdef DUST_LIVE
      sprintf(dustString, "%10llu  ", tot_dust);
#else
      sprintf(dustString, "%s", "");
#endif

#ifdef PMGRID
      if(All.PM_Ti_endstep == All.Ti_Current)
        {
          fprintf(FdTimebin, "PM-Step. Total: %10llu  %10llu  %s  %s\n", tot_grav, tot_sph, tracerString, dustString);
        }
      else
#endif
        {
          fprintf(FdTimebin, "Total active:   %10llu  %10llu  %s  %s\n", tot_grav, tot_sph, tracerString, dustString);
        }
      fprintf(FdTimebin, "\n");
      myflush(FdTimebin);

#ifdef CHEMISTRY
      printf("Abundances elec: %g, HM: %g, H2I: %g, H2II: %g\n", SphP[1].elec, SphP[1].HM, SphP[1].H2I, SphP[1].H2II);
#endif


#ifdef DARKENERGY
      if(All.ComovingIntegrationOn == 1)
        {
          double hubble_a;

          hubble_a = hubble_function(All.Time);
          fprintf(FdDE, "%d %g %e ", All.NumCurrentTiStep, All.Time, hubble_a);
#ifndef TIMEDEPDE
          fprintf(FdDE, "%e ", All.DarkEnergyParam);
#else
          fprintf(FdDE, "%e %e ", get_wa(All.Time), DarkEnergy_a(All.Time));
#endif
#ifdef TIMEDEPGRAV
          fprintf(FdDE, "%e %e", dHfak(All.Time), dGfak(All.Time));
#endif
          fprintf(FdDE, "\n");
          myflush(FdDE);
        }
#endif
    }

#if defined (VS_TURB) || defined (AB_TURB)
  log_turb_temp();
#endif

#ifdef RESTART_DEBUG
  log_restart_debug();
#endif

  TIMER_STOP(CPU_LOGS);
}

void init_cpu_log(void)
{
  int i = 0;

#define TIMER_STRUCT
#include "timer.h"


  for(i = 0; i < CPU_LAST; i++)
    {
      if(Timer_data[i].parent >= 0)
        Timer_data[i].depth = Timer_data[Timer_data[i].parent].depth + 1;
      else
        Timer_data[i].depth = 0;
    }

  for(i = 0; i < CPU_LAST; i++)
    {
      All.CPU_Sum[i] = 0.;
      CPU_Step[i] = 0.;
    }


  TimerStackPos = 0;
  TimerStack[0] = CPU_MISC;

  CPUThisRun = 0.;

#ifdef VTUNE_INSTRUMENT
  vtune_domain =__itt_domain_create("Task Domain");
#endif

  WallclockTime = second();
  StartOfRun = second();
}



/*! \brief Write the FdBalance and FdCPU files.
 *
 * At each time step this function writes on to two log-files.
 * In FdBalance, it outputs in a graphical way the amount of 
 * time spent in the various parts of the code, while
 * in FdCPU it writes information about the cpu-time consumption
 * of the various modules.
 */
void write_cpu_log(void)
{
  double max_CPU_Step[CPU_LAST], avg_CPU_Step[CPU_LAST], summed_CPU_Step[CPU_LAST];
  double t0, t1, tsum;
  double avg_total = 0;
  double local_total = 0;
  double max_total = 0;
  int i;


  TIMER_START(CPU_LOGS);

  for(i = 0; i < CPU_LAST; i++)
    {
      local_total += CPU_Step[i];
    }

  MPI_Reduce(CPU_Step, max_CPU_Step, CPU_LAST, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_total, &max_total, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(CPU_Step, avg_CPU_Step, CPU_LAST, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {

      /* sum up cpu items into groups */
      for(i = 0; i < CPU_LAST; i++)
        {
          summed_CPU_Step[i] = avg_CPU_Step[i];
        }
      for(i = CPU_LAST - 1; i > CPU_ALL; i--)
        {
          if(Timer_data[i].parent >= 0)
            {
              summed_CPU_Step[Timer_data[i].parent] += summed_CPU_Step[i];
            }
        }

      /* calc averages, update All.CPU_Sum */
      for(i = 0; i < CPU_LAST; i++)
        {
          avg_CPU_Step[i] /= NTask;
          avg_total += avg_CPU_Step[i];

          summed_CPU_Step[i] /= NTask;
          All.CPU_Sum[i] += summed_CPU_Step[i];
        }

      /* create balance.txt string */
      char cpu_String[CPU_STRING_LEN + 1];
      put_symbol(cpu_String, 0., 1.0, '-');

      for(i = 1, tsum = 0.0; i < CPU_LAST; i++)
        {
          if(max_CPU_Step[i] > 0 && Timer_data[i].symb != 0 && Timer_data[i].symbImbal != 0)
            {
              t0 = tsum;
              t1 = tsum + avg_CPU_Step[i] * (avg_CPU_Step[i] / max_CPU_Step[i]);
              put_symbol(cpu_String, t0 / avg_total, t1 / avg_total, Timer_data[i].symb);
              tsum += t1 - t0;

              t0 = tsum;
              t1 = tsum + avg_CPU_Step[i] * ((max_CPU_Step[i] - avg_CPU_Step[i]) / max_CPU_Step[i]);
              put_symbol(cpu_String, t0 / avg_total, t1 / avg_total, Timer_data[i].symbImbal);
              tsum += t1 - t0;
            }
        }
      //put_symbol(cpu_String, tsum / max_total, 1.0, '-');
      fprintf(FdBalance, "Step=%7d  sec=%10.3f Nsync-grv=%10llu Nsync-hyd=%10llu  %s\n", All.NumCurrentTiStep, max_total, All.GlobalNSynchronizedGravity, All.GlobalNSynchronizedHydro, cpu_String);
      myflush(FdBalance);



      if(All.CPU_TimeBinCountMeasurements[All.HighestActiveTimeBin] == NUMBER_OF_MEASUREMENTS_TO_RECORD)
        {
          All.CPU_TimeBinCountMeasurements[All.HighestActiveTimeBin]--;
          memmove(&All.CPU_TimeBinMeasurements[All.HighestActiveTimeBin][0], &All.CPU_TimeBinMeasurements[All.HighestActiveTimeBin][1], (NUMBER_OF_MEASUREMENTS_TO_RECORD - 1) * sizeof(double));
        }

      All.CPU_TimeBinMeasurements[All.HighestActiveTimeBin][All.CPU_TimeBinCountMeasurements[All.HighestActiveTimeBin]++] = max_total;


#ifdef OUTPUT_CPU_CSV
      fprintf(FdCPUCSV, "%d, %g, %d, %d, %d, ", All.NumCurrentTiStep, All.Time, NTask, All.MultipleDomains, All.HighestActiveTimeBin);
#endif

      fprintf(FdCPU, "Step %d, Time: %g, CPUs: %d, MultiDomains: %d, HighestActiveTimeBin: %d\n", All.NumCurrentTiStep, All.Time, NTask, All.MultipleDomains, All.HighestActiveTimeBin);
      fprintf(FdCPU, "                          diff               cumulative\n");

      for(i = 0; i < CPU_LAST; i++)
        {
          fprintf(FdCPU, "%*s%*s%10.2f  %5.1f%% %10.2f  %*s%5.1f%%\n", 2 * Timer_data[i].depth, "",
                  -20 + 2 * Timer_data[i].depth,
                  Timer_data[i].longname, summed_CPU_Step[i],
                  summed_CPU_Step[i] / summed_CPU_Step[CPU_ALL] * 100., All.CPU_Sum[i], 5 * Timer_data[i].depth, "", All.CPU_Sum[i] / All.CPU_Sum[CPU_ALL] * 100.);
#ifdef OUTPUT_CPU_CSV
          fprintf(FdCPUCSV, "%f, %f, %f, ", summed_CPU_Step[i], All.CPU_Sum[i], All.CPU_Sum[i] / All.CPU_Sum[CPU_ALL] * 100.);
#endif

        }

      fprintf(FdCPU, "\n");
      myflush(FdCPU);

#ifdef OUTPUT_CPU_CSV
      fprintf(FdCPUCSV, "\n");
      myflush(FdCPUCSV);
#endif

    }

  for(i = 0; i < CPU_LAST; i++)
    CPU_Step[i] = 0.;

  CPUThisRun = timediff(StartOfRun, second());

  TIMER_STOP(CPU_LOGS);
}


/*! \brief Fill the cpu balance string representing the cpu usage in a graphical way
 * 
 * This function fills a fraction, specified by the parameters t0 and t1,
 * of the array string with the debug symbol given by c.
 *
 * \param string string to fill
 * \param t0 initial position of the symbol in the array as a fraction of its maximum dimension 
 * \param t1 final position of the symbol in the array as a fraction of its maximum dimension 
 * \param c symbol to be put on string
 */
void put_symbol(char *string, double t0, double t1, char c)
{
  int i, j;

  i = (int) (t0 * CPU_STRING_LEN + 0.5);
  j = (int) (t1 * CPU_STRING_LEN);

  if(i < 0)
    i = 0;
  if(j >= CPU_STRING_LEN)
    j = CPU_STRING_LEN;

  while(i <= j)
    string[i++] = c;

  string[CPU_STRING_LEN] = 0;
}
