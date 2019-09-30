/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        
 * \date        10/2014
 * \author      Christine Simpson
 * \brief
 * \details     Create star particles from cell sfr
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <math.h>
#include "../allvars.h"
#include "../proto.h"

#if defined(LOCAL_FEEDBACK) && defined(USE_SFR) 
#if defined(EXTERNALSHEARBOX_KSRATE_RANDOM) || defined(EXTERNALSHEARBOX_MIXED_INJECTION)

void create_dummy_particles(void)
{
  CPU_Step[CPU_MISC] += measure_time();
    
  double SNRate = All.LocalFeedbackSNRate;
  double dt = (All.HighestActiveTimeBin ? (((int) 1) << All.HighestActiveTimeBin) : 0) * All.Timebase_interval;
  if(dt == 0.0)
    return;
  
  double fg = All.ShearBoxFg;
  double mu = All.ShearBoxMu;

#ifdef EXTERNALSHEARBOX_KSRATE_UPDATE_PARAM
  double Sigma = All.ShearBoxSigmaNow;
  double bscale = All.ShearBoxB;
#else
  double Sigma = All.ShearBoxSigma0;
  double bscale = 61.0 * (fg / 0.1 / mu) / (Sigma / 10.0);
#endif

#ifdef EXTERNALSHEARBOX_MIXED_INJECTION
  double prob_factor = (1.0/EXTERNALSHEARBOX_MIXED_INJECTION - 1.0);

  int idx,i;
  double sum_sfr = 0.0;
  /* loop over all active cells and particles */
    
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      /* just consider gas cells */
      if(P[i].Type == 0)
        {
          /* skip cells that have been swallowed or eliminated */
          if(P[i].Mass == 0 && P[i].ID == 0)
            continue;
          sum_sfr += SphP[i].Sfr;
	}
    }
  double totpeaksfrrate;
  MPI_Allreduce(&sum_sfr, &totpeaksfrrate, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  double prob = prob_factor*totpeaksfrrate* dt * SNRate / 100.0;
  mpi_printf("CMS_FEEDBACK: create particles: time = %f prob = %f prob_factor = %f sfr = %g dt = %g %d %d\n",All.Time,prob,prob_factor,totpeaksfrrate,dt,All.HighestActiveTimeBin,All.HighestOccupiedTimeBin);
#else /* not mixed injection */

  double sigma_star = 2.5e-4 * pow(Sigma, 1.4) * (All.UnitTime_in_s / 3.15569e7) * (1.9891e33 / All.UnitMass_in_g) * (All.UnitLength_in_cm / 3.0857e21) * (All.UnitLength_in_cm / 3.0857e21);   /*KS relation in code units */
  double totsfrrate = sigma_star * boxSize_X * boxSize_Y;
  double total_sm = totsfrrate * dt;

  double prob = sigma_star * boxSize_X * boxSize_Y * dt * SNRate / 100.0;

#endif

  int N_SN = 0;

  if(All.MaxID == 0)            /* MaxID not calculated yet */
    calculate_maxid();

  if(ThisTask == 0)
    {
      while(prob > 1)
        {
          N_SN += 1;
          prob -= 1.0;
        }
      double rand = get_random_number();
      if(rand < prob)
        N_SN += 1;

      double xpos, ypos, zpos, F;
      int istar, isn;
      for(isn = 0; isn < N_SN; isn++)
        {
          xpos = boxSize_X * get_random_number();
          ypos = boxSize_Y * get_random_number();
          F = get_random_number();
          zpos = bscale * atanh(2.0 * F - 1.0) + boxSize_Z / 2.0;

          istar = NumPart + isn;
          P[istar] = P[0];
          P[istar].Type = 4;
          P[istar].Pos[0] = xpos;
          P[istar].Pos[1] = ypos;
          P[istar].Pos[2] = zpos;

          P[istar].Mass = 0.0;
          P[istar].SofteningType = All.SofteningTypeOfPartType[P[istar].Type];
          int k;
          for(k = 0; k < 3; k++)
            P[istar].Vel[k] = 0.0;

          P[istar].ID = All.MaxID + isn + 1;
          timebin_add_particle(&TimeBinsGravity, istar, 0, P[istar].TimeBinGrav, TimeBinSynchronized[P[istar].TimeBinGrav]);
        }
      NumPart += N_SN;
    }

  int tot_N_SN;
  MPI_Allreduce(&N_SN, &tot_N_SN, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  mpi_printf("CMS_FEEDBACK: Creating dummy particles %d\n", tot_N_SN);

  All.TotNumPart += tot_N_SN;
  All.MaxID += tot_N_SN;

#ifndef  EXTERNALSHEARBOX_MIXED_INJECTION
  /* Print to sfr file */
  if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
    {                           /* only do this for full steps */
      double total_sum_mass_stars = tot_N_SN * 100.0 / SNRate;
      double rate_in_msunperyear = totsfrrate * (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);
      double cum_mass_stars = -1.0;
      if(ThisTask == 0)
        {
          fprintf(FdSfr, "%14e %14e %14e %14e %14e %14e\n", All.Time, total_sm, totsfrrate, rate_in_msunperyear, total_sum_mass_stars, cum_mass_stars);
          myflush(FdSfr);
        }
    }
#endif
}

#endif
#endif

#if defined(LOCAL_FEEDBACK) && defined(USE_SFR) && !defined(EXTERNALSHEARBOX_KSRATE_RANDOM)

void create_star_particles(void)
{
  CPU_Step[CPU_MISC] += measure_time();
  
  int i;

  double sum_sfr = 0.0;
  double sum_sm = 0.0;
  double sum_star_mass = 0.0;

  double dt, prob, rand, deltaM;
  int count_spawned = 0, count_converted = 0, count_feedback = 0;
  double MinimumStarMass;
  int idx;
  /* loop over all active cells and particles */
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      /* just consider gas cells */
      if(P[i].Type == 0)
        {
          /* skip cells that have been swallowed or eliminated */
          if(P[i].Mass == 0 && P[i].ID == 0)
            continue;

          dt = (All.HighestActiveTimeBin ? (((int) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;

	  double factor = 1.0;
#ifdef EXTERNALSHEARBOX_MIXED_INJECTION
	  factor =  1.0/EXTERNALSHEARBOX_MIXED_INJECTION;
#endif
          deltaM = SphP[i].Sfr *factor* dt;

          sum_sfr += SphP[i].Sfr*factor;
          sum_sm += deltaM;



#ifdef LOCAL_FEEDBACK_PARTICLES
	  if(P[i].Mass < 2.0*All.TargetGasMass)
	    MinimumStarMass = P[i].Mass;
	  else
	    MinimumStarMass = All.TargetGasMass;

	  prob = deltaM/MinimumStarMass;
	  rand = get_random_number();

    	  if(rand < prob)
	    {
	      if(P[i].Mass == MinimumStarMass){
		  convert_cell_into_star(i,All.Time);
		  timebin_remove_particle(&TimeBinsHydro, idx, P[i].TimeBinHydro);
		  Stars_converted++;
		  count_converted++;
		  sum_star_mass += P[i].Mass;
		}else{
		  spawn_star_from_cell(i,All.Time,NumPart+count_spawned,MinimumStarMass);
		  count_spawned++;
		  sum_star_mass += MinimumStarMass;
		} 
	    }
#endif
        }
    }

#ifdef LOCAL_FEEDBACK_PARTICLES
  /* Calculate ids of spawned particles and print out number of created particles */
  int countall_spawned, countall_converted;
  MPI_Allreduce(&count_spawned, &countall_spawned, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&count_converted, &countall_converted, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(countall_spawned > 0)
    {
      int *list;

      if(All.MaxID == 0)        /* MaxID not calculated yet */
        calculate_maxid();

      list = mymalloc("list", NTask * sizeof(int));

      MPI_Allgather(&count_spawned, 1, MPI_INT, list, 1, MPI_INT, MPI_COMM_WORLD);

      MyIDType newid = All.MaxID + 1;

      for(i = 0; i < ThisTask; i++)
        newid += list[i];

      myfree(list);

      for(i = 0; i < count_spawned; i++)
        {
          P[NumPart + i].ID = newid++;
	  //          printf("LOCAL_FEEDBACK: spawning star from cell; Star: time = %f; id = %d; mass = %f; height = %f\n", All.Time, P[NumPart + i].ID, P[NumPart + i].Mass, abs(P[NumPart + i].Pos[2] - 5000.0));
        }
      All.MaxID += countall_spawned;

    }

  if(countall_spawned > 0 || countall_converted > 0)
    {
      All.TotNumPart += countall_spawned;
      All.TotNumGas -= countall_converted;
      NumPart += count_spawned;
    }

  int countall = countall_converted + countall_spawned;

  if(countall > 0)
    {
      //    ngb_recompute_nodes();
      mpi_printf("LOCAL_FEEDBACK: Created %d particles at time %f\n", countall, All.Time);
    }
#endif

  /* Print to sfr file */
  if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
    {                           /* only do this for full steps */
      double totsfrrate, total_sm, total_sum_mass_stars;
      MPI_Allreduce(&sum_sfr, &totsfrrate, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Reduce(&sum_sm, &total_sm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#ifdef LOCAL_FEEDBACK_PARTICLES
      MPI_Reduce(&sum_star_mass, &total_sum_mass_stars, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
      total_sum_mass_stars = 0.0;
#endif
      double rate_in_msunperyear = totsfrrate * (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);
      double cum_mass_stars = -1.0;
      if(ThisTask == 0)
        {
          fprintf(FdSfr, "%14e %14e %14e %14e %14e %14e\n", All.Time, total_sm, totsfrrate, rate_in_msunperyear, total_sum_mass_stars, cum_mass_stars);
          myflush(FdSfr);
        }
    }
}
#endif
