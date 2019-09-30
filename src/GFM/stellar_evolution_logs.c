/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/GFM/stellar_evolution_logs.c
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


#ifdef GFM_STELLAR_EVOLUTION

#if defined(FM_STAR_FEEDBACK) && defined(OUTPUT_STELLAR_FEEDBACK)
/* global variables for feedback log files */
static MyFloat global_tot_feed_egy = 0.0;
static MyFloat global_th_feed_egy = 0.0;
static MyFloat global_kin_feed_egy = 0.0;
static MyFloat global_SNII_Num = 0.0;

static MyFloat feed_egy_converted_cells = 0.0;
static MyFloat th_feed_egy_converted_cells = 0.0;
static MyFloat kin_feed_egy_converted_cells = 0.0;
#endif

void output_stellar_evolution_statistics(void)
{
  int i, j;

  MyFloat metal_mass_gas[GFM_N_CHEM_ELEMENTS + 1], metal_mass_stars[GFM_N_CHEM_ELEMENTS + 1], metal_mass_total[GFM_N_CHEM_ELEMENTS + 1], global_metal_mass[GFM_N_CHEM_ELEMENTS + 1];
  MyFloat SNIaRate = 0.0, global_SNIaRate = 0.0;
  MyFloat SNIIRate = 0.0, global_SNIIRate = 0.0;

  /* metal data */
  for(i = 0; i < GFM_N_CHEM_ELEMENTS + 1; i++)
    {
      metal_mass_gas[i] = 0.0;
      metal_mass_stars[i] = 0.0;
    }

  for(j = 0; j < NumPart; j++)
    {
      if(P[j].Type == 0)        /* gas particle */
        {
          for(i = 0; i < GFM_N_CHEM_ELEMENTS; i++)
            metal_mass_gas[i] += SphP[j].MassMetals[i];

          metal_mass_gas[GFM_N_CHEM_ELEMENTS] += P[j].Mass;
        }
      if(P[j].Type == 4 && P[j].Mass > 0 && STP(j).BirthTime > 0)       /* star particle */
        {
          for(i = 0; i < GFM_N_CHEM_ELEMENTS; i++)
            metal_mass_stars[i] += STP(j).MassMetals[i];

          metal_mass_stars[GFM_N_CHEM_ELEMENTS] += P[j].Mass;
        }
    }

  for(i = 0; i < GFM_N_CHEM_ELEMENTS; i++)
    metal_mass_total[i] = metal_mass_gas[i] + metal_mass_stars[i];

  metal_mass_total[GFM_N_CHEM_ELEMENTS] = metal_mass_gas[GFM_N_CHEM_ELEMENTS] + metal_mass_stars[GFM_N_CHEM_ELEMENTS];

  /* output total metal mass in gas particles */
  MPI_Reduce(metal_mass_gas, global_metal_mass, GFM_N_CHEM_ELEMENTS + 1, MPI_MYFLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      fprintf(FdMetalsGas, "%e", All.Time);
      for(i = 0; i < GFM_N_CHEM_ELEMENTS + 1; i++)
        fprintf(FdMetalsGas, " %e", global_metal_mass[i]);
      fprintf(FdMetalsGas, "\n");
      myflush(FdMetalsGas);
    }

  /* output total metal mass in star particles */
  MPI_Reduce(metal_mass_stars, global_metal_mass, GFM_N_CHEM_ELEMENTS + 1, MPI_MYFLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      fprintf(FdMetalsStars, "%e", All.Time);
      for(i = 0; i < GFM_N_CHEM_ELEMENTS + 1; i++)
        fprintf(FdMetalsStars, " %e", global_metal_mass[i]);
      fprintf(FdMetalsStars, "\n");
      myflush(FdMetalsStars);
    }

  /* output total metal mass in all particles */
  MPI_Reduce(metal_mass_total, global_metal_mass, GFM_N_CHEM_ELEMENTS + 1, MPI_MYFLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      fprintf(FdMetalsTot, "%e", All.Time);
      for(i = 0; i < GFM_N_CHEM_ELEMENTS + 1; i++)
        fprintf(FdMetalsTot, " %e", global_metal_mass[i]);
      fprintf(FdMetalsTot, "\n");
      myflush(FdMetalsTot);
    }


  /* SN data */
  for(i = 0; i < N_star; i++)
    {
      SNIaRate += StarP[i].SNIaRate;
      SNIIRate += StarP[i].SNIIRate;
    }

  MPI_Reduce(&SNIaRate, &global_SNIaRate, 1, MPI_MYFLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&SNIIRate, &global_SNIIRate, 1, MPI_MYFLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      fprintf(FdSN, "%e %e %e\n", All.Time, global_SNIaRate, global_SNIIRate);
      myflush(FdSN);
    }
}


#if defined(FM_STAR_FEEDBACK) && defined(OUTPUT_STELLAR_FEEDBACK)
void output_stellar_feedback_statistics(void)
{
  mpi_printf("GFM_STELLAR_EVOLUTION: Writing feedback statistics.\n");

  int i;
  /* star variables */
  MyFloat tot_feed_egy_time_step = 0.0, global_tot_feed_egy_time_step = 0.0;
  MyFloat th_feed_egy_time_step = 0.0, global_th_feed_egy_time_step = 0.0;
  MyFloat kin_feed_egy_time_step = 0.0, global_kin_feed_egy_time_step = 0.0;
  MyFloat SNII_Num_time_step = 0.0, global_SNII_Num_time_step = 0.0;

  /* gas variables */
  MyFloat tot_gas_feed_egy = 0.0, global_tot_gas_feed_egy = 0.0;
  MyFloat th_gas_feed_egy = 0.0, global_th_gas_feed_egy = 0.0;
  MyFloat kin_gas_feed_egy = 0.0, global_kin_gas_feed_egy = 0.0;

  /* stellar feedback data (only active star particles are counted) */
  for(i = 0; i < Nstar; i++)
    {
      tot_feed_egy_time_step += STP(StarParticle[i].index).FeedbackEnergy;
      th_feed_egy_time_step += STP(StarParticle[i].index).FeedbackThEnergy;
      kin_feed_egy_time_step += STP(StarParticle[i].index).FeedbackKinEnergy;
      SNII_Num_time_step += STP(StarParticle[i].index).SNII_Num;
    }

  MPI_Reduce(&tot_feed_egy_time_step, &global_tot_feed_egy_time_step, 1, MPI_MYFLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&th_feed_egy_time_step, &global_th_feed_egy_time_step, 1, MPI_MYFLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&kin_feed_egy_time_step, &global_kin_feed_egy_time_step, 1, MPI_MYFLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&SNII_Num_time_step, &global_SNII_Num_time_step, 1, MPI_MYFLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      global_tot_feed_egy += global_tot_feed_egy_time_step;
      global_th_feed_egy += global_th_feed_egy_time_step;
      global_kin_feed_egy += global_kin_feed_egy_time_step;
      global_SNII_Num += global_SNII_Num_time_step;

      fprintf(FdFeedback, "%e %e %e %e %e %e %e %e %e\n", All.Time, global_SNII_Num_time_step, global_SNII_Num,
              global_tot_feed_egy_time_step, global_th_feed_egy_time_step, global_kin_feed_egy_time_step, global_tot_feed_egy, global_th_feed_egy, global_kin_feed_egy);
      myflush(FdFeedback);
    }

  /* cumulative feedback energy received by gas cells (all cells, included those converted into stars) */
  for(i = 0; i < NumGas; i++)
    {

      /* do not count cells converted into stars (this is done after this loop) */
      if(P[i].Type != 0)
        continue;

      /* do not count derefined cells (their energy is already assigned to neighbours) */
      if(P[i].Mass > 0 && P[i].ID != 0)
        {
          tot_gas_feed_egy += SphP[i].TotEgyFeed;
          th_gas_feed_egy += SphP[i].IntEgyFeed;
          kin_gas_feed_egy += SphP[i].KinEgyFeed;
        }
    }

  /* add the energy of gas cells converted into stars */
  tot_gas_feed_egy += feed_egy_converted_cells;
  th_gas_feed_egy += th_feed_egy_converted_cells;
  kin_gas_feed_egy += kin_feed_egy_converted_cells;

  MPI_Reduce(&tot_gas_feed_egy, &global_tot_gas_feed_egy, 1, MPI_MYFLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&th_gas_feed_egy, &global_th_gas_feed_egy, 1, MPI_MYFLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&kin_gas_feed_egy, &global_kin_gas_feed_egy, 1, MPI_MYFLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      fprintf(FdGasFeed, "%e %e %e %e\n", All.Time, global_tot_gas_feed_egy, global_th_gas_feed_egy, global_kin_gas_feed_egy);

      myflush(FdGasFeed);
    }
}

/* takes into account the feedback energy of cells that are converted into star particles */
void update_cumulative_feedback_energy_of_converted_cells(int index)
{
  feed_egy_converted_cells += SphP[index].TotEgyFeed;
  th_feed_egy_converted_cells += SphP[index].IntEgyFeed;
  kin_feed_egy_converted_cells += SphP[index].KinEgyFeed;
}

#endif

#endif
