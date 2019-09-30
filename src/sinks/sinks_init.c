/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/sinks/sinks_init.c
 * \date        01/2013
 * \author      Thomas Greif
 * \brief       Sink particles
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 */

#include "../allvars.h"
#include "../proto.h"

int NSinks, TotNSinks;
struct SinksAux_struct *SinksAux;

void sinks_begrun(void)
{
  CPU_Step[CPU_MISC] += measure_time();

  MPI_Bcast(&SKD, sizeof(struct SKD_struct), MPI_BYTE, 0, MPI_COMM_WORLD);

  SKD.Flag = 0;
  SKD.AccRad2 = pow(SKD.AccRad, 2);

  CPU_Step[CPU_SINKS] += measure_time();
}


void sinks_get_num_sinks(void)
{
  int i; 
  SKD.Flag = 0;
  SKD.NumSinks = 0;

  if(!SKD.Flag)
    {
      for(i = 0; i < NumPart; i++)
        if(P[i].Type == 5 && P[i].Mass > 0)
          SKD.NumSinks++;

      MPI_Allreduce(&SKD.NumSinks, &SKD.TotNumSinks, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    }

  SKD.Flag = 1;
}

void sinks_get_active_sinks()
{
  int i;
  NSinks = 0;

  for(int idx = 0; idx < TimeBinsSinksAccretion.NActiveParticles; idx++)
    {
      i = TimeBinsSinksAccretion.ActiveParticleList[idx];
      if(i < 0)
        continue;

     if(P[i].Ti_Current != All.Ti_Current)
        terminate("how can this be?");

      if(P[i].Type == 5 && P[i].Mass > 0)
        {
          SinksAux[NSinks].SinksAuxID = i;
          SinksAux[NSinks].MinTimeBin = TIMEBINS;
          NSinks++;
        }
    }

  if(NSinks > SKD.NumSinks)
    terminate("Number of active sinks %d larger than total sink number %d", NSinks, SKD.NumSinks);

  MPI_Allreduce(&NSinks, &TotNSinks, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
}

void sinks_begin()
{
  SinksAux = mymalloc("SinksAux", SKD.NumSinks * sizeof(struct SinksAux_struct));
}

void sinks_end()
{
  myfree(SinksAux);
}

