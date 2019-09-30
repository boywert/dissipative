/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/windtunnel.c
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"
#include "voronoi.h"


#if defined(WINDTUNNEL) && defined(WINDTUNNEL_EXTERNAL_SOURCE)


static int NWindTable;

static struct wind_table
{
  double time, rho, vel;
} *WindTable;


void interpolate_from_wind_table(double t, double *rho, double *vel)
{
  int binlow, binhigh, binmid;

  if(t < WindTable[0].time || t > WindTable[NWindTable - 1].time)
    {
      terminate("time outside of wind interpolation table");
    }

  binlow = 0;
  binhigh = NWindTable - 1;

  while(binhigh - binlow > 1)
    {
      binmid = (binhigh + binlow) / 2;
      if(t < WindTable[binmid].time)
        binhigh = binmid;
      else
        binlow = binmid;
    }

  double dt = WindTable[binhigh].time - WindTable[binlow].time;

  if(dt == 0)
    terminate("dt=0");

  double u = (t - WindTable[binlow].time) / dt;

  *rho = (1 - u) * WindTable[binlow].rho + u * WindTable[binhigh].rho;
  *vel = (1 - u) * WindTable[binlow].vel + u * WindTable[binhigh].vel;
}



void read_windtunnel_file(void)
{
  FILE *fd;
  double k, p;

  if(!(fd = fopen(All.WindTunnelExternalSourceFile, "r")))
    {
      char buf[1000];
      sprintf(buf, "can't read file '%s' with windtunnel data on task %d\n", buf, ThisTask);
      terminate(buf);
    }

  NWindTable = 0;
  do
    {
      double t, rho, v;
      if(fscanf(fd, " %lg %lg %lg ", &t, &rho, &v) == 3)
        NWindTable++;
      else
        break;
    }
  while(1);

  fclose(fd);

  mpi_printf("found %d rows in input wind table\n", NWindTable);


  WindTable = mymalloc("WindTable", NWindTable * sizeof(struct wind_table));

  if(!(fd = fopen(All.WindTunnelExternalSourceFile, "r")))
    {
      char buf[1000];
      sprintf(buf, "can't read file '%s' with windtunnel data on task %d\n", buf, ThisTask);
      terminate(buf);
    }


  NWindTable = 0;
  do
    {
      double t, rho, v;
      if(fscanf(fd, " %lg %lg %lg ", &t, &rho, &v) == 3)
        {
          WindTable[NWindTable].time = t;
          WindTable[NWindTable].rho = rho;
          WindTable[NWindTable].vel = v;
          NWindTable++;
        }
      else
        break;
    }
  while(1);

  fclose(fd);

  /* note: we'll assume that this file is sorted according to time */
}


#endif
