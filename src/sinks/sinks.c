/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/sinks/sinks.c
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


void sinks(void)
{
  double t0, t1, dt, dt_max;

  CPU_Step[CPU_MISC] += measure_time();

  if(All.Time == All.TimeBegin || SKD.AccRad == 0)
    return;

  if(SKD.NHThresh)
    {
      t0 = second();

      sinks_set_constants();

      sinks_get_num_sinks();
      mpi_printf("SINKS: Got sink numbers %d\n", SKD.TotNumSinks);

      sinks_begin();

      sinks_get_active_sinks();
      mpi_printf("SINKS: Got active sink numbers %d\n", TotNSinks);

      sinks_dmass();

      sinks_accrete(); // Debugging, 2016

      sinks_end();

#ifndef SGCHEM
      if(!SKD.TotNumSinks)
#endif
        sinks_create();

      t1 = second();

      dt = timediff(t0, t1);

      MPI_Allreduce(&dt, &dt_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

      mpi_printf("SINKS: Done! Took %g seconds\n", dt_max);
    }

  CPU_Step[CPU_SINKS] += measure_time();

  //endrun();
}
