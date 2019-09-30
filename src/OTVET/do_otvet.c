/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/OTVET/do_otvet.c
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


#include "../allvars.h"
#include "otvet_proto.h"


#ifdef OTVET

int do_otvet(void)
{
  int i, no, p;
  /* do only for highest time bin */
  if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
    {

      double timeeach = 0, timeall = 0, tstart = 0, tend = 0;
      All.otvet_Radiation_Ti_endstep = All.Ti_Current;

      for(i = 0; i < NumGas; i++)
        {
          SphP[i].Hsml = SphP[i].MaxDelaunayRadius;
        }
      density();                /* this changes SphP.Hsml value to the classical SPH value for gas */

      /* WARNING: this could change P[i].Soft of 
         star particles to their Hsml (like in Gadget). Not yet implemented this way */
      mpi_printf("OTVET: done with density\n");

      otvet_eddington();
      mpi_printf("OTVET: done with Eddington\n");


#if defined(OTVET_SCATTER_SOURCE) && defined(EDDINGTON_TENSOR_STARS) && !defined(GFM)
      otvet_stardensity();
      mpi_printf("OTVET: done with star_density\n");
#endif

      otvet_star_lum();
      mpi_printf("OTVET: done with star_lum\n");

      /***** evolve the transport of radiation *****/
      mpi_printf("OTVET: start radtransfer...\n");

      tstart = second();

      otvet_radtransfer();
      mpi_printf("OTVET: done with CG method\n");

      tend = second();
      timeeach = timediff(tstart, tend);
      MPI_Allreduce(&timeeach, &timeall, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#ifndef OTVET_SILENT
      mpi_printf("OTVET: time consumed is %g \n", timeall);
      mpi_printf("OTVET: done with radtransfer! \n");
#endif

      All.otvet_Radiation_Ti_begstep = All.otvet_Radiation_Ti_endstep;
    }


  /* chemistry updated at sub-stepping as well */
  otvet_update_chemistry();
  mpi_printf("OTVET: done with update_chem\n");

  /* radiative pressure */
#ifdef RADPRESS_OPT_THICK
  mpi_printf("RADPRESS_OPT_THICK: starting...\n");
  radpressthick();
  mpi_printf("RADPRESS_OPT_THICK: done.\n");
#endif

#ifdef RADPRESS_OPT_THIN
  mpi_printf("RADPRESS_OPT_THIN: starting...\n");
  radpressthin();
  mpi_printf("RADPRESS_OPT_THIN: done.\n");
#endif

#ifdef OTVET_CHEMISTRY_PS2009
  if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
    otvet_write_stats();
#endif

  return (0);
}
#endif
