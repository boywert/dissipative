/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/sinks/sinks_aux.c
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


void sinks_set_constants(void)
{
  set_cosmo_factors_for_current_time();

  SKD.DistFac = pow(All.cf_atime / All.HubbleParam, 2);

  SKD.NHFac = HYDROGEN_MASSFRAC * All.UnitDensity_in_cgs * All.cf_a3inv * pow(All.HubbleParam, 2) / PROTONMASS;

#ifdef SGCHEM
  SKD.DistFac = 1.0;
  SKD.NHFac = HYDROGEN_MASSFRAC * All.UnitDensity_in_cgs / PROTONMASS;

  mpi_printf("Set sink constants %d %g %g\n", ThisTask, SKD.NHFac, All.UnitDensity_in_cgs / PROTONMASS);

#endif


}


void write_sink_data(int xaxis, int yaxis, int zaxis, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
{
  int i;
  char buf[500];
  FILE *fd = 0;
  sprintf(buf, "%s/sinks_%03d", All.OutputDir, RestartSnapNum);

  fd = open_file(buf);

  int num_sinks = 0;
  double sink_pos[2];

  for(i = 0; i < SKD.TotNumSinks; i++)
    {
      if(SKD.SINK[i].Pos[xaxis] >= xmin && SKD.SINK[i].Pos[xaxis] < xmax
         && SKD.SINK[i].Pos[yaxis] >= ymin && SKD.SINK[i].Pos[yaxis] < ymax && SKD.SINK[i].Pos[zaxis] >= zmin && SKD.SINK[i].Pos[zaxis] < zmax)
        num_sinks++;
    }

  my_fwrite(&num_sinks, sizeof(int), 1, fd);

  for(i = 0; i < SKD.TotNumSinks; i++)
    {
      if(SKD.SINK[i].Pos[xaxis] >= xmin && SKD.SINK[i].Pos[xaxis] < xmax
         && SKD.SINK[i].Pos[yaxis] >= ymin && SKD.SINK[i].Pos[yaxis] < ymax && SKD.SINK[i].Pos[zaxis] >= zmin && SKD.SINK[i].Pos[zaxis] < zmax)
        {
          sink_pos[0] = (SKD.SINK[i].Pos[xaxis] - xmin) / (xmax - xmin);
          sink_pos[1] = (SKD.SINK[i].Pos[yaxis] - ymin) / (ymax - ymin);

          my_fwrite(&sink_pos, sizeof(double), 2, fd);
          my_fwrite(&SKD.SINK[i].Mass, sizeof(double), 1, fd);
        }
    }

  fclose(fd);
}
