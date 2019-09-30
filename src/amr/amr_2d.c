/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/amr/amr_2d.c
 * \date        10/2015
 * \author      Shy Genel
 * \brief       image dump based on AMR mesh
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
#include <gmp.h>

#include "../allvars.h"
#include "../proto.h"
#include "amr.h"


#if defined(AMR) && defined(TWODIMS) && !defined(ONEDIMS) && !defined(LONG_X) && !defined(LONG_Y)       /* will only be compiled in 2D case */

#if defined(VORONOI_FIELD_DUMP_PIXELS_X) && defined(VORONOI_FIELD_DUMP_PIXELS_Y)
void do_special_dump(int num, int gradients_flag)
{
  CPU_Step[CPU_MISC] += measure_time();

  char buf[1000];
  int pixels_x, pixels_y;
  float *dens, *denssum, rho_L;
  float *dp;

#ifdef TRACER_FIELD
  float *tracer, *tracersum, tracer_L;
#endif
#ifdef TRACER_MC
  float *tracerMC, *tracerMCsum, tracerMC_L;
#endif
#ifdef VORONOI_VELOCITY_FIELD_2D
  float *velx, *velxsum, velx_L;
  float *vely, *velysum, vely_L;
#endif

  FILE *fd = 0;
  double l_dx, l_dy;
  int i, j;
  double minposx, maxposx, minposy, maxposy;
  int minpixx, maxpixx, minpixy, maxpixy;

  if(gradients_flag == 1)
    sprintf(buf, "%s/density_field_%03d", All.OutputDir, num);
  else if(gradients_flag == 0)
    sprintf(buf, "%s/density_field_nograds_%03d", All.OutputDir, num);
  else
    terminate("gradients_flag != 1 && gradients_flag != 0");

  mpi_printf("we start to generate a special dump... gradients_flag=%d\n", gradients_flag);

  pixels_x = VORONOI_FIELD_DUMP_PIXELS_X;
  pixels_y = VORONOI_FIELD_DUMP_PIXELS_Y;

  if(ThisTask == 0)
    {
      if(!(fd = fopen(buf, "w")))
        {
          char buf[1000];
          sprintf(buf, "can't open file `%s' for writing snapshot.\n", buf);
          terminate(buf);
        }

      my_fwrite(&pixels_x, sizeof(int), 1, fd);
      my_fwrite(&pixels_y, sizeof(int), 1, fd);
    }

  dens = mymalloc("dens", pixels_x * pixels_y * sizeof(float));
  denssum = mymalloc("denssum", pixels_x * pixels_y * sizeof(float));

  for(i = 0, dp = dens; i < pixels_x; i++)
    for(j = 0; j < pixels_y; j++)
      *dp++ = 0;

#ifdef TRACER_FIELD
  tracer = mymalloc("tracer", pixels_x * pixels_y * sizeof(float));
  tracersum = mymalloc("tracersum", pixels_x * pixels_y * sizeof(float));

  for(i = 0, dp = tracer; i < pixels_x; i++)
    for(j = 0; j < pixels_y; j++)
      *dp++ = 0;
#endif

#ifdef TRACER_MC
  tracerMC = mymalloc("tracerMC", pixels_x * pixels_y * sizeof(float));
  tracerMCsum = mymalloc("tracerMCsum", pixels_x * pixels_y * sizeof(float));

  for(i = 0, dp = tracerMC; i < pixels_x; i++)
    for(j = 0; j < pixels_y; j++)
      *dp++ = 0;
#endif

#ifdef VORONOI_VELOCITY_FIELD_2D
  velx = mymalloc("velx", pixels_x * pixels_y * sizeof(float));
  velxsum = mymalloc("velxsum", pixels_x * pixels_y * sizeof(float));
  vely = mymalloc("vely", pixels_x * pixels_y * sizeof(float));
  velysum = mymalloc("velysum", pixels_x * pixels_y * sizeof(float));

  for(i = 0, dp = velx; i < pixels_x; i++)
    for(j = 0; j < pixels_y; j++)
      *dp++ = 0;
  for(i = 0, dp = vely; i < pixels_x; i++)
    for(j = 0; j < pixels_y; j++)
      *dp++ = 0;
#endif

  for(int li = 0; li < NumGas; li++)
    {
      minposx = P[li].Pos[0] - pow(SphP[li].Volume, 1.0 / 2) / 2;
      maxposx = P[li].Pos[0] + pow(SphP[li].Volume, 1.0 / 2) / 2;
      minposy = P[li].Pos[1] - pow(SphP[li].Volume, 1.0 / 2) / 2;
      maxposy = P[li].Pos[1] + pow(SphP[li].Volume, 1.0 / 2) / 2;

      minpixx = (int) ceil(minposx / All.BoxSize * (pixels_x - 1));
      maxpixx = (int) floor(maxposx / All.BoxSize * (pixels_x - 1));
      minpixy = (int) ceil(minposy / All.BoxSize * (pixels_y - 1));
      maxpixy = (int) floor(maxposy / All.BoxSize * (pixels_y - 1));

      for(i = minpixx; i <= maxpixx; i++)
        {
          for(j = minpixy; j <= maxpixy; j++)
            {
              l_dx = ((double) i / (pixels_x - 1) * All.BoxSize) - P[li].Pos[0];
              l_dy = ((double) j / (pixels_y - 1) * All.BoxSize) - P[li].Pos[1];

#ifdef PERIODIC
#if !defined(REFLECTIVE_X)
              if(l_dx < -boxHalf_X)
                l_dx += boxSize_X;
              if(l_dx > boxHalf_X)
                l_dx -= boxSize_X;
#endif
#if !defined(REFLECTIVE_Y)
              if(l_dy < -boxHalf_Y)
                l_dy += boxSize_Y;
              if(l_dy > boxHalf_Y)
                l_dy -= boxSize_Y;
#endif
#endif

              if(gradients_flag == 1)
                {
                  rho_L = SphP[li].Density + SphP[li].Grad.drho[0] * l_dx + SphP[li].Grad.drho[1] * l_dy;
#ifdef TRACER_FIELD
                  tracer_L = SphP[li].Tracer + SphP[li].Grad.dtracer[0] * l_dx + SphP[li].Grad.dtracer[1] * l_dy;
#endif
#ifdef TRACER_MC
                  tracerMC_L = (SphP[li].Density + SphP[li].Grad.drho[0] * l_dx + SphP[li].Grad.drho[1] * l_dy) * get_number_of_tracers(li) * All.ReferenceTracerMCMass / P[li].Mass;
#endif
#ifdef VORONOI_VELOCITY_FIELD_2D
                  velx_L = SphP[li].Momentum[0] / P[li].Mass + SphP[li].Grad.dvel[0][0] * l_dx + SphP[li].Grad.dvel[0][1] * l_dy;
                  vely_L = SphP[li].Momentum[1] / P[li].Mass + SphP[li].Grad.dvel[1][0] * l_dx + SphP[li].Grad.dvel[1][1] * l_dy;
#endif
                }
              else if(gradients_flag == 0)
                {
                  rho_L = SphP[li].Density;
#ifdef TRACER_FIELD
                  tracer_L = SphP[li].Tracer;
#endif
#ifdef TRACER_MC
                  tracerMC_L = SphP[li].Density * get_number_of_tracers(li) * All.ReferenceTracerMCMass / P[li].Mass;
#endif
#ifdef VORONOI_VELOCITY_FIELD_2D
                  velx_L = SphP[li].Momentum[0] / P[li].Mass;
                  vely_L = SphP[li].Momentum[1] / P[li].Mass;
#endif
                }
              else
                terminate("gradients_flag != 1 && gradients_flag != 0");

              dens[i * pixels_y + j] = rho_L;

#ifdef TRACER_FIELD
              tracer[i * pixels_y + j] = tracer_L;
#endif
#ifdef TRACER_MC
              tracerMC[i * pixels_y + j] = tracerMC_L;
#endif
#ifdef VORONOI_VELOCITY_FIELD_2D
              velx[i * pixels_y + j] = velx_L;
              vely[i * pixels_y + j] = vely_L;
#endif
            }
        }
    }

  MPI_Reduce(dens, denssum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#ifdef TRACER_FIELD
  MPI_Reduce(tracer, tracersum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
#ifdef TRACER_MC
  MPI_Reduce(tracerMC, tracerMCsum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
#ifdef VORONOI_VELOCITY_FIELD_2D
  MPI_Reduce(velx, velxsum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(vely, velysum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#endif


  if(ThisTask == 0)
    {
      my_fwrite(denssum, sizeof(float), pixels_x * pixels_y, fd);
#ifdef TRACER_FIELD
      my_fwrite(tracersum, sizeof(float), pixels_x * pixels_y, fd);
#endif
#ifdef TRACER_MC
      my_fwrite(tracerMCsum, sizeof(float), pixels_x * pixels_y, fd);
#endif
#ifdef VORONOI_VELOCITY_FIELD_2D
      my_fwrite(velxsum, sizeof(float), pixels_x * pixels_y, fd);
      my_fwrite(velysum, sizeof(float), pixels_x * pixels_y, fd);
#endif

      fclose(fd);
    }

#ifdef VORONOI_VELOCITY_FIELD_2D
  myfree(velysum);
  myfree(vely);
  myfree(velxsum);
  myfree(velx);
#endif
#ifdef TRACER_MC
  myfree(tracerMCsum);
  myfree(tracerMC);
#endif
#ifdef TRACER_FIELD
  myfree(tracersum);
  myfree(tracer);
#endif

  myfree(denssum);
  myfree(dens);

  CPU_Step[CPU_MAKEIMAGES] += measure_time();
}
#endif

#endif
