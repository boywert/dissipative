/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/rt/rt_optical_depth.c
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

#include "../allvars.h"
#include "../proto.h"
#include "../voronoi.h"

#ifdef RT_SHORT_CHARACTERISTICS

static struct columnexch
{
  double Density;
  double nHI;
  double ColumnDensity0;
  double ColumnDensity1;
  double RayThroughCellDistance;
  double RayEntryDistance;
  double RayExitDistance;
  double SourceDistance;
  char Flag_Computed;
}
 *ColumnExch;

static char *Flag_Computed, *Flag_Visited;
static int *Upwind_DP;
static double *RayThroughCellDistance;
static double *RayEntryDistance;
static double *RayExitDistance;



/* this function calculates the column 
 * densities to the given source 
 */
void rt_calc_column_density(double *source_pos)
{
  Upwind_DP = mymalloc("Upwind_DP", NumGas * sizeof(int));
  RayThroughCellDistance = mymalloc("RayThroughCellDistance", NumGas * sizeof(double));
  RayEntryDistance = mymalloc("RayEntryDistance", NumGas * sizeof(double));
  RayExitDistance = mymalloc("RayExitDistance", NumGas * sizeof(double));


  rt_find_upstream_ray_intersections(source_pos);

  mpi_printf("upstream cells found.\n");

  rt_find_short_characteristics_column_densities(source_pos);

  mpi_printf("column densities found.\n");


  myfree(RayExitDistance);
  myfree(RayEntryDistance);
  myfree(RayThroughCellDistance);
  myfree(Upwind_DP);
}





void rt_find_short_characteristics_column_densities(double *source_pos)
{
  int i, iter = 0, nleft, nleft_tot;

  Flag_Computed = mymalloc("Flag_Computed", NumGas * sizeof(char));
  Flag_Visited = mymalloc("Flag_Visited", NumGas * sizeof(char));
  ColumnExch = (struct columnexch *) mymalloc("ColumnExch", Mesh_nimport * sizeof(struct columnexch));

  for(i = 0; i < NumGas; i++)
    Flag_Computed[i] = 0;

  do
    {
      rt_exchange_column_densities();

      for(i = 0; i < NumGas; i++)
        Flag_Visited[i] = 0;

      for(i = 0; i < NumGas; i++)
        {
          if(Flag_Visited[i] == 0)
            rt_try_to_compute_column(i);
        }

      for(i = 0, nleft = 0; i < NumGas; i++)
        {
          if(Flag_Computed[i] == 0)
            nleft++;
        }

      MPI_Allreduce(&nleft, &nleft_tot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      mpi_printf("RT optical depth calculation: iter=%d,  left are %d cells\n", iter, nleft_tot);

      iter++;

      if(iter > 100)
        terminate("too many iterations");
    }
  while(nleft_tot > 0);


  /* now convert to hydrogen atom column density */
  for(i = 0; i < NumGas; i++)
    {
      SphP[i].ColumnDensity0 *= HYDROGEN_MASSFRAC / PROTONMASS * All.UnitMass_in_g / All.HubbleParam;
      SphP[i].ColumnDensity1 *= HYDROGEN_MASSFRAC / PROTONMASS * All.UnitMass_in_g / All.HubbleParam;
    }

  myfree(ColumnExch);
  myfree(Flag_Visited);
  myfree(Flag_Computed);
}



void rt_try_to_compute_column(int i)
{
  int dp_upwind, cell_upwind;

  Flag_Visited[i] = 1;

  dp_upwind = Upwind_DP[i];

  if(dp_upwind == -1)
    terminate("invalid upwind cell found");


  if(dp_upwind == -2)           /* this is the cell with the source */
    {
      SphP[i].ColumnDensity0 = 0.0;
      SphP[i].ColumnDensity1 = SphP[i].Density * SphP[i].nHI * RayThroughCellDistance[i];

      Flag_Computed[i] = 1;
    }
  else
    {
      if(DP[dp_upwind].task == ThisTask && DP[dp_upwind].index < NumGas)
        {
          cell_upwind = DP[dp_upwind].index;

          if(!Flag_Computed[cell_upwind])
            rt_try_to_compute_column(cell_upwind);

          if(Flag_Computed[cell_upwind])
            {
              SphP[i].ColumnDensity0 = SphP[cell_upwind].ColumnDensity1 + (RayEntryDistance[i] - RayExitDistance[cell_upwind]) * SphP[cell_upwind].Density * SphP[cell_upwind].nHI;

              SphP[i].ColumnDensity1 = SphP[i].ColumnDensity0 + SphP[i].Density * SphP[i].nHI * RayThroughCellDistance[i];

              Flag_Computed[i] = 1;
            }
        }
      else
        {
          if(DP[dp_upwind].task == ThisTask)
            terminate("this case should not occur");


          if(DP[dp_upwind].task != ThisTask)    /* the upstream cell is a ghost cell from another processor */
            {
              cell_upwind = DP[dp_upwind].index;

              if(ColumnExch[cell_upwind].Flag_Computed)
                {
                  SphP[i].ColumnDensity0 = ColumnExch[cell_upwind].ColumnDensity1 +
                    (RayEntryDistance[i] - ColumnExch[cell_upwind].RayExitDistance) * ColumnExch[cell_upwind].Density * ColumnExch[cell_upwind].nHI;

                  SphP[i].ColumnDensity1 = SphP[i].ColumnDensity0 + SphP[i].Density * SphP[i].nHI * RayThroughCellDistance[i];

                  Flag_Computed[i] = 1;
                }
            }
        }
    }
}





/* this function finds for every cell the upstream cell towards the source, and the
   point at which the ray from the cell to the source goes into the upstream cell
*/
void rt_find_upstream_ray_intersections(double *source_pos)
{
  int i, j, p1, p2, n;
  double *CurrentEntry_t, *CurrentExit_t;

  CurrentEntry_t = mymalloc("CurrentEntry_t", NumGas * sizeof(double));
  CurrentExit_t = mymalloc("CurrentExit_t", NumGas * sizeof(double));

  for(i = 0; i < NumGas; i++)
    {
      CurrentEntry_t[i] = MAX_REAL_NUMBER;
      CurrentExit_t[i] = -MAX_REAL_NUMBER;
      Upwind_DP[i] = -1;

      SphP[i].SourceDistance = sqrt(pow(SphP[i].Center[0] - source_pos[0], 2) + pow(SphP[i].Center[1] - source_pos[1], 2) + pow(SphP[i].Center[2] - source_pos[2], 2));
    }

  for(i = 0; i < Nvf; i++)      /* loop over all cell faces */
    {
      for(j = 0; j < 2; j++)
        {
          if(j == 0)
            {
              p1 = VF[i].p1;
              p2 = VF[i].p2;
            }
          else
            {
              p1 = VF[i].p2;
              p2 = VF[i].p1;
            }

          if(DP[p1].index < 0 || DP[p2].index < 0)
            continue;

          if(DP[p1].task == ThisTask && DP[p1].index < NumGas)  /* is a local cell */
            {
              n = DP[p1].index; /* index of SPH particle */

              double dx = DP[p2].x - DP[p1].x;
              double dy = DP[p2].y - DP[p1].y;
              double dz = DP[p2].z - DP[p1].z;

              double sx = source_pos[0] - SphP[n].Center[0];
              double sy = source_pos[1] - SphP[n].Center[1];
              double sz = source_pos[2] - SphP[n].Center[2];

              double sd = sx * dx + sy * dy + sz * dz;
              double mx = 0.5 * (DP[p2].x + DP[p1].x) - SphP[n].Center[0];
              double my = 0.5 * (DP[p2].y + DP[p1].y) - SphP[n].Center[1];
              double mz = 0.5 * (DP[p2].z + DP[p1].z) - SphP[n].Center[2];

              double md = mx * dx + my * dy + mz * dz;

              if(sd > 0)
                {
                  double t = md / sd;

                  if(t < 0)
                    terminate("t < 0");

                  if(t < CurrentEntry_t[n] || Upwind_DP[n] < 0) /* new ray entry point */
                    {
                      CurrentEntry_t[n] = t;
                      Upwind_DP[n] = p2;
                    }
                }
              else if(sd < 0)
                {
                  double t = md / sd;

                  if(t > 0)
                    terminate("t > 0");

                  if(t > CurrentExit_t[n])      /* new ray exit point */
                    {
                      CurrentExit_t[n] = t;
                    }
                }
            }
        }
    }

  for(i = 0; i < NumGas; i++)
    {
      RayEntryDistance[i] = (1 - CurrentEntry_t[i]) * SphP[i].SourceDistance;
      RayExitDistance[i] = (1 - CurrentExit_t[i]) * SphP[i].SourceDistance;

      if(RayEntryDistance[i] < 0 || SphP[i].SourceDistance == 0)
        {
          Upwind_DP[i] = -2;    /* this cell contains the source */

          RayEntryDistance[i] = 0;
          RayExitDistance[i] = pow(SphP[i].Volume / (4.0 / 3 * M_PI), 1.0 / 3); /* use fiducial cell radius */
        }

      RayThroughCellDistance[i] = RayExitDistance[i] - RayEntryDistance[i];

      SphP[i].v_shell = 4 * M_PI / 3.0 * (pow(RayExitDistance[i], 3) - pow(RayEntryDistance[i], 3));
    }

  myfree(CurrentExit_t);
  myfree(CurrentEntry_t);
}


void rt_exchange_column_densities(void)
{
  int listp;
  struct columnexch *tmpColumnExch;
  int j, p, task, off;
  int ngrp, sendTask, recvTask, place;

  tmpColumnExch = (struct columnexch *) mymalloc("tmpColumnExch", Mesh_nexport * sizeof(struct columnexch));

  /* prepare data for export */
  for(j = 0; j < NTask; j++)
    Mesh_Send_count[j] = 0;

  for(p = 0; p < NumGas; p++)
    {
      if(P[p].Type == 0)
        {
          listp = List_P[p].firstexport;
          while(listp >= 0)
            {
              if((task = ListExports[listp].origin) != ThisTask)
                {
                  place = ListExports[listp].index;
                  off = Mesh_Send_offset[task] + Mesh_Send_count[task]++;

                  tmpColumnExch[off].Density = SphP[place].Density;
                  tmpColumnExch[off].nHI = SphP[place].nHI;
                  tmpColumnExch[off].ColumnDensity0 = SphP[place].ColumnDensity0;
                  tmpColumnExch[off].ColumnDensity1 = SphP[place].ColumnDensity1;
                  tmpColumnExch[off].RayThroughCellDistance = RayThroughCellDistance[place];
                  tmpColumnExch[off].RayEntryDistance = RayEntryDistance[place];
                  tmpColumnExch[off].RayExitDistance = RayExitDistance[place];
                  tmpColumnExch[off].SourceDistance = SphP[place].SourceDistance;
                  tmpColumnExch[off].Flag_Computed = Flag_Computed[place];
                }
              listp = ListExports[listp].nextexport;
            }
        }
    }

  /* exchange data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Mesh_Send_count[recvTask] > 0 || Mesh_Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&tmpColumnExch[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask]
                           * sizeof(struct columnexch), MPI_BYTE, recvTask, TAG_DENS_A,
                           &ColumnExch[Mesh_Recv_offset[recvTask]], Mesh_Recv_count[recvTask] * sizeof(struct columnexch), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  myfree(tmpColumnExch);
}





#endif
