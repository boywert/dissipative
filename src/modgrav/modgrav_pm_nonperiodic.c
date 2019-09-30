/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/modgrav/modgrav_pm_nonperiodic.c
 * \date        MM/YYYY
 * \author
 * \brief       routines for adding the effective mass to the nonperiodic PM-grid
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


/*!
 *  \brief routines for adding the effective mass to the PM-grids
 */

#include "../allvars.h"
#include "../proto.h"

#if defined(PMGRID) && (!defined(PERIODIC) || defined(PLACEHIGHRESREGION) || defined(PLACEHIGHRESREGION) || defined(GRAVITY_NOT_PERIODIC))

#if defined(LONG_X) || defined(LONG_Y) || defined (LONG_Z)
#error "LONG_X/Y/Z not supported for the non-periodic FFT gravity code"
#endif

#if defined(GRAVITY_TALLBOX)
#error "GRAVITY_TALLBOX not supported for the non-periodic FFT gravity code"
#endif

#ifndef GRIDBOOST
#define GRIDBOOST 2
#endif

#define  GRID  (GRIDBOOST*PMGRID)
#define  GRIDz (GRID/2 + 1)
#define  GRID2 (2*GRIDz)

#if (GRID > 1024)
typedef long long large_array_offset;   /* use a larger data type in this case so that we can always address all cells of the 3D grid with a single index */
#else
typedef unsigned int large_array_offset;
#endif


void modgrav_add_amr_nodes_to_nonperiodic_grid(int grnr, fft_real *rhogrid, double to_slab_fac, fft_plan myplan,  int fftsize)
{
  int no, i, k;
  int mintask, maxtask, minslab, maxslab, actual_task;
  int num_nodes_to_exchange = 0;
  int alloc_size;
  size_t bytes;

  int minx[3], maxx[3], ix[3];

  double overlap[3];
  double pm_cell_size = All.TotalMeshSize[grnr] / GRID;
  double node_vol, mass_for_node;

  large_array_offset offset;

  int max_export_nodes = (int) ceil(2.5 * Tree_MaxNodes);

  struct pm_amr_data *exchange_nodes;
  if(!(exchange_nodes = (struct pm_amr_data *) mymalloc("exchange_nodes", bytes = max_export_nodes * sizeof(struct pm_amr_data))))
    {
      printf("failed to allocate memory for `exchange_nodes' (%g MB), presently allocated=%g MB\n", bytes / (1024.0 * 1024.0), AllocatedBytes / (1024.0 * 1024.0));
      terminate("Something wrong");
    }

  to_slab_fac = GRID / All.TotalMeshSize[grnr];

  alloc_size = max_export_nodes;

  /* first loop, prepare nodes for export */
  for(no = All.MaxPart; no < All.MaxPart + Tree_NumNodes; no++)
    if(Nodes[no].mg_bitflag.amr_node == 1)
      {
        if(no < Tree_FirstNonTopLevelNode)    /* top nodes are also exported, this simplifies the code and should not have a big impact on performance */
          if((no - All.MaxPart) % NTask != ThisTask)    /* makes sure that topnodes are only exported by one task */
            continue;

        /*pos = Nodes[no].s_eff; */

        if(Nodes[no].center[0] + 0.5 * Nodes[no].len < All.Corner[grnr][0] || Nodes[no].center[0] - 0.5 * Nodes[no].len >= All.UpperCorner[grnr][0])
          continue;
        if(Nodes[no].center[1] + 0.5 * Nodes[no].len < All.Corner[grnr][1] || Nodes[no].center[1] - 0.5 * Nodes[no].len >= All.UpperCorner[grnr][1])
          continue;
        if(Nodes[no].center[2] + 0.5 * Nodes[no].len < All.Corner[grnr][2] || Nodes[no].center[2] - 0.5 * Nodes[no].len >= All.UpperCorner[grnr][2])
          continue;

        minslab = (int) floor(to_slab_fac * (Nodes[no].center[0] - All.Corner[grnr][0] - 0.5 * Nodes[no].len) + 0.5);
        if(minslab < 0)
          minslab = 0;
        else if(minslab >= GRID)
          minslab = GRID - 1;
        maxslab = (int) floor(to_slab_fac * (Nodes[no].center[0] - All.Corner[grnr][0] + 0.5 * Nodes[no].len) + 0.5);
        if(maxslab >= GRID)
          maxslab = GRID - 1;
        else if(maxslab < 0)
          maxslab = 0;

        if(minslab < 0 || maxslab >= GRID || minslab > maxslab)
          {
            printf("task: %i, no = %i, minslab = %i, maxslab = %i \n", ThisTask, no, minslab, maxslab);
            fflush(stdout);
            terminate("Something wrong with node export.");
          }

        mintask = myplan.slab_to_task[minslab];
        maxtask = myplan.slab_to_task[maxslab];

        if(minslab < 0 || maxslab >= GRID)
          {
            printf("task: %i, no = %i, mintask =%i, maxtask = %i, slab to task done.\n", ThisTask, no, mintask, maxtask);
            fflush(stdout);
          }

        /*
         *if (mintask>0 && maxtask==0) maxtask = NTask;
         */
        for(i = mintask; i <= maxtask; ++i)
          {
            actual_task = i;

            if(myplan.slabs_x_per_task[actual_task] > 0)
              {
                if(num_nodes_to_exchange > max_export_nodes)
                  {
                    printf("num_nodes_to_exchange = %i, max_export_nodes = %i \n", num_nodes_to_exchange, max_export_nodes);
                    terminate("too little memory allocated for exchange_nodes");
                  }

                exchange_nodes[num_nodes_to_exchange].task = actual_task;
                exchange_nodes[num_nodes_to_exchange].center[0] = Nodes[no].center[0] - All.Corner[grnr][0];
                exchange_nodes[num_nodes_to_exchange].center[1] = Nodes[no].center[1] - All.Corner[grnr][1];
                exchange_nodes[num_nodes_to_exchange].center[2] = Nodes[no].center[2] - All.Corner[grnr][2];
                exchange_nodes[num_nodes_to_exchange].len = Nodes[no].len;
                exchange_nodes[num_nodes_to_exchange].eff_mass = Nodes[no].mg.eff.eff_mass;

                ++num_nodes_to_exchange;
              }
          }
      }

  mpi_distribute_items_to_tasks((void *) exchange_nodes, 0, &num_nodes_to_exchange, &alloc_size, sizeof(struct pm_amr_data), TAG_MG_PM);

  /* second loop, putting imported nodes on PM grid */
  for(i = 0; i < num_nodes_to_exchange; ++i)
    {
      mass_for_node = 0.0;

      if(exchange_nodes[i].task != ThisTask)
        terminate("something is wrong here!");

      node_vol = exchange_nodes[i].len * exchange_nodes[i].len * exchange_nodes[i].len;

      for(k=0; k<3; k++)
        {
          minx[k] = (int) floor(to_slab_fac*(exchange_nodes[i].center[k] - 0.5*exchange_nodes[i].len) + 0.5);
          maxx[k] = (int) floor(to_slab_fac*(exchange_nodes[i].center[k] + 0.5*exchange_nodes[i].len) + 0.5);

          if(minx[k] < 0)
            minx[k] = 0;
          if(maxx[k] >= GRID)
            maxx[k] = GRID - 1;

          if(maxx[k] < minx[k])
            terminate("Something wrong");
        }

      for(ix[0] = minx[0]; ix[0] <= maxx[0]; ++ix[0])
        if(ix[0] >= myplan.first_slab_x_of_task[ThisTask] && ix[0] < myplan.first_slab_x_of_task[ThisTask] + myplan.slabs_x_per_task[ThisTask])
          for(ix[1] = minx[1]; ix[1] <= maxx[1]; ++ix[1])
            for(ix[2] = minx[2]; ix[2] <= maxx[2]; ++ix[2])
              {
                for(k=0; k<3; k++)
                  {
                    overlap[k] = fabs(exchange_nodes[i].center[k] - ix[k] * pm_cell_size);

                    if(overlap[k] > 0.5 * All.BoxSize)
                      terminate("seems strange");

                    overlap[k] = 0.5 * exchange_nodes[i].len + 0.5 * pm_cell_size - overlap[k];  /* calculate overlap in x[k] */

                    if(exchange_nodes[i].len < overlap[k])
                      overlap[k] = exchange_nodes[i].len;
                    if(pm_cell_size < overlap[k])
                      overlap[k] = pm_cell_size;
                  }

                offset = ((large_array_offset) GRID2) * (GRID * (ix[0] - myplan.first_slab_x_of_task[ThisTask]) + ix[1]) + ix[2];
                if(offset < 0 || offset > fftsize)
                  {
                    printf("Task %i, offset %ui, fftsize %i, first_slab_of_task[ThisTask] %i, actual_x = [%i, %i, %i]\n",
                           ThisTask, offset, fftsize, myplan.first_slab_x_of_task[ThisTask], ix[0], ix[1], ix[2]);
                    terminate("incorrect offset");
                  }

                rhogrid[offset] += exchange_nodes[i].eff_mass * overlap[0] * overlap[1] * overlap[2] / node_vol;
                mass_for_node += exchange_nodes[i].eff_mass * overlap[0] * overlap[1] * overlap[2] / node_vol;
              }

      if(fabs(mass_for_node) > fabs(1.000001 * exchange_nodes[i].eff_mass))
        {
          printf("mass_for_node %f, exchange_nodes[i].eff_mass %f\n", mass_for_node, exchange_nodes[i].eff_mass);
          printf("minx %i, maxx %i, miny %i, maxy %i, minz %i, maxz %i\n", minx[0], maxx[0], minx[1], maxx[1], minx[2], maxx[2]);
          printf("first slab %i, last slab %i \n", myplan.first_slab_x_of_task[ThisTask], myplan.first_slab_x_of_task[ThisTask] + myplan.slabs_x_per_task[ThisTask] - 1);
          terminate("something wrong in mass assignment");
        }
    }
  myfree(exchange_nodes);
}
#endif
