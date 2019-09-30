/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/modgrav/modgrav_pm_periodic.c
 * \date        MM/YYYY
 * \author
 * \brief       routines that add the effective mass of the AMR nodes to the periodic PM grid
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

#if defined(PMGRID) && defined(PERIODIC)

#define GRIDX (PMGRID * STRETCHX * DBX + DBX_EXTRA)
#define GRIDY (PMGRID * STRETCHY * DBY + DBY_EXTRA)
#define GRIDZ (PMGRID * STRETCHZ * DBZ + DBZ_EXTRA)

#define  GRIDz (GRIDZ/2 + 1)
#define  GRID2 (2*GRIDz)

#if (GRIDX > 1024) || (GRIDY > 1024) || (GRIDZ > 1024)
typedef long long large_array_offset;   /* use a larger data type in this case so that we can always address all cells of the 3D grid with a single index */
#else
typedef unsigned int large_array_offset;
#endif

void modgrav_add_amr_nodes_to_periodic_grid(fft_real *rhogrid, MyFloat to_slab_fac, int fftsize, fft_plan myplan)
{
  int no, i;
  int mintask, maxtask, minslab, maxslab, actual_task;
  int num_nodes_to_exchange = 0;
  int alloc_size;
  size_t bytes;

  int minx, maxx, miny, maxy, minz, maxz;
  int ix, iy, iz;
  int actual_x, actual_y, actual_z;

  double dx, dy, dz;
  double pm_cell_size = All.BoxSize / PMGRID;
  double node_vol, mass_for_node;

  int max_export_nodes = (int) ceil(1.1 * Tree_MaxNodes);

  struct pm_amr_data *exchange_nodes;
  if(!(exchange_nodes = (struct pm_amr_data *) mymalloc("exchange_nodes", bytes = max_export_nodes * sizeof(struct pm_amr_data))))
    {
      printf("failed to allocate memory for `exchange_nodes' (%g MB), presently allocated=%g MB\n", bytes / (1024.0 * 1024.0), AllocatedBytes / (1024.0 * 1024.0));
      terminate("Something wrong");
    }

  alloc_size = max_export_nodes;

  /* first loop, prepare nodes for export */
  for(no = All.MaxPart; no < All.MaxPart + Tree_NumNodes; no++)
    if(Nodes[no].mg_bitflag.amr_node == 1)
      {
        if(no < Tree_FirstNonTopLevelNode)    /* top nodes are also exported, this simplifies the code and should not have a big impact on performance */
          if((no - All.MaxPart) % NTask != ThisTask)    /* makes sure that topnodes are only exported by one task */
            continue;

        minslab = (int) floor(to_slab_fac * (Nodes[no].center[0] - 0.5 * Nodes[no].len) + 0.5);
        if(minslab == PMGRID)
          minslab = 0;
        maxslab = (int) floor(to_slab_fac * (Nodes[no].center[0] + 0.5 * Nodes[no].len) + 0.5);
        if(maxslab == PMGRID)
          maxslab = 0;

        mintask = myplan.slab_to_task[minslab];
        maxtask = myplan.slab_to_task[maxslab];

        if(mintask > 0 && maxtask == 0)
          maxtask = NTask;

        for(i = mintask; i <= maxtask; ++i)
          {
            if(i == NTask)
              actual_task = 0;
            else
              actual_task = i;

            if(myplan.slabs_x_per_task[actual_task] > 0)
              {
                if(num_nodes_to_exchange > max_export_nodes)
                  {
                    printf("num_nodes_to_exchange = %i, max_export_nodes = %i \n", num_nodes_to_exchange, max_export_nodes);
                    terminate("too little memory allocated for exchange_nodes");
                  }

                exchange_nodes[num_nodes_to_exchange].task = actual_task;
                exchange_nodes[num_nodes_to_exchange].center[0] = Nodes[no].center[0];
                exchange_nodes[num_nodes_to_exchange].center[1] = Nodes[no].center[1];
                exchange_nodes[num_nodes_to_exchange].center[2] = Nodes[no].center[2];
                exchange_nodes[num_nodes_to_exchange].len = Nodes[no].len;
                exchange_nodes[num_nodes_to_exchange].eff_mass = Nodes[no].mg.eff.eff_mass;
                ++num_nodes_to_exchange;
              }
          }
      }

  /* exchange of node data */
  mpi_distribute_items_to_tasks((void *) exchange_nodes, 0, &num_nodes_to_exchange, &alloc_size, sizeof(struct pm_amr_data), TAG_MG_PM);

  /* second loop, putting imported nodes on PM grid */
  for(i = 0; i < num_nodes_to_exchange; ++i)
    {
      mass_for_node = 0.0;

      if(exchange_nodes[i].task != ThisTask)
        terminate("something is wrong here!");

      node_vol = exchange_nodes[i].len * exchange_nodes[i].len * exchange_nodes[i].len;

      minx = (int) floor(to_slab_fac*(exchange_nodes[i].center[0] - 0.5*exchange_nodes[i].len) + 0.5);
      maxx = (int) floor(to_slab_fac*(exchange_nodes[i].center[0] + 0.5*exchange_nodes[i].len) + 0.5);
      miny = (int) floor(to_slab_fac*(exchange_nodes[i].center[1] - 0.5*exchange_nodes[i].len) + 0.5);
      maxy = (int) floor(to_slab_fac*(exchange_nodes[i].center[1] + 0.5*exchange_nodes[i].len) + 0.5);
      minz = (int) floor(to_slab_fac*(exchange_nodes[i].center[2] - 0.5*exchange_nodes[i].len) + 0.5);
      maxz = (int) floor(to_slab_fac*(exchange_nodes[i].center[2] + 0.5*exchange_nodes[i].len) + 0.5);

      for(ix = minx; ix <= maxx; ++ix)
        {
          if(ix == PMGRID)
            actual_x = 0;
          else
            actual_x = ix;

          if(actual_x >= myplan.first_slab_x_of_task[ThisTask] && actual_x < myplan.first_slab_x_of_task[ThisTask] + myplan.slabs_x_per_task[ThisTask])
            {
              dx = fabs(exchange_nodes[i].center[0] - ix * pm_cell_size);
              if(dx > 0.5 * All.BoxSize)
                terminate("seems strange");

              dx = 0.5 * exchange_nodes[i].len + 0.5 * pm_cell_size - dx;  /* calculate overlap in x */
              if(exchange_nodes[i].len < dx)
                dx = exchange_nodes[i].len;
              if(pm_cell_size < dx)
                dx = pm_cell_size;

              for(iy = miny; iy <= maxy; ++iy)
                {
                  if(iy == PMGRID)
                    actual_y = 0;
                  else
                    actual_y = iy;

                  dy = fabs(exchange_nodes[i].center[1] - iy * pm_cell_size);
                  if(dy > 0.5 * All.BoxSize)
                    terminate("seems strange");

                  dy = 0.5 * exchange_nodes[i].len + 0.5 * pm_cell_size - dy;      /* calculate overlap in y */
                  if(exchange_nodes[i].len < dy)
                    dy = exchange_nodes[i].len;
                  if(pm_cell_size < dy)
                    dy = pm_cell_size;

                  for(iz = minz; iz <= maxz; ++iz)
                    {
                      if(iz == PMGRID)
                        actual_z = 0;
                      else
                        actual_z = iz;

                      dz = fabs(exchange_nodes[i].center[2] - iz * pm_cell_size);
                      if(dz > 0.5 * All.BoxSize)
                        terminate("seems strange");

                      dz = 0.5 * exchange_nodes[i].len + 0.5 * pm_cell_size - dz;  /* calculate overlap in z */
                      if(exchange_nodes[i].len < dz)
                        dz = exchange_nodes[i].len;
                      if(pm_cell_size < dz)
                        dz = pm_cell_size;

                      large_array_offset offset = ((large_array_offset) GRID2) * (GRIDZ * (actual_x - myplan.first_slab_x_of_task[ThisTask]) + actual_y) + actual_z;
                      if(offset < 0 || offset > fftsize)
                        terminate("incorrect offset");

                      rhogrid[offset] += exchange_nodes[i].eff_mass * dx * dy * dz / node_vol;
                      mass_for_node += exchange_nodes[i].eff_mass * dx * dy * dz / node_vol;

                    }
                }
            }
        }

      if(fabs(mass_for_node) > fabs(1.000001 * exchange_nodes[i].eff_mass))
        {
          printf("mass_for_node %f, exchange_nodes[i].eff_mass %f\n", mass_for_node, exchange_nodes[i].eff_mass);
          printf("minx %i, maxx %i, miny %i, maxy %i, minz %i, maxz %i\n", minx, maxx, miny, maxy, minz, maxz);
          printf("first slab %i, last slab %i \n", myplan.first_slab_x_of_task[ThisTask], myplan.first_slab_x_of_task[ThisTask] + myplan.slabs_x_per_task[ThisTask] - 1);
          terminate("something wrong in mass assignment");
        }
    }

  myfree(exchange_nodes);
}
#endif
