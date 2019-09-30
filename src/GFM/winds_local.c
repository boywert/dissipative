/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/GFM/winds_local.c
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
#include "../voronoi.h"



#ifdef GFM_WINDS_LOCAL

#ifdef GFM_WINDS
#error "GFM_WINDS_LOCAL not allowed together with GFM_WINDS"
#endif

void create_winds_local(void)
{
  int i;
  int local_winds = 0, global_winds;
  int local_excess = 0, global_excess;
  int local_received = 0, global_received;
  double windenergy_local_sum_should = 0, windenergy_global_sum_should;
  double windenergy_local_sum_is = 0, windenergy_global_sum_is;
  double vel, prob;
  double halo_sigma;

  WindEnergy_Should = 0.0;
  WindEnergy_Is = 0.0;

  for(i = 0; i < NumGas; i++)
    {
      halo_sigma = GFM_FAC_SIGMA * SphP[i].w.DMVelDisp;
      vel = (All.VariableWindVelFactor * halo_sigma) > All.MinWindVel ? (All.VariableWindVelFactor * halo_sigma) : All.MinWindVel;
      if(vel > 0 && SphP[i].WindEnergyReceived > 0)
        {
          if(P[i].Mass == 0 && P[i].ID == 0)
            terminate("GFM_WINDS_LOCAL: cell gone i=%d WindEnergyReceived=%g\n", i, SphP[i].WindEnergyReceived);

          windenergy_local_sum_should += SphP[i].WindEnergyReceived;
          prob = SphP[i].WindEnergyReceived / (0.5 * P[i].Mass * vel * vel);
          if(prob > 1)
            local_excess++;
          if(prob > 0)
            local_received++;
          if(get_random_number() < prob)
            {
              windenergy_local_sum_is += 0.5 * P[i].Mass * vel * vel;
              gfm_add_wind(i, vel, 0);
              local_winds++;
            }
        }
    }


  MPI_Reduce(&local_winds, &global_winds, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_excess, &global_excess, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_received, &global_received, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  MPI_Allreduce(&windenergy_local_sum_should, &windenergy_global_sum_should, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&windenergy_local_sum_is, &windenergy_global_sum_is, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  WindEnergy_Should += windenergy_global_sum_should;
  WindEnergy_Is += windenergy_global_sum_is;

  mpi_printf("GFM_WINDS_LOCAL: generated %d/%d wind particles of which %d have excess energy\n", global_winds, global_received, global_excess);
}
#endif /* GFM_WINDS_LOCAL */
