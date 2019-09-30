/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/cosmic_rays_shock_acceleration.c
 * \date        09/2015
 * \author      C. Pfrommer
 * \brief       Simple recipes for CR acceleration at shocks
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 */

#include "../allvars.h"
#include "../proto.h"
#include "string.h"

#ifdef COSMIC_RAYS_SHOCK_ACCELERATION
/** \brief Inject cosmic rays at shocks in the simulations
 *
 *  We assume instantaneous deposition of cosmic rays whenever the
 *  on-the-fly shock finder detects a shock. The amount of energy
 *  deposited in cosmic rays may depend on Mach number and magnetic
 *  obliquity (future)
 *    
 */

static double getUthermEstimate( int cell )
{
  double Utherm;
  if(TimeBinSynchronized[P[cell].TimeBinHydro])
    {
      Utherm = SphP[cell].Utherm;
    }
  else
    {
      double Vel[3];
      for(int k=0; k<3; k++) Vel[k] = SphP[cell].Momentum[k] / P[cell].Mass;
      Utherm = (SphP[cell].Energy / P[cell].Mass - 0.5 * (Vel[0]*Vel[0] + Vel[1]*Vel[1] + Vel[2]*Vel[2])) / (All.cf_atime * All.cf_atime);
#ifdef MHD
      double B[3];
      for(int k=0; k<3; k++) B[k] = SphP[cell].BConserved[k] / SphP[cell].Volume;
      Utherm -= 0.5 * (B[0]*B[0] + B[1]*B[1] + B[2]*B[2]) / (P[cell].Mass/SphP[cell].Volume) / All.cf_atime;
#endif
    }
  return Utherm;
}

static double lost_injected_energy_sum = 0.;

void accelerate_cosmic_rays_at_shock(shock_ray * first_ray, shock_ray * last_ray)
{
  //send_count[i]: number of elements which have to be send to task i
  unsigned int *send_count = (unsigned int *) mymalloc("send_count", NTask * sizeof(unsigned int));
  memset(send_count, 0, NTask * sizeof(unsigned int));

  //receive_count[i]: number of elements which have to be received from task i
  unsigned int *receive_count = (unsigned int *) mymalloc("receive_count", NTask * sizeof(unsigned int));
  memset(receive_count, 0, NTask * sizeof(unsigned int));

  shock_ray *current_ray = first_ray;

  double eth_tot_pre_shock, total_energy, dissipated_energy, dt_cell;
  int i;
  int cell = 0;
  int task = 0;

  int export_counter = 0;

  double lost_injected_energy = 0.;
  double lost_injected_energy_all = 0.;
  double acceleration_efficiency = 0;

#ifdef COSMIC_RAYS_MAGNETIC_OBLIQUITY
  double critical_obliquity = M_PI / 4.;
  double delta = M_PI / 18.;
#endif

  //count number of elements for export
  while(current_ray != last_ray)
    {
      if(current_ray->current_cell == -1 && SDATA(current_ray->original_cell).ShockSurface)
        {
          cell = current_ray->original_cell;

          if(SphP[cell].Machnumber > All.CriticalMachnumber)
            {
              for(i = 0; i < current_ray->cell_count; i++)
                {
                  task = current_ray->post_shock_cells_tasks[i];

                  if(task != ThisTask)
                    {
                      export_counter++;
                    }
                }
            }
        }
      current_ray++;
    }

  //memory for export ray indices
  int *cells_export = (int *) mymalloc("cells_export", export_counter * sizeof(int));
  int *to_task = (int *) mymalloc("to_task", export_counter * sizeof(int));
  double *diss_per_tot = (double *) mymalloc("diss_per_tot", export_counter * sizeof(double));
  double *tot_pre_shock = (double *) mymalloc("tot_pre_shock", export_counter * sizeof(double));
#ifdef COSMIC_RAYS_MAGNETIC_OBLIQUITY
  double *acc_efficiencies = (double *) mymalloc("acc_efficiencies", export_counter * sizeof(double));
#endif

  export_counter = 0;
  current_ray = first_ray;

  //second loop, caluclate stuff!
  while(current_ray != last_ray)
    {
      if(current_ray->current_cell == -1 && SDATA(current_ray->original_cell).ShockSurface)
        {
          cell = current_ray->original_cell;

          if(SphP[cell].Machnumber > All.CriticalMachnumber)
            {
#ifdef COSMIC_RAYS_MAGNETIC_OBLIQUITY
              acceleration_efficiency = All.AccelerationEfficiency * 0.5 * (tanh(-(current_ray->magnetic_obliquity_pre_shock - critical_obliquity) / delta) + 1);
#else
              acceleration_efficiency = All.AccelerationEfficiency;
#endif

              dt_cell = (P[cell].TimeBinHydro ? (((integertime) 1) << P[cell].TimeBinHydro) : 0) * All.Timebase_interval / All.cf_hubble_a;
              dissipated_energy = SphP[cell].EnergyDissipation * dt_cell;
              eth_tot_pre_shock = current_ray->eth_tot_pre_shock;

              assert(current_ray->post_shock_cells[current_ray->cell_count] == -1 || current_ray->cell_count == 3);
              assert(current_ray->eth_tot_pre_shock > 0);
              assert(current_ray->cell_count >= 1);
              assert(current_ray->post_shock_cells[current_ray->cell_count - 1] != -1);

              total_energy = current_ray->post_shock_eth_tot_sum - current_ray->eth_tot_pre_shock * (current_ray->cell_count + 1);

              double dEnergy_tmp = (injection_weight(cell) - eth_tot_pre_shock) * acceleration_efficiency * dissipated_energy / total_energy / All.cf_atime;    /* in physical units */
              
              double dEnergy = dmax(0., dmin(dEnergy_tmp, (getUthermEstimate(cell) - All.MinEgySpec) * P[cell].Mass));
              lost_injected_energy += dEnergy_tmp - dEnergy;

              /* convert energies to comoving units and inject CR Energy from the thermal pool */
              SphP[cell].Utherm -= dEnergy / P[cell].Mass;
              SphP[cell].Energy -= dEnergy * All.cf_atime * All.cf_atime;
              SphP[cell].CR_Energy += dEnergy * All.cf_atime;

              assert(SphP[cell].Energy > 0);

              for(i = 0; i < current_ray->cell_count; i++)
                {
                  task = current_ray->post_shock_cells_tasks[i];

                  if(task != ThisTask)
                    {
                      send_count[task]++;
                      cells_export[export_counter] = current_ray->post_shock_cells[i];
                      to_task[export_counter] = task;
                      diss_per_tot[export_counter] = dissipated_energy / total_energy;
                      tot_pre_shock[export_counter] = eth_tot_pre_shock;
#ifdef COSMIC_RAYS_MAGNETIC_OBLIQUITY
                      acc_efficiencies[export_counter] = acceleration_efficiency;
#endif
                      export_counter++;
                    }
                  else
                    {
                      cell = current_ray->post_shock_cells[i];

                      dEnergy_tmp = (injection_weight(cell) - eth_tot_pre_shock) * acceleration_efficiency * dissipated_energy / total_energy / All.cf_atime;
                      
                      dEnergy = dmax(0., dmin(dEnergy_tmp, (getUthermEstimate(cell) - All.MinEgySpec) * P[cell].Mass));
                      lost_injected_energy += dEnergy_tmp - dEnergy;

                      SphP[cell].Utherm -= dEnergy / P[cell].Mass;
                      SphP[cell].Energy -= dEnergy * All.cf_atime * All.cf_atime;
                      SphP[cell].CR_Energy += dEnergy * All.cf_atime;

                      if(SphP[cell].Energy <= 0)
                        {
                          print_particle_info( cell );
                          printf( "dEnergy=%g, dUtherm=%g, atime=%g\n", -dEnergy / P[cell].Mass, -dEnergy, All.cf_atime );
                        }

                      assert(SphP[cell].Energy > 0);
                    }
                }
            }
        }
      current_ray++;
    }

  //parallelization
  struct data_exch
  {
    int cell;
    double tot_pre_shock;
    double diss_per_tot_energy;
#ifdef COSMIC_RAYS_MAGNETIC_OBLIQUITY
    double acc_efficiency;
#endif
  } *data_import, *data_export;


  MPI_Alltoall(send_count, 1, MPI_INT, receive_count, 1, MPI_INT, MPI_COMM_WORLD);

  unsigned int nimport = 0;
  unsigned int nexport = 0;
  unsigned int send_offset[NTask];
  unsigned int receive_offset[NTask];
  unsigned int j;

  //calculate nimport/nexport, send_offset/receive_offset
  for(j = 0, nimport = 0, nexport = 0, receive_offset[0] = 0, send_offset[0] = 0; j < NTask; j++)
    {
      nimport += receive_count[j];
      nexport += send_count[j];

      if(j > 0)
        {
          send_offset[j] = send_offset[j - 1] + send_count[j - 1];
          receive_offset[j] = receive_offset[j - 1] + receive_count[j - 1];
        }
    }

  //prepare data for export
  for(j = 0; j < NTask; j++)
    {
      send_count[j] = 0;
    }

  unsigned int off;

  data_export = (struct data_exch *) mymalloc("data_export", nexport * sizeof(struct data_exch));
  data_import = (struct data_exch *) mymalloc("data_import", nimport * sizeof(struct data_exch));

  for(j = 0; j < export_counter; j++)
    {
      off = send_offset[to_task[j]] + send_count[to_task[j]]++;
      data_export[off].cell = cells_export[j];
      data_export[off].diss_per_tot_energy = diss_per_tot[j];
      data_export[off].tot_pre_shock = tot_pre_shock[j];
#ifdef COSMIC_RAYS_MAGNETIC_OBLIQUITY
      data_export[off].acc_efficiency = acc_efficiencies[j];
#endif
    }

  int ngrp, recvTask;

  // exchange data
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(send_count[recvTask] > 0 || receive_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&data_export[send_offset[recvTask]], send_count[recvTask] * sizeof(struct data_exch), MPI_BYTE, recvTask, 0,
                           &data_import[receive_offset[recvTask]], receive_count[recvTask] * sizeof(struct data_exch), MPI_BYTE, recvTask, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  //inject into cells on other cores
  for(j = 0; j < nimport; j++)
    {
#ifdef COSMIC_RAYS_MAGNETIC_OBLIQUITY
      acceleration_efficiency = data_import[j].acc_efficiency;
#else
      acceleration_efficiency = All.AccelerationEfficiency;
#endif

      double dEnergy_tmp = (injection_weight(data_import[j].cell) - data_import[j].tot_pre_shock) * acceleration_efficiency * data_import[j].diss_per_tot_energy / All.cf_atime;
      double dEnergy = dmax(0., dmin(dEnergy_tmp, (SphP[data_import[j].cell].Utherm - All.MinEgySpec) * P[data_import[j].cell].Mass));
      lost_injected_energy += dEnergy_tmp - dEnergy;

      SphP[data_import[j].cell].Utherm -= dEnergy / P[data_import[j].cell].Mass;
      SphP[data_import[j].cell].Energy -= dEnergy * All.cf_atime * All.cf_atime;
      SphP[data_import[j].cell].CR_Energy += dEnergy * All.cf_atime;

      assert(SphP[cell].Energy > 0);
    }

  myfree(data_import);
  myfree(data_export);
#ifdef COSMIC_RAYS_MAGNETIC_OBLIQUITY
  myfree(acc_efficiencies);
#endif
  myfree(tot_pre_shock);
  myfree(diss_per_tot);
  myfree(to_task);
  myfree(cells_export);
  myfree(receive_count);
  myfree(send_count);

  update_primitive_variables();

  MPI_Reduce(&lost_injected_energy, &lost_injected_energy_all, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  lost_injected_energy_sum += lost_injected_energy_all;
  mpi_printf("COSMIC_RAYS: lost injected energy E = %g cumulative E = %g.\n", lost_injected_energy_all, lost_injected_energy_sum);

}
#endif
