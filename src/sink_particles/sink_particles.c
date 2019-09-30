
#include "../allvars.h"
#include "../proto.h"

void sink_particles(void)
{
  double t0, t1;

  CPU_Step[CPU_MISC] += measure_time();

  if(All.Time == All.TimeBegin || SinkAccretionRadius == 0)
    return;

  t0 = second();


  /* get the information about all the sinks that are present 
     in the simulation
  */
  int num_sinks;
  num_sinks = get_all_sink_particle_info(1);
  mpi_printf("SINK_PARTICLES: number of sinks on all tasks %d \n", num_sinks);

#ifndef SINK_PARTICLE_TEST
  /* Allow the current sinks to accrete if present
  */
  mpi_printf("SINK_PARTICLES: Entering accretion function... NSinksAllTasks %d \n", NSinksAllTasks); 
  double mass_accreted = 0;
  if(NSinksAllTasks > 0)
    mass_accreted = accrete_onto_sink_particles();
  mpi_printf("SINK_PARTICLES: mass accreted by all sinks this timestep %g \n", mass_accreted);
#endif

  /* Allow new sinks to form 
  */
  int success;
  success = create_sink_particles();
  if(success <= 0  &&  ThisTask == 0) 
    printf("SINK_PARTICLES: success flag %d \n", success);
#ifndef SINK_PARTICLE_TEST
  if(success > 0)
    {
      /* New sink present, so we need to updated the SinkP lists on all tasks
         before letting the new sink accrete mass
      */
      num_sinks = get_all_sink_particle_info(1);
      if(ThisTask == 0)
        printf("SINK_PARTICLES: number of sinks on all tasks after creation %d \n", num_sinks);
      mass_accreted = accrete_onto_sink_particles();
      mpi_printf("SINK_PARTICLES: mass accreted after sink creation %g \n", mass_accreted);
    }
#endif

#ifdef DUMP_SINK_PARTICLE_INFO 
  /* Dump the sink particle data to a file
  */
  if(NSinksAllTasks > 0)
    dump_sink_particle_info(success);
#endif
  

  MPI_Barrier(MPI_COMM_WORLD);

  t1 = second();

  CPU_Step[CPU_SINKS] += measure_time();
}

