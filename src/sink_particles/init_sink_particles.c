#include "../allvars.h"
#include "../proto.h"

/* PCC - 02.01.2014
   	Takes care of the initialization of the ITA-style sink
        particles. Called early in AREPO, with different modes
        depending on where it's called (before/after ICs/snapshot)
        is read
*/
void init_sink_particles(int mode)
  {
 
    int LastTaskNewSink;
    int sink_histories_success;
    char msg[100];

    mpi_printf("SINK_PARTICLES: Initializing the global variables. In mode %d \n", mode);    

    if (mode == 0)
     {
       /* This branch is called from begrun.c
       */
    
       NSinksAllTasksOld = 0;
       NSinksAllTasks = 0;
       NSinkBufferSize = 1000;
       LastTaskNewSink = 0;
       SinksFormedSinceLastDomain = 0;
    

       /* Broadcast the variables set in the parameter file to the other CPUs
        */
       MPI_Bcast(&SinkCreationDensityCodeUnits, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
       MPI_Bcast(&SinkFormationRadius, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
       MPI_Bcast(&SinkEvolutionDumpRateYears, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);


       SinkFormationRadiusSquared = SinkFormationRadius * SinkFormationRadius;
       SinkAccretionRadius = SinkFormationRadius;
       SinkSofteningLength = SinkAccretionRadius * 0.2;
       SinkEvolutionDumpRateYears *= 31557600.0 / All.UnitTime_in_s;
       SinkTestsGiveupDensityCodeUnits = SinkCreationDensityCodeUnits*100;
       mpi_printf("SINK_PARTICLES: Creation density %g  and radius %g \n", SinkCreationDensityCodeUnits, SinkFormationRadius); 
       mpi_printf("SINK_PARTICLES: Evolution dump rate in code units %g \n", SinkEvolutionDumpRateYears);
       mpi_printf("SINK_PARTICLES: ... variables set!\n");
#ifdef DEBUG_SINK_PARTICLES
       printf("SINK PARTICLES: debug formation radius %g acc rad %g \n", SinkFormationRadius, SinkAccretionRadius);
#endif
     }
   else
     {
       /* This branch is called from init.c, once the particle structures have been read, etc.
       */
       SinksLastEvolutionDumpTime = All.Time; // this ensures that we get a dump on startup
       /* Read the histories of the sink particles... init() only gets called if RestartFlag !=1
          so it's safe to call it here (it won't overwrite the correct histories in the restart files)
       */
#ifdef DUMP_SINK_PARTICLE_INFO
       sink_histories_success = open_old_sink_file_and_read_sink_data();
#endif
     }

  }

