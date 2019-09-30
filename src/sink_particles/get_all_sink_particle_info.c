#include "../allvars.h"
#include "../proto.h"

/* PCC - 19.12.2013
	This function loads the data for all sinks in the simulation into
	external structure SinkP. Returns the number of sinks.
        The current version tries out a different model of the MPI communication
        than is common for our previous sink particle algorithms. Instead of 
        having the usual send/rec setup, where all tasks communicate the 
        information in their buffers, we've adopted a Bcast setup, where only those
        tasks that have sinks, Bcast their information to the others. The rational
        is that for many applications where these sinks will be used, only a handful
        of tasks will hold active sinks. In the future, we should allow this function
    	to implement the old communication model if it's more appropriate.
*/
int get_all_sink_particle_info(int mode)
{
  int i, ii, j;
  int *NSinksEachTask;
  int *NOffset;
  int NSinksThisTask;
  int itask, ioffset;
  int NSinkCount;
  int SinkList[1000];  

  struct unique_to_sinks
    {
      double FormationMass;
      double FormationTime;
      int FormationOrder;
      int ID;
    } *store_SinkP;


  /* Find the number of sinks on this Task and their location in
     the particle list 
  */
  for(i = NSinksThisTask = 0; i < NumPart; i++)
    {
      if(P[i].Type==5) 
        {
          SinkList[NSinksThisTask] = i; 
          NSinksThisTask++;
          printf("SINK_PARTICLES: in get_sinks, found sink particle with ID %d \n", P[i].ID);
          /* XXX-pcc: need to be careful here that we don't run out of space in the
             the SinkList array! Fix this later.
          */
        }
    }


  /* Now create an array that holds the number of sinks on all Tasks 
  */
  NSinksEachTask = mymalloc("NSinksEachTask", NTask * sizeof(int));
  MPI_Allgather(&NSinksThisTask, 1, MPI_INT, NSinksEachTask, 1, MPI_INT, MPI_COMM_WORLD );


  /* Get the total number of sinks in simulation, and work out the offsets
     for the communication of the sink data.
  */ 
  for(i = NSinksAllTasks = 0; i < NTask; i++)
    NSinksAllTasks += NSinksEachTask[i];
  
  if(mode == 0)
    {
      myfree(NSinksEachTask);
      return(NSinksAllTasks);
    }


  if(NSinksAllTasks == 0)
    {
      mpi_printf("SINK_PARTICLES: No sinks in domain yet!\n");
      myfree(NSinksEachTask);
      return(NSinksAllTasks);
    }
  mpi_printf("SINK_PARTICLES: %d sinks present in domain \n", NSinksAllTasks);

  NOffset = mymalloc("NOffset", NTask * sizeof(int));
  NOffset[0] = 0;
  for(i = 1; i < NTask; i++)
    NOffset[i] = NOffset[i-1] + NSinksEachTask[i-1];


  /* Check that we still have space in the SinkP array. In the future we could be smart
     and allow the code to dynamically reallocate. This would involve first freeing the
     'NOffset' array, and so would require the function to be called again. 
  */
  if(NSinksAllTasks > NSinkBufferSize)
    {
      mpi_printf("Ran out of sink buffer, currently set at %d \n", NSinkBufferSize);
      endrun();
    }

  /* Now make a copy of the SinkP array since we're going to overwriting it in the loop below
  */
  store_SinkP = mymalloc("store_SinkP", NSinksAllTasks * sizeof(struct unique_to_sinks));
  for(i = 0; i < NSinksAllTasks; i++)
    {
      store_SinkP[i].FormationMass = SinkP[i].FormationMass;
      store_SinkP[i].FormationTime = SinkP[i].FormationTime;
      store_SinkP[i].FormationOrder = SinkP[i].FormationOrder;
      store_SinkP[i].ID = SinkP[i].ID;
    }
#ifdef DEBUG_SINK_PARTICLES
  printf("SINK_PARTICLES: GET -- just set the store_Sink array. TASK %d \n", ThisTask);
  mpi_printf("SINK_PARTICLES: GET -- Sink formass %g formtime %g \n", store_SinkP[0].FormationMass, store_SinkP[0].FormationTime);
#endif

  /* In the loop below, we do the sink communication between the Tasks.
     Looping around each Task, we allow those containing sinks to Bcast
     their data to the other tasks.
  */
  
  for(itask = 0; itask < NTask; itask++)  
    {
      if(NSinksEachTask[itask] > 0) 
        {
          /* Ready the communication buffer for this set of sinks */
          export_SinkP = (struct global_sink_particle_data *) mymalloc("export_SinkP", NSinksEachTask[itask] * sizeof(struct global_sink_particle_data));   
          /* If this Task has sinks, it might be the one broadcasting!
              If so, then it needs to fill the communication buffer 
          */
          if(itask==ThisTask)
            {
              for(i = 0; i < NSinksThisTask; i++)
                {
                  j = SinkList[i];
                  export_SinkP[i].Pos[0] = P[j].Pos[0];
                  export_SinkP[i].Pos[1] = P[j].Pos[1];
                  export_SinkP[i].Pos[2] = P[j].Pos[2];
                  export_SinkP[i].Vel[0] = P[j].Vel[0];
                  export_SinkP[i].Vel[1] = P[j].Vel[1];
                  export_SinkP[i].Vel[2] = P[j].Vel[2];
                  export_SinkP[i].Accel[0] = P[j].GravAccel[0];
                  export_SinkP[i].Accel[1] = P[j].GravAccel[1];
                  export_SinkP[i].Accel[2] = P[j].GravAccel[2];
                  export_SinkP[i].Mass = P[j].Mass;
                  export_SinkP[i].ID = P[j].ID;
                  export_SinkP[i].HomeTask = ThisTask;
                  export_SinkP[i].Index = j;
                  for(ii = 0; ii < NSinksAllTasks; ii++)
                    {
                      /* Keep track of particle properties that are ONLY recorded on the 
                         SinkP array, such as formation mass and the order in which the
                         sinks formed. Could also keep age and stellar properties here too.
                         XXX Warning! This is a N^2 operation. If we have many sinks then this
                         may become a bottleneck and will need to be rethought (i.e. with a hash)  
                      */
                      if(P[j].ID == store_SinkP[ii].ID)
                        {
                          export_SinkP[i].FormationMass = store_SinkP[ii].FormationMass;
                          export_SinkP[i].FormationTime = store_SinkP[ii].FormationTime;
                          export_SinkP[i].FormationOrder = store_SinkP[ii].FormationOrder;
                        }
                    }
                }
            }
          /* Do the Bcast to communicate the sink data */
          MPI_Bcast(export_SinkP, NSinksEachTask[itask] * sizeof(struct global_sink_particle_data), MPI_BYTE, itask, MPI_COMM_WORLD);
          /* Now copy the communication buffer to the SinkP array using the offset info */
          for(i = 0; i < NSinksEachTask[itask]; i++)
            {
              ioffset = NOffset[itask] + i;
              SinkP[ioffset].Pos[0] = export_SinkP[i].Pos[0];
              SinkP[ioffset].Pos[1] = export_SinkP[i].Pos[1];
              SinkP[ioffset].Pos[2] = export_SinkP[i].Pos[2];
              SinkP[ioffset].Vel[0] = export_SinkP[i].Vel[0];
              SinkP[ioffset].Vel[1] = export_SinkP[i].Vel[1];
              SinkP[ioffset].Vel[2] = export_SinkP[i].Vel[2];
              SinkP[ioffset].Accel[0] = export_SinkP[i].Accel[0];
              SinkP[ioffset].Accel[1] = export_SinkP[i].Accel[1];
              SinkP[ioffset].Accel[2] = export_SinkP[i].Accel[2];
              SinkP[ioffset].Mass = export_SinkP[i].Mass;
              SinkP[ioffset].FormationMass = export_SinkP[i].FormationMass;
              SinkP[ioffset].FormationTime = export_SinkP[i].FormationTime;
              SinkP[ioffset].ID = export_SinkP[i].ID;
              SinkP[ioffset].HomeTask = export_SinkP[i].HomeTask;
              SinkP[ioffset].Index = export_SinkP[i].Index; 
              SinkP[ioffset].FormationOrder = export_SinkP[i].FormationOrder;
            }
          /* Clear up the memory */
          myfree(export_SinkP);
        }
    }

  /* De-allocate anything used above (needs to be done in the reverse order
     in which they were allocated.) 
  */
  myfree(store_SinkP);
  myfree(NOffset);
  myfree(NSinksEachTask);

  /* Print the "sink densities". Might use these for free-fall time checks
     in the sink creation algorithm
  */
#ifdef DEBUG_SINK_PARTICLES_VERBOSE
  for(i = 0; i < NSinksAllTasks; i++)
    {
         mpi_printf("SINK_PARTICLES - GET: HomeTask %d ID %d pos %e %e %e mass %e form_mass %e form order %d \n", SinkP[i].HomeTask, SinkP[i].ID, SinkP[i].Pos[0], SinkP[i].Pos[1], SinkP[i].Pos[2], SinkP[i].Mass, SinkP[i].FormationMass, SinkP[i].FormationOrder);
    }
#endif

  /* Should be finished! Update the 'old' number of sinks in domain, and 
     return the number of sinks to the calling function 
   */
  NSinksAllTasksOld = NSinksAllTasks;
  return(NSinksAllTasks);
}

