#include "../allvars.h"
#include "../proto.h"

/* PCC - 08.01.2014
         This function deals with the creation of new sinks. 
*/
int create_sink_particles(void)
{
  int sink_creation_success;
  int idx, i, j;
  int potential_peak, still_a_peak;
  int candidate_index, candidate_task, candidate_idx;
  double dx, dy, dz, dist;
  double dvx, dvy, dvz;
  double dax, day, daz;
  double egrav, esupport;
  double density_potential;
  double dist_to_sink_boundary, tenc, tff_candidate, vrad;
  double density_min_in_region;
  struct search_buffer
    {
      double density;
      int rank;
    } search_in, search_out;

  struct candidate_info
    {
      double Pos[3];
      double Vel[3];
      double Accel[3];
      double Mass;
      double Utherm;
      double Potential;
      double Density;
    } candidate;

  struct sink_environment
    {
      double mass;
      double ekin;
      double etherm;
      double divv;
      double diva;
    } environment_in, environment_out;

#ifdef SINK_PARTICLE_TEST
  if (NSinksAllTasks > 0)
    {
      mpi_printf("Already have a sink so exiting \n");
      return(-1); 
    }
#endif

  /* Initialize
  */

  sink_creation_success = -1; 

  /* Find our candidate sink for this timestep on THIS task 
     Candidate needs to be above the sink creation threshold
     density, and be the lowest potential on this Task
  */
  search_in.density = -1;
  search_in.rank = ThisTask;
  candidate_index = -1;
  candidate_idx = -1;
  density_potential = 1;
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;
      if(P[i].Type > 0 || P[i].ID == 0 || P[i].Mass == 0)
        continue;
#ifndef SINK_PARTICLE_TEST
      if(SphP[i].Density > SinkCreationDensityCodeUnits)
#endif
        {
          if(P[i].Potential < density_potential)
            {
              search_in.density = SphP[i].Density;
              candidate_index = i;
              candidate_idx = idx;
              density_potential = P[i].Potential;
            }
        }
    }

#ifdef DEBUG_SINK_PARTICLES_VERBOSE
        printf("SINK_PARTICLES: Task %d max density %g potential %g \n", ThisTask, search_in.density, density_potential);
#endif


  /* Check that the candidate isn't inside the accretion radius of another sink.
     This can be done BEFORE Bcast of candidate, since we have SinkP.
  */
  if(NSinksAllTasks > 0)
    {
      for(i = 0; i < NSinksAllTasks; i++)
        {
          dx = GRAVITY_NEAREST_X(SinkP[i].Pos[0] - P[candidate_index].Pos[0]);
          dy = GRAVITY_NEAREST_Y(SinkP[i].Pos[1] - P[candidate_index].Pos[1]);
          dz = GRAVITY_NEAREST_Z(SinkP[i].Pos[2] - P[candidate_index].Pos[2]);
          dist = dx*dx + dy*dy + dz*dz;
          if(dist < SinkFormationRadiusSquared)
            {
#ifdef DEBUG_SINK_PARTICLES_VERBOSE
              mpi_printf("SINK_PARTICLES: Sink too close on Task %d !\n", ThisTask);
#endif
              search_in.density = -1;
              candidate_index = -1; /* this is never actually used in this case, as we break immedidately*/
              break;
            }
        } 
    }



  /* Our candidate now needs to be the DENSEST of the possible potential minima 
  */
  MPI_Allreduce(&search_in, &search_out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);  

  /* If there are no candidates (no particle above rho_crit and in isolation), then we return.
  */
  if(search_out.density < 0)
    {
      mpi_printf("SINK_PARTICLES: density max %g too low or near existing sinks \n", search_out.density);
      return(sink_creation_success);
    }

  /* Check whether the particle is a real potential minimum. Bcast the position and density
     to all tasks, and let them work out whether they have a particle within the accretion
     radius that has a lower potential. They might have denser particles, but that's fine.
     Actually, we Bcast other information here too, that will be used if below, if the 
     check is passed. Might want to consider whether this is a good idea. I think it is,
     since we're anyway Bcasting a small amount of information, and one Bcast will be
     significantly cheaper than two.
  */
  candidate_task = search_out.rank;
  if(ThisTask == candidate_task)
    {
      candidate.Pos[0] = P[candidate_index].Pos[0];
      candidate.Pos[1] = P[candidate_index].Pos[1];
      candidate.Pos[2] = P[candidate_index].Pos[2];
      candidate.Vel[0] = P[candidate_index].Vel[0];
      candidate.Vel[1] = P[candidate_index].Vel[1];
      candidate.Vel[2] = P[candidate_index].Vel[2];
      candidate.Accel[0] = P[candidate_index].GravAccel[0];
      candidate.Accel[1] = P[candidate_index].GravAccel[1];
      candidate.Accel[2] = P[candidate_index].GravAccel[2];
      candidate.Mass = P[candidate_index].Mass;
      candidate.Utherm = SphP[candidate_index].Utherm;
      candidate.Potential = P[candidate_index].Potential;
      candidate.Density = SphP[candidate_index].Density;
    }

  MPI_Bcast(&candidate, 13, MPI_DOUBLE, candidate_task, MPI_COMM_WORLD);

  density_min_in_region = candidate.Density;  
  potential_peak = 1; /* Starts off 'true' */
  for(i = 0; i < NumPart; i++)
    {
      /* Since we've already checked for the sinks here, we're just looking for the
         gas. DM doesn't count here 
      */
      if(P[i].Type != 0)
        continue;
      dx = GRAVITY_NEAREST_X(candidate.Pos[0] - P[i].Pos[0]);
      dy = GRAVITY_NEAREST_Y(candidate.Pos[1] - P[i].Pos[1]);
      dz = GRAVITY_NEAREST_Z(candidate.Pos[2] - P[i].Pos[2]);
      dist = dx * dx  +  dy * dy  +  dz * dz; 
      if (dist < SinkFormationRadiusSquared)  
        {
          // could be the same particle!
          if((ThisTask == candidate_task)  &&  (i == candidate_index))
            continue;
          // Might not be a real 'peak' in potential space
          if (SphP[i].Density < density_min_in_region)
            density_min_in_region = SphP[i].Density;
          if (P[i].Potential < candidate.Potential) 
            {
#ifdef DEBUG_SINK_PARTICLES
              printf("SINK_PARTICLES: Candidate failed potential check. ThisTask %d Candidate task %d \n", ThisTask, candidate_task);
              printf("SINK_PARTICLES: Candidate failed potential check. This pot %g and rho %g and TimeBin %d. Candidate pot %g and density %g \n", P[i].Potential,SphP[i].Density, P[i].TimeBinHydro, candidate.Potential, candidate.Density);
#endif
              potential_peak = 0;
              break; 
            }
        }
    } 

  MPI_Allreduce(&potential_peak, &still_a_peak, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 

  if(ThisTask == candidate_task)
    printf("SINK_PARTICLES: Considering particle %d for sink. Potential %g Density %g accels %g %g %g \n", P[candidate_index].ID, P[candidate_index].Potential, SphP[candidate_index].Density, P[candidate_index].GravAccel[0], P[candidate_index].GravAccel[1], P[candidate_index].GravAccel[2]);

  if(still_a_peak < NTask  &&  candidate.Density < SinkTestsGiveupDensityCodeUnits)
    {
      sink_creation_success = -2;
      return(sink_creation_success);
    }

#ifndef SINK_PARTICLE_TEST

  /* Now check that the region inside the candidate's interaction radius is
     indeed both bound (energy check) and collapsing (divv & diva checks). 
     Once again, each CPU will run in parallel. Each task bins its gas particles
     radially, centred on the candidate sink. We then do a MPI_All sum to get
     the full radial profile. This is then used to make an energy profile, which
     will be used to determine whether this gas cell can become a sink. If so,
     the sink details are then added to the SinkP array immediately, and obviously
     a new particle is created on the host task with P.Type = 5. 
     The newly-formed sink is then given a chance to accrete 
  */

  memset(&environment_in, 0, sizeof(struct sink_environment));

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(P[i].Mass == 0 || P[i].ID == 0)
        continue;
      dx = GRAVITY_NEAREST_X(P[i].Pos[0] - candidate.Pos[0]);
      dy = GRAVITY_NEAREST_Y(P[i].Pos[1] - candidate.Pos[1]);
      dz = GRAVITY_NEAREST_Z(P[i].Pos[2] - candidate.Pos[2]);
      dist = dx * dx  +  dy * dy  +  dz * dz;
      if(dist > SinkFormationRadiusSquared)
        continue;
      if(ThisTask == candidate_task  &&  i == candidate_index) 
        continue;
      dist = sqrt(dist);
      dvx = P[i].Vel[0] - candidate.Vel[0];
      dvy = P[i].Vel[1] - candidate.Vel[1];
      dvz = P[i].Vel[2] - candidate.Vel[2];
      dax = P[i].GravAccel[0] - candidate.Accel[0];
      day = P[i].GravAccel[1] - candidate.Accel[1];
      daz = P[i].GravAccel[2] - candidate.Accel[2];
      /* Now work out the partial sums (mass weighting!) for ekin, etherm 
         divv, and diva
      */
      environment_in.mass += P[i].Mass;
      environment_in.ekin += 0.5 * (dvx * dvx  +  dvy * dvy  +  dvz * dvz) * P[i].Mass;
      environment_in.etherm += SphP[i].Utherm * P[i].Mass;
      environment_in.divv += P[i].Mass * (dvx * dx  +  dvy * dy  +  dvz * dz) / dist;
      environment_in.diva += P[i].Mass * (dax * dx  +  day * dy  +  daz * dz) / dist;
    }


  /* Do the MPI_All.. Sum to get the full values.
  */
  memset(&environment_out, 0, sizeof(struct sink_environment));
  MPI_Allreduce(&environment_in, &environment_out, 5, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);



  /* Now get what we need from the sums. Need to add the candiate to the mass and 
     the thermal energy calculation. We don't need to add it to the kinetic energy though
  */
  environment_out.mass += candidate.Mass;
  environment_out.etherm += (candidate.Mass * candidate.Utherm); 
  environment_out.divv /= environment_out.mass;
  environment_out.diva /= environment_out.mass; 
  egrav = environment_out.mass * environment_out.mass / SinkFormationRadius; 
  esupport = environment_out.ekin + environment_out.etherm;


  /* Is this candidate the centre of a collapsing structure? If so, we have a sink,
     if not, we return with the appropriate value of the success flag.
  */
  if (candidate.Density > SinkTestsGiveupDensityCodeUnits) 
    {
       mpi_printf("SINK_PARTICLES: Candidate density of %g is now very dense. Forcing sink formation! \n", candidate.Density);
    }
  else
    {
      sink_creation_success = -3;
      if((egrav <= 2.*esupport)  ||  (environment_out.divv >= 0)  ||  (environment_out.diva >= 0))
        {
          mpi_printf("SINK_PARTICLES: Candidate failed the energy checks: egrav %g esupport %g divv %g diva %g \n", egrav, esupport, environment_out.divv, environment_out.diva);
          return(sink_creation_success);
        }
    }
 
#ifdef SINK_PARTICLE_FREE_FALL_TEST 
  /* Final test! Check that the candidate isn't going to fall into another sink's accretion radius
     before its is able to collapse locally: i.e. tff(candidate) < interaction time
     This test is useful if you have discs forming around existing sinks, and there's nothing to
     stabilise the inner disc, or if you have very dense clusters.
  */
  if (NSinksAllTasks > 0  &&  candidate.Density < SinkTestsGiveupDensityCodeUnits)
    {
#ifdef DEBUG_SINK_PARTICLES
      mpi_printf("SINK_PARTICLES--FORM: checking the free fall time criterion \n");
#endif
      for(i = 0; i < NSinksAllTasks; i++)
        {
          dx = GRAVITY_NEAREST_X(candidate.Pos[0] - SinkP[i].Pos[0]);
          dy = GRAVITY_NEAREST_Y(candidate.Pos[1] - SinkP[i].Pos[1]);
          dz = GRAVITY_NEAREST_Z(candidate.Pos[2] - SinkP[i].Pos[2]);
          dist = sqrt(dx*dx + dy*dy + dz*dz);
          dist_to_sink_boundary  = dist - SinkFormationRadius;     
          if (dist_to_sink_boundary > 0)          
            {
              dvx = candidate.Vel[0] - SinkP[i].Vel[0];
              dvy = candidate.Vel[1] - SinkP[i].Vel[1];
              dvz = candidate.Vel[2] - SinkP[i].Vel[2];          
              vrad = (dvx*dx + dvy*dy + dvz*dz) / dist;
#ifdef DEBUG_SINK_PARTICLES
              mpi_printf("SINK_PARTICLES--FORM: free-fall check dist %g dist to sink %g vrad %g \n", dist, dist_to_sink_boundary, vrad);
#endif
              if (vrad <= 0)
                {
                  tenc = dist_to_sink_boundary / (-vrad);
                  tff_candidate = sqrt(3 * M_PI / 32.0/ All.G/ candidate.Density);
#ifdef DEBUG_SINK_PARTICLES
                  mpi_printf("SINK_PARTICLES--FORM: free fall check tff %g tenc %g \n", tff_candidate, tenc);
#endif
                  if (tenc <  tff_candidate)
                    {
#ifdef DEBUG_SINK_PARTICLES
                      mpi_printf("SINK_PARTICLES--FORM: candidate will fall into sink before formation on Task %d \n", ThisTask);
                      mpi_printf("SINK_PARTICLES--FORM: dist_to_sink_boundary %g tenc %g tff_candidate %g \n", dist_to_sink_boundary, tenc, tff_candidate);
#endif
                      sink_creation_success = -4;
                      return(sink_creation_success);
                    }
                }
            }
          else
            { 
              mpi_printf("SINK_PARTICLES--FORM: candidate is inside another sink on Task %d! This shouldn't happen! Aborting\n", ThisTask);
              endrun();
            }
        }
    } 
#endif

#endif

  /* Joy! Bang a drum...
  */ 
  mpi_printf("SINK_PARTICLES: SINK CREATION SUCCESSFUL!\n");
  mpi_printf("SINK_PARTICLES: CREATION egrav %g esupport %g divv %g diva %g \n", egrav, esupport, environment_out.divv, environment_out.diva);
  mpi_printf("SINK_PARTICLES: CREATION formation mass %g cell mass %g ", environment_out.mass, candidate.Mass);

  sink_creation_success = 1;


  /* Update the SinkP array. The ID and Index are only needed on the host Task, so they
     don't need to be Bcast to the other Tasks.
  */
  for(i = 0; i < 3; i++)
    SinkP[NSinksAllTasks].Pos[i] = candidate.Pos[i];
  for(i = 0; i < 3; i++)
    SinkP[NSinksAllTasks].Vel[i] = candidate.Vel[i];
  for(i = 0; i < 3; i++)
    SinkP[NSinksAllTasks].Accel[i] = candidate.Accel[i];
  SinkP[NSinksAllTasks].Mass = candidate.Mass;
  SinkP[NSinksAllTasks].FormationMass = environment_out.mass; /* Used in early accretion*/
  SinkP[NSinksAllTasks].FormationTime = All.Time;
  SinkP[NSinksAllTasks].HomeTask = candidate_task;
  if(ThisTask == candidate_task)
    {
      SinkP[NSinksAllTasks].ID = P[candidate_index].ID;
      SinkP[NSinksAllTasks].Index = NumPart;
    }
  else
    {
      SinkP[NSinksAllTasks].ID = -1; 
      SinkP[NSinksAllTasks].Index = -1; 
    }
  SinkP[NSinksAllTasks].FormationOrder = NSinksAllTasks + 1;

  NSinksAllTasks++;
  NSinksThisTask++;

#ifdef SINK_PARTICLE_TEST
  /* Set a random amount of mass to see how potential reacts....*/
  SinkP[NSinksAllTasks].Mass = 10.;
  SinkP[NSinksAllTasks].FormationMass = 10.;
  if(ThisTask == candidate_task)
    P[candidate_index].Mass = 10.;
#endif

  /* Convert cell to particle 
  */
  add_sink_to_particle_structure(candidate_idx, candidate_index, candidate_task);


  /* Keep track of which task was responsible for forming this new sink... Might be useful!
  */
  LastTaskNewSink = candidate_task;  

  return(sink_creation_success);
}





/* This function turns the gas cell ('particle' Type 0) into a sink particle (Type 5).
   The code is based on the technique in star_formation.c (functions make_star and 
   convert_cell_into_star). Unlike the old version, where we made a new Type 5 particle, 
   and then ripped the Type 0 particle out of the mesh, here we simply change its type!
   We also free up the memory that the particle is taking in the mesh, and remove it from
   the list of particle in the hydro time bin. 
*/
void add_sink_to_particle_structure(int idx, int candidate_index, int candidate_task)
  {
    mpi_printf("SINK_PARTICLES: Adding particle to P array and cleaning up! \n");

    if(ThisTask == candidate_task) 
      {   
        P[candidate_index].Type = 5;
        P[candidate_index].SofteningType = All.SofteningTypeOfPartType[P[candidate_index].Type];
#ifdef INDIVIDUAL_GRAVITY_SOFTENING
        if(((1 << P[candidate_index].Type) & (INDIVIDUAL_GRAVITY_SOFTENING)))
          {
            P[candidate_index].SofteningType = get_softening_type_from_mass(P[candidate_index].Mass);
            printf("SINK_PARTICLES: Doing the weird mass softening thing... \n");
          } 
#endif

        printf("SINK_PARTICLES: Task %d converting particle with ID %d to sink. the idx is %d \n", ThisTask, P[candidate_index].ID, idx);

  /* Can remove the candidate gas cell from the mesh now
  */
        timebin_remove_particle(&TimeBinsHydro, idx, P[candidate_index].TimeBinHydro);
#ifdef VORONOI_DYNAMIC_UPDATE
        voronoi_remove_connection(candidate_index);
#endif
        printf("SINK_PARTICLES: Finshed the update of the P-array \n");
     }
   SinksFormedSinceLastDomain++;
 }
