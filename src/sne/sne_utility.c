/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sne/sne_utility.c
 * \date        MM/YYYY
 * \author      Robin Tress, Rowan Smith, Andre Bubel 
 * \brief       This file holds all the minor functions used throughout the SNe feedback routine workflow
 * \details     Please contact the authors at robin.tress@uni-heidelberg.de
 *              and rowan.smith@manchester.ac.uk before using this to avoid overlapping 
 *              of projects. And please report any issue you may encounter by using this 
 *              routine.
 * 
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include "../allvars.h"
#include "../proto.h"

/*! \brief Initialize variables for the SNe module
 *
 * Note that e.g. the log file descriptor is initialized in other parts of the code.
 *
 * We only initialize only the PRNG here and use the gsl_rng_ranlxd2 PRNG algorithm
 * with the supplied seed via the SNESeed parameter.
 */
void sne_init(void)
{
  sne_rng = gsl_rng_alloc(gsl_rng_ranlxd2);
  gsl_rng_set(sne_rng, All.SNESeed);
}

/*! \brief Clean up variables for the SNe module
 *
 * Note that e.g. the log file descriptor is cleaned up in other parts of the code.
 */
void sne_destroy(void)
{
  gsl_rng_free(sne_rng);

//  if(All.SNEInjectionCriterion == 4)
//    myfree_movable(SNEPDF_radial_CDF);

} 

/*! \brief initialize here the discrete CDF of the radial distribution of the SNe.
 *         This is necessary if the user selects SNEInjectionCriterion == 4, in which case
 *         the cylindrical radius is found by inverting this CDF (in inverse_radial_CDF())
 */
void sne_initialize_radial_CDF()
{
  double radial_PDF[SNEPDF_NBINS];
  SNEPDF_radial_CDF = mymalloc_movable(&SNEPDF_radial_CDF, "SNEPDF_radial_CDF", SNEPDF_NBINS * sizeof(double));

  int i;

  double Rm = 12. * KILOPARSEC / All.UnitLength_in_cm;
  double Rd = 1.5 * KILOPARSEC / All.UnitLength_in_cm;

  double dR = SNEPDF_RMAX / SNEPDF_NBINS;

  /*define radial(cylindrical) PDF (not normalized)*/
  for(i=0; i<SNEPDF_NBINS; i++)
    {
      double R = i * dR;
      radial_PDF[i] = R * exp(-Rm/R - R/Rd); //Radial distribution of the Galactic H2 disc taken from McMillan 2017
    }

  /*find radial CDF*/
  SNEPDF_radial_CDF[0] = 0.;

  for(i=1; i<SNEPDF_NBINS; i++)
    SNEPDF_radial_CDF[i] = SNEPDF_radial_CDF[i-1] + (radial_PDF[i-1] + radial_PDF[i]) * dR / 2.;

  /*normalize shit*/
  for(i=0; i<SNEPDF_NBINS; i++)
    SNEPDF_radial_CDF[i] /= SNEPDF_radial_CDF[SNEPDF_NBINS-1];
}

/*! \brief given a normalized random number x this function returns 
 *         a variable that follows a distribution whose CDF is initialized in sne_initialize_radial_CDF()
 */
double inverse_radial_CDF(double x)
{
  int ind = 0;
  int ind1 = SNEPDF_NBINS-1;

  while((ind1-ind) > 1)
    {
      int ind3 = (ind+ind1) / 2;
      double xx = SNEPDF_radial_CDF[ind3] - x;
      if(xx < 0)
        ind = ind3;
      else
        ind1 = ind3;
    }

  double dR = SNEPDF_RMAX / SNEPDF_NBINS;
  double x1 = SNEPDF_radial_CDF[ind];
  double x2 = SNEPDF_radial_CDF[ind+1];
  double fx1 = ind * dR;
  double fx2 = (ind+1) * dR;

  double dfdx = (fx2 - fx1) / (x2 - x1);

  return fx1 + dfdx * (x - x1);
}

/*! \brief
 *        
 */
double inverse_z_CDF(double x)
{
  double x1 = 2.*x -1.;
  double zd = 0.045 * KILOPARSEC / All.UnitLength_in_cm;

  return zd * log((1+x1)/(1-x1)); /*Corresponds to a density distribution of 1/(cosh(z/(2*zd)))**2, 
                                    which is the z distribution of the H2 gas of the Galaxy, taken from McMillan2017*/
}

/*! \brief A number drawn from a exponential distribution with rate 'rate'.
 *
 * See Knuth, The Art of Computer Programming 3.4.1 (D)
 */
double random_exponential_distributed(double rate)
{
  double r = gsl_rng_uniform(sne_rng);
  return -log(1.0 - r) / rate;
}

/*! \brief  finds the temperature of the cell of index idx 
 *
 */
double find_temperature_of_cell(int idx)
{
#ifdef SGCHEM
  double yn = (SphP[idx].Density * All.UnitDensity_in_cgs) / ((1.0 + 4.0 * ABHE) * PROTONMASS);
  double gas_energy_density = SphP[idx].Utherm * SphP[idx].Density * All.UnitEnergy_in_cgs / pow(All.UnitLength_in_cm,3);

#if CHEMISTRYNETWORK == 1
  double yntot = (1.0 + ABHE - SphP[idx].TracAbund[IH2] + SphP[idx].TracAbund[IHP] + SphP[i].TracAbund[IHEP] + 2.*SphP[i].TracAbund[IHEPP]) * yn;
#elif CHEMISTRYNETWORK == 15
  double yntot = (1.0 + ABHE - SphP[idx].TracAbund[IH2] + SphP[idx].TracAbund[IHP] + SphP[i].TracAbund[IHEP]) * yn;
#else
  double yntot = (1.0 + ABHE - SphP[idx].TracAbund[IH2] + SphP[idx].TracAbund[IHP]) * yn;
#endif

  double gas_temp = GAMMA_MINUS1 * gas_energy_density / (yntot * BOLTZMANN);
#else
  double gas_temp = 0.;
#endif

  return gas_temp;
}

/*! \brief finds the total mass, volume, the mean density, internal energy, temperature and the velocity of the center of mass
 *         of given cells which in principle can reside on any processor. The indices of these cells on the local processor 
 *         are given in indices[]
 *
 */
void find_gas_properties_of_injection_region(int local_n, int indices[], double *total_mass, double *total_volume, double *mean_density, double *mean_utherm, double *mean_temp, double vel_CM[3])
{

  double local_mass = 0.0f;
  double local_volume = 0.0f;
  double local_utherm = 0.0f;
  double local_temp = 0.0f;
  double vel_CM_x, local_v_CM_x = 0.0f;
  double vel_CM_y, local_v_CM_y = 0.0f;
  double vel_CM_z, local_v_CM_z = 0.0f;

  for(int i=0; i < local_n; i++)
    {
      int idx = indices[i];
      local_mass   += P[idx].Mass;
      local_volume += SphP[idx].Volume;
      local_utherm += P[idx].Mass * SphP[idx].Utherm;
      local_v_CM_x += P[idx].Mass * P[idx].Vel[0];
      local_v_CM_y += P[idx].Mass * P[idx].Vel[1];
      local_v_CM_z += P[idx].Mass * P[idx].Vel[2];

      double gas_temp = find_temperature_of_cell(idx);
      local_temp += P[idx].Mass * gas_temp;
    }

  MPI_Allreduce(&local_mass, total_mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);//combine everything into a single Allreduce
  MPI_Allreduce(&local_volume, total_volume, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&local_v_CM_x, &vel_CM_x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&local_v_CM_y, &vel_CM_y, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&local_v_CM_z, &vel_CM_z, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&local_utherm, mean_utherm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&local_temp, mean_temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  vel_CM[0] = vel_CM_x / *total_mass;
  vel_CM[1] = vel_CM_y / *total_mass;
  vel_CM[2] = vel_CM_z / *total_mass;

  *mean_utherm /= *total_mass;
  *mean_temp /= *total_mass;

  *mean_density = *total_mass / *total_volume;

}

/*! \brief smooths out the mass, internal energy and velocity of a given amount of cells.
 *
 *  \param local_n numer of particles affected on the local processor
 *  \param indices[] indices of these particles
 *  \param mean_density new density that these particles have to be given
 *  \param mean_utherm specific internal energy that these particles have to be given
 *  \param vel_CM[3] velocity to be assigned to these particles
 *
 */
void distribute_mass(int local_n, int indices[], double mean_density, double mean_utherm, double vel_CM[3])
{
  int i;
  
  for(i = 0; i < local_n; i++)
    {
      int idx = indices[i];
      SphP[idx].Density = mean_density;
      P[idx].Mass = SphP[idx].Density * SphP[idx].Volume;
      
      P[idx].Vel[0] = vel_CM[0];
      P[idx].Vel[1] = vel_CM[1];
      P[idx].Vel[2] = vel_CM[2];
      
      SphP[idx].Momentum[0] = P[idx].Mass * P[idx].Vel[0];
      SphP[idx].Momentum[1] = P[idx].Mass * P[idx].Vel[1];
      SphP[idx].Momentum[2] = P[idx].Mass * P[idx].Vel[2];
      
      SphP[idx].Utherm = mean_utherm;
      
      double abs_Mom2 = pow(SphP[idx].Momentum[0],2) + pow(SphP[idx].Momentum[1],2) + pow(SphP[idx].Momentum[2],2);
      SphP[idx].Energy = SphP[idx].Utherm * P[idx].Mass + 0.5*abs_Mom2 / P[idx].Mass;

#ifdef MAXSCALARS
      for(int k = 0; k < N_Scalar; k++)
        *(MyFloat *) (((char *) (&SphP[idx])) + scalar_elements[k].offset_mass) = *(MyFloat *) (((char *) (&SphP[idx])) + scalar_elements[k].offset) * P[idx].Mass;
#endif
 
    }

}

/*! \brief Tries to find the radius containing All.SNEMinimalParticleNumber Sph particles and 
 *         stores the indices of these cells into a variable to be used outside of this function.
 *
 *         Initially guesses the radius based on the mean volume of the cells and then adjusts the 
 *         radius based on a bisection scheme until exactly All.SNEMinimalParticleNumber are contained 
 *         within it.
 *
 *  \param center[3] center around which the cells have to be found
 *  \param *radius pointer to the variable where the computed radius has to be stored to
 *  \param *local_n this will hold the number of local cells found 
 *  \param *total_n this will hold the total cells found on any processor, this should result equal to All.SNEMinimalParticleNumber
 *  \param **indices here we will store the indices of the local particles found
 *
 */
void find_injection_cells(double center[3], double *radius, int *local_n, int *total_n, int **indices)
{
  mpi_printf("\t\t Trying to find radius containing %d cells\n", All.SNEMinimalParticleNumber);
  int ntot = 0;
  double r1, r2, tmp_radius;

  int cycle = 0;

  r1 = 0.;
  while(1)
    {
      if(cycle == 0)  //first radius estimate based on the mean volume of all the cells
        {
          r2 = pow(All.SNEMinimalParticleNumber * (boxSize_X * boxSize_Y * boxSize_Z) / (All.TotNumGas), (1./3.)) / 2.;
          find_particles_within_a_sphere(center, r2, local_n, &ntot, indices);

          while(ntot <  All.SNEMinimalParticleNumber)
            {
              myfree(*indices);
              r1 = r2;
              r2 *= 2.;
              find_particles_within_a_sphere(center, r2, local_n, &ntot, indices);
            }
        }
      else
        {
          tmp_radius = (r1 + r2)/2.;

          find_particles_within_a_sphere(center, tmp_radius, local_n, &ntot, indices);

          if(ntot < All.SNEMinimalParticleNumber)
            r1 = tmp_radius;
          else
            r2 = tmp_radius;
        }

      if (ntot == All.SNEMinimalParticleNumber)
        break;
      else
        {
          cycle++;
          myfree(*indices);
        }

    }

  *radius = r2;
  *total_n = ntot;

  mpi_printf("\t\t Done, using radius %g pc\n", r2 * All.UnitLength_in_cm / PARSEC);

}

/*! \brief Finds the particles within a sphere
 *  
 *  remember to deallocate indices after this function is called
 */
void find_particles_within_a_sphere(double center[3], double radius, int *local_n, int *total_n, int **indices)
{
//only search on local threads
//kind of naive use of the function ngb_treefind_variable_threads(), but it seems to work
  int n, j;
  int thread_id = get_thread_num();

  Thread[thread_id].R2list = (double *) mymalloc("R2list", NumPart * sizeof(double));
  Thread[thread_id].Ngblist = (int *) mymalloc("Ngblist", NumPart * sizeof(int));

  int nfound =  ngb_treefind_variable_threads(center, radius, -1, MODE_LOCAL_PARTICLES, thread_id, 1, NULL);

  int index[nfound];
  double r2;
  *local_n = 0;

  for(int n = 0; n < nfound; n++)
    {
      j = Thread[thread_id].Ngblist[n];

      //r2 = Thread[thread_id].R2list[n];
      r2 = pow(center[0] - P[j].Pos[0], 2) + pow(center[1] - P[j].Pos[1], 2) + pow(center[2] - P[j].Pos[2], 2);
      if((r2 < radius*radius) && (P[j].ID != 0) && (P[j].Mass > 0) && (P[j].Type == 0))
        {
          index[(*local_n)++] = j;
        }
    }

  myfree(Thread[thread_id].Ngblist);
  myfree(Thread[thread_id].R2list);
  *indices = (int *) mymalloc("indices", *local_n * sizeof(int));

  memcpy(*indices, &index, *local_n * sizeof(int));
  MPI_Allreduce(local_n, total_n, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

}

#ifdef INJECT_TRACER_INTO_SN
/*! \brief injects MC tracer particles equally distributed into the injection region
 *
 */
int inject_tracer_particles(int local_n, int indices[], int total_n, enum SNE_type sne_type)
{
  if (All.SNETracerBitmask & sne_type) 
    {
      int i,j;

      int tracers_per_cell = All.SNETracersForEachSn / total_n;

      if( (All.SNETracersForEachSn % total_n) >= (All.SNETracersForEachSn/2.))
        tracers_per_cell += 1;

      int total_tr_injected = (total_n * tracers_per_cell);

      MyIDType id_num = MaxTracerID + (ThisTask * total_tr_injected) + 1;

      for(i = 0; i < local_n; i++)
        {
          for(j = 0; j < tracers_per_cell; j++)
            add_new_tracer(indices[i], id_num + j);

          id_num += tracers_per_cell;
        }

      All.N_alltracer_global += total_tr_injected;
      MaxTracerID += NTask * total_tr_injected; 

      if(local_n != 0)
        printf("INJECT_TRACER_INTO_SN: Task %d injected count = %d tracers.\n", ThisTask, local_n*tracers_per_cell);

      return total_tr_injected;
    }
  else
    return 0;
}
#endif

#ifdef INJECT_TRACER_INTO_SN
/*! \brief Adds a new tracer particle to cell
 *  
 *  \param p P_index of destination parent
 *  \param newid id to assign to the new tracer
 */
void add_new_tracer(int p, MyIDType newid)
{
  int itr;
  itr = get_free_tracer_slot();
#ifdef TRACER_MC_CHECKS
  TracerLinkedList[itr].ParentID=P[p].ID;
#endif

  TracerLinkedList[itr].ID = newid;
  add_tracer_to_parent(p, itr);

  if(TracerLinkedList[itr].ID < 0)
    terminate("INJECT_TRACER_INTO_SN: tracer id too large");

#ifdef TRACER_MC_NUM_FLUID_QUANTITIES
  for(int l = 0; l < TRACER_MC_NUM_FLUID_QUANTITIES; l++)
    TracerLinkedList[itr].fluid_quantities[l] = 0.0;
#endif
}
#endif /* INJECT_TRACER_INTO_SN */


/*! \brief Logs a supernova occurrence
 *    
 * All values are in internal units
 */   
void sne_log(double sne_pos[3], int n, int injection_scheme, double radius,
             double mean_density, double mean_temperature
#ifdef INJECT_TRACER_INTO_SN
             , int n_tracer 
#endif             
             )
{
  if(ThisTask == 0)
    {
#ifdef INJECT_TRACER_INTO_SN
      fprintf(FdSNe, "%e, %e, %e, %e, %d, %d, %e, %e, %e, %d\n", All.Time, sne_pos[0], sne_pos[1], sne_pos[2], n, injection_scheme, radius, mean_density, mean_temperature, n_tracer);
#else
      fprintf(FdSNe, "%e, %e, %e, %e, %d, %d, %e, %e, %e\n", All.Time, sne_pos[0], sne_pos[1], sne_pos[2], n, injection_scheme, radius, mean_density, mean_temperature);
#endif
      myflush(FdSNe);
    }
}

