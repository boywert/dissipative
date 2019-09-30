/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sne/sne_injection_criteria.c
 * \date        MM/YYYY
 * \author      Robin Tress, Rowan Smith, Andre Bubel
 * \brief       This file holds the functions that control where and when a new Supernova has 
 *              to be injected based on the injection criterion used
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

/*! \brief based on the injection scheme used this function decides if it is
 *         time for a new Supernova
 *
 *   returns 1 if it is time for a new SN
 *
 */
enum SNE_type is_it_time_for_a_new_sn()
{
  enum SNE_type kaboom;
  switch(All.SNEInjectionCriterion)
    {
      case 1:
        kaboom = sne_only_once_at_center_TIMING();
        break;
      case 2:
        kaboom = sne_random_TIMING();
        break;
      case 3:
        kaboom = sne_at_cluster_position_TIMING();   //sne at the position of a single cluster particle representing a stellar cluster
        break;
      case 4:
        kaboom = sne_random_in_thin_disc_TIMING();
        break;
      case 5:
        kaboom = sne_following_a_given_general_distribution_TIMING();
        break;
      case 6:
        kaboom = sne_at_cluster_position_and_random_TIMING();
        break;

      default:
        mpi_printf("Wrong injection criterion chosen, no supernova will be created\n");
        kaboom = NO_SNE;
        break;
    }

  return kaboom;
}

/*! \brief based on the injection scheme used this function decides where the SN has to be positioned.
 *         It first finds a candidate position based on the injection scheme used. Then, with find_injection_cells(), finds the 
 *         cells around this position where the energy has later to be deposited (the indices of these are stored into indices). 
 *         If nobody complains about the position (complaints can be raised by is_this_position_viable()) the function returns, 
 *         otherwise a new candidate is chosen.
 *
 *  \param sne_pos[3] variable where the position of the SN will be stored
 *  \param *radius pointer to the variable where the radius of the injection region will be stored
 *  \param *local_n will store the numer of injection particles found on the local processor
 *  \param *total_n will store the total number on any processor of injectin particles
 *  \param **indices will hold the indices of the local injection cells found
 *
 */
void determine_position_and_injection_region_of_sn(double sne_pos[3], enum SNE_type sne_type, double *radius, int *local_n, int *total_n, int **indices)
{
  int viablePosition;

  while(1)
    {
      switch(sne_type)
        {
          case SNE_ONLY_ONCE:
            sne_only_once_at_center_POSITIONING(sne_pos);
            break;
          case SNE_RANDOM:
            sne_random_POSITIONING(sne_pos);
            break;
          case SNE_CLUSTER:
            sne_at_cluster_position_POSITIONING(sne_pos);
            break;
          case SNE_IN_DISC:
            sne_random_in_thin_disc_POSITIONING(sne_pos);
            break;
          case SNE_GENERAL_DISTRIBUTION:
            sne_following_a_given_general_distribution_POSITIONING(sne_pos);
            break;

          default:
            terminate("Wrong injection criterion chosen\n");
        }

      /*find injection region*/
      find_injection_cells(sne_pos, radius, local_n, total_n, indices);
      /*is the position and the injection region fine? or should I choose a different position*/
      viablePosition = is_this_position_viable(sne_pos, *radius);
      if (viablePosition)
        {
          mpi_printf("\t found a viable position \n");
          return;
        }
      else
        myfree(*indices);
    }
}

/*! \brief based on the injeciton scheme used this function will ask if 
 *         the chosen position for the SN is fine.
 *
 *  returns 1 if the position is ok, 0 otherwise.
 *
 *  \param sne_pos[3] position of the supernova
 *  \param radius radius of the injection region
 *
 */
int is_this_position_viable(double sne_pos[3], double radius)
{
  int viable_position = 1;
  switch(All.SNEInjectionCriterion)
    {
      case 1: //sne_random
#ifndef PERIODIC
      //checks if sphere touches the box boundaries 

        if   ( ((sne_pos[0] - radius) < 0.0)         || ((sne_pos[1] - radius) < 0.0)         || ((sne_pos[2] - radius) < 0.0)
          || ( ((sne_pos[0] + radius) > All.BoxSize) || ((sne_pos[1] + radius) > All.BoxSize) || ((sne_pos[2] + radius) > All.BoxSize) )
          viable_position = 0;
        else viable_position = 1;
#endif
        break;

      case 3:
        viable_position = sne_random_in_thin_disc_VALIDATE_POSITION(sne_pos, radius);
        break;

      default:
        viable_position = 1; 
    }
  return viable_position;
}

/*! \brief decides when it is time for a new supernova
 *         given the injection criterion All.SNEInjectionCriterion = 0 (Single SN at t = 0 at the center of the box) 
 *         i.e. at t=0
 *
 *  returns 1 if it is time for a new SN
 */
enum SNE_type sne_only_once_at_center_TIMING()
{
  // Flag: Have we run yet?
  static int run_once = 0;
  // Only run once and only if time == 0.0
  if(run_once == 0)
    run_once = 1;
  else
    return NO_SNE;

  if(All.Time > 0.0)
    return NO_SNE;

  return SNE_ONLY_ONCE;
}

/*! \brief decides where to inject the new supernova
 *         given the injection criterion All.SNEInjectionCriterion = 0 (Single SN at t = 0 at the center of the box)
 *         i.e. at the center of the box.
 *
 */
void sne_only_once_at_center_POSITIONING(double sne_pos[3])
{
  sne_pos[0] = boxHalf_X;
  sne_pos[1] = boxHalf_Y;
  sne_pos[2] = boxHalf_Z;
}

/*! \brief decides when it is time for a new supernova
 *         given the injection criterion All.SNEInjectionCriterion = 1 (SNe randomly distributed in space, with a given rate in time).
 *         i.e. a new SN is due if we exeed a given time wich is updated based on a distribution peaked around a given rate.
 *
 */
enum SNE_type sne_random_TIMING()
{
  double rate = 1. / ( All.SNEPeriodInYears * SEC_PER_YEAR / All.UnitTime_in_s);

  /* if we are restarting from a snapshot, skip the SNe at time 0.0*/
  if(RestartFlag == 2 && SNERandomNextTime == 0.0)
    SNERandomNextTime = All.Time + random_exponential_distributed(rate);

  if(All.Time >= SNERandomNextTime)
    {
      SNERandomNextTime += random_exponential_distributed(rate);
      mpi_printf("SNe: Next supernova around t = %g \n", SNERandomNextTime);
      return SNE_RANDOM;
    }

  return NO_SNE;
}

/*! \brief decides where to inject the new supernova
 *         given the injection criterion All.SNEInjectionCriterion = 1 (SNe randomly distributed in space, with a given rate in time).
 *         I.e. randomly within the domain.
 *
 */
void sne_random_POSITIONING(double sne_pos[3])
{
  sne_pos[0] = boxSize_X * gsl_rng_uniform(sne_rng);
  sne_pos[1] = boxSize_Y * gsl_rng_uniform(sne_rng);
  sne_pos[2] = boxSize_Z * gsl_rng_uniform(sne_rng);
}

/*! \brief 
 *
 */
enum SNE_type sne_at_cluster_position_TIMING()
{
#ifdef CLUSTERED_SNE
  /*supernovae are equally spaced in time*/
  double deltat = (All.SNEClusterTfin - All.SNEClusterTinit) / (All.SNENumber - 1);
  /*number of supernovae that already exploded*/
  int nt = ceil((All.Time - All.SNEClusterTinit) / deltat);
 
  /*if not already set, set the next time for supernova injection*/
  if(SNEClusterNextTime == 0.0)
      SNEClusterNextTime = All.SNEClusterTinit + ((nt < 0) ? 0. : nt) * deltat;

  if(All.Time >= SNEClusterNextTime && All.Time < (All.SNEClusterTfin + deltat))
    {
      SNEClusterNextTime += deltat;
      return SNE_CLUSTER;
    }

#else
  terminate("If you want SN feedback from a single stellar cluster recompile the code with CLUSTERED_SNE, otherwise change injection criterion.\n");
#endif

  return NO_SNE;
}

/*! \brief 
 *
 */
void sne_at_cluster_position_POSITIONING(double sne_pos[3])
{
  double sink_pos[3];
      
  //CAREFUL: will work only if you have a single cluster and if that is only on a single task, how it should be
  if(NumPart != NumGas) //i.e. the cluster is on this task 
    {
      sink_pos[0] = P[NumPart-1].Pos[0];  //the sink particle index is the last in the structure
      sink_pos[1] = P[NumPart-1].Pos[1];
      sink_pos[2] = P[NumPart-1].Pos[2];
      
      for(int task = 0; task < NTask; task++)
        if(task != ThisTask)
          MPI_Send(&sink_pos, 3, MPI_DOUBLE, task, 0, MPI_COMM_WORLD);
    }
  else
    MPI_Recv(&sink_pos, 3, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  sne_pos[0] = sink_pos[0];
  sne_pos[1] = sink_pos[1];
  sne_pos[2] = sink_pos[2];
  
}

/*! \brief 
 *
 */
enum SNE_type sne_random_in_thin_disc_TIMING()
{
  if (sne_random_TIMING())
    return SNE_IN_DISC;

  return NO_SNE;
}

/*! \brief 
 *
 */
void sne_random_in_thin_disc_POSITIONING(double sne_pos[3])
{

  sne_pos[0] = boxSize_X * gsl_rng_uniform(sne_rng);
  sne_pos[1] = boxSize_Y * gsl_rng_uniform(sne_rng);

  double sne_height= log((100.*PARSEC/All.UnitLength_in_cm)/gsl_rng_uniform(sne_rng));

  double sign=1.0;
  if (gsl_rng_uniform(sne_rng) > 0.5) sign=-1.0;
    sne_pos[2] = boxHalf_Z + sign*sne_height;
        
}

/*! \brief 
 *
 */
int sne_random_in_thin_disc_VALIDATE_POSITION(double sne_pos[3], double radius)
{
  double rdisc_inner=4.0*KILOPARSEC/All.UnitLength_in_cm;
  double rdisc_outer=10.0*KILOPARSEC/All.UnitLength_in_cm;

  /*Add the bubble size to the supernova central position*/
  double xminus = sne_pos[0] - radius;
  double xplus =  sne_pos[0] + radius;
  double yminus = sne_pos[1] - radius;
  double yplus =  sne_pos[1] + radius;
  double zminus = sne_pos[2] - radius;
  double zplus =  sne_pos[2] + radius;
  double galactic_radius = sqrt(pow(sne_pos[0] - boxHalf_X, 2) + pow(sne_pos[1] - boxHalf_Y, 2));

  /* Ensure SNe only go off within central 500 pc in the z-direction for computational reasons- this should be rare anyway due to the exponential function.*/

  double snzmin = boxHalf_Z - 250.*PARSEC/All.UnitLength_in_cm;
  double snzmax = boxHalf_Z + 250.*PARSEC/All.UnitLength_in_cm;

  if ( xminus < 0.0 || yminus < 0.0 || zminus < snzmin || xplus > boxSize_X || yplus > boxSize_Y || zplus > snzmax || galactic_radius < rdisc_inner*1.05 || galactic_radius > rdisc_outer*0.95)
    return 0; /*Rejected*/
  else
    return 1; /*Accepted*/

}

/*! \brief 
 *
 */
enum SNE_type sne_following_a_given_general_distribution_TIMING()
{
  if (sne_random_TIMING())
    return SNE_GENERAL_DISTRIBUTION;

  return NO_SNE;
}

/*! \brief positions of the supernovae follow a general distribution.
 *         The CDF of the radial and z distribution are defined elsewhere
 *         (see inverse_radial_CDF() and inverse_z_CDF() in sne_utility.c)
 */
void sne_following_a_given_general_distribution_POSITIONING(double sne_pos[3])
{
  // is this the first call of this function? then initialize the inverse_radial_CDF 
  static int first_call = 1;
  if(first_call)
    {
      first_call = 0;
      sne_initialize_radial_CDF();
    } 

  double rr = inverse_radial_CDF(gsl_rng_uniform(sne_rng));
  double theta = gsl_rng_uniform(sne_rng) * 2. * M_PI;
  double z = inverse_z_CDF(gsl_rng_uniform(sne_rng));

  sne_pos[0] = rr * cos(theta) + boxHalf_X;
  sne_pos[1] = rr * sin(theta) + boxHalf_Y;
  sne_pos[2] = z + boxHalf_Z;
}

enum SNE_type sne_at_cluster_position_and_random_TIMING()
{
  if (sne_at_cluster_position_TIMING())
    return SNE_CLUSTER;

  if (sne_random_TIMING())
    return SNE_GENERAL_DISTRIBUTION;

  return NO_SNE;
}
