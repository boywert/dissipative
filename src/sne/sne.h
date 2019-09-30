/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sne/sne.h
 * \date        MM/YYYY
 * \author      Robin Tress, Rowan Smith, Andre Bubel
 * \brief
 * \details     Please contact the authors at robin.tress@uni-heidelberg.de
 *              and rowan.smith@manchester.ac.uk before using this to avoid overlapping 
 *              of projects. And please report any issue you may encounter by using this 
 *              routine.
 *
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */


#include <gsl/gsl_rng.h>
#include <stdio.h>

extern gsl_rng *sne_rng;      /**< a random number generator used for sne concerns, seeded the same on all tasks */
extern FILE *FdSNe;  /** Supernova logfile **/

extern double SNERandomNextTime;
#ifdef CLUSTERED_SNE 
extern double SNEClusterNextTime;
#endif

/*variables and parameters needed for SNEInjectionCriterion = 4 
 * i.e. to define the cumulative distribution function of the radial supernova distribution*/
#define SNEPDF_NBINS 1000
#define SNEPDF_RMAX 20. * KILOPARSEC / All.UnitLength_in_cm
double *SNEPDF_radial_CDF;

enum SNE_type
{
   NO_SNE = 0,

   SNE_ONLY_ONCE = 1,
   SNE_RANDOM = 2,
   SNE_CLUSTER = 4,
   SNE_IN_DISC = 8,
   SNE_GENERAL_DISTRIBUTION = 16
};


