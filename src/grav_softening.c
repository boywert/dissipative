/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/gravtree.c
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
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"
#include "domain.h"


/*! \file grav_softening.c
 *  \brief routines for setting the gravitational softening lengths
 *
 */



/*! \brief This function sets the (comoving) softening length of all particle
 *  types in the table All.SofteningTable[...].
 *
 *  A check is performed that the physical
 *  softening length is bounded by the Softening-MaxPhys values.
 */
void set_softenings(void)
{
  int i;

  if(All.ComovingIntegrationOn)
    {
      for(i = 0; i < NSOFTTYPES; i++)
        if(All.SofteningComoving[i] * All.Time > All.SofteningMaxPhys[i])
          All.SofteningTable[i] = All.SofteningMaxPhys[i] / All.Time;
        else
          All.SofteningTable[i] = All.SofteningComoving[i];
    }
  else
    {
      for(i = 0; i < NSOFTTYPES; i++)
        All.SofteningTable[i] = All.SofteningComoving[i];
    }

#ifdef ADAPTIVE_HYDRO_SOFTENING
  for(i = 0; i < NSOFTTYPES_HYDRO; i++)
    All.SofteningTable[i + NSOFTTYPES] = All.MinimumComovingHydroSoftening * pow(All.AdaptiveHydroSofteningSpacing, i);

  if(All.AdaptiveHydroSofteningSpacing < 1)
    terminate("All.AdaptiveHydroSofteningSpacing < 1");

#ifdef MULTIPLE_NODE_SOFTENING
  /* we check that type=0 has its own slot 0 in the softening types, so that only gas masses are stored there */
  if(All.SofteningTypeOfPartType[0] != 0)
    terminate("All.SofteningTypeOfPartType[0] != 0");

  for(i = 1; i < NTYPES; i++)
    if(All.SofteningTypeOfPartType[i] == All.SofteningTypeOfPartType[0])
      terminate("i=%d: All.SofteningTypeOfPartType[i] == All.SofteningTypeOfPartType[0]", i);
#endif

#endif

#ifdef REDUCE_SOFTENINGS
  /* get position of RG core and companion */
  double pos_rgcore[3], pos_companion[3];
  double vel_rgcore[3], vel_companion[3];
  double m_rgcore, m_companion;
  int j;
  for(j = 0; j < 3; j++)
    pos_rgcore[j] = pos_companion[j] = vel_rgcore[j] = vel_companion[j] = 0.0;
  m_rgcore = 0;
  m_companion = 0;

  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type == 0)
        continue;
      if(P[i].ID == ID_RGCORE)
        for(j = 0; j < 3; j++)
          pos_rgcore[j] = P[i].Pos[j];
      if(P[i].ID == ID_RGCORE + 1)
        for(j = 0; j < 3; j++)
          pos_companion[j] = P[i].Pos[j];
    }

  /* communicate everything */
  double buf[3];
  MPI_Allreduce(pos_rgcore, buf, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for(j = 0; j < 3; j++)
    pos_rgcore[j] = buf[j];
  MPI_Allreduce(pos_companion, buf, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for(j = 0; j < 3; j++)
    pos_companion[j] = buf[j];

  double distance;
  distance = sqrt((pos_rgcore[0] - pos_companion[0]) * (pos_rgcore[0] - pos_companion[0]) +
                  (pos_rgcore[1] - pos_companion[1]) * (pos_rgcore[1] - pos_companion[1]) + (pos_rgcore[2] - pos_companion[2]) * (pos_rgcore[2] - pos_companion[2]));
  /* Check if force softening is larger than 1/5 of distance 
   * of RG core to companion */
  double reduce_factor = 2.5;
  if(2.8 * All.SofteningTable[1] > distance / reduce_factor)
    {
      All.SofteningTable[1] = distance / reduce_factor / 2.8;
      All.SofteningComoving[1] = distance / reduce_factor / 2.8;
      mpi_printf("SOFTENINGS: Reducing Type 1 softenings to %g\n", All.SofteningTable[1]);
    }
  if(2.8 * All.SofteningTable[2] > distance / reduce_factor)
    {
      All.SofteningTable[2] = distance / reduce_factor / 2.8;
      All.SofteningComoving[2] = distance / reduce_factor / 2.8;
      mpi_printf("SOFTENINGS: Reducing Type 2 softenings to %g\n", All.SofteningTable[2]);
    }
#endif

  for(i = 0; i < NSOFTTYPES + NSOFTTYPES_HYDRO; i++)
    All.ForceSoftening[i] = 2.8 * All.SofteningTable[i];

  All.ForceSoftening[NSOFTTYPES + NSOFTTYPES_HYDRO] = 0;        /* important - this entry is actually used */

}

#ifdef ADAPTIVE_HYDRO_SOFTENING
int get_softeningtype_for_hydro_cell(int i)
{
  double soft = All.GasSoftFactor * get_cell_radius(i);

  if(soft <= All.ForceSoftening[NSOFTTYPES])
    return NSOFTTYPES;

  int k = 0.5 + log(soft / All.ForceSoftening[NSOFTTYPES]) / log(All.AdaptiveHydroSofteningSpacing);
  if(k >= NSOFTTYPES_HYDRO)
    k = NSOFTTYPES_HYDRO - 1;

  return NSOFTTYPES + k;
}
#endif



/*! \brief Returns the default softening length (differening from the cut-off scale by 2.8) for particle type 'type'.
 *
 * \param i the index of the local particle
 * \return the softening length of particle i
 */
double get_default_softening_of_particletype(int type)
{
  return All.SofteningTable[All.SofteningTypeOfPartType[type]];
}



#ifdef INDIVIDUAL_GRAVITY_SOFTENING
int get_softening_type_from_mass(double mass)
{
  int i, min_type = -1;
  double eps = get_desired_softening_from_mass(mass);
  double min_dln = MAX_FLOAT_NUMBER;

#if defined(MULTIPLE_NODE_SOFTENING) && defined(ADAPTIVE_HYDRO_SOFTENING)
  i = 1;
#else
  i = 0;
#endif

  for(; i < NSOFTTYPES; i++)
    {
      if(All.ForceSoftening[i] > 0)
        {
          double dln = fabs(log(eps) - log(All.ForceSoftening[i]));

          if(dln < min_dln)
            {
              min_dln = dln;
              min_type = i;
            }
        }
    }
  if(min_type < 0)
    terminate("min_type < 0  mass=%g  eps=%g   All.AvgType1Mass=%g  All.ForceSoftening[1]=%g", mass, eps, All.AvgType1Mass, All.ForceSoftening[1]);

  return min_type;
}

/*! \brief Returns the softening length of Type 1
 * particles depending on the particle mass.
 *
 * \param mass particle mass
 * \return softening length for a type 1 particle of mass #mass
 */
double get_desired_softening_from_mass(double mass)
{
  if(mass <= All.AvgType1Mass)
    return 2.8 * All.SofteningComoving[1];
  else
    return 2.8 * All.SofteningComoving[1] * pow(mass / All.AvgType1Mass, 1.0 / 3);
}

/*! \brief Initializes the mass dependent softening calculation for Type 1 particles
 *
 * The average mass of Type 1 particles is calculated.
 */
void init_individual_softenings(void)
{
  int i, ndm;
  double mass, masstot;
  long long ndmtot;

  for(i = 0, ndm = 0, mass = 0; i < NumPart; i++)
    if(P[i].Type == 1)
      {
        ndm++;
        mass += P[i].Mass;
      }
  sumup_large_ints(1, &ndm, &ndmtot);
  MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  All.AvgType1Mass = masstot / ndmtot;

  mpi_printf("INIT: AvgType1Mass = %g\n", All.AvgType1Mass);

  for(i = 0; i < NumPart; i++)
    {
      if(((1 << P[i].Type) & (INDIVIDUAL_GRAVITY_SOFTENING)))
        P[i].SofteningType = get_softening_type_from_mass(P[i].Mass);
    }
}
#endif
