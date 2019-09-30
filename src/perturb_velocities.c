/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/perturb_velocities.c
 * \date
 * \author
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include "allvars.h"
#include "proto.h"

#ifdef PERTURB_VELOCITIES

void perturb_velocities()
{
  int idx, i;

  mpi_printf("PERTURBVELOCITIES: Perturbing velocities by %g\n", All.VelocityPerturbation);

  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      SphP[i].Momentum[0] *= 1.0 + All.VelocityPerturbation * (2.0 * get_random_number() - 1.0);
      SphP[i].Momentum[1] *= 1.0 + All.VelocityPerturbation * (2.0 * get_random_number() - 1.0);
      SphP[i].Momentum[2] *= 1.0 + All.VelocityPerturbation * (2.0 * get_random_number() - 1.0);
    }

}

#endif
