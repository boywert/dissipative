/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/riemann_rosunov.c
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

#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>

#include "allvars.h"
#include "proto.h"
#include "voronoi.h"


#if defined(RIEMANN_ROSUNOV) || (defined(RIEMANN_HLLC) && defined(DG))

#if (defined(RIEMANN_HLLC) && !defined(DG)) || defined(RIEMANN_HLLD) || defined(RIEMANN_HLL) || defined(RIEMANN_GAMMA)
#error option RIEMANN_ROSUNOV is incompatible with options RIEMANN_HLLC, RIEMANN_HLLD, RIEMANN_HLL \
       and RIEMANN_GAMMA. Only one Riemann solver can be chosen among the above options. If none   \
       of them is selected, the exact Riemann solver will be used.
#endif

static void rosunov_get_fluxes_from_state(struct state *st, struct fluxes *flux)
{
  flux->mass = st->rho * st->velx;
  flux->momentum[0] = st->rho * st->velx * st->velx + st->press;
  flux->momentum[1] = st->rho * st->velx * st->vely;
  flux->momentum[2] = st->rho * st->velx * st->velz;

  st->Energy = st->press / GAMMA_MINUS1 + 0.5 * st->rho * (st->velx * st->velx + st->vely * st->vely + st->velz * st->velz);
  flux->energy = (st->Energy + st->press) * st->velx;
}

static void face_get_rosunov_state(const struct state *st_L, const struct state *st_R, struct state_face *st_face, const struct fluxes *flux_L, const struct fluxes *flux_R, double S)
{
  st_face->rho = 0.5 * (flux_L->mass - flux_R->mass + S * (st_L->rho + st_R->rho)) / S;

  st_face->velx = 0.5 * (flux_L->momentum[0] - flux_R->momentum[0] + S * (st_L->rho * st_L->velx + st_R->rho * st_R->velx)) / (S * st_face->rho);
  st_face->vely = 0.5 * (flux_L->momentum[1] - flux_R->momentum[1] + S * (st_L->rho * st_L->vely + st_R->rho * st_R->vely)) / (S * st_face->rho);
  st_face->velz = 0.5 * (flux_L->momentum[2] - flux_R->momentum[2] + S * (st_L->rho * st_L->velz + st_R->rho * st_R->velz)) / (S * st_face->rho);

  double face_energy = 0.5 * (flux_L->energy - flux_R->energy + S * (st_L->Energy + st_R->Energy)) / S;
  st_face->press = GAMMA_MINUS1 * (face_energy - 0.5 * st_face->rho * (st_face->velx * st_face->velx + st_face->vely * st_face->vely + st_face->velz * st_face->velz));
}

static void get_rosunov_fluxes(const struct state *st_L, const struct state *st_R, const struct fluxes *flux_L, const struct fluxes *flux_R, struct fluxes *rosunov_flux, double S)
{
  rosunov_flux->mass = 0.5 * (flux_L->mass + flux_R->mass - S * (st_R->rho - st_L->rho));

  rosunov_flux->momentum[0] = 0.5 * (flux_L->momentum[0] + flux_R->momentum[0] - S * (st_R->rho * st_R->velx - st_L->rho * st_L->velx));

  rosunov_flux->momentum[1] = 0.5 * (flux_L->momentum[1] + flux_R->momentum[1] - S * (st_R->rho * st_R->vely - st_L->rho * st_L->vely));

  rosunov_flux->momentum[2] = 0.5 * (flux_L->momentum[2] + flux_R->momentum[2] - S * (st_R->rho * st_R->velz - st_L->rho * st_L->velz));

  rosunov_flux->energy = 0.5 * (flux_L->energy + flux_R->energy - S * (st_R->Energy - st_L->Energy));
}

double godunov_flux_3d_rosunov(struct state *st_L, struct state *st_R, struct state_face *st_face, struct fluxes *flux)
{

  if(st_L->rho > 0 && st_R->rho > 0)
    {
      double S_plus;
      struct fluxes flux_L, flux_R;

      st_L->csnd = sqrt(GAMMA * st_L->press / st_L->rho);
      st_R->csnd = sqrt(GAMMA * st_R->press / st_R->rho);

      /* first calculate the wave speeds */
      S_plus = dmax(dmax(fabs(st_L->velx - st_L->csnd), fabs(st_R->velx - st_R->csnd)), dmax(fabs(st_L->velx + st_L->csnd), fabs(st_R->velx + st_R->csnd)));

      /* compute fluxes for the left and right states */
      rosunov_get_fluxes_from_state(st_L, &flux_L);
      rosunov_get_fluxes_from_state(st_R, &flux_R);

      /* set the primitive variables at the face (taken from star region of the solver) */
      face_get_rosunov_state(st_L, st_R, st_face, &flux_L, &flux_R, S_plus);

      /* compute Rosunov fluxes */
      get_rosunov_fluxes(st_L, st_R, &flux_L, &flux_R, flux, S_plus);
    }
  else
    {
      printf("Left:  st_L->press=%g st_L->rho=%g  st_L->velx=%g\n", st_L->press, st_L->rho, st_L->velx);
      printf("Right: st_R->press=%g st_R->rho=%g  st_R->velx=%g\n", st_R->press, st_R->rho, st_R->velx);
      terminate("density is zero\n");
      return 0;
    }

  return st_face->press;
}

#endif
