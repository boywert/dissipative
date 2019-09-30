/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/riemann_hll.c
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


#if defined(RIEMANN_HLL) || (defined(RIEMANN_HLLC) && defined(DG))

#if (defined(RIEMANN_HLLC) && !defined(DG)) || defined(RIEMANN_HLLD) || defined(RIEMANN_ROSUNOV) || defined(RIEMANN_GAMMA)
#error option RIEMANN_HLL is incompatible with options RIEMANN_HLLC, RIEMANN_HLLD, RIEMANN_ROSUNOV \
       and RIEMANN_GAMMA. Only one Riemann solver can be chosen among the above options. If none   \
       of them is selected, the exact Riemann solver will be used.
#endif

static void hll_get_fluxes_from_state(struct state *st, struct fluxes *flux)
{
  flux->mass = st->rho * st->velx;
  flux->momentum[0] = st->rho * st->velx * st->velx + st->press;
  flux->momentum[1] = st->rho * st->velx * st->vely;
  flux->momentum[2] = st->rho * st->velx * st->velz;

  st->Energy = st->press / GAMMA_MINUS1 + 0.5 * st->rho * (st->velx * st->velx + st->vely * st->vely + st->velz * st->velz);
  flux->energy = (st->Energy + st->press) * st->velx;
}

static void face_get_hll_state(const struct state *st_L, const struct state *st_R, struct state_face *st_face, const struct fluxes *flux_L, const struct fluxes *flux_R, double S_L, double S_R)
{
  if(S_L >= 0.0)
    {
      st_face->press = st_L->press;
      st_face->rho = st_L->rho;
      st_face->velx = st_L->velx;
      st_face->vely = st_L->vely;
      st_face->velz = st_L->velz;
    }
  else if(S_R <= 0.0)
    {
      st_face->press = st_R->press;
      st_face->rho = st_R->rho;
      st_face->velx = st_R->velx;
      st_face->vely = st_R->vely;
      st_face->velz = st_R->velz;
    }
  else
    {
      double fac = 1.0 / (S_R - S_L);
      st_face->rho = fac * (flux_L->mass - flux_R->mass - S_L * st_L->rho + S_R * st_R->rho);

      st_face->velx = fac * (flux_L->momentum[0] - flux_R->momentum[0] - S_L * st_L->rho * st_L->velx + S_R * st_R->rho * st_R->velx) / st_face->rho;
      st_face->vely = fac * (flux_L->momentum[1] - flux_R->momentum[1] - S_L * st_L->rho * st_L->vely + S_R * st_R->rho * st_R->vely) / st_face->rho;
      st_face->velz = fac * (flux_L->momentum[2] - flux_R->momentum[2] - S_L * st_L->rho * st_L->velz + S_R * st_R->rho * st_R->velz) / st_face->rho;

      double face_energy = fac * (flux_L->energy - flux_R->energy - S_L * st_L->Energy + S_R * st_R->Energy);
      st_face->press = GAMMA_MINUS1 * (face_energy - 0.5 * st_face->rho * (st_face->velx * st_face->velx + st_face->vely * st_face->vely + st_face->velz * st_face->velz));
    }
}

static void get_hll_fluxes(const struct state *st_L, const struct state *st_R, const struct fluxes *flux_L, const struct fluxes *flux_R, struct fluxes *hll_flux, double S_L, double S_R)
{
  /* upwinding */
  S_L = dmin(S_L, 0.0);
  S_R = dmax(S_R, 0.0);

  double fac = 1.0 / (S_R - S_L);

  hll_flux->mass = fac * (S_R * flux_L->mass - S_L * flux_R->mass + S_R * S_L * (st_R->rho - st_L->rho));

  hll_flux->momentum[0] = fac * (S_R * flux_L->momentum[0] - S_L * flux_R->momentum[0] + S_R * S_L * (st_R->rho * st_R->velx - st_L->rho * st_L->velx));

  hll_flux->momentum[1] = fac * (S_R * flux_L->momentum[1] - S_L * flux_R->momentum[1] + S_R * S_L * (st_R->rho * st_R->vely - st_L->rho * st_L->vely));

  hll_flux->momentum[2] = fac * (S_R * flux_L->momentum[2] - S_L * flux_R->momentum[2] + S_R * S_L * (st_R->rho * st_R->velz - st_L->rho * st_L->velz));

  hll_flux->energy = fac * (S_R * flux_L->energy - S_L * flux_R->energy + S_R * S_L * (st_R->Energy - st_L->Energy));
}

double godunov_flux_3d_hll(struct state *st_L, struct state *st_R, struct state_face *st_face, struct fluxes *flux)
{
  if(st_L->rho > 0 && st_R->rho > 0)
    {
      double S_L, S_R;
      struct fluxes flux_L, flux_R;

      st_L->csnd = sqrt(GAMMA * st_L->press / st_L->rho);
      st_R->csnd = sqrt(GAMMA * st_R->press / st_R->rho);

      /* first calculate the wave speeds */
      S_L = dmin(st_L->velx - st_L->csnd, st_R->velx - st_R->csnd);
      S_R = dmax(st_L->velx + st_L->csnd, st_R->velx + st_R->csnd);

      /* compute fluxes for the left and right states */
      hll_get_fluxes_from_state(st_L, &flux_L);
      hll_get_fluxes_from_state(st_R, &flux_R);

      /* set the primitive variables at the face (taken from solution computed by the solver) */
      face_get_hll_state(st_L, st_R, st_face, &flux_L, &flux_R, S_L, S_R);

      /* compute HLL fluxes */
      get_hll_fluxes(st_L, st_R, &flux_L, &flux_R, flux, S_L, S_R);
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
