/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/riemann_hlld.c
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


#if defined(RIEMANN_HLLD)

static double hlld_get_fast_wave(struct state *st);
static void hlld_get_fluxes_from_state(struct state *st_face, struct fluxes *flux, double *st_ptot);
static void hlld_get_star(struct state *st_star, struct state *st, double S, double S_M, double ptot, double ptot_star);
static void hlld_get_fluxes_star(struct state *st_A, struct state *st_A_star, struct fluxes *flux_A, double S_A, struct fluxes *flux);
static void hlld_get_starstar_L(struct state *st_star_L, struct state *st_star_R, struct state *st_starstar);
static void hlld_get_starstar_R(struct state *st_star_L, struct state *st_star_R, struct state *st_starstar);
static void hlld_get_starstar(struct state *st_star_L, struct state *st_star_R, struct state *st_starstar, struct state *st_star_A, double sign);
static void hlld_get_fluxes_starstar(struct state *st_A, struct state *st_A_star, struct state *st_A_starstar, struct fluxes *flux_A, double S_A, double S_A_star, struct fluxes *flux);
static void hll_get_star(struct state *st_star, struct fluxes *flux_L, struct fluxes *flux_R, struct state *st_L, struct state *st_R, double S_L, double S_R);
static void hll_get_flux(struct fluxes *flux, struct fluxes *flux_L, struct fluxes *flux_R, struct state *st_L, struct state *st_R, double S_L, double S_R);
static void lax_get_flux(struct fluxes *flux, struct fluxes *flux_L, struct fluxes *flux_R, struct state *st_L, struct state *st_R, double S);

static inline int state_and_flux_valid(const struct state *st, const struct fluxes *flux)
{
  return (st->press >= 0) && gsl_finite(st->press) && gsl_finite(flux->energy);
}

double godunov_flux_3d_hlld(struct state *st_L, struct state *st_R, double *vel_face, struct state_face *st_face, struct fluxes *flux)
{
  struct state st_Lstar, st_Rstar, st_star;
  struct state st_Lstarstar, st_Rstarstar;
  struct state *st_middle;
  double Bx;
  double cf_L, cf_R;
  double S, S_L, S_R, S_M, S_L_star, S_R_star;
  double ptot_L, ptot_R;

  S_R_star = S_L_star = S_M = 0.;

  if(st_L->rho > 0 && st_R->rho > 0)
    {
      Bx = 0.5 * (st_L->Bx + st_R->Bx);
      flux->B[0] = 0.;

      st_L->Bx = Bx;
      st_R->Bx = Bx;
      st_face->Bx = Bx;

      /* get wave speeds first */
      cf_L = hlld_get_fast_wave(st_L);
      cf_R = hlld_get_fast_wave(st_R);

      S = dmax(dmax(fabs(st_L->velx - cf_L), fabs(st_R->velx - cf_R)), dmax(fabs(st_L->velx + cf_L), fabs(st_R->velx + cf_R)));

      S_L = dmin(st_L->velx - cf_L, st_R->velx - cf_R);
      S_R = dmax(st_L->velx + cf_L, st_R->velx + cf_R);

      if(S_L >= 0)
        {
          st_middle = st_L;
          hlld_get_fluxes_from_state(st_L, flux, NULL);
        }
      else if(S_R <= 0)
        {
          st_middle = st_R;
          hlld_get_fluxes_from_state(st_R, flux, NULL);
        }
      else
        {
          // stars are needed
          struct fluxes flux_R, flux_L;

          hlld_get_fluxes_from_state(st_L, &flux_L, &ptot_L);
          hlld_get_fluxes_from_state(st_R, &flux_R, &ptot_R);

          S_M = ((S_R - st_R->velx) * st_R->rho * st_R->velx - (S_L - st_L->velx) * st_L->rho * st_L->velx - ptot_R + ptot_L) / ((S_R - st_R->velx) * st_R->rho - (S_L - st_L->velx) * st_L->rho);

          double ptot_star =
            ((S_R - st_R->velx) * st_R->rho * ptot_L - (S_L - st_L->velx) * st_L->rho * ptot_R +
             st_L->rho * st_R->rho * (S_R - st_R->velx) * (S_L - st_L->velx) * (st_R->velx - st_L->velx)) / ((S_R - st_R->velx) * st_R->rho - (S_L - st_L->velx) * st_L->rho);

          hlld_get_star(&st_Lstar, st_L, S_L, S_M, ptot_L, ptot_star);
          hlld_get_star(&st_Rstar, st_R, S_R, S_M, ptot_R, ptot_star);

          S_L_star = S_M - fabs(st_L->Bx) / sqrt(st_Lstar.rho);
          S_R_star = S_M + fabs(st_R->Bx) / sqrt(st_Rstar.rho);

          if(S_L_star >= 0 || (Bx == 0 && S_M >= 0))    // we already know: S_L <= 0 
            {
              st_middle = &st_Lstar;
              hlld_get_fluxes_star(st_L, &st_Lstar, &flux_L, S_L, flux);
            }
          else if(S_R_star <= 0 || (Bx == 0))   // we already know: S_R >= 0
            {
              st_middle = &st_Rstar;
              hlld_get_fluxes_star(st_R, &st_Rstar, &flux_R, S_R, flux);
            }
          else
            {
              // double stars are needed
              if(S_M >= 0)      // we already know: S_L_star <= 0)
                {
                  st_middle = &st_Lstarstar;
                  hlld_get_starstar_L(&st_Lstar, &st_Rstar, &st_Lstarstar);
                  hlld_get_fluxes_starstar(st_L, &st_Lstar, &st_Lstarstar, &flux_L, S_L, S_L_star, flux);
                }
              else              // we already know: S_R_star >= 0 and S_M <= 0
                {
                  st_middle = &st_Rstarstar;
                  hlld_get_starstar_R(&st_Lstar, &st_Rstar, &st_Rstarstar);
                  hlld_get_fluxes_starstar(st_R, &st_Rstar, &st_Rstarstar, &flux_R, S_R, S_R_star, flux);
                }
            }
        }
    }
  else
    {
      printf("Left:  st_L->press=%g st_L->rho=%g  st_L->velx=%g\n", st_L->press, st_L->rho, st_L->velx);
      printf("Right: st_R->press=%g st_R->rho=%g  st_R->velx=%g\n", st_R->press, st_R->rho, st_R->velx);
      terminate("density is zero\n");
      return 0;
    }

  if(!state_and_flux_valid(st_middle, flux))
    {
      /* HLLD did not work => use HLL instead */
      struct fluxes flux_R, flux_L;

      hlld_get_fluxes_from_state(st_L, &flux_L, NULL);
      hlld_get_fluxes_from_state(st_R, &flux_R, NULL);

      hll_get_star(&st_star, &flux_L, &flux_R, st_L, st_R, S_L, S_R);
      hll_get_flux(flux, &flux_L, &flux_R, st_L, st_R, S_L, S_R);

      st_middle = &st_star;

      if(!state_and_flux_valid(st_middle, flux))
        {
          /* HLL did not work, use lax-friedrich flux instead */
          lax_get_flux(flux, &flux_L, &flux_R, st_L, st_R, S);

          st_star.press = 0.5 * (st_L->press + st_R->press);
        }
    }

  st_face->rho = st_middle->rho;
  st_face->velx = st_middle->velx;
  st_face->vely = st_middle->vely;
  st_face->velz = st_middle->velz;
  st_face->press = st_middle->press;
  st_face->By = st_middle->By;
  st_face->Bz = st_middle->Bz;

  if(!state_and_flux_valid(st_middle, flux))
    {
      printf("M: rho=%g, v=(%g,%g,%g), p=%g, B=(%g,%g,%g)\n", st_middle->rho, st_middle->velx + vel_face[0],
             st_middle->vely + vel_face[1], st_middle->velz + vel_face[2], st_middle->press, st_middle->Bx, st_middle->By, st_middle->Bz);
      printf("S_L=%g, S_L_star=%g, S_M=%g, S_R_star=%g, S_R=%g, cf_L=%g, cf_R=%g\n", S_L, S_L_star, S_M, S_R_star, S_R, cf_L, cf_R);
    }

  return st_middle->press;
}


static double hlld_get_fast_wave(struct state *st)
{
#ifdef VARIABLE_GAMMA
  double gamma = st->gammaC;
#else
  double gamma = GAMMA;
#endif

  double gPress = gamma * st->press;
#ifdef COSMIC_RAYS
  gPress += All.GammaCR * st->crPressure;
#endif

  /* TODO: add cosmic ray pressure here? */

  double Bsqr = st->Bx * st->Bx + st->By * st->By + st->Bz * st->Bz;
  double gpb2 = gPress + Bsqr;

  return sqrt(0.5 / st->rho * (gpb2 + sqrt(gpb2 * gpb2 - 4. * gPress * st->Bx * st->Bx)));
}

static void hlld_get_fluxes_from_state(struct state *st, struct fluxes *flux, double *st_ptot)
{
#ifdef VARIABLE_GAMMA
  double gamma = st->gammaE;
#else
  double gamma = GAMMA;
#endif
  double gamma_minus1 = gamma - 1.;

#ifdef COSMIC_RAYS
  double cr_press = st->crPressure;
#else
  double cr_press = 0.;
#endif

  flux->mass = st->rho * st->velx;
  double Bsqr = st->Bx * st->Bx + st->By * st->By + st->Bz * st->Bz;
  flux->momentum[0] = st->rho * st->velx * st->velx + st->press + 0.5 * Bsqr - st->Bx * st->Bx + cr_press;
  flux->momentum[1] = st->rho * st->velx * st->vely - st->Bx * st->By;
  flux->momentum[2] = st->rho * st->velx * st->velz - st->Bx * st->Bz;

  flux->B[1] = st->By * st->velx - st->Bx * st->vely;
  flux->B[2] = st->Bz * st->velx - st->Bx * st->velz;

  double etot = st->press / gamma_minus1 + 0.5 * st->rho * (st->velx * st->velx + st->vely * st->vely + st->velz * st->velz) + 0.5 * Bsqr;
  double ptot = st->press + 0.5 * Bsqr + cr_press;

  flux->energy = (etot + ptot) * st->velx - st->Bx * (st->velx * st->Bx + st->vely * st->By + st->velz * st->Bz);

  st->Energy = etot;
  if(st_ptot)
    *st_ptot = ptot;
}

static void hlld_get_star(struct state *st_star, struct state *st, double S, double S_M, double ptot, double ptot_star)
{
  st_star->rho = st->rho * (S - st->velx) / (S - S_M);
  st_star->velx = S_M;
  st_star->vely = st->vely - st->Bx * st->By * (S_M - st->velx) / (st->rho * (S - st->velx) * (S - S_M) - st->Bx * st->Bx);
  st_star->velz = st->velz - st->Bx * st->Bz * (S_M - st->velx) / (st->rho * (S - st->velx) * (S - S_M) - st->Bx * st->Bx);

  st_star->Bx = st->Bx;
  st_star->By = st->By * (st->rho * (S - st->velx) * (S - st->velx) - st->Bx * st->Bx) / (st->rho * (S - st->velx) * (S - S_M) - st->Bx * st->Bx);
  st_star->Bz = st->Bz * (st->rho * (S - st->velx) * (S - st->velx) - st->Bx * st->Bx) / (st->rho * (S - st->velx) * (S - S_M) - st->Bx * st->Bx);

  st_star->Energy = ((S - st->velx) * st->Energy - ptot * st->velx + ptot_star * S_M
                     + st->Bx * (st->velx * st->Bx + st->vely * st->By + st->velz * st->Bz - st_star->velx * st->Bx - st_star->vely * st_star->By - st_star->velz * st_star->Bz)) / (S - S_M);

  st_star->press = ptot_star - 0.5 * (st_star->Bx * st_star->Bx + st_star->By * st_star->By + st_star->Bz * st_star->Bz);
}

static void hlld_get_fluxes_star(struct state *st_A, struct state *st_A_star, struct fluxes *flux_A, double S_A, struct fluxes *flux)
{
  flux->mass = flux_A->mass - S_A * (st_A->rho - st_A_star->rho);

  flux->momentum[0] = flux_A->momentum[0] - S_A * (st_A->rho * st_A->velx - st_A_star->rho * st_A_star->velx);
  flux->momentum[1] = flux_A->momentum[1] - S_A * (st_A->rho * st_A->vely - st_A_star->rho * st_A_star->vely);
  flux->momentum[2] = flux_A->momentum[2] - S_A * (st_A->rho * st_A->velz - st_A_star->rho * st_A_star->velz);

  flux->B[1] = flux_A->B[1] - S_A * (st_A->By - st_A_star->By);
  flux->B[2] = flux_A->B[2] - S_A * (st_A->Bz - st_A_star->Bz);

  flux->energy = flux_A->energy - S_A * (st_A->Energy - st_A_star->Energy);
}

static void hlld_get_starstar_L(struct state *st_star_L, struct state *st_star_R, struct state *st_starstar)
{
  hlld_get_starstar(st_star_L, st_star_R, st_starstar, st_star_L, -1.0);
}

static void hlld_get_starstar_R(struct state *st_star_L, struct state *st_star_R, struct state *st_starstar)
{
  hlld_get_starstar(st_star_L, st_star_R, st_starstar, st_star_R, 1.0);
}

static void hlld_get_starstar(struct state *st_star_L, struct state *st_star_R, struct state *st_starstar, struct state *st_star_A, double sign)
{
  double sBx = st_star_A->Bx < 0 ? -1.0 : 1.0;

  double sqLrho = sqrt(st_star_L->rho);
  double sqRrho = sqrt(st_star_R->rho);

  st_starstar->rho = st_star_A->rho;

  st_starstar->velx = st_star_L->velx;  /* == st_star_R->velx == S_M */
  st_starstar->vely = ((sqLrho * st_star_L->vely) + (sqRrho * st_star_R->vely) + (st_star_R->By - st_star_L->By) * sBx) / (sqLrho + sqRrho);
  st_starstar->velz = ((sqLrho * st_star_L->velz) + (sqRrho * st_star_R->velz) + (st_star_R->Bz - st_star_L->Bz) * sBx) / (sqLrho + sqRrho);

  st_starstar->Bx = st_star_A->Bx;
  st_starstar->By = ((sqLrho * st_star_R->By) + (sqRrho * st_star_L->By) + sqLrho * sqRrho * (st_star_R->vely - st_star_L->vely) * sBx) / (sqLrho + sqRrho);
  st_starstar->Bz = ((sqLrho * st_star_R->Bz) + (sqRrho * st_star_L->Bz) + sqLrho * sqRrho * (st_star_R->velz - st_star_L->velz) * sBx) / (sqLrho + sqRrho);

  st_starstar->Energy = st_star_A->Energy + sign * sqrt(st_star_A->rho) * sBx
    * (st_star_A->velx * st_star_A->Bx + st_star_A->vely * st_star_A->By + st_star_A->velz * st_star_A->Bz
       - st_starstar->velx * st_star_A->Bx - st_starstar->vely * st_starstar->By - st_starstar->velz * st_starstar->Bz);

  st_starstar->press = st_star_A->press;
}

static void hlld_get_fluxes_starstar(struct state *st_A, struct state *st_A_star, struct state *st_A_starstar, struct fluxes *flux_A, double S_A, double S_A_star, struct fluxes *flux)
{
  flux->mass = flux_A->mass + S_A_star * st_A_starstar->rho - (S_A_star - S_A) * st_A_star->rho - S_A * st_A->rho;

  flux->momentum[0] = flux_A->momentum[0] + S_A_star * st_A_starstar->rho * st_A_starstar->velx - (S_A_star - S_A) * st_A_star->rho * st_A_star->velx - S_A * st_A->rho * st_A->velx;
  flux->momentum[1] = flux_A->momentum[1] + S_A_star * st_A_starstar->rho * st_A_starstar->vely - (S_A_star - S_A) * st_A_star->rho * st_A_star->vely - S_A * st_A->rho * st_A->vely;
  flux->momentum[2] = flux_A->momentum[2] + S_A_star * st_A_starstar->rho * st_A_starstar->velz - (S_A_star - S_A) * st_A_star->rho * st_A_star->velz - S_A * st_A->rho * st_A->velz;

  flux->B[1] = flux_A->B[1] + S_A_star * st_A_starstar->By - (S_A_star - S_A) * st_A_star->By - S_A * st_A->By;
  flux->B[2] = flux_A->B[2] + S_A_star * st_A_starstar->Bz - (S_A_star - S_A) * st_A_star->Bz - S_A * st_A->Bz;

  flux->energy = flux_A->energy + S_A_star * st_A_starstar->Energy - (S_A_star - S_A) * st_A_star->Energy - S_A * st_A->Energy;
}

static void hll_get_star(struct state *st_star, struct fluxes *flux_L, struct fluxes *flux_R, struct state *st_L, struct state *st_R, double S_L, double S_R)
{
#ifdef VARIABLE_GAMMA
  st_star->gammaE = 0.5 * (st_L->gammaE + st_R->gammaE);
  double gamma = st_star->gammaE;
#else
  double gamma = GAMMA;
#endif
  double gamma_minus1 = gamma - 1.;

  double fac = 1.0 / (S_R - S_L);

  st_star->rho = fac * (S_R * st_R->rho - S_L * st_L->rho - flux_R->mass + flux_L->mass);

  st_star->velx = fac * (S_R * st_R->rho * st_R->velx - S_L * st_L->rho * st_L->velx - flux_R->momentum[0] + flux_L->momentum[0]) / st_star->rho;
  st_star->vely = fac * (S_R * st_R->rho * st_R->vely - S_L * st_L->rho * st_L->vely - flux_R->momentum[1] + flux_L->momentum[1]) / st_star->rho;
  st_star->velz = fac * (S_R * st_R->rho * st_R->velz - S_L * st_L->rho * st_L->velz - flux_R->momentum[2] + flux_L->momentum[2]) / st_star->rho;

  st_star->Energy = fac * (S_R * st_R->Energy - S_L * st_L->Energy - flux_R->energy + flux_L->energy);

  st_star->Bx = st_R->Bx;       /* == st_L->Bx */
  st_star->By = fac * (S_R * st_R->By - S_L * st_L->By - flux_R->B[1] + flux_L->B[1]);
  st_star->Bz = fac * (S_R * st_R->Bz - S_L * st_L->Bz - flux_R->B[2] + flux_L->B[2]);

  st_star->press = gamma_minus1 * (st_star->Energy
                                   - 0.5 * st_star->rho * (st_star->velx * st_star->velx +
                                                           st_star->vely * st_star->vely +
                                                           st_star->velz * st_star->velz) - 0.5 * (st_star->Bx * st_star->Bx + st_star->By * st_star->By + st_star->Bz * st_star->Bz));
}

static void hll_get_flux(struct fluxes *flux, struct fluxes *flux_L, struct fluxes *flux_R, struct state *st_L, struct state *st_R, double S_L, double S_R)
{
  double fac = 1.0 / (S_R - S_L);

  flux->mass = fac * (S_R * flux_L->mass - S_L * flux_R->mass + S_R * S_L * (st_R->rho - st_L->rho));

  flux->momentum[0] = fac * (S_R * flux_L->momentum[0] - S_L * flux_R->momentum[0] + S_R * S_L * (st_R->rho * st_R->velx - st_L->rho * st_L->velx));
  flux->momentum[1] = fac * (S_R * flux_L->momentum[1] - S_L * flux_R->momentum[1] + S_R * S_L * (st_R->rho * st_R->vely - st_L->rho * st_L->vely));
  flux->momentum[2] = fac * (S_R * flux_L->momentum[2] - S_L * flux_R->momentum[2] + S_R * S_L * (st_R->rho * st_R->velz - st_L->rho * st_L->velz));

  flux->energy = fac * (S_R * flux_L->energy - S_L * flux_R->energy + S_R * S_L * (st_R->Energy - st_L->Energy));

  flux->B[1] = fac * (S_R * flux_L->B[1] - S_L * flux_R->B[1] + S_R * S_L * (st_R->By - st_L->By));
  flux->B[2] = fac * (S_R * flux_L->B[2] - S_L * flux_R->B[2] + S_R * S_L * (st_R->Bz - st_L->Bz));
}

static void lax_get_flux(struct fluxes *flux, struct fluxes *flux_L, struct fluxes *flux_R, struct state *st_L, struct state *st_R, double S)
{

  flux->mass = 0.5 * (flux_L->mass + flux_R->mass) - 0.5 * S * (st_R->rho - st_L->rho);

  flux->momentum[0] = 0.5 * (flux_L->momentum[0] + flux_R->momentum[0]) - 0.5 * S * (st_R->rho * st_R->velx - st_L->rho * st_L->velx);
  flux->momentum[1] = 0.5 * (flux_L->momentum[1] + flux_R->momentum[1]) - 0.5 * S * (st_R->rho * st_R->vely - st_L->rho * st_L->vely);
  flux->momentum[2] = 0.5 * (flux_L->momentum[2] + flux_R->momentum[2]) - 0.5 * S * (st_R->rho * st_R->velz - st_L->rho * st_L->velz);

  flux->energy = 0.5 * (flux_L->energy + flux_R->energy) - 0.5 * S * (st_R->Energy - st_L->Energy);

  flux->B[1] = 0.5 * (flux_L->B[1] + flux_R->B[1]) - 0.5 * S * (st_R->By - st_L->By);
  flux->B[2] = 0.5 * (flux_L->B[2] + flux_R->B[2]) - 0.5 * S * (st_R->Bz - st_L->Bz);
}

#endif
