/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/riemann_gamma.c
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

#if defined(RIEMANN_GAMMA)

#if defined(RIEMANN_HLLC) || defined(RIEMANN_HLLD) || defined(RIEMANN_ROSUNOV) || defined(RIEMANN_HLL)
#error option RIEMANN_GAMMA is incompatible with options RIEMANN_HLLC, RIEMANN_HLLD, RIEMANN_ROSUNOV \
       and RIEMANN_HLL. Only one Riemann solver can be chosen among the above options. If none of    \
       them is selected, the exact Riemann solver will be used.
#endif

#define TOL 1.0e-8
#define SMALLP 1.0e-20

double godunov_flux_3d_gamma(struct state *st_L, struct state *st_R, struct state_face *st_face)
{
#ifdef ISOTHERM_EQS
#error "ISOTHERM_EQS does not work with RIEMANN_GAMMA"
#endif
  if(st_L->rho > 0 && st_R->rho > 0)
    {
      double Press_star, Vel_star, W_L, W_R;

      st_L->csnd = sqrt(st_L->gammaC * st_L->press / st_L->rho);
      st_R->csnd = sqrt(st_R->gammaC * st_R->press / st_R->rho);

      riemann_gamma(st_L, st_R, &Press_star, &Vel_star, &W_L, &W_R);
      sample_solution_3d_gamma(0.0,     /* S=x/t */
                               st_L, st_R, Press_star, Vel_star, W_L, W_R, st_face);
      return Press_star;
    }
  else
    {
      terminate("density is zero\n");
      return 0;
    }
}

void get_mach_numbers_gamma(struct state *st_L, struct state *st_R, double Press, double *Mach_L, double *Mach_R)
{
#if (defined ENTROPY_MACH_THRESHOLD) || (defined GODUNOV_STATS)
  if(Press <= st_L->press)      /* left fan */
    {
      st_L->mach = 0;
    }
  else                          /* left shock */
    {
      double pml = Press / st_L->press;

      double GammaG2 = (st_L->gammaC + 1.0) / (2.0 * st_L->gammaC);
      double GammaG1 = (st_L->gammaC - 1.0) / (2.0 * st_L->gammaC);

      st_L->mach = sqrt(GammaG2 * pml + GammaG1);
    }

  if(Press > st_R->press)       /* right shock */
    {
      double pmr = Press / st_R->press;

      double GammaG2 = (st_R->gammaC + 1.0) / (2.0 * st_R->gammaC);
      double GammaG1 = (st_R->gammaC - 1.0) / (2.0 * st_R->gammaC);

      st_R->mach = sqrt(GammaG2 * pmr + GammaG1);
    }
  else
    {
      st_R->mach = 0;
    }
#endif
}


void sample_solution_3d_gamma(double S, struct state *st_L, struct state *st_R, double Press_star, double Vel_star, double W_L, double W_R, struct state_face *st_face)
{
  double Rho, Vel, Press, GammaE, GammaC, W, Csnd;
  double Rho_star, Csnd_star, Sign, We, We_star, gamfac, GammaE_star;
  double fac1, fac2;

  if(S <= Vel_star)             /* sample point is left of contact */
    {
      Rho = st_L->rho;
      Vel = st_L->velx;
      Press = st_L->press;
      GammaE = st_L->gammaE;
      GammaC = st_L->gammaC;

      st_face->vely = st_L->vely;
      st_face->velz = st_L->velz;
#ifdef MAXSCALARS
      st_face->scalars = st_L->scalars;
#endif

      Csnd = st_L->csnd;
      W = W_L;
      Sign = 1.0;
    }
  else
    {
      Rho = st_R->rho;
      Vel = st_R->velx;
      Press = st_R->press;
      GammaE = st_R->gammaE;
      GammaC = st_R->gammaC;

      st_face->vely = st_R->vely;
      st_face->velz = st_R->velz;
#ifdef MAXSCALARS
      st_face->scalars = st_R->scalars;
#endif

      Csnd = st_R->csnd;
      W = W_R;
      Sign = -1.0;
    }

  Rho_star = 1.0 / (1.0 / Rho - (Press_star - Press) / (W * W));
  Csnd_star = sqrt(GammaC * Press_star / Rho_star);

  if(Press_star >= Press)
    {
      We = W / Rho - Sign * Vel;
      We_star = We;
    }
  else
    {
      We = Csnd - Sign * Vel;
      We_star = Csnd_star - Sign * Vel_star;
    }

  gamfac = (1.0 - GammaE / GammaC) * (GammaE - 1.0);
  GammaE_star = GammaE + 2.0 * gamfac * (Press_star - Press) / (Press_star + Press);

  Rho = dmax(Rho, 1e-5);
  Rho_star = dmax(Rho_star, 1e-5);

  if(We < 0.)
    {
      st_face->rho = Rho;
      st_face->velx = Vel;
      st_face->press = Press;
      st_face->gammaE = GammaE;
    }
  else if(We_star >= 0.)
    {
      st_face->rho = Rho_star;
      st_face->velx = Vel_star;
      st_face->press = Press_star;
      st_face->gammaE = GammaE_star;
    }
  else
    {
      fac1 = 0.5 * (1. + (We + We_star) / dmax(dmax(We - We_star, We + We_star), 1e-5));
      fac2 = 1. - fac1;

      st_face->rho = fac1 * Rho_star + fac2 * Rho;
      st_face->velx = fac1 * Vel_star + fac2 * Vel;
      st_face->press = fac1 * Press_star + fac2 * Press;
      st_face->gammaE = fac1 * GammaE_star + fac2 * GammaE;
    }

  if(!gsl_finite(st_face->gammaE))
    printf("GammaE=%g, GammaE_star=%g, fac1=%g, fac2=%g, We=%g, We_star=%g, W_L=%g, W_R=%g, Rho=%g\n", GammaE, GammaE_star, fac1, fac2, We, We_star, W_L, W_R, Rho);
}


int riemann_gamma(struct state *st_L, struct state *st_R, double *Press, double *Vel, double *W_L, double *W_R)
{
  double p, pold, pnew, udiff, udiffold;

  double GammaEav = st_L->gammaE + st_R->gammaE;
  double GammaCav = st_L->gammaC + st_R->gammaC;
  double GammaFac = (1. - GammaEav / GammaCav) * (GammaEav - 2.);

  double GammaEmin = dmin(st_L->gammaE, st_R->gammaE);
  double GammaEmax = dmax(st_L->gammaE, st_R->gammaE);

  double critVel = 2.0 / (st_L->gammaC - 1.0) * st_L->csnd + 2.0 / (st_R->gammaC - 1.0) * st_R->csnd - (st_R->velx - st_L->velx);
  if(critVel < 0)
    {
      /* results will be crap (vacuum), but go on... */
    }

  pold = dmax(st_L->press + (st_R->press - st_L->press - st_R->csnd * (st_R->velx - st_L->velx)) * (st_L->csnd / (st_L->csnd + st_R->csnd)), SMALLP);

  calcW(pold, GammaFac, st_L, GammaEmin, GammaEmax, W_L);
  calcW(pold, GammaFac, st_R, GammaEmin, GammaEmax, W_R);

  udiff = (st_L->velx - (pold - st_L->press) / (*W_L)) - (st_R->velx + (pold - st_R->press) / (*W_R));
  p = dmax(st_L->press + (st_R->press - st_L->press - *W_R * (st_R->velx - st_L->velx)) * (*W_L / (*W_L + *W_R)), SMALLP);
  pnew = p;

  int iter = 0;

  do                            /* newton-raphson scheme */
    {
      calcW(p, GammaFac, st_L, GammaEmin, GammaEmax, W_L);
      calcW(p, GammaFac, st_R, GammaEmin, GammaEmax, W_R);

      udiffold = udiff;
      udiff = (st_L->velx - (p - st_L->press) / (*W_L)) - (st_R->velx + (p - st_R->press) / (*W_R));

      if(udiff - udiffold == 0.)
        break;

      pnew = dmax(p - udiff * (p - pold) / (udiff - udiffold), SMALLP);
      pold = p;
      p = pnew;

      iter++;
    }
  while(2. * fabs((p - pold) / (p + pold)) > TOL && iter < MAXITER);

  /* prepare output values */
  *Press = p;
  *Vel = 0.5 * ((st_L->velx - (p - st_L->press) / (*W_L)) + (st_R->velx + (p - st_R->press) / (*W_R)));

  if(!gsl_finite(*W_L) || !gsl_finite(*W_R))
    printf
      ("W_L=%g, W_R=%g, Press=%g, pold=%g, Vel=%g, critVel=%g, v_L=%g, v_R=%g, P_L=%g, P_R=%g, csnd=%g|%g, GammaC=%g|%g, Gamma=%g|%g|%g\n",
       *W_L, *W_R, *Press, pold, *Vel, critVel, st_L->velx, st_R->velx, st_L->press, st_R->press, st_L->csnd, st_R->csnd, st_L->gammaC, st_R->gammaC, GammaFac, GammaEmin, GammaEmax);

  if(iter >= MAXITER)
    terminate("failed convergence in Riemann solver");

  return 1;
}

void calcW(double Pstar, double GammaFac, struct state *st, double GammaEmin, double GammaEmax, double *W)
{
  double pdiff = Pstar - st->press;
  double psum = Pstar + st->press;

  /* GammaFac = ( 1. - GammaEav / GammaCav ) * ( GammaEav - 1. ) */
  double gammaE_star = st->gammaE + GammaFac * pdiff / psum;

  gammaE_star = dmax(GammaEmin, dmin(GammaEmax, gammaE_star));

  double sqgame = sqrt(0.5 * (st->gammaE - 1.0) / st->gammaE);

  if(fabs(pdiff / psum) < 1e-20)
    {
      *W = st->csnd;
    }
  else
    {
      *W = sqrt(pdiff * (Pstar + 0.5 * (gammaE_star - 1.0) * psum) * st->rho / (Pstar - (gammaE_star - 1.0) / (st->gammaE - 1.0) * st->press));
    }
  *W = dmax(*W, sqgame * st->csnd);
}

#endif
