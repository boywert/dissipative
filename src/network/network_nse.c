/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/network/network_nse.c
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

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_const_cgs.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multiroots.h>

#include "arepoconfig.h"
#include "../allvars.h"

#if defined (NUCLEAR_NETWORK) && defined (NETWORK_NSE)

#include "network_nse.h"

#include "../helm_eos.h"

#include "network.h"
#include "utilities.h"

const static double conv = 1.602177e-12 * 1.0e3 * 6.0221367e23; /* eV2erg * 1.0e3 [keV] * avogadro */

const static double network_nse_parttemp[24] = { 1.0e8, 1.5e8, 2.0e8, 3.0e8, 4.0e8, 5.0e8, 6.0e8, 7.0e8,
  8.0e8, 9.0e8, 1.0e9, 1.5e9, 2.0e9, 2.5e9, 3.0e9, 3.5e9,
  4.0e9, 4.5e9, 5.0e9, 6.0e9, 7.0e9, 8.0e9, 9.0e9, 1.0e10
};

static int network_nse_efunc(const gsl_vector * x, void *params, gsl_vector * f);
static void guessmu(double ye, double *mu_n, double *mu_p);
static void calcprefac(double temp, double rho, struct network_data *nd, struct network_workspace *nw);
static void calcyi(double temp, double rho, double ye, double mu_n, double mu_p, double *yi, double *x, struct network_data *nd, struct network_workspace *nw);
static int get_ye_rate_energy(double time, const double *ye, double *dyedt, void *params);
static int get_ye_rate(double t, const double *ye, double *dyedt, void *params);
static double get_ye_rate_simple(double temp, double rho, double ye, struct network_data *nd, struct network_workspace *nw);
static double network_nse_iterfunc(double temp, void *pp);

#ifdef NETWORK_SCREENING
static double screening_factor(double Gamma_e, double mu_p_coul, int i, struct network_data *nd, struct network_workspace *nw);
static void screening_params(double rho, double ye, double T, double *Gamma_e, double *mu_p_coul);
#endif /* NETWORK_SCREENING */

int network_nse_init(char *speciesfile, char *ratesfile, char *partfile, char *massesfile, char *weakratesfile, struct network_data *nd, struct network_workspace *nw)
{
  network_init_onlyweak(speciesfile, ratesfile, partfile, massesfile, weakratesfile, nd);

  int k;
  for(k = 0; k < NUM_THREADS; k++)
    network_workspace_init(nd, &nw[k]);

  nw->prefac = malloc(nd->nuc_count * sizeof(double));
  return 0;
}

void network_nse_deinit(struct network_data *nd, struct network_workspace *nw)
{
  free(nw->prefac);
  if(nd->initialized)
    network_workspace_deinit(nw);
  if(nw->initialized)
    network_deinit(nd);
}

int network_nse(double temp, double rho, double ye, double *x, struct network_data *nd, struct network_workspace *nw)
{
  int iter, i;
  double mu_n, mu_p, dmu_n, dmu_p;
  double y[2], y2[2];
  double jac[4];                /* order: a11, a12, a21, a22 */
  double det, kt;
  double xsum, nn, ne;
#ifdef NETWORK_SCREENING
  double Gamma_e, mu_p_coul;
#endif /* NETWORK_SCREENING */

  /* preparations */
  calcprefac(temp, rho, nd, nw);
  guessmu(ye, &mu_n, &mu_p);

  /* newton raphson */
  iter = 0;
  while(iter < NSE_MAXITER)
    {
      calcyi(temp, rho, ye, mu_n, mu_p, y, NULL, nd, nw);
      if(dmax(fabs(y[0]), fabs(y[1])) < 1e-12)
        break;

      dmu_n = fabs(mu_n) * NSE_DIFFVAR;
      if(dmu_n == 0.0)
        dmu_n = NSE_DIFFVAR;
      calcyi(temp, rho, ye, mu_n + dmu_n, mu_p, y2, NULL, nd, nw);
      for(i = 0; i < 2; i++)
        jac[i * 2] = (y2[i] - y[i]) / dmu_n;

      dmu_p = fabs(mu_p) * NSE_DIFFVAR;
      if(dmu_p == 0.0)
        dmu_p = NSE_DIFFVAR;
      calcyi(temp, rho, ye, mu_n, mu_p + dmu_p, y2, NULL, nd, nw);
      for(i = 0; i < 2; i++)
        jac[i * 2 + 1] = (y2[i] - y[i]) / dmu_p;

      det = 1.0 / (jac[0] * jac[3] - jac[1] * jac[2]);
      dmu_n = det * (y[0] * jac[3] - y[1] * jac[1]);
      dmu_p = det * (y[1] * jac[0] - y[0] * jac[2]);

      mu_n -= dmu_n;
      mu_p -= dmu_p;
      iter++;
    }

  kt = 1.0 / (GSL_CONST_CGS_BOLTZMANN * temp);

#ifdef NETWORK_SCREENING
  screening_params(rho, ye, temp, &Gamma_e, &mu_p_coul);
#endif /* NETWORK_SCREENING */

  for(i = 0; i < nd->nuc_count; i++)
    {
      x[i] = nw->prefac[i] * exp(kt * (nd->nucdata[i].nz * mu_p + nd->nucdata[i].nn * mu_n + nd->nucdata[i].q));
#ifdef NETWORK_SCREENING
      x[i] *= screening_factor(Gamma_e, mu_p_coul, i, nd, nw);
#endif /* NETWORK_SCREENING */
    }
  if(iter < NSE_MAXITER)
    {
      return 0;
    }
  else
    {
      xsum = 0.0;
      ne = 0.0;
      nn = 0.0;
      for(i = 0; i < nd->nuc_count; i++)
        {
          xsum += x[i];
          ne += nd->nucdata[i].nz * x[i] / nd->nucdata[i].m;
          nn += nd->nucdata[i].na * x[i] / nd->nucdata[i].m;
        }

      myprintf("NSE not converged for T=%13.6e, rho=%13.6e, ye=%13.6e, sum_xi=%13.6e, ye=%13.6e\n", temp, rho, ye, xsum, ne / nn);
      return 1;
    }
}

void network_nse_integrate_ye(double energy, double rho, double x[], double dt, double *dedt, double *temp, struct network_data *nd, struct network_workspace *nw)
{
  double x_old[EOS_NSPECIES];
  int i;

  for(i = 0; i < nd->nuc_count; i++)
    {
      x_old[i] = x[i];
    }

  *temp = network_nse_iter(energy, rho, x, nd, nw);

  double ye = 0.0;
  for(i = 0; i < nd->nuc_count; i++)
    {
      ye += x[i] / nd->nucdata[i].na * nd->nucdata[i].nz;
    }

  double dyedt = get_ye_rate_simple(*temp, rho, ye, nd, nw);

  ye += dyedt * dt;
  /* we ignore the energy change due to the change of ye here */
  network_nse(*temp, rho, ye, x, nd, nw);

  /* calculate change of mass fractions and energy release */
  *dedt = 0.0;
  for(i = 0; i < nd->nuc_count; i++)
    {
      double dxdt = (x[i] - x_old[i]) / dt;
      *dedt -= dxdt / nd->nucdata[i].na * nd->nucdata[i].exm;
    }
  *dedt *= conv;
}

double network_nse_iter(double energy, double rho, double x[], struct network_data *nd, struct network_workspace *nw)
{
  struct nse_params_iter params;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  gsl_function f;
  double ye;
  int i;

  ye = 0.0;
  for(i = 0; i < nd->nuc_count; i++)
    {
      ye += x[i] / nd->nucdata[i].na * nd->nucdata[i].nz;
    }

  params.nd = nd;
  params.nw = nw;
  params.energy = energy;
  params.rho = rho;
  params.ye = ye;
  params.x = x;

  f.function = &network_nse_iterfunc;
  f.params = &params;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc(T);
  gsl_root_fsolver_set(s, &f, 4e9, 2e10);

  int iter = 0;
  int status;
  do
    {
      status = gsl_root_fsolver_iterate(s);

      double temp_lo = gsl_root_fsolver_x_lower(s);
      double temp_hi = gsl_root_fsolver_x_upper(s);

      status = gsl_root_test_interval(temp_lo, temp_hi, 0., 1e-4);

      iter++;
    }
  while(status == GSL_CONTINUE && iter < NSE_MAXITER);

  if(iter < NSE_MAXITER)
    {
      double temp = gsl_root_fsolver_root(s);

      network_nse(temp, rho, ye, x, nd, nw);
      return temp;
    }
  else
    {
      return GSL_FAILURE;
    }
}

static double network_nse_iterfunc(double temp, void *pp)
{
  struct nse_params_iter *params = (struct nse_params_iter *) pp;
  struct network_data *nd = params->nd;
  struct network_workspace *nw = params->nw;
  double energy = params->energy;
  double rho = params->rho;
  double ye = params->ye;
  double *x = params->x;

  double x_new[EOS_NSPECIES];

  network_nse(temp, rho, ye, x_new, nd, nw);

  double delta_energy = 0;
  int i;
  for(i = 0; i < nd->nuc_count; i++)
    {
      double dx = x_new[i] - x[i];
      delta_energy -= dx / nd->nucdata[i].na * nd->nucdata[i].exm;
    }
  delta_energy *= conv;

  double energy_new = energy + delta_energy;

  double temp_new = temp;
  struct eos_result res;
  eos_calc_egiven(rho, x_new, energy_new, &temp_new, &res);

  return (temp - temp_new) / temp;
}

int network_nse_egiven(double energy, double rho, double ye, double *xnuc, double *tempguess, struct network_data *nd, struct network_workspace *nw)
{
  struct network_nse_eparams params;
  double mu_n, mu_p, temp;
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;
  gsl_multiroot_function f;
  gsl_vector *x;
  int status, iter, i;
  double xsum, nn, ne, kt;
#ifdef NETWORK_SCREENING
  double Gamma_e, mu_p_coul;
#endif /* NETWORK_SCREENING */

  /* no NSE */
  if(rho < 1e7)
    return -1;

  guessmu(ye, &mu_n, &mu_p);

  params.nd = nd;
  params.nw = nw;
  params.rho = rho;
  params.energy = energy;
  params.ye = ye;

  if(tempguess == 0 || *tempguess == 0)
    temp = 7e9;
  else
    temp = *tempguess;

  f.f = &network_nse_efunc;
  f.n = 3;
  f.params = &params;

  x = gsl_vector_alloc(3);
  gsl_vector_set(x, 0, mu_n);
  gsl_vector_set(x, 1, mu_p);
  gsl_vector_set(x, 2, temp / 1e14);

  //T = gsl_multiroot_fsolver_dnewton;
  T = gsl_multiroot_fsolver_hybrids;
  s = gsl_multiroot_fsolver_alloc(T, 3);
  gsl_multiroot_fsolver_set(s, &f, x);

  iter = 0;

  do
    {
      /*
         printf( "iter=%d\n x=(%25.16e|%25.16e|%25.16e)\n dx=(%25.16e|%25.16e|%25.16e)\n f=(%25.16e|%25.16e|%25.16e)\n", iter,
         gsl_vector_get( s->x, 0 ), gsl_vector_get( s->x, 1 ), gsl_vector_get( s->x, 2 ),
         gsl_vector_get( s->dx, 0 ), gsl_vector_get( s->dx, 1 ), gsl_vector_get( s->dx, 2 ),
         gsl_vector_get( s->f, 0 ), gsl_vector_get( s->f, 1 ), gsl_vector_get( s->f, 2 ) );
       */
      status = gsl_multiroot_fsolver_iterate(s);

      if(status)
        break;

      status = gsl_multiroot_test_residual(s->f, 1e-10);

      iter++;
    }
  while(status == GSL_CONTINUE && iter < NSE_MAXITER);
  mu_n = gsl_vector_get(s->x, 0);
  mu_p = gsl_vector_get(s->x, 1);
  temp = gsl_vector_get(s->x, 2) * 1e14;

  gsl_multiroot_fsolver_free(s);
  gsl_vector_free(x);

  if(status == GSL_ENOMEM || status == GSL_EBADFUNC || status == GSL_ENOPROG)
    {
      myprintf("NSE solver failed for T=%13.6e, rho=%13.6e, ye=%13.6e, energy=%25.13e\n", temp, rho, ye, energy);
      return status;
    }

  kt = 1.0 / (GSL_CONST_CGS_BOLTZMANN * temp);

#ifdef NETWORK_SCREENING
  screening_params(rho, ye, temp, &Gamma_e, &mu_p_coul);
#endif /* NETWORK_SCREENING */

  for(i = 0; i < nd->nuc_count; i++)
    {
      xnuc[i] = nw->prefac[i] * exp(kt * (nd->nucdata[i].nz * mu_p + nd->nucdata[i].nn * mu_n + nd->nucdata[i].q));
#ifdef NETWORK_SCREENING
      xnuc[i] *= screening_factor(Gamma_e, mu_p_coul, i, nd, nw);
#endif /* NETWORK_SCREENING */
    }

  if(tempguess)
    *tempguess = temp;

  if(iter < NSE_MAXITER)
    {
      return 0;
    }
  else
    {
      xsum = 0.0;
      ne = 0.0;
      nn = 0.0;
      for(i = 0; i < nd->nuc_count; i++)
        {
          xsum += xnuc[i];
          ne += nd->nucdata[i].nz * xnuc[i] / nd->nucdata[i].m;
          nn += nd->nucdata[i].na * xnuc[i] / nd->nucdata[i].m;
        }

      myprintf("NSE not converged for T=%13.6e, rho=%13.6e, ye=%13.6e, sum_xi=%13.6e, ye=%13.6e, energy=%25.13e\n", temp, rho, ye, xsum, ne / nn, energy);

      return GSL_FAILURE;
    }
}

static int network_nse_efunc(const gsl_vector * x, void *params, gsl_vector * f)
{
  struct network_data *nd = ((struct network_nse_eparams *) params)->nd;
  struct network_workspace *nw = ((struct network_nse_eparams *) params)->nw;
  double rho = ((struct network_nse_eparams *) params)->rho;
  double energy = ((struct network_nse_eparams *) params)->energy;
  double ye = ((struct network_nse_eparams *) params)->ye;
  double abar, zbar, xsum;
  struct helm_eos_cache cache;

  double *xnuc = malloc(nd->nuc_count * sizeof(double));
  struct eos_result res;

  double mu_n = gsl_vector_get(x, 0);
  double mu_p = gsl_vector_get(x, 1);
  double temp = gsl_vector_get(x, 2) * 1e14;

  double yi[2];

  temp = dmin(dmax(temp, 1e9), 2e10);
  mu_n = dmin(dmax(mu_n, -1e-4), 1e-4);
  mu_p = dmin(dmax(mu_p, -1e-4), 1e-4);

  calcprefac(temp, rho, nd, nw);

  calcyi(temp, rho, ye, mu_n, mu_p, yi, xnuc, nd, nw);
  gsl_vector_set(f, 0, yi[0]);
  gsl_vector_set(f, 1, yi[1]);

  xsum = 0;
  abar = 0;
  zbar = 0;
  int i;
  for(i = 0; i < nd->nuc_count; i++)
    {
      xsum += xnuc[i];
      abar += xnuc[i] / nd->nucdata[i].na;
      zbar += xnuc[i] / nd->nucdata[i].na * nd->nucdata[i].nz;
    }

  zbar = zbar / abar;
  abar = xsum / abar;

  helm_eos_update_cache(rho, abar, zbar, &cache);
  if(eos_calc_tgiven_azbar(rho, &cache, temp, &res, 1))
    {
      gsl_vector_set(f, 2, energy);
      return GSL_FAILURE;
    }

  gsl_vector_set(f, 2, (res.e.v - energy) / energy);

  free(xnuc);

  return GSL_SUCCESS;
}

#ifdef NETWORK_SCREENING
static double screening_factor(double Gamma_e, double mu_p_coul, int i, struct network_data *nd, struct network_workspace *nw)
{
  /* screening as mentioned in Seitenzahl et al. Atomic and Nuclear Data Tables 95 (2009) 96-114 eq. (14) */
  const double A1 = -0.9052, A2 = 0.6322, A3 = -sqrt(3.0) / 2.0 - A1 / sqrt(A2);
  double Gamma_i;

  if(nd->nucdata[i].nz != 0.0)
    {
      Gamma_i = Gamma_e * pow(nd->nucdata[i].nz, 5.0 / 3.0);
      return exp(-(A1 * (sqrt(Gamma_i * (A2 + Gamma_i)) - A2 * log(sqrt(Gamma_i / A2) + sqrt(1.0 + Gamma_i / A2))) + 2.0 * A3 * (sqrt(Gamma_i) - atan(sqrt(Gamma_i)))) + nd->nucdata[i].nz * mu_p_coul);
    }
  else
    return 1.0;
}

static void screening_params(double rho, double ye, double T, double *Gamma_e, double *mu_p_coul)
{
  const double A1 = -0.9052, A2 = 0.6322, A3 = -sqrt(3.0) / 2.0 - A1 / sqrt(A2);

  *Gamma_e = gsl_pow_2(ELECTRON_CHARGE_ESU) / GSL_CONST_CGS_BOLTZMANN / T * pow(4.0 / 3.0 * M_PI * rho * ye / GSL_CONST_CGS_UNIFIED_ATOMIC_MASS, 1.0 / 3.0);
  *mu_p_coul = (A1 * (sqrt(*Gamma_e * (A2 + *Gamma_e)) - A2 * log(sqrt(*Gamma_e / A2) + sqrt(1.0 + *Gamma_e / A2))) + 2.0 * A3 * (sqrt(*Gamma_e) - atan(sqrt(*Gamma_e))));
}
#endif /* NETWORK_SCREENING */

static void guessmu(double ye, double *mu_n, double *mu_p)
{
  if(ye > 0.55)
    {
      *mu_n = -2.7e-5;
      *mu_p = -1.1e-6;
    }
  else if(ye > 0.53)
    {
      *mu_n = -0.22e-4;
      *mu_p = -0.06e-4;
    }
  else if(ye > 0.46)
    {
      *mu_n = -0.19e-4;
      *mu_p = -0.08e-4;
    }
  else if(ye > 0.44)
    {
      *mu_n = -0.145e-4;
      *mu_p = -0.135e-4;
    }
  else if(ye > 0.42)
    {
      *mu_n = -0.11e-4;
      *mu_p = -0.18e-4;
    }
  else
    {
      *mu_n = -0.08e-4;
      *mu_p = -0.22e-4;
    }
}

static void calcprefac(double temp, double rho, struct network_data *nd, struct network_workspace *nw)
{
  int i;
  network_part(temp, nd, nw);

  for(i = 0; i < nd->nuc_count; i++)
    {
      nw->prefac[i] =
        nd->nucdata[i].m / rho * (2.0 * nd->nucdata[i].spin +
                                  1.0) * nw->gg[i].v * pow(2.0 * M_PI * nd->nucdata[i].m * GSL_CONST_CGS_BOLTZMANN * temp / (GSL_CONST_CGS_PLANCKS_CONSTANT_H * GSL_CONST_CGS_PLANCKS_CONSTANT_H), 1.5);
    }
}

static void calcyi(double temp, double rho, double ye, double mu_n, double mu_p, double *yi, double *x, struct network_data *nd, struct network_workspace *nw)
{
  int i;
  double kt, xi;
#ifdef NETWORK_SCREENING
  double Gamma_e, mu_p_coul;
#endif /* NETWORK_SCREENING */

  yi[0] = -1.0;
  yi[1] = 0.0;

  kt = 1.0 / (GSL_CONST_CGS_BOLTZMANN * temp);
#ifdef NETWORK_SCREENING
  screening_params(rho, ye, temp, &Gamma_e, &mu_p_coul);
#endif /* NETWORK_SCREENING */
  for(i = 0; i < nd->nuc_count; i++)
    {
      xi = nw->prefac[i] * exp(kt * (nd->nucdata[i].nz * mu_p + nd->nucdata[i].nn * mu_n + nd->nucdata[i].q));
#ifdef NETWORK_SCREENING
      xi *= screening_factor(Gamma_e, mu_p_coul, i, nd, nw);
#endif /* NETWORK_SCREENING */
      yi[0] += xi;
      yi[1] += GSL_CONST_CGS_UNIFIED_ATOMIC_MASS / nd->nucdata[i].m * ((ye - 1) * nd->nucdata[i].nz + ye * nd->nucdata[i].nn) * xi;

      if(x)
        x[i] = xi;
    }
}


void network_nse_integrate_traj(struct network_data *nd, struct network_workspace *nw, struct network_solver_trajectory *traj)
{
  double rho, energy, temp, ye, dt;
  double lasttime, lasttemp, lastye;
  // variables for the GSL ODE solver
  gsl_odeiv_system ode;
  gsl_odeiv_evolve *evolve;
  gsl_odeiv_control *control;
  gsl_odeiv_step *step;
  struct nse_params_ene params;
  int status;
  double maxtime;

  evolve = gsl_odeiv_evolve_alloc(1);
  control = gsl_odeiv_control_standard_new(1e-3, 1e-3, 0.5, 0.5);
  // step = gsl_odeiv_step_alloc(gsl_odeiv_step_rk8pd, 1);
  step = gsl_odeiv_step_alloc(gsl_odeiv_step_rk4, 1);
  ode.function = &get_ye_rate_energy;
  ode.jacobian = NULL;
  ode.dimension = 1;

  temp = 0;

  params.nd = nd;
  params.nw = nw;
  params.traj = traj;
  params.temp = temp;
  ode.params = &params;

  ye = 0.0;
  int i;
  for(i = 0; i < nd->nuc_count; i++)
    {
      ye += traj->x[i] / nd->nucdata[i].na * nd->nucdata[i].nz;
    }

  dt = 1e-2;

  maxtime = traj->maxtime;
  while(traj->time < maxtime)
    {
      lasttime = traj->time;
      lasttemp = temp;
      lastye = ye;

      dt = dmin(dt, 0.1);       // we do not accept a timestep larger than 0.1s
      status = gsl_odeiv_evolve_apply(evolve, control, step, &ode, &traj->time, maxtime, &dt, &ye);

      if(status == 666)
        {
          dt *= 0.1;
          myprintf("step failed, reducing timestep to %g\n", dt);
          continue;
        }

      rho = params.rho;
      temp = params.temp;

      if(traj->mintemp && temp < traj->mintemp)
        {
          if(temp < 0.9 * traj->mintemp)
            {
              maxtime = lasttime + (traj->mintemp - lasttemp) / (temp - lasttemp) * (traj->time - lasttime);

              traj->time = lasttime;
              params.temp = lasttemp;
              ye = lastye;
              dt *= 0.5;
            }
          else
            {
              break;
            }
        }
    }

  // get the new abundances for the new value of ye
  network_solver_interpolate_trajectory(traj, traj->time, &rho, &energy);
  network_nse_egiven(energy, rho, ye, traj->x, NULL, nd, nw);

  gsl_odeiv_evolve_free(evolve);
  gsl_odeiv_control_free(control);
  gsl_odeiv_step_free(step);
}

static int get_ye_rate_energy(double time, const double *ye, double *dyedt, void *params)
{
  struct network_data *nd;
  struct network_workspace *nw;
  struct network_solver_trajectory *traj;
  double *y, *dydt;
  double energy, rho, temp, yetemp;

  nd = ((struct nse_params_ene *) params)->nd;
  nw = ((struct nse_params_ene *) params)->nw;
  traj = ((struct nse_params_ene *) params)->traj;
  temp = ((struct nse_params_ene *) params)->temp;

  if(*ye < 0 || *ye > 1.0)
    return 666;

  yetemp = dmin(dmax(0, *ye), 1.0);

  network_solver_interpolate_trajectory(traj, time, &rho, &energy);
  if(network_nse_egiven(energy, rho, yetemp, traj->x, &temp, nd, nw))
    return 666;

  if(isnan(temp))
    return 666;

  y = malloc(nd->n_matrix * sizeof(double));
  dydt = malloc(nd->n_matrix * sizeof(double));

  int i;
  for(i = 0; i < nd->nuc_count; i++)
    {
      y[i] = traj->x[i] / nd->nucdata[i].na;
    }
#ifdef NETWORK_VARIABLE
  y[nd->iTemp] = temp;
#endif

  network_getrhs(rho, temp, y, 0, nd, nw, dydt, nw->nsd->rhs_deriv);

  *dyedt = 0.0;
  for(i = 0; i < nd->nuc_count; i++)
    {
      *dyedt += nd->nucdata[i].nz * dydt[i];
    }

  free(y);
  free(dydt);

  ((struct nse_params_ene *) params)->rho = rho;
  ((struct nse_params_ene *) params)->temp = temp;

  return GSL_SUCCESS;
}

void network_nse_integrate_traj_temp(struct network_data *nd, struct network_workspace *nw, struct network_solver_trajectory *traj)
{
  double rho, energy, temp, ye, dt;
  // variables for the GSL ODE solver
  gsl_odeiv_system ode;
  gsl_odeiv_evolve *evolve;
  gsl_odeiv_control *control;
  gsl_odeiv_step *step;
  struct nse_params params;
  struct eos_result res;

  evolve = gsl_odeiv_evolve_alloc(1);
  control = gsl_odeiv_control_standard_new(1e-6, 1e-6, 0.5, 0.5);
  step = gsl_odeiv_step_alloc(gsl_odeiv_step_rk8pd, 1);
  ode.function = &get_ye_rate;
  ode.jacobian = NULL;
  ode.dimension = 1;
  params.nd = nd;
  params.nw = nw;
  params.temp = 0;
  params.rho = 0;
  ode.params = &params;

  ye = 0.0;
  int i;
  for(i = 0; i < nd->nuc_count; i++)
    {
      ye += traj->x[i] / nd->nucdata[i].na * nd->nucdata[i].nz;
    }

  temp = 7e9;
  network_solver_interpolate_trajectory(traj, traj->time, &rho, &energy);
  network_nse_egiven(energy, rho, ye, traj->x, &temp, nd, nw);

  dt = 1e-2;
  while(traj->time < traj->maxtime)
    {
      network_solver_interpolate_trajectory(traj, traj->time, &rho, &energy);
      eos_calc_egiven(rho, traj->x, energy, &temp, &res);

      params.temp = temp;
      params.rho = rho;

      if((traj->mintemp && temp < traj->mintemp) || (traj->maxtemp && temp > traj->maxtemp))
        {
          break;
        }

      gsl_odeiv_evolve_apply(evolve, control, step, &ode, &traj->time, traj->maxtime, &dt, &ye);
    }

  // get the new abundances for the new value of ye
  network_solver_interpolate_trajectory(traj, traj->time, &rho, &energy);
  network_nse_egiven(energy, rho, ye, traj->x, NULL, nd, nw);

  gsl_odeiv_evolve_free(evolve);
  gsl_odeiv_control_free(control);
  gsl_odeiv_step_free(step);
}

void network_nse_integrate(double temp, double rho, double x[], double dt, double *dedt, struct network_data *nd, struct network_workspace *nw, int verbose)
{
  double ye;
  double *new_x, dx;
  double t, sub_dt;
  int i;
  /* variables for the GSL ODE solver */
  gsl_odeiv_system ode;
  gsl_odeiv_evolve *evolve;
  gsl_odeiv_control *control;
  gsl_odeiv_step *step;
  struct nse_params params;

  evolve = gsl_odeiv_evolve_alloc(1);
  control = gsl_odeiv_control_standard_new(1e-6, 1e-6, 0.5, 0.5);
  step = gsl_odeiv_step_alloc(gsl_odeiv_step_rk8pd, 1);
  ode.function = &get_ye_rate;
  ode.jacobian = NULL;
  ode.dimension = 1;
  params.nd = nd;
  params.nw = nw;
  params.temp = temp;
  params.rho = rho;
  ode.params = &params;

  new_x = malloc(nd->nuc_count * sizeof(double));

  ye = 0.0;
  for(i = 0; i < nd->nuc_count; i++)
    {
      ye += x[i] / nd->nucdata[i].na * nd->nucdata[i].nz;
    }

  if(verbose)
    myprintf("Integrating NSE for temp=%g, rho=%g, ye=%g, xHe=%g, xC12=%g, xO16=%g, xNi56=%g\n", temp, rho, ye, x[2], x[3], x[5], x[62]);

  t = 0.0;
  sub_dt = dt;
  while(t < dt)
    {
      if(verbose)
        myprintf("t = %g: time step in NSE for the weak rates %g of %g\n", t, sub_dt, dt);
      if(gsl_odeiv_evolve_apply(evolve, control, step, &ode, &t, dt, &sub_dt, &ye) != GSL_SUCCESS)
        {
          if(verbose)
            myprintf("step was too large: reducing by a factor of 2\n");
          sub_dt *= 0.5;
        }
    }

  /* get the new abundances for the new value of ye */
  network_nse(temp, rho, ye, new_x, nd, nw);

  if(verbose)
    myprintf("t = %g: finished integration, new ye=%g, xHe=%g, xC12=%g, xO16=%g, xNi56=%g\n", t, ye, new_x[2], new_x[3], new_x[5], new_x[62]);

  /* calculate change of mass fractions and energy release */
  *dedt = 0.0;
  for(i = 0; i < nd->nuc_count; i++)
    {
      dx = (new_x[i] - x[i]) / dt;
      *dedt -= dx / nd->nucdata[i].na * nd->nucdata[i].exm;
      x[i] = new_x[i];
    }
  *dedt *= conv;

  if(verbose)
    myprintf("Energy released: %g\n", *dedt);

  gsl_odeiv_evolve_free(evolve);
  gsl_odeiv_control_free(control);
  gsl_odeiv_step_free(step);
  free(new_x);
}

static double get_ye_rate_simple(double temp, double rho, double ye, struct network_data *nd, struct network_workspace *nw)
{
  double *y, *dydt;
  int i;

  y = malloc(nd->n_matrix * sizeof(double));
  dydt = malloc(nd->n_matrix * sizeof(double));

  network_nse(temp, rho, ye, y, nd, nw);        /* the returned y is a mass fraction, we change it to a number fraction immediately */

  for(i = 0; i < nd->nuc_count; i++)
    {
      y[i] /= nd->nucdata[i].na;
    }
#ifdef NETWORK_VARIABLE
  y[nd->iTemp] = temp;
#endif

  network_getrhs(rho, temp, y, 0, nd, nw, dydt, nw->nsd->rhs_deriv);

  double dyedt = 0.0;
  for(i = 0; i < nd->nuc_count; i++)
    {
      dyedt += nd->nucdata[i].nz * dydt[i];
    }

  free(y);
  free(dydt);

  return dyedt;
}

static int get_ye_rate(double t, const double *ye, double *dyedt, void *params)
{
  struct network_data *nd;
  struct network_workspace *nw;
  double temp, rho;

  if(*ye < 0.0 || *ye > 1.0)
    {
      /* probably the time step was too large */
      return GSL_SUCCESS + 1;   /* return a value different from GSL_SUCCESS */
    }

  nd = ((struct nse_params *) params)->nd;
  nw = ((struct nse_params *) params)->nw;
  temp = ((struct nse_params *) params)->temp;
  rho = ((struct nse_params *) params)->rho;

  /* we do not need the parameter t */
  (void) t;

  *dyedt = get_ye_rate_simple(temp, rho, *ye, nd, nw);

  return GSL_SUCCESS;
}
#endif /* NUCLEAR_NETWORK */
