/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/tgchem/tgchem_step.c
 * \date        01/2013
 * \author      Thomas Greif
 * \brief       Primordial chemistry and cooling network
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 */

#ifdef TGCHEM_TEST
#include "tgchem_test.h"
#else
#include "../allvars.h"
#include "../proto.h"
#endif

void tgchem_bisection(int mode, double *val, double min_frac, double max_frac, double *aux);
double tgchem_eq_func(int mode, double *val, double init_val, double *aux, double init_aux);
double tgchem_quadratic(int branch, double a, double b, double c);


void tgchem_step(double dt)
{
  int i, flag_break, flag_repeat;
  double dt_step, dt_step_sec, dt_sum, dt_sum_sec, tout;
  double t0, t1;
  double abh2, abhii, energy;

  dt_step = TGCHEM_TOL_STEPSIZE * dt;

  dt_sum = 0.;

  flag_break = 0;

  while(!flag_break)
    {
      if(dt_sum + dt_step >= dt)
        {
          dt_step = dmax(dt - dt_sum, 0);

          flag_break = 1;
        }

      t0 = second();

      if(TGCD.EqHIIFlag)
        {
          abhii = NV_Ith_S(TGCD.Species, 2);

          tgchem_bisection(1, &abhii, 1e-2, 1e2, 0);

          NV_Ith_S(TGCD.Species, 2) = abhii;

          if(tgchem_check_species())
            {
              printf("Something wrong with species! Mode 1: %d %d %d %g %g %g\n", TGCD.Task, TGCD.Index, TGCD.ID, NV_Ith_S(TGCD.Species, 0), NV_Ith_S(TGCD.Species, 1), NV_Ith_S(TGCD.Species, 2));

              terminate("");
            }
        }

      if(TGCD.EqH2Flag)
        {
          abh2 = NV_Ith_S(TGCD.Species, 1);
          energy = NV_Ith_S(TGCD.Species, TGCHEM_NUM_ABUNDANCES);

          tgchem_bisection(0, &energy, 1e-1, 1e1, &abh2);

          NV_Ith_S(TGCD.Species, 1) = abh2;
          NV_Ith_S(TGCD.Species, TGCHEM_NUM_ABUNDANCES) = energy;

          if(tgchem_check_species())
            {
              printf("Something wrong with species! Mode 2: %d %d %d %g %g %g\n", TGCD.Task, TGCD.Index, TGCD.ID, NV_Ith_S(TGCD.Species, 0), NV_Ith_S(TGCD.Species, 1), NV_Ith_S(TGCD.Species, 2));

              terminate("");
            }
        }

      t1 = second();

      TGCD.DtEq += timediff(t0, t1);

      if(TGCD.EqHIIFlag)
        break;

      dt_step_sec = dt_step;
      dt_sum_sec = 0.;

      while(dt_sum_sec < dt_step)
        {
          for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
            NV_Ith_S(TGCD.SpeciesSave, i) = NV_Ith_S(TGCD.Species, i);

          flag_repeat = CVodeReInit(TGCD.CVODEMem, 0, TGCD.Species);

          if(!flag_repeat)
            {
              t0 = second();

              flag_repeat = CVode(TGCD.CVODEMem, dt_step_sec, TGCD.Species, &tout, CV_NORMAL);

              t1 = second();

              TGCD.DtNEq += timediff(t0, t1);
            }

          if(!flag_repeat)
            flag_repeat = tgchem_check_species();

          if(flag_repeat)
            {
              for(i = 0; i < TGCHEM_NUM_SPECIES; i++)
                NV_Ith_S(TGCD.Species, i) = NV_Ith_S(TGCD.SpeciesSave, i);

              dt_step_sec /= 2.;
            }
          else
            {
              dt_sum_sec += dt_step_sec;
              dt_step_sec = dt_step - dt_sum_sec;

              TGCD.NumSubSteps++;
            }
        }

      dt_sum += dt_step;
    }
}


void tgchem_bisection(int mode, double *val, double min_frac, double max_frac, double *aux)
{
  double init_val, init_aux, min_val, max_val, old_val, func;

  init_val = *val;
  init_aux = *aux;

  min_val = min_frac * (*val);
  max_val = max_frac * (*val);

  if(mode)
    max_val = dmin(max_val, TGCD.AbMax[2]);

  while(*val != min_val && *val != max_val)
    {
      old_val = *val;

      func = tgchem_eq_func(mode, val, init_val, aux, init_aux);

      if(func <= 0)
        max_val = *val;
      else
        min_val = *val;

      *val = (min_val + max_val) / 2.;

      if(dabs(*val - old_val) < TGCHEM_TOL_EQ * old_val)
        break;
    }
}


double tgchem_eq_func(int mode, double *val, double init_val, double *aux, double init_aux)
{
  double abh, K0, K1, quad, func;

  if(!mode)
    {
      abh = dmax(1. - NV_Ith_S(TGCD.Species, 2), 0.);

      NV_Ith_S(TGCD.Species, TGCHEM_NUM_ABUNDANCES) = *val;
    }
  else
    {
      abh = 1.;

      NV_Ith_S(TGCD.Species, 2) = *val;
    }

  tgchem_compute_vars(0, TGCD.Species);

  TGCHEM_TEMP();

  if(!mode)
    {
      TGCHEM_EQRATE(0);

      K0 = TGCD.EqRate[0] / TGCD.NH;

      quad = tgchem_quadratic(1., 4., -(4. * abh + K0), abh * abh);

      *aux = dmin(dmax(quad, 0.), TGCD.AbMax[1]);

      func = init_val + TGCHEM_CHI_H2 * (*aux - init_aux) * TGCD.NH - (*val);
    }
  else
    {
      TGCHEM_EQRATE(1);

      K1 = TGCD.EqRate[1] / TGCD.NH;

      quad = tgchem_quadratic(0., 1., K1, -K1 * abh);

      func = dmin(dmax(quad, 0.), TGCD.AbMax[2]) - (*val);
    }

  return func;
}


double tgchem_quadratic(int branch, double a, double b, double c)
{
  double b2, d, quad;

  b2 = b * b;

  d = 4. * a * c;

  if(dabs(d) < 1e-10 * b2)
    {
      quad = d / 2. / dabs(b);

      if(!branch)
        quad *= -1.;
    }
  else
    {
      quad = sqrt(b2 - d);

      if(branch)
        quad *= -1.;

      quad -= b;
    }

  quad /= 2. * a;

  return quad;
}
