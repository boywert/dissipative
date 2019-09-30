/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/fld/fld.c
 * \date        09/2014
 * \author      Andreas Bauer (andreas.bauer@h-its.org)
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "../allvars.h"
#include "../proto.h"
#include "fld_proto.h"


#ifdef FLD


#define ACCURACY 1e-8
#define MAX_ITER 100

int fld(void)
{
  int i;

  TIMER_START(CPU_FLD) set_cosmo_factors_for_current_time();

#ifdef FLD_CONES
  fld_get_vectors();
#endif

  /* do only for highest time bin */
  if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
    {
      double timeeach = 0, timeall = 0, tstart = 0, tend = 0;
      All.fld_Radiation_Ti_endstep = All.Ti_Current;
      if(All.fld_Radiation_Ti_endstep - All.fld_Radiation_Ti_begstep == 0)
        return 0;

      double dt = (All.fld_Radiation_Ti_endstep - All.fld_Radiation_Ti_begstep) * All.Timebase_interval;

      double (*gamma_new)[NumGas] = (double (*)[NumGas]) mymalloc("gamma_new", sizeof(double) * NumGas * FLD_NCONES);
      double *u_new = (double *) mymalloc("u_new", sizeof(double) * NumGas);
      double *u_old = (double *) mymalloc("u_old", sizeof(double) * NumGas);


      tstart = second();

#ifdef FLD_CONES
      fld_update_n_gamma();
#endif

      exchange_primitive_variables();
      calculate_gradients();

      fld_update_kappa();
      //fld_compute_flux_limiter();

      fld_source(dt);

      for(i = 0; i < NumGas; i++)
        {
#ifndef FLD_CONES
          gamma_new[0][i] = SphP[i].n_gamma;
#else
          int j;
          for(j = 0; j < FLD_NCONES; j++)
            {
              gamma_new[j][i] = SphP[i].gammas[j];
            }
#endif

#ifndef FLD_MARSHAK
          double mu = 2.33;
          double u = SphP[i].Utherm * All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;
          double temp = GAMMA_MINUS1 / BOLTZMANN * u * PROTONMASS * mu;
          u_new[i] = temp;
          u_old[i] = temp;
#else
          u_new[i] = SphP[i].Utherm;
          u_old[i] = SphP[i].Utherm;
#endif
        }

      exchange_primitive_variables();

      //double tau = fld_sigma(gamma, u_new, dt) - 1.;

      double tau = 0.;

      int iter = 0;
      double change = 0, res = 0.;
      do
        {
          if(iter >= MAX_ITER)
            terminate("failed to converge\n");

          int cone;
          for(cone = 0; cone < FLD_NCONES; cone++)
            {
#ifndef FLD_HYPRE_IJ1
              fld_setup_matrix(cone, dt);
#endif
              fld_compute_coeff(cone, dt, gamma_new[cone], u_new, u_old, tau);

              TIMER_STOP(CPU_FLD) TIMER_START(CPU_FLD_MATRIX) fld_radtransfer(gamma_new[cone], dt);
            TIMER_STOP(CPU_FLD_MATRIX) TIMER_START(CPU_FLD)}

#ifndef FLD_NO_TEMP_UPDATE
          res = fld_update_state(dt, gamma_new[0], u_new, u_old, tau, &change);
#endif

          mpi_printf("FLD: finished PTC iteration %d: res %g, change %g, tau %g\n", iter, res, change, tau);

          tau = 0.5 * tau;

          iter++;
        }
      while(res > ACCURACY && change > ACCURACY);

      fld_update_gas(dt, gamma_new, u_new);

#ifdef FLD_CONES
      fld_update_n_gamma();
#endif

      exchange_primitive_variables();
      calculate_gradients();

#ifndef FLD_NO_TEMP_UPDATE
      fld_explicit_term(dt);
#endif

#ifndef FLD_SILENT
      tend = second();
      timeeach = timediff(tstart, tend);
      MPI_Allreduce(&timeeach, &timeall, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      mpi_printf("FLD: time consumed is %g \n", timeall);
#endif

      All.fld_Radiation_Ti_begstep = All.fld_Radiation_Ti_endstep;

      myfree(u_old);
      myfree(u_new);
      myfree(gamma_new);
    }

  TIMER_STOP(CPU_FLD) return 0;
}

#ifdef FLD_CONES
void fld_update_n_gamma()
{
  int i;

  for(i = 0; i < NumGas; i++)
    {
      SphP[i].n_gamma = 0.;

      int j;
      for(j = 0; j < FLD_NCONES; j++)
        {
          SphP[i].n_gamma += SphP[i].gammas[j];
        }
    }
}
#endif


void fld_compute_coeff(int cone, double dt, double *gamma_new, double *u_new, double *u_old, double tau)
{
  int i;

  double c_light = CLIGHT / All.UnitVelocity_in_cm_per_s;

#ifndef FLD_MARSHAK
  double mu = 2.33;
  double c_v = 1. / GAMMA_MINUS1 * BOLTZMANN / All.UnitPressure_in_cgs * All.UnitDensity_in_cgs / PROTONMASS / mu;
  double a_r = RAD_CONST / All.UnitEnergy_in_cgs * All.UnitLength_in_cm * All.UnitLength_in_cm * All.UnitLength_in_cm;
#endif

#ifdef FLD_CONES
  a_r *= 1. / (FLD_NCONES);
#endif

  double sigma = 1. + tau;

  for(i = 0; i < NumGas; i++)
    {
#ifndef FLD_NO_TEMP_UPDATE

#ifndef FLD_MARSHAK
      double planck = a_r * pow(u_new[i], 4);
      double dPlanck = 4. / u_new[i] * planck;

      double rho_cv = SphP[i].Density * c_v;
      double a = dt * c_light * SphP[i].Kappa_P;
#else
      double planck = All.Epsilon * u_new[i];
      double dPlanck = All.Epsilon;

      double rho_cv = SphP[i].Density;
      double a = dt * c_light * SphP[i].Kappa_P;
#endif
      double delta = 1. / (rho_cv * sigma + a * dPlanck);
      double f = a * dPlanck * delta;

#ifdef FLD_HYPRE_IJ1
      SphP[i].w = a - f * a + tau;
#else
      SphP[i].w[2 * NUMDIMS] += a - f * a + tau;
#endif

#ifdef FLD_MARSHAK
      if(P[i].ID >= FLD_LOWER_BOUNDARY_MINID && P[i].ID <= FLD_LOWER_BOUNDARY_MAXID)
        {
          SphP[i].w[2 * NUMDIMS] += dt * c_light * 2 * amr_area[Mesh.DP[i].level] / (4. + 3 * SphP[i].Kappa_P * amr_area[Mesh.DP[i].level]) / SphP[i].Volume;
        }
#endif

#ifndef FLD_CONES
      SphP[i].b = SphP[i].n_gamma + tau * gamma_new[i] + a * (planck) + f * (rho_cv * (u_old[i] - u_new[i]) - a * planck);
#else
      SphP[i].b = SphP[i].gammas[cone] + tau * gamma_new[i] + a * (planck) + f * (rho_cv * (u_old[i] - u_new[i]) - a * planck);

      if(i == 1000)
        {
          printf("%g %g %g %g %g\n", SphP[i].b, SphP[i].gammas[cone], a, planck, f);
        }
#endif

#else

#ifndef FLD_CONES
      SphP[i].b = SphP[i].n_gamma;
#else
      SphP[i].b = SphP[i].gammas[cone];
#endif

#endif
    }

}


double fld_update_state(double dt, double *gamma_new, double *u_new, double *u_old, double tau, double *change)
{
  int i;

  double res_u = 0., res_nl = 0.;
  double global_res_u, global_res_nl, global_res_change;

  double res_change = 0.;

  double c_light = CLIGHT / All.UnitVelocity_in_cm_per_s;

#ifndef FLD_MARSHAK
  double mu = 2.33;
  double c_v = 1. / GAMMA_MINUS1 * BOLTZMANN / All.UnitPressure_in_cgs * All.UnitDensity_in_cgs / PROTONMASS / mu;
  double a_r = RAD_CONST / All.UnitEnergy_in_cgs * All.UnitLength_in_cm * All.UnitLength_in_cm * All.UnitLength_in_cm;
#endif

  double sigma = 1. + tau;

  for(i = 0; i < NumGas; i++)
    {
      double gamma;
#ifndef FLD_CONES
      gamma = gamma_new[i];
#else
      gamma = 0.;
      int j;
      for(j = 0; j < FLD_NCONES; j++)
        {
          gamma += gamma_new[i + j * NumGas];
        }
#endif

#ifndef FLD_MARSHAK
      double planck = a_r * pow(u_new[i], 4);
      double dPlanck = 4. / u_new[i] * planck;

      double rho_cv = SphP[i].Density * c_v;
      double a = dt * c_light * SphP[i].Kappa_P;
#else
      double planck = All.Epsilon * u_new[i];
      double dPlanck = All.Epsilon;

      double rho_cv = SphP[i].Density;
      double a = dt * c_light * SphP[i].Kappa_P;
#endif
      double delta = 1. / (rho_cv * sigma + a * dPlanck);

      double du = delta * rho_cv * (u_old[i] - u_new[i]) - a * delta * (planck - gamma);

      u_new[i] = u_new[i] + du;

      res_change += fabs(SphP[i].Volume * rho_cv * du);

      res_u += fabs(SphP[i].Volume * rho_cv * u_new[i]);

#ifndef FLD_MARSHAK
      planck = a_r * pow(u_new[i], 4);
#else
      planck = All.Epsilon * u_new[i];
#endif
      res_nl += fabs(SphP[i].Volume * (rho_cv * (u_old[i] - u_new[i]) - a * (planck - gamma)));
    }

  MPI_Allreduce(&res_u, &global_res_u, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&res_nl, &global_res_nl, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&res_change, &global_res_change, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  *change = global_res_change;

  return global_res_nl / global_res_u;
}

void fld_update_gas(double dt, double gamma_new[FLD_NCONES][NumGas], MyFloat * u_new)
{
  int i;

  double local_source_min = MAX_REAL_NUMBER, local_source_max = MIN_REAL_NUMBER;
  double global_source_min, global_source_max;

  int idx;
  /* loop over all active cells and particles */
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];

      if(P[i].Type == 0 && P[i].Mass != 0 && P[i].ID != 0)
        {

          if(!gsl_finite(gamma_new[0][i]) || !gsl_finite(u_new[i]))
            {
              terminate("FLD: bad radiation update P[i].Mass=%g dtime=%g ID=%d Lambda=%g Kappa_diff=%g i=%d", P[i].Mass, dt, P[i].ID, SphP[i].Lambda, SphP[i].Kappa_diff, i);
            }


#ifndef FLD_CONES
          SphP[i].n_gamma = gamma_new[0][i];
#else
          int j;
          for(j = 0; j < FLD_NCONES; j++)
            {
              SphP[i].gammas[j] = gamma_new[j][i];
            }
#endif

          double source = 0.;
#ifndef FLD_NO_TEMP_UPDATE

#ifndef FLD_MARSHAK
          double mu = 2.33;
          double c_v = 1. / GAMMA_MINUS1 * BOLTZMANN / All.UnitPressure_in_cgs * All.UnitDensity_in_cgs / PROTONMASS / mu;
          source = u_new[i] * c_v * SphP[i].Density;
          SphP[i].Utherm = u_new[i] * c_v;
#else
          source = u_new[i] * SphP[i].Density;
          SphP[i].Utherm = u_new[i];
#endif

          SphP[i].Energy = SphP[i].Utherm * P[i].Mass + 0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);
#endif


          if(source < local_source_min)
            local_source_min = source;
          if(source > local_source_max)
            local_source_max = source;
        }
    }

  MPI_Reduce(&local_source_min, &global_source_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_source_max, &global_source_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  mpi_printf("total source range = (%g - %g) \n", global_source_min, global_source_max);

}

void fld_explicit_term(double dt)
{
  int i, k;

  double dvel;
  double local_ekin_min = MAX_REAL_NUMBER, local_ekin_max = MIN_REAL_NUMBER;
  double global_ekin_min, global_ekin_max;
  double ekin_before, ekin_after, ekin_ratio;

  double c_light = CLIGHT / All.UnitVelocity_in_cm_per_s;

  int idx;
  /* loop over all active cells and particles */
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(P[i].Type == 0 && P[i].Mass != 0 && P[i].ID != 0)
        {
          dt = (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;

#ifndef FLD_MARSHAK
          // momentum update
          // subtract kin. energy
          SphP[i].Energy -= 0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);
          ekin_before = 0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);

#ifdef FLD_TEST_BOUNDARY
          if(P[i].ID >= FLD_LOWER_BOUNDARY_MINID && P[i].ID <= FLD_LOWER_BOUNDARY_MAXID)
            {
              SphP[i].Grad.dngamma[0] = 0.;
              SphP[i].Grad.dngamma[1] = -All.fld_Flux / c_light / SphP[i].Lambda * SphP[i].Kappa_R;
              SphP[i].Grad.dngamma[2] = 0.;
            }
#endif
          // do the kick for gas cells
          for(k = 0; k < 2; k++)
            {
              dvel = -dt * SphP[i].Lambda * SphP[i].Grad.dngamma[k] / SphP[i].Density;
              // add velocity
              P[i].Vel[k] += dvel;
              // add momentum
              SphP[i].Momentum[k] += P[i].Mass * dvel;
              if(!gsl_finite(dvel))
                terminate("FLD: bad radiation pressure dvel=%g P[i].Mass=%g dtime=%g ID=%d Grad[k]=%g Lambda=%g k=%d", dvel, P[i].Mass, dt, P[i].ID, SphP[i].Grad.dngamma[k], SphP[i].Lambda, k);
            }
          // add kin. energy
          SphP[i].Energy += 0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);
          ekin_after = 0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);

          ekin_ratio = (ekin_before > 0) ? (ekin_after / ekin_before) : 1.0;

          if(ekin_ratio < local_ekin_min)
            local_ekin_min = ekin_ratio;

          if(ekin_ratio > local_ekin_max)
            local_ekin_max = ekin_ratio;
#endif

          MyFloat rel1 =
            dt * SphP[i].Lambda * (2. * SphP[i].Kappa_P / SphP[i].Kappa_R - 1.) * (P[i].Vel[0] * SphP[i].Grad.dngamma[0] + P[i].Vel[1] * SphP[i].Grad.dngamma[1] +
                                                                                   P[i].Vel[2] * SphP[i].Grad.dngamma[2]);

          SphP[i].Energy += rel1 * SphP[i].Volume;
          SphP[i].n_gamma -= rel1;

#ifdef FLD_CONES
          for(k = 0; k < FLD_NCONES; k++)
            {
              SphP[i].gammas[k] -= rel1 / FLD_NCONES;
            }
#endif

          SphP[i].Utherm = (SphP[i].Energy / P[i].Mass - 0.5 * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]));
        }
    }

  MPI_Reduce(&local_ekin_min, &global_ekin_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_ekin_max, &global_ekin_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  mpi_printf("kin. energy ratio = (%g - %g) \n", global_ekin_min, global_ekin_max);
}

void fld_update_kappa()
{
  int i;
  MyFloat minKappa = MAX_REAL_NUMBER, maxKappa = MIN_REAL_NUMBER;
  MyFloat global_minKappa, global_maxKappa;

  MyFloat kappa_R;
  MyFloat kappa_P;

  for(i = 0; i < NumGas; i++)
    {
#ifndef FLD_CONST_KAPPA
      MyFloat mu = 2.33;
      MyFloat u = SphP[i].Utherm * All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;
      MyFloat temp = GAMMA_MINUS1 / BOLTZMANN * u * PROTONMASS * mu;

      kappa_R = 0.0316 * pow(temp / 10., 2) * All.UnitMass_in_g / All.UnitLength_in_cm / All.UnitLength_in_cm;
      kappa_R *= SphP[i].Density;
      kappa_P = 0.1 * pow(temp / 10., 2) * All.UnitMass_in_g / All.UnitLength_in_cm / All.UnitLength_in_cm;
      kappa_P *= SphP[i].Density;
#else
      kappa_P = All.Kappa_P * SphP[i].Density;
      kappa_R = All.Kappa_R * SphP[i].Density;
#endif

      SphP[i].Kappa_R = kappa_R;
      SphP[i].Kappa_P = kappa_P;

      SphP[i].Kappa_diff = 1. / kappa_R;

      if(kappa_P < minKappa)
        minKappa = kappa_P;
      if(kappa_P > maxKappa)
        maxKappa = kappa_P;
    }

  MPI_Reduce(&minKappa, &global_minKappa, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&maxKappa, &global_maxKappa, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  mpi_printf("FLD: minKappa_P %g maxKappa_P %g\n", global_minKappa, global_maxKappa);
}

void fld_update_kappa_P(MyFloat * u_new)
{
  int i;

  MyFloat kappa_P;

  for(i = 0; i < NumGas; i++)
    {
#ifndef FLD_CONST_KAPPA
      MyFloat temp = u_new[i];

      kappa_P = 0.1 * pow(temp / 10., 2) * All.UnitMass_in_g / All.UnitLength_in_cm / All.UnitLength_in_cm;
      kappa_P *= SphP[i].Density;
#else
      kappa_P = All.Kappa_P * SphP[i].Density;
#endif

      SphP[i].Kappa_P = kappa_P;
    }
}

void fld_compute_flux_limiter(int cone)
{
  int i;
  double c_light = CLIGHT / All.UnitVelocity_in_cm_per_s;

  for(i = 0; i < NumGas; i++)
    {
#ifdef FLD_TEST_BOUNDARY
      if(P[i].ID >= FLD_LOWER_BOUNDARY_MINID && P[i].ID <= FLD_LOWER_BOUNDARY_MAXID)
        {
          SphP[i].Grad.dngamma[0] = 0.;
          SphP[i].Grad.dngamma[1] = -All.fld_Flux / c_light / SphP[i].Lambda * SphP[i].Kappa_R;
          SphP[i].Grad.dngamma[2] = 0.;
        }
#endif

      /* now calculate flux limiter */
#ifndef FLD_MARSHAK
      if(SphP[i].n_gamma > 0)
        {
          double beta = 1e-4;
#ifndef FLD_CONES
          double R = (sqrt(SphP[i].Grad.dngamma[0] * SphP[i].Grad.dngamma[0] +
                           SphP[i].Grad.dngamma[1] * SphP[i].Grad.dngamma[1] + SphP[i].Grad.dngamma[2] * SphP[i].Grad.dngamma[2]) + beta) / (SphP[i].n_gamma * SphP[i].Kappa_R);
#else
          double R = (sqrt(SphP[i].Grad.dngamma[0] * SphP[i].Grad.dngamma[0] +
                           SphP[i].Grad.dngamma[1] * SphP[i].Grad.dngamma[1] + SphP[i].Grad.dngamma[2] * SphP[i].Grad.dngamma[2]) + beta) / (SphP[i].n_gamma * SphP[i].Kappa_R) * FLD_NCONES;
#endif

          if(All.ComovingIntegrationOn)
            R /= All.Time;


          //SphP[i].Lambda = (2 + R) / (6 + 3 * R + R * R);
          SphP[i].Lambda = 1. / R * (1. / tanh(R) - 1. / R);

          if(SphP[i].Lambda < 1e-100)
            SphP[i].Lambda = 0;

          SphP[i].R2 = SphP[i].Lambda + SphP[i].Lambda * SphP[i].Lambda * R * R;
        }
      else
        {
          SphP[i].Lambda = 1. / 3.;
          SphP[i].R2 = 0.;      //FIXME what is the correct value?
        }

#else
      SphP[i].Lambda = 1. / 3;
#endif

      SphP[i].Kappa_diff = c_light * SphP[i].Lambda / SphP[i].Kappa_R;

#ifdef FLD_CONES
      SphP[i].Kappa_diff *= FLD_NCONES;
#endif
    }
}

double fld_limit_coeff(double *grad, double gamma, double kappa)
{
  double lambda;

#ifndef FLD_MARSHAK
  if(gamma > 0)
    {
      double beta = 1e-4;

      double R = (sqrt(grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]) + beta) / gamma * kappa;

      if(All.ComovingIntegrationOn)
        R /= All.Time;


      //SphP[i].Lambda = (2 + R) / (6 + 3 * R + R * R);
      lambda = 1. / R * (1. / tanh(R) - 1. / R);
    }
  else
    {
      lambda = 1. / 3.;
    }
#else
  lambda = 1. / 3;
#endif

  return lambda;
}

double fld_sigma(double *gamma, double *u_new, double dt)
{
  int i;

  double sigmaMax1 = 1.;
  double sigmaMax2 = 1.;
  double sigmaMax = 1.;


  double c_light = CLIGHT / All.UnitVelocity_in_cm_per_s;

  double mu = 2.33;
  double c_v = 1. / GAMMA_MINUS1 * BOLTZMANN / All.UnitPressure_in_cgs * All.UnitDensity_in_cgs / PROTONMASS / mu;
  double a_r = RAD_CONST / All.UnitEnergy_in_cgs * All.UnitLength_in_cm * All.UnitLength_in_cm * All.UnitLength_in_cm;


  for(i = 0; i < NumGas; i++)
    {
      double sigma1 = 1.;

      double planck = a_r * pow(u_new[i], 4);
      double dPlanck = 4. / u_new[i] * planck;

      double a = dt * c_light * SphP[i].Kappa_P;

      double b = 0.5 * (a * planck + a * dPlanck / SphP[i].Density / c_v * gamma[i]);
      double c = a * dPlanck / SphP[i].Density / c_v * (a * planck) + a * dPlanck * (-a * planck / SphP[i].Density / c_v);

      if(gamma[i] + 2 * b + c < 0.)
        {
          if(gamma[i] > 0.)
            {
              sigma1 = (sqrt(b * b - gamma[i] * c) - b) / gamma[i];
            }
          else
            {
              sigma1 = -c / (2. * b);
            }

        }

      sigmaMax1 = fmax(sigmaMax1, sigma1);
    }

  mpi_printf("sigma %g sigma1 %g sigma2 %g\n", sigmaMax, sigmaMax1, sigmaMax2);

  return sigmaMax1;
}

void fld_source(double dt)
{
  int i;

  //source terms at lower boundary
  face *VF = Mesh.VF;
  point *DP = Mesh.DP;
  for(i = 0; i < Mesh.Nvf; i++)
    {
      struct geometry geom;
      if(face_get_normals(&Mesh, i, &geom))
        continue;

#ifdef FLD_TEST_BOUNDARY
      if(DP[VF[i].p1].ID >= FLD_LOWER_BOUNDARY_MINID && DP[VF[i].p1].ID <= FLD_LOWER_BOUNDARY_MAXID && DP[VF[i].p1].ID == DP[VF[i].p2].ID)
        {
          int particle = DP[VF[i].p1].index;
          if(particle >= NumGas)
            particle -= NumGas;

          double source;
#ifdef FLD_MARSHAK
          source = All.fld_Flux * VF[i].area / SphP[particle].Volume * dt * 8. / (4. + 3. * All.Kappa_P * VF[i].area);
#else
          source = All.fld_Flux * VF[i].area / SphP[particle].Volume * dt;
#endif

#ifndef FLD_CONES
          SphP[particle].n_gamma += source;
#else
          int cone;
          for(cone = 0; cone < FLD_NCONES; cone++)
            {
              SphP[particle].gammas[cone] += source / FLD_NCONES;
            }
#endif

        }
      else if(DP[VF[i].p1].ID >= FLD_UPPER_BOUNDARY_MINID && DP[VF[i].p1].ID <= FLD_UPPER_BOUNDARY_MAXID && DP[VF[i].p1].ID == DP[VF[i].p2].ID)
        {
          int particle = DP[VF[i].p1].index;
          if(particle >= NumGas)
            particle -= NumGas;

          double r_ij = 0.5;
          double source = dt * SphP[particle].Kappa_diff * VF[particle].area / SphP[particle].Volume / r_ij * All.fld_n_gamma;
          SphP[particle].n_gamma += source;
        }
#endif

    }
}

#ifndef FLD_HYPRE_IJ1
void fld_setup_matrix(int cone, double dt)
{
  double c_light = CLIGHT / All.UnitVelocity_in_cm_per_s;

  int i;

#ifdef FLD_CONES
  int indices[4][6] = { {0, 4, 5, 6, 2, 3}, {4, 1, 2, 3, 7, 8}, {2, 4, 0, 1, 5, 7}, {4, 3, 6, 8, 0, 1} };
#endif

  for(i = 0; i < NumGas; i++)
    {
      double dx_ij, dy_ij, dz_ij, r_ij;
      double Kappa_i, Kappa_j;
      double area_i;
      double volume;
      double x_i, x_j, y_i, y_j, z_i, z_j;

      double gamma_i, gamma_j;

      area_i = amr_area[Mesh.DP[i].level];
      volume = SphP[i].Volume;
      Kappa_i = SphP[i].Kappa_diff;

      x_i = P[i].Pos[0];
      y_i = P[i].Pos[1];
      z_i = P[i].Pos[2];

      double diag = 1.;

      double grad_i[3];

      grad_i[0] = SphP[i].Grad.dngamma[0];
      grad_i[1] = SphP[i].Grad.dngamma[1];
      grad_i[2] = SphP[i].Grad.dngamma[2];

#ifndef FLD_CONES
      gamma_i = SphP[i].n_gamma;
#else
      gamma_i = SphP[i].gammas[cone];
#endif

      int j;

#ifdef FLD_CONES
      for(j = 0; j < 9; j++)
        SphP[i].w[j] = 0.;
#else
      for(j = 0; j < 5; j++)
        SphP[i].w[j] = 0.;
#endif

      for(j = 0; j < 2 * NUMDIMS; j++)
        {
          int other = Mesh.DP[i].neighbors[j];

          if(other == i)
            {
              SphP[i].w[j] = 0;
              continue;
            }

          int dp = other;

          double grad_j[3];

          if(other < Ngb_MaxPart)
            {
              Kappa_j = SphP[other].Kappa_diff;

              grad_j[0] = SphP[other].Grad.dngamma[0];
              grad_j[1] = SphP[other].Grad.dngamma[1];
              grad_j[2] = SphP[other].Grad.dngamma[2];

#ifndef FLD_CONES
              gamma_j = SphP[other].n_gamma;
#else
              gamma_j = SphP[other].gammas[cone];
#endif
            }
          else if(other >= Ngb_MaxPart + Mesh.nodes_total)
            {
              dp = other - Mesh.nodes_total;
              Kappa_j = PrimExch[Mesh.DP[dp].index].Kappa_diff;

              grad_j[0] = GradExch[other].dngamma[0];
              grad_j[1] = GradExch[other].dngamma[1];
              grad_j[2] = GradExch[other].dngamma[2];

#ifndef FLD_CONES
              gamma_j = PrimExch[other].n_gamma;
#else
              terminate("blubb");
#endif
            }
          else
            {
              terminate("bad mesh");
            }

          x_j = Mesh.DP[dp].x;
          y_j = Mesh.DP[dp].y;
          z_j = Mesh.DP[dp].z;


          dx_ij = nearest_x(x_j - x_i);
          dy_ij = nearest_y(y_j - y_i);
          dz_ij = nearest_z(z_j - z_i);


          r_ij = sqrt(dx_ij * dx_ij + dy_ij * dy_ij + dz_ij * dz_ij);
          double area_ij = area_i;




          double grad[3];

          grad[0] = 0.5 * (grad_i[0] + grad_j[0]);
          grad[1] = 0.5 * (grad_i[1] + grad_j[1]);
          grad[2] = 0.5 * (grad_i[2] + grad_j[2]);

#ifdef FLD_CONES
          double n[3];

          /*double x = 0.5 * (x_i+x_j);
             double y = 0.5 * (y_i+y_j);
             grad[0] = x-1.0078125;
             grad[1] = y-1.0078125;
             grad[2] = 0.; */

          fld_limit_cone(n, grad, cone);

          /*n[0] = (x-1.0078125);
             n[1] = (y-1.0078125);
             n[2] = 0.;

             double len = sqrt(n[0] * n[0] + n[1] * n[1]);

             n[0] /= len;
             n[1] /= len; */

#ifdef FLD_ANISOTROPIC_CIRCULAR
          double x = 0.5 * (x_i + x_j);
          double y = 0.5 * (y_i + y_j);

          double rx = x - 0.5 * boxSize_X;
          double ry = y - 0.5 * boxSize_Y;
          double alpha = atan2(ry, rx);
          acos(rx / sqrt(rx * rx + ry * ry));
          n[0] = -sin(alpha);
          n[1] = cos(alpha);

          if(rx * rx + ry * ry > 0.95 * 0.95)
            {
              double fac = exp(-(sqrt(rx * rx + ry * ry) - 0.95) * 20);
              n[0] *= fac;
              n[1] *= fac;
            }

          if(rx * rx + ry * ry < 0.05 * 0.05)
            {
              double fac = exp((sqrt(rx * rx + ry * ry) - 0.05) * 20);
              n[0] *= fac;
              n[1] *= fac;
            }

#endif
#endif

          if((Kappa_i + Kappa_j) > 0)
            {
              double coeff_mean_ij = 2.0 * (Kappa_i * Kappa_j) / (Kappa_i + Kappa_j);

              double gamma = 0.5 * (gamma_i + gamma_j);

              double lambda = fld_limit_coeff(grad, gamma, coeff_mean_ij);;
              coeff_mean_ij *= c_light * lambda;

#ifdef FLD_CONES
              coeff_mean_ij *= FLD_NCONES;
#endif

              SphP[i].Lambda = lambda;

              //coeff_mean_ij = 0.01;
              double w = coeff_mean_ij * area_ij * dt / volume / r_ij;


#ifdef FLD_CONES
              double rdotn = dx_ij * n[0] + dy_ij * n[1] + dz_ij * n[2];
              w *= rdotn / (amr_length[Mesh.DP[i].level]);

              double n0, n1;
              if(j < 2)
                {
                  n0 = n[0];
                  n1 = n[1];
                }
              else if(j < 4)
                {
                  n0 = n[1];
                  n1 = n[0];
                }
              else
                {
                  terminate("blubb");
                }

              /*if(n0 < 0)
                 {
                 n0 = - n0;
                 n1 = - n1;

                 w = -w;
                 } */

              int *ind = indices[j];

              SphP[i].w[ind[0]] += n0 * w;
              SphP[i].w[ind[1]] += -n0 * w;

              SphP[i].w[ind[2]] += n1 / 4. * w;
              SphP[i].w[ind[3]] += -n1 / 4. * w;

              SphP[i].w[ind[4]] += n1 / 4. * w;
              SphP[i].w[ind[5]] += -n1 / 4. * w;

              if(!gsl_finite(SphP[i].w[ind[0]]) || !gsl_finite(SphP[i].w[ind[1]]) || !gsl_finite(SphP[i].w[ind[2]]) || !gsl_finite(SphP[i].w[ind[3]]) || !gsl_finite(SphP[i].w[ind[4]])
                 || !gsl_finite(SphP[i].w[ind[5]]))
                {
                  terminate("blubb");
                }
#else
              SphP[i].w[j] = -w;
              diag += w;
#endif
            }
          else
            {
              SphP[i].w[j] = 0;
            }
        }

#ifdef FLD_TEST_BOUNDARY
      if(P[i].ID >= FLD_UPPER_BOUNDARY_MINID && P[i].ID <= FLD_UPPER_BOUNDARY_MAXID)
        {
          diag += Kappa_i * area_i * dt / volume / 0.5;
        }
#endif

      SphP[i].w[2 * NUMDIMS] += diag;
    }
}
#endif


#ifdef FLD_CONES
double rt_vec[FLD_NCONES][3];

void fld_limit_cone(double out[], double grad[], int cone)
{
#ifdef TWODIMS
  double phimax = 2. * 0.5 * (M_PI / FLD_NCONES);
#else
  double phimax = sqrt((4 * M_PI / FLD_NCONES) / M_PI);
#endif

  double len = sqrt(grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]);

  if(len > 0)
    {
      double n[3], g[3], m[3], gg[3];

      g[0] = -grad[0] / len;
      g[1] = -grad[1] / len;
      g[2] = -grad[2] / len;

      n[0] = rt_vec[cone][0];
      n[1] = rt_vec[cone][1];
      n[2] = rt_vec[cone][2];

      double cosphi = g[0] * n[0] + g[1] * n[1] + g[2] * n[2];
      double phi = acos(cosphi);


      if(phi > M_PI / 2.)
        {
          n[0] = -n[0];
          n[1] = -n[1];
          n[2] = -n[2];

          cosphi = -cosphi;
          phi = acos(cosphi);
        }

      if(phi > phimax)
        {
          gg[0] = n[1] * g[2] - n[2] * g[1];
          gg[1] = n[2] * g[0] - n[0] * g[2];
          gg[2] = n[0] * g[1] - n[1] * g[0];

          m[0] = gg[1] * n[2] - gg[2] * n[1];
          m[1] = gg[2] * n[0] - gg[0] * n[2];
          m[2] = gg[0] * n[1] - gg[1] * n[0];

          double mcoeff = sin(phimax);
          double ncoeff = cos(phimax);

          double sq = 1. / sqrt(1. - cosphi * cosphi);

          out[0] = mcoeff * sq * m[0] + ncoeff * n[0];
          out[1] = mcoeff * sq * m[1] + ncoeff * n[1];
          out[2] = mcoeff * sq * m[2] + ncoeff * n[2];
#ifdef TWODIMS
          out[2] = 0;
#endif
          return;
        }

      out[0] = g[0];
      out[1] = g[1];
      out[2] = g[2];

#ifdef TWODIMS
      out[2] = 0;
#endif

      return;
    }
  else
    {
      out[0] = rt_vec[cone][0];
      out[1] = rt_vec[cone][1];
      out[2] = rt_vec[cone][2];
    }
}

void fld_get_vectors(void)
{
  int i;

#ifdef TWODIMS
  for(i = 0; i < FLD_NCONES; i++)
    {
      rt_vec[i][0] = cos((M_PI / FLD_NCONES) * (i + 0.5));
      rt_vec[i][1] = sin((M_PI / FLD_NCONES) * (i + 0.5));
      rt_vec[i][2] = 0;

      mpi_printf("%d | %g %g %g\n", i, rt_vec[i][0], rt_vec[i][1], rt_vec[i][2]);
    }
#else
  double vec[3];

  for(i = 0; i < FLD_NCONES; i++)
    {
      terminate("TODO initialize only half sphere");
      pix2vec_ring(CONES_NSIDE, i, &vec[0]);

      rt_vec[i][0] = vec[0];
      rt_vec[i][1] = vec[1];
      rt_vec[i][2] = vec[2];

      mpi_printf("%d | %g %g %g\n", i, rt_vec[i][0], rt_vec[i][1], rt_vec[i][2]);
    }
#endif
}
#endif


#ifdef FLD_UPPER_BOUNDARY_MINID
#endif
#ifdef FLD_UPPER_BOUNDARY_MAXID
#endif
#ifdef FLD_LOWER_BOUNDARY_MINID
#endif
#ifdef FLD_LOWER_BOUNDARY_MAXID
#endif



#endif
