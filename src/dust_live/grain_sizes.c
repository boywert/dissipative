/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dust_live/grain_sizes.c
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../proto.h"
#include "../allvars.h"

#ifdef DUST_LIVE
#ifdef DL_GRAIN_BINS

/* Overlap between [a_1, b_1] and [a_2, b_2]. */
double interval_overlap(double a_1, double b_1, double a_2, double b_2)
{
  return dmax(0.0, dmin(b_1, b_2) - dmax(a_1, a_2));
}

double x_1_ij(int i, int j, double adot, double dt)
{
  /* Use special handling for the boundary bins. */
  if(j == -1)
    {
      return GSD.Edges[i];
    }
  else if(j == DL_GRAIN_BINS)
    {
      return dmax(GSD.Edges[i], GSD.Edges[DL_GRAIN_BINS] - adot*dt);
    }
  else
    {
      return dmax(GSD.Edges[i], GSD.Edges[j] - adot*dt);
    }
}

double x_2_ij(int i, int j, double adot, double dt)
{
  /* Use special handling for the boundary bins. */
  if(j == -1)
    {
      return dmin(GSD.Edges[i+1], GSD.Edges[0] - adot*dt);
    }
  else if(j == DL_GRAIN_BINS)
    {
      return GSD.Edges[i+1];
    }
  else
    {
      return dmin(GSD.Edges[i+1], GSD.Edges[j+1] - adot*dt);
    }
}

double ind_overlap(double x_1, double x_2)
{
  return (x_2 >= x_1) ? 1.0 : 0.0;
}

enum gsd_integral_type
{
  GSD_NUMBER,
  GSD_RADIUS,
  GSD_MASS
};

#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
double radius_helper(double a, double adot, double dt, double mid)
{
  return a*a*a/3.0 + (adot*dt - mid)*a*a/2.0 - mid*adot*dt*a;
}

double mass_helper(double a, double adot, double dt, double mid)
{
   double res = pow(a, 5.0)/5.0;
   res += pow(a, 4.0)/4.0 * (3.0*adot*dt - mid);
   res += a*a*a * adot*dt * (adot*dt - mid);
   res += a*a*(adot*dt)*(adot*dt)/2.0 * (adot*dt - 3.0*mid);
   res -= a * adot*adot*adot * dt*dt*dt * mid;
   return res;
}
#endif

double integrate_i_to_j(int p, int i, int j, double adot, double dt, enum gsd_integral_type type)
{
    double x_1 = x_1_ij(i, j, adot, dt);
    double x_2 = x_2_ij(i, j, adot, dt);
    double is_overlap = ind_overlap(x_1, x_2);
    /* When tracking grains that evolve below a_min, need special handling to
     * ensure that we don't integrate contributions from grains that have a +
     * \dot{a} dt <= 0.0. */
    if(j == -1)
      {
        double a_star = -adot * dt;
        if((x_1 <= a_star) && (a_star <= x_2))
          x_1 = a_star;
        else if((x_1 <= a_star) && (x_2 <= a_star))
          is_overlap = 0.0;
      }

    double num;
    if(type == GSD_NUMBER)
      {
        num = DTP(p).NumGrains[i] * (x_2 - x_1) / GSD.Widths[i];
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
        num += DTP(p).BinSlopes[i] * ((x_2*x_2 - x_1*x_1)/2.0 - GSD.Midpoints[i] * (x_2 - x_1));
#endif
      }
    /* Caution: this is total radius of all grains that end up in bin j.
     * Divide elsewhere by the number of grains that end up in bin j if
     * you want an average radius. */
    else if(type == GSD_RADIUS)
      {
        num = DTP(p).NumGrains[i] * ((x_2*x_2 - x_1*x_1)/2.0 + adot*dt*(x_2 - x_1)) / GSD.Widths[i];
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
        num += DTP(p).BinSlopes[i] * (radius_helper(x_2, adot, dt, GSD.Midpoints[i]) - radius_helper(x_1, adot, dt, GSD.Midpoints[i]));
#endif
      }
    else if(type == GSD_MASS)
      {
        num = DTP(p).NumGrains[i] / GSD.Widths[i] * (pow(x_2 + adot*dt, 4.0) - pow(x_1 + adot*dt, 4.0)) / 4.0;
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
        num += DTP(p).BinSlopes[i] * (mass_helper(x_2, adot, dt, GSD.Midpoints[i]) - mass_helper(x_1, adot, dt, GSD.Midpoints[i]));
#endif
        num *= 4.0*M_PI/3.0 * GSD.InternalDensity;
      }
    else
      {
        terminate("DUST_LIVE: Invalid gsd_integral_type!\n");
      }

    if(is_overlap * num < 0.0)
      {
        num = 0.0;
      }
    return is_overlap * num;
}

#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
double bin_avg_mass(int j, double num_grains, double slope)
{
  double res = GSD.AvgMasses[j];
  /* Correction from piecewise-linear form.  Default to piecewise-constant
   * result if no grains are present. */
  if(num_grains > 0.0)
    {
      double edge_fac1 = (pow(GSD.Edges[j+1], 5.0) - pow(GSD.Edges[j], 5.0))/5.0;
      double edge_fac2 = (pow(GSD.Edges[j+1], 4.0) - pow(GSD.Edges[j], 4.0))*GSD.Midpoints[j]/4.0;
      res += 4.0*M_PI/3.0 * GSD.InternalDensity * slope/num_grains * (edge_fac1 - edge_fac2);
    }
  return res;
}
#endif

#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
double bin_slope_from_avg_a(int j, double avg_a, double num_grains)
{
  double num = num_grains*avg_a - num_grains/2.0 * (pow(GSD.Edges[j+1], 2.0) - pow(GSD.Edges[j], 2.0)) / GSD.Widths[j];
  double denom = (pow(GSD.Edges[j+1], 3.0) - pow(GSD.Edges[j], 3.0))/3.0 - GSD.Midpoints[j] * (pow(GSD.Edges[j+1], 2.0) - pow(GSD.Edges[j], 2.0))/2.0;
  return num/denom;
}

double bin_slope_from_mass(int j, double bin_mass, double num_grains)
{
  double num = bin_mass - num_grains * M_PI*GSD.InternalDensity/3.0 * (pow(GSD.Edges[j+1], 4.0) - pow(GSD.Edges[j], 4.0)) / GSD.Widths[j];
  double denom = (pow(GSD.Edges[j+1], 5.0) - pow(GSD.Edges[j], 5.0))/5.0 - GSD.Midpoints[j] * (pow(GSD.Edges[j+1], 4.0) - pow(GSD.Edges[j], 4.0))/4.0;
  denom *= 4.0*M_PI*GSD.InternalDensity/3.0;
  return num/denom;
}

void slope_limit_bin(int k, int j, double bin_mass)
{
  double N_j = DustParticle[k].NewNumGrains[j];
  double s_j = DustParticle[k].NewBinSlopes[j];
  double lhs_val = N_j/GSD.Widths[j] + s_j*(GSD.Edges[j] - GSD.Midpoints[j]);
  double rhs_val = N_j/GSD.Widths[j] + s_j*(GSD.Edges[j+1] - GSD.Midpoints[j]);
  if((lhs_val >= 0.0) && (rhs_val >= 0.0))
    {
      /* No slope-limiting necessary since grain size distribution always non-negative. */
      return;
    }

  double a, b, c, d;
  if(lhs_val < 0.0)
    {
      b = GSD.Edges[j] - GSD.Midpoints[j];
    }
  else if(rhs_val < 0.0)
    {
      b = GSD.Edges[j+1] - GSD.Midpoints[j];
    }

  a = 1.0 / GSD.Widths[j];
  c = M_PI*GSD.InternalDensity/3.0 * (pow(GSD.Edges[j+1], 4.0) - pow(GSD.Edges[j], 4.0)) / GSD.Widths[j];
  d = (pow(GSD.Edges[j+1], 5.0) - pow(GSD.Edges[j], 5.0))/5.0 - GSD.Midpoints[j] * (pow(GSD.Edges[j+1], 4.0) - pow(GSD.Edges[j], 4.0))/4.0;
  d *= 4.0*M_PI*GSD.InternalDensity/3.0;

  double det = a*d - b*c;
  if(det == 0.0)
    {
      terminate("DUST_LIVE: Unable to slope limit grain size distribution bin!\n");
    }

  DustParticle[k].NewNumGrains[j] = -b*bin_mass / det;
  DustParticle[k].NewBinSlopes[j] = a*bin_mass / det;

  if(DustParticle[k].NewNumGrains[j] < 0.0)
    {
      DustParticle[k].NewNumGrains[j] = 0.0;
      DustParticle[k].NewBinSlopes[j] = 0.0;
    }
}
#endif

/* j should be 0 or DL_GRAIN_BINS-1 */
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
void rebin_boundary_mass(int k, int j, double boundary_mass, double bin_mass)
{
  double a, b, c, d;
  double new_mass = boundary_mass + bin_mass;
  double c_fac = 1.0/2.0 * (pow(GSD.Edges[j+1], 2.0) - pow(GSD.Edges[j], 2.0)) / GSD.Widths[j];

  a = M_PI*GSD.InternalDensity/3.0 * (pow(GSD.Edges[j+1], 4.0) - pow(GSD.Edges[j], 4.0)) / GSD.Widths[j];
  b = (pow(GSD.Edges[j+1], 5.0) - pow(GSD.Edges[j], 5.0))/5.0 - GSD.Midpoints[j] * (pow(GSD.Edges[j+1], 4.0) - pow(GSD.Edges[j], 4.0))/4.0;
  b *= 4.0*M_PI*GSD.InternalDensity/3.0;

  d = (pow(GSD.Edges[j+1], 3.0) - pow(GSD.Edges[j], 3.0))/3.0 - GSD.Midpoints[j] * (pow(GSD.Edges[j+1], 2.0) - pow(GSD.Edges[j], 2.0))/2.0;

  double atot = DustParticle[k].NewNumGrains[j]*c_fac + DustParticle[k].NewBinSlopes[j]*d;
  double N_rebin = 0.0;
  if(j == 0)
    {
      N_rebin = boundary_mass / (4.0*M_PI*GSD.InternalDensity/3.0 * pow(GSD.Edges[0], 3.0));
      atot += N_rebin * GSD.Edges[0];
    }
  else if(j == DL_GRAIN_BINS-1)
    {
      N_rebin = boundary_mass / (4.0*M_PI*GSD.InternalDensity/3.0 * pow(GSD.Edges[DL_GRAIN_BINS], 3.0));
      atot += N_rebin * GSD.Edges[DL_GRAIN_BINS];
    }
  else
    {
      terminate("DUST_LIVE: Trying to rebin non-boundary mass!\n");
    }

  if(DustParticle[k].NewNumGrains[j] + N_rebin == 0.0)
    {
      /* There are no grains, so keep things the same. */
      return;
    }
  double new_a = atot / (DustParticle[k].NewNumGrains[j] + N_rebin);
  c = c_fac - new_a;
  double det = a*d - b*c;
  if(det == 0.0)
    {
      terminate("DUST_LIVE: Unable to rebin boundary mass!\n");
    }

  DustParticle[k].NewNumGrains[j] = (d*new_mass) / det;
  DustParticle[k].NewBinSlopes[j] = (-c*new_mass) / det;
}
#else
void rebin_boundary_mass(int k, int j, double boundary_mass)
{
  DustParticle[k].NewNumGrains[j] += boundary_mass / GSD.AvgMasses[j];
}
#endif

#if defined(DL_GROWTH) || defined(DL_SPUTTERING)
void update_grain_sizes(int k, double adot, double dt, double dt_code)
{
  int p = DustParticle[k].index;
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
  double bin_masses[DL_GRAIN_BINS];
#endif

  for(int j = 0; j < DL_GRAIN_BINS; j++)
    {
      DustParticle[k].NewNumGrains[j] = 0.0;
      for(int i = 0; i < DL_GRAIN_BINS; i++)
        {
          DustParticle[k].NewNumGrains[j] += integrate_i_to_j(p, i, j, adot, dt, GSD_NUMBER);
        }
    }

#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
  for(int j = 0; j < DL_GRAIN_BINS; j++)
    {
      bin_masses[j] = 0.0;
      for(int i = 0; i < DL_GRAIN_BINS; i++)
        {
          bin_masses[j] += integrate_i_to_j(p, i, j, adot, dt, GSD_MASS);
        }

      if(DustParticle[k].NewNumGrains[j] > 0.0)
        {
          DustParticle[k].NewBinSlopes[j] = bin_slope_from_mass(j, bin_masses[j], DustParticle[k].NewNumGrains[j]);
        }
      else
        {
          /* There are no grains in the bin, so we can fix the slope arbitrarily. */
          DustParticle[k].NewBinSlopes[j] = 0.0;
        }
    }
#endif

#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
  for(int j = 0; j < DL_GRAIN_BINS; j++)
    {
      /* Conserves mass. */
      slope_limit_bin(k, j, bin_masses[j]);
    }
#endif

  /* Special handling of boundary bins. */
  MyFloat M_m1 = 0.0;
  MyFloat M_N = 0.0;
  for(int i = 0; i < DL_GRAIN_BINS; i++)
    {
      M_m1 += integrate_i_to_j(p, i, -1, adot, dt, GSD_MASS);
      M_N += integrate_i_to_j(p, i, DL_GRAIN_BINS, adot, dt, GSD_MASS);
    }

#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
#if (DL_GRAIN_BINS == 1)
  /* If we're using only one bin, both boundary bins get rebinned in the same
   * place. */
  rebin_boundary_mass(k, 0, M_m1 + M_N, bin_masses[0]);
#else
  rebin_boundary_mass(k, 0, M_m1, bin_masses[0]);
  rebin_boundary_mass(k, DL_GRAIN_BINS-1, M_N, bin_masses[DL_GRAIN_BINS-1]);
#endif
#else
  rebin_boundary_mass(k, 0, M_m1);
  rebin_boundary_mass(k, DL_GRAIN_BINS-1, M_N);
#endif

#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
  /* Conserves mass. */
#if (DL_GRAIN_BINS == 1)
  slope_limit_bin(k, 0, bin_masses[0] + M_m1 + M_N);
#else
  slope_limit_bin(k, 0, bin_masses[0] + M_m1);
  slope_limit_bin(k, DL_GRAIN_BINS-1, bin_masses[DL_GRAIN_BINS-1] + M_N);
#endif
#endif

  /* Keep track of expected change in mass overall.  We need to store this for
   * later in case there aren't enough metals in surrounding gas. */
  DustParticle[k].DeltaMassExpected = 0.0;
  for(int j = 0; j < DL_GRAIN_BINS; j++)
    {
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
      double new_avg = bin_avg_mass(j, DustParticle[k].NewNumGrains[j], DustParticle[k].NewBinSlopes[j]);
      double old_avg = bin_avg_mass(j, DTP(p).NumGrains[j], DTP(p).BinSlopes[j]);
#else
      double new_avg = GSD.AvgMasses[j];
      double old_avg = GSD.AvgMasses[j];
#endif
      double delta_M_j = DustParticle[k].NewNumGrains[j]*new_avg - DTP(p).NumGrains[j]*old_avg;
      DustParticle[k].DeltaMassExpected += delta_M_j;

      double old_M = DTP(p).OrigMass;
      if((delta_M_j != 0.0) && (old_M > 0.0))
        {
          double tau = fabs(old_M * dt_code / delta_M_j);
          DTP(p).BinMassChgTau = dmin(DTP(p).BinMassChgTau, tau);
        }
    }
}
#endif

#ifdef DL_SNE_DESTRUCTION
void update_grain_sizes_sne(int k, double dt, double dt_code)
{
  int p = DustParticle[k].index;
  double sn_prefactor = DustParticle[k].LocalSNPrefactor;
  if((DTP(p).DustDensity * All.cf_a3inv) > (10.0 * All.PhysDensThresh))
    {
      /* Stronger dust destruction due to high local dust density. */
      sn_prefactor *= ((DTP(p).DustDensity * All.cf_a3inv) / (10.0 * All.PhysDensThresh));
    }

  double num_grains[DL_GRAIN_BINS], bin_masses[DL_GRAIN_BINS], new_masses[DL_GRAIN_BINS];
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
  double bin_slopes[DL_GRAIN_BINS];
#endif

  /* After handling grain growth and sputtering, account for any updates due to
   * supernova destruction.  At the end, NewNumGrains and NewBinSlopes include
   * contributions from all these processes. */
  for(int j = 0; j < DL_GRAIN_BINS; j++)
    {
      num_grains[j] = DTP(p).NumGrains[j];
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
      bin_slopes[j] = DTP(p).BinSlopes[j];
      bin_masses[j] = DTP(p).NumGrains[j] * bin_avg_mass(j, DTP(p).NumGrains[j], DTP(p).BinSlopes[j]);
#else
      bin_masses[j] = DTP(p).NumGrains[j] * GSD.AvgMasses[j];
#endif
      new_masses[j] = 0.0;
    }

  for(int j = 0; j < DL_GRAIN_BINS; j++)
    {
      double dNj = 0.0;
      for(int i = 0; i < DL_GRAIN_BINS; i++)
        {
          dNj += num_grains[i] * GSD.XiFrac[j][i] * GSD.Widths[j];
        }
      /* Take the new number of grains and subtract the starting number to get
       * the change. */
      dNj -= num_grains[j];

      double dNj_dt = sn_prefactor * dNj;
      DustParticle[k].NewNumGrains[j] = num_grains[j] + dNj_dt * dt;
      if(DustParticle[k].NewNumGrains[j] < 0.0)
        DustParticle[k].NewNumGrains[j] = 0.0;
    }

#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
  for(int j = 0; j < DL_GRAIN_BINS; j++)
    {
      double dMj = 0.0;
      for(int i = 0; i < DL_GRAIN_BINS; i++)
        {
          dMj += num_grains[i] * GSD.XiFrac[j][i] * 4.0*M_PI/3.0 * GSD.InternalDensity * (pow(GSD.Edges[j+1], 4) - pow(GSD.Edges[j], 4)) / 4.0;
        }
      dMj -= bin_masses[j];

      double dMj_dt = sn_prefactor * dMj;
      new_masses[j] = bin_masses[j] + dMj_dt * dt;
      /* If we should have zero mass in a bin, also set the number of grains to
       * zero. */
      if(new_masses[j] < 0.0)
        {
          new_masses[j] = 0.0;
          DustParticle[k].NewNumGrains[j] = 0.0;
        }

      if(DustParticle[k].NewNumGrains[j] > 0.0)
        {
          DustParticle[k].NewBinSlopes[j] = bin_slope_from_mass(j, new_masses[j], DustParticle[k].NewNumGrains[j]);
        }
      else
        {
          /* There are no grains in the bin, so we can fix the slope arbitrarily. */
          DustParticle[k].NewBinSlopes[j] = 0.0;
        }
    }
#else
  /* For the piecewise constant method, all grains in a bin are assumed to have
   * the same mass. */
  for(int j = 0; j < DL_GRAIN_BINS; j++)
    {
      new_masses[j] = DustParticle[k].NewNumGrains[j] * GSD.AvgMasses[j];
    }
#endif

#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
  for(int j = 0; j < DL_GRAIN_BINS; j++)
    {
      /* Conserves mass. */
      slope_limit_bin(k, j, new_masses[j]);
    }
#endif

  /* Do not zero DeltaMassExpected to start, because it was initialized by
   * growth and sputtering calculations in update_grain_sizes(). */
  DustParticle[k].DeltaMassExpected = 0.0;
  for(int j = 0; j < DL_GRAIN_BINS; j++)
    {
      DustParticle[k].DeltaMassExpected += (new_masses[j] - bin_masses[j]);

      double delta_M_j = new_masses[j] - bin_masses[j];
      double old_M = DTP(p).OrigMass;
      if((delta_M_j != 0.0) && (old_M > 0.0))
        {
          double tau = fabs(old_M * dt_code / delta_M_j);
          DTP(p).BinMassChgTau = dmin(DTP(p).BinMassChgTau, tau);
        }
    }
}
#endif

#if defined(DL_GROWTH) || defined(DL_SPUTTERING) || defined(DL_SNE_DESTRUCTION)
void correct_grain_conserved_quantities(int k)
{
  int p = DustParticle[k].index;

  double mass_frac = 0.0;
  if (DustParticle[k].DeltaMassExpected > 0.0)
    {
      /* We may not have been able to accrete the expected mass from gas. */
      mass_frac = DustParticle[k].DeltaMassActual / DustParticle[k].DeltaMassExpected;
    }
  else
    {
      /* We definitely are able to return the expected mass to gas. */
      mass_frac = 1.0;
    }

  for(int j = 0; j < DL_GRAIN_BINS; j++)
    {
      /* If we're trying to accrete mass and there weren't enough metals in the
       * surrounding gas, we only gain a fraction of the desired number of
       * grains. */
      double grain_diff = DustParticle[k].NewNumGrains[j] - DTP(p).NumGrains[j];
      DTP(p).NumGrains[j] += grain_diff * mass_frac;

#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
      /* This similar procedure for the slope ensures that the dust particle
       * changes in mass according to the mass transfer with gas. */
      double slope_diff = DustParticle[k].NewBinSlopes[j] - DTP(p).BinSlopes[j];
      DTP(p).BinSlopes[j] += slope_diff * mass_frac;
#endif
    }

  /* The Delta quantities below were calculated in the kernel over gas cells
   * and may be positive or negative depending on the sign of
   * DeltaMassExpected. */
  double m_new = P[p].Mass + DustParticle[k].DeltaMassActual;
  for(int l = 0; l < 3; l++)
    {
      double p_old = P[p].Mass * P[p].Vel[l];
      double p_new = p_old + DustParticle[k].DeltaMomentum[l];
      if(m_new > 0.0)
        P[p].Vel[l] = p_new / m_new;
    }
  P[p].Mass = m_new;
}
#endif

#if defined(DL_SHATTERING) || defined(DL_COAGULATION)
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
/* Results from Integrate[(x+y)^2 * x^3 * (Ni+si*(x-aic)) * (Nk+sk*(y-akc)), {x, ai, aip1}, {y, ak, akp1}]
 *
 * Note that for simplicity in Mathematica, N = NumGrains / GSD.Widths and that
 * 4pi/3 * rho_grain is omitted from the integral prefactor. */
double I_1_ik(int p, int i, int k)
{
  double ai = GSD.Edges[i];
  double aip1 = GSD.Edges[i+1];
  double aic = GSD.Midpoints[i];
  double ak = GSD.Edges[k];
  double akp1 = GSD.Edges[k+1];
  double akc = GSD.Midpoints[k];
  double Ni = DTP(p).NumGrains[i] / GSD.Widths[i];
  double Nk = DTP(p).NumGrains[k] / GSD.Widths[k];
  double si = DTP(p).BinSlopes[i];
  double sk = DTP(p).BinSlopes[k];

  double res = -((ak - akp1)*(-(pow(ai,4)*(360*pow(ai,3)*si*(2*Nk + (ak - 2*akc + akp1)*sk) +
              105*(Ni - aic*si)*(3*pow(ak,3)*sk + pow(ak,2)*(4*Nk - 4*akc*sk + 3*akp1*sk) + ak*akp1*(4*Nk - 4*akc*sk + 3*akp1*sk) +
                 pow(akp1,2)*(4*Nk - 4*akc*sk + 3*akp1*sk)) + 140*pow(ai,2)*
               (3*Ni*(2*Nk + (ak - 2*akc + akp1)*sk) + si*(6*ak*Nk + 6*akp1*Nk + 4*pow(ak,2)*sk - 6*ak*akc*sk + 4*ak*akp1*sk - 6*akc*akp1*sk + 4*pow(akp1,2)*sk -
                    3*aic*(2*Nk + (ak - 2*akc + akp1)*sk))) + 84*ai*(3*pow(ak,3)*si*sk + pow(ak,2)*(4*Nk*si + (8*Ni - 8*aic*si - 4*akc*si + 3*akp1*si)*sk) +
                 ak*(4*Ni*(3*Nk - 3*akc*sk + 2*akp1*sk) + si*(-4*aic*(3*Nk - 3*akc*sk + 2*akp1*sk) + akp1*(4*Nk - 4*akc*sk + 3*akp1*sk))) +
                 akp1*(4*Ni*(3*Nk - 3*akc*sk + 2*akp1*sk) + si*(-4*aic*(3*Nk - 3*akc*sk + 2*akp1*sk) + akp1*(4*Nk - 4*akc*sk + 3*akp1*sk)))))) +
         pow(aip1,4)*(360*pow(aip1,3)*si*(2*Nk + (ak - 2*akc + akp1)*sk) +
            105*(Ni - aic*si)*(3*pow(ak,3)*sk + pow(ak,2)*(4*Nk - 4*akc*sk + 3*akp1*sk) + ak*akp1*(4*Nk - 4*akc*sk + 3*akp1*sk) +
               pow(akp1,2)*(4*Nk - 4*akc*sk + 3*akp1*sk)) + 140*pow(aip1,2)*
             (3*Ni*(2*Nk + (ak - 2*akc + akp1)*sk) + si*(6*ak*Nk + 6*akp1*Nk + 4*pow(ak,2)*sk - 6*ak*akc*sk + 4*ak*akp1*sk - 6*akc*akp1*sk + 4*pow(akp1,2)*sk -
                  3*aic*(2*Nk + (ak - 2*akc + akp1)*sk))) + 84*aip1*(3*pow(ak,3)*si*sk + pow(ak,2)*(4*Nk*si + (8*Ni - 8*aic*si - 4*akc*si + 3*akp1*si)*sk) +
               ak*(4*Ni*(3*Nk - 3*akc*sk + 2*akp1*sk) + si*(-4*aic*(3*Nk - 3*akc*sk + 2*akp1*sk) + akp1*(4*Nk - 4*akc*sk + 3*akp1*sk))) +
               akp1*(4*Ni*(3*Nk - 3*akc*sk + 2*akp1*sk) + si*(-4*aic*(3*Nk - 3*akc*sk + 2*akp1*sk) + akp1*(4*Nk - 4*akc*sk + 3*akp1*sk)))))))/5040.;

  return 4.0/3.0 * M_PI * GSD.InternalDensity * res; /* internal mass * um^2 */
}

/* Results from Integrate[(x+y)^2 * (Nk+sk*(x-akc)) * (Nj+sj*(y-ajc)), {x, ak, akp1}, {y, aj, ajp1}]
 *
 * As in I_1_ik, in Mathematica we use N = NumGrains / GSD.Widths. */
double I_2_kj(int p, int k, int j)
{
  double Nk = DTP(p).NumGrains[k] / GSD.Widths[k];
  double Nj = DTP(p).NumGrains[j] / GSD.Widths[j];
#ifdef DL_SHATTERING_DETAILED_INTEGRALS
  double ak = GSD.Edges[k];
  double akp1 = GSD.Edges[k+1];
  double akc = GSD.Midpoints[k];
  double aj = GSD.Edges[j];
  double ajp1 = GSD.Edges[j+1];
  double ajc = GSD.Midpoints[j];
  double sk = DTP(p).BinSlopes[k];
  double sj = DTP(p).BinSlopes[j];

  double res = (9*pow(aj,4)*(ak - akp1)*sj*(2*Nk + (ak - 2*akc + akp1)*sk) +
       6*aj*(ak - akp1)*(Nj - ajc*sj)*(3*pow(ak,3)*sk + pow(ak,2)*(4*Nk - 4*akc*sk + 3*akp1*sk) + ak*akp1*(4*Nk - 4*akc*sk + 3*akp1*sk) +
          pow(akp1,2)*(4*Nk - 4*akc*sk + 3*akp1*sk)) + 4*pow(aj,3)*(ak - akp1)*
        (3*Nj*(2*Nk + (ak - 2*akc + akp1)*sk) + sj*(6*ak*Nk + 6*akp1*Nk + 4*pow(ak,2)*sk - 6*ak*akc*sk + 4*ak*akp1*sk - 6*akc*akp1*sk + 4*pow(akp1,2)*sk -
             3*ajc*(2*Nk + (ak - 2*akc + akp1)*sk))) + 3*pow(aj,2)*(3*pow(ak,4)*sj*sk + 12*pow(ak,2)*(Nj - ajc*sj)*(Nk - akc*sk) +
          4*pow(ak,3)*(Nk*sj + (2*Nj - (2*ajc + akc)*sj)*sk) + pow(akp1,2)*
           (-4*Nj*(3*Nk - 3*akc*sk + 2*akp1*sk) + sj*(akp1*(-4*Nk + 4*akc*sk - 3*akp1*sk) + 4*ajc*(3*Nk - 3*akc*sk + 2*akp1*sk)))) -
       ajp1*(9*pow(ajp1,3)*(ak - akp1)*sj*(2*Nk + (ak - 2*akc + akp1)*sk) +
          6*(ak - akp1)*(Nj - ajc*sj)*(3*pow(ak,3)*sk + pow(ak,2)*(4*Nk - 4*akc*sk + 3*akp1*sk) + ak*akp1*(4*Nk - 4*akc*sk + 3*akp1*sk) +
             pow(akp1,2)*(4*Nk - 4*akc*sk + 3*akp1*sk)) + 4*pow(ajp1,2)*(ak - akp1)*
           (3*Nj*(2*Nk + (ak - 2*akc + akp1)*sk) + sj*(6*ak*Nk + 6*akp1*Nk + 4*pow(ak,2)*sk - 6*ak*akc*sk + 4*ak*akp1*sk - 6*akc*akp1*sk + 4*pow(akp1,2)*sk -
                3*ajc*(2*Nk + (ak - 2*akc + akp1)*sk))) + 3*ajp1*(3*pow(ak,4)*sj*sk + 12*pow(ak,2)*(Nj - ajc*sj)*(Nk - akc*sk) +
             4*pow(ak,3)*(Nk*sj + (2*Nj - (2*ajc + akc)*sj)*sk) + pow(akp1,2)*
              (-4*Nj*(3*Nk - 3*akc*sk + 2*akp1*sk) + sj*(akp1*(-4*Nk + 4*akc*sk - 3*akp1*sk) + 4*ajc*(3*Nk - 3*akc*sk + 2*akp1*sk))))))/72.;
  return res; /* um^2 */
#else
  /* Results from piecewise constant approximation of the above integral, via
   * Mathematica for Integrate[(x+y)^2 * (Nk) * (Nj), {x,
   * ak, akp1}, {y, aj, ajp1}].  I_2_kj_Prefac contains the portion of the
   * integral apart from the Nk and Nj dependence. */
  double res = Nj*Nk*GSD.I_2_kj_Prefac[k][j];

  return res; /* um^2 */
#endif
}
#endif

#ifdef DL_SHATTERING
/* Equation 11 in Hirashita+ (2009) */
double sigma_helper(double M, double s_fac)
{
  return 0.30 * pow(s_fac + 1.0/M - 0.11, 0.13) / (s_fac + 1.0/M - 1.0);
}

void correct_fragment_range(double *af_min, double *af_max)
{
  /* Need to ensure size range of fragments lies within allowable range.
   * Intersect [af_min, af_max] with [a_min, a_max]. */
  double new_max = dmin(*af_max, GSD.Edges[DL_GRAIN_BINS]);
  double new_min = dmax(*af_min, GSD.Edges[0]);
  /* If the expected fragment size range does not overlap the allowable range,
   * just put all grains in the smallest bin. */
  if(new_max > new_min)
    {
      *af_max = new_max;
      *af_min = new_min;
    }
  else
    {
      *af_max = GSD.Edges[1];
      *af_min = GSD.Edges[0];
    }
}

/* Return the index corresponding to the grain size bin in which a grain of
 * mass m (in internal units) belongs. */
int bin_index_from_mass(double m)
{
  /* Follows gsl_histogram_find(). */
  int upper = DL_GRAIN_BINS;
  int lower = 0;
  /* If a lies outside of the grain size range, it will be put in the closest
   * bin. */
  double a = pow(3.0*m/(4.0*M_PI*GSD.InternalDensity), 1.0/3.0);

  while(upper - lower > 1)
    {
      int mid = (upper + lower) / 2;

      if(a >= GSD.Edges[mid])
        {
          lower = mid;
        }
      else
        {
          upper = mid;
        }
    }

  return lower;
}

/* Updates mass of grain fragments for size i produced by collisions of grains
 * of size k and j.  Follows notation in Section 2.3 of Hirashita+ (2009).
 *
 * Arguments:
 *   p: index into particle array
 *   k, j: indices of colliding bins
 *   v_kj: relative collision velocity, in um/s
 *   delta_bins: array of change in dust mass in internal mass units in each bin during time-step
 *   norm_fac: dimensionless prefactor (including dt) multiplying m_shat^kj(i) to give dMi in internal mass units
 */
void m_shat_ikj(int p, int k, int j, double v_kj, double *delta_bins, double norm_fac)
{
  /* Adopt consistent notation that k is the larger grain. */
  if(j > k)
    {
      int kj = j;
      j = k;
      k = kj;
    }

  double R_script = 1.0; /* dimensionless */
  double c_0 = (5.0+1.8)/2.0; /* km/s */
  double c_0_cms = c_0 * 1.0e5; /* cm/s */
  double c_0_int = c_0 * 1.0e9; /* um/s */
  double s_fac = (1.2+1.9)/2.0; /* dimensionless */
  double M_r = v_kj / c_0_int; /* dimensionless */

  double P_1 = (3.0e11 + 4.0e10) / 2.0; /* dyn/cm^2 */
  double P_v = (5.4e12 + 5.8e12) / 2.0; /* dyn/cm^2 */
  double phi_1 = P_1 / (All.GrainDensity * c_0_cms * c_0_cms); /* dimensionless */
  double M_1 = (2.0 * phi_1) / (1.0 + sqrt(1.0 + 4.0*s_fac*phi_1)); /* dimensionless */
  double sigma_1 = sigma_helper(M_1, s_fac); /* dimensionless */
  double sigma_1i = sigma_helper(M_r / (1.0 + R_script), s_fac); /*dimensionless */

  double M_shock = GSD.AvgMasses[j] * (1.0+2.0*R_script) / pow(1.0+R_script, 9.0/16.0) * pow(M_r*M_r / (sigma_1*M_1*M_1), 8.0/9.0) / pow(sigma_1i, 1.0/9.0); /* internal mass */

  double m_frag, v_cat;
  double af_min, af_max;

  /* First, handle fragmented mass from the projectile.  Hirashita+ (2009),
   * Section 2.3.1. */
  m_frag = GSD.AvgMasses[j]; /* internal mass */
  v_cat = c_0_int * pow(GSD.AvgMasses[k] / ((1.0 + 2.0*R_script) * GSD.AvgMasses[j]), 9.0/16.0) * sqrt(sigma_1) * pow(sigma_1i, 1.0/16.0) * (1.0 + R_script) * M_1; /* um/s */
  af_max = 0.22 * GSD.Midpoints[j] * v_cat/v_kj; /* um */
  af_min = 0.03 * af_max; /* um */
  correct_fragment_range(&af_min, &af_max);

  /* Distribute mass of fragmented projectile into resulting bins. */
  for(int i = 0; i < DL_GRAIN_BINS; i++)
    {
      /* Hirashita+ (2009), Equation 14. */
      double lim2 = dmin(af_max, GSD.Edges[i+1]); /* um */
      double lim1 = dmax(af_min, GSD.Edges[i]); /* um */
      /* Only if the fragmented size distribution overlaps bin i does mass
       * enter bin i. */
      if(lim2 > lim1)
        {
          double res_proj = m_frag * (pow(lim2, 0.7) - pow(lim1, 0.7)) / (pow(af_max, 0.7) - pow(af_min, 0.7)); /* internal mass */
          delta_bins[i] += norm_fac * res_proj;
        }
    }

  /* Next, handle fragmented mass from the target.  Hirashita+ (2009), Section
   * 2.3.2. */
  if(M_shock > 0.5*GSD.AvgMasses[k])
    {
      m_frag = GSD.AvgMasses[k]; /* internal mass */
      /* Note k and j are interchanged from above. */
      v_cat = c_0_int * pow(GSD.AvgMasses[j] / ((1.0 + 2.0*R_script) * GSD.AvgMasses[k]), 9.0/16.0) * sqrt(sigma_1) * pow(sigma_1i, 1.0/16.0) * (1.0 + R_script) * M_1; /* um/s */

      af_max = 0.22 * GSD.Midpoints[k] * v_cat/v_kj; /* um */
      af_min = 0.03 * af_max; /* um */
    }
  else
    {
      double M_ej = 0.4 * M_shock; /* internal mass */
      m_frag = M_ej; /* internal mass */
      double z_fac = 3.4;
      double af_max3 = (M_ej / GSD.InternalDensity) * 3.0 / (16.0*M_PI) * (z_fac+1.0) / (z_fac*z_fac*z_fac*(z_fac-2.0)); /* um^3 */
      af_max = pow(af_max3, 1.0/3.0); /* um */
      af_min = af_max * pow(P_1 / P_v, 1.47); /* um */

      /* Entire target has not fragmented.  Need to put leftover portion in
       * proper bin. */
      double m_left = GSD.AvgMasses[k] - M_ej; /* internal mass */
      int ii = bin_index_from_mass(m_left);
      delta_bins[ii] += norm_fac * m_left;
    }

  correct_fragment_range(&af_min, &af_max);

  /* Distribute mass of fragmented target into resulting bins. */
  for(int i = 0; i < DL_GRAIN_BINS; i++)
    {
      /* Hirashita+ (2009), Equation 14. */
      double lim2 = dmin(af_max, GSD.Edges[i+1]); /* um */
      double lim1 = dmax(af_min, GSD.Edges[i]); /* um */
      /* Only if the fragmented size distribution overlaps bin i does mass
       * enter bin i. */
      if(lim2 > lim1)
        {
          double res_targ = m_frag * (pow(lim2, 0.7) - pow(lim1, 0.7)) / (pow(af_max, 0.7) - pow(af_min, 0.7)); /* internal mass */
          delta_bins[i] += norm_fac * res_targ;
        }
    }
}
#endif

#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
void bin_state_post_shattering(int l, int i, double old_mass, double delta_mass)
{
  int p = DustParticle[l].index;

  double a, b, c, d;
  double new_mass = old_mass + delta_mass;
  double c_fac = 1.0/2.0 * (pow(GSD.Edges[i+1], 2.0) - pow(GSD.Edges[i], 2.0)) / GSD.Widths[i];

  a = M_PI*GSD.InternalDensity/3.0 * (pow(GSD.Edges[i+1], 4.0) - pow(GSD.Edges[i], 4.0)) / GSD.Widths[i];
  b = (pow(GSD.Edges[i+1], 5.0) - pow(GSD.Edges[i], 5.0))/5.0 - GSD.Midpoints[i] * (pow(GSD.Edges[i+1], 4.0) - pow(GSD.Edges[i], 4.0))/4.0;
  b *= 4.0*M_PI*GSD.InternalDensity/3.0;

  d = (pow(GSD.Edges[i+1], 3.0) - pow(GSD.Edges[i], 3.0))/3.0 - GSD.Midpoints[i] * (pow(GSD.Edges[i+1], 2.0) - pow(GSD.Edges[i], 2.0))/2.0;

  double old_atot = DTP(p).NumGrains[i]*c_fac + DTP(p).BinSlopes[i]*d;
  double new_a;
  if(delta_mass > 0.0)
    {
      double a_shat_i = 2.3/1.3 * (pow(GSD.Edges[i+1], -1.3) - pow(GSD.Edges[i], -1.3)) / (pow(GSD.Edges[i+1], -2.3) - pow(GSD.Edges[i], -2.3));
      double m_shat_i = 4.0*M_PI*GSD.InternalDensity/3.0 * -2.3/0.7 * (pow(GSD.Edges[i+1], 0.7) - pow(GSD.Edges[i], 0.7)) / (pow(GSD.Edges[i+1], -2.3) - pow(GSD.Edges[i], -2.3));
      double delta_N_i = delta_mass / m_shat_i;
      new_a = (old_atot + delta_N_i * a_shat_i) / (DTP(p).NumGrains[i] + delta_N_i);
    }
  else
    {
      /* Grains maintain their average size.  We should only lose mass if the
       * old number of grains is positive, but to be safe we revert to the
       * midpoint if necessary. */
      new_a = (DTP(p).NumGrains[i] > 0.0) ? old_atot / DTP(p).NumGrains[i] : GSD.Midpoints[i];
    }

  c = c_fac - new_a;
  double det = a*d - b*c;
  if(det == 0.0)
    {
      terminate("DUST_LIVE: Unable to solve for bin state after shattering!\n");
    }

  DustParticle[l].NewNumGrains[i] = (d*new_mass) / det;
  DustParticle[l].NewBinSlopes[i] = (-c*new_mass) / det;
}
#else
void bin_state_post_shattering(int l, int i, double old_mass, double delta_mass)
{
  DustParticle[l].NewNumGrains[i] = (old_mass + delta_mass) / GSD.AvgMasses[i];
}
#endif

#ifdef DL_COAGULATION
/* Updates masses of grains of size i in collisions between grains of sizes k
 * and j. */
void m_coag_ikj(int p, int k, int j, double *delta_bins, double norm_fac)
{
  double m_coag = GSD.AvgMasses[k] + GSD.AvgMasses[j];
  /* Could speed this up with some binary search if it actually matters.  If
   * the resulting mass doesn't end up in the N-1 smallest bins, put it in the
   * last bin. */
  for(int i = 0; i < DL_GRAIN_BINS-1; i++)
    {
      if((m_coag >= GSD.EdgeMasses[i]) && (m_coag < GSD.EdgeMasses[i+1]))
        {
          delta_bins[i] += norm_fac * m_coag;
          return;
        }
    }

  delta_bins[DL_GRAIN_BINS-1] += norm_fac * m_coag;
}
#endif

#if defined(DL_SHATTERING) || defined(DL_COAGULATION)
enum vel_type
{
  VEL_CNM,
  VEL_WIM
};

/* Returns velocity of a grain in bin i for given ISM phase type, in units of
 * um/s. */
double grain_vel(int i, enum vel_type type)
{
  if(type == VEL_CNM)
    {
      return GSD.VelocitiesCNM[i];
    }
  else if(type == VEL_WIM)
    {
      return GSD.VelocitiesWIM[i];
    }
  else
    {
      terminate("DUST_LIVE: Invalid vel_type!\n");
    }
}

double v_rel(int i1, int i2, enum vel_type type)
{
  double v1 = grain_vel(i1, type);
  double v2 = grain_vel(i2, type);
  double cos_theta = 2.0*get_random_number() - 1.0;

  return sqrt(v1*v1 + v2*v2 - 2.0*v1*v2*cos_theta);
}

double v_rel_eff(int i1, int i2, double cold_frac)
{
  return cold_frac * v_rel(i1, i2, VEL_CNM) + (1.0-cold_frac) * v_rel(i1, i2, VEL_WIM);
}

/* Perform collisional shattering or coagulation update for dust particle with
 * local index k over time-step dt in Gyr.  Here, type must be one of
 * GSD_SHATTERING or GSD_COAGULATION. */
void update_grain_sizes_shattering(int l, double dt, double dt_code, enum gsd_collision_type type)
{
  int p = DustParticle[l].index;

  double V_d = P[p].Mass / (DTP(p).DustDensity * All.cf_a3inv); /* internal length^3 */
  double V_d_um3 = V_d * pow(All.UnitLength_in_cm/All.HubbleParam * 1.0e4, 3.0); /* um^3 */
  double dt_s = dt * SEC_PER_GIGAYEAR;
  double mass_bins[DL_GRAIN_BINS], delta_bins[DL_GRAIN_BINS];
  for(int i = 0; i < DL_GRAIN_BINS; i++)
    {
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
      mass_bins[i] = bin_avg_mass(i, DTP(p).NumGrains[i], DTP(p).BinSlopes[i]) * DTP(p).NumGrains[i];
#else
      mass_bins[i] = GSD.AvgMasses[i] * DTP(p).NumGrains[i];
#endif
      delta_bins[i] = 0.0;
    }


  /* Precompute relative velocities.  To ensure that total mass is conserved,
   * it's important to manually ensure v_rels[i][j] = v_rels[j][i], which may
   * not otherwise be the case because of stochastic collision angles. */
  double v_rels[DL_GRAIN_BINS][DL_GRAIN_BINS];
  for(int i = 0; i < DL_GRAIN_BINS; i++)
    {
      for(int j = 0; j <= i; j++)
        {
          v_rels[i][j] = v_rel_eff(i, j, DTP(p).LocalCloudFrac); /* um/s */
          v_rels[j][i] = v_rels[i][j];
        }
    }

  /* Mass loss due to collisions with other grains. */
  for(int i = 0; i < DL_GRAIN_BINS; i++)
    {
      for(int k = 0; k < DL_GRAIN_BINS; k++)
        {
          double v_ik = v_rels[i][k]; /* um/s */
          double ind = 0.0;
#ifdef DL_SHATTERING
          if(type == GSD_SHATTERING)
            {
              ind = (v_ik > GSD.VelShat) ? 1.0 : 0.0;
            }
#endif
#ifdef DL_COAGULATION
          if(type == GSD_COAGULATION)
            {
              ind = (v_ik < GSD.VelCoag[i][k]) ? 1.0 : 0.0;
            }
#endif
          /* Skip update if velocity not within range. */
          if(ind == 0.0)
            {
              continue;
            }
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
          /* Can approximate I_1_ik ~ <m_i> * I_2_ik to ensure derivatives
           * analytically sum to zero. */
          /* double dMi_dt = -M_PI * v_ik * ind * I_1_kj(p, i, k) / V_d_um3; */
          double dMi_dt = -M_PI * v_ik * ind * I_2_kj(p, i, k) * GSD.AvgMasses[i] / V_d_um3; /* internal mass / s */
#else
          double dMi_dt = -M_PI * v_ik * ind * pow(GSD.Midpoints[i] + GSD.Midpoints[k], 2.0) * DTP(p).NumGrains[i] * DTP(p).NumGrains[k] * GSD.AvgMasses[i] / V_d_um3; /* internal mass / s */
#endif
          delta_bins[i] += dMi_dt * dt_s; /* internal mass */
        }
    }

  /* Mass gain due to collisions from two other grains.  We loop over the grain
   * sizes of colliding grains and distribute the resulting mass into bins in
   * m_shat_ikj. */
  for(int k = 0; k < DL_GRAIN_BINS; k++)
    {
      /* I believe Hirashita+ (2009) was double-counting things, by letting j <
       * DL_GRAIN_BINS. */
      for(int j = 0; j < DL_GRAIN_BINS; j++)
        {
          double v_kj = v_rels[k][j]; /* um/s */
          double ind = 0.0;
#ifdef DL_SHATTERING
          if(type == GSD_SHATTERING)
            {
              ind = (v_kj > GSD.VelShat) ? 1.0 : 0.0;
            }
#endif
#ifdef DL_COAGULATION
          if(type == GSD_COAGULATION)
            {
              ind = (v_kj < GSD.VelCoag[k][j]) ? 1.0 : 0.0;
            }
#endif
          /* Skip update if velocity not within range. */
          if(ind == 0.0)
            {
              continue;
            }
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
          double norm_fac = M_PI * v_kj * ind * I_2_kj(p, k, j) / V_d_um3 * dt_s / 2.0; /* dimensionless */
#else
          double norm_fac = M_PI * v_kj * ind * pow(GSD.Midpoints[k] + GSD.Midpoints[j], 2.0) * DTP(p).NumGrains[k] * DTP(p).NumGrains[j] / V_d_um3 * dt_s / 2.0; /* dimensionless */
#endif
#ifdef DL_SHATTERING
          if(type == GSD_SHATTERING)
            {
              m_shat_ikj(p, k, j, v_kj, delta_bins, norm_fac);
            }
#endif
#ifdef DL_COAGULATION
          if(type == GSD_COAGULATION)
            {
              m_coag_ikj(p, k, j, delta_bins, norm_fac);
            }
#endif
        }
    }

  for(int j = 0; j < DL_GRAIN_BINS; j++)
    {
      double delta_M_j = delta_bins[j];
      double old_M = DTP(p).OrigMass;
      if((delta_M_j != 0.0) && (old_M > 0.0))
        {
          double tau = fabs(old_M * dt_code / delta_M_j);
          DTP(p).BinMassChgTau = dmin(DTP(p).BinMassChgTau, tau);
        }
    }

  /* Ensure no bin loses too much mass. */
  for(int i = 0; i < DL_GRAIN_BINS; i++)
    {
      if(mass_bins[i] + delta_bins[i] < 0.0)
        {
          delta_bins[i] = -mass_bins[i];
        }
    }

  /* Using the mass changes in all bins, may need to rescale things to ensure
   * dust particle's mass is conserved. */
  double delta_mass = 0.0;
  for(int i = 0; i < DL_GRAIN_BINS; i++)
    {
      delta_mass += delta_bins[i];
    }

  if(delta_mass > 0.0)
    {
      double delta_mass_sub = 0.0;
      for(int i = 0; i < DL_GRAIN_BINS; i++)
        {
          delta_mass_sub += (delta_bins[i] > 0.0) ? delta_bins[i] : 0.0;
        }
      for(int i = 0; i < DL_GRAIN_BINS; i++)
        {
          if(delta_bins[i] > 0.0)
            {
              delta_bins[i] -= delta_mass * delta_bins[i] / delta_mass_sub;
            }
        }
    }
  else if(delta_mass < 0.0)
    {
      double delta_mass_sub = 0.0;
      for(int i = 0; i < DL_GRAIN_BINS; i++)
        {
          delta_mass_sub += (delta_bins[i] < 0.0) ? delta_bins[i] : 0.0;
        }
      for(int i = 0; i < DL_GRAIN_BINS; i++)
        {
          if(delta_bins[i] < 0.0)
            {
              delta_bins[i] += fabs(delta_mass * delta_bins[i] / delta_mass_sub);
            }
        }
    }

#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
  /* Apply average grain size heuristic to impose second constraint and
   * determine new number of grains and bin slope. */
  for(int i = 0; i < DL_GRAIN_BINS; i++)
    {
      bin_state_post_shattering(l, i, mass_bins[i], delta_bins[i]);
      slope_limit_bin(l, i, mass_bins[i] + delta_bins[i]);
    }
#else
  /* Can get new number of grains just from mass updates. */
  for(int i = 0; i < DL_GRAIN_BINS; i++)
    {
      bin_state_post_shattering(l, i, mass_bins[i], delta_bins[i]);
    }
#endif

  /* Because there's no mass transfer to or from gas, we don't need to worry
   * about limiting these terms. */
  for(int i = 0; i < DL_GRAIN_BINS; i++)
    {
      DTP(p).NumGrains[i] = DustParticle[l].NewNumGrains[i];
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
      DTP(p).BinSlopes[i] = DustParticle[l].NewBinSlopes[i];
#endif
    }
}
#endif

#endif
#endif
#endif
