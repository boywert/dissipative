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

/*! \file grav_external.c
 *  \brief special gravity routines for external forces
 *
 */

#ifdef EXTERNALGRAVITY
void gravity_external(void)
{
  int idx, i;

  TIMER_START(CPU_TREE);

#ifdef THERMAL_INSTABILITY_TEST
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      P[i].GravAccel[0] = 0.0;
      P[i].GravAccel[1] = 1.0;
      P[i].GravAccel[2] = 0.0;
    }

#endif

#ifdef MRT
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      P[i].GravAccel[0] = 0.0;
      P[i].GravAccel[1] = -1.46e-6 * pow(All.UnitTime_in_s, 2.0) / All.UnitLength_in_cm ;
      P[i].GravAccel[2] = 0.0;
    }
#endif

#ifdef DG_TEST_PROBLEM
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      double acc[3];
      dg_acceleration(SphP[i].Center[0], SphP[i].Center[1], SphP[i].Center[2], acc);
      P[i].GravAccel[0] = acc[0];
      P[i].GravAccel[1] = acc[1];
      P[i].GravAccel[2] = acc[2];
    }
#endif

#ifndef SPIRAL
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      P[i].ExtPotential = 0;    /* this is only needed for calculating total energy */
    }
#endif

#ifdef EXTERNALGY
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      P[i].GravAccel[1] += EXTERNALGY;
#ifdef EVALPOTENTIAL
      P[i].Potential += -(EXTERNALGY) * P[i].Pos[1];
#endif
      P[i].ExtPotential += -(EXTERNALGY) * P[i].Pos[1];
    }
#endif

#ifdef EXTERNALDISKPOTENTIAL
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      double dx, dy, r;

      dx = P[i].Pos[0] - boxHalf_X;
      dy = P[i].Pos[1] - boxHalf_Y;

      r = sqrt(dx * dx + dy * dy);

      double Sigma0 = 1.0 / (2 * M_PI);

      double y = r / (2);

      double dphidR = externaldisk_dphidR(r, NULL);
      double pot = externaldisk_potential(r);

      P[i].GravAccel[0] -= dphidR * dx / r;
      P[i].GravAccel[1] -= dphidR * dy / r;
#ifdef EVALPOTENTIAL
      P[i].Potential += pot;
#endif
      P[i].ExtPotential += pot;
    }
#endif

#ifdef EXTERNALSHEARBOX
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      double dz, b;

      dz = P[i].Pos[2] - boxHalf_Z;


      double Sigma0 = All.ShearBoxSigma0;
      double fg = All.ShearBoxFg;
      double mu = All.ShearBoxMu;


      b = 61.0 * (fg / 0.1 / mu) / (Sigma0 / 10.0);

      double f;

#ifdef SELFGRAVITY
      f = 1.0 / (1.0 / fg - 1.0);
#else
      f = fg;
#endif

      P[i].GravAccel[2] -= 2.0 * M_PI * All.G * Sigma0 * tanh(dz / b) / f;
      double pot = 2.0 * M_PI * All.G * b * Sigma0 * log(cosh(dz / b)) / f;

#ifdef EVALPOTENTIAL
      P[i].Potential += pot;
#endif
      P[i].ExtPotential += pot;
    }
#endif

#ifdef GRAVITY_TABLE 
  double sp_acc[3];

  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      double xi = P[i].Pos[0] - boxHalf_X;
      double yi = P[i].Pos[1] - boxHalf_Y;
      double zi = P[i].Pos[2] - boxHalf_Z;

      sp_acc[0] = 0.0;
      sp_acc[1] = 0.0;
      sp_acc[2] = 0.0;

      /* contributions from a multi component Galaxy model defined in an external file */

      grav_table_find_grav_acceleration(xi, yi, zi, sp_acc);
      
      P[i].GravAccel[0] = sp_acc[0];
      P[i].GravAccel[1] = sp_acc[1];
      P[i].GravAccel[2] = sp_acc[2];
    }
#endif

#ifdef GALPOT
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      double xi = P[i].Pos[0] - boxHalf_X;
      double yi = P[i].Pos[1] - boxHalf_Y;
      double zi = P[i].Pos[2] - boxHalf_Z;
      //double zi = 0.;

      double dPhidx[3]; 
      galpot_dPhidx(xi, yi, zi, All.Time, &dPhidx[0]);
      
      P[i].GravAccel[0] = -dPhidx[0];
      P[i].GravAccel[1] = -dPhidx[1];
      P[i].GravAccel[2] = -dPhidx[2];
    }
#endif

#ifdef EXTERNALSHEETY
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;
      if(P[i].Type == 0)
        P[i].GravAccel[1] = (-2 * M_PI * tanh((SphP[i].Center[1] - 1.0) / 0.1));
    }
#endif


#ifdef STATICISO
  {
    double r, m;
    double dx, dy, dz;

    for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
      {
        i = TimeBinsGravity.ActiveParticleList[idx];
        if(i < 0)
          continue;

#ifdef AXISYMMETRY
        dx = P[i].Pos[0];
        dy = P[i].Pos[1] - boxHalf_Y;
        dz = 0;

#ifdef CELL_CENTER_GRAVITY
        if(P[i].Type == 0)
        {
          dx = SphP[i].Center[0];
          dy = SphP[i].Center[1] - boxHalf_Y;
        }
#endif /* CELL_CENTER_GRAVITY */

#else

        dx = P[i].Pos[0] - boxHalf_X;
        dy = P[i].Pos[1] - boxHalf_Y;
        dz = P[i].Pos[2] - boxHalf_Z;

#ifdef CELL_CENTER_GRAVITY
        if(P[i].Type == 0) /* use center of mass, not mesh generating point in gas cell */
        {
          dx = SphP[i].Center[0] - boxHalf_X;
          dy = SphP[i].Center[1] - boxHalf_Y;
          dz = SphP[i].Center[2] - boxHalf_Z;
        }
#endif /* CELL_CENTER_GRAVITY */

#endif
        r = sqrt(dx * dx + dy * dy + dz * dz);

        if(r > ISO_R200)
          m = ISO_M200;
        else
          m = ISO_M200 * r / ISO_R200;

#ifdef ISO_FRACTION
        m *= ISO_FRACTION;
#endif

        if(r > 0)
          {
            P[i].GravAccel[0] += -All.G * m * dx / r / (r * r + ISO_Eps * ISO_Eps);
            P[i].GravAccel[1] += -All.G * m * dy / r / (r * r + ISO_Eps * ISO_Eps);
            P[i].GravAccel[2] += -All.G * m * dz / r / (r * r + ISO_Eps * ISO_Eps);
          }
      }
  }
#endif


#ifdef GROWING_DISK_POTENTIAL
  {
    double mdisk, dx, dy, dz, r, z, aR, az;

    growing_disk_init();

    mdisk = get_disk_mass(All.Time);

    for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
      {
        i = TimeBinsGravity.ActiveParticleList[idx];
        if(i < 0)
          continue;

        dx = P[i].Pos[0] - boxHalf_X;
        dy = P[i].Pos[1] - boxHalf_Y;
        dz = P[i].Pos[2] - boxHalf_Z;

        r = sqrt(dx * dx + dy * dy);
        z = fabs(dz);

        get_disk_forces(r, z, &aR, &az);

        aR *= mdisk;
        az *= mdisk;

        if(r > 0)
          {
            P[i].GravAccel[0] += -dx / r * aR;
            P[i].GravAccel[1] += -dy / r * aR;
            P[i].GravAccel[2] += -dz / z * az;
          }
      }
  }
#endif



#ifdef STATICNFW
  {
    double r, m;
    double dx, dy, dz;

    for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
      {
        i = TimeBinsGravity.ActiveParticleList[idx];
        if(i < 0)
          continue;

        dx = P[i].Pos[0] - boxHalf_X;
        dy = P[i].Pos[1] - boxHalf_Y;
        dz = P[i].Pos[2] - boxHalf_Z;

#ifdef CELL_CENTER_GRAVITY
        if(P[i].Type == 0) /* use center of mass, not mesh generating point in gas cell */
        {
          dx = SphP[i].Center[0] - boxHalf_X;
          dy = SphP[i].Center[1] - boxHalf_Y;
          dz = SphP[i].Center[2] - boxHalf_Z;
        }
#endif /* CELL_CENTER_GRAVITY */

        r = sqrt(dx * dx + dy * dy + dz * dz);
        m = enclosed_mass(r);
#ifdef NFW_DARKFRACTION
        m *= NFW_DARKFRACTION;
#endif
        if(r > 0)
          {
            P[i].GravAccel[0] += -All.G * m * dx / (r * r * r);
            P[i].GravAccel[1] += -All.G * m * dy / (r * r * r);
            P[i].GravAccel[2] += -All.G * m * dz / (r * r * r);
          }
#if defined (VS_TURB) || defined (AB_TURB)
        P[i].ExtPotential += get_turb_pot(P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
#endif
      }
  }
#endif



#ifdef STATICHQ
  {
    double r, m, a;
    double dx, dy, dz;

    for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
      {
        i = TimeBinsGravity.ActiveParticleList[idx];
        if(i < 0)
          continue;

        dx = P[i].Pos[0] - boxHalf_X;
        dy = P[i].Pos[1] - boxHalf_Y;
        dz = P[i].Pos[2] - boxHalf_Z;

#ifdef CELL_CENTER_GRAVITY
        if(P[i].Type == 0) /* use center of mass, not mesh generating point in gas cell */
        {
          dx = SphP[i].Center[0] - boxHalf_X;
          dy = SphP[i].Center[1] - boxHalf_Y;
          dz = SphP[i].Center[2] - boxHalf_Z;
        }
#endif /* CELL_CENTER_GRAVITY */

        r = sqrt(dx * dx + dy * dy + dz * dz);

        a = pow(All.G * HQ_M200 / (100 * All.Hubble * All.Hubble), 1.0 / 3) / HQ_C * sqrt(2 * (log(1 + HQ_C) - HQ_C / (1 + HQ_C)));

        m = HQ_M200 * pow(r / (r + a), 2);
#ifdef HQ_DARKFRACTION
        m *= HQ_DARKFRACTION;
#endif
        if(r > 0)
          {
            P[i].GravAccel[0] += -All.G * m * dx / (r * r * r);
            P[i].GravAccel[1] += -All.G * m * dy / (r * r * r);
            P[i].GravAccel[2] += -All.G * m * dz / (r * r * r);
          }
      }
  }
#endif


#ifdef CONSTANT_GRAVITY

  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      P[i].GravAccel[0] += 0.0;
      P[i].GravAccel[1] += -1.0;
      P[i].GravAccel[2] += 0.0;
    }

#endif


#ifdef CENTRAL_MASS_POTENTIAL
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      double fac, wp;
      double dx, dy, r, r2;
      double dz = 0;
      double h, h_inv, h3_inv, u;

      h = All.SofteningCentral * 2.8;

      dx = P[i].Pos[0] - boxHalf_X;
      dy = P[i].Pos[1] - boxHalf_Y;
#ifndef TWODIMS
      dz = P[i].Pos[2] - boxHalf_Z;
#endif

      r2 = dx * dx + dy * dy + dz * dz;
      r = sqrt(r2);

      //using spline softening
      if(r >= h)
        {
          fac = 1 / (r2 * r);
          wp = -1 / r;
        }
      else
        {
          h_inv = 1.0 / h;
          h3_inv = h_inv * h_inv * h_inv;
          u = r * h_inv;

          if(u < 0.5)
            {
              fac = h3_inv * (10.666666666667 + u * u * (32.0 * u - 38.4));
              wp = h_inv * (-2.8 + u * u * (5.333333333333 + u * u * (6.4 * u - 9.6)));
            }
          else
            {
              fac = h3_inv * (21.333333333333 - 48.0 * u + 38.4 * u * u - 10.666666666667 * u * u * u - 0.066666666667 / (u * u * u));
              wp = h_inv * (-3.2 + 0.066666666667 / u + u * u * (10.666666666667 + u * (-16.0 + u * (9.6 - 2.133333333333 * u))));
            }
        }

      P[i].GravAccel[0] -= All.G * All.CentralMass * fac * dx;
      P[i].GravAccel[1] -= All.G * All.CentralMass * fac * dy;
#ifndef TWODIMS
      P[i].GravAccel[2] -= All.G * All.CentralMass * fac * dz;
#endif

#ifdef EVALPOTENTIAL
      P[i].Potential += All.G * All.CentralMass * wp;
#endif
      P[i].ExtPotential += All.G * All.CentralMass * wp;

#ifdef SPECIAL_BOUNDARY
      if(P[i].ID <= -3)
        {
          P[i].GravAccel[0] = 0.0;
          P[i].GravAccel[1] = 0.0;
          P[i].GravAccel[2] = 0.0;
          continue;
        }
#endif
    }
#endif



#ifdef STAR_PLANET_POTENTIAL
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      double dx, dy, r;
      double first_mass, second_mass, dx_2, dy_2, r_2, soft;
      double indirect;

      soft = All.SofteningPlanet * 2.8;
      first_mass = (1 - All.MassRatio);
      second_mass = All.MassRatio;

      if(All.Time <= All.PlanetGrowthTime)
        second_mass *= sin(0.5 * M_PI * All.Time / All.PlanetGrowthTime) * sin(0.5 * M_PI * All.Time / All.PlanetGrowthTime);

      dx = P[i].Pos[0] - boxHalf_X;
      dy = P[i].Pos[1] - boxHalf_Y;
      dx_2 = P[i].Pos[0] - boxHalf_X - cos(All.Time);
      dy_2 = P[i].Pos[1] - boxHalf_Y - sin(All.Time);
      indirect = All.G * second_mass * (dx * cos(All.Time) + dy * sin(All.Time));

      r = sqrt(dx * dx + dy * dy);
      r_2 = sqrt(dx_2 * dx_2 + dy_2 * dy_2 + soft * soft);

      //or using spline softening
      r_2 = sqrt(dx_2 * dx_2 + dy_2 * dy_2);
      double u_2, fac_2, wp_2, soft_inv, soft3_inv;
      if(r_2 >= soft)
        {
          fac_2 = 1 / (r_2 * r_2 * r_2);
          wp_2 = -1 / r_2;
        }
      else
        {
          soft_inv = 1.0 / soft;
          soft3_inv = soft_inv * soft_inv * soft_inv;
          u_2 = r_2 * soft_inv;

          if(u_2 < 0.5)
            {
              fac_2 = soft3_inv * (10.666666666667 + u_2 * u_2 * (32.0 * u_2 - 38.4));
              wp_2 = soft_inv * (-2.8 + u_2 * u_2 * (5.333333333333 + u_2 * u_2 * (6.4 * u_2 - 9.6)));
            }
          else
            {
              fac_2 = soft3_inv * (21.333333333333 - 48.0 * u_2 + 38.4 * u_2 * u_2 - 10.666666666667 * u_2 * u_2 * u_2 - 0.066666666667 / (u_2 * u_2 * u_2));
              wp_2 = soft_inv * (-3.2 + 0.066666666667 / u_2 + u_2 * u_2 * (10.666666666667 + u_2 * (-16.0 + u_2 * (9.6 - 2.133333333333 * u_2))));
            }
        }

      double pot = -All.G * first_mass / r + All.G * second_mass * wp_2 + indirect;

#ifdef EVALPOTENTIAL
      P[i].Potential += pot;
#endif
      P[i].ExtPotential += pot;

#ifdef SPECIAL_BOUNDARY
      if(P[i].ID <= -3)
        {
          P[i].GravAccel[0] = 0.0;
          P[i].GravAccel[1] = 0.0;
          P[i].GravAccel[2] = 0.0;
          continue;
        }
#endif

      P[i].GravAccel[0] -= All.G * first_mass / r / r / r * dx + All.G * second_mass * fac_2 * dx_2 + All.G * second_mass * cos(All.Time);

      P[i].GravAccel[1] -= All.G * first_mass / r / r / r * dy + All.G * second_mass * fac_2 * dy_2 + All.G * second_mass * sin(All.Time);

      P[i].GravAccel[2] = 0;
    }
#endif



#ifdef BINARY_POTENTIAL
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      double dx_1, dy_1, dr_1, dx_2, dy_2, dr_2;
      double first_mass, second_mass;

      double soft;
      double soft1 = All.BinarySoftening * 2.8;
      double soft2 = soft1;

      double indirect = 0.0;

      double mu = All.BinaryMassRatio / (1.0 + All.BinaryMassRatio);
      first_mass = (1.0 - mu);
      second_mass = mu;
      if(All.Time < All.BinaryGrowthTime)
        second_mass *= sin(0.5 * M_PI * All.Time / All.BinaryGrowthTime) * sin(0.5 * M_PI * All.Time / All.BinaryGrowthTime);

      if(P[i].ID == 5)
        second_mass = 0;

      double x, y, vx, vy;
      circumstellar_solve_kepler(All.Time + M_PI, All.BinaryEccentricity, &x, &y, &vx, &vy);

      double x1, x2, y1, y2;
      if(All.BinaryBarycentricCoord)
        {
          x1 = boxHalf_X + mu * x;
          y1 = boxHalf_Y + mu * y;

          x2 = boxHalf_X - (1.0 - mu) * x;
          y2 = boxHalf_Y - (1.0 - mu) * y;
        }
      else
        {
          x1 = boxHalf_X;
          y1 = boxHalf_Y;

          x2 = boxHalf_X - x;
          y2 = boxHalf_Y - y;
        }
#ifdef CIRCUMSTELLAR_WBOUNDARIES
      if(All.BinaryBarycentricCoord == 0)
        soft1 = 0.0;
#endif

      dx_1 = P[i].Pos[0] - x1;
      dy_1 = P[i].Pos[1] - y1;
      dx_2 = P[i].Pos[0] - x2;
      dy_2 = P[i].Pos[1] - y2;
      double r_2 = sqrt(x2 * x2 + y2 * y2);

      if(All.BinaryBarycentricCoord)
        indirect = 0;
      else
        indirect = All.G * second_mass * (dx_1 * x2 + dy_1 * y2) / r_2 / r_2 / r_2;

      dr_1 = sqrt(dx_1 * dx_1 + dy_1 * dy_1);
      dr_2 = sqrt(dx_2 * dx_2 + dy_2 * dy_2);

      int k;
      double r, u, fac, wp, soft_inv, soft3_inv;
      double wp_1, wp_2, fac_1, fac_2;
      for(k = 0; k < 2; k++)
        {
          if(k == 0)
            {
              r = dr_1;
              soft = soft1;
            }
          else
            {
              r = dr_2;
              soft = soft2;
            }

          if(r >= soft)
            {
              fac = 1 / (r * r * r);
              wp = -1 / r;
            }
          else
            {
              soft_inv = 1.0 / soft;
              soft3_inv = soft_inv * soft_inv * soft_inv;
              u = r * soft_inv;

              if(u < 0.5)
                {
                  fac = soft3_inv * (10.666666666667 + u * u * (32.0 * u - 38.4));
                  wp = soft_inv * (-2.8 + u * u * (5.333333333333 + u * u * (6.4 * u - 9.6)));
                }
              else
                {
                  fac = soft3_inv * (21.333333333333 - 48.0 * u + 38.4 * u * u - 10.666666666667 * u * u * u - 0.066666666667 / (u * u * u));
                  wp = soft_inv * (-3.2 + 0.066666666667 / u + u * u * (10.666666666667 + u * (-16.0 + u * (9.6 - 2.133333333333 * u))));
                }
            }

          if(k == 0)
            {
              wp_1 = wp;
              fac_1 = fac;
            }
          else
            {
              wp_2 = wp;
              fac_2 = fac;
            }

        }

      double pot = -All.G * first_mass * wp_1 - All.G * second_mass * wp_2 + indirect;

#ifdef EVALPOTENTIAL
      P[i].Potential += pot;
#endif
      P[i].ExtPotential += pot;

#ifdef SPECIAL_BOUNDARY
      if(P[i].ID <= -3)
        {
          P[i].GravAccel[0] = 0.0;
          P[i].GravAccel[1] = 0.0;
          P[i].GravAccel[2] = 0.0;
          continue;
        }
#endif


      P[i].GravAccel[0] -= All.G * first_mass * fac_1 * dx_1 + All.G * second_mass * fac_2 * dx_2;

      P[i].GravAccel[1] -= All.G * first_mass * fac_1 * dy_1 + All.G * second_mass * fac_2 * dy_2;

      P[i].GravAccel[2] = 0;

      if(All.BinaryBarycentricCoord == 0)
        {
          P[i].GravAccel[0] -= All.G * second_mass * x2 / r_2 / r_2 / r_2;

          P[i].GravAccel[1] -= All.G * second_mass * y2 / r_2 / r_2 / r_2;
        }

    }
#endif

#if defined(CIRCUMSTELLAR) && defined(GRAVITY_FROM_STARS_PLANETS_ONLY)
  circumstellar_calc_gravity_from_stars_planets_only();
#endif

  TIMER_STOP(CPU_TREE);
}
#endif



#ifdef ONEDIMS_SPHERICAL
void gravity_monopole_1d_spherical()
{
  printf("Doing 1D gravity...\n");

  int i;
  double msum = All.CoreMass;

  for(i = 0; i < NumGas; i++)
    {
      double r0;
      if(i > 0)
        r0 = 0.5 * (P[i].Pos[0] + P[i - 1].Pos[0]);
      else
        r0 = All.CoreRadius;
      double dm = 4. / 3. * M_PI * (SphP[i].Center[0] * SphP[i].Center[0] * SphP[i].Center[0] - r0 * r0 * r0) * SphP[i].Density;
      double rad = SphP[i].Center[0];

      P[i].GravAccel[0] = -(msum + dm) * All.G / (rad * rad);

#ifdef EVALPOTENTIAL
      P[i].Potential = -(msum + dm) * All.G / rad;
#endif

      msum += P[i].Mass;

      P[i].GravAccel[1] = 0;
      P[i].GravAccel[2] = 0;
    }

  printf("... 1D gravity done.\n");
}
#endif





#ifdef STATICNFW
/*! auxiliary function for static NFW potential
 */
double enclosed_mass(double R)
{
  /* Eps is in units of Rs !!!! */

  if(R > Rs * NFW_C)
    R = Rs * NFW_C;

  return fac * 4 * M_PI * RhoCrit * Dc *
    (-(Rs * Rs * Rs * (1 - NFW_Eps + log(Rs) - 2 * NFW_Eps * log(Rs) + NFW_Eps * NFW_Eps * log(NFW_Eps * Rs)))
     / ((NFW_Eps - 1) * (NFW_Eps - 1)) +
     (Rs * Rs * Rs * (Rs - NFW_Eps * Rs - (2 * NFW_Eps - 1) * (R + Rs) * log(R + Rs) + NFW_Eps * NFW_Eps * (R + Rs) * log(R + NFW_Eps * Rs))) / ((NFW_Eps - 1) * (NFW_Eps - 1) * (R + Rs)));
}

#ifdef DG_EXTERNAL_ACCELERATION
void dg_acceleration(double x, double y, double z, double* acc)
{
	double r, M, a, dx, dy, dz;
	dx = x - boxHalf_X;
	dy = y - boxHalf_Y;
	dz = z - boxHalf_Z;
	r = sqrt( dx * dx + dy * dy + dz * dz );
	M = enclosed_mass(r);
    a = -All.G * M / r / r / r;

    acc[0] = a * dx;
    acc[1] = a * dy;
    acc[2] = a * dz;
}
#endif
#endif





#ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE
void calc_exact_gravity_for_particle_type(void)
{
  int i, idx;
#ifdef EXACT_GRAVITY_REACTION
  double *accx, *accy, *accz;
  accx = (double *)mymalloc("accx", All.MaxPartSpecial * sizeof(double));
  accy = (double *)mymalloc("accy", All.MaxPartSpecial * sizeof(double));
  accz = (double *)mymalloc("accz", All.MaxPartSpecial * sizeof(double));
#ifdef EVALPOTENTIAL
  double *pot;
  pot = (double *)mymalloc("pot", All.MaxPartSpecial * sizeof(double));
#endif
  int n;
  for(n = 0; n < All.MaxPartSpecial; n++)
    {
      accx[n] = accy[n] = accz[n] = 0.0;
#ifdef EVALPOTENTIAL
      pot[n] = 0.0;
#endif
    }
#endif

  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      double fac, wp;
      double dx, dy, dz, r, r2;
      double h, h_inv, h3_inv, u;
      int k;

      /* set softening to corresponding particle's softening length */
      h = All.ForceSoftening[All.SofteningTypeOfPartType[EXACT_GRAVITY_FOR_PARTICLE_TYPE]];

      for(k = 0; k < All.MaxPartSpecial; k++)
        {
          if(PartSpecialListGlobal[k].ID == P[i].ID)
            continue;
#ifdef REFINEMENT_AROUND_DM
          /* reset softening length to the length for this particle */
          h = DMPartListGlobal[k].softening;
#endif

          dx = P[i].Pos[0] - PartSpecialListGlobal[k].pos[0];
          dy = P[i].Pos[1] - PartSpecialListGlobal[k].pos[1];
          dz = P[i].Pos[2] - PartSpecialListGlobal[k].pos[2];

          r2 = dx * dx + dy * dy + dz * dz;
          r = sqrt(r2);

          //using spline softening
          if(r >= h)
            {
              fac = 1 / (r2 * r);
              wp = -1 / r;
            }
          else
            {
              h_inv = 1.0 / h;
              h3_inv = h_inv * h_inv * h_inv;
              u = r * h_inv;

              if(u < 0.5)
                {
                  fac = h3_inv * (10.666666666667 + u * u * (32.0 * u - 38.4));
                  wp = h_inv * (-2.8 + u * u * (5.333333333333 + u * u * (6.4 * u - 9.6)));
                }
              else
                {
                  fac = h3_inv * (21.333333333333 - 48.0 * u + 38.4 * u * u - 10.666666666667 * u * u * u - 0.066666666667 / (u * u * u));
                  wp = h_inv * (-3.2 + 0.066666666667 / u + u * u * (10.666666666667 + u * (-16.0 + u * (9.6 - 2.133333333333 * u))));
                }
            }

          P[i].GravAccel[0] -= All.G * PartSpecialListGlobal[k].mass * fac * dx;
          P[i].GravAccel[1] -= All.G * PartSpecialListGlobal[k].mass * fac * dy;
          P[i].GravAccel[2] -= All.G * PartSpecialListGlobal[k].mass * fac * dz;

#ifdef EVALPOTENTIAL
          P[i].Potential += All.G * PartSpecialListGlobal[k].mass * wp;
#endif
#ifdef EXACT_GRAVITY_REACTION
          /* avoid double counting */
          if(P[i].Type != EXACT_GRAVITY_FOR_PARTICLE_TYPE)
            {
              accx[k] += All.G * P[i].Mass * fac * dx;
              accy[k] += All.G * P[i].Mass * fac * dy;
              accz[k] += All.G * P[i].Mass * fac * dz;
#ifdef EVALPOTENTIAL
              pot[k] += All.G * P[i].Mass * wp;
#endif
            }
#endif
        }
    }
#ifdef EXACT_GRAVITY_REACTION
  double * buf = (double *) mymalloc("buf", All.MaxPartSpecial * sizeof(double));

  MPI_Allreduce(accx, buf, All.MaxPartSpecial, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for(n = 0; n < All.MaxPartSpecial; n++)
    accx[n] = buf[n];
  MPI_Allreduce(accy, buf, All.MaxPartSpecial, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for(n = 0; n < All.MaxPartSpecial; n++)
    accy[n] = buf[n];
  MPI_Allreduce(accz, buf, All.MaxPartSpecial, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for(n = 0; n < All.MaxPartSpecial; n++)
    accz[n] = buf[n];
#ifdef EVALPOTENTIAL
  MPI_Allreduce(pot, buf, All.MaxPartSpecial, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for(n = 0; n < All.MaxPartSpecial; n++)
    pot[n] = buf[n];
#endif
  myfree(buf);

  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;
      for(n = 0; n < All.MaxPartSpecial; n++)
        {
          if(PartSpecialListGlobal[n].ID == P[i].ID)
            {
              P[i].GravAccel[0] += accx[n];
              P[i].GravAccel[1] += accy[n];
              P[i].GravAccel[2] += accz[n];
#ifdef EVALPOTENTIAL
              P[i].Potential += pot[n];
#endif
            }
        }
    }

#ifdef EVALPOTENTIAL
  myfree(pot);
#endif
  myfree(accz);
  myfree(accy);
  myfree(accx);
#endif
}
#endif



#ifdef  EXACT_GRAVITY_FOR_PARTICLE_TYPE
void special_particle_create_list()
{
  struct special_particle_data *SpecialPartList;
  SpecialPartList = (struct special_particle_data *) mymalloc("SpecialPartList", All.MaxPartSpecial * sizeof(struct special_particle_data));


  int i, j, nsrc, nimport, ngrp;
  for(i = 0, nsrc = 0; i < NumPart; i++)
    {
      if(P[i].Type == EXACT_GRAVITY_FOR_PARTICLE_TYPE)
        {
          SpecialPartList[nsrc].ID = P[i].ID;

          SpecialPartList[nsrc].pos[0] = P[i].Pos[0];
          SpecialPartList[nsrc].pos[1] = P[i].Pos[1];
          SpecialPartList[nsrc].pos[2] = P[i].Pos[2];

          SpecialPartList[nsrc++].mass = P[i].Mass;
        }
    }

  for(j = 0; j < NTask; j++)
    Send_count[j] = nsrc;

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = 0;
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  /* exchange particle data */
  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&SpecialPartList[Send_offset[recvTask]],
                           Send_count[recvTask] * sizeof(struct special_particle_data), MPI_BYTE,
                           recvTask, TAG_DENS_A,
                           &PartSpecialListGlobal[Recv_offset[recvTask]],
                           Recv_count[recvTask] * sizeof(struct special_particle_data), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  myfree(SpecialPartList);

}

void special_particle_update_list()
{
  struct special_particle_data *SpecialPartList;
  SpecialPartList = (struct special_particle_data *) mymalloc("SpecialPartList", All.MaxPartSpecial * sizeof(struct special_particle_data));

  int i, j, nsrc, nimport, ngrp;
  for(i = 0, nsrc = 0; i < NumPart; i++)
    {
      if(P[i].Type == EXACT_GRAVITY_FOR_PARTICLE_TYPE)
        {
          SpecialPartList[nsrc].ID = P[i].ID;

          SpecialPartList[nsrc].pos[0] = P[i].Pos[0];
          SpecialPartList[nsrc].pos[1] = P[i].Pos[1];
          SpecialPartList[nsrc].pos[2] = P[i].Pos[2];

          SpecialPartList[nsrc++].mass = P[i].Mass;
        }
    }

  for(j = 0; j < NTask; j++)
    Send_count[j] = nsrc;

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = 0;
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  /* exchange particle data */
  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&SpecialPartList[Send_offset[recvTask]],
                           Send_count[recvTask] * sizeof(struct special_particle_data), MPI_BYTE,
                           recvTask, TAG_DENS_A,
                           &PartSpecialListGlobal[Recv_offset[recvTask]],
                           Recv_count[recvTask] * sizeof(struct special_particle_data), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  myfree(SpecialPartList);

}
#endif

#if defined(ACCRETE_ONTO_CENTRAL_POTENTIAL) && defined(CENTRAL_MASS_POTENTIAL)
void accrete_onto_central_potential()
{
  int idx, i;
  double dx, dy, dz, r2, r, racc2, vr, new_racc, nma;
  double racc_global;
  double dm_local = 0.0;
  double dm_global;
  int count_local = 0;
  int count_global;

  mpi_printf("CENTRALPOTENTIAL: Accreting all cells within = %g\n", All.CentralAccretionRadius);
  racc2 = All.CentralAccretionRadius * All.CentralAccretionRadius;
  new_racc = 2 * sqrt(racc2);

  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      dx = P[i].Pos[0] - boxHalf_X;
      dy = P[i].Pos[1] - boxHalf_Y;
      dz = P[i].Pos[2] - boxHalf_Z;

      r2 = dx * dx + dy * dy + dz * dz;

      // Accretion criterion based on density and position
      if(r2 < racc2 && P[i].Mass > 0.1 * All.TargetGasMass)
      {
        count_local += 1;

        P[i].Mass *= 0.5;
        dm_local += P[i].Mass;

#ifdef MHD
        double Emag = 0.5 * (SphP[i].B[0] * SphP[i].B[0] + SphP[i].B[1] * SphP[i].B[1] + SphP[i].B[2] * SphP[i].B[2]) * SphP[i].Volume * All.cf_atime;
        SphP[i].Energy -= Emag;
#endif
        SphP[i].Energy *= 0.5;
#ifdef MHD
        SphP[i].Energy += Emag;
#endif

        SphP[i].Momentum[0] *= 0.5;
        SphP[i].Momentum[1] *= 0.5;
        SphP[i].Momentum[2] *= 0.5;

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
        SphP[i].Entropy *= 0.5;
#endif
#ifdef MAXSCALARS
        for(int s = 0; s < N_Scalar; s++)
          *(MyFloat *) (((char *) (&SphP[i])) + scalar_elements[s].offset_mass) *= 0.5;
#endif

      }
      else if(All.HighestOccupiedTimeBin == All.HighestActiveTimeBin && r2 < 4.0 * racc2)
      {
        r = sqrt(r2);
        vr = (dx * P[i].Vel[0] + dy * P[i].Vel[1] + dz * P[i].Vel[2]) / r;
        nma = vr / SphP[i].Csnd;
        if(nma > -2.0 && r < new_racc)
          new_racc = r;
      }
    }

    if(All.HighestOccupiedTimeBin == All.HighestActiveTimeBin)
    {
      MPI_Allreduce(&new_racc, &racc_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      if(racc_global > All.CentralAccretionRadius)
      {
        All.CentralAccretionRadius = racc_global;
        mpi_printf("CENTRALPOTENTIAL: Adjusting new accretion radius to %g\n",  All.CentralAccretionRadius);
      }
    }

  // Sum accreted mass
  MPI_Allreduce(&count_local, &count_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&dm_local, &dm_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  All.CentralMass += dm_global;

  mpi_printf("CENTRALPOTENTIAL: Accreted %g from %i cells onto central mass potential, total = %g\n", dm_global, count_global, All.CentralMass );
}
#endif
