/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/rt/rt_CGmethod.c
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "../allvars.h"
#include "../proto.h"

#if defined(RT_ADVECT) && defined(RT_CGMETHOD)

#define MAX_ITER 100
#define ACCURACY 1.0e-2
#define EPSILON 1.0e-5
#define tiny 1e-37


face *VF;
point *DP;

static double *XVec;
static double *QVec, *DVec, *Residue, *Zvec;
static double *Kappa, *Diag, *Diag2;
static double c_light, dt;

void rt_cgmethod(tessellation * T, double dt)
{
  int i, j, iter;
  double alpha_cg, beta, delta_old, delta_new, sum, old_sum, min_diag, glob_min_diag, max_diag, glob_max_diag;
  double nH, mp_inv;
  double rel, res, maxrel, glob_maxrel;
  double DQ;

  c_light = CLIGHT / All.UnitVelocity_in_cm_per_s;
  mp_inv = 1.0 / (PROTONMASS / All.UnitMass_in_g * All.HubbleParam);

  set_cosmo_factors_for_current_time();
  if(All.ComovingIntegrationOn)
    {
      /* in comoving case, timestep is dloga at this point. Convert to dt */
      dt /= All.cf_hubble_a;
    }

  XVec = (double *) mymalloc("XVec", NumGas * sizeof(double));
  QVec = (double *) mymalloc("QVec", NumGas * sizeof(double));
  DVec = (double *) mymalloc("DVec", NumGas * sizeof(double));
  Kappa = (double *) mymalloc("Kappa", NumGas * sizeof(double));
  Residue = (double *) mymalloc("Residue", NumGas * sizeof(double));
  Diag = (double *) mymalloc("Diag", NumGas * sizeof(double));
  Zvec = (double *) mymalloc("Zvec", NumGas * sizeof(double));
  Diag2 = (double *) mymalloc("Diag2", NumGas * sizeof(double));

  /* loop over all directions */
  for(i = 0; i < RT_N_DIR; i++)
    {
      /* initialization for the CG method */
      for(j = 0; j < NumGas; j++)
        if(P[j].Type == 0)
          {
            QVec[j] = 0;
            DVec[j] = 0;
            Residue[j] = 0;
            Diag[j] = 0;
            Zvec[j] = 0;
            Diag2[j] = 0;

            XVec[j] = SphP[j].DensPhot[i];

            nH = (HYDROGEN_MASSFRAC * SphP[j].Density) * mp_inv;
            Kappa[j] = All.cf_a3inv * (SphP[j].nHI + tiny) * nH;

            if(All.ComovingIntegrationOn)
              Kappa[j] *= All.Time;
          }

      radtransfer_matrix_multiply(T, dt, XVec, Residue, Diag);

      /* Let's take the diagonal matrix elements as Jacobi preconditioner */

      for(j = 0, min_diag = MAX_REAL_NUMBER, max_diag = -MAX_REAL_NUMBER; j < NumGas; j++)
        if(P[j].Type == 0)
          {
            Residue[j] = SphP[j].DensPhot[i] - Residue[j];

            /* note: in principle we would have to substract the w_ii term, but this is always zero */
            if(Diag[j] < min_diag)
              min_diag = Diag[j];
            if(Diag[j] > max_diag)
              max_diag = Diag[j];

            Zvec[j] = Residue[j] / Diag[j];
            DVec[j] = Zvec[j];
          }

      MPI_Allreduce(&min_diag, &glob_min_diag, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&max_diag, &glob_max_diag, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

      delta_new = radtransfer_vector_multiply(Zvec, Residue);
      delta_old = delta_new;

      old_sum = radtransfer_vector_sum(XVec);

      if(ThisTask == 0)
        {
          printf("Begin RT_N_DIR %d/%d\n", i + 1, RT_N_DIR);
          printf("\nBegin CG iteration\nold |x|=%g, min-diagonal=%g, max-diagonal=%g\n", old_sum, glob_min_diag, glob_max_diag);
        }

      /* begin the CG method iteration */
      iter = 0;

      do
        {
          radtransfer_matrix_multiply(T, dt, DVec, QVec, Diag2);

          DQ = radtransfer_vector_multiply(DVec, QVec);
          if(DQ == 0)
            alpha_cg = 0;
          else
            alpha_cg = delta_new / DQ;


          for(j = 0, maxrel = 0; j < NumGas; j++)
            {
              XVec[j] += alpha_cg * DVec[j];
              Residue[j] -= alpha_cg * QVec[j];

              Zvec[j] = Residue[j] / Diag[j];

              rel = fabs(alpha_cg * DVec[j]) / (XVec[j] + tiny);
              if(rel > maxrel)
                maxrel = rel;
            }

          delta_old = delta_new;
          delta_new = radtransfer_vector_multiply(Zvec, Residue);

          sum = radtransfer_vector_sum(XVec);
          res = radtransfer_vector_sum(Residue);

          if(delta_old)
            beta = delta_new / delta_old;
          else
            beta = 0;

          for(j = 0; j < NumGas; j++)
            DVec[j] = Zvec[j] + beta * DVec[j];

          MPI_Allreduce(&maxrel, &glob_maxrel, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

          if(ThisTask == 0)
            {
              printf("radtransfer: iter=%3d  |res|/|x|=%12.6g  maxrel=%12.6g  |x|=%12.6g | res|=%12.6g\n", iter, res / sum, glob_maxrel, sum, res);
              myflush(stdout);
            }

          iter++;
        }
      while((res > ACCURACY * sum && iter < MAX_ITER) || iter < 2);

      if(ThisTask == 0)
        printf("\n");

      /* update the intensity */
      for(j = 0; j < NumGas; j++)
        if(P[j].Type == 0)
          {
            if(XVec[j] < 0)
              XVec[j] = 0;
            SphP[j].DensPhot[i] = XVec[j];
          }
    }

  myfree(Diag2);
  myfree(Zvec);
  myfree(Diag);
  myfree(Residue);
  myfree(Kappa);
  myfree(DVec);
  myfree(QVec);
  myfree(XVec);

}

/* internal product of two vectors */
double radtransfer_vector_multiply(double *a, double *b)
{
  int i;
  double sum, sumall;

  for(i = 0, sum = 0; i < NumGas; i++)
    if(P[i].Type == 0)
      sum += a[i] * b[i];

  MPI_Allreduce(&sum, &sumall, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return sumall;
}


double radtransfer_vector_sum(double *a)
{
  int i;
  double sum, sumall;

  for(i = 0, sum = 0; i < NumGas; i++)
    if(P[i].Type == 0)
      sum += fabs(a[i]);

  MPI_Allreduce(&sum, &sumall, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return sumall;
}


/* this function computes the vector b(out) given the vector x(in) such as Ax = b, where A is a matrix */
void radtransfer_matrix_multiply(tessellation * T, double dt, double *in, double *out, double *sum)
{
  int i, j, iside;
  double nx, ny, nz, nn;
  double dx, dy, dz, fac;
  int ri, li, target, other;
  double densphot;
  double velPhot[3];
  int point;
  double dist;
  double ainv, c_light;
  double volume_inv, area;

  VF = T->VF;
  DP = T->DP;

  c_light = CLIGHT / All.UnitVelocity_in_cm_per_s;

  if(All.ComovingIntegrationOn)
    ainv = 1.0 / All.Time;
  else
    ainv = 1.0;

  for(i = 0; i < T->Nvf; i++)
    {
      if(rt_face_get_normals(i))
        continue;

      if(rt_check_responsibility_of_this_task(VF[i].p1, VF[i].p2))
        continue;

      li = DP[VF[i].p1].index;
      ri = DP[VF[i].p2].index;

      if(li < 0 || ri < 0)
        continue;

      if(li >= NumGas && DP[VF[i].p1].task == ThisTask)
        li -= NumGas;

      if(ri >= NumGas && DP[VF[i].p2].task == ThisTask)
        ri -= NumGas;

      /* normal vector pointing to "right" state */
      nx = DP[VF[i].p2].x - DP[VF[i].p1].x;
      ny = DP[VF[i].p2].y - DP[VF[i].p1].y;
      nz = DP[VF[i].p2].z - DP[VF[i].p1].z;

      nn = sqrt(nx * nx + ny * ny + nz * nz);
      nx /= nn;
      ny /= nn;
      nz /= nn;

      velPhot[0] = rt_vec[j][0];
      velPhot[1] = rt_vec[j][1];
      velPhot[2] = rt_vec[j][2];

      dist = sqrt(velPhot[0] * velPhot[0] + velPhot[1] * velPhot[1] + velPhot[2] * velPhot[2]);

      velPhot[0] *= (CLIGHT / All.UnitVelocity_in_cm_per_s) / dist;
      velPhot[1] *= (CLIGHT / All.UnitVelocity_in_cm_per_s) / dist;
      velPhot[2] *= (CLIGHT / All.UnitVelocity_in_cm_per_s) / dist;

#ifdef TWODIMS
      double arclength = 2 * M_PI / RT_N_DIR;
      if(iside == 1)
        {
          nx = -nx;
          ny = -ny;
        }

      double vx = ny * velPhot[0] - nx * velPhot[1];
      double vy = nx * velPhot[0] + ny * velPhot[1];
      double phi = atan2(vy, vx);
      double alpha = phi - 0.5 * arclength;
      double beta = phi + 0.5 * arclength;
      if(beta < 0 && alpha < 0)
        fac = 0;
      else
        {
          if(beta > M_PI)
            beta = M_PI;
          if(alpha < 0 && beta > 0)
            alpha = 0;
          if(beta < 0)
            terminate("a");
          double p1 = velPhot[1] * nx - velPhot[0] * ny;
          double p2 = velPhot[0] * nx + velPhot[1] * ny;

          alpha -= phi;
          beta -= phi;
          fac = (p1 * (cos(beta) - cos(alpha)) + p2 * (sin(beta) - sin(alpha))) / arclength;
        }

      if(iside == 1)
        {
          nx = -nx;
          ny = -ny;
          fac = -fac;
        }
#else
      fac = velPhot[0] * nx + velPhot[1] * ny + velPhot[2] * nz;
#endif

      for(iside = 0; iside < 2; iside++)
        {
          if(iside == 0)
            {
              target = ri;
              other = li;
              fac = +fac;
              point = VF[i].p2;
            }
          else
            {
              target = li;
              other = ri;
              fac = -fac;
              point = VF[i].p1;
            }

          if(DP[point].task == ThisTask)
            volume_inv = 1 / SphP[target].Volume;
          else
            volume_inv = 1 / PrimExch[target].Volume;

          area = VF[i].area;

          if(fac > 0)
            sum[target] -= dt * volume_inv * c_light * ainv * area * fac;

          if(fac < 0)
            out[target] -= dt * volume_inv * c_light * ainv * area * fac * in[other];
        }

    }

  for(i = 0; i < NumGas; i++)
    if(P[i].Type == 0)
      {
        sum[i] += 1.0 + dt * c_light * ainv * Kappa[i];

        out[i] += in[i] * sum[i];
      }
}

#endif
