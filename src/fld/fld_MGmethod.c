/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/fld/fld_MGmethod.c
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

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "../allvars.h"
#include "../proto.h"
#include "fld_proto.h"



#ifdef FLD
#ifdef FLD_MG

#define MAX_ITER 10000
#define ACCURACY 1.0e-4

static double c_light, dt;

double fld_vector_sum(double *a);
void fld_smooth(int level, double *XVec);
double fld_residuum(int level, double *XVec);

void fld_multigrid(int level, double *XVec);
void fld_restrict(int level);
void fld_prolongate(int level, double *XVec);

void fld_output_line(double *XVec, FILE * f_rad2);
void fld_output(int level, FILE * f_rad2);
void fld_output_xval(int level, double *XVal, FILE * f_rad2);
void fld_output_b(int level, FILE * r_rad4);

double fld_sum_xval(int level, double *XVal);
double fld_sum_b(int level);

static FILE *f_rad = 0;
static FILE *f_rad1 = 0;
static FILE *f_rad2 = 0;
static FILE *f_rad3 = 0;
static FILE *f_rad4 = 0;

void fld_reorder_lists()
{
  int level;

  for(level = 1; level <= amr_maxlevel; level++)
    {
      int father = amr_lastinlevel[level - 1];

      int current = Ngb_Nodes[father].u.suns[0];
      int last = -1;

      int start0 = current;
      int start1 = current;

      do
        {
          start1 = current;
          do
            {
              if(last > -1)
                {
                  if(last >= Ngb_MaxPart)
                    {
                      Ngb_Nodes[last].nextinlevel = current;
                    }
                  else
                    {
                      Mesh.DP[last].nextinlevel = current;
                    }
                }
              else
                {
                  amr_lastinlevel[level] = current;
                }
              if(last >= Ngb_MaxPart)
                {
                  Ngb_Nodes[current].previnlevel = last;
                }
              else
                {
                  Mesh.DP[current].previnlevel = last;
                }

              last = current;
              if(last >= Ngb_MaxPart)
                {
                  current = Ngb_Nodes[current].neighbors[1];
                }
              else
                {
                  current = Mesh.DP[current].neighbors[1];
                }
            }
          while(current != start1 && current != last);

#ifdef TWODIMS
          if(last >= Ngb_MaxPart)
            {
              current = Ngb_Nodes[start1].neighbors[3];
            }
          else
            {
              current = Mesh.DP[start1].neighbors[3];
            }
#endif
        }
      while(current != start0 && current != start1);


      if(last >= Ngb_MaxPart)
        Ngb_Nodes[last].nextinlevel = -1;
      else
        Mesh.DP[last].nextinlevel = -1;
    }
}


void fld_radtransfer(double *XVec, double dt)
{
  int i;
  int iter = 0;
  double sum, res;

  c_light = CLIGHT / All.UnitVelocity_in_cm_per_s;
  dt = (All.fld_Radiation_Ti_endstep - All.fld_Radiation_Ti_begstep) * All.Timebase_interval;

  if(f_rad == 0)
    f_rad = fopen("rad_mg.txt", "w");

  if(f_rad1 == 0)
    f_rad1 = fopen("rad_out1.txt", "w");

  if(f_rad2 == 0)
    f_rad2 = fopen("rad_out2.txt", "w");

  if(f_rad3 == 0)
    f_rad3 = fopen("rad_out3.txt", "w");

  if(f_rad4 == 0)
    f_rad4 = fopen("rad_out4.txt", "w");

  for(i = 0; i < NumGas; i++)
    if(P[i].Type == 0)
      {
        XVec[i] = SphP[i].n_gamma;
      }

  double sum2 = fld_vector_sum(XVec);

  do
    {
      if(iter >= MAX_ITER)
        terminate("failed to converge\n");

#ifdef FLD_MG_GS
      fld_smooth(amr_maxlevel, XVec);
#else
      for(i = Ngb_MaxPart; i < Ngb_MaxPart + Ngb_NumNodes; i++)
        {
          XVec[i] = 0.;
        }
      fld_multigrid(amr_maxlevel, XVec);
#endif

      sum = fld_vector_sum(XVec);
      res = fld_residuum(amr_maxlevel, XVec);

#ifndef FLD_SILENT
      if(ThisTask == 0)
        {
          printf("FLD: radtransfer: iter=%3d  |res|/|x|=%12.6g |x|=%12.6g | res|=%12.6g\n\n", iter, res / sum, sum, res);
          fflush(stdout);
        }
#endif
      iter++;

      fld_output_line(XVec, f_rad1);

      fprintf(f_rad, "%g\n", res);

    }
  while((res > ACCURACY * sum) || iter < 2);


  sum = fld_vector_sum(XVec);
  printf("sum before: %g; sum after %g; defect: %g\n", sum2, sum, sum - sum2);

  fflush(f_rad);
  fflush(f_rad2);
  fflush(f_rad3);
  fflush(f_rad4);
}


void fld_multigrid(int level, double *XVec)
{
  if(level > 2)
    {
      fld_output_b(level, f_rad4);

      double res0 = fld_residuum(level, XVec);
      fld_output(level, f_rad2);
      fld_output_xval(level, XVec, f_rad3);

      fld_smooth(level, XVec);

      double res1 = fld_residuum(level, XVec);
      fld_output(level, f_rad2);
      fld_output_xval(level, XVec, f_rad3);

      fld_restrict(level);

      double b = fld_sum_b(level - 1);

      int gamma = 0;
      for(gamma = 0; gamma < 1; gamma++)
        {
          fld_multigrid(level - 1, XVec);
        }

      double x = fld_sum_xval(level - 1, XVec);

      fld_prolongate(level, XVec);

      fld_output(level, f_rad3);

      double res2 = fld_residuum(level, XVec);
      fld_output(level, f_rad2);
      fld_output_xval(level, XVec, f_rad3);

      fld_smooth(level, XVec);

      double res3 = fld_residuum(level, XVec);
      fld_output(level, f_rad2);
      fld_output_xval(level, XVec, f_rad3);



      printf("at level %d: initial res: %12.6g after smooth: %12.6g after MG iteration: %12.6g after smooth: %12.6g; b %g x %g\n", level, res0, res1, res2, res3, b, x);
    }
  else
    {
      /*fld_output_b(level, f_rad4);

         double res0 = fld_residuum(level, XVec);

         int node = amr_lastinlevel[level];

         XVec[node] = Ngb_Nodes[node].hydro.b;

         double res1 = fld_residuum(level, XVec);

         printf("at level %d: initial res: %g after smooth: %g, b-vector: %g kappa %g\n",level, res0,res1, Ngb_Nodes[node].hydro.b,  Ngb_Nodes[node].hydro.kappa);
       */

      /*int i;
         for(i = 0; i < 200; i++)
         {
         fld_smooth(level, XVec);
         } */

      int iter = 0;
      double sum;
      double res;

      do
        {
          if(iter >= MAX_ITER)
            terminate("failed to converge\n");

          fld_smooth(level, XVec);

          sum = fld_vector_sum(XVec);
          res = fld_residuum(level, XVec);
          iter++;
        }
      while((res > ACCURACY * sum) || iter < 2);
      printf("at level %d: final res %12.6g final sum %12.6g iterations %d\n", level, res, sum, iter);

    }
}

double fld_vector_sum(double *a)
{
  int i;
  double sum, sumall;

  for(i = 0, sum = 0; i < NumGas; i++)
    if(P[i].Type == 0)
      sum += fabs(a[i]);

  MPI_Allreduce(&sum, &sumall, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return sumall;
}

void fld_smooth(int level, double *XVec)
{
  int node = amr_lastinlevel[level];

  while(node >= 0)
    {
      double b;

      int *neighbors;

      double *ws;

      if(node < Ngb_MaxPart)
        {
          b = SphP[node].n_gamma;

          neighbors = Mesh.DP[node].neighbors;

          ws = SphP[node].w;
        }
      else if(node >= Ngb_MaxPart && node < Ngb_MaxPart + Ngb_NumNodes)
        {
          b = Ngb_Nodes[node].hydro.b;

          neighbors = Ngb_Nodes[node].neighbors;

          ws = Ngb_Nodes[node].hydro.w;
        }
      else
        {
          terminate("blubb");
        }

      double sum = 0.;
      double diag = 0.;

      int i;
      for(i = 0; i < 2 * NUMDIMS; i++)
        {
          int other = neighbors[i];

          if(other == node)
            continue;

          {
            double w = ws[i];

            sum += w * XVec[other];
          }

          diag = ws[2 * NUMDIMS];
        }

      XVec[node] = (b - sum) / (diag);


      if(node >= Ngb_MaxPart)
        node = Ngb_Nodes[node].nextinlevel;
      else
        node = Mesh.DP[node].nextinlevel;
    }
}

double fld_residuum(int level, double *XVec)
{
  double result = 0.;

  int node = amr_lastinlevel[level];
  while(node >= 0)
    {
      double sum = 0.;

      double b;

      int *neighbors;

      double *ws;

      if(node < Ngb_MaxPart)
        {
          b = SphP[node].n_gamma;

          ws = SphP[node].w;
          neighbors = Mesh.DP[node].neighbors;
        }
      else if(node >= Ngb_MaxPart && node < Ngb_MaxPart + Ngb_NumNodes)
        {
          b = Ngb_Nodes[node].hydro.b;

          ws = Ngb_Nodes[node].hydro.w;

          neighbors = Ngb_Nodes[node].neighbors;
        }
      else
        {
          terminate("blubb");
        }


      int i;
      for(i = 0; i < 2 * NUMDIMS; i++)
        {
          int other = neighbors[i];

          if(other == node)
            continue;

          {
            double w = ws[i];

            sum += w * XVec[other];

          }

        }

      double diag = ws[2 * NUMDIMS];

      double res = b - (diag * XVec[node] + sum);

      if(node < Ngb_MaxPart)
        {
          SphP[node].residuum = res;
        }
      else if(node >= Ngb_MaxPart && node < Ngb_MaxPart + Ngb_NumNodes)
        {
          Ngb_Nodes[node].hydro.residuum = res;
        }
      result += fabs(res);


      if(node >= Ngb_MaxPart)
        node = Ngb_Nodes[node].nextinlevel;
      else
        node = Mesh.DP[node].nextinlevel;
    }

  return result;
}

inline double fld_resval(int node)
{
  if(node < Ngb_MaxPart)
    {
      return SphP[node].residuum;
    }
  else if(node >= Ngb_MaxPart && node < Ngb_MaxPart + Ngb_NumNodes)
    {
      return Ngb_Nodes[node].hydro.residuum;
    }
  return 0.;
}

inline double fld_reskappa(int node)
{
  if(node < Ngb_MaxPart)
    {
      return SphP[node].Kappa_diff;
    }
  else if(node >= Ngb_MaxPart && node < Ngb_MaxPart + Ngb_NumNodes)
    {
      return Ngb_Nodes[node].hydro.kappa;
    }
  return 0.;
}

inline double fld_getw(int node, int dir, int side)
{
  int child;
  if(dir == -1)
    {
      child = Ngb_Nodes[Ngb_Nodes[node].neighbors[0]].u.suns[1];
    }
  else if(dir == 0)
    {
      child = Ngb_Nodes[node].u.suns[0];
    }
  else if(dir == 1)
    {
      child = Ngb_Nodes[node].u.suns[1];
    }

  if(child < Ngb_MaxPart)
    {
      return SphP[child].w[side];
    }
  else if(node >= Ngb_MaxPart && node < Ngb_MaxPart + Ngb_NumNodes)
    {
      return Ngb_Nodes[child].hydro.w[side];
    }
  return 0.;
}

void fld_comp_w(int level)
{
  int node = amr_lastinlevel[level];
  while(node >= 0)
    {

      double dx_ij, dy_ij, dz_ij, r_ij;
      double Kappa_i, Kappa_j;
      double area_i;
      double volume;
      double x_i, x_j, y_i, y_j, z_i, z_j;

      area_i = amr_area[level];
      volume = amr_volume[level];

      Kappa_i = Ngb_Nodes[node].hydro.kappa;

      x_i = Ngb_Nodes[node].Center[0];
      y_i = Ngb_Nodes[node].Center[1];
      z_i = Ngb_Nodes[node].Center[2];

      double diag = 1.;


      int j;
      for(j = 0; j < 2 * NUMDIMS; j++)
        {
          int other = Ngb_Nodes[node].neighbors[j];

          if(other == node)
            Ngb_Nodes[node].hydro.w[j] = 0;


          Kappa_j = Ngb_Nodes[other].hydro.kappa;

          x_j = Ngb_Nodes[other].Center[0];
          y_j = Ngb_Nodes[other].Center[1];
          z_j = Ngb_Nodes[other].Center[2];


          dx_ij = nearest_x(x_i - x_j);
          dy_ij = nearest_y(y_i - y_j);
          dz_ij = nearest_z(z_i - z_j);



          r_ij = sqrt(dx_ij * dx_ij + dy_ij * dy_ij + dz_ij * dz_ij);
          double area_ij = area_i;

          if((Kappa_i + Kappa_j) > 0)
            {
              double coeff_mean_ij = 2.0 * (Kappa_i * Kappa_j) / (Kappa_i + Kappa_j);
              double w = (coeff_mean_ij * area_ij) * dt / volume / r_ij;

              Ngb_Nodes[node].hydro.w[j] = -w;
              diag += w;

            }
          else
            {
              Ngb_Nodes[node].hydro.w[j] = 0;
            }
        }

      Ngb_Nodes[node].hydro.w[2 * NUMDIMS] = diag;
      node = Ngb_Nodes[node].nextinlevel;
    }
}

void fld_restrict(int level)
{
  int node = amr_lastinlevel[level - 1];
  while(node >= 0)
    {
#ifdef TWODIMS
      Ngb_Nodes[node].hydro.b = 1. / 4. * fld_resval(Ngb_Nodes[node].u.suns[2]) +
        1. / 8. * (fld_resval(Ngb_Nodes[Ngb_Nodes[node].neighbors[0]].u.suns[3]) + fld_resval(Ngb_Nodes[Ngb_Nodes[node].neighbors[3]].u.suns[0]) + fld_resval(Ngb_Nodes[node].u.suns[0]) +
                   fld_resval(Ngb_Nodes[node].u.suns[3])) + 1. / 16. * (fld_resval(Ngb_Nodes[Ngb_Nodes[Ngb_Nodes[node].neighbors[0]].neighbors[3]].u.suns[1]) +
                                                                        fld_resval(Ngb_Nodes[Ngb_Nodes[node].neighbors[3]].u.suns[1]) + fld_resval(Ngb_Nodes[Ngb_Nodes[node].neighbors[0]].u.suns[1]) +
                                                                        fld_resval(Ngb_Nodes[node].u.suns[1]));

      Ngb_Nodes[node].hydro.kappa = 1. / 4. * fld_reskappa(Ngb_Nodes[node].u.suns[2]) +
        1. / 8. * (fld_reskappa(Ngb_Nodes[Ngb_Nodes[node].neighbors[0]].u.suns[3]) + fld_reskappa(Ngb_Nodes[Ngb_Nodes[node].neighbors[3]].u.suns[0]) + fld_reskappa(Ngb_Nodes[node].u.suns[0]) +
                   fld_reskappa(Ngb_Nodes[node].u.suns[3])) + 1. / 16. * (fld_reskappa(Ngb_Nodes[Ngb_Nodes[Ngb_Nodes[node].neighbors[0]].neighbors[3]].u.suns[1]) +
                                                                          fld_reskappa(Ngb_Nodes[Ngb_Nodes[node].neighbors[3]].u.suns[1]) +
                                                                          fld_reskappa(Ngb_Nodes[Ngb_Nodes[node].neighbors[0]].u.suns[1]) + fld_reskappa(Ngb_Nodes[node].u.suns[1]));

#else
      Ngb_Nodes[node].hydro.b = 1. / 2. * fld_resval(Ngb_Nodes[node].u.suns[0]) + 1. / 4. * (fld_resval(Ngb_Nodes[Ngb_Nodes[node].neighbors[0]].u.suns[1]) + fld_resval(Ngb_Nodes[node].u.suns[1]));

      //Ngb_Nodes[node].hydro.kappa = 1./2. * fld_reskappa(Ngb_Nodes[node].u.suns[0]) +
      //1./4. *( fld_reskappa(Ngb_Nodes[Ngb_Nodes[node].neighbors[0]].u.suns[1]) + fld_reskappa(Ngb_Nodes[node].u.suns[1]) );

      Ngb_Nodes[node].hydro.w[0] = 1. / 8. * fld_getw(node, -1, 2) + 1. / 2. * fld_getw(node, -1, 1);
      Ngb_Nodes[node].hydro.w[1] = 1. / 8. * fld_getw(node, 1, 2) + 1. / 2. * fld_getw(node, 1, 0);
      Ngb_Nodes[node].hydro.w[2] = 3. / 4. * fld_getw(node, 0, 2) + 1. / 2. * fld_getw(node, 0, 0) + 1. / 2. * fld_getw(node, 0, 1);

      /*if(node == amr_lastinlevel[level-1])
         {
         printf("Node w[0] %g w[1] %g w[2] %g\n", Ngb_Nodes[node].hydro.w[0], Ngb_Nodes[node].hydro.w[1], Ngb_Nodes[node].hydro.w[2]);
         } */

#endif

      if(node >= Ngb_MaxPart)
        node = Ngb_Nodes[node].nextinlevel;
      else
        node = Mesh.DP[node].nextinlevel;
    }

#ifdef TWODIMS
  fld_comp_w(level - 1);
#endif
}

inline void fld_addval(double val, int node, double *XVec)
{
  if(node < Ngb_MaxPart)
    {
      SphP[node].residuum += val;
    }
  else if(node >= Ngb_MaxPart && node < Ngb_MaxPart + Ngb_NumNodes)
    {
      Ngb_Nodes[node].hydro.residuum += val;
    }

  XVec[node] += val;
}

void fld_prolongate(int level, double *XVec)
{
  int node = amr_lastinlevel[level];
  while(node >= 0)
    {

      if(node < Ngb_MaxPart)
        {
          SphP[node].residuum = 0.;
        }
      else if(node >= Ngb_MaxPart && node < Ngb_MaxPart + Ngb_NumNodes)
        {
          Ngb_Nodes[node].hydro.residuum = 0;
        }

      if(node >= Ngb_MaxPart)
        node = Ngb_Nodes[node].nextinlevel;
      else
        node = Mesh.DP[node].nextinlevel;
    }

  node = amr_lastinlevel[level - 1];
  while(node >= 0)
    {
#ifdef TWODIMS
      fld_addval(XVec[node], Ngb_Nodes[node].u.suns[2], XVec);

      fld_addval(1. / 2. * XVec[node], Ngb_Nodes[Ngb_Nodes[node].neighbors[0]].u.suns[3], XVec);
      fld_addval(1. / 2. * XVec[node], Ngb_Nodes[Ngb_Nodes[node].neighbors[3]].u.suns[0], XVec);
      fld_addval(1. / 2. * XVec[node], Ngb_Nodes[node].u.suns[0], XVec);
      fld_addval(1. / 2. * XVec[node], Ngb_Nodes[node].u.suns[3], XVec);

      fld_addval(1. / 4. * XVec[node], Ngb_Nodes[Ngb_Nodes[Ngb_Nodes[node].neighbors[0]].neighbors[3]].u.suns[1], XVec);
      fld_addval(1. / 4. * XVec[node], Ngb_Nodes[Ngb_Nodes[node].neighbors[0]].u.suns[1], XVec);
      fld_addval(1. / 4. * XVec[node], Ngb_Nodes[Ngb_Nodes[node].neighbors[3]].u.suns[1], XVec);
      fld_addval(1. / 4. * XVec[node], Ngb_Nodes[node].u.suns[1], XVec);
#else
      fld_addval(XVec[node], Ngb_Nodes[node].u.suns[0], XVec);

      fld_addval(1. / 2. * XVec[node], Ngb_Nodes[Ngb_Nodes[node].neighbors[0]].u.suns[1], XVec);

      fld_addval(1. / 2. * XVec[node], Ngb_Nodes[node].u.suns[1], XVec);
#endif

      if(node >= Ngb_MaxPart)
        node = Ngb_Nodes[node].nextinlevel;
      else
        node = Mesh.DP[node].nextinlevel;
    }
}

void fld_output(int level, FILE * f_rad2)
{
  int node = amr_lastinlevel[level];
  while(node >= 0)
    {

      double res = fld_resval(node);

      fprintf(f_rad2, "%g ", res);

      if(node >= Ngb_MaxPart)
        node = Ngb_Nodes[node].nextinlevel;
      else
        node = Mesh.DP[node].nextinlevel;
    }

  fprintf(f_rad2, "\n");
}

void fld_output_xval(int level, double *XVal, FILE * f_rad2)
{
  int node = amr_lastinlevel[level];
  while(node >= 0)
    {

      double res = XVal[node];

      fprintf(f_rad2, "%g ", res);

      if(node >= Ngb_MaxPart)
        node = Ngb_Nodes[node].nextinlevel;
      else
        node = Mesh.DP[node].nextinlevel;
    }

  fprintf(f_rad2, "\n");
}

double fld_sum_xval(int level, double *XVal)
{
  int node = amr_lastinlevel[level];
  double sum = 0.;
  while(node >= 0)
    {

      double res = XVal[node];
      sum += res;

      if(node >= Ngb_MaxPart)
        node = Ngb_Nodes[node].nextinlevel;
      else
        node = Mesh.DP[node].nextinlevel;
    }

  return sum;
}

double fld_sum_b(int level)
{
  int node = amr_lastinlevel[level];
  double sum = 0.;
  while(node >= 0)
    {

      if(node < Ngb_MaxPart)
        {
          sum += SphP[node].n_gamma;
        }
      else
        {
          sum += Ngb_Nodes[node].hydro.b;
        }


      if(node >= Ngb_MaxPart)
        node = Ngb_Nodes[node].nextinlevel;
      else
        node = Mesh.DP[node].nextinlevel;
    }

  return sum;
}

void fld_output_b(int level, FILE * r_rad4)
{
  int node = amr_lastinlevel[level];
  double out = 0.;
  while(node >= 0)
    {

      if(node < Ngb_MaxPart)
        {
          out = SphP[node].n_gamma;
        }
      else
        {
          out = Ngb_Nodes[node].hydro.b;
        }
      fprintf(f_rad4, "%g ", out);

      if(node >= Ngb_MaxPart)
        node = Ngb_Nodes[node].nextinlevel;
      else
        node = Mesh.DP[node].nextinlevel;
    }

  fprintf(f_rad4, "\n");

}

void fld_output_line(double *XVec, FILE * f_rad2)
{
  int i;

  double line[512];

  for(i = 0; i < NumGas; i++)
    {
      if(P[i].ID > 1000000000 && P[i].ID < 2000000000)
        {
          if((P[i].ID - 1000000000) % 512 == 256)
            {
              int element = (P[i].ID - 1000000000) / 512;
              line[element] = XVec[i];
            }
        }

      else if(P[i].ID == 256)
        {
          line[0] = XVec[i];
        }

      else if(P[i].ID == 2000000000 + 256 + 512 * 511)
        {
          line[511] = XVec[i];
        }

    }
  if(All.NumCurrentTiStep == 1)
    {
      for(i = 0; i < 512; i++)
        {
          fprintf(f_rad2, "%g ", line[i]);
        }
      fprintf(f_rad2, "\n");
    }
}

#endif
#endif
