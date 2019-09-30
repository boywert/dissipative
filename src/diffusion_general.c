/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/conduction.c
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

#include "arepoconfig.h"

#ifdef DIFFUSION

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "voronoi.h"
#include "proto.h"

/*! \file conduction.c
*  \brief Computes isotropic/anisotropic diffusion bases on an implicit diffusion solver
*  \brief for general quantities. Based on implementation by V.Springel and V.Pauz
*/

#define MAX_ITER_DIFFUSION 500
#define ACCURACY_DIFFUSION 1e-6

static void diffusion_exchange_in_vector(double *in, double *diffusion_coeff);
static double diffusion_vector_multiply(double *a, double *b);
static void diffusion_matrix_multiply(double *in, double *out, double theta, double *diffusion_coeff);

static struct inexch_isotropic
{
  double in;
  double Kappa;
  double Volume;
  double Mass;
}
 *InExch_Isotropic;

void diffuse(double *data, double *diffusion_coeff, double dt)
{
  double *data_old, *data_old_lhs, *residual, *dvec, *qvec;
  double delta0, delta1, change, max_change, max_change_all, alpha, beta;
  double theta = 0.5;           /* Crank Nicolson */
  int i, iter;

  data_old = (double *) mymalloc("diffusion_data_old", NumGas * sizeof(double));
  data_old_lhs = (double *) mymalloc("diffusion_data_old_lhs", NumGas * sizeof(double));
  residual = (double *) mymalloc("diffusion_residual", NumGas * sizeof(double));
  dvec = (double *) mymalloc("diffusion_dvec", NumGas * sizeof(double));
  qvec = (double *) mymalloc("diffusion_dqec", NumGas * sizeof(double));

  for(i = 0; i < NumGas; i++)
    if(P[i].Type == 0 && P[i].ID != 0 && P[i].Mass > 0)
      {
        data_old[i] = data[i];
        diffusion_coeff[i] *= dt;
      }

  diffusion_matrix_multiply(data_old, data_old_lhs, theta - 1., diffusion_coeff);
  diffusion_matrix_multiply(data_old, residual, theta, diffusion_coeff);

  for(i = 0; i < NumGas; i++)
    if(P[i].Type == 0 && P[i].ID != 0 && P[i].Mass > 0)
      {
        residual[i] = data_old_lhs[i] - residual[i];
        dvec[i] = residual[i];
      }

  delta0 = diffusion_vector_multiply(residual, residual);
  delta1 = delta0;

  iter = 0;
  max_change_all = 1. + ACCURACY_DIFFUSION;
  //  while(iter < MAX_ITER_DIFFUSION && max_change_all > ACCURACY_DIFFUSION && delta1 > 0)
  while(iter < MAX_ITER_DIFFUSION && max_change_all > ACCURACY_DIFFUSION && delta1 > 1e-300)
    {
      diffusion_matrix_multiply(dvec, qvec, theta, diffusion_coeff);
      alpha = delta1 / diffusion_vector_multiply(dvec, qvec);

      max_change = 0;
      for(i = 0; i < NumGas; i++)
        {
          data[i] += alpha * dvec[i];
          residual[i] -= alpha * qvec[i];

          /* calculate maximum residual */
          change = fabs(residual[i] / data[i]);
          //change = fabs( alpha * dvec[i] / data[i] );
          if(change > max_change)
            max_change = change;
        }

      delta0 = delta1;
      delta1 = diffusion_vector_multiply(residual, residual);

      beta = delta1 / delta0;

      double energy_sum_all;
      double energy_sum = 0;
      for(i = 0; i < NumGas; i++)
        {
          dvec[i] = residual[i] + beta * dvec[i];
          if(P[i].Type == 0 && P[i].ID != 0 && P[i].Mass > 0)
            {
              if(data[i] < 0.)
                {
                  printf("Changing CR Energy in Cell ID=%d from %g erg/g to %g erg/g.\n", P[i].ID, data[i], SphP[i].Utherm * 1e-10);
                  data[i] = SphP[i].Utherm * 1e-10;
                }
              energy_sum += data[i] * SphP[i].Volume;
            }
        }

      MPI_Allreduce(&max_change, &max_change_all, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      MPI_Reduce(&energy_sum, &energy_sum_all, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

      mpi_printf("Diffusion iter=%d, delta1=%g, delta1/delta0=%g, maximum residual=%g.\n", iter, delta1, delta1 / delta0, max_change_all);
      mpi_printf("   Total CR Energy = %g erg.\n", energy_sum_all * All.UnitEnergy_in_cgs);

      iter++;
    }

  for(i = 0; i < NumGas; i++)
    if(P[i].Type == 0 && P[i].ID != 0 && P[i].Mass > 0)
      {
        diffusion_coeff[i] /= dt;
      }

  myfree(qvec);
  myfree(dvec);
  myfree(residual);
  myfree(data_old_lhs);
  myfree(data_old);
}

double diffusion_vector_multiply(double *a, double *b)
{
  int i;
  double sum, sumall;

  sum = 0;
  for(i = 0; i < NumGas; i++)
    if(P[i].Type == 0 && P[i].ID != 0 && P[i].Mass > 0)
      sum += a[i] * b[i];

  MPI_Allreduce(&sum, &sumall, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return sumall;
}

void diffusion_matrix_multiply(double *in, double *out, double theta, double *diffusion_coeff)
{
  int i;
  double in_i, in_j;
  double dx_ij, dy_ij, dz_ij, r_ij;
  double Kappa_i, Kappa_j;
  double kappa_mean_ij = 0.0;
  double area_ij = 0.0;

  InExch_Isotropic = (struct inexch_isotropic *) mymalloc("InExch_Isotropic", Mesh_nimport * sizeof(struct inexch_isotropic));

  diffusion_exchange_in_vector(in, diffusion_coeff);

  for(i = 0; i < NumGas; i++)
    {
      if(P[i].Type == 0 && P[i].ID != 0 && P[i].Mass > 0)
        {
          Kappa_i = diffusion_coeff[i];
          in_i = in[i];

          double out_sum = 0.0;
          double alpha = 0.0;
          int q = SphP[i].first_connection;

          while(q >= 0)
            {
              int dp = DC[q].dp_index;
              int vf = DC[q].vf_index;
              int particle = Mesh.DP[dp].index;

              if(particle < 0)
                {
                  q = DC[q].next;
                  continue;
                }

              if(particle >= NumGas && Mesh.DP[dp].task == ThisTask)
                particle -= NumGas;

              if(Mesh.DP[dp].task == ThisTask)
                {
                  Kappa_j = diffusion_coeff[particle];
                  in_j = in[particle];
                }
              else
                {
                  Kappa_j = InExch_Isotropic[particle].Kappa;
                  in_j = InExch_Isotropic[particle].in;
                }

              dx_ij = P[i].Pos[0] - Mesh.DP[dp].x;
              dy_ij = P[i].Pos[1] - Mesh.DP[dp].y;
              dz_ij = P[i].Pos[2] - Mesh.DP[dp].z;
              r_ij = sqrt(dx_ij * dx_ij + dy_ij * dy_ij + dz_ij * dz_ij);

              area_ij = Mesh.VF[vf].area;

              if((Kappa_i + Kappa_j) > 0)
                {
                  kappa_mean_ij = 2.0 * (Kappa_i * Kappa_j) / (Kappa_i + Kappa_j);

                  double w = (kappa_mean_ij * area_ij) / r_ij;
                  alpha += w;
                  out_sum += w * (in_j - in_i);
                }

              if(q == SphP[i].last_connection)
                break;

              q = DC[q].next;
            }

          out[i] = in[i] - theta * out_sum / SphP[i].Volume;
        }
    }                           /* for loop ends */

  myfree(InExch_Isotropic);
}

void diffusion_exchange_in_vector(double *in, double *diffusion_coeff)
{
  int listp;
  int j, p, task, off;
  int ngrp, recvTask, place;    /* there was earlier the following defined: sendTask */

  struct inexch_isotropic *tmpInExch;
  tmpInExch = (struct inexch_isotropic *) mymalloc("tmpInExch", Mesh_nexport * sizeof(struct inexch_isotropic));

  /* prepare data for export */
  for(j = 0; j < NTask; j++)
    Mesh_Send_count[j] = 0;

  for(p = 0; p < NumGas; p++)
    {
      if(P[p].Type == 0)
        {
          listp = List_P[p].firstexport;

          while(listp >= 0)
            {
              if((task = ListExports[listp].origin) != ThisTask)
                {
                  place = ListExports[listp].index;
                  off = Mesh_Send_offset[task] + Mesh_Send_count[task]++;

                  tmpInExch[off].in = in[place];
                  tmpInExch[off].Kappa = diffusion_coeff[place];
                  tmpInExch[off].Volume = SphP[place].Volume;
                  tmpInExch[off].Mass = P[place].Mass;
                }

              listp = ListExports[listp].nextexport;
            }
        }
    }

  /* exchange data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Mesh_Send_count[recvTask] > 0 || Mesh_Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&tmpInExch[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask] *
                           sizeof(struct inexch_isotropic), MPI_BYTE, recvTask, TAG_DENS_A, &InExch_Isotropic[Mesh_Recv_offset[recvTask]],
                           Mesh_Recv_count[recvTask] * sizeof(struct inexch_isotropic), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  myfree(tmpInExch);
}

#endif
