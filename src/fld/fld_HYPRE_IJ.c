/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/fld/fld_HYPRE.c
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
#if defined(FLD_HYPRE_IJ1) || defined(FLD_HYPRE_IJ2)

#include "_hypre_utilities.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"


#define MAX_ITER 100
#define ACCURACY 1.0e-8

#define IND_MAX 1000

#ifdef FLD_HYPRE_IJ1
void fld_set_coeff(HYPRE_IJMatrix * A, HYPRE_IJVector * x, HYPRE_IJVector * b, int *fld_offsets, double dt)
{

  int i;

  int *rows = (int *) mymalloc("rows", sizeof(int) * NumGas);
  double *xval = (double *) mymalloc("x", sizeof(double) * NumGas);
  double *bval = (double *) mymalloc("b", sizeof(double) * NumGas);

  for(i = 0; i < NumGas; i++)
    {
      rows[i] = fld_offsets[ThisTask] + i;
      xval[i] = 0.;

      bval[i] = SphP[i].b;
    }

  HYPRE_IJVectorSetValues(*b, NumGas, rows, bval);
  HYPRE_IJVectorSetValues(*x, NumGas, rows, xval);

  myfree(bval);
  myfree(xval);
  myfree(rows);

  int *cols = mymalloc("cols", sizeof(int) * IND_MAX);
  double *vals = mymalloc("vals", sizeof(double) * IND_MAX);

  for(i = 0; i < NumGas; i++)
    {
      cols[0] = fld_offsets[ThisTask] + i;
      vals[0] = SphP[i].w;

      int ind = 1;

      MyFloat volume = SphP[i].Volume;
      MyFloat Kappa_i = SphP[i].Kappa_diff;

      MyFloat x_i = P[i].Pos[0];
      MyFloat y_i = P[i].Pos[1];
      MyFloat z_i = P[i].Pos[2];

      MyFloat diag = 1.;

      int q = SphP[i].first_connection;

      while(q >= 0)
        {
          int dp = DC[q].dp_index;
          int vf = DC[q].vf_index;
          int other = Mesh.DP[dp].index;

          double area_ij = Mesh.VF[vf].area;

          if(other < 0)
            {
              if(q == SphP[i].last_connection)
                break;
              q = DC[q].next;
              continue;
            }

          if(other >= NumGas && Mesh.DP[dp].task == ThisTask)
            other -= NumGas;

          if(other == i)
            {
#ifdef FLD_TEST_BOUNDARY
              if(P[i].ID >= FLD_UPPER_BOUNDARY_MINID && P[i].ID <= FLD_UPPER_BOUNDARY_MAXID)
                {
                  diag += Kappa_i * area_ij * dt / volume / 0.5;
                }
#endif
              if(q == SphP[i].last_connection)
                break;
              q = DC[q].next;
              continue;
            }

          MyFloat Kappa_j;
          if(Mesh.DP[dp].task == ThisTask)
            {
              Kappa_j = SphP[other].Kappa_diff;
            }
          else
            {
              Kappa_j = PrimExch[other].Kappa_diff;
            }

          MyFloat x_j = Mesh.DP[dp].x;
          MyFloat y_j = Mesh.DP[dp].y;
          MyFloat z_j = Mesh.DP[dp].z;


          MyFloat dx_ij = nearest_x(x_i - x_j);
          MyFloat dy_ij = nearest_y(y_i - y_j);
          MyFloat dz_ij = nearest_z(z_i - z_j);



          MyFloat r_ij = sqrt(dx_ij * dx_ij + dy_ij * dy_ij + dz_ij * dz_ij);

          if((Kappa_i + Kappa_j) > 0)
            {
              double coeff_mean_ij = 2.0 * (Kappa_i * Kappa_j) / (Kappa_i + Kappa_j);
              double w = coeff_mean_ij * area_ij * dt / volume / r_ij;

              cols[ind] = fld_offsets[Mesh.DP[dp].task] + Mesh.DP[dp].originalindex;
              vals[ind] = -w;

              ind++;

              if(ind == IND_MAX)
                terminate("a cell with too many connections found");

              diag += w;
            }


          if(q == SphP[i].last_connection)
            break;

          q = DC[q].next;
        }

      vals[0] += diag;

      int item = fld_offsets[ThisTask] + i;
      HYPRE_IJMatrixSetValues(*A, 1, &ind, &item, cols, vals);
    }

  myfree(vals);
  myfree(cols);
}
#endif

#ifdef FLD_HYPRE_IJ2
void fld_set_coeff(HYPRE_IJMatrix * A, HYPRE_IJVector * x, HYPRE_IJVector * b, int *fld_offsets, double dt)
{

  int i;

  int *rows = (int *) mymalloc("rows", sizeof(int) * NumGas);
  double *xval = (double *) mymalloc("x", sizeof(double) * NumGas);
  double *bval = (double *) mymalloc("b", sizeof(double) * NumGas);

  for(i = 0; i < NumGas; i++)
    {
      rows[i] = fld_offsets[ThisTask] + i;
      xval[i] = 0.;

      bval[i] = SphP[i].b;
    }

  HYPRE_IJVectorSetValues(*b, NumGas, rows, bval);
  HYPRE_IJVectorSetValues(*x, NumGas, rows, xval);

  myfree(bval);
  myfree(xval);
  myfree(rows);

  int *cols = mymalloc("cols", sizeof(int) * IND_MAX);
  double *vals;

  for(i = 0; i < NumGas; i++)
    {
      cols[0] = fld_offsets[ThisTask] + i;


      vals = SphP[i].w;

      cols[0] = Mesh.DP[i].neighbors[0];
      cols[1] = Mesh.DP[i].neighbors[1];
      cols[2] = Mesh.DP[i].neighbors[2];
      cols[3] = Mesh.DP[i].neighbors[3];
      cols[4] = i;
      cols[5] = Mesh.DP[Mesh.DP[i].neighbors[0]].neighbors[2];
      cols[6] = Mesh.DP[Mesh.DP[i].neighbors[0]].neighbors[3];
      cols[7] = Mesh.DP[Mesh.DP[i].neighbors[1]].neighbors[2];
      cols[8] = Mesh.DP[Mesh.DP[i].neighbors[1]].neighbors[3];

      int ind = 9;

      int item = fld_offsets[ThisTask] + i;
      HYPRE_IJMatrixSetValues(*A, 1, &ind, &item, cols, vals);
    }

  myfree(cols);
}
#endif

void fld_get_xval(HYPRE_IJVector * x, double *XVec, int *fld_offsets)
{

  int *rows = (int *) mymalloc("rows", sizeof(int) * NumGas);

  int i;
  for(i = 0; i < NumGas; i++)
    {
      rows[i] = fld_offsets[ThisTask] + i;
    }

  HYPRE_IJVectorGetValues(*x, NumGas, rows, XVec);

  myfree(rows);

}

void fld_radtransfer(double *XVec, double dt)
{
  HYPRE_IJMatrix A;
  HYPRE_ParCSRMatrix parcsr_A;
  HYPRE_IJVector b;
  HYPRE_ParVector par_b;
  HYPRE_IJVector x;
  HYPRE_ParVector par_x;

  HYPRE_Solver solver;


  int *fld_offsets = (int *) mymalloc("fld_offsets", sizeof(int) * (NTask + 1));

  MPI_Allgather(&NumGas, 1, MPI_INT, fld_offsets, 1, MPI_INT, MPI_COMM_WORLD);

  int i;
  int cumm = 0;
  for(i = 0; i <= NTask; i++)
    {
      int tmp = cumm;
      cumm += fld_offsets[i];
      fld_offsets[i] = tmp;
    }

  HYPRE_IJMatrixCreate(MPI_COMM_WORLD, fld_offsets[ThisTask], fld_offsets[ThisTask + 1] - 1, fld_offsets[ThisTask], fld_offsets[ThisTask + 1] - 1, &A);
  HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
  HYPRE_IJMatrixInitialize(A);

  HYPRE_IJVectorCreate(MPI_COMM_WORLD, fld_offsets[ThisTask], fld_offsets[ThisTask + 1] - 1, &b);
  HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(b);

  HYPRE_IJVectorCreate(MPI_COMM_WORLD, fld_offsets[ThisTask], fld_offsets[ThisTask + 1] - 1, &x);
  HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(x);

  fld_set_coeff(&A, &x, &b, fld_offsets, dt);

  HYPRE_IJMatrixAssemble(A);
  HYPRE_IJMatrixGetObject(A, (void **) &parcsr_A);

  HYPRE_IJVectorAssemble(b);
  HYPRE_IJVectorGetObject(b, (void **) &par_b);

  HYPRE_IJVectorAssemble(x);
  HYPRE_IJVectorGetObject(x, (void **) &par_x);

  int num_iterations;
  double final_res_norm;

  /* Create solver */
  HYPRE_BoomerAMGCreate(&solver);

  /* Set some parameters (See Reference Manual for more parameters) */
  HYPRE_BoomerAMGSetPrintLevel(solver, 3);      /* print solve info + parameters */
  HYPRE_BoomerAMGSetCoarsenType(solver, 6);     /* Falgout coarsening */
  HYPRE_BoomerAMGSetRelaxType(solver, 3);       /* G-S/Jacobi hybrid relaxation */
  HYPRE_BoomerAMGSetNumSweeps(solver, 1);       /* Sweeeps on each level */
  HYPRE_BoomerAMGSetMaxLevels(solver, 20);      /* maximum number of levels */
  HYPRE_BoomerAMGSetTol(solver, ACCURACY);      /* conv. tolerance */
  HYPRE_BoomerAMGSetMaxIter(solver, MAX_ITER);

  /* Now setup and solve! */
  HYPRE_BoomerAMGSetup(solver, parcsr_A, par_b, par_x);
  HYPRE_BoomerAMGSolve(solver, parcsr_A, par_b, par_x);

  /* Run info - needed logging turned on */
  HYPRE_BoomerAMGGetNumIterations(solver, &num_iterations);
  HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &final_res_norm);

  /* Destroy solver */
  HYPRE_BoomerAMGDestroy(solver);

  fld_get_xval(&x, XVec, fld_offsets);

  HYPRE_IJMatrixDestroy(A);
  HYPRE_IJVectorDestroy(b);
  HYPRE_IJVectorDestroy(x);

  myfree(fld_offsets);


  if(final_res_norm > ACCURACY)
    terminate("HYPRE failed to converge");

  mpi_printf("HYPRE: iter %d, res_norm %g\n", num_iterations, final_res_norm);
}

#endif
#endif
