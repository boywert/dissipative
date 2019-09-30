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
#ifdef FLD_HYPRE

#include "_hypre_utilities.h"
#include "HYPRE_struct_ls.h"

#define MAX_ITER 100
#define ACCURACY 1.0e-8

#ifdef FLD_CONES
#define NCOEFF 9
#else
#define NCOEFF 5
#endif



void fld_set_matrix_coeff(double *coeff, double *b, int node, int size)
{
  int i;
  int j;
  int current = node;
  int start1;

  while(current > Ngb_MaxPart)
    {
      current = Ngb_Nodes[current].u.suns[0];
    }

  for(i = 0; i < size; i++)
    {
      start1 = current;
      for(j = 0; j < size; j++)
        {
          int k;

          for(k = 0; k < NCOEFF; k++)
            {
              coeff[NCOEFF * (i * size + j) + k] = SphP[current].w[k];
              if(!gsl_finite(SphP[current].w[k]))
                {
                  terminate("bad matrix coeff: %d %d %d %d %g", current, k, i, j, SphP[current].w[k]);
                }
            }

          b[i * size + j] = SphP[current].b;

          current = Mesh.DP[current].neighbors[1];
        }
      current = Mesh.DP[start1].neighbors[3];
    }
}

void fld_get_xval(double *xtmp, double *XVec, int node, int size)
{
  int i;
  int j;
  int current = node;
  int start1;

  while(current > Ngb_MaxPart)
    {
      current = Ngb_Nodes[current].u.suns[0];
    }

  for(i = 0; i < size; i++)
    {
      start1 = current;
      for(j = 0; j < size; j++)
        {

          XVec[current] = xtmp[i * size + j];
          current = Mesh.DP[current].neighbors[1];
        }
      current = Mesh.DP[start1].neighbors[3];
    }

}

void fld_radtransfer(double *XVec, double dt)
{
  int sizeX = (1 << Mesh.maxlevel) * boxSize_X / DomainLen;
  int sizeY = (1 << Mesh.maxlevel) * boxSize_Y / DomainLen;
  //printf("sizeX %d, sizeY %d\n", sizeX, sizeY);

  int ilower[2], iupper[2];
  int sizeBox[2] = { sizeX, sizeY };

  ilower[0] = 0;
  ilower[1] = 0;
  iupper[0] = sizeBox[0] - 1;
  iupper[1] = sizeBox[1] - 1;

  HYPRE_StructGrid grid;
  HYPRE_StructStencil stencil;
  HYPRE_StructMatrix A;
  HYPRE_StructVector b;
  HYPRE_StructVector x;
  HYPRE_StructSolver solver;

  int num_iterations;
  double final_res_norm;

  HYPRE_StructGridCreate(MPI_COMM_WORLD, 2, &grid);

  int i;
  for(i = 0; i < NTopleaves; i++)
    {
      if(DomainTask[i] == ThisTask)
        {
          int node = Ngb_DomainNodeIndex[i];

          if(Ngb_Nodes[node].u.suns[0] == -1)
            continue;

          int size = 1 << (Mesh.maxlevel - Ngb_Nodes[node].level);

          ilower[0] = (Ngb_Nodes[node].Center[0] - amr_length[Ngb_Nodes[node].level] / 2.) / amr_length[Mesh.maxlevel];
          ilower[1] = (Ngb_Nodes[node].Center[1] - amr_length[Ngb_Nodes[node].level] / 2.) / amr_length[Mesh.maxlevel];
          iupper[0] = ilower[0] + size - 1;
          iupper[1] = ilower[1] + size - 1;

          HYPRE_StructGridSetExtents(grid, ilower, iupper);
        }
    }

  int periodicity[2];
#ifdef REFLECTIVE_X
  periodicity[0] = 0;
#else
  periodicity[0] = sizeBox[0];
#endif
#ifdef REFLECTIVE_Y
  periodicity[1] = 0;
#else
  periodicity[1] = sizeBox[1];
#endif
  HYPRE_StructGridSetPeriodic(grid, periodicity);
  HYPRE_StructGridAssemble(grid);


  HYPRE_StructStencilCreate(2, NCOEFF, &stencil);
  int entry;
  int offsets[13][2] = { {-1, 0}, {1, 0}, {0, -1}, {0, 1}, {0, 0}, {-1, -1}, {-1, 1}, {1, -1}, {1, 1}, {-2, 0}, {2, 0}, {0, -2}, {0, 2} };
  for(entry = 0; entry < NCOEFF; entry++)
    {
      HYPRE_StructStencilSetElement(stencil, entry, offsets[entry]);
    }


  HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, stencil, &A);
  HYPRE_StructMatrixInitialize(A);

  HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &b);
  HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &x);
  HYPRE_StructVectorInitialize(b);
  HYPRE_StructVectorInitialize(x);

  int stencil_indices[NCOEFF];
  int j;
  for(j = 0; j < NCOEFF; j++)
    {
      stencil_indices[j] = j;
    }

  for(i = 0; i < NTopleaves; i++)
    {
      if(DomainTask[i] == ThisTask)
        {
          int node = Ngb_DomainNodeIndex[i];

          if(Ngb_Nodes[node].u.suns[0] == -1)
            continue;

          int size = 1 << (Mesh.maxlevel - Ngb_Nodes[node].level);

          ilower[0] = (Ngb_Nodes[node].Center[0] - amr_length[Ngb_Nodes[node].level] / 2.) / amr_length[Mesh.maxlevel];
          ilower[1] = (Ngb_Nodes[node].Center[1] - amr_length[Ngb_Nodes[node].level] / 2.) / amr_length[Mesh.maxlevel];
          iupper[0] = ilower[0] + size - 1;
          iupper[1] = ilower[1] + size - 1;

          double *xtmp = mymalloc("xtmp", size * size * sizeof(double));
          double *btmp = mymalloc("btmp", size * size * sizeof(double));
          double *coeff = mymalloc("coeff", NCOEFF * size * size * sizeof(double));

          memset(xtmp, 0, size * size * sizeof(double));
          fld_set_matrix_coeff(coeff, btmp, node, size);

          HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, NCOEFF, stencil_indices, coeff);
          HYPRE_StructVectorSetBoxValues(b, ilower, iupper, btmp);
          HYPRE_StructVectorSetBoxValues(x, ilower, iupper, xtmp);

          myfree(coeff);
          myfree(btmp);
          myfree(xtmp);
        }
    }

  HYPRE_StructMatrixAssemble(A);
  HYPRE_StructVectorAssemble(b);
  HYPRE_StructVectorAssemble(x);

#ifdef HYPRE_PCG
  HYPRE_StructPCGCreate(MPI_COMM_WORLD, &solver);
  HYPRE_StructPCGSetMaxIter(solver, MAX_ITER);
  HYPRE_StructPCGSetTol(solver, ACCURACY);
  HYPRE_StructPCGSetTwoNorm(solver, 1);
  HYPRE_StructPCGSetRelChange(solver, 0);
  HYPRE_StructPCGSetPrintLevel(solver, 1);
  HYPRE_StructPCGSetLogging(solver, 2);


  HYPRE_StructSolver precond;
  HYPRE_StructSMGCreate(MPI_COMM_WORLD, &precond);
  HYPRE_StructSMGSetMemoryUse(precond, 0);
  HYPRE_StructSMGSetMaxIter(precond, 1);
  HYPRE_StructSMGSetTol(precond, 0.0);
  HYPRE_StructSMGSetZeroGuess(precond);
  HYPRE_StructSMGSetNumPreRelax(precond, 1);
  HYPRE_StructSMGSetNumPostRelax(precond, 1);
  HYPRE_StructPCGSetPrecond(solver, HYPRE_StructSMGSolve, HYPRE_StructSMGSetup, precond);

  HYPRE_StructPCGSetup(solver, A, b, x);
  HYPRE_StructPCGSolve(solver, A, b, x);


  HYPRE_StructPCGGetNumIterations(solver, &num_iterations);
  HYPRE_StructPCGGetFinalRelativeResidualNorm(solver, &final_res_norm);


  HYPRE_StructPCGDestroy(solver);

#else

  HYPRE_StructSMGCreate(MPI_COMM_WORLD, &solver);
  HYPRE_StructSMGSetMemoryUse(solver, 0);
  HYPRE_StructSMGSetMaxIter(solver, MAX_ITER);
  HYPRE_StructSMGSetTol(solver, ACCURACY);
  HYPRE_StructSMGSetRelChange(solver, 0);
  HYPRE_StructSMGSetNumPreRelax(solver, 1);
  HYPRE_StructSMGSetNumPostRelax(solver, 1);
  HYPRE_StructSMGSetLogging(solver, 1);

  HYPRE_StructSMGSetup(solver, A, b, x);
  HYPRE_StructSMGSolve(solver, A, b, x);

  HYPRE_StructSMGGetNumIterations(solver, &num_iterations);
  HYPRE_StructSMGGetFinalRelativeResidualNorm(solver, &final_res_norm);

  HYPRE_StructSMGDestroy(solver);
#endif

  for(i = 0; i < NTopleaves; i++)
    {
      if(DomainTask[i] == ThisTask)
        {
          int node = Ngb_DomainNodeIndex[i];

          if(Ngb_Nodes[node].u.suns[0] == -1)
            continue;

          int size = 1 << (Mesh.maxlevel - Ngb_Nodes[node].level);

          ilower[0] = (Ngb_Nodes[node].Center[0] - amr_length[Ngb_Nodes[node].level] / 2.) / amr_length[Mesh.maxlevel];
          ilower[1] = (Ngb_Nodes[node].Center[1] - amr_length[Ngb_Nodes[node].level] / 2.) / amr_length[Mesh.maxlevel];
          iupper[0] = ilower[0] + size - 1;
          iupper[1] = ilower[1] + size - 1;

          double *xtmp = mymalloc("xtmp", size * size * sizeof(double));

          HYPRE_StructVectorGetBoxValues(x, ilower, iupper, xtmp);
          fld_get_xval(xtmp, XVec, node, size);

          myfree(xtmp);
        }
    }

  HYPRE_StructGridDestroy(grid);
  HYPRE_StructStencilDestroy(stencil);
  HYPRE_StructMatrixDestroy(A);
  HYPRE_StructVectorDestroy(b);
  HYPRE_StructVectorDestroy(x);

  if(final_res_norm > ACCURACY)
    terminate("HYPRE failed to converge: iter %d, res_norm %g\n", num_iterations, final_res_norm);

  mpi_printf("HYPRE: iter %d, res_norm %g\n", num_iterations, final_res_norm);
}

#endif
#endif
