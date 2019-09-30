/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/MHD_non_ideal/implicit_ohmic_diffusion_new.c
 * \date        03/2016
 * \author      Federico Marinacci
 * \brief        
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 */

/*
 * ! \file implicit_ohmic_diffusion.c 
 *   \brief main driver for an isotropic diffusion solver
 * 
 */


#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <gsl/gsl_math.h>

#include "../allvars.h"
#include "../proto.h"
#include "../domain.h"
#include "../voronoi.h"

#ifdef IMPLICIT_OHMIC_DIFFUSION

#ifdef MHD_CT

#ifndef OHM_CRANK_NICHOLSON

#define MAXTETRA 200
#define EPSILON 1.0e-10
#define MAX_ITER 100

static struct inexch
{
  double Eta;
  double A;
  double SurfaceArea;
  double Pos[3];
} *InExch;


static MyFloat *Eta;
static MyFloat *A0;

static void do_joule_heating(double dt);

void ohmic_diffusion()
{
  mpi_printf("NON-IDEAL MHD: Doing ohmic diffusion for this time step\n");

  TIMER_START(CPU_OHM);

  double dt;

  Eta = (MyFloat *) mymalloc("Eta", NumGas * sizeof(MyFloat));
  A0 = (MyFloat *) mymalloc("A0", NumGas * sizeof(MyFloat));

  int numcycles;
  dt = get_time_step(&numcycles);
  mpi_printf("total time step = %g\tnumcycles = %d\n\n\n\n", dt, numcycles);

  mpi_printf("NON-IDEAL MHD: performing Joule heating \n");
  do_joule_heating(dt);
  mpi_printf("NON-IDEAL MHD: Joule heating done\n");

  InExch = (struct inexch *) mymalloc("InExch", Mesh_nimport * sizeof(struct inexch));

  for(int direction = 0; direction < 3; direction++)
    {
      mpi_printf("NON-IDEAL MHD: direction = %d\n", direction);

      initialize_eta(direction);

      diffusion_exchange_vector(direction);

      prepare_linear_iterations(dt, direction);

      mpi_printf("NON-IDEAL MHD: linear iterations done\n");
    }

  update_B_and_Energy();

  myfree(InExch);
  myfree(A0);
  myfree(Eta);

  TIMER_STOP(CPU_OHM);

  mpi_printf("NON-IDEAL MHD: Done\n");
}

void update_B_and_Energy()
{
  int i, k;
  struct pv_update_data pvd;
 
  if(All.ComovingIntegrationOn)
    {
      pvd.atime = All.Time;
      pvd.hubble_a = hubble_function(All.Time);
      pvd.a3inv = 1 / (All.Time * All.Time * All.Time);
    }
  else
    pvd.atime = pvd.hubble_a = pvd.a3inv = 1.0;

#ifdef ONLY_OHMIC_DIFFUSION
  for(i = 0; i < NumGas; i++)
    for(k = 0; k < DIMS; k++)
      SphP[i].AConserved[k] = SphP[i].A[k] * SphP[i].Volume;

  correct_ctr_b();
#else
  double *dEmag;
  dEmag = (double*) mymalloc("dEmag", sizeof(double) * NumGas);

  for(i = 0; i < NumGas; i++)
    {
      for(k = 0; k < DIMS; k++)
        SphP[i].AConserved[k] = SphP[i].A[k] * SphP[i].Volume;

      dEmag[i] = -0.5 * (SphP[i].BConserved[0] * SphP[i].BConserved[0] + SphP[i].BConserved[1] * SphP[i].BConserved[1] + SphP[i].BConserved[2] * SphP[i].BConserved[2]);
    }

  correct_ctr_b();

  for(i = 0; i < NumGas; i++)
    {
      dEmag[i] += 0.5 * (SphP[i].BConserved[0] * SphP[i].BConserved[0] + SphP[i].BConserved[1] * SphP[i].BConserved[1] + SphP[i].BConserved[2] * SphP[i].BConserved[2]);

      SphP[i].Energy += dEmag[i] / SphP[i].Volume * pvd.atime;

      update_internal_energy(P, SphP, i, &pvd);
      set_pressure_of_cell_internal(P, SphP, i);
    }

  myfree(dEmag);
#endif
}

double get_time_step(int *numcycles)
{
  double dt;
  dt = (All.ohmdiffusion_Ti_endstep - All.ohmdiffusion_Ti_begstep) * All.Timebase_interval;
  dt *= All.cf_atime / All.cf_time_hubble_a;

  *numcycles = 1;

  mpi_printf("dt = %le\n", dt);

  return dt;
}

void initialize_eta(int direction)
{
  int i;
  for(i = 0; i < NumGas; i++)
    {
      A0[i] = SphP[i].A[direction];
      Eta[i] = All.OhmicDiffusionCoefficient / All.cf_atime * All.cf_atime;
    }
}

void init_ohm_conductivity(void)
{
  mpi_printf("NON-IDEAL MHD: Conductivity set to %g (int units)\n", All.OhmicDiffusionCoefficient);
}


void prepare_linear_iterations(double dt, int direction)
{
  linear_iterations(dt, direction);
}

void set_coeff(HYPRE_IJMatrix * A, HYPRE_IJVector * x, HYPRE_IJVector * b, int *fld_offsets, double dt, int direction) 
{
  point *DP = Mesh.DP;
  face *VF = Mesh.VF;

  int i;
  int *cols = mymalloc("cols", sizeof(int) * MAXTETRA);
  double *vals = mymalloc("vals", sizeof(double) * MAXTETRA);

  for(i = 0; i < NumGas; i++)
    {
      cols[0] = fld_offsets[ThisTask] + i;
      vals[0] = 0.0;
      double sum = 0.0;
      double sfarea;
      double eta_q;
      int ind = 1;
      int q = SphP[i].first_connection;

      while(q >= 0)
        {
          int dp = DC[q].dp_index;
          int vf = DC[q].vf_index;
          int particle = DP[dp].index;

          if(particle < 0)
            {
              q = DC[q].next;
              continue;
            }

          if(DP[dp].task == ThisTask)
            {
              if(particle >= NumGas)
                particle -= NumGas;

              eta_q = Eta[particle];
              sfarea = SphP[particle].SurfaceArea;
            }
          else
            {
              eta_q = InExch[particle].Eta;
              sfarea = InExch[particle].SurfaceArea;
            }

          if(VF[vf].area < 1.0e-5 * SphP[i].SurfaceArea || VF[vf].area < 1.0e-5 * sfarea)
            {
              if(q == SphP[i].last_connection)
                break;

              q = DC[q].next;
              continue;
            }

          double p0[3], p1[3], dr, dx, dy, dz;

          int dp_other = VF[vf].p1;
          if(dp_other == dp)
            dp_other = VF[vf].p2;

          p0[0] = DP[dp].x;
          p0[1] = DP[dp].y;
          p0[2] = DP[dp].z;
          p1[0] = DP[dp_other].x;
          p1[1] = DP[dp_other].y;
          p1[2] = DP[dp_other].z;

          int originalindex;

          if(DP[dp].task == ThisTask)
            {
              if(particle < NumGas)
                originalindex = particle;
              else
                originalindex = particle - NumGas;
            }
          else
            originalindex = DP[dp].originalindex;

          dx = p1[0] - p0[0];
          dy = p1[1] - p0[1];
          dz = p1[2] - p0[2];

          dr = sqrt(dx * dx + dy * dy + dz * dz);

          double eta_ij = 0.5 * (Eta[i] + eta_q) / SphP[i].Volume; 
          double dval = -dt * eta_ij * VF[vf].area / dr;
	  sum += dval;

          cols[ind] = fld_offsets[DP[dp].task] + originalindex;
          vals[ind] = dval;
          ind++;

          if(ind >= MAXTETRA)
            terminate("No space left for storing matrix elements");

          if(q == SphP[i].last_connection)
            break;

          q = DC[q].next;
        }

      vals[0] = 1.0 - sum;
      int item = cols[0];
      HYPRE_IJMatrixSetValues(*A, 1, &ind, &item, cols, vals);
    }

  myfree(vals);
  myfree(cols);

  int *rows = (int *) mymalloc("rows", sizeof(int) * NumGas);
  double *xval = (double *) mymalloc("x", sizeof(double) * NumGas);
  double *bval = (double *) mymalloc("b", sizeof(double) * NumGas);

  for(i = 0; i < NumGas; i++)
    {
      rows[i] = fld_offsets[ThisTask] + i;
      xval[i] = SphP[i].A[direction];
      bval[i] = A0[i];
    }

  HYPRE_IJVectorSetValues(*b, NumGas, rows, bval);
  HYPRE_IJVectorSetValues(*x, NumGas, rows, xval);

  myfree(bval);
  myfree(xval);
  myfree(rows);
}

void linear_iterations(double dt, int direction)
{
  HYPRE_IJMatrix A;
  HYPRE_ParCSRMatrix parcsr_A;
  HYPRE_IJVector b;
  HYPRE_ParVector par_b;
  HYPRE_IJVector x;
  HYPRE_ParVector par_x;
  HYPRE_Solver solver;
  HYPRE_Solver precond;

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
  set_coeff(&A, &x, &b, fld_offsets, dt, direction);


  HYPRE_IJMatrixAssemble(A);
  HYPRE_IJMatrixGetObject(A, (void **) &parcsr_A);
  HYPRE_IJVectorAssemble(b);
  HYPRE_IJVectorGetObject(b, (void **) &par_b);
  HYPRE_IJVectorAssemble(x);
  HYPRE_IJVectorGetObject(x, (void **) &par_x);

  int num_iterations;
  double final_res_norm;

  /* Create solver */
  HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, &solver);

  /* Set Parameters */
  HYPRE_ParCSRGMRESSetTol(solver, EPSILON);
  HYPRE_ParCSRGMRESSetMaxIter(solver, MAX_ITER);
  HYPRE_ParCSRGMRESSetPrintLevel(solver, 4);
  HYPRE_ParCSRGMRESSetTol(solver, EPSILON); /* conv. tolerance */
  HYPRE_ParCSRGMRESSetMaxIter(solver, MAX_ITER);

  /* use BoomerAMG as preconditioner */
  HYPRE_BoomerAMGCreate(&precond);
  HYPRE_BoomerAMGSetCoarsenType(precond, 10);
  HYPRE_BoomerAMGSetStrongThreshold(precond, 0.25);
  HYPRE_BoomerAMGSetTol(precond, 0.0);
  HYPRE_BoomerAMGSetPrintLevel(precond, 1);
  HYPRE_BoomerAMGSetMaxIter(precond, 2);

  /* set the preconditioner */
  HYPRE_ParCSRGMRESSetPrecond(solver, HYPRE_BoomerAMGSolve, HYPRE_BoomerAMGSetup, precond);
  HYPRE_ParCSRGMRESSetup(solver, parcsr_A, par_b, par_x);

  HYPRE_ParCSRGMRESSolve(solver, parcsr_A, par_b, par_x);
  HYPRE_ParCSRGMRESGetNumIterations(solver, &num_iterations);
  HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);

  /* Destroy solver */
  HYPRE_ParCSRGMRESDestroy(solver);
  HYPRE_BoomerAMGDestroy(precond);
  get_xval(&x, fld_offsets, direction);
  HYPRE_IJMatrixDestroy(A);
  HYPRE_IJVectorDestroy(b);
  HYPRE_IJVectorDestroy(x);

  myfree(fld_offsets);

  if(final_res_norm > EPSILON)
    terminate("HYPRE failed to converge");

  mpi_printf("HYPRE: iter %d, res_norm %g\n", num_iterations, final_res_norm);
}

void get_xval(HYPRE_IJVector * x, int *fld_offsets, int direction)
{
  int *rows = (int *) mymalloc("rows", sizeof(int) * NumGas);
  double *XVec = (double *) mymalloc("XVec", sizeof(double) * NumGas);
  int i;

  for(i = 0; i < NumGas; i++)
    rows[i] = fld_offsets[ThisTask] + i;

  HYPRE_IJVectorGetValues(*x, NumGas, rows, XVec);

  //discard boundary cells from the solution
  for(i = 0; i < NumGas; i++)
    SphP[i].A[direction]= XVec[i];

  myfree(XVec);
  myfree(rows);
}

void diffusion_exchange_vector(int direction)
{
  int listp;
  int j, p, task, off;
  int ngrp, recvTask, place;    /* there was earlier the following defined: sendTask */
  struct inexch *tmpInExch;
  tmpInExch = (struct inexch *) mymalloc("tmpInExch", Mesh_nexport * sizeof(struct inexch));

  /* prepare data for export */
  for(j = 0; j < NTask; j++)
    Mesh_Send_count[j] = 0;

  for(j = 0; j < NumGasInMesh; j++)
    {
      p = List_InMesh[j];
      listp = List_P[p].firstexport;
      while(listp >= 0)
        {
          if((task = ListExports[listp].origin) != ThisTask)
            {
              place = ListExports[listp].index;
              off = Mesh_Send_offset[task] + Mesh_Send_count[task]++;

              tmpInExch[off].Eta = Eta[place];
              tmpInExch[off].A = SphP[place].A[direction];

              for(int k = 0; k < DIMS; k++)
                tmpInExch[off].Pos[k] = P[place].Pos[k];

              tmpInExch[off].SurfaceArea = SphP[place].SurfaceArea;
            }
          listp = ListExports[listp].nextexport;
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
              MPI_Sendrecv(&tmpInExch[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask] * sizeof(struct inexch), MPI_BYTE, recvTask, TAG_DENS_A, &InExch[Mesh_Recv_offset[recvTask]], Mesh_Recv_count[recvTask] * sizeof(struct inexch), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  myfree(tmpInExch);
}

void do_joule_heating(double dt)
{
#ifdef ONLY_OHMIC_DIFFUSION
  return;
#else
  struct pv_update_data pvd;

  if(All.ComovingIntegrationOn)
    {
      pvd.atime = All.Time;
    }
  else
    pvd.atime = 1.0;

  for(int i = 0; i < NumGas; i++)
    {
      double Jsq = SphP[i].CurlB[0] * SphP[i].CurlB[0] + SphP[i].CurlB[1] * SphP[i].CurlB[1] + SphP[i].CurlB[2] * SphP[i].CurlB[2];
      double du = dt * All.OhmicDiffusionCoefficient * Jsq / (SphP[i].Density * pvd.atime * pvd.atime);
      SphP[i].Utherm += du;
      SphP[i].Energy += pvd.atime * pvd.atime * du * P[i].Mass;
    }
#endif
}

#endif  /* OHM_CRANK_NICHOLSON */

#endif  /* MHD_CT */

#endif /* IMPLICIT_OHMIC_DIFFUSION */

