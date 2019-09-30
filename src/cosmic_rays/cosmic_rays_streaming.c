/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/cosmic_rays_streaming.c
 * \date        02/2016
 * \author      R. Pakmor
 * \brief       (anisotropic) Cosmic Ray streaming
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 */

#include "../allvars.h"
#include "../proto.h"
#include "string.h"
#include "sys/stat.h"

#ifdef COSMIC_RAYS_STREAMING

#include <gsl/gsl_linalg.h>

#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "cosmic_rays.h"

#define STREAMING_MAX_ITER 100
#define STREAMING_MAX_ITER_NONLINEAR 500
#define STREAMING_ACCURACY 1.0e-6
#define STREAMING_MAX_COLS 1000

/* structs first ... */
static struct diff_face_data {
  int cornerFirst;
  int cornerCount;
  double bfld[3];
  double btot;
  double rho;
  double chi;
  double grad[3];

  int active;
  double dt;
  double area;
  double tanh;
  double vel_alfven;
  double gradb;
  double crEnergyDensity;
  
  double nx, ny, nz;
  double mx, my, mz;
  double px, py, pz;

  double failWeight;
} *diff_face_data;

static struct corner_list {
  int index;
  double weight;
} *corner_list;

static struct corner_data {
  int active;
  double matrix[NUMDIMS+1][NUMDIMS+1];
  
  double CREnergyDensity;
  double CREnergyDensityGrad[3];
  double bfld[3];
  double rho;
  double chi;

  int tetra;
  int fail;
} *corner_data;

static struct crflux_list_data {
  int task;
  int index;
  double dCREnergy;
} *CRFluxList;

#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
static double *refCenter;
static int *pIsReflective;
#endif

static int cornerCount;
static void add_corner( int tetra, int iface );
static int get_timebin_of_point( int point );
static double get_surface_area_of_cell( int point );

static void report_total_cr_energy();
static void prepare_stuff();
static void compute_least_squares_matrix_at_corners();
static void compute_geometry_of_interface( int iface );
static void compute_quantities_at_interface( int iface );
static void compute_quantities_at_corners();
static int crflux_list_data_compare(const void *a, const void *b);
static void streaming_explicit();
static void streaming_implicit();
static void streaming_implicit_extract_solution();
static void free_stuff();

struct column_data {
  int index;
  double value;
  int next_column;
  int task;
};

struct matrix_data {
  int local_row_count;
  int *local_row_index;  /* row index of a particle */
  int *local_part_index; /* particle index of a row */
  int *first_column;
  int *last_column;
  int *offsets;
  int *imported_particle_indizes;
  
  struct column_data *column_data;
  int cd_count;
  
  int *external_index; /* stores global index of external particles imported to PrimExch */
};

struct hypre_data {
  int nRows;
  int *columnCounts;
  int *rows;
  int *cols;
  double *vals;
  double *bval;
  double *xval;
  int *tasks;
  int initialized;
};

struct solver_data {
  double *residuals;
  double *initial_cr_energy_density;
  double *last_cr_energy_density;
  double *current_cr_energy_density;
  double *delta_cr_energy_density;
  double norm_initial_cr_energy_density;
  double norm_initial_residual;
  double norm_last_residual;
  double *g;
};

struct external_element {
  int task;
  int index;
  int column_index;
  int originaltask;
  double value;
};

struct external_elements {
  int N_external_elements;
  int Nmax_external_elements;
  struct external_element *elements;
};

struct cr_energy_log {
  double sum;
  double min, max;
};

double *CR_EnergyExch;

static void exchange_cr_energy();
static void allocate_rows_offsets_and_indizes( struct matrix_data *md );
static void compute_rows_offsets_and_indizes( struct matrix_data *md );
static void free_hypre_data( struct hypre_data *hd );
static void free_rows_offsets_and_indizes( struct matrix_data *md );
static void set_matrix_coefficients( struct matrix_data *md, struct hypre_data *hd, struct solver_data *sd );
static void point_get_center( int p, double *Center );
static void add_matrix_element( int point_row, int point_column, double val, struct matrix_data *md, struct external_elements *ee );
static void exchange_row_indizes( struct matrix_data *md );
static void apply_implicit_fluxes( HYPRE_IJVector* x, struct matrix_data *md );
static int streaming_implicit_solve_matrix( struct matrix_data *md, struct hypre_data *hd, int useMultigridPreconditionier, MPI_Comm diffComm, double *res_norm, struct solver_data *sd );
static void streaming_implicit_compute_energy_log( struct cr_energy_log *lg, struct matrix_data *md );
static void compute_and_check_residuals( struct matrix_data *md, struct hypre_data *hd, MPI_Comm diffComm, int *success );
static void streaming_implicit_init_solver( struct matrix_data *md, struct hypre_data *hd, struct solver_data *sd );
static void streaming_implicit_free_solver( struct solver_data *sd );
static void streaming_implicit_compute_residuals( struct matrix_data *md, struct hypre_data *hd, struct solver_data *sd );
static int streaming_implicit_update_solution( struct matrix_data *md, struct hypre_data *hd, struct solver_data *sd, int iter, double lin_res_norm );
static double compute_abs( struct matrix_data *md, struct hypre_data *hd, double * val );
static double compute_error( struct matrix_data *md, struct hypre_data *hd, struct solver_data *sd );

// add #ifdef COSMIC_RAYS_STREAMING_CONSTANT_CHI

void cosmic_rays_do_streaming(void)
{
  TIMER_START(CPU_CR_STREAMING);
  
  int CountAll;
  MPI_Allreduce( &TimeBinsHydro.NActiveParticles, &CountAll, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
  if(CountAll == 0)
    {
      mpi_printf( "COSMIC_RAYS STREAMING: No gas cells active, skipping.\n" );
      return;
    }
  else
    mpi_printf( "COSMIC_RAYS STREAMING: Doing CR streaming for %d active cells.\n", CountAll );
  
  update_primitive_variables();
  exchange_primitive_variables();
  report_total_cr_energy();
  
  prepare_stuff();

#ifdef COSMIC_RAYS_STREAMING_EXPLICIT
  streaming_explicit();
#else
  streaming_implicit();
#endif
  
  free_stuff();
  update_primitive_variables();
  report_total_cr_energy();
  
  TIMER_STOP(CPU_CR_STREAMING);
}

void exchange_cr_energy()
{
  int listp;
  int i, j, p, task, off;
  int ngrp, recvTask, place;
  
  struct crDataExch {
    double CR_Energy;
    double CR_Chi;
  } *tmpDataExch, *tmpDataRecv;

  tmpDataExch = (struct crDataExch*) mymalloc("tmpDataExch", Mesh_nexport * sizeof(struct crDataExch));

  /* prepare data for export */
  for(j = 0; j < NTask; j++)
    Mesh_Send_count[j] = 0;

  for(i = 0; i < NumGasInMesh; i++)
    {
      p = List_InMesh[i];

      listp = List_P[p].firstexport;
      while(listp >= 0)
        {
          if((task = ListExports[listp].origin) != ThisTask)
            {
              place = ListExports[listp].index;
              off = Mesh_Send_offset[task] + Mesh_Send_count[task]++;
              
              tmpDataExch[off].CR_Energy = SphP[place].CR_Energy;
              tmpDataExch[off].CR_Chi = SphP[place].CR_Chi;
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
              tmpDataRecv = (struct crDataExch*) mymalloc("tmpDataRecv", Mesh_Recv_count[recvTask] * sizeof(struct crDataExch));

              /* get the values */
              MPI_Sendrecv(&tmpDataExch[Mesh_Send_offset[recvTask]],
                           Mesh_Send_count[recvTask] * sizeof(struct crDataExch), MPI_BYTE,
                           recvTask, TAG_DENS_A, tmpDataRecv, Mesh_Recv_count[recvTask] * sizeof(struct crDataExch), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

              for(i = 0; i < Mesh_Recv_count[recvTask]; i++)
                {
                  CR_EnergyExch[Mesh_Recv_offset[recvTask] + i] = tmpDataRecv[i].CR_Energy;
                  PrimExch[Mesh_Recv_offset[recvTask] + i].CR_Chi = tmpDataRecv[i].CR_Chi;
                }

              myfree(tmpDataRecv);
            }
        }
    }

  myfree(tmpDataExch);
}

void report_total_cr_energy()
{
  if(All.HighestActiveTimeBin != All.HighestOccupiedTimeBin)
    return;

  double CR_Energy, CR_EnergyAll;
  double CR_EDensMin, CR_EDensMax;
  double CR_EDensMinAll, CR_EDensMaxAll;
  
  CR_Energy = 0;
  CR_EDensMin = +MAX_DOUBLE_NUMBER;
  CR_EDensMax = -MAX_DOUBLE_NUMBER;

  int cell;
  for(cell = 0; cell < NumGas; cell++)
    {
      CR_Energy += SphP[cell].CR_Energy;
      double CR_EnergyDensity = SphP[cell].CR_Energy / SphP[cell].Volume;
      if(CR_EnergyDensity < CR_EDensMin)
        CR_EDensMin = CR_EnergyDensity;
      if(CR_EnergyDensity > CR_EDensMax)
        CR_EDensMax = CR_EnergyDensity;
    }
  
  MPI_Reduce( &CR_Energy, &CR_EnergyAll, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
  MPI_Reduce( &CR_EDensMin, &CR_EDensMinAll, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
  MPI_Reduce( &CR_EDensMax, &CR_EDensMaxAll, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
  mpi_printf( "COSMIC_RAYS: Total cosmic ray energy: %g, CR energy density min=%e, max=%e\n", CR_EnergyAll, CR_EDensMinAll, CR_EDensMaxAll );
}

void prepare_stuff()
{
#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
  refCenter = (double *)mymalloc( "refCenter", Mesh.Ndp * 3 * sizeof(double) );
  pIsReflective = (int*)mymalloc( "pIsReflective", Mesh.Ndp * sizeof(int) );
  memset( pIsReflective, 0, Mesh.Ndp*sizeof(int) );
  // 0: no interesting boundary, 1: reflective boundary, 2: outflow boundary
  
  {
    // we go through all the interfaces, look for boundary points and calculate their mirrored center of mass
    int iface;
    for(iface = 0; iface < Mesh.Nvf; iface++)
      {
        int p1 = Mesh.VF[iface].p1;
        int p2 = Mesh.VF[iface].p2;
        
        if(Mesh.DP[p1].task == ThisTask && Mesh.DP[p2].task == ThisTask)
          {
            int part1 = Mesh.DP[p1].index;
            int part2 = Mesh.DP[p2].index;
            
            if((part1 >= NumGas || part2 >= NumGas) && (part1 < NumGas || part2 < NumGas))
              {
                // we found a boundary point
                if(part1 >= NumGas) part1 -= NumGas;
                if(part2 >= NumGas) part2 -= NumGas;
                
                if(Mesh.VF[iface].area < 1e-5 * dmin(SphP[part1].SurfaceArea,SphP[part2].SurfaceArea))
                  continue;

                if(P[part1].ID == P[part2].ID)
                  {
                    // it is a reflective of outflow boundary, i.e. not a periodic boundary
                    int pCell, pGhost;
                    
                    if(Mesh.DP[p1].index >= NumGas)
                      {
                        pCell = p2;
                        pGhost = p1;
                      }
                    else if (Mesh.DP[p2].index >= NumGas)
                      {
                        pCell = p1;
                        pGhost = p2;
                      }
                    else
                      {
                        // just to get rid of the warning
                        pCell = pGhost = -1;
                      }
                    
                    int particle = Mesh.DP[pCell].index;

                    // calculate normal vector of the interface
                    double nx = Mesh.DP[pGhost].x - P[particle].Pos[0];
                    double ny = Mesh.DP[pGhost].y - P[particle].Pos[1];
                    double nz = Mesh.DP[pGhost].z - P[particle].Pos[2];
                    
                    // perpendicular on the surface
                    double nn = sqrt(nx*nx + ny*ny + nz*nz);
                    nx /= nn;
                    ny /= nn;
                    nz /= nn;
                    double fx = (SphP[particle].Center[0] - Mesh.VF[iface].cx);
                    double fy = (SphP[particle].Center[1] - Mesh.VF[iface].cy);
                    double fz = (SphP[particle].Center[2] - Mesh.VF[iface].cz);
                    double ff = (fx * nx + fy * ny + fz * nz);
                    
                    double px = SphP[particle].Center[0] - ff * nx;
                    double py = SphP[particle].Center[1] - ff * ny;
                    double pz = SphP[particle].Center[2] - ff * nz;
                    
                    refCenter[pGhost*3+0] = 2. * px - SphP[particle].Center[0];
                    refCenter[pGhost*3+1] = 2. * py - SphP[particle].Center[1];
                    refCenter[pGhost*3+2] = 2. * pz - SphP[particle].Center[2];
                    
                    int image_flags = Mesh.DP[pGhost].image_flags;
#if defined(REFLECTIVE_X)
                    if((image_flags & REFL_X_FLAGS) && (image_flags & OUTFLOW_X))
                      pIsReflective[pGhost] = 2;
                    else
                      pIsReflective[pGhost] = 1;
#endif
#if defined(REFLECTIVE_Y)
                    if((image_flags & REFL_Y_FLAGS) && (image_flags & OUTFLOW_Y))
                      pIsReflective[pGhost] = 2;
                    else
                      pIsReflective[pGhost] = 1;
#endif
#if defined(REFLECTIVE_Z)
                    if((image_flags & REFL_Z_FLAGS) && (image_flags & OUTFLOW_Z))
                      pIsReflective[pGhost] = 2;
                    else
                      pIsReflective[pGhost] = 1;
#endif
                  }
              }
          }
      }
  }
#endif
  
  /* import CR Energy for external cells */
  CR_EnergyExch = (double*) mymalloc( "CR_EnergyExch", sizeof(double) * Mesh_nimport );
  
  /* for all interfaces, make a list of all cells sharing a corner with the interface 
    and save their center positions */
  
  diff_face_data = (struct diff_face_data*)mymalloc( "diff_face_data", Mesh.Nvf*sizeof(struct diff_face_data) );  
  corner_list = (struct corner_list*)mymalloc( "corner_list", Mesh.Nvf*20*sizeof(struct corner_list) );
  corner_data = (struct corner_data*)mymalloc( "corner_data", Mesh.Ndt*sizeof(struct corner_data) );
  
  /* compute matrices and magnetic fields at corners */
  int icorner;
  for(icorner = 0; icorner < Mesh.Ndt; icorner++)
    corner_data[icorner].active = 0;
  
  /* flag interfaces that are needed and compute their timestep */
  int iface;
  for(iface = 0; iface < Mesh.Nvf; iface++)
    {
      struct diff_face_data *fd = &diff_face_data[iface];
      fd->active = 0;
      
      int p1 = Mesh.VF[iface].p1;
      int p2 = Mesh.VF[iface].p2;
      
      if(p1 < 0 || p2 < 0)
        continue;
      
      int timebin1 = get_timebin_of_point( p1 );
      int timebin2 = get_timebin_of_point( p2 );
      
      if(!TimeBinSynchronized[timebin1] && !TimeBinSynchronized[timebin2])
        continue;
      
      if(Mesh.DP[p1].ID < Mesh.DP[p2].ID)
        {
          if(TimeBinSynchronized[timebin1])
            {
              /* lower ID is active, its task is responsible */
              if(Mesh.DP[p1].task != ThisTask || Mesh.DP[p1].index >= NumGas)
                continue;
            }
          else
            {
              /* only higher ID is active, its task is responsible */
              if(Mesh.DP[p2].task != ThisTask || Mesh.DP[p2].index >= NumGas)
                continue;
            }
        }
      else if (Mesh.DP[p1].ID > Mesh.DP[p2].ID)
        {
          if(TimeBinSynchronized[timebin2])
            {
              /* lower ID is active, its task is responsible */
              if(Mesh.DP[p2].task != ThisTask || Mesh.DP[p2].index >= NumGas)
                continue;
            }
          else
            {
              /* only higher ID is active, its task is responsible */
              if(Mesh.DP[p1].task != ThisTask || Mesh.DP[p1].index >= NumGas)
                continue;
            }
        }
      else
        {
          /* interface with at least one ghost point */
          
          if(Mesh.DP[p1].task != ThisTask && Mesh.DP[p2].task != ThisTask)
            continue;
          
          if(Mesh.DP[p1].task != Mesh.DP[p2].task)
            terminate( "This should not happen, I think..." );
          
          if(Mesh.DP[p1].index >= NumGas && Mesh.DP[p2].index >= NumGas)
            continue;
          
          int p;
          if(Mesh.DP[p1].index < NumGas)
            p = p1;
          else
            p = p2;
          
          if(!TimeBinSynchronized[P[Mesh.DP[p].index].TimeBinHydro])
            continue;
        }
      
      double surfacearea = dmax( get_surface_area_of_cell(p1), get_surface_area_of_cell(p2) );
      if(Mesh.VF[iface].area <= 1e-5 * surfacearea)
        continue;
      
      /* if we make it here, the interface needs to be included and we are responsible for it */
      fd->active = 1;
      
      /* now lets get its timestep */
      int timeBin = imin( timebin1, timebin2 );
      if(timeBin == 0)
        fd->dt = 0;
      else
        fd->dt = (((integertime) 1) << timeBin) * All.Timebase_interval / All.cf_hubble_a;
    }

  cornerCount = 0;
  
#ifndef ONEDIMS
#ifdef TWODIMS
  const int edge_start[3] = { 1, 2, 0 };
  const int edge_end[3] = { 2, 0, 1 };
#else
  const int edge_start[6] = { 0, 0, 0, 1, 1, 2 };
  const int edge_end[6] = { 1, 2, 3, 2, 3, 3 };
  const int edge_opposite[6] = { 3, 1, 2, 3, 0, 1 };
  const int edge_nexttetra[6] = { 2, 3, 1, 0, 2, 0 };
#endif
#endif
  
  for(iface = 0; iface < Mesh.Nvf; iface++)
    {
      struct diff_face_data *fd = &diff_face_data[iface];
      int p1 = Mesh.VF[iface].p1;
      int p2 = Mesh.VF[iface].p2;
      
      if(p1 < 0 || p2 < 0 || fd->active == 0)
        {
          fd->area = 0;
          fd->active = 0;
          continue;
        }

      fd->area = Mesh.VF[iface].area;
      fd->cornerFirst = -1;
      fd->cornerCount = 0;

#ifdef ONEDIMS
      terminate( "not implemented at the moment." );
#else
#ifdef TWODIMS
      int tt = Mesh.VF[iface].dt_index;
      tetra *t = &Mesh.DT[tt];
      
      int nr;
      for(nr = 0; nr < 3; nr++)
        {
          int start_index = t->p[edge_start[nr]];
          int end_index = t->p[edge_end[nr]];
          
          if((start_index == p1 && end_index == p2) || (start_index == p2 && end_index == p1))
            break;
        }
      
      int qq = t->t[nr];
      
      fd->cornerFirst = cornerCount;

      add_corner( tt, iface );
      add_corner( qq, iface );

      for(icorner = fd->cornerFirst; icorner < fd->cornerFirst + fd->cornerCount; icorner++)
        corner_list[icorner].weight = 0.5;
#else /* 3D */
      int tt = Mesh.VF[iface].dt_index;
      tetra *t = &Mesh.DT[tt];
      
      fd->cornerFirst = cornerCount;

      int nr;
      for(nr = 0; nr < 6; nr++)
        {
          int start_index = t->p[edge_start[nr]];
          int end_index = t->p[edge_end[nr]];
          
          if((start_index == p1 && end_index == p2) || (start_index == p2 && end_index == p1))
            break;
        }
      
      tetra *prev, *next;      
      int i, j, k, l, m, ii, jj, kk, ll, nn;
      
      i = edge_start[nr];
      j = edge_end[nr];
      k = edge_opposite[nr];
      l = edge_nexttetra[nr];

      prev = t;

      do
        {
          nn = prev->t[l];
          next = &Mesh.DT[nn];

          add_corner( nn, iface );

          for(m = 0, ll = ii = jj = -1; m < 4; m++)
            {
              if(next->p[m] == prev->p[k])
                ll = m;
              if(next->p[m] == prev->p[i])
                ii = m;
              if(next->p[m] == prev->p[j])
                jj = m;
            }

          if(ll < 0 || ii < 0 || jj < 0)
            terminate("inconsistency");

          kk = 6 - (ll + ii + jj);

          prev = next;

          i = ii;
          l = ll;
          j = jj;
          k = kk;
        }
      while(next != t);
      
      /* fix weights */
#endif
#endif
    }

  compute_least_squares_matrix_at_corners();
  
  for(iface = 0; iface < Mesh.Nvf; iface++)
    {
      compute_geometry_of_interface( iface );

      compute_quantities_at_interface( iface );
    }
  
  MPI_Barrier( MPI_COMM_WORLD );
  mpi_printf( "COSMIC_RAYS: STREAMING: Preparations done.\n" );  
}

int get_timebin_of_point( int point )
{
  int particle = Mesh.DP[point].index;
  
  if(Mesh.DP[point].task == ThisTask)
    {
      if(particle >= NumGas)
        particle -= NumGas;
      return P[particle].TimeBinHydro;
    }
  else
    return PrimExch[particle].TimeBinHydro;
}

double get_surface_area_of_cell( int point )
{
  int particle = Mesh.DP[point].index;
  
  if(Mesh.DP[point].task == ThisTask)
    {
      if(particle >= NumGas)
        particle -= NumGas;
      return SphP[particle].SurfaceArea;
    }
  else
    return PrimExch[particle].SurfaceArea;
}

void add_corner( int tt, int iface )
{
  if(cornerCount >= Mesh.Nvf*20)
    terminate( "urg" );
  
  corner_list[cornerCount].index = tt;
  cornerCount++;
  
  corner_data[tt].active = 1;
  
  struct diff_face_data *fd = &diff_face_data[iface];
  fd->cornerCount++;
}

void compute_least_squares_matrix_at_corners()
{
  int failCount = 0;
  int activeCount = 0;
  
  int totMoves = 0;
  
  int icorner;
  for(icorner = 0; icorner < Mesh.Ndt; icorner++)
    {
      struct corner_data *cd = &corner_data[icorner];

      cd->fail = 0;

      if(!cd->active)
        continue;
      
      activeCount++;

      cd->tetra = icorner;
      
      double weights[NUMDIMS+1];
      int boundary = 0;

      int row, col, k;
      for(row = 0; row < NUMDIMS+1; row++)
        {
          int point = Mesh.DT[cd->tetra].p[row];

          if(point < 0)
            {
              for(k = 0; k < NUMDIMS+1; k++)
                for(col = 0; col < NUMDIMS+1; col++)
                  cd->matrix[k][col] = 0;
              
              boundary = 1;
              break;
            }

          double Center[3];
          if(!TimeBinSynchronized[Mesh.DP[point].timebin])
            {
              Center[0] = Mesh.DP[point].x;
              Center[1] = Mesh.DP[point].y;
              Center[2] = Mesh.DP[point].z;
            }
          else
            {
              int particle = Mesh.DP[point].index;

              if(particle >= NumGas && Mesh.DP[point].task == ThisTask)
                particle -= NumGas;

          
              if(Mesh.DP[point].task == ThisTask)
                {
                  for(k = 0; k < 3; k++)
                    Center[k] = SphP[particle].Center[k];
#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
                  if(pIsReflective[point])
                    for(k = 0; k < 3; k++)
                      Center[k] = refCenter[point*3+k];
#endif
                }
              else
                {
                  for(k = 0; k < 3; k++)
                    Center[k] = PrimExch[particle].Center[k];
                }
            }

          double weight = 1.;
          weights[row] = weight;

          double dist[NUMDIMS];
          double xtmp; dist[0] = NEAREST_X( Center[0] - Mesh.DTC[icorner].cx );
#if NUMDIMS > 1
          double ytmp; dist[1] = NEAREST_Y( Center[1] - Mesh.DTC[icorner].cy );
#endif
#if NUMDIMS > 2
          double ztmp; dist[2] = NEAREST_Z( Center[2] - Mesh.DTC[icorner].cz );
#endif

          cd->matrix[row][0] = weight;
          for(col = 0; col < NUMDIMS; col++)
            cd->matrix[row][col+1] = dist[col] * weight;
        }

      if(boundary)
        continue;

      double matrixT[(NUMDIMS+1)][(NUMDIMS+1)];
      double matrix[(NUMDIMS+1)*(NUMDIMS+1)];
      for(row = 0; row < NUMDIMS+1; row++)
        for(col = 0; col < NUMDIMS+1; col++)
          {
            matrixT[row][col] = cd->matrix[col][row];
            int idx = row*(NUMDIMS+1) + col;
            matrix[idx] = 0;
            for(k = 0; k < NUMDIMS+1; k++)
              matrix[idx] += cd->matrix[k][row] * cd->matrix[k][col];        
          }

      int s;
      gsl_matrix_view m = gsl_matrix_view_array( matrix, NUMDIMS+1, NUMDIMS+1 );
      gsl_permutation *perm = gsl_permutation_alloc( NUMDIMS+1 );
      gsl_linalg_LU_decomp( &m.matrix, perm, &s );

      if( gsl_linalg_LU_det( &m.matrix, s ) != 0. )
        {
          double matrix_inverse[(NUMDIMS+1)*(NUMDIMS+1)];
          gsl_matrix_view minv = gsl_matrix_view_array( matrix_inverse, NUMDIMS+1, NUMDIMS+1 );
          gsl_linalg_LU_invert( &m.matrix, perm, &minv.matrix );
          gsl_permutation_free( perm );

          for(row = 0; row < NUMDIMS+1; row++)
            for(col = 0; col < NUMDIMS+1; col++)
              {
                cd->matrix[col][row] = 0;
                for(k = 0; k < NUMDIMS+1; k++)
                  cd->matrix[col][row] += matrix_inverse[row*(NUMDIMS+1)+k] * matrixT[k][col];

                cd->matrix[col][row] *= weights[row];
              }
        }
      else
        {
          for(row = 0; row < NUMDIMS+1; row++)
            {
              cd->matrix[row][0] = 1./(NUMDIMS+1);
              for(col = 1; col < NUMDIMS+1; col++)
                cd->matrix[row][col] = 0.;
            }
        }

      for(row = 0; row < NUMDIMS+1; row++)
        {
#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
          if(cd->matrix[row][0] < -0.01 || pIsReflective[Mesh.DT[cd->tetra].p[row]])
#else
          if(cd->matrix[row][0] < -0.01)
#endif
            {
              cd->fail = 1;
              break;
            }
        }

      if(cd->fail)
        failCount++;
      
      for(k = 0; k < 3; k++)
        cd->bfld[k] = 0;
      cd->rho = 0;
      cd->chi = 0;
      
      double Bmin[3], Bmax[3];
      
      for(k = 0; k < 3; k++)
        {
          Bmin[k] = +MAX_DOUBLE_NUMBER;
          Bmax[k] = -MAX_DOUBLE_NUMBER;
        }

      int p;
      for(p = 0; p < NUMDIMS+1; p++)
        {
          int point = Mesh.DT[cd->tetra].p[p];
          int particle = Mesh.DP[point].index;
  
          if(particle >= NumGas && Mesh.DP[point].task == ThisTask)
            particle -= NumGas;

          double *B, Rho, Chi;
          if(Mesh.DP[point].task == ThisTask)
            {
              B = SphP[particle].B;
              Rho = SphP[particle].Density;
              Chi = SphP[particle].CR_Chi;
            }
          else
            {
              B = PrimExch[particle].B;
              Rho = PrimExch[particle].Density;
              Chi = PrimExch[particle].CR_Chi;
            }
      
          struct corner_data *cd = &corner_data[icorner];          
          for(k = 0; k < 3; k++)
            {
              cd->bfld[k] += B[k] * cd->matrix[p][0];
              
              if(B[k] < Bmin[k])
                Bmin[k] = B[k];
              if(B[k] > Bmax[k])
                Bmax[k] = B[k];
            }
          cd->rho += Rho * cd->matrix[p][0];
          cd->chi += Chi * cd->matrix[p][0];
        }
      
      for(k = 0; k < 3; k++)
        {
          if(cd->bfld[k] > Bmax[k])
            cd->bfld[k] = Bmax[k];
          if(cd->bfld[k] < Bmin[k])
            cd->bfld[k] = Bmin[k];
        }


    }
  
  int failCountSum, activeCountSum, totMovesSum;
  MPI_Reduce( &failCount, &failCountSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
  MPI_Reduce( &activeCount, &activeCountSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
  MPI_Reduce( &totMoves, &totMovesSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
  
  mpi_printf( "COSMIC_RAYS: STREAMING: %d out of %d active corners were bad, totMoves=%d.\n", failCountSum, activeCountSum, totMovesSum );
}

static double get_vector_distance( double a[3], double b[3] )
{
  double dx = a[0] - b[0];
  double dy = a[1] - b[1];
  double dz = a[2] - b[2];
  
  return sqrt( dx*dx + dy*dy + dz*dz );
}

static double get_area_triangle( double p1[3], double p2[3], double p3[3] )
{
  double a = get_vector_distance( p1, p2 );
  double b = get_vector_distance( p2, p3 );
  double c = get_vector_distance( p3, p1 );
      
  double s = 0.5 * (a + b + c);
  double prod = s * (s-a) * (s-b) * (s-c);
  
  if(prod < 0.)
    return 0.; // can be a tiny bit below zero for degenerate triangles because of roundoff errors
  else
    return sqrt( prod );
}

void compute_geometry_of_interface( int iface )
{
  struct diff_face_data *fd = &diff_face_data[iface];
  
  if(!fd->active)
    return;

  fd->failWeight = 0;

  int p1 = Mesh.VF[iface].p1;
  int p2 = Mesh.VF[iface].p2;
  
  double nx = Mesh.DP[p2].x - Mesh.DP[p1].x;
  double ny = Mesh.DP[p2].y - Mesh.DP[p1].y;
  double nz = Mesh.DP[p2].z - Mesh.DP[p1].z;
  double nn = sqrt( nx*nx + ny*ny + nz*nz );
  
  nx /= nn;
  ny /= nn;
  nz /= nn;
  
  fd->nx = nx;
  fd->ny = ny;
  fd->nz = nz;

  // need an ortonormal basis
  if(fd->nx != 0 || fd->ny != 0)
    {
      fd->mx = -fd->ny;
      fd->my = fd->nx;
      fd->mz = 0;
    }
  else
    {
      fd->mx = 1;
      fd->my = 0;
      fd->mz = 0;
    }

  double mm = sqrt(fd->mx * fd->mx + fd->my * fd->my + fd->mz * fd->mz);
  fd->mx /= mm;
  fd->my /= mm;
  fd->mz /= mm;

  fd->px = fd->ny * fd->mz - fd->nz * fd->my;
  fd->py = fd->nz * fd->mx - fd->nx * fd->mz;
  fd->pz = fd->nx * fd->my - fd->ny * fd->mx;

#ifdef TWODIMS
  int icorner;
  for(icorner = 0; icorner < fd->cornerCount; icorner++)
    {
      struct corner_list *cl = &corner_list[fd->cornerFirst+icorner];
      cl->weight = 0.5;
    }
#else
  // compute weights to get value at center from values at corners
  double cold[3];
  cold[0] = 0.5 * (Mesh.DTC[corner_list[fd->cornerFirst+fd->cornerCount-1].index].cx + Mesh.DTC[corner_list[fd->cornerFirst].index].cx );
  cold[1] = 0.5 * (Mesh.DTC[corner_list[fd->cornerFirst+fd->cornerCount-1].index].cy + Mesh.DTC[corner_list[fd->cornerFirst].index].cy );
  cold[2] = 0.5 * (Mesh.DTC[corner_list[fd->cornerFirst+fd->cornerCount-1].index].cz + Mesh.DTC[corner_list[fd->cornerFirst].index].cz );
  
  double fc[3];
  fc[0] = Mesh.VF[iface].cx;
  fc[1] = Mesh.VF[iface].cy;
  fc[2] = Mesh.VF[iface].cz;

  int icorner;
  for(icorner = 0; icorner < fd->cornerCount; icorner++)
    {
      struct corner_list *cl = &corner_list[fd->cornerFirst+icorner];
      
      double p[3];
      p[0] = Mesh.DTC[cl->index].cx;
      p[1] = Mesh.DTC[cl->index].cy;
      p[2] = Mesh.DTC[cl->index].cz;
      
      double area = get_area_triangle( cold, fc, p );
      
      int icurr = fd->cornerFirst + icorner;
      int inext = fd->cornerFirst + ((icorner+1) % fd->cornerCount);

      cold[0] = 0.5 * (Mesh.DTC[corner_list[icurr].index].cx + Mesh.DTC[corner_list[inext].index].cx );
      cold[1] = 0.5 * (Mesh.DTC[corner_list[icurr].index].cy + Mesh.DTC[corner_list[inext].index].cy );
      cold[2] = 0.5 * (Mesh.DTC[corner_list[icurr].index].cz + Mesh.DTC[corner_list[inext].index].cz );

      area += get_area_triangle( cold, fc, p );
      cl->weight = area / Mesh.VF[iface].area;
    }

  double bx = 0;
  double by = 0;
  double bz = 0;

  double grady = 0.;
  double gradz = 0.;
  for(icorner = 0; icorner < fd->cornerCount; icorner++)
    {
      struct corner_list *cl = &corner_list[fd->cornerFirst+icorner];
      struct corner_data *cd = &corner_data[cl->index];
      
      int k;
      double grad_corner[3];
      for(k = 0; k < NUMDIMS; k++)
        grad_corner[k] = cd->CREnergyDensityGrad[k];
      for(k = NUMDIMS; k < 3; k++)
        grad_corner[k] = 0.;
      
      grady += (grad_corner[0] * fd->mx + grad_corner[1] * fd->my + grad_corner[2] * fd->mz) * cl->weight;
      gradz += (grad_corner[0] * fd->px + grad_corner[1] * fd->py + grad_corner[2] * fd->pz) * cl->weight;
      
      bx += (cd->bfld[0] * fd->nx + cd->bfld[1] * fd->ny + cd->bfld[2] * fd->nz) * cl->weight;
      by += (cd->bfld[0] * fd->mx + cd->bfld[1] * fd->my + cd->bfld[2] * fd->mz) * cl->weight;
      bz += (cd->bfld[0] * fd->px + cd->bfld[1] * fd->py + cd->bfld[2] * fd->pz) * cl->weight;
    }
  
  
  /*
  // orient coordinate system in interface along gradient
  double mx = grady * fd->mx + gradz * fd->px;
  double my = grady * fd->my + gradz * fd->py;
  double mz = grady * fd->mz + gradz * fd->pz;
  */
  
  double mx = by * fd->mx + bz * fd->px;
  double my = by * fd->my + bz * fd->py;
  double mz = bz * fd->mz + bz * fd->pz;
  
  mm = sqrt(mx*mx + my*my + mz*mz);
  if(mm > 0)
    {
      fd->mx = mx / mm;
      fd->my = my / mm;
      fd->mz = mz / mm;

      fd->px = fd->ny * fd->mz - fd->nz * fd->my;
      fd->py = fd->nz * fd->mx - fd->nx * fd->mz;
      fd->pz = fd->nx * fd->my - fd->ny * fd->mx;
    }
#endif

  for(icorner = 0; icorner < fd->cornerCount; icorner++)
    {
      struct corner_list *cl = &corner_list[fd->cornerFirst+icorner];
      struct corner_data *cd = &corner_data[cl->index];

      if(cd->fail)
        fd->failWeight += cl->weight;
    }
}

void compute_quantities_at_interface( int iface )
{
  struct diff_face_data *fd = &diff_face_data[iface];
  
  if(!fd->active)
    return;
  
  int k;
  for(k = 0; k < 3; k++)
    fd->bfld[k] = 0;
  fd->rho = 0;
  
  fd->chi = 0;
  
  int icorner;
  for(icorner = 0; icorner < fd->cornerCount; icorner++)
    {
      struct corner_list *cl = &corner_list[fd->cornerFirst+icorner];
      struct corner_data *cd = &corner_data[cl->index];

      fd->bfld[0] += (cd->bfld[0] * fd->nx + cd->bfld[1] * fd->ny + cd->bfld[2] * fd->nz) * cl->weight;
      fd->bfld[1] += (cd->bfld[0] * fd->mx + cd->bfld[1] * fd->my + cd->bfld[2] * fd->mz) * cl->weight;
      fd->bfld[2] += (cd->bfld[0] * fd->px + cd->bfld[1] * fd->py + cd->bfld[2] * fd->pz) * cl->weight;
      
      fd->rho += cd->rho * cl->weight;
      fd->chi += cd->chi * cl->weight;
    }

  // normalize b-field
  fd->btot = sqrt( fd->bfld[0]*fd->bfld[0] + fd->bfld[1]*fd->bfld[1] + fd->bfld[2]*fd->bfld[2] );
  if(fd->btot > 0.)
    for(k = 0; k < 3; k++)
      fd->bfld[k] /= fd->btot;
}

void compute_cr_energy_and_grad_at_corners()
{
  int icorner, k;
  for(icorner = 0; icorner < Mesh.Ndt; icorner++)
    {
      struct corner_data *cd = &corner_data[icorner];
      
      if(!cd->active)
        continue;

      cd->CREnergyDensity = 0;
      for(k = 0; k < NUMDIMS; k++)
        cd->CREnergyDensityGrad[k] = 0;
      
      tetra *tt = &Mesh.DT[cd->tetra];
      int p;
      for(p = 0; p < NUMDIMS+1; p++)
        {
          int point = tt->p[p];
          
          if(point < 0)
            break;
          
          int particle = Mesh.DP[point].index;
          
          if(particle >= NumGas && Mesh.DP[point].task == ThisTask)
            particle -= NumGas;
          
          double crEnergyDensity;
          if(Mesh.DP[point].task == ThisTask)
            {
              // for reflective boundaries the value is constant across the interface, for outflow boundaries we
              // set it to zero. alternatively one could set it to a constant background value here
#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
              if(pIsReflective[point] == 2)
                crEnergyDensity = 0;
              else
#endif
                crEnergyDensity = SphP[particle].CR_Energy / SphP[particle].Volume;
            }
          else
            crEnergyDensity = CR_EnergyExch[particle] / PrimExch[particle].Volume;
          
          cd->CREnergyDensity += crEnergyDensity * cd->matrix[p][0];

          for(k = 0; k < NUMDIMS; k++)
            cd->CREnergyDensityGrad[k] += crEnergyDensity * cd->matrix[p][k+1];          
        }
    }
}

double compute_chi( int i )
{
#ifdef COSMIC_RAYS_STREAMING_GLOBAL_CHI
  return All.CR_Chi;
#else
  double L = 0.1 * All.BoxSize; // set to r200 of object for cosmological runs
  double chi = SphP[i].CR_Pressure * 2. * get_cell_radius(i) / (L*L);
  
  if(chi < 1e-30)
    chi = 1e-30;
  
  return chi;
#endif
}

double get_timestep_streaming( int i )
{
  double rad = get_cell_radius( i );
  return rad * rad * SphP[i].CR_Chi / SphP[i].CR_Pressure;
}

void streaming_compute_effective_velocities()
{
  exchange_cr_energy();
  
  compute_cr_energy_and_grad_at_corners();
  
  int iface;
  for(iface = 0; iface < Mesh.Nvf; iface++)
    {
      struct diff_face_data *fd = &diff_face_data[iface];
      
      if(!fd->active)
        continue;
      
      double grad[3];
      grad[0] = 0;
      grad[1] = 0;
      grad[2] = 0; // this gradient always points perpendicular to the magnetic field
      
      fd->crEnergyDensity = 0;
      
      /* compute gradient, grad[0] is the normal gradient, grad[1] and grad[2] the gradient in the interface */
      int icorner;
      for(icorner = 0; icorner < fd->cornerCount; icorner++)
        {
          struct corner_list *cl = &corner_list[fd->cornerFirst+icorner];
          struct corner_data *cd = &corner_data[cl->index];
          
          int k;
          double grad_corner[3];
          for(k = 0; k < NUMDIMS; k++)
            grad_corner[k] = cd->CREnergyDensityGrad[k];
          for(k = NUMDIMS; k < 3; k++)
            grad_corner[k] = 0.;

          double gx = (grad_corner[0] * fd->nx + grad_corner[1] * fd->ny + grad_corner[2] * fd->nz);
          double gy = (grad_corner[0] * fd->mx + grad_corner[1] * fd->my + grad_corner[2] * fd->mz);
          double gz = (grad_corner[0] * fd->px + grad_corner[1] * fd->py + grad_corner[2] * fd->pz);
          
          grad[0] += gx * cl->weight;
          grad[1] += gy * cl->weight;
          grad[2] += gz * cl->weight;
          
          fd->crEnergyDensity += cd->CREnergyDensity * cl->weight;
        }
      
      /* compute streaming velocity */
      fd->tanh = tanh( (fd->bfld[0] * grad[0] + fd->bfld[1] * grad[1] + fd->bfld[2] * grad[2]) / fd->chi );
      fd->vel_alfven = fd->btot / sqrt( fd->rho );
      fd->gradb = fd->bfld[0] * grad[0] + fd->bfld[1] * grad[1] + fd->bfld[2] * grad[2];
      
      if(gsl_isnan(fd->tanh * fd->vel_alfven))
        {
          printf( "svel=%g, bfld=%g,%g,%g, grad=%g,%g,%g, btot=%g, rho=%g\n", fd->tanh * fd->vel_alfven, fd->bfld[0], fd->bfld[1], fd->bfld[2], grad[0], grad[1], grad[2], fd->btot, fd->rho );
          terminate( "HAR" );
        }
    }
}

void streaming_explicit()
{
#ifdef COSMIC_RAYS_STREAMING_EXPLICIT
  mpi_printf( "COSMIC_RAYS: STREAMING: Doing explicit streaming.\n" );
#endif
  streaming_compute_effective_velocities();
  
  int iface;
  
  int NCRflux, MaxNCRflux;
  MaxNCRflux = Mesh.Indi.AllocFacNflux;
  NCRflux = 0;
  CRFluxList = mymalloc_movable(&CRFluxList, "FluxList", MaxNCRflux * sizeof(struct crflux_list_data));
  
  for(iface = 0; iface < Mesh.Nvf; iface++)
    {
      struct diff_face_data *fd = &diff_face_data[iface];
      
      if(!fd->active)
        continue;

      double flux = fd->tanh * fd->vel_alfven * fd->bfld[0] * fd->crEnergyDensity;

#ifndef COSMIC_RAYS_STREAMING_EXPLICIT
      flux *= -1;
#endif

      if(flux != 0.)
        {
          int k;
          for(k = 0; k < 2; k++)
            {
              double dir = +1. - k * 2.; /* switches between -1 and +1 */
              int q = (k == 0 ? Mesh.VF[iface].p1 : Mesh.VF[iface].p2);
              int p = Mesh.DP[q].index;

              if(Mesh.DP[q].task == ThisTask)
                {
                  if(Mesh.DP[q].index >= NumGas)     /* this is a local ghost point */
                    {
                      if(Mesh.DP[Mesh.VF[iface].p1].ID == Mesh.DP[Mesh.VF[iface].p2].ID)    /* this may happen for reflective points */
                        continue;
                      p -= NumGas;
                    }

                  SphP[p].CR_Energy += dir * fd->area * fd->dt * flux;
                }
              else
                {
                  if(NCRflux >= MaxNCRflux)
                    {
                      Mesh.Indi.AllocFacNflux *= ALLOC_INCREASE_FACTOR;
                      MaxNCRflux = Mesh.Indi.AllocFacNflux;
                      CRFluxList = myrealloc_movable(CRFluxList, MaxNCRflux * sizeof(struct crflux_list_data));

                      if(NCRflux >= MaxNCRflux)
                        terminate("NCRflux >= MaxNCRflux");
                    }

                  CRFluxList[NCRflux].task = Mesh.DP[q].task;
                  CRFluxList[NCRflux].index = Mesh.DP[q].originalindex;
                  CRFluxList[NCRflux].dCREnergy = dir * fd->area * fd->dt * flux;
                  NCRflux++;
                }
            }
        }
    }
  
  /* exchange and apply fluxes */
  int i, j, nimport, ngrp, recvTask;
  mysort(CRFluxList, NCRflux, sizeof(struct crflux_list_data), crflux_list_data_compare);

  for(j = 0; j < NTask; j++)
    Send_count[j] = 0;

  for(i = 0; i < NCRflux; i++)
    Send_count[CRFluxList[i].task]++;

  if(Send_count[ThisTask] > 0)
    terminate("Send_count[ThisTask]");

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  struct crflux_list_data *CRFluxListGet = (struct crflux_list_data *) mymalloc("FluxListGet", nimport * sizeof(struct crflux_list_data));
  
  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              MPI_Sendrecv(&CRFluxList[Send_offset[recvTask]],
                           Send_count[recvTask] * sizeof(struct crflux_list_data), MPI_BYTE,
                           recvTask, TAG_DENS_A,
                           &CRFluxListGet[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct crflux_list_data), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }
  
  for(i = 0; i < nimport; i++)
    {
      int p = CRFluxListGet[i].index;
      
      SphP[p].CR_Energy += CRFluxListGet[i].dCREnergy;
    }
  
  myfree(CRFluxListGet);
  myfree(CRFluxList);
}

int crflux_list_data_compare(const void *a, const void *b)
{
  if(((struct crflux_list_data *) a)->task < (((struct crflux_list_data *) b)->task))
    return -1;

  if(((struct crflux_list_data *) a)->task > (((struct crflux_list_data *) b)->task))
    return +1;

  return 0;
}

void free_stuff()
{
  /* free temporary arrays */
  myfree(corner_data);
  myfree(corner_list);
  myfree(diff_face_data);
  myfree(CR_EnergyExch);
#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
  myfree(pIsReflective);
  myfree(refCenter);
#endif
}

void streaming_implicit()
{
  mpi_printf( "COSMIC_RAYS: STREAMING: Doing implicit streaming.\n" );
  
  struct matrix_data md;
  struct hypre_data hd;
  struct cr_energy_log logB, logE;
  struct solver_data sd;
  
  hd.initialized = 0;
  
  allocate_rows_offsets_and_indizes( &md );
  compute_rows_offsets_and_indizes( &md );
  
  streaming_implicit_init_solver( &md, &hd, &sd );
  
  streaming_implicit_compute_energy_log( &logB, &md );
  
  int needSplit = 0;
  MPI_Comm diffComm;
  double lin_res_norm = 0;
  
  int iter;
  for(iter = 0; iter <= STREAMING_MAX_ITER_NONLINEAR; iter++)
    {
      if(streaming_implicit_update_solution( &md, &hd, &sd, iter, lin_res_norm ))
        break;
      
      if(iter == STREAMING_MAX_ITER_NONLINEAR)
        terminate( "COSMIC_RAYS: STREAMING: Nonlinear iteration did not converge after %d iterations, stopping.\n", iter );
      
      set_matrix_coefficients( &md, &hd, &sd );
      
      if(iter == 0)
        {
          int task;
          for(task = 0; task < NTask; task++)
            if(md.offsets[task] > md.offsets[task+1]-1)
              needSplit++;

          if(needSplit == NTask)
            {
              free_hypre_data( &hd );
              streaming_implicit_free_solver( &sd );
              free_rows_offsets_and_indizes( &md );
              return;
            }

          mpi_printf( "COSMIC_RAYS: STREAMING: Using %d of %d cores for implicit streaming.\n", NTask-needSplit, NTask );
          
          if(needSplit)
            MPI_Comm_split( MPI_COMM_WORLD, md.local_row_count > 0 ? 1 : 2, ThisTask, &diffComm );
          else
            diffComm = MPI_COMM_WORLD;
        }
      
      streaming_implicit_solve_matrix( &md, &hd, 1, diffComm, &lin_res_norm, &sd );
  }
  
  streaming_implicit_compute_energy_log( &logE, &md );
  
  if(iter > 0)
    free_hypre_data( &hd );
  streaming_implicit_free_solver( &sd );
  free_rows_offsets_and_indizes( &md );
  
  if(needSplit)
    {
      MPI_Barrier( diffComm );
      MPI_Comm_free( &diffComm );
    }
  
  mpi_printf( "COSMIC_RAYS: STREAMING: Total CR Energy (min/max density) of all participating cells changed from %g (%g,%g) to %g (%g,%g)\n",
    logB.sum, logB.min, logB.max, logE.sum, logE.min, logE.max );
  
  MPI_Barrier( MPI_COMM_WORLD );
}

void streaming_implicit_init_solver( struct matrix_data *md, struct hypre_data *hd, struct solver_data *sd )
{
  sd->residuals = (double*)mymalloc( "residuals", md->local_row_count*sizeof(double) );
  sd->initial_cr_energy_density = (double*)mymalloc( "initial_cr_energy_density", md->local_row_count*sizeof(double) );
  sd->last_cr_energy_density = (double*)mymalloc( "last_cr_energy_density", md->local_row_count*sizeof(double) );
  sd->current_cr_energy_density = (double*)mymalloc( "current_cr_energy_density", md->local_row_count*sizeof(double) );
  sd->delta_cr_energy_density = (double*)mymalloc( "delta_cr_energy_density", md->local_row_count*sizeof(double) );
  sd->g = (double*)mymalloc( "g", md->local_row_count*sizeof(double) );
  
  double norm_initial_cr_energy_density = 0;
  for(int row = 0; row < md->local_row_count; row++)
    {
      int cell = md->local_part_index[row];
      sd->initial_cr_energy_density[row] = SphP[cell].CR_Energy / SphP[cell].Volume;
      norm_initial_cr_energy_density += sd->initial_cr_energy_density[row] * sd->initial_cr_energy_density[row];
      sd->last_cr_energy_density[row] = 0;
      sd->delta_cr_energy_density[row] = 0;
    }
  
  MPI_Allreduce( &norm_initial_cr_energy_density, &sd->norm_initial_cr_energy_density, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  sd->norm_initial_cr_energy_density = sqrt( sd->norm_initial_cr_energy_density );
}

void streaming_implicit_free_solver( struct solver_data *sd )
{
  myfree( sd->g );
  myfree( sd->delta_cr_energy_density );
  myfree( sd->current_cr_energy_density );
  myfree( sd->last_cr_energy_density );
  myfree( sd->initial_cr_energy_density );
  myfree( sd->residuals );
}

void streaming_implicit_compute_residuals( struct matrix_data *md, struct hypre_data *hd, struct solver_data *sd )
{
  /* we assume that CR_Energy is set to the value we want to calculate the residual for */
  streaming_explicit(); /* also updates streaming velocities */
  
  for(int row = 0; row < md->local_row_count; row++)
    {
      int cell = md->local_part_index[row];
      sd->residuals[row] = SphP[cell].CR_Energy / SphP[cell].Volume - sd->initial_cr_energy_density[row];
    }
}

int streaming_implicit_update_solution( struct matrix_data *md, struct hypre_data *hd, struct solver_data *sd, int iter, double lin_res_norm )
{
  /* store the current cosmic ray energy somewhere else, as we will use CR_Energy for other purposes */
  for(int row = 0; row < md->local_row_count; row++)
    {
      int cell = md->local_part_index[row];
      sd->current_cr_energy_density[row] = SphP[cell].CR_Energy / SphP[cell].Volume;
      
      /* update chi */
      SphP[cell].CR_Chi = compute_chi( cell );
    }

  if(iter == 0)
    {
      /* only do initializations */
      streaming_implicit_compute_residuals( md, hd, sd );
      
      double norm_residual = 0;
      for(int row = 0; row < md->local_row_count; row++)
        {
          int cell = md->local_part_index[row];
          sd->last_cr_energy_density[row] = sd->current_cr_energy_density[row];
          
          /* restore old value */
          SphP[cell].CR_Energy = sd->current_cr_energy_density[row] * SphP[cell].Volume;
          
          norm_residual += sd->residuals[row] * sd->residuals[row];
        }
      
      double norm_residual_tot;
      MPI_Allreduce( &norm_residual, &norm_residual_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
      sd->norm_initial_residual = sqrt( norm_residual_tot );
      sd->norm_last_residual = 0.5 * norm_residual_tot;
      
      if(sd->norm_initial_residual == 0.)
        return 1;
      else
        return 0;
    }

  double local_slope = 0;
  for(int row = 0; row < md->local_row_count; row++)
    {
      local_slope += sd->g[row] * sd->delta_cr_energy_density[row];
    }
  
  double slope;
  MPI_Allreduce( &local_slope, &slope, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  
  double stepsize = 1.0;
  double stepsize_last = 0.;
  double norm = 0;
  double norm_last = 0;
  while(1)
    {
      if(stepsize < 1e-8)
        break;
      
      for(int row = 0; row < md->local_row_count; row++)
        {
          int cell = md->local_part_index[row];
          SphP[cell].CR_Energy = (sd->current_cr_energy_density[row] + stepsize * sd->delta_cr_energy_density[row]) * SphP[cell].Volume;
        }
      
      streaming_implicit_compute_residuals( md, hd, sd );
      
      double local_norm = 0;
      for(int row = 0; row < md->local_row_count; row++)
        {
          local_norm += sd->residuals[row] * sd->residuals[row];
        }
      
      local_norm *= 0.5;
      MPI_Allreduce( &local_norm, &norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
      
      if(norm < sd->norm_last_residual + 1e-4 * stepsize * slope)
        break;
      
      if(stepsize == 1.0)
        {
          stepsize_last = stepsize;
          stepsize = dmax( 0.1 * stepsize, -slope/(2.*(norm-sd->norm_last_residual-slope)) );
          norm_last = norm;
        }
      else
        {
          double rhs1 = norm - sd->norm_last_residual - stepsize * slope;
          double rhs2 = norm_last - sd->norm_last_residual - stepsize_last * slope;
          
          double a = (rhs1/(stepsize*stepsize) - rhs2/(stepsize_last*stepsize_last)) / (stepsize-stepsize_last);
          double b = (-stepsize_last*rhs1/(stepsize*stepsize) + stepsize*rhs2/(stepsize_last*stepsize_last)) / (stepsize-stepsize_last);
          
          stepsize_last = stepsize;
          norm_last = norm;
          
          if(a == 0.)
            {
              stepsize = dmax( 0.1 * stepsize, dmin( 0.5 * stepsize, -slope/(2.*b) ) );
            }
          else
            {
              double stepsize_new;
              double disc = b*b - 3.*a*slope;
              
              if(disc < 0.)
                stepsize_new = 0.5 * stepsize;
              else if (b <= 0.)
                stepsize_new = (sqrt(disc) - b) / (3.*a);
              else
                stepsize_new = -slope/(b+sqrt(disc));
              
              stepsize = dmax( 0.1 * stepsize, dmin( 0.5 * stepsize, stepsize_new ) );
            }
        }
    }
  
  double norm_difference = 0;
  for(int row = 0; row < md->local_row_count; row++)
    {
      int cell = md->local_part_index[row];
      
      double current_cr_energy_density = sd->current_cr_energy_density[row] + stepsize * sd->delta_cr_energy_density[row];
      SphP[cell].CR_Energy = current_cr_energy_density * SphP[cell].Volume;
      
      double difference = current_cr_energy_density - sd->last_cr_energy_density[row];
      norm_difference += difference * difference;
      
      sd->last_cr_energy_density[row] = current_cr_energy_density;
    }
  
  double norm_difference_tot;
  MPI_Allreduce( &norm_difference, &norm_difference_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  norm_difference_tot = sqrt( norm_difference_tot );
  
  sd->norm_last_residual = norm;
  double norm_residual_tot = sqrt( 2. * norm );
  
  /* check size of the residual and progress as stopping condition */
  if(iter == 0)
    {
      sd->norm_initial_residual = norm_residual_tot;
      return 0;
    }
  else
    {
      struct cr_energy_log logC;
      streaming_implicit_compute_energy_log( &logC, md );
      
      mpi_printf( "COSMIC_RAYS: STREAMING iter=%02d, res norm=%8.5e, rel. res norm=%8.5e, rel change=%8.5e, stepsize=%8.5e, lin. res norm=%8.5e, energies = %8.5e, %8.5e, %8.5e\n", iter, norm_residual_tot,
        norm_residual_tot / sd->norm_initial_residual, norm_difference_tot / sd->norm_initial_cr_energy_density, stepsize, lin_res_norm,
        logC.sum, logC.min, logC.max );
      
      double rel_norm = norm_residual_tot / sd->norm_initial_residual;
      double rel_change = norm_difference_tot / sd->norm_initial_cr_energy_density;
      if(rel_norm < 1e-3 || rel_change < 1e-6)
        {
          return 1;
        }
      else
        return 0;
    }
}

int streaming_implicit_solve_matrix( struct matrix_data *md, struct hypre_data *hd, int useMultigridPreconditionier, MPI_Comm diffComm, double *res_norm, struct solver_data *sd )
{
  *res_norm = -1;
  
  if(hd->rows > 0)
    {
      /* initialize, run, and free the HYPRE matrix solver */
      HYPRE_IJMatrix A;
      HYPRE_ParCSRMatrix parcsr_A;
      HYPRE_IJVector b;
      HYPRE_ParVector par_b;
      HYPRE_IJVector x;
      HYPRE_ParVector par_x;

      HYPRE_Solver solver, precond;
      
      int ilower = md->offsets[ThisTask];
      int iupper = md->offsets[ThisTask+1]-1;
      
      HYPRE_IJMatrixCreate(diffComm, ilower, iupper, ilower, iupper, &A);
      HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
      HYPRE_IJMatrixInitialize(A);

      HYPRE_IJVectorCreate(diffComm, ilower, iupper, &b);
      HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
      HYPRE_IJVectorInitialize(b);
      
      HYPRE_IJVectorCreate(diffComm, ilower, iupper, &x);
      HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
      HYPRE_IJVectorInitialize(x);
      
      HYPRE_IJVectorSetValues(b, hd->nRows, hd->rows, hd->bval);
      HYPRE_IJVectorSetValues(x, hd->nRows, hd->rows, hd->xval);
      HYPRE_IJMatrixSetValues(A, hd->nRows, hd->columnCounts, hd->rows, hd->cols, hd->vals);
      
      HYPRE_IJMatrixAssemble(A);
      HYPRE_IJMatrixGetObject(A, (void**) &parcsr_A);

      HYPRE_IJVectorAssemble(b);
      HYPRE_IJVectorGetObject(b, (void **) &par_b);

      HYPRE_IJVectorAssemble(x);
      HYPRE_IJVectorGetObject(x, (void **) &par_x);
      
      if(useMultigridPreconditionier == 2)
        {
          /* we failed, lets write out the matrix and stop */
          mkdir( "matrix", 02755 );

          HYPRE_IJMatrixPrint(A, "matrix/IJ.out.A");
          HYPRE_IJVectorPrint(b, "matrix/IJ.out.b");
          HYPRE_IJVectorPrint(x, "matrix/IJ.out.x");
          
          MPI_Barrier( diffComm );
          terminate( "HYPRE did not converge." );
        }
      
      int nIter;
      double final_res_norm;

      HYPRE_ParCSRFlexGMRESCreate(diffComm, &solver);
      HYPRE_ParCSRFlexGMRESSetTol(solver, STREAMING_ACCURACY);
      HYPRE_ParCSRFlexGMRESSetPrintLevel(solver, 0); // 3
      HYPRE_ParCSRFlexGMRESSetLogging(solver, 1);
      
      if(useMultigridPreconditionier)
        {
          HYPRE_ParCSRFlexGMRESSetKDim(solver, 500);
          HYPRE_ParCSRFlexGMRESSetMaxIter(solver, 500);
        }
      else
        {
          HYPRE_ParCSRFlexGMRESSetKDim(solver, 200);
          HYPRE_ParCSRFlexGMRESSetMaxIter(solver, 200);
        }

      if(useMultigridPreconditionier)
        {
          HYPRE_BoomerAMGCreate(&precond);
#ifdef TWODIMS
          HYPRE_BoomerAMGSetCoarsenType(precond, 6);
          HYPRE_BoomerAMGSetStrongThreshold(precond, 0.25);
#else
          HYPRE_BoomerAMGSetCoarsenType(precond, 6);
          HYPRE_BoomerAMGSetInterpType(precond, 6); // test 6/7 here??
          HYPRE_BoomerAMGSetStrongThreshold(precond, 0.25); // increase to 0.9?
          HYPRE_BoomerAMGSetPMaxElmts(precond, 4);
#endif
          HYPRE_BoomerAMGSetRelaxType(precond, 18);
          HYPRE_BoomerAMGSetNumSweeps(precond, 1);
          HYPRE_BoomerAMGSetTol(precond, 0.0);
          HYPRE_BoomerAMGSetPrintLevel(precond, 0); // 1
          HYPRE_BoomerAMGSetMaxIter(precond, 2);

          HYPRE_ParCSRFlexGMRESSetPrecond(solver,
                                     HYPRE_BoomerAMGSolve,
                                     HYPRE_BoomerAMGSetup,
                                     precond);
        }
        
      HYPRE_ParCSRFlexGMRESSetup(solver, parcsr_A, par_b, par_x);
      
      HYPRE_ParCSRFlexGMRESSolve(solver, parcsr_A, par_b, par_x);
      
      HYPRE_ParCSRFlexGMRESGetNumIterations(solver, &nIter);
      HYPRE_ParCSRFlexGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
      
      *res_norm = final_res_norm;
  
      if(useMultigridPreconditionier)
        HYPRE_BoomerAMGDestroy(precond);
        
      HYPRE_ParCSRFlexGMRESDestroy(solver);

      streaming_implicit_extract_solution(&x, md, sd);
      /* writes the result back to CR_Energy */
      
      HYPRE_IJMatrixDestroy(A);
      HYPRE_IJVectorDestroy(b);
      HYPRE_IJVectorDestroy(x);
    }

  double error = compute_error( md, hd, sd );
  double b = compute_abs( md, hd, hd->bval );

  double res = error / b;

  int success = res < STREAMING_ACCURACY;  
  if((res < 1e-6) && !success)
    {
      mpi_printf( "COSMIC_RAYS: STREAMING: Residual is larger than %g, but still smaller than %g, so we accept it and continue.\n", STREAMING_ACCURACY, 1e-6 );
      success = 1;
    }
  
  return success;
}

double compute_abs( struct matrix_data *md, struct hypre_data *hd, double * val )
{
  double sum = 0;
  
  int row;
  for(row = 0; row < hd->nRows; row++)
    {
      sum += val[row] * val[row];
    }
  
  double totsum;
  MPI_Allreduce( &sum, &totsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  
  return sqrt(totsum);
}

double compute_error( struct matrix_data *md, struct hypre_data *hd, struct solver_data *sd )
{
  double *residuals = (double*)mymalloc( "residuals", hd->nRows * sizeof(double) );
  
  struct residual_exchange_request {
    int task;
    int index;
    int row;
    double val;
  } *res_exch_req_send_unordered, *res_exch_req_send, *res_exch_req_get;
  
  struct residual_exchange_result {
    int row;
    double val;
    double x;
    double residual;
  } *res_exch_res_send, *res_exch_res_get;
  
  res_exch_req_send_unordered = (struct residual_exchange_request *)mymalloc( "res_exch_req_send_unordered", md->cd_count*sizeof(struct residual_exchange_request) );
  int Nsend = 0;
  
  int row;
  int colCount = 0;
  for(row = 0; row < hd->nRows; row++)
    {
      residuals[row] = -hd->bval[row];
      
      int col;
      for(col = 0; col < hd->columnCounts[row]; col++)
        {
          int task = hd->tasks[colCount];
          
          if(task == ThisTask)
            {
              int particle = md->local_part_index[hd->cols[colCount] - md->offsets[task]];
              residuals[row] += hd->vals[colCount] * SphP[particle].CR_Energy / SphP[particle].Volume;
              sd->g[row] += hd->vals[colCount] * sd->residuals[row];
            }
          else
            {
              res_exch_req_send_unordered[Nsend].task = task;
              res_exch_req_send_unordered[Nsend].index = hd->cols[colCount] - md->offsets[task];
              res_exch_req_send_unordered[Nsend].row = row;
              res_exch_req_send_unordered[Nsend].val = hd->vals[colCount];
              Nsend++;
            }
          colCount++;
        }
    }
  
  int sendCount[NTask];
  int s, t;
  for(t = 0; t < NTask; t++)
    sendCount[t] = 0;
  for(s = 0; s < Nsend; s++)
    sendCount[res_exch_req_send_unordered[s].task]++;
  
  int sendOffset[NTask];
  sendOffset[0] = 0;
  for(t = 1; t < NTask; t++)
    sendOffset[t] = sendOffset[t-1] + sendCount[t-1];
  
  for(t = 0; t < NTask; t++)
    sendCount[t] = 0;
  
  res_exch_req_send = (struct residual_exchange_request *)mymalloc( "res_exch_req_send", Nsend*sizeof(struct residual_exchange_request) );
  
  for(s = 0; s < Nsend; s++)
    {
      int task = res_exch_req_send_unordered[s].task;
      int index = sendOffset[task] + sendCount[task];
      sendCount[task]++;
      res_exch_req_send[index] = res_exch_req_send_unordered[s];
    }
  
  int recvCount[NTask];
  MPI_Alltoall(sendCount, 1, MPI_INT, recvCount, 1, MPI_INT, MPI_COMM_WORLD);
  
  int recvOffset[NTask];
  recvOffset[0] = 0;
  for(t = 1; t < NTask; t++)
    recvOffset[t] = recvOffset[t-1] + recvCount[t-1];
  
  int nimport = recvOffset[NTask-1] + recvCount[NTask-1];
  
  res_exch_req_get = (struct residual_exchange_request *)mymalloc( "res_exch_req_get", nimport*sizeof(struct residual_exchange_request) );
  
  int ngrp;
  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(sendCount[recvTask] > 0 || recvCount[recvTask] > 0)
            {
              MPI_Sendrecv(&res_exch_req_send[sendOffset[recvTask]],
                           sendCount[recvTask] * sizeof(struct residual_exchange_request), MPI_BYTE,
                           recvTask, TAG_DENS_A,
                           &res_exch_req_get[recvOffset[recvTask]], recvCount[recvTask] * sizeof(struct residual_exchange_request), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }
  
  res_exch_res_send = (struct residual_exchange_result *)mymalloc( "res_exch_res_send", nimport*sizeof(struct residual_exchange_result) );
  
  int i;
  for(i = 0; i < nimport; i++)
    {
      res_exch_res_send[i].row = res_exch_req_get[i].row;
      res_exch_res_send[i].val = res_exch_req_get[i].val;
      
      int particle = md->local_part_index[res_exch_req_get[i].index];
      res_exch_res_send[i].x = SphP[particle].CR_Energy / SphP[particle].Volume;
      res_exch_res_send[i].residual = sd->residuals[res_exch_req_get[i].index];
    }
  
  res_exch_res_get = (struct residual_exchange_result *)mymalloc( "res_exch_res_get", nimport*sizeof(struct residual_exchange_result) );
  
  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(sendCount[recvTask] > 0 || recvCount[recvTask] > 0)
            {
              MPI_Sendrecv(&res_exch_res_send[recvOffset[recvTask]],
                           recvCount[recvTask] * sizeof(struct residual_exchange_result), MPI_BYTE,
                           recvTask, TAG_DENS_A,
                           &res_exch_res_get[sendOffset[recvTask]], sendCount[recvTask] * sizeof(struct residual_exchange_result), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }
  
  for(i = 0; i < Nsend; i++)
    {
      residuals[ res_exch_res_get[i].row ] += res_exch_res_get[i].val * res_exch_res_get[i].x;
      sd->g[ res_exch_res_get[i].row ] += res_exch_res_get[i].val * res_exch_res_get[i].residual;
    }
  
  myfree(res_exch_res_get);
  myfree(res_exch_res_send);
  myfree(res_exch_req_get);
  myfree(res_exch_req_send);
  myfree(res_exch_req_send_unordered);
  
  double error = 0;
  for(row = 0; row < hd->nRows; row++)
    {
      error += residuals[row] * residuals[row];
    }
  
  myfree( residuals );
  
  double toterror;
  MPI_Allreduce( &error, &toterror, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  
  return sqrt(toterror);
}

void allocate_rows_offsets_and_indizes( struct matrix_data *md )
{
  md->local_row_index = (int*) mymalloc( "local_row_index", sizeof(int) * NumGas );
  int icell;
  for(icell = 0; icell < NumGas; icell++)
    md->local_row_index[icell] = -1;
  
  md->local_part_index = (int*) mymalloc( "local_part_index", sizeof(int) * NumGas );
  md->first_column = (int*) mymalloc( "first_column", sizeof(int) * NumGas );
  md->last_column = (int*) mymalloc( "last_column", sizeof(int) * NumGas );
  
  md->column_data = (struct column_data*)mymalloc( "column_data", sizeof(struct column_data) * NumGas * 20 * NUMDIMS );
  md->cd_count = 0;
  
  md->imported_particle_indizes = (int*) mymalloc( "imported_particle_indizes", sizeof(int) * Mesh_nimport );
  for(icell = 0; icell < Mesh_nimport; icell++)
    md->imported_particle_indizes[icell] = -1;
  
  md->offsets = (int*) mymalloc( "offsets", sizeof(int) * (NTask+1) );
  
  md->local_row_count = 0;
}

struct external_cells {
  int task;
  int index;
};

static int external_cells_compare(const void *a, const void *b)
{
  if(((struct external_cells *) a)->task < (((struct external_cells *) b)->task))
    return -1;

  if(((struct external_cells *) a)->task > (((struct external_cells *) b)->task))
    return +1;

  return 0;
}

void compute_rows_offsets_and_indizes( struct matrix_data *md )
{
  /* look through all interfaces, mark cells that are needed */
  int N_external_cells = 0;
  int Nmax_external_cells = Mesh.Indi.AllocFacNflux;
  struct external_cells *external_cells;
  external_cells = (struct external_cells*)mymalloc_movable( &external_cells, "external_cells", sizeof(struct external_cells) * Nmax_external_cells );
  
  int iface;
  for(iface = 0; iface < Mesh.Nvf; iface++)
    {
      struct diff_face_data *fd = &diff_face_data[iface];
      
      if(!fd->active)
        continue;
      
      int k;
      for(k = 0; k < 2; k++)
        {
          int p = (k == 0 ? Mesh.VF[iface].p1 : Mesh.VF[iface].p2);
          
          if(Mesh.DP[p].task == ThisTask)
            {
              int particle = Mesh.DP[p].index;
              if(particle >= NumGas) particle -= NumGas;
          
              if(md->local_row_index[particle] == -1)
                {
                  /* we still have to add this */
                  md->local_row_index[particle] = md->local_row_count;
                  md->local_part_index[md->local_row_count] = particle;
                  md->first_column[md->local_row_count] = -1;
                  md->last_column[md->local_row_count] = -1;
                  md->local_row_count++;
                }
            }
          else
            {
              if(N_external_cells >= Nmax_external_cells)
                {
                  Nmax_external_cells *= ALLOC_INCREASE_FACTOR;
                  external_cells = myrealloc_movable(external_cells, Nmax_external_cells * sizeof(struct external_cells));
                }
          
              external_cells[N_external_cells].task = Mesh.DP[p].task;
              external_cells[N_external_cells].index = Mesh.DP[p].originalindex;
              N_external_cells++;
            }
        }

      int icorner;
      for(icorner = 0; icorner < fd->cornerCount; icorner++)
        {
          struct corner_list *cl = &corner_list[fd->cornerFirst+icorner];
          struct corner_data *cd = &corner_data[cl->index];
          tetra *tt = &Mesh.DT[cd->tetra];
  
          int icol;
          for(icol = 0; icol < NUMDIMS+1; icol++)
            {
              int p = tt->p[icol];

              if(Mesh.DP[p].task == ThisTask)
                {
                  int particle = Mesh.DP[p].index;
                  if(particle >= NumGas) particle -= NumGas;
      
                  if(md->local_row_index[particle] == -1)
                    {
                      /* we still have to add this */
                      md->local_row_index[particle] = md->local_row_count;
                      md->local_part_index[md->local_row_count] = particle;
                      md->first_column[md->local_row_count] = -1;
                      md->last_column[md->local_row_count] = -1;
                      md->local_row_count++;
                    }
                }
              else
                {
                  if(N_external_cells >= Nmax_external_cells)
                    {
                      Nmax_external_cells *= ALLOC_INCREASE_FACTOR;
                      external_cells = myrealloc_movable(external_cells, Nmax_external_cells * sizeof(struct external_cells));
                    }
      
                  external_cells[N_external_cells].task = Mesh.DP[p].task;
                  external_cells[N_external_cells].index = Mesh.DP[p].originalindex;
                  N_external_cells++;
                }
            }
        }
    }
  
  /* exchange list of external cells that need to be included */
  int i, j, ngrp, recvTask, nimport;
  mysort(external_cells, N_external_cells, sizeof(struct external_cells), external_cells_compare);

  for(j = 0; j < NTask; j++)
    Send_count[j] = 0;

  for(i = 0; i < N_external_cells; i++)
    Send_count[external_cells[i].task]++;

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  struct external_cells *external_cells_get = (struct external_cells *) mymalloc("external_cells_get", nimport * sizeof(struct external_cells));

  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              MPI_Sendrecv(&external_cells[Send_offset[recvTask]],
                           Send_count[recvTask] * sizeof(struct external_cells), MPI_BYTE,
                           recvTask, TAG_DENS_A,
                           &external_cells_get[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct external_cells), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  /* add cells that are needed from extern but are not already in our local list */
  for(i = 0; i < nimport; i++)
    {
      int particle = external_cells_get[i].index;
      
      if(md->local_row_index[particle] == -1)
        {
          /* we still have to add this */
          md->local_row_index[particle] = md->local_row_count;
          md->local_part_index[md->local_row_count] = particle;
          md->first_column[md->local_row_count] = -1;
          md->last_column[md->local_row_count] = -1;
          md->local_row_count++;
        }
    }

  myfree(external_cells_get);
  myfree(external_cells);
  
  /* compute global offsets */
  md->offsets[NTask] = 0;

  MPI_Allgather( &md->local_row_count, 1, MPI_INT, md->offsets, 1, MPI_INT, MPI_COMM_WORLD );

  int sum = 0;
  for(i = 0; i <= NTask; i++)
    {
      int tmp = sum;
      sum += md->offsets[i];
      md->offsets[i] = tmp;
    }
  
  exchange_row_indizes( md );
}

void exchange_row_indizes( struct matrix_data *md )
{
  int listp;
  int i, j, p, task, off;
  int ngrp, recvTask, place;
  int *tmpIndizesExch, *tmpIndizesRecv;

  tmpIndizesExch = (int *) mymalloc("tmpIndizesExch", Mesh_nexport * sizeof(int));

  /* prepare data for export */
  for(j = 0; j < NTask; j++)
    Mesh_Send_count[j] = 0;

  for(i = 0; i < NumGasInMesh; i++)
    {
      p = List_InMesh[i];

      listp = List_P[p].firstexport;
      while(listp >= 0)
        {
          if((task = ListExports[listp].origin) != ThisTask)
            {
              place = ListExports[listp].index;
              off = Mesh_Send_offset[task] + Mesh_Send_count[task]++;
              
              if(md->local_row_index[place] == -1)
                tmpIndizesExch[off] = -1;
              else
                tmpIndizesExch[off] = md->local_row_index[place] + md->offsets[ThisTask];
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
              tmpIndizesRecv = (int *) mymalloc("tmpIndizesRecv", Mesh_Recv_count[recvTask] * sizeof(int));

              /* get the values */
              MPI_Sendrecv(&tmpIndizesExch[Mesh_Send_offset[recvTask]],
                           Mesh_Send_count[recvTask] * sizeof(int), MPI_BYTE,
                           recvTask, TAG_DENS_A, tmpIndizesRecv, Mesh_Recv_count[recvTask] * sizeof(int), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

              for(i = 0; i < Mesh_Recv_count[recvTask]; i++)
                {
                  if(Mesh_Recv_offset[recvTask] + i >= Mesh_nimport)
                    terminate( "blurb" );
                  
                  md->imported_particle_indizes[Mesh_Recv_offset[recvTask] + i] = tmpIndizesRecv[i];
                }

              myfree(tmpIndizesRecv);
            }
        }
    }

  myfree(tmpIndizesExch);
}

void free_hypre_data( struct hypre_data *hd )
{
  if(hd->nRows > 0)
    {
      myfree(hd->tasks);
      myfree(hd->xval);
      myfree(hd->bval);
      myfree(hd->vals);
      myfree(hd->cols);
      myfree(hd->rows);
      myfree(hd->columnCounts);      
    }
}

void free_rows_offsets_and_indizes( struct matrix_data *md )
{
  /* free temporary arrays */
  myfree( md->offsets );
  myfree( md->imported_particle_indizes );
  myfree( md->column_data );
  myfree( md->last_column );
  myfree( md->first_column );
  myfree( md->local_part_index );
  myfree( md->local_row_index );
}

static int external_elements_compare(const void *a, const void *b)
{
  if(((struct external_element *) a)->task < (((struct external_element *) b)->task))
    return -1;

  if(((struct external_element *) a)->task > (((struct external_element *) b)->task))
    return +1;

  return 0;
}

void set_matrix_coefficients( struct matrix_data *md, struct hypre_data *hd, struct solver_data *sd )
{
  struct external_elements ee;
  
  ee.N_external_elements = 0;
  ee.Nmax_external_elements = Mesh.Indi.AllocFacNflux;
  ee.elements = (struct external_element*)mymalloc_movable( &ee.elements, "external_elements", sizeof(struct external_element) * ee.Nmax_external_elements );
  
  md->cd_count = 0;
  
  int row;
  for(row = 0; row < md->local_row_count; row++)
    {
      int particle = md->local_part_index[row];
      
      struct column_data *cd = &md->column_data[md->cd_count];
      md->first_column[particle] = md->cd_count;
      md->last_column[particle] = md->cd_count;
      md->cd_count++;
      
      if(md->cd_count >= NumGas * 20 * NUMDIMS)
        terminate( "baaad" );
      
      cd->next_column = -1;
      cd->value = 1.0;
      cd->index = row + md->offsets[ThisTask];
      cd->task = ThisTask;
    }
  
  /* go through all interfaces and compute matrix entries */
  int iface;
  for(iface = 0; iface < Mesh.Nvf; iface++)
    {
      struct diff_face_data *fd = &diff_face_data[iface];
      
      if(!fd->active)
        continue;
      
#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
      if((pIsReflective[Mesh.VF[iface].p1] == 1) || (pIsReflective[Mesh.VF[iface].p2] == 1))
        continue;
#endif
      
      // we use full normal gradients
      int icorner;
      for(icorner = 0; icorner < fd->cornerCount; icorner++)
        {
          struct corner_list *cl = &corner_list[fd->cornerFirst+icorner];
          struct corner_data *cd = &corner_data[cl->index];
          tetra *tt = &Mesh.DT[cd->tetra];

          if(cd->fail || cl->weight == 0.)
            continue;

          int icol;
          for(icol = 0; icol < NUMDIMS+1; icol++)
            {
              int col = tt->p[icol];
              
              double val = -fd->area * fd->dt * fd->vel_alfven * fd->bfld[0] * fd->tanh * cd->matrix[icol][0] * cl->weight;
              
              double fac = -fd->area * fd->dt * fd->vel_alfven * fd->bfld[0] * fd->crEnergyDensity * (1. - fd->tanh * fd->tanh) * cl->weight / fd->chi;
#ifdef TWODIMS
              val += fac * (cd->matrix[icol][1] * fd->nx + cd->matrix[icol][2] * fd->ny) * fd->bfld[0];
              val += fac * (cd->matrix[icol][1] * fd->mx + cd->matrix[icol][2] * fd->my) * fd->bfld[1];
#else
              val += fac * (cd->matrix[icol][1] * fd->nx + cd->matrix[icol][2] * fd->ny + cd->matrix[icol][3] * fd->nz) * fd->bfld[0];
              val += fac * (cd->matrix[icol][1] * fd->mx + cd->matrix[icol][2] * fd->my + cd->matrix[icol][3] * fd->mz) * fd->bfld[1];
              val += fac * (cd->matrix[icol][1] * fd->px + cd->matrix[icol][2] * fd->py + cd->matrix[icol][3] * fd->pz) * fd->bfld[2];
#endif

#ifndef COSMIC_RAYS_STREAMING_GLOBAL_CHI
              val += fd->area * fd->dt * fd->vel_alfven * fd->bfld[0] * (1. - fd->tanh * fd->tanh) * fd->gradb / fd->chi * cd->matrix[icol][0] * cl->weight;
#endif
              add_matrix_element( Mesh.VF[iface].p1, col, +val, md, &ee );
              add_matrix_element( Mesh.VF[iface].p2, col, -val, md, &ee );
            }
        }
    }
  
  /* exchange matrix entries for external rows */
  int i, j, nimport, ngrp, recvTask;
  mysort(ee.elements, ee.N_external_elements, sizeof(struct external_element), external_elements_compare);

  for(j = 0; j < NTask; j++)
    Send_count[j] = 0;

  for(i = 0; i < ee.N_external_elements; i++)
    Send_count[ee.elements[i].task]++;

  if(Send_count[ThisTask] > 0)
    terminate("Send_count[ThisTask]");

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  struct external_element *external_elements_get = (struct external_element *) mymalloc("external_elements_get", nimport * sizeof(struct external_element));
  
  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              MPI_Sendrecv(&ee.elements[Send_offset[recvTask]],
                           Send_count[recvTask] * sizeof(struct external_element), MPI_BYTE,
                           recvTask, TAG_DENS_A,
                           &external_elements_get[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct external_element), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }
 
  for(i = 0; i < nimport; i++)
    {
      int row = external_elements_get[i].index;
      int column = external_elements_get[i].column_index;
      int particle = md->local_part_index[row - md->offsets[ThisTask]];
      int originaltask = external_elements_get[i].originaltask;
      
      struct column_data *cd = NULL;
      
      int entry = md->first_column[particle];
      while(entry >= 0)
        {
          cd = &md->column_data[entry];
          
          if(cd->index == column)
            break;
          
          entry = cd->next_column;
        }
      
      if(!cd || cd->index != column)
        {
          if(md->cd_count >= NumGas * 20 * NUMDIMS)
            terminate( "list too long" );
      
          cd = &md->column_data[md->cd_count];
          if(md->first_column[particle] == -1)
            md->first_column[particle] = md->cd_count;
          else
            md->column_data[md->last_column[particle]].next_column = md->cd_count;

          md->last_column[particle] = md->cd_count;
          md->cd_count++;
      
          cd->index = column;
          cd->value = external_elements_get[i].value;
          cd->next_column = -1;
          cd->task = originaltask;
        }
      else
        {
          cd->value += external_elements_get[i].value;
          
          if(cd->task != originaltask)
            {
              printf( "ThisTask=%d, cd->task=%d, originaltask=%d, entry=%d\n", ThisTask, cd->task, originaltask, entry );
              terminate( "failed" );
            }
        }
    }
  
  myfree( external_elements_get );
  myfree( ee.elements );
  
  /* check if we contribute */
  if(md->local_row_count > 0)
    {
      /* build matrix and send data to HYPRE */
      if(!hd->initialized)
        {
          hd->nRows = md->local_row_count;
      
          hd->columnCounts = (int*)mymalloc( "columnCounts", md->local_row_count * sizeof(int) );
          hd->rows = (int*)mymalloc( "rows", md->local_row_count * sizeof(int) );
          hd->cols = (int*)mymalloc( "cols", md->cd_count * sizeof(int) );
  
          hd->vals = (double*)mymalloc( "vals", md->cd_count * sizeof(double) );
          hd->bval = (double*)mymalloc( "bval", md->local_row_count * sizeof(double) );
          hd->xval = (double*)mymalloc( "xval", md->local_row_count * sizeof(double) );
          hd->tasks = (int*)mymalloc( "point", md->cd_count * sizeof(int) );
          
          hd->initialized = 1;
        }
  
      int colCount = 0;
      for(row = 0; row < md->local_row_count; row++)
        {
          int particle = md->local_part_index[row];

          hd->bval[row] = -sd->residuals[row];
          //hd->xval[row] = 0;
          hd->xval[row] = sd->delta_cr_energy_density[row];
          hd->rows[row] = row + md->offsets[ThisTask];
  
          hd->columnCounts[row] = 0;
          
          double fac = 1.; /* changing this changes the weighting for the global residual!!! */
          
          int column_index = md->first_column[particle];
          while( column_index >= 0 )
            {
              struct column_data *cd = &md->column_data[column_index];
              
              hd->columnCounts[row]++;
              hd->cols[colCount] = cd->index;
              hd->vals[colCount] = cd->value * fac;
              hd->tasks[colCount] = cd->task;
              colCount++;
              
              column_index = cd->next_column;
            }

          hd->bval[row] *= fac;
        }
    }
  else
    {
      hd->nRows = 0;
      hd->columnCounts = 0;
      hd->rows = 0;
      hd->cols = 0;
      hd->vals = 0;
      hd->bval = 0;
      hd->xval = 0;
    }
}

void streaming_implicit_compute_energy_log( struct cr_energy_log *lg, struct matrix_data *md )
{
  double CR_Energy, CR_EnergyAll;
  double CR_EDensMin, CR_EDensMax;
  double CR_EDensMinAll, CR_EDensMaxAll;
  
  CR_Energy = 0;
  CR_EDensMin = +MAX_DOUBLE_NUMBER;
  CR_EDensMax = -MAX_DOUBLE_NUMBER;

  int row;
  for(row = 0; row < md->local_row_count; row++)
    {
      int cell = md->local_part_index[row];
      
      CR_Energy += SphP[cell].CR_Energy;
      double CR_EnergyDensity = SphP[cell].CR_Energy / SphP[cell].Volume;
      if(CR_EnergyDensity < CR_EDensMin)
        CR_EDensMin = CR_EnergyDensity;
      if(CR_EnergyDensity > CR_EDensMax)
        CR_EDensMax = CR_EnergyDensity;
    }
  
  MPI_Reduce( &CR_Energy, &CR_EnergyAll, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
  MPI_Reduce( &CR_EDensMin, &CR_EDensMinAll, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
  MPI_Reduce( &CR_EDensMax, &CR_EDensMaxAll, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
  
  if(ThisTask == 0)
    {
      lg->sum = CR_EnergyAll;
      lg->min = CR_EDensMinAll;
      lg->max = CR_EDensMaxAll;
    }
}

void point_get_center( int p, double *Center )
{
  if(!TimeBinSynchronized[Mesh.DP[p].timebin])
    {
      Center[0] = Mesh.DP[p].x;
      Center[1] = Mesh.DP[p].y;
      Center[2] = Mesh.DP[p].z;
      return;
    }
  
  int particle = Mesh.DP[p].index;
  int j;
  if(Mesh.DP[p].task == ThisTask)
    {
      if(particle >= NumGas)
        particle -= NumGas;

      for(j = 0; j < 3; j++) 
        Center[j] = SphP[particle].Center[j];

#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
      if(pIsReflective[p])
        for(j = 0; j < 3; j++)
          Center[j] = refCenter[p*3+j];
#endif
    }
  else
    {
      for(j = 0; j < 3; j++)
        Center[j] = PrimExch[particle].Center[j];
    }
}

void add_matrix_element( int point_row, int point_column, double val, struct matrix_data *md, struct external_elements *ee )
{
  int column;
  
  if(Mesh.DP[point_column].task == ThisTask)
    {
      int particle = Mesh.DP[point_column].index;
      if(particle >= NumGas)
        particle -= NumGas;
        
      column = md->local_row_index[particle] + md->offsets[ThisTask];

      if(column < 0 || column > md->offsets[NTask]-1)
        {
          printf( "column fail: %d, last offset: %d, particle: %d, local_row_index: %d, offset: %d\n", column, md->offsets[NTask], particle, md->local_row_index[particle], md->offsets[ThisTask] );
          terminate( "Z" );
        }
    }
  else
    {
      int particle = Mesh.DP[point_column].index;
      
      if(particle < 0)
        terminate( "A" );
      if(particle >= Mesh_nimport)
        terminate( "B" );
      
      column = md->imported_particle_indizes[particle];
      
      if(column < 0 || column > md->offsets[NTask]-1)
        {
          printf( "column fail: %d, last offset: %d\n", column, md->offsets[NTask] );
          terminate( "C" );
        }
    }

  if(Mesh.DP[point_row].task == ThisTask)
    {
      /* add to local list */
      int particle = Mesh.DP[point_row].index;
      if(particle >= NumGas)
        particle -= NumGas;
      
      int entry = md->first_column[particle];
      while(entry >= 0)
        {
          struct column_data *cd = &md->column_data[entry];
          
          if(cd->index == column)
            break;
          
          entry = cd->next_column;
        }
      
      if(entry == -1)
        {
          /* need new entry */
          if(md->cd_count >= NumGas * 20 * NUMDIMS)
            terminate( "list too long" );

          struct column_data *cd = &md->column_data[md->cd_count];
          if(md->first_column[particle] == -1)
            md->first_column[particle] = md->cd_count;
          else
            md->column_data[md->last_column[particle]].next_column = md->cd_count;

          md->last_column[particle] = md->cd_count;
          md->cd_count++;
      
          cd->index = column;
          cd->value = val / SphP[particle].Volume;
          cd->next_column = -1;
          cd->task = Mesh.DP[point_column].task;
        }
      else
        {
          struct column_data *cd = &md->column_data[entry];
          cd->value += val / SphP[particle].Volume;
          
          if(cd->task != Mesh.DP[point_column].task)
            terminate( "shit" );
        }
    }
  else
    {
      /* create entry to send to other task */
      if(ee->N_external_elements >= ee->Nmax_external_elements)
        {
          ee->Nmax_external_elements *= ALLOC_INCREASE_FACTOR;
          ee->elements = myrealloc_movable( ee->elements, sizeof(struct external_element) * ee->Nmax_external_elements );
        }
      
      struct external_element *el = &ee->elements[ee->N_external_elements];
      ee->N_external_elements++;
      
      el->task = Mesh.DP[point_row].task;
      el->index = md->imported_particle_indizes[Mesh.DP[point_row].index];
      el->column_index = column;
      el->value = val / PrimExch[Mesh.DP[point_row].index].Volume;
      el->originaltask = Mesh.DP[point_column].task;
    }
}

void streaming_implicit_extract_solution( HYPRE_IJVector* x, struct matrix_data *md, struct solver_data *sd )
{
  /* extract solution from HYPRE */
  int rows[md->local_row_count];
  double xval[md->local_row_count];

  for(int row = 0; row <md->local_row_count; row++)
    rows[row] = row + md->offsets[ThisTask];

  HYPRE_IJVectorGetValues(*x,md->local_row_count, rows, xval);
  
  for(int row = 0; row < md->local_row_count; row++)
    {
      sd->delta_cr_energy_density[row] = xval[row];
      
      //int particle = md->local_part_index[row];
      //SphP[particle].CR_Energy += xval[row] * SphP[particle].Volume;
    }
}

#endif
