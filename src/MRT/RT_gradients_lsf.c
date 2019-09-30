/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/voronoi_gradients.c
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

#include "../allvars.h"
#include "../proto.h"
#include "RT.h"
#include "RT_proto.h"



#ifdef MRT
#if defined(VORONOI_DYNAMIC_UPDATE) || defined(AMR_CONNECTIONS)


static double *minvalues, *maxvalues;

static void limit_gradients_RT();

static double boundaryX(double dx);
static double boundaryY(double dy);
static double boundaryZ(double dz);
static void inline add_row(double X[NUMDIMS][NUMDIMS], double y[NUMDIMS], int source_row, double fac, int target_row)
{
  y[target_row] += fac * y[source_row];

  for(int i = 0; i < NUMDIMS; i++)
    {
      X[target_row][i] += fac * X[source_row][i];
    }
}


/* note, we now here that X is symmetric, and that we can pivot on the diagonal elements 
 * NOTE: the matrix X and the vector y are modified by this routine.
 */
static void solve_matrix_problem(double X[NUMDIMS][NUMDIMS], double y[NUMDIMS], double grad[NUMDIMS])
 {
#if NUMDIMS==2
  int perm[NUMDIMS];

  if(fabs(X[0][0]) > fabs(X[1][1]))
    {
      perm[0] = 0;
      perm[1] = 1;
    }
  else
    {
      perm[0] = 1;
      perm[1] = 0;
    }

  add_row(X, y, perm[0], -X[perm[1]][perm[0]] / X[perm[0]][perm[0]], perm[1]);

  grad[perm[1]] = y[perm[1]] / X[perm[1]][perm[1]];
  grad[perm[0]] = (y[perm[0]] - X[perm[0]][perm[1]] * grad[perm[1]]) / X[perm[0]][perm[0]];

#else

  int perm[NUMDIMS];

  if(fabs(X[2][2]) > fabs(X[1][1]) && fabs(X[2][2]) > fabs(X[0][0]))
    {
      perm[0] = 2;
      perm[1] = 0;
      perm[2] = 1;
    }
  else if(fabs(X[1][1]) > fabs(X[0][0]))
    {
      perm[0] = 1;
      perm[1] = 0;
      perm[2] = 2;
    }
  else
    {
      perm[0] = 0;
      perm[1] = 1;
      perm[2] = 2;
    }

  add_row(X, y, perm[0], -X[perm[1]][perm[0]] / X[perm[0]][perm[0]], perm[1]);
  add_row(X, y, perm[0], -X[perm[2]][perm[0]] / X[perm[0]][perm[0]], perm[2]);

  if(fabs(X[perm[1]][perm[1]]) < fabs(X[perm[2]][perm[2]]))
    {
      int p = perm[1];
      perm[1] = perm[2];
      perm[2] = p;
    }

  add_row(X, y, perm[1], -X[perm[2]][perm[1]] / X[perm[1]][perm[1]], perm[2]);

  grad[perm[2]] = y[perm[2]] / X[perm[2]][perm[2]];
  grad[perm[1]] = (y[perm[1]] - X[perm[1]][perm[2]] * grad[perm[2]]) / X[perm[1]][perm[1]];
  grad[perm[0]] = (y[perm[0]] - X[perm[0]][perm[1]] * grad[perm[1]] - X[perm[0]][perm[2]] * grad[perm[2]]) / X[perm[0]][perm[0]];

#endif
}



void calculate_gradients_RT(void)
{
  //  TIMER_START(CPU_GRADIENTS);

  mpi_printf("RT: Calculating Gradients...\n");

  minvalues = mymalloc("gradmin", NumGas * N_Grad_RT * sizeof(double));
  maxvalues = mymalloc("gradmax", NumGas * N_Grad_RT * sizeof(double));

  struct matrix_vec_data
  {
    double X[NUMDIMS][NUMDIMS];    /* input matrix */
    double y[NUMDIMS];             /* input vector */
    double grad[NUMDIMS];          /* output */
  }
  *mdata;

  mdata = mymalloc("mdata", N_Grad_RT * sizeof(struct matrix_vec_data));

  double *Value =  mymalloc("Value", N_Grad_RT * sizeof(double));


  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      for(int k = 0; k < N_Grad_RT; k++)
        {
          minvalues[i * N_Grad_RT + k] = +MAX_REAL_NUMBER;
          maxvalues[i * N_Grad_RT + k] = -MAX_REAL_NUMBER;

	  Value[k] = *(MyFloat *) (((char *) (&SphP[i])) + grad_elements_RT[k].offset);
        }


      MyDouble *Center = SphP[i].Center;

      /* reset matrix and vector to 0 */
      memset(mdata, 0, N_Grad_RT * sizeof(struct matrix_vec_data));

#ifdef REFLECTIVE_X
      int OutFlowX = 0;
#endif
#ifdef REFLECTIVE_Y
      int OutFlowY = 0;
#endif
#ifdef REFLECTIVE_Z
      int OutFlowZ = 0;
#endif

      int q = SphP[i].first_connection;

      while(q >= 0)
        {
          int dp = DC[q].dp_index;
          int vf = DC[q].vf_index;
          int particle = Mesh.DP[dp].index;

          if(particle < 0)
            {
              /* cell has been removed */
              q = DC[q].next;
              continue;
            }

          if(Mesh.VF[vf].area > 1e-10 * SphP[i].SurfaceArea)
            {
              MyDouble *CenterOther, Mirror[3];

              if(particle >= NumGas && Mesh.DP[dp].task == ThisTask)
                particle -= NumGas;

#ifdef REFLECTIVE_X
              if((Mesh.DP[dp].image_flags & REFL_X_FLAGS) && (Mesh.DP[dp].image_flags & OUTFLOW_X))
              OutFlowX = 1;
#endif
#ifdef REFLECTIVE_Y
              if((Mesh.DP[dp].image_flags & REFL_Y_FLAGS) && (Mesh.DP[dp].image_flags & OUTFLOW_Y))
              OutFlowY = 1;
#endif
#ifdef REFLECTIVE_Z
              if((Mesh.DP[dp].image_flags & REFL_Z_FLAGS) && (Mesh.DP[dp].image_flags & OUTFLOW_Z))
              OutFlowZ = 1;
#endif

              if(Mesh.DP[dp].task == ThisTask)
                {
#ifndef VORONOI_STATIC_MESH
                  if(P[particle].Ti_Current != All.Ti_Current)
                    terminate("surprise! we don't expect this here anymore");
#endif

                  if(P[particle].ID == P[i].ID)
                    {
                      /* mirrored cell, we have to mirror the Center */

                      /* calculate normal vector of the interface */
                      double nx = Mesh.DP[dp].x - P[i].Pos[0];
                      double ny = Mesh.DP[dp].y - P[i].Pos[1];
                      double nz = Mesh.DP[dp].z - P[i].Pos[2];

                      /* perpendicular on the surface */
                      double nn = sqrt(nx * nx + ny * ny + nz * nz);
                      nx /= nn;
                      ny /= nn;
                      nz /= nn;
                      double fx = (Center[0] - Mesh.VF[vf].cx);
                      double fy = (Center[1] - Mesh.VF[vf].cy);
                      double fz = (Center[2] - Mesh.VF[vf].cz);
                      double ff = (fx * nx + fy * ny + fz * nz);

                      double px = Center[0] - ff * nx;
                      double py = Center[1] - ff * ny;
                      double pz = Center[2] - ff * nz;

                      Mirror[0] = 2. * px - Center[0];
                      Mirror[1] = 2. * py - Center[1];
                      Mirror[2] = 2. * pz - Center[2];
                      CenterOther = Mirror;
                    }
                  else
                    CenterOther = SphP[particle].Center;
                }
              else
                CenterOther = PrimExch[particle].Center;

              double norm[3];
              norm[0] = boundaryX(CenterOther[0] - Center[0]);
              norm[1] = boundaryY(CenterOther[1] - Center[1]);
              norm[2] = boundaryZ(CenterOther[2] - Center[2]);

              double dist = sqrt(norm[0] * norm[0] + norm[1] * norm[1] + norm[2] * norm[2]);
              double distinv = 1.0 / dist;
              norm[0] *= distinv;
              norm[1] *= distinv;
              norm[2] *= distinv;

              double weight = Mesh.VF[vf].area;


              for(int k = 0; k < N_Grad_RT; k++)
                {
                  double ValueOther;

                  if(Mesh.DP[dp].task == ThisTask)
		      ValueOther = *(MyFloat *) (((char *) (&SphP[particle])) + grad_elements_RT[k].offset);
		  else
                      ValueOther = *(MyFloat *) (((char *) (&RTPrimExch[particle])) + grad_elements_RT[k].offset_exch);
                    

#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
                      correct_for_reflective_boundaries(&ValueOther, Value[k], grad_elements_RT[k].type, &Mesh.DP[dp].image_flags);
#endif

                  double fac = weight * (ValueOther - Value[k]) / dist;

                  for(int ia = 0; ia < NUMDIMS; ia++)
                    {
                      mdata[k].y[ia] += fac * norm[ia];

                      for(int ib = 0; ib < NUMDIMS; ib++)
                        mdata[k].X[ia][ib] += weight * norm[ia] * norm[ib];
                    }

                  if(ValueOther < minvalues[i * N_Grad_RT + k])
                    minvalues[i * N_Grad_RT + k] = ValueOther;

                  if(ValueOther > maxvalues[i * N_Grad_RT + k])
                    maxvalues[i * N_Grad_RT + k] = ValueOther;
                }
            }

          if(q == SphP[i].last_connection)
            break;

          q = DC[q].next;
        }

      for(int k = 0; k < N_Grad_RT; k++)
        {
          solve_matrix_problem(mdata[k].X, mdata[k].y, mdata[k].grad);

          MySingle *data = (MySingle *) (((char *) (&(SphP[i].RTGrad))) + grad_elements_RT[k].offset_grad);
          for(int j = 0; j < NUMDIMS; j++)
            data[j] = mdata[k].grad[j];
          for(int j = NUMDIMS; j < 3; j++)
            data[j] = 0.;

#ifdef REFLECTIVE_X
          if(OutFlowX)
            data[0] = 0;
#endif
#ifdef REFLECTIVE_Y
          if(OutFlowY)
            data[1] = 0;
#endif
#ifdef REFLECTIVE_Z
          if(OutFlowZ)
            data[2] = 0;
#endif
        }
    }

  myfree(Value);
  myfree(mdata);


#ifndef UNLIMITED_GRADIENTS
  limit_gradients_RT();
#endif

  myfree(maxvalues);
  myfree(minvalues);


//  TIMER_STOP(CPU_GRADIENTS);
}



void limit_gradients_RT(void)
{

  mpi_printf("RT: Limiting gradients...\n");

  point *DP = Mesh.DP;
  face *VF = Mesh.VF;

  for(int i = 0; i < Mesh.Nvf; i++)
    {
      if(DP[VF[i].p1].index < 0 || DP[VF[i].p2].index < 0)
        continue;
      for(int j = 0; j < 2; j++)
        {
          point *p;
#ifdef GRADIENT_LIMITER_DUFFELL
          point *pother;
#endif
          if(j == 0)
            {
              p = &DP[VF[i].p1];
#ifdef GRADIENT_LIMITER_DUFFELL
              pother = &DP[VF[i].p2];
#endif
            }
          else
            {
              p = &DP[VF[i].p2];
#ifdef GRADIENT_LIMITER_DUFFELL
              pother = &DP[VF[i].p1];
#endif
            }
          
          if(p->task == ThisTask && p->index >= 0 && p->index < NumGas)
            {
              int q = p->index;
              if(TimeBinSynchronized[P[q].TimeBinHydro])
                {
                  double d[3];
                  d[0] = VF[i].cx - SphP[q].Center[0];
                  d[1] = VF[i].cy - SphP[q].Center[1];
                  d[2] = VF[i].cz - SphP[q].Center[2];
#if !defined(REFLECTIVE_X)
                  double xtmp;
                  d[0] = NEAREST_X(d[0]);
#endif
#if !defined(REFLECTIVE_Y)
                  double ytmp;
                  d[1] = NEAREST_Y(d[1]);
#endif
#if !defined(REFLECTIVE_Z)
                  double ztmp;
                  d[2] = NEAREST_Z(d[2]);
#endif
                  double value;
                  MySingle *data;
                  if(VF[i].area > 1.0e-10 * SphP[q].SurfaceArea)
                    {
                      for(int k = 0; k < N_Grad_RT; k++)
                        {
			  value = *(MyFloat *) (((char *) (&SphP[q])) + grad_elements_RT[k].offset);
                          
                          data = (MySingle *) (((char *) (&(SphP[q].RTGrad))) + grad_elements_RT[k].offset_grad);
                          
#ifndef GRADIENT_LIMITER_DUFFELL
			  limit_gradient(d, value, minvalues[q * N_Grad_RT + k], maxvalues[q * N_Grad_RT + k], data);
#else 
                          double ValueOther;
                          
                          int particle = pother->index;
                          if(particle >= NumGas && pother->task == ThisTask)
                            particle -= NumGas;

                          if(pother->task == ThisTask)
			    ValueOther = *(MyFloat *) (((char *) (&SphP[particle])) + grad_elements_RT[k].offset);
			  else
			    ValueOther = *(MyFloat *) (((char *) (&RTPrimExch[particle])) + grad_elements_RT[k].offset_exch);

                          limit_gradient(d, value, ValueOther, ValueOther, data);
#endif
                        }
                    }
                }
            }
        }
    }
}

double boundaryX(double dx)
{
#if !defined(REFLECTIVE_X)
  if(dx < -boxHalf_X)
    dx += boxSize_X;
  if(dx > boxHalf_X)
    dx -= boxSize_X;
#endif
  return dx;
}

double boundaryY(double dy)
{
#if !defined(REFLECTIVE_Y)
  if(dy < -boxHalf_Y)
    dy += boxSize_Y;
  if(dy > boxHalf_Y)
    dy -= boxSize_Y;
#endif
  return dy;
}

double boundaryZ(double dz)
{
#if !defined(REFLECTIVE_Z)
  if(dz < -boxHalf_Z)
    dz += boxSize_Z;
  if(dz > boxHalf_Z)
    dz -= boxSize_Z;
#endif
  return dz;
}

#endif
#endif
