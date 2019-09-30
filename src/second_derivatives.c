/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/second_derivatives.c
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

#include "allvars.h"
#include "proto.h"
#include "voronoi.h"

#ifdef SECOND_DERIVATIVES

static int N_Hess = 0;

static struct hessian_elements
{
  size_t offset_grad;           /* offset of the gradient of the quantity in the SphP struct */
  size_t offset_grad_exch;      /* offset of the gradient in the PrimExch struct */
  size_t offset_hessian;        /* offset in the hessian_data struct */
#ifdef SLOPE_LIMIT_HESSIANS
  MyFloat *min_value, *max_value;
#endif
  MyFloat gradient0[3], gradient1[3];
}
hessian_elements[MAXGRADIENTS];
#endif

void init_hessians()
{
#ifdef SECOND_DERIVATIVES
  hessian_init(SphP[0].Grad.drho, GradExch[0].drho, &SphP[0].Hessian.ddrho[0][0]);

  hessian_init(SphP[0].Grad.dvel[0], GradExch[0].dvel[0], &SphP[0].Hessian.ddvelx[0][0]);
  hessian_init(SphP[0].Grad.dvel[1], GradExch[0].dvel[1], &SphP[0].Hessian.ddvely[0][0]);
  hessian_init(SphP[0].Grad.dvel[2], GradExch[0].dvel[2], &SphP[0].Hessian.ddvelz[0][0]);

  hessian_init(SphP[0].Grad.dpress, GradExch[0].dpress, &SphP[0].Hessian.ddpress[0][0]);
/*
#ifdef RT_ADVECT
  hessian_init(SphP[0].Grad.ddensphot, GradExch[0].ddensphot, SphP[0].Hessian.ddensphot[0][0]);
#endif

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
  hessian_init(SphP[0].Grad.dA, GradExch[0].dA, &SphP[0].Hessian.ddA[0][0]);
#endif

#ifdef VARIABLE_GAMMA
  hessian_init(SphP[0].Grad.dgammaE, GradExch[0].dgammaE, &SphP[0].Hessian.ddgammaE[0][0]);
  hessian_init(SphP[0].Grad.dammaC, GradExch[0].dgammaC, &SphP[0].Hessian.ddgammaC[0][0]);
#endif

#ifdef MAXSCALARS
  int k;
  MyFloat *addr;

  for(k = 0; k < N_Scalar; k++)
    {
      addr = (MyFloat *) (((char *) (&SphP[0])) + scalar_elements[k].offset);
      hessian_init(addr, GradExch[0].dscalars[k], &SphP[0].Hessian.dscalars[k][0][0]);
    }
#endif

#ifdef TRACER_FIELD
  hessian_init(SphP[0].Grad.dtracer, GradExch[0].dtracer, &SphP[0].Hessian.ddtracer[0][0]);
#endif
*/

  mpi_printf("%d/%d Hessians used.\n", N_Hess, MAXGRADIENTS);

#endif
}


void hessian_init(MyFloat * addr_grad, MyFloat * addr_grad_exch, double *addr_hessian)
{
#ifdef SECOND_DERIVATIVES
  if(N_Hess == MAXGRADIENTS)
    {
      mpi_printf("Failed to register hessian, maximum of %d already reached\n", MAXGRADIENTS);
      terminate("MAXGRADIENTS reached");
    }


  hessian_elements[N_Hess].offset_grad = ((char *) addr_grad) - ((char *) &SphP[0].Grad);
  hessian_elements[N_Hess].offset_grad_exch = ((char *) addr_grad_exch) - ((char *) &GradExch[0]);
  hessian_elements[N_Hess].offset_hessian = ((char *) addr_hessian) - ((char *) &(SphP[0].Hessian));

  N_Hess++;
#endif
}


void calculate_green_gauss_hessian(tessellation * T)
/*We compute the Hessian matrix at the center of a cell
using the previously computed cell-centered gradients
in the neighboring cells*/
{
#ifdef SECOND_DERIVATIVES
  CPU_Step[CPU_MISC] += measure_time();

  face *VF = T->VF;
  point *DP = T->DP;

  int idx, i, j, k, l, t0, t1, q0, q1, s0, s1;
  MySingle *gradient0, *gradient1;

  MyFloat value;
  double fac0, fac1, n[3], nn, d[3], *data;
  double direct = 0, transpose = 0;


  double facA, facB, c[3];

#ifdef SLOPE_LIMIT_HESSIANS
  for(k = 0; k < N_Hess; k++)
    {
      hessian_elements[k].min_value = mymalloc("gradmin", 3 * NumGas * sizeof(double));
      hessian_elements[k].max_value = mymalloc("gradmax", 3 * NumGas * sizeof(double));
    }
#endif

#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
  unsigned int *image_flags0, *image_flags1;
#endif

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      for(k = 0; k < N_Hess; k++)
        {
          data = (double *) (((char *) (&(SphP[i].Hessian))) + hessian_elements[k].offset_hessian);
          for(l = 0; l < 3; l++)
            {
              for(j = 0; j < 3; j++)
                data[l * 3 + j] = 0;

#ifdef SLOPE_LIMIT_HESSIANS
              hessian_elements[k].min_value[3 * i + l] = +MAX_REAL_NUMBER;
              hessian_elements[k].max_value[3 * i + l] = -MAX_REAL_NUMBER;
#endif
            }
        }
    }


  for(i = 0; i < T->Nvf; i++)
    {
      point *p1 = &DP[VF[i].p1];
      point *p0 = &DP[VF[i].p2];

      if(p0->index < 0 || p1->index < 0)
        continue;

      n[0] = p1->x - p0->x;
      n[1] = p1->y - p0->y;
      n[2] = p1->z - p0->z;

      nn = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);

      for(j = 0; j < 3; j++)
        n[j] /= nn;

      c[0] = VF[i].cx - 0.5 * (p1->x + p0->x);
      c[1] = VF[i].cy - 0.5 * (p1->y + p0->y);
      c[2] = VF[i].cz - 0.5 * (p1->z + p0->z);

      /* one of the sides */
      q0 = p0->index;
      s0 = p1->index;
      t0 = p1->task;
#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
      image_flags0 = &(p1->image_flags);
#endif

      /* the other side */
      q1 = p1->index;
      s1 = p0->index;
      t1 = p0->task;
#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
      image_flags1 = &(p0->image_flags);
#endif

      for(k = 0; k < N_Hess; k++)
        {
          for(j = 0; j < 3; j++)
            {
              hessian_elements[k].gradient0[j] = 0;
              hessian_elements[k].gradient1[j] = 0;
            }
        }



      /* let's get the physical quantities of interest for the particle one one of the sides */
      if(t0 == ThisTask)
        {
          if(s0 >= NumGas)
            s0 -= NumGas;

          if(s0 >= 0)
            {
              for(k = 0; k < N_Hess; k++)
                {
                  gradient0 = (MySingle *) (((char *) (&(SphP[s0].Grad))) + hessian_elements[k].offset_grad);
                  hessian_elements[k].gradient0[0] = gradient0[0];
                  hessian_elements[k].gradient0[1] = gradient0[1];
                  hessian_elements[k].gradient0[2] = gradient0[2];
                }
            }
        }
      else
        {
          for(k = 0; k < N_Hess; k++)
            {
              gradient0 = (MySingle *) (((char *) (&GradExch[s0])) + hessian_elements[k].offset_grad_exch);
              hessian_elements[k].gradient0[0] = gradient0[0];
              hessian_elements[k].gradient0[1] = gradient0[1];
              hessian_elements[k].gradient0[2] = gradient0[2];
            }
        }


      /* let's get the physical quantities of interest for the particle on the other side */
      if(t1 == ThisTask)
        {
          if(s1 >= NumGas)
            s1 -= NumGas;

          if(s1 >= 0)
            {
              for(k = 0; k < N_Hess; k++)
                {
                  gradient1 = (MySingle *) (((char *) (&SphP[s1].Grad)) + hessian_elements[k].offset_grad);
                  hessian_elements[k].gradient1[0] = gradient1[0];
                  hessian_elements[k].gradient1[1] = gradient1[1];
                  hessian_elements[k].gradient1[2] = gradient1[2];
                }

            }
        }
      else
        {
          for(k = 0; k < N_Hess; k++)
            {
              gradient1 = (MySingle *) (((char *) (&GradExch[s1])) + hessian_elements[k].offset_grad_exch);
              hessian_elements[k].gradient1[0] = gradient1[0];
              hessian_elements[k].gradient1[1] = gradient1[1];
              hessian_elements[k].gradient1[2] = gradient1[2];
            }
        }

      /*We EXTRAPOLATE the velocity gradients in case of reflective or non-slip boundaries.
         Although gradients are IDENTICAL on both sides of the surface, we do this because
         gradients do not change continuously accross boundaries like primitive variables do. */
#if defined(REFLECTIVE_X)
      if((*image_flags0 & REFL_X_FLAGS) && !(*image_flags0 & OUTFLOW_X))
        for(j = 0; j < 3; j++)
          {
            hessian_elements[1].gradient0[j] *= 2;
            if((*image_flags0 & REFL_X_FLAGS) && (*image_flags0) && (REFLECTIVE_X == 3))
              {
                hessian_elements[2].gradient0[j] *= 2;
                hessian_elements[3].gradient0[j] *= 2;
              }
          }
#endif
#if defined(REFLECTIVE_Y)
      if((*image_flags0 & REFL_Y_FLAGS) && !(*image_flags0 & OUTFLOW_Y))
        for(j = 0; j < 3; j++)
          {
            hessian_elements[2].gradient0[j] *= 2;
            if((*image_flags0 & REFL_Y_FLAGS) && (*image_flags0) && (REFLECTIVE_Y == 3))
              {
                hessian_elements[1].gradient0[j] *= 2;
                hessian_elements[3].gradient0[j] *= 2;
              }
          }
#endif
#if defined(REFLECTIVE_Z)
      if((*image_flags0 & REFL_Z_FLAGS) && !(*image_flags0 & OUTFLOW_Z))
        for(j = 0; j < 3; j++)
          {
            hessian_elements[3].gradient0[j] *= 2;
            if((*image_flags0 & REFL_Z_FLAGS) && (*image_flags0) && (REFLECTIVE_Z == 3))
              {
                hessian_elements[1].gradient0[j] *= 2;
                hessian_elements[2].gradient0[j] *= 2;
              }
          }
#endif

#if defined(REFLECTIVE_X)
      if((*image_flags1 & REFL_X_FLAGS) && !(*image_flags1 & OUTFLOW_X))
        for(j = 0; j < 3; j++)
          {
            hessian_elements[1].gradient1[j] *= 2;
            if((*image_flags1 & REFL_X_FLAGS) && (*image_flags0) && (REFLECTIVE_X == 3))
              {
                hessian_elements[2].gradient1[j] *= 2;
                hessian_elements[3].gradient1[j] *= 2;
              }
          }
#endif
#if defined(REFLECTIVE_Y)
      if((*image_flags1 & REFL_Y_FLAGS) && !(*image_flags1 & OUTFLOW_Y))
        for(j = 0; j < 3; j++)
          {
            hessian_elements[2].gradient1[j] *= 2;
            if((*image_flags1 & REFL_Y_FLAGS) && (*image_flags0) && (REFLECTIVE_Y == 3))
              {
                hessian_elements[1].gradient1[j] *= 2;
                hessian_elements[3].gradient1[j] *= 2;
              }
          }
#endif
#if defined(REFLECTIVE_Z)
      if((*image_flags1 & REFL_Z_FLAGS) && !(*image_flags1 & OUTFLOW_Z))
        for(j = 0; j < 3; j++)
          {
            hessian_elements[3].gradient1[j] *= 2;
            if((*image_flags1 & REFL_Z_FLAGS) && (*image_flags0) && (REFLECTIVE_Z == 3))
              {
                hessian_elements[1].gradient1[j] *= 2;
                hessian_elements[2].gradient1[j] *= 2;
              }
          }
#endif


#ifdef SPECIAL_BOUNDARY
      MyIDType id1, id0;

      int j, l, k;

      id0 = p1->ID;
      id1 = p0->ID;

      if((id1 == -1 && id0 == -2) || (id0 == -1 && id1 == -2))
        {

          /*if id0 is in the surface side and id1 is in the fluid side */
          if(id1 == -1 && id0 == -2)
            {
              for(k = 0; k < N_Hess; k++)
                {
                  if(k != 1 && k != 2 && k != 3)
                    for(l = 0; l < 3; l++)
                      hessian_elements[k].gradient0[l] = hessian_elements[k].gradient1[l];

                  /*NOTE: The following is only for Hessian estimation, not for linear reconstruction! */

                  if(k == 1 || k == 2 || k == 3)
                    for(l = 0; l < 3; l++)      /*no-slip condition. Extrapolate all velocity gradients */
                      {
#ifdef VISCOSITY
                        hessian_elements[k].gradient0[l] = 2 * hessian_elements[k].gradient1[l];
#else
                        hessian_elements[k].gradient0[l] = hessian_elements[k].gradient1[l]
                          + hessian_elements[k].gradient1[0] * n[0] + hessian_elements[k].gradient1[1] * n[1] + hessian_elements[k].gradient1[2] * n[2];
#endif
                      }
                }

            }
          else
            {
              for(k = 0; k < N_Hess; k++)
                {
                  if(k != 1 && k != 2 && k != 3)
                    for(l = 0; l < 3; l++)
                      hessian_elements[k].gradient1[l] = hessian_elements[k].gradient0[l];

                  /*NOTE: The following is only for Hessian estimation, not for linear reconstruction! */

                  if(k == 1 || k == 2 || k == 3)
                    for(l = 0; l < 3; l++)      /*no-slip condition. Extrapolate all velocity gradients */
                      {
#ifdef VISCOSITY
                        hessian_elements[k].gradient1[l] = 2 * hessian_elements[k].gradient0[l];
#else
                        /* extrapolate only the normal component of the gradient */
                        hessian_elements[k].gradient1[l] = hessian_elements[k].gradient0[l]
                          + hessian_elements[k].gradient0[0] * n[0] + hessian_elements[k].gradient0[1] * n[1] + hessian_elements[k].gradient0[2] * n[2];
#endif
                      }


                }

            }
        }
#endif




      /* if the cell q0 is a local particle, construct the Hessian estimate, and the minmax values */
      if(p0->task == ThisTask && q0 >= 0 && q0 < NumGas)
        {
          if(TimeBinSynchronized[P[q0].TimeBinHydro])
            {
              fac0 = 0.5 * VF[i].area / SphP[q0].Volume;
              facA = VF[i].area / (nn * SphP[q0].Volume);

              for(k = 0; k < N_Hess; k++)
                {
                  data = (double *) (((char *) (&(SphP[q0].Hessian))) + hessian_elements[k].offset_hessian);
                  for(l = 0; l < 3; l++)
                    {
                      for(j = l; j < 3; j++)
                        {
                          direct = fac0 * n[j] * (hessian_elements[k].gradient0[l] +
                                                  hessian_elements[k].gradient1[l]) + facA * c[j] * (hessian_elements[k].gradient0[l] - hessian_elements[k].gradient1[l]);

                          transpose = fac0 * n[l] * (hessian_elements[k].gradient0[j] +
                                                     hessian_elements[k].gradient1[j]) + facA * c[l] * (hessian_elements[k].gradient0[j] - hessian_elements[k].gradient1[j]);


                          data[l * 3 + j] += (direct + transpose) / 2.0;

                        }
                    }

#ifdef SLOPE_LIMIT_HESSIANS
                  if(VF[i].area > 1.0e-5 * SphP[q1].SurfaceArea)
                    for(l = 0; l < 3; l++)
                      {
                        if(hessian_elements[k].max_value[3 * q0 + l] < hessian_elements[k].gradient0[l])
                          hessian_elements[k].max_value[3 * q0 + l] = hessian_elements[k].gradient0[l];
                        if(hessian_elements[k].min_value[3 * q0 + l] > hessian_elements[k].gradient0[l])
                          hessian_elements[k].min_value[3 * q0 + l] = hessian_elements[k].gradient0[l];
                      }
#endif

                }
            }
        }

      /* if the cell q1 is a local particle, construct the gradient estimate, and the minmax values */
      if(p1->task == ThisTask && q1 >= 0 && q1 < NumGas)
        {
          if(TimeBinSynchronized[P[q1].TimeBinHydro])
            {
              fac1 = -0.5 * VF[i].area / SphP[q1].Volume;
              facB = VF[i].area / (nn * SphP[q1].Volume);

              for(k = 0; k < N_Hess; k++)
                {
                  data = (double *) (((char *) (&(SphP[q1].Hessian))) + hessian_elements[k].offset_hessian);
                  for(l = 0; l < 3; l++)
                    {
                      for(j = l; j < 3; j++)
                        {
                          direct = fac1 * n[j] * (hessian_elements[k].gradient1[l] +
                                                  hessian_elements[k].gradient0[l]) + facB * c[j] * (hessian_elements[k].gradient1[l] - hessian_elements[k].gradient0[l]);

                          transpose = fac1 * n[l] * (hessian_elements[k].gradient1[j] +
                                                     hessian_elements[k].gradient0[j]) + facB * c[l] * (hessian_elements[k].gradient1[j] - hessian_elements[k].gradient0[j]);


                          data[l * 3 + j] += (direct + transpose) / 2.0;

                        }
                    }

#ifdef SLOPE_LIMIT_HESSIANS
                  if(VF[i].area > 1.0e-5 * SphP[q1].SurfaceArea)
                    for(l = 0; l < 3; l++)
                      {
                        if(hessian_elements[k].max_value[3 * q1 + l] < hessian_elements[k].gradient1[l])
                          hessian_elements[k].max_value[3 * q1 + l] = hessian_elements[k].gradient1[l];
                        if(hessian_elements[k].min_value[3 * q1 + l] > hessian_elements[k].gradient1[l])
                          hessian_elements[k].min_value[3 * q1 + l] = hessian_elements[k].gradient1[l];
                      }
#endif

                }
            }
        }
    }



  /* let's now implement a slope limitation if appropriate */
#ifdef SLOPE_LIMIT_HESSIANS
  for(i = 0; i < T->Nvf; i++)
    {
      point *p;
      int q;

      if(DP[VF[i].p1].index < 0 || DP[VF[i].p2].index < 0)
        continue;

      for(j = 0; j < 2; j++)
        {
          if(j == 0)
            p = &DP[VF[i].p1];
          else
            p = &DP[VF[i].p2];

          if(p->task == ThisTask && p->index >= 0 && p->index < NumGas)
            {
              q = p->index;
              if(TimeBinSynchronized[P[q].TimeBinHydro])
                {
                  d[0] = VF[i].cx - SphP[q].Center[0];
                  d[1] = VF[i].cy - SphP[q].Center[1];
                  d[2] = VF[i].cz - SphP[q].Center[2];

#ifdef PERIODIC
#if !defined(REFLECTIVE_X)
                  if(d[0] < -boxHalf_X)
                    d[0] += boxSize_X;
                  if(d[0] > boxHalf_X)
                    d[0] -= boxSize_X;
#endif
#if !defined(REFLECTIVE_Y)
                  if(d[1] < -boxHalf_Y)
                    d[1] += boxSize_Y;
                  if(d[1] > boxHalf_Y)
                    d[1] -= boxSize_Y;
#endif
#if !defined(REFLECTIVE_Z)
                  if(d[2] < -boxHalf_Z)
                    d[2] += boxSize_Z;
                  if(d[2] > boxHalf_Z)
                    d[2] -= boxSize_Z;
#endif
#endif

                  if(VF[i].area > 1.0e-5 * SphP[q].SurfaceArea)
                    {
                      for(k = 0; k < N_Hess; k++)
                        {
                          for(l = 0; l < 3; l++)
                            {
                              MySingle *gradient = (MySingle *) (((char *) (&(SphP[q].Grad))) + hessian_elements[k].offset_grad);
                              value = gradient[l];

                              data = (double *) (((char *) (&(SphP[q].Hessian))) + hessian_elements[k].offset_hessian);
                              data += 3 * l;

                              limit_gradient(d, value, hessian_elements[k].min_value[3 * q + l], hessian_elements[k].max_value[3 * q + l], data);
                            }
                        }
                    }
                }
            }
        }
    }

  for(k = N_Hess - 1; k >= 0; k--)
    {
      myfree(hessian_elements[k].max_value);
      myfree(hessian_elements[k].min_value);
    }
#endif

  CPU_Step[CPU_HESSIAN] += measure_time();
#endif
}
