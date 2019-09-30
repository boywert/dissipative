/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/rt/rt_gradients.c
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

#include "../allvars.h"
#include "../proto.h"
#include "../voronoi.h"

#ifdef RT_ADVECT

static int RT_N_Grad;
#define GRADIENT_TYPE_PHOTONS  1000

static struct rt_grad_elements
{
  size_t offset;                /* offset of the quantity in the SphP struct */
  size_t offset_exch;           /* offset of the quantity in the PrimExch struct */
  size_t offset_grad;           /* offset in the grad_data struct */
  double *min_value, *max_value;
  double value0, value1;
#if !defined(RT_HEALPIX_NSIDE)
  int sourceid0, sourceid1;
  size_t offset_sourceid, offset_exch_sourceid;
#endif
} rt_grad_elements[RT_MAXGRADIENTS];


void rt_init_gradients()
{
  int k;

#ifndef RT_HEALPIX_NSIDE
  for(k = 0; k < RT_N_DIR; k++)
    rt_gradient_init(&SphP[0].DensPhot[k], &RTPrimExch[0].DensPhot[k], SphP[0].rt_Grad.ddensphot[k], &SphP[0].SourceID[k], &RTPrimExch[0].SourceID[k], GRADIENT_TYPE_PHOTONS + k);
#else
  for(k = 0; k < RT_N_DIR; k++)
    rt_gradient_init(&SphP[0].DensPhot[k], &RTPrimExch[0].DensPhot[k], SphP[0].rt_Grad.ddensphot[k], NULL, NULL, GRADIENT_TYPE_PHOTONS + k);
#endif

  mpi_printf("INIT: %d/%d Gradients used for RT.\n", RT_N_Grad, RT_MAXGRADIENTS);
}

void rt_gradient_init(MyFloat * addr, MyFloat * addr_exch, double *addr_grad, int *addr_sourceid, int *addr_sourceid_exch, int type)
{
  if(RT_N_Grad == RT_MAXGRADIENTS)
    {
      mpi_printf("Failed to register gradient, maximum of %d already reached\n", MAXGRADIENTS);
      terminate("MAXGRADIENTS reached");
    }

  /* basic structure is SphP */
  rt_grad_elements[RT_N_Grad].offset = ((char *) addr) - ((char *) &SphP[0]);

  rt_grad_elements[RT_N_Grad].offset_exch = ((char *) addr_exch) - ((char *) &RTPrimExch[0]);
  rt_grad_elements[RT_N_Grad].offset_grad = ((char *) addr_grad) - ((char *) &(SphP[0].rt_Grad));

#if !defined(RT_HEALPIX_NSIDE)
  rt_grad_elements[RT_N_Grad].offset_sourceid = ((char *) addr_sourceid) - ((char *) &SphP[0]);
  rt_grad_elements[RT_N_Grad].offset_exch_sourceid = ((char *) addr_sourceid_exch) - ((char *) &RTPrimExch[0]);
#endif

  RT_N_Grad++;
}

void rt_calculate_green_gauss_gradients(void)
{
  CPU_Step[CPU_MISC] += measure_time();

  point *DP = Mesh.DP;
  face *VF = Mesh.VF;

  int idx, i, j, k, t0, t1, q0, q1, s0, s1;
  MyFloat value;
  double fac0, fac1, n[3], nn, d[3], *data;

  double facA, facB, c[3];

  for(k = 0; k < RT_N_Grad; k++)
    {
      rt_grad_elements[k].min_value = mymalloc("rt_gradmin", NumGas * sizeof(double));
      rt_grad_elements[k].max_value = mymalloc("rt_gradmax", NumGas * sizeof(double));
    }

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      for(k = 0; k < RT_N_Grad; k++)
        {
          rt_grad_elements[k].min_value[i] = +MAX_REAL_NUMBER;
          rt_grad_elements[k].max_value[i] = -MAX_REAL_NUMBER;
        }

      for(k = 0; k < RT_N_Grad; k++)
        {
          data = (double *) (((char *) (&(SphP[i].rt_Grad))) + rt_grad_elements[k].offset_grad);
          for(j = 0; j < 3; j++)
            data[j] = 0;
        }
    }

  for(i = 0; i < Mesh.Nvf; i++)
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

      /* the other side */
      q1 = p1->index;
      s1 = p0->index;
      t1 = p0->task;

      for(k = 0; k < RT_N_Grad; k++)
        {
          rt_grad_elements[k].value0 = 0;
          rt_grad_elements[k].value1 = 0;
        }


      /* let's get the physical quantities of interest for the particle one one of the sides */
      if(t0 == ThisTask)
        {
          if(s0 >= NumGas)
            s0 -= NumGas;

          if(s0 >= 0)
            {
              for(k = 0; k < RT_N_Grad; k++)
                {
                  rt_grad_elements[k].value0 = *(MyFloat *) (((char *) (&SphP[s0])) + rt_grad_elements[k].offset);

#if !defined(RT_HEALPIX_NSIDE)
                  rt_grad_elements[k].sourceid0 = *(int *) (((char *) (&SphP[s0])) + rt_grad_elements[k].offset_sourceid);
#endif
                }
            }
        }
      else
        {
          for(k = 0; k < RT_N_Grad; k++)
            {
              rt_grad_elements[k].value0 = *(MyFloat *) (((char *) (&RTPrimExch[s0])) + rt_grad_elements[k].offset_exch);

#if !defined(RT_HEALPIX_NSIDE)
              rt_grad_elements[k].sourceid0 = *(int *) (((char *) (&RTPrimExch[s0])) + rt_grad_elements[k].offset_exch_sourceid);
#endif
            }

        }

      /* let's get the physical quantities of interest for the particle on the other side */
      if(t1 == ThisTask)
        {
          if(s1 >= NumGas)
            s1 -= NumGas;

          if(s1 >= 0)
            {
              for(k = 0; k < RT_N_Grad; k++)
                {
                  rt_grad_elements[k].value1 = *(MyFloat *) (((char *) (&SphP[s1])) + rt_grad_elements[k].offset);

#if !defined(RT_HEALPIX_NSIDE)
                  rt_grad_elements[k].sourceid1 = *(int *) (((char *) (&SphP[s1])) + rt_grad_elements[k].offset_sourceid);
#endif
                }
            }
        }
      else
        {
          for(k = 0; k < RT_N_Grad; k++)
            {
              rt_grad_elements[k].value1 = *(MyFloat *) (((char *) (&RTPrimExch[s1])) + rt_grad_elements[k].offset_exch);

#if !defined(RT_HEALPIX_NSIDE)
              rt_grad_elements[k].sourceid1 = *(int *) (((char *) (&RTPrimExch[s1])) + rt_grad_elements[k].offset_exch_sourceid);
#endif
            }
        }

      /* if the cell q0 is a local particle, construct the gradient estimate, and the minmax values */
      if(p0->task == ThisTask && q0 >= 0 && q0 < NumGas)
        {
          if(TimeBinSynchronized[P[q0].TimeBinHydro])
            {
              fac0 = 0.5 * VF[i].area / SphP[q0].Volume;
              facA = VF[i].area / (nn * SphP[q0].Volume);

              for(k = 0; k < RT_N_Grad; k++)
                {
                  double value0;

                  value0 = rt_grad_elements[k].value0;

                  data = (double *) (((char *) (&(SphP[q0].rt_Grad))) + rt_grad_elements[k].offset_grad);

#if !defined(RT_HEALPIX_NSIDE)
                  /* here we need to see whether we find matching source IDs */
                  int sourceid0 = rt_grad_elements[k].sourceid0;
                  int sourceid1 = rt_grad_elements[k].sourceid1;

                  if(sourceid0 != sourceid1 && !(sourceid0 == 1000 || sourceid1 == 1000))
                    {
                      value0 = rt_grad_elements[k].value1;      /* use this if no match is found */
                      int dir, kk = k;
                      if(sourceid0 > sourceid1)
                        dir = -1;
                      else
                        dir = 1;
                      while(sourceid0 != sourceid1)
                        {
                          kk += dir;
                          if(kk < 0 || kk >= RT_N_Grad)
                            break;
                          sourceid0 = rt_grad_elements[kk].sourceid0;
                          if(sourceid0 == sourceid1)
                            {
                              value0 = rt_grad_elements[kk].value0;
                              break;
                            }
                        }
                    }
#endif
                  for(j = 0; j < 3; j++)
                    {
                      data[j] += fac0 * n[j] * value0;
                      data[j] += facA * c[j] * (value0 - rt_grad_elements[k].value1);
                    }

                  if(VF[i].area > 1.0e-5 * SphP[q0].SurfaceArea)
                    {
                      if(rt_grad_elements[k].max_value[q0] < value0)
                        rt_grad_elements[k].max_value[q0] = value0;

                      if(rt_grad_elements[k].min_value[q0] > value0)
                        rt_grad_elements[k].min_value[q0] = value0;
                    }
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

              for(k = 0; k < RT_N_Grad; k++)
                {
                  double value1;

                  value1 = rt_grad_elements[k].value1;

#if !defined(RT_HEALPIX_NSIDE)
                  /* here we need to see whether we find matching source IDs */
                  int sourceid0 = rt_grad_elements[k].sourceid0;
                  int sourceid1 = rt_grad_elements[k].sourceid1;

                  if(sourceid0 != sourceid1 && !(sourceid0 == 1000 || sourceid1 == 1000))
                    {
                      value1 = rt_grad_elements[k].value0;      /* use this if no match is found */
                      int dir, kk = k;
                      if(sourceid1 > sourceid0)
                        dir = -1;
                      else
                        dir = 1;
                      while(sourceid0 != sourceid1)
                        {
                          kk += dir;
                          if(kk < 0 || kk >= RT_N_Grad)
                            break;
                          sourceid1 = rt_grad_elements[kk].sourceid1;
                          if(sourceid0 == sourceid1)
                            {
                              value1 = rt_grad_elements[kk].value1;
                              break;
                            }
                        }
                    }
#endif
                  data = (double *) (((char *) (&(SphP[q1].rt_Grad))) + rt_grad_elements[k].offset_grad);
                  for(j = 0; j < 3; j++)
                    {
                      data[j] += fac1 * n[j] * value1;
                      data[j] += facB * c[j] * (value1 - rt_grad_elements[k].value0);
                    }

                  if(VF[i].area > 1.0e-5 * SphP[q1].SurfaceArea)
                    {
                      if(rt_grad_elements[k].max_value[q1] < value1)
                        rt_grad_elements[k].max_value[q1] = value1;

                      if(rt_grad_elements[k].min_value[q1] > value1)
                        rt_grad_elements[k].min_value[q1] = value1;
                    }
                }
            }
        }
    }

#ifdef RT_HEALPIX_NSIDE
  /* calculate a gradient for the total radiation field (without slope limiting), used for determining the direction of radiation transport */
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      for(k = 0; k < 3; k++)
        SphP[i].rt_Grad.ddensphot_unlimited[k] = 0;
      for(j = 0; j < RT_N_DIR; j++)
        {
          for(k = 0; k < 3; k++)
            SphP[i].rt_Grad.ddensphot_unlimited[k] += SphP[i].rt_Grad.ddensphot[j][k];
        }
    }
#endif

  /* let's now implement a slope limitation if appropriate */

  for(i = 0; i < Mesh.Nvf; i++)
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
                  if(d[0] < -boxHalf_X)
                    d[0] += boxSize_X;
                  if(d[0] > boxHalf_X)
                    d[0] -= boxSize_X;

                  if(d[1] < -boxHalf_Y)
                    d[1] += boxSize_Y;
                  if(d[1] > boxHalf_Y)
                    d[1] -= boxSize_Y;

                  if(d[2] < -boxHalf_Z)
                    d[2] += boxSize_Z;
                  if(d[2] > boxHalf_Z)
                    d[2] -= boxSize_Z;
#endif

                  if(VF[i].area > 1.0e-5 * SphP[q].SurfaceArea)
                    {
                      for(k = 0; k < RT_N_Grad; k++)
                        {
                          value = *(MyFloat *) (((char *) (&SphP[q])) + rt_grad_elements[k].offset);

                          data = (double *) (((char *) (&(SphP[q].rt_Grad))) + rt_grad_elements[k].offset_grad);

                          limit_gradient(d, value, rt_grad_elements[k].min_value[q], rt_grad_elements[k].max_value[q], data);

                        }
                    }
                }
            }
        }
    }

  for(k = RT_N_Grad - 1; k >= 0; k--)
    {
      myfree(rt_grad_elements[k].max_value);
      myfree(rt_grad_elements[k].min_value);
    }

  CPU_Step[CPU_GRADIENTS] += measure_time();
}

#endif
