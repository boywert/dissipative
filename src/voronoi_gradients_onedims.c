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

#include "allvars.h"
#include "proto.h"

#if defined(VORONOI) || defined(AMR)

#if defined(ONEDIMS)

#ifdef OUTPUT_DIVVEL
static void compute_divvel();
#endif

double getValue(int i, int k)
{
  if((grad_elements[k].type == GRADIENT_TYPE_VELX) || (grad_elements[k].type == GRADIENT_TYPE_VELY) || (grad_elements[k].type == GRADIENT_TYPE_VELZ))
    return *(MyFloat *) (((char *) (&P[i])) + grad_elements[k].offset);
  else
    return *(MyFloat *) (((char *) (&SphP[i])) + grad_elements[k].offset);
}

void calculate_gradients(void)
{
  CPU_Step[CPU_MISC] += measure_time();

  printf("Calculating 1D gradients...\n");

  int idx, i, k;
#pragma omp parallel for private(idx,i,k)
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      for(k = 0; k < N_Grad; k++)
        {
          double Value = getValue(i, k);
          double Pos = P[i].Pos[0];

#if defined (ONEDIMS_SPHERICAL) || defined (REFLECTIVE_X)
          if(i == 0 || i == NumGas - 1)
            {
              MySingle *data = (MySingle *) (((char *) (&(SphP[i].Grad))) + grad_elements[k].offset_grad);
              memset(data, 0, 3 * sizeof(MySingle));
#ifdef TVD_SLOPE_LIMITER
              data = (MySingle *) (((char *) (&(SphP[i].GradUl))) + grad_elements[k].offset_grad);
              memset(data, 0, 3 * sizeof(MySingle));
#endif
              continue;
            }
#endif
          /* if we get here, we have periodic boundary conditions or are not at the boundaries */
          double ValueL, ValueR;

          if(i == 0)
            ValueL = getValue(NumGas - 1, k);
          else
            ValueL = getValue(i - 1, k);

          if(i == NumGas - 1)
            ValueR = getValue(0, k);
          else
            ValueR = getValue(i + 1, k);

          double PosL = Mesh.DP[i - 1].x;
          double PosR = Mesh.DP[i + 1].x;

          //double grad = ((Pos - PosL) * (Value - ValueL) + (Pos - PosR) * (Value - ValueR)) / ((Pos - PosL) * (Pos - PosL) + (Pos - PosR) * (Pos - PosR));
          double grad = (ValueL - ValueR) / (PosL - PosR);

          MySingle *data = (MySingle *) (((char *) (&(SphP[i].Grad))) + grad_elements[k].offset_grad);
          data[0] = grad;
          data[1] = 0;
          data[2] = 0;

#ifdef TVD_SLOPE_LIMITER
          MySingle *dataUL = (MySingle *) (((char *) (&(SphP[i].GradUl))) + grad_elements[k].offset_grad);

          dataUL[0] = grad;
          dataUL[1] = data[1];
          dataUL[2] = data[2];
#endif

#ifdef NUCLEAR_NETWORK
          if(grad_elements[k].type == GRADIENT_TYPE_DENSITY)
            {
              int j;
              for(j = 0; j < 3; j++)
                SphP[i].Grad.drhoU[j] = SphP[i].Grad.drho[j];
            }
          
          if(grad_elements[k].type == GRADIENT_TYPE_PRESSURE)
            {
              int j;
              for(j = 0; j < 3; j++)
                SphP[i].Grad.dpressU[j] = SphP[i].Grad.dpress[j];
            }
#endif

#ifndef UNLIMITED_GRADIENTS
          double ValueMin = dmin(ValueL, ValueR);
          double ValueMax = dmax(ValueL, ValueR);

          if(Value + grad * (PosL - Pos) < ValueMin)
            {
              if(ValueMin < Value)
                grad = (ValueMin - Value) / (PosL - Pos);
              else
                grad = 0.;
            }

          if(Value + grad * (PosL - Pos) > ValueMax)
            {
              if(ValueMax > Value)
                grad = (ValueMax - Value) / (PosL - Pos);
              else
                grad = 0.;
            }

          if(Value + grad * (PosR - Pos) < ValueMin)
            {
              if(ValueMin < Value)
                grad = (ValueMin - Value) / (PosR - Pos);
              else
                grad = 0.;
            }

          if(Value + grad * (PosR - Pos) > ValueMax)
            {
              if(ValueMax > Value)
                grad = (ValueMax - Value) / (PosR - Pos);
              else
                grad = 0.;
            }

          data[0] = grad;
#endif
        }
    }

#ifdef OUTPUT_DIVVEL
  compute_divvel();
#endif

  CPU_Step[CPU_GRADIENTS] += measure_time();
}

#ifdef OUTPUT_DIVVEL
void compute_divvel()
{
  face *VF = Mesh.VF;
  double VelxL, VelxR;

  int idx, i;
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(i == 0)
        {
#if defined (ONEDIMS_SPHERICAL) || defined (REFLECTIVE_X)
          VelxL = P[i].Vel[0];
#else
          VelxL = P[NumGas - 1].Vel[0];
#endif
        }
      else
        VelxL = P[i - 1].Vel[0];

      if(i == NumGas - 1)
        {
#if defined (ONEDIMS_SPHERICAL) || defined (REFLECTIVE_X)
          VelxR = P[i].Vel[0];
#else
          VelxR = P[0].Vel[0];
#endif
        }
      else
        VelxR = P[i + 1].Vel[0];

      SphP[i].DivVel = 0.5 * (VF[i].area * VelxR - VF[i - 1].area * VelxL) / SphP[i].Volume;
    }
}
#endif

#endif

#endif
