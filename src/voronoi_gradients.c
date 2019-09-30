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

#if defined(VORONOI) || (defined(AMR) && !defined(AMR_GRADIENTS))

#if defined(GRADIENTS_GREEN_GAUSS) && !defined(ONEDIMS)

static double atime, hubble_a;

void calculate_gradients(void)
{
#ifdef DG
  return;
#endif

  CPU_Step[CPU_MISC] += measure_time();

  point *DP = Mesh.DP;
  face *VF = Mesh.VF;

  int idx, i, j, k, t0, t1, q0, q1, s0, s1;
  MyFloat value;
  double fac0, fac1, n[3], nn, d[3], *curl;
  MySingle *data;

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
  double M0 = 0, M1 = 0;
#endif
  double facA, facB, c[3];

#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
  unsigned int *image_flags0, *image_flags1;
#endif

  if(All.ComovingIntegrationOn)
    {
      atime = All.Time;
      hubble_a = hubble_function(atime);
    }
  else
    {
      atime = 1.0;
      hubble_a = 0.0;
    }

  for(k = 0; k < N_Grad; k++)
    {
      grad_elements[k].min_value = mymalloc("gradmin", NumGas * sizeof(double));
      grad_elements[k].max_value = mymalloc("gradmax", NumGas * sizeof(double));
    }

  curl = mymalloc("curl", NumGas * 3 * sizeof(double));

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;
#if defined (OUTPUT_DIVVEL) || defined(TGCHEM) || defined(SGCHEM)
      SphP[i].DivVel = 0;
#endif

#ifdef DEREFINE_GENTLY
      SphP[i].DoNotDerefFlag = 0;
#endif

      for(k = 0; k < N_Grad; k++)
        {
          grad_elements[k].min_value[i] = +MAX_REAL_NUMBER;
          grad_elements[k].max_value[i] = -MAX_REAL_NUMBER;
        }

#ifdef SPECIAL_BOUNDARY
      SphP[i].MinDistBoundaryCell = +MAX_REAL_NUMBER;
#endif

      for(k = 0; k < N_Grad; k++)
        {
          data = (MySingle *) (((char *) (&(SphP[i].Grad))) + grad_elements[k].offset_grad);
          for(j = 0; j < 3; j++)
            data[j] = 0;
        }

      for(j = 0; j < 3; j++)
        curl[3 * i + j] = 0;
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

      for(k = 0; k < N_Grad; k++)
        {
          grad_elements[k].value0 = 0;
          grad_elements[k].value1 = 0;
        }

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
      M0 = M1 = 0;
#endif

      /* let's get the physical quantities of interest for the particle on one of the sides */
      if(t0 == ThisTask)
        {
          if(s0 >= NumGas)
            s0 -= NumGas;

          if(s0 >= 0)
            {
              for(k = 0; k < N_Grad; k++)
                {
                  if((grad_elements[k].type == GRADIENT_TYPE_VELX) || (grad_elements[k].type == GRADIENT_TYPE_VELY) || (grad_elements[k].type == GRADIENT_TYPE_VELZ))
                    grad_elements[k].value0 = *(MyFloat *) (((char *) (&P[s0])) + grad_elements[k].offset);
                  else
                    grad_elements[k].value0 = *(MyFloat *) (((char *) (&SphP[s0])) + grad_elements[k].offset);
                }
#ifdef USE_ENTROPY_FOR_COLD_FLOWS
              M0 = P[s0].Mass;
#endif
            }
        }
      else
        {
          for(k = 0; k < N_Grad; k++)
            {
              grad_elements[k].value0 = *(MyFloat *) (((char *) (&PrimExch[s0])) + grad_elements[k].offset_exch);
            }

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
          M0 = PrimExch[s0].Mass;
#endif
        }

      /* let's get the physical quantities of interest for the particle on the other side */
      if(t1 == ThisTask)
        {
          if(s1 >= NumGas)
            s1 -= NumGas;

          if(s1 >= 0)
            {
              for(k = 0; k < N_Grad; k++)
                {
                  if((grad_elements[k].type == GRADIENT_TYPE_VELX) || (grad_elements[k].type == GRADIENT_TYPE_VELY) || (grad_elements[k].type == GRADIENT_TYPE_VELZ))
                    grad_elements[k].value1 = *(MyFloat *) (((char *) (&P[s1])) + grad_elements[k].offset);
                  else
                    grad_elements[k].value1 = *(MyFloat *) (((char *) (&SphP[s1])) + grad_elements[k].offset);

                }
#ifdef USE_ENTROPY_FOR_COLD_FLOWS
              M1 = P[s1].Mass;
#endif
            }
        }
      else
        {
          for(k = 0; k < N_Grad; k++)
            {
              grad_elements[k].value1 = *(MyFloat *) (((char *) (&PrimExch[s1])) + grad_elements[k].offset_exch);
            }

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
          M1 = PrimExch[s1].Mass;
#endif
        }

      /* convert velocity components to peculiar velocities */
      if(All.ComovingIntegrationOn)
        {
          GVelx->value0 /= atime;
          GVelx->value1 /= atime;
          GVely->value0 /= atime;
          GVely->value1 /= atime;
          GVelz->value0 /= atime;
          GVelz->value1 /= atime;
        }

#if defined(REFLECTIVE_X)
      if((*image_flags0 & REFL_X_FLAGS) && !(*image_flags0 & OUTFLOW_X))
        GVelx->value0 *= -1;
      if((*image_flags0 & REFL_X_FLAGS) && (*image_flags0) && (REFLECTIVE_X == 3))
        {
          GVely->value0 *= -1;
          GVelz->value0 *= -1;
        }
#endif
#if defined(REFLECTIVE_Y)
      if((*image_flags0 & REFL_Y_FLAGS) && !(*image_flags0 & OUTFLOW_Y))
        GVely->value0 *= -1;
      if((*image_flags0 & REFL_Y_FLAGS) && (*image_flags0) && (REFLECTIVE_Y == 3))
        {
          GVelx->value0 *= -1;
          GVelz->value0 *= -1;
        }
#endif
#if defined(REFLECTIVE_Z)
      if((*image_flags0 & REFL_Z_FLAGS) && !(*image_flags0 & OUTFLOW_Z))
        GVelz->value0 *= -1;
      if((*image_flags0 & REFL_Z_FLAGS) && (*image_flags0) && (REFLECTIVE_Z == 3))
        {
          GVelx->value0 *= -1;
          GVely->value0 *= -1;
        }
#endif

#if defined(REFLECTIVE_X)
      if((*image_flags1 & REFL_X_FLAGS) && !(*image_flags1 & OUTFLOW_X))
        GVelx->value1 *= -1;
      if((*image_flags1 & REFL_X_FLAGS) && (*image_flags0) && (REFLECTIVE_X == 3))
        {
          GVely->value1 *= -1;
          GVelz->value1 *= -1;
        }
#endif
#if defined(REFLECTIVE_Y)
      if((*image_flags1 & REFL_Y_FLAGS) && !(*image_flags1 & OUTFLOW_Y))
        GVely->value1 *= -1;
      if((*image_flags1 & REFL_Y_FLAGS) && (*image_flags0) && (REFLECTIVE_Y == 3))
        {
          GVelx->value1 *= -1;
          GVelz->value1 *= -1;
        }
#endif
#if defined(REFLECTIVE_Z)
      if((*image_flags1 & REFL_Z_FLAGS) && !(*image_flags1 & OUTFLOW_Z))
        GVelz->value1 *= -1;
      if((*image_flags1 & REFL_Z_FLAGS) && (*image_flags0) && (REFLECTIVE_Z == 3))
        {
          GVelx->value1 *= -1;
          GVely->value1 *= -1;
        }
#endif



#if defined(BOUNDARY_REFL_FLUIDSIDE_MINID) && defined(BOUNDARY_REFL_FLUIDSIDE_MAXID) && defined(BOUNDARY_REFL_SOLIDSIDE_MINID) && defined(BOUNDARY_REFL_SOLIDSIDE_MAXID)
      {
        MyIDType id1, id0;

        double m[3], vx, vy, vz, vp0, vp1;

        double velGas0[3], velGas1[3];
        double dt;

        int j = 0;

        id0 = p1->ID;
        id1 = p0->ID;

        if((id0 >= BOUNDARY_REFL_FLUIDSIDE_MINID && id0 < BOUNDARY_REFL_FLUIDSIDE_MAXID && id1 >= BOUNDARY_REFL_SOLIDSIDE_MINID
            && id1 < BOUNDARY_REFL_SOLIDSIDE_MAXID) || (id1 >= BOUNDARY_REFL_FLUIDSIDE_MINID && id1 < BOUNDARY_REFL_FLUIDSIDE_MAXID
                                                        && id0 >= BOUNDARY_REFL_SOLIDSIDE_MINID && id0 < BOUNDARY_REFL_SOLIDSIDE_MAXID))
          {
            m[0] = 0.5 * (p1->x + p0->x);
            m[1] = 0.5 * (p1->y + p0->y);
            m[2] = 0.5 * (p1->z + p0->z);

            velGas0[0] = GVelx->value0;
            velGas0[1] = GVely->value0;
            velGas0[2] = GVelz->value0;

            velGas1[0] = GVelx->value1;
            velGas1[1] = GVely->value1;
            velGas1[2] = GVelz->value1;

            /*we project against the local normal of the solid/artificial
               surface that points from p0(index s1) to p1(index s0) */
            vp0 = velGas0[0] * n[0] + velGas0[1] * n[1] + velGas0[2] * n[2];
            vp1 = velGas1[0] * n[0] + velGas1[1] * n[1] + velGas1[2] * n[2];

            /*if id0 is in the solid side and id1 is in the fluid side */
            if(id1 >= BOUNDARY_REFL_FLUIDSIDE_MINID && id1 < BOUNDARY_REFL_FLUIDSIDE_MAXID && id0 >= BOUNDARY_REFL_SOLIDSIDE_MINID && id0 < BOUNDARY_REFL_SOLIDSIDE_MAXID)
              {
                for(k = 0; k < N_Grad; k++)
                  {
                    if(grad_elements[k].type != GRADIENT_TYPE_VELX && grad_elements[k].type != GRADIENT_TYPE_VELY && grad_elements[k].type != GRADIENT_TYPE_VELZ)
                      grad_elements[k].value0 = grad_elements[k].value1;
                  }

                velGas0[0] = velGas1[0] - 2 * vp1 * n[0];
                velGas0[1] = velGas1[1] - 2 * vp1 * n[1];
                velGas0[2] = velGas1[2] - 2 * vp1 * n[2];
              }
            else
              {
                for(k = 0; k < N_Grad; k++)
                  {
                    if(grad_elements[k].type != GRADIENT_TYPE_VELX && grad_elements[k].type != GRADIENT_TYPE_VELY && grad_elements[k].type != GRADIENT_TYPE_VELZ)
                      grad_elements[k].value1 = grad_elements[k].value0;
                  }

                velGas1[0] = velGas0[0] - 2 * vp0 * n[0];
                velGas1[1] = velGas0[1] - 2 * vp0 * n[1];
                velGas1[2] = velGas0[2] - 2 * vp0 * n[2];
              }

            /*we rewrite the physical quantities */

            GVelx->value0 = velGas0[0];
            GVely->value0 = velGas0[1];
            GVelz->value0 = velGas0[2];

            GVelx->value1 = velGas1[0];
            GVely->value1 = velGas1[1];
            GVelz->value1 = velGas1[2];
          }
      }
#endif




#ifdef SPECIAL_BOUNDARY

      MyIDType id1, id0;

      double m[3], vx, vy, vz, vp0, vp1;

      double velGas0[3], velGas1[3];
      double dt;

      int j = 0;

      id0 = p1->ID;
      id1 = p0->ID;

      if(id0 == -1 && id1 > 0)  /* id1 is a regular cell */
        {
          if(p0->task == ThisTask && p0->index >= 0 && p0->index < NumGas)
            {
              if((nn < SphP[p0->index].MinDistBoundaryCell) && (SphP[p0->index].MinDistBoundaryCell > 0))
                {
                  SphP[p0->index].MinDistBoundaryCell = nn;
                }
              if(get_cell_radius(p0->index) < 0.25 * get_cell_radius(p1->index))
                SphP[p0->index].MinDistBoundaryCell = 0;
            }
        }

      if(id1 == -1 && id0 > 0)  /* id0 is a regular cell */
        {
          if(p1->task == ThisTask && p1->index >= 0 && p1->index < NumGas)
            {
              if((nn < SphP[p1->index].MinDistBoundaryCell) && (SphP[p1->index].MinDistBoundaryCell > 0))
                {
                  SphP[p1->index].MinDistBoundaryCell = nn;
                }
              if(get_cell_radius(p1->index) < 0.25 * get_cell_radius(p0->index))
                SphP[p1->index].MinDistBoundaryCell = 0;
            }
        }
      if(id1 == -2 && id0 > 0)  /* id0 is already next to an ID=-2 particle */
        if(p1->task == ThisTask && p1->index >= 0 && p1->index < NumGas)
          SphP[p1->index].MinDistBoundaryCell = 0;
      if(id0 == -2 && id1 > 0)  /* id1 is already next to an ID=-2 particle */
        if(p0->task == ThisTask && p0->index >= 0 && p0->index < NumGas)
          SphP[p0->index].MinDistBoundaryCell = 0;



      if((id1 == -1 && id0 == -2) || (id0 == -1 && id1 == -2))
        {
          m[0] = 0.5 * (p1->x + p0->x);
          m[1] = 0.5 * (p1->y + p0->y);
          m[2] = 0.5 * (p1->z + p0->z);

          velGas0[0] = GVelx->value0;
          velGas0[1] = GVely->value0;
          velGas0[2] = GVelz->value0;

          velGas1[0] = GVelx->value1;
          velGas1[1] = GVely->value1;
          velGas1[2] = GVelz->value1;

          /* let's get the velocity of the boundary */
          dt = 0.5 * ((P[q0].TimeBinHydro ? (((integertime) 1) << P[q0].TimeBinHydro) : 0) * All.Timebase_interval +
                      (P[q1].TimeBinHydro ? (((integertime) 1) << P[q1].TimeBinHydro) : 0) * All.Timebase_interval);
          boundary_get_velocity(m[0], m[1], m[2], &vx, &vy, &vz, dt);

          /*and we boost to the frame of the moving contour */
          velGas0[0] -= vx;
          velGas0[1] -= vy;
          velGas0[2] -= vz;
          velGas1[0] -= vx;
          velGas1[1] -= vy;
          velGas1[2] -= vz;

          /*we project against the local normal of the solid/artificial
             surface that points from p0(index s1) to p1(index s0) */
          vp0 = velGas0[0] * n[0] + velGas0[1] * n[1] + velGas0[2] * n[2];
          vp1 = velGas1[0] * n[0] + velGas1[1] * n[1] + velGas1[2] * n[2];

          /*if id0 is in the surface side and id1 is in the fluid side */
          if(id1 == -1 && id0 == -2)
            {

              if(All.SpecialBoundaryType == 1)
                {
                  velGas0[0] = velGas1[0] - (vp0 + vp1) * n[0];
                  velGas0[1] = velGas1[1] - (vp0 + vp1) * n[1];
                  velGas0[2] = velGas1[2] - (vp0 + vp1) * n[2];

                }
              if(All.SpecialBoundaryType == 2)
                {
                  /*enforce no-slip condition (reflect all velocity components) */
                  velGas0[0] = -velGas1[0];
                  velGas0[1] = -velGas1[1];
                  velGas0[2] = -velGas1[2];
                }
              if(All.SpecialBoundaryType == 3 || All.SpecialBoundaryType == 4)
                {
                  velGas0[0] = velGas1[0];
                  velGas0[1] = velGas1[1];
                  velGas0[2] = velGas1[2];
                }
              for(k = 0; k < N_Grad; k++)
                {
                  if(grad_elements[k].type != GRADIENT_TYPE_VELX && grad_elements[k].type != GRADIENT_TYPE_VELY && grad_elements[k].type != GRADIENT_TYPE_VELZ)
                    grad_elements[k].value0 = grad_elements[k].value1;
                }
            }
          else
            {

              if(All.SpecialBoundaryType == 1)
                {
                  velGas1[0] = velGas0[0] - (vp0 + vp1) * n[0];
                  velGas1[1] = velGas0[1] - (vp0 + vp1) * n[1];
                  velGas1[2] = velGas0[2] - (vp0 + vp1) * n[2];
                }
              if(All.SpecialBoundaryType == 2)
                {
                  /*enforce no-slip condition */
                  velGas1[0] = -velGas0[0];
                  velGas1[1] = -velGas0[1];
                  velGas1[2] = -velGas0[2];
                }
              if(All.SpecialBoundaryType == 3 || All.SpecialBoundaryType == 4)
                {
                  velGas1[0] = velGas0[0];
                  velGas1[1] = velGas0[1];
                  velGas1[2] = velGas0[2];
                }

              for(k = 0; k < N_Grad; k++)
                {
                  if(grad_elements[k].type != GRADIENT_TYPE_VELX && grad_elements[k].type != GRADIENT_TYPE_VELY && grad_elements[k].type != GRADIENT_TYPE_VELZ)
                    grad_elements[k].value1 = grad_elements[k].value0;
                }

            }

          velGas0[0] += vx;
          velGas0[1] += vy;
          velGas0[2] += vz;
          velGas1[0] += vx;
          velGas1[1] += vy;
          velGas1[2] += vz;

          /*we rewrite the physical quantities */

          GVelx->value0 = velGas0[0];
          GVely->value0 = velGas0[1];
          GVelz->value0 = velGas0[2];

          GVelx->value1 = velGas1[0];
          GVely->value1 = velGas1[1];
          GVelz->value1 = velGas1[2];

        }
#endif






#ifdef COFFEE_PROBLEM
      MyIDType id1, id2;

      id1 = DP[VF[i].p1].ID;
      id2 = DP[VF[i].p2].ID;
      double m[2], vx, vy, vp0, vp1, velGas0[2], velGas1[2];

      if((id1 >= 10000000 && id1 < 20000000 && id2 >= 20000000) || (id2 >= 10000000 && id2 < 20000000 && id1 >= 20000000))
        {
          m[0] = 0.5 * (DP[VF[i].p1].x + DP[VF[i].p2].x);
          m[1] = 0.5 * (DP[VF[i].p1].y + DP[VF[i].p2].y);

          coffee_get_velocity(m[0], m[1], &vx, &vy);

          velGas0[0] = GVelx->value0;
          velGas0[1] = GVely->value0;

          velGas1[0] = GVelx->value1;
          velGas1[1] = GVely->value1;

          velGas0[0] -= vx;
          velGas0[1] -= vy;
          velGas1[0] -= vx;
          velGas1[1] -= vy;

          vp0 = velGas0[0] * n[0] + velGas0[1] * n[1];
          vp1 = velGas1[0] * n[0] + velGas1[1] * n[1];

          if(id2 >= 10000000 && id2 < 20000000 && id1 >= 20000000)
            {
              velGas0[0] = velGas1[0] - (vp0 + vp1) * n[0];
              velGas0[1] = velGas1[1] - (vp0 + vp1) * n[1];

              for(k = 0; k < N_Grad; k++)
                {
                  if(grad_elements[k].type != GRADIENT_TYPE_VELX && grad_elements[k].type != GRADIENT_TYPE_VELY)
                    grad_elements[k].value0 = grad_elements[k].value1;
                }
            }
          else
            {
              velGas1[0] = velGas0[0] - (vp0 + vp1) * n[0];
              velGas1[1] = velGas0[1] - (vp0 + vp1) * n[1];

              for(k = 0; k < N_Grad; k++)
                {
                  if(grad_elements[k].type != GRADIENT_TYPE_VELX && grad_elements[k].type != GRADIENT_TYPE_VELY)
                    grad_elements[k].value1 = grad_elements[k].value0;
                }

            }

          velGas0[0] += vx;
          velGas0[1] += vy;
          velGas1[0] += vx;
          velGas1[1] += vy;

          GVelx->value0 = velGas0[0];
          GVely->value0 = velGas0[1];

          GVelx->value1 = velGas1[0];
          GVely->value1 = velGas1[1];
        }
#endif

      /* if the cell q0 is a local particle, construct the gradient estimate, and the minmax values */
      if(p0->task == ThisTask && q0 >= 0 && q0 < NumGas)
        {
          if(TimeBinSynchronized[P[q0].TimeBinHydro])
            {
              fac0 = 0.5 * VF[i].area / SphP[q0].Volume;
              facA = VF[i].area / (nn * SphP[q0].Volume);

              for(k = 0; k < N_Grad; k++)
                {
                  double value0;

                  if(grad_elements[k].type == GRADIENT_TYPE_VELX)
                    value0 = grad_elements[k].value0 + nn * n[0] * atime * hubble_a;
                  else if(grad_elements[k].type == GRADIENT_TYPE_VELY)
                    value0 = grad_elements[k].value0 + nn * n[1] * atime * hubble_a;
                  else if(grad_elements[k].type == GRADIENT_TYPE_VELZ)
                    value0 = grad_elements[k].value0 + nn * n[2] * atime * hubble_a;
                  else
                    value0 = grad_elements[k].value0;

                  data = (MySingle *) (((char *) (&(SphP[q0].Grad))) + grad_elements[k].offset_grad);

                  for(j = 0; j < 3; j++)
                    {
                      data[j] += fac0 * n[j] * value0;
                      data[j] += facA * c[j] * (value0 - grad_elements[k].value1);
                    }

                  /* force gradient to be flat if we are at a outflow boundary */
#ifdef REFLECTIVE_X
                  if(((*image_flags0 & REFL_X_FLAGS) && (*image_flags0 & OUTFLOW_X)) || ((*image_flags1 & REFL_X_FLAGS) && (*image_flags1 & OUTFLOW_X)))
                    data[0] = 0;
#endif
#ifdef REFLECTIVE_Y
                  if(((*image_flags0 & REFL_Y_FLAGS) && (*image_flags0 & OUTFLOW_Y)) || ((*image_flags1 & REFL_Y_FLAGS) && (*image_flags1 & OUTFLOW_Y)))
                    data[1] = 0;
#endif
#ifdef REFLECTIVE_Z
                  if(((*image_flags0 & REFL_Z_FLAGS) && (*image_flags0 & OUTFLOW_Z)) || ((*image_flags1 & REFL_Z_FLAGS) && (*image_flags1 & OUTFLOW_Z)))
                    data[2] = 0;
#endif

                  if(VF[i].area > 1.0e-5 * SphP[q0].SurfaceArea)
                    {
                      if(grad_elements[k].max_value[q0] < value0)
                        grad_elements[k].max_value[q0] = value0;

                      if(grad_elements[k].min_value[q0] > value0)
                        grad_elements[k].min_value[q0] = value0;
                    }
                }

              /* now also construct an estimate of the curl */
              curl[3 * q0 + 0] += fac0 * (n[1] * GVelz->value0 - n[2] * GVely->value0);
              curl[3 * q0 + 1] += fac0 * (n[2] * GVelx->value0 - n[0] * GVelz->value0);
              curl[3 * q0 + 2] += fac0 * (n[0] * GVely->value0 - n[1] * GVelx->value0);

              /* now add the contribution to the velocity divergence */
#if defined (OUTPUT_DIVVEL) || defined(TGCHEM) || defined(SGCHEM)
              SphP[q0].DivVel += fac0 * (n[0] * GVelx->value0 + n[1] * GVely->value0 + n[2] * GVelz->value0);
              SphP[q0].DivVel += facA * (c[0] * (GVelx->value0 - GVelx->value1) + c[1] * (GVely->value0 - GVely->value1) + c[2] * (GVelz->value0 - GVelz->value1));
#endif

            }
        }

      /* if the cell q1 is a local particle, construct the gradient estimate, and the minmax values */
      if(p1->task == ThisTask && q1 >= 0 && q1 < NumGas)
        {
          if(TimeBinSynchronized[P[q1].TimeBinHydro])
            {
              fac1 = -0.5 * VF[i].area / SphP[q1].Volume;
              facB = VF[i].area / (nn * SphP[q1].Volume);

              for(k = 0; k < N_Grad; k++)
                {
                  double value1;

                  if(grad_elements[k].type == GRADIENT_TYPE_VELX)
                    value1 = grad_elements[k].value1 - nn * n[0] * atime * hubble_a;
                  else if(grad_elements[k].type == GRADIENT_TYPE_VELY)
                    value1 = grad_elements[k].value1 - nn * n[1] * atime * hubble_a;
                  else if(grad_elements[k].type == GRADIENT_TYPE_VELZ)
                    value1 = grad_elements[k].value1 - nn * n[2] * atime * hubble_a;
                  else
                    value1 = grad_elements[k].value1;

                  data = (MySingle *) (((char *) (&(SphP[q1].Grad))) + grad_elements[k].offset_grad);
                  for(j = 0; j < 3; j++)
                    {
                      data[j] += fac1 * n[j] * value1;
                      data[j] += facB * c[j] * (value1 - grad_elements[k].value0);
                    }

                  /* force the gradient to be flat if we are at outflow boundary */
#ifdef REFLECTIVE_X
                  if(((*image_flags0 & REFL_X_FLAGS) && (*image_flags0 & OUTFLOW_X)) || ((*image_flags1 & REFL_X_FLAGS) && (*image_flags1 & OUTFLOW_X)))
                    data[0] = 0;
#endif
#ifdef REFLECTIVE_Y
                  if(((*image_flags0 & REFL_Y_FLAGS) && (*image_flags0 & OUTFLOW_Y)) || ((*image_flags1 & REFL_Y_FLAGS) && (*image_flags1 & OUTFLOW_Y)))
                    data[1] = 0;
#endif
#ifdef REFLECTIVE_Z
                  if(((*image_flags0 & REFL_Z_FLAGS) && (*image_flags0 & OUTFLOW_Z)) || ((*image_flags1 & REFL_Z_FLAGS) && (*image_flags1 & OUTFLOW_Z)))
                    data[2] = 0;
#endif

                  if(VF[i].area > 1.0e-5 * SphP[q1].SurfaceArea)
                    {
                      if(grad_elements[k].max_value[q1] < value1)
                        grad_elements[k].max_value[q1] = value1;

                      if(grad_elements[k].min_value[q1] > value1)
                        grad_elements[k].min_value[q1] = value1;
                    }
                }
              /* now also construct an estimate of the curl */
              curl[3 * q1 + 0] += fac1 * (n[1] * GVelz->value1 - n[2] * GVely->value1);
              curl[3 * q1 + 1] += fac1 * (n[2] * GVelx->value1 - n[0] * GVelz->value1);
              curl[3 * q1 + 2] += fac1 * (n[0] * GVely->value1 - n[1] * GVelx->value1);

              /* now add the contribution to the velocity divergence */
#if defined (OUTPUT_DIVVEL) || defined(TGCHEM) || defined(SGCHEM)
              SphP[q1].DivVel += fac1 * (n[0] * GVelx->value1 + n[1] * GVely->value1 + n[2] * GVelz->value1);
              SphP[q1].DivVel += facB * (c[0] * (GVelx->value1 - GVelx->value0) + c[1] * (GVely->value1 - GVely->value0) + c[2] * (GVelz->value1 - GVelz->value0));
#endif
            }
        }
    }


#ifdef DEREFINE_GENTLY
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

#ifdef SINKS
      /* cells that are inside the sink accretion radius are eligible for derefinement */
      if(SphP[i].InAccrRadius)
        {
          SphP[i].InAccrRadius = 0;
          continue;
        }
#endif

      if(SphP[i].Density > GENTLE_DEREFINE_FACTOR * GDensity->min_value[i])
        SphP[i].DoNotDerefFlag = 1;

      if(SphP[i].Density < GDensity->max_value[i] / GENTLE_DEREFINE_FACTOR)
        SphP[i].DoNotDerefFlag = 1;

      if(SphP[i].Utherm > GENTLE_DEREFINE_FACTOR * GUtherm->min_value[i])
        SphP[i].DoNotDerefFlag = 1;

      if(SphP[i].Utherm < GUtherm->max_value[i] / GENTLE_DEREFINE_FACTOR)
        SphP[i].DoNotDerefFlag = 1;

      if(P[i].Vel[0] / atime > GENTLE_DEREFINE_FACTOR * GVelx->min_value[i])
        SphP[i].DoNotDerefFlag = 1;

      if(P[i].Vel[0] / atime < GVelx->max_value[i] / GENTLE_DEREFINE_FACTOR)
        SphP[i].DoNotDerefFlag = 1;

      if(P[i].Vel[1] / atime > GENTLE_DEREFINE_FACTOR * GVely->min_value[i])
        SphP[i].DoNotDerefFlag = 1;

      if(P[i].Vel[1] / atime < GVely->max_value[i] / GENTLE_DEREFINE_FACTOR)
        SphP[i].DoNotDerefFlag = 1;

      if(P[i].Vel[2] / atime > GENTLE_DEREFINE_FACTOR * GVelz->min_value[i])
        SphP[i].DoNotDerefFlag = 1;

      if(P[i].Vel[2] / atime < GVelz->max_value[i] / GENTLE_DEREFINE_FACTOR)
        SphP[i].DoNotDerefFlag = 1;

      if(SphP[i].Pressure > GENTLE_DEREFINE_FACTOR * GPressure->min_value[i])
        SphP[i].DoNotDerefFlag = 1;

      if(SphP[i].Pressure < GPressure->max_value[i] / GENTLE_DEREFINE_FACTOR)
        SphP[i].DoNotDerefFlag = 1;
    }
#endif

#ifdef TVD_SLOPE_LIMITER
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;
      SphP[i].GradUl = SphP[i].Grad;
    }
#endif

#ifdef NUCLEAR_NETWORK
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      for(j = 0; j < 3; j++)
        {
          SphP[i].Grad.drhoU[j] = SphP[i].Grad.drho[j];
          SphP[i].Grad.dpressU[j] = SphP[i].Grad.dpress[j];
        }
    }
#endif

#ifdef MHD
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;
      SphP[i].CurlB[0] = SphP[i].Grad.dB[2][1] - SphP[i].Grad.dB[1][2];
      SphP[i].CurlB[1] = SphP[i].Grad.dB[0][2] - SphP[i].Grad.dB[2][0];
      SphP[i].CurlB[2] = SphP[i].Grad.dB[1][0] - SphP[i].Grad.dB[0][1];
    }
#endif

#ifndef UNLIMITED_GRADIENTS
#if defined(SHOCK_FINDER_BEFORE_OUTPUT) || defined(SHOCK_FINDER_ON_THE_FLY)
  if(SfVars.limit_gradients)
    {
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
                          for(k = 0; k < N_Grad; k++)
                            {
                              if(grad_elements[k].type == GRADIENT_TYPE_FLD)
                                continue;

                              if((grad_elements[k].type == GRADIENT_TYPE_VELX) || (grad_elements[k].type == GRADIENT_TYPE_VELY) || (grad_elements[k].type == GRADIENT_TYPE_VELZ))
                                {
                                  value = *(MyFloat *) (((char *) (&P[q])) + grad_elements[k].offset);
                                  value /= atime;
                                }
                              else
                                value = *(MyFloat *) (((char *) (&SphP[q])) + grad_elements[k].offset);

                              data = (MySingle *) (((char *) (&(SphP[q].Grad))) + grad_elements[k].offset_grad);

#ifdef SPECIAL_BOUNDARY
                              //if(p->ID != -1 && p->ID != -2)
#endif
                              //if((grad_elements[k].type != GRADIENT_TYPE_AX) && (grad_elements[k].type != GRADIENT_TYPE_AY) && (grad_elements[k].type != GRADIENT_TYPE_AZ) && (grad_elements[k].type != GRADIENT_TYPE_ASHIFT))                         
                              limit_gradient(d, value, grad_elements[k].min_value[q], grad_elements[k].max_value[q], data);

#ifdef WINDTUNNEL

                              MyIDType id;
                              int l = 0;
                              id = p->ID;
                              MyFloat *gradient;
                              gradient = data;
                              /*set the gradients to zero within the buffer cells */
                              if((id >= BOUNDARY_INFLOWOUTFLOW_MINID) && (id < BOUNDARY_INFLOWOUTFLOW_MAXID))
                                {
                                  for(l = 0; l < 3; l++)
                                    data[l] = 0.0;
                                }
#endif
#ifdef SPECIAL_BOUNDARY
                              if(All.SpecialBoundaryType == 3 || All.SpecialBoundaryType == 4)
                                {
                                  MyIDType id;
                                  int l = 0;
                                  id = p->ID;
                                  MyFloat *gradient;
                                  gradient = data;

                                  /*set the normal components of the gradients to zero on both sides
                                     of the transparent boundary */
                                  if(id == -1 || id == -2)
                                    {

                                      gradient = data;
                                      for(l = 0; l < 3; l++)
                                        data[l] = gradient[l] - (gradient[0] * n[0] + gradient[1] * n[1] + gradient[2] * n[2]) * n[l];
                                    }
                                }
#endif

                            }

#ifndef DISABLE_VELOCITY_CSND_SLOPE_LIMITING
                          /* let's now limit the overall size of the velocity gradient */

                          MySingle *grad_vx = (MySingle *) (((char *) (&(SphP[q].Grad))) + GVelx->offset_grad);
                          MySingle *grad_vy = (MySingle *) (((char *) (&(SphP[q].Grad))) + GVely->offset_grad);
                          MySingle *grad_vz = (MySingle *) (((char *) (&(SphP[q].Grad))) + GVelz->offset_grad);

                          limit_vel_gradient(d, grad_vx, grad_vy, grad_vz, get_sound_speed(q));
#endif
                        }
                    }
                }
            }
        }
#if defined(SHOCK_FINDER_BEFORE_OUTPUT) || defined(SHOCK_FINDER_ON_THE_FLY)
    }
#endif

#ifdef REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED
  /* compute magnitude of curl */
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;
      SphP[i].CurlVel = sqrt(curl[3 * i + 0] * curl[3 * i + 0] + curl[3 * i + 1] * curl[3 * i + 1] + curl[3 * i + 2] * curl[3 * i + 2]);
    }
#endif

#if defined(VS_TURB) || defined(AB_TURB) || defined(CALCULATE_QUANTITIES_IN_POSTPROCESS)
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      mpi_printf("CQIP: Entered Here\n") ;
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      mpi_printf("CQIP: Entered Here\n") ;

      SphP[i].Vorticity[0] = curl[3 * i + 0];
      SphP[i].Vorticity[1] = curl[3 * i + 1];
      SphP[i].Vorticity[2] = curl[3 * i + 2];
    }
#endif

#endif


#ifdef NO_SCALAR_GRADIENTS
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      for(k = 0; k < N_Scalar; k++)
        {
          SphP[i].Grad.dscalars[k][0] = 0;
          SphP[i].Grad.dscalars[k][1] = 0;
          SphP[i].Grad.dscalars[k][2] = 0;
        }
    }
#endif



  myfree(curl);

  for(k = N_Grad - 1; k >= 0; k--)
    {
      myfree(grad_elements[k].max_value);
      myfree(grad_elements[k].min_value);
    }

  CPU_Step[CPU_GRADIENTS] += measure_time();
}

void limit_vel_gradient(double *d, MySingle * grad_vx, MySingle * grad_vy, MySingle * grad_vz, double csnd)
{
  int i;

#define VEL_GRADIENT_LIMIT_FAC 1.0

  if(All.ComovingIntegrationOn)
    {
      grad_vx[0] -= atime * hubble_a;
      grad_vy[1] -= atime * hubble_a;
      grad_vz[2] -= atime * hubble_a;
    }

  double dvx = fabs(grad_vx[0] * d[0] + grad_vx[1] * d[1] + grad_vx[2] * d[2]);
  double dvy = fabs(grad_vy[0] * d[0] + grad_vy[1] * d[1] + grad_vy[2] * d[2]);
  double dvz = fabs(grad_vz[0] * d[0] + grad_vz[1] * d[1] + grad_vz[2] * d[2]);

  if(dvx > VEL_GRADIENT_LIMIT_FAC * csnd)
    {
      double fac = VEL_GRADIENT_LIMIT_FAC * csnd / dvx;
      for(i = 0; i < 3; i++)
        {
          grad_vx[i] *= fac;
        }
    }

  if(dvy > VEL_GRADIENT_LIMIT_FAC * csnd)
    {
      double fac = VEL_GRADIENT_LIMIT_FAC * csnd / dvy;
      for(i = 0; i < 3; i++)
        {
          grad_vy[i] *= fac;
        }
    }
  if(dvz > VEL_GRADIENT_LIMIT_FAC * csnd)
    {
      double fac = VEL_GRADIENT_LIMIT_FAC * csnd / dvz;
      for(i = 0; i < 3; i++)
        {
          grad_vz[i] *= fac;
        }
    }

  if(All.ComovingIntegrationOn)
    {
      grad_vx[0] += atime * hubble_a;
      grad_vy[1] += atime * hubble_a;
      grad_vz[2] += atime * hubble_a;
    }
}

void limit_gradient(double *d, double phi, double min_phi, double max_phi, MySingle * dphi)
{
  double dp, fac;

  dp = dphi[0] * d[0] + dphi[1] * d[1] + dphi[2] * d[2];

  if(dp > 0)
    {
      if(phi + dp > max_phi)
        {
          if(max_phi > phi)
            fac = (max_phi - phi) / dp;
          else
            fac = 0;

          if(fac < 0 || fac > 1)
            terminate("fac=%g\ndp=%g max_phi=%g phi=%g", fac, dp, max_phi, phi);

          dphi[0] *= fac;
          dphi[1] *= fac;
          dphi[2] *= fac;
        }
    }
  else if(dp < 0)
    {
      if(phi + dp < min_phi)
        {
          if(min_phi < phi)
            fac = (min_phi - phi) / dp;
          else
            fac = 0;

          if(fac < 0 || fac > 1)
            terminate("fac=%g\ndp=%g max_phi=%g phi=%g", fac, dp, max_phi, phi);

          dphi[0] *= fac;
          dphi[1] *= fac;
          dphi[2] *= fac;
        }
    }
}

#ifdef OUTPUT_VERTEX_VELOCITY_DIVERGENCE
void calculate_vertex_velocity_divergence(void)
{
  voronoi_update_ghost_velvertex();

  int idx, i, j, t0, t1, q0, q1, s0, s1;

  double fac0, fac1, n[3], nn, velGas0[3] = { 0, 0, 0 }, velGas1[3] =
  {
  0, 0, 0};

  double facA, facB, c[3];

#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
  unsigned int *image_flags0, *image_flags1;
#endif

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      SphP[i].DivVelVertex = 0;
    }

  for(i = 0; i < Nvf; i++)
    {
      point *p1 = &DP[VF[i].p1];
      point *p0 = &DP[VF[i].p2];

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

      /* let's get the physical quantities of interest for the particle one one of the  sides */
      if(t0 == ThisTask)
        {
          if(s0 >= NumGas)
            s0 -= NumGas;

          if(s0 >= 0)
            {
              velGas0[0] = SphP[s0].VelVertex[0];
              velGas0[1] = SphP[s0].VelVertex[1];
              velGas0[2] = SphP[s0].VelVertex[2];
            }
        }
      else
        {
          velGas0[0] = PrimExch[s0].VelVertex[0];
          velGas0[1] = PrimExch[s0].VelVertex[1];
          velGas0[2] = PrimExch[s0].VelVertex[2];
        }

      /* let's get the physical quantities of interest for the particle on the other side */
      if(t1 == ThisTask)
        {
          if(s1 >= NumGas)
            s1 -= NumGas;

          if(s1 >= 0)
            {
              velGas1[0] = SphP[s1].VelVertex[0];
              velGas1[1] = SphP[s1].VelVertex[1];
              velGas1[2] = SphP[s1].VelVertex[2];
            }
        }
      else
        {
          velGas1[0] = PrimExch[s1].VelVertex[0];
          velGas1[1] = PrimExch[s1].VelVertex[1];
          velGas1[2] = PrimExch[s1].VelVertex[2];
        }

#if defined(REFLECTIVE_X)
      if((*image_flags0 & REFL_X_FLAGS) && !(*image_flags0 & OUTFLOW_X))
        velGas0[0] *= -1;
      if((*image_flags0 & REFL_X_FLAGS) && (*image_flags0) && (REFLECTIVE_X == 3))
        {
          velGas0[1] *= -1;
          velGas0[2] *= -1;
        }
#endif
#if defined(REFLECTIVE_Y)
      if((*image_flags0 & REFL_Y_FLAGS) && !(*image_flags0 & OUTFLOW_Y))
        velGas0[1] *= -1;
      if((*image_flags0 & REFL_Y_FLAGS) && (*image_flags0) && (REFLECTIVE_Y == 3))
        {
          velGas0[0] *= -1;
          velGas0[2] *= -1;
        }
#endif
#if defined(REFLECTIVE_Z)
      if((*image_flags0 & REFL_Z_FLAGS) && !(*image_flags0 & OUTFLOW_Z))
        velGas0[2] *= -1;
      if((*image_flags0 & REFL_Z_FLAGS) && (*image_flags0) && (REFLECTIVE_Z == 3))
        {
          velGas0[0] *= -1;
          velGas0[1] *= -1;
        }
#endif

#if defined(REFLECTIVE_X)
      if((*image_flags1 & REFL_X_FLAGS) && !(*image_flags1 & OUTFLOW_X))
        velGas1[0] *= -1;
      if((*image_flags1 & REFL_X_FLAGS) && (*image_flags0) && (REFLECTIVE_X == 3))
        {
          velGas1[1] *= -1;
          velGas1[2] *= -1;
        }
#endif
#if defined(REFLECTIVE_Y)
      if((*image_flags1 & REFL_Y_FLAGS) && !(*image_flags1 & OUTFLOW_Y))
        velGas1[1] *= -1;
      if((*image_flags1 & REFL_Y_FLAGS) && (*image_flags0) && (REFLECTIVE_Y == 3))
        {
          velGas1[0] *= -1;
          velGas1[2] *= -1;
        }
#endif
#if defined(REFLECTIVE_Z)
      if((*image_flags1 & REFL_Z_FLAGS) && !(*image_flags1 & OUTFLOW_Z))
        velGas1[2] *= -1;
      if((*image_flags1 & REFL_Z_FLAGS) && (*image_flags0) && (REFLECTIVE_Z == 3))
        {
          velGas1[0] *= -1;
          velGas1[1] *= -1;
        }
#endif

      /* if the cell q0 is a local particle, construct the gradient estimate, and the minmax values */
      if(p0->task == ThisTask && q0 >= 0 && q0 < NumGas)
        {
          if(TimeBinSynchronized[P[q0].TimeBinHydro])
            {
              fac0 = 0.5 * VF[i].area / SphP[q0].Volume;
              facA = VF[i].area / (nn * SphP[q0].Volume);
              SphP[q0].DivVelVertex += fac0 * (n[0] * velGas0[0] + n[1] * velGas0[1] + n[2] * velGas0[2]);
              SphP[q0].DivVelVertex += facA * (c[0] * (velGas0[0] - velGas1[0]) + c[1] * (velGas0[1] - velGas1[1]) + c[2] * (velGas0[2] - velGas1[2]));
            }
        }

      /* if the cell q1 is a local particle, construct the gradient estimate, and the minmax values */
      if(p1->task == ThisTask && q1 >= 0 && q1 < NumGas)
        {
          if(TimeBinSynchronized[P[q1].TimeBinHydro])
            {
              fac1 = -0.5 * VF[i].area / SphP[q1].Volume;
              facB = VF[i].area / (nn * SphP[q1].Volume);

              SphP[q1].DivVelVertex += fac1 * (n[0] * velGas1[0] + n[1] * velGas1[1] + n[2] * velGas1[2]);
              SphP[q1].DivVelVertex += facB * (c[0] * (velGas1[0] - velGas0[0]) + c[1] * (velGas1[1] - velGas0[1]) + c[2] * (velGas1[2] - velGas0[2]));
            }
        }
    }
}
#endif

#endif

#endif /* VORONOI */
