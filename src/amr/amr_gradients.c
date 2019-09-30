/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/amr/amr_gradients.c
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

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "../allvars.h"
#include "../proto.h"

#if defined(AMR) && defined(AMR_GRADIENTS)

#ifdef TVD_SLOPE_LIMITER2
#if !defined(TVD_SLOPE_LIMITER_VANLEER) && !defined(TVD_SLOPE_LIMITER_SUPERBEE) && !defined(TVD_SLOPE_LIMITER_ALBADA) && !defined(TVD_SLOPE_LIMITER_MINBEE)
#error "Choose a slope limiter"
#endif
#endif

MyFloat amr_gradient_get_value(int node, int direction, int gradient);

static double atime, hubble_a;

void calculate_gradients(void)
{
  TIMER_START(CPU_GRADIENTS) if(All.TotNumGas == 0)
    return;

  int level;

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

#if !defined(UNLIMITED_GRADIENTS) && !defined(TVD_SLOPE_LIMITER2)
  int k;
  for(k = 0; k < N_Grad; k++)
    {
      grad_elements[k].min_value = mymalloc("gradmin", NumGas * sizeof(double));
      grad_elements[k].max_value = mymalloc("gradmax", NumGas * sizeof(double));
    }
#endif

  for(level = Mesh.minlevel; level <= Mesh.maxlevel; level++)
    {
      int node;

      node = Mesh.lastinlevel[level];
      while(node >= 0)
        {
          if(node < Ngb_MaxPart)
            {
              int dir;

#if !defined(UNLIMITED_GRADIENTS) && !defined(TVD_SLOPE_LIMITER2)
              int k;
              for(k = 0; k < N_Grad; k++)
                {
                  grad_elements[k].min_value[node] = +MAX_REAL_NUMBER;
                  grad_elements[k].max_value[node] = -MAX_REAL_NUMBER;

                  MyFloat value = amr_gradient_get_value(node, -1, k);

                  if(All.ComovingIntegrationOn && (grad_elements[k].type == GRADIENT_TYPE_VELX || grad_elements[k].type == GRADIENT_TYPE_VELY || grad_elements[k].type == GRADIENT_TYPE_VELZ))
                    {
                      value /= atime;
                      value /= atime;
                    }

                  if(grad_elements[k].max_value[node] < value)
                    {
                      grad_elements[k].max_value[node] = value;
                    }

                  if(grad_elements[k].min_value[node] > value)
                    {
                      grad_elements[k].min_value[node] = value;
                    }

                }
#endif

              for(dir = 0; dir < 3; dir++)
                {
                  int k;
                  for(k = 0; k < N_Grad; k++)
                    {
                      double *data = (double *) (((char *) (&(SphP[node].Grad))) + grad_elements[k].offset_grad);

#ifdef TVD_SLOPE_LIMITER
                      double *dataUl = (double *) (((char *) (&(SphP[node].GradUl))) + grad_elements[k].offset_grad);
#endif

                      if(dir < NUMDIMS)
                        {
                          MyFloat lvalue = amr_gradient_get_value(node, 2 * dir, k);
                          MyFloat rvalue = amr_gradient_get_value(node, 2 * dir + 1, k);

                          if(All.ComovingIntegrationOn && (grad_elements[k].type == GRADIENT_TYPE_VELX || grad_elements[k].type == GRADIENT_TYPE_VELY || grad_elements[k].type == GRADIENT_TYPE_VELZ))
                            {
                              lvalue /= atime;
                              rvalue /= atime;
                            }

                          data[dir] = (rvalue - lvalue) / (2 * amr_length[Mesh.DP[node].level]);

#if !defined(UNLIMITED_GRADIENTS) && !defined(TVD_SLOPE_LIMITER2)
                          if(grad_elements[k].max_value[node] < lvalue)
                            {
                              grad_elements[k].max_value[node] = lvalue;
                            }

                          if(grad_elements[k].min_value[node] > lvalue)
                            {
                              grad_elements[k].min_value[node] = lvalue;
                            }

                          if(grad_elements[k].max_value[node] < rvalue)
                            {
                              grad_elements[k].max_value[node] = rvalue;
                            }

                          if(grad_elements[k].min_value[node] > rvalue)
                            {
                              grad_elements[k].min_value[node] = rvalue;
                            }
#endif

#ifdef TVD_SLOPE_LIMITER2

                          MyFloat value = amr_gradient_get_value(node, -1, k);

                          if(All.ComovingIntegrationOn && (grad_elements[k].type == GRADIENT_TYPE_VELX || grad_elements[k].type == GRADIENT_TYPE_VELY || grad_elements[k].type == GRADIENT_TYPE_VELZ))
                            {
                              value /= atime;
                              value /= atime;
                            }

                          double rf = (value - lvalue) / (rvalue - value);

                          double psi;

                          if(rf != rf)
                            {
                              psi = 0.;
                            }
                          else
                            {
#ifdef TVD_SLOPE_LIMITER_VANLEER
                              if(rf <= 0.)
                                {
                                  psi = 0.;
                                }
                              else
                                {
                                  double psi_R = 2. / (1. + rf);
                                  psi = dmin(2. * rf / (1. + rf), psi_R);
                                }
#endif

#ifdef TVD_SLOPE_LIMITER_SUPERBEE
                              if(rf <= 0.)
                                {
                                  psi = 0.;
                                }
                              else if(rf <= 0.5)
                                {
                                  psi = 2. * rf;
                                }
                              else if(rf <= 1.)
                                {
                                  psi = 1.;
                                }
                              else
                                {
                                  double psi_R = 2. / (1. + rf);
                                  psi = dmin(dmin(rf, psi_R), 2.);
                                }
#endif

#ifdef TVD_SLOPE_LIMITER_ALBADA
                              if(rf <= 0.)
                                {
                                  psi = 0.;
                                }
                              else
                                {
                                  double psi_R = 2. / (1. + rf);
                                  psi = dmin((rf * (1. + rf)) / (1. + rf * rf), psi_R);
                                }
#endif

#ifdef TVD_SLOPE_LIMITER_MINBEE
                              if(rf <= 0.)
                                {
                                  psi = 0.;
                                }
                              else if(rf <= 1.)
                                {
                                  psi = rf;
                                }
                              else
                                {
                                  double psi_R = 2. / (1. + rf);
                                  psi = dmin(1., psi_R);
                                }
#endif
                            }
                          /* limiter is applied */ ;
                          data[dir] *= psi;
#endif
                        }
                      else
                        {
                          data[dir] = 0.;
                        }
#ifdef TVD_SLOPE_LIMITER
                      dataUl[dir] = data[dir];
#endif
                    }
                }

#ifdef TVD_SLOPE_LIMITER2
#ifndef DISABLE_VELOCITY_CSND_SLOPE_LIMITING
              for(dir = 0; dir < NUMDIMS; dir++)
                {
                  int side;
                  for(side = -1; side <= 1; side += 2)
                    {
                      double d[3];
                      d[0] = 0.;
                      d[1] = 0.;
                      d[2] = 0.;

                      d[dir] += side * amr_length[Mesh.DP[node].level + 1];

                      /* let's now limit the overall size of the velocity gradient */
                      double *grad_vx = (double *) (((char *) (&(SphP[node].Grad))) + GVelx->offset_grad);
                      double *grad_vy = (double *) (((char *) (&(SphP[node].Grad))) + GVely->offset_grad);
                      double *grad_vz = (double *) (((char *) (&(SphP[node].Grad))) + GVelz->offset_grad);

                      limit_vel_gradient(d, grad_vx, grad_vy, grad_vz, get_sound_speed(node));
                    }
                }
#endif
#endif
            }

          if(node >= Ngb_MaxPart)
            node = Ngb_Nodes[node].nextinlevel;
          else
            node = Mesh.DP[node].nextinlevel;
        }

      exchange_gradients();
    }

#if !defined(UNLIMITED_GRADIENTS) && !defined(TVD_SLOPE_LIMITER2)
  for(level = Mesh.minlevel; level <= Mesh.maxlevel; level++)
    {
      int node;

      node = Mesh.lastinlevel[level];
      while(node >= 0)
        {
          if(node < Ngb_MaxPart)
            {


              int dir;
              for(dir = 0; dir < NUMDIMS; dir++)
                {

                  int side;
                  for(side = -1; side <= 1; side += 2)
                    {
                      double d[3];
                      d[0] = 0.;
                      d[1] = 0.;
                      d[2] = 0.;

                      d[dir] += side * amr_length[Mesh.DP[node].level + 1];

                      for(k = 0; k < N_Grad; k++)
                        {
                          double value;
                          double *data;
                          if((grad_elements[k].type == GRADIENT_TYPE_VELX) || (grad_elements[k].type == GRADIENT_TYPE_VELY) || (grad_elements[k].type == GRADIENT_TYPE_VELZ))
                            {
                              value = *(MyFloat *) (((char *) (&P[node])) + grad_elements[k].offset);
                              value /= atime;
                            }
                          else
                            {
                              value = *(MyFloat *) (((char *) (&SphP[node])) + grad_elements[k].offset);
                            }

                          data = (double *) (((char *) (&(SphP[node].Grad))) + grad_elements[k].offset_grad);


                          limit_gradient(d, value, grad_elements[k].min_value[node], grad_elements[k].max_value[node], data);

                        }

#ifndef DISABLE_VELOCITY_CSND_SLOPE_LIMITING
                      /* let's now limit the overall size of the velocity gradient */

                      double *grad_vx = (double *) (((char *) (&(SphP[node].Grad))) + GVelx->offset_grad);
                      double *grad_vy = (double *) (((char *) (&(SphP[node].Grad))) + GVely->offset_grad);
                      double *grad_vz = (double *) (((char *) (&(SphP[node].Grad))) + GVelz->offset_grad);

                      limit_vel_gradient(d, grad_vx, grad_vy, grad_vz, get_sound_speed(node));
#endif
                    }
                }
            }

          if(node >= Ngb_MaxPart)
            node = Ngb_Nodes[node].nextinlevel;
          else
            node = Mesh.DP[node].nextinlevel;
        }

      //exchange_gradients();
    }

  for(k = N_Grad - 1; k >= 0; k--)
    {
      myfree(grad_elements[k].max_value);
      myfree(grad_elements[k].min_value);
    }
#endif
TIMER_STOP(CPU_GRADIENTS)}


MyFloat amr_gradient_get_value(int node, int direction, int k)
{
  MyFloat value;
  MyFloat *grad = NULL;
  int node_other;
  if(direction >= 0)
    {
      node_other = Mesh.DP[node].neighbors[direction];
    }
  else
    {
      node_other = node;
    }

  if(node_other < Ngb_MaxPart)
    {
      if((grad_elements[k].type == GRADIENT_TYPE_VELX) || (grad_elements[k].type == GRADIENT_TYPE_VELY) || (grad_elements[k].type == GRADIENT_TYPE_VELZ))
        value = *(MyFloat *) (((char *) (&P[node_other])) + grad_elements[k].offset);
      else
        value = *(MyFloat *) (((char *) (&SphP[node_other])) + grad_elements[k].offset);

      grad = (MyFloat *) (((char *) (&SphP[node_other].Grad)) + grad_elements[k].offset_grad);
    }
  else if(node_other < Ngb_MaxPart + Mesh.nodes_total)
    {
      if(grad_elements[k].type == GRADIENT_TYPE_DENSITY)
        {
          value = Ngb_Nodes[node_other].hydro.mass / amr_volume[Ngb_Nodes[node_other].level];
        }
      else if(grad_elements[k].type == GRADIENT_TYPE_VELX)
        {
          value = Ngb_Nodes[node_other].hydro.momentum[0] / Ngb_Nodes[node_other].hydro.mass;
        }
      else if(grad_elements[k].type == GRADIENT_TYPE_VELY)
        {
          value = Ngb_Nodes[node_other].hydro.momentum[1] / Ngb_Nodes[node_other].hydro.mass;
        }
      else if(grad_elements[k].type == GRADIENT_TYPE_VELZ)
        {
          value = Ngb_Nodes[node_other].hydro.momentum[2] / Ngb_Nodes[node_other].hydro.mass;
        }
      else if(grad_elements[k].type == GRADIENT_TYPE_PRESSURE)
        {
          value = GAMMA_MINUS1 * Ngb_Nodes[node_other].hydro.etherm / amr_volume[Ngb_Nodes[node_other].level];
        }
      else
        {
          terminate("not supported");
        }

      return value;
    }
  else
    {
      node_other -= Mesh.nodes_total;
      int i = Mesh.DP[node_other].index;
      value = *(MyFloat *) (((char *) (&PrimExch[i])) + grad_elements[k].offset_exch);
      grad = (MyFloat *) (((char *) (&SphP[i].Grad)) + grad_elements[k].offset_grad);
    }

  if(Mesh.DP[node_other].level < Mesh.DP[node].level)
    {
      assert(grad);

      MyFloat dx, dy, dz;
      MyFloat pos[3];

      pos[0] = Mesh.DP[node].x;
      pos[1] = Mesh.DP[node].y;
      pos[2] = Mesh.DP[node].z;

      int shift_dir = direction / 2;
      int shift = 2*(direction % 2) - 1;
      pos[shift_dir] += shift * amr_length[Mesh.DP[node].level];

      dx = pos[0] - Mesh.DP[node_other].x;
      dy = pos[1] - Mesh.DP[node_other].y;
      dz = pos[2] - Mesh.DP[node_other].z;

      if(dx > boxSize_X / 2)
        dx -= boxSize_X;
      else if(dx < -boxSize_X / 2)
        dx += boxSize_X;

      if(dy > boxSize_Y / 2)
        dy -= boxSize_Y;
      else if(dy < -boxSize_Y / 2)
        dy += boxSize_Y;

      if(dz > boxSize_Z / 2)
        dz -= boxSize_Z;
      else if(dz < -boxSize_Z / 2)
        dz += boxSize_Z;

      value += grad[0] * dx + grad[1] * dy + grad[2] * dz;
    }

  return value;

}

void limit_vel_gradient(double *d, double *grad_vx, double *grad_vy, double *grad_vz, double csnd)
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

void limit_gradient(double *d, double phi, double min_phi, double max_phi, double *dphi)
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
#endif
