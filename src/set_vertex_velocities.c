/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/set_vertex_velocities.c
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

#ifdef ONEDIMS_SPHERICAL
static void validate_vertex_velocities_1d();
#endif

void set_vertex_velocities(void)
{
  TIMER_START(CPU_SET_VERTEXVELS);

  int idx, i, j;
  double dt;

#if defined (VORONOI_STATIC_MESH) || defined (AMR) || defined (NOHYDRO)
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      for(j = 0; j < 3; j++)
        SphP[i].VelVertex[j] = 0;
    }
  TIMER_STOP(CPU_SET_VERTEXVELS);
  return;
#endif

#ifdef DMFIXED
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;
      if(P[i].Type == 1)
        for(j = 0; j < 3; j++)
          P[i].Vel[j] = 0;
    }
#endif

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

#ifdef MESHRELAX
      for(j = 0; j < 3; j++)
        SphP[i].VelVertex[j] = 0;
#else
      for(j = 0; j < 3; j++)
        SphP[i].VelVertex[j] = P[i].Vel[j];     /* make cell velocity equal to fluid's velocity */
#endif

#ifdef GENERAL_RELATIVITY

      // note that the primitive velocity is not the coordinate velocity but the velocity measured by an Eulerian observer
      // hence we need to get the coordinate velocity

      double xmet = SphP[i].Center[0];
      double ymet = SphP[i].Center[1];
      double zmet = SphP[i].Center[2];

      double alp, psi, bx, by, bz;
  
      get_metric_general_relativity(xmet,ymet,zmet,&alp,&psi,&bx,&by,&bz);

      SphP[i].VelVertex[0] = SphP[i].VelVertex[0] - bx / alp;
      SphP[i].VelVertex[1] = SphP[i].VelVertex[1] - by / alp;
      SphP[i].VelVertex[2] = SphP[i].VelVertex[2] - bz / alp;     
      
#endif

      double acc[3];

      /*  the actual time-step of particle */
      integertime ti_step = P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0;
      dt = ti_step * All.Timebase_interval;
      dt /= All.cf_hubble_a;    /* this gives the actual timestep: dt = dloga/ (adot/a) */


      /* now let's add the gradient of the pressure force */

      if(SphP[i].Density > 0)
        {
          acc[0] = -SphP[i].Grad.dpress[0] / SphP[i].Density;
          acc[1] = -SphP[i].Grad.dpress[1] / SphP[i].Density;
          acc[2] = -SphP[i].Grad.dpress[2] / SphP[i].Density;

#ifdef MHD
      /* we also add the acceleration due to the Lorentz force */
          acc[0] += (SphP[i].CurlB[1] * SphP[i].B[2] - SphP[i].CurlB[2] * SphP[i].B[1]) / SphP[i].Density; 
          acc[1] += (SphP[i].CurlB[2] * SphP[i].B[0] - SphP[i].CurlB[0] * SphP[i].B[2]) / SphP[i].Density; 
          acc[2] += (SphP[i].CurlB[0] * SphP[i].B[1] - SphP[i].CurlB[1] * SphP[i].B[0]) / SphP[i].Density; 

#endif

          SphP[i].VelVertex[0] += 0.5 * dt * acc[0];
          SphP[i].VelVertex[1] += 0.5 * dt * acc[1];
          SphP[i].VelVertex[2] += 0.5 * dt * acc[2];
        }

#ifdef ACTIVE_CELL_SPIN
      SphP[i].VelVertex[0] += SphP[i].CenterOffset[1] * SphP[i].Omega[2] - SphP[i].CenterOffset[2] * SphP[i].Omega[1];
      SphP[i].VelVertex[1] += SphP[i].CenterOffset[2] * SphP[i].Omega[0] - SphP[i].CenterOffset[0] * SphP[i].Omega[2];
      SphP[i].VelVertex[2] += SphP[i].CenterOffset[0] * SphP[i].Omega[1] - SphP[i].CenterOffset[1] * SphP[i].Omega[0];
#endif
    }

#if (defined(BOUNDARY_STICKY_MINID) && defined(BOUNDARY_STICKY_MAXID)) || defined(STICKYFLAGS)
  sticky_boundary_vertex_velocities();
#endif

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

#ifdef REGULARIZE_MESH_CM_DRIFT

      double dx, dy, dz, d, fraction;

      dx = nearest_x(P[i].Pos[0] - SphP[i].Center[0]);
      dy = nearest_y(P[i].Pos[1] - SphP[i].Center[1]);
      dz = nearest_z(P[i].Pos[2] - SphP[i].Center[2]);

      /*  the actual time-step of particle */
      dt = (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;
      dt /= All.cf_hubble_a;    /* this is dt, the actual timestep  */

      double cellrad = get_cell_radius(i);

#if !defined(REGULARIZE_MESH_FACE_ANGLE)
      /* if there is a density gradient, use a center that is displaced slightly in the direction of the gradient.
       * This makes sure that the Lloyd scheme does not simply iterate towards cells of equal volume, instead
       * we keep cells of roughly equal mass.
       */
      double dgrad = sqrt(SphP[i].Grad.drho[0] * SphP[i].Grad.drho[0] + SphP[i].Grad.drho[1] * SphP[i].Grad.drho[1] + SphP[i].Grad.drho[2] * SphP[i].Grad.drho[2]);

      if(dgrad > 0)
        {
          double scale = SphP[i].Density / dgrad;
          double tmp = 3 * cellrad + scale;
          double x = (tmp - sqrt(tmp * tmp - 8 * cellrad * cellrad)) / 4;

          if(x < 0.25 * cellrad)
            {
              dx = nearest_x(P[i].Pos[0] - (SphP[i].Center[0] + x * SphP[i].Grad.drho[0] / dgrad));
              dy = nearest_y(P[i].Pos[1] - (SphP[i].Center[1] + x * SphP[i].Grad.drho[1] / dgrad));
              dz = nearest_z(P[i].Pos[2] - (SphP[i].Center[2] + x * SphP[i].Grad.drho[2] / dgrad));
            }
        }
#endif


      d = sqrt(dx * dx + dy * dy + dz * dz);

      fraction = 0;

#if !defined(REGULARIZE_MESH_FACE_ANGLE)
      if(d > 0.75 * All.CellShapingFactor * cellrad && dt > 0)
        {
          if(d > All.CellShapingFactor * cellrad)
            fraction = All.CellShapingSpeed;
          else
            fraction = All.CellShapingSpeed * (d - 0.75 * All.CellShapingFactor * cellrad) / (0.25 * All.CellShapingFactor * cellrad);
        }
#else
      if(SphP[i].MaxFaceAngle > 0.75 * All.CellMaxAngleFactor && dt > 0)
        {
          if(SphP[i].MaxFaceAngle > All.CellMaxAngleFactor)
            fraction = All.CellShapingSpeed;
          else
            fraction = All.CellShapingSpeed * (SphP[i].MaxFaceAngle - 0.75 * All.CellMaxAngleFactor) / (0.25 * All.CellMaxAngleFactor);
        }
#endif

      if(d > 0 && fraction > 0)
        {
          double v;
#ifdef REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED

          v = All.cf_atime * get_sound_speed(i);

#if defined(SELFGRAVITY) || defined(EXTERNALGRAVITY) || defined(EXACT_GRAVITY_FOR_PARTICLE_TYPE)
          /* calculate gravitational velocity scale */
          double ax, ay, az, ac, vgrav;
#ifdef HIERARCHICAL_GRAVITY
          ax = SphP[i].FullGravAccel[0];
          ay = SphP[i].FullGravAccel[1];
          az = SphP[i].FullGravAccel[2];
#else
          ax = P[i].GravAccel[0];
          ay = P[i].GravAccel[1];
          az = P[i].GravAccel[2];
#endif
#ifdef PMGRID
          ax += P[i].GravPM[0];
          ay += P[i].GravPM[1];
          az += P[i].GravPM[2];
#endif
          ac = sqrt(ax * ax + ay * ay + az * az);
          vgrav = 4 * sqrt(All.cf_atime * cellrad * ac);
          if(v < vgrav)
            v = vgrav;
#endif

          double vcurl = cellrad * SphP[i].CurlVel;
          if(v < vcurl)
            v = vcurl;

#else
          v = All.cf_atime * All.cf_atime * d / dt;     /* use fiducial velocity */

          double vel = sqrt(P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);
          double vmax = dmax(All.cf_atime * get_sound_speed(i), vel);
          if(v > vmax)
            v = vmax;
#endif

#ifdef REFINEMENT_SPLIT_CELLS
          double proj = SphP[i].SepVector[0] * dx + SphP[i].SepVector[1] * dy + SphP[i].SepVector[2] * dz;

          if(proj != 0)
            {
              dx = proj * SphP[i].SepVector[0];
              dy = proj * SphP[i].SepVector[1];
              dz = proj * SphP[i].SepVector[2];
            }

          SphP[i].SepVector[0] = 0;
          SphP[i].SepVector[1] = 0;
          SphP[i].SepVector[2] = 0;
#endif

#ifdef REFINEMENT_MERGE_PAIRS
          /* do not do any grid regularization when we are actually trying to move two
             mesh generating points very close together */
          if(SphP[i].DerefPartnerId > 0)
            fraction = 0;
#endif

          SphP[i].VelVertex[0] += fraction * v * (-dx / d);
          SphP[i].VelVertex[1] += fraction * v * (-dy / d);
          SphP[i].VelVertex[2] += fraction * v * (-dz / d);
        }
#endif

#if defined(BOUNDARY_INFLOWOUTFLOW_MINID) && defined(BOUNDARY_INFLOWOUTFLOW_MAXID)
      if(P[i].ID >= BOUNDARY_INFLOWOUTFLOW_MINID && P[i].ID < BOUNDARY_INFLOWOUTFLOW_MAXID)
        {
          SphP[i].VelVertex[0] = 0;
          SphP[i].VelVertex[1] = 0;
          SphP[i].VelVertex[2] = 0;
        }
#endif

#if defined(BOUNDARY_STICKY_MINID) && defined(BOUNDARY_STICKY_MAXID)
      if(P[i].ID >= BOUNDARY_STICKY_MINID && P[i].ID < BOUNDARY_STICKY_MAXID)
        {
          SphP[i].VelVertex[0] = 0;
          SphP[i].VelVertex[1] = 0;
          SphP[i].VelVertex[2] = 0;
        }
#endif
#ifdef STICKYFLAGS

      if(SphP[i].StickyFlag > 0)
        {
          SphP[i].VelVertex[0] = 0;
          SphP[i].VelVertex[1] = 0;
          SphP[i].VelVertex[2] = 0;
        }
      
      /*if outflow boundaries are in use, make non-sticky cells within sticky layer
         stationary */
      double boundary_dist;
      boundary_dist = fmax(boxSize_X, fmax(boxSize_Y, boxSize_Z));

#if (REFLECTIVE_X == 2)
      boundary_dist = fmin(boundary_dist, fmin(boxSize_X - P[i].Pos[0], P[i].Pos[0]));
#endif
#if (REFLECTIVE_Y == 2)
      boundary_dist = fmin(boundary_dist, fmin(boxSize_Y - P[i].Pos[1], P[i].Pos[1]));
#endif
#if (REFLECTIVE_Z == 2)
      boundary_dist = fmin(boundary_dist, fmin(boxSize_Z - P[i].Pos[2], P[i].Pos[2]));
#endif

      if(boundary_dist < All.StickyLayerMaxDist)
        {
          SphP[i].VelVertex[0] = 0;
          SphP[i].VelVertex[1] = 0;
          SphP[i].VelVertex[2] = 0;
        }
#endif

#if defined(BOUNDARY_REFL_FLUIDSIDE_MINID) && defined(BOUNDARY_REFL_FLUIDSIDE_MAXID)
      if(P[i].ID >= BOUNDARY_REFL_FLUIDSIDE_MINID && P[i].ID < BOUNDARY_REFL_FLUIDSIDE_MAXID)
        {
          SphP[i].VelVertex[0] = 0;
          SphP[i].VelVertex[1] = 0;
          SphP[i].VelVertex[2] = 0;
        }
#endif

#if defined(BOUNDARY_REFL_SOLIDSIDE_MINID) && defined(BOUNDARY_REFL_SOLIDSIDE_MAXID)
      if(P[i].ID >= BOUNDARY_REFL_SOLIDSIDE_MINID && P[i].ID < BOUNDARY_REFL_SOLIDSIDE_MAXID)
        {
          SphP[i].VelVertex[0] = 0;
          SphP[i].VelVertex[1] = 0;
          SphP[i].VelVertex[2] = 0;
        }
#endif

#ifdef SPECIAL_BOUNDARY
      boundary_overide_velocities(P, SphP, i);
#endif

      for(j = NUMDIMS; j < 3; j++)
        SphP[i].VelVertex[j] = 0;       /* vertex velocities for unused dimensions set to zero */
    }

#ifdef COFFEE_PROBLEM
  coffee_overide_velocities();
#endif

#ifdef REGULARIZE_MESH_LLOYD
  set_vertex_velocities_lloyd();
#endif

#ifdef REGULARIZE_MESH_SMOOTH
  set_vertex_velocities_smooth();
#endif

#ifdef OUTPUT_VERTEX_VELOCITY_DIVERGENCE
  voronoi_exchange_primitive_variables();
  calculate_vertex_velocity_divergence();
#endif

#ifdef REFINEMENT_MERGE_PAIRS
  voronoi_derefinement_pairs_velvertex_corrections();
#endif


#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
  validate_vertex_velocities();
#endif

#ifdef ONEDIMS_SPHERICAL
  validate_vertex_velocities_1d();
#endif

  TIMER_STOP(CPU_SET_VERTEXVELS);
}

#ifdef ONEDIMS_SPHERICAL
void validate_vertex_velocities_1d()
{
  double dt = (P[0].TimeBinHydro ? (((integertime) 1) << P[0].TimeBinHydro) : 0) * All.Timebase_interval;
  if(P[0].Pos[0] + dt * SphP[0].VelVertex[0] < All.CoreRadius)
    SphP[0].VelVertex[0] = 0.;
}
#endif

#if (defined(BOUNDARY_STICKY_MINID) && defined(BOUNDARY_STICKY_MAXID)) || defined(STICKYFLAGS)
void sticky_boundary_vertex_velocities(void)
{
  int i, k, particle, p1, p2;

  point *DP = Mesh.DP;
  face *VF = Mesh.VF;

  for(i = 0; i < Mesh.Nvf; i++)
    {
      for(k = 0; k < 2; k++)
        {
          if(k == 0)
            {
              p1 = VF[i].p1;
              p2 = VF[i].p2;
            }
          else
            {
              p1 = VF[i].p2;
              p2 = VF[i].p1;
            }

          if(DP[p1].task != ThisTask)
            continue;

          particle = DP[p1].index;

          if(particle < 0 || particle >= NumGas)
            continue;

          if(!TimeBinSynchronized[P[particle].TimeBinHydro])
            continue;

	  
#if defined(BOUNDARY_STICKY_MINID) && defined(BOUNDARY_STICKY_MAXID)
          MyIDType ID1 = DP[p1].ID;
          MyIDType ID2 = DP[p2].ID;

          if((ID1 >= BOUNDARY_STICKY_MINID && ID1 < BOUNDARY_STICKY_MAXID) || (ID2 >= BOUNDARY_STICKY_MINID && ID2 < BOUNDARY_STICKY_MAXID))
            {
              SphP[particle].VelVertex[0] = 0;
              SphP[particle].VelVertex[1] = 0;
              SphP[particle].VelVertex[2] = 0;
            }

#endif

#ifdef STICKYFLAGS
          
	  int index1 = DP[p1].index;
	  int index2 = DP[p2].index;

          if((SphP[index1].StickyFlag > 0) || (SphP[index2].StickyFlag > 0))
            {
              SphP[particle].VelVertex[0] = 0;
              SphP[particle].VelVertex[1] = 0;
              SphP[particle].VelVertex[2] = 0;
            }
#endif

        }
    }
}
#endif


#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
/* in case we have reflecting boundaries, make sure that we do not drift points beyond boundary */
void validate_vertex_velocities(void)
{
  int idx, i;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      integertime ti_step = P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0;
      double dt_drift;

      if(All.ComovingIntegrationOn)
        dt_drift = get_drift_factor(All.Ti_Current, All.Ti_Current + ti_step);
      else
        dt_drift = ti_step * All.Timebase_interval;

#if defined(REFLECTIVE_X)
      if((P[i].Pos[0] + dt_drift * SphP[i].VelVertex[0]) < 0 || (P[i].Pos[0] + dt_drift * SphP[i].VelVertex[0]) >= boxSize_X)
        SphP[i].VelVertex[0] = 0;
#endif
#if defined(REFLECTIVE_Y)
      if((P[i].Pos[1] + dt_drift * SphP[i].VelVertex[1]) < 0 || (P[i].Pos[1] + dt_drift * SphP[i].VelVertex[1]) >= boxSize_Y)
        SphP[i].VelVertex[1] = 0;
#endif
#if defined(REFLECTIVE_Z)
      if((P[i].Pos[2] + dt_drift * SphP[i].VelVertex[2]) < 0 || (P[i].Pos[2] + dt_drift * SphP[i].VelVertex[2]) >= boxSize_Z)
        SphP[i].VelVertex[2] = 0;
#endif
    }
}
#endif
