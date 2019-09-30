/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/boundaries.c
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

#ifdef SPECIAL_BOUNDARY
void boundary_overide_velocities(struct particle_data *localP, struct sph_particle_data *localSphP, int p)
{
  double vx, vy, vz;
  double dt;
  double weight = 0;

#ifdef WINDTUNNEL
  if(localP[p].ID >= BOUNDARY_INFLOWOUTFLOW_MINID && localP[p].ID < BOUNDARY_INFLOWOUTFLOW_MAXID)
    {
      localSphP[p].VelVertex[0] = 0;
      localSphP[p].VelVertex[1] = 0;
      localSphP[p].VelVertex[2] = 0;
      return;
    }
#endif
  if((localP[p].ID == -1) || (localP[p].ID == -2) || (localP[p].ID <= -3))
    {
      if(localP[p].ID <= -3)
        {
          localSphP[p].VelVertex[0] = 0;
          localSphP[p].VelVertex[1] = 0;
          localSphP[p].VelVertex[2] = 0;
          return;
        }

      if((localP[p].ID == -1) || (localP[p].ID == -2))
        {
          dt = (localP[p].TimeBinHydro ? (((integertime) 1) << localP[p].TimeBinHydro) : 0) * All.Timebase_interval;
          boundary_get_velocity(localP[p].Pos[0], localP[p].Pos[1], localP[p].Pos[2], &vx, &vy, &vz, dt);

          localSphP[p].VelVertex[0] = vx;
          localSphP[p].VelVertex[1] = vy;
          localSphP[p].VelVertex[2] = vz;
        }
    }
  else
    {
      dt = (localP[p].TimeBin ? (((integertime) 1) << localP[p].TimeBin) : 0) * All.Timebase_interval;
      boundary_get_velocity(localP[p].Pos[0], localP[p].Pos[1], localP[p].Pos[2], &vx, &vy, &vz, dt);

      weight = usr_defined_get_damping_weight(p);


      localSphP[p].VelVertex[0] = localSphP[p].VelVertex[0] * (1.0 - weight) + weight * vx;
      localSphP[p].VelVertex[1] = localSphP[p].VelVertex[1] * (1.0 - weight) + weight * vy;
      localSphP[p].VelVertex[2] = localSphP[p].VelVertex[2] * (1.0 - weight) + weight * vz;

    }
  return;
}

double usr_defined_get_damping_weight(int i)
{
#ifdef CIRCUMSTELLAR_WBOUNDARIES
  double x, y, r;
  double weight;

  x = P[i].Pos[0];
  y = P[i].Pos[1];

  x -= boxSize_X / 2.0;
  y -= boxSize_Y / 2.0;

  r = sqrt(x * x + y * y);

  if((r < All.inner_radius) || (r > 1.5 * All.inner_radius))
    weight = 0;
  else
    weight = (1.5 * All.inner_radius - r) * (1.5 * All.inner_radius - r) / (0.5 * All.inner_radius) / (0.5 * All.inner_radius);

  return weight;
#endif
  return 0;
}

void boundary_get_velocity(double x, double y, double z, double *vx, double *vy, double *vz, double dt)
{
  switch (All.SpecialBoundaryMotion)
    {
      // Static special boundary
    case 0:
      *vx = 0.0;
      *vy = 0.0;
      *vz = 0.0;
      break;

      // Motion along the x-direction
    case 1:
      boundary_x_motion(x, y, z, vx, vy, vz);
      break;

      // Motion along the y-direction
    case 2:
      boundary_y_motion(x, y, z, vx, vy, vz);
      break;

      // Motion along the z-direction
    case 3:
      boundary_z_motion(x, y, z, vx, vy, vz);
      break;

      // Circular Motion (e.g spoon)
    case 4:
      boundary_circular_motion(x, y, z, vx, vy, vz);
      break;

      // Concentric Boundaries (e.g Taylor Couette)
    case 5:
      boundary_double_concentric_motion(x, y, z, vx, vy, vz);
      break;

      // Disk Boundaries
    case 6:
      boundary_disk_motion(x, y, z, vx, vy, vz, dt);
      break;

    default:
      terminate("invalid boundary motion specified");
      break;
    }
}



void boundary_x_motion(double x, double y, double z, double *vx, double *vy, double *vz)
{
  *vx = All.SpecialBoundarySpeed;
  *vy = 0.0;
  *vz = 0.0;
}


void boundary_y_motion(double x, double y, double z, double *vx, double *vy, double *vz)
{
  *vx = 0.0;
  *vy = All.SpecialBoundarySpeed;
  *vz = 0.0;
}


void boundary_z_motion(double x, double y, double z, double *vx, double *vy, double *vz)
{
  *vx = 0.0;
  *vy = 0.0;
  *vz = All.SpecialBoundarySpeed;
}


void boundary_circular_motion(double x, double y, double z, double *vx, double *vy, double *vz)
{
  double dt, r, phi, omega;

  x -= 0.5;
  y -= 0.5;

  phi = atan2(y, x);
  r = sqrt(x * x + y * y);

  /*omega = All.BoundaryVel; */

  omega = 2 * M_PI / 5;         /* angular velocity */

  /*  the actual time-step of the particle */
  dt = (P[0].TimeBinHydro ? (((integertime) 1) << P[0].TimeBinHydro) : 0) * All.Timebase_interval;

  if(P[0].TimeBinHydro == 0)
    {
      dt = 1.0e-5;
    }

  phi += omega * dt;

  *vx = (r * cos(phi) - x) / dt;
  *vy = (r * sin(phi) - y) / dt;
  *vz = 0.0;
}


void boundary_double_concentric_motion(double x, double y, double z, double *vx, double *vy, double *vz)
{
#ifdef COAXIAL_BOUNDARIES
  double dt, r, phi, omega = 0;

  /*  the actual time-step of the particle */
  dt = (P[0].TimeBinHydro ? (((integertime) 1) << P[0].TimeBinHydro) : 0) * All.Timebase_interval;

  if(P[0].TimeBinHydro == 0)
    {
      dt = 1.0e-5;
    }

  x -= boxSize_X / 2.0;
  y -= boxSize_Y / 2.0;

  phi = atan2(y, x);
  r = sqrt(x * x + y * y);

  if(r >= 0.7 * All.outer_radius)
    omega = All.omega_out;

  if(r <= 1.5 * All.inner_radius)
    omega = All.omega_in;

  phi += omega * dt;

  *vx = (r * cos(phi) - x) / dt;
  *vy = (r * sin(phi) - y) / dt;
  *vz = 0.0;
#endif
}

void boundary_disk_motion(double x, double y, double z, double *vx, double *vy, double *vz, double dt)
{
#ifdef CIRCUMSTELLAR_WBOUNDARIES
  double r, phi, omega_sq = 0;
  double CentralMass, boundaryradius;

  if(dt < All.MinSizeTimestep)
    {
      dt = All.MinSizeTimestep;
    }
  //dt = All.Timebase_interval;

  //mpi_printf("boundary dt %g\n",dt);

  x -= boxSize_X / 2.0;
  y -= boxSize_Y / 2.0;

  phi = atan2(y, x);
  r = sqrt(x * x + y * y);
#ifdef CENTRAL_MASS_POTENTIAL
  CentralMass = All.CentralMass;
#else
#ifdef BINARY_POTENTIAL
  CentralMass = 1.0;            // - All.BinaryMassRatio / (1 + All.BinaryMassRatio);
#else
#ifdef STAR_PLANET_POTENTIAL
  CentralMass = 1.0 - All.MassRatio;
#else
  CentralMass = 1.0;
#endif
#endif
#endif


  if((r - All.inner_radius) * (r - All.inner_radius) < (r - All.outer_radius) * (r - All.outer_radius))
    boundaryradius = All.inner_radius;
  else
    boundaryradius = All.outer_radius;

  omega_sq = All.G * CentralMass / boundaryradius / boundaryradius / boundaryradius;

#ifdef BINARY_POTENTIAL
  omega_sq += 0.75 * All.G * CentralMass / boundaryradius / boundaryradius / boundaryradius / boundaryradius / boundaryradius
    * All.BinaryMassRatio / (1 + All.BinaryMassRatio) / (1 + All.BinaryMassRatio);
#endif

  phi += sqrt(omega_sq) * dt;

  *vx = (r * cos(phi) - x) / dt;
  *vy = (r * sin(phi) - y) / dt;
  *vz = 0.0;
#endif
}

void get_boundary_cell_state(struct state *state_inside, struct state *state_outside, int orient)
{

  switch (All.SpecialBoundaryType)
    {

      // Reflective Slip Boundary
    case 1:
      outside_state_reflective(state_inside, state_outside);
      break;

      // No-Slip Boundary
    case 2:
      outside_state_noslip(state_inside, state_outside);
      break;

      // Subsonic Outflow Non-Reflective Boundary
    case 3:
      outside_state_nrbc(state_inside, state_outside);
      break;

      // Subsonic Ouflow Open Boundary
    case 4:
      outside_state_pressure(state_inside, state_outside);
      break;

      // Diode-type Ouflow Open Boundary
    case 5:
      outside_state_diode(state_inside, state_outside, orient);
      break;

    default:
      terminate("invalid boundary type specified");
      break;
    }

}


void outside_state_reflective(struct state *state_inside, struct state *state_outside)
{
  state_outside->velx = -state_inside->velx;
  state_outside->vely = state_inside->vely;
  state_outside->velz = state_inside->velz;

  state_outside->rho = state_inside->rho;
  state_outside->press = state_inside->press;
#ifdef TRACER_FIELD
  state_outside->tracer = state_inside->tracer;
#endif
#ifdef VISCOSITY
  int l, k;
  for(l = 0; l < 3; l++)
    for(k = 0; k < 3; k++)
      {
        state_inside->grad->dvel[k][l] = 0.0;
        state_outside->grad->dvel[k][l] = 0.0;
      }
#endif
}


void outside_state_noslip(struct state *state_inside, struct state *state_outside)
{
  state_outside->velx = -state_inside->velx;
  state_outside->vely = -state_inside->vely;
  state_outside->velz = -state_inside->velz;

  state_outside->rho = state_inside->rho;
  state_outside->press = state_inside->press;
#ifdef TRACER_FIELD
  state_outside->tracer = state_inside->tracer;
#endif

#ifdef VISCOSITY
  int l, k;
  for(l = 0; l < 3; l++)
    for(k = 0; k < 3; k++)
      state_outside->grad->dvel[k][l] = state_inside->grad->dvel[k][l];
#endif

}


void outside_state_nrbc(struct state *state_inside, struct state *state_outside)
{
  state_outside->velx = state_inside->velx;
  state_outside->vely = state_inside->vely;
  state_outside->velz = state_inside->velz;

  state_outside->rho = state_inside->rho;
  state_outside->press = state_inside->press;
#ifdef TRACER_FIELD
  state_outside->tracer = state_inside->tracer;
#endif
}



void outside_state_pressure(struct state *state_inside, struct state *state_outside)
{
  state_outside->velx = state_inside->velx;
  state_outside->vely = state_inside->vely;
  state_outside->velz = state_inside->velz;

  state_outside->rho = state_inside->rho;

#ifdef TRACER_FIELD
  state_outside->tracer = state_inside->tracer;
#endif
}

void outside_state_diode(struct state *state_inside, struct state *state_outside, int orient)
{
  state_outside->velx = orient * state_inside->velx;
  state_outside->vely = state_inside->vely;
  state_outside->velz = state_inside->velz;

  state_outside->rho = state_inside->rho;

#ifdef TRACER_FIELD
  state_outside->tracer = state_inside->tracer;
#endif
}

#endif
