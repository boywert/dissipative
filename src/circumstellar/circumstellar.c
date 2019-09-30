/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/circumstellar/circumstellar.c
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

#include "../allvars.h"
#include "../proto.h"

#include "./circumstellar_proto.h"

#ifdef CIRCUMSTELLAR
#if (defined(CIRCUMSTELLAR_IRRADIATION) || defined(ALPHA_VISCOSITY) || defined(CIRCUMSTELLAR_REFINEMENTS)) && !defined (EXTERNALGRAVITY)
void source_particle_create_list()
{
  struct source_particle_data *SourcePartList;
  SourcePartList = (struct source_particle_data *) mymalloc("SourcePartList", All.MaxPartSources * sizeof(struct source_particle_data));


  int i, j, nsrc, nimport, ngrp;
  for(i = 0, nsrc = 0; i < NumPart; i++)
    {
      if(P[i].Type == 4 || P[i].Type == 5)
        {
          SourcePartList[nsrc].SofteningType = P[i].SofteningType;
          SourcePartList[nsrc].SourceID = P[i].ID;
          SourcePartList[nsrc].Mass = P[i].Mass;

          SourcePartList[nsrc].Pos[0] = P[i].Pos[0];
          SourcePartList[nsrc].Pos[1] = P[i].Pos[1];
          SourcePartList[nsrc].Pos[2] = P[i].Pos[2];

          if(P[i].Mass > 0.08 * SOLAR_MASS / All.UnitMass_in_g)
            SourcePartList[nsrc++].SurfaceTemp = SOLAR_EFF_TEMP * sqrt(P[i].Mass * All.UnitMass_in_g / SOLAR_MASS);
          else
            SourcePartList[nsrc++].SurfaceTemp = 0;
        }
    }

  for(j = 0; j < NTask; j++)
    Send_count[j] = nsrc;

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = 0;
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  /* exchange particle data */
  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&SourcePartList[Send_offset[recvTask]],
                           Send_count[recvTask] * sizeof(struct source_particle_data), MPI_BYTE,
                           recvTask, TAG_DENS_A,
                           &SourcePartListGlobal[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct source_particle_data), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  myfree(SourcePartList);

}

void source_particle_update_list()
{
  struct source_particle_data *SourcePartList;
  SourcePartList = (struct source_particle_data *) mymalloc("SourcePartList", All.MaxPartSources * sizeof(struct source_particle_data));

  int i, j, nsrc, nimport, ngrp;
  for(i = 0, nsrc = 0; i < NumPart; i++)
    {
      if(P[i].Type == 4 || P[i].Type == 5)
        {
          SourcePartList[nsrc].SofteningType = P[i].SofteningType;
          SourcePartList[nsrc].SourceID = P[i].ID;
          SourcePartList[nsrc].Mass = P[i].Mass;

          SourcePartList[nsrc].Pos[0] = P[i].Pos[0];
          SourcePartList[nsrc].Pos[1] = P[i].Pos[1];
          SourcePartList[nsrc].Pos[2] = P[i].Pos[2];
          if(P[i].Mass > 0.08 * SOLAR_MASS / All.UnitMass_in_g)
            SourcePartList[nsrc++].SurfaceTemp = SOLAR_EFF_TEMP * sqrt(P[i].Mass * All.UnitMass_in_g / SOLAR_MASS);
          else
            SourcePartList[nsrc++].SurfaceTemp = 0;
        }
    }

  for(j = 0; j < NTask; j++)
    Send_count[j] = nsrc;

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = 0;
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  /* exchange particle data */
  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&SourcePartList[Send_offset[recvTask]],
                           Send_count[recvTask] * sizeof(struct source_particle_data), MPI_BYTE,
                           recvTask, TAG_DENS_A,
                           &SourcePartListGlobal[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct source_particle_data), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  myfree(SourcePartList);

}


double get_circumstellar_distance(int i)
{
  double dx, dy, dz, r, r2;
  double min_r2 = +MAX_REAL_NUMBER;
  int k;
  for(k = 0; k < All.TotPartSources; k++)
    {
      if(SourcePartListGlobal[k].SurfaceTemp == 0)
        continue;
      dx = P[i].Pos[0] - SourcePartListGlobal[k].Pos[0];
      dy = P[i].Pos[1] - SourcePartListGlobal[k].Pos[1];
      dz = P[i].Pos[2] - SourcePartListGlobal[k].Pos[2];

      r2 = dx * dx + dy * dy + dz * dz;
      if(r2 < min_r2)
        min_r2 = r2;
    }
  r = sqrt(min_r2);
  return r;
}
#endif

#ifdef CIRCUMSTELLAR_IRRADIATION
double circumstellar_DoCoolingHeating(int i, double dtime)
{
  if(All.Time == All.TimeBegin)
    return 0;
  if(SphP[i].Density < All.LimitUBelowThisDensity)
    return 0;

  double eq_utherm = circumstellar_irradiation_utherm(i);
  double delta = eq_utherm - SphP[i].Utherm;
  return delta;
}

double circumstellar_irradiation_utherm(int i)
{
  double wp;
  double dx, dy, dz, r, r2;
  double h, h_inv, h3_inv, u;
  int k;
  double min_r2 = +MAX_REAL_NUMBER;
  int ParentStarIndex = 0;
  double temp4, temp, utherm;

  for(k = 0; k < All.TotPartSources; k++)
    {
      if(SourcePartListGlobal[k].SurfaceTemp == 0)
        continue;
      dx = P[i].Pos[0] - SourcePartListGlobal[k].Pos[0];
      dy = P[i].Pos[1] - SourcePartListGlobal[k].Pos[1];
      dz = P[i].Pos[2] - SourcePartListGlobal[k].Pos[2];

      r2 = dx * dx + dy * dy + dz * dz;
      if(r2 < min_r2)
        {
          min_r2 = r2;
          ParentStarIndex = k;
        }
    }
  h = All.ForceSoftening[SourcePartListGlobal[ParentStarIndex].SofteningType];
  h_inv = 1.0 / h;
  r = sqrt(min_r2);
  r2 = min_r2;

  //using spline softening
  if(r >= h)
    {
      wp = 1 / r;
    }
  else
    {
      h3_inv = h_inv * h_inv * h_inv;
      u = r * h_inv;

      if(u < 0.5)
        wp = -h_inv * (-2.8 + u * u * (5.333333333333 + u * u * (6.4 * u - 9.6)));
      else
        wp = -h_inv * (-3.2 + 0.066666666667 / u + u * u * (10.666666666667 + u * (-16.0 + u * (9.6 - 2.133333333333 * u))));

    }

  temp4 = pow(SourcePartListGlobal[ParentStarIndex].SurfaceTemp, 4) * wp * wp / h_inv / h_inv / 2.8 / 2.8;
  temp = All.IrradiationTempScaling * pow(temp4, 0.25);

  if(temp > 0)
    utherm = circumstellar_convert_temp_to_u(temp);
  else
    utherm = SphP[i].Utherm;

  return utherm;

}

double circumstellar_convert_temp_to_u(double T)
{
  double X = 0.7;
  double Y = 0.28;
  double Z = 0.02;
  double meanweight = 1.0 / (0.5 * X + Y / 4.002602 + 0.5 / 17.0 * Z);
  return BOLTZMANN * T / GAMMA_MINUS1 / meanweight / PROTONMASS / All.UnitEnergy_in_cgs * All.UnitMass_in_g;
}
#endif


double get_circumstellar_alpha_viscosity(double x, double y, double z, double rho, double press)
{
#if defined(VISCOSITY) && defined(ALPHA_VISCOSITY)

  double dx = 0, dy = 0, dz = 0;
  double r, r2;
  double OmegaKep, Csnd;

#ifndef EXTERNALGRAVITY
  double h, h_inv, h3_inv, u;
  double fac;
  int k;
  double min_r2 = +MAX_REAL_NUMBER;
  int ParentStarIndex = 0;
  for(k = 0; k < All.TotPartSources; k++)
    {
      if(SourcePartListGlobal[k].SurfaceTemp == 0)
        continue;
      dx = x - SourcePartListGlobal[k].Pos[0];
      dy = y - SourcePartListGlobal[k].Pos[1];
      dz = z - SourcePartListGlobal[k].Pos[2];

      r2 = dx * dx + dy * dy + dz * dz;
      if(r2 < min_r2)
        {
          min_r2 = r2;
          ParentStarIndex = k;
        }
    }

  r = sqrt(r2);
  h = All.ForceSoftning[SourcePartListGlobal[ParentStarIndex].SofteningType];

  //using spline softening
  if(r >= h)
    fac = 1 / (r2 * r);
  else
    {
      h_inv = 1.0 / h;
      h3_inv = h_inv * h_inv * h_inv;
      u = r * h_inv;
      if(u < 0.5)
        fac = h3_inv * (10.666666666667 + u * u * (32.0 * u - 38.4));
      else
        fac = h3_inv * (21.333333333333 - 48.0 * u + 38.4 * u * u - 10.666666666667 * u * u * u - 0.066666666667 / (u * u * u));
    }

  OmegaKep = sqrt(All.G * SourcePartListGlobal[ParentStarIndex].Mass * fac);
#else
  dx = x - 0.5 * boxSize_X;
  dy = y - 0.5 * boxSize_Y;
#ifndef TWODIMS
  dz = z - 0.5 * boxSize_Z;
#endif

  r2 = dx * dx + dy * dy + dz * dz;
  r = sqrt(r2);

  OmegaKep = sqrt(All.G / r2 / r);
#endif

#ifndef ISOTHERM_EQS
  Csnd = GAMMA * press / rho;
#else
  Csnd = All.IsoSoundSpeed;
#endif

  return All.AlphaCoefficient * Csnd * Csnd / OmegaKep;

#endif
  return 0;
}


#ifdef CIRCUMSTELLAR_WBOUNDARIES
void do_circumstellar_disk_update(struct particle_data *localP, struct sph_particle_data *localSphP, int p)
{

  double x, y, radius, cosphi, sinphi;
  double Tin, Tout, damping, lambda;
  double DRin, DRout;
  double omega_keplerian, velx0, vely0, velz0;
  double dens0, press0;
  double CentralMass;
  double dt, aphys;
  double delta_rho, delta_velx, delta_vely, delta_velz, delta_press;

#ifdef CENTRAL_MASS_POTENTIAL
  CentralMass = All.CentralMass;
#else
#ifdef BINARY_POTENTIAL
  CentralMass = 1.0 - (All.BinaryMassRatio) / (1 + All.BinaryMassRatio);
#else
  CentralMass = 1.0;
#endif
#endif


  if(All.EvanescentBoundaryStrength <= 0.01)
    return;

  x = localP[p].Pos[0] - boxSize_X / 2.0;
  y = localP[p].Pos[1] - boxSize_Y / 2.0;
  radius = sqrt(x * x + y * y);
  cosphi = x / radius;
  sinphi = y / radius;
  omega_keplerian = sqrt(All.G * CentralMass / radius / radius / radius);

  DRin = (All.inner_radius) * 5.0;
  DRout = (All.outer_radius) * 0.84;
  if((radius > All.outer_radius) || (radius < All.inner_radius))
    return;
  if((radius > DRin) && (radius < DRout))
    return;

  Tin = 2.0 * M_PI * All.inner_radius * sqrt(All.inner_radius / All.G / CentralMass);
  Tout = 2.0 * M_PI * All.outer_radius * sqrt(All.outer_radius / All.G / CentralMass);


  dt = get_timestep(p, &aphys, 0) * All.Timebase_interval;      /* TODO call get_timestep_hydro now */
  //dt = get_cell_radius(p) / get_sound_speed(p);
  //dt = get_cell_radius(p) * sqrt(radius/ All.G / CentralMass);
  //dt = (localP[p].TimeBinHydro ? (((integertime) 1) << localP[p].TimeBinHydro) : 0) * All.Timebase_interval;

  //These are the ideal steady-state quanitities
  dens0 = All.CircumstellarBoundaryDensity;
  velx0 = -sinphi * radius * omega_keplerian;
  vely0 = cosphi * radius * omega_keplerian;
  velz0 = 0;
  press0 = All.CircumstellarBoundaryDensity * All.IsoSoundSpeed * All.IsoSoundSpeed;
#ifdef LOCALLY_ISOTHERM
  press0 = All.CircumstellarBoundaryDensity * All.IsoSoundSpeed * All.IsoSoundSpeed / radius;
#endif

  if(radius < DRin)
    {
      damping = (DRin - radius) / (DRin - All.inner_radius);
      //damping = sin(0.5 * M_PI * (DRin - radius)/(DRin - All.inner_radius));
      //damping = exp(-8 * (radius - All.inner_radius) * (radius - All.inner_radius)/(DRin - All.inner_radius) / (DRin - All.inner_radius));
      lambda = damping * damping * dt / Tin;
      //lambda = 0;
      //lambda = damping * dt / Tin;
    }
  if(radius > DRout)
    {
      damping = (radius - DRout) / (All.outer_radius - DRout);
      //damping = sin(0.5 * M_PI * (radius - DRout)/(All.outer_radius - DRout));
      //damping = exp(-8 * (radius - All.outer_radius) * (radius - All.outer_radius)/(DRout - All.outer_radius) / (DRout - All.outer_radius));

      lambda = damping * damping * dt / Tout;
      //lambda = damping * dt / Tout;
    }

  //lambda = damping * damping * get_cell_radius(p) / radius;
  //lambda = damping *  get_cell_radius(p) / radius;

  lambda *= All.EvanescentBoundaryStrength;

  delta_rho = (dens0 - localSphP[p].Density) * lambda;
  delta_velx = (velx0 - localP[p].Vel[0]) * lambda;
  delta_vely = (vely0 - localP[p].Vel[1]) * lambda;
  delta_velz = (velz0 - localP[p].Vel[2]) * lambda;
  delta_press = (press0 - localSphP[p].Pressure) * lambda;

  /*we update the primitive variables */
  //localSphP[p].Density = (localSphP[p].Density + lambda * dens0) / (1.0 + lambda);

  //localP[p].Vel[0] = (localP[p].Vel[0] + lambda * velx0) / (1.0 + lambda);
  //localP[p].Vel[1] = (localP[p].Vel[1] + lambda * vely0) / (1.0 + lambda);
  //localP[p].Vel[2] = (localP[p].Vel[2] + lambda * velz0) / (1.0 + lambda);

  //localSphP[p].Pressure = (localSphP[p].Pressure + lambda * press0) / (1.0 + lambda);

  localSphP[p].Density += delta_rho;

  localP[p].Vel[0] += delta_velx;
  localP[p].Vel[1] += delta_vely;
  localP[p].Vel[2] += delta_velz;

  localSphP[p].Pressure += delta_press;

  /*Now we correct the conserved quantities */
  localP[p].Mass = localSphP[p].Density * localSphP[p].Volume;
  localSphP[p].Momentum[0] = localP[p].Mass * localP[p].Vel[0];
  localSphP[p].Momentum[1] = localP[p].Mass * localP[p].Vel[1];
  localSphP[p].Momentum[2] = localP[p].Mass * localP[p].Vel[2];

#ifndef ISOTHERM_EQS
  localSphP[p].Utherm = localSphP[p].Pressure / localSphP[p].Density / (GAMMA - 1);
#endif

}
#endif


#endif
