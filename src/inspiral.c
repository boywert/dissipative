/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/inspiral.c
 * \date        06/2016
 * \author      R. Pakmor
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

#ifdef INSPIRAL

struct wd {
  double center[3];
  double vel[3];
  double mass;
};

static void get_wd_properties( struct wd *wd, int ipass )
{
  double mass = 0;
  double center[3], mom[3];
  for(int k=0; k<3; k++)
    {
      center[k] = 0;
      mom[k] = 0;
    }
  
  for(int i=0; i<NumGas; i++)
    {
      if(P[i].Type != 0 || (P[i].Mass == 0 && P[i].ID == 0))
        continue;
      
      double cmass = P[i].Mass * SphP[i].PScalars[ipass];
      mass += cmass;
      for(int k=0; k<3; k++)
        {
          center[k] += SphP[i].Center[k] * cmass;
          mom[k] += P[i].Vel[k] * cmass;
        }
    }
  
  double mtot, centertot[3], momtot[3];
  MPI_Allreduce(&mass, &mtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(center, centertot, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(mom, momtot, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  wd->mass = mtot;
  if(wd->mass > 0)
    {
      for(int k=0; k<3; k++)
        {
          wd->center[k] = centertot[k] / wd->mass;
          wd->vel[k] = momtot[k] / wd->mass;
        }
    }
  else
    {
      for(int k=0; k<3; k++)
        {
          wd->center[k] = 0;
          wd->vel[k] = 0;
        }
    }
}

static void inspiral_compute_acceleration()
{
  struct wd wd1, wd2;
  get_wd_properties( &wd1, 0 );
  get_wd_properties( &wd2, 1 );

  double mtot = wd1.mass + wd2.mass;

  double dx = wd1.center[0] - wd2.center[0];
  double dy = wd1.center[1] - wd2.center[1];
  double dz = wd1.center[2] - wd2.center[2];
  double dist = sqrt( dx*dx + dy*dy + dz*dz );
  
  double fac = All.G * All.InspiralVelocity / (2.*dist*dist*mtot);
  
  double m12 = wd1.mass * wd1.mass;
  double m22 = wd2.mass * wd2.mass;
  
  double v12 = wd1.vel[0]*wd1.vel[0] + wd1.vel[1]*wd1.vel[1] + wd1.vel[2]*wd1.vel[2];
  double v22 = wd2.vel[0]*wd2.vel[0] + wd2.vel[1]*wd2.vel[1] + wd2.vel[2]*wd2.vel[2];
  
  for(int k=0; k<3; k++)
    {
      All.WD1_Acc[k] = - m22 * fac * wd1.vel[k] / v12;
      All.WD2_Acc[k] = - m12 * fac * wd2.vel[k] / v22;
    }
  
  if(ThisTask == 0)
    {
      char buf[1000];
      sprintf(buf, "%sbinary.txt", All.OutputDir);
      
      FILE *fp;
      if(All.Time == 0)
        fp = fopen( buf, "w" );
      else
        fp = fopen( buf, "a" );
      fprintf( fp, "%g %g %g %g %g %g %g %g %g %g\n", All.Time, dist, wd1.mass, wd2.mass, wd1.center[0], wd1.center[1], wd1.center[2],
        wd2.center[0], wd2.center[1], wd2.center[2] );
      fclose(fp);
    }
}

static void apply_friction_force()
{
  int idx;
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      int i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      double dt_cell = 0.5 * (P[i].TimeBinGrav ? (((integertime) 1) << P[i].TimeBinGrav) : 0) * All.Timebase_interval;
      
      if(P[i].Type == 0)
        {
          double ekin_old = 0.5 * (SphP[i].Momentum[0]*SphP[i].Momentum[0] + SphP[i].Momentum[1]*SphP[i].Momentum[1] + 
            SphP[i].Momentum[2]*SphP[i].Momentum[2]) / P[i].Mass;
          
          for(int k=0; k<3; k++)
            {
              double dvel = dt_cell * (All.WD1_Acc[k] * SphP[i].PScalars[0] + All.WD2_Acc[k] * SphP[i].PScalars[1]);
              P[i].Vel[k] += dvel;
              SphP[i].Momentum[k] += dvel * P[i].Mass;
            }
          
          double ekin_new = 0.5 * (SphP[i].Momentum[0]*SphP[i].Momentum[0] + SphP[i].Momentum[1]*SphP[i].Momentum[1] + 
            SphP[i].Momentum[2]*SphP[i].Momentum[2]) / P[i].Mass;
          
          SphP[i].Energy += ekin_new - ekin_old;
        }
    }
}

void do_inspiral_source_terms_first_half()
{
  apply_friction_force();
}

void do_inspiral_source_terms_second_half()
{
  if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
    inspiral_compute_acceleration();
  
  apply_friction_force();
}

#endif
