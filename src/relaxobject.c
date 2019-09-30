/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/relaxobject.c
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

#include "arepoconfig.h"

#ifdef RELAXOBJECT

#include "allvars.h"
#include "proto.h"

#ifdef RELAXOBJECT_COOLING2
double *wdR, *wdTemp;
int length;
#endif

void relaxobject()
{
  double relaxFac;
  double dt;
  double ekin_old, ekin_new;
  double dvelx, dvely, dvelz;
  int idx, i;

#ifdef RELAXOBJECT_BINARY
  relaxFac = 1. / All.RelaxBaseFac;
#else
  if(All.Time < 0.2 * All.TimeMax)
    {
      relaxFac = 1. / All.RelaxBaseFac;
    }
  else if(All.Time > 0.8 * All.TimeMax)
    {
      relaxFac = 0.;
    }
  else
    {
      relaxFac = 1. / (All.RelaxBaseFac * pow(10., (All.Time - 0.2 * All.TimeMax) / (0.6 * All.TimeMax) * 3.));
    }
#endif
#ifdef RELAX_RUNTIME
  /* the runtime is assumed to be 10 times the dynamical or acoustic time */
  double dyntime = All.TimeMax / 10.;
  /* damping timescales */
  double tau1 = dyntime / 10.;
  double tau2 = dyntime / 1.;
  double t1 = 2.0  * dyntime;
  double t2 = 5.0  * dyntime;

  if (All.Time < t1)
    relaxFac = 1/tau1;
  else if (All.Time < t2)
    {
      double tmp = (All.Time - t1) / (t2 - t1);
      relaxFac = 1. / tau1 * pow(tau1/tau2, tmp);
    }
  else
    relaxFac = 0.;
#endif

  update_primitive_variables();

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(relaxFac > 0.)
        {
          dt = (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;

          ekin_old = .5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);

          dvelx = -dt * relaxFac * P[i].Vel[0];
          dvely = -dt * relaxFac * P[i].Vel[1];
          dvelz = -dt * relaxFac * P[i].Vel[2];

          P[i].Vel[0] += dvelx;
          SphP[i].Momentum[0] += dvelx * P[i].Mass;
          P[i].Vel[1] += dvely;
          SphP[i].Momentum[1] += dvely * P[i].Mass;
          P[i].Vel[2] += dvelz;
          SphP[i].Momentum[2] += dvelz * P[i].Mass;

          ekin_new = .5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);

#ifdef RELAXOBJECT_COOLING
          if(All.Time < 0.2 * All.TimeMax)
            {
              struct eos_result res;
#ifndef RELAXOBJECT_COOLING2
              eos_calc_tgiven(SphP[i].Density, SphP[i].Composition, All.RelaxTemperature, &res);

              SphP[i].Utherm = res.e.v;
#endif

#ifdef RELAXOBJECT_COOLING2
 //used for WD+HeShell relaxation
              double Temp;
              int k;
              if ( P[i].Pos[0] <= All.ShellBaseRadius)
                {
                  Temp = All.TempCore;
                  eos_calc_tgiven( SphP[i].Density, SphP[i].Composition, Temp, &res );
                  SphP[i].Utherm = res.e.v;
                }
// /*              if ( P[i].Pos[0] >= All.ShellBaseRadius) // for small linear Temperature decrease
                {
                  Temp = All.TempShell - abs(P[i].Pos[0] - All.ShellBaseRadius) / abs(All.ShellBaseRadius - All.RTempMin) * abs(All.TempShell - All.TempMin); // for a linear interpolation of the temperature
                } d/

              if ( P[i].Pos[0] < wdR[0])
                {
                  Temp = All.TempCore; 
                }

              for (k=1; k<length; k++)
                {
                  if ( ( P[i].Pos[0] < wdR[k] ) && ( P[i].Pos[0] > wdR[k-1] ) )
                    {
//                    Temp = wdTemp[k];
                      Temp = wdTemp[k-1] + (wdTemp[k-1] - wdTemp[k]) / (wdR[k-1] - wdR[k]) * (P[i].Pos[0] - wdR[k-1]); // linear interpolation between known positions and temperatures
//                      if (Temp < All.BackgroundTemp)
//                        {
//                          Temp = All.BackgroundTemp;
//                        }
                      break;
                    }
                }
              if ( P[i].Pos[0] > wdR[length-1])
                {
                  Temp = All.BackgroundTemp; //wdTemp[length-1];
                }

//            Temp = 5e7;*/
//              eos_calc_tgiven( SphP[i].Density, SphP[i].Composition, Temp, &res );

//              SphP[i].Utherm = res.e.v;         
#endif
            }

          SphP[i].Energy = ekin_new + SphP[i].Utherm * P[i].Mass;
#else
          SphP[i].Energy += ekin_new - ekin_old;
          EgyInjection -= ekin_new - ekin_old;  /* we remove energy from the system... */
#endif
          /* no need to update pressure here, as Utherm does not change */
        }
    }

  update_primitive_variables();
}

#endif

#ifdef RELAXOBJECT_BINARY

static void compute_orbital_period();
static void apply_orbital_acceleration();

void do_binary_source_terms_first_half()
{
  apply_orbital_acceleration();
}

void do_binary_source_terms_second_half()
{
  if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
    compute_orbital_period();

  apply_orbital_acceleration();
}

void compute_orbital_period()
{
  double pos0[3], pos1[3], mass0, mass1;
  int idx, i, k;

  mass0 = 0;
  mass1 = 0;
  for(k = 0; k < 3; k++)
    {
      pos0[k] = 0;
      pos1[k] = 0;
    }

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      mass0 += P[i].Mass * SphP[i].PScalars[0];
      mass1 += P[i].Mass * SphP[i].PScalars[1];
      for(k = 0; k < 3; k++)
        {
          pos0[k] += P[i].Mass * SphP[i].PScalars[0] * (SphP[i].Center[k] - 0.5 * All.BoxSize);
          pos1[k] += P[i].Mass * SphP[i].PScalars[1] * (SphP[i].Center[k] - 0.5 * All.BoxSize);
        }
    }

  double m0tot, m1tot;
  double p0tot[3], p1tot[3];

  MPI_Allreduce(&mass0, &m0tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&mass1, &m1tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  MPI_Allreduce(pos0, p0tot, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(pos1, p1tot, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  double a = fabs(p0tot[0] / m0tot - p1tot[0] / m1tot);
  All.Omega = sqrt(All.G * (m0tot + m1tot) / (a * a * a));
}

void apply_orbital_acceleration()
{
  int idx;
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      int i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      double dt_cell = 0.5 * (P[i].TimeBinGrav ? (((integertime) 1) << P[i].TimeBinGrav) : 0) * All.Timebase_interval;

      double *pos;
      if(P[i].Type == 0)
        pos = SphP[i].Center;
      else
        pos = P[i].Pos;

      double accel[3];
      accel[0] = All.Omega * All.Omega * (pos[0] - 0.5 * All.BoxSize) + 2. * All.Omega * P[i].Vel[1];
      accel[1] = All.Omega * All.Omega * (pos[1] - 0.5 * All.BoxSize) - 2. * All.Omega * P[i].Vel[0];
      accel[2] = 0;

      double dvel[3];
      int k;
      for(k = 0; k < 3; k++)
        dvel[k] = accel[k] * dt_cell;

      if(P[i].Type == 0)
        {
          for(k = 0; k < 3; k++)
            SphP[i].Momentum[k] += dvel[k] * P[i].Mass;
        }
      for(k = 0; k < 3; k++)
        P[i].Vel[k] += dvel[k];
    }
}

#endif

#ifdef RELAXOBJECT_COOLING2
void load_temperature_profil()
{
  FILE *fp;
  double rad, temp;
  int  len, j;

  fp = fopen("TemperatureProfil.txt","r");
  if (!fp)
    {
      printf("File could not be opened.\n");
      terminate("blurb");
    }
  fscanf(fp, "%x\n", &len);
//printf("length %d\n", length);
  wdR = (double*) malloc( len * sizeof(double) );
  wdTemp = (double*) malloc( len * sizeof(double) );

  for (j=0; j< len; j++)
    {
      fscanf(fp, "%lf %lf\n", &rad, &temp);
//    printf("%f %f\n", rad, temp);
      wdR[j] = rad;
      wdTemp[j] = temp;
    }
  fclose(fp);
  length = len;
}
#endif
