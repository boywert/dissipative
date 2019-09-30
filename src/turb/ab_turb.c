/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/turb/ab_turb.c
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

#ifdef AB_TURB


double st_grn(void);
void st_init_ouseq(void);
void st_calc_phases(void);
void st_compute_injection(void);


//Ornstein-Uhlenbeck variables
double StOUVar;
double *StOUPhases;
gsl_rng *StRng;


//forcing field in fourie space
double *StAmpl;
double *StAka;                  //phases (real part)
double *StAkb;                  //phases (imag part)
double *StMode;
int StNModes;


integertime StTPrev;
double StSolWeightNorm;


void init_turb(void)
{
  int ikx, iky, ikz;
  double kx, ky, kz, k;
  double ampl;

  int ikxmax = All.BoxSize * All.StKmax / 2. / M_PI;
#ifndef ONEDIMS
  int ikymax = All.BoxSize * All.StKmax / 2. / M_PI;

#ifndef TWODIMS
  int ikzmax = All.BoxSize * All.StKmax / 2. / M_PI;
#else
  int ikzmax = 0;
#endif

#else
  int ikymax = 0;
  int ikzmax = 0;
#endif

  StNModes = 0;
  for(ikx = 0; ikx <= ikxmax; ikx++)
    {
      kx = 2. * M_PI * ikx / All.BoxSize;
      for(iky = 0; iky <= ikymax; iky++)
        {
          ky = 2. * M_PI * iky / All.BoxSize;
          for(ikz = 0; ikz <= ikzmax; ikz++)
            {
              kz = 2. * M_PI * ikz / All.BoxSize;
              k = sqrt(kx * kx + ky * ky + kz * kz);
              if(k >= All.StKmin && k <= All.StKmax)
                {
#if NUMDIMS ==1
                  StNModes += 1;
#endif
#if NUMDIMS == 2
                  StNModes += 2;
#endif
#if NUMDIMS == 3
                  StNModes += 4;
#endif
                }
            }
        }
    }

  if(ThisTask == 0)
    mpi_printf("TURB: Using %d modes, %d %d %d\n", StNModes, ikxmax, ikymax, ikzmax);

  StMode = (double *) mymalloc_movable(&StMode, "StModes", StNModes * 3 * sizeof(double));
  StAka = (double *) mymalloc_movable(&StAka, "StAka", StNModes * 3 * sizeof(double));
  StAkb = (double *) mymalloc_movable(&StAkb, "StAkb", StNModes * 3 * sizeof(double));
  StAmpl = (double *) mymalloc_movable(&StAmpl, "StAmpl", StNModes * sizeof(double));
  StOUPhases = (double *) mymalloc_movable(StOUPhases, "StOUPhases", StNModes * 6 * sizeof(double));

  StOUVar = sqrt(All.StEnergy / All.StDecay);
  double kc = 0.5 * (All.StKmin + All.StKmax);
  double amin = 0.;

#if NUMDIMS == 3
  StSolWeightNorm = sqrt(3.0 / 3.0) * sqrt(3.0) * 1.0 / sqrt(1.0 - 2.0 * All.StSolWeight + 3.0 * All.StSolWeight * All.StSolWeight);
#endif
#if NUMDIMS == 2
  StSolWeightNorm = sqrt(3.0 / 2.0) * sqrt(3.0) * 1.0 / sqrt(1.0 - 2.0 * All.StSolWeight + 2.0 * All.StSolWeight * All.StSolWeight);
#endif
#if NUMDIMS == 1
  StSolWeightNorm = sqrt(3.0 / 1.0) * sqrt(3.0) * 1.0 / sqrt(1.0 - 2.0 * All.StSolWeight + 1.0 * All.StSolWeight * All.StSolWeight);
#endif

  StNModes = 0;

  for(ikx = 0; ikx <= ikxmax; ikx++)
    {
      kx = 2. * M_PI * ikx / All.BoxSize;

      for(iky = 0; iky <= ikymax; iky++)
        {
          ky = 2. * M_PI * iky / All.BoxSize;

          for(ikz = 0; ikz <= ikzmax; ikz++)
            {
              kz = 2. * M_PI * ikz / All.BoxSize;

              k = sqrt(kx * kx + ky * ky + kz * kz);
              if(k >= All.StKmin && k <= All.StKmax)
                {
                  if(All.StSpectForm == 0)
                    {
                      ampl = 1.;
                    }
                  else if(All.StSpectForm == 1)
                    {
                      ampl = 4.0 * (amin - 1.0) / ((All.StKmax - All.StKmin) * (All.StKmax - All.StKmin)) * ((k - kc) * (k - kc)) + 1.0;
                    }
                  else if(All.StSpectForm == 2)
                    {
                      ampl = pow(All.StKmin, 5. / 3) / pow(k, 5. / 3);
                    }
                  else if(All.StSpectForm == 3)
                    {
                      ampl = pow(All.StKmin, 2.) / pow(k, 2.);
                    }
                  else
                    {
                      terminate("unkown spectral form");
                    }

                  StAmpl[StNModes] = ampl;
                  StMode[3 * StNModes + 0] = kx;
                  StMode[3 * StNModes + 1] = ky;
                  StMode[3 * StNModes + 2] = kz;
                  if(ThisTask == 0)
                    mpi_printf("TURB: Mode: %d, ikx=%d, iky=%d, ikz=%d, kx=%f, ky=%f, kz=%f, ampl=%f\n", StNModes, ikx, iky, ikz, kx, ky, kz, ampl);

                  StNModes++;

#if NUMDIMS > 1
                  StAmpl[StNModes] = ampl;
                  StMode[3 * StNModes + 0] = kx;
                  StMode[3 * StNModes + 1] = -ky;
                  StMode[3 * StNModes + 2] = kz;
                  if(ThisTask == 0)
                    mpi_printf("TURB: Mode: %d, ikx=%d, iky=%d, ikz=%d, kx=%f, ky=%f, kz=%f, ampl=%f\n", StNModes, ikx, -iky, ikz, kx, -ky, kz, ampl);

                  StNModes++;

#if NUMDIMS > 2
                  StAmpl[StNModes] = ampl;
                  StMode[3 * StNModes + 0] = kx;
                  StMode[3 * StNModes + 1] = ky;
                  StMode[3 * StNModes + 2] = -kz;
                  if(ThisTask == 0)
                    mpi_printf("TURB: Mode: %d, ikx=%d, iky=%d, ikz=%d, kx=%f, ky=%f, kz=%f, ampl=%f\n", StNModes, ikx, iky, -ikz, kx, ky, -kz, ampl);

                  StNModes++;

                  StAmpl[StNModes] = ampl;
                  StMode[3 * StNModes + 0] = kx;
                  StMode[3 * StNModes + 1] = -ky;
                  StMode[3 * StNModes + 2] = -kz;
                  if(ThisTask == 0)
                    mpi_printf("TURB: Mode: %d, ikx=%d, iky=%d, ikz=%d, kx=%f, ky=%f, kz=%f, ampl=%f\n", StNModes, ikx, -iky, -ikz, kx, -ky, -kz, ampl);

                  StNModes++;
#endif
#endif
                }
            }
        }
    }

  StTPrev = -1;


  StRng = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(StRng, All.StSeed);

  st_init_ouseq();
  st_calc_phases();

  mpi_printf("TURB: calling set_turb_ampl in init_turb\n");
  set_turb_ampl();

  StTPrev = All.Ti_Current;

  init_static_nfw();
  mpi_printf("TURB: init done.\n");
}

void st_init_ouseq(void)
{
  int i;

  for(i = 0; i < 6 * StNModes; i++)
    {
      StOUPhases[i] = st_grn() * StOUVar;
    }
}

void st_update_ouseq(void)
{
  int i;
  double damping = exp(-All.StDtFreq / All.StDecay);

  for(i = 0; i < 6 * StNModes; i++)
    {
      StOUPhases[i] = StOUPhases[i] * damping + StOUVar * sqrt(1. - damping * damping) * st_grn();
    }
}

double st_grn(void)
{
  double r0 = gsl_rng_uniform(StRng);
  double r1 = gsl_rng_uniform(StRng);

  return sqrt(2. * log(1. / r0)) * cos(2. * M_PI * r1);
}



void st_calc_phases(void)
{
  int i, j;
  for(i = 0; i < StNModes; i++)
    {
      double ka = 0.;
      double kb = 0.;
      double kk = 0.;

      int dim = NUMDIMS;

      for(j = 0; j < dim; j++)
        {
          kk += StMode[3 * i + j] * StMode[3 * i + j];
          ka += StMode[3 * i + j] * StOUPhases[6 * i + 2 * j + 1];
          kb += StMode[3 * i + j] * StOUPhases[6 * i + 2 * j + 0];
        }
      for(j = 0; j < dim; j++)
        {
          double diva = StMode[3 * i + j] * ka / kk;
          double divb = StMode[3 * i + j] * kb / kk;
          double curla = StOUPhases[6 * i + 2 * j + 0] - divb;
          double curlb = StOUPhases[6 * i + 2 * j + 1] - diva;

          StAka[3 * i + j] = All.StSolWeight * curla + (1. - All.StSolWeight) * divb;
          StAkb[3 * i + j] = All.StSolWeight * curlb + (1. - All.StSolWeight) * diva;
        }
    }
}

void set_turb_ampl(void)
{
  TIMER_START(CPU_TURB_UPDATE) int i;

  double delta = (All.Ti_Current - StTPrev) * All.Timebase_interval;
  if(delta >= All.StDtFreq)
    {
      if(delta > 0)
        {
          mpi_printf("TURB: updating dudt_*\n");
          for(i = 0; i < NumGas; i++)
            {
              SphP[i].DuDt_diss = SphP[i].EgyDiss / P[i].Mass / delta;
              SphP[i].EgyDiss = 0;

              SphP[i].DuDt_drive = SphP[i].EgyDrive / P[i].Mass / delta;
              SphP[i].EgyDrive = 0;

#ifdef DG
              int j;

              for(j = 0; j < NOF_BASE_FUNCTIONS; j++)
                {
                  SphP[i].WeightsDeDt[j] = SphP[i].WeightsDiss[j] / delta;
                  SphP[i].WeightsDiss[j] = 0;
                }

#endif
            }
        }
      st_update_ouseq();
      st_calc_phases();

      StTPrev = StTPrev + All.StDtFreq / All.Timebase_interval;

      mpi_printf("TURB: Updated turbulent stirring field at time %f.\n", StTPrev * All.Timebase_interval);
    }

TIMER_STOP(CPU_TURB_UPDATE)}

#ifndef DG
void add_turb_accel()
{
  int idx, i, j, m;
  double acc[3];

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      double fx = 0;
      double fy = 0;
      double fz = 0;

      for(m = 0; m < StNModes; m++)     //calc force
        {
          double kxx = StMode[3 * m + 0] * P[i].Pos[0];
          double kyy = StMode[3 * m + 1] * P[i].Pos[1];
          double kzz = StMode[3 * m + 2] * P[i].Pos[2];
          double kdotx = kxx + kyy + kzz;
          double ampl = StAmpl[m];

          double realt = cos(kdotx);
          double imagt = sin(kdotx);

          fx += ampl * (StAka[3 * m + 0] * realt - StAkb[3 * m + 0] * imagt);
          fy += ampl * (StAka[3 * m + 1] * realt - StAkb[3 * m + 1] * imagt);
          fz += ampl * (StAka[3 * m + 2] * realt - StAkb[3 * m + 2] * imagt);
        }

      fx *= 2. * All.StAmplFac * StSolWeightNorm;
      fy *= 2. * All.StAmplFac * StSolWeightNorm;
      fz *= 2. * All.StAmplFac * StSolWeightNorm;

      if(P[i].Mass > 0.)
        {
          acc[0] = fx;
          acc[1] = fy;
#ifndef TWODIMS
          acc[2] = fz;
#else
          acc[2] = 0;
#endif

          for(j = 0; j < 3; j++)
            SphP[i].TurbAccel[j] = acc[j];
        }
    }

  mpi_printf("TURB: Finished turbulent accel computation.\n");
}

#else
void dg_acceleration(double x, double y, double z, double *acc)
{

  TIMER_START(CPU_TURB_FORCE) double fx = 0.;
  double fy = 0.;
#ifndef TWODIMS
  double fz = 0.;
#endif

  int m;
  for(m = 0; m < StNModes; m++) //calc force
    {
      double kxx = StMode[3 * m + 0] * x;
      double kyy = StMode[3 * m + 1] * y;
#ifndef TWODIMS
      double kzz = StMode[3 * m + 2] * z;
#else
      double kzz = 0.;
#endif

      double kdotx = kxx + kyy + kzz;
      double ampl = StAmpl[m];

      double realt = cos(kdotx);
      double imagt = sin(kdotx);

      fx += ampl * (StAka[3 * m + 0] * realt - StAkb[3 * m + 0] * imagt);
      fy += ampl * (StAka[3 * m + 1] * realt - StAkb[3 * m + 1] * imagt);
#ifndef TWODIMS
      fz += ampl * (StAka[3 * m + 2] * realt - StAkb[3 * m + 2] * imagt);
#endif
    }

  fx *= 2. * All.StAmplFac * StSolWeightNorm;
  fy *= 2. * All.StAmplFac * StSolWeightNorm;
#ifndef TWODIMS
  fz *= 2. * All.StAmplFac * StSolWeightNorm;
#endif

  acc[0] = fx;
  acc[1] = fy;
#ifndef TWODIMS
  acc[2] = fz;
#endif
TIMER_STOP(CPU_TURB_FORCE)}
#endif

#endif
