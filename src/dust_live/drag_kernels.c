/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dust_live/drag_kernels.c
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
#include <gsl/gsl_math.h>

#include "../allvars.h"
#include "../proto.h"

#ifdef DUST_LIVE

typedef struct
{
  MyDouble Pos[3];
  MyFloat Hsml;
  MyFloat TotNgbMass;
  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

typedef struct
{
  MyFloat LocalGasDensity;
  MyFloat LocalGasVelocity[3];
  MyFloat LocalSoundSpeed;
  MyFloat LocalGasAccel[3];
  MyFloat LocalGradP[3];
#ifdef DL_GRAIN_BINS
  MyFloat LocalGasDensityH;
  MyFloat LocalGasZ;
  MyFloat LocalGasTemp;
#ifdef DL_SNE_DESTRUCTION
  MyFloat LocalSNPrefactor;
#endif
#if defined(DL_SNE_DESTRUCTION) || defined(DL_SHATTERING) || defined(DL_COAGULATION)
  MyFloat LocalCloudFrac;
#endif
#endif
} data_out;

static data_out *DataResult, *DataOut;

static void particle2in(data_in *in, int i, int firstnode)
{
  in->Pos[0] = P[DustParticle[i].index].Pos[0];
  in->Pos[1] = P[DustParticle[i].index].Pos[1];
  in->Pos[2] = P[DustParticle[i].index].Pos[2];
  in->Hsml = DTP(DustParticle[i].index).Hsml;
  in->TotNgbMass = DustParticle[i].TotNgbMass;
  in->Firstnode = firstnode;
}

static void out2particle(data_out *out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES)
    {
      DTP(DustParticle[i].index).LocalGasDensity = out->LocalGasDensity;
      DTP(DustParticle[i].index).LocalSoundSpeed = out->LocalSoundSpeed;
#ifdef DL_GRAIN_BINS
      DustParticle[i].LocalGasDensityH = out->LocalGasDensityH;
      DustParticle[i].LocalGasZ = out->LocalGasZ;
      DustParticle[i].LocalGasTemp = out->LocalGasTemp;
#ifdef DL_SNE_DESTRUCTION
      DustParticle[i].LocalSNPrefactor = out->LocalSNPrefactor;
#endif
#if defined(DL_SNE_DESTRUCTION) || defined(DL_SHATTERING) || defined(DL_COAGULATION)
      DTP(DustParticle[i].index).LocalCloudFrac = out->LocalCloudFrac;
#endif
#endif
      for(int k = 0; k < 3; k++)
        {
          DTP(DustParticle[i].index).LocalGasVelocity[k] = out->LocalGasVelocity[k];
          DustParticle[i].LocalGasAccel[k] = out->LocalGasAccel[k];
          DustParticle[i].LocalGradP[k] = out->LocalGradP[k];
        }
    }
  else
    {
      DTP(DustParticle[i].index).LocalGasDensity += out->LocalGasDensity;
      DTP(DustParticle[i].index).LocalSoundSpeed += out->LocalSoundSpeed;
#ifdef DL_GRAIN_BINS
      DustParticle[i].LocalGasDensityH += out->LocalGasDensityH;
      DustParticle[i].LocalGasZ += out->LocalGasZ;
      DustParticle[i].LocalGasTemp += out->LocalGasTemp;
#ifdef DL_SNE_DESTRUCTION
      DustParticle[i].LocalSNPrefactor += out->LocalSNPrefactor;
#endif
#if defined(DL_SNE_DESTRUCTION) || defined(DL_SHATTERING) || defined(DL_COAGULATION)
      DTP(DustParticle[i].index).LocalCloudFrac += out->LocalCloudFrac;
#endif
#endif
      for(int k = 0; k < 3; k++)
        {
          DTP(DustParticle[i].index).LocalGasVelocity[k] += out->LocalGasVelocity[k];
          DustParticle[i].LocalGasAccel[k] += out->LocalGasAccel[k];
          DustParticle[i].LocalGradP[k] += out->LocalGradP[k];
        }
    }
}

#include "../generic_comm_helpers2.h"

static void kernel_local(void)
{
  int i;
#ifdef GENERIC_ASYNC
  int flag = 0;
#endif

#pragma omp parallel private(idx)
  {
    int j, threadid = get_thread_num();
#ifdef GENERIC_ASYNC
    int count = 0;
#endif

    for(j = 0; j < NTask; j++)
      Thread[threadid].Exportflag[j] = -1;

    while(1)
      {
        if(Thread[threadid].ExportSpace < MinSpace)
          break;

#ifdef GENERIC_ASYNC
        if(threadid == 0)
          {
            if((count & POLLINGINTERVAL) == 0)
              if(generic_polling_primary(count, Ndust))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        i = NextParticle++;

        if(i >= Ndust)
          break;

        drag_kernel_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
      }
  }
}

static void kernel_imported(void)
{
  /* now do the particles that were sent to us */
  int i, cnt = 0;
#pragma omp parallel private(i)
  {
    int threadid = get_thread_num();
#ifdef GENERIC_ASYNC
    int count = 0;
#endif

    while(1)
      {
#pragma omp atomic capture
        i = cnt++;

        if(i >= Nimport)
          break;

#ifdef GENERIC_ASYNC
        if(threadid == 0)
          {
            if((count & POLLINGINTERVAL) == 0)
              generic_polling_secondary();
          }

        count++;
#endif

        drag_kernel_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

void drag_kernel(void)
{
  long long ntot;

  sumup_large_ints(1, &Ndust, &ntot);
  if(ntot == 0)
    return;

  mpi_printf("DUST_LIVE: Using gas cells for density estimates.\n");

  generic_set_MaxNexport();

  double t0 = second();

  generic_comm_pattern(Ndust, kernel_local, kernel_imported);

  double t1 = second();

  mpi_printf("DUST_LIVE: Done! Calculation took %g sec\n", timediff(t0, t1));
}

int drag_kernel_evaluate(int target, int mode, int thread_id)
{
  int numnodes, *firstnode;
  data_in local, *in;
  data_out out;

  double h, h2, h3;
  double dx, dy, dz, r2;
  int k;

  double weight_fac, wk, hinv, hinv3;
#ifndef GFM_TOPHAT_KERNEL
  double u;
#endif

#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif

  MyDouble *pos;
  memset(&out, 0, sizeof(data_out));

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      in = &local;

      numnodes = 1;
      firstnode = NULL;
    }
  else
    {
      in = &DataGet[target];
      generic_get_numnodes(target, &numnodes, &firstnode);
    }

  pos = in->Pos;
  h = in->Hsml;
  MyFloat totngbmass = in->TotNgbMass;

  h2 = h * h;
  h3 = h * h * h;

  hinv = 1.0 / h;
#ifndef  TWODIMS
  hinv3 = hinv * hinv * hinv;
#else
  hinv3 = hinv * hinv / boxSize_Z;
#endif

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, thread_id, numnodes, firstnode);

  for(int n = 0; n < nfound; n++)
    {
      int j = Thread[thread_id].Ngblist[n];

      if(P[j].Mass > 0 && P[j].ID != 0) /* skip cells that have been swallowed or dissolved */
        {
          dx = NGB_PERIODIC_LONG_X(pos[0] - P[j].Pos[0]);
          dy = NGB_PERIODIC_LONG_Y(pos[1] - P[j].Pos[1]);
          dz = NGB_PERIODIC_LONG_Z(pos[2] - P[j].Pos[2]);

          r2 = dx * dx + dy * dy + dz * dz;

          if(r2 < h2)
            {
#ifndef GFM_TOPHAT_KERNEL
              u = sqrt(r2) * hinv;
              if(u < 0.5)
                wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
              else
                wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

              weight_fac = NORM_COEFF * P[j].Mass * wk * h3 / totngbmass;

#else
              wk = hinv3 / NORM_COEFF;
              weight_fac = P[j].Mass / totngbmass;
#endif

              out.LocalGasDensity += wk * P[j].Mass;
              out.LocalSoundSpeed += weight_fac * get_sound_speed(j);
#ifdef DL_GRAIN_BINS
              out.LocalGasDensityH += wk * SphP[j].MassMetals[element_index_Hydrogen];
              out.LocalGasZ += weight_fac * SphP[j].MassMetallicity / P[j].Mass;
              // TODO: do this properly, check for config conflicts
              double yHe = (1.0 - HYDROGEN_MASSFRAC) / (4.0 * HYDROGEN_MASSFRAC);
#ifdef COOLING
              double mu = (1.0 + 4.0*yHe) / (1.0 + yHe + SphP[j].Ne);
#else
              double mu = (1.0 + 4.0*yHe) / (1.0 + yHe);
#endif
              double SphP_temp = (mu * PROTONMASS) / BOLTZMANN * GAMMA_MINUS1 * SphP[j].Utherm * All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;
#ifdef USE_SFR
              if (SphP[j].Sfr > 0.0)
                SphP_temp = 1.0e4;
#endif
              out.LocalGasTemp += weight_fac * SphP_temp;

#ifdef DL_SNE_DESTRUCTION
              double n_cgs = SphP[j].Density * All.cf_a3inv * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam / PROTONMASS; /* cm^-3 */
              /* Just to avoid potentially small, negative values in MassMetallicity. */
              double Z_mod = dmax(SphP[j].MassMetallicity/P[j].Mass / 0.0127, 0.0) + 0.039;
              double M_swept = 1535.0 * pow(n_cgs, -0.202) * pow(Z_mod, -0.289); /* M_sun */
              M_swept *= (SOLAR_MASS / All.UnitMass_in_g * All.HubbleParam); /* internal mass units */
              double sn_prefactor = M_swept * SphP[j].SNRate / P[j].Mass; /* want to be 1/Gyr */
              out.LocalSNPrefactor += weight_fac * sn_prefactor;
#endif

#if defined(DL_SNE_DESTRUCTION) || defined(DL_SHATTERING) || defined(DL_COAGULATION)
              out.LocalCloudFrac += weight_fac * SphP[j].CloudFrac;
#endif
#endif
              for(k = 0; k < 3; k++)
                {
                  out.LocalGasVelocity[k] += weight_fac * P[j].Vel[k];
                  out.LocalGasAccel[k] += weight_fac * P[j].GravAccel[k];
#ifdef PMGRID
                  out.LocalGasAccel[k] += weight_fac * P[j].GravPM[k];
#endif
                  out.LocalGradP[k] += weight_fac * SphP[j].Grad.dpress[k];
                }
            }
        }
    }

  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#endif
