/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dust_live/dust_transfer_kernels.c
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
#ifdef DL_GRAIN_BINS

typedef struct
{
  MyDouble Pos[3];
  MyDouble Vel[3];
  MyFloat Hsml;
  MyFloat NormSph;
  MyFloat DeltaMassExpected;
  MyFloat MetalFractions[GFM_N_CHEM_ELEMENTS];
  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

typedef struct
{
  MyFloat DeltaMassActual;
  MyFloat DeltaMetalMasses[GFM_N_CHEM_ELEMENTS];
  MyFloat DeltaMomentum[3];
} data_out;

static data_out *DataResult, *DataOut;

static void particle2in(data_in *in, int i, int firstnode)
{
  in->Pos[0] = P[DustParticle[i].index].Pos[0];
  in->Pos[1] = P[DustParticle[i].index].Pos[1];
  in->Pos[2] = P[DustParticle[i].index].Pos[2];
  in->Vel[0] = P[DustParticle[i].index].Vel[0];
  in->Vel[1] = P[DustParticle[i].index].Vel[1];
  in->Vel[2] = P[DustParticle[i].index].Vel[2];
  in->Hsml = DTP(DustParticle[i].index).Hsml;
  in->NormSph = DustParticle[i].NormSph;
  in->DeltaMassExpected = DustParticle[i].DeltaMassExpected;
  for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
    {
      in->MetalFractions[k] = DTP(DustParticle[i].index).MetalFractions[k];
    }
  in->Firstnode = firstnode;
}

static void out2particle(data_out *out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES)
    {
      DustParticle[i].DeltaMassActual = out->DeltaMassActual;
      for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
        {
          DustParticle[i].DeltaMetalMasses[k] = out->DeltaMetalMasses[k];
        }
      for(int k = 0; k < 3; k++)
        {
          DustParticle[i].DeltaMomentum[k] = out->DeltaMomentum[k];
        }
    }
  else
    {
      DustParticle[i].DeltaMassActual += out->DeltaMassActual;
      for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
        {
          DustParticle[i].DeltaMetalMasses[k] += out->DeltaMetalMasses[k];
        }
      for(int k = 0; k < 3; k++)
        {
          DustParticle[i].DeltaMomentum[k] += out->DeltaMomentum[k];
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

        dust_transfer_kernel_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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

        dust_transfer_kernel_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

void dust_transfer_kernel(void)
{
  long long ntot;

  sumup_large_ints(1, &Ndust, &ntot);
  if(ntot == 0)
    return;

  mpi_printf("DUST_LIVE: Performing dust transfer update using local gas cells.\n");

  generic_set_MaxNexport();

  double t0 = second();

  generic_comm_pattern(Ndust, kernel_local, kernel_imported);

  double t1 = second();

  mpi_printf("DUST_LIVE: Done! Calculation took %g sec\n", timediff(t0, t1));
}

int dust_transfer_kernel_evaluate(int target, int mode, int thread_id)
{
  int numnodes, *firstnode;
  data_in local, *in;
  data_out out;

  double h, h2, h3;
  double dx, dy, dz, r2;

//#ifndef GFM_TOPHAT_KERNEL
//  double wk, u, hinv, hinv3;
//#endif
  double weight_fac, wk, hinv, hinv3;
#ifndef GFM_TOPHAT_KERNEL
  double u;
#endif

#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif

  out.DeltaMassActual = 0.0;
  for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
    {
      out.DeltaMetalMasses[k] = 0.0;
    }
  for(int k = 0; k < 3; k++)
    {
      out.DeltaMomentum[k] = 0.0;
    }

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

  MyDouble *pos = in->Pos;
  MyDouble *vel = in->Vel;
  h = in->Hsml;
  MyFloat normsph = in->NormSph;
  MyFloat delta_mass_expected = in->DeltaMassExpected;
  MyFloat *metal_fractions = in->MetalFractions;

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

              weight_fac = SphP[j].Volume * wk / normsph;
#else
              wk = hinv3 / NORM_COEFF;
              weight_fac = SphP[j].Volume / normsph;
#endif

              /* We use volume-weighted kernel interpolation over the gas
               * cells, because mass-weighted interpolation would run into
               * issues as mass is transferred to and from the gas. */
              if(delta_mass_expected > 0.0)
                {
                  double avail_metals[GFM_N_CHEM_ELEMENTS];
                  for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
                    {
                      if((elem_can_be_dust(k)) && (SphP[j].MassMetals[k] > 0.0))
                        avail_metals[k] = SphP[j].MassMetals[k];
                      else
                        avail_metals[k] = 0.0;
                    }

                  double avail_mass = 0.0;
                  for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
                    {
                      avail_mass += avail_metals[k];
                    }
                  double avail_fracs[GFM_N_CHEM_ELEMENTS];
                  for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
                    {
                      avail_fracs[k] = (avail_mass > 0.0 ? avail_metals[k]/avail_mass : 0.0);
                    }

                  double desired_mass = weight_fac * delta_mass_expected;
                  /* Prevent against any small negative values in SphP.MassMetals. */
                  double actual_mass = dmax(dmin(desired_mass, avail_mass), 0.0);
                  double dm = actual_mass;

                  double E_kin = 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] + SphP[j].Momentum[2] * SphP[j].Momentum[2]) / P[j].Mass;
                  double Utherm = (SphP[j].Energy - E_kin) / P[j].Mass;

                  P[j].Mass -= dm;
                  out.DeltaMassActual += dm;
                  SphP[j].MassMetallicity -= dm;
                  SphP[j].Metallicity = SphP[j].MassMetallicity / P[j].Mass; /* update primitives manually */
                  for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
                    {
                      double dm_metal = dm * avail_fracs[k];
                      SphP[j].MassMetals[k] -= dm_metal;
                      /* This check seems to really help prevent spurious
                       * negative metal masses. */
                      if(SphP[j].MassMetals[k] < 0.0)
                        {
                          SphP[j].MassMetals[k] = 0.0;
                        }
                      SphP[j].MetalsFraction[k] = SphP[j].MassMetals[k] / P[j].Mass; /* update primitives manually */
                      out.DeltaMetalMasses[k] += dm_metal;
                    }
                  for(int l = 0; l < 3; l++)
                    {
                      /* Assume velocity of gas is unchanged. */
                      SphP[j].Momentum[l] -= dm * P[j].Vel[l];
                      out.DeltaMomentum[l] += dm * P[j].Vel[l];
                    }
                  E_kin = 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] + SphP[j].Momentum[2] * SphP[j].Momentum[2]) / P[j].Mass;
                  /* New total energy keeps the old Utherm and adds in the new kinetic energy. */
                  SphP[j].Energy = Utherm * P[j].Mass + E_kin;
                }

              /* Dust particle is enriching surrounding gas. */
              if(delta_mass_expected < 0.0)
                {
                  double dm = weight_fac * fabs(delta_mass_expected);
                  double E_kin = 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] + SphP[j].Momentum[2] * SphP[j].Momentum[2]) / P[j].Mass;
                  double Utherm = (SphP[j].Energy - E_kin) / P[j].Mass;
                  P[j].Mass += dm;
                  out.DeltaMassActual -= dm;
                  SphP[j].MassMetallicity += dm;
                  SphP[j].Metallicity = SphP[j].MassMetallicity / P[j].Mass; /* update primitives manually */
                  for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
                    {
                      // Assume that a dust particle returns mass in different
                      // elements such that their relative ratios in dust don't
                      // change.
                      SphP[j].MassMetals[k] += dm * metal_fractions[k];
                      SphP[j].MetalsFraction[k] = SphP[j].MassMetals[k] / P[j].Mass; /* update primitives manually */
                      out.DeltaMetalMasses[k] -= dm * metal_fractions[k];
                    }
                  for(int l = 0; l < 3; l++)
                    {
                      double dmom = dm * vel[l];
                      SphP[j].Momentum[l] += dmom;
                      out.DeltaMomentum[l] -= dmom;
                    }
                  E_kin = 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] + SphP[j].Momentum[2] * SphP[j].Momentum[2]) / P[j].Mass;
                  /* New total energy keeps the old Utherm and adds in the new kinetic energy. */
                  SphP[j].Energy = Utherm * P[j].Mass + E_kin;
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
#endif
