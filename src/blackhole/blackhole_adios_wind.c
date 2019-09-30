/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/blackhole/blackhole.c
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

/*! \file blackhole_feedback.c
 *  \brief distribute feedback energy
 */


#ifdef BLACK_HOLES
#ifdef BH_ADIOS_WIND


void blackhole_update_wind_affected_cells(void)
{
  int idx, i, j;

  mpi_printf("BH_ADIOS_WIND: Updating momenta of affected cells.\n");

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Mass == 0 && P[i].ID == 0)
        continue;

      if(P[i].Type == 0)
        {
          SphP[i].Energy -= 0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);

          for(j = 0; j < 3; j++)        /* do the kick for gas cells */
            {
              P[i].Vel[j] += SphP[i].Injected_BH_Wind_Momentum[j] / P[i].Mass;
              SphP[i].Momentum[j] += SphP[i].Injected_BH_Wind_Momentum[j];

              SphP[i].Injected_BH_Wind_Momentum[j] = 0;
            }

          SphP[i].Energy += 0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);
        }
    }
}



#ifndef BH_ADIOS_RANDOMIZED


static int blackhole_evaluate(int target, int mode, int threadid);
static int phase;


/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble Pos[3];
  MyFloat BH_Hsml;

  MyFloat BH_WindEnergy;

  MyFloat Asum;
  MyFloat Bsum;

  MyFloat Msum;
  MyFloat Qsum[3];
#ifdef BH_ADIOS_WIND_DIRECTIONAL
  MyFloat WindDir[3];
#endif

  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

/* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in * in, int i, int firstnode)
{
  for(int k = 0; k < 3; k++)
    in->Pos[k] = P[i].Pos[k];

  in->BH_Hsml = BPP(i).BH_Hsml;

  if(phase > 0)
    {
      in->Msum = BPP(i).Msum;

      for(int k = 0; k < 3; k++)
        in->Qsum[k] = BPP(i).Qsum[k];
    }

  if(phase == 2)
    {
      in->Asum = BPP(i).Asum;
      in->Bsum = BPP(i).Bsum;

      in->BH_WindEnergy = BPP(i).BH_WindEnergy;
    }

#ifdef BH_ADIOS_WIND_DIRECTIONAL
  for(int k = 0; k < 3; k++)
    in->WindDir[k] = BPP(i).WindDir[k];
#endif

  in->Firstnode = firstnode;
}




/* local data structure that holds results acquired on remote processors */
typedef struct
{
  MyFloat Asum;
  MyFloat Bsum;

  MyFloat Msum;
  MyFloat Qsum[3];
} data_out;

static data_out *DataResult, *DataOut;



/* routine to store or combine result data */
static void out2particle(data_out * out, int i, int mode)
{
  int j;

  if(phase == 0)
    {
      if(mode == MODE_LOCAL_PARTICLES)  /* initial store */
        {
          BPP(i).Msum = out->Msum;
          for(j = 0; j < 3; j++)
            BPP(i).Qsum[j] = out->Qsum[j];
        }
      else
        {
          BPP(i).Msum += out->Msum;

          for(j = 0; j < 3; j++)
            BPP(i).Qsum[j] += out->Qsum[j];
        }
    }
  else if(phase == 1)
    {
      if(mode == MODE_LOCAL_PARTICLES)  /* initial store */
        {
          BPP(i).Asum = out->Asum;
          BPP(i).Bsum = out->Bsum;
        }
      else
        {
          BPP(i).Asum += out->Asum;
          BPP(i).Bsum += out->Bsum;
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
#pragma omp parallel private(i)
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
              if(generic_polling_primary(count, TimeBinsBHAccretion.NActiveParticles))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        i = NextParticle++;

        if(i >= TimeBinsBHAccretion.NActiveParticles)
          break;

        int idx = TimeBinsBHAccretion.ActiveParticleList[i];

        if(idx < 0)
          continue;

        if(BPP(idx).BH_WindEnergy == 0)
          continue;

        blackhole_evaluate(idx, MODE_LOCAL_PARTICLES, threadid);
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

        blackhole_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

/* This function distributes the feedback energy and look for gas cells that we might stochastically swallow
 */
void blackhole_blow_wind(void)
{
  mpi_printf("BH_ADIOS_WIND: Start assigning BH momentum feedback\n");
  double t0 = second();

#ifdef BH_ADIOS_WIND_DIRECTIONAL
  for(int i = 0; i < TimeBinsBHAccretion.NActiveParticles; i++)
     {
       int idx = TimeBinsBHAccretion.ActiveParticleList[i];
       if(idx < 0)
         continue;

       double theta = acos(2 * get_random_number() - 1);
       double phi = 2 * M_PI * get_random_number();
       BPP(idx).WindDir[0] = sin(theta) * cos(phi);
       BPP(idx).WindDir[1] = sin(theta) * sin(phi);
       BPP(idx).WindDir[2] = cos(theta);
     }
#endif

  generic_set_MaxNexport();

  for(phase = 0; phase < 3; phase++)
    {
      generic_comm_pattern(TimeBinsBHAccretion.NActiveParticles, kernel_local, kernel_imported);
    }

  for(int i = 0; i < TimeBinsBHAccretion.NActiveParticles; i++)
    {
      int idx = TimeBinsBHAccretion.ActiveParticleList[i];
      if(idx < 0)
        continue;

      BPP(idx).BH_WindEnergy = 0;
    }

  double t1 = second();
  mpi_printf("BH_ADIOS_WIND: Done assigning BH momentum feedback\n", timediff(t0, t1));
}




static int blackhole_evaluate(int target, int mode, int threadid)
{
  int numnodes, *firstnode, j, n;
  MyDouble *pos;
  double dx, dy, dz, h_i, h_i2, r2, r, u, hinv, hinv3, wk, fac = 0.0, vel[3];
#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif

  data_in local, *in;
  data_out out;;

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
  h_i = in->BH_Hsml;

  out.Msum = 0;
  out.Asum = 0;
  out.Bsum = 0;
  for(j = 0; j < 3; j++)
    out.Qsum[j] = 0;

  if(phase == 2)
    {
      double windenergy = in->BH_WindEnergy;

      if(in->Asum > 0)
        fac = (-in->Bsum + sqrt(in->Bsum * in->Bsum + 4 * in->Asum * windenergy * All.cf_atime * All.cf_atime)) / (2 * in->Asum);
    }

  h_i2 = h_i * h_i;

  int nfound = ngb_treefind_variable_threads(pos, h_i, target, mode, threadid, numnodes, firstnode);

  for(n = 0; n < nfound; n++)
    {
      j = Thread[threadid].Ngblist[n];

      if(P[j].Mass > 0)
        {
          dx = NEAREST_X(P[j].Pos[0] - pos[0]);
          dy = NEAREST_Y(P[j].Pos[1] - pos[1]);
          dz = NEAREST_Z(P[j].Pos[2] - pos[2]);

          r2 = dx * dx + dy * dy + dz * dz;

          if(r2 < h_i2 && P[j].Mass > 0 && P[j].ID != 0)
            {
              if(P[j].Type == 0)
                {
                  r = sqrt(r2);
                  hinv = 1 / h_i;
                  hinv3 = hinv * hinv * hinv;

                  u = r * hinv;

                  if(u < 0.5)
                    wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
                  else
                    wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

                  if(r > 0)
                    {
                      dx /= r;
                      dy /= r;
                      dz /= r;

#ifdef BH_ADIOS_WIND_DIRECTIONAL
                      double dircos = cos(dx * in->WindDir[0] + dy * in->WindDir[1] + dz * in->WindDir[2]);
                      wk *= dircos * dircos;
#endif

                      if(phase > 0)
                        {
                          dx -= in->Qsum[0] / in->Msum;
                          dy -= in->Qsum[1] / in->Msum;
                          dz -= in->Qsum[2] / in->Msum;
                        }

                      if(phase == 0)
                        {
                          out.Msum += P[j].Mass * wk;
                          out.Qsum[0] += P[j].Mass * wk * dx;
                          out.Qsum[1] += P[j].Mass * wk * dy;
                          out.Qsum[2] += P[j].Mass * wk * dz;
                        }
                      else if(phase == 1)
                        {
                    	  vel[0] = (P[j].Vel[0] + SphP[j].Injected_BH_Wind_Momentum[0] / P[j].Mass);
                    	  vel[1] = (P[j].Vel[1] + SphP[j].Injected_BH_Wind_Momentum[1] / P[j].Mass);
                    	  vel[2] = (P[j].Vel[2] + SphP[j].Injected_BH_Wind_Momentum[2] / P[j].Mass);

                          out.Asum += 0.5 * P[j].Mass * wk * wk * (dx * dx + dy * dy + dz * dz);
                          out.Bsum += P[j].Mass * wk * (vel[0] * dx + vel[1] * dy + vel[2] * dz);
                        }
                      else if(phase == 2)
                        {
                          SphP[j].Injected_BH_Wind_Momentum[0] += P[j].Mass * wk * fac * dx;
                          SphP[j].Injected_BH_Wind_Momentum[1] += P[j].Mass * wk * fac * dy;
                          SphP[j].Injected_BH_Wind_Momentum[2] += P[j].Mass * wk * fac * dz;
                        }
                    }
                }
            }
        }
    }

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#endif
#endif
#endif
