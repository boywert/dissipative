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


#ifdef BH_ADIOS_RANDOMIZED


static int blackhole_evaluate(int target, int mode, int threadid);



/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble Pos[3];
  MyFloat BH_Hsml;
  MyFloat BH_WindEnergy;
  MyFloat BH_Density;
  MyFloat WindDir[3];

  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

/* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in * in, int i, int firstnode)
{
  for(int k = 0; k < 3; k++)
    in->Pos[k] = P[i].Pos[k];

  in->BH_Hsml = BPP(i).BH_Hsml;
  in->BH_WindEnergy = BPP(i).BH_WindEnergy;
  in->BH_Density = BPP(i).BH_Density;

  for(int k = 0; k < 3; k++)
    in->WindDir[k] = BPP(i).WindDir[k];

  in->Firstnode = firstnode;
}




/* local data structure that holds results acquired on remote processors */
typedef struct
{
   int dummy;
} data_out;

static data_out *DataResult, *DataOut;



/* routine to store or combine result data */
static void out2particle(data_out * out, int i, int mode)
{

}


#include "../generic_comm_helpers2.h"


static int blackhole_wind_test_for_injection(int idx)
{
#if (defined(GFM_WINDS_VARIABLE) && (GFM_WINDS_VARIABLE==1)) || defined(GFM_WINDS_LOCAL)
  double egycmp = All.RadioFeedbackReiorientationFactor * 0.5 * pow(BPP(idx).BH_DMVelDisp, 2) * BPP(idx).BH_Density * (4.0 * M_PI/3.0) * pow(BPP(idx).BH_Hsml, 3);
#else
  double egycmp = All.RadioFeedbackReiorientationFactor * BPP(idx).BH_U * BPP(idx).BH_Density * (4.0 * M_PI/3.0) * pow(BPP(idx).BH_Hsml, 3);
#endif

  if(BPP(idx).BH_WindEnergy >= egycmp)
    return 1;
  else
    return 0;
}


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

        if(blackhole_wind_test_for_injection(idx))
          {
            /* do this not necessarily every BH timestep */
            double theta = acos(2 * get_random_number() - 1);
            double phi = 2 * M_PI * get_random_number();
            BPP(idx).WindDir[0] = sin(theta) * cos(phi);
            BPP(idx).WindDir[1] = sin(theta) * sin(phi);
            BPP(idx).WindDir[2] = cos(theta);

            blackhole_evaluate(idx, MODE_LOCAL_PARTICLES, threadid);
          }
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

  generic_set_MaxNexport();

  generic_comm_pattern(TimeBinsBHAccretion.NActiveParticles, kernel_local, kernel_imported);


  for(int i = 0; i < TimeBinsBHAccretion.NActiveParticles; i++)
    {
      int idx = TimeBinsBHAccretion.ActiveParticleList[i];
      if(idx < 0)
        continue;

      if(blackhole_wind_test_for_injection(idx))
        BPP(idx).BH_WindEnergy = 0;
    }

  double t1 = second();
  mpi_printf("BH_ADIOS_WIND: Done assigning BH momentum feedback\n", timediff(t0, t1));
}




static int blackhole_evaluate(int target, int mode, int threadid)
{
  int numnodes, *firstnode;

  data_in local, *in;
  data_out out;

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
  double h_i = in->BH_Hsml;
  double windenergy = in->BH_WindEnergy;
  double dens = in->BH_Density;

  double hinv = 1.0 / h_i;
  double hinv3 = hinv * hinv * hinv;

  int nfound = ngb_treefind_variable_threads(pos, h_i, target, mode, threadid, numnodes, firstnode);

  for(int n = 0; n < nfound; n++)
    {
      int j = Thread[threadid].Ngblist[n];

      if(P[j].Mass > 0 && P[j].ID != 0)
        {
          double r2 = Thread[threadid].R2list[n];
          double r = sqrt(r2);
          double u = r * hinv;

          double wk;

          if(u < 0.5)
            wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
          else
            wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

          double fac = All.cf_atime * P[j].Mass * sqrt(2.0 * windenergy * wk / dens);

          SphP[j].Injected_BH_Wind_Momentum[0] += fac * in->WindDir[0];
          SphP[j].Injected_BH_Wind_Momentum[1] += fac * in->WindDir[1];
          SphP[j].Injected_BH_Wind_Momentum[2] += fac * in->WindDir[2];
        }
    }

  out.dummy = 0;

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
