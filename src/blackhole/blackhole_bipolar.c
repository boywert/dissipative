/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/blackhole/blackhole_bipolar.c
 * \date        12/2015
 * \author      Mike Curtis
 * \brief       Routines for black hole bipolar implementation 
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

#ifdef BH_BIPOLAR_FEEDBACK
/*! Structure for communication during the density computation. Holds data that is sent to other processors.
 */

int blackhole_bipolar_evaluate(int target, int mode, int threadid);

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble Pos[3];
  MyFloat BH_Hsml;
  MyFloat BH_BipolarJ[3];
  int BH_BipolarColdDisk;

  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;


 /* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in * in, int i, int firstnode)
{
  in->Pos[0] = P[i].Pos[0];
  in->Pos[1] = P[i].Pos[1];
  in->Pos[2] = P[i].Pos[2];
  in->BH_Hsml = BPP(i).BH_Hsml;

  in->BH_BipolarJ[0] = BPP(i).BH_BipolarJ[0];
  in->BH_BipolarJ[1] = BPP(i).BH_BipolarJ[1];
  in->BH_BipolarJ[2] = BPP(i).BH_BipolarJ[2];
  in->BH_BipolarColdDisk = BPP(i).BH_BipolarColdDisk;

  in->Firstnode = firstnode;
}

/* local data structure that holds results acquired on remote processors */
typedef struct
{
  MyFloat BH_BipolarSum;
  MyFloat BH_Density;
  MyFloat BH_U;
  MyFloat BH_VolSum;
} data_out;

static data_out *DataResult, *DataOut;

 /* routine to store or combine result data */
static void out2particle(data_out * out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES)	/* initial store */
    {
      BPP(i).BH_BipolarSum = out->BH_BipolarSum;
      BPP(i).BH_Density = out->BH_Density;
      BPP(i).BH_U = out->BH_U;
      BPP(i).BH_VolSum = out->BH_VolSum;
    }
  else				/* merge */
    {
      BPP(i).BH_BipolarSum += out->BH_BipolarSum;
      BPP(i).BH_Density += out->BH_Density;
      BPP(i).BH_VolSum += out->BH_VolSum;
      BPP(i).BH_U += out->BH_U;
    }
}

#include "../generic_comm_helpers2.h"

static void kernel_local(void)
{
  int idx;
#ifdef GENERIC_ASYNC
  int flag = 0;
#endif
#pragma omp parallel private(i, idx)
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
	idx = NextParticle++;

	if(idx >= TimeBinsBHAccretion.NActiveParticles)
	  break;

	int i = TimeBinsBHAccretion.ActiveParticleList[idx];
	if(i < 0)
	  continue;

	if(blackhole_isactive(i))
	  blackhole_bipolar_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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

	blackhole_bipolar_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

void blackhole_bipolar(void)
{
  if(TimeBinsBHAccretion.GlobalNActiveParticles == 0)
    return;

  int nleft;
  long long ntot;

  generic_set_MaxNexport();
  printf("BLACK_HOLES: Beginning black hole bipolar calculation.");

  /* we will repeat the whole thing for those black holes where we didn't find enough neighbours */
  do
    {
      printf("BLACK_HOLES: Beginning coms...");
      generic_comm_pattern(TimeBinsBHAccretion.NActiveParticles, kernel_local, kernel_imported);


      /* do final operations on results */
      nleft = 0;

      for(int idx = 0; idx < TimeBinsBHAccretion.NActiveParticles; idx++)
	{
	  int i = TimeBinsBHAccretion.ActiveParticleList[idx];
	  if(i < 0)
	    continue;

	  if(BPP(i).BH_Density > 0)
	    {
	      BPP(i).BH_U /= BPP(i).BH_Density;
	      BPP(i).BH_Density /= BPP(i).BH_VolSum;
	    }

	  if(BPP(i).BH_Density == 0 && BPP(i).BH_BipolarColdDisk == 1)
	    {
	      /* Didn't find any cells in the cold disk - need to repeat.
	       * This is very rare so haven't optimised for now, but could just follow inactive marking used in blackhole_density.c */
	      BPP(i).BH_BipolarColdDisk = 0;
	      nleft++;
	      warn("Didn't find any cells in the cold disk. Repeating without.");
	    }
	}

      sumup_large_ints(1, &nleft, &ntot);
    }
  while(ntot > 0);
}


/*! This function represents the core of the blackhole density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
static int blackhole_bipolar_evaluate(int target, int mode, int threadid)
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

  MyFloat bipolar_sum = 0;
  MyDouble *bipolar_j = in->BH_BipolarJ;
  int cold_disk = in->BH_BipolarColdDisk;
  MyFloat rho = 0, rho_ker = 0, smooth_u = 0, vol_weight = 0;

  double xtmp, ytmp, ztmp;
  MyDouble dx, dy, dz;

  MyDouble *pos = in->Pos;
  MyFloat h = in->BH_Hsml;
  MyFloat hinv = 1.0 / h;
  MyFloat hinv3 = hinv * hinv * hinv;

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, threadid, numnodes, firstnode);

  for(int n = 0; n < nfound; n++)
    {
      int j = Thread[threadid].Ngblist[n];

      if(P[j].Mass > 0 && P[j].ID != 0)
	{
	  dx = NGB_PERIODIC_LONG_X(pos[0] - P[j].Pos[0]);
	  dy = NGB_PERIODIC_LONG_Y(pos[1] - P[j].Pos[1]);
	  dz = NGB_PERIODIC_LONG_Z(pos[2] - P[j].Pos[2]);

	  MyDouble r2 = dx * dx + dy * dy + dz * dz;
	  MyFloat mass_j = P[j].Mass;

	  if(r2 < h * h)
	    {
	      MyDouble r = sqrt(r2);
	      MyDouble u = r * hinv;

	      MyDouble wk;

	      if(u < 0.5)
		{
		  wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
		}
	      else
		{
		  wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
		}

	      if(is_cell_in_bipolar_cone(j, pos, bipolar_j))
		{
		  bipolar_sum += P[j].Mass * wk;
		}
	      else
		{
		  if(is_cell_in_bipolar_cold_disk(j, cold_disk))
		    {
		      rho += mass_j;
		      rho_ker += (mass_j * wk);
		      vol_weight += SphP[j].Volume;
		      smooth_u += (mass_j * SphP[j].Utherm);
		    }
		}
	    }
	}
    }

  out.BH_BipolarSum = bipolar_sum;
  out.BH_Density = rho;
  out.BH_U = smooth_u;
  out.BH_VolSum = vol_weight;

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}


/* take a cell id and check whether it is sitting in a blackhole bipolar cone */
int is_cell_in_bipolar_cone(MyIDType i, MyDouble * pos, MyDouble * bipolar_j)
{
  double r[3], r2;

  for(int j = 0; j < 3; j++)
    {
      r[j] = P[i].Pos[j] - pos[j];
      r2 += r[j] * r[j];
    }

  double dot = 0;

  for(int j = 0; j < 3; j++)
    {
      dot += r[j] * bipolar_j[j];
    }

  dot /= sqrt(r2);

  if(fabs(dot) > cos(All.BHBipolarTheta))
    {
      return 1;
    }

  return 0;
}

int is_cell_in_bipolar_cold_disk(MyIDType i, int ColdDisk)
{
  double u_to_temp_fac =
    (4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC))) * PROTONMASS / BOLTZMANN * GAMMA_MINUS1 * All.UnitEnergy_in_cgs /
    All.UnitMass_in_g;

  // if we have a significant cold component, we estimate from the disk
  if(ColdDisk)
    {
      if(SphP[i].Utherm * u_to_temp_fac < All.BHBipolarColdTemp)
	{
	  return 1;
	}
      else
	{
	  return 0;
	}
    }

  return 1;
}


#endif
