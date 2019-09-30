/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dust_live/dust_util.c
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
#include "../proto.h"
#include "../allvars.h"

#ifdef DUST_LIVE

void start_dust(void)
{
  TIMER_START(CPU_DUST);

  DustParticle = mymalloc("DustParticle", N_dust * sizeof(struct dust_particle));

  Ndust = 0;
  for(int idx = 0; idx < TimeBinsDust.NActiveParticles; idx++)
    {
      int i = TimeBinsDust.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Ti_Current != All.Ti_Current)
        {
          terminate("how can this be?");
        }

      if((P[i].Type == DUST_LIVE) && (P[i].Mass > 0))
        {
          DustParticle[Ndust].index = i;
          DustParticle[Ndust].active_idx = idx;
          DustParticle[Ndust].NumNgb = 0.0;
          DustParticle[Ndust].NormSph = 0.0;
          DustParticle[Ndust].TotNgbMass = 0.0;
          DustParticle[Ndust].Dhsmlrho = 0.0;
          for(int k = 0; k < 3; k++)
            {
              DustParticle[Ndust].LocalGasAccel[k] = 0.0;
              DustParticle[Ndust].LocalGradP[k] = 0.0;
            }
          Ndust++;
        }
    }

}

void end_dust(void)
{
  myfree(DustParticle);

  TIMER_STOP(CPU_DUST);
}

#ifdef DL_GRAIN_BINS
#if defined(DL_SNE_DESTRUCTION) || defined(DL_SHATTERING) || defined(DL_COAGULATION)
void begin_shattering(void)
{
  int idx, i;

  Nforces = 0;
  TargetList = mymalloc("TargetList", (NumPart + Tree_NumPartImported) * sizeof(int));

  for(idx = 0; idx < TimeBinsDust.NActiveParticles; idx++)
    {
      i = TimeBinsDust.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if((P[i].Type == DUST_LIVE) && (Tree_Task_list[i] == ThisTask) && (P[i].Mass > 0.0))
        {
          TargetList[Nforces++] = i;
        }
    }

  for(i = 0; i < Tree_NumPartImported; i++)
#ifndef HIERARCHICAL_GRAVITY
    if(Tree_Points[i].ActiveFlag)
#endif
      if((Tree_Points[i].Type == DUST_LIVE) && (Tree_Points[i].Mass > 0.0))
        {
          TargetList[Nforces++] = i + Tree_ImportedNodeOffset;
        }

  PShatter = (struct shatter_data_p *) mymalloc("PShatter", Nforces * sizeof(struct shatter_data_p));

  memset(&PShatter[0], 0, Nforces * sizeof(struct shatter_data_p));
}

void end_shattering(void)
{
  myfree(PShatter);
  myfree(TargetList);
}

static struct shatterimported_data
{
  int index;
  MyFloat DustDensity;
  MyFloat DustNumNgb;
  MyFloat DustHsml;
} *Tree_ShatterImported;

void exchange_shattering_results(void)
{
  int idx, i, j, k, nexport, nimport;
  int ngrp, recvTask, n;
  int ncount, ncount1, ncount2;
  struct shatterimported_data *tmp_results;

  /* first all active particles -> update all properties that are always calculated; i.e. those that also need updates even without a scattering happening */
  for(i = 0, ncount = 0; i < Tree_NumPartImported; i++)
#ifndef HIERARCHICAL_GRAVITY
    if(Tree_Points[i].ActiveFlag)
#endif
      if((Tree_Points[i].Type == DUST_LIVE) && (Tree_Points[i].Mass > 0.0))
        ncount++;

  Tree_ShatterImported = mymalloc("Tree_ShatterImported", ncount * sizeof(struct shatterimported_data));

  memset(Tree_ShatterImported, 0, ncount * sizeof(struct shatterimported_data));

  for(idx = 0, ncount1 = 0; idx < TimeBinsDust.NActiveParticles; idx++)
    {
      i = TimeBinsDust.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if((P[i].Type == DUST_LIVE) && (Tree_Task_list[i] == ThisTask) && (P[i].Mass > 0.0))
        {
          DTP(i).DustNumNgb = PShatter[ncount1].DustNumNgb;
          DTP(i).DustHsml = PShatter[ncount1].DustHsml;
          DTP(i).DustDensity = PShatter[ncount1].DustDensity;
          ncount1++;
        }
    }

  for(i = 0, ncount2 = 0; i < Tree_NumPartImported; i++)
#ifndef HIERARCHICAL_GRAVITY
    if(Tree_Points[i].ActiveFlag)
#endif
    if((Tree_Points[i].Type == DUST_LIVE) && (Tree_Points[i].Mass > 0.0))
      {
        Tree_ShatterImported[ncount2].DustNumNgb = PShatter[ncount1].DustNumNgb;
        Tree_ShatterImported[ncount2].DustHsml = PShatter[ncount1].DustHsml;
        Tree_ShatterImported[ncount2].DustDensity = PShatter[ncount1].DustDensity;
        ncount1++;
        ncount2++;
      }


  for(j = 0; j < NTask; j++)
    Recv_count[j] = 0;

  for(i = 0, n = 0, k = 0; i < NTask; i++)
    for(j = 0; j < Force_Recv_count[i]; j++, n++)
      {
#ifndef HIERARCHICAL_GRAVITY
        if(Tree_Points[n].ActiveFlag)
#endif
          if((Tree_Points[n].Type == DUST_LIVE) && (Tree_Points[n].Mass > 0.0))
            {
              Tree_ShatterImported[k].index = Tree_Points[n].index;
              Recv_count[i]++;
              k++;
            }
      }

  MPI_Alltoall(Recv_count, 1, MPI_INT, Send_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nexport = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nexport += Send_count[j];
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  tmp_results = mymalloc("tmp_results", nexport * sizeof(struct shatterimported_data));
  memset(tmp_results, -1, nexport * sizeof(struct shatterimported_data));

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          MPI_Sendrecv(&Tree_ShatterImported[Recv_offset[recvTask]],
                       Recv_count[recvTask] * sizeof(struct shatterimported_data), MPI_BYTE, recvTask, TAG_FOF_A,
                       &tmp_results[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct shatterimported_data), MPI_BYTE, recvTask, TAG_FOF_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

  for(i = 0; i < nexport; i++)
    {
      int target = tmp_results[i].index;
      DTP(target).DustDensity = tmp_results[i].DustDensity;
      DTP(target).DustNumNgb = tmp_results[i].DustNumNgb;
      DTP(target).DustHsml = tmp_results[i].DustHsml;
    }

  myfree(tmp_results);
  myfree(Tree_ShatterImported);
}
#endif
#endif

#ifdef DL_GRAIN_BINS
int elem_can_be_dust(int k)
{
  if((k == element_index_Carbon) || (k == element_index_Oxygen) || (k == element_index_Magnesium) || (k == element_index_Silicon) || (k == element_index_Iron))
    {
      return 1;
    }
  /* If GFM_NORMALIZED_METAL_ADVECTION is on, OtherMetals are assumed to not
   * deplete onto dust. */
  return 0;
}
#endif

#endif
