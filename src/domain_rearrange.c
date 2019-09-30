/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/domain_rearrange.c
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
#include <strings.h>
#include <math.h>


#include "allvars.h"
#include "proto.h"
#include "domain.h"
#include "voronoi.h"

void domain_rearrange_particle_sequence(void)
{
#if defined(USE_SFR) || defined(SINK_PARTICLES)
#ifdef USE_SFR
  if(Stars_converted)
#else
  if(SinksFormedSinceLastDomain)
#endif
    {
      struct particle_data psave;
      peanokey key;

      for(int i = 0; i < NumGas; i++)
        if(P[i].Type != 0)       /*If not a gas particle, swap to the end of the list */
          {
            psave = P[i];
            key = Key[i];

            P[i] = P[NumGas - 1];
            SphP[i] = SphP[NumGas - 1];
            Key[i] = Key[NumGas - 1];

            P[NumGas - 1] = psave;
            Key[NumGas - 1] = key;

#ifdef GFM
            if(P[NumGas - 1].Type == 4)
              StarP[P[NumGas - 1].AuxDataID].PID = NumGas - 1;
#endif
#ifdef BLACK_HOLES
            if(P[NumGas - 1].Type == 5)
              BHP[P[NumGas - 1].AuxDataID].PID = NumGas - 1;
#endif
#ifdef DUST_LIVE
            if(P[NumGas - 1].Type == DUST_LIVE)
              DustP[P[NumGas - 1].AuxDataID].PID = NumGas - 1;
#endif
            NumGas--;
            i--;
          }
      /*Now we have rearranged the particles,
       *we don't need to do it again unless there are more stars*/
#ifdef USE_SFR 
      Stars_converted = 0;
#endif
#ifdef SINK_PARTICLES
      SinksFormedSinceLastDomain = 0;
#endif
    }
#endif



#if defined(BLACK_HOLES) || defined(REFINEMENT_MERGE_CELLS) || defined(GFM_WINDS) || defined(GFM_WINDS_LOCAL) || defined(SINKS) || defined(SINK_PARTICLES) || defined(DUST_LIVE)
  int i, count_elim, count_gaselim, count_BHelim, count_windelim, count_dustelim;

  count_elim = 0;
  count_gaselim = 0;
  count_BHelim = 0;
  count_windelim = 0;
  count_dustelim = 0;

  for(i = 0; i < NumPart; i++)
#ifndef DUST_LIVE
    if((P[i].Mass == 0 && P[i].ID == 0) || (P[i].Type == 4 && P[i].Mass == 0))
#else
    if((P[i].Mass == 0 && P[i].ID == 0) || (P[i].Type == 4 && P[i].Mass == 0) || (P[i].Type == DUST_LIVE && P[i].Mass == 0))
#endif
      {
#ifdef TRACER_MC
        if(P[i].NumberOfTracers > 0)
          terminate("have found a particle with no mass but tracers\n");
#endif

        if(P[i].Type == 0)
          {
            P[i] = P[NumGas - 1];
            SphP[i] = SphP[NumGas - 1];
            Key[i] = Key[NumGas - 1];

            P[NumGas - 1] = P[NumPart - 1];
            Key[NumGas - 1] = Key[NumPart - 1];

#ifdef GFM
            if(P[NumGas - 1].Type == 4)
              StarP[P[NumGas - 1].AuxDataID].PID = NumGas - 1;
#endif
#ifdef BLACK_HOLES
            if(P[NumGas - 1].Type == 5)
              BHP[P[NumGas - 1].AuxDataID].PID = NumGas - 1;
#endif
#ifdef DUST_LIVE
            if(P[NumGas - 1].Type == DUST_LIVE)
              DustP[P[NumGas - 1].AuxDataID].PID = NumGas - 1;
#endif
            NumGas--;
            count_gaselim++;
          }
#ifdef GFM
        else if(P[i].Type == 4)
          {
            StarP[P[i].AuxDataID] = StarP[N_star - 1];
            P[StarP[N_star - 1].PID].AuxDataID = P[i].AuxDataID;

            if(i < NumPart - 1)
              {
                P[i] = P[NumPart - 1];
                Key[i] = Key[NumPart - 1];

                if(P[i].Type == 4)
                  StarP[P[i].AuxDataID].PID = i;
#ifdef BLACK_HOLES
                if(P[i].Type == 5)
                  BHP[P[i].AuxDataID].PID = i;
#endif
#ifdef DUST_LIVE
                if(P[i].Type == DUST_LIVE)
                  DustP[P[i].AuxDataID].PID = i;
#endif

              }
            N_star--;
            count_windelim++;
          }
#endif /* GFM */
#ifdef BLACK_HOLES
        else if(P[i].Type == 5)
          {
            BHP[P[i].AuxDataID] = BHP[NumBHs - 1];
            P[BHP[NumBHs - 1].PID].AuxDataID = P[i].AuxDataID;

            if(i < NumPart - 1)
              {
                P[i] = P[NumPart - 1];
                Key[i] = Key[NumPart - 1];

#ifdef GFM
                if(P[i].Type == 4)
                  StarP[P[i].AuxDataID].PID = i;
#endif
                if(P[i].Type == 5)
                  BHP[P[i].AuxDataID].PID = i;
#ifdef DUST_LIVE
                if(P[i].Type == DUST_LIVE)
                  DustP[P[i].AuxDataID].PID = i;
#endif

              }

            NumBHs--;
            count_BHelim++;
          }
#endif /* BLACK_HOLES */
#ifdef DUST_LIVE
        else if(P[i].Type == DUST_LIVE)
          {
            DustP[P[i].AuxDataID] = DustP[N_dust - 1];
            P[DustP[N_dust - 1].PID].AuxDataID = P[i].AuxDataID;

            if(i < NumPart - 1)
              {
                P[i] = P[NumPart - 1];
                Key[i] = Key[NumPart - 1];

#ifdef GFM
                if(P[i].Type == 4)
                  StarP[P[i].AuxDataID].PID = i;
#endif
#ifdef BLACK_HOLES
                if(P[i].Type == 5)
                  BHP[P[i].AuxDataID].PID = i;
#endif
                if(P[i].Type == DUST_LIVE)
                  DustP[P[i].AuxDataID].PID = i;

              }

            N_dust--;
            count_dustelim++;
          }
#endif /* DUST_LIVE */
        NumPart--;
        i--;

        count_elim++;
      }

  int count[5] = {count_elim, count_gaselim, 0, 0, 0}, tot[5], nelem = 2;
#ifdef BLACK_HOLES
  count[2] = count_BHelim;
  nelem = 3;
#endif
#ifdef GFM
  count[3] = count_windelim;
  nelem = 4;
#endif
#ifdef DUST_LIVE
  count[4] = count_dustelim;
  nelem = 5;
#endif

  MPI_Allreduce(count, tot, nelem, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
#ifndef DUST_LIVE
      printf("DOMAIN: Eliminated %d derefined/swallowed gas cells, merged away %d black holes, removed %d recoupled wind particles.\n",
            tot[1], tot[2], tot[3]);
#else
      printf("DOMAIN: Eliminated %d derefined/swallowed gas cells, merged away %d black holes, removed %d recoupled wind particles, removed %d dust particles.\n",
            tot[1], tot[2], tot[3], tot[4]);
#endif
      myflush(stdout);
    }

  All.TotNumPart -= tot[0];
  All.TotNumGas -= tot[1];
#ifdef BLACK_HOLES
  All.TotNumBHs -= tot[2];
#endif
#ifdef GFM
  All.TotN_star -= tot[3];
#endif
#ifdef DUST_LIVE
  All.TotN_dust -= tot[4];
#endif
#endif
}
