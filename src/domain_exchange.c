/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/domain_exchange.c
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



void domain_resize_storage(int count_get, int count_get_sph, int option_flag)
{
  int load = NumPart + count_get;
  int sphload = NumGas + count_get_sph;
  int loc_data[2] = {load, sphload}, res[2];

  MPI_Allreduce(loc_data, res, 2, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  int max_load = res[0];
  int max_sphload = res[1];

  if(max_load > (1.0 - ALLOC_TOLERANCE) * All.MaxPart || max_load < (1.0 - 3 * ALLOC_TOLERANCE) * All.MaxPart)
    {
      All.MaxPart = max_load / (1.0 - 2 * ALLOC_TOLERANCE);
      reallocate_memory_maxpart();

      if(option_flag == 1)
        Key = (peanokey *) myrealloc_movable(Key, sizeof(peanokey) * All.MaxPart);
    }

  if(max_sphload >= (1.0 - ALLOC_TOLERANCE) * All.MaxPartSph || max_sphload < (1.0 - 3 * ALLOC_TOLERANCE) * All.MaxPartSph)
    {
      All.MaxPartSph = max_sphload / (1.0 - 2 * ALLOC_TOLERANCE);
      if(option_flag == 2)
        {
          if(All.MaxPartSph > Ngb_MaxPart)
            ngb_treemodifylength(All.MaxPartSph - Ngb_MaxPart);
        }
      reallocate_memory_maxpartsph();
    }
}

#ifdef GFM
void domain_resize_storage_stars(int count_get_star)
{
  int max_starload, starload = N_star + count_get_star;
  MPI_Allreduce(&starload, &max_starload, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  if(max_starload >= All.MaxPartStar - 0.5 * ALLOC_STARBH_ROOM * (All.TotNumGas / NTask))
    {
      All.MaxPartStar = max_starload + ALLOC_STARBH_ROOM * (All.TotNumGas / NTask);
      reallocate_memory_maxpartstar();
    }
}
#endif

#ifdef BLACK_HOLES
void domain_resize_storage_blackholes(int count_get_BHs)
{
  int max_BHsload, BHsload = NumBHs + count_get_BHs;
  MPI_Allreduce(&BHsload, &max_BHsload, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  if(max_BHsload >= All.MaxPartBHs - 0.5 * ALLOC_STARBH_ROOM * (All.TotNumGas / NTask))
    {
      All.MaxPartBHs = max_BHsload + ALLOC_STARBH_ROOM * (All.TotNumGas / NTask);
      reallocate_memory_maxpartBHs();
    }
}
#endif

#ifdef DUST_LIVE
void domain_resize_storage_dust(int count_get_dust)
{
  int max_dustload, dustload = N_dust + count_get_dust;
  MPI_Allreduce(&dustload, &max_dustload, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  if(max_dustload >= All.MaxPartDust - 0.5 * ALLOC_STARBH_ROOM * (All.TotNumGas / NTask))
    {
      All.MaxPartDust = max_dustload + ALLOC_STARBH_ROOM * (All.TotNumGas / NTask);
      reallocate_memory_maxpartdust();
    }
}
#endif

#ifdef TRACER_MC
void domain_resize_storage_tracer(int count_get_tracer)
{
  int max_N_tracer, N_tracerload = N_tracer + count_get_tracer;
  MPI_Allreduce(&N_tracerload, &max_N_tracer, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  if(max_N_tracer >= All.MaxPartTracer - sqrt(All.MaxPartTracer))
    {
      int delta = max_N_tracer - All.MaxPartTracer + 2 * sqrt(All.MaxPartTracer);
      reallocate_memory_maxparttracer(delta);
    }
}
#endif



void domain_exchange(void)
{
  double t0 = second();

  int count_togo = 0, count_togo_sph = 0, count_get = 0, count_get_sph = 0;
  int *count, *count_sph, *offset, *offset_sph;
  int *count_recv, *count_recv_sph, *offset_recv, *offset_recv_sph;
  int i, n, no, target;
  struct particle_data *partBuf;
  struct sph_particle_data *sphBuf;

#ifdef TRACER_MC
  int next;
  int count_togo_tracer = 0, count_get_tracer = 0;
  int *count_tracer, *offset_tracer;
  int *count_recv_tracer, *offset_recv_tracer;
  struct tracer_linked_list *tracerBuf, *tracerRecvBuf;
#endif
#ifdef GFM
  int count_togo_star = 0, count_get_star = 0;
  int *count_star, *offset_star;
  int *count_recv_star, *offset_recv_star;
  struct star_particle_data *starBuf;
  int *starbuf_offset, *starbuf_recv_offset;
#endif
#ifdef BLACK_HOLES
  int count_togo_BHs = 0, count_get_BHs = 0;
  int *count_BHs, *offset_BHs;
  int *count_recv_BHs, *offset_recv_BHs;
  struct bh_particle_data *BHsBuf;
  int *bhbuf_offset, *bhbuf_recv_offset;
#endif
#ifdef DUST_LIVE
  int count_togo_dust = 0, count_get_dust = 0;
  int *count_dust, *offset_dust;
  int *count_recv_dust, *offset_recv_dust;
  struct dust_particle_data *dustBuf;
  int *dustbuf_offset, *dustbuf_recv_offset;
#endif
  peanokey *keyBuf;

  long long sumtogo = 0;

  for(i = 0; i < NTask; i++)
    sumtogo += toGo[i];

  sumup_longs(1, &sumtogo, &sumtogo);

#ifdef TRACER_MC_CHECKS
  check_tracer_lists();
#endif

  count = (int *) mymalloc_movable(&count, "count", NTask * sizeof(int));
  count_sph = (int *) mymalloc_movable(&count_sph, "count_sph", NTask * sizeof(int));
  offset = (int *) mymalloc_movable(&offset, "offset", NTask * sizeof(int));
  offset_sph = (int *) mymalloc_movable(&offset_sph, "offset_sph", NTask * sizeof(int));
  count_recv = (int *) mymalloc_movable(&count_recv, "count_recv", NTask * sizeof(int));
  count_recv_sph = (int *) mymalloc_movable(&count_recv_sph, "count_recv_sph", NTask * sizeof(int));
  offset_recv = (int *) mymalloc_movable(&offset_recv, "offset_recv", NTask * sizeof(int));
  offset_recv_sph = (int *) mymalloc_movable(&offset_recv_sph, "offset_recv_sph", NTask * sizeof(int));

#ifdef TRACER_MC
  count_tracer = (int *) mymalloc_movable(&count_tracer, "count_tracer", NTask * sizeof(int));
  offset_tracer = (int *) mymalloc_movable(&offset_tracer, "offset_tracer", NTask * sizeof(int));
  count_recv_tracer = (int *) mymalloc_movable(&count_recv_tracer, "count_recv_tracer", NTask * sizeof(int));
  offset_recv_tracer = (int *) mymalloc_movable(&offset_recv_tracer, "offset_recv_tracer", NTask * sizeof(int));
#endif
#ifdef GFM
  count_star = (int *) mymalloc_movable(&count_star, "count_star", NTask * sizeof(int));
  offset_star = (int *) mymalloc_movable(&offset_star, "offset_star", NTask * sizeof(int));
  count_recv_star = (int *) mymalloc_movable(&count_recv_star, "count_recv_star", NTask * sizeof(int));
  offset_recv_star = (int *) mymalloc_movable(&offset_recv_star, "offset_recv_star", NTask * sizeof(int));
  starbuf_offset = (int *) mymalloc_movable(&starbuf_offset, "starbuf_offset", NTask * sizeof(int));
  starbuf_recv_offset = (int *) mymalloc_movable(&starbuf_recv_offset, "starbuf_recv_offset", NTask * sizeof(int));
#endif
#ifdef BLACK_HOLES
  count_BHs = (int *) mymalloc_movable(&count_BHs, "count_BHs", NTask * sizeof(int));
  offset_BHs = (int *) mymalloc_movable(&offset_BHs, "offset_BHs", NTask * sizeof(int));
  count_recv_BHs = (int *) mymalloc_movable(&count_recv_BHs, "count_recv_BHs", NTask * sizeof(int));
  offset_recv_BHs = (int *) mymalloc_movable(&offset_recv_BHs, "offset_recv_BHs", NTask * sizeof(int));
  bhbuf_offset = (int *) mymalloc_movable(&bhbuf_offset, "bhbuf_offset", NTask * sizeof(int));
  bhbuf_recv_offset = (int *) mymalloc_movable(&bhbuf_recv_offset, "bhbuf_recv_offset", NTask * sizeof(int));
#endif
#ifdef DUST_LIVE
  count_dust = (int *) mymalloc_movable(&count_dust, "count_dust", NTask * sizeof(int));
  offset_dust = (int *) mymalloc_movable(&offset_dust, "offset_dust", NTask * sizeof(int));
  count_recv_dust = (int *) mymalloc_movable(&count_recv_dust, "count_recv_dust", NTask * sizeof(int));
  offset_recv_dust = (int *) mymalloc_movable(&offset_recv_dust, "offset_recv_dust", NTask * sizeof(int));
  dustbuf_offset = (int *) mymalloc_movable(&dustbuf_offset, "dustbuf_offset", NTask * sizeof(int));
  dustbuf_recv_offset = (int *) mymalloc_movable(&dustbuf_recv_offset, "dustbuf_recv_offset", NTask * sizeof(int));
#endif

  int prec_offset;
  int *decrease;

  decrease = (int *) mymalloc_movable(&decrease, "decrease", NTask * sizeof(int));

  for(i = 1, offset_sph[0] = 0, decrease[0] = 0; i < NTask; i++)
    {
      offset_sph[i] = offset_sph[i - 1] + toGoSph[i - 1];
      decrease[i] = toGoSph[i - 1];
    }

  prec_offset = offset_sph[NTask - 1] + toGoSph[NTask - 1];

#ifdef TRACER_MC
  offset_tracer[0] = 0;
  for(i = 1; i < NTask; i++)
    offset_tracer[i] = offset_tracer[i - 1] + toGoTracer[i - 1];
#endif

#ifdef GFM
  offset_star[0] = prec_offset;
  for(i = 1; i < NTask; i++)
    {
      offset_star[i] = offset_star[i - 1] + toGoStar[i - 1];
      decrease[i] += toGoStar[i - 1];
    }
  prec_offset = offset_star[NTask - 1] + toGoStar[NTask - 1];
#endif

#ifdef BLACK_HOLES
  offset_BHs[0] = prec_offset;
  for(i = 1; i < NTask; i++)
    {
      offset_BHs[i] = offset_BHs[i - 1] + toGoBHs[i - 1];
      decrease[i] += toGoBHs[i - 1];
    }
  prec_offset = offset_BHs[NTask - 1] + toGoBHs[NTask - 1];
#endif

#ifdef DUST_LIVE
  offset_dust[0] = prec_offset;
  for(i = 1; i < NTask; i++)
    {
      offset_dust[i] = offset_dust[i - 1] + toGoDust[i - 1];
      decrease[i] += toGoDust[i - 1];
    }
  prec_offset = offset_dust[NTask - 1] + toGoDust[NTask - 1];
#endif

  offset[0] = prec_offset;
  for(i = 1; i < NTask; i++)
    offset[i] = offset[i - 1] + (toGo[i - 1] - decrease[i]);

  myfree(decrease);

  for(i = 0; i < NTask; i++)
    {
      count_togo += toGo[i];
      count_togo_sph += toGoSph[i];
      count_get += toGet[i];
      count_get_sph += toGetSph[i];
#ifdef TRACER_MC
      count_togo_tracer += toGoTracer[i];
      count_get_tracer += toGetTracer[i];
#endif
#ifdef GFM
      count_togo_star += toGoStar[i];
      count_get_star += toGetStar[i];
#endif
#ifdef BLACK_HOLES
      count_togo_BHs += toGoBHs[i];
      count_get_BHs += toGetBHs[i];
#endif
#ifdef DUST_LIVE
      count_togo_dust += toGoDust[i];
      count_get_dust += toGetDust[i];
#endif
    }


  partBuf = (struct particle_data *) mymalloc_movable(&partBuf, "partBuf", count_togo * sizeof(struct particle_data));
  sphBuf = (struct sph_particle_data *) mymalloc_movable(&sphBuf, "sphBuf", count_togo_sph * sizeof(struct sph_particle_data));
#ifdef TRACER_MC
  tracerBuf = (struct tracer_linked_list *) mymalloc_movable(&tracerBuf, "tracerBuf", count_togo_tracer * sizeof(struct tracer_linked_list));
  tracerRecvBuf = (struct tracer_linked_list *) mymalloc_movable(&tracerRecvBuf, "tracerRecvBuf", count_get_tracer * sizeof(struct tracer_linked_list));
#endif
#ifdef GFM
  starBuf = (struct star_particle_data *) mymalloc_movable(&starBuf, "starBuf", count_togo_star * sizeof(struct star_particle_data));
#endif
#ifdef BLACK_HOLES
  BHsBuf = (struct bh_particle_data *) mymalloc_movable(&BHsBuf, "BHsBuf", count_togo_BHs * sizeof(struct bh_particle_data));
#endif
#ifdef DUST_LIVE
  dustBuf = (struct dust_particle_data *) mymalloc_movable(&dustBuf, "dustBuf", count_togo_dust * sizeof(struct dust_particle_data));
#endif

  keyBuf = (peanokey *) mymalloc_movable(&keyBuf, "keyBuf", count_togo * sizeof(peanokey));


  for(i = 0; i < NTask; i++)
    {
      count[i] = count_sph[i] = 0;
#ifdef TRACER_MC
      count_tracer[i] = 0;
#endif
#ifdef GFM
      count_star[i] = 0;
#endif
#ifdef BLACK_HOLES
      count_BHs[i] = 0;
#endif
#ifdef DUST_LIVE
      count_dust[i] = 0;
#endif
    }


  for(n = 0; n < NumPart; n++)
    {
      no = 0;

      peanokey mask = ((peanokey) 7) << (3 * (BITS_PER_DIMENSION - 1));
      int shift = 3 * (BITS_PER_DIMENSION - 1);

      while(topNodes[no].Daughter >= 0)
        {
          no = topNodes[no].Daughter + (int) ((Key[n] & mask) >> shift);
          mask >>= 3;
          shift -= 3;
        }

      no = topNodes[no].Leaf;

      target = DomainTask[no];

#ifdef TRACER_MC
      P[n].OriginTask = ThisTask;
#endif

      if(target != ThisTask)
        {
          /* copy this particle into the exchange buffer */
          if(P[n].Type == 0)
            {
              partBuf[offset_sph[target] + count_sph[target]] = P[n];
#ifdef TRACER_MC
              partBuf[offset_sph[target] + count_sph[target]].TracerHead = count_tracer[target];
#endif
              keyBuf[offset_sph[target] + count_sph[target]] = Key[n];
              sphBuf[offset_sph[target] + count_sph[target]] = SphP[n];
              count_sph[target]++;
            }
#ifdef GFM
          else if(P[n].Type == 4)
            {
              partBuf[offset_star[target] + count_star[target]] = P[n];
#ifdef TRACER_MC
              partBuf[offset_star[target] + count_star[target]].TracerHead = count_tracer[target];
#endif
              keyBuf[offset_star[target] + count_star[target]] = Key[n];
              starBuf[offset_star[target] - offset_star[0] + count_star[target]] = StarP[P[n].AuxDataID];
              count_star[target]++;
            }
#endif
#ifdef BLACK_HOLES
          else if(P[n].Type == 5)
            {
              partBuf[offset_BHs[target] + count_BHs[target]] = P[n];
#ifdef TRACER_MC
              partBuf[offset_BHs[target] + count_BHs[target]].TracerHead = count_tracer[target];
#endif
              keyBuf[offset_BHs[target] + count_BHs[target]] = Key[n];
              BHsBuf[offset_BHs[target] - offset_BHs[0] + count_BHs[target]] = BHP[P[n].AuxDataID];
              count_BHs[target]++;
            }
#endif
#ifdef DUST_LIVE
          else if(P[n].Type == DUST_LIVE)
            {
              partBuf[offset_dust[target] + count_dust[target]] = P[n];
#ifdef TRACER_MC
              partBuf[offset_dust[target] + count_dust[target]].TracerHead = count_tracer[target];
#endif
              keyBuf[offset_dust[target] + count_dust[target]] = Key[n];
              dustBuf[offset_dust[target] - offset_dust[0] + count_dust[target]] = DustP[P[n].AuxDataID];
              count_dust[target]++;
            }
#endif
          else
            {
              partBuf[offset[target] + count[target]] = P[n];
#ifdef TRACER_MC
              partBuf[offset[target] + count[target]].TracerHead = count_tracer[target];
#endif
              keyBuf[offset[target] + count[target]] = Key[n];
              count[target]++;
            }

#ifdef TRACER_MC
          next = P[n].TracerHead;

          while(next >= 0)      /* fill tracerBuf export buffer and tag entries for removal */
            {
              tracerBuf[offset_tracer[target] + count_tracer[target]++] = TracerLinkedList[next];
              remove_tracer_from_parent(n, next);
              release_tracer_slot(next);

              next = TracerLinkedList[next].Next;
            }

#endif /* TRACER_MC */


          if(P[n].Type == 0)
            {
              P[n] = P[NumGas - 1];
              P[NumGas - 1] = P[NumPart - 1];
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

              Key[n] = Key[NumGas - 1];
              Key[NumGas - 1] = Key[NumPart - 1];

              SphP[n] = SphP[NumGas - 1];

              NumGas--;
            }
#ifdef GFM
          else if(P[n].Type == 4)
            {
              StarP[P[n].AuxDataID] = StarP[N_star - 1];
              P[StarP[N_star - 1].PID].AuxDataID = P[n].AuxDataID;

              if(n < NumPart - 1)
                {
                  P[n] = P[NumPart - 1];
                  Key[n] = Key[NumPart - 1];
#ifdef GFM
                  if(P[n].Type == 4)
                    StarP[P[n].AuxDataID].PID = n;
#endif
#ifdef BLACK_HOLES
                  if(P[n].Type == 5)
                    BHP[P[n].AuxDataID].PID = n;
#endif
#ifdef DUST_LIVE
                  if(P[n].Type == DUST_LIVE)
                    DustP[P[n].AuxDataID].PID = n;
#endif
                }

              N_star--;
            }
#endif
#ifdef BLACK_HOLES
          else if(P[n].Type == 5)
            {
              BHP[P[n].AuxDataID] = BHP[NumBHs - 1];
              P[BHP[NumBHs - 1].PID].AuxDataID = P[n].AuxDataID;

              if(n < NumPart - 1)
                {
                  P[n] = P[NumPart - 1];
                  Key[n] = Key[NumPart - 1];
#ifdef GFM
                  if(P[n].Type == 4)
                    StarP[P[n].AuxDataID].PID = n;
#endif
#ifdef BLACK_HOLES
                  if(P[n].Type == 5)
                    BHP[P[n].AuxDataID].PID = n;
#endif
#ifdef DUST_LIVE
                  if(P[n].Type == DUST_LIVE)
                    DustP[P[n].AuxDataID].PID = n;
#endif
                }

              NumBHs--;
            }
#endif
#ifdef DUST_LIVE
          else if(P[n].Type == DUST_LIVE)
            {
              DustP[P[n].AuxDataID] = DustP[N_dust - 1];
              P[DustP[N_dust - 1].PID].AuxDataID = P[n].AuxDataID;

              if(n < NumPart - 1)
                {
                  P[n] = P[NumPart - 1];
                  Key[n] = Key[NumPart - 1];
#ifdef GFM
                  if(P[n].Type == 4)
                    StarP[P[n].AuxDataID].PID = n;
#endif
#ifdef BLACK_HOLES
                  if(P[n].Type == 5)
                    BHP[P[n].AuxDataID].PID = n;
#endif
#ifdef DUST_LIVE
                  if(P[n].Type == DUST_LIVE)
                    DustP[P[n].AuxDataID].PID = n;
#endif
                }

              N_dust--;
            }
#endif
          else
            {
              P[n] = P[NumPart - 1];
              Key[n] = Key[NumPart - 1];
#ifdef GFM
              if(P[n].Type == 4)
                StarP[P[n].AuxDataID].PID = n;
#endif
#ifdef BLACK_HOLES
              if(P[n].Type == 5)
                BHP[P[n].AuxDataID].PID = n;
#endif
#ifdef DUST_LIVE
              if(P[n].Type == DUST_LIVE)
                DustP[P[n].AuxDataID].PID = n;
#endif
            }

          NumPart--;
          n--;

        }                       /* target != ThisTask */
    }                           /* n < NumPart */



  /**** now resize the storage for the P[] and SphP[] arrays if needed ****/
  domain_resize_storage(count_get, count_get_sph, 1);

#ifdef GFM
  domain_resize_storage_stars(count_get_star);
#endif

#ifdef BLACK_HOLES
  domain_resize_storage_blackholes(count_get_BHs);
#endif

#ifdef DUST_LIVE
  domain_resize_storage_dust(count_get_dust);
#endif

#ifdef TRACER_MC
  domain_resize_storage_tracer(count_get_tracer);
#endif


  /*****  space has been created, now can do the actual exchange *****/

  int count_totget = count_get_sph;
#ifdef GFM
  count_totget += count_get_star;
#endif
#ifdef BLACK_HOLES
  count_totget += count_get_BHs;
#endif
#ifdef DUST_LIVE
  count_totget += count_get_dust;
#endif

  if(count_totget)
    {
      memmove(P + NumGas + count_totget, P + NumGas, (NumPart - NumGas) * sizeof(struct particle_data));
      memmove(Key + NumGas + count_totget, Key + NumGas, (NumPart - NumGas) * sizeof(peanokey));
#ifdef GFM
      for(i = 0; i < N_star; i++)
        {
          StarP[i].PID += count_totget;
        }
#endif
#ifdef BLACK_HOLES
      for(i = 0; i < NumBHs; i++)
        {
          BHP[i].PID += count_totget;
        }
#endif
#ifdef DUST_LIVE
      for(i = 0; i < N_dust; i++)
        {
          DustP[i].PID += count_totget;
        }
#endif
    }


  for(i = 0; i < NTask; i++)
    {
      count_recv_sph[i] = toGetSph[i];
      count_recv[i] = toGet[i] - toGetSph[i];
#ifdef GFM
      count_recv_star[i] = toGetStar[i];
      count_recv[i] -= toGetStar[i];
#endif
#ifdef BLACK_HOLES
      count_recv_BHs[i] = toGetBHs[i];
      count_recv[i] -= toGetBHs[i];
#endif
#ifdef DUST_LIVE
      count_recv_dust[i] = toGetDust[i];
      count_recv[i] -= toGetDust[i];
#endif

#ifdef TRACER_MC
      count_recv_tracer[i] = toGetTracer[i];
#endif
    }


  int prec_count;
  for(i = 1, offset_recv_sph[0] = NumGas; i < NTask; i++)
    offset_recv_sph[i] = offset_recv_sph[i - 1] + count_recv_sph[i - 1];
  prec_count = NumGas + count_get_sph;

#ifdef TRACER_MC
  offset_recv_tracer[0] = 0;
  for(i = 1; i < NTask; i++)
    offset_recv_tracer[i] = offset_recv_tracer[i - 1] + count_recv_tracer[i - 1];
#endif

#ifdef GFM
  offset_recv_star[0] = prec_count;
  for(i = 1; i < NTask; i++)
    offset_recv_star[i] = offset_recv_star[i - 1] + count_recv_star[i - 1];
  prec_count += count_get_star;

  for(target = 0; target < NTask; target++)
    {
      starbuf_offset[target] = offset_star[target] - offset_sph[NTask - 1] - count_sph[NTask - 1];
      starbuf_recv_offset[target] = N_star + offset_recv_star[target] - offset_recv_sph[NTask - 1] - count_recv_sph[NTask - 1];
    }
#endif

#ifdef BLACK_HOLES
  offset_recv_BHs[0] = prec_count;
  for(i = 1; i < NTask; i++)
    offset_recv_BHs[i] = offset_recv_BHs[i - 1] + count_recv_BHs[i - 1];
  prec_count += count_get_BHs;

  for(target = 0; target < NTask; target++)
    {
#ifdef GFM
      bhbuf_offset[target] = offset_BHs[target] - offset_star[NTask - 1] - count_star[NTask - 1];
      bhbuf_recv_offset[target] = NumBHs + offset_recv_BHs[target] - offset_recv_star[NTask - 1] - count_recv_star[NTask - 1];
#else
      bhbuf_offset[target] = offset_BHs[target] - offset_sph[NTask - 1] - count_sph[NTask - 1];
      bhbuf_recv_offset[target] = NumBHs + offset_recv_BHs[target] - offset_recv_sph[NTask - 1] - count_recv_sph[NTask - 1];
#endif
    }

#endif

#ifdef DUST_LIVE
  offset_recv_dust[0] = prec_count;
  for(i = 1; i < NTask; i++)
    offset_recv_dust[i] = offset_recv_dust[i - 1] + count_recv_dust[i - 1];
  prec_count += count_get_dust;

  for(target = 0; target < NTask; target++)
    {
#ifdef BLACK_HOLES
      dustbuf_offset[target] = offset_dust[target] - offset_BHs[NTask - 1] - count_BHs[NTask - 1];
      dustbuf_recv_offset[target] = N_dust + offset_recv_dust[target] - offset_recv_BHs[NTask - 1] - count_recv_BHs[NTask - 1];
#elif defined(GFM)
      dustbuf_offset[target] = offset_dust[target] - offset_star[NTask - 1] - count_star[NTask - 1];
      dustbuf_recv_offset[target] = N_dust + offset_recv_dust[target] - offset_recv_star[NTask - 1] - count_recv_star[NTask - 1];
#else
      dustbuf_offset[target] = offset_dust[target] - offset_sph[NTask - 1] - count_sph[NTask - 1];
      dustbuf_recv_offset[target] = N_dust + offset_recv_dust[target] - offset_recv_sph[NTask - 1] - count_recv_sph[NTask - 1];
#endif
    }
#endif

  offset_recv[0] = NumPart - NumGas + prec_count;

  for(i = 1; i < NTask; i++)
    offset_recv[i] = offset_recv[i - 1] + count_recv[i - 1];

#ifndef USE_MPIALLTOALLV_IN_DOMAINDECOMP

  int ngrp;
#ifdef NO_ISEND_IRECV_IN_DOMAIN /* synchronous communication */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      target = ThisTask ^ ngrp;

      if(target < NTask)
        {
          if(count_sph[target] > 0 || count_recv_sph[target] > 0)
            {
              MPI_Sendrecv(partBuf + offset_sph[target], count_sph[target] * sizeof(struct particle_data),
                           MPI_BYTE, target, TAG_PDATA_SPH,
                           P + offset_recv_sph[target], count_recv_sph[target] * sizeof(struct particle_data), MPI_BYTE, target, TAG_PDATA_SPH, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

              MPI_Sendrecv(sphBuf + offset_sph[target], count_sph[target] * sizeof(struct sph_particle_data),
                           MPI_BYTE, target, TAG_SPHDATA,
                           SphP + offset_recv_sph[target], count_recv_sph[target] * sizeof(struct sph_particle_data), MPI_BYTE, target, TAG_SPHDATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

              MPI_Sendrecv(keyBuf + offset_sph[target], count_sph[target] * sizeof(peanokey),
                           MPI_BYTE, target, TAG_KEY_SPH, Key + offset_recv_sph[target], count_recv_sph[target] * sizeof(peanokey), MPI_BYTE, target, TAG_KEY_SPH, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

#ifdef TRACER_MC
          if(count_tracer[target] > 0 || count_recv_tracer[target] > 0)
            {
              MPI_Sendrecv(tracerBuf + offset_tracer[target], count_tracer[target] * sizeof(struct tracer_linked_list),
                           MPI_BYTE, target, TAG_TRACERDATA,
                           tracerRecvBuf + offset_recv_tracer[target],
                           count_recv_tracer[target] * sizeof(struct tracer_linked_list), MPI_BYTE, target, TAG_TRACERDATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
#endif
#ifdef GFM
          if(count_star[target] > 0 || count_recv_star[target] > 0)
            {
              MPI_Sendrecv(partBuf + offset_star[target], count_star[target] * sizeof(struct particle_data),
                           MPI_BYTE, target, TAG_PDATA_STAR,
                           P + offset_recv_star[target], count_recv_star[target] * sizeof(struct particle_data), MPI_BYTE, target, TAG_PDATA_STAR, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

              MPI_Sendrecv(starBuf + offset_star[target] - offset_sph[NTask - 1] - count_sph[NTask - 1],
                           count_star[target] * sizeof(struct star_particle_data), MPI_BYTE, target, TAG_STARDATA,
                           StarP + N_star + offset_recv_star[target] - offset_recv_sph[NTask - 1] - count_recv_sph[NTask - 1],
                           count_recv_star[target] * sizeof(struct star_particle_data), MPI_BYTE, target, TAG_STARDATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

              MPI_Sendrecv(keyBuf + offset_star[target], count_star[target] * sizeof(peanokey),
                           MPI_BYTE, target, TAG_KEY_STAR,
                           Key + offset_recv_star[target], count_recv_star[target] * sizeof(peanokey), MPI_BYTE, target, TAG_KEY_STAR, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
#endif
#ifdef BLACK_HOLES
          if(count_BHs[target] > 0 || count_recv_BHs[target] > 0)
            {
              MPI_Sendrecv(partBuf + offset_BHs[target], count_BHs[target] * sizeof(struct particle_data),
                           MPI_BYTE, target, TAG_PDATA_BH,
                           P + offset_recv_BHs[target], count_recv_BHs[target] * sizeof(struct particle_data), MPI_BYTE, target, TAG_PDATA_BH, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

#ifdef GFM
              MPI_Sendrecv(BHsBuf + offset_BHs[target] - offset_star[NTask - 1] - count_star[NTask - 1],
                           count_BHs[target] * sizeof(struct bh_particle_data), MPI_BYTE, target, TAG_BHDATA,
                           BHP + NumBHs + offset_recv_BHs[target] - offset_recv_star[NTask - 1] - count_recv_star[NTask - 1],
                           count_recv_BHs[target] * sizeof(struct bh_particle_data), MPI_BYTE, target, TAG_BHDATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#else
              MPI_Sendrecv(BHsBuf + offset_BHs[target] - offset_sph[NTask - 1] - count_sph[NTask - 1],
                           count_BHs[target] * sizeof(struct bh_particle_data), MPI_BYTE, target, TAG_BHDATA,
                           BHP + NumBHs + offset_recv_BHs[target] - offset_recv_sph[NTask - 1] - count_recv_sph[NTask - 1],
                           count_recv_BHs[target] * sizeof(struct bh_particle_data), MPI_BYTE, target, TAG_BHDATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#endif

              MPI_Sendrecv(keyBuf + offset_BHs[target], count_BHs[target] * sizeof(peanokey),
                           MPI_BYTE, target, TAG_KEY_BH, Key + offset_recv_BHs[target], count_recv_BHs[target] * sizeof(peanokey), MPI_BYTE, target, TAG_KEY_BH, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
#endif
#ifdef DUST_LIVE
          if(count_dust[target] > 0 || count_recv_dust[target] > 0)
            {
              MPI_Sendrecv(partBuf + offset_dust[target], count_dust[target] * sizeof(struct particle_data),
                           MPI_BYTE, target, TAG_PDATA_DUST,
                           P + offset_recv_dust[target], count_recv_dust[target] * sizeof(struct particle_data), MPI_BYTE, target, TAG_PDATA_DUST, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

#ifdef BLACK_HOLES
              MPI_Sendrecv(dustBuf + offset_dust[target] - offset_BHs[NTask - 1] - count_BHs[NTask - 1],
                           count_dust[target] * sizeof(struct dust_particle_data), MPI_BYTE, target, TAG_DUSTDATA,
                           DustP + N_dust + offset_recv_dust[target] - offset_recv_BHs[NTask - 1] - count_recv_BHs[NTask - 1],
                           count_recv_dust[target] * sizeof(struct dust_particle_data), MPI_BYTE, target, TAG_DUSTDATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#elif defined(GFM)
              MPI_Sendrecv(dustBuf + offset_dust[target] - offset_star[NTask - 1] - count_star[NTask - 1],
                           count_dust[target] * sizeof(struct dust_particle_data), MPI_BYTE, target, TAG_DUSTDATA,
                           DustP + N_dust + offset_recv_dust[target] - offset_recv_star[NTask - 1] - count_recv_star[NTask - 1],
                           count_recv_dust[target] * sizeof(struct dust_particle_data), MPI_BYTE, target, TAG_DUSTDATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#else
              MPI_Sendrecv(dustBuf + offset_dust[target] - offset_sph[NTask - 1] - count_sph[NTask - 1],
                           count_dust[target] * sizeof(struct dust_particle_data), MPI_BYTE, target, TAG_DUSTDATA,
                           DustP + N_dust + offset_recv_dust[target] - offset_recv_sph[NTask - 1] - count_recv_sph[NTask - 1],
                           count_recv_dust[target] * sizeof(struct dust_particle_data), MPI_BYTE, target, TAG_DUSTDATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#endif

              MPI_Sendrecv(keyBuf + offset_dust[target], count_dust[target] * sizeof(peanokey),
                           MPI_BYTE, target, TAG_KEY_DUST, Key + offset_recv_dust[target], count_recv_dust[target] * sizeof(peanokey), MPI_BYTE, target, TAG_KEY_DUST, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
#endif

          if(count[target] > 0 || count_recv[target] > 0)
            {
              MPI_Sendrecv(partBuf + offset[target], count[target] * sizeof(struct particle_data),
                           MPI_BYTE, target, TAG_PDATA, P + offset_recv[target], count_recv[target] * sizeof(struct particle_data), MPI_BYTE, target, TAG_PDATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

              MPI_Sendrecv(keyBuf + offset[target], count[target] * sizeof(peanokey),
                           MPI_BYTE, target, TAG_KEY, Key + offset_recv[target], count_recv[target] * sizeof(peanokey), MPI_BYTE, target, TAG_KEY, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

#else
  /* asynchronous communication */

  MPI_Request *requests = (MPI_Request *) mymalloc_movable(&requests, "requests", 30 * NTask * sizeof(MPI_Request));
  int n_requests = 0;

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      target = ThisTask ^ ngrp;

      if(target < NTask)
        {
          if(count_recv_sph[target] > 0)
            {
              MPI_Irecv(P + offset_recv_sph[target], count_recv_sph[target] * sizeof(struct particle_data), MPI_BYTE, target, TAG_PDATA_SPH, MPI_COMM_WORLD, &requests[n_requests++]);

              MPI_Irecv(SphP + offset_recv_sph[target], count_recv_sph[target] * sizeof(struct sph_particle_data), MPI_BYTE, target, TAG_SPHDATA, MPI_COMM_WORLD, &requests[n_requests++]);

              MPI_Irecv(Key + offset_recv_sph[target], count_recv_sph[target] * sizeof(peanokey), MPI_BYTE, target, TAG_KEY_SPH, MPI_COMM_WORLD, &requests[n_requests++]);
            }

#ifdef TRACER_MC
          if(count_recv_tracer[target] > 0)
            {
              MPI_Irecv(tracerRecvBuf + offset_recv_tracer[target],
                        count_recv_tracer[target] * sizeof(struct tracer_linked_list), MPI_BYTE, target, TAG_TRACERDATA, MPI_COMM_WORLD, &requests[n_requests++]);
            }
#endif
#ifdef GFM
          if(count_recv_star[target] > 0)
            {
              MPI_Irecv(P + offset_recv_star[target], count_recv_star[target] * sizeof(struct particle_data), MPI_BYTE, target, TAG_PDATA_STAR, MPI_COMM_WORLD, &requests[n_requests++]);

              MPI_Irecv(StarP + N_star + offset_recv_star[target] - offset_recv_sph[NTask - 1] - count_recv_sph[NTask - 1],
                        count_recv_star[target] * sizeof(struct star_particle_data), MPI_BYTE, target, TAG_STARDATA, MPI_COMM_WORLD, &requests[n_requests++]);

              MPI_Irecv(Key + offset_recv_star[target], count_recv_star[target] * sizeof(peanokey), MPI_BYTE, target, TAG_KEY_STAR, MPI_COMM_WORLD, &requests[n_requests++]);
            }
#endif
#ifdef BLACK_HOLES
          if(count_recv_BHs[target] > 0)
            {
              MPI_Irecv(P + offset_recv_BHs[target], count_recv_BHs[target] * sizeof(struct particle_data), MPI_BYTE, target, TAG_PDATA_BH, MPI_COMM_WORLD, &requests[n_requests++]);

              MPI_Irecv(Key + offset_recv_BHs[target], count_recv_BHs[target] * sizeof(peanokey), MPI_BYTE, target, TAG_KEY_BH, MPI_COMM_WORLD, &requests[n_requests++]);

#ifdef GFM
              MPI_Irecv(BHP + NumBHs + offset_recv_BHs[target] - offset_recv_star[NTask - 1] -
                        count_recv_star[NTask - 1], count_recv_BHs[target] * sizeof(struct bh_particle_data), MPI_BYTE, target, TAG_BHDATA, MPI_COMM_WORLD, &requests[n_requests++]);
#else
              MPI_Irecv(BHP + NumBHs + offset_recv_BHs[target] - offset_recv_sph[NTask - 1] -
                        count_recv_sph[NTask - 1], count_recv_BHs[target] * sizeof(struct bh_particle_data), MPI_BYTE, target, TAG_BHDATA, MPI_COMM_WORLD, &requests[n_requests++]);
#endif
            }
#endif
#ifdef DUST_LIVE
          if(count_recv_dust[target] > 0)
            {
              MPI_Irecv(P + offset_recv_dust[target], count_recv_dust[target] * sizeof(struct particle_data), MPI_BYTE, target, TAG_PDATA_DUST, MPI_COMM_WORLD, &requests[n_requests++]);

              MPI_Irecv(Key + offset_recv_dust[target], count_recv_dust[target] * sizeof(peanokey), MPI_BYTE, target, TAG_KEY_DUST, MPI_COMM_WORLD, &requests[n_requests++]);

#ifdef BLACK_HOLES
              MPI_Irecv(DustP + N_dust + offset_recv_dust[target] - offset_recv_BHs[NTask - 1] -
                        count_recv_BHs[NTask - 1], count_recv_dust[target] * sizeof(struct dust_particle_data), MPI_BYTE, target, TAG_DUSTDATA, MPI_COMM_WORLD, &requests[n_requests++]);
#elif defined(GFM)
              MPI_Irecv(DustP + N_dust + offset_recv_dust[target] - offset_recv_star[NTask - 1] -
                        count_recv_star[NTask - 1], count_recv_dust[target] * sizeof(struct dust_particle_data), MPI_BYTE, target, TAG_DUSTDATA, MPI_COMM_WORLD, &requests[n_requests++]);
#else
              MPI_Irecv(DustP + N_dust + offset_recv_dust[target] - offset_recv_sph[NTask - 1] -
                        count_recv_sph[NTask - 1], count_recv_dust[target] * sizeof(struct dust_particle_data), MPI_BYTE, target, TAG_DUSTDATA, MPI_COMM_WORLD, &requests[n_requests++]);
#endif
            }
#endif
          if(count_recv[target] > 0)
            {
              MPI_Irecv(P + offset_recv[target], count_recv[target] * sizeof(struct particle_data), MPI_BYTE, target, TAG_PDATA, MPI_COMM_WORLD, &requests[n_requests++]);

              MPI_Irecv(Key + offset_recv[target], count_recv[target] * sizeof(peanokey), MPI_BYTE, target, TAG_KEY, MPI_COMM_WORLD, &requests[n_requests++]);
            }
        }
    }

  MPI_Barrier(MPI_COMM_WORLD);  /* not really necessary, but this will guarantee that all receives are
                                   posted before the sends, which helps the stability of MPI on
                                   bluegene, and perhaps some mpich1-clusters */


  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      target = ThisTask ^ ngrp;

      if(target < NTask)
        {
          if(count_sph[target] > 0)
            {
              MPI_Isend(partBuf + offset_sph[target], count_sph[target] * sizeof(struct particle_data), MPI_BYTE, target, TAG_PDATA_SPH, MPI_COMM_WORLD, &requests[n_requests++]);

              MPI_Isend(sphBuf + offset_sph[target], count_sph[target] * sizeof(struct sph_particle_data), MPI_BYTE, target, TAG_SPHDATA, MPI_COMM_WORLD, &requests[n_requests++]);

              MPI_Isend(keyBuf + offset_sph[target], count_sph[target] * sizeof(peanokey), MPI_BYTE, target, TAG_KEY_SPH, MPI_COMM_WORLD, &requests[n_requests++]);
            }

#ifdef TRACER_MC
          if(count_tracer[target] > 0)
            {
              MPI_Isend(tracerBuf + offset_tracer[target], count_tracer[target] * sizeof(struct tracer_linked_list), MPI_BYTE, target, TAG_TRACERDATA, MPI_COMM_WORLD, &requests[n_requests++]);
            }
#endif
#ifdef GFM
          if(count_star[target] > 0)
            {
              MPI_Isend(partBuf + offset_star[target], count_star[target] * sizeof(struct particle_data), MPI_BYTE, target, TAG_PDATA_STAR, MPI_COMM_WORLD, &requests[n_requests++]);

              MPI_Isend(starBuf + offset_star[target] - offset_sph[NTask - 1] - count_sph[NTask - 1],
                        count_star[target] * sizeof(struct star_particle_data), MPI_BYTE, target, TAG_STARDATA, MPI_COMM_WORLD, &requests[n_requests++]);

              MPI_Isend(keyBuf + offset_star[target], count_star[target] * sizeof(peanokey), MPI_BYTE, target, TAG_KEY_STAR, MPI_COMM_WORLD, &requests[n_requests++]);
            }
#endif
#ifdef BLACK_HOLES
          if(count_BHs[target] > 0)
            {
              MPI_Isend(partBuf + offset_BHs[target], count_BHs[target] * sizeof(struct particle_data), MPI_BYTE, target, TAG_PDATA_BH, MPI_COMM_WORLD, &requests[n_requests++]);

#ifdef GFM
              MPI_Isend(BHsBuf + offset_BHs[target] - offset_star[NTask - 1] - count_star[NTask - 1],
                        count_BHs[target] * sizeof(struct bh_particle_data), MPI_BYTE, target, TAG_BHDATA, MPI_COMM_WORLD, &requests[n_requests++]);
#else
              MPI_Isend(BHsBuf + offset_BHs[target] - offset_sph[NTask - 1] - count_sph[NTask - 1],
                        count_BHs[target] * sizeof(struct bh_particle_data), MPI_BYTE, target, TAG_BHDATA, MPI_COMM_WORLD, &requests[n_requests++]);
#endif
              MPI_Isend(keyBuf + offset_BHs[target], count_BHs[target] * sizeof(peanokey), MPI_BYTE, target, TAG_KEY_BH, MPI_COMM_WORLD, &requests[n_requests++]);
            }
#endif
#ifdef DUST_LIVE
          if(count_dust[target] > 0)
            {
              MPI_Isend(partBuf + offset_dust[target], count_dust[target] * sizeof(struct particle_data), MPI_BYTE, target, TAG_PDATA_DUST, MPI_COMM_WORLD, &requests[n_requests++]);

#ifdef BLACK_HOLES
              MPI_Isend(dustBuf + offset_dust[target] - offset_BHs[NTask - 1] - count_BHs[NTask - 1],
                        count_dust[target] * sizeof(struct dust_particle_data), MPI_BYTE, target, TAG_DUSTDATA, MPI_COMM_WORLD, &requests[n_requests++]);
#elif defined(GFM)
              MPI_Isend(dustBuf + offset_dust[target] - offset_star[NTask - 1] - count_star[NTask - 1],
                        count_dust[target] * sizeof(struct dust_particle_data), MPI_BYTE, target, TAG_DUSTDATA, MPI_COMM_WORLD, &requests[n_requests++]);
#else
              MPI_Isend(dustBuf + offset_dust[target] - offset_sph[NTask - 1] - count_sph[NTask - 1],
                        count_dust[target] * sizeof(struct dust_particle_data), MPI_BYTE, target, TAG_DUSTDATA, MPI_COMM_WORLD, &requests[n_requests++]);
#endif
              MPI_Isend(keyBuf + offset_dust[target], count_dust[target] * sizeof(peanokey), MPI_BYTE, target, TAG_KEY_DUST, MPI_COMM_WORLD, &requests[n_requests++]);
            }
#endif

          if(count[target] > 0)
            {
              MPI_Isend(partBuf + offset[target], count[target] * sizeof(struct particle_data), MPI_BYTE, target, TAG_PDATA, MPI_COMM_WORLD, &requests[n_requests++]);

              MPI_Isend(keyBuf + offset[target], count[target] * sizeof(peanokey), MPI_BYTE, target, TAG_KEY, MPI_COMM_WORLD, &requests[n_requests++]);
            }
        }
    }

  MPI_Waitall(n_requests, requests, MPI_STATUSES_IGNORE);
  myfree(requests);
#endif

#else /* begins block of myMPI_Alltoallv communications */

  myMPI_Alltoallv(partBuf, count_sph, offset_sph, P, count_recv_sph, offset_recv_sph, sizeof(struct particle_data), 0, MPI_COMM_WORLD);

  myMPI_Alltoallv(sphBuf, count_sph, offset_sph, SphP, count_recv_sph, offset_recv_sph, sizeof(struct sph_particle_data), 0, MPI_COMM_WORLD);

  myMPI_Alltoallv(keyBuf, count_sph, offset_sph, Key, count_recv_sph, offset_recv_sph, sizeof(peanokey), 0, MPI_COMM_WORLD);

  myMPI_Alltoallv(partBuf, count, offset, P, count_recv, offset_recv, sizeof(struct particle_data), 0, MPI_COMM_WORLD);

  myMPI_Alltoallv(keyBuf, count, offset, Key, count_recv, offset_recv, sizeof(peanokey), 0, MPI_COMM_WORLD);

#ifdef TRACER_MC
  myMPI_Alltoallv(tracerBuf, count_tracer, offset_tracer, tracerRecvBuf, count_recv_tracer, offset_recv_tracer, sizeof(struct tracer_linked_list), 0, MPI_COMM_WORLD);
#endif

#ifdef GFM
  myMPI_Alltoallv(partBuf, count_star, offset_star, P, count_recv_star, offset_recv_star, sizeof(struct particle_data), 0, MPI_COMM_WORLD);

  myMPI_Alltoallv(starBuf, count_star, starbuf_offset, StarP, count_recv_star, starbuf_recv_offset, sizeof(struct star_particle_data), 0, MPI_COMM_WORLD);

  myMPI_Alltoallv(keyBuf, count_star, offset_star, Key, count_recv_star, offset_recv_star, sizeof(peanokey), 0, MPI_COMM_WORLD);
#endif

#ifdef BLACK_HOLES
  myMPI_Alltoallv(partBuf, count_BHs, offset_BHs, P, count_recv_BHs, offset_recv_BHs, sizeof(struct particle_data), 0, MPI_COMM_WORLD);

  myMPI_Alltoallv(BHsBuf, count_BHs, bhbuf_offset, BHP, count_recv_BHs, bhbuf_recv_offset, sizeof(struct bh_particle_data), 0, MPI_COMM_WORLD);

  myMPI_Alltoallv(keyBuf, count_BHs, offset_BHs, Key, count_recv_BHs, offset_recv_BHs, sizeof(peanokey), 0, MPI_COMM_WORLD);
#endif

#ifdef DUST_LIVE
  myMPI_Alltoallv(partBuf, count_dust, offset_dust, P, count_recv_dust, offset_recv_dust, sizeof(struct particle_data), 0, MPI_COMM_WORLD);

  myMPI_Alltoallv(dustBuf, count_dust, dustbuf_offset, DustP, count_recv_dust, dustbuf_recv_offset, sizeof(struct dust_particle_data), 0, MPI_COMM_WORLD);

  myMPI_Alltoallv(keyBuf, count_dust, offset_dust, Key, count_recv_dust, offset_recv_dust, sizeof(peanokey), 0, MPI_COMM_WORLD);
#endif

#endif /* close block of myMPI_Alltoallv communications */


#ifdef GFM
  for(i = 0; i < count_get_star; i++)
    {
      StarP[N_star + i].PID = offset_recv_star[0] + i;
      P[offset_recv_star[0] + i].AuxDataID = N_star + i;
    }
#endif
#ifdef BLACK_HOLES
  for(i = 0; i < count_get_BHs; i++)
    {
      BHP[NumBHs + i].PID = offset_recv_BHs[0] + i;
      P[offset_recv_BHs[0] + i].AuxDataID = NumBHs + i;
    }
#endif
#ifdef DUST_LIVE
  for(i = 0; i < count_get_dust; i++)
    {
      DustP[N_dust + i].PID = offset_recv_dust[0] + i;
      P[offset_recv_dust[0] + i].AuxDataID = N_dust + i;
    }
#endif

  NumPart += count_get;
  NumGas += count_get_sph;

#ifdef GFM
  N_star += count_get_star;
#endif
#ifdef BLACK_HOLES
  NumBHs += count_get_BHs;
#endif
#ifdef DUST_LIVE
  N_dust += count_get_dust;
#endif

#ifdef TRACER_MC
  /* reattach imported tracers with imported parents */

  for(int i = 0; i < NumPart; i++)
    if(P[i].OriginTask != ThisTask)
      {
        if(P[i].NumberOfTracers)
          {
            int count = P[i].NumberOfTracers;
            int task = P[i].OriginTask;
            int off = offset_recv_tracer[task] + P[i].TracerHead;

            P[i].NumberOfTracers = 0;
            P[i].TracerHead = -1;

            for(int j = 0; j < count; j++)
              {
                int itr = get_free_tracer_slot();

                TracerLinkedList[itr] = tracerRecvBuf[off++];

                add_tracer_to_parent(i, itr);
              }
          }
        else
          {
            P[i].TracerHead = -1;
          }
      }
#endif /* TRACER_MC */


  myfree(keyBuf);
#ifdef DUST_LIVE
  myfree(dustBuf);
#endif
#ifdef BLACK_HOLES
  myfree(BHsBuf);
#endif
#ifdef GFM
  myfree(starBuf);
#endif
#ifdef TRACER_MC
  myfree(tracerRecvBuf);
  myfree(tracerBuf);
#endif
  myfree(sphBuf);
  myfree(partBuf);
#ifdef DUST_LIVE
  myfree(dustbuf_recv_offset);
  myfree(dustbuf_offset);
  myfree(offset_recv_dust);
  myfree(count_recv_dust);
  myfree(offset_dust);
  myfree(count_dust);
#endif
#ifdef BLACK_HOLES
  myfree(bhbuf_recv_offset);
  myfree(bhbuf_offset);
  myfree(offset_recv_BHs);
  myfree(count_recv_BHs);
  myfree(offset_BHs);
  myfree(count_BHs);
#endif
#ifdef GFM
  myfree(starbuf_recv_offset);
  myfree(starbuf_offset);
  myfree(offset_recv_star);
  myfree(count_recv_star);
  myfree(offset_star);
  myfree(count_star);
#endif
#ifdef TRACER_MC
  myfree(offset_recv_tracer);
  myfree(count_recv_tracer);
  myfree(offset_tracer);
  myfree(count_tracer);
#endif
  myfree(offset_recv_sph);
  myfree(offset_recv);
  myfree(count_recv_sph);
  myfree(count_recv);
  myfree(offset_sph);
  myfree(offset);
  myfree(count_sph);
  myfree(count);





  double t1 = second();
  mpi_printf("DOMAIN: exchange of %lld particles done. (took %g sec)\n", sumtogo, timediff(t0, t1));


#ifdef TRACER_MC_CHECKS
  MPI_Barrier(MPI_COMM_WORLD);
  for(int i = 0; i < NumPart; i++)
    {
      if(P[i].Type == 1 && (P[i].TracerHead >= 0 || P[i].NumberOfTracers > 0))
        terminate("i=%d P[i].Type=%d  P[i].TracerHead=%d  P[i].NumberOfTracers=%d", i, P[i].Type, P[i].TracerHead, P[i].NumberOfTracers);
    }
  check_tracer_lists();
#endif
}
