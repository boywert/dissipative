/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/fof/fof_distribute.c
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

#include "../allvars.h"
#include "../proto.h"
#include "../domain.h"
#include "fof.h"
#include "../subfind/subfind.h"

#ifdef FOF



/* This function redistributes the particles according to what is stored in
 * PS[].TargetTask, and PS[].TargetIndex.
 */
void fof_subfind_exchange(MPI_Comm Communicator)
{
  int nimport, nexport;
  int i, j, n, type, ngrp, target;
  int max_load, max_loadsph, load;
  struct particle_data *partBuf;
  struct subfind_data *subBuf;
  struct sph_particle_data *sphBuf;


  int CommThisTask, CommNTask;

  MPI_Comm_size(Communicator, &CommNTask);
  MPI_Comm_rank(Communicator, &CommThisTask);


#ifdef GFM
  int max_loadstar;
  int star_nimport, star_nexport;
  struct star_particle_data *starBuf;
  int *star_send_count = mymalloc_movable(&star_send_count, "star_send_count", sizeof(int) * CommNTask);
  int *star_recv_count = mymalloc_movable(&star_recv_count, "star_recv_count", sizeof(int) * CommNTask);
  int *star_send_offset = mymalloc_movable(&star_send_offset, "star_send_offset", sizeof(int) * CommNTask);
  int *star_recv_offset = mymalloc_movable(&star_recv_offset, "star_recv_offset", sizeof(int) * CommNTask);
#endif
#ifdef BLACK_HOLES
  int max_loadbh;
  int bh_nimport, bh_nexport;
  struct bh_particle_data *BHsBuf;
  int *bh_send_count = mymalloc_movable(&bh_send_count, "bh_send_count", sizeof(int) * CommNTask);
  int *bh_recv_count = mymalloc_movable(&bh_recv_count, "bh_recv_count", sizeof(int) * CommNTask);
  int *bh_send_offset = mymalloc_movable(&bh_send_offset, "bh_send_offset", sizeof(int) * CommNTask);
  int *bh_recv_offset = mymalloc_movable(&bh_recv_offset, "bh_recv_offset", sizeof(int) * CommNTask);
#endif
#ifdef TRACER_MC
  int tracer_nimport, tracer_nexport;
  struct tracer_linked_list *tracerSendBuf, *tracerRecvBuf;
  int *tracer_send_count = mymalloc_movable(&tracer_send_count, "tracer_send_count", sizeof(int) * CommNTask);
  int *tracer_recv_count = mymalloc_movable(&tracer_recv_count, "tracer_recv_count", sizeof(int) * CommNTask);
  int *tracer_send_offset = mymalloc_movable(&tracer_send_offset, "tracer_send_offset", sizeof(int) * CommNTask);
  int *tracer_recv_offset = mymalloc_movable(&tracer_recv_offset, "tracer_recv_offset", sizeof(int) * CommNTask);
#endif

  int old_AllMaxPart = All.MaxPart;
  int old_AllMaxPartSph = All.MaxPartSph;
#ifdef GFM
  int old_AllMaxPartStar = All.MaxPartStar;
#endif
#ifdef BLACK_HOLES
  int old_AllMaxPartBHs = All.MaxPartBHs;
#endif


  for(type = 0; type < NTYPES; type++)
    {
      size_t ExportSpace = 0.5 * (FreeBytes); /* we will try to grab at most half of the still available memory  */
      size_t PartSpace = sizeof(struct particle_data) + sizeof(struct subfind_data) + sizeof(struct sph_particle_data);
#ifdef GFM
      PartSpace += sizeof(struct star_particle_data);
#endif
#ifdef BLACK_HOLES
      PartSpace += sizeof(struct bh_particle_data);
#endif
      if(PartSpace > ExportSpace)
	terminate("seems like we have insufficient storage, PartSpace=%lld ExportSpace=%lld", (long long) PartSpace, (long long) ExportSpace);

      int glob_flag = 0;

      do
	{
	  for(n = 0; n < CommNTask; n++)
	    {
	      Send_count[n] = 0;
#ifdef GFM
	      star_send_count[n] = 0;
#endif
#ifdef BLACK_HOLES
	      bh_send_count[n] = 0;
#endif
#ifdef TRACER_MC
	      tracer_send_count[n] = 0;
#endif
	    }

	  ptrdiff_t AvailableSpace = ExportSpace; /* this must be a type that can become negative */

	  for(n = 0; n < NumPart; n++)
	    {
	      if(AvailableSpace < 0)
		break;

	      if(P[n].Type == type && PS[n].TargetTask != CommThisTask)
		{
		  target = PS[n].TargetTask;

		  if(target < 0 || target >= CommNTask)
		    terminate("n=%d targettask=%d", n, target);

		  AvailableSpace -= PartSpace;

		  Send_count[target]++;
#ifdef GFM
		  if(P[n].Type == 4)
		    star_send_count[target]++;
#endif
#ifdef BLACK_HOLES
		  if(P[n].Type == 5)
		    bh_send_count[target]++;
#endif
#ifdef TRACER_MC
		  int next = P[n].TracerHead;
		  while(next >= 0)
		    {
		      tracer_send_count[target]++;
		      next = TracerLinkedList[next].Next;
		    }
#endif
		}
	    }

	  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, Communicator);
#ifdef GFM
	  MPI_Alltoall(star_send_count, 1, MPI_INT, star_recv_count, 1, MPI_INT, Communicator);
#endif
#ifdef BLACK_HOLES
	  MPI_Alltoall(bh_send_count, 1, MPI_INT, bh_recv_count, 1, MPI_INT, Communicator);
#endif
#ifdef TRACER_MC
	  MPI_Alltoall(tracer_send_count, 1, MPI_INT, tracer_recv_count, 1, MPI_INT, Communicator);
#endif

	  for(j = 0, nimport = 0, nexport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < CommNTask; j++)
	    {
	      nexport += Send_count[j];
	      nimport += Recv_count[j];

	      if(j > 0)
		{
		  Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
		  Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
		}
	    }

#ifdef GFM
	  for(j = 0, star_nimport = 0, star_nexport = 0, star_recv_offset[0] = 0, star_send_offset[0] = 0; j < CommNTask; j++)
	    {
	      star_nexport += star_send_count[j];
	      star_nimport += star_recv_count[j];

	      if(j > 0)
		{
		  star_send_offset[j] = star_send_offset[j - 1] + star_send_count[j - 1];
		  star_recv_offset[j] = star_recv_offset[j - 1] + star_recv_count[j - 1];
		}
	    }
#endif
#ifdef BLACK_HOLES
	  for(j = 0, bh_nimport = 0, bh_nexport = 0, bh_recv_offset[0] = 0, bh_send_offset[0] = 0; j < CommNTask; j++)
	    {
	      bh_nexport += bh_send_count[j];
	      bh_nimport += bh_recv_count[j];

	      if(j > 0)
		{
		  bh_send_offset[j] = bh_send_offset[j - 1] + bh_send_count[j - 1];
		  bh_recv_offset[j] = bh_recv_offset[j - 1] + bh_recv_count[j - 1];
		}
	    }
#endif
#ifdef TRACER_MC
	  for(j = 0, tracer_nimport = 0, tracer_nexport = 0, tracer_recv_offset[0] = 0, tracer_send_offset[0] = 0; j < CommNTask; j++)
	    {
	      tracer_nexport += tracer_send_count[j];
	      tracer_nimport += tracer_recv_count[j];

	      if(j > 0)
		{
		  tracer_send_offset[j] = tracer_send_offset[j - 1] + tracer_send_count[j - 1];
		  tracer_recv_offset[j] = tracer_recv_offset[j - 1] + tracer_recv_count[j - 1];
		}
	    }
#endif

	  /* for resize */
	  load = (NumPart + nimport - nexport);
	  MPI_Allreduce(&load, &max_load, 1, MPI_INT, MPI_MAX, Communicator);

	  if(type == 0)
	    {
	      load = (NumGas + nimport - nexport);
	      MPI_Allreduce(&load, &max_loadsph, 1, MPI_INT, MPI_MAX, Communicator);
	    }

#ifdef GFM
	  if(type == 4)
	    {
	      load = (N_star + star_nimport - star_nexport);
	      MPI_Allreduce(&load, &max_loadstar, 1, MPI_INT, MPI_MAX, Communicator);
	    }
#endif

#ifdef BLACK_HOLES
	  if(type == 5)
	    {
	      load = (NumBHs + bh_nimport - bh_nexport);
	      MPI_Allreduce(&load, &max_loadbh, 1, MPI_INT, MPI_MAX, Communicator);
	    }
#endif

#ifdef TRACER_MC
	  load = (N_tracer + tracer_nimport - tracer_nexport);
	  int max_loadtracer;
	  MPI_Allreduce(&load, &max_loadtracer, 1, MPI_INT, MPI_MAX, Communicator);
#endif

	  partBuf = (struct particle_data *) mymalloc_movable(&partBuf, "partBuf", nexport * sizeof(struct particle_data));
	  subBuf = (struct subfind_data *) mymalloc_movable(&subBuf, "subBuf", nexport * sizeof(struct subfind_data));
	  if(type == 0)
	    sphBuf = (struct sph_particle_data *) mymalloc_movable(&sphBuf, "sphBuf", nexport * sizeof(struct sph_particle_data));
#ifdef GFM
	  if(type == 4)
	    starBuf = (struct star_particle_data *) mymalloc_movable(&starBuf, "starBuf", star_nexport * sizeof(struct star_particle_data));
#endif
#ifdef BLACK_HOLES
	  if(type == 5)
	    BHsBuf = (struct bh_particle_data *) mymalloc_movable(&BHsBuf, "BHsBuf", bh_nexport * sizeof(struct bh_particle_data));
#endif
#ifdef TRACER_MC
	  tracerSendBuf = (struct tracer_linked_list *) mymalloc_movable(&tracerSendBuf, "tracerSendBuf",
	      tracer_nexport * sizeof(struct tracer_linked_list));
	  tracerRecvBuf = (struct tracer_linked_list *) mymalloc_movable(&tracerRecvBuf, "tracerRecvBuf",
	      tracer_nimport * sizeof(struct tracer_linked_list));
#endif

	  for(i = 0; i < CommNTask; i++)
	    {
	      Send_count[i] = 0;
#ifdef GFM
	      star_send_count[i] = 0;
#endif
#ifdef BLACK_HOLES
	      bh_send_count[i] = 0;
#endif
#ifdef TRACER_MC
	      tracer_send_count[i] = 0;
#endif
	    }


#ifdef TRACER_MC
	  for(n = 0; n < NumPart; n++)
	    P[n].OriginTask = CommThisTask;
#endif

	  AvailableSpace = ExportSpace; /* this must be allowed to become negative */

	  int nstay = 0;
          int delta_numpart = 0;
          int delta_numgas = 0;

	  for(n = 0; n < NumPart; n++)
	    {
	      if(AvailableSpace < 0)
		break;

 	      if(P[n].Type == type && PS[n].TargetTask != CommThisTask)
		{
		  target = PS[n].TargetTask;

		  AvailableSpace -= PartSpace;

		  partBuf[Send_offset[target] + Send_count[target]] = P[n];
		  subBuf[Send_offset[target] + Send_count[target]] = PS[n];

		  if(P[n].Type == 0)
		    {
		      sphBuf[Send_offset[target] + Send_count[target]] = SphP[n];
		      delta_numgas++;
		    }

#ifdef TRACER_MC
		  partBuf[Send_offset[target] + Send_count[target]].TracerHead = tracer_send_count[target];
#endif

#ifdef GFM
		  if(P[n].Type == 4)
		    {
		      if(N_star < 0)
			terminate("n=%d type=%d N_star=%d NumPart=%d\n", n, type, N_star, NumPart);

		      starBuf[star_send_offset[target] + star_send_count[target]++] = StarP[P[n].AuxDataID];
		      StarP[P[n].AuxDataID] = StarP[N_star - 1];
		      P[StarP[N_star - 1].PID].AuxDataID = P[n].AuxDataID;
		      N_star--;
		    }
#endif
#ifdef BLACK_HOLES
		  if(P[n].Type == 5)
		    {
		      BHsBuf[bh_send_offset[target] + bh_send_count[target]++] = BHP[P[n].AuxDataID];
		      BHP[P[n].AuxDataID] = BHP[NumBHs - 1];
		      P[BHP[NumBHs - 1].PID].AuxDataID = P[n].AuxDataID;
		      NumBHs--;
		    }
#endif

#ifdef TRACER_MC
		  int next = P[n].TracerHead;

		  while(next >= 0) /* fill tracerSendBuf export buffer and tag entries for removal */
		    {
		      int q = tracer_send_offset[target] + tracer_send_count[target]++;
		      tracerSendBuf[q] = TracerLinkedList[next];
		      remove_tracer_from_parent(n, next);
		      release_tracer_slot(next);

		      /* we must reform the Next, Prev connectivity on the remote task so we override them here,
			   use next=-1 to mark the end of the tracer chain for this parent,
			  and use prev for the Type of the parent */

		      next = TracerLinkedList[next].Next;
		    }
#endif
		  Send_count[target]++;
		  delta_numpart++;
		}
 	      else
 	        {
 	          if(nstay != n)
 	            {
 	              /* now move P[n] to P[nstay] */

 	              P[nstay] = P[n];
 	              PS[nstay] = PS[n];

		      if(P[nstay].Type == 0)
			SphP[nstay] = SphP[n];

#ifdef GFM
		      if(P[nstay].Type == 4)
			StarP[P[nstay].AuxDataID].PID = nstay;
#endif
#ifdef BLACK_HOLES
		      if(P[nstay].Type == 5)
			BHP[P[nstay].AuxDataID].PID = nstay;
#endif
		    }

		  nstay++;
 	        }
	    }

	  if(delta_numgas > 0)
	    if(delta_numpart != delta_numgas)
	      terminate("delta_numpart=%d != delta_numgas=%d", delta_numpart, delta_numgas);


	  /* now close gap (if present) */
	  memmove(P + nstay, P + nstay + delta_numpart, (NumPart - (nstay + delta_numpart)) * sizeof(struct particle_data));
	  memmove(PS + nstay, PS + nstay + delta_numpart, (NumPart - (nstay + delta_numpart)) * sizeof(struct subfind_data));

	  if(delta_numgas > 0)
	    if(NumGas - (nstay + delta_numgas) > 0)
	      memmove(SphP + nstay, SphP + nstay + delta_numpart, (NumGas - (nstay + delta_numgas)) * sizeof(struct sph_particle_data));
	    

          NumPart -= delta_numpart;
          NumGas -= delta_numgas;

#ifdef GFM
	  for(i = nstay; i < NumPart; i++)
	    if(P[i].Type == 4)
	      StarP[P[i].AuxDataID].PID = i;
#endif
#ifdef BLACK_HOLES
          for(i = nstay; i < NumPart; i++)
            if(P[i].Type == 5)
              BHP[P[i].AuxDataID].PID = i;
#endif

	  /* do resize, but only increase arrays!! (otherwise data in ActiveParticleList etc. gets lost */
	  if(max_load > (1.0 - ALLOC_TOLERANCE) * All.MaxPart)
	    {
	      All.MaxPart = max_load / (1.0 - 2 * ALLOC_TOLERANCE);
	      reallocate_memory_maxpart();
	      PS = (struct subfind_data *) myrealloc_movable(PS, All.MaxPart * sizeof(struct subfind_data));
	    }

	  if(type == 0)
	    {
	      if(max_loadsph > (1.0 - ALLOC_TOLERANCE) * All.MaxPartSph)
		{
		  All.MaxPartSph = max_loadsph / (1.0 - 2 * ALLOC_TOLERANCE);
		  reallocate_memory_maxpartsph();
		}
	    }

#ifdef GFM
	  if(type == 4)
	    {
	      if(max_loadstar >= All.MaxPartStar - 0.5 * ALLOC_STARBH_ROOM * (All.TotNumGas / NTask))
		{
		  All.MaxPartStar = max_loadstar + ALLOC_STARBH_ROOM * (All.TotNumGas / NTask);
		  reallocate_memory_maxpartstar();
		}
	    }
#endif

#ifdef BLACK_HOLES
	  if(type == 5)
	    {
	      if(max_loadbh >= All.MaxPartBHs - 0.5 * ALLOC_STARBH_ROOM * (All.TotNumGas / NTask))
		{
		  All.MaxPartBHs = max_loadbh + ALLOC_STARBH_ROOM * (All.TotNumGas / NTask);
		  reallocate_memory_maxpartBHs();
		}
	    }
#endif

#ifdef TRACER_MC
	  if(max_loadtracer >= All.MaxPartTracer - sqrt(All.MaxPartTracer))
	    {
	      int delta = max_loadtracer - All.MaxPartTracer + 2 * sqrt(All.MaxPartTracer);
	      reallocate_memory_maxparttracer(delta);
	    }
#endif

	  /* create a gap behind the existing gas particles where we will insert the incoming particles */
	  memmove(P + NumGas + nimport, P + NumGas, (NumPart - NumGas) * sizeof(struct particle_data));
	  memmove(PS + NumGas + nimport, PS + NumGas, (NumPart - NumGas) * sizeof(struct subfind_data));

#ifdef GFM
	  for(i = 0; i < N_star; i++)
	    if(StarP[i].PID >= NumGas)
	      StarP[i].PID += nimport;
#endif
#ifdef BLACK_HOLES
	  for(i = 0; i < NumBHs; i++)
	    if(BHP[i].PID >= NumGas)
	      BHP[i].PID += nimport;
#endif

	  for(i = 0; i < CommNTask; i++)
	    Recv_offset[i] += NumGas;

	  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	    {
	      target = CommThisTask ^ ngrp;

	      if(target < CommNTask)
		{
		  if(Send_count[target] > 0 || Recv_count[target] > 0)
		    {
		      MPI_Sendrecv(partBuf + Send_offset[target], Send_count[target] * sizeof(struct particle_data), MPI_BYTE, target, TAG_PDATA,
			  P + Recv_offset[target], Recv_count[target] * sizeof(struct particle_data), MPI_BYTE, target, TAG_PDATA, Communicator,
			  MPI_STATUS_IGNORE);

		      MPI_Sendrecv(subBuf + Send_offset[target], Send_count[target] * sizeof(struct subfind_data),
		      MPI_BYTE, target, TAG_KEY, PS + Recv_offset[target], Recv_count[target] * sizeof(struct subfind_data), MPI_BYTE, target,
			  TAG_KEY, Communicator, MPI_STATUS_IGNORE);

		      if(type == 0)
			MPI_Sendrecv(sphBuf + Send_offset[target], Send_count[target] * sizeof(struct sph_particle_data), MPI_BYTE, target,
			TAG_SPHDATA, SphP + Recv_offset[target], Recv_count[target] * sizeof(struct sph_particle_data), MPI_BYTE, target, TAG_SPHDATA,
			    Communicator, MPI_STATUS_IGNORE);
#ifdef GFM
		      if(type == 4)
			MPI_Sendrecv(starBuf + star_send_offset[target], star_send_count[target] * sizeof(struct star_particle_data), MPI_BYTE,
			    target,
			    TAG_STARDATA, StarP + star_recv_offset[target] + N_star, star_recv_count[target] * sizeof(struct star_particle_data),
			    MPI_BYTE, target, TAG_STARDATA, Communicator, MPI_STATUS_IGNORE);
#endif
#ifdef BLACK_HOLES
		      if(type == 5)
			MPI_Sendrecv(BHsBuf + bh_send_offset[target], bh_send_count[target] * sizeof(struct bh_particle_data), MPI_BYTE, target,
			TAG_BHDATA, BHP + bh_recv_offset[target] + NumBHs, bh_recv_count[target] * sizeof(struct bh_particle_data), MPI_BYTE, target,
			    TAG_BHDATA, Communicator, MPI_STATUS_IGNORE);
#endif

#ifdef TRACER_MC
		      MPI_Sendrecv(tracerSendBuf + tracer_send_offset[target], tracer_send_count[target] * sizeof(struct tracer_linked_list),
		      MPI_BYTE, target, TAG_TRACERDATA, tracerRecvBuf + tracer_recv_offset[target],
			  tracer_recv_count[target] * sizeof(struct tracer_linked_list), MPI_BYTE, target, TAG_TRACERDATA, Communicator,
			  MPI_STATUS_IGNORE);
#endif
		    }
		}
	    }

#ifdef GFM
	  if(type == 4)
	    {
	      int i_tmp = 0;

	      for(i = 0; i < nimport; i++)
		if(P[NumGas + i].Type == 4)
		  {
		    StarP[N_star + i_tmp].PID = NumGas + i;
		    P[NumGas + i].AuxDataID = N_star + i_tmp;
		    i_tmp++;
		  }

	      N_star += star_nimport;
	    }
#endif
#ifdef BLACK_HOLES
	  if(type == 5)
	    {
	      int i_tmp = 0;

	      for(i = 0; i < nimport; i++)
		if(P[NumGas + i].Type == 5)
		  {
		    BHP[NumBHs + i_tmp].PID = NumGas + i;
		    P[NumGas + i].AuxDataID = NumBHs + i_tmp;
		    i_tmp++;
		  }

	      NumBHs += bh_nimport;
	    }
#endif

	  if(type == 0)
	    NumGas += nimport;

	  NumPart += nimport;

#ifdef TRACER_MC
	  /* reattach imported tracers with imported parents */

	  for(int i = 0; i < NumPart; i++)
	    if(P[i].OriginTask != CommThisTask)
	      {
		if(P[i].NumberOfTracers)
		  {
		    int count = P[i].NumberOfTracers;
		    int task = P[i].OriginTask;
		    int off = tracer_recv_offset[task] + P[i].TracerHead;

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

#ifdef TRACER_MC
	  myfree(tracerRecvBuf);
	  myfree(tracerSendBuf);
#endif
#ifdef BLACK_HOLES
	  if(type == 5)
	    myfree_movable(BHsBuf);
#endif
#ifdef GFM
	  if(type == 4)
	    myfree_movable(starBuf);
#endif
	  if(type == 0)
	    myfree_movable(sphBuf);

	  myfree_movable(subBuf);
	  myfree_movable(partBuf);


	  int loc_flag = 0;
	  if(AvailableSpace < 0)
	    loc_flag = 1;

	  MPI_Allreduce(&loc_flag, &glob_flag, 1, MPI_INT, MPI_SUM, Communicator);
	  if(glob_flag > 0 && CommThisTask == 0)
	    {
	      printf("FOF-DISTRIBUTE: Need to cycle in particle exchange due to memory shortage. type=%d glob_flag=%d ThisTask=%d CommThisTask=%d   PartSpace=%lld  ExportSpace=%lld\n",
		      type, glob_flag, ThisTask, CommThisTask, (long long) PartSpace, (long long) ExportSpace);
	      fflush(stdout);
	    }
	}
      while(glob_flag);
    }



  /* if there was a temporary memory shortage during the exchange, we may had to increase the maximum allocations. Go back to smaller values again if possible */

  load = NumPart;
  MPI_Allreduce(&load, &max_load, 1, MPI_INT, MPI_MAX, Communicator);
  max_load = max_load / (1.0 - 2 * ALLOC_TOLERANCE);
  if(max_load < old_AllMaxPart)
    max_load = old_AllMaxPart;
  if(max_load != All.MaxPart)
    {
      All.MaxPart = max_load;
      reallocate_memory_maxpart();
      PS = (struct subfind_data *) myrealloc_movable(PS, All.MaxPart * sizeof(struct subfind_data));
    }

  load = NumGas;
  MPI_Allreduce(&load, &max_loadsph, 1, MPI_INT, MPI_MAX, Communicator);
  max_loadsph = max_loadsph / (1.0 - 2 * ALLOC_TOLERANCE);
  if(max_loadsph < old_AllMaxPartSph)
    max_loadsph = old_AllMaxPartSph;
  if(max_loadsph != All.MaxPartSph)
    {
      All.MaxPartSph = max_loadsph;
      reallocate_memory_maxpartsph();
    }

#ifdef GFM
  load = N_star;
  MPI_Allreduce(&load, &max_loadstar, 1, MPI_INT, MPI_MAX, Communicator);
  max_loadstar += ALLOC_STARBH_ROOM * (All.TotNumGas / NTask);
  if(max_loadstar < old_AllMaxPartStar)
    max_loadstar = old_AllMaxPartStar;
  if(max_loadstar != All.MaxPartStar)
    {
      All.MaxPartStar = max_loadstar;
      reallocate_memory_maxpartsph();
    }
#endif

#ifdef BLACK_HOLES
  load = NumBHs;
  MPI_Allreduce(&load, &max_loadbh, 1, MPI_INT, MPI_MAX, Communicator);
  max_loadbh += ALLOC_STARBH_ROOM * (All.TotNumGas / NTask);
  if(max_loadbh < old_AllMaxPartBHs)
    max_loadbh = old_AllMaxPartBHs;
  if(max_loadbh != All.MaxPartBHs)
    {
      All.MaxPartBHs = max_loadbh;
      reallocate_memory_maxpartBHs();
    }
#endif


  /* finally, let's also address the desired local order according to PS[].TargetIndex */

  struct fof_local_sort_data *mp;
  int *Id;

  if(NumGas)
    {
      mp = (struct fof_local_sort_data *) mymalloc("mp", sizeof(struct fof_local_sort_data) * NumGas);
      Id = (int *) mymalloc("Id", sizeof(int) * NumGas);

      for(i = 0; i < NumGas; i++)
        {
          mp[i].index = i;
          mp[i].targetindex = PS[i].TargetIndex;
        }

      qsort(mp, NumGas, sizeof(struct fof_local_sort_data), fof_compare_local_sort_data_targetindex);

      for(i = 0; i < NumGas; i++)
        Id[mp[i].index] = i;

      reorder_gas(Id);

      for(i = 0; i < NumGas; i++)
        Id[mp[i].index] = i;

      fof_reorder_PS(Id, 0, NumGas);

      myfree(Id);
      myfree(mp);
    }

  if(NumPart - NumGas > 0)
    {
      mp = (struct fof_local_sort_data *) mymalloc("mp", sizeof(struct fof_local_sort_data) * (NumPart - NumGas));
      mp -= NumGas;

      Id = (int *) mymalloc("Id", sizeof(int) * (NumPart - NumGas));
      Id -= NumGas;

      for(i = NumGas; i < NumPart; i++)
        {
          mp[i].index = i;
          mp[i].targetindex = PS[i].TargetIndex;
        }

      qsort(mp + NumGas, NumPart - NumGas, sizeof(struct fof_local_sort_data), fof_compare_local_sort_data_targetindex);

      for(i = NumGas; i < NumPart; i++)
        Id[mp[i].index] = i;

      reorder_particles(Id);

      for(i = NumGas; i < NumPart; i++)
        Id[mp[i].index] = i;

      fof_reorder_PS(Id, NumGas, NumPart);


      Id += NumGas;
      myfree(Id);
      mp += NumGas;
      myfree(mp);
    }

#ifdef TRACER_MC
  myfree_movable(tracer_recv_offset);
  myfree_movable(tracer_send_offset);
  myfree_movable(tracer_recv_count);
  myfree_movable(tracer_send_count);
#endif

#ifdef BLACK_HOLES
  myfree_movable(bh_recv_offset);
  myfree_movable(bh_send_offset);
  myfree_movable(bh_recv_count);
  myfree_movable(bh_send_count);
#endif
#ifdef GFM
  myfree_movable(star_recv_offset);
  myfree_movable(star_send_offset);
  myfree_movable(star_recv_count);
  myfree_movable(star_send_count);
#endif
}


void fof_reorder_PS(int *Id, int Nstart, int N)
{
  int i;
  struct subfind_data PSsave, PSsource;
  int idsource, idsave, dest;

  for(i = Nstart; i < N; i++)
    {
      if(Id[i] != i)
        {
          PSsource = PS[i];

          idsource = Id[i];
          dest = Id[i];

          do
            {
              PSsave = PS[dest];
              idsave = Id[dest];

              PS[dest] = PSsource;
              Id[dest] = idsource;

              if(dest == i)
                break;

              PSsource = PSsave;
              idsource = idsave;

              dest = idsource;
            }
          while(1);
        }
    }
}


#endif
