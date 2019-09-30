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
#include <math.h>
#include <gsl/gsl_math.h>

#include "../allvars.h"
#include "../proto.h"


/*! \file mpi_util.c
 *  \brief a wrapper around MPI_Alltoallv that can deal with data in individual sends that are very big
 */


void myMPI_Alltoallv(void *sendb, size_t * sendcounts, size_t * sdispls, void *recvb, size_t * recvcounts, size_t * rdispls, int len, int big_flag, MPI_Comm comm)
{
  char *sendbuf = (char *) sendb;
  char *recvbuf = (char *) recvb;

  if(big_flag == 0)
    {
      int ntask;
      MPI_Comm_size(comm, &ntask);

      int *scount = (int *) mymalloc("scount", ntask * sizeof(int));
      int *rcount = (int *) mymalloc("rcount", ntask * sizeof(int));
      int *soff = (int *) mymalloc("soff", ntask * sizeof(int));
      int *roff = (int *) mymalloc("roff", ntask * sizeof(int));

      for(int i = 0; i < ntask; i++)
        {
          scount[i] = sendcounts[i] * len;
          rcount[i] = recvcounts[i] * len;
          soff[i] = sdispls[i] * len;
          roff[i] = rdispls[i] * len;
        }

      MPI_Alltoallv(sendbuf, scount, soff, MPI_BYTE, recvbuf, rcount, roff, MPI_BYTE, comm);

      myfree(roff);
      myfree(soff);
      myfree(rcount);
      myfree(scount);
    }
  else
    {
      /* here we definitely have some large messages. We default to the
       * pair-wise protocoll, which should be most robust anyway.
       */

      int ntask, thistask;
      MPI_Comm_size(comm, &ntask);
      MPI_Comm_rank(comm, &thistask);

      for(int ngrp = 0; ngrp < (1 << PTask); ngrp++)
        {
          int target = thistask ^ ngrp;

          if(target < ntask)
            {
              if(sendcounts[target] > 0 || recvcounts[target] > 0)
                myMPI_Sendrecv(sendbuf + sdispls[target] * len, sendcounts[target] * len, MPI_BYTE, target, TAG_PDATA + ngrp,
                               recvbuf + rdispls[target] * len, recvcounts[target] * len, MPI_BYTE, target, TAG_PDATA + ngrp, comm, MPI_STATUS_IGNORE);
            }
        }
    }
}



void my_int_MPI_Alltoallv(void *sendb, int *sendcounts, int *sdispls, void *recvb, int *recvcounts, int *rdispls, int len, int big_flag, MPI_Comm comm)
{
  char *sendbuf = (char *) sendb;
  char *recvbuf = (char *) recvb;

  if(big_flag == 0)
    {
      int i, ntask;
      MPI_Comm_size(comm, &ntask);

      int *scount = (int *) mymalloc("scount", ntask * sizeof(int));
      int *rcount = (int *) mymalloc("rcount", ntask * sizeof(int));
      int *soff = (int *) mymalloc("soff", ntask * sizeof(int));
      int *roff = (int *) mymalloc("roff", ntask * sizeof(int));

      for(i = 0; i < ntask; i++)
        {
          scount[i] = sendcounts[i] * len;
          rcount[i] = recvcounts[i] * len;
          soff[i] = sdispls[i] * len;
          roff[i] = rdispls[i] * len;
        }

      MPI_Alltoallv(sendbuf, scount, soff, MPI_BYTE, recvbuf, rcount, roff, MPI_BYTE, comm);

      myfree(roff);
      myfree(soff);
      myfree(rcount);
      myfree(scount);
    }
  else
    {
      /* here we definitely have some large messages. We default to the
       * pair-wise protocoll, which should be most robust anyway.
       */

      int ntask, thistask;
      MPI_Comm_size(comm, &ntask);
      MPI_Comm_rank(comm, &thistask);

      for(int ngrp = 0; ngrp < (1 << PTask); ngrp++)
        {
          int target = thistask ^ ngrp;

          if(target < ntask)
            {
              if(sendcounts[target] > 0 || recvcounts[target] > 0)
                myMPI_Sendrecv(sendbuf + sdispls[target] * len, sendcounts[target] * len, MPI_BYTE, target, TAG_PDATA + ngrp,
                               recvbuf + rdispls[target] * len, recvcounts[target] * len, MPI_BYTE, target, TAG_PDATA + ngrp, comm, MPI_STATUS_IGNORE);
            }
        }
    }
}







/* this implements a sparse pattern where only MPI_COMM_WORLDunications of data elements occur that are non-zero, followed by
 * an asynchronous barrier.
 * It is means to act as a potential replacement for the frequent construct:
 * MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);
 */
void mySparse_MPI_Alltoall_OneInt(int *send_count, int *recv_count)
{
  static int MyTagOffset = 0;
  MyTagOffset++;
  if(MyTagOffset > 30000)
    MyTagOffset = 0;

  /* variables for non-blacking barrier */
  int nLevels = my_fls(NTask - 1);
  int received_levels = 0, sent_levels = 0;
  int *stagelist = (int *) mymalloc("stagelist", nLevels * sizeof(int));
  for(int j = 0; j < nLevels; j++)
    stagelist[j] = j;

  MPI_Request *level_requests = (MPI_Request *) mymalloc("level_requests", nLevels * sizeof(MPI_Request));
  MPI_Request *requests = (MPI_Request *) mymalloc("requests", NTask * sizeof(MPI_Request));
  int n_requests = 0;

  /* post our send requests */
  for(int ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int j = ThisTask ^ ngrp;

      if(j < NTask)
        if(send_count[j] != 0)
          MPI_Issend(&send_count[j], 1, MPI_INT, j, TAG_N + MyTagOffset, MPI_COMM_WORLD, &requests[n_requests++]);
    }

  /* make the default result equal to zero */
  memset(recv_count, 0, NTask * sizeof(int));

  int barrier_active = 0;

  while(1)
    {
      int flag;
      MPI_Status status;

      MPI_Iprobe(MPI_ANY_SOURCE, TAG_N + MyTagOffset, MPI_COMM_WORLD, &flag, &status);

      if(flag)
        {
          int source = status.MPI_SOURCE;
          int tag = status.MPI_TAG;

          MPI_Recv(&recv_count[source], 1, MPI_INT, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

      MPI_Iprobe(MPI_ANY_SOURCE, TAG_BARRIER + MyTagOffset, MPI_COMM_WORLD, &flag, &status);

      if(flag)
        {
          int source = status.MPI_SOURCE;
          int tag = status.MPI_TAG;

          int stage;
          MPI_Recv(&stage, 1, MPI_INT, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          received_levels |= (1 << stage);
        }

      if(barrier_active)
        {
          for(int stage = 0; stage < nLevels; stage++)
            if(!(sent_levels & (1 << stage)))
              {
                int mask = ((1 << stage) - 1);

                if((mask & received_levels) == mask)
                  {
                    sent_levels |= (1 << stage);

                    int target = (ThisTask + (1 << stage)) % NTask;

                    MPI_Issend(&stagelist[stage], 1, MPI_INT, target, TAG_BARRIER + MyTagOffset, MPI_COMM_WORLD, &level_requests[stage]);
                  }
              }
        }
      else
        {
          MPI_Testall(n_requests, requests, &flag, MPI_STATUSES_IGNORE);

          if(flag)
            {
              barrier_active = 1;
              n_requests = 0;
            }
        }

      if(received_levels == ((1 << nLevels) - 1) && sent_levels == ((1 << nLevels) - 1))
        break;
    }

  MPI_Waitall(nLevels, level_requests, MPI_STATUSES_IGNORE);  /* as we are going to free stagelist */

  myfree(requests);
  myfree(level_requests);
  myfree(stagelist);
}


#ifdef PERFORMANCE_TEST_SPARSE_MPI_ALLTOALL

void test_mpi_alltoall_performance(void)
{
  Send_count = (int *) mymalloc("Send_count", sizeof(int) * NTask);
  Recv_count = (int *) mymalloc("Recv_count", sizeof(int) * NTask);

  for(int rep = 0; rep < 3; rep++)
    {
      if(rep == 0)
	{
	  for(int i = 0; i < NTask; i++)
	    Send_count[i] = i + 1;
	}
      else if(rep == 1)
	{
	  for(int i = 0; i < NTask; i++)
	    Send_count[i] = 0;

	  for(int i = 0; i < NTask / 8; i++)
	    {
	      int j = get_random_number() * NTask;
	      Send_count[j] = 10;
	    }
	}
      else if(rep == 2)
	{
	  for(int i = 0; i < NTask; i++)
	    Send_count[i] = 0;

	  if(get_random_number() < 0.2)
	    {
	      for(int i = 0; i < 10; i++)
		{
		  int j = get_random_number() * NTask;
		  Send_count[j] = 10;
		}
	    }
	}

      mpi_printf("\n");

      MPI_Barrier(MPI_COMM_WORLD);

      double t0 = second();

      for(int i = 0; i < 1000; i++)
	MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

      double t1 = second();

      mpi_printf("REP=%d:  Ordinary MPI_Alltoall took %g sec\n", rep, timediff(t0, t1));

      MPI_Barrier(MPI_COMM_WORLD);


      MPI_Barrier(MPI_COMM_WORLD);

      t0 = second();

      for(int i = 0; i < 1000; i++)
	mySparse_MPI_Alltoall_OneInt(Send_count, Recv_count);

      t1 = second();

      mpi_printf("REP=%d:  Sparse MPI_Alltoall took %g sec\n", rep, timediff(t0, t1));

      MPI_Barrier(MPI_COMM_WORLD);

      mpi_printf("\n");
    }

  endrun();
}


#endif





