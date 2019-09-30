/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/mpi_utils/hypercube_allgatherv.c
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

#ifdef MPI_HYPERCUBE_ALLGATHERV

#define TAG 100

int MPI_hypercube_Allgatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int *recvcount, int *displs, MPI_Datatype recvtype, MPI_Comm comm)
{
  int ntask, thistask, ptask, ngrp, size_sendtype, size_recvtype;
  MPI_Status status;

  MPI_Comm_rank(comm, &thistask);
  MPI_Comm_size(comm, &ntask);

  MPI_Type_size(sendtype, &size_sendtype);
  MPI_Type_size(recvtype, &size_recvtype);

  for(ptask = 0; ntask > (1 << ptask); ptask++);

  for(ngrp = 1; ngrp < (1 << ptask); ngrp++)
    {
      int recvtask = thistask ^ ngrp;

      if(recvtask < ntask)
        MPI_Sendrecv(sendbuf, sendcount, sendtype, recvtask, TAG, recvbuf + displs[recvtask] * size_recvtype, recvcount[recvtask], recvtype, recvtask, TAG, comm, &status);
    }

  if(sendbuf != recvbuf + displs[thistask] * size_recvtype)
    memcpy(recvbuf + displs[thistask] * size_recvtype, sendbuf, sendcount * size_sendtype);

  return 0;
}

#endif
