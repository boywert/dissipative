/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/cuda_util.c
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
#include "string.h"
#include <unistd.h>

#include "allvars.h"
#include "proto.h"


#ifdef CUDA
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cudaProfiler.h>
#include "cuda_util.h"


static CUdevice cuDev;
static CUcontext cuCtx;


void __cudaSafeCall(cudaError_t err, const char *file, const int line)
{

  if(cudaSuccess != err)
    {
      terminate("cudaSafeCall() failed at %s:%i : %s\n", file, line, cudaGetErrorString(err));
    }

  return;
}

void __cuDrvSafeCall(CUresult err, const char *file, const int line)
{

  if(err != CUDA_SUCCESS)
    {
      terminate("cuDrvSafeCall() failed at %s:%i with error code %d\n", file, line, err);
    }

  return;
}

void __cudaCheckError(const char *file, const int line)
{

  cudaError_t err = cudaGetLastError();
  if(cudaSuccess != err)
    {
      terminate("cudaCheckError() failed at %s:%i : %s.\n", file, line, cudaGetErrorString(err));
    }

  // More careful checking. However, this will affect performance.
  // Comment if not needed.
  /*err = cudaThreadSynchronize();
     if (cudaSuccess != err)
     {
     terminate("cudaCheckError() with sync failed at %s:%i : %s.\n", file, line, cudaGetErrorString(err));
     }
   */
  return;
}

typedef struct host
{
  char name[MPI_MAX_PROCESSOR_NAME];
  int rank;
} host;

int cuda_cmphost(const void *p1, const void *p2)
{
  int i;
  if((i = strcmp(((struct host *) p1)->name, ((struct host *) p2)->name)) != 0)
    return i;
  if(((struct host *) p1)->rank < (((struct host *) p2)->rank))
    return -1;

  if(((struct host *) p1)->rank > (((struct host *) p2)->rank))
    return +1;

  return 0;
}

int cuda_selectDevice()
{
  char host_name[MPI_MAX_PROCESSOR_NAME];
  size_t namelen = MPI_MAX_PROCESSOR_NAME;

  //MPI_Get_processor_name(host_name,&namelen);
  gethostname(host_name, namelen);

  host hosts[NTask];
  host this_host;

  strcpy(this_host.name, host_name);
  this_host.rank = ThisTask;

  MPI_Allgather(&this_host, sizeof(host), MPI_BYTE, hosts, sizeof(host), MPI_BYTE, MPI_COMM_WORLD);

  qsort(hosts, NTask, sizeof(host), cuda_cmphost);

  int count = 0;
  int i;
  for(i = 0; i < NTask; i++)
    {

      if(i > 0 && !strcmp(hosts[i - 1].name, hosts[i].name))
        {
          count++;
        }
      else
        {
          count = 0;
        }
      if(hosts[i].rank == ThisTask)
        {
          break;
        }

    }

  printf("Assigning device %d  to process %d on node %s\n", count, ThisTask, host_name);

  return count;

}

void cuda_init()
{
  CuDrvSafeCall(cuInit(0));

  CuDrvSafeCall(cuDeviceGet(&cuDev, cuda_selectDevice()));
  CuDrvSafeCall(cuCtxCreate(&cuCtx, 0, cuDev));
}



void cuda_finish()
{
  CuDrvSafeCall(cuCtxDestroy(cuCtx));
}

#endif
