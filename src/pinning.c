/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/pinning.c
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
#include <unistd.h>
#include <sys/types.h>
#include <sys/syscall.h>
#include "allvars.h"
#include "proto.h"


#ifdef IMPOSE_PINNING

#include <hwloc.h>
#include <hwloc/bitmap.h>

#define MAX_CORES 4096

static int flag_pinning_error = 0;

static hwloc_cpuset_t cpuset, cpuset_after_MPI_init;
static hwloc_topology_t topology;
static int topodepth;
static int sockets;
static int cores;
static int pus;
static int hyperthreads_per_core;


void get_core_set(void)
{
  cpuset = hwloc_bitmap_alloc();
  hwloc_get_proc_cpubind(topology, getpid(), cpuset, 0);
}

void detect_topology(void)
{
  unsigned depth;

  /* Allocate and initialize topology object. */
  hwloc_topology_init(&topology);

  /* Perform the topology detection. */
  hwloc_topology_load(topology);

  /* Get some additional topology information
     in case we need the topology depth later. */
  topodepth = hwloc_topology_get_depth(topology);

  depth = hwloc_get_type_depth(topology, HWLOC_OBJ_SOCKET);

  if(depth == HWLOC_TYPE_DEPTH_UNKNOWN)
    sockets = -1;
  else
    sockets = hwloc_get_nbobjs_by_depth(topology, depth);

  depth = hwloc_get_type_depth(topology, HWLOC_OBJ_CORE);

  if(depth == HWLOC_TYPE_DEPTH_UNKNOWN)
    cores = -1;
  else
    cores = hwloc_get_nbobjs_by_depth(topology, depth);

  depth = hwloc_get_type_depth(topology, HWLOC_OBJ_PU);

  if(depth == HWLOC_TYPE_DEPTH_UNKNOWN)
    pus = -1;
  else
    pus = hwloc_get_nbobjs_by_depth(topology, depth);
}

void pin_to_core_set(void)
{
  int i, num_threads, thread;
  char buf[MAX_CORES + 1];
  char *p = getenv("OMP_NUM_THREADS");
  if(p)
    num_threads = atoi(p);
  else
    num_threads = 1;

#ifdef NUM_THREADS
  if(num_threads != NUM_THREADS)
    {
      terminate("You have activated NUM_THREADS, but the value of the environment variable OMP_NUM_THREADS=%d is not equal to NUM_THREADS=%d!", num_threads, NUM_THREADS);
    }
#endif

  mpi_printf("\n\n");
  mpi_printf("PINNING: We have %d sockets, %d physical cores and %d logical cores on the first MPI-task's node.\n", sockets, cores, pus);
  if(cores <= 0 || sockets <= 0 || pus <= 0)
    {
      mpi_printf("PINNING: The topology cannot be recognized. We refrain from any pinning attempt.\n");
      flag_pinning_error = 1;
      return;
    }

  hyperthreads_per_core = pus / cores;

  if(hyperthreads_per_core < 1)
    terminate("Need at least one logical thread per physical core\n");

  if(pus > cores)
    mpi_printf("PINNING: Looks like %d hyperthreads per physical core are in principle possible.\n", hyperthreads_per_core);

  cpuset_after_MPI_init = hwloc_bitmap_alloc();
  hwloc_get_proc_cpubind(topology, getpid(), cpuset_after_MPI_init, 0);

  if(!hwloc_bitmap_isequal(cpuset, cpuset_after_MPI_init))
    mpi_printf("PINNING: Apparently, the MPI library set some pinning itself. We'll override this.\n");


  int id, available_pus = 0;

  for(id = hwloc_bitmap_first(cpuset); id != -1; id = hwloc_bitmap_next(cpuset, id))
    available_pus++;

  mpi_printf("PINNING: Looks like %d logical cores are available\n", available_pus);

  if(available_pus == pus)
    mpi_printf("PINNING: Looks like all available logical cores are at our disposal.\n");
  else
    {
      if(available_pus >= 1)
        {
          mpi_printf("PINNING: Looks like allready before start of the code, a tight binding was imposed.\n");
#ifdef IMPOSE_PINNING_OVERRIDE_MODE
          for(id = 0; id < pus; id++)
            hwloc_bitmap_set(cpuset, id);
          available_pus = pus;
          mpi_printf("PINNING: We are overridung this and make all %d available to us.\n", available_pus);
#else
          mpi_printf("PINNING: We refrain from any pinning attempt ourselves. (This can be changed by setting USE_PINNING_OVERRIDE_MODE.)\n");
          flag_pinning_error = 1;
          return;
#endif
        }
    }

  for(i = 0; i < pus && i < MAX_CORES; i++)
    if(hwloc_bitmap_isset(cpuset, i))
      buf[i] = '1';
    else
      buf[i] = '-';
  buf[pus] = 0;

  mpi_printf("PINNING: Available logical cores on first node:  %s\n", buf);

  int pus_per_task = available_pus / TasksInThisNode;

  mpi_printf("PINNING: %d logical cores are available per MPI Task.\n", pus_per_task);

  if(pus_per_task <= 0)
    terminate("Need at least one logical core per MPI task for pinning to make sense.  available_pus=%d TasksInThisNode=%d\n", available_pus, TasksInThisNode);



  int depth, cid, cores_before, id_this, id_found, count;
  hwloc_obj_t obj;
  hwloc_cpuset_t cpuset_core;

  /* go through all logical cores in sequence of proximity */
  depth = hwloc_get_type_depth(topology, HWLOC_OBJ_PU);

  for(cid = 0, cores_before = 0; cores_before < RankInThisNode * pus_per_task && cid < pus; cid++)
    {
      obj = hwloc_get_obj_by_depth(topology, depth, cid);

      cpuset_core = hwloc_bitmap_dup(obj->cpuset);
      if(hwloc_bitmap_isincluded(cpuset_core, cpuset))
        {
          cores_before++;
        }
      hwloc_bitmap_free(cpuset_core);
    }

  int pus_per_thread, skip;

  if(pus_per_task > NUM_THREADS)
    pus_per_thread = pus_per_task / NUM_THREADS;
  else
    pus_per_thread = 1;

  /* cid should now be the logical index of the first PU for this MPI task */
  for(thread = 0, id_this = id_found = cid, count = 0; thread < NUM_THREADS; thread++)
    {
      obj = hwloc_get_obj_by_depth(topology, depth, id_found);
      cpuset_thread[thread] = hwloc_bitmap_dup(obj->cpuset);

      for(skip = 0; skip < pus_per_thread; skip++)
        {
          id_this++;
          count++;

          id_found = -1;
          if(count >= pus_per_task)
            {
              id_this = cid;
              count = 0;
            }
          do
            {
              obj = hwloc_get_obj_by_depth(topology, depth, id_this);
              cpuset_core = hwloc_bitmap_dup(obj->cpuset);
              if(hwloc_bitmap_isincluded(cpuset_core, cpuset))
                {
                  id_found = id_this;
                }
              else
                {
                  id_this++;
                  if(id_this >= pus)
                    terminate("id_this >= pus");
                }
              hwloc_bitmap_free(cpuset_core);
            }
          while(id_found < 0);
        }
    }

  hwloc_set_proc_cpubind(topology, getpid(), cpuset_thread[0], HWLOC_CPUBIND_PROCESS);
}


void report_pinning(void)
{
  int i;
  char buf[MAX_CORES + 1];

  if(flag_pinning_error)
    return;

  hwloc_get_cpubind(topology, cpuset, 0);

  for(i = 0; i < pus && i < MAX_CORES; i++)
    if(hwloc_bitmap_isset(cpuset, i))
      buf[i] = '1';
    else
      buf[i] = '-';
  buf[pus] = 0;

  for(i = 0; i < NTask; i++)
    {
      if(ThisTask == i && ThisNode == 0)
        printf("PINNING: Node=%4d: Task=%04d:                   %s\n", ThisNode, ThisTask, buf);
      fflush(stdout);
      MPI_Barrier(MPI_COMM_WORLD);
    }
}


#if (NUM_THREADS > 1)

void set_pinning_openmp_threads(void)
{
  int threadid;

  if(flag_pinning_error)
    return;

#pragma omp parallel private(threadid)
  {
    threadid = omp_get_thread_num();

    hwloc_set_cpubind(topology, cpuset_thread[threadid], HWLOC_CPUBIND_THREAD);
  }
}


void report_pinning_openmp_threads(void)
{
  int i, task, tid;
  char buf[MAX_CORES + 1];
  hwloc_cpuset_t cpuset_thr;

  if(flag_pinning_error)
    return;

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  for(task = 0; task < NTask; task++)
    {
#pragma omp parallel private(tid, cpuset_thr, i, buf)
      {
        tid = omp_get_thread_num();
        cpuset_thr = hwloc_bitmap_alloc();
        hwloc_get_cpubind(topology, cpuset_thr, HWLOC_CPUBIND_THREAD);

        for(i = 0; i < pus && i < MAX_CORES; i++)
          if(hwloc_bitmap_isset(cpuset_thr, i))
            buf[i] = '1';
          else
            buf[i] = '-';
        buf[pus] = 0;

        if(ThisTask == task && ThisNode == 0)
          printf("PINNING: Node=%4d: Task=%04d OpenMP-Thread= %02d: %s\n", ThisNode, ThisTask, tid, buf);
        fflush(stdout);

        hwloc_bitmap_free(cpuset_thr);
      }

      MPI_Barrier(MPI_COMM_WORLD);
    }
}
#endif


#endif
