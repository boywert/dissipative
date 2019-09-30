/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/system.c
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
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <unistd.h>
#include <signal.h>
#include <gsl/gsl_rng.h>


#include "allvars.h"
#include "proto.h"



void subdivide_evenly(int N, int pieces, int index, int *first, int *count)
{
  int avg = (N - 1) / pieces + 1;
  int exc = pieces * avg - N;
  int indexlastsection = pieces - exc;

  if(index < indexlastsection)
    {
      *first = index * avg;
      *count = avg;
    }
  else
    {
      *first = index * avg - (index - indexlastsection);
      *count = avg - 1;
    }
}



#ifdef CHUNKING
void permutate_chunks_in_list(int ncount, int *list)
{
#define WALK_N_PIECES  32       /*!< Number of sets, the chunks are divided into */
#define WALK_N_SIZE    500      /*!< Number of particles per chunk */

  int nchunk;                   /*!< Number of chunk sets used */
  int nchunksize;               /*!< Size of each chunk */
  int currentchunk;             /*!< Chunk set currently processed */
  int nextparticle;

  if(ncount > WALK_N_PIECES * WALK_N_SIZE)
    {
      nchunk = WALK_N_PIECES;
      nchunksize = WALK_N_SIZE;
    }
  else
    {
      nchunk = 1;
      nchunksize = ncount;
    }

  currentchunk = 0;

  int *chunked_TargetList = (int *) mymalloc("chunked_TargetList", ncount * sizeof(int));
  int n, i;
  for(n = 0, nextparticle = 0; n < ncount; n++)
    {
      i = nextparticle;

      chunked_TargetList[n] = list[i];
      if(i < ncount)
        {
          nextparticle++;

          if((nextparticle % nchunksize) == 0)
            nextparticle += (nchunk - 1) * nchunksize;

          if(nextparticle >= ncount)
            {
              currentchunk++;
              if(currentchunk < nchunk)
                nextparticle = currentchunk * nchunksize;
            }
        }
    }

  for(n = 0; n < ncount; n++)
    list[n] = chunked_TargetList[n];

  myfree(chunked_TargetList);
}
#endif



int get_thread_num(void)
{
#if (NUM_THREADS > 1)           /* This enables OpenMP */
  return omp_get_thread_num();
#else
  return 0;
#endif
}

static struct node_data
{
  int task, this_node, first_task_in_this_node;
  int first_index, rank_in_node, tasks_in_node;
  char name[MPI_MAX_PROCESSOR_NAME];
}
loc_node, *list_of_nodes;

int system_compare_hostname(const void *a, const void *b)
{
  int cmp = strcmp(((struct node_data *) a)->name, ((struct node_data *) b)->name);

  if(cmp == 0)
    {
      if(((struct node_data *) a)->task < ((struct node_data *) b)->task)
        cmp = -1;
      else
        cmp = +1;
    }

  return cmp;
}

int system_compare_first_task(const void *a, const void *b)
{
  if(((struct node_data *) a)->first_task_in_this_node < ((struct node_data *) b)->first_task_in_this_node)
    return -1;

  if(((struct node_data *) a)->first_task_in_this_node > ((struct node_data *) b)->first_task_in_this_node)
    return +1;

  if(((struct node_data *) a)->task < ((struct node_data *) b)->task)
    return -1;

  if(((struct node_data *) a)->task > ((struct node_data *) b)->task)
    return +1;

  return 0;
}


int system_compare_task(const void *a, const void *b)
{
  if(((struct node_data *) a)->task < ((struct node_data *) b)->task)
    return -1;

  if(((struct node_data *) a)->task > ((struct node_data *) b)->task)
    return +1;

  return 0;
}



void determine_compute_nodes(void)
{
  int len, nodes, i, no, rank, first_index;

  MPI_Get_processor_name(loc_node.name, &len);
  loc_node.task = ThisTask;

  list_of_nodes = malloc(sizeof(struct node_data) * NTask);     /* Note: Internal memory allocation routines are not yet available when this function is called */

  MPI_Allgather(&loc_node, sizeof(struct node_data), MPI_BYTE, list_of_nodes, sizeof(struct node_data), MPI_BYTE, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      FILE *fd;
      if(!(fd = fopen("uses-machines.txt", "w")))
        terminate("can't write file with used machines");
      for(i = 0; i < NTask; i++)
        fprintf(fd, "%5d  %s\n", list_of_nodes[i].task, list_of_nodes[i].name);
      fclose(fd);
    }

  qsort(list_of_nodes, NTask, sizeof(struct node_data), system_compare_hostname);

  list_of_nodes[0].first_task_in_this_node = list_of_nodes[0].task;

  for(i = 1, nodes = 1; i < NTask; i++)
    {
      if(strcmp(list_of_nodes[i].name, list_of_nodes[i - 1].name) != 0)
        {
          list_of_nodes[i].first_task_in_this_node = list_of_nodes[i].task;
          nodes++;
        }
      else
        list_of_nodes[i].first_task_in_this_node = list_of_nodes[i - 1].first_task_in_this_node;
    }

  qsort(list_of_nodes, NTask, sizeof(struct node_data), system_compare_first_task);

  for(i = 0; i < NTask; i++)
    list_of_nodes[i].tasks_in_node = 0;

  for(i = 0, no = 0, rank = 0, first_index = 0; i < NTask; i++)
    {
      if(i ? list_of_nodes[i].first_task_in_this_node != list_of_nodes[i - 1].first_task_in_this_node : 0)
        {
          no++;
          rank = 0;
          first_index = i;
        }

      list_of_nodes[i].first_index = first_index;
      list_of_nodes[i].this_node = no;
      list_of_nodes[i].rank_in_node = rank++;
      list_of_nodes[first_index].tasks_in_node++;
    }

  int max_count = 0;
  int min_count = (1 << 30);

  for(i = 0; i < NTask; i++)
    {
      list_of_nodes[i].tasks_in_node = list_of_nodes[list_of_nodes[i].first_index].tasks_in_node;

      if(list_of_nodes[i].tasks_in_node > max_count)
        max_count = list_of_nodes[i].tasks_in_node;
      if(list_of_nodes[i].tasks_in_node < min_count)
        min_count = list_of_nodes[i].tasks_in_node;
    }

  qsort(list_of_nodes, NTask, sizeof(struct node_data), system_compare_task);

  TasksInThisNode = list_of_nodes[ThisTask].tasks_in_node;
  RankInThisNode = list_of_nodes[ThisTask].rank_in_node;

  ThisNode = list_of_nodes[ThisTask].this_node;

  NumNodes = nodes;
  MinTasksPerNode = min_count;
  MaxTasksPerNode = max_count;

  free(list_of_nodes);
}





void allreduce_sparse_double_sum(double *loc, double *glob, int N)
{
  int i, j, n, loc_first_n, nimport, nexport, task, ngrp;

  int *send_count = mymalloc("send_count", sizeof(int) * NTask);
  int *recv_count = mymalloc("recv_count", sizeof(int) * NTask);
  int *send_offset = mymalloc("send_offset", sizeof(int) * NTask);
  int *recv_offset = mymalloc("recv_offset", sizeof(int) * NTask);
  int *blocksize = mymalloc("blocksize", sizeof(int) * NTask);

  int blk = N / NTask;
  int rmd = N - blk * NTask;    /* remainder */
  int pivot_n = rmd * (blk + 1);


  for(task = 0, loc_first_n = 0; task < NTask; task++)
    {
      if(task < rmd)
        blocksize[task] = blk + 1;
      else
        blocksize[task] = blk;

      if(task < ThisTask)
        loc_first_n += blocksize[task];
    }

  double *loc_data = mymalloc("loc_data", blocksize[ThisTask] * sizeof(double));
  memset(loc_data, 0, blocksize[ThisTask] * sizeof(double));

  for(j = 0; j < NTask; j++)
    send_count[j] = 0;

  /* find for each non-zero element the processor where it should go for being summed */
  for(n = 0; n < N; n++)
    {
      if(loc[n] != 0)
        {
          if(n < pivot_n)
            task = n / (blk + 1);
          else
            task = rmd + (n - pivot_n) / blk;   /* note: if blk=0, then this case can not occur */

          send_count[task]++;
        }
    }

  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nexport = 0, nimport = 0, recv_offset[0] = 0, send_offset[0] = 0; j < NTask; j++)
    {
      nexport += send_count[j];
      nimport += recv_count[j];
      if(j > 0)
        {
          send_offset[j] = send_offset[j - 1] + send_count[j - 1];
          recv_offset[j] = recv_offset[j - 1] + recv_count[j - 1];
        }
    }

  struct ind_data
  {
    int n;
    double val;
  }
   *export_data, *import_data;

  export_data = mymalloc("export_data", nexport * sizeof(struct ind_data));
  import_data = mymalloc("import_data", nimport * sizeof(struct ind_data));

  for(j = 0; j < NTask; j++)
    send_count[j] = 0;

  for(n = 0; n < N; n++)
    {
      if(loc[n] != 0)
        {
          if(n < pivot_n)
            task = n / (blk + 1);
          else
            task = rmd + (n - pivot_n) / blk;   /* note: if blk=0, then this case can not occur */

          int index = send_offset[task] + send_count[task]++;
          export_data[index].n = n;
          export_data[index].val = loc[n];
        }
    }

  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)    /* note: here we also have a transfer from each task to itself (for ngrp=0) */
    {
      int recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        if(send_count[recvTask] > 0 || recv_count[recvTask] > 0)
          MPI_Sendrecv(&export_data[send_offset[recvTask]], send_count[recvTask] * sizeof(struct ind_data), MPI_BYTE,
                       recvTask, TAG_DENS_B, &import_data[recv_offset[recvTask]], recv_count[recvTask] * sizeof(struct ind_data), MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

  for(i = 0; i < nimport; i++)
    {
      int j = import_data[i].n - loc_first_n;

      if(j < 0 || j >= blocksize[ThisTask])
        terminate("j=%d < 0 || j>= blocksize[ThisTask]=%d", j, blocksize[ThisTask]);

      loc_data[j] += import_data[i].val;
    }

  myfree(import_data);
  myfree(export_data);

  /* now share the cost data across all processors */
  int *bytecounts = (int *) mymalloc("bytecounts", sizeof(int) * NTask);
  int *byteoffset = (int *) mymalloc("byteoffset", sizeof(int) * NTask);

  for(task = 0; task < NTask; task++)
    bytecounts[task] = blocksize[task] * sizeof(double);

  for(task = 1, byteoffset[0] = 0; task < NTask; task++)
    byteoffset[task] = byteoffset[task - 1] + bytecounts[task - 1];

  MPI_Allgatherv(loc_data, bytecounts[ThisTask], MPI_BYTE, glob, bytecounts, byteoffset, MPI_BYTE, MPI_COMM_WORLD);

  myfree(byteoffset);
  myfree(bytecounts);

  myfree(loc_data);
  myfree(blocksize);
  myfree(recv_offset);
  myfree(send_offset);
  myfree(recv_count);
  myfree(send_count);
}



void allreduce_sparse_imin(int *loc, int *glob, int N)
{
  int i, j, n, loc_first_n, nimport, nexport, task, ngrp;

  int *send_count = mymalloc("send_count", sizeof(int) * NTask);
  int *recv_count = mymalloc("recv_count", sizeof(int) * NTask);
  int *send_offset = mymalloc("send_offset", sizeof(int) * NTask);
  int *recv_offset = mymalloc("recv_offset", sizeof(int) * NTask);
  int *blocksize = mymalloc("blocksize", sizeof(int) * NTask);

  int blk = N / NTask;
  int rmd = N - blk * NTask;    /* remainder */
  int pivot_n = rmd * (blk + 1);


  for(task = 0, loc_first_n = 0; task < NTask; task++)
    {
      if(task < rmd)
        blocksize[task] = blk + 1;
      else
        blocksize[task] = blk;

      if(task < ThisTask)
        loc_first_n += blocksize[task];
    }

  int *loc_data = mymalloc("loc_data", blocksize[ThisTask] * sizeof(int));
  for(i = 0; i < blocksize[ThisTask]; i++)
    {
      loc_data[i] = INT_MAX;
    }

  for(j = 0; j < NTask; j++)
    send_count[j] = 0;

  /* find for each non-zero element the processor where it should go for being summed */
  for(n = 0; n < N; n++)
    {
      if(loc[n] != 0)
        {
          if(n < pivot_n)
            task = n / (blk + 1);
          else
            task = rmd + (n - pivot_n) / blk;   /* note: if blk=0, then this case can not occur */

          send_count[task]++;
        }
    }

  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nexport = 0, nimport = 0, recv_offset[0] = 0, send_offset[0] = 0; j < NTask; j++)
    {
      nexport += send_count[j];
      nimport += recv_count[j];
      if(j > 0)
        {
          send_offset[j] = send_offset[j - 1] + send_count[j - 1];
          recv_offset[j] = recv_offset[j - 1] + recv_count[j - 1];
        }
    }

  struct ind_data
  {
    int n;
    int val;
  }
   *export_data, *import_data;

  export_data = mymalloc("export_data", nexport * sizeof(struct ind_data));
  import_data = mymalloc("import_data", nimport * sizeof(struct ind_data));

  for(j = 0; j < NTask; j++)
    send_count[j] = 0;

  for(n = 0; n < N; n++)
    {
      if(loc[n] != 0)
        {
          if(n < pivot_n)
            task = n / (blk + 1);
          else
            task = rmd + (n - pivot_n) / blk;   /* note: if blk=0, then this case can not occur */

          int index = send_offset[task] + send_count[task]++;
          export_data[index].n = n;
          export_data[index].val = loc[n];
        }
    }

  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)    /* note: here we also have a transfer from each task to itself (for ngrp=0) */
    {
      int recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        if(send_count[recvTask] > 0 || recv_count[recvTask] > 0)
          MPI_Sendrecv(&export_data[send_offset[recvTask]], send_count[recvTask] * sizeof(struct ind_data), MPI_BYTE,
                       recvTask, TAG_DENS_B, &import_data[recv_offset[recvTask]], recv_count[recvTask] * sizeof(struct ind_data), MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

  for(i = 0; i < nimport; i++)
    {
      int j = import_data[i].n - loc_first_n;

      if(j < 0 || j >= blocksize[ThisTask])
        terminate("j=%d < 0 || j>= blocksize[ThisTask]=%d", j, blocksize[ThisTask]);

      loc_data[j] = imin(loc_data[j], import_data[i].val);
    }

  myfree(import_data);
  myfree(export_data);

  /* now share the cost data across all processors */
  int *bytecounts = (int *) mymalloc("bytecounts", sizeof(int) * NTask);
  int *byteoffset = (int *) mymalloc("byteoffset", sizeof(int) * NTask);

  for(task = 0; task < NTask; task++)
    bytecounts[task] = blocksize[task] * sizeof(int);

  for(task = 1, byteoffset[0] = 0; task < NTask; task++)
    byteoffset[task] = byteoffset[task - 1] + bytecounts[task - 1];

  MPI_Allgatherv(loc_data, bytecounts[ThisTask], MPI_BYTE, glob, bytecounts, byteoffset, MPI_BYTE, MPI_COMM_WORLD);

  myfree(byteoffset);
  myfree(bytecounts);

  myfree(loc_data);
  myfree(blocksize);
  myfree(recv_offset);
  myfree(send_offset);
  myfree(recv_count);
  myfree(send_count);
}






/* this function is a wrapper function for selecting the sort-function to use, and it returns the elapsed CPU time */
double mysort(void *base, size_t nel, size_t width, int (*compar) (const void *, const void *))
{
  double t0, t1;

  t0 = second();

  qsort(base, nel, width, compar);

  t1 = second();

  return timediff(t0, t1);
}

double dabs(double a)
{
  if(a < 0)
    return -a;
  else
    return a;
}

double dmax(double a, double b)
{
  if(a > b)
    return a;
  else
    return b;
}

size_t smax(size_t a, size_t b)
{
  if(a > b)
    return a;
  else
    return b;
}

double dmin(double a, double b)
{
  if(a < b)
    return a;
  else
    return b;
}

double max_array(double *a, int num_elements)
{
  int i;
  double max = -DBL_MAX;
  for(i = 0; i < num_elements; i++)
    {
      if(a[i] > max)
        {
          max = a[i];
        }
    }
  return (max);
}

int imax(int a, int b)
{
  if(a > b)
    return a;
  else
    return b;
}

int imin(int a, int b)
{
  if(a < b)
    return a;
  else
    return b;
}

int myflush(FILE * fstream)
{
#ifdef REDUCE_FLUSH
  /* do nothing */
  return 0;
#else
  return fflush(fstream);
#endif
}

int flush_everything(void)
{
#ifndef REDUCE_FLUSH
  return 0;
#else
  if(ThisTask == 0)
    {
      if((CPUThisRun - All.FlushLast) < All.FlushCpuTimeDiff)
        {
          return 0;
        }
      else
        {
          All.FlushLast = CPUThisRun;
        }
    }
  else
    {
      return 0;
    }
#endif

  mpi_printf("Flushing...\n");

  fflush(FdDomain);
  fflush(FdMemory);
  fflush(FdTimings);
  fflush(FdInfo);
  fflush(FdTimebin);
  fflush(FdBalance);
  fflush(FdCPU);
  fflush(FdEnergy);

#ifdef LOCAL_FEEDBACK
  fflush(FdLocalFeedback);
#endif

#ifdef OUTPUT_CPU_CSV
  fflush(FdCPUCSV);
#endif

#ifdef BLACK_HOLES
  fflush(FdBlackHoles);
  fflush(FdBlackHolesDetails);
  fflush(FdBlackHolesMergers);
#ifdef BH_SPIN_EVOLUTION
  fflush(FdBlackHolesSpin);
#endif
#ifdef BH_BIPOLAR_FEEDBACK
  fflush(FdBlackHolesBipolar);
#endif
#endif

#ifdef DARKENERGY
  fflush(FdDE);
#endif

#ifdef GFM_STELLAR_EVOLUTION
  fflush(FdMetalsGas);
  fflush(FdMetalsStars);
  fflush(FdMetalsTot);
  fflush(FdSN);
#endif

#if defined(FM_STAR_FEEDBACK) && defined(OUTPUT_STELLAR_FEEDBACK)
  fflush(FdFeedback);
  fflush(FdGasFeed);
#endif

#ifdef RT_ADVECT
  fflush(FdRad);
#endif

#ifdef USE_SFR
  fflush(FdSfr);
#endif

#ifdef BINARYLOG
  fflush(FdBinary);
#endif

#ifdef DG_CON_VARS_SUM_TO_FILE
  fflush(FdAngularMomentumDG);
  fflush(FdMassDG);
  fflush(FdEnergyDG);
#endif

#ifdef DG
  fflush(FdInfoDG);
#endif

#ifdef NUCLEAR_NETWORK
  fflush(FdNetwork);
#endif

#ifdef GENERAL_RELATIVITY
  fflush(FdGR);
#endif

#ifdef COSMIC_RAYS
  fflush(FdCREnergy);
#endif

#if defined(DUST_LIVE) && defined(DL_PRODUCTION)
  fflush(FdDust);
#endif

  return 1;
}

#ifdef DEBUG
#include <fenv.h>
void enable_core_dumps_and_fpu_exceptions(void)
{
#ifdef DEBUG_ENABLE_FPU_EXCEPTIONS
  /* enable floating point exceptions */

  extern int feenableexcept(int __excepts);
  feenableexcept(FE_DIVBYZERO | FE_INVALID);

  /* Note: FPU exceptions appear not to work properly
   * when the Intel C-Compiler for Linux is used
   */
#endif

  /* set core-dump size to infinity */
  struct rlimit rlim;
  getrlimit(RLIMIT_CORE, &rlim);
  rlim.rlim_cur = RLIM_INFINITY;
  setrlimit(RLIMIT_CORE, &rlim);

  /* MPICH catches the signales SIGSEGV, SIGBUS, and SIGFPE....
   * The following statements reset things to the default handlers,
   * which will generate a core file.
   */
  signal(SIGSEGV, SIG_DFL);
  signal(SIGBUS, SIG_DFL);
  signal(SIGFPE, SIG_DFL);
  signal(SIGINT, SIG_DFL);
}
#endif


void my_gsl_error_handler(const char *reason, const char *file, int line, int gsl_errno)
{
  terminate("GSL has reported an error: reason='%s', error handler called from file '%s', line %d, with error code %d", reason, file, line, gsl_errno);
}



double get_random_number(void)
{
  return gsl_rng_uniform(random_generator);
}


double get_random_number_aux(void)
{
  return gsl_rng_uniform(random_generator_aux);
}


/* returns the number of cpu-ticks in seconds that
 * have elapsed. (or the wall-clock time)
 */
double second(void)
{
  return MPI_Wtime();

  /*
   * possible alternative:
   *
   * return ((double) clock()) / CLOCKS_PER_SEC;
   *
   * but note: on AIX and presumably many other 32bit systems,
   * clock() has only a resolution of 10ms=0.01sec
   */
}

double measure_time(void)       /* strategy: call this at end of functions to account for time in this function, and before another (nontrivial) function is called */
{
  double t, dt;

  t = second();
  dt = t - WallclockTime;
  WallclockTime = t;

  return dt;
}

/* returns the time difference between two measurements
 * obtained with second(). The routine takes care of the
 * possible overflow of the tick counter on 32bit systems.
 */
double timediff(double t0, double t1)
{
  double dt;

  dt = t1 - t0;

  if(dt < 0)                    /* overflow has occured (for systems with 32bit tick counter) */
    {
#ifdef WALLCLOCK
      dt = 0;
#else
      dt = t1 + pow(2, 32) / CLOCKS_PER_SEC - t0;
#endif
    }

  return dt;
}



void minimum_large_ints(int n, long long *src, long long *res)
{
  if(src == res)
    {
      /* we need a buffer */
      long long buf[n];
      memcpy(buf, src, n * sizeof(long long));
      MPI_Allreduce(buf, res, n, MPI_LONG_LONG_INT, MPI_MIN, MPI_COMM_WORLD);
    }
  else
    MPI_Allreduce(src, res, n, MPI_LONG_LONG_INT, MPI_MIN, MPI_COMM_WORLD);
}


void sumup_large_ints_comm(int n, int *src, long long *res, MPI_Comm comm)
{
  long long lsrc[n];

  for(int i = 0; i < n; i++)
    lsrc[i] = src[i];

  MPI_Allreduce(lsrc, res, n, MPI_LONG_LONG_INT, MPI_SUM, comm);
}


void sumup_large_ints(int n, int *src, long long *res)
{
  sumup_large_ints_comm(n, src, res, MPI_COMM_WORLD);
}

void sumup_longs(int n, long long *src, long long *res)
{
  if(src == res)
    {
      /* we need a buffer */
      long long buf[n];
      memcpy(buf, src, n * sizeof(long long));
      MPI_Allreduce(buf, res, n, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
    }
  else
    MPI_Allreduce(src, res, n, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
}


void sumup_floats(int n, float *x, float *res)
{
  terminate("Remove this statement if you really want to use this function.");
  int i, j, p;
  float *numlist;

  double min_FreeBytes_glob, FreeBytes_local = 1.0 * FreeBytes;
  MPI_Allreduce(&FreeBytes_local, &min_FreeBytes_glob, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  int sum_chunksize = (int) (min_FreeBytes_glob / sizeof(float) / NTask);
  int sum_pieces = n / sum_chunksize;
  int sum_restsize = n % sum_chunksize;

  if(sum_chunksize == 0)
    terminate("min_FreeBytes_glob too small - not enough memory for sumup_floats.\n");

  for(j = 0; j < n; j++)
    res[j] = 0;

  for(p = 0; p < sum_pieces; p++)
    {
      numlist = (float *) mymalloc("numlist", NTask * sum_chunksize * sizeof(float));
      MPI_Allgather(x + p * sum_chunksize, sum_chunksize, MPI_FLOAT, numlist, sum_chunksize, MPI_FLOAT, MPI_COMM_WORLD);

      for(i = 0; i < NTask; i++)
        for(j = 0; j < sum_chunksize; j++)
          res[p * sum_chunksize + j] += numlist[i * sum_chunksize + j];
      myfree(numlist);
    }

  if(sum_restsize > 0)
    {
      numlist = (float *) mymalloc("numlist", NTask * sum_restsize * sizeof(float));
      MPI_Allgather(x + sum_pieces * sum_chunksize, sum_restsize, MPI_FLOAT, numlist, sum_restsize, MPI_FLOAT, MPI_COMM_WORLD);

      for(i = 0; i < NTask; i++)
        for(j = 0; j < sum_restsize; j++)
          res[sum_pieces * sum_chunksize + j] += numlist[i * sum_restsize + j];
      myfree(numlist);
    }
}

void sumup_doubles(int n, double *x, double *res)
{
  terminate("Remove this statement if you really want to use this function.");
  int i, j, p;
  double *numlist;

  double min_FreeBytes_glob, FreeBytes_local = 1.0 * FreeBytes;
  MPI_Allreduce(&FreeBytes_local, &min_FreeBytes_glob, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  int sum_chunksize = (int) (min_FreeBytes_glob / sizeof(double) / NTask);
  int sum_pieces = n / sum_chunksize;
  int sum_restsize = n % sum_chunksize;

  if(sum_chunksize == 0)
    terminate("min_FreeBytes_glob too small - not enough memory for sumup_doubles.\n");

  for(j = 0; j < n; j++)
    res[j] = 0;

  for(p = 0; p < sum_pieces; p++)
    {
      numlist = (double *) mymalloc("numlist", NTask * sum_chunksize * sizeof(double));
      MPI_Allgather(x + p * sum_chunksize, sum_chunksize, MPI_DOUBLE, numlist, sum_chunksize, MPI_DOUBLE, MPI_COMM_WORLD);

      for(i = 0; i < NTask; i++)
        for(j = 0; j < sum_chunksize; j++)
          res[p * sum_chunksize + j] += numlist[i * sum_chunksize + j];
      myfree(numlist);
    }

  if(sum_restsize > 0)
    {
      numlist = (double *) mymalloc("numlist", NTask * sum_restsize * sizeof(double));
      MPI_Allgather(x + sum_pieces * sum_chunksize, sum_restsize, MPI_DOUBLE, numlist, sum_restsize, MPI_DOUBLE, MPI_COMM_WORLD);

      for(i = 0; i < NTask; i++)
        for(j = 0; j < sum_restsize; j++)
          res[sum_pieces * sum_chunksize + j] += numlist[i * sum_restsize + j];
      myfree(numlist);
    }
}


size_t sizemax(size_t a, size_t b)
{
  if(a < b)
    return b;
  else
    return a;
}


void report_VmRSS(void)
{
  pid_t my_pid;
  FILE *fd;
  char buf[1024];

  my_pid = getpid();

  sprintf(buf, "/proc/%d/status", my_pid);

  if((fd = fopen(buf, "r")))
    {
      while(1)
        {
          if(fgets(buf, 500, fd) != buf)
            break;

          if(strncmp(buf, "VmRSS", 5) == 0)
            {
              printf("ThisTask=%d: %s", ThisTask, buf);
            }
          if(strncmp(buf, "VmSize", 6) == 0)
            {
              printf("ThisTask=%d: %s", ThisTask, buf);
            }
        }
      fclose(fd);
    }
}




long long report_comittable_memory(long long *MemTotal, long long *Committed_AS, long long *SwapTotal, long long *SwapFree)
{
  FILE *fd;
  char buf[1024];

  if((fd = fopen("/proc/meminfo", "r")))
    {
      while(1)
        {
          if(fgets(buf, 500, fd) != buf)
            break;

          if(bcmp(buf, "MemTotal", 8) == 0)
            {
              *MemTotal = atoll(buf + 10);
            }
          if(strncmp(buf, "Committed_AS", 12) == 0)
            {
              *Committed_AS = atoll(buf + 14);
            }
          if(strncmp(buf, "SwapTotal", 9) == 0)
            {
              *SwapTotal = atoll(buf + 11);
            }
          if(strncmp(buf, "SwapFree", 8) == 0)
            {
              *SwapFree = atoll(buf + 10);
            }
        }
      fclose(fd);
    }

  return (*MemTotal - *Committed_AS);
}



void check_maxmemsize_setting(void)
{
  int errflag = 0, errflag_tot;

  if(All.MaxMemSize > (MemoryOnNode / 1024.0 / TasksInThisNode) && RankInThisNode == 0)
    {
      printf("On node '%s', we have %d MPI ranks and at most %g MB available. This is not enough space for MaxMemSize = %g MB\n",
             loc_node.name, TasksInThisNode, MemoryOnNode / 1024.0, (double) All.MaxMemSize);
      errflag = 1;
      fflush(stdout);
    }

  MPI_Allreduce(&errflag, &errflag_tot, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#ifndef __OSX__
  if(errflag_tot)
    endrun();
#endif
}


void mpi_report_committable_memory(void)
{
  long long *sizelist, maxsize[6], minsize[6];
  double avgsize[6];
  int i, imem, mintask[6], maxtask[6];
  long long Mem[6];
  char label[512];

  Mem[0] = report_comittable_memory(&Mem[1], &Mem[2], &Mem[3], &Mem[4]);
  Mem[5] = Mem[1] - Mem[0];

  MemoryOnNode = Mem[1];

  for(imem = 0; imem < 6; imem++)
    {
      sizelist = (long long *) malloc(NTask * sizeof(long long));
      MPI_Allgather(&Mem[imem], sizeof(long long), MPI_BYTE, sizelist, sizeof(long long), MPI_BYTE, MPI_COMM_WORLD);

      for(i = 1, mintask[imem] = 0, maxtask[imem] = 0, maxsize[imem] = minsize[imem] = sizelist[0], avgsize[imem] = sizelist[0]; i < NTask; i++)
        {
          if(sizelist[i] > maxsize[imem])
            {
              maxsize[imem] = sizelist[i];
              maxtask[imem] = i;
            }
          if(sizelist[i] < minsize[imem])
            {
              minsize[imem] = sizelist[i];
              mintask[imem] = i;
            }
          avgsize[imem] += sizelist[i];
        }

      free(sizelist);
    }

  if(ThisTask == 0)
    {
      printf("\n-------------------------------------------------------------------------------------------------------------------------\n");
      for(imem = 0; imem < 6; imem++)
        {
          switch (imem)
            {
            case 0:
              sprintf(label, "AvailMem");
              break;
            case 1:
              sprintf(label, "Total Mem");
              break;
            case 2:
              sprintf(label, "Committed_AS");
              break;
            case 3:
              sprintf(label, "SwapTotal");
              break;
            case 4:
              sprintf(label, "SwapFree");
              break;
            case 5:
              sprintf(label, "AllocMem");
              break;
            }
          printf
            ("%s:\t Largest = %10.2f Mb (on task=%4d), Smallest = %10.2f Mb (on task=%4d), Average = %10.2f Mb\n",
             label, maxsize[imem] / (1024.0), maxtask[imem], minsize[imem] / (1024.0), mintask[imem], avgsize[imem] / (1024.0 * NTask));
        }
      printf("-------------------------------------------------------------------------------------------------------------------------\n");
    }

  char name[MPI_MAX_PROCESSOR_NAME];

  if(ThisTask == maxtask[2])
    {
      int len;
      MPI_Get_processor_name(name, &len);
    }

  MPI_Bcast(name, MPI_MAX_PROCESSOR_NAME, MPI_BYTE, maxtask[2], MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      printf("Task=%d has the maximum commited memory and is host: %s\n", maxtask[2], name);
      printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    }

  fflush(stdout);
}


/* Find the first bit set in the argument  */
int my_ffsll(peanokey i)
{
  int res = 0;

  while(i > 0xffffffff)
    {
      res += 32;
      i >>= 32;
    }

  return res + ffs(i);
}


/* the following function appears in the linux kernel */
int my_fls(int x)
{
  int r = 32;

  if(!x)
    return 0;
  if(!(x & 0xffff0000u))
    {
      x <<= 16;
      r -= 16;
    }
  if(!(x & 0xff000000u))
    {
      x <<= 8;
      r -= 8;
    }
  if(!(x & 0xf0000000u))
    {
      x <<= 4;
      r -= 4;
    }
  if(!(x & 0xc0000000u))
    {
      x <<= 2;
      r -= 2;
    }
  if(!(x & 0x80000000u))
    {
      x <<= 1;
      r -= 1;
    }
  return r;
}
