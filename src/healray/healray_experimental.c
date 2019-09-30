/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/healray/healray_experimental.c
 * \date        01/2013
 * \author      Thomas Greif
 * \brief       Adaptive ray-tracing
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 */

#include "../allvars.h"
#include "../proto.h"

int compare_idx(const void *a, const void *b);


void healray_merge()
{
  int i, j, first_i;

  int tot_num_merge, num_merge = 0;
  int all_rays;
  mesh_search_data *searchdata;

  int lvl_res = 0;

  for(i = 0; i < HRD.NumRays; i++)
    {
      if(HRD.HRRD[i].lvl > lvl_res)
        HRD.HRRD[i].hpidx_orig = HRD.HRRD[i].hpidx / (1 << 2 * (HRD.HRRD[i].lvl - lvl_res));
      else if(HRD.HRRD[i].lvl == lvl_res)
        HRD.HRRD[i].hpidx_orig = HRD.HRRD[i].hpidx;
      else
        HRD.HRRD[i].hpidx_orig = -1;
    }

  qsort((char *) HRD.HRRD, HRD.NumRays, sizeof(struct HRRD_struct), compare_idx);

  first_i = 0;

  for(i = 1; i < HRD.NumRays; i++)
    {
      if(HRD.HRRD[i].idx == HRD.HRRD[first_i].idx)
        {
          if(HRD.HRRD[first_i].hpidx_orig >= 0)
            {
              if(HRD.HRRD[i].hpidx_orig == HRD.HRRD[first_i].hpidx_orig)
                {
                  HRD.HRRD[i].iter = -1;

                  HRD.HRRD[first_i].lvl = lvl_res;
                  HRD.HRRD[first_i].hpidx = HRD.HRRD[i].hpidx_orig;

                  for(j = 0; j < 3; j++)
                    HRD.HRRD[first_i].pos[j] = P[HRD.HRRD[i].idx].Pos[j];

                  pix2vec_nest(1 << HRD.HRRD[first_i].lvl, HRD.HRRD[first_i].hpidx, HRD.HRRD[first_i].dir);

                  HRD.HRRD[first_i].energy += HRD.HRRD[i].energy;

                  num_merge++;
                }
              else
                first_i = i;
            }
          else
            first_i = i;
        }
      else
        first_i = i;
    }

  for(i = 0; i < HRD.NumRays; i++)
    {
      if(HRD.HRRD[i].iter < 0)
        {
          if(i != HRD.NumRays - 1)
            memcpy(HRD.HRRD + i, HRD.HRRD + HRD.NumRays - 1, sizeof(struct HRRD_struct));

          i--;

          HRD.NumRays--;
        }
    }

  //for(i = 0; i < HRD.NumRays; i++)
  //printf("ray 2: %d %d %d %d %d %d\n", ThisTask, HRD.HRRD[i].idx, HRD.HRRD[i].lvl, HRD.HRRD[i].hpidx, HRD.HRRD[i].hpidx_orig, HRD.HRRD[i].iter);

  //MPI_Barrier(MPI_COMM_WORLD);

  /*
     dx = HRD.HRRD[i].dir[0] * HRD.HRRD[j].dir[0];
     dy = HRD.HRRD[i].dir[1] * HRD.HRRD[j].dir[1];
     dz = HRD.HRRD[i].dir[2] * HRD.HRRD[j].dir[2];

     dphi = acos(dx + dy + dz);

     mpi_printf("match 1: %d %d %d %d %d %d %g %g %g\n", ThisTask, HRD.HRRD[i].task, HRD.HRRD[i].id, HRD.HRRD[i].hpidx, HRD.HRRD[i].iter, HRD.HRRD[i].idx, HRD.HRRD[i].dir[0], HRD.HRRD[i].dir[1], HRD.HRRD[i].dir[2]);

     mpi_printf("match 2: %d %d %d %d %d %d %g %g %g\n", ThisTask, HRD.HRRD[j].task, HRD.HRRD[j].id, HRD.HRRD[j].hpidx, HRD.HRRD[j].iter, HRD.HRRD[j].idx, HRD.HRRD[j].dir[0], HRD.HRRD[j].dir[1], HRD.HRRD[j].dir[2]);

     mpi_printf("dphi = %g\n", dphi);
   */
  /*
     double dir_old[3];

     for(i = 0; i < HRD.NumRays; i++)
     {
     if(HRD.HRRD[i].lvl < 0)
     {
     HRD.HRRD[i].lvl *= -1;

     for(j = 0; j < 3; j++)
     dir_old[j] = HRD.HRRD[i].dir[j];

     pix2vec_nest(1 << HRD.HRRD[i].lvl, HRD.HRRD[i].hpidx, HRD.HRRD[i].dir);

     for(j = 0; j < 3; j++)
     HRD.HRRD[i].pos[j] = HRD.HRRD[i].pos[j] + HRD.HRRD[i].len * (HRD.HRRD[i].dir[j] - dir_old[j]);
     }
     else if(HRD.HRRD[i].iter < 0)
     {
     if(i != HRD.NumRays - 1)
     memcpy(HRD.HRRD + i, HRD.HRRD + HRD.NumRays - 1, sizeof(struct HRRD_struct));

     i--;

     HRD.NumRays--;
     }

     //HRD.HRRD[i].energy *= nside2npix(HRD.HRRD[i].nside);
     }

     searchdata = mymalloc_movable(&searchdata, "searchdata", HRD.NumRays * sizeof(mesh_search_data));

     for(i = 0; i < HRD.NumRays; i++)
     {
     searchdata[i].Pos[0] = HRD.HRRD[i].pos[0];
     searchdata[i].Pos[1] = HRD.HRRD[i].pos[1];
     searchdata[i].Pos[2] = HRD.HRRD[i].pos[2];
     searchdata[i].u.hsmlguess = get_cell_radius(HRD.HRRD[i].idx);
     }

     find_nearest_meshpoint_global(searchdata, HRD.NumRays, 1, 0);

     mpi_distribute_items_from_search(searchdata, HRD.HRRD, &HRD.NumRays, &HRD.MaxNumRays, sizeof(struct HRRD_struct), TAG_DENS_A, offsetof(struct HRRD_struct, task), offsetof(struct HRRD_struct, idx));

     myfree_movable(searchdata);
   */
  MPI_Allreduce(&HRD.NumRays, &all_rays, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&num_merge, &tot_num_merge, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  mpi_printf("Iteration %d done, tot num merge: %d, rays after merge: %d \n", HRD.Iter, tot_num_merge, all_rays);
}


void healray_mic_test()
{
  char *env;
  int i, nthreads;
  int tot_nrays, nrays_start[NUM_THREADS], nrays[NUM_THREADS];
  double t0, t1, dt;

  t0 = second();

#ifdef HEALRAY_MIC
  env = getenv("HEALRAY_MIC_NUM_THREADS");
#else
  env = getenv("HEALRAY_NUM_THREADS");
#endif

  nthreads = atoi(env);

  omp_set_num_threads(nthreads);
  omp_set_dynamic(0);

  mpi_printf("val = %d\n", nthreads);

  for(i = tot_nrays = 0; i < nthreads; i++)
    {
      nrays_start[i] = tot_nrays;

      nrays[i] = HRD.NumRays / nthreads;

      if(i < HRD.NumRays % nthreads)
        nrays[i] += 1;

      tot_nrays += nrays[i];

      //printf("rays: %d %d %d %d %d %d\n", ThisTask, t, HRD.NumRays, nrays[t], nrays_start[t], nrays_start[t] + nrays[t] - 1);
    }

  //if(ThisTask == 0)
  //printf("rays: %d %d %d %d %d\n", ThisTask, t, nrays[t], nrays_start[t], nrays_start[t] + nrays[t] - 1);

#ifdef HEALRAY_MIC
#pragma offload target(mic:0)
#endif
  {
    printf("Hello! %d %d\n", ThisTask, nthreads);

#pragma omp parallel
    {
      int i, j, tid = omp_get_thread_num();
      int totnum = 2;
      double bla;

      //bla = malloc(totnum * sizeof(double));

      //if(ThisTask == 0)
      //printf("rays: %d %d %d %d %d\n", ThisTask, tid, nrays[tid], nrays_start[tid], nrays_start[tid] + nrays[tid] - 1);

      //printf("Hello! %d %d %d\n", ThisTask, nthreads, tid);

      for(i = nrays_start[tid]; i < nrays_start[tid] + nrays[tid]; i++)
        for(j = 0; j < totnum; j++)
          bla = exp(-j);

      //free(bla);
    }
  }

  t1 = second();

  dt = timediff(t0, t1);

  mpi_printf("dt = %g\n", dt);

  endrun();
}


int compare_idx(const void *a, const void *b)
{
  if(((struct HRRD_struct *) a)->idx < (((struct HRRD_struct *) b)->idx))
    return -1;

  if(((struct HRRD_struct *) a)->idx > (((struct HRRD_struct *) b)->idx))
    return +1;

  if(((struct HRRD_struct *) a)->idx == (((struct HRRD_struct *) b)->idx))
    {
      if(((struct HRRD_struct *) a)->hpidx_orig < (((struct HRRD_struct *) b)->hpidx_orig))
        return -1;

      if(((struct HRRD_struct *) a)->hpidx_orig > (((struct HRRD_struct *) b)->hpidx_orig))
        return +1;
    }

  return 0;
}
