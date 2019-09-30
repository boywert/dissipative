/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/healray/healray_init.c
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

int healray_intpointer_compare(const void *a, const void *b);
int healray_mesh_search_compare_task(const void *a, const void *b);


void healray_comm(int mode)
{
  void *data;
  char *data_copy;
  int i, iter, size, size_hrd, offset_task, offset_idx, nimport, **perm;
  size_t orig_pos;
  mesh_search_data *searchdata, **perm_search;

  if(mode)
    {
      searchdata = mymalloc_movable(&searchdata, "searchdata", HRD.NumRays * sizeof(mesh_search_data));

      for(i = 0; i < HRD.NumRays; i++)
        {
          searchdata[i].Pos[0] = HRD.HRRD[i].pos[0];
          searchdata[i].Pos[1] = HRD.HRRD[i].pos[1];
          searchdata[i].Pos[2] = HRD.HRRD[i].pos[2];
          searchdata[i].u.hsmlguess = get_cell_radius(HRD.HRRD[i].idx);
        }

      find_nearest_meshpoint_global(searchdata, HRD.NumRays, 1, 0);
    }

  size_hrd = sizeof(struct HRRD_struct);
  offset_task = offsetof(struct HRRD_struct, task);
  offset_idx = offsetof(struct HRRD_struct, idx);

  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;

  for(i = 0; i < HRD.NumRays; i++)
    {
      if(mode)
        {
          HRD.HRRD[i].task = searchdata[i].Task;
          HRD.HRRD[i].idx = searchdata[i].u.Index;
        }

      Send_count[HRD.HRRD[i].task]++;
    }

  nimport = mpi_calculate_offsets(Send_count, Send_offset, Recv_count, Recv_offset, 0);

  if(HRD.MaxNumRays < nimport)
    {
      HRD.HRRD = myrealloc_movable(HRD.HRRD, (size_t) nimport * size_hrd);
      HRD.RayEnergyLine = myrealloc_movable(HRD.RayEnergyLine, (size_t) nimport * HRD.RayNumLines * sizeof(double));
      HRD.RayEnergyLineInit = myrealloc_movable(HRD.RayEnergyLineInit, (size_t) nimport * HRD.RayNumLines * sizeof(double));
      HRD.RayEnergy = myrealloc_movable(HRD.RayEnergy, (size_t) nimport * HRD.RayNumBins * sizeof(float));

      HRD.MaxNumRays = nimport;
    }

  if(mode)
    perm_search = mymalloc("perm_search", HRD.NumRays * sizeof(*perm_search));
  else
    perm = mymalloc("perm", HRD.NumRays * sizeof(*perm));

  for(i = 0; i < HRD.NumRays; ++i)
    {
      if(mode)
        perm_search[i] = &searchdata[i];
      else
        perm[i] = (int *) ((char *) HRD.HRRD + (size_t) i * size_hrd + offset_task);
    }

  if(mode)
    mysort(perm_search, HRD.NumRays, sizeof(*perm_search), healray_mesh_search_compare_task);
  else
    mysort(perm, HRD.NumRays, sizeof(*perm), healray_intpointer_compare);

  for(iter = 0; iter < 4; iter++)
    {
      if(iter == 0)
        {
          data = HRD.RayEnergyLine;
          size = HRD.RayNumLines * sizeof(double);
        }

      if(iter == 1)
        {
          data = HRD.RayEnergyLineInit;
          size = HRD.RayNumLines * sizeof(double);
        }

      if(iter == 2)
        {
          data = HRD.RayEnergy;
          size = HRD.RayNumBins * sizeof(float);
        }

      if(iter == 3)
        {
          data = HRD.HRRD;
          size = size_hrd;
        }

      data_copy = mymalloc_movable(&data_copy, "data_copy", (size_t) HRD.NumRays * size);

      for(i = 0; i < HRD.NumRays; ++i)
        {
          if(mode)
            orig_pos = perm_search[i] - searchdata;
          else
            orig_pos = ((char *) perm[i] - ((char *) HRD.HRRD + offset_task)) / size_hrd;

          memcpy(data_copy + size * (size_t) i, (char *) data + size * orig_pos, size);
        }

      mpi_exchange_buffers(data_copy, Send_count, Send_offset, data, Recv_count, Recv_offset, size, TAG_DENS_A, 1);

      myfree_movable(data_copy);
    }

  HRD.NumRays = nimport;

  if(mode)
    {
      myfree(perm_search);

      myfree_movable(searchdata);
    }
  else
    myfree(perm);
}


void healray_split_comm()
{
  void *data;
  char *data_copy;
  int i, iter, task, size, size_hrd, nimport;
  size_t orig_pos;
  mesh_search_data *searchdata, **perm;

  searchdata = mymalloc_movable(&searchdata, "searchdata", HRD.TmpNumRays * sizeof(mesh_search_data));

  for(i = 0; i < HRD.TmpNumRays; i++)
    {
      searchdata[i].Pos[0] = HRD.TmpHRRD[i].pos[0];
      searchdata[i].Pos[1] = HRD.TmpHRRD[i].pos[1];
      searchdata[i].Pos[2] = HRD.TmpHRRD[i].pos[2];
      searchdata[i].u.hsmlguess = get_cell_radius(HRD.TmpHRRD[i].idx);
    }

  find_nearest_meshpoint_global(searchdata, HRD.TmpNumRays, 1, 0);

  //mpi_distribute_items_from_search(searchdata, HRD.TmpHRRD, &HRD.TmpNumRays, &HRD.TmpMaxNumRays, sizeof(struct HRRD_struct), TAG_DENS_A, offsetof(struct HRRD_struct, task), offsetof(struct HRRD_struct, idx));

  size_hrd = sizeof(struct HRRD_struct);

  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;

  for(i = 0; i < HRD.TmpNumRays; i++)
    {
      task = searchdata[i].Task;

      Send_count[task]++;

      HRD.TmpHRRD[i].task = task;
      HRD.TmpHRRD[i].idx = searchdata[i].u.Index;
    }

  nimport = mpi_calculate_offsets(Send_count, Send_offset, Recv_count, Recv_offset, 0);

  if(HRD.TmpMaxNumRays < nimport)
    {
      HRD.TmpHRRD = myrealloc_movable(HRD.TmpHRRD, (size_t) nimport * size_hrd);
      HRD.TmpRayEnergyLine = myrealloc_movable(HRD.TmpRayEnergyLine, (size_t) nimport * HRD.RayNumLines * sizeof(double));
      HRD.TmpRayEnergyLineInit = myrealloc_movable(HRD.TmpRayEnergyLineInit, (size_t) nimport * HRD.RayNumLines * sizeof(double));
      HRD.TmpRayEnergy = myrealloc_movable(HRD.TmpRayEnergy, (size_t) nimport * HRD.RayNumBins * sizeof(float));

      HRD.TmpMaxNumRays = nimport;
    }

  perm = mymalloc("perm", HRD.TmpNumRays * sizeof(*perm));

  for(i = 0; i < HRD.TmpNumRays; ++i)
    perm[i] = &searchdata[i];

  mysort(perm, HRD.TmpNumRays, sizeof(*perm), healray_mesh_search_compare_task);

  for(iter = 0; iter < 4; iter++)
    {
      if(iter == 0)
        {
          data = HRD.TmpRayEnergyLine;
          size = HRD.RayNumLines * sizeof(double);
        }

      if(iter == 1)
        {
          data = HRD.TmpRayEnergyLineInit;
          size = HRD.RayNumLines * sizeof(double);
        }

      if(iter == 2)
        {
          data = HRD.TmpRayEnergy;
          size = HRD.RayNumBins * sizeof(float);
        }

      if(iter == 3)
        {
          data = HRD.TmpHRRD;
          size = size_hrd;
        }

      data_copy = mymalloc("data_copy", (size_t) HRD.TmpNumRays * size);

      for(i = 0; i < HRD.TmpNumRays; ++i)
        {
          orig_pos = perm[i] - searchdata;

          memcpy(data_copy + size * (size_t) i, (char *) data + size * orig_pos, size);
        }

      mpi_exchange_buffers(data_copy, Send_count, Send_offset, data, Recv_count, Recv_offset, size, TAG_DENS_A, 1);

      myfree(data_copy);
    }

  HRD.TmpNumRays = nimport;

  myfree(perm);

  myfree_movable(searchdata);
}


/*
void healray_alt_comm()
{
  int nrays_comm, nrays_max, nrays_left, tot_nrays_left, nrays_init, nrays_send, nrays_recv;
  double comm_mbytes, min_comm_mbytes;

  comm_mbytes = 0.1 * (double) FreeBytes / (1024 * 1024);

  MPI_Allreduce(&comm_mbytes, &min_comm_mbytes, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  if(min_comm_mbytes < 5.)
    mpi_printf("HEALRAY: Warning: very little space left for ray communication! (%g MB)", min_comm_mbytes);

  nrays_comm = nrays_max = imax((int) (comm_mbytes / ((double) sizeof(struct HRRD_struct) / (1024 * 1024))), 1);

  HRD.TmpHRRD = mymalloc_movable(&HRD.TmpHRRD, "HRD.TmpHRRD", nrays_comm * sizeof(struct HRRD_struct));

  nrays_left = HRD.NumRays;

  MPI_Allreduce(&nrays_left, &tot_nrays_left, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  while(tot_nrays_left)
    {
      nrays_send = imin(nrays_comm, nrays_left);

      memcpy(HRD.TmpHRRD, HRD.HRRD, nrays_send * sizeof(struct HRRD_struct));

      memmove(HRD.HRRD, HRD.HRRD + nrays_send, (HRD.NumRays - nrays_send) * sizeof(struct HRRD_struct));

      nrays_recv = nrays_send;

      mpi_distribute_items_to_tasks(HRD.TmpHRRD, offsetof(struct HRRD_struct, task), &nrays_recv, &nrays_max, sizeof(struct HRRD_struct), TAG_DENS_A);

      HRD.NumRays += nrays_recv - nrays_send;

      while(HRD.NumRays > HRD.MaxNumRays)
	{
	  HRD.MaxNumRays = ALLOC_INCREASE_FACTOR * HRD.MaxNumRays + 1;

	  HRD.HRRD = myrealloc_movable(HRD.HRRD, HRD.MaxNumRays * sizeof(struct HRRD_struct));
	}

      memcpy(HRD.HRRD + HRD.NumRays - nrays_recv, HRD.TmpHRRD, nrays_recv * sizeof(struct HRRD_struct));

      nrays_left -= nrays_send;

      MPI_Allreduce(&nrays_left, &tot_nrays_left, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    }

  myfree_movable(HRD.TmpHRRD);
}
*/


int healray_intpointer_compare(const void *a, const void *b)
{
  if((**(int **) a) < (**(int **) b))
    return -1;

  if((**(int **) a) > (**(int **) b))
    return +1;

  return 0;
}


int healray_mesh_search_compare_task(const void *a, const void *b)
{
  if((*(mesh_search_data **) a)->Task < (*(mesh_search_data **) b)->Task)
    return -1;

  if((*(mesh_search_data **) a)->Task > (*(mesh_search_data **) b)->Task)
    return +1;

  return 0;
}
