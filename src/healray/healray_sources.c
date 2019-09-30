/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/healray/healray_sources.c
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


void healray_read_sources()
{
  char fname[MAXLEN_PATH], buf[MAXLEN_PATH];
  int i, read_ok, read_ID;
  FILE *file;

  if(ThisTask == 0)
    {
      sprintf(fname, "data/%s", HRD.RaySourceFile);

      if(!(file = fopen(fname, "r")))
        {
          sprintf(buf, "Can't open HEALRAY sources file `%s'!", fname);
          terminate(buf);
        }

      while(1)
        {
          read_ok = fscanf(file, "%d\n", &read_ID);

          if(read_ok >= 0)
            HRD.TotNumSources++;
          else
            break;
        }
    }

  MPI_Bcast(&HRD.TotNumSources, 1, MPI_INT, 0, MPI_COMM_WORLD);

  HRSL = mymalloc_movable(&HRSL, "HRSL", HRD.TotNumSources * sizeof(struct HRSL_struct));

  if(ThisTask == 0)
    {
      rewind(file);

      for(i = 0; i < HRD.TotNumSources; i++)
        fscanf(file, "%d\n", &HRSL[i].id);

      fclose(file);
    }

  MPI_Bcast(HRSL, HRD.TotNumSources * sizeof(struct HRSL_struct), MPI_BYTE, 0, MPI_COMM_WORLD);

  healray_update_sources();
}


void healray_update_sources()
{
  int i, j, k;
  double pos[3], hsml;

  for(i = 0; i < HRD.TotNumSources; i++)
    {
      for(j = 0; j < 3; j++)
        pos[j] = 0;

      for(j = 0; j < NumGas; j++)
        if(P[j].ID == HRSL[i].id)
          {
            for(k = 0; k < 3; k++)
              pos[k] = P[j].Pos[k];

            hsml = get_cell_radius(j);
          }

      MPI_Allreduce(pos, HRSL[i].pos, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&hsml, &HRSL[i].hsml, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
}


void healray_init_sources()
{
  int i, j, k;
  mesh_search_data *searchdata;

  healray_update_sources();

  HRD.NumSources = 0;

  HRD.HRSD = mymalloc_movable(&HRD.HRSD, "HRD.HRSD", HRD.TotNumSources * sizeof(struct HRSD_struct));

  for(i = 0; i < HRD.TotNumSources; i++)
    {
      for(j = 0; j < NumGas; j++)
        if(P[j].ID == HRSL[i].id)
          {
            HRD.HRSD[HRD.NumSources].idx = j;
            HRD.HRSD[HRD.NumSources].task = ThisTask;

            for(k = 0; k < 3; k++)
              HRD.HRSD[HRD.NumSources].pos[k] = P[j].Pos[k];

            HRD.NumSources++;
          }
    }

  searchdata = mymalloc_movable(&searchdata, "searchdata", HRD.NumSources * sizeof(mesh_search_data));

  for(i = 0; i < HRD.NumSources; i++)
    {
      searchdata[i].Pos[0] = HRD.HRSD[i].pos[0];
      searchdata[i].Pos[1] = HRD.HRSD[i].pos[1];
      searchdata[i].Pos[2] = HRD.HRSD[i].pos[2];
      searchdata[i].u.hsmlguess = 0.1;
    }

  find_nearest_meshpoint_global(searchdata, HRD.NumSources, 1, 0);

  mpi_distribute_items_from_search(searchdata, HRD.HRSD, &HRD.NumSources, &HRD.TotNumSources, sizeof(struct HRSD_struct), TAG_DENS_A,
                                   offsetof(struct HRSD_struct, task), offsetof(struct HRSD_struct, idx));

  myfree_movable(searchdata);

  for(i = 0; i < HRD.NumSources; i++)
    {
      if(HRD.NumRays + HRD.NumRaysInit > HRD.MaxNumRays)
        {
          HRD.MaxNumRays = ALLOC_INCREASE_FACTOR * HRD.MaxNumRays + HRD.NumRaysInit;

          HRD.HRRD = myrealloc_movable(HRD.HRRD, HRD.MaxNumRays * sizeof(struct HRRD_struct));
        }

      for(j = 0; j < HRD.NumRaysInit; j++)
        {
          HRD.HRRD[HRD.NumRays].iter = 0;
          HRD.HRRD[HRD.NumRays].task = ThisTask;
          HRD.HRRD[HRD.NumRays].task_orig = HRD.HRRD[HRD.NumRays].task;
          HRD.HRRD[HRD.NumRays].idx = HRD.HRSD[i].idx;
          HRD.HRRD[HRD.NumRays].idx_prev = -1;
          HRD.HRRD[HRD.NumRays].idx_orig = HRD.HRRD[HRD.NumRays].idx;
          HRD.HRRD[HRD.NumRays].lvl = HRD.RayLvlInit;
          HRD.HRRD[HRD.NumRays].hpidx = j;
          HRD.HRRD[HRD.NumRays].id = 0;
          HRD.HRRD[HRD.NumRays].len = 0;

          for(k = 0; k < 3; k++)
            {
              HRD.HRRD[HRD.NumRays].pos[k] = HRD.HRSD[i].pos[k];
              HRD.HRRD[HRD.NumRays].pos_old[k] = HRD.HRRD[HRD.NumRays].pos[k];
            }

          pix2vec_nest(1 << HRD.RayLvlInit, j, HRD.HRRD[HRD.NumRays].dir);

          HRD.HRRD[HRD.NumRays].energy = HRD.SourceLuminosity / HRD.NumRaysInit / HRD.EnergyUnit;
          HRD.HRRD[HRD.NumRays].energy_init = HRD.HRRD[HRD.NumRays].energy;
          HRD.HRRD[HRD.NumRays].energy_init_tot = HRD.HRRD[HRD.NumRays].energy_init * HRD.NumRaysInit;

          HRD.NumRays++;
        }
    }
}
