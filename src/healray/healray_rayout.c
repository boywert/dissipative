/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/healray/healray_rayout.c
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

struct frac_struct
{
  int task;
  int idx;
  double energy_frac;
} *frac;


void healray_init_rayout()
{
  int i, j, flag_redo;
  long tot_num_rays, nprev, ray_idx;

  long num_rays = HRD.NumRays;

  MPI_Allreduce(&num_rays, &tot_num_rays, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

  if(tot_num_rays < HRD.RayNumTrace)
    HRD.RayNumTrace = tot_num_rays;

  if(HRD.RayNumTrace)
    {
      int *num_rays_list = mymalloc_movable(&num_rays_list, "num_rays_list", NTask * sizeof(int));

      MPI_Allgather(&HRD.NumRays, 1, MPI_INT, num_rays_list, 1, MPI_INT, MPI_COMM_WORLD);

      long prev_num_rays = 0;

      for(i = 0; i < ThisTask; i++)
        prev_num_rays += num_rays_list[i];

      long *id_list = mymalloc_movable(&id_list, "id_list", num_rays_list[ThisTask] * sizeof(long));

      int k = 0;
      long loc_ray_idx = -1;

      for(i = 0; i < NumGas; i++)
        if(HRD.EnergyInit[i] >= 0)
          {
#ifdef TGSET
            if(ThisTask == TGD.NHMaxTask && i == TGD.NHMaxIdx)
              loc_ray_idx = prev_num_rays + k;
#endif
            for(j = 0; j < HRD.NumRaysInit; j++)
              id_list[k++] = (long) P[i].ID * HRD.NumRaysInit + j;
          }

      if(k != num_rays_list[ThisTask])
        terminate("HEALRAY: counter does not match ray list entry!");

      MPI_Allreduce(&loc_ray_idx, &ray_idx, 1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD);

      for(i = 0; i < HRD.RayNumThreads; i++)
        {
          memset(HRD.TraceEnergyFrac[i], 0, HRD.TraceNumIter * HRD.RayNumTrace * sizeof(double));
#ifdef TGCHEM
          if(!TGCD.ChemMode && HRD.RayMultiFreqMode)
            memset(HRD.TraceEnergyFracFreq[i], 0, HRD.TraceNumIter * HRD.RayNumTrace * HRD.RayNumBins * sizeof(double));
#endif
        }

      long *target_list = mymalloc_movable(&target_list, "target_list", HRD.RayNumTrace * sizeof(long));
      int *task_list = mymalloc_movable(&task_list, "task_list", HRD.RayNumTrace * sizeof(int));
      int *idx_list = mymalloc_movable(&idx_list, "idx_list", HRD.RayNumTrace * sizeof(int));

      if(ThisTask == 0)
        for(i = 0; i < HRD.RayNumTrace; i++)
          {
            do
              {
                flag_redo = 0;
#ifdef TGSET
                target_list[i] = ray_idx + get_random_number() * HRD.NumRaysInit;
#endif
                //target_list[i] = get_random_number() * tot_num_rays;

                for(j = 0; j < i; j++)
                  if(target_list[i] == target_list[j])
                    flag_redo = 1;

                if(!flag_redo)
                  {
                    j = 0;
                    nprev = 0;

                    while(target_list[i] >= nprev + num_rays_list[j])
                      nprev += num_rays_list[j++];

                    task_list[i] = j;
                    idx_list[i] = target_list[i] - nprev;
                  }
              }
            while(flag_redo);
          }

      MPI_Bcast(target_list, HRD.RayNumTrace, MPI_LONG, 0, MPI_COMM_WORLD);
      MPI_Bcast(task_list, HRD.RayNumTrace, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(idx_list, HRD.RayNumTrace, MPI_INT, 0, MPI_COMM_WORLD);

      for(i = 0; i < HRD.RayNumTrace; i++)
        {
          if(ThisTask == task_list[i])
            HRD.TraceList[i] = id_list[idx_list[i]];

          MPI_Bcast(&HRD.TraceList[i], 1, MPI_LONG, task_list[i], MPI_COMM_WORLD);

          if(HRD.RayDebugMode)
            mpi_printf("id: %ld %d %ld %d %d %ld\n", tot_num_rays, HRD.RayNumTrace, target_list[i], task_list[i], idx_list[i], HRD.TraceList[i]);
        }

      myfree_movable(idx_list);
      myfree_movable(task_list);
      myfree_movable(target_list);
      myfree_movable(id_list);
      myfree_movable(num_rays_list);

      if(ThisTask == 0)
        {
          char buf[MAXLEN_PATH];
          int tgchem_flag;

          sprintf(buf, "data/healray_ray.dat");

          if(!(HRD.TraceFile = fopen(buf, "w")))
            terminate("Could not open file!\n");
#ifdef TGCHEM
          if(!TGCD.ChemMode)
            tgchem_flag = 1;
#else
          tgchem_flag = 0;
#endif
          fwrite(&tgchem_flag, sizeof(int), 1, HRD.TraceFile);
          fwrite(&HRD.RayMultiFreqMode, sizeof(int), 1, HRD.TraceFile);
          fwrite(&HRD.RayNumTrace, sizeof(int), 1, HRD.TraceFile);
          fwrite(&HRD.TraceNumIter, sizeof(int), 1, HRD.TraceFile);
#ifdef TGCHEM
          if(!TGCD.ChemMode && HRD.RayMultiFreqMode)
            {
              fwrite(&HRD.RayNumLines, sizeof(int), 1, HRD.TraceFile);
              fwrite(&HRD.RayNumFreq, sizeof(int), 1, HRD.TraceFile);
              fwrite(&HRD.NuDRange, sizeof(double), 1, HRD.TraceFile);
              fwrite(HRD.NuDPos, sizeof(double), HRD.RayNumFreq, HRD.TraceFile);
            }
#endif
        }
    }
}


void healray_finish_rayout()
{
  if(HRD.RayNumTrace)
    {
      int i, j, k, idx;

      int nfrac = 0;
      int max_nfrac = 0;

      int bins = HRD.TraceNumIter * HRD.RayNumTrace;

      int *last_iter = mymalloc_movable(&last_iter, "last_iter", HRD.RayNumTrace * sizeof(int));
      double *energy_frac = mymalloc_movable(&energy_frac, "energy_frac", bins * sizeof(double));

      frac = mymalloc_movable(&frac, "frac", 0);

      memset(energy_frac, 0, bins * sizeof(double));

      for(i = 0; i < bins; i++)
        for(j = 0; j < HRD.RayNumThreads; j++)
          energy_frac[i] += *((double *) HRD.TraceEnergyFrac[j] + i);

      for(i = 0; i < bins; i++)
        if(energy_frac[i])
          {
            if(nfrac + 1 > max_nfrac)
              {
                max_nfrac = ALLOC_INCREASE_FACTOR * max_nfrac + 1;

                frac = myrealloc_movable(frac, max_nfrac * sizeof(struct frac_struct));
              }

            frac[nfrac].task = 0;
            frac[nfrac].idx = i;
            frac[nfrac].energy_frac = energy_frac[i];

            nfrac++;
          }

      mpi_distribute_items_to_tasks(frac, offsetof(struct frac_struct, task), &nfrac, &max_nfrac, sizeof(struct frac_struct), TAG_DENS_A);

      memset(energy_frac, 0, bins * sizeof(double));

      if(ThisTask == 0)
        for(i = 0; i < nfrac; i++)
          energy_frac[frac[i].idx] = frac[i].energy_frac;

      if(ThisTask == 0)
        {
          for(i = 0; i < HRD.RayNumTrace; i++)
            {
              for(j = 0; j < HRD.TraceNumIter; j++)
                if(energy_frac[j * HRD.RayNumTrace + i])
                  last_iter[i] = j;

              if(HRD.RayDebugMode)
                mpi_printf("trace: %d %ld %d %g\n", i, HRD.TraceList[i], last_iter[i], energy_frac[last_iter[i] * HRD.RayNumTrace + i]);
            }

          fwrite(last_iter, sizeof(int), HRD.RayNumTrace, HRD.TraceFile);

          for(i = 0; i < HRD.RayNumTrace; i++)
            for(j = 0; j < HRD.TraceNumIter; j++)
              fwrite(&energy_frac[j * HRD.RayNumTrace + i], sizeof(double), 1, HRD.TraceFile);
        }
#ifdef TGCHEM
      double *energy_frac_freq;

      if(!TGCD.ChemMode && HRD.RayMultiFreqMode)
        {
          bins *= HRD.RayNumBins;

          energy_frac_freq = mymalloc_movable(&energy_frac_freq, "energy_frac_freq", bins * sizeof(double));

          memset(energy_frac_freq, 0, bins * sizeof(double));

          for(i = 0; i < bins; i++)
            for(j = 0; j < HRD.RayNumThreads; j++)
              energy_frac_freq[i] += *((double *) HRD.TraceEnergyFracFreq[j] + i);

          for(i = nfrac = max_nfrac = 0; i < bins; i++)
            if(energy_frac_freq[i])
              {
                if(nfrac + 1 > max_nfrac)
                  {
                    max_nfrac = ALLOC_INCREASE_FACTOR * max_nfrac + 1;

                    frac = myrealloc_movable(frac, max_nfrac * sizeof(struct frac_struct));
                  }

                frac[nfrac].task = 0;
                frac[nfrac].idx = i;
                frac[nfrac].energy_frac = energy_frac_freq[i];

                nfrac++;
              }

          mpi_distribute_items_to_tasks(frac, offsetof(struct frac_struct, task), &nfrac, &max_nfrac, sizeof(struct frac_struct), TAG_DENS_A);

          memset(energy_frac_freq, 0, bins * sizeof(double));

          if(ThisTask == 0)
            for(i = 0; i < nfrac; i++)
              energy_frac_freq[frac[i].idx] = frac[i].energy_frac;

          if(ThisTask == 0)
            for(i = 0; i < HRD.RayNumTrace; i++)
              for(j = 0; j < HRD.TraceNumIter; j++)
                for(k = 0; k < HRD.RayNumBins; k++)
                  {
                    idx = j * HRD.RayNumTrace * HRD.RayNumBins + i * HRD.RayNumBins + k;

                    fwrite(&energy_frac_freq[idx], sizeof(double), 1, HRD.TraceFile);
                  }
        }

      if(!TGCD.ChemMode && HRD.RayMultiFreqMode)
        myfree_movable(energy_frac_freq);
#endif

      myfree_movable(frac);
      myfree_movable(energy_frac);
      myfree_movable(last_iter);

      if(ThisTask == 0)
        fclose(HRD.TraceFile);
    }
}
