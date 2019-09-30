/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/healray/healray.c
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

void healray_merge();
void healray_split();
void healray_advance();
void healray_get_alpha_and_source(int idx, double *alpha, double *source);
void healray_add_outray(int tid, int i);


void healray()
{
  int do_ray_step, tot_iter, max_num_rays;
  long long tot_num_rays, all_num_rays;
  double avg_num_rays, imba, avg_imba;
  double t0, t1, t2, t3, dt, dt_init, dt_advance, dt_imba, dt_comm, dt_rest, dt_sum;

  CPU_Step[CPU_MISC] += measure_time();

  do_ray_step = 0;

  if(HRD.SourceFlag >= 0)
    if(!HRD.RayTimeFac || All.NumCurrentTiStep % HRD.RayTimeFac == 0)
      do_ray_step = 1;

  if(do_ray_step)
    {
      t0 = second();

      mpi_printf("HEALRAY: Doing healray step...\n");

      healray_init_step();

      HRD.Loop = -1;
      tot_iter = -1;

      all_num_rays = 0;
      avg_imba = 0.;

      do
        {
          HRD.Loop++;

          t2 = second();

          if(HRD.RayMode)
            healray_init_rays_alt();
          else
            healray_init_rays();

          t3 = second();

          HRD.DtInit += timediff(t2, t3);

          HRD.Iter = -1;

          do
            {
              HRD.Iter++;
              tot_iter++;

              if(0)
                healray_merge();

              if(HRD.RaySplitFac)
                healray_split();

              sumup_large_ints(1, &HRD.NumRays, &tot_num_rays);

              MPI_Allreduce(&HRD.NumRays, &max_num_rays, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

              avg_num_rays = (double) tot_num_rays / NTask;
              all_num_rays += tot_num_rays;

              if(avg_num_rays)
                imba = max_num_rays / avg_num_rays;
              else
                imba = 1.;

              avg_imba += tot_num_rays * imba;

              healray_advance();

              sumup_large_ints(1, &HRD.NumRays, &tot_num_rays);
            }
          while(tot_num_rays);
        }
      while(HRD.TotNumEmitters);

      healray_finish_step();

      if(all_num_rays)
        avg_imba /= all_num_rays;
      else
        avg_imba = 1.;

      t1 = second();

      dt = timediff(t0, t1);

      MPI_Allreduce(&HRD.DtInit, &dt_init, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&HRD.DtAdvance, &dt_advance, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&HRD.DtImba, &dt_imba, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&HRD.DtComm, &dt_comm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&dt, &dt_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      dt_init /= NTask;
      dt_advance /= NTask;
      dt_imba /= NTask;
      dt_comm /= NTask;
      dt_sum /= NTask;

      dt_rest = dmax(dt_sum - (dt_init + dt_advance + dt_comm + dt_imba), 0.);
#ifndef TGCHEM
      if(HRD.RayNumSources)
        mpi_printf("HEALRAY: Error = %g\n", HRD.Error);
#endif
      mpi_printf("HEALRAY: Done! %d loops, %d iterations, %g rays, %g operations, (%g/%g) advances/comms per ray\n", HRD.Loop + 1, tot_iter,
                 (double) HRD.TotInitNumRays, (double) HRD.NumOps, HRD.AdvancesPerRay, HRD.CommsPerRay);
      mpi_printf("HEALRAY: Ray Imbalance: %g, Init: %g, Advance: %g, Imba: %g, Comm: %g, Rest: %g, Total: %g seconds\n\n", avg_imba, dt_init, dt_advance, dt_imba, dt_comm, dt_rest, dt_sum);
#ifndef TGCHEM
      if(!HRD.RayTimeFac)
        {
          produce_dump();

          endrun();
        }
#endif
    }

  CPU_Step[CPU_HEALRAY] += measure_time();
}


void healray_split()
{
  int i, j, k, l, idx, rmray_flag, flag, glob_flag = 1;
  long new_num_rays, tot_new_num_rays;
  double r, d, cell_area, ray_area, dir[3];

  while(glob_flag)
    {
      HRD.TmpNumRays = HRD.TmpMaxNumRays = flag = 0;

      HRD.TmpHRRD = mymalloc_movable(&HRD.TmpHRRD, "TmpHRRD", HRD.TmpMaxNumRays * sizeof(struct HRRD_struct));
      HRD.TmpRayEnergyLine = mymalloc_movable(&HRD.TmpRayEnergyLine, "TmpRayEnergyLine", HRD.TmpMaxNumRays * HRD.RayNumLines * sizeof(double));
      HRD.TmpRayEnergyLineInit = mymalloc_movable(&HRD.TmpRayEnergyLineInit, "TmpRayEnergyLineInit", HRD.TmpMaxNumRays * HRD.RayNumLines * sizeof(double));
      HRD.TmpRayEnergy = mymalloc_movable(&HRD.TmpRayEnergy, "TmpRayEnergy", HRD.TmpMaxNumRays * HRD.RayNumBins * sizeof(float));

      for(i = 0; i < HRD.NumRays; i++)
        {
          idx = HRD.HRRD[i].idx;

          r = get_cell_radius(idx);

          cell_area = M_PI * r * r;

          d = dmax(HRD.HRRD[i].len, r);

          ray_area = 4 * M_PI * d * d / (12 * (1 << 2 * HRD.HRRD[i].lvl));

          if(HRD.RaySplitFac * ray_area > cell_area)
            {
              flag = 1;

              if(HRD.TmpNumRays + 4 > HRD.TmpMaxNumRays)
                {
                  HRD.TmpMaxNumRays = ALLOC_INCREASE_FACTOR * HRD.TmpMaxNumRays + 4;

                  HRD.TmpHRRD = myrealloc_movable(HRD.TmpHRRD, HRD.TmpMaxNumRays * sizeof(struct HRRD_struct));
                  HRD.TmpRayEnergyLine = myrealloc_movable(HRD.TmpRayEnergyLine, HRD.TmpMaxNumRays * HRD.RayNumLines * sizeof(double));
                  HRD.TmpRayEnergyLineInit = myrealloc_movable(HRD.TmpRayEnergyLineInit, HRD.TmpMaxNumRays * HRD.RayNumLines * sizeof(double));
                  HRD.TmpRayEnergy = myrealloc_movable(HRD.TmpRayEnergy, HRD.TmpMaxNumRays * HRD.RayNumBins * sizeof(float));
                }

              for(j = 0; j < 4; j++)
                {
                  memcpy(HRD.TmpHRRD + HRD.TmpNumRays, HRD.HRRD + i, sizeof(struct HRRD_struct));
                  memcpy(HRD.TmpRayEnergyLine + HRD.TmpNumRays * HRD.RayNumLines, HRD.RayEnergyLine + i * HRD.RayNumLines, HRD.RayNumLines * sizeof(double));
                  memcpy(HRD.TmpRayEnergyLineInit + HRD.TmpNumRays * HRD.RayNumLines, HRD.RayEnergyLineInit + i * HRD.RayNumLines, HRD.RayNumLines * sizeof(double));
                  memcpy(HRD.TmpRayEnergy + HRD.TmpNumRays * HRD.RayNumBins, HRD.RayEnergy + i * HRD.RayNumBins, HRD.RayNumBins * sizeof(float));

                  HRD.TmpHRRD[HRD.TmpNumRays].lvl++;
                  HRD.TmpHRRD[HRD.TmpNumRays].hpidx = (HRD.HRRD[i].hpidx << 2) + j;

                  if(j > 0)
                    {
                      HRD.TmpHRRD[HRD.TmpNumRays].id = -1;
                      HRD.TmpHRRD[HRD.TmpNumRays].num_advances = 0;
                      HRD.TmpHRRD[HRD.TmpNumRays].num_comms = 0;
                    }

                  pix2vec_nest(1 << HRD.TmpHRRD[HRD.TmpNumRays].lvl, HRD.TmpHRRD[HRD.TmpNumRays].hpidx, dir);

                  for(k = 0; k < 3; k++)
                    {
                      HRD.TmpHRRD[HRD.TmpNumRays].dir[k] = 0;

                      for(l = 0; l < 3; l++)
                        HRD.TmpHRRD[HRD.TmpNumRays].dir[k] += HRD.HRRD[i].rot[k * 3 + l] * dir[l];

                      HRD.TmpHRRD[HRD.TmpNumRays].pos_old[k] = HRD.TmpHRRD[HRD.TmpNumRays].pos[k];
                      HRD.TmpHRRD[HRD.TmpNumRays].pos[k] = HRD.HRRD[i].pos[k] + HRD.HRRD[i].len * (HRD.TmpHRRD[HRD.TmpNumRays].dir[k] - HRD.HRRD[i].dir[k]);
                    }

                  HRD.TmpHRRD[HRD.TmpNumRays].energy = HRD.HRRD[i].energy / 4.;
                  HRD.TmpHRRD[HRD.TmpNumRays].energy_init = HRD.HRRD[i].energy_init / 4.;

                  for(k = 0; k < HRD.RayNumLines; k++)
                    {
                      HRD.TmpRayEnergyLine[HRD.TmpNumRays * HRD.RayNumLines + k] = HRD.RayEnergyLine[i * HRD.RayNumLines + k] / 4.;
                      HRD.TmpRayEnergyLineInit[HRD.TmpNumRays * HRD.RayNumLines + k] = HRD.RayEnergyLineInit[i * HRD.RayNumLines + k] / 4.;
                    }

                  for(k = 0; k < HRD.RayNumBins; k++)
                    HRD.TmpRayEnergy[HRD.TmpNumRays * HRD.RayNumBins + k] = HRD.RayEnergy[i * HRD.RayNumBins + k] / 4.;

                  HRD.TmpNumRays++;
                }

              if(i != HRD.NumRays - 1)
                {
                  memcpy(HRD.HRRD + i, HRD.HRRD + HRD.NumRays - 1, sizeof(struct HRRD_struct));
                  memcpy(HRD.RayEnergyLine + i * HRD.RayNumLines, HRD.RayEnergyLine + (HRD.NumRays - 1) * HRD.RayNumLines, HRD.RayNumLines * sizeof(double));
                  memcpy(HRD.RayEnergyLineInit + i * HRD.RayNumLines, HRD.RayEnergyLineInit + (HRD.NumRays - 1) * HRD.RayNumLines, HRD.RayNumLines * sizeof(double));
                  memcpy(HRD.RayEnergy + i * HRD.RayNumBins, HRD.RayEnergy + (HRD.NumRays - 1) * HRD.RayNumBins, HRD.RayNumBins * sizeof(float));
                }

              HRD.NumRays--;

              i--;
            }
        }

      MPI_Allreduce(&flag, &glob_flag, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      healray_split_comm();

      if(!HRD.RayMode)
        {
          HRD.NumEscFrac = 0;

          HRD.MaxNumEscFrac = HRD.TmpNumRays;

          HRD.EscFrac = mymalloc_movable(&HRD.EscFrac, "HRD.EscFrac", HRD.MaxNumEscFrac * sizeof(struct EscFrac_struct));
        }

      for(i = 0; i < HRD.TmpNumRays; i++)
        {
          rmray_flag = 0;

          idx = HRD.TmpHRRD[i].idx;

          if(HRD.TmpHRRD[i].pos[0] < 0 || HRD.TmpHRRD[i].pos[0] >= boxSize_X)
            rmray_flag = 1;

          if(HRD.TmpHRRD[i].pos[1] < 0 || HRD.TmpHRRD[i].pos[1] >= boxSize_Y)
            rmray_flag = 1;

          if(HRD.TmpHRRD[i].pos[2] < 0 || HRD.TmpHRRD[i].pos[2] >= boxSize_Z)
            rmray_flag = 1;

          if(dabs(P[idx].Pos[0] - HRD.TmpHRRD[i].pos_old[0]) > boxHalf_X)
            rmray_flag = 1;

          if(dabs(P[idx].Pos[1] - HRD.TmpHRRD[i].pos_old[1]) > boxHalf_Y)
            rmray_flag = 1;

          if(dabs(P[idx].Pos[2] - HRD.TmpHRRD[i].pos_old[2]) > boxHalf_Z)
            rmray_flag = 1;

          if(rmray_flag)
            {
              if(!HRD.RayMode)
                {
                  HRD.TotEnergyEsc += HRD.TmpHRRD[i].energy;

                  HRD.EscFrac[HRD.NumEscFrac].task = HRD.TmpHRRD[i].task_orig;
                  HRD.EscFrac[HRD.NumEscFrac].idx = HRD.TmpHRRD[i].idx_orig;
                  HRD.EscFrac[HRD.NumEscFrac].energy_frac = HRD.TmpHRRD[i].energy / HRD.TmpHRRD[i].energy_init_tot;

                  HRD.NumEscFrac++;
                }

              if(i != HRD.TmpNumRays - 1)
                {
                  memcpy(HRD.TmpHRRD + i, HRD.TmpHRRD + HRD.TmpNumRays - 1, sizeof(struct HRRD_struct));
                  memcpy(HRD.TmpRayEnergyLine + i * HRD.RayNumLines, HRD.TmpRayEnergyLine + (HRD.TmpNumRays - 1) * HRD.RayNumLines, HRD.RayNumLines * sizeof(double));
                  memcpy(HRD.TmpRayEnergyLineInit + i * HRD.RayNumLines, HRD.TmpRayEnergyLineInit + (HRD.TmpNumRays - 1) * HRD.RayNumLines, HRD.RayNumLines * sizeof(double));
                  memcpy(HRD.TmpRayEnergy + i * HRD.RayNumBins, HRD.TmpRayEnergy + (HRD.TmpNumRays - 1) * HRD.RayNumBins, HRD.RayNumBins * sizeof(float));
                }

              i--;

              HRD.TmpNumRays--;
            }
        }

      if(!HRD.RayMode)
        {
          mpi_distribute_items_to_tasks(HRD.EscFrac, offsetof(struct EscFrac_struct, task), &HRD.NumEscFrac, &HRD.MaxNumEscFrac, sizeof(struct EscFrac_struct), TAG_DENS_A);

          for(i = 0; i < HRD.NumEscFrac; i++)
            {
              if(HRD.EscFrac[i].task != ThisTask)
                terminate("EscFrac tasks don't agree!");

              if(HRD.EscFrac[i].idx < 0 || HRD.EscFrac[i].idx >= NumGas)
                terminate("EscFrac idx wrong!");

              HRD.EnergyFrac[HRD.EscFrac[i].idx] += HRD.EscFrac[i].energy_frac;
            }

          myfree_movable(HRD.EscFrac);
        }

      if(HRD.TmpNumRays)
        {
          HRD.MaxNumRays = HRD.NumRays + HRD.TmpNumRays;

          HRD.HRRD = myrealloc_movable(HRD.HRRD, HRD.MaxNumRays * sizeof(struct HRRD_struct));
          HRD.RayEnergyLine = myrealloc_movable(HRD.RayEnergyLine, HRD.MaxNumRays * HRD.RayNumLines * sizeof(double));
          HRD.RayEnergyLineInit = myrealloc_movable(HRD.RayEnergyLineInit, HRD.MaxNumRays * HRD.RayNumLines * sizeof(double));
          HRD.RayEnergy = myrealloc_movable(HRD.RayEnergy, HRD.MaxNumRays * HRD.RayNumBins * sizeof(float));

          memcpy(HRD.HRRD + HRD.NumRays, HRD.TmpHRRD, HRD.TmpNumRays * sizeof(struct HRRD_struct));
          memcpy(HRD.RayEnergyLine + HRD.NumRays * HRD.RayNumLines, HRD.TmpRayEnergyLine, HRD.TmpNumRays * HRD.RayNumLines * sizeof(double));
          memcpy(HRD.RayEnergyLineInit + HRD.NumRays * HRD.RayNumLines, HRD.TmpRayEnergyLineInit, HRD.TmpNumRays * HRD.RayNumLines * sizeof(double));
          memcpy(HRD.RayEnergy + HRD.NumRays * HRD.RayNumBins, HRD.TmpRayEnergy, HRD.TmpNumRays * HRD.RayNumBins * sizeof(float));

          HRD.NumRays = HRD.MaxNumRays;
        }

      new_num_rays = HRD.NumRays + 3 * HRD.TmpNumRays;

      MPI_Allreduce(&new_num_rays, &tot_new_num_rays, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

      HRD.TotInitNumRays += tot_new_num_rays;

      myfree_movable(HRD.TmpRayEnergy);
      myfree_movable(HRD.TmpRayEnergyLineInit);
      myfree_movable(HRD.TmpRayEnergyLine);
      myfree_movable(HRD.TmpHRRD);
    }
}


void healray_advance()
{
  int i, tot_nrays;
  double t0, t1;
  point *DP = Mesh.DP;

  t0 = second();

  for(i = tot_nrays = 0; i < HRD.RayNumThreads; i++)
    {
      HRD.ThreadNumRays[i] = HRD.NumRays / HRD.RayNumThreads;

      if(i < HRD.NumRays % HRD.RayNumThreads)
        HRD.ThreadNumRays[i] += 1;

      HRD.ThreadPrevNumRays[i] = tot_nrays;

      tot_nrays += HRD.ThreadNumRays[i];
    }

  if(tot_nrays != HRD.NumRays)
    terminate("HEALRAY: Ray numbers don't agree!");

  if(!HRD.RayMode)
    {
      HRD.NumEscFrac = 0;

      HRD.MaxNumEscFrac = HRD.NumRays;

      HRD.EscFrac = mymalloc_movable(&HRD.EscFrac, "HRD.EscFrac", HRD.MaxNumEscFrac * sizeof(struct EscFrac_struct));
    }

  if(0)
    healray_mic_test();

#pragma omp parallel private(i)
  {
    int j, tid, idx, last_idx, next_edge, rmray_flag;
    long hpidx;
    double r, len, nh, alpha, denergy;

    double tau, exp_tau, one_minus_exp_tau, energy_ave, source, heat_rate;
#ifdef TGCHEM
    int k, l, tidx, temp_idx, temp_idx_orig;
    double abh2, abhii, fac, mu, temp, dtemp;
    double pos[3], dr[3], dr2, dv[3], vrel, shift_temp, shift_vrel;
    double new_energy, tot_emiss, TotEmiss, TotCross;

    double *Emiss = malloc(HRD.RayNumLines * sizeof(double));
    double *Cross = malloc(HRD.RayNumLines * sizeof(double));
    double *ReducExp = malloc(HRD.RayNumLines * sizeof(double));
    double *NuDExp = malloc(HRD.RayNumFreq * sizeof(double));
    double *EnergyReducFac = malloc(HRD.RayNumBins * sizeof(double));
    double *DEnergyTable = malloc(HRD.RayNumBins * sizeof(double));
#endif
    tid = omp_get_thread_num();

    for(i = HRD.ThreadPrevNumRays[tid]; i < HRD.ThreadPrevNumRays[tid] + HRD.ThreadNumRays[tid]; i++)
      {
        idx = HRD.HRRD[i].idx;

        next_edge = find_next_voronoi_cell(&Mesh, idx, HRD.HRRD[i].pos, HRD.HRRD[i].dir, HRD.HRRD[i].idx_prev, &len);

        if(next_edge < 0 || len < 0)
          {
            printf("HEALRAY: Ray ID = %ld, next_edge = %d, len = %g\n", HRD.HRRD[i].id, next_edge, len);

            terminate("HEALRAY: Something is wrong in find_next_voronoi_cell!");
          }

        if(HRD.RayMode)
          {
            vec2pix_nest(1 << HRD.RayLvlInit, HRD.HRRD[i].dir, &hpidx);

            HRD.PartRayCount[idx * HRD.NumRaysInit + hpidx]++;
          }

        r = get_cell_radius(idx);

        nh = HRD.RhoToNHCgs * SphP[idx].Density;

        if(HRD.RayMode)
          {
            source = 1.;

            healray_get_alpha_and_source(idx, &alpha, &source);

            tau = alpha * len * HRD.LenToCgs;

            if(tau)
              {
                exp_tau = exp(-tau);

                one_minus_exp_tau = (1. - exp_tau);

                energy_ave = HRD.HRRD[i].energy * one_minus_exp_tau / tau + source;     // * (tau - one_minus_exp_tau) / tau;

                HRD.HRRD[i].energy = HRD.HRRD[i].energy * exp_tau + source * one_minus_exp_tau;

                heat_rate = energy_ave * alpha;
#pragma omp atomic
                HRD.Energy[idx] += tau * heat_rate;
#pragma omp atomic
                HRD.EnergyAux[idx] += tau;
              }
          }
        else
          {
#ifdef TGCHEM
            if(!TGCD.ChemMode)
              {
                if(HRD.RayDebugMode)
                  {
                    temp_idx_orig = HRD.HRRD[i].temp_idx_orig;

                    dtemp = HRD.HRRD[i].temp_init - TGCD.TempTable[temp_idx_orig];

                    TotEmiss = TGCD.H2TotEmissTable[temp_idx_orig] + dtemp * TGCD.H2DTotEmissTable[temp_idx_orig];

                    tot_emiss = 0;

                    if(HRD.RayMultiFreqMode)
                      {
                        for(j = 0; j < HRD.RayNumLines; j++)
                          {
                            tidx = TGCD.H2LineIdxTable[j * TGCHEM_NUM_TEMP + temp_idx_orig] * TGCHEM_NUM_TEMP + temp_idx_orig;

                            Emiss[j] = TGCD.H2EmissTable[tidx] + dtemp * TGCD.H2DEmissTable[tidx];

                            tot_emiss += Emiss[j];
                          }

                        for(j = 0; j < HRD.RayNumLines; j++)
                          Emiss[j] *= TotEmiss / tot_emiss;
                      }
                  }

                abh2 = SphP[idx].Abund[1];
                abhii = SphP[idx].Abund[2];

                fac = abh2 * nh * HRD.LenToCgs * len;

                mu = (1. + 4. * HE_ABUND) / (1. + HE_ABUND - abh2 + abhii);

                temp = (SphP[idx].Gamma - 1.) * mu * HRD.UThermToTemp * SphP[idx].Utherm;

                temp_idx = imin((int) (log10(temp / TGCHEM_TEMP_MIN) / TGCHEM_LOG_DTEMP), TGCHEM_NUM_TEMP - 1);

                dtemp = temp - TGCD.TempTable[temp_idx];

                TotCross = TGCD.H2TotCrossTable[temp_idx] + dtemp * TGCD.H2DTotCrossTable[temp_idx];

                if(HRD.RayMultiFreqMode)
                  {
                    for(j = 0; j < HRD.RayNumLines; j++)
                      {
                        tidx = TGCD.H2LineIdxTable[j * TGCHEM_NUM_TEMP + HRD.HRRD[i].temp_idx_orig] * TGCHEM_NUM_TEMP + temp_idx;

                        Cross[j] = TGCD.H2CrossTable[tidx] + dtemp * TGCD.H2DCrossTable[tidx];
                      }

                    for(j = 0; j < HRD.RayNumLines; j++)
                      ReducExp[j] = -fac * Cross[j];

                    for(j = dr2 = 0; j < 3; j++)
                      {
                        dr[j] = P[idx].Pos[j] - HRD.HRRD[i].pos_init[j];
                        dv[j] = P[idx].Vel[j] - HRD.HRRD[i].vel_init[j];

                        dr2 += dr[j] * dr[j];
                      }

                    for(j = vrel = 0; j < 3; j++)
                      {
                        if(dr2)
                          dr[j] /= sqrt(dr2);

                        vrel += dv[j] * dr[j];
                      }

                    vrel *= HRD.VelToCgs;

                    if(HRD.RayMultiFreqMode == 3 || HRD.RayMultiFreqMode == 4)
                      shift_temp = 1;
                    else
                      shift_temp = sqrt(HRD.HRRD[i].temp_init / temp);

                    if(HRD.RayMultiFreqMode == 2 || HRD.RayMultiFreqMode == 4)
                      shift_vrel = 0;
                    else
                      shift_vrel = vrel / sqrt(BOLTZMANN * temp / PROTONMASS);

                    for(j = 0; j < HRD.RayNumFreq; j++)
                      NuDExp[j] = HRD.NuDPos[j] * shift_temp - shift_vrel;

                    for(j = 0; j < HRD.RayNumFreq; j++)
                      NuDExp[j] = exp(-NuDExp[j] * NuDExp[j]);

                    for(j = l = 0; j < HRD.RayNumLines; j++)
                      for(k = 0; k < HRD.RayNumFreq; k++)
                        EnergyReducFac[l++] = exp(ReducExp[j] * NuDExp[k]);

                    for(j = 0; j < HRD.RayNumBins; j++)
                      DEnergyTable[j] = HRD.RayEnergy[i * HRD.RayNumBins + j] * (1. - EnergyReducFac[j]);

                    for(j = 0; j < HRD.RayNumBins; j++)
                      HRD.RayEnergy[i * HRD.RayNumBins + j] -= DEnergyTable[j];

                    for(j = denergy = 0; j < HRD.RayNumBins; j++)
                      denergy += DEnergyTable[j];

                    if(HRD.RayDebugMode)
                      {
                        printf("reduc: %g %g %g %g %g %g %d %g %d %g %g %g %g\n", nh, abh2, len, HRD.LenToCgs, fac, HRD.HRRD[i].temp_init,
                               temp_idx_orig, temp, temp_idx, dtemp, shift_temp, vrel, shift_vrel);

                        for(j = 0; j < HRD.RayNumLines; j++)
                          {
                            tidx = TGCD.H2LineIdxTable[j * TGCHEM_NUM_TEMP + temp_idx_orig] * TGCHEM_NUM_TEMP + temp_idx_orig;

                            printf("reduc line: %d %d %g %g %g %g %g %g %g\n", j, TGCD.H2LineIdxTable[j * TGCHEM_NUM_TEMP + temp_idx_orig], TotEmiss,
                                   TGCD.H2EmissTable[tidx], Emiss[j], TotCross, TGCD.H2CrossTable[tidx], Cross[j], ReducExp[j]);
                          }

                        for(j = 0; j < HRD.RayNumFreq; j++)
                          printf("reduc nu: %d %g %g %g\n", j, HRD.NuDPos[j], HRD.NuDPos[j] * shift_temp - shift_vrel, NuDExp[j]);

                        for(j = l = 0; j < HRD.RayNumLines; j++)
                          for(k = 0; k < HRD.RayNumFreq; k++)
                            {
                              printf("reduc bin: %d %d %d %g %g\n", l, j, k, HRD.NuDFac[k], EnergyReducFac[l]);

                              l++;
                            }

                        for(j = l = 0; j < HRD.RayNumLines; j++)
                          {
                            HRD.RayEnergyLine[i * HRD.RayNumLines + j] = 0.;

                            for(k = 0; k < HRD.RayNumFreq; k++)
                              HRD.RayEnergyLine[i * HRD.RayNumLines + j] += HRD.RayEnergy[i * HRD.RayNumBins + l++];

                            printf("reduc absorption: %d %g\n", j, HRD.RayEnergyLine[i * HRD.RayNumLines + j] / HRD.RayEnergyLineInit[i * HRD.RayNumLines + j]);
                          }
                      }

                    for(j = new_energy = 0; j < HRD.RayNumBins; j++)
                      new_energy += HRD.RayEnergy[i * HRD.RayNumBins + j];

                    denergy = dmax(HRD.HRRD[i].energy - new_energy, 0.);
                  }
                else
                  denergy = HRD.HRRD[i].energy * (1. - exp(-fac * TotCross));
              }
            else
              denergy = 0;
#else
            denergy = HRD.HRRD[i].energy * (1. - exp(-HRD.Alpha * len * HRD.LenToCgs));
#endif
            HRD.HRRD[i].energy -= denergy;
#pragma omp atomic
            HRD.Energy[idx] += denergy;
          }

        rmray_flag = 0;

        if(HRD.HRRD[i].energy <= HRD.RayMinEgyFrac * HRD.HRRD[i].energy_init)
          rmray_flag = 1;
        else
          {
            if(HRD.HRRD[i].pos[0] < 0. || HRD.HRRD[i].pos[0] >= boxSize_X)
              rmray_flag = 1;

            if(HRD.HRRD[i].pos[1] < 0. || HRD.HRRD[i].pos[1] >= boxSize_Y)
              rmray_flag = 1;

            if(HRD.HRRD[i].pos[2] < 0. || HRD.HRRD[i].pos[2] >= boxSize_Z)
              rmray_flag = 1;
          }

        healray_add_outray(tid, i);

        if(HRD.RayDebugMode && rmray_flag)
          {
#ifdef TGCHEM
            if(!TGCD.ChemMode)
              {
                printf("absorption: %d %d %g %g %g %g %g %g %g %g %d %g %g %g %g %g %g %g\n", HRD.HRRD[i].idx_orig, HRD.HRRD[i].idx, nh,
                       HRD.LenToCgs * len, fac, abh2, abhii, mu, temp, dtemp, temp_idx, TotEmiss, TotCross, fac * TotCross, HRD.HRRD[i].energy_init,
                       denergy, HRD.HRRD[i].energy, HRD.HRRD[i].energy / HRD.HRRD[i].energy_init);

                if(HRD.RayMultiFreqMode)
                  for(j = 0; j < HRD.RayNumBins; j++)
                    printf("absorption line: %d %g\n", j, HRD.RayEnergy[i * HRD.RayNumBins + j] / HRD.HRRD[i].energy_init);
              }
#else
            printf("absorption: %d %g %g %g %g %g\n", HRD.HRRD[i].idx_orig, HRD.Alpha, HRD.LenToCgs * len, HRD.HRRD[i].energy_init, denergy, HRD.HRRD[i].energy / HRD.HRRD[i].energy_init);
#endif
          }

        if(rmray_flag)
          {
            if(!HRD.RayMode)
              {
#pragma omp atomic
                HRD.TotEnergyEsc += HRD.HRRD[i].energy;

                HRD.EscFrac[HRD.NumEscFrac].task = HRD.HRRD[i].task_orig;
                HRD.EscFrac[HRD.NumEscFrac].idx = HRD.HRRD[i].idx_orig;
                HRD.EscFrac[HRD.NumEscFrac].energy_frac = HRD.HRRD[i].energy / HRD.HRRD[i].energy_init_tot;
#pragma omp atomic
                HRD.NumEscFrac++;
              }

            if(i != HRD.ThreadPrevNumRays[tid] + HRD.ThreadNumRays[tid] - 1)
              {
                last_idx = HRD.ThreadPrevNumRays[tid] + HRD.ThreadNumRays[tid] - 1;

                memcpy(HRD.HRRD + i, HRD.HRRD + last_idx, sizeof(struct HRRD_struct));
                memcpy(HRD.RayEnergyLine + i * HRD.RayNumLines, HRD.RayEnergyLine + last_idx * HRD.RayNumLines, HRD.RayNumLines * sizeof(double));
                memcpy(HRD.RayEnergyLineInit + i * HRD.RayNumLines, HRD.RayEnergyLineInit + last_idx * HRD.RayNumLines, HRD.RayNumLines * sizeof(double));
                memcpy(HRD.RayEnergy + i * HRD.RayNumBins, HRD.RayEnergy + last_idx * HRD.RayNumBins, HRD.RayNumBins * sizeof(float));
              }

            HRD.ThreadNumRays[tid]--;
#pragma omp atomic
            HRD.LocNumAdvances += HRD.HRRD[i].num_advances + 1;
#pragma omp atomic
            HRD.LocNumComms += HRD.HRRD[i].num_comms;

            i--;
          }
        else
          {
            HRD.HRRD[i].iter++;
            HRD.HRRD[i].task = DC[next_edge].task;

            if(HRD.HRRD[i].task != ThisTask)
              {
                HRD.HRRD[i].idx_prev = -1;
                HRD.HRRD[i].num_comms++;
              }
            else
              HRD.HRRD[i].idx_prev = HRD.HRRD[i].idx;

            HRD.HRRD[i].idx = DC[next_edge].index;
            HRD.HRRD[i].num_advances++;

            for(j = 0; j < 3; j++)
              {
                HRD.HRRD[i].pos_old[j] = HRD.HRRD[i].pos[j];
                HRD.HRRD[i].pos[j] += len * HRD.HRRD[i].dir[j];
              }

            HRD.HRRD[i].len += len;

            if(!HRD.RaySplitFac && HRD.HRRD[i].task == ThisTask)
              i--;
          }
#pragma omp atomic
        HRD.LocNumOps += HRD.RayNumBins;
      }
#ifdef TGCHEM
    free(DEnergyTable);
    free(EnergyReducFac);
    free(NuDExp);
    free(ReducExp);
    free(Cross);
    free(Emiss);
#endif
  }

  for(i = HRD.NumRays = 0; i < HRD.RayNumThreads; i++)
    {
      memmove(HRD.HRRD + HRD.NumRays, HRD.HRRD + HRD.ThreadPrevNumRays[i], HRD.ThreadNumRays[i] * sizeof(struct HRRD_struct));
      memmove(HRD.RayEnergyLine + HRD.NumRays * HRD.RayNumLines, HRD.RayEnergyLine + HRD.ThreadPrevNumRays[i] * HRD.RayNumLines, HRD.ThreadNumRays[i] * HRD.RayNumLines * sizeof(double));
      memmove(HRD.RayEnergyLineInit + HRD.NumRays * HRD.RayNumLines, HRD.RayEnergyLineInit + HRD.ThreadPrevNumRays[i] * HRD.RayNumLines, HRD.ThreadNumRays[i] * HRD.RayNumLines * sizeof(double));
      memmove(HRD.RayEnergy + HRD.NumRays * HRD.RayNumBins, HRD.RayEnergy + HRD.ThreadPrevNumRays[i] * HRD.RayNumBins, HRD.ThreadNumRays[i] * HRD.RayNumBins * sizeof(float));

      HRD.NumRays += HRD.ThreadNumRays[i];
    }

  t1 = second();

  HRD.DtAdvance += timediff(t0, t1);

  t0 = second();

  MPI_Barrier(MPI_COMM_WORLD);

  t1 = second();

  HRD.DtImba += timediff(t0, t1);

  t0 = second();

  healray_comm(0);

  if(!HRD.RayMode)
    mpi_distribute_items_to_tasks(HRD.EscFrac, offsetof(struct EscFrac_struct, task), &HRD.NumEscFrac, &HRD.MaxNumEscFrac, sizeof(struct EscFrac_struct), TAG_DENS_A);

  t1 = second();

  HRD.DtComm += timediff(t0, t1);

  if(!HRD.RayMode)
    {
      for(i = 0; i < HRD.NumEscFrac; i++)
        HRD.EnergyFrac[HRD.EscFrac[i].idx] += HRD.EscFrac[i].energy_frac;

      myfree_movable(HRD.EscFrac);
    }
}


void healray_get_alpha_and_source(int idx, double *alpha, double *source)
{
#ifdef TGCHEM
  int temp_idx;
  double nh, abh2, abhii, mu, temp, dtemp, cross, emiss;

  if(!TGCD.ChemMode)
    {
      nh = HRD.RhoToNHCgs * SphP[idx].Density;

      abh2 = SphP[idx].Abund[1];
      abhii = SphP[idx].Abund[2];

      mu = (1. + 4. * HE_ABUND) / (1. + HE_ABUND - abh2 + abhii);

      temp = (SphP[idx].Gamma - 1.) * mu * HRD.UThermToTemp * SphP[idx].Utherm;

      temp_idx = imin((int) (log10(temp / TGCHEM_TEMP_MIN) / TGCHEM_LOG_DTEMP), TGCHEM_NUM_TEMP - 1);

      dtemp = temp - TGCD.TempTable[temp_idx];

      cross = TGCD.H2TotCrossTable[temp_idx] + dtemp * TGCD.H2DTotCrossTable[temp_idx];

      *alpha = abh2 * nh * cross;

      emiss = TGCD.H2TotEmissTable[temp_idx] + dtemp * TGCD.H2DTotEmissTable[temp_idx];

      *source = emiss / cross;
    }
  else
    {
      *alpha = 0.;

      *source = 0.;
    }
#else
  *alpha = HRD.Alpha;

  *source = 1.;
#endif
}


void healray_add_outray(int tid, int i)
{
  int j, k, idx;

  for(j = 0; j < HRD.RayNumTrace; j++)
    if(HRD.HRRD[i].id == HRD.TraceList[j])
      {
        if(HRD.HRRD[i].iter + 1 > HRD.TraceNumIter)
          terminate("HEALRAY: HRD.TraceNumIter not chosen large enough!");

        idx = HRD.HRRD[i].iter * HRD.RayNumTrace + j;

        *((double *) HRD.TraceEnergyFrac[tid] + idx) = HRD.HRRD[i].energy / HRD.HRRD[i].energy_init;
#ifdef TGCHEM
        if(!TGCD.ChemMode && HRD.RayMultiFreqMode)
          for(k = 0; k < HRD.RayNumBins; k++)
            {
              idx = HRD.HRRD[i].iter * HRD.RayNumTrace * HRD.RayNumBins + j * HRD.RayNumBins + k;

              *((double *) HRD.TraceEnergyFracFreq[tid] + idx) += HRD.RayEnergy[i * HRD.RayNumBins + k] / HRD.HRRD[i].energy_init;
            }
#endif
      }
}
