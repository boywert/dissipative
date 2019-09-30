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

void healray_alloc();
void healray_init_mesh();
void healray_init_variables();
void healray_init_emitters();
void healray_init_emitters_alt();
void healray_find_intersection(double *l0, double *dir, double *sec);
int healray_idx_list_compare(const void *a, const void *b);

struct idx_list_struct
{
  int idx;
  double rnd;
} *idx_list;



void healray_begrun()
{
  CPU_Step[CPU_MISC] += measure_time();

  MPI_Bcast(&HRD, sizeof(struct HRD_struct), MPI_BYTE, 0, MPI_COMM_WORLD);

  if(!strncmp(HRD.RaySourceFile, "-1", MAXLEN_PATH))
    HRD.SourceFlag = -1;
  else
    {
#if defined (REFLECTIVE_X) || defined (REFLECTIVE_Y) || defined (REFLECTVIE_Z)
      terminate("HEALRAY currently does not work with reflective boundary conditions!");
#endif

#ifndef VORONOI_DYNAMIC_UPDATE
      terminate("HEALRAY requires VORONOI_DYNAMIC_UPDATE!");
#endif
      if(!strncmp(HRD.RaySourceFile, "0", MAXLEN_PATH))
        HRD.SourceFlag = 0;
      else
        HRD.SourceFlag = 1;

      HRD.TraceNumIter = 1000;

      HRD.RayNumLines = imax(1, HRD.RayNumLines);
      HRD.RayNumFreq = imax(1, HRD.RayNumFreq);

      if(HRD.RayNumFreq == 1)
        HRD.RayMultiFreqMode = 0;

      HRD.RayNumBins = HRD.RayNumLines * HRD.RayNumFreq;
    }

  CPU_Step[CPU_HEALRAY] += measure_time();
}


void healray_init()
{
  CPU_Step[CPU_MISC] += measure_time();

  if(HRD.SourceFlag > 0)
    healray_read_sources();

  CPU_Step[CPU_HEALRAY] += measure_time();
}


void healray_init_step()
{
  double t0, t1;
#ifndef FORCE_EQUAL_TIMESTEPS
  healray_init_mesh();
#endif
  t0 = second();

  healray_alloc();

  healray_init_variables();

  if(HRD.SourceFlag)
    healray_init_sources();
  else
    {
      if(HRD.RayMode)
        healray_init_emitters_alt();
      else
        healray_init_emitters();
    }

  healray_init_rayout();

  t1 = second();

  HRD.DtInit += timediff(t0, t1);
}


void healray_init_mesh()
{
  int i;

  mpi_printf("HEALRAY: preparing mesh for ray-tracing\n\n");

  free_mesh();

  HRD.SaveTimeBin = mymalloc_movable(&HRD.SaveTimeBin, "HRD.SaveTimeBin", NumPart * sizeof(short int));

  for(i = 0; i < NumPart; i++)
    {
      HRD.SaveTimeBin[i] = P[i].TimeBin;
      P[i].TimeBin = 0;
    }

  for(i = 0; i < TIMEBINS; i++)
    {
      HRD.SaveTimeBinActive[i] = TimeBinSynchronized[i];
      TimeBinSynchronized[i] = 0;
    }

  for(i = 0; i < NumGas; i++)
    P[i].TimeBin = TIMEBINS - 1;

  TimeBinSynchronized[TIMEBINS - 1] = 1;

  CPU_Step[CPU_HEALRAY] += measure_time();

  reconstruct_timebins();

  create_mesh();

  mesh_setup_exchange();
}


void healray_alloc()
{
  int i;

  HRD.ThreadNumRays = mymalloc_movable(&HRD.ThreadNumRays, "HRD.ThreadNumRays", HRD.RayNumThreads * sizeof(int));
  HRD.ThreadPrevNumRays = mymalloc_movable(&HRD.ThreadPrevNumRays, "HRD.ThreadPrevNumRays", HRD.RayNumThreads * sizeof(int));

  HRD.NuDPos = mymalloc_movable(&HRD.NuDPos, "HRD.NuDPos", HRD.RayNumFreq * sizeof(double));
  HRD.NuDFac = mymalloc_movable(&HRD.NuDFac, "HRD.NuDFac", HRD.RayNumFreq * sizeof(double));

  if(HRD.RayNumTrace)
    {
      HRD.TraceList = mymalloc_movable(&HRD.TraceList, "HRD.TraceList", HRD.RayNumTrace * sizeof(long));
      HRD.TraceEnergyFrac = mymalloc_movable(&HRD.TraceEnergyFrac, "HRD.TraceEnergyFrac", HRD.RayNumThreads * sizeof(void *));

      for(i = 0; i < HRD.RayNumThreads; i++)
        {
          HRD.TraceEnergyFrac[i] = mymalloc_movable(&HRD.TraceEnergyFrac[i], "HRD.TraceEnergyFrac[i]", HRD.TraceNumIter * HRD.RayNumTrace * sizeof(double));
        }
#ifdef TGCHEM
      if(!TGCD.ChemMode && HRD.RayMultiFreqMode)
        {
          HRD.TraceEnergyFracFreq = mymalloc_movable(&HRD.TraceEnergyFracFreq, "HRD.TraceEnergyFracFreq", HRD.RayNumThreads * sizeof(void *));

          for(i = 0; i < HRD.RayNumThreads; i++)
            HRD.TraceEnergyFracFreq[i] = mymalloc_movable(&HRD.TraceEnergyFracFreq[i], "HRD.TraceEnergyFracFreq[i]", HRD.TraceNumIter * HRD.RayNumTrace * HRD.RayNumBins * sizeof(double));
        }
#endif
    }

  if(HRD.RayNumSources)
    {
      HRD.SourceIDList = mymalloc_movable(&HRD.SourceIDList, "HRD.SourceIDList", HRD.RayNumSources * sizeof(int));
      HRD.SourcePosList = mymalloc_movable(&HRD.SourcePosList, "HRD.SourcePosList", 3 * HRD.RayNumSources * sizeof(double));
      HRD.SourceEgyList = mymalloc_movable(&HRD.SourceEgyList, "HRD.SourceEgyList", HRD.RayNumSources * sizeof(double));
    }

  if(HRD.RayMode)
    HRD.IdxList = mymalloc_movable(&HRD.IdxList, "HRD.IdxList", NumGas * sizeof(int));
  else
    HRD.EnergyInit = mymalloc_movable(&HRD.EnergyInit, "HRD.EnergyInit", NumGas * sizeof(double));

  HRD.Energy = mymalloc_movable(&HRD.Energy, "HRD.Energy", NumGas * sizeof(double));

  if(HRD.RayMode)
    HRD.EnergyAux = mymalloc_movable(&HRD.EnergyAux, "HRD.EnergyAux", NumGas * sizeof(double));
  else
    HRD.EnergyFrac = mymalloc_movable(&HRD.EnergyFrac, "HRD.EnergyFrac", NumGas * sizeof(double));
}


void healray_init_variables()
{
  int i;
  double tot_nud = 0.;

  omp_set_num_threads(HRD.RayNumThreads);

  set_cosmo_factors_for_current_time();

  HRD.NumRaysInit = 12 * (1 << 2 * HRD.RayLvlInit);
  //HRD.NumRaysInit = 1;

  HRD.EnergyUnit = 1.;
  HRD.UThermToTemp = PROTONMASS / BOLTZMANN * All.UnitEnergy_in_cgs / All.UnitMass_in_g;
  HRD.LenToCgs = All.cf_atime / All.HubbleParam * All.UnitLength_in_cm;
  HRD.VolToCgs = pow(HRD.LenToCgs, 3);
  HRD.RhoToNHCgs = HYDROGEN_MASSFRAC / PROTONMASS * All.UnitDensity_in_cgs * All.cf_a3inv * pow(All.HubbleParam, 2);
  HRD.VelToCgs = All.UnitVelocity_in_cm_per_s / All.cf_atime;

  HRD.InitEnergy = 1. / HRD.EnergyUnit;
  HRD.Alpha = 5. * All.BoxSize / All.cf_atime * All.HubbleParam / All.UnitLength_in_cm;

  HRD.NuDRange = 3.;
  HRD.NuDWidth = (2. * HRD.NuDRange / HRD.RayNumFreq);

  for(i = 0; i < HRD.RayNumFreq; i++)
    {
      HRD.NuDPos[i] = HRD.NuDRange * ((2. * i + 1.) / HRD.RayNumFreq - 1.);
      HRD.NuDFac[i] = exp(-pow(HRD.NuDPos[i], 2.));

      tot_nud += HRD.NuDFac[i];
    }

  for(i = 0; i < HRD.RayNumFreq; i++)
    HRD.NuDFac[i] /= tot_nud;

  if(!HRD.RayMode)
    memset(HRD.EnergyInit, -1, NumGas * sizeof(double));

  memset(HRD.Energy, 0., NumGas * sizeof(double));

  if(HRD.RayMode)
    memset(HRD.EnergyAux, 0., NumGas * sizeof(double));
  else
    memset(HRD.EnergyFrac, -1, NumGas * sizeof(double));

  if(!HRD.RayMode)
    {
      HRD.TotEnergyInit = 0.;
      HRD.TotEnergyEsc = 0.;
    }

  HRD.TotInitNumRays = 0;

  HRD.LocNumOps = 0;
  HRD.LocNumAdvances = 0;
  HRD.LocNumComms = 0;

  HRD.DtInit = 0.;
  HRD.DtAdvance = 0.;
  HRD.DtImba = 0.;
  HRD.DtComm = 0.;
}


void healray_init_emitters()
{
  char buf[200];
  int i, j, k, max_size, flag_redo, tot_num_gas;
  int nrays_per_emitter, num_emitters, max_num_emitters, num_iter, max_num_iter;
  double nh, rnd_min, rnd_max, *source_pos_list;
  double energy, emitted_energy, tot_emitted_energy, injected_energy, tot_injected_energy;
  FILE *file;
#ifdef TGCHEM
  int temp_idx;
  double fac, abh2, abhii, mu, temp, dtemp, H2TotEmiss;
#endif

  //HRD.NumRaysInit = 1;

  MPI_Allreduce(&NumGas, &tot_num_gas, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(HRD.RaySplitFac)
    nrays_per_emitter = HRD.RaySplitFac * tot_num_gas / NTask;
  else
    nrays_per_emitter = HRD.NumRaysInit;

  max_size = sizeof(struct HRRD_struct) + 2 * HRD.RayNumLines * sizeof(double) + HRD.RayNumBins * sizeof(float);

  max_num_emitters = (double) FreeBytes / max_size / nrays_per_emitter / HRD_MEM_SAFETY_FAC;

  MPI_Allreduce(&max_num_emitters, &HRD.MaxNumEmitters, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

  HRD.MaxNumEmitters = 5000;

  if(!HRD.MaxNumEmitters)
    terminate("HEALRAY: Not enough memory even for a single emitter!");

  HRD.MaxNumRays = HRD.MaxNumEmitters * HRD.NumRaysInit;

  HRD.HRRD = mymalloc_movable(&HRD.HRRD, "HRD.HRRD", HRD.MaxNumRays * sizeof(struct HRRD_struct));
  HRD.RayEnergyLine = mymalloc_movable(&HRD.RayEnergyLine, "HRD.RayEnergyLine", HRD.MaxNumRays * HRD.RayNumLines * sizeof(double));
  HRD.RayEnergyLineInit = mymalloc_movable(&HRD.RayEnergyLineInit, "HRD.RayEnergyLineInit", HRD.MaxNumRays * HRD.RayNumLines * sizeof(double));
  HRD.RayEnergy = mymalloc_movable(&HRD.RayEnergy, "HRD.RayEnergy", HRD.MaxNumRays * HRD.RayNumBins * sizeof(float));

  HRD.NumRays = 0;

  num_emitters = 0;

  emitted_energy = 0.;
  injected_energy = 0.;

  if(HRD.RayNumSources)
    {
      rnd_min = -1.;
      rnd_max = 1.;

      source_pos_list = mymalloc_movable(&source_pos_list, "source_pos_list", 3 * HRD.RayNumSources * sizeof(double));

      memset(source_pos_list, 0., 3 * HRD.RayNumSources * sizeof(double));

      if(ThisTask == 0)
        {
          for(i = 0; i < HRD.RayNumSources; i++)
            {
              do
                {
                  flag_redo = 0;

                  HRD.SourceIDList[i] = imax(get_random_number() * tot_num_gas, 1);
                  //HRD.SourceIDList[i] = 5662285;

                  HRD.SourceEgyList[i] = pow(10., rnd_min + get_random_number() * (rnd_max - rnd_min)) * HRD.InitEnergy;

                  for(j = 0; j < i; j++)
                    if(HRD.SourceIDList[i] == HRD.SourceIDList[j])
                      flag_redo = 1;
                }
              while(flag_redo);
            }

          sprintf(buf, "sator/data/healray_sources.dat");

          if(!(file = fopen(buf, "w")))
            terminate("HEALRAY: Could not open file!\n");

          fwrite(&HRD.RayNumSources, sizeof(int), 1, file);
          fwrite(&HRD.InitEnergy, sizeof(double), 1, file);
          fwrite(&HRD.EnergyUnit, sizeof(double), 1, file);
          fwrite(&HRD.Alpha, sizeof(double), 1, file);
          fwrite(HRD.SourceIDList, sizeof(int), HRD.RayNumSources, file);
          fwrite(HRD.SourceEgyList, sizeof(double), HRD.RayNumSources, file);

          fclose(file);
        }

      MPI_Bcast(HRD.SourceIDList, HRD.RayNumSources, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(HRD.SourceEgyList, HRD.RayNumSources, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

  for(i = 0; i < NumGas; i++)
    {
      if(P[i].ID == 0 && P[i].Mass == 0)
        continue;

      if(HRD.RayNumSources)
        {
          for(j = 0; j < HRD.RayNumSources; j++)
            if(P[i].ID == HRD.SourceIDList[j])
              {
                for(k = 0; k < 3; k++)
                  source_pos_list[3 * j + k] = P[i].Pos[k];

                HRD.EnergyInit[i] = HRD.SourceEgyList[j];

                num_emitters++;

                HRD.NumRays += HRD.NumRaysInit;
              }
        }
      else
        {
          nh = HRD.RhoToNHCgs * SphP[i].Density;

          //nh = get_random_number();
          //if(P[i].ID == 1023449619)
          //if(nh > 9e10)
          //if(nh < 0.1)
          //if(ThisTask == TGD.NHMaxTask && i == TGD.NHMaxIdx)
          //if(ThisTask == 180 && i == 5074)
          {
#ifdef TGCHEM
            if(!TGCD.ChemMode)
              {
                abh2 = SphP[i].Abund[1];
                abhii = SphP[i].Abund[2];

                fac = abh2 * nh * SphP[i].Volume * HRD.VolToCgs / HRD.EnergyUnit;

                mu = (1. + 4. * HE_ABUND) / (1. + HE_ABUND - abh2 + abhii);

                temp = (SphP[i].Gamma - 1.) * mu * HRD.UThermToTemp * SphP[i].Utherm;

                temp_idx = imin((int) (log10(temp / TGCHEM_TEMP_MIN) / TGCHEM_LOG_DTEMP), TGCHEM_NUM_TEMP - 1);

                dtemp = temp - TGCD.TempTable[temp_idx];

                H2TotEmiss = TGCD.H2TotEmissTable[temp_idx] + dtemp * TGCD.H2DTotEmissTable[temp_idx];

                energy = fac * H2TotEmiss;
              }
            else
              energy = 1.;
#else
            energy = 1.;
#endif
            emitted_energy += energy;

            if(get_random_number() < HRD.RayPartFrac)
              {
                HRD.EnergyInit[i] = energy;

                injected_energy += energy;

                num_emitters++;

                HRD.NumRays += HRD.NumRaysInit;
              }
          }
        }
    }

  if(HRD.RayNumSources)
    {
      MPI_Allreduce(source_pos_list, HRD.SourcePosList, 3 * HRD.RayNumSources, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      myfree_movable(source_pos_list);
    }
  else
    {
      MPI_Allreduce(&emitted_energy, &tot_emitted_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&injected_energy, &tot_injected_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      if(tot_injected_energy)
        HRD.EnergyIncreaseFac = tot_emitted_energy / tot_injected_energy;
      else
        HRD.EnergyIncreaseFac = 0;
    }

  MPI_Allreduce(&num_emitters, &HRD.TotNumEmitters, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(num_emitters)
    {
      num_iter = num_emitters / HRD.MaxNumEmitters;

      if(num_emitters > HRD.MaxNumEmitters && num_emitters % HRD.MaxNumEmitters)
        num_iter++;
    }
  else
    num_iter = 0;

  MPI_Allreduce(&num_iter, &max_num_iter, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  if(max_num_iter > 1)
    mpi_printf("HEALRAY: Warning: %d extra loops needed!\n", max_num_iter - 1);
}


void healray_init_emitters_alt()
{
  int i, tot_num_gas;

  MPI_Allreduce(&NumGas, &tot_num_gas, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  HRD.LimitNumRays = imax((int) (2 * pow(tot_num_gas, 2. / 3.) / NTask) * HRD.NumRaysInit, HRD.NumRaysInit);

  HRD.HRRD = mymalloc_movable(&HRD.HRRD, "HRD.HRRD", HRD.LimitNumRays * sizeof(struct HRRD_struct));
  HRD.RayEnergyLine = mymalloc_movable(&HRD.RayEnergyLine, "HRD.RayEnergyLine", HRD.LimitNumRays * HRD.RayNumLines * sizeof(double));
  HRD.RayEnergyLineInit = mymalloc_movable(&HRD.RayEnergyLineInit, "HRD.RayEnergyLineInit", HRD.LimitNumRays * HRD.RayNumLines * sizeof(double));
  HRD.RayEnergy = mymalloc_movable(&HRD.RayEnergy, "HRD.RayEnergy", HRD.LimitNumRays * HRD.RayNumBins * sizeof(float));
  HRD.PartRayCount = mymalloc_movable(&HRD.PartRayCount, "HRD.PartRayCount", NumGas * HRD.NumRaysInit * sizeof(int));

  memset(HRD.PartRayCount, 0, NumGas * HRD.NumRaysInit * sizeof(float));

  idx_list = mymalloc("idx_list", NumGas * sizeof(struct idx_list_struct));

  for(i = 0; i < NumGas; i++)
    {
      idx_list[i].idx = i;

      idx_list[i].rnd = get_random_number();
    }

  mysort(idx_list, NumGas, sizeof(struct idx_list_struct), healray_idx_list_compare);

  for(i = 0; i < NumGas; i++)
    HRD.IdxList[i] = idx_list[i].idx;

  myfree(idx_list);

  HRD.PIdx = 0;
}


void healray_init_rays()
{
  int i, j, k, l, num_emitters, tot_num_emitters;
  long long init_num_rays;
  double phi, theta, psi, dir[3], energy_init, tot_energy_init;
#ifdef TGCHEM
  int m, tidx, temp_idx;
  double nh, fac, abh2, abhii, mu, temp, dtemp;
  double tot_emiss, tot_emiss_test, H2TotEmiss;
  double *Emiss;
#endif

  num_emitters = 0;

  HRD.NumRays = 0;

#ifdef TGCHEM
  Emiss = mymalloc_movable(&Emiss, "Emiss", HRD.RayNumLines * sizeof(double));
#endif
  energy_init = 0.;

  for(i = 0; i < NumGas; i++)
    if(HRD.EnergyInit[i] >= 0)
      {
        phi = get_random_number() * 2. * M_PI;
        theta = get_random_number() * 2. * M_PI;
        psi = get_random_number() * 2. * M_PI;

        for(j = 0; j < HRD.NumRaysInit; j++)
          {
            HRD.HRRD[HRD.NumRays].iter = 0;
            HRD.HRRD[HRD.NumRays].task = ThisTask;
            HRD.HRRD[HRD.NumRays].task_orig = HRD.HRRD[HRD.NumRays].task;
            HRD.HRRD[HRD.NumRays].idx = i;
            HRD.HRRD[HRD.NumRays].idx_prev = -1;
            HRD.HRRD[HRD.NumRays].idx_orig = HRD.HRRD[HRD.NumRays].idx;
            HRD.HRRD[HRD.NumRays].lvl = HRD.RayLvlInit;
            HRD.HRRD[HRD.NumRays].hpidx = j;
            HRD.HRRD[HRD.NumRays].num_advances = 0;
            HRD.HRRD[HRD.NumRays].num_comms = 0;
            HRD.HRRD[HRD.NumRays].id = (long) P[i].ID * HRD.NumRaysInit + j;
            HRD.HRRD[HRD.NumRays].len = 0.;

            for(k = 0; k < 3; k++)
              {
                HRD.HRRD[HRD.NumRays].pos[k] = P[i].Pos[k];
                HRD.HRRD[HRD.NumRays].pos_old[k] = HRD.HRRD[HRD.NumRays].pos[k];
              }

            HRD.HRRD[HRD.NumRays].rot[0] = cos(theta) * cos(psi);
            HRD.HRRD[HRD.NumRays].rot[1] = cos(phi) * sin(psi) + sin(phi) * sin(theta) * cos(psi);
            HRD.HRRD[HRD.NumRays].rot[2] = sin(phi) * sin(psi) - cos(phi) * sin(theta) * cos(psi);
            HRD.HRRD[HRD.NumRays].rot[3] = -cos(theta) * sin(psi);
            HRD.HRRD[HRD.NumRays].rot[4] = cos(phi) * cos(psi) - sin(phi) * sin(theta) * sin(psi);
            HRD.HRRD[HRD.NumRays].rot[5] = sin(phi) * cos(psi) + cos(phi) * sin(theta) * sin(psi);
            HRD.HRRD[HRD.NumRays].rot[6] = sin(theta);
            HRD.HRRD[HRD.NumRays].rot[7] = -sin(phi) * cos(theta);
            HRD.HRRD[HRD.NumRays].rot[8] = cos(phi) * cos(theta);

            pix2vec_nest(1 << HRD.RayLvlInit, j, dir);

            for(k = 0; k < 3; k++)
              {
                HRD.HRRD[HRD.NumRays].dir[k] = 0.;

                for(l = 0; l < 3; l++)
                  HRD.HRRD[HRD.NumRays].dir[k] += HRD.HRRD[HRD.NumRays].rot[k * 3 + l] * dir[l];
              }
#ifdef TGCHEM
            if(!TGCD.ChemMode)
              {
                for(k = 0; k < 3; k++)
                  {
                    HRD.HRRD[HRD.NumRays].pos_init[k] = P[i].Pos[k];
                    HRD.HRRD[HRD.NumRays].vel_init[k] = P[i].Vel[k];
                  }

                nh = HRD.RhoToNHCgs * SphP[i].Density;

                abh2 = SphP[i].Abund[1];
                abhii = SphP[i].Abund[2];

                fac = abh2 * nh * SphP[i].Volume * HRD.VolToCgs / HRD.NumRaysInit / HRD.EnergyUnit;

                mu = (1. + 4. * HE_ABUND) / (1. + HE_ABUND - abh2 + abhii);

                temp = (SphP[i].Gamma - 1.) * mu * HRD.UThermToTemp * SphP[i].Utherm;

                HRD.HRRD[HRD.NumRays].temp_init = temp;

                temp_idx = imin((int) (log10(temp / TGCHEM_TEMP_MIN) / TGCHEM_LOG_DTEMP), TGCHEM_NUM_TEMP - 1);

                HRD.HRRD[HRD.NumRays].temp_idx_orig = temp_idx;

                dtemp = temp - TGCD.TempTable[temp_idx];

                H2TotEmiss = TGCD.H2TotEmissTable[temp_idx] + dtemp * TGCD.H2DTotEmissTable[temp_idx];

                tot_emiss = tot_emiss_test = 0;

                if(HRD.RayMultiFreqMode)
                  {
                    for(k = 0; k < HRD.RayNumLines; k++)
                      {
                        tidx = TGCD.H2LineIdxTable[k * TGCHEM_NUM_TEMP + temp_idx] * TGCHEM_NUM_TEMP + temp_idx;

                        Emiss[k] = TGCD.H2EmissTable[tidx] + dtemp * TGCD.H2DEmissTable[tidx];
                      }

                    for(k = 0; k < HRD.RayNumLines; k++)
                      tot_emiss += Emiss[k];

                    for(k = 0; k < HRD.RayNumLines; k++)
                      Emiss[k] *= H2TotEmiss / tot_emiss;

                    for(k = m = 0; k < HRD.RayNumLines; k++)
                      for(l = 0; l < HRD.RayNumFreq; l++)
                        HRD.RayEnergy[HRD.NumRays * HRD.RayNumBins + m++] = HRD.EnergyIncreaseFac * fac * Emiss[k] * HRD.NuDFac[l];

                    if(HRD.RayDebugMode)
                      for(k = m = 0; k < HRD.RayNumLines; k++)
                        {
                          tot_emiss_test += Emiss[k];

                          HRD.RayEnergyLineInit[HRD.NumRays * HRD.RayNumLines + k] = 0.;

                          for(l = 0; l < HRD.RayNumFreq; l++)
                            HRD.RayEnergyLineInit[HRD.NumRays * HRD.RayNumLines + k] += HRD.RayEnergy[HRD.NumRays * HRD.RayNumBins + m++];

                          printf("line emiss: %d %d %d %g %g\n", k, temp_idx, TGCD.H2LineIdxTable[k * TGCHEM_NUM_TEMP + temp_idx],
                                 TGCD.H2EmissTable[TGCD.H2LineIdxTable[k * TGCHEM_NUM_TEMP + temp_idx] * TGCHEM_NUM_TEMP + temp_idx], HRD.RayEnergyLineInit[HRD.NumRays * HRD.RayNumLines + k] / fac);
                        }
                  }

                HRD.HRRD[HRD.NumRays].energy = HRD.EnergyIncreaseFac * fac * H2TotEmiss;

                if(HRD.RayDebugMode)
                  printf("cell: %g %g %g %g %g %g %g %g %d %g %g %g %g %g\n", nh, abh2, abhii, SphP[i].Volume, HRD.VolToCgs, fac, mu, temp, temp_idx,
                         dtemp, H2TotEmiss, tot_emiss, tot_emiss_test, HRD.HRRD[HRD.NumRays].energy);
              }
            else
              HRD.HRRD[HRD.NumRays].energy = 1.;
#else
            HRD.HRRD[HRD.NumRays].energy = HRD.EnergyInit[i] / HRD.NumRaysInit;
#endif
            HRD.HRRD[HRD.NumRays].energy_init = HRD.HRRD[HRD.NumRays].energy;
            HRD.HRRD[HRD.NumRays].energy_init_tot = HRD.HRRD[HRD.NumRays].energy_init * HRD.NumRaysInit;

            energy_init += HRD.HRRD[HRD.NumRays].energy;

            HRD.NumRays++;
          }

        HRD.EnergyInit[i] = -1;
        HRD.EnergyFrac[i] = 0.;

        num_emitters++;

        if(num_emitters == HRD.MaxNumEmitters)
          break;
      }
#ifdef TGCHEM
  myfree(Emiss);
#endif
  sumup_large_ints(1, &HRD.NumRays, &init_num_rays);

  HRD.TotInitNumRays += init_num_rays;

  MPI_Allreduce(&num_emitters, &tot_num_emitters, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  HRD.TotNumEmitters -= tot_num_emitters;

  MPI_Allreduce(&energy_init, &tot_energy_init, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  HRD.TotEnergyInit += tot_energy_init;
}


void healray_init_rays_alt()
{
  int i, j, k, pidx, num_emitters;
  long long init_num_rays;
  double phi, theta, psi, dir[3], rot[9];

  HRD.MaxNumRays = HRD.LimitNumRays;

  HRD.HRRD = myrealloc_movable(HRD.HRRD, HRD.MaxNumRays * sizeof(struct HRRD_struct));
  HRD.RayEnergyLine = myrealloc_movable(HRD.RayEnergyLine, HRD.MaxNumRays * HRD.RayNumLines * sizeof(double));
  HRD.RayEnergyLineInit = myrealloc_movable(HRD.RayEnergyLineInit, HRD.MaxNumRays * HRD.RayNumLines * sizeof(double));
  HRD.RayEnergy = myrealloc_movable(HRD.RayEnergy, HRD.MaxNumRays * HRD.RayNumBins * sizeof(float));

  HRD.NumRays = 0;

  int max_part = NumGas;

  if(ThisTask < NTask)
    while(HRD.NumRays + HRD.NumRaysInit <= HRD.MaxNumRays && HRD.PIdx < max_part)
      {
        if(!(P[pidx].ID == 0 && P[pidx].Mass == 0.))
          {
            phi = get_random_number() * 2. * M_PI;
            theta = get_random_number() * 2. * M_PI;
            psi = get_random_number() * 2. * M_PI;

            pidx = HRD.IdxList[HRD.PIdx];

            rot[0] = cos(theta) * cos(psi);
            rot[1] = cos(phi) * sin(psi) + sin(phi) * sin(theta) * cos(psi);
            rot[2] = sin(phi) * sin(psi) - cos(phi) * sin(theta) * cos(psi);
            rot[3] = -cos(theta) * sin(psi);
            rot[4] = cos(phi) * cos(psi) - sin(phi) * sin(theta) * sin(psi);
            rot[5] = sin(phi) * cos(psi) + cos(phi) * sin(theta) * sin(psi);
            rot[6] = sin(theta);
            rot[7] = -sin(phi) * cos(theta);
            rot[8] = cos(phi) * cos(theta);

            for(i = 0; i < HRD.NumRaysInit; i++)
              if(!HRD.PartRayCount[pidx * HRD.NumRaysInit + i])
                {
                  HRD.HRRD[HRD.NumRays].iter = 0;
                  HRD.HRRD[HRD.NumRays].task = ThisTask;
                  HRD.HRRD[HRD.NumRays].task_orig = HRD.HRRD[HRD.NumRays].task;
                  HRD.HRRD[HRD.NumRays].idx = pidx;
                  HRD.HRRD[HRD.NumRays].idx_prev = -1;
                  HRD.HRRD[HRD.NumRays].lvl = HRD.RayLvlInit;
                  HRD.HRRD[HRD.NumRays].num_advances = 0;
                  HRD.HRRD[HRD.NumRays].num_comms = 0;
                  HRD.HRRD[HRD.NumRays].id = (long) P[pidx].ID * HRD.NumRaysInit + i;
                  HRD.HRRD[HRD.NumRays].len = 0.;

                  for(j = 0; j < 9; j++)
                    HRD.HRRD[HRD.NumRays].rot[j] = rot[j];

                  pix2vec_nest(1 << HRD.RayLvlInit, i, dir);

                  for(j = 0; j < 3; j++)
                    {
                      HRD.HRRD[HRD.NumRays].dir[j] = 0.;

                      for(k = 0; k < 3; k++)
                        HRD.HRRD[HRD.NumRays].dir[j] += HRD.HRRD[HRD.NumRays].rot[j * 3 + k] * dir[k];
                    }

                  vec2pix_nest(1 << HRD.RayLvlInit, HRD.HRRD[HRD.NumRays].dir, &HRD.HRRD[HRD.NumRays].hpidx);

                  healray_find_intersection(P[pidx].Pos, HRD.HRRD[HRD.NumRays].dir, HRD.HRRD[HRD.NumRays].pos);

                  for(j = 0; j < 3; j++)
                    HRD.HRRD[HRD.NumRays].pos_old[j] = HRD.HRRD[HRD.NumRays].pos[j];

                  HRD.HRRD[HRD.NumRays].energy = 0.;

                  HRD.NumRays++;
                }
          }

        HRD.PIdx++;
      }

  healray_comm(1);

  sumup_large_ints(1, &HRD.NumRays, &init_num_rays);

  HRD.TotInitNumRays += init_num_rays;

  num_emitters = max_part - HRD.PIdx;

  MPI_Allreduce(&num_emitters, &HRD.TotNumEmitters, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
}


void healray_find_intersection(double *l0, double *dir, double *sec)
{
  int i, j, axis;
  double p0[3], n[3], l_dot_n, p_dot_n, d, min_d;

  min_d = MAX_DOUBLE_NUMBER;

  for(i = 0; i < 6; i++)
    {
      axis = i / 2;

      for(j = 0; j < 3; j++)
        p0[j] = n[j] = 0;

      p0[axis] = (i % 2) * All.BoxSize;
      n[axis] = All.BoxSize;

      l_dot_n = 0.;

      for(j = 0; j < 3; j++)
        l_dot_n += -(*(dir + j)) * n[j];

      p_dot_n = 0;

      for(j = 0; j < 3; j++)
        p_dot_n += (p0[j] - *(l0 + j)) * n[j];

      d = p_dot_n / l_dot_n;

      if(d >= 0 && d < min_d)
        {
          min_d = d;

          for(j = 0; j < 3; j++)
            sec[j] = l0[j] + d * -(*(dir + j));
        }
    }

  for(j = 0; j < 3; j++)
    sec[j] = dmin(dmax(sec[j], 0.), (1. - 1e-10) * All.BoxSize);
}


int healray_idx_list_compare(const void *a, const void *b)
{
  if(((struct idx_list_struct *) a)->rnd < ((struct idx_list_struct *) b)->rnd)
    return -1;

  if(((struct idx_list_struct *) a)->rnd > ((struct idx_list_struct *) b)->rnd)
    return +1;

  return 0;
}
