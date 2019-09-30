/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        
 * \date        10/2014
 * \author      Christine Simpson
 * \brief
 * \details     Injects thermal & kinetic feedback in sphere
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <math.h>
#include "../allvars.h"
#include "../proto.h"

#if defined(LOCAL_FEEDBACK) && defined(EXTERNALSHEARBOX)

static int calc;                /* determines phase of calculation */
static int find_feedback_cells_evaluate(int target, int mode, int thread_id);

/* communication structures */
typedef struct
{
  MyDouble Pos[3];
  MyDouble Vel[3];
  MyDouble Volume;

  MyDouble SNEnergy;
  MyDouble fkin;
  MyDouble SNMassReturn;
  MyDouble deltap;
  MyDouble total_vol;
  MyDouble h;

  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;


/* temporary result structures */

static struct temp_celldata
{
  char flag_feedback;

  MyDouble total_vol;

  MyDouble SNEnergy;
  MyDouble fkin;
  MyDouble SNMassReturn;
  MyDouble deltap;
  MyDouble h;

  MyDouble aterm1, b, c, tot_ekin_b;
  MyDouble resid[3];
  int Ncells;

  /* If feedback is centered on particles, record info about closest hydro cell */
#ifdef LOCAL_FEEDBACK_PARTICLES
  MyDouble NgbDistance;
  MyDouble NgbVolume;
#endif

}
 *TempCellData;



/* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in * in, int i, int firstnode)
{
  int k;

#ifdef LOCAL_FEEDBACK_PARTICLES
  for(k = 0; k < 3; k++)
    in->Pos[k] = P[i].Pos[k];

  if(calc == 0)
    if(TempCellData[i].Ncells == 0)
      in->Volume = 0.0;
    else
      in->Volume = TempCellData[i].NgbVolume;

#else
  for(k = 0; k < 3; k++)
    in->Pos[k] = SphP[i].Center[k];

  in->Volume = SphP[i].Volume;
#endif
  for(k = 0; k < 3; k++)
    in->Vel[k] = P[i].Vel[k];


  in->total_vol = TempCellData[i].total_vol;
  in->h = TempCellData[i].h;
  if(calc > 0)
    {
      in->SNEnergy = TempCellData[i].SNEnergy;
      in->fkin = TempCellData[i].fkin;
      in->SNMassReturn = TempCellData[i].SNMassReturn;
      in->deltap = TempCellData[i].deltap;
    }
  in->Firstnode = firstnode;
}


typedef struct
{
  MyDouble total_vol;
  MyDouble aterm1, b, c, tot_ekin_b;
  MyDouble resid[3];
  int Ncells;

#ifdef LOCAL_FEEDBACK_PARTICLES
  MyDouble NgbDistance;
  MyDouble NgbVolume;
  MyDouble NgbRadius;
#endif

} data_out;

static data_out *DataResult, *DataOut;

/* routine to store or combine result data */
static void out2particle(data_out * out, int i, int mode)
{
  int k;
  if(mode == MODE_LOCAL_PARTICLES)      /* initial store */
    {
      if(calc == 0)
        {
          TempCellData[i].total_vol = out->total_vol;
          TempCellData[i].Ncells = out->Ncells;
        }
      if(calc == -1)
        {
#ifdef LOCAL_FEEDBACK_PARTICLES
          TempCellData[i].NgbDistance = out->NgbDistance;
          TempCellData[i].NgbVolume = out->NgbVolume;
#endif
        }
      if(calc == 1)
        {
          TempCellData[i].aterm1 = out->aterm1;
          TempCellData[i].b = out->b;
          TempCellData[i].c = out->c;
          TempCellData[i].tot_ekin_b = out->tot_ekin_b;
        }

      if(calc == 2)
        {
          for(k = 0; k < 3; k++)
            TempCellData[i].resid[k] = out->resid[k];
        }
    }
  else                          /* combine */
    {

      if(calc == 0)
        {
          TempCellData[i].total_vol += out->total_vol;
          TempCellData[i].Ncells += out->Ncells;

        }
      if(calc == -1)
#ifdef LOCAL_FEEDBACK_PARTICLES
        if(out->NgbDistance < TempCellData[i].NgbDistance)
          {
            TempCellData[i].NgbDistance = out->NgbDistance;
            TempCellData[i].NgbVolume = out->NgbVolume;
            TempCellData[i].h = fmax(3.0 * pow(3.0 * out->NgbVolume / 4.0 / M_PI, 1.0 / 3.0), TempCellData[i].h);
          }
#endif

      if(calc == 1)
        {
          TempCellData[i].aterm1 += out->aterm1;
          TempCellData[i].b += out->b;
          TempCellData[i].c += out->c;
          TempCellData[i].tot_ekin_b += out->tot_ekin_b;
        }
      if(calc == 2)
        {
          for(k = 0; k < 3; k++)
            TempCellData[i].resid[k] += out->resid[k];
        }
    }
}

#include "../generic_comm_helpers2.h"

#ifdef LOCAL_FEEDBACK_PARTICLES

void create_shear_box_particles(void)
{
  CPU_Step[CPU_MISC] += measure_time();

  //mpi_printf("CMS_FEEDBACK: Starting Starformation and Feedback Routine %d\n",TimeBinsHydro.NActiveParticles);

  int i;
  double pi = 3.14159265358979323846;
  double rho, sfr, dz, dt, prob, sech, deltaM;
  int count_spawned = 0, count_converted = 0, count_feedback = 0;

  double SNEnergy = All.LocalFeedbackSNEnergy * 1e51 / All.UnitVelocity_in_cm_per_s / All.UnitVelocity_in_cm_per_s / All.UnitMass_in_g;
  double SNRate = All.LocalFeedbackSNRate;


  double StarParticleCreationTime;
  double MinimumStarMass = 0.5 * All.TargetGasMass;

  double Sigma0 = All.ShearBoxSigma0;
  double fg = All.ShearBoxFg;
  double mu = All.ShearBoxMu;
  double bscale = 61.0 * (fg / 0.1 / mu) / (Sigma0 / 10.0);

#ifndef LOCAL_FEEDBACK_KSRATE
  double SFEff = All.LocalFeedbackSFEff;
  double SFDenThresh = All.LocalFeedbackSFDenThresh * PROTONMASS * All.UnitLength_in_cm * All.UnitLength_in_cm * All.UnitLength_in_cm / All.UnitMass_in_g;
  double density_threshold = SFDenThresh;
#else
  double sigma_star = 2.5e-10 * pow(Sigma0, 1.4);
#endif

  double SNTimeDelay = All.LocalFeedbackSNTimeDelay;
  double SNTimeSpread = All.LocalFeedbackSNTimeSpread;

#ifndef LOCAL_FEEDBACK_KSRATE
  if(All.Time < SNTimeDelay)
    density_threshold = 0.0;
#endif

  int idx;
  /* loop over all active cells and particles */
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      /* just consider gas cells */
      if(P[i].Type == 0)
        {

          /* skip cells that have been swallowed or eliminated */
          if(P[i].Mass == 0 && P[i].ID == 0)
            continue;

          dz = P[i].Pos[2] - boxHalf_Z;
          sech = 1.0 / cosh(dz / bscale);
          rho = Sigma0 * sech * sech / 2.0 / bscale;

          dt = (All.HighestActiveTimeBin ? (((int) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;

#ifdef LOCAL_FEEDBACK_KSRATE    /* calculate sfr based on the KS relation */

          sfr = sigma_star * P[i].Mass / Sigma0;
          dt *= All.UnitTime_in_Megayears * 1e6;
          deltaM = sfr * dt;

#else /* calculate sfr based on local dynamical time */

          if(SphP[i].Density < density_threshold)
            continue;

          double tdyn, tot_rho;

          tot_rho = SphP[i].Density + (1.0 / fg - 1.0) * rho;
          tdyn = sqrt(3.0 * M_PI / 32.0 / All.G / tot_rho);
          sfr = SFEff * P[i].Mass / tdyn;
          deltaM = sfr * dt;

#endif
          prob = deltaM / MinimumStarMass;
          SphP[i].Sfr = sfr;
          double rand = get_random_number();
          //printf("CMS_FEEDBACK: prob = %f; dt = %f\n",prob,dt);
          if(rand < prob)
            {
              if(P[i].Mass < 2.0 * MinimumStarMass)
                {
                  convert_cell_into_star(i, All.Time);
                  timebin_remove_particle(&TimeBinsHydro, idx, P[i].TimeBinHydro);
                  printf("CMS_FEEDBACK: converting cell to star; time = %f; id = %d; sfr = %f; prob = %f; mass = %f; height = %f\n", All.Time, P[i].ID, sfr, prob, P[i].Mass, dz);
                  Stars_converted++;
                  count_converted++;
                }
              else
                {
                  printf("CMS_FEEDBACK: spawning star from cell; Cell: time = %f; id = %d; sfr = %f; prob = %f; mass = %f; height = %f\n", All.Time, P[i].ID, sfr, prob, P[i].Mass, dz);
                  spawn_star_from_cell(i, All.Time, NumPart + count_spawned, MinimumStarMass);
                  count_spawned++;
                }
            }
        }
    }

  /* Calculate ids of spawned particles and print out number of created particles */
  int countall_spawned, countall_converted;
  MPI_Allreduce(&count_spawned, &countall_spawned, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&count_converted, &countall_converted, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(countall_spawned > 0)
    {
      int *list;

      if(All.MaxID == 0)        /* MaxID not calculated yet */
        calculate_maxid();

      list = mymalloc("list", NTask * sizeof(int));

      MPI_Allgather(&count_spawned, 1, MPI_INT, list, 1, MPI_INT, MPI_COMM_WORLD);

      MyIDType newid = All.MaxID + 1;

      for(i = 0; i < ThisTask; i++)
        newid += list[i];

      myfree(list);

      for(i = 0; i < count_spawned; i++)
        {
          P[NumPart + i].ID = newid++;
          printf("CMS_FEEDBACK: spawning star from cell; Star: time = %f; id = %d; mass = %f; height = %f\n", All.Time, P[NumPart + i].ID, P[NumPart + i].Mass, abs(P[NumPart + i].Pos[2] - 5000.0));
        }
      All.MaxID += countall_spawned;

    }

  if(countall_spawned > 0 || countall_converted > 0)
    {
      All.TotNumPart += countall_spawned;
      All.TotNumGas -= countall_converted;
      NumPart += count_spawned;
    }

  int countall = countall_converted + countall_spawned;

  if(countall > 0)
    {
      //    ngb_recompute_nodes();
      mpi_printf("CMS_FEEDBACK: Created %d particles at time %f\n", countall, All.Time);
    }

}
#endif





static void kernel_local(void)
{
  int i, idx;
#ifdef GENERIC_ASYNC
  int flag = 0;
#endif

#pragma omp parallel private(i, idx)
  {
    int j, threadid = get_thread_num();
#ifdef GENERIC_ASYNC
    int count = 0;
#endif

    for(j = 0; j < NTask; j++)
      Thread[threadid].Exportflag[j] = -1;

    while(1)
      {
        if(Thread[threadid].ExportSpace < MinSpace)
          break;

#ifdef GENERIC_ASYNC
        if(threadid == 0)
          {
            if((count & POLLINGINTERVAL) == 0)
              if(generic_polling_primary(count, TimeBinsGravity.NActiveParticles))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        idx = NextParticle++;
        if(idx >= TimeBinsGravity.NActiveParticles)
          break;

        i = TimeBinsGravity.ActiveParticleList[idx];
        if(i < 0)
          c ontinue;

        if(TempCellData[i].flag_feedback == 1)
          find_feedback_cells_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
      }
  }
}


static void kernel_imported(void)
{
  /* now do the particles that were sent to us */
  int i, cnt = 0;
#pragma omp parallel private(i)
  {
    int threadid = get_thread_num();
#ifdef GENERIC_ASYNC
    int count = 0;
#endif

    while(1)
      {
#pragma omp atomic capture
        i = cnt++;

        if(i >= Nimport)
          break;

#ifdef GENERIC_ASYNC
        if(threadid == 0)
          {
            if((count & POLLINGINTERVAL) == 0)
              generic_polling_secondary();
          }

        count++;
#endif

        //printf("CMS_FEEDBACK[%d]: id = %d; total_vol = %f\n",ThisTask,P[i].ID,TempCellData[i].total_vol);
        find_feedback_cells_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}


/* Routine that injects feedback either around previously created star particles, or if star particles
   are not used, instantaneously injected based on local gas properties.
   If instantaneous injection is used, the injection rate is computed either based on the KS rate or
   based on the local dynamical time.
   The total amount of energy injected per SN event, and the fraction of that energy that takes a kinetic 
   form, is set by run time parameters. */

void inject_shear_box_sphere_feedback(void)
{
  CPU_Step[CPU_MISC] += measure_time();

  mpi_printf("CMS_FEEDBACK: Starting Starformation and Feedback Routine\n");

  int i;
  double pi = 3.14159265358979323846;
  double rho, sfr, dz, dt, prob, sech, deltaM;

  double SNEnergy = All.LocalFeedbackSNEnergy * 1e51 / All.UnitVelocity_in_cm_per_s / All.UnitVelocity_in_cm_per_s / All.UnitMass_in_g;
  double SNRate = All.LocalFeedbackSNRate;

  double StarParticleCreationTime;
  double MinimumStarMass = 0.5 * All.TargetGasMass;

  double Sigma0 = All.ShearBoxSigma0;
  double fg = All.ShearBoxFg;
  double mu = All.ShearBoxMu;
  double bscale = 61.0 * (fg / 0.1 / mu) / (Sigma0 / 10.0);
  double sigma_star = 2.5e-10 * pow(Sigma0, 1.4);

#ifndef LOCAL_FEEDBACK_KSRATE
  double SFEff = All.LocalFeedbackSFEff;
  double SFDenThresh = All.LocalFeedbackSFDenThresh * PROTONMASS * All.UnitLength_in_cm * All.UnitLength_in_cm * All.UnitLength_in_cm / All.UnitMass_in_g;
  double density_threshold = SFDenThresh;
#endif

#ifdef LOCAL_FEEDBACK_PARTICLES
  double SNTimeDelay_max = All.LocalFeedbackSNTimeDelay;
  double SNTimeDelay = All.LocalFeedbackSNTimeDelay;

  double SNTimeSpread_max = All.LocalFeedbackSNTimeSpread;
  double SNTimeSpread = All.LocalFeedbackSNTimeSpread;
  double tform, Age, N_SN;
#endif

  TempCellData = mymalloc("TemCellData", NumPart * sizeof(struct temp_celldata));

  int count_spawned = 0, count_converted = 0, count_feedback = 0;

  /*** Let's first flag the cells or paritcles that inject a feedback event ***/
  /* loop over all active cells and particles */
  int idx;
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      TempCellData[i].flag_feedback = 0;
      TempCellData[i].Ncells = 0;
#ifdef LOCAL_FEEDBACK_PARTICLES /* Feedback events created from star particles that were previously created */
      if(P[i].Type == 4)
        {
          N_SN = P[i].Mass * SNRate / 100.0 - P[i].FeedbackDone;
          if(N_SN < 0.0)
            continue;

          tform = P[i].StellarAge;
          Age = All.Time - tform;

          dt = (All.HighestActiveTimeBin ? (((int) 1) << P[i].TimeBinGrav) : 0) * All.Timebase_interval;

          prob = 0.5 * (1.0 + erf((Age + dt - SNTimeDelay) / sqrt(2) / SNTimeSpread)) - 0.5 * (1.0 + erf((Age - SNTimeDelay) / sqrt(2) / SNTimeSpread));
          prob *= N_SN;

          /* force a feedback event when Age ~ SNTimeDelay */
          if(P[i].FeedbackDone == 0 && fabs(Age - SNTimeDelay) < 2.0 * dt)
            prob = 1.0;

          double rand = get_random_number();

          if(rand < prob)
            {
              TempCellData[i].flag_feedback = 1;
              P[i].FeedbackDone += 1;
              printf("CMS_FEEDBACK: Feedback particle ID = %d; time = %f; rand = %f; prob = %f; N_SN = %f; Age = %f; Age-SNTD = %f; |Age-SNTD| = %f; SNTD = %f\n",
                     P[i].ID, All.Time, rand, prob, N_SN, Age, Age - SNTimeDelay, fabs(Age - SNTimeDelay), All.LocalFeedbackSNTimeDelay);
              printf("CMS_FEEDBACK: center x,y,z,t: %f %f %f %f\n", P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], All.Time);
            }
        }

#else /* inject feedback without particles */

      /* just consider gas cells */
      if(P[i].Type == 0)
        {

          /* skip cells that have been swallowed or eliminated */
          if(P[i].Mass == 0 && P[i].ID == 0)
            continue;

          /* Compute sfr */

          dz = P[i].Pos[2] - boxHalf_Z;
          sech = 1.0 / cosh(dz / bscale);
          rho = Sigma0 * sech * sech / 2.0 / bscale;
          double rand = get_random_number();
          dt = (All.HighestActiveTimeBin ? (((int) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;

#ifdef LOCAL_FEEDBACK_KSRATE

          sfr = sigma_star * P[i].Mass / Sigma0;
          dt *= All.UnitTime_in_Megayears * 1e6;
          deltaM = sfr * dt;
          prob = deltaM * SNRate / 100.0;

#else /*use local dyn time to calc sfr */

          if(SphP[i].Density < density_threshold)
            continue;

          double tdyn, tot_rho;

          tot_rho = SphP[i].Density + (1.0 / fg - 1.0) * rho;
          tdyn = sqrt(3.0 * M_PI / 32.0 / All.G / tot_rho);
          sfr = SFEff * P[i].Mass / tdyn;
          deltaM = sfr * dt;
#endif
          prob = deltaM * SNRate / 100.0;

#ifdef USE_SFR
          SphP[i].Sfr = sfr;
#endif
          /* select a cell randomly with a certain rate */
          if(rand < prob)
            {
              TempCellData[i].flag_feedback = 1;
              printf("CMS_FEEDBACK: center x,y,z,t: %f %f %f %f\n", SphP[i].Center[0], SphP[i].Center[1], SphP[i].Center[2], All.Time);
            }
        }
#endif
    }


  /*** Now, let's go for each of the selected cells three times over all neighbours in a search sphere.
       on iteration -1: Find closest gas cell
       on iteration 0: Just determine total volume of the neighbours cells in the aperture
       on iteration 1: Compute momentum we will add to the gas
       on iteration 2: Inject energy, mass and momentum
  ***/

  generic_set_MaxNexport();

  int j, nexport, nimport, place, ngrp, recvTask, ndone_flag, ndone, dummy;
  int repeat_calc0 = 0;
  double cell_radius, h, h2;

  for(calc = -1; calc < 3; calc++)
    {
      //printf("CMS_FEEDBACK[%d]: starting calc %d\n",ThisTask,calc);
#ifndef LOCAL_FEEDBACK_PARTICLES
      if(calc == -1)
        continue;
#endif
      /* let's do some preparatory calculations */

      int idx;
      /* loop over all active cells and particles */
      for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
        {
          i = TimeBinsGravity.ActiveParticleList[idx];
          if(i < 0)
            continue;

          if(TempCellData[i].flag_feedback)
            {

#ifdef LOCAL_FEEDBACK_PARTICLES

              /*if(calc == 1 && repeat_calc0 == 0 && TempCellData[i].Ncells < 27){
                 repeat_calc0 += 1;
                 calc = 0;
                 } */

              /* on first pass, tie radius to softening length */
              if(calc == -1)
                cell_radius = 25.0 * All.MinimumComovingHydroSoftening;

              /* on subsequent passes, tie radius to resolution of nearest gas cell */
              if(calc > -1)
                cell_radius = pow(TempCellData[i].NgbVolume * 3.0 / 4.0 / pi, 1.0 / 3.0);

#else /*cell centered feedback */
              /* determine feedback radius */
              cell_radius = pow(SphP[i].Volume * 3.0 / 4.0 / pi, 1.0 / 3.0);
#endif

              TempCellData[i].h = 3.0 * cell_radius;
              h2 = TempCellData[i].h * TempCellData[i].h;
              printf("CMS_FEEDBACK: calculated h = %f\n", TempCellData[i].h);

              printf("CMS_FEEDBACK[%d]: id = %d; calc = %d; repeat_calc0 = %d; tot_vol = %f; Ncells = %d\n", ThisTask, P[i].ID, calc, repeat_calc0, TempCellData[i].total_vol, TempCellData[i].Ncells);

#ifdef LOCAL_FEEDBACK_PARTICLES
              if(calc >= 1)
                if(TempCellData[i].NgbVolume == 0.0)
                  {
                    printf("CMS_FEEDBACK[%d]: id = %d; NgbDistance = %g; NgbVolume = %g\n", ThisTask, P[i].ID, TempCellData[i].NgbDistance, TempCellData[i].NgbVolume);
                    terminate("zero radius feedback sphere\n");
                  }

#endif
              if(calc == 1)
                {
#ifdef LOCAL_KINETIC
                  TempCellData[i].fkin = All.LocalFeedbackKineticInjectionFraction;
#else
                  TempCellData[i].fkin = 0.0;
#endif
                  TempCellData[i].SNEnergy = SNEnergy;
                  TempCellData[i].SNMassReturn = All.LocalFeedbackSNMassReturn;
                }

              if(calc == 2)
                {
                  double a = TempCellData[i].aterm1 - TempCellData[i].tot_ekin_b - TempCellData[i].fkin * TempCellData[i].SNEnergy;
                  double b = TempCellData[i].b;
                  double c = TempCellData[i].c;
                  TempCellData[i].deltap = (-b + sqrt(b * b - 4.0 * a * c)) / 2.0 / c;

                  if(a == 0.0)
                    TempCellData[i].deltap = 0.0;

                  count_feedback++;
                  printf("CMS_FEEDBACK: a,b,c,deltap: %f, %f, %f, %f\n", a, b, c, TempCellData[i].deltap);

                }
            }
        }                       //end loop over particles


      /* this drives the processing of local and imported particles, and all the communication */
      generic_comm_pattern(TimeBinsGravity.NActiveParticles, kernel_local, kernel_imported);

    }


  /* loop over all active cells and particles */
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(TempCellData[i].flag_feedback)
        {
          TempCellData[i].flag_feedback = 0;

          double resid = sqrt(TempCellData[i].resid[0] * TempCellData[i].resid[0] + TempCellData[i].resid[1] * TempCellData[i].resid[1] + TempCellData[i].resid[2] * TempCellData[i].resid[2]);
          printf("CMS_FEEDBACK:    Cell %d: resid/deltap %f \n", P[i].ID, resid / TempCellData[i].deltap);
        }
    }



  /* Print statement */
  //printf("CMS_FEEDBACK: Injected %d events\n",count);
  int countall_feedback;
  MPI_Allreduce(&count_feedback, &countall_feedback, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(countall_feedback > 0)
    {
      mpi_printf("CMS_FEEDBACK: Injected %d events\n", countall_feedback);
      //ngb_recompute_nodes();
    }

  myfree(TempCellData);

}


/*! This function represents the core of the star density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
static int find_feedback_cells_evaluate(int target, int mode, int thread_id)
{
  int j, n, k;
  double h2;
  double dx, dy, dz, r2;
  double pi = 3.14159265358979323846;
  /* mode is */
#ifdef LOCAL_FEEDBACK_PARTICLES
  double density_threshold = 10.0;
  density_threshold *= PROTONMASS * All.UnitLength_in_cm * All.UnitLength_in_cm * All.UnitLength_in_cm / All.UnitMass_in_g;
#endif

#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif

  MyDouble *pos, *frame_vel, target_vol, vol;
  double local_resid[3];
  double cell_radius;

  double total_vol, fkin, SNEnergy, SNMassReturn, deltap;
  double aterm1, a, b, c, tot_ekin_b;
  /* determine input parameters for particle */

  int numnodes, *firstnode;
  data_in local, *in;
  data_out out;




  if(mode == MODE_LOCAL_PARTICLES)      /* local particle */
    {
      particle2in(&local, target, 0);
      in = &local;

      numnodes = 1;
      firstnode = NULL;
    }
  else                          /* imported particle */
    {
      in = &DataGet[target];
      generic_get_numnodes(target, &numnodes, &firstnode);
    }

  pos = in->Pos;
  frame_vel = in->Vel;
  double h = in->h;
  if(calc > 0)
    {
      total_vol = in->total_vol;
      fkin = in->fkin;
      SNEnergy = in->SNEnergy;
      SNMassReturn = in->SNMassReturn;
    }

  if(calc == 2)
    deltap = in->deltap;

#ifndef LOCAL_FEEDBACK_PARTICLES
  target_vol = in->Volume;
#endif



  printf("CMS_FEEDBACK[%d]: calc = %d; h = %f; pos = %f %f %f\n", ThisTask, calc, h, pos[0], pos[1], pos[2]);

  double cell_aterm1, cell_b, cell_c, cell_tot_ekin_b;
  double local_vol_sum, local_aterm1, local_b, local_c, local_tot_ekin_b;
  int local_Ncells, ingb;
  double ngb_distance2;


  ngb_distance2 = boxSize_X * boxSize_Z;
  ingb = -1;
  if(calc == 0)
    {
      local_vol_sum = 0.0;
      local_Ncells = 0;

    }

  if(calc == 1)
    {
      local_aterm1 = 0.0;
      local_b = 0.0;
      local_c = 0.0;
      local_tot_ekin_b = 0.0;
    }

  if(calc == 2)
    {
      for(k = 0; k < 3; k++)
        local_resid[k] = 0.0;
    }

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, thread_id, numnodes, firstnode);
  printf("CMS_FEEDBACK[%d]: nfound: %d\n", ThisTask, nfound);
  if(mode == MODE_LOCAL_PARTICLES)
    while(nfound == 0)
      {
        printf("CMS_FEEDBACK[%d]: Increasing h: %g\n", ThisTask, h * 2.0);
        h *= 2.0;
        nfound = ngb_treefind_variable_threads(pos, h, target, mode, thread_id, numnodes, firstnode);
      }

  h2 = h * h;
  for(n = 0; n < nfound; n++)
    {
      j = Thread[thread_id].Ngblist[n];

      if(P[j].Mass > 0 && P[j].ID != 0 && P[j].Type == 0)
        {
          dx = pos[0] - SphP[j].Center[0];
          dy = pos[1] - SphP[j].Center[1];
          dz = pos[2] - SphP[j].Center[2];

#ifdef PERIODIC
#if !defined(REFLECTIVE_X)
          if(dx < -boxHalf_X)
            dx += boxSize_X;
          if(dx > boxHalf_X)
            dx -= boxSize_X;
#endif
#if !defined(REFLECTIVE_Y)
          if(dy < -boxHalf_Y)
            dy += boxSize_Y;
          if(dy > boxHalf_Y)
            dy -= boxSize_Y;
#endif
#if !defined(REFLECTIVE_Z)
          if(dz < -boxHalf_Z)
            dz += boxSize_Z;
          if(dz > boxHalf_Z)
            dz -= boxSize_Z;
#endif
#endif

          r2 = dx * dx + dy * dy + dz * dz;

          if(r2 < ngb_distance2 && calc == -1)
            {
              printf("CMS_FEEDBACK[%d]: reassigning nearest particle %e %e\n", ThisTask, r2, ngb_distance2);
              ngb_distance2 = r2;
              ingb = j;
            }

          if(r2 < h2 && calc > -1)
            {
              //printf("CMS_FEEDBACK[%d]: r2 h2 = %f %f\n",ThisTask,r2,h2);
              /* compute total volume for normalization purposes */

              if(calc == 0)
                {
                  local_vol_sum += SphP[j].Volume;
                  local_Ncells += 1;
                  printf("CMS_FEEDBACK[%d]: Calculation0; Cell ID: %d; Cell Vol: %f; Tot Cell Vol: %f\n", ThisTask, P[j].ID, SphP[j].Volume, local_vol_sum);
                }
              /* compute kinetic energy injection constants */
              if(calc == 1)
                {
                  double momentum[3], mass;

                  mass = P[j].Mass;
                  for(k = 0; k < 3; k++)
                    momentum[k] = SphP[j].Momentum[k] - mass * frame_vel[k];

                  vol = SphP[j].Volume;

                  double fvol = vol / total_vol;
                  double fvol_cell_centered = vol / (total_vol - target_vol);
                  double pvec[3];
                  printf("CMS_FEEDBACK[%d]: Calculation1; Cell ID: %d; total_vol: %f; fvol: %f\n", ThisTask, P[j].ID, total_vol, fvol);
                  if(r2 > 0.0)
                    {
                      double r = sqrt(r2);
                      pvec[0] = dx / r;
                      pvec[1] = dy / r;
                      pvec[2] = dz / r;

                      double pb2 = (momentum[0] * momentum[0] + momentum[1] * momentum[1] + momentum[2] * momentum[2]);
                      double dM = SNMassReturn;

                      cell_tot_ekin_b = 0.5 * pb2 / mass;
                      local_tot_ekin_b += cell_tot_ekin_b;

                      cell_aterm1 = 0.5 * pb2 / (mass + dM * fvol);
                      local_aterm1 += cell_aterm1;

                      cell_b = (momentum[0] * pvec[0] + momentum[1] * pvec[1] + momentum[2] * pvec[2]) * fvol_cell_centered / (mass + dM * fvol);
                      local_b += cell_b;
                      cell_c = 0.5 * fvol_cell_centered * fvol_cell_centered / (mass + dM * fvol);
                      local_c += cell_c;
                      //printf("CMS_FEEDBACK[%d]: Local %d: %e, %e, %e, %e\n",ThisTask,j,local_aterm1,local_b,local_c,local_tot_ekin_b);
                      //printf("CMS_FEEDBACK[%d]: Cell %d: %e, %e, %e, %e\n",ThisTask,P[j].ID,cell_aterm1,cell_b,cell_c,cell_tot_ekin_b);
                      //printf("CMS_FEEDBACK[%d]: Cell %d: %e %e %e; %e %e %e\n",ThisTask,P[j].ID,pvec[0],pvec[1],pvec[2],momentum[0],momentum[1],momentum[2]);
                      //printf("CMS_FEEDBACK[%d]: Cell %d: %e, %e %e\n",ThisTask,P[j].ID,frame_vel[0],frame_vel[1],frame_vel[2]);
                    }
                }

              /* add feedback to cells */
              if(calc == 2)
                {
                  double momentum[3], mass;
                  vol = SphP[j].Volume;

                  double pvec[3];
                  if(r2 > 0.0)
                    {
                      double r = sqrt(r2);
                      pvec[0] = dx / r;
                      pvec[1] = dy / r;
                      pvec[2] = dz / r;
                    }
                  else
                    {
                      pvec[0] = 0.0;
                      pvec[1] = 0.0;
                      pvec[2] = 0.0;
                    }

                  double fvol = vol / total_vol;
                  double fvol_cell_centered = vol / (total_vol - target_vol);
                  mass = P[j].Mass;
                  for(k = 0; k < 3; k++)
                    momentum[k] = SphP[j].Momentum[k] - mass * frame_vel[k];

                  double pb2 = (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] + SphP[j].Momentum[2] * SphP[j].Momentum[2]);

                  double kin_b = 0.5 * pb2 / P[j].Mass;

                  P[j].Mass += SNMassReturn * fvol;
                  for(k = 0; k < 3; k++)
                    momentum[k] += deltap * pvec[k] * fvol_cell_centered;

                  SphP[j].Energy += SNEnergy * (1.0 - fkin) * fvol;

                  mass = P[j].Mass;
                  for(k = 0; k < 3; k++)
                    SphP[j].Momentum[k] = momentum[k] + mass * frame_vel[k];

                  //printf("CMS_FEEDBACK[%d]: Cell %d: %f, %f, %f %f %f, %f, %f\n",ThisTask,P[j].ID,fvol,fvol_cell_centered,pvec[0],pvec[1],pvec[2],deltap,total_vol);
                  struct pv_update_data pvd;
                  pvd.atime = pvd.hubble_a = pvd.a3inv = 1.0;
                  update_primitive_variables_single(P, SphP, j, &pvd);

                  pb2 = (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] + SphP[j].Momentum[2] * SphP[j].Momentum[2]);

                  double kin_a = 0.5 * pb2 / P[j].Mass;

                  double dkin = kin_a - kin_b;

                  SphP[j].Energy += dkin;

                  update_internal_energy(P, SphP, j, &pvd);
                  set_pressure_of_cell(j);
                  SphP[j].Csnd = get_sound_speed(j);

                  for(k = 0; k < 3; k++)
                    local_resid[k] += deltap * pvec[k] * fvol;
                  //printf("CMS_FEEDBACK[%d]: Cell %d: %f, %f, %f %f %f\n",ThisTask,P[j].ID, SNEnergy*(1.0-fkin)*fvol,SNMassReturn*fvol,dkin,deltap,total_vol);
                }
            }                   //end r2 < h2
        }
    }

  /* now store the results at the appropriate place */

  if(calc == -1)
    {
#ifdef LOCAL_FEEDBACK_PARTICLES
      out.NgbDistance = sqrt(ngb_distance2);
      if(ingb >= 0)
        out.NgbVolume = SphP[ingb].Volume;
      else
        out.NgbVolume = 0.0;

      printf("CMS_FEEDBACK[%d]: id = %d; nearest cell has volume %e at distance %f\n", ThisTask, P[ingb].ID, SphP[ingb].Volume, sqrt(ngb_distance2));
#endif
    }
  if(calc == 0)
    {
      out.total_vol = local_vol_sum;
      out.Ncells = local_Ncells;

    }

  if(calc == 1)
    {
      out.aterm1 = local_aterm1;
      out.b = local_b;
      out.c = local_c;
      out.tot_ekin_b = local_tot_ekin_b;
    }

  if(calc == 2)
    {
      out.resid[0] = local_resid[0];
      out.resid[1] = local_resid[1];
      out.resid[2] = local_resid[2];
      out.Ncells = nfound;
    }

  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#endif
