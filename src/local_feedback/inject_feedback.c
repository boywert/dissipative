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

#if defined(LOCAL_FEEDBACK)

#include "../parallel_logs.h"

static int calc;                /* determines phase of calculation */
static int find_feedback_cells_evaluate(int target, int mode, int thread_id);

static struct plog_data plog;

/* communication structures */
typedef struct
{
  MyDouble Pos[3];
  MyDouble Vel[3];
  MyDouble Volume;

  MyDouble SNEnergy;
  MyDouble fkin;
  MyDouble fcr;
  MyDouble SNMassReturn;
  MyDouble SNMetalReturn;
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
  MyDouble fcr;
  MyDouble SNMassReturn;
  MyDouble SNMetalReturn;
  MyDouble deltap;
  MyDouble h, h_h, h_l, dh;

  MyDouble aterm1, b, c, tot_ekin_b;
  MyDouble resid[3];
  int Ncells;

  /* If feedback is centered on particles, record info about closest hydro cell */

  MyDouble NgbDistance;
  MyDouble NgbVolume;
  MyDouble NgbVel[3];


}
 *TempCellData;

/* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in * in, int i, int firstnode)
{
  int k;

#if defined(LOCAL_FEEDBACK_PARTICLES) || defined(EXTERNALSHEARBOX_KSRATE_RANDOM) || defined(EXTERNALSHEARBOX_MIXED_INJECTION)
  for(k = 0; k < 3; k++)
    in->Pos[k] = P[i].Pos[k];

  if(TempCellData[i].Ncells == 0)
    in->Volume = 0.0;
  else
    in->Volume = TempCellData[i].NgbVolume;

#else
  for(k = 0; k < 3; k++)
    in->Pos[k] = SphP[i].Center[k];

  in->Volume = SphP[i].Volume;
#endif

#if defined(EXTERNALSHEARBOX_KSRATE_RANDOM) || defined(EXTERNALSHEARBOX_MIXED_INJECTION)
  for(k = 0; k < 3; k++)
    in->Vel[k] = TempCellData[i].NgbVel[k];
#else
  for(k = 0; k < 3; k++)
    in->Vel[k] = P[i].Vel[k];
#endif

  in->total_vol = TempCellData[i].total_vol;
  in->h = TempCellData[i].h;
  if(calc > 0)
    {
      in->SNEnergy = TempCellData[i].SNEnergy;
      in->fkin = TempCellData[i].fkin;
      in->fcr = TempCellData[i].fcr;
      in->SNMassReturn = TempCellData[i].SNMassReturn;
      in->deltap = TempCellData[i].deltap;
    }
  in->Firstnode = firstnode;
}


typedef struct
{
  MyDouble total_vol;
  MyDouble aterm1, b, c, tot_ekin_b, h;
  MyDouble resid[3];
  int Ncells;

  MyDouble NgbDistance;
  MyDouble NgbVolume;
  MyDouble NgbRadius;
  MyDouble NgbVel[3];

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
#if defined(LOCAL_FEEDBACK_PARTICLES) || defined(EXTERNALSHEARBOX_KSRATE_RANDOM) || defined(EXTERNALSHEARBOX_MIXED_INJECTION)
          TempCellData[i].NgbDistance = out->NgbDistance;
          TempCellData[i].NgbVolume = out->NgbVolume;
          for(k = 0; k < 3; k++)
            TempCellData[i].NgbVel[k] = out->NgbVel[k];
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
#if defined(LOCAL_FEEDBACK_PARTICLES) || defined(EXTERNALSHEARBOX_KSRATE_RANDOM) || defined(EXTERNALSHEARBOX_MIXED_INJECTION)
          if(out->NgbDistance < TempCellData[i].NgbDistance)
            {
              TempCellData[i].NgbDistance = out->NgbDistance;
              TempCellData[i].NgbVolume = out->NgbVolume;
              for(k = 0; k < 3; k++)
                TempCellData[i].NgbVel[k] = out->NgbVel[k];
            }
#endif
        }
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

/* Routine that injects feedback either around previously created star particles, or if star particles
   are not used, instantaneously injected based on local gas properties.
   If instantaneous injection is used, the injection rate is computed either based on the KS rate or
   based on the local dynamical time.
   The total amount of energy injected per SN event, and the fraction of that energy that takes a kinetic 
   form, is set by run time parameters. */

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
          continue;

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
        find_feedback_cells_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }

}

void inject_feedback(void)
{
  CPU_Step[CPU_MISC] += measure_time();
  mpi_printf("CMS_FEEDBACK: Starting Starformation and Feedback Routine\n");

  int i;
  double rho, sfr, dz, dt, prob, sech, deltaM;
  double rand;
  double SNEnergy = All.LocalFeedbackSNEnergy * 1e51 / All.UnitVelocity_in_cm_per_s / All.UnitVelocity_in_cm_per_s / All.UnitMass_in_g;
  double SNRate = All.LocalFeedbackSNRate;

#ifdef LOCAL_FEEDBACK_PARTICLES
  double StarParticleCreationTime;
  double MinimumStarMass = 0.5 * All.TargetGasMass;
  double SNTimeDelay_max = All.LocalFeedbackSNTimeDelay;
  double SNTimeDelay = All.LocalFeedbackSNTimeDelay;
  double SNTimeSpread_max = All.LocalFeedbackSNTimeSpread;
  double SNTimeSpread = All.LocalFeedbackSNTimeSpread;
  double tform, Age, N_SN;
#endif

  TempCellData = mymalloc("TemCellData", NumPart * sizeof(struct temp_celldata));

  parallel_log_allocate( &plog, 0 );

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

          rand = get_random_number();

          if(rand < prob)
            {
              TempCellData[i].flag_feedback = 1;
              P[i].FeedbackDone += 1;
	    }
        }

#else /* inject feedback without real star particles */

#if defined(EXTERNALSHEARBOX_KSRATE_RANDOM)  || defined(EXTERNALSHEARBOX_MIXED_INJECTION) /* inject feedback around dummy particles */
      if(P[i].Type == 4)
        {
          TempCellData[i].flag_feedback = 1;
          pl_fprintf( &plog, "CMS_FEEDBACK[%d]: center t,x,y,z,ty: %f %f %f %f %d\n", 
		      ThisTask, All.Time, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2],P[i].Type);
	  printf( "CMS_FEEDBACK[%d]: center t,x,y,z,ty: %f %f %f %f %d\n", 
		  ThisTask, All.Time, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2],P[i].Type);

        }
#endif

/* inject feedback centered on cells (no particles!) */
#ifndef EXTERNALSHEARBOX_KSRATE_RANDOM /* inject feedback centered on cells (no particles!) */

      /* just consider gas cells */
      if(P[i].Type == 0){
	
	/* skip cells that have been swallowed or eliminated */
	if(P[i].Mass == 0 && P[i].ID == 0)
	  continue;
		
	dt = (All.HighestActiveTimeBin ? (((int) 1) << All.HighestActiveTimeBin) : 0) * All.Timebase_interval;
	deltaM = SphP[i].Sfr*dt;
	prob = deltaM*SNRate/100.0; 
	rand = get_random_number();
	
	/* select a cell randomly with a certain rate */
	if(rand < prob){
	  TempCellData[i].flag_feedback = 1;
	  pl_fprintf( &plog, "CMS_FEEDBACK[%d]: center t,x,y,z,ty: %f %f %f %f %d\n",
		      ThisTask, All.Time, SphP[i].Center[0],SphP[i].Center[1],SphP[i].Center[2],P[i].Type);
	  printf( "CMS_FEEDBACK[%d]: center t,x,y,z: %f %f %f %f %d\n",
		  ThisTask, All.Time, SphP[i].Center[0],SphP[i].Center[1],SphP[i].Center[2],P[i].Type);
	}
      }
#endif
#endif
    }
  mpi_printf("CMS_FEEDBACK: Finished selecting explosions\n");
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

  for(calc = 0; calc < 3; calc++)
    {
      int idx;
      int enlarge_radius = 1;
      /* loop over all active cells and particles */
      for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
        {
          i = TimeBinsGravity.ActiveParticleList[idx];
          if(i < 0)
            continue;

          if(TempCellData[i].flag_feedback)
            {
	      /* on first pass, tie radius to softening length */
              if(calc == 0)
                {
#if defined(LOCAL_FEEDBACK_PARTICLES) || defined(EXTERNALSHEARBOX_KSRATE_RANDOM) || defined(EXTERNALSHEARBOX_MIXED_INJECTION)
		  if(P[i].Type == 4)
		    cell_radius = 25.0 * All.MinimumComovingHydroSoftening;
		  else
		    cell_radius = pow(SphP[i].Volume * 3.0 / 4.0 / M_PI, 1.0 / 3.0);
#else
                  cell_radius = pow(SphP[i].Volume * 3.0 / 4.0 / M_PI, 1.0 / 3.0);
#endif
                  TempCellData[i].h = 3.0 * cell_radius;
                  TempCellData[i].h_l = 0.0;
                  TempCellData[i].h_h = 0.0;
                  TempCellData[i].dh = 0.5 * cell_radius;

                }
              h2 = TempCellData[i].h * TempCellData[i].h;
              
#ifdef LOCAL_FEEDBACK_PARTICLES
              if(calc >= 1)
                if(TempCellData[i].h == 0.0)
		  terminate("CMS_FEEDBACK[%d]: zero radius feedback sphere: id = %d; NgbDistance = %g; NgbVolume = %g\n", 
			    ThisTask, P[i].ID, TempCellData[i].NgbDistance, TempCellData[i].NgbVolume);                  
#endif
              if(calc == 1)
                {
#ifdef LOCAL_KINETIC
                  TempCellData[i].fkin = All.LocalFeedbackKineticInjectionFraction;
#else
                  TempCellData[i].fkin = 0.0;
#endif
#if defined(COSMIC_RAYS) && !defined(COSMIC_RAYS_SHOCK_ACCELERATION)
                  TempCellData[i].fcr = All.LocalFeedbackCRInjectionFraction;
#endif

	      TempCellData[i].SNEnergy = SNEnergy;
	      TempCellData[i].SNMassReturn = All.LocalFeedbackSNMassReturn;
#ifdef METALS
	      TempCellData[i].SNMetalReturn = 1.0;
#else
	      TempCellData[i].SNMetalReturn = 0.0;
#endif

	    }
	    
	    if(calc == 2){
	      double a = TempCellData[i].aterm1 - TempCellData[i].tot_ekin_b - TempCellData[i].fkin*TempCellData[i].SNEnergy;
	      double b = TempCellData[i].b;
	      double c = TempCellData[i].c;
	      TempCellData[i].deltap = (-b + sqrt(b*b - 4.0*a*c))/2.0/c;
	      
	      if(a == 0.0)
		TempCellData[i].deltap = 0.0;
	      
	      count_feedback++;
	    }
	  }
      } //end loop over particles

      while(enlarge_radius){	

	generic_comm_pattern(TimeBinsGravity.NActiveParticles, kernel_local, kernel_imported);

	enlarge_radius = 0;
	if(calc == 0)
	  {
	    for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
	      {
		i = TimeBinsGravity.ActiveParticleList[idx];
		if(i < 0)
		  continue;

		if(TempCellData[i].flag_feedback && (TempCellData[i].Ncells < 32 || TempCellData[i].Ncells > 32))
		  {
		    enlarge_radius = 1;
		    
		    if(TempCellData[i].Ncells < 32)
		      {
			TempCellData[i].h_l = TempCellData[i].h;
			TempCellData[i].h += TempCellData[i].dh;
			
		      }
		    if(TempCellData[i].Ncells > 32)
		      {
			TempCellData[i].h_h = TempCellData[i].h;
			TempCellData[i].h -= TempCellData[i].dh;
		      }
		    
		    if(TempCellData[i].h_h > 0.0 && TempCellData[i].h_l > 0.0)
                        TempCellData[i].dh = 0.5 * (TempCellData[i].h_h - TempCellData[i].h_l);
		   
		  }
	      }
	  }
	MPI_Allreduce(MPI_IN_PLACE, &enlarge_radius, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      }

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
          //pl_fprintf( &plog, "CMS_FEEDBACK:    Cell %d: resid/deltap %f \n", P[i].ID, resid / TempCellData[i].deltap);
#ifdef EXTERNALSHEARBOX_KSRATE_RANDOM
          P[i].ID = 0;
          P[i].Mass = 0;
          timebin_remove_particle(&TimeBinsGravity, idx, P[i].TimeBinGrav);
#endif
        }
    }

/* Print statement */
  int countall_feedback;
  MPI_Allreduce(&count_feedback, &countall_feedback, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(countall_feedback > 0)
    mpi_printf("CMS_FEEDBACK: Injected %d events\n", countall_feedback);
  
  parallel_log_clear( &plog, PLOG_LOCAL_FEEDBACK );
  parallel_log_free( &plog );

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
 
  double total_vol,fkin,fcr,SNEnergy,SNMassReturn,deltap,SNMetalReturn;
  double aterm1,a,b,c,tot_ekin_b;
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
      fcr = in->fcr;
      SNEnergy = in->SNEnergy;
      SNMassReturn = in->SNMassReturn;
      SNMetalReturn = in->SNMetalReturn;
    }

  if(calc == 2)
    deltap = in->deltap;

#if !defined(LOCAL_FEEDBACK_PARTICLES) && !defined(EXTERNALSHEARBOX_KSRATE_RANDOM) && !defined(EXTERNALSHEARBOX_MIXED_INJECTION)
  target_vol = in->Volume;
#endif

  double cell_aterm1, cell_b, cell_c, cell_tot_ekin_b;
  double local_vol_sum, local_aterm1, local_b, local_c, local_tot_ekin_b;
  int local_Ncells, ingb;
  double ngb_distance2;
#if defined(LOCAL_FEEDBACK_PARTICLES) || defined(EXTERNALSHEARBOX_KSRATE_RANDOM) || defined(EXTERNALSHEARBOX_MIXED_INJECTION)
  ngb_distance2 = boxSize_X * boxSize_Z;
  ingb = -1;
#endif

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
#if defined(LOCAL_FEEDBACK_PARTICLES) || defined(EXTERNALSHEARBOX_KSRATE_RANDOM) || defined(EXTERNALSHEARBOX_MIXED_INJECTION)
          if(r2 < ngb_distance2 && calc == 0)
            {
	      ngb_distance2 = r2;
              ingb = j;
            }
#endif
	  if(r2 < h2)
	    {
	      /* compute total volume for normalization purposes */
	      
	      if(calc == 0){
		local_vol_sum += SphP[j].Volume;
		local_Ncells += 1;
	      }
	      /* compute kinetic energy injection constants */
	      if(calc == 1){
		double momentum[3], mass;
		
		mass = P[j].Mass;
		for(k = 0; k < 3; k++)
		  momentum[k] = SphP[j].Momentum[k] - mass*frame_vel[k];
		
		vol = SphP[j].Volume;
	    
		double fvol = vol / total_vol;
		double fvol_cell_centered = vol / (total_vol - target_vol);
		double pvec[3];
		
		if(r2 > 0.0){
		  double r = sqrt(r2);
		  pvec[0] = dx/r;
		  pvec[1] = dy/r;
		  pvec[2] = dz/r;
			        
		  double pb2 = (momentum[0] * momentum[0] + momentum[1] * momentum[1] + momentum[2] * momentum[2]);
		  double dM = SNMassReturn;
	      
		  cell_tot_ekin_b = 0.5 * pb2 / mass;
		  local_tot_ekin_b += cell_tot_ekin_b;
		  
		  cell_aterm1 = 0.5 * pb2 / (mass + dM * fvol);
		  local_aterm1 += cell_aterm1;
		  
		  cell_b = (momentum[0] * pvec[0] + momentum[1] * pvec[1] + momentum[2] * pvec[2])*fvol_cell_centered / (mass + dM * fvol);
		  local_b += cell_b;
		  cell_c = 0.5 * fvol_cell_centered * fvol_cell_centered / (mass + dM * fvol);
		  local_c += cell_c;
		}
	      }
	  
	      /* add feedback to cells */
	      if(calc == 2){		  
		double momentum[3], mass;
		vol = SphP[j].Volume;
		
		double pvec[3];
		if(r2 > 0.0){
		  double r = sqrt(r2);
		  pvec[0] = dx/r;
		  pvec[1] = dy/r;
		  pvec[2] = dz/r;
		}else{
		  pvec[0] = 0.0;
		  pvec[1] = 0.0;
		  pvec[2] = 0.0;
		}
		  
		double fvol = vol/total_vol;
		double fvol_cell_centered = vol / (total_vol - target_vol);
		mass = P[j].Mass;
		for(k = 0; k < 3; k++)
		  momentum[k] = SphP[j].Momentum[k] - mass*frame_vel[k];
		
		double pb2 = (SphP[j].Momentum[0] * SphP[j].Momentum[0] + 
			      SphP[j].Momentum[1] * SphP[j].Momentum[1] + 
			      SphP[j].Momentum[2] * SphP[j].Momentum[2]);
	    
		double kin_b = 0.5*pb2/P[j].Mass;
		
#ifdef METALS
		SNMetalReturn = 1.0;
		double metal_mass = SphP[j].Metallicity*P[j].Mass;
		metal_mass += SNMetalReturn*fvol;
		SphP[j].MassMetallicity = metal_mass;
		SphP[j].Metallicity = metal_mass/P[j].Mass;
		assert(SphP[j].Metallicity >= 0);
		P[j].Metallicity = SphP[j].Metallicity;
		assert(P[j].Metallicity >= 0);
#endif
		
		P[j].Mass += SNMassReturn*fvol;
		for(k = 0; k < 3; k++)
		  momentum[k] += deltap*pvec[k]*fvol_cell_centered;
		
		double fraction_of_thermal_energy;
		fraction_of_thermal_energy = 1.0 - fkin;
#if defined(COSMIC_RAYS) && !defined(COSMIC_RAYS_SHOCK_ACCELERATION)
		fraction_of_thermal_energy -= fcr;
		SphP[j].CR_Energy += SNEnergy * fcr * fvol;
#endif
		SphP[j].Energy += SNEnergy*fraction_of_thermal_energy*fvol;
		
		mass = P[j].Mass;
		for(k = 0; k < 3; k++)
		  SphP[j].Momentum[k] = momentum[k]+mass*frame_vel[k];
		
		struct pv_update_data pvd;
		pvd.atime = pvd.hubble_a = pvd.a3inv = 1.0;
		update_primitive_variables_single(P, SphP, j, &pvd);
		
		pb2 = (SphP[j].Momentum[0] * SphP[j].Momentum[0] + 
		       SphP[j].Momentum[1] * SphP[j].Momentum[1] + 
		       SphP[j].Momentum[2] * SphP[j].Momentum[2]);
		
		double kin_a = 0.5*pb2/P[j].Mass;
		double dkin = kin_a - kin_b;
		
		SphP[j].Energy += dkin;
		
		update_internal_energy(P, SphP, j, &pvd);
		set_pressure_of_cell(j);
		SphP[j].Csnd = get_sound_speed(j);
	    
		for(k = 0; k < 3; k++)
		  local_resid[k] += deltap*pvec[k]*fvol;
	      }      
	    } //end r2 < h2
	}
    }
  
  /* now store the results at the appropriate place */
  if(calc == 0)
    {
#if defined(LOCAL_FEEDBACK_PARTICLES) || defined(EXTERNALSHEARBOX_KSRATE_RANDOM) || defined(EXTERNALSHEARBOX_MIXED_INJECTION)
      out.NgbDistance = sqrt(ngb_distance2);
      if(ingb >= 0)
	{
	  out.NgbVolume = SphP[ingb].Volume;
	  for(k = 0; k < 3; k++)
	    out.NgbVel[k] = P[ingb].Vel[k];
	}else
	{
	  out.NgbVolume = 0.0;
	  for(k = 0; k < 3; k++)
	    out.NgbVel[k] = 0.0;
	}
#else
      out.NgbDistance = 0.0;
      out.NgbVolume = target_vol;
      for(k = 0; k < 3; k++)
	out.NgbVel[k] = frame_vel[k];
#endif
      out.Ncells = local_Ncells;
      out.total_vol = local_vol_sum;
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
