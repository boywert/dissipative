/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/voronoi_exchange.c
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

#include "allvars.h"
#include "proto.h"
#include "voronoi.h"

#ifdef VORONOI

struct data_primexch_compare
{
  int rank, task, index;
} *SortPrimExch, *SortPrimExch2;



void mesh_setup_exchange(void)
{
  if(All.TotNumGas == 0)
    return;

  TIMER_START(CPU_MESH_EXCHANGE);

  int listp;
  struct indexexch
  {
    int task, index;
  } *tmpIndexExch, *IndexExch;
  int i, j, p, task, off, count;
  int ngrp, recvTask, place;


  for(j = 0; j < NTask; j++)
    Mesh_Send_count[j] = 0;

  for(i = 0; i < NumGasInMesh; i++)
    {
      p = List_InMesh[i];

      listp = List_P[p].firstexport;
      while(listp >= 0)
        {
          if(ListExports[listp].origin != ThisTask)
            {
              Mesh_Send_count[ListExports[listp].origin]++;
            }
          listp = ListExports[listp].nextexport;
        }
    }

  MPI_Alltoall(Mesh_Send_count, 1, MPI_INT, Mesh_Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, Mesh_nimport = 0, Mesh_nexport = 0, Mesh_Recv_offset[0] = 0, Mesh_Send_offset[0] = 0; j < NTask; j++)
    {
      Mesh_nimport += Mesh_Recv_count[j];
      Mesh_nexport += Mesh_Send_count[j];

      if(j > 0)
        {
          Mesh_Send_offset[j] = Mesh_Send_offset[j - 1] + Mesh_Send_count[j - 1];
          Mesh_Recv_offset[j] = Mesh_Recv_offset[j - 1] + Mesh_Recv_count[j - 1];
        }
    }

  IndexExch = (struct indexexch *) mymalloc("IndexExch", Mesh_nimport * sizeof(struct indexexch));
  tmpIndexExch = (struct indexexch *) mymalloc("tmpIndexExch", Mesh_nexport * sizeof(struct indexexch));

  /* prepare data for export */
  for(j = 0; j < NTask; j++)
    Mesh_Send_count[j] = 0;

  for(i = 0; i < NumGasInMesh; i++)
    {
      p = List_InMesh[i];

      listp = List_P[p].firstexport;
      while(listp >= 0)
        {
          if((task = ListExports[listp].origin) != ThisTask)
            {
              place = ListExports[listp].index;
              off = Mesh_Send_offset[task] + Mesh_Send_count[task]++;

              tmpIndexExch[off].task = ThisTask;
              tmpIndexExch[off].index = place;
            }
          listp = ListExports[listp].nextexport;
        }
    }

  /* exchange data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Mesh_Send_count[recvTask] > 0 || Mesh_Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&tmpIndexExch[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask]
                           * sizeof(struct indexexch), MPI_BYTE, recvTask, TAG_DENS_A,
                           &IndexExch[Mesh_Recv_offset[recvTask]], Mesh_Recv_count[recvTask] * sizeof(struct indexexch), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  myfree(tmpIndexExch);

  /* now we need to associate the imported data with the points stored in the DP[] array */

  SortPrimExch = (struct data_primexch_compare *) mymalloc("SortPrimExch", Mesh_nimport * sizeof(struct data_primexch_compare));

  for(i = 0; i < Mesh_nimport; i++)
    {
      SortPrimExch[i].rank = i;
      SortPrimExch[i].task = IndexExch[i].task;
      SortPrimExch[i].index = IndexExch[i].index;
    }

  /* let sort the data according to task and index */
  mysort(SortPrimExch, Mesh_nimport, sizeof(struct data_primexch_compare), compare_primexch);

  SortPrimExch2 = (struct data_primexch_compare *) mymalloc("SortPrimExch2", Mesh.Ndp * sizeof(struct data_primexch_compare));

  for(i = 0, count = 0; i < Mesh.Ndp; i++)
    {
      if(Mesh.DP[i].task != ThisTask)
        {
          SortPrimExch2[count].rank = i;
          SortPrimExch2[count].task = Mesh.DP[i].task;
          SortPrimExch2[count].index = Mesh.DP[i].index;
          count++;
        }
    }

  /* let sort according to task and index */
  mysort(SortPrimExch2, count, sizeof(struct data_primexch_compare), compare_primexch);

  /* count can be larger than nimport because a foreigh particle can appear
     multiple times on the local domain, due to periodicity */

  for(i = 0, j = 0; i < count; i++)
    {
      if(SortPrimExch2[i].task != SortPrimExch[j].task || SortPrimExch2[i].index != SortPrimExch[j].index)
        j++;

      if(j >= Mesh_nimport)
        terminate("j >= Mesh_nimport");

      Mesh.DP[SortPrimExch2[i].rank].index = SortPrimExch[j].rank;      /* note: this change is now permanent and available for next exchange */
    }

  myfree(SortPrimExch2);
  myfree(SortPrimExch);
  myfree(IndexExch);



  /* allocate structures needed to exchange the actual information for ghost cells */
  PrimExch = (struct primexch *) mymalloc_movable(&PrimExch, "PrimExch", Mesh_nimport * sizeof(struct primexch));
  GradExch = (struct grad_data *) mymalloc_movable(&GradExch, "GradExch", Mesh_nimport * sizeof(struct grad_data));
#ifdef TVD_SLOPE_LIMITER
  GradExchUl = (struct grad_data *) mymalloc_movable(&GradExchUl, "GradExchUl", Mesh_nimport * sizeof(struct grad_data));
#endif

#if defined(RT_ADVECT) || defined(MRT)
  RTPrimExch = (struct rt_primexch *) mymalloc_movable(&RTPrimExch, "RTPrimExch", Mesh_nimport * sizeof(struct rt_primexch));
  RTGradExch = (struct rt_grad_data *) mymalloc_movable(&RTGradExch, "RTGradExch", Mesh_nimport * sizeof(struct rt_grad_data));
#endif

#ifdef SECOND_DERIVATIVES
  HessianExch = (struct hessian_data *) mymalloc_movable(&HessianExch, "HessianExch", Mesh_nimport * sizeof(struct hessian_data));
#endif


  TIMER_STOP(CPU_MESH_EXCHANGE);
}

void exchange_primitive_variables(void)
{
  if(All.TotNumGas == 0)
    return;

  TIMER_START(CPU_MESH_EXCHANGE);

  int listp;
  struct primexch *tmpPrimExch;
  int i, j, p, task, off;
  int ngrp, recvTask, place;
#if defined(DVR_RENDER) || defined(TGCHEM) || defined(MHD_CT)
  int k;
#endif

  tmpPrimExch = (struct primexch *) mymalloc("tmpPrimExch", Mesh_nexport * sizeof(struct primexch));

  /* prepare data for export */
  for(j = 0; j < NTask; j++)
    Mesh_Send_count[j] = 0;

  for(i = 0; i < NumGasInMesh; i++)
    {
      p = List_InMesh[i];

      listp = List_P[p].firstexport;
      while(listp >= 0)
        {
          if((task = ListExports[listp].origin) != ThisTask)
            {
              place = ListExports[listp].index;
              off = Mesh_Send_offset[task] + Mesh_Send_count[task]++;

              tmpPrimExch[off].Volume = SphP[place].Volume;

              tmpPrimExch[off].Density = SphP[place].Density;
#if defined(DEREFINE_GENTLY) || defined(CONDUCTION_SATURATION) || (defined(COSMIC_RAYS) && defined(SHOCK_FINDER_ON_THE_FLY) || defined(SHOCK_FINDER_BEFORE_OUTPUT)) || defined(SPECIAL_RELATIVITY) || defined(GENERAL_RELATIVITY) || defined(NON_LINEAR_SLOPE_LIMITERS) || defined(CALCULATE_QUANTITIES_IN_POSTPROCESS)
              tmpPrimExch[off].Utherm = SphP[place].Utherm;
#endif

#ifdef DVR_RENDER
              for(k = 0; k < DVR_NUM_FIELDS; k++)
                tmpPrimExch[off].DvrFields[k] = SphP[place].DvrFields[k];
#endif
              tmpPrimExch[off].Pressure = SphP[place].Pressure;
#ifdef COSMIC_RAYS
              tmpPrimExch[off].CR_Pressure = SphP[place].CR_Pressure;
              tmpPrimExch[off].CR_SpecificEnergy = SphP[place].CR_SpecificEnergy;
#ifdef COSMIC_RAYS_STREAMING
              tmpPrimExch[off].CR_Chi = SphP[place].CR_Chi;
#endif
#endif
#ifdef TGCHEM
              for(k = 0; k < TGCHEM_NUM_ABUNDANCES; k++)
                tmpPrimExch[off].Abund[k] = SphP[place].Abund[k];

              tmpPrimExch[off].Gamma = SphP[place].Gamma;
#endif

#ifdef MRT_LSF_GRADIENTS
	      tmpPrimExch[off].HI = SphP[place].HI ;
	      tmpPrimExch[off].HII = SphP[place].HII ;
	      tmpPrimExch[off].Ne = SphP[place].Ne ;
#ifdef MRT_INCLUDE_HE
	      tmpPrimExch[off].HeI = SphP[place].HeI ;
	      tmpPrimExch[off].HeII = SphP[place].HeII ;
	      tmpPrimExch[off].HeIII = SphP[place].HeIII ;
#endif
	  
#endif

#ifdef MHD
              tmpPrimExch[off].B[0] = SphP[place].B[0];
              tmpPrimExch[off].B[1] = SphP[place].B[1];
              tmpPrimExch[off].B[2] = SphP[place].B[2];
#ifdef MHD_POWELL
              tmpPrimExch[off].DivB = SphP[place].DivB;
#endif
#ifdef MHD_CT
              for(k = 0; k < 3; k++)
                {
                  tmpPrimExch[off].A[k] = SphP[place].A[k];
                }
              tmpPrimExch[off].TimeLastBUpdate = SphP[place].TimeLastBUpdate;
#endif
#ifdef NON_IDEAL_MHD
#if defined(OHMIC_DIFFUSION) || defined(IMPLICIT_OHMIC_DIFFUSION) || defined(AMBIPOLAR_DIFFUSION)
              tmpPrimExch[off].CurlB[0] = SphP[place].CurlB[0];
              tmpPrimExch[off].CurlB[1] = SphP[place].CurlB[1];
              tmpPrimExch[off].CurlB[2] = SphP[place].CurlB[2];
#endif
#endif
#endif
#ifdef VARIABLE_GAMMA
              tmpPrimExch[off].GammaC = SphP[place].GammaC;
              tmpPrimExch[off].GammaE = SphP[place].GammaE;
#endif
#ifdef USE_ENTROPY_FOR_COLD_FLOWS
              tmpPrimExch[off].A = SphP[place].A;
#endif
#if defined(USE_ENTROPY_FOR_COLD_FLOWS) || (defined(TRACER_MC) && defined(TWODIMS)) || defined(COSMIC_RAYS_SHOCK_ACCELERATION)
              tmpPrimExch[off].Mass = P[place].Mass;
#endif

              tmpPrimExch[off].OldMass = SphP[place].OldMass;
              tmpPrimExch[off].SurfaceArea = SphP[place].SurfaceArea;
              tmpPrimExch[off].ActiveArea = SphP[place].ActiveArea;
              tmpPrimExch[off].TimeBinHydro = P[place].TimeBinHydro;

#ifdef MAXSCALARS
              for(j = 0; j < N_Scalar; j++)
                tmpPrimExch[off].Scalars[j] = *(MyFloat *) (((char *) (&SphP[place])) + scalar_elements[j].offset);
#endif
#ifdef GFM_CHEMTAGS
              for(j = 0; j < GFM_N_CHEM_TAGS; j++)
                tmpPrimExch[off].MassMetalsChemTagsFraction[j] = SphP[place].MassMetalsChemTagsFraction[j];
#endif

#ifdef TRACER_FIELD
              tmpPrimExch[off].Tracer = SphP[place].Tracer;
#endif
#if defined(TRACER_MC) && defined(TWODIMS)
              tmpPrimExch[off].tracerMC_num = get_number_of_tracers(place);
#endif
#ifdef ATOMIC_DM
              tmpPrimExch[off].Ne = SphP[place].Ne;
#endif


#ifdef FLD
              tmpPrimExch[off].n_gamma = SphP[place].n_gamma;
              tmpPrimExch[off].Kappa_diff = SphP[place].Kappa_diff;
              tmpPrimExch[off].R2 = SphP[place].R2;
#endif

              tmpPrimExch[off].TimeLastPrimUpdate = SphP[place].TimeLastPrimUpdate;

#if (defined(SHOCK_FINDER_POST_PROCESSING) || defined(SHOCK_FINDER_ON_THE_FLY) || defined(SHOCK_FINDER_BEFORE_OUTPUT)) && defined(USE_SFR)
              tmpPrimExch[off].Sfr = SphP[place].Sfr;
#endif

#ifdef SHOCK_FINDER_ON_THE_FLY
              tmpPrimExch[off].Divvel = SphP[place].Divvel;
              tmpPrimExch[off].ShockZone = SphP[place].ShockZone;
#ifdef COSMIC_RAYS
              tmpPrimExch[off].CRpseudoTemperature = SphP[place].CRpseudoTemperature;
#else
              tmpPrimExch[off].Temperature = SphP[place].Temperature;
#endif
#endif
              for(j = 0; j < 3; j++)
                {
                  tmpPrimExch[off].VelGas[j] = P[place].Vel[j];
#if defined(USE_ENTROPY_FOR_COLD_FLOWS) || defined(MHD_CT) || defined(REGULARIZE_MESH_LLOYD) || defined(REGULARIZE_MESH_SMOOTH)
                  tmpPrimExch[off].VelVertex[j] = SphP[place].VelVertex[j];
#endif
                  tmpPrimExch[off].Center[j] = SphP[place].Center[j];
#ifdef OUTPUT_CELL_SPIN
                  tmpPrimExch[off].CenterOffset[j] = SphP[place].CenterOffset[j];
#endif
#ifdef ACTIVE_CELL_SPIN
                  tmpPrimExch[off].Omega[j] = SphP[place].Omega[j];
#endif
#ifdef SHOCK_FINDER_ON_THE_FLY
                  tmpPrimExch[off].ShockDir[j] = SphP[place].ShockDir[j];
#endif
                }
              tmpPrimExch[off].Csnd = get_sound_speed(place);
            }
          listp = ListExports[listp].nextexport;
        }
    }


  /* exchange data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Mesh_Send_count[recvTask] > 0 || Mesh_Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&tmpPrimExch[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask]
                           * sizeof(struct primexch), MPI_BYTE, recvTask, TAG_DENS_A,
                           &PrimExch[Mesh_Recv_offset[recvTask]], Mesh_Recv_count[recvTask] * sizeof(struct primexch), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  myfree(tmpPrimExch);

  TIMER_STOP(CPU_MESH_EXCHANGE);
}

void exchange_primitive_variables_and_gradients(void)
{
  if(All.TotNumGas == 0)
    return;

  TIMER_START(CPU_MESH_EXCHANGE);

  int listp;
  struct grad_data *tmpGradExch;
  struct primexch *tmpPrimExch;
#ifdef TVD_SLOPE_LIMITER
  struct grad_data *tmpGradExchUl;
#endif
  int i, j, p, task, off;
  int ngrp, recvTask, place;
#if defined(DVR_RENDER) || defined(TGCHEM) || defined(MHD_CT)
  int k;
#endif

  tmpPrimExch = (struct primexch *) mymalloc("tmpPrimExch", Mesh_nexport * sizeof(struct primexch));
  tmpGradExch = (struct grad_data *) mymalloc("tmpGradExch", Mesh_nexport * sizeof(struct grad_data));
#ifdef TVD_SLOPE_LIMITER
  tmpGradExchUl = (struct grad_data *) mymalloc("tmpGradExchUl", Mesh_nexport * sizeof(struct grad_data));
#endif


  /* prepare data for export */
  for(j = 0; j < NTask; j++)
    Mesh_Send_count[j] = 0;

  for(i = 0; i < NumGasInMesh; i++)
    {
      p = List_InMesh[i];

      /* in case previous steps already lowered the Mass, update OldMass to yield together with metallicity vector conservative estimate of metal mass of each species contained in cell */
      if(P[p].Mass < SphP[p].OldMass)
        SphP[p].OldMass = P[p].Mass;


      listp = List_P[p].firstexport;
      while(listp >= 0)
        {
          if((task = ListExports[listp].origin) != ThisTask)
            {
              place = ListExports[listp].index;
              off = Mesh_Send_offset[task] + Mesh_Send_count[task]++;

              tmpPrimExch[off].Volume = SphP[place].Volume;
              tmpPrimExch[off].Density = SphP[place].Density;
#if defined(DEREFINE_GENTLY) || defined(CONDUCTION_SATURATION) || (defined(COSMIC_RAYS) && defined(SHOCK_FINDER_ON_THE_FLY) || defined(SHOCK_FINDER_BEFORE_OUTPUT)) || defined(SPECIAL_RELATIVITY) || defined(GENERAL_RELATIVITY) || defined(NON_LINEAR_SLOPE_LIMITERS) || defined(CALCULATE_QUANTITIES_IN_POSTPROCESS)
              tmpPrimExch[off].Utherm = SphP[place].Utherm;
#endif

#ifdef DVR_RENDER
              for(k = 0; k < DVR_NUM_FIELDS; k++)
                tmpPrimExch[off].DvrFields[k] = SphP[place].DvrFields[k];
#endif
              tmpPrimExch[off].Pressure = SphP[place].Pressure;

#ifdef COSMIC_RAYS
              tmpPrimExch[off].CR_Pressure = SphP[place].CR_Pressure;
              tmpPrimExch[off].CR_SpecificEnergy = SphP[place].CR_SpecificEnergy;
#ifdef COSMIC_RAYS_STREAMING
              tmpPrimExch[off].CR_Chi = SphP[place].CR_Chi;
#endif
#endif

#ifdef TGCHEM
              for(k = 0; k < TGCHEM_NUM_ABUNDANCES; k++)
                tmpPrimExch[off].Abund[k] = SphP[place].Abund[k];

              tmpPrimExch[off].Gamma = SphP[place].Gamma;
#endif


#ifdef MRT_LSF_GRADIENTS
	      tmpPrimExch[off].HI = SphP[place].HI ;
	      tmpPrimExch[off].HII = SphP[place].HII ;
	      tmpPrimExch[off].Ne = SphP[place].Ne ;
#ifdef MRT_INCLUDE_HE
	      tmpPrimExch[off].HeI = SphP[place].HeI ;
	      tmpPrimExch[off].HeII = SphP[place].HeII ;
	      tmpPrimExch[off].HeIII = SphP[place].HeIII ;
#endif

	  
#endif

#ifdef MHD
              tmpPrimExch[off].B[0] = SphP[place].B[0];
              tmpPrimExch[off].B[1] = SphP[place].B[1];
              tmpPrimExch[off].B[2] = SphP[place].B[2];
#ifdef MHD_POWELL
              tmpPrimExch[off].DivB = SphP[place].DivB;
#endif
#ifdef MHD_CT
              for(k = 0; k < 3; k++)
                {
                  tmpPrimExch[off].A[k] = SphP[place].A[k];
                }
              tmpPrimExch[off].TimeLastBUpdate = SphP[place].TimeLastBUpdate;
#endif
#ifdef NON_IDEAL_MHD
#if defined(OHMIC_DIFFUSION) || defined(IMPLICIT_OHMIC_DIFFUSION) || defined(AMBIPOLAR_DIFFUSION)
              tmpPrimExch[off].CurlB[0] = SphP[place].CurlB[0];
              tmpPrimExch[off].CurlB[1] = SphP[place].CurlB[1];
              tmpPrimExch[off].CurlB[2] = SphP[place].CurlB[2];
#endif
#endif
#endif
#ifdef VARIABLE_GAMMA
              tmpPrimExch[off].GammaC = SphP[place].GammaC;
              tmpPrimExch[off].GammaE = SphP[place].GammaE;
#endif
#ifdef USE_ENTROPY_FOR_COLD_FLOWS
              tmpPrimExch[off].A = SphP[place].A;
#endif
#if defined(USE_ENTROPY_FOR_COLD_FLOWS) || (defined(TRACER_MC) && defined(TWODIMS)) || defined(COSMIC_RAYS_SHOCK_ACCELERATION)
              tmpPrimExch[off].Mass = P[place].Mass;
#endif
              tmpPrimExch[off].OldMass = SphP[place].OldMass;
              tmpPrimExch[off].SurfaceArea = SphP[place].SurfaceArea;
              tmpPrimExch[off].ActiveArea = SphP[place].ActiveArea;

              tmpPrimExch[off].TimeBinHydro = P[place].TimeBinHydro;

#ifdef MAXSCALARS
              for(j = 0; j < N_Scalar; j++)
                tmpPrimExch[off].Scalars[j] = *(MyFloat *) (((char *) (&SphP[place])) + scalar_elements[j].offset);
#endif
#ifdef GFM_CHEMTAGS
              for(j = 0; j < GFM_N_CHEM_TAGS; j++)
                tmpPrimExch[off].MassMetalsChemTagsFraction[j] = SphP[place].MassMetalsChemTagsFraction[j];
#endif
#ifdef TRACER_FIELD
              tmpPrimExch[off].Tracer = SphP[place].Tracer;
#endif
#if defined(TRACER_MC) && defined(TWODIMS)
              tmpPrimExch[off].tracerMC_num = get_number_of_tracers(place);
#endif
#ifdef ATOMIC_DM
              tmpPrimExch[off].Ne = SphP[place].Ne;
#endif

#ifdef FLD
              tmpPrimExch[off].n_gamma = SphP[place].n_gamma;
              tmpPrimExch[off].Kappa_diff = SphP[place].Kappa_diff;
              tmpPrimExch[off].R2 = SphP[place].R2;
#endif

              tmpPrimExch[off].TimeLastPrimUpdate = SphP[place].TimeLastPrimUpdate;

#if (defined(SHOCK_FINDER_POST_PROCESSING) || defined(SHOCK_FINDER_ON_THE_FLY) || defined(SHOCK_FINDER_BEFORE_OUTPUT)) && defined(USE_SFR)
              tmpPrimExch[off].Sfr = SphP[place].Sfr;
#endif

#ifdef SHOCK_FINDER_ON_THE_FLY
              tmpPrimExch[off].Divvel = SphP[place].Divvel;
              tmpPrimExch[off].ShockZone = SphP[place].ShockZone;
#ifdef COSMIC_RAYS
              tmpPrimExch[off].CRpseudoTemperature = SphP[place].CRpseudoTemperature;
#else
              tmpPrimExch[off].Temperature = SphP[place].Temperature;
#endif
#endif
              for(j = 0; j < 3; j++)
                {
                  tmpPrimExch[off].VelGas[j] = P[place].Vel[j];
                  tmpPrimExch[off].Center[j] = SphP[place].Center[j];
                  tmpPrimExch[off].VelVertex[j] = SphP[place].VelVertex[j];
#ifdef OUTPUT_CELL_SPIN
                  tmpPrimExch[off].CenterOffset[j] = SphP[place].CenterOffset[j];
#endif
#ifdef ACTIVE_CELL_SPIN
                  tmpPrimExch[off].Omega[j] = SphP[place].Omega[j];
#endif
#ifdef SHOCK_FINDER_ON_THE_FLY
                  tmpPrimExch[off].ShockDir[j] = SphP[place].ShockDir[j];
#endif
                }

              tmpGradExch[off] = SphP[place].Grad;

#ifdef TVD_SLOPE_LIMITER
              tmpGradExchUl[off] = SphP[place].GradUl;
#endif

              tmpPrimExch[off].Csnd = get_sound_speed(place);
            }
          listp = ListExports[listp].nextexport;
        }
    }


  /* exchange data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Mesh_Send_count[recvTask] > 0 || Mesh_Recv_count[recvTask] > 0)
            {
              /* exchange the data */
              MPI_Sendrecv(&tmpPrimExch[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask]
                           * sizeof(struct primexch), MPI_BYTE, recvTask, TAG_DENS_A,
                           &PrimExch[Mesh_Recv_offset[recvTask]], Mesh_Recv_count[recvTask] * sizeof(struct primexch), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

              MPI_Sendrecv(&tmpGradExch[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask]
                           * sizeof(struct grad_data), MPI_BYTE, recvTask, TAG_HYDRO_A,
                           &GradExch[Mesh_Recv_offset[recvTask]], Mesh_Recv_count[recvTask] * sizeof(struct grad_data), MPI_BYTE, recvTask, TAG_HYDRO_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

#ifdef TVD_SLOPE_LIMITER
              MPI_Sendrecv(&tmpGradExchUl[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask]
                           * sizeof(struct grad_data), MPI_BYTE, recvTask, TAG_HYDRO_A,
                           &GradExchUl[Mesh_Recv_offset[recvTask]], Mesh_Recv_count[recvTask] * sizeof(struct grad_data), MPI_BYTE, recvTask, TAG_HYDRO_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#endif
            }
        }
    }

#ifdef TVD_SLOPE_LIMITER
  myfree(tmpGradExchUl);
#endif
  myfree(tmpGradExch);
  myfree(tmpPrimExch);

  TIMER_STOP(CPU_MESH_EXCHANGE);

  /* note: because the sequence is the same as before, we don't have to do the sorts again */
}

int compare_primexch(const void *a, const void *b)
{
  if(((struct data_primexch_compare *) a)->task < ((struct data_primexch_compare *) b)->task)
    return -1;

  if(((struct data_primexch_compare *) a)->task > ((struct data_primexch_compare *) b)->task)
    return +1;

  if(((struct data_primexch_compare *) a)->index < ((struct data_primexch_compare *) b)->index)
    return -1;

  if(((struct data_primexch_compare *) a)->index > ((struct data_primexch_compare *) b)->index)
    return +1;

  return 0;
}


#ifdef OUTPUT_VERTEX_VELOCITY_DIVERGENCE
void voronoi_update_ghost_velvertex(void)
{
  CPU_Step[CPU_MISC] += measure_time();

  int listp;
  int i, j, p, task, off;
  int ngrp, recvTask, place;
  struct velvertex_data
  {
    MyFloat VelVertex[3];
  }
   *tmpVelVertexExch, *tmpVelVertexRecv;

  tmpVelVertexExch = (struct velvertex_data *) mymalloc("tmpVelVertexExch", Mesh_nexport * sizeof(struct velvertex_data));

  /* prepare data for export */
  for(j = 0; j < NTask; j++)
    Mesh_Send_count[j] = 0;

  for(i = 0; i < NumGasInMesh; i++)
    {
      p = List_InMesh[i];

      listp = List_P[p].firstexport;
      while(listp >= 0)
        {
          if((task = ListExports[listp].origin) != ThisTask)
            {
              place = ListExports[listp].index;
              off = Mesh_Send_offset[task] + Mesh_Send_count[task]++;

              for(j = 0; j < 3; j++)
                {
                  tmpVelVertexExch[off].VelVertex[j] = SphP[place].VelVertex[j];
                }
            }
          listp = ListExports[listp].nextexport;
        }
    }

  /* exchange data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Mesh_Send_count[recvTask] > 0 || Mesh_Recv_count[recvTask] > 0)
            {
              tmpVelVertexRecv = (struct velvertex_data *) mymalloc("tmpVelVertexRecv", Mesh_Recv_count[recvTask] * sizeof(struct velvertex_data));

              /* get the values */
              MPI_Sendrecv(&tmpVelVertexExch[Mesh_Send_offset[recvTask]],
                           Mesh_Send_count[recvTask] * sizeof(struct velvertex_data), MPI_BYTE,
                           recvTask, TAG_DENS_A, tmpVelVertexRecv, Mesh_Recv_count[recvTask] * sizeof(struct velvertex_data), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

              for(i = 0; i < Mesh_Recv_count[recvTask]; i++)
                {
                  for(j = 0; j < 3; j++)
                    {
                      PrimExch[Mesh_Recv_offset[recvTask] + i].VelVertex[j] = tmpVelVertexExch[i].VelVertex[j];
                    }
                }

              myfree(tmpVelVertexRecv);
            }
        }
    }

  myfree(tmpVelVertexExch);

  CPU_Step[CPU_SET_VERTEXVELS] += measure_time();
}
#endif


#ifdef SECOND_DERIVATIVES
void voronoi_exchange_hessians(void)
{
  CPU_Step[CPU_MISC] += measure_time();
  int listp;
  struct hessian_data *tmpHessianExch;
  int j, p, task, off;
  int ngrp, recvTask, place;

  tmpHessianExch = (struct hessian_data *) mymalloc("tmpHessianExch", Mesh_nexport * sizeof(struct hessian_data));

  /* prepare data for export */
  for(j = 0; j < NTask; j++)
    Mesh_Send_count[j] = 0;

  for(i = 0; i < NumGasInMesh; i++)
    {
      p = List_InMesh[i];

      listp = List_P[p].firstexport;
      while(listp >= 0)
        {
          if((task = ListExports[listp].origin) != ThisTask)
            {
              place = ListExports[listp].index;
              off = Mesh_Send_offset[task] + Mesh_Send_count[task]++;

              tmpHessianExch[off] = SphP[place].Hessian;
            }
          listp = ListExports[listp].nextexport;
        }
    }

  /* exchange data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Mesh_Send_count[recvTask] > 0 || Mesh_Recv_count[recvTask] > 0)
            {
              /* exchange the data */

              MPI_Sendrecv(&tmpHessianExch[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask]
                           * sizeof(struct hessian_data), MPI_BYTE, recvTask, TAG_HYDRO_A,
                           &HessianExch[Mesh_Recv_offset[recvTask]], Mesh_Recv_count[recvTask] * sizeof(struct hessian_data), MPI_BYTE, recvTask, TAG_HYDRO_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  myfree(tmpHessianExch);

  CPU_Step[CPU_MESH_EXCHANGE] += measure_time();
}
#endif

#endif /* VORONOI */
