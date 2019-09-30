/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/tracer_mc.c
 * \date        07/2012
 * \author      Shy Genel, Dylan Nelson, Mark Vogelsberger
 * \brief       Monte Carlo Tracer Particles
 * \details     
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>

#include "allvars.h"
#include "proto.h"
#include "voronoi.h"

#ifdef TRACER_MC

/*! \brief Process exchange of tracers requiring remote communication.
 *  Standard MPI communication technique.
 */
int exchange_tracer_flux_list(void)
{
  int i, j, nimport, nexport, ngrp, recvTask;

  for(j = 0; j < NTask; j++)
    Send_count[j] = 0;

  for(i = 0; i < Ntracerflux; i++)
    Send_count[TracerFluxListIn[i].task]++;

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, nexport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nimport += Recv_count[j];
      nexport += Send_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  struct tracer_flux_list_data *TracerFluxListSend = (struct tracer_flux_list_data *) mymalloc("TracerFluxListSend", nexport * sizeof(struct tracer_flux_list_data));
  TracerFluxListGet = (struct tracer_flux_list_data *) mymalloc("TracerFluxListGet", nimport * sizeof(struct tracer_flux_list_data));

  for(j = 0; j < NTask; j++)
    Send_count[j] = 0;

  for(i = 0; i < Ntracerflux; i++)
    {
      int task = TracerFluxListIn[i].task;
      int ind = Send_offset[task] + Send_count[task]++;

      TracerFluxListSend[ind] = TracerFluxListIn[i];
    }

  /* exchange particle data */
  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&TracerFluxListSend[Send_offset[recvTask]],
                           Send_count[recvTask] * sizeof(struct tracer_flux_list_data), MPI_BYTE,
                           recvTask, TAG_DENS_A,
                           &TracerFluxListGet[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct tracer_flux_list_data), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  /* process imported tracers */
  for(i = 0; i < nimport; i++)
    {
      int p = TracerFluxListGet[i].index;

      int itr = get_free_tracer_slot();

      TracerLinkedList[itr].ID = TracerFluxListGet[i].ID;
#ifdef TRACER_MC_NUM_FLUID_QUANTITIES
      for(int j = 0; j < TRACER_MC_NUM_FLUID_QUANTITIES; j++)
        TracerLinkedList[itr].fluid_quantities[j] = TracerFluxListGet[i].fluid_quantities[j];
#if (TRACER_MC_EXCHANGE_DISTANCE)
      if(P[p].Type == 0)
        TracerLinkedList[itr].fluid_quantities[TracerMCExchangeDistanceIndex] += get_cell_radius(p);
#endif
#if (TRACER_MC_EXCHANGE_DISTANCE_ERROR)
      if(P[p].Type == 0)
        {
          TracerLinkedList[itr].fluid_quantities[TracerMCExchangeDistanceErrorIndex] +=
            get_cell_radius(p) * (sqrt(TracerLinkedList[itr].fluid_quantities[TracerMCExchangeCounterIndex]) - sqrt(TracerLinkedList[itr].fluid_quantities[TracerMCExchangeCounterIndex] - 1));
        }
#endif
#endif
#ifdef TRACER_MC_CHECKS
      TracerLinkedList[itr].ParentID = TracerFluxListGet[i].ParentID;
#endif

      if(p < 0 || p >= NumPart)
        terminate("p=%d NumPart=%d\n", p, NumPart);

      add_tracer_to_parent(p, itr);
    }

  myfree(TracerFluxListGet);
  myfree(TracerFluxListSend);

  return nimport;
}

/*! \brief Prepare for tracer exchange
 *  
 *  \param nexport_tracer room required in export buffer, set to N_tracer if unknown
 */
void start_MC_tracer(int nexport_tracer)
{
  /* count how many tracers need to be exported */
  Ntracerflux = 0;

  /* export buffer */
  TracerFluxListIn = mymalloc("TracerFluxListIn", nexport_tracer * sizeof(struct tracer_flux_list_data));
}

/*! \brief Carry out tracer exchange and attach to new local parents. */
void finish_MC_tracer(void)
{
  exchange_tracer_flux_list();

  myfree(TracerFluxListIn);

  TracerFluxListIn = NULL;
}

/*! \brief Add a single tracer to the TracerFluxList for remote exchange, and remove it from the 
 *         local TracerLinkedList at the same time.
 *
 * \param p_orig P_index of original parent (use -1 to indicate it was unattached previously).
 * \param tr_ind Index in TracerLinkedList of the tracer to move.
 * \param pother_task Task number where pother_ID potential parent exists.
 * \param pother_index Pindex of the potential parent on the task given by pother_task.
 * \param pother_ID Particle/cell ID of the potential parent.
 */
void add_tracer_to_TFL(int p_orig, int tr_ind, int pother_task, int pother_index, MyIDType pother_ID)
{
  /* fill export structure */
  TracerFluxListIn[Ntracerflux].task = pother_task;
  TracerFluxListIn[Ntracerflux].index = pother_index;
#ifdef TRACER_MC_CHECKS
  TracerFluxListIn[Ntracerflux].ParentID = pother_ID;
#endif

  TracerFluxListIn[Ntracerflux].ID = TracerLinkedList[tr_ind].ID;

#ifdef TRACER_MC_NUM_FLUID_QUANTITIES
  for(int j = 0; j < TRACER_MC_NUM_FLUID_QUANTITIES; j++)
    TracerFluxListIn[Ntracerflux].fluid_quantities[j] = TracerLinkedList[tr_ind].fluid_quantities[j];
#endif

  /* remove from existing parent, delete entry in TLL */
  if(p_orig >= 0)
    remove_tracer_from_parent(p_orig, tr_ind);

  release_tracer_slot(tr_ind);

  Ntracerflux++;
}

/*! \brief Start (allocate for) re-attachment of snapshot loaded MC tracers to their parents. */
void load_MC_tracer_start(void)
{
  if(RestartFlag < 2)
    return;

  // allocate temporary buffer to hold ParentIDs for reattachment
  tracer_cellids = (MyIDType *) mymalloc("tracer_cellids", sizeof(MyIDType) * All.MaxPartTracer);
  memset(tracer_cellids, 0, sizeof(MyIDType) * All.MaxPartTracer);
}

/*! \brief Sort comparator function for load_MC_tracer_reattach(). */
int P_ID_cmp(const void *a, const void *b)
{
  int ia = *(int *) a;
  int ib = *(int *) b;

  return P[ia].ID < P[ib].ID ? -1 : P[ia].ID > P[ib].ID;
}

/*! \brief Binary search helper for load_MC_tracer_reattach().
 * \param sorted_ids Pointer to array of MyIDType entries, must be pre-sorted in ascending order.
 * \param N Size of sorted_ids.
 * \param search_id Single value to locate within sorted_ids array. If an exact match is found, 
 *   the index in sorted_ids is returned, otherwise -1 is returned.
 */
inline int bsearch_tracer(MyIDType *sorted_ids, int N, MyIDType search_id)
{
  int L = 0;
  int R = N - 1;
  int m = 0;

  while(L <= R)
    {
      m = floor(0.5*(L+R));

      if(sorted_ids[m] == search_id)
        break;
      if(sorted_ids[m] < search_id)
        L = m + 1;
      if(sorted_ids[m] > search_id)
        R = m - 1;
    }

  if(sorted_ids[m] == search_id)
    return m;

  return -1;
}

/*! \brief Carry out (possibly parallel) re-attachment of snapshot loaded MC tracers to their parents.
 */
void load_MC_tracer_reattach(void)
{
  if(RestartFlag < 2)
    return;

  mpi_printf("TRACER_MC: load_MC_tracer_reattach()...\n");

  int i, k, sorted_par_index, ngrp, recvTask;
  int count_nonlocal = 0, count_incoming = 0, max_count_nonlocal = 0;

  int *P_ID_indices, *ParentTask, *ParentIndex, *ParentIndex_Recv, *ParentIndex_Remote;
  MyIDType *P_ID_sorted, *ParentID_Recv;

  // (1) allocate and init
  P_ID_indices = (int *) mymalloc("P_ID_indices", sizeof(int) * NumPart);
  P_ID_sorted  = (MyIDType *) mymalloc("P_ID_sorted", sizeof(MyIDType) * NumPart);

  for(i = 0; i < NumPart; i++)
    {
      P[i].TracerHead = -1;
      P[i].NumberOfTracers = 0;
      P_ID_indices[i] = i;
    }

  // (2) get sorted list of P[].IDs and sort indices
  qsort(P_ID_indices, NumPart, sizeof(int), P_ID_cmp);

  for(i = 0; i < NumPart; i++)
    P_ID_sorted[i] = P[P_ID_indices[i]].ID;

  // (3) attach those with local parents, use TTL.Next == -2 as a flag for those remaining unattached
  for(i = 0; i < N_tracer; i++)
    {
      int sorted_par_index = bsearch_tracer(P_ID_sorted, NumPart, tracer_cellids[i]);
      if(sorted_par_index >= 0)
        {
          add_tracer_to_parent(P_ID_indices[sorted_par_index], i);
        }
      else
        {
          TracerLinkedList[i].Next = -2;
          count_nonlocal++;
        }
    }

  // (4) compact tracer_cellids to those which are still unattached
  for(i = 0, k = 0; i < N_tracer; i++)
    if(TracerLinkedList[i].Next == -2)
      tracer_cellids[k++] = tracer_cellids[i];

  // (5) determine maximum count_nonlocal among tasks, allocate arrays for communication
  memset(Send_count, 0, sizeof(int) * NTask);
  memset(Recv_count, 0, sizeof(int) * NTask);

  Send_count[ThisTask] = count_nonlocal;
  MPI_Allreduce(&Send_count[0], &Recv_count[0], NTask, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  for(i = 0; i < NTask; i++)
    if(Recv_count[i] > max_count_nonlocal)
      max_count_nonlocal = Recv_count[i];

  ParentTask = (int *) mymalloc("ParentTask", sizeof(int) * count_nonlocal);
  ParentIndex = (int *) mymalloc("ParentIndex", sizeof(int) * count_nonlocal);

  ParentID_Recv = (MyIDType *) mymalloc("ParentID_Recv", sizeof(MyIDType) * max_count_nonlocal);
  ParentIndex_Recv = (int *) mymalloc("ParentIndex_Recv", sizeof(int) * max_count_nonlocal);
  ParentIndex_Remote = (int *) mymalloc("ParentIndex_Remote", sizeof(int) * max_count_nonlocal);

  for(i = 0; i < count_nonlocal; i++)
    ParentTask[i] = -1;

  // (6) hypercube swap orphaned ParentIDs
  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask >= NTask || recvTask == ThisTask)
        continue;

      // (6a) exchange ParentIDs (nonlocal)
      MPI_Sendrecv(&tracer_cellids[0],
                   Recv_count[ThisTask] * sizeof(MyIDType), // = count_nonlocal * sizeof(MyIDType)
                   MPI_BYTE,
                   recvTask, 
                   TAG_DENS_A,
                   &ParentID_Recv[0], 
                   Recv_count[recvTask] * sizeof(MyIDType), 
                   MPI_BYTE, 
                   recvTask, 
                   TAG_DENS_A, 
                   MPI_COMM_WORLD, 
                   MPI_STATUS_IGNORE);

      // (6b) check if any receieved ParentIDs are local
      for(i = 0; i < max_count_nonlocal; i++)
        ParentIndex_Remote[i] = -1;

      for(i = 0; i < Recv_count[recvTask]; i++)
      {
        sorted_par_index = bsearch_tracer(P_ID_sorted, NumPart, ParentID_Recv[i]);

        // (6c) fill ParentID_Index array with Pindex on ThisTask in these cases
        if(sorted_par_index >= 0)
          {
            ParentIndex_Remote[i] = P_ID_indices[sorted_par_index];
            count_incoming++;
          }
      }

      // (6d) swap ParentID_Task/Index results
      MPI_Sendrecv(&ParentIndex_Remote[0],
                   Recv_count[recvTask] * sizeof(int),
                   MPI_BYTE,
                   recvTask, 
                   TAG_DENS_A,
                   &ParentIndex_Recv[0], 
                   Recv_count[ThisTask] * sizeof(int), // = count_nonlocal * sizeof(int)
                   MPI_BYTE, 
                   recvTask, 
                   TAG_DENS_A, 
                   MPI_COMM_WORLD, 
                   MPI_STATUS_IGNORE);

      // (6e) check received task/index array for valid entries and save them
      for(i = 0; i < count_nonlocal; i++)
        if(ParentIndex_Recv[i] >= 0)
          {
            if(ParentTask[i] != -1)
              terminate("[%d] TRACER_MC: attempt to overwrite ParentID_Task[%d] = %d to %d\n",
                ThisTask, i, ParentTask[i], recvTask);

            ParentTask[i] = recvTask;
            ParentIndex[i] = ParentIndex_Recv[i];
          }
    }

  myfree(ParentIndex_Remote);
  myfree(ParentIndex_Recv);
  myfree(ParentID_Recv);

  // (7) initalize a TFL, populate, and exchange
  int N_tracer_orig = N_tracer; /* decremented when we add to TFL */

  start_MC_tracer(count_nonlocal);

  for(i = 0, k = 0; i < N_tracer_orig; i++)
    {
      if(TracerLinkedList[i].Next != -2)
        continue; /* skip loaded tracers which were attached to local parents */

      if(ParentTask[k] < 0 || ParentTask[k] >= NTask || ParentTask[k] == ThisTask)
        terminate("[%d] TRACER_MC: nonlocal[i=%d] has ParentTask = %d ParentID=%lld\n",
          ThisTask, k, ParentTask[k], (long long)tracer_cellids[k]);

#ifdef TRACER_MC_CHECKS
      if(tracer_cellids[k] != TracerLinkedList[i].ParentID)
        terminate("[%d] TRACER_MC: mismatch in export i=%d k=%d\n", ThisTask, i, k);
#endif

      add_tracer_to_TFL(-1, i, ParentTask[k], ParentIndex[k], tracer_cellids[k]);
      k++;
    }

  finish_MC_tracer();
  
  // (8) finish and verify
  myfree(ParentIndex);
  myfree(ParentTask);
  myfree(P_ID_sorted);
  myfree(P_ID_indices);
  myfree(tracer_cellids);

#ifdef TRACER_MC_CHECKS
  check_tracer_lists();
#endif

}

/*! \brief for each tracer, update its tracked properties.
 * Fluid quantities are updated according to the current (end of timestep) 
 * values of its parent gas cell. These should be set on the quantities of interest for a given simulation.
 */
void record_tracer_parent_fluid_properties(void)
{
#ifdef TRACER_MC_NUM_FLUID_QUANTITIES
  int idx, i, parent_ind;
  double mu;
#if (TRACER_MC_TMAX) || (TRACER_MC_ENTMAX)
  double fluid_val;
#endif

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      parent_ind = TimeBinsHydro.ActiveParticleList[idx];
      if(parent_ind < 0)
        continue;

      if(P[parent_ind].Type == 0)
        {
          i = P[parent_ind].TracerHead;
          while(i != -1)
            {

              /* temperature [Kelvin] */
#ifdef COOLING
              mu = 4.0 / (1.0 + 3.0 * HYDROGEN_MASSFRAC + 4.0 * HYDROGEN_MASSFRAC * SphP[parent_ind].Ne);
#else
              mu = 4.0 / (1.0 + 3.0 * HYDROGEN_MASSFRAC);
#endif
#if (TRACER_MC_TMAX)
              fluid_val = GAMMA_MINUS1 / BOLTZMANN * SphP[parent_ind].Utherm * PROTONMASS * mu * All.UnitEnergy_in_cgs / All.UnitMass_in_g;
#ifdef USE_SFR
              if(SphP[parent_ind].Sfr == 0 && fluid_val > TracerLinkedList[i].fluid_quantities[TracerMCTmaxIndex])
#else
              if(fluid_val > TracerLinkedList[i].fluid_quantities[TracerMCTmaxIndex])
#endif
                {
                  TracerLinkedList[i].fluid_quantities[TracerMCTmaxIndex] = fluid_val;
#if (TRACER_MC_TMAX_TIME)
                  TracerLinkedList[i].fluid_quantities[TracerMCTmaxTimeIndex] = All.Time;
#endif
#if (TRACER_MC_TMAX_RHO)
                  TracerLinkedList[i].fluid_quantities[TracerMCTmaxRhoIndex] = SphP[parent_ind].Density;
#endif
                }
#endif
              /* density [code units] */
#if (TRACER_MC_RHOMAX)
              if(SphP[parent_ind].Density > TracerLinkedList[i].fluid_quantities[TracerMCRhomaxIndex])
                {
                  TracerLinkedList[i].fluid_quantities[TracerMCRhomaxIndex] = SphP[parent_ind].Density;
#if (TRACER_MC_RHOMAX_TIME)
                  TracerLinkedList[i].fluid_quantities[TracerMCRhomaxTimeIndex] = All.Time;
#endif
                }
#endif
              /* mach number from Riemann solver */
#if (TRACER_MC_MACHMAX)
              if(SphP[parent_ind].MaxMach > TracerLinkedList[i].fluid_quantities[TracerMCMachmaxIndex])
                TracerLinkedList[i].fluid_quantities[TracerMCMachmaxIndex] = SphP[parent_ind].MaxMach;
#endif
              /* entropy (note units - physical, not comoving) [(P/a^3)/(rho/a^3)^gamma] */
#if (TRACER_MC_ENTMAX)
              fluid_val = (SphP[parent_ind].Pressure * All.cf_a3inv) / pow(SphP[parent_ind].Density * All.cf_a3inv, GAMMA);
              if(fluid_val > TracerLinkedList[i].fluid_quantities[TracerMCEntmaxIndex])
                {
                  TracerLinkedList[i].fluid_quantities[TracerMCEntmaxIndex] = fluid_val;
#if (TRACER_MC_ENTMAX_TIME)
                  TracerLinkedList[i].fluid_quantities[TracerMCEntmaxTimeIndex] = All.Time;
#endif
                }
#endif

              /* mach number from on-the-fly shock finder */
#if (TRACER_MC_SHOCKMACHNUM_MAX)
              if(SphP[parent_ind].Machnumber > TracerLinkedList[i].fluid_quantities[TracerMCShockMachMaxIndex])
                TracerLinkedList[i].fluid_quantities[TracerMCShockMachMaxIndex] = SphP[parent_ind].Machnumber;
#endif

              i = TracerLinkedList[i].Next;
            }
        }
    }
#endif
}

/*! \brief get number of tracers in cell
 *
 * \param i P_index
 */
int get_number_of_tracers(int i)
{
#ifndef TRACER_MC_CHECKS
  return P[i].NumberOfTracers;
#else
  int next = P[i].TracerHead;
  int count = 0;

  while(next >= 0)
    {
      if(P[i].ID != TracerLinkedList[next].ParentID)
        terminate("i=%d P[i].Type=%d P[i].ID=%lld TracerLinkedList[next].ParentID=%lld\n", 
          i, P[i].Type, (long long)P[i].ID, (long long)TracerLinkedList[next].ParentID);

      count++; /* increment count regardless of if ID==0 or not */
      next = TracerLinkedList[next].Next;

      if(count > N_tracer)
        terminate("[%d] TRACER_MC: count %d > N_tracer %d (i=%d Type=%d)", ThisTask, count, N_tracer, i, P[i].Type);
    }

  if(count != P[i].NumberOfTracers)
    terminate("Mismatch found. i=%d type=%d count=%d P[i].TracerHead=%d P[i].NumberOfTracers=%d", i, P[i].Type, count, P[i].TracerHead, P[i].NumberOfTracers);

  return count;
#endif
}

/*! \brief Get global maximum number of tracers attached to any one parent.
 */
int get_max_number_of_tracers()
{
  int i, tmp;
  int max_local = 0, max_global = 0;

  for(i = 0; i < NumPart; i++)
    {
      tmp = get_number_of_tracers(i);
      if(tmp > max_local)
        max_local = tmp;
    }

  MPI_Allreduce(&max_local, &max_global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  return max_global;
}

/*! \brief get global number of tracers for all particle types
 *
 *  \param flag set flag=-1 to return for all particle types (global total), 
 *         or set flag>=0 to return just for that particle type (global subtotal)
 */
long long get_total_number_of_tracers(int flag)
{
  int i;
  int num_local = 0;
  long long num_global = 0;

  for(i = 0; i < NumPart; i++)
    {
      if(flag < 0 || P[i].Type == flag)
        num_local += get_number_of_tracers(i);
    }

  sumup_large_ints(1, &num_local, &num_global);
  return num_global;
}


/*! \brief restore connectivity of the tracer linked list broken by collective subfind 
 *  during particle unbinding
 */
void restore_tracer_connectivity(void)
{
  int i;

  for(i = 0; i < NumPart; i++)
    {
      int next = P[i].TracerHead;
      while(next >= 0)
        {
          int q = TracerLinkedList[next].Next;

          if(q < 0)
            break;

          TracerLinkedList[q].Prev = next;
          next = q;
        }
    }
}

/*! \brief Monte Carlo exchange of tracers from a local cell/particle to a particle/cell that may be local or remote.
 *
 * \param p P_index of original parent.
 * \param pother_task Task number where pother_ID potential parent exists.
 * \param pother_index Pindex of the potential parent on the task given by pother_task.
 * \param pother_ID Particle/cell ID of the potential parent.
 * \param prob probability whether or not to move each child tracer.
 */
int consider_moving_tracers(int p, int pother_task, int pother_index, MyIDType pother_ID, double prob)
{
  int nmoved = 0;
  int next = P[p].TracerHead;

  if(prob < 0.0 || prob > 1.0)
    terminate("TRACER_MC: bad prob=%g\n", prob);

  while(next >= 0)
    {
      int move = next;
      next = TracerLinkedList[next].Next;

      double rndnum = get_random_number();

      if(rndnum < prob)
        {

#ifdef TRACER_MC_NUM_FLUID_QUANTITIES
#if (TRACER_MC_EXCHANGE_COUNTER)
          TracerLinkedList[move].fluid_quantities[TracerMCExchangeCounterIndex] += 1;
#endif
#if (TRACER_MC_LAST_STAR_TIME)
          if(P[p].Type == 4)    /* moving from a star or wind particle */
            {
              if(StarP[P[p].AuxDataID].BirthTime == 0)  /* recoupling wind */
                TracerLinkedList[move].fluid_quantities[TracerMCLastStarTimeIndex] = -All.Time;
              else if(StarP[P[p].AuxDataID].BirthTime > 0)      /* stellar mass loss */
                TracerLinkedList[move].fluid_quantities[TracerMCLastStarTimeIndex] = All.Time;
              else              /* if BirthTime < 0 */
                terminate("returning tracer from active wind particle to gas: should not happen");
            }
#endif
#if (TRACER_MC_EXCHANGE_DISTANCE)
          if(P[p].Type == 0)
            TracerLinkedList[move].fluid_quantities[TracerMCExchangeDistanceIndex] += get_cell_radius(p);
#endif
#if (TRACER_MC_EXCHANGE_DISTANCE_ERROR)
          if(P[p].Type == 0)
            {
              TracerLinkedList[move].fluid_quantities[TracerMCExchangeDistanceErrorIndex] +=
                get_cell_radius(p) *
                (sqrt(TracerLinkedList[move].fluid_quantities[TracerMCExchangeCounterIndex]) - 
                 sqrt(TracerLinkedList[move].fluid_quantities[TracerMCExchangeCounterIndex] - 1));
            }
#endif
#endif

          add_tracer_to_TFL(p, move, pother_task, pother_index, pother_ID);
          nmoved++;
        }
    }

  return nmoved;
}

/*! \brief Monte Carlo exchange of exactly one tracer, chosen at random from a specified
 *  parent particle/cell. Note: Destination cell can be either local or remote.
 *
 * \param p P_index of original parent.
 * \param pother_task Task number where pother_ID potential parent exists.
 * \param pother_index Pindex of the potential parent on the task given by pother_task.
 * \param pother_ID Particle/cell ID of the potential parent.
 */
void move_one_tracer(int p, int pother_task, int pother_index, MyIDType pother_ID)
{
  int next = P[p].TracerHead;

  int id_to_move = P[p].NumberOfTracers * get_random_number();
  int cnt = 0;

  while(next >= 0)
    {
      int move = next;
      next = TracerLinkedList[next].Next;

      if(cnt == id_to_move)
        {

#ifdef TRACER_MC_NUM_FLUID_QUANTITIES
#if (TRACER_MC_EXCHANGE_COUNTER)
          TracerLinkedList[move].fluid_quantities[TracerMCExchangeCounterIndex] += 1;
#endif
#if (TRACER_MC_LAST_STAR_TIME)
          if(P[p].Type == 4)    /* moving from a star or wind particle */
            {
              if(StarP[P[p].AuxDataID].BirthTime > 0)   /* stellar mass loss */
                TracerLinkedList[move].fluid_quantities[TracerMCLastStarTimeIndex] = All.Time;
              /* if source is a wind particle, this is because of chemical stripping,
                 in which case the tracer should not be considered to have come from a wind,
                 since the mass is stripped at the same time as the wind launch itself.
                 hence, we restore its original TracerMCLastStarTimeIndex value */
              else if(StarP[P[p].AuxDataID].BirthTime < 0)
                TracerLinkedList[move].fluid_quantities[TracerMCLastStarTimeIndex] -= 3 * All.TimeMax;
              else              /* if BirthTime == 0 */
                terminate("moving tracer from a stellar particle with age=0, should be done while recoupling, not here");
            }
#endif
#endif

          add_tracer_to_TFL(p, move, pother_task, pother_index, pother_ID);
          break;
        }

      cnt++;
    }
}

/*! \brief Monte Carlo exchange of tracers between two cells/particles that are certainly local. 
 * 
 *  \param p_from P_index to move tracers from
 *  \param p_to P_index to move tracers to
 *  \param prob probability whether or not to move each child tracer
 */
int consider_moving_tracers_local(int p_from, int p_to, double prob)
{
  int Nmoved = 0;
  int next = P[p_from].TracerHead;

  if(prob < 0.0 || prob > 1.0)
    terminate("TRACER_MC: bad prob=%g\n", prob);

  while(next >= 0)
    {
      int move = next;
      next = TracerLinkedList[next].Next;

      double rndnum = get_random_number();

      if(rndnum < prob)         /* move tracer */
        {
          move_tracer_between_parents(p_from, p_to, move);
          Nmoved++;
        }
    }

  return Nmoved;
}

/*! \brief Return free index in TLL from the heap and increment N_tracer. */
int get_free_tracer_slot(void)
{
  if(N_tracer >= All.MaxPartTracer)
    terminate("N_tracer=%d >= All.MaxPartTracer=%d", N_tracer, All.MaxPartTracer);

  return TracerLinkedListHeap[N_tracer++];
}

/*! \brief Free an entry in the TLL by adding its index to the heap, and decrement N_tracer.
 * 
 * \param itr Index of TLL to free. 
 */
void release_tracer_slot(int itr)
{
  if(N_tracer <= 0 || itr < 0)
    terminate("N_tracer=%d itr=%d", N_tracer, itr);

  TracerLinkedListHeap[--N_tracer] = itr;
}


/*! \brief Remove tracer from parent. Your responsibility to add it elsewhere.
 *
 *  \param p P_index of the original parent
 *  \param itracer TLL_index of tracer to remove
 */
void remove_tracer_from_parent(int p, int itracer)
{
  if(TracerLinkedList[itracer].Prev < 0)        /* first tracer */
    P[p].TracerHead = TracerLinkedList[itracer].Next;
  else
    TracerLinkedList[TracerLinkedList[itracer].Prev].Next = TracerLinkedList[itracer].Next;

  if(TracerLinkedList[itracer].Next != -1)
    TracerLinkedList[TracerLinkedList[itracer].Next].Prev = TracerLinkedList[itracer].Prev;

#ifdef TRACER_MC_CHECKS
  if(P[p].ID != TracerLinkedList[itracer].ParentID)
    terminate("P[p=%d].ID=%lld != TracerLinkedList[itracer=%d].ParentID=%lld\n", p, (long long) P[p].ID, itracer, (long long) TracerLinkedList[itracer].ParentID);
#endif

  P[p].NumberOfTracers--;
}

/*! \brief Add tracer to parent. Your responsibility to remove it elsewhere.
 *
 *  \param p P_index of destination parent
 *  \param itracer TLL_index of tracer to add
 */
void add_tracer_to_parent(int p, int itracer)
{
  int next = P[p].TracerHead;

  P[p].TracerHead = itracer;
  TracerLinkedList[itracer].Next = next;
  TracerLinkedList[itracer].Prev = -p - 1;      /* Prev = -Pindex - 1 */

#ifdef TRACER_MC_CHECKS
  if(P[p].ID != TracerLinkedList[itracer].ParentID)
    terminate("P[p=%d].ID=%lld != TracerLinkedList[itracer=%d].ParentID=%lld\n", p, (long long) P[p].ID, itracer, (long long) TracerLinkedList[itracer].ParentID);
#endif

  if(next >= 0)
    TracerLinkedList[next].Prev = itracer;

  P[p].NumberOfTracers++;
}

/*! \brief move a tracer from one parent to another, both particles local (any particle type)
 *
 *  \param p_from P_index of priginal parent
 *  \param p_to P_index of destination parent
 *  \param itracer TLL_index of tracer to move
 */
void move_tracer_between_parents(int p_from, int p_to, int itracer)
{
  if(P[p_from].Type == 4)       /* moving from a star or wind particle */
    terminate("the origin should always be a gas cell here");

  /* remove itracer from p */
  remove_tracer_from_parent(p_from, itracer);

#ifdef TRACER_MC_CHECKS
  TracerLinkedList[itracer].ParentID = P[p_to].ID;
#endif

  /* add itracer to pother */
  add_tracer_to_parent(p_to, itracer);

#ifdef TRACER_MC_NUM_FLUID_QUANTITIES

#if (TRACER_MC_EXCHANGE_COUNTER)
  TracerLinkedList[itracer].fluid_quantities[TracerMCExchangeCounterIndex] += 1;
#endif
#ifdef GFM
#if (TRACER_MC_LAST_STAR_TIME)
  if(P[p_to].Type == 4)
    {
      if(StarP[P[p_to].AuxDataID].BirthTime < 0)        /* moving to a wind particle */
        TracerLinkedList[itracer].fluid_quantities[TracerMCLastStarTimeIndex] += 3 * All.TimeMax;
      else                      /* moving to a star */
        TracerLinkedList[itracer].fluid_quantities[TracerMCLastStarTimeIndex] = 2 * All.TimeMax;
    }
#endif
#if (TRACER_MC_WIND_COUNTER)
  if(P[p_to].Type == 4)
    if(StarP[P[p_to].AuxDataID].BirthTime < 0)  /* moving to a wind particle */
      TracerLinkedList[itracer].fluid_quantities[TracerMCWindCounterIndex] += 1;
#endif
#endif /* GFM */
#if (TRACER_MC_EXCHANGE_DISTANCE)
  if(P[p_from].Type == 0)
    TracerLinkedList[itracer].fluid_quantities[TracerMCExchangeDistanceIndex] += get_cell_radius(p_from);
  if(P[p_to].Type == 0)
    TracerLinkedList[itracer].fluid_quantities[TracerMCExchangeDistanceIndex] += get_cell_radius(p_to);
#endif
#if (TRACER_MC_EXCHANGE_DISTANCE_ERROR)
  if(P[p_from].Type == 0)
    {
      TracerLinkedList[itracer].fluid_quantities[TracerMCExchangeDistanceErrorIndex] +=
        get_cell_radius(p_from) * (sqrt(TracerLinkedList[itracer].fluid_quantities[TracerMCExchangeCounterIndex]) - sqrt(TracerLinkedList[itracer].fluid_quantities[TracerMCExchangeCounterIndex] - 1));
    }
  if(P[p_to].Type == 0)
    {
      TracerLinkedList[itracer].fluid_quantities[TracerMCExchangeDistanceErrorIndex] +=
        get_cell_radius(p_to) * (sqrt(TracerLinkedList[itracer].fluid_quantities[TracerMCExchangeCounterIndex]) - sqrt(TracerLinkedList[itracer].fluid_quantities[TracerMCExchangeCounterIndex] - 1));
    }
#endif

#endif
}

#ifdef TRACER_MC_CHECKS
/*! \brief Helper used in check_tracer_lists().
 *
 *  \param p_id P_index of cell/particle to check
 */
int check_tracer_children(int p_id)
{
  int local_count = 0;
  int next = P[p_id].TracerHead;

  if(next < 0)                  /* no tracer children in this parent */
    return 0;

  while(next >= 0)              /* traverse forwards to end */
    {
      if(P[p_id].ID != TracerLinkedList[next].ParentID)
        terminate("p_id=%d P[p_id].Type=%d  P[p_id].NumberOfTracers=%d   P[p_id].ID=%lld  next=%d  TracerLinkedList[next].ParentID=%lld   parent ID mismatch", p_id, P[p_id].Type,
                  P[p_id].NumberOfTracers, (long long) P[p_id].ID, next, (long long) TracerLinkedList[next].ParentID);

      next = TracerLinkedList[next].Next;
      local_count++;
    }

  return local_count;
}

/*! \brief Verify integrity of tracer linked list.
 */
void check_tracer_lists(void)
{
  int i;
  int tot_tr_count = 0;
  int max_N_tracer, N_tracerload = N_tracer;
  long long all_tracers;

  MPI_Allreduce(&N_tracerload, &max_N_tracer, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  sumup_large_ints(1, &N_tracerload, &all_tracers);
  mpi_printf("TRACER_MC: max_N_tracer=%d total-tracers=%lld All.MaxPartTracer=%d\n", max_N_tracer, all_tracers, All.MaxPartTracer);

  /* check: global Monte Carlo tracer count must be conserved */
  if(All.N_alltracer_global != get_total_number_of_tracers(-1))
    terminate("TRACER_MC: N_alltracer_global is not conserved! N_alltracer_global=%lld but current number of tracers is %lld", All.N_alltracer_global, get_total_number_of_tracers(-1));

  /* check: all particle types with tracers have good Prev and Next connectivity */
  for(i = 0; i < NumPart; i++)
    tot_tr_count += check_tracer_children(i);

  /* check: total number of tracers in linked list equals the total local number */
  if(tot_tr_count != N_tracer)
    terminate("TRACER_MC: check_tracer_list: tot_tr_count %d != N_tracer %d", tot_tr_count, N_tracer);
}
#endif

#endif /* TRACER_MC */
