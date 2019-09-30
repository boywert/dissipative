/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/fm_star_formation/radpressthin.c
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


#include "../allvars.h"
#include "../proto.h"


#ifdef RADPRESS_OPT_THIN

int radpressthin_treeevaluate(int i, int mode, int *nexport, int *nsend_local);

extern int *TargetList;
extern int Nforces;

/*structures for radpressthin tensor*/
struct radpressthindata_in
{
  int NodeList[NODELISTLENGTH];
  MyDouble Pos[3];
  unsigned char SofteningType;
}
 *RadpressthinDataIn, *RadpressthinDataGet;


struct radpressthindata_out
{
  MyFloat RadPress[3];
}
 *RadpressthinDataResult, *RadpressthinDataOut;

struct radpressthinresultsimported_data
{
  int index;
  MyFloat RadPress[3];
}
 *RadpressthinResultsImported;

double PI = 3.141592653589;
/* Walk tree to compute contribution of radiation to single gas cell
 * and arrangement of particle communication */
void radpressthin(void)
{
  int idx, i, j, k, ngrp, dummy;
  int recvTask, nexport, nimport, place, ndone, ndone_flag;
  int target, ncount, n;
  double dt, dtime, dvel, lum_per_mass;



  /* allocate buffers to arrange communication */

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(data_index) + sizeof(struct data_nodelist) +
                                             sizeof(struct radpressthindata_in) +
                                             sizeof(struct radpressthindata_out) + sizemax(sizeof(struct radpressthindata_in), sizeof(struct radpressthindata_out))));
  DataIndexTable = (data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(data_index));
  DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));



  /* Create list of targets. We do this here to simplify the treatment of the two possible sources of points */

  TargetList = mymalloc("TargetList", (NumPart + Tree_NumPartImported) * sizeof(int));
  Tree_ResultIndexList = mymalloc("Tree_ResultIndexList", Tree_NumPartImported * sizeof(int));

  Nforces = 0;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Type == 0 && Tree_Task_list[i] == ThisTask)
        TargetList[Nforces++] = i;
    }

  for(i = 0, ncount = 0; i < Tree_NumPartImported; i++)
#ifndef HIERARCHICAL_GRAVITY
    if(Tree_Points[i].ActiveFlag)
#endif
    if(Tree_Points[i].Type == 0)
      {
        Tree_ResultIndexList[i] = ncount++;
        TargetList[Nforces++] = i + Tree_ImportedNodeOffset;
      }

  RadpressthinResultsImported = mymalloc("RadpressthinResultsImported", ncount * sizeof(struct radpressthinresultsimported_data));

  /* set radiation pressure to zero, will be incremented below */
  for(i = 0; i < ncount; i++)
    for(k = 0; k < 3; k++)
      RadpressthinResultsImported[i].RadPress[k] = 0.0;

  for(i = 0; i < NumGas; i++)
    for(k = 0; k < 3; k++)
      SphP[i].RadPress[k] = 0;


  i = 0;                        /* beginn with this index */

  do
    {
      for(j = 0; j < NTask; j++)
        {
          Send_count[j] = 0;
          Exportflag[j] = -1;
        }


      /* do local particles and prepare export list */
      for(nexport = 0; i < Nforces; i++)
        {
          if(radpressthin_treeevaluate(i, 0, &nexport, Send_count) < 0)
            break;
        }

      mysort_dataindex(DataIndexTable, nexport, sizeof(data_index), data_index_compare);

      MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

      for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
        {
          nimport += Recv_count[j];

          if(j > 0)
            {
              Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
              Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
            }
        }

      RadpressthinDataGet = (struct radpressthindata_in *) mymalloc("RadpressthinDataGet", nimport * sizeof(struct radpressthindata_in));
      RadpressthinDataIn = (struct radpressthindata_in *) mymalloc("RadpressthinDataIn", nexport * sizeof(struct radpressthindata_in));

      /* prepare particle data for export */

      for(j = 0; j < nexport; j++)
        {
          place = DataIndexTable[j].Index;
          target = TargetList[place];

          if(target < NumPart)
            {
              for(k = 0; k < 3; k++)
                RadpressthinDataIn[j].Pos[k] = P[target].Pos[k];
              RadpressthinDataIn[j].SofteningType = P[target].SofteningType;
            }
          else
            {
              target -= Tree_ImportedNodeOffset;

              for(k = 0; k < 3; k++)
                RadpressthinDataIn[j].Pos[k] = Tree_Points[target].Pos[k];
              RadpressthinDataIn[j].SofteningType = Tree_Points[target].SofteningType;
            }

          memcpy(RadpressthinDataIn[j].NodeList, DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
        }



      /* exchange particle data */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
          recvTask = ThisTask ^ ngrp;
          if(recvTask < NTask)
            if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
              /* get the particles */
              MPI_Sendrecv(&RadpressthinDataIn[Send_offset[recvTask]],
                           Send_count[recvTask] * sizeof(struct radpressthindata_in), MPI_BYTE,
                           recvTask, TAG_HYDRO_A,
                           &RadpressthinDataGet[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct radpressthindata_in), MPI_BYTE, recvTask, TAG_HYDRO_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

      myfree(RadpressthinDataIn);
      RadpressthinDataResult = (struct radpressthindata_out *) mymalloc("RadpressthinDataResult", nimport * sizeof(struct radpressthindata_out));
      RadpressthinDataOut = (struct radpressthindata_out *) mymalloc("RadpressthinDataOut", nexport * sizeof(struct radpressthindata_out));


      /* now do the particles that were sent to us */
      for(j = 0; j < nimport; j++)
        radpressthin_treeevaluate(j, 1, &dummy, &dummy);

      /* check whether this is the last iteration */
      if(i >= Nforces)
        ndone_flag = 1;
      else
        ndone_flag = 0;

      MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      /* get the result */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
          recvTask = ThisTask ^ ngrp;
          if(recvTask < NTask)
            if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
              MPI_Sendrecv(&RadpressthinDataResult[Recv_offset[recvTask]],
                           Recv_count[recvTask] * sizeof(struct radpressthindata_out),
                           MPI_BYTE, recvTask, TAG_HYDRO_B,
                           &RadpressthinDataOut[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct radpressthindata_out), MPI_BYTE, recvTask, TAG_HYDRO_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
      /* add the result to the local particles */
      for(j = 0; j < nexport; j++)
        {
          place = DataIndexTable[j].Index;
          target = TargetList[place];

          if(target < NumPart)
            for(k = 0; k < 3; k++)
              SphP[target].RadPress[k] += RadpressthinDataOut[j].RadPress[k];
          else
            {
              int idx = Tree_ResultIndexList[target - Tree_ImportedNodeOffset];

              for(k = 0; k < 3; k++)
                RadpressthinResultsImported[idx].RadPress[k] += RadpressthinDataOut[j].RadPress[k];
            }
        }

      myfree(RadpressthinDataOut);
      myfree(RadpressthinDataResult);
      myfree(RadpressthinDataGet);
    }
  while(ndone < NTask);


  /* now communicate the results in RadpressthinResultsImported */

  for(j = 0; j < NTask; j++)
    Recv_count[j] = 0;

  for(i = 0, n = 0, k = 0; i < NTask; i++)
    for(j = 0; j < Mesh_Recv_count[i]; j++, n++)
      {
#ifndef HIERARCHICAL_GRAVITY
        if(Tree_Points[n].ActiveFlag)
#endif
        if(Tree_Points[n].Type == 0)
          {
            RadpressthinResultsImported[k].index = Tree_Points[n].index;
            Recv_count[i]++;
            k++;
          }
      }

  MPI_Alltoall(Recv_count, 1, MPI_INT, Send_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nexport = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nexport += Send_count[j];
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  struct radpressthinresultsimported_data *tmp_results = mymalloc("tmp_results", nexport * sizeof(struct radpressthinresultsimported_data));
  memset(tmp_results, -1, nexport * sizeof(struct radpressthinresultsimported_data));

  /* exchange  data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          MPI_Sendrecv(&RadpressthinResultsImported[Recv_offset[recvTask]],
                       Recv_count[recvTask] * sizeof(struct radpressthinresultsimported_data), MPI_BYTE, recvTask, TAG_FOF_A,
                       &tmp_results[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct radpressthinresultsimported_data), MPI_BYTE, recvTask, TAG_FOF_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

  for(i = 0; i < nexport; i++)
    {
      int target = tmp_results[i].index;
      for(k = 0; k < 3; k++)
        SphP[target].RadPress[k] = tmp_results[i].RadPress[k];
    }

  myfree(tmp_results);

  myfree(RadpressthinResultsImported);
  myfree(Tree_ResultIndexList);
  myfree(TargetList);

  myfree(DataNodeList);
  myfree(DataIndexTable);

  double local_ekin_min = MAX_REAL_NUMBER, local_ekin_max = MIN_REAL_NUMBER;
  double global_ekin_min, global_ekin_max;
  double ekin_before, ekin_after, ekin_ratio;
  double nH, radius_cell, nHI;

  /* [HIcross_sec] = UnitLength^2 */
  double HIcross_sec = 6.3e-18 * (1.0 / All.UnitLength_in_cm / All.UnitLength_in_cm * All.HubbleParam * All.HubbleParam);

  /* [c_light] = UnitVelocity */
  double c_light = CLIGHT / All.UnitVelocity_in_cm_per_s;


/* [All.IonizingLumPerSolarMass] = photons s^-1 Msun^-1; */
#ifdef RADPRESS_OPT_THIN_LUMPERMASS
  lum_per_mass = (All.UnitMass_in_g / SOLAR_MASS / All.HubbleParam * All.UnitTime_in_s / All.HubbleParam) * All.IonizingLumPerSolarMass;
#else
  lum_per_mass = (All.UnitTime_in_s / All.HubbleParam) * All.IonizingLumPerSolarMass;
#endif

  double energy_per_photon = 13.6 * ELECTRONVOLT_IN_ERGS * (All.HubbleParam / All.UnitEnergy_in_cgs);
  /* [All.IonizingLumPerSolarMass] = erg s^-1 Msun^-1; [lum_per_mass] = UnitEnergy UnitTime^-1 UnitMass^-1 */
  lum_per_mass *= energy_per_photon;


  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      /* dt = (All.otvet_Radiation_Ti_endstep - All.otvet_Radiation_Ti_begstep) * All.Timebase_interval; */
      dt = (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;
      /* in comoving case, timestep is dloga at this point. Convert to dt */
      if(All.ComovingIntegrationOn)
        dt /= All.cf_hubble_a;


      dtime = dt;

      /* 3D number of atoms density in code units */
      nH = HYDROGEN_MASSFRAC * SphP[i].Density / PROTONMASS * All.UnitMass_in_g / All.HubbleParam;      /* LVS shall I multiply by a3inv?? */
      nHI = SphP[i].HI * nH;

      /* subtract kin. energy */
      SphP[i].Energy -= 0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);
      ekin_before = 0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);

      /* do the kick for gas cells */
      for(k = 0; k < 3; k++)
        {

          /* [SphP[i].RadPress[k]] = UnitEnergy UnitLength^-2 */
          SphP[i].RadPress[k] *= lum_per_mass / 4. / PI / c_light;

          SphP[i].RadPress[k] *= nHI * HIcross_sec / SphP[i].Density;   /* in nHI I use results of OTVET */

          dvel = (SphP[i].RadPress[k] * dtime);

          /* add velocity */
          P[i].Vel[k] += dvel;

          /* add momentum */
          SphP[i].Momentum[k] += P[i].Mass * dvel;

          if(!gsl_finite(dvel))
            terminate("RADPRESS_OPT_THIN: bad radiation pressure dvel=%g P[i].Mass=%g SphP[i].RadPress[k]=%g dtime=%g ID=%d", dvel, P[i].Mass, SphP[i].RadPress[k], dtime, P[i].ID);
        }

      /* add kin. energy */
      SphP[i].Energy += 0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);
      ekin_after = 0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);

      ekin_ratio = (ekin_before > 0) ? (ekin_after / ekin_before) : 1.0;

      if(ekin_ratio < local_ekin_min)
        local_ekin_min = ekin_ratio;

      if(ekin_ratio > local_ekin_max)
        local_ekin_max = ekin_ratio;
    }

  MPI_Reduce(&local_ekin_min, &global_ekin_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_ekin_max, &global_ekin_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  mpi_printf("RADPRESS_OPT_THIN: kin. energy ratio = (%g - %g) \n", global_ekin_min, global_ekin_max);
}


/*! This routine computes the radiation pressure in optically thin approx.
 *  for a given local
 *  particle, or for a particle in the communication buffer. Depending on
 *  the value of TypeOfOpeningCriterion, either the geometrical BH
 *  cell-opening criterion, or the `relative' opening criterion is used.
 */
int radpressthin_treeevaluate(int i, int mode, int *nexport, int *nsend_local)
{
  struct NODE *nop = 0;
  int k, no, nodesinlist, ninteractions, task, ptype;
  int target;
  int listindex = 0;
  double r2, dx = 0, dy = 0, dz = 0, r;
  double pos_x, pos_y, pos_z, h_i, h_j, hmax;
  double young_stellar_mass, age_in_Gyr, r2_invfac;
  MyDouble RadPress[3] = { 0, 0, 0 };
#if defined(PERIODIC) && !defined(GRAVITY_NOT_PERIODIC)
  double xtmp, ytmp, ztmp;
#endif
  int nexport_save;
  nexport_save = *nexport;
  double young_stellar_mass_accum = 0;

  ninteractions = 0;
  nodesinlist = 0;
  ptype = 0;                    /* we only deal with gas particles */


  if(mode == 0)
    {
      target = TargetList[i];

      if(target < NumPart)
        {
          pos_x = Tree_Pos_list[3 * target + 0];
          pos_y = Tree_Pos_list[3 * target + 1];
          pos_z = Tree_Pos_list[3 * target + 2];
          h_i = All.ForceSoftening[P[target].SofteningType];
        }
      else
        {
          int n = target - Tree_ImportedNodeOffset;

          pos_x = Tree_Points[n].Pos[0];
          pos_y = Tree_Points[n].Pos[1];
          pos_z = Tree_Points[n].Pos[2];
          h_i = All.ForceSoftening[Tree_Points[n].SofteningType];
        }
    }
  else
    {
      target = i;
      pos_x = RadpressthinDataGet[target].Pos[0];
      pos_y = RadpressthinDataGet[target].Pos[1];
      pos_z = RadpressthinDataGet[target].Pos[2];
      h_i = All.ForceSoftening[RadpressthinDataGet[target].SofteningType];
    }

  if(mode == 0)
    no = Tree_MaxPart;          /* root node */
  else
    {
      nodesinlist++;
      no = RadpressthinDataGet[target].NodeList[0];
      no = Nodes[no].u.d.nextnode;      /* open it */
    }

  while(no >= 0)
    {
      while(no >= 0)
        {
          if(no < Tree_MaxPart) /* single particle */
            {
              /* the index of the node is the index of the particle */
              /* observe the sign */

              if(P[no].Type == 4)
                {
                  dx = GRAVITY_NEAREST_X(Tree_Pos_list[3 * no + 0] - pos_x);    /* x_star - x_gas ==> -dx is radial */
                  dy = GRAVITY_NEAREST_Y(Tree_Pos_list[3 * no + 1] - pos_y);
                  dz = GRAVITY_NEAREST_Z(Tree_Pos_list[3 * no + 2] - pos_z);
                  young_stellar_mass = 0.;

#ifdef GFM
                  age_in_Gyr = get_time_difference_in_Gyr(STP(no).BirthTime, All.Ti_Current);
#else
                  age_in_Gyr = 0.;
#endif

                  if(age_in_Gyr < All.RadPressure_MaxAge)
                    young_stellar_mass = P[no].Mass;
                }
              else
                {
                  dx = 0;
                  dy = 0;
                  dz = 0;
                  young_stellar_mass = 0;
                }

              r2 = dx * dx + dy * dy + dz * dz;
              h_j = All.ForceSoftening[P[no].SofteningType];

              if(h_j > h_i)
                hmax = h_j;
              else
                hmax = h_i;

              r = sqrt(r2 + hmax * hmax);

              no = Nextnode[no];

            }
          else if(no < Tree_MaxPart + Tree_MaxNodes)    /* internal node */
            {
              if(mode == 1)
                {
                  if(no < Tree_FirstNonTopLevelNode)    /* we reached a top-level node again, which means that we are done with the branch */
                    {
                      no = -1;
                      continue;
                    }
                }

              nop = &Nodes[no];
              young_stellar_mass = nop->young_stellar_mass;

              dx = GRAVITY_NEAREST_X(nop->young_stellar_s[0] - pos_x);
              dy = GRAVITY_NEAREST_Y(nop->young_stellar_s[1] - pos_y);
              dz = GRAVITY_NEAREST_Z(nop->young_stellar_s[2] - pos_z);

              r2 = dx * dx + dy * dy + dz * dz;

              /* we have an  internal node. Need to check opening criterion */

              if(All.ErrTolTheta)       /* check Barnes-Hut opening criterion */
                {
                  if(nop->len * nop->len > r2 * All.ErrTolTheta * All.ErrTolTheta)
                    {
                      /* open cell */
                      no = nop->u.d.nextnode;
                      continue;
                    }
                }
              else
                {
                  /* check in addition whether we lie inside the cell */
                  if(fabs(nop->center[0] - pos_x) < 0.60 * nop->len)
                    {
                      if(fabs(nop->center[1] - pos_y) < 0.60 * nop->len)
                        {
                          if(fabs(nop->center[2] - pos_z) < 0.60 * nop->len)
                            {
                              no = nop->u.d.nextnode;
                              continue;
                            }
                        }
                    }
                }

              h_j = nop->u.d.maxsoft;

              if(h_j > h_i)
                {
                  if(r2 < h_j * h_j)
                    {
                      /* open cell */
                      no = nop->u.d.nextnode;
                      continue;
                    }
                  hmax = h_j;
                }
              else
                hmax = h_i;

              r = sqrt(r2 + hmax * hmax);

              no = nop->u.d.sibling;
            }

          else if(no >= Tree_ImportedNodeOffset)        /* point from imported nodelist */
            {
              int n = no - Tree_ImportedNodeOffset;


              if(Tree_Points[n].Type == 4)
                {
                  dx = GRAVITY_NEAREST_X(Tree_Points[n].Pos[0] - pos_x);
                  dy = GRAVITY_NEAREST_Y(Tree_Points[n].Pos[1] - pos_y);
                  dz = GRAVITY_NEAREST_Z(Tree_Points[n].Pos[2] - pos_z);
                  young_stellar_mass = 0.;

#ifdef GFM
                  age_in_Gyr = get_time_difference_in_Gyr(Tree_Points[n].BirthTime, All.Time);
#else
                  age_in_Gyr = 0.;
#endif
                  if(age_in_Gyr < All.RadPressure_MaxAge)
                    young_stellar_mass = Tree_Points[n].Mass;

                }
              else
                {
                  dx = 0;
                  dy = 0;
                  dz = 0;
                  young_stellar_mass = 0;
                }


              r2 = dx * dx + dy * dy + dz * dz;
              h_j = All.ForceSoftening[Tree_Points[n].SofteningType];

              if(h_j > h_i)
                hmax = h_j;
              else
                hmax = h_i;

              r = sqrt(r2 + hmax * hmax);

              no = Nextnode[no - Tree_MaxNodes];
            }
          else                  /* pseudo particle */
            {
              if(mode == 0)
                {
                  if(Exportflag[task = DomainNewTask[no - (Tree_MaxPart + Tree_MaxNodes)]] != i)
                    {
                      Exportflag[task] = i;
                      Exportnodecount[task] = NODELISTLENGTH;
                    }

                  if(Exportnodecount[task] == NODELISTLENGTH)
                    {
                      if(*nexport >= All.BunchSize)
                        {
                          /* out of buffer space. Need to discard work for this particle and interrupt */
                          *nexport = nexport_save;
                          if(nexport_save == 0)
                            terminate("17998"); /* in this case, the buffer is too small to process even a single particle */
                          for(task = 0; task < NTask; task++)
                            nsend_local[task] = 0;
                          for(no = 0; no < nexport_save; no++)
                            nsend_local[DataIndexTable[no].Task]++;
                          return -1;
                        }
                      Exportnodecount[task] = 0;
                      Exportindex[task] = *nexport;
                      DataIndexTable[*nexport].Task = task;
                      DataIndexTable[*nexport].Index = i;
                      DataIndexTable[*nexport].IndexGet = *nexport;
                      *nexport = *nexport + 1;
                      nsend_local[task]++;
                    }

                  DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]++] = DomainNodeIndex[no - (Tree_MaxPart + Tree_MaxNodes)];

                  if(Exportnodecount[task] < NODELISTLENGTH)
                    DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]] = -1;
                }

              no = Nextnode[no - Tree_MaxNodes];
              continue;
            }


          if(young_stellar_mass > 0)
            ninteractions++;


          r2_invfac = 1.0 / (r2 + hmax * hmax);

          RadPress[0] += young_stellar_mass * r2_invfac * (-dx) / r;    /* -dx is radially outwards from star */
          RadPress[1] += young_stellar_mass * r2_invfac * (-dy) / r;
          RadPress[2] += young_stellar_mass * r2_invfac * (-dz) / r;

          young_stellar_mass_accum += young_stellar_mass;
        }
      if(mode == 1)
        {
          listindex++;
          if(listindex < NODELISTLENGTH)
            {
              no = RadpressthinDataGet[target].NodeList[listindex];
              if(no >= 0)
                {
                  nodesinlist++;
                  no = Nodes[no].u.d.nextnode;  /* open it */
                }
            }
        }
    }

/* IF no weighted by mass, all sources emmit at the rate in the parameter file */
#ifndef RADPRESS_OPT_THIN_LUMPERMASS
  if(young_stellar_mass_accum > 0)
    for(k = 0; k < 3; k++)
      RadPress[k] /= young_stellar_mass_accum;
#endif





  /* store result at the proper place */
  if(mode == 0)
    {
      if(target < NumPart)
        {
          for(k = 0; k < 3; k++)
            SphP[target].RadPress[k] = RadPress[k];
        }
      else
        {
          int idx = Tree_ResultIndexList[target - Tree_ImportedNodeOffset];
          for(k = 0; k < 3; k++)
            RadpressthinResultsImported[idx].RadPress[k] = RadPress[k];

        }
    }
  else
    {
      for(k = 0; k < 3; k++)
        RadpressthinDataResult[target].RadPress[k] = RadPress[k];

      *nexport = nodesinlist;
    }


  return ninteractions;
}

#endif
