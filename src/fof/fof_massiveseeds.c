/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/fof/fof_massiveseeds.c
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
#include <sys/stat.h>
#include <sys/types.h>
#include <gsl/gsl_math.h>
#include <inttypes.h>

#include "../allvars.h"
#include "../proto.h"

#ifdef FOF
#ifdef BLACK_HOLES
#ifdef MASSIVE_SEEDS

static struct densdata_in
{
  MyDouble Pos[3];
  MyFloat BH_Hsml;
  MyIDType ID;
  int NodeList[NODELISTLENGTH];
}
 *DensDataIn, *DensDataGet;

static struct densdata_out
{
  MyFloat Ngb;
  MyFloat Mass;
  MyFloat CoM[3];
}
 *DensDataResult, *DensDataOut;

/* find a predefined number of gas Ngbs of each BH seed, add their mass to the seed mass, and mark them for deletion */
void blackhole_massiveseeds(int nseed, int *seed_indices)
{
  MyFloat *Left, *Right;
  int i, j, k, ndone, ndone_flag, npleft, dummy, iter = 0;
  int ngrp, recvTask, place, nexport, nimport;
  long long ntot;
  double *mold;
  MyFloat *vxold, *vyold, *vzold;

  mold = (double *) mymalloc("mold", nseed * sizeof(double));
  vxold = (MyFloat *) mymalloc("vxold", nseed * sizeof(MyFloat));
  vyold = (MyFloat *) mymalloc("vyold", nseed * sizeof(MyFloat));
  vzold = (MyFloat *) mymalloc("vzold", nseed * sizeof(MyFloat));

  Left = (MyFloat *) mymalloc("Left", NumPart * sizeof(MyFloat));
  Right = (MyFloat *) mymalloc("Right", NumPart * sizeof(MyFloat));

  memset(mold, 0, nseed * sizeof(double));
  memset(vxold, 0, nseed * sizeof(MyFloat));
  memset(vyold, 0, nseed * sizeof(MyFloat));
  memset(vzold, 0, nseed * sizeof(MyFloat));

  for(i = 0; i < nseed; i++)
    Left[seed_indices[i]] = Right[seed_indices[i]] = 0;

  for(i = 0; i < nseed; i++)
    {
      mold[i] = P[seed_indices[i]].Mass;
      vxold[i] = P[seed_indices[i]].Vel[0];
      vyold[i] = P[seed_indices[i]].Vel[1];
      vzold[i] = P[seed_indices[i]].Vel[2];
    }

  Ngblist = (int *) mymalloc("Ngblist", NumPart * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(data_index) + sizeof(struct data_nodelist) +
                                             sizeof(struct densdata_in) + sizeof(struct densdata_out) + sizemax(sizeof(struct densdata_in), sizeof(struct densdata_out))));
  DataIndexTable = (data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(data_index));
  DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

  do
    {
      i = 0;                    /* begin with this index */

      for(k = 0; k < NumGas; k++)
        if(P[k].ID != 0)
          SphP[k].SwallowID = 0;

      do
        {
          for(j = 0; j < NTask; j++)
            {
              Send_count[j] = 0;
              Exportflag[j] = -1;
            }

          /* do local particles and prepare export list */

          for(nexport = 0; i < nseed; i++)
            {
              int target = seed_indices[i];
              if(blackhole_isactive(target))
                {
                  if(blackhole_gasngb_evaluate(target, 0, &nexport, Send_count) < 0)
                    break;
                }
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

          DensDataGet = (struct densdata_in *) mymalloc("DensDataGet", nimport * sizeof(struct densdata_in));
          DensDataIn = (struct densdata_in *) mymalloc("DensDataIn", nexport * sizeof(struct densdata_in));

          /* prepare particle data for export */
          for(j = 0; j < nexport; j++)
            {
              place = DataIndexTable[j].Index;

              DensDataIn[j].Pos[0] = P[place].Pos[0];
              DensDataIn[j].Pos[1] = P[place].Pos[1];
              DensDataIn[j].Pos[2] = P[place].Pos[2];
              DensDataIn[j].BH_Hsml = BPP(place).BH_Hsml;
              DensDataIn[j].ID = P[place].ID;
              memcpy(DensDataIn[j].NodeList, DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
            }

          /* exchange particle data */
          for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
            {
              recvTask = ThisTask ^ ngrp;
              if(recvTask < NTask)
                if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                  /* get the particles */
                  MPI_Sendrecv(&DensDataIn[Send_offset[recvTask]],
                               Send_count[recvTask] * sizeof(struct densdata_in), MPI_BYTE,
                               recvTask, TAG_DENS_A,
                               &DensDataGet[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct densdata_in), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

          myfree(DensDataIn);
          DensDataResult = (struct densdata_out *) mymalloc("DensDataResult", nimport * sizeof(struct densdata_out));
          DensDataOut = (struct densdata_out *) mymalloc("DensDataOut", nexport * sizeof(struct densdata_out));

          /* now do the particles that were sent to us */
          for(j = 0; j < nimport; j++)
            blackhole_gasngb_evaluate(j, 1, &dummy, &dummy);

          if(i >= nseed)
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
                  MPI_Sendrecv(&DensDataResult[Recv_offset[recvTask]],
                               Recv_count[recvTask] * sizeof(struct densdata_out),
                               MPI_BYTE, recvTask, TAG_DENS_B,
                               &DensDataOut[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct densdata_out), MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

          /* add the result to the local particles */
          for(j = 0; j < nexport; j++)
            {
              place = DataIndexTable[j].Index;
              BPP(place).BH_NumNgb += DensDataOut[j].Ngb;
              P[place].Mass += DensDataOut[j].Mass;
              P[place].Vel[0] += DensDataOut[j].CoM[0];
              P[place].Vel[1] += DensDataOut[j].CoM[1];
              P[place].Vel[2] += DensDataOut[j].CoM[2];
            }

          myfree(DensDataOut);
          myfree(DensDataResult);
          myfree(DensDataGet);
        }
      while(ndone < NTask);

      for(j = 0; j < nseed; j++)
        {
          i = seed_indices[j];

          if(blackhole_isactive(i))
            BPP(i).BH_NumNgb /= All.ReferenceGasPartMass;
        }

      for(j = 0, npleft = 0; j < nseed; j++)
        {
          i = seed_indices[j];

          if(blackhole_isactive(i))
            {
              /* now check whether we had enough neighbours */

              if(BPP(i).BH_NumNgb < (All.DesNumNgbSeed - All.MaxNumNgbDeviation) || (BPP(i).BH_NumNgb > (All.DesNumNgbSeed + All.MaxNumNgbDeviation) && BPP(i).BH_Hsml > (1.01 * All.MinGasHsml)))
                {
                  /* need to redo this particle */
                  npleft++;

                  if(Left[i] > 0 && Right[i] > 0)
                    if((Right[i] - Left[i]) < 1.0e-3 * Left[i])
                      {
                        /* this one should be ok */
                        npleft--;
                        P[i].TimeBinGrav = -P[i].TimeBinGrav - 1;       /* Mark as inactive */
                        P[i].Mass += mold[j];
                        P[i].Vel[0] += mold[j] * vxold[j];
                        P[i].Vel[1] += mold[j] * vyold[j];
                        P[i].Vel[2] += mold[j] * vzold[j];
                        P[i].Vel[0] /= P[i].Mass;
                        P[i].Vel[1] /= P[i].Mass;
                        P[i].Vel[2] /= P[i].Mass;

                        continue;
                      }

                  if(BPP(i).BH_NumNgb < (All.DesNumNgbSeed - All.MaxNumNgbDeviation))
                    Left[i] = dmax(BPP(i).BH_Hsml, Left[i]);
                  else
                    {
                      if(Right[i] != 0)
                        {
                          if(BPP(i).BH_Hsml < Right[i])
                            Right[i] = BPP(i).BH_Hsml;
                        }
                      else
                        Right[i] = BPP(i).BH_Hsml;
                    }

                  if(iter >= MAXITER - 10)
                    {
                      printf
                        ("i=%d task=%d ID=%d Hsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g\n   pos=(%g|%g|%g)\n",
                         i, ThisTask, (int) P[i].ID, BPP(i).BH_Hsml, Left[i], Right[i], (double) BPP(i).BH_NumNgb, Right[i] - Left[i], P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
                      myflush(stdout);
                    }

                  if(Right[i] > 0 && Left[i] > 0)
                    BPP(i).BH_Hsml = pow(0.5 * (pow(Left[i], 3) + pow(Right[i], 3)), 1.0 / 3);
                  else
                    {
                      if(Right[i] == 0 && Left[i] == 0)
                        terminate("Right[i] = Left[i] = 0 in the determination of BH_Hsml");    /* can't occur */

                      if(Right[i] == 0 && Left[i] > 0)
                        BPP(i).BH_Hsml *= 1.26;

                      if(Right[i] > 0 && Left[i] == 0)
                        BPP(i).BH_Hsml /= 1.26;
                    }

                  if(BPP(i).BH_Hsml < All.MinGasHsml)
                    BPP(i).BH_Hsml = All.MinGasHsml;

                  if(Left[i] > All.SeedMaxAccretionRadius)
                    {
                      /* this will stop the search for a new BH smoothing length in the next iteration */
                      BPP(i).BH_Hsml = Left[i] = Right[i] = All.SeedMaxAccretionRadius;
                    }
                }
              else
                {
                  P[i].TimeBinGrav = -P[i].TimeBinGrav - 1;     /* Mark as inactive */
                  P[i].Mass += mold[j];
                  P[i].Vel[0] += mold[j] * vxold[j];
                  P[i].Vel[1] += mold[j] * vyold[j];
                  P[i].Vel[2] += mold[j] * vzold[j];
                  P[i].Vel[0] /= P[i].Mass;
                  P[i].Vel[1] /= P[i].Mass;
                  P[i].Vel[2] /= P[i].Mass;

                }
            }
        }

      sumup_large_ints(1, &npleft, &ntot);

      if(ntot > 0)
        {
          iter++;

          if(iter > 0)
            mpi_printf("BLACK_HOLE SEEDS: blackhole ngb iteration %d: need to repeat for %lld particles.\n", iter, ntot);

          if(iter > MAXITER)
            terminate("failed to converge in neighbour iteration in fof_massiveseeds()\n");
        }

      remove_swallowed_gas(nseed, seed_indices);
    }
  while(ntot > 0);

  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);
  myfree(Right);
  myfree(Left);
  myfree(vzold);
  myfree(vyold);
  myfree(vxold);
  myfree(mold);


  /* mark as active again */
  for(j = 0; j < nseed; j++)
    {
      i = seed_indices[j];
      if(P[i].TimeBinGrav < 0)
        P[i].TimeBinGrav = -P[i].TimeBinGrav - 1;

#ifdef INDIVIDUAL_GRAVITY_SOFTENING
      if(((1 << P[i].Type) & (INDIVIDUAL_GRAVITY_SOFTENING)))
        P[i].SofteningType = get_softening_type_from_mass(P[i].Mass);
#endif
    }

  for(j = 0; j < NumGas; j++)
    if(P[j].ID != 0)
      SphP[j].SwallowID = 0;

}

/*! This function represents the core of the blackhole density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
int blackhole_gasngb_evaluate(int target, int mode, int *nexport, int *nsend_local)
{
  int j, n, k;
  int startnode, numngb, numngb_inbox, listindex = 0;
  double h, h2, hinv, hinv3, hinv4;
  double wk;
  double dx, dy, dz, r, r2, u;
  MyFloat weighted_numngb;
  MyDouble *pos;
  MyFloat accreted_mass = 0, center_of_mass[3];
  MyIDType id;

  weighted_numngb = 0;

  for(k = 0; k < 3; k++)
    center_of_mass[k] = 0;

  if(mode == 0)
    {
      pos = P[target].Pos;
      h = BPP(target).BH_Hsml;
      id = P[target].ID;
    }
  else
    {
      pos = DensDataGet[target].Pos;
      h = DensDataGet[target].BH_Hsml;
      id = DensDataGet[target].ID;
    }

  h2 = h * h;
  hinv = 1.0 / h;
#ifndef  TWODIMS
  hinv3 = hinv * hinv * hinv;
#else
  hinv3 = hinv * hinv / boxSize_Z;
#endif
  hinv4 = hinv3 * hinv;

  if(mode == 0)
    {
      startnode = Ngb_MaxPart;  /* root node */
    }
  else
    {
      startnode = DensDataGet[target].NodeList[0];
      startnode = Ngb_Nodes[startnode].u.d.nextnode;    /* open it */
    }

  numngb = 0;

  while(startnode >= 0)
    {
      while(startnode >= 0)
        {
          numngb_inbox = ngb_treefind_variable(pos, h, target, &startnode, mode, nexport, nsend_local);

          if(numngb_inbox < 0)
            return -1;

          for(n = 0; n < numngb_inbox; n++)
            {
              j = Ngblist[n];

              if(P[j].Type == 0 && P[j].ID != 0)
                {
                  dx = pos[0] - P[j].Pos[0];
                  dy = pos[1] - P[j].Pos[1];
                  dz = pos[2] - P[j].Pos[2];

#ifdef PERIODIC                 /*  now find the closest image in the given box size  */
                  if(dx > boxHalf_X)
                    dx -= boxSize_X;
                  if(dx < -boxHalf_X)
                    dx += boxSize_X;
                  if(dy > boxHalf_Y)
                    dy -= boxSize_Y;
                  if(dy < -boxHalf_Y)
                    dy += boxSize_Y;
                  if(dz > boxHalf_Z)
                    dz -= boxSize_Z;
                  if(dz < -boxHalf_Z)
                    dz += boxSize_Z;
#endif
                  r2 = dx * dx + dy * dy + dz * dz;

                  if(r2 < h2)
                    {
                      numngb++;

                      r = sqrt(r2);

                      u = r * hinv;

                      if(u < 0.5)
                        wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
                      else
                        wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

                      weighted_numngb += (NORM_COEFF * P[j].Mass * wk / hinv3); /* 4.0/3 * PI = 4.188790204786 */

                      accreted_mass += P[j].Mass;

                      for(k = 0; k < 3; k++)
                        center_of_mass[k] += P[j].Mass * P[j].Vel[k];

                      if(SphP[j].SwallowID != 0)
                        terminate("SphP.SwallowID = %llu BH_Hsml = %g\n", (long long) SphP[j].SwallowID, h);


                      SphP[j].SwallowID = id;

                    }
                }
            }
        }

      if(mode == 1)
        {
          listindex++;
          if(listindex < NODELISTLENGTH)
            {
              startnode = DensDataGet[target].NodeList[listindex];
              if(startnode >= 0)
                startnode = Ngb_Nodes[startnode].u.d.nextnode;  /* open it */
            }
        }
    }

  if(mode == 0)
    {
      P[target].Mass = accreted_mass;
      BPP(target).BH_NumNgb = weighted_numngb;
      for(k = 0; k < 3; k++)
        P[target].Vel[k] = center_of_mass[k];
    }
  else
    {
      DensDataResult[target].Mass = accreted_mass;
      DensDataResult[target].Ngb = weighted_numngb;
      for(k = 0; k < 3; k++)
        DensDataResult[target].CoM[k] = center_of_mass[k];
    }

  return 0;
}

void remove_swallowed_gas(int nseed, int *seed_indices)
{
  int i, j, k, nseed_remove = 0;
  int *seed_indices_remove;
  int total_nseed_remove = 0;

  for(j = 0; j < nseed; j++)
    {
      i = seed_indices[j];
      if(P[i].TimeBinGrav < 0)
        {
          nseed_remove++;
        }
    }

  seed_indices_remove = mymalloc("seed_indices_remove", nseed_remove * sizeof(int));

  for(j = 0, k = 0; j < nseed; j++)
    {
      i = seed_indices[j];
      if(P[i].TimeBinGrav < 0)
        {
          seed_indices_remove[k] = i;
          k++;
        }
    }

  MPI_Allreduce(&nseed_remove, &total_nseed_remove, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(total_nseed_remove > 0)
    {
      int *common_nseed, *disp;
      MyIDType *seed_ID, *tot_seed_ID;

      seed_ID = mymalloc("seed_ID", nseed_remove * sizeof(MyIDType));
      tot_seed_ID = mymalloc("tot_seed_ID", total_nseed_remove * sizeof(MyIDType));

      for(j = 0; j < nseed_remove; j++)
        seed_ID[j] = 0;

      for(j = 0; j < total_nseed_remove; j++)
        tot_seed_ID[j] = 0;

      for(j = 0; j < nseed_remove; j++)
        {
          i = seed_indices_remove[j];
          seed_ID[j] = P[i].ID;
        }

      common_nseed = mymalloc("common_nseed", NTask * sizeof(int));
      disp = mymalloc("disp", NTask * sizeof(int));

      MPI_Allgather(&nseed_remove, 1, MPI_INT, common_nseed, 1, MPI_INT, MPI_COMM_WORLD);

      for(k = 1, disp[0] = 0; k < NTask; k++)
        disp[k] = disp[k - 1] + common_nseed[k - 1];

#ifndef LONGIDS
      MPI_Allgatherv(seed_ID, nseed_remove, MPI_UNSIGNED, tot_seed_ID, common_nseed, disp, MPI_UNSIGNED, MPI_COMM_WORLD);
#else
      MPI_Allgatherv(seed_ID, nseed_remove, MPI_UNSIGNED_LONG_LONG, tot_seed_ID, common_nseed, disp, MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);
#endif

      for(j = 0; j < total_nseed_remove; j++)
        {
          for(k = 0; k < NumGas; k++)
            {
              if(SphP[k].SwallowID == tot_seed_ID[j])
                {
                  P[k].Mass = 0;
                  P[k].ID = 0;
                  SphP[k].SwallowID = 0;
#ifdef VORONOI_DYNAMIC_UPDATE
                  voronoi_remove_connection(k);
#endif
                }
            }
        }

      myfree(disp);
      myfree(common_nseed);
      myfree(tot_seed_ID);
      myfree(seed_ID);
    }
  myfree(seed_indices_remove);
}

#endif
#endif
#endif
