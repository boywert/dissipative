/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/amr/amr_refinement.c
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

#include "../allvars.h"
#include "../proto.h"

#ifdef AMR_REMAP

static int distribute_evaluate(int target, int mode, int *nexport, int *nsend_local);
static void distribute_sph_particles(void);

static int spread_evaluate(int target, int mode, int *nexport, int *nsend_local);
static void spread_sph_particles(void);

static MyFloat *NumNgb, *DhsmlDensityFactor;
MyIDType amr_IDNew = 0;

int amr_remap(void)
{
  int i, j, k, no, numnodes;
  long long ngas_count_all_old;
  double vol, voltot, mgas, mtot;
  int flag_all, flag = 0;

  mpi_printf("\n\nAMR: Remaping particles to AMR mesh\n\n");

  ngas_count_all_old = All.TotNumGas;

  ngb_treefree();
  amr_amrtree = 1;

  domain_free();
  domain_Decomposition();       /* do new domain decomposition, will also make a new chained-list of synchronized particles */

  //build a minimal ngb tree
  ngb_treeallocate();
  ngb_treebuild(0);

  int *topNodeLevel = mymalloc_movable(&topNodeLevel, "topNodeLevel", sizeof(int) * NTopleaves);

  for(i = 0; i < NTopleaves; i++)
    {
      if(DomainTask[i] == ThisTask)
        {
          int index = Ngb_DomainNodeIndex[i];
          topNodeLevel[i] = imax(Ngb_Nodes[index].level, All.MinRefLevel);
        }
    }

  //FIXME ensure AMR mesh validity between top node patches

  int count_leaves = 0, count_leaves_all;
  for(i = 0; i < NTopleaves; i++)
    {
      if(DomainTask[i] == ThisTask)
        {
          int index = Ngb_DomainNodeIndex[i];
          count_leaves += 1 << (3 * (topNodeLevel[i] + 1 - Ngb_Nodes[index].level));
        }
    }

  MPI_Allreduce(&count_leaves, &count_leaves_all, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  mpi_printf("AMR_REMAP: starting with base grid of %d cells\n\n", count_leaves_all);

  /* determine maximum ID */
  MyIDType maxid, newid, *tmp;
  int *list;

  for(i = 0, maxid = 0; i < NumPart; i++)
    if(P[i].ID > maxid)
      maxid = P[i].ID;

  tmp = mymalloc("tmp", NTask * sizeof(MyIDType));

  MPI_Allgather(&maxid, sizeof(MyIDType), MPI_BYTE, tmp, sizeof(MyIDType), MPI_BYTE, MPI_COMM_WORLD);

  for(i = 0; i < NTask; i++)
    if(tmp[i] > maxid)
      maxid = tmp[i];

  myfree(tmp);
  // maxid is now the total maximum ID number of all particles

  list = mymalloc("list", NTask * sizeof(int));

  MPI_Allgather(&count_leaves, 1, MPI_INT, list, 1, MPI_INT, MPI_COMM_WORLD);

  newid = maxid + 1;
  for(i = 0; i < ThisTask; i++)
    newid += list[i];

  myfree(list);

  // newid is now the maxid+total of count_leaves over all previous tasks
  amr_IDNew = maxid + 1;        /* old gas particles will have IDs below this */

  if((NumGas + count_leaves >= All.MaxPartSph) || (NumPart + count_leaves >= All.MaxPart))
    flag = 1;

  MPI_Allreduce(&flag, &flag_all, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  /*Increase storage for newly added gas particles */
  if(flag_all)
    domain_resize_storage(count_leaves, count_leaves, 0);

  int count_leaves_old = count_leaves;
  count_leaves = 0;
  int ind;

  for(ind = 0; ind < NTopleaves; ind++)
    {
      int index = Ngb_DomainNodeIndex[ind];

      if(DomainTask[ind] == ThisTask)
        {
          int cellsPerDim = 1 << ((topNodeLevel[ind] + 1 - Ngb_Nodes[index].level));

          double corner[3];
          for(i = 0; i < 3; i++)
            {
              corner[i] = Ngb_Nodes[index].Center[i] - amr_length[Ngb_Nodes[index].level + 1] + amr_length[topNodeLevel[ind] + 2];
            }

          for(i = 0; i < cellsPerDim; i++)
            {
              for(j = 0; j < cellsPerDim; j++)
                {
                  for(k = 0; k < cellsPerDim; k++)
                    {

                      P[NumGas].Pos[0] = corner[0] + i * amr_length[topNodeLevel[ind] + 1];
                      P[NumGas].Pos[1] = corner[1] + j * amr_length[topNodeLevel[ind] + 1];
                      P[NumGas].Pos[2] = corner[2] + k * amr_length[topNodeLevel[ind] + 1];

                      P[NumGas].TimeBin = 0;

                      P[NumGas].Type = 0;
                      P[NumGas].SofteningType = 0;
                      // this would seem to put the new ID at the right spot
                      P[NumGas].ID = newid++;
                      SphP[NumGas].Volume = amr_volume[topNodeLevel[ind] + 1];

                      NumGas++;
                      NumPart++;
                      count_leaves++;
                    }
                }
            }
        }
    }
  assert(count_leaves == count_leaves_old);
  myfree(topNodeLevel);

  int refiter;

  do
    {
      for(i = 0; i < NumGas; i++)
        {
          if(P[i].ID < amr_IDNew)
            continue;

          P[i].Vel[0] = 0;
          P[i].Vel[1] = 0;
          P[i].Vel[2] = 0;

          P[i].Mass = 0;
          P[i].TimeBin = 0;

#ifdef METALS
          SphP[i].MassMetallicity = 0.0;
#endif

#ifdef REFINEMENT_RPS
          SphP[i].RPSGalaxyMass = 0.0;
#endif

#ifdef MHD
          SphP[i].B[0] = 0;
          SphP[i].B[1] = 0;
          SphP[i].B[2] = 0;
          SphP[i].divB = 0;
#endif

          SphP[i].Utherm = 0;
          SphP[i].Energy = 0;
          SphP[i].Momentum[0] = 0;
          SphP[i].Momentum[1] = 0;
          SphP[i].Momentum[2] = 0;
        }

      ngb_treefree();

      ngb_treeallocate();
      ngb_treebuild(NumGas);

      distribute_sph_particles();
      spread_sph_particles();

      ngb_recompute_nodes();

      create_mesh();

      refiter = 0;

      do
        {
          do_derefinements_and_refinements();
          ngb_recompute_nodes();
          refiter++;
        }
      while(amr_Refined > 0 || amr_Derefined > 0);
    }
  while(refiter > 1);

  int count_elim = 0, count_elim_all;
  for(i = 0; i < NumGas; i++)
    if(P[i].Type == 0)
      {
        if(P[i].ID <= maxid)
          {
            // remove particle i by swapping in the last sph particle
            // and then swap the last particle to that spot
            P[i] = P[NumGas - 1];
            P[NumGas - 1] = P[NumPart - 1];

            SphP[i] = SphP[NumGas - 1];

            NumPart--;
            NumGas--;
            i--;

            count_elim++;
          }
        else
          {
            if(P[i].Mass > 0)
              {
                SphP[i].Utherm = SphP[i].Energy / P[i].Mass;
                P[i].Vel[0] = SphP[i].Momentum[0] / P[i].Mass;
                P[i].Vel[1] = SphP[i].Momentum[1] / P[i].Mass;
                P[i].Vel[2] = SphP[i].Momentum[2] / P[i].Mass;

#ifdef METALS
                SphP[i].Metallicity = SphP[i].MassMetallicity / P[i].Mass;
#endif
              }
            else
              {
#ifdef METALS
                SphP[i].Metallicity = 0.0;
#endif
              }
          }
      }

  MPI_Allreduce(&count_elim, &count_elim_all, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  sumup_large_ints(1, &NumPart, &All.TotNumPart);
  sumup_large_ints(1, &NumGas, &All.TotNumGas);

  mpi_printf("\nADD BACKGROUND GRID: count_elim_all=%d  IDNew=%d\n", count_elim_all, amr_IDNew);
  mpi_printf("ADD BACKGROUND GRID: added particles=%d  (task 0: NumGas=%d)\n", count_leaves_all - count_elim_all, NumGas);
  mpi_printf("ADD BACKGROUND GRID: new particle number=%d\n", All.TotNumPart);
  mpi_printf("ADD BACKGROUND GRID: new gas particle number=%d\n\n", All.TotNumGas);

  update_primitive_variables();

  savepositions(0, 0);

  return 0;
}

/*****************************/

static struct distributedata_in
{
  MyDouble Pos[3];
  MyFloat Vel[3];
  MyFloat Hsml;
  int NodeList[NODELISTLENGTH];
} *DistributeDataIn, *DistributeDataGet;

static struct distributedata_out
{
  MyFloat Rho;
  MyFloat DhsmlDensity;
  MyFloat Ngb;
} *DistributeDataResult, *DistributeDataOut;

static void distribute_sph_particles(void)
{
  MyFloat *Left, *Right;
  int i, j, ndone, ndone_flag, npleft, dummy, iter = 0;
  int ngrp, recvTask, place, nexport, nimport;
  long long ntot;
  double desnumngb, desnumngbdev;

  mpi_printf("ADD BACKGROUND GRID: distribution of fluid quantities in a SPH-like fashion\n");
  mpi_printf("ADD BACKGROUND GRID: finding the normalization factors\n");

  /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
  FirstActiveParticle = 0;

  for(i = 0; i < NumGas; i++)
    {
      P[i].TimeBin = 0;
      NextActiveParticle[i] = i + 1;
    }

  NextActiveParticle[NumGas - 1] = -1;

  NumNgb = (MyFloat *) mymalloc("NumNgb", NumPart * sizeof(MyFloat));
  DhsmlDensityFactor = (MyFloat *) mymalloc("DhsmlDensityFactor", NumPart * sizeof(MyFloat));

  Left = (MyFloat *) mymalloc("Left", NumPart * sizeof(MyFloat));
  Right = (MyFloat *) mymalloc("Right", NumPart * sizeof(MyFloat));

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      Left[i] = Right[i] = 0;

    }

  /* allocate buffers to arrange communication */

  Ngblist = (int *) mymalloc("Ngblist", NumPart * sizeof(int));

  All.BunchSize = (int) ((All.BufferSize * 1024 * 1024)
                         / (sizeof(data_index) + sizeof(struct data_nodelist) + sizeof(struct distributedata_in) +
                            sizeof(struct distributedata_out) + sizemax(sizeof(struct distributedata_in), sizeof(struct distributedata_out))));
  DataIndexTable = (data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(data_index));
  DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

  desnumngb = All.DesNumNgb;

  i = FirstActiveParticle;      /* begin with this index */

  do
    {
      for(j = 0; j < NTask; j++)
        {
          Send_count[j] = 0;
          Exportflag[j] = -1;
        }

      /* do local particles and prepare export list */

      for(nexport = 0; i >= 0; i = NextActiveParticle[i])
        {
          if(P[i].Type == 0 && P[i].ID < amr_IDNew && P[i].TimeBin == 0 && P[i].ID != 0)
            {
              if(distribute_evaluate(i, 0, &nexport, Send_count) < 0)
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

      DistributeDataGet = (struct distributedata_in *) mymalloc("DistributeDataGet", nimport * sizeof(struct distributedata_in));
      DistributeDataIn = (struct distributedata_in *) mymalloc("DistributeDataIn", nexport * sizeof(struct distributedata_in));

      /* prepare particle data for export */
      for(j = 0; j < nexport; j++)
        {
          place = DataIndexTable[j].Index;

          DistributeDataIn[j].Pos[0] = P[place].Pos[0];
          DistributeDataIn[j].Pos[1] = P[place].Pos[1];
          DistributeDataIn[j].Pos[2] = P[place].Pos[2];
          DistributeDataIn[j].Hsml = SphP[place].Hsml;

          memcpy(DistributeDataIn[j].NodeList, DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
        }

      /* exchange particle data */

      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
          recvTask = ThisTask ^ ngrp;

          if(recvTask < NTask)
            {
              if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                {
                  /* get the particles */
                  MPI_Sendrecv(&DistributeDataIn[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct distributedata_in), MPI_BYTE, recvTask,
                               TAG_DENS_A, &DistributeDataGet[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct distributedata_in), MPI_BYTE, recvTask,
                               TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }

      myfree(DistributeDataIn);
      DistributeDataResult = (struct distributedata_out *) mymalloc("DistributeDataResult", nimport * sizeof(struct distributedata_out));
      DistributeDataOut = (struct distributedata_out *) mymalloc("DistributeDataResult", nexport * sizeof(struct distributedata_out));

      /* now do the particles that were sent to us */

      for(j = 0; j < nimport; j++)
        distribute_evaluate(j, 1, &dummy, &dummy);

      if(i < 0)
        ndone_flag = 1;
      else
        ndone_flag = 0;

      MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      /* get the result */

      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
          recvTask = ThisTask ^ ngrp;
          if(recvTask < NTask)
            {
              if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                {
                  /* send the results */
                  MPI_Sendrecv(&DistributeDataResult[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct distributedata_out), MPI_BYTE, recvTask,
                               TAG_DENS_B, &DistributeDataOut[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct distributedata_out), MPI_BYTE, recvTask,
                               TAG_DENS_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }

        }

      for(j = 0; j < nexport; j++)
        {
          place = DataIndexTable[j].Index;

          NumNgb[place] += DistributeDataOut[j].Ngb;

          if(P[place].Type == 0)
            {
              SphP[place].Density += DistributeDataOut[j].Rho;
              DhsmlDensityFactor[place] += DistributeDataOut[j].DhsmlDensity;
            }
        }

      myfree(DistributeDataOut);
      myfree(DistributeDataResult);
      myfree(DistributeDataGet);
    }
  while(ndone < NTask);

  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);
  myfree(Right);
  myfree(Left);
  myfree(DhsmlDensityFactor);
  myfree(NumNgb);

  mpi_printf("ADD BACKGROUND GRID: done\n");
}

static int distribute_evaluate(int target, int mode, int *nexport, int *nsend_local)
{
  int j, n;
  int startnode, numngb, numngb_inbox, listindex = 0;
  double h, h2, hinv, hinv3, hinv4;
  MyFloat rho;
  double wk, dwk;
  double dx, dy, dz, r, r2, u;
  MyFloat weighted_numngb;
  MyFloat dhsmlrho;
  MyDouble *pos;

  rho = weighted_numngb = dhsmlrho = 0;

  if(mode == 0)
    {
      pos = P[target].Pos;
      h = SphP[target].Hsml;
    }
  else
    {
      pos = DistributeDataGet[target].Pos;
      h = DistributeDataGet[target].Hsml;
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
      startnode = DistributeDataGet[target].NodeList[0];
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

              if(P[j].Mass == 0)
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
                        {
                          wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
                          dwk = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
                        }
                      else
                        {
                          wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
                          dwk = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
                        }

                      double vol;

                      vol = SphP[j].Volume;

                      rho += FLT(vol * wk);

                      weighted_numngb += FLT(NORM_COEFF * wk / hinv3);  /* 4.0/3 * PI = 4.188790204786 */

                      dhsmlrho += FLT(-vol * (NUMDIMS * hinv * wk + u * dwk));
                    }
                }
            }
        }

      if(mode == 1)
        {
          listindex++;
          if(listindex < NODELISTLENGTH)
            {
              startnode = DistributeDataGet[target].NodeList[listindex];
              if(startnode >= 0)
                startnode = Ngb_Nodes[startnode].u.d.nextnode;  /* open it */
            }
        }
    }

  if(mode == 0)
    {
      NumNgb[target] = weighted_numngb;
      SphP[target].Density = rho;
      DhsmlDensityFactor[target] = dhsmlrho;
    }
  else
    {
      DistributeDataResult[target].Rho = rho;
      DistributeDataResult[target].Ngb = weighted_numngb;
      DistributeDataResult[target].DhsmlDensity = dhsmlrho;
    }

  return 0;
}

/******************/

static struct spreaddata_in
{
  MyDouble Pos[3];
  MyFloat Hsml;
  MyFloat SumWeight;

  MyFloat Mass;
  MyFloat Momentum[3];
  MyFloat Etherm;

#ifdef METALS
  MyFloat Metallicity;
#endif

#ifdef REFINEMENT_RPS
  MyFloat RPSGalaxyMass;
#endif

#ifdef MHD
  MyFloat B[3];
#endif

  int NodeList[NODELISTLENGTH];
} *SpreadDataIn, *SpreadDataGet;

static void spread_sph_particles(void)
{
  int i, j, ndone, ndone_flag, dummy;
  int ngrp, recvTask, place, nexport, nimport;

  /* allocate buffers to arrange communication */

  Ngblist = (int *) mymalloc("Ngblist", NumPart * sizeof(int));

  All.BunchSize = (int) ((All.BufferSize * 1024 * 1024) / (sizeof(data_index) + sizeof(struct data_nodelist) + sizeof(struct spreaddata_in) + sizeof(struct spreaddata_in)));
  DataIndexTable = (data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(data_index));
  DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

  mpi_printf("ADD BACKGROUND GRID: distributing the fluid quantities\n");

  i = FirstActiveParticle;      /* begin with this index */

  do
    {
      for(j = 0; j < NTask; j++)
        {
          Send_count[j] = 0;
          Exportflag[j] = -1;
        }

      /* do local particles and prepare export list */

      for(nexport = 0; i >= 0; i = NextActiveParticle[i])
        {
          if(P[i].Type == 0 && P[i].ID < amr_IDNew && P[i].ID != 0)
            {
              if(spread_evaluate(i, 0, &nexport, Send_count) < 0)
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

      SpreadDataGet = (struct spreaddata_in *) mymalloc("SpreadDataGet", nimport * sizeof(struct spreaddata_in));
      SpreadDataIn = (struct spreaddata_in *) mymalloc("SpreadDataIn", nexport * sizeof(struct spreaddata_in));

      /* prepare particle data for export */
      for(j = 0; j < nexport; j++)
        {
          place = DataIndexTable[j].Index;

          SpreadDataIn[j].Pos[0] = P[place].Pos[0];
          SpreadDataIn[j].Pos[1] = P[place].Pos[1];
          SpreadDataIn[j].Pos[2] = P[place].Pos[2];
          SpreadDataIn[j].Hsml = SphP[place].Hsml;
          SpreadDataIn[j].SumWeight = SphP[place].Density;
          SpreadDataIn[j].Mass = P[place].Mass;
          SpreadDataIn[j].Etherm = SphP[place].Utherm * P[place].Mass;
          SpreadDataIn[j].Momentum[0] = P[place].Vel[0] * P[place].Mass;
          SpreadDataIn[j].Momentum[1] = P[place].Vel[1] * P[place].Mass;
          SpreadDataIn[j].Momentum[2] = P[place].Vel[2] * P[place].Mass;
#ifdef MHD
          SpreadDataIn[j].B[0] = SphP[place].B[0];
          SpreadDataIn[j].B[1] = SphP[place].B[1];
          SpreadDataIn[j].B[2] = SphP[place].B[2];
#endif
#ifdef METALS
          SpreadDataIn[j].Metallicity = SphP[place].Metallicity;
#endif
#ifdef REFINEMENT_RPS
          SpreadDataIn[j].RPSGalaxyMass = SphP[place].RPSGalaxyMass;
#endif

          memcpy(SpreadDataIn[j].NodeList, DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
        }

      /* exchange particle data */

      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
          recvTask = ThisTask ^ ngrp;

          if(recvTask < NTask)
            {
              if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                {
                  /* get the particles */
                  MPI_Sendrecv(&SpreadDataIn[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct spreaddata_in), MPI_BYTE, recvTask, TAG_DENS_A,
                               &SpreadDataGet[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct spreaddata_in), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }

      myfree(SpreadDataIn);

      /* now do the particles that were sent to us */

      for(j = 0; j < nimport; j++)
        spread_evaluate(j, 1, &dummy, &dummy);

      if(i < 0)
        ndone_flag = 1;
      else
        ndone_flag = 0;

      MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      myfree(SpreadDataGet);
    }
  while(ndone < NTask);

#ifdef MHD
  /* now divide the B field in each cell by the weight (sum of the wk's,
     which we stored in SphP.divB */
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    if(P[i].Type == 0 && P[i].ID < amr_IDNew && P[i].ID != 0)
      {
        for(j = 0; j < 3; j++)
          if(SphP[i].divB > 0)
            SphP[i].B[j] /= SphP[i].divB;
      }
#endif

  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);

  mpi_printf("ADD BACKGROUND GRID: done\n");
}

static int spread_evaluate(int target, int mode, int *nexport, int *nsend_local)
{
  int j, n;
  int startnode, numngb, numngb_inbox, listindex = 0;
  double h, h2, hinv, hinv3;
  double wk;
  double dx, dy, dz, r, r2, u;
  MyDouble *pos;
  double mass, sumweight, etherm, momentum[3], weight;
#ifdef MHD
  double B[3];
#endif

#ifdef METALS
  MyFloat metallicity;
#endif
#ifdef REFINEMENT_RPS
  MyFloat RPSGalaxyMass;
#endif

  if(mode == 0)
    {
      pos = P[target].Pos;
      h = SphP[target].Hsml;

      sumweight = SphP[target].Density;
      mass = P[target].Mass;
      etherm = SphP[target].Utherm * P[target].Mass;
      momentum[0] = P[target].Vel[0] * P[target].Mass;
      momentum[1] = P[target].Vel[1] * P[target].Mass;
      momentum[2] = P[target].Vel[2] * P[target].Mass;
#ifdef METALS
      metallicity = SphP[target].Metallicity;
#endif
#ifdef REFINEMENT_RPS
      RPSGalaxyMass = SphP[target].RPSGalaxyMass;
#endif
#ifdef MHD
      B[0] = SphP[target].B[0];
      B[1] = SphP[target].B[1];
      B[2] = SphP[target].B[2];
#endif
    }
  else
    {
      pos = SpreadDataGet[target].Pos;
      h = SpreadDataGet[target].Hsml;

      sumweight = SpreadDataGet[target].SumWeight;
      mass = SpreadDataGet[target].Mass;
      etherm = SpreadDataGet[target].Etherm;
      momentum[0] = SpreadDataGet[target].Momentum[0];
      momentum[1] = SpreadDataGet[target].Momentum[1];
      momentum[2] = SpreadDataGet[target].Momentum[2];
#ifdef METALS
      metallicity = SpreadDataGet[target].Metallicity;
#endif
#ifdef REFINEMENT_RPS
      RPSGalaxyMass = SpreadDataGet[target].RPSGalaxyMass;
#endif
#ifdef MHD
      B[0] = SpreadDataGet[target].B[0];
      B[1] = SpreadDataGet[target].B[1];
      B[2] = SpreadDataGet[target].B[2];
#endif
    }

  h2 = h * h;
  hinv = 1.0 / h;
#ifndef  TWODIMS
  hinv3 = hinv * hinv * hinv;
#else
  hinv3 = hinv * hinv / boxSize_Z;
#endif

  if(mode == 0)
    {
      startnode = Ngb_MaxPart;  /* root node */
    }
  else
    {
      startnode = SpreadDataGet[target].NodeList[0];
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

              if(P[j].ID >= amr_IDNew)
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

                      double vol;

                      vol = SphP[j].Volume;

                      weight = FLT(vol * wk) / sumweight;

                      P[j].Mass += mass * weight;
                      SphP[j].Energy += etherm * weight;
                      SphP[j].Momentum[0] += momentum[0] * weight;
                      SphP[j].Momentum[1] += momentum[1] * weight;
                      SphP[j].Momentum[2] += momentum[2] * weight;
#ifdef METALS
                      SphP[j].MassMetallicity += mass * weight * metallicity;
#endif
#ifdef REFINEMENT_RPS
                      SphP[j].RPSGalaxyMass += RPSGalaxyMass * weight;
#endif

#ifdef MHD
                      /* B is an intrinsic quantity, so we want to assign the kernel-weighted
                         mean value to a cell rather than redistribute as done for, e.g., mass */
                      SphP[j].B[0] += wk * B[0];
                      SphP[j].B[1] += wk * B[1];
                      SphP[j].B[2] += wk * B[2];
                      SphP[j].divB += wk;
#endif
                    }
                }
            }
        }

      if(mode == 1)
        {
          listindex++;
          if(listindex < NODELISTLENGTH)
            {
              startnode = SpreadDataGet[target].NodeList[listindex];
              if(startnode >= 0)
                startnode = Ngb_Nodes[startnode].u.d.nextnode;  /* open it */
            }
        }
    }

  return 0;
}

void init_aux_fields(void)
{
#ifdef METALS
  int i;
  for(i = 0; i < NumGas; i++)
    {
      SphP[i].Metallicity = P[i].Metallicity;   /* read from IC */
      SphP[i].MassMetallicity = SphP[i].Metallicity * P[i].Mass;
    }
#endif
}

#endif
