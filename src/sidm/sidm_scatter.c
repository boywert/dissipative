/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/sidm/sidm_neighbors.c
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
#include "sidm_vars.h"

#ifdef SIDM



typedef struct
{
  MyDouble Pos[3];
  MyFloat Vel[3];
  MyDouble Mass;
  MyDouble Hsml;
  MyIDType ScatterID;
  unsigned char State;
  int Reaction;
  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

typedef struct
{
  float Vx, Vy, Vz;
  unsigned char State;
  double Mass;
  int ScattersInStep[SIDM_REACTIONS];
  double Ekin_before[SIDM_REACTIONS], Ekin_after[SIDM_REACTIONS], Scattered_Mass[SIDM_REACTIONS];
} data_out;

static data_out *DataResult, *DataOut;





static void particle2in(data_in * in, int i, int firstnode)
{
  int k;

  int idx = TargetList[i];

  if(idx < NumPart)
    {
      in->State = P[idx].sidm_State;
      in->Mass = P[idx].Mass;
      for(k = 0; k < 3; k++)
        {
          in->Pos[k] = P[idx].Pos[k];
          in->Vel[k] = P[idx].Vel[k];
        }
    }
  else
    {
      idx -= Tree_ImportedNodeOffset;
      in->State = Tree_Points[idx].sidm_State;
      in->Mass = Tree_Points[idx].Mass;
      for(k = 0; k < 3; k++)
        {
          in->Pos[k] = Tree_Points[idx].Pos[k];
          in->Vel[k] = Tree_Points[idx].Vel[k];
        }
    }
  in->Hsml = PSIDM[i].Hsml;
  in->ScatterID = PSIDM[i].ScatterID;
  in->Reaction = PSIDM[i].ScatterReaction;
  in->Firstnode = firstnode;
}

static void out2particle(data_out * out, int i, int mode)
{
  int target = TargetList[i];
  int totScattersInStep;
  unsigned char reaction;

  if(target > NumPart)
    target = NumPart + target - Tree_ImportedNodeOffset;

  if(mode == MODE_LOCAL_PARTICLES)      /* initial store */
    {
      for(totScattersInStep = 0, reaction = 0; reaction < SIDM_REACTIONS; reaction++)
        {
          TSIDM[target].ScattersInStep[reaction] = out->ScattersInStep[reaction];
          if(out->ScattersInStep[reaction] > 0)
            TSIDM[target].ScatterReaction = reaction;
          totScattersInStep += out->ScattersInStep[reaction];
          if(TSIDM[target].ScattersInStep[reaction] > 1)
            terminate("SIDM: multiple scatter -> reaction=%d   out->ScattersInStep[reaction]=%d  target=%d  TSIDM[target].ScattersInStep[target] = %d\n", reaction, out->ScattersInStep[reaction],
                      target, TSIDM[target].ScattersInStep[reaction]);
        }
      if(totScattersInStep > 0)
        {
          TSIDM[target].NewVel[0] = out->Vx;
          TSIDM[target].NewVel[1] = out->Vy;
          TSIDM[target].NewVel[2] = out->Vz;
          TSIDM[target].NewState  = out->State;
          TSIDM[target].NewMass = out->Mass;
          for (reaction = 0; reaction < SIDM_REACTIONS; reaction++)
            {  
              SIDM_Ekin_before_total[reaction] += out->Ekin_before[reaction];
              SIDM_Ekin_after_total[reaction] += out->Ekin_after[reaction];
              SIDM_Scattered_Mass[reaction] += out->Scattered_Mass[reaction];
            }
        }
    }
  else                          /* combine */
    {
      for(totScattersInStep = 0, reaction = 0; reaction < SIDM_REACTIONS; reaction++)
        {
          TSIDM[target].ScattersInStep[reaction] += out->ScattersInStep[reaction];
          if(out->ScattersInStep[reaction] > 0)
            TSIDM[target].ScatterReaction = reaction;
          totScattersInStep += out->ScattersInStep[reaction];
          if(TSIDM[target].ScattersInStep[reaction] > 1)
            terminate("SIDM: mutliple scatter -> reaction=%d   out->ScattersInStep[reaction]=%d  target=%d  TSIDM[target].ScattersInStep[reaction] = %d\n", reaction, out->ScattersInStep[reaction],
                      target, TSIDM[target].ScattersInStep[reaction]);
        }

      if(totScattersInStep > 0)
        {
          TSIDM[target].NewVel[0] = out->Vx;
          TSIDM[target].NewVel[1] = out->Vy;
          TSIDM[target].NewVel[2] = out->Vz;
          TSIDM[target].NewState  = out->State;
          TSIDM[target].NewMass = out->Mass;
          for (reaction = 0; reaction < SIDM_REACTIONS; reaction++)
            {
              SIDM_Ekin_before_total[reaction] += out->Ekin_before[reaction];
              SIDM_Ekin_after_total[reaction] += out->Ekin_after[reaction];
              SIDM_Scattered_Mass[reaction] += out->Scattered_Mass[reaction];
            }
        }
    }
}



#include "../generic_comm_helpers2.h"



static void kernel_local(void)
{
  int idx;
#ifdef GENERIC_ASYNC
  int flag = 0;
#endif

#pragma omp parallel private(idx)
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
              if(generic_polling_primary(count, Nforces))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        idx = NextParticle++;

        if(idx >= Nforces)
          break;

        if(PSIDM[idx].ShouldScatterInStep[PSIDM[idx].ScatterReaction] > 0)
          sidm_Scatter_evaluate(idx, MODE_LOCAL_PARTICLES, threadid);
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

        sidm_Scatter_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}


void sidm_Scatter(void)
{
  mpi_printf("SIDM: Scatter\n");

  int i;
  unsigned char reaction;

  for(i = 0; i < (NumPart + Tree_NumPartImported); i++)
    for(reaction = 0; reaction < SIDM_REACTIONS; reaction++)
      if((TSIDM[i].ScattersInStep[reaction] > 1) || (TSIDM[i].ScattersInStep[reaction] < -1))
        terminate("SIDM: wrong scatter i=%d  reaction=%d  TSIDM[i].ScattersInStep[reaction]=%d  (NumPart + Tree_NumPartImported)=%d\n", i, reaction, TSIDM[i].ScattersInStep[reaction],
                  (NumPart + Tree_NumPartImported));


  generic_set_MaxNexport();

  generic_comm_pattern(Nforces, kernel_local, kernel_imported);

  mpi_printf("SIDM: done with scatter\n");
}


void sidm_Scatter_evaluate(int target, int mode, int threadid)
{
  int numnodes, *firstnode;
  data_in local, *in;
  data_out out;
  int k;
  int no;
#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif
  double h, h2;
  MyDouble *pos;
  MyFloat *vel;
  MyDouble dx, dy, dz;
  int didScatterInStep = 0;
  MyDouble mass;
  MyIDType scatterID;
  unsigned char pstate;
  unsigned char reaction;
  scatter_process_data_in sprdata_in;
  scatter_process_data_out sprdata_out;
  double Ekin_before = 0.0, Ekin_after = 0.0, Scattered_Mass = 0.0;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      in = &local;

      numnodes = 1;
      firstnode = NULL;
    }
  else
    {
      in = &DataGet[target];
      generic_get_numnodes(target, &numnodes, &firstnode);
    }



  memset(&out, 0, sizeof(data_out));

  pos = in->Pos;
  vel = in->Vel;
  mass = in->Mass;
  h = in->Hsml;
  scatterID = in->ScatterID;
  reaction = in->Reaction;
  pstate = in->State;

  h2 = h * h;


  for(k = 0; k < numnodes; k++)
    {
      if(mode == MODE_LOCAL_PARTICLES)
        {
          no = Tree_MaxPart;    /* root node */
        }
      else
        {
          no = firstnode[k];
          no = Nodes[no].u.d.nextnode;  /* open it */
        }


      while(no >= 0)
        {
          if(no < Tree_MaxPart) /* single particle */
            {



              if(P[no].ID == scatterID)
                {
                  if(!((1 << P[no].Type) & (SIDM)))
                    terminate("SIDM: no SIDM scatter partner");

                  Ekin_before = SIDM_avel * SIDM_avel * (0.5 * mass * (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]) + 0.5 * P[no].Mass * (P[no].Vel[0] * P[no].Vel[0] + P[no].Vel[1] * P[no].Vel[1] + P[no].Vel[2] * P[no].Vel[2]));

                  sprdata_in.Pos1[0] = pos[0];
                  sprdata_in.Pos1[1] = pos[1];
                  sprdata_in.Pos1[2] = pos[2];
                  sprdata_in.Pos2[0] = Tree_Pos_list[3 * no + 0];
                  sprdata_in.Pos2[1] = Tree_Pos_list[3 * no + 1];
                  sprdata_in.Pos2[2] = Tree_Pos_list[3 * no + 2];
                  sprdata_in.Vel1[0] = vel[0];
                  sprdata_in.Vel1[1] = vel[1];
                  sprdata_in.Vel1[2] = vel[2];
                  sprdata_in.Vel2[0] = P[no].Vel[0];
                  sprdata_in.Vel2[1] = P[no].Vel[1];
                  sprdata_in.Vel2[2] = P[no].Vel[2];
                  sprdata_in.Mass1 = mass;
                  sprdata_in.Mass2 = P[no].Mass;
                  sprdata_in.Reaction = reaction;
                  sprdata_in.State1 = pstate;
                  sprdata_in.State2 = P[no].sidm_State;
                  sidm_evaluate_scatter_process(&sprdata_in, &sprdata_out, P[no].ID);

                  /* assign new velocity */
                  TSIDM[no].NewVel[0] = sprdata_out.Vel2[0];
                  TSIDM[no].NewVel[1] = sprdata_out.Vel2[1];
                  TSIDM[no].NewVel[2] = sprdata_out.Vel2[2];
                  TSIDM[no].NewState  = sprdata_out.OutState2;
                  TSIDM[no].NewMass = sprdata_out.Mass2;
                  didScatterInStep++;
                  TSIDM[no].ScattersInStep[reaction]++;
                  TSIDM[no].ScatterReaction = reaction;
                  if(TSIDM[no].ScattersInStep[reaction] > 1)
                    terminate("SIDM: multiple scatter -> reaction=%d  TSIDM[no].ScattersInStep[reaction]=%d  no=%d  scatterID=%llu\n", reaction, TSIDM[no].ScattersInStep[reaction], no, (unsigned long long)scatterID);

                }

              no = Nextnode[no];
            }
          else if(no < Tree_MaxPart + Tree_MaxNodes)    /* internal node */
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                {
                  if(no < Tree_FirstNonTopLevelNode)    /* we reached a top-level node again, which means that we are done with the branch */
                    break;

                }

              struct NODE *current = &Nodes[no];

              no = current->u.d.sibling;        /* in case the node can be discarded */

              double dist = h + 0.5 * current->len;
              dx = NGB_PERIODIC_LONG_X(current->center[0] - pos[0]);
              if(dx > dist)
                continue;
              dy = NGB_PERIODIC_LONG_Y(current->center[1] - pos[1]);
              if(dy > dist)
                continue;
              dz = NGB_PERIODIC_LONG_Z(current->center[2] - pos[2]);
              if(dz > dist)
                continue;
              /* now test against the minimal sphere enclosing everything */
              dist += FACT1 * current->len;
              if(dx * dx + dy * dy + dz * dz > dist * dist)
                continue;

              no = current->u.d.nextnode;       /* ok, we need to open the node */

            }


          else if(no >= Tree_ImportedNodeOffset)        /* point from imported nodelist */
            {
              int n = no - Tree_ImportedNodeOffset;

              if(Tree_Points[n].sidm_ID == scatterID)
                {

                  if(!((1 << (Tree_Points[n].Type)) & (SIDM)))
                    terminate("SIDM: no SIDM scatter partner");

                  Ekin_before = SIDM_avel * SIDM_avel * (0.5 * mass * (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]) + 0.5 * Tree_Points[n].Mass * (Tree_Points[n].Vel[0] * Tree_Points[n].Vel[0] + Tree_Points[n].Vel[1] * Tree_Points[n].Vel[1] + Tree_Points[n].Vel[2] * Tree_Points[n].Vel[2]));

                  sprdata_in.Pos1[0] = pos[0];
                  sprdata_in.Pos1[1] = pos[1];
                  sprdata_in.Pos1[2] = pos[2];
                  sprdata_in.Pos2[0] = Tree_Points[n].Pos[0];
                  sprdata_in.Pos2[1] = Tree_Points[n].Pos[1];
                  sprdata_in.Pos2[2] = Tree_Points[n].Pos[2];
                  sprdata_in.Vel1[0] = vel[0];
                  sprdata_in.Vel1[1] = vel[1];
                  sprdata_in.Vel1[2] = vel[2];
                  sprdata_in.Vel2[0] = Tree_Points[n].Vel[0];
                  sprdata_in.Vel2[1] = Tree_Points[n].Vel[1];
                  sprdata_in.Vel2[2] = Tree_Points[n].Vel[2];
                  sprdata_in.Mass1 = mass;
                  sprdata_in.Mass2 = Tree_Points[n].Mass;
                  sprdata_in.Reaction = reaction;
                  sprdata_in.State1 = pstate;
                  sprdata_in.State2 = Tree_Points[n].sidm_State;
                  sidm_evaluate_scatter_process(&sprdata_in, &sprdata_out, Tree_Points[n].sidm_ID);

                  /* assign new velocity */
                  TSIDM[NumPart + n].NewVel[0] = sprdata_out.Vel2[0];
                  TSIDM[NumPart + n].NewVel[1] = sprdata_out.Vel2[1];
                  TSIDM[NumPart + n].NewVel[2] = sprdata_out.Vel2[2];
                  TSIDM[NumPart + n].NewState  = sprdata_out.OutState2;
                  TSIDM[NumPart + n].NewMass = sprdata_out.Mass2;
                  didScatterInStep++;
                  TSIDM[NumPart + n].ScattersInStep[reaction]++;
                  TSIDM[NumPart + n].ScatterReaction = reaction;
                  if(TSIDM[NumPart + n].ScattersInStep[reaction] > 1)
                    terminate("SIDM: multiple scatter -> reaction=%d  TSIDM[NumPart + n].ScattersInStep[reaction]=%d  n=%d  scatterID=%llu    NumPart=%d   (NumPart + Tree_NumPartImported)=%d\n",
                              reaction, TSIDM[NumPart + n].ScattersInStep[reaction], n, (unsigned long long)scatterID, NumPart, (NumPart + Tree_NumPartImported));
                }

              no = Nextnode[no - Tree_MaxNodes];

            }
          else                  /* pseudo particle */
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                terminate("mode == MODE_IMPORTED_PARTICLES");

              if(target >= 0)
                tree_treefind_export_node_threads(no, target, threadid);

              no = Nextnode[no - Tree_MaxNodes];
              continue;
            }
        }

    }

  Ekin_after = SIDM_avel * SIDM_avel * (0.5 * sprdata_out.Mass1 * (sprdata_out.Vel1[0] * sprdata_out.Vel1[0] + sprdata_out.Vel1[1] * sprdata_out.Vel1[1] + sprdata_out.Vel1[2] * sprdata_out.Vel1[2]) +
               0.5 * sprdata_out.Mass2 * (sprdata_out.Vel2[0] * sprdata_out.Vel2[0] + sprdata_out.Vel2[1] * sprdata_out.Vel2[1] + sprdata_out.Vel2[2] * sprdata_out.Vel2[2]));
  Scattered_Mass = sprdata_out.Mass1 + sprdata_out.Mass2;
  out.Vx = sprdata_out.Vel1[0];
  out.Vy = sprdata_out.Vel1[1];
  out.Vz = sprdata_out.Vel1[2];
  out.State = sprdata_out.OutState1;
  out.Mass = sprdata_out.Mass1;
  out.Ekin_before[reaction] = Ekin_before;
  out.Ekin_after[reaction] = Ekin_after;
  out.Scattered_Mass[reaction] = Scattered_Mass;
  out.ScattersInStep[reaction] = didScatterInStep;

  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;
}


#endif
