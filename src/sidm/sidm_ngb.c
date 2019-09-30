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
  MyDouble Hsml;
  double Mass;
  unsigned char State;
  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

typedef struct
{
  int Ngb;
  ngb_entry E[SIDM_MAX_NGBS];
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

  in->Firstnode = firstnode;
}

static void out2particle(data_out * out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES)      /* initial store */
    {
      memcpy(&PSIDM[LISIDM[i].List2].ngb_Entry[PSIDM[LISIDM[i].List2].ngb_Offset], &(out->E[0]), out->Ngb * sizeof(ngb_entry));
      PSIDM[LISIDM[i].List2].ngb_Offset = out->Ngb;
    }
  else                          /* combine */
    {
      memcpy(&PSIDM[LISIDM[i].List2].ngb_Entry[PSIDM[LISIDM[i].List2].ngb_Offset], &(out->E[0]), out->Ngb * sizeof(ngb_entry));
      PSIDM[LISIDM[i].List2].ngb_Offset += out->Ngb;
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
          sidm_NgbList_evaluate(idx, MODE_LOCAL_PARTICLES, threadid);
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

        sidm_NgbList_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}


void sidm_NgbList(void)
{
  mpi_printf("SIDM: Finding NGB list\n");

  generic_set_MaxNexport();

  generic_comm_pattern(Nforces, kernel_local, kernel_imported);

  mpi_printf("SIDM: done with NGB list\n");
}


void sidm_NgbList_evaluate(int target, int mode, int threadid)
{
  int numnodes, *firstnode;
  data_in local, *in;
  data_out out;
  int k;
  int no;
#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif
  double wk, h, h2, hinv, hinv3;
  double r, r2, u;
  double in_mass_1, in_mass_2;
  MyDouble phys_rho;
  MyDouble *pos;
  MyFloat *vel;
  MyDouble phys_rel_velx, phys_rel_vely, phys_rel_velz;
  MyDouble phys_rel_vel;
  MyDouble dx, dy, dz;
  unsigned char pstate, state;
  unsigned char reaction;
  double Ekin;
  int retval;

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
  h = in->Hsml;
  in_mass_1 = in->Mass;
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
              dx = GRAVITY_NEAREST_X(Tree_Pos_list[3 * no + 0] - pos[0]);
              dy = GRAVITY_NEAREST_Y(Tree_Pos_list[3 * no + 1] - pos[1]);
              dz = GRAVITY_NEAREST_Z(Tree_Pos_list[3 * no + 2] - pos[2]);

              r2 = dx * dx + dy * dy + dz * dz;


              if((r2 < h2) && (r2 > 0) && ((1 << P[no].Type) & (SIDM)))
                {
                  hinv = 1. / h;
                  hinv3 = hinv * hinv * hinv;

                  r = sqrt(r2);

                  u = r * hinv;

                  if(u < 0.5)
                    {
                      wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
                    }
                  else
                    {
                      wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
                    }



                  state = P[no].sidm_State;

                  in_mass_2 = P[no].Mass;

                  phys_rho = in_mass_2 * wk;
                  phys_rho *= SIDM_arho;
                  phys_rel_velx = SIDM_avel * (P[no].Vel[0] - vel[0]);
                  phys_rel_vely = SIDM_avel * (P[no].Vel[1] - vel[1]);
                  phys_rel_velz = SIDM_avel * (P[no].Vel[2] - vel[2]);
                  phys_rel_vel = sqrt(phys_rel_velx * phys_rel_velx + phys_rel_vely * phys_rel_vely + phys_rel_velz * phys_rel_velz);

                  Ekin = SIDM_avel * SIDM_avel * (0.5 * in_mass_1 * (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]) + 0.5 * in_mass_2 * (P[no].Vel[0] * P[no].Vel[0] + P[no].Vel[1] * P[no].Vel[1] + P[no].Vel[2] * P[no].Vel[2]));

                  for(reaction = 0; reaction < SIDM_REACTIONS; reaction++)
                    if((pstate == SMSIDM[reaction].In1) && (state == SMSIDM[reaction].In2))
                      out.E[out.Ngb].P0j_half[reaction] = sidm_scatter_P(phys_rho, phys_rel_vel, Ekin, reaction, &retval); 

                  out.E[out.Ngb].NgbIDs = P[no].ID;
                  out.E[out.Ngb].Distance = r;
                  out.E[out.Ngb].State = P[no].sidm_State;

                  out.Ngb++;


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

              dx = GRAVITY_NEAREST_X(Tree_Points[n].Pos[0] - pos[0]);
              dy = GRAVITY_NEAREST_Y(Tree_Points[n].Pos[1] - pos[1]);
              dz = GRAVITY_NEAREST_Z(Tree_Points[n].Pos[2] - pos[2]);

              r2 = dx * dx + dy * dy + dz * dz;

              if((r2 < h2) && (r2 > 0) && ((1 << (Tree_Points[n].Type)) & (SIDM)))
                {
                  hinv = 1. / h;
                  hinv3 = hinv * hinv * hinv;

                  r = sqrt(r2);

                  u = r * hinv;

                  if(u < 0.5)
                    {
                      wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
                    }
                  else
                    {
                      wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
                    }


                  state = Tree_Points[n].sidm_State;

                  in_mass_2 = Tree_Points[n].Mass;

                  phys_rho = in_mass_2 * wk;
                  phys_rho *= SIDM_arho;
                  phys_rel_velx = SIDM_avel * (Tree_Points[n].Vel[0] - vel[0]);
                  phys_rel_vely = SIDM_avel * (Tree_Points[n].Vel[1] - vel[1]);
                  phys_rel_velz = SIDM_avel * (Tree_Points[n].Vel[2] - vel[2]);
                  phys_rel_vel = sqrt(phys_rel_velx * phys_rel_velx + phys_rel_vely * phys_rel_vely + phys_rel_velz * phys_rel_velz);

                  Ekin = SIDM_avel * SIDM_avel * (0.5 * in_mass_1 * (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]) + 0.5 * in_mass_2 * (Tree_Points[n].Vel[0] * Tree_Points[n].Vel[0] + Tree_Points[n].Vel[1] * Tree_Points[n].Vel[1] + Tree_Points[n].Vel[2] * Tree_Points[n].Vel[2]));

                  for(reaction = 0; reaction < SIDM_REACTIONS; reaction++)
                    if((pstate == SMSIDM[reaction].In1) && (state == SMSIDM[reaction].In2))
                      out.E[out.Ngb].P0j_half[reaction] = sidm_scatter_P(phys_rho, phys_rel_vel, Ekin, reaction, &retval); 

                  out.E[out.Ngb].NgbIDs = Tree_Points[n].sidm_ID;
                  out.E[out.Ngb].Distance = r;
                  out.E[out.Ngb].State = Tree_Points[n].sidm_State;
                  out.Ngb++;

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


  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;


}


#endif
