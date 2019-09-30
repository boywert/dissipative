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
  /* your fields go here */
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
  /* your fields go here */
  MyDouble Density[SIDM_STATES];
  MyDouble VelDisp[SIDM_STATES], Vx[SIDM_STATES], Vy[SIDM_STATES], Vz[SIDM_STATES];
  MyDouble PSum[SIDM_REACTIONS];
  int Ngb, NgbStates[SIDM_STATES];
  int EnergyForbidden[SIDM_REACTIONS];
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
  unsigned char state;
  unsigned char reaction;
  if(mode == MODE_LOCAL_PARTICLES)      /* initial store */
    {
      for(state = 0; state < SIDM_STATES; state++)
        {
          PSIDM[i].Density[state] = out->Density[state];
          PSIDM[i].Vx[state] = out->Vx[state];
          PSIDM[i].Vy[state] = out->Vy[state];
          PSIDM[i].Vz[state] = out->Vz[state];
          PSIDM[i].VelDisp[state] = out->VelDisp[state];
          PSIDM[i].NumNgbState[state] = out->NgbStates[state];
        }
      for(reaction = 0; reaction < SIDM_REACTIONS; reaction++)
        {
          PSIDM[i].PSum[reaction] = out->PSum[reaction];
          PSIDM[i].EnergyForbidden[reaction] = out->EnergyForbidden[reaction];
        }

      PSIDM[i].NumNgb = out->Ngb;
    }
  else                          /* combine */
    {
      for(state = 0; state < SIDM_STATES; state++)
        {
          PSIDM[i].Density[state] += out->Density[state];
          PSIDM[i].Vx[state] += out->Vx[state];
          PSIDM[i].Vy[state] += out->Vy[state];
          PSIDM[i].Vz[state] += out->Vz[state];
          PSIDM[i].VelDisp[state] += out->VelDisp[state];
          PSIDM[i].NumNgbState[state] += out->NgbStates[state];
        }
      for(reaction = 0; reaction < SIDM_REACTIONS; reaction++)
        {
          PSIDM[i].PSum[reaction] += out->PSum[reaction];
          PSIDM[i].EnergyForbidden[reaction] += out->EnergyForbidden[reaction];
        }

      PSIDM[i].NumNgb += out->Ngb;
    }
}



#include "../generic_comm_helpers2.h"

static unsigned char *Todo;

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

        if(Todo[idx] > 0)
          sidm_findHsml_evaluate(idx, MODE_LOCAL_PARTICLES, threadid);
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

        sidm_findHsml_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}


void sidm_findHsml(void)
{
  int idx, i;
  int iter = 0, npleft;
  MyDouble *Left, *Right;
  double t0, t1;
  long long ntot;

  mpi_printf("SIDM: Finding Hsml values\n");

  Left = mymalloc("Left", Nforces * sizeof(MyDouble));
  Right = mymalloc("Right", Nforces * sizeof(MyDouble));
  Todo = mymalloc("Todo", Nforces * sizeof(unsigned char));

  int nforces = 0;
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(((1 << P[i].Type) & (SIDM)) && Tree_Task_list[i] == ThisTask)
        {
          Left[nforces] = 0;
          Right[nforces] = 0;
          Todo[nforces] = 1;
          PSIDM[nforces].ID = P[i].ID;
          PSIDM[nforces++].Hsml = P[i].sidm_Hsml;
        }
    }

  for(i = 0; i < Tree_NumPartImported; i++)
#ifndef HIERARCHICAL_GRAVITY
    if(Tree_Points[i].ActiveFlag)
#endif
    if(((1 << (Tree_Points[i].Type)) & (SIDM)))
      {
        Left[nforces] = 0;
        Right[nforces] = 0;
        Todo[nforces] = 1;
        PSIDM[nforces].ID = Tree_Points[i].sidm_ID;
        PSIDM[nforces++].Hsml = Tree_Points[i].sidm_Hsml;
      }


  generic_set_MaxNexport();


  do
    {

      t0 = second();



      generic_comm_pattern(Nforces, kernel_local, kernel_imported);



      /* do final operations on results */
      for(i = 0, npleft = 0; i < Nforces; i++)
        {

          if(Todo[i])
            {
              if((PSIDM[i].NumNgb < (All.SIDMDesNumNgb - All.SIDMMaxNumNgbDeviation)) || (PSIDM[i].NumNgb > (All.SIDMDesNumNgb + All.SIDMMaxNumNgbDeviation)))
                {
                  /* need to redo this particle */
                  npleft++;

                  if(PSIDM[i].NumNgb < All.SIDMDesNumNgb - All.SIDMMaxNumNgbDeviation)
                    Left[i] = dmax(PSIDM[i].Hsml, Left[i]);
                  else
                    {
                      if(Right[i] != 0)
                        {
                          if(PSIDM[i].Hsml < Right[i])
                            Right[i] = PSIDM[i].Hsml;
                        }
                      else
                        Right[i] = PSIDM[i].Hsml;
                    }

                  if(iter >= MAXITER - 10)
                    {
                      printf("SIDM: i=%d task=%d Hsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g\n", i, ThisTask, PSIDM[i].Hsml, Left[i], Right[i], (double) PSIDM[i].NumNgb, Right[i] - Left[i]);
                      fflush(stdout);
                    }

                  if(Right[i] > 0 && Left[i] > 0)
                    PSIDM[i].Hsml = pow(0.5 * (pow(Left[i], 3) + pow(Right[i], 3)), 1.0 / 3);
                  else
                    {
                      if(Right[i] == 0 && Left[i] == 0)
                        terminate("BAD");       /* can't occur */

                      if(Right[i] == 0 && Left[i] > 0)
                        PSIDM[i].Hsml *= 1.26;

                      if(Right[i] > 0 && Left[i] == 0)
                        PSIDM[i].Hsml /= 1.26;
                    }

                }
              else
                {
                  Todo[i] = 0;
                  if(PSIDM[i].NumNgb > SIDM_MAX_NGBS || PSIDM[i].NumNgb < SIDM_MIN_NGBS)
                    terminate("SIDM: numngb = %d", PSIDM[i].NumNgb);
                }

            }

        }

      sumup_large_ints(1, &npleft, &ntot);

      t1 = second();

      if(ntot > 0)
        {
          iter++;

          if(iter > 0 && ThisTask == 0)
            {
              printf("SIDM: ngb iteration %d: need to repeat for %llu particles. (took %g sec)\n", iter, (unsigned long long) ntot, timediff(t0, t1));
              fflush(stdout);
            }

          if(iter > MAXITER)
            {
              printf("SIDM: failed to converge in neighbour iteration in density()\n");
              fflush(stdout);
              terminate("BAD");
            }
        }
    }
  while(ntot > 0);


  myfree(Todo);
  myfree(Right);
  myfree(Left);

  mpi_printf("SIDM: done with HSML values\n");

}


void sidm_findHsml_evaluate(int target, int mode, int threadid)
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

                  out.Density[state] += (in_mass_2 * wk);
                  out.Vx[state] += P[no].Vel[0];
                  out.Vy[state] += P[no].Vel[1];
                  out.Vz[state] += P[no].Vel[2];
                  out.VelDisp[state] += P[no].Vel[0] * P[no].Vel[0] + P[no].Vel[1] * P[no].Vel[1] + P[no].Vel[2] * P[no].Vel[2];

                  phys_rho = in_mass_2 * wk;
                  phys_rho *= SIDM_arho;
                  phys_rel_velx = SIDM_avel * (P[no].Vel[0] - vel[0]);
                  phys_rel_vely = SIDM_avel * (P[no].Vel[1] - vel[1]);
                  phys_rel_velz = SIDM_avel * (P[no].Vel[2] - vel[2]);
                  phys_rel_vel = sqrt(phys_rel_velx * phys_rel_velx + phys_rel_vely * phys_rel_vely + phys_rel_velz * phys_rel_velz);

                  Ekin = SIDM_avel * SIDM_avel * (0.5 * in_mass_1 * (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]) + 0.5 * in_mass_2 * (P[no].Vel[0] * P[no].Vel[0] + P[no].Vel[1] * P[no].Vel[1] + P[no].Vel[2] * P[no].Vel[2]));

                  int cflag = 0;
                  for(reaction = 0; reaction < SIDM_REACTIONS; reaction++)      
                    if((pstate == SMSIDM[reaction].In1) && (state == SMSIDM[reaction].In2))       
                      {
                        out.PSum[reaction] += sidm_scatter_P(phys_rho, phys_rel_vel, Ekin, reaction, &retval);
		        if (retval == -1) 
                          out.EnergyForbidden[reaction]++;
                        cflag++;
                      }
                  if(cflag == 0)
                    terminate("SIDM: scatter states problem %d\n", cflag);

                  out.Ngb++;
                  out.NgbStates[state]++;

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

                  out.Density[state] += (in_mass_2 * wk);
                  out.Vx[state] += Tree_Points[n].Vel[0];
                  out.Vy[state] += Tree_Points[n].Vel[1];
                  out.Vz[state] += Tree_Points[n].Vel[2];
                  out.VelDisp[state] += Tree_Points[n].Vel[0] * Tree_Points[n].Vel[0] + Tree_Points[n].Vel[1] * Tree_Points[n].Vel[1] + Tree_Points[n].Vel[2] * Tree_Points[n].Vel[2];

                  phys_rho = in_mass_2 * wk;
                  phys_rho *= SIDM_arho;
                  phys_rel_velx = SIDM_avel * (Tree_Points[n].Vel[0] - vel[0]);
                  phys_rel_vely = SIDM_avel * (Tree_Points[n].Vel[1] - vel[1]);
                  phys_rel_velz = SIDM_avel * (Tree_Points[n].Vel[2] - vel[2]);
                  phys_rel_vel = sqrt(phys_rel_velx * phys_rel_velx + phys_rel_vely * phys_rel_vely + phys_rel_velz * phys_rel_velz);


                  Ekin = SIDM_avel * SIDM_avel * (0.5 * in_mass_1 * (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]) + 0.5 * in_mass_2 * (Tree_Points[n].Vel[0] * Tree_Points[n].Vel[0] + Tree_Points[n].Vel[1] * Tree_Points[n].Vel[1] + Tree_Points[n].Vel[2] * Tree_Points[n].Vel[2])); 

                  int cflag = 0;
                  for(reaction = 0; reaction < SIDM_REACTIONS; reaction++)      
                    if((pstate == SMSIDM[reaction].In1) && (state == SMSIDM[reaction].In2))       
                      {
                        out.PSum[reaction] += sidm_scatter_P(phys_rho, phys_rel_vel, Ekin, reaction, &retval);
                        if (retval == -1)
                          out.EnergyForbidden[reaction]++;
                        cflag++;
                      }
                  if(cflag == 0)
                    terminate("SIDM: scatter states problem %d\n", cflag);

                  out.Ngb++;
                  out.NgbStates[state]++;


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
