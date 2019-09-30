/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/subfind/subfind_density.c
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
#include <sys/stat.h>
#include <sys/types.h>

#include "../allvars.h"
#include "../proto.h"

#ifdef SUBFIND

#include "../fof/fof.h"
#include "subfind.h"


static char *Todo;
static int *DM_NumNgb;
#ifdef SUBFIND_CALC_MORE
static MyFloat *Vx, *Vy, *Vz;
#endif

static int subfind_density_evaluate(int target, int mode, int threadid);


/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble Pos[3];
  MyFloat Hsml;

  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;


 /* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in * in, int i, int firstnode)
{
#ifdef CELL_CENTER_GRAVITY
  if(P[i].Type == 0)
    {
      in->Pos[0] = SphP[i].Center[0];
      in->Pos[1] = SphP[i].Center[1];
      in->Pos[2] = SphP[i].Center[2];
    }
  else
#endif
    {
      in->Pos[0] = P[i].Pos[0];
      in->Pos[1] = P[i].Pos[1];
      in->Pos[2] = P[i].Pos[2];
    }
  in->Hsml = PS[i].Hsml;

  in->Firstnode = firstnode;
}

/* local data structure that holds results acquired on remote processors */
typedef struct
{
  int Ngb;
  MyFloat Rho;
#ifdef SUBFIND_CALC_MORE
  MyFloat VelDisp, Vx, Vy, Vz, RhoDM;
#endif
} data_out;

static data_out *DataResult, *DataOut;




 /* routine to store or combine result data */
static void out2particle(data_out * out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES)      /* initial store */
    {
      DM_NumNgb[i] = out->Ngb;
      PS[i].Density = out->Rho;
#ifdef SUBFIND_CALC_MORE
      Vx[i] = out->Vx;
      Vy[i] = out->Vy;
      Vz[i] = out->Vz;
      PS[i].SubfindVelDisp = out->VelDisp;
      PS[i].SubfindDMDensity = out->RhoDM;
#endif
    }
  else                          /* combine */
    {
      DM_NumNgb[i] += out->Ngb;
      PS[i].Density += out->Rho;
#ifdef SUBFIND_CALC_MORE
      Vx[i] += out->Vx;
      Vy[i] += out->Vy;
      Vz[i] += out->Vz;
      PS[i].SubfindVelDisp += out->VelDisp;
      PS[i].SubfindDMDensity += out->RhoDM;
#endif
    }
}


#include "../generic_comm_helpers2.h"



static void kernel_local(void)
{
  int i;
#ifdef GENERIC_ASYNC
  int flag = 0;
#endif

#pragma omp parallel private(i)
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
              if(generic_polling_primary(count, NumPart))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        i = NextParticle++;

        if(i >= NumPart)
          break;

        if(Todo[i])
          subfind_density_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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

        subfind_density_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

double subfind_density(int mode)
{
  long long ntot;
  int i, npleft, iter;
  MyFloat *Left, *Right;
  double t0, t1, tstart, tend;

  if(mode == FIND_SMOOTHING_LENGTHS)
    mpi_printf("SUBFIND: finding smoothing length for all particles\n");
  else
    mpi_printf("SUBFIND: finding total densities around all particles\n");

  tstart = second();


  int HsmlFlag = 0;

#ifdef SUBFIND_CALC_MORE
  HsmlFlag = 1;                 /* in this case, calculate densities for all particles, not only those in groups */
#endif


  DM_NumNgb = (int *) mymalloc_movable(&DM_NumNgb, "DM_NumNgb", sizeof(int) * NumPart);
  Left = (MyFloat *) mymalloc_movable(&Left, "Left", sizeof(MyFloat) * NumPart);
  Right = (MyFloat *) mymalloc_movable(&Right, "Right", sizeof(MyFloat) * NumPart);
  Todo = (char *) mymalloc_movable(&Todo, "Todo", sizeof(char) * NumPart);

#ifdef SUBFIND_CALC_MORE
  Vx = (MyFloat *) mymalloc("Vx", sizeof(MyFloat) * NumPart);
  Vy = (MyFloat *) mymalloc("Vy", sizeof(MyFloat) * NumPart);
  Vz = (MyFloat *) mymalloc("Vz", sizeof(MyFloat) * NumPart);
#endif


  generic_set_MaxNexport();

  for(i = 0; i < NumPart; i++)
    {
      Left[i] = Right[i] = 0;
      DM_NumNgb[i] = 0;
      Todo[i] = 1;
      if((PS[i].GrNr >= TotNgroups) && (HsmlFlag == 0)) // particle not in groups
        Todo[i] = 0;

#ifdef REFINEMENT_HIGH_RES_GAS
      if((PS[i].GrNr >= TotNgroups) && (P[i].Type == 4 || P[i].Type == 5))      // particle a star or BH but not in group
        Todo[i] = 0;

      if(P[i].Type != 0 && P[i].Type != 1 && P[i].Type != 4 && P[i].Type != 5)
        Todo[i] = 0;
      if(P[i].Type == 0)
        if(SphP[i].AllowRefinement == 0)
          Todo[i] = 0;
#endif

#ifdef TRACER_PARTICLE
      if(P[i].Type == TRACER_PARTICLE)  /* tracer particle */
        Todo[i] = 0;
#endif

      PS[i].Density = 0;
#ifdef SUBFIND_CALC_MORE
      PS[i].SubfindHsml = 0;
      PS[i].SubfindDensity = 0;
      PS[i].SubfindDMDensity = 0;
      PS[i].SubfindVelDisp = 0;
#endif
    }

  iter = 0;

  /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
  do
    {
      t0 = second();

      generic_comm_pattern(NumPart, kernel_local, kernel_imported);


      /* do final operations on results */
      for(i = 0, npleft = 0; i < NumPart; i++)
        {
          /* now check whether we had enough neighbours */

          if(Todo[i] && mode == FIND_SMOOTHING_LENGTHS)
            {
              if(abs(DM_NumNgb[i] - All.DesNumNgb) > All.MaxNumNgbDeviation && ((Right[i] - Left[i]) > 1.0e-4 * Left[i] || Left[i] == 0 || Right[i] == 0))
                {
                  /* need to redo this particle */
                  npleft++;

                  if(DM_NumNgb[i] < All.DesNumNgb)
                    Left[i] = (MyFloat) dmax(PS[i].Hsml, Left[i]);
                  else
                    {
                      if(Right[i] != 0)
                        {
                          if(PS[i].Hsml < Right[i])
                            Right[i] = PS[i].Hsml;
                        }
                      else
                        Right[i] = PS[i].Hsml;
                    }

                  if(iter >= MAXITER - 10)
                    {
                      printf("SUBFIND: i=%d task=%d ID=%d Hsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g\n   pos=(%g|%g|%g)\n", i, ThisTask,
                             (int) P[i].ID, PS[i].Hsml, Left[i], Right[i], (double) DM_NumNgb[i], Right[i] - Left[i], P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
                      myflush(stdout);
                    }

                  if(Right[i] > 0 && Left[i] > 0)
                    PS[i].Hsml = (MyFloat) pow(0.5 * (pow(Left[i], 3) + pow(Right[i], 3)), 1.0 / 3);
                  else
                    {
                      if(Right[i] == 0 && Left[i] == 0)
                        terminate("can't occur");

                      if(Right[i] == 0 && Left[i] > 0)
                        PS[i].Hsml *= 1.26;

                      if(Right[i] > 0 && Left[i] == 0)
                        PS[i].Hsml /= 1.26;
                    }
                }
              else
                Todo[i] = 0;
            }
        }

      sumup_large_ints(1, &npleft, &ntot);

      t1 = second();

      if(ntot > 0 && mode == FIND_SMOOTHING_LENGTHS)
        {
          iter++;

          if(iter > 0)
            mpi_printf("SUBFIND: ngb iteration %2d: need to repeat for %15lld particles. (took %g sec)\n", iter, ntot, timediff(t0, t1));

          if(iter > MAXITER)
            terminate("failed to converge in neighbour iteration in density()\n");
        }
    }
  while(ntot > 0);

#ifdef SUBFIND_CALC_MORE
  double vel_to_phys;

  vel_to_phys = 1.0 / All.cf_atime;

  for(i = 0; i < NumPart; i++)
    {
      Vx[i] /= DM_NumNgb[i];
      Vy[i] /= DM_NumNgb[i];
      Vz[i] /= DM_NumNgb[i];
      PS[i].SubfindVelDisp /= DM_NumNgb[i];
      PS[i].SubfindVelDisp = vel_to_phys * sqrt(PS[i].SubfindVelDisp - Vx[i] * Vx[i] - Vy[i] * Vy[i] - Vz[i] * Vz[i]);
    }
#endif

#ifdef SUBFIND_CALC_MORE
  myfree_movable(Vz);
  myfree_movable(Vy);
  myfree_movable(Vx);
#endif
  myfree_movable(Todo);
  myfree_movable(Right);
  myfree_movable(Left);
  myfree_movable(DM_NumNgb);


#ifdef SUBFIND_CALC_MORE
  for(i = 0; i < NumPart; i++)
    {
      PS[i].SubfindHsml = PS[i].Hsml;
      PS[i].SubfindDensity = PS[i].Density;
    }
#endif

  tend = second();
  return timediff(tstart, tend);
}



/*! This function represents the core of the SPH density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */

static int subfind_density_evaluate(int target, int mode, int threadid)
{
  int k, numnodes, *firstnode, type;
  double hsml;
  double rhosum = 0;
  MyDouble *pos;
  int numngb = 0, no, p;
  struct NODE *current;
  double dx, dy, dz, r2, mass;
  double h2, hinv, hinv3, r, u, wk;
#ifdef PERIODIC
  MyDouble xtmp, ytmp, ztmp;
#endif
#ifdef SUBFIND_CALC_MORE
  double vxsum = 0, vysum = 0, vzsum = 0, v2sum = 0, rhodmsum = 0;
#endif

  data_in local, *target_data;
  data_out out;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      target_data = &local;

      numnodes = 1;
      firstnode = NULL;
    }
  else
    {
      target_data = &DataGet[target];

      generic_get_numnodes(target, &numnodes, &firstnode);
    }

  pos = target_data->Pos;
  hsml = target_data->Hsml;

  h2 = hsml * hsml;
  hinv = 1.0 / hsml;
  hinv3 = hinv * hinv * hinv;

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
              p = no;
              no = Nextnode[no];

              dx = FOF_NEAREST_LONG_X(Tree_Pos_list[3 * p + 0] - pos[0]);
              if(dx > hsml)
                continue;
              dy = FOF_NEAREST_LONG_Y(Tree_Pos_list[3 * p + 1] - pos[1]);
              if(dy > hsml)
                continue;
              dz = FOF_NEAREST_LONG_Z(Tree_Pos_list[3 * p + 2] - pos[2]);
              if(dz > hsml)
                continue;

              if((r2 = (dx * dx + dy * dy + dz * dz)) > hsml * hsml)
                continue;

              mass = P[p].Mass;
              type = P[p].Type;
            }
          else if(no < Tree_MaxPart + Tree_MaxNodes)    /* internal node */
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                {
                  if(no < Tree_FirstNonTopLevelNode)    /* we reached a top-level node again, which means that we are done with the branch */
                    break;
                }

              current = &Nodes[no];

              no = current->u.d.sibling;        /* in case the node can be discarded */

              double dist = hsml + 0.5 * current->len;

              dx = (MyFloat) FOF_NEAREST_LONG_X(current->center[0] - pos[0]);
              if(dx > dist)
                continue;
              dy = (MyFloat) FOF_NEAREST_LONG_Y(current->center[1] - pos[1]);
              if(dy > dist)
                continue;
              dz = (MyFloat) FOF_NEAREST_LONG_Z(current->center[2] - pos[2]);
              if(dz > dist)
                continue;
              /* now test against the minimal sphere enclosing everything */
              dist += FACT1 * current->len;
              if(dx * dx + dy * dy + dz * dz > dist * dist)
                continue;

              no = current->u.d.nextnode;       /* ok, we need to open the node */
              continue;
            }
          else if(no >= Tree_ImportedNodeOffset)        /* point from imported nodelist */
            {
              int n = no - Tree_ImportedNodeOffset;
              no = Nextnode[no - Tree_MaxNodes];

              dx = FOF_NEAREST_LONG_X(Tree_Points[n].Pos[0] - pos[0]);
              if(dx > hsml)
                continue;
              dy = FOF_NEAREST_LONG_Y(Tree_Points[n].Pos[1] - pos[1]);
              if(dy > hsml)
                continue;
              dz = FOF_NEAREST_LONG_Z(Tree_Points[n].Pos[2] - pos[2]);
              if(dz > hsml)
                continue;

              if((r2 = (dx * dx + dy * dy + dz * dz)) > hsml * hsml)
                continue;

              mass = Tree_Points[n].Mass;
              type = Tree_Points[n].Type;

              p = -1;
            }
          else                  /* pseudo particle */
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                terminate("can't be");

              if(target >= 0)   /* if no target is given, export will not occur */
                tree_treefind_export_node_threads(no, target, threadid);

              no = Nextnode[no - Tree_MaxNodes];
              continue;
            }



          if((1 << type) & (FOF_PRIMARY_LINK_TYPES))
            {
              numngb++;

#ifdef SUBFIND_CALC_MORE
              if(p < 0)
                terminate("this should not occur");

              vxsum += P[p].Vel[0];
              vysum += P[p].Vel[1];
              vzsum += P[p].Vel[2];
              v2sum += P[p].Vel[0] * P[p].Vel[0] + P[p].Vel[1] * P[p].Vel[1] + P[p].Vel[2] * P[p].Vel[2];
#endif
            }

          if(((1 << type) & (FOF_PRIMARY_LINK_TYPES)) || ((1 << type) & (FOF_SECONDARY_LINK_TYPES)))
            if(r2 < h2)
              {
                r = sqrt(r2);

                u = r * hinv;

                if(u < 0.5)
                  wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
                else
                  wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

                rhosum += mass * wk;

#ifdef SUBFIND_CALC_MORE
                if((1 << type) & (FOF_PRIMARY_LINK_TYPES))
                  rhodmsum += mass * wk;
#endif
              }
        }
    }

  out.Ngb = numngb;
  out.Rho = rhosum;
#ifdef SUBFIND_CALC_MORE
  out.Vx = vxsum;
  out.Vy = vysum;
  out.Vz = vzsum;
  out.VelDisp = v2sum;
  out.RhoDM = rhodmsum;
#endif

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}




void subfind_density_hsml_guess(void)   /* set the initial guess for the smoothing length */
{
  int i;
  double hsml_prev = 0;

  for(i = 0; i < NumPart; i++)
    {
#ifdef TRACER_PARTICLE
      if(P[i].Type == TRACER_PARTICLE)
        {
          PS[i].Density = -1.0; // exclude tracers from subfind_serial
          continue;
        }
#endif

      int no, p;

      if((1 << P[i].Type) & (FOF_PRIMARY_LINK_TYPES))
        {
          no = Father[i];

          while(8 * All.DesNumNgb * P[i].Mass > Nodes[no].u.d.mass && Nodes[no].len == 0)
            {
              p = Nodes[no].u.d.father;

              if(p < 0)
                break;

              no = p;
            }

          PS[i].Hsml = hsml_prev = (pow(3.0 / (4 * M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 3) * Nodes[no].len);

          if(PS[i].Hsml == 0)
            {
              printf("Hsml=0 task=%d i=%d no=%d Nodes[no].len=%g Nodes[no].u.d.mass=%g P[i].Mass=%g type=%d ID=%llu  pos=(%g|%g|%g)\n", ThisTask, i,
                     no, Nodes[no].len, Nodes[no].u.d.mass, P[i].Mass, P[i].Type, (long long) P[i].ID, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
              terminate("zero hsml guess\n");
            }
        }
      else
        {
          if(hsml_prev)
            PS[i].Hsml = hsml_prev;
          else
            PS[i].Hsml = All.SofteningTable[P[i].SofteningType];
        }
    }
}

#endif
