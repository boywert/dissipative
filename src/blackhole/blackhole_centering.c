/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/blackhole/blackhole_centering.c
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

#if defined(BLACK_HOLES) && defined(BH_NEW_CENTERING)

static struct bh_aux_data
{
  MyFloat Ngb;
  MyFloat RhoTot;
  MyFloat Vel[3];
  MyFloat MinPot;
  MyFloat MinPotPos[3];
  MyFloat MinPotGravAccel[3];
#ifdef PMGRID
  MyFloat MinPotGravPM[3];
#endif
  MyFloat Left;
  MyFloat Right;
  MyFloat Dhsmlrho;
}
 *BHPaux;


#ifdef BLACK_HOLES
#define BPPaux(i) BHPaux[P[(i)].AuxDataID]
#endif


 static struct treepoint_auxdata
 {
   MyFloat Potential;
   MyFloat Vel[3];
   MyFloat GravAccel[3];
 #ifdef PMGRID
   MyFloat GravPM[3];
 #endif
 } *TreeAux_Points;


static int blackhole_centering_evaluate(int target, int mode, int threadid);
static void blackhole_centering_provide_treeexported_potential_values(void);


/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble Pos[3];
  MyFloat  Hsml;
  MyFloat  MassOfBH;

  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;


 /* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in * in, int i, int firstnode)
{
  for(int k = 0; k < 3; k++)
    in->Pos[k] = P[i].Pos[k];

  in->Hsml = BPP(i).HsmlCentering;
  in->MassOfBH = P[i].Mass;

  in->Firstnode = firstnode;
}

/* local data structure that holds results acquired on remote processors */
typedef struct
{
  MyFloat Ngb;
  MyFloat RhoTot;
  MyFloat Dhsmlrho;
  MyFloat Vel[3];
  MyFloat MinPot;
  MyFloat MinPotPos[3];
  MyFloat MinPotGravAccel[3];
#ifdef PMGRID
  MyFloat MinPotGravPM[3];
#endif

} data_out;

static data_out *DataResult, *DataOut;



 /* routine to store or combine result data */
static void out2particle(data_out * out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES)      /* initial store */
    {
      BPPaux(i).Ngb = out->Ngb;
      BPPaux(i).RhoTot = out->RhoTot;
      BPPaux(i).Dhsmlrho = out->Dhsmlrho;
      BPPaux(i).MinPot = out->MinPot;
      for(int k = 0; k < 3; k++)
        {
          BPPaux(i).Vel[k] = out->Vel[k];
          BPPaux(i).MinPotPos[k] = out->MinPotPos[k];
          BPPaux(i).MinPotGravAccel[k] = out->MinPotGravAccel[k];
#ifdef PMGRID
          BPPaux(i).MinPotGravPM[k] = out->MinPotGravPM[k];
#endif
        }
    }
  else                          /* merge */
    {
      BPPaux(i).Ngb += out->Ngb;
      BPPaux(i).RhoTot += out->RhoTot;
      BPPaux(i).Dhsmlrho += out->Dhsmlrho;
      for(int k = 0; k < 3; k++)
        BPPaux(i).Vel[k] += out->Vel[k];

      if(out->MinPot < BPPaux(i).MinPot)
        {
          BPPaux(i).MinPot = out->MinPot;
          for(int k = 0; k < 3; k++)
            {
              BPPaux(i).MinPotPos[k] = out->MinPotPos[k];
              BPPaux(i).MinPotGravAccel[k] = out->MinPotGravAccel[k];
#ifdef PMGRID
              BPPaux(i).MinPotGravPM[k] = out->MinPotGravPM[k];
#endif
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
              if(generic_polling_primary(count, TimeBinsBHAccretion.NActiveParticles))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        idx = NextParticle++;

        if(idx >= TimeBinsBHAccretion.NActiveParticles)
          break;

        int i = TimeBinsBHAccretion.ActiveParticleList[idx];
        if(i < 0)
          continue;

        if(blackhole_isactive(i))
          blackhole_centering_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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

        blackhole_centering_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}


void blackhole_centering(void)
{
  if(TimeBinsBHAccretion.GlobalNActiveParticles == 0)
    return;

  TIMER_START(CPU_BH_DENSITY);

  mpi_printf("BH_CENTERING: Begin finding region around BHs, Tree_NumPartImported=%d\n", Tree_NumPartImported);

  BHPaux = mymalloc("BHPaux", All.MaxPartBHs * sizeof(struct bh_aux_data));
  TreeAux_Points = (struct treepoint_auxdata *) mymalloc("Treeaux_Points", Tree_NumPartImported * sizeof(struct treepoint_auxdata));

  blackhole_centering_provide_treeexported_potential_values();

  for(int idx = 0; idx < TimeBinsBHAccretion.NActiveParticles; idx++)
    {
      int i = TimeBinsBHAccretion.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(blackhole_isactive(i))
        {
          BPPaux(i).Left = BPPaux(i).Right = 0;
          if(BPP(i).HsmlCentering == 0)
            BPP(i).HsmlCentering = BPP(i).BH_Hsml;
        }
    }

  generic_set_MaxNexport();

  int iter = 0;
  long long ntot;

  /* we will repeat the whole thing for those black holes where we didn't find the right number of neighbors */
  do
    {
      generic_comm_pattern(TimeBinsBHAccretion.NActiveParticles, kernel_local, kernel_imported);

      /* do final operations on results */
      int npleft = 0;
      for(int idx = 0; idx < TimeBinsBHAccretion.NActiveParticles; idx++)
        {
          int i = TimeBinsBHAccretion.ActiveParticleList[idx];
          if(i < 0)
            continue;

          if(blackhole_isactive(i))
            {
              if(BPPaux(i).RhoTot > 0)
                {
                  BPPaux(i).Vel[0] /= BPPaux(i).RhoTot;
                  BPPaux(i).Vel[1] /= BPPaux(i).RhoTot;
                  BPPaux(i).Vel[2] /= BPPaux(i).RhoTot;
                }

              /* now check whether we had the right number of neighbors */

              if(BPPaux(i).Ngb < 0.95 * All.BlackHoleCenteringMassMultiplier || BPPaux(i).Ngb > 1.05 * All.BlackHoleCenteringMassMultiplier)
                {
                  /* need to redo this particle */
                  npleft++;

                  if(BPPaux(i).Ngb > 0)
                    {
                      BPPaux(i).Dhsmlrho *= BPP(i).HsmlCentering / (NUMDIMS * BPPaux(i).Ngb / (NORM_COEFF * pow(BPP(i).HsmlCentering, 3)));

                      if(BPPaux(i).Dhsmlrho > -0.9)    /* note: this would be -1 if only a single particle at zero lag is found */
                        BPPaux(i).Dhsmlrho = 1 / (1 + BPPaux(i).Dhsmlrho);
                      else
                        BPPaux(i).Dhsmlrho = 1;
                    }
                  else
                    BPPaux(i).Dhsmlrho = 1;

                  if(BPPaux(i).Left > 0 && BPPaux(i).Right > 0)
                    if((BPPaux(i).Right - BPPaux(i).Left) < 1.0e-3 * BPPaux(i).Left && BPPaux(i).Ngb > 0)
                      {
                        /* this one should be ok */
                        npleft--;
                        P[i].TimeBinHydro = -P[i].TimeBinHydro - 1;     /* Mark as inactive */
                        continue;
                      }

                  if(BPPaux(i).Ngb < 0.95 * All.BlackHoleCenteringMassMultiplier)
                    BPPaux(i).Left = dmax(BPP(i).HsmlCentering, BPPaux(i).Left);
                  else
                    {
                      if(BPPaux(i).Right != 0)
                        {
                          if(BPP(i).HsmlCentering < BPPaux(i).Right)
                            BPPaux(i).Right = BPP(i).HsmlCentering;
                        }
                      else
                        BPPaux(i).Right = BPP(i).HsmlCentering;
                    }

                  if(BPPaux(i).Right > 0 && BPPaux(i).Left > 0)
                    BPP(i).HsmlCentering = pow(0.5 * (pow(BPPaux(i).Left, 3) + pow(BPPaux(i).Right, 3)), 1.0 / 3);
                  else
                    {
                      if(BPPaux(i).Right == 0 && BPPaux(i).Left == 0)
                        terminate("BLACK_HOLES: Right[i] = Left[i] = 0 in the computation of HsmlCentering"); /* can't occur */

                      if(BPPaux(i).Right == 0 && BPPaux(i).Left > 0)
                        {
                          double fac = 1.26;

                          if(fabs(BPPaux(i).Ngb - All.BlackHoleCenteringMassMultiplier) < 0.5 * All.BlackHoleCenteringMassMultiplier)
                            {
                              fac = 1 - (BPPaux(i).Ngb - All.BlackHoleCenteringMassMultiplier) / (NUMDIMS * BPPaux(i).Ngb) * BPPaux(i).Dhsmlrho;

                              if(fac > 1.26)
                                fac = 1.26;
                            }

                          BPP(i).HsmlCentering *= fac;
                        }

                      if(BPPaux(i).Right > 0 && BPPaux(i).Left == 0)
                        {
                          double fac = 1 / 1.26;

                          if(fabs(BPPaux(i).Ngb - All.BlackHoleCenteringMassMultiplier) < 0.5 * All.BlackHoleCenteringMassMultiplier)
                            {
                              fac = 1 - (BPPaux(i).Ngb - All.BlackHoleCenteringMassMultiplier) / (NUMDIMS * BPPaux(i).Ngb) * BPPaux(i).Dhsmlrho;

                              if(fac < 1 / 1.26)
                                fac = 1 / 1.26;
                            }
                          BPP(i).HsmlCentering *= fac;
                        }
                    }
                }
              else
                P[i].TimeBinHydro = -P[i].TimeBinHydro - 1;     /* Mark as inactive */
            }
        }

      sumup_large_ints(1, &npleft, &ntot);

      if(ntot > 0)
        {
          iter++;

          if(iter > 0)
            mpi_printf("BH_CENTERING: blackhole ngb iteration %d: need to repeat for %lld particles.\n", iter, ntot);

          if(iter > MAXITER)
            terminate("BH_CENTERING: failed to converge in neighbor iteration in blackhole_centering()\n");
        }
    }
  while(ntot > 0);


  /* mark as active again */
  for(int idx = 0; idx < TimeBinsBHAccretion.NActiveParticles; idx++)
    {
      int i = TimeBinsBHAccretion.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].TimeBinHydro < 0)
        P[i].TimeBinHydro = -P[i].TimeBinHydro - 1;
    }


  /* Now carry out the repositioning */

  for(int idx = 0; idx < TimeBinsBHAccretion.NActiveParticles; idx++)
    {
      int i = TimeBinsBHAccretion.ActiveParticleList[idx];
      if(i < 0)
        continue;

#ifdef PERIODIC
      double xtmp, ytmp, ztmp;
#endif

      double dx = GRAVITY_NEAREST_X(BPPaux(i).MinPotPos[0] - P[i].Pos[0]);
      double dy = GRAVITY_NEAREST_Y(BPPaux(i).MinPotPos[1] - P[i].Pos[1]);
      double dz = GRAVITY_NEAREST_Z(BPPaux(i).MinPotPos[2] - P[i].Pos[2]);

      double dvx = BPPaux(i).Vel[0] - P[i].Vel[0];
      double dvy = BPPaux(i).Vel[1] - P[i].Vel[1];
      double dvz = BPPaux(i).Vel[2] - P[i].Vel[2];

      double r = sqrt(dx * dx + dy * dy + dz * dz);
      double v = sqrt(dvx * dvx + dvy * dvy + dvz * dvz);

      fprintf(FdBlackHolesRepos, "%g %lld %g   %g   %g\n", All.Time, (long long)P[i].ID, BPP(i).BH_Mass, r, All.Time * v);

      for(int k = 0; k < 3; k++)
        {
          P[i].Pos[k] = BPPaux(i).MinPotPos[k];
          P[i].Vel[k] = BPPaux(i).Vel[k];
          P[i].GravAccel[k] = BPPaux(i).MinPotGravAccel[k];
#ifdef PMGRID
          P[i].GravPM[k] = BPPaux(i).MinPotGravPM[k];
#endif
        }
    }

  myfree(TreeAux_Points);
  myfree(BHPaux);

  mpi_printf("BH_CENTERING: Done.\n");

  TIMER_STOP(CPU_BH_DENSITY);
}



#ifdef HIERARCHICAL_GRAVITY
#define INDEX(idx) (TimeBinsGravity.ActiveParticleList[idx])
#else
#define INDEX(idx) (idx)
#endif


static void blackhole_centering_provide_treeexported_potential_values(void)
{
  struct treepoint_auxdata *export_Tree_Points = (struct treepoint_auxdata *) mymalloc("export_Tree_Points", Tree_NumPartExported * sizeof(struct treepoint_auxdata));

  for(int j = 0; j < NTask; j++)
    Force_Send_count[j] = 0;

  for(int idx = 0; idx < NTreeInsert; idx++)        /* prepare particle data to be copied to other tasks */
    {
      int i = INDEX(idx);
      if(i < 0)
        continue;

#ifdef TRACER_PARTICLE
      if(P[i].Type == TRACER_PARTICLE)
        continue;
#endif

      int task = Tree_Task_list[i];

      if(task != ThisTask)
        {
          int n = Force_Send_offset[task] + Force_Send_count[task]++;

          export_Tree_Points[n].Potential = P[i].Potential;
          for(int k = 0; k < 3; k++)
            {
              export_Tree_Points[n].Vel[k] = P[i].Vel[k];
              export_Tree_Points[n].GravAccel[k] = P[i].GravAccel[k];
#ifdef PMGRID
              export_Tree_Points[n].GravPM[k] = P[i].GravPM[k];
#endif
            }

        }
    }

  /* exchange  data */
  for(int ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        if(Force_Send_count[recvTask] > 0 || Force_Recv_count[recvTask] > 0)
          MPI_Sendrecv(&export_Tree_Points[Force_Send_offset[recvTask]], Force_Send_count[recvTask] * sizeof(struct treepoint_auxdata), MPI_BYTE, recvTask, TAG_DENS_A,
                       &TreeAux_Points[Force_Recv_offset[recvTask]], Force_Recv_count[recvTask] * sizeof(struct treepoint_auxdata), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

  myfree(export_Tree_Points);
}




static int blackhole_centering_evaluate(int target, int mode, int threadid)
{
  int numnodes, *firstnode;
#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
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

  MyDouble *pos = target_data->Pos;
  double h = target_data->Hsml;
  double h2 = h * h;
  double hinv = 1.0 / h;
  double hinv3 = hinv * hinv * hinv;
  double hinv4 = hinv * hinv3;

  double vel[3] = {0, 0, 0};
  double rhotot = 0;
  double weighted_numngb = 0;
  double dhsmlrho = 0;

  MyFloat minpotpos[3] = {0, 0, 0 };
  MyFloat minpot = BHPOTVALUEINIT;
  MyFloat gravaccel[3] = {0, 0, 0};
#ifdef PMGRID
  MyFloat gravpm[3] = {0, 0, 0};
#endif

  for(int k = 0; k < numnodes; k++)
    {
      int no;

      if(mode == MODE_LOCAL_PARTICLES)
        {
          no = Tree_MaxPart; /* root node */
        }
      else
        {
          no = firstnode[k];
          no = Nodes[no].u.d.nextnode; /* open it */
        }

      while(no >= 0)
        {
          if(no < Tree_MaxPart) /* single particle */
            {
              double dx = GRAVITY_NEAREST_X(Tree_Pos_list[3 * no + 0] - pos[0]);
              double dy = GRAVITY_NEAREST_Y(Tree_Pos_list[3 * no + 1] - pos[1]);
              double dz = GRAVITY_NEAREST_Z(Tree_Pos_list[3 * no + 2] - pos[2]);

              double r2 = dx * dx + dy * dy + dz * dz;

              if(r2 < h2)
                {
                  /* minpot */
                  if(P[no].Potential < minpot)
                    {
                      minpot = P[no].Potential;

                      for(int k = 0; k < 3; k++)
                        {
                          minpotpos[k] = Tree_Pos_list[3 * no + k];
                          gravaccel[k] = P[no].GravAccel[k];
#ifdef PMGRID
                          gravpm[k] = P[no].GravPM[k];
#endif
                        }
                    }

                  double r = sqrt(r2);
                  double u = r * hinv;

                  double wk, dwk;
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

                  rhotot += P[no].Mass * wk;
                  vel[0] += P[no].Mass * wk * P[no].Vel[0];
                  vel[1] += P[no].Mass * wk * P[no].Vel[1];
                  vel[2] += P[no].Mass * wk * P[no].Vel[2];

                  double ngb_eff = P[no].Mass / target_data->MassOfBH;

                  weighted_numngb += (NORM_COEFF * ngb_eff * wk / hinv3); /* 4.0/3 * PI = 4.188790204786 */
                  dhsmlrho += (-ngb_eff * (NUMDIMS * hinv * wk + u * dwk));
                }
              no = Nextnode[no];
            }
          else if(no < Tree_MaxPart + Tree_MaxNodes) /* internal node */
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                {
                  if(no < Tree_FirstNonTopLevelNode) /* we reached a top-level node again, which means that we are done with the branch */
                    break;
                }

              struct NODE *current = &Nodes[no];

              no = current->u.d.sibling; /* in case the node can be discarded */

              double dist = h + 0.5 * current->len;
              double dx = NGB_PERIODIC_LONG_X(current->center[0] - pos[0]);
              if(dx > dist)
                continue;
              double dy = NGB_PERIODIC_LONG_Y(current->center[1] - pos[1]);
              if(dy > dist)
                continue;
              double dz = NGB_PERIODIC_LONG_Z(current->center[2] - pos[2]);
              if(dz > dist)
                continue;
              /* now test against the minimal sphere enclosing everything */
              dist += FACT1 * current->len;
              if(dx * dx + dy * dy + dz * dz > dist * dist)
                continue;

              no = current->u.d.nextnode; /* ok, we need to open the node */
            }
          else if(no >= Tree_ImportedNodeOffset) /* point from imported nodelist */
            {
              int n = no - Tree_ImportedNodeOffset;

              double dx = GRAVITY_NEAREST_X(Tree_Points[n].Pos[0] - pos[0]);
              double dy = GRAVITY_NEAREST_Y(Tree_Points[n].Pos[1] - pos[1]);
              double dz = GRAVITY_NEAREST_Z(Tree_Points[n].Pos[2] - pos[2]);

              double r2 = dx * dx + dy * dy + dz * dz;

              if(r2 < h2)
                {
                  /* minpot */
                  if(TreeAux_Points[n].Potential < minpot)
                    {
                      minpot = TreeAux_Points[n].Potential;

                      for(int k = 0; k < 3; k++)
                        {
                          minpotpos[k] = Tree_Points[n].Pos[k];
                          gravaccel[k] = TreeAux_Points[n].GravAccel[k];
#ifdef PMGRID
                          gravpm[k] = TreeAux_Points[n].GravPM[k];
#endif
                        }
                    }

                  double r = sqrt(r2);
                  double u = r * hinv;

                  double wk, dwk;
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

                  rhotot += Tree_Points[n].Mass * wk;
                  vel[0] += Tree_Points[n].Mass * wk * TreeAux_Points[n].Vel[0];
                  vel[1] += Tree_Points[n].Mass * wk * TreeAux_Points[n].Vel[1];
                  vel[2] += Tree_Points[n].Mass * wk * TreeAux_Points[n].Vel[2];

                  double ngb_eff = Tree_Points[n].Mass / target_data->MassOfBH;

                  weighted_numngb += (NORM_COEFF * ngb_eff * wk / hinv3); /* 4.0/3 * PI = 4.188790204786 */
                  dhsmlrho += (-ngb_eff * (NUMDIMS * hinv * wk + u * dwk));
                }
              no = Nextnode[no - Tree_MaxNodes];
            }
          else /* pseudo particle */
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

  out.Ngb = weighted_numngb;
  out.RhoTot = rhotot;
  out.Dhsmlrho = dhsmlrho;
  out.MinPot = minpot;

  for(int k = 0; k < 3; k++)
    {
      out.MinPotPos[k] = minpotpos[k];
      out.Vel[k] = vel[k];
      out.MinPotGravAccel[k] = gravaccel[k];
#ifdef PMGRID
      out.MinPotGravPM[k] = gravpm[k];
#endif
    }

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#endif
