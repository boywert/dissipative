/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/blackhole/blackhole_mergers.c
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

#ifdef BLACK_HOLES

static int blackhole_tree_merger_evaluate(int i, int mode, int threadid);

struct bh_properties
{
  MyDouble Mass;
  MyFloat BHMass;
  MyFloat BH_CumMass_QM;
  MyFloat BH_CumEgy_QM;
  MyFloat BH_CumMass_RM;
  MyFloat BH_CumEgy_RM;
  MyFloat Momentum[3];
#ifdef BH_RECOIL_KICK
  MyFloat Pos[3];
  MyFloat Vel[3];
  MyFloat Spin[3];
#endif
#ifdef DRAINGAS
  MyDouble DrainBucketMass;
#endif
#ifdef BH_NF_RADIO
  MyFloat BH_Mdot_radio;
  MyFloat BH_Mdot_quasar;
  MyFloat BH_RadioEgyFeedback;
  MyFloat BH_HaloVvir;
  MyFloat BH_XrayLum;
  MyFloat BH_RadioLum;
#endif
#ifdef BH_BUBBLES
  MyFloat BHMass_bubbles;
  MyFloat BHMass_ini;
#endif
  int BH_CountProgs;
};


static struct bhresultsimported_data
{
  struct bh_properties Prop;
  int index;
}
 *Tree_BhngbResultsImported;

static struct bhdeleted_data
{
  int index;
#ifdef TRACER_MC
  int rindex;                   /* P_index of BH that does the swallowing (recepient of victim's tracers) */
  int rtask;
  MyIDType rID;
#endif
}
 *Tree_BhngbDeleted;

static int N_BH_swallowed_local, N_BH_swallowed_imported;
static char *treeBHexportflag;



static struct bh_properties *BH_accreted;



#ifdef TRACER_MC
static int *treeBHexportRtask;   /**< task for target parent as above */
static int *treeBHexportRindex; /**< PID for target parent when local BH swallows remote tree BH */
static MyIDType *treeBHexportRID; /**< PID for target parent when local BH swallows remote tree BH */
#endif



/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble Pos[3];
  MyDouble Mass;
  MyFloat BH_Hsml;
  MyIDType ID;
#ifdef TRACER_MC
  int rtask;                    /* task which holds the swallowing BH */
  int rindex;                   /* P_index of the BH on rtask */
#endif

  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;


 /* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in * in, int i, int firstnode)
{
  int k;

  if(i < NumPart)
    {
      for(k = 0; k < 3; k++)
        in->Pos[k] = P[i].Pos[k];

      in->Mass = P[i].Mass;
      in->BH_Hsml = BPP(i).BH_Hsml;
      in->ID = P[i].ID;
#ifdef TRACER_MC
      in->rtask = ThisTask;
      in->rindex = i;
#endif
    }
  else
    {
      i -= Tree_ImportedNodeOffset;

      for(k = 0; k < 3; k++)
        in->Pos[k] = Tree_Points[i].Pos[k];

      in->Mass = Tree_Points[i].Mass;
      in->BH_Hsml = TBPP(i).BH_Hsml;
      in->ID = TBPP(i).ID;
#ifdef TRACER_MC
      in->rtask = Tree_Points[i].origTask;
      in->rindex = Tree_Points[i].index;
#endif
    }

  in->Firstnode = firstnode;
}


 /* local data structure that holds results acquired on remote processors */
typedef struct
{
  struct bh_properties Prop;
} data_out;

static data_out *DataResult, *DataOut;


 /* routine to store or combine result data */
static void out2particle(data_out * out, int i, int mode)
{
  int k;

  if(mode == MODE_LOCAL_PARTICLES)      /* initial store */
    {
      if(i < NumPart)
        {
          if(P[i].AuxDataID >= NumBHs)
            terminate("BLACK_HOLES: P[target(=%d)].AuxDataID(=%lld) >= NumBHs=%d", i, (long long) P[i].AuxDataID, NumBHs);

          BH_accreted[P[i].AuxDataID] = out->Prop;
        }
      else
        {
          int idx = Tree_ResultIndexList[i - Tree_ImportedNodeOffset];

          Tree_BhngbResultsImported[idx].Prop = out->Prop;
        }
    }
  else                          /* combine */
    {
      if(i < NumPart)
        {
          if(P[i].AuxDataID >= NumBHs)
            terminate("BLACK_HOLES: P[i(=%d)].AuxDataID(=%lld) >= NumBHs=%d", i, (long long) P[i].AuxDataID, NumBHs);

          BH_accreted[P[i].AuxDataID].Mass += out->Prop.Mass;
          BH_accreted[P[i].AuxDataID].BHMass += out->Prop.BHMass;
          BH_accreted[P[i].AuxDataID].BH_CumMass_QM += out->Prop.BH_CumMass_QM;
          BH_accreted[P[i].AuxDataID].BH_CumEgy_QM += out->Prop.BH_CumEgy_QM;
          BH_accreted[P[i].AuxDataID].BH_CumMass_RM += out->Prop.BH_CumMass_RM;
          BH_accreted[P[i].AuxDataID].BH_CumEgy_RM += out->Prop.BH_CumEgy_RM;

          for(k = 0; k < 3; k++)
            BH_accreted[P[i].AuxDataID].Momentum[k] += out->Prop.Momentum[k];

#ifdef BH_RECOIL_KICK
          for(k = 0; k < 3; k++)
          {
            BH_accreted[P[i].AuxDataID].Pos[k] += out->Prop.Pos[k];
            BH_accreted[P[i].AuxDataID].Vel[k] += out->Prop.Vel[k];
#ifdef BH_SPIN_EVOLUTION
            BH_accreted[P[i].AuxDataID].Spin[k] += out->Prop.Spin[k];
#endif
          }
#endif

#ifdef DRAINGAS
          BH_accreted[P[i].AuxDataID].DrainBucketMass += out->Prop.DrainBucketMass;
#endif
#ifdef BH_BUBBLES
          BH_accreted[P[i].AuxDataID].BHMass_bubbles += out->Prop.BHMass_bubbles;
          BH_accreted[P[i].AuxDataID].BHMass_ini += out->Prop.BHMass_ini;
#endif
#ifdef BH_NF_RADIO
          BH_accreted[P[i].AuxDataID].BH_RadioEgyFeedback += out->Prop.BH_RadioEgyFeedback;
          BH_accreted[P[i].AuxDataID].BH_Mdot_radio = dmax(BH_accreted[P[i].AuxDataID].BH_Mdot_radio, out->Prop.BH_Mdot_radio);
          BH_accreted[P[i].AuxDataID].BH_Mdot_quasar = dmax(BH_accreted[P[i].AuxDataID].BH_Mdot_quasar, out->Prop.BH_Mdot_quasar);
          BH_accreted[P[i].AuxDataID].BH_HaloVvir = dmax(BH_accreted[P[i].AuxDataID].BH_HaloVvir, out->Prop.BH_HaloVvir);
          BH_accreted[P[i].AuxDataID].BH_XrayLum = dmax(BH_accreted[P[i].AuxDataID].BH_XrayLum, out->Prop.BH_XrayLum);
          BH_accreted[P[i].AuxDataID].BH_RadioLum = dmax(BH_accreted[P[i].AuxDataID].BH_RadioLum, out->Prop.BH_RadioLum);
#endif
          BH_accreted[P[i].AuxDataID].BH_CountProgs += out->Prop.BH_CountProgs;
        }
      else
        {
          int idx = Tree_ResultIndexList[i - Tree_ImportedNodeOffset];

          Tree_BhngbResultsImported[idx].Prop.Mass += out->Prop.Mass;
          Tree_BhngbResultsImported[idx].Prop.BHMass += out->Prop.BHMass;
          Tree_BhngbResultsImported[idx].Prop.BH_CumMass_QM += out->Prop.BH_CumMass_QM;
          Tree_BhngbResultsImported[idx].Prop.BH_CumEgy_QM += out->Prop.BH_CumEgy_QM;
          Tree_BhngbResultsImported[idx].Prop.BH_CumMass_RM += out->Prop.BH_CumMass_RM;
          Tree_BhngbResultsImported[idx].Prop.BH_CumEgy_RM += out->Prop.BH_CumEgy_RM;

          for(k = 0; k < 3; k++)
            Tree_BhngbResultsImported[idx].Prop.Momentum[k] += out->Prop.Momentum[k];
#ifdef BH_RECOIL_KICK
          for(k = 0; k < 3; k++)
          {
        	Tree_BhngbResultsImported[idx].Prop.Pos[k] += out->Prop.Pos[k];
        	Tree_BhngbResultsImported[idx].Prop.Vel[k] += out->Prop.Vel[k];
#ifdef BH_SPIN_EVOLUTION
        	Tree_BhngbResultsImported[idx].Prop.Spin[k] += out->Prop.Spin[k];
#endif
          }
#endif
#ifdef DRAINGAS
          Tree_BhngbResultsImported[idx].Prop.DrainBucketMass += out->Prop.DrainBucketMass;
#endif
#ifdef BH_BUBBLES
          Tree_BhngbResultsImported[idx].Prop.BHMass_bubbles += out->Prop.BHMass_bubbles;
          Tree_BhngbResultsImported[idx].Prop.BHMass_ini += out->Prop.BHMass_ini;
#endif
#ifdef BH_NF_RADIO
          Tree_BhngbResultsImported[idx].Prop.BH_RadioEgyFeedback += out->Prop.BH_RadioEgyFeedback;
          Tree_BhngbResultsImported[idx].Prop.BH_HaloVvir = dmax(Tree_BhngbResultsImported[idx].Prop.BH_HaloVvir, out->Prop.BH_HaloVvir);
          Tree_BhngbResultsImported[idx].Prop.BH_Mdot_radio = dmax(Tree_BhngbResultsImported[idx].Prop.BH_Mdot_radio, out->Prop.BH_Mdot_radio);
          Tree_BhngbResultsImported[idx].Prop.BH_Mdot_quasar = dmax(Tree_BhngbResultsImported[idx].Prop.BH_Mdot_quasar, out->Prop.BH_Mdot_quasar);
          Tree_BhngbResultsImported[idx].Prop.BH_XrayLum = dmax(Tree_BhngbResultsImported[idx].Prop.BH_XrayLum, out->Prop.BH_XrayLum);
          Tree_BhngbResultsImported[idx].Prop.BH_RadioLum = dmax(Tree_BhngbResultsImported[idx].Prop.BH_RadioLum, out->Prop.BH_RadioLum);
#endif
          Tree_BhngbResultsImported[idx].Prop.BH_CountProgs += out->Prop.BH_CountProgs;
        }
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
              if(generic_polling_primary(count, Nforces))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        i = NextParticle++;

        if(i >= Nforces)
          break;

        int idx = TargetList[i];

        blackhole_tree_merger_evaluate(idx, MODE_LOCAL_PARTICLES, threadid);
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

        blackhole_tree_merger_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

#ifdef BH_RECOIL_KICK


/* Routine to compute cross product of two vectors */
void cross_product(MyFloat * A, MyFloat * B, MyFloat * C)
{
  C[0] = A[1]*B[2] - A[2]*B[1];
  C[1] = A[2]*B[0] - A[0]*B[2];
  C[2] = A[0]*B[1] - A[1]*B[0];
}

/* Routine to apply recoil kick to a BH remnant just after merger. It assumes that m1<m2.
 * It can compute two contributions for the recoil kick, e.g. due to mass asymmetry and due to spinning
 * BHs. The latter requires the flag BH_SPIN_EVOLUTION to be enabled.
 */
void apply_recoil_kick(MyDouble Mass1, MyDouble Mass2, MyFloat * Pos1, MyFloat * Pos2, MyFloat * Vel1, MyFloat * Vel2, MyFloat * Spin1, MyFloat * Spin2, MyFloat * RecoilVel, MyFloat * FinalSpin)
{
  int k;
  double Q, Eta;
  double r[3] = {Pos1[0]-Pos2[0], Pos1[1]-Pos2[1], Pos1[2]-Pos2[2]};
  double v[3] = {Vel1[0]-Vel2[0], Vel1[1]-Vel2[1], Vel1[2]-Vel2[2]};
  double Rmag = 0, Vmag = 0, Lmag = 0;
  double Vm, Vper, Vpar;
  double eccentricity;
  double e1[3], e2[3], ez[3];
  double cosmofac = All.cf_atime;
  double Spin1z = 0, Spin1xy = 0, Spin2z = 0, Spin2xy = 0, Spin12xyMag = 0;
  double Theta = 0;

  /* Parameters for recoil velocity*/
  double A = 1.2e4; 	// [km/s]	Gonzalez, et al. 2007
  double B = -0.93; 	// [-]		Gonzalez, et al. 2007
  double H = 7.3e3;  	// [km/s]	Campanelli, et al. 2007
  double Xi = M_PI_2;	// [-]		Sijacki, et al. 2009
  double K = 6.0e4;		// [km/s]	Campanelli, et al. 2007
  double Theta0 = 0.184;// [-]		Campanelli, et al. 2007

  Q = Mass1/Mass2;
  Eta = Q/((1 + Q)*(1 + Q));

  /* Building coordinate system */
  cross_product(r, v, ez);
  for(k = 0; k<3; k++)
  {
	  Rmag += r[k]*r[k];
	  Vmag += v[k]*v[k];
	  Lmag += ez[k]*ez[k];
  }
  Rmag = sqrt(Rmag);
  Vmag = sqrt(Vmag);
  Lmag = sqrt(Lmag);
  for(k = 0; k<3; k++)
  {
    e1[k] = r[k]/Rmag;
    ez[k] = ez[k]/Lmag;
  }
  cross_product(ez, e1, e2);

  /* Contribution due to mass asymmetry */
  Vm = A*Eta*Eta*sqrt(1 - 4*Eta)*(1 + B*Eta);

  for(k = 0; k<3; k++)
    RecoilVel[k] = Vm*e1[k];

#ifdef BH_SPIN_EVOLUTION
  /* Recoil kick on the orbital plane, using spin components along angular momentum */
  for(k = 0; k<3; k++)
  {
    Spin1z += Spin1[k]*ez[k];
    Spin2z += Spin2[k]*ez[k];
  }
  Vper = H*Eta*Eta/(1+Q)*(Spin2z - Q*Spin1z);

  /* Recoil kick along angular momentum, using spin components on the orbital plane */
  for(k = 0; k<3; k++)
  {
    Spin1xy += (Spin1[k] - Spin1z*ez[k])*(Spin1[k] - Spin1z*ez[k]);
    Spin2xy += (Spin2[k] - Spin2z*ez[k])*(Spin2[k] - Spin2z*ez[k]);
    Spin12xyMag += (Spin1[k] - Spin2[k] - (Spin1z - Spin2z)*ez[k])*(Spin1[k] - Spin2[k] - (Spin1z - Spin2z)*ez[k]);
    Theta += (Spin1[k] - Spin2[k] - (Spin1z - Spin2z)*ez[k])*e1[k];
  }
  Spin1xy = sqrt(Spin1xy);
  Spin2xy = sqrt(Spin2xy);
  Theta = acos(Theta/sqrt(Spin12xyMag));

  Vpar = K*cos(Theta-Theta0)*Eta*Eta / (1+Q)*(Spin2xy - Q*Spin1xy );

  /* Total spin contribution */
  for(k = 0; k<3; k++)
    RecoilVel[k] += Vper*(cos(Xi)*e1[k] + sin(Xi)*e2[k]) + Vpar*ez[k];

  /* Final spin of the remnant */
  double Lspin = 0;
  double cosa1a2 = 0, cosa1l = 0, cosa2l = 0;
  double a1 = 0, a2 = 0;
  /* Parameters for remnant's spin */
  double s4 = -0.129;	//[-]		Fanidakis, et al. 2010, Rezzolla, et al. 2008
  double s5 = -0.384;	//[-]		Fanidakis, et al. 2010, Rezzolla, et al. 2008
  double t0 = -2.686;	//[-]		Fanidakis, et al. 2010, Rezzolla, et al. 2008
  double t2 = -3.454;	//[-]		Fanidakis, et al. 2010, Rezzolla, et al. 2008
  double t3 = 2.353;	//[-]		Fanidakis, et al. 2010, Rezzolla, et al. 2008

  for(k = 0; k<3; k++)
  {
    a1 += Spin1[k]*Spin1[k];
    a2 += Spin2[k]*Spin2[k];
    cosa1a2 += Spin1[k]*Spin2[k];
    cosa1l += Spin1[k]*ez[k];
    cosa2l += Spin2[k]*ez[k];
  }
  a1 = sqrt(a1);
  a2 = sqrt(a2);
  cosa1a2 = cosa1a2/(a1*a2);
  cosa1l = cosa1l/a1;
  cosa2l = cosa2l/a2;

  Lspin = s4/pow(1 + Q*Q,2)*(a2*a2 + a1*a1*Q*Q*Q*Q + 2*a1*a2*Q*Q*cosa1a2) +
		  ((s5*Eta + t0 + 2)/( 1 + Q*Q ))*(a2*cosa2l + a1*Q*Q*cosa1l) +
		  2*sqrt(3) + t2*Eta + t3*Eta*Eta;

  for(k = 0; k<3; k++)
    FinalSpin[k] = 1/((1+Q)*(1+Q))*( Spin2[k] + Q*Q*Spin1[k] + Lspin*Q*ez[k] );

#endif

  /* Correction due to orbital eccentricity */
  eccentricity = sqrt(1 + 2*(0.5*Vmag*Vmag - All.G*(Mass1 + Mass2)/Rmag)*Lmag*Lmag / (All.G*All.G*(Mass1 + Mass2)*(Mass1 + Mass2)));
  for(k = 0; k<3; k++)
    RecoilVel[k] = cosmofac*(1 + 0*eccentricity)*RecoilVel[k]*1e5/All.UnitVelocity_in_cm_per_s;

  //DEBUG
  printf("BLACK_HOLES R: Applying recoil kick: Vm=%e\t Vper=%e\t Vpar=%e\t ecc=%e.\n",Vm, Vper, Vpar, eccentricity);
  printf("BLACK_HOLES R: Applying recoil kick: Vx=%e\t Vy=%e\t Vz=%e\n", RecoilVel[0], RecoilVel[1], RecoilVel[2]);
  printf("BLACK_HOLES R: Masses: [%e, %e]\n", Mass1, Mass2);
  printf("BLACK_HOLES R: R1 components: [%e, %e, %e]\n",Pos1[0], Pos1[1], Pos1[2]);
  printf("BLACK_HOLES R: R2 components: [%e, %e, %e]\n",Pos2[0], Pos2[1], Pos2[2]);
  printf("BLACK_HOLES R: R components: [%e, %e, %e]\n",r[0], r[1], r[2]);
  printf("BLACK_HOLES R: V1 components: [%e, %e, %e]\n",Vel1[0], Vel1[1], Vel1[2]);
  printf("BLACK_HOLES R: V2 components: [%e, %e, %e]\n",Vel2[0], Vel2[1], Vel2[2]);
  printf("BLACK_HOLES R: V components: [%e, %e, %e]\n",v[0], v[1], v[2]);
  printf("BLACK_HOLES R: S1 components: [%e, %e, %e]\n",Spin1[0], Spin1[1], Spin1[2]);
  printf("BLACK_HOLES R: S2 components: [%e, %e, %e]\n",Spin2[0], Spin2[1], Spin2[2]);
  printf("BLACK_HOLES R: Sf components: [%e, %e, %e]\n",FinalSpin[0], FinalSpin[1], FinalSpin[2]);
}
#endif


/* This routine uses the gravitational tree to search in the BH_Hsml neighborhood of
 * an active BH for other BHs that are to be merged, and carries out the mergers.
 * Only BHs with SwallowID==0 can swallow other BHs because they are guaranteed to survive themselves.
 */
void blackhole_do_mergers(void)
{
  int idx, i, j, k, n, ncount, nexport, nimport;
  int Ntot_BH_swallowed_local, Ntot_BH_swallowed_imported;
  int ngrp, recvTask;

  TIMER_START(CPU_BH_MERGERS);

  mpi_printf("BLACK_HOLES: Begin BH mergers.\n");

  treeBHexportflag = mymalloc_clear("treeBHexportflag", Tree_NumBHImported * sizeof(char));

#ifdef TRACER_MC_CHECKS
  long long check_total_tracers = get_total_number_of_tracers(-1);
  long long check_total_bh_tracers_rear = get_total_number_of_tracers(5);
  long long check_total_star_tracers_rear = get_total_number_of_tracers(4);
#endif

#ifdef TRACER_MC
  treeBHexportRindex = mymalloc_clear("treeBHexportRindex", Tree_NumBHImported * sizeof(int));
  treeBHexportRtask = mymalloc_clear("treeBHexportRtask", Tree_NumBHImported * sizeof(int));
  treeBHexportRID = mymalloc_clear("treeBHexportRID", Tree_NumBHImported * sizeof(MyIDType));
#endif

  N_BH_swallowed_local = 0;
  N_BH_swallowed_imported = 0;

  /* allocate temporary variables */
  BH_accreted = mymalloc_clear("BH_accreted", NumBHs * sizeof(struct bh_properties));

  generic_set_MaxNexport();

  /* Create list of is. We do this here to simplify the treatment of the two possible locations of source points */
  TargetList = mymalloc("TargetList", (NumPart + Tree_NumPartImported) * sizeof(int));
  Tree_ResultIndexList = mymalloc("Tree_ResultIndexList", Tree_NumPartImported * sizeof(int));

  Nforces = 0;
  for(idx = 0; idx < TimeBinsBHAccretion.NActiveParticles; idx++)
    {
      i = TimeBinsBHAccretion.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(Tree_Task_list[i] == ThisTask)
        if(BPP(i).SwallowID == 0)
          TargetList[Nforces++] = i;
    }

  for(i = 0, ncount = 0; i < Tree_NumPartImported; i++)
#ifndef HIERARCHICAL_GRAVITY
    if(Tree_Points[i].ActiveFlag)
#endif
      if(Tree_Points[i].Type == 5)
        if(TBPP(i).SwallowID == 0)
          {
            Tree_ResultIndexList[i] = ncount++;
            TargetList[Nforces++] = i + Tree_ImportedNodeOffset;
          }


#ifdef TRACER_MC
  start_MC_tracer(N_tracer);    /* allocate buffer for tracer exchange */
#endif


  Tree_BhngbResultsImported = mymalloc_clear("Tree_BhngbResultsImported", ncount * sizeof(struct bhresultsimported_data));


  generic_comm_pattern(Nforces, kernel_local, kernel_imported);


  /* now communicate the results in Tree_BhngbResultsImported */

  for(j = 0; j < NTask; j++)
    Recv_count[j] = 0;

  for(i = 0, n = 0, k = 0; i < NTask; i++)
    for(j = 0; j < Force_Recv_count[i]; j++, n++)
      {
#ifndef HIERARCHICAL_GRAVITY
        if(Tree_Points[n].ActiveFlag)
#endif
          if(Tree_Points[n].Type == 5)
            if(TBPP(n).SwallowID == 0)
              {
                Tree_BhngbResultsImported[k].index = Tree_Points[n].index;
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

  struct bhresultsimported_data *tmp_results = mymalloc("tmp_results", nexport * sizeof(struct bhresultsimported_data));
  memset(tmp_results, -1, nexport * sizeof(struct bhresultsimported_data));

  /* exchange data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          MPI_Sendrecv(&Tree_BhngbResultsImported[Recv_offset[recvTask]],
                       Recv_count[recvTask] * sizeof(struct bhresultsimported_data), MPI_BYTE, recvTask, TAG_FOF_A,
                       &tmp_results[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct bhresultsimported_data), MPI_BYTE, recvTask, TAG_FOF_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

  for(i = 0; i < nexport; i++)
    {
      int target = tmp_results[i].index;
      if(P[target].AuxDataID >= NumBHs)
        terminate("BLACK_HOLES: P[target(=%d)].AuxDataID(=%lld) >= NumBHs=%d", target, (long long) P[target].AuxDataID, NumBHs);

      BH_accreted[P[target].AuxDataID] = tmp_results[i].Prop;
    }

  myfree(tmp_results);
  myfree(Tree_BhngbResultsImported);

  /* now in case we have swallowed remote black holes out of Tree_Points, need to erase them locally too */
  Tree_BhngbDeleted = mymalloc("Tree_BhngbDeleted", sizeof(struct bhdeleted_data) * N_BH_swallowed_imported);

  for(j = 0; j < NTask; j++)
    Recv_count[j] = 0;

  for(i = 0, n = 0, k = 0; i < NTask; i++)
    for(j = 0; j < Force_Recv_count[i]; j++, n++)
      {
        if(Tree_Points[n].Type  == 5 && treeBHexportflag[Tree_Points[n].AuxDataIndex] == 1)
          {
            if(k >= N_BH_swallowed_imported)
              terminate("BLACK_HOLES: k >= N_BH_swallowed k=%d N_BH_swallowed_imported=%d n=%d i=%d Tree_NumPartImported=%d", k, N_BH_swallowed_imported, n, i, Tree_NumPartImported);

            Tree_BhngbDeleted[k].index = Tree_Points[n].index;
#ifdef TRACER_MC
            /* we will delete P[Tree_Points[n].index] on some other task, and at that time send
               its child tracers back to this task. note the local BH that swallowed this remote
               BH so that we can attach these imported tracers. */
            Tree_BhngbDeleted[k].rindex = treeBHexportRindex[Tree_Points[n].AuxDataIndex];
            Tree_BhngbDeleted[k].rtask = treeBHexportRtask[Tree_Points[n].AuxDataIndex];
            Tree_BhngbDeleted[k].rID = treeBHexportRID[Tree_Points[n].AuxDataIndex];
#endif
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

  struct bhdeleted_data *tmp_deleted = mymalloc("tmp_deleted", nexport * sizeof(struct bhdeleted_data));

  /* exchange data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          MPI_Sendrecv(&Tree_BhngbDeleted[Recv_offset[recvTask]],
                       Recv_count[recvTask] * sizeof(struct bhdeleted_data), MPI_BYTE, recvTask, TAG_FOF_A,
                       &tmp_deleted[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct bhdeleted_data), MPI_BYTE, recvTask, TAG_FOF_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

  for(i = 0; i < nexport; i++)
    {
      int no = tmp_deleted[i].index;
      if(P[no].Type != 5)
        terminate("BLACK_HOLES: P[no].Type != 5 when deleting BHs");

      /* deleting this black hole */
      int bin = P[no].TimeBinHydro;

      TimeBin_BH_mass[bin] -= BPP(no).BH_Mass;
      TimeBin_BH_dynamicalmass[bin] -= P[no].Mass;
      TimeBin_BH_Mdot[bin] -= BPP(no).BH_Mdot;
      if(BPP(no).BH_Mass > 0)
        TimeBin_BH_Medd[bin] -= BPP(no).BH_Mdot / BPP(no).BH_Mass;

#ifdef TRACER_MC
      consider_moving_tracers(no, tmp_deleted[i].rtask, tmp_deleted[i].rindex, tmp_deleted[i].rID, 1.0);
#endif

      P[no].Mass = 0;
      P[no].ID = 0;

      BPP(no).BH_Mass = 0;
      BPP(no).BH_Mdot = 0;

#ifdef BH_BUBBLES
      BPP(no).BH_Mass_bubbles = 0;
      BPP(no).BH_Mass_ini = 0;
#endif
      //FIXME set zero for new radio mode
    }

  myfree(tmp_deleted);
  myfree(Tree_BhngbDeleted);


#ifdef TRACER_MC
  finish_MC_tracer();
#endif

  myfree(Tree_ResultIndexList);
  myfree(TargetList);


  /* now update momentum of BH */
  for(idx = 0; idx < TimeBinsBHAccretion.NActiveParticles; idx++)
    {
      n = TimeBinsBHAccretion.ActiveParticleList[idx];
      if(n < 0)
        continue;

      if(P[n].AuxDataID >= NumBHs)
        terminate("BLACK_HOLES: P[n(=%d)].AuxDataID(=%lld) >= NumBHs=%d", n, (long long) P[n].AuxDataID, NumBHs);

      if(BH_accreted[P[n].AuxDataID].Mass > 0)
        {
#ifdef BH_RECOIL_KICK
    	  double RecoilVel[3] = {};
    	  double FinalSpin[3] = {};
    	  double M1, M2, Pos1[3], Pos2[3], Vel1[3], Vel2[3];
    	  double Spin1[3] = {0}, Spin2[3] = {0};
    	  double SpinParameter = 0;
    	  MyFloat CurrentTime;

		  M1 = P[n].Mass;
		  M2 = BH_accreted[P[n].AuxDataID].Mass;
		  for(k = 0; k < 3; k++)
		  {
	        Pos1[k] = P[n].Pos[k];
	        Pos2[k] = BH_accreted[P[n].AuxDataID].Pos[k];
	        Vel1[k] = P[n].Vel[k];
	        Vel2[k] = BH_accreted[P[n].AuxDataID].Vel[k];
		  }
#ifdef BH_SPIN_EVOLUTION
		  for(k = 0; k < 3; k++)
		  {
		    Spin1[k] = BPP(n).BH_SpinParameter*BPP(n).BH_SpinOrientation[k];
		    Spin2[k] = BH_accreted[P[n].AuxDataID].Spin[k];
		  }
#endif

    	  if(P[n].Mass < BH_accreted[P[n].AuxDataID].Mass)
    	    apply_recoil_kick(M1, M2, Pos1, Pos2, Vel1, Vel2, Spin1, Spin2, RecoilVel, FinalSpin);
    	  else
    	    apply_recoil_kick(M2, M1, Pos2, Pos1, Vel2, Vel1, Spin2, Spin1, RecoilVel, FinalSpin);

          for(k = 0; k < 3; k++)
        	P[n].Pos[k] = (P[n].Mass*P[n].Pos[k] + BH_accreted[P[n].AuxDataID].Mass*BH_accreted[P[n].AuxDataID].Pos[k]) / (P[n].Mass + BH_accreted[P[n].AuxDataID].Mass);

#ifdef BH_SPIN_EVOLUTION
          for(k = 0; k < 3; k++)
            SpinParameter += FinalSpin[k]*FinalSpin[k];
          SpinParameter = sqrt(SpinParameter);
          BPP(n).BH_SpinParameter = SpinParameter;
          for(k = 0; k < 3; k++)
        	BPP(n).BH_SpinOrientation[k] = FinalSpin[k]/SpinParameter;

          if(All.ComovingIntegrationOn)
        	CurrentTime = All.TimeBegin * exp(P[n].Ti_Current * All.Timebase_interval);
          else
            CurrentTime = All.TimeBegin + P[n].Ti_Current * All.Timebase_interval;
          CurrentTime = get_time_difference_in_Gyr(All.TimeBegin, CurrentTime);

          blackhole_new_accretion_episode(n, CurrentTime, BPP(n).BH_Mass + BH_accreted[P[n].AuxDataID].BHMass);
#endif

#endif

          for(k = 0; k < 3; k++)
            P[n].Vel[k] = (P[n].Vel[k] * P[n].Mass + BH_accreted[P[n].AuxDataID].Momentum[k]) / (P[n].Mass + BH_accreted[P[n].AuxDataID].Mass);

#ifdef BH_RECOIL_KICK
          for(k = 0; k < 3; k++)
        	  P[n].Vel[k] = P[n].Vel[k] + RecoilVel[k];
#endif

          P[n].Mass += BH_accreted[P[n].AuxDataID].Mass;

#ifdef INDIVIDUAL_GRAVITY_SOFTENING
          if(((1 << P[n].Type) & (INDIVIDUAL_GRAVITY_SOFTENING)))
            P[n].SofteningType = get_softening_type_from_mass(P[n].Mass);
#endif

          BPP(n).BH_Mass += BH_accreted[P[n].AuxDataID].BHMass;
          BPP(n).BH_CumMass_QM += BH_accreted[P[n].AuxDataID].BH_CumMass_QM;
          BPP(n).BH_CumEgy_QM += BH_accreted[P[n].AuxDataID].BH_CumEgy_QM;
          BPP(n).BH_CumMass_RM += BH_accreted[P[n].AuxDataID].BH_CumMass_RM;
          BPP(n).BH_CumEgy_RM += BH_accreted[P[n].AuxDataID].BH_CumEgy_RM;

#ifdef DRAINGAS
          BPP(n).DrainBucketMass += BH_accreted[P[n].AuxDataID].DrainBucketMass;
#endif
          BPP(n).BH_CountProgs += BH_accreted[P[n].AuxDataID].BH_CountProgs;

#ifdef BH_BUBBLES
          BPP(n).BH_Mass_bubbles += BH_accreted[P[n].AuxDataID].BHMass_bubbles;
          BPP(n).BH_Mass_ini += BH_accreted[P[n].AuxDataID].BHMass_ini;
#endif
#ifdef BH_NF_RADIO
          BPP(n).BH_RadioEgyFeedback += BH_accreted[P[n].AuxDataID].BH_RadioEgyFeedback;
          BPP(n).BH_HaloVvir = dmax(BPP(n).BH_HaloVvir, BH_accreted[P[n].AuxDataID].BH_HaloVvir);
          BPP(n).BH_Mdot_radio = dmax(BPP(n).BH_Mdot_radio, BH_accreted[P[n].AuxDataID].BH_Mdot_radio);
          BPP(n).BH_Mdot_quasar = dmax(BPP(n).BH_Mdot_quasar, BH_accreted[P[n].AuxDataID].BH_Mdot_quasar);
          BPP(n).BH_XrayLum = dmax(BPP(n).BH_XrayLum, BH_accreted[P[n].AuxDataID].BH_XrayLum);
          BPP(n).BH_RadioLum = dmax(BPP(n).BH_RadioLum, BH_accreted[P[n].AuxDataID].BH_RadioLum);
#endif
        }

      if(P[n].Mass == 0 && P[n].ID == 0)
        timebin_remove_particle(&TimeBinsBHAccretion, idx, P[n].TimeBinHydro);
    }


  timebin_cleanup_list_of_active_particles(&TimeBinsGravity);

  MPI_Reduce(&N_BH_swallowed_imported, &Ntot_BH_swallowed_imported, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&N_BH_swallowed_local, &Ntot_BH_swallowed_local, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  myfree(BH_accreted);

#ifdef TRACER_MC
  myfree(treeBHexportRID);
  myfree(treeBHexportRtask);
  myfree(treeBHexportRindex);
#endif

#ifdef TRACER_MC_CHECKS
  if(check_total_tracers != get_total_number_of_tracers(-1))
    terminate("TRACER_MC: unconserved tracers during BLACK_HOLES (merger): total number of tracers BEFORE = %lld, AFTER = %lld\n", check_total_tracers, get_total_number_of_tracers(-1));
  if(check_total_bh_tracers_rear != get_total_number_of_tracers(5))
    terminate("TRACER_MC: strange, global number of BH tracers changed in MERGER, BEFORE = %lld, AFTER = %lld\n", check_total_bh_tracers_rear, get_total_number_of_tracers(5));
  if(check_total_star_tracers_rear != get_total_number_of_tracers(4))
    terminate("TRACER_MC: strange, global number of STAR tracers changed in bh MERGER, BEFORE = %lld, AFTER = %lld\n", check_total_star_tracers_rear, get_total_number_of_tracers(4));
#endif

  myfree(treeBHexportflag);
  mpi_printf("BLACK_HOLES: Black hole merging done: %d (%d/%d) BH particles swallowed\n", Ntot_BH_swallowed_local + Ntot_BH_swallowed_imported, Ntot_BH_swallowed_local, Ntot_BH_swallowed_imported);

  TIMER_STOP(CPU_BH_MERGERS);
}






static int blackhole_tree_merger_evaluate(int target, int mode, int threadid)
{
  int k, no, numnodes, *firstnode;
  double h, h2;
#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif
  double dx, dy, dz, r2;
  MyIDType id;
  MyDouble *pos;
  MyDouble mass;

  data_in local, *in;
  data_out out;

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

  pos = in->Pos;
  mass = in->Mass;
  h = in->BH_Hsml;
  id = in->ID;

  h2 = h * h;

  struct bh_properties accreted;
  memset(&accreted, 0, sizeof(struct bh_properties));

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

              if(r2 < h2)
                {
                  if(P[no].Type == 5)
                    {
                      if((BPP(no).SwallowID == id) && (id != 0) && (P[no].ID != 0))     /* we have a black hole merger */
                        {
                          fprintf(FdBlackHolesMergers, "%d %g %llu %g %llu %g\n", ThisTask, All.Time, (long long) id, mass, (long long) P[no].ID, BPP(no).BH_Mass);
                          myflush(FdBlackHolesMergers);

                          accreted.Mass += P[no].Mass;
                          accreted.BHMass += BPP(no).BH_Mass;
                          accreted.BH_CumMass_QM += BPP(no).BH_CumMass_QM;
                          accreted.BH_CumEgy_QM += BPP(no).BH_CumEgy_QM;
                          accreted.BH_CumMass_RM += BPP(no).BH_CumMass_RM;
                          accreted.BH_CumEgy_RM += BPP(no).BH_CumEgy_RM;
                          for(k = 0; k < 3; k++)
                            accreted.Momentum[k] += P[no].Mass * P[no].Vel[k];
#ifdef BH_RECOIL_KICK
                          for(k = 0; k < 3; k++)
                          {
                        	accreted.Pos[k] += P[no].Pos[k];
                        	accreted.Vel[k] += P[no].Vel[k];
#ifdef BH_SPIN_EVOLUTION
                        	accreted.Spin[k] += BPP(no).BH_SpinParameter*BPP(no).BH_SpinOrientation[k];
#endif
                          }
#endif
                          accreted.BH_CountProgs += BPP(no).BH_CountProgs;

#ifdef TRACER_MC
                          consider_moving_tracers(no, in->rtask, in->rindex, in->ID, 1.0);
#endif

                          int bin = P[no].TimeBinHydro;

                          TimeBin_BH_mass[bin] -= BPP(no).BH_Mass;
                          TimeBin_BH_dynamicalmass[bin] -= P[no].Mass;
                          TimeBin_BH_Mdot[bin] -= BPP(no).BH_Mdot;
                          if(BPP(no).BH_Mass > 0)
                            TimeBin_BH_Medd[bin] -= BPP(no).BH_Mdot / BPP(no).BH_Mass;

                          P[no].Mass = 0;
                          P[no].ID = 0;

                          BPP(no).BH_Mass = 0;
                          BPP(no).BH_Mdot = 0;

#ifdef DRAINGAS
                          accreted.DrainBucketMass += BPP(no).DrainBucketMass;
#endif

#ifdef BH_BUBBLES
                          accreted.BHMass_bubbles += BPP(no).BH_Mass_bubbles;
                          accreted.BHMass_ini += BPP(no).BH_Mass_ini;
                          BPP(no).BH_Mass_bubbles = 0;
                          BPP(no).BH_Mass_ini = 0;
#endif

#ifdef BH_NF_RADIO
                          accreted.BH_RadioEgyFeedback += BPP(no).BH_RadioEgyFeedback;
                          accreted.BH_Mdot_quasar = dmax(accreted.BH_Mdot_quasar, BPP(no).BH_Mdot_quasar);
                          accreted.BH_Mdot_radio = dmax(accreted.BH_Mdot_radio, BPP(no).BH_Mdot_radio);
                          accreted.BH_HaloVvir = dmax(accreted.BH_HaloVvir, BPP(no).BH_HaloVvir);
                          accreted.BH_XrayLum = dmax(accreted.BH_XrayLum, BPP(no).BH_XrayLum);
                          accreted.BH_RadioLum = dmax(accreted.BH_RadioLum, BPP(no).BH_RadioLum);
#endif

                          N_BH_swallowed_local++;
                        }
                    }
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

              if(r2 < h2)
                {
                  if(Tree_Points[n].Type == 5)   /* we have a potential black hole merger */
                    {
                      if((TBPP(n).SwallowID == id) && (id != 0) && (TBPP(n).ID != 0))   /* we have a black hole merger */
                        {
                          fprintf(FdBlackHolesMergers, "%d %g %llu %g %llu %g\n", ThisTask, All.Time, (long long) id, mass, (long long) TBPP(n).ID, TBPP(n).BH_Mass);
                          myflush(FdBlackHolesMergers);

                          accreted.Mass += Tree_Points[n].Mass;
                          accreted.BHMass += TBPP(n).BH_Mass;
                          accreted.BH_CumMass_QM += TBPP(n).BH_CumMass_QM;
                          accreted.BH_CumEgy_QM += TBPP(n).BH_CumEgy_QM;
                          accreted.BH_CumMass_RM += TBPP(n).BH_CumMass_RM;
                          accreted.BH_CumEgy_RM += TBPP(n).BH_CumEgy_RM;

                          for(k = 0; k < 3; k++)
                            accreted.Momentum[k] += Tree_Points[n].Mass * TBPP(n).Vel[k];
#ifdef BH_RECOIL_KICK
                          for(k = 0; k < 3; k++)
                          {
                        	accreted.Pos[k] += Tree_Points[n].Pos[k];
                        	accreted.Vel[k] += TBPP(n).Vel[k];
#ifdef BH_SPIN_EVOLUTION
                        	accreted.Spin[k] += TBPP(n).BH_SpinParameter*TBPP(n).BH_SpinOrientation[k];
#endif
                          }
#endif
                          accreted.BH_CountProgs += TBPP(n).BH_CountProgs;
#ifdef DRAINGAS
                          accreted.DrainBucketMass += TBPP(n).DrainBucketMass;
#endif
#ifdef BH_BUBBLES
                          accreted.BHMass_bubbles += TBPP(n).BH_Mass_bubbles;
                          accreted.BHMass_ini += TBPP(n).BH_Mass_ini;
#endif
#ifdef BH_NF_RADIO
                          accreted.BH_RadioEgyFeedback += TBPP(n).BH_RadioEgyFeedback;
                          accreted.BH_Mdot_quasar = dmax(accreted.BH_Mdot_quasar, TBPP(n).BH_Mdot_quasar);
                          accreted.BH_Mdot_radio = dmax(accreted.BH_Mdot_radio, TBPP(n).BH_Mdot_radio);
                          accreted.BH_HaloVvir = dmax(accreted.BH_HaloVvir, TBPP(n).BH_HaloVvir);
                          accreted.BH_XrayLum = dmax(accreted.BH_XrayLum, TBPP(n).BH_XrayLum);
                          accreted.BH_RadioLum = dmax(accreted.BH_RadioLum, TBPP(n).BH_RadioLum);
#endif

                          TBPP(n).ID = 0;
                          treeBHexportflag[Tree_Points[n].AuxDataIndex] = 1;

#ifdef TRACER_MC
                          /* 3 cases are all handled by setting id,rindex,rtask previously:
                             (a) mode=0 (local P search): we found a Tree_Points BH to swallow,
                             request this found BH be deleted later, tracers come back to this task
                             (b) mode=0 (local Tree search): we found a Tree_Points BH to swallow,
                             request this found BH be deleted later, but tracers go to different task
                             (c) mode=1 (remote search, P or Tree) we found a Tree_Points BH to swallow,
                             request this found BH be deleted later, but tracers go to different task */
                          treeBHexportRindex[Tree_Points[n].AuxDataIndex] = in->rindex;
                          treeBHexportRtask[Tree_Points[n].AuxDataIndex] = in->rtask;
                          treeBHexportRID[Tree_Points[n].AuxDataIndex] = in->ID;
#endif

                          Tree_Points[n].Mass = 0;

                          N_BH_swallowed_imported++;
                        }
                    }
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

  out.Prop = accreted;

  /* store result at the proper place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target].Prop = accreted;

  return 0;
}

#endif