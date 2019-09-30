/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/blackhole/blackhole_swallowgas.c
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


#if defined(BLACK_HOLES)
 

static int blackhole_evaluate_swallow(int target, int mode, int threadid);

static int NumGas_swallowed, Ntot_gas_swallowed;
static MyDouble *BH_accreted_Mass;
static MyFloat *BH_accreted_momentum;


/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble Pos[3];
  MyFloat BH_Hsml;
  MyIDType ID;
  MyFloat BH_Mass;
#ifdef DRAINGAS
  MyFloat Mdot;
  MyFloat Dt;
  MyDouble Mass;
  MyIDType DrainID;
  MyDouble DrainBucketMass;
#if (DRAINGAS == 3)
  MyDouble BH_VolSum;
#endif
#endif
#ifdef TRACER_MC
  int rtask;
  int rindex;
#endif
#ifdef BH_BIPOLAR_FEEDBACK
  MyDouble BH_BipolarJ[3];
#endif

  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;


 /* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in * in, int i, int firstnode)
{
  in->Pos[0] = P[i].Pos[0];
  in->Pos[1] = P[i].Pos[1];
  in->Pos[2] = P[i].Pos[2];
  in->BH_Hsml = BPP(i).BH_Hsml;
  in->ID = P[i].ID;
  in->BH_Mass = BPP(i).BH_Mass;
#ifdef DRAINGAS
  in->Mdot = BPP(i).BH_Mdot;
  in->Dt = (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval / All.cf_hubble_a;
  in->Mass = P[i].Mass;
  in->DrainID = BPP(i).DrainID;
  in->DrainBucketMass = BPP(i).DrainBucketMass;
#if (DRAINGAS == 3)
  in->BH_VolSum = BPP(i).BH_VolSum;
#endif
#endif
#ifdef TRACER_MC
  in->rtask = ThisTask;
  in->rindex = i;
#endif

  in->Firstnode = firstnode;
}

/* local data structure that holds results acquired on remote processors */
typedef struct
{
  MyDouble Mass;
  MyFloat AccretedMomentum[3];
#ifdef DRAINGAS
  MyDouble DrainBucketMass;
#endif
} data_out;

static data_out *DataResult, *DataOut;


 /* routine to store or combine result data */
static void out2particle(data_out * out, int i, int mode)
{
  int k;

  if(mode == MODE_LOCAL_PARTICLES)      /* initial store */
    {
      if(P[i].AuxDataID >= NumBHs)
        terminate("BLACK_HOLES: P[i(=%d)].AuxDataID >= NumBHs=%d", i, NumBHs);

      BH_accreted_Mass[P[i].AuxDataID] = out->Mass;
      for(k = 0; k < 3; k++)
        BH_accreted_momentum[3 * P[i].AuxDataID + k] = out->AccretedMomentum[k];
#ifdef DRAINGAS
      BPP(i).DrainBucketMass = out->DrainBucketMass;
#endif
    }
  else                          /* combine */
    {
      if(P[i].AuxDataID >= NumBHs)
        terminate("BLACK_HOLES: P[i(=%d)].AuxDataID >= NumBHs=%d", i, NumBHs);

      BH_accreted_Mass[P[i].AuxDataID] += out->Mass;

      for(k = 0; k < 3; k++)
        BH_accreted_momentum[3 * P[i].AuxDataID + k] += out->AccretedMomentum[k];

#ifdef DRAINGAS
      BPP(i).DrainBucketMass += out->DrainBucketMass;
#endif
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

        if(BPP(i).SwallowID == 0)
          blackhole_evaluate_swallow(i, MODE_LOCAL_PARTICLES, threadid);
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

        blackhole_evaluate_swallow(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}


void blackhole_swallow_gas(void)
{
  int idx, n;

  mpi_printf("BLACK_HOLES: Begin BH gas swallowing.\n");
  double t0 = second();

  /* Now do the swallowing of particles */
  NumGas_swallowed = 0;

  /* allocate temporary variables */
  BH_accreted_Mass = (MyDouble *) mymalloc("BH_accreted_Mass", NumBHs * sizeof(MyDouble));
  memset(BH_accreted_Mass, 0, NumBHs * sizeof(MyDouble));
  BH_accreted_momentum = (MyFloat *) mymalloc("BH_accreted_momentum", 3 * NumBHs * sizeof(MyFloat));
  memset(BH_accreted_momentum, 0, 3 * NumBHs * sizeof(MyFloat));

  generic_set_MaxNexport();

#ifdef TRACER_MC
  start_MC_tracer(N_tracer);    /* allocate buffer for tracer exchange */
#endif

  generic_comm_pattern(TimeBinsBHAccretion.NActiveParticles, kernel_local, kernel_imported);


#ifdef TRACER_MC
  finish_MC_tracer();
#endif


  for(idx = 0; idx < TimeBinsBHAccretion.NActiveParticles; idx++)
    {
      n = TimeBinsBHAccretion.ActiveParticleList[idx];
      if(n < 0)
        continue;

      if(P[n].AuxDataID >= NumBHs)
        terminate("BLACK_HOLES: P[n(=%d)].AuxDataID(=%lld) >= NumBHs=%d", n, (long long) P[n].AuxDataID, NumBHs);

      if(BH_accreted_Mass[P[n].AuxDataID] > 0)
        {
#if !defined(BH_FRICTION) && !defined(BH_NEW_CENTERING)
          int k;
          for(k = 0; k < 3; k++)
            P[n].Vel[k] = (P[n].Vel[k] * P[n].Mass + BH_accreted_momentum[3 * P[n].AuxDataID + k]) / (P[n].Mass + BH_accreted_Mass[P[n].AuxDataID]);
#endif
          P[n].Mass += BH_accreted_Mass[P[n].AuxDataID];

#ifdef INDIVIDUAL_GRAVITY_SOFTENING
          if(((1 << P[n].Type) & (INDIVIDUAL_GRAVITY_SOFTENING)))
            P[n].SofteningType = get_softening_type_from_mass(P[n].Mass);
#endif
        }
    }

  MPI_Reduce(&NumGas_swallowed, &Ntot_gas_swallowed, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  myfree(BH_accreted_momentum);
  myfree(BH_accreted_Mass);

  double t1 = second();
  mpi_printf("BLACK_HOLES: Accretion done: %d gas particles swallowed/drained. (took %g sec)\n", Ntot_gas_swallowed, timediff(t0, t1));
}


static int blackhole_evaluate_swallow(int target, int mode, int threadid)
{
  int numnodes, *firstnode, j, k, n;
  MyIDType id;
  MyDouble accreted_mass;
  MyFloat accreted_momentum[3];
  MyDouble *pos;
  MyFloat h_i, bh_mass;

#ifdef DRAINGAS
  MyDouble mass, dt, mdot, dmass, fac;
  MyIDType drainID;
  MyDouble diff_drainBucketMass = 0;
#if (DRAINGAS == 3)
  double xtmp, ytmp, ztmp;
  MyDouble volsum, dx, dy, dz, r2, r, u, wk, h2, hinv, hinv3;
#else
  MyDouble drainBucketMass;
#endif
#endif

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
  h_i = in->BH_Hsml;
  id = in->ID;
  bh_mass = in->BH_Mass;
#ifdef DRAINGAS
  mass = in->Mass;
  mdot = in->Mdot;
  dt = in->Dt;
  drainID = in->DrainID;
#if (DRAINGAS == 3)
  volsum = in->BH_VolSum;
#else
  drainBucketMass = in->DrainBucketMass;
#endif
#endif
#ifdef BH_BIPOLAR_FEEDBACK
  MyDouble *bipolar_j = in->BH_BipolarJ;
#endif

  accreted_mass = 0;
  accreted_momentum[0] = accreted_momentum[1] = accreted_momentum[2] = 0;

#if (DRAINGAS == 3)
  h2 = h_i * h_i;
  hinv = 1.0 / h_i;
  hinv3 = hinv * hinv * hinv;
#endif


  int nfound = ngb_treefind_variable_threads(pos, h_i, target, mode, threadid, numnodes, firstnode);

  for(n = 0; n < nfound; n++)
    {
      j = Thread[threadid].Ngblist[n];

      /* gas cell -> drain or swallow */
      if(P[j].Type == 0)
        {
#ifdef DRAINGAS
          /* continuous draining; note: draining only starts once subgrid mass >= dynamical mass */
#if (DRAINGAS == 3)
          dx = NGB_PERIODIC_LONG_X(pos[0] - P[j].Pos[0]);
          dy = NGB_PERIODIC_LONG_Y(pos[1] - P[j].Pos[1]);
          dz = NGB_PERIODIC_LONG_Z(pos[2] - P[j].Pos[2]);

          r2 = dx * dx + dy * dy + dz * dz;

          if(r2 < h2 && P[j].Mass > 0 && P[j].ID != 0 && (bh_mass >= mass))
            {

#ifdef BH_BIPOLAR_FEEDBACK
              if (!is_cell_in_bipolar_cone(j,pos,bipolar_j)) {
                continue;
              }
#endif

              dmass = (1. - All.BlackHoleRadiativeEfficiency) * mdot * dt;

              /* predicted final dynamical mass for the BH: if it is still lower than the subgrid   */
              /* mass, the mdot is tempararily increased to gently bring the two masses at the same */
              /* value. Note that this does not affect the feedback energy released by the BH which */
              /* is computed from the subgrid mass only                                             */
              if((bh_mass - mass) > dmass)
                dmass = 0.5 * (bh_mass - mass + dmass);

              r = sqrt(r2);

              u = r * hinv;

              if(u < 0.5)
                wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
              else
                wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

              dmass *= wk * SphP[j].Volume / volsum;

              if(volsum <= 0.0)
                terminate
                  ("The BH seems not active, but it tries to accrete mass nevertheless. volsum=%g, BH ID=%llu, BH hsml=%g, BH drainID=%llu, celld=%g",
                   volsum, (unsigned long long) id, h_i, (unsigned long long) drainID, r);

              fac = (P[j].Mass - dmass) / P[j].Mass;

              if(fac <= 0.01)
                {
                  warn("Trying to accrete more than the mass available from a cell: j=%d, fac=%g, P[j].Mass=%g, dmass=%g", j, fac, P[j].Mass, dmass);

                  fac = 0.01;
                  dmass = (1 - fac) * P[j].Mass;
                }

              accreted_mass += dmass;

#ifdef TRACER_MC
              consider_moving_tracers(j, in->rtask, in->rindex, in->ID, 1 - fac);
#endif

              for(k = 0; k < 3; k++)
                accreted_momentum[k] += dmass * P[j].Vel[k];

              P[j].Mass *= fac;
              SphP[j].Energy *= fac;
              SphP[j].Momentum[0] *= fac;
              SphP[j].Momentum[1] *= fac;
              SphP[j].Momentum[2] *= fac;

#ifdef MAXSCALARS
              for(int s = 0; s < N_Scalar; s++)     /* Note, the changes in MATERIALS, HIGHRESGASMASS, etc., are treated as part of the Scalars */
                *(MyFloat *) (((char *) (&SphP[j])) + scalar_elements[s].offset_mass) *= fac;
#endif
#ifdef GFM_CHEMTAGS
              for(int k = 0; k < GFM_N_CHEM_TAGS; k++)
                SphP[j].MassMetalsChemTags[k] *= fac;
#endif

              /* count drained as swallowed for statistics */
              NumGas_swallowed++;
            }
#else
          if((drainID == P[j].ID) && (bh_mass >= mass))
            {
              dmass = (1. - All.BlackHoleRadiativeEfficiency) * mdot * dt;

              /* predicted final dynamical mass for the BH: if it is still lower than the subgrid   */
              /* mass, the mdot is tempararily increased to gently bring the two masses at the same */
              /* value. Note that this does not affect the feedback energy released by the BH which */
              /* is computed from the subgrid mass only                                             */
              if((bh_mass - mass) > dmass)
                dmass = 0.5 * (bh_mass - mass + dmass);

              /* fill the bucket if the cell mass is too low for the mass that should be drained */
              if(dmass > 0.9 * P[j].Mass)
                {
                  diff_drainBucketMass += (dmass - 0.9 * P[j].Mass);
                  warn
                    ("BLACK_HOLES DRAINGAS (All.Time=%g): filling bucket for BH: id=%llu mass=%g/%g drainBucketMass=%g; from cell: P[j].ID=%llu j=%d P[j].Mass=%g; dmass=%g diff_drainBucketMass=%g\n",
                     All.Time, (unsigned long long) id, mass, bh_mass, drainBucketMass, (unsigned long long) P[j].ID, j, P[j].Mass, dmass, diff_drainBucketMass);
                  dmass = 0.9 * P[j].Mass;
                }
              /* if bucket not empty, try to empty it */
              else
                {
                  if(drainBucketMass > 0)
                    {
                      if(0.9 * P[j].Mass > dmass + drainBucketMass)
                        {
                          warn
                            ("BLACK_HOLES DRAINGAS (All.Time=%g): fully emptying bucket of BH: id=%llu mass=%g/%g drainBucketMass=%g; and cell: P[j].ID=%llu j=%d P[j].Mass=%g; dmass=%g\n",
                             All.Time, (unsigned long long) id, mass, bh_mass, drainBucketMass, (unsigned long long) P[j].ID, j, P[j].Mass, dmass);
                          dmass += drainBucketMass;
                          diff_drainBucketMass += -drainBucketMass;
                        }
                      else
                        {
                          diff_drainBucketMass += -(0.9 * P[j].Mass - dmass);
                          warn
                            ("BLACK_HOLES DRAINGAS (All.Time=%g): partially emptying bucket of BH: id=%llu mass=%g/%g drainBucketMass=%g; with cell: P[j].ID=%llu j=%d P[j].Mass=%g; dmass=%g diff_drainBucketMass=%g\n",
                             All.Time, (unsigned long long) id, mass, bh_mass, drainBucketMass, (unsigned long long) P[j].ID, j, P[j].Mass, dmass, diff_drainBucketMass);
                          dmass = 0.9 * P[j].Mass;
                        }
                    }
                }

              fac = (P[j].Mass - dmass) / P[j].Mass;
              accreted_mass += dmass;

#ifdef TRACER_MC
              consider_moving_tracers(j, in->rtask, in->rindex, in->ID, 1 - fac);
#endif

              for(k = 0; k < 3; k++)
                accreted_momentum[k] += dmass * P[j].Vel[k];

              P[j].Mass *= fac;
              SphP[j].Energy *= fac;
              SphP[j].Momentum[0] *= fac;
              SphP[j].Momentum[1] *= fac;
              SphP[j].Momentum[2] *= fac;

#ifdef MAXSCALARS
              for(int s = 0; s < N_Scalar; s++)     /* Note, the changes in MATERIALS, HIGHRESGASMASS, etc., are treated as part of the Scalars */
                *(MyFloat *) (((char *) (&SphP[j])) + scalar_elements[s].offset_mass) *= fac;
#endif
#ifdef GFM_CHEMTAGS
              for(int k = 0; k < GFM_N_CHEM_TAGS; k++)
                SphP[j].MassMetalsChemTags[k] *= fac;
#endif
              /* count drained as swallowed for statistics */
              NumGas_swallowed++;
            }
#endif
#endif /* end DRAINGAS */
        }
    }


  out.Mass = accreted_mass;
  for(k = 0; k < 3; k++)
    out.AccretedMomentum[k] = accreted_momentum[k];
#ifdef DRAINGAS
  out.DrainBucketMass = diff_drainBucketMass;
#endif

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;
  return 0;
}

#endif
