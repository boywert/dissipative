/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/blackhole/blackhole_density.c
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
/*! Structure for communication during the density computation. Holds data that is sent to other processors.
 */

#define MAX_CONTRIBUTION_OF_SINGLE_CELL 2.5

#if defined(BH_USE_ALFVEN_SPEED_IN_BONDI) && !defined(MHD)
#error "BH_USE_ALFVEN_SPEED_IN_BONDI requires MHD"
#endif


static int blackhole_density_evaluate(int target, int mode, int threadid);
static MyFloat *Dhsmlrho;

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble Pos[3];
  MyFloat BH_Hsml;

#ifdef BH_BIPOLAR_FEEDBACK
  MyFloat BH_Vel[3];
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
#ifdef BH_BIPOLAR_FEEDBACK
  in->BH_Vel[0] = P[i].Vel[0];
  in->BH_Vel[1] = P[i].Vel[1];
  in->BH_Vel[2] = P[i].Vel[2];
#endif

  in->Firstnode = firstnode;
}

/* local data structure that holds results acquired on remote processors */
typedef struct
{
  MyFloat Rho;
  MyFloat VolSum;
  MyFloat Dhsmlrho;
  MyFloat SmoothedU;
  MyFloat GasVel[3];
  MyFloat Ngb;
  MyFloat DtGasNeighbor;
#ifdef DRAINGAS
  MyFloat NearestDist;
  MyIDType DrainID;
  MyFloat CellDensity;
#endif
#ifdef BH_USE_ALFVEN_SPEED_IN_BONDI
  MyFloat Bpress;
#endif
#if (defined(GFM_WINDS_VARIABLE) && (GFM_WINDS_VARIABLE==1)) || defined(GFM_WINDS_LOCAL)
  MyFloat SmoothedDMveldisp;
#endif
#ifdef BH_BIPOLAR_FEEDBACK
  MyFloat BH_BipolarJ[3];
  MyFloat BH_BipolarColdMass;
  MyFloat BH_BipolarColdFraction;
#endif
#ifdef BH_SPIN_EVOLUTION
  MyFloat GasAngMom[3];
#endif
} data_out;

static data_out *DataResult, *DataOut;



 /* routine to store or combine result data */
static void out2particle(data_out * out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES)      /* initial store */
    {
      BPP(i).BH_NumNgb = out->Ngb;
      BPP(i).BH_Density = out->Rho;
      BPP(i).BH_VolSum = out->VolSum;
      BPP(i).BH_U = out->SmoothedU;
      Dhsmlrho[i] = out->Dhsmlrho;
      BPP(i).BH_SurroundingGasVel[0] = out->GasVel[0];
      BPP(i).BH_SurroundingGasVel[1] = out->GasVel[1];
      BPP(i).BH_SurroundingGasVel[2] = out->GasVel[2];
#ifdef DRAINGAS
      BPP(i).NearestDist = out->NearestDist;
      BPP(i).DrainID = out->DrainID;
      BPP(i).CellDensity = out->CellDensity;
#endif
      BPP(i).BH_DtGasNeighbor = out->DtGasNeighbor;
#ifdef BH_USE_ALFVEN_SPEED_IN_BONDI
      BPP(i).BH_Bpress = out->Bpress;
#endif
#if (defined(GFM_WINDS_VARIABLE) && (GFM_WINDS_VARIABLE==1)) || defined(GFM_WINDS_LOCAL)
      BPP(i).BH_DMVelDisp = out->SmoothedDMveldisp;
#endif
#ifdef BH_BIPOLAR_FEEDBACK
      for(int k=0; k<3; k++)
        {
          BPP(i).BH_BipolarJ[k] = out->BH_BipolarJ[k];
        }
      BPP(i).BH_BipolarColdMass = out->BH_BipolarColdMass;
      BPP(i).BH_BipolarColdFraction = out->BH_BipolarColdFraction;
#endif
#ifdef BH_SPIN_EVOLUTION
      BPP(i).BH_AngMomGasCells[0] = out->GasAngMom[0];
      BPP(i).BH_AngMomGasCells[1] = out->GasAngMom[1];
      BPP(i).BH_AngMomGasCells[2] = out->GasAngMom[2];
#endif

    }
  else                          /* merge */
    {
      BPP(i).BH_NumNgb += out->Ngb;
      BPP(i).BH_Density += out->Rho;
      BPP(i).BH_VolSum += out->VolSum;
      BPP(i).BH_U += out->SmoothedU;
      Dhsmlrho[i] += out->Dhsmlrho;
      BPP(i).BH_SurroundingGasVel[0] += out->GasVel[0];
      BPP(i).BH_SurroundingGasVel[1] += out->GasVel[1];
      BPP(i).BH_SurroundingGasVel[2] += out->GasVel[2];
#ifdef DRAINGAS
      if(out->NearestDist < BPP(i).NearestDist)
        {
          BPP(i).NearestDist = out->NearestDist;
          BPP(i).DrainID = out->DrainID;
          BPP(i).CellDensity = out->CellDensity;
        }
#endif
      if(BPP(i).BH_DtGasNeighbor > out->DtGasNeighbor)
        BPP(i).BH_DtGasNeighbor = out->DtGasNeighbor;
#ifdef BH_USE_ALFVEN_SPEED_IN_BONDI
      BPP(i).BH_Bpress += out->Bpress;
#endif
#if (defined(GFM_WINDS_VARIABLE) && (GFM_WINDS_VARIABLE==1)) || defined(GFM_WINDS_LOCAL)
      BPP(i).BH_DMVelDisp += out->SmoothedDMveldisp;
#endif
#ifdef BH_BIPOLAR_FEEDBACK
      for(int k=0; k<3; k++)
        {
          BPP(i).BH_BipolarJ[k] += out->BH_BipolarJ[k];
        }
      BPP(i).BH_BipolarColdMass += out->BH_BipolarColdMass;
      BPP(i).BH_BipolarColdFraction += out->BH_BipolarColdFraction;
#endif
#ifdef BH_SPIN_EVOLUTION
      BPP(i).BH_AngMomGasCells[0] += out->GasAngMom[0];
      BPP(i).BH_AngMomGasCells[1] += out->GasAngMom[1];
      BPP(i).BH_AngMomGasCells[2] += out->GasAngMom[2];
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

        if(blackhole_isactive(i))
          blackhole_density_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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

        blackhole_density_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

void blackhole_density(void)
{
  if(TimeBinsBHAccretion.GlobalNActiveParticles == 0)
    return;

  TIMER_START(CPU_BH_DENSITY);

  MyFloat *Left, *Right;
  int idx, i, npleft, iter = 0;
  long long ntot;

  Left = (MyFloat *) mymalloc("Left", NumPart * sizeof(MyFloat));
  Right = (MyFloat *) mymalloc("Right", NumPart * sizeof(MyFloat));
  Dhsmlrho = (MyFloat *) mymalloc("Dhsmlrho", NumPart * sizeof(MyFloat));

  for(i = 0; i < NumBHs; i++)
    if(P[BHP[i].PID].ID != 0)
      BHP[i].SwallowID = 0;

  for(idx = 0; idx < TimeBinsBHAccretion.NActiveParticles; idx++)
    {
      i = TimeBinsBHAccretion.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(blackhole_isactive(i))
        {
          Left[i] = Right[i] = 0;
#ifdef DRAINGAS
          BPP(i).NearestDist = MAX_REAL_NUMBER;
          BPP(i).DrainID = 0;
          BPP(i).CellDensity = 0;
#endif
        }
    }

  generic_set_MaxNexport();

  /* we will repeat the whole thing for those black holes where we didn't find enough neighbours */
  do
    {
      generic_comm_pattern(TimeBinsBHAccretion.NActiveParticles, kernel_local, kernel_imported);


      /* do final operations on results */
      npleft = 0;
      for(idx = 0; idx < TimeBinsBHAccretion.NActiveParticles; idx++)
        {
          i = TimeBinsBHAccretion.ActiveParticleList[idx];
          if(i < 0)
            continue;

          if(blackhole_isactive(i))
            {
              if(BPP(i).BH_Density > 0)
                {
                  BPP(i).BH_Pressure = GAMMA_MINUS1 * BPP(i).BH_U / BPP(i).BH_VolSum;
                  BPP(i).BH_U /= BPP(i).BH_Density;
                  BPP(i).BH_SurroundingGasVel[0] /= BPP(i).BH_Density;
                  BPP(i).BH_SurroundingGasVel[1] /= BPP(i).BH_Density;
                  BPP(i).BH_SurroundingGasVel[2] /= BPP(i).BH_Density;

#if (defined(GFM_WINDS_VARIABLE) && (GFM_WINDS_VARIABLE==1)) || defined(GFM_WINDS_LOCAL)
                  BPP(i).BH_DMVelDisp /= BPP(i).BH_Density;
#endif
#ifdef BH_BIPOLAR_FEEDBACK
                  /* We calculate whether there is a cold component present, and normalise the net AM vector */
                  BPP(i).BH_BipolarColdFraction /= BPP(i).BH_BipolarColdMass;
                  if(BPP(i).BH_BipolarColdFraction > All.BHBipolarColdFraction)
                    {
                      printf("BLACK_HOLES: COLD DISK present, fraction: %g \n",BPP(i).BH_BipolarColdFraction);
                      BPP(i).BH_BipolarColdDisk = 1;
                    } else {
                      printf("BLACK_HOLES: NO cold disk present, fraction: %g \n",BPP(i).BH_BipolarColdFraction);
                      BPP(i).BH_BipolarColdDisk = 0;
                    }
                  float jnorm = 0;
                  for (int k=0; k<3; k++) 
                    {
                      jnorm += BPP(i).BH_BipolarJ[k]*BPP(i).BH_BipolarJ[k];
                    }
                  jnorm = sqrt(jnorm);
                  for (int k=0; k<3; k++) 
                    {
                      BPP(i).BH_BipolarJ[k] /= jnorm;
                    }
                  printf("BLACK_HOLES: AM Direction = (%g,%g,%g)\n", BPP(i).BH_BipolarJ[0],BPP(i).BH_BipolarJ[1],BPP(i).BH_BipolarJ[2]);
#endif

                  BPP(i).BH_Density /= BPP(i).BH_VolSum;

#ifdef BH_USE_ALFVEN_SPEED_IN_BONDI
                  BPP(i).BH_Bpress /= BPP(i).BH_VolSum;
#endif
                }

              /* now check whether we had enough neighbours */

              if(BPP(i).BH_NumNgb < (All.DesNumNgbBlackHole - All.MaxNumNgbDeviationBlackHole) || (BPP(i).BH_NumNgb > (All.DesNumNgbBlackHole + All.MaxNumNgbDeviationBlackHole)))
                {
                  /* need to redo this particle */
                  npleft++;

                  if(BPP(i).BH_NumNgb > 0)
                    {
                      Dhsmlrho[i] *= BPP(i).BH_Hsml / (NUMDIMS * BPP(i).BH_NumNgb / (NORM_COEFF * pow(BPP(i).BH_Hsml, 3)));

                      if(Dhsmlrho[i] > -0.9)    /* note: this would be -1 if only a single particle at zero lag is found */
                        Dhsmlrho[i] = 1 / (1 + Dhsmlrho[i]);
                      else
                        Dhsmlrho[i] = 1;
                    }
                  else
                    Dhsmlrho[i] = 1;

                  if(Left[i] > 0 && Right[i] > 0)
                    if((Right[i] - Left[i]) < 1.0e-3 * Left[i] && BPP(i).BH_NumNgb > 0)
                      {
                        /* this one should be ok */
                        npleft--;
                        P[i].TimeBinHydro = -P[i].TimeBinHydro - 1;     /* Mark as inactive */
                        continue;
                      }

                  if(BPP(i).BH_NumNgb < (All.DesNumNgbBlackHole - All.MaxNumNgbDeviationBlackHole))
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
                        ("i=%d task=%d ID=%lld Hsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g\n   pos=(%g|%g|%g)\n",
                         i, ThisTask, (long long) P[i].ID, BPP(i).BH_Hsml, Left[i], Right[i], (double) BPP(i).BH_NumNgb, Right[i] - Left[i], P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
                      myflush(stdout);
                    }

                  if(Right[i] > 0 && Left[i] > 0)
                    BPP(i).BH_Hsml = pow(0.5 * (pow(Left[i], 3) + pow(Right[i], 3)), 1.0 / 3);
                  else
                    {
                      if(Right[i] == 0 && Left[i] == 0)
                        terminate("BLACK_HOLES: Right[i] = Left[i] = 0 in the computation of BH_Hsml"); /* can't occur */

                      if(Right[i] == 0 && Left[i] > 0)
                        {
                          double fac = 1.26;

                          if(fabs(BPP(i).BH_NumNgb - All.DesNumNgbBlackHole) < 0.5 * All.DesNumNgbBlackHole)
                            {
                              fac = 1 - (BPP(i).BH_NumNgb - All.DesNumNgbBlackHole) / (NUMDIMS * BPP(i).BH_NumNgb) * Dhsmlrho[i];

                              if(fac > 1.26)
                                fac = 1.26;
                            }

                          BPP(i).BH_Hsml *= fac;
                        }

                      if(Right[i] > 0 && Left[i] == 0)
                        {
                          double fac = 1 / 1.26;

                          if(fabs(BPP(i).BH_NumNgb - All.DesNumNgbBlackHole) < 0.5 * All.DesNumNgbBlackHole)
                            {
                              fac = 1 - (BPP(i).BH_NumNgb - All.DesNumNgbBlackHole) / (NUMDIMS * BPP(i).BH_NumNgb) * Dhsmlrho[i];

                              if(fac < 1 / 1.26)
                                fac = 1 / 1.26;
                            }
                          BPP(i).BH_Hsml *= fac;
                        }
                    }

                  if(Left[i] > All.BlackHoleMaxAccretionRadius)
                    {
                      /* this will stop the search for a new BH smoothing length in the next iteration */
                      BPP(i).BH_Hsml = Left[i] = Right[i] = All.BlackHoleMaxAccretionRadius;
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
            mpi_printf("BLACK_HOLES: blackhole ngb iteration %d: need to repeat for %lld particles.\n", iter, ntot);

          if(iter > MAXITER)
            terminate("BLACK_HOLES: failed to converge in neighbour iteration in blackhole_density()\n");
        }
    }
  while(ntot > 0);


#ifdef DRAINGAS
  for(idx = 0; idx < TimeBinsBHAccretion.NActiveParticles; idx++)
    {
      i = TimeBinsBHAccretion.ActiveParticleList[idx];
      if(i < 0)
        continue;
      if(blackhole_isactive(i))
        {
      if(BPP(i).NearestDist == MAX_REAL_NUMBER || BPP(i).DrainID == 0)
        terminate
          ("BLACK_HOLES DRAINGAS: no primary cell found: ID=%lld Left=%g Right=%g NearestDist=%g DrainID=%lld rho=%g hsml=%g NumNgb=%g Mdot=%g",
           (long long) P[i].ID, Left[i], Right[i], BPP(i).NearestDist, (long long) BPP(i).DrainID, BPP(i).BH_Density, BPP(i).BH_Hsml, BPP(i).BH_NumNgb, BPP(i).BH_Mdot);
        }
    }
#endif

  myfree(Dhsmlrho);
  myfree(Right);
  myfree(Left);


  /* mark as active again */
  for(idx = 0; idx < TimeBinsBHAccretion.NActiveParticles; idx++)
    {
      i = TimeBinsBHAccretion.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].TimeBinHydro < 0)
        P[i].TimeBinHydro = -P[i].TimeBinHydro - 1;
    }

  TIMER_STOP(CPU_BH_DENSITY);
}


/*! This function represents the core of the blackhole density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
static int blackhole_density_evaluate(int target, int mode, int threadid)
{
  int numnodes, *firstnode;
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


  MyDouble *pos = in->Pos;
  double h = in->BH_Hsml;

  double hinv = 1.0 / h;
  double hinv3 = hinv * hinv * hinv;
  double hinv4 = hinv3 * hinv;

  double dt_neighbor = MAX_REAL_NUMBER;
  double vol_weight = 0;
  double rho = 0;
  double weighted_numngb = 0;
  double dhsmlrho = 0;
  double smoothU = 0;
  double gasvel[3] = { 0, 0, 0 };
#if (defined(GFM_WINDS_VARIABLE) && (GFM_WINDS_VARIABLE==1)) || defined(GFM_WINDS_LOCAL)
  double smooth_DMveldisp = 0;
#endif

#ifdef DRAINGAS
  MyFloat nearestDist = MAX_REAL_NUMBER;
  MyIDType drainID = 0;
  MyFloat cellDensity = 0;
#endif
#ifdef BH_USE_ALFVEN_SPEED_IN_BONDI
  MyFloat Bpress = 0;
#endif

#ifdef BH_BIPOLAR_FEEDBACK
  MyFloat dr[3];
  MyFloat bipolar_j[3];
  bipolar_j[0] = bipolar_j[1] = bipolar_j[2] = 0;
  MyFloat *bh_vel = in->BH_Vel;
  MyFloat total_mass = 0;
  MyFloat cold_fraction = 0;
  double u_to_temp_fac = (4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC))) * PROTONMASS / BOLTZMANN * GAMMA_MINUS1 * All.UnitEnergy_in_cgs / All.UnitMass_in_g;
#endif

#ifdef BH_SPIN_EVOLUTION
  MyFloat gas_angular_momentum[3] =  { 0, 0, 0 };
#endif

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, threadid, numnodes, firstnode);

  for(int n = 0; n < nfound; n++)
    {
      int j = Thread[threadid].Ngblist[n];

      if(P[j].Mass > 0 && P[j].ID != 0)
        {
          double r2 = Thread[threadid].R2list[n];

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

          double mass_j = P[j].Mass;

          rho += (mass_j * wk);

          vol_weight += SphP[j].Volume * wk;

          double ngb_j = mass_j / All.ReferenceGasPartMass;

          if(ngb_j > MAX_CONTRIBUTION_OF_SINGLE_CELL)
            ngb_j = MAX_CONTRIBUTION_OF_SINGLE_CELL;

          weighted_numngb += (NORM_COEFF * ngb_j * wk / hinv3); /* 4.0/3 * PI = 4.188790204786 */
          dhsmlrho += (-ngb_j * (NUMDIMS * hinv * wk + u * dwk));

          smoothU += (mass_j * wk * SphP[j].Utherm);
          gasvel[0] += (mass_j * wk * P[j].Vel[0]);
          gasvel[1] += (mass_j * wk * P[j].Vel[1]);
          gasvel[2] += (mass_j * wk * P[j].Vel[2]);

#ifdef BH_SPIN_EVOLUTION
          gas_angular_momentum[0] += SphP[j].Momentum[1]*(P[j].Pos[2]-pos[2]) - SphP[j].Momentum[2]*(P[j].Pos[1]-pos[1]);
          gas_angular_momentum[1] += SphP[j].Momentum[2]*(P[j].Pos[0]-pos[0]) - SphP[j].Momentum[0]*(P[j].Pos[2]-pos[2]);
          gas_angular_momentum[2] += SphP[j].Momentum[0]*(P[j].Pos[1]-pos[1]) - SphP[j].Momentum[1]*(P[j].Pos[0]-pos[0]);
#endif

#if (defined(GFM_WINDS_VARIABLE) && (GFM_WINDS_VARIABLE==1)) || defined(GFM_WINDS_LOCAL)
          smooth_DMveldisp += (mass_j * wk * SphP[j].w.DMVelDisp);
#endif

#ifdef BH_USE_ALFVEN_SPEED_IN_BONDI
          Bpress += 0.5 * (SphP[j].B[0] * SphP[j].B[0] + SphP[j].B[1] * SphP[j].B[1] + SphP[j].B[2] * SphP[j].B[2]) * SphP[j].Volume * wk;
#endif

#ifdef BH_BIPOLAR_FEEDBACK
          for(int k=0; k<3; k++)
            dr[k] = P[j].Pos[k] - pos[k];

          /* Here we want the mass weighted averaged AM. */
          bipolar_j[0] += mass_j * mass_j * (dr[1] * (SphP[j].Momentum[2]/mass_j - bh_vel[2]) - dr[2] * (SphP[j].Momentum[1]/mass_j - bh_vel[1]));
          bipolar_j[1] += mass_j * mass_j * (dr[2] * (SphP[j].Momentum[0]/mass_j - bh_vel[0]) - dr[0] * (SphP[j].Momentum[2]/mass_j - bh_vel[2]));
          bipolar_j[2] += mass_j * mass_j * (dr[0] * (SphP[j].Momentum[1]/mass_j - bh_vel[1]) - dr[1] * (SphP[j].Momentum[0]/mass_j - bh_vel[0]));

          total_mass += mass_j;
          if(SphP[j].Utherm * u_to_temp_fac < All.BHBipolarColdTemp)
            cold_fraction += mass_j;
#endif

          double csnd = get_sound_speed(j);

          if(SphP[j].Density < 1.0e-10)
            csnd = 0;

#ifdef VORONOI_STATIC_MESH
          csnd += sqrt(P[j].Vel[0] * P[j].Vel[0] + P[j].Vel[1] * P[j].Vel[1] + P[j].Vel[2] * P[j].Vel[2]);
#else
          csnd += sqrt((P[j].Vel[0] - SphP[j].VelVertex[0]) * (P[j].Vel[0] - SphP[j].VelVertex[0]) + (P[j].Vel[1]
                                                                                                      - SphP[j].VelVertex[1]) * (P[j].Vel[1] - SphP[j].VelVertex[1]) + (P[j].Vel[2] -
                                                                                                                                                                        SphP[j].VelVertex[2]) *
                       (P[j].Vel[2] - SphP[j].VelVertex[2]));
#endif

          double rad = pow(SphP[j].Volume / (4.0 / 3 * M_PI), 1.0 / 3);

          if(csnd <= 0)
            csnd = 1.0e-30;

          double dt_courant = rad / csnd;

          if(dt_neighbor > dt_courant)
            dt_neighbor = dt_courant;
#ifdef DRAINGAS
          if((r < nearestDist) && (P[j].Mass != 0) && (P[j].ID != 0))
            {
              drainID = P[j].ID;
              nearestDist = r;
              cellDensity = SphP[j].Density;
            }
#endif
        }
    }


  out.Rho = rho;
  out.Ngb = weighted_numngb;
  out.Dhsmlrho = dhsmlrho;
  out.GasVel[0] = gasvel[0];
  out.GasVel[1] = gasvel[1];
  out.GasVel[2] = gasvel[2];
  out.SmoothedU = smoothU;
  out.DtGasNeighbor = dt_neighbor;
  out.VolSum = vol_weight;
#ifdef DRAINGAS
  out.NearestDist = nearestDist;
  out.DrainID = drainID;
  out.CellDensity = cellDensity;
#endif
#ifdef BH_USE_ALFVEN_SPEED_IN_BONDI
  out.Bpress = Bpress;
#endif
#if (defined(GFM_WINDS_VARIABLE) && (GFM_WINDS_VARIABLE==1)) || defined(GFM_WINDS_LOCAL)
  out.SmoothedDMveldisp = smooth_DMveldisp;
#endif
#ifdef BH_BIPOLAR_FEEDBACK
  for(int k=0; k<3; k++)
    out.BH_BipolarJ[k] = bipolar_j[k];

  out.BH_BipolarColdFraction = cold_fraction;
  out.BH_BipolarColdMass = total_mass;
#endif
#ifdef BH_SPIN_EVOLUTION
  out.GasAngMom[0] = gas_angular_momentum[0];
  out.GasAngMom[1] = gas_angular_momentum[1];
  out.GasAngMom[2] = gas_angular_momentum[2];
#endif

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}


int blackhole_isactive(int n)
{
  if(P[n].TimeBinHydro < 0)
    return 0;

  if(P[n].Type == 5)
    return 1;

  return 0;
}

#endif
