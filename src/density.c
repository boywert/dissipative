/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/density.c
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

#include "allvars.h"
#include "proto.h"
#include "domain.h"


/*! Structure for communication during the density computation. Holds data that is sent to other processors.
 */
static int density_evaluate(int target, int mode, int threadid);


static MyFloat *NumNgb, *DhsmlDensityFactor;
#ifdef FIX_SPH_PARTICLES_AT_IDENTICAL_COORDINATES
static MyFloat *MinDist;
#endif

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble Pos[3];
  MyFloat Hsml;
#ifdef FIX_SPH_PARTICLES_AT_IDENTICAL_COORDINATES
  MyIDType ID;
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
  in->Hsml = SphP[i].Hsml;
#ifdef FIX_SPH_PARTICLES_AT_IDENTICAL_COORDINATES
  in->ID = P[i].ID;
#endif

  in->Firstnode = firstnode;
}

 /* local data structure that holds results acquired on remote processors */
typedef struct
{
  MyFloat Rho;
  MyFloat DhsmlDensity;
  MyFloat Ngb;
#ifdef FIX_SPH_PARTICLES_AT_IDENTICAL_COORDINATES
  MyFloat MinDist;
#endif

#ifdef CONDUCTION_SATURATION    //VITALI
  MyFloat GradEntr[3];
#endif

#ifdef OTVET_FLUXLIMITER
  MyFloat Grad_ngamma[3][OT_N_BINS];
#endif
} data_out;

static data_out *DataResult, *DataOut;



 /* routine to store or combine result data */
static void out2particle(data_out * out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES)      /* initial store */
    {
      NumNgb[i] = out->Ngb;
      if(P[i].Type == 0)
        {
          SphP[i].Density = out->Rho;
          DhsmlDensityFactor[i] = out->DhsmlDensity;
#ifdef FIX_SPH_PARTICLES_AT_IDENTICAL_COORDINATES
          MinDist[i] = out->MinDist;
#endif
#ifdef OTVET_FLUXLIMITER
          int k, ii;
          for(k = 0; k < OT_N_BINS; k++)
            for(ii = 0; ii < 3; ii++)
              SphP[i].Grad_ngamma[ii][k] = out->Grad_ngamma[ii][k];
#endif
        }
    }
  else                          /* combine */
    {
      NumNgb[i] += out->Ngb;
      if(P[i].Type == 0)
        {
          SphP[i].Density += out->Rho;
          DhsmlDensityFactor[i] += out->DhsmlDensity;
#ifdef FIX_SPH_PARTICLES_AT_IDENTICAL_COORDINATES
          if(MinDist[i] > out->MinDist)
            MinDist[i] = out->MinDist;
#endif
#ifdef OTVET_FLUXLIMITER
          int k, ii;
          for(k = 0; k < OT_N_BINS; k++)
            for(ii = 0; ii < 3; ii++)
              SphP[i].Grad_ngamma[ii][k] += out->Grad_ngamma[ii][k];
#endif
        }
    }
}


#include "generic_comm_helpers2.h"

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
              if(generic_polling_primary(count, TimeBinsHydro.NActiveParticles))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        idx = NextParticle++;

        if(idx >= TimeBinsHydro.NActiveParticles)
          break;

        int i = TimeBinsHydro.ActiveParticleList[idx];
        if(i < 0)
          continue;

        if(density_isactive(i))
          density_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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

        density_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}




/*! \file density.c
 *  \brief SPH density computation and smoothing length determination
 *
 *  This file contains the "first SPH loop", where the SPH densities and some
 *  auxiliary quantities are computed.  There is also functionality that
 *  corrects the smoothing length if needed.
 */


/*! This function computes the local density for each active SPH particle, the
 * number of neighbours in the current smoothing radius, and the divergence
 * and rotation of the velocity field.  The pressure is updated as well.  If a
 * particle with its smoothing region is fully inside the local domain, it is
 * not exported to the other processors. The function also detects particles
 * that have a number of neighbours outside the allowed tolerance range. For
 * these particles, the smoothing length is adjusted accordingly, and the
 * density() computation is called again.  Note that the smoothing length is
 * not allowed to fall below the lower bound set by MinGasHsml (this may mean
 * that one has to deal with substantially more than normal number of
 * neighbours.)
 */

static MyFloat *NumNgb, *DhsmlDensityFactor;
#ifdef FIX_SPH_PARTICLES_AT_IDENTICAL_COORDINATES
static MyFloat *MinDist;
#endif

void density(void)
{
  MyFloat *Left, *Right;
  int idx, i, npleft, iter = 0;
  long long ntot;
  double desnumngb, t0, t1;
#ifdef OTVET_FLUXLIMITER
  int k, ii;
#endif

  CPU_Step[CPU_MISC] += measure_time();

  NumNgb = (MyFloat *) mymalloc("NumNgb", NumPart * sizeof(MyFloat));
  DhsmlDensityFactor = (MyFloat *) mymalloc("DhsmlDensityFactor", NumPart * sizeof(MyFloat));
  Left = (MyFloat *) mymalloc("Left", NumPart * sizeof(MyFloat));
  Right = (MyFloat *) mymalloc("Right", NumPart * sizeof(MyFloat));

#ifdef FIX_SPH_PARTICLES_AT_IDENTICAL_COORDINATES
  MinDist = (MyFloat *) mymalloc("MinDist", NumPart * sizeof(MyFloat));
#endif

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(density_isactive(i))
        {
          Left[i] = Right[i] = 0;
        }
    }

  generic_set_MaxNexport();

  desnumngb = All.DesNumNgb;



  /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
  do
    {
      t0 = second();

      generic_comm_pattern(TimeBinsHydro.NActiveParticles, kernel_local, kernel_imported);

      /* do final operations on results */
      for(idx = 0, npleft = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
        {
          i = TimeBinsHydro.ActiveParticleList[idx];
          if(i < 0)
            continue;

          if(density_isactive(i))
            {
              if(P[i].Type == 0)
                {
                  if(SphP[i].Density > 0)
                    {
                      DhsmlDensityFactor[i] *= SphP[i].Hsml / (NUMDIMS * SphP[i].Density);
                      if(DhsmlDensityFactor[i] > -0.9)  /* note: this would be -1 if only a single particle at zero lag is found */
                        DhsmlDensityFactor[i] = 1 / (1 + DhsmlDensityFactor[i]);
                      else
                        DhsmlDensityFactor[i] = 1;

#ifdef OTVET_FLUXLIMITER
                      for(k = 0; k < OT_N_BINS; k++)
                        {
                          SphP[i].Grad_ngamma[0][k] /= SphP[i].Density;
                          SphP[i].Grad_ngamma[1][k] /= SphP[i].Density;
                          SphP[i].Grad_ngamma[2][k] /= SphP[i].Density;
                        }
#endif
                    }
                }

              if(NumNgb[i] < (desnumngb - All.MaxNumNgbDeviation) || NumNgb[i] > (desnumngb + All.MaxNumNgbDeviation))
                {
                  /* need to redo this particle */
                  npleft++;

                  if(Left[i] > 0 && Right[i] > 0)
                    if((Right[i] - Left[i]) < 1.0e-3 * Left[i])
                      {
                        /* this one should be ok */
                        npleft--;
                        P[i].TimeBinHydro = -P[i].TimeBinHydro - 1;     /* Mark as inactive */
                        continue;
                      }

                  if(NumNgb[i] < (desnumngb - All.MaxNumNgbDeviation))
                    Left[i] = dmax(SphP[i].Hsml, Left[i]);
                  else
                    {
                      if(Right[i] != 0)
                        {
                          if(SphP[i].Hsml < Right[i])
                            Right[i] = SphP[i].Hsml;
                        }
                      else
                        Right[i] = SphP[i].Hsml;
                    }

                  if(iter >= MAXITER - 10)
                    {
                      printf("i=%d task=%d ID=%d Hsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g\n   pos=(%g|%g|%g)\n",
                             i, ThisTask, (int) P[i].ID, SphP[i].Hsml, Left[i], Right[i], (float) NumNgb[i], Right[i] - Left[i], P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
                      myflush(stdout);
                    }

                  if(Right[i] > 0 && Left[i] > 0)
                    SphP[i].Hsml = pow(0.5 * (pow(Left[i], 3) + pow(Right[i], 3)), 1.0 / 3);
                  else
                    {
                      if(Right[i] == 0 && Left[i] == 0)
                        terminate("should not occur");

                      if(Right[i] == 0 && Left[i] > 0)
                        {
                          if(P[i].Type == 0 && fabs(NumNgb[i] - desnumngb) < 0.5 * desnumngb)
                            {
                              /*
                                 fac = 1 - (SphP[i].NumNgb -
                                 desnumngb) / (NUMDIMS * SphP[i].NumNgb) * SphP[i].DhsmlDensityFactor;

                                 if(fac < 1.26)
                                 SphP[i].Hsml *= fac;
                                 else
                               */

                              SphP[i].Hsml *= 1.26;
                            }
                          else
                            SphP[i].Hsml *= 1.26;
                        }

                      if(Right[i] > 0 && Left[i] == 0)
                        {
                          if(P[i].Type == 0 && fabs(NumNgb[i] - desnumngb) < 0.5 * desnumngb)
                            {
                              /*
                                 fac = 1 - (SphP[i].NumNgb -
                                 desnumngb) / (NUMDIMS * SphP[i].NumNgb) * SphP[i].DhsmlDensityFactor;

                                 if(fac > 1 / 1.26)
                                 SphP[i].Hsml *= fac;
                                 else
                                 SphP[i].Hsml /= 1.26;
                               */

                              SphP[i].Hsml /= 1.26;
                            }
                          else
                            SphP[i].Hsml /= 1.26;
                        }
                    }

                }
              else
                P[i].TimeBinHydro = -P[i].TimeBinHydro - 1;     /* Mark as inactive */
            }
        }

      sumup_large_ints(1, &npleft, &ntot);

      t1 = second();

      if(ntot > 0)
        {
          iter++;

          if(iter > 0)
            mpi_printf("DENSITY: ngb iteration %3d: need to repeat for %12lld particles. (took %g sec)\n", iter, ntot, timediff(t0, t1));

          if(iter > MAXITER)
            terminate("failed to converge in neighbour iteration in density()\n");
        }
    }
  while(ntot > 0);


#ifdef FIX_SPH_PARTICLES_AT_IDENTICAL_COORDINATES

#if defined(REFLECTIVE_X) && defined(REFLECTIVE_Y) && defined(REFLECTIVE_Z)

  int count2 = 0;
  int countall2 = 0;

  for(i = 0; i < NumGas; i++)
    {
      /*
       * If the distance to the border of a particle is too small,
       * then the ghost particle will be too close to this particle.
       * Therefore we shift the particle in this case into the direction of the box center.
       */
      if(distance_to_border(i) < 0.5 * 0.001 * SphP[i].Hsml)
        {
          count2++;

          double dir[3];

          dir[0] = boxSize_X * 0.5 - P[i].Pos[0];
          dir[1] = boxSize_Y * 0.5 - P[i].Pos[1];
          dir[2] = boxSize_Z * 0.5 - P[i].Pos[2];

          double n = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
          //note: it's not possible that the operand of sqrt is zero here.

          dir[0] /= n;
          dir[1] /= n;
          dir[2] /= n;

          P[i].Pos[0] += 0.05 * SphP[i].Hsml * dir[0];
          P[i].Pos[1] += 0.05 * SphP[i].Hsml * dir[1];
          P[i].Pos[2] += 0.05 * SphP[i].Hsml * dir[2];
        }
    }

  MPI_Allreduce(&count2, &countall2, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  mpi_printf("\nFOUND %d particles extremely close to the reflective boundary. Fixing this. \n\n", countall2);
#endif

  int count = 0, countall;

  for(i = 0; i < NumGas; i++)
    if(MinDist[i] < 0.001 * SphP[i].Hsml)
      count++;

  MPI_Allreduce(&count, &countall, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(countall)
    {
      mpi_printf("\nFOUND %d SPH particles with an extremely close neighbor. Fixing this. \n\n", countall);

      for(i = 0; i < NumGas; i++)
        if(MinDist[i] < 0.001 * SphP[i].Hsml)
          {
            double theta = acos(2 * get_random_number() - 1);
            double phi = 2 * M_PI * get_random_number();

            P[i].Pos[0] += 0.1 * SphP[i].Hsml * sin(theta) * cos(phi);
            P[i].Pos[1] += 0.1 * SphP[i].Hsml * sin(theta) * sin(phi);
            P[i].Pos[2] += 0.1 * SphP[i].Hsml * cos(theta);
          }
    }

#endif



#ifdef FIX_SPH_PARTICLES_AT_IDENTICAL_COORDINATES
  myfree(MinDist);
#endif
  myfree(Right);
  myfree(Left);
  myfree(DhsmlDensityFactor);
  myfree(NumNgb);


  /* mark as active again */
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].TimeBinHydro < 0)
        P[i].TimeBinHydro = -P[i].TimeBinHydro - 1;
    }

  /* collect some timing information */

  CPU_Step[CPU_INIT] += measure_time();

}


/*! This function represents the core of the SPH density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
static int density_evaluate(int target, int mode, int threadid)
{
  int j, n;
  int numngb, numnodes, *firstnode;
  double h, h2, hinv, hinv3, hinv4;
  MyFloat rho;
  double wk, dwk;
  double dx, dy, dz, r, r2, u, mass_j;
  MyFloat weighted_numngb;
  MyFloat dhsmlrho;
  MyDouble *pos;
#ifdef OTVET_FLUXLIMITER
  int k;
  double mj_dwk_r;
  MyFloat grad_ngamma[3][OT_N_BINS];
#endif

#ifdef FIX_SPH_PARTICLES_AT_IDENTICAL_COORDINATES
  MyFloat mindist = MAX_REAL_NUMBER;
  MyIDType ID;
#endif
#ifdef OTVET_FLUXLIMITER
  for(k = 0; k < OT_N_BINS; k++)
    for(j = 0; j < 3; j++)
      grad_ngamma[j][k] = 0;
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
  h = target_data->Hsml;
#ifdef FIX_SPH_PARTICLES_AT_IDENTICAL_COORDINATES
  ID = target_data->ID;
#endif

  h2 = h * h;
  hinv = 1.0 / h;
#ifndef  TWODIMS
  hinv3 = hinv * hinv * hinv;
#else
  hinv3 = hinv * hinv / boxSize_Z;
#endif
  hinv4 = hinv3 * hinv;

  numngb = 0;
  rho = weighted_numngb = dhsmlrho = 0;

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, threadid, numnodes, firstnode);

  for(n = 0; n < nfound; n++)
    {
      j = Thread[threadid].Ngblist[n];

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

          mass_j = P[j].Mass;

          rho += FLT(mass_j * wk);

          weighted_numngb += FLT(NORM_COEFF * wk / hinv3);      /* 4.0/3 * PI = 4.188790204786 */

          dhsmlrho += FLT(-mass_j * (NUMDIMS * hinv * wk + u * dwk));

#ifdef FIX_SPH_PARTICLES_AT_IDENTICAL_COORDINATES
          if(ID != P[j].ID && mindist > r)
            mindist = r;
#endif

#ifdef OTVET_FLUXLIMITER
          if(r > 0)
            mj_dwk_r = mass_j * dwk / r;
          for(k = 0; k < OT_N_BINS; k++)
            {
              grad_ngamma[0][k] += mj_dwk_r * dx * SphP[j].n_gamma[k];
              grad_ngamma[1][k] += mj_dwk_r * dy * SphP[j].n_gamma[k];
              grad_ngamma[2][k] += mj_dwk_r * dz * SphP[j].n_gamma[k];
            }
#endif
        }
    }

  out.Rho = rho;
  out.Ngb = weighted_numngb;
  out.DhsmlDensity = dhsmlrho;
#ifdef FIX_SPH_PARTICLES_AT_IDENTICAL_COORDINATES
  out.MinDist = mindist;
#endif
#ifdef OTVET_FLUXLIMITER
  for(k = 0; k < OT_N_BINS; k++)
    for(j = 0; j < 3; j++)
      out.Grad_ngamma[j][k] = grad_ngamma[j][k];
#endif

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}



int density_isactive(int n)
{
  if(P[n].TimeBinHydro < 0)
    return 0;

  if(P[n].Type == 0)
    return 1;

  return 0;
}
