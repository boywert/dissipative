/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/powerspec_vel.c
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
#include <string.h>
#include <gsl/gsl_errno.h>

#include "allvars.h"
#include "proto.h"

#if defined(PMGRID) && defined(PERIODIC) && defined(VEL_POWERSPEC)

#define  PMGRIDz (PMGRID/2 + 1)
#define  PMGRID2 (2*PMGRIDz)

#if (PMGRID > 1024)
typedef long long large_array_offset;
#else
typedef unsigned int large_array_offset;
#endif

static double *powerspec_vel_nearest_distance, *powerspec_vel_nearest_hsml;

#define BINS_PS  2000           /* number of bins for power spectrum computation */
struct powerspectrum {
  double K0, K1;
  double binfac;
  long long CountModes[BINS_PS];
  double SumPower[BINS_PS];
  double Power[BINS_PS];
  double Kbin[BINS_PS];
  double LostPower;
};

static void powerspec_vel_calc_and_bin_spectrum( fft_real *field, struct powerspectrum *ps, fft_plan *myplan, int flag);
static void powerspec_vel_collect_and_save( struct powerspectrum *ps, char *name, int RestartSnapNum );
static double powerspec_vel_obtain_fields( fft_plan *myplan );
static int powerspec_vel_find_nearest_evaluate(int target, int mode, int thread_id);

typedef struct data_in
{
  MyDouble Pos[3];
  double Hsml;
  int Firstnode;
} data_in;

static double CMVel[3];

static data_in *DataIn, *DataGet;

#ifdef MHD
static fft_real *bfield[3];
#endif
static fft_real *velfield[3];
static fft_real *velrhofield[3];
static fft_real *densityfield, *workspace;

static fft_plan myplan;

static void particle2in(struct data_in *in, int place, int firstnode)
{
  int xx = place / (PMGRID * PMGRID);
  int yy = (place - xx * PMGRID * PMGRID) / PMGRID;
  int zz = (place - xx * PMGRID * PMGRID - yy * PMGRID);
  xx += myplan.slabstart_x;

#ifndef VEL_POWERSPEC_BOX
  double x = (xx + 0.5) / PMGRID * All.BoxSize;
  double y = (yy + 0.5) / PMGRID * All.BoxSize;
  double z = (zz + 0.5) / PMGRID * All.BoxSize;
#else
  double x = All.PSBoxMinX + (xx + 0.5) * All.PSCellSize;
  double y = All.PSBoxMinY + (yy + 0.5) * All.PSCellSize;
  double z = All.PSBoxMinZ + (zz + 0.5) * All.PSCellSize;
#endif

  in->Pos[0] = x;
  in->Pos[1] = y;
  in->Pos[2] = z;
  in->Hsml = powerspec_vel_nearest_hsml[place];

  in->Firstnode = firstnode;
}

typedef struct data_out
{
  double Distance;
#ifdef MHD
  MyDouble B[3];
#endif
  MyDouble Vel[3];
  MyDouble Density;
} data_out;

static data_out *DataResult, *DataOut;

static void out2particle(struct data_out *out, int place, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES)      /* initial store */
    {
      powerspec_vel_nearest_distance[place] = sqrt(out->Distance);

      int i = place / (PMGRID * PMGRID);
      int j = (place - i * PMGRID * PMGRID) / PMGRID;
      int k = (place - i * PMGRID * PMGRID - j * PMGRID);
      int ip = PMGRID2 * (PMGRID * i + j) + k;

#ifdef MHD
      for(int m = 0; m < 3; m++)
        bfield[m][ip] = out->B[m];
#endif	

      for(int m = 0; m < 3; m++)
        velfield[m][ip] = out->Vel[m];

      for(int m = 0; m < 3; m++)
        velrhofield[m][ip] = sqrt(out->Density) * out->Vel[m];
    }
  else                          /* combine */
    {
      if(out->Distance < powerspec_vel_nearest_distance[place])
        {
          powerspec_vel_nearest_distance[place] = out->Distance;

          int ii = place / (PMGRID * PMGRID);
          int jj = (place - ii * PMGRID * PMGRID) / PMGRID;
          int kk = (place - ii * PMGRID * PMGRID - jj * PMGRID);
          int ip = PMGRID2 * (PMGRID * ii + jj) + kk;

#ifdef MHD
          for(int m=0; m < 3; m++)
            bfield[m][ip] = out->B[m];
#endif
          for(int m=0; m < 3; m++)
            velfield[m][ip] = out->Vel[m];
          
          for(int m=0; m < 3; m++)
            velrhofield[m][ip] = sqrt(out->Density) * out->Vel[m];
        }
    }
}

#include "generic_comm_helpers2.h"

static large_array_offset Ncount;

static void kernel_local(void)
{
  large_array_offset i;
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
              if(generic_polling_primary(count, Ncount))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        i = NextParticle++;

        if(i >= Ncount)
          break;

        if(powerspec_vel_nearest_distance[i] > 1.0e29)
          powerspec_vel_find_nearest_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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

        powerspec_vel_find_nearest_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

void powerspec_vel( int RestartSnapNum )
{
#ifdef MHD
  mpi_printf("Starting velocity and magnetic field  computation, t=%g.\n", All.Time );
#else
  mpi_printf("Starting velocity powerspec computation, t=%g.\n", All.Time );
#endif

#ifdef VEL_POWERSPEC_BOX
  All.PSBoxMinX = All.PSCenterX - 2. * All.PSRadius;
  All.PSBoxMinY = All.PSCenterY - 2. * All.PSRadius;
  All.PSBoxMinZ = All.PSCenterZ - 2. * All.PSRadius;
  All.PSCellSize = 4. * All.PSRadius / PMGRID;
#endif

  double tstart = second();

  densityfield = (fft_real *) mymalloc("densityfield", PMGRID2 * sizeof(fft_real));
  workspace = (fft_real *) mymalloc("workspace", PMGRID2 * sizeof(fft_real));
  
#ifdef DOUBLEPRECISION_FFTW
  int alignflag = 0;
#else
  /* for single precision, the start of our FFT columns is presently only guaranteed to be 8-byte aligned */
  int alignflag = FFTW_UNALIGNED;
#endif

  int stride = PMGRIDz;

  int ndim[1] = { PMGRID };     /* dimension of the 1D transforms */
  /* Set up the FFTW plan  */
  myplan.forward_plan_zdir = FFTW(plan_many_dft_r2c) (1, ndim, 1, densityfield, 0, 1, PMGRID2, (fft_complex *) workspace, 0, 1, PMGRIDz, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  myplan.forward_plan_ydir = FFTW(plan_many_dft) (1, ndim, 1, (fft_complex *) densityfield, 0, stride, PMGRIDz * PMGRID, (fft_complex *) workspace, 0, stride, PMGRIDz * PMGRID, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  myplan.forward_plan_xdir = FFTW(plan_many_dft) (1, ndim, 1, (fft_complex *) densityfield, 0, stride, PMGRIDz * PMGRID, (fft_complex *) workspace, 0, stride, PMGRIDz * PMGRID, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  myfree(workspace);
  myfree(densityfield);

  my_slab_based_fft_init(&myplan, PMGRID, PMGRID, PMGRID);

  size_t maxfftsize = imax(myplan.largest_x_slab * PMGRID, myplan.largest_y_slab * PMGRID) * ((size_t) PMGRID2);
  
  /* allocate the memory to hold the FFT fields */
  for(int m=0; m<3; m++) {
    velfield[m] = (fft_real *) mymalloc("velgrid", maxfftsize * sizeof(fft_real));
    memset(velfield[m], 0, maxfftsize * sizeof(fft_real));

    velrhofield[m] = (fft_real *) mymalloc("velrhogrid", maxfftsize * sizeof(fft_real));
    memset(velrhofield[m], 0, maxfftsize * sizeof(fft_real));

#ifdef MHD
    bfield[m] = (fft_real *) mymalloc("bgrid", maxfftsize * sizeof(fft_real));
    memset(bfield[m], 0, maxfftsize * sizeof(fft_real));
#endif
  }

  workspace = (fft_real *) mymalloc("workspace", maxfftsize * sizeof(fft_real));
  
  double tdiff = powerspec_vel_obtain_fields( &myplan );
  mpi_printf( "Mapping to grid took %g seconds.\n", tdiff );
  
  struct powerspectrum ps;
  
#ifdef VEL_POWERSPEC_BOX
  ps.K0 = 2 * M_PI / (4.*All.PSRadius);
#else
  ps.K0 = 2 * M_PI / All.BoxSize;  /* minimum k */
#endif
  
  ps.K1 = ps.K0 * PMGRID / 2;         /* maximum k */
  ps.binfac = BINS_PS / (log(ps.K1) - log(ps.K0));
  
  for(int m=0; m<3; m++)
    powerspec_vel_calc_and_bin_spectrum(velfield[m], &ps, &myplan, m == 0);
  powerspec_vel_collect_and_save( &ps, "vel", RestartSnapNum );

  for(int m=0; m<3; m++)
    powerspec_vel_calc_and_bin_spectrum(velrhofield[m], &ps, &myplan, m == 0);
  powerspec_vel_collect_and_save( &ps, "velrho", RestartSnapNum );

  for(int m=0; m<3; m++)
    powerspec_vel_calc_and_bin_spectrum(bfield[m], &ps, &myplan, m == 0);
  powerspec_vel_collect_and_save( &ps, "B", RestartSnapNum );

  myfree(workspace);
  for(int m=2; m>=0; m--) {
    myfree(bfield[m]);
    myfree(velrhofield[m]);
    myfree(velfield[m]);
  }
  
  double tend = second();

  if(ThisTask == 0)
    {
      printf("End power spectrum, took %g seconds\n", timediff(tstart, tend));
      myflush(stdout);
    }
}

void powerspec_vel_calc_and_bin_spectrum( fft_real *field, struct powerspectrum *ps, fft_plan *myplan, int flag)
{
  my_slab_based_fft(myplan, &field[0], &workspace[0], 1);
  
  fft_complex *fft_of_field = (fft_complex *) & field[0];
  
  if(flag) {
    for(int i = 0; i < BINS_PS; i++) {
      ps->SumPower[i] = 0;
      ps->CountModes[i] = 0;
    }
    ps->LostPower = 0;
  }
  
#pragma omp parallel for private(y, z)
  for(int x = 0; x < PMGRID; x++)
    for(int y = myplan->slabstart_y; y < myplan->slabstart_y + myplan->nslab_y; y++)
      for(int z = 0; z < PMGRID; z++)
        {
          int zz = z;
          if(z >= PMGRID / 2 + 1)
            zz = PMGRID - z;

          double kx, ky, kz;
          if(x > PMGRID / 2)
            kx = x - PMGRID;
          else
            kx = x;
          if(y > PMGRID / 2)
            ky = y - PMGRID;
          else
            ky = y;
          if(z > PMGRID / 2)
            kz = z - PMGRID;
          else
            kz = z;

          double k2 = kx * kx + ky * ky + kz * kz;

          large_array_offset ip = ((large_array_offset) PMGRIDz) * (PMGRID * (y - myplan->slabstart_y) + x) + zz;

          double po = (fft_of_field[ip][0] * fft_of_field[ip][0] + fft_of_field[ip][1] * fft_of_field[ip][1]) / pow(PMGRID, 6);

          if(k2 > 0)
            {
              if(k2 < (PMGRID / 2.0) * (PMGRID / 2.0))
                {
                  double k = sqrt(k2) * ps->K0;

                  if(k >= ps->K0 && k < ps->K1)
                    {
                      int bin = log(k / ps->K0) * ps->binfac;

                      ps->SumPower[bin] += po;

                      if(flag)
                        ps->CountModes[bin] += 1;
                    }
                  else
                    ps->LostPower += po;
                }
            }
        }
}

void powerspec_vel_collect_and_save( struct powerspectrum *ps, char *name, int RestartSnapNum )
{
  long long int *countbuf = (long long int *) mymalloc("countbuf", NTask * BINS_PS * sizeof(long long));
  double *powerbuf = (double *) mymalloc("powerbuf", NTask * BINS_PS * sizeof(double));

  MPI_Allgather(ps->CountModes, BINS_PS * sizeof(long long), MPI_BYTE, countbuf, BINS_PS * sizeof(long long), MPI_BYTE, MPI_COMM_WORLD);

  for(int i = 0; i < BINS_PS; i++)
    {
      ps->CountModes[i] = 0;
      for(int n = 0; n < NTask; n++)
        ps->CountModes[i] += countbuf[n * BINS_PS + i];
    }

  MPI_Allgather(ps->SumPower, BINS_PS * sizeof(double), MPI_BYTE, powerbuf, BINS_PS * sizeof(double), MPI_BYTE, MPI_COMM_WORLD);

  for(int i = 0; i < BINS_PS; i++)
    {
      ps->SumPower[i] = 0;
      for(int n = 0; n < NTask; n++)
        ps->SumPower[i] += powerbuf[n * BINS_PS + i];
    }

  myfree(powerbuf);
  myfree(countbuf);

  double prefac = 0.5 * pow( 2. * M_PI / ps->K0, 3. );
  double totPower = 0;
  for(int i = 0; i < BINS_PS; i++)
    {
      ps->Kbin[i] = exp((i + 0.5) / ps->binfac + log(ps->K0));
      /*
      if(ps->CountModes[i] > 0)
        ps->Power[i] = ps->SumPower[i] / ps->CountModes[i];
      else
        ps->Power[i] = 0;
      */

      double KMin = exp((i + 0.) / ps->binfac + log(ps->K0));
      double KMax = exp((i + 1.) / ps->binfac + log(ps->K0));
      
      ps->Power[i] = prefac * ps->SumPower[i] / (KMax - KMin); /* total energy per 1D binsize */
      totPower += prefac * ps->SumPower[i];
    }

  double LostPower;
  ps->LostPower *= prefac;
  MPI_Reduce(&ps->LostPower, &LostPower, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  mpi_printf( "Total energy in spectrum %s: %g (code units), K0=%g, LostPower=%g, sum=%g.\n", name, totPower, ps->K0, LostPower, totPower+LostPower );
  
  if(ThisTask == 0)
    {
      FILE *fd;
      char fname[500], buf[500];
      
      sprintf(fname, "%s/powerspec_%s_%03d.txt", All.OutputDir, name, RestartSnapNum);
      if(!(fd = fopen(fname, "w")))
        {
          sprintf(buf, "can't open file `%s`\n", fname);
          terminate(buf);
        }

      fprintf(fd, "%g\n", All.Time);
      fprintf(fd, "%g\n", ps->K0 );
      fprintf(fd, "%g\n", ps->K1 );
      fprintf(fd, "%d\n", (int)(PMGRID));
      fprintf(fd, "%d\n", (int)(BINS_PS));

      for(int i = 0; i < BINS_PS; i++)
        {
          fprintf(fd, "%g %g %g %g\n", ps->Kbin[i], ps->Power[i], (double)ps->CountModes[i], ps->SumPower[i]);
        }

      fclose(fd);
    }
}

/* this function determines the velocity fields by using the nearest cell's values 
 */
double powerspec_vel_obtain_fields( fft_plan *myplan )
{
  double tstart = second();

  mpi_printf("Start finding nearest gas-particle for mesh-cell centers (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

#ifdef VEL_POWERSPEC_BOX
  mpi_printf("Limiting to radius=%g (code units), %g (physical units).\n", All.PSRadius, All.PSRadius * All.cf_atime / All.HubbleParam );
#endif

  Ncount = ((large_array_offset) myplan->nslab_x) * (PMGRID * PMGRID); /* number of grid points on the local slab */

  powerspec_vel_nearest_distance = (double *) mymalloc("powerspec_vel_nearest_distance", sizeof(double) * Ncount);
  powerspec_vel_nearest_hsml = (double *) mymalloc("powerspec_vel_nearest_hsml", sizeof(double) * Ncount);

  double mass = 0;
  double cmvel[3];
  for(int m=0; m<3; m++)
    cmvel[m] = 0;

  int gasCount = 0;
  for(int n=0; n<NumGas; n++)
    {
      double dx = P[n].Pos[0] - All.PSCenterX;
      double dy = P[n].Pos[1] - All.PSCenterY;
      double dz = P[n].Pos[2] - All.PSCenterZ;

      double r2 = dx*dx + dy*dy + dz*dz;
      if(r2 < All.PSRadius*All.PSRadius)
        {
          gasCount++;
          
          mass += P[n].Mass;
          for(int m=0; m<3; m++)
            cmvel[m] += P[n].Mass * P[n].Vel[m];
        }
    }
  
  double massAll;
  MPI_Allreduce(&mass, &massAll, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(cmvel, CMVel, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for(int m=0; m<3; m++)
    CMVel[m] /= massAll;

  mpi_printf( "CMVel: %g, %g, %g\n", CMVel[0], CMVel[1], CMVel[2] );

  double eKin = 0;
#ifdef MHD
  double eMag = 0;
#endif

  for(int n=0; n<NumGas; n++)
    {
      double dx = P[n].Pos[0] - All.PSCenterX;
      double dy = P[n].Pos[1] - All.PSCenterY;
      double dz = P[n].Pos[2] - All.PSCenterZ;

      double r2 = dx*dx + dy*dy + dz*dz;
      if(r2 < All.PSRadius*All.PSRadius)
        {
          double vx = P[n].Vel[0] - CMVel[0];
          double vy = P[n].Vel[1] - CMVel[1];
          double vz = P[n].Vel[2] - CMVel[2];
          eKin += 0.5 * P[n].Mass * (vx*vx + vy*vy + vz*vz);
#ifdef MHD
          eMag += 0.5 * SphP[n].Volume * (SphP[n].B[0]*SphP[n].B[0] + SphP[n].B[1]*SphP[n].B[1] + SphP[n].B[1]*SphP[n].B[1]);
#endif
        }
    }
  
  double eKinAll;
  MPI_Reduce(&eKin, &eKinAll, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  mpi_printf( "Total kinetic energy in radius: %g (code units), %g (physical units).\n", eKinAll, eKinAll * All.cf_a2inv / All.HubbleParam );
  
#ifdef MHD
  double eMagAll;
  MPI_Reduce(&eMag, &eMagAll, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  mpi_printf( "Total magnetic energy in radius: %g (code units), %g (physical units).\n", eMagAll, eMagAll / All.cf_atime / All.HubbleParam );
#endif
  
  long long gasCountAll;
  sumup_large_ints(1, &gasCount, &gasCountAll);

#ifdef VEL_POWERSPEC_BOX
  double hsml = 4. * All.PSRadius / pow(gasCountAll, 1.0 / 3);
#else
  double hsml = All.BoxSize / pow(All.TotNumGas, 1.0 / 3);
#endif

  mpi_printf( "Initial hsml=%g (code units), %g (physical), cells=%lld.\n", hsml, hsml * All.cf_atime / All.HubbleParam, gasCountAll );

  for(large_array_offset n = 0; n < Ncount; n++)
    {
      powerspec_vel_nearest_distance[n] = 1.0e30;
      powerspec_vel_nearest_hsml[n] = hsml;
    }

  generic_set_MaxNexport();

  long long ntot;
  int iter = 0;
  do {
    generic_comm_pattern(Ncount, kernel_local, kernel_imported);

    long long npleft;
    /* do final operations on results */
    for(large_array_offset i = 0, npleft = 0; i < Ncount; i++)
      {
        if(powerspec_vel_nearest_distance[i] > 1.0e29)
          {
            /* need to redo this particle */
            npleft++;
            powerspec_vel_nearest_hsml[i] *= 2.0;
          }
         else
          {
            powerspec_vel_nearest_distance[i] = 0;    /* we not continue to search for this particle */
          }
       }
    sumup_longs(1, &npleft, &ntot);
    if(ntot > 0)
      {
        iter++;
        if(iter > 0 && ThisTask == 0)
          {
            printf("powerpec_vel nearest iteration %d: need to repeat for %lld particles.\n", iter, ntot);
            myflush(stdout);
          }

        if(iter > MAXITER)
          terminate("failed to converge");
      }
    }
  while(ntot > 0);

  myfree(powerspec_vel_nearest_hsml);
  myfree(powerspec_vel_nearest_distance);
  
#ifdef VEL_POWERSPEC_BOX
  double cellsize = All.PSCellSize;
#else
  double cellsize = All.BoxSize / PMGRID;
#endif
  double vol = cellsize * cellsize * cellsize;
  
  double ekin = 0;
  for(large_array_offset n = 0; n < Ncount; n++)
    {
      for(int m=0; m < 3; m++)
        ekin += velrhofield[m][n] * velrhofield[m][n];
    }
  
  ekin *= 0.5 * vol;
  double ekin_all;
  MPI_Reduce(&ekin, &ekin_all, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
  mpi_printf( "Total kinetic energy on grid: %g (code units), %g (physical units).\n", ekin_all, ekin_all * All.cf_a2inv / All.HubbleParam );

#ifdef MHD
  double emag = 0;
  for(large_array_offset n = 0; n < Ncount; n++)
    {
      for(int m=0; m < 3; m++)
        emag += bfield[m][n] * bfield[m][n];
    }
  
  emag *= 0.5 * vol;
  double emag_all;
  MPI_Reduce(&emag, &emag_all, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
  mpi_printf( "Total magnetic energy on grid: %g (code units), %g (physical units).\n", emag_all, emag_all / All.cf_atime / All.HubbleParam );
#endif

  mpi_printf("Done computing fields.\n");

  double tend = second();
  return timediff(tstart, tend);
}

int powerspec_vel_find_nearest_evaluate(int target, int mode, int thread_id)
{
  int numnodes, *firstnode;
  double xtmp, ytmp, ztmp;

  struct data_in local, *target_data;
  struct data_out out;

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
  
#ifdef VEL_POWERSPEC_BOX
  if(mode == MODE_LOCAL_PARTICLES){
    double dx = NGB_PERIODIC_LONG_X(pos[0] - All.PSCenterX);
    double dy = NGB_PERIODIC_LONG_Y(pos[1] - All.PSCenterY);
    double dz = NGB_PERIODIC_LONG_Z(pos[2] - All.PSCenterZ);
    
    double r2 = dx*dx + dy*dy + dz*dz;
    if(r2 > All.PSRadius*All.PSRadius)
      {
        out.Distance = 0;
#ifdef MHD
        for(int m = 0; m < 3; m++)
          out.B[m] = 0;
#endif
        for(int m = 0; m < 3; m++)
          out.Vel[m] = 0;
        out.Density = 0;
        
        return 0;
      }
  }
#endif

  int index = -1;
  double r2max = 1.0e60;
  out.Distance = sqrt(r2max);

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, thread_id, numnodes, firstnode);
  
  if(nfound < 0)
    return -1;

  for(int n = 0; n < nfound; n++)
    {
      int j = Thread[thread_id].Ngblist[n];

      double dx = NGB_PERIODIC_LONG_X(pos[0] - P[j].Pos[0]);
      double dy = NGB_PERIODIC_LONG_Y(pos[1] - P[j].Pos[1]);
      double dz = NGB_PERIODIC_LONG_Z(pos[2] - P[j].Pos[2]);

      double r2 = dx * dx + dy * dy + dz * dz;
      if(r2 < r2max && r2 < h * h)
        {
          index = j;
          r2max = r2;
        }
    }

  if(index >= 0)
    {
      if(index >= NumGas)
        terminate("index >= NumGas");

      out.Distance = sqrt(r2max);

#ifdef MHD
      for(int m = 0; m < 3; m++)
        out.B[m] = SphP[index].B[m];
#endif
      for(int m = 0; m < 3; m++)
        out.Vel[m] = P[index].Vel[m] - CMVel[m];
      out.Density = SphP[index].Density;
    }
  else
    {
      out.Distance = 1.e60;
    }

  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#endif
