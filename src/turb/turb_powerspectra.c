/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/turb/turb_powerspectra.c
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
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "../allvars.h"
#include "../proto.h"


#if defined(POWERSPEC_GRID) && defined(PERIODIC) && (defined(VS_TURB) || defined(AB_TURB))

#define  POWERSPEC_GRIDz (POWERSPEC_GRID/2 + 1)
#define  POWERSPEC_GRID2 (2*POWERSPEC_GRIDz)


void powerspec_turb_calc_and_bin_spectrum(fft_real * field, int flag);
void powerspec_dump_field(fft_real * density, int num);


#if (POWERSPEC_GRID > 1024)
typedef long long large_array_offset;
#else
typedef unsigned int large_array_offset;
#endif


static size_t maxfftsize;

static fft_plan myplan;

#ifdef MHD
static fft_real *bfield[3];
#endif
static fft_real *velfield[3];
static fft_real *vorticityfield[3];
static fft_real *velrhofield[3];
static fft_real *dis1field;
static fft_real *dis2field;
static fft_real *densityfield;

static fft_real *randomfield;
static fft_real *workspace;

static float *RandomValue;

static fft_complex *fft_of_field;

static double *powerspec_turb_nearest_distance, *powerspec_turb_nearest_hsml;


#define BINS_PS  2000           /* number of bins for power spectrum computation */

static long long CountModes[BINS_PS];
static double SumPower[BINS_PS];
static double Power[BINS_PS];
static double Kbin[BINS_PS];
static double K0, K1;
static double binfac;
#ifdef MHD
static double b_disp[3];
#endif
static double vel_disp[3];
static double velrho_disp[3];
static double empty_disp[3] = { 0, 0, 0 };

#ifdef DG
#define MAX_LEVEL_DIFF 4

double (*dg_base_function_level)[NOF_BASE_FUNCTIONS];
int GridLevel;
int Offsets[MAX_LEVEL_DIFF];
#endif

static int particletype;

typedef struct data_in
{
  MyDouble Pos[3];
  MyFloat Hsml;
  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

static void particle2in(struct data_in *in, int place, int firstnode)
{
  int xx = place / (POWERSPEC_GRID * POWERSPEC_GRID);
  int yy = (place - xx * POWERSPEC_GRID * POWERSPEC_GRID) / POWERSPEC_GRID;
  int zz = (place - xx * POWERSPEC_GRID * POWERSPEC_GRID - yy * POWERSPEC_GRID);
  xx += myplan.slabstart_x;

  double x = (xx + 0.5) / POWERSPEC_GRID * All.BoxSize;
  double y = (yy + 0.5) / POWERSPEC_GRID * All.BoxSize;
  double z = (zz + 0.5) / POWERSPEC_GRID * All.BoxSize;

  in->Pos[0] = x;
  in->Pos[1] = y;
  in->Pos[2] = z;
  in->Hsml = powerspec_turb_nearest_hsml[place];

  in->Firstnode = firstnode;
}

typedef struct data_out
{
  double Distance;
#ifdef MHD
  MyDouble B[3];
#endif
  MyDouble Vel[3];
  MyDouble Vorticity[3];
  MyDouble Density;
  MyDouble DuDt_diss;
  MyDouble DuDt_drive;
  MyDouble RandomValue;
#ifdef TRACER_MC
  int NtracerMC;
#endif
} data_out;

static data_out *DataResult, *DataOut;

static void out2particle(struct data_out *out, int place, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES)      /* initial store */
    {
      powerspec_turb_nearest_distance[place] = sqrt(out->Distance);

      int i = place / (POWERSPEC_GRID * POWERSPEC_GRID);
      int j = (place - i * POWERSPEC_GRID * POWERSPEC_GRID) / POWERSPEC_GRID;
      int k = (place - i * POWERSPEC_GRID * POWERSPEC_GRID - j * POWERSPEC_GRID);
      int ip = POWERSPEC_GRID2 * (POWERSPEC_GRID * i + j) + k;

      int m;

#ifdef MHD
      for(m = 0; m < 3; m++)
        {
          bfield[m][ip] = out->B[m];
        }
#endif	

      for(m = 0; m < 3; m++)
        {
          velfield[m][ip] = out->Vel[m];
        }

      for(m = 0; m < 3; m++)
        {
          velrhofield[m][ip] = sqrt(out->Density) * out->Vel[m];
        }
#ifdef TRACER_MC
      if(particletype == TRACER_MC)
        for(m = 0; m < 3; m++)
          velrhofield[m][ip] *= sqrt(out->NtracerMC / (double) All.TracerMCPerCell);
#endif

      for(m = 0; m < 3; m++)
        {
          vorticityfield[m][ip] = out->Vorticity[m];
        }


      if(out->DuDt_diss >= 0)
        {
          dis1field[ip] = sqrt(out->DuDt_diss);
          dis2field[ip] = 0;
        }
      else
        {
          dis1field[ip] = 0;
          dis2field[ip] = sqrt(-out->DuDt_diss);
        }

      randomfield[ip] = out->RandomValue;
      densityfield[ip] = out->Density;

#ifdef TRACER_MC
      if(particletype == TRACER_MC)
        densityfield[ip] *= sqrt(out->NtracerMC / (double) All.TracerMCPerCell);
#endif
    }
  else                          /* combine */
    {
      if(out->Distance < powerspec_turb_nearest_distance[place])
        {
          powerspec_turb_nearest_distance[place] = out->Distance;

          int ii = place / (POWERSPEC_GRID * POWERSPEC_GRID);
          int jj = (place - ii * POWERSPEC_GRID * POWERSPEC_GRID) / POWERSPEC_GRID;
          int kk = (place - ii * POWERSPEC_GRID * POWERSPEC_GRID - jj * POWERSPEC_GRID);
          int ip = POWERSPEC_GRID2 * (POWERSPEC_GRID * ii + jj) + kk;

#ifdef MHD
          bfield[0][ip] = out->B[0];
          bfield[1][ip] = out->B[1];
          bfield[2][ip] = out->B[2];
#endif

          velfield[0][ip] = out->Vel[0];
          velfield[1][ip] = out->Vel[1];
          velfield[2][ip] = out->Vel[2];

          velrhofield[0][ip] = sqrt(out->Density) * out->Vel[0];
          velrhofield[1][ip] = sqrt(out->Density) * out->Vel[1];
          velrhofield[2][ip] = sqrt(out->Density) * out->Vel[2];

          vorticityfield[0][ip] = out->Vorticity[0];
          vorticityfield[1][ip] = out->Vorticity[1];
          vorticityfield[2][ip] = out->Vorticity[2];



          if(out->DuDt_diss >= 0)
            {
              dis1field[ip] = sqrt(out->DuDt_diss);
              dis2field[ip] = 0;
            }
          else
            {
              dis1field[ip] = 0;
              dis2field[ip] = sqrt(-out->DuDt_diss);
            }

          randomfield[ip] = out->RandomValue;
          densityfield[ip] = out->Density;
        }
    }
}

#include "../generic_comm_helpers2.h"

static large_array_offset Ncount;

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

        if(powerspec_turb_nearest_distance[i] > 1.0e29)
          powerspec_turb_find_nearest_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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

        powerspec_turb_find_nearest_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}


void powersepc_turb_init()
{
  /* temporarily allocate some arrays to make sure that out-of-place plans are created */
  densityfield = (fft_real *) mymalloc("densityfield", POWERSPEC_GRID2 * sizeof(fft_real));
  workspace = (fft_real *) mymalloc("workspace", POWERSPEC_GRID2 * sizeof(fft_real));

#ifdef DOUBLEPRECISION_FFTW
  int alignflag = 0;
#else
  /* for single precision, the start of our FFT columns is presently only guaranteed to be 8-byte aligned */
  int alignflag = FFTW_UNALIGNED;
#endif

  int stride = POWERSPEC_GRIDz;

  int ndim[1] = { POWERSPEC_GRID };     /* dimension of the 1D transforms */
  /* Set up the FFTW plan  */
  myplan.forward_plan_zdir =
    FFTW(plan_many_dft_r2c) (1, ndim, 1, densityfield, 0, 1, POWERSPEC_GRID2, (fft_complex *) workspace, 0, 1, POWERSPEC_GRIDz, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  myplan.forward_plan_ydir =
    FFTW(plan_many_dft) (1, ndim, 1, (fft_complex *) densityfield, 0, stride, POWERSPEC_GRIDz * POWERSPEC_GRID, (fft_complex *) workspace, 0, stride, POWERSPEC_GRIDz * POWERSPEC_GRID, FFTW_FORWARD,
                         FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  myplan.forward_plan_xdir =
    FFTW(plan_many_dft) (1, ndim, 1, (fft_complex *) densityfield, 0, stride, POWERSPEC_GRIDz * POWERSPEC_GRID, (fft_complex *) workspace, 0, stride, POWERSPEC_GRIDz * POWERSPEC_GRID, FFTW_FORWARD,
                         FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  myfree(workspace);
  myfree(densityfield);


  my_slab_based_fft_init(&myplan, POWERSPEC_GRID, POWERSPEC_GRID, POWERSPEC_GRID);

  maxfftsize = imax(myplan.largest_x_slab * POWERSPEC_GRID, myplan.largest_y_slab * POWERSPEC_GRID) * ((size_t) POWERSPEC_GRID2);


#ifdef DG
  int l, i, x, y, z;

  int points = 0;
  for(l = 0; l < MAX_LEVEL_DIFF; l++)
    {
      points += 1 << (3 * l);
    }

  GridLevel = 0;
  int tmp = POWERSPEC_GRID;
  while(tmp > 1)
    {
      GridLevel++;
      tmp /= 2;
    }

  dg_base_function_level = (double (*)[NOF_BASE_FUNCTIONS]) mymalloc("dg_weights_level", points * NOF_BASE_FUNCTIONS * sizeof(double));

  int offset = 0;
  for(l = 0; l < MAX_LEVEL_DIFF; l++)
    {
      Offsets[l] = offset;
      double xx = -1. + 1. / (1 << (l));
      for(x = 0; x < (1 << l); x++)
        {
          double yy = -1. + 1. / (1 << (l));
          for(y = 0; y < (1 << l); y++)
            {
              double zz = -1. + 1. / (1 << (l));
              for(z = 0; z < (1 << l); z++)
                {
                  for(i = 0; i < Nof_base_functions; i++)
                    {
                      dg_base_function_level[offset][i] = base_function(i, xx, yy, zz);
                    }

                  offset++;
                  zz += 1. / (1 << (l - 1));
                }
              yy += 1. / (1 << (l - 1));
            }
          xx += 1. / (1 << (l - 1));
        }
    }
#endif
}

void powerspec_turb(int filenr, int type)
{
  int i;
  char fname[1000];

  particletype = type;

  mpi_printf("TURB: start turbulent powerspec computation for type=%d (filenr=%d)\n", type, filenr);

  double tstart = second();

  /* allocate the memory to hold the FFT fields */
#ifdef MHD
  bfield[0] = (fft_real *) mymalloc("bfield[0]", maxfftsize * sizeof(fft_real));
  bfield[1] = (fft_real *) mymalloc("bfield[1]", maxfftsize * sizeof(fft_real));
  bfield[2] = (fft_real *) mymalloc("bfield[2]", maxfftsize * sizeof(fft_real));
#endif

  velfield[0] = (fft_real *) mymalloc("velfield[0]", maxfftsize * sizeof(fft_real));
  velfield[1] = (fft_real *) mymalloc("velfield[1]", maxfftsize * sizeof(fft_real));
  velfield[2] = (fft_real *) mymalloc("velfield[2]", maxfftsize * sizeof(fft_real));

  velrhofield[0] = (fft_real *) mymalloc("velrhofield[0]", maxfftsize * sizeof(fft_real));
  velrhofield[1] = (fft_real *) mymalloc("velrhofield[1]", maxfftsize * sizeof(fft_real));
  velrhofield[2] = (fft_real *) mymalloc("velrhofield[2]", maxfftsize * sizeof(fft_real));

  vorticityfield[0] = (fft_real *) mymalloc("vorticityfield[0]", maxfftsize * sizeof(fft_real));
  vorticityfield[1] = (fft_real *) mymalloc("vorticityfield[1]", maxfftsize * sizeof(fft_real));
  vorticityfield[2] = (fft_real *) mymalloc("vorticityfield[2]", maxfftsize * sizeof(fft_real));

  dis1field = (fft_real *) mymalloc("dis1field", maxfftsize * sizeof(fft_real));
  dis2field = (fft_real *) mymalloc("dis2field", maxfftsize * sizeof(fft_real));
  randomfield = (fft_real *) mymalloc("randomfield", maxfftsize * sizeof(fft_real));

  densityfield = (fft_real *) mymalloc("densityfield", maxfftsize * sizeof(fft_real));

  workspace = (fft_real *) mymalloc("workspace", maxfftsize * sizeof(fft_real));

#ifdef MHD
  memset(bfield[0], 0, maxfftsize * sizeof(fft_real));
  memset(bfield[1], 0, maxfftsize * sizeof(fft_real));
  memset(bfield[2], 0, maxfftsize * sizeof(fft_real));
#endif

  memset(velfield[0], 0, maxfftsize * sizeof(fft_real));
  memset(velfield[1], 0, maxfftsize * sizeof(fft_real));
  memset(velfield[2], 0, maxfftsize * sizeof(fft_real));

  memset(velrhofield[0], 0, maxfftsize * sizeof(fft_real));
  memset(velrhofield[1], 0, maxfftsize * sizeof(fft_real));
  memset(velrhofield[2], 0, maxfftsize * sizeof(fft_real));

  memset(vorticityfield[0], 0, maxfftsize * sizeof(fft_real));
  memset(vorticityfield[1], 0, maxfftsize * sizeof(fft_real));
  memset(vorticityfield[2], 0, maxfftsize * sizeof(fft_real));

  memset(dis1field, 0, maxfftsize * sizeof(fft_real));
  memset(dis2field, 0, maxfftsize * sizeof(fft_real));
  memset(randomfield, 0, maxfftsize * sizeof(fft_real));

  memset(densityfield, 0, maxfftsize * sizeof(fft_real));

  RandomValue = (float *) mymalloc("RndField", NumGas * sizeof(float));

  gsl_rng *random_gen = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(random_gen, 42 + ThisTask);       /* start-up seed */

  for(i = 0; i < NumGas; i++)
    RandomValue[i] = gsl_ran_gaussian(random_gen, 1.0);

  powerspec_turb_obtain_fields(type);

  powerspec_turb_calc_dispersion();

  //powerspec_dump_field(densityfield, filenr);


  /* Now compute the power spectrum of the B field */
  fname[0] = 0;
  for(i = 0; i < BINS_PS; i++)
    {
      SumPower[i] = 0;
      CountModes[i] = 0;
    }

#ifdef MHD
  powerspec_turb_calc_and_bin_spectrum(bfield[0], 1); /* only here the modes are counted */
  powerspec_turb_calc_and_bin_spectrum(bfield[1], 0);
  powerspec_turb_calc_and_bin_spectrum(bfield[2], 0);

  powerspec_turb_collect();

  if(type == 0)
    sprintf(fname, "%s/powerspec_B_%03d.txt", All.OutputDir, filenr);
  if(fname[0] != 0)
    powerspec_turb_save(fname, b_disp);
#endif

  /* Now compute the power spectrum of the velocities */
  fname[0] = 0;
  for(i = 0; i < BINS_PS; i++)
    {
      SumPower[i] = 0;
      CountModes[i] = 0;
    }

  powerspec_turb_calc_and_bin_spectrum(velfield[0], 1); /* only here the modes are counted */
  powerspec_turb_calc_and_bin_spectrum(velfield[1], 0);
  powerspec_turb_calc_and_bin_spectrum(velfield[2], 0);

  powerspec_turb_collect();

  if(type == 0)
    sprintf(fname, "%s/powerspec_vel_%03d.txt", All.OutputDir, filenr);
#ifdef TRACER_PARTICLE
  if(type == TRACER_PARTICLE)
    sprintf(fname, "%s/powerspec_TRPART_vel_%03d.txt", All.OutputDir, filenr);
#endif
  if(fname[0] != 0)
    powerspec_turb_save(fname, vel_disp);


  /* now compute the power spectrum of the sqrt(rho)-weighted veloicty */
  fname[0] = 0;
  for(i = 0; i < BINS_PS; i++)
    {
      SumPower[i] = 0;
      CountModes[i] = 0;
    }

  powerspec_turb_calc_and_bin_spectrum(velrhofield[0], 1);
  powerspec_turb_calc_and_bin_spectrum(velrhofield[1], 0);
  powerspec_turb_calc_and_bin_spectrum(velrhofield[2], 0);

  powerspec_turb_collect();

  if(type == 0)
    sprintf(fname, "%s/powerspec_velrho_%03d.txt", All.OutputDir, filenr);
#ifdef TRACER_PARTICLE
  if(type == TRACER_PARTICLE)
    sprintf(fname, "%s/powerspec_TRPART_velrho_%03d.txt", All.OutputDir, filenr);
#endif
#ifdef TRACER_MC
  if(type == TRACER_MC)
    sprintf(fname, "%s/powerspec_TRMC_velrho_%03d.txt", All.OutputDir, filenr);
#endif
  if(fname[0] != 0)
    powerspec_turb_save(fname, velrho_disp);


  /* now compute the power spectrum of the vorticity */
  fname[0] = 0;
  for(i = 0; i < BINS_PS; i++)
    {
      SumPower[i] = 0;
      CountModes[i] = 0;
    }

  powerspec_turb_calc_and_bin_spectrum(vorticityfield[0], 1);
  powerspec_turb_calc_and_bin_spectrum(vorticityfield[1], 0);
  powerspec_turb_calc_and_bin_spectrum(vorticityfield[2], 0);

  powerspec_turb_collect();

  if(type == 0)
    sprintf(fname, "%s/powerspec_vorticity_%03d.txt", All.OutputDir, filenr);
#ifdef TRACER_PARTICLE
  if(type == TRACER_PARTICLE)
    sprintf(fname, "%s/powerspec_TRPART_vorticity_%03d.txt", All.OutputDir, filenr);
#endif
  if(fname[0] != 0)
    powerspec_turb_save(fname, velrho_disp);


  /* Now compute the power spectrum of the dissipation1 */
  fname[0] = 0;
  for(i = 0; i < BINS_PS; i++)
    {
      SumPower[i] = 0;
      CountModes[i] = 0;
    }

  powerspec_turb_calc_and_bin_spectrum(dis1field, 1);

  powerspec_turb_collect();

  if(type == 0)
    sprintf(fname, "%s/powerspec_dis1_%03d.txt", All.OutputDir, filenr);
  if(fname[0] != 0)
    powerspec_turb_save(fname, empty_disp);


  /* Now compute the power spectrum of the dissipation2 */
  fname[0] = 0;
  for(i = 0; i < BINS_PS; i++)
    {
      SumPower[i] = 0;
      CountModes[i] = 0;
    }

  powerspec_turb_calc_and_bin_spectrum(dis2field, 1);

  powerspec_turb_collect();

  if(type == 0)
    sprintf(fname, "%s/powerspec_dis2_%03d.txt", All.OutputDir, filenr);
  if(fname[0] != 0)
    powerspec_turb_save(fname, empty_disp);


  /* Now compute the power spectrum of the random field */
  fname[0] = 0;
  for(i = 0; i < BINS_PS; i++)
    {
      SumPower[i] = 0;
      CountModes[i] = 0;
    }

  powerspec_turb_calc_and_bin_spectrum(randomfield, 1);

  powerspec_turb_collect();

  if(type == 0)
    sprintf(fname, "%s/powerspec_random_%03d.txt", All.OutputDir, filenr);
  if(fname[0] != 0)
    powerspec_turb_save(fname, empty_disp);


  /* Now compute the power spectrum of the density field */
  fname[0] = 0;
  for(i = 0; i < BINS_PS; i++)
    {
      SumPower[i] = 0;
      CountModes[i] = 0;
    }

  powerspec_turb_calc_and_bin_spectrum(densityfield, 1);

  powerspec_turb_collect();

  if(type == 0)
    sprintf(fname, "%s/powerspec_density_%03d.txt", All.OutputDir, filenr);
#ifdef TRACER_PARTICLE
  if(type == TRACER_PARTICLE)
    sprintf(fname, "%s/powerspec_TRPART_density_%03d.txt", All.OutputDir, filenr);
#endif
#ifdef TRACER_MC
  if(type == TRACER_MC)
    sprintf(fname, "%s/powerspec_TRMC_density_%03d.txt", All.OutputDir, filenr);
#endif
  if(fname[0] != 0)
    powerspec_turb_save(fname, empty_disp);


  /* finish */
  myfree(RandomValue);

  myfree(workspace);
  myfree(densityfield);
  myfree(randomfield);
  myfree(dis2field);
  myfree(dis1field);
  myfree(vorticityfield[2]);
  myfree(vorticityfield[1]);
  myfree(vorticityfield[0]);
  myfree(velrhofield[2]);
  myfree(velrhofield[1]);
  myfree(velrhofield[0]);
  myfree(velfield[2]);
  myfree(velfield[1]);
  myfree(velfield[0]);
#ifdef MHD
  myfree(bfield[2]);
  myfree(bfield[1]);
  myfree(bfield[0]);
#endif
  double tend = second();

  mpi_printf("TURB: end turbulent power spectra  took %g seconds\n", timediff(tstart, tend));
}



void powerspec_turb_calc_and_bin_spectrum(fft_real * field, int flag)
{
  double k2, kx, ky, kz;
  int x, y, z, zz;

  K0 = 2 * M_PI / All.BoxSize;  /* minimum k */
  K1 = K0 * POWERSPEC_GRID / 2; /* maximum k */
  binfac = BINS_PS / (log(K1) - log(K0));

  /* Do the FFT of field */

  my_slab_based_fft(&myplan, &field[0], &workspace[0], 1);

  fft_of_field = (fft_complex *) & field[0];

#pragma omp parallel for private(y, z)
  for(x = 0; x < POWERSPEC_GRID; x++)
    for(y = myplan.slabstart_y; y < myplan.slabstart_y + myplan.nslab_y; y++)
      for(z = 0; z < POWERSPEC_GRID; z++)
        {
          zz = z;
          if(z >= POWERSPEC_GRID / 2 + 1)
            zz = POWERSPEC_GRID - z;

          if(x > POWERSPEC_GRID / 2)
            kx = x - POWERSPEC_GRID;
          else
            kx = x;
          if(y > POWERSPEC_GRID / 2)
            ky = y - POWERSPEC_GRID;
          else
            ky = y;
          if(z > POWERSPEC_GRID / 2)
            kz = z - POWERSPEC_GRID;
          else
            kz = z;

          k2 = kx * kx + ky * ky + kz * kz;

          large_array_offset ip = ((large_array_offset) POWERSPEC_GRIDz) * (POWERSPEC_GRID * (y - myplan.slabstart_y) + x) + zz;

          double po = (fft_of_field[ip][0] * fft_of_field[ip][0] + fft_of_field[ip][1] * fft_of_field[ip][1]) / pow(POWERSPEC_GRID, 6);

          if(k2 > 0)
            {
              if(k2 < (POWERSPEC_GRID / 2.0) * (POWERSPEC_GRID / 2.0))
                {
                  double k = sqrt(k2) * 2 * M_PI / All.BoxSize;

                  if(k >= K0 && k < K1)
                    {
                      int bin = log(k / K0) * binfac;

                      SumPower[bin] += po;

                      if(flag)
                        CountModes[bin] += 1;
                    }
                }
            }
        }
}



void powerspec_turb_collect(void)
{
  int i, n;
  long long int *countbuf = (long long int *) mymalloc("countbuf", NTask * BINS_PS * sizeof(long long));
  double *powerbuf = (double *) mymalloc("powerbuf", NTask * BINS_PS * sizeof(double));

  MPI_Allgather(CountModes, BINS_PS * sizeof(long long), MPI_BYTE, countbuf, BINS_PS * sizeof(long long), MPI_BYTE, MPI_COMM_WORLD);

  for(i = 0; i < BINS_PS; i++)
    {
      CountModes[i] = 0;
      for(n = 0; n < NTask; n++)
        CountModes[i] += countbuf[n * BINS_PS + i];
    }

  MPI_Allgather(SumPower, BINS_PS * sizeof(double), MPI_BYTE, powerbuf, BINS_PS * sizeof(double), MPI_BYTE, MPI_COMM_WORLD);

  for(i = 0; i < BINS_PS; i++)
    {
      SumPower[i] = 0;
      for(n = 0; n < NTask; n++)
        SumPower[i] += powerbuf[n * BINS_PS + i];
    }

  myfree(powerbuf);
  myfree(countbuf);

  for(i = 0; i < BINS_PS; i++)
    {
      Kbin[i] = exp((i + 0.5) / binfac + log(K0));

      if(CountModes[i] > 0)
        Power[i] = SumPower[i] / CountModes[i];
      else
        Power[i] = 0;
    }
}



void powerspec_turb_save(char *fname, double *disp)
{
  FILE *fd;
  char buf[500];
  int i;

  if(ThisTask == 0)
    {
      if(!(fd = fopen(fname, "w")))
        {
          sprintf(buf, "can't open file `%s`\n", fname);
          terminate(buf);
        }

      fprintf(fd, "%g\n", All.Time);
      i = POWERSPEC_GRID;
      fprintf(fd, "%d\n", i);
      i = BINS_PS;
      fprintf(fd, "%d\n", i);

      fprintf(fd, "%g\n", disp[0]);
      fprintf(fd, "%g\n", disp[1]);
      fprintf(fd, "%g\n", disp[2]);

      for(i = 0; i < BINS_PS; i++)
        {
          fprintf(fd, "%g %g %g %g\n", Kbin[i], Power[i], (double) CountModes[i], SumPower[i]);
        }

      fclose(fd);
    }
}


/* this function determines the velocity fields by using the nearest cell's values 
 */
double powerspec_turb_obtain_fields(int type)
{
  double tstart = second();
  long long ntot, npleft, i;
  int iter;

  mpi_printf("TURB: Start finding nearest gas-particle for mesh-cell centers (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

  large_array_offset n;
  Ncount = ((large_array_offset) myplan.nslab_x) * (POWERSPEC_GRID * POWERSPEC_GRID);   /* number of grid points on the local slab */


  powerspec_turb_nearest_distance = (double *) mymalloc("powerspec_turb_nearest_distance", sizeof(double) * Ncount);
  powerspec_turb_nearest_hsml = (double *) mymalloc("powerspec_turb_nearest_hsml", sizeof(double) * Ncount);

  for(n = 0; n < Ncount; n++)
    {
      powerspec_turb_nearest_distance[n] = 1.0e30;
      powerspec_turb_nearest_hsml[n] = All.BoxSize / pow(All.TotNumGas, 1.0 / 3);
    }

  generic_set_MaxNexport();

  iter = 0;
  do {
  generic_comm_pattern(Ncount, kernel_local, kernel_imported);

  /* do final operations on results */
  for(i = 0, npleft = 0; i < Ncount; i++)
    {
      if(powerspec_turb_nearest_distance[i] > 1.0e29)
        {
          /* need to redo this particle */
          npleft++;
          powerspec_turb_nearest_hsml[i] *= 2.0;
        }
       else
        {
          powerspec_turb_nearest_distance[i] = 0;    /* we not continue to search for this particle */
        }
     }
    sumup_longs(1, &npleft, &ntot);
    if(ntot > 0)
      {
        iter++;
        if(iter > 0 && ThisTask == 0)
          {
            printf("TURB: powerpec_vel nearest iteration %d: need to repeat for %lld particles.\n", iter, ntot);
            myflush(stdout);
          }

        if(iter > MAXITER)
          terminate("failed to converge");
      }
    }
  while(ntot > 0);

  myfree(powerspec_turb_nearest_hsml);
  myfree(powerspec_turb_nearest_distance);

  mpi_printf("done finding velocity field\n");

  double tend = second();
  return timediff(tstart, tend);
}

static inline unsigned long long ngb_double_to_int(double d)
{
  union
  {
    double d;
    unsigned long long ull;
  } u;
  u.d = d;
  return (u.ull & 0xFFFFFFFFFFFFFllu);
}

int powerspec_turb_find_nearest_evaluate(int target, int mode, int thread_id)
{
  int j, n;
  int m;
  int numnodes, *firstnode;
#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif

  double h, r2max;
  double dx, dy, dz, r2;
  MyDouble *pos;

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

  pos = target_data->Pos;
  h = target_data->Hsml;

  int index = -1;
  r2max = 1.0e60;
  out.Distance = sqrt(r2max);

#ifndef AMR
  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, thread_id, numnodes, firstnode);
#else
  int nfound = amr_treefind_single_threads(pos, target, mode, thread_id, numnodes, firstnode);
#endif

  if(nfound < 0)
    return -1;

  for(n = 0; n < nfound; n++)
    {
      j = Thread[thread_id].Ngblist[n];

      dx = NGB_PERIODIC_LONG_X(pos[0] - P[j].Pos[0]);
      dy = NGB_PERIODIC_LONG_Y(pos[1] - P[j].Pos[1]);
      dz = NGB_PERIODIC_LONG_Z(pos[2] - P[j].Pos[2]);

      r2 = dx * dx + dy * dy + dz * dz;
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

#ifndef DG
#ifdef MHD
      for(m = 0; m < 3; m++)
        out.B[m] = SphP[index].B[m];
#endif
      for(m = 0; m < 3; m++)
        out.Vel[m] = P[index].Vel[m];
      for(m = 0; m < 3; m++)
        out.Vorticity[m] = SphP[index].Vorticity[m];
      out.Density = SphP[index].Density;
#ifdef TRACER_MC
      out.NtracerMC = get_number_of_tracers(index);
#endif
      out.DuDt_diss = SphP[index].DuDt_diss;
      out.DuDt_drive = SphP[index].DuDt_drive;
      out.RandomValue = RandomValue[index];
#else
      int level_diff = GridLevel - Mesh.DP[j].level;

      assert(level_diff >= 0 && level_diff < MAX_LEVEL_DIFF);

      unsigned long long xxb = ngb_double_to_int(((pos[0] - P[j].Pos[0] + amr_length[Mesh.DP[j].level + 1]) / amr_length[Mesh.DP[j].level]) + 1.0);
      unsigned long long yyb = ngb_double_to_int(((pos[1] - P[j].Pos[1] + amr_length[Mesh.DP[j].level + 1]) / amr_length[Mesh.DP[j].level]) + 1.0);
      unsigned long long zzb = ngb_double_to_int(((pos[2] - P[j].Pos[2] + amr_length[Mesh.DP[j].level + 1]) / amr_length[Mesh.DP[j].level]) + 1.0);

      unsigned long long mask = 0;
      int i;
      for(i = 0; i < level_diff; i++)
        {
          mask |= ((unsigned long long) 1) << (51 - i);
        }


      int shift_x = (52 - 3 * level_diff);
      int shift_y = (52 - 2 * level_diff);
      int shift_z = (52 - 1 * level_diff);

      int point = Offsets[level_diff] + ((xxb & mask) >> shift_x) + ((yyb & mask) >> shift_y) + ((zzb & mask) >> shift_z);

      out.Density = 0.;
      for(m = 0; m < 3; m++)
        out.Vel[m] = 0.;

      out.DuDt_diss = 0.;


      out.DuDt_drive = SphP[index].DuDt_drive;
      out.RandomValue = RandomValue[index];

      //loop over base functions
      int l;
      for(l = 0; l < Nof_base_functions; l++)
        {
          out.Density += SphP[index].Weights[l][0] * dg_base_function_level[point][l];
          out.DuDt_diss += SphP[index].WeightsDiss[l] * dg_base_function_level[point][l];
          for(m = 0; m < 3; m++)
            out.Vel[m] += SphP[index].Weights[l][m + 1] * dg_base_function_level[point][l];
        }
      for(m = 0; m < 3; m++)
        out.Vel[m] /= out.Density;

      out.DuDt_diss /= out.Density;
#endif
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



void powerspec_turb_calc_dispersion(void)
{
  int dim, i, j, k;

#ifdef MHD
  for(dim = 0; dim < 3; dim++)
    {
      double vsum = 0, vsum_all, vmean, vdisp = 0, vdisp_all;

      for(i = 0; i < myplan.nslab_x; i++)
        for(j = 0; j < POWERSPEC_GRID; j++)
          for(k = 0; k < POWERSPEC_GRID; k++)
            {
              int ip = POWERSPEC_GRID2 * (POWERSPEC_GRID * i + j) + k;

              vsum += bfield[dim][ip];
            }

      MPI_Allreduce(&vsum, &vsum_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      vmean = vsum_all / pow(POWERSPEC_GRID, 3);

      for(i = 0; i < myplan.nslab_x; i++)
        for(j = 0; j < POWERSPEC_GRID; j++)
          for(k = 0; k < POWERSPEC_GRID; k++)
            {
              int ip = POWERSPEC_GRID2 * (POWERSPEC_GRID * i + j) + k;

              bfield[dim][ip] -= vmean;
            }

      for(i = 0; i < myplan.nslab_x; i++)
        for(j = 0; j < POWERSPEC_GRID; j++)
          for(k = 0; k < POWERSPEC_GRID; k++)
            {
              int ip = POWERSPEC_GRID2 * (POWERSPEC_GRID * i + j) + k;

              vdisp += bfield[dim][ip] * bfield[dim][ip];
            }

      MPI_Allreduce(&vdisp, &vdisp_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      b_disp[dim] = vdisp_all / pow(POWERSPEC_GRID, 3);
    }
#endif

  for(dim = 0; dim < 3; dim++)
    {
      double vsum = 0, vsum_all, vmean, vdisp = 0, vdisp_all;

      for(i = 0; i < myplan.nslab_x; i++)
        for(j = 0; j < POWERSPEC_GRID; j++)
          for(k = 0; k < POWERSPEC_GRID; k++)
            {
              int ip = POWERSPEC_GRID2 * (POWERSPEC_GRID * i + j) + k;

              vsum += velfield[dim][ip];
            }

      MPI_Allreduce(&vsum, &vsum_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      vmean = vsum_all / pow(POWERSPEC_GRID, 3);

      for(i = 0; i < myplan.nslab_x; i++)
        for(j = 0; j < POWERSPEC_GRID; j++)
          for(k = 0; k < POWERSPEC_GRID; k++)
            {
              int ip = POWERSPEC_GRID2 * (POWERSPEC_GRID * i + j) + k;

              velfield[dim][ip] -= vmean;
            }

      for(i = 0; i < myplan.nslab_x; i++)
        for(j = 0; j < POWERSPEC_GRID; j++)
          for(k = 0; k < POWERSPEC_GRID; k++)
            {
              int ip = POWERSPEC_GRID2 * (POWERSPEC_GRID * i + j) + k;

              vdisp += velfield[dim][ip] * velfield[dim][ip];
            }

      MPI_Allreduce(&vdisp, &vdisp_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      vel_disp[dim] = vdisp_all / pow(POWERSPEC_GRID, 3);
    }


  for(dim = 0; dim < 3; dim++)
    {
      double vsum = 0, vsum_all, vmean, vdisp = 0, vdisp_all;

      for(i = 0; i < myplan.nslab_x; i++)
        for(j = 0; j < POWERSPEC_GRID; j++)
          for(k = 0; k < POWERSPEC_GRID; k++)
            {
              int ip = POWERSPEC_GRID2 * (POWERSPEC_GRID * i + j) + k;

              vsum += velrhofield[dim][ip];
            }

      MPI_Allreduce(&vsum, &vsum_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      vmean = vsum_all / pow(POWERSPEC_GRID, 3);

      for(i = 0; i < myplan.nslab_x; i++)
        for(j = 0; j < POWERSPEC_GRID; j++)
          for(k = 0; k < POWERSPEC_GRID; k++)
            {
              int ip = POWERSPEC_GRID2 * (POWERSPEC_GRID * i + j) + k;

              velrhofield[dim][ip] -= vmean;
            }

      for(i = 0; i < myplan.nslab_x; i++)
        for(j = 0; j < POWERSPEC_GRID; j++)
          for(k = 0; k < POWERSPEC_GRID; k++)
            {
              int ip = POWERSPEC_GRID2 * (POWERSPEC_GRID * i + j) + k;

              vdisp += velrhofield[dim][ip] * velrhofield[dim][ip];
            }

      MPI_Allreduce(&vdisp, &vdisp_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      velrho_disp[dim] = vdisp_all / pow(POWERSPEC_GRID, 3);

    }
}

#endif
