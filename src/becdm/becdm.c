/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/becdm/becdm.c
 * \date        02/2017
 * \author      Philip Mocz
 * \brief       includes functions for BECDM pseudo-spectral algorithm
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>

#include "../allvars.h"
#include "../proto.h"


#if defined(PMGRID) && defined(PERIODIC) && defined(DOUBLEPRECISION_FFTW) && defined(BECDM)

#if (PMGRID > 1024)
typedef long long large_array_offset;
#else
typedef unsigned int large_array_offset;
#endif

static fft_complex *densityfield, *psifield, *workspace;

static fft_plan myplan;

static unsigned int *psi_send_count;       //!< psi_send_count[i]: number of elements which have to be send to task i
static unsigned int *psi_recv_count;       //!< psi_recv_count[i]: number of elements which have to be received from task i

typedef struct
{
  double psiRe;
  double psiIm;
  large_array_offset ip;
} psiType;


/** Potential kick operator
 *  psi = exp(-1.i * dt/2 * m/hbar * V).*psi;
 * psi is psi_comoving = a^(3/2) psi_physical
 *  dt_gravkick is int_{t}^{t+dt} dt/a
 **/
void do_becdm_potential_kick(int i, double dt_gravkick) {
  double psiRe = P[i].PsiRe;
  double psiIm = P[i].PsiIm;
  double cosx = cos(dt_gravkick/2. * All.mAxion / All.hbar * P[i].Potential);
  double sinx = sin(dt_gravkick/2. * All.mAxion / All.hbar * P[i].Potential);
  P[i].PsiRe = cosx * psiRe +  sinx * psiIm;  // evolves as  psi = exp(-i dt V) * psi
  P[i].PsiIm = cosx * psiIm -  sinx * psiRe;
}


/** Kinetic drift operator
 *  % kinetic - drift
 *  psi = fftn(psi);
 *  psi = exp(dt * (hbar/m) * (-1.i*((2*pi/Lbox)^2*double(kSq))/2)) .*psi;
 *  psi = ifftn(psi);
 **/
void do_becdm_kinetic_drift(void) {
  assert(PMGRID % NTask == 0); // assume # processors divides PMGRID evenly, for fft optimality

  //memory for the send and receive counts
  psi_send_count = (unsigned int *) mymalloc("psi_send_count", NTask * sizeof(unsigned int));
  psi_recv_count = (unsigned int *) mymalloc("psi_recv_count", NTask * sizeof(unsigned int));

  /* [1] fft setup */
  double k2, kx, ky, kz, K0;

  mpi_printf("Starting becdm kinetic fft computation\n");

  double tstart = second();

  /* allocate the memory to hold the FFT fields */
  densityfield = (fft_complex *) mymalloc("densityfield", PMGRID * sizeof(fft_complex));
  workspace = (fft_complex *) mymalloc("workspace", PMGRID * sizeof(fft_complex));

  int alignflag = 0;
  int stride = PMGRID;
  int ndim[1] = { PMGRID };     /* dimension of the 1D transforms */

  /* Set up the FFTW */
  myplan.forward_plan_zdir = FFTW(plan_many_dft) (1, ndim, 1, densityfield, 0, 1, PMGRID, workspace, 0, 1, PMGRID, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);
  myplan.forward_plan_ydir = FFTW(plan_many_dft) (1, ndim, 1, densityfield, 0, stride, PMGRID * PMGRID, workspace, 0, stride, PMGRID * PMGRID, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);
  myplan.forward_plan_xdir = FFTW(plan_many_dft) (1, ndim, 1, densityfield, 0, stride, PMGRID * PMGRID, workspace, 0, stride, PMGRID * PMGRID, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  myplan.backward_plan_xdir = FFTW(plan_many_dft) (1, ndim, 1, densityfield, 0, stride, PMGRID * PMGRID, workspace, 0, stride, PMGRID * PMGRID, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);
  myplan.backward_plan_ydir = FFTW(plan_many_dft) (1, ndim, 1, densityfield, 0, stride, PMGRID * PMGRID, workspace, 0, stride, PMGRID * PMGRID, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);
  myplan.backward_plan_zdir = FFTW(plan_many_dft) (1, ndim, 1, densityfield, 0, 1, PMGRID, workspace, 0, 1, PMGRID, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  myfree(workspace);
  myfree(densityfield);

  my_slab_based_fft_init(&myplan, PMGRID, PMGRID, PMGRID);

  size_t maxfftsize = imax(myplan.largest_x_slab * PMGRID, myplan.largest_y_slab * PMGRID) * ((size_t) PMGRID);

  /* allocate the memory to hold the FFT fields */
  psifield = (fft_complex *) mymalloc("psigrid", maxfftsize * sizeof(fft_complex));
  memset(psifield, 0, maxfftsize * sizeof(fft_complex));

  workspace = (fft_complex *) mymalloc("workspace", maxfftsize * sizeof(fft_complex));
  /* [2] make grid of psi for fft */
  particle2fftgrid(&myplan);

  /* [3] fourier transform it */

  K0 = 2.0 * M_PI / All.BoxSize;  /* minimum k */
  //K1 = K0 * PMGRID / 2.0;         /* maximum k */
  
  /* Do the FFT of the psi_field */
  my_slab_based_fft_c2c(&myplan, &psifield[0], &workspace[0], 1);
  fft_complex *fft_of_psifield = (fft_complex *) & psifield[0];

  /*  get timestep dt */
  double dt_drift;
  integertime ti_step, itstart, itend;
  int bin = -1;
  for(int p = 0; p < NumPart; p++) {
    if(P[p].Type == 1) {
      bin = P[p].TimeBinGrav;  // XXX nicer way to grab global time?
      assert(bin >= 0);
      break;
    }
  }
  ti_step = bin ? (((integertime) 1) << bin) : 0;
  itstart = All.Ti_begstep[bin];     /* beginning of step */
  itend = itstart + ti_step;     /* end of step */
  
  if(All.ComovingIntegrationOn)
    dt_drift = get_drift_factor(itstart, itend);  // dt_drift is int_{t}^{t+dt} dt/a^2
  else
    dt_drift = (itend - itstart) * All.Timebase_interval;

  /* [4] apply k^2 */
#pragma omp parallel for private(y, z)
  for(int x = 0; x < PMGRID; x++)
    for(int y = myplan.slabstart_y; y < myplan.slabstart_y + myplan.nslab_y; y++)
      for(int z = 0; z < PMGRID; z++)
        {
          kx = (x < PMGRID / 2) ? x : x - PMGRID;
          ky = (y < PMGRID / 2) ? y : y - PMGRID;
          kz = (z < PMGRID / 2) ? z : z - PMGRID;
  
          k2 = kx * kx + ky * ky + kz * kz;
  
          large_array_offset ip = ((large_array_offset) PMGRID) * (PMGRID * (y - myplan.slabstart_y) + x) + z;
  
          double psihatRe = fft_of_psifield[ip][0];
          double psihatIm = fft_of_psifield[ip][1];
          double cosx = cos(dt_drift * All.hbar / All.mAxion * K0 * K0 * 0.5 * k2);
          double sinx = sin(dt_drift * All.hbar / All.mAxion * K0 * K0 * 0.5 * k2);
          fft_of_psifield[ip][0] = (cosx * psihatRe +  sinx * psihatIm) / pow(PMGRID,3);  // evolves as  psi = exp(dt * (hbar/m) * (-1.i*((2*pi/Lbox)^2*double(kSq))/2)) * psi
          fft_of_psifield[ip][1] = (cosx * psihatIm -  sinx * psihatRe) / pow(PMGRID,3);  // also add in normalizaiton ofr ifft
        }

  /* [5] ifft the result */
  my_slab_based_fft_c2c(&myplan, &psifield[0], &workspace[0], -1);

  /* [6] get result back to particles */
  fftgrid2particle(&myplan);

  /* [7] cleanup */

  myfree(workspace);
  myfree(psifield);
  
  myfree(myplan.first_slab_y_of_task);
  myfree(myplan.slabs_y_per_task);
  myfree(myplan.first_slab_x_of_task);
  myfree(myplan.slabs_x_per_task);
  myfree(myplan.slab_to_task);
  
  myfree(psi_recv_count);
  myfree(psi_send_count);

  double tend = second();
  
  if(ThisTask == 0) {
    printf("end becdm fft. took %g seconds\n", timediff(tstart, tend));
    myflush(stdout);
  } 

  /* [8] Update the density/masses of the static grid cells */
  for(int p = 0; p < NumPart; p++) {
    if(P[p].Type == 1)
      P[p].Mass = (pow(P[p].PsiRe,2) + pow(P[p].PsiIm,2)) * pow(All.BoxSize/((double) PMGRID), 3);
  } 
}






void particle2fftgrid(fft_plan *myplan) {
  int i, p, x, y, z, correspTask, correspSlabStartY;
  
  // [1] add locals to psifield[][]
  for(p = 0; p < NumPart; p++) {
    if(P[p].Type == 1) {
      x = (int) (P[p].Pos[0] / ( All.BoxSize / PMGRID ) - 0.25);
      y = (int) (P[p].Pos[1] / ( All.BoxSize / PMGRID ) - 0.25);
      z = (int) (P[p].Pos[2] / ( All.BoxSize / PMGRID ) - 0.25);
        
      correspTask = y / (PMGRID/NTask);  
      correspSlabStartY = correspTask * (PMGRID/NTask);  //correspSlabStartY = myplan->slabstart_y; // also true, if local
      large_array_offset ip = ((large_array_offset) PMGRID) * (PMGRID * (y - correspSlabStartY) + x) + z;
      
      if(correspTask == ThisTask) {
        psifield[ip][0] = P[p].PsiRe;
        psifield[ip][1] = P[p].PsiIm;
      }
    }
  }

  // [2] send off non-locals
  /* count how many we have on each task */
  for(i = 0; i < NTask; i++) {
    psi_send_count[i] = 0;
    psi_recv_count[i] = 0;
  }

  for(p = 0; p < NumPart; p++) {
    if(P[p].Type == 1) {
      x = (int) (P[p].Pos[0] / ( All.BoxSize / PMGRID ) - 0.25);
      y = (int) (P[p].Pos[1] / ( All.BoxSize / PMGRID ) - 0.25);
      z = (int) (P[p].Pos[2] / ( All.BoxSize / PMGRID ) - 0.25);
        
      correspTask = y / (PMGRID/NTask);  
     
      if(correspTask != ThisTask)
        psi_send_count[correspTask]++;
    }
  }

  MPI_Alltoall(psi_send_count, 1, MPI_INT, psi_recv_count, 1, MPI_INT, MPI_COMM_WORLD);
     
  /* send off --*/
  unsigned int nimport = 0;
  unsigned int nexport = 0;
  unsigned int send_offset[NTask];
  unsigned int receive_offset[NTask];
  unsigned int j;

  //calculate nimport/nexport, send_offset/receive_offset
  for(j = 0, send_offset[0] = 0, receive_offset[0] = 0; j < NTask; j++) {
    nimport += psi_recv_count[j];
    nexport += psi_send_count[j];

    if(j > 0) {
      send_offset[j] = send_offset[j - 1] + psi_send_count[j - 1];
      receive_offset[j] = receive_offset[j - 1] + psi_recv_count[j - 1];
    }
  }

  psiType *psi_in = (psiType *) mymalloc("psi_in", nimport * sizeof(psiType));
  psiType *psi_out = (psiType *) mymalloc("psi_out", nexport * sizeof(psiType));

  //prepare data for export
  for(j = 0; j < NTask; j++)
    psi_send_count[j] = 0;

  unsigned int off;

  for(p = 0; p < NumPart; p++) {
    if(P[p].Type == 1) {
      x = (int) (P[p].Pos[0] / ( All.BoxSize / PMGRID ) - 0.25);
      y = (int) (P[p].Pos[1] / ( All.BoxSize / PMGRID ) - 0.25);
      z = (int) (P[p].Pos[2] / ( All.BoxSize / PMGRID ) - 0.25);
        
      correspTask = y / (PMGRID/NTask);  
      correspSlabStartY = correspTask * (PMGRID/NTask);
      large_array_offset ip = ((large_array_offset) PMGRID) * (PMGRID * (y - correspSlabStartY) + x) + z;
      
      if(correspTask != ThisTask) {
        off = send_offset[correspTask] + psi_send_count[correspTask]++;
        psi_out[off].psiRe = P[p].PsiRe;
        psi_out[off].psiIm = P[p].PsiIm;
        psi_out[off].ip = ip;
      }
    }
  }

  //exchange data
  unsigned int ngrp, receive_task;
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++) {
    //send_task = ThisTask;
    receive_task = ThisTask ^ ngrp;

    if(receive_task < NTask) {
      if(psi_send_count[receive_task] > 0 || psi_recv_count[receive_task] > 0) {
          MPI_Sendrecv(&psi_out[send_offset[receive_task]], psi_send_count[receive_task] * sizeof(psiType), MPI_BYTE, receive_task, 0,
                        &psi_in[receive_offset[receive_task]], psi_recv_count[receive_task] * sizeof(psiType), MPI_BYTE, receive_task, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
    }
  }

  // [3] loop over received data and add to psifield[][]
  //store the received psi
  for(j = 0; j < nimport; j++) {
    psifield[psi_in[j].ip][0] = psi_in[j].psiRe;
    psifield[psi_in[j].ip][1] = psi_in[j].psiIm;
  }

  myfree(psi_out);
  myfree(psi_in);
  
}








void fftgrid2particle(fft_plan *myplan) {
  int i, p, x, y, z, correspTask, correspSlabStartY;

  // [1] add psifield[][] local particles
  for(p = 0; p < NumPart; p++) {
    if(P[p].Type == 1) {
      x = (int) (P[p].Pos[0] / ( All.BoxSize / PMGRID ) - 0.25);
      y = (int) (P[p].Pos[1] / ( All.BoxSize / PMGRID ) - 0.25);
      z = (int) (P[p].Pos[2] / ( All.BoxSize / PMGRID ) - 0.25);
        
      correspTask = y / (PMGRID/NTask);  
      correspSlabStartY = correspTask * (PMGRID/NTask);  //correspSlabStartY = myplan->slabstart_y; // also true, if local
      large_array_offset ip = ((large_array_offset) PMGRID) * (PMGRID * (y - correspSlabStartY) + x) + z;
      
      if(correspTask == ThisTask) {
        P[p].PsiRe = psifield[ip][0];
        P[p].PsiIm = psifield[ip][1];
      }
    }
  }

  // [2] send non-local data (i.e., requests for the fftgrid data by indix ip)
  /* count how many we have on each task */
  for(i = 0; i < NTask; i++) {
    psi_send_count[i] = 0;
    psi_recv_count[i] = 0;
  }

  for(p = 0; p < NumPart; p++) {
    if(P[p].Type == 1) {
      x = (int) (P[p].Pos[0] / ( All.BoxSize / PMGRID ) - 0.25);
      y = (int) (P[p].Pos[1] / ( All.BoxSize / PMGRID ) - 0.25);
      z = (int) (P[p].Pos[2] / ( All.BoxSize / PMGRID ) - 0.25);
        
      correspTask = y / (PMGRID/NTask);  
     
      if(correspTask != ThisTask)
        psi_send_count[correspTask]++;
    }
  }

  MPI_Alltoall(psi_send_count, 1, MPI_INT, psi_recv_count, 1, MPI_INT, MPI_COMM_WORLD);
     
  /* send off --*/
  unsigned int nimport = 0;
  unsigned int nexport = 0;
  unsigned int send_offset[NTask];
  unsigned int receive_offset[NTask];
  unsigned int j;

  //calculate nimport/nexport, send_offset/receive_offset
  for(j = 0, send_offset[0] = 0, receive_offset[0] = 0; j < NTask; j++) {
    nimport += psi_recv_count[j];
    nexport += psi_send_count[j];

    if(j > 0) {
      send_offset[j] = send_offset[j - 1] + psi_send_count[j - 1];
      receive_offset[j] = receive_offset[j - 1] + psi_recv_count[j - 1];
    }
  }

  psiType *psi_in = (psiType *) mymalloc("psi_in", nimport * sizeof(psiType));
  psiType *psi_out = (psiType *) mymalloc("psi_out", nexport * sizeof(psiType));

  //prepare data for export
  for(j = 0; j < NTask; j++)
    psi_send_count[j] = 0;

  unsigned int off;

  for(p = 0; p < NumPart; p++) {
    if(P[p].Type == 1) {
      x = (int) (P[p].Pos[0] / ( All.BoxSize / PMGRID ) - 0.25);
      y = (int) (P[p].Pos[1] / ( All.BoxSize / PMGRID ) - 0.25);
      z = (int) (P[p].Pos[2] / ( All.BoxSize / PMGRID ) - 0.25);
        
      correspTask = y / (PMGRID/NTask);  
      correspSlabStartY = correspTask * (PMGRID/NTask);
      large_array_offset ip = ((large_array_offset) PMGRID) * (PMGRID * (y - correspSlabStartY) + x) + z;
      
      if(correspTask != ThisTask) {
        off = send_offset[correspTask] + psi_send_count[correspTask]++;
        psi_out[off].psiRe = P[p].PsiRe;
        psi_out[off].psiIm = P[p].PsiIm;
        psi_out[off].ip = ip;
      }
    }
  }

  //exchange data
  unsigned int ngrp, receive_task;
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++) {
    //send_task = ThisTask;
    receive_task = ThisTask ^ ngrp;

    if(receive_task < NTask) {
      if(psi_send_count[receive_task] > 0 || psi_recv_count[receive_task] > 0) {
          MPI_Sendrecv(&psi_out[send_offset[receive_task]], psi_send_count[receive_task] * sizeof(psiType), MPI_BYTE, receive_task, 0,
                        &psi_in[receive_offset[receive_task]], psi_recv_count[receive_task] * sizeof(psiType), MPI_BYTE, receive_task, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
    }
  }

  // [3] loop over received data and grab the corresponding fftgrid[][] data into psi_in
  //store the received psi
  for(j = 0; j < nimport; j++) {
    psi_in[j].psiRe = psifield[psi_in[j].ip][0];
    psi_in[j].psiIm = psifield[psi_in[j].ip][1];
  }

  // [4] exchange psi_in back to local processors
  //exchange data
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++) {
    receive_task = ThisTask ^ ngrp;

    if(receive_task < NTask) {
      if(psi_recv_count[receive_task] > 0 || psi_send_count[receive_task] > 0) {
          MPI_Sendrecv(&psi_in[receive_offset[receive_task]], psi_recv_count[receive_task] * sizeof(psiType), MPI_BYTE, receive_task, 0,
                        &psi_out[send_offset[receive_task]], psi_send_count[receive_task] * sizeof(psiType), MPI_BYTE, receive_task, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
    }
  }
  
  // [5] add grabbed data to particles
  for(j = 0; j < NTask; j++)
    psi_send_count[j] = 0;
  
  for(p = 0; p < NumPart; p++) {
    if(P[p].Type == 1) {
      x = (int) (P[p].Pos[0] / ( All.BoxSize / PMGRID ) - 0.25);
      y = (int) (P[p].Pos[1] / ( All.BoxSize / PMGRID ) - 0.25);
      z = (int) (P[p].Pos[2] / ( All.BoxSize / PMGRID ) - 0.25);
        
      correspTask = y / (PMGRID/NTask);  
      
      if(correspTask != ThisTask) {
        off = send_offset[correspTask] + psi_send_count[correspTask]++;
        P[p].PsiRe = psi_out[off].psiRe;
        P[p].PsiIm = psi_out[off].psiIm;
      }
    }
  }
  
  myfree(psi_out);
  myfree(psi_in);  
  
}


#endif
