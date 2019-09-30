/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/turb/mhd_seedpspec.c
 * \date        01/2016
 * \author      Philip Mocz
 * \brief       Sets power-spectrum based tangled initial magnetic fields (3D)
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_randist.h>

#include "../allvars.h"
#include "../proto.h"

#ifdef MHD_SEEDPSPEC

gsl_rng *StRng;

//Fourier Space variables
double *StAhatRe;
double *StAhatIm;
double *StMode;
int StNModes;
double Nk;


/** Initialize the Fourier Modes*/
void init_pspec_bfld() {
 
  int ikx, iky, ikz, i, j;
  double kx, ky, kz, k;
  double ampl;

  int ikxmax = All.B_spec_kmax;//128;//256;
  int ikymax = ikxmax;
  int ikzmax = ikxmax;
  double StKmin = 1.;
  double StKmax = 1.*ikxmax;
  Nk = 2.0 * StKmax;

  StRng = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(StRng, 42);

  StNModes = 0;
  for(ikx = 0; ikx <= ikxmax; ikx++) {
      kx = 1. * ikx;
      for(iky = 0; iky <= ikymax; iky++) {
          ky = 1. * iky;
          for(ikz = 0; ikz <= ikzmax; ikz++) {
              kz = 1. * ikz;
              k = sqrt(kx * kx + ky * ky + kz * kz);
              if(k >= StKmin && k <= StKmax) {              
                  StNModes += 4;
              }
          }
      }
  }

  StMode   = (double *) mymalloc_movable(&StMode,   "StModes",  StNModes * 3 * sizeof(double));
  StAhatRe = (double *) mymalloc_movable(&StAhatRe, "StAhatRe", StNModes * 3 * sizeof(double));
  StAhatIm = (double *) mymalloc_movable(&StAhatIm, "StAhatIm", StNModes * 3 * sizeof(double));

  // set k vectors
  StNModes = 0;
  for(ikx = 0; ikx <= ikxmax; ikx++) {
      kx = 1. * ikx;
      for(iky = 0; iky <= ikymax; iky++) {
          ky = 1. * iky;
          for(ikz = 0; ikz <= ikzmax; ikz++) {
              kz = 1. * ikz;
              k = sqrt(kx * kx + ky * ky + kz * kz);
              if(k >= StKmin && k <= StKmax) {
                  StMode[3 * StNModes + 0] = kx;
                  StMode[3 * StNModes + 1] = ky;
                  StMode[3 * StNModes + 2] = kz;
                  StNModes++;
                  
                  StMode[3 * StNModes + 0] = kx;
                  StMode[3 * StNModes + 1] = -ky;
                  StMode[3 * StNModes + 2] = kz;
                  StNModes++;
                  
                  StMode[3 * StNModes + 0] = kx;
                  StMode[3 * StNModes + 1] = ky;
                  StMode[3 * StNModes + 2] = -kz;
                  StNModes++;
                  
                  StMode[3 * StNModes + 0] = kx;
                  StMode[3 * StNModes + 1] = -ky;
                  StMode[3 * StNModes + 2] = -kz;
                  StNModes++;
              }
          }
      }
  }
  
  // set Ahat
  for(i = 0; i < StNModes; i++)
    {
      kx = StMode[3 * i + 0];
      ky = StMode[3 * i + 1];
      kz = StMode[3 * i + 2];
      k = sqrt(kx * kx + ky * ky + kz * kz);
      ampl = 1.0 / sqrt(2.0/M_PI) * sqrt(All.B_pspec_ampl/3.0 * pow(k,All.B_pspec_slope-2.0) * exp(-pow(k/All.B_pspec_kcut,4))) * pow(Nk,3) / (2.0 * M_PI);
      
      if(All.B_pspec_helical == 0.0)  //  no helicity
        {
          for(j = 0; j < 3; j++)
            {
              double phase = gsl_rng_uniform(StRng);
              double amplP = gsl_ran_gaussian(StRng, 1.0);
              amplP *= ampl;
              StAhatRe[3 * i + j] = amplP * cos(2. * M_PI * phase);
              StAhatIm[3 * i + j] = amplP * sin(2. * M_PI * phase);
            }
        }
      else if(All.B_pspec_helical == 1.0)  //  max helicity
        {
          double zhat_x = .1;
          double zhat_y = .4;
          double zhat_z = 0.911043357914430;
          
          //e1hat = cross(kfull,zhat) / norm( cross(kfull,zhat) );
          //e2hat = cross(kfull, cross(kfull,zhat)) / norm(  cross(kfull, cross(kfull,zhat)) );
          //eP = (e1hat + 1i * e2hat) / sqrt(2);
          double e1hat[3], e2hat[3], ePRe[3], ePIm[3], nm;
          e1hat[0] = ky * zhat_z - kz * zhat_y;
          e1hat[1] = kz * zhat_x - kx * zhat_z;
          e1hat[2] = kx * zhat_y - ky * zhat_x;
          nm = sqrt(e1hat[0]*e1hat[0] + e1hat[1]*e1hat[1] + e1hat[2]*e1hat[2]);
          e1hat[0] /= nm*sqrt(2.);
          e1hat[1] /= nm*sqrt(2.);
          e1hat[2] /= nm*sqrt(2.);
          
          e2hat[0] = ky * e1hat[2] - kz * e1hat[1];
          e2hat[1] = kz * e1hat[0] - kx * e1hat[2];
          e2hat[2] = kx * e1hat[1] - ky * e1hat[0]; 
          nm = sqrt(e2hat[0]*e2hat[0] + e2hat[1]*e2hat[1] + e2hat[2]*e2hat[2]);
          e2hat[0] /= nm*sqrt(2.);
          e2hat[1] /= nm*sqrt(2.);
          e2hat[2] /= nm*sqrt(2.); 
           
          for(j = 0; j < 3; j++)
            {
              double phase = gsl_rng_uniform(StRng);
              double amplP = gsl_ran_gaussian(StRng, 1.0);
              amplP *= sqrt(3) * ampl; 
              StAhatRe[3 * i + j] = amplP * e1hat[j];
              StAhatIm[3 * i + j] = amplP * e2hat[j]
            }
        }
      else
        {
          terminate("unknown helicity setting");
        }
    }
  mpi_printf("B-fld power spectra init done.\n");
}

void free_pspec_bfld() {
  myfree_movable(StAhatIm);
  myfree_movable(StAhatRe);
  myfree_movable(StMode);
  mpi_printf("B-fld power spectra cells done.\n");
}


/** Set B-field of cell i from fourier modes*/
void set_pspec_bfld(int i) {
  int m, j;
  double B_value[3];
  B_value[0] = 0;
  B_value[1] = 0;
  B_value[2] = 0;
  
  for(m = 0; m < StNModes; m++)
  {
    
    double kx = StMode[3 * m + 0];
    double ky = StMode[3 * m + 1];
    double kz = StMode[3 * m + 2];
    double kdotx = kx * P[i].Pos[0] / boxSize_X  +  ky * P[i].Pos[1] / boxSize_Y  +  kz * P[i].Pos[2] / boxSize_Z;
    
    double realt = cos(2. * M_PI * kdotx);
    double imagt = sin(2. * M_PI * kdotx);

    //for(j = 0; j < 3; j++) {
    //  A_value[j] += 2.0 * (StAhatRe[3 * m + j] * realt - StAhatIm[3 * m + j] * imagt) / pow(Nk,3); // contribution from k and -k
    //}
    B_value[0] += (ky*StAhatIm[3 * m + 2] - kz*StAhatIm[3 * m + 1]) * realt + (ky*StAhatRe[3 * m + 2] - kz*StAhatRe[3 * m + 1]) * imagt;
    B_value[1] += (kz*StAhatIm[3 * m + 0] - kx*StAhatIm[3 * m + 2]) * realt + (kz*StAhatRe[3 * m + 0] - kx*StAhatRe[3 * m + 2]) * imagt;
    B_value[2] += (kx*StAhatIm[3 * m + 1] - ky*StAhatIm[3 * m + 0]) * realt + (kx*StAhatRe[3 * m + 1] - ky*StAhatRe[3 * m + 0]) * imagt;
  }
  
  double bfac = 1. / (sqrt(All.UnitMass_in_g / All.UnitLength_in_cm) / (All.UnitTime_in_s / All.HubbleParam));
  for(j = 0; j < 3; j++) {
    SphP[i].BConserved[j] = (-2.0 * (2.0 * M_PI) / pow(Nk,3)) * B_value[j] * SphP[i].Volume * bfac;// / sqrt(4. * M_PI);   /* convert Gauss-cgs to heavyside - lorentz done later!*/
    SphP[i].B[j] = SphP[i].BConserved[j] / SphP[i].Volume;
  }
}


void set_pscpec_ICs() {
  init_pspec_bfld();
  for(int i = 0; i < NumGas; i++)
      set_pspec_bfld(i);
  free_pspec_bfld();
}  


#endif



