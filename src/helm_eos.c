/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/helm_eos.c
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

#include "arepoconfig.h"
#include "proto.h"

#ifdef EOS_DEGENERATE

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <assert.h>

#include <gsl/gsl_const.h>
#include <gsl/gsl_math.h>

#include "helm_eos.h"

#define ELECTRON_CHARGE_ESU (4.80320427e-10)

#ifndef HELM_EOS_MAXITER
#define HELM_EOS_MAXITER 100
#endif // HELM_EOS_MAXITER

#ifndef HELM_EOS_EPS
#define HELM_EOS_EPS 1.0e-13
#endif // HELM_EOS_EPS

#define IMAX 271
#define JMAX 101

typedef double helm_eos_table_entry[IMAX][JMAX];

struct helm_eos_table
{
  // number of entries
  int ntemp, nrho;

  // density and temperature ranges
  double ltempMin, ltempMax;
  double lrhoMin, lrhoMax;
  double ltempDelta;
  double lrhoDelta;
  double tempDelta;
  double rhoDelta;
  double tempMin, tempMax;
  double rhoMin, rhoMax;
  double temp[JMAX];
  double rho[IMAX];


  // d means derivative w. r. t. density; t means derivative w. r. t. temperature

  // Helmholtz free energy
  helm_eos_table_entry f;
  helm_eos_table_entry fd;
  helm_eos_table_entry ft;
  helm_eos_table_entry fdd;
  helm_eos_table_entry ftt;
  helm_eos_table_entry fdt;
  helm_eos_table_entry fddt;
  helm_eos_table_entry fdtt;
  helm_eos_table_entry fddtt;

  // pressure derivative w. r. t. density
  helm_eos_table_entry dpdf;
  helm_eos_table_entry dpdfd;
  helm_eos_table_entry dpdft;
  helm_eos_table_entry dpdfdd;
  helm_eos_table_entry dpdftt;
  helm_eos_table_entry dpdfdt;

  // chemical potential
  helm_eos_table_entry ef;
  helm_eos_table_entry efd;
  helm_eos_table_entry eft;
  helm_eos_table_entry efdd;
  helm_eos_table_entry eftt;
  helm_eos_table_entry efdt;

  // number density
  helm_eos_table_entry xf;
  helm_eos_table_entry xfd;
  helm_eos_table_entry xft;
  helm_eos_table_entry xfdd;
  helm_eos_table_entry xftt;
  helm_eos_table_entry xfdt;


  // species information
  int nspecies;
  double *na;
  double *nai;
  double *nz;
};

// global pointer to the EOS table
static struct helm_eos_table *helm_eos_table;


// internal functions
#ifndef _GNU_SOURCE
//static double exp10(double);
static inline double exp10(double x)
{
  return exp(log(10) * x);
}
#endif // _GNU_SOURCE

// quintic hermite polynomials

// psi0 and its derivatives
static inline double psi0(double z)
{
  return z * z * z * (z * (-6.0 * z + 15.0) - 10.0) + 1.0;
}

static inline double dpsi0(double z)
{
  return z * z * (z * (-30.0 * z + 60.0) - 30.0);
}

static inline double ddpsi0(double z)
{
  return z * (z * (-120.0 * z + 180.0) - 60.0);
}

// psi1 and its derivatives
static inline double psi1(double z)
{
  return z * (z * z * (z * (-3.0 * z + 8.0) - 6.0) + 1.0);
}

static inline double dpsi1(double z)
{
  return z * z * (z * (-15.0 * z + 32.0) - 18.0) + 1.0;
}

static inline double ddpsi1(double z)
{
  return z * (z * (-60.0 * z + 96.0) - 36.0);
}

// psi2 and its derivatives
static inline double psi2(double z)
{
  return 0.5 * z * z * (z * (z * (-z + 3.0) - 3.0) + 1.0);
}

static inline double dpsi2(double z)
{
  return 0.5 * z * (z * (z * (-5.0 * z + 12.0) - 9.0) + 2.0);
}

static inline double ddpsi2(double z)
{
  return 0.5 * (z * (z * (-20.0 * z + 36.0) - 18.0) + 2.0);
}

// biquintic hermite polynomial
static inline double h5(double fi[36], double w0t, double w1t, double w2t, double w0mt, double w1mt, double w2mt, double w0d, double w1d, double w2d, double w0md, double w1md, double w2md)
{
  return fi[0] * w0d * w0t + fi[1] * w0md * w0t + fi[2] * w0d * w0mt + fi[3] * w0md * w0mt +
    fi[4] * w0d * w1t + fi[5] * w0md * w1t + fi[6] * w0d * w1mt + fi[7] * w0md * w1mt + fi[8] * w0d * w2t +
    fi[9] * w0md * w2t + fi[10] * w0d * w2mt + fi[11] * w0md * w2mt + fi[12] * w1d * w0t +
    fi[13] * w1md * w0t + fi[14] * w1d * w0mt + fi[15] * w1md * w0mt + fi[16] * w2d * w0t +
    fi[17] * w2md * w0t + fi[18] * w2d * w0mt + fi[19] * w2md * w0mt + fi[20] * w1d * w1t +
    fi[21] * w1md * w1t + fi[22] * w1d * w1mt + fi[23] * w1md * w1mt + fi[24] * w2d * w1t +
    fi[25] * w2md * w1t + fi[26] * w2d * w1mt + fi[27] * w2md * w1mt + fi[28] * w1d * w2t +
    fi[29] * w1md * w2t + fi[30] * w1d * w2mt + fi[31] * w1md * w2mt + fi[32] * w2d * w2t + fi[33] * w2md * w2t + fi[34] * w2d * w2mt + fi[35] * w2md * w2mt;
}

// cubic hermite polynomial
// psi0 and its derivatives
static inline double xpsi0(double z)
{
  return z * z * (2.0 * z - 3.0) + 1.0;
}

static inline double xdpsi0(double z)
{
  return z * (6.0 * z - 6.0);
}

// psi1 and its derivatives
static inline double xpsi1(double z)
{
  return z * (z * (z - 2.0) + 1.0);
}

static inline double xdpsi1(double z)
{
  return z * (3.0 * z - 4.0) + 1.0;
}

// bicubic hermite polynomial
static inline double h3(double fi[16], double w0t, double w1t, double w0mt, double w1mt, double w0d, double w1d, double w0md, double w1md)
{
  return fi[0] * w0d * w0t + fi[1] * w0md * w0t + fi[2] * w0d * w0mt + fi[3] * w0md * w0mt +
    fi[4] * w0d * w1t + fi[5] * w0md * w1t + fi[6] * w0d * w1mt + fi[7] * w0md * w1mt + fi[8] * w1d * w0t +
    fi[9] * w1md * w0t + fi[10] * w1d * w0mt + fi[11] * w1md * w0mt + fi[12] * w1d * w1t + fi[13] * w1md * w1t + fi[14] * w1d * w1mt + fi[15] * w1md * w1mt;
}

static inline void azbar(const double xnuc[], double *abar, double *zbar)
{
  const double *nai = helm_eos_table->nai, *nz = helm_eos_table->nz;
  int nspecies = helm_eos_table->nspecies;
  double xsum = 0.0;

  // compute abar, zbar
  *abar = 0.0;
  *zbar = 0.0;

  {
    int i;
    for(i = 0; i < nspecies; i++)
      {
        const double ymass = xnuc[i] * nai[i];
        *abar += ymass;
        *zbar += nz[i] * ymass;
        xsum += xnuc[i];
      }
  }

  *abar = xsum / *abar;
  *zbar = *zbar / xsum * *abar;
}

// different contributions to the EOS
// the five-element arrays correspond to the actual quantity and its derivatives w. r. t. density, temperature, abar, and zbar

// radiation
static int __attribute__ ((unused)) helm_eos_rad(double rho, double temp, double abar, double zbar, double prad[5], double erad[5], double srad[5])
{
  const double rhoi = 1.0 / rho, tempi = 1.0 / temp;
  // prad
  prad[0] = 4.0 / 3.0 * GSL_CONST_CGS_STEFAN_BOLTZMANN_CONSTANT / GSL_CONST_CGS_SPEED_OF_LIGHT * gsl_pow_4(temp);
  // dprad dd
  prad[1] = 0.0;
  // dprad dt
  prad[2] = 4.0 * prad[0] * tempi;
  // dprad da
  prad[3] = 0.0;
  // dprad dz
  prad[4] = 0.0;

  // erad
  erad[0] = 3.0 * prad[0] * rhoi;
  // derad dd
  erad[1] = -erad[0] * rhoi;
  // derad dt
  erad[2] = 3.0 * prad[2] * rhoi;
  // derad da
  erad[3] = 0.0;
  // derad dz
  erad[4] = 0.0;

  // srad
  srad[0] = (prad[0] * rhoi + erad[0]) * tempi;
  // dsrad dd
  srad[1] = ((prad[1] - prad[0] * rhoi) * rhoi + erad[1]) * tempi;
  // dsrad dt
  srad[2] = (prad[2] * rhoi + erad[2] - srad[0]) * tempi;
  // dsrad da
  srad[3] = 0.0;
  // dsrad dz
  srad[4] = 0.0;

  // independent of zbar and abar; tell the compiler that this is intentional
  (void) zbar;
  (void) abar;

  return 0;
}


// ion section
static int __attribute__ ((unused)) helm_eos_ion(double rho, double temp, double ltemp, const struct helm_eos_cache *cache, double abar, double zbar, double pion[5], double eion[5], double sion[5])
{
  const double kt = GSL_CONST_CGS_BOLTZMANN * temp;
  const double ytot = cache->ytot;
  const double xni = cache->xni;
  const double dxnidd = cache->dxnidd;
  const double dxnida = cache->dxnida;
  //const double s = (2.0 * M_PI * GSL_CONST_CGS_UNIFIED_ATOMIC_MASS * GSL_CONST_CGS_BOLTZMANN) / (GSL_CONST_CGS_PLANCKS_CONSTANT_H * GSL_CONST_CGS_PLANCKS_CONSTANT_H) * temp;
  //const double y = log(abar * abar * sqrt(abar) / (rho * GSL_CONST_NUM_AVOGADRO) * s * sqrt(s));
  const double y = cache->ywot + cache->lswot15 * ltemp;

  // pion
  pion[0] = xni * kt;
  // dpion dd
  pion[1] = dxnidd * kt;
  // dpion dt
  pion[2] = xni * GSL_CONST_CGS_BOLTZMANN;
  // dpion da
  pion[3] = dxnida * kt;
  // dpion dz
  pion[4] = 0.0;

  // eion
  eion[0] = 1.5 * pion[0] / rho;
  // deion dd
  eion[1] = (1.5 * pion[1] - eion[0]) / rho;
  // deion dt
  eion[2] = 1.5 * pion[2] / rho;
  // deion da
  eion[3] = 1.5 * pion[3] / rho;
  // deion dz
  eion[4] = 0.0;

  // sion
  sion[0] = (pion[0] / rho + eion[0]) / temp + GSL_CONST_CGS_BOLTZMANN * GSL_CONST_NUM_AVOGADRO * ytot * y;
  // dsion dd
  sion[1] = (pion[1] / rho - pion[0] / (rho * rho) + eion[1]) / temp - GSL_CONST_CGS_BOLTZMANN * GSL_CONST_NUM_AVOGADRO * ytot / rho;
  // dsion dt
  sion[2] = (pion[2] / rho + eion[2]) / temp - (pion[0] / rho + eion[0]) / (temp * temp) + 1.5 * GSL_CONST_CGS_BOLTZMANN * GSL_CONST_NUM_AVOGADRO * ytot / temp;
  // dsion da
  sion[3] = (pion[3] / rho + eion[3]) / temp + GSL_CONST_CGS_BOLTZMANN * GSL_CONST_NUM_AVOGADRO * ytot * ytot * (2.5 - y);
  // dsion dz
  sion[4] = 0.0;

  // independent of zbar; tell the compiler, this is intentional
  (void) zbar;

  return 0;
}

// electron-positron
static int
  __attribute__ ((unused)) helm_eos_ele(double rho, double temp, double l10temp, double pele[5],
                                        double eele[5], double sele[5], double etaele[5], double xne[5], const struct helm_eos_cache *cache, const int only_e)
{
  const double ye = cache->ye;  // electron number fraction
  const double ytot = cache->ytot;
  const double din = cache->din;
  int iat, jat;                 // temperature and density indices in the table
  double fi[36];                // cache for the table values
  // temperature and density deltas
  double dth, dt2, dti, dt2i;
  double dd, dd2, ddi, dd2i;
  // various differences
  double xt, xd, mxt, mxd;
  // the six density and six temperature basis functions
  double si0t, si1t, si2t, si0mt, si1mt, si2mt;
  double si0d, si1d, si2d, si0md, si1md, si2md;
  // derivatives of the weight functions
  double dsi0t, dsi1t, dsi2t, dsi0mt, dsi1mt, dsi2mt;
  double dsi0d, dsi1d, dsi2d, dsi0md, dsi1md, dsi2md;
  // second derivatives
  double ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt;
  double ddsi0d, ddsi1d, ddsi2d, ddsi0md, ddsi1md, ddsi2md;
  // free energy and its derivatives
  double free_en, df_d, df_t, df_tt, df_dt, df_dd;
  // scratch variables
  double x, s;

  if(rho < helm_eos_table->rhoMin || rho > helm_eos_table->rhoMax)
    {
      fprintf(stderr, "density (%g) out of table [%g:%g]\n", rho, helm_eos_table->rhoMin, helm_eos_table->rhoMax);
      return -1;
    }
  if(temp < helm_eos_table->tempMin || temp > helm_eos_table->tempMax)
    {
      fprintf(stderr, "temperature (%g) out of table [%g:%g]\n", temp, helm_eos_table->tempMin, helm_eos_table->tempMax);
      return -1;
    }

  // hash locate temperature and density
  jat = (l10temp - helm_eos_table->ltempMin) / helm_eos_table->ltempDelta;
  iat = (cache->ldin - helm_eos_table->lrhoMin) / helm_eos_table->lrhoDelta;

  if(jat < 0 || jat >= helm_eos_table->ntemp)
    {
      fprintf(stderr, "temperature (%g) off table [%d:%d]\n", temp, 0, helm_eos_table->ntemp);
      return -1;
    }
  if(iat < 0 || iat >= helm_eos_table->nrho)
    {
      fprintf(stderr, "density (%g) off table [%d:%d], %g %g %g %g %d\n", rho, 0, helm_eos_table->nrho, cache->abar, cache->zbar, cache->ldin, cache->ye, iat);
      return -1;
    }

  // compute temperature and density deltas
  dth = helm_eos_table->temp[jat + 1] - helm_eos_table->temp[jat];
  dt2 = dth * dth;
  dti = 1.0 / dth;
  dt2i = 1.0 / dt2;
  dd = helm_eos_table->rho[iat + 1] - helm_eos_table->rho[iat];
  dd2 = dd * dd;
  ddi = 1.0 / dd;
  dd2i = 1.0 / dd2;

  // access the table locations only once
  fi[0] = helm_eos_table->f[iat][jat];
  fi[1] = helm_eos_table->f[iat + 1][jat];
  fi[2] = helm_eos_table->f[iat][jat + 1];
  fi[3] = helm_eos_table->f[iat + 1][jat + 1];
  fi[4] = helm_eos_table->ft[iat][jat];
  fi[5] = helm_eos_table->ft[iat + 1][jat];
  fi[6] = helm_eos_table->ft[iat][jat + 1];
  fi[7] = helm_eos_table->ft[iat + 1][jat + 1];
  fi[8] = helm_eos_table->ftt[iat][jat];
  fi[9] = helm_eos_table->ftt[iat + 1][jat];
  fi[10] = helm_eos_table->ftt[iat][jat + 1];
  fi[11] = helm_eos_table->ftt[iat + 1][jat + 1];
  fi[12] = helm_eos_table->fd[iat][jat];
  fi[13] = helm_eos_table->fd[iat + 1][jat];
  fi[14] = helm_eos_table->fd[iat][jat + 1];
  fi[15] = helm_eos_table->fd[iat + 1][jat + 1];
  fi[16] = helm_eos_table->fdd[iat][jat];
  fi[17] = helm_eos_table->fdd[iat + 1][jat];
  fi[18] = helm_eos_table->fdd[iat][jat + 1];
  fi[19] = helm_eos_table->fdd[iat + 1][jat + 1];
  fi[20] = helm_eos_table->fdt[iat][jat];
  fi[21] = helm_eos_table->fdt[iat + 1][jat];
  fi[22] = helm_eos_table->fdt[iat][jat + 1];
  fi[23] = helm_eos_table->fdt[iat + 1][jat + 1];
  fi[24] = helm_eos_table->fddt[iat][jat];
  fi[25] = helm_eos_table->fddt[iat + 1][jat];
  fi[26] = helm_eos_table->fddt[iat][jat + 1];
  fi[27] = helm_eos_table->fddt[iat + 1][jat + 1];
  fi[28] = helm_eos_table->fdtt[iat][jat];
  fi[29] = helm_eos_table->fdtt[iat + 1][jat];
  fi[30] = helm_eos_table->fdtt[iat][jat + 1];
  fi[31] = helm_eos_table->fdtt[iat + 1][jat + 1];
  fi[32] = helm_eos_table->fddtt[iat][jat];
  fi[33] = helm_eos_table->fddtt[iat + 1][jat];
  fi[34] = helm_eos_table->fddtt[iat][jat + 1];
  fi[35] = helm_eos_table->fddtt[iat + 1][jat + 1];

  // various differences
  xt = dmax(0.0, (temp - helm_eos_table->temp[jat]) * dti);
  xd = dmax(0.0, (din - helm_eos_table->rho[iat]) * ddi);
  mxt = 1.0 - xt;
  mxd = 1.0 - xd;

  // the six density and six temperature basis functions
  si0t = psi0(xt);
  si1t = psi1(xt) * dth;
  si2t = psi2(xt) * dt2;

  si0mt = psi0(mxt);
  si1mt = -psi1(mxt) * dth;
  si2mt = psi2(mxt) * dt2;

  si0d = psi0(xd);
  si1d = psi1(xd) * dd;
  si2d = psi2(xd) * dd2;

  si0md = psi0(mxd);
  si1md = -psi1(mxd) * dd;
  si2md = psi2(mxd) * dd2;

  // derivatives of the weight functions
  dsi0t = dpsi0(xt) * dti;
  dsi1t = dpsi1(xt);
  dsi2t = dpsi2(xt) * dth;

  dsi0mt = -dpsi0(mxt) * dti;
  dsi1mt = dpsi1(mxt);
  dsi2mt = -dpsi2(mxt) * dth;

  dsi0d = dpsi0(xd) * ddi;
  dsi1d = dpsi1(xd);
  dsi2d = dpsi2(xd) * dd;

  dsi0md = -dpsi0(mxd) * ddi;
  dsi1md = dpsi1(mxd);
  dsi2md = -dpsi2(mxd) * dd;

  // second derivatives of the weight functions
  ddsi0t = ddpsi0(xt) * dt2i;
  ddsi1t = ddpsi1(xt) * dti;
  ddsi2t = ddpsi2(xt);

  ddsi0mt = ddpsi0(mxt) * dt2i;
  ddsi1mt = -ddpsi1(mxt) * dti;
  ddsi2mt = ddpsi2(mxt);

  ddsi0d = ddpsi0(xd) * dd2i;
  ddsi1d = ddpsi1(xd) * ddi;
  ddsi2d = ddpsi2(xd);

  ddsi0md = ddpsi0(mxd) * dd2i;
  ddsi1md = -ddpsi1(mxd) * ddi;
  ddsi2md = ddpsi2(mxd);

  // the free energy
  free_en = h5(fi, si0t, si1t, si2t, si0mt, si1mt, si2mt, si0d, si1d, si2d, si0md, si1md, si2md);

  // derivative with respect to temperature
  df_t = h5(fi, dsi0t, dsi1t, dsi2t, dsi0mt, dsi1mt, dsi2mt, si0d, si1d, si2d, si0md, si1md, si2md);

  // derivative with respect to temperature**2
  df_tt = h5(fi, ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt, si0d, si1d, si2d, si0md, si1md, si2md);

  if(!only_e)
    {
      // derivative with respect to density
      df_d = h5(fi, si0t, si1t, si2t, si0mt, si1mt, si2mt, dsi0d, dsi1d, dsi2d, dsi0md, dsi1md, dsi2md);

      // derivative with respect to density**2
      df_dd = h5(fi, si0t, si1t, si2t, si0mt, si1mt, si2mt, ddsi0d, ddsi1d, ddsi2d, ddsi0md, ddsi1md, ddsi2md);

      // derivative with respect to temperature and density
      df_dt = h5(fi, dsi0t, dsi1t, dsi2t, dsi0mt, dsi1mt, dsi2mt, dsi0d, dsi1d, dsi2d, dsi0md, dsi1md, dsi2md);


      // now get the pressure derivative with density, chemical potential, and
      // electron positron number densities
      // get the interpolation weight functions
      si0t = xpsi0(xt);
      si1t = xpsi1(xt) * dth;

      si0mt = xpsi0(mxt);
      si1mt = -xpsi1(mxt) * dth;

      si0d = xpsi0(xd);
      si1d = xpsi1(xd) * dd;

      si0md = xpsi0(mxd);
      si1md = -xpsi1(mxd) * dd;


      // derivatives of weight functions
      dsi0t = xdpsi0(xt) * dti;
      dsi1t = xdpsi1(xt);

      dsi0mt = -xdpsi0(mxt) * dti;
      dsi1mt = xdpsi1(mxt);

      dsi0d = xdpsi0(xd) * ddi;
      dsi1d = xdpsi1(xd);

      dsi0md = -xdpsi0(mxd) * ddi;
      dsi1md = xdpsi1(mxd);


      // look in the pressure derivative only once
      fi[0] = helm_eos_table->dpdf[iat][jat];
      fi[1] = helm_eos_table->dpdf[iat + 1][jat];
      fi[2] = helm_eos_table->dpdf[iat][jat + 1];
      fi[3] = helm_eos_table->dpdf[iat + 1][jat + 1];
      fi[4] = helm_eos_table->dpdft[iat][jat];
      fi[5] = helm_eos_table->dpdft[iat + 1][jat];
      fi[6] = helm_eos_table->dpdft[iat][jat + 1];
      fi[7] = helm_eos_table->dpdft[iat + 1][jat + 1];
      fi[8] = helm_eos_table->dpdfd[iat][jat];
      fi[9] = helm_eos_table->dpdfd[iat + 1][jat];
      fi[10] = helm_eos_table->dpdfd[iat][jat + 1];
      fi[11] = helm_eos_table->dpdfd[iat + 1][jat + 1];
      fi[12] = helm_eos_table->dpdfdt[iat][jat];
      fi[13] = helm_eos_table->dpdfdt[iat + 1][jat];
      fi[14] = helm_eos_table->dpdfdt[iat][jat + 1];
      fi[15] = helm_eos_table->dpdfdt[iat + 1][jat + 1];

      //pressure derivative with density
      pele[1] = dmax(0.0, ye * h3(fi, si0t, si1t, si0mt, si1mt, si0d, si1d, si0md, si1md));



      //look in the electron chemical potential table only once
      fi[0] = helm_eos_table->ef[iat][jat];
      fi[1] = helm_eos_table->ef[iat + 1][jat];
      fi[2] = helm_eos_table->ef[iat][jat + 1];
      fi[3] = helm_eos_table->ef[iat + 1][jat + 1];
      fi[4] = helm_eos_table->eft[iat][jat];
      fi[5] = helm_eos_table->eft[iat + 1][jat];
      fi[6] = helm_eos_table->eft[iat][jat + 1];
      fi[7] = helm_eos_table->eft[iat + 1][jat + 1];
      fi[8] = helm_eos_table->efd[iat][jat];
      fi[9] = helm_eos_table->efd[iat + 1][jat];
      fi[10] = helm_eos_table->efd[iat][jat + 1];
      fi[11] = helm_eos_table->efd[iat + 1][jat + 1];
      fi[12] = helm_eos_table->efdt[iat][jat];
      fi[13] = helm_eos_table->efdt[iat + 1][jat];
      fi[14] = helm_eos_table->efdt[iat][jat + 1];
      fi[15] = helm_eos_table->efdt[iat + 1][jat + 1];


      //electron chemical potential etaele
      etaele[0] = h3(fi, si0t, si1t, si0mt, si1mt, si0d, si1d, si0md, si1md);

      //derivative with respect to density
      x = h3(fi, si0t, si1t, si0mt, si1mt, dsi0d, dsi1d, dsi0md, dsi1md);
      etaele[1] = ye * x;

      //derivative with respect to temperature
      etaele[2] = h3(fi, dsi0t, dsi1t, dsi0mt, dsi1mt, si0d, si1d, si0md, si1md);

      //derivative with respect to abar and zbar
      etaele[3] = -x * din * ytot;
      etaele[4] = x * rho * ytot;



      //look in the number density table only once
      fi[0] = helm_eos_table->xf[iat][jat];
      fi[1] = helm_eos_table->xf[iat + 1][jat];
      fi[2] = helm_eos_table->xf[iat][jat + 1];
      fi[3] = helm_eos_table->xf[iat + 1][jat + 1];
      fi[4] = helm_eos_table->xft[iat][jat];
      fi[5] = helm_eos_table->xft[iat + 1][jat];
      fi[6] = helm_eos_table->xft[iat][jat + 1];
      fi[7] = helm_eos_table->xft[iat + 1][jat + 1];
      fi[8] = helm_eos_table->xfd[iat][jat];
      fi[9] = helm_eos_table->xfd[iat + 1][jat];
      fi[10] = helm_eos_table->xfd[iat][jat + 1];
      fi[11] = helm_eos_table->xfd[iat + 1][jat + 1];
      fi[12] = helm_eos_table->xfdt[iat][jat];
      fi[13] = helm_eos_table->xfdt[iat + 1][jat];
      fi[14] = helm_eos_table->xfdt[iat][jat + 1];
      fi[15] = helm_eos_table->xfdt[iat + 1][jat + 1];

      //electron + positron number densities
      xne[0] = h3(fi, si0t, si1t, si0mt, si1mt, si0d, si1d, si0md, si1md);

      //derivative with respect to density
      x = dmax(0.0, h3(fi, si0t, si1t, si0mt, si1mt, dsi0d, dsi1d, dsi0md, dsi1md));
      xne[1] = ye * x;

      //derivative with respect to temperature
      xne[2] = h3(fi, dsi0t, dsi1t, dsi0mt, dsi1mt, si0d, si1d, si0md, si1md);

      //derivative with respect to abar and zbar
      xne[3] = -x * din * ytot;
      xne[4] = x * rho * ytot;



      //the desired electron-positron thermodynamic quantities

      //dpepdd at high temperatures and low densities is below the
      //floating point limit of the subtraction of two large terms.
      //since dpresdd doesn't enter the maxwell relations at all, use the
      //bicubic interpolation done above instead of this one
      x = din * din;
      pele[0] = x * df_d;
      pele[2] = x * df_dt;
      //pele[1]  = ye * (x * df_dd + 2.0 * din * df_d);
      s = pele[1] / ye - 2.0 * din * df_d;
      pele[3] = -ytot * (2.0 * pele[0] + s * din);
      pele[4] = rho * ytot * (2.0 * din * df_d + s);

      x = ye * ye;
      sele[0] = -df_t * ye;
      sele[2] = -df_tt * ye;
      sele[1] = -df_dt * x;
      sele[3] = ytot * (ye * df_dt * din - sele[0]);
      sele[4] = -ytot * (ye * df_dt * rho + df_t);


      eele[0] = ye * free_en + temp * sele[0];
      eele[2] = temp * sele[2];
      eele[1] = x * df_d + temp * sele[1];
      eele[3] = -ye * ytot * (free_en + df_d * din) + temp * sele[3];
      eele[4] = ytot * (free_en + ye * df_d * rho) + temp * sele[4];
    }
  else
    {
      eele[0] = ye * free_en + temp * (-df_t) * ye;
      eele[2] = temp * (-df_tt) * ye;
    }

  return 0;
}

#ifdef EOS_COULOMB_CORRECTIONS_SMOOTH
// butterworth low-pass and first derivative by ^_^ : beware the butterbomb
struct Filter{
  double g; // gain
  double dgdf; // d(gain)/d(frequency)
};

//struct Filter butterworth();
static void butterworth(double freq, double cfreq, int n, struct Filter *result) {
  //    butterworth filter input: frequency, central frequency and order of
  //    filter (freq, cfreq, n), returns: gain and d(gain)/d(frequnency)

  result->g = 1.0 + gsl_pow_int((freq / cfreq),2 * n);
  result->g = 1.0 / result->g;
  result->dgdf = -gsl_pow_2(result->g) * 2 * n / gsl_pow_int(cfreq,2 * n) * gsl_pow_int(freq,(2 * n - 1));
}
#endif

// Coulomb corrections
static int __attribute__ ((unused)) helm_eos_coul(double rho, double temp, double abar, double zbar, double pcoul[5], double ecoul[5], double scoul[5])
{
  // fitting parameters
  const double a1 = -0.898004, b1 = 0.96786, c1 = 0.220703, d1 = -0.86097, e1 = 2.5269, a2 = 0.29561, b2 = 1.9885, c2 = 0.288675;
  //const double ye = zbar / abar; // electron number fraction
  const double ytot = 1.0 / abar;
  const double kt = GSL_CONST_CGS_BOLTZMANN * temp;
  double xni, dxnidd, dxnida;
  double pion, dpiondd, dpiondt, dpionda, dpiondz;
  double s, dsdd, dsda;
  double lami, lamidd, lamida;  // average ion serperation
  double plasg, plasgdd, plasgda, plasgdt, plasgdz;     // plasma coupling parameter

  xni = GSL_CONST_NUM_AVOGADRO * ytot * rho;
  dxnidd = GSL_CONST_NUM_AVOGADRO * ytot;
  dxnida = -xni * ytot;

  pion = xni * kt;
  dpiondd = dxnidd * kt;
  dpiondt = xni * GSL_CONST_CGS_BOLTZMANN;
  dpionda = dxnida * kt;
  dpiondz = 0.0;

  s = 4.0 / 3.0 * M_PI * xni;
  dsdd = 4.0 / 3.0 * M_PI * dxnidd;
  dsda = 4.0 / 3.0 * M_PI * dxnida;

  lami = 1.0 / cbrt(s);
  lamidd = -lami / 3.0 * dsdd / s;
  lamida = -lami / 3.0 * dsda / s;

  plasg = gsl_pow_2(ELECTRON_CHARGE_ESU * zbar) / (kt * lami);
  plasgdd = -plasg / lami * lamidd;
  plasgda = -plasg / lami * lamida;
  plasgdt = -plasg / temp;
  plasgdz = 2.0 * plasg / zbar;

  if(plasg >= 1.0)
    {
      double x = sqrt(sqrt(plasg)), y = GSL_CONST_NUM_AVOGADRO * GSL_CONST_CGS_BOLTZMANN * ytot;
      ecoul[0] = y * temp * (a1 * plasg + b1 * x + c1 / x + d1);
      pcoul[0] = rho * ecoul[0] / 3.0;
      scoul[0] = -y * (3.0 * b1 * x - 5.0 * c1 / x + d1 * (log(plasg) - 1.0) - e1);

      y = GSL_CONST_NUM_AVOGADRO * GSL_CONST_CGS_BOLTZMANN * temp * ytot * (a1 + 0.25 / plasg * (b1 * x - c1 / x));
      ecoul[1] = y * plasgdd;
      ecoul[2] = y * plasgdt + ecoul[0] / temp;
      ecoul[3] = y * plasgda - ecoul[0] / abar;
      ecoul[4] = y * plasgdz;

      y = rho / 3.0;
      pcoul[1] = ecoul[0] / 3.0 + y * ecoul[1];
      pcoul[2] = y * ecoul[2];
      pcoul[3] = y * ecoul[3];
      pcoul[4] = y * ecoul[4];

      y = -GSL_CONST_NUM_AVOGADRO * GSL_CONST_CGS_BOLTZMANN / (abar * plasg) * (0.75 * b1 * x + 1.25 * c1 / x + d1);
      scoul[1] = y * plasgdd;
      scoul[2] = y * plasgdt;
      scoul[3] = y * plasgda - scoul[0] / abar;
      scoul[4] = y * plasgdz;
    }
  else
    {
      const double x = plasg * sqrt(plasg), y = pow(plasg, b2), z = c2 * x - a2 / 3.0 * y;

      pcoul[0] = -pion * z;
      ecoul[0] = 3.0 * pcoul[0] / rho;
      scoul[0] = -GSL_CONST_NUM_AVOGADRO * GSL_CONST_CGS_BOLTZMANN / abar * (c2 * x - a2 * (b2 - 1.0) / b2 * y);

      s = 1.5 * c2 * x / plasg - a2 * b2 / 3.0 * y / plasg;
      pcoul[1] = -dpiondd * z - pion * s * plasgdd;
      pcoul[2] = -dpiondt * z - pion * s * plasgdt;
      pcoul[3] = -dpionda * z - pion * s * plasgda;
      pcoul[4] = -dpiondz * z - pion * s * plasgdz;

      s = 3.0 / rho;
      ecoul[1] = s * pcoul[1] - ecoul[0] / rho;
      ecoul[2] = s * pcoul[2];
      ecoul[3] = s * pcoul[3];
      ecoul[4] = s * pcoul[4];

      s = -GSL_CONST_NUM_AVOGADRO * GSL_CONST_CGS_BOLTZMANN / (abar * plasg) * (1.5 * c2 * x - a2 * (b2 - 1.0) * y);
      scoul[1] = s * plasgdd;
      scoul[2] = s * plasgdt;
      scoul[3] = s * plasgda - scoul[0] / abar;
      scoul[4] = s * plasgdz;
    }

#ifdef EOS_COULOMB_CORRECTIONS_SMOOTH
  // butterworth bomb proofing by ^_^ : beware the butterbomb
  struct Filter tfilter, dfilter;
  double gain, dgaindt, dgaindd;

  butterworth(log10(temp) - 4.5, 3.0, 12, &tfilter);
  butterworth(log10(rho) + 1.0, 6.0, 12, &dfilter);

  // derivatives (and conversion from logarithmic derivative)
  gain = (1.0 - tfilter.g * dfilter.g);
  dgaindt = -dfilter.g * tfilter.dgdf / temp / log(10.);
  dgaindd = -tfilter.g * dfilter.dgdf / rho / log(10.);

  // straight up gain
  pcoul[0] = pcoul[0] * gain;
  ecoul[0] = ecoul[0] * gain;
  scoul[0] = scoul[0] * gain;

  // derivatives via chain rule
  pcoul[1] = gain * pcoul[1] + pcoul[0] * dgaindd;
  pcoul[2] = gain * pcoul[2] + pcoul[0] * dgaindt;
  pcoul[3] = gain * pcoul[3];
  pcoul[4] = gain * pcoul[4];

  ecoul[1] = gain * ecoul[1] + ecoul[0] * dgaindd;
  ecoul[2] = gain * ecoul[2] + ecoul[0] * dgaindt;
  ecoul[3] = gain * ecoul[3];
  ecoul[4] = gain * ecoul[4];

  scoul[1] = gain * scoul[1] + scoul[0] * dgaindd;
  scoul[2] = gain * scoul[2] + scoul[0] * dgaindt;
  scoul[3] = gain * scoul[3];
  scoul[4] = gain * scoul[4];
#endif
  return 0;
}


int eos_init(const char *datafile, const char *speciesfile)
{
  FILE *file;

  file = fopen(datafile, "r");
  if(file == NULL)
    {
      perror("error opening EOS table file");
      fprintf(stderr, "the filname was `%s'\n", datafile);
      return -1;
    }

  helm_eos_table = malloc(sizeof(struct helm_eos_table));
  if(helm_eos_table == NULL)
    {
      perror("could not allocate memory for the EOS table");
      fclose(file);
      return -1;
    }

  helm_eos_table->ntemp = JMAX;
  helm_eos_table->nrho = IMAX;

  helm_eos_table->ltempMin = 3.0;
  helm_eos_table->ltempMax = 13.0;
  helm_eos_table->lrhoMin = -12.0;
  helm_eos_table->lrhoMax = 15.0;

  helm_eos_table->tempMin = exp10(helm_eos_table->ltempMin);
  helm_eos_table->tempMax = exp10(helm_eos_table->ltempMax);
  helm_eos_table->rhoMin = exp10(helm_eos_table->lrhoMin);
  helm_eos_table->rhoMax = exp10(helm_eos_table->lrhoMax);

  helm_eos_table->ltempDelta = (helm_eos_table->ltempMax - helm_eos_table->ltempMin) / (double) (helm_eos_table->ntemp - 1);
  helm_eos_table->lrhoDelta = (helm_eos_table->lrhoMax - helm_eos_table->lrhoMin) / (double) (helm_eos_table->nrho - 1);


  // read the Helmholtz free energy table and its derivatives
  {
    int i, j;
    for(j = 0; j < helm_eos_table->ntemp; j++)
      {
        helm_eos_table->temp[j] = exp10(helm_eos_table->ltempMin + j * helm_eos_table->ltempDelta);
        for(i = 0; i < helm_eos_table->nrho; i++)
          {
            helm_eos_table->rho[i] = exp10(helm_eos_table->lrhoMin + i * helm_eos_table->lrhoDelta);
            if(fscanf
               (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &helm_eos_table->f[i][j],
                &helm_eos_table->fd[i][j], &helm_eos_table->ft[i][j], &helm_eos_table->fdd[i][j],
                &helm_eos_table->ftt[i][j], &helm_eos_table->fdt[i][j], &helm_eos_table->fddt[i][j], &helm_eos_table->fdtt[i][j], &helm_eos_table->fddtt[i][j]) != 9)
              {
                fprintf(stderr, "error reading the Helmholtz free energy table at i = %d, j = %d\n", i, j);
                fclose(file);
                free(helm_eos_table);
                return -1;
              }
          }
      }
  }

  // read the pressure derivative table
  {
    int i, j;
    for(j = 0; j < helm_eos_table->ntemp; j++)
      {
        for(i = 0; i < helm_eos_table->nrho; i++)
          {
            if(fscanf(file, "%lf %lf %lf %lf", &helm_eos_table->dpdf[i][j], &helm_eos_table->dpdfd[i][j], &helm_eos_table->dpdft[i][j], &helm_eos_table->dpdfdt[i][j]) != 4)
              {
                fprintf(stderr, "error reading the pressure derivative table at i = %d, j = %d\n", i, j);
                fclose(file);
                free(helm_eos_table);
                return -1;
              }
          }
      }
  }

  // read the electron chemical potential table
  {
    int i, j;
    for(j = 0; j < helm_eos_table->ntemp; j++)
      {
        for(i = 0; i < helm_eos_table->nrho; i++)
          {
            if(fscanf(file, "%lf %lf %lf %lf", &helm_eos_table->ef[i][j], &helm_eos_table->efd[i][j], &helm_eos_table->eft[i][j], &helm_eos_table->efdt[i][j]) != 4)
              {
                fprintf(stderr, "error reading the electron chemical potential table at i = %d, j = %d\n", i, j);
                fclose(file);
                free(helm_eos_table);
                return -1;
              }
          }
      }
  }

  // read the number density table
  {
    int i, j;
    for(j = 0; j < helm_eos_table->ntemp; j++)
      {
        for(i = 0; i < helm_eos_table->nrho; i++)
          {
            if(fscanf(file, "%lf %lf %lf %lf", &helm_eos_table->xf[i][j], &helm_eos_table->xfd[i][j], &helm_eos_table->xft[i][j], &helm_eos_table->xfdt[i][j]) != 4)
              {
                fprintf(stderr, "error reading the number density table at i = %d, j = %d\n", i, j);
                fclose(file);
                free(helm_eos_table);
                return -1;
              }
          }
      }
  }

  fclose(file);

  // read the species file
  file = fopen(speciesfile, "r");
  if(file == NULL)
    {
      perror("cannot open species file");
      fprintf(stderr, "The filename was `%s'.\n", speciesfile);
      free(helm_eos_table);
      fclose(file);
      return -1;
    }

  if(fscanf(file, "%d", &helm_eos_table->nspecies) != 1)
    {
      fprintf(stderr, "error in species file format");
      free(helm_eos_table);
      fclose(file);
      return -1;
    }

  helm_eos_table->na = malloc(helm_eos_table->nspecies * sizeof(double));
  helm_eos_table->nai = malloc(helm_eos_table->nspecies * sizeof(double));
  helm_eos_table->nz = malloc(helm_eos_table->nspecies * sizeof(double));

  if(helm_eos_table->na == NULL || helm_eos_table->nai == NULL || helm_eos_table->nz == NULL)
    {
      perror("error allocating species table");
      if(!helm_eos_table->na)
        free(helm_eos_table->na);
      if(!helm_eos_table->nai)
        free(helm_eos_table->na);
      if(!helm_eos_table->nz)
        free(helm_eos_table->nz);
      free(helm_eos_table);
      fclose(file);
      return -1;
    }

  {
    int i;
    for(i = 0; i < helm_eos_table->nspecies; i++)
      {
        if(fscanf(file, "%*5s %lf %lf", helm_eos_table->na + i, helm_eos_table->nz + i) != 2)
          {
            fprintf(stderr, "error in species file at element %d\n", i);
            free(helm_eos_table->na);
            free(helm_eos_table->nai);
            free(helm_eos_table->nz);
            free(helm_eos_table);
            fclose(file);
            return -1;
          }
        helm_eos_table->nai[i] = 1.0 / helm_eos_table->na[i];
      }
  }

  fclose(file);
  return 0;
}

void eos_deinit(void)
{
  free(helm_eos_table->na);
  free(helm_eos_table->nai);
  free(helm_eos_table->nz);
  free(helm_eos_table);
}

int eos_calc_tgiven(double rho, const double xnuc[], double temp, struct eos_result *res)
{
  struct helm_eos_cache cache;

  azbar(xnuc, &cache.abar, &cache.zbar);
  helm_eos_update_cache(rho, cache.abar, cache.zbar, &cache);

  return eos_calc_tgiven_azbar(rho, &cache, temp, res, 0);
}

int eos_calc_tgiven_onlye(double rho, const double xnuc[], double temp, struct eos_result *res)
{
  struct helm_eos_cache cache;

  azbar(xnuc, &cache.abar, &cache.zbar);
  helm_eos_update_cache(rho, cache.abar, cache.zbar, &cache);

  return eos_calc_tgiven_azbar(rho, &cache, temp, res, 1);
}

int eos_calc_tgiven_azbar(double rho, const struct helm_eos_cache *cache, double temp, struct eos_result *res, const int only_e)
{
  // arrays for the different contributions and their derivatives
  // all are initialized to zero
  double prad[5] = { 0 }, pion[5] =
  {
  0}, pele[5] =
  {
  0}, pcoul[5] =
  {
  0};
  double erad[5] = { 0 }, eion[5] =
  {
  0}, eele[5] =
  {
  0}, ecoul[5] =
  {
  0};
  double srad[5] = { 0 }, sion[5] =
  {
  0}, sele[5] =
  {
  0}, scoul[5] =
  {
  0};
  // electron chemical potential and electron + positron number density
  double etaele[5] = { 0 }, xne[5] =
  {
  0};
  const size_t offsets[5] = { offsetof(struct eos_value, v), offsetof(struct eos_value, drho), offsetof(struct eos_value, dtemp),
    offsetof(struct eos_value, dabar), offsetof(struct eos_value, dzbar)
  };
  const double abar = cache->abar, zbar = cache->zbar;
  double ltemp = log(temp);
  double *data;

  if(helm_eos_ion(rho, temp, ltemp, cache, abar, zbar, pion, eion, sion) != 0)
    return -1;
#ifndef EOS_IDEAL
  if(helm_eos_rad(rho, temp, abar, zbar, prad, erad, srad) != 0)
    return -1;
  if(helm_eos_ele(rho, temp, ltemp * (1.0 / log(10)), pele, eele, sele, etaele, xne, cache, only_e) != 0)
    return -1;
#ifdef EOS_COULOMB_CORRECTIONS
  if(helm_eos_coul(rho, temp, abar, zbar, pcoul, ecoul, scoul) != 0)
    return -1;
#ifndef EOS_COULOMB_CORRECTIONS_SMOOTH
  //If Coulomb corrections result in P or E falling below zero, reset Coulomb corrections to zero.
  if(prad[0] + pion[0] + pele[0] + pcoul[0] <= 0.0 && only_e == 0 || erad[0] + eion[0] + eele[0] + ecoul[0] <= 0.0)
    {
      pcoul[0] = 0.0;
      pcoul[1] = 0.0;
      pcoul[2] = 0.0;
      pcoul[3] = 0.0;
      pcoul[4] = 0.0;
      ecoul[0] = 0.0;
      ecoul[1] = 0.0;
      ecoul[2] = 0.0;
      ecoul[3] = 0.0;
      ecoul[4] = 0.0;
      scoul[0] = 0.0;
      scoul[1] = 0.0;
      scoul[2] = 0.0;
      scoul[3] = 0.0;
      scoul[4] = 0.0;
    }
#endif // ! EOS_COULOMB_CORRECTIONS_SMOOTH
#endif // EOS_COULOMB_CORRECTIONS
#endif // ! EOS_IDEAL

  if(res != NULL)
    res->temp = temp;

  if(res != NULL && !only_e)
    {
      double x;

      {
        int i;
        for(i = 0; i < 5; i++)
          {
            *(double *) ((char *) &res->p + offsets[i]) = prad[i] + pion[i] + pele[i] + pcoul[i];
            *(double *) ((char *) &res->e + offsets[i]) = erad[i] + eion[i] + eele[i] + ecoul[i];
            *(double *) ((char *) &res->s + offsets[i]) = srad[i] + sion[i] + sele[i] + scoul[i];
            *(double *) ((char *) &res->etaele + offsets[i]) = etaele[i];
            *(double *) ((char *) &res->nep + offsets[i]) = xne[i];
          }
      }

      res->cv = res->e.dtemp;
      res->chit = temp / res->p.v * res->p.dtemp;
      res->chid = res->p.drho * rho / res->p.dtemp;
      x = res->p.v / rho * res->chit / (temp * res->cv);
      res->gamma_3 = x + 1.0;
      res->gamma_1 = res->chit * x + res->chid;
      res->nabla_ad = x / res->gamma_1;
      res->gamma_2 = 1.0 / (1.0 - res->nabla_ad);
      res->cp = res->cv * res->gamma_1 / res->chid;
      res->sound = GSL_CONST_CGS_SPEED_OF_LIGHT * sqrt(res->gamma_1 / (1.0 + (res->e.v + gsl_pow_2(GSL_CONST_CGS_SPEED_OF_LIGHT)) * rho / res->p.v));

      res->abar = abar;
      res->zbar = zbar;
    }

  if(res != NULL && only_e)
    {
      {
        int i;
        for(i = 0; i < 5; i++)
          {
            data = &res->e.v;
            data[i] = erad[i] + eion[i] + eele[i] + ecoul[i];
          }
      }
    }

  return 0;
}

int eos_calc_egiven(double rho, const double xnuc[], double e, double *tempguess, struct eos_result *res)
{
  int iter;                     // number of Newton-Raphson iterations
  double ni, ne, nn;            // number densities
  double _temp, _tempold;       // temperature, pressure, and energy during the iteration
  _Bool have_to_free = false;
  struct helm_eos_cache cache;

  azbar(xnuc, &cache.abar, &cache.zbar);
  helm_eos_update_cache(rho, cache.abar, cache.zbar, &cache);

  ni = 1.0 / cache.abar * rho * GSL_CONST_NUM_AVOGADRO;
  nn = rho * GSL_CONST_NUM_AVOGADRO;
  ne = cache.zbar * ni;

  // use negative temperatures to show that we need to make a guess for the temperature
  if(tempguess == 0 || *tempguess <= 0.0)
    {
      // guess using an ideal gas
      _temp = 2.0 / 3.0 * e * rho / (ni + ne) / GSL_CONST_CGS_BOLTZMANN;
    }
  else
    {
      _temp = *tempguess;
    }

#ifndef EOS_IDEAL
  _temp = dmax(dmin(_temp, helm_eos_table->tempMax), helm_eos_table->tempMin);
#endif // ! EOS_IDEAL

  if(res == NULL)
    {
      res = malloc(sizeof(struct eos_result));
      have_to_free = true;
    }

  _tempold = 0.0;
  for(iter = 0; iter < HELM_EOS_MAXITER; iter++)
    {
      /* update only energy on res */
      if(eos_calc_tgiven_azbar(rho, &cache, _temp, res, 1) != 0)
        return -1;

      // check if we are converged already
      if(fabs(res->e.v - e) <= (HELM_EOS_EPS * e))
        break;

      // not converged; compute the next step
      _tempold = _temp;
      _temp = _temp - (res->e.v - e) / res->e.dtemp;

#ifndef EOS_IDEAL
      // check if we left the table
      if(_temp <= helm_eos_table->tempMin && _tempold <= helm_eos_table->tempMin)
        {
          _temp = helm_eos_table->tempMin;
          break;
        }
      if(_temp >= helm_eos_table->tempMax && _tempold >= helm_eos_table->tempMax)
        {
          _temp = helm_eos_table->tempMax;
          break;
        }

      _temp = dmax(dmin(_temp, helm_eos_table->tempMax), helm_eos_table->tempMin);
#endif // ! EOS_IDEAL
    }

  if(iter >= HELM_EOS_MAXITER)
    {
      fprintf(stderr, "Newton-Raphson in function `%s' did not converge.\n", __func__);

      if(tempguess && *tempguess > 0.)
        {
          if(have_to_free)
            free(res);
          return -2;
        }
    }

  if(tempguess)
    *tempguess = _temp;

  /* update all quantities on res */
  if(eos_calc_tgiven_azbar(rho, &cache, _temp, res, 0) != 0)
    {
      if(have_to_free)
        free(res);
      return -1;
    }

  /* FIXME: dpdr should be checked again */
  /* *dpdr = res->p.drho + _temp * gsl_pow_2(res->p.dtemp / rho) / res->e.dtemp;
     assert(*dpdr > 0); */

  if(have_to_free)
    {
      free(res);
    }

  return 0;
}

int eos_calc_egiven_y(double rho, const double y[], double e, double *tempguess, struct eos_result *res)
{
  const double *na = helm_eos_table->na;
  int nspecies = helm_eos_table->nspecies;
  double *xnuc;
  int result;

  xnuc = malloc(nspecies * sizeof(double));

  {
    int i;
    for(i = 0; i < nspecies; i++)
      xnuc[i] = y[i] * na[i];
  }
  result = eos_calc_egiven(rho, xnuc, e, tempguess, res);

  free(xnuc);
  return result;
}

int eos_calc_ptgiven(double p, const double xnuc[], double temp, double *rho, struct eos_result *res)
{
  int iter;                     // number of Newton-Raphson iterations
  double abar, zbar;            // mean nuclear mass and charge
  double _rho, _rhoold;         // temperature and pressure during the iteration
  double dpdd;                  // derivative of pressure w. r. t. density
  _Bool have_to_free = false;

  azbar(xnuc, &abar, &zbar);

  // use negative densities to show that we need to make a guess for the density
  if(*rho < 0.0)
    {
      // guess using an ideal gas
      _rho = abar / (1.0 + zbar) * p / (GSL_CONST_NUM_AVOGADRO * GSL_CONST_CGS_BOLTZMANN * temp);
    }
  else
    {
      _rho = *rho;
    }

  if(res == NULL)
    {
      res = malloc(sizeof(struct eos_result));
      have_to_free = true;
    }

  _rhoold = 0.0;
  for(iter = 0; iter < HELM_EOS_MAXITER; iter++)
    {
      if(eos_calc_tgiven(_rho, xnuc, temp, res) != 0)
        return -1;

      // check if we are converged already
      if(fabs(res->p.v - p) <= (HELM_EOS_EPS * p))
        break;

      // not converged; compute the next step
      _rhoold = _rho;
      dpdd = res->p.drho;
      _rho = _rho - (res->p.v - p) / dpdd;

      // check if we left the table
      if(_rho <= helm_eos_table->rhoMin && _rhoold <= helm_eos_table->rhoMin)
        {
          _rho = helm_eos_table->rhoMin;
          break;
        }
      if(_rho >= helm_eos_table->rhoMax && _rhoold >= helm_eos_table->rhoMax)
        {
          _rho = helm_eos_table->rhoMax;
          break;
        }

      _rho = dmax(dmin(_rho, helm_eos_table->rhoMax), helm_eos_table->rhoMin);
    }

  if(iter >= HELM_EOS_MAXITER)
    {
      fprintf(stderr, "Newton-Raphson in function `%s' did not converge.\n", __func__);
      return -1;
    }
#ifndef NDEBUG
  printf("needed %d iterations\n", iter);
#endif

  if(eos_calc_tgiven(_rho, xnuc, temp, res) != 0)
    return -1;

  *rho = _rho;

  if(have_to_free)
    {
      free(res);
    }

  return 0;
}

int helm_eos_update_cache(double rho, double abar, double zbar, struct helm_eos_cache *cache)
{
  cache->abar = abar;
  cache->zbar = zbar;
  cache->ytot = 1.0 / cache->abar;
  cache->ye = zbar * cache->ytot;
  cache->din = rho * cache->ye;
  cache->ldin = log10(cache->din);
  cache->xni = GSL_CONST_NUM_AVOGADRO * cache->ytot * rho;
  cache->dxnidd = GSL_CONST_NUM_AVOGADRO * cache->ytot;
  cache->dxnida = -cache->xni * cache->ytot;
  cache->lswot15 = 1.5 * log((2.0 * M_PI * GSL_CONST_CGS_UNIFIED_ATOMIC_MASS * GSL_CONST_CGS_BOLTZMANN) / (GSL_CONST_CGS_PLANCKS_CONSTANT_H * GSL_CONST_CGS_PLANCKS_CONSTANT_H));
  cache->ywot = log(cache->abar * cache->abar * sqrt(cache->abar) / (rho * GSL_CONST_NUM_AVOGADRO));

  return 0;
}

#endif /* EOS_DEGENERATE */
