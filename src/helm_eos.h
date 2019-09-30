/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/helm_eos.h
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

#ifndef HELM_EOS_H
#define HELM_EOS_H

#include "arepoconfig.h"

#ifdef EOS_DEGENERATE

struct eos_value
{
  double v;                     /* value */
  double drho;                  /* derivative with density */
  double dtemp;                 /* derivative with temperature */
  double dabar;                 /* derivative with abar */
  double dzbar;                 /* derivative with zbar */
};

struct eos_result
{
  double temp;
  struct eos_value p;           /*  pressure */
  struct eos_value e;           /*  specific energy */
  struct eos_value s;           /*  specific entropy */
  struct eos_value etaele;      /*  degeneracy parameter: electron chemical potential / (k_B * T) */
  struct eos_value nep;         /*  electron + positron number density */
  double cv, cp;                /*  specific heat at constant volume, pressure */
  double chit, chid;            /*  temperature and density exponents from Cox & Giuli */
  double gamma_1, gamma_2, gamma_3;     /*  gammas from Cox & Giuli */
  double nabla_ad;              /*  nabla adiabatic */
  double sound;                 /*  relativistic speed of sound */

  double abar, zbar;            /*  mean mass number and charge number */
};

/* cached values that only depend on composition and density
 * huge time saver for the egiven iterations
 */
struct helm_eos_cache
{
  double abar, zbar;            /* mean nuclear mass and charge */
  double ye;                    /* zbar / abar */
  double din, ldin;             /* rho * ye */
  double ytot;                  /* 1.0 / abar; */
  double xni;                   /* GSL_CONST_NUM_AVOGADRO * ytot * rho */
  double dxnidd;                /* GSL_CONST_NUM_AVOGADRO * ytot */
  double dxnida;                /* - xni * ytot */
  double lswot15;               /* 1.5 * log((2.0 * M_PI * GSL_CONST_CGS_UNIFIED_ATOMIC_MASS * GSL_CONST_CGS_BOLTZMANN) / (GSL_CONST_CGS_PLANCKS_CONSTANT_H * GSL_CONST_CGS_PLANCKS_CONSTANT_H) * temp) without the temp factor */
  double ywot;                  /* log(abar * abar * sqrt(abar) / (rho * GSL_CONST_NUM_AVOGADRO)) */
};


int eos_init(const char *datafile, const char *speciesfile);
void eos_deinit(void);

int eos_calc_tgiven(double rho, const double xnuc[], double temp, struct eos_result *res);
int eos_calc_tgiven_onlye(double rho, const double xnuc[], double temp, struct eos_result *res);
int eos_calc_tgiven_azbar(double rho, const struct helm_eos_cache *cache, double temp, struct eos_result *res, const int only_e);
int eos_calc_egiven(double rho, const double xnuc[], double e, double *tempguess, struct eos_result *res);
int eos_calc_egiven_y(double rho, const double y[], double e, double *tempguess, struct eos_result *res);
int eos_calc_ptgiven(double p, const double xnuc[], double temp, double *rho, struct eos_result *res);

int helm_eos_update_cache(double rho, double abar, double zbar, struct helm_eos_cache *cache);

#endif /* EOS_DEGENERATE */

#endif /* HELM_EOS_H */
