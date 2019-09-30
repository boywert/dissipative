/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/opal_eos.h
 * \date        02/2014
 * \author      Sebastian Ohlmann <sohlmann@astro.uni-wuerzburg.de>
 * \brief       header file for OPAL equation of state
 * \details     implementation of the OPAL equation of state in C using
 *    interpolation of the freely available data
 *
 *
 * \par Major modifications and contributions:
 *
 * - 10.02.2014 Integration in Arepo
 * - 21.04.2016 extension to larger area; blending of ideal gas to OPAL
 */
#ifndef OPAL_EOS_H
#define OPAL_EOS_H

#include "arepoconfig.h"

#ifdef EOS_OPAL

/* #define OPAL_DEBUG */
#undef OPAL_DEBUG

/* physical constants */
//#define SIGMA 5.670400e-5  /* erg /(s cm**2 K**4); 2010 CODATA value */
#define SIGMA 5.67051e-5        /* erg /(s cm**2 K**4); old value */
//#define SPEED_OF_LIGHT 29979245800 /* cm / s; exact */
#define SPEED_OF_LIGHT 29979245000      /* cm / s; old value */
#define SIGMAC SIGMA / SPEED_OF_LIGHT
#define BOLTZMANN_K 1.3806504e-16       /* erg/K; 2010 CODATA value */
#define UAMU 1.660538782e-24    /* g; 2010 CODATA value */
//#define GAS_CONSTANT 8.31447248e7 /* = k/u; unit: erg/(K g) */
#define GAS_CONSTANT 8.314510e7 /* = k/u; unit: erg/(K g) old value */
#define HBAR 1.054571628e-27 /* unit: erg s; 2010 CODATA value */
#define PI 3.141592653589793238
#define MASS_E 9.10938215e-28 /* unit: g; 2010 CODATA value */

/* definitions for eos results */
#define EOS_ENERGY   (1 << 0)
#define EOS_PRESSURE (1 << 1)
#define EOS_ENTROPY  (1 << 2)
#define EOS_DEDRHO   (1 << 3)
#define EOS_CV       (1 << 4)
#define EOS_CHIR     (1 << 5)
#define EOS_CHIT     (1 << 6)
#define EOS_GAMMA1   (1 << 7)
#define EOS_GAMMA2   (1 << 8)
#define EOS_GAMMA3   (1 << 9)
#define EOS_RADIATION (1 << 10)
#define EOS_ALL (EOS_ENERGY | EOS_PRESSURE | EOS_ENTROPY | EOS_DEDRHO |   \
    EOS_CV | EOS_CHIR | EOS_CHIT | EOS_GAMMA1 | EOS_GAMMA2 | EOS_GAMMA3 | \
    EOS_RADIATION)

/* eos table definitions */
#define NRHO 169
#define NT 197
#define NX 5

#define IP_LINEAR 1
#define IP_QUADRATIC 2
#define INTERPOLATION_METHOD IP_QUADRATIC
//#define INTERPOLATION_METHOD IP_LINEAR

/* global table object? */
#define OPAL_GLOBAL

/* maximum number of Newton iterations */
#define OPAL_NEWTON_MAX 70
#define OPAL_NEWTON_ACCURACY 1e-10

typedef double eostable[NT][NRHO][NX]; /**< table for eos variables */

/** @brief EOS result struct.
 * It contains variables for all EOS variables.
 */
struct opal_eos_result
{
  /* eos variables */
  double p;      /**< pressure in dyn/cm^2                             */
  double e;      /**< energy in erg/g                                  */
  double s;      /**< entropy in erg/(g K)                             */
  double dedrho; /**< dE/dRHO at constant T6                           */
  double cv;     /**< specific heat, dE/dT6 at constant V in erg/(K g) */
  double chir;   /**< dlogP/dlogRho at constant T6.
                  *   Cox and Guil1 eq 9.82                            */
  double chit;   /**< dlogP/dlogT6 at constant Rho.
                  *   Cox and Guil1 eq 9.81                            */
  double gamma1; /**< gamma1. Eqs. 9.88 Cox and Guili.                 */
  double gamma2; /**< gamma2, Eqs. 9.88 Cox and Guili                  */
  double gamma3; /**< gamma3 is computed from gamma1 and gamma2        */
};

/**
 * @brief cache structure for faster computation
 */
struct opal_eos_cache
{
  double x, rho, t;                  /**< values for evaluation */
  int ix, irho, it;                  /**< corresponding low indices */
  int mx, mrho, mt;                  /**< indices directly near values */
#if INTERPOLATION_METHOD == IP_LINEAR
  double slopex, sloperho, slopet;   /**< pre-computed slopes for linear interpolation */
#elif INTERPOLATION_METHOD == IP_QUADRATIC
  double Ax1, Bx1, Cx1, Arho1, Brho1, Crho1, At1, Bt1, Ct1; /**< pre-computed constants for quadratic interpolation 1 */
  double Ax2, Bx2, Cx2, Arho2, Brho2, Crho2, At2, Bt2, Ct2; /**< pre-computed constants for quadratic interpolation 2 */
  double slopex, sloperho, slopet;   /**< pre-computed slopes for overlapping quadratics */
#endif
  double mint;                       /**< minimum T for current rho */
  int valid;                         /**< flag, set to 0 if current values in range of table */
};


/**
 * @brief Main EOS struct.
 * It contains the data tables and the X, rho and T values.
 *
 * EOS tables have the following dimensions:
 * - 169 points in rho
 * - 197 points in T
 * -   5 points in X
 */
struct opal_eos_table
{
  double z; /**< metallicity */

  /* eos variables */
  eostable p;      /**< pressure in megabars (10**12dyne/cm**2).
                    *   stored is P/P0, where P0=T6*rho             */
  eostable e;      /**< energy in 10**12 ergs/gm. Zero is zero T6.
                    *   stored is E/T6                              */
  eostable s;      /**< entropy in units of energy/T6
                    *   i.e. in 10**6 ergs/(g K)                    */
  eostable dedrho; /**< dE/dRHO at constant T6                      */
  eostable cv;     /**< specific heat, dE/dT6 at constant V         */
  eostable chir;   /**< dlogP/dlogRho at constant T6.
                    *   Cox and Guil1 eq 9.82                       */
  eostable chit;   /**< dlogP/dlogT6 at conxtant Rho.
                    *   Cox and Guil1 eq 9.81                       */
  eostable gamma1; /**< gamma1. Eqs. 9.88 Cox and Guili.            */
  eostable gamma2; /**< stored is gamma2/(gamma2-1).
                    *   Eqs. 9.88 Cox and Guili                     */


  double x[NX];     /**< X values (H mass fraction) */
  double rho[NRHO]; /**< rho values (in g/cm^-3) */
  double t[NT];     /**< T values (in 10^6 K) */

  int numt[NRHO];    /**< number of T points for each rho row */
  double minx;       /**< minimum X */
  double maxx;       /**< maximum X */
  double deltax;     /**< spacing in X */
  double minrho;     /**< minimum rho */
  double maxrho;     /**< maximum rho */
  double mint[NRHO]; /**< minimum T for each rho row */
  double maxt;       /**< maximum T */

  struct opal_eos_cache *cache; /**< eos cache object */
};


/* function declarations */
struct opal_eos_table *opaleos_init_local(const char *filename);
void opaleos_deinit_local(struct opal_eos_table *opal_table);
#ifdef OPAL_GLOBAL
int opaleos_init_global(const char *filename);
void opaleos_deinit_global();
void opaleos_get_rho_limits(double *lower, double *upper);
#endif
int opaleos(struct opal_eos_table *opal_table, double x, double rho, double t,
    int mode, struct opal_eos_result *result);
int opaleos_all_tgiven_local(struct opal_eos_table *opal_table, double x, double rho,
    double t, struct opal_eos_result *result);
int opaleos_all_egiven_local(struct opal_eos_table *opal_table, double x, double rho,
    double e, double *tguess, struct opal_eos_result *result);
int opaleos_egiven_local(struct opal_eos_table *opal_table, double x, double rho,
    double e, double *tguess, int mode, struct opal_eos_result *result);
int opaleos_all_pgiven_local(struct opal_eos_table *opal_table, double x, double rho,
    double p, double *tguess, struct opal_eos_result *result);
int opaleos_pgiven_local(struct opal_eos_table *opal_table, double x, double rho,
    double p, double *tguess, int mode, struct opal_eos_result *result);
#ifdef OPAL_GLOBAL
int opaleos_all_tgiven_global(double x, double rho, double t,
    struct opal_eos_result *result);
int opaleos_all_egiven_global(double x, double rho, double e, double *tguess,
    struct opal_eos_result *result);
int opaleos_egiven_global(double x, double rho, double e, double *tguess,
    int mode, struct opal_eos_result *result);
int opaleos_all_pgiven_global(double x, double rho, double p, double *tguess,
    struct opal_eos_result *result);
int opaleos_pgiven_global(double x, double rho, double p, double *tguess,
    int mode, struct opal_eos_result *result);
#endif
int get_eostable_value(struct opal_eos_table *opal_table, eostable var, double x, double rho, double t, double *result);

#ifdef OPAL_GLOBAL
#define opaleos_init opaleos_init_global
#define opaleos_deinit opaleos_deinit_global
#define opaleos_all_tgiven opaleos_all_tgiven_global
#define opaleos_all_egiven opaleos_all_egiven_global
#define opaleos_egiven opaleos_egiven_global
#define opaleos_all_pgiven opaleos_all_pgiven_global
#define opaleos_pgiven opaleos_pgiven_global
#else
#define opaleos_init opaleos_init_local
#define opaleos_deinit opaleos_deinit_local
#define opaleos_all_tgiven opaleos_all_tgiven_local
#define opaleos_all_egiven opaleos_all_egiven_local
#define opaleos_egiven opaleos_egiven_local
#define opaleos_all_pgiven opaleos_all_pgiven_local
#define opaleos_pgiven opaleos_pgiven_local
#endif

#endif /* EOS_OPAL */

#endif /* OPAL_EOS_H */
