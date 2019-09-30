/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/opal_eos.h
 * \date        02/2014
 * \author      Sebastian Ohlmann <sohlmann@astro.uni-wuerzburg.de>
 * \brief       source file for OPAL equation of state
 * \details     implementation of the OPAL equation of state in C using
 *    interpolation of the freely available data
 *
 *
 * \par Major modifications and contributions:
 *
 * - 10.02.2014 Integration in Arepo
 * - 15.07.2015 Improvement in convergence, added pgiven routines
 * - 21.04.2016 Enlarged range of EOS; includes blending to ideal EOS outside of
 *     OPAL tables
 * - 02.08.2016 Corrected EOS: now minimum T is 1.87e3 K, maximum T 2e8 K; only for
 *    rho < 1e14 g/cc, ideal gas + radiation is used; otherwise, T is set to value at
 *    boundary of table; converges now everywhere for egiven mode
 */

#include "arepoconfig.h"
#include "proto.h"

#ifdef EOS_OPAL

#undef USE_ZLIB

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef USE_ZLIB
# include <zlib.h>
#endif
#ifdef USE_ZLIB
# define fclose gzclose
# define fgets(x, y, z) gzgets(z, x, y)
#endif

#include "opal_eos.h"

#define LEFT 0
#define RIGHT 2
#define MIDDLE 1

#ifdef OPAL_GLOBAL
static struct opal_eos_table *opal_table; /**< main eos table object */
#endif


/* function declarations */
static inline double interpolate_linear(double x, double *xs, double *ys);
static inline double interpolate_linear_slope(double slope, double *ys);
static inline double interpolate_quadratic(double x, double *xs, double *ys);
static inline double interpolate_quadratic_coeff(double A, double B, double C, double *ys);
static inline void binary_search(double x, double *list, int len, int points, int *low, int *m);
static inline int linear_search_asc(double x, double *list, int len);
int checkmode(int mode);
double mean_molecular_weight(double xh, double z);
double mean_molecular_weight_neutral(double xh, double z);
double entropy_constant(double xh, double z, int neutral);
double get_temp_hion(double rho);
int update_cache(struct opal_eos_table *opal_table, double x, double rho, double t);
void compute_quadratic_coeff(double x, double *xs, double *A, double *B, double *C);
/* simple min/max functions */
/*
static inline int imax(int i1, int i2) { return i1 > i2 ? i1 : i2; }
static inline int imin(int i1, int i2) { return i1 < i2 ? i1 : i2; }
static inline double dmax(double i1, double i2) { return i1 > i2 ? i1 : i2; }
static inline double dmin(double i1, double i2) { return i1 < i2 ? i1 : i2; }
*/
int get_opal_result(struct opal_eos_table *opal_table, double x, double rho, double t,
    int mode, struct opal_eos_result *result);
int get_ideal_result(struct opal_eos_table *opal_table, double x, double rho, double t,
    int mode, struct opal_eos_result *result);
int check_limits(double rho, double t, double *F);
void blend_results(struct opal_eos_result *result, struct opal_eos_result *result_tmp,
    double F);
double outerborder(double rho);
double innerborder(double rho);
double get_blending_factor(double rho, double t);

/* constants for borders of OPAL region */
const double rho1 = 3e-4;
const double rho2 = 100;
const double rho3 = 3e5;
const double rho4 = 7e6;
const double t1 = 1.87e3;
const double t2 = 2.5e6;
const double t3 = 1.3e7;
const double rho1_i = 3e-4*0.9;
const double rho2_i = 100*0.9;
const double rho3_i = 3e5*0.9;
const double rho4_i = 7e6*0.9;
const double t1_i = 1.87e3*1.1;
const double t2_i = 2.5e6*1.1;
const double t3_i = 1.3e7*1.1;

/**
 * @brief initialization routine.
 * allocate table, open data file and read table from data file
 *
 * @param filename the name of the EOS data file
 *
 * @return the eos table struct
 */
struct opal_eos_table *opaleos_init_local(const char *filename)
{
  /* allocate table */
  struct opal_eos_table *opal_table = (struct opal_eos_table *) malloc(sizeof(struct opal_eos_table));

  if(opal_table == NULL)
    {
#ifdef OPAL_DEBUG
      fprintf(stderr, "Error allocating memory for EOS tables.\n");
#endif
      return NULL;
    }

  /* allocate cache */
  opal_table->cache = (struct opal_eos_cache *) malloc(sizeof(struct opal_eos_cache));

  if(opal_table->cache == NULL)
    {
#ifdef OPAL_DEBUG
      fprintf(stderr, "Error allocating memory for EOS cache.\n");
#endif
      return NULL;
    }
  /* initialize with negative values */
  opal_table->cache->x = -1.0;
  opal_table->cache->rho = -1.0;
  opal_table->cache->t = -1.0;
  opal_table->cache->mint = -1.0;

  /* open data file */
#ifndef USE_ZLIB
  FILE *datafile = fopen(filename, "r");
#else
  gzFile datafile = gzopen(filename, "r");
#endif

  if(datafile == NULL)
    {
#ifdef OPAL_DEBUG
      fprintf(stderr, "Error opening EOS data file %s.\n", filename);
#endif
      return NULL;
    }

  /* read in eos data */
  int ix, irho, it;
  int jrho;
  int n;
  char line[200];
  double dum;

  /* Reading EOS tables from data file */
  for(ix = 0; ix < NX; ix++)
    {
      /* read first line of new X block */
      if(fgets(line, 200, datafile) == NULL)
        {
#ifdef OPAL_DEBUG
          fprintf(stderr, "Error reading line.\n");
#endif
          fclose(datafile);
          return NULL;
        }
      n = sscanf(line, " X=%lf Z=%lf", &opal_table->x[ix], &opal_table->z);
      /* consistency checks */
      if(n < 2)
        {
#ifdef OPAL_DEBUG
          fprintf(stderr, "Error reading X line, read only %d elements.\n", n);
#endif
          fclose(datafile);
          return NULL;
        }
      /*printf("Reading values for X=%f.\n", opal_table->x[ix]); */
      /* 2 empty lines */
      fgets(line, 200, datafile);
      fgets(line, 200, datafile);

      for(irho = 0; irho < NRHO; irho++)
        {
          /* read first line of rho block */
          if(fgets(line, 200, datafile) == NULL)
            {
#ifdef OPAL_DEBUG
              fprintf(stderr, "Error reading line.\n");
#endif
              fclose(datafile);
              return NULL;
            }
          n = sscanf(line, "%d %d %lf %lf density(g/cc)= %lf", &jrho, &opal_table->numt[irho], &dum, &dum, &opal_table->rho[irho]);
          /* consistency checks */
          if(n < 5)
            {
#ifdef OPAL_DEBUG
              fprintf(stderr, "Error reading rho line, read only %d elements.\n", n);
#endif
              fclose(datafile);
              return NULL;
            }
          if((irho != jrho - 1) || (opal_table->numt[irho] > NT))
            {
#ifdef OPAL_DEBUG
              fprintf(stderr, "Error: inconsistent file format.\n");
#endif
              fclose(datafile);
              return NULL;
            }
          /* 2 empty lines */
          fgets(line, 200, datafile);
          fgets(line, 200, datafile);

          for(it = 0; it < opal_table->numt[irho]; it++)
            {
              /* read lines of T block */
              if(fgets(line, 200, datafile) == NULL)
                {
#ifdef OPAL_DEBUG
                  fprintf(stderr, "Error reading line.\n");
#endif
                  fclose(datafile);
                  return NULL;
                }
              /* read 12 values per line */
              n = sscanf(line, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
                         &opal_table->t[it], &dum, &dum, &opal_table->p[it][irho][ix],
                         &opal_table->e[it][irho][ix], &opal_table->s[it][irho][ix],
                         &opal_table->dedrho[it][irho][ix], &opal_table->cv[it][irho][ix],
                         &opal_table->chir[it][irho][ix], &opal_table->chit[it][irho][ix], &opal_table->gamma1[it][irho][ix], &opal_table->gamma2[it][irho][ix]);
              /* consistency checks */
              if(n < 12)
                {
#ifdef OPAL_DEBUG
                  fprintf(stderr, "Error reading T line, read only %d elements.\n", n);
#endif
                  fclose(datafile);
                  return NULL;
                }
            }
          /* 3 empty lines */
          fgets(line, 200, datafile);
          fgets(line, 200, datafile);
          fgets(line, 200, datafile);
        }
      /* final line with two zeros */
      fgets(line, 200, datafile);
    }

  /* close data file */
  fclose(datafile);

  /* get data boundaries */
  opal_table->minx = opal_table->x[0];
  opal_table->maxx = opal_table->x[NX - 1];
  opal_table->deltax = (opal_table->maxx - opal_table->minx) / (double) (NX - 1);
  opal_table->minrho = opal_table->rho[0];
  opal_table->maxrho = opal_table->rho[NRHO - 1];
  opal_table->maxt = opal_table->t[0];
  /* T min boundary depends on density */
  for(irho = 0; irho < NRHO; irho++)
    {
      opal_table->mint[irho] = opal_table->t[opal_table->numt[irho] - 1];
    }
  /* check boundaries */
  /*
     printf("X   - min: %e, max: %e\n", opal_table->minx, opal_table->maxx);
     printf("Rho - min: %e, max: %e\n", opal_table->minrho, opal_table->maxrho);
     printf("T   - min: %e, max: %e\n", opal_table->mint[0], opal_table->maxt);
   */

  printf("Read EOS tables from data file %s for metallicity %.3f.\n", filename, opal_table->z);
  return opal_table;
}

/**
 * @brief deinitialization routine.
 *
 * free EOS table data
 */
void opaleos_deinit_local(struct opal_eos_table *opal_table)
{
  free(opal_table->cache);
  free(opal_table);
}

#ifdef OPAL_GLOBAL
/**
 * @brief Initialization routine if global table is used
 *
 * @param filename Path to table file
 *
 * @return  <0 if errror
 */
int opaleos_init_global(const char *filename)
{
  /* read in file with task 0 and send to other tasks */
  if(ThisTask == 0)
    {
      opal_table = opaleos_init_local(filename);
    }
  else
    {
      /* allocate table */
      opal_table = (struct opal_eos_table *) malloc(sizeof(struct opal_eos_table));

      if(opal_table == NULL)
        {
#ifdef OPAL_DEBUG
          fprintf(stderr, "Error allocating memory for EOS table.\n");
#endif
          return -1;
        }

      /* allocate cache */
      opal_table->cache = (struct opal_eos_cache *) malloc(sizeof(struct opal_eos_cache));

      if(opal_table->cache == NULL)
        {
#ifdef OPAL_DEBUG
          fprintf(stderr, "Error allocating memory for EOS cache.\n");
#endif
          return -1;
        }
      /* initialize with negative values */
      opal_table->cache->x = -1.0;
      opal_table->cache->rho = -1.0;
      opal_table->cache->t = -1.0;
      opal_table->cache->mint = -1.0;
    }
  /* now broadcast table data */
  MPI_Bcast(&opal_table->z, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Bcast(opal_table->p, NT * NRHO * NX, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(opal_table->e, NT * NRHO * NX, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(opal_table->s, NT * NRHO * NX, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(opal_table->dedrho, NT * NRHO * NX, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(opal_table->cv, NT * NRHO * NX, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(opal_table->chir, NT * NRHO * NX, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(opal_table->chit, NT * NRHO * NX, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(opal_table->gamma1, NT * NRHO * NX, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(opal_table->gamma2, NT * NRHO * NX, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Bcast(opal_table->x, NX, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(opal_table->rho, NRHO, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(opal_table->t, NT, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Bcast(opal_table->numt, NRHO, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&opal_table->minx, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&opal_table->maxx, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&opal_table->deltax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&opal_table->minrho, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&opal_table->maxrho, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(opal_table->mint, NRHO, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&opal_table->maxt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  return 0;
}

/**
 * @brief deinitialization routine for global table
 */
void opaleos_deinit_global()
{
  opaleos_deinit_local(opal_table);
}

/**
 * @brief get limits of table in density
 *
 * @param lower lower border
 * @param upper upper border
 */
void opaleos_get_rho_limits(double *lower, double *upper)
{
  *lower = opal_table->minrho;
  *upper = opal_table->maxrho;
}
#endif

/**
 * @brief main driver routine for the OPAL EOS
 *
 * @param opal_table the opal table struct
 * @param x H mass fraction
 * @param t temperature in K
 * @param rho density in g/cm^3
 * @param mode bit mask deciding which values to retrieve
 * @param result pointer to double array where results are stored
 *
 * @result negative if error occured, otherwise zero
 */
int opaleos(struct opal_eos_table *opal_table, double x, double rho, double t, int mode, struct opal_eos_result *result)
{
  int err = 0;
  int do_opal = 1; /* > determines if input values are in opal range; 1: OPAL, 0: ideal, -1: blend */
  double res;
  double prad, erad, srad, p0, gamma1_0, gamma2_0, gamma3_0, gamma1_tot, gamma2_tot, gamma3_tot;
  double t1, t2;
  double F; /**> blending factor */

  /* interpolation is done in T6 */
  double t6 = t * 1e-6;

  /* check consistency of EOS mode */
  mode = checkmode(mode);

  /* check limits */
  do_opal = check_limits(rho, t, &F);

  /* check range in rho, T, perhaps also Q?
   * compare to blending range of MESA
   * SCVH tables?
   * FreeEOS?
   * blend to Helmholtz? Create helmholtz tables to lower temperatures, as
   *   ideal gas plus radiation?
   * extend to lower densities < 1e-14: HELM?/ideal gas?
   */

  if (do_opal == 1)
    {
      /* update cache and return with error if value not in table */
      if(update_cache(opal_table, x, rho, t6) < 0)
        {
          if ((x < opal_table->minx) || (x > opal_table->maxx))
            return -1;
        }
      err = get_opal_result(opal_table, x, rho, t, mode, result);
    }
  else if (do_opal == 0) /* out of opal range, do ideal gas */
    {
      err = get_ideal_result(opal_table, x, rho, t, mode, result);
    }
  else if (do_opal == -1) /* in border region; blend OPAL and ideal gas*/
    {
      //struct opal_eos_result *result_tmp = malloc(sizeof(struct opal_eos_result));
      /* blend */
      err = get_opal_result(opal_table, x, rho, t, mode, result);
      //err = get_ideal_result(opal_table, x, rho, t, mode, result_tmp);
      //blend_results(result, result_tmp, F);
      //free(result_tmp);
    }

  /* radiation corrections */
  if(mode & EOS_RADIATION)
    {
      /* quantities for radiation; in cgs units */
      prad = 4.0 / 3.0 * SIGMAC * t * t * t * t;
      erad = 3.0 * prad / rho;
      //srad = 4.0 / 3.0 * erad / t;
      srad = 16.0 / 3.0 * SIGMAC * t * t * t / rho;

      /* store quantities without radiation */
      if(mode & EOS_CHIR)
        p0 = result->p;
      if((mode & EOS_GAMMA1) || (mode & EOS_GAMMA2) || (mode & EOS_GAMMA3))
        {
          gamma3_0 = 1.0 + result->p * result->chit / (t * rho * result->cv);
          gamma1_0 = result->chir + result->chit * (gamma3_0 - 1.0);
          gamma2_0 = gamma1_0 / (gamma1_0 - gamma3_0 + 1.0);
        }

      /* update quantities */
      if(mode & EOS_PRESSURE)
        result->p += prad;
      if(mode & EOS_ENERGY)
        result->e += erad;
      if(mode & EOS_ENTROPY)
        result->s += srad;
      if(mode & EOS_DEDRHO)
        result->dedrho -= erad / rho;
      if(mode & EOS_CV)
        result->cv += 3.0 * srad;
      if(mode & EOS_CHIR)
        result->chir = p0 * result->chir / result->p;
      if(mode & EOS_CHIT)
        result->chit = (p0 * result->chit + prad * 4.0) / result->p;
      if((mode & EOS_GAMMA1) || (mode & EOS_GAMMA2) || (mode & EOS_GAMMA3))
        {
          gamma3_tot = 1.0 + result->p * result->chit / (t * rho * result->cv);
          gamma1_tot = result->chir + result->chit * (gamma3_tot - 1.0);
          gamma2_tot = gamma1_tot / (gamma1_tot - gamma3_tot + 1.0);
        }
      if(mode & EOS_GAMMA1)
        result->gamma1 += gamma1_tot - gamma1_0;
      if(mode & EOS_GAMMA2)
        result->gamma2 += gamma2_tot - gamma2_0;
      if(mode & EOS_GAMMA3)
        result->gamma3 += gamma3_tot - gamma3_0;
    }

  return 0;
}

/**
 * @brief routine for retrieving the tabulated OPAL values
 *
 * @param opal_table the opal table struct
 * @param x H mass fraction
 * @param t temperature in K
 * @param rho density in g/cm^3
 * @param mode bit mask deciding which values to retrieve
 * @param result pointer to double array where results are stored
 *
 * @result negative if error occured, otherwise zero
 */
int get_opal_result(struct opal_eos_table *opal_table, double x, double rho, double t,
    int mode, struct opal_eos_result *result)
{
  int err;
  double res;
  double t6;
  err = 0;
  t6 = t/1e6;
  /* retrieve eos variables depending on mode */
  if(mode & EOS_PRESSURE)
    {
      err += get_eostable_value(opal_table, opal_table->p, x, rho, t6, &res);
      /* stored is P/P0, where P0=T6*rho, in 10^12 dyn/cm^2 */
      /* return pressure in dyn/cm^2 */
      result->p = res * t6 * rho * 1e12;;
    }
  if(mode & EOS_ENERGY)
    {
      err += get_eostable_value(opal_table, opal_table->e, x, rho, t6, &res);
      /* stored is E/T6 in 10^12 erg/g */
      /* return energy in erg/g */
      result->e = res * t6 * 1e12;
    }
  if(mode & EOS_ENTROPY)
    {
      err += get_eostable_value(opal_table, opal_table->s, x, rho, t6, &res);
      /* stored in 10^6 erg/(g K) */
      /* return entropy in erg/(g K) */
      result->s = res * 1e6;
    }
  if(mode & EOS_DEDRHO)
    {
      err += get_eostable_value(opal_table, opal_table->dedrho, x, rho, t6, &res);
      result->dedrho = res;
    }
  if(mode & EOS_CV)
    {
      err += get_eostable_value(opal_table, opal_table->cv, x, rho, t6, &res);
      /* cv is normalized to 3/2 for ideal, non-degenerate, fully ionized gas */
      /* formula for ideal gas: cv = 3/2 * R / mu */
      result->cv = res * GAS_CONSTANT / mean_molecular_weight(x, opal_table->z);
    }
  if(mode & EOS_CHIR)
    {
      err += get_eostable_value(opal_table, opal_table->chir, x, rho, t6, &res);
      result->chir = res;
    }
  if(mode & EOS_CHIT)
    {
      err += get_eostable_value(opal_table, opal_table->chit, x, rho, t6, &res);
      result->chit = res;
    }
  if(mode & EOS_GAMMA1)
    {
      err += get_eostable_value(opal_table, opal_table->gamma1, x, rho, t6, &res);
      result->gamma1 = res;
    }
  if(mode & EOS_GAMMA2)
    {
      err += get_eostable_value(opal_table, opal_table->gamma2, x, rho, t6, &res);
      /* stored is gamma2/(gamma2-1) */
      result->gamma2 = res / (res - 1.0);
    }
  if(mode & EOS_GAMMA3)
    {
      /* compute from gamma1 and gamma2 */
      result->gamma3 = 1.0 + result->gamma1 - result->gamma1 / result->gamma2;
    }

  /* check for errors */
  if(err < 0)
    {
      fprintf(stderr, "Error obtaining EOS values for x=%e, rho=%e, T=%e.\n", x, rho, t);
      return -1;
    }
  return 0;
}

/**
 * @brief routine for computing the ideal EOS
 *
 * @param x H mass fraction
 * @param t temperature in K
 * @param rho density in g/cm^3
 * @param mode bit mask deciding which values to retrieve
 * @param result pointer to double array where results are stored
 *
 * @result negative if error occured, otherwise zero
 */
int get_ideal_result(struct opal_eos_table *opal_table, double x, double rho, double t,
    int mode, struct opal_eos_result *result)
{
  double t6;
  t6 = t/1e6;

  /* check T: neutral for T<T(X(H I)==0.5), ionized otherwise */
  double txh1 = get_temp_hion(rho);
  double mu;
  if (t6*1e6 < txh1)
    mu = mean_molecular_weight_neutral(x, opal_table->z);
  else
    mu = mean_molecular_weight(x, opal_table->z);
  result->p = rho * GAS_CONSTANT / mu * t6 * 1e6;
  result->e = 1.5 * GAS_CONSTANT / mu * t6 * 1e6;
  result->s = GAS_CONSTANT/mu * log(pow(t6*1e6, 1.5)/rho);
  if (t6*1e6 < txh1)
    result->s += entropy_constant(x, opal_table->z, 1);
  else
    result->s += entropy_constant(x, opal_table->z, 0);
  result->s = dmax(0, result->s);
  result->dedrho = 0;
  result->cv = 1.5 * GAS_CONSTANT / mu;
  result->chir = 1;
  result->chit = 1;
  result->gamma1 = 5./3.;
  result->gamma2 = 5./3.;
  result->gamma3 = 5./3.;
  return 0;
}
/**
 * @brief driver routine for retrieving all eos values given T
 *
 * @param opal_table the eos table struct
 * @param x hydrogen mass fraction
 * @param rho density
 * @param t temperature
 * @param result eos result structure
 *
 * @return <0 if error occured, 0 otherwise
 */
int opaleos_all_tgiven_local(struct opal_eos_table *opal_table, double x, double rho, double t, struct opal_eos_result *result)
{
  /*
  */
  double t6, t1, t2;
  t6 = t * 1e-6;
  update_cache(opal_table, x, rho, t6);
  //t1 = opal_table->cache->mint * 1e6 * 1.0001;
  t1 =1.87e3;
  t2 = opal_table->maxt * 1e6 * 0.9999;
  //if (t < t1) t = t1;
  if (t < t1 && rho > opal_table->minrho) t = t1;
  if (t > t2) t = t2;
  t6 = t * 1e-6;
  update_cache(opal_table, x, rho, t6);
  return opaleos(opal_table, x, rho, t, EOS_ALL, result);
}

/**
 * @brief driver routine for retrieving all eos values given T
 *  for global tables
 *
 * @param x hydrogen mass fraction
 * @param rho density
 * @param t temperature
 * @param result eos result structure
 *
 * @return <0 if error occured, 0 otherwise
 */
#ifdef OPAL_GLOBAL
int opaleos_all_tgiven_global(double x, double rho, double t, struct opal_eos_result *result)
{
  /*
  */
  double t6, t1, t2;
  t6 = t * 1e-6;
  update_cache(opal_table, x, rho, t6);
  //t1 = opal_table->cache->mint * 1e6 * 1.0001;
  t1 =1.87e3;
  t2 = opal_table->maxt * 1e6 * 0.9999;
  //if (t < t1) t = t1;
  if (t < t1 && rho > opal_table->minrho) t = t1;
  if (t > t2) t = t2;
  t6 = t * 1e-6;
  update_cache(opal_table, x, rho, t6);
  return opaleos(opal_table, x, rho, t, EOS_ALL, result);
}
#endif


/**
 * @brief driver routine for retrieving all eos values given the energy
 *
 * @param opal_table the eos table struct
 * @param x hydrogen mass fraction
 * @param rho density
 * @param e internal energy
 * @param tguess guess for temperature, if <0 make guess for ideal gas
 * @param result eos result structure
 *
 * @return <0 if error occured, 0 otherwise
 */
int opaleos_all_egiven_local(struct opal_eos_table *opal_table, double x, double rho, double e, double *tguess, struct opal_eos_result *result)
{
  return opaleos_egiven_local(opal_table, x, rho, e, tguess, EOS_ALL, result);
}

/**
 * @brief driver routine for retrieving certain eos values given the energy
 *
 * @param opal_table the eos table struct
 * @param x hydrogen mass fraction
 * @param rho density
 * @param e internal energy
 * @param tguess guess for temperature, if <0 make guess for ideal gas
 * @param mode only retrieve values given by mode
 * @param result eos result structure
 *
 * @return <0 if error occured, 0 otherwise
 */
int opaleos_egiven_local(struct opal_eos_table *opal_table, double x, double rho, double e, double *tguess, int mode, struct opal_eos_result *result)
{
  double t, told;               /* temperature */
  double de;                    /* change in energy */
  int i, err;
  int num_out;

  /* if no guess given, use ideal gas with radiation */
  if(*tguess < 0.0)
    {
      double A = 3.0 / 2.0 * GAS_CONSTANT / mean_molecular_weight(x, opal_table->z);
      double B = 4.0 * SIGMAC / rho;
      double T0i = e / A;       /* ideal gas */
      double Ti = T0i * (1.0 - 1.0 / (4.0 + A / (B * T0i * T0i * T0i)));        /* 1st order correction */
      double T0r = sqrt(sqrt(e / B));   /* radiation */
      double Tr = T0r * (1.0 - 1.0 / (1.0 + 4.0 * B * T0i * T0i * T0i / A));    /* 1st order correction */
      /*
         printf("T0i: %e, Ti: %e, T0r: %e, Tr: %e\n", T0i, Ti, T0r, Tr);
         printf("ei: %e, er: %e\n", Ti * A, B * Tr * Tr * Tr * Tr);
         printf("ei,tr: %e, er,ti: %e\n", Tr * A, B * Ti * Ti * Ti * Ti);
       */
      double ei = A * Ti + B * Ti * Ti * Ti * Ti;
      double er = A * Tr + B * Tr * Tr * Tr * Tr;

      /* guess for temperature with smallest difference */
      if(fabs(e - ei) < fabs(e - er))
        t = Ti;
      else
        t = Tr;
    }
  else
    t = *tguess;

  /* check bounds for density, if out of bounds move inside -> for robustness! */
  /* TODO disable?
  if (rho > opal_table->maxrho)
    rho = opal_table->maxrho;
  else if (rho < opal_table->minrho)
    rho = opal_table->minrho;
  */

  /* Newton loop */
  /* counter for iterations outside ranges */
  num_out = 0;
  for(i = 0; i < OPAL_NEWTON_MAX; i++)
    {
      err = opaleos(opal_table, x, rho, t, EOS_ENERGY | EOS_CV | EOS_RADIATION, result);
      de = (result->e - e);
      told = t;
      t -= de / result->cv;
      /* check boundaries */
      /* check bounds */
      /* old version :
         if(t > opal_table->maxt * 1e6)
         t = 0.5 * (told + opal_table->maxt * 1e6);
         else if(t < opal_table->cache->mint * 1e6)
         t = 0.5 * (told + opal_table->cache->mint * 1e6);
       */
      if(t > opal_table->maxt * 1e6)
        {
          t = opal_table->maxt * 1e6;
          num_out++;
        }
      else if(t < 1.87e3)
        {
          t = 1.87e3;
          num_out++;
        }
      /*
      else if(t < opal_table->cache->mint * 1e6)
        {
          t = opal_table->cache->mint * 1e6;
          num_out++;
        }
        */

      /* converged to desired accuracy? */
      if(fabs(de) <= OPAL_NEWTON_ACCURACY * e)
        break;

      /* out of range? */
      if(err < 0)               /* error in EOS */
        {
#ifdef OPAL_DEBUG
          fprintf(stderr, "Error in iteration %d for X=%e, rho=%e, e=%e, t=%e, mint=%e, maxt=%e\n", i, x, rho, e, t, opal_table->cache->mint, opal_table->maxt);
//          return -1;
#endif
        }

      /* check number of iterations outside ranges, if too large, do bisection */
      if (num_out > 8)
        break;
#ifdef OPAL_DEBUG
      fprintf(stderr, "Iteration %d for X=%e, rho=%e, e=%e => t=%e, e_n=%e, de=%e\n", i, x, rho, e, t, result->e, de);
#endif
    }

  /* do bisection if newton has not converged */
  if(num_out > 5 || i >= OPAL_NEWTON_MAX)
    {
      num_out = 0;
#ifdef OPAL_DEBUG
      fprintf(stderr, "Switching to bisection for X=%e, rho=%e, e=%e\n", x, rho, e);
#endif
      double t1, t2, e1, e2, ttry;

      /* get energy at boundaries */
      //update_cache(opal_table, x, rho, opal_table->maxt);
      //t1 = opal_table->cache->mint * 1e6 * 1.0001;
      t1 = 1.87e3;
      err = opaleos(opal_table, x, rho, t1, EOS_ENERGY | EOS_RADIATION, result);
      e1 = result->e;
#ifdef OPAL_DEBUG
      fprintf(stderr, "X=%e, rho=%e, e=%e, t=%e\n", x, rho, result->e, t1);
#endif
      t2 = opal_table->maxt * 1e6 * 0.9999;
      err = opaleos(opal_table, x, rho, t2, EOS_ENERGY | EOS_RADIATION, result);
      e2 = result->e;
#ifdef OPAL_DEBUG
      fprintf(stderr, "X=%e, rho=%e, e=%e, t=%e\n", x, rho, result->e, t2);
#endif
      /* check range of energy */
      if(e < e1 || e > e2)
        {
          if(e < e1)
            t = t1;
          if(e > e2)
            t = t2;
#ifdef OPAL_DEBUG
          fprintf(stderr, "Error in EOS: e out of range for X=%e, rho=%e, e=%e => t=%e.\n", x, rho, e, t);
#endif
          *tguess = t;
          opaleos(opal_table, x, rho, t, mode, result);
          return -1;
        }

      /* check bounds, if out of bounds move inside */
      /*
      if (t > opal_table->maxt * 1e6)
        t = pow(10, log10(opal_table->maxt * 1e6)
            - 0.1 * log10(opal_table->maxt / opal_table->cache->mint));
      else if (t < opal_table->cache->mint * 1e6)
        t = pow(10, log10(opal_table->cache->mint * 1e6)
            + 0.1 * log10(opal_table->maxt / opal_table->cache->mint));
            */

      /* find better initial bracketing */
      err = opaleos(opal_table, x, rho, t, EOS_ENERGY | EOS_RADIATION, result);
      de = (result->e - e);
      if(de * (e2 - e) > 0)
        {
          t2 = t;
          //ttry = fmax(t/3., opal_table->cache->mint * 1e6);
          ttry = fmax(t/3., 100);
        }
      else
        {
          t1 = t;
          ttry = fmin(t * 3., opal_table->maxt * 1e6);
        }
      err = opaleos(opal_table, x, rho, ttry, EOS_ENERGY | EOS_RADIATION, result);
#ifdef OPAL_DEBUG
      fprintf(stderr, "X=%e, rho=%e, e=%e, t=%e\n\n\n", x, rho, result->e, ttry);
#endif
      de = (result->e - e);
      if(de * (e2 - e) > 0)
        t2 = ttry;
      else
        t1 = ttry;

      /* now set start value for loop */
      t = 0.5 * (t1 + t2);

      /* bisection loop */
      for(i = 0; i < OPAL_NEWTON_MAX; i++)
        {
          err = opaleos(opal_table, x, rho, t, EOS_ENERGY | EOS_RADIATION, result);
          de = (result->e - e);
          if(de * (e2 - e) > 0)
            {
              t2 = t;
            }
          else if(de * (e1 - e) > 0)
            {
              t1 = t;
            }
          else
#ifdef OPAL_DEBUG
            fprintf(stderr, "Error in iteration %d for X=%e, rho=%e, e=%e, t=%e, mint=%e, maxt=%e\n", i, x, rho, e, t, opal_table->cache->mint, opal_table->maxt);
#endif

          told = t;
          t = 0.5 * (t1 + t2);
          /* checking boundaries not necessary because of bisection */

          /* converged to desired accuracy? */
          if(fabs(de) <= OPAL_NEWTON_ACCURACY * e)
            break;

          /* bisection stopped at border? */
          if(fabs(told - t) <= OPAL_NEWTON_ACCURACY * t)
            num_out++;

          if(num_out > 10)
            break;

          /* out of range? */
          if(err < 0)           /* error in EOS */
            {
#ifdef OPAL_DEBUG
              fprintf(stderr, "Error in iteration %d for X=%e, rho=%e, e=%e, t=%e, mint=%e, maxt=%e\n", i, x, rho, e, t, opal_table->cache->mint, opal_table->maxt);
#endif
              //          return -1;
            }
#ifdef OPAL_DEBUG
          fprintf(stderr, "Iteration %d for X=%e, rho=%e, e=%e => t=%e, e_n=%e, de=%e\n", i, x, rho, e, t, result->e, de);
#endif
        }
    }                           /* end of bisection loop */


  if(i >= OPAL_NEWTON_MAX || num_out > 5)
    {
#ifdef OPAL_DEBUG
      fprintf(stderr, "Error in EOS: Newton-Raphson or bisection not converged for X=%e, rho=%e, e=%e => t=%e.\n", x, rho, e, t);
#endif
      *tguess = t;
      opaleos(opal_table, x, rho, t, mode, result);
      return -1;
    }

  *tguess = t;
  return opaleos(opal_table, x, rho, t, mode, result);
}

#ifdef OPAL_GLOBAL
/**
 * @brief driver routine for retrieving all eos values given the energy
 *
 * @param x hydrogen mass fraction
 * @param rho density
 * @param e internal energy
 * @param tguess guess for temperature, if <0 make guess for ideal gas
 * @param result eos result structure
 *
 * @return <0 if error occured, 0 otherwise
 */
int opaleos_all_egiven_global(double x, double rho, double e, double *tguess, struct opal_eos_result *result)
{
  return opaleos_all_egiven_local(opal_table, x, rho, e, tguess, result);
}

int opaleos_egiven_global(double x, double rho, double e, double *tguess, int mode, struct opal_eos_result *result)
{
  return opaleos_egiven_local(opal_table, x, rho, e, tguess, mode, result);
}
#endif

/**
 * @brief driver routine for retrieving all eos values given the pressure
 *
 * @param opal_table the eos table struct
 * @param x hydrogen mass fraction
 * @param rho density
 * @param p pressure
 * @param tguess guess for temperature, if <0 make guess for ideal gas
 * @param result eos result structure
 *
 * @return <0 if error occured, 0 otherwise
 */
int opaleos_all_pgiven_local(struct opal_eos_table *opal_table, double x, double rho, double p, double *tguess, struct opal_eos_result *result)
{
  return opaleos_pgiven_local(opal_table, x, rho, p, tguess, EOS_ALL, result);
}

/**
 * @brief driver routine for retrieving certain eos values given the energy
 *
 * @param opal_table the eos table struct
 * @param x hydrogen mass fraction
 * @param rho density
 * @param p pressure
 * @param tguess guess for temperature, if <0 make guess for ideal gas
 * @param mode only retrieve values given by mode
 * @param result eos result structure
 *
 * @return <0 if error occured, 0 otherwise
 */
int opaleos_pgiven_local(struct opal_eos_table *opal_table, double x, double rho, double p, double *tguess, int mode, struct opal_eos_result *result)
{
  double t, told;               /* temperature */
  double dp;                    /* change in pressure */
  double egas, erad;
  int i, err;
  double t1, t2, p1, p2, ttry;

  /* get pressure at boundaries */
  //update_cache(opal_table, x, rho, opal_table->maxt);
  //t1 = opal_table->cache->mint * 1e6 * 1.0001;
  t1 = 10;
  err = opaleos(opal_table, x, rho, t1, EOS_PRESSURE | EOS_RADIATION, result);
  p1 = result->p;
#ifdef OPAL_DEBUG
  fprintf(stderr, "X=%e, rho=%e, p=%e, t=%e, chit=%e\n", x, rho, result->p, t1, result->chit);
#endif
  t2 = opal_table->maxt * 1e6 * 0.9999;
  err = opaleos(opal_table, x, rho, t2, EOS_PRESSURE | EOS_RADIATION, result);
  p2 = result->p;
#ifdef OPAL_DEBUG
  fprintf(stderr, "X=%e, rho=%e, p=%e, t=%e, chit=%e\n\n\n", x, rho, result->p, t2, result->chit);
#endif
  /* check range of pressure */
  if(p < p1 || p > p2)
    {
      if(p < p1)
        t = t1;
      if(p > p2)
        t = t2;
#ifdef OPAL_DEBUG
      fprintf(stderr, "Error in EOS: p out of range for X=%e, rho=%e, p=%e => t=%e.\n", x, rho, p, t);
#endif
      *tguess = t;
      opaleos(opal_table, x, rho, t, mode, result);
      return -1;
    }

  /* if no guess given, use ideal gas with radiation */
  if(*tguess < 0.0)
    {
      double A = 3.0 / 2.0 * GAS_CONSTANT / mean_molecular_weight(x, opal_table->z);
      double B = 4.0 * SIGMAC / rho;
      egas = 3. / 2. * p / rho;
      double T0i = egas / A;    /* ideal gas */
      double Ti = T0i * (1.0 - 1.0 / (4.0 + A / (B * T0i * T0i * T0i)));        /* 1st order correction */
      erad = 3. * p / rho;
      double T0r = sqrt(sqrt(erad / B));        /* radiation */
      double Tr = T0r * (1.0 - 1.0 / (1.0 + 4.0 * B * T0i * T0i * T0i / A));    /* 1st order correction */
#ifdef OPAL_DEBUG
      printf("T0i: %e, Ti: %e, T0r: %e, Tr: %e\n", T0i, Ti, T0r, Tr);
      printf("ei: %e, er: %e\n", Ti * A, B * Tr * Tr * Tr * Tr);
      printf("ei,tr: %e, er,ti: %e\n", Tr * A, B * Ti * Ti * Ti * Ti);
      printf("pi: %e, pr: %e\n", Ti * A * rho * 2. / 3., B * Tr * Tr * Tr * Tr * rho / 3.);
#endif
      double pi = A * Ti * rho * 2. / 3. + B * Ti * Ti * Ti * Ti * rho / 3.;
      double pr = A * Tr * rho * 2. / 3. + B * Tr * Tr * Tr * Tr * rho / 3.;

      /* guess for temperature with smallest difference */
      if(fabs(p - pi) < fabs(p - pr))
        t = Ti;
      else
        t = Tr;
    }
  else
    t = *tguess;

  /* check bounds, if out of bounds move inside */
  /*
  if (t > opal_table->maxt * 1e6)
    t = pow(10, log10(opal_table->maxt * 1e6)
        - 0.1 * log10(opal_table->maxt / opal_table->cache->mint));
  else if (t < opal_table->cache->mint * 1e6)
    t = pow(10, log10(opal_table->cache->mint * 1e6)
        + 0.1 * log10(opal_table->maxt / opal_table->cache->mint));
  */

  /* find better initial bracketing */
  err = opaleos(opal_table, x, rho, t, EOS_PRESSURE | EOS_RADIATION, result);
  dp = (result->p - p);
  if(dp * (p2 - p) > 0)
    {
      t2 = t;
      //ttry = fmax(t/3., opal_table->cache->mint * 1e6);
      ttry = fmax(t/3., 10);
    }
  else
    {
      t1 = t;
      ttry = fmin(t * 3., opal_table->maxt * 1e6);
    }
  err = opaleos(opal_table, x, rho, ttry, EOS_PRESSURE | EOS_RADIATION, result);
#ifdef OPAL_DEBUG
  fprintf(stderr, "X=%e, rho=%e, p=%e, t=%e\n\n\n", x, rho, result->p, ttry);
#endif
  dp = (result->p - p);
  if(dp * (p2 - p) > 0)
    t2 = ttry;
  else
    t1 = ttry;

  /* now set start value for loop */
  t = 0.5 * (t1 + t2);

  /* bisection loop */
  for(i = 0; i < OPAL_NEWTON_MAX; i++)
    {
      err = opaleos(opal_table, x, rho, t, EOS_PRESSURE | EOS_RADIATION, result);
      dp = (result->p - p);
      if(dp * (p2 - p) > 0)
        {
          t2 = t;
        }
      else if(dp * (p1 - p) > 0)
        {
          t1 = t;
        }
      else
#ifdef OPAL_DEBUG
        fprintf(stderr, "Error in iteration %d for X=%e, rho=%e, p=%e, t=%e, mint=%e, maxt=%e\n", i, x, rho, p, t, opal_table->cache->mint, opal_table->maxt);
#endif

      told = t;
      t = 0.5 * (t1 + t2);
      /* checking boundaries not necessary because of bisection */

      /* converged to desired accuracy? */
      if(fabs(dp) <= OPAL_NEWTON_ACCURACY * p)
        break;

      /* out of range? */
      if(err < 0)               /* error in EOS */
        {
#ifdef OPAL_DEBUG
          fprintf(stderr, "Error in iteration %d for X=%e, rho=%e, p=%e, t=%e, mint=%e, maxt=%e\n", i, x, rho, p, t, opal_table->cache->mint, opal_table->maxt);
#endif
//          return -1;
        }
#ifdef OPAL_DEBUG
      fprintf(stderr, "Iteration %d for X=%e, rho=%e, p=%e => t=%e, dt=%e, p_n=%e, dp=%e, chit=%e\n", i, x, rho, p, t, dp / (result->p / t * result->chit), result->p, dp, result->chit);
#endif
    }

  if(i >= OPAL_NEWTON_MAX)
    {
#ifdef OPAL_DEBUG
      fprintf(stderr, "Error in EOS: Newton-Raphson not converged for X=%e, rho=%e, p=%e => t=%e.\n", x, rho, p, t);
#endif
      *tguess = t;
      opaleos(opal_table, x, rho, t, mode, result);
      return -1;
    }

  *tguess = t;
  return opaleos(opal_table, x, rho, t, mode, result);
}

#ifdef OPAL_GLOBAL
/**
 * @brief driver routine for retrieving all eos values given the pressure
 *
 * @param x hydrogen mass fraction
 * @param rho density
 * @param p pressure
 * @param tguess guess for temperature, if <0 make guess for ideal gas
 * @param result eos result structure
 *
 * @return <0 if error occured, 0 otherwise
 */
int opaleos_all_pgiven_global(double x, double rho, double p, double *tguess, struct opal_eos_result *result)
{
  return opaleos_all_pgiven_local(opal_table, x, rho, p, tguess, result);
}

int opaleos_pgiven_global(double x, double rho, double p, double *tguess, int mode, struct opal_eos_result *result)
{
  return opaleos_pgiven_local(opal_table, x, rho, p, tguess, mode, result);
}
#endif


/**
 * @brief retrieve value of certain variable
 *
 * @param opal_table the opal table struct
 * @param var EOS variable to interpolate
 * @param x H mass fraction
 * @param rho density in g/cm^3
 * @param t temperature in 10^6 K
 * @param result interpolated value
 *
 * @result negative if error occured, otherwise zero
 */
int get_eostable_value(struct opal_eos_table *opal_table, eostable var, double x, double rho, double t, double *result)
{
  int i, j;

  /* update cache */
  if(update_cache(opal_table, x, rho, t) < 0)
    return -1;

#if INTERPOLATION_METHOD == IP_LINEAR
  /* interpolate in x first */
  double varrt[2][2];
  /* use cached slopes for faster interpolation */
  for(i = 0; i < 2; i++)
    for(j = 0; j < 2; j++)
      varrt[i][j] = interpolate_linear_slope(opal_table->cache->slopex, var[opal_table->cache->it + i][opal_table->cache->irho + j] + opal_table->cache->ix);

  /* interpolate in rho */
  double vart[2];
  for(i = 0; i < 2; i++)
    vart[i] = interpolate_linear_slope(opal_table->cache->sloperho, varrt[i]);

  /* interpolate in T */
  //*result = interpolate_linear(t, opal_table->t+it, vart);
  *result = interpolate_linear_slope(opal_table->cache->slopet, vart);
#elif INTERPOLATION_METHOD == IP_QUADRATIC
  /* Determine position */
  int rhoposition;              /* 0: left, 1: middle, 2:right */
  if((opal_table->cache->mrho == 0) || (opal_table->cache->mt >= opal_table->numt[opal_table->cache->mrho + 1] - 1))
    rhoposition = LEFT;
  /* right border */
  else if(opal_table->cache->mrho >= (NRHO - 2))
    rhoposition = RIGHT;
  /* in the middle */
  else
    rhoposition = MIDDLE;

  /* max indices for right T border */
  int imaxt1 = imin(opal_table->numt[opal_table->cache->irho],
                    imin(opal_table->numt[opal_table->cache->irho + 1],
                         opal_table->numt[opal_table->cache->irho + 2])) - 1;
  int imaxt2 = imin(opal_table->numt[opal_table->cache->irho + 1],
                    imin(opal_table->numt[opal_table->cache->irho + 2],
                         opal_table->numt[opal_table->cache->irho + 3])) - 1;

  /* conditions for variable T border */
  int topleft, topright, bottomleft, bottomright, left, right, top, bottom;

  /* when do we need the corresponding interpolation? */
  topleft = !((rhoposition == RIGHT) || (opal_table->cache->mt == opal_table->cache->it + 2)
              //|| (imaxt1 < opal_table->cache->it + 2)
    );
  bottomleft = !((rhoposition == RIGHT) || (opal_table->cache->mt == opal_table->cache->it)
                 //|| (imaxt1 < opal_table->cache->it + 3)
    );
  topright = !((rhoposition == LEFT) || (opal_table->cache->mt == opal_table->cache->it + 2) || (imaxt2 < opal_table->cache->it + 2));
  bottomright = !((rhoposition == LEFT) || (opal_table->cache->mt == opal_table->cache->it) || (imaxt2 < opal_table->cache->it + 3));

  left = topleft || bottomleft;
  right = topright || bottomright;

  top = topleft || topright;
  bottom = bottomleft || bottomright;

  if((!left) && (!right))
    return -1;



  /* interpolate in x first */
  double varrt[4][4], varrt1[4][4], varrt2[4][4];
  /* use cached slopes for faster interpolation */
  /* left border */
  if(opal_table->cache->mx == 0)
    {
      for(i = 0; i < 4; i++)
        {
          if(i == 0 && !top)
            continue;
          if(i == 4 && !bottom)
            continue;
          for(j = 0; j < 4; j++)
            {
              if(j == 0)
                {
                  if(!left)
                    continue;
                  if(i == 0 && !topleft)
                    continue;
                  if(i == 4 && !bottomleft)
                    continue;
                }
              if(j == 4)
                {
                  if(!right)
                    continue;
                  if(i == 0 && !topright)
                    continue;
                  if(i == 4 && !bottomright)
                    continue;
                }
              varrt[i][j] = interpolate_quadratic_coeff(opal_table->cache->Ax1,
                                                        opal_table->cache->Bx1, opal_table->cache->Cx1, var[opal_table->cache->it + i][opal_table->cache->irho + j] + opal_table->cache->ix);
            }
        }
    }
  /* right border */
  else if(opal_table->cache->mx >= (NX - 2))
    {
      for(i = 0; i < 4; i++)
        {
          if(i == 0 && !top)
            continue;
          if(i == 4 && !bottom)
            continue;
          for(j = 0; j < 4; j++)
            {
              if(j == 0)
                {
                  if(!left)
                    continue;
                  if(i == 0 && !topleft)
                    continue;
                  if(i == 4 && !bottomleft)
                    continue;
                }
              if(j == 4)
                {
                  if(!right)
                    continue;
                  if(i == 0 && !topright)
                    continue;
                  if(i == 4 && !bottomright)
                    continue;
                }
              varrt[i][j] = interpolate_quadratic_coeff(opal_table->cache->Ax2,
                                                        opal_table->cache->Bx2, opal_table->cache->Cx2, var[opal_table->cache->it + i][opal_table->cache->irho + j] + opal_table->cache->ix + 1);
            }
        }
    }
  /* in the middle */
  else
    {
      for(i = 0; i < 4; i++)
        {
          if(i == 0 && !top)
            continue;
          if(i == 4 && !bottom)
            continue;
          for(j = 0; j < 4; j++)
            {
              if(j == 0)
                {
                  if(!left)
                    continue;
                  if(i == 0 && !topleft)
                    continue;
                  if(i == 4 && !bottomleft)
                    continue;
                }
              if(j == 4)
                {
                  if(!right)
                    continue;
                  if(i == 0 && !topright)
                    continue;
                  if(i == 4 && !bottomright)
                    continue;
                }
              varrt1[i][j] = interpolate_quadratic_coeff(opal_table->cache->Ax1,
                                                         opal_table->cache->Bx1, opal_table->cache->Cx1, var[opal_table->cache->it + i][opal_table->cache->irho + j] + opal_table->cache->ix);
              varrt2[i][j] = interpolate_quadratic_coeff(opal_table->cache->Ax2,
                                                         opal_table->cache->Bx2, opal_table->cache->Cx2, var[opal_table->cache->it + i][opal_table->cache->irho + j] + opal_table->cache->ix + 1);
              varrt[i][j] = opal_table->cache->slopex * varrt1[i][j] + (1.0 - opal_table->cache->slopex) * varrt2[i][j];
            }
        }
    }

  /* interpolate in rho */
  double vart1[4], vart2[4];
  /* left border */
  if(rhoposition == LEFT)
    {
      for(i = 0; i < 4; i++)
        {
          if(i == 0 && !topleft)
            continue;
          if(i == 4 && !bottomleft)
            continue;
          vart1[i] = interpolate_quadratic_coeff(opal_table->cache->Arho1, opal_table->cache->Brho1, opal_table->cache->Crho1, varrt[i]);
        }
    }
  /* right border */
  else if(rhoposition == RIGHT)
    {
      for(i = 0; i < 4; i++)
        {
          if(i == 0 && !topright)
            continue;
          if(i == 4 && !bottomright)
            continue;
          vart2[i] = interpolate_quadratic_coeff(opal_table->cache->Arho2, opal_table->cache->Brho2, opal_table->cache->Crho2, varrt[i] + 1);
        }
    }
  /* in the middle */
  else
    {
      for(i = 0; i < 4; i++)
        {
          if(i == 0 && !top)
            continue;
          if(i == 4 && !bottom)
            continue;
          vart1[i] = interpolate_quadratic_coeff(opal_table->cache->Arho1, opal_table->cache->Brho1, opal_table->cache->Crho1, varrt[i]);
          vart2[i] = interpolate_quadratic_coeff(opal_table->cache->Arho2, opal_table->cache->Brho2, opal_table->cache->Crho2, varrt[i] + 1);
        }
    }

  /* interpolate in T, depending on position, final overlapping again in rho */
  /* top left, bottom left, top right, bottom right (for rho left to right, T bottom to top) */
  double var1 = 0., var2 = 0., var3 = 0., var4 = 0.;
  double var12 = 0., var34 = 0.;        /* left, right */

  /* interpolations */
  if(topleft)
    var1 = interpolate_quadratic_coeff(opal_table->cache->At1, opal_table->cache->Bt1, opal_table->cache->Ct1, vart1);
  if(bottomleft)
    var2 = interpolate_quadratic_coeff(opal_table->cache->At2, opal_table->cache->Bt2, opal_table->cache->Ct2, vart1 + 1);

  if(topright)
    var3 = interpolate_quadratic_coeff(opal_table->cache->At1, opal_table->cache->Bt1, opal_table->cache->Ct1, vart2);
  if(bottomright)
    var4 = interpolate_quadratic_coeff(opal_table->cache->At2, opal_table->cache->Bt2, opal_table->cache->Ct2, vart2 + 1);

  /* overlap in temperature */
  if(topleft && bottomleft)
    var12 = opal_table->cache->slopet * var1 + (1.0 - opal_table->cache->slopet) * var2;
  else if(topleft)
    var12 = var1;
  else if(bottomleft)
    var12 = var2;

  if(topright && bottomright)
    var34 = opal_table->cache->slopet * var3 + (1.0 - opal_table->cache->slopet) * var4;
  else if(topright)
    var34 = var3;
  else if(bottomright)
    var34 = var4;

  /* overlap in density */
  if(left && right)
    *result = opal_table->cache->sloperho * var12 + (1.0 - opal_table->cache->sloperho) * var34;
  else if(left)
    *result = var12;
  else if(right)
    *result = var34;

#endif
  return 0;
}

/**
 * @brief linear interpolation
 *
 * @param x value to interpolate at
 * @param xs pointer to array x values for interpolation
 * @param ys pointer to array y values for interpolation
 *
 * @return interpolated value
 */
double interpolate_linear(double x, double *xs, double *ys)
{
  return ys[0] + (x - xs[0]) * (ys[1] - ys[0]) / (xs[1] - xs[0]);
}

/**
 * @brief linear interpolation using given slope
 *
 * @param slope given slope
 * @param ys pointer to array y values for interpolation
 *
 * @return interpolated value
 */
double interpolate_linear_slope(double slope, double *ys)
{
  return ys[0] + slope * (ys[1] - ys[0]);
}


/**
 * @brief quadratic interpolation
 *
 * @param x value to interpolate at
 * @param xs pointer to array x values for interpolation
 * @param ys pointer to array y values for interpolation
 *
 * @return interpolated value
 */
double interpolate_quadratic(double x, double *xs, double *ys)
{
  double A = (x - xs[1]) * (x - xs[2]) / (xs[0] - xs[1]) / (xs[0] - xs[2]);
  double B = (x - xs[0]) * (x - xs[2]) / (xs[1] - xs[0]) / (xs[1] - xs[2]);
  double C = (x - xs[0]) * (x - xs[1]) / (xs[2] - xs[0]) / (xs[2] - xs[1]);
  return A * ys[0] + B * ys[1] + C * ys[2];
}

/**
 * @brief quadratic interpolation using given coefficients
 *
 * @param A,B,C given coefficients
 * @param ys y values for interpolation
 *
 * @return interpolated value
 */
double interpolate_quadratic_coeff(double A, double B, double C, double *ys)
{
  return A * ys[0] + B * ys[1] + C * ys[2];
}

/**
 * @brief binary search in lists in arbitrary order
 *
 * @param x value to search for
 * @param list list to search in
 * @param len length of list
 * @param points number of points for interpolating polynomials
 *        (2 for linear, 3 for quadratic etc)
 * @param low lowest index in array with x centered in a subrange of points points
 * @param m index with list[index] <= x < list[index+1] (ascending)
 *  or with list[index] < x <= list[index+1] (descending)
 */
void binary_search(double x, double *list, int len, int points, int *low, int *m)
{
  int ilow = 0, ihigh = len - 1, imid;
  int ascend = (list[0] <= list[len - 1]);
  while(ilow < ihigh - 1)
    {
      /* imid = (ilow + ihigh) / 2; */
      imid = (ilow + ihigh) >> 1;
      if((list[imid] <= x) == ascend)
        ilow = imid;
      else
        ihigh = imid;
    }
  /* now ilow is index with list[index] <= x < list[index+1] (ascending)
   *  or with list[index] < x <= list[index+1] (descending) */
  *m = ilow;

  /* low: min 0, max len-points, else m-floor((points-2)/2) */
  *low = ilow - ((points - 2) >> 1);
  *low = fmax(0, fmin(len - points, *low));
}

/**
 * @brief linear search in lists in ascending order
 *
 * @param x value to search for
 * @param list list to search in
 * @param len length of list
 *
 * @return index with list[index] <= x < list[index+1]
 */
int linear_search_asc(double x, double *list, int len)
{
  int ilow = 0;
  while((list[ilow] > x) && (ilow < len))
    ilow++;
  if(ilow < len)
    return ilow;
  else
    return -1;
}

/**
 * @brief checks if EOS mode is consistent
 *
 * @param mode EOS mode
 *
 * @return 0 if ok, -1 if error
 */
int checkmode(int mode)
{
  /* for gamma3, the other gammas are needed */
  if ((mode & EOS_GAMMA3) && ( !(mode & EOS_GAMMA1) || !(mode & EOS_GAMMA2) ))
    mode = mode | EOS_GAMMA1 | EOS_GAMMA2;
  if(mode & EOS_RADIATION)
    {
      if ((mode & EOS_GAMMA2) && !(mode & EOS_GAMMA1)) mode = mode | EOS_GAMMA1;
      if ((mode & EOS_GAMMA1) && !(mode & EOS_GAMMA3)) mode = mode | EOS_GAMMA3;
      if ((mode & EOS_GAMMA3) && ( !(mode & EOS_CHIT) || !(mode & EOS_CHIR) ||
            !(mode & EOS_CV) ))
        mode = mode | EOS_CHIT | EOS_CHIR | EOS_CV;
      if ((mode & EOS_CHIT) && !(mode & EOS_PRESSURE)) mode = mode | EOS_PRESSURE;
      if ((mode & EOS_CHIR) && !(mode & EOS_PRESSURE)) mode = mode | EOS_PRESSURE;
    }
  return mode;
}


/**
 * @brief computes mean molecular weight of composition
 *
 * @param xh H mass fraction
 * @param z metallicity
 *
 * @return mean molecular weight
 */
double mean_molecular_weight(double xh, double z)
{
  /* atomic data: H, He, C, N, O, Ne */
  const double zi[6] = { 1.0, 2.0, 6.0, 7.0, 8.0, 10.0 };
                                                       /**< number of protons */
  const double mui[6] = { 1.0079, 4.0026, 12.011, 14.0067,
    15.9994, 20.179
  };                                                     /**< molecular weight */
  const double xi0[6] = { 0.0, 0.0, 0.247137766, 0.0620782,
    0.52837118, 0.1624188
  };                                     /**< initial mass fraction of metals */
  double muz; /**< molecular weight of metals */
  double xi[7]; /**< mass fractions */
  double mu; /**< mean molecular weight */
  int i;

  /* compute molecular weight of metals */
  for(muz = 0.0, i = 2; i < 6; i++)
    {
      muz += xi0[i] * mui[i];
    }

  /* compute mass fractions */
  xi[0] = xh;
  xi[1] = 1.0 - xh - z;
  for(i = 2; i < 6; i++)
    {
      xi[i] = xi0[i] * z * mui[i] / muz;
    }

  /* compute mean molecular weight */
  mu = 0.0;
  for(i = 0; i < 6; i++)
    {
      mu += (1.0 + zi[i]) * xi[i] / mui[i];
    }
  return 1.0/mu;
}

/**
 * @brief computes mean molecular weight of composition for neutral gas
 *
 * @param xh H mass fraction
 * @param z metallicity
 *
 * @return mean molecular weight for neutral gas
 */
double mean_molecular_weight_neutral(double xh, double z)
{
  /* atomic data: H, He, C, N, O, Ne */
  const double zi[6] = {1.0, 2.0, 6.0, 7.0, 8.0, 10.0};/**< number of protons */
  const double mui[6] = {1.0079, 4.0026, 12.011, 14.0067,
    15.9994, 20.179};                                    /**< molecular weight */
  const double xi0[6] = {0.0, 0.0, 0.247137766, 0.0620782,
    0.52837118, 0.1624188};              /**< initial mass fraction of metals */
  double muz; /**< molecular weight of metals */
  double xi[7]; /**< mass fractions */
  double mu; /**< mean molecular weight */
  int i;

  /* compute molecular weight of metals */
  for (muz = 0.0, i = 2; i < 6; i++)
    {
      muz += xi0[i] * mui[i];
    }

  /* compute mass fractions */
  xi[0] = xh;
  xi[1] = 1.0 - xh - z;
  for (i = 2; i < 6; i++)
    {
      xi[i] = xi0[i] * z * mui[i] / muz;
    }

  /* compute mean molecular weight */
  mu = 0.0;
  for (i = 0; i < 6; i++)
    {
      mu +=  xi[i] / mui[i];
    }
  return 1.0/mu;
}

/**
 * @brief computes constant for entropy computation
 *
 * @param xh H mass fraction
 * @param z metallicity
 * @param neutral either 1 for neutral gas or 0 for ionized gas
 *
 * @return constant to add to entropy
 */
double entropy_constant(double xh, double z, int neutral)
{
  /* atomic data: H, He, C, N, O, Ne */
  const double zi[6] = {1.0, 2.0, 6.0, 7.0, 8.0, 10.0};/**< number of protons */
  const double mui[6] = {1.0079, 4.0026, 12.011, 14.0067,
    15.9994, 20.179};                                    /**< molecular weight */
  const double xi0[6] = {0.0, 0.0, 0.247137766, 0.0620782,
    0.52837118, 0.1624188};              /**< initial mass fraction of metals */
  const double gii[6] = {2.0, 2.0, 1.0, 3.0, 1.0, 1.0};
                  /**< statistical weights of ionized atoms */
  const double gi0[6] = {1.0, 1.0, 1.0, 4.0, 5.0, 1.0};
                  /**< statistical weights of neutral atoms */
  double muz; /**< molecular weight of metals */
  double xi[7]; /**< mass fractions */
  double mue; /**< electron fraction */
  double mu; /**< mean molecular weight */
  double result; /**< final result */
  int i;

  /* compute molecular weight of metals */
  for (muz = 0.0, i = 2; i < 6; i++)
    {
      muz += xi0[i] * mui[i];
    }

  /* compute mass fractions */
  xi[0] = xh;
  xi[1] = 1.0 - xh - z;
  for (i = 2; i < 6; i++)
    {
      xi[i] = xi0[i] * z * mui[i] / muz;
    }

  /* compute mean molecular weight */
  mu = 0.0;
  if (neutral == 1)
    {
      for (i = 0; i < 6; i++)
          mu +=  xi[i] / mui[i];
    }
  else
    {
      for (i = 0; i < 6; i++)
          mu += (1.0 + zi[i]) * xi[i] / mui[i];
    }
  mu = 1.0/mu;

  /* first part of constant */
  result = 2.5 * GAS_CONSTANT / mu;
  if (neutral == 1)
    {
      /* contribution by ions */
      for (i = 0; i < 6; i++)
        {
          result += GAS_CONSTANT * xi[i]/mui[i] * log(
              pow(mui[i]*UAMU*BOLTZMANN_K/(2.*PI*HBAR*HBAR), 1.5)
              * gi0[i]*mui[i]*UAMU/xi[i]);
        }
      /* no contribution by electrons since gas is neutral */
    }
  else
    {
      /* contribution by ions */
      for (i = 0; i < 6; i++)
        {
          result += GAS_CONSTANT * xi[i]/mui[i] * log(
              pow(mui[i]*UAMU*BOLTZMANN_K/(2.*PI*HBAR*HBAR), 1.5)
              * gii[i]*mui[i]*UAMU/xi[i]);
        }
      /* contribution by electrons */
      mue = 0;
      for (i = 0; i < 6; i++)
          mue += zi[i]*xi[i]/mui[i];

      result += GAS_CONSTANT * mue * log(
              pow(MASS_E*BOLTZMANN_K/(2.*PI*HBAR*HBAR), 1.5)
              * 2.0*UAMU/mue);
    }

  return result;
}

/**
 * @brief Get temperature for which X(H I)=0.5
 *
 * @param rho density
 *
 * @return temperature from fit
 */
double get_temp_hion(double rho)
{
  /* values from broken power law fit to the line
   * X(H I)=0.5 in the rho-T plane
   */
  double a1 = 4.4555e4;
  double alpha1 = 6.4315e-2;
  double alpha2 = 0.64903;
  double rhobr = 3.8612e-4;

  if (rho < 1e-14) rho = 1e-14;

  if (rho < rhobr)
    return a1*pow(rho, alpha1);
  else
    return a1*pow(rhobr, alpha1)*pow(rho/rhobr, alpha2);
}

/**
 * @brief Updates values in cache for faster interpolation
 *
 * @param opal_table the opal table struct
 * @param x,rho,t EOS parameters
 *
 * @return -1 if error occurred, 0 otherwise
 */
int update_cache(struct opal_eos_table *opal_table, double x, double rho, double t)
{
  /* check if update needed */
  if (fabs(opal_table->cache->x   - x)   < 1e-15 &&
      fabs(opal_table->cache->rho - rho) < 1e-15 &&
      fabs(opal_table->cache->t   - t)   < 1e-15 &&
      opal_table->cache->valid)
    return 0;

  /* set to non valid first */
  opal_table->cache->valid = 0;

  int num_points = 0;
#if INTERPOLATION_METHOD == IP_LINEAR
  num_points = 2;
#elif INTERPOLATION_METHOD == IP_QUADRATIC
  num_points = 4;
#endif

  /* check boundaries and get indices if values are different */
  if(fabs(opal_table->cache->x - x) > 1e-15)
    {
      if((x < opal_table->minx) || (x > opal_table->maxx))
        {
#ifdef OPAL_DEBUG
          fprintf(stderr, "Error: X out of range. Value: %e, limits: %e to %e.\n", x, opal_table->minx, opal_table->maxx);
#endif
          return -1;
        }
      /* table is equidistant in X */
      binary_search(x, opal_table->x, NX, num_points, &opal_table->cache->ix, &opal_table->cache->mx);
      /*printf("ix: %d, x: %f %f %f\n", ix, x, opal_table->x[ix], opal_table->x[ix+1]); */
    }

  if(fabs(opal_table->cache->rho - rho) > 1e-15)
    {
      if((rho < opal_table->minrho) || (rho > opal_table->maxrho))
        {
#ifdef OPAL_DEBUG
          fprintf(stderr, "Error: rho out of range. Value: %e, limits: %e to %e.\n", rho, opal_table->minrho, opal_table->maxrho);
#endif
          return -1;
        }
      /* table is in ascending order for rho */
      binary_search(rho, opal_table->rho, NRHO, num_points, &opal_table->cache->irho, &opal_table->cache->mrho);
      /*printf("irho: %d, rho: %f %f %f\n", irho, rho, opal_table->rho[irho], opal_table->rho[irho+1]); */
    }

  /* check also if temperature limits changed with rho! */
  double mint = dmax(opal_table->mint[opal_table->cache->mrho],
                     opal_table->mint[opal_table->cache->mrho + 1]);
  //double mint = opal_table->mint[opal_table->cache->mrho];
  if((fabs(opal_table->cache->t - t) > 1e-15) || (fabs(opal_table->cache->mint - mint) > 1e-15))
    {
      if((t < mint) || (t > opal_table->maxt))
        {
#ifdef OPAL_DEBUG
          fprintf(stderr, "Error: T out of range. Value: %e, limits: %e to %e \
for x=%e, rho=%e.\n", t * 1e6, mint * 1e6, opal_table->maxt * 1e6, x, rho);
#endif
          opal_table->cache->mint = mint;
          return -1;
        }
      /* table is in descending order for T */
      int minnum = imin(opal_table->numt[opal_table->cache->mrho],
                        opal_table->numt[opal_table->cache->mrho + 1]);
      //int minnum = opal_table->numt[opal_table->cache->mrho];
      binary_search(t, opal_table->t, minnum, num_points, &opal_table->cache->it, &opal_table->cache->mt);
      if(fabs(t - opal_table->maxt * 1e6) < 1e-15)
        opal_table->cache->it = opal_table->cache->mt = 0;
    }

  /* pre-compute interpolation quantities */
#if INTERPOLATION_METHOD == IP_LINEAR
  /* set slopes */
  if(fabs(opal_table->cache->x - x) > 1e-15)
    opal_table->cache->slopex = (x - opal_table->x[opal_table->cache->ix]) / (opal_table->x[opal_table->cache->ix + 1] - opal_table->x[opal_table->cache->ix]);
  if(fabs(opal_table->cache->rho - rho) > 1e-15)
    opal_table->cache->sloperho = (rho - opal_table->rho[opal_table->cache->irho]) / (opal_table->rho[opal_table->cache->irho + 1] - opal_table->rho[opal_table->cache->irho]);
  if((fabs(opal_table->cache->t - t) > 1e-15) || (fabs(opal_table->cache->mint - mint) > 1e-15))
    opal_table->cache->slopet = (t - opal_table->t[opal_table->cache->it]) / (opal_table->t[opal_table->cache->it + 1] - opal_table->t[opal_table->cache->it]);
#elif INTERPOLATION_METHOD == IP_QUADRATIC
  /* update coefficients only if values change */
  if(fabs(opal_table->cache->x - x) > 1e-15)
    {
      /* set A, B and C coefficients */
      compute_quadratic_coeff(x, opal_table->x + opal_table->cache->ix, &opal_table->cache->Ax1, &opal_table->cache->Bx1, &opal_table->cache->Cx1);
      compute_quadratic_coeff(x, opal_table->x + opal_table->cache->ix + 1, &opal_table->cache->Ax2, &opal_table->cache->Bx2, &opal_table->cache->Cx2);
      /* set slopes for overlapping quadratics */
      opal_table->cache->slopex = (opal_table->x[opal_table->cache->ix + 2] - x) / (opal_table->x[opal_table->cache->ix + 2] - opal_table->x[opal_table->cache->ix + 1]);
    }
  if(fabs(opal_table->cache->rho - rho) > 1e-15)
    {
      /* set A, B and C coefficients */
      compute_quadratic_coeff(rho, opal_table->rho + opal_table->cache->irho, &opal_table->cache->Arho1, &opal_table->cache->Brho1, &opal_table->cache->Crho1);
      compute_quadratic_coeff(rho, opal_table->rho + opal_table->cache->irho + 1, &opal_table->cache->Arho2, &opal_table->cache->Brho2, &opal_table->cache->Crho2);
      /* set slopes for overlapping quadratics */
      opal_table->cache->sloperho = (opal_table->rho[opal_table->cache->irho + 2] - rho) / (opal_table->rho[opal_table->cache->irho + 2] - opal_table->rho[opal_table->cache->irho + 1]);
    }
  if((fabs(opal_table->cache->t - t) > 1e-15) || (fabs(opal_table->cache->mint - mint) > 1e-15))
    {
      /* set A, B and C coefficients */
      compute_quadratic_coeff(t, opal_table->t + opal_table->cache->it, &opal_table->cache->At1, &opal_table->cache->Bt1, &opal_table->cache->Ct1);
      compute_quadratic_coeff(t, opal_table->t + opal_table->cache->it + 1, &opal_table->cache->At2, &opal_table->cache->Bt2, &opal_table->cache->Ct2);
      /* set slopes for overlapping quadratics */
      opal_table->cache->slopet = (opal_table->t[opal_table->cache->it + 2] - t) / (opal_table->t[opal_table->cache->it + 2] - opal_table->t[opal_table->cache->it + 1]);
    }
#endif

  opal_table->cache->x = x;
  opal_table->cache->rho = rho;
  opal_table->cache->t = t;
  opal_table->cache->mint = mint;


  /* for debugging
     printf("irho: %d, mrho: %d; it: %d, mt: %d; numt[mrho]: %d, numt[mrho+1]: %d, numt[mrho+2]: %d;\n",
     opal_table->cache->irho, opal_table->cache->mrho,
     opal_table->cache->it, opal_table->cache->mt,
     opal_table->numt[opal_table->cache->mrho],
     opal_table->numt[opal_table->cache->mrho + 1],
     opal_table->numt[opal_table->cache->mrho + 2]);
     printf("t: %e, t[mt]: %e, t[mt+1]: %e, t[mt+2]: %e, mint[mrho+1]: %e\n",
     opal_table->cache->t, opal_table->t[opal_table->cache->mt],
     opal_table->t[opal_table->cache->mt + 1],
     opal_table->t[opal_table->cache->mt + 2],
     opal_table->mint[opal_table->cache->mrho + 1]);
   */

  /* everything fine, set to valid */
  opal_table->cache->valid = 1;
  return 0;
}

/**
 * @brief pre-compute coefficients for quadratic interpolation
 *
 * @param x value to interpolate at
 * @param xs table x values
 * @param A,B,C the coefficients are stored here
 */
void compute_quadratic_coeff(double x, double *xs, double *A, double *B, double *C)
{
  *A = (x - xs[1]) * (x - xs[2]) / (xs[0] - xs[1]) / (xs[0] - xs[2]);
  *B = (x - xs[0]) * (x - xs[2]) / (xs[1] - xs[0]) / (xs[1] - xs[2]);
  *C = (x - xs[0]) * (x - xs[1]) / (xs[2] - xs[0]) / (xs[2] - xs[1]);
}

/**
 * @brief check limits of OPAL EOS
 *
 * @param rho density to check for
 * @param t temperature to check for
 * @param F distance to border is returned
 *
 * @result 1 - in OPAL range; 0 - ideal gas only; -1 - in blending region;
 */
int check_limits(double rho, double t, double *F)
{
  const double blendfactor = 0.01; /* 5 per cent for blending region */
  double upper = 1.0 - blendfactor;
  double lower = 1.0 + blendfactor;
  const double rholow = 1e-14;
  const double rhohigh = 1e7;
  const double tlow = 1.87e3;
  const double thigh = 2e8;
  double touter, tinner;
  /* check if out of limits of OPAL EOS first */
  if (rho < rholow || rho > rhohigh || t < tlow || t > thigh)
    return 0; /* do ideal gas only */
  touter = outerborder(rho);
  if (t < touter)
    return 0;
  /* check if inside OPAL region */
  tinner = innerborder(rho);
  if (rho > lower*rholow && rho < upper*rhohigh && t > tinner && t < upper*thigh)
    return 1;
  /* now in blending region; compute blending factor F*/
  if (t < tinner)
    *F = get_blending_factor(rho, t);
  if (t > upper*thigh)
    *F = (t - upper*thigh)/(blendfactor*thigh);
  if (rho > upper*rhohigh)
    *F = (rho - upper*rhohigh)/(blendfactor*rhohigh);
  if (rho < lower*rholow)
    *F = (lower*rholow - rho)/(blendfactor*rholow);
  *F = dmax(0.0, dmin(1.0, *F));

  return -1;
}

/**
 * @brief straight line through (rho1,t1) and (rho2,t2) in log space
 *
 * @param rho density value to compute temperature for
 * @param rho1,t1 first point in plane
 * @param rho2,t2 secon point in plane
 *
 * @result returns the temperature corresponding to the density rho
 */
double line(double rho, double rho1, double t1, double rho2, double t2)
{
  return pow(10, log10(t1) + log10(t2/t1)/log10(rho2/rho1) * log10(rho/rho1));
}

/**
 * @brief distance to straight line through (rho1,t1) and (rho2,t2) in log space
 *
 * @param rho density value to compute distance for
 * @param t temperature value to compute distance for
 * @param rho1,t1 first point in plane
 * @param rho2,t2 secon point in plane
 *
 * @result returns the distance
 */
double linedist(double rho, double t, double rho1, double t1, double rho2, double t2)
{
  double m, b;
  double rhol, tl;
  m = log10(t2/t1)/log10(rho2/rho1);
  b = log10(t1) - log10(t2/t1)/log10(rho2/rho1) * log10(rho1);
  rhol = (rho + m*(t - b))/(1.0 + m*m);
  tl = m*rhol + b;
  return sqrt((rhol - rho)*(rhol - rho) + (tl - t)*(tl - t));
}

/**
 * @brief outer border of OPAL EOS; end of blending region
 *
 * @param rho density
 *
 * @return temperature of outer border
 */
double outerborder(double rho)
{
  if (rho < rho1)
    return t1;
  else if (rho < rho2)
    return line(rho, rho1, t1, rho2, t2);
  else if (rho < rho3)
    return t2;
  else if (rho < rho4)
    return line(rho, rho3, t2, rho4, t3);
  else
    return t3;
}

/**
 * @brief inner border of OPAL EOS; start of blending region
 *
 * @param rho density
 *
 * @return temperature of inner border
 */
double innerborder(double rho)
{
  if (rho < rho1_i)
    return t1_i;
  else if (rho < rho2_i)
    return line(rho, rho1_i, t1_i, rho2_i, t2_i);
  else if (rho < rho3_i)
    return t2_i;
  else if (rho < rho4_i)
    return line(rho, rho3_i, t2_i, rho4_i, t3_i);
  else
    return t3_i;
}

/**
 * @brief compute blending factor in transition zone between OPAL and ideal EOS
 *
 * @param rho density
 * @param t temperature
 *
 * @result blending factor F, between 0 and 1
 */
double get_blending_factor(double rho, double t)
{
  double d1, d2;
  if (rho < rho1)
    {
      d1 = t1_i - t;
      d2 = t - t1;
    }
  else if (rho < rho2)
    {
      d1 = linedist(rho, t, rho1_i, t1_i, rho2_i, t2_i);
      d2 = linedist(rho, t, rho1, t1, rho2, t2);
    }
  else if (rho < rho3)
    {
      d1 = t2_i - t;
      d2 = t - t2;
    }
  else if (rho < rho4)
    {
      d1 = linedist(rho, t, rho3_i, t2_i, rho4_i, t3_i);
      d2 = linedist(rho, t, rho3, t2, rho4, t3);
    }
  else
    {
      d1 = t3_i - t;
      d2 = t - t3;
    }
  return d1 / (d1 + d2);
}

/**
 * @brief blend results of OPAL and ideal EOS
 *
 * @param result first result structure
 * @param result_tmp second result structure
 * @param F blending parameter, between 0 and 1
 */
void blend_results(struct opal_eos_result *result, struct opal_eos_result *result_tmp,
    double F)
{
  double S;
  S = 0.5*(1-cos(F*PI));
  result->p = result->p * S + (1.0 - S) * result_tmp->p;
  //printf("F: %e, S: %e, p: %e %e\n", F, S, result->p, result_tmp->p);
  result->e = result->e * S + (1.0 - S) * result_tmp->e;
  result->s = result->s * S + (1.0 - S) * result_tmp->s;
  result->dedrho = result->dedrho * S + (1.0 - S) * result_tmp->dedrho;
  result->cv = result->cv * S + (1.0 - S) * result_tmp->cv;
  result->chir = result->chir * S + (1.0 - S) * result_tmp->chir;
  result->chit = result->chit * S + (1.0 - S) * result_tmp->chit;
  result->gamma1 = result->gamma1 * S + (1.0 - S) * result_tmp->gamma1;
  result->gamma2 = result->gamma2 * S + (1.0 - S) * result_tmp->gamma2;
  result->gamma3 = result->gamma3 * S + (1.0 - S) * result_tmp->gamma3;
}
#endif /* EOS_OPAL */
