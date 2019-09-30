/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/forcetree_ewald.c
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
#include <time.h>

#include "allvars.h"
#include "proto.h"

/*! \file forcetree_ewald.c
 *  \brief code for Ewald correction
 *
 *  This file contains the computation of the Ewald correction table
 */
#if !defined(PMGRID) && defined(SELFGRAVITY) && !defined(GRAVITY_NOT_PERIODIC) && defined(PERIODIC) && !defined(ONEDIMS_SPHERICAL)

#include <gsl/gsl_sf_bessel.h>

/* variables for Ewald correction lookup table */
MyFloat Ewd_fcorrx[ENX + 1][ENY + 1][ENZ + 1];
MyFloat Ewd_fcorry[ENX + 1][ENY + 1][ENZ + 1];
MyFloat Ewd_fcorrz[ENX + 1][ENY + 1][ENZ + 1];
MyFloat Ewd_potcorr[ENX + 1][ENY + 1][ENZ + 1];
double Ewd_fac_intp;

typedef struct
{
  int resx, resy, resz, varsize, ewaldtype;
} ewald_header;


/*! \brief This function initializes tables with the correction force and the
 *  correction potential due to the periodic images of a point mass located
 *  at the origin.
 *
 *  These corrections are obtained by Ewald summation. (See for example
 *  Hernquist, Bouchet, Suto, ApJS, 1991, 75, 231) The correction fields
 *  are used to obtain the full periodic force if periodic boundaries
 *  combined with the pure tree algorithm are used. For the TreePM
 *  algorithm, the Ewald correction is not used.
 *
 *  The correction terms are computed by ewald_psi() and ewald_force() and stored in
 *  the arrays Ewd_fcorrx, Ewd_fcorry, Ewd_fcorrz and Ewd_potcorr.
 *
 *  The correction fields are stored on disk once they are computed. If a
 *  corresponding file is found, they are loaded from disk to speed up the
 *  initialization. The Ewald summation issrc/gravtree_forcetest.c done in parallel, i.e. the
 *  processors share the work to compute the tables if needed.
 */
void ewald_init(void)
{
  int recomputeflag = 0;
  double force[3];
  char buf[200];
  FILE *fd;

  mpi_printf("EWALD: initialize Ewald correction...\n");

#ifdef LONG_X
  if(LONG_X != (int) (LONG_X))
    terminate("LONG_X must be an integer");
#endif

#ifdef LONG_Y
  if(LONG_Y != (int) (LONG_Y))
    terminate("LONG_Y must be an integer");
#endif

#ifdef LONG_Z
  if(LONG_Z != (int) (LONG_Z))
    terminate("LONG_Z must be an integer");
#endif

  sprintf(buf, "ewald_table_%d_%d_%d.dat", ENX, ENY, ENZ);

  if(ThisTask == 0)
    {
      if((fd = fopen(buf, "r")))
        {
          mpi_printf("\nEWALD: reading Ewald tables from file `%s'\n", buf);

          ewald_header tabh;
          my_fread(&tabh, sizeof(ewald_header), 1, fd);

#ifndef GRAVITY_TALLBOX
          int ewaldtype = -1;
#else
          int ewaldtype = GRAVITY_TALLBOX + 1;
#endif
          if(tabh.resx != ENX || tabh.resy != ENY || tabh.resz != ENZ || tabh.varsize != sizeof(MyFloat) || tabh.ewaldtype != ewaldtype)
            {
              mpi_printf("\nEWALD: something's wrong with this table file. Discarding it.\n");
              recomputeflag = 1;
            }
          else
            {
              my_fread(Ewd_fcorrx, sizeof(MyFloat), (ENX + 1) * (ENY + 1) * (ENZ + 1), fd);
              my_fread(Ewd_fcorry, sizeof(MyFloat), (ENX + 1) * (ENY + 1) * (ENZ + 1), fd);
              my_fread(Ewd_fcorrz, sizeof(MyFloat), (ENX + 1) * (ENY + 1) * (ENZ + 1), fd);
              my_fread(Ewd_potcorr, sizeof(MyFloat), (ENX + 1) * (ENY + 1) * (ENZ + 1), fd);

              recomputeflag = 0;
            }
          fclose(fd);
        }
      else
        recomputeflag = 1;
    }

  MPI_Bcast(&recomputeflag, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if(recomputeflag)
    {
      mpi_printf("\nEWALD: No usable Ewald tables in file `%s' found. Recomputing them...\n", buf);

      /* ok, let's recompute things. Actually, we do that in parallel. */

      int size = (ENX + 1) * (ENY + 1) * (ENZ + 1);
      int first, count;

      subdivide_evenly(size, NTask, ThisTask, &first, &count);

      for(int n = first; n < first + count; n++)
        {
          // n = i * (ENY+1)*(ENZ + 1) + j * (ENZ + 1) + k;

          int i = n / ((ENY + 1) * (ENZ + 1));
          int j = (n - i * (ENY + 1) * (ENZ + 1)) / (ENZ + 1);
          int k = (n - i * (ENY + 1) * (ENZ + 1) - j * (ENZ + 1));


          if(ThisTask == 0)
            {
              if(((n - first) % (count / 20)) == 0)
                {
                  printf("%4.1f percent done\n", (n - first) / (count / 100.0));
                  myflush(stdout);
                }
            }

          double xx = 0.5 * DBX * STRETCHX * ((double) i) / ENX;
          double yy = 0.5 * DBY * STRETCHY * ((double) j) / ENY;
          double zz = 0.5 * DBZ * STRETCHZ * ((double) k) / ENZ;

#ifndef GRAVITY_TALLBOX
          Ewd_potcorr[i][j][k] = ewald_psi(xx, yy, zz);

          ewald_force(xx, yy, zz, force);

          Ewd_fcorrx[i][j][k] = force[0];
          Ewd_fcorry[i][j][k] = force[1];
          Ewd_fcorrz[i][j][k] = force[2];
#else
          switch (GRAVITY_TALLBOX)
            {
            case 0:
              Ewd_potcorr[i][j][k] = ewald_psi(yy, zz, xx);

              ewald_force(yy, zz, xx, force);

              Ewd_fcorrx[i][j][k] = force[2];
              Ewd_fcorry[i][j][k] = force[0];
              Ewd_fcorrz[i][j][k] = force[1];
              break;
            case 1:
              Ewd_potcorr[i][j][k] = ewald_psi(xx, zz, yy);

              ewald_force(xx, zz, yy, force);

              Ewd_fcorrx[i][j][k] = force[0];
              Ewd_fcorry[i][j][k] = force[2];
              Ewd_fcorrz[i][j][k] = force[1];
              break;
            case 2:
              Ewd_potcorr[i][j][k] = ewald_psi(xx, yy, zz);

              ewald_force(xx, yy, zz, force);

              Ewd_fcorrx[i][j][k] = force[0];
              Ewd_fcorry[i][j][k] = force[1];
              Ewd_fcorrz[i][j][k] = force[2];
              break;
            }
#endif
        }

      int *recvcnts = (int *) mymalloc("recvcnts", NTask * sizeof(int));
      int *recvoffs = (int *) mymalloc("recvoffs", NTask * sizeof(int));

      for(int i = 0; i < NTask; i++)
        {
          int off, cnt;
          subdivide_evenly(size, NTask, i, &off, &cnt);
          recvcnts[i] = cnt * sizeof(MyFloat);
          recvoffs[i] = off * sizeof(MyFloat);
        }

      MPI_Allgatherv(MPI_IN_PLACE, size * sizeof(MyFloat), MPI_BYTE, Ewd_fcorrx, recvcnts, recvoffs, MPI_BYTE, MPI_COMM_WORLD);
      MPI_Allgatherv(MPI_IN_PLACE, size * sizeof(MyFloat), MPI_BYTE, Ewd_fcorry, recvcnts, recvoffs, MPI_BYTE, MPI_COMM_WORLD);
      MPI_Allgatherv(MPI_IN_PLACE, size * sizeof(MyFloat), MPI_BYTE, Ewd_fcorrz, recvcnts, recvoffs, MPI_BYTE, MPI_COMM_WORLD);
      MPI_Allgatherv(MPI_IN_PLACE, size * sizeof(MyFloat), MPI_BYTE, Ewd_potcorr, recvcnts, recvoffs, MPI_BYTE, MPI_COMM_WORLD);

      myfree(recvoffs);
      myfree(recvcnts);

      mpi_printf("\nEWALD: writing Ewald tables to file `%s'\n", buf);
      if(ThisTask == 0)
        {
          if((fd = fopen(buf, "w")))
            {
              ewald_header tabh;
              tabh.resx = ENX;
              tabh.resy = ENY;
              tabh.resz = ENZ;
              tabh.varsize = sizeof(MyFloat);
#ifndef GRAVITY_TALLBOX
              tabh.ewaldtype = -1;
#else
              tabh.ewaldtype = GRAVITY_TALLBOX + 1;
#endif
              my_fwrite(&tabh, sizeof(ewald_header), 1, fd);

              my_fwrite(Ewd_fcorrx, sizeof(MyFloat), (ENX + 1) * (ENY + 1) * (ENZ + 1), fd);
              my_fwrite(Ewd_fcorry, sizeof(MyFloat), (ENX + 1) * (ENY + 1) * (ENZ + 1), fd);
              my_fwrite(Ewd_fcorrz, sizeof(MyFloat), (ENX + 1) * (ENY + 1) * (ENZ + 1), fd);
              my_fwrite(Ewd_potcorr, sizeof(MyFloat), (ENX + 1) * (ENY + 1) * (ENZ + 1), fd);
              fclose(fd);
            }
        }
    }
  else
    {
      /* here we got them from disk */
      int len = (ENX + 1) * (ENY + 1) * (ENZ + 1) * sizeof(MyFloat);

      MPI_Bcast(Ewd_fcorrx, len, MPI_BYTE, 0, MPI_COMM_WORLD);
      MPI_Bcast(Ewd_fcorry, len, MPI_BYTE, 0, MPI_COMM_WORLD);
      MPI_Bcast(Ewd_fcorrz, len, MPI_BYTE, 0, MPI_COMM_WORLD);
      MPI_Bcast(Ewd_potcorr, len, MPI_BYTE, 0, MPI_COMM_WORLD);
    }


  /* now scale things to the boxsize that is actually used */

  Ewd_fac_intp = 2 * EN / All.BoxSize;

  for(int i = 0; i <= ENX; i++)
    for(int j = 0; j <= ENY; j++)
      for(int k = 0; k <= ENZ; k++)
        {
          Ewd_potcorr[i][j][k] /= All.BoxSize;
          Ewd_fcorrx[i][j][k] /= All.BoxSize * All.BoxSize;
          Ewd_fcorry[i][j][k] /= All.BoxSize * All.BoxSize;
          Ewd_fcorrz[i][j][k] /= All.BoxSize * All.BoxSize;
        }

  mpi_printf("EWALD: Initialization of periodic boundaries finished.\n");
}


/*! \brief This function looks up the correction force due to the infinite number
 *  of periodic particle/node images.
 *
 *  We here use trilinear interpolation
 *  to get it from the precomputed tables, which contain one octant
 *  around the target particle at the origin. The other octants are
 *  obtained from it by exploiting the symmetry properties.
 *
 *  \param dx x component of the distance between the two particles
 *  \param dx y component of the distance between the two particles
 *  \param dx z component of the distance between the two particles
 *  \param fper pointer to array containing the correction force
 */
void ewald_corr(double dx, double dy, double dz, double *fper)
{
  int signx, signy, signz;
  int i, j, k;
  double u, v, w;
  double f1, f2, f3, f4, f5, f6, f7, f8;

  if(dx < 0)
    {
      dx = -dx;
      signx = +1;
    }
  else
    signx = -1;
  if(dy < 0)
    {
      dy = -dy;
      signy = +1;
    }
  else
    signy = -1;
  if(dz < 0)
    {
      dz = -dz;
      signz = +1;
    }
  else
    signz = -1;
  u = dx * Ewd_fac_intp;
  i = (int) u;
  if(i >= ENX)
    i = ENX - 1;
  u -= i;
  v = dy * Ewd_fac_intp;
  j = (int) v;
  if(j >= ENY)
    j = ENY - 1;
  v -= j;
  w = dz * Ewd_fac_intp;
  k = (int) w;
  if(k >= ENZ)
    k = ENZ - 1;
  w -= k;
  f1 = (1 - u) * (1 - v) * (1 - w);
  f2 = (1 - u) * (1 - v) * (w);
  f3 = (1 - u) * (v) * (1 - w);
  f4 = (1 - u) * (v) * (w);
  f5 = (u) * (1 - v) * (1 - w);
  f6 = (u) * (1 - v) * (w);
  f7 = (u) * (v) * (1 - w);
  f8 = (u) * (v) * (w);
  fper[0] = signx * (Ewd_fcorrx[i][j][k] * f1 +
                     Ewd_fcorrx[i][j][k + 1] * f2 +
                     Ewd_fcorrx[i][j + 1][k] * f3 +
                     Ewd_fcorrx[i][j + 1][k + 1] * f4 + Ewd_fcorrx[i + 1][j][k] * f5 + Ewd_fcorrx[i + 1][j][k + 1] * f6 + Ewd_fcorrx[i + 1][j + 1][k] * f7 + Ewd_fcorrx[i + 1][j + 1][k + 1] * f8);
  fper[1] =
    signy * (Ewd_fcorry[i][j][k] * f1 + Ewd_fcorry[i][j][k + 1] * f2 + Ewd_fcorry[i][j + 1][k] * f3 + Ewd_fcorry[i][j + 1][k + 1] * f4 +
             Ewd_fcorry[i + 1][j][k] * f5 + Ewd_fcorry[i + 1][j][k + 1] * f6 + Ewd_fcorry[i + 1][j + 1][k] * f7 + Ewd_fcorry[i + 1][j + 1][k + 1] * f8);
  fper[2] =
    signz * (Ewd_fcorrz[i][j][k] * f1 + Ewd_fcorrz[i][j][k + 1] * f2 + Ewd_fcorrz[i][j + 1][k] * f3 + Ewd_fcorrz[i][j + 1][k + 1] * f4 +
             Ewd_fcorrz[i + 1][j][k] * f5 + Ewd_fcorrz[i + 1][j][k + 1] * f6 + Ewd_fcorrz[i + 1][j + 1][k] * f7 + Ewd_fcorrz[i + 1][j + 1][k + 1] * f8);
}



/*! \brief This function looks up the correction potential due to the infinite
 *  number of periodic particle/node images.
 *
 *  We here use tri-linear interpolation to get it from the precomputed
 *  table, which contains one octant around the target particle at the
 *  origin. The other octants are obtained from it by exploiting symmetry
 *  properties.
 *
 *  \param dx x component of the distance between the two particles
 *  \param dx y component of the distance between the two particles
 *  \param dx z component of the distance between the two particles
 *  \return the correction potential
 */
double ewald_pot_corr(double dx, double dy, double dz)
{
  int i, j, k;
  double u, v, w;
  double f1, f2, f3, f4, f5, f6, f7, f8;

  if(dx < 0)
    dx = -dx;
  if(dy < 0)
    dy = -dy;
  if(dz < 0)
    dz = -dz;
  u = dx * Ewd_fac_intp;
  i = (int) u;
  if(i >= ENX)
    i = ENX - 1;
  u -= i;
  v = dy * Ewd_fac_intp;
  j = (int) v;
  if(j >= ENY)
    j = ENY - 1;
  v -= j;
  w = dz * Ewd_fac_intp;
  k = (int) w;
  if(k >= ENZ)
    k = ENZ - 1;
  w -= k;
  f1 = (1 - u) * (1 - v) * (1 - w);
  f2 = (1 - u) * (1 - v) * (w);
  f3 = (1 - u) * (v) * (1 - w);
  f4 = (1 - u) * (v) * (w);
  f5 = (u) * (1 - v) * (1 - w);
  f6 = (u) * (1 - v) * (w);
  f7 = (u) * (v) * (1 - w);
  f8 = (u) * (v) * (w);
  return Ewd_potcorr[i][j][k] * f1 +
    Ewd_potcorr[i][j][k + 1] * f2 +
    Ewd_potcorr[i][j + 1][k] * f3 +
    Ewd_potcorr[i][j + 1][k + 1] * f4 + Ewd_potcorr[i + 1][j][k] * f5 + Ewd_potcorr[i + 1][j][k + 1] * f6 + Ewd_potcorr[i + 1][j + 1][k] * f7 + Ewd_potcorr[i + 1][j + 1][k + 1] * f8;
}


#ifdef GRAVITY_TALLBOX

typedef struct
{
  double R;
  double rs;
  double z;
} param_data;

double bessel_int(double k, void *par)
{

  param_data *p = (param_data *) par;

  double R = p->R;
  double rs = p->rs;
  double z = p->z;

  return exp(-k * z) * gsl_sf_bessel_J0(k * R) * exp(-k * k * rs * rs);
}

#endif


/*! \brief This function computes the potential correction term by means of Ewald
 *  summation.
 *
 *  \param x, y, z contains the distance vector for which the correction
 *  term should be computed
 *  \return the correction term
 */
double ewald_psi(double x, double y, double z)
{
  static int printed = 0;

  double r = sqrt(x * x + y * y + z * z);

  if(r == 0)
    return 0;

#ifndef GRAVITY_TALLBOX

  double lmin = imin(imin(STRETCHX, STRETCHY), STRETCHZ);
  double alpha = 3.0 / lmin;

  const int nmax = 4;

  double sum1 = 0;
  for(int nx = -nmax; nx <= nmax; nx++)
    for(int ny = -nmax; ny <= nmax; ny++)
      for(int nz = -nmax; nz <= nmax; nz++)
        {
          double dx = x - nx * STRETCHX;
          double dy = y - ny * STRETCHY;
          double dz = z - nz * STRETCHZ;
          double r = sqrt(dx * dx + dy * dy + dz * dz);
          sum1 += erfc(alpha * r) / r;
        }

  double alpha2 = alpha * alpha;

  int nxmax = (int) (2 * alpha * (STRETCHX / lmin) + 0.5);
  int nymax = (int) (2 * alpha * (STRETCHY / lmin) + 0.5);
  int nzmax = (int) (2 * alpha * (STRETCHZ / lmin) + 0.5);

  if(printed == 0)
    {
      mpi_printf("EWALD: potential tab: nxmax=%d nymax=%d nzmax=%d\n", nxmax, nymax, nzmax);
      printed = 1;
    }

  double sum2 = 0.0;


  for(int nx = -nxmax; nx <= nxmax; nx++)
    for(int ny = -nymax; ny <= nymax; ny++)
      for(int nz = -nzmax; nz <= nzmax; nz++)
        {
          double kx = (2.0 * M_PI / (STRETCHX)) * nx;
          double ky = (2.0 * M_PI / (STRETCHY)) * ny;
          double kz = (2.0 * M_PI / (STRETCHZ)) * nz;
          double k2 = kx * kx + ky * ky + kz * kz;
          if(k2 > 0)
            {
              double kdotx = (x * kx + y * ky + z * kz);
              sum2 += 4.0 * M_PI / (k2 * STRETCHX * STRETCHY * STRETCHZ) * exp(-k2 / (4.0 * alpha2)) * cos(kdotx);
            }
        }
  
  double psi = /*-2.83729 + */  M_PI / (alpha * alpha * STRETCHX * STRETCHY * STRETCHZ) - sum1 - sum2 + 1.0 / r;

#else
  /* in the tallbox case, the third dimension, z, is assumed to be the non-periodic one */

  double lmin = imin(BOXX, BOXY);
  double alpha = 2.0 / lmin;

  double sum1 = 0.0;

  const int nmax = 4;

  for(int nx = -nmax; nx <= nmax; nx++)
    for(int ny = -nmax; ny <= nmax; ny++)
      {
        double dx = x - nx * BOXX;
        double dy = y - ny * BOXY;
        double r = sqrt(dx * dx + dy * dy + z * z);
        sum1 += erfc(alpha * r) / r;
      }

  int nxmax = (int) (2 * alpha * (BOXX / lmin) + 0.5);
  int nymax = (int) (2 * alpha * (BOXY / lmin) + 0.5);

  if(printed == 0)
    {
      mpi_printf("EWALD: potential tab: nxmax=%d nymax=%d\n", nxmax, nymax);
      printed = 1;
    }

  double alpha2 = alpha * alpha;

  double sum2 = 0.0;

  for(int nx = -nxmax; nx <= nxmax; nx++)
    for(int ny = -nymax; ny <= nymax; ny++)
      {
        if(nx != 0 || ny != 0)
          {
            double kx = (2.0 * M_PI / BOXX) * nx;
            double ky = (2.0 * M_PI / BOXY) * ny;
            double k2 = kx * kx + ky * ky;
            double k = sqrt(k2);

            sum2 += cos(kx * x + ky * y) * (exp(k * z) * erfc(k / (2 * alpha) + alpha * z) + exp(-k * z) * erfc(k / (2 * alpha) - alpha * z)) / k;
          }
      }

  sum2 *= M_PI / (BOXX * BOXY);

  r = sqrt(x * x + y * y + z * z);

  double psi = 2.0 * alpha / sqrt(M_PI) + 2 * sqrt(M_PI) / (BOXX * BOXY) * (exp(-alpha2 * z * z) / alpha + sqrt(M_PI) * z * erf(alpha * z)) - sum1 - sum2 + 1.0 / r;

#endif

  return psi;
}


/*! \brief This function computes the force correction term (difference between full
 *  force of infinite lattice and nearest image) by Ewald summation.
 *
 *  \param x, y, z contains the distance vector for which the correction
 *  force should be computed
 *  \param force  array will containing the correction force
 */
void ewald_force(double x, double y, double z, double force[3])
{
  static int printed = 0;
  for(int i = 0; i < 3; i++)
    force[i] = 0;
  double r2 = x * x + y * y + z * z;

  if(r2 == 0)
    return;

#ifndef GRAVITY_TALLBOX

  double lmin = imin(imin(STRETCHX, STRETCHY), STRETCHZ);
  double alpha = 2.0 / lmin;
  double alpha2 = alpha * alpha;

  double r3inv = 1.0 / (r2 * sqrt(r2));

  force[0] += r3inv * x;
  force[1] += r3inv * y;
  force[2] += r3inv * z;

  const int nmax = 4;

  for(int nx = -nmax; nx <= nmax; nx++)
    for(int ny = -nmax; ny <= nmax; ny++)
      for(int nz = -nmax; nz <= nmax; nz++)
        {
          double dx = x - nx * STRETCHX;
          double dy = y - ny * STRETCHY;
          double dz = z - nz * STRETCHZ;
          double r2 = dx * dx + dy * dy + dz * dz;
          double r = sqrt(r2);
          double val = erfc(alpha * r) + 2.0 * alpha * r / sqrt(M_PI) * exp(-alpha2 * r2);
          double val2 = val / (r2 * r);

          force[0] -= dx * val2;
          force[1] -= dy * val2;
          force[2] -= dz * val2;
        }


  int nxmax = (int) (2 * alpha * (STRETCHX / lmin) + 0.5);
  int nymax = (int) (2 * alpha * (STRETCHY / lmin) + 0.5);
  int nzmax = (int) (2 * alpha * (STRETCHZ / lmin) + 0.5);

  if(printed == 0)
    {
      mpi_printf("EWALD: force tab: nxmax=%d nymax=%d nzmax=%d\n", nxmax, nymax, nzmax);
      printed = 1;
    }

  for(int hx = -nxmax; hx <= nxmax; hx++)
    for(int hy = -nymax; hy <= nymax; hy++)
      for(int hz = -nzmax; hz <= nzmax; hz++)
        {
          double h2 = hx * hx + hy * hy + hz * hz;
          if(h2 > 0)
            {
              double hdotx = x * hx + y * hy + z * hz;
              double val = 2.0 / h2 * exp(-M_PI * M_PI * h2 / alpha2) * sin(2.0 * M_PI * hdotx);

              force[0] -= hx * val;
              force[1] -= hy * val;
              force[2] -= hz * val;
            }
        }
#else
  /* this is the case with periodicity only in two dimensions */

  double lmin = imin(BOXX, BOXY);
  double alpha = 2.0 / lmin;
  double alpha2 = alpha * alpha;

  double r3inv = 1.0 / (r2 * sqrt(r2));

  force[0] += r3inv * x;
  force[1] += r3inv * y;
  force[2] += r3inv * z;

  const int nmax = 4;

  for(int nx = -nmax; nx <= nmax; nx++)
    for(int ny = -nmax; ny <= nmax; ny++)
      {
        double dx = x - nx * BOXX;
        double dy = y - ny * BOXY;
        double dz = z;
        double r2 = dx * dx + dy * dy + dz * dz;
        double r = sqrt(r2);
        double val = erfc(alpha * r) + 2.0 * alpha * r / sqrt(M_PI) * exp(-alpha2 * r2);
        double val2 = val / (r2 * r);

        force[0] -= dx * val2;
        force[1] -= dy * val2;
        force[2] -= dz * val2;
      }

  int nxmax = (int) (2 * alpha * (BOXX / lmin) + 0.5);
  int nymax = (int) (2 * alpha * (BOXY / lmin) + 0.5);

  if(printed == 0)
    {
      mpi_printf("EWALD: force tab: nxmax=%d nymax=%d\n", nxmax, nymax);
      printed = 1;
    }

  for(int nx = -nxmax; nx <= nxmax; nx++)
    for(int ny = -nymax; ny <= nymax; ny++)
      {
        if(nx != 0 || ny != 0)
          {
            double kx = (2.0 * M_PI / BOXX) * nx;
            double ky = (2.0 * M_PI / BOXY) * ny;
            double k2 = kx * kx + ky * ky;
            double k = sqrt(k2);

            double val = M_PI / (BOXX * BOXY) * sin(kx * x + ky * y) * (exp(k * z) * erfc(k / (2 * alpha) + alpha * z) + exp(-k * z) * erfc(k / (2 * alpha) - alpha * z)) / k;

            force[0] += -kx * val;
            force[1] += -ky * val;
            force[2] += M_PI / (BOXX * BOXY) * cos(kx * x + ky * y) * (exp(k * z) * erfc(k / (2 * alpha) + alpha * z) - exp(-k * z) * erfc(k / (2 * alpha) - alpha * z));
          }
      }

  force[2] += -2.0 * M_PI / (BOXX * BOXY) * erf(alpha * z);
#endif
}

#endif
