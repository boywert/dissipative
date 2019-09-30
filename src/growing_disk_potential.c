/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/growing_disk_potential.c
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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stddef.h>

#include "allvars.h"
#include "proto.h"


#ifdef GROWING_DISK_POTENTIAL

static double comp_Dphi_z_disk_exact(double RR, double zz);
static double comp_Dphi_R_disk_exact(double RR, double zz);
static double intz_di(double k);
static double intz_di_abs(double);
static double intR_di(double k);
static double intR_di_abs(double k);
static double bessj0(double x);
static double bessj1(double x);
static double qromb(double (*func) (double), double a, double b);
static double *vector(long nl, long nh);
static void free_vector(double *v, long nl, long nh);



static double Rmin, Rmax, H, M_DISK;
static int Nbin, init_flag = 0;
static double *force_R, *force_z;
static double R, z;



void get_disk_forces(double RR, double zz, double *f_R, double *f_z)
{
  int i, j;
  double u, v;

  if(RR <= Rmin)
    {
      i = 0;
      u = 0;
    }
  else if(RR >= Rmax)
    {
      i = Nbin - 1;
      u = 1;
    }
  else
    {
      u = log(RR / Rmin) / log(Rmax / Rmin) * Nbin;
      i = (int) (u);
      u -= i;
    }

  if(zz <= Rmin)
    {
      j = 0;
      v = 0;
    }
  else if(zz >= Rmax)
    {
      j = Nbin - 1;
      v = 1;
    }
  else
    {
      v = log(zz / Rmin) / log(Rmax / Rmin) * Nbin;
      j = (int) (v);
      v -= j;
    }

  *f_R =
    tanh(RR / 0.01) * ((force_R[j * (Nbin + 1) + i]) * (1 - u) * (1 - v) +
                       (force_R[j * (Nbin + 1) + (i + 1)]) * (u) * (1 - v) + (force_R[(j + 1) * (Nbin + 1) + i]) * (1 - u) * (v) + (force_R[(j + 1) * (Nbin + 1) + (i + 1)]) * (u) * (v));

  *f_z =
    tanh(zz / 0.01) * ((force_z[j * (Nbin + 1) + i]) * (1 - u) * (1 - v) +
                       (force_z[j * (Nbin + 1) + (i + 1)]) * (u) * (1 - v) + (force_z[(j + 1) * (Nbin + 1) + i]) * (1 - u) * (v) + (force_z[(j + 1) * (Nbin + 1) + (i + 1)]) * (u) * (v));
}



double get_disk_mass(double time)
{
  return 5.8 * sqrt(time / 2.0);
}


void growing_disk_init(void)
{
  int i, j;

  if(init_flag)
    return;

  if(ThisTask == 0)
    printf("initializing disk grid\n");

  Nbin = 120;
  Rmin = 0.01;
  Rmax = 240.0;

  force_R = mymalloc("force_R", (Nbin + 1) * (Nbin + 1) * sizeof(double));
  force_z = mymalloc("force_z", (Nbin + 1) * (Nbin + 1) * sizeof(double));

  H = 5.0;                      /* disk scale length */
  M_DISK = 1.0;

  for(i = 0; i <= Nbin; i++)
    for(j = 0; j <= Nbin; j++)
      {
        double RR = exp((log(Rmax / Rmin) * i) / Nbin + log(Rmin));
        double zz = exp((log(Rmax / Rmin) * j) / Nbin + log(Rmin));

        force_z[j * (Nbin + 1) + i] = comp_Dphi_z_disk_exact(RR, zz);
        force_R[j * (Nbin + 1) + i] = comp_Dphi_R_disk_exact(RR, zz);

        double fR = All.G * M_DISK / pow(RR * RR + zz * zz, 1.5) * RR;
        double fz = All.G * M_DISK / pow(RR * RR + zz * zz, 1.5) * zz;

        if(sqrt(RR * RR + zz * zz) > 15.0 * H)
          {
            force_z[j * (Nbin + 1) + i] = fz;
            force_R[j * (Nbin + 1) + i] = fR;
          }
      }

  if(ThisTask == 0)
    printf("done with initializing disk grid\n");

  init_flag = 1;
}





static double comp_Dphi_z_disk_exact(double RR, double zz)
{
  double dphiz;
  double Sigma0;
  double in1, in2, in3, bb;

  R = RR;
  z = zz;
  Sigma0 = (M_DISK) / (2 * M_PI * H * H);

  in1 = qromb(intz_di, 0, 2 / H);

  bb = 2;
  do
    {
      in2 = qromb(intz_di, bb / H, (bb + 2) / H);
      in3 = qromb(intz_di_abs, bb / H, (bb + 2) / H);
      in1 += in2;
      bb += 2;
    }
  while(fabs(in3 / in1) > 1e-3);

  dphiz = 2 * M_PI * All.G * Sigma0 * H * H * (in1);

  return dphiz;
}


static double comp_Dphi_R_disk_exact(double RR, double zz)
{
  double dphiR;
  double Sigma0;
  double in1, in2, in3, bb;


  R = RR;
  z = zz;
  Sigma0 = (M_DISK) / (2 * M_PI * H * H);

  in1 = qromb(intR_di, 0, 2 / H);

  bb = 2;
  do
    {
      in2 = qromb(intR_di, bb / H, (bb + 2) / H);
      in3 = qromb(intR_di_abs, bb / H, (bb + 2) / H);
      in1 += in2;
      bb += 2;
    }
  while(fabs(in3 / in1) > 1e-3);

  dphiR = 2 * M_PI * All.G * Sigma0 * H * H * (in1);

  return dphiR;
}



static double intz_di(double k)
{
  if(z > 0)
    return (bessj0(k * R) * k * exp(-z * k) / pow(1 + k * k * H * H, 1.5));
  else
    return (-bessj0(k * R) * k * exp(z * k) / pow(1 + k * k * H * H, 1.5));
}


static double intz_di_abs(double k)
{
  if(z > 0)
    return fabs(bessj0(k * R) * k * exp(-z * k) / pow(1 + k * k * H * H, 1.5));
  else
    return fabs(-bessj0(k * R) * k * exp(z * k) / pow(1 + k * k * H * H, 1.5));
}


static double intR_di(double k)
{
  if(z >= 0)
    return bessj1(k * R) * k * exp(-z * k) / pow(1 + k * k * H * H, 1.5);
  else
    return bessj1(k * R) * k * exp(z * k) / pow(1 + k * k * H * H, 1.5);
}

static double intR_di_abs(double k)
{
  if(z >= 0)
    return fabs(bessj1(k * R) * k * exp(-z * k) / pow(1 + k * k * H * H, 1.5));
  else
    return fabs(bessj1(k * R) * k * exp(z * k) / pow(1 + k * k * H * H, 1.5));
}

/*
static double mass_cumulative_disk(double R)
{
  return M_DISK * (1 - (1 + R / H) * exp(-R / H));
}
*/



/*** some legacy functions from numerical recipies *****/


static void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
static double trapzd(double (*func) (double), double a, double b, int n);

#define EPS 1.0e-6
#define JMAX 40
#define JMAXP (JMAX+1)
#define K 5

static double qromb(double (*func) (double), double a, double b)
{
  double ss, dss;
  double s[JMAXP], h[JMAXP + 1];
  int j;

  h[1] = 1.0;
  for(j = 1; j <= JMAX; j++)
    {
      s[j] = trapzd(func, a, b, j);
      if(j >= K)
        {
          polint(&h[j - K], &s[j - K], K, 0.0, &ss, &dss);
          if(fabs(dss) <= EPS * fabs(ss))
            return ss;
        }
      h[j + 1] = 0.25 * h[j];
    }
  terminate("Too many steps in routine qromb");
  return 0.0;
}

#undef EPS
#undef JMAX
#undef JMAXP
#undef K

#define FUNC(x) ((*func)(x))

double trapzd(double (*func) (double), double a, double b, int n)
{
  double x, tnm, sum, del;
  static double s;
  int it, j;

  if(n == 1)
    {
      return (s = 0.5 * (b - a) * (FUNC(a) + FUNC(b)));
    }
  else
    {
      for(it = 1, j = 1; j < n - 1; j++)
        it <<= 1;
      tnm = it;
      del = (b - a) / tnm;
      x = a + 0.5 * del;
      for(sum = 0.0, j = 1; j <= it; j++, x += del)
        sum += FUNC(x);
      s = 0.5 * (s + (b - a) * sum / tnm);
      return s;
    }
}

#undef FUNC

#define NRANSI

void polint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
  int i, m, ns = 1;
  double den, dif, dift, ho, hp, w;
  double *c, *d;

  dif = fabs(x - xa[1]);
  c = vector(1, n);
  d = vector(1, n);
  for(i = 1; i <= n; i++)
    {
      if((dift = fabs(x - xa[i])) < dif)
        {
          ns = i;
          dif = dift;
        }
      c[i] = ya[i];
      d[i] = ya[i];
    }
  *y = ya[ns--];
  for(m = 1; m < n; m++)
    {
      for(i = 1; i <= n - m; i++)
        {
          ho = xa[i] - x;
          hp = xa[i + m] - x;
          w = c[i + 1] - d[i];
          if((den = ho - hp) == 0.0)
            terminate("Error in routine polint");
          den = w / den;
          d[i] = hp * den;
          c[i] = ho * den;
        }
      *y += (*dy = (2 * ns < (n - m) ? c[ns + 1] : d[ns--]));
    }
  free_vector(d, 1, n);
  free_vector(c, 1, n);
}

#undef NRANSI

#define NR_END 1
#define FREE_ARG char*

/* allocate a double vector with subscript range v[nl..nh] */
static double *vector(long nl, long nh)
{
  double *v;

  v = (double *) malloc((size_t) ((nh - nl + 1 + NR_END) * sizeof(double)));
  if(!v)
    terminate("allocation failure in vector()");
  return v - nl + NR_END;
}

static void free_vector(double *v, long nl, long nh)
/* free a double vector allocated with vector() */
{
  free((FREE_ARG) (v + nl - NR_END));
}


static double bessj0(double x)
{
  double ax, z;
  double xx, y, ans, ans1, ans2;

  if((ax = fabs(x)) < 8.0)
    {
      y = x * x;
      ans1 = 57568490574.0 + y * (-13362590354.0 + y * (651619640.7 + y * (-11214424.18 + y * (77392.33017 + y * (-184.9052456)))));
      ans2 = 57568490411.0 + y * (1029532985.0 + y * (9494680.718 + y * (59272.64853 + y * (267.8532712 + y * 1.0))));
      ans = ans1 / ans2;
    }
  else
    {
      z = 8.0 / ax;
      y = z * z;
      xx = ax - 0.785398164;
      ans1 = 1.0 + y * (-0.1098628627e-2 + y * (0.2734510407e-4 + y * (-0.2073370639e-5 + y * 0.2093887211e-6)));
      ans2 = -0.1562499995e-1 + y * (0.1430488765e-3 + y * (-0.6911147651e-5 + y * (0.7621095161e-6 - y * 0.934935152e-7)));
      ans = sqrt(0.636619772 / ax) * (cos(xx) * ans1 - z * sin(xx) * ans2);
    }
  return ans;
}


static double bessj1(double x)
{
  double ax, z;
  double xx, y, ans, ans1, ans2;

  if((ax = fabs(x)) < 8.0)
    {
      y = x * x;
      ans1 = x * (72362614232.0 + y * (-7895059235.0 + y * (242396853.1 + y * (-2972611.439 + y * (15704.48260 + y * (-30.16036606))))));
      ans2 = 144725228442.0 + y * (2300535178.0 + y * (18583304.74 + y * (99447.43394 + y * (376.9991397 + y * 1.0))));
      ans = ans1 / ans2;
    }
  else
    {
      z = 8.0 / ax;
      y = z * z;
      xx = ax - 2.356194491;
      ans1 = 1.0 + y * (0.183105e-2 + y * (-0.3516396496e-4 + y * (0.2457520174e-5 + y * (-0.240337019e-6))));
      ans2 = 0.04687499995 + y * (-0.2002690873e-3 + y * (0.8449199096e-5 + y * (-0.88228987e-6 + y * 0.105787412e-6)));
      ans = sqrt(0.636619772 / ax) * (cos(xx) * ans1 - z * sin(xx) * ans2);
      if(x < 0.0)
        ans = -ans;
    }
  return ans;
}

#endif
