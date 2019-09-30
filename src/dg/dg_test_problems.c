/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dg/dg_test_problems.c
 * \date        10/2014
 * \author		Kevin Schaal
 * \brief		initial conditions for various test problems
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

//3d done

#include "../allvars.h"
#include "../proto.h"
#include "dg_test_problems.h"

#ifdef DG_TEST_PROBLEM

#ifndef DG

#ifdef DG_EXTERNAL_ACCELERATION
#ifndef EXTERNALGRAVITY
#error Switch on the flag EXTERNALGRAVITY!
#endif
#endif

#ifndef NO_INI_PROJECTION

static double cell_density(int cell, double x, double y, double z)
{
  return ic_density(x, y, z);
}

static double cell_mom_x_density(int cell, double x, double y, double z)
{
  return ic_density(x, y, z) * ic_velocity_x(x, y, z);
}

static double cell_mom_y_density(int cell, double x, double y, double z)
{
  return ic_density(x, y, z) * ic_velocity_y(x, y, z);
}

static double cell_mom_z_density(int cell, double x, double y, double z)
{
  return ic_density(x, y, z) * ic_velocity_z(x, y, z);
}

static double cell_thermal_energy_density(int cell, double x, double y, double z)
{
  double rho_u = ic_pressure(x, y, z) / (GAMMA - 1);

  return rho_u;
}

#endif

void arepo_set_initial_conditions()
{

#ifndef MESHRELAX_DENSITY_IN_INPUT
  terminate("Error in arepo_set_initial_conditions: compile with MESHRELAX_DENSITY_IN_INPUT!\n");
#endif

  int i;

#ifndef NO_INI_PROJECTION
  double mass;
  double volume;
  double px, py, pz;
  double etherm;
#else
  double x;
  double y;
  double z;
#endif

  for(i = 0; i < NumGas; i++)
    {

#ifndef NO_INI_PROJECTION
      volume = SphP[i].Volume;

      //calculate mass of the cell
      mass = integrate(i, &cell_density);
      px = integrate(i, &cell_mom_x_density);
      py = integrate(i, &cell_mom_y_density);
      pz = integrate(i, &cell_mom_z_density);
      etherm = integrate(i, &cell_thermal_energy_density);
      mpi_printf("cell %d integrated.\n", i);

      P[i].Mass = mass / volume;        // (density in input!)
      P[i].Vel[0] = px / mass;
      P[i].Vel[1] = py / mass;
      P[i].Vel[2] = pz / mass;
      SphP[i].Utherm = etherm / mass;
#else
      x = SphP[i].Center[0];
      y = SphP[i].Center[1];
      z = SphP[i].Center[2];

      P[i].Mass = ic_density(x, y, z);
      P[i].Vel[0] = ic_velocity_x(x, y, z);
      P[i].Vel[1] = ic_velocity_y(x, y, z);
      P[i].Vel[2] = ic_velocity_z(x, y, z);
      SphP[i].Utherm = ic_pressure(x, y, z) / (ic_density(x, y, z) * (GAMMA - 1));
#endif
    }

  // update_primitive_variables();
  // exchange_primitive_variables();
  // calculate_gradients();
  // exchange_primitive_variables_and_gradients();
}

#endif

#ifdef SHU_OSHER_TUBE

//box: [0,20], gamma=1.4, tmax=1.8;
//periodic boundaries (shift everything 14 to the right:
//discontinuity: -4 -> 10
//sin wave:  5x -> 5(x-14)

//static double lambda = 5;
//static double epsilon = 0.2;

double ic_pressure(double x, double y, double z)
{
  if(x < 10)
    {
      return 10.33333;
    }
  else
    {
      return 1;
    }
}

double ic_density(double x, double y, double z)
{
  if(x < 10)
    {
      return 3.857143;
    }
  else
    {
      return 1 + 0.2 * sin(5 * (x - 14));
    }
}



double ic_velocity_x(double x, double y, double z)
{
  if(x < 10)
    {
      return 2.629369;
    }
  else
    {
      return 0;
    }
}

double ic_velocity_y(double x, double y, double z)
{
  return 0;
}

double ic_velocity_z(double x, double y, double z)
{
  return 0;
}


#endif

#ifdef SOUND_WAVE_SKEW

//box: [0,1], gamma=5./3., tmax=sqrt(2);

static double p0 = 3. / 5;
static double rho0 = 1;
static double A = 1e-6;
static double K = 2 * M_PI * sqrt(2);



static double ic_density_old(double x, double y, double z)
{
  return rho0 * (1. + A * sin(K * x));
}

double ic_pressure(double x, double y, double z)
{

  x -= 0.5;
  y -= 0.5;

  double x_old = x;

  x = x / sqrt(2) + y / sqrt(2);
  y = -x_old / sqrt(2) + y / sqrt(2);

  x += 0.5;
  y += 0.5;

  return p0 * pow(ic_density_old(x, y, z) / rho0, GAMMA);
}

double ic_density(double x, double y, double z)
{

  x -= 0.5;
  y -= 0.5;

  double x_old = x;

  x = x / sqrt(2) + y / sqrt(2);
  y = -x_old / sqrt(2) + y / sqrt(2);

  x += 0.5;
  y += 0.5;

  return rho0 * (1. + A * sin(K * x));
}

static double ic_velocity_x_old(double x, double y, double z)
{
  double cs = sqrt(GAMMA * p0 / rho0);
  return A * cs * sin(K * x);
}

double ic_velocity_x(double x, double y, double z)
{
  x -= 0.5;
  y -= 0.5;

  double x_old = x;

  x = x / sqrt(2) + y / sqrt(2);
  y = -x_old / sqrt(2) + y / sqrt(2);

  x += 0.5;
  y += 0.5;

  return ic_velocity_x_old(x, y, z) / sqrt(2);
}

double ic_velocity_y(double x, double y, double z)
{
  x -= 0.5;
  y -= 0.5;

  double x_old = x;

  x = x / sqrt(2) + y / sqrt(2);
  y = -x_old / sqrt(2) + y / sqrt(2);

  x += 0.5;
  y += 0.5;

  return ic_velocity_x_old(x, y, z) / sqrt(2);
}

double ic_velocity_z(double x, double y, double z)
{
  return 0;
}

#endif

#ifdef SOUND_WAVE

//box: [0,1], gamma=5./3., tmax=1

static double p0 = 3. / 5;
static double rho0 = 1;
static double A = 1e-6;

double ic_pressure(double x, double y, double z)
{

  return p0 * pow(ic_density(x, y, z) / rho0, GAMMA);
}

double ic_density(double x, double y, double z)
{

  return rho0 * (1. + A * sin(2 * M_PI * x));
}

double ic_velocity_x(double x, double y, double z)
{
  double cs = sqrt(GAMMA * p0 / rho0);
  return A * cs * sin(2 * M_PI * x);
}

double ic_velocity_y(double x, double y, double z)
{
  return 0;
}

double ic_velocity_z(double x, double y, double z)
{
  return 0;
}

#endif

#ifdef KELVIN_HELMHOLTZ_SINGLE_LAYER

//box: [0,1], gamma=5./3.

double offset = 0.;

double ic_pressure(double x, double y, double z)
{
  return 2.5;
}

double ic_density(double x, double y, double z)
{
  if(y < 0.5 + offset)          //bottom stripe
    {
      return 2;
    }
  else                          //top stripe
    {
      return 1;
    }
}

double ic_velocity_x(double x, double y, double z)
{
  if(y < 0.5 + offset)          //bottom stripe
    {
      return 0.5;
    }
  else                          //top stripe
    {
      return -0.5;
    }
}

double ic_velocity_y(double x, double y, double z)
{
  double omega = 0.1;
  double sigma = 0.05 / sqrt(2);

  return omega * sin(4 * M_PI * x) * exp(-(y - (0.5 + offset)) * (y - (0.5 + offset)) * 0.5 / (sigma * sigma));
}

double ic_velocity_z(double x, double y, double z)
{
  return 0;
}


#endif


#ifdef KELVIN_HELMHOLTZ

//box: [0,1], gamma=5./3.

double offset = 0.;

double ic_pressure(double x, double y, double z)
{
  return 2.5;
}

double ic_density(double x, double y, double z)
{
  if(y < 0.25 + offset)         //bottom stripe
    {
      return 1;
    }
  else if(y < 0.75 + offset)    //middle stripe
    {
      return 2;
    }
  else                          //top stripe
    {
      return 1;
    }
}

double ic_velocity_x(double x, double y, double z)
{
  if(y < 0.25 + offset)         //bottom stripe
    {
      return -0.5;
    }
  else if(y < 0.75 + offset)    //middle stripe
    {
      return 0.5;
    }
  else                          //top stripe
    {
      return -0.5;
    }
}

double ic_velocity_y(double x, double y, double z)
{
  double omega = 0.1;
  double sigma = 0.05 / sqrt(2);

  return omega * sin(4 * M_PI * x) * (exp(-(y - (0.25 + offset)) * (y - (0.25 + offset)) * 0.5 / (sigma * sigma)) + exp(-(y - (0.75 + offset)) * (y - (0.75 + offset)) * 0.5 / (sigma * sigma)));
}

double ic_velocity_z(double x, double y, double z)
{
  return 0;
}

#endif

#ifdef SEDOV

//box: [0,1] (periodic), gamma=5./3., tmax=0.1;

double ic_pressure(double x, double y, double z)
{
  double dl = amr_length[Mesh.DP[0].level];
  double p0 = 1e-6;


  //center
  double cx = 0.5;
  double cy = 0.5;
#ifndef TWODIMS
  double cz = 0.5;
#endif

#ifdef SEDOV_KERNEL
#ifdef TWODIMS
  r = sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy));
#else
  double r = sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy) + (z - cz) * (z - cz));
#endif
  double sigma = 0.4 * dl;

  if(r < dl)
    {
      return 1. / pow(dl, NUMDIMS) * exp(-r * r / (2 * sigma * sigma));
    }
  else
    {
      return p0;
    }

#elif defined(SEDOV_SINGLE_CELL)
  cx = cx + 0.5 * dl;
  cy = cy + 0.5 * dl;
#ifndef TWODIMS
  cz = cz + 0.5 * dl;
#endif

#ifdef TWODIMS
  if(fabs(x - cx) < dl * 0.5 && fabs(y - cy) < dl * 0.5)
    {
      return (GAMMA - 1) / (dl * dl);
    }
  else
    {
      return p0;
    }
#else

  if(fabs(x - cx) < dl * 0.5 && fabs(y - cy) < dl * 0.5 && fabs(z - cz) < dl * 0.5)
    {
      return (GAMMA - 1) / (dl * dl * dl);
    }
  else
    {
      return p0;
    }
#endif


#else
#ifdef TWODIMS
  if(fabs(x - cx) < dl && fabs(y - cy) < dl)
    {
      return 1. / 4. * (GAMMA - 1) / (dl * dl); //factor 1./4.: 4 cells
    }
  else
    {
      return p0;
    }
#else

  if(fabs(x - cx) < dl && fabs(y - cy) < dl && fabs(z - cz) < dl)
    {
      return 1. / 8. * (GAMMA - 1) / (dl * dl * dl);    //factor 1./8.: 8 cells
    }
  else
    {
      return p0;
    }
#endif

#endif

}

double ic_density(double x, double y, double z)
{
  return 1;
}

double ic_velocity_x(double x, double y, double z)
{
  return 0;
}

double ic_velocity_y(double x, double y, double z)
{
  return 0;
}

double ic_velocity_z(double x, double y, double z)
{
  return 0;
}
#endif


#ifdef SHOCK_TUBE

//box: [0,1] (reflective), gamma=1.4, tmax=0.228;

double ic_pressure(double x, double y, double z)
{
  if(x < 0.5)
    {
      return 1;
    }
  else
    {
      return 0.1;
    }
}

double ic_density(double x, double y, double z)
{
  if(x < 0.5)
    {
      return 1;
    }
  else
    {
      return 0.125;
    }
}

double ic_velocity_x(double x, double y, double z)
{
  return 0;
}

double ic_velocity_y(double x, double y, double z)
{
  return 0;
}

double ic_velocity_z(double x, double y, double z)
{
  return 0;
}
#endif

#ifdef SHOCK_TUBE_Z

//box: [0,1] (reflective), gamma=1.4, tmax=0.228;

double ic_pressure(double x, double y, double z)
{
  if(z < 0.5)
    {
      return 1;
    }
  else
    {
      return 0.1;
    }
}

double ic_density(double x, double y, double z)
{
  if(z < 0.5)
    {
      return 1;
    }
  else
    {
      return 0.125;
    }
}

double ic_velocity_x(double x, double y, double z)
{
  return 0;
}

double ic_velocity_y(double x, double y, double z)
{
  return 0;
}

double ic_velocity_z(double x, double y, double z)
{
  return 0;
}
#endif

#ifdef SHOCK_TUBE_SKEW

//box: [0,2] (reflective), gamma=1.4, tmax=0.228;

double ic_pressure(double x, double y, double z)
{
  if(x + y < 2)
    {
      return 1;
    }
  else
    {
      return 0.1;
    }
}

double ic_density(double x, double y, double z)
{
  if(x + y < 2)
    {
      return 1;
    }
  else
    {
      return 0.125;
    }
}

double ic_velocity_x(double x, double y, double z)
{
  return 0;
}

double ic_velocity_y(double x, double y, double z)
{
  return 0;
}

double ic_velocity_z(double x, double y, double z)
{
  return 0;
}
#endif

#ifdef YEE_VORTEX

//box: [0,10], gamma=1.4
//amr suggestion: targetslope=0.01 for 32^2

static double center_x = 5;
static double center_y = 5;
static double beta = 5;

double ic_pressure(double x, double y, double z)
{
  return pow(ic_density(x, y, z), GAMMA);
}

double ic_density(double x, double y, double z)
{
  double r = sqrt(pow(x - center_x, 2) + pow(y - center_y, 2));
  return pow(1 - (GAMMA - 1) * beta * beta / (8 * GAMMA * M_PI * M_PI) * exp(1 - r * r), 1 / (GAMMA - 1));
}

double ic_velocity_x(double x, double y, double z)
{
  double boost = 0;
  double r = sqrt(pow(x - center_x, 2) + pow(y - center_y, 2));
  return beta / (2 * M_PI) * exp((1 - r * r) * 0.5) * (-(y - center_y)) + boost;
}

double ic_velocity_y(double x, double y, double z)
{
  double r = sqrt(pow(x - center_x, 2) + pow(y - center_y, 2));
  return beta / (2 * M_PI) * exp((1 - r * r) * 0.5) * (+(x - center_x));
}

double ic_velocity_z(double x, double y, double z)
{
  return 0;
}

#endif


#ifdef WHATEVER

double ic_pressure(double x, double y, double z)
{
  return 1;
}

double ic_density(double x, double y, double z)
{

  return 4 * x * x + 5;
}

double ic_velocity_x(double x, double y, double z)
{
  return 0;
}

double ic_velocity_y(double x, double y, double z)
{
  return 0;
}

double ic_velocity_z(double x, double y, double z)
{
  return 0;
}

#ifdef DG_EXTERNAL_ACCELERATION

void dg_acceleration(double x, double y, double z, double *acc)
{
  acc[0] = 0;
  acc[1] = 0;
  acc[2] = 0;
}

#endif

#endif

#ifdef RAYLEIGH_TAYLOR

static double g = -0.1;

//box: [0,0.5]x[0,2], gamma=1.4

#ifndef DG_EXTERNAL_ACCELERATION
#error Switch on the external acceleration flag!
#endif

double ic_pressure(double x, double y, double z)
{
  double P0 = 2.5;

  return P0 + g * (y - 1) * ic_density(x, y, z);
}

double ic_density(double x, double y, double z)
{
  if(y > 1)
    {
      return 2;
    }
  else
    {
      return 1;
    }
}

double ic_velocity_x(double x, double y, double z)
{
  return 0;
}

double ic_velocity_y(double x, double y, double z)
{
  double w0 = 0.0025;

  return w0 * (1 - cos(4 * M_PI * x)) * (1 - cos(M_PI * y));
}

double ic_velocity_z(double x, double y, double z)
{
  return 0;
}

double ic_velocity_z(double x, double y, double z)
{
  return 0;
}

void dg_acceleration(double x, double y, double z, double *acc)
{
  acc[0] = 0;
  acc[1] = g;
  acc[2] = 0;
}
#endif

#ifdef KEPLERIAN_DISK

#ifndef DG_EXTERNAL_ACCELERATION
#error Switch on the external acceleration flag!
#endif

//box: [0,10]x[0,10], gamma=5./3.

//center
static double rmin = 0.5;
static double rmax = 2;
static double rim = 0.1;
static double epsilon = 0.25;

double ic_pressure(double x, double y, double z)
{
  return 1e-6;
}

double ic_density(double x, double y, double z)
{
  double cx = 0.5 * boxSize_X;
  double cy = 0.5 * boxSize_Y;

  double rho_zero = 1e-6;
  double rho_disk = 1;

  //coords in center frame
  double xc = x - cx;
  double yc = y - cy;

  double rsquare = xc * xc + yc * yc;
  double r = sqrt(rsquare);

  double rimhalf = rim * 0.5;

  double a = 2 * (rho_disk - rho_zero) / (rim * rim);

  if(r < rmin - rimhalf)
    {
      return rho_zero;
    }
  else if(r < rmin + rimhalf)
    {
      //return (rho_zero - rho_disk) * 0.5 * cos(M_PI / rim * (r - (rmin - rimhalf))) + (rho_zero + rho_disk) * 0.5;
      return (rho_disk - rho_zero) / rim * (r - (rmin - rimhalf)) + rho_zero;

      //parabolas
      if(r < rmin)
        {
          return a * pow((r - (rmin - rimhalf)), 2) + rho_zero;
        }
      else
        {
          return -a * pow((r - (rmin + rimhalf)), 2) + rho_disk;
        }
    }
  else if(r < rmax - rimhalf)
    {
      return rho_disk;
    }
  else if(r < rmax + rimhalf)
    {
      //return (rho_disk - rho_zero) * 0.5 * cos(M_PI / rim * (r - (rmax - rimhalf))) + (rho_disk + rho_zero) * 0.5;
      return (rho_zero - rho_disk) / rim * (r - (rmax - rimhalf)) + rho_disk;

      //parabolas
      if(r < rmax)
        {
          return -a * pow((r - (rmax - rimhalf)), 2) + rho_disk;
        }
      else
        {
          return a * pow((r - (rmax + rimhalf)), 2) + rho_zero;
        }
    }
  else
    {
      return rho_zero;
    }
}

double ic_velocity_x(double x, double y, double z)
{
  double cx = 0.5 * boxSize_X;
  double cy = 0.5 * boxSize_Y;

  //coords in center frame
  double xc = x - cx;
  double yc = y - cy;
  double r = sqrt(xc * xc + yc * yc);
  double v = 1. / sqrt(r);

  return -v * yc / r;
}

double ic_velocity_y(double x, double y, double z)
{
  double cx = 0.5 * boxSize_X;
  double cy = 0.5 * boxSize_Y;

  //coords in center frame
  double xc = x - cx;
  double yc = y - cy;
  double r = sqrt(xc * xc + yc * yc);
  double v = 1. / sqrt(r);

  return v * xc / r;
}

double ic_velocity_z(double x, double y, double z)
{
  return 0;
}

static double dg_acceleration_x(double x, double y, double z)
{
  double cx = 0.5 * boxSize_X;
  double cy = 0.5 * boxSize_Y;

  //coords in center frame
  double xc = x - cx;
  double yc = y - cy;
  double rsquare = xc * xc + yc * yc;
  double r = sqrt(rsquare);

  if(r < rmin - 0.5 * rim)
    {
      return -xc / (r * (rsquare + epsilon * epsilon));
    }
  else
    {
      return -xc / (r * rsquare);
    }
}

static double dg_acceleration_y(double x, double y, double z)
{
  double cx = 0.5 * boxSize_X;
  double cy = 0.5 * boxSize_Y;

  //coords in center frame
  double xc = x - cx;
  double yc = y - cy;
  double rsquare = xc * xc + yc * yc;
  double r = sqrt(rsquare);

  if(r < rmin - 0.5 * rim)
    {
      return -yc / (r * (rsquare + epsilon * epsilon));
    }
  else
    {
      return -yc / (r * rsquare);
    }
}

static double dg_acceleration_z(double x, double y, double z)
{
  return 0;
}

void dg_acceleration(double x, double y, double z, double *acc)
{
  acc[0] = dg_acceleration_x(x, y, z);
  acc[1] = dg_acceleration_y(x, y, z);
  acc[2] = dg_acceleration_z(x, y, z);
}

#endif

#ifdef SQUARE_ADVECTION

//periodic box: [0,1]x[0,1], gamma=1.4, time max: 10, res: 64^2

double ic_pressure(double x, double y, double z)
{
  return 2.5;
}

double ic_density(double x, double y, double z)
{
  if(x < 0.25 || x > 0.75 || y < 0.25 || y > 0.75)
    {
      return 1;
    }
  else
    {
      return 4;
    }
}

double ic_velocity_x(double x, double y, double z)
{
  return 100;
}

double ic_velocity_y(double x, double y, double z)
{
  return 50;
}

double ic_velocity_z(double x, double y, double z)
{
  return 0;
}

#endif

#ifdef DG_TURBULENCE

#ifndef DG_EXTERNAL_ACCELERATION
#error Switch on the external acceleration flag!
#endif


double ic_pressure(double x, double y, double z)
{
  return ic_density(x, y, z) * All.IsoSoundSpeed * All.IsoSoundSpeed;
}

double ic_density(double x, double y, double z)
{
  return 1.;
}

double ic_velocity_x(double x, double y, double z)
{
  return 0.;
}

double ic_velocity_y(double x, double y, double z)
{
  return 0.;
}

double ic_velocity_z(double x, double y, double z)
{
  return 0;
}


#endif

#endif /* DG_TEST_PROBLEM */
