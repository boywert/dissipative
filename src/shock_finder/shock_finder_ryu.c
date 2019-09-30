/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/shock_finder/shock_finder_ryu.c
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

#include "../allvars.h"
#include "../proto.h"


#ifdef SHOCK_FINDER_RYU


/*
 * shock finder following Ryu 2003 (TAJ, 593:599-610),
 * implemented for test purposes only!
 *
 * use this finder only for a simulation on a fixed cartesian grid (flag VORONOI_STATIC_GRID).
 * However, for this shock finder, Arepo has to be compiled with VORONOI_DYNAMIC_UPDATE (without the flag VORONOI_STATIC_GRID).
 *
 */

static void shock_finder_dir(double dir_p[3], double dir_m[3]);
static int is_in_shock_zone(int k, double dir_p[3], double dir_m[3]);
static int is_in_shock_surface(int k, double dir_p[3], double dir_m[3]);
static double calculate_machnumber(int k, double dir_p[3], double dir_m[3]);


void shock_finder_ryu()
{
  double x_plus[3] = { 1, 0, 0 };
  double x_minus[3] = { -1, 0, 0 };

  double y_plus[3] = { 0, 1, 0 };
  double y_minus[3] = { 0, -1, 0 };

  double z_plus[3] = { 0, 0, 1 };
  double z_minus[3] = { 0, 0, -1 };

  shock_finder_dir(x_plus, x_minus);
  shock_finder_dir(y_plus, y_minus);
  shock_finder_dir(z_plus, z_minus);
}

void shock_finder_dir(double dir_p[3], double dir_m[3])
{
  int i;
  double calculated_machnumber;

  mpi_printf("Calculating shock zones (direction (%.1f, %.1f, %.1f))\n", dir_p[0], dir_p[1], dir_p[2]);

  for(i = 0; i < NumGas; i++)
    {
      //check whether the particle is in the shock zone and set flag correspondingly
      SphP[i].ShockZone = is_in_shock_zone(i, dir_p, dir_m);

      if(SphP[i].ShockZone)
        {
          SphP[i].ShockZoneGlobal = 1;
        }
    }

  mpi_printf("Calculating shock surface (direction (%.1f, %.1f, %.1f))\n", dir_p[0], dir_p[1], dir_p[2]);

  for(i = 0; i < NumGas; i++)
    {
      //check whether the particle is part of the shock surface and set flag
      SphP[i].ShockSurface = is_in_shock_surface(i, dir_p, dir_m);

      if(SphP[i].ShockSurface)
        {
          SphP[i].ShockSurfaceGlobal = 1;
        }
    }

  mpi_printf("Calculating Mach numbers (direction (%.1f, %.1f, %.1f))\n", dir_p[0], dir_p[1], dir_p[2]);

  for(i = 0; i < NumGas; i++)
    {
      //calculate the Machnumber if the cell is in the shock surface
      calculated_machnumber = calculate_machnumber(i, dir_p, dir_m);

      if(calculated_machnumber > SphP[i].Machnumber)
        {
          SphP[i].Machnumber = calculated_machnumber;
        }
    }

  //reset
  for(i = 0; i < NumGas; i++)
    {
      SphP[i].ShockZone = 0;
      SphP[i].ShockSurface = 0;
    }
}

static int is_in_shock_zone(int k, double dir_p[3], double dir_m[3])
{

  VPRINTF("particle: (%f, %f, %f)\n", P[k].Pos[0], P[k].Pos[1], P[k].Pos[2]);

  int particle_index_p = find_neighbouring_cell(k, dir_p);
  int particle_index_m = find_neighbouring_cell(k, dir_m);

  double rho_p = SphP[particle_index_p].Density;
  double rho_m = SphP[particle_index_m].Density;
  double p_p = SphP[particle_index_p].Pressure;
  double p_m = SphP[particle_index_m].Pressure;
  double vol_p = SphP[particle_index_p].Volume;
  double vol_m = SphP[particle_index_m].Volume;

  double grad_T = (p_p / rho_p - p_m / rho_m) / SfVars.R;       //temperature gradient
  //double grad_S = SfVars.c_v * (log(p_p/pow(rho_p, GAMMA)) * rho_p * vol_p - log(p_m/pow(rho_m, GAMMA)) * rho_m * vol_m); //enropy gradient
  double grad_S = log(p_p / pow(rho_p, GAMMA)) - log(p_m / pow(rho_m, GAMMA));  //gradient of the entropic function
  double grad_log_T = log(p_p / rho_p) - log(p_m / rho_m);

  double div_v = divergence(k); //divergence of the velocity

  VPRINTF("\tplus neighbour: (%f, %f, %f), density: %f, pressure: %f\n", P[particle_index_p].Pos[0], P[particle_index_p].Pos[1], P[particle_index_p].Pos[2], rho_p, p_p);
  VPRINTF("\tminus neighbour: (%f, %f, %f), density: %f, pressure: %f\n", P[particle_index_m].Pos[0], P[particle_index_m].Pos[1], P[particle_index_m].Pos[2], rho_m, p_m);
  VPRINTF("\t\tTemperature gradient: %f\n", grad_T);
  VPRINTF("\t\tEntropy gradient: %f\n", grad_S);
  VPRINTF("\t\tLog T gradient: %f\n", grad_log_T);

  //Delta T x Delta S > 0 => maybe a shock zone
  if(!(grad_T * grad_S > 0))
    {
      return 0;
    }

  //div_v < 0 => maybe a shock zone
  if(!(div_v < 0))
    {
      return 0;
    }

  // |Delta log T|>=0.11 => jump of a M=1.3 shock or higher
  if(!(fabs(grad_log_T) >= 0.11))
    {
      return 0;
    }


  return 1;
}

static int is_in_shock_surface(int k, double dir_p[3], double dir_m[3])
{
  if(!SphP[k].ShockZone)
    {
      return 0;
    }

  int particle_index_p = find_neighbouring_cell(k, dir_p);
  int particle_index_m = find_neighbouring_cell(k, dir_m);


  double div_v = divergence(k); //divergence of the velocity
  double div_v_p = divergence(particle_index_p);        //divergence of the neighbour in plus direction
  double div_v_m = divergence(particle_index_m);        //divergence of the neighbour in plus direction


  if(div_v < div_v_p && div_v < div_v_m)
    {
      return 1;
    }
  else
    {
      return 0;
    }
}

double calculate_machnumber(int k, double dir_p[3], double dir_m[3])
{

  if(!SphP[k].ShockSurface)
    {
      return 0;
    }

  VPRINTF("particle: (%f, %f, %f)\n", P[k].Pos[0], P[k].Pos[1], P[k].Pos[2]);


  int plus_neighbour = find_neighbouring_cell(k, dir_p);
  int minus_neighbour = find_neighbouring_cell(k, dir_m);
  int neighbour;

  while(SphP[plus_neighbour].ShockZone) //search for the first cell in plus direction which is not in the shock zone
    {
      neighbour = plus_neighbour;

      plus_neighbour = find_neighbouring_cell(plus_neighbour, dir_p);

      if(plus_neighbour == neighbour)
        {
          return 0;             //we are at a reflective boundary
        }
    }

  while(SphP[minus_neighbour].ShockZone)        //search for the first cell in minus direction which is not in the shock zone
    {
      neighbour = minus_neighbour;

      minus_neighbour = find_neighbouring_cell(minus_neighbour, dir_m);

      if(minus_neighbour == neighbour)
        {
          return 0;             //we are at a reflective boundary
        }
    }

  VPRINTF("\tplus neighbour: (%f, %f, %f)\n", P[plus_neighbour].Pos[0], P[plus_neighbour].Pos[1], P[plus_neighbour].Pos[2]);
  VPRINTF("\tminus neighbour: (%f, %f, %f)\n", P[minus_neighbour].Pos[0], P[minus_neighbour].Pos[1], P[minus_neighbour].Pos[2]);


  double T_p = SphP[plus_neighbour].Pressure / (SphP[plus_neighbour].Density * SfVars.R);
  double T_l = SphP[minus_neighbour].Pressure / (SphP[minus_neighbour].Density * SfVars.R);

  double T2;                    //post-shock temperature
  double T1;                    //pre-shock temperature

  if(T_p > T_l)
    {
      T2 = T_p;
      T1 = T_l;
    }
  else
    {
      T2 = T_l;
      T1 = T_p;
    }

  VPRINTF("\t\tpost shock T: %f\n", T2);
  VPRINTF("\t\tpre shock T: %f\n", T1);

  double machnumber = calculate_machnumber_T_jump(T1, T2);

  VPRINTF("\t\t\tmachnumber: %f\n", machnumber);

  return machnumber;
}

#endif
