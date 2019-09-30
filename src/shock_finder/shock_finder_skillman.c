/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/shock_finder/shock_finder_skillman.c
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


#ifdef SHOCK_FINDER_SKILLMAN

/*
 * shock finder following Skillman 2008 (TAJ, 689:1063-1077)
 * implemented for test purposes only!
 *
 * Skillman shock finder was made for the AMR code Enzo, this is a generalization for Arepo
 *
 */


static int is_in_shock_zone(int k);
static int is_in_shock_surface(int k);
static void flag_holes(int k);
static void flag_neighbours(int k);
static int is_in_shock_surface2(int k);

void shock_finder_skillman()
{
  int i;

  for(i = 0; i < NumGas; i++)
    {
      //check whether the particle is in the shock zone and set flag correspondingly
      SphP[i].ShockZone = is_in_shock_zone(i);
    }

  for(i = 0; i < NumGas; i++)
    {
      //check whether the particle is part of the shock surface and set flag
      //this call also calculates the machnumber
      SphP[i].ShockSurface = is_in_shock_surface(i);

    }
}




//check whether a cell is in the shock zone
int is_in_shock_zone(int k)
{
  VPRINTF("particle: (%f, %f, %f)\n", P[k].Pos[0], P[k].Pos[1], P[k].Pos[2]);


  double div_v = divergence(k); //divergence of the velocity
  double rho = SphP[k].Density;
  double p = SphP[k].Pressure;
  double T = p / (rho * SfVars.R);      //Temperature

  //Temperature gradient
  double grad_T_x = 1. / SfVars.R * (SphP[k].Grad.dpress[0] * 1. / rho - p / (rho * rho) * SphP[k].Grad.drho[0]);
  double grad_T_y = 1. / SfVars.R * (SphP[k].Grad.dpress[1] * 1. / rho - p / (rho * rho) * SphP[k].Grad.drho[1]);
  double grad_T_z = 1. / SfVars.R * (SphP[k].Grad.dpress[2] * 1. / rho - p / (rho * rho) * SphP[k].Grad.drho[2]);

  VPRINTF("\tTemperature gradient: (%f, %f, %f)\n", grad_T_x, grad_T_y, grad_T_z);

  //Entropy gradient
  double grad_S_x = 1. / pow(rho, GAMMA - 1) * grad_T_x + T * (1 - GAMMA) * pow(rho, -GAMMA) * SphP[k].Grad.drho[0];
  double grad_S_y = 1. / pow(rho, GAMMA - 1) * grad_T_y + T * (1 - GAMMA) * pow(rho, -GAMMA) * SphP[k].Grad.drho[1];
  double grad_S_z = 1. / pow(rho, GAMMA - 1) * grad_T_z + T * (1 - GAMMA) * pow(rho, -GAMMA) * SphP[k].Grad.drho[2];

  VPRINTF("\tEntropy gradient: (%f, %f, %f)\n", grad_S_x, grad_S_y, grad_S_z);


  //div_v < 0 => maybe a shock zone
  if(!(div_v < 0))
    {
      return 0;
    }

  //grad(T)*grad(S) > 0 => maybe a shock zone
  if(!(grad_T_x * grad_S_x + grad_T_y * grad_S_y + grad_T_z * grad_S_z > 0))
    {
      return 0;
    }

  //grad(T)*grad(rho) > 0 => maybe a shock zone
  if(!(grad_T_x * SphP[k].Grad.drho[0] + grad_T_y * SphP[k].Grad.drho[1] + grad_T_z * SphP[k].Grad.drho[2] > 0))
    {
      return 0;
    }

  return 1;
}

//check whether a cell is in the shock surfacen and calculate the machnumber
int is_in_shock_surface(int k)
{
  LPRINTF("particle: (%f, %f, %f)\n", P[k].Pos[0], P[k].Pos[1], P[k].Pos[2]);
  VPRINTF("particle: (%f, %f, %f)\n", P[k].Pos[0], P[k].Pos[1], P[k].Pos[2]);


  if(!SphP[k].ShockZone)
    {
      LPRINTF("not in shock zone\n");
      return 0;
    }

  double div_v = divergence(k); //divergence of the velocity
  double rho = SphP[k].Density;
  double p = SphP[k].Pressure;


  //pre-shock and post-shock temperatures
  double T_preshock;
  double T_postshock;


  //Temperature gradient
  double grad_T[3];
  grad_T[0] = 1. / SfVars.R * (SphP[k].Grad.dpress[0] * 1. / rho - p / (rho * rho) * SphP[k].Grad.drho[0]);
  grad_T[1] = 1. / SfVars.R * (SphP[k].Grad.dpress[1] * 1. / rho - p / (rho * rho) * SphP[k].Grad.drho[1]);
  grad_T[2] = 1. / SfVars.R * (SphP[k].Grad.dpress[2] * 1. / rho - p / (rho * rho) * SphP[k].Grad.drho[2]);

  //search for cells along the temperature gradient
  int current = k;
  int new;
  double *pos = SphP[k].Center;
  double endpoint[3];


  if(grad_T[0] * grad_T[0] + grad_T[1] * grad_T[1] + grad_T[2] * grad_T[2] == 0)
    {
      LPRINTF("denominator is zero\n");
      return 0;
    }

  //plus direction
  while(SphP[current].ShockZone)
    {
      VPRINTF("\tray-pos: (%f, %f, %f)", pos[0], pos[1], pos[2]);
      VPRINTF("\tray-dir: (%.14f, %.14f, %.14f)\n", grad_T[0], grad_T[1], grad_T[2]);

      new = find_neighbouring_cell_pos(current, pos, grad_T, -1, endpoint);

      if(new == current)        //we are at a reflective boundary
        {
          LPRINTF("found same cell\n");
          return 0;
        }
      else
        {
          current = new;
        }

      pos = endpoint;

      VPRINTF("\tneighbour: (%f, %f, %f)\n", P[current].Pos[0], P[current].Pos[1], P[current].Pos[2]);


      if(divergence(current) < div_v && SphP[current].ShockZone)
        {

          LPRINTF("other particle has lower div\n");
          LPRINTF("other: (%f, %f)\n", P[current].Pos[0], P[current].Pos[1]);
          return 0;
        }

    }

  T_postshock = SphP[current].Pressure / (SphP[current].Density * SfVars.R);

  current = k;
  pos = SphP[k].Center;
  grad_T[0] *= -1;
  grad_T[1] *= -1;
  grad_T[2] *= -1;


  //minus direction
  while(SphP[current].ShockZone)
    {
      VPRINTF("\tray-pos: (%f, %f, %f)", pos[0], pos[1], pos[2]);
      VPRINTF("\tray-dir: (%.14f, %.14f, %.14f)\n", grad_T[0], grad_T[1], grad_T[2]);

      new = find_neighbouring_cell_pos(current, pos, grad_T, -1, endpoint);

      if(new == current)        //we are at a reflective boundary
        {
          LPRINTF("found same cell\n");
          return 0;
        }
      else
        {
          current = new;
        }

      pos = endpoint;

      VPRINTF("\tneighbour: (%f, %f, %f)\n", P[current].Pos[0], P[current].Pos[1], P[current].Pos[2]);


      if(divergence(current) < div_v && SphP[current].ShockZone)
        {

          LPRINTF("other particle has lower div\n");
          LPRINTF("other: (%f, %f)\n", P[current].Pos[0], P[current].Pos[1]);
          return 0;
        }
    }

  T_preshock = SphP[current].Pressure / (SphP[current].Density * SfVars.R);

  if(T_postshock < T_preshock)
    {
      LPRINTF("jump in wrong direction\n");
      return 0;
    }

  VPRINTF("\t\tpostshock temperature: %f\n", T_postshock);
  VPRINTF("\t\tpreshock temperature: %f\n", T_preshock);



  // |Delta log T|>=0.11 => jump of a M=1.3 shock or higher
  if(!(fabs(log(T_postshock) - log(T_preshock)) >= 0.11))
    {
      LPRINTF("jump too low\n");
      return 0;
    }

  SphP[k].Machnumber = calculate_machnumber_T_jump(T_preshock, T_postshock);

  VPRINTF("\t\t\tMachnumber: %f\n", SphP[k].Machnumber);


  return 1;
}

#endif
