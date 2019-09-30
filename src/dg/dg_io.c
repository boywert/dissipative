/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dg/dg_io.c
 * \date        08/2015
 * \author      Kevin Schaal
 * \brief       IO routines for the Galerkin (DG) module
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include "../allvars.h"
#include "../proto.h"

#ifdef DG
void io_func_dgw0(int particle, int components, void* buffer, int mode)
{
  int i;

  if(mode == 0)
    {
      MyOutputFloat* out_buffer = buffer;
      for(i = 0; i < components; i++)
        {
          out_buffer[i] = SphP[particle].Weights[i][0];
        }
    }
  else
    {
      MyInputFloat* in_buffer = buffer;
      for(i = 0; i < components; i++)
        {
          SphP[particle].Weights[i][0] = in_buffer[i];
        }
    }
}

void io_func_dgw1(int particle, int components, void* buffer, int mode)
{
  int i;

  if(mode == 0)
    {
      MyOutputFloat* out_buffer = buffer;
      for(i = 0; i < components; i++)
        {
          out_buffer[i] = SphP[particle].Weights[i][1];
        }
    }
  else
    {
      MyInputFloat* in_buffer = buffer;
      for(i = 0; i < components; i++)
        {
          SphP[particle].Weights[i][1] = in_buffer[i];
        }
    }
}

void io_func_dgw2(int particle, int components, void* buffer, int mode)
{
  int i;

  if(mode == 0)
    {
      MyOutputFloat* out_buffer = buffer;
      for(i = 0; i < components; i++)
        {
          out_buffer[i] = SphP[particle].Weights[i][2];
        }
    }
  else
    {
      MyInputFloat* in_buffer = buffer;
      for(i = 0; i < components; i++)
        {
          SphP[particle].Weights[i][2] = in_buffer[i];
        }
    }
}

void io_func_dgw3(int particle, int components, void* buffer, int mode)
{
  int i;

  if(mode == 0)
    {
      MyOutputFloat* out_buffer = buffer;
      for(i = 0; i < components; i++)
        {
          out_buffer[i] = SphP[particle].Weights[i][3];
        }
    }
  else
    {
      MyInputFloat* in_buffer = buffer;
      for(i = 0; i < components; i++)
        {
          SphP[particle].Weights[i][3] = in_buffer[i];
        }
    }
}

void io_func_dgw4(int particle, int components, void* buffer, int mode)
{
  int i;

  if(mode == 0)
    {
      MyOutputFloat* out_buffer = buffer;
      for(i = 0; i < components; i++)
        {
          out_buffer[i] = SphP[particle].Weights[i][4];
        }
    }
  else
    {
      MyInputFloat* in_buffer = buffer;
      for(i = 0; i < components; i++)
        {
          SphP[particle].Weights[i][4] = in_buffer[i];
        }
    }
}

#ifdef DG_SET_IC_FROM_AVERAGES
void load_weights_from_averages()
{
  if(RestartFlag != 1)
  {
    int i;
    double dens;
    mpi_printf("DG_IO: calculating weights from density, velocity and internal energy in ICs \n");
    for(i=0;i<NumGas;i++)
    {
      dens = SphP[i].Density;
      SphP[i].Weights[0][0] = dens;
      SphP[i].Weights[0][1] = dens*P[i].Vel[0];
      SphP[i].Weights[0][2] = dens*P[i].Vel[1];
      SphP[i].Weights[0][3] = dens*P[i].Vel[2];
      SphP[i].Weights[0][4] = dens*SphP[i].Utherm;
      SphP[i].Weights[0][4] += 0.5*dens*P[i].Vel[0]*P[i].Vel[0];
      SphP[i].Weights[0][4] += 0.5*dens*P[i].Vel[1]*P[i].Vel[1];
      SphP[i].Weights[0][4] += 0.5*dens*P[i].Vel[2]*P[i].Vel[2];
    }
  }
}
#endif
#endif

#ifdef OUTPUT_DG_ACCELERATION
void io_func_dg_accel(int particle, int components, void* buffer, int mode)
{
  double acc[3];
  MyOutputFloat* out_buffer = buffer;

  dg_acceleration(SphP[particle].Center[0], SphP[particle].Center[1], SphP[particle].Center[2], acc);

  out_buffer[0] = acc[0];
  out_buffer[1] = acc[1];
  out_buffer[2] = acc[2];
}
#endif

#ifdef OUTPUT_DG_ANGULAR_MOMENTUM
void io_func_dg_angular_momentum(int particle, int components, void* buffer, int mode)
{
  MyOutputFloat* out_buffer = buffer;
  out_buffer[0] = angular_momentum(SphP[particle].Weights, particle);
}
#endif

#ifdef OUTPUT_DG_SPIN
void io_func_dg_spin(int particle, int components, void* buffer, int mode)
{
  MyOutputFloat* out_buffer = buffer;
  out_buffer[0] = spin(SphP[particle].Weights, particle);
}
#endif

#ifdef OUTPUT_DG_TIMESTEP
void io_func_dg_timestep(int particle, int components, void* buffer, int mode)
{
  MyOutputFloat* out_buffer = buffer;
  out_buffer[0] = dg_time_step_cell(particle);
}
#endif

#ifdef OUTPUT_DG_TEMPERATURE
void io_func_dgw_temperature(int particle, int components, void* buffer, int mode)
{
  int i;
  MyOutputFloat* out_buffer = buffer;

  double temperature_weights[NOF_BASE_FUNCTIONS];
  calc_temperature_weights(SphP[particle].Weights, temperature_weights);

  for(i = 0; i < components; i++)
    {
      out_buffer[i] = temperature_weights[i];
    }
}
#endif

#ifdef OUTPUT_DG_U
void io_func_dgw_u(int particle, int components,  void* buffer, int mode)
{
  int i;
  MyOutputFloat* out_buffer = buffer;

  double u_weights[NOF_BASE_FUNCTIONS];
  calc_u_weights(SphP[particle].Weights, u_weights);

  for(i = 0; i < components; i++)
    {
      out_buffer[i] = u_weights[i];
    }
}
#endif

#ifdef OUTPUT_DG_L1_NORM
void io_func_dg_norm(int particle, int components, void* buffer, int mode)
{
  MyOutputFloat* out_buffer = buffer;
  out_buffer[0] = calc_L1_norm(particle);
}
#endif
