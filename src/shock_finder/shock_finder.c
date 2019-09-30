/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/shock_finder/shock_finder.c
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

#include <stdio.h>
#include <hdf5.h>


#include "../allvars.h"
#include "../proto.h"

#if defined(SHOCK_FINDER_POST_PROCESSING) || defined(SHOCK_FINDER_BEFORE_OUTPUT) || defined(SHOCK_FINDER_ON_THE_FLY)

//throw compiler errors
#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
#ifndef SKIP_BORDER
#error For shock finding with reflective boundaries compile with flag SKIP_BORDER!
#endif
#endif

struct ShockFinderVariables SfVars = {.limit_gradients = 1};


#ifndef SHOCK_FINDER_ON_THE_FLY
ShockData *SData;
ShockData *SDataExch;
#endif

//main function used for the post-processing shock finder
#ifdef SHOCK_FINDER_POST_PROCESSING

static void shock_finder_warnings();
static void shock_finder_info();

void shock_finder()
{
  TIMER_START(CPU_SHOCK_FINDER);

  //prepare state for shock finder
  update_primitive_variables();
  exchange_primitive_variables();
  calculate_gradients(); // unlimited gradients are calculated if the recommended option UNLIMITED_GRADIENTS is used.
  exchange_primitive_variables_and_gradients();

  mpi_printf("\n\n---Entering shock finder---\n\n");

  shock_finder_warnings();

  shock_finder_info();

  mpi_printf("Initializing shock finder\n");
  shock_finder_alloc_memory();
  reset_shock_data();
  print_ini_info();

  calculate_temperatures();
  exchange_shock_data();
  calculate_gradT();
  calculate_divvel();
  exchange_shock_data();

#ifdef SHOCK_FINDER_AREPO
  mpi_printf("Using shock finder for Arepo.\n");
  shock_finder_arepo();
#endif

#ifdef SHOCK_FINDER_RYU
  mpi_printf("Using shock finder following Ryu.\n");
  shock_finder_ryu();
#endif

#ifdef SHOCK_FINDER_SKILLMAN
  mpi_printf("Using shock finder following Skillman.\n");
  shock_finder_skillman();
#endif

  exchange_shock_data();

  TIMER_START(CPU_SF_OUTPUT);

  mpi_printf("Data hdf5 output.\n");

  shock_finder_output();

  TIMER_STOP(CPU_SF_OUTPUT);

  mpi_printf("Finishing shock finder.\n");
  shock_finder_free_memory();

  TIMER_STOP(CPU_SHOCK_FINDER);
}
#else

//main function used for the on-the-fly shock finder and shock finder before output
void shock_finder_on_the_fly()
{
  TIMER_START(CPU_SHOCK_FINDER);

  double tstart = second();

  mpi_printf("\nSHOCKS: Start shock finding algorithm\n");

#ifndef SHOCK_FINDER_ON_THE_FLY
  shock_finder_alloc_memory();
#endif

  reset_shock_data();

  SfVars.limit_gradients = 0;   //need unlimited gradients (density gradient)
  update_primitive_variables(); //needed for calculate gradients
  exchange_primitive_variables();       //needed for calculate gradients
  calculate_gradients();        //needed for calculate_divvel and density gradient

  calculate_temperatures();     //needed for the shock finder

  exchange_primitive_variables_and_gradients(); //needed for calculate_divvel, also exchanges temperatures for on-the-fly shock finding
  exchange_shock_data();        //exchange temperatures

  calculate_gradT();
  calculate_divvel();

//  exchange_shock_data(); //not needed by now
  shock_finder_arepo();
//  exchange_shock_data(); //not needed by now

  //calculate the normal (limited) gradients for the hydro scheme
  SfVars.limit_gradients = 1;
  calculate_gradients();

  calculate_additional_shock_data();

  exchange_primitive_variables_and_gradients();

#ifndef SHOCK_FINDER_ON_THE_FLY
  shock_finder_free_memory();
#endif

  double tend = second();

  mpi_printf("SHOCKS: Shock finding done, took %g secs\n\n", timediff(tstart, tend));

  TIMER_STOP(CPU_SHOCK_FINDER);
}
#endif

//calculate additional variables
void calculate_additional_shock_data()
{

  int idx, i;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

#ifndef COSMIC_RAYS //The energy dissipation with cosmic rays is set in set_machnumber()
      SphP[i].EnergyDissipation = SDATA(i).ShockSurfaceArea * generated_thermal_energy_flux(i);
      //multiply by 1./a * UnitMass_in_g / UnitLength_in_cm * UnitVelocity_in_cm_per_s^3 in order to get erg/s (h factors cancel)
#endif

#ifdef SHOCK_FINDER_BEFORE_OUTPUT_MORE
      SphP[i].ShockSurfaceArea = SDATA(i).ShockSurfaceArea;
      SphP[i].ShockDir[0] = SDATA(i).ShockDir[0];
      SphP[i].ShockDir[1] = SDATA(i).ShockDir[1];
      SphP[i].ShockDir[2] = SDATA(i).ShockDir[2];
      SphP[i].Divvel = SDATA(i).Divvel;

      if(SDATA(i).ShockSurface)
        {
          SphP[i].ZoneFlag = 2;
        }
      else if(SDATA(i).ShockZone)
        {
          SphP[i].ZoneFlag = 1;
        }
      else
        {
          SphP[i].ZoneFlag = 0;
        }
      SphP[i].RhopreShock = SDATA(i).RhopreShock;
      SphP[i].RhopostShock = SDATA(i).RhopostShock;
      SphP[i].VpostShock[0] = SDATA(i).VpostShock[0];
      SphP[i].VpostShock[1] = SDATA(i).VpostShock[1];
      SphP[i].VpostShock[2] = SDATA(i).VpostShock[2];
      SphP[i].VpreShock[0] = SDATA(i).VpreShock[0];
      SphP[i].VpreShock[1] = SDATA(i).VpreShock[1];
      SphP[i].VpreShock[2] = SDATA(i).VpreShock[2];
#ifndef COSMIC_RAYS
      SphP[i].Temperature = SDATA(i).Temperature;
      SphP[i].Tgrad[0] = SDATA(i).Tgrad[0];
      SphP[i].Tgrad[1] = SDATA(i).Tgrad[1];
      SphP[i].Tgrad[2] = SDATA(i).Tgrad[2];
      SphP[i].CpreShock = SDATA(i).CpreShock;
      SphP[i].PpreShock = SDATA(i).PpreShock;
      SphP[i].TpreShock = SDATA(i).TpreShock;
      SphP[i].PpostShock = SDATA(i).PpostShock;
#endif
#endif
    }
}


//exchange shock data
void exchange_shock_data()
{
#ifdef SHOCK_FINDER_ON_THE_FLY
  return;
#else
  if(All.TotNumGas == 0)
    return;

  int listp;
  ShockData *tmpSDataExch;
  int i, j, p, task, off;
  int ngrp, recvTask, place;
  int k;

  tmpSDataExch = (ShockData *) mymalloc("tmpSDataExch", Mesh_nexport * sizeof(ShockData));

  /* prepare data for export */
  for(j = 0; j < NTask; j++)
    Mesh_Send_count[j] = 0;

  for(i = 0; i < NumGasInMesh; i++)
    {
      p = List_InMesh[i];

      if(P[p].Type == 0)
        {
          listp = List_P[p].firstexport;
          while(listp >= 0)
            {
              if((task = ListExports[listp].origin) != ThisTask)
                {
                  place = ListExports[listp].index;
                  off = Mesh_Send_offset[task] + Mesh_Send_count[task]++;

                  tmpSDataExch[off].ShockSurfaceArea = SDATA(place).ShockSurfaceArea;
                  tmpSDataExch[off].ShockZone = SDATA(place).ShockZone;
                  tmpSDataExch[off].ShockSurface = SDATA(place).ShockSurface;
                  tmpSDataExch[off].Divvel = SDATA(place).Divvel;

                  for(k = 0; k < 3; k++)
                    {
                      tmpSDataExch[off].ShockDir[k] = SDATA(place).ShockDir[k];
                    }

                  tmpSDataExch[off].RhopreShock = SDATA(place).RhopreShock;
                  tmpSDataExch[off].RhopostShock = SDATA(place).RhopostShock;
                  tmpSDataExch[off].VpreShock[0] = SDATA(place).VpreShock[0];
                  tmpSDataExch[off].VpreShock[1] = SDATA(place).VpreShock[1];
                  tmpSDataExch[off].VpreShock[2] = SDATA(place).VpreShock[2];

                  tmpSDataExch[off].VpostShock[0] = SDATA(place).VpostShock[0];
                  tmpSDataExch[off].VpostShock[1] = SDATA(place).VpostShock[1];
                  tmpSDataExch[off].VpostShock[2] = SDATA(place).VpostShock[2];
#ifdef COSMIC_RAYS
                  tmpSDataExch[off].CRpseudoTemperature = SDATA(place).CRpseudoTemperature;
                  tmpSDataExch[off].CRpseudoTgrad[0] = SDATA(place).CRpseudoTgrad[0];
                  tmpSDataExch[off].CRpseudoTgrad[1] = SDATA(place).CRpseudoTgrad[1];
                  tmpSDataExch[off].CRpseudoTgrad[2] = SDATA(place).CRpseudoTgrad[2];
#else
                  tmpSDataExch[off].CpreShock = SDATA(place).CpreShock;
                  tmpSDataExch[off].PpreShock = SDATA(place).PpreShock;
                  tmpSDataExch[off].PpostShock = SDATA(place).PpostShock;
                  tmpSDataExch[off].TpreShock = SDATA(place).TpreShock;
                  tmpSDataExch[off].Temperature = SDATA(place).Temperature;
                  tmpSDataExch[off].Tgrad[0] = SDATA(place).Tgrad[0];
                  tmpSDataExch[off].Tgrad[1] = SDATA(place).Tgrad[1];
                  tmpSDataExch[off].Tgrad[2] = SDATA(place).Tgrad[2];
#endif
                }
              listp = ListExports[listp].nextexport;
            }
        }
    }

  /* exchange data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Mesh_Send_count[recvTask] > 0 || Mesh_Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&tmpSDataExch[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask]
                           * sizeof(ShockData), MPI_BYTE, recvTask, TAG_SHOCK_DATA,
                           &SDataExch[Mesh_Recv_offset[recvTask]], Mesh_Recv_count[recvTask] * sizeof(ShockData), MPI_BYTE, recvTask, TAG_SHOCK_DATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  myfree(tmpSDataExch);
#endif
}

//print
void print_ini_info()
{
  mpi_printf("\tNumber of tasks: %d\n", NTask);
  mpi_printf_each("\tTask %d: %d gas particles!\n", ThisTask, NumGas);
  long long int num_gas_global;
  long long int num_gas_local = (long long int) NumGas;
  MPI_Allreduce(&num_gas_local, &num_gas_global, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
  mpi_printf("\tTotal number of gas particles: %lld\n", num_gas_global);

  SfVars.nof_gas_part_total = num_gas_global;
}

//reset variables
void reset_shock_data()
{
  SfVars.nof_gas_part_total = 0;
  SfVars.nof_shockzone_part_total = 0;
  SfVars.nof_shocksurface_part_total = 0;
  SfVars.nof_shockzone_part_local = 0;
  SfVars.nof_shocksurface_part_local = 0;

  int idx, i;
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      SphP[i].Machnumber = 0;
      SphP[i].EnergyDissipation = 0;

      SDATA(i).ShockSurfaceArea = 0;
      SDATA(i).ShockDir[0] = 0;
      SDATA(i).ShockDir[1] = 0;
      SDATA(i).ShockDir[2] = 0;
      SDATA(i).ShockZone = 0;
      SDATA(i).ShockSurface = 0;
      SDATA(i).Divvel = 0;
      SDATA(i).RhopreShock = 0;
      SDATA(i).RhopostShock = 0;
      SDATA(i).VpostShock[0] = 0;
      SDATA(i).VpostShock[1] = 0;
      SDATA(i).VpostShock[2] = 0;
      SDATA(i).VpreShock[0] = 0;
      SDATA(i).VpreShock[1] = 0;
      SDATA(i).VpreShock[2] = 0;

#ifdef COSMIC_RAYS
      SDATA(i).CRpseudoTemperature = 0;
      SDATA(i).CRpseudoTgrad[0] = 0;
      SDATA(i).CRpseudoTgrad[1] = 0;
      SDATA(i).CRpseudoTgrad[2] = 0;
#else
      SDATA(i).CpreShock = 0;
      SDATA(i).PpreShock = 0;
      SDATA(i).PpostShock = 0;
      SDATA(i).TpreShock = 0;
      SDATA(i).Temperature = 0;
      SDATA(i).Tgrad[0] = 0;
      SDATA(i).Tgrad[1] = 0;
      SDATA(i).Tgrad[2] = 0;
#endif
    }
}

//allocate memory
void shock_finder_alloc_memory()
{
#ifndef SHOCK_FINDER_ON_THE_FLY
  SData = (ShockData *) mymalloc_movable(&SData, "SData", NumGas * sizeof(ShockData));
  SDataExch = (ShockData *) mymalloc_movable(&SDataExch, "SDataExch", Mesh_nimport * sizeof(ShockData));
#endif
}

//free memory
void shock_finder_free_memory()
{
#ifndef SHOCK_FINDER_ON_THE_FLY
  myfree(SDataExch);
  myfree(SData);
#endif
}

void calculate_divvel()
{

  double atime, hubble_a;

  point *DP = Mesh.DP;
  face *VF = Mesh.VF;

  int idx, i, j, k, t0, t1, q0, q1, s0, s1;
  double fac0, fac1, n[3], nn;
  MySingle *data;

  double facA, facB, c[3];

#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
  unsigned int *image_flags0, *image_flags1;
#endif

  if(All.ComovingIntegrationOn)
    {
      atime = All.Time;
      hubble_a = hubble_function(atime);
    }
  else
    {
      atime = 1.0;
      hubble_a = 0.0;
    }

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      SDATA(i).Divvel = 0;


      for(k = 0; k < N_Grad; k++)
        {
          if((grad_elements[k].type == GRADIENT_TYPE_VELX) || (grad_elements[k].type == GRADIENT_TYPE_VELY) || (grad_elements[k].type == GRADIENT_TYPE_VELZ))
            {
              data = (MySingle *) (((char *) (&(SphP[i].Grad))) + grad_elements[k].offset_grad);
              for(j = 0; j < 3; j++)
                data[j] = 0;
            }
        }

    }

  for(i = 0; i < Mesh.Nvf; i++)
    {
      point *p1 = &DP[VF[i].p1];
      point *p0 = &DP[VF[i].p2];

      if(p0->index < 0 || p1->index < 0)
        continue;

      n[0] = p1->x - p0->x;
      n[1] = p1->y - p0->y;
      n[2] = p1->z - p0->z;

      nn = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);

      for(j = 0; j < 3; j++)
        n[j] /= nn;

      c[0] = VF[i].cx - 0.5 * (p1->x + p0->x);
      c[1] = VF[i].cy - 0.5 * (p1->y + p0->y);
      c[2] = VF[i].cz - 0.5 * (p1->z + p0->z);

      /* one of the sides */
      q0 = p0->index;
      s0 = p1->index;
      t0 = p1->task;
#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
      image_flags0 = &(p1->image_flags);
#endif

      /* the other side */
      q1 = p1->index;
      s1 = p0->index;
      t1 = p0->task;
#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
      image_flags1 = &(p0->image_flags);
#endif

      for(k = 0; k < N_Grad; k++)
        {
          if((grad_elements[k].type == GRADIENT_TYPE_VELX) || (grad_elements[k].type == GRADIENT_TYPE_VELY) || (grad_elements[k].type == GRADIENT_TYPE_VELZ))
            {
              grad_elements[k].value0 = 0;
              grad_elements[k].value1 = 0;
            }
        }

      /* let's get the physical quantities of interest for the particle on one of the sides */
      if(t0 == ThisTask)
        {
          if(s0 >= NumGas)
            s0 -= NumGas;

          if(s0 >= 0)
            {
              for(k = 0; k < N_Grad; k++)
                {
                  if((grad_elements[k].type == GRADIENT_TYPE_VELX) || (grad_elements[k].type == GRADIENT_TYPE_VELY) || (grad_elements[k].type == GRADIENT_TYPE_VELZ))
                    {
                      grad_elements[k].value0 = *(MyFloat *) (((char *) (&P[s0])) + grad_elements[k].offset);
                    }
                }
            }
        }
      else
        {
          for(k = 0; k < N_Grad; k++)
            {
              if((grad_elements[k].type == GRADIENT_TYPE_VELX) || (grad_elements[k].type == GRADIENT_TYPE_VELY) || (grad_elements[k].type == GRADIENT_TYPE_VELZ))
                {
                  grad_elements[k].value0 = *(MyFloat *) (((char *) (&PrimExch[s0])) + grad_elements[k].offset_exch);
                }
            }
        }

      /* let's get the physical quantities of interest for the particle on the other side */
      if(t1 == ThisTask)
        {
          if(s1 >= NumGas)
            s1 -= NumGas;

          if(s1 >= 0)
            {
              for(k = 0; k < N_Grad; k++)
                {
                  if((grad_elements[k].type == GRADIENT_TYPE_VELX) || (grad_elements[k].type == GRADIENT_TYPE_VELY) || (grad_elements[k].type == GRADIENT_TYPE_VELZ))
                    {
                      grad_elements[k].value1 = *(MyFloat *) (((char *) (&P[s1])) + grad_elements[k].offset);
                    }
                }
            }
        }
      else
        {
          for(k = 0; k < N_Grad; k++)
            {
              if((grad_elements[k].type == GRADIENT_TYPE_VELX) || (grad_elements[k].type == GRADIENT_TYPE_VELY) || (grad_elements[k].type == GRADIENT_TYPE_VELZ))
                {
                  grad_elements[k].value1 = *(MyFloat *) (((char *) (&PrimExch[s1])) + grad_elements[k].offset_exch);
                }
            }
        }

      /* convert velocity components to peculiar velocities */
      if(All.ComovingIntegrationOn)
        {
          GVelx->value0 /= atime;
          GVelx->value1 /= atime;
          GVely->value0 /= atime;
          GVely->value1 /= atime;
          GVelz->value0 /= atime;
          GVelz->value1 /= atime;
        }

#if defined(REFLECTIVE_X)
      if((*image_flags0 & REFL_X_FLAGS) && !(*image_flags0 & OUTFLOW_X))
        GVelx->value0 *= -1;
      if((*image_flags0 & REFL_X_FLAGS) && (*image_flags0) && (REFLECTIVE_X == 3))
        {
          GVely->value0 *= -1;
          GVelz->value0 *= -1;
        }
#endif
#if defined(REFLECTIVE_Y)
      if((*image_flags0 & REFL_Y_FLAGS) && !(*image_flags0 & OUTFLOW_Y))
        GVely->value0 *= -1;
      if((*image_flags0 & REFL_Y_FLAGS) && (*image_flags0) && (REFLECTIVE_Y == 3))
        {
          GVelx->value0 *= -1;
          GVelz->value0 *= -1;
        }
#endif
#if defined(REFLECTIVE_Z)
      if((*image_flags0 & REFL_Z_FLAGS) && !(*image_flags0 & OUTFLOW_Z))
        GVelz->value0 *= -1;
      if((*image_flags0 & REFL_Z_FLAGS) && (*image_flags0) && (REFLECTIVE_Z == 3))
        {
          GVelx->value0 *= -1;
          GVely->value0 *= -1;
        }
#endif

#if defined(REFLECTIVE_X)
      if((*image_flags1 & REFL_X_FLAGS) && !(*image_flags1 & OUTFLOW_X))
        GVelx->value1 *= -1;
      if((*image_flags1 & REFL_X_FLAGS) && (*image_flags0) && (REFLECTIVE_X == 3))
        {
          GVely->value1 *= -1;
          GVelz->value1 *= -1;
        }
#endif
#if defined(REFLECTIVE_Y)
      if((*image_flags1 & REFL_Y_FLAGS) && !(*image_flags1 & OUTFLOW_Y))
        GVely->value1 *= -1;
      if((*image_flags1 & REFL_Y_FLAGS) && (*image_flags0) && (REFLECTIVE_Y == 3))
        {
          GVelx->value1 *= -1;
          GVelz->value1 *= -1;
        }
#endif
#if defined(REFLECTIVE_Z)
      if((*image_flags1 & REFL_Z_FLAGS) && !(*image_flags1 & OUTFLOW_Z))
        GVelz->value1 *= -1;
      if((*image_flags1 & REFL_Z_FLAGS) && (*image_flags0) && (REFLECTIVE_Z == 3))
        {
          GVelx->value1 *= -1;
          GVely->value1 *= -1;
        }
#endif


      /* if the cell q0 is a local particle, construct the gradient estimate */
      if(p0->task == ThisTask && q0 >= 0 && q0 < NumGas)
        {
          if(TimeBinSynchronized[P[q0].TimeBinHydro])
            {
              fac0 = 0.5 * VF[i].area / SphP[q0].Volume;
              facA = VF[i].area / (nn * SphP[q0].Volume);

              for(k = 0; k < N_Grad; k++)
                {
                  if(!(grad_elements[k].type == GRADIENT_TYPE_VELX) || (grad_elements[k].type == GRADIENT_TYPE_VELY) || (grad_elements[k].type == GRADIENT_TYPE_VELZ))
                    {
                      continue;
                    }

                  double value0;

                  if(grad_elements[k].type == GRADIENT_TYPE_VELX)
                    value0 = grad_elements[k].value0 + nn * n[0] * atime * hubble_a;
                  else if(grad_elements[k].type == GRADIENT_TYPE_VELY)
                    value0 = grad_elements[k].value0 + nn * n[1] * atime * hubble_a;
                  else if(grad_elements[k].type == GRADIENT_TYPE_VELZ)
                    value0 = grad_elements[k].value0 + nn * n[2] * atime * hubble_a;
                  else
                    value0 = 0;

                  data = (MySingle *) (((char *) (&(SphP[q0].Grad))) + grad_elements[k].offset_grad);

                  for(j = 0; j < 3; j++)
                    {
                      data[j] += fac0 * n[j] * value0;
                      data[j] += facA * c[j] * (value0 - grad_elements[k].value1);
                    }

                  /* force gradient to be flat if we are at a outflow boundary */
#ifdef REFLECTIVE_X
                  if(((*image_flags0 & REFL_X_FLAGS) && (*image_flags0 & OUTFLOW_X)) || ((*image_flags1 & REFL_X_FLAGS) && (*image_flags1 & OUTFLOW_X)))
                    data[0] = 0;
#endif
#ifdef REFLECTIVE_Y
                  if(((*image_flags0 & REFL_Y_FLAGS) && (*image_flags0 & OUTFLOW_Y)) || ((*image_flags1 & REFL_Y_FLAGS) && (*image_flags1 & OUTFLOW_Y)))
                    data[1] = 0;
#endif
#ifdef REFLECTIVE_Z
                  if(((*image_flags0 & REFL_Z_FLAGS) && (*image_flags0 & OUTFLOW_Z)) || ((*image_flags1 & REFL_Z_FLAGS) && (*image_flags1 & OUTFLOW_Z)))
                    data[2] = 0;
#endif
                }

              /* now add the contribution to the velocity divergence */

              SDATA(q0).Divvel += fac0 * (n[0] * GVelx->value0 + n[1] * GVely->value0 + n[2] * GVelz->value0);
              SDATA(q0).Divvel += facA * (c[0] * (GVelx->value0 - GVelx->value1) + c[1] * (GVely->value0 - GVely->value1) + c[2] * (GVelz->value0 - GVelz->value1));
            }
        }

      /* if the cell q1 is a local particle, construct the gradient estimate */
      if(p1->task == ThisTask && q1 >= 0 && q1 < NumGas)
        {
          if(TimeBinSynchronized[P[q1].TimeBinHydro])
            {
              fac1 = -0.5 * VF[i].area / SphP[q1].Volume;
              facB = VF[i].area / (nn * SphP[q1].Volume);

              for(k = 0; k < N_Grad; k++)
                {
                  if(!(grad_elements[k].type == GRADIENT_TYPE_VELX) || (grad_elements[k].type == GRADIENT_TYPE_VELY) || (grad_elements[k].type == GRADIENT_TYPE_VELZ))
                    {
                      continue;
                    }
                  double value1;


                  if(grad_elements[k].type == GRADIENT_TYPE_VELX)
                    value1 = grad_elements[k].value1 - nn * n[0] * atime * hubble_a;
                  else if(grad_elements[k].type == GRADIENT_TYPE_VELY)
                    value1 = grad_elements[k].value1 - nn * n[1] * atime * hubble_a;
                  else if(grad_elements[k].type == GRADIENT_TYPE_VELZ)
                    value1 = grad_elements[k].value1 - nn * n[2] * atime * hubble_a;
                  else
                    value1 = 0;

                  data = (MySingle *) (((char *) (&(SphP[q1].Grad))) + grad_elements[k].offset_grad);
                  for(j = 0; j < 3; j++)
                    {
                      data[j] += fac1 * n[j] * value1;
                      data[j] += facB * c[j] * (value1 - grad_elements[k].value0);
                    }

                  /* force the gradient to be flat if we are at outflow boundary */
#ifdef REFLECTIVE_X
                  if(((*image_flags0 & REFL_X_FLAGS) && (*image_flags0 & OUTFLOW_X)) || ((*image_flags1 & REFL_X_FLAGS) && (*image_flags1 & OUTFLOW_X)))
                    data[0] = 0;
#endif
#ifdef REFLECTIVE_Y
                  if(((*image_flags0 & REFL_Y_FLAGS) && (*image_flags0 & OUTFLOW_Y)) || ((*image_flags1 & REFL_Y_FLAGS) && (*image_flags1 & OUTFLOW_Y)))
                    data[1] = 0;
#endif
#ifdef REFLECTIVE_Z
                  if(((*image_flags0 & REFL_Z_FLAGS) && (*image_flags0 & OUTFLOW_Z)) || ((*image_flags1 & REFL_Z_FLAGS) && (*image_flags1 & OUTFLOW_Z)))
                    data[2] = 0;
#endif
                }

              SDATA(q1).Divvel += fac1 * (n[0] * GVelx->value1 + n[1] * GVely->value1 + n[2] * GVelz->value1);
              SDATA(q1).Divvel += facB * (c[0] * (GVelx->value1 - GVelx->value0) + c[1] * (GVely->value1 - GVely->value0) + c[2] * (GVelz->value1 - GVelz->value0));
            }
        }
    }
}



/*!
 * calculate the gradient of the temperature for every cell;
 * in the case of cosmic rays the effective temperature is calculated
 */
void calculate_gradT()
{
  point *DP = Mesh.DP;
  face *VF = Mesh.VF;

  int idx, i, j, k, t0, t1, q0, q1, s0, s1;

  double fac0, fac1, n[3], nn;

  double facA, facB, c[3];

  double *grad_min_value;
  double *grad_max_value;

  double gradT_value0;
  double gradT_value1;

#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
  unsigned int *image_flags0, *image_flags1;
#endif

  grad_min_value = mymalloc("grad_min_value", NumGas * sizeof(double));
  grad_max_value = mymalloc("grad_max_value", NumGas * sizeof(double));

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      grad_min_value[i] = +MAX_REAL_NUMBER;
      grad_max_value[i] = -MAX_REAL_NUMBER;

#ifdef COSMIC_RAYS
      SDATA(i).CRpseudoTgrad[0] = 0;
      SDATA(i).CRpseudoTgrad[1] = 0;
      SDATA(i).CRpseudoTgrad[2] = 0;
#else
      SDATA(i).Tgrad[0] = 0;
      SDATA(i).Tgrad[1] = 0;
      SDATA(i).Tgrad[2] = 0;
#endif
    }

  for(i = 0; i < Mesh.Nvf; i++)
    {
      point *p1 = &DP[VF[i].p1];
      point *p0 = &DP[VF[i].p2];

      if(p0->index < 0 || p1->index < 0)
        continue;

      n[0] = p1->x - p0->x;
      n[1] = p1->y - p0->y;
      n[2] = p1->z - p0->z;

      nn = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);

      for(j = 0; j < 3; j++)
        n[j] /= nn;

      c[0] = VF[i].cx - 0.5 * (p1->x + p0->x);
      c[1] = VF[i].cy - 0.5 * (p1->y + p0->y);
      c[2] = VF[i].cz - 0.5 * (p1->z + p0->z);

      /* one of the sides */
      q0 = p0->index;
      s0 = p1->index;
      t0 = p1->task;
#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
      image_flags0 = &(p1->image_flags);
#endif

      /* the other side */
      q1 = p1->index;
      s1 = p0->index;
      t1 = p0->task;
#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
      image_flags1 = &(p0->image_flags);
#endif

      gradT_value0 = 0;
      gradT_value1 = 0;


      /* let's get the physical quantities of interest for the particle on one of the sides */
      if(t0 == ThisTask)
        {
          if(s0 >= NumGas)
            s0 -= NumGas;

          if(s0 >= 0)
            {
#ifdef COSMIC_RAYS
              gradT_value0 = SDATA(s0).CRpseudoTemperature;
#else
              gradT_value0 = SDATA(s0).Temperature;
#endif
            }
        }
      else
        {
#ifdef COSMIC_RAYS
          gradT_value0 = SDATAEXCH(s0).CRpseudoTemperature;
#else
          gradT_value0 = SDATAEXCH(s0).Temperature;
#endif
        }

      /* let's get the physical quantities of interest for the particle on the other side */
      if(t1 == ThisTask)
        {
          if(s1 >= NumGas)
            s1 -= NumGas;

          if(s1 >= 0)
            {
#ifdef COSMIC_RAYS
              gradT_value1 = SDATA(s1).CRpseudoTemperature;
#else
              gradT_value1 = SDATA(s1).Temperature;
#endif
            }
        }
      else
        {
          for(k = 0; k < N_Grad; k++)
            {
#ifdef COSMIC_RAYS
              gradT_value1 = SDATAEXCH(s1).CRpseudoTemperature;
#else
              gradT_value1 = SDATAEXCH(s1).Temperature;
#endif
            }
        }


      /* if the cell q0 is a local particle, construct the gradient estimate, and the minmax values */
      if(p0->task == ThisTask && q0 >= 0 && q0 < NumGas)
        {
          if(TimeBinSynchronized[P[q0].TimeBinHydro])
            {
              fac0 = 0.5 * VF[i].area / SphP[q0].Volume;
              facA = VF[i].area / (nn * SphP[q0].Volume);

              for(j = 0; j < 3; j++)
                {
#ifdef COSMIC_RAYS
                  SDATA(q0).CRpseudoTgrad[j] += fac0 * n[j] * gradT_value0;
                  SDATA(q0).CRpseudoTgrad[j] += facA * c[j] * (gradT_value0 - gradT_value1);
#else
                  SDATA(q0).Tgrad[j] += fac0 * n[j] * gradT_value0;
                  SDATA(q0).Tgrad[j] += facA * c[j] * (gradT_value0 - gradT_value1);
#endif
                }

              /* force gradient to be flat if we are at a outflow boundary */
#ifdef REFLECTIVE_X
              if(((*image_flags0 & REFL_X_FLAGS) && (*image_flags0 & OUTFLOW_X)) || ((*image_flags1 & REFL_X_FLAGS) && (*image_flags1 & OUTFLOW_X)))
                {
#ifdef COSMIC_RAYS
                  SDATA(q0).CRpseudoTgrad[0] = 0;
#else
                  SDATA(q0).Tgrad[0] = 0;
#endif
                }
#endif
#ifdef REFLECTIVE_Y
              if(((*image_flags0 & REFL_Y_FLAGS) && (*image_flags0 & OUTFLOW_Y)) || ((*image_flags1 & REFL_Y_FLAGS) && (*image_flags1 & OUTFLOW_Y)))
                {
#ifdef COSMIC_RAYS
                  SDATA(q0).CRpseudoTgrad[1] = 0;
#else
                  SDATA(q0).Tgrad[1] = 0;
#endif
                }
#endif
#ifdef REFLECTIVE_Z
              if(((*image_flags0 & REFL_Z_FLAGS) && (*image_flags0 & OUTFLOW_Z)) || ((*image_flags1 & REFL_Z_FLAGS) && (*image_flags1 & OUTFLOW_Z)))
                {
#ifdef COSMIC_RAYS
                  SDATA(q0).CRpseudoTgrad[2] = 0;
#else
                  SDATA(q0).Tgrad[2] = 0;
#endif
                }
#endif

              if(VF[i].area > 1.0e-5 * SphP[q0].SurfaceArea)
                {
                  if(grad_max_value[q0] < gradT_value0)
                    grad_max_value[q0] = gradT_value0;

                  if(grad_min_value[q0] > gradT_value0)
                    grad_min_value[q0] = gradT_value0;
                }
            }
        }

      /* if the cell q1 is a local particle, construct the gradient estimate, and the minmax values */
      if(p1->task == ThisTask && q1 >= 0 && q1 < NumGas)
        {
          if(TimeBinSynchronized[P[q1].TimeBinHydro])
            {
              fac1 = -0.5 * VF[i].area / SphP[q1].Volume;
              facB = VF[i].area / (nn * SphP[q1].Volume);

              for(j = 0; j < 3; j++)
                {
#ifdef COSMIC_RAYS
                  SDATA(q1).CRpseudoTgrad[j] += fac1 * n[j] * gradT_value1;
                  SDATA(q1).CRpseudoTgrad[j] += facB * c[j] * (gradT_value1 - gradT_value0);
#else
                  SDATA(q1).Tgrad[j] += fac1 * n[j] * gradT_value1;
                  SDATA(q1).Tgrad[j] += facB * c[j] * (gradT_value1 - gradT_value0);
#endif
                }

              /* force the gradient to be flat if we are at outflow boundary */
#ifdef REFLECTIVE_X
              if(((*image_flags0 & REFL_X_FLAGS) && (*image_flags0 & OUTFLOW_X)) || ((*image_flags1 & REFL_X_FLAGS) && (*image_flags1 & OUTFLOW_X)))
                {
#ifdef COSMIC_RAYS
                  SDATA(q1).CRpseudoTgrad[0] = 0;
#else
                  SDATA(q1).Tgrad[0] = 0;
#endif
                }
#endif
#ifdef REFLECTIVE_Y
              if(((*image_flags0 & REFL_Y_FLAGS) && (*image_flags0 & OUTFLOW_Y)) || ((*image_flags1 & REFL_Y_FLAGS) && (*image_flags1 & OUTFLOW_Y)))
                {
#ifdef COSMIC_RAYS
                  SDATA(q1).CRpseudoTgrad[1] = 0;
#else
                  SDATA(q1).Tgrad[1] = 0;
#endif
                }
#endif
#ifdef REFLECTIVE_Z
              if(((*image_flags0 & REFL_Z_FLAGS) && (*image_flags0 & OUTFLOW_Z)) || ((*image_flags1 & REFL_Z_FLAGS) && (*image_flags1 & OUTFLOW_Z)))
                {
#ifdef COSMIC_RAYS
                  SDATA(q1).CRpseudoTgrad[2] = 0;
#else
                  SDATA(q1).Tgrad[2] = 0;
#endif
                }
#endif

              if(VF[i].area > 1.0e-5 * SphP[q1].SurfaceArea)
                {
                  if(grad_max_value[q1] < gradT_value1)
                    grad_max_value[q1] = gradT_value1;

                  if(grad_min_value[q1] > gradT_value1)
                    grad_min_value[q1] = gradT_value1;
                }
            }
        }
    }

  myfree(grad_max_value);
  myfree(grad_min_value);
}

//calculate the temperature for each cell
void calculate_temperatures()
{
  int idx, i;

  double meanweight;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

#ifdef COOLING
      meanweight = 4. / (3 * HYDROGEN_MASSFRAC + 1 + 4 * HYDROGEN_MASSFRAC * SphP[i].Ne);
#else
      meanweight = 0.5882352941176471;  //fully ionized
#endif

#ifdef COSMIC_RAYS
      SDATA(i).CRpseudoTemperature = (SphP[i].Pressure + SphP[i].CR_Pressure) / SphP[i].Density * All.UnitEnergy_in_cgs / All.UnitMass_in_g * meanweight * PROTONMASS / BOLTZMANN;
#else
      SDATA(i).Temperature = SphP[i].Utherm * All.UnitEnergy_in_cgs / All.UnitMass_in_g * GAMMA_MINUS1 * meanweight * PROTONMASS / BOLTZMANN;
#endif
    }
}

//calculate the shock surface area of a cell
double shock_surface_area(int k)
{
  if(SDATA(k).ShockSurface == 0)
    {
      return 0.0;
    }

#ifdef SURFACE_EXACT
#ifdef TETRA_INDEX_IN_FACE
  return cross_section_plane_cell(k, k, SphP[k].Center, SphP[k].ShockDir);
#else
  terminate("termination in shock_finder.c: shock_surface_area: flag TETRA_INDEX_IN_FACE is missing in Config.sh");
#endif
#endif

#ifdef SURFACE_SPHERE_APPROX
#ifdef ONEDIMS
  return 1;
#endif
#ifdef TWODIMS
  return 2 * sqrt(SphP[k].Volume / M_PI);
#else
  return pow(sqrt(M_PI) * 3. / 4. * SphP[k].Volume, 2. / 3.);
#endif
#endif

#ifdef SURFACE_CUBE_APPROX
#ifdef ONEDIMS
  return 1;
#endif
#ifdef TWODIMS
  terminate("termination in shock_finder.c: shock_surface_area: surface cube approximation in 2D not yet implemented!\n");
#else
  return 1.19 * pow(SphP[k].Volume, 2. / 3.);
#endif
#endif

#ifdef SURFACE_ANGLE_APPROX
#ifdef ONEDIMS
  return 1;
#elif TWODIMS
  terminate("termination in shock_finder.c: shock_surface_area: surface angle approximation in 2D not implemented. Use flag SURFACE_SPHERE_APPROX!\n");
#else
  return 1.074 * pow(SphP[k].MaxFaceAngle, 0.4378) * pow(SphP[k].Volume, 2. / 3.);
#endif
#endif

  assert(0);
  return 0;
}



//density compression factor, needed for delta_M
static inline double dcfR(double M)
{
  return GAMMA_PLUS1 / (GAMMA_MINUS1 + 2 / (M * M));
}

//shock thermalization efficiency, function of the Mach number (Kang et al. 2007)
double delta_M(double M)
{
  assert(M >= 1);

  double result = 2 / (GAMMA * GAMMA_MINUS1 * M * M * dcfR(M)) * ((2 * GAMMA * M * M - GAMMA_MINUS1) / GAMMA_PLUS1 - pow(dcfR(M), GAMMA));

  if(result <= 0)               //can happen for small M due to numerical round-off errors
    {
      return 0;
    }
  else
    {
      return result;
    }
}

//calculate the energy flux (energy per time and area) of a shocked cell
double generated_thermal_energy_flux(int k)
{
#ifdef COSMIC_RAYS
  terminate("For shock finding with cosmic rays call general_thermal_energy_flux_cosmic_rays!\n");
  return 0;
#else

  if(!SDATA(k).ShockSurface)
    {
      return 0;
    }
  else
    {
      double gte_flux = delta_M(SphP[k].Machnumber) * 0.5 * SDATA(k).RhopreShock * pow(SDATA(k).CpreShock * SphP[k].Machnumber, 3);     //rhopreshock: comoving, note: no factors for cpreshock

      return gte_flux;
    }
#endif
}

//calculate the dissipated energy in erg per s
//warning: you might get an overflow with floats for cosmological simulations (ediss > 1e38)
double dissipated_energy_in_erg_per_s(int k)
{
  if(All.Time == 0)
    {
      return 0;
    }

  double surface_physical = SDATA(k).ShockSurfaceArea * (All.Time * All.Time) * All.UnitLength_in_cm * All.UnitLength_in_cm;
  double eflux_physical = generated_thermal_energy_flux(k) / (All.Time * All.Time * All.Time) * All.UnitMass_in_g / pow(All.UnitLength_in_cm,
                                                                                                                        3) * pow(All.UnitVelocity_in_cm_per_s, 3);
  return surface_physical * eflux_physical;     //note h factors cancel
}


//Rankine-Hugoniot Temperature jump condition
double calculate_machnumber_T_jump(double T_preshock, double T_postshock)
{
  myassert(T_preshock > 0);

  double solution1;
  double solution2;

  //solve_quadratic_equation(5, 14 - 16 * T_postshock / T_preshock, -3, &solution1, &solution2);        //gamma = 5/3
  solve_quadratic_equation(2 * GAMMA * (GAMMA - 1), 4 * GAMMA - (GAMMA - 1) * (GAMMA - 1) - (GAMMA + 1) * (GAMMA + 1) * T_postshock / T_preshock, -2 * (GAMMA - 1), &solution1, &solution2);    //general gamma

  return sqrt(solution1);
}

//Rankine-Hugoniot Density jump condition
double calculate_machnumber_rho_jump(double rho_preshock, double rho_postshock)
{
  myassert(rho_preshock > 0);

  double jump = rho_postshock/rho_preshock;

  myassert(jump<(GAMMA+1.)/(GAMMA-1.));

  return sqrt(2*jump/((GAMMA+1)-jump*(GAMMA-1)));
}

//Rankine-Hugoniot Pressure jump condition
double calculate_machnumber_p_jump(double p_preshock, double p_postshock)
{
  myassert(p_preshock > 0);

  //return sqrt((p_postshock / p_preshock * 4 + 1) / 5);        //gamma = 5/3
  return sqrt((GAMMA + 1) * p_postshock / (2 * GAMMA * p_preshock) + (GAMMA - 1) / (2 * GAMMA));        //general gamma
}

//velocity-jump condition, uses the difference instead of the ratio, note: dv has to be the PHYSICAL peculiar velocity (P[i].Vel / a)
double calculate_machnumber_vdiff(double dv, double c_preshock)
{
  myassert(dv <= 0);
  myassert(c_preshock > 0);     //c_preshock == 0: infinite Mach number

  double solution_plus;
  double solution_minus;

  solve_quadratic_equation_pq(dv / c_preshock * (GAMMA + 1) * 0.5, -1, &solution_plus, &solution_minus);

  return solution_plus;
}

//jump condition for the entropy
double calculate_machnumber_S_jump(double S_preshock, double S_postshock, double guess)
{
  //Solve equation with Newton-Raphson iterations
  //

  myassert(S_preshock > 0);
  double jump = S_postshock / S_preshock;

  static unsigned int max_iterations = 20;
  double tolerance = 1e-8;

  double m = guess;             //current mach number
  double m_next, numerator, denominator, s_1, s_2, ds_1, ds_2;

  unsigned int i;

  for(i = 0; i < max_iterations; i++)
    {

      //jump = s_1 * s_2
      s_1 = (2 * GAMMA * m * m - GAMMA + 1) / (GAMMA + 1);
      s_2 = pow((GAMMA - 1 + 2 / (m * m)) / (GAMMA + 1), GAMMA);

      //the derivatives of s_1 and s_2
      ds_1 = (4 * GAMMA * m) / (GAMMA + 1);
      ds_2 = GAMMA * pow((GAMMA - 1 + 2 / (m * m)) / (GAMMA + 1), GAMMA - 1) * (-4 / (m * m * m)) / (GAMMA + 1);

      numerator = s_1 * s_2 - jump;
      denominator = ds_1 * s_2 + s_1 * ds_2;

      myassert(denominator > 0);

      m_next = m - numerator / denominator;

      if(2 * fabs((m_next - m) / (m_next + m)) < tolerance)     //converged
        {
          m = m_next;
          break;
        }
      else                      //not yet converged
        {
          m = m_next;
        }
    }

  if(i == max_iterations)
    {
      terminate("Error in calculate_machnumber_S_jump: divergence!\n");
    }

  return m;
}

#ifdef COSMIC_RAYS
/*!
 * Calculate the Mach number for shock finding with cosmic rays using the mass and momentum equation (default).
 * This implementation is only stable for significant jumps.
 *
 * \param pth1 the pre shock thermal pressure
 * \param pth2 the post shock thermal pressure
 * \param pcr1 the pre shock cosmic ray pressure
 * \param pcr2 the post shock cosmic ray pressure
 * \param rho1 the pre shock density
 * \param rho2 the post shock density
 */
double calculate_machnumber_cosmic_rays(double pth1, double pth2, double pcr1, double pcr2, double rho1, double rho2)
{
  double xs = rho2 / rho1;
  double gammaEff1 = calculate_gamma_effective(pcr1, pth1);
  double Msquare = ((pth2 + pcr2) / (pth1 + pcr1) - 1) * xs / (gammaEff1 * (xs - 1));

  return sqrt(Msquare);
}

/*!
 * Calculate the Mach number for shock finding with cosmic rays using the mass and energy equation.
 * This implementation is only stable for significant jumps.
 *
 * \param pth1 the pre shock thermal pressure
 * \param pth2 the post shock thermal pressure
 * \param pcr1 the pre shock cosmic ray pressure
 * \param pcr2 the post shock cosmic ray pressure
 * \param rho1 the pre shock density
 * \param rho2 the post shock density
 */
double calculate_machnumber_cosmic_rays_from_energy_equation(double pth1, double pth2, double pcr1, double pcr2, double rho1, double rho2)
{
  assert(rho1 > 0);
  assert(rho2 > 0);

  double xs = rho2 / rho1;
  double gammaEff1 = calculate_gamma_effective(pcr1, pth1);
  double gammaEff2 = calculate_gamma_effective(pcr2, pth2);
  double ptot1 = pcr1 + pth1; 
  double ptot2 = pcr2 + pth2; 
  double c1_squared = gammaEff1 * ptot1 / rho1;
  double c2_squared = gammaEff2 * ptot2 / rho2;
  double eps1 = pth1 / GAMMA_MINUS1 + pcr1 / (All.GammaCR - 1.);
  double eps2 = pth2 / GAMMA_MINUS1 + pcr2 / (All.GammaCR - 1.);

  assert(ptot1 > 0);
  assert(ptot2 > 0);

  double Msquare = 2. * xs * xs / (xs * xs - 1.) * ( (eps2 / ptot2 + 1.) * c2_squared / c1_squared / gammaEff2 - (eps1 / ptot1 + 1.) / gammaEff1  );

  return sqrt(Msquare);
}

/*!
 * Calculate the Mach number for shock finding with cosmic rays.
 * This implementation is only stable for different gammas.
 *
 * \param p1th the pre shock thermal pressure
 * \param p2th the post shock thermal pressure
 * \param p1cr the pre shock cosmic ray pressure
 * \param p2cr the post shock cosmic ray pressure
 * \param gamma1 the pre shock gamma defined by the eos (gamma=p/e+1)
 * \param gamma2 the post shock gamma defined by the eos (gamma=p/e+1)
 */
double calculate_machnumber_cosmic_rays_zone(double pth1, double pth2, double pcr1, double pcr2, double gamma1, double gamma2)
{
  assert(pcr1 + pth1 > 0);

  double y = (pth2 + pcr2) / (pth1 + pcr1);

  double gamma_eff1 = (All.GammaCR * pcr1 + GAMMA * pth1) / (pth1 + pcr1);

  double Y = ((gamma2 + 1) * y + gamma2 - 1) * (gamma1 - 1);

  double Msquare = ((y - 1) * Y) / (gamma_eff1 * (Y - ((gamma1 + 1) + (gamma1 - 1) * y) * (gamma2 - 1)));

  return sqrt(Msquare);
}

/*!
 * Calculate the generated thermal energy flux for shock finding with cosmic rays
 *
 * \param eth1 the pre shock thermal energy density (energy / volume)
 * \param eth2 the post shock thermal energy density (energy / volume)
 * \param ecr1 the pre shock cosmic ray energy density (energy / volume)
 * \param ecr2 the post shock cosmic ray energy density (energy / volume)
 * \param rho1 the pre shock density
 * \param rho2 the post shock density
 * \param M the Mach number
 * \param c1 the pre shock sound speed
 */
double generated_thermal_energy_flux_cosmic_rays(double eth1, double eth2, double ecr1, double ecr2, double rho1, double rho2, double M, double c1)
{
  double xs = rho2 / rho1;

  return (eth2 + ecr2 - eth1 * pow(xs, GAMMA) - ecr1 * pow(xs, All.GammaCR)) * M * c1 / xs;
}

/*!
 * Calculate the effective gamma
 *
 * \param pcr the cosmic ray pressure
 * \param pth the thermal pressure
 */
double calculate_gamma_effective(double pcr, double pth)
{
  return (All.GammaCR * pcr + GAMMA * pth) / (pcr + pth);
}

/*!
 * Calculate the sound speed in the case of cosmic rays
 *
 * \param pcr: cosmic ray pressure
 * \param pth: thermal pressure
 * \param rho: density
 * \return the sound speed
 *
 */
double sound_speed_cosmic_rays(double pcr, double pth, double rho)
{
  return sqrt(calculate_gamma_effective(pcr, pth) * (pcr + pth) / rho);
}

#ifdef COSMIC_RAYS_SHOCK_ACCELERATION

/*!
 * Calculate the angle of the magnetic shock obliquity from the shock direction and B field.
 */

double magnetic_shock_obliquity(double *sdir, double *bfield)
{
  double b_norm = sqrt(bfield[0] * bfield[0] + bfield[1] * bfield[1] + bfield[2] * bfield[2]);

  //maximum efficiency for vanishing B field.
  if(b_norm == 0)
    {
      return 0;
    }

  double theta = (sdir[0] * bfield[0] + sdir[1] * bfield[1] + sdir[2] * bfield[2]) / b_norm;

  theta = fmin(theta, fabs(M_PI - theta));

  return acos(theta);
}


/*!
 * Returns the weight for CR injection
 */
double injection_weight(int cell)
{
#ifdef CR_INJECTION_THERM
  return (SphP[cell].Pressure / GAMMA_MINUS1) * SphP[cell].Volume;
#else
  return ((SphP[cell].Pressure / GAMMA_MINUS1) + (SphP[cell].CR_Pressure / (All.GammaCR - 1.))) * SphP[cell].Volume;
#endif
}

/*!
 * Returns the weight for CR injection from PrimExch
 */
double injection_weight_exch(int cell)
{
#ifdef CR_INJECTION_THERM
  return (PrimExch[cell].Pressure / GAMMA_MINUS1) * PrimExch[cell].Volume;
#else
  return ((PrimExch[cell].Pressure / GAMMA_MINUS1) + (PrimExch[cell].CR_Pressure / (All.GammaCR - 1.))) * PrimExch[cell].Volume;
#endif
}
#endif
#endif

/*! Resets the position of a ray in the periodic box
 *  if it crosses the box.
 *
 *  \param pos The position of the ray
 *  \param particle_index The index of the found particle
 */
static void restore_pos_periodic(double pos[3], int particle_index)
{
  double pos_particle[3];
  pos_particle[0] = P[particle_index].Pos[0];
  pos_particle[1] = P[particle_index].Pos[1];
  pos_particle[2] = P[particle_index].Pos[2];

  double dx = pos[0] - pos_particle[0];
  double dy = pos[1] - pos_particle[1];
  double dz = pos[2] - pos_particle[2];

  if(dx > boxHalf_X)
    {
      pos[0] -= boxSize_X;
    }
  else if(dx < -boxHalf_X)
    {
      pos[0] += boxSize_X;
    }

  if(dy > boxHalf_Y)
    {
      pos[1] -= boxSize_Y;
    }
  else if(dy < -boxHalf_Y)
    {
      pos[1] += boxSize_Y;
    }

  if(dz > boxHalf_Z)
    {
      pos[2] -= boxSize_Z;
    }
  else if(dz < -boxHalf_Z)
    {
      pos[2] += boxSize_Z;
    }


}

//! Finds the next Voronoi cell starting at the center of a cell, single task version
/*!
 *  \param k The index of current cell
 *  \param dir The direction in which to search the next cell
 *  \result The index of the found cell
 */
int find_neighbouring_cell(int k, double dir[3])
{
  double length;

  //delaunay connection
  int dc = find_next_voronoi_cell(&Mesh, k, P[k].Pos, dir, -1, &length);

  //delaunay point index
  int dp_index = DC[dc].dp_index;

  //particle indices
  int particle_index = Mesh.DP[dp_index].index;

  //handle borders
  if(particle_index >= NumGas && Mesh.DP[dp_index].task == ThisTask)
    {
      particle_index -= NumGas;
    }

  return particle_index;
}


//! Finds the next Voronoi cell starting at the center of a cell, parallel version
/*!
 *  \param k The index of current cell
 *  \param dir The direction in which to search the next cell
 *  \param task Output value, the task of the found cell
 *  \param original_index Output value, the hydro cell index on the original task
 *  \param edge Output value, the found connection (put this in DC)
 *  \param boundary_reached, indicates whether the boundary is reached and stores the image flag in this case
 *  \result The index of the found cell
 */
int find_neighbouring_cell_para(int k, double dir[3], int *task, int *original_index, int *edge, int *boundary_reached)
{
  double length;
  *boundary_reached = 0;

  //connection
  *edge = find_next_voronoi_cell(&Mesh, k, P[k].Pos, dir, -1, &length);

  //delaunay point index
  int dp_index = DC[*edge].dp_index;

  //particle indices
  int particle_index = Mesh.DP[dp_index].index;

  //check whether the found Delaunay point is imported
  if(Mesh.DP[dp_index].task != ThisTask)
    {
      *task = Mesh.DP[dp_index].task;
      *original_index = Mesh.DP[dp_index].originalindex;
      return particle_index;
    }

  //handle borders
  if(particle_index >= NumGas && Mesh.DP[dp_index].task == ThisTask)
    {
      *boundary_reached = Mesh.DP[dp_index].image_flags;
      particle_index -= NumGas;
    }

  *task = ThisTask;
  *original_index = -1;
  return particle_index;
}

//! Finds the next Voronoi cell, single task version
/*!
 *  \param k The index of the current cell
 *  \param start_point The point where the ray starts
 *  \param dir The direction of the ray
 *  \param previous The cell from which the ray came (needed for stability reasons in find_next_voronoi_cell)
 *  \param end_point Output value, the intersection of the ray with the next cell
 *  \result The index of the neighbouring cell
 */
int find_neighbouring_cell_pos(int k, double start_point[3], double dir[3], int previous, double *end_point)
{
  double dir_l = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
  myassert(dir_l > 0);

  double length;

#ifndef DOUBLEPRECISION
  MyDouble start_point_md[3];

  start_point_md[0] = start_point[0];
  start_point_md[1] = start_point[1];
  start_point_md[2] = start_point[2];

  int dc = find_next_voronoi_cell(&Mesh, k, start_point_md, dir, previous, &length);
#else
  //delaunay connection
  int dc = find_next_voronoi_cell(&Mesh, k, start_point, dir, previous, &length);
#endif

  //delaunay point index
  int dp_index = DC[dc].dp_index;

  //particle indices
  int particle_index = Mesh.DP[dp_index].index;

  //handle borders
  if(particle_index >= NumGas && Mesh.DP[dp_index].task == ThisTask)
    {
      particle_index -= NumGas;

      restore_pos_periodic(end_point, particle_index);
    }

  //calculate end point
  end_point[0] = start_point[0] + dir[0] * length / dir_l;
  end_point[1] = start_point[1] + dir[1] * length / dir_l;
  end_point[2] = start_point[2] + dir[2] * length / dir_l;

  return particle_index;
}



//! Finds the next Voronoi cell, parallel version
/*!
 *  \param k The index of the current cell
 *  \param start_point The point where the ray starts
 *  \param dir The direction of the ray
 *  \param previous The cell from which the ray came (needed for stability reasons in find_next_voronoi_cell)
 *  \param end_point Output value, the intersection of the ray with the next cell
 *  \param task Output value, the task of the found cell
 *  \param original_index Output value, the hydro cell index on the original task
 *  \param edge Output value, the found connection (put this in DC)
 *  \param boundary_reached, indicates whether the boundary is reached and stores the image flag in this case
 *  \result The index of the neighbouring cell
 */
int find_neighbouring_cell_pos_para(int k, double start_point[3], double dir[3], int previous, int task_of_previous, double *end_point, int *task, int *original_index, int *edge,
                                    int *boundary_reached)
{
  double dir_l = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
  myassert(dir_l > 0);

  double length;
  *boundary_reached = 0;

#ifndef DOUBLEPRECISION
  MyDouble start_point_md[3];

  start_point_md[0] = start_point[0];
  start_point_md[1] = start_point[1];
  start_point_md[2] = start_point[2];

  *edge = find_next_voronoi_cell2(&Mesh, k, start_point_md, dir, previous, task_of_previous, &length);
#else
  //connection
  *edge = find_next_voronoi_cell2(&Mesh, k, start_point, dir, previous, task_of_previous, &length);
#endif

  //delaunay point index
  int dp_index = DC[*edge].dp_index;

  //particle indices
  int particle_index = Mesh.DP[dp_index].index;

  //calculate end point
  end_point[0] = start_point[0] + dir[0] * length / dir_l;
  end_point[1] = start_point[1] + dir[1] * length / dir_l;
  end_point[2] = start_point[2] + dir[2] * length / dir_l;

  //check whether the found Delaunay point is imported
  if(Mesh.DP[dp_index].task != ThisTask)
    {
      *task = Mesh.DP[dp_index].task;
      *original_index = Mesh.DP[dp_index].originalindex;
      return particle_index;
    }

  //handle borders
  if(particle_index >= NumGas && Mesh.DP[dp_index].task == ThisTask)
    {
      *boundary_reached = Mesh.DP[dp_index].image_flags;;
      particle_index -= NumGas;

      restore_pos_periodic(end_point, particle_index);
    }

  *task = ThisTask;
  *original_index = -1;
  return particle_index;
}


//a stable way of solving a quadratic equation (ax^2+bx+c=0); solution plus/minus corresponds to the abc-formula solution
void solve_quadratic_equation(double a, double b, double c, double *solution_plus, double *solution_minus)
{
  myassert(a != 0);

  double p = -b / a;
  double q = c / a;

  if(p > 0)
    {
      if(a > 0)
        {
          *solution_plus = p * 0.5 + sqrt(p * p * 0.25 - q);
          *solution_minus = q / (*solution_plus);
        }
      else
        {
          *solution_minus = p * 0.5 + sqrt(p * p * 0.25 - q);
          *solution_plus = q / (*solution_minus);
        }
    }

  else                          //p<=0
    {
      if(a > 0)
        {
          *solution_minus = p * 0.5 - sqrt(p * p * 0.25 - q);
          *solution_plus = q / (*solution_minus);
        }
      else
        {
          *solution_plus = p * 0.5 - sqrt(p * p * 0.25 - q);
          *solution_minus = q / (*solution_plus);
        }
    }
}

//a stable way of solving a quadratic equation with the p-q formula ( x^2 + px + q = 0)
//note: the plus/minus solution of the pq-formula does not always correspond to the plus/minus solution of the abc-formula, the solutions may be swaped
void solve_quadratic_equation_pq(double p, double q, double *solution_plus, double *solution_minus)
{
  p *= -1;

  if(p > 0)
    {
      *solution_plus = p * 0.5 + sqrt(p * p * 0.25 - q);
      *solution_minus = q / (*solution_plus);
    }

  else                          //p<=0
    {
      *solution_minus = p * 0.5 - sqrt(p * p * 0.25 - q);
      *solution_plus = q / (*solution_minus);

    }
}

double distance_to_border_x(int cell)
{
  double d1 = boxSize_X - P[cell].Pos[0];
  assert(d1 > 0);

  double d2 = P[cell].Pos[0];

  return fmin(d1, d2);
}

double distance_to_border_y(int cell)
{
  double d1 = boxSize_Y - P[cell].Pos[1];
  assert(d1 > 0);

  double d2 = P[cell].Pos[1];

  return fmin(d1, d2);
}

double distance_to_border_z(int cell)
{
  double d1 = boxSize_Z - P[cell].Pos[2];
  assert(d1 > 0);

  double d2 = P[cell].Pos[2];

  return fmin(d1, d2);
}

#ifdef SHOCK_FINDER_POST_PROCESSING

static void shock_finder_info()
{
  mpi_printf("Total particle number: %d\n", All.TotNumPart);

#ifdef SHOCK_DIR_GRAD_T
  mpi_printf("Using the temperature gradient for calculating the shock direction.\n");
#endif

#ifdef LOCAL_DIVERGENCE
  mpi_printf("Using local divergence\n");
#endif

#ifdef UNLIMITED_GRADIENTS
  mpi_printf("Using unlimited gradients\n");
#endif

  mpi_printf("\n");
}

static void shock_finder_warnings()
{
  //check parameter file values
  if(All.RayStepsMax <= 0)
    {
      terminate("RayStepsMax has to be greater equal zero! Set it to 100 for no limitation.\n");
    }

  if(All.NumFilesPerOutput <= 0 || (NTask % All.NumFilesPerOutput != 0))
    {
      terminate("The number of tasks has to be a multiple of NumFilesPerOutput, change NumFilesPerOutput in param.txt!\n");
    }

  int shock_finder_counter = 0;
  int shock_dir_counter = 0;
  int shock_jump_counter = 0;
  int zone_jump_counter = 0;
  int surface_approx_counter = 0;

#ifdef SHOCK_FINDER_RYU
  shock_finder_counter++;
#endif

#ifdef SHOCK_FINDER_SKILLMAN
  shock_finder_counter++;
#endif

#ifdef SHOCK_FINDER_AREPO
  shock_finder_counter++;
#endif

#if (defined(SHOCK_FINDER_RYU) || defined(SHOCK_FINDER_SKILLMAN))
  if(NTask > 1)
    {
      terminate("Chosen shockfinder only works on one processor.\n");
    }

#endif

  if(shock_finder_counter != 1)
    {
      terminate("Activate at least/only one shock finder in Config.sh");
    }

#if (defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z))
  mpi_printf("WARNING: Shockfinder does not work directly at the reflective boundaries!\n");

#if defined(SHOCK_FINDER_AREPO) && !defined(SKIP_BORDER)
  terminate("Arepo shock finder works only for periodic boundaries so far!\n");
#endif
#endif

#ifdef SHOCK_DIR_GRAD_T
  shock_dir_counter++;
#endif

#ifdef SHOCK_DIR_GRAD_P
  shock_dir_counter++;
#endif

#ifdef SHOCK_DIR_GRAD_S
  shock_dir_counter++;
#endif

  if(shock_dir_counter != 1)
    {
      terminate("Activate at least/only one flag for calculating the shock direction in Config.sh");
    }

#ifdef SHOCK_JUMP_V
  shock_jump_counter++;
#endif

#ifdef SHOCK_JUMP_P
  shock_jump_counter++;
#endif

#ifdef SHOCK_JUMP_T
  shock_jump_counter++;
#endif

#ifdef SHOCK_JUMP_S
  shock_jump_counter++;
#endif


  if(shock_jump_counter != 1)
    {
      terminate("Activate at least/only one flag for calculating the shock jump in Config.sh");
    }

#ifdef ZONE_JUMP_V
  zone_jump_counter++;
#endif

#ifdef ZONE_JUMP_P
  zone_jump_counter++;
#endif

#ifdef ZONE_JUMP_S
  zone_jump_counter++;
#endif

#ifdef ZONE_JUMP_T
  zone_jump_counter++;
#endif

  if(zone_jump_counter < 1)
    {
      terminate("Activate at least/only one flag for calculating the zone jump in Config.sh");
    }

#ifdef SURFACE_SPHERE_APPROX
  surface_approx_counter++;
#endif

#ifdef SURFACE_CUBE_APPROX
  surface_approx_counter++;
#endif

#ifdef SURFACE_EXACT
  surface_approx_counter++;
#endif

#ifdef SURFACE_ANGLE_APPROX
  surface_approx_counter++;
#endif

  if(surface_approx_counter != 1)
    {
      terminate("Activate at least/only one flag for calculating surface area of shocked cells in Config.sh");
    }
}

//number of particles already written
static unsigned int Offset = 0;

//buffer size in number of particles before sending to Task 0 and writing to disk
static unsigned int Buffer_size = 100000;

//writes out a vector quantity, e.g. coordinates (x1,y1,z1,x2,y2,z2,...)
static void write_float_buffer_vector(hid_t file_id, unsigned int nof_particles, float *buffer, const char *dataset_name)
{
  //identifier
  herr_t status;

  //open dataset and get dataspace
  hid_t dataset = my_H5Dopen(file_id, dataset_name);
  hid_t filespace = my_H5Dget_space(dataset, dataset_name);

  //velocity file hyperslab
  hsize_t file_offset[2];
  hsize_t file_count[2];
  file_offset[0] = Offset;
  file_offset[1] = 0;
  file_count[0] = nof_particles;
  file_count[1] = 3;

  status = my_H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_offset, NULL, file_count, NULL);

  //velocity memory hyperslab
  hsize_t mem[1];
  mem[0] = nof_particles * 3;

  hid_t memspace = my_H5Screate_simple(1, mem, NULL);
  hsize_t mem_offset[1];
  hsize_t mem_count[1];
  mem_offset[0] = 0;
  mem_count[0] = nof_particles * 3;

  status = my_H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_offset, NULL, mem_count, NULL);

  //write
  status = my_H5Dwrite(dataset, H5T_NATIVE_FLOAT, memspace, filespace, H5P_DEFAULT, buffer, dataset_name);

  myassert(status >= 0);

  //close handles
  status = my_H5Sclose(memspace, H5S_SIMPLE);
  status = my_H5Sclose(filespace, H5S_SIMPLE);
  status = my_H5Dclose(dataset, dataset_name);
}

static void write_float_buffer(hid_t file_id, unsigned int nof_particles, float *buffer, const char *dataset_name)
{
  //identifier
  herr_t status;

  //open dataset and get dataspace
  hid_t dataset = my_H5Dopen(file_id, dataset_name);
  hid_t filespace = my_H5Dget_space(dataset, dataset_name);

  //file hyperslab
  hsize_t file_offset[1];
  hsize_t file_count[1];
  file_offset[0] = Offset;
  file_count[0] = nof_particles;

  status = my_H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_offset, NULL, file_count, NULL);

  //memory hyperslab
  hsize_t mem[1];
  mem[0] = nof_particles;

  hid_t memspace = my_H5Screate_simple(1, mem, NULL);
  hsize_t mem_offset[1];
  hsize_t mem_count[1];
  mem_offset[0] = 0;
  mem_count[0] = nof_particles;

  status = my_H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_offset, NULL, mem_count, NULL);

  //write
  status = my_H5Dwrite(dataset, H5T_NATIVE_FLOAT, memspace, filespace, H5P_DEFAULT, buffer, dataset_name);

  myassert(status >= 0);

  //close handles
  status = my_H5Sclose(memspace, H5S_SIMPLE);
  status = my_H5Sclose(filespace, H5S_SIMPLE);
  status = my_H5Dclose(dataset, dataset_name);

}

static void write_zoneflag_buffer(hid_t file_id, unsigned int nof_particles, int *zoneflag_buffer)
{
  //identifier
  herr_t status;

  //open dataset and get dataspace
  hid_t zoneflag_dataset = my_H5Dopen(file_id, "/Zoneflag");
  hid_t zoneflag_filespace = my_H5Dget_space(zoneflag_dataset, "/Zoneflag");

  //zoneflag file hyperslab
  hsize_t zoneflag_file_offset[1];
  hsize_t zoneflag_file_count[1];
  zoneflag_file_offset[0] = Offset;
  zoneflag_file_count[0] = nof_particles;

  status = my_H5Sselect_hyperslab(zoneflag_filespace, H5S_SELECT_SET, zoneflag_file_offset, NULL, zoneflag_file_count, NULL);

  //zoneflag memory hyperslab
  hsize_t zoneflag_mem[1];
  zoneflag_mem[0] = nof_particles;

  hid_t zoneflag_memspace = my_H5Screate_simple(1, zoneflag_mem, NULL);
  hsize_t zoneflag_mem_offset[1];
  hsize_t zoneflag_mem_count[1];
  zoneflag_mem_offset[0] = 0;
  zoneflag_mem_count[0] = nof_particles;

  status = my_H5Sselect_hyperslab(zoneflag_memspace, H5S_SELECT_SET, zoneflag_mem_offset, NULL, zoneflag_mem_count, NULL);

  //write
  status = my_H5Dwrite(zoneflag_dataset, H5T_NATIVE_INT, zoneflag_memspace, zoneflag_filespace, H5P_DEFAULT, zoneflag_buffer, "/Zoneflag");

  myassert(status >= 0);

  //close handles
  status = my_H5Sclose(zoneflag_memspace, H5S_SIMPLE);
  status = my_H5Sclose(zoneflag_filespace, H5S_SIMPLE);
  status = my_H5Dclose(zoneflag_dataset, "/Zoneflag");
}

static void write_id_buffer(hid_t file_id, unsigned int nof_particles, unsigned long long int *id_buffer)
{
  //identifier
  herr_t status;

  //open dataset and get dataspace
  hid_t id_dataset = my_H5Dopen(file_id, "/ID");
  hid_t id_filespace = my_H5Dget_space(id_dataset, "/ID");

  //id file hyperslab
  hsize_t id_file_offset[1];
  hsize_t id_file_count[1];
  id_file_offset[0] = Offset;
  id_file_count[0] = nof_particles;

  status = my_H5Sselect_hyperslab(id_filespace, H5S_SELECT_SET, id_file_offset, NULL, id_file_count, NULL);

  //id memory hyperslab
  hsize_t id_mem[1];
  id_mem[0] = nof_particles;

  hid_t id_memspace = my_H5Screate_simple(1, id_mem, NULL);
  hsize_t id_mem_offset[1];
  hsize_t id_mem_count[1];
  id_mem_offset[0] = 0;
  id_mem_count[0] = nof_particles;

  status = my_H5Sselect_hyperslab(id_memspace, H5S_SELECT_SET, id_mem_offset, NULL, id_mem_count, NULL);

  //write
  status = my_H5Dwrite(id_dataset, H5T_NATIVE_ULLONG, id_memspace, id_filespace, H5P_DEFAULT, id_buffer, "/ID");

  myassert(status >= 0);

  //close handles
  status = my_H5Sclose(id_memspace, H5S_SIMPLE);
  status = my_H5Sclose(id_filespace, H5S_SIMPLE);
  status = my_H5Dclose(id_dataset, "/ID");
}




#define WRITE_BUFFERS(size) \
write_float_buffer_vector(file_id, size, SfVars.coord_buffer, "/Coordinates"); \
write_float_buffer_vector(file_id, size, SfVars.velocity_buffer, "/Velocities"); \
SDIR(write_float_buffer_vector(file_id, size, SfVars.shockdir_buffer, "/ShockDirection");) \
write_float_buffer_vector(file_id, size, SfVars.vpre_buffer, "/PreShockVelocity"); \
write_float_buffer_vector(file_id, size, SfVars.vpost_buffer, "/PostShockVelocity"); \
FANGLE(write_float_buffer(file_id, size, SfVars.fangle_buffer, "/MaxFaceAngle");) \
write_float_buffer(file_id, size, SfVars.machnum_buffer, "/Machnumber"); \
write_float_buffer(file_id, size, SfVars.surface_buffer, "/Surface"); \
write_float_buffer(file_id, size, SfVars.volume_buffer, "/Volume"); \
write_float_buffer(file_id, size, SfVars.density_buffer, "/Density"); \
write_float_buffer(file_id, size, SfVars.genUflux_buffer, "/GeneratedInternalEnergyFlux"); \
write_float_buffer(file_id, size, SfVars.cpre_buffer, "/PreShockSoundSpeed"); \
write_float_buffer(file_id, size, SfVars.rhopre_buffer, "/PreShockDensity"); \
write_float_buffer(file_id, size, SfVars.ppre_buffer, "/PreShockPressure"); \
write_float_buffer(file_id, size, SfVars.rhopost_buffer, "/PostShockDensity"); \
write_float_buffer(file_id, size, SfVars.ppost_buffer, "/PostShockPressure"); \
write_float_buffer(file_id, size, SfVars.tpre_buffer, "/PreShockTemperature"); \
write_float_buffer(file_id, size, SfVars.intenergy_buffer, "/InternalEnergy"); \
write_float_buffer(file_id, size, SfVars.temperature_buffer, "/Temperature"); \
write_zoneflag_buffer(file_id, size, SfVars.zoneflag_buffer); \
write_id_buffer(file_id, size, SfVars.id_buffer);

#define SEND_BUFFERS(size) \
MPI_Send(SfVars.coord_buffer, size * 3, MPI_FLOAT, q, 0, MPI_COMM_WORLD); \
MPI_Send(SfVars.velocity_buffer, size * 3, MPI_FLOAT, q, 0, MPI_COMM_WORLD); \
SDIR(MPI_Send(SfVars.shockdir_buffer, size * 3, MPI_FLOAT, q, 0, MPI_COMM_WORLD);) \
MPI_Send(SfVars.vpre_buffer, size * 3, MPI_FLOAT, q, 0, MPI_COMM_WORLD); \
MPI_Send(SfVars.vpost_buffer, size * 3, MPI_FLOAT, q, 0, MPI_COMM_WORLD); \
FANGLE(MPI_Send(SfVars.fangle_buffer, size, MPI_FLOAT, q, 0, MPI_COMM_WORLD);) \
MPI_Send(SfVars.machnum_buffer, size, MPI_FLOAT, q, 0, MPI_COMM_WORLD); \
MPI_Send(SfVars.surface_buffer, size, MPI_FLOAT, q, 0, MPI_COMM_WORLD); \
MPI_Send(SfVars.volume_buffer, size, MPI_FLOAT, q, 0, MPI_COMM_WORLD); \
MPI_Send(SfVars.density_buffer, size, MPI_FLOAT, q, 0, MPI_COMM_WORLD); \
MPI_Send(SfVars.genUflux_buffer, size, MPI_FLOAT, q, 0, MPI_COMM_WORLD); \
MPI_Send(SfVars.cpre_buffer, size, MPI_FLOAT, q, 0, MPI_COMM_WORLD); \
MPI_Send(SfVars.rhopre_buffer, size, MPI_FLOAT, q, 0, MPI_COMM_WORLD); \
MPI_Send(SfVars.ppre_buffer, size, MPI_FLOAT, q, 0, MPI_COMM_WORLD); \
MPI_Send(SfVars.rhopost_buffer, size, MPI_FLOAT, q, 0, MPI_COMM_WORLD); \
MPI_Send(SfVars.ppost_buffer, size, MPI_FLOAT, q, 0, MPI_COMM_WORLD); \
MPI_Send(SfVars.tpre_buffer, size, MPI_FLOAT, q, 0, MPI_COMM_WORLD); \
MPI_Send(SfVars.intenergy_buffer, size, MPI_FLOAT, q, 0, MPI_COMM_WORLD); \
MPI_Send(SfVars.temperature_buffer, size, MPI_FLOAT, q, 0, MPI_COMM_WORLD); \
MPI_Send(SfVars.zoneflag_buffer, size, MPI_INT, q, 0, MPI_COMM_WORLD); \
MPI_Send(SfVars.id_buffer, Buffer_size, MPI_UNSIGNED_LONG_LONG, q, 0, MPI_COMM_WORLD);


#define RECEIVE_BUFFERS(status) \
MPI_Recv(SfVars.coord_buffer, Buffer_size * 3, MPI_FLOAT, n, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); \
MPI_Recv(SfVars.velocity_buffer, Buffer_size * 3, MPI_FLOAT, n, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); \
SDIR(MPI_Recv(SfVars.shockdir_buffer, Buffer_size * 3, MPI_FLOAT, n, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);) \
MPI_Recv(SfVars.vpre_buffer, Buffer_size * 3, MPI_FLOAT, n, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); \
MPI_Recv(SfVars.vpost_buffer, Buffer_size * 3, MPI_FLOAT, n, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); \
FANGLE(MPI_Recv(SfVars.fangle_buffer, Buffer_size, MPI_FLOAT, n, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);) \
MPI_Recv(SfVars.machnum_buffer, Buffer_size, MPI_FLOAT, n, 0, MPI_COMM_WORLD, status); \
MPI_Recv(SfVars.surface_buffer, Buffer_size, MPI_FLOAT, n, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); \
MPI_Recv(SfVars.volume_buffer, Buffer_size, MPI_FLOAT, n, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); \
MPI_Recv(SfVars.density_buffer, Buffer_size, MPI_FLOAT, n, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); \
MPI_Recv(SfVars.genUflux_buffer, Buffer_size, MPI_FLOAT, n, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); \
MPI_Recv(SfVars.cpre_buffer, Buffer_size, MPI_FLOAT, n, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); \
MPI_Recv(SfVars.rhopre_buffer, Buffer_size, MPI_FLOAT, n, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); \
MPI_Recv(SfVars.ppre_buffer, Buffer_size, MPI_FLOAT, n, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); \
MPI_Recv(SfVars.rhopost_buffer, Buffer_size, MPI_FLOAT, n, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); \
MPI_Recv(SfVars.ppost_buffer, Buffer_size, MPI_FLOAT, n, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); \
MPI_Recv(SfVars.tpre_buffer, Buffer_size, MPI_FLOAT, n, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); \
MPI_Recv(SfVars.intenergy_buffer, Buffer_size, MPI_FLOAT, n, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); \
MPI_Recv(SfVars.temperature_buffer, Buffer_size, MPI_FLOAT, n, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); \
MPI_Recv(SfVars.zoneflag_buffer, Buffer_size, MPI_INT, n, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); \
MPI_Recv(SfVars.id_buffer, Buffer_size, MPI_UNSIGNED_LONG_LONG, n, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


//!send data from task n to task q and write to disk
/*!
 * \param file_id The file handle
 * \param n The sending task
 * \param q The receiving/writing task
 */
static void send_write_data(hid_t file_id, unsigned int n, unsigned int q)
{
  set_cosmo_factors_for_current_time();

  if(!(ThisTask == q || ThisTask == n)) //only task q and n should enter this function
    {
      return;
    }

  unsigned int buffer_counter = 0;
  unsigned int nof_sends;       //number of sends

  if(ThisTask != q || n == q)   //n==q: task n, only write!
    {
      if(n != q)
        {
          //send to task q the number of sends
          nof_sends = ceil((double) NumGas / Buffer_size);

          MPI_Send(&nof_sends, 1, MPI_UNSIGNED, q, 0, MPI_COMM_WORLD);
        }


      unsigned int i;

      for(i = 0; i < NumGas; i++)
        {
          SfVars.coord_buffer[buffer_counter * 3] = P[i].Pos[0];
          SfVars.coord_buffer[buffer_counter * 3 + 1] = P[i].Pos[1];
          SfVars.coord_buffer[buffer_counter * 3 + 2] = P[i].Pos[2];

          SfVars.velocity_buffer[buffer_counter * 3] = P[i].Vel[0] * sqrt(All.cf_a3inv);
          SfVars.velocity_buffer[buffer_counter * 3 + 1] = P[i].Vel[1] * sqrt(All.cf_a3inv);
          SfVars.velocity_buffer[buffer_counter * 3 + 2] = P[i].Vel[2] * sqrt(All.cf_a3inv);

          SDIR(SfVars.shockdir_buffer[buffer_counter * 3] = SDATA(i).ShockDir[0];
          SfVars.shockdir_buffer[buffer_counter * 3 + 1] = SDATA(i).ShockDir[1];
          SfVars.shockdir_buffer[buffer_counter * 3 + 2] = SDATA(i).ShockDir[2];)

          SfVars.vpre_buffer[buffer_counter * 3] = SDATA(i).VpreShock[0]*sqrt(All.cf_a3inv);
          SfVars.vpre_buffer[buffer_counter * 3 + 1] = SDATA(i).VpreShock[1]*sqrt(All.cf_a3inv);
          SfVars.vpre_buffer[buffer_counter * 3 + 2] = SDATA(i).VpreShock[2]*sqrt(All.cf_a3inv);

          SfVars.vpost_buffer[buffer_counter * 3] = SDATA(i).VpostShock[0]*sqrt(All.cf_a3inv);
          SfVars.vpost_buffer[buffer_counter * 3 + 1] = SDATA(i).VpostShock[1]*sqrt(All.cf_a3inv);
          SfVars.vpost_buffer[buffer_counter * 3 + 2] = SDATA(i).VpostShock[2]*sqrt(All.cf_a3inv);

          FANGLE(SfVars.fangle_buffer[buffer_counter] = SphP[i].MaxFaceAngle;)
          SfVars.machnum_buffer[buffer_counter] = SphP[i].Machnumber;
          SfVars.surface_buffer[buffer_counter] = SDATA(i).ShockSurfaceArea;
          SfVars.volume_buffer[buffer_counter] = SphP[i].Volume;
          SfVars.density_buffer[buffer_counter] = SphP[i].Density;
          SfVars.genUflux_buffer[buffer_counter] = generated_thermal_energy_flux(i);
          SfVars.cpre_buffer[buffer_counter] = SDATA(i).CpreShock;
          SfVars.rhopre_buffer[buffer_counter] = SDATA(i).RhopreShock;
          SfVars.ppre_buffer[buffer_counter] = SDATA(i).PpreShock;
          SfVars.rhopost_buffer[buffer_counter] = SDATA(i).RhopostShock;
          SfVars.ppost_buffer[buffer_counter] = SDATA(i).PpostShock;
          SfVars.tpre_buffer[buffer_counter] = SDATA(i).TpreShock;
          SfVars.intenergy_buffer[buffer_counter] = SphP[i].Utherm;
          SfVars.temperature_buffer[buffer_counter] = SDATA(i).Temperature;
          SfVars.id_buffer[buffer_counter] = (unsigned long long int) P[i].ID;

          if(SDATA(i).ShockSurface)
            {
              SfVars.zoneflag_buffer[buffer_counter] = 2;
            }
          else if(SDATA(i).ShockZone)
            {
              SfVars.zoneflag_buffer[buffer_counter] = 1;
            }
          else
            {
              SfVars.zoneflag_buffer[buffer_counter] = 0;
            }

          buffer_counter++;

          if(buffer_counter == Buffer_size)     //write buffer or send buffer
            {
              if(n == q)
                {
                  WRITE_BUFFERS(Buffer_size)

                  Offset += Buffer_size;
                }
              else
                {
                  SEND_BUFFERS(Buffer_size)
                }

                  buffer_counter = 0;
            }
        }

      //send/write remaining partial buffer
      if(buffer_counter != 0)
        {
          if(n == q)
            {
              WRITE_BUFFERS(buffer_counter)

              Offset += buffer_counter;

              buffer_counter = 0;
            }
          else
            {
              SEND_BUFFERS(buffer_counter)
            }
        }
    }
  else                          //task q
    {
      MPI_Recv(&nof_sends, 1, MPI_UNSIGNED, n, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      unsigned int k;

      for(k = 0; k < nof_sends - 1; k++)
        {
          RECEIVE_BUFFERS(MPI_STATUS_IGNORE)

          WRITE_BUFFERS(Buffer_size)

          Offset += Buffer_size;
        }

      //receive remaining partial buffer
      MPI_Status status;
      int count;

      RECEIVE_BUFFERS(&status)

      MPI_Get_count(&status, MPI_FLOAT, &count);

      WRITE_BUFFERS(count)

      Offset += count;
    }
}


//!get for each task the number of particles which have to be written to disk
static int get_particle_numbers()
{
  int tasks_group = NTask / All.NumFilesPerOutput;

  if(ThisTask % tasks_group == 0)       //its a write task
    {
      int result = NumGas;

      unsigned int k;
      int buffer;

      for(k = 1; k < tasks_group; k++)
        {
          MPI_Recv(&buffer, 1, MPI_INT, ThisTask + k, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          result += buffer;
        }

      return result;
    }
  else                          //its not a write task
    {

      int write_task = ThisTask - (ThisTask % tasks_group);

      MPI_Send(&NumGas, 1, MPI_INT, write_task, 0, MPI_COMM_WORLD);

      return 0;
    }
}


void shock_finder_output()
{
  //memory reservation
  SfVars.coord_buffer = (float *) mymalloc("coord_buffer", 3 * Buffer_size * sizeof(float));
  SfVars.velocity_buffer = (float *) mymalloc("velocity_buffer", 3 * Buffer_size * sizeof(float));
  SDIR(SfVars.shockdir_buffer = (float *) mymalloc("shockdir_buffer", 3 * Buffer_size * sizeof(float));)
  SfVars.vpre_buffer = (float *) mymalloc("vpre_buffer", 3 * Buffer_size * sizeof(float));
  SfVars.vpost_buffer = (float *) mymalloc("vpost_buffer", 3 * Buffer_size * sizeof(float));
  FANGLE(SfVars.fangle_buffer = (float *) mymalloc("fangle_buffer", Buffer_size * sizeof(float));)
  SfVars.machnum_buffer = (float *) mymalloc("machnum_buffer", Buffer_size * sizeof(float));
  SfVars.surface_buffer = (float *) mymalloc("surface_buffer", Buffer_size * sizeof(float));
  SfVars.volume_buffer = (float *) mymalloc("volume_buffer", Buffer_size * sizeof(float));
  SfVars.density_buffer = (float *) mymalloc("density_buffer", Buffer_size * sizeof(float));
  SfVars.genUflux_buffer = (float *) mymalloc("genUflux_buffer", Buffer_size * sizeof(float));
  SfVars.cpre_buffer = (float *) mymalloc("cpre_buffer", Buffer_size * sizeof(float));
  SfVars.rhopre_buffer = (float *) mymalloc("rhopre_buffer", Buffer_size * sizeof(float));
  SfVars.ppre_buffer = (float *) mymalloc("ppre_buffer", Buffer_size * sizeof(float));
  SfVars.rhopost_buffer = (float *) mymalloc("rhopost_buffer", Buffer_size * sizeof(float));
  SfVars.ppost_buffer = (float *) mymalloc("ppost_buffer", Buffer_size * sizeof(float));
  SfVars.tpre_buffer = (float *) mymalloc("tpre_buffer", Buffer_size * sizeof(float));
  SfVars.intenergy_buffer = (float *) mymalloc("intenergy_buffer", Buffer_size * sizeof(float));
  SfVars.temperature_buffer = (float *) mymalloc("temperature_buffer", Buffer_size * sizeof(float));
  SfVars.zoneflag_buffer = (int *) mymalloc("zoneflag_buffer", Buffer_size * sizeof(int));
  SfVars.id_buffer = (unsigned long long int *) mymalloc("id_buffer", Buffer_size * sizeof(unsigned long long int));

  //identifiers
  hid_t file_id = 0;
  herr_t status;

  int num_particles;
  //create a new hdf5 file
  char filename[1024];

  //number of particles which have to be written to disk (for this task)
  num_particles = get_particle_numbers();

  int tasks_group = NTask / All.NumFilesPerOutput;

  if(ThisTask % tasks_group == 0)       // it's a write task
    {
      if(All.NumFilesPerOutput == 1)
        {
#ifdef SHOCK_SUBVOLUME
          sprintf(filename, "%s/shocks%d_%03d.hdf5", All.OutputDir, All.ShockSubvolumeNum, RestartSnapNum);
#else
          sprintf(filename, "%s/shocks_%03d.hdf5", All.OutputDir, RestartSnapNum);
#endif
        }
      else if(NTask == All.NumFilesPerOutput)
        {
          int filenum = ThisTask;
          sprintf(filename, "%s/shocks_%03d_%03d.hdf5", All.OutputDir, RestartSnapNum, filenum);
        }
      else
        {
          int filenum = ThisTask / tasks_group;
          sprintf(filename, "%s/shocks_%03d_%03d.hdf5", All.OutputDir, RestartSnapNum, filenum);
        }


      file_id = my_H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      hid_t dataspace_id, dataset_id, attribute_id, group_id;
      hsize_t dims1d[1];
      hsize_t dims2d[2];

      dims1d[0] = num_particles;
      dims2d[0] = num_particles;
      dims2d[1] = 3;


      /*
       *  Header
       */


      group_id = my_H5Gcreate(file_id, "/Header", 0);

      //Number of files
      dataspace_id = my_H5Screate(H5S_SCALAR);
      attribute_id = my_H5Acreate(group_id, "NumFiles", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT);
      status = my_H5Awrite(attribute_id, H5T_NATIVE_INT, &All.NumFilesPerOutput, "NumFiles");
      status = my_H5Aclose(attribute_id, "NumFiles");
      status = my_H5Sclose(dataspace_id, H5S_SCALAR);

      //Time
      dataspace_id = my_H5Screate(H5S_SCALAR);
      attribute_id = my_H5Acreate(group_id, "Time", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
      status = my_H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &All.Time, "Time");
      status = my_H5Aclose(attribute_id, "Time");
      status = my_H5Sclose(dataspace_id, H5S_SCALAR);

      //Redshift
      double redshift = 1.0 / All.Time - 1;
      dataspace_id = my_H5Screate(H5S_SCALAR);
      attribute_id = my_H5Acreate(group_id, "Redshift", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
      status = my_H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &redshift, "Redshift");
      status = my_H5Aclose(attribute_id, "Redshift");
      status = my_H5Sclose(dataspace_id, H5S_SCALAR);

      //Boxsize
      dataspace_id = my_H5Screate(H5S_SCALAR);
      attribute_id = my_H5Acreate(group_id, "BoxSize", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
      status = my_H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &All.BoxSize, "BoxSize");
      status = my_H5Aclose(attribute_id, "BoxSize");
      status = my_H5Sclose(dataspace_id, H5S_SCALAR);

      //number of gas cells
      dataspace_id = my_H5Screate(H5S_SCALAR);
      attribute_id = my_H5Acreate(group_id, "NumGasParticles", H5T_NATIVE_LLONG, dataspace_id, H5P_DEFAULT);
      status = my_H5Awrite(attribute_id, H5T_NATIVE_LLONG, &SfVars.nof_gas_part_total, "NumGasParticles");
      status = my_H5Aclose(attribute_id, "NumGasParticles");
      status = my_H5Sclose(dataspace_id, H5S_SCALAR);

      //number of gas cells in this hdf5 file
      dataspace_id = my_H5Screate(H5S_SCALAR);
      attribute_id = my_H5Acreate(group_id, "NumGasParticlesThisFile", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT);
      status = my_H5Awrite(attribute_id, H5T_NATIVE_INT, &num_particles, "NumGasParticlesThisFile");
      status = my_H5Aclose(attribute_id, "NumGasParticlesThisFile");
      status = my_H5Sclose(dataspace_id, H5S_SCALAR);

      //number of cells in the shock zone
      dataspace_id = my_H5Screate(H5S_SCALAR);
      attribute_id = my_H5Acreate(group_id, "NumShockZoneCells", H5T_NATIVE_LLONG, dataspace_id, H5P_DEFAULT);
      status = my_H5Awrite(attribute_id, H5T_NATIVE_LLONG, &SfVars.nof_shockzone_part_total, "NumShockZoneCells");
      status = my_H5Aclose(attribute_id, "NumShockZoneCells");
      status = my_H5Sclose(dataspace_id, H5S_SCALAR);

      //number of cells in the shock surface
      dataspace_id = my_H5Screate(H5S_SCALAR);
      attribute_id = my_H5Acreate(group_id, "NumShockSurfaceCells", H5T_NATIVE_LLONG, dataspace_id, H5P_DEFAULT);
      status = my_H5Awrite(attribute_id, H5T_NATIVE_LLONG, &SfVars.nof_shocksurface_part_total, "NumShockSurfaceCells");
      status = my_H5Aclose(attribute_id, "NumShockSurfaceCells");
      status = my_H5Sclose(dataspace_id, H5S_SCALAR);

      //Omega 0
      dataspace_id = my_H5Screate(H5S_SCALAR);
      attribute_id = my_H5Acreate(group_id, "Omega0", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
      status = my_H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &All.Omega0, "Omega0");
      status = my_H5Aclose(attribute_id, "Omega0");
      status = my_H5Sclose(dataspace_id, H5S_SCALAR);

      //Omega Lambda
      dataspace_id = my_H5Screate(H5S_SCALAR);
      attribute_id = my_H5Acreate(group_id, "OmegaLambda", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
      status = my_H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &All.OmegaLambda, "OmegaLambda");
      status = my_H5Aclose(attribute_id, "OmegaLambda");
      status = my_H5Sclose(dataspace_id, H5S_SCALAR);

      //Hubble parameter
      dataspace_id = my_H5Screate(H5S_SCALAR);
      attribute_id = my_H5Acreate(group_id, "HubbleParam", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
      status = my_H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &All.HubbleParam, "HubbleParam");
      status = my_H5Aclose(attribute_id, "HubbleParam");
      status = my_H5Sclose(dataspace_id, H5S_SCALAR);


      //close header
      status = my_H5Gclose(group_id, "/Header");


      /*
       * Coordinate dataset
       */

      dataspace_id = my_H5Screate_simple(2, dims2d, NULL);
      dataset_id = my_H5Dcreate(file_id, "/Coordinates", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);
      status = my_H5Dclose(dataset_id, "/Coordinates");
      status = my_H5Sclose(dataspace_id, H5S_SIMPLE);

      /*
       * velocity dataset
       */

      dataspace_id = my_H5Screate_simple(2, dims2d, NULL);
      dataset_id = my_H5Dcreate(file_id, "/Velocities", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);
      status = my_H5Dclose(dataset_id, "/Velocities");
      status = my_H5Sclose(dataspace_id, H5S_SIMPLE);

      /*
       * shock direction dataset
       */

      SDIR(dataspace_id = my_H5Screate_simple(2, dims2d, NULL);
      dataset_id = my_H5Dcreate(file_id, "/ShockDirection", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);
      status = my_H5Dclose(dataset_id, "/ShockDirection");
      status = my_H5Sclose(dataspace_id, H5S_SIMPLE);)

      /*
       * pre-shock velocity dataset
       */

      dataspace_id = my_H5Screate_simple(2, dims2d, NULL);
      dataset_id = my_H5Dcreate(file_id, "/PreShockVelocity", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);
      status = my_H5Dclose(dataset_id, "/PreShockVelocity");
      status = my_H5Sclose(dataspace_id, H5S_SIMPLE);

      /*
       * post-shock velocity dataset
       */

      dataspace_id = my_H5Screate_simple(2, dims2d, NULL);
      dataset_id = my_H5Dcreate(file_id, "/PostShockVelocity", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);
      status = my_H5Dclose(dataset_id, "/PostShockVelocity");
      status = my_H5Sclose(dataspace_id, H5S_SIMPLE);

      /*
       * Maximum face angle dataset
       */

      FANGLE(dataspace_id = my_H5Screate_simple(1, dims1d, NULL);
      dataset_id = my_H5Dcreate(file_id, "/MaxFaceAngle", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);
      status = my_H5Dclose(dataset_id, "/MaxFaceAngle");
      status = my_H5Sclose(dataspace_id, H5S_SIMPLE);)

      /*
       * Machnumber dataset
       */

      dataspace_id = my_H5Screate_simple(1, dims1d, NULL);
      dataset_id = my_H5Dcreate(file_id, "/Machnumber", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);
      status = my_H5Dclose(dataset_id, "/Machnumber");
      status = my_H5Sclose(dataspace_id, H5S_SIMPLE);

      /*
       * Surface dataset
       */

      dataspace_id = my_H5Screate_simple(1, dims1d, NULL);
      dataset_id = my_H5Dcreate(file_id, "/Surface", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);
      status = my_H5Dclose(dataset_id, "/Surface");
      status = my_H5Sclose(dataspace_id, H5S_SIMPLE);

      /*
       * Volume dataset
       */

      dataspace_id = my_H5Screate_simple(1, dims1d, NULL);
      dataset_id = my_H5Dcreate(file_id, "/Volume", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);
      status = my_H5Dclose(dataset_id, "/Volume");
      status = my_H5Sclose(dataspace_id, H5S_SIMPLE);

      /*
       * Density dataset
       */

      dataspace_id = my_H5Screate_simple(1, dims1d, NULL);
      dataset_id = my_H5Dcreate(file_id, "/Density", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);
      status = my_H5Dclose(dataset_id, "/Density");
      status = my_H5Sclose(dataspace_id, H5S_SIMPLE);

      /*
       * generated internal energy dataset
       */

      dataspace_id = my_H5Screate_simple(1, dims1d, NULL);
      dataset_id = my_H5Dcreate(file_id, "/GeneratedInternalEnergyFlux", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);
      status = my_H5Dclose(dataset_id, "/GeneratedInternalEnergyFlux");
      status = my_H5Sclose(dataspace_id, H5S_SIMPLE);

      /*
       * pre-shock sound speed dataset
       */

      dataspace_id = my_H5Screate_simple(1, dims1d, NULL);
      dataset_id = my_H5Dcreate(file_id, "/PreShockSoundSpeed", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);
      status = my_H5Dclose(dataset_id, "/PreShockSoundSpeed");
      status = my_H5Sclose(dataspace_id, H5S_SIMPLE);

      /*
       * pre-shock density dataset
       */

      dataspace_id = my_H5Screate_simple(1, dims1d, NULL);
      dataset_id = my_H5Dcreate(file_id, "/PreShockDensity", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);
      status = my_H5Dclose(dataset_id, "/PreShockDensity");
      status = my_H5Sclose(dataspace_id, H5S_SIMPLE);

      /*
       * pre-shock pressure dataset
       */

      dataspace_id = my_H5Screate_simple(1, dims1d, NULL);
      dataset_id = my_H5Dcreate(file_id, "/PreShockPressure", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);
      status = my_H5Dclose(dataset_id, "/PreShockPressure");
      status = my_H5Sclose(dataspace_id, H5S_SIMPLE);

      /*
       * post-shock density dataset
       */

      dataspace_id = my_H5Screate_simple(1, dims1d, NULL);
      dataset_id = my_H5Dcreate(file_id, "/PostShockDensity", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);
      status = my_H5Dclose(dataset_id, "/PostShockDensity");
      status = my_H5Sclose(dataspace_id, H5S_SIMPLE);

      /*
       * post-shock pressure dataset
       */

      dataspace_id = my_H5Screate_simple(1, dims1d, NULL);
      dataset_id = my_H5Dcreate(file_id, "/PostShockPressure", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);
      status = my_H5Dclose(dataset_id, "/PostShockPressure");
      status = my_H5Sclose(dataspace_id, H5S_SIMPLE);

      /*
       * pre-shock temperature dataset
       */

      dataspace_id = my_H5Screate_simple(1, dims1d, NULL);
      dataset_id = my_H5Dcreate(file_id, "/PreShockTemperature", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);
      status = my_H5Dclose(dataset_id, "/PreShockTemperature");
      status = my_H5Sclose(dataspace_id, H5S_SIMPLE);

      /*
       * Internal Energy dataset
       */

      dataspace_id = my_H5Screate_simple(1, dims1d, NULL);
      dataset_id = my_H5Dcreate(file_id, "/InternalEnergy", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);
      status = my_H5Dclose(dataset_id, "/InternalEnergy");
      status = my_H5Sclose(dataspace_id, H5S_SIMPLE);

      /*
       * Temperature dataset
       */

      dataspace_id = my_H5Screate_simple(1, dims1d, NULL);
      dataset_id = my_H5Dcreate(file_id, "/Temperature", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);
      status = my_H5Dclose(dataset_id, "/Temperature");
      status = my_H5Sclose(dataspace_id, H5S_SIMPLE);


      /*
       * Zoneflag dataset
       */

      dataspace_id = my_H5Screate_simple(1, dims1d, NULL);
      dataset_id = my_H5Dcreate(file_id, "/Zoneflag", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT);
      status = my_H5Dclose(dataset_id, "/Zoneflag");
      status = my_H5Sclose(dataspace_id, H5S_SIMPLE);

      /*
       * ID dataset
       */


      dataspace_id = my_H5Screate_simple(1, dims1d, NULL);
      dataset_id = my_H5Dcreate(file_id, "/ID", H5T_NATIVE_ULLONG, dataspace_id, H5P_DEFAULT);
      status = my_H5Dclose(dataset_id, "/ID");
      status = my_H5Sclose(dataspace_id, H5S_SIMPLE);


      myassert(status >= 0);
    }

  int tasks_per_group = NTask / All.NumFilesPerOutput;

  unsigned int k;
  unsigned int l;
  unsigned int m;

  for(k = 0; k < All.NumFilesPerOutput; k++)
    {
      for(l = 0; l < tasks_per_group; l++)
        {
          m = k * tasks_per_group + l;
          send_write_data(file_id, m, k * tasks_per_group);
        }
    }

  //close hdf5 file handle
  if(ThisTask % (NTask / All.NumFilesPerOutput) == 0)
    {
      status = my_H5Fclose(file_id, filename);
    }

  //free memory
  myfree(SfVars.id_buffer);
  myfree(SfVars.zoneflag_buffer);
  myfree(SfVars.temperature_buffer);
  myfree(SfVars.intenergy_buffer);
  myfree(SfVars.tpre_buffer);
  myfree(SfVars.ppost_buffer);
  myfree(SfVars.rhopost_buffer);
  myfree(SfVars.ppre_buffer);
  myfree(SfVars.rhopre_buffer);
  myfree(SfVars.cpre_buffer);
  myfree(SfVars.genUflux_buffer);
  myfree(SfVars.density_buffer);
  myfree(SfVars.volume_buffer);
  myfree(SfVars.surface_buffer);
  myfree(SfVars.machnum_buffer);
  FANGLE(myfree(SfVars.fangle_buffer);)
  myfree(SfVars.vpost_buffer);
  myfree(SfVars.vpre_buffer);
  SDIR(myfree(SfVars.shockdir_buffer);)
  myfree(SfVars.velocity_buffer);
  myfree(SfVars.coord_buffer);
}

#endif



#endif
