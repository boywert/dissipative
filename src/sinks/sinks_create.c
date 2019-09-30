/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/sinks/sinks_create.c
 * \date        01/2013
 * \author      Thomas Greif
 * \brief       Sink particles
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 */

#include "../allvars.h"
#include "../proto.h"

void sinks_find_candidate(void);
void sinks_get_accreted_particles(void);
void sinks_get_center_of_mass(void);
void sinks_add_sink(void);
int sinks_energycheck(void);


void sinks_create(void)
{

  sinks_find_candidate();

  mpi_printf("SINKS: task %d SKD task %d SKD index %d SKD ID %d SKD NHMax %g SKD NThresh %g SKD pos %g|%g|%g\n", ThisTask, SKD.Task, SKD.Index, SKD.ID, SKD.NHMax, SKD.NHThresh, SKD.Pos[0], SKD.Pos[1], SKD.Pos[2]);

  if(SKD.NHMax > SKD.NHThresh)
    {
      sinks_get_accreted_particles();

#ifdef SGCHEM
      if(sinks_energycheck() > 0)
        {
          sinks_add_sink();

          mpi_printf("SINKS: Created a new sink.\n");
        }

#else

      sinks_get_center_of_mass();

      sinks_add_sink();

      mpi_printf("SINKS: Created a new sink.\n");

#endif
    }
}


void sinks_find_candidate(void)
{
  int idx, i, j, flag, totflag;
  double r2_min, nh;

  struct
  {
    double nh_max;
    int task;

  } local, global;


  flag = 0;

  local.nh_max = 0.0;
  r2_min = 0.0;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

        nh = SKD.NHFac * SphP[i].Density;

        if(flag == 0 || nh > local.nh_max)
          {
            flag = 1;

            local.nh_max = nh;

            SKD.Index = i;
            SKD.ID = P[i].ID;

            for(j = 0; j < 3; j++)
              SKD.Pos[j] = P[i].Pos[j];
          }
    }

  MPI_Allreduce(&flag, &totflag, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(totflag)
    {
      local.task = ThisTask;

      MPI_Allreduce(&local, &global, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

      SKD.Task = global.task;
      SKD.NHMax = global.nh_max;

      MPI_Bcast(&SKD.Index, 1, MPI_INT, SKD.Task, MPI_COMM_WORLD);
      MPI_Bcast(&SKD.ID, 1, MPI_INT, SKD.Task, MPI_COMM_WORLD);
      MPI_Bcast(SKD.Pos, 3, MPI_DOUBLE, SKD.Task, MPI_COMM_WORLD);
    }
  else
    SKD.NHMax = 0;
}


void sinks_get_accreted_particles(void)
{
  int idx, i, j;
  double dx, dy, dz, xtmp, ytmp, ztmp, r2;

  SKD.NumAcc = 0;
  SKD.MaxNumAcc = SKD_INIT_MAX_NUM_ACC;

  SKD.ACC = mymalloc_movable(&SKD.ACC, "SKD.ACC", SKD.MaxNumAcc * sizeof(struct ACC_struct));

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].ID == 0 && P[i].Mass == 0)
        continue;

      dx = NGB_PERIODIC_LONG_X(P[i].Pos[0] - SKD.Pos[0]);
      dy = NGB_PERIODIC_LONG_Y(P[i].Pos[1] - SKD.Pos[1]);
      dz = NGB_PERIODIC_LONG_Z(P[i].Pos[2] - SKD.Pos[2]);

      r2 = SKD.DistFac * (dx * dx + dy * dy + dz * dz);

      if(r2 < SKD.AccRad2 && TimeBinSynchronized[P[i].TimeBinHydro] > 0 && P[i].TimeBinHydro == All.LowestActiveTimeBin)
        {
          if(SKD.NumAcc >= SKD.MaxNumAcc)
            {
              SKD.MaxNumAcc = ALLOC_INCREASE_FACTOR * SKD.MaxNumAcc + SKD_INIT_MAX_NUM_ACC;

              SKD.ACC = myrealloc_movable(SKD.ACC, SKD.MaxNumAcc * sizeof(struct ACC_struct));
            }

          SKD.ACC[SKD.NumAcc].Task = SKD.Task;
          SKD.ACC[SKD.NumAcc].Index = i;
          SKD.ACC[SKD.NumAcc].Active = TimeBinSynchronized[P[i].TimeBinHydro];

          for(j = 0; j < 3; j++)
            {
              SKD.ACC[SKD.NumAcc].Pos[j] = P[i].Pos[j];
              SKD.ACC[SKD.NumAcc].Vel[j] = P[i].Vel[j];
            }

          SKD.ACC[SKD.NumAcc].Mass = P[i].Mass;
#ifdef SGCHEM
          SKD.ACC[SKD.NumAcc].Utherm = SphP[i].Utherm;
          for(j = 0; j < 3; j++)
            SKD.ACC[SKD.NumAcc].Accel[j] = P[i].GravAccel[j];
#endif

#ifndef SGCHEM

          P[i].Mass = 0;
          P[i].ID = 0;

#ifdef VORONOI_DYNAMIC_UPDATE
          if(P[i].Type == 0)
            voronoi_remove_connection(i);
#endif
          timebin_remove_particle(&TimeBinsHydro, idx, P[i].TimeBinHydro);
#endif
          SKD.NumAcc++;
        }
    }

#ifndef SGCHEM
  timebin_cleanup_list_of_active_particles(&TimeBinsGravity);
#endif

  mpi_distribute_items_to_tasks(SKD.ACC, offsetof(struct ACC_struct, Task), &SKD.NumAcc, &SKD.MaxNumAcc, sizeof(struct ACC_struct), TAG_DENS_A);

  if(ThisTask == SKD.Task)
    printf("acc: %d %d\n", SKD.NumAcc, SKD.MaxNumAcc);
}


void sinks_get_center_of_mass(void)
{
  int i, j;

  if(ThisTask == SKD.Task)
    {
      for(i = 0; i < 3; i++)
        SKD.CoMPos[i] = SKD.CoMVel[i] = 0;

      SKD.Mass = 0;

      for(i = 0; i < SKD.NumAcc; i++)
        {
          for(j = 0; j < 3; j++)
            {
              SKD.CoMPos[j] += SKD.ACC[i].Mass * SKD.ACC[i].Pos[j];
              SKD.CoMVel[j] += SKD.ACC[i].Mass * SKD.ACC[i].Vel[j];
            }

          SKD.Mass += SKD.ACC[i].Mass;
        }

      for(i = 0; i < 3; i++)
        {
          SKD.CoMPos[i] /= SKD.Mass;
          SKD.CoMVel[i] /= SKD.Mass;
        }
    }

  myfree_movable(SKD.ACC);
}


void sinks_add_sink(void)
{
  int i, igas, isink;

  domain_resize_storage(1, 0, 2);

  if(ThisTask == SKD.Task)
    {
      igas = SKD.Index;
      isink = NumPart;

      if(igas >= isink)
        terminate("This is not possible! isink %d < igas %d", isink, igas);

      P[isink] = P[igas];

      for(i = 0; i < 3; i++)
        {
          P[isink].Pos[i] = SKD.CoMPos[i];
          P[isink].Vel[i] = SKD.CoMVel[i];
        }

      P[isink].Mass = SKD.Mass;
      P[isink].Type = 5;
      P[isink].SofteningType = All.SofteningTypeOfPartType[P[isink].Type];

#ifdef INDIVIDUAL_GRAVITY_SOFTENING
      if(((1 << P[isink].Type) & (INDIVIDUAL_GRAVITY_SOFTENING)))
        P[isink].SofteningType = get_softening_type_from_mass(P[isink].Mass);
#endif

      P[isink].ID = SKD.ID;

      P[isink].TimeBinSink = P[isink].TimeBinHydro;

      timebin_add_particle(&TimeBinsGravity, isink, igas, P[isink].TimeBinGrav, TimeBinSynchronized[P[isink].TimeBinGrav]); // Marinacci, 2016
      timebin_add_particle(&TimeBinsSinksAccretion, isink, -1, P[isink].TimeBinHydro, TimeBinSynchronized[P[isink].TimeBinHydro]);

      NumPart++;

      printf("SINKS: Mass %g ID %d\n", P[isink].Mass, P[isink].ID);
    }

  SKD.TotNumSinks++;
  All.TotNumPart++;
}


#ifdef SGCHEM

int sinks_energycheck(void)
{
  int i, j;
  double e_grav, e_kin, e_therm, divv, diva;
  double dx, dy, dz, xtmp, ytmp, ztmp;
  double r, dvx, dvy, dvz, r2;
  double delx, dely, delz, delr, soft, wp, uu;
  double dax, day, daz;
  double coAcc[3];
  int Etest, Etest_loc;

  Etest = 0;
  Etest_loc = 0;
  soft = get_default_softening_of_particletype(5);

  if(ThisTask == SKD.Task)
    {

      /* find the com */

      for(i = 0; i < 3; i++)
        SKD.CoMPos[i] = SKD.CoMVel[i] = 0;

      SKD.Mass = 0;

      for(i = 0; i < SKD.NumAcc; i++)
        {
          for(j = 0; j < 3; j++)
            {
              SKD.CoMPos[j] += SKD.ACC[i].Mass * SKD.ACC[i].Pos[j];
              SKD.CoMVel[j] += SKD.ACC[i].Mass * SKD.ACC[i].Vel[j];
              coAcc[j] += SKD.ACC[i].Mass * SKD.ACC[i].Accel[j];
            }

          SKD.Mass += SKD.ACC[i].Mass;
        }

      for(i = 0; i < 3; i++)
        {
          SKD.CoMPos[i] /= SKD.Mass;
          SKD.CoMVel[i] /= SKD.Mass;
          coAcc[i] /= SKD.Mass;
        }

      printf("com: %g %g %g %g %g %g %g\n", SKD.CoMPos[0], SKD.CoMPos[1], SKD.CoMPos[2], SKD.CoMVel[0], SKD.CoMVel[1], SKD.CoMVel[2], SKD.Mass);

      /* Check the energies */


      e_grav = e_kin = e_therm = divv = diva = 0.0;



      for(i = 0; i < SKD.NumAcc - 1; i++)
        {
          dx = SKD.ACC[i].Pos[0] - SKD.CoMPos[0];
          dy = SKD.ACC[i].Pos[1] - SKD.CoMPos[1];
          dz = SKD.ACC[i].Pos[2] - SKD.CoMPos[2];
          r = sqrt(dx * dx + dy * dy + dz * dz);
          dvx = SKD.ACC[i].Vel[0] - SKD.CoMVel[0];
          dvy = SKD.ACC[i].Vel[1] - SKD.CoMVel[1];
          dvz = SKD.ACC[i].Vel[2] - SKD.CoMVel[2];
          dax = SKD.ACC[i].Accel[0] - coAcc[0];
          day = SKD.ACC[i].Accel[1] - coAcc[1];
          daz = SKD.ACC[i].Accel[2] - coAcc[2];


          e_kin += 0.5 * SKD.ACC[i].Mass * (dvx * dvx + dvy * dvy + dvz * dvz);

          e_therm += SKD.ACC[i].Utherm * SKD.ACC[i].Mass;


          if(r > 0)
            {
              divv += SKD.ACC[i].Mass * (dvx * dx + dvy * dy + dvz * dz) / r;
              diva += SKD.ACC[i].Mass * (dax * dx + day * dy + daz * dz) / r;
            }

          for(j = i + 1; j < SKD.NumAcc; j++)
            {
              delx = SKD.ACC[i].Pos[0] - SKD.ACC[j].Pos[0];
              dely = SKD.ACC[i].Pos[1] - SKD.ACC[j].Pos[1];
              delz = SKD.ACC[i].Pos[2] - SKD.ACC[j].Pos[2];
              delr = sqrt(delx * delx + dely * dely + delz * delz);


              wp = 1.0 / delr;

              /*
                 if (delr > soft) 
                 wp=1.0/delr;
                 else
                 {
                 uu = delr / soft;

                 if(uu < 0.5)
                 wp = -1.0 / soft * (-2.8 + uu * uu * (5.333333333333 + uu * uu * (6.4 * uu - 9.6)));
                 else
                 wp =-1.0 / soft * (-3.2 + 0.066666666667 / uu + uu * uu * (10.666666666667 + uu * (-16.0 + uu * (9.6 - 2.133333333333 * uu))));
                 } */

              e_grav += All.G * SKD.ACC[i].Mass * SKD.ACC[j].Mass * wp;

            }


        }

      /*printf("Sink Check: %d particles, egrav = %g, ekin = %g, etherm = %g, divv = %g diva = %g\n", SKD.NumAcc, e_grav, e_kin, e_therm, divv, diva); */

      if(e_grav > 2 * e_therm && e_grav > e_kin + e_therm && divv < 0.0 && SKD.NumAcc > 7 && diva < 0.0)
        Etest_loc = 1;
      else
        Etest_loc = 0;
    }

  /* Here we need to broadcast the value of energy check */

  MPI_Allreduce(&Etest_loc, &Etest, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  /*printf("MPI TEST; Task %d  Etest = %d \n",ThisTask,Etest); */


/*If the energy check is passed the accreted particles should be removed from the mesh*/

  if(Etest > 0)
    {
      for(i = 0; i < NumGas; i++)
        {
          if(P[i].ID == 0 && P[i].Mass == 0)
            continue;

          dx = NGB_PERIODIC_LONG_X(P[i].Pos[0] - SKD.Pos[0]);
          dy = NGB_PERIODIC_LONG_Y(P[i].Pos[1] - SKD.Pos[1]);
          dz = NGB_PERIODIC_LONG_Z(P[i].Pos[2] - SKD.Pos[2]);

          r2 = SKD.DistFac * (dx * dx + dy * dy + dz * dz);

          if(r2 < SKD.AccRad2)
            {
              printf("Removing particle i= %d ID= %d\n", i, P[i].ID);
              P[i].Mass = 0;
              P[i].ID = 0;

#ifdef VORONOI_DYNAMIC_UPDATE
              voronoi_remove_connection(i);
#endif

            }
        }
    }

  myfree_movable(SKD.ACC);

  return Etest;

}
#endif
