/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/blackhole/blackhole.c
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
#include <gsl/gsl_math.h>

#include "../allvars.h"
#include "../proto.h"

/*! \file blackhole_feedback.c
 *  \brief distribute feedback energy
 */


#ifdef BLACK_HOLES

#if defined(BH_THERMALFEEDBACK) || defined(BH_THERMALFEEDBACK_ACC)

static int blackhole_evaluate(int target, int mode, int threadid);

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble Pos[3];
  MyFloat Density;
  MyFloat BH_Hsml;
#ifdef BH_THERMALFEEDBACK
  MyFloat Energy;
#endif
#ifdef BH_THERMALFEEDBACK_ACC
  MyFloat BH_AccEnergy;
  MyFloat BH_AccTime;
#endif
#ifdef BH_BIPOLAR_FEEDBACK
  MyFloat BH_BipolarSum;
  MyFloat BH_BipolarJ[3];
#endif

  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

/* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in * in, int i, int firstnode)
{
  int k;
  for(k = 0; k < 3; k++)
    in->Pos[k] = P[i].Pos[k];

  in->BH_Hsml = BPP(i).BH_Hsml;
  in->Density = BPP(i).BH_Density * BPP(i).BH_VolSum;   /* to get back the plain sum over m_j W_j */

#ifdef BH_THERMALFEEDBACK
  in->Energy = BPP(i).BH_ThermEnergy;
#endif

#ifdef BH_THERMALFEEDBACK_ACC
  in->BH_AccEnergy = BPP(i).BH_AccEnergy;
  in->BH_AccTime = BPP(i).BH_AccTime;
#endif

#ifdef BH_BIPOLAR_FEEDBACK
  in->BH_BipolarSum = BPP(i).BH_BipolarSum;
  for(k = 0; k < 3; k++)
    {
      in->BH_BipolarJ[k] = BPP(i).BH_BipolarJ[k];
    }
#endif

  in->Firstnode = firstnode;
}




/* local data structure that holds results acquired on remote processors */
typedef struct
{
  char dummy;
} data_out;

static data_out *DataResult, *DataOut;



/* routine to store or combine result data */
static void out2particle(data_out * out, int i, int mode)
{

}


#include "../generic_comm_helpers2.h"


static void kernel_local(void)
{
  int i;
#ifdef GENERIC_ASYNC
  int flag = 0;
#endif

#pragma omp parallel private(i)
  {
    int j, threadid = get_thread_num();
#ifdef GENERIC_ASYNC
    int count = 0;
#endif

    for(j = 0; j < NTask; j++)
      Thread[threadid].Exportflag[j] = -1;

    while(1)
      {
        if(Thread[threadid].ExportSpace < MinSpace)
          break;

#ifdef GENERIC_ASYNC
        if(threadid == 0)
          {
            if((count & POLLINGINTERVAL) == 0)
              if(generic_polling_primary(count, TimeBinsBHAccretion.NActiveParticles))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        i = NextParticle++;

        if(i >= TimeBinsBHAccretion.NActiveParticles)
          break;

        int idx = TimeBinsBHAccretion.ActiveParticleList[i];

        if(idx < 0)
          continue;

        if(BPP(idx).BH_ThermEnergy == 0)
          continue;

        blackhole_evaluate(idx, MODE_LOCAL_PARTICLES, threadid);
      }
  }
}


static void kernel_imported(void)
{
  /* now do the particles that were sent to us */
  int i, cnt = 0;
#pragma omp parallel private(i)
  {
    int threadid = get_thread_num();
#ifdef GENERIC_ASYNC
    int count = 0;
#endif

    while(1)
      {
#pragma omp atomic capture
        i = cnt++;

        if(i >= Nimport)
          break;

#ifdef GENERIC_ASYNC
        if(threadid == 0)
          {
            if((count & POLLINGINTERVAL) == 0)
              generic_polling_secondary();
          }

        count++;
#endif

        blackhole_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }

}

/* This function distributes the feedback energy and look for gas cells that we might stochastically swallow
 */
void blackhole_assign_feedback(void)
{
  mpi_printf("BH_THERMALFEEDBACK: Start assigning BH feedback\n");
  double t0 = second();

  generic_set_MaxNexport();

  /** Let's first spread the feedback energy, and determine which particles may be swallowed and by whom */

  generic_comm_pattern(TimeBinsBHAccretion.NActiveParticles, kernel_local, kernel_imported);


  for(int i = 0; i < TimeBinsBHAccretion.NActiveParticles; i++)
    {
      int idx = TimeBinsBHAccretion.ActiveParticleList[i];
      if(idx < 0)
        continue;

      BPP(idx).BH_ThermEnergy = 0;
    }

  double t1 = second();
  mpi_printf("BH_THERMALFEEDBACK: Done assigning BH thermal feedback (%g sec)\n", timediff(t0, t1));
}




static int blackhole_evaluate(int target, int mode, int threadid)
{
  int numnodes, *firstnode, j, n;
  MyDouble *pos;
  MyFloat rho;
  double dx, dy, dz, h_i, h_i2, r2, r, u, hinv, hinv3, wk;
#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif

#ifdef BH_THERMALFEEDBACK_ACC
  double accenergy = 0;
  double energy_to_temp_fac =
    (4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC))) * PROTONMASS / BOLTZMANN * GAMMA_MINUS1 * All.UnitEnergy_in_cgs / All.UnitMass_in_g / (All.ReferenceGasPartMass * All.DesNumNgbBlackHole);
  double delta_temp;
  double delta_time = 0;
#endif

  data_in local, *target_data;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      target_data = &local;

      numnodes = 1;
      firstnode = NULL;
    }
  else
    {
      target_data = &DataGet[target];

      generic_get_numnodes(target, &numnodes, &firstnode);
    }

  pos = target_data->Pos;
  rho = target_data->Density;

  h_i = target_data->BH_Hsml;

#ifdef BH_THERMALFEEDBACK
  double energy = target_data->Energy;
#endif

#ifdef BH_THERMALFEEDBACK_ACC
  accenergy = target_data->BH_AccEnergy;
  delta_time = target_data->BH_AccTime;
  delta_temp = energy_to_temp_fac * accenergy;
#endif

#ifdef BH_BIPOLAR_FEEDBACK
  rho = target_data->BH_BipolarSum;
  MyDouble *bipolar_j = target_data->BH_BipolarJ;
#endif

  h_i2 = h_i * h_i;


  int nfound = ngb_treefind_variable_threads(pos, h_i, target, mode, threadid, numnodes, firstnode);

  for(n = 0; n < nfound; n++)
    {
      j = Thread[threadid].Ngblist[n];

      if(P[j].Mass > 0)
        {
          dx = NGB_PERIODIC_LONG_X(pos[0] - P[j].Pos[0]);
          dy = NGB_PERIODIC_LONG_Y(pos[1] - P[j].Pos[1]);
          dz = NGB_PERIODIC_LONG_Z(pos[2] - P[j].Pos[2]);

          r2 = dx * dx + dy * dy + dz * dz;

          if(r2 < h_i2 && P[j].Mass > 0 && P[j].ID != 0)
            {
              if(P[j].Type == 0)
                {

#ifdef BH_BIPOLAR_FEEDBACK
                  if(!is_cell_in_bipolar_cone(j,pos,bipolar_j))
                    {
                      /* we only inject feedback into the cone */
                      continue;
                    }
#endif
                  r = sqrt(r2);
                  hinv = 1 / h_i;
                  hinv3 = hinv * hinv * hinv;

                  u = r * hinv;

                  if(u < 0.5)
                    wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
                  else
                    wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

                  if(rho > 0.0)
                    {

#ifdef BH_THERMALFEEDBACK_ACC
                      if(blackhole_check_duty_cycle(delta_temp, delta_time))
                        SphP[j].Injected_BH_Energy += accenergy * P[j].Mass * wk / rho;
#endif
#ifdef BH_THERMALFEEDBACK
                      SphP[j].Injected_BH_Energy += energy * P[j].Mass * wk / rho;
#endif
                    }
                }
            }
        }
    }

  return 0;
}

#endif



#ifdef BH_THERMALFEEDBACK_ACC
void blackhole_accumulate_energy(void)
{
  int idx, i;
  double meddington, energy, mdot, bh_mass, dt;

  for(idx = 0; idx < TimeBinsBHAccretion.NActiveParticles; idx++)
    {
      i = TimeBinsBHAccretion.ActiveParticleList[idx];
      if(i < 0)
        continue;

      mdot = BPP(i).BH_Mdot;
      bh_mass = BPP(i).BH_Mass;
      dt = (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval / All.cf_hubble_a;

      meddington = blackhole_mdot_eddington(bh_mass);

      energy = All.BlackHoleFeedbackFactor * All.BlackHoleRadiativeEfficiency * mdot * dt * pow(CLIGHT / All.UnitVelocity_in_cm_per_s, 2);

      if(blackhole_get_mode(mdot, meddington) == BH_RADIO_MODE)
        energy = 0.0;

      BPP(i).BH_AccEnergy += energy;

      BPP(i).BH_AccTime += dt / All.HubbleParam * All.UnitTime_in_s / SEC_PER_GIGAYEAR; /* increase time in Gyrs */
    }
}


void blackhole_reset_energy(void)
{
  int idx i;

  double delta_temp, delta_time;
  double energy_to_temp_fac =
    (4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC))) * PROTONMASS / BOLTZMANN * GAMMA_MINUS1 * All.UnitEnergy_in_cgs / All.UnitMass_in_g / (All.ReferenceGasPartMass * All.DesNumNgbBlackHole);

  for(idx = 0; idx < TimeBinsBHAccretion.NActiveParticles; idx++)
    {
      i = TimeBinsBHAccretion.ActiveParticleList[idx];
      if(i < 0)
        continue;

      delta_temp = energy_to_temp_fac * BPP(i).BH_AccEnergy;
      delta_time = BPP(i).BH_AccTime;
      if(blackhole_check_duty_cycle(delta_temp, delta_time))
        {
          BPP(i).BH_AccEnergy = 0.0;
          BPP(i).BH_AccTime = 0.0;
        }
    }
}


int blackhole_check_duty_cycle(double delta_temp, double delta_time)
{
  int flag1, flag2;

  if(delta_time >= All.BlackholeDeltaTime)
    flag1 = 1;
  else
    flag1 = 0;

  if(delta_temp >= All.BlackholeDeltaTemp)
    flag2 = 1;
  else
    flag2 = 0;

  if(flag1 * flag2 != 0)
    return 1;
  else
    return 0;
}

#endif /* end BH_THERMALFEEDBACK_ACC */


#ifdef BH_DRAG
/** This function is a place holder for some sub-grid modelling of dynamical friction
 *  on the black hole.
 *  However, the current implementation just accounts for momentum transfer from the
 *  accreted gas (in the sub-grid sense, via the instantaneous BH_Mdot, which not necessarily
 *  the actual mass accreted) to the black hole. This momentum transfer is already implemented
 *  in blackhole_swallowgas.c whenever gas mass is actually transferred to the black hole.
 *  Therefore, in its current form, BH_DRAG should NOT be used.
 */
void blackhole_drag_force(void)
{
  int idx, n;
  double fac, dt;

  for(idx = 0; idx < TimeBinsBHAccretion.NActiveParticles; idx++)
    {
      n = TimeBinsBHAccretion.ActiveParticleList[idx];
      if(n < 0)
        continue;

      dt = (P[n].TimeBinHydro ? (((integertime) 1) << P[n].TimeBinHydro) : 0) * All.Timebase_interval / All.cf_hubble_a;
      if(BPP(n).BH_Mass > 0)
        {
          fac = BPP(n).BH_Mdot * dt / BPP(n).BH_Mass;

          if(fac > 1)
            fac = 1;

          /* TODO: unclear what to do here when we use DECOUPLE_TIMESTEPS, but this part of the code seems to be broken anyway? */
          if(dt > 0)
            {
              int k;
              for(k = 0; k < 3; k++)
                P[n].GravAccel[k] += -All.cf_atime * All.cf_atime * fac / dt * (P[n].Vel[k] - BPP(n).BH_SurroundingGasVel[k]) / All.cf_atime;
            }
        }
    }
}
#endif


#endif
