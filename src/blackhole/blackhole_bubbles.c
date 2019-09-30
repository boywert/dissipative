/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/blackhole/blackhole_bubbles.c
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

/*! \file blackhole_bubbles.c
 *  \brief routines for placing bubbles in radio mode
 */

#ifdef BH_BUBBLES
#ifdef BH_MAGNETIC_BUBBLES
#ifndef MHD
#error BH_MAGNETIC_BUBBLES requires MHD
#endif

#ifdef BH_MAGNETIC_DIPOLAR_BUBBLES
static double set_dipole_strength(double bubbleVol, double soft, double magneticEnergy);
#endif
#endif

#define MEASURE_BUBBLE_CONTENT  0
#define PLACE_BUBBLE_FEEDBACK   1

#define TOLERANCE_FACTOR        1.001

static void bh_assign_default_bubble_size_and_position(void);
static void bh_work_on_bubbles(int actioncode);
static int blackhole_place_bubble_evaluate(int target, int mode, int threadid, int actioncode);
static void bh_rescale_bubble_size_and_position(void);
static void bh_final_bubble_preparations(void);


static double local_energy_bubbles, local_mass_bubbles;

static int action;
static double u_to_temp_fac;

static int num_bubbles;

static struct localbubble_data
{
  MyDouble pos[3];
  MyDouble bh_pos[3];
  double dir[3];
  double bh_dmass;
  MyIDType bh_id;

  double radius;                /* bubble radius */
  double original_radius;
  double min_radius;
  double bubbleEnergy;
  double distance;

  double E_bubble;
  double Mass_bubble;
  double Volume_bubble;
  int numngb;

#ifdef BH_MAGNETIC_BUBBLES
#ifdef BH_MAGNETIC_DIPOLAR_BUBBLES
  double m[3];
  double mstrength;
#endif
#endif
} *radio_bubble;

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  struct localbubble_data bubble;

  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

/* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in * in, int i, int firstnode)
{
  in->bubble = radio_bubble[i];

  in->Firstnode = firstnode;
}

/* local data structure that holds results acquired on remote processors */
typedef struct
{
  double E_bubble;
  double Mass_bubble;
  double Volume_bubble;
  double min_radius;
  int numngb;
} data_out;

static data_out *DataResult, *DataOut;

/* routine to store or combine result data */
static void out2particle(data_out * out, int i, int mode)
{
  if(action == MEASURE_BUBBLE_CONTENT)
    {
      if(mode == MODE_LOCAL_PARTICLES)  /* initial store */
        {
          radio_bubble[i].min_radius = out->min_radius;
          radio_bubble[i].E_bubble = out->E_bubble;
          radio_bubble[i].Mass_bubble = out->Mass_bubble;
          radio_bubble[i].Volume_bubble = out->Volume_bubble;
          radio_bubble[i].numngb = out->numngb;
        }
      else                      /* combine */
        {
          if(radio_bubble[i].min_radius > out->min_radius)
            radio_bubble[i].min_radius = out->min_radius;

          radio_bubble[i].E_bubble += out->E_bubble;
          radio_bubble[i].Mass_bubble += out->Mass_bubble;
          radio_bubble[i].Volume_bubble += out->Volume_bubble;
          radio_bubble[i].numngb += out->numngb;
        }
    }
}

#include "../generic_comm_helpers2.h"

void blackhole_do_bubbles(int num_activebh)
{
  int idx;
  double tot_energy_bubbles, tot_mass_bubbles;
  local_energy_bubbles = local_mass_bubbles = 0.0;

  int total_num_activebh = 0;

  /* note: assuming FULL ionization */
  u_to_temp_fac = (4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC))) * PROTONMASS / BOLTZMANN * GAMMA_MINUS1 * All.UnitEnergy_in_cgs / All.UnitMass_in_g;

  MPI_Allreduce(&num_activebh, &total_num_activebh, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  mpi_printf("BH_BUBBLES: The total number of active BHs is: %d\n", total_num_activebh);

  if(total_num_activebh > 0)
    {
      radio_bubble = mymalloc("radio_bubble", num_activebh * sizeof(struct localbubble_data));
      num_bubbles = 0;

      for(idx = 0; idx < TimeBinsBHAccretion.NActiveParticles; idx++)
        {
          int n = TimeBinsBHAccretion.ActiveParticleList[idx];
          if(n < 0)
            continue;

          if(BPP(n).BH_Mass_bubbles > 0 && BPP(n).BH_Mass_bubbles > All.BlackHoleRadioTriggeringFactor * BPP(n).BH_Mass_ini)
            {
              radio_bubble[num_bubbles].bh_dmass = BPP(n).BH_Mass_bubbles - BPP(n).BH_Mass_ini;

              BPP(n).BH_CumEgy_RM += All.RadioFeedbackFactor * All.BlackHoleRadiativeEfficiency * radio_bubble[num_bubbles].bh_dmass * pow(CLIGHT / All.UnitVelocity_in_cm_per_s, 2);

              BPP(n).BH_Mass_ini = BPP(n).BH_Mass;
              BPP(n).BH_Mass_bubbles = BPP(n).BH_Mass;

              radio_bubble[num_bubbles].bh_pos[0] = P[n].Pos[0];
              radio_bubble[num_bubbles].bh_pos[1] = P[n].Pos[1];
              radio_bubble[num_bubbles].bh_pos[2] = P[n].Pos[2];
              radio_bubble[num_bubbles].bh_id = P[n].ID;

              if(radio_bubble[num_bubbles].bh_dmass > 0)
                num_bubbles++;
            }
        }

      bh_assign_default_bubble_size_and_position();

      bh_work_on_bubbles(MEASURE_BUBBLE_CONTENT);

      bh_rescale_bubble_size_and_position();

      bh_work_on_bubbles(MEASURE_BUBBLE_CONTENT);

      bh_final_bubble_preparations();

      bh_work_on_bubbles(PLACE_BUBBLE_FEEDBACK);

      myfree(radio_bubble);
    }


  MPI_Reduce(&local_energy_bubbles, &tot_energy_bubbles, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_mass_bubbles, &tot_mass_bubbles, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  mpi_printf("BH_BUBBLES: injected a total energy of %g into a total bubble mass of %g\n", tot_energy_bubbles, tot_mass_bubbles);
}

static void bh_rescale_bubble_size_and_position(void)
{
  int i, j;

  for(i = 0; i < num_bubbles; i++)
    {
      /* calculate the ICM density inside the bubble */
      double ICMDensity = radio_bubble[i].Mass_bubble / (4.0 * M_PI / 3.0 * pow(radio_bubble[i].radius, 3));
      ICMDensity = ICMDensity / (pow(All.cf_atime, 3)); /*now physical */

      /* calculate parameters of rescaled bubble */
      radio_bubble[i].min_radius = MAX_DOUBLE_NUMBER;
      radio_bubble[i].radius *= pow((radio_bubble[i].bubbleEnergy * All.DefaultICMDensity / (All.BubbleEnergy * ICMDensity)), 1. / 5.);
      radio_bubble[i].distance *= pow((radio_bubble[i].bubbleEnergy * All.DefaultICMDensity / (All.BubbleEnergy * ICMDensity)), 1. / 5.);
      radio_bubble[i].original_radius = radio_bubble[i].radius;

      /*switch to comoving if it is assumed that Rbub should be constant with redshift */
      /* BubbleRadius = BubbleRadius / All.cf_atime;
         BubbleDistance = BubbleDistance / All.cf_atime; */

      for(j = 0; j < 3; j++)
        radio_bubble[i].pos[j] = radio_bubble[i].distance * radio_bubble[i].dir[j] + radio_bubble[i].bh_pos[j];

      radio_bubble[i].numngb = 0;
    }
}

static void bh_assign_default_bubble_size_and_position(void)
{
  int i, j;
  for(i = 0; i < num_bubbles; i++)
    {
      radio_bubble[i].min_radius = MAX_DOUBLE_NUMBER;
      radio_bubble[i].radius = All.BubbleRadius;
      radio_bubble[i].original_radius = radio_bubble[i].radius;
      radio_bubble[i].distance = pow(get_random_number(), 1. / 3.) * All.BubbleDistance;

      /*switch to comoving if it is assumed that Rbub should be constant with redshift */
      /* BubbleDistance = All.BubbleDistance / All.cf_atime;
         BubbleRadius = All.BubbleRadius / All.cf_atime; */

      radio_bubble[i].bubbleEnergy = All.RadioFeedbackFactor * All.BlackHoleRadiativeEfficiency * radio_bubble[i].bh_dmass * All.UnitMass_in_g / All.HubbleParam * pow(CLIGHT, 2);      /*in cgs units */

      AGNEnergyM_Should += radio_bubble[i].bubbleEnergy * All.HubbleParam / All.UnitEnergy_in_cgs;      /* in internal units */


      double phi = 2 * M_PI * get_random_number();
      double theta = acos(2 * get_random_number() - 1);

      radio_bubble[i].dir[0] = sin(theta) * cos(phi);
      radio_bubble[i].dir[1] = sin(theta) * sin(phi);
      radio_bubble[i].dir[2] = cos(theta);

      for(j = 0; j < 3; j++)
        radio_bubble[i].pos[j] = radio_bubble[i].distance * radio_bubble[i].dir[j] + radio_bubble[i].bh_pos[j];

      radio_bubble[i].numngb = 0;
    }
}

static void bh_final_bubble_preparations(void)
{
#ifdef BH_MAGNETIC_BUBBLES
#ifdef BH_MAGNETIC_DIPOLAR_BUBBLES

  int i;
  for(i = 0; i < num_bubbles; i++)
    {
      /* energy we want to inject in total */
      double dE = (radio_bubble[i].bubbleEnergy * All.HubbleParam / All.UnitEnergy_in_cgs);
      /* note: if this applies, a difference between AGNEnergyM_Should and AGNEnergyM_Is will appear */
      if(u_to_temp_fac * dE / radio_bubble[i].Mass_bubble > 5.0e9)
        dE = 5.0e9 * radio_bubble[i].Mass_bubble / u_to_temp_fac;

      double soft = get_default_softening_of_particletype(0) * All.cf_atime;    /* current physical gas softening length */
      radio_bubble[i].mstrength = set_dipole_strength(radio_bubble[i].Volume_bubble / All.cf_a3inv, soft, All.MagneticEnergyFraction * dE);     /* in physical units */

      double phi = 2 * M_PI * get_random_number();
      double theta = acos(2 * get_random_number() - 1);

      /* transform B field in comoving units (CHECK ME!!!) */
      radio_bubble[i].mstrength *= All.cf_atime * All.cf_atime;
      soft /= All.cf_atime;

      radio_bubble[i].m[0] = radio_bubble[i].mstrength * sin(theta) * cos(phi);
      radio_bubble[i].m[1] = radio_bubble[i].mstrength * sin(theta) * sin(phi);
      radio_bubble[i].m[2] = radio_bubble[i].mstrength * cos(theta);
    }

#endif
#endif
}



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
              if(generic_polling_primary(count, num_bubbles))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        i = NextParticle++;

        if(i >= num_bubbles)
          break;

        /* Don't reprocess a bubble whose content has been measured */
        if((action == MEASURE_BUBBLE_CONTENT) && (radio_bubble[i].numngb > 0))
          continue;

        blackhole_place_bubble_evaluate(i, MODE_LOCAL_PARTICLES, threadid, action);
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

        blackhole_place_bubble_evaluate(i, MODE_IMPORTED_PARTICLES, threadid, action);
      }
  }
}

static void bh_work_on_bubbles(int actioncode)
{
  int i;
  int nleft, nleft_all, iter = 0;

  action = actioncode;

  generic_set_MaxNexport();

  do
    {
      generic_comm_pattern(num_bubbles, kernel_local, kernel_imported);

      if(actioncode == MEASURE_BUBBLE_CONTENT)
        {
          nleft = 0;

          for(i = 0; i < num_bubbles; i++)
            {
              if(radio_bubble[i].numngb == 0)
                {
                  radio_bubble[i].radius *= 1.26;
                  nleft++;
                }
              else if(radio_bubble[i].radius > radio_bubble[i].original_radius)
                {
                  radio_bubble[i].radius = TOLERANCE_FACTOR * radio_bubble[i].min_radius;
                  radio_bubble[i].original_radius = radio_bubble[i].radius;

                  /* force processing of the bubble to recompute its properties */
                  radio_bubble[i].numngb = 0;
                  nleft++;
                }
            }

          MPI_Allreduce(&nleft, &nleft_all, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        }
      else
        nleft_all = 0;

      iter++;
      if(iter > MAXITER)
        terminate("BH_BUBBLES: failed to converge in neighbour iteration in bh_bubble()");
    }
  while(nleft_all > 0);
}

static int blackhole_place_bubble_evaluate(int target, int mode, int threadid, int actioncode)
{
  int j, n, numnodes, *firstnode;

  data_in local, *in;
  data_out out;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      in = &local;

      numnodes = 1;
      firstnode = NULL;
    }
  else
    {
      in = &DataGet[target];

      generic_get_numnodes(target, &numnodes, &firstnode);
    }

  MyDouble *pos = in->bubble.pos;
  double h = in->bubble.radius;
  double totMass_bubble = in->bubble.Mass_bubble;
  double bubbleEnergy = in->bubble.bubbleEnergy;



  out.E_bubble = 0;
  out.Mass_bubble = 0;
  out.Volume_bubble = 0;
  out.numngb = 0;
  out.min_radius = MAX_DOUBLE_NUMBER;


  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, threadid, numnodes, firstnode);

  for(n = 0; n < nfound; n++)
    {
      j = Thread[threadid].Ngblist[n];

      if(P[j].Type == 0 && P[j].Mass > 0)
        {
          if(actioncode == MEASURE_BUBBLE_CONTENT)
            {
              out.numngb++;
              out.E_bubble += SphP[j].Utherm * P[j].Mass;
              out.Mass_bubble += P[j].Mass;
              out.Volume_bubble += SphP[j].Volume;
#ifdef PERIODIC
              double xtmp, ytmp, ztmp;
#endif
              double dx = NGB_PERIODIC_LONG_X(pos[0] - P[j].Pos[0]);
              double dy = NGB_PERIODIC_LONG_Y(pos[1] - P[j].Pos[1]);
              double dz = NGB_PERIODIC_LONG_Z(pos[2] - P[j].Pos[2]);
              double r2 = dx * dx + dy * dy + dz * dz;
              double r = sqrt(r2);
              if(r < out.min_radius)
                out.min_radius = r;
            }
          else if(actioncode == PLACE_BUBBLE_FEEDBACK)
            {
#ifndef BH_MAGNETIC_BUBBLES
              /* energy we want to inject in this particle */
              double dE = ((bubbleEnergy * All.HubbleParam / All.UnitEnergy_in_cgs) / totMass_bubble) * P[j].Mass;

              /* note: if this applies, a difference between AGNEnergyM_Should and AGNEnergyM_Is will appear */
              if(u_to_temp_fac * dE / P[j].Mass > 5.0e9)
                dE = 5.0e9 * P[j].Mass / u_to_temp_fac;

              SphP[j].Utherm += dE / P[j].Mass;
              SphP[j].Energy += All.cf_atime * All.cf_atime * dE;
#else
#ifndef BH_MAGNETIC_DIPOLAR_BUBBLES
              double totVolume_bubble = in->bubble.Volume_bubble;

              /* energy we want to inject in total */
              double dE = (bubbleEnergy * All.HubbleParam / All.UnitEnergy_in_cgs);

              /* note: if this applies, a difference between AGNEnergyM_Should and AGNEnergyM_Is will appear */
              if(u_to_temp_fac * dE / totMass_bubble > 5.0e9)
                dE = 5.0e9 * totMass_bubble / u_to_temp_fac;

              /* thermal energy we want to inject in this particle */
              double dEtherm = (1.0 - All.MagneticEnergyFraction) * dE * P[j].Mass / totMass_bubble;
              SphP[j].Utherm += dEtherm;

              /* magnetic energy we want to inject in this particle */
              double dEmagnetic = All.MagneticEnergyFraction * dE * SphP[j].Volume / totVolume_bubble;

              double dBConserved = sqrt(2.0 * dEmagnetic * All.cf_atime / SphP[j].Volume) * SphP[j].Volume;
              double B = sqrt(SphP[j].BConserved[0] * SphP[j].BConserved[0] + SphP[j].BConserved[1] * SphP[j].BConserved[1] + SphP[j].BConserved[2] * SphP[j].BConserved[2]);
              double dir[3];
              dir[0] = SphP[j].BConserved[0] / B;
              dir[1] = SphP[j].BConserved[1] / B;
              dir[2] = SphP[j].BConserved[2] / B;
              SphP[j].BConserved[0] += dBConserved * dir[0];
              SphP[j].BConserved[1] += dBConserved * dir[1];
              SphP[j].BConserved[2] += dBConserved * dir[2];

              /* reset value of dE so that statistics come out right (for the next particle dE will be recomputed) */
              dE = dEtherm + dEmagnetic;
              /* update particle total energy */
              SphP[j].Energy += All.cf_atime * All.cf_atime * dE;
#else
              /* energy we want to inject in this particle */
              double dE = ((bubbleEnergy * All.HubbleParam / All.UnitEnergy_in_cgs) / totMass_bubble) * P[j].Mass;
              double *m = in->bubble.m;
              double mstrength = in->bubble.mstrength;

              /* note: if this applies, a difference between AGNEnergyM_Should and AGNEnergyM_Is will appear     */
              if(u_to_temp_fac * dE / P[j].Mass > 5.0e9)
                dE = 5.0e9 * P[j].Mass / u_to_temp_fac;

              SphP[j].Utherm += (1.0 - All.MagneticEnergyFraction) * dE / P[j].Mass;
              SphP[j].Energy += All.cf_atime * All.cf_atime * dE;

#ifdef PERIODIC
              double xtmp, ytmp, ztmp;
#endif
              double dx = NGB_PERIODIC_LONG_X(pos[0] - P[j].Pos[0]);
              double dy = NGB_PERIODIC_LONG_Y(pos[1] - P[j].Pos[1]);
              double dz = NGB_PERIODIC_LONG_Z(pos[2] - P[j].Pos[2]);
              double r2 = dx * dx + dy * dy + dz * dz;

              double soft = get_default_softening_of_particletype(0) * All.cf_atime;

              double r = sqrt(r2);
              if(r > 0.0)
                {
                  double dir[3];
                  dir[0] = -dx / r;
                  dir[1] = -dy / r;
                  dir[2] = -dz / r;

                  double mdotn = m[0] * dir[0] + m[1] * dir[1] + m[2] * dir[2];
                  SphP[j].BConserved[0] += (3.0 * mdotn * dir[0] - m[0]) / (r * r2 + soft * soft * soft) * SphP[j].Volume;
                  SphP[j].BConserved[1] += (3.0 * mdotn * dir[1] - m[1]) / (r * r2 + soft * soft * soft) * SphP[j].Volume;
                  SphP[j].BConserved[2] += (3.0 * mdotn * dir[2] - m[2]) / (r * r2 + soft * soft * soft) * SphP[j].Volume;
                }
              else
                {
                  double dir[3];
                  /* B field is aligned with the magnetic moment vector at the very center of the dipole */
                  dir[0] = m[0] / mstrength;
                  dir[1] = m[1] / mstrength;
                  dir[2] = m[2] / mstrength;

                  double mdotn = m[0] * dir[0] + m[1] * dir[1] + m[2] * dir[2];
                  SphP[j].BConserved[0] += (3.0 * mdotn * dir[0] - m[0]) / (soft * soft * soft) * SphP[j].Volume;
                  SphP[j].BConserved[1] += (3.0 * mdotn * dir[1] - m[1]) / (soft * soft * soft) * SphP[j].Volume;
                  SphP[j].BConserved[2] += (3.0 * mdotn * dir[2] - m[2]) / (soft * soft * soft) * SphP[j].Volume;
                }
#endif
#endif
              AGNEnergyM_Is += dE;
              local_energy_bubbles += dE;
              local_mass_bubbles += P[j].Mass;
            }
        }
    }


  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#ifdef BH_MAGNETIC_BUBBLES
#ifdef BH_MAGNETIC_DIPOLAR_BUBBLES
double set_dipole_strength(double bubbleVol, double soft, double magneticEnergy)
{
  double bubbleRadiusCubed = 3.0 / (4.0 * M_PI) * bubbleVol;
  double softCubed = soft * soft * soft;
  double m = sqrt(3. * magneticEnergy * softCubed * (bubbleRadiusCubed + softCubed) / bubbleRadiusCubed);

  /* such that the field is in Lorentz-Heaviside units already */
  return m / sqrt(4.0 * M_PI);
}
#endif
#endif

#endif /* end of BH_BUBBLE */
