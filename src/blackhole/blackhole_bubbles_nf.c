/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/blackhole/blackhole_bubbles_nf.c
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


#ifdef BH_NF_RADIO

#define MAX_BUBBLES_PER_ACTIVE_HOLE 15
#define MIN_UFAC_FOR_BUBBLE_HEATING 0.1

static int blackhole_place_nf_bubbles_evaluate(int target, int mode, int threadid, int actioncode);

static struct localbubble_data
{
  MyDouble pos[3];
  double radius;
  double u;
  double umin;
  int pindex;
  double mass;
}
 *radio_bubble;


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
  double mass;
} data_out;

static data_out *DataResult, *DataOut;


 /* routine to store or combine result data */
static void out2particle(data_out * out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES)      /* initial store */
    {
      radio_bubble[i].mass = out->mass;
    }
  else                          /* combine */
    {
      radio_bubble[i].mass += out->mass;
    }
}


#include "../generic_comm_helpers2.h"




void blackhole_do_bubbles_nf(int num_activebh)
{
  int total_num_activebh, idx, i, n, l, num_bubbles, rep;

  MPI_Allreduce(&num_activebh, &total_num_activebh, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  mpi_printf("BH_NF_RADIO: The total number of active BHs is: %d\n", total_num_activebh);

  if(total_num_activebh > 0)
    {
      double num_placed = 0, total_placed;

      num_bubbles = MAX_BUBBLES_PER_ACTIVE_HOLE * num_activebh; /* this is the potential number of bubbles we place */

      radio_bubble = mymalloc("radio_bubble", num_bubbles * sizeof(struct localbubble_data));

      l = 0;
      for(idx = 0; idx < TimeBinsBHAccretion.NActiveParticles; idx++)
        {
          n = TimeBinsBHAccretion.ActiveParticleList[idx];
          if(n < 0)
            continue;

          if(BPP(n).BH_HaloVvir > 0)
            {
              double egy_thresh = blackhole_get_bubble_energy_thresh(n);

              if(BPP(n).BH_RadioEgyFeedback > egy_thresh && egy_thresh > 0)
                {
                  int idoff = 0;
                  double u = 3.0 / 4 * BPP(n).BH_HaloVvir * BPP(n).BH_HaloVvir;

                  double rvir_phys = BPP(n).BH_HaloVvir / (10 * All.cf_H);      /* physical */
                  double rvir = (rvir_phys / All.cf_atime);     /* comoving */

                  for(rep = 0; rep < MAX_BUBBLES_PER_ACTIVE_HOLE; rep++)
                    {
                      /* determine the distance of the bubble to the black hole */
                      double dist = get_random_number() * All.RadioRelativeMaxDist * rvir;

                      /* calculate bubble position */
                      double phi = 2 * M_PI * get_random_number();
                      double theta = acos(2 * get_random_number() - 1);

                      radio_bubble[l].pos[0] = P[n].Pos[0] + dist * sin(theta) * cos(phi);
                      radio_bubble[l].pos[1] = P[n].Pos[1] + dist * sin(theta) * sin(phi);
                      radio_bubble[l].pos[2] = P[n].Pos[2] + dist * cos(theta);
                      radio_bubble[l].radius = All.RadioRelativeBubbleSize * rvir;
                      radio_bubble[l].u = All.RadioRelativeBubbleEnergy * u;
                      radio_bubble[l].umin = MIN_UFAC_FOR_BUBBLE_HEATING * u;

                      radio_bubble[l].pindex = n;

                      l++;
                    }
                }
            }
        }

      if(l != num_bubbles)
        terminate("l != num_bubbles");

      /* first, measure just the mass in each of the bubbles */
      blackhole_place_nf_bubbles(num_bubbles, 0);

      for(i = 0; i < num_activebh; i++)
        {
          double mtot = 0;

          l = i * MAX_BUBBLES_PER_ACTIVE_HOLE;  /* first bubble index for this hole */

          for(rep = 0; rep < MAX_BUBBLES_PER_ACTIVE_HOLE; rep++)
            mtot += radio_bubble[l + rep].mass;

          if(mtot * radio_bubble[l].u > BPP(radio_bubble[l].pindex).BH_RadioEgyFeedback)
            {
              /* we do not need to use all bubbles */
              double egy_thresh = blackhole_get_bubble_energy_thresh(radio_bubble[l].pindex);

              for(rep = 0; rep < MAX_BUBBLES_PER_ACTIVE_HOLE; rep++)
                {
                  if(BPP(radio_bubble[l + rep].pindex).BH_RadioEgyFeedback > egy_thresh && radio_bubble[l + rep].mass > 0)
                    {
                      double egy = radio_bubble[l + rep].mass * radio_bubble[l + rep].u;

                      BPP(radio_bubble[l + rep].pindex).BH_RadioEgyFeedback -= egy;
                      BPP(radio_bubble[l + rep].pindex).BH_CumEgy_RM += egy;
                      AGNEnergyM_Should += egy;
                      num_placed++;
                    }
                  else
                    radio_bubble[l + rep].u = 0;
                }
            }
          else if(mtot > 0)
            {
              /* There is so much energy that even placing MAX_BUBBLES_PER_ACTIVE_HOLE bubbles with the nominal energy increase will not get rid off it all.
               * In this situation, we increase the amount of energy that is given per bubble.
               */

              double unew = BPP(radio_bubble[l].pindex).BH_RadioEgyFeedback / mtot;

              for(rep = 0; rep < MAX_BUBBLES_PER_ACTIVE_HOLE; rep++)
                {
                  radio_bubble[l + rep].u = unew;
                  double egy = radio_bubble[l + rep].mass * unew;

                  BPP(radio_bubble[l + rep].pindex).BH_RadioEgyFeedback -= egy;
                  BPP(radio_bubble[l + rep].pindex).BH_CumEgy_RM += egy;
                  AGNEnergyM_Should += egy;
                  num_placed++;
                }
            }
          else
            {
              /* we put no energy as there is no target mass in the fiducial bubble positions */

              for(rep = 0; rep < MAX_BUBBLES_PER_ACTIVE_HOLE; rep++)
                radio_bubble[l + rep].u = 0;
            }
        }

      /* now place the actual energy */
      blackhole_place_nf_bubbles(num_bubbles, 1);

      myfree(radio_bubble);

      MPI_Allreduce(&num_placed, &total_placed, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      mpi_printf("BH_NF_RADIO: The average number of bubbles placed per active BH: %g\n", total_placed / total_num_activebh);
    }
}

double blackhole_get_radio_efficiency(double vvir)
{
  double R = 0;

  if(vvir > 0)
    {
      double time_current = All.Time;

      if(All.ComovingIntegrationOn)
        {
          All.Time = 1.0;
          set_cosmo_factors_for_current_time();
          IonizeParams();
        }

      double Q = 2.5;

      double RefVelDisp = 180.0 * (1.0e5 / All.UnitVelocity_in_cm_per_s);
      double RefBHMass = 1.0e6 * SOLAR_MASS / (All.HubbleParam * All.UnitMass_in_g);

      double MB = RefBHMass * pow(vvir / RefVelDisp, 4);

      double u = 3.0 / 4 * pow(vvir, 2);

      /* choose a high fiducial physical density to avoid that we pick up a compton cooling contribution */
      double dens = 1.0e12 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

#ifdef GFM_COOLING_METAL
      double XHe_primordial = 1 - HYDROGEN_MASSFRAC;
      double XHe_solar = 4.0 * GFM_SOLAR_HE_ABUNDANCE / (1.0 + 4.0 * GFM_SOLAR_HE_ABUNDANCE + GFM_SOLAR_METALLICITY);
      double XHe = XHe_primordial + (XHe_solar - XHe_primordial) * All.RadioModeMetallicityInSolar;
      double XH = 1 - XHe - All.RadioModeMetallicityInSolar * GFM_SOLAR_METALLICITY;

      update_gas_state(dens, XH, All.RadioModeMetallicityInSolar * GFM_SOLAR_METALLICITY);
#if defined(GFM_AGN_RADIATION) || defined(GFM_UVB_CORRECTIONS)
      update_radiation_state(dens, XH, 0);
#endif
#endif

      double ne = 1.0;
      double tcool = GetCoolingTime(u, dens, &ne);

      double mdot = 2 * M_PI * Q * GAMMA_MINUS1 * tcool * dens * All.G * MB * pow(All.RadioModeMachnumber, 1.5);

      /* convert mdot to g/s */

      mdot *= All.UnitMass_in_g / All.UnitTime_in_s;

      double Lradio = All.BlackHoleRadiativeEfficiency * mdot * CLIGHT * CLIGHT;

      double Tvir = 35.9 * pow(vvir * All.UnitVelocity_in_cm_per_s / 1.0e5, 2); /* now in Kelvin */

      /*  Pratt et al. 2009, REXCESS */
      double Lx = 7.26e44 * pow(Tvir / (5.0 * 1.16e7), 2.62);   /* in ergs/sec */

      R = Lradio / Lx;

      if(All.ComovingIntegrationOn)
        {
          All.Time = time_current;
          set_cosmo_factors_for_current_time();
          IonizeParams();
        }
    }

  return R;
}

/* minimum bubble energy */
double blackhole_get_bubble_energy_thresh(int n)
{
  double rhocore = 1.0e4 * All.OmegaBaryon * 3 * All.cf_H * All.cf_H / (8 * M_PI * All.G);
  double u = 3.0 / 4 * BPP(n).BH_HaloVvir * BPP(n).BH_HaloVvir;
  double rvir = BPP(n).BH_HaloVvir / (10 * All.cf_H);

  double egythresh = All.RadioRelativeBubbleEnergy * u * rhocore * 4 * M_PI / 3.0 * pow(All.RadioRelativeBubbleSize * rvir, 3.0);

  return egythresh;
}


static int num_bubbles, actioncode;

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

        blackhole_place_nf_bubbles_evaluate(i, MODE_LOCAL_PARTICLES, threadid, actioncode);
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

        blackhole_place_nf_bubbles_evaluate(i, MODE_IMPORTED_PARTICLES, threadid, actioncode);
      }
  }
}



void blackhole_place_nf_bubbles(int num_bubbles_loc, int actioncode_loc)
{
  num_bubbles = num_bubbles_loc;
  actioncode = actioncode_loc;

  generic_set_MaxNexport();

  generic_comm_pattern(num_bubbles, kernel_local, kernel_imported);
}


static int blackhole_place_nf_bubbles_evaluate(int target, int mode, int threadid, int actioncode)
{
  int j, n;
  int numnodes, *firstnode;
  double h, h2, dx, dy, dz, r2, u, bubble_mass, umin, dens;

  double eos_dens_threshold = All.PhysDensThresh;
#ifdef MODIFIED_EOS
  eos_dens_threshold *= All.FactorDensThresh;
#endif

#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif
  MyDouble *pos;

  data_in local, *target_data;
  data_out out;

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

  pos = target_data->bubble.pos;
  h = target_data->bubble.radius;
  u = target_data->bubble.u;
  umin = target_data->bubble.umin;

  bubble_mass = 0;

  h2 = h * h;

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, threadid, numnodes, firstnode);

  for(n = 0; n < nfound; n++)
    {
      j = Thread[threadid].Ngblist[n];

      dx = NGB_PERIODIC_LONG_X(pos[0] - P[j].Pos[0]);
      dy = NGB_PERIODIC_LONG_Y(pos[1] - P[j].Pos[1]);
      dz = NGB_PERIODIC_LONG_Z(pos[2] - P[j].Pos[2]);

      r2 = dx * dx + dy * dy + dz * dz;

      dens = SphP[j].Density * All.cf_a3inv;

      if(r2 < h2 && P[j].Mass > 0 && P[j].ID != 0 && dens < eos_dens_threshold && SphP[j].Utherm > umin)
        {
          if(actioncode == 0)
            {
              bubble_mass += P[j].Mass;
            }
          else
            {
              double dE = u * P[j].Mass;
              SphP[j].Utherm += u;
              SphP[j].Energy += All.cf_atime * All.cf_atime * dE;
              AGNEnergyM_Is += dE;
            }
        }
    }

  out.mass = bubble_mass;

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}


#endif
