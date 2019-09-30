/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/GFM/winds.c
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

#include "../allvars.h"
#include "../proto.h"
#include "../voronoi.h"


#if defined(GFM_WINDS) || defined(GFM_WINDS_LOCAL)

static int recouple_wind_particles_evaluate(int target, int mode, int threadid);

static int Nwind;

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyIDType CellID;              /* cell to recouple with */
  MyDouble Pos[3];              /* wind particle position */
  MyFloat Vel[3];               /* wind particle velocity */
  MyDouble Mass;                /* wind particle mass */
  MyFloat Hsml;                 /* hsml for at least one ngb */
  MyFloat Utherm;
#ifdef GFM_STELLAR_EVOLUTION
  MyFloat MassMetals[GFM_N_CHEM_ELEMENTS];
  MyFloat Metallicity;
#endif
#ifdef GFM_DUST
  MyFloat InitialMetalFractions[GFM_N_CHEM_ELEMENTS];
  MyFloat InitialDustFractions[GFM_DUST_N_CHANNELS][GFM_N_CHEM_ELEMENTS];
#endif
#ifdef GFM_CHEMTAGS
  MyFloat MassMetalsChemTags[GFM_N_CHEM_TAGS];
#endif
#ifdef REFINEMENT_HIGH_RES_GAS
  MyFloat HighResMass;
#endif

  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;



 /* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in * in, int i, int firstnode)
{
  int k;

  in->CellID = WindParticle[i].CellID;
  for(k = 0; k < 3; k++)
    {
      in->Pos[k] = P[WindParticle[i].index].Pos[k];
      in->Vel[k] = P[WindParticle[i].index].Vel[k];
    }
  in->Mass = P[WindParticle[i].index].Mass;
  in->Hsml = STP(WindParticle[i].index).Hsml;
  in->Utherm = STP(WindParticle[i].index).Utherm;
#ifdef GFM_STELLAR_EVOLUTION
  for(k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
    in->MassMetals[k] = STP(WindParticle[i].index).MassMetals[k];
  in->Metallicity = STP(WindParticle[i].index).Metallicity;
#endif
#ifdef GFM_DUST
  for(k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
    {
      in->InitialMetalFractions[k] = STP(WindParticle[i].index).InitialMetalFractions[k];
      int chan;
      for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
        {
          in->InitialDustFractions[chan][k] = STP(WindParticle[i].index).InitialDustFractions[chan][k];
        }
    }
#endif
#ifdef GFM_CHEMTAGS
  for(k = 0; k < GFM_N_CHEM_TAGS; k++)
    in->MassMetalsChemTags[k] = STP(WindParticle[i].index).MassMetalsChemTags[k];
#endif
#ifdef REFINEMENT_HIGH_RES_GAS
  in->HighResMass = STP(WindParticle[i].index).HighResMass;
#endif

  in->Firstnode = firstnode;
}

 /* local data structure that holds results acquired on remote processors */
typedef struct
{
#ifdef TRACER_MC
  int attached_to_task;
  int attached_to_index;
#ifdef TRACER_MC_CHECKS
  MyIDType attached_to_ID;
#endif
#else
  char dummy;
#endif
} data_out;

static data_out *DataResult, *DataOut;



  /* routine to store or combine result data */
static void out2particle(data_out * out, int i, int mode)
{
#ifdef TRACER_MC
  if(out->attached_to_task >= 0)
    {
      if(mode == MODE_LOCAL_PARTICLES)  /* initial store */
        {
          if(out->attached_to_task != ThisTask)
            terminate("out->attached_to_task != ThisTask");

#ifdef TRACER_MC_CHECKS
          consider_moving_tracers(WindParticle[i].index, out->attached_to_task, out->attached_to_index, out->attached_to_ID, 1.0);
#else
          consider_moving_tracers(WindParticle[i].index, out->attached_to_task, out->attached_to_index, 0, 1.0);
#endif
        }
      else                      /* merge */
        {
#ifdef TRACER_MC_CHECKS
          consider_moving_tracers(WindParticle[i].index, out->attached_to_task, out->attached_to_index, out->attached_to_ID, 1.0);
#else
          consider_moving_tracers(WindParticle[i].index, out->attached_to_task, out->attached_to_index, 0, 1.0);
#endif
        }
    }
#endif
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
              if(generic_polling_primary(count, Nwind))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        i = NextParticle++;

        if(i >= Nwind)
          break;

        if(STP(WindParticle[i].index).BirthTime == 0)
          recouple_wind_particles_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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

        recouple_wind_particles_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}



/* counter for recoupled local wind particles */
static int local_recoupled;
/* total mass of recoupled local wind particles */
static double local_mass_recoupled;


void recouple_wind_particles(int nwind, int *ret_recoupled, double *ret_mass_recoupled)
{
  Nwind = nwind;

  long long ntot;
  sumup_large_ints(1, &Nwind, &ntot);
  if(ntot == 0)
    {
      mpi_printf("GFM_WINDS: Nothing to do.\n");
      *ret_recoupled = 0;
      *ret_mass_recoupled = 0;
      return;
    }

  generic_set_MaxNexport();


#ifdef TRACER_MC
  start_MC_tracer(N_tracer);    /* allocate buffer for tracer exchange */
#endif

  local_recoupled = 0;
  local_mass_recoupled = 0;

  generic_comm_pattern(Nwind, kernel_local, kernel_imported);


#ifdef TRACER_MC
  finish_MC_tracer();
#endif


  *ret_recoupled = local_recoupled;
  *ret_mass_recoupled = local_mass_recoupled;

  mpi_printf("GFM_WINDS: done with recoupling\n");
}


static int recouple_wind_particles_evaluate(int target, int mode, int threadid)
{
  int j, n, numnodes, *firstnode;
  MyDouble *pos, mass;
  MyFloat *vel, h, Utherm;
#ifdef REFINEMENT_HIGH_RES_GAS
  MyFloat highresmass;
#endif
#ifdef GFM_STELLAR_EVOLUTION
  MyFloat *massmetals, metallicity;
  int k;
#endif
#ifdef GFM_DUST
  int chan;
  MyFloat InitialMetalFractions[GFM_N_CHEM_ELEMENTS];
  MyFloat InitialDustFractions[GFM_DUST_N_CHANNELS][GFM_N_CHEM_ELEMENTS];
#endif
#ifdef GFM_CHEMTAGS
  MyFloat *massmetalschemtags;
#endif

  MyIDType cellID;
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

  cellID = in->CellID;
  pos = in->Pos;
  vel = in->Vel;
  mass = in->Mass;
  h = in->Hsml;
  Utherm = in->Utherm;
#ifdef GFM_STELLAR_EVOLUTION
  massmetals = in->MassMetals;
  metallicity = in->Metallicity;
#endif
#ifdef GFM_DUST
  memcpy(InitialMetalFractions, in->InitialMetalFractions, sizeof(MyFloat) * GFM_N_CHEM_ELEMENTS);
  memcpy(InitialDustFractions, in->InitialDustFractions, sizeof(MyFloat) * GFM_DUST_N_CHANNELS * GFM_N_CHEM_ELEMENTS);
#endif
#ifdef GFM_CHEMTAGS
  massmetalschemtags = in->MassMetalsChemTags;
#endif
#ifdef REFINEMENT_HIGH_RES_GAS
  highresmass = in->HighResMass;
#endif
#ifdef TRACER_MC
  out.attached_to_task = -1;
#endif

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, threadid, numnodes, firstnode);

  for(n = 0; n < nfound; n++)
    {
      j = Thread[threadid].Ngblist[n];

      /* nearest ngb cell? */
      if(P[j].ID == cellID)
        {
          local_recoupled++;
          local_mass_recoupled += mass;

          double inj_mass, inj_thermalenergy, inj_mom[3];
          inj_mass = mass;
          inj_mom[0] = mass * vel[0];
          inj_mom[1] = mass * vel[1];
          inj_mom[2] = mass * vel[2];
          inj_thermalenergy = All.cf_atime * All.cf_atime * mass * Utherm;
          gfm_inject_into_cell(j, inj_mass, inj_thermalenergy, inj_mom);

#ifdef GFM_STELLAR_EVOLUTION
#ifndef GFM_DUST
          for(k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
            SphP[j].MassMetals[k] += massmetals[k];
          SphP[j].MassMetallicity += mass * metallicity;
#else /* ifdef GFM_DUST */
          /* Use the relative makeup of the ISM when the wind was */
          /* launched to determine the recoupling into gas-phase */
          /* versus dust.  Since the wind's metallicity included */
          /* gas-phase and dust, subtract the dust contribution */
          /* when updating the SphP mass metallicity. */
          SphP[j].MassMetallicity += mass * metallicity;
          for(k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
            {
              SphP[j].MassMetals[k] += massmetals[k] * InitialMetalFractions[k];
              for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
                {
                  SphP[j].MassMetalsDust[chan][k] += massmetals[k] * InitialDustFractions[chan][k];
                  SphP[j].MassMetallicity -= massmetals[k] * InitialDustFractions[chan][k];
                }
            }
#endif
#endif
#ifdef GFM_CHEMTAGS
          for(k = 0; k < GFM_N_CHEM_TAGS; k++)
            SphP[j].MassMetalsChemTags[k] += massmetalschemtags[k];
#endif
#ifdef REFINEMENT_HIGH_RES_GAS
          SphP[j].HighResMass += highresmass;
#endif

#ifdef TRACER_MC
          out.attached_to_task = ThisTask;
          out.attached_to_index = j;
#ifdef TRACER_MC_CHECKS
          out.attached_to_ID = P[j].ID;
#endif
#else
          out.dummy = 0;
#endif
        }
    }

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}



#endif
