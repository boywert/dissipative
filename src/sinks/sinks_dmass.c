/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/sinks/sinks_accrete.c
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

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble Pos[3];
  signed char SinkStep;

  int Firstnode;
} data_in;

static data_in *DataGet, *DataIn;

/* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in *in, int i, int firstnode)
{
  for(int k = 0; k < 3; k++)
    in->Pos[k] = P[SinksAux[i].SinksAuxID].Pos[k];

  in->SinkStep = P[SinksAux[i].SinksAuxID].TimeBinGrav;

  in->Firstnode = firstnode;
}

/* local data structure that holds results acquired on remote processors */
typedef struct
{
  MyFloat AccretedMass;
  MyFloat TotAccretedMass;
  signed char MinTimeBin;
} data_out;

static data_out *DataResult, *DataOut;

 /* routine to store or combine result data */
static void out2particle(data_out *out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES)      /* initial store */
    {
      SinksAux[i].dMass = out->AccretedMass; 
      SinksAux[i].MassNorm = out->TotAccretedMass; 
      SinksAux[i].MinTimeBin = out->MinTimeBin;
    }
  else                          /* combine */
    {
      SinksAux[i].dMass += out->AccretedMass;
      SinksAux[i].MassNorm += out->TotAccretedMass; 

      if(out->MinTimeBin < SinksAux[i].MinTimeBin)
        SinksAux[i].MinTimeBin = out->MinTimeBin;
    }
}

static int sinks_dmass_evaluate(int target, int mode, int threadid);

#include "../generic_comm_helpers2.h"

static void kernel_local(void)
{
  int i;
#pragma omp parallel private(i)
    {
      int j, threadid = get_thread_num();

      for(j = 0; j < NTask; j++)
        Thread[threadid].Exportflag[j] = -1;

      while(1)
        {
          if(Thread[threadid].ExportSpace < MinSpace)
            break;

#pragma omp atomic capture
          i = NextParticle++;

          if(i >= NSinks)
            break;

          sinks_dmass_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
        }
    }
}

static void kernel_imported(void)
{
  /* now do the particles that were sent to us */
  int i, count = 0;
#pragma omp parallel private(i)
  {
    int threadid = get_thread_num();

    while(1)
      {
#pragma omp atomic capture
        i = count++;

        if(i >= Nimport)
          break;

        sinks_dmass_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

void sinks_dmass(void)
{
  generic_set_MaxNexport();

  generic_comm_pattern(NSinks, kernel_local, kernel_imported);

  mpi_printf("SINKS: computed mass accretion rate.\n");
}

static int sinks_dmass_evaluate(int target, int mode, int threadid)
{
  int j, n, numnodes, *firstnode;
  double rad;
  signed char min_time_bin = TIMEBINS;
  signed char sinkstep;
  MyFloat accreted_mass, tot_accreted_mass;
  MyDouble *pos;

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

  pos = in->Pos;
  sinkstep = in->SinkStep;

  rad = SKD.AccRad / sqrt(SKD.DistFac);

  accreted_mass = 0;
  tot_accreted_mass = 0;

  int nfound = ngb_treefind_variable_threads(pos, rad, target, mode, threadid, numnodes, firstnode);

  for(n = 0; n < nfound; n++)
    {
      j = Thread[threadid].Ngblist[n];

      if(P[j].Type == 0)
        {
          if(P[j].Mass > 0 && P[j].ID != 0)
            {
              signed char binstep = P[j].TimeBinHydro;
              
              if(binstep < sinkstep)
                binstep = sinkstep;

              integertime dti_step = (binstep ? (((integertime) 1) << binstep) : 0);

              if(dti_step > 0)
                accreted_mass += P[j].Mass / dti_step;

              if(P[j].TimeBinHydro < min_time_bin)
                min_time_bin = P[j].TimeBinHydro;

              /* mark the cell as eligible for derefinement */
              SphP[j].InAccrRadius = 1;

              tot_accreted_mass += P[j].Mass;
            }
        }
    }

  out.AccretedMass = accreted_mass;
  out.TotAccretedMass = tot_accreted_mass;
  out.MinTimeBin = min_time_bin;

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}
