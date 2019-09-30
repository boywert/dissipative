/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/cosmic_rays.c
 * \date        11/2013
 * \author      R. Pakmor & C. Pfrommer
 * \brief       A simple cosmic ray implementation
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 */

#include "../allvars.h"
#include "../proto.h"
#include "string.h"

#ifdef COSMIC_RAYS_SN_INJECTION
/** \brief Inject cosmic rays for all newly formed stars
 *
 *  We assume instantanious deposition of cosmic rays from exploding
 *  core-collapse supernovae into the surrounding medium when a star
 *  is formed in the simulations. The amount of energy deposited in
 *  cosmic rays is assumed to scale linearly with the amount of mass
 *  that is converted into stars.
 *    
 */

int *NgbCount;
MyFloat *NumNgb;
MyFloat *NormSph;
double CREnergy_Injected;

static void cosmic_rays_inject(int target, int mode, int thread_id);

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble Pos[3];
  MyFloat Hsml;
  MyFloat Mass;
  MyFloat NormSph;

  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

 /* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in * in, int i, int firstnode)
{
  in->Pos[0] = P[i].Pos[0];
  in->Pos[1] = P[i].Pos[1];
  in->Pos[2] = P[i].Pos[2];

  in->Mass = P[i].Mass;
  in->Hsml = P[i].Hsml;
  in->NormSph = NormSph[i];

  in->Firstnode = firstnode;
}


 /* local data structure that holds results acquired on remote processors */
typedef struct
{
  MyFloat EnergyInjected;
} data_out;

static data_out *DataResult, *DataOut;


 /* routine to store or combine result data */
static void out2particle(data_out * out, int i, int mode)
{
  CREnergy_Injected += out->EnergyInjected;
}

#include "../generic_comm_helpers2.h"


static void kernel_local(void)
{
  int idx;
#ifdef GENERIC_ASYNC
  int flag = 0;
#endif
  /* do local particles */
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
              if(generic_polling_primary(count, TimeBinsGravity.NActiveParticles))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        idx = NextParticle++;

        if(idx >= TimeBinsGravity.NActiveParticles)
          break;

        int i = TimeBinsGravity.ActiveParticleList[idx];
        if(i < 0)
          continue;

        if(P[i].Type != 4)
          continue;

#ifdef GFM
	if(STP(i).BirthTime < 0.)
	  continue;
#endif

	/* don't inject into low res star particles */
	if(P[i].Mass > 10. * All.TargetGasMass)
	  continue;

        if(P[i].CRInjection == 1)
          {
            cosmic_rays_inject(i, MODE_LOCAL_PARTICLES, threadid);
          }
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

        cosmic_rays_inject(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}


void deposit_cosmic_rays_from_supernovae(void)
{
  int idx, i;

  NgbCount = (int *) mymalloc("NgbCount", NumPart * sizeof(int));
  NumNgb = (MyFloat *) mymalloc("NumNgb", NumPart * sizeof(MyFloat));
  NormSph = (MyFloat *) mymalloc("NormSph", NumPart * sizeof(MyFloat));

  cosmic_rays_find_ngbs(NgbCount, NumNgb, NormSph);

  double MassStars = 0;
  CREnergy_Injected = 0;

  generic_set_MaxNexport();

  generic_comm_pattern(TimeBinsGravity.NActiveParticles, kernel_local, kernel_imported);

  myfree(NormSph);
  myfree(NumNgb);
  myfree(NgbCount);

  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0 || P[i].Type != 4)
        continue;

      if(P[i].CRInjection == 1)
        {
          P[i].CRInjection = 0;
          MassStars += P[i].Mass;
        }
    }

  double totCREnergy_Injected, totMassStars;
  MPI_Reduce(&CREnergy_Injected, &totCREnergy_Injected, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&MassStars, &totMassStars, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0 && totMassStars > 0)
    {
      mpi_printf("COSMIC_RAYS: Injected a total of %g erg in cosmic rays for a total of %g solar masses of newly formed stars.\n",
                 totCREnergy_Injected * All.UnitEnergy_in_cgs, totMassStars * All.UnitMass_in_g / SOLAR_MASS);
    }
}

void cosmic_rays_inject(int target, int mode, int thread_id)
{
  int j, n;
  int numnodes, *firstnode;
  double h, h2;
  double wk, hinv, hinv3;
  double u;
  double dx, dy, dz, r, r2;
  MyDouble *pos, normSph;
  double mass, e_inj;
#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif
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
  h = in->Hsml;
  mass = in->Mass;
  normSph = in->NormSph;
  h2 = h * h;

  hinv = 1.0 / h;
#ifndef TWODIMS
  hinv3 = hinv * hinv * hinv;
#else
  hinv3 = hinv * hinv / boxSize_Z;
#endif

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, thread_id, numnodes, firstnode);

  for(n = 0; n < nfound; n++)
    {
      j = Thread[thread_id].Ngblist[n];

      if(P[j].Mass > 0 && P[j].ID != 0)
        {
          dx = NGB_PERIODIC_LONG_X(pos[0] - P[j].Pos[0]);
          dy = NGB_PERIODIC_LONG_Y(pos[1] - P[j].Pos[1]);
          dz = NGB_PERIODIC_LONG_Z(pos[2] - P[j].Pos[2]);

          r2 = dx * dx + dy * dy + dz * dz;

          if(r2 < h2)
            {
              r = sqrt(r2);
              u = r * hinv;

              if(u < 0.5)
                wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
              else
                wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

              double weight_fac = SphP[j].Volume * wk / normSph;

              e_inj = All.CREnergyInputPerSolarMassOfStarFormation / All.UnitEnergy_in_cgs * mass * All.UnitMass_in_g / SOLAR_MASS * weight_fac;
              SphP[j].CR_Energy += e_inj * All.cf_atime;
              All.TotalCREnergyInjected += e_inj;
            }
        }
    }

  out.EnergyInjected = 0;

  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;
}
#endif
