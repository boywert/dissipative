/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/fof/fof_bh.c
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
#include <sys/stat.h>
#include <sys/types.h>
#include <gsl/gsl_math.h>
#include <inttypes.h>

#include "../allvars.h"
#include "../proto.h"
#include "../domain.h"
#include "fof.h"
#include "../subfind/subfind.h"


#ifdef FOF
#ifdef BLACK_HOLES

void fof_make_black_holes(void)
{
  int idx, i, j, n, ntot;
  int nexport, nimport, recvTask, level;
  int *import_indices, *export_indices;
  double csnd = 0, rad = 0, dt_courant = 0;
  char msg[200];

  for(n = 0; n < NTask; n++)
    Send_count[n] = 0;

  for(i = 0; i < Ngroups; i++)
    {
      if(Group[i].LenType[1] * All.massDMpart >= (All.Omega0 - All.OmegaBaryon) / All.Omega0 * All.MinFoFMassForNewSeed)
        if(Group[i].LenType[5] == 0)
          {
            if(Group[i].index_maxdens >= 0)
              if(P[Group[i].index_maxdens].Mass != 0 || P[Group[i].index_maxdens].ID != 0)
                Send_count[Group[i].task_maxdens]++;
          }
    }

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = nexport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nexport += Send_count[j];
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  import_indices = mymalloc("import_indices", nimport * sizeof(int));
  export_indices = mymalloc("export_indices", nexport * sizeof(int));

  for(n = 0; n < NTask; n++)
    Send_count[n] = 0;

  for(i = 0; i < Ngroups; i++)
    {
      if(Group[i].LenType[1] * All.massDMpart >= (All.Omega0 - All.OmegaBaryon) / All.Omega0 * All.MinFoFMassForNewSeed)
        if(Group[i].LenType[5] == 0)
          {
            if(Group[i].index_maxdens >= 0)
              if(P[Group[i].index_maxdens].Mass != 0 || P[Group[i].index_maxdens].ID != 0)
                export_indices[Send_offset[Group[i].task_maxdens] + Send_count[Group[i].task_maxdens]++] = Group[i].index_maxdens;
          }
    }

  memcpy(&import_indices[Recv_offset[ThisTask]], &export_indices[Send_offset[ThisTask]], Send_count[ThisTask] * sizeof(int));

  for(level = 1; level < (1 << PTask); level++)
    {
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
        MPI_Sendrecv(&export_indices[Send_offset[recvTask]],
                     Send_count[recvTask] * sizeof(int),
                     MPI_BYTE, recvTask, TAG_FOF_E, &import_indices[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(int), MPI_BYTE, recvTask, TAG_FOF_E, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

  MPI_Allreduce(&nimport, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  mpi_printf("BLACK_HOLES: Making %d new black hole particles\n", ntot);

  All.TotNumBHs += ntot;

  for(n = 0; n < nimport; n++)
    {
      if(P[import_indices[n]].Type != 0)
        {
          sprintf(msg, "Particle n = %d import_indices[n] = %d is not a gas particle\n", n, import_indices[n]);
          terminate(msg);
        }

      P[import_indices[n]].Type = 5;    /* make it a black hole particle */
      P[import_indices[n]].SofteningType = All.SofteningTypeOfPartType[5];

#ifdef INDIVIDUAL_GRAVITY_SOFTENING
      if(((1 << P[import_indices[n]].Type) & (INDIVIDUAL_GRAVITY_SOFTENING)))
        P[import_indices[n]].SofteningType = get_softening_type_from_mass(P[import_indices[n]].Mass);
#endif

      if(NumBHs + 1 > All.MaxPartBHs)
        {
          sprintf(msg, "On Task=%d with NumPart=%d I tried to make too many BH. Sorry, no space left...(All.MaxPartBHs=%d)\n", ThisTask, NumPart, All.MaxPartBHs);
          terminate(msg);
        }

      /* clear all fields to zero in the new BH structure */
      memset(&BHP[NumBHs], 0, sizeof(struct bh_particle_data));

      P[import_indices[n]].AuxDataID = NumBHs;
#ifdef MEASURE_POTMIN_AROUND_BH
      BHP[NumBHs].BH_MinPot = 0;        /* set potential to zero to avoid first repositioning step */
#endif
      BHP[NumBHs].PID = import_indices[n];
      NumBHs++;

#ifdef STELLARAGE
      P[import_indices[n]].StellarAge = All.Time;
#endif
      BPP(import_indices[n]).BH_Mass = All.SeedBlackHoleMass;

      BPP(import_indices[n]).BH_Mdot = 0;
      BPP(import_indices[n]).BH_CumMass_QM = 0;
      BPP(import_indices[n]).BH_CumEgy_QM = 0;
      BPP(import_indices[n]).BH_CountProgs = 1;

      /* the new blackhole is not active for draining gas now. Mdot is zero anyways */
      BPP(import_indices[n]).SwallowID = P[import_indices[n]].ID;

      /* assign Courant timestep of the cell to BH_DtGasNeighbor */

      rad = csnd = 0;

      rad = get_cell_radius(import_indices[n]);
      csnd = get_sound_speed(import_indices[n]);
#ifdef VORONOI_STATIC_MESH
      csnd += sqrt(P[import_indices[n]].Vel[0] * P[import_indices[n]].Vel[0] + P[import_indices[n]].Vel[1] * P[import_indices[n]].Vel[1] + P[import_indices[n]].Vel[2] * P[import_indices[n]].Vel[2]);
#endif

      if(csnd <= 0)
        csnd = 1.0e-30;

      dt_courant = rad / csnd;

#ifdef TREE_BASED_TIMESTEPS
      if(dt_courant > SphP[import_indices[n]].CurrentMaxTiStep)
        dt_courant = SphP[import_indices[n]].CurrentMaxTiStep;
#endif

      BPP(import_indices[n]).BH_DtGasNeighbor = dt_courant;

      BPP(import_indices[n]).BH_Hsml = get_default_softening_of_particletype(5);
#ifdef BH_BUBBLES
      BPP(import_indices[n]).BH_Mass_bubbles = All.SeedBlackHoleMass;
      BPP(import_indices[n]).BH_Mass_ini = All.SeedBlackHoleMass;
#endif

#ifdef DRAINGAS
      BPP(import_indices[n]).NearestDist = MAX_REAL_NUMBER;
#endif

#ifdef BH_FRICTION
      BPP(import_indices[n]).BH_MinPotTime = -1;
      BPP(import_indices[n]).BH_MinPotTime_Previous = -1;
      BPP(import_indices[n]).BH_MinPotCumAvgTime = 0;
#endif

#ifdef BH_NF_RADIO
      BPP(import_indices[n]).BH_RadioEgyFeedback = 0;
      BPP(import_indices[n]).BH_Mdot_quasar = 0;
      BPP(import_indices[n]).BH_Mdot_radio = 0;
      BPP(import_indices[n]).BH_HaloVvir = 0;
      BPP(import_indices[n]).BH_XrayLum = 0;
      BPP(import_indices[n]).BH_RadioLum = 0;
      BPP(import_indices[n]).BH_CumMass_QM = 0;
      BPP(import_indices[n]).BH_CumEgy_QM = 0;
#endif

#ifdef BH_SPIN_EVOLUTION
      double Phi, CosTheta, SinTheta;

	  BPP(import_indices[n]).BH_SpinParameter = All.BHInitialSpin;
	  BPP(import_indices[n]).BlackHoleRadiativeEfficiency = All.BlackHoleRadiativeEfficiency;
	  BPP(import_indices[n]).BH_FlagOngAccEpis = 0;

	  Phi = 2*M_PI*get_random_number();
	  CosTheta = 1-2*get_random_number();
	  SinTheta = sqrt(1-CosTheta*CosTheta);
	  BPP(import_indices[n]).BH_SpinOrientation[0] = SinTheta*cos(Phi);
	  BPP(import_indices[n]).BH_SpinOrientation[1] = SinTheta*sin(Phi);
	  BPP(import_indices[n]).BH_SpinOrientation[2] = CosTheta;
#endif

#ifdef USE_SFR
      Stars_converted++;
#endif

#ifdef VORONOI_DYNAMIC_UPDATE
      voronoi_remove_connection(import_indices[n]);
#endif

      timebin_add_particle(&TimeBinsBHAccretion, import_indices[n], -1, P[import_indices[n]].TimeBinHydro, TimeBinSynchronized[P[import_indices[n]].TimeBinHydro]);
    }

  /* remove the cells that we converted to black holes from the list of active cells */
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      n = TimeBinsHydro.ActiveParticleList[idx];
      if(n < 0)
        continue;

      if(P[n].Type == 5)
        timebin_remove_particle(&TimeBinsHydro, idx, P[n].TimeBinHydro);
    }


  All.TotNumGas -= ntot;

#ifdef MASSIVE_SEEDS
  blackhole_massiveseeds(nimport, import_indices);
#endif

  myfree(export_indices);
  myfree(import_indices);
}
#endif /* BLACK_HOLES */
#endif
