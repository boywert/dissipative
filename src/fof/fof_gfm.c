/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/fof/fof_gfm.c
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

/*! \file fof.c
 *  \brief parallel FoF group finder
 */

#ifdef FOF



#if (defined(GFM_WINDS_VARIABLE) && (GFM_WINDS_VARIABLE==0)) || defined(GFM_BIPOLAR_WINDS) || defined(GFM_AGN_RADIATION) || defined(MASSIVE_SEEDS_MERGER)  || defined(BH_NF_RADIO)

struct group_mass_MinID
{
  MyIDType MinID;
  MyFloat mass;
#ifdef GFM_BIPOLAR_WINDS
#if (GFM_BIPOLAR_WINDS == 3)
  MyFloat DensGasAngMomentum[3];
#else
  MyFloat Vel[3];
  MyFloat GravAcc[3];
#endif
#endif
#if defined(BH_NF_RADIO)
  MyIDType ID_Min_BH_Potential;
#endif
#ifdef BH_NF_RADIO
  MyFloat XrayLum;
#endif
};

int compare_group_mass_ID(const void *a, const void *b)
{
  if(((struct group_mass_MinID *) a)->MinID < (((struct group_mass_MinID *) b)->MinID))
    return -1;

  if(((struct group_mass_MinID *) a)->MinID > (((struct group_mass_MinID *) b)->MinID))
    return +1;

  return 0;
}

void fof_assign_HostHaloMass(void)      /* assigns mass of host FoF group to SphP[].w.HostHaloMass for SPH particles and/or to BHP[].HostHaloMass for BH particles */
{
  int i, j, k, start, lenloc, nimport;
  struct group_mass_MinID *required_groups, *groups_to_export;
#ifdef GFM_BIPOLAR_WINDS
  int l;
#endif

  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;
  for(i = 0; i < NgroupsExt; i++)       /* loop over all groups for which at least one particle is on this task */
    Send_count[FOF_GList[i].MinIDTask]++;       /* its FoF group properties are stored on Task = MinIDTask */

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  qsort(FOF_GList, NgroupsExt, sizeof(struct fof_group_list), fof_compare_FOF_GList_MinIDTask_MinID);

  required_groups = (struct group_mass_MinID *) mymalloc("required_groups", NgroupsExt * sizeof(struct group_mass_MinID));

  for(i = 0; i < NgroupsExt; i++)
    required_groups[i].MinID = FOF_GList[i].MinID;

  groups_to_export = (struct group_mass_MinID *) mymalloc("groups_to_export", nimport * sizeof(struct group_mass_MinID));

  int slen = sizeof(struct group_mass_MinID);
  for(i = 0; i < NTask; i++)
    {
      Send_count[i] *= slen;
      Send_offset[i] *= slen;
      Recv_count[i] *= slen;
      Recv_offset[i] *= slen;
    }

  MPI_Alltoallv(required_groups, Send_count, Send_offset, MPI_BYTE, groups_to_export, Recv_count, Recv_offset, MPI_BYTE, MPI_COMM_WORLD);

  for(i = 0; i < NTask; i++)
    {
      Send_count[i] /= slen;
      Send_offset[i] /= slen;
      Recv_count[i] /= slen;
      Recv_offset[i] /= slen;
    }


  for(j = 0, start = 0; j < NTask; j++)
    {
      i = 0;
      k = 0;

      while(i < Recv_count[j] && k < Ngroups)
        {
          if(groups_to_export[start].MinID == Group[k].MinID)
            {
              groups_to_export[start].mass = Group[k].Mass;
#ifdef GFM_BIPOLAR_WINDS
              for(l = 0; l < 3; l++)
                {
#if (GFM_BIPOLAR_WINDS == 3)
                  groups_to_export[start].DensGasAngMomentum[l] = Group[k].DensGasAngMomentum[l];
#else
                  groups_to_export[start].Vel[l] = Group[k].Vel[l];
                  groups_to_export[start].GravAcc[l] = Group[k].GravAcc[l];
#endif
                }
#endif
#if defined(BH_NF_RADIO)
              groups_to_export[start].ID_Min_BH_Potential = Group[k].ID_Min_BH_Potential;
#endif
#ifdef BH_NF_RADIO
              groups_to_export[start].XrayLum = Group[k].XrayLum;
#endif
              i++;
              k++;
              start++;
            }
          else
            k++;
        }
    }
  if(start != nimport)
    terminate("start != nimport");

  for(i = 0; i < NTask; i++)
    {
      Send_count[i] *= slen;
      Send_offset[i] *= slen;
      Recv_count[i] *= slen;
      Recv_offset[i] *= slen;
    }

  MPI_Alltoallv(groups_to_export, Recv_count, Recv_offset, MPI_BYTE, required_groups, Send_count, Send_offset, MPI_BYTE, MPI_COMM_WORLD);

  for(i = 0; i < NTask; i++)
    {
      Send_count[i] /= slen;
      Send_offset[i] /= slen;
      Recv_count[i] /= slen;
      Recv_offset[i] /= slen;
    }


  myfree(groups_to_export);

  qsort(required_groups, NgroupsExt, sizeof(struct group_mass_MinID), compare_group_mass_ID);

  for(i = 0; i < NumGas; i++)
    {
#if defined(GFM_WINDS_VARIABLE) && (GFM_WINDS_VARIABLE==0)
      SphP[i].w.HostHaloMass = 0;
#endif
#ifdef GFM_BIPOLAR_WINDS
      for(j = 0; j < 3; j++)
        {
#if (GFM_BIPOLAR_WINDS == 3)
          SphP[i].DensGasAngMomentum[j] = 0;
#else
          SphP[i].GroupVel[j] = 0;
          SphP[i].GroupGravAcc[j] = 0;
#endif
        }
#endif
    }

#ifdef BH_NF_RADIO
  /* clear radio-mode activity for all BHs. At most one BH per halo will be "switched on" again in this routine. */
  for(i = 0; i < NumBHs; i++)
    {
      BHP[i].BH_HaloVvir = 0;
      BHP[i].BH_XrayLum = 0;
      BHP[i].BH_Mdot_radio = 0;
      BHP[i].BH_RadioLum = 0;
    }
#endif

#if defined(BLACK_HOLES) && (defined(GFM_AGN_RADIATION) || defined(MASSIVE_SEEDS_MERGER))
  for(i = 0; i < NumBHs; i++)
    BHP[i].HostHaloMass = 0.0;
#endif


  for(i = 0, start = 0; i < NgroupsExt; i++)
    {
      while(FOF_PList[start].MinID < required_groups[i].MinID)
        {
          start++;
          if(start > NumPart)
            terminate("start > NumPart");
        }

      if(FOF_PList[start].MinID != required_groups[i].MinID)
        terminate("FOF_PList[start].MinID != required_groups[i].MinID");

      for(lenloc = 0; start + lenloc < NumPart;)
        if(FOF_PList[start + lenloc].MinID == required_groups[i].MinID)
          {
            if(P[FOF_PList[start + lenloc].Pindex].Type == 0)
              {
#if defined(GFM_WINDS_VARIABLE) && (GFM_WINDS_VARIABLE==0)
                SphP[FOF_PList[start + lenloc].Pindex].w.HostHaloMass = required_groups[i].mass;
#endif
#ifdef GFM_BIPOLAR_WINDS
                for(j = 0; j < 3; j++)
                  {
#if (GFM_BIPOLAR_WINDS == 3)
                    SphP[FOF_PList[start + lenloc].Pindex].DensGasAngMomentum[j] = required_groups[i].DensGasAngMomentum[j];
#else
                    SphP[FOF_PList[start + lenloc].Pindex].GroupVel[j] = required_groups[i].Vel[j];
                    SphP[FOF_PList[start + lenloc].Pindex].GroupGravAcc[j] = required_groups[i].GravAcc[j];
#endif
                  }
#endif
              }

#if defined(BH_NF_RADIO)
            if(P[FOF_PList[start + lenloc].Pindex].Type == 5)
              BPP(FOF_PList[start + lenloc].Pindex).ID_Min_BH_Potential = required_groups[i].ID_Min_BH_Potential;
#endif

#ifdef BH_NF_RADIO
            if(P[FOF_PList[start + lenloc].Pindex].Type == 5)
              if(P[FOF_PList[start + lenloc].Pindex].ID == BPP(FOF_PList[start + lenloc].Pindex).ID_Min_BH_Potential)
                {
                  double vvir = pow(10 * All.G * All.cf_H * required_groups[i].mass, 1.0 / 3);

                  BPP(FOF_PList[start + lenloc].Pindex).BH_HaloVvir = vvir;
                  BPP(FOF_PList[start + lenloc].Pindex).BH_XrayLum = required_groups[i].XrayLum;

                  BPP(FOF_PList[start + lenloc].Pindex).BH_RadioLum = (All.Hubble / All.cf_H) * blackhole_get_radio_efficiency(vvir) * BPP(FOF_PList[start + lenloc].Pindex).BH_XrayLum;

                  BPP(FOF_PList[start + lenloc].Pindex).BH_Mdot_radio = blackhole_get_mdot_radio_from_radiolum(FOF_PList[start + lenloc].Pindex);
                }
#endif


#if defined (BLACK_HOLES) && (defined(GFM_AGN_RADIATION) || defined(MASSIVE_SEEDS_MERGER))
            if(P[FOF_PList[start + lenloc].Pindex].Type == 5)
              BPP(FOF_PList[start + lenloc].Pindex).HostHaloMass = required_groups[i].mass;
#endif
            lenloc++;
          }
        else
          break;

      start += lenloc;
    }

  myfree(required_groups);

  qsort(FOF_GList, NgroupsExt, sizeof(struct fof_group_list), fof_compare_FOF_GList_MinID);     /* restore original order */
}
#endif

#endif
