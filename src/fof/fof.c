/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/fof/fof.c
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


static MyIDType *MinID;
static int *Head, *Len, *Next, *Tail, *MinIDTask;

/* If called with num == -1 as argument, only FOF is carried out and no group catalogs are saved to disk.
 * If num >= 0, the code will store the group/subgroup catalogs, and bring the particles into output order.
 * In this case, the calling routine (which is normally savepositions()) will need to free PS[] and bring
 * the particles back into the original order, as well as reestablished the mesh.
 */
void fof_fof(int num)
{
  int i, start, lenloc, largestgroup;
  double t0, t1, cputime;

  TIMER_START(CPU_FOF);

  mpi_printf("FOF: Begin to compute FoF group catalogue...  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

  if(num >= 0 && RestartFlag != 3 && RestartFlag != 6)
    {
      /* let's discard an existing mesh - we do this here to reduce the peak memory usage, even at the price of
       * having to recreate it later */
      free_mesh();
    }

  if(RestartFlag != 6)
    {
      ngb_treefree();

      domain_free();
    }

#ifdef ADD_GROUP_PROPERTIES
  int ngroups_cat = get_number_of_groups_in_catalogue(num);
  int nsubgroups_cat = get_number_of_subgroups_in_catalogue(num);
  subfind_add_grp_props_read_catalogue(num, ngroups_cat, nsubgroups_cat);
#endif

  domain_Decomposition();

  ngb_treeallocate();
  ngb_treebuild(NumGas);

  /* check */
  for(i = 0; i < NumPart; i++)
    if((P[i].Mass == 0 && P[i].ID == 0) || (P[i].Type == 4 && P[i].Mass == 0))
      terminate("this should not happen");

  /* this structure will hold auxiliary information for each particle, needed only during group finding */
  PS = (struct subfind_data *) mymalloc_movable(&PS, "PS", All.MaxPart * sizeof(struct subfind_data));

  memset(PS, 0, NumPart * sizeof(struct subfind_data));

  /* First, we save the original location of the particles, in order to be able to revert to this layout later on */
  for(i = 0; i < NumPart; i++)
    {
      PS[i].OriginTask = ThisTask;
      PS[i].OriginIndex = i;
    }

  fof_OldMaxPart = All.MaxPart;
  fof_OldMaxPartSph = All.MaxPartSph;
#ifdef GFM
  fof_OldMaxPartStar = All.MaxPartStar;
#endif
#ifdef BLACK_HOLES
  fof_OldMaxPartBHs = All.MaxPartBHs;
#endif
#ifdef DUST_LIVE
  fof_OldMaxPartDust = All.MaxPartDust;
#endif

  LinkL = fof_get_comoving_linking_length();

  mpi_printf("FOF: Comoving linking length: %g    (presently allocated=%g MB)\n", LinkL, AllocatedBytes / (1024.0 * 1024.0));

  MinID = (MyIDType *) mymalloc("MinID", NumPart * sizeof(MyIDType));
  MinIDTask = (int *) mymalloc("MinIDTask", NumPart * sizeof(int));

  Head = (int *) mymalloc("Head", NumPart * sizeof(int));
  Len = (int *) mymalloc("Len", NumPart * sizeof(int));
  Next = (int *) mymalloc("Next", NumPart * sizeof(int));
  Tail = (int *) mymalloc("Tail", NumPart * sizeof(int));

#ifdef HIERARCHICAL_GRAVITY
  timebin_make_list_of_active_particles_up_to_timebin(&TimeBinsGravity, All.HighestOccupiedTimeBin);
#endif

  force_treeallocate(NumPart, All.MaxPart);
  force_treebuild(NumPart, 0, 1, All.HighestOccupiedTimeBin);

#if defined(SUBFIND) || (defined(GFM_WINDS_VARIABLE) && (GFM_WINDS_VARIABLE==1)) || defined(GFM_WINDS_LOCAL)
  subfind_density_hsml_guess();
#endif

  /* initialize link-lists */
  for(i = 0; i < NumPart; i++)
    {
      Head[i] = Tail[i] = i;
      Len[i] = 1;
      Next[i] = -1;
      MinID[i] = P[i].ID;
      MinIDTask[i] = ThisTask;
    }

#ifndef ADD_GROUP_PROPERTIES

  /* call routine to find primary groups */
  cputime = fof_find_groups(MinID, Head, Len, Next, Tail, MinIDTask);
  mpi_printf("FOF: group finding took = %g sec\n", cputime);

#endif

#ifdef FOF_SECONDARY_LINK_TARGET_TYPES
  myfree(Father);
  myfree(Nextnode);
#ifdef BLACK_HOLES
  myfree(Tree_AuxBH_Points);
#endif
  myfree(Tree_Points);

  /* now rebuild the tree with all the types selected as secondary link targets */
  force_treebuild(NumPart, 0, 2, All.HighestOccupiedTimeBin);
#endif

#ifdef HIERARCHICAL_GRAVITY
  timebin_make_list_of_active_particles_up_to_timebin(&TimeBinsGravity, All.HighestActiveTimeBin);
#endif

#ifndef ADD_GROUP_PROPERTIES

  /* call routine to attach secondary particles/cells to primary groups */
  cputime = fof_find_nearest_dmparticle(MinID, Head, Len, Next, Tail, MinIDTask);


  mpi_printf("FOF: attaching gas and star particles to nearest dm particles took = %g sec\n", cputime);

#endif

  /* calculate velocity dispersion etc. */
#if (defined(GFM_WINDS_VARIABLE) && (GFM_WINDS_VARIABLE==1)) || defined(GFM_WINDS_LOCAL) || (defined(GFM_VARIABLE_IMF) && (GFM_VARIABLE_IMF==0))
  if(num == -2)
    {
      /* at this point, TotNgroups/PS[].GrNr are not defined yet. The following selects just gas particles
         for processing in subfind_density() */

      TotNgroups = 1;
      for(i = 0; i < NumPart; i++)
        {
          if(P[i].Type == 0)
            PS[i].GrNr = 0;
          else
            PS[i].GrNr = 1;
        }

      cputime = subfind_density(0);
      mpi_printf("GFM_WINDS_VARIABLE/GFM_WINDS_LOCAL/GFM_VARIABLE_IMF: velocity dispersion calculation took %g sec\n", cputime);
      for(i = 0; i < NumPart; i++)
        {
          if(P[i].Type == 0)
            SphP[i].w.DMVelDisp = PS[i].SubfindVelDisp;
        }
    }
#endif

  myfree(Father);
  myfree(Nextnode);
#ifdef BLACK_HOLES
  myfree(Tree_AuxBH_Points);
#endif
  myfree(Tree_Points);
  force_treefree();

  myfree(Tail);
  myfree(Next);
  myfree(Len);

  t0 = second();

  FOF_PList = (struct fof_particle_list *) mymalloc_movable(&FOF_PList, "FOF_PList", NumPart * sizeof(struct fof_particle_list));

  for(i = 0; i < NumPart; i++)
    {
#ifdef ADD_GROUP_PROPERTIES
      FOF_PList[i].MinID = P[i].MinID;
      FOF_PList[i].MinIDTask = P[i].MinIDTask;
      FOF_PList[i].OriginalGrNr = P[i].OriginalGrNr;
#else
      FOF_PList[i].MinID = MinID[Head[i]];
      FOF_PList[i].MinIDTask = MinIDTask[Head[i]];
#endif
      FOF_PList[i].Pindex = i;
    }

  myfree_movable(Head);
  myfree_movable(MinIDTask);
  myfree_movable(MinID);





  FOF_GList = (struct fof_group_list *) mymalloc_movable(&FOF_GList, "FOF_GList", sizeof(struct fof_group_list) * NumPart);

  fof_compile_catalogue();

  t1 = second();
  mpi_printf("FOF: compiling local group data and catalogue took = %g sec\n", timediff(t0, t1));

  MPI_Allreduce(&Ngroups, &TotNgroups, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  sumup_large_ints(1, &Nids, &TotNids);

  if(TotNgroups > 0)
    {
      int largestloc = 0;

      for(i = 0; i < NgroupsExt; i++)
        if(FOF_GList[i].LocCount + FOF_GList[i].ExtCount > largestloc)
          largestloc = FOF_GList[i].LocCount + FOF_GList[i].ExtCount;
      MPI_Allreduce(&largestloc, &largestgroup, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    }
  else
    largestgroup = 0;

  mpi_printf("FOF: Total number of FOF groups with at least %d particles: %d\n", FOF_GROUP_MIN_LEN, TotNgroups);
  mpi_printf("FOF: Largest FOF group has %d particles.\n", largestgroup);
  mpi_printf("FOF: Total number of particles in FOF groups: %lld\n", TotNids);

  t0 = second();

  MaxNgroups = 2 * imax(NgroupsExt, TotNgroups / NTask + 1);

  Group = (struct group_properties *) mymalloc_movable(&Group, "Group", sizeof(struct group_properties) * MaxNgroups);

  mpi_printf("FOF: group properties are now allocated.. (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

  for(i = 0, start = 0; i < NgroupsExt; i++)
    {
      while(FOF_PList[start].MinID < FOF_GList[i].MinID)
        {
          start++;
          if(start > NumPart)
            terminate("start > NumPart");
        }

      if(FOF_PList[start].MinID != FOF_GList[i].MinID)
        terminate("ID mismatch");

      for(lenloc = 0; start + lenloc < NumPart;)
        if(FOF_PList[start + lenloc].MinID == FOF_GList[i].MinID)
          lenloc++;
        else
          break;

      Group[i].MinID = FOF_GList[i].MinID;
      Group[i].MinIDTask = FOF_GList[i].MinIDTask;

      fof_compute_group_properties(i, start, lenloc);

      start += lenloc;
    }

  fof_exchange_group_data();


  fof_finish_group_properties();

  t1 = second();
  mpi_printf("FOF: computation of group properties took = %g sec\n", timediff(t0, t1));


#if (GFM_BIPOLAR_WINDS == 3)
  fof_spin_measurement();
#endif


#if (defined(GFM_WINDS_VARIABLE) && (GFM_WINDS_VARIABLE==0)) || defined(GFM_BIPOLAR_WINDS) || defined(GFM_AGN_RADIATION) || defined(MASSIVE_SEEDS_MERGER) || defined(BH_NF_RADIO)
  if(num < 0)
    fof_assign_HostHaloMass();
#endif

#ifdef BLACK_HOLES
  if(num < 0)
    fof_make_black_holes();
#endif

  fof_assign_group_numbers();

  mpi_printf("FOF: Finished computing FoF groups.  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));



  myfree_movable(FOF_GList);
  myfree_movable(FOF_PList);

#ifdef SUBFIND
  if(num >= 0)
    {
      TIMER_STOP(CPU_FOF);

      subfind(num);

      TIMER_START(CPU_FOF);
    }
#else
  Nsubgroups = 0;
  TotNsubgroups = 0;
  if(num >= 0)
    {
      TIMER_STOP(CPU_FOF);
      TIMER_START(CPU_SNAPSHOT);

      fof_save_groups(num);

      TIMER_STOP(CPU_SNAPSHOT);
      TIMER_START(CPU_FOF);
    }
#endif

  myfree_movable(Group);

  mpi_printf("FOF: All FOF related work finished.  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

#ifndef FOF_STOREIDS
  if(num >= 0)
    {
      TIMER_STOP(CPU_FOF);
      TIMER_START(CPU_SNAPSHOT);

      /* now distribute the particles into output order */
      t0 = second();
      fof_prepare_output_order();
      fof_subfind_exchange(MPI_COMM_WORLD);     /* distribute particles such that FOF groups will appear in coherent way in snapshot files */
      t1 = second();
      mpi_printf("FOF: preparing output order of particles took %g sec\n", timediff(t0, t1));

      TIMER_STOP(CPU_SNAPSHOT);
      TIMER_START(CPU_FOF);
    }
  else
    myfree(PS);
#else
  myfree(PS);
#endif

  TIMER_STOP(CPU_FOF);
}


void fof_prepare_output_order(void)
{
  int i, off, ntype[NTYPES];

  struct data_aux_sort *aux_sort = (struct data_aux_sort *) mymalloc("aux_sort", sizeof(struct data_aux_sort) * NumPart);

  for(i = 0; i < NTYPES; i++)
    ntype[i] = 0;

  for(i = 0; i < NumPart; i++)
    {
      aux_sort[i].OriginTask = ThisTask;
      aux_sort[i].OriginIndex = i;
      aux_sort[i].GrNr = PS[i].GrNr;
#ifdef SUBFIND
      aux_sort[i].SubNr = PS[i].SubNr;
      aux_sort[i].DM_BindingEnergy = PS[i].BindingEnergy;
#endif
      aux_sort[i].Type = P[i].Type;
#ifdef FOF_FUZZ_SORT_BY_NEAREST_GROUP
      aux_sort[i].key = PS[i].GroupNr;
#else
      aux_sort[i].ID = P[i].ID;
#endif
#if defined(RECOMPUTE_POTENTIAL_IN_SNAPSHOT) || defined(CALCULATE_QUANTITIES_IN_POSTPROCESS)
      aux_sort[i].FileOrder = P[i].FileOrder;
#endif

      ntype[P[i].Type]++;
    }

  qsort(aux_sort, NumPart, sizeof(struct data_aux_sort), fof_compare_aux_sort_Type);

  if(RestartFlag == 18 || RestartFlag == 19)
    {
#if defined(RECOMPUTE_POTENTIAL_IN_SNAPSHOT) || defined(CALCULATE_QUANTITIES_IN_POSTPROCESS)
      for(i = 0, off = 0; i < NTYPES; off += ntype[i], i++)
	parallel_sort(aux_sort + off, ntype[i], sizeof(struct data_aux_sort), fof_compare_aux_sort_FileOrder);
#endif
    }
  else
    {
      for(i = 0, off = 0; i < NTYPES; off += ntype[i], i++)
	parallel_sort(aux_sort + off, ntype[i], sizeof(struct data_aux_sort), fof_compare_aux_sort_GrNr);
    }

  for(i = 0; i < NumPart; i++)
    {
      aux_sort[i].TargetTask = ThisTask;
      aux_sort[i].TargetIndex = i;
    }


  /* now bring back into starting order */
  parallel_sort(aux_sort, NumPart, sizeof(struct data_aux_sort), fof_compare_aux_sort_OriginTask_OriginIndex);

  for(i = 0; i < NumPart; i++)
    {
      PS[i].TargetTask = aux_sort[i].TargetTask;
      PS[i].TargetIndex = aux_sort[i].TargetIndex;
    }

  myfree(aux_sort);
}

/* calculate linking length based on mean particle separation */
double fof_get_comoving_linking_length(void)
{
  int i, ndm;
  long long ndmtot;
  double mass, masstot, rhodm;

  for(i = 0, ndm = 0, mass = 0; i < NumPart; i++)
    if(((1 << P[i].Type) & (FOF_PRIMARY_LINK_TYPES)))
      {
        ndm++;
        mass += P[i].Mass;
      }
  sumup_large_ints(1, &ndm, &ndmtot);
  MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  rhodm = (All.Omega0 - All.OmegaBaryon) * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

  return FOF_LINKLENGTH * pow(masstot / ndmtot / rhodm, 1.0 / 3);
}




void fof_compile_catalogue(void)
{
  int i, j, start, nimport, ngrp, recvTask;
  struct fof_group_list *get_FOF_GList;

  /* sort according to MinID */
  mysort(FOF_PList, NumPart, sizeof(struct fof_particle_list), fof_compare_FOF_PList_MinID);

  for(i = 0; i < NumPart; i++)
    {
      FOF_GList[i].MinID = FOF_PList[i].MinID;
      FOF_GList[i].MinIDTask = FOF_PList[i].MinIDTask;
      if(FOF_GList[i].MinIDTask == ThisTask)
        {
          FOF_GList[i].LocCount = 1;
          FOF_GList[i].ExtCount = 0;
        }
      else
        {
          FOF_GList[i].LocCount = 0;
          FOF_GList[i].ExtCount = 1;
        }

#ifdef ADD_GROUP_PROPERTIES
      FOF_GList[i].OriginalGrNr = FOF_PList[i].OriginalGrNr;
#endif

#ifdef TRACER_PARTICLE
      /* count tracer particles independently so we can enforce FOF_GROUP_MIN_LEN against only non-tracers */
      if(P[FOF_PList[i].Pindex].Type == TRACER_PARTICLE)
        {
          if(FOF_GList[i].MinIDTask == ThisTask)
            {
              FOF_GList[i].LocTrCount = 1;
              FOF_GList[i].ExtTrCount = 0;
            }
          else
            {
              FOF_GList[i].LocTrCount = 0;
              FOF_GList[i].ExtTrCount = 1;
            }
        }
      else
        {
          FOF_GList[i].LocTrCount = 0;
          FOF_GList[i].ExtTrCount = 0;
        }
#endif
    }

  /* eliminate duplicates in FOF_GList with respect to MinID */

  if(NumPart)
    NgroupsExt = 1;
  else
    NgroupsExt = 0;

  for(i = 1, start = 0; i < NumPart; i++)
    {
      if(FOF_GList[i].MinID == FOF_GList[start].MinID)
        {
          FOF_GList[start].LocCount += FOF_GList[i].LocCount;
          FOF_GList[start].ExtCount += FOF_GList[i].ExtCount;
#ifdef TRACER_PARTICLE
          FOF_GList[start].LocTrCount += FOF_GList[i].LocTrCount;
          FOF_GList[start].ExtTrCount += FOF_GList[i].ExtTrCount;
#endif
        }
      else
        {
          start = NgroupsExt;
          FOF_GList[start] = FOF_GList[i];
          NgroupsExt++;
        }
    }

  /* sort the remaining ones according to task */
  mysort(FOF_GList, NgroupsExt, sizeof(struct fof_group_list), fof_compare_FOF_GList_MinIDTask);

  /* count how many we have of each task */
  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;
  for(i = 0; i < NgroupsExt; i++)
    Send_count[FOF_GList[i].MinIDTask]++;

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      if(j == ThisTask)         /* we will not exchange the ones that are local */
        Recv_count[j] = 0;
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  get_FOF_GList = (struct fof_group_list *) mymalloc("get_FOF_GList", nimport * sizeof(struct fof_group_list));

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the group info */
              MPI_Sendrecv(&FOF_GList[Send_offset[recvTask]],
                           Send_count[recvTask] * sizeof(struct fof_group_list), MPI_BYTE,
                           recvTask, TAG_DENS_A,
                           &get_FOF_GList[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct fof_group_list), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  for(i = 0; i < nimport; i++)
    get_FOF_GList[i].MinIDTask = i;


  /* sort the groups according to MinID */
  mysort(FOF_GList, NgroupsExt, sizeof(struct fof_group_list), fof_compare_FOF_GList_MinID);
  mysort(get_FOF_GList, nimport, sizeof(struct fof_group_list), fof_compare_FOF_GList_MinID);

  /* merge the imported ones with the local ones */
  for(i = 0, start = 0; i < nimport; i++)
    {
      while(FOF_GList[start].MinID < get_FOF_GList[i].MinID)
        {
          start++;
          if(start >= NgroupsExt)
            terminate("start >= NgroupsExt");
        }

      if(get_FOF_GList[i].LocCount != 0)
        terminate("start >= NgroupsExt");

      if(FOF_GList[start].MinIDTask != ThisTask)
        terminate("FOF_GList[start].MinIDTask != ThisTask");

      if(FOF_GList[start].MinID != get_FOF_GList[i].MinID)
        terminate("FOF_GList[start].MinID != get_FOF_GList[i].MinID start=%d i=%d FOF_GList[start].MinID=%llu get_FOF_GList[i].MinID=%llu\n", start,
                  i, (long long) FOF_GList[start].MinID, (long long) get_FOF_GList[i].MinID);

      FOF_GList[start].ExtCount += get_FOF_GList[i].ExtCount;
#ifdef TRACER_PARTICLE
      FOF_GList[start].ExtTrCount += get_FOF_GList[i].ExtTrCount;

      if(FOF_GList[start].LocTrCount > FOF_GList[start].LocCount || FOF_GList[start].ExtTrCount > FOF_GList[start].ExtCount)
        terminate("Error: FOF_GList tracer counts exceed total counts!");
#endif
    }

  /* copy the size information back into the list, to inform the others */
  for(i = 0, start = 0; i < nimport; i++)
    {
      while(FOF_GList[start].MinID < get_FOF_GList[i].MinID)
        {
          start++;
          if(start >= NgroupsExt)
            terminate("start >= NgroupsExt");
        }

      get_FOF_GList[i].ExtCount = FOF_GList[start].ExtCount;
      get_FOF_GList[i].LocCount = FOF_GList[start].LocCount;
#ifdef TRACER_PARTICLE
      get_FOF_GList[i].ExtTrCount = FOF_GList[start].ExtTrCount;
      get_FOF_GList[i].LocTrCount = FOF_GList[start].LocTrCount;
#endif
    }

  /* sort the imported/exported list according to MinIDTask */
  mysort(get_FOF_GList, nimport, sizeof(struct fof_group_list), fof_compare_FOF_GList_MinIDTask);
  mysort(FOF_GList, NgroupsExt, sizeof(struct fof_group_list), fof_compare_FOF_GList_MinIDTask);


  for(i = 0; i < nimport; i++)
    get_FOF_GList[i].MinIDTask = ThisTask;

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the group info */
              MPI_Sendrecv(&get_FOF_GList[Recv_offset[recvTask]],
                           Recv_count[recvTask] * sizeof(struct fof_group_list), MPI_BYTE,
                           recvTask, TAG_DENS_A,
                           &FOF_GList[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct fof_group_list), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  myfree(get_FOF_GList);

  /* eliminate all groups that are too small, and count local groups */
  for(i = 0, Ngroups = 0, Nids = 0; i < NgroupsExt; i++)
    {
      if(FOF_GList[i].LocCount + FOF_GList[i].ExtCount
#ifdef TRACER_PARTICLE
         - FOF_GList[i].LocTrCount - FOF_GList[i].ExtTrCount
#endif
         < FOF_GROUP_MIN_LEN)
        {
          FOF_GList[i] = FOF_GList[NgroupsExt - 1];
          NgroupsExt--;
          i--;
        }
      else
        {
          if(FOF_GList[i].MinIDTask == ThisTask)
            {
              Ngroups++;
              Nids += FOF_GList[i].LocCount + FOF_GList[i].ExtCount;
            }
        }
    }

  /* sort the group list according to MinID */
  mysort(FOF_GList, NgroupsExt, sizeof(struct fof_group_list), fof_compare_FOF_GList_MinID);
}


void fof_assign_group_numbers(void)
{
  int i, j, ngr, start, lenloc;
  long long totNids;
  double t0, t1;

  mpi_printf("FOF: start assigning group numbers\n");

  t0 = second();

  /* assign group numbers (at this point, both Group and FOF_GList are sorted by MinID) */
  for(i = 0; i < NgroupsExt; i++)
    {
      FOF_GList[i].LocCount += FOF_GList[i].ExtCount;   /* total length */
      FOF_GList[i].ExtCount = ThisTask; /* original task */
    }

  parallel_sort(FOF_GList, NgroupsExt, sizeof(struct fof_group_list), fof_compare_FOF_GList_LocCountTaskDiffMinID);


  for(i = 0, ngr = 0; i < NgroupsExt; i++)
    {
      if(FOF_GList[i].ExtCount == FOF_GList[i].MinIDTask)
        ngr++;

      FOF_GList[i].GrNr = ngr - 1;
    }

  MPI_Allgather(&ngr, 1, MPI_INT, Send_count, 1, MPI_INT, MPI_COMM_WORLD);


  /* count how many groups there are on earlier CPUs */
  long long ngr_sum;
  for(j = 0, ngr_sum = 0; j < ThisTask; j++)
    ngr_sum += Send_count[j];

  for(i = 0; i < NgroupsExt; i++)
    FOF_GList[i].GrNr += ngr_sum;

  sumup_large_ints(1, &ngr, &ngr_sum);
  if(ngr_sum != TotNgroups)
    {
      printf("ngr_sum=%d\n", (int) ngr_sum);
      terminate("inconsistency");
    }

  /* bring the group list back into the original order */
  parallel_sort(FOF_GList, NgroupsExt, sizeof(struct fof_group_list), fof_compare_FOF_GList_ExtCountMinID);

  /* Assign the group numbers to the group properties array */
  for(i = 0, start = 0; i < Ngroups; i++)
    {
      while(FOF_GList[start].MinID < Group[i].MinID)
        {
          start++;
          if(start >= NgroupsExt)
            terminate("start >= NgroupsExt");
        }
      Group[i].GrNr = FOF_GList[start].GrNr;
    }


  /* sort the groups according to group-number */
  parallel_sort(Group, Ngroups, sizeof(struct group_properties), fof_compare_Group_GrNr);

  for(i = 0; i < NumPart; i++)
    PS[i].GrNr = TotNgroups + 1;        /* this marks all particles that are not in any group */

  for(i = 0, start = 0, Nids = 0; i < NgroupsExt; i++)
    {
      while(FOF_PList[start].MinID < FOF_GList[i].MinID)
        {
          start++;
          if(start > NumPart)
            terminate("start > NumPart");
        }

      if(FOF_PList[start].MinID != FOF_GList[i].MinID)
        terminate("FOF_PList[start=%d].MinID=%lld != FOF_GList[i=%d].MinID=%lld", start, (long long) FOF_PList[start].MinID, i, (long long) FOF_GList[i].MinID);

      for(lenloc = 0; start + lenloc < NumPart;)
        if(FOF_PList[start + lenloc].MinID == FOF_GList[i].MinID)
          {
            PS[FOF_PList[start + lenloc].Pindex].GrNr = FOF_GList[i].GrNr;
            Nids++;
            lenloc++;
          }
        else
          break;

      start += lenloc;
    }

  sumup_large_ints(1, &Nids, &totNids);

  if(totNids != TotNids)
    {
      char buf[1000];
      sprintf(buf, "Task=%d Nids=%d totNids=%d TotNids=%d\n", ThisTask, Nids, (int) totNids, (int) TotNids);
      terminate(buf);
    }

  t1 = second();

  mpi_printf("FOF: Assigning of group numbers took = %g sec\n", timediff(t0, t1));
}


//---------------------------------
void fof_compute_group_properties(int gr, int start, int len)
{
  int j, k, index, type, start_index = FOF_PList[start].Pindex;
  double xyz[3];

  // declare + initialize
//#ifdef GFM_STELLAR_EVOLUTION
//  MyFloat *GroupMassMetallicity = NULL, *GroupMassMetals = NULL;
//#endif

  Group[gr].Len = 0;
  double gr_Mass = 0;
#ifdef USE_SFR
  double gr_Sfr = 0;
#endif
#ifdef GFM_STELLAR_EVOLUTION
  double *GroupMassMetallicity = NULL, *GroupMassMetals = NULL; // FIXME
  double gr_GasMassMetallicity = 0, gr_StellarMassMetallicity = 0;
  double gr_GasMassMetals[GFM_N_CHEM_ELEMENTS], gr_StellarMassMetals[GFM_N_CHEM_ELEMENTS];
  for(k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
    {
      gr_GasMassMetals[k] = 0;
      gr_StellarMassMetals[k] = 0;
    }
#ifdef GFM_DUST
  double *GroupMassDustMetallicity = NULL;
  double gr_GasMassDustMetallicity = 0;
  int chan = 0;
#endif
#endif
#ifdef BLACK_HOLES
  double gr_BH_Mass = 0, gr_BH_Mdot = 0, gr_MaxDens = 0;
  Group[gr].index_maxdens = Group[gr].task_maxdens = -1;
#ifdef BH_NF_RADIO
  double gr_Min_BH_Potential = MAX_FLOAT_NUMBER;
  Group[gr].ID_Min_BH_Potential = 0;
  double gr_XrayLum = 0, gr_RadioLum = 0;
#endif
#endif
#ifdef GFM_WINDS
  double gr_WindMass = 0;
#endif
  double gr_CM[3], gr_Vel[3];
#ifdef GFM_BIPOLAR_WINDS
#if (GFM_BIPOLAR_WINDS != 3)
  double gr_GravAcc[3];
#endif
#endif
  for(k = 0; k < 3; k++)
    {
      gr_CM[k] = 0;
      gr_Vel[k] = 0;
#ifdef GFM_BIPOLAR_WINDS
#if (GFM_BIPOLAR_WINDS != 3)
      gr_GravAcc[k] = 0;
#endif
#endif
      Group[gr].FirstPos[k] = P[start_index].Pos[k];
    }

#if (GFM_BIPOLAR_WINDS == 3)
  double gr_MinPotential = MAX_FLOAT_NUMBER;
#endif

  double gr_MassType[NTYPES];
  for(k = 0; k < NTYPES; k++)
    {
      Group[gr].LenType[k] = 0;
      gr_MassType[k] = 0;
    }


  // calculate
  for(k = 0; k < len; k++)
    {
      index = FOF_PList[start + k].Pindex;

      Group[gr].Len++;
      gr_Mass += P[index].Mass;
      type = P[index].Type;

#ifdef TRACER_MC
      Group[gr].LenType[TRACER_MC] += get_number_of_tracers(index);
#endif

#ifdef GFM_WINDS_SAVE_PARTTYPE
      /* count wind as new particle type for LenType since we save these separately */
      if(P[index].Type == 4 && STP(index).BirthTime <= 0)
        Group[gr].LenType[GFM_WINDS_SAVE_PARTTYPE]++;
      else
#endif
      Group[gr].LenType[type]++;

#ifdef GFM_WINDS
      /* count wind as gas for mass, but not for LenType since we use this to construct offset tables */
      if(P[index].Type == 4 && STP(index).BirthTime <= 0)
        type = 0;
#endif
      gr_MassType[type] += P[index].Mass;

#ifdef BH_NF_RADIO
      if(P[index].Type == 0)
        gr_XrayLum += get_cooling_luminosity(index);
#endif

#ifdef USE_SFR
      if(P[index].Type == 0)
        gr_Sfr += SphP[index].Sfr;
#endif
#ifdef GFM_STELLAR_EVOLUTION
      if(P[index].Type == 0)
        {
          GroupMassMetallicity = &(gr_GasMassMetallicity);
          GroupMassMetals = gr_GasMassMetals;
        }
      if(P[index].Type == 4)
        {
          GroupMassMetallicity = &(gr_StellarMassMetallicity);
          GroupMassMetals = gr_StellarMassMetals;
        }
#ifdef GFM_WINDS
      if(P[index].Type == 4 && STP(index).BirthTime <= 0)
        {
          GroupMassMetallicity = &(gr_GasMassMetallicity);
          GroupMassMetals = gr_GasMassMetals;
        }
#endif
      if(P[index].Type == 0)
        {
#if (GFM_STELLAR_EVOLUTION==1)
          *GroupMassMetallicity += SphP[index].Metallicity * P[index].Mass;
#else
          *GroupMassMetallicity += SphP[index].MassMetallicity;
#endif
          for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
            {
#if (GFM_STELLAR_EVOLUTION==1)
              GroupMassMetals[j] += SphP[index].MetalsFraction[j] * P[index].Mass;
#else
              GroupMassMetals[j] += SphP[index].MassMetals[j];
#endif
            }
        }
      if(P[index].Type == 4)
        {
          *GroupMassMetallicity += STP(index).Metallicity * P[index].Mass;
          for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
            {
              GroupMassMetals[j] += STP(index).MassMetals[j];
            }
        }
#endif
#if (defined(GFM_STELLAR_EVOLUTION)) && (defined(GFM_DUST))
      if(P[index].Type == 0)
        {
          GroupMassDustMetallicity = &(gr_GasMassDustMetallicity);
        }
#ifdef GFM_WINDS
      if(P[index].Type == 4 && STP(index).BirthTime <= 0)
        {
          GroupMassDustMetallicity = &(gr_GasMassDustMetallicity);
        }
#endif
      if(P[index].Type == 0)
        {
          for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
            for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
              *GroupMassDustMetallicity += SphP[index].MetalsDustFraction[chan][j] * P[index].Mass;
        }
      if(P[index].Type == 4 && STP(index).BirthTime <= 0)
        {
          for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
            for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
              *GroupMassDustMetallicity += STP(index).MassMetals[j] * STP(index).InitialDustFractions[chan][j];
        }
#endif

#ifdef BLACK_HOLES
      if(P[index].Type == 5)
        {
          gr_BH_Mdot += BPP(index).BH_Mdot;
          gr_BH_Mass += BPP(index).BH_Mass;
#ifdef BH_NF_RADIO
          if(P[index].Potential < gr_Min_BH_Potential)
            {
              gr_Min_BH_Potential = P[index].Potential;
              Group[gr].ID_Min_BH_Potential = P[index].ID;
            }
#endif
        }
      if(P[index].Type == 0)
        {
          if(SphP[index].Density > gr_MaxDens)
            {
              gr_MaxDens = SphP[index].Density;
              Group[gr].index_maxdens = index;
              Group[gr].task_maxdens = ThisTask;
            }
        }
#endif

#ifdef GFM_WINDS
      if(P[index].Type == 4 && STP(index).BirthTime <= 0)
        gr_WindMass += P[index].Mass;
#endif

#if (GFM_BIPOLAR_WINDS == 3)
      if(P[index].Potential < gr_MinPotential)
        {
          gr_MinPotential = P[index].Potential;
          Group[gr].Pos_MinPotential[0] = P[index].Pos[0];
          Group[gr].Pos_MinPotential[1] = P[index].Pos[1];
          Group[gr].Pos_MinPotential[2] = P[index].Pos[2];
        }
#endif

      for(j = 0; j < 3; j++)
        {
          xyz[j] = P[index].Pos[j];
#ifdef PERIODIC
          xyz[j] = fof_periodic(xyz[j] - P[start_index].Pos[j]);
#endif
          gr_CM[j] += P[index].Mass * xyz[j];
          gr_Vel[j] += P[index].Mass * P[index].Vel[j];
#ifdef GFM_BIPOLAR_WINDS
#if  (GFM_BIPOLAR_WINDS != 3)
          gr_GravAcc[j] += P[index].Mass * P[index].GravAccel[j];
#endif
#endif
        }
    }

  // put values into group struct
  Group[gr].Mass = gr_Mass;
#ifdef USE_SFR
  Group[gr].Sfr = gr_Sfr;
#endif

#ifdef GFM_STELLAR_EVOLUTION
  Group[gr].GasMassMetallicity = gr_GasMassMetallicity;
  Group[gr].StellarMassMetallicity = gr_StellarMassMetallicity;
  for(k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
    {
      Group[gr].GasMassMetals[k] = gr_GasMassMetals[k];
      Group[gr].StellarMassMetals[k] = gr_StellarMassMetals[k];
    }
#ifdef GFM_DUST
  Group[gr].GasMassDustMetallicity = gr_GasMassDustMetallicity;
#endif
#endif

#ifdef BLACK_HOLES
  Group[gr].BH_Mass = gr_BH_Mass;
  Group[gr].BH_Mdot = gr_BH_Mdot;
  Group[gr].MaxDens = gr_MaxDens;
#ifdef BH_NF_RADIO
  Group[gr].Min_BH_Potential = gr_Min_BH_Potential;
  Group[gr].XrayLum = gr_XrayLum;
  Group[gr].RadioLum = gr_RadioLum;
#endif
#endif

#ifdef GFM_WINDS
  Group[gr].WindMass = gr_WindMass;
#endif
  for(k = 0; k < 3; k++)
    {
      Group[gr].CM[k] = gr_CM[k];
      Group[gr].Vel[k] = gr_Vel[k];
#ifdef GFM_BIPOLAR_WINDS
#if (GFM_BIPOLAR_WINDS != 3)
      Group[gr].GravAcc[k] = gr_GravAcc[k];
#endif
#endif
    }

#if (GFM_BIPOLAR_WINDS == 3)
  Group[gr].MinPotential = gr_MinPotential;
#endif

  for(k = 0; k < NTYPES; k++)
    Group[gr].MassType[k] = gr_MassType[k];


}

//---------------------------------
#if defined(BH_NF_RADIO)
double get_cooling_luminosity(int i)
{
  double lum = 0;
  double eos_dens_threshold = All.PhysDensThresh;
#ifdef MODIFIED_EOS
  eos_dens_threshold *= All.FactorDensThresh;
#endif

  double rho = SphP[i].Density * All.cf_a3inv;

  //if(SphP[i].Sfr == 0)
  if(rho < eos_dens_threshold)
    {
      double ne = SphP[i].Ne;

#if defined(GFM_AGN_RADIATION) || defined(GFM_UVB_CORRECTIONS)
#ifdef GFM_AGN_RADIATION
      update_radiation_state(rho, SphP[i].MetalsFraction[element_index_Hydrogen], SphP[i].AGNBolIntensity);
#else
      update_radiation_state(rho, SphP[i].MetalsFraction[element_index_Hydrogen], 0);
#endif
#endif

#ifdef GFM_COOLING_METAL
      update_gas_state(rho, SphP[i].MetalsFraction[element_index_Hydrogen], SphP[i].Metallicity);
#endif

      double tcool = GetCoolingTime(SphP[i].Utherm, rho, &ne);

      if(tcool > 0)
        lum = P[i].Mass * SphP[i].Utherm / tcool;
    }

  return lum;
}
#endif


void fof_exchange_group_data(void)
{
  struct group_properties *get_Group;
  int i, j, ngrp, recvTask, nimport, start;
  double xyz[3];

  /* sort the groups according to task */
  mysort(Group, NgroupsExt, sizeof(struct group_properties), fof_compare_Group_MinIDTask);

  /* count how many we have of each task */
  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;
  for(i = 0; i < NgroupsExt; i++)
    Send_count[FOF_GList[i].MinIDTask]++;

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      if(j == ThisTask)         /* we will not exchange the ones that are local */
        Recv_count[j] = 0;
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  get_Group = (struct group_properties *) mymalloc("get_Group", sizeof(struct group_properties) * nimport);

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the group data */
              MPI_Sendrecv(&Group[Send_offset[recvTask]],
                           Send_count[recvTask] * sizeof(struct group_properties), MPI_BYTE,
                           recvTask, TAG_DENS_A,
                           &get_Group[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct group_properties), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  /* sort the groups again according to MinID */
  mysort(Group, NgroupsExt, sizeof(struct group_properties), fof_compare_Group_MinID);
  mysort(get_Group, nimport, sizeof(struct group_properties), fof_compare_Group_MinID);

  /* now add in the partial imported group data to the main ones */
  for(i = 0, start = 0; i < nimport; i++)
    {
      while(Group[start].MinID < get_Group[i].MinID)
        {
          start++;
          if(start >= NgroupsExt)
            terminate("start >= NgroupsExt");
        }

      Group[start].Len += get_Group[i].Len;
      Group[start].Mass += get_Group[i].Mass;

      for(j = 0; j < NTYPES; j++)
        {
          Group[start].LenType[j] += get_Group[i].LenType[j];
          Group[start].MassType[j] += get_Group[i].MassType[j];
        }

#ifdef USE_SFR
      Group[start].Sfr += get_Group[i].Sfr;
#endif
#ifdef BH_NF_RADIO
      Group[start].XrayLum += get_Group[i].XrayLum;
#endif

#ifdef GFM_STELLAR_EVOLUTION
      Group[start].GasMassMetallicity += get_Group[i].GasMassMetallicity;
      Group[start].StellarMassMetallicity += get_Group[i].StellarMassMetallicity;
      for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
        {
          Group[start].GasMassMetals[j] += get_Group[i].GasMassMetals[j];
          Group[start].StellarMassMetals[j] += get_Group[i].StellarMassMetals[j];
        }
#ifdef GFM_DUST
      Group[start].GasMassDustMetallicity += get_Group[i].GasMassDustMetallicity;
#endif
#endif
#ifdef BLACK_HOLES
      Group[start].BH_Mdot += get_Group[i].BH_Mdot;
      Group[start].BH_Mass += get_Group[i].BH_Mass;
      if(get_Group[i].MaxDens > Group[start].MaxDens)
        {
          Group[start].MaxDens = get_Group[i].MaxDens;
          Group[start].index_maxdens = get_Group[i].index_maxdens;
          Group[start].task_maxdens = get_Group[i].task_maxdens;
        }

#if defined(BH_NF_RADIO)
      if(get_Group[i].Min_BH_Potential < Group[start].Min_BH_Potential)
        {
          Group[start].Min_BH_Potential = get_Group[i].Min_BH_Potential;
          Group[start].ID_Min_BH_Potential = get_Group[i].ID_Min_BH_Potential;
        }
#endif

#endif

#ifdef GFM_WINDS
      Group[start].WindMass += get_Group[i].WindMass;
#endif

#if (GFM_BIPOLAR_WINDS == 3)
      if(get_Group[i].MinPotential < Group[start].MinPotential)
        {
          Group[start].MinPotential = get_Group[i].MinPotential;
          Group[start].Pos_MinPotential[0] = get_Group[i].Pos_MinPotential[0];
          Group[start].Pos_MinPotential[1] = get_Group[i].Pos_MinPotential[1];
          Group[start].Pos_MinPotential[2] = get_Group[i].Pos_MinPotential[2];
        }
#endif

      for(j = 0; j < 3; j++)
        {
          xyz[j] = get_Group[i].CM[j] / get_Group[i].Mass;
#ifdef PERIODIC
          xyz[j] = fof_periodic(xyz[j] + get_Group[i].FirstPos[j] - Group[start].FirstPos[j]);
#endif
          Group[start].CM[j] += get_Group[i].Mass * xyz[j];
          Group[start].Vel[j] += get_Group[i].Vel[j];
#ifdef GFM_BIPOLAR_WINDS
#if (GFM_BIPOLAR_WINDS != 3)
          Group[start].GravAcc[j] += get_Group[i].GravAcc[j];
#endif
#endif
        }
    }

  myfree(get_Group);
}

void fof_finish_group_properties(void)
{
  double cm[3];
  int i, j, ngr;

  for(i = 0; i < NgroupsExt; i++)
    {
      if(Group[i].MinIDTask == ThisTask)
        {
          for(j = 0; j < 3; j++)
            {
              Group[i].Vel[j] /= Group[i].Mass;
#ifdef GFM_BIPOLAR_WINDS
#if  (GFM_BIPOLAR_WINDS != 3)
              Group[i].GravAcc[j] /= Group[i].Mass;
#endif
#endif
              cm[j] = Group[i].CM[j] / Group[i].Mass;
#ifdef PERIODIC
              cm[j] = fof_periodic_wrap(cm[j] + Group[i].FirstPos[j]);
#endif
              Group[i].CM[j] = cm[j];
            }

#ifdef BH_NF_RADIO
          if(Group[i].LenType[5] > 0)   /* only if a BH is actually present */
            {
              double vvir = pow(10.0 * All.G * All.cf_H * Group[i].Mass, 1.0 / 3);
              Group[i].RadioLum = (All.Hubble / All.cf_H) * blackhole_get_radio_efficiency(vvir) * Group[i].XrayLum;
            }
          else
            Group[i].RadioLum = 0;
#endif

        }
    }

  /* eliminate the non-local groups */
  for(i = 0, ngr = NgroupsExt; i < ngr; i++)
    {
      if(Group[i].MinIDTask != ThisTask)
        {
          Group[i] = Group[ngr - 1];
          i--;
          ngr--;
        }
    }

  if(ngr != Ngroups)
    terminate("ngr != Ngroups");

  mysort(Group, Ngroups, sizeof(struct group_properties), fof_compare_Group_MinID);
}



double fof_periodic(double x)
{
  if(x >= 0.5 * All.BoxSize)
    x -= All.BoxSize;
  if(x < -0.5 * All.BoxSize)
    x += All.BoxSize;
  return x;
}


double fof_periodic_wrap(double x)
{
  while(x >= All.BoxSize)
    x -= All.BoxSize;
  while(x < 0)
    x += All.BoxSize;
  return x;
}


#endif /* of FOF */
