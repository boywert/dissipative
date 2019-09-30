/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/subfind/subfind_fof_reprocess.c
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
#include "../fof/fof.h"
#include "subfind.h"

#ifdef HAVE_HDF5
#include <hdf5.h>


#ifdef SUBFIND_EXTENDED_PROPERTIES
void subfind_add_grp_props_calc_fof_angular_momentum(int num, int ngroups_cat)
{
  mpi_printf("FOF: Begin Angular Momentum Calculation for FOF Groups.\n");

  /* assign target CPUs for the particles in groups */
  /* the particles not in groups will be distributed such that a uniform particle load results */
  double t0 = second();
  int *count_loc_task = mymalloc_clear("count_loc_task", NTask * sizeof(int));
  int *count_task = mymalloc("count_task", NTask * sizeof(int));
  int *count_free = mymalloc("count_free", NTask * sizeof(int));
  int count_loc_free = 0;

  for(int i = 0; i < NumPart; i++)
    {
      if(PS[i].GrNr < 0)
        terminate("PS[i].GrNr=%d", PS[i].GrNr);

      if(PS[i].GrNr < TotNgroups)       /* particle is in a group */
        {
          if(PS[i].GrNr < Ncollective)  /* we are in a collective group */
            PS[i].TargetTask = ProcAssign[PS[i].GrNr].FirstTask + (i % ProcAssign[PS[i].GrNr].NTask);
          else
            PS[i].TargetTask = ((PS[i].GrNr - Ncollective) % (NTask - NprocsCollective)) + NprocsCollective;

          if(PS[i].TargetTask < 0 || PS[i].TargetTask >= NTask)
            terminate("PS[i].TargetTask=%d PS[i].GrNr=%d", PS[i].TargetTask, PS[i].GrNr);

          count_loc_task[PS[i].TargetTask]++;
        }
      else
        count_loc_free++;

      PS[i].TargetIndex = 0;    /* unimportant here */
    }

  MPI_Allgather(&count_loc_free, 1, MPI_INT, count_free, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Allreduce(count_loc_task, count_task, NTask, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  long long sum = 0;
  for(int i = 0; i < NTask; i++)
    sum += count_task[i] + count_free[i];

  int maxload = (sum + NTask - 1) / NTask;
  for(int i = 0; i < NTask; i++)
    {
      count_task[i] = maxload - count_task[i];  /* this is the amount that can fit on this task */
      if(count_task[i] < 0)
        count_task[i] = 0;
    }

  int current_task = 0;

  for(int i = 0; i < ThisTask; i++)
    {
      while(count_free[i] > 0 && current_task < NTask)
        {
          if(count_free[i] < count_task[current_task])
            {
              count_task[current_task] -= count_free[i];
              count_free[i] = 0;
            }
          else
            {
              count_free[i] -= count_task[current_task];
              count_task[current_task] = 0;
              current_task++;
            }
        }
    }

  for(int i = 0; i < NumPart; i++)
    {
      if(PS[i].GrNr >= TotNgroups)      /* particle not in a group. Can in principle stay but we move it such that a good load balance is obtained. */
        {
          while(count_task[current_task] == 0 && current_task < NTask - 1)
            current_task++;

          PS[i].TargetTask = current_task;      /* particle not in any group, move it here so that uniform load is achieved */
          count_task[current_task]--;
        }
    }

  myfree(count_free);
  myfree(count_task);
  myfree(count_loc_task);

  double balance = subfind_get_particle_balance();
  mpi_printf("SUBFIND: particle balance=%g\n", balance);

  /* distribute particles such that groups are completely on the CPU(s) that do the corresponding group(s) */
  fof_subfind_exchange(MPI_COMM_WORLD);
  double t1 = second();
  mpi_printf("SUBFIND: subfind_exchange() took %g sec\n", timediff(t0, t1));

  balance = subfind_get_particle_balance();
  mpi_printf("SUBFIND: particle balance for AM processing=%g\n", balance);

  /* we can now split the communicator to give each collectively treated group its own processor set */
  MPI_Comm_split(MPI_COMM_WORLD, CommSplitColor, ThisTask, &SubComm);
  MPI_Comm_size(SubComm, &SubNTask);
  MPI_Comm_rank(SubComm, &SubThisTask);
  SubTagOffset = TagOffset;

  /* here the execution paths for collective groups and serial groups branch. The collective CPUs work in small sets that each
   * deal with one large group. The serial CPUs each deal with several halos by themselves
   */
  if(CommSplitColor < Ncollective)      /* we are one of the CPUs that does a collective group */
    {
      /* we now apply a collective version of subfind to the group split across the processors belonging to communicator SubComm
       * The relevant group is the one stored in Group[0] on SubThisTask==0.
       */
      subfind_fof_calc_am_collective(num, ngroups_cat);
    }
  else
    {
      /* now let us sort according to GrNr and Density. This step will temporarily break the association with SphP[] and other arrays! */
      submp = (struct submp_data *) mymalloc("submp", sizeof(struct submp_data) * NumPart);
      for(int i = 0; i < NumPart; i++)
        {
          PS[i].OldIndex = i;
          submp[i].index = i;
          submp[i].GrNr = PS[i].GrNr;
          submp[i].DM_Density = PS[i].Density;
        }
      qsort(submp, NumPart, sizeof(struct submp_data), subfind_compare_submp_GrNr_DM_Density);
      subfind_reorder_according_to_submp();
      myfree(submp);

      /* now we have the particles in each group consecutively */
      if(SubThisTask == 0)
        printf("SUBFIND-SERIAL: Start to do AM for %d small groups with serial subfind algorithm on %d processors (root-node=%d)\n", TotNgroups - Ncollective, SubNTask, ThisTask);

      /* we now apply a serial version of subfind to the local groups */

      t0 = second();
      for(int gr = 0, offset = 0; gr < Ngroups; gr++)
        {
          if(((Group[gr].GrNr - Ncollective) % (NTask - NprocsCollective)) + NprocsCollective == ThisTask)
            offset = subfind_fof_calc_am_serial(gr, offset, num, ngroups_cat);
          else
            terminate("how come that we have this group number?");
        }

      MPI_Barrier(SubComm);
      t1 = second();
      if(SubThisTask == 0)
        printf("SUBFIND-SERIAL: processing AM of serial groups took %g sec\n", timediff(t0, t1));

      /* undo local rearrangement that made groups consecutive. After that, the association of SphP[] will be correct again */
      submp = (struct submp_data *) mymalloc("submp", sizeof(struct submp_data) * NumPart);
      for(int i = 0; i < NumPart; i++)
        {
          submp[i].index = i;
          submp[i].OldIndex = PS[i].OldIndex;
        }
      qsort(submp, NumPart, sizeof(struct submp_data), subfind_compare_submp_OldIndex);
      subfind_reorder_according_to_submp();
      myfree(submp);
    }

  /* free the communicator */
  MPI_Comm_free(&SubComm);

  /* distribute particles back to original CPU */
  t0 = second();
  for(int i = 0; i < NumPart; i++)
    {
      PS[i].TargetTask = PS[i].OriginTask;
      PS[i].TargetIndex = PS[i].OriginIndex;
    }

  fof_subfind_exchange(MPI_COMM_WORLD);
  t1 = second();
  if(ThisTask == 0)
    printf("SUBFIND: subfind_exchange() (for return to original CPU after AM)  took %g sec\n", timediff(t0, t1));

  mpi_printf("FOF: Angular Momentum Calculation for FOF Groups finished successfully.\n");

#ifdef ADD_GROUP_PROPERTIES
  myfree_movable(subgroupidlist);
  myfree_movable(groupidlist);
  myfree_movable(SubGroupMostBoundID);
  myfree_movable(GroupCat);
#endif
}
#endif


#ifdef ADD_GROUP_PROPERTIES     // used with RestartFlag = 3 only

static int *SubGroupGrNr;
static int *SubGroupLenType;


void subfind_add_grp_props_read_catalogue(int num, int ngroups_cat, int nsubgroups_cat)
{
  int i;

  // allocate structs for group properties from group catalogue
  GroupCat = (struct group_catalogue *) mymalloc_movable(&GroupCat, "GroupCat", ngroups_cat * sizeof(struct group_catalogue));
  SubGroupMostBoundID = (MyIDType *) mymalloc_movable(&SubGroupMostBoundID, "SubGroupMostBoundID", nsubgroups_cat * sizeof(MyIDType));

  SubGroupPos = (MyFloat *) mymalloc_movable(&SubGroupPos, "SubGroupPos", nsubgroups_cat * 3 * sizeof(MyFloat));

  // read in the group properties from group catalogue (especially identifier IDs) into above structs
  get_catalogue_prop(num);

  assign_group_numbers_based_on_catalogue(ngroups_cat, nsubgroups_cat);


  /* make list so that GroupCat can be accessed in increasing order of FirstSubMostBoundID */
  groupidlist = (struct idsortlist *) mymalloc_movable(&groupidlist, "grouplist", ngroups_cat * sizeof(struct idsortlist));
  subgroupidlist = (struct idsortlist *) mymalloc_movable(&subgroupidlist, "subgrouplist", nsubgroups_cat * sizeof(struct idsortlist));

  for(i = 0; i < ngroups_cat; i++)      // do this only for groups with subs!!
    {
      groupidlist[i].MostBoundID = GroupCat[i].FirstSubMostBoundID;
      groupidlist[i].index = i;
    }

  for(i = 0; i < nsubgroups_cat; i++)
    {
      subgroupidlist[i].MostBoundID = SubGroupMostBoundID[i];
      subgroupidlist[i].index = i;
    }

  /* sort according to MostBoundID */
  qsort(groupidlist, ngroups_cat, sizeof(struct idsortlist), subfind_fof_compare_MostBoundID);
  qsort(subgroupidlist, nsubgroups_cat, sizeof(struct idsortlist), subfind_fof_compare_MostBoundID);
}



void subfind_add_grp_props_finalize(int num, int ngroups_cat, int nsubgroups_cat)
{
  // append additional group properties
  double t0 = second();

  if(ThisTask == 0)
    GroupAll = (struct group_properties *) mymalloc("GroupAll", TotNgroups * sizeof(struct group_properties));

  mpi_printf("FOF REPROCESS: Task 0 has to collect and store locally %d groups corresponding to %lld Byte.\n", TotNgroups, TotNgroups * sizeof(struct group_properties));

  fof_collect_groups();

  if(ThisTask == 0)
    {
      fof_append_group_properties(num, ngroups_cat);
      myfree(GroupAll);
    }

  double t1 = second();
  mpi_printf("FOF/SUBFIND: Group catalogue extension (FOF) saved. took = %g sec\n", timediff(t0, t1));
  mpi_printf("FOF REPROCESS: done.\n");

  // append additional subhalo properties
  t0 = second();

  if(ThisTask == 0)
    SubGroupAll = (struct subgroup_properties *) mymalloc("SubGroupAll", TotNsubgroups * sizeof(struct subgroup_properties));

  mpi_printf("SUBFIND REPROCESS: Task 0 has to collect and store locally %d subgroups corresponding to %lld Byte.\n", TotNsubgroups, TotNsubgroups * sizeof(struct subgroup_properties));

  fof_collect_subgroups();

  if(ThisTask == 0)
    {
      fof_append_subgroup_properties(num, nsubgroups_cat);
      myfree(SubGroupAll);
    }

  t1 = second();
  mpi_printf("FOF/SUBFIND: Group catalogue extension (Subhaloes) saved. took = %g sec\n", timediff(t0, t1));
  mpi_printf("SUBFIND REPROCESS: done.\n");

  endrun();
}



int get_number_of_group_catalogue_files(int snapnr)
{
  int numfiles, flag_doubleprecision;
  hid_t hdf5_file, hdf5_headergrp, hdf5_attribute;

  char fname[1000];

  if(All.NumFilesPerSnapshot > 1)
    sprintf(fname, "%sgroups_%03d/%s_%03d.%d.hdf5", All.OutputDir, snapnr, "fof_subhalo_tab", snapnr, 0);
  else
    sprintf(fname, "%s%s_%03d.hdf5", All.OutputDir, "fof_subhalo_tab", snapnr);

  hdf5_file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if(hdf5_file < 0)
    terminate("cannot read initial conditions file %s", fname);
  hdf5_headergrp = H5Gopen(hdf5_file, "/Header");

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "FlagDoubleprecision");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &flag_doubleprecision);       /* read statement */
  H5Aclose(hdf5_attribute);

#ifdef INPUT_IN_DOUBLEPRECISION
  if(flag_doubleprecision == 0)
    terminate("\nProblem: Code compiled with INPUT_IN_DOUBLEPRECISION, but input files are in single precision!\n");
#else
  if(flag_doubleprecision)
    terminate("\nProblem: Code not compiled with INPUT_IN_DOUBLEPRECISION, but input files are in double precision!\n");
#endif

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumFiles");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &numfiles);   /* read statement */
  H5Aclose(hdf5_attribute);

  H5Gclose(hdf5_headergrp);
  H5Fclose(hdf5_file);

  return numfiles;
}

int get_number_of_groups_in_file(int snapnr, int fnr)
{
  int numgroups, flag_doubleprecision;
  hid_t hdf5_file, hdf5_headergrp, hdf5_attribute;

  char fname[1000];

  if(All.NumFilesPerSnapshot > 1)
    sprintf(fname, "%sgroups_%03d/%s_%03d.%d.hdf5", All.OutputDir, snapnr, "fof_subhalo_tab", snapnr, fnr);
  else
    sprintf(fname, "%s%s_%03d.hdf5", All.OutputDir, "fof_subhalo_tab", snapnr);

  hdf5_file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if(hdf5_file < 0)
    terminate("cannot read initial conditions file %s", fname);
  hdf5_headergrp = H5Gopen(hdf5_file, "/Header");

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "FlagDoubleprecision");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &flag_doubleprecision);       /* read statement */
  H5Aclose(hdf5_attribute);

#ifdef INPUT_IN_DOUBLEPRECISION
  if(flag_doubleprecision == 0)
    terminate("\nProblem: Code compiled with INPUT_IN_DOUBLEPRECISION, but input files are in single precision!\n");
#else
  if(flag_doubleprecision)
    terminate("\nProblem: Code not compiled with INPUT_IN_DOUBLEPRECISION, but input files are in double precision!\n");
#endif

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Ngroups_ThisFile");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &numgroups);  /* read statement */
  H5Aclose(hdf5_attribute);

  H5Gclose(hdf5_headergrp);
  H5Fclose(hdf5_file);

  return numgroups;
}

int get_number_of_subgroups_in_file(int snapnr, int fnr)        // only for subfind_fof_reprocess.c
{
  int nsubgroups, flag_doubleprecision;
  hid_t hdf5_file, hdf5_headergrp, hdf5_attribute;

  char fname[1000];

  if(All.NumFilesPerSnapshot > 1)
    sprintf(fname, "%sgroups_%03d/%s_%03d.%d.hdf5", All.OutputDir, snapnr, "fof_subhalo_tab", snapnr, fnr);
  else
    sprintf(fname, "%s%s_%03d.hdf5", All.OutputDir, "fof_subhalo_tab", snapnr);

  hdf5_file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if(hdf5_file < 0)
    terminate("cannot read initial conditions file %s", fname);
  hdf5_headergrp = H5Gopen(hdf5_file, "/Header");

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "FlagDoubleprecision");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &flag_doubleprecision);       /* read statement */
  H5Aclose(hdf5_attribute);

#ifdef INPUT_IN_DOUBLEPRECISION
  if(flag_doubleprecision == 0)
    terminate("\nProblem: Code compiled with INPUT_IN_DOUBLEPRECISION, but input files are in single precision!\n");
#else
  if(flag_doubleprecision)
    terminate("\nProblem: Code not compiled with INPUT_IN_DOUBLEPRECISION, but input files are in double precision!\n");
#endif

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Nsubgroups_ThisFile");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &nsubgroups); /* read statement */
  H5Aclose(hdf5_attribute);

  H5Gclose(hdf5_headergrp);
  H5Fclose(hdf5_file);

  return nsubgroups;
}

int get_number_of_groups_in_catalogue(int snapnr)
{
  int fnr, number_of_files, ngroups_total = 0;
  number_of_files = get_number_of_group_catalogue_files(snapnr);        // fct def in subfind_so_reprocess.c 

  for(fnr = 0; fnr < number_of_files; fnr++)
    ngroups_total += get_number_of_groups_in_file(snapnr, fnr); // fct def in subfind_so_reprocess.c

  return ngroups_total;
}

int get_number_of_subgroups_in_catalogue(int snapnr)
{
  int fnr, number_of_files, nsubgroups_total = 0;
  number_of_files = get_number_of_group_catalogue_files(snapnr);        // fct def in subfind_so_reprocess.c 

  for(fnr = 0; fnr < number_of_files; fnr++)
    nsubgroups_total += get_number_of_subgroups_in_file(snapnr, fnr);   // fct def in subfind_so_reprocess.c

  return nsubgroups_total;
}

void get_catalogue_prop(int snapnr)
{
  int i;
  double t0, t1;
  t0 = second();

  // this is happening on every task 
  int fnr, number_of_files, ngroups_total = 0, nsubgroups_total = 0;
  number_of_files = get_number_of_group_catalogue_files(snapnr);        // fct def in subfind_so_reprocess.c
  int ngroups_in_file[number_of_files], nsubgroups_in_file[number_of_files];
  int group_range[number_of_files + 1], subgroup_range[number_of_files + 1];

  calculate_maxid();
  if(All.MaxID + 1 < All.MaxID)
    terminate("FOF REPROCESS: we will have an overflow in MaxID");      // then we have to think of something smarter

  for(fnr = 0; fnr < number_of_files; fnr++)
    {
      ngroups_in_file[fnr] = get_number_of_groups_in_file(snapnr, fnr); // fct def in subfind_so_reprocess.c
      ngroups_total += ngroups_in_file[fnr];
      nsubgroups_in_file[fnr] = get_number_of_subgroups_in_file(snapnr, fnr);   // fct def in subfind_so_reprocess.c
      nsubgroups_total += nsubgroups_in_file[fnr];
    }

  group_range[0] = 0;
  subgroup_range[0] = 0;
  for(fnr = 0; fnr < number_of_files; fnr++)
    {
      group_range[fnr + 1] = group_range[fnr] + ngroups_in_file[fnr];
      subgroup_range[fnr + 1] = subgroup_range[fnr] + nsubgroups_in_file[fnr];
    }

  // get group cat Nsubs and FirstSub into GroupCat
  fof_read_group_prop(snapnr, number_of_files, &group_range[0]);        // group_range instead of ngroups as we need it as boundaries for GroupCat 

  // get group cat SubMostBoundID of FirstSub into GroupcCat
  SubGroupGrNr = mymalloc_movable(&SubGroupGrNr, "SubGroupGrNr", nsubgroups_total * sizeof(int));
  memset(SubGroupGrNr, 0, nsubgroups_total * sizeof(int));

  SubGroupLenType = mymalloc_movable(&SubGroupLenType, "SubGroupLenType", 6 * nsubgroups_total * sizeof(int));

  fof_read_subgroup_prop(snapnr, number_of_files, &subgroup_range[0]);

  for(i = 0; i < ngroups_total; i++)
    {
      if(GroupCat[i].Nsubs > 0)
        {
          GroupCat[i].FirstSubMostBoundID = SubGroupMostBoundID[GroupCat[i].FirstSub];
          GroupCat[i].GrNr = SubGroupGrNr[GroupCat[i].FirstSub];
        }
      else
        {
          GroupCat[i].FirstSubMostBoundID = All.MaxID + 1;
          GroupCat[i].GrNr = All.MaxID + 1;
        }
    }

  t1 = second();
  mpi_printf("FOF/SUBFIND REPROCESS: Group catalogue properties imported. took = %g sec\n", timediff(t0, t1));
}


void fof_read_group_prop(int snapnr, int number_of_files, int *group_range)
{
  int i, j, rank, fnr, bnr;
  hid_t hdf5_file, hdf5_grp, hdf5_dataset, hdf5_dataspace_in_memory, hdf5_datatype, hdf5_dataspace_in_file;
  hsize_t dims[2];
  enum fof_subfind_iofields blocknr;
  char fname[1000];

  for(fnr = 0; fnr < number_of_files; fnr++)
    {
      if(All.NumFilesPerSnapshot > 1)
        sprintf(fname, "%sgroups_%03d/%s_%03d.%d.hdf5", All.OutputDir, snapnr, "fof_subhalo_tab", snapnr, fnr);
      else
        sprintf(fname, "%s%s_%03d.hdf5", All.OutputDir, "fof_subhalo_tab", snapnr);

      hdf5_file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
      if(hdf5_file < 0)
        terminate("cannot read initial conditions file %s", fname);
      hdf5_grp = H5Gopen(hdf5_file, "/Group");

      for(bnr = 0; bnr < 1000; bnr++)
        {
          blocknr = (enum fof_subfind_iofields) bnr;

          if((blocknr == IO_FOF_NSUBS) || (blocknr == IO_FOF_FIRSTSUB) || (blocknr == IO_FOF_LENTYPE))
            {
              char buf[1000];
              fof_subfind_get_dataset_name(blocknr, buf);

              hdf5_dataset = H5Dopen(hdf5_grp, buf);

              dims[0] = group_range[fnr + 1] - group_range[fnr];        // ngroups[fnr];
              dims[1] = fof_subfind_get_values_per_blockelement(blocknr);
              if(dims[1] == 1)
                rank = 1;
              else
                rank = 2;

              hdf5_datatype = H5Tcopy(H5T_NATIVE_INT);
              int *tmp = mymalloc_clear("tmp", dims[0] * dims[1] * sizeof(int));

              hdf5_dataspace_in_memory = H5Screate_simple(rank, dims, NULL);
              hdf5_dataspace_in_file = H5Screate_simple(rank, dims, NULL);

              H5Dread(hdf5_dataset, hdf5_datatype, hdf5_dataspace_in_memory, hdf5_dataspace_in_file, H5P_DEFAULT, tmp);

              for(i = group_range[fnr], j = 0; i < group_range[fnr + 1]; i++)
                {
                  switch (blocknr)
                    {
                    case IO_FOF_NSUBS:
                      GroupCat[i].Nsubs = tmp[j];
                      j++;
                      break;
                    case IO_FOF_FIRSTSUB:
                      GroupCat[i].FirstSub = tmp[j];
                      j++;
                      break;
                    case IO_FOF_LENTYPE:
                      for(int count = 0; count < 6; count++, j++)
                        {
                          GroupCat[i].GroupLenType[count] = tmp[j];
                        }
                      break;

                    default:
                      break;
                    }
                }

              myfree(tmp);

              H5Dclose(hdf5_dataset);
            }

          if((blocknr == IO_FOF_R_MEAN200) || (blocknr == IO_FOF_R_TOPHAT200) || (blocknr == IO_FOF_R_CRIT200)
              || (blocknr == IO_FOF_R_CRIT500))
            {
              char buf[1000];
              fof_subfind_get_dataset_name(blocknr, buf);

              hdf5_dataset = H5Dopen(hdf5_grp, buf);

              dims[0] = group_range[fnr + 1] - group_range[fnr];
              dims[1] = fof_subfind_get_values_per_blockelement(blocknr);
              if(dims[1] == 1)
                rank = 1;
              else
                rank = 2;

#ifdef INPUT_IN_DOUBLEPRECISION
              hdf5_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
              double *tmp = mymalloc_clear("tmp", dims[0] * dims[1] * sizeof(double));
#else
              hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
              float *tmp = mymalloc_clear("tmp", dims[0] * dims[1] * sizeof(float));
#endif

              hdf5_dataspace_in_memory = H5Screate_simple(rank, dims, NULL);
              hdf5_dataspace_in_file = H5Screate_simple(rank, dims, NULL);

              H5Dread(hdf5_dataset, hdf5_datatype, hdf5_dataspace_in_memory, hdf5_dataspace_in_file, H5P_DEFAULT, tmp);

              for(i = group_range[fnr], j = 0; i < group_range[fnr + 1]; i++)
                {
                  switch(blocknr)
                    {
                  case IO_FOF_R_MEAN200:
                    GroupCat[i].R_Mean200 = tmp[j++];
                    break;
                  case IO_FOF_R_TOPHAT200:
                    GroupCat[i].R_TopHat200 = tmp[j++];
                    break;
                  case IO_FOF_R_CRIT200:
                    GroupCat[i].R_Crit200 = tmp[j++];
                    break;
                  case IO_FOF_R_CRIT500:
                    GroupCat[i].R_Crit500 = tmp[j++];
                    break;
                  default:
                    break;
                    }
                }

              myfree(tmp);

              H5Dclose(hdf5_dataset);
            }


          if(blocknr == IO_FOF_LASTENTRY)
            break;
        }

      H5Gclose(hdf5_grp);
      H5Fclose(hdf5_file);

    }
}

void fof_read_subgroup_prop(int snapnr, int number_of_files, int *subgroup_range)
{
  int i, j, rank, fnr, bnr;
  hid_t hdf5_file, hdf5_grp, hdf5_dataset, hdf5_dataspace_in_memory, hdf5_datatype, hdf5_dataspace_in_file;
  hsize_t dims[2];
  enum fof_subfind_iofields blocknr;
  char fname[1000];

  for(fnr = 0; fnr < number_of_files; fnr++)
    {
      if(All.NumFilesPerSnapshot > 1)
        sprintf(fname, "%sgroups_%03d/%s_%03d.%d.hdf5", All.OutputDir, snapnr, "fof_subhalo_tab", snapnr, fnr);
      else
        sprintf(fname, "%s%s_%03d.hdf5", All.OutputDir, "fof_subhalo_tab", snapnr);

      hdf5_file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

      if(hdf5_file < 0)
        terminate("cannot read initial conditions file %s", fname);
      hdf5_grp = H5Gopen(hdf5_file, "/Subhalo");

      for(bnr = 0; bnr < 1000; bnr++)
        {
          blocknr = (enum fof_subfind_iofields) bnr;

          if(blocknr == IO_SUB_IDMOSTBOUND)
            {
              char buf[1000];
              fof_subfind_get_dataset_name(blocknr, buf);

              hdf5_dataset = H5Dopen(hdf5_grp, buf);

              dims[0] = subgroup_range[fnr + 1] - subgroup_range[fnr];
              dims[1] = fof_subfind_get_values_per_blockelement(blocknr);
              if(dims[1] == 1)
                rank = 1;
              else
                rank = 2;

#ifdef LONGIDS                  // changes needed for SPECIAL_BOUNDARY switch!!
              hdf5_datatype = H5Tcopy(H5T_NATIVE_ULLONG);       // Types stimmen - checked.
#else
              hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT);
#endif
              MyIDType *tmp = mymalloc_clear("tmp", dims[0] * dims[1] * sizeof(MyIDType));

              hdf5_dataspace_in_memory = H5Screate_simple(rank, dims, NULL);
              hdf5_dataspace_in_file = H5Screate_simple(rank, dims, NULL);

              H5Dread(hdf5_dataset, hdf5_datatype, hdf5_dataspace_in_memory, hdf5_dataspace_in_file, H5P_DEFAULT, tmp);

              for(i = subgroup_range[fnr], j = 0; i < subgroup_range[fnr + 1]; i++)
                {
                  SubGroupMostBoundID[i] = tmp[j];
                  j++;
                }

              myfree(tmp);

              H5Dclose(hdf5_dataset);
            }

          if(blocknr == IO_SUB_GRNR)
            {
              char buf[1000];
              fof_subfind_get_dataset_name(blocknr, buf);

              hdf5_dataset = H5Dopen(hdf5_grp, buf);

              dims[0] = subgroup_range[fnr + 1] - subgroup_range[fnr];
              dims[1] = fof_subfind_get_values_per_blockelement(blocknr);
              if(dims[1] == 1)
                rank = 1;
              else
                rank = 2;

              hdf5_datatype = H5Tcopy(H5T_NATIVE_INT);
              int *tmp = mymalloc_clear("tmp", dims[0] * dims[1] * sizeof(int));

              hdf5_dataspace_in_memory = H5Screate_simple(rank, dims, NULL);
              hdf5_dataspace_in_file = H5Screate_simple(rank, dims, NULL);

              H5Dread(hdf5_dataset, hdf5_datatype, hdf5_dataspace_in_memory, hdf5_dataspace_in_file, H5P_DEFAULT, tmp);

              for(i = subgroup_range[fnr], j = 0; i < subgroup_range[fnr + 1]; i++)
                {
                  SubGroupGrNr[i] = tmp[j];
                  j++;
                }

              myfree(tmp);

              H5Dclose(hdf5_dataset);
            }

          if(blocknr == IO_SUB_POS)
            {
              char buf[1000];
              fof_subfind_get_dataset_name(blocknr, buf);

              hdf5_dataset = H5Dopen(hdf5_grp, buf);

              dims[0] = subgroup_range[fnr + 1] - subgroup_range[fnr];
              dims[1] = fof_subfind_get_values_per_blockelement(blocknr);
              if(dims[1] == 1)
                rank = 1;
              else
                rank = 2;

              hdf5_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
              double *tmp = mymalloc_clear("tmp", dims[0] * dims[1] * sizeof(double));

              hdf5_dataspace_in_memory = H5Screate_simple(rank, dims, NULL);
              hdf5_dataspace_in_file = H5Screate_simple(rank, dims, NULL);

              H5Dread(hdf5_dataset, hdf5_datatype, hdf5_dataspace_in_memory, hdf5_dataspace_in_file, H5P_DEFAULT, tmp);

              for(int i = subgroup_range[fnr], k = 0; i < subgroup_range[fnr + 1]; i++)
                {
                  for(int j = 0; j < 3; j++)
                    SubGroupPos[i * 3 + j] = tmp[k++];
                }

              myfree(tmp);

              H5Dclose(hdf5_dataset);
            }

          if(blocknr == IO_SUB_LENTYPE)
            {
              char buf[1000];
              fof_subfind_get_dataset_name(blocknr, buf);

              hdf5_dataset = H5Dopen(hdf5_grp, buf);

              dims[0] = subgroup_range[fnr + 1] - subgroup_range[fnr];
              dims[1] = fof_subfind_get_values_per_blockelement(blocknr);
              if(dims[1] == 1)
                rank = 1;
              else
                rank = 2;

              hdf5_datatype = H5Tcopy(H5T_NATIVE_INT);
              int *tmp = mymalloc_clear("tmp", dims[0] * dims[1] * sizeof(int));

              hdf5_dataspace_in_memory = H5Screate_simple(rank, dims, NULL);
              hdf5_dataspace_in_file = H5Screate_simple(rank, dims, NULL);

              H5Dread(hdf5_dataset, hdf5_datatype, hdf5_dataspace_in_memory, hdf5_dataspace_in_file, H5P_DEFAULT, tmp);

              for(int i = subgroup_range[fnr], k = 0; i < subgroup_range[fnr + 1]; i++)
                {
                  for(int j = 0; j < 6; j++)
                    SubGroupLenType[i * 6 + j] = tmp[k++];
                }

              myfree(tmp);

              H5Dclose(hdf5_dataset);
            }

          if(blocknr == IO_FOF_LASTENTRY)
            break;

        }

      H5Gclose(hdf5_grp);
      H5Fclose(hdf5_file);
    }
}




int subfind_compare_originaltask_originalindex(const void *a, const void *b)
{
  if(((struct particle_data *) a)->OriginalTask < ((struct particle_data *) b)->OriginalTask)
    return -1;

  if(((struct particle_data *) a)->OriginalTask > ((struct particle_data *) b)->OriginalTask)
    return +1;

  if(((struct particle_data *) a)->OriginalIndex < ((struct particle_data *) b)->OriginalIndex)
    return -1;

  if(((struct particle_data *) a)->OriginalIndex > ((struct particle_data *) b)->OriginalIndex)
    return +1;

  return 0;
}

int subfind_compare_originalgrnr_minid(const void *a, const void *b)
{
  if(((struct particle_data *) a)->OriginalGrNr < ((struct particle_data *) b)->OriginalGrNr)
    return -1;

  if(((struct particle_data *) a)->OriginalGrNr > ((struct particle_data *) b)->OriginalGrNr)
    return +1;

  int type0 = (((struct particle_data *) a)->Type & FOF_PRIMARY_LINK_TYPES) ? 0 : 1;
  int type1 = (((struct particle_data *) a)->Type & FOF_PRIMARY_LINK_TYPES) ? 0 : 1;

  if(type0 < type1)
    return -1;

  if(type0 > type1)
    return +1;

  if(((struct particle_data *) a)->ID < ((struct particle_data *) b)->ID)
    return -1;

  if(((struct particle_data *) a)->ID > ((struct particle_data *) b)->ID)
    return +1;

  return 0;
}


void assign_group_numbers_based_on_catalogue(int ngroups_cat, int nsubgroups_cat)
{
  /* store original order of particle data */
  for(int i = 0; i < NumPart; i++)
    {
      P[i].OriginalIndex = i;
      P[i].OriginalTask = ThisTask;

      P[i].OriginalGrNr = 2000000000;
      P[i].OriginalSubNr = 2000000000;
    }

  /* sort by Fileorder (implies sort by type as well) */
  parallel_sort(P, NumPart, sizeof(struct particle_data), subfind_compare_FileOrder);

  for(int i = 0; i < NumPart - 1; i++)
    if(P[i].Type > P[i + 1].Type)
      terminate("P[i].Type > P[i+1].Type");


  int ntype_loc[6], *ntype_all;
  long long ntype_tot[6];

  for(int i = 0; i < 6; i++)
    ntype_loc[i] = 0;

  for(int i = 0; i < NumPart; i++)
    ntype_loc[P[i].Type]++;

  ntype_all = mymalloc("ntype_all", 6 * NTask * sizeof(int));
  MPI_Allgather(ntype_loc, 6, MPI_INT, ntype_all, 6, MPI_INT, MPI_COMM_WORLD);

  for(int i = 0; i < 6; i++)
    ntype_tot[i] = 0;

  for(int i = 0; i < NTask; i++)
    for(int j = 0; j < 6; j++)
      ntype_tot[j] += ntype_all[i * 6 + j];

  mpi_printf("FOF/SUBFIND REPROCESS: Begin FOF GrNr assignment\n");

  /* assign FOF numbers */
  for(int type = 0; type < 6; type++)
    {
      long long nprevious = 0;
      long long nprevious_range0 = 0;

      for(int task = 0; task < ThisTask; task++)
        nprevious_range0 += ntype_all[task * 6 + type];

      int kmin = 0;
      int kmax = 0;

      for(int t = 0; t < type; t++)
        kmin += ntype_all[ThisTask * 6 + t];

      kmax = kmin + ntype_all[ThisTask * 6 + type];

      int gr = 0;

      for(int k = kmin; k < kmax; k++)
        {
          while(gr < ngroups_cat && nprevious + GroupCat[gr].GroupLenType[type] <= nprevious_range0 + (k - kmin))
            {
              nprevious += GroupCat[gr].GroupLenType[type];
              gr++;
            }

          if(gr >= ngroups_cat)
            break;

          if(P[k].Type != type)
            terminate("gr=%d P[k=%d].Type=%d != type=%d   nprevious_types=%lld", gr, k, P[k].Type, type, nprevious);

          P[k].OriginalGrNr = gr;
        }
    }


  mpi_printf("FOF/SUBFIND REPROCESS: Begin SubhaloNr assignment\n");

  /* assign subgroup numbers */
  for(int type = 0; type < 6; type++)
    {
      long long nprevious = 0;
      long long nprevious_range0 = 0;

      for(int task = 0; task < ThisTask; task++)
        nprevious_range0 += ntype_all[task * 6 + type];

      int kmin = 0;
      int kmax = 0;

      for(int t = 0; t < type; t++)
        kmin += ntype_all[ThisTask * 6 + t];

      kmax = kmin + ntype_all[ThisTask * 6 + type];

      int gr = 0;
      int subnr = 0;
      long long nprevious_subs = nprevious;

      for(int k = kmin; k < kmax; k++)
        {
          while(gr < ngroups_cat && nprevious + GroupCat[gr].GroupLenType[type] <= nprevious_range0 + (k - kmin))
            {
              nprevious += GroupCat[gr].GroupLenType[type];
              gr++;

              subnr = 0;
              nprevious_subs = nprevious;
            }

          if(gr >= ngroups_cat)
            break;

          while(subnr < GroupCat[gr].Nsubs && nprevious_subs + SubGroupLenType[(GroupCat[gr].FirstSub + subnr) * 6 + type] <= nprevious_range0 + (k - kmin))
            {
              nprevious_subs += SubGroupLenType[(GroupCat[gr].FirstSub + subnr) * 6 + type];
              subnr++;
            }

          if(subnr < GroupCat[gr].Nsubs)
            {
              P[k].OriginalSubNr = GroupCat[gr].FirstSub + subnr;

              if(P[k].Type != type)
                terminate("P[k=%d].Type=%d != type=%d", k, P[k].Type, type);
            }
        }
    }

  mpi_printf("FOF/SUBFIND REPROCESS: Done with SubhaloNr assignment\n");


  /* assign MinID and MinIDTask such that previous FOF groups will be exactly recovered */
  for(int i = 0; i < NumPart; i++)
    {
      P[i].MinID = P[i].ID;
      P[i].MinIDTask = P[i].OriginalTask;
    }

  for(int grnr = 0; grnr < ngroups_cat; grnr++)
    {
      GroupCat[grnr].MinID = 0;
      GroupCat[grnr].MinIDTask = 0;
    }

  parallel_sort(P, NumPart, sizeof(struct particle_data), subfind_compare_originalgrnr_minid);

  mpi_printf("FOF/SUBFIND REPROCESS: Done with sorting according to originalgrnr\n");

  for(int i = 0; i < 6; i++)
    ntype_loc[i] = 0;

  for(int i = 0; i < NumPart; i++)
    ntype_loc[P[i].Type]++;

  MPI_Allgather(ntype_loc, 6, MPI_INT, ntype_all, 6, MPI_INT, MPI_COMM_WORLD);

  for(int i = 0; i < 6; i++)
    ntype_tot[i] = 0;

  for(int i = 0; i < NTask; i++)
    for(int j = 0; j < 6; j++)
      ntype_tot[j] += ntype_all[i * 6 + j];


  long long nprevious = 0;
  for(int task = 0; task < ThisTask; task++)
    for(int type = 0; type < 6; type++)
      nprevious += ntype_all[task * 6 + type];

  int gr = 0;
  int grlen = 0;
  for(int type = 0; type < 6; type++)
    grlen += GroupCat[gr].GroupLenType[type];

  long long nprevious_grp = 0;

  int count = 0;

  for(int k = 0; k < NumPart; k++)
    {
      while(gr < ngroups_cat && nprevious_grp + grlen <= nprevious + k)
        {
          nprevious_grp += grlen;
          gr++;

          if(gr >= ngroups_cat)
            break;

          grlen = 0;
          for(int type = 0; type < 6; type++)
            grlen += GroupCat[gr].GroupLenType[type];
        }

      if(gr >= ngroups_cat)
        break;

      if(nprevious_grp == nprevious + k)
        {
          GroupCat[gr].MinID = P[k].MinID;
          GroupCat[gr].MinIDTask = P[k].MinIDTask;
          count++;
        }
    }


  MPI_Allreduce(MPI_IN_PLACE, &count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(count != ngroups_cat)
    terminate("count=%d  ngroups_cat=%d", count, ngroups_cat);


  mpi_printf("FOF/SUBFIND REPROCESS: Done with determination of MinID/MinIDTask\n");

  MPI_Allreduce(MPI_IN_PLACE, GroupCat, ngroups_cat * sizeof(struct group_catalogue), MPI_BYTE, MPI_BOR, MPI_COMM_WORLD);


  /* and now distribute the values to the particles */
  nprevious = 0;
  for(int task = 0; task < ThisTask; task++)
    for(int type = 0; type < 6; type++)
      nprevious += ntype_all[task * 6 + type];

  gr = 0;
  grlen = 0;
  for(int type = 0; type < 6; type++)
    grlen += GroupCat[gr].GroupLenType[type];

  nprevious_grp = 0;

  for(int k = 0; k < NumPart; k++)
    {
      while(gr < ngroups_cat && nprevious_grp + grlen <= nprevious + k)
        {
          nprevious_grp += grlen;
          gr++;

          if(gr >= ngroups_cat)
            break;

          grlen = 0;
          for(int type = 0; type < 6; type++)
            grlen += GroupCat[gr].GroupLenType[type];
        }

      if(gr >= ngroups_cat)
        break;

      P[k].MinID = GroupCat[gr].MinID;
      P[k].MinIDTask = GroupCat[gr].MinIDTask;
    }

  /* revert to original order */
  parallel_sort(P, NumPart, sizeof(struct particle_data), subfind_compare_originaltask_originalindex);

  myfree(ntype_all);

  mpi_printf("FOF/SUBFIND REPROCESS: Done with restoring to original order\n");
}




void fof_append_group_properties(int snapnr, int NgroupsCat)
{
  int rank, fnr, bnr, ngroups, ngroups_already_done;
  int number_of_files, start = 0;
  int bytes_per_blockelement;
  hid_t hdf5_file, hdf5_grp;
  hid_t hdf5_dataset, hdf5_datatype, hdf5_dataspace_memory, hdf5_dataspace_in_file;
  hsize_t dims[2];
  enum fof_subfind_iofields blocknr;
  char fname[1000];
  int NgroupsDiff;

  if(NgroupsCat >= TotNgroups)
    NgroupsDiff = TotNgroups;
  else
    NgroupsDiff = NgroupsCat;   // in the case when the group catalogue is shorter i simply ignore all groups with higer GrNr - done.
  printf("FOF REPROCESS: TotNgroups %d, NgroupsCat %d\n", TotNgroups, NgroupsCat);

  number_of_files = get_number_of_group_catalogue_files(snapnr);        // fct def in subfind_so_reprocess.c

  for(fnr = 0; fnr < number_of_files; fnr++)
    {
      ngroups = get_number_of_groups_in_file(snapnr, fnr);      // fct def in subfind_so_reprocess.c

      if(ThisTask == 0)
        Group = &GroupAll[start];

      start += ngroups;
      NgroupsDiff -= ngroups;

      printf("FOF REPROCESS: determining groups in file %d: %d, NgroupsDiff: %d\n", fnr, ngroups, NgroupsDiff);

      // this is the case when the group catalogue is longer then the number of groups determined now
      if(NgroupsDiff < 0)
        {
          ngroups += NgroupsDiff;       // should be relevant in the last file only - contains number to fill up group cat exactly, all groups eith higher number simply ignored (but these are groups that could not be matched anyway)
          // if its not the last file: terminate and think of something smarter.
          if(fnr < (number_of_files - 1))
            terminate("FOF REPROCESS: run out of groups in file %d out of %d files in total", fnr, number_of_files);
        }

      if(All.NumFilesPerSnapshot > 1)
        sprintf(fname, "%sgroups_%03d/%s_%03d.%d.hdf5", All.OutputDir, snapnr, "fof_subhalo_tab", snapnr, fnr);
      else
        sprintf(fname, "%s%s_%03d.hdf5", All.OutputDir, "fof_subhalo_tab", snapnr);

      hdf5_file = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);
      if(hdf5_file < 0)
        terminate("cannot read initial conditions file %s", fname);


      if(ThisTask == 0)
        {
          // open
          if(All.SnapFormat == 3)
            {
              printf("FOF REPROCESS: extending group catalogue: '%s' (file %d of %d)\n", fname, fnr, All.NumFilesPerSnapshot);
              hdf5_grp = H5Gopen(hdf5_file, "/Group");
            }
          else
            {
              printf("can't open file `%s' to extend snapshot, as file not in hd5f format.\n", fname);
              terminate("file open error");
            }

          // write
          for(bnr = 0; bnr < 1000; bnr++)
            {
              blocknr = (enum fof_subfind_iofields) bnr;

              if(blocknr == IO_FOF_LASTENTRY)
                break;

              if(fof_additional_properties(blocknr))    // 1 if its an additional property
                {
                  bytes_per_blockelement = fof_subfind_get_bytes_per_blockelement(blocknr);

                  if(ngroups > 0)
                    {
                      char buf[1000];

                      fof_subfind_get_dataset_name(blocknr, buf);
                      printf("FOF/SUBFIND: writing block %d (%s)...\n", blocknr, buf);

                      switch (fof_subfind_get_datatype(blocknr))
                        {
                        default:
                          hdf5_datatype = H5Tcopy(H5T_NATIVE_INT);
                          break;
                        case 1:
#ifdef OUTPUT_IN_DOUBLEPRECISION
                          hdf5_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
#else
                          hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
#endif
                          break;
                        case 2:
                          hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT64);
                          break;
                        }

                      CommBuffer = mymalloc("CommBuffer", bytes_per_blockelement * ngroups);    // exactly adjusted

                      ngroups_already_done = 0; // in this file! not whole snapshot

                      fof_subfind_fill_write_buffer(blocknr, &ngroups_already_done, ngroups);

                      dims[0] = ngroups;
                      dims[1] = fof_subfind_get_values_per_blockelement(blocknr);
                      if(dims[1] == 1)
                        rank = 1;
                      else
                        rank = 2;

                      fof_subfind_get_dataset_name(blocknr, buf);

                      hdf5_dataspace_in_file = H5Screate_simple(rank, dims, NULL);      // creates space in file

                      if((hdf5_dataset = H5Dcreate(hdf5_grp, buf, hdf5_datatype, hdf5_dataspace_in_file, H5P_DEFAULT)) < 0)
                        hdf5_dataset = H5Dopen(hdf5_grp, buf);

                      hdf5_dataspace_memory = H5Screate_simple(rank, dims, NULL);

                      H5Dwrite(hdf5_dataset, hdf5_datatype, hdf5_dataspace_memory, hdf5_dataspace_in_file, H5P_DEFAULT, CommBuffer);

                      myfree(CommBuffer);

                      H5Sclose(hdf5_dataspace_memory);
                      H5Dclose(hdf5_dataset);
                      H5Sclose(hdf5_dataspace_in_file);
                      H5Tclose(hdf5_datatype);
                    }
                }

            }

          H5Gclose(hdf5_grp);
        }

      H5Fclose(hdf5_file);
    }
}

void fof_append_subgroup_properties(int snapnr, int NsubgroupsCat)
{
  int rank, fnr, bnr, nsubgroups, nsubgroups_already_done;
  int number_of_files, start = 0;
  int bytes_per_blockelement;
  hid_t hdf5_file, hdf5_grp;
  hid_t hdf5_dataset, hdf5_datatype, hdf5_dataspace_memory, hdf5_dataspace_in_file;
  hsize_t dims[2];
  enum fof_subfind_iofields blocknr;
  char fname[1000];
  int NsubgroupsDiff;

  if(NsubgroupsCat >= TotNsubgroups)
    NsubgroupsDiff = TotNsubgroups;
  else
    NsubgroupsDiff = NsubgroupsCat;     // in the case when the group catalogue is shorter i simply ignore all groups with higer GrNr - done.
  printf("SUBFIND REPROCESS: TotNsubgroups %d, NsubgroupsCat %d\n", TotNsubgroups, NsubgroupsCat);

  number_of_files = get_number_of_group_catalogue_files(snapnr);        // fct def in subfind_so_reprocess.c

  for(fnr = 0; fnr < number_of_files; fnr++)
    {
      nsubgroups = get_number_of_subgroups_in_file(snapnr, fnr);        // fct def in subfind_so_reprocess.c

      if(ThisTask == 0)
        SubGroup = &SubGroupAll[start];

      start += nsubgroups;
      NsubgroupsDiff -= nsubgroups;

      printf("SUBFIND REPROCESS: determining subhaloes in file %d: %d, NsubgroupsDiff: %d\n", fnr, nsubgroups, NsubgroupsDiff);

      // this is the case when the group catalogue is longer then the number of groups determined now
      if(NsubgroupsDiff < 0)
        {
          nsubgroups += NsubgroupsDiff; // should be relevant in the last file only - contains number to fill up group cat exactly, all groups eith higher number simply ignored (but these are groups that could not be matched anyway)
          // if its not the last file: terminate and think of something smarter.
          if(fnr < (number_of_files - 1))
            terminate("SUBFIND REPROCESS: run out of subhaloes in file %d out of %d files in total", fnr, number_of_files);
        }

      if(All.NumFilesPerSnapshot > 1)
        sprintf(fname, "%sgroups_%03d/%s_%03d.%d.hdf5", All.OutputDir, snapnr, "fof_subhalo_tab", snapnr, fnr);
      else
        sprintf(fname, "%s%s_%03d.hdf5", All.OutputDir, "fof_subhalo_tab", snapnr);

      hdf5_file = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);
      if(hdf5_file < 0)
        terminate("cannot read initial conditions file %s", fname);


      if(ThisTask == 0)
        {
          // open
          if(All.SnapFormat == 3)
            {
              printf("SUBFIND REPROCESS: extending group catalogue: '%s' (file %d of %d)\n", fname, fnr, All.NumFilesPerSnapshot);
              hdf5_grp = H5Gopen(hdf5_file, "/Subhalo");
            }
          else
            {
              printf("can't open file `%s' to extend snapshot, as file not in hd5f format.\n", fname);
              terminate("file open error");
            }

          // write
          for(bnr = 0; bnr < 1000; bnr++)
            {
              blocknr = (enum fof_subfind_iofields) bnr;

              if(blocknr == IO_FOF_LASTENTRY)
                break;

              if(sub_additional_properties(blocknr))    // 1 if its an additional property
                {
                  bytes_per_blockelement = fof_subfind_get_bytes_per_blockelement(blocknr);

                  if(nsubgroups > 0)
                    {
                      char buf[1000];

                      fof_subfind_get_dataset_name(blocknr, buf);
                      printf("FOF/SUBFIND: writing block %d (%s)...\n", blocknr, buf);

                      switch (fof_subfind_get_datatype(blocknr))
                        {
                        default:
                          hdf5_datatype = H5Tcopy(H5T_NATIVE_INT);
                          break;
                        case 1:
#ifdef OUTPUT_IN_DOUBLEPRECISION
                          hdf5_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
#else
                          hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
#endif
                          break;
                        case 2:
                          hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT64);
                          break;
                        }

                      CommBuffer = mymalloc("CommBuffer", bytes_per_blockelement * nsubgroups); // exactly adjusted

                      nsubgroups_already_done = 0;      // in this file! not whole snapshot

                      fof_subfind_fill_write_buffer(blocknr, &nsubgroups_already_done, nsubgroups);

                      dims[0] = nsubgroups;
                      dims[1] = fof_subfind_get_values_per_blockelement(blocknr);
                      if(dims[1] == 1)
                        rank = 1;
                      else
                        rank = 2;

                      fof_subfind_get_dataset_name(blocknr, buf);

                      hdf5_dataspace_in_file = H5Screate_simple(rank, dims, NULL);      // creates space in file

                      if((hdf5_dataset = H5Dcreate(hdf5_grp, buf, hdf5_datatype, hdf5_dataspace_in_file, H5P_DEFAULT)) < 0)
                        hdf5_dataset = H5Dopen(hdf5_grp, buf);

                      hdf5_dataspace_memory = H5Screate_simple(rank, dims, NULL);

                      H5Dwrite(hdf5_dataset, hdf5_datatype, hdf5_dataspace_memory, hdf5_dataspace_in_file, H5P_DEFAULT, CommBuffer);

                      myfree(CommBuffer);

                      H5Sclose(hdf5_dataspace_memory);
                      H5Dclose(hdf5_dataset);
                      H5Sclose(hdf5_dataspace_in_file);
                      H5Tclose(hdf5_datatype);
                    }
                }

            }

          H5Gclose(hdf5_grp);
        }

      H5Fclose(hdf5_file);
    }
}

void fof_collect_groups(void)
{
  int i;
  int *count = NULL;

  mpi_printf("FOF REPROCESS: Collecting groups...\n");

  if(ThisTask == 0)
    count = mymalloc_clear("count", NTask * sizeof(int));

  MPI_Gather(&Ngroups, 1, MPI_INT, count, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      char *buf = (char *) GroupAll;

      memcpy(buf, Group, Ngroups * sizeof(struct group_properties));
      buf += Ngroups * sizeof(struct group_properties);

      MPI_Status status;

      for(i = 1; i < NTask; i++)
        {
          if(count[i])
            MPI_Recv(buf,  count[i] * sizeof(struct group_properties), MPI_BYTE, i, TAG_N + i, MPI_COMM_WORLD, &status);

          buf += count[i] * sizeof(struct group_properties);
	}
    }
  else
    {
      if(Ngroups)
        MPI_Ssend(Group, Ngroups * sizeof(struct group_properties), MPI_BYTE, 0, TAG_N + ThisTask, MPI_COMM_WORLD);
    }

  if(ThisTask == 0)
    myfree(count);

  mpi_printf("FOF REPROCESS: Collecting groups done.\n");
}

void fof_collect_subgroups(void)
{
  int i;
  int *count = NULL;

  mpi_printf("FOF REPROCESS: Collecting subgroups...\n");

  if(ThisTask == 0)
    count = mymalloc_clear("count", NTask * sizeof(int));

  MPI_Gather(&Nsubgroups, 1, MPI_INT, count, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      char *buf = (char *) SubGroupAll;

      memcpy(buf, SubGroup, Nsubgroups * sizeof(struct subgroup_properties));
      buf += Nsubgroups * sizeof(struct subgroup_properties);

      MPI_Status status;

      for(i = 1; i < NTask; i++)
        {
          if(count[i])
            MPI_Recv(buf, count[i] * sizeof(struct subgroup_properties), MPI_BYTE, i, TAG_N + i, MPI_COMM_WORLD, &status);

          buf += count[i] * sizeof(struct subgroup_properties);
        }
    }
  else
    {
      if(Nsubgroups)
        MPI_Ssend(SubGroup, Nsubgroups * sizeof(struct subgroup_properties), MPI_BYTE, 0, TAG_N + ThisTask, MPI_COMM_WORLD);
    }

  if(ThisTask == 0)
    myfree(count);

  mpi_printf("FOF REPROCESS: Collecting subgroups done.\n");
}

int fof_additional_properties(enum fof_subfind_iofields blocknr)
{
  int present = 0;

  switch (blocknr)
    {
    default:
      present = 0;
      break;

      /* add new properties here   */
    case IO_FOF_J:
    case IO_FOF_JDM:
    case IO_FOF_JGAS:
    case IO_FOF_JSTARS:
    case IO_FOF_CMFRAC:
    case IO_FOF_CMFRACTYPE:
    case IO_FOF_EKIN:
    case IO_FOF_ETHR:
    case IO_FOF_EPOT:

    case IO_FOF_J_MEAN200:
    case IO_FOF_JDM_MEAN200:
    case IO_FOF_JGAS_MEAN200:
    case IO_FOF_JSTARS_MEAN200:
    case IO_FOF_MASSTYPE_MEAN200:
    case IO_FOF_LENTYPE_MEAN200:
    case IO_FOF_CMFRAC_MEAN200:
    case IO_FOF_CMFRACTYPE_MEAN200:
    case IO_FOF_J_CRIT200:
    case IO_FOF_JDM_CRIT200:
    case IO_FOF_JGAS_CRIT200:
    case IO_FOF_JSTARS_CRIT200:
    case IO_FOF_MASSTYPE_CRIT200:
    case IO_FOF_LENTYPE_CRIT200:
    case IO_FOF_CMFRAC_CRIT200:
    case IO_FOF_CMFRACTYPE_CRIT200:
    case IO_FOF_J_TOPHAT200:
    case IO_FOF_JDM_TOPHAT200:
    case IO_FOF_JGAS_TOPHAT200:
    case IO_FOF_JSTARS_TOPHAT200:
    case IO_FOF_MASSTYPE_TOPHAT200:
    case IO_FOF_LENTYPE_TOPHAT200:
    case IO_FOF_CMFRAC_TOPHAT200:
    case IO_FOF_CMFRACTYPE_TOPHAT200:
    case IO_FOF_J_CRIT500:
    case IO_FOF_JDM_CRIT500:
    case IO_FOF_JGAS_CRIT500:
    case IO_FOF_JSTARS_CRIT500:
    case IO_FOF_MASSTYPE_CRIT500:
    case IO_FOF_LENTYPE_CRIT500:
    case IO_FOF_CMFRAC_CRIT500:
    case IO_FOF_CMFRACTYPE_CRIT500:
    case IO_FOF_EPOT_CRIT200:
    case IO_FOF_EKIN_CRIT200:
    case IO_FOF_ETHR_CRIT200:
    case IO_FOF_EPOT_MEAN200:
    case IO_FOF_EKIN_MEAN200:
    case IO_FOF_ETHR_MEAN200:
    case IO_FOF_EPOT_TOPHAT200:
    case IO_FOF_EKIN_TOPHAT200:
    case IO_FOF_ETHR_TOPHAT200:
    case IO_FOF_EPOT_CRIT500:
    case IO_FOF_EKIN_CRIT500:
    case IO_FOF_ETHR_CRIT500:

      present = 1;
      break;
    }

  return present;
}

int sub_additional_properties(enum fof_subfind_iofields blocknr)
{
  int present = 0;

  switch (blocknr)
    {
    default:
      present = 0;
      break;

      /* add new properties here  */
    case IO_SUB_EKIN:
    case IO_SUB_ETHR:
    case IO_SUB_EPOT:
    case IO_SUB_J:
    case IO_SUB_JDM:
    case IO_SUB_JGAS:
    case IO_SUB_JSTARS:
    case IO_SUB_JINHALFRAD:
    case IO_SUB_JDMINHALFRAD:
    case IO_SUB_JGASINHALFRAD:
    case IO_SUB_JSTARSINHALFRAD:
    case IO_SUB_JINRAD:
    case IO_SUB_JDMINRAD:
    case IO_SUB_JGASINRAD:
    case IO_SUB_JSTARSINRAD:
    case IO_SUB_CMFRAC:
    case IO_SUB_CMFRACTYPE:
    case IO_SUB_CMFRACINHALFRAD:
    case IO_SUB_CMFRACTYPEINHALFRAD:
    case IO_SUB_CMFRACINRAD:
    case IO_SUB_CMFRACTYPEINRAD:
      present = 1;
      break;
    }

  return present;
}


int subfind_compare_FileOrder(const void *a, const void *b)
{
  if(((struct particle_data *) a)->FileOrder < ((struct particle_data *) b)->FileOrder)
    return -1;

  if(((struct particle_data *) a)->FileOrder > ((struct particle_data *) b)->FileOrder)
    return +1;

  return 0;
}

#endif
#endif
