/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/subfind/subfind_sort_kernels.c
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
#include <unistd.h>
#include <gsl/gsl_rng.h>

#include "../allvars.h"
#include "../proto.h"
#include "../fof/fof.h"
#include "../domain.h"
#include "subfind.h"


#ifdef SUBFIND

int subfind_compare_procassign_GrNr(const void *a, const void *b)
{
  if(((struct proc_assign_data *) a)->GrNr < ((struct proc_assign_data *) b)->GrNr)
    return -1;

  if(((struct proc_assign_data *) a)->GrNr > ((struct proc_assign_data *) b)->GrNr)
    return +1;

  return 0;
}

int subfind_compare_submp_GrNr_DM_Density(const void *a, const void *b)
{
  if(((struct submp_data *) a)->GrNr < ((struct submp_data *) b)->GrNr)
    return -1;

  if(((struct submp_data *) a)->GrNr > ((struct submp_data *) b)->GrNr)
    return +1;

#ifdef ADD_GROUP_PROPERTIES
  if(((struct submp_data *) a)->OriginalSubNr < ((struct submp_data *) b)->OriginalSubNr)
    return -1;

  if(((struct submp_data *) a)->OriginalSubNr > ((struct submp_data *) b)->OriginalSubNr)
    return +1;
#endif

  if(((struct submp_data *) a)->DM_Density > ((struct submp_data *) b)->DM_Density)
    return -1;

  if(((struct submp_data *) a)->DM_Density < ((struct submp_data *) b)->DM_Density)
    return +1;

  return 0;
}

int subfind_compare_submp_OldIndex(const void *a, const void *b)
{
  if(((struct submp_data *) a)->OldIndex < ((struct submp_data *) b)->OldIndex)
    return -1;

  if(((struct submp_data *) a)->OldIndex > ((struct submp_data *) b)->OldIndex)
    return +1;

  return 0;
}

int subfind_compare_ID_list(const void *a, const void *b)
{
  if(((struct id_list *) a)->GrNr < ((struct id_list *) b)->GrNr)
    return -1;

  if(((struct id_list *) a)->GrNr > ((struct id_list *) b)->GrNr)
    return +1;

  if(((struct id_list *) a)->SubNr < ((struct id_list *) b)->SubNr)
    return -1;

  if(((struct id_list *) a)->SubNr > ((struct id_list *) b)->SubNr)
    return +1;

  if(((struct id_list *) a)->Type < ((struct id_list *) b)->Type)
    return -1;

  if(((struct id_list *) a)->Type > ((struct id_list *) b)->Type)
    return +1;

  if(((struct id_list *) a)->BindingEgy < ((struct id_list *) b)->BindingEgy)
    return -1;

  if(((struct id_list *) a)->BindingEgy > ((struct id_list *) b)->BindingEgy)
    return +1;

  return 0;
}

int subfind_compare_SubGroup_GrNr_SubNr(const void *a, const void *b)
{
  if(((struct subgroup_properties *) a)->GrNr < ((struct subgroup_properties *) b)->GrNr)
    return -1;

  if(((struct subgroup_properties *) a)->GrNr > ((struct subgroup_properties *) b)->GrNr)
    return +1;

  if(((struct subgroup_properties *) a)->SubNr < ((struct subgroup_properties *) b)->SubNr)
    return -1;

  if(((struct subgroup_properties *) a)->SubNr > ((struct subgroup_properties *) b)->SubNr)
    return +1;

  return 0;
}


int subfind_compare_dist_rotcurve(const void *a, const void *b)
{
  if(((sort_r2list *) a)->r < ((sort_r2list *) b)->r)
    return -1;

  if(((sort_r2list *) a)->r > ((sort_r2list *) b)->r)
    return +1;

  return 0;
}

#if defined(MHD) && defined(ADD_MAGNETIC_GROUP_PROPERTIES)
int subfind_compare_rlist_mhd(const void *a, const void *b)
{
  if(((struct rlist_mhd *) a)->r < ((struct rlist_mhd *) b)->r)
    return -1;

  if(((struct rlist_mhd *) a)->r > ((struct rlist_mhd *) b)->r)
    return +1;

  return 0;
}
#endif

int subfind_compare_binding_energy(const void *a, const void *b)
{
  if(*((double *) a) > *((double *) b))
    return -1;

  if(*((double *) a) < *((double *) b))
    return +1;

  return 0;
}


int subfind_compare_serial_candidates_boundlength(const void *a, const void *b)
{
  if(((struct cand_dat *) a)->bound_length > ((struct cand_dat *) b)->bound_length)
    return -1;

  if(((struct cand_dat *) a)->bound_length < ((struct cand_dat *) b)->bound_length)
    return +1;

  if(((struct cand_dat *) a)->rank < ((struct cand_dat *) b)->rank)
    return -1;

  if(((struct cand_dat *) a)->rank > ((struct cand_dat *) b)->rank)
    return +1;

  return 0;
}

int subfind_compare_serial_candidates_rank(const void *a, const void *b)
{
  if(((struct cand_dat *) a)->rank < ((struct cand_dat *) b)->rank)
    return -1;

  if(((struct cand_dat *) a)->rank > ((struct cand_dat *) b)->rank)
    return +1;

  if(((struct cand_dat *) a)->len > ((struct cand_dat *) b)->len)
    return -1;

  if(((struct cand_dat *) a)->len < ((struct cand_dat *) b)->len)
    return +1;

  return 0;
}

int subfind_compare_serial_candidates_subnr(const void *a, const void *b)
{
  if(((struct cand_dat *) a)->subnr < ((struct cand_dat *) b)->subnr)
    return -1;

  if(((struct cand_dat *) a)->subnr > ((struct cand_dat *) b)->subnr)
    return +1;

  return 0;
}





int subfind_compare_coll_candidates_subnr(const void *a, const void *b)
{
  if(((struct coll_cand_dat *) a)->subnr < ((struct coll_cand_dat *) b)->subnr)
    return -1;

  if(((struct coll_cand_dat *) a)->subnr > ((struct coll_cand_dat *) b)->subnr)
    return +1;

  return 0;
}

int subfind_compare_coll_candidates_nsubs(const void *a, const void *b)
{
  if(((struct coll_cand_dat *) a)->nsub < ((struct coll_cand_dat *) b)->nsub)
    return -1;

  if(((struct coll_cand_dat *) a)->nsub > ((struct coll_cand_dat *) b)->nsub)
    return +1;

  return 0;
}

int subfind_compare_coll_candidates_boundlength(const void *a, const void *b)
{
  if(((struct coll_cand_dat *) a)->bound_length > ((struct coll_cand_dat *) b)->bound_length)
    return -1;

  if(((struct coll_cand_dat *) a)->bound_length < ((struct coll_cand_dat *) b)->bound_length)
    return +1;

  if(((struct coll_cand_dat *) a)->rank < ((struct coll_cand_dat *) b)->rank)
    return -1;

  if(((struct coll_cand_dat *) a)->rank > ((struct coll_cand_dat *) b)->rank)
    return +1;

  return 0;
}

int subfind_compare_coll_candidates_rank(const void *a, const void *b)
{
  if(((struct coll_cand_dat *) a)->rank < ((struct coll_cand_dat *) b)->rank)
    return -1;

  if(((struct coll_cand_dat *) a)->rank > ((struct coll_cand_dat *) b)->rank)
    return +1;

  if(((struct coll_cand_dat *) a)->len > ((struct coll_cand_dat *) b)->len)
    return -1;

  if(((struct coll_cand_dat *) a)->len < ((struct coll_cand_dat *) b)->len)
    return +1;

  return 0;
}

int subfind_fof_compare_ID(const void *a, const void *b)
{
  if(*((MyIDType *) a) < *((MyIDType *) b))
    return -1;

  if(*((MyIDType *) a) > *((MyIDType *) b))
    return +1;

  return 0;
}

#ifdef ADD_GROUP_PROPERTIES
int subfind_fof_compare_MostBoundID(const void *a, const void *b)
{
  if(((struct idsortlist *) a)->MostBoundID < ((struct idsortlist *) b)->MostBoundID)
    return -1;

  if(((struct idsortlist *) a)->MostBoundID > ((struct idsortlist *) b)->MostBoundID)
    return +1;

  return 0;
}


int subfind_compare_SubGroup_SubNr(const void *a, const void *b)
{
  if(((struct subgroup_properties *) a)->SubNr < ((struct subgroup_properties *) b)->SubNr)
    return -1;

  if(((struct subgroup_properties *) a)->SubNr > ((struct subgroup_properties *) b)->SubNr)
    return +1;

  return 0;
}
#endif

#endif
