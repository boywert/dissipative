/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/subfind/subfind_vars.c
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

#include "../allvars.h"

#ifdef SUBFIND

#include "../fof/fof.h"
#include "subfind.h"
#include "../domain.h"


double SubDomainCorner[3], SubDomainCenter[3], SubDomainLen, SubDomainFac;
double SubDomainInverseLen, SubDomainBigFac;

MyDouble GrCM[3];

int GrNr;
int NumPartGroup;

MPI_Comm SubComm;
int CommSplitColor;
int SubNTask, SubThisTask;
int SubTagOffset;

struct topnode_data *SubTopNodes;
struct local_topnode_data *Sub_LocTopNodes;

double SubTreeAllocFactor;

#if defined(SUBFIND) && defined(SUBFIND_EXTENDED_PROPERTIES)
int *NodeGrNr;
#endif

int *SubDomainTask;
int *SubDomainNodeIndex;
int *SubNextnode;
int SubNTopleaves;
int SubNTopnodes;

int SubTree_MaxPart;
int SubTree_NumNodes;
int SubTree_MaxNodes;
int SubTree_FirstNonTopLevelNode;
int SubTree_NumPartImported;
int SubTree_NumPartExported;
int SubTree_ImportedNodeOffset;
int SubTree_NextFreeNode;
struct NODE *SubNodes;
struct ExtNODE *SubExtNodes;
int *SubTree_ResultIndexList;
int *SubTree_Task_list;
unsigned long long *SubTree_IntPos_list;
MyDouble *SubTree_Pos_list;


int Ncollective;
int NprocsCollective;
int MaxNsubgroups = 0;
int MaxNgbs;
int MaxSerialGroupLen;

r2type *R2list;

int NumPaux;

struct paux_data *Paux;
struct proc_assign_data *ProcAssign;
struct subgroup_properties *SubGroup;
struct nearest_r2_data *R2Loc;
struct nearest_ngb_data *NgbLoc;
struct submp_data *submp;
struct cand_dat *candidates;
struct coll_cand_dat *coll_candidates;

#ifdef ADD_GROUP_PROPERTIES
struct subgroup_properties *SubGroupAll;
MyIDType *SubGroupMostBoundID;
MyFloat *SubGroupPos;
struct idsortlist *subgroupidlist;
#endif

#endif
