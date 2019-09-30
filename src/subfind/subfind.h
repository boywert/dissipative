/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/subfind/subfind.h
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

#ifndef SUBFIND_H
#define SUBFIND_H

#include "../allvars.h"
#include "../domain.h"



#define FIND_SMOOTHING_LENGTHS   0
#define FIND_TOTAL_DENSITIES     1

#define SUBFIND_SO_POT_CALCULATION_PARTICLE_NUMBER  10000

#define SUBFIND_GAL_RADIUS_FAC    2.0   /* for subfind metal calculation */

#if defined(SUBFIND) && defined(SUBFIND_EXTENDED_PROPERTIES)
extern int *NodeGrNr;
#endif

extern int GrNr;
extern int NumPartGroup;

extern struct topnode_data *SubTopNodes;
extern struct local_topnode_data *Sub_LocTopNodes;

extern int *SubDomainTask;
extern int *SubDomainNodeIndex;
extern int *SubNextnode;
extern int SubNTopleaves;
extern int SubNTopnodes;

extern int SubTree_MaxPart;
extern int SubTree_NumNodes;
extern int SubTree_MaxNodes;
extern int SubTree_FirstNonTopLevelNode;
extern int SubTree_NumPartImported;
extern int SubTree_NumPartExported;
extern int SubTree_ImportedNodeOffset;
extern int SubTree_NextFreeNode;
extern MyDouble *SubTree_Pos_list;
extern struct NODE *SubNodes;
extern struct ExtNODE *SubExtNodes;

extern double SubTreeAllocFactor;

extern int *SubTree_ResultIndexList;
extern int *SubTree_Task_list;
extern unsigned long long *SubTree_IntPos_list;


extern double SubDomainCorner[3], SubDomainCenter[3], SubDomainLen, SubDomainFac;
extern double SubDomainInverseLen, SubDomainBigFac;

extern MyDouble GrCM[3];

extern int Ncollective;
extern int NprocsCollective;
extern int MaxNsubgroups;
extern int MaxNgbs;
extern int MaxSerialGroupLen;
extern r2type *R2list;

extern int CommSplitColor;
extern MPI_Comm SubComm;

extern int SubNTask, SubThisTask;
extern int SubTagOffset;

#ifdef ADD_GROUP_PROPERTIES
extern MyIDType *SubGroupMostBoundID;
extern MyFloat  *SubGroupPos;
#endif

extern struct proc_assign_data
{
  int GrNr;
  int Len;
  int FirstTask;
  int NTask;
}
 *ProcAssign;

extern struct subgroup_properties
{
  int Len;
  int LenType[NTYPES];
  int GrNr;
  int SubNr;
  int SubParent;
  MyIDType SubMostBoundID;
  MyFloat Mass;
  MyFloat MassType[NTYPES];
  MyFloat SubVelDisp;
  MyFloat SubVmax;
  MyFloat SubVmaxRad;
  MyFloat SubHalfMassRad;
  MyFloat SubHalfMassRadType[NTYPES];
  MyFloat SubMassInRad;
  MyFloat SubMassInRadType[NTYPES];
  MyFloat SubMassInHalfRad;
  MyFloat SubMassInHalfRadType[NTYPES];
  MyFloat SubMassInMaxRad;
  MyFloat SubMassInMaxRadType[NTYPES];
  MyFloat Pos[3];
  MyFloat CM[3];
  MyFloat Vel[3];
  MyFloat Spin[3];

#ifdef MHD
  MyFloat Bfld_Halo, Bfld_Disk;
#endif

#ifdef SUBFIND_EXTENDED_PROPERTIES
  MyFloat Ekin, Epot, Ethr;
  MyFloat J[3], Jdm[3], Jgas[3], Jstars[3], CMFrac, CMFracType[NTYPES];
  MyFloat J_inRad[3], Jdm_inRad[3], Jgas_inRad[3], Jstars_inRad[3], CMFrac_inRad, CMFracType_inRad[NTYPES];
  MyFloat J_inHalfRad[3], Jdm_inHalfRad[3], Jgas_inHalfRad[3], Jstars_inHalfRad[3], CMFrac_inHalfRad, CMFracType_inHalfRad[NTYPES];
#endif

#ifdef USE_SFR
  MyFloat Sfr, SfrInRad, SfrInHalfRad, SfrInMaxRad, GasMassSfr;
#endif
#ifdef GFM_STELLAR_EVOLUTION
  MyFloat GasMassMetallicity;
  MyFloat GasMassMetallicityHalfRad;
  MyFloat GasMassMetallicityMaxRad;
  MyFloat GasMassMetals[GFM_N_CHEM_ELEMENTS];
  MyFloat GasMassMetalsHalfRad[GFM_N_CHEM_ELEMENTS];
  MyFloat GasMassMetalsMaxRad[GFM_N_CHEM_ELEMENTS];
  MyFloat StellarMassMetallicity;
  MyFloat StellarMassMetallicityHalfRad;
  MyFloat StellarMassMetallicityMaxRad;
  MyFloat StellarMassMetals[GFM_N_CHEM_ELEMENTS];
  MyFloat StellarMassMetalsHalfRad[GFM_N_CHEM_ELEMENTS];
  MyFloat StellarMassMetalsMaxRad[GFM_N_CHEM_ELEMENTS];
  MyFloat GasMassMetallicitySfr;
  MyFloat GasMassMetalsSfr[GFM_N_CHEM_ELEMENTS];
  MyFloat GasMassMetallicitySfrWeighted;
  MyFloat GasMassMetalsSfrWeighted[GFM_N_CHEM_ELEMENTS];
#ifdef GFM_DUST
  MyFloat GasMassDustMetallicity;
  MyFloat GasMassDustMetallicityHalfRad;
  MyFloat GasMassDustMetallicityMaxRad;
  MyFloat GasMassDustMetallicitySfr;
  MyFloat GasMassDustMetallicitySfrWeighted;
#endif
#endif
#ifdef BLACK_HOLES
  MyFloat BH_Mass;
  MyFloat BH_Mdot;
#endif
#ifdef GFM_WINDS
  MyFloat WindMass;
#endif
#ifdef SUBFIND_MEASURE_H2MASS
  MyFloat H2_Mass;
#endif
#ifdef GFM_STELLAR_PHOTOMETRICS
  MyFloat Magnitude_U, Magnitude_B, Magnitude_V, Magnitude_K;
  MyFloat Magnitude_g, Magnitude_r, Magnitude_i, Magnitude_z;
  MyFloat SurfaceBrightnessLimitRad;
  MyFloat SubMassInPhotRad;
#endif
} *SubGroup;


#ifdef ADD_GROUP_PROPERTIES
extern struct subgroup_properties *SubGroupAll;
#endif

#if defined(MHD) && defined(ADD_MAGNETIC_GROUP_PROPERTIES)
struct rlist_mhd
{
  double r;
  double b_egy;
};
#endif


extern struct nearest_r2_data
{
  double dist[2];
}
 *R2Loc;

extern struct nearest_ngb_data
{
  long long index[2];
  int count;
}
 *NgbLoc;

extern int NumPaux;

extern struct paux_data
 {
    int TaskOfGr;
    int LocGrIndex;
    unsigned char Type;
    unsigned char SofteningType;
    MyDouble Pos[3];
    MyDouble  Mass;
 }
 *Paux;

extern struct submp_data
{
  int index;
  int GrNr;
  int OldIndex;
  MyFloat DM_Density;
#ifdef ADD_GROUP_PROPERTIES
  int OriginalSubNr;
#endif
}
 *submp;

extern struct cand_dat
{
  int head;
  int len;
  int nsub;
  int rank, subnr, parent;
  int bound_length;
}
 *candidates;

extern struct coll_cand_dat
{
  long long head;
  long long rank;
  int len;
  int nsub;
  int subnr, parent;
  int bound_length;
}
 *coll_candidates;

typedef struct
{
  double rho;
#ifdef SUBFIND_CALC_MORE
  double vx, vy, vz;
  double v2;
#endif
} SubDMData;

void subfind_determine_sub_halo_properties(struct unbind_data *d, int num, struct subgroup_properties *subgroup, int grnr, int subnr, int parallel_flag, int nsubgroups_cat);
int subfind_ngb_treefind_density(MyDouble searchcenter[3], double hsml, int target, int *startnode, int mode, int *exportflag, int *exportnodecount, int *exportindex, SubDMData * sub_dm_data);
int subfind_treefind_collective_export_node_threads(int no, int i, int thread_id);
void subfind_domain_do_local_refine(int n, int *list);
void assign_group_numbers_based_on_catalogue(int ngroups_cat, int nsubgroups_cat);
int subfind_compare_rlist_mhd(const void *a, const void *b);
#endif
