/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/fof/fof.h
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

#ifndef FOF_H
#define FOF_H

#include "../allvars.h"

extern int Ngroups, NgroupsExt, MaxNgroups, TotNgroups, Nsubgroups, TotNsubgroups;
extern int Nids;
extern long long TotNids;

extern int fof_OldMaxPart;
extern int fof_OldMaxPartSph;
#ifdef GFM
extern int fof_OldMaxPartStar;
#endif
#ifdef BLACK_HOLES
extern int fof_OldMaxPartBHs;
#endif
#ifdef DUST_LIVE
extern int fof_OldMaxPartDust;
#endif


extern double LinkL;

extern unsigned char *flag_node_inside_linkinglength;

#define BITFLAG_INSIDE_LINKINGLENGTH    1


#ifndef FOF_SECONDARY_LINK_TARGET_TYPES
#define FOF_SECONDARY_LINK_TARGET_TYPES   FOF_PRIMARY_LINK_TYPES
#endif


extern struct group_properties
{
  int Len;
  MyIDType MinID;
  MyIDType MinIDTask;
  int GrNr;
  int LenType[NTYPES];
  MyFloat MassType[NTYPES];
  MyFloat Mass;
  MyDouble CM[3];
  MyFloat Vel[3];
  MyDouble Pos[3];
#ifdef FOF_FUZZ_SORT_BY_NEAREST_GROUP
  unsigned long long FuzzOffsetType[NTYPES];
#endif

#ifdef GFM_BIPOLAR_WINDS
#if (GFM_BIPOLAR_WINDS == 3)
  MyFloat DensGasMass;
  MyFloat DensGasCenter[3];
  MyFloat DensGasMomentum[3];
  MyFloat DensGasAngMomentum[3];
  MyFloat MinPotential;
  MyDouble Pos_MinPotential[3];
#else
  MyFloat GravAcc[3];
#endif
#endif

  MyDouble FirstPos[3];
#ifdef USE_SFR
  MyFloat Sfr;
#endif
#ifdef GFM_STELLAR_EVOLUTION
  MyFloat GasMassMetallicity;
  MyFloat GasMassMetals[GFM_N_CHEM_ELEMENTS];
  MyFloat StellarMassMetallicity;
  MyFloat StellarMassMetals[GFM_N_CHEM_ELEMENTS];
#ifdef GFM_DUST
  MyFloat GasMassDustMetallicity;
#endif
#endif
#ifdef BLACK_HOLES
  MyFloat BH_Mass;
  MyFloat BH_Mdot;
  MyFloat MaxDens;
  int index_maxdens, task_maxdens;
#ifdef BH_NF_RADIO
  MyFloat Min_BH_Potential;
  MyIDType ID_Min_BH_Potential;
  MyFloat XrayLum;
  MyFloat RadioLum;
#endif
#endif

#ifdef GFM_WINDS
  MyFloat WindMass;
#endif

#ifdef SUBFIND
  int TargetTask;               /* primary CPU responsible for this group */
  int Nsubs;
  int FirstSub;
  MyFloat M_Mean200, R_Mean200;
  MyFloat M_Crit200, R_Crit200;
  MyFloat M_Crit500, R_Crit500;
  MyFloat M_TopHat200, R_TopHat200;
#ifdef SUBFIND_EXTENDED_PROPERTIES
  MyFloat J_Mean200[3], JDM_Mean200[3], JGas_Mean200[3], JStars_Mean200[3], MassType_Mean200[NTYPES], CMFrac_Mean200, CMFracType_Mean200[NTYPES];
  MyFloat J_Crit200[3], JDM_Crit200[3], JGas_Crit200[3], JStars_Crit200[3], MassType_Crit200[NTYPES], CMFrac_Crit200, CMFracType_Crit200[NTYPES];
  MyFloat J_Crit500[3], JDM_Crit500[3], JGas_Crit500[3], JStars_Crit500[3], MassType_Crit500[NTYPES], CMFrac_Crit500, CMFracType_Crit500[NTYPES];
  MyFloat J_TopHat200[3], JDM_TopHat200[3], JGas_TopHat200[3], JStars_TopHat200[3], MassType_TopHat200[NTYPES], CMFrac_TopHat200, CMFracType_TopHat200[NTYPES];
  int LenType_Mean200[NTYPES], LenType_Crit200[NTYPES], LenType_Crit500[NTYPES], LenType_TopHat200[NTYPES];
  MyFloat J[3], JDM[3], JGas[3], JStars[3], CMFrac, CMFracType[NTYPES];
  MyFloat Ekin, Epot, Ethr;
  MyFloat Ekin_Crit200, Epot_Crit200, Ethr_Crit200;
  MyFloat Ekin_Crit500, Epot_Crit500, Ethr_Crit500;
  MyFloat Ekin_Mean200, Epot_Mean200, Ethr_Mean200;
  MyFloat Ekin_TopHat200, Epot_TopHat200, Ethr_TopHat200;
#endif
#endif

} *Group;


#ifdef ADD_GROUP_PROPERTIES
extern struct group_properties *GroupAll;

extern struct group_catalogue
{
  int GrNr;
  int Nsubs;
  int FirstSub;
  int GroupLenType[6];
  MyIDType MinID;
  int MinIDTask;
  MyIDType FirstSubMostBoundID;
#ifdef ADD_GROUP_PROPERTIES
  MyFloat R_Crit200;
  MyFloat R_Mean200;
  MyFloat R_TopHat200;
  MyFloat R_Crit500;
#endif
}
 *GroupCat;

extern struct idsortlist
{
  MyIDType MostBoundID;
  int index;
}
 *groupidlist, *subgroupidlist;

#endif


struct data_aux_sort
{
  int OriginTask, OriginIndex;
  int TargetTask, TargetIndex;
  int GrNr;
  int Type;
#ifdef FOF_FUZZ_SORT_BY_NEAREST_GROUP
  peanokey key;
#else
  MyIDType ID;
#endif
#if defined(RECOMPUTE_POTENTIAL_IN_SNAPSHOT) || defined(CALCULATE_QUANTITIES_IN_POSTPROCESS)
  MyIDType FileOrder;
#endif
#ifdef SUBFIND
  int SubNr;
  MyFloat DM_BindingEnergy;
#endif
};


extern struct fof_particle_list
{
  MyIDType MinID;
  int MinIDTask;
  int Pindex;
#ifdef ADD_GROUP_PROPERTIES
  int OriginalGrNr;
#endif
}
 *FOF_PList;


extern struct fof_group_list
{
  MyIDType MinID;
  int MinIDTask;
  int LocCount;
  int ExtCount;
#ifdef TRACER_PARTICLE
  int LocTrCount;
  int ExtTrCount;
#endif
  int GrNr;
#ifdef ADD_GROUP_PROPERTIES
  int OriginalGrNr;
#endif
}
 *FOF_GList;


extern struct id_list
{
  MyIDType ID;
  int GrNr;
  int Type;
#ifdef SUBFIND
  int SubNr;
  MyFloat BindingEgy;
#endif
}
 *ID_list;


extern struct bit_flags
{
  unsigned char Nonlocal:2, MinIDChanged:2, Marked:2, Changed:2;
} *Flags;

struct fof_local_sort_data
{
  int targetindex;
  int index;
};


extern struct fof_subfind_header
{
  int Ngroups;
  int Nsubgroups;
  int Nids;
  int TotNgroups;
  int TotNsubgroups;
  long long TotNids;
  int num_files;
  double time;
  double redshift;
  double HubbleParam;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  int flag_doubleprecision;
}
catalogue_header;


enum fof_subfind_iofields
{
  IO_FOF_LEN,
  IO_FOF_MTOT,
  IO_FOF_POS,
  IO_FOF_CM,
  IO_FOF_POSMINPOT,
  IO_FOF_VEL,
  IO_FOF_LENTYPE,
  IO_FOF_MASSTYPE,
  IO_FOF_SFR,
  IO_FOF_GASMETAL,
  IO_FOF_STARMETAL,
  IO_FOF_GASMETALELEMENTS,
  IO_FOF_STARMETALELEMENTS,
  IO_FOF_GASDUSTMETAL,
  IO_FOF_BHMASS,
  IO_FOF_BHMDOT,
  IO_FOF_WINDMASS,
  IO_FOF_XRAYLUM,
  IO_FOF_RADIOLUM,

  IO_FOF_M_MEAN200,
  IO_FOF_R_MEAN200,
  IO_FOF_M_CRIT200,
  IO_FOF_R_CRIT200,
  IO_FOF_M_TOPHAT200,
  IO_FOF_R_TOPHAT200,
  IO_FOF_M_CRIT500,
  IO_FOF_R_CRIT500,

#ifdef SUBFIND_EXTENDED_PROPERTIES
  IO_FOF_J_MEAN200,
  IO_FOF_JDM_MEAN200,
  IO_FOF_JGAS_MEAN200,
  IO_FOF_JSTARS_MEAN200,
  IO_FOF_MASSTYPE_MEAN200,
  IO_FOF_LENTYPE_MEAN200,
  IO_FOF_CMFRAC_MEAN200,
  IO_FOF_CMFRACTYPE_MEAN200,
  IO_FOF_J_CRIT200,
  IO_FOF_JDM_CRIT200,
  IO_FOF_JGAS_CRIT200,
  IO_FOF_JSTARS_CRIT200,
  IO_FOF_MASSTYPE_CRIT200,
  IO_FOF_LENTYPE_CRIT200,
  IO_FOF_CMFRAC_CRIT200,
  IO_FOF_CMFRACTYPE_CRIT200,
  IO_FOF_J_TOPHAT200,
  IO_FOF_JDM_TOPHAT200,
  IO_FOF_JGAS_TOPHAT200,
  IO_FOF_JSTARS_TOPHAT200,
  IO_FOF_MASSTYPE_TOPHAT200,
  IO_FOF_LENTYPE_TOPHAT200,
  IO_FOF_CMFRAC_TOPHAT200,
  IO_FOF_CMFRACTYPE_TOPHAT200,
  IO_FOF_J_CRIT500,
  IO_FOF_JDM_CRIT500,
  IO_FOF_JGAS_CRIT500,
  IO_FOF_JSTARS_CRIT500,
  IO_FOF_MASSTYPE_CRIT500,
  IO_FOF_LENTYPE_CRIT500,
  IO_FOF_CMFRAC_CRIT500,
  IO_FOF_CMFRACTYPE_CRIT500,
  IO_FOF_J,
  IO_FOF_JDM,
  IO_FOF_JGAS,
  IO_FOF_JSTARS,
  IO_FOF_CMFRAC,
  IO_FOF_CMFRACTYPE,
  IO_FOF_EKIN,
  IO_FOF_ETHR,
  IO_FOF_EPOT,
  IO_FOF_EPOT_CRIT200,
  IO_FOF_EKIN_CRIT200,
  IO_FOF_ETHR_CRIT200,
  IO_FOF_EPOT_MEAN200,
  IO_FOF_EKIN_MEAN200,
  IO_FOF_ETHR_MEAN200,
  IO_FOF_EPOT_TOPHAT200,
  IO_FOF_EKIN_TOPHAT200,
  IO_FOF_ETHR_TOPHAT200,
  IO_FOF_EPOT_CRIT500,
  IO_FOF_EKIN_CRIT500,
  IO_FOF_ETHR_CRIT500,
#endif


  IO_FOF_NSUBS,
  IO_FOF_FIRSTSUB,
  IO_FOF_FUZZOFFTYPE,
  IO_FOF_DENSLVEC,

  IO_SUB_LEN,
  IO_SUB_MTOT,
  IO_SUB_POS,
  IO_SUB_VEL,
  IO_SUB_LENTYPE,
  IO_SUB_MASSTYPE,
  IO_SUB_CM,
  IO_SUB_SPIN,
  IO_SUB_BFLD_HALO,
  IO_SUB_BFLD_DISK,

#ifdef SUBFIND_EXTENDED_PROPERTIES
  IO_SUB_EKIN,
  IO_SUB_ETHR,
  IO_SUB_EPOT,
  IO_SUB_J,
  IO_SUB_JDM,
  IO_SUB_JGAS,
  IO_SUB_JSTARS,
  IO_SUB_JINHALFRAD,
  IO_SUB_JDMINHALFRAD,
  IO_SUB_JGASINHALFRAD,
  IO_SUB_JSTARSINHALFRAD,
  IO_SUB_JINRAD,
  IO_SUB_JDMINRAD,
  IO_SUB_JGASINRAD,
  IO_SUB_JSTARSINRAD,
  IO_SUB_CMFRAC,
  IO_SUB_CMFRACTYPE,
  IO_SUB_CMFRACINHALFRAD,
  IO_SUB_CMFRACTYPEINHALFRAD,
  IO_SUB_CMFRACINRAD,
  IO_SUB_CMFRACTYPEINRAD,
#endif

  IO_SUB_VELDISP,
  IO_SUB_VMAX,
  IO_SUB_VMAXRAD,
  IO_SUB_HALFMASSRAD,
  IO_SUB_HALFMASSRADTYPE,
  IO_SUB_MASSINRAD,
  IO_SUB_MASSINHALFRAD,
  IO_SUB_MASSINMAXRAD,
  IO_SUB_MASSINRADTYPE,
  IO_SUB_MASSINHALFRADTYPE,
  IO_SUB_MASSINMAXRADTYPE,
  IO_SUB_IDMOSTBOUND,
  IO_SUB_GRNR,
  IO_SUB_PARENT,
  IO_SUB_SFR,
  IO_SUB_SFRINRAD,
  IO_SUB_SFRINHALFRAD,
  IO_SUB_SFRINMAXRAD,
  IO_SUB_GASMETAL,
  IO_SUB_GASMETALHALFRAD,
  IO_SUB_GASMETALMAXRAD,
  IO_SUB_GASMETALSFR,
  IO_SUB_GASMETALSFRWEIGHTED,
  IO_SUB_STARMETAL,
  IO_SUB_STARMETALHALFRAD,
  IO_SUB_STARMETALMAXRAD,
  IO_SUB_GASMETALELEMENTS,
  IO_SUB_GASMETALELEMENTSHALFRAD,
  IO_SUB_GASMETALELEMENTSMAXRAD,
  IO_SUB_GASMETALELEMENTSSFR,
  IO_SUB_GASMETALELEMENTSSFRWEIGHTED,
  IO_SUB_STARMETALELEMENTS,
  IO_SUB_STARMETALELEMENTSHALFRAD,
  IO_SUB_STARMETALELEMENTSMAXRAD,
  IO_SUB_GASDUSTMETAL,
  IO_SUB_GASDUSTMETALHALFRAD,
  IO_SUB_GASDUSTMETALMAXRAD,
  IO_SUB_GASDUSTMETALSFR,
  IO_SUB_GASDUSTMETALSFRWEIGHTED,
  IO_SUB_BHMASS,
  IO_SUB_BHMDOT,
  IO_SUB_WINDMASS,
  IO_SUB_H2MASS,
  IO_SUB_STELLARPHOTOMETRICS,
  IO_SUB_STELLARPHOTOMETRICSRAD,
  IO_SUB_STELLARPHOTOMETRICSMASSINRAD,
  IO_FOFSUB_IDS,
  IO_FOF_LASTENTRY
};


int fof_subfind_blockpresent(enum fof_subfind_iofields blocknr);
int fof_subfind_get_datatype(enum fof_subfind_iofields blocknr);
int fof_subfind_get_bytes_per_blockelement(enum fof_subfind_iofields blocknr);
int fof_subfind_get_particles_in_block(enum fof_subfind_iofields blocknr);
void fof_subfind_get_dataset_name(enum fof_subfind_iofields blocknr, char *label);
void fof_subfind_get_Tab_IO_Label(enum fof_subfind_iofields blocknr, char *label);
int fof_subfind_get_dataset_group(enum fof_subfind_iofields blocknr);
void fof_subfind_fill_write_buffer(enum fof_subfind_iofields blocknr, int *startindex, int pc);
int fof_subfind_get_values_per_blockelement(enum fof_subfind_iofields blocknr);


#ifdef ADD_GROUP_PROPERTIES
int fof_additional_properties(enum fof_subfind_iofields blocknr);
int sub_additional_properties(enum fof_subfind_iofields blocknr);
#endif

#endif
