/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/forcetree.h
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

#ifndef FORCETREE_H
#define FORCETREE_H

#ifndef INLINE_FUNC
#ifdef INLINE
#define INLINE_FUNC inline
#else
#define INLINE_FUNC
#endif
#endif


typedef struct
{
  MyDouble Pos[3];
#ifdef TREE_RAD_VEL
  MyFloat Vel[3];
  MyFloat Vth2;
#endif
  float OldAcc;
  unsigned char Type;
  unsigned char SofteningType;

  int Firstnode;
} gravdata_in;


typedef struct
 {
   MyFloat Acc[3];
 #ifdef EVALPOTENTIAL
   MyFloat Potential;
 #endif
 #ifdef OUTPUTGRAVINTERACTIONS
   int GravInteractions;
 #endif
 #ifdef TREE_RAD
   MyFloat Projection[NPIX];
 #ifdef TREE_RAD_H2
   MyFloat ProjectionH2[NPIX];
 #endif
 #ifdef TREE_RAD_CO
   MyFloat ProjectionCO[NPIX];
 #endif
 #endif

 #ifdef RADCOOL
   double Phios, Phins;
 #ifdef RADCOOL_HOTHALO
   double PhiT6, PhiT7, PhiT8;
 #endif
 #endif
#ifdef MODGRAV_EFF_MASS
   MyFloat ModgravAcc[3];
#endif

 } gravdata_out;


#ifdef LONG_X
#define STRETCHX (LONG_X)
#else
#define STRETCHX 1
#endif

#ifdef LONG_Y
#define STRETCHY (LONG_Y)
#else
#define STRETCHY 1
#endif

#ifdef LONG_Z
#define STRETCHZ (LONG_Z)
#else
#define STRETCHZ 1
#endif

#ifdef GRAVITY_TALLBOX

#if (GRAVITY_TALLBOX==0)
#define DBX 2
#define DBX_EXTRA 6
#define BOXX STRETCHY
#define BOXY STRETCHZ
#else
#define DBX 1
#define DBX_EXTRA 0
#endif

#if (GRAVITY_TALLBOX==1)
#define DBY 2
#define DBY_EXTRA 6
#define BOXX STRETCHX
#define BOXY STRETCHZ
#else
#define DBY 1
#define DBY_EXTRA 0
#endif

#if (GRAVITY_TALLBOX==2)
#define DBZ 2
#define DBZ_EXTRA 6
#define BOXX STRETCHX
#define BOXY STRETCHY
#else
#define DBZ 1
#define DBZ_EXTRA 0
#endif

#else

#define DBX 1
#define DBY 1
#define DBZ 1
#define DBX_EXTRA 0
#define DBY_EXTRA 0
#define DBZ_EXTRA 0
#endif



/*! length of lock-up table for short-range force kernel in TreePM algorithm */
#define NTAB 127



#if defined(SELFGRAVITY) && !defined(GRAVITY_NOT_PERIODIC) && defined(PERIODIC)


#define EN 64



#define ENX (DBX*STRETCHX*EN)
#define ENY (DBY*STRETCHY*EN)
#define ENZ (DBZ*STRETCHZ*EN)




extern MyFloat Ewd_fcorrx[ENX + 1][ENY + 1][ENZ + 1];
extern MyFloat Ewd_fcorry[ENX + 1][ENY + 1][ENZ + 1];
extern MyFloat Ewd_fcorrz[ENX + 1][ENY + 1][ENZ + 1];
extern MyFloat Ewd_potcorr[ENX + 1][ENY + 1][ENZ + 1];
extern double Ewd_fac_intp;

extern int NTreeInsert;

#endif


#define MAX_TREE_LEVEL        30
#define MAX_TREE_ALLOC_FACTOR 30.0

#define TAKE_NSLOTS_IN_ONE_GO      32


#define MAX_IMPACT_BEFORE_OPTIMIZATION 1.03


#define BITFLAG_TOPLEVEL                    0
#define BITFLAG_DEPENDS_ON_LOCAL_MASS       1
#define BITFLAG_DEPENDS_ON_EXTERN_MASS      2
#define BITFLAG_INTERNAL_TOPLEVEL           6
#define BITFLAG_MULTIPLEPARTICLES           7
#define BITFLAG_CONTAINS_GAS                10


#define BITFLAG_MASK  ((1<< BITFLAG_CONTAINS_GAS) + (1 << BITFLAG_MULTIPLEPARTICLES))



static inline unsigned long long force_double_to_int(double d)
{
  union
  {
    double d;
    unsigned long long ull;
  } u;
  u.d = d;
  return (u.ull & 0xFFFFFFFFFFFFFllu);
}

static inline double force_int_to_double(unsigned long long x)
{
  union
  {
    double d;
    unsigned long long ull;
  } u;
  u.d = 1.0;
  u.ull |= x;
  return u.d;
}

int tree_treefind_export_node_threads(int no, int target, int thread_id);
int force_treebuild(int npart, int optimized_domain_mapping, int insert_only_primary, int timebin);
int force_treebuild_construct(int npart, int optimized_domain_mapping, int insert_only_primary, int timebin);
int force_treebuild_insert_single_point(int i, unsigned long long *intpos, int th, unsigned char level);
int force_create_empty_nodes(int no, int topnode, int bits, int x, int y, int z);
void force_insert_pseudo_particles(void);
void force_update_node_recursive(int no, int sib, int father, int *last);
void force_exchange_topleafdata(void);
void force_treeupdate_toplevel(int no, int topnode, int bits, int x, int y, int z);
void force_treeallocate(int maxpart, int maxindex);
void force_treefree(void);
void dump_particles(void);
int force_add_empty_nodes(void);
void force_short_range_init(void);
int force_treeevaluate(gravdata_in *in, gravdata_out *out, int target, int mode, int thread_id, int numnodes, int *firstnode, int measure_cost_flag);
int force_treeevaluate_shortrange(gravdata_in *in, gravdata_out *out, int target, int mode, int thread_id, int numnodes, int *firstnode, int measure_cost_flag);
int force_treeevaluate_ewald_correction(int i, int mode, int thread_id);
int force_treeevaluate_direct(int target, int mode);
void force_assign_cost_values(void);
void force_update_node_recursive_sse(int no, int sib, int father, int *last);
void force_optimize_domain_mapping(void);
double force_get_current_balance(double *impact);
void force_get_global_cost_for_leavenodes(int nexport);
void forcetree_update_exported_potential_values(void);
void forcetest_ewald_init(void);

#endif



#ifdef RADCOOL_HOTHALO
double calculate_HH_temperature(double intU
#ifdef COOLING
                                , double intNe
#endif
  );
#endif
