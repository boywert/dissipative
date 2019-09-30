/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/amr/amr.h
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

#ifndef AMR_H
#define AMR_H

#define AMR_NONE (0)
#define AMR_REFINE (1<<0)
#define AMR_SMOOTH (1<<1)
#define AMR_REFM (1<<2) //cell already added to list of refinement candidates

#define AMR_REF1_LEFT (1<<3)
#define AMR_REF1_RIGHT (1<<4)
#define AMR_REF1 (AMR_REF1_LEFT | AMR_REF1_RIGHT)

#define AMR_REF2_FRONT_LEFT (1<<5)
#define AMR_REF2_FRONT_RIGHT (1<<6)
#define AMR_REF2_BACK_LEFT (1<<7)
#define AMR_REF2_BACK_RIGHT (1<<8)
#define AMR_REF2_FRONT (AMR_REF2_FRONT_LEFT | AMR_REF2_FRONT_RIGHT)
#define AMR_REF2_BACK (AMR_REF2_BACK_LEFT | AMR_REF2_BACK_RIGHT)
#define AMR_REF2_LEFT (AMR_REF2_FRONT_LEFT | AMR_REF2_BACK_LEFT)
#define AMR_REF2_RIGHT (AMR_REF2_FRONT_RIGHT | AMR_REF2_BACK_RIGHT)
#define AMR_REF2 (AMR_REF2_FRONT | AMR_REF2_BACK)

#define AMR_REF3_TOP_FRONT_LEFT (1<<9)
#define AMR_REF3_TOP_FRONT_RIGHT (1<<10)
#define AMR_REF3_TOP_BACK_LEFT (1<<11)
#define AMR_REF3_TOP_BACK_RIGHT (1<<12)
#define AMR_REF3_BOTTOM_FRONT_LEFT (1<<13)
#define AMR_REF3_BOTTOM_FRONT_RIGHT (1<<14)
#define AMR_REF3_BOTTOM_BACK_LEFT (1<<15)
#define AMR_REF3_BOTTOM_BACK_RIGHT (1<<16)
#define AMR_REF3_TOP_FRONT (AMR_REF3_TOP_FRONT_LEFT | AMR_REF3_TOP_FRONT_RIGHT)
#define AMR_REF3_TOP_BACK (AMR_REF3_TOP_BACK_LEFT | AMR_REF3_TOP_BACK_RIGHT)
#define AMR_REF3_TOP (AMR_REF3_TOP_FRONT | AMR_REF3_TOP_BACK)
#define AMR_REF3_BOTTOM_FRONT (AMR_REF3_BOTTOM_FRONT_LEFT | AMR_REF3_BOTTOM_FRONT_RIGHT)
#define AMR_REF3_BOTTOM_BACK (AMR_REF3_BOTTOM_BACK_LEFT | AMR_REF3_BOTTOM_BACK_RIGHT)
#define AMR_REF3_BOTTOM (AMR_REF3_BOTTOM_FRONT | AMR_REF3_BOTTOM_BACK)
#define AMR_REF3 (AMR_REF3_TOP | AMR_REF3_BOTTOM)
#define AMR_REF3_LEFT ( AMR_REF3_TOP_FRONT_LEFT | AMR_REF3_TOP_BACK_LEFT | AMR_REF3_BOTTOM_FRONT_LEFT | AMR_REF3_BOTTOM_BACK_LEFT)
#define AMR_REF3_RIGHT ( AMR_REF3_TOP_FRONT_RIGHT | AMR_REF3_TOP_BACK_RIGHT | AMR_REF3_BOTTOM_FRONT_RIGHT | AMR_REF3_BOTTOM_BACK_RIGHT)
#define AMR_REF3_FRONT (AMR_REF3_BOTTOM_FRONT | AMR_REF3_TOP_FRONT)
#define AMR_REF3_BACK (AMR_REF3_BOTTOM_BACK | AMR_REF3_TOP_BACK)
#define AMR_REF3_FRONT_LEFT (AMR_REF3_TOP_FRONT_LEFT | AMR_REF3_BOTTOM_FRONT_LEFT)
#define AMR_REF3_FRONT_RIGHT (AMR_REF3_TOP_FRONT_RIGHT | AMR_REF3_BOTTOM_FRONT_RIGHT)
#define AMR_REF3_BACK_LEFT (AMR_REF3_TOP_BACK_LEFT | AMR_REF3_BOTTOM_BACK_LEFT)
#define AMR_REF3_BACK_RIGHT (AMR_REF3_TOP_BACK_RIGHT | AMR_REF3_BOTTOM_BACK_RIGHT)

#define AMR_REF_ANY (AMR_REF1 | AMR_REF2 | AMR_REF3)


#define AMR_REFINE_EXPAND_LEFT (1<<17)
#define AMR_REFINE_EXPAND_RIGHT (1<<18)
#define AMR_REFINE_EXPAND_FRONT (1<<19)
#define AMR_REFINE_EXPAND_BACK (1<<20)
#define AMR_REFINE_EXPAND_TOP (1<<21)
#define AMR_REFINE_EXPAND_BOTTOM (1<<22)

#define AMR_REFINE_EXPAND (AMR_REFINE_EXPAND_LEFT | AMR_REFINE_EXPAND_RIGHT | AMR_REFINE_EXPAND_FRONT | AMR_REFINE_EXPAND_BACK | AMR_REFINE_EXPAND_TOP | AMR_REFINE_EXPAND_BOTTOM)

#define AMR_REFINE_SMOOTH_A (1<<23)
#define AMR_REFINE_SMOOTH_B (1<<24)

#define AMR_REFINE_CRITERION (1<<25)

#define AMR_REFINE_ABORT (1<<26)


#define AMR_MAX_REFLEVEL    40  //to be on the safe side with double precision

#define MIN_ALLOC_NUMBER       1000
#define ALLOC_INCREASE_FACTOR  1.1

#ifndef REFINEMENT
#define AMR_STATIC_MESH
#endif


extern double amr_length[];
extern double amr_area[];
extern double amr_volume[];

typedef struct
{
  int task, index;
}
list_export_data;


/*! This is the AMR version of DP
 */
typedef struct point
{
  double x, y, z;               // The 3-space position of the point

  MyIDType ID;
  int task;                     // The MPI task owning this cell
  int index;                    // The hydro quantity index of the cell
  int originalindex;
  int timebin;
#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z) || defined(AMR_CONNECTIONS)
  unsigned int image_flags;
#endif
  int father; //points to the father node
  int level; // amr refinement level

  int neighbors[2*NUMDIMS];
  int nextinlevel;
  int previnlevel;
}
point;

typedef struct amr_node_data
{
  MyFloat mass;
  MyFloat momentum[3];
  MyFloat etherm;

#ifdef DG
  MyDouble Weights[NOF_BASE_FUNCTIONS][5];
#endif
  
#ifdef FLD
#ifdef FLD_MG
  MyFloat b;
  MyFloat residuum;
  MyFloat kappa;
  MyFloat w[2*NUMDIMS+1];
#endif
#endif
}
amr_node_data;

typedef struct point_exchange
{
  double x, y, z;               // The 3-space position of the point
  MyIDType ID;
  int timebin;
  int originalindex;
  int level;
} point_exchange;


typedef struct node_exchange
{
  double x, y, z;               // The 3-space position of the point
  int originalindex;
  int level;
} node_exchange;

typedef struct individual_alloc_data
{
  double AllocFacNdp;
  double AllocFacNdt;           //FIXME remove??
  double AllocFacNnodes;
  double AllocFacNgnodes;

  double AllocFacNvf;
  double AllocFacNinlist;
  double AllocFacN_DP_Buffer;
  double AllocFacNflux;
  double AllocFacNradinflux;

#ifdef AMR_CONNECTIONS
  double AllocFacNvc;
#endif
}
mesh_alloc_facs;

typedef struct tessellation_data
{
  int LastDP;                   /* index of the last non ghost DP */
  int Ndp;                      /* number of points */
  int MaxNdp;                   /* maximum number of points */
  point *DP;                    /* points */

  int Nvf;                      /* number of faces */
  int MaxNvf;                   /* maximum number of faces */
  face *VF;                     /* faces */


  //struct particle_data *p_node_data;
  //struct sph_particle_data *sph_node_data;
  /* num, max for these arrays are the same as for amr_nodes */

  list_export_data *expList;
  int NexpList, MaxNexpList;

  list_export_data *expListNode;
  int NexpListNode, MaxNexpListNode;

  int Node_nimport, Node_nexport, *Node_Send_offset, *Node_Send_count, *Node_Recv_count, *Node_Recv_offset;

  mesh_alloc_facs Indi;
  int allocated;

  int Nghost_dp;
  int Nghost_nodes;
  int nodes_total;

  int* original_index;
  int* original_task;

  int lastinlevel[AMR_MAX_REFLEVEL + 1];

  int minlevel;
  int maxlevel;

#ifdef REFINEMENT
  int* refcells;
  int Nrefine;

  int* derefcand;
  int Nderefine;

  int Refined;
  int Derefined;

  MyIDType IDNew;

  int* Refflag;
#endif
}
tessellation;

extern tessellation Mesh;



extern int amr_amrtree;

#if defined(LONG_X) || defined(LONG_Y) ||defined(LONG_Z)
extern int amr_last_long_level[3];
#endif


#ifdef AMR_CONNECTIONS
typedef struct connection_data
{
  int task;
  int index;
  int image_flags;
  int next;

  int dp_index; /* this seems to be needed always the way voronoi_makeimage is implemented at the moment */
  int vf_index; /* index to the corresponding face */
  MyIDType ID;
}
connection;


extern int Nvc;                        /* number of connections */
extern int MaxNvc;                     /* maximum number of connections */
extern int Largest_Nvc;
extern connection *DC;                 /* Connections */
extern int FirstUnusedConnection;
#endif

#endif /* AMR_H_ */
