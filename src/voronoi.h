/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/voronoi.h
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

#ifndef HAVE_H_VORONOI
#define HAVE_H_VORONOI

#ifdef VORONOI

#include <gmp.h>

#define STACKSIZE_TETRA        10000
#define MIN_ALLOC_NUMBER       1000
#define ALLOC_INCREASE_FACTOR  1.1
#define ALLOC_DECREASE_FACTOR  0.7
#define MAX_VORONOI_ITERATIONS 500


#define GENTLE_DEREFINE_FACTOR 1.2

#define USEDBITS 52

#if USEDBITS > 31
typedef signed long long int IntegerMapType;
void MY_mpz_set_si(mpz_t dest, signed long long int val);
void MY_mpz_mul_si(mpz_t prod, mpz_t mult, signed long long int val);
void MY_mpz_sub_ui(mpz_t prod, mpz_t mult, unsigned long long int val);
#else
typedef signed long int IntegerMapType;
#define MY_mpz_set_si mpz_set_si
#define MY_mpz_mul_si mpz_mul_si
#define MY_mpz_sub_ui mpz_sub_ui
#endif

#define DOUBLE_to_VORONOIINT(y)   ((IntegerMapType)(((*((long long *) &y)) & 0xFFFFFFFFFFFFFllu) >> (52 - USEDBITS)))

#ifdef NOINLINE
#define inline
#endif




/*
    Prerequisites for this function:
      sizeof(double)==sizeof(unsigned long long)
      doubles must be stored according to IEEE 754
*/
static inline IntegerMapType double_to_voronoiint(double d)
{
  union
  {
    double d;
    unsigned long long ull;
  } u;
  u.d = d;
  return (u.ull & 0xFFFFFFFFFFFFFllu) >> (52 - USEDBITS);
}

static inline double mask_voronoi_int(double x)
{
  union
  {
    double d;
    unsigned long long ull;
  } u;
  u.d = x;
  u.ull = u.ull & (~((1llu << (52 - USEDBITS)) - 1));
  return u.d;
}


#ifndef TWODIMS

#define EDGE_0 1                /* points 0-1 */
#define EDGE_1 2                /* points 0-2 */
#define EDGE_2 4                /* points 0-3 */
#define EDGE_3 8                /* points 1-2 */
#define EDGE_4 16               /* points 1-3 */
#define EDGE_5 32               /* points 2-3 */
#define EDGE_ALL 63

#else

#define EDGE_0 1                /* points 1-2 */
#define EDGE_1 2                /* points 0-2 */
#define EDGE_2 4                /* points 0-1 */
#define EDGE_ALL 7

#endif


#define HSML_INCREASE_FACTOR 1.3


#ifdef TWODIMS                  /* will only be compiled in 2D case */
#define DIMS 2
#else
#define DIMS 3
#endif



typedef struct
{
  double x, y, z;               // The 3-space position of the point
  MyIDType ID;
  int task;                     // The MPI task owning this cell
  int index;                    // The hydro quantity index of the cell
  int originalindex, timebin;
#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z) || defined(VORONOI_DYNAMIC_UPDATE)
  unsigned int image_flags;
#endif
#ifndef OPTIMIZE_MEMORY_USAGE
  double xx, yy, zz;
  IntegerMapType ix, iy, iz;
#endif
#ifdef DOUBLE_STENCIL
  MyFloat Hsml;
  int first_connection;
  int last_connection;
  char flag_primary_triangle;
#endif
}
point;


typedef struct tetra_data
{
  int p[DIMS + 1];              /* oriented tetrahedron points */
  int t[DIMS + 1];              /* adjacent tetrahedrons, always opposite to corresponding point */
  unsigned char s[DIMS + 1];    /* gives the index of the point in the adjacent tetrahedron that
                                   lies opposite to the common face */

  /* Note: if t[0] == -1, the tetrahedron has been deleted */
}
tetra;

typedef struct tetra_center_data
{
#ifndef OPTIMIZE_MEMORY_USAGE
  double cx, cy, cz;            /* describes circumcircle center */
#else
  MyFloat cx, cy, cz;
#endif
}
tetra_center;

typedef struct tri_data
{
  double p[DIMS + 1][DIMS];
  int owner;
}
triangle;

extern unsigned char *Edge_visited;

extern struct list_export_data
{
  unsigned int image_bits;
  int origin, index;
  int nextexport;
}
 *ListExports;

extern int Ninlist, MaxNinlist;


extern struct area_list_data
{
  int task, index;
  double darea;
}
 *AreaList;

extern int Narea, MaxNarea;

extern int NumGasInMesh;
extern int *List_InMesh;

extern struct list_P_data
{
  int firstexport, currentexport;

} *List_P;



typedef struct connection_data
{
  int task;
  int index;
  int image_flags;
  int next;

  int dp_index;                 /* this seems to be needed always the way voronoi_makeimage is implemented at the moment */
  int vf_index;                 /* index to the corresponding face */
#if defined(TETRA_INDEX_IN_FACE) || defined(MONOTONE_CONDUCTION)
  int dt_index;
#endif
  MyIDType ID;
}
connection;


/** This structure contains the points where a line segment intersects
    the tetrahedron faces and the internal voronoi faces. Is returned
    by calc_voronoi_intersections(). */
typedef struct intersection_list_data
{
  double s;                     /// the distance from the entry point (fraction of whole segment)
  point p;                      /// the intersection point
  int indA, indB;               /// the indices of the tetra points (0-4) defining the face
} intersection_list;


extern int CountInSphereTests, CountInSphereTestsExact;
extern int CountConvexEdgeTest, CountConvexEdgeTestExact;
extern int CountFlips, Count_1_to_3_Flips2d, Count_2_to_4_Flips2d;
extern int Count_1_to_4_Flips, Count_2_to_3_Flips, Count_3_to_2_Flips, Count_4_to_4_Flips;
extern int Count_EdgeSplits, Count_FaceSplits;
extern int Count_InTetra, Count_InTetraExact;
extern int Largest_N_DP_Buffer;

extern int Ninlist, MaxNinlist;


typedef struct individual_alloc_data
{
  double AllocFacNdp;
  double AllocFacNdt;
  double AllocFacNvf;
  double AllocFacNinlist;
  double AllocFacN_DP_Buffer;
  double AllocFacNflux;
  double AllocFacNradinflux;
#if defined(VORONOI_DYNAMIC_UPDATE)
  double AllocFacNvc;
#endif
}
mesh_alloc_facs;



typedef struct tessellation_data
{
  int Ndp;                      /* number of delaunay points */
  int MaxNdp;                   /* maximum number of delaunay points */
  point *DP;                    /* delaunay points */

  int Ndt;
  int MaxNdt;                   /* number of delaunary tetrahedra */
  tetra *DT;                    /* Delaunay tetrahedra */
  tetra_center *DTC;            /* circumcenters of delaunay tetrahedra */
  char *DTF;


  int Nvf;                      /* number of Voronoi faces */
  int MaxNvf;                   /* maximum number of Voronoi faces */
  face *VF;                     /* Voronoi faces */

  mesh_alloc_facs Indi;
}
tessellation;


extern tessellation Mesh, DeRefMesh;





extern int DPinfinity;

#ifdef VORONOI_DYNAMIC_UPDATE

extern int Nvc;                 /* number of connections */
extern int MaxNvc;              /* maximum number of connections */
extern int Largest_Nvc;
extern connection *DC;          /* Connections */
extern int FirstUnusedConnection;

#endif

void voronoi_derefinement_pairs_velvertex_corrections(void);

extern double CentralOffsetX, CentralOffsetY, CentralOffsetZ, ConversionFac;

void check_for_cut(double pp[3][3], int xaxis, int yaxis, int zaxis, double zval);
void drift_voronoi_face_half_a_step(void);
void recompute_circumcircles_and_faces(void);

void reprocess_edge_faces(tetra * t, int nr);

int derefine_add_point_and_split_tri(int q, triangle * trilist, int n, int max_n, double vol);

void make_3d_voronoi_listfaces_check_for_cut(tessellation * T, int tt, int nr, int xaxis, int yaxis, int zaxis, double zval);
void calculate_volume_changes(void);
void derefine_refine_process_edge(tessellation * T, double *vol, int tt, int nr);
void derefine_refine_compute_volumes(double *vol);

int derefine_refine_get_triangles(tessellation * T, int tt, int nr, point * dtip, triangle * trilist, int ntri, int max_n_tri);

void save_mass_flux_list(void);
void do_gravity_massflux_based_gravitywork(void);

void voronoi_save_old_gravity_forces(void);

void calculate_volume_changes_and_correct(void);

void create_mesh(void);
void mesh_setup_exchange(void);
void free_mesh(void);
void free_mesh_structures_not_needed_for_derefinement_refinement(void);
void free_all_remaining_mesh_structures(void);

void make_2d_voronoi_image(int num, int pixels_x, int pixels_y);
void make_2d_voronoi_image_zoomed(tessellation * T, int num, int pixels_x, int pixels_y, double xmin, double xmax, double ymin, double ymax);

void apply_area_list(void);
int area_list_data_compare(const void *a, const void *b);

void write_voronoi_mesh(tessellation * T, char *fname, int writeTask, int lastTask);
void finalize_output_of_voronoi_geometry(void);
void prepare_output_of_voronoi_geometry(void);
void initialize_and_create_first_tetra(tessellation * T);
void compute_voronoi_faces_and_volumes(void);
void get_line_segments(int sphp_index, int dp_index, double *segments, unsigned int *nof_elements, unsigned int max_elements);
double cross_section_plane_cell(int sphp_index, int dp_index, double *center, double *n);
void intersections_plane_cell(int sphp_index, int dp_index, double *center, double *n, double *polygon, unsigned int *nof_elements);
void intersection_plane_grid(double *center, double *n, const char *filename);
void process_edge_faces_and_volumes(tessellation * T, int tt, int nr);

int insert_point(tessellation * T, int pp, int ttstart);


void make_an_edge_split(tessellation * T, int tt0, int edge_nr, int count, int pp, int *ttlist);

void make_a_face_split(tessellation * T, int tt0, int face_nr, int pp, int tt1, int tt2, int qq1, int qq2);


double calculate_tetra_volume(point * p0, point * p1, point * p2, point * p3);
void make_a_4_to_4_flip(tessellation * T, int tt, int tip_index, int edge_nr);
double get_tri_volume(int i, triangle * trilist);


void make_a_1_to_4_flip(tessellation * T, int pp, int tt0, int tt1, int tt2, int tt3);

void make_a_3_to_2_flip(tessellation * T, int tt0, int tt1, int tt2, int tip, int edge, int bottom);

void make_a_2_to_3_flip(tessellation * T, int tt0, int tip, int tt1, int bottom, int qq, int tt2);


int get_tetra(tessellation * T, point * p, int *moves, int ttstart, int *flag, int *edgeface_nr);

int InTetra(tessellation * T, int tt, point * pp, int *edgeface_nr, int *nexttetra);
double InSphere(point * p0, point * p1, point * p2, point * p3, point * p);
void update_circumcircle(tessellation * T, int tt);
int test_tetra_orientation(point * p0, point * p1, point * p2, point * p3);
int test_intersect_triangle(point * p0, point * p1, point * p2, point * q, point * s);
double deter4(point * p0, point * p1, point * p2, point * p3);
double deter3(point * p0, point * p1, point * p2);
double deter4_orient(point * p0, point * p1, point * p2, point * p3);
double determinante3(double *a, double *b, double *c);
double deter_special(double *a, double *b, double *c, double *d);
int voronoi_ghost_search_alternative(tessellation * T);
int voronoi_exchange_evaluate(tessellation * T, int target, int mode, int *nexport, int *nsend_local);
int ngb_treefind_voronoi(tessellation * T, MyDouble searchcenter[3], MyFloat hsml, int target, int origin, int *startnode, int mode, int *nexport, int *nsend_local, int id);
void compute_circumcircles(tessellation * T);
int compute_max_delaunay_radius(void);
void check_for_min_distance(tessellation * T);
void check_links(tessellation * T);
void check_orientations(tessellation * T);
void check_tetras(tessellation * T, int npoints);
int voronoi_get_local_particles(void);
void check_for_vertex_crossings(void);
void voronoi_calculate_gravity_work_from_potential(int mode);
int convex_edge_test(tessellation * T, int tt, int tip, int *edgenr);
void voronoi_update_ghost_potential(void);

void do_hydro_calculations(void);
void update_cells_with_fluxes(void);
#ifdef RT_ADVECT
void rt_gradient_init(MyFloat * addr, MyFloat * addr_exch, double *addr_grad, int *addr_sourceid, int *addr_sourceid_exch, int type);
#endif
void calculate_gradients(void);
#ifdef RT_ADVECT
void rt_calculate_green_gauss_gradients(void);
#endif
void half_step_evolution(void);
void limit_gradient(double *d, double phi, double min_phi, double max_phi, MySingle *dphi);

#ifdef SECOND_DERIVATIVES
void hessian_init(MyFloat * addr_grad, MyFloat * addr_grad_exch, double *addr_hessian);
void calculate_green_gauss_hessian(tessellation * T);
void voronoi_exchange_hessians(void);
#endif

#if defined(SECOND_DERIVATIVES) && defined(RECONSTRUCT_GRADIENTS)
void face_time_advance_gradients(struct grad_data *delta_grad, struct state *st);
void face_space_extrapolate_gradients(struct grad_data *delta_grad, struct state *st);
void face_extrapolate_gradient(double *delta, double *hessian, double dx, double dy, double dz);
void face_add_gradient_extrapolations(struct state *st, struct grad_data *delta_grad_time, struct grad_data *delta_grad_space);
#endif

#if defined(VISCOSITY) || defined(THERMAL_CONDUCTION) || defined(TRACER_DIFFUSION)
void face_get_gradients(struct state *st_L, struct state *st_R, struct state_face *st_face, struct fluxes *flux);
#endif

#ifdef VISCOSITY
double local_get_dynvisc_coefficient(double x, double y, double z, double rho, double press);
double local_get_bulkvisc_coefficient(double x, double y, double z, double rho, double press);
double get_alpha_viscosity(double x, double y, double z, double rho, double press);
void face_get_viscous_fluxes(struct state_face *st_face, struct fluxes *flux, struct geometry *geom, double dyn_visc, double bulk_visc);
#ifdef SECOND_DERIVATIVES
void face_extrapolate_viscous_kick(struct state *st, double dyn_visc, double bulk_visc);
#endif
#endif

#ifdef THERMAL_CONDUCTION
void face_get_conduction_fluxes(struct state_face *st_face, struct fluxes *flux, struct geometry *geom);
#endif

#ifdef TRACER_DIFFUSION
void face_get_scalar_diffusion_fluxes(struct state_face *st_face, struct fluxes *flux, struct geometry *geom);
#endif



void exchange_primitive_variables(void);
void exchange_primitive_variables_and_gradients(void);
#ifdef RT_ADVECT
void rt_voronoi_exchange_primitive_variables(void);
void rt_voronoi_exchange_primitive_variables_and_gradients(void);
#endif
void voronoi_exchange_vertex_velocities(void);

int compare_primexch(const void *a, const void *b);

void apply_axisymmetric_source_terms(tessellation * T);


/* 2D voronoi routines */
void check_edge_and_flip_if_needed(tessellation * T, int ip, int it);
int get_triangle(tessellation * T, int pp, int *moves, int *degenerate_flag, int ttstart);
int InTriangle(point * p0, point * p1, point * p2, point * p);
double InCircle(point * p0, point * p1, point * p2, point * p);
double v2d_deter3(point * p0, point * p1, point * p2);
void make_a_1_to_3_flip(tessellation * T, int pp, int tt0, int tt1, int tt2);


double test_triangle_orientation(tessellation * T, int pp0, int pp1, int pp2);
void get_circle_center(point * a, point * b, point * c, double *x, double *y);
void do_special_dump(int num, int gradients_flag);

void make_a_2_to_4_flip(tessellation * T, int pp, int tt0, int tt1, int tt2, int tt3, int i0, int j0);

void dump_points(tessellation * T);

void set_integers_for_pointer(point * p);
#if !defined(ONEDIMS)
#ifndef OPTIMIZE_MEMORY_USAGE
static inline void set_integers_for_point(tessellation * T, int pp)
{
  point *p = &T->DP[pp];
  set_integers_for_pointer(p);
}
#else
static inline void get_integers_for_point(point * p, IntegerMapType ixyz[], double xyz[])
{
  xyz[0] = (p->x - CentralOffsetX) * ConversionFac + 1.0;
  xyz[1] = (p->y - CentralOffsetY) * ConversionFac + 1.0;
  xyz[2] = (p->z - CentralOffsetZ) * ConversionFac + 1.0;

  ixyz[0] = double_to_voronoiint(xyz[0]);
  ixyz[1] = double_to_voronoiint(xyz[1]);
  ixyz[2] = double_to_voronoiint(xyz[2]);

  xyz[0] = mask_voronoi_int(xyz[0]);
  xyz[1] = mask_voronoi_int(xyz[1]);
  xyz[2] = mask_voronoi_int(xyz[2]);
}
#endif
#else
void set_integers_for_point(tessellation * T, int pp);
#endif

/* quick function to compare a point to the infinity point */
static inline int isInfinity(point * p)
{
  //return gsl_isinf(p->x);       // DPinfinity has xyz set to INFINITY
  return p->x == MAX_DOUBLE_NUMBER;
}


int solve_linear_equations(double *m, double *res);


void check_triangles(tessellation * T, int npoints);


int InCircle_Quick(tessellation * T, int pp0, int pp1, int pp2, int pp);
int InCircle_Errorbound(tessellation * T, int pp0, int pp1, int pp2, int pp);

int InCircle_Exact(tessellation * T, int pp0, int pp1, int pp2, int pp);

int Orient2d_Exact(tessellation * T, int pp0, int pp1, int pp2);
int Orient2d_Quick(tessellation * T, int pp0, int pp1, int pp2);


void compute_axisymmetric_geometry_factors(tessellation * T);

int FindTriangle(tessellation * T, int tt, int pp, int *degnerate_flag, int *nexttetra);

int InSphere_Exact(point * p0, point * p1, point * p2, point * p3, point * p);
int InSphere_Quick(point * p0, point * p1, point * p2, point * p3, point * p);
int InSphere_Errorbound(point * p0, point * p1, point * p2, point * p3, point * p);
int InSphere_Gauss(point * p0, point * p1, point * p2, point * p3, point * p);

int Orient3d_Exact(point * p0, point * p1, point * p2, point * p3);
int Orient3d_Quick(point * p0, point * p1, point * p2, point * p3);
int Orient3d(point * p0, point * p1, point * p2, point * p3);

void make_3d_voronoi_listfaces(tessellation * T, int num, int xaxis, int yaxis, int zaxis, double zval);

int count_undecided_tetras(tessellation * T);
int ngb_treefind_ghost_search(tessellation * T, MyDouble searchcenter[3], MyDouble refpos[3], MyFloat hsml, MyFloat maxdist, int target, int origin, int mode, int thread_id, int numnodes,
                              int *firstnode);
int voronoi_ghost_search_evaluate(tessellation * T, int target, int mode, int q, int thread_id);

int voronoi_ghost_search(tessellation * T);

int find_next_voronoi_cell(tessellation * T, int cell, MyDouble p0[3], double dir[3], int previous, double *length);
int find_next_voronoi_cell2(tessellation * T, int cell, MyDouble p0[3], double dir[3], int previous, int task_of_previous, double *length);
double assert_contains(tessellation * T, int cell, MyDouble p0[3]);
void calc_delaunay_intersections(tessellation * T, int tt, point * pp0, point * pp1, intersection_list * list, int *nlist);
int decode_intersection_list(int i, intersection_list * list);

double distance_to_border(int cell);

#ifdef DVR_RENDER
void voronoi_make_new_tessellation(void);
void voronoi_restore_old_tessellation(void);
#endif

#endif /* VORONOI */

#endif /* HAVE_H_VORONOI */
