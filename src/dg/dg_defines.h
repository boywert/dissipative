#ifdef TWODIMS
#define DG3D(...)
#define NOF_BASE_FUNCTIONS ((DEGREE_K + 1) * (DEGREE_K + 2) / 2)
#define NOF_INNER_QUAD_POINTS ((DEGREE_K + 1) * (DEGREE_K + 1))
#define NOF_OUTER_QUAD_POINTS (DEGREE_K + 1)
#define NOF_OUTER_FINE_QUAD_POINTS (2*NOF_OUTER_QUAD_POINTS)
#define NOF_INTERFACES 4
#define ONE_OR_TWO 1.
#define TWO_OR_FOUR 2.
#define DG_PROJ_NORM 4. //speedup: change to 0.25


#else //THREEDIMS
#define DG3D(...) __VA_ARGS__
#define NOF_BASE_FUNCTIONS ((DEGREE_K + 1) * (DEGREE_K + 2) * (DEGREE_K + 3) / 6)
#define NOF_INNER_QUAD_POINTS ((DEGREE_K + 1) * (DEGREE_K + 1) * (DEGREE_K + 1))
#define NOF_OUTER_QUAD_POINTS ((DEGREE_K + 1) * (DEGREE_K + 1))
#define NOF_OUTER_FINE_QUAD_POINTS (4*NOF_OUTER_QUAD_POINTS)
#define NOF_INTERFACES 6
#define ONE_OR_TWO 2.
#define TWO_OR_FOUR 4.
#define DG_PROJ_NORM 8. //speedup: change to 0.125
#endif

//general
#define NOF_QUAD_POINTS_1D (DEGREE_K + 1)
#define NOF_LOBATTO_POINTS_1D (ceil(0.5 * (DEGREE_K + 3)))

#if defined(REFINEMENT_SPLIT_CELLS) || defined(REFINEMENT_MERGE_CELLS)
#define DG_REFINEMENT
#endif

//order of the scheme
#define DG_ORDER (DEGREE_K + 1)

