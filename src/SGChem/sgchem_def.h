#ifndef SG_HEADER_FLAG
#define SG_HEADER_FLAG

/* Param. definitions for SG chemistry */

#ifndef ABHE
#if CHEMISTRYNETWORK == 1
#define ABHE 0.079
#else
#define ABHE 0.1
#endif
#endif

/* Default values from Network 15; for the species used in the other networks, we
 * need to set the values explicitly below, to avoid CPP complaining about macro
 * redefinition
 */
#define ICHX 3
#define IOHX 4
#define IHCOP 6
#define IMP   8
#define ICATOM  11
#define IOATOM  12
#define IMATOM  13

#ifdef SGCHEM_VARIABLE_Z
#define SGCHEM_NUM_ELEMS 5
#else
#define SGCHEM_NUM_ELEMS 0
#endif

#if CHEMISTRYNETWORK == 1
#define SGCHEM_NUM_SPECIES 6
#define SGCHEM_NUM_ADVECTED_SPECIES 9
#define SGCHEM_NUM_THERMAL_RATES 12
#define IH2   0
#define IHP   1
#define IDP   2
#define IHD   3
#define IHEP  4
#define IHEPP 5
#define IHATOM  6
#define IDATOM  7
#define IHEATOM 8
/* Not used, but need definitions */
#define ICP -1
#define ICO -1
#endif

#if CHEMISTRYNETWORK == 4
#define SGCHEM_NUM_SPECIES 2
#define SGCHEM_NUM_ADVECTED_SPECIES 3
#define SGCHEM_NUM_THERMAL_RATES 18
#define IH2 0
#define IHP 1
#define IHATOM 2
/* Not used, but need definitions */
#define ICP  -1
#define ICO  -1
#define IHEP -1
#define IHEATOM -1
#endif

#if CHEMISTRYNETWORK == 5
#define SGCHEM_NUM_SPECIES 3
#define SGCHEM_NUM_ADVECTED_SPECIES 5
#define SGCHEM_NUM_THERMAL_RATES 21
#define IH2 0
#define IHP 1
#define ICO 2
#define IHATOM 3
#define ICP 4
/* Not used, but need definitions */
#define IHEP -1
#define IHEATOM -1
#endif

#if CHEMISTRYNETWORK == 15
#define SGCHEM_NUM_SPECIES 9
#define SGCHEM_NUM_ADVECTED_SPECIES 14
#define SGCHEM_NUM_THERMAL_RATES 22
#define IH2 0
#define IHP 1
#define ICP 2
#define ICO  5
#define IHEP  7
#define IHATOM  9
#define IHEATOM 10
/* For the rest, we use the defaults above */
#ifdef MCMA
#define NELEM_CMA 3
#define NSPEC_CMA 10
#endif
#endif

/* Temperature always goes after species */
#define ITMP SGCHEM_NUM_SPECIES

#endif
