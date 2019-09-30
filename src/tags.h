/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/tags.h
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

#define TAG_N             10    /*!< Various tags used for labelling MPI messages */
#define TAG_HEADER        11
#define TAG_PDATA         12
#define TAG_SPHDATA       13
#define TAG_KEY           14
#define TAG_DMOM          15
#define TAG_NODELEN       16
#define TAG_HMAX          17
#define TAG_GRAV_A        18
#define TAG_GRAV_B        19
#define TAG_DIRECT_A      20
#define TAG_DIRECT_B      21
#define TAG_HYDRO_A       22
#define TAG_HYDRO_B       23
#define TAG_NFORTHISTASK  24
#define TAG_PERIODIC_A    25
#define TAG_PERIODIC_B    26
#define TAG_PERIODIC_C    27
#define TAG_PERIODIC_D    28
#define TAG_NONPERIOD_A   29
#define TAG_NONPERIOD_B   30
#define TAG_NONPERIOD_C   31
#define TAG_NONPERIOD_D   32
#define TAG_POTENTIAL_A   33
#define TAG_POTENTIAL_B   34
#define TAG_DENS_A        35
#define TAG_DENS_B        36
#define TAG_LOCALN        37
#define TAG_BH_A          38
#define TAG_BH_B          39
#define TAG_SMOOTH_A      40
#define TAG_SMOOTH_B      41
#define TAG_ENRICH_A      42
#define TAG_CONDUCT_A     43
#define TAG_CONDUCT_B     44
#define TAG_FOF_A         45
#define TAG_FOF_B         46
#define TAG_FOF_C         47
#define TAG_FOF_D         48
#define TAG_FOF_E         49
#define TAG_FOF_F         50
#define TAG_FOF_G         51
#define TAG_HOTNGB_A      52
#define TAG_HOTNGB_B      53
#define TAG_GRAD_A        54
#define TAG_GRAD_B        55

// removed ifdefs here. we should never share tags. /Patrik
#define TAG_SE            56
#define TAG_METDATA       57

#define TAG_SEARCH_A      58
#define TAG_SEARCH_B      59

#define TAG_INJECT_A      61

#define TAG_PDATA_SPH     70
#define TAG_KEY_SPH       71

#define TAG_PDATA_STAR    72
#define TAG_STARDATA      73
#define TAG_KEY_STAR      74

#define TAG_PDATA_BH      75
#define TAG_BHDATA        76
#define TAG_KEY_BH        77

#define TAG_TRACERDATA    78

#define TAG_GRAVCOST_A    79
#define TAG_GRAVCOST_B    80

#ifdef OTVET
#define TAG_OTVET_A       81
#define TAG_OTVET_B       82
#endif

#define TAG_PM_FOLD       83
#define TAG_SHOCK_DATA    84

#define TAG_BARRIER       85
#define TAG_PART_DATA     95
#define TAG_NODE_DATA    105
#define TAG_RESULTS      115

#define TAG_PLOGS        116

#define TAG_MG_PM       117

#define TAG_PDATA_DUST    120
#define TAG_DUSTDATA      121
#define TAG_KEY_DUST      122
