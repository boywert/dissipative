/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/mesh.h
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


#ifndef MESH_H
#define MESH_H

#ifdef DG
#include "dg/dg_defines.h"
#endif

#define SCALAR_TYPE_PASSIVE   0 /* only advection */
#define SCALAR_TYPE_SPECIES   1 /* species are normalised to guarantee sum{species}=1 */
#define SCALAR_TYPE_NORMALIZE 2 /* the same normalisation factor as for species is applied, but no contribution to sum{species} */
#ifdef SGCHEM
#define SCALAR_TYPE_CELEM     3 /* Only used when SGCHEM_VARIABLE_Z is defined */
#define SCALAR_TYPE_OELEM     4 /* Only used when SGCHEM_VARIABLE_Z is defined */
#define SCALAR_TYPE_MELEM     5 /* Only used when SGCHEM_VARIABLE_Z is defined */
#define SCALAR_TYPE_HYDROGEN  6
#define SCALAR_TYPE_H2        7
#if CHEMISTRYNETWORK == 1
#define SCALAR_TYPE_HD        8
#define SCALAR_TYPE_DEUTERIUM 9
#define SCALAR_TYPE_HELIUM    10
#elif CHEMISTRYNETWORK == 15
#define SCALAR_TYPE_CARBON    8
#define SCALAR_TYPE_HELIUM    9
#define SCALAR_TYPE_METAL     10
#define SCALAR_TYPE_MCMA      11
#else
#define SCALAR_TYPE_CARBON    8
#endif
#endif


#define REFL_X_FLAGS 115043766
#define REFL_Y_FLAGS 132379128
#define REFL_Z_FLAGS 134217216

#define OUTFLOW_X (1<<27)
#define OUTFLOW_Y (1<<28)
#define OUTFLOW_Z (1<<29)

#if defined MAXSCALARS
extern struct scalar_elements
{
  int type;                     /* scalar type, determines whether a normalization is applied */
  size_t offset;                /* offset of the primitive quantity in the SphP struct */
  size_t offset_mass;           /* offset of the conserved quantity in the SphP struct */
}
scalar_elements[MAXSCALARS];

extern struct scalar_index
{
#ifdef REFINEMENT_HIGH_RES_GAS
  int HighResMass;
#endif
#ifdef GFM_STELLAR_EVOLUTION
  int Metallicity;
#endif
#ifdef COSMIC_RAYS
  int CR_Energy;
#endif
#ifdef MHD_THERMAL_ENERGY_SWITCH
  int Etherm;
#endif
} ScalarIndex;

extern int N_Scalar;            /* number of registered scalars */
#endif

#define GRADIENT_TYPE_NORMAL   0
#define GRADIENT_TYPE_VELX     1
#define GRADIENT_TYPE_VELY     2
#define GRADIENT_TYPE_VELZ     3
#define GRADIENT_TYPE_DENSITY  4
#define GRADIENT_TYPE_PRESSURE 5
#define GRADIENT_TYPE_UTHERM   6
#define GRADIENT_TYPE_AX       7
#define GRADIENT_TYPE_AY       8
#define GRADIENT_TYPE_AZ       9
#define GRADIENT_TYPE_FLD     10
#define GRADIENT_TYPE_RTF     11

extern struct grad_elements
{
  int type;                     /* gradient type, ensures special treatment for velocities and speed of sound */
  size_t offset;                /* offset of the quantity in the SphP struct */
  size_t offset_exch;           /* offset of the quantity in the PrimExch struct */
  size_t offset_grad;           /* offset in the grad_data struct */
  double *min_value, *max_value;
  double value0, value1;
} grad_elements[MAXGRADIENTS], *GDensity, *GVelx, *GVely, *GVelz, *GPressure, *GUtherm
#ifdef MRT
  , grad_elements_RT[MAXRTGRADIENTS] 
#endif
  ;

extern int N_Grad;              /* number of gradients to be calculated */

#ifdef MRT
extern int N_Grad_RT;
#endif

extern struct grad_data
{
  MySingle drho[3];
#if defined(DEREFINE_GENTLY) || defined(CONDUCTION_SATURATION) || defined(SPECIAL_RELATIVITY) || defined(GENERAL_RELATIVITY) || defined(NON_LINEAR_SLOPE_LIMITERS) || defined(CALCULATE_QUANTITIES_IN_POSTPROCESS)
  MySingle dutherm[3];
#endif

#ifdef DVR_RENDER
  MySingle dDvrFields[DVR_NUM_FIELDS][3];
#endif
  MySingle dvel[3][3];
  MySingle dpress[3];

#if defined(MRT) && defined(MRT_LSF_GRADIENTS)
  MySingle dFN[MRT_BINS][3][3] ;
  MySingle dmodFN[MRT_BINS][3] ;
  MySingle dDensPhot[MRT_BINS][3] ;
#ifdef MRT_TIME_EXTRAPOLATION
  MySingle dPT[MRT_BINS][3][3][3] ;
#endif
#endif

#ifdef COSMIC_RAYS
  MySingle dcrPressure[3];
#endif

#ifdef NUCLEAR_NETWORK
  MySingle drhoU[3];
  MySingle dpressU[3];
#endif

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
  MySingle dA[3];
#endif
#ifdef TGCHEM
  MySingle dgamma[3];
#endif
#ifdef VARIABLE_GAMMA
  MySingle dgammaE[3];
  MySingle dgammaC[3];
#endif
#ifdef MHD
  MySingle dB[3][3];
#endif
#ifdef MHD_CT
	MySingle dA[3][3];
#endif
#ifdef MAXSCALARS
  MySingle dscalars[MAXSCALARS][3];
#endif
#ifdef TRACER_FIELD
  MySingle dtracer[3];
#endif

#ifdef FLD
	MyFloat dngamma[3];
#endif
#if defined(VORONOI_PROJ_TAU)
  MySingle dtemp[3];
#endif
}
 *GradExch;

#ifdef TVD_SLOPE_LIMITER
extern struct grad_data *GradExchUl;
#endif

#ifdef RT_ADVECT
extern struct rt_grad_data
{
  double ddensphot[RT_N_DIR][3];
#ifdef RT_HEALPIX_NSIDE
  double ddensphot_unlimited[3];
#else
  MyFloat sourcepos[RT_N_DIR][3];
  int sourceid;
#endif
}
 *RTGradExch;


extern struct rt_primexch
{
  double DensPhot[RT_N_DIR];
#ifndef RT_HEALPIX_NSIDE
  int SourceID[RT_N_DIR];
  MyFloat SourcePos[RT_N_DIR][3];
#endif
}
 *RTPrimExch;
#endif


#ifdef MRT
extern struct rt_grad_data
{
  MySingle dFN[MRT_BINS][3][3] ;
  MySingle dmodFN[MRT_BINS][3] ;
  MySingle dDensPhot[MRT_BINS][3] ;
#ifdef MRT_TIME_EXTRAPOLATION
  MySingle dPT[MRT_BINS][3][3][3] ;
#endif
}
  *RTGradExch;

extern struct rt_primexch
{
  double PT[MRT_BINS][3][3] ;
  double RT_F[MRT_BINS][3] ;
  double modFN[MRT_BINS] ;
  double FN[MRT_BINS][3] ;
  double DensPhot[MRT_BINS] ;
  double OldCons_DensPhot[MRT_BINS] ;
}
  *RTPrimExch ;
#endif

extern struct primexch
{
  double Volume;
  MyFloat Density;
#if defined(DEREFINE_GENTLY) || defined(CONDUCTION_SATURATION) || (defined(COSMIC_RAYS) && defined(SHOCK_FINDER_ON_THE_FLY) || defined(SHOCK_FINDER_BEFORE_OUTPUT)) || defined(SPECIAL_RELATIVITY) || defined(GENERAL_RELATIVITY) || defined(NON_LINEAR_SLOPE_LIMITERS) || defined(CALCULATE_QUANTITIES_IN_POSTPROCESS)
  MyFloat Utherm;
#endif

#ifdef DVR_RENDER
  MyFloat DvrFields[DVR_NUM_FIELDS];
#endif
  MyFloat VelGas[3];
  MyFloat VelVertex[3];
#ifdef TGCHEM
  double Abund[TGCHEM_NUM_ABUNDANCES];
  double Gamma;
#endif
#ifdef VARIABLE_GAMMA
  MyFloat GammaE;
  MyFloat GammaC;
#endif
#ifdef MHD
  MyFloat B[3];
#ifdef MHD_POWELL
  MyFloat DivB;
#endif
#ifdef MHD_CT
  MyFloat A[3];
  double TimeLastBUpdate;
#endif
  MyFloat CurlB[3];
#endif
  MyFloat Pressure;
#ifdef USE_ENTROPY_FOR_COLD_FLOWS
  MyFloat A;
#endif
#if defined(USE_ENTROPY_FOR_COLD_FLOWS) || (defined(TRACER_MC) && defined(TWODIMS)) || defined(COSMIC_RAYS_SHOCK_ACCELERATION)
  MyFloat Mass;
#endif
#ifdef MAXSCALARS
  MyFloat Scalars[MAXSCALARS];
#endif
#ifdef GFM_CHEMTAGS
  MyFloat MassMetalsChemTagsFraction[GFM_N_CHEM_TAGS];
#endif
#ifdef TRACER_FIELD
  MyFloat Tracer;
#ifdef SECOND_DERIVATIVES
  struct grad_data Grad;
#endif
#endif
#if defined(TRACER_MC) && defined(TWODIMS)
  int tracerMC_num;
#endif
#ifdef DG
	MyDouble Weights[NOF_BASE_FUNCTIONS][5];
#endif
#ifdef COSMIC_RAYS
  double CR_Pressure;
  double CR_SpecificEnergy;
#ifdef COSMIC_RAYS_STREAMING
  double CR_Chi;
#endif
#endif

#if defined(MRT) && defined(MRT_LSF_GRADIENTS)
  double PT[MRT_BINS][3][3] ;
  double RT_F[MRT_BINS][3] ;
  double modFN[MRT_BINS] ;
  double FN[MRT_BINS][3] ;
  double DensPhot[MRT_BINS] ;
  double OldCons_DensPhot[MRT_BINS] ;
  double HI, HII, Ne, HeI, HeII, HeIII ;

#ifdef MRT_IR_PHOTON_TRAPPING
  double Trapped_DensPhot[UV_BINS] ;
#endif
  int flag ;
#endif

#ifdef OUTPUT_CELL_SPIN
  double CenterOffset[3];
#endif
#ifdef ACTIVE_CELL_SPIN
  double Omega[3];
#endif

#ifdef FLD
      MyFloat n_gamma;
      MyFloat Kappa_diff;
      MyFloat R2;
#ifdef FLD_CONES
      MyFloat gammas[FLD_NCONES];
#endif
#endif

  double TimeLastPrimUpdate;

  MyDouble Center[3];
  MyFloat OldMass;
  MySingle Csnd;
  MySingle SurfaceArea;
  MySingle ActiveArea;
  /*  int task, index; */
  short int TimeBinHydro;
#if (defined(SHOCK_FINDER_POST_PROCESSING) || defined(SHOCK_FINDER_ON_THE_FLY) || defined(SHOCK_FINDER_BEFORE_OUTPUT)) && defined(USE_SFR)
  MySingle Sfr;
#endif
#ifdef SHOCK_FINDER_ON_THE_FLY
  MyFloat Divvel;
  MyFloat ShockDir[3];
  char ShockZone;
#ifdef COSMIC_RAYS
  MyFloat CRpseudoTemperature;
#else
  MyFloat Temperature;
#endif
#endif
#if defined(VORONOI_PROJ_TAU) && !defined(SHOCK_FINDER_ON_THE_FLY)
  MyFloat Temperature;
#endif
}
 *PrimExch;

#ifdef REFINEMENT
extern struct refdata
{
#ifdef REFINEMENT_VOLUME_LIMIT
  double Volume;
#endif
  short int TimeBinHydro;
}
 *RefExch;
#endif

#ifdef SECOND_DERIVATIVES
extern struct hessian_data
{
  double ddrho[3][3];
  double ddvelx[3][3];
  double ddvely[3][3];
  double ddvelz[3][3];
  double ddpress[3][3];
#if defined(TRACER_FIELD) && defined(TRACER_DIFFUSION)
  double ddtracer[3][3];
#endif
}
 *HessianExch;
#endif

typedef struct face_data
{
  int p1, p2;
#ifdef REFINEMENT_MERGE_CELLS
  int t, nr;                    /* delaunay tetra and edge number that generated this face */
#endif
#ifdef OPTIMIZE_MEMORY_USAGE
  MyFloat area;
  MyFloat cx, cy, cz;           /* center-of-mass of face */
#else
  double area;
  double cx, cy, cz;            /* center-of-mass of face */
#endif
#ifdef VORONOI_BACKUP_RESTORE_FACE_AREAS
  double area_backup;
#endif
#ifdef TETRA_INDEX_IN_FACE
  int dt_index;
#endif


#ifdef MONOTONE_CONDUCTION
#ifndef SEMI_IMPLICIT_TI
  double Flux1p1, Flux2p1 ;
  double Flux1p2, Flux2p2 ;
#endif
#endif
#ifdef SEMI_IMPLICIT_TI
  double akssp1, akssp2 ;
  int kidp1, kidp2 ;
  int ktaskp1, ktaskp2 ;
  double F2p1, F2p2 ;
  double akssp1p2, akssp2p1 ;
  double F2p1p2, F2p2p1 ;
  int originalindexp1, originalindexp2 ;
  int doface ;
  double alpha_fac ;
#endif
#ifdef IMPLICIT_OHMIC_DIFFUSION
  double akssp1, akssp2 ;
  int kidp1, kidp2 ;
  int ktaskp1, ktaskp2 ;
  double F2p1, F2p2 ;
  double akssp1p2, akssp2p1 ;
  double F2p1p2, F2p2p1 ;
  int originalindexp1, originalindexp2 ;
  int doface ;
#endif

}
face;


/*!< left or right state of a face */
struct state
{
  double dx, dy, dz;
  double dt_half;
  short int timeBin;

  double rho;
#if defined(SPECIAL_RELATIVITY)  || defined(GENERAL_RELATIVITY)
  double utherm;
#endif
  double velx, vely, velz;
  double press;
  double oldmass;
  double surfacearea;
  double activearea;
  double volume;

  MyFloat velGas[3];
  MyFloat velVertex[3];
  struct grad_data *grad;
#ifdef TVD_SLOPE_LIMITER
  struct grad_data *gradul;
#endif
#ifdef MRT
  struct rt_grad_data *rtgrad;
#endif

  double csnd;
#ifdef TGCHEM
  double pcabund[TGCHEM_NUM_ABUNDANCES];
  double gamma;
#endif
#ifdef VARIABLE_GAMMA
  double gammaE;
  double gammaC;
#endif
  double Energy;
#ifdef MHD
  double Bx, By, Bz;
#ifdef MHD_CT
	double Ax, Ay, Az;
#endif  
#ifdef MHD_POWELL
  double divB;
#endif
  double CurlB[3];
#ifdef NON_IDEAL_MHD
#if defined(OHMIC_DIFFUSION)
  MyFloat dcp[3];
#endif
#endif
#endif
#ifdef COSMIC_RAYS
  double crPressure;
#endif
  
#ifdef MRT_LSF_GRADIENTS
  double DensPhot[MRT_BINS] ;
  double RT_F[MRT_BINS][3] ;
  double FN[MRT_BINS][3] ;
  double modFN[MRT_BINS] ;
  double PT[MRT_BINS][3][3] ;
  double OldCons_DensPhot[MRT_BINS] ;
  int flag ;
  double nHI, nHII, ne, nHeI, nHeII, nHeIII ;

#ifdef MRT_TIME_EXTRAPOLATION
  double divF[MRT_BINS] ;
#endif
#ifdef MRT_IR_PHOTON_TRAPPING
  double Trapped_DensPhot[UV_BINS] ;
#endif
#endif
  
#ifdef USE_ENTROPY_FOR_COLD_FLOWS
  double A;
#endif
#if defined(ENTROPY_MACH_THRESHOLD) || defined(GODUNOV_STATS) || (TRACER_MC_MACHMAX) || (TRACER_PART_MACHMAX)
  double mach;
#endif
#ifdef MAXSCALARS
  double scalars[MAXSCALARS];
#endif
#ifdef GFM_CHEMTAGS
  double chemtagsfraction[GFM_N_CHEM_TAGS];
#endif
#if defined(TRACER_FIELD) || !defined(MESHRELAX)
  double tracer;
#endif
  MyIDType ID;
#ifdef SECOND_DERIVATIVES
  struct hessian_data *hessian;
#endif
#ifdef LOCALLY_ISOTHERM_DISK
  double localSoundSpeed;
  int inDiskFlag;
#endif
#ifdef ACTIVE_CELL_SPIN
  double dxL, dyL, dzL;
  double wx, wy, wz;
  double dvx, dvy, dvz;
#endif
#ifdef ONEDIMS_SPHERICAL
  double radius;
#endif

  double dtExtrapolation;

#ifdef MHD_CT
  double dtBExtrapolation;
#endif

#ifdef FLD
  double n_gamma;
  double R2;
#endif
};

/*!< state on a face determined by riemann solver */
extern struct state_face
{
  double rho;
#if defined(SPECIAL_RELATIVITY) || defined(GENERAL_RELATIVITY)
  double utherm;
#endif
  double velx, vely, velz;
  double press;
#ifdef MHD
  double Bx, By, Bz;
#endif
#ifdef MHD_CT
	double Ax, Ay, Az;
#endif
#ifdef TGCHEM
  double pcabund[TGCHEM_NUM_ABUNDANCES];
  double gamma;
#endif
#ifdef VARIABLE_GAMMA
  double gammaE;
#endif
#ifdef USE_ENTROPY_FOR_COLD_FLOWS
  double A;
#endif
#ifdef MAXSCALARS
  double *scalars;
#endif
#if defined(TRACER_FIELD) || !defined(MESHRELAX)
  double tracer;
#endif
#ifdef TRACER_DIFFUSION
  double dConservedTracer[3];
#endif
#ifdef THERMAL_CONDUCTION
  double dTemp[3];
#endif
#ifdef VISCOSITY
  double vel_grad[3][3];        /*velocity gradients at the interface */
#endif
#ifdef LOCALLY_ISOTHERM_DISK
  double localSoundSpeed;
#endif
#ifdef FLD
      double n_gamma;
      double R2;  
#endif
}
state_face;

/*!< flux through a face */
extern struct fluxes
{
  double mass;
  double momentum[3];
  double energy;

#ifdef MRT_LSF_GRADIENTS
  double DensPhot[MRT_BINS] ;
  double RT_F[MRT_BINS][3] ;
  double nHI, nHII, ne, nHeI, nHeII, nHeIII ;
#ifdef MRT_IR_PHOTON_TRAPPING
  double Trapped_DensPhot[IR_BINS] ;
#endif
#endif

#ifdef MHD
  double B[3];
#endif
#ifdef MHD_CT
	double A[3];
#endif
#ifdef USE_ENTROPY_FOR_COLD_FLOWS
  double entropy;
#endif
#ifdef MAXSCALARS
  double scalars[MAXSCALARS];
#endif
#if defined(TRACER_FIELD) || !defined(MESHRELAX)
  double tracer;
#endif
#ifdef TGCHEM
  double pcabund[TGCHEM_NUM_ABUNDANCES];
#endif
#ifdef FLD
  double dFLD;
#endif
#ifdef GFM_CHEMTAGS
  double chemtags[GFM_N_CHEM_TAGS];
#endif
}
fluxes, diffusionfluxes;


extern struct geometry
{
  double nn;
  double nx, ny, nz;
  double mx, my, mz;
  double px, py, pz;
  double cx, cy, cz;
}
geom;

struct pv_update_data
{
  double atime;
  double hubble_a;
  double a3inv;

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
  int count_keep_entropy, count_update_entropy;
#endif
};
#endif //MESH_H

struct fvs_stat
{
  int count_disable_extrapolation;
#ifdef COSMIC_RAYS
  int count_CR_limiter;
#endif
};
