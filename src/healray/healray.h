/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/healray/healray.h
 * \date        01/2013
 * \author      Thomas Greif
 * \brief       Adaptive ray-tracing
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 */

#include "chealpix.h"

#define HRD_MEM_SAFETY_FAC 10.
#define HRD_RECOMBINATION_COEFFICIENT 2.6e-13

struct HRSD_struct
{
  int idx;
  int task;
  double pos[3];
};

struct HRRD_struct
{
  int iter;
  int task;
  int task_init;
  int task_orig;
  int idx;
  int idx_prev;
  int idx_init;
  int idx_orig;
  int lvl;
  long hpidx;
  long hpidx_orig;
  int num_advances;
  int num_comms;
  long id;
  double len;
  double pos[3];
  double pos_old[3];
  double pos_init[3];
  double rot[9];
  double dir[3];
  double energy;
  double energy_init;
  double energy_init_tot;
#ifdef TGCHEM
  int temp_idx_orig;
  //double pos_init[3];
  double vel_init[3];
  double temp_init;
#endif
};

struct EscFrac_struct
{
  int task;
  int idx;
  double energy_frac;
};

extern struct HRSL_struct
{
  int id;
  double pos[3];
  double hsml;

} *HRSL;

extern struct HRD_struct
{
  // Parameters
  char RaySourceFile[MAXLEN_PATH];
  int RayMode;
  int RayIOMode;
  int RayNumThreads;
  int RayDebugMode;
  int RayTimeFac;
  int RayMultiFreqMode;
  int RayNumSources;
  int RayLvlInit;
  int RayNumLines;
  int RayNumFreq;
  double RaySplitFac;
  double RayPartFrac;
  double RayMinEgyFrac;

  // Constants
  int SourceFlag;
  int TraceNumIter;
  int RayNumBins;
  int NumThreads;
  int NumRaysInit;
  double EnergyUnit;
  double UThermToTemp;
  double RayTemp;
  double SourceLuminosity;
#ifdef TGCHEM
  double HIPhotonEnergy;
#endif
  double LenToCgs;
  double VolToCgs;
  double RhoToNHCgs;
  double VelToCgs;
  double EnergyIncreaseFac;
  double NuDRange;
  double NuDWidth;
  double *NuDPos;
  double *NuDFac;

  // Counters
  int Loop;
  int Iter;

  // Sanity check
  double TotEnergyEsc;
  double TotEnergyInit;

  // Timing
  double DtInit;
  double DtAdvance;
  double DtImba;
  double DtComm;

  // Logging
  long LocNumOps;
  long LocNumAdvances;
  long LocNumComms;
  long NumOps;
  long NumAdvances;
  long NumComms;
  double AdvancesPerRay;
  double CommsPerRay;

  // Test run
  double InitEnergy;
  double Alpha;
  double Error;
  int *SourceIDList;
  double *SourcePosList;
  double *SourceEgyList;

  // Ray output
  int RayNumTrace;
  long *TraceList;
  void **TraceEnergyFrac;
  void **TraceEnergyFracFreq;
  FILE *TraceFile;

  // Mesh;
  short int *SaveTimeBin;
  int SaveTimeBinActive[TIMEBINS];

  // Sources
  int NumSources;
  int TotNumSources;
  struct HRSD_struct *HRSD;

  // Rays
  int NumRays;
  int MaxNumRays;
  int LimitNumRays;
  int TmpNumRays;
  int TmpMaxNumRays;
  long long TotInitNumRays;
  int *ThreadNumRays;
  int *ThreadPrevNumRays;
  struct HRRD_struct *HRRD;
  struct HRRD_struct *TmpHRRD;
  double *RayEnergyLine;
  double *TmpRayEnergyLine;
  double *RayEnergyLineInit;
  double *TmpRayEnergyLineInit;
  float *RayEnergy;
  float *TmpRayEnergy;

  // Gas
  int PIdx;
  int MaxNumEmitters;
  int TotNumEmitters;
  int NumEscFrac;
  int MaxNumEscFrac;
  int *IdxList;
  int *PartRayCount;
  int *ThreadNumEscFrac;
  double *EnergyInit;
  double *Energy;
  double *EnergyFrac;
  double *EnergyAux;
  struct EscFrac_struct *EscFrac;

} HRD;
