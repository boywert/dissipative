/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/tgchem/tgchem.h
 * \date        01/2013
 * \author      Thomas Greif
 * \brief       Primordial chemistry and cooling network 
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 */

#define ELECTRON_VOLT 1.60219e-12
#define TEMP_CMB 2.725
#define TGCHEM_CHI_HI (13.6 * ELECTRON_VOLT)
#define TGCHEM_CHI_H2 (4.48 * ELECTRON_VOLT)
#define TGCHEM_CHI_HM (3.53 * ELECTRON_VOLT)

#define TGCHEM_NUM_ABUNDANCES 3
#define TGCHEM_NUM_SPECIES (TGCHEM_NUM_ABUNDANCES + 1)
#define TGCHEM_NUM_EQ 2
#define TGCHEM_NUM_CHEM 7
#define TGCHEM_NUM_COOL 4
#define TGCHEM_NUM_RATES 6
#define TGCHEM_NUM_COOLING 5
#define TGCHEM_NUM_PHOTO 1
#define TGCHEM_NUM_TEMP 5000

#define TGCHEM_TEMP_MIN 1e1
#define TGCHEM_TEMP_MAX 1e8
#define TGCHEM_LOG_DTEMP (log10(TGCHEM_TEMP_MAX / TGCHEM_TEMP_MIN) / TGCHEM_NUM_TEMP)

#define TGCHEM_TOL 1e-4
#define TGCHEM_TOL_EQ 1e-7
#define TGCHEM_TOL_STEPSIZE 0.2
#define TGCHEM_TOL_ABHM 1e-20
#define TGCHEM_TOL_ABH2 1e-20
#define TGCHEM_TOL_ABHII 1e-20
#define TGCHEM_TOL_ABENERGY 0.

#define TGCHEM_MAX_RATE 1e100

#define TGCHEM_H2_ENERGY_SA 8
#define TGCHEM_H2_ENERGY_SB 5
#define TGCHEM_H2_NUM_V 3
#define TGCHEM_H2_NUM_J 20
#define TGCHEM_H2_NUM_T 3
#define TGCHEM_H2_TOT_NUM_LINES (TGCHEM_H2_NUM_V * TGCHEM_H2_NUM_J * TGCHEM_H2_NUM_V * TGCHEM_H2_NUM_T)
#define TGCHEM_H2_NUM_COLUMN 200
#define TGCHEM_H2_COLUMN_MIN 1e21
#define TGCHEM_H2_COLUMN_MAX 1e30
#define TGCHEM_H2_LOG_DCOLUMN (log10(TGCHEM_H2_COLUMN_MAX / TGCHEM_H2_COLUMN_MIN) / TGCHEM_H2_NUM_COLUMN)

#define TGCHEM_OPAC_NUM_RHO 15
#define TGCHEM_OPAC_NUM_TEMP 29

#define TGCHEM_TEMP() {TGCD.TempIndex = imin((int) (log10(TGCD.Temp / TGCHEM_TEMP_MIN) / TGCHEM_LOG_DTEMP), TGCHEM_NUM_TEMP - 1); TGCD.DTemp = TGCD.Temp - TGCD.TempTable[TGCD.TempIndex];}
#define TGCHEM_EQRATE(x) {TGCD.EqRate[x] = TGCD.EqTable[x * TGCHEM_NUM_TEMP + TGCD.TempIndex] + TGCD.DTemp * TGCD.DEqTable[x * TGCHEM_NUM_TEMP + TGCD.TempIndex];}
#define TGCHEM_CHEMRATE(x) {TGCD.ChemRate[x] = TGCD.ChemTable[x * TGCHEM_NUM_TEMP + TGCD.TempIndex] + TGCD.DTemp * TGCD.DChemTable[x * TGCHEM_NUM_TEMP + TGCD.TempIndex];}
#define TGCHEM_COOLRATE(x) {TGCD.CoolRate[x] = TGCD.CoolTable[x * TGCHEM_NUM_TEMP + TGCD.TempIndex] + TGCD.DTemp * TGCD.DCoolTable[x * TGCHEM_NUM_TEMP + TGCD.TempIndex];}

#define TGCHEM_H2_COLUMN() {TGCD.H2ColumnIndex = imin((int) (log10(TGCD.H2Column / TGCHEM_H2_COLUMN_MIN) / TGCHEM_H2_LOG_DCOLUMN), TGCHEM_H2_NUM_COLUMN - 1); TGCD.H2DColumn = TGCD.H2Column - TGCD.H2ColumnTable[TGCD.H2ColumnIndex];}
#define TGCHEM_H2_SOBEMISS() {int idx1 = TGCD.TempIndex * TGCHEM_H2_NUM_COLUMN + TGCD.H2ColumnIndex; int idx2 = dmin((TGCD.TempIndex + 1), TGCHEM_NUM_TEMP - 1) * TGCHEM_H2_NUM_COLUMN + TGCD.H2ColumnIndex; double emiss1 = TGCD.H2SobEmissTable[idx1] + TGCD.H2DColumn * TGCD.H2DSobEmissTable[idx1]; double emiss2 = TGCD.H2SobEmissTable[idx2] + TGCD.H2DColumn * TGCD.H2DSobEmissTable[idx2]; TGCD.H2SobEmiss = emiss1 + TGCD.DTemp * (emiss2 - emiss1);}

extern struct TGCD_struct
{
  // Input parameters
  int ChemMode;
  int ChemIOMode;
  int ChemRndMode;
  int ChemH2Mode;
  double ChemH2Thresh;
  double ChemHIITrans;
  double ChemHIIThresh;
  double ChemInitAbH2;
  double ChemInitAbHII;
  double ChemJ21;

  // Global parameters
  int DebugFlag;
  double RedShift;
  double TempCMB;
  double ConvertEnergy;

  // Global test parameters
  double TestStepSize;
  double CollapseFac;
  double NHStart;
  double NHEnd;

  // Constants
  double H2OptThickConst;
  double H2OptThickNHThresh;
  double CIEOptThickNHThresh;
  double LWDissRate;
  double HIPhotonEnergy;

  // Step Constants
  int EqH2Flag;
  int EqHIIFlag;

  // Logging
  int NumEq;
  int NumNeq;
  int NumRateCalls;
  int MaxNumRateCalls;
  long TotNumRateCalls;
  int NumSubSteps;
  int MaxNumSubSteps;
  long TotNumSubSteps;
  long NumCellsDone;
  double DtEq;
  double DtNEq;

  // Variables
  int Task;
  int Index;
  int ID;
  int TempIndex;
  int H2ColumnIndex;
  double IntDensity;
  double Density;
  double NH;
  double TransPower;
  double DivVel;
  double Temp;
  double DTemp;
  double H2Column;
  double H2DColumn;
  double H2SobEmiss;
  double H2EscFrac;
  double NTot;
  double Mu;
  double Energy;
  double Gamma;
  double Csnd;
  double VThermal;
  double JeansLength;

  // CVODE
  void *CVODEMem;
  N_Vector Species;
  N_Vector SpeciesSave;
  N_Vector SpeciesStepSave;
  N_Vector DSpecies;
  N_Vector SpeciesTol;

  // Abundances
  double AbHM;
  double AbH2;
  double AbHI;
  double AbHII;
  double AbE;
  double AbMax[TGCHEM_NUM_ABUNDANCES];

  // Tables
  double *TempTable;
  double *EqTable;
  double *DEqTable;
  double *ChemTable;
  double *DChemTable;
  double *CoolTable;
  double *DCoolTable;

  // Rates
  double EqRate[TGCHEM_NUM_EQ];
  double ChemRate[TGCHEM_NUM_CHEM];
  double CoolRate[TGCHEM_NUM_COOL];
  double CrRate[TGCHEM_NUM_ABUNDANCES][TGCHEM_NUM_RATES];
  double DsRate[TGCHEM_NUM_ABUNDANCES][TGCHEM_NUM_RATES];
  double ChemicalRate[TGCHEM_NUM_ABUNDANCES][TGCHEM_NUM_RATES];
  double TotCrRate[TGCHEM_NUM_ABUNDANCES];
  double TotDsRate[TGCHEM_NUM_ABUNDANCES];
  double TotChemicalRate[TGCHEM_NUM_ABUNDANCES];
  double CoolingRate[TGCHEM_NUM_COOLING];
  double PhRate[TGCHEM_NUM_PHOTO];
  double ChemicalRateSum;
  double CoolingRateSum;
  double HydroHeatRate;
  double HRHeatRate;
  double HeatRate;

  // H2 data
  double H2Energy[TGCHEM_H2_NUM_V * TGCHEM_H2_NUM_J];
  double *H2SpontCoeff;
  double *H2ColumnTable;
  double *H2SobEmissTable;
  double *H2DSobEmissTable;
  double *H2TestSobXVal;
  double *H2TestFitEscFrac;
  double *H2TestSobEscFrac;

  double *H2EmissTable;
  double *H2DEmissTable;
  double *H2TotEmissTable;
  double *H2DTotEmissTable;
  double *H2CrossTable;
  double *H2DCrossTable;
  double *H2TotCrossTable;
  double *H2DTotCrossTable;
#ifdef HEALRAY
  int *H2LineIdxTable;
#endif

  // Opacity
  double Opac;
  double *OpacTable;

} TGCD;
