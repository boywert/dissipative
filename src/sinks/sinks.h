/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/sinks/sinks.h
 * \date        01/2013
 * \author      Thomas Greif
 * \brief       Sink particles
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 */

#define SKD_INIT_MAX_NUM_ACC 10


struct ACC_struct
{
  int Task;
  int Index;
  int Active;
  double Pos[3];
  double Vel[3];
  double Mass;
#ifdef SGCHEM
  double Utherm;
  double Accel[3];
#endif
};

struct SINK_struct
{
  double Pos[3];
  double Mass;
};

extern struct SKD_struct
{
  int Task;
  int Index;
  int Flag;
  int NumAcc;
  int MaxNumAcc;
  int NumSinks;
  int TotNumSinks;
  MyIDType ID;
  double NHThresh;
  double AccRad;
  double AccRad2;
  double DistFac;
  double NHFac;
  double NHMax;
  double Pos[3];
  double Mass;
  double CoMPos[3];
  double CoMVel[3];

  struct ACC_struct *ACC;
  struct SINK_struct *SINK;

} SKD;

extern struct SinksAux_struct
{
  int SinksAuxID;
  signed char MinTimeBin;
  double dMass;
  double MassNorm;
} *SinksAux;

extern int NSinks, TotNSinks;
