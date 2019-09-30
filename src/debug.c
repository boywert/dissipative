/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/debug.c
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

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"

/** \file debug.c
 *  \brief Print relevant information about a particle / face to standard output for debug. 
 *
 *  The functions contained in this file are mostly called when a condition, that 
 *  causes the abort of the run, is met. In that case, the information about the state of 
 *  the particle / face which triggered that condition is printed to the standard output. 
 */


/** \brief Print information relative to a particle / cell to standard output.
 *
 *  \param i particle / cell index
 */
void print_particle_info(int i)
{
  printf("Task=%d, ID=%llu, Type=%d, TimeBinGrav=%d, TimeBinHydro=%d, Mass=%g, pos=%g|%g|%g, vel=%g|%g|%g\n", ThisTask,
         (unsigned long long) P[i].ID, P[i].Type, P[i].TimeBinGrav, P[i].TimeBinHydro, P[i].Mass, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], P[i].Vel[0], P[i].Vel[1], P[i].Vel[2]);
#ifdef PMGRID
  printf("GravAccel=%g|%g|%g, GravPM=%g|%g|%g, Soft=%g, SoftType=%d, OldAcc=%g\n", P[i].GravAccel[0], P[i].GravAccel[1], P[i].GravAccel[2], P[i].GravPM[0],
         P[i].GravPM[1], P[i].GravPM[2], All.ForceSoftening[P[i].SofteningType], P[i].SofteningType, P[i].OldAcc);
#else
  printf("GravAccel=%g|%g|%g, Soft=%g, SoftType=%d, OldAcc=%g\n", P[i].GravAccel[0], P[i].GravAccel[1], P[i].GravAccel[2], All.ForceSoftening[P[i].SofteningType], P[i].SofteningType, P[i].OldAcc);
#endif
  if(P[i].Type == 0)
    {
#ifndef COSMIC_RAYS
      printf("Vol=%g, rad=%g, rho=%g, p=%g,u=%g, velVertex=%g|%g|%g, csnd=%g\n",
             SphP[i].Volume, get_cell_radius(i), SphP[i].Density, SphP[i].Pressure, SphP[i].Utherm, SphP[i].VelVertex[0], SphP[i].VelVertex[1], SphP[i].VelVertex[2], get_sound_speed(i));
#else
      printf("Vol=%g, rad=%g, rho=%g, p=%g, p_CR=%g, u=%g, velVertex=%g|%g|%g, csnd=%g\n",
             SphP[i].Volume, get_cell_radius(i), SphP[i].Density, SphP[i].Pressure, SphP[i].CR_Pressure, SphP[i].Utherm,
             SphP[i].VelVertex[0], SphP[i].VelVertex[1], SphP[i].VelVertex[2], get_sound_speed(i));
#endif
      printf("Center-Pos=%g|%g|%g\n", SphP[i].Center[0] - P[i].Pos[0], SphP[i].Center[1] - P[i].Pos[1], SphP[i].Center[2] - P[i].Pos[2]);
#ifndef MHD
      printf("Mom=%g|%g|%g, Energy=%g, EInt=%g, EKin=%g\n", SphP[i].Momentum[0], SphP[i].Momentum[1],
             SphP[i].Momentum[2], SphP[i].Energy, SphP[i].Utherm * P[i].Mass,
             0.5 * P[i].Mass * ((SphP[i].Momentum[0] / P[i].Mass) * (SphP[i].Momentum[0] / P[i].Mass) +
                                (SphP[i].Momentum[1] / P[i].Mass) * (SphP[i].Momentum[1] / P[i].Mass) + (SphP[i].Momentum[2] / P[i].Mass) * (SphP[i].Momentum[2] / P[i].Mass)));
#else
#ifndef COSMIC_RAYS
      printf("Mom=%g|%g|%g, Energy=%g, EInt=%g, EKin=%g, EB=%g\n", SphP[i].Momentum[0], SphP[i].Momentum[1],
             SphP[i].Momentum[2], SphP[i].Energy, SphP[i].Utherm * P[i].Mass,
             0.5 * P[i].Mass * ((SphP[i].Momentum[0] / P[i].Mass) * (SphP[i].Momentum[0] / P[i].Mass) +
                                (SphP[i].Momentum[1] / P[i].Mass) * (SphP[i].Momentum[1] / P[i].Mass) +
                                (SphP[i].Momentum[2] / P[i].Mass) * (SphP[i].Momentum[2] / P[i].Mass)),
             0.5 * SphP[i].Volume * (SphP[i].B[0] * SphP[i].B[0] + SphP[i].B[1] * SphP[i].B[1] + SphP[i].B[2] * SphP[i].B[2]));
#else
      printf("Mom=%g|%g|%g, Energy=%g, CR_Energy=%g, EInt=%g, EKin=%g, EB=%g\n", SphP[i].Momentum[0], SphP[i].Momentum[1],
             SphP[i].Momentum[2], SphP[i].Energy, SphP[i].CR_Energy, SphP[i].Utherm * P[i].Mass,
             0.5 * P[i].Mass * ((SphP[i].Momentum[0] / P[i].Mass) * (SphP[i].Momentum[0] / P[i].Mass) +
                                (SphP[i].Momentum[1] / P[i].Mass) * (SphP[i].Momentum[1] / P[i].Mass) +
                                (SphP[i].Momentum[2] / P[i].Mass) * (SphP[i].Momentum[2] / P[i].Mass)),
             0.5 * SphP[i].Volume * (SphP[i].B[0] * SphP[i].B[0] + SphP[i].B[1] * SphP[i].B[1] + SphP[i].B[2] * SphP[i].B[2]));
#endif
#endif
#ifdef GFM_COOLING_METAL
      printf("H Mass Frac=%g|He Mass Frac=%g|Z Mass Frac=%g\n", SphP[i].MetalsFraction[element_index("Hydrogen")], SphP[i].MetalsFraction[element_index("Helium")], SphP[i].Metallicity);
      printf("H Mass =%g|He Mass=%g|Z Mass=%g\n", SphP[i].MassMetals[element_index("Hydrogen")], SphP[i].MassMetals[element_index("Helium")], SphP[i].MassMetallicity);
#endif
#if defined(EOS_DEGENERATE) || defined(EOS_OPAL)
      printf("gammaC=%g, gammaE=%g, temp=%g\n", SphP[i].GammaC, SphP[i].GammaE, SphP[i].EOSTemperature);
#endif
#ifdef MHD
      double err = pow(SphP[i].Volume,
                       1. / 3.) * fabs(SphP[i].DivB) / sqrt(SphP[i].B[0] * SphP[i].B[0] + SphP[i].B[1] * SphP[i].B[1] + SphP[i].B[2] * SphP[i].B[2]);
      printf("B=%g|%g|%g, divb=%g, err=%g\n", SphP[i].B[0], SphP[i].B[1], SphP[i].B[2], SphP[i].DivB, err);
#endif

#ifdef MRT
      printf("N=%g|Fx=%g|Fy=%g|Fz=%g\n", SphP[i].DensPhot[0], SphP[i].RT_F[0][0], SphP[i].RT_F[0][1], SphP[i].RT_F[0][2]) ;
#ifdef MRT_IR
      printf("Px=%g|Py=%g|Pz=%g|Kappa_P=%g|Kappa_R=%g\n", SphP[i].Momentum[0], SphP[i].Momentum[0], SphP[i].Momentum[0], SphP[i].KappaIR_P[0], SphP[i].KappaIR_R[0]) ;
#endif
#endif

#ifdef OUTPUT_CELL_SPIN
      printf("Spin=%g|%g|%g.\n", SphP[i].Spin[0], SphP[i].Spin[1], SphP[i].Spin[2]);
#ifdef ACTIVE_CELL_SPIN
      printf("Omega=%g|%g|%g.\n", SphP[i].Omega[0], SphP[i].Omega[1], SphP[i].Omega[2]);
      int k;
      for(k = 0; k < 3; k++)
        printf("MomentOfInertia: %6g %6g %6g\n", SphP[i].MomentOfInertia[k][0], SphP[i].MomentOfInertia[k][1], SphP[i].MomentOfInertia[k][2]);
      printf("MomentOfInertiaSphere: %g\n", 2. / 5. * P[i].Mass * get_cell_radius(i) * get_cell_radius(i) / SphP[i].Density);
#endif
#endif

#ifdef TREE_BASED_TIMESTEPS
      printf("ID=%llu SphP[p].CurrentMaxTiStep=%g\n", (unsigned long long) P[i].ID, SphP[i].CurrentMaxTiStep);
#endif
    }

#ifdef GFM_STELLAR_EVOLUTION
  if(P[i].Type == 4 && StarP[P[i].AuxDataID].BirthTime > 0)     /* it's a star */
    {
      int k;
      double age = get_time_difference_in_Gyr(StarP[P[i].AuxDataID].BirthTime, All.Time);

      printf("AuxDataID=%llu, InitialMass=%g, Metallicity=%g,", (unsigned long long) P[i].AuxDataID, StarP[P[i].AuxDataID].InitialMass, StarP[P[i].AuxDataID].Metallicity);
      printf(" BirthTime=%g, Age=%g, SNIaRate=%g, SNIIRate=%g\n", StarP[P[i].AuxDataID].BirthTime, age, StarP[P[i].AuxDataID].SNIaRate, StarP[P[i].AuxDataID].SNIIRate);

      for(k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
        printf("%s mass=%g, mass fraction=%g\n", ElementNames[k], StarP[P[i].AuxDataID].MassMetals[k], StarP[P[i].AuxDataID].MassMetals[k] / P[i].Mass);
    }
#endif
#ifdef BLACK_HOLES
  if(P[i].Type == 5)            /* it's a black hole */
    {
      printf("SwallowID=%llu, AuxDataID=%llu, NumBHs=%d, Mass=%g, Mdot=%g\n", (unsigned long long) BPP(i).SwallowID, (unsigned long long) P[i].AuxDataID, NumBHs, BPP(i).BH_Mass, BPP(i).BH_Mdot);
#ifdef DRAINGAS
      printf("NearestDist=%g DrainID=%llu rho=%g hsml=%g NumNgb=%g\n", BPP(i).NearestDist, (unsigned long long) BPP(i).DrainID, BPP(i).BH_Density, BPP(i).BH_Hsml, BPP(i).BH_NumNgb);
#endif
#if (DRAINGAS==3)
      printf("VolSum=%g\n", BPP(i).BH_VolSum);
#endif
    }
#endif
#ifdef DUST_LIVE
  if(P[i].Type == DUST_LIVE)
    {
      printf("hsml=%g, GasVel=%g|%g|%g, csnd=%g, DragAccel=%g|%g|%g, t_s=%g\n", DTP(i).Hsml, DTP(i).LocalGasVelocity[0], DTP(i).LocalGasVelocity[1], DTP(i).LocalGasVelocity[2], DTP(i).LocalSoundSpeed, DTP(i).DragAccel[0], DTP(i).DragAccel[1], DTP(i).DragAccel[2], DTP(i).StoppingTime);
#ifdef DL_GRAIN_BINS
      for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
        printf("k=%d metal_frac=%g\n", k, DTP(i).MetalFractions[k]);

      for(int k = 0; k < DL_GRAIN_BINS; k++)
        {
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
          printf("bin k=%d NumGrains=%g BinSlopes=%g\n", k, DTP(i).NumGrains[k], DTP(i).BinSlopes[k]);
#else
          printf("bin k=%d NumGrains=%g\n", k, DTP(i).NumGrains[k]);
#endif
        }

#if defined(DL_SNE_DESTRUCTION) || defined(DL_SHATTERING) || defined(DL_COAGULATION)
      printf("dust_hsml=%g, DustDensity=%g, LocalCloudFrac=%g\n", DTP(i).DustHsml, DTP(i).DustDensity, DTP(i).LocalCloudFrac);
#endif
#endif
    }
#endif
}

/** \brief Print information relative to a particle / cell to standard output given its ID.
 *  *
 *   *  \param ID particle / cell ID
 *    */
void print_particle_info_from_ID(MyIDType ID)
{
  int i;
  for(i = 0; i < NumPart; i++)
    if(P[i].ID == ID)
      print_particle_info(i);
}

/** \brief Print information relative to the left or right state of a face to standard output.
 *
 *  \param st structure containing the left or right state of a face.
 */
void print_state_info(struct state *st)
{
  printf("Task=%d, ID=%llu rho=%g, p=%g, vel=%g|%g|%g, velVertex=%g|%g|%g\n", ThisTask, (unsigned long long) st->ID,
         st->rho, st->press, st->velx, st->vely, st->velz, st->velVertex[0], st->velVertex[1], st->velVertex[2]);
  printf("dx=%g, dy=%g, dz=%g, dt_half=%g\n", st->dx, st->dy, st->dz, st->dt_half);
  printf("timeBin=%d, volume=%g, activearea=%g, surfacearea=%g, csnd=%g\n", st->timeBin, st->volume, st->activearea, st->surfacearea, st->csnd);
#if defined(EOS_DEGENERATE) || defined(EOS_OPAL)
  printf("gammaC=%g, gammaE=%g\n", st->gammaC, st->gammaE);
#endif
#ifdef MHD
  printf("B=%g|%g|%g\n", st->Bx, st->By, st->Bz);
#endif
}

/** \brief Print information relative to state the of a face as determined by the Riemman solver to standard output.
 *
 *  \param st structure containing the state of a face after the solution of the Riemann problem
 */
void print_state_face_info(struct state_face *st)
{
  printf("rho=%g, p=%g, vel=%g|%g|%g\n", st->rho, st->press, st->velx, st->vely, st->velz);
}
