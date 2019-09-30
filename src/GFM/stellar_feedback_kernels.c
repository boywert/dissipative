/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/GFM/stellar_feedback_kernels.c
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

#include "../allvars.h"
#include "../proto.h"
#include "stellar_feedback_kernels.h"

int is_doing_stellar_feedback(int i)
{
#ifdef FM_STAR_FEEDBACK
  if(StarParticle[i].TotalEnergyReleased > 0)
    return 1;
#endif
#ifdef GFM_STELLAR_FEEDBACK
  if(StarParticle[i].SNIaEnergyReleased > 0 || StarParticle[i].AGBMomentumReleased > 0)
    return 1;
#endif
#ifdef GFM_WINDS_LOCAL
  if(StarParticle[i].WindEnergyReleased > 0)
    return 1;
#endif

  return 0;
}

#ifdef GFM_STELLAR_FEEDBACK
void GFM_stellar_feedback(int target, int mode, int thread_id, int numnodes, int *firstnode, data_in * in)
{
  MyDouble inj_mom[3], *pos, SNIaEnergyReleased, AGBMomentumReleased;
  MyFloat normsph;
  double weight_fac, h, h2;
  double dr[3], r2, r;

  pos = in->Pos;
  h = in->Hsml;
  normsph = in->NormSph;
  SNIaEnergyReleased = in->SNIaEnergyReleased * All.cf_atime * All.cf_atime;
  AGBMomentumReleased = in->AGBMomentumReleased * All.cf_atime;

#ifndef GFM_TOPHAT_KERNEL
  double wk;
  MyFloat hinv = 1.0 / h;
#ifndef TWODIMS
  MyFloat hinv3 = hinv * hinv * hinv;
#else
  MyFloat hinv3 = hinv * hinv / boxSize_Z;
#endif
#endif

#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif

  h2 = h * h;

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, thread_id, numnodes, firstnode);

  for(int n = 0; n < nfound; n++)
    {
      int j = Thread[thread_id].Ngblist[n];

      if(P[j].Mass > 0 && P[j].ID != 0) /* skip cells that have been swallowed or dissolved */
        {
          dr[0] = NEAREST_X(pos[0] - P[j].Pos[0]);
          dr[1] = NEAREST_Y(pos[1] - P[j].Pos[1]);
          dr[2] = NEAREST_Z(pos[2] - P[j].Pos[2]);

          r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];

          if(r2 < h2)
            {
              r = sqrt(r2);

#ifndef GFM_TOPHAT_KERNEL
              double u = r * hinv;

              if(u < 0.5)
                wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
              else
                wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

              weight_fac = SphP[j].Volume * wk / normsph;
#else
              weight_fac = SphP[j].Volume / normsph;
#endif

              /* injected feedback energy */
              MyDouble de_feedback = weight_fac * SNIaEnergyReleased;
              /* injected feedback momentum */
              MyDouble dp_feedback = weight_fac * AGBMomentumReleased;

              inj_mom[0] = -dp_feedback * dr[0] / r;
              inj_mom[1] = -dp_feedback * dr[1] / r;
              inj_mom[2] = -dp_feedback * dr[2] / r;

              /* this accounts both for thermal and kinetic energy from feedback */
              SphP[j].Energy += de_feedback;
              /* momentum due to stellar mass return is injected in radial direction */
              SphP[j].Momentum[0] += inj_mom[0];
              SphP[j].Momentum[1] += inj_mom[1];
              SphP[j].Momentum[2] += inj_mom[2];

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
              SphP[j].Utherm =
                (SphP[j].Energy -
                 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] +
                        SphP[j].Momentum[2] * SphP[j].Momentum[2]) / P[j].Mass) / P[j].Mass / (All.cf_atime * All.cf_atime);
              SphP[j].A = (GAMMA - 1.0) * SphP[j].Utherm / pow(SphP[j].Density * All.cf_a3inv, GAMMA - 1);
              SphP[j].Entropy = log(SphP[j].A) * P[j].Mass;
#endif
            }
        }
    }
}
#endif

#ifdef GFM_WINDS_LOCAL
void GFM_winds_local(int target, int mode, int thread_id, int numnodes, int *firstnode, data_in *in)
{
  MyDouble *pos, WindEnergyReleased;
  MyFloat normsph;
  double weight_fac, h, h2;
  double dr[3], r2, r;

  pos = in->Pos;
  h = in->Hsml;
  normsph = in->NormSph;
  WindEnergyReleased = in->WindEnergyReleased;

#ifndef GFM_TOPHAT_KERNEL
  double wk;
  MyFloat hinv = 1.0 / h;
#ifndef TWODIMS
  MyFloat hinv3 = hinv * hinv * hinv;
#else
  MyFloat hinv3 = hinv * hinv / boxSize_Z;
#endif
#endif

#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif

  h2 = h * h;

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, thread_id, numnodes, firstnode);

  for(int n = 0; n < nfound; n++)
    {
      int j = Thread[thread_id].Ngblist[n];

      if(P[j].Mass > 0 && P[j].ID != 0) /* skip cells that have been swallowed or dissolved */
        {
          dr[0] = NEAREST_X(pos[0] - P[j].Pos[0]);
          dr[1] = NEAREST_Y(pos[1] - P[j].Pos[1]);
          dr[2] = NEAREST_Z(pos[2] - P[j].Pos[2]);

          r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];

          if(r2 < h2)
            {
              r = sqrt(r2);

#ifndef GFM_TOPHAT_KERNEL
              double u = r * hinv;

              if(u < 0.5)
                wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
              else
                wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

              weight_fac = SphP[j].Volume * wk / normsph;
#else
              weight_fac = SphP[j].Volume / normsph;
#endif
              SphP[j].WindEnergyReceived += weight_fac * WindEnergyReleased;

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
              SphP[j].Utherm =
                (SphP[j].Energy -
                 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] +
                        SphP[j].Momentum[2] * SphP[j].Momentum[2]) / P[j].Mass) / P[j].Mass / (All.cf_atime * All.cf_atime);
              SphP[j].A = (GAMMA - 1.0) * SphP[j].Utherm / pow(SphP[j].Density * All.cf_a3inv, GAMMA - 1);
              SphP[j].Entropy = log(SphP[j].A) * P[j].Mass;
#endif
            }
        }
    }
}
#endif

#ifdef FM_STAR_FEEDBACK

#ifdef DELAYED_COOLING
void delayed_cooling(int target, int mode, int thread_id, int numnodes, int *firstnode, data_in *in, data_out *out)
{
  int n_cells_not_cooling = 0;
  double weight_fac, dr[3], r2, r, h, h2;

  MyDouble TotalEnergyReleased, TotalKinEnergyReleased;
  MyFloat normsph_feed, BlastRadius, CoolShutoffTime;
  MyDouble de_feedback, *pos;

  pos = in->Pos;
  h = in->Hsml;

  TotalEnergyReleased = in->TotalEnergyReleased * All.cf_atime * All.cf_atime;
  TotalKinEnergyReleased = All.EtaKineticEnergy * TotalEnergyReleased;

  normsph_feed = in->NormSphFeedback;
  BlastRadius = in->BlastRadius;
  CoolShutoffTime = in->CoolShutoffTime;

#ifndef GFM_TOPHAT_KERNEL
  double wk;
  MyFloat blast_inv = 1.0 / BlastRadius;
#ifndef TWODIMS
  MyFloat blast_inv3 = blast_inv * blast_inv * blast_inv;
#else
  MyFloat blast_inv3 = blast_inv * blast_inv / boxSize_Z;
#endif
#endif

#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif

  h2 = h * h;

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, thread_id, numnodes, firstnode);

  for(int n = 0; n < nfound; n++)
    {
      int j = Thread[thread_id].Ngblist[n];

      if(P[j].Mass > 0 && P[j].ID != 0) /* skip cells that have been swallowed or dissolved */
        {
          dr[0] = NEAREST_X(pos[0] - P[j].Pos[0]);
          dr[1] = NEAREST_Y(pos[1] - P[j].Pos[1]);
          dr[2] = NEAREST_Z(pos[2] - P[j].Pos[2]);

          r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];

          if(r2 < h2)
            {
              r = sqrt(r2);
#ifdef FM_STAR_FEEDBACK
              if(TotalEnergyReleased > 0)
                {
                  /* disable cooling for the particle if there are SN explosions */
                  if(r <= BlastRadius)
                    {
#ifndef GFM_TOPHAT_KERNEL
                      double u = r * blast_inv;

                      if(u < 0.5)
                        wk = blast_inv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
                      else
                        wk = blast_inv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

                      weight_fac = SphP[j].Volume * wk / normsph_feed;
#else
                      weight_fac = SphP[j].Volume / normsph_feed;
#endif
                      de_feedback = weight_fac * TotalEnergyReleased;
                      SphP[j].Energy += de_feedback;
                      set_cooling_shutoff_time(j, CoolShutoffTime);

#ifdef OUTPUT_STELLAR_FEEDBACK
                      /* save cumulative feedback energy for debugging */
                      SphP[j].TotEgyFeed += de_feedback;
                      SphP[j].IntEgyFeed += (1.0 - All.EtaKineticEnergy) * de_feedback;
                      SphP[j].KinEgyFeed += weight_fac * TotalKinEnergyReleased;
#endif
                      n_cells_not_cooling++;
                    }
                }
#endif

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
              SphP[j].Utherm =
                (SphP[j].Energy -
                 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] +
                        SphP[j].Momentum[2] * SphP[j].Momentum[2]) / P[j].Mass) / P[j].Mass / (All.cf_atime * All.cf_atime);
              SphP[j].A = (GAMMA - 1.0) * SphP[j].Utherm / pow(SphP[j].Density * All.cf_a3inv, GAMMA - 1);
              SphP[j].Entropy = log(SphP[j].A) * P[j].Mass;
#endif
            }
        }
    }

  out->kicked_cells = n_cells_not_cooling;
  out->BlastRadius = BlastRadius;
  out->ShutoffTime = CoolShutoffTime;
}
#endif

#ifdef DELAYED_COOLING_TURB
void delayed_cooling_turbulence(int target, int mode, int thread_id, int numnodes, int *firstnode, data_in *in)
{
  double weight_fac, dr[3], r2, r, h, h2;
  MyDouble de_feedback, dturb_feedback, *pos;
  MyDouble TotalEnergyReleased, TotalKinEnergyReleased;
  MyFloat normsph;

  pos = in->Pos;
  h = in->Hsml;

  TotalEnergyReleased = in->TotalEnergyReleased * All.cf_atime * All.cf_atime;
  TotalKinEnergyReleased = All.EtaKineticEnergy * TotalEnergyReleased;

  normsph = in->NormSph;

#ifndef GFM_TOPHAT_KERNEL
  double wk;
  MyFloat hinv = 1.0 / h;
#ifndef TWODIMS
  MyFloat hinv3 = hinv * hinv * hinv;
#else
  MyFloat hinv3 = hinv * hinv / boxSize_Z;
#endif
#endif

#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif

  h2 = h * h;

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, thread_id, numnodes, firstnode);

  for(int n = 0; n < nfound; n++)
    {
      int j = Thread[thread_id].Ngblist[n];

      if(P[j].Mass > 0 && P[j].ID != 0) /* skip cells that have been swallowed or dissolved */
        {
          dr[0] = NEAREST_X(pos[0] - P[j].Pos[0]);
          dr[1] = NEAREST_Y(pos[1] - P[j].Pos[1]);
          dr[2] = NEAREST_Z(pos[2] - P[j].Pos[2]);

          r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];

          if(r2 < h2)
            { 
              r = sqrt(r2);
              if(TotalEnergyReleased > 0)
                {
#ifndef GFM_TOPHAT_KERNEL
                  double u = r * hinv;

                  if(u < 0.5)
                    wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
                  else
                    wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

                  weight_fac = SphP[j].Volume * wk / normsph;
#else
                  weight_fac = SphP[j].Volume / normsph;
#endif

                  /* injected feedback energy */
                  de_feedback = weight_fac * TotalEnergyReleased;
                  dturb_feedback = weight_fac * TotalKinEnergyReleased;
                  /* this accounts both for thermal and kinetic energy from feedback */
                  SphP[j].Energy += de_feedback;
                  /* even for cosmological run this is a physical energy (is that ok?) */
                  //SphP[j].MassUturb += de_feedback / All.cf_atime / All.cf_atime;
                  SphP[j].MassUturb += dturb_feedback / All.cf_atime / All.cf_atime;

#ifdef OUTPUT_STELLAR_FEEDBACK
                  /* save cumulative feedback energy for debugging */
                  SphP[j].TotEgyFeed += de_feedback;
                  SphP[j].IntEgyFeed += (1.0 - All.EtaKineticEnergy) * de_feedback;
                  SphP[j].KinEgyFeed += weight_fac * TotalKinEnergyReleased;
#endif
                }

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
              SphP[j].Utherm =
                (SphP[j].Energy -
                 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] +
                        SphP[j].Momentum[2] * SphP[j].Momentum[2]) / P[j].Mass) / P[j].Mass / (All.cf_atime * All.cf_atime);
              SphP[j].A = (GAMMA - 1.0) * SphP[j].Utherm / pow(SphP[j].Density * All.cf_a3inv, GAMMA - 1);
              SphP[j].Entropy = log(SphP[j].A) * P[j].Mass;
#endif
            }
        }
    }
}
#endif

#ifdef NON_STOCHASTIC_MOMENTUM_FEEDBACK
void momentum_feedback(int target, int mode, int thread_id, int numnodes, int *firstnode, data_in *in, data_out *out)
{
  double weight_fac, dr[3], r2, r, h, h2, h3;
  MyFloat normsph;
  MyDouble inj_mom[3], *pos;
  MyDouble TotalEnergyReleased, TotalKinEnergyReleased;
  MyDouble dptot = 0;

#ifdef INJECT_INTO_SINGLE_CELL
  MyIDType ClosestNeighbourID = 0;
#endif
#ifdef FM_MASS_WEIGHT_SN
  MyFloat totngbmass;
#endif

  pos = in->Pos;
  h = in->Hsml;
  normsph = in->NormSph;

#ifdef FM_MASS_WEIGHT_SN
  totngbmass = in->TotNgbMass;
#endif

  TotalEnergyReleased = in->TotalEnergyReleased * All.cf_atime * All.cf_atime;
  TotalKinEnergyReleased = All.EtaKineticEnergy * TotalEnergyReleased;

#ifdef INJECT_INTO_SINGLE_CELL
  ClosestNeighbourID = in->ClosestNeighbourID;
#endif

#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif

  h2 = h * h;

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, thread_id, numnodes, firstnode);

  for(int n = 0; n < nfound; n++)
    {
      int j = Thread[thread_id].Ngblist[n];

      if(P[j].Mass > 0 && P[j].ID != 0) /* skip cells that have been swallowed or dissolved */
        {
          dr[0] = NEAREST_X(pos[0] - P[j].Pos[0]);
          dr[1] = NEAREST_Y(pos[1] - P[j].Pos[1]);
          dr[2] = NEAREST_Z(pos[2] - P[j].Pos[2]);

          r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];

          if(r2 < h2)
            {
              r = sqrt(r2);
              if(TotalEnergyReleased > 0)
                {
#ifdef INJECT_INTO_SINGLE_CELL
                  if(ClosestNeighbourID == 0)
                    terminate("No neighbours found for this stellar particle");

                  if(P[j].ID == ClosestNeighbourID)
                    {
                      /* Compute internal energy */
                      double EKin = 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] + SphP[j].Momentum[2] * SphP[j].Momentum[2]) / P[j].Mass;
                      SphP[j].Energy -= EKin;

                      weight_fac = 1.; 
                      MyDouble de_feedback = weight_fac * TotalEnergyReleased;
                      MyDouble dp_feedback = sqrt(2.0 * weight_fac * P[j].Mass * TotalKinEnergyReleased);

                      inj_mom[0] = -dp_feedback * dr[0] / r;
                      inj_mom[1] = -dp_feedback * dr[1] / r;
                      inj_mom[2] = -dp_feedback * dr[2] / r;

                      /* this accounts for the thermal energy from feedback */
                      SphP[j].Energy += (1.0 - All.EtaKineticEnergy) * de_feedback;

                      /* momentum is injected in radial direction */
                      SphP[j].Momentum[0] += inj_mom[0];
                      SphP[j].Momentum[1] += inj_mom[1];
                      SphP[j].Momentum[2] += inj_mom[2];

                      /* Compute new kinetic energy and add to the internal energy */
                      EKin = 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] + SphP[j].Momentum[2] * SphP[j].Momentum[2]) / P[j].Mass;
                      SphP[j].Energy += EKin;

#ifdef OUTPUT_STELLAR_FEEDBACK
                      /* save cumulative feedback energy for debugging */
                      SphP[j].TotEgyFeed += de_feedback;
                      SphP[j].IntEgyFeed += (1.0 - All.EtaKineticEnergy) * de_feedback;
                      SphP[j].KinEgyFeed += weight_fac * TotalKinEnergyReleased;
#endif
                    }

#else

#if defined(FM_MASS_WEIGHT_SN) || !defined(GFM_TOPHAT_KERNEL)
                  double wk;
                  MyFloat hinv = 1.0 / h;
#ifndef TWODIMS
                  MyFloat hinv3 = hinv * hinv * hinv;
#else
                  MyFloat hinv3 = hinv * hinv / boxSize_Z;
#endif
                  h3 = 1.0 / hinv3;
#endif
                  double mass_j = P[j].Mass;

#if defined(FM_MASS_WEIGHT_SN) || !defined(GFM_TOPHAT_KERNEL)
                  double u = r * hinv;

                  if(u < 0.5)
                    wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
                  else
                    wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

#ifdef FM_MASS_WEIGHT_SN  
                  weight_fac = NORM_COEFF * P[j].Mass * wk * h3 / totngbmass;
		  mass_j *= NORM_COEFF * wk * h3;
#else
                  weight_fac = SphP[j].Volume * wk / normsph;
#endif

#else
                  weight_fac = SphP[j].Volume / normsph;
#endif


                  /* Compute internal energy */
                  double EKin = 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] + SphP[j].Momentum[2] * SphP[j].Momentum[2]) / P[j].Mass;
                  SphP[j].Energy -= EKin;

                  MyDouble de_feedback = weight_fac * TotalEnergyReleased;
                  MyDouble dp_feedback = sqrt(2.0 * weight_fac * mass_j * TotalKinEnergyReleased);
                  dptot += dp_feedback;

                  /* this accounts for the thermal energy from feedback */
                  SphP[j].Energy += (1.0 - All.EtaKineticEnergy) * de_feedback;

#if (FM_STAR_FEEDBACK_KICK_TYPE == 0)
                  double theta = acos(2 * get_random_number() - 1);
                  double phi = 2 * M_PI * get_random_number();
                  inj_mom[0] = -dp_feedback * sin(theta) * cos(phi);
                  inj_mom[1] = -dp_feedback * sin(theta) * sin(phi);
                  inj_mom[2] = -dp_feedback * cos(theta);
#endif
#if (FM_STAR_FEEDBACK_KICK_TYPE == 1)
                  inj_mom[0] = -dp_feedback * dr[0] / r;
                  inj_mom[1] = -dp_feedback * dr[1] / r;
                  inj_mom[2] = -dp_feedback * dr[2] / r;
#endif
#if (FM_STAR_FEEDBACK_KICK_TYPE == 2)
                  double dir[3];
                  double norm;
                  dir[0] = P[j].GravAccel[1] * SphP[j].Momentum[2] - P[j].GravAccel[2] * SphP[j].Momentum[1];
                  dir[1] = P[j].GravAccel[2] * SphP[j].Momentum[0] - P[j].GravAccel[0] * SphP[j].Momentum[2];
                  dir[2] = P[j].GravAccel[0] * SphP[j].Momentum[1] - P[j].GravAccel[1] * SphP[j].Momentum[0];
                  norm = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);

                  if(get_random_number() < 0.5)
                    norm = -norm;

                  /* what happens if norm is zero ??? */
                  inj_mom[0] = -dp_feedback * dir[0] / norm;
                  inj_mom[1] = -dp_feedback * dir[1] / norm;
                  inj_mom[2] = -dp_feedback * dir[2] / norm;
#endif

                  /* momentum is injected in radial direction */
                  SphP[j].Momentum[0] += inj_mom[0];
                  SphP[j].Momentum[1] += inj_mom[1];
                  SphP[j].Momentum[2] += inj_mom[2];

                  /* Compute new kinetic energy and add to the internal energy */
                  EKin = 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] + SphP[j].Momentum[2] * SphP[j].Momentum[2]) / P[j].Mass;
                  SphP[j].Energy += EKin;

#ifdef OUTPUT_STELLAR_FEEDBACK
                  /* save cumulative feedback energy for debugging */
                  SphP[j].TotEgyFeed += de_feedback;
                  SphP[j].IntEgyFeed += (1.0 - All.EtaKineticEnergy) * de_feedback;
                  SphP[j].KinEgyFeed += weight_fac * TotalKinEnergyReleased;
#endif

#endif
                }
#ifdef USE_ENTROPY_FOR_COLD_FLOWS
              SphP[j].Utherm =
                (SphP[j].Energy -
                 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] +
                        SphP[j].Momentum[2] * SphP[j].Momentum[2]) / P[j].Mass) / P[j].Mass / (All.cf_atime * All.cf_atime);
              SphP[j].A = (GAMMA - 1.0) * SphP[j].Utherm / pow(SphP[j].Density * All.cf_a3inv, GAMMA - 1);
              SphP[j].Entropy = log(SphP[j].A) * P[j].Mass;
#endif
            }
        }
    }

  out->TotalMomentumInjected = dptot;
}
#endif

#ifdef DIRECT_MOMENTUM_INJECTION_FEEDBACK
void direct_momentum_feedback(int target, int mode, int thread_id, int numnodes, int *firstnode, data_in *in, data_out *out)
{
  double weight_fac, dr[3], r2, r, h, h2;
#if defined(FM_MASS_WEIGHT_SN) || !defined(GFM_TOPHAT_KERNEL)
  double h3;
#endif
  MyFloat normsph;
  MyDouble inj_mom[3], *pos;
  MyDouble TotalEnergyReleased, TotalMomentumReleased;
  MyDouble dptot = 0;

#ifdef FM_MASS_WEIGHT_SN
  MyFloat totngbmass;
#endif

  pos = in->Pos;
  h = in->Hsml;
  normsph = in->NormSph;

#ifdef FM_MASS_WEIGHT_SN
  totngbmass = in->TotNgbMass;
#endif

  TotalEnergyReleased = in->TotalEnergyReleased * All.cf_atime * All.cf_atime;
  TotalMomentumReleased = in->TotalMomentumReleased * All.cf_atime * All.cf_atime;

#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif

  h2 = h * h;

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, thread_id, numnodes, firstnode);

  for(int n = 0; n < nfound; n++)
    {
      int j = Thread[thread_id].Ngblist[n];

      if(P[j].Mass > 0 && P[j].ID != 0) /* skip cells that have been swallowed or dissolved */
        {
          dr[0] = NEAREST_X(pos[0] - P[j].Pos[0]);
          dr[1] = NEAREST_Y(pos[1] - P[j].Pos[1]);
          dr[2] = NEAREST_Z(pos[2] - P[j].Pos[2]);

          r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];

          if(r2 < h2)
            {
              r = sqrt(r2);
              if(TotalEnergyReleased > 0)
                {
#if defined(FM_MASS_WEIGHT_SN) || !defined(GFM_TOPHAT_KERNEL)
                  double wk;
                  MyFloat hinv = 1.0 / h;
#ifndef TWODIMS
                  MyFloat hinv3 = hinv * hinv * hinv;
#else
                  MyFloat hinv3 = hinv * hinv / boxSize_Z;
#endif
                  h3 = 1.0 / hinv3;
#endif

#if defined(FM_MASS_WEIGHT_SN) || !defined(GFM_TOPHAT_KERNEL)
                  double u = r * hinv;

                  if(u < 0.5)
                    wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
                  else
                    wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

#ifdef FM_MASS_WEIGHT_SN  
                  weight_fac = NORM_COEFF * P[j].Mass * wk * h3 / totngbmass;
#else
                  weight_fac = SphP[j].Volume * wk / normsph;
#endif

#else
                  weight_fac = SphP[j].Volume / normsph;
#endif

                  /* Compute internal energy */
                  double EKin = 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] + SphP[j].Momentum[2] * SphP[j].Momentum[2]) / P[j].Mass;
                  SphP[j].Energy -= EKin;

                  MyDouble de_feedback = weight_fac * TotalEnergyReleased;
                  MyDouble dp_feedback = weight_fac * TotalMomentumReleased;
                  dptot += dp_feedback;

                  /* this accounts for the thermal energy from feedback */
                  SphP[j].Energy += (1.0 - All.EtaKineticEnergy) * de_feedback;

#if (FM_STAR_FEEDBACK_KICK_TYPE == 0)
                  double theta = acos(2 * get_random_number() - 1);
                  double phi = 2 * M_PI * get_random_number();
                  inj_mom[0] = -dp_feedback * sin(theta) * cos(phi);
                  inj_mom[1] = -dp_feedback * sin(theta) * sin(phi);
                  inj_mom[2] = -dp_feedback * cos(theta);
#endif
#if (FM_STAR_FEEDBACK_KICK_TYPE == 1)
                  inj_mom[0] = -dp_feedback * dr[0] / r;
                  inj_mom[1] = -dp_feedback * dr[1] / r;
                  inj_mom[2] = -dp_feedback * dr[2] / r;
#endif
#if (FM_STAR_FEEDBACK_KICK_TYPE == 2)
                  double dir[3];
                  double norm;
                  dir[0] = P[j].GravAccel[1] * SphP[j].Momentum[2] - P[j].GravAccel[2] * SphP[j].Momentum[1];
                  dir[1] = P[j].GravAccel[2] * SphP[j].Momentum[0] - P[j].GravAccel[0] * SphP[j].Momentum[2];
                  dir[2] = P[j].GravAccel[0] * SphP[j].Momentum[1] - P[j].GravAccel[1] * SphP[j].Momentum[0];
                  norm = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);

                  if(get_random_number() < 0.5)
                    norm = -norm;

                  /* what happens if norm is zero ??? */
                  inj_mom[0] = -dp_feedback * dir[0] / norm;
                  inj_mom[1] = -dp_feedback * dir[1] / norm;
                  inj_mom[2] = -dp_feedback * dir[2] / norm;
#endif

                  /* momentum is injected in radial direction */
                  SphP[j].Momentum[0] += inj_mom[0];
                  SphP[j].Momentum[1] += inj_mom[1];
                  SphP[j].Momentum[2] += inj_mom[2];

                  /* Compute new kinetic energy and add to the internal energy */
                  EKin = 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] + SphP[j].Momentum[2] * SphP[j].Momentum[2]) / P[j].Mass;
                  SphP[j].Energy += EKin;

#ifdef OUTPUT_STELLAR_FEEDBACK
                  /* save cumulative feedback energy for debugging */
                  SphP[j].TotEgyFeed += de_feedback;
                  SphP[j].IntEgyFeed += (1.0 - All.EtaKineticEnergy) * de_feedback;
                  SphP[j].KinEgyFeed += All.EtaKineticEnergy * de_feedback;
#endif
                }
#ifdef USE_ENTROPY_FOR_COLD_FLOWS
              SphP[j].Utherm =
                (SphP[j].Energy -
                 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] +
                        SphP[j].Momentum[2] * SphP[j].Momentum[2]) / P[j].Mass) / P[j].Mass / (All.cf_atime * All.cf_atime);
              SphP[j].A = (GAMMA - 1.0) * SphP[j].Utherm / pow(SphP[j].Density * All.cf_a3inv, GAMMA - 1);
              SphP[j].Entropy = log(SphP[j].A) * P[j].Mass;
#endif
            }
        }
    }

  out->TotalMomentumInjected = dptot;
}
#endif


#ifdef FM_SN_COOLING_RADIUS_BOOST
void cooling_radius_momentum_feedback(int target, int mode, int thread_id, int numnodes, int *firstnode, data_in *in, data_out *out)
{
  int n, j, k, n_int_SN, nfound;
  double weight_fac, dr[3], r2, r, h, h2,  normsph;
#ifdef FM_MASS_WEIGHT_SN
  double  wk, u, h3, hinv, hinv3;
#endif
  double n0, z0;                /* properties of the local ISM for cooling radius calc */
  double cool_rad, mass_enc, momentum_boost, n_SN, e_sn_shock, e_sn_shock_resolved = 0;
#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif

  MyDouble *pos;
  MyDouble TotalEnergyReleased, TotalMomentumReleased, TotalMassReleased, dp, du;
  MyDouble dptot = 0;

  pos = in->Pos;
  h = in->Hsml;

  normsph = in->NormSph;
#ifdef FM_MASS_WEIGHT_SN
  MyFloat totngbmass = in->TotNgbMass;
#endif
  n_SN = in->n_SNII;

  n_int_SN=floor(n_SN);  /* rounded down integer number of SN in this timestep */
  n_SN -= n_int_SN;
  if(get_random_number() < n_SN) n_int_SN++;

  if(n_int_SN==0) return;

  TotalMassReleased     = n_int_SN * All.one_SNe_mass;
  TotalEnergyReleased   = n_int_SN * All.one_SNe_energy;
  TotalMomentumReleased = n_int_SN * All.one_SNe_mass * All.SNe_velocity;

  /* ToDo:  include mean molecular weight? */
  n0 = in->LocISMdens / PROTONMASS * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;                                 /* local ISM density (code units) */
  z0 = in->LocISMZdens / in->LocISMdens;     /* local ISM metallicity          */

  /* cooling radius taken from pg 599 of Hopkins+ (2014MNRAS.445..581H) */
  /* orig calculation of cooling radius taken from Cioffi+ (1988ApJ...334..252C) in weak external pressure limit*/
  cool_rad = (28.0*3.086e+18/All.UnitLength_in_cm) * pow(All.FeedbackEfficiency, 0.286) * pow(n0,-0.43) * pow(z0+0.01,-0.18);
  /* ToDo:  Do we need a hubble factor */
  mass_enc = 4.18879*(in->LocISMdens)*cool_rad*cool_rad*cool_rad;
  e_sn_shock = 0.5 * All.one_SNe_mass * n_int_SN * All.SNe_velocity * All.SNe_velocity;                              /* in code units */

  MyFloat mom_at_cool_rad = TotalMomentumReleased * sqrt(1.0 + mass_enc / TotalMassReleased  );
  MyFloat egy_at_cool_rad = TotalEnergyReleased - 0.5 * mom_at_cool_rad * mom_at_cool_rad / mass_enc;               /* this is a total energy */

  h2 = h * h;

  nfound = ngb_treefind_variable_threads(pos, h, target, mode, thread_id, numnodes, firstnode);

  for(n = 0; n < nfound; n++)
    {
      j = Thread[thread_id].Ngblist[n];
      if(P[j].Mass > 0 && P[j].ID != 0) /* skip cells that have been swallowed or dissolved */
        {
          dr[0] = NEAREST_X(pos[0] - P[j].Pos[0]);
          dr[1] = NEAREST_Y(pos[1] - P[j].Pos[1]);
          dr[2] = NEAREST_Z(pos[2] - P[j].Pos[2]);
          r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];

          if(r2 < h2)
            {
              r = sqrt(r2);
              if(TotalEnergyReleased > 0)
                {
#ifdef FM_MASS_WEIGHT_SN
                  hinv = 1.0 / h;
                  hinv3 = hinv * hinv * hinv;
                  h3 = 1.0 / hinv3;
                  u = r * hinv;
                  if(u < 0.5)  wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
                  else  wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
                  weight_fac = NORM_COEFF * P[j].Mass * wk * h3 / totngbmass;
#else
                  weight_fac = SphP[j].Volume / normsph;
#endif
                  if(h>cool_rad)
                  {
                      dp = mom_at_cool_rad * weight_fac;
                      du = egy_at_cool_rad * (cool_rad*cool_rad*cool_rad*cool_rad*cool_rad*cool_rad)/(h*h*h*h*h*h) * weight_fac / P[j].Mass;
                  } else {
                      // TODO:  Fix this for partially resolved cooling radius limit.  
                      dp = mom_at_cool_rad * weight_fac;
                      du = egy_at_cool_rad * weight_fac / P[j].Mass;
                  }

                  MyFloat f_startup = All.Time * All.Time / (0.01*0.01 + All.Time * All.Time);

                  if(r>0.1){
                    du = 0;
                    dp = 0;
                  }
                  for(k=0;k<3;k++) SphP[j].Energy -= 0.5 * SphP[j].Momentum[k] * SphP[j].Momentum[k] / P[j].Mass;

                  SphP[j].Energy      += du * f_startup * P[j].Mass;
                  SphP[j].Utherm      += du * f_startup;                // probably useless
                  for(k=0;k<3;k++) {
                        SphP[j].Momentum[k] -= dp * f_startup * dr[k] / r;
                  }
                  dptot = dp;

                  for(k=0;k<3;k++) SphP[j].Energy += 0.5 * SphP[j].Momentum[k] * SphP[j].Momentum[k] / P[j].Mass;
                }
            }
        }
    }
  out->TotalMomentumInjected = dptot;
}
#endif

#if !defined(DELAYED_COOLING) && !defined(DELAYED_COOLING_TURB) && !defined(NON_STOCHASTIC_MOMENTUM_FEEDBACK) && !defined(DIRECT_MOMENTUM_INJECTION_FEEDBACK)
void stochastic_momentum_feedback(int target, int mode, int thread_id, int numnodes, int *firstnode, data_in *in, data_out *out)
{
  int n_cells_kicked = 0;
  double dr[3], r2, r, h, h2;

  MyDouble dp_feedback, deltaEKin, TotalEnergyReleased, TotalKinEnergyReleased, *pos;
  MyFloat normsph, prob, velocity, tot_mass;
  MyDouble mass_kicked = 0, sum_deltaEKin = 0;
  MyDouble sum_deltaMomentum[3] = {0.0, 0.0, 0.0};
  MyDouble de_feedback, EKin;
  MyDouble inj_mom[3], dp_tot = 0;

  pos = in->Pos;
  h = in->Hsml;
  normsph = in->NormSph;
  TotalEnergyReleased = in->TotalEnergyReleased * All.cf_atime * All.cf_atime;
  TotalKinEnergyReleased = All.EtaKineticEnergy * TotalEnergyReleased;
  tot_mass = in->NumNgb * All.ReferenceGasPartMass;

  velocity = sqrt(All.G * tot_mass / h);
  velocity *= 2.;    /* factor of 2 to be on the safe side. LVS: check cosmology factor?? */

  velocity = dmax(All.FeedbackVelocity, velocity) * All.cf_atime;

  /* couldn't find any neighbour (it can happen in zoom runs in the low-low res region) */
  if(tot_mass <= 0)
    prob = 0.0;
  else
    {
      prob = 2.0 * TotalKinEnergyReleased / (velocity * velocity * tot_mass);
      velocity *= dmax(1.0, sqrt(prob) + 1.0e-6);
      prob = 2.0 * TotalKinEnergyReleased / (velocity * velocity * tot_mass);
    }

  if(prob > 1.0)
    terminate("Total mass within kernel is less than mass to be kicked for stochastic feedback: tot_mass %g, kick_mass %g", tot_mass, prob * tot_mass);

#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif

  h2 = h * h;

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, thread_id, numnodes, firstnode);

  for(int n = 0; n < nfound; n++)
    {
      int j = Thread[thread_id].Ngblist[n];

      if(P[j].Mass > 0 && P[j].ID != 0) /* skip cells that have been swallowed or dissolved */
        {
          dr[0] = NEAREST_X(pos[0] - P[j].Pos[0]);
          dr[1] = NEAREST_Y(pos[1] - P[j].Pos[1]);
          dr[2] = NEAREST_Z(pos[2] - P[j].Pos[2]);

          r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];

          if(r2 < h2)
            {
              r = sqrt(r2);
              if(TotalEnergyReleased > 0)
                {
                  if(prob > 1.0)
                    terminate("Total mass within kernel is less than mass to be kicked for stochastic feedback: tot_mass %g, kick_mass %g", tot_mass, prob * tot_mass);

                  /* do feedback */
                  if(get_random_number() < prob)
                    {
                      /* injected feedback momentum */
                      dp_feedback = P[j].Mass * velocity;
                      dp_tot += dp_feedback;
#if (FM_STAR_FEEDBACK_KICK_TYPE == 0)
                      double theta = acos(2 * get_random_number() - 1);
                      double phi = 2 * M_PI * get_random_number();
                      inj_mom[0] = -dp_feedback * sin(theta) * cos(phi);
                      inj_mom[1] = -dp_feedback * sin(theta) * sin(phi);
                      inj_mom[2] = -dp_feedback * cos(theta);
#endif
#if (FM_STAR_FEEDBACK_KICK_TYPE == 1)
                      inj_mom[0] = -dp_feedback * dr[0] / r;
                      inj_mom[1] = -dp_feedback * dr[1] / r;
                      inj_mom[2] = -dp_feedback * dr[2] / r;
#endif
#if (FM_STAR_FEEDBACK_KICK_TYPE == 2)
                      double dir[3];
                      double norm;
                      dir[0] = P[j].GravAccel[1] * SphP[j].Momentum[2] - P[j].GravAccel[2] * SphP[j].Momentum[1];
                      dir[1] = P[j].GravAccel[2] * SphP[j].Momentum[0] - P[j].GravAccel[0] * SphP[j].Momentum[2];
                      dir[2] = P[j].GravAccel[0] * SphP[j].Momentum[1] - P[j].GravAccel[1] * SphP[j].Momentum[0];
                      norm = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);

                      if(get_random_number() < 0.5)
                        norm = -norm;

                      /* what happens if norm is zero ??? */
                      inj_mom[0] = -dp_feedback * dir[0] / norm;
                      inj_mom[1] = -dp_feedback * dir[1] / norm;
                      inj_mom[2] = -dp_feedback * dir[2] / norm;
#endif
                      /* Compute internal energy */
                      EKin = 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] + SphP[j].Momentum[2] * SphP[j].Momentum[2]) / P[j].Mass;
                      SphP[j].Energy -= EKin;
                      deltaEKin = -EKin;

                      /* this accounts for the thermal energy from feedback */
                      de_feedback = 0.5 * dp_feedback * dp_feedback / (P[j].Mass * All.EtaKineticEnergy);
                      SphP[j].Energy += (1.0 - All.EtaKineticEnergy) * de_feedback;

                      /* momentum is injected in radial direction */
                      SphP[j].Momentum[0] += inj_mom[0];
                      SphP[j].Momentum[1] += inj_mom[1];
                      SphP[j].Momentum[2] += inj_mom[2];

                      /* Compute new kinetic energy and add to the internal energy */
                      EKin = 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] + SphP[j].Momentum[2] * SphP[j].Momentum[2]) / P[j].Mass;
                      SphP[j].Energy += EKin;

                      //SphP[j].IsKicked = 1;

                      deltaEKin += EKin;
                      sum_deltaEKin += deltaEKin;

                      sum_deltaMomentum[0] += inj_mom[0];
                      sum_deltaMomentum[1] += inj_mom[1];
                      sum_deltaMomentum[2] += inj_mom[2];

                      mass_kicked += P[j].Mass;

                      n_cells_kicked++;

#ifdef OUTPUT_STELLAR_FEEDBACK
                      /* save cumulative feedback energy for debugging */
                      SphP[j].TotEgyFeed += de_feedback;
                      SphP[j].IntEgyFeed += (1.0 - All.EtaKineticEnergy) * de_feedback;
                      SphP[j].KinEgyFeed += deltaEKin;      // CHECK ME!!!
#endif
                    }
                }
#ifdef USE_ENTROPY_FOR_COLD_FLOWS
              SphP[j].Utherm =
                (SphP[j].Energy -
                 0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] +
                        SphP[j].Momentum[2] * SphP[j].Momentum[2]) / P[j].Mass) / P[j].Mass / (All.cf_atime * All.cf_atime);
              SphP[j].A = (GAMMA - 1.0) * SphP[j].Utherm / pow(SphP[j].Density * All.cf_a3inv, GAMMA - 1);
              SphP[j].Entropy = log(SphP[j].A) * P[j].Mass;
#endif
            }

        }
    }

  out->kicked_cells = n_cells_kicked;
  out->mass_to_kick = prob * tot_mass;
  out->mass_kicked = mass_kicked;
  out->kick_vel = velocity;
  out->deltaEKin = sum_deltaEKin;
  out->deltaMomentum[0] = sum_deltaMomentum[0];
  out->deltaMomentum[1] = sum_deltaMomentum[1];
  out->deltaMomentum[2] = sum_deltaMomentum[2];
  out->TotalMomentumInjected = dp_tot;
}
#endif

#endif
