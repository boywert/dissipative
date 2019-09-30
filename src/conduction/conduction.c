/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/conduction.c
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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <gsl/gsl_math.h>

#include "../allvars.h"
#include "../proto.h"
#include "../domain.h"
#include "../voronoi.h"


/*! \file conduction.c
 *  \brief main driver for an anisotropic diffusion solver
 *
 *  This file contains the code for a monotonicity preserving anisotrpic diffusion solver. This implementation is based on the method outlined in Gao & Wu (2013) - Journal of Computational Physics 10/2013; 250:308-331. DOI: 10.1016/j.jcp.2013.05.013 .  
 */



#ifdef MONOTONE_CONDUCTION


//static double atime, hubble_a ;

#define MaxNoc 100
#define MaxTetra 10*MaxNoc
#define epsilon 1e-15
#define MINANGLE 1e-5
#define gamma 0.0


#define cm (All.HubbleParam/All.UnitLength_in_cm)
#define g  (All.HubbleParam/All.UnitMass_in_g)
#define s  (All.HubbleParam/All.UnitTime_in_s)
#define erg (g*cm*cm/(s*s))
#define keV (1.602e-9*erg)
#define deg 1.0
#define m_p (PROTONMASS * g)
#define k_B (BOLTZMANN * erg / deg)


static struct inexch
{
  double Kappa;
  double Utherm;
#ifdef CONDUCTION_ANISOTROPIC
  double B[DIMS];
#endif
}
 *InExch;

static struct fluxexch
{
  double Flux1[MaxNoc];
  double Flux2[MaxNoc];
  int ID[MaxNoc];
}
 *FluxExch;

#ifndef TWODIMS
static struct savedtetradata
{
  double y_sigma_p1[3];
  double y_sigma_p2[3];
  double y_sigma_p3[3];
  double mdls_ysigma_p1, mdls_ysigma_p2, mdls_ysigma_p3;
  double u1, u2, u3;
  int p1Noc, p2Noc, p3Noc;
}
 *datatetra;
#endif

static MyFloat *Kappa;

void monotone_conduction()
{
  mpi_printf("MONOTONE CONDUCTION: Doing Conduction For This Time Step\n");
  double dt;
  Kappa = (MyFloat *) mymalloc("Kappa", NumGas * sizeof(MyFloat));

  dt = get_time_step();
  if(dt > 0)
    {
      initialize_kappa();
#ifdef TWODIMS
      calculate_unidirectional_fluxes_2D();
#else
      calculate_unidirectional_fluxes_3D();
#endif
      calculate_incremental_u(dt);
    }
  myfree(Kappa);
}


double get_time_step()
{
  double dt;
  int idx, i;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i > 0)
        break;
    }

  dt = (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;
  return dt;
}


void initialize_kappa()
{
  int i;
  for(i = 0; i < NumGas; i++)
    Kappa[i] = All.ConductionCoeff / SphP[i].Volume;

  for(i = 0; i < Mesh.Nvf; i++)
    {
      Mesh.VF[i].Flux1p1 = 0.0;
      Mesh.VF[i].Flux2p1 = 0.0;
      Mesh.VF[i].Flux1p2 = 0.0;
      Mesh.VF[i].Flux2p2 = 0.0;
    }
}


void init_conductivity(void)
{
  All.ConductionCoeff = 1.0 * All.ConductionEfficiency; /* Only for test and only if CONDUCTION_CONSTANT IS ON!!! */

#ifndef CONDUCTION_CONSTANT
  double coulomb_log, meanweight;

  meanweight = m_p * 4.0 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));   /* assuming full ionization */
  coulomb_log = 37.8;           /* accordin1g to Sarazin's book */

  /* Kappa_Spitzer definition taken from Zakamska & Narayan 2003
   * ( ApJ 582:162-169, Eq. (5) )
   */
  All.ConductionCoeff = (1.84e-5 / coulomb_log * pow(meanweight / k_B * GAMMA_MINUS1, 2.5) * erg / (s * deg * cm));

  /* Note: Because we replace \nabla(T) in the conduction equation with
   * \nable(u), our conduction coefficient is not the usual kappa, but
   * rather kappa*(gamma-1)*mu/kB. We therefore need to multiply with
   * another factor of (meanweight / k_B * GAMMA_MINUS1).
   */

  printf("unit conv1 = %le\n", meanweight / k_B * GAMMA_MINUS1);
  All.ConductionCoeff *= meanweight / k_B * GAMMA_MINUS1;

  /* The conversion of ConductionCoeff between internal units and cgs
   * units involves one factor of 'h'. We take care of this here.
   */
  All.ConductionCoeff /= All.HubbleParam;

  /* Include the fraction factor f, which is called ConductionEfficiency,
   * of the classical Spitzer conductivity - see (ApJ 582:162-169, Eq.(5))
   */
  All.ConductionCoeff *= All.ConductionEfficiency;

#ifdef CONDUCTION_SATURATION
  All.ElectronFreePathFactor = 8 * pow(3.0, 1.5) * pow(GAMMA_MINUS1, 2) / pow(3 + 5 * HYDROGEN_MASSFRAC, 2)
    / (1 + HYDROGEN_MASSFRAC) / sqrt(M_PI) / coulomb_log * pow(PROTONMASS, 3) / pow(ELECTRONCHARGE, 4)
    / (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam) * pow(All.UnitPressure_in_cgs / All.UnitDensity_in_cgs, 2);

  /* If the above value is multiplied with u^2/rho in code units (with rho being the physical density), then
   * one gets the electrong mean free path in centimeter. Since we want to compare this with another length
   * scale in code units, we now add an additional factor to convert back to code units.
   */
  All.ElectronFreePathFactor *= All.HubbleParam / All.UnitLength_in_cm;
#endif
#endif /* CONDUCTION_CONSTANT */

  mpi_printf("CONDUCTION: Spitzer conductivity set to %g (int units)\n", All.ConductionCoeff);
}


#ifdef TWODIMS

void calculate_unidirectional_fluxes_2D()
{
  point *DP = Mesh.DP;
  face *VF = Mesh.VF;
  int i;

  double n[DIMS];

  int i1, j1, k1;

  int Noc;

  double lTnK[MaxNoc][DIMS], y_sigma[MaxNoc][DIMS];
  double lambda_K, lambda_L;

  double Lambda_K[DIMS][DIMS], Lambda_L[DIMS][DIMS];

  double Kappa_K, modB, B_K[DIMS], u_K, Kappa_L, B_L[DIMS], u_L, mdls, orthdist;
  double omega_K, omega_L, u[MaxNoc], zycoeff1, zycoeff2, area_kl[MaxNoc];


  double mdls_ltnk[MaxNoc];


  double mdls_ysigma[MaxNoc];
  int face_id[MaxNoc];

  double diffL[DIMS][DIMS];
  double dL[DIMS];


  InExch = (struct inexch *) mymalloc("InExch", Mesh_nimport * sizeof(struct inexch));

  conduction_exchange_vector();

  int dp = -1;
  int vf = -1;
  int particle = -1;

  for(i = 0; i < NumGas; i++)
    {
      Kappa_K = Kappa[i];

      if((SphP[i].B[0] == 0.0) && (SphP[i].B[1] == 0.0))
        {
          Lambda_K[0][0] = Kappa_K;
          Lambda_K[0][1] = 0.0;
          Lambda_K[1][0] = 0.0;
          Lambda_K[1][1] = Kappa_K;
        }
      else
        {

          modB = sqrt(square(SphP[i].B[0]) + square(SphP[i].B[1]));

          if(modB == 0.0)
            modB = 1.0;

          B_K[0] = SphP[i].B[0] / modB;
          B_K[1] = SphP[i].B[1] / modB;

          Lambda_K[0][0] = Kappa_K * B_K[0] * B_K[0];
          Lambda_K[0][1] = Kappa_K * B_K[0] * B_K[1];
          Lambda_K[1][0] = Kappa_K * B_K[1] * B_K[0];
          Lambda_K[1][1] = Kappa_K * B_K[1] * B_K[1];
        }
      u_K = SphP[i].Utherm;

      /*Do loop over neighbours of i */

      int q = SphP[i].first_connection;

      Noc = -1;

      while(q >= 0)
        {
          dp = DC[q].dp_index;
          vf = DC[q].vf_index;
          particle = Mesh.DP[dp].index;

          if(particle < 0)
            {
              q = DC[q].next;
              continue;
            }

          if(VF[vf].area > 1e-10 * SphP[i].SurfaceArea)
            {
              Noc++;
              face_id[Noc] = vf;

              if(Mesh.DP[dp].task == ThisTask)
                {

                  if(particle >= NumGas)
                    particle -= NumGas;

                  Kappa_L = Kappa[particle];

                  B_L[0] = SphP[particle].B[0];
                  B_L[1] = SphP[particle].B[1];

                  u_L = SphP[particle].Utherm;
                }
              else
                {
                  Kappa_L = InExch[particle].Kappa;

                  B_L[0] = InExch[particle].B[0];
                  B_L[1] = InExch[particle].B[1];

                  u_L = InExch[particle].Utherm;
                }


              if((B_L[0] == 0.0) && (B_L[1] == 0.0))
                {
                  Lambda_L[0][0] = Kappa_L;
                  Lambda_L[0][1] = 0.0;
                  Lambda_L[1][0] = 0.0;
                  Lambda_L[1][1] = Kappa_L;
                }
              else
                {

                  modB = sqrt(square(B_L[0]) + square(B_L[1]));

                  if(modB == 0.0)
                    modB = 1.0;

                  B_L[0] /= modB;
                  B_L[1] /= modB;

                  Lambda_L[0][0] = Kappa_L * B_L[0] * B_L[0];
                  Lambda_L[0][1] = Kappa_L * B_L[0] * B_L[1];
                  Lambda_L[1][0] = Kappa_L * B_L[1] * B_L[0];
                  Lambda_L[1][1] = Kappa_L * B_L[1] * B_L[1];
                }

              n[0] = DP[dp].x - P[i].Pos[0];
              n[1] = DP[dp].y - P[i].Pos[1];

              //rintf("n x = %le\t y = %le\n", n[0], n[1]) ;


              mdls = sqrt(n[0] * n[0] + n[1] * n[1]);


              n[0] /= mdls;
              n[1] /= mdls;


              lambda_K = square(B_K[0] * n[0] + B_K[1] * n[1]);
              lambda_L = square(B_L[0] * n[0] + B_L[1] * n[1]);



              int c1, c2;
              for(c1 = 0; c1 < 2; c1++)
                for(c2 = 0; c2 < 2; c2++)
                  diffL[c1][c2] = Lambda_K[c1][c2] - Lambda_L[c1][c2];

              dL[0] = diffL[0][0] * n[0] + diffL[0][1] * n[1];
              dL[1] = diffL[1][0] * n[0] + diffL[1][1] * n[1];


              if((lambda_K == 0.0) && (lambda_L == 0.0))
                {
                  omega_K = 0.5;
                  omega_L = 0.5;
                }
              else
                {
                  omega_K = lambda_K / (lambda_K + lambda_L);
                  omega_L = 1.0 - omega_K;
                }


              orthdist = mdls / 2.0;

              area_kl[Noc] = VF[vf].area;


              lTnK[Noc][0] = Lambda_K[0][0] * n[0] + Lambda_K[0][1] * n[1];
              lTnK[Noc][1] = Lambda_K[1][0] * n[0] + Lambda_K[1][1] * n[1];

              //printf("ltnk x = %le\t y = %le\n", lTnK[Noc][0], lTnK[Noc][1]) ;


              mdls_ltnk[Noc] = sqrt(square(lTnK[Noc][0]) + square(lTnK[Noc][1]));

              if(mdls_ltnk[Noc] == 0.0)
                mdls_ltnk[Noc] = 1.0;

              lTnK[Noc][0] /= mdls_ltnk[Noc];
              lTnK[Noc][1] /= mdls_ltnk[Noc];

              y_sigma[Noc][0] = omega_L * (DP[dp].x - P[i].Pos[0]);     //+ orthdist*dL[0]/(lambda_K+lambda_L);
              y_sigma[Noc][1] = omega_L * (DP[dp].y - P[i].Pos[1]);     //+ orthdist*dL[1]/(lambda_K+lambda_L); 

              mdls_ysigma[Noc] = sqrt(square(y_sigma[Noc][0]) + square(y_sigma[Noc][1]));

              y_sigma[Noc][0] /= mdls_ysigma[Noc];
              y_sigma[Noc][1] /= mdls_ysigma[Noc];


              u[Noc] = omega_K * u_K + omega_L * u_L;

            }
          if(q == SphP[i].last_connection)
            break;

          q = DC[q].next;

        }

      double a, b;
      double x1, y1, x2, y2, x, y;
      double denom;

      for(i1 = 0; i1 <= Noc; i1++)
        {

          x = lTnK[i1][0];
          y = lTnK[i1][1];

          int numberidx = -1;
          for(j1 = 0; j1 <= Noc; j1++)
            {

              x1 = y_sigma[j1][0];
              y1 = y_sigma[j1][1];

              for(k1 = 0; k1 <= Noc; k1++)
                {
                  if(k1 == j1)
                    continue;

                  x2 = y_sigma[k1][0];
                  y2 = y_sigma[k1][1];


                  denom = y2 * x1 - x2 * y1;
                  if(mod(denom) < MINANGLE)
                    {
                      a = -1.0;
                      b = -1.0;
                    }
                  else
                    {
                      a = (x * y2 - x2 * y) / (y2 * x1 - x2 * y1);
                      b = (y * x1 - x * y1) / (y2 * x1 - x2 * y1);
                    }

                  if(a >= 0 && b >= 0)
                    {
                      numberidx = 1;
                      break;
                    }
                }
              if(numberidx == 1)
                break;
            }

          if((a < 0.0) && (b < 0.0))
            {
              printf("ltnk x = %le\ty = %le\n", lTnK[i1][0], lTnK[i1][1]);
              printf("a = %le\tb=%le\n", a, b);
              terminate("|||| Fatal  Error one of a,b < 0\n");
            }

#ifdef WG15_INTPL
          zycoeff1 = (area_kl[i1] * mdls_ltnk[i1] * a / mdls_ysigma[j1] + area_kl[i1] * mdls_ltnk[i1] * b / mdls_ysigma[k1]) * u_K;
          zycoeff2 = area_kl[i1] * mdls_ltnk[i1] * a * u[j1] / mdls_ysigma[j1] + area_kl[i1] * mdls_ltnk[i1] * b * u[k1] / mdls_ysigma[k1];
#else
          zycoeff1 = area_kl[i1] * mdls_ltnk[i1] * a * (u_K - u[j1]) / mdls_ysigma[j1];
          zycoeff2 = area_kl[i1] * mdls_ltnk[i1] * b * (u_K - u[k1]) / mdls_ysigma[k1];
#endif
          int vfidx = face_id[i1];
          int to = VF[vfidx].p1;
          int index = DP[to].ID;

          double F1, F2;
#ifdef WG15_INTPL
          F1 = zycoeff1;
          F2 = zycoeff2;
#else
          if(i1 == j1)
            {
              F1 = (1.0 - gamma) * zycoeff1;
              F2 = gamma * zycoeff1 + zycoeff2;
            }
          else if(i1 == k1)
            {
              F1 = (1.0 - gamma) * zycoeff2;
              F2 = gamma * zycoeff2 + zycoeff1;
            }
          else
            {
              F1 = 0.0;
              F2 = zycoeff1 + zycoeff2;
            }
          //      if((zycoeff1*(u_K - u[j1]) < 0.0) || (zycoeff2 * (u_K - u[k1]) < 0.0))
          //terminate("|||| Fatal  Error one of a,b,c < 0\n") ;
#endif
          if(index == P[i].ID)
            {
              VF[vfidx].Flux1p1 = F1;
              VF[vfidx].Flux2p1 = F2;
            }
          else
            {
              VF[vfidx].Flux1p2 = F1;
              VF[vfidx].Flux2p2 = F2;
            }
        }
    }
  myfree(InExch);
}

#else

void calculate_unidirectional_fluxes_3D()
{
  point *DP = Mesh.DP;
  face *VF = Mesh.VF;
  tetra *DT = Mesh.DT;
  int i;


  double n[DIMS];
  int i1, j1;
  int Noc;

  double lTnK[MaxNoc][DIMS];

  const int edge_start[6] = { 0, 0, 0, 1, 1, 2 };
  const int edge_end[6] = { 1, 2, 3, 2, 3, 3 };
  const int edge_opposite[6] = { 3, 1, 2, 3, 0, 1 };
  const int edge_nexttetra[6] = { 2, 3, 1, 0, 2, 0 };

  double Lambda_K[DIMS][DIMS];

  double Kappa_K, modB, modBL, B_K[DIMS], u_K, Kappa_L, B_L[DIMS], u_L, mdls;
  double omega_K, omega_L, u, zycoeff1, zycoeff2, area_kl[MaxNoc];
  double zycoeff3;

  double mdls_ltnk[MaxNoc];

  int tetraindexes[MaxTetra];

  double lambda_K, lambda_L;

  int face_id[MaxNoc];

  double px1, py1, pz1;

  int particle_id[MaxNoc];

  InExch = (struct inexch *) mymalloc("InExch", Mesh_nimport * sizeof(struct inexch));

  datatetra = (struct savedtetradata *) mymalloc("InExch", MaxTetra * sizeof(struct savedtetradata));

  conduction_exchange_vector();

  int dp = -1;
  int vf = -1;
  int particle = -1;

  for(i = 0; i < NumGas; i++)
    {
      Kappa_K = Kappa[i];

      if((SphP[i].B[0] == 0.0) && (SphP[i].B[1] == 0.0) && (SphP[i].B[2] == 0.0))
        {
          Lambda_K[0][0] = Kappa_K;
          Lambda_K[0][1] = 0.0;
          Lambda_K[1][0] = 0.0;
          Lambda_K[1][1] = Kappa_K;
          Lambda_K[0][2] = 0.0;
          Lambda_K[1][2] = 0.0;
          Lambda_K[2][0] = 0.0;
          Lambda_K[2][1] = 0.0;
          Lambda_K[2][2] = Kappa_K;
        }
      else
        {
          B_K[0] = SphP[i].B[0];
          B_K[1] = SphP[i].B[1];
          B_K[2] = SphP[i].B[2];

          modB = sqrt(square(SphP[i].B[0]) + square(SphP[i].B[1]) + square(SphP[i].B[2]));

          if(modB == 0.0)
            modB = 1.0;

          B_K[0] /= modB;
          B_K[1] /= modB;
          B_K[2] /= modB;


          Lambda_K[0][0] = Kappa_K * B_K[0] * B_K[0];
          Lambda_K[0][1] = Kappa_K * B_K[0] * B_K[1];
          Lambda_K[0][2] = Kappa_K * B_K[0] * B_K[2];
          Lambda_K[1][0] = Kappa_K * B_K[1] * B_K[0];
          Lambda_K[1][1] = Kappa_K * B_K[1] * B_K[1];
          Lambda_K[1][2] = Kappa_K * B_K[1] * B_K[2];
          Lambda_K[2][0] = Kappa_K * B_K[2] * B_K[0];
          Lambda_K[2][1] = Kappa_K * B_K[2] * B_K[1];
          Lambda_K[2][2] = Kappa_K * B_K[2] * B_K[2];
        }

      u_K = SphP[i].Utherm;

      /*Do loop over neighbours of i */

      int q = SphP[i].first_connection;


      int tottetra = 0;
      Noc = -1;
      while(q >= 0)
        {
          dp = DC[q].dp_index;
          vf = DC[q].vf_index;
          particle = DP[dp].index;


          if(particle < 0)
            {
              q = DC[q].next;
              continue;
            }

          Noc++;

          face_id[Noc] = vf;

          particle_id[Noc] = particle;

          if(DP[dp].task == ThisTask)
            {

              if(particle >= NumGas)
                particle -= NumGas;

              Kappa_L = Kappa[particle];

              B_L[0] = SphP[particle].B[0];
              B_L[1] = SphP[particle].B[1];
              B_L[2] = SphP[particle].B[2];

              u_L = SphP[particle].Utherm;
            }
          else
            {
              Kappa_L = InExch[particle].Kappa;

              B_L[0] = InExch[particle].B[0];
              B_L[1] = InExch[particle].B[1];
              B_L[2] = InExch[particle].B[2];

              u_L = InExch[particle].Utherm;
            }

          modBL = sqrt(B_L[0] * B_L[0] + B_L[1] * B_L[1] + B_L[2] * B_L[2]);

          if(modBL == 0.0)
            modBL = 1.0;

          B_L[0] /= modBL;
          B_L[1] /= modBL;
          B_L[2] /= modBL;



          n[0] = DP[dp].x - P[i].Pos[0];
          n[1] = DP[dp].y - P[i].Pos[1];
          n[2] = DP[dp].z - P[i].Pos[2];

          mdls = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);


          n[0] /= mdls;
          n[1] /= mdls;
          n[2] /= mdls;


          lambda_K = Kappa_K * square(B_K[0] * n[0] + B_K[1] * n[1] + B_K[2] * n[2]);
          lambda_L = Kappa_L * square(B_L[0] * n[0] + B_L[1] * n[1] + B_L[2] * n[2]);


          if((lambda_K == 0.0) || (lambda_L == 0.0))
            {
              omega_K = 0.5;
              omega_L = 0.5;
            }
          else
            {
              omega_K = lambda_K / (lambda_K + lambda_L);
              omega_L = 1.0 - omega_K;
            }



          area_kl[Noc] = VF[vf].area;


          lTnK[Noc][0] = Lambda_K[0][0] * n[0] + Lambda_K[0][1] * n[1] + Lambda_K[0][2] * n[2];
          lTnK[Noc][1] = Lambda_K[1][0] * n[0] + Lambda_K[1][1] * n[1] + Lambda_K[1][2] * n[2];
          lTnK[Noc][2] = Lambda_K[2][0] * n[0] + Lambda_K[2][1] * n[1] + Lambda_K[2][2] * n[2];

          mdls_ltnk[Noc] = sqrt(square(lTnK[Noc][0]) + square(lTnK[Noc][1]) + square(lTnK[Noc][2]));


          if(mdls_ltnk[Noc] == 0.0)
            mdls_ltnk[Noc] = 1.0;


          lTnK[Noc][0] /= mdls_ltnk[Noc];
          lTnK[Noc][1] /= mdls_ltnk[Noc];
          lTnK[Noc][2] /= mdls_ltnk[Noc];

          u = omega_K * u_K + omega_L * u_L;


          px1 = omega_L * (DP[dp].x - P[i].Pos[0]);
          py1 = omega_L * (DP[dp].y - P[i].Pos[1]);
          pz1 = omega_L * (DP[dp].z - P[i].Pos[2]);

          mdls = sqrt(px1 * px1 + py1 * py1 + pz1 * pz1);
          if(mdls == 0.0)
            mdls = 1.0;

          px1 /= mdls;
          py1 /= mdls;
          pz1 /= mdls;


          int tt = VF[vf].dt_index;
          tetra *t = &DT[tt];

          //      int p1 = i ;
          //int p2 = particle_id[Noc] ;
          int p1, p2;
          if(DP[VF[vf].p1].ID == P[i].ID)
            {
              p1 = VF[vf].p1;
              p2 = VF[vf].p2;
            }
          else
            {
              p2 = VF[vf].p1;
              p1 = VF[vf].p2;
            }

          int nr;

          int idp1, idp2, idp3;
          int taskp2, taskp3;

          double px, py, pz;



          for(nr = 0; nr < 6; nr++)
            {
              int start_index = t->p[edge_start[nr]];
              int end_index = t->p[edge_end[nr]];

              if((start_index == p1 && end_index == p2) || (start_index == p2 && end_index == p1))
                break;
            }


          tetra *this, *next;
          int i2, j2, k2, l2, m2, ii2, jj2, kk2, ll2, nn2;

          i2 = edge_start[nr];
          j2 = edge_end[nr];
          k2 = edge_opposite[nr];
          l2 = edge_nexttetra[nr];

          this = t;

          do
            {
              nn2 = this->t[l2];
              next = &DT[nn2];

              int breakindex = 1;
              for(i1 = 0; i1 < tottetra; i1++)
                {
                  if(tt == tetraindexes[i1])
                    {
                      breakindex = -1;
                      break;
                    }
                }


              if(breakindex == 1)
                {
                  tetraindexes[tottetra] = tt;

                  datatetra[tottetra].y_sigma_p1[0] = px1;
                  datatetra[tottetra].y_sigma_p1[1] = py1;
                  datatetra[tottetra].y_sigma_p1[2] = pz1;

                  datatetra[tottetra].mdls_ysigma_p1 = mdls;
                  datatetra[tottetra].u1 = u;
                  idp1 = particle_id[Noc];
                  datatetra[tottetra].p1Noc = idp1;



                  idp2 = DP[this->p[l2]].index;
                  taskp2 = DP[this->p[l2]].task;
                  datatetra[tottetra].p2Noc = idp2;

                  if(taskp2 == ThisTask)
                    {

                      if(idp2 >= NumGas)
                        idp2 -= NumGas;

                      Kappa_L = Kappa[idp2];

                      B_L[0] = SphP[idp2].B[0];
                      B_L[1] = SphP[idp2].B[1];
                      B_L[2] = SphP[idp2].B[2];

                      u_L = SphP[idp2].Utherm;
                    }
                  else
                    {
                      Kappa_L = InExch[idp2].Kappa;

                      B_L[0] = InExch[idp2].B[0];
                      B_L[1] = InExch[idp2].B[1];
                      B_L[2] = InExch[idp2].B[2];

                      u_L = InExch[idp2].Utherm;
                    }

                  modBL = sqrt(B_L[0] * B_L[0] + B_L[1] * B_L[1] + B_L[2] * B_L[2]);


                  if(modBL == 0.0)
                    modBL = 1.0;


                  B_L[0] /= modBL;
                  B_L[1] /= modBL;
                  B_L[2] /= modBL;


                  lambda_K = Kappa_K * square(B_K[0] * n[0] + B_K[1] * n[1] + B_K[2] * n[2]);
                  lambda_L = Kappa_L * square(B_L[0] * n[0] + B_L[1] * n[1] + B_L[2] * n[2]);


                  if((lambda_K == 0.0) || (lambda_L == 0.0))
                    {
                      omega_K = 0.5;
                      omega_L = 0.5;
                    }
                  else
                    {
                      omega_K = lambda_K / (lambda_K + lambda_L);
                      omega_L = 1.0 - omega_K;
                    }



                  datatetra[tottetra].u2 = omega_K * u_K + omega_L * u_L;

                  px = omega_L * (DP[this->p[l2]].x - P[i].Pos[0]);
                  py = omega_L * (DP[this->p[l2]].y - P[i].Pos[1]);
                  pz = omega_L * (DP[this->p[l2]].z - P[i].Pos[2]);

                  datatetra[tottetra].mdls_ysigma_p2 = sqrt(px * px + py * py + pz * pz);
                  if(datatetra[tottetra].mdls_ysigma_p2 == 0.0)
                    datatetra[tottetra].mdls_ysigma_p2 = 1.0;

                  datatetra[tottetra].y_sigma_p2[0] = px / datatetra[tottetra].mdls_ysigma_p2;
                  datatetra[tottetra].y_sigma_p2[1] = py / datatetra[tottetra].mdls_ysigma_p2;
                  datatetra[tottetra].y_sigma_p2[2] = pz / datatetra[tottetra].mdls_ysigma_p2;


                  idp3 = DP[this->p[k2]].index;
                  taskp3 = DP[this->p[k2]].task;
                  datatetra[tottetra].p3Noc = idp3;

                  if(taskp3 == ThisTask)
                    {

                      if(idp3 >= NumGas)
                        idp3 -= NumGas;

                      Kappa_L = Kappa[idp3];

                      B_L[0] = SphP[idp3].B[0];
                      B_L[1] = SphP[idp3].B[1];
                      B_L[2] = SphP[idp3].B[2];

                      u_L = SphP[idp3].Utherm;
                    }
                  else
                    {
                      Kappa_L = InExch[idp3].Kappa;

                      B_L[0] = InExch[idp3].B[0];
                      B_L[1] = InExch[idp3].B[1];
                      B_L[2] = InExch[idp3].B[2];

                      u_L = InExch[idp3].Utherm;
                    }

                  modBL = sqrt(B_L[0] * B_L[0] + B_L[1] * B_L[1] + B_L[2] * B_L[2]);

                  if(modBL == 0.0)
                    modBL = 1.0;

                  B_L[0] /= modBL;
                  B_L[1] /= modBL;
                  B_L[2] /= modBL;


                  lambda_K = Kappa_K * square(B_K[0] * n[0] + B_K[1] * n[1] + B_K[2] * n[2]);
                  lambda_L = Kappa_L * square(B_L[0] * n[0] + B_L[1] * n[1] + B_L[2] * n[2]);


                  if((lambda_K == 0.0) || (lambda_L == 0.0))
                    {
                      omega_K = 0.5;
                      omega_L = 0.5;
                    }
                  else
                    {
                      omega_K = lambda_K / (lambda_K + lambda_L);
                      omega_L = 1.0 - omega_K;
                    }


                  datatetra[tottetra].u3 = omega_K * u_K + omega_L * u_L;


                  px = omega_L * (DP[this->p[k2]].x - P[i].Pos[0]);
                  py = omega_L * (DP[this->p[k2]].y - P[i].Pos[1]);
                  pz = omega_L * (DP[this->p[k2]].z - P[i].Pos[2]);

                  datatetra[tottetra].mdls_ysigma_p3 = sqrt(px * px + py * py + pz * pz);
                  if(datatetra[tottetra].mdls_ysigma_p3 == 0.0)
                    datatetra[tottetra].mdls_ysigma_p3 = 1.0;

                  datatetra[tottetra].y_sigma_p3[0] = px / datatetra[tottetra].mdls_ysigma_p3;
                  datatetra[tottetra].y_sigma_p3[1] = py / datatetra[tottetra].mdls_ysigma_p3;
                  datatetra[tottetra].y_sigma_p3[2] = pz / datatetra[tottetra].mdls_ysigma_p3;

                  tottetra++;
                }


              for(m2 = 0, ll2 = ii2 = jj2 = -1; m2 < 4; m2++)
                {
                  if(next->p[m2] == this->p[k2])
                    ll2 = m2;
                  if(next->p[m2] == this->p[i2])
                    ii2 = m2;
                  if(next->p[m2] == this->p[j2])
                    jj2 = m2;
                }

              if(ll2 < 0 || ii2 < 0 || jj2 < 0)
                terminate("inconsistency");


              kk2 = 6 - (ll2 + ii2 + jj2);

              this = next;
              tt = nn2;

              i2 = ii2;
              l2 = ll2;
              j2 = jj2;
              k2 = kk2;
            }
          while(next != t);

          if(q == SphP[i].last_connection)
            break;

          q = DC[q].next;
        }


      double a, b, c;
      double x1, y1, x2, y2, x3, y3, x, y, z, z1, z2, z3;


      for(i1 = 0; i1 <= Noc; i1++)
        {
          if(VF[face_id[i1]].area > 1e-10 * SphP[i].SurfaceArea)
            {

              x = lTnK[i1][0];
              y = lTnK[i1][1];
              z = lTnK[i1][2];


              for(j1 = 0; j1 < tottetra; j1++)
                {

                  x1 = datatetra[j1].y_sigma_p1[0];
                  y1 = datatetra[j1].y_sigma_p1[1];
                  z1 = datatetra[j1].y_sigma_p1[2];

                  x2 = datatetra[j1].y_sigma_p2[0];
                  y2 = datatetra[j1].y_sigma_p2[1];
                  z2 = datatetra[j1].y_sigma_p2[2];

                  x3 = datatetra[j1].y_sigma_p3[0];
                  y3 = datatetra[j1].y_sigma_p3[1];
                  z3 = datatetra[j1].y_sigma_p3[2];

                  double denom = (z1 * (x2 * y3 - y2 * x3) + z2 * (x3 * y1 - y3 * x1) + z3 * (y2 * x1 - x2 * y1));

                  if(mod(denom) < MINANGLE)
                    {
                      c = -1.0;
                    }
                  else
                    {
                      a = (z * (y3 * x2 - y2 * x3) + z2 * (y * x3 - x * y3) + z3 * (y2 * x - y * x2)) / denom;
                      b = (z * (y1 * x3 - x1 * y3) + z1 * (x * y3 - y * x3) + z3 * (y * x1 - x * y1)) / denom;
                      c = (z * (x1 * y2 - x2 * y1) + z1 * (x2 * y - x * y2) + z2 * (y1 * x - x1 * y)) / denom;
                    }

                  if(a > -MINANGLE && b > -MINANGLE && c > -MINANGLE)
                    break;
                }

              if(mod(a) < MINANGLE)
                a = 0.0;
              if(mod(b) < MINANGLE)
                b = 0.0;
              if(mod(c) < MINANGLE)
                c = 0.0;

              if(a < 0.0 || b < 0.0 || c < 0.0)
                {
                  printf("ltnk x = %le\ty=%le\tz=%le\n", x, y, z);
                  printf("x1 x = %le\ty=%le\tz=%le\n", x1, y1, z1);
                  printf("x2 x = %le\ty=%le\tz=%le\n", x2, y2, z2);
                  printf("x3 x = %le\ty=%le\tz=%le\n", x3, y3, z3);
                  printf("a=%le\tb=%le\te=%le\n", a, b, c);
                  terminate("Fatal Error : One of a,b,c < 0.0");
                }


              zycoeff1 = area_kl[i1] * mdls_ltnk[i1] * a * (u_K - datatetra[j1].u1) / datatetra[j1].mdls_ysigma_p1;
              zycoeff2 = area_kl[i1] * mdls_ltnk[i1] * b * (u_K - datatetra[j1].u2) / datatetra[j1].mdls_ysigma_p2;
              zycoeff3 = area_kl[i1] * mdls_ltnk[i1] * c * (u_K - datatetra[j1].u3) / datatetra[j1].mdls_ysigma_p3;


              int vfidx = face_id[i1];
              int to = VF[vfidx].p1;
              int index = DP[to].ID;

              double F1, F2;

              if(particle_id[i1] == datatetra[j1].p1Noc)
                {
                  F1 = zycoeff1;
                  F2 = zycoeff2 + zycoeff3;
                }
              else if(particle_id[i1] == datatetra[j1].p2Noc)
                {
                  F1 = zycoeff2;
                  F2 = zycoeff1 + zycoeff3;
                }
              else if(particle_id[i1] == datatetra[j1].p3Noc)
                {
                  F1 = zycoeff3;
                  F2 = zycoeff1 + zycoeff2;
                }
              else
                {
                  F1 = 0.0;
                  F2 = zycoeff1 + zycoeff2 + zycoeff3;
                }

              if((zycoeff1 * (u_K - datatetra[j1].u1) < 0.0) || (zycoeff2 * (u_K - datatetra[j1].u2) < 0.0) || (zycoeff3 * (u_K - datatetra[j1].u3) < 0.0))
                terminate("|||| Fatal  Error one of a,b,c < 0\n\n");


              if(index == P[i].ID)
                {
                  VF[vfidx].Flux1p1 = F1;
                  VF[vfidx].Flux2p1 = F2;
                }
              else
                {
                  VF[vfidx].Flux1p2 = F1;
                  VF[vfidx].Flux2p2 = F2;
                }
            }
        }
    }
  myfree(datatetra);
  myfree(InExch);
}

#endif


void calculate_incremental_u(double dt)
{
  int i;

  point *DP = Mesh.DP;
  face *VF = Mesh.VF;

  double F1_K, F2_K, F1_L, F2_L, mu_K, mu_L, mup_K, mup_L;

  double u_K, u_L;

  FluxExch = (struct fluxexch *) mymalloc("fluxExch", Mesh_nimport * sizeof(struct fluxexch));

  flux_exchange_vector();

  InExch = (struct inexch *) mymalloc("InExch", Mesh_nimport * sizeof(struct inexch));

  conduction_exchange_vector();


  for(i = 0; i < NumGas; i++)
    {
      double Ftot = 0.0;
      int q = SphP[i].first_connection;

      u_K = SphP[i].Utherm;
      while(q >= 0)
        {
          int dp = DC[q].dp_index;
          int vf = DC[q].vf_index;
          int particle = Mesh.DP[dp].index;

          if(particle < 0)
            {
              q = DC[q].next;
              continue;
            }


          if(VF[vf].area > 1e-10 * SphP[i].SurfaceArea)
            {


              int to = VF[vf].p1;
              int index = DP[to].ID;



              if(index == P[i].ID)
                {
                  F1_K = Mesh.VF[vf].Flux1p1;
                  F2_K = Mesh.VF[vf].Flux2p1;
                }
              else
                {
                  F1_K = Mesh.VF[vf].Flux1p2;
                  F2_K = Mesh.VF[vf].Flux2p2;
                }


              if(Mesh.DP[dp].task == ThisTask)
                {
                  if(particle < NumGas)
                    {
                      u_L = SphP[particle].Utherm;
                      if(index == P[i].ID)
                        {
                          F1_L = Mesh.VF[vf].Flux1p2;
                          F2_L = Mesh.VF[vf].Flux2p2;
                        }
                      else
                        {
                          F1_L = Mesh.VF[vf].Flux1p1;
                          F2_L = Mesh.VF[vf].Flux2p1;
                        }
                    }
                  else
                    {
                      particle -= NumGas;

                      u_L = SphP[particle].Utherm;

                      int q_L = SphP[particle].first_connection;

                      while(q_L >= 0)
                        {
                          int dp_L = DC[q_L].dp_index;
                          int vf_L = DC[q_L].vf_index;
                          int particle_L = Mesh.DP[dp_L].index;

                          if(particle_L < 0)
                            {
                              q_L = DC[q_L].next;
                              continue;
                            }



                          int to_L = VF[vf_L].p1;
                          int index_L = DP[to_L].ID;

                          if(Mesh.DP[dp_L].task == ThisTask)
                            {
                              if(particle_L >= NumGas)
                                particle_L -= NumGas;

                              if(particle_L == i)
                                {
                                  if(index_L == P[particle].ID)
                                    {
                                      F1_L = Mesh.VF[vf_L].Flux1p1;
                                      F2_L = Mesh.VF[vf_L].Flux2p1;
                                    }
                                  else
                                    {
                                      F1_L = Mesh.VF[vf_L].Flux1p2;
                                      F2_L = Mesh.VF[vf_L].Flux2p2;
                                    }
                                  break;
                                }
                            }
                          if(q_L == SphP[particle].last_connection)
                            break;

                          q_L = DC[q_L].next;
                        }
                    }
                }
              else
                {
                  int number;
                  u_L = InExch[particle].Utherm;
                  for(number = 0; number < MaxNoc; number++)
                    {
                      if(FluxExch[particle].ID[number] == P[i].ID)
                        break;
                    }


                  F1_L = FluxExch[particle].Flux1[number];
                  F2_L = FluxExch[particle].Flux2[number];


                }
#ifdef WG15_INTPL
              if((F2_K == 0.0) && (F2_L == 0.0))
                {
                  mu_K = 0.5;
                  mu_L = 0.5;
                }
              else
                {
                  mu_K = mod(F2_L) / (mod(F2_L) + mod(F2_K));
                  mu_L = 1.0 - mu_K;
                }

              double B_sigma, B_sigma_plus, B_sigma_minus;

              B_sigma = mu_L * F2_L - mu_K * F2_K;
              B_sigma_plus = (mod(B_sigma) + B_sigma) / 2.0;
              B_sigma_minus = (mod(B_sigma) - B_sigma) / 2.0;

              Ftot += mu_K * F1_K - mu_L * F1_L + mu_K * (1.0 - sign(F2_L * F2_K)) * MC(B_sigma_plus, B_sigma_minus);
              //Ftot += mu_K*F1_K - mu_L*F1_L - mu_K*(1.0 - sign(F2_L*F2_K))*F2_K  + mu_L*(1.0 - sign(F2_L*F2_K))*F2_L;
#else
              mu_K = (mod(F2_L) + epsilon) / (mod(F2_L) + mod(F2_K) + 2.0 * epsilon);
              mu_L = (mod(F2_K) + epsilon) / (mod(F2_L) + mod(F2_K) + 2.0 * epsilon);
              mup_K = mod(F2_L) / (mod(F2_L) + mod(F2_K) + 2.0 * epsilon);
              mup_L = mod(F2_K) / (mod(F2_L) + mod(F2_K) + 2.0 * epsilon);

              Ftot += mu_K * F1_K - mu_L * F1_L + mup_K * (1.0 - sign(F2_L * F2_K)) * F2_K;
              //Ftot += mu_K*F1_K - mu_L*F1_L +  mup_K*(1.0 - sign(F2_L*F2_K))*MC(F2_K,F2_L) ;        
#endif
              //Ftot += mu_K*F1_K - mu_L*F1_L + mup_K*(1.0 - sign(F2_L*F2_K))*MC(F2_K,F2_L) ;
              //Ftot += mu_K*F1_K - mu_L*F1_L + mup_K*MC(F2_K,F2_L) ; 
            }

          if(q == SphP[i].last_connection)
            break;

          q = DC[q].next;

        }
      SphP[i].Energy += -dt * Ftot * P[i].Mass;
    }
  myfree(InExch);
  myfree(FluxExch);
}



void conduction_exchange_vector()
{
  int listp;
  int j, p, task, off;
  int ngrp, recvTask, place;    /* there was earlier the following defined: sendTask */


  struct inexch *tmpInExch;
  tmpInExch = (struct inexch *) mymalloc("tmpInExch", Mesh_nexport * sizeof(struct inexch));

  /* prepare data for export */
  for(j = 0; j < NTask; j++)
    Mesh_Send_count[j] = 0;

  for(p = 0; p < NumGas; p++)
    {
      if(P[p].Type == 0)
        {
          listp = List_P[p].firstexport;

          while(listp >= 0)
            {
              if((task = ListExports[listp].origin) != ThisTask)
                {
                  place = ListExports[listp].index;
                  off = Mesh_Send_offset[task] + Mesh_Send_count[task]++;

                  tmpInExch[off].Kappa = Kappa[place];
                  tmpInExch[off].Utherm = SphP[place].Utherm;
#ifdef CONDUCTION_ANISOTROPIC
                  tmpInExch[off].B[0] = SphP[place].B[0];
                  tmpInExch[off].B[1] = SphP[place].B[1];
                  tmpInExch[off].B[2] = SphP[place].B[2];
#endif
                }

              listp = ListExports[listp].nextexport;
            }
        }
    }

  /* exchange data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Mesh_Send_count[recvTask] > 0 || Mesh_Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&tmpInExch[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask] *
                           sizeof(struct inexch), MPI_BYTE, recvTask, TAG_DENS_A, &InExch[Mesh_Recv_offset[recvTask]],
                           Mesh_Recv_count[recvTask] * sizeof(struct inexch), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  myfree(tmpInExch);
}


void flux_exchange_vector()
{
  int listp;
  int j, p, task, off;
  int ngrp, recvTask, place;    /* there was earlier the following defined: sendTask */

  point *DP = Mesh.DP;
  face *VF = Mesh.VF;


  struct fluxexch *tmpFluxExch;
  tmpFluxExch = (struct fluxexch *) mymalloc("tmpFluxExch", Mesh_nexport * sizeof(struct fluxexch));

  /* prepare data for export */
  for(j = 0; j < NTask; j++)
    Mesh_Send_count[j] = 0;

  for(p = 0; p < NumGas; p++)
    {
      if(P[p].Type == 0)
        {
          listp = List_P[p].firstexport;

          while(listp >= 0)
            {
              if((task = ListExports[listp].origin) != ThisTask)
                {
                  place = ListExports[listp].index;
                  off = Mesh_Send_offset[task] + Mesh_Send_count[task]++;

                  int init;

                  for(init = 0; init < MaxNoc; init++)
                    {
                      tmpFluxExch[off].Flux1[init] = 0.0;
                      tmpFluxExch[off].Flux2[init] = 0.0;
                      tmpFluxExch[off].ID[init] = -1;
                    }


                  int q = SphP[place].first_connection;

                  int Noc = -1;

                  while(q >= 0)
                    {
                      int dp = DC[q].dp_index;
                      int vf = DC[q].vf_index;
                      int particle = Mesh.DP[dp].index;


                      int nto = VF[vf].p1;
                      int nindex = DP[nto].ID;


                      if(particle < 0)
                        {
                          q = DC[q].next;
                          continue;
                        }

                      if(VF[vf].area > 1e-10 * SphP[place].SurfaceArea)
                        {


                          Noc++;
                          if(DP[dp].task != ThisTask)
                            {
                              tmpFluxExch[off].ID[Noc] = DP[dp].ID;
                              if(nindex == P[place].ID)
                                {
                                  tmpFluxExch[off].Flux1[Noc] = VF[vf].Flux1p1;
                                  tmpFluxExch[off].Flux2[Noc] = VF[vf].Flux2p1;
                                }
                              else
                                {
                                  tmpFluxExch[off].Flux1[Noc] = VF[vf].Flux1p2;
                                  tmpFluxExch[off].Flux2[Noc] = VF[vf].Flux2p2;
                                }
                            }
                        }

                      if(q == SphP[place].last_connection)
                        break;

                      q = DC[q].next;
                    }

                }

              listp = ListExports[listp].nextexport;
            }


        }
    }

  /* exchange data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Mesh_Send_count[recvTask] > 0 || Mesh_Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&tmpFluxExch[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask] *
                           sizeof(struct fluxexch), MPI_BYTE, recvTask, TAG_DENS_A, &FluxExch[Mesh_Recv_offset[recvTask]],
                           Mesh_Recv_count[recvTask] * sizeof(struct fluxexch), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  myfree(tmpFluxExch);
}




double square(double x)
{
  return x * x;
}

double mod(double x)
{
  if(x < 0)
    return -x;
  else
    return x;
}

double MC(double a, double b)
{
  return minmod(2.0 * minmod(a, b), (a + b) / 2.0);
}

double minmod(double a, double b)
{
  if((a > 0.0) && (b > 0.0))
    return min(a, b);
  else if((a < 0.0) && (b < 0.0))
    return max(a, b);
  else
    return 0;
}

double min(double a, double b)
{
  if(a <= b)
    return a;
  else
    return b;
}

double max(double a, double b)
{
  if(a >= b)
    return a;
  else
    return b;
}

double sign(double x)
{
  if(x > 0)
    return 1.0;
  else if(x == 0.0)
    return 0.0;
  else
    return -1.0;
}

#endif /*MONOTONE_CONDUCTION */
