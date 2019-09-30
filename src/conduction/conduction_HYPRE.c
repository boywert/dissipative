/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/conduction/conduction_HYPRE.c
 * \date        MM/YYYY
 * \author      Rahul Kannan
 * \brief        
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 */


/*! \file conduction_HYPRE.c
 *  \brief main driver for an anisotropic/isotropic diffusion solver
 *
 *  This file contains the code for a extremum preserving anisotrpic diffusion solver. This implementation is based on the method outlined in Gao & Wu (2013) - Journal of Computational Physics 10/2013; 250:308-331. DOI: 10.1016/j.jcp.2013.05.013 .  

Implicit Time Integration - Non-linear system ( M(U)U = F(U) ) solved using  Picard iterations with  GMRES and multigrid preconditioner (using HYPRE library) for the linear system.  
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





#ifdef MONOTONE_CONDUCTION
#ifdef IMPLICIT_TI

#ifdef TWODIMS
#define MaxNoc 100
#else
#define MaxNoc 1000
#endif

#define MaxTetra 10*MaxNoc
#define epsilon 1e-25
#define MINANGLE 1e-5
#define gamma 0.0

#define epsilon_non_min 1e-2
#define epsilon_non_max 1e-2

//#define epsilon_non_min 1e-3
//#define epsilon_non_max 1e-5


#ifdef CONDUCTION_ISOTROPIC
#define epsilon_u_err 1e-15
#else
#ifndef SEMI_IMPLICIT_TI
#define epsilon_u_err 1e-4
//#define epsilon_u_err 1e-8
#else
#define epsolin_u_err 1e-10
#endif
#endif

#ifdef CONDUCTION_ISOTROPIC
#define epsilon_lin 1e-12
#else
#ifndef SEMI_IMPLICIT_TI
#define epsilon_lin 5e-4
//#define epsilon_lin 1e-8
#else
#define epsilon_lin 1e-8
#endif
#endif

#define MAX_ITER 100

#ifndef UMC_CORRECTOR
#define PICARD_BETA 1.0
#endif

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
#ifdef CONDUCTION_ANISOTROPIC
  double B[DIMS];
#endif
}
 *InExch;

static struct fluxexch
{
  double alss[MaxNoc];
  MyIDType ID[MaxNoc];
}
 *FluxExch;

static struct uexch
{
  double Utherm;
}
 *UExch;

static struct f2exch
{
  double F2[MaxNoc];
}
 *F2Exch;


static struct conductionfacedata
{
  int dp_kidp1, dp_kidp2;
  int dp1_kidp1, dp1_kidp2;
  int dp2_kidp1, dp2_kidp2;
  int dp3_kidp1, dp3_kidp2;

  double dp_akssp1, dp_akssp2;
  double dp1_akssp1, dp1_akssp2;
  double dp2_akssp1, dp2_akssp2;
  double dp3_akssp1, dp3_akssp2;

  double dp_akssp1p2, dp_akssp2p1;

  double F2p1, F2p2, F2p1p2, F2p2p1;

}
 *ConductionFaceData;

#ifndef TWODIMS
static struct savedtetradata
{
  double y_sigma_p1[3];
  double y_sigma_p2[3];
  double y_sigma_p3[3];
  double mdls_ysigma_p1, mdls_ysigma_p2, mdls_ysigma_p3;
  double omega1, omega2, omega3;
  int p1Noc, p2Noc, p3Noc;
}
 *datatetra;
#else
static struct savedtetradata
{
  double y_sigma_p1[2];
  double y_sigma_p2[2];
  double mdls_ysigma_p1, mdls_ysigma_p2;
  double omega1, omega2;
  int p1Noc, p2Noc;
  int p1task, p2task;
}
 *datatetra;
#endif

static MyFloat *Kappa;
static MyFloat *U0;
static MyFloat *Uminus1;
static MyFloat u_error;
static double mgm1;


void monotone_conduction()
{
  mpi_printf("CONDUCTION: Doing Conduction For This Time Step\n");
  double dt;
  Kappa = (MyFloat *) mymalloc("Kappa", NumGas * sizeof(MyFloat));

  ConductionFaceData = (struct conductionfacedata *) mymalloc("ConductionFaceData", Mesh.Nvf * sizeof(struct conductionfacedata));

  U0 = (MyFloat *) mymalloc("U0", NumGas * sizeof(MyFloat));
  Uminus1 = (MyFloat *) mymalloc("Uminus1", NumGas * sizeof(MyFloat));


  initialize_kappa();

  InExch = (struct inexch *) mymalloc("InExch", Mesh_nimport * sizeof(struct inexch));

  conduction_exchange_vector();

  dt = get_time_step();

  if(dt > 0)
    {
#ifdef TWODIMS
      calculate_unidirectional_fluxes_2D();
#else
      calculate_unidirectional_fluxes_3D();
#endif
      nonlinear_iterations(dt);
      update_U_and_Energy();
    }
  myfree(InExch);
  myfree(Uminus1);
  myfree(U0);
  myfree(ConductionFaceData);
  myfree(Kappa);

  mpi_printf("CONDUCTION: DONE\n");
}



void update_U_and_Energy()
{
  int i;

  for(i = 0; i < NumGas; i++)
    {
      double du = SphP[i].Utherm - U0[i];       /* [Temp] = erg/g   */
      SphP[i].Utherm = U0[i];
      SphP[i].Energy += All.cf_atime * All.cf_atime * du * P[i].Mass;   /* [Energy] = erg     */
      set_pressure_of_cell(i);
    }
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
  //#else
  //dt = (All.Conduction_Ti_endstep - All.Conduction_Ti_begstep) * All.Timebase_interval;
  //dt *= All.cf_atime / All.cf_time_hubble_a;
  //#endif
  return dt;
}


void initialize_kappa()
{
  int i;

  for(i = 0; i < NumGas; i++)
    {
      Kappa[i] = All.ConductionCoeff / SphP[i].Volume;
      U0[i] = SphP[i].Utherm;
      Uminus1[i] = SphP[i].Utherm;
    }

  for(i = 0; i < Mesh.Nvf; i++)
    {
      ConductionFaceData[i].dp_kidp1 = -1;
      ConductionFaceData[i].dp1_kidp1 = -1;
      ConductionFaceData[i].dp2_kidp1 = -1;
      ConductionFaceData[i].dp3_kidp1 = -1;

      ConductionFaceData[i].dp_kidp2 = -1;
      ConductionFaceData[i].dp1_kidp2 = -1;
      ConductionFaceData[i].dp2_kidp2 = -1;
      ConductionFaceData[i].dp3_kidp2 = -1;

      ConductionFaceData[i].dp_akssp1 = 0.0;
      ConductionFaceData[i].dp1_akssp1 = 0.0;
      ConductionFaceData[i].dp2_akssp1 = 0.0;
      ConductionFaceData[i].dp3_akssp1 = 0.0;

      ConductionFaceData[i].dp_akssp2 = 0.0;
      ConductionFaceData[i].dp1_akssp2 = 0.0;
      ConductionFaceData[i].dp2_akssp2 = 0.0;
      ConductionFaceData[i].dp3_akssp2 = 0.0;

      ConductionFaceData[i].dp_akssp1p2 = 0.0;

      ConductionFaceData[i].dp_akssp2p1 = 0.0;

      ConductionFaceData[i].F2p1 = 0.0;
      ConductionFaceData[i].F2p2 = 0.0;
      ConductionFaceData[i].F2p1p2 = 0.0;
      ConductionFaceData[i].F2p2p1 = 0.0;

    }
}



void init_conductivity(void)
{
  All.ConductionCoeff = 1.0 * All.ConductionEfficiency; /* Only for test and only if CONDUCTION_CONSTANT IS ON!!! */

  /*  double meanweight;

     meanweight = m_p * 4.0 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));   /* assuming full ionization */
  /*All.ConductionCoeff *= meanweight / k_B * GAMMA_MINUS1;
   */
  mpi_printf("MONOTONE CONDUCTION: Conduction Coeff = %le\n", All.ConductionCoeff);

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
  tetra *DT = Mesh.DT;
  int i;

  double n[DIMS];

  int i1, j1;

  int Noc;

  double lTnK[MaxNoc][DIMS];
  double lambda_K, lambda_L;

  double Lambda_K[DIMS][DIMS], Lambda_L[DIMS][DIMS];

  double Kappa_K, modB, B_K[DIMS], Kappa_L, B_L[DIMS], mdls, orthdist;
  double omega_K, omega_L, area_kl[MaxNoc];


  int tetraindexes[MaxTetra];

  double mdls_ltnk[MaxNoc];

  int taskp1, taskp2;

  double modBL;

  int particle_id[MaxNoc];
  int face_id[MaxNoc];

  int lktask[MaxNoc];


  const int edge_start[3] = { 1, 2, 0 };
  const int edge_end[3] = { 2, 0, 1 };
  const int edge_nexttetra[3] = { 0, 1, 2 };



  datatetra = (struct savedtetradata *) mymalloc("datatetra", MaxTetra * sizeof(struct savedtetradata));


  int tottetra;

  int dp = -1;
  int vf = -1;
  int particle = -1;

  for(i = 0; i < NumGas; i++)
    {

      Kappa_K = Kappa[i];

#ifdef CONDUCTION_ISOTROPIC
      Lambda_K[0][0] = Kappa_K;
      Lambda_K[0][1] = 0.0;
      Lambda_K[1][0] = 0.0;
      Lambda_K[1][1] = Kappa_K;
#endif

#ifdef CONDUCTION_ANISOTROPIC
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
#endif
      /*Do loop over neighbours of i */

      int q = SphP[i].first_connection;

      tottetra = 0;
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

          Noc++;
          face_id[Noc] = vf;
          particle_id[Noc] = particle;
          lktask[Noc] = Mesh.DP[dp].task;

          if(DP[dp].task == ThisTask)
            {
              if(particle >= NumGas)
                particle -= NumGas;

              Kappa_L = Kappa[particle];
#ifdef CONDUCTION_ANISOTROPIC
              B_L[0] = SphP[particle].B[0];
              B_L[1] = SphP[particle].B[1];
#endif
            }
          else
            {
              Kappa_L = InExch[particle].Kappa;
#ifdef CONDUCTION_ANISOTROPIC
              B_L[0] = InExch[particle].B[0];
              B_L[1] = InExch[particle].B[1];
#endif
            }

#ifdef CONDUCTION_ISOTROPIC
          Lambda_L[0][0] = Kappa_L;
          Lambda_L[0][1] = 0.0;
          Lambda_L[1][0] = 0.0;
          Lambda_L[1][1] = Kappa_L;
#endif
#ifdef CONDUCTION_ANISOTROPIC
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
#endif

          n[0] = DP[dp].x - P[i].Pos[0];
          n[1] = DP[dp].y - P[i].Pos[1];

          mdls = sqrt(n[0] * n[0] + n[1] * n[1]);


          n[0] /= mdls;
          n[1] /= mdls;

#ifdef CONDUCTION_ISOTROPIC
          omega_K = 0.5;
          omega_L = 0.5;
#endif
#ifdef CONDUCTION_ANISOTROPIC
          lambda_K = Kappa_K * square(B_K[0] * n[0] + B_K[1] * n[1]);
          lambda_L = Kappa_L * square(B_L[0] * n[0] + B_L[1] * n[1]);

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

#endif
          orthdist = mdls / 2.0;

          area_kl[Noc] = VF[vf].area;


          lTnK[Noc][0] = Lambda_K[0][0] * n[0] + Lambda_K[0][1] * n[1];
          lTnK[Noc][1] = Lambda_K[1][0] * n[0] + Lambda_K[1][1] * n[1];


          mdls_ltnk[Noc] = sqrt(square(lTnK[Noc][0]) + square(lTnK[Noc][1]));

          if(mdls_ltnk[Noc] == 0.0)
            mdls_ltnk[Noc] = 1.0;

          lTnK[Noc][0] /= mdls_ltnk[Noc];
          lTnK[Noc][1] /= mdls_ltnk[Noc];

          double px1, py1;

          px1 = omega_L * (DP[dp].x - P[i].Pos[0]);
          py1 = omega_L * (DP[dp].y - P[i].Pos[1]);


          mdls = sqrt(square(px1) + square(py1));

          px1 /= mdls;
          py1 /= mdls;


          int tt = VF[vf].dt_index;
          tetra *t = &DT[tt];

          //          int p1 = i ;
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

          int idp1, idp2;

          double px, py;

          for(nr = 0; nr < 3; nr++)
            {
              int start_index = t->p[edge_start[nr]];
              int end_index = t->p[edge_end[nr]];

              if((start_index == p1 && end_index == p2) || (start_index == p2 && end_index == p1))
                break;
            }


          tetra *this, *next;
          int i2, j2, l2, m2, ii2, jj2, ll2, nn2;

          int numbertetra = 1;

          i2 = edge_start[nr];
          j2 = edge_end[nr];
          l2 = edge_nexttetra[nr];

          this = t;
          do
            {
              nn2 = this->t[l2];
              next = &DT[nn2];

              numbertetra++;
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
                  datatetra[tottetra].omega1 = omega_K;
                  datatetra[tottetra].mdls_ysigma_p1 = mdls;

                  idp1 = particle_id[Noc];
                  datatetra[tottetra].p1Noc = idp1;
                  taskp1 = lktask[Noc];
                  datatetra[tottetra].p1task = taskp1;

                  idp2 = DP[this->p[l2]].index;
                  taskp2 = DP[this->p[l2]].task;

                  datatetra[tottetra].p2Noc = idp2;
                  datatetra[tottetra].p2task = taskp2;


                  if(taskp2 == ThisTask)
                    {

                      if(idp2 >= NumGas)
                        idp2 -= NumGas;

                      Kappa_L = Kappa[idp2];
#ifdef CONDUCTION_ANISOTROPIC
                      B_L[0] = SphP[idp2].B[0];
                      B_L[1] = SphP[idp2].B[1];
#endif
                    }
                  else
                    {
                      Kappa_L = InExch[idp2].Kappa;
#ifdef CONDUCTION_ANISOTROPIC
                      B_L[0] = InExch[idp2].B[0];
                      B_L[1] = InExch[idp2].B[1];
#endif
                    }

#ifdef CONDUCTION_ISOTROPIC
                  omega_K = 0.5;
                  omega_L = 0.5;
#endif
#ifdef CONDUCTION_ANISOTROPIC
                  modBL = sqrt(B_L[0] * B_L[0] + B_L[1] * B_L[1]);


                  B_L[0] /= modBL;
                  B_L[1] /= modBL;

                  if(modBL == 0.0)
                    modBL = 1.0;

                  lambda_K = Kappa_K * square(B_K[0] * n[0] + B_K[1] * n[1]);
                  lambda_L = Kappa_L * square(B_L[0] * n[0] + B_L[1] * n[1]);


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
#endif


                  px = omega_L * (DP[this->p[l2]].x - P[i].Pos[0]);
                  py = omega_L * (DP[this->p[l2]].y - P[i].Pos[1]);

                  datatetra[tottetra].omega2 = omega_K;
                  datatetra[tottetra].mdls_ysigma_p2 = sqrt(px * px + py * py);
                  if(datatetra[tottetra].mdls_ysigma_p2 == 0.0)
                    datatetra[tottetra].mdls_ysigma_p2 = 1.0;

                  datatetra[tottetra].y_sigma_p2[0] = px / datatetra[tottetra].mdls_ysigma_p2;
                  datatetra[tottetra].y_sigma_p2[1] = py / datatetra[tottetra].mdls_ysigma_p2;

                  tottetra++;
                }

              for(m2 = 0, ii2 = jj2 = -1; m2 < 3; m2++)
                {
                  if(next->p[m2] == this->p[j2])
                    jj2 = m2;
                  if(next->p[m2] == this->p[i2])
                    ii2 = m2;
                }

              if(ii2 < 0 || jj2 < 0)
                terminate("inconsistency");


              ll2 = 3 - (ii2 + jj2);

              this = next;
              tt = nn2;

              i2 = ii2;
              l2 = ll2;
              j2 = jj2;

            }
          while(numbertetra != 3);

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


          for(j1 = 0; j1 <= tottetra; j1++)
            {

              x1 = datatetra[j1].y_sigma_p1[0];
              y1 = datatetra[j1].y_sigma_p1[1];

              x2 = datatetra[j1].y_sigma_p2[0];
              y2 = datatetra[j1].y_sigma_p2[1];

              denom = y2 * x1 - x2 * y1;
              if(mod(denom) < MINANGLE)
                {
                  a = -1.0;
                  b = -1.0;
                }
              else
                {
                  a = (x * y2 - x2 * y) / denom;
                  b = (y * x1 - x * y1) / denom;
                }

              //              if(a>=0 && b>=0)
              //break ;

              if((a > -MINANGLE) && (b > -MINANGLE))
                break;
            }

          if(a >= -MINANGLE && a < 0.0)
            a = 0.0;

          if(b >= -MINANGLE && b < 0.0)
            b = 0.0;
          /*
             if(c>=-MINANGLE && c<0.0)
             c = 0.0 ;
           */


          if((a < 0.0) || (b < 0.0))
            {
              printf("ltnk x = %le\ty = %le\n", lTnK[i1][0], lTnK[i1][1]);
              printf("a = %le\tb=%le\n", a, b);
              terminate("|||| Fatal  Error one of a,b < 0\n");
            }


          double aKL1, aKL2;

          /*
             if(a<MINANGLE)
             a = 0.0 ;
             if(b<MINANGLE)
             b = 0.0 ;
           */
          aKL1 = area_kl[i1] * mdls_ltnk[i1] * a / datatetra[j1].mdls_ysigma_p1;
          aKL2 = area_kl[i1] * mdls_ltnk[i1] * b / datatetra[j1].mdls_ysigma_p2;



          int vfidx = face_id[i1];
          int to = VF[vfidx].p1;
          int index = DP[to].ID;

          int number1;


          double akss, aksp1, aksp2;

          if(particle_id[i1] == datatetra[j1].p1Noc)
            {
              akss = aKL1;
              aksp1 = aKL2;
              aksp2 = 0.0;
              number1 = 1;
            }
          else if(particle_id[i1] == datatetra[j1].p2Noc)
            {
              akss = aKL2;
              aksp1 = aKL1;
              aksp2 = 0.0;
              number1 = 2;
            }
          else
            {
              akss = 0.0;
              aksp1 = aKL1;
              aksp2 = aKL2;
              number1 = -1;
            }

          if(index == P[i].ID)
            {
              if(number1 == 1)
                {
                  ConductionFaceData[vfidx].Kid = datatetra[j1].p1Noc;
                  ConductionFaceData[vfidx].Kid1 = datatetra[j1].p2Noc;
                  ConductionFaceData[vfidx].Kid2 = -1;
                  ConductionFaceData[vfidx].Kid3 = -1;
                  ConductionFaceData[vfidx].akss = akss;
                  ConductionFaceData[vfidx].aksp1 = aksp1;
                  ConductionFaceData[vfidx].aksp2 = aksp2;
                  ConductionFaceData[vfidx].aksp3 = 0.0;
                  ConductionFaceData[vfidx].omegaKs = datatetra[j1].omega1;
                  ConductionFaceData[vfidx].omegaKsp1 = datatetra[j1].omega2;
                  ConductionFaceData[vfidx].omegaKsp2 = 0.0;
                  ConductionFaceData[vfidx].omegaKsp3 = 0.0;
                  ConductionFaceData[vfidx].ktask = datatetra[j1].p1task;
                  ConductionFaceData[vfidx].ktask1 = datatetra[j1].p2task;
                  ConductionFaceData[vfidx].ktask2 = -1;
                  ConductionFaceData[vfidx].ktask3 = -1;
                }
              else if(number1 == 2)
                {
                  ConductionFaceData[vfidx].Kid = datatetra[j1].p2Noc;
                  ConductionFaceData[vfidx].Kid1 = datatetra[j1].p1Noc;;
                  ConductionFaceData[vfidx].Kid2 = -1;
                  ConductionFaceData[vfidx].Kid3 = -1;
                  ConductionFaceData[vfidx].aksp1 = aksp1;
                  ConductionFaceData[vfidx].aksp2 = aksp2;
                  ConductionFaceData[vfidx].akss = akss;
                  ConductionFaceData[vfidx].aksp3 = 0.0;
                  ConductionFaceData[vfidx].omegaKsp1 = datatetra[j1].omega1;
                  ConductionFaceData[vfidx].omegaKsp2 = 0.0;
                  ConductionFaceData[vfidx].omegaKs = datatetra[j1].omega2;
                  ConductionFaceData[vfidx].omegaKsp3 = 0.0;
                  ConductionFaceData[vfidx].ktask1 = datatetra[j1].p1task;
                  ConductionFaceData[vfidx].ktask2 = -1;
                  ConductionFaceData[vfidx].ktask = datatetra[j1].p2task;
                  ConductionFaceData[vfidx].ktask3 = -1;
                }
              else
                {
                  ConductionFaceData[vfidx].Kid = -1;
                  ConductionFaceData[vfidx].Kid1 = datatetra[j1].p1Noc;
                  ConductionFaceData[vfidx].Kid2 = datatetra[j1].p2Noc;;
                  ConductionFaceData[vfidx].Kid3 = -1;
                  ConductionFaceData[vfidx].aksp1 = aksp1;
                  ConductionFaceData[vfidx].aksp2 = aksp2;
                  ConductionFaceData[vfidx].akss = akss;
                  ConductionFaceData[vfidx].aksp3 = 0.0;
                  ConductionFaceData[vfidx].omegaKsp1 = datatetra[j1].omega1;
                  ConductionFaceData[vfidx].omegaKsp2 = datatetra[j1].omega2;
                  ConductionFaceData[vfidx].omegaKs = 0.0;;
                  ConductionFaceData[vfidx].omegaKsp3 = 0.0;
                  ConductionFaceData[vfidx].ktask1 = datatetra[j1].p1task;
                  ConductionFaceData[vfidx].ktask2 = datatetra[j1].p2task;
                  ConductionFaceData[vfidx].ktask = -1;
                  ConductionFaceData[vfidx].ktask3 = -1;
                }
            }
          else
            {
              if(number1 == 1)
                {
                  ConductionFaceData[vfidx].Lid = datatetra[j1].p1Noc;
                  ConductionFaceData[vfidx].Lid1 = datatetra[j1].p2Noc;
                  ConductionFaceData[vfidx].Lid2 = -1;
                  ConductionFaceData[vfidx].Lid3 = -1;
                  ConductionFaceData[vfidx].alss = akss;
                  ConductionFaceData[vfidx].alsp1 = aksp1;
                  ConductionFaceData[vfidx].alsp2 = aksp2;
                  ConductionFaceData[vfidx].alsp3 = 0.0;
                  ConductionFaceData[vfidx].omegaLs = datatetra[j1].omega1;
                  ConductionFaceData[vfidx].omegaLsp1 = datatetra[j1].omega2;
                  ConductionFaceData[vfidx].omegaLsp2 = 0.0;
                  ConductionFaceData[vfidx].omegaLsp3 = 0.0;
                  ConductionFaceData[vfidx].ltask = datatetra[j1].p1task;
                  ConductionFaceData[vfidx].ltask1 = datatetra[j1].p2task;
                  ConductionFaceData[vfidx].ltask2 = -1;
                  ConductionFaceData[vfidx].ltask3 = -1;

                }
              else if(number1 == 2)
                {
                  ConductionFaceData[vfidx].Lid = datatetra[j1].p2Noc;
                  ConductionFaceData[vfidx].Lid1 = datatetra[j1].p1Noc;;
                  ConductionFaceData[vfidx].Lid2 = -1;
                  ConductionFaceData[vfidx].Lid3 = -1;
                  ConductionFaceData[vfidx].alsp1 = aksp1;
                  ConductionFaceData[vfidx].alsp2 = aksp2;
                  ConductionFaceData[vfidx].alss = akss;
                  ConductionFaceData[vfidx].alsp3 = 0.0;
                  ConductionFaceData[vfidx].omegaLsp1 = datatetra[j1].omega1;
                  ConductionFaceData[vfidx].omegaLsp2 = 0.0;
                  ConductionFaceData[vfidx].omegaLs = datatetra[j1].omega2;
                  ConductionFaceData[vfidx].omegaLsp3 = 0.0;
                  ConductionFaceData[vfidx].ltask1 = datatetra[j1].p1task;
                  ConductionFaceData[vfidx].ltask2 = -1;
                  ConductionFaceData[vfidx].ltask = datatetra[j1].p2task;
                  ConductionFaceData[vfidx].ltask3 = -1;
                }
              else
                {
                  ConductionFaceData[vfidx].Lid = -1;
                  ConductionFaceData[vfidx].Lid1 = datatetra[j1].p1Noc;
                  ConductionFaceData[vfidx].Lid2 = datatetra[j1].p2Noc;;
                  ConductionFaceData[vfidx].Lid3 = -1;
                  ConductionFaceData[vfidx].alsp1 = aksp1;
                  ConductionFaceData[vfidx].alsp2 = aksp2;
                  ConductionFaceData[vfidx].alss = akss;
                  ConductionFaceData[vfidx].alsp3 = 0.0;
                  ConductionFaceData[vfidx].omegaLsp1 = datatetra[j1].omega1;
                  ConductionFaceData[vfidx].omegaLsp2 = datatetra[j1].omega2;
                  ConductionFaceData[vfidx].omegaLs = 0.0;;
                  ConductionFaceData[vfidx].omegaLsp3 = 0.0;
                  ConductionFaceData[vfidx].ltask1 = datatetra[j1].p1task;
                  ConductionFaceData[vfidx].ltask2 = datatetra[j1].p2task;
                  ConductionFaceData[vfidx].ltask = -1;
                  ConductionFaceData[vfidx].ltask3 = -1;
                }
            }
        }
    }
  myfree(datatetra);
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



  double Kappa_K, modB, modBL, B_K[DIMS], Kappa_L, B_L[DIMS], mdls;
  double omega_K, omega_L, area_kl[MaxNoc];

  double mdls_ltnk[MaxNoc];

  int tetraindexes[MaxTetra];

  double lambda_K, lambda_L;

  int face_id[MaxNoc];

  double px1, py1, pz1;

  int particle_id[MaxNoc];

  datatetra = (struct savedtetradata *) mymalloc("datatetra", MaxTetra * sizeof(struct savedtetradata));


  int dp = -1;
  int vf = -1;
  int particle = -1;

  for(i = 0; i < NumGas; i++)
    {

      //      printf("i=%d\tx=%lf\ty=%lf\tz=%lf\tutherm=%lf\n", i, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], SphP[i].Utherm) ;

      Kappa_K = Kappa[i];

#ifdef CONDUCTION_ISOTROPIC
      Lambda_K[0][0] = Kappa_K;
      Lambda_K[0][1] = 0.0;
      Lambda_K[1][0] = 0.0;
      Lambda_K[1][1] = Kappa_K;
      Lambda_K[0][2] = 0.0;
      Lambda_K[1][2] = 0.0;
      Lambda_K[2][0] = 0.0;
      Lambda_K[2][1] = 0.0;
      Lambda_K[2][2] = Kappa_K;
#endif
#ifdef CONDUCTION_ANISOTROPIC
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
#endif
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

          particle_id[Noc] = dp;

          if(DP[dp].task == ThisTask)
            {

              if(particle >= NumGas)
                particle -= NumGas;

              Kappa_L = Kappa[particle];
#ifdef CONDUCTION_ANISOTROPIC
              B_L[0] = SphP[particle].B[0];
              B_L[1] = SphP[particle].B[1];
              B_L[2] = SphP[particle].B[2];
#endif
            }
          else
            {
              Kappa_L = InExch[particle].Kappa;
#ifdef CONDUCTION_ANISOTROPIC
              B_L[0] = InExch[particle].B[0];
              B_L[1] = InExch[particle].B[1];
              B_L[2] = InExch[particle].B[2];
#endif
            }

#ifdef CONDUCTION_ANISOTROPIC

          modBL = sqrt(B_L[0] * B_L[0] + B_L[1] * B_L[1] + B_L[2] * B_L[2]);

          if(modBL == 0.0)
            modBL = 1.0;

          B_L[0] /= modBL;
          B_L[1] /= modBL;
          B_L[2] /= modBL;

#endif

          n[0] = DP[dp].x - P[i].Pos[0];
          n[1] = DP[dp].y - P[i].Pos[1];
          n[2] = DP[dp].z - P[i].Pos[2];

          mdls = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);


          n[0] /= mdls;
          n[1] /= mdls;
          n[2] /= mdls;

#ifdef CONDUCTION_ISOTROPIC
          omega_K = 0.5;
          omega_L = 0.5;
#endif
#ifdef CONDUCTION_ANISOTROPIC
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
#endif


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

                  idp1 = particle_id[Noc];
                  datatetra[tottetra].p1Noc = idp1;
                  datatetra[tottetra].omega1 = omega_K;

                  idp2 = DP[this->p[l2]].index;
                  taskp2 = DP[this->p[l2]].task;
                  datatetra[tottetra].p2Noc = this->p[l2];

                  if(taskp2 == ThisTask)
                    {
                      if(idp2 >= NumGas)
                        idp2 -= NumGas;

                      Kappa_L = Kappa[idp2];
#ifdef CONDUCTION_ANISOTROPIC
                      B_L[0] = SphP[idp2].B[0];
                      B_L[1] = SphP[idp2].B[1];
                      B_L[2] = SphP[idp2].B[2];
#endif
                    }
                  else
                    {
                      Kappa_L = InExch[idp2].Kappa;
#ifdef CONDUCTION_ANISOTROPIC
                      B_L[0] = InExch[idp2].B[0];
                      B_L[1] = InExch[idp2].B[1];
                      B_L[2] = InExch[idp2].B[2];
#endif
                    }
#ifdef CONDUCTION_ANISOTROPIC
                  modBL = sqrt(B_L[0] * B_L[0] + B_L[1] * B_L[1] + B_L[2] * B_L[2]);


                  if(modBL == 0.0)
                    modBL = 1.0;

                  B_L[0] /= modBL;
                  B_L[1] /= modBL;
                  B_L[2] /= modBL;
#endif

#ifdef CONDUCTION_ISOTROPIC
                  omega_K = 0.5;
                  omega_L = 0.5;
#endif
#ifdef CONDUCTION_ANISOTROPIC
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

#endif
                  datatetra[tottetra].omega2 = omega_K;


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
                  datatetra[tottetra].p3Noc = this->p[k2];

                  if(taskp3 == ThisTask)
                    {
                      if(idp3 >= NumGas)
                        idp3 -= NumGas;

                      Kappa_L = Kappa[idp3];
#ifdef CONDUCTION_ANISOTROPIC
                      B_L[0] = SphP[idp3].B[0];
                      B_L[1] = SphP[idp3].B[1];
                      B_L[2] = SphP[idp3].B[2];
#endif
                    }
                  else
                    {
                      Kappa_L = InExch[idp3].Kappa;
#ifdef CONDUCTION_ANISOTROPIC
                      B_L[0] = InExch[idp3].B[0];
                      B_L[1] = InExch[idp3].B[1];
                      B_L[2] = InExch[idp3].B[2];
#endif
                    }
#ifdef CONDUCTION_ANISOTROPIC
                  modBL = sqrt(B_L[0] * B_L[0] + B_L[1] * B_L[1] + B_L[2] * B_L[2]);

                  if(modBL == 0.0)
                    modBL = 1.0;

                  B_L[0] /= modBL;
                  B_L[1] /= modBL;
                  B_L[2] /= modBL;
#endif

#ifdef CONDUCTION_ISOTROPIC
                  omega_K = 0.5;
                  omega_L = 0.5;
#endif
#ifdef CONDUCTION_ANISOTROPIC
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
#endif
                  datatetra[tottetra].omega3 = omega_K;


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

                  if((a > -MINANGLE) && (b > -MINANGLE) && (c > -MINANGLE))
                    break;
                }



              /*
                 if(i == 1)
                 {
                 printf("i = %d \t particle = %d\n", i, DP[particle_id[j1]].index) ;
                 printf("ltnk x = %le\t y = %le\t z = %le\n", x, y, z) ;
                 printf("p1 x = %le\t y = %le\t z = %le\n", x1, y1, z1) ;
                 printf("p2 x = %le\t y = %le\t z = %le\n", x2, y2, z2) ;
                 printf("p3 x = %le\t y = %le\t z = %le\n", x3, y3, z3) ;
                 printf("a,b,c = %le , %le , %le\n", a,b,c) ;
                 } 
               */

              if(a >= -MINANGLE && a < MINANGLE)
                a = 0.0;

              if(b >= -MINANGLE && b < MINANGLE)
                b = 0.0;

              if(c >= -MINANGLE && c < MINANGLE)
                c = 0.0;


              if(a < 0.0 || b < 0.0 || c < 0.0)
                {
                  printf("ID = %d\n", P[i].ID);
                  printf("NumGas = %d\n", NumGas);
                  printf("Particle Num = %d\n", i);
                  printf("task = %d\n", DP[particle_id[i1]].task);
                  printf("Noc=%d\ttottetra=%d\ti1=%d\tj1=%d\n", Noc, tottetra, i1, j1);
                  printf("x = %le\ty = %le\tz = %le\n", x, y, z);
                  printf("x1 = %le\ty1 = %le\tz1 = %le\n", x1, y1, z1);
                  printf("x2 = %le\ty2 = %le\tz2 = %le\n", x2, y2, z2);
                  printf("x3 = %le\ty3 = %le\tz3 = %le\n", x3, y3, z3);
                  printf("a=%le\tb=%le\tc=%le\n", a, b, c);
                  terminate("|||| Fatal  Error one of a,b,c < 0\n\n");
                }

              double aKL1, aKL2, aKL3;


              aKL1 = area_kl[i1] * mdls_ltnk[i1] * a / datatetra[j1].mdls_ysigma_p1;
              aKL2 = area_kl[i1] * mdls_ltnk[i1] * b / datatetra[j1].mdls_ysigma_p2;
              aKL3 = area_kl[i1] * mdls_ltnk[i1] * c / datatetra[j1].mdls_ysigma_p3;


              aKL1 = aKL1 * (1.0 - datatetra[j1].omega1);
              aKL2 = aKL2 * (1.0 - datatetra[j1].omega2);
              aKL3 = aKL3 * (1.0 - datatetra[j1].omega3);

              /*
                 if(i == 1)
                 {
                 printf("aKL1 = %le\t aKL2 = %le \t aKL3 = %le \n", aKL1, aKL2, aKL3) ;
                 printf("particle_id = %d \t datatetra j1 p1 = %d \t p2  = %d \t p3 = %d \n", particle_id[i1], datatetra[j1].p1Noc, datatetra[j1].p2Noc, datatetra[j1].p3Noc) ;
                 }
               */
              int vfidx = face_id[i1];
              int to = VF[vfidx].p1;
              int index = DP[to].ID;

              int dp_kid, dp1_kid, dp2_kid, dp3_kid;
              double dp_akss, dp1_akss, dp2_akss, dp3_akss;


              if(particle_id[i1] == datatetra[j1].p1Noc)
                {
                  dp_akss = aKL1;
                  dp1_akss = aKL2;
                  dp2_akss = aKL3;
                  dp3_akss = 0.0;
                  dp_kid = datatetra[j1].p1Noc;
                  dp1_kid = datatetra[j1].p2Noc;
                  dp2_kid = datatetra[j1].p3Noc;
                  dp3_kid = -1;
                }
              else if(particle_id[i1] == datatetra[j1].p2Noc)
                {
                  dp_akss = aKL2;
                  dp1_akss = aKL3;
                  dp2_akss = aKL1;
                  dp3_akss = 0.0;
                  dp_kid = datatetra[j1].p2Noc;
                  dp1_kid = datatetra[j1].p3Noc;
                  dp2_kid = datatetra[j1].p1Noc;
                  dp3_kid = -1;
                }
              else if(particle_id[i1] == datatetra[j1].p3Noc)
                {
                  dp_akss = aKL3;
                  dp1_akss = aKL1;
                  dp2_akss = aKL2;
                  dp3_akss = 0.0;
                  dp_kid = datatetra[j1].p3Noc;
                  dp1_kid = datatetra[j1].p1Noc;
                  dp2_kid = datatetra[j1].p2Noc;
                  dp3_kid = -1;
                }
              else
                {
                  dp_akss = 0.0;
                  dp1_akss = aKL1;
                  dp2_akss = aKL2;
                  dp3_akss = aKL3;
                  dp_kid = -1;
                  dp1_kid = datatetra[j1].p1Noc;
                  dp2_kid = datatetra[j1].p2Noc;
                  dp3_kid = datatetra[j1].p3Noc;
                }

              /*
                 if(i == 1)
                 {
                 printf("dp_kid = %d\t dp1_kid = %d\t dp2_kid = %d\t dp3_kid = %d\n", dp_kid, dp1_kid, dp2_kid, dp3_kid) ;
                 printf("dp_akss = %le\t dp1_akss = %le\t dp2_akss = %le\t dp3_akss = %le\n", dp_akss, dp1_akss, dp2_akss, dp3_akss) ;
                 }
               */

              if(index == P[i].ID)
                {
                  ConductionFaceData[vfidx].dp_kidp1 = dp_kid;
                  ConductionFaceData[vfidx].dp1_kidp1 = dp1_kid;
                  ConductionFaceData[vfidx].dp2_kidp1 = dp2_kid;
                  ConductionFaceData[vfidx].dp3_kidp1 = dp3_kid;

                  ConductionFaceData[vfidx].dp_akssp1 = dp_akss;
                  ConductionFaceData[vfidx].dp1_akssp1 = dp1_akss;
                  ConductionFaceData[vfidx].dp2_akssp1 = dp2_akss;
                  ConductionFaceData[vfidx].dp3_akssp1 = dp3_akss;
                }
              else
                {
                  ConductionFaceData[vfidx].dp_kidp2 = dp_kid;
                  ConductionFaceData[vfidx].dp1_kidp2 = dp1_kid;
                  ConductionFaceData[vfidx].dp2_kidp2 = dp2_kid;
                  ConductionFaceData[vfidx].dp3_kidp2 = dp3_kid;

                  ConductionFaceData[vfidx].dp_akssp2 = dp_akss;
                  ConductionFaceData[vfidx].dp1_akssp2 = dp1_akss;
                  ConductionFaceData[vfidx].dp2_akssp2 = dp2_akss;
                  ConductionFaceData[vfidx].dp3_akssp2 = dp3_akss;
                }
            }
        }
    }
  myfree(datatetra);
}

#endif

void nonlinear_iterations(double deltat)
{
  double error, error_init;

  FluxExch = (struct fluxexch *) mymalloc("FluxExch", Mesh_nimport * sizeof(struct fluxexch));
  flux_exchange_vector();

  double dt = deltat;

  int niter = 0;


  double errorminus1;

  error_init = calculate_nonlinear_error(dt);

#if defined(CONDUCTION_ANISOTROPIC)

  error = error_init;
  errorminus1 = 10.0 * error;


  mpi_printf("CONDUCTION: Initial error - non linear iterations = %le\n", error);

  while(niter < MAX_ITER)
    {
      if((error / error_init > epsilon_non_max) && (u_error > epsilon_u_err))
        {
          if(niter > 20 && (error / error_init < epsilon_non_min))
            break;

          //if(error > errorminus1) 
          //break ;

          errorminus1 = error;
          error = calculate_nonlinear_error(dt);
          mpi_printf("CONDUCTION: NON LINEAR ITERATION = %d \t relative error = %le\t initial error = %le\t u_error = %le\n", niter, error / error_init, error_init, u_error);
        }
      else
        break;

      niter++;
    }

  MPI_Barrier(MPI_COMM_WORLD);

  if(niter >= MAX_ITER)
    {
      if(error / error_init < epsilon_non_min)
        {
          mpi_printf("CONDUCTION: error/errorinit = %le\n", error / error_init);
          mpi_printf("CONDUCTION: Weak Convergence\n");
        }
      else
        {
          mpi_printf("CONDUCTION: error/errorinit = %le\n", error / error_init);
          mpi_printf("CONDUCTION: Non linear iterations failed to converge!\n");
        }
    }
#endif
  myfree(FluxExch);
}


double calculate_nonlinear_error(double dt)
{
  UExch = (struct uexch *) mymalloc("UExch", Mesh_nimport * sizeof(struct uexch));
  F2Exch = (struct f2exch *) mymalloc("F2Exch", Mesh_nimport * sizeof(struct f2exch));

  double error;
  MPI_Barrier(MPI_COMM_WORLD);
  U_exchange_vector();

  calculate_F2();
  MPI_Barrier(MPI_COMM_WORLD);
  F2_exchange_vector();


  set_initial_coeff();

  error = linear_iterations(dt);

  myfree(F2Exch);
  myfree(UExch);

  return error;
}



void calculate_F2()
{
  point *DP = Mesh.DP;
  face *VF = Mesh.VF;
  int i;

  int dp_kid, dp1_kid, dp2_kid, dp3_kid;
  double dp_akss, dp1_akss, dp2_akss, dp3_akss;

  double F2;
  double u_L1, u_L2, u_L3;

  for(i = 0; i < NumGas; i++)
    {
      int q = SphP[i].first_connection;
      double u_K = SphP[i].Utherm;

      while(q >= 0)
        {
          int dp = DC[q].dp_index;
          int vf = DC[q].vf_index;
          int particle = DP[dp].index;


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
                  dp_kid = ConductionFaceData[vf].dp_kidp1;
                  dp1_kid = ConductionFaceData[vf].dp1_kidp1;
                  dp2_kid = ConductionFaceData[vf].dp2_kidp1;
                  dp3_kid = ConductionFaceData[vf].dp3_kidp1;
                  dp_akss = ConductionFaceData[vf].dp_akssp1;
                  dp1_akss = ConductionFaceData[vf].dp1_akssp1;
                  dp2_akss = ConductionFaceData[vf].dp2_akssp1;
                  dp3_akss = ConductionFaceData[vf].dp3_akssp1;
                }
              else
                {
                  dp_kid = ConductionFaceData[vf].dp_kidp2;
                  dp1_kid = ConductionFaceData[vf].dp1_kidp2;
                  dp2_kid = ConductionFaceData[vf].dp2_kidp2;
                  dp3_kid = ConductionFaceData[vf].dp3_kidp2;
                  dp_akss = ConductionFaceData[vf].dp_akssp2;
                  dp1_akss = ConductionFaceData[vf].dp1_akssp2;
                  dp2_akss = ConductionFaceData[vf].dp2_akssp2;
                  dp3_akss = ConductionFaceData[vf].dp3_akssp2;
                }

              u_L1 = get_Utherm(dp1_kid);

              u_L2 = get_Utherm(dp2_kid);

              u_L3 = get_Utherm(dp3_kid);

              F2 = dp1_akss * (u_K - u_L1) + dp2_akss * (u_K - u_L2) + dp3_akss * (u_K - u_L3);

              if(index == P[i].ID)
                ConductionFaceData[vf].F2p1 = F2;
              else
                ConductionFaceData[vf].F2p2 = F2;
            }
          if(q == SphP[i].last_connection)
            break;

          q = DC[q].next;
        }
    }
}


void set_initial_coeff()
{

  int i;
  point *DP = Mesh.DP;
  face *VF = Mesh.VF;

  double u_K, u_L;

  double dp_akss;
  double dp_alss;
  int F2_K, F2_L;

  for(i = 0; i < NumGas; i++)
    {


      int q = SphP[i].first_connection;

      u_K = SphP[i].Utherm;


      while(q >= 0)
        {
          int dp = DC[q].dp_index;
          int vf = DC[q].vf_index;
          int particle = DP[dp].index;

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
                  F2_K = ConductionFaceData[vf].F2p1;
                  dp_akss = ConductionFaceData[vf].dp_akssp1;
                }
              else
                {
                  F2_K = ConductionFaceData[vf].F2p2;
                  dp_akss = ConductionFaceData[vf].dp_akssp2;
                }

              if(Mesh.DP[dp].task == ThisTask)
                {
                  if(particle < NumGas)
                    {
                      if(index == P[i].ID)
                        {
                          F2_L = ConductionFaceData[vf].F2p2;
                          dp_alss = ConductionFaceData[vf].dp_akssp2;
                        }
                      else
                        {
                          F2_L = ConductionFaceData[vf].F2p1;
                          dp_alss = ConductionFaceData[vf].dp_akssp1;
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
                                      F2_L = ConductionFaceData[vf_L].F2p1;
                                      dp_alss = ConductionFaceData[vf_L].dp_akssp1;
                                    }
                                  else
                                    {
                                      F2_L = ConductionFaceData[vf_L].F2p2;
                                      dp_alss = ConductionFaceData[vf_L].dp_akssp2;
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
                  for(number = 0; number < MaxNoc; number++)
                    {
                      if(FluxExch[particle].ID[number] == P[i].ID)
                        break;
                    }
                  F2_L = F2Exch[particle].F2[number];
                  dp_alss = FluxExch[particle].alss[number];
                }


              if(index == P[i].ID)
                {
                  ConductionFaceData[vf].F2p1p2 = F2_L;
                  ConductionFaceData[vf].dp_akssp1p2 = dp_alss;
                }
              else
                {
                  ConductionFaceData[vf].F2p2p1 = F2_L;
                  ConductionFaceData[vf].dp_akssp2p1 = dp_alss;
                }
            }
          if(q == SphP[i].last_connection)
            break;

          q = DC[q].next;
        }
    }
}



double get_Utherm(int dp)
{
  int index = Mesh.DP[dp].index;
  double utherm;
  if(dp == -1)
    {
      utherm = 0.0;
    }
  else
    {
      if(Mesh.DP[dp].task == ThisTask)
        {
          if(index >= NumGas)
            index -= NumGas;

          utherm = SphP[index].Utherm;
        }
      else
        utherm = UExch[index].Utherm;
    }
  return utherm;
}


double set_coeff(HYPRE_IJMatrix * A, HYPRE_IJVector * x, HYPRE_IJVector * b, int *fld_offsets, double dt)
{
  int i;
  point *DP = Mesh.DP;
  face *VF = Mesh.VF;

  int *rows = (int *) mymalloc("rows", sizeof(int) * NumGas);
  double *xval = (double *) mymalloc("x", sizeof(double) * NumGas);
  double *bval = (double *) mymalloc("b", sizeof(double) * NumGas);

  for(i = 0; i < NumGas; i++)
    {
      rows[i] = fld_offsets[ThisTask] + i;
      xval[i] = SphP[i].Utherm;
      bval[i] = U0[i];
    }

  HYPRE_IJVectorSetValues(*b, NumGas, rows, bval);
  HYPRE_IJVectorSetValues(*x, NumGas, rows, xval);

  myfree(bval);
  myfree(xval);
  myfree(rows);


  double u_K;

  int dp_kid, dp1_kid, dp2_kid, dp3_kid;
  double dp_akss, dp1_akss, dp2_akss, dp3_akss, dp_alss;
  double F2_K, F2_L;


  double globtotsum;

  int *idxs = (int *) mymalloc("idxs", sizeof(int) * MaxTetra);


  int *cols = mymalloc("cols", sizeof(int) * MaxTetra);
  double *vals = mymalloc("vals", sizeof(double) * MaxTetra);


  double totsum = 0.0;
  for(i = 0; i < NumGas; i++)
    {
      cols[0] = fld_offsets[ThisTask] + i;
      vals[0] = 0.0;

      int iii;
      for(iii = 0; iii < MaxTetra; iii++)
        {
          idxs[iii] = -2;
          vals[iii] = 0.0;
        }



      double sum = 0.0;

      double sum2 = 0.0;
      int ind = 1;
      int q = SphP[i].first_connection;

      u_K = SphP[i].Utherm;

      while(q >= 0)
        {
          int dp = DC[q].dp_index;
          int vf = DC[q].vf_index;
          int particle = DP[dp].index;

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
                  dp_kid = ConductionFaceData[vf].dp_kidp1;
                  dp1_kid = ConductionFaceData[vf].dp1_kidp1;
                  dp2_kid = ConductionFaceData[vf].dp2_kidp1;
                  dp3_kid = ConductionFaceData[vf].dp3_kidp1;

                  dp_akss = ConductionFaceData[vf].dp_akssp1;
                  dp1_akss = ConductionFaceData[vf].dp1_akssp1;
                  dp2_akss = ConductionFaceData[vf].dp2_akssp1;
                  dp3_akss = ConductionFaceData[vf].dp3_akssp1;

                  dp_alss = ConductionFaceData[vf].dp_akssp1p2;

                  F2_K = ConductionFaceData[vf].F2p1;
                  F2_L = ConductionFaceData[vf].F2p1p2;
                }
              else
                {
                  dp_kid = ConductionFaceData[vf].dp_kidp2;
                  dp1_kid = ConductionFaceData[vf].dp1_kidp2;
                  dp2_kid = ConductionFaceData[vf].dp2_kidp2;
                  dp3_kid = ConductionFaceData[vf].dp3_kidp2;

                  dp_akss = ConductionFaceData[vf].dp_akssp2;
                  dp1_akss = ConductionFaceData[vf].dp1_akssp2;
                  dp2_akss = ConductionFaceData[vf].dp2_akssp2;
                  dp3_akss = ConductionFaceData[vf].dp3_akssp2;

                  dp_alss = ConductionFaceData[vf].dp_akssp2p1;

                  F2_K = ConductionFaceData[vf].F2p2;
                  F2_L = ConductionFaceData[vf].F2p2p1;
                }

              double u_L, u_L1, u_L2, u_L3, mu_K, mu_L, mup_K, alpha, beta;

              u_L = get_Utherm(dp);

              u_L1 = get_Utherm(dp1_kid);

              u_L2 = get_Utherm(dp2_kid);

              u_L3 = get_Utherm(dp3_kid);

              mu_K = (mod(F2_L) + epsilon) / (mod(F2_K) + mod(F2_L) + 2.0 * epsilon);
              mu_L = 1.0 - mu_K;

              mup_K = mod(F2_L) / (mod(F2_K) + mod(F2_L) + 2.0 * epsilon);


              alpha = mu_K * dp_akss + mu_L * dp_alss;

              //beta = mup_K*(1.0 - erf(F2_K*F2_L/4.0)) ; 
              beta = mup_K * (1.0 - sign(F2_K * F2_L));
              int count;
              double dval;
              int idx;



              /*
                 if(i == 0 || i == 1)
                 printf("i = %d\t particle = %d\t dp_akss = %le\t dp_alss  = %le\t alpha = %le\n", i, particle, dp_akss, dp_alss, alpha) ;
               */
              for(count = 0; count < ind; count++)
                if(dp == idxs[count])
                  break;
              idx = count;
              idxs[count] = dp;
              cols[idx] = fld_offsets[DP[dp].task] + calculate_cols(dp);
              dval = -dt * alpha;
              vals[idx] += dval;
              sum += dval;
              sum2 += dval * u_L;
              if(idx == ind)
                ind++;

              for(count = 0; count < ind; count++)
                if(dp1_kid == idxs[count])
                  break;
              idx = count;
              idxs[count] = dp1_kid;
              cols[idx] = fld_offsets[DP[dp1_kid].task] + calculate_cols(dp1_kid);
              dval = -dt * beta * dp1_akss;
              vals[idx] += dval;
              sum += dval;
              sum2 += dval * u_L1;
              if(idx == ind)
                ind++;




              for(count = 0; count < ind; count++)
                if(dp2_kid == idxs[count])
                  break;
              idxs[count] = dp2_kid;
              idx = count;
              cols[idx] = fld_offsets[DP[dp2_kid].task] + calculate_cols(dp2_kid);
              dval = -dt * beta * dp2_akss;
              vals[idx] += dval;
              sum += dval;
              sum2 += dval * u_L2;
              if(idx == ind)
                ind++;


              if(dp3_kid != -1)
                {
                  for(count = 0; count < ind; count++)
                    if(dp3_kid == idxs[count])
                      break;
                  idxs[count] = dp3_kid;
                  idx = count;
                  cols[idx] = fld_offsets[DP[dp3_kid].task] + calculate_cols(dp3_kid);
                  dval = -dt * beta * dp3_akss;
                  vals[idx] += dval;
                  sum += dval;
                  sum2 += dval * u_L3;
                  if(idx == ind)
                    ind++;
                }
            }
          if(q == SphP[i].last_connection)
            break;

          q = DC[q].next;
        }

      vals[0] = 1 - sum;
      sum2 += vals[0] * u_K - U0[i];
      totsum += square(sum2);


      int item = fld_offsets[ThisTask] + i;

      HYPRE_IJMatrixSetValues(*A, 1, &ind, &item, cols, vals);
    }
  MPI_Allreduce(&totsum, &globtotsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  myfree(vals);
  myfree(cols);
  myfree(idxs);

  return sqrt(globtotsum);
}



int calculate_cols(int dp)
{
  int idx;
  int index = Mesh.DP[dp].index;
  if(Mesh.DP[dp].task == ThisTask)
    {
      if(index >= NumGas)
        index -= NumGas;

      idx = index;
    }
  else
    idx = Mesh.DP[dp].originalindex;

  return idx;
}


double linear_iterations(double dt)
{

  double error;

  HYPRE_IJMatrix A;
  HYPRE_ParCSRMatrix parcsr_A;
  HYPRE_IJVector b;
  HYPRE_ParVector par_b;
  HYPRE_IJVector x;
  HYPRE_ParVector par_x;

  HYPRE_Solver solver;
  HYPRE_Solver precond;

  int *fld_offsets = (int *) mymalloc("fld_offsets", sizeof(int) * (NTask + 1));

  MPI_Allgather(&NumGas, 1, MPI_INT, fld_offsets, 1, MPI_INT, MPI_COMM_WORLD);

  int i;
  int cumm = 0;
  for(i = 0; i <= NTask; i++)
    {
      int tmp = cumm;
      cumm += fld_offsets[i];
      fld_offsets[i] = tmp;
    }

  HYPRE_IJMatrixCreate(MPI_COMM_WORLD, fld_offsets[ThisTask], fld_offsets[ThisTask + 1] - 1, fld_offsets[ThisTask], fld_offsets[ThisTask + 1] - 1, &A);
  HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
  HYPRE_IJMatrixInitialize(A);

  HYPRE_IJVectorCreate(MPI_COMM_WORLD, fld_offsets[ThisTask], fld_offsets[ThisTask + 1] - 1, &b);
  HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(b);

  HYPRE_IJVectorCreate(MPI_COMM_WORLD, fld_offsets[ThisTask], fld_offsets[ThisTask + 1] - 1, &x);
  HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(x);


  error = set_coeff(&A, &x, &b, fld_offsets, dt);

  MPI_Barrier(MPI_COMM_WORLD);

  HYPRE_IJMatrixAssemble(A);
  HYPRE_IJMatrixGetObject(A, (void **) &parcsr_A);

  HYPRE_IJVectorAssemble(b);
  HYPRE_IJVectorGetObject(b, (void **) &par_b);

  HYPRE_IJVectorAssemble(x);
  HYPRE_IJVectorGetObject(x, (void **) &par_x);

  int num_iterations;
  double final_res_norm;

  /*
     char fname[200] ;
     sprintf(fname, "MatrixA.dat") ;
     HYPRE_IJMatrixPrint(A, fname) ;
     terminate("Matrix printed") ;
   */

  /* Create solver */
  HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, &solver);

  /*Set Parameters */
  HYPRE_ParCSRGMRESSetTol(solver, epsilon_lin);
  HYPRE_ParCSRGMRESSetMaxIter(solver, MAX_ITER);
  HYPRE_ParCSRGMRESSetPrintLevel(solver, 3);

  HYPRE_ParCSRGMRESSetTol(solver, epsilon_lin); /* conv. tolerance */
  HYPRE_ParCSRGMRESSetMaxIter(solver, MAX_ITER);

  /* use BoomerAMG as preconditioner */
  HYPRE_BoomerAMGCreate(&precond);
  HYPRE_BoomerAMGSetCoarsenType(precond, 10);
  HYPRE_BoomerAMGSetStrongThreshold(precond, 0.25);
  HYPRE_BoomerAMGSetTol(precond, 0.0);
  HYPRE_BoomerAMGSetPrintLevel(precond, 1);
  HYPRE_BoomerAMGSetMaxIter(precond, 2);

  /* set the preconditioner */
  MPI_Barrier(MPI_COMM_WORLD);
  HYPRE_ParCSRGMRESSetPrecond(solver, HYPRE_BoomerAMGSolve, HYPRE_BoomerAMGSetup, precond);



  HYPRE_ParCSRGMRESSetup(solver, parcsr_A, par_b, par_x);


  MPI_Barrier(MPI_COMM_WORLD);




  HYPRE_ParCSRGMRESSolve(solver, parcsr_A, par_b, par_x);
  HYPRE_ParCSRGMRESGetNumIterations(solver, &num_iterations);
  HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);

  /* Destroy solver */
  HYPRE_ParCSRGMRESDestroy(solver);
  HYPRE_BoomerAMGDestroy(precond);
  get_xval(&x, fld_offsets);

  HYPRE_IJMatrixDestroy(A);
  HYPRE_IJVectorDestroy(b);
  HYPRE_IJVectorDestroy(x);

  myfree(fld_offsets);


  if(final_res_norm > epsilon_lin)
    mpi_printf("HYPRE failed to converge\n");
  //terminate("HYPRE failed to converge") ;

  mpi_printf("HYPRE: iter %d, res_norm %g\n", num_iterations, final_res_norm);
  return error;
}


void get_xval(HYPRE_IJVector * x, int *fld_offsets)
{

  int *rows = (int *) mymalloc("rows", sizeof(int) * NumGas);
  double *XVec = (double *) mymalloc("XVec", sizeof(double) * NumGas);
  int i;

  for(i = 0; i < NumGas; i++)
    rows[i] = fld_offsets[ThisTask] + i;

  HYPRE_IJVectorGetValues(*x, NumGas, rows, XVec);

#ifndef UMC_CORRECTOR
  double rel_diff_sum = 0.0;
  double glob_sum;
  double norm = 0.0;
  double glob_norm;

  for(i = 0; i < NumGas; i++)
    {
      rel_diff_sum += square(SphP[i].Utherm - (PICARD_BETA * XVec[i] + (1.0 - PICARD_BETA) * SphP[i].Utherm));
      norm += square(SphP[i].Utherm);
      SphP[i].Utherm = PICARD_BETA * XVec[i] + (1.0 - PICARD_BETA) * SphP[i].Utherm;
    }

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Allreduce(&rel_diff_sum, &glob_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  MPI_Allreduce(&norm, &glob_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


  u_error = sqrt(glob_sum) / sqrt(glob_norm);

#else

  double c1, c1_norm, c2, c2_norm, c1c2dot;

  double sumdot;
  double sumdot_global;
  double c1_global;
  double c2_global;
  double rel_diff = 0.0;
  double rel_diff_global;
  double norm = 0.0;
  double glob_norm;

  sumdot = 0.0;
  c1_norm = 0.0;
  c2_norm = 0.0;


  for(i = 0; i < NumGas; i++)
    {
      c1 = XVec[i] - SphP[i].Utherm;
      c2 = SphP[i].Utherm - Uminus1[i];

      c1c2dot = c1 * c2;
      sumdot += c1c2dot;
      c1_norm += square(c1);
      c2_norm += square(c2);
    }


  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Allreduce(&sumdot, &sumdot_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&c1_norm, &c1_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&c2_norm, &c2_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if(c1_global == 0.0)
    c1_global = 1.0;

  if(c2_global == 0.0)
    c2_global = 1.0;

  double theta_gamma;

  theta_gamma = acos(sumdot_global / sqrt(c1_global * c2_global));

  double mu_gamma;


  if(theta_gamma <= M_PI / 8.0)
    {
      if(mgm1 == 0.3)
        mu_gamma = 0.5;
      else if(mgm1 == 0.5)
        mu_gamma = 1.0;
      else
        mu_gamma = 1.2;
    }
  else if((theta_gamma > M_PI / 8.0) && (theta_gamma <= 5.0 * M_PI / 20.0))
    {
      if(mgm1 == 0.3)
        mu_gamma = 0.5;
      else
        mu_gamma = 1.0;
    }
  else if((theta_gamma > 10.0 * M_PI / 20.0) && (theta_gamma <= 18.0 * M_PI / 20.0))
    mu_gamma = 0.5;
  else
    mu_gamma = 0.3;



  mpi_printf("CONDUCTION: theta_gamma = %le\n", theta_gamma);
  mpi_printf("CONDUCTION: mu_gamma = %le\n", mu_gamma);
  mpi_printf("CONDUCTION: mgm1 = %le\n", mgm1);

  mgm1 = mu_gamma;

  for(i = 0; i < NumGas; i++)
    {
      rel_diff += square(XVec[i] - SphP[i].Utherm);
      norm += square(SphP[i].Utherm);
      Uminus1[i] = SphP[i].Utherm;
      SphP[i].Utherm = SphP[i].Utherm + mu_gamma * (XVec[i] - SphP[i].Utherm);
    }

  MPI_Allreduce(&rel_diff, &rel_diff_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&norm, &glob_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  u_error = sqrt(rel_diff_global) / sqrt(glob_norm);

#endif

  myfree(XVec);
  myfree(rows);

}




void U_exchange_vector()
{
  int listp;
  int j, p, task, off;
  int ngrp, recvTask, place;    /* there was earlier the following defined: sendTask */


  struct uexch *tmpUExch;
  tmpUExch = (struct uexch *) mymalloc("tmpUExch", Mesh_nexport * sizeof(struct uexch));

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

                  tmpUExch[off].Utherm = SphP[place].Utherm;
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
              MPI_Sendrecv(&tmpUExch[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask] *
                           sizeof(struct uexch), MPI_BYTE, recvTask, TAG_DENS_A, &UExch[Mesh_Recv_offset[recvTask]],
                           Mesh_Recv_count[recvTask] * sizeof(struct uexch), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  myfree(tmpUExch);
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
#ifdef CONDUCTION_ANISOTROPIC
                  tmpInExch[off].B[0] = SphP[place].B[0];
                  tmpInExch[off].B[1] = SphP[place].B[1];
#ifndef TWODIMS
                  tmpInExch[off].B[2] = SphP[place].B[2];
#endif
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
                      tmpFluxExch[off].alss[init] = 0.0;
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
                                tmpFluxExch[off].alss[Noc] = ConductionFaceData[vf].dp_akssp1;
                              else
                                tmpFluxExch[off].alss[Noc] = ConductionFaceData[vf].dp_akssp2;
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




void F2_exchange_vector()
{
  int listp;
  int j, p, task, off;
  int ngrp, recvTask, place;    /* there was earlier the following defined: sendTask */

  point *DP = Mesh.DP;
  face *VF = Mesh.VF;


  struct f2exch *tmpF2Exch;
  tmpF2Exch = (struct f2exch *) mymalloc("tmpF2Exch", Mesh_nexport * sizeof(struct f2exch));

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
                    tmpF2Exch[off].F2[init] = 0.0;

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
                              if(nindex == P[place].ID)
                                tmpF2Exch[off].F2[Noc] = ConductionFaceData[vf].F2p1;
                              else
                                tmpF2Exch[off].F2[Noc] = ConductionFaceData[vf].F2p2;
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
              MPI_Sendrecv(&tmpF2Exch[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask] *
                           sizeof(struct f2exch), MPI_BYTE, recvTask, TAG_DENS_A, &F2Exch[Mesh_Recv_offset[recvTask]],
                           Mesh_Recv_count[recvTask] * sizeof(struct f2exch), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  myfree(tmpF2Exch);
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

double compute_timestep(double deltat, int i, int N)
{
  double dN = ((double) (N));
  double di = ((double) (i + 1));

  double sts = deltat / dN / dN * di * 2.0 * dN / (dN - 1.0);
  return sts;
}

#endif /*IMPLICIT_TI */
#endif /*MONOTONE_CONDUCTION */
