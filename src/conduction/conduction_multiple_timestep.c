/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/conduction/conduction_multiple_timestep.c
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


/*! \file conduction_multiple_timestep.c
 *  \brief main driver for an anisotropic/isotropic diffusion solver
 *
 *  This file contains the code for a extremum preserving anisotrpic diffusion solver. This implementation is based on the method outlined in Gao & Wu (2013) - Journal of Computational Physics 10/2013; 250:308-331. DOI: 10.1016/j.jcp.2013.05.013 .  

 Semi implicit time integration - Split the Non-linear system ( M(U)U = F(U) ) into explicit and implicit parts as in Sharma & Hammett 2011, Journal of Computational Physics 10.1016/j.jcp.2011.03.009 solved using  GMRES and multigrid preconditioner (using HYPRE library).  
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
#ifdef SEMI_IMPLICIT_TI
#ifdef MULTIPLE_TIME_STEPPING

#ifndef DOUBLE_STENCIL
terminate("MONOTONE CONDUCTION: MULTIPLE_TIME_STEPPING requires DOUBLE_STENCIL option") ;
#endif



#ifdef TWODIMS
#define MaxNoc 100
#else
#define MaxNoc 500
#endif

#define MaxFluxExch 100


#define MaxTetra MaxNoc
#define epsilon 1e-25
#define MINANGLE 1e-5
#define gamma 0.0

#ifdef TWODIMS
#define epsilon_lin 1e-13
#else
#define epsilon_lin 1e-8
#endif
#define MAX_ITER 100

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
  double SurfaceArea;
#ifdef CONDUCTION_ANISOTROPIC
  double B[DIMS];
#endif
}
 *InExch;

static struct fluxexch
{
  double F2[MaxFluxExch];
  double alss[MaxFluxExch];
  MyIDType ID[MaxFluxExch];
}
 *FluxExch;

#ifndef TWODIMS
static struct savedtetradata
{
  double y_sigma_p1[3];
  double y_sigma_p2[3];
  double y_sigma_p3[3];
  double mdls_ysigma_p1, mdls_ysigma_p2, mdls_ysigma_p3;
  double omega1, omega2, omega3;
  int p1Noc, p2Noc, p3Noc;
  int p1task, p2task, p3task;
  double u1, u2, u3;
  int originalindex1, originalindex2, originalindex3;
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
  //int p1task, p2task ;
  double u1, u2;
}
 *datatetra;
#endif

static MyFloat *Kappa;
static MyFloat *U0;
//static MyFloat *Uminus1 ;
//static MyFloat u_error ;
//static double mgm1 ;


static MyFloat *Uinit;
static MyFloat *U_con;


void monotone_conduction()
{
  mpi_printf("CONDUCTION: Doing Conduction For This Time Step\n");
  TIMER_START(CPU_CONDUCTION) ;
  double dt_init, dt, dt_sum;
  Kappa = (MyFloat *) mymalloc("Kappa", NumGas * sizeof(MyFloat));

  U0 = (MyFloat *) mymalloc("U0", NumGas * sizeof(MyFloat));
  U_con = (MyFloat *) mymalloc("Uminus1", NumGas * sizeof(MyFloat));


  Uinit = (MyFloat *) mymalloc("Uinit", NumGas * sizeof(MyFloat));

  int numcycles, contsubcycle;

  dt_init = get_time_step(&numcycles);

  mpi_printf("total time step = %g\tnumcycles = %d\n", dt_init, numcycles) ;

  initialize_u_init();
  int i=0;

  contsubcycle = -1 ;

  dt_sum = 0.0 ;
  if(dt_init > 0)
    {
      while(contsubcycle<0)
	{
	  mpi_printf("CONDUCTION: Total loops = %d, no. of this loop = %d\n", numcycles, i + 1);

	  i++ ;

	  if(i>500)
	    contsubcycle = 1 ;

	  initialize_kappa();
	  
	  dt = get_subcycle_timestep() ;
	  if(dt_sum+dt >= dt_init)
	    {
	      dt = dt_init - dt_sum ;
	      contsubcycle = 1 ;
	    }
	  dt_sum += dt ;
	  mpi_printf("CONDUCTION: Total time Step = %le \t This time step  = %le \t Sum dt = %le \n", dt_init, dt, dt_sum) ;	  
#ifdef TWODIMS
	  calculate_unidirectional_fluxes_2D();
#else
	  calculate_unidirectional_fluxes_3D();
#endif
	  prepare_linear_iterations(dt);
	}
    }
  update_U_and_Energy();
  myfree(Uinit);
  myfree(U_con);
  myfree(U0);
  myfree(Kappa);

  TIMER_STOP(CPU_CONDUCTION) ;

  mpi_printf("CONDUCTION: DONE\n");
}



void initialize_u_init()
{
  int i;
  for(i = 0; i < NumGas; i++)
    U_con[i] = SphP[i].Utherm;
}


void update_U_and_Energy()
{
  int i;
  double du;
  for(i = 0; i < NumGas; i++)
    {
      du = SphP[i].Utherm - U_con[i];   /* [Temp] = erg/g   */
      SphP[i].Utherm = U_con[i];
      SphP[i].Energy += All.cf_atime * All.cf_atime * du * P[i].Mass;
      set_pressure_of_cell(i);
    }
}


double get_time_step(int *numcycles)
{
  double dt;
  dt = (All.conduction_Ti_endstep - All.conduction_Ti_begstep) * All.Timebase_interval;
  dt *= All.cf_atime / All.cf_time_hubble_a;

  double temp = dt / All.dt_conduction;
  if(temp == round(temp))
    *numcycles = ((int) (round(temp)));
  else
    *numcycles = ((int) (round(dt / All.dt_conduction))) + 1;
  mpi_printf("dt = %le\tdt_conduction = %le\n", dt, All.dt_conduction) ;
  return dt;
}


double get_subcycle_timestep()
{
  double dt, dt_all, dt_min;
  int i ;
  double chi ;
  dt_min = MAX_DOUBLE_NUMBER ;
  for(i = 0 ; i < NumGas; i++)
    {
      chi = Kappa[i]*SphP[i].Volume ;
      double rad = get_cell_radius(i) ;

      if(Kappa[i] == 0.0)
	dt = MAX_DOUBLE_NUMBER ;
      else
	dt = 4.0 * rad * rad / chi ;


      if(dt<dt_min)
	dt_min = dt ;
    }
  MPI_Allreduce(&dt_min, &dt_all, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD) ;
  return dt_all ;
}

void initialize_kappa()
{
  int i, idx;

#ifdef RESTRICT_KAPPA
  double Kappa_max = All.Chi_max_int/All.cf_atime/All.cf_atime ;
#endif

  //for(i = 0; i < NumGas; i++)
  //{
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      Uinit[i] = SphP[i].Utherm;
      U0[i] = SphP[i].Utherm;

#ifdef USE_SFR
      if(SphP[i].Sfr>0.0)
	{
	  Kappa[i] = 0.0 ;
	  continue ;
	}
#endif
#ifdef CONDUCTION_CONSTANT
      Kappa[i] = All.ConductionCoeff / SphP[i].Volume;
#else
      Kappa[i] = All.ConductionCoeff * pow(SphP[i].Utherm, 2.5) / P[i].Mass;
#ifdef CONDUCTION_SATURATION
      double electron_free_path = All.ElectronFreePathFactor * SphP[i].Utherm * SphP[i].Utherm / (SphP[i].Density * All.cf_a3inv);
      double mod_dU = sqrt(SphP[i].Grad.dutherm[0] * SphP[i].Grad.dutherm[0] + SphP[i].Grad.dutherm[1] * SphP[i].Grad.dutherm[1] + SphP[i].Grad.dutherm[2] * SphP[i].Grad.dutherm[2]);
      if(mod_dU != 0.0)
        {
          double temp_scale_length = All.Time * fabs(SphP[i].Utherm) / mod_dU;

          Kappa[i] /= (1 + 4.2 * electron_free_path / temp_scale_length);
        }
#endif
#endif

      if(All.ComovingIntegrationOn)
        Kappa[i] *= All.cf_atime;

#ifdef RESTRICT_KAPPA
      if(Kappa[i]>Kappa_max/SphP[i].Volume)
	Kappa[i] = Kappa_max/SphP[i].Volume ; 
#endif
      //mpi_printf("CONDUCTION: All cf time = %le\t Kappa = %le\n", All.cf_atime, Kappa[i]) ;

      
    }

  for(i = 0; i < Mesh.Nvf; i++)
    {
      Mesh.VF[i].akssp1 = 0.0;
      Mesh.VF[i].akssp2 = 0.0;
      Mesh.VF[i].akssp1p2 = 0.0;
      Mesh.VF[i].akssp2p1 = 0.0;
      Mesh.VF[i].F2p1 = 0.0;
      Mesh.VF[i].F2p2 = 0.0;
      Mesh.VF[i].F2p1p2 = 0.0;
      Mesh.VF[i].F2p2p1 = 0.0;
    }
}



void init_conductivity(void)
{
  All.ConductionCoeff = 1.0 * All.ConductionEfficiency; /* Only for test and only if CONDUCTION_CONSTANT IS ON!!! */

  /*  double meanweight;

     meanweight = m_p * 4.0 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));    assuming full ionization */
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

  // printf("unit conv1 = %le\n", meanweight / k_B * GAMMA_MINUS1);
  All.ConductionCoeff *= meanweight / k_B * GAMMA_MINUS1;

  /* The conversion of ConductionCoeff between internal units and cgs
   * units involves one factor of 'h'. We take care of this here.
   */
  //All.ConductionCoeff /= All.HubbleParam;

  /* Include the fraction factor f, which is called ConductionEfficiency,
   * of the classical Spitzer conductivity - see (ApJ 582:162-169, Eq.(5))
   */
  All.ConductionCoeff *= All.ConductionEfficiency;

  mpi_printf("MONOTONE CONDUCTION: Conduction Coeff = %le\n", All.ConductionCoeff);


#ifdef RESTRICT_KAPPA
  All.Chi_max_int = All.MaxDiffusivity*cm*cm/s ;
  /* Take care of h factor*/
  //All.Chi_max_int *= All.HubbleParam ;
#endif


#ifdef CONDUCTION_SATURATION

    double meanweight_cgs = PROTONMASS * 4.0 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));
  All.ElectronFreePathFactor = pow(3.0, 1.5)*GAMMA_MINUS1*GAMMA_MINUS1*meanweight_cgs*meanweight_cgs*meanweight_cgs*(3+HYDROGEN_MASSFRAC)/(8.0*sqrt(M_PI)*pow(ELECTRONCHARGE, 4)* coulomb_log) ;

  All.ElectronFreePathFactor /= erg*erg*cm*cm/g/g/g ; 
 
  //      All.ElectronFreePathFactor = 8 * pow(3.0, 1.5) * pow(GAMMA_MINUS1, 2) / pow(3 + 5 * HYDROGEN_MASSFRAC, 2) / (1 + HYDROGEN_MASSFRAC) / sqrt(M_PI) / coulomb_log * pow(PROTONMASS, 3) / pow(ELECTRONCHARGE, 4)  / (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam) * pow(All.UnitPressure_in_cgs / All.UnitDensity_in_cgs, 2);

  /* If the above value is multiplied with u^2/rho in code units (with rho being the physical density), then
   * one gets the electrong mean free path in centimeter. Since we want to compare this with another length
   * scale in code units, we now add an additional factor to convert back to code units.
   */
  //All.ElectronFreePathFactor *= All.HubbleParam / All.UnitLength_in_cm;

  mpi_printf("MONOTONE CONDUCTION: Electron Free Path  Coeff = %le\n", All.ElectronFreePathFactor);

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


  const int edge_start[3] = { 1, 2, 0 };
  const int edge_end[3] = { 2, 0, 1 };
  const int edge_nexttetra[3] = { 0, 1, 2 };


  double Lambda_K[DIMS][DIMS];

  double Kappa_K, modB, B_K[DIMS], Kappa_L, B_L[DIMS], mdls;
  double omega_K, omega_L, area_kl[MaxNoc];


  int tetraindexes[MaxTetra];

  double mdls_ltnk[MaxNoc];

  double sfarea;

  double modBL;

  int particle_id[MaxNoc];
  int face_id[MaxNoc];


  double u_K, u_L, u;

  InExch = (struct inexch *) mymalloc("InExch", Mesh_nimport * sizeof(struct inexch));

  datatetra = (struct savedtetradata *) mymalloc("datatetra", MaxTetra * sizeof(struct savedtetradata));

  conduction_exchange_vector();

  int tottetra;

  int dp = -1;
  int vf = -1;
  int particle = -1;

  //for(i = 0; i < NumGas; i++)
  //{
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
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

      u_K = SphP[i].Utherm;

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
          particle_id[Noc] = dp;

          if(DP[dp].task == ThisTask)
            {
              if(particle >= NumGas)
                particle -= NumGas;

              Kappa_L = Kappa[particle];
              u_L = SphP[particle].Utherm;
              sfarea = SphP[particle].SurfaceArea;
	      if(TimeBinSynchronized[DP[dp].timebin] != 1)
		

#ifdef CONDUCTION_ANISOTROPIC
              B_L[0] = SphP[particle].B[0];
              B_L[1] = SphP[particle].B[1];
#endif
            }
          else
            {
              Kappa_L = InExch[particle].Kappa;
              u_L = InExch[particle].Utherm;
              sfarea = InExch[particle].SurfaceArea;
#ifdef CONDUCTION_ANISOTROPIC
              B_L[0] = InExch[particle].B[0];
              B_L[1] = InExch[particle].B[1];
#endif
            }


          if((VF[vf].area > 1e-10 * SphP[i].SurfaceArea) || (VF[vf].area > 1e-10 * sfarea))
            VF[vf].doface = 1;
          else
            VF[vf].doface = 0;

#ifdef CONDUCTION_ANISOTROPIC
          modBL = sqrt(B_L[0] * B_L[0] + B_L[1] * B_L[1]);

          if(modBL == 0.0)
            modBL = 1.0;

          B_L[0] /= modBL;
          B_L[1] /= modBL;
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
          area_kl[Noc] = VF[vf].area;


          lTnK[Noc][0] = Lambda_K[0][0] * n[0] + Lambda_K[0][1] * n[1];
          lTnK[Noc][1] = Lambda_K[1][0] * n[0] + Lambda_K[1][1] * n[1];


          mdls_ltnk[Noc] = sqrt(square(lTnK[Noc][0]) + square(lTnK[Noc][1]));

          if(mdls_ltnk[Noc] == 0.0)
            mdls_ltnk[Noc] = 1.0;

          lTnK[Noc][0] /= mdls_ltnk[Noc];
          lTnK[Noc][1] /= mdls_ltnk[Noc];

          u = omega_K * u_K + omega_L * u_L;

          double px1, py1;

          if(omega_L == 0.0)
            {
              omega_L = MINANGLE;
              omega_K = 1.0 - MINANGLE;
            }

          px1 = omega_L * (DP[dp].x - P[i].Pos[0]);
          py1 = omega_L * (DP[dp].y - P[i].Pos[1]);


          mdls = sqrt(square(px1) + square(py1));


          if(mdls == 0.0)
            mdls = 1.0;


          px1 /= mdls;
          py1 /= mdls;


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

          int idp1, idp2;
          int taskp1, taskp2;
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
                  datatetra[tottetra].u1 = u;

                  idp2 = DP[this->p[l2]].index;
                  taskp2 = DP[this->p[l2]].task;

                  datatetra[tottetra].p2Noc = this->p[l2];


                  if(taskp2 == ThisTask)
                    {

                      if(idp2 >= NumGas)
                        idp2 -= NumGas;

                      Kappa_L = Kappa[idp2];
                      u_L = SphP[idp2].Utherm;
#ifdef CONDUCTION_ANISOTROPIC
                      B_L[0] = SphP[idp2].B[0];
                      B_L[1] = SphP[idp2].B[1];
#endif
                    }
                  else
                    {
                      Kappa_L = InExch[idp2].Kappa;
                      u_L = InExch[idp2].Utherm;
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

                  if(omega_L == 0.0)
                    {
                      omega_K = 1.0 - MINANGLE;
                      omega_L = MINANGLE;
                    }

                  px = omega_L * (DP[this->p[l2]].x - P[i].Pos[0]);
                  py = omega_L * (DP[this->p[l2]].y - P[i].Pos[1]);

                  datatetra[tottetra].omega2 = omega_K;
                  datatetra[tottetra].mdls_ysigma_p2 = sqrt(px * px + py * py);
                  if(datatetra[tottetra].mdls_ysigma_p2 == 0.0)
                    datatetra[tottetra].mdls_ysigma_p2 = 1.0;


                  datatetra[tottetra].u2 = omega_K * u_K + omega_L * u_L;
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

          if(VF[face_id[i1]].doface == 1)
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


                  if((a > -MINANGLE) && (b > -MINANGLE))
                    break;
                }

              if(a >= -MINANGLE && a < 0.0)
                a = 0.0;

              if(b >= -MINANGLE && b < 0.0)
                b = 0.0;


              if((a < 0.0) || (b < 0.0))
                {
                  printf("i = %d\n", i);
                  printf("Noc = %d i1 = %d tottetra = %d j1 = %d \n", Noc, i1, tottetra, j1);
                  printf("ltnk x = %le\ty = %le\n", x, y);
                  printf("1 - x = %le \t y = %le \n", x1, y1);
                  printf("2 - x = %le \t y = %le \n", x2, y2);
                  printf("a = %le\tb=%le\n", a, b);
                  terminate("|||| Fatal  Error one of a,b < 0\n");
                }


              double aKL1, aKL2;
              double zycoeff1, zycoeff2;

              aKL1 = area_kl[i1] * mdls_ltnk[i1] * a / datatetra[j1].mdls_ysigma_p1;
              aKL2 = area_kl[i1] * mdls_ltnk[i1] * b / datatetra[j1].mdls_ysigma_p2;

              zycoeff1 = area_kl[i1] * mdls_ltnk[i1] * a * (u_K - datatetra[j1].u1) / datatetra[j1].mdls_ysigma_p1;
              zycoeff2 = area_kl[i1] * mdls_ltnk[i1] * b * (u_K - datatetra[j1].u2) / datatetra[j1].mdls_ysigma_p2;

              int vfidx = face_id[i1];
              int to = VF[vfidx].p1;
              MyIDType index = DP[to].ID;


              int kid;
              double F2, akss;


              if(particle_id[i1] == datatetra[j1].p1Noc)
                {
                  akss = aKL1 * (1.0 - datatetra[j1].omega1);
                  F2 = zycoeff2;
                  kid = datatetra[j1].p1Noc;
                }
              else if(particle_id[i1] == datatetra[j1].p2Noc)
                {
                  akss = aKL2 * (1.0 - datatetra[j1].omega2);
                  F2 = zycoeff1;
                  kid = datatetra[j1].p2Noc;
                }
              else
                {
                  akss = 0.0;
                  F2 = zycoeff1 + zycoeff2;
                  kid = -1;
                }


              if(index == P[i].ID)
                {
                  VF[vfidx].kidp1 = kid;
                  VF[vfidx].F2p1 = F2;
                  VF[vfidx].akssp1 = akss;
                }
              else
                {
                  VF[vfidx].kidp2 = kid;
                  VF[vfidx].F2p2 = F2;
                  VF[vfidx].akssp2 = akss;
                }
            }
        }
    }
  myfree(datatetra);
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

  double Kappa_K, modB, modBL, B_K[DIMS], Kappa_L, B_L[DIMS], mdls;
  double omega_K, omega_L, area_kl[MaxNoc];

  double mdls_ltnk[MaxNoc];

  int tetraindexes[MaxTetra];

  double lambda_K, lambda_L;

  int face_id[MaxNoc];

  double px1, py1, pz1;

  double sfarea;

  int particle_id[MaxNoc];

  double u_K, u_L, u;

  InExch = (struct inexch *) mymalloc("InExch", Mesh_nimport * sizeof(struct inexch));

  datatetra = (struct savedtetradata *) mymalloc("datatetra", MaxTetra * sizeof(struct savedtetradata));

  conduction_exchange_vector();

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
          particle_id[Noc] = dp;

          if(DP[dp].task == ThisTask)
            {

              if(particle >= NumGas)
                particle -= NumGas;

              Kappa_L = Kappa[particle];
              u_L = SphP[particle].Utherm;
              sfarea = SphP[particle].SurfaceArea;

#ifdef CONDUCTION_ANISOTROPIC
              B_L[0] = SphP[particle].B[0];
              B_L[1] = SphP[particle].B[1];
              B_L[2] = SphP[particle].B[2];
#endif
            }
          else
            {
              Kappa_L = InExch[particle].Kappa;
              u_L = InExch[particle].Utherm;
              sfarea = InExch[particle].SurfaceArea;

#ifdef CONDUCTION_ANISOTROPIC
              B_L[0] = InExch[particle].B[0];
              B_L[1] = InExch[particle].B[1];
              B_L[2] = InExch[particle].B[2];
#endif
            }


	  if((VF[vf].area > 1e-10 * SphP[i].SurfaceArea) || (VF[vf].area > 1e-10 * sfarea))
            VF[vf].doface = 1;
          else
            VF[vf].doface = 0;

#ifdef USE_SFR
	  if(Kappa_K == 0.0 || Kappa_L == 0.0)
	    VF[vf].doface = 0 ;
#endif
      


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


          if(omega_L == 0.0)
            {
              omega_L = MINANGLE;
              omega_K = 1.0 - MINANGLE;
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
                  datatetra[tottetra].u1 = u;

                  idp2 = DP[this->p[l2]].index;
                  taskp2 = DP[this->p[l2]].task;
                  datatetra[tottetra].p2Noc = this->p[l2];

                  if(taskp2 == ThisTask)
                    {
                      if(idp2 >= NumGas)
                        idp2 -= NumGas;

                      Kappa_L = Kappa[idp2];
                      u_L = SphP[idp2].Utherm;

#ifdef CONDUCTION_ANISOTROPIC
                      B_L[0] = SphP[idp2].B[0];
                      B_L[1] = SphP[idp2].B[1];
                      B_L[2] = SphP[idp2].B[2];
#endif
                    }
                  else
                    {
                      Kappa_L = InExch[idp2].Kappa;
                      u_L = InExch[idp2].Utherm;

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


                  if(omega_L == 0.0)
                    {
                      omega_L = MINANGLE;
                      omega_K = 1.0 - MINANGLE;
                    }



                  datatetra[tottetra].u2 = omega_K * u_K + omega_L * u_L;
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
                      u_L = SphP[idp3].Utherm;

#ifdef CONDUCTION_ANISOTROPIC
                      B_L[0] = SphP[idp3].B[0];
                      B_L[1] = SphP[idp3].B[1];
                      B_L[2] = SphP[idp3].B[2];
#endif
                    }
                  else
                    {
                      Kappa_L = InExch[idp3].Kappa;
                      u_L = InExch[idp3].Utherm;

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


                  if(omega_L == 0.0)
                    {
                      omega_L = MINANGLE;
                      omega_K = 1.0 - MINANGLE;
                    }

                  datatetra[tottetra].u3 = omega_K * u_K + omega_L * u_L;

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

          if(VF[face_id[i1]].doface == 1)
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

              if(a >= -MINANGLE && a < MINANGLE)
                a = 0.0;

              if(b >= -MINANGLE && b < MINANGLE)
                b = 0.0;

              if(c >= -MINANGLE && c < MINANGLE)
                c = 0.0;

              double aKL1, aKL2, aKL3;

              double zycoeff1, zycoeff2, zycoeff3;

              aKL1 = area_kl[i1] * mdls_ltnk[i1] * a / datatetra[j1].mdls_ysigma_p1;
              aKL2 = area_kl[i1] * mdls_ltnk[i1] * b / datatetra[j1].mdls_ysigma_p2;
              aKL3 = area_kl[i1] * mdls_ltnk[i1] * c / datatetra[j1].mdls_ysigma_p3;

              zycoeff1 = area_kl[i1] * mdls_ltnk[i1] * a * (u_K - datatetra[j1].u1) / datatetra[j1].mdls_ysigma_p1;
              zycoeff2 = area_kl[i1] * mdls_ltnk[i1] * b * (u_K - datatetra[j1].u2) / datatetra[j1].mdls_ysigma_p2;
              zycoeff3 = area_kl[i1] * mdls_ltnk[i1] * c * (u_K - datatetra[j1].u3) / datatetra[j1].mdls_ysigma_p3;

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


              int vfidx = face_id[i1];
              int to = VF[vfidx].p1;
              MyIDType index = DP[to].ID;

              int kid;
              double F2, akss;

              if(particle_id[i1] == datatetra[j1].p1Noc)
                {
                  akss = aKL1 * (1.0 - datatetra[j1].omega1);
                  F2 = zycoeff2 + zycoeff3;
                  kid = datatetra[j1].p1Noc;
                }
              else if(particle_id[i1] == datatetra[j1].p2Noc)
                {
                  akss = aKL2 * (1.0 - datatetra[j1].omega2);
                  F2 = zycoeff1 + zycoeff3;
                  kid = datatetra[j1].p2Noc;
                }
              else if(particle_id[i1] == datatetra[j1].p3Noc)
                {
                  akss = aKL3 * (1.0 - datatetra[j1].omega3);
                  F2 = zycoeff1 + zycoeff2;
                  kid = datatetra[j1].p3Noc;
                }
              else
                {
                  akss = 0.0;
                  F2 = zycoeff1 + zycoeff2 + zycoeff3;
                  kid = -1;
                }



              if(index == P[i].ID)
                {
                  VF[vfidx].kidp1 = kid;
                  VF[vfidx].F2p1 = F2;
                  VF[vfidx].akssp1 = akss;
                }
              else
                {
                  VF[vfidx].kidp2 = kid;
                  VF[vfidx].F2p2 = F2;
                  VF[vfidx].akssp2 = akss;
                }
            }
        }
    }

  myfree(datatetra);
  myfree(InExch);
}

#endif

void prepare_linear_iterations(double dt)
{
  set_initial_coeff(dt);
  linear_iterations(dt);
}




void set_initial_coeff(double dt)
{

  FluxExch = (struct fluxexch *) mymalloc("FluxExch", Mesh_nimport * sizeof(struct fluxexch));
  flux_exchange_vector();


  int i;
  point *DP = Mesh.DP;
  face *VF = Mesh.VF;

  double u_K, u_L;

  double F2_L, F2_K;
  double akss, alss;

  for(i = 0; i < NumGas; i++)
    {
      double sum = 0.0;

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

          if(VF[vf].doface == 1)
            {
              int to = VF[vf].p1;
              MyIDType index = DP[to].ID;

              if(index == P[i].ID)
                {
                  F2_K = Mesh.VF[vf].F2p1;
                  akss = Mesh.VF[vf].akssp1;
                }
              else
                {
                  F2_K = Mesh.VF[vf].F2p2;
                  akss = Mesh.VF[vf].akssp2;
                }

              if(Mesh.DP[dp].task == ThisTask)
                {
                  if(particle < NumGas)
                    {
                      if(index == P[i].ID)
                        {
                          F2_L = Mesh.VF[vf].F2p2;
                          alss = Mesh.VF[vf].akssp2;
                        }
                      else
                        {
                          F2_L = Mesh.VF[vf].F2p1;
                          alss = Mesh.VF[vf].akssp1;
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
                          MyIDType index_L = DP[to_L].ID;

                          if(Mesh.DP[dp_L].task == ThisTask)
                            {
                              if(particle_L >= NumGas)
                                particle_L -= NumGas;

                              if(particle_L == i)
                                {
                                  if(index_L == P[particle].ID)
                                    {
                                      F2_L = Mesh.VF[vf_L].F2p1;
                                      alss = Mesh.VF[vf_L].akssp1;
                                    }
                                  else
                                    {
                                      F2_L = Mesh.VF[vf_L].F2p2;
                                      alss = Mesh.VF[vf_L].akssp2;
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
                  for(number = 0; number < MaxFluxExch; number++)
                    {
                      if(FluxExch[particle].ID[number] == P[i].ID)
                        break;
                    }
                  F2_L = FluxExch[particle].F2[number];
                  alss = FluxExch[particle].alss[number];
                }

              double mup_K, beta;


              mup_K = mod(F2_L) / (mod(F2_K) + mod(F2_L) + 2.0 * epsilon);

              if(index == P[i].ID)
                {
                  Mesh.VF[vf].F2p1p2 = F2_L;
                  Mesh.VF[vf].akssp1p2 = alss;
                }
              else
                {
                  Mesh.VF[vf].F2p2p1 = F2_L;
                  Mesh.VF[vf].akssp2p1 = alss;
                }

              beta = mup_K * (1.0 - sign(F2_K * F2_L));

              sum += -dt * beta * F2_K;

            }
          if(q == SphP[i].last_connection)
            break;

          q = DC[q].next;
        }
      U0[i] += sum;
    }
  myfree(FluxExch);
}


void set_coeff(HYPRE_IJMatrix * A, HYPRE_IJVector * x, HYPRE_IJVector * b, int *fld_offsets, double dt)
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


  int Kid;
  double akss, alss, F2_L, F2_K;


  int *cols = mymalloc("cols", sizeof(int) * MaxTetra);
  double *vals = mymalloc("vals", sizeof(double) * MaxTetra);

  for(i = 0; i < NumGas; i++)
    {
      cols[0] = fld_offsets[ThisTask] + i;
      vals[0] = 0.0;

      double sum = 0.0;

      int ind = 1;
      int q = SphP[i].first_connection;


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

          if(VF[vf].doface == 1)
            {
              int to = VF[vf].p1;
              MyIDType index = DP[to].ID;

              if(index == P[i].ID)
                {
                  Kid = Mesh.VF[vf].kidp1;
                  akss = Mesh.VF[vf].akssp1;
                  F2_K = Mesh.VF[vf].F2p1;
                  alss = Mesh.VF[vf].akssp1p2;
                  F2_L = Mesh.VF[vf].F2p1p2;
                }
              else
                {
                  Kid = Mesh.VF[vf].kidp2;
                  akss = Mesh.VF[vf].akssp2;
                  F2_K = Mesh.VF[vf].F2p2;
                  alss = Mesh.VF[vf].akssp2p1;
                  F2_L = Mesh.VF[vf].F2p2p1;
                }

              double mu_K, mu_L;
              double alpha;

              mu_K = (mod(F2_L) + epsilon) / (mod(F2_K) + mod(F2_L) + 2.0 * epsilon);
              mu_L = 1.0 - mu_K;


              alpha = mu_K * akss + mu_L * alss;

              double dval;
              int originalindex;
              if(DP[dp].task == ThisTask)
                {
                  if(particle < NumGas)
                    originalindex = particle;
                  else
                    originalindex = particle - NumGas;
                }
              else
                originalindex = DP[dp].originalindex;


              cols[ind] = fld_offsets[DP[dp].task] + originalindex;
              dval = -dt * alpha;
              sum += dval;
              vals[ind] = dval;
              ind++;
            }
          if(q == SphP[i].last_connection)
            break;

          q = DC[q].next;
        }

      vals[0] = 1.0 - sum;



      int item = cols[0];
      HYPRE_IJMatrixSetValues(*A, 1, &ind, &item, cols, vals);
    }
  myfree(vals);
  myfree(cols);
}




void linear_iterations(double dt)
{

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


  set_coeff(&A, &x, &b, fld_offsets, dt);

  //MPI_Barrier(MPI_COMM_WORLD) ;

  HYPRE_IJMatrixAssemble(A);
  HYPRE_IJMatrixGetObject(A, (void **) &parcsr_A);

  HYPRE_IJVectorAssemble(b);
  HYPRE_IJVectorGetObject(b, (void **) &par_b);

  HYPRE_IJVectorAssemble(x);
  HYPRE_IJVectorGetObject(x, (void **) &par_x);

  int num_iterations;
  double final_res_norm;


  /*char fname[200] ;
     sprintf(fname, "MatrixA.dat") ;
     HYPRE_IJMatrixPrint(A, fname) ;
     MPI_Barrier(MPI_COMM_WORLD) ;
     terminate("Printed Matrix") ;

     char fname1[200] ;
     sprintf(fname1, "Vectorx.dat") ;
     HYPRE_IJVectorPrint(x, fname1) ;
     MPI_Barrier(MPI_COMM_WORLD) ;


     char fname2[200] ;
     sprintf(fname2, "Vectorb.dat") ;
     HYPRE_IJVectorPrint(b, fname2) ;
     MPI_Barrier(MPI_COMM_WORLD) ;


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
  //MPI_Barrier(MPI_COMM_WORLD);
  HYPRE_ParCSRGMRESSetPrecond(solver, HYPRE_BoomerAMGSolve, HYPRE_BoomerAMGSetup, precond);



  HYPRE_ParCSRGMRESSetup(solver, parcsr_A, par_b, par_x);


  //MPI_Barrier(MPI_COMM_WORLD);




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
    terminate("HYPRE failed to converge");

  mpi_printf("HYPRE: iter %d, res_norm %g\n", num_iterations, final_res_norm);
}


void get_xval(HYPRE_IJVector * x, int *fld_offsets)
{

  int *rows = (int *) mymalloc("rows", sizeof(int) * NumGas);
  double *XVec = (double *) mymalloc("XVec", sizeof(double) * NumGas);
  int i;

  for(i = 0; i < NumGas; i++)
    rows[i] = fld_offsets[ThisTask] + i;

  HYPRE_IJVectorGetValues(*x, NumGas, rows, XVec);

  for(i = 0; i < NumGas; i++)
    SphP[i].Utherm = XVec[i];

  myfree(XVec);
  myfree(rows);

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
                  tmpInExch[off].SurfaceArea = SphP[place].SurfaceArea;
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

                  for(init = 0; init < MaxFluxExch; init++)
                    {
                      tmpFluxExch[off].F2[init] = 0.0;
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

                      if(VF[vf].doface == 1)
                        {



                          if(DP[dp].task != ThisTask)
                            {

                              Noc++;

                              if(Noc >= MaxFluxExch)
                                terminate("Found a Cell with more than %d number of neighbours on other tasks\n", MaxFluxExch);

                              tmpFluxExch[off].ID[Noc] = DP[dp].ID;
                              if(nindex == P[place].ID)
                                {
                                  tmpFluxExch[off].alss[Noc] = VF[vf].akssp1;
                                  tmpFluxExch[off].F2[Noc] = VF[vf].F2p1;
                                }
                              else
                                {
                                  tmpFluxExch[off].alss[Noc] = VF[vf].akssp2;
                                  tmpFluxExch[off].F2[Noc] = VF[vf].F2p2;
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

#endif /*MULTIPLE_TIME_STEPPING*/
#endif /*SEMI_IMPLICIT_TI */
#endif /*MONOTONE_CONDUCTION */
