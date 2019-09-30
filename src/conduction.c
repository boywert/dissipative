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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "voronoi.h"
#include "proto.h"

/*! \file conduction.c
*  \brief Computes isotropic/anisotropic conduction bases on an implicit diffusion solver
*/

#ifdef CONDUCTION

#if !defined(CONDUCTION_ANISOTROPIC) && !defined(CONDUCTION_ISOTROPIC)
#error CONDUCTION requires either CONDUCTION_ISOTROPIC or CONDUCTION_ANISOTROPIC
#endif

#ifdef CONDUCTION_ISOTROPIC
#if defined(CONDUCTION_ANISOTROPIC)
#error CONDUCTION_ISOTROPIC cannot be used with CONDUCTION_ANISOTROPIC
#endif
#endif

#ifdef CONDUCTION_ANISOTROPIC
#if !defined(MHD)
#error CONDUCTION_ANISOTROPIC requires MHD
#endif
#endif


#define cm (All.HubbleParam/All.UnitLength_in_cm)
#define g  (All.HubbleParam/All.UnitMass_in_g)
#define s  (All.HubbleParam/All.UnitTime_in_s)
#define erg (g*cm*cm/(s*s))
#define keV (1.602e-9*erg)
#define deg 1.0
#define m_p (PROTONMASS * g)
#define k_B (BOLTZMANN * erg / deg)

#define MAX_COND_ITER 400
#define COND_ITER_ACCURACY 1.0e-10

static struct inexch
{
  double in;
  double Kappa;
  double Volume;
  double Mass;

#ifdef CONDUCTION_ANISOTROPIC
  double B[3];
  double gradU[3];
  double res_sym;
  double egy_sym;
#endif
}
 *InExch;

static double *Temp, *TempOld;
static double *Residual, *DVec, *QVec;
static MyFloat *Kappa;
static unsigned short PCG_bool = 0;



/* 
 * we will use the conjugate gradient method to compute a solution
 * of the implicitly formulate diffusion equation 
 * 
 * Note: the conduction equation we solve is really formulated with u instead of T, i.e.
 * the factor (gamma-1)*mu*mp/k_B that converts from T to u is implicitely absorbed in a
 * redefinition of kappa 
 */

void conduction_exchange_in_vector(double *in);
void conduction_temperature_gradient_calculation(double *in);


void init_spitzer_conductivity(void)
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


void conduction(void)
{
  int i, iter;
  double delta0, delta1, alpha, beta, dt, rel_change, loc_max_rel_change, glob_max_rel_change;
  double sumnew, sumold, sumtransfer, sumnew_tot, sumold_tot, sumtransfer_tot, du;

#ifdef CONDUCTION_SATURATION
  double electron_free_path, temp_scale_length;
#endif

  mpi_printf("Start thermal conduction...\n");

  Temp = (double *) mymalloc("Temp", NumGas * sizeof(double));  /* [Temp] = erg/g    */
  TempOld = (double *) mymalloc("TempOld", NumGas * sizeof(double));    /* [TempOld] = erg/g */
  Residual = (double *) mymalloc("Residual", NumGas * sizeof(double));
  DVec = (double *) mymalloc("DVec", NumGas * sizeof(double));
  QVec = (double *) mymalloc("QVec", NumGas * sizeof(double));
  Kappa = (MyFloat *) mymalloc("Kappa", NumGas * sizeof(MyFloat));

#ifdef CONDUCTION_CRANK_NICOLSON
  double *TempOld_lhs;
  TempOld_lhs = (double *) mymalloc("TempOld_lhs", NumGas * sizeof(double));
#endif

#ifdef CONDUCTION_PCG
  double *Helpvar;
  Helpvar = (double *) mymalloc("Helvar", NumGas * sizeof(double));
#endif

  /* Actually diag is not needed for the CG-method without Jacobipreconditioner      */
  double *diag;
  diag = (double *) mymalloc("diag", NumGas * sizeof(double));


  dt = (All.Conduction_Ti_endstep - All.Conduction_Ti_begstep) * All.Timebase_interval;
  dt *= All.cf_atime / All.cf_time_hubble_a;


  mpi_printf("dt=%g\n", dt);

  /* First, let's compute the thermal energies per unit mass and conductivities for all particles */
  TIMER_START(CPU_CONDUCTION);

  for(i = 0; i < NumGas; i++)
    {
      if(P[i].Type == 0)
        {
          /* this gives the thermal energy per unit mass for particle i */
          Temp[i] = TempOld[i] = SphP[i].Utherm;

#ifdef CONDUCTION_CONSTANT
          Kappa[i] = All.ConductionCoeff;
#else
          Kappa[i] = All.ConductionCoeff * pow(TempOld[i], 2.5);
#ifdef CONDUCTION_SATURATION
          electron_free_path = All.ElectronFreePathFactor * Temp[i] * Temp[i] / (SphP[i].Density * All.cf_a3inv);
          temp_scale_length = All.Time * fabs(SphP[i].Utherm) / sqrt(SphP[i].Grad.dutherm[0] * SphP[i].Grad.dutherm[0] +
                                                                     SphP[i].Grad.dutherm[1] * SphP[i].Grad.dutherm[1] + SphP[i].Grad.dutherm[2] * SphP[i].Grad.dutherm[2]);

          Kappa[i] /= (1 + 4.2 * electron_free_path / temp_scale_length);
#endif /* CONDUCTION_SATURATION */
#endif /* CONDUCTION_CONSTANT */

#ifdef SFR
          if(SphP[i].d.Density * All.cf_a3inv >= All.PhysDensThresh)
            Kappa[i] = 0;
#endif /* SFR */

          /* we'll factor the timestep into the conductivities, for simplicity */
          Kappa[i] *= dt;
        }
    }
  /* Let's start the Conjugate Gradient Algorithm */
  /*
   *    CG Algorithmus for solving A * x = b
   *    or with Crank-Nicolson:          A * x = G * b = c
   *                                                                    x = u^{n+1}
   *                                                                    b = u^{n}
   *                                    
   *    Initialvalue:   x_(0) \in \R^n
   *                                    r_(0) := b - A*x_(0)
   *                                    d_(0) := r_(0)
   *
   *  For i >= 0:               
   *                            alpha_i = |r_(i)|^2 / <d_(i), A*d_(i)>
   *                            x_(i+1) = x_(i) + alpha_i * d_(i)
   *                            r_(i+1) = r_(i) - alpha_i * A*d_(i)
   *
   *                            beta_(i) = |r_(i+1)|^2 / |r_(i)|^2
   *                            d_(i+1) = r_(i+1) + beta_(i) * d_(i)
   */

  /*
   *    Preconditioner CG Algorithmus for solving A * x = b
   *    or with Crank-Nicolson:          A * x = G * b = c
   *                                                                    x = u^{n+1}
   *                                                                    b = u^{n}
   *                                    
   *    Initialvalue:   x_(0) \in \R^n
   *                                    r_(0) := b - A*x_(0)
   *                                    h_(0) := C^-1 * r_0
   *                                    d_(0) := h_(0)
   *
   *  For i >= 0:               
   *                            alpha_i = < r_(i), h_(i) > / <d_(i), A*d_(i)>
   *                            x_(i+1) = x_(i) + alpha_i * d_(i)
   *                            r_(i+1) = r_(i) - alpha_i * A*d_(i)
   *
   *                                    h_(i+1) = C^-1 * r_(i+1)                                        //this is new
   *                                    
   *                            beta_(i) = < r_(i+1), h_(i+1) > / < r_(i), h_(i) >
   *                            d_(i+1) = h_(i+1) + beta_(i) * d_(i)
   */

  /* Initialization */
#ifdef CONDUCTION_CRANK_NICOLSON
  double theta = 0.5;
#else
  double theta = 1.0;
#endif /* CONDUCTION_CRANK_NICOLSON */

#ifdef CONDUCTION_ANISOTROPIC
  conduction_temperature_gradient_calculation(TempOld);
#endif /* CONDUCTION_ANISOTROPIC */

#ifdef CONDUCTION_CRANK_NICOLSON
  /*      [TempOld_lhs] = erg     
   *              [TempOld] = erg / g
   *              [A] = g
   *              G * b = c
   */
  conduction_matrix_multiply(TempOld, TempOld_lhs, diag, -theta);
#endif /* CONDUCTION_CRANK_NICOLSON */

  conduction_matrix_multiply(TempOld, Residual, diag, theta);

  for(i = 0; i < NumGas; i++)
    {
      if(P[i].Type == 0)
        {
#ifdef CONDUCTION_CRANK_NICOLSON
          /*  r_(0) = c - A * x_(0)   = G * b - Ax_(0)
           * [Residual] = erg
           */
          Residual[i] = TempOld_lhs[i] - Residual[i];
#else
          /*  r_(0) = b - A * x_(0)   */
          /* [Residal] = erg 
           */
          Residual[i] = P[i].Mass * TempOld[i] - Residual[i];
#endif /* CONDUCTION_CRANK_NICOLSON */

#ifdef CONDUCTION_PCG
          /*  h_0 = C^-1 * r_0                    */
          /*  Jacobi Preconditioner:   [Helpvar] = [C^-1]*[Residual] = erg / g         */

          /*
           * you have two options for using the preconditioner: *
           *    a)      you can use PCG_bool and implement in conduction_matrix_multiply the
           *                    operation h_(i+1) = C^{-1} * r_(i+1)
           *
           *    b)      you can use conduction_PCG(Residual, Helpvar) with a extra function,
           *                    where you can implement a more sophisticated preconditioner.
           *                     
           *    c)      At the moment the Jacobi-preconditioner is available. 
           *                    The function conduction_matrix_multiply calculate the diagonal elements
           *                    of the Matrix A. The Jacobi-preconditioner is C^{-1}_{ij} = \delta_{ij} / A_ij.
           *                    
           *
           *    At the moment option c) is used in case of CONDUCTION_PCG on.
           *
           */

          /*    option a)
           *    PCG_bool = 1;
           *    conduction_matrix_multiply(Residual, Helpvar, theta);
           *    PCG_bool = 0;
           */

          /*
           *    option b)       
           *    conduction_PCG(Residual, Helpvar);
           */

          /* [DVec] = erg / g     */
          Helpvar[i] = Residual[i] * diag[i];
          DVec[i] = Helpvar[i];
#else
          /*      d_(0) = b - A * x_(0) or d_(0) = c - A * x_(0)  */
          /*      [DVec] = erg                                    */
          DVec[i] = Residual[i];
#endif /* CONDUCTION_PCG */
        }
    }

#ifdef CONDUCTION_PCG
  /*      delta1 = <r_(i), h_(i)>         */
  /*      [delta1] = erg * erg / g        */
  delta1 = conduction_vector_multiply(Residual, Helpvar);
#else
  /*      delta1 = |r_(i)|^2                      */
  /*      [delta1] = erg * erg    */
  delta1 = conduction_vector_multiply(Residual, Residual);
#endif /* CONDUCTION_PCG */

  delta0 = delta1;              /* [delta0] = [delta1] */

  iter = 0;                     /* iteration counter */
  glob_max_rel_change = 1 + COND_ITER_ACCURACY; /* to make sure that we enter the iteration */

  /*      (P)CG Loop      */
  while(iter < MAX_COND_ITER && glob_max_rel_change > COND_ITER_ACCURACY && delta1 > 0)
    {
      /*      DVec[i] = d_(i) and QVec -> A * DVec = A * d_(i)        */
      /*      CG:     [DVec] = erg,           [QVec] = erg *g,        [A_ij] = g      */
      /*      PCG:    [DVec] = erg/g,         [QVec] = erg,           [A_ij] = g      */
      conduction_matrix_multiply(DVec, QVec, diag, theta);

      /*      CG:     [alpha] = 1/g   or PCG: [alpha] = 1     */
      alpha = delta1 / conduction_vector_multiply(DVec, QVec);

      for(i = 0, loc_max_rel_change = 0; i < NumGas; i++)
        {
          /*      x_(i+1) = x_(i) + alpha_i * d_(i)               */
          /* CG:  [Temp] = [alpha]*[DVec] = 1/g * erg = erg/g     */
          /* PCG: [Temp] = [alpha]*[DVec] = 1 * erg/g = erg/g     */
          Temp[i] += alpha * DVec[i];

          /*  r_(i+1) = r_(i) - alpha_i * Ad_(i)          */
          /*  CG:     [Residual] = 1/g*erg*g = erg        */
          /*  PCG:    [Residual] = 1 * erg = erg          */
          Residual[i] -= alpha * QVec[i];

          rel_change = alpha * DVec[i] / Temp[i];

          if(loc_max_rel_change < rel_change)
            loc_max_rel_change = rel_change;
        }

      /* delta_0 = <r_(i), r_(i) >   or  delta_0 = <r_(i), h_(i) > */
      delta0 = delta1;

#ifdef CONDUCTION_PCG
      /*      h_(i+1) = C^{-1} * r_(i+1)      */

      /* PCG_bool = 1;
       * conduction_matrix_multiply(Residual, Helpvar, theta);
       * PCG_bool = 0;
       */
      for(i = 0; i < NumGas; i++)
        Helpvar[i] = Residual[i] * diag[i];

      /*      
       * option b)
       * conduction_PCG(Residual, Helpvar);
       */

      /* delta1 = <r_(i+1), h_(i+1)>  */
      delta1 = conduction_vector_multiply(Residual, Helpvar);
#else
      /* delta1 = <r_(i+1), r_(i+1)>  */
      delta1 = conduction_vector_multiply(Residual, Residual);
#endif /* CONDUCTION_PCG */

      beta = delta1 / delta0;   /* [beta] = 1 */

      /*      d_(i+1) = r_(i+1) + beta_(i) * d_(i)    */
      /*      [DVec] = erg                            */
      for(i = 0; i < NumGas; i++)
#ifdef	CONDUCTION_PCG
        DVec[i] = Helpvar[i] + beta * DVec[i];
#else
        DVec[i] = Residual[i] + beta * DVec[i];
#endif /* CONDUCTION_PCG */

      iter++;

      MPI_Allreduce(&loc_max_rel_change, &glob_max_rel_change, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

      mpi_printf("conduction iter=%d  delta1=%g delta1/delta0=%g  max-rel-change=%g\n", iter, delta1, delta1 / delta0, glob_max_rel_change);
    }
  /* Now we have the solution vector in Temp[] */
  /* assign it to the entropies, and update the pressure */
  for(i = 0, sumnew = sumold = sumtransfer = 0; i < NumGas; i++)
    {
      if(P[i].Type == 0)
        {
          sumnew += P[i].Mass * Temp[i];        /*  [sumnew] = erg  */
          sumold += P[i].Mass * TempOld[i];     /*  [sumold] = erg  */
          sumtransfer += P[i].Mass * fabs(Temp[i] - TempOld[i]);
          du = Temp[i] - SphP[i].Utherm;        /* [Temp] = erg/g   */
          SphP[i].Utherm += du; /* [Utherm] = erg/g   */
          SphP[i].Energy += All.cf_atime * All.cf_atime * du * P[i].Mass;       /* [Energy] = erg     */
          set_pressure_of_cell(i);
        }
    }

  MPI_Allreduce(&sumnew, &sumnew_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&sumold, &sumold_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&sumtransfer, &sumtransfer_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  TIMER_STOP(CPU_CONDUCTION);

  mpi_printf("\nconduction finished. energy_before=%g erg, energy_after=%g erg, rel-change=%g rel-transfer=%g\n\n",
             sumold_tot, sumnew_tot, (sumnew_tot - sumold_tot) / sumold_tot, sumtransfer_tot / sumold_tot);

  myfree(diag);

#ifdef CONDUCTION_PCG
  myfree(Helpvar);
#endif

#ifdef CONDUCTION_CRANK_NICOLSON
  myfree(TempOld_lhs);
#endif

  myfree(Kappa);
  myfree(QVec);
  myfree(DVec);
  myfree(Residual);
  myfree(TempOld);
  myfree(Temp);
}


double conduction_vector_multiply(double *a, double *b)
{
  /* scalar product of vector a and b */
  int i;
  double sum, sumall;

  for(i = 0, sum = 0; i < NumGas; i++)
    if(P[i].Type == 0)
      sum += a[i] * b[i];

  MPI_Allreduce(&sum, &sumall, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return sumall;
}


#ifdef CONDUCTION_PCG
void conduction_PCG(double *in, double *out)
{
  /*  h = C^-1 * r */
#ifdef CONDUCTION_ISOTROPIC
  conduction_PCG_isotropic(in, out);
#endif

#ifdef CONDUCTION_ANISOTROPIC
  conduction_PCG_anisotropic(in, out);
#endif
}

#ifdef CONDUCTION_ISOTROPIC
void conduction_PCG_isotropic(double *in, double *out)
{
  /*      Here you can implement an arbitrary Preconditioner
   *      We use for the first time the Jacobi Preconditioner 
   *      A = D + L + R
   *      C = D   => C^{-1}_ii = D^{-1}_ii
   *      
   *      out = C^{-1} * in
   *
   */
  terminate("here you can insert a more complicated preconditioner\n");

  /*
     int i;
     double Kappa_i, Kappa_j;
     double kappa_mean_ij;
     double alpha = 0.0;
     double dx_ij, dy_ij, dz_ij;
     double r_ij;
     double area_ij;

     for(i = 0; i < NumGas; i++)
     {
     if(P[i].Type == 0)
     {

     Kappa_i = Kappa[i];

     int q = SphP[i].first_connection;
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

     if(particle >= NumGas && Mesh.DP[dp].task == ThisTask)
     particle -= NumGas;
     if(Mesh.DP[dp].task == ThisTask)
     {
     Kappa_j = Kappa[particle];
     }
     else
     {
     Kappa_j = InExch[particle].Kappa;
     }
     dx_ij = P[i].Pos[0] - Mesh.DP[dp].x;
     dy_ij = P[i].Pos[1] - Mesh.DP[dp].y;
     dz_ij = P[i].Pos[2] - Mesh.DP[dp].z;

     r_ij = sqrt(dx_ij * dx_ij + dy_ij * dy_ij + dz_ij * dz_ij);
     area_ij = Mesh.VF[vf].area;

     if((Kappa_i + Kappa_j) > 0)
     {
     kappa_mean_ij = 2.0 * (Kappa_i * Kappa_j) / (Kappa_i + Kappa_j);
     alpha += (kappa_mean_ij * area_ij) / r_ij; // [alpha] = g 
     }

     if(q == SphP[i].last_connection)
     break;

     q = DC[q].next;
     }

     #ifdef CONDUCTION_CRANK_NICOLSON
     out[i] = in[i] / ( P[i].Mass + 0.5 * alpha  );
     #else
     out[i] = in[i] / ( P[i].Mass +  alpha  );
     #endif

     }  // if ends 
     }
   */
}
#endif /* CONDUCTION_ISOTROPIC */

#ifdef CONDUCTION_ANISOTROPIC
void conduction_PCG_anisotropic(double *in, double *out)
{
  terminate("here you can insert a more complicated preconditioner\n");
}
#endif /* CONDUCTION_ANISOTROPIC */

#endif /* CONDUCTION_PCG */


void conduction_matrix_multiply(double *in, double *out, double *diag, double theta)
{
#ifdef CONDUCTION_ISOTROPIC
  conduction_matrix_multiply_isotropic(in, out, diag, theta);
#endif

#ifdef CONDUCTION_ANISOTROPIC
  conduction_matrix_multiply_anisotropic(in, out, diag, theta);
#endif
}


#ifdef CONDUCTION_ISOTROPIC
void conduction_matrix_multiply_isotropic(double *in, double *out, double *diag, double theta)
{
  /* [in] = erg/g, [out] = [in] * g = erg    */
  int i;
  double in_i, in_j;
  double dx_ij, dy_ij, dz_ij, r_ij;
  double Kappa_i, Kappa_j;
  double kappa_mean_ij = 0.0;
  double area_ij = 0.0;

  InExch = (struct inexch *) mymalloc("InExch", Mesh_nimport * sizeof(struct inexch));

  conduction_exchange_in_vector(in);

  for(i = 0; i < NumGas; i++)
    {
      if(P[i].Type == 0)
        {
          Kappa_i = Kappa[i];
          in_i = in[i];

          double out_sum = 0.0;
          double alpha = 0.0;
          int q = SphP[i].first_connection;

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

              if(particle >= NumGas && Mesh.DP[dp].task == ThisTask)
                particle -= NumGas;

              if(Mesh.DP[dp].task == ThisTask)
                {
                  Kappa_j = Kappa[particle];
                  in_j = in[particle];
                }
              else
                {
                  Kappa_j = InExch[particle].Kappa;
                  in_j = InExch[particle].in;
                }

              dx_ij = P[i].Pos[0] - Mesh.DP[dp].x;
              dy_ij = P[i].Pos[1] - Mesh.DP[dp].y;
              dz_ij = P[i].Pos[2] - Mesh.DP[dp].z;
              r_ij = sqrt(dx_ij * dx_ij + dy_ij * dy_ij + dz_ij * dz_ij);

              area_ij = Mesh.VF[vf].area;

              if((Kappa_i + Kappa_j) > 0)
                {
                  kappa_mean_ij = 2.0 * (Kappa_i * Kappa_j) / (Kappa_i + Kappa_j);

                  /*if(PCG_bool == 1)     {
                     alpha += (kappa_mean_ij * area_ij) / r_ij;
                     }
                     else 
                     {
                     double w = (kappa_mean_ij * area_ij) / r_ij;  //[w] = g      
                     out_sum += w * (in_j - in_i);
                     } */

                  alpha += (kappa_mean_ij * area_ij) / r_ij;
                  double w = (kappa_mean_ij * area_ij) / r_ij;  /* [w] = g */
                  out_sum += w * (in_j - in_i);
                }

              if(q == SphP[i].last_connection)
                break;

              q = DC[q].next;
            }

          /*    if(PCG_bool == 1) {
             //out[i] = in[i] / (P[i].Mass + theta * alpha);
             out[i] = in[i];
             }
             else{
             out[i] = P[i].Mass * in[i] - theta * out_sum;
             }
             }  *//* if condition (P[i].Type == 0) ends   */

          out[i] = P[i].Mass * in[i] - theta * out_sum;
          diag[i] = 1.0 / (P[i].Mass + 0.5 * alpha);
        }
    }                           /* for loop ends */

  myfree(InExch);
}
#endif /* CONDUCTION_ISOTROPIC */

#ifdef CONDUCTION_ANISOTROPIC
void conduction_temperature_gradient_calculation(double *in)
{
  int i;
  double volume_i;
  double in_u_i = 0.0;
  double in_u_k = 0.0;
  double grad_i[3] = { 0, 0, 0 };

  face *VF = Mesh.VF;
  point *DP = Mesh.DP;

  InExch = (struct inexch *) mymalloc("InExch", Mesh_nimport * sizeof(struct inexch));

  conduction_exchange_in_vector(in);

  for(i = 0; i < NumGas; i++)
    {
      if(P[i].Type == 0)
        {
          in_u_i = in[i];
          volume_i = SphP[i].Volume;
          grad_i[0] = 0;
          grad_i[1] = 0;
          grad_i[2] = 0;

          int q_k = SphP[i].first_connection;

          while(q_k >= 0)
            {
              int dp_k = DC[q_k].dp_index;
              int vf_k = DC[q_k].vf_index;
              int particle_k = DP[dp_k].index;

              if(particle_k < 0)
                {
                  q_k = DC[q_k].next;
                  continue;
                }

              if(particle_k >= NumGas && DP[dp_k].task == ThisTask)
                particle_k -= NumGas;

              if(DP[dp_k].task == ThisTask)
                in_u_k = in[particle_k];
              else
                in_u_k = InExch[particle_k].in;

              double n[3], c[3];
              n[0] = P[i].Pos[0] - DP[dp_k].x;
              n[1] = P[i].Pos[1] - DP[dp_k].y;
              n[2] = P[i].Pos[2] - DP[dp_k].z;
              double nn = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);

              c[0] = VF[vf_k].cx - 0.5 * (P[i].Pos[0] + DP[dp_k].x);
              c[1] = VF[vf_k].cy - 0.5 * (P[i].Pos[1] + DP[dp_k].y);
              c[2] = VF[vf_k].cz - 0.5 * (P[i].Pos[2] + DP[dp_k].z);

              grad_i[0] += VF[vf_k].area * (c[0] * (in_u_k - in_u_i) - 0.5 * (in_u_k + in_u_i) * n[0]) / nn;
              grad_i[1] += VF[vf_k].area * (c[1] * (in_u_k - in_u_i) - 0.5 * (in_u_k + in_u_i) * n[1]) / nn;
              grad_i[2] += VF[vf_k].area * (c[2] * (in_u_k - in_u_i) - 0.5 * (in_u_k + in_u_i) * n[2]) / nn;

              if(q_k == SphP[i].last_connection)
                break;

              q_k = DC[q_k].next;
            }

          grad_i[0] /= volume_i;
          grad_i[1] /= volume_i;
          grad_i[2] /= volume_i;

          SphP[i].gradU[0] = grad_i[0];
          SphP[i].gradU[1] = grad_i[1];
          SphP[i].gradU[2] = grad_i[2];
        }
    }

  myfree(InExch);
}


void conduction_matrix_multiply_anisotropic(double *in, double *out, double *diag, double theta)
{
  int i;
  double B_i[3], B_k[3];
  double in_u_i, in_u_k;
  double Kappa_i, Kappa_k;
  double grad_i[3], grad_k[3];
  double kappa_mean_ik;
  double alpha = 0.0;

  point *DP = Mesh.DP;

  InExch = (struct inexch *) mymalloc("InExch", Mesh_nimport * sizeof(struct inexch));

  conduction_exchange_in_vector(in);

  /* prepare auxiliary arrays for symmetrization */
  double *res = mymalloc("res", NumGas * sizeof(double));

  for(i = 0; i < NumGas; i++)
    res[i] = 0;

  for(i = 0; i < NumGas; i++)
    {
      if(P[i].Type == 0)
        {
          in_u_i = in[i];
          Kappa_i = Kappa[i];

          B_i[0] = SphP[i].B[0];
          B_i[1] = SphP[i].B[1];
          B_i[2] = SphP[i].B[2];

          grad_i[0] = SphP[i].gradU[0];
          grad_i[1] = SphP[i].gradU[1];
          grad_i[2] = SphP[i].gradU[2];

          alpha = 0.0;

          /* do main loop over the neighbours of i */
          int q_k = SphP[i].first_connection;

          while(q_k >= 0)
            {
              int dp_k = DC[q_k].dp_index;      //Initialisation of Delaunay points
              int vf_k = DC[q_k].vf_index;      //Initialisation of Voronoi faces
              int particle_k = DP[dp_k].index;

              if(particle_k < 0)
                {
                  q_k = DC[q_k].next;
                  continue;
                }

              if(particle_k >= NumGas && DP[dp_k].task == ThisTask)
                particle_k -= NumGas;

              if(DP[dp_k].task == ThisTask)
                {
                  in_u_k = in[particle_k];
                  Kappa_k = Kappa[particle_k];

                  B_k[0] = SphP[particle_k].B[0];
                  B_k[1] = SphP[particle_k].B[1];
                  B_k[2] = SphP[particle_k].B[2];

                  grad_k[0] = SphP[particle_k].gradU[0];
                  grad_k[1] = SphP[particle_k].gradU[1];
                  grad_k[2] = SphP[particle_k].gradU[2];
                }
              else
                {
                  in_u_k = InExch[particle_k].in;
                  Kappa_k = InExch[particle_k].Kappa;

                  B_k[0] = InExch[particle_k].B[0];
                  B_k[1] = InExch[particle_k].B[1];
                  B_k[2] = InExch[particle_k].B[2];

                  grad_k[0] = InExch[particle_k].gradU[0];
                  grad_k[1] = InExch[particle_k].gradU[1];
                  grad_k[2] = InExch[particle_k].gradU[2];
                }

              /* get average B-field */
              double B[3];
              B[0] = 0.5 * (B_i[0] + B_k[0]);
              B[1] = 0.5 * (B_i[1] + B_k[1]);
              B[2] = 0.5 * (B_i[2] + B_k[2]);

              /* get average temperature gradient direction */
              double grad[3];
              grad[0] = 0.5 * (grad_i[0] + grad_k[0]);
              grad[1] = 0.5 * (grad_i[1] + grad_k[1]);
              grad[2] = 0.5 * (grad_i[2] + grad_k[2]);

              /* calculate outward normal */
              double n[3];
              n[0] = DP[dp_k].x - P[i].Pos[0];
              n[1] = DP[dp_k].y - P[i].Pos[1];
              n[2] = DP[dp_k].z - P[i].Pos[2];
              /* length */
              double nn = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);

              double area_ik = Mesh.VF[vf_k].area;

              if((Kappa_i + Kappa_k) > 0)
                kappa_mean_ik = 2.0 * (Kappa_i * Kappa_k) / (Kappa_i + Kappa_k);
              else
                kappa_mean_ik = 0.0;

              double grad_dot_B = (grad[0] * B[0] + grad[1] * B[1] + grad[2] * B[2]);
              double grad_dot_n = (grad[0] * n[0] + grad[1] * n[1] + grad[2] * n[2]);
              double B_dot_n = (B[0] * n[0] + B[1] * n[1] + B[2] * n[2]);
              double B_dot_B = (B[0] * B[0] + B[1] * B[1] + B[2] * B[2]);

              if(B_dot_B > MIN_FLOAT_NUMBER && fabs(grad_dot_n) > MIN_FLOAT_NUMBER)
                {
                  double A_ik = area_ik / nn * kappa_mean_ik * B_dot_n * grad_dot_B / (grad_dot_n * B_dot_B);

                  if(A_ik > 0)
                    alpha += A_ik;

                  /* A_ik < 0 can destroy monotonicity - therefore we use the cos^2 approximation */
                  if(A_ik < 0)
                    {
                      A_ik = area_ik / nn * kappa_mean_ik * B_dot_n * B_dot_n / (nn * nn);
                      alpha += A_ik;
                    }

                  if(PCG_bool == 1)
                    res[i] += A_ik;
                  else
                    res[i] += A_ik * (in_u_k - in_u_i);
                }

              if(q_k == SphP[i].last_connection)
                break;

              q_k = DC[q_k].next;
            }
        }
    }

  for(i = 0; i < NumGas; i++)
    {
      if(PCG_bool == 1)
        out[i] = in[i] / (P[i].Mass + theta * res[i]);
      else
        {
          out[i] = P[i].Mass * in[i] - theta * res[i];
          diag[i] = 1.0 / (P[i].Mass + 0.5 * alpha);
        }
    }

  myfree(res);
  myfree(InExch);
}
#endif /* CONDUCTION_ANISOTROPIC */

void conduction_exchange_in_vector(double *in)
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

                  tmpInExch[off].in = in[place];
                  tmpInExch[off].Kappa = Kappa[place];
                  tmpInExch[off].Volume = SphP[place].Volume;
                  tmpInExch[off].Mass = P[place].Mass;
#ifdef CONDUCTION_ANISOTROPIC
                  tmpInExch[off].B[0] = SphP[place].B[0];
                  tmpInExch[off].B[1] = SphP[place].B[1];
                  tmpInExch[off].B[2] = SphP[place].B[2];

                  tmpInExch[off].gradU[0] = SphP[place].gradU[0];
                  tmpInExch[off].gradU[1] = SphP[place].gradU[1];
                  tmpInExch[off].gradU[2] = SphP[place].gradU[2];
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
#endif /* CONDUCTION */
