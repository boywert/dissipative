/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/mm.c
 * \date        MM/YYYY
 * \author      R. Grand
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

#ifdef TURBULENT_METALDIFFUSION

/*! \file mm.c
*  \brief Computes turbulent metal diffusion on an implicit diffusion solver
*/



#define cm (All.HubbleParam/All.UnitLength_in_cm)
#define g  (All.HubbleParam/All.UnitMass_in_g)
#define s  (All.HubbleParam/All.UnitTime_in_s)
#define erg (g*cm*cm/(s*s))
#define keV (1.602e-9*erg)
#define deg 1.0
#define m_p (PROTONMASS * g)
#define k_B (BOLTZMANN * erg / deg)

#define MAX_MM_ITER 400
#define MM_ITER_ACCURACY 1.0e-10

static struct inexch
{
  double in;
  double Kappa;
  double Volume;
  double Mass;

}
 *InExch;

static double *mfield, *mfield_old;
static double *Residual, *DVec, *QVec;
static MyFloat *Kappa;



void mm_gradient(double *in)
{
  int i;
  double volume_i;
  double in_i = 0.0;
  double in_k = 0.0;
  double grad_i[3] = { 0, 0, 0 };

  face *VF = Mesh.VF;
  point *DP = Mesh.DP;

  InExch = (struct inexch *) mymalloc("InExch", Mesh_nimport * sizeof(struct inexch));

  mm_exchange_in_vector(in);

  for(i = 0; i < NumGas; i++)
    {
      if(P[i].Type == 0)
        {
          in_i = in[i];
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
                in_k = in[particle_k];
              else
                in_k = InExch[particle_k].in;

              double n[3], c[3];
              n[0] = P[i].Pos[0] - DP[dp_k].x;
              n[1] = P[i].Pos[1] - DP[dp_k].y;
              n[2] = P[i].Pos[2] - DP[dp_k].z;
              double nn = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);

              c[0] = VF[vf_k].cx - 0.5 * (P[i].Pos[0] + DP[dp_k].x);
              c[1] = VF[vf_k].cy - 0.5 * (P[i].Pos[1] + DP[dp_k].y);
              c[2] = VF[vf_k].cz - 0.5 * (P[i].Pos[2] + DP[dp_k].z);

              grad_i[0] += VF[vf_k].area * (c[0] * (in_k - in_i) - 0.5 * (in_k + in_i) * n[0]) / nn;
              grad_i[1] += VF[vf_k].area * (c[1] * (in_k - in_i) - 0.5 * (in_k + in_i) * n[1]) / nn;
              grad_i[2] += VF[vf_k].area * (c[2] * (in_k - in_i) - 0.5 * (in_k + in_i) * n[2]) / nn;

              if(q_k == SphP[i].last_connection)
                break;

              q_k = DC[q_k].next;
            }

          grad_i[0] /= volume_i;
          grad_i[1] /= volume_i;
          grad_i[2] /= volume_i;

          SphP[i].gradm[0] = grad_i[0];
          SphP[i].gradm[1] = grad_i[1];
          SphP[i].gradm[2] = grad_i[2];
        }
    }

  myfree(InExch);
}


void calculate_shear_tensor(void)
{

  int i,k,l;

  for(i = 0; i < NumGas; i++)
    {
      if(P[i].Type == 0)
        {
	  /* Calculate the shear tensor */
	  for(l = 0; l < 3; l++)
	    {
	      for(k = 0; k < 3; k++)
		{
		  SphP[i].st[k][l] = 0.5 * ( SphP[i].Grad.dvel[k][l] + SphP[i].Grad.dvel[l][k] );
		  /*if(i<10){
		    mpi_printf("\nSphP.st = %g\n",SphP[i].st[k][l]);
		    mpi_printf("\nSphP.graddvel = %g\n",SphP[i].Grad.dvel[k][l]);
		    }*/
		}
	    }
	  
	}
    }

}



void turbulent_metal_mixing(void)
{

  int i,iter,l,k;
  double delta0, delta1, alpha, beta, rel_change, loc_max_rel_change, glob_max_rel_change;
  double sumnew, sumold, sumtransfer, sumnew_tot, sumold_tot, sumtransfer_tot, dm;
  double tr;
  double tracefreenorm;
  double cfac = 0.05;
  double stmp[3][3];

  if(All.dt_metaldiff==0)
    return;

  All.metaldiff_kappamax = 0.;

  Kappa = (MyFloat *) mymalloc("Kappa", NumGas * sizeof(MyFloat));

  calculate_shear_tensor();

  for(i = 0; i < NumGas; i++)
    {
      if(P[i].Type == 0)
        {

	  tracefreenorm = 0.;
	  tr = 0.;
	  for(l = 0; l < 3; l++)
	    tr += SphP[i].st[l][l];

	  tr /= 3.;

	  for(l = 0; l < 3; l++)
	    {

	      for(k = 0; k < 3; k++)
		{
		  stmp[k][l] = 0.5 * (SphP[i].st[k][l] + SphP[i].st[l][k]);

		    if(k == l)
		      stmp[k][l] -= tr;
		}
	    }

	  for(l = 0; l < 3; l++)
	    {
	      for(k = 0; k < 3; k++)
		{
		  SphP[i].st[k][l] = stmp[k][l];

		  tracefreenorm += SphP[i].st[k][l]*SphP[i].st[k][l];
		}
	    }

	  Kappa[i] = pow(tracefreenorm, 0.5) * pow(SphP[i].Volume * ( 3. / 4. * M_PI ), 2./3.) * cfac;
	  if(Kappa[i] > All.metaldiff_kappamax)
	    All.metaldiff_kappamax = Kappa[i];

	}
    }

  /*mpi_printf("\nAll.dt_metaldiff=%g\n",All.dt_metaldiff);
    mpi_printf("Hubble fac=%g\n",(All.cf_atime / All.cf_time_hubble_a));*/

  for(i = 0; i < NumGas; i++) // fold dt into difusion coefficient
    {
      if(P[i].Type == 0)
        {
	  Kappa[i] *= ( All.dt_metaldiff * All.cf_atime / All.cf_time_hubble_a );
	  /*All.dt_metaldiff = 0.1;// test values
	    Kappa[i] = 0.2;*/
	}
    }

  for(k=1; k<GFM_N_CHEM_ELEMENTS; k++)
    {
      tr=0.;
      tracefreenorm=0.;

      mfield = (double *) mymalloc("mfield", NumGas * sizeof(double)); 
      mfield_old =  (double *) mymalloc("mfield_old", NumGas * sizeof(double));
      Residual = (double *) mymalloc("Residual", NumGas * sizeof(double));
      DVec = (double *) mymalloc("DVec", NumGas * sizeof(double));
      QVec = (double *) mymalloc("QVec", NumGas * sizeof(double));
      double *diag;
      diag = (double *) mymalloc("diag", NumGas * sizeof(double));

      for(i = 0; i < NumGas; i++)
	{
	  if(P[i].Type == 0)
	    {
	      mfield[i] = mfield_old[i] = SphP[i].MassMetals[k];
	    }
	}

      mm_gradient(mfield);

      double *mfield_old_lhs;
      mfield_old_lhs = (double *) mymalloc("mfield_old_lhs", NumGas * sizeof(double));

      double theta = 0.5;

      mm_matrix_multiply(mfield_old, mfield_old_lhs, diag, -theta);

      mm_matrix_multiply(mfield_old, Residual, diag, theta);

      for(i = 0; i < NumGas; i++)
	{
	  if(P[i].Type == 0)
	    {
	      Residual[i] = mfield_old_lhs[i] - Residual[i];
	      
	      DVec[i] = Residual[i];
	    }
	}

      delta1 = mm_vector_multiply(Residual, Residual);

      delta0 = delta1;

      iter = 0;                     /* iteration counter */
      glob_max_rel_change = 1 + MM_ITER_ACCURACY;

      while(iter < 200 && glob_max_rel_change > 1e-8 && delta1 > 0)
	{
	  mm_matrix_multiply(DVec, QVec, diag, theta);

	  alpha = delta1 / mm_vector_multiply(DVec, QVec);

	  for(i = 0, loc_max_rel_change = 0; i < NumGas; i++)
	    {
	      if(mfield[i] == 0.)
		continue;

	      mfield[i] += alpha * DVec[i];
	      
	      Residual[i] -= alpha * QVec[i];

	      rel_change = alpha * DVec[i] / mfield[i];
	   
	      if(loc_max_rel_change < rel_change)
		loc_max_rel_change = rel_change;
	    }

	  delta0 = delta1;

	  delta1 = mm_vector_multiply(Residual, Residual);

	  beta = delta1 / delta0;

	  for(i = 0; i < NumGas; i++)
	    DVec[i] = Residual[i] + beta * DVec[i];

	  iter++;

	  MPI_Allreduce(&loc_max_rel_change, &glob_max_rel_change, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	}
      mpi_printf("TURBULENT_METALDIFFUSION: Finished loop with: iter=%d  delta1=%g delta1/delta0=%g  max-rel-change=%g\n", iter, delta1, delta1 / delta0, glob_max_rel_change);


      for(i = 0, sumnew = sumold = sumtransfer = 0; i < NumGas; i++)
	{
	  if(P[i].Type == 0)
	    {
	      sumnew += P[i].Mass * mfield[i];        
	      sumold += P[i].Mass * mfield_old[i];    
	      sumtransfer += P[i].Mass * fabs(mfield[i] - mfield_old[i]);
	      dm = mfield[i] - SphP[i].MassMetals[k];
	      SphP[i].MassMetals[k] += dm; 
	      if(SphP[i].MassMetals[k] < 0.)
		SphP[i].MassMetals[k] = 1e-8;
	      SphP[i].MetalsFraction[k] = SphP[i].MassMetals[k] / P[i].Mass;
	    }
	}

      MPI_Allreduce(&sumnew, &sumnew_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&sumold, &sumold_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&sumtransfer, &sumtransfer_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      mpi_printf("TURBULENT_METALDIFFUSION: finished for element %d. mass_before=%g, mass_after=%g, rel-change=%g rel-transfer=%g\n",k,sumold_tot, sumnew_tot, (sumnew_tot - sumold_tot) / sumold_tot, sumtransfer_tot / sumold_tot);


      myfree(mfield_old_lhs);
      myfree(diag);

      myfree(QVec);
      myfree(DVec);
      myfree(Residual);
      myfree(mfield_old);
      myfree(mfield);

    }

  myfree(Kappa);

}


double mm_vector_multiply(double *a, double *b)
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


void mm_matrix_multiply(double *in, double *out, double *diag, double theta)
{

  mm_matrix_multiply_isotropic(in, out, diag, theta);

}


void mm_matrix_multiply_isotropic(double *in, double *out, double *diag, double theta)
{
  int i;
  double in_i, in_j;
  double dx_ij, dy_ij, dz_ij, r_ij;
  double Kappa_i, Kappa_j;
  double kappa_mean_ij = 0.0;
  double area_ij = 0.0;

  InExch = (struct inexch *) mymalloc("InExch", Mesh_nimport * sizeof(struct inexch));

  mm_exchange_in_vector(in);

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

                  alpha += (kappa_mean_ij * area_ij) / r_ij;
                  double w = (kappa_mean_ij * area_ij) / r_ij; 
                  out_sum += w * (in_j - in_i);
                }

              if(q == SphP[i].last_connection)
                break;

              q = DC[q].next;
            }


          out[i] = P[i].Mass * in[i] - theta * out_sum;
          diag[i] = 1.0 / (P[i].Mass + 0.5 * alpha);
        }
    }                           /* for loop ends */

  myfree(InExch);
}




void mm_exchange_in_vector(double *in)
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
                  tmpInExch[off].in = in[place];
                  tmpInExch[off].Volume = SphP[place].Volume;
                  tmpInExch[off].Mass = P[place].Mass;
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
#endif
