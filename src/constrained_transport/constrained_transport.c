/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/constrained_transport/constrained_transport.c
 * \date        06/2015
 * \author      Philip Mocz
 * \brief       includes functions for constrained transport algorithm
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY 3D version of Constrained Transport of Delaunay triangles w/ vector potential 
 * - we evolve the variable Avec * voronoi cell radius (double AvecR[3])
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>

#include "../allvars.h"
#include "../proto.h"
#include "../voronoi.h"
#include "../domain.h"

#ifdef MHD
#ifdef MHD_CT

extern const int edge_start[6];
extern const int edge_end[6];
extern const int edge_opposite[6];
extern const int edge_nexttetra[6];


/** set ICs for A_periodic (A(x,t) = A_periodic(x,t) + A_meanfield(x)) determined by MHD_CT_IC */
void set_A_ICs(int ic_flag)
{
  int i, k;
  double A[3];
  for(i = 0; i < NumGas; i++)
    {
      get_init_A(SphP[i].Center[0], SphP[i].Center[1], SphP[i].Center[2], &A[0], &A[1], &A[2], ic_flag);
      for(k = 0; k < 3; k++)
        {
          SphP[i].A[k] = A[k];
          SphP[i].AConserved[k] = SphP[i].A[k] * SphP[i].Volume;
        }
      SphP[i].TimeLastBUpdate = All.Time;
    }
  correct_ctr_b();
}


/** IC helper - get initial A_periodic (vector potential) as function of position */
void get_init_A(double x, double y, double z, double *Ax, double *Ay, double *Az, int ic_flag)
{
  if(ic_flag == 0 || ic_flag == 1)
    {                           // B_non-mean-field = 0
      *Ax = 0.;
      *Ay = 0.;
      *Az = 0.;
    }
  else if(ic_flag == 2)
    {                           // Orszag-Tang
      *Ax = 0.;
      *Ay = 0.;
      *Az = cos(4 * M_PI * x) / (4 * M_PI * sqrt(4 * M_PI)) + cos(2 * M_PI * y) / (2 * M_PI * sqrt(4 * M_PI));
    }
  else if(ic_flag == 3)
    {                           // Field Loop
      *Ax = 0.;
      *Ay = 0.;
      *Az = 0.;
      if(sqrt((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5)) < 0.3)
        {
          *Az = 0.001 / sqrt(4. * M_PI) * (0.3 - sqrt((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5)));
        }
    }
  else if(ic_flag == 4)
    {                           // Alfven wave along z
      double k = 2. * M_PI;
      *Ax =  cos(k * z) / k;
      *Ay = -sin(k * z) / k;
      *Az =  0.;
    }
  else if(ic_flag == 5)
    {                           // 1D gaussian B field (for diffusion test only!)
      *Ax = 0.;
      *Az = 0.;

      //forces the periodicity of the vector potential
      if(x >= boxSize_X)
        x -= boxSize_X;

      if(x > 0.5 * boxSize_X)
        {
          double xprime = x - 0.75 * boxSize_X;
          *Ay = 0.5 * erf(xprime / sqrt(4 * All.TimeBegin));
        }
      else
        {
          double xprime = x - 0.25 * boxSize_X;
          *Ay = 0.5 * erf(-xprime / sqrt(4 * All.TimeBegin));
        }
    }
  else if(ic_flag == 6)
    {                           // 2D gaussian B field (for diffusion test only!)

      //forces the periodicity of the vector potential
      if(x >= boxSize_X)
        x -= boxSize_X;

      if(y >= boxSize_Y)
        y -= boxSize_Y;

      double xprime = x - 0.5 * boxSize_X;
      double yprime = y - 0.5 * boxSize_Y;
      double r2 = xprime*xprime + yprime*yprime;

      *Ax = -yprime * exp(-r2 / sqrt(4 * All.TimeBegin)) / (2 * M_PI * r2);
      *Ay =  xprime * exp(-r2 / sqrt(4 * All.TimeBegin)) / (2 * M_PI * r2);
      *Az =  0.;
    }
  else if(ic_flag == 7)
    {                           // Harris current sheet 
      double B_0 = 1.;
      double delta = 0.1;
      double off = log(cosh((0.25 * boxSize_Y) / delta)) + log(cosh((-0.25 * boxSize_Y / delta)));
      *Ax = 0.0;
      *Ay = 0.0;

      /* to ensure periodicity */
      if(y >= boxSize_Y)
        y -= boxSize_Y;

      if(y > 0.5 * boxSize_Y)
        {
          double yprime = y - 0.75 * boxSize_Y;
          *Az = B_0 * delta * log(cosh(yprime / delta));
        }
      else
        {
          double yprime = y - 0.25 * boxSize_Y;
          *Az = -B_0 * delta * (log(cosh(yprime / delta)) - off);
        }
    }
  else
    {
      terminate("Specified MHD_CT_IC non-existant. Code your MHD ICs here!\n");
    }
}




/** update A - right after finite_volume_solver is called; or in mhd.c 2-step Strange-split source term */
void update_A(int i, double dt)
{
  int k;
  double V[3], B[3], E[3];

  for(k = 0; k < 3; k++)
    {
      V[k] = P[i].Vel[k];
      B[k] = SphP[i].B[k];
    }

  // add in cosmological factors
  if(All.ComovingIntegrationOn)
    {                           // variables are rho_c, rho_cw, calE, B_c, w, uth, A_c: B_c = nabla_x cross A_c
      double atime = All.Time;
      // convert to peculiar velocity. (note: dt is already in cosmological units)
      for(k = 0; k < 3; k++)
        {
          V[k] /= atime;
        }
    }

  // calculate E field: E = - V x B
  E[0] = -(V[1] * B[2] - V[2] * B[1]);
  E[1] =  (V[0] * B[2] - V[2] * B[0]);
  E[2] = -(V[0] * B[1] - V[1] * B[0]);

  for(k = 0; k < 3; k++) {
    //SphP[i].A[k] += dt * -E[k];  // AConserved is the main variable being updated
    SphP[i].AConserved[k] += dt * -E[k] * SphP[i].Volume;
  }

}




#if !defined(TWODIMS) && !defined(MHD_CT_CTR_B_ALT)
/** get center B-fields from A - at the end of each main loop iteration */
void correct_ctr_b(void)
{
  int idx, pi;
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      pi = TimeBinsHydro.ActiveParticleList[idx];
      if(pi < 0)
        continue;

      int k;
      for(k = 0; k < 3; k++)
        {
          SphP[pi].A[k] = SphP[pi].AConserved[k] / SphP[pi].Volume;
        }
      SphP[pi].TimeLastBUpdate = All.Time;
    }

  exchange_primitive_variables();       // to ensure we have all neighboring A

  // recover the B-fields: A -> tetra face phi's -> tetra B's -> cell ctr B's

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      pi = TimeBinsHydro.ActiveParticleList[idx];
      if(pi < 0)
        continue;

      int edge = SphP[pi].first_connection;
      int last_edge = SphP[pi].last_connection;

      int vf = DC[edge].vf_index;
      int dp_index = Mesh.VF[vf].p1;
      if(dp_index == DC[edge].dp_index)
        dp_index = Mesh.VF[vf].p2;

      point *dp_ref = &(Mesh.DP[dp_index]);

      point dp_ref_cm;
      switch_to_CM(&dp_ref, &dp_ref_cm);

      double Bx_running = 0;
      double By_running = 0;
      double Bz_running = 0;
      double normalization = 0;
      double voronoi_face_area = Mesh.VF[vf].area;

      //loop over all the connections around a cell
      while(1)
        {
          int dq_index = DC[edge].dp_index;

          //one of the tetrahedras around the Delaunay connection
          int tt = DC[edge].dt_index;
          tetra *t = &Mesh.DT[tt];

          //find the local index of the edge
          int nr = 6;
          int e, dp_start_index, dp_end_index;

          for(e = 0; e < 6; e++)
            {
              dp_start_index = t->p[edge_start[e]];
              dp_end_index = t->p[edge_end[e]];
              if((dp_start_index == dp_index && dp_end_index == dq_index) || (dp_start_index == dq_index && dp_end_index == dp_index))
                {
                  nr = e;
                  break;
                }
            }

          //ensure that the local edge index has been found
          assert(nr != 6);

          //already set: t,tt,nr
          int i, j, k, l, m, ii, jj, kk, ll, nn;
          tetra *prev, *next;

          //tetra_center *nextc;
          i = edge_start[nr];
          j = edge_end[nr];
          k = edge_opposite[nr];
          l = edge_nexttetra[nr];

          prev = t;

          // loop over all the tetra around connection
          do
            {
              nn = prev->t[l];
              next = &Mesh.DT[nn];
              //nextc = &Mesh.DTC[nn];

              // calculate (outward) phi's of the tetra faces from the A's
              point *p0 = &(Mesh.DP[prev->p[0]]);
              point *p1 = &(Mesh.DP[prev->p[1]]);
              point *p2 = &(Mesh.DP[prev->p[2]]);
              point *p3 = &(Mesh.DP[prev->p[3]]);

              point p0_cpy, p1_cpy, p2_cpy, p3_cpy;
              switch_to_CM(&p0, &p0_cpy);       // this is better for refinement
              switch_to_CM(&p1, &p1_cpy);
              switch_to_CM(&p2, &p2_cpy);
              switch_to_CM(&p3, &p3_cpy);

              double phi[3];
              phi[0] = get_phi(p1, p2, p3, dp_ref);
              phi[1] = get_phi(p0, p3, p2, dp_ref);
              phi[2] = get_phi(p0, p1, p3, dp_ref);

              // calculate tetra B field & add contribution to cell
              double Area[3][3];
              get_face_area(p1, p2, p3, &Area[0][0], &Area[0][1], &Area[0][2]);
              get_face_area(p0, p3, p2, &Area[1][0], &Area[1][1], &Area[1][2]);
              get_face_area(p0, p1, p3, &Area[2][0], &Area[2][1], &Area[2][2]);
              double rBx, rBy, rBz;
              solve3x3(Area, phi, &rBx, &rBy, &rBz);

              double tetra_vol = calculate_tetra_volume(p0, p1, p2, p3);        // more geom. precision if calculate exact intersection vol, but not needed
              double weight = tetra_vol;
              if(gsl_finite(rBx) && gsl_finite(rBy) && gsl_finite(rBz) && gsl_finite(weight) && weight > 0)
                {
                  Bx_running += weight * rBx;
                  By_running += weight * rBy;
                  Bz_running += weight * rBz;
                  normalization += weight;
                }

              for(m = 0, ll = ii = jj = -1; m < 4; m++)
                {
                  if(next->p[m] == prev->p[k])
                    ll = m;
                  if(next->p[m] == prev->p[i])
                    ii = m;
                  if(next->p[m] == prev->p[j])
                    jj = m;
                }

              if(ll < 0 || ii < 0 || jj < 0)
                terminate("inconsistency");

              kk = 6 - (ll + ii + jj);

              prev = next;

              i = ii;
              l = ll;
              j = jj;
              k = kk;

            }
          while(next != t);

          if(edge == last_edge)
            break;

          edge = DC[edge].next;

        }                       //end of while loop
      SphP[pi].B[0] = Bx_running / normalization;
      SphP[pi].B[1] = By_running / normalization;
      SphP[pi].B[2] = Bz_running / normalization;
      // add the mean field
      // mean B-field is B_c  (B = B_c a^-2) for cosmological runs
      double bfac = 1. / (sqrt(All.UnitMass_in_g / All.UnitLength_in_cm) / (All.UnitTime_in_s / All.HubbleParam));
      SphP[pi].B[0] += All.CT_mean_Bx * bfac / sqrt(4. * M_PI);  /* convert Gauss-cgs to heavyside - lorentz */
      SphP[pi].B[1] += All.CT_mean_By * bfac / sqrt(4. * M_PI); 
      SphP[pi].B[2] += All.CT_mean_Bz * bfac / sqrt(4. * M_PI); 
      //update Bconserved too
      SphP[pi].BConserved[0] = SphP[pi].B[0] * SphP[pi].Volume;
      SphP[pi].BConserved[1] = SphP[pi].B[1] * SphP[pi].Volume;
      SphP[pi].BConserved[2] = SphP[pi].B[2] * SphP[pi].Volume;
    }

}
#else

/** get center B-fields from integral formulation of curl of A */
void correct_ctr_b(void) 
{  
  int idx, i;
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++) {
    i = TimeBinsHydro.ActiveParticleList[idx];
    if(i < 0)
      continue;  
      
    int k;
    for(k=0; k<3; k++) {
      SphP[i].A[k] = SphP[i].AConserved[k] / SphP[i].Volume;
    }  
    SphP[i].TimeLastBUpdate = All.Time;  
  }

  exchange_primitive_variables();  // to ensure we have all neighboring A
  
  // recover the B-fields from A: B*vol = sum area cross A

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++) {
    i = TimeBinsHydro.ActiveParticleList[idx];
    if(i < 0)
      continue;

    int edge = SphP[i].first_connection;
    int last_edge = SphP[i].last_connection;

    int vf = DC[edge].vf_index;
    int dp_index = Mesh.VF[vf].p1;
    if(dp_index == DC[edge].dp_index)
	    dp_index = Mesh.VF[vf].p2;

    point * dp_ref = &(Mesh.DP[dp_index]);

    double Bx_running = 0.;
    double By_running = 0.;
    double Bz_running = 0.; 
      
	  //loop over all the connections around a cell
	  while(1) {
	    int dq_index = DC[edge].dp_index;
      point *dq = &Mesh.DP[dq_index];
      int j = dq->index;
      if(j >= NumGas && dq->task == ThisTask)  j -= NumGas;
      
      vf = DC[edge].vf_index;
      double face_cx = Mesh.VF[vf].cx;  
      double face_cy = Mesh.VF[vf].cy;  
      double face_cz = Mesh.VF[vf].cz; 
      double face_area = Mesh.VF[vf].area;

      // extrapolate A to CM of face
      double Ap[3], Aq[3];
      Ap[0] = SphP[i].A[0];
      Ap[1] = SphP[i].A[1];
      Ap[2] = SphP[i].A[2];
      if(dq->task == ThisTask) {
        Aq[0] = SphP[j].A[0];
        Aq[1] = SphP[j].A[1];
        Aq[2] = SphP[j].A[2];
      } 
      else {
        Aq[0] = PrimExch[j].A[0];
        Aq[1] = PrimExch[j].A[1];
        Aq[2] = PrimExch[j].A[2];
      } 
    
      double faceAxArea, faceAyArea, faceAzArea;
      faceAxArea = 0.5*face_area*(Ap[0]+Aq[0]);
      faceAyArea = 0.5*face_area*(Ap[1]+Aq[1]);
      faceAzArea = 0.5*face_area*(Ap[2]+Aq[2]);
    
      double noutx, nouty, noutz, norm;
      noutx = dq->x - dp_ref->x;
      nouty = dq->y - dp_ref->y;
      noutz = dq->z - dp_ref->z;
      norm = sqrt(noutx*noutx + nouty*nouty + noutz*noutz);
      noutx /= norm;
      nouty /= norm;
      noutz /= norm;
      Bx_running += nouty*faceAzArea - faceAyArea*noutz;
      By_running += noutz*faceAxArea - faceAzArea*noutx;
      Bz_running += noutx*faceAyArea - faceAxArea*nouty;
      double coutx, couty, coutz, cnorm;
      coutx = face_cx - 0.5*(dq->x + dp_ref->x);
      couty = face_cy - 0.5*(dq->y + dp_ref->y);
      coutz = face_cz - 0.5*(dq->z + dp_ref->z);
      coutx /= norm;
      couty /= norm;
      coutz /= norm;
      double faceDAxArea, faceDAyArea, faceDAzArea;
      faceDAxArea = face_area*(Ap[0]-Aq[0]);
      faceDAyArea = face_area*(Ap[1]-Aq[1]);
      faceDAzArea = face_area*(Ap[2]-Aq[2]);
      Bx_running -= couty*faceDAzArea - faceDAyArea*coutz;
      By_running -= coutz*faceDAxArea - faceDAzArea*coutx;
      Bz_running -= coutx*faceDAyArea - faceDAxArea*couty;
    
      
	    if(edge == last_edge)
	      break;
      
	    edge = DC[edge].next;
	
    }			//end of while loop
    SphP[i].B[0] = Bx_running / SphP[i].Volume;
    SphP[i].B[1] = By_running / SphP[i].Volume;
    SphP[i].B[2] = Bz_running / SphP[i].Volume;
    SphP[i].BConserved[0] = SphP[i].B[0] * SphP[i].Volume;
    SphP[i].BConserved[1] = SphP[i].B[1] * SphP[i].Volume;
    SphP[i].BConserved[2] = SphP[i].B[2] * SphP[i].Volume;
    // add the mean field
    // mean B-field is B_c  (B = B_c a^-2) for cosmological runs
    double bfac = 1. / (sqrt(All.UnitMass_in_g / All.UnitLength_in_cm) / (All.UnitTime_in_s / All.HubbleParam));
    SphP[i].B[0] += All.CT_mean_Bx * bfac / sqrt(4. * M_PI);  /* convert Gauss-cgs to heavyside - lorentz */
    SphP[i].B[1] += All.CT_mean_By * bfac / sqrt(4. * M_PI); 
    SphP[i].B[2] += All.CT_mean_Bz * bfac / sqrt(4. * M_PI); 
    SphP[i].BConserved[0] = SphP[i].B[0] * SphP[i].Volume;
    SphP[i].BConserved[1] = SphP[i].B[1] * SphP[i].Volume;
    SphP[i].BConserved[2] = SphP[i].B[2] * SphP[i].Volume;
  }  

}
#endif













/* get face area - outward vector -- a helper function */
void get_face_area(point * dp0, point * dp1, point * dp2, double *Ax, double *Ay, double *Az)
{
  double vecA[3], vecB[3];
  vecA[0] = dp1->x - dp0->x;
  vecA[1] = dp1->y - dp0->y;
  vecA[2] = dp1->z - dp0->z;
  vecB[0] = dp2->x - dp0->x;
  vecB[1] = dp2->y - dp0->y;
  vecB[2] = dp2->z - dp0->z;
  *Ax = 0.5 * (vecA[1] * vecB[2] - vecA[2] * vecB[1]);
  *Ay = -0.5 * (vecA[0] * vecB[2] - vecA[2] * vecB[0]);
  *Az = 0.5 * (vecA[0] * vecB[1] - vecA[1] * vecB[0]);
}



/* get phi - outward vector -- a helper function */
double get_phi(point * dp0, point * dp1, point * dp2, point * dp_ref)
{
  double phi, A[3][3], p[3][3];

  int i0, i1, i2, k;
  i0 = dp0->index;
  i1 = dp1->index;
  i2 = dp2->index;
  if(i0 >= NumGas && dp0->task == ThisTask)
    i0 -= NumGas;
  if(i1 >= NumGas && dp1->task == ThisTask)
    i1 -= NumGas;
  if(i2 >= NumGas && dp2->task == ThisTask)
    i2 -= NumGas;

  // get the A's
  if(dp0->task == ThisTask)
    {
      for(k = 0; k < 3; k++)
        A[0][k] = SphP[i0].A[k];
    }
  else
    {
      for(k = 0; k < 3; k++)
        A[0][k] = PrimExch[i0].A[k];
    }

  if(dp1->task == ThisTask)
    {
      for(k = 0; k < 3; k++)
        A[1][k] = SphP[i1].A[k];
    }
  else
    {
      for(k = 0; k < 3; k++)
        A[1][k] = PrimExch[i1].A[k];
    }

  if(dp2->task == ThisTask)
    {
      for(k = 0; k < 3; k++)
        A[2][k] = SphP[i2].A[k];
    }
  else
    {
      for(k = 0; k < 3; k++)
        A[2][k] = PrimExch[i2].A[k];
    }

  //time_extrapolate_A_if_needed(A[0], dp0, dp_ref);
  //time_extrapolate_A_if_needed(A[1], dp1, dp_ref);
  //time_extrapolate_A_if_needed(A[2], dp2, dp_ref);

  // calculate phi
  phi = 0.5 * ((A[1][0] + A[0][0]) * (dp1->x - dp0->x) + (A[1][1] + A[0][1]) * (dp1->y - dp0->y) + (A[1][2] + A[0][2]) * (dp1->z - dp0->z));
  phi += 0.5 * ((A[2][0] + A[1][0]) * (dp2->x - dp1->x) + (A[2][1] + A[1][1]) * (dp2->y - dp1->y) + (A[2][2] + A[1][2]) * (dp2->z - dp1->z));
  phi += 0.5 * ((A[0][0] + A[2][0]) * (dp0->x - dp2->x) + (A[0][1] + A[2][1]) * (dp0->y - dp2->y) + (A[0][2] + A[2][2]) * (dp0->z - dp2->z));
  return phi;
}


/** Solve 3x3 matrix A * x = b equation. A = [A00, A01, A02 \\ A10 , ...]*/
void solve3x3(double A[3][3], double b[3], double *x, double *y, double *z)
{
  double det = A[0][2] * A[1][1] * A[2][0] - A[0][1] * A[1][2] * A[2][0] - A[0][2] * A[1][0] * A[2][1] + A[0][0] * A[1][2] * A[2][1] + A[0][1] * A[1][0] * A[2][2] - A[0][0] * A[1][1] * A[2][2];
  *x = ((A[1][2] * A[2][1] - A[1][1] * A[2][2]) * b[0] + (-A[0][2] * A[2][1] + A[0][1] * A[2][2]) * b[1] + (A[0][2] * A[1][1] - A[0][1] * A[1][2]) * b[2]) / det;
  *y = ((-A[1][2] * A[2][0] + A[1][0] * A[2][2]) * b[0] + (A[0][2] * A[2][0] - A[0][0] * A[2][2]) * b[1] + (-A[0][2] * A[1][0] + A[0][0] * A[1][2]) * b[2]) / det;
  *z = ((A[1][1] * A[2][0] - A[1][0] * A[2][1]) * b[0] + (-A[0][1] * A[2][0] + A[0][0] * A[2][1]) * b[1] + (A[0][1] * A[1][0] - A[0][0] * A[1][1]) * b[2]) / det;
}




/** Get fluxes due to mesh motion - upwind scheme */
void get_A_fluxes(int k, int q, int qother, double face_normal_vel, double face_velx, struct state *st_L, struct state *st_R, struct fluxes *flux)
{
  // face_velx is in face frame
  int d;
  double faceA[3];

  // get face states
  if(face_normal_vel < 0)
    {                           // face state is L
      faceA[0] = st_L->Ax;
      faceA[1] = st_L->Ay;
      faceA[2] = st_L->Az;
    }
  else
    {                           // face state is R
      faceA[0] = st_R->Ax;
      faceA[1] = st_R->Ay;
      faceA[2] = st_R->Az;
    }

  // moving mesh advection terms
  for(d = 0; d < 3; d++)
    {
      flux->A[d] = -face_normal_vel * faceA[d];
    }
}



void switch_to_CM(point ** p0, point * p1)
{
  *p1 = **p0;
  *p0 = p1;

  int i = p1->index;
  if(i >= NumGas && p1->task == ThisTask)
    i -= NumGas;

  double x0, y0, z0;
  x0 = p1->x;
  y0 = p1->y;
  z0 = p1->z;

  if(p1->task == ThisTask)
    {
      p1->x = SphP[i].Center[0];
      p1->y = SphP[i].Center[1];
      p1->z = SphP[i].Center[2];
    }
  else
    {
      p1->x = PrimExch[i].Center[0];
      p1->y = PrimExch[i].Center[1];
      p1->z = PrimExch[i].Center[2];
    }

  if(p1->x - x0 < -0.5 * boxSize_X)
    p1->x += boxSize_X;
  if(p1->x - x0 >= 0.5 * boxSize_X)
    p1->x -= boxSize_X;
  if(p1->y - y0 < -0.5 * boxSize_Y)
    p1->y += boxSize_Y;
  if(p1->y - y0 >= 0.5 * boxSize_Y)
    p1->y -= boxSize_Y;
  if(p1->z - z0 < -0.5 * boxSize_Z)
    p1->z += boxSize_Z;
  if(p1->z - z0 >= 0.5 * boxSize_Z)
    p1->z -= boxSize_Z;
}

/** Apply source terms to A using update_A() */
void do_mhd_ct_source_terms()
{
  int idx, i;
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      double hubble_a = 1.0;
      double atime = 1.0;
      if(All.ComovingIntegrationOn)
        {
          atime = All.Time;
          hubble_a = hubble_function(atime);
        }
      double dt_cell = 0.5 * (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval / hubble_a;     /* half the timestep of the cell */
      dt_cell /= atime;         // dA / dt = -(1/a) E_c

      update_A(i, dt_cell);
    }
}




void time_extrapolate_A_if_needed(double A[3], point * dp, point * dp_ref)
{
  int i, i_ref, k;
  double tlast, tlast_ref;
  double V[3], B[3], E[3];
  double hubble_a = 1.0;
  double atime = 1.0;
  if(All.ComovingIntegrationOn)
    {
      atime = All.Time;
      hubble_a = hubble_function(atime);
    }

  i = dp->index;
  if(i >= NumGas && dp->task == ThisTask)
    i -= NumGas;

  i_ref = dp_ref->index;
  if(i_ref >= NumGas && dp_ref->task == ThisTask)
    i_ref -= NumGas;

  if(dp->task == ThisTask)
    {
      tlast = SphP[i].TimeLastBUpdate;
    }
  else
    {
      tlast = PrimExch[i].TimeLastBUpdate;
    }

  if(dp_ref->task == ThisTask)
    {
      tlast_ref = SphP[i_ref].TimeLastBUpdate;
    }
  else
    {
      tlast_ref = PrimExch[i_ref].TimeLastBUpdate;
    }

  double dt_extrapolate = tlast_ref - tlast;

  if(dt_extrapolate > 0)
    {
      if(dp->task == ThisTask)
        {
          for(k = 0; k < 3; k++)
            {
              V[k] = P[i].Vel[k] - SphP[i].VelVertex[k];
              B[k] = SphP[i].B[k];
              V[k] /= atime;    // convert to peculiar velocity
            }
        }
      else
        {
          for(k = 0; k < 3; k++)
            {
              V[k] = PrimExch[i].VelGas[k] - PrimExch[i].VelVertex[k];
              B[k] = PrimExch[i].B[k];
              V[k] /= atime;    // convert to peculiar velocity
            }
        }
      // calculate E field: E = - V x B
      E[0] = -(V[1] * B[2] - V[2] * B[1]);
      E[1] = (V[0] * B[2] - V[2] * B[0]);
      E[2] = -(V[0] * B[1] - V[1] * B[0]);

      for(k = 0; k < 3; k++)
        A[k] += -dt_extrapolate * E[k];
    }

}



int not_both_active(point *dp_ref, point *dp0) {
  int i_ref, i0;

  i_ref = dp_ref->index;
  if(i_ref >= NumGas && dp_ref->task == ThisTask)
    i_ref -= NumGas;
    
  i0 = dp0->index;
  if(i0 >= NumGas && dp0->task == ThisTask)
    i0 -= NumGas;
    
  if((dp_ref->task == ThisTask) & (dp0->task == ThisTask)) {
    if(SphP[i_ref].TimeLastBUpdate != SphP[i0].TimeLastBUpdate)   return 1;
  }

  if((dp_ref->task != ThisTask) & (dp0->task == ThisTask)) {
    if(PrimExch[i_ref].TimeLastBUpdate != SphP[i0].TimeLastBUpdate)   return 1;
  }  
  
  if((dp_ref->task == ThisTask) & (dp0->task != ThisTask)) {
    if(SphP[i_ref].TimeLastBUpdate != PrimExch[i0].TimeLastBUpdate)   return 1;
  } 

  if((dp_ref->task != ThisTask) & (dp0->task != ThisTask)) {
    if(PrimExch[i_ref].TimeLastBUpdate != PrimExch[i0].TimeLastBUpdate)   return 1;
  } 
  
  return 0;
}




#endif /* #ifdef MHD_CTA */
#endif /* #ifdef MHD */
