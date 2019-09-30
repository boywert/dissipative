/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/constrained_transport/regularize_mesh_opt.c
 * \date        05/2016
 * \author      Philip Mocz
 * \brief       strongly centroidal lloyd mesh regularization and smooth Duffell regularization for set_vertex_velocities.c
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../allvars.h"
#include "../proto.h"
#include "../voronoi.h"



/* Duffell regularization */
#ifdef REGULARIZE_MESH_SMOOTH
void set_vertex_velocities_smooth(void)
{
  int idx, i, j, iter;
  
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;
      for(j = 0; j < 3; j++)
        SphP[i].VelVertex[j] = P[i].Vel[j];     /* make cell velocity equal to fluid's velocity */ // XXX
    }
    
  exchange_primitive_variables();   // exchange VelVertex   

  /* iterative algorithm - find avg vertex velocity of neighbors, and avg with self*/
  for(iter = 0; iter < 8; iter++)
    {
      mpi_printf("CM iteration: %d.\n", iter);
      point *DP = Mesh.DP;
      
      // calculate avg vertex velocity of neighbors
      for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
        {
          i = TimeBinsHydro.ActiveParticleList[idx];
          if(i < 0)
            continue;
            
          for(j = 0; j < 3; j++)
            SphP[i].VelVertexAvg[j] = 0.0;
            
          SphP[i].VelVertexAvgNorm = 0.0;
          
          if((P[i].Mass == 0) && (P[i].ID == 0))
            continue;           /* skip cells that have been swallowed or dissolved */
            
          int q = SphP[i].first_connection;
          while(q >= 0)
            {
              int dp = DC[q].dp_index;
              int particle = DP[dp].index;
              int vf = DC[q].vf_index;    
      
              if(particle < 0)
                {
                  q = DC[q].next;
                  continue;
                }
      
              if(particle >= NumGas && Mesh.DP[dp].task == ThisTask)
                particle -= NumGas;
      
              if(Mesh.DP[dp].task==ThisTask) {
                for(j = 0; j < 3; j++)
                  SphP[i].VelVertexAvg[j] += SphP[particle].VelVertex[j] * Mesh.VF[vf].area;
              } else  {
                for(j = 0; j < 3; j++)
                  SphP[i].VelVertexAvg[j] += PrimExch[particle].VelVertex[j] * Mesh.VF[vf].area;
              }  
              SphP[i].VelVertexAvgNorm += Mesh.VF[vf].area;
      
              if(q == SphP[i].last_connection)
                break;
      
              q = DC[q].next;
            }
            
          if(SphP[i].VelVertexAvgNorm==0)
            terminate("norm is 0"); 
               
        }

      // avg with self
      for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
        {
          i = TimeBinsHydro.ActiveParticleList[idx];
          if(i < 0)
            continue;
          for(j = 0; j < 3; j++)
            SphP[i].VelVertexAvg[j] /= SphP[i].VelVertexAvgNorm;
          for(j = 0; j < 3; j++)
            SphP[i].VelVertex[j] = 0.5*SphP[i].VelVertex[j] + 0.5*SphP[i].VelVertexAvg[j];
        }

      exchange_primitive_variables();   // exchange VelVertex     

    }

}

#endif




/* Quick and dirty temporary fix for Llyod reg since the original module was broken with the time stepping change - only works for single tstep, 1 processor*/
#ifdef REGULARIZE_MESH_LLOYD

void set_vertex_velocities_lloyd() {
  int i,j,idx, r;
  double dt_cell;
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      for(j = 0; j < 3; j++)
        SphP[i].VelVertex[j] = P[i].Vel[j];
        
      dt_cell = (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;  
     }
     
  if(dt_cell >0) 
  for(r = 0; r < 3; r++) {
    mpi_printf("MESH_REGULARIZATION LLOYD: step %d, %f\n", r, dt_cell);
    free_mesh();
    drift_all_particles();
    create_mesh();
    mesh_setup_exchange();
    
    for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
      {
        i = TimeBinsHydro.ActiveParticleList[idx];
        if(i < 0)
          continue;
        for(j = 0; j < 3; j++)
          SphP[i].VelVertex[j] = -SphP[i].VelVertex[j];
       }
    
    drift_all_particles();
    
    for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
      {
        i = TimeBinsHydro.ActiveParticleList[idx];
        if(i < 0)
          continue;
        for(j = 0; j < 3; j++)
          SphP[i].VelVertex[j] = (SphP[i].Center[j] - P[i].Pos[j]) / dt_cell;
       }
    
  }
  
  free_mesh();
  create_mesh();
  mesh_setup_exchange();
  exchange_primitive_variables_and_gradients();
}

#endif



/* This was BROKEN when the time stepping scheme changed*/
#if defined(REGULARIZE_MESH_LLOYD) && !defined(REGULARIZE_MESH_LLOYD) //_BROKEN
extern unsigned char *Edge_visited;
#ifdef TWODIMS
extern const int edge_start[3];
extern const int edge_end[3];
#else
extern const int edge_start[6];
extern const int edge_end[6];
extern const int edge_opposite[6];
extern const int edge_nexttetra[6];
#endif


void predict_DP(point * p_in, point * p_out, double dt)
{
  p_out->task = p_in->task;
  p_out->index = p_in->index;
  p_out->x = p_in->x;
  p_out->y = p_in->y;
  p_out->z = p_in->z;
  int p_idx = p_out->index;
  if(p_idx >= NumGas)
    p_idx -= NumGas;
  if(p_idx >= 0)
    {
      p_out->x += dt * SphP[p_idx].VelVertex[0];
      p_out->y += dt * SphP[p_idx].VelVertex[1];
      p_out->z += dt * SphP[p_idx].VelVertex[2];
    }
}

#ifdef TWODIMS
void vv_tricircumcenter(double a[2], double b[2], double c[2], double circumcenter[2])
{
  double xba, yba, xca, yca;
  double balength, calength;
  double denominator;
  double xcirca, ycirca;

  /* Use coordinates relative to point `a' of the triangle. */
  xba = b[0] - a[0];
  yba = b[1] - a[1];
  xca = c[0] - a[0];
  yca = c[1] - a[1];
  /* Squares of lengths of the edges incident to `a'. */
  balength = xba * xba + yba * yba;
  calength = xca * xca + yca * yca;

  /* Calculate the denominator of the formulae. */
  denominator = 0.5 / (xba * yca - yba * xca);

  /* Calculate offset (from `a') of circumcenter. */
  xcirca = (yca * balength - yba * calength) * denominator;
  ycirca = (xba * calength - xca * balength) * denominator;
  circumcenter[0] = xcirca + a[0];
  circumcenter[1] = ycirca + a[1];
}

void predict_DTC(tessellation * T, tetra * t_in, tetra_center * tc_out, double dt)
{
  point *DP = T->DP;
  point *dp1 = &DP[t_in->p[0]];
  point *dp2 = &DP[t_in->p[1]];
  point *dp3 = &DP[t_in->p[2]];
  point dp1_p, dp2_p, dp3_p;
  double a[2], b[2], c[2], circumcenter[2];

  predict_DP(dp1, &dp1_p, dt);
  predict_DP(dp2, &dp2_p, dt);
  predict_DP(dp3, &dp3_p, dt);
  a[0] = dp1_p.x;
  a[1] = dp1_p.y;
  b[0] = dp2_p.x;
  b[1] = dp2_p.y;
  c[0] = dp3_p.x;
  c[1] = dp3_p.y;

  vv_tricircumcenter(a, b, c, circumcenter);
  tc_out->cx = circumcenter[0];
  tc_out->cy = circumcenter[1];
}

void process_edge_faces_and_volumes_predict(tessellation * T, int tt, int nr)
{
  int side;
  for(side = 0; side < 2; ++side)
    {
      int i, j, qq;
      double dt_drift;
      double nx, ny;
      double sx, sy;
      double hx, hy;
      double dvol, h;

      tetra *DT = T->DT;
      point *DP = T->DP;

      tetra *t = &DT[tt];
      //tetra_center *DTC = T->DTC;

      i = edge_start[nr];
      j = edge_end[nr];

      point dpi_obj;
      point dpj_obj;
      point *dpi = &dpi_obj;
      point *dpj = &dpj_obj;

      if(side == 0)
        {
          /*  the actual time-step of particle */
          integertime ti_step = P[DP[t->p[i]].index].TimeBin ? (((integertime) 1) << P[DP[t->p[i]].index].TimeBin) : 0;
          if(All.ComovingIntegrationOn)
            dt_drift = get_drift_factor(All.Ti_Current, All.Ti_Current + ti_step);
          else
            dt_drift = ti_step * All.Timebase_interval;
        }
      else
        {
          /*  the actual time-step of particle */
          integertime ti_step = P[DP[t->p[j]].index].TimeBin ? (((integertime) 1) << P[DP[t->p[j]].index].TimeBin) : 0;
          if(All.ComovingIntegrationOn)
            dt_drift = get_drift_factor(All.Ti_Current, All.Ti_Current + ti_step);
          else
            dt_drift = ti_step * All.Timebase_interval;
        }

      predict_DP(&DP[t->p[i]], dpi, dt_drift);
      predict_DP(&DP[t->p[j]], dpj, dt_drift);

      qq = t->t[nr];
      tetra *q = &DT[qq];

      tetra_center DTC_tt, DTC_qq;
      predict_DTC(T, t, &DTC_tt, dt_drift);
      predict_DTC(T, q, &DTC_qq, dt_drift);

      if(side == 0)
        {
          Edge_visited[tt] |= (1 << nr);
          Edge_visited[qq] |= (1 << (t->s[nr]));
        }

      double f_cx, f_cy, f_area;

      f_cx = 0.5 * (DTC_tt.cx + DTC_qq.cx);
      f_cy = 0.5 * (DTC_tt.cy + DTC_qq.cy);

      nx = DTC_tt.cx - DTC_qq.cx;
      ny = DTC_tt.cy - DTC_qq.cy;

      f_area = sqrt(nx * nx + ny * ny);

      hx = 0.5 * (dpi->x - dpj->x);
      hy = 0.5 * (dpi->y - dpj->y);

      h = sqrt(hx * hx + hy * hy);
      dvol = 0.5 * f_area * h;

      if(side == 0)
        if(dpi->task == ThisTask && dpi->index >= 0 && dpi->index < NumGas)
          {
            if(TimeBinSynchronized[P[dpi->index].TimeBin])
              {
                SphP[dpi->index].Volume_predict += dvol;

                /* let's now compute the center-of-mass of the pyramid at the bottom top */
                sx = (2.0 / 3) * f_cx + (1.0 / 3) * dpi->x;
                sy = (2.0 / 3) * f_cy + (1.0 / 3) * dpi->y;

                SphP[dpi->index].Center_predict[0] += dvol * sx;
                SphP[dpi->index].Center_predict[1] += dvol * sy;
              }
          }

      if(side == 1)
        if(dpj->task == ThisTask && dpj->index >= 0 && dpj->index < NumGas)
          {
            if(TimeBinSynchronized[P[dpj->index].TimeBin])
              {
                SphP[dpj->index].Volume_predict += dvol;

                /* let's now compute the center-of-mass of the pyramid on top */
                sx = (2.0 / 3) * f_cx + (1.0 / 3) * dpj->x;
                sy = (2.0 / 3) * f_cy + (1.0 / 3) * dpj->y;

                SphP[dpj->index].Center_predict[0] += dvol * sx;
                SphP[dpj->index].Center_predict[1] += dvol * sy;
              }
          }
    }
}
#endif

#if !defined(TWODIMS) && !defined(ONEDIMS)
void vv_tetcircumcenter(double a[3], double b[3], double c[3], double d[3], double circumcenter[3])
{
  double xba, yba, zba, xca, yca, zca, xda, yda, zda;
  double balength, calength, dalength;
  double xcrosscd, ycrosscd, zcrosscd;
  double xcrossdb, ycrossdb, zcrossdb;
  double xcrossbc, ycrossbc, zcrossbc;
  double denominator;
  double xcirca, ycirca, zcirca;

  /* Use coordinates relative to point `a' of the tetrahedron. */
  xba = b[0] - a[0];
  yba = b[1] - a[1];
  zba = b[2] - a[2];
  xca = c[0] - a[0];
  yca = c[1] - a[1];
  zca = c[2] - a[2];
  xda = d[0] - a[0];
  yda = d[1] - a[1];
  zda = d[2] - a[2];
  /* Squares of lengths of the edges incident to `a'. */
  balength = xba * xba + yba * yba + zba * zba;
  calength = xca * xca + yca * yca + zca * zca;
  dalength = xda * xda + yda * yda + zda * zda;
  /* Cross products of these edges. */
  xcrosscd = yca * zda - yda * zca;
  ycrosscd = zca * xda - zda * xca;
  zcrosscd = xca * yda - xda * yca;
  xcrossdb = yda * zba - yba * zda;
  ycrossdb = zda * xba - zba * xda;
  zcrossdb = xda * yba - xba * yda;
  xcrossbc = yba * zca - yca * zba;
  ycrossbc = zba * xca - zca * xba;
  zcrossbc = xba * yca - xca * yba;

  /* Calculate the denominator of the formulae. */
  denominator = 0.5 / (xba * xcrosscd + yba * ycrosscd + zba * zcrosscd);

  /* Calculate offset (from `a') of circumcenter. */
  xcirca = (balength * xcrosscd + calength * xcrossdb + dalength * xcrossbc) * denominator;
  ycirca = (balength * ycrosscd + calength * ycrossdb + dalength * ycrossbc) * denominator;
  zcirca = (balength * zcrosscd + calength * zcrossdb + dalength * zcrossbc) * denominator;
  circumcenter[0] = xcirca + a[0];
  circumcenter[1] = ycirca + a[1];
  circumcenter[2] = zcirca + a[2];
}


void predict_DTC(tessellation * T, tetra * t_in, tetra_center * tc_out, double dt)
{
  point *DP = T->DP;
  point *dp1 = &DP[t_in->p[0]];
  point *dp2 = &DP[t_in->p[1]];
  point *dp3 = &DP[t_in->p[2]];
  point *dp4 = &DP[t_in->p[3]];
  point dp1_p, dp2_p, dp3_p, dp4_p;
  double a[3], b[3], c[3], d[3], circumcenter[3];

  predict_DP(dp1, &dp1_p, dt);
  predict_DP(dp2, &dp2_p, dt);
  predict_DP(dp3, &dp3_p, dt);
  predict_DP(dp4, &dp4_p, dt);
  a[0] = dp1_p.x;
  a[1] = dp1_p.y;
  a[2] = dp1_p.z;
  b[0] = dp2_p.x;
  b[1] = dp2_p.y;
  b[2] = dp2_p.z;
  c[0] = dp3_p.x;
  c[1] = dp3_p.y;
  c[2] = dp3_p.z;
  d[0] = dp4_p.x;
  d[1] = dp4_p.y;
  d[2] = dp4_p.z;

  vv_tetcircumcenter(a, b, c, d, circumcenter);
  tc_out->cx = circumcenter[0];
  tc_out->cy = circumcenter[1];
  tc_out->cz = circumcenter[2];
}


void process_edge_faces_and_volumes_predict(tessellation * T, int tt, int nr)
{
  int side;
  for(side = 0; side < 2; ++side)
    {
      int i, j, k, l, m, ii, jj, kk, ll, nn, count, nr_next, p1, p2;
      double dt_drift;
      tetra *prev, *next;
      tetra_center prevc_obj, nextc_obj;
      tetra_center *prevc = &prevc_obj;
      tetra_center *nextc = &nextc_obj;
      double ax, ay, az;
      double bx, by, bz;
      double cx, cy, cz;
      double nx, ny, nz;
      double sx, sy, sz;
      double hhx, hhy, hhz;
      double darea, dvol, h;
      double f_cx, f_cy, f_cz, f_area;

      tetra *DT = T->DT;
      point *DP = T->DP;
      tetra_center *DTC = T->DTC;

      tetra *t = &DT[tt];

      i = edge_start[nr];
      j = edge_end[nr];
      k = edge_opposite[nr];
      l = edge_nexttetra[nr];

      if(side == 1)
        Edge_visited[tt] |= (1 << nr);

      p1 = t->p[i];
      p2 = t->p[j];

      point dpi_obj;
      point dpj_obj;
      point *dpi = &dpi_obj;
      point *dpj = &dpj_obj;

      if(side == 0)
        {
          /*  the actual time-step of particle */
          integertime ti_step = P[DP[p1].index].TimeBin ? (((integertime) 1) << P[DP[p1].index].TimeBin) : 0;
          if(All.ComovingIntegrationOn)
            dt_drift = get_drift_factor(All.Ti_Current, All.Ti_Current + ti_step);
          else
            dt_drift = ti_step * All.Timebase_interval;
        }
      else
        {
          /*  the actual time-step of particle */
          integertime ti_step = P[DP[p2].index].TimeBin ? (((integertime) 1) << P[DP[p2].index].TimeBin) : 0;
          if(All.ComovingIntegrationOn)
            dt_drift = get_drift_factor(All.Ti_Current, All.Ti_Current + ti_step);
          else
            dt_drift = ti_step * All.Timebase_interval;
        }

      predict_DP(&DP[p1], dpi, dt_drift);
      predict_DP(&DP[p2], dpj, dt_drift);

      f_area = 0;
      f_cx = 0;
      f_cy = 0;
      f_cz = 0;

      hhx = 0.5 * (dpi->x - dpj->x);
      hhy = 0.5 * (dpi->y - dpj->y);
      hhz = 0.5 * (dpi->z - dpj->z);

      h = sqrt(hhx * hhx + hhy * hhy + hhz * hhz);

      tetra_center DTC_tt;
      predict_DTC(T, t, &DTC_tt, dt_drift);
      cx = DTC_tt.cx;
      cy = DTC_tt.cy;
      cz = DTC_tt.cz;

      count = 0;

      prev = t;
      predict_DTC(T, prev, prevc, dt_drift);
      do
        {
          nn = prev->t[l];
          next = &DT[nn];
          predict_DTC(T, next, nextc, dt_drift);

          if(prev != t && next != t)
            {
              ax = prevc->cx - cx;
              ay = prevc->cy - cy;
              az = prevc->cz - cz;

              bx = nextc->cx - cx;
              by = nextc->cy - cy;
              bz = nextc->cz - cz;

              nx = ay * bz - az * by;
              ny = az * bx - ax * bz;
              nz = ax * by - ay * bx;

              sx = nextc->cx + prevc->cx + cx;
              sy = nextc->cy + prevc->cy + cy;
              sz = nextc->cz + prevc->cz + cz;

              darea = 0.5 * sqrt(nx * nx + ny * ny + nz * nz);
              f_area += darea;

              darea *= (1.0 / 3);

              f_cx += darea * sx;
              f_cy += darea * sy;
              f_cz += darea * sz;
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


          /* need to determine the edge number to be able to flag it */
          if(side == 1)
            for(nr_next = 0; nr_next < 6; nr_next++)
              if((edge_start[nr_next] == ii && edge_end[nr_next] == jj) || (edge_start[nr_next] == jj && edge_end[nr_next] == ii))
                {
                  if((Edge_visited[nn] & (1 << nr_next)) && next != t)
                    terminate("inconsistency");

                  Edge_visited[nn] |= (1 << nr_next);
                  break;
                }

          prev = next;
          prevc->cx = nextc->cx;
          prevc->cy = nextc->cy;
          prevc->cz = nextc->cz;
          i = ii;
          l = ll;
          j = jj;
          k = kk;

          count++;

          if(count > 1000)
            terminate("count is too large");
        }
      while(next != t);

      i = edge_start[nr];
      j = edge_end[nr];

      if(f_area)
        {
          f_cx /= f_area;
          f_cy /= f_area;
          f_cz /= f_area;
        }

      dvol = (1.0 / 3) * f_area * h;

      if(side == 0)
        if(dpi->task == ThisTask && dpi->index >= 0 && dpi->index < NumGas)
          {
            if(TimeBinSynchronized[P[dpi->index].TimeBin])
              {
                SphP[dpi->index].Volume_predict += dvol;

                /* let's now compute the center-of-mass of the pyramid at the bottom top */
                sx = 0.75 * f_cx + 0.25 * dpi->x;
                sy = 0.75 * f_cy + 0.25 * dpi->y;
                sz = 0.75 * f_cz + 0.25 * dpi->z;

                SphP[dpi->index].Center_predict[0] += dvol * sx;
                SphP[dpi->index].Center_predict[1] += dvol * sy;
                SphP[dpi->index].Center_predict[2] += dvol * sz;

              }
          }

      if(side == 1)
        if(dpj->task == ThisTask && dpj->index >= 0 && dpj->index < NumGas)
          {
            if(TimeBinSynchronized[P[dpj->index].TimeBin])
              {
                SphP[dpj->index].Volume_predict += dvol;

                /* let's now compute the center-of-mass of the pyramid on top */
                sx = 0.75 * f_cx + 0.25 * dpj->x;
                sy = 0.75 * f_cy + 0.25 * dpj->y;
                sz = 0.75 * f_cz + 0.25 * dpj->z;

                SphP[dpj->index].Center_predict[0] += dvol * sx;
                SphP[dpj->index].Center_predict[1] += dvol * sy;
                SphP[dpj->index].Center_predict[2] += dvol * sz;

              }
          }
    }
}
#endif


#ifdef TWODIMS
void process_edge_faces_and_volumes_predict_exact(tessellation * T, double *vol, double *cmxv, double *cmyv, double *cmzv, int tt, int nr)
{
  int i, j, qq, p1, p2;
  tetra *q;
  double nx, ny;
  double sx, sy;
  double hx, hy;
  double dvol, h;

  tetra *DT = T->DT;
  point *DP = T->DP;

  tetra *t = &DT[tt];
  tetra_center *DTC = T->DTC;

  i = edge_start[nr];
  j = edge_end[nr];

  point *dpi = &DP[t->p[i]];
  point *dpj = &DP[t->p[j]];

  qq = t->t[nr];
  q = &DT[qq];

  Edge_visited[tt] |= (1 << nr);
  Edge_visited[qq] |= (1 << (t->s[nr]));

  p1 = t->p[i];
  p2 = t->p[j];

  double f_cx, f_cy, f_area;

  f_cx = 0.5 * (DTC[tt].cx + DTC[qq].cx);
  f_cy = 0.5 * (DTC[tt].cy + DTC[qq].cy);

  nx = DTC[tt].cx - DTC[qq].cx;
  ny = DTC[tt].cy - DTC[qq].cy;

  f_area = sqrt(nx * nx + ny * ny);

  hx = 0.5 * (dpi->x - dpj->x);
  hy = 0.5 * (dpi->y - dpj->y);

  h = sqrt(hx * hx + hy * hy);
  dvol = 0.5 * f_area * h;

  if(p1 >= 0 && p1 < DeRefMesh.Ndp)
    {
      vol[p1] += dvol;

      /* let's now compute the center-of-mass of the pyramid at the bottom top */
      sx = (2.0 / 3) * f_cx + (1.0 / 3) * dpi->x;
      sy = (2.0 / 3) * f_cy + (1.0 / 3) * dpi->y;

      cmxv[p1] += dvol * sx;
      cmyv[p1] += dvol * sy;
    }

  if(p2 >= 0 && p2 < DeRefMesh.Ndp)
    {
      vol[p2] += dvol;

      /* let's now compute the center-of-mass of the pyramid at the bottom top */
      sx = (2.0 / 3) * f_cx + (1.0 / 3) * dpj->x;
      sy = (2.0 / 3) * f_cy + (1.0 / 3) * dpj->y;

      cmxv[p2] += dvol * sx;
      cmyv[p2] += dvol * sy;
    }

}
#endif

#if !defined(TWODIMS) && !defined(ONEDIMS)
void process_edge_faces_and_volumes_predict_exact(tessellation * T, double *vol, double *cmxv, double *cmyv, double *cmzv, int tt, int nr)
{
  int i, j, k, l, m, ii, jj, kk, ll, nn, count, nr_next, p1, p2;
  tetra *prev, *next;
  tetra_center *prevc, *nextc;
  double ax, ay, az;
  double bx, by, bz;
  double cx, cy, cz;
  double nx, ny, nz;
  double sx, sy, sz;
  double hhx, hhy, hhz;
  double darea, dvol, h;
  double f_cx, f_cy, f_cz, f_area;

  tetra *DT = T->DT;
  point *DP = T->DP;
  tetra_center *DTC = T->DTC;

  tetra *t = &DT[tt];

  i = edge_start[nr];
  j = edge_end[nr];
  k = edge_opposite[nr];
  l = edge_nexttetra[nr];

  Edge_visited[tt] |= (1 << nr);

  p1 = t->p[i];
  p2 = t->p[j];

  f_area = 0;
  f_cx = 0;
  f_cy = 0;
  f_cz = 0;

  hhx = 0.5 * (DP[p1].x - DP[p2].x);
  hhy = 0.5 * (DP[p1].y - DP[p2].y);
  hhz = 0.5 * (DP[p1].z - DP[p2].z);

  h = sqrt(hhx * hhx + hhy * hhy + hhz * hhz);

  cx = DTC[tt].cx;
  cy = DTC[tt].cy;
  cz = DTC[tt].cz;

  count = 0;

  prev = t;
  prevc = &DTC[tt];
  do
    {
      nn = prev->t[l];
      next = &DT[nn];
      nextc = &DTC[nn];

      if(prev != t && next != t)
        {
          ax = prevc->cx - cx;
          ay = prevc->cy - cy;
          az = prevc->cz - cz;

          bx = nextc->cx - cx;
          by = nextc->cy - cy;
          bz = nextc->cz - cz;

          nx = ay * bz - az * by;
          ny = az * bx - ax * bz;
          nz = ax * by - ay * bx;

          sx = nextc->cx + prevc->cx + cx;
          sy = nextc->cy + prevc->cy + cy;
          sz = nextc->cz + prevc->cz + cz;

          darea = 0.5 * sqrt(nx * nx + ny * ny + nz * nz);
          f_area += darea;

          darea *= (1.0 / 3);

          f_cx += darea * sx;
          f_cy += darea * sy;
          f_cz += darea * sz;
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


      /* need to determine the edge number to be able to flag it */

      for(nr_next = 0; nr_next < 6; nr_next++)
        if((edge_start[nr_next] == ii && edge_end[nr_next] == jj) || (edge_start[nr_next] == jj && edge_end[nr_next] == ii))
          {
            if((Edge_visited[nn] & (1 << nr_next)) && next != t)
              terminate("inconsistency");

            Edge_visited[nn] |= (1 << nr_next);
            break;
          }

      prev = next;
      prevc = nextc;
      i = ii;
      l = ll;
      j = jj;
      k = kk;

      count++;

      if(count > 1000)
        terminate("count is too large");
    }
  while(next != t);

  i = edge_start[nr];
  j = edge_end[nr];

  if(f_area)
    {
      f_cx /= f_area;
      f_cy /= f_area;
      f_cz /= f_area;
    }

  dvol = (1.0 / 3) * f_area * h;

  if(p1 >= 0 && p1 < DeRefMesh.Ndp)
    {
      vol[p1] += dvol;

      /* let's now compute the center-of-mass of the pyramid at the bottom top */
      sx = 0.75 * f_cx + 0.25 * DP[p1].x;
      sy = 0.75 * f_cy + 0.25 * DP[p1].y;
      sz = 0.75 * f_cz + 0.25 * DP[p1].z;

      cmxv[p1] += dvol * sx;
      cmyv[p1] += dvol * sy;
      cmzv[p1] += dvol * sz;
    }

  if(p2 >= 0 && p2 < DeRefMesh.Ndp)
    {
      vol[p2] += dvol;

      /* let's now compute the center-of-mass of the pyramid at the bottom top */
      sx = 0.75 * f_cx + 0.25 * DP[p2].x;
      sy = 0.75 * f_cy + 0.25 * DP[p2].y;
      sz = 0.75 * f_cz + 0.25 * DP[p2].z;

      cmxv[p2] += dvol * sx;
      cmyv[p2] += dvol * sy;
      cmzv[p2] += dvol * sz;
    }
}
#endif

void compute_CM_predict_exact(double *vol, double *cmxv, double *cmyv, double *cmzv)
{
  int i, bit, nr;

  for(i = 0; i < DeRefMesh.Ndp; i++)
    {
      vol[i] = 0;
      cmxv[i] = 0;
      cmyv[i] = 0;
      cmzv[i] = 0;
    }

  Edge_visited = mymalloc_movable(&Edge_visited, "Edge_visited", DeRefMesh.Ndt * sizeof(unsigned char));

  for(i = 0; i < DeRefMesh.Ndt; i++)
    Edge_visited[i] = 0;

  for(i = 0; i < DeRefMesh.Ndt; i++)
    {
      if(DeRefMesh.DT[i].t[0] < 0)      /* deleted ? */
        continue;

      bit = 1;
      nr = 0;

      while(Edge_visited[i] != EDGE_ALL)
        {
          if((Edge_visited[i] & bit) == 0)
            process_edge_faces_and_volumes_predict_exact(&DeRefMesh, vol, cmxv, cmyv, cmzv, i, nr);

          bit <<= 1;
          nr++;
        }
    }

  myfree(Edge_visited);
}


void compute_CM_predict(void)
{
  int i, bit, nr;

  Edge_visited = mymalloc_movable(&Edge_visited, "Edge_visited", Mesh.Ndt * sizeof(unsigned char));

  for(i = 0; i < Mesh.Ndt; i++)
    Edge_visited[i] = 0;

  for(i = 0; i < Mesh.Ndt; i++)
    {
      if(Mesh.DT[i].t[0] < 0)   /* deleted ? */
        continue;

      bit = 1;
      nr = 0;

      while(Edge_visited[i] != EDGE_ALL)
        {
          if((Edge_visited[i] & bit) == 0)
            process_edge_faces_and_volumes_predict(&Mesh, i, nr);

          bit <<= 1;
          nr++;
        }
    }

  myfree(Edge_visited);
}

static int *ref_SphP_dp_index;


void set_vertex_velocities_lloyd(void)
{

  int idx, i, j, k, iter;
  double dt, d_anchor;

  /* set sph particle to delaunay particle map */
  ref_SphP_dp_index = mymalloc_movable(&ref_SphP_dp_index, "ref_SphP_dp_index", NumGas * sizeof(int));
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;
      ref_SphP_dp_index[i] = -1;
    }
  for(i = 0; i < Mesh.Ndp; i++)
    if(Mesh.DP[i].task == ThisTask)
      {
        int li = Mesh.DP[i].index;
        if(li >= 0 && li < NumGas)
          if(ref_SphP_dp_index[li] < 0)
            ref_SphP_dp_index[li] = i;  /* only guaranteed to be set for active cells */
      }

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      /*  the actual time-step of particle */
      integertime ti_step = P[i].TimeBin ? (((integertime) 1) << P[i].TimeBin) : 0;
      if(All.ComovingIntegrationOn)
        dt = get_drift_factor(All.Ti_Current, All.Ti_Current + ti_step);
      else
        dt = ti_step * All.Timebase_interval;

      for(j = 0; j < 3; j++)
        SphP[i].VelVertex[j] = P[i].Vel[j];     /* make cell velocity equal to fluid's velocity */
      /* now let's add the gradient of the pressure force */
      double acc[3];
      if(SphP[i].Density > 0)
        {
          acc[0] = -SphP[i].Grad.dpress[0] / SphP[i].Density;
          acc[1] = -SphP[i].Grad.dpress[1] / SphP[i].Density;
          acc[2] = -SphP[i].Grad.dpress[2] / SphP[i].Density;
          SphP[i].VelVertex[0] += 0.5 * dt * acc[0];
          SphP[i].VelVertex[1] += 0.5 * dt * acc[1];
          SphP[i].VelVertex[2] += 0.5 * dt * acc[2];
        }
      for(j = 0; j < 3; j++)
        SphP[i].Center_anchor[j] = P[i].Pos[j] + dt * SphP[i].VelVertex[j];
    }

  /* iterative Lloyd's algorithm */
  for(iter = 0; iter < 5; ++iter)
    {
      mpi_printf("CM iteration: %d.\n", iter);
      for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
        {
          i = TimeBinsHydro.ActiveParticleList[idx];
          if(i < 0)
            continue;
          SphP[i].Volume_predict = 0;
          for(j = 0; j < 3; j++)
            SphP[i].Center_predict[j] = 0;
        }

      compute_CM_predict();

      /* if cell connection changes too much, we need to compute the predicted CM exactly */
      for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
        {
          i = TimeBinsHydro.ActiveParticleList[idx];
          if(i < 0)
            continue;
          if((P[i].Mass == 0) && (P[i].ID == 0))
            continue;           /* skip cells that have been swallowed or dissolved */

          /* if cell connection changes too much, we need to compute the predicted CM exactly */
          if(fabs(SphP[i].Volume_predict - SphP[i].Volume) > 0.2 * SphP[i].Volume ||
             pow(4. * fabs(SphP[i].Center_predict[0] / SphP[i].Volume_predict - P[i].Pos[0]), NUMDIMS) > SphP[i].Volume ||
             pow(4. * fabs(SphP[i].Center_predict[1] / SphP[i].Volume_predict - P[i].Pos[1]), NUMDIMS) > SphP[i].Volume ||
             pow(4. * fabs(SphP[i].Center_predict[2] / SphP[i].Volume_predict - P[i].Pos[2]), NUMDIMS) > SphP[i].Volume)
            {
              printf("making a more precise CM prediction...\n");

              /*  the actual time-step of particle */
              integertime ti_step = P[i].TimeBin ? (((integertime) 1) << P[i].TimeBin) : 0;
              if(All.ComovingIntegrationOn)
                dt = get_drift_factor(All.Ti_Current, All.Ti_Current + ti_step);
              else
                dt = ti_step * All.Timebase_interval;
              SphP[i].Volume_predict = 0;
              for(j = 0; j < 3; j++)
                SphP[i].Center_predict[j] = 0;

              point *DP = Mesh.DP;
              point insert_pt;
              int q = SphP[i].first_connection;
            /**** create the voronoi cell of i as an auxiliary mesh */
              initialize_and_create_first_tetra(&DeRefMesh);
              DeRefMesh.DTC = mymalloc_movable(&DeRefMesh.DTC, "DeRefDTC", DeRefMesh.MaxNdt * sizeof(tetra_center));
              DeRefMesh.DTF = mymalloc_movable(&DeRefMesh.DTF, "DeRefDTF", DeRefMesh.MaxNdt * sizeof(char));
              for(k = 0; k < DeRefMesh.Ndt; k++)
                DeRefMesh.DTF[k] = 0;
              int tlast = 0;
              while(q >= 0)
                {
                  int dp = DC[q].dp_index;
                  int particle = DP[dp].index;

                  if(particle < 0)
                    {
                      q = DC[q].next;
                      continue;
                    }

                  if(particle >= NumGas)
                    particle -= NumGas;
                  //
                  if(DeRefMesh.Ndp + 2 >= DeRefMesh.MaxNdp)
                    {
                      DeRefMesh.Indi.AllocFacNdp *= ALLOC_INCREASE_FACTOR;
                      DeRefMesh.MaxNdp = DeRefMesh.Indi.AllocFacNdp;
                      DeRefMesh.DP -= 5;
                      DeRefMesh.DP = myrealloc_movable(DeRefMesh.DP, (DeRefMesh.MaxNdp + 5) * sizeof(point));
                      DeRefMesh.DP += 5;
                    }

                  insert_pt = Mesh.DP[dp];
                  insert_pt.x += dt * SphP[particle].VelVertex[0];
                  insert_pt.y += dt * SphP[particle].VelVertex[1];
                  insert_pt.z += dt * SphP[particle].VelVertex[2];
                  DeRefMesh.DP[DeRefMesh.Ndp] = insert_pt;

#ifndef OPTIMIZE_MEMORY_USAGE
                  set_integers_for_point(&DeRefMesh, DeRefMesh.Ndp);
#endif
                  tlast = insert_point(&DeRefMesh, DeRefMesh.Ndp, tlast);

                  DeRefMesh.Ndp++;
                  //

                  if(q == SphP[i].last_connection)
                    break;

                  q = DC[q].next;
                }

              /* now add also the point i itself  */
              insert_pt = Mesh.DP[ref_SphP_dp_index[i]];
              insert_pt.x += dt * SphP[i].VelVertex[0];
              insert_pt.y += dt * SphP[i].VelVertex[1];
              insert_pt.z += dt * SphP[i].VelVertex[2];
              DeRefMesh.DP[DeRefMesh.Ndp] = insert_pt;

#ifndef OPTIMIZE_MEMORY_USAGE
              set_integers_for_point(&DeRefMesh, DeRefMesh.Ndp);
#endif
              tlast = insert_point(&DeRefMesh, DeRefMesh.Ndp, tlast);
              DeRefMesh.Ndp++;

              /* compute circumcircles */
              compute_circumcircles(&DeRefMesh);

              double *Volume = mymalloc("Volume", DeRefMesh.Ndp * sizeof(double));
              double *CMxV = mymalloc("CMxV", DeRefMesh.Ndp * sizeof(double));
              double *CMyV = mymalloc("CMyV", DeRefMesh.Ndp * sizeof(double));
              double *CMzV = mymalloc("CMzV", DeRefMesh.Ndp * sizeof(double));

              /* compute CMs */
              compute_CM_predict_exact(Volume, CMxV, CMyV, CMzV);
              SphP[i].Volume_predict = Volume[DeRefMesh.Ndp - 1];
              SphP[i].Center_predict[0] = CMxV[DeRefMesh.Ndp - 1];
              SphP[i].Center_predict[1] = CMyV[DeRefMesh.Ndp - 1];
              SphP[i].Center_predict[2] = CMzV[DeRefMesh.Ndp - 1];
              myfree(CMzV);
              myfree(CMyV);
              myfree(CMxV);
              myfree(Volume);

              myfree(DeRefMesh.DTF);
              myfree(DeRefMesh.DTC);
              DeRefMesh.DTC = NULL;
              myfree(DeRefMesh.DT);
              myfree(DeRefMesh.DP - 5);
              myfree(DeRefMesh.VF);

            }
        }

      for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
        {
          i = TimeBinsHydro.ActiveParticleList[idx];
          if(i < 0)
            continue;
          /*  the actual time-step of particle */
          integertime ti_step = P[i].TimeBin ? (((integertime) 1) << P[i].TimeBin) : 0;
          if(All.ComovingIntegrationOn)
            dt = get_drift_factor(All.Ti_Current, All.Ti_Current + ti_step);
          else
            dt = ti_step * All.Timebase_interval;
          if((dt > 0) && (SphP[i].Volume_predict > 0))
            {
              for(j = 0; j < 3; j++)
                SphP[i].Center_predict[j] /= SphP[i].Volume_predict;
              d_anchor =
                sqrt(pow(SphP[i].Center_predict[0] - SphP[i].Center_anchor[0], 2) + pow(SphP[i].Center_predict[1] - SphP[i].Center_anchor[1], 2) +
                     pow(SphP[i].Center_predict[2] - SphP[i].Center_anchor[2], 2));
              /* limit displacement to a fraction of the sound speed motion XXX check cosmological units XX */
              for(j = 0; j < 3; j++)
                {
                  SphP[i].Center_predict[j] = SphP[i].Center_anchor[j] + (SphP[i].Center_predict[j] - SphP[i].Center_anchor[j]) * fmin(1.0, 0.05 * dt * All.cf_atime * get_sound_speed(i) / d_anchor);
                }
#ifdef REFINEMENT_MERGE_PAIRS
              /* do not do any grid regularization when we are actually trying to move two
                 mesh generating points very close together */
              if(SphP[i].DerefPartnerId <= 0)
                {
                  SphP[i].VelVertex[0] = nearest_x(SphP[i].Center_predict[0] - P[i].Pos[0]) / dt;
                  SphP[i].VelVertex[1] = nearest_y(SphP[i].Center_predict[1] - P[i].Pos[1]) / dt;
                  SphP[i].VelVertex[2] = nearest_z(SphP[i].Center_predict[2] - P[i].Pos[2]) / dt;
                }
#else
              SphP[i].VelVertex[0] = nearest_x(SphP[i].Center_predict[0] - P[i].Pos[0]) / dt;
              SphP[i].VelVertex[1] = nearest_y(SphP[i].Center_predict[1] - P[i].Pos[1]) / dt;
              SphP[i].VelVertex[2] = nearest_z(SphP[i].Center_predict[2] - P[i].Pos[2]) / dt;
#endif
            }
        }

      exchange_primitive_variables();   // exchange VelVertex     

    }
  myfree(ref_SphP_dp_index);
}

#endif
