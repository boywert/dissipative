/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/rt/rt_advect.c
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../allvars.h"
#include "../proto.h"
#include "../voronoi.h"

#ifdef RT_ADVECT

static struct flux_list_data
{
  int task, index, dirindex;
  double dPhotons;
}
 *FluxList;

int Nflux, MaxNflux;


#ifndef RT_HEALPIX_NSIDE

static struct radinflux_list_data
{
  int task, index, sourceid;
  MyFloat sourcepos[3];
  double dPhotons;
}
 *RadinFluxList;

struct diff_list
{
  int sourceid;
  double dPhotons;
  MyFloat sourcepos[3];
}
 *diff;

int Nradinflux, MaxNradinflux;
#endif

face *VF;
point *DP;

void rt_apply_radiation_flux(void);
int rt_check_responsibility_of_this_task(int p1, int p2);
void pix2vec_ring(long nside, long ipix, double *vec);
int rt_face_get_normals(int i);

void rt_advect_radiation(tessellation * T, double dt)
{
  CPU_Step[CPU_MISC] += measure_time();
  int i, j, iside, processflag;
  double nx, ny, nz, nn;
  double dx, dy, dz, fac;
  int ri, li, give;
  double densphot;
  double velPhot[3];
  struct rt_grad_data *g;
  double face_dt, delta_densphot;
  double dt_half, prefac = 1.0;
  double flux_phot;
  double sum, phot_before, phot_after, xtmp, ytmp, ztmp;
  int point;
  double dist;
#ifndef RT_HEALPIX_NSIDE
  MyFloat *sourcepos;
  int sourceid;
#endif

  set_cosmo_factors_for_current_time();

  VF = T->VF;
  DP = T->DP;

  for(i = 0, sum = 0; i < NumGas; i++)
    for(j = 0; j < RT_N_DIR; j++)
      {
        sum += SphP[i].Photons[j];
      }

  MPI_Allreduce(&sum, &phot_before, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


#ifndef RT_HEALPIX_NSIDE
  MaxNradinflux = T->Indi.AllocFacNradinflux;
  Nradinflux = 0;
  RadinFluxList = mymalloc_movable(&RadinFluxList, "RadinFluxList", MaxNradinflux * sizeof(struct radinflux_list_data));
#endif

  MaxNflux = T->Indi.AllocFacNflux;
  Nflux = 0;
  FluxList = mymalloc_movable(&FluxList, "FluxList", MaxNflux * sizeof(struct flux_list_data));

  for(i = 0; i < T->Nvf; i++)
    {
      face_dt = 0;              /* the default is that this face is not active */

      if(rt_face_get_normals(i))
        continue;

      if(rt_check_responsibility_of_this_task(VF[i].p1, VF[i].p2))
        continue;

      li = DP[VF[i].p1].index;
      ri = DP[VF[i].p2].index;

      if(li < 0 || ri < 0)
        continue;

      if(li >= NumGas && DP[VF[i].p1].task == ThisTask)
        li -= NumGas;

      if(ri >= NumGas && DP[VF[i].p2].task == ThisTask)
        ri -= NumGas;

      /* normal vector pointing to "right" state */
      nx = DP[VF[i].p2].x - DP[VF[i].p1].x;
      ny = DP[VF[i].p2].y - DP[VF[i].p1].y;
      nz = DP[VF[i].p2].z - DP[VF[i].p1].z;

      nn = sqrt(nx * nx + ny * ny + nz * nz);
      nx /= nn;
      ny /= nn;
      nz /= nn;


      /* calculate average direction vector from the two photon field gradients on both sides */
#ifdef RT_HEALPIX_NSIDE
      struct rt_grad_data *g1, *g2;

      if(DP[VF[i].p1].task == ThisTask)
        g1 = &SphP[li].rt_Grad;
      else
        g1 = &RTGradExch[li];

      if(DP[VF[i].p2].task == ThisTask)
        g2 = &SphP[ri].rt_Grad;
      else
        g2 = &RTGradExch[ri];

      double dir[3];

      dir[0] = -0.5 * (g1->ddensphot_unlimited[0] + g2->ddensphot_unlimited[0]);
      dir[1] = -0.5 * (g1->ddensphot_unlimited[1] + g2->ddensphot_unlimited[1]);
      dir[2] = -0.5 * (g1->ddensphot_unlimited[2] + g2->ddensphot_unlimited[2]);

      double len = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
#endif

      for(j = 0; j < RT_N_DIR; j++)
        {
          for(iside = 0; iside < 2; iside++)    /* consider both sides of the face in turn */
            {
              if(iside == 0)
                {
                  give = li;
                  point = VF[i].p1;
                }
              else
                {
                  give = ri;
                  point = VF[i].p2;
                }

              /* get the gradients and photon densitites of the used source side */
              if(DP[point].task == ThisTask)
                {
                  densphot = SphP[give].DensPhot[j];
                  g = &SphP[give].rt_Grad;
#ifndef RT_HEALPIX_NSIDE
                  sourcepos = &SphP[give].SourcePos[j][0];
                  sourceid = SphP[give].SourceID[j];
#endif
                }
              else
                {
                  densphot = RTPrimExch[give].DensPhot[j];
                  g = &RTGradExch[give];
#ifndef RT_HEALPIX_NSIDE
                  sourcepos = &RTPrimExch[give].SourcePos[j][0];
                  sourceid = RTPrimExch[give].SourceID[j];
#endif
                }


              /* calculate the photon velocity vector in the direction of propagation */
#ifdef RT_HEALPIX_NSIDE
              velPhot[0] = rt_vec[j][0];
              velPhot[1] = rt_vec[j][1];
              velPhot[2] = rt_vec[j][2];

              if(len > 0)       /* if there is a defined gradient direction, we take that, but projected into the allowed cone */
                {
                  double d[3], gg[3], m[3], n[3];

                  gg[0] = dir[0] / len;
                  gg[1] = dir[1] / len;
                  gg[2] = dir[2] / len;

                  d[0] = rt_vec[j][0];
                  d[1] = rt_vec[j][1];
                  d[2] = rt_vec[j][2];

                  n[0] = d[1] * gg[2] - d[2] * gg[1];
                  n[1] = d[2] * gg[0] - d[0] * gg[2];
                  n[2] = d[0] * gg[1] - d[1] * gg[0];

                  double alpha = acos(gg[0] * d[0] + gg[1] * d[1] + gg[2] * d[2]);
#ifdef TWODIMS
                  double alphamax = 0.5 * (2 * M_PI / RT_N_DIR);
#else
                  double alphamax = sqrt((4 * M_PI / RT_N_DIR) / M_PI);
#endif

                  if(alpha > alphamax)
                    {
                      m[0] = n[1] * d[2] - n[2] * d[1];
                      m[1] = n[2] * d[0] - n[0] * d[2];
                      m[2] = n[0] * d[1] - n[1] * d[0];

                      double mcoeff = sin(alphamax);
                      double dcoeff = cos(alphamax);

                      gg[0] = mcoeff * m[0] + dcoeff * d[0];
                      gg[1] = mcoeff * m[1] + dcoeff * d[1];
                      gg[2] = mcoeff * m[2] + dcoeff * d[2];
#ifdef TWODIMS
                      gg[2] = 0;
#endif
                    }

                  velPhot[0] = gg[0];
                  velPhot[1] = gg[1];
                  velPhot[2] = gg[2];
                }


#else
              velPhot[0] = VF[i].cx - sourcepos[0];
              velPhot[1] = VF[i].cy - sourcepos[1];
              velPhot[2] = VF[i].cz - sourcepos[2];
#endif

              prefac = 1.0;

#ifndef RT_HEALPIX_NSIDE
              if(j == 0)        /* this the diffuse component */
                {
                  if(iside == 0)
                    {
                      velPhot[0] = nx;
                      velPhot[1] = ny;
                      velPhot[2] = nz;
                    }
                  else
                    {
                      velPhot[0] = -nx;
                      velPhot[1] = -ny;
                      velPhot[2] = -nz;
                    }

                  prefac = 0.5 / NUMDIMS;
                  densphot *= prefac;
                }
#endif

              dist = sqrt(velPhot[0] * velPhot[0] + velPhot[1] * velPhot[1] + velPhot[2] * velPhot[2]);

              velPhot[0] *= (CLIGHT / All.UnitVelocity_in_cm_per_s) / dist;
              velPhot[1] *= (CLIGHT / All.UnitVelocity_in_cm_per_s) / dist;
              velPhot[2] *= (CLIGHT / All.UnitVelocity_in_cm_per_s) / dist;

#ifndef RT_HEALPIX_NSIDE
              fac = velPhot[0] * nx + velPhot[1] * ny + velPhot[2] * nz;
#else
              /* now do a more accurate integral of the discretized solid angle over the face normal */
#ifdef TWODIMS
              double arclength = 2 * M_PI / RT_N_DIR;
              if(iside == 1)
                {
                  nx = -nx;
                  ny = -ny;
                }

              double vx = ny * velPhot[0] - nx * velPhot[1];
              double vy = nx * velPhot[0] + ny * velPhot[1];
              double phi = atan2(vy, vx);
              double alpha = phi - 0.5 * arclength;
              double beta = phi + 0.5 * arclength;
              if(beta < 0 && alpha < 0)
                fac = 0;
              else
                {
                  if(beta > M_PI)
                    beta = M_PI;
                  if(alpha < 0 && beta > 0)
                    alpha = 0;
                  if(beta < 0)
                    terminate("a");
                  double p1 = velPhot[1] * nx - velPhot[0] * ny;
                  double p2 = velPhot[0] * nx + velPhot[1] * ny;

                  alpha -= phi;
                  beta -= phi;
                  fac = (p1 * (cos(beta) - cos(alpha)) + p2 * (sin(beta) - sin(alpha))) / arclength;
                }

              if(iside == 1)
                {
                  nx = -nx;
                  ny = -ny;
                  fac = -fac;
                }
#else
              fac = velPhot[0] * nx + velPhot[1] * ny + velPhot[2] * nz;        /* 3D case */
#endif
#endif
              /* we only need to consider a side if it is on the upwind side */
              if(iside == 0)
                {
                  if(fac <= 0)
                    continue;
                }
              else
                {
                  if(fac >= 0)
                    continue;
                }

              /* interpolation vector for the used state */
              if(DP[point].task == ThisTask)
                {
                  dx = NEAREST_X(VF[i].cx - SphP[give].Center[0]);
                  dy = NEAREST_Y(VF[i].cy - SphP[give].Center[1]);
                  dz = NEAREST_Z(VF[i].cz - SphP[give].Center[2]);
                }
              else
                {
                  dx = NEAREST_X(VF[i].cx - PrimExch[give].Center[0]);
                  dy = NEAREST_Y(VF[i].cy - PrimExch[give].Center[1]);
                  dz = NEAREST_Z(VF[i].cz - PrimExch[give].Center[2]);
                }

              /* set the actual time-step for the face */
              face_dt = dt / All.cf_hubble_a;

              /* compute the half-step prediction times */
              dt_half = 0.5 * face_dt;

              /* calculate the extrapolated used state, predicted forward in time by half a time-step */
              if(densphot > 0)
                {
                  delta_densphot = -dt_half / All.cf_atime * (velPhot[0] * g->ddensphot[j][0] + velPhot[1]
                                                              * g->ddensphot[j][1] + velPhot[2] * g->ddensphot[j][2]) + g->ddensphot[j][0] * dx + g->ddensphot[j][1] * dy + g->ddensphot[j][2] * dz;

                  delta_densphot *= prefac;
                }
              else
                delta_densphot = 0;

              if(densphot + delta_densphot < 0)
                delta_densphot = 0;

              densphot += delta_densphot;

              /* compute net flux with dot-product of outward normal and area of face */
              /* multiplication with area and time-step comes later */

              flux_phot = densphot * fac / All.cf_atime;

              /* now apply the flux to update the conserved states of the cells */
              if(face_dt > 0 && flux_phot != 0) /* selects active faces */
                {
                  int k, p, q;
                  double dir;
                  double fac = face_dt * VF[i].area;

                  for(k = 0; k < 2; k++)
                    {
                      if(k == 0)
                        {
                          q = VF[i].p1;
                          p = DP[q].index;
                          dir = -fac;
                        }
                      else
                        {
                          q = VF[i].p2;
                          p = DP[q].index;
                          dir = +fac;
                        }

                      if(DP[q].task == ThisTask)
                        {
                          if(DP[q].index >= NumGas)     /* this is a local ghost point */
                            p -= NumGas;

                          /* note: this will be executed if P[p] is a local point */

#ifdef RT_HEALPIX_NSIDE
                          processflag = 1;
#else
                          if(dir * flux_phot < 0 || j == 0)
                            processflag = 1;    // outgoing
                          else
                            processflag = 0;    // incoming
#endif

                          if(processflag == 1)
                            SphP[p].Photons[j] += dir * flux_phot;

#ifndef RT_HEALPIX_NSIDE
                          if(processflag == 0)
                            {
                              if(Nradinflux >= MaxNradinflux)
                                {
                                  T->Indi.AllocFacNradinflux *= ALLOC_INCREASE_FACTOR;
                                  MaxNradinflux = T->Indi.AllocFacNradinflux;
                                  RadinFluxList = myrealloc_movable(RadinFluxList, MaxNradinflux * sizeof(struct radinflux_list_data));
                                  if(Nradinflux >= MaxNradinflux)
                                    terminate("Nradinflux >= MaxNradinflux");
                                }

                              RadinFluxList[Nradinflux].task = ThisTask;
                              RadinFluxList[Nradinflux].index = p;
                              RadinFluxList[Nradinflux].dPhotons = dir * flux_phot;
                              RadinFluxList[Nradinflux].sourceid = sourceid;
                              memcpy(RadinFluxList[Nradinflux].sourcepos, sourcepos, 3 * sizeof(MyFloat));
                              Nradinflux++;
                            }
#endif
                        }
                      else
                        {
                          /* here we have a foreign ghost point */
                          if(DP[q].originalindex < 0)
                            terminate("should not happen");
#ifdef RT_HEALPIX_NSIDE
                          processflag = 1;
#else
                          if(dir * flux_phot < 0 || j == 0)
                            processflag = 1;
                          else
                            processflag = 0;
#endif

                          if(processflag == 1)
                            {
                              if(Nflux >= MaxNflux)
                                {
                                  T->Indi.AllocFacNflux *= ALLOC_INCREASE_FACTOR;
                                  MaxNflux = T->Indi.AllocFacNflux;
                                  FluxList = myrealloc_movable(FluxList, MaxNflux * sizeof(struct flux_list_data));
                                  if(Nflux >= MaxNflux)
                                    terminate("Nflux >= MaxNflux");
                                }

                              FluxList[Nflux].task = DP[q].task;
                              FluxList[Nflux].index = DP[q].originalindex;
                              FluxList[Nflux].dirindex = j;
                              FluxList[Nflux].dPhotons = dir * flux_phot;
                              Nflux++;
                            }
#ifndef RT_HEALPIX_NSIDE
                          if(processflag == 0)
                            {
                              if(Nradinflux >= MaxNradinflux)
                                {
                                  T->Indi.AllocFacNradinflux *= ALLOC_INCREASE_FACTOR;
                                  MaxNradinflux = T->Indi.AllocFacNradinflux;
                                  RadinFluxList = myrealloc_movable(RadinFluxList, MaxNradinflux * sizeof(struct radinflux_list_data));
                                  if(Nradinflux >= MaxNradinflux)
                                    terminate("Nradinflux >= MaxNradinflux");
                                }

                              RadinFluxList[Nradinflux].task = DP[q].task;
                              RadinFluxList[Nradinflux].index = DP[q].originalindex;
                              RadinFluxList[Nradinflux].dPhotons = dir * flux_phot;
                              RadinFluxList[Nradinflux].sourceid = sourceid;
                              memcpy(RadinFluxList[Nradinflux].sourcepos, sourcepos, 3 * sizeof(MyFloat));
                              Nradinflux++;
                            }
#endif
                        }
                    }
                }
            }                   /* end loop over the two sides */
        }                       /* end loop over RT_N_DIR */
    }                           /* end of big loop over all faces */

  /* now exchange the flux-list and apply it when needed */

  rt_apply_radiation_flux();
  myfree(FluxList);

#ifndef RT_HEALPIX_NSIDE
  rt_select_new_brightest_sources();
  myfree(RadinFluxList);
#endif

  for(i = 0, sum = 0; i < NumGas; i++)
    for(j = 0; j < RT_N_DIR; j++)
      {
        SphP[i].DensPhot[j] = SphP[i].Photons[j] / SphP[i].Volume;
        sum += SphP[i].Photons[j];
      }

  MPI_Allreduce(&sum, &phot_after, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if(ThisTask == 0)
    printf("Photons before and after transfer: %g %g  diff=%g\n", phot_before, phot_after, phot_after - phot_before);

  if(fabs(phot_after - phot_before) > 0.05 * phot_before)
    terminate("fabs(phot_after - phot_before) > 0.05 * phot_before\n");

  for(i = 0; i < NumGas; i++)
    for(j = 0; j < RT_N_DIR; j++)
      {
#ifdef RT_ALLOW_ABSORBING_CELLS
        if(P[i].ID >= 1000000000)
          {
            SphP[i].Photons[j] = 0;
            SphP[i].DensPhot[j] = 0;
          }
#endif
        if(SphP[i].Photons[j] < -1.0e-6)
          {
            printf("neg phot %g ID=%d\n", SphP[i].Photons[j], P[i].ID);
            terminate("neg phot\n");
            SphP[i].Photons[j] = 0;
            SphP[i].DensPhot[j] = 0;
          }
      }

  CPU_Step[CPU_FLUXES] += measure_time();
}


#ifndef RT_HEALPIX_NSIDE




void rt_select_new_brightest_sources(void)
{
#define MAX_COUNT_DIFF 200

  int i, j, k, count, ngrp, nimport, pindex;

  /* first, let's bring all influxes to the correct target task */

  mysort(RadinFluxList, Nradinflux, sizeof(struct radinflux_list_data), radinflux_compare_task);

  for(j = 0; j < NTask; j++)
    Send_count[j] = 0;

  for(i = 0; i < Nradinflux; i++)
    Send_count[RadinFluxList[i].task]++;

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  struct radinflux_list_data *RadinFluxListGet = (struct radinflux_list_data *) mymalloc("RadinFluxListGet", nimport * sizeof(struct radinflux_list_data));

  /* exchange particle data */
  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&RadinFluxList[Send_offset[recvTask]],
                           Send_count[recvTask] * sizeof(struct radinflux_list_data), MPI_BYTE,
                           recvTask, TAG_DENS_A,
                           &RadinFluxListGet[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct radinflux_list_data), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  /* now let's sort the incoming sources according to index, and then sourceid */
  mysort(RadinFluxListGet, nimport, sizeof(struct radinflux_list_data), radinflux_compare_index_sourceid);

  /* we should now go over each of the particle indices, and first construct the cumulative flux into the cell */
  /* among the different ones we can then select the brightest ones */
  diff = mymalloc("diff", MAX_COUNT_DIFF * sizeof(struct diff_list));

  for(i = 0; i < nimport;)
    {
      pindex = RadinFluxListGet[i].index;       /* particle index */

      count = 0;
      diff[count].sourceid = RadinFluxListGet[i].sourceid;
      memcpy(diff[count].sourcepos, RadinFluxListGet[i].sourcepos, 3 * sizeof(MyFloat));
      diff[count].dPhotons = 0;
      count++;

      for(; i < nimport; i++)
        {
          if(RadinFluxListGet[i].index != pindex)
            break;

          if(RadinFluxListGet[i].sourceid == diff[count - 1].sourceid)
            diff[count - 1].dPhotons += RadinFluxListGet[i].dPhotons;
          else
            {
              if(count >= MAX_COUNT_DIFF)
                terminate("count >= MAX_COUNT_DIFF");

              diff[count].sourceid = RadinFluxListGet[i].sourceid;
              diff[count].dPhotons = RadinFluxListGet[i].dPhotons;

              memcpy(diff[count].sourcepos, RadinFluxListGet[i].sourcepos, 3 * sizeof(MyFloat));
              count++;
            }
        }

      /* let's now sort the different ones according to total input flux */
      mysort(diff, count, sizeof(struct diff_list), rt_compare_difflist_dphotons);

      if(count > RT_N_DIR - 1)
        {
          /* the weaker ones we assign to the diffuse flux */
          for(j = RT_N_DIR - 1; j < count; j++)
            SphP[pindex].Photons[0] += diff[j].dPhotons;

          count = RT_N_DIR - 1;
        }

      /* append the present sources that have photons left in the cell */

      for(j = 1; j < RT_N_DIR; j++)
        if(SphP[pindex].Photons[j] > 0)
          {
            if(count >= MAX_COUNT_DIFF)
              terminate("count >= MAX_COUNT_DIFF");

            diff[count].sourceid = SphP[pindex].SourceID[j];
            diff[count].dPhotons = SphP[pindex].Photons[j];
            memcpy(diff[count].sourcepos, &SphP[pindex].SourcePos[j][0], 3 * sizeof(MyFloat));
            count++;
          }

      /* let's now sort the different ones according to source-id, and then sum up equal source ids */
      mysort(diff, count, sizeof(struct diff_list), rt_compare_difflist_sourceid);

      for(k = 0, j = 0; k < count;)
        {
          diff[j] = diff[k];
          k++;
          j++;

          for(; k < count; k++)
            {
              if(diff[j - 1].sourceid == diff[k].sourceid)
                diff[j - 1].dPhotons += diff[k].dPhotons;
              else
                {
                  diff[j] = diff[k];
                  j++;
                }
            }
        }
      count = j;                /* new count of different ones */

      /* let's now sort again the different ones according to total photon number */

      mysort(diff, count, sizeof(struct diff_list), rt_compare_difflist_dphotons);

      if(count > RT_N_DIR - 1)
        {
          /* the weaker ones we assign to the diffuse flux */
          for(j = RT_N_DIR - 1; j < count; j++)
            SphP[pindex].Photons[0] += diff[j].dPhotons;

          count = RT_N_DIR - 1;
        }

      /* let's sort the remaining most luminous ones according to sourceid, which is adopted in the gradient estimation */
      mysort(diff, count, sizeof(struct diff_list), rt_compare_difflist_sourceid);

      for(j = 0; j < RT_N_DIR - 1; j++)
        {
          if(j < count && diff[j].dPhotons > 0)
            {
              SphP[pindex].Photons[1 + j] = diff[j].dPhotons;
              SphP[pindex].SourceID[1 + j] = diff[j].sourceid;
              memcpy(&SphP[pindex].SourcePos[1 + j][0], diff[j].sourcepos, 3 * sizeof(MyFloat));
            }
          else
            SphP[pindex].Photons[1 + j] = 0;
        }
    }

  myfree(diff);
  myfree(RadinFluxListGet);
}

#endif

double rt_get_advect_tistep(void)
{
  double rad, radmin, globradmin;
  int idx, i;

  set_cosmo_factors_for_current_time();

  radmin = MAX_REAL_NUMBER;
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      rad = get_cell_radius(i);
      if(rad < radmin)
        radmin = rad;
    }

  MPI_Allreduce(&radmin, &globradmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  double v = (CLIGHT / All.UnitVelocity_in_cm_per_s) / All.cf_atime;
  double dt = All.CourantFac * globradmin / v;

  return dt * All.cf_hubble_a;
}



int radflux_list_data_compare(const void *a, const void *b)
{
  if(((struct flux_list_data *) a)->task < (((struct flux_list_data *) b)->task))
    return -1;

  if(((struct flux_list_data *) a)->task > (((struct flux_list_data *) b)->task))
    return +1;

  return 0;
}

void rt_apply_radiation_flux(void)
{
  int i, j, p, nimport, ngrp, sendTask, recvTask;

  /* now exchange the flux-list and apply it when needed */

  mysort(FluxList, Nflux, sizeof(struct flux_list_data), radflux_list_data_compare);

  for(j = 0; j < NTask; j++)
    Send_count[j] = 0;

  for(i = 0; i < Nflux; i++)
    Send_count[FluxList[i].task]++;

  if(Send_count[ThisTask] > 0)
    terminate("Send_count[ThisTask]");

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  struct flux_list_data *FluxListGet = (struct flux_list_data *) mymalloc("FluxListGet", nimport * sizeof(struct flux_list_data));

  /* exchange particle data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&FluxList[Send_offset[recvTask]],
                           Send_count[recvTask] * sizeof(struct flux_list_data), MPI_BYTE,
                           recvTask, TAG_DENS_A,
                           &FluxListGet[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct flux_list_data), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  /* apply the fluxes */

  for(i = 0; i < nimport; i++)
    {
      p = FluxListGet[i].index;
      j = FluxListGet[i].dirindex;

      SphP[p].Photons[j] += FluxListGet[i].dPhotons;
    }


  myfree(FluxListGet);
}

int rt_check_responsibility_of_this_task(int p1, int p2)
{
  int low_p, high_p;

  if(DP[p1].ID < DP[p2].ID)
    {
      low_p = p1;
      high_p = p2;
    }
  else
    {
      low_p = p2;
      high_p = p1;
    }

  /* we need to check whether the one with the lower ID is a local particle */
  if(DP[low_p].task == ThisTask && DP[low_p].index < NumGas)
    return 0;

  return -1;                    /* we can skip this face on the local task */
}

int rt_face_get_normals(int i)
{
  int li, ri;
  double surface, surface_l, surface_r;
  int present_left, present_right;

  li = DP[VF[i].p1].index;
  ri = DP[VF[i].p2].index;

  if(li < 0 || ri < 0)
    return -1;

  if(li >= NumGas && DP[VF[i].p1].task == ThisTask)
    li -= NumGas;

  if(ri >= NumGas && DP[VF[i].p2].task == ThisTask)
    ri -= NumGas;

  if(DP[VF[i].p1].task == ThisTask)
    surface_l = SphP[li].SurfaceArea;
  else
    surface_l = PrimExch[li].SurfaceArea;

  if(DP[VF[i].p2].task == ThisTask)
    surface_r = SphP[ri].SurfaceArea;
  else
    surface_r = PrimExch[ri].SurfaceArea;

  if(surface_r > surface_l)
    surface = 1.0e-5 * surface_r;
  else
    surface = 1.0e-5 * surface_l;

  present_left = present_right = 0;

  /* if the area of this face is negligible compared to the surface
     of the larger cell, skip it */
  if(DP[VF[i].p1].task == ThisTask && DP[VF[i].p1].index < NumGas)
    if(VF[i].area > surface)
      present_left = 1;

  if(DP[VF[i].p2].task == ThisTask && DP[VF[i].p2].index < NumGas)
    if(VF[i].area > surface)
      present_right = 1;

  if(present_left == 0 && present_right == 0)
    {
      VF[i].area = 0;
      return -1;
    }

  return 0;
}

/** \bug This routine still uses the incorrect cell location routine based
    on Delaunay tetrahedrons. */
int rt_get_cell(tessellation * T, double x, double y, double z)
{
  double dist2max, dist2;
  int Ndp = T->Ndp;
  int MaxNdp = T->MaxNdp;
  int pp = Ndp, ret, moves, tt, i, ind, q, li;
  int ttrow = 0;

  point *DP = T->DP;
  tetra *DT = T->DT;
  point *p = &DP[pp];

  while((Ndp >= (MaxNdp - 5)))
    {
      T->Indi.AllocFacNdp *= ALLOC_INCREASE_FACTOR;
      MaxNdp = T->Indi.AllocFacNdp;
#ifdef VERBOSE
      printf("Task=%d: increase memory allocation, MaxNdp=%d Indi.AllocFacNdp=%g\n", ThisTask, MaxNdp, T->Indi.AllocFacNdp);
#endif
      DP -= 5;
      DP = myrealloc_movable(DP, (MaxNdp + 5) * sizeof(point));
      DP += 5;

      if(Ndp >= (MaxNdp - 5) && NumGas == 0)
        terminate("(Ndp >= (MaxNdp - 5)");
    }

#ifndef TWODIMS
  while(DT[ttrow].t[0] < 0 || DT[ttrow].p[0] == DPinfinity || DT[ttrow].p[1] == DPinfinity || DT[ttrow].p[2] == DPinfinity || DT[ttrow].p[3] == DPinfinity)
    ttrow++;
#else
  while(DT[ttrow].t[0] < 0 || DT[ttrow].p[0] == DPinfinity || DT[ttrow].p[1] == DPinfinity || DT[ttrow].p[2] == DPinfinity)
    ttrow++;
#endif

  p->x = x;
  p->y = y;
  p->z = z;

  set_integers_for_point(T, pp);

  /* get the tetrahedron that this point contains */
#ifdef TWODIMS
  tt = get_triangle(T, pp, &moves, &ret, ttrow);
#else
  tt = get_tetra(T, p, &moves, ttrow, &ret, &ret);
#endif
  tetra *t = &DT[tt];




  for(i = 0, ind = -1, dist2max = 1.0e30; i < DIMS + 1; i++)
    {
      double dx = x - DP[t->p[i]].x;
      double dy = y - DP[t->p[i]].y;
#ifndef TWODIMS
      double dz = z - DP[t->p[i]].z;
      dist2 = dx * dx + dy * dy + dz * dz;
#else
      dist2 = dx * dx + dy * dy;
#endif
      if(dist2 < dist2max)
        {
          dist2max = dist2;
          ind = i;
        }
    }

  q = -1;

  if(ind >= 0)
    if(DP[t->p[ind]].task == ThisTask)
      {
        li = DP[t->p[ind]].index;

        if(li >= 0 && li < NumGas)
          q = li;
      }

  return q;
}


void rt_init_sourceid()
{
#ifndef RT_HEALPIX_NSIDE
  int i, k;

  for(i = 0; i < NumGas; i++)
    for(k = 1; k < RT_N_DIR; k++)
      SphP[i].SourceID[k] = 1000;
#endif
}

#ifndef RT_HEALPIX_NSIDE
int radinflux_compare_index_sourceid(const void *a, const void *b)
{
  if(((struct radinflux_list_data *) a)->index < (((struct radinflux_list_data *) b)->index))
    return -1;

  if(((struct radinflux_list_data *) a)->index > (((struct radinflux_list_data *) b)->index))
    return +1;

  if(((struct radinflux_list_data *) a)->sourceid < (((struct radinflux_list_data *) b)->sourceid))
    return -1;

  if(((struct radinflux_list_data *) a)->sourceid > (((struct radinflux_list_data *) b)->sourceid))
    return +1;

  return 0;
}

int radinflux_compare_task(const void *a, const void *b)
{
  if(((struct radinflux_list_data *) a)->task < (((struct radinflux_list_data *) b)->task))
    return -1;

  if(((struct radinflux_list_data *) a)->task > (((struct radinflux_list_data *) b)->task))
    return +1;

  return 0;
}


int rt_compare_difflist_dphotons(const void *a, const void *b)
{
  if(((struct diff_list *) a)->dPhotons < (((struct diff_list *) b)->dPhotons))
    return +1;

  if(((struct diff_list *) a)->dPhotons > (((struct diff_list *) b)->dPhotons))
    return -1;

  return 0;
}


int rt_compare_difflist_sourceid(const void *a, const void *b)
{
  if(((struct diff_list *) a)->sourceid < (((struct diff_list *) b)->sourceid))
    return -1;

  if(((struct diff_list *) a)->sourceid > (((struct diff_list *) b)->sourceid))
    return +1;

  return 0;
}


#endif


#ifdef RT_HEALPIX_NSIDE
void rt_get_vectors(void)
{
  int i;

#ifdef TWODIMS
  for(i = 0; i < RT_N_DIR; i++)
    {
      rt_vec[i][0] = cos((2 * M_PI / RT_N_DIR) * (i));
      rt_vec[i][1] = sin((2 * M_PI / RT_N_DIR) * (i));
      rt_vec[i][2] = 0;

      mpi_printf("%d | %g %g %g\n", i, rt_vec[i][0], rt_vec[i][1], rt_vec[i][2]);
    }
#else
  double vec[3];

  for(i = 0; i < RT_N_DIR; i++)
    {
      pix2vec_ring(RT_HEALPIX_NSIDE, i, &vec[0]);

      rt_vec[i][0] = vec[0];
      rt_vec[i][1] = vec[1];
      rt_vec[i][2] = vec[2];

      mpi_printf("%d | %g %g %g\n", i, rt_vec[i][0], rt_vec[i][1], rt_vec[i][2]);
    }
#endif
}
#endif

#endif //RT_ADVECT
