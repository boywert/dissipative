/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/voronoi_proj.c
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
#include <stddef.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"
#include "voronoi.h"
#include "voronoi_proj.h"

#ifdef OPACITIES
#include "opacities/opacities_combined.h"
opacitydata opdata;
#endif

#ifdef VORONOI_PROJ
#ifndef VORONOI_DYNAMIC_UPDATE
#error VORONOI_PROJ requires VORONOI_DYNAMIC_UPDATE
#endif
#endif

#ifdef VORONOI_PROJ_TAU
#include "opal_eos.h"
#endif

#ifdef VORONOI_PROJ
static void voronoi_proj_write_file(char *name, struct projection_data *pdata, double *data);
static void voronoi_proj_exchange_rays(struct projection_data *pdata);
static void voronoi_proj_update_column_integral(int cell, double len, int pix, double delta_center[3], struct projection_data *pdata);
static int voronoi_proj_advance_rays_for_one_cell(struct projection_data *pdata);

int voronoi_proj_setup(int argc, char **argv, struct projection_data *pdata)
{
  int k;

  if(argc < 15)
    {
      printf("Not enough parameters, aborting.\n");
      return 0;
    }

  pdata->nx = atoi(argv[4]);
  pdata->ny = atoi(argv[5]);

  pdata->center[0] = atof(argv[6]);
  pdata->center[1] = atof(argv[7]);
  pdata->center[2] = atof(argv[8]);

  pdata->dir[0] = atof(argv[9]);
  pdata->dir[1] = atof(argv[10]);
  pdata->dir[2] = atof(argv[11]);

  pdata->boxx = atof(argv[12]);
  pdata->boxy = atof(argv[13]);
  pdata->boxz = atof(argv[14]);

  pdata->pflag = 0;

  double dd = sqrt(pdata->dir[0] * pdata->dir[0] + pdata->dir[1] * pdata->dir[1] + pdata->dir[2] * pdata->dir[2]);
  for(k = 0; k < 3; k++)
    pdata->dir[k] /= dd;

  if(pdata->dir[0] != 0 || pdata->dir[1] != 0)
    {
      pdata->normal1[0] = -pdata->dir[1];
      pdata->normal1[1] = pdata->dir[0];
      pdata->normal1[2] = 0;
    }
  else
    {
      pdata->normal1[0] = 1;
      pdata->normal1[1] = 0;
      pdata->normal1[2] = 0;
    }

  double nn = sqrt(pdata->normal1[0] * pdata->normal1[0] + pdata->normal1[1] * pdata->normal1[1] + pdata->normal1[2] * pdata->normal1[2]);
  for(k = 0; k < 3; k++)
    pdata->normal1[k] /= nn;

  pdata->normal2[0] = pdata->dir[1] * pdata->normal1[2] - pdata->dir[2] * pdata->normal1[1];
  pdata->normal2[1] = pdata->dir[2] * pdata->normal1[0] - pdata->dir[0] * pdata->normal1[2];
  pdata->normal2[2] = pdata->dir[0] * pdata->normal1[1] - pdata->dir[1] * pdata->normal1[0];

  pdata->dx = pdata->boxx / pdata->nx;
  pdata->dy = pdata->boxy / pdata->ny;

  pdata->npixels = pdata->nx * pdata->ny;

#ifdef OPACITIES
  char * dir = "./";
  initialize_opacities(&opdata, dir);
#endif

  return 1;
}

int voronoi_proj_init(struct projection_data *pdata, int pflag)
{
  pdata->Nray = 0;
  pdata->MaxNray = pdata->npixels;
  pdata->Ray = mymalloc("Ray", pdata->MaxNray * sizeof(ray_data));

  pdata->rho = mymalloc("rho", pdata->npixels * sizeof(double));
#ifdef MHD
  pdata->bsqr = mymalloc("bsqr", pdata->npixels * sizeof(double));
#endif
  pdata->ptot = mymalloc("ptot", pdata->npixels * sizeof(double));
#ifdef GFM_DUST
  pdata->metal_rho = mymalloc("metal_rho", pdata->npixels * sizeof(double));
  pdata->dust_rho = mymalloc("dust_rho", pdata->npixels * sizeof(double));
  pdata->dust_lum = mymalloc("dust_lum", pdata->npixels * sizeof(double));
#endif
#ifdef VORONOI_PROJ_TAU
  pdata->tau = mymalloc("tau", pdata->npixels * sizeof(double));
  pdata->oldk = mymalloc("oldk", pdata->npixels * sizeof(double));
  pdata->tps = mymalloc("tps", pdata->npixels * sizeof(double));
  pdata->rhops = mymalloc("rhops", pdata->npixels * sizeof(double));
  pdata->len = mymalloc("len", pdata->npixels * sizeof(double));
  pdata->finished = mymalloc("len", pdata->npixels * sizeof(int));
#endif

  memset(pdata->rho, 0, pdata->npixels * sizeof(double));
#ifdef MHD
  memset(pdata->bsqr, 0, pdata->npixels * sizeof(double));
#endif
  memset(pdata->ptot, 0, pdata->npixels * sizeof(double));
#ifdef GFM_DUST
  memset(pdata->metal_rho, 0, pdata->npixels * sizeof(double));
  memset(pdata->dust_rho, 0, pdata->npixels * sizeof(double));
  memset(pdata->dust_lum, 0, pdata->npixels * sizeof(double));
#endif
#ifdef VORONOI_PROJ_TAU
  memset(pdata->tau, 0, pdata->npixels * sizeof(double));
  memset(pdata->oldk, 0, pdata->npixels * sizeof(double));
  memset(pdata->tps, 0, pdata->npixels * sizeof(double));
  memset(pdata->rhops, 0, pdata->npixels * sizeof(double));
  memset(pdata->len, 0, pdata->npixels * sizeof(double));
  memset(pdata->finished, 0, pdata->npixels * sizeof(int));
#endif

  /* initialize rays */
  int i, k, x, y;

  pdata->pflag = pflag;

  // first all tasks create the rays that lie in their domain
  for(x = 0; x < pdata->nx; x++)
    for(y = 0; y < pdata->nx; y++)
      {
        // set position and direction
        for(k = 0; k < 3; k++)
          {
            if(!pdata->pflag)
              {
                pdata->Ray[pdata->Nray].pos[k] = pdata->center[k]
                  - 0.5 * pdata->boxz * pdata->dir[k] + (-0.5 * pdata->boxx + (x + 0.5) * pdata->dx) * pdata->normal1[k] + (-0.5 * pdata->boxy + (y + 0.5) * pdata->dy) * pdata->normal2[k];

                pdata->Ray[pdata->Nray].dir[k] = pdata->dir[k];
              }
            else
              {
                pdata->Ray[pdata->Nray].pos[k] = pdata->center[k]
                  - 0.5 * pdata->boxx * pdata->normal1[k] + (-0.5 * pdata->boxy + (x + 0.5) * pdata->dx) * pdata->normal2[k] + (-0.5 * pdata->boxz + (y + 0.5) * pdata->dy) * pdata->dir[k];

                pdata->Ray[pdata->Nray].dir[k] = pdata->normal1[k];
              }
          }

        pdata->Ray[pdata->Nray].len = 0;
        pdata->Ray[pdata->Nray].target_len = pdata->boxz;
        pdata->Ray[pdata->Nray].pixel = y * pdata->nx + x;
        pdata->Ray[pdata->Nray].index = -1;
        pdata->Ray[pdata->Nray].prev = -1;

        // now check domain
        peanokey key = position_to_peanokey(pdata->Ray[pdata->Nray].pos);
        int no = peanokey_to_topnode(key);
        pdata->Ray[pdata->Nray].task = DomainTask[no];
        if(pdata->Ray[pdata->Nray].task == ThisTask)
          {
            // if ray is on our domain, keep it. otherwise it will be
            // overwritten by the next
            pdata->Nray++;
          }
      }

  mesh_search_data *searchdata;

  // determine the index and task of the mesh cells containing the
  // rays
  searchdata = mymalloc("searchdata", pdata->Nray * sizeof(mesh_search_data));

  // fill search array
  for(i = 0; i < pdata->Nray; i++)
    {
      searchdata[i].Pos[0] = pdata->Ray[i].pos[0];
      searchdata[i].Pos[1] = pdata->Ray[i].pos[1];
      searchdata[i].Pos[2] = pdata->Ray[i].pos[2];
    }

  find_nearest_meshpoint_global(searchdata, pdata->Nray, 0, 1);

  /* set ray index/task and send rays to correct tasks */
  mpi_distribute_items_from_search(searchdata, pdata->Ray, &pdata->Nray, &pdata->MaxNray, sizeof(ray_data), TAG_DENS_A, offsetof(ray_data, task), offsetof(ray_data, index));

  myfree(searchdata);

  return 0;
}

int voronoi_proj_run(struct projection_data *pdata)
{
  MPI_Barrier(MPI_COMM_WORLD);
  printf("Task %d has %d rays\n", ThisTask, pdata->Nray);
  MPI_Barrier(MPI_COMM_WORLD);

  measure_time();

  int i;
  // double check the positions are correct after exchanging
  for(i = 0; i < pdata->Nray; i++)
    {
      assert_contains(&Mesh, pdata->Ray[i].index, pdata->Ray[i].pos);
    }

#ifdef VORONOI_PROJ_TAU
  double *tauex = 0;
  tauex = mymalloc("tauex", pdata->npixels * sizeof(double));
  double *lenex = 0;
  lenex = mymalloc("lenex", pdata->npixels * sizeof(double));
#endif
  int left_this_task, rays_left;
  i = 0;
  do
    {
      left_this_task = voronoi_proj_advance_rays_for_one_cell(pdata);
      voronoi_proj_exchange_rays(pdata);
      ++i;
#ifdef VORONOI_PROJ_TAU
      int j;
      /* ensure that optical depth is communicated for all rays */
      memset(tauex, 0, pdata->npixels * sizeof(double));
      MPI_Allreduce(pdata->tau, tauex, pdata->npixels, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      for (j = 0; j < pdata->npixels; j++)
        pdata->tau[j] = tauex[j];
      memset(lenex, 0, pdata->npixels * sizeof(double));
      MPI_Allreduce(pdata->len, lenex, pdata->npixels, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      for (j = 0; j < pdata->npixels; j++)
        pdata->len[j] = lenex[j];
#endif
      MPI_Allreduce(&left_this_task, &rays_left, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      mpi_printf("iteration %d rays left: %d        \r", i, rays_left);
    }
  while(rays_left);

#ifdef VORONOI_PROJ_TAU
  /* set all quantities from non-finished rays to 0 on this task */
  int j;
  for (j = 0; j < pdata->npixels; j++)
    {
      if (!pdata->finished[j])
        {
          pdata->tau[j] = 0.0;
          pdata->tps[j] = 0.0;
          pdata->rhops[j] = 0.0;
          pdata->len[j] = 0.0;
        }
    }
#endif

#ifdef VORONOI_PROJ_TAU
  myfree(lenex);
  myfree(tauex);
#endif

  double *rho = 0, *ptot = 0;
#ifdef MHD
  double *bsqr = 0;
#endif
#ifdef GFM_DUST
  double *metal_rho = 0, *dust_rho = 0, *dust_lum = 0;
#endif
#ifdef VORONOI_PROJ_TAU
  double *tau = 0, *tps = 0, *rhops = 0, *len = 0;
#endif
  if(ThisTask == 0)
    {
      rho = mymalloc("rho", pdata->npixels * sizeof(double));
#ifdef MHD
      bsqr = mymalloc("bsqr", pdata->npixels * sizeof(double));
#endif
      ptot = mymalloc("ptot", pdata->npixels * sizeof(double));
#ifdef GFM_DUST
      metal_rho = mymalloc("metal_rho", pdata->npixels * sizeof(double));
      dust_rho = mymalloc("dust_rho", pdata->npixels * sizeof(double));
      dust_lum = mymalloc("dust_lum", pdata->npixels * sizeof(double));
#endif
#ifdef VORONOI_PROJ_TAU
      tau = mymalloc("tau", pdata->npixels * sizeof(double));
      tps = mymalloc("tps", pdata->npixels * sizeof(double));
      rhops = mymalloc("rhops", pdata->npixels * sizeof(double));
      len = mymalloc("len", pdata->npixels * sizeof(double));
#endif
    }

  MPI_Reduce(pdata->rho, rho, pdata->npixels, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#ifdef MHD
  MPI_Reduce(pdata->bsqr, bsqr, pdata->npixels, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
  MPI_Reduce(pdata->ptot, ptot, pdata->npixels, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#ifdef GFM_DUST
  MPI_Reduce(pdata->metal_rho, metal_rho, pdata->npixels, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(pdata->dust_rho, dust_rho, pdata->npixels, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(pdata->dust_lum, dust_lum, pdata->npixels, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
#ifdef VORONOI_PROJ_TAU
  /* get max of quantities from all cores to core 0 -> only on the cores where the rays finished,
   * quantities are > 0 */
  MPI_Reduce(pdata->tau,   tau,   pdata->npixels, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(pdata->tps,   tps,   pdata->npixels, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(pdata->rhops, rhops, pdata->npixels, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(pdata->len,   len,   pdata->npixels, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
#endif

  double t = measure_time();
  mpi_printf("Integration done after %d iterations and %f s\n", i, t);

  if(ThisTask == 0)
    {
      voronoi_proj_write_file("rho", pdata, rho);
#ifdef MHD
      voronoi_proj_write_file("bsqr", pdata, bsqr);
#endif
      voronoi_proj_write_file("ptot", pdata, ptot);
#ifdef GFM_DUST
      voronoi_proj_write_file("metal_rho", pdata, metal_rho);
      voronoi_proj_write_file("dust_rho", pdata, dust_rho);
      voronoi_proj_write_file("dust_lum", pdata, dust_lum);
#endif
#ifdef VORONOI_PROJ_TAU
      voronoi_proj_write_file("tau", pdata, tau);
      voronoi_proj_write_file("tps", pdata, tps);
      voronoi_proj_write_file("rhops", pdata, rhops);
      voronoi_proj_write_file("len", pdata, len);
#endif

#ifdef VORONOI_PROJ_TAU
      myfree(len);
      myfree(rhops);
      myfree(tps);
      myfree(tau);
#endif
#ifdef GFM_DUST
      myfree(dust_lum);
      myfree(dust_rho);
      myfree(metal_rho);
#endif
      myfree(ptot);
#ifdef MHD
      myfree(bsqr);
#endif
      myfree(rho);
    }

  return 0;
}

void voronoi_proj_write_file(char *name, struct projection_data *pdata, double *data)
{
  FILE *fd;
  char fname[1000];

  int res[2];
  double box[3];

  res[0] = pdata->nx;
  res[1] = pdata->ny;
  box[0] = pdata->boxx;
  box[1] = pdata->boxy;
  box[2] = pdata->boxz;

  if(!pdata->pflag)
    sprintf(fname, "%s/vproj_%s_%03d", All.OutputDir, name, RestartSnapNum);
  else
    sprintf(fname, "%s/vprojp_%s_%03d", All.OutputDir, name, RestartSnapNum);

  printf("Writing file %s.\n", fname);

  fd = fopen(fname, "w");
  fwrite(res, sizeof(int), 2, fd);
  fwrite(box, sizeof(double), 3, fd);
  fwrite(data, sizeof(double), pdata->npixels, fd);
  fclose(fd);
}

int voronoi_proj_deinit(struct projection_data *pdata)
{
#ifdef VORONOI_PROJ_TAU
  myfree(pdata->finished);
  myfree(pdata->len);
  myfree(pdata->rhops);
  myfree(pdata->tps);
  myfree(pdata->oldk);
  myfree(pdata->tau);
#endif
#ifdef GFM_DUST
  myfree(pdata->dust_lum);
  myfree(pdata->dust_rho);
  myfree(pdata->metal_rho);
#endif
  myfree(pdata->ptot);
#ifdef MHD
  myfree(pdata->bsqr);
#endif
  myfree(pdata->rho);
  myfree(pdata->Ray);

  return 0;
}

void voronoi_proj_exchange_rays(struct projection_data *pdata)
{
  int i;
  for(i = 0; i < pdata->Nray; i++)
    {
      if(pdata->Ray[i].task != ThisTask)
        pdata->Ray[i].prev = -1;
    }

  mpi_distribute_items_to_tasks(pdata->Ray, offsetof(ray_data, task), &pdata->Nray, &pdata->MaxNray, sizeof(*pdata->Ray), TAG_DENS_A);
}

void voronoi_proj_update_column_integral(int cell, double len, int pix, double delta_center[3], struct projection_data *pdata)
{
#ifdef PROJ_WRITE_RAYS
  char fname[200];
  FILE *fd;
#endif
  //SphP[sph_idx].Grad.drho[0] * delta_center[0] + SphP[sph_idx].Grad.drho[1] * delta_center[1] + SphP[sph_idx].Grad.drho[2] * delta_center[2];
#ifdef MHD
  double bsqr = (SphP[cell].B[0] * SphP[cell].B[0] + SphP[cell].B[1] * SphP[cell].B[1] + SphP[cell].B[2] * SphP[cell].B[2]);
#endif

  pdata->rho[pix] += len * SphP[cell].Density;
#ifdef MHD
  pdata->bsqr[pix] += len * bsqr;
  pdata->ptot[pix] += len * (SphP[cell].Pressure + 0.5 * bsqr);
#endif
#ifdef GFM_DUST
  pdata->metal_rho[pix] += len * SphP[cell].Density * SphP[cell].Metallicity;

  double dustZ = 0.0;
  int chan, k_elem;
  for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
    {
      for(k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
        {
          dustZ += SphP[cell].MetalsDustFraction[chan][k_elem];
        }
    }
  pdata->dust_rho[pix] += len * SphP[cell].Density * dustZ;

  double mu = 4.0 / (1. + 3. * HYDROGEN_MASSFRAC + 4. * HYDROGEN_MASSFRAC * SphP[cell].Ne);

  double temp_in_K = GAMMA_MINUS1 * SphP[cell].Utherm / BOLTZMANN * All.UnitEnergy_in_cgs / All.UnitMass_in_g * mu * PROTONMASS;

  double dust_to_gas_ratio;
  char dust_cool; 

  get_CoolingDustState(cell, &dust_to_gas_ratio, &dust_cool);

  double nHcgs = HYDROGEN_MASSFRAC * SphP[cell].Density * All.UnitDensity_in_cgs / PROTONMASS  * All.HubbleParam*All.HubbleParam;
  double vol = SphP[cell].Volume * All.UnitLength_in_cm * All.UnitLength_in_cm * All.UnitLength_in_cm / (All.HubbleParam*All.HubbleParam*All.HubbleParam);
  if (temp_in_K >0)
    pdata->dust_lum[pix] += get_CoolingDustRate(log10(temp_in_K), dust_to_gas_ratio, mu, dust_cool) * nHcgs * nHcgs * vol; 
#endif
#ifdef VORONOI_PROJ_TAU
  if (!pdata->finished[pix])
    {
      double rho = SphP[cell].Density + SphP[cell].Grad.drho[0] * delta_center[0] + SphP[cell].Grad.drho[1] * delta_center[1] + SphP[cell].Grad.drho[2] * delta_center[2];
      rho = dmax(1e-14, rho);
      /* use cell value */
      //double t = SphP[cell].Temperature;
      /* use T gradient */
      double t = SphP[cell].Temperature + SphP[cell].Grad.dtemp[0] * delta_center[0] + SphP[cell].Grad.dtemp[1] * delta_center[1] + SphP[cell].Grad.dtemp[2] * delta_center[2];
      t = dmax(1e3, t);
      double kappa = opacities(rho, t, &opdata);
      double tauold = pdata->tau[pix];
      pdata->tau[pix] += 0.5 * len * (pdata->oldk[pix] + kappa * rho);
      pdata->oldk[pix] = kappa * rho;
      pdata->len[pix] += len;
      if (pdata->tau[pix] > 1.0)
        {
          pdata->tps[pix] = t;
          pdata->rhops[pix] = rho;
          pdata->finished[pix] = 1;
          /*
          // interpolate linealy back to tau=1
          // x = x1 + (x2-x1) / (tau2-tau1) * (tau-tau1)
          double deltalen = len - len / (pdata->tau[pix] - tauold) * (1.0 - tauold);
          int j;
          for(j = 0; j < 3; j++)
            {
              delta_center[j] = pdata->Ray[pix].pos[j] - (0.5 * len + deltalen) * pdata->Ray[pix].dir[j] - SphP[cell].Center[j];
            }
          //printf("Ray %d reached tau = 1\n", pix);
          pdata->tps[pix] = SphP[cell].Temperature + SphP[cell].Grad.dtemp[0] * delta_center[0] + SphP[cell].Grad.dtemp[1] * delta_center[1] + SphP[cell].Grad.dtemp[2] * delta_center[2];
          pdata->rhops[pix] = SphP[cell].Density + SphP[cell].Grad.drho[0] * delta_center[0] + SphP[cell].Grad.drho[1] * delta_center[1] + SphP[cell].Grad.drho[2] * delta_center[2];
          */
        }
#ifdef PROJ_WRITE_RAYS
      if(!pdata->pflag)
        sprintf(fname, "ray_%06d", pix);
      else
        sprintf(fname, "rayp_%06d", pix);

      fd = fopen(fname, "a");
      fprintf(fd, "%e %e %e %e %e\n", rho, t, kappa, pdata->tau[pix], pdata->len[pix]);
      fclose(fd);
#endif
    }
#endif
}

int voronoi_proj_advance_rays_for_one_cell(struct projection_data *pdata)
{
#if !defined(VORONOI) || !defined(VORONOI_DYNAMIC_UPDATE)
  terminate("Ray tracing only works with VORONOI and DYNAMIC_UPDATE enabled.");
#else
  int i, j, nleft = 0;

  for(i = 0; i < pdata->Nray; i++)
    {
      if(pdata->Ray[i].len >= pdata->Ray[i].target_len)
        // ray done
        continue;

      /* this is the index of the cell in which we currently are */
      int cell = pdata->Ray[i].index;
      double len;

      int next_edge = find_next_voronoi_cell(&Mesh, cell, pdata->Ray[i].pos, pdata->Ray[i].dir,
                                             pdata->Ray[i].prev, &len);

      if(pdata->Ray[i].len + len > pdata->Ray[i].target_len)
        {
          /* The next cell boundary is past the target length. truncate
             propagation and mark ray as done. do not change the
             cell/task. set prev to -1 since we are inside the cell */
          len = pdata->Ray[i].target_len - pdata->Ray[i].len;
          pdata->Ray[i].len = pdata->Ray[i].target_len;
          pdata->Ray[i].prev = -1;
        }
      else
        {
          /* we will enter next cell. update cell info */
          pdata->Ray[i].task = DC[next_edge].task;
          pdata->Ray[i].prev = pdata->Ray[i].index;
          pdata->Ray[i].index = DC[next_edge].index;
          pdata->Ray[i].len += len;
          ++nleft;
        }

      /* now determine the interpolated value at midpoint of segment */
      double delta_center[3];

#ifndef VORONOI_PROJ_SUBSTEPS
      for(j = 0; j < 3; j++)
        {
          delta_center[j] = pdata->Ray[i].pos[j] + 0.5 * len * pdata->Ray[i].dir[j] - SphP[cell].Center[j];
          pdata->Ray[i].pos[j] += len * pdata->Ray[i].dir[j];
        }

      voronoi_proj_update_column_integral(cell, len, pdata->Ray[i].pixel, delta_center, pdata);
#else
      double fac = 3.0;
      if (pdata->tau[pdata->Ray[i].pixel] > 0.001 && pdata->tau[pdata->Ray[i].pixel] < 1.0)
        fac = 10.0;
      double sublen = len / (VORONOI_PROJ_SUBSTEPS * fac);
      int k;
      for(k = 1; k < VORONOI_PROJ_SUBSTEPS * fac + 1; k++)
        {
          for(j = 0; j < 3; j++)
            {
              delta_center[j] = pdata->Ray[i].pos[j] + 0.5 * sublen * pdata->Ray[i].dir[j] - SphP[cell].Center[j];
              pdata->Ray[i].pos[j] += sublen * pdata->Ray[i].dir[j];
            }

          voronoi_proj_update_column_integral(cell, sublen, pdata->Ray[i].pixel, delta_center, pdata);
        }
      if (pdata->finished[pdata->Ray[i].pixel])
        {
          /* mark ray as finished */
          pdata->Ray[i].len = pdata->Ray[i].target_len;
          pdata->Ray[i].prev = -1;
        }
#endif

      // see if we need to wrap ray around boundary. (This is hoaky,
      // because if the point struct is changed it will break.)
      periodic_wrap_point(pdata->Ray[i].pos, &Mesh.DP[DC[next_edge].dp_index].x);
    }
  return nleft;
#endif
}
#endif
