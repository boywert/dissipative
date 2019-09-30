/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dvr_render/dvr_render.c
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
#include <sys/stat.h>
#include <sys/types.h>

#include "../allvars.h"
#include "../proto.h"
#include "../voronoi.h"
#include "dvr_render.h"

#ifdef DVR_RENDER
#ifndef VORONOI_DYNAMIC_UPDATE
#error DVR_RENDER requires VORONOI_DYNAMIC_UPDATE
#endif
#endif

#ifdef DVR_RENDER


/* select kernel for smoothing */
//#define DVR_TOPHAT


static int dvr_do_next_camera_setup(FILE * fdcam, dvr_cam_data * dvrdata);
static int dvr_init(dvr_cam_data * dvrdata);
static int dvr_advance_rays_for_one_cell(dvr_cam_data * dvrdata);
static int dvr_run(dvr_cam_data * dvrdata);
static int dvr_deinit(dvr_cam_data * dvrdata);
static void dvr_write_file(char *name, dvr_cam_data * dvrdata, double *data, int field);
static void dvr_exchange_rays(dvr_cam_data * dvrdata);
static void dvr_update_render_integral(int cell, double len, dvr_ray_data * ray, dvr_cam_data * dvrdata, int field);

/* camera path file */
static FILE *fdcam = NULL;

/* entries in transfer function table */
static int transRGBEntries[DVR_NUM_FIELDS];

/* transfer function table */
static transRGB_data transRGBTable[DVR_NUM_FIELDS][DVR_NUM_TRANS];

/* set render field data type for casting */
typedef MyFloat RenderDataType;


static inline void renderValueSphP(int index, int field, double *rgba_voxel_data)
{
  rgba_voxel_data[0] = SphP[index].DvrFields[field];
  rgba_voxel_data[1] = SphP[index].DvrFields[DVR_RENDER_RHO];
}

static inline void renderValuePrimExch(int index, int field, double *rgba_voxel_data)
{
  rgba_voxel_data[0] = PrimExch[index].DvrFields[field];
  rgba_voxel_data[1] = PrimExch[index].DvrFields[DVR_RENDER_RHO];
}


/* kernel for smoothing */
static inline float kernel(float u)
{
  if(u < 0.5)
    return (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
  if(u < 1.0)
    return KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
  return 0;
}


/* mark frustum for mesh generation */
void mark_frustum_for_tessellation(dvr_cam_data * dvrdata)
{
  double off_vector[3];
  int i;
  double delta_x, delta_y, delta_z;
  int local_cells = 0;
  long long tot_cells;
#ifndef DVR_TESS_ALL
  double distscale = 0;
#endif

  for(i = 0; i < NumGas; i++)
    {
      /* off vector */
      off_vector[0] = P[i].Pos[0] - dvrdata->camPos[0];
      off_vector[1] = P[i].Pos[1] - dvrdata->camPos[1];
      off_vector[2] = P[i].Pos[2] - dvrdata->camPos[2];

      /* project into camera coordinate system (x: along viewing direction) note: all dvrdata vectors are normalized */
      delta_x = dvrdata->dir[0] * off_vector[0] + dvrdata->dir[1] * off_vector[1] + dvrdata->dir[2] * off_vector[2];
      delta_y = fabs(dvrdata->upVec[0] * off_vector[0] + dvrdata->upVec[1] * off_vector[1] + dvrdata->upVec[2] * off_vector[2]);
      delta_z = fabs(dvrdata->normalVec[0] * off_vector[0] + dvrdata->normalVec[1] * off_vector[1] + dvrdata->normalVec[2] * off_vector[2]);

      SphP[i].DoMesh = 1;

#ifndef DVR_TESS_ALL
      if(delta_x > (1.0 + DVR_FRUST_BUFF_ZONE) * dvrdata->far || delta_x < (1.0 - DVR_FRUST_BUFF_ZONE) * dvrdata->near)
        {
          SphP[i].DoMesh = 0;
          continue;
        }

      distscale = (delta_x + dvrdata->near) / dvrdata->near;

      if(2 * delta_y / dvrdata->height > (1.0 + DVR_FRUST_BUFF_ZONE) * distscale)
        {
          SphP[i].DoMesh = 0;
          continue;
        }

      if(2 * delta_z / dvrdata->width > (1.0 + DVR_FRUST_BUFF_ZONE) * distscale)
        {
          SphP[i].DoMesh = 0;
          continue;
        }
#endif
      local_cells++;
    }

  sumup_large_ints(1, &local_cells, &tot_cells);

  mpi_printf("DVR_RENDER: marked cells in render frustum = %lld\n", tot_cells);
}


/* read in transfer function table */
void dvr_load_transfer_function(char *fname, int field)
{
  FILE *fdtrans;
  double dummy;
  int entries;

  transRGBEntries[field] = 0;
  if(!(fdtrans = fopen(fname, "r")))
    terminate("DVR_RENDER: Cannot read transfer function file `%s'\n", fname);
  while(1)
    {
      if(fscanf(fdtrans, "%lf %lf %lf %lf", &dummy, &dummy, &dummy, &dummy) == EOF)
        break;
      transRGBEntries[field]++;
    }
  fclose(fdtrans);

  if(transRGBEntries[field] > DVR_NUM_TRANS)
    terminate("DVR_RENDER: too many entries for transfer function");

  fdtrans = fopen(fname, "r");
  for(entries = 0; entries < transRGBEntries[field]; entries++)
    {
      fscanf(fdtrans, "%lf %lf %lf %lf", &(transRGBTable[field][entries].Val), &(transRGBTable[field][entries].Red), &(transRGBTable[field][entries].Green), &(transRGBTable[field][entries].Blue));
      mpi_printf("DVR_RENDER: transRGBTable: %d/%d %g %g %g %g\n", entries, transRGBEntries[field], transRGBTable[field][entries].Val,
                 transRGBTable[field][entries].Red, transRGBTable[field][entries].Green, transRGBTable[field][entries].Blue);
    }

  fclose(fdtrans);
}


#ifdef DVR_RENDER_SMOOTH
/* make scalar field smoother by taking some average over connected neighbor cells */
static void smooth_cell(int cell, double raypos[3], int field, double *rgba_voxel_data)
{
  int q, channel, iter;
  double tmp_rgba_voxel_data[2];
  float dx, dy, dz;
  float wk;
  float tot_norm, volume;
  float r2;

#ifndef DVR_TOPHAT
  float hinv, hsml2;
#endif

  tot_norm = 0.0;
  rgba_voxel_data[0] = rgba_voxel_data[1] = 0.0;

#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif

  dx = NEAREST_X(raypos[0] - SphP[cell].Center[0]);
  dy = NEAREST_Y(raypos[1] - SphP[cell].Center[1]);
  dz = NEAREST_Z(raypos[2] - SphP[cell].Center[2]);

#ifndef DVR_TOPHAT
  hsml2 = dx * dx + dy * dy + dz * dz;
#endif

  for(iter = 0; iter < 2; iter++)
    {
      if(iter == 1)
        {
          dx = NEAREST_X(raypos[0] - SphP[cell].Center[0]);
          dy = NEAREST_Y(raypos[1] - SphP[cell].Center[1]);
          dz = NEAREST_Z(raypos[2] - SphP[cell].Center[2]);

          r2 = dx * dx + dy * dy + dz * dz;

#ifndef DVR_TOPHAT
          wk = kernel(sqrtf(r2) * hinv);
#else
          wk = 1.0;
#endif
          renderValueSphP(cell, field, tmp_rgba_voxel_data);
          volume = SphP[cell].Volume;
          rgba_voxel_data[0] += tmp_rgba_voxel_data[0] * volume * wk;
          rgba_voxel_data[1] += tmp_rgba_voxel_data[1] * volume * wk;
          tot_norm += volume * wk;
        }

      q = SphP[cell].first_connection;
      while(q >= 0)
        {
          int dp = DC[q].dp_index;
          int part = Mesh.DP[dp].index;

          if(part < 0)
            {
              if(q == SphP[cell].last_connection)
                break;
              q = DC[q].next;
              continue;
            }

          if(part >= NumGas && Mesh.DP[dp].task == ThisTask)
            part -= NumGas;

          if(DC[q].task == ThisTask)
            {
              dx = NEAREST_X(raypos[0] - SphP[part].Center[0]);
              dy = NEAREST_Y(raypos[1] - SphP[part].Center[1]);
              dz = NEAREST_Z(raypos[2] - SphP[part].Center[2]);
              if(iter == 1)
                {
                  renderValueSphP(part, field, tmp_rgba_voxel_data);
                  volume = SphP[part].Volume;
                }
            }
          else
            {
              dx = NEAREST_X(raypos[0] - PrimExch[part].Center[0]);
              dy = NEAREST_Y(raypos[1] - PrimExch[part].Center[1]);
              dz = NEAREST_Z(raypos[2] - PrimExch[part].Center[2]);
              if(iter == 1)
                {
                  renderValuePrimExch(part, field, tmp_rgba_voxel_data);
                  volume = PrimExch[part].Volume;
                }
            }

          r2 = dx * dx + dy * dy + dz * dz;

          if(iter == 0)
            {
#ifndef DVR_TOPHAT
              if(r2 > hsml2)
                hsml2 = r2;
#endif
            }
          else
            {
#ifndef DVR_TOPHAT
              wk = kernel(sqrtf(r2) * hinv);
#else
              wk = 1.0;
#endif
              rgba_voxel_data[0] += tmp_rgba_voxel_data[0] * volume * wk;
              rgba_voxel_data[1] += tmp_rgba_voxel_data[1] * volume * wk;
              tot_norm += volume * wk;
            }

          if(q == SphP[cell].last_connection)
            break;

          q = DC[q].next;
        }
#ifndef DVR_TOPHAT
      hinv = 1.0f / sqrtf(hsml2);
#endif
    }
  tot_norm = 1.0f / tot_norm;
  rgba_voxel_data[0] *= tot_norm;
  rgba_voxel_data[1] *= tot_norm;
}
#endif


/* voxel data -> RGBA assignment (transfer function evaluation) */
static void transferFunctions(double *rgba_voxel_data, double len, int field, color_data * color)
{
  double tau;
  int entries;
  double value;

  /* ALPHA: opacity transfer function */
  tau = (rgba_voxel_data[1] * len >= 0) ? (rgba_voxel_data[1] * len) : (0);
  color->Alpha = 1.0 - exp(-All.DvrTauScaleFactor * tau);       //All.DvrTauScaleFactor * tau + All.DvrTauFloor;      /* ~ 1.0 - exp(-All.DvrTauScaleFactor * tau) */

  /* RGB: transfer functions */
  /* default black */
  color->Red = 0.0;
  color->Green = 0.0;
  color->Blue = 0.0;

  if(rgba_voxel_data[0] > 0)
    {
      value = log10(rgba_voxel_data[0]);

      for(entries = 1; entries < transRGBEntries[field]; entries++)
        if(transRGBTable[field][entries - 1].Val < value && transRGBTable[field][entries].Val > value)
          break;

      if(entries >= 1 && entries < transRGBEntries[field])
        {
          double fac = (value - transRGBTable[field][entries - 1].Val) / (transRGBTable[field][entries].Val - transRGBTable[field][entries - 1].Val);
          color->Red = transRGBTable[field][entries - 1].Red + fac * (transRGBTable[field][entries].Red - transRGBTable[field][entries - 1].Red);
          color->Green = transRGBTable[field][entries - 1].Green + fac * (transRGBTable[field][entries].Green - transRGBTable[field][entries - 1].Green);
          color->Blue = transRGBTable[field][entries - 1].Blue + fac * (transRGBTable[field][entries].Blue - transRGBTable[field][entries - 1].Blue);
        }
    }
}

void dvr_render_main(void)
{
  double deltat = get_time_difference_in_Gyr(All.TimeBegin, All.Time);
  int newframenumber = (int) (deltat / All.DvrRenderTimeIntverallInGyr);

  mpi_printf("DVR_RENDER: check... timediff=%g  ratio=%g\n", deltat, deltat / All.DvrRenderTimeIntverallInGyr);
  if(newframenumber == All.DvrFrameFileNum)
    mpi_printf("DVR_RENDER: check... no rendering (%d|%d)\n", newframenumber, All.DvrFrameFileNum);
  if(newframenumber == All.DvrFrameFileNum + 1)
    {
      mpi_printf("DVR_RENDER: check... time to render (%d|%d)\n", newframenumber, All.DvrFrameFileNum);
      dvr_render_loop("campath.txt", All.DvrPixelX, All.DvrPixelY, 1);
    }
}

/* main render loop, walk through camera file and shoot rays from near field plane */
int dvr_render_loop(char *fname, int pixelx, int pixely, int mode)
{
  dvr_cam_data dvrdata;

  dvrdata.nx = pixelx;
  dvrdata.ny = pixely;

  if(mode == 1)
    fname = "campath.txt";

  if(fdcam == NULL)
    {
      if(!(fdcam = fopen(fname, "r")))
        terminate("DVR_RENDER: Cannot read camera setup file `%s'\n", fname);

      mpi_printf("DVR_RENDER: opened camera path file `%s'\n", fname);
      int i;
      double dummy;
      for(i = 0; i < All.DvrFrameFileNum; i++)
        if(fscanf(fdcam, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy) == EOF)
          terminate("DVR_RENDER: camera path file corrupted.");
      mpi_printf("DVR_RENDER: jumped to camera position %d\n", All.DvrFrameFileNum);
    }


  if(mode == 0)
    {
      while(1)
        {
          if(dvr_do_next_camera_setup(fdcam, &dvrdata) < 0)
            break;
          All.DvrFrameFileNum++;
        }
    }
  if(mode == 1)
    {
      dvr_do_next_camera_setup(fdcam, &dvrdata);
      All.DvrFrameFileNum++;
    }

  if(mode == 0)
    {
      fclose(fdcam);
    }

  return 0;
}


static void setup_transfer_function(void)
{
  int field;

  for(field = 0; field < DVR_NUM_FIELDS; field++)
    {
      if(field == DVR_RENDER_TEMP)
        dvr_load_transfer_function("dvr_render/transRGB_temp.txt", DVR_RENDER_TEMP);

      if(field == DVR_RENDER_MET)
        dvr_load_transfer_function("dvr_render/transRGB_met.txt", DVR_RENDER_MET);

      if(field == DVR_RENDER_RHO)
        dvr_load_transfer_function("dvr_render/transRGB_rho.txt", DVR_RENDER_RHO);

      if(field == DVR_RENDER_XRAY)
        dvr_load_transfer_function("dvr_render/transRGB_xray.txt", DVR_RENDER_XRAY);

      if(field >= DVR_RENDER_ELEMENT_0 && field <= DVR_RENDER_ELEMENT_8)
        {
          int off = field - DVR_RENDER_ELEMENT_0;
          char buf[MAXLEN_PATH];
          sprintf(buf, "dvr_render/transRGB_element_%d.txt", off);
          dvr_load_transfer_function(buf, field);
        }

#ifdef GFM_DUST
      if(field == DVR_RENDER_DUST_RHO)
        {
          dvr_load_transfer_function("dvr_render/transRGB_dust_rho.txt", DVR_RENDER_DUST_RHO);
        }
      if(field == DVR_RENDER_DUST_MET)
        {
          dvr_load_transfer_function("dvr_render/transRGB_dust_met.txt", DVR_RENDER_DUST_MET);
        }
#endif
    }
}

/* render next camera setup */
static int dvr_do_next_camera_setup(FILE * fdcam, dvr_cam_data * dvrdata)
{
  int k;

  if(fscanf(fdcam, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
            &(dvrdata->camPos[0]), &(dvrdata->camPos[1]), &(dvrdata->camPos[2]),
            &(dvrdata->camTo[0]), &(dvrdata->camTo[1]), &(dvrdata->camTo[2]),
            &(dvrdata->width), &(dvrdata->height), &(dvrdata->near), &(dvrdata->far), &(dvrdata->upVec[0]), &(dvrdata->upVec[1]), &(dvrdata->upVec[2])) == EOF)
    return -1;

  /* direction vector of view axis */
  for(k = 0; k < 3; k++)
    dvrdata->dir[k] = dvrdata->camTo[k] - dvrdata->camPos[k];
  double dd = sqrt(dvrdata->dir[0] * dvrdata->dir[0] + dvrdata->dir[1] * dvrdata->dir[1] + dvrdata->dir[2] * dvrdata->dir[2]);
  for(k = 0; k < 3; k++)
    dvrdata->dir[k] /= dd;

  /* normalize up vector */
  double nn = sqrt(dvrdata->upVec[0] * dvrdata->upVec[0] + dvrdata->upVec[1] * dvrdata->upVec[1] + dvrdata->upVec[2] * dvrdata->upVec[2]);
  for(k = 0; k < 3; k++)
    dvrdata->upVec[k] /= nn;

  dvrdata->normalVec[0] = dvrdata->dir[1] * dvrdata->upVec[2] - dvrdata->dir[2] * dvrdata->upVec[1];
  dvrdata->normalVec[1] = dvrdata->dir[2] * dvrdata->upVec[0] - dvrdata->dir[0] * dvrdata->upVec[2];
  dvrdata->normalVec[2] = dvrdata->dir[0] * dvrdata->upVec[1] - dvrdata->dir[1] * dvrdata->upVec[0];

  dvrdata->dx = dvrdata->width / dvrdata->nx;
  dvrdata->dy = dvrdata->height / dvrdata->ny;

  dvrdata->npixels = dvrdata->nx * dvrdata->ny;


  mpi_printf("DVR_RENDER: -------------------------- > RENDER FRAME NUMBER=%04d\n", All.DvrFrameFileNum);
  mpi_printf("DVR_RENDER: pixelx/pixely=%d/%d\n", dvrdata->nx, dvrdata->ny);
  mpi_printf("DVR_RENDER: camera position=(%g|%g|%g)\n", dvrdata->camPos[0], dvrdata->camPos[1], dvrdata->camPos[2]);
  mpi_printf("DVR_RENDER: camera lookat position=(%g|%g|%g)\n", dvrdata->camTo[0], dvrdata->camTo[1], dvrdata->camTo[2]);
  mpi_printf("DVR_RENDER: width/height=%g/%g\n", dvrdata->width, dvrdata->height);
  mpi_printf("DVR_RENDER: near/far=%g/%g\n", dvrdata->near, dvrdata->far);
  mpi_printf("DVR_RENDER: normalized viewing vector=(%g|%g|%g)\n", dvrdata->dir[0], dvrdata->dir[1], dvrdata->dir[2]);
  mpi_printf("DVR_RENDER: upVec=(%g|%g|%g)\n", dvrdata->upVec[0], dvrdata->upVec[1], dvrdata->upVec[2]);
  mpi_printf("DVR_RENDER: normalVec=(%g|%g|%g)\n", dvrdata->normalVec[0], dvrdata->normalVec[1], dvrdata->normalVec[2]);


  double temp, u, ne, rho, dens;
#ifdef GFM_COOLING_METAL
  double mu;
#endif

  int i;

  /* new mesh */
  mark_frustum_for_tessellation(dvrdata);

  voronoi_make_new_tessellation();

  /* set fields */
  for(i = 0; i < NumGas; i++)
    {
      if(P[i].Mass <= 0)
        continue;

      //DVR_RENDER_TEMP
#if (DVR_RENDER_TEMP < DVR_NUM_FIELDS)
      u = SphP[i].Utherm * All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;
#ifdef COOLING
      ne = SphP[i].Ne;
      rho = SphP[i].Density * All.cf_a3inv;
      dens = rho * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
#ifdef GFM_COOLING_METAL
      mu = (4. / (1. + 3. * SphP[i].MetalsFraction[element_index_Hydrogen] + 4. * SphP[i].MetalsFraction[element_index_Hydrogen] * ne) * PROTONMASS);
      update_gas_state(rho, SphP[i].MetalsFraction[element_index_Hydrogen], SphP[i].Metallicity);
#endif
      temp = convert_u_to_temp(u, dens, &ne);
#else
      temp = u;
#endif
      SphP[i].DvrFields[DVR_RENDER_TEMP] = temp;
#endif

      //DVR_RENDER_MET
#ifdef GFM_STELLAR_EVOLUTION
#if (DVR_RENDER_MET < DVR_NUM_FIELDS)
      SphP[i].DvrFields[DVR_RENDER_MET] = SphP[i].MassMetallicity / P[i].Mass;
#endif
#endif

      //DVR_RENDER_RHO
#if (DVR_RENDER_RHO < DVR_NUM_FIELDS)
      SphP[i].DvrFields[DVR_RENDER_RHO] = P[i].Mass / SphP[i].Volume;
#endif

      //DVR_RENDER_XRAY
#if (DVR_RENDER_XRAY < DVR_NUM_FIELDS)
      if(temp > DVR_XRAY_TEMP)
        SphP[i].DvrFields[DVR_RENDER_XRAY] = 1.2e-24 * P[i].Mass / All.HubbleParam * All.UnitMass_in_g * dens / (mu * mu) * sqrt(temp / 1.1e7);
      else
        SphP[i].DvrFields[DVR_RENDER_XRAY] = 0.0;
#endif

      //DVR_RENDER_ELEMENT_0-DVR_RENDER_ELEMENT_8
#ifdef GFM_STELLAR_EVOLUTION
#if (DVR_RENDER_ELEMENT_0 < DVR_NUM_FIELDS)
      SphP[i].DvrFields[DVR_RENDER_ELEMENT_0] = SphP[i].MassMetals[DVR_RENDER_ELEMENT_0 - DVR_RENDER_ELEMENT_0] / P[i].Mass;
#endif
#if (DVR_RENDER_ELEMENT_1 < DVR_NUM_FIELDS)
      SphP[i].DvrFields[DVR_RENDER_ELEMENT_1] = SphP[i].MassMetals[DVR_RENDER_ELEMENT_1 - DVR_RENDER_ELEMENT_0] / P[i].Mass;
#endif
#if (DVR_RENDER_ELEMENT_2 < DVR_NUM_FIELDS)
      SphP[i].DvrFields[DVR_RENDER_ELEMENT_2] = SphP[i].MassMetals[DVR_RENDER_ELEMENT_2 - DVR_RENDER_ELEMENT_0] / P[i].Mass;
#endif
#if (DVR_RENDER_ELEMENT_3 < DVR_NUM_FIELDS)
      SphP[i].DvrFields[DVR_RENDER_ELEMENT_3] = SphP[i].MassMetals[DVR_RENDER_ELEMENT_3 - DVR_RENDER_ELEMENT_0] / P[i].Mass;
#endif
#if (DVR_RENDER_ELEMENT_4 < DVR_NUM_FIELDS)
      SphP[i].DvrFields[DVR_RENDER_ELEMENT_4] = SphP[i].MassMetals[DVR_RENDER_ELEMENT_4 - DVR_RENDER_ELEMENT_0] / P[i].Mass;
#endif
#if (DVR_RENDER_ELEMENT_5 < DVR_NUM_FIELDS)
      SphP[i].DvrFields[DVR_RENDER_ELEMENT_5] = SphP[i].MassMetals[DVR_RENDER_ELEMENT_5 - DVR_RENDER_ELEMENT_0] / P[i].Mass;
#endif
#if (DVR_RENDER_ELEMENT_6 < DVR_NUM_FIELDS)
      SphP[i].DvrFields[DVR_RENDER_ELEMENT_6] = SphP[i].MassMetals[DVR_RENDER_ELEMENT_6 - DVR_RENDER_ELEMENT_0] / P[i].Mass;
#endif
#if (DVR_RENDER_ELEMENT_7 < DVR_NUM_FIELDS)
      SphP[i].DvrFields[DVR_RENDER_ELEMENT_7] = SphP[i].MassMetals[DVR_RENDER_ELEMENT_7 - DVR_RENDER_ELEMENT_0] / P[i].Mass;
#endif
#if (DVR_RENDER_ELEMENT_8 < DVR_NUM_FIELDS)
      SphP[i].DvrFields[DVR_RENDER_ELEMENT_8] = SphP[i].MassMetals[DVR_RENDER_ELEMENT_8 - DVR_RENDER_ELEMENT_0] / P[i].Mass;
#endif
#endif

#ifdef GFM_DUST
      double dustZ = 0.0;
      int chan, k_elem;
      for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
        {
          for(k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
            {
              dustZ += SphP[i].MassMetalsDust[chan][k_elem] / P[i].Mass;
            }
        }

#if (DVR_RENDER_DUST_RHO < DVR_NUM_FIELDS)
      SphP[i].DvrFields[DVR_RENDER_DUST_RHO] = dustZ * P[i].Mass / SphP[i].Volume;
#endif
#if (DVR_RENDER_DUST_MET < DVR_NUM_FIELDS)
      SphP[i].DvrFields[DVR_RENDER_DUST_MET] = dustZ;
#endif
#endif
    }

  /* gradients */
  exchange_primitive_variables();
  calculate_gradients();
  exchange_primitive_variables_and_gradients();

  /* render */
  CPU_Step[CPU_MISC] += measure_time();

  setup_transfer_function();
  dvr_init(dvrdata);
  dvr_run(dvrdata);
  dvr_deinit(dvrdata);

  CPU_Step[CPU_MAKEIMAGES] += measure_time();

  /* restore old time line */
  voronoi_restore_old_tessellation();

  return 0;
}


static int dvr_init(dvr_cam_data * dvrdata)
{
  int i, k, x, y;
  double dd;
  int field;

  dvrdata->Nray = 0;
  dvrdata->MaxNray = dvrdata->npixels;
  dvrdata->Ray = mymalloc("Ray", dvrdata->MaxNray * sizeof(dvr_ray_data));

  for(i = 0; i < dvrdata->MaxNray; i++)
    for(field = 0; field < DVR_NUM_FIELDS; field++)
      dvrdata->Ray[i].Red[field] = dvrdata->Ray[i].Green[field] = dvrdata->Ray[i].Blue[field] = dvrdata->Ray[i].Opacity[field] = 0.0;

  /* all tasks create the rays that lie in their domain */
  for(x = 0; x < dvrdata->nx; x++)
    for(y = 0; y < dvrdata->ny; y++)
      {
        dd = 0;
        for(k = 0; k < 3; k++)
          {
            /* rays start at near field plane -> we do front-to-back alpha composite rendering */
            dvrdata->Ray[dvrdata->Nray].pos[k] = dvrdata->camPos[k] +   /* camera */
              dvrdata->near * dvrdata->dir[k] + /* center of near plane */
              (-0.5 * dvrdata->width + (x + 0.5) * dvrdata->dx) * dvrdata->normalVec[k] +       /* offset left-right  */
              (-0.5 * dvrdata->height + (y + 0.5) * dvrdata->dy) * dvrdata->upVec[k];   /* offset bottom-top */

#ifndef DVR_RENDER_ORTHOGONAL
            /* perspective projection: vector from camera pos to near plane */
            dvrdata->Ray[dvrdata->Nray].dir[k] = dvrdata->Ray[dvrdata->Nray].pos[k] - dvrdata->camPos[k];
#else
            /* orthogonal projection: vector parallel to viewing vector */
            dvrdata->Ray[dvrdata->Nray].dir[k] = dvrdata->dir[k];
#endif
            /* ray dir norm */
            dd += dvrdata->Ray[dvrdata->Nray].dir[k] * dvrdata->Ray[dvrdata->Nray].dir[k];
          }

        /* normalize ray direction vector */
        for(k = 0; k < 3; k++)
          dvrdata->Ray[dvrdata->Nray].dir[k] /= sqrt(dd);

        dvrdata->Ray[dvrdata->Nray].len = 0;

        /* scalar product of dir vectors */
        double cos_angle = (dvrdata->dir[0] * dvrdata->Ray[dvrdata->Nray].dir[0] + dvrdata->dir[1] * dvrdata->Ray[dvrdata->Nray].dir[1] + dvrdata->dir[2] * dvrdata->Ray[dvrdata->Nray].dir[2]);

        /* distance from near to far plane */
        dvrdata->Ray[dvrdata->Nray].target_len = (dvrdata->far - dvrdata->near) / cos_angle;

/*
        mpi_printf("CHECK SETUP: %d %d: cos=%g len=%g pos=(%g|%g|%g) dir=(%g|%g|%g)\n",
                x, y, cos_angle, dvrdata->Ray[dvrdata->Nray].target_len,
                dvrdata->Ray[dvrdata->Nray].pos[0], dvrdata->Ray[dvrdata->Nray].pos[1], dvrdata->Ray[dvrdata->Nray].pos[2],
                dvrdata->Ray[dvrdata->Nray].dir[0], dvrdata->Ray[dvrdata->Nray].dir[1], dvrdata->Ray[dvrdata->Nray].dir[2]);
*/
        dvrdata->Ray[dvrdata->Nray].pixel = y * dvrdata->nx + x;
        dvrdata->Ray[dvrdata->Nray].index = -1;
        dvrdata->Ray[dvrdata->Nray].prev = -1;

        /* check local ray */
        peanokey key = position_to_peanokey(dvrdata->Ray[dvrdata->Nray].pos);
        int no = peanokey_to_topnode(key);
        dvrdata->Ray[dvrdata->Nray].task = DomainTask[no];
        if(dvrdata->Ray[dvrdata->Nray].task == ThisTask)
          dvrdata->Nray++;
      }


  mesh_search_data *searchdata;

  searchdata = mymalloc("searchdata", dvrdata->Nray * sizeof(mesh_search_data));

  for(i = 0; i < dvrdata->Nray; i++)
    {
      searchdata[i].Pos[0] = dvrdata->Ray[i].pos[0];
      searchdata[i].Pos[1] = dvrdata->Ray[i].pos[1];
      searchdata[i].Pos[2] = dvrdata->Ray[i].pos[2];
      searchdata[i].u.hsmlguess = 2.0 * get_default_softening_of_particletype(0);
    }


  int tot_rays;
  MPI_Reduce(&dvrdata->Nray, &tot_rays, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  mpi_printf("DVR_RENDER: total number of rays distributed among tasks = %d\n", tot_rays);

  find_nearest_meshpoint_global(searchdata, dvrdata->Nray, 1, 1);
  mpi_distribute_items_from_search(searchdata, dvrdata->Ray, &dvrdata->Nray, &dvrdata->MaxNray, sizeof(dvr_ray_data), TAG_DENS_A, offsetof(dvr_ray_data, task), offsetof(dvr_ray_data, index));
  myfree(searchdata);
  return 0;
}

static int dvr_run(dvr_cam_data * dvrdata)
{
  int i, left_this_task, rays_left, field;
  double t0, t1, t2, t1_tot, t2_tot;
  /* pixel arrays to collect ray data */
  double *localRed = NULL, *Red = NULL;
  double *localGreen = NULL, *Green = NULL;
  double *localBlue = NULL, *Blue = NULL;
  double *localOpacity = NULL, *Opacity = NULL;


  MPI_Barrier(MPI_COMM_WORLD);

  // double check the positions are correct after exchanging
  for(i = 0; i < dvrdata->Nray; i++)
    assert_contains(&Mesh, dvrdata->Ray[i].index, dvrdata->Ray[i].pos);

  t1_tot = 0;
  t2_tot = 0;
  i = 0;
  do
    {
      t0 = second();
      left_this_task = dvr_advance_rays_for_one_cell(dvrdata);
      t1 = second();
      dvr_exchange_rays(dvrdata);
      t2 = second();
      ++i;

      t1_tot += t1 - t0;
      t2_tot += t2 - t1;
      MPI_Allreduce(&left_this_task, &rays_left, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      mpi_printf("DVR_RENDER: iteration %d rays left: %d  t1-t0=%g t2-t1=%g\n", i, rays_left, t1 - t0, t2 - t1);
    }
  while(rays_left);

  mpi_printf("DVR_RENDER: t1_tot=%g t2_tot=%g\n", t1_tot, t2_tot);


  for(field = 0; field < DVR_NUM_FIELDS; field++)
    {
      /* need to collect integrated rays, do this collectively for simplicity */
      localRed = mymalloc("localRed", dvrdata->npixels * sizeof(double));
      localGreen = mymalloc("localGreen", dvrdata->npixels * sizeof(double));
      localBlue = mymalloc("localBlue", dvrdata->npixels * sizeof(double));
      localOpacity = mymalloc("localOpacity", dvrdata->npixels * sizeof(double));

      memset(localRed, 0, dvrdata->npixels * sizeof(double));
      memset(localGreen, 0, dvrdata->npixels * sizeof(double));
      memset(localBlue, 0, dvrdata->npixels * sizeof(double));
      memset(localOpacity, 0, dvrdata->npixels * sizeof(double));

      /* contruct local values; if final ray segment is not on task entry will be zero */
      for(i = 0; i < dvrdata->Nray; i++)
        {
          localRed[dvrdata->Ray[i].pixel] = dvrdata->Ray[i].Red[field];
          localGreen[dvrdata->Ray[i].pixel] = dvrdata->Ray[i].Green[field];
          localBlue[dvrdata->Ray[i].pixel] = dvrdata->Ray[i].Blue[field];
          localOpacity[dvrdata->Ray[i].pixel] = dvrdata->Ray[i].Opacity[field];
        }

      if(ThisTask == 0)
        {
          Red = mymalloc("Red", dvrdata->npixels * sizeof(double));
          Green = mymalloc("Green", dvrdata->npixels * sizeof(double));
          Blue = mymalloc("Blue", dvrdata->npixels * sizeof(double));
          Opacity = mymalloc("Opacity", dvrdata->npixels * sizeof(double));
        }


      MPI_Reduce(localRed, Red, dvrdata->npixels, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(localGreen, Green, dvrdata->npixels, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(localBlue, Blue, dvrdata->npixels, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(localOpacity, Opacity, dvrdata->npixels, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

      mpi_printf("DVR_RENDER: writing files for field %d...\n", field);
      if(ThisTask == 0)
        {
          dvr_write_file("Red", dvrdata, Red, field);
          dvr_write_file("Green", dvrdata, Green, field);
          dvr_write_file("Blue", dvrdata, Blue, field);
          dvr_write_file("Opacity", dvrdata, Opacity, field);

          myfree(Opacity);
          myfree(Blue);
          myfree(Green);
          myfree(Red);
        }

      mpi_printf("DVR_RENDER: files written.\n");
      myfree(localOpacity);
      myfree(localBlue);
      myfree(localGreen);
      myfree(localRed);
    }
  return 0;
}

void dvr_write_file(char *name, dvr_cam_data * dvrdata, double *data, int field)
{
  FILE *fd;
  char fname[MAXLEN_PATH], ext[MAXLEN_PATH];

  float *float_data = mymalloc("float_data", dvrdata->npixels * sizeof(float));
  int i;
  for(i = 0; i < dvrdata->npixels; i++)
    float_data[i] = (float) data[i];

#if (DVR_RENDER==0)
  sprintf(ext, "%03d_%04d", All.SnapshotFileCount, All.DvrFrameFileNum);
#endif
#if (DVR_RENDER==1)
  sprintf(ext, "%04d", All.DvrFrameFileNum);
#endif

  mkdir(All.DvrOutputDir, 02755);

  if(field == DVR_RENDER_TEMP)
    sprintf(fname, "%s/dvr_render_temp_%s_%s", All.DvrOutputDir, name, ext);
  if(field == DVR_RENDER_MET)
    sprintf(fname, "%s/dvr_render_met_%s_%s", All.DvrOutputDir, name, ext);
  if(field == DVR_RENDER_RHO)
    sprintf(fname, "%s/dvr_render_rho_%s_%s", All.DvrOutputDir, name, ext);
  if(field == DVR_RENDER_XRAY)
    sprintf(fname, "%s/dvr_render_xray_%s_%s", All.DvrOutputDir, name, ext);
  if(field >= DVR_RENDER_ELEMENT_0 && field <= DVR_RENDER_ELEMENT_8)
    sprintf(fname, "%s/dvr_render_element_%d_%s_%s", All.DvrOutputDir, field - DVR_RENDER_ELEMENT_0, name, ext);
#ifdef GFM_DUST
  if(field == DVR_RENDER_DUST_RHO)
    sprintf(fname, "%s/dvr_render_dust_rho_%s_%s", All.DvrOutputDir, name, ext);
  if(field == DVR_RENDER_DUST_MET)
    sprintf(fname, "%s/dvr_render_dust_met_%s_%s", All.DvrOutputDir, name, ext);
#endif

  fd = fopen(fname, "w");
  fwrite(&(dvrdata->nx), sizeof(int), 1, fd);
  fwrite(&(dvrdata->ny), sizeof(int), 1, fd);
  fwrite(&(dvrdata->width), sizeof(double), 1, fd);
  fwrite(&(dvrdata->height), sizeof(double), 1, fd);
  fwrite(&(dvrdata->near), sizeof(double), 1, fd);
  fwrite(&(dvrdata->far), sizeof(double), 1, fd);
  fwrite(float_data, sizeof(float), dvrdata->npixels, fd);
  fclose(fd);

  myfree(float_data);
}

static int dvr_deinit(dvr_cam_data * dvrdata)
{
  myfree(dvrdata->Ray);

  return 0;
}

static void dvr_exchange_rays(dvr_cam_data * dvrdata)
{
  int i;
  for(i = 0; i < dvrdata->Nray; i++)
    {
      if(dvrdata->Ray[i].task != ThisTask)
        dvrdata->Ray[i].prev = -1;
    }
  mpi_distribute_items_to_tasks(dvrdata->Ray, offsetof(dvr_ray_data, task), &dvrdata->Nray, &dvrdata->MaxNray, sizeof(*dvrdata->Ray), TAG_DENS_A);
}


/* F2B composite render integral */
static void dvr_update_render_integral(int cell, double len, dvr_ray_data * ray, dvr_cam_data * dvrdata, int field)
{
  color_data color;
  double rgba_voxel_data[2];

  if(ray->Opacity[field] >= 1.0)
    return;

#ifdef DVR_RENDER_SMOOTH
  smooth_cell(cell, ray->pos, field, rgba_voxel_data);
#else
  renderValueSphP(cell, field, rgba_voxel_data);
#endif

  /* evaluate transfer function */
  transferFunctions(rgba_voxel_data, len, field, &color);

  /* front-to-back recursion step: color -> local voxel value; dvrdata -> accumulated data */
  ray->Red[field] += color.Red * color.Alpha * (1.0 - ray->Opacity[field]);
  ray->Green[field] += color.Green * color.Alpha * (1.0 - ray->Opacity[field]);
  ray->Blue[field] += color.Blue * color.Alpha * (1.0 - ray->Opacity[field]);
  ray->Opacity[field] += color.Alpha * (1.0 - ray->Opacity[field]);
}

static int dvr_advance_rays_for_one_cell(dvr_cam_data * dvrdata)
{
  int i, j, nleft = 0;
  int prevtask = -1, nexttask = -1;
  double len;
  int cell, next_edge;
  int field;

  for(i = 0; i < dvrdata->Nray; i++)
    {
      if(dvrdata->Ray[i].len >= dvrdata->Ray[i].target_len)
        continue;

      /* this is the index of the cell in which we currently are */
      cell = dvrdata->Ray[i].index;

      next_edge = find_next_voronoi_cell(&Mesh, cell, dvrdata->Ray[i].pos, dvrdata->Ray[i].dir, dvrdata->Ray[i].prev, &len);

      if((dvrdata->Ray[i].len + len > dvrdata->Ray[i].target_len))
        {
          /*
             ray termination is reached:
             -next cell is past the target length
             -opacity >= threshold (note that we do front-to-back rendering)
           */
          len = dvrdata->Ray[i].target_len - dvrdata->Ray[i].len;
          dvrdata->Ray[i].len = dvrdata->Ray[i].target_len;
          dvrdata->Ray[i].prev = -1;
        }
      else
        {
          /* we will enter next cell. update cell info */
          prevtask = dvrdata->Ray[i].task, nexttask = DC[next_edge].task;
          dvrdata->Ray[i].task = DC[next_edge].task;
          dvrdata->Ray[i].prev = dvrdata->Ray[i].index;
          dvrdata->Ray[i].index = DC[next_edge].index;
          dvrdata->Ray[i].len += len;
          ++nleft;
        }

      /* advance ray */
      for(j = 0; j < 3; j++)
        dvrdata->Ray[i].pos[j] += len * dvrdata->Ray[i].dir[j];

#ifdef DVR_STAY_IN_BOX
      /* do not leave simulation box */
      if((dvrdata->Ray[i].pos[0] < 0.0) || (dvrdata->Ray[i].pos[0] > All.BoxSize) ||
         (dvrdata->Ray[i].pos[1] < 0.0) || (dvrdata->Ray[i].pos[1] > All.BoxSize) || (dvrdata->Ray[i].pos[2] < 0.0) || (dvrdata->Ray[i].pos[2] > All.BoxSize))
        dvrdata->Ray[i].len = 1.01 * dvrdata->Ray[i].target_len;
#endif

/*
      if (dvrdata->Ray[i].pixel==100)
        {
          printf("CHECK ThisTask=%d pos=(%g|%g|%g) len=%g/%g tasks:%d->%d opacity=%g RGB=(%g|%g|%g) rho=%g\n",
                  ThisTask,
                  dvrdata->Ray[i].pos[0], dvrdata->Ray[i].pos[1], dvrdata->Ray[i].pos[2],
                  dvrdata->Ray[i].len, dvrdata->Ray[i].target_len, prevtask, nexttask,
                  dvrdata->Ray[i].Opacity, dvrdata->Ray[i].Red, dvrdata->Ray[i].Green, dvrdata->Ray[i].Blue, SphP[cell].Density);
          fflush(stdout);
        }
*/
      for(field = 0; field < DVR_NUM_FIELDS; field++)
        dvr_update_render_integral(cell, len, &(dvrdata->Ray[i]), dvrdata, field);
      //TODO: periodic wrapping
    }
  return nleft;
}
#endif
