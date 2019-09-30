/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dvr_render/dvr_render.h
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

#ifndef DVR_RENDER_H
#define DVR_RENDER_H

#include "arepoconfig.h"
#include "../allvars.h"

#ifdef DVR_RENDER

#define DVR_RENDER_TEMP 0
#define DVR_RENDER_MET  1
#define DVR_RENDER_RHO  2
#define DVR_RENDER_XRAY 3
#define DVR_RENDER_ELEMENT_0 4
#define DVR_RENDER_ELEMENT_1 5
#define DVR_RENDER_ELEMENT_2 6
#define DVR_RENDER_ELEMENT_3 7
#define DVR_RENDER_ELEMENT_4 8
#define DVR_RENDER_ELEMENT_5 9
#define DVR_RENDER_ELEMENT_6 10
#define DVR_RENDER_ELEMENT_7 11
#define DVR_RENDER_ELEMENT_8 12

#ifndef DVR_NUM_FIELDS
#define DVR_NUM_FIELDS  13
#endif

#ifdef GFM_DUST
#undef DVR_NUM_FIELDS
#define DVR_RENDER_DUST_RHO 13
#define DVR_RENDER_DUST_MET 14
#define DVR_NUM_FIELDS 15
#endif

#if (DVR_NUM_FIELDS < 3)
#error "DVR_NUM_FIELDS < 3 not allowed"
#endif

#define DVR_NUM_TRANS 200

#define DVR_FRUST_BUFF_ZONE 0.5
#define DVR_TESS_ALL

#define DVR_XRAY_TEMP 1e5

/** Struct holding the data for the rays used for making projections. */
typedef struct
{
  /** Index of the primary Voronoi cell in which the ray is located. */
  int index;
  /** Index of the previous cell the particle traversed. */
  int prev;
  /** Task of the primary Voronoi cell in which the ray is located. */
  int task;
  /** The index in the pixel array this ray is associated with. */
  int pixel;
  /** The length the ray has traveled. */
  double len;
  /** The total length the ray will travel when done. */
  double target_len;
  /** The current location of the ray. */
  MyDouble pos[3];
  /** The direction the ray is travelling in. */
  double dir[3];

  float Red[DVR_NUM_FIELDS], Green[DVR_NUM_FIELDS], Blue[DVR_NUM_FIELDS];
  float Opacity[DVR_NUM_FIELDS];
} dvr_ray_data;

typedef struct
{
  /* cam data */
  double camPos[3];
  double camTo[3];
  double dir[3];
  double upVec[3];
  double normalVec[3];
  double width, height;
  double near, far;

  /* image data */
  int nx, ny, npixels, pflag;
  double dx, dy;

  /* ray data */
  int Nray, MaxNray;
  dvr_ray_data *Ray;
} dvr_cam_data;

typedef struct
{
  double Red, Green, Blue, Alpha;
} color_data;

typedef struct
{
  double Val;
  double Red, Green, Blue;
} transRGB_data;

int dvr_render_loop(char *fname, int pixelx, int pixely, int mode);
void dvr_load_transfer_function(char *fname, int field);
#endif

#endif /* DVR_RENDER_PROJ_H */
