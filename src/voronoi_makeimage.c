/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/voronoi_makeimage.c
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

#include <sys/types.h>
#include <sys/stat.h>


#include "allvars.h"
#include "proto.h"
#include "voronoi.h"


#ifndef GFM_N_CHEM_ELEMENTS
#define GFM_N_CHEM_ELEMENTS 0
#endif

/** Extracts the limits and resolutions for the image slice and
    projection routines. These are largely the same for both
    projection and slicing, so this means less duplication, and that
    the slicing routine gets the same centering options as the
    projection ones. */
int get_image_limits(int argc, char **argv, int RestartFlag,
                     int *pixels_x, int *pixels_y, int *xaxis, int *yaxis, int *zaxis, double *xmin, double *xmax, double *ymin, double *ymax, double *zmin, double *zmax, int *weight_flag)
{
#ifdef TGSET
  if(Argc != 13 && Argc != 16)
    {
      mpi_printf
        ("\nwrong parameters found: Expecting:\n Call with: <ParameterFile> %d <SnapNum> <pixelsX> <pixelsY> <axisX> <axisY> <axisZ> <center_flag> <weight_flag>\n <width_flag> <width> <center_x> <center_y> <center_z> (The last three only needed for center_flag = 1)\n\n",
         RestartFlag);
      return (0);
    }
#else
  if(!(((RestartFlag == 4) && (argc == 14)) || ((RestartFlag == 5) && ((argc == 15) || (argc == 16)))
#ifdef SUBBOX_SNAPSHOTS
       || ((RestartFlag == 5) && ((argc == 16) || (argc == 17)))
#endif
       || ((RestartFlag == 8) && (argc == 15))))
    {
      mpi_printf
        ("\nwrong parameters found: Expecting:\n Call with: <ParameterFile> %d <SnapNum> <pixelsX> <pixelsY> <axisX> <axisY> <axisZ>  <xmin> <xmax> <ymin> <ymax> %s\n",
         RestartFlag, RestartFlag == 4 ? "<zval>" : "<zmin> <zmax>");
      if(RestartFlag == 5)
        mpi_printf("Optionally also with <weight_flag> as last argument (default is 0)\n\n");
      return (0);
    }
#endif

  *pixels_x = atoi(argv[4]);
  *pixels_y = atoi(argv[5]);

  *xaxis = atoi(argv[6]);
  *yaxis = atoi(argv[7]);
  *zaxis = atoi(argv[8]);

#ifdef TGSET
  tgset_get_image_limits(argv, xaxis, yaxis, zaxis, xmin, xmax, ymin, ymax, zmin, zmax, weight_flag);
#else
  *xmin = atof(argv[9]);
  *xmax = atof(argv[10]);
  *ymin = atof(argv[11]);
  *ymax = atof(argv[12]);
  *zmin = atof(argv[13]);

  if(RestartFlag == 5 || RestartFlag == 8)
    {
#ifndef SUBBOX_SNAPSHOTS
      *zmax = atof(argv[14]);
      if(Argc == 16)
        {
          *weight_flag = atoi(Argv[15]);
          mpi_printf("VORONOI_MAKEIMAGE: weight_flag = %d\n", *weight_flag);
        }
#else
      *zmax = atof(argv[14]);
      if(Argc == 17)
        {
          *weight_flag = atoi(Argv[16]);
          mpi_printf("VORONOI_MAKEIMAGE: weight_flag = %d\n", *weight_flag);
        }
#endif
    }
#endif

  return -1;
}

#ifdef VORONOI

#if defined(VORONOI_IMAGES_FOREACHSNAPSHOT) || defined(VORONOI_FREQUENT_IMAGES)

void make_voronoi_image(int num)
{
#if !defined(TWODIMS) && !defined (ONEDIMS)

  make_3d_voronoi_projected_image(num, 1, All.PicXpixels, All.PicYpixels, All.PicXaxis, All.PicYaxis, All.PicZaxis, All.PicXmin, All.PicXmax, All.PicYmin, All.PicYmax, All.PicZmin, All.PicZmax, 0);

#ifdef VORONOI_NOGRADS
  make_3d_voronoi_projected_image(num, 0, All.PicXpixels, All.PicYpixels, All.PicXaxis, All.PicYaxis, All.PicZaxis, All.PicXmin, All.PicXmax, All.PicYmin, All.PicYmax, All.PicZmin, All.PicZmax, 0);
#endif
#endif
}

void make_voronoi_image_slice(int num)
{
#if !defined(TWODIMS) && !defined (ONEDIMS)
  make_3d_voronoi_slice_image(num, 1, All.PicXpixels, All.PicYpixels, All.PicXaxis, All.PicYaxis, All.PicZaxis, All.PicXmin, All.PicXmax, All.PicYmin, All.PicYmax, All.PicZmin);
#ifdef VORONOI_NOGRADS
  make_3d_voronoi_slice_image(num, 0, All.PicXpixels, All.PicYpixels, All.PicXaxis, All.PicYaxis, All.PicZaxis, All.PicXmin, All.PicXmax, All.PicYmin, All.PicYmax, All.PicZmin);
#endif
#endif
}

#endif

#if !defined(TWODIMS) && !defined(ONEDIMS)      /* will only be compiled in 3D case */

#define INSIDE_EPS   1.0e-6
#define MAX_COUNT_MOVES 500000


extern const int edge_start[6];
extern const int edge_end[6];
extern const int edge_opposite[6];
extern const int edge_nexttetra[6];


/** Tests the line segment between ppstart and ppend for intersections
    with the tetrahedron tt. Count is set to the number of
    intersections found. face, edge, and corner is set to the
    intersecting elements. orientations is a int[6] array with the
    orientations of the points. (If the line segment intersects more
    than one face, only one of them will be returned. What do we do in
    that case?) */
void voronoi_probe_intersections(tessellation * T, int tt, point * pstart, point * pend, int *count, int *f, int *edge, int *corner, int *orientations)
{
  tetra *DT = T->DT;
  point *DP = T->DP;

  int o1, o2, o3, o4, o5, o6;

  tetra *t = &DT[tt];

  point *p0 = &DP[t->p[0]];
  point *p1 = &DP[t->p[1]];
  point *p2 = &DP[t->p[2]];
  point *p3 = &DP[t->p[3]];

  o1 = Orient3d(p0, p1, pend, pstart);
  o2 = Orient3d(p1, p2, pend, pstart);
  o3 = Orient3d(p2, p0, pend, pstart);
  o4 = Orient3d(p1, p3, pend, pstart);
  o5 = Orient3d(p3, p2, pend, pstart);
  o6 = Orient3d(p3, p0, pend, pstart);

  /*
     o1 = Orient3d_Exact(p0, p1, pend, pstart);
     o2 = Orient3d_Exact(p1, p2, pend, pstart);
     o3 = Orient3d_Exact(p2, p0, pend, pstart);
     o4 = Orient3d_Exact(p1, p3, pend, pstart);
     o5 = Orient3d_Exact(p3, p2, pend, pstart);
     o6 = Orient3d_Exact(p3, p0, pend, pstart);
   */

  *f = -1;
  *edge = -1;
  *corner = -1;
  *count = 0;

  if(o4 == 1 && o5 == 1 && o2 == -1)
    {
      *f = 0;
      *count = *count + 1;
    }
  if(o3 == -1 && o5 == -1 && o6 == 1)
    {
      *f = 1;
      *count = *count + 1;
    }
  if(o4 == -1 && o1 == -1 && o6 == -1)
    {
      *f = 2;
      *count = *count + 1;
    }
  if(o1 == 1 && o2 == 1 && o3 == 1)
    {
      *f = 3;
      *count = *count + 1;
    }


  /* edge between face 0 and face 1 (is edge between corners 2 & 3) */
  if((o4 == 1 && o5 == 0 && o2 == -1 && o3 == -1 && o5 == 0 && o6 == 1) ||
     (o4 == 1 && o5 == 0 && o2 == -1 && o3 == 0 && o5 == 0 && o6 == 0) || (o4 == 0 && o5 == 0 && o2 == 0 && o3 == -1 && o5 == 0 && o6 == 1))
    {
      *edge = 5;
      *count = *count + 1;
    }

  /* edge between face 0 and face 2 (is edge between corners 1 & 3) */
  if((o4 == 0 && o5 == 1 && o2 == -1 && o4 == 0 && o1 == -1 && o6 == -1) ||
     (o4 == 0 && o5 == 1 && o2 == -1 && o4 == 0 && o1 == 0 && o6 == 0) || (o4 == 0 && o5 == 0 && o2 == 0 && o4 == 0 && o1 == -1 && o6 == -1))
    {
      *edge = 4;
      *count = *count + 1;
    }

  /* edge between face 0 and face 3 (is edge between corners 1 & 2) */
  if((o4 == 1 && o5 == 1 && o2 == 0 && o1 == 1 && o2 == 0 && o3 == 1) ||
     (o4 == 1 && o5 == 1 && o2 == 0 && o1 == 0 && o2 == 0 && o3 == 0) || (o4 == 0 && o5 == 0 && o2 == 0 && o1 == 1 && o2 == 0 && o3 == 1))
    {
      *edge = 3;
      *count = *count + 1;
    }

  /* edge between face 1 and face 2 (is edge between corners 0 & 3) */
  if((o3 == -1 && o5 == -1 && o6 == 0 && o4 == -1 && o1 == -1 && o6 == 0) ||
     (o3 == -1 && o5 == -1 && o6 == 0 && o4 == 0 && o1 == 0 && o6 == 0) || (o3 == 0 && o5 == 0 && o6 == 0 && o4 == -1 && o1 == -1 && o6 == 0))
    {
      *edge = 2;
      *count = *count + 1;
    }

  /* edge between face 1 and face 3 (is edge between corners 0 & 2) */
  if((o3 == 0 && o5 == -1 && o6 == 1 && o1 == 1 && o2 == 1 && o3 == 0) ||
     (o3 == 0 && o5 == -1 && o6 == 1 && o1 == 0 && o2 == 0 && o3 == 0) || (o3 == 0 && o5 == 0 && o6 == 0 && o1 == 1 && o2 == 1 && o3 == 0))
    {
      *edge = 1;
      *count = *count + 1;
    }

  /* edge between face 2 and face 3 (is edge between corners 0 & 1) */
  if((o4 == -1 && o1 == 0 && o6 == -1 && o1 == 0 && o2 == 1 && o3 == 1) ||
     (o4 == -1 && o1 == 0 && o6 == -1 && o1 == 0 && o2 == 0 && o3 == 0) || (o4 == 0 && o1 == 0 && o6 == 0 && o1 == 0 && o2 == 1 && o3 == 1))
    {
      *edge = 0;
      *count = *count + 1;
    }

  /* corner 0 (true if all adjacent edges are cut */

  if(o3 == 0 && o5 == -1 && o6 == 0 && o4 == -1 && o1 == 0 && o2 == 1)
    {
      *corner = 0;
      *count = *count + 1;
    }

  /* corner 1 (true if all adjacent edges are cut */
  if(o4 == 0 && o1 == 0 && o6 == -1 && o2 == 0 && o3 == 1 && o5 == 1)
    {
      *corner = 1;
      *count = *count + 1;
    }

  /* corner 2 (true if all adjacent edges are cut */
  if(o3 == 0 && o5 == 0 && o6 == 1 && o1 == 1 && o2 == 0 && o4 == 1)
    {
      *corner = 2;
      *count = *count + 1;
    }

  /* corner 3 (true if all adjacent edges are cut */
  if(o5 == 0 && o6 == 0 && o4 == 0 && o1 == -1 && o2 == -1 && o3 == -1)
    {
      *corner = 3;
      *count = *count + 1;
    }

  orientations[0] = o1;
  orientations[1] = o2;
  orientations[2] = o3;
  orientations[3] = o4;
  orientations[4] = o5;
  orientations[5] = o6;
}


/** Tests whether the line defined by [pstart, pend] has an
  intersection with the tetrahedron tt. IF yes, 1 is returned and the
  coordinates of the point ppexit are updated with the coordinates of
  the exit point.  If no (or if it intersects a corner?) 0 is returned
  and ppexit is updated to the index of the corner point. If the tetra
  contains infinity points and the segment did not cross a valid
  (finity) face, 0 is returned and ppexit is unset. previous_tetra
  indicates a previously visited tetrahedron in case a ray is followed
  from pstart to pend. */
/* OLD? pnew is either the point onto the face entering the next tetra
   given by nexttetra, or equal to p if we have arrived */
int image_get_next_tetra(tessellation * T,
                         /// The DT index of the tetrahedron to test
                         int tt,
                         /// The DP index of the starting point
                         point * pstart,
                         /// The DP index of the end point
                         point * pend,
                         /// [out] The DT index of the next tetrahedron
                         int *nexttetra,
                         /// [out] DP index of the point in which the exit coordinates will be placed
                         point * pexit, int *previous_tetra)
{
  tetra *DT = T->DT;
  point *DP = T->DP;
  tetra *t = &DT[tt];

  double w, u, v, s = 0;
  double denom;
  int status;

  //int pp0, pp1, pp2, pp3;
  point *p0, *p1, *p2, *p3;

  /*
     pp0 = t->p[0];
     pp1 = t->p[1];
     pp2 = t->p[2];
     pp3 = t->p[3];
   */

  p0 = &DP[t->p[0]];
  p1 = &DP[t->p[1]];
  p2 = &DP[t->p[2]];
  p3 = &DP[t->p[3]];

  // This is not a valid end condition! We could go through the valid face!
  // XXX not true anymore, is it?
#ifndef SUNRISE
  if(isInfinity(p0) || isInfinity(p1) || isInfinity(p2) || isInfinity(p3))
    {
      printf("we are in a tetraeder with in infinity point. tetra=%d\n", tt);
      terminate("tetraeder with in infinity point");
    }
#endif

  int face, edge, corner, count, orientations[6];

  voronoi_probe_intersections(T, tt, pstart, pend, &count, &face, &edge, &corner, orientations);



  if(count != 1)
    {
      if(count == 0)
        // no intersection found, that means (I think) we're in a
        // tetra with infinity points and did not cross a valid face
        return 0;

      printf("count=%d\n", count);

      printf("orientations: %d %d %d %d %d %d\n", orientations[0], orientations[1], orientations[2], orientations[3], orientations[4], orientations[5]);
      printf("tetra: %d   adjacent ones %d %d %d %d\n", tt, t->t[0], t->t[1], t->t[2], t->t[3]);
      printf("points: %ld %ld %ld %ld\n", p0 - DP, p1 - DP, p2 - DP, p3 - DP);
      printf("face-cuts:  %d\n", face);
      printf("edge-cuts:  %d\n", edge);
      printf("corner-cuts: %d\n", corner);


      printf("\n");

      printf("p0:     %g %g %g\n", p0->x - pstart->x, p0->y - pstart->y, p0->z - pstart->z);
      printf("p1:     %g %g %g\n", p1->x - pstart->x, p1->y - pstart->y, p1->z - pstart->z);
      printf("p2:     %g %g %g\n", p2->x - pstart->x, p2->y - pstart->y, p2->z - pstart->z);
      printf("p3:     %g %g %g\n", p3->x - pstart->x, p3->y - pstart->y, p3->z - pstart->z);
      printf("\n");
      printf("pstart: %g %g %g\n", pstart->x, pstart->y, pstart->z);
      printf("pend:   %g %g %g\n", pend->x, pend->y, pend->z);

      terminate("we have not found one exit intersection with the tetrahedron");
    }

#ifndef OPTIMIZE_MEMORY_USAGE
  double ax = p1->xx - p0->xx;
  double ay = p1->yy - p0->yy;
  double az = p1->zz - p0->zz;

  double bx = p2->xx - p0->xx;
  double by = p2->yy - p0->yy;
  double bz = p2->zz - p0->zz;

  double cx = p3->xx - p0->xx;
  double cy = p3->yy - p0->yy;
  double cz = p3->zz - p0->zz;

  double qx = pend->xx - p0->xx;
  double qy = pend->yy - p0->yy;
  double qz = pend->zz - p0->zz;
#else
  double ax, ay, az, bx, by, bz, cx, cy, cz, qx, qy, qz;
  IntegerMapType pA_ixyz[3], pB_ixyz[3];
  double pA_xyz[3], pB_xyz[3];

  get_integers_for_point(p0, pA_ixyz, pA_xyz);

  get_integers_for_point(p1, pB_ixyz, pB_xyz);
  ax = pB_xyz[0] - pA_xyz[0];
  ay = pB_xyz[1] - pA_xyz[1];
  az = pB_xyz[2] - pA_xyz[2];

  get_integers_for_point(p2, pB_ixyz, pB_xyz);
  bx = pB_xyz[0] - pA_xyz[0];
  by = pB_xyz[1] - pA_xyz[1];
  bz = pB_xyz[2] - pA_xyz[2];

  get_integers_for_point(p3, pB_ixyz, pB_xyz);
  cx = pB_xyz[0] - pA_xyz[0];
  cy = pB_xyz[1] - pA_xyz[1];
  cz = pB_xyz[2] - pA_xyz[2];

  get_integers_for_point(pend, pB_ixyz, pB_xyz);
  qx = pB_xyz[0] - pA_xyz[0];
  qy = pB_xyz[1] - pA_xyz[1];
  qz = pB_xyz[2] - pA_xyz[2];
#endif

  double mv_data[] = { ax, bx, cx, qx, ay, by, cy, qy, az, bz, cz, qz };
  double xend[3];


  status = solve_linear_equations(mv_data, xend);

  if(status < 0)
    {
      /*
         printf("warning. Potentially inaccurate solution of linear system\n");
         printf("a:     %g %g %g\n", ax, ay, az);
         printf("b:     %g %g %g\n", bx, by, bz);
         printf("b:     %g %g %g\n", cx, cy, cz);
         printf("q:     %g %g %g\n", qx, qy, qz);
         printf("vol=%g test_orient=%d exact_orient=%d\n",
         calculate_tetra_volume(p0, p1, p2, p3),
         test_tetra_orientation(p0, p1, p2, p3),
         Orient3d_Exact(p0, p1, p2, p3));
       */
      //      terminate("status < 0");
    }


#ifndef OPTIMIZE_MEMORY_USAGE
  double qstartx = pstart->xx - p0->xx;
  double qstarty = pstart->yy - p0->yy;
  double qstartz = pstart->zz - p0->zz;
#else
  get_integers_for_point(pstart, pB_ixyz, pB_xyz);

  double qstartx = pB_xyz[0] - pA_xyz[0];
  double qstarty = pB_xyz[1] - pA_xyz[1];
  double qstartz = pB_xyz[2] - pA_xyz[2];
#endif

  double mv_start[] = { ax, bx, cx, qstartx, ay, by, cy, qstarty, az, bz, cz, qstartz };
  double xstart[3];

  status = solve_linear_equations(mv_start, xstart);

  if(status < 0)
    {
      /*
         printf("warning. Potentially inaccurate solution of linear system\n");
         printf("a:     %g %g %g\n", ax, ay, az);
         printf("b:     %g %g %g\n", bx, by, bz);
         printf("b:     %g %g %g\n", cx, cy, cz);
         printf("q:     %g %g %g\n", qx, qy, qz);
         printf("vol=%g test_orient=%d exact_orient=%d\n",
         calculate_tetra_volume(p0, p1, p2, p3),
         test_tetra_orientation(p0, p1, p2, p3),
         Orient3d_Exact(p0, p1, p2, p3));
       */
      //      terminate("status < 0");
    }

  *nexttetra = *previous_tetra;

  if(face >= 0 && face <= 3)
    {
      if(face == 0)             /* cut with face 0 */
        {
          denom = (xend[0] + xend[1] + xend[2] - (xstart[0] + xstart[1] + xstart[2]));
          if(denom != 0)
            s = (1 - (xstart[0] + xstart[1] + xstart[2])) / denom;
          else
            terminate("denom==0");
        }

      if(face == 1)             /* cut with face 1 */
        {
          denom = (xstart[0] - xend[0]);

          if(denom != 0)
            s = xstart[0] / denom;
          else
            terminate("denom==0");
        }


      if(face == 2)             /* cut with face 2 */
        {
          denom = (xstart[1] - xend[1]);

          if(denom != 0)
            s = xstart[1] / denom;
          else
            terminate("denom == 0");
        }


      if(face == 3)             /* cut with face 3 */
        {
          denom = (xstart[2] - xend[2]);

          if(denom != 0)
            s = xstart[2] / denom;
          else
            terminate("denom == 0");
        }
      *nexttetra = t->t[face];
    }

  if(edge >= 0 && face <= 5)
    {
      int o1, o2, o3, o4, o5, o6;

      o1 = orientations[0];
      o2 = orientations[1];
      o3 = orientations[2];
      o4 = orientations[3];
      o5 = orientations[4];
      o6 = orientations[5];

      if(edge == 0)
        {
          if(o4 != 0 || o1 != 0 || o6 != 0)
            {
              denom = (xstart[1] - xend[1]);
              if(denom != 0)
                s = xstart[1] / denom;
              else
                terminate("denom == 0");
            }
          else
            {
              denom = (xstart[2] - xend[2]);
              if(denom != 0)
                s = xstart[2] / denom;
              else
                terminate("denom == 0");
            }
        }
      else if(edge == 1)
        {
          if(o3 != 0 || o5 != 0 || o6 != 0)
            {
              denom = (xstart[0] - xend[0]);
              if(denom != 0)
                s = xstart[0] / denom;
              else
                terminate("denom == 0");
            }
          else
            {
              denom = (xstart[2] - xend[2]);
              if(denom != 0)
                s = xstart[2] / denom;
              else
                {
                  printf("xstart[0]-xend[0]=%g\n", xstart[0] - xend[0]);
                  printf("xstart[1]-xend[1]=%g\n", xstart[1] - xend[1]);
                  printf("xstart[2]-xend[2]=%g\n", xstart[2] - xend[2]);
                  terminate("denom == 0");
                }
            }
        }
      else if(edge == 2)
        {
          if(o3 != 0 || o5 != 0 || o6 != 0)
            {
              denom = (xstart[0] - xend[0]);
              if(denom != 0)
                s = xstart[0] / denom;
              else
                terminate("denom == 0");
            }
          else
            {
              denom = (xstart[1] - xend[1]);
              if(denom != 0)
                s = xstart[1] / denom;
              else
                terminate("denom == 0");
            }
        }
      else if(edge == 3)
        {
          if(o4 != 0 || o5 != 0 || o2 != 0)
            {
              denom = (xend[0] + xend[1] + xend[2] - (xstart[0] + xstart[1] + xstart[2]));
              if(denom != 0)
                s = (1 - (xstart[0] + xstart[1] + xstart[2])) / denom;
              else
                terminate("denom == 0");
            }
          else
            {
              denom = (xstart[2] - xend[2]);
              if(denom != 0)
                s = xstart[2] / denom;
              else
                terminate("denom == 0");
            }
        }
      else if(edge == 4)
        {
          if(o4 != 0 || o5 != 0 || o2 != 0)
            {
              denom = (xend[0] + xend[1] + xend[2] - (xstart[0] + xstart[1] + xstart[2]));
              if(denom != 0)
                s = (1 - (xstart[0] + xstart[1] + xstart[2])) / denom;
              else
                terminate("denom == 0");
            }
          else
            {
              denom = (xstart[1] - xend[1]);
              if(denom != 0)
                s = xstart[1] / denom;
              else
                terminate("denom == 0");
            }
        }
      else if(edge == 5)
        {
          if(o4 != 0 || o5 != 0 || o2 != 0)
            {
              denom = (xend[0] + xend[1] + xend[2] - (xstart[0] + xstart[1] + xstart[2]));
              if(denom != 0)
                s = (1 - (xstart[0] + xstart[1] + xstart[2])) / denom;
              else
                terminate("denom == 0");
            }
          else
            {
              denom = (xstart[0] - xend[0]);
              if(denom != 0)
                s = xstart[0] / denom;
              else
                terminate("denom==0");
            }
        }

      /* now we circle around the edge to find the correct next tetrahedron */

      int i, j, k, l, m, ii, jj, kk, ll, pp, nn, iter;
      tetra *prev, *next;

      i = edge_start[edge];
      j = edge_end[edge];
      k = edge_opposite[edge];
      l = edge_nexttetra[edge];

      iter = 0;

      pp = tt;
      prev = &DT[pp];
      do
        {
          nn = prev->t[l];
          next = &DT[nn];

          //      printf("nn=%d tt=%d *previous_tetra=%d\n", nn, tt, *previous_tetra);

          if(nn != tt && nn != (*previous_tetra))
            {
              int face2, edge2, corner2, count2, orientations2[6];

              voronoi_probe_intersections(T, nn, pstart, pend, &count2, &face2, &edge2, &corner2, orientations2);
              /*
                 printf("count2=%d\n", count2);
                 printf("orientations2: %d %d %d %d %d %d\n", orientations2[0], orientations2[1], orientations2[2], orientations2[3], orientations2[4], orientations2[5]);
               */
              if(count2 > 1)
                terminate("count2 > 1");

              if(count2 == 1)
                {
                  *nexttetra = nn;
                  break;
                }
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
          pp = nn;
          i = ii;
          l = ll;
          j = jj;
          k = kk;

          iter++;

          if(iter > 1000)
            terminate("iter is too large");
        }
      while(nn != tt);

    }


  if(corner >= 0 && corner <= 3)
    {
      *pexit = DP[t->p[corner]];

      if(pow(pexit->x - pstart->x, 2) + pow(pexit->y - pstart->y, 2) + pow(pexit->z - pstart->z, 2) > pow(pend->x - pstart->x, 2) + pow(pend->y - pstart->y, 2) + pow(pend->z - pstart->z, 2)
#ifdef SUNRISE
         || t->p[corner] == DPinfinity
#endif
        )
        {
          *pexit = *pend;
          return 0;
        }

      /* now we look among all tetras that share the point for a suitable next tetrahedron */

      int i;
      int pp = t->p[corner];

      for(i = 0; i < T->Ndt; i++)
        {
          if(DT[i].p[0] == pp || DT[i].p[1] == pp || DT[i].p[2] == pp || DT[i].p[3] == pp)
            if(i != tt && i != (*previous_tetra))
              {
                if(DT[i].t[0] < 0)      /* don't consider deleted ones */
                  continue;

                int face2, edge2, corner2, count2, orientations2[6];

                voronoi_probe_intersections(T, i, pstart, pend, &count2, &face2, &edge2, &corner2, orientations2);
                /*
                   printf("count2=%d\n", count2);
                   printf("orientations2: %d %d %d %d %d %d\n", orientations2[0], orientations2[1], orientations2[2], orientations2[3], orientations2[4], orientations2[5]);
                 */
                if(count2 > 1)
                  terminate("count2 > 1");

                if(count2 == 1)
                  {
                    *nexttetra = i;
                    break;
                  }
              }
        }
    }


  if(face >= 0 || edge >= 0)
    {
      if(s >= 1)
        {
          *pexit = *pend;
          return 0;
        }

      if(status >= 0)
        {
          u = xstart[0] + s * (xend[0] - xstart[0]);
          v = xstart[1] + s * (xend[1] - xstart[1]);
          w = xstart[2] + s * (xend[2] - xstart[2]);

#ifndef OPTIMIZE_MEMORY_USAGE
          pexit->xx = u * ax + v * bx + w * cx + p0->xx;
          pexit->yy = u * ay + v * by + w * cy + p0->yy;
          pexit->zz = u * az + v * bz + w * cz + p0->zz;
          pexit->x = (pexit->xx - 1.0) / ConversionFac + CentralOffsetX;
          pexit->y = (pexit->yy - 1.0) / ConversionFac + CentralOffsetY;
          pexit->z = (pexit->zz - 1.0) / ConversionFac + CentralOffsetZ;
#else
          double pexit_xyz[3];
          pexit_xyz[0] = u * ax + v * bx + w * cx + pA_xyz[0];
          pexit_xyz[1] = u * ay + v * by + w * cy + pA_xyz[1];
          pexit_xyz[2] = u * az + v * bz + w * cz + pA_xyz[2];

          pexit->x = (pexit_xyz[0] - 1.0) / ConversionFac + CentralOffsetX;
          pexit->y = (pexit_xyz[1] - 1.0) / ConversionFac + CentralOffsetY;
          pexit->z = (pexit_xyz[2] - 1.0) / ConversionFac + CentralOffsetZ;
#endif
        }
    }

  if(*nexttetra == *previous_tetra)
    {
      printf("orientations: %d %d %d %d %d %d\n", orientations[0], orientations[1], orientations[2], orientations[3], orientations[4], orientations[5]);

      printf("face-cuts:  %d\n", face);
      printf("edge-cuts:  %d\n", edge);
      printf("corner-cuts: %d\n", corner);


      printf("\n");
      printf("p0:     %g %g %g\n", p0->x - pstart->x, p0->y - pstart->y, p0->z - pstart->z);
      printf("p1:     %g %g %g\n", p1->x - pstart->x, p1->y - pstart->y, p1->z - pstart->z);
      printf("p2:     %g %g %g\n", p2->x - pstart->x, p2->y - pstart->y, p2->z - pstart->z);
      printf("p3:     %g %g %g\n", p3->x - pstart->x, p3->y - pstart->y, p3->z - pstart->z);
      printf("\n");
      printf("pstart: %g %g %g\n", pstart->x, pstart->y, pstart->z);
      printf("pend:   %g %g %g\n", pend->x, pend->y, pend->z);

      terminate("nexttetra == previous_tetra");
    }

  return 1;
}



/* if this is defined, the definition is in voronoi_makeimage_new.c */
#ifndef VORONOI_NEW_IMAGE
void make_3d_voronoi_projected_image(int num, int gradients_flag, int pixels_x, int pixels_y, int xaxis,
                                     int yaxis, int zaxis, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, int weight_flag)
{
  CPU_Step[CPU_MISC] += measure_time();
  tessellation *T = &Mesh;

  char buf[1000], msg[1000];
  float *dens, *denssum, *dp, *temp, *tempsum, *tp, *weight, *weightsum, *wp;
#if defined(COOLING) && !defined(GRACKLE)
  float *szy, *szysum, *sp;
#endif
#ifdef GFM_STELLAR_EVOLUTION
  float *metal, *metalsum, *mp;
#endif
#ifdef GFM_AGN_RADIATION
  float *agnbol, *agnbolsum, *ap;
#endif
#ifdef TRACER_MC
  float *tracernum, *tracernumsum, *trp, *trweight, *trweightsum, *trwp;
#endif
#ifdef CHEM_IMAGE
  float *dust, *xH2, *xHP, *xCO, *dustsum, *h2sum, *hpsum, *cosum, *dst, *dh2, *dhp, *dco;
#endif

  FILE *fd = 0;
  //int pp, ppend, ppold, ppnew, ppstart;
  point p, pold, pnew, pend, pstart;
  int tt0, tt, ttstart, ttrow, next_tetra, previous_tetra;
  double sigma, sigmatemp, sigmaweight;
#if defined(COOLING) && !defined(GRACKLE)
  double sigmaszy;
#endif
#ifdef GFM_STELLAR_EVOLUTION
  double sigmametal;
#endif
#ifdef GFM_AGN_RADIATION
  double sigmaagnbol;
#endif
#ifdef TRACER_MC
  double sigmatracernum, sigmatrweight;
#endif
#ifdef CHEM_IMAGE
  double sigmadust, sigmah2, sigmahp, sigmaco;
#endif
  int i, j, moves, ret, count_moves;

  mpi_printf("Warning: Generating projected image without VORONOI_NEW_IMAGE.\nThe resulting column densities are known to be incorrect.\n");
  fflush(stdout);

  mpi_printf("we start to project the image... gradients_flag=%d\n", gradients_flag);

  T->DTF = mymalloc_movable(&T->DTF, "DTF", T->MaxNdt * sizeof(char));
  tetra *DT = T->DT;

  for(i = 0; i < T->Ndt; i++)
    T->DTF[i] = 0;

  compute_circumcircles(T);

  /*
     xmin += All.BoxSize / 2;
     xmax += All.BoxSize / 2;
     ymin += All.BoxSize / 2;
     ymax += All.BoxSize / 2;
     zmin += All.BoxSize / 2;
     zmax += All.BoxSize / 2;
   */

  // check that we have enough space in DP for the points we
  // need. Note that we do NOT increase Ndp -- this means that this
  // code is not reentrant, though we can pass separate DP arrays to
  // different threads
  while((T->Ndp >= (T->MaxNdp - 5)))
    {
      T->Indi.AllocFacNdp *= ALLOC_INCREASE_FACTOR;
      T->MaxNdp = T->Indi.AllocFacNdp;
#ifdef VERBOSE
      printf("Task=%d: increase memory allocation, MaxNdp=%d Indi.AllocFacNdp=%g\n", ThisTask, T->MaxNdp, T->Indi.AllocFacNdp);
#endif
      T->DP -= 5;
      T->DP = myrealloc_movable(T->DP, (T->MaxNdp + 5) * sizeof(point));
      T->DP += 5;

      if(T->Ndp >= (T->MaxNdp - 5) && NumGas == 0)
        terminate("(Ndp >= (MaxNdp - 5)");
    }

  if(gradients_flag == 1)
    sprintf(buf, "%s/proj_density_field_%03d", All.OutputDir, num);
  else if(gradients_flag == 0)
    sprintf(buf, "%s/proj_density_field_nograds_%03d", All.OutputDir, num);
  else
    terminate("gradients_flag != 0 && gradients_flag != 1");

  if(ThisTask == 0)
    {
      if(!(fd = fopen(buf, "w")))
        {
          sprintf(msg, "can't open file `%s' for writing snapshot.\n", buf);
          terminate(msg);
        }

      my_fwrite(&pixels_x, sizeof(int), 1, fd);
      my_fwrite(&pixels_y, sizeof(int), 1, fd);
    }

  // allocate images and zero them out
  dens = mymalloc("dens", pixels_x * pixels_y * sizeof(float));
  denssum = mymalloc("denssum", pixels_x * pixels_y * sizeof(float));

  temp = mymalloc("temp", pixels_x * pixels_y * sizeof(float));
  tempsum = mymalloc("tempsum", pixels_x * pixels_y * sizeof(float));

  weight = mymalloc("weight", pixels_x * pixels_y * sizeof(float));
  weightsum = mymalloc("weightsum", pixels_x * pixels_y * sizeof(float));

#if defined(COOLING) && !defined(GRACKLE)
  szy = mymalloc("szy", pixels_x * pixels_y * sizeof(float));
  szysum = mymalloc("szysum", pixels_x * pixels_y * sizeof(float));
#endif

#ifdef GFM_STELLAR_EVOLUTION
  metal = mymalloc("metal", pixels_x * pixels_y * sizeof(float));
  metalsum = mymalloc("metalsum", pixels_x * pixels_y * sizeof(float));
#endif

#ifdef GFM_AGN_RADIATION
  agnbol = mymalloc("agnbol", pixels_x * pixels_y * sizeof(float));
  agnbolsum = mymalloc("agnbolsum", pixels_x * pixels_y * sizeof(float));
#endif

#ifdef TRACER_MC
  tracernum = mymalloc("tracernum", pixels_x * pixels_y * sizeof(float));
  tracernumsum = mymalloc("tracernumsum", pixels_x * pixels_y * sizeof(float));

  trweight = mymalloc("trweight", pixels_x * pixels_y * sizeof(float));
  trweightsum = mymalloc("trweightsum", pixels_x * pixels_y * sizeof(float));
#endif

#ifdef CHEM_IMAGE
  dust = mymalloc("dust", pixels_x * pixels_y * sizeof(float));
  dustsum = mymalloc("dustsum", pixels_x * pixels_y * sizeof(float));

  xH2 = mymalloc("xH2", pixels_x * pixels_y * sizeof(float));
  h2sum = mymalloc("h2sum", pixels_x * pixels_y * sizeof(float));

  xHP = mymalloc("xHP", pixels_x * pixels_y * sizeof(float));
  hpsum = mymalloc("hpsum", pixels_x * pixels_y * sizeof(float));

  xCO = mymalloc("xCO", pixels_x * pixels_y * sizeof(float));
  cosum = mymalloc("cosum", pixels_x * pixels_y * sizeof(float));
#endif

  for(i = 0, dp = dens; i < pixels_x; i++)
    for(j = 0; j < pixels_y; j++)
      *dp++ = 0;

  for(i = 0, tp = temp; i < pixels_x; i++)
    for(j = 0; j < pixels_y; j++)
      *tp++ = 0;

  for(i = 0, wp = weight; i < pixels_x; i++)
    for(j = 0; j < pixels_y; j++)
      *wp++ = 0;

#if defined(COOLING) && !defined(GRACKLE)
  for(i = 0, sp = szy; i < pixels_x; i++)
    for(j = 0; j < pixels_y; j++)
      *sp++ = 0;
#endif

#ifdef GFM_STELLAR_EVOLUTION
  for(i = 0, mp = metal; i < pixels_x; i++)
    for(j = 0; j < pixels_y; j++)
      *mp++ = 0;
#endif

#ifdef GFM_AGN_RADIATION
  for(i = 0, ap = agnbol; i < pixels_x; i++)
    for(j = 0; j < pixels_y; j++)
      *ap++ = 0;
#endif

#ifdef TRACER_MC
  for(i = 0, trp = tracernum; i < pixels_x; i++)
    for(j = 0; j < pixels_y; j++)
      *trp++ = 0;

  for(i = 0, trwp = trweight; i < pixels_x; i++)
    for(j = 0; j < pixels_y; j++)
      *trwp++ = 0;
#endif

#ifdef CHEM_IMAGE
  for(i = 0, dst = dust; i < pixels_x; i++)
    for(j = 0; j < pixels_y; j++)
      *dst++ = 0;

  for(i = 0, dh2 = xH2; i < pixels_x; i++)
    for(j = 0; j < pixels_y; j++)
      *dh2++ = 0;

  for(i = 0, dhp = xHP; i < pixels_x; i++)
    for(j = 0; j < pixels_y; j++)
      *dhp++ = 0;

  for(i = 0, dco = xCO; i < pixels_x; i++)
    for(j = 0; j < pixels_y; j++)
      *dco++ = 0;
#endif

  // Get a suitable start tetrahedron?
  ttrow = 0;
  while(DT[ttrow].t[0] < 0 || DT[ttrow].p[0] == DPinfinity || DT[ttrow].p[1] == DPinfinity || DT[ttrow].p[2] == DPinfinity || DT[ttrow].p[3] == DPinfinity)
    ttrow++;

  // Loop over x (columns)
  for(i = 0; i < pixels_x; i++)
    {
      ttstart = ttrow;

      // Loop over y (rows)
      for(j = 0; j < pixels_y; j++)
        {
          // The logic below sets the points p and pend correctly
          // depending on which axes we have chosen.
          if(xaxis == 0 && yaxis == 1)
            {
              p.x = (i + 0.5) / pixels_x * (xmax - xmin) + xmin;
              p.y = (j + 0.5) / pixels_y * (ymax - ymin) + ymin;
              p.z = zmin;
              pend.x = p.x;
              pend.y = p.y;
              pend.z = zmax;
            }
          else if(xaxis == 1 && yaxis == 0)
            {
              p.x = (j + 0.5) / pixels_y * (ymax - ymin) + ymin;
              p.y = (i + 0.5) / pixels_x * (xmax - xmin) + xmin;
              p.z = zmin;
              pend.x = p.x;
              pend.y = p.y;
              pend.z = zmax;
            }
          else if(xaxis == 0 && yaxis == 2)
            {
              p.x = (i + 0.5) / pixels_x * (xmax - xmin) + xmin;
              p.y = zmin;
              p.z = (j + 0.5) / pixels_y * (ymax - ymin) + ymin;
              pend.x = p.x;
              pend.y = zmax;
              pend.z = p.z;
            }
          else if(xaxis == 2 && yaxis == 0)
            {
              p.x = (j + 0.5) / pixels_y * (ymax - ymin) + ymin;
              p.y = zmin;
              p.z = (i + 0.5) / pixels_x * (xmax - xmin) + xmin;
              pend.x = p.x;
              pend.y = zmax;
              pend.z = p.z;
            }
          else if(xaxis == 1 && yaxis == 2)
            {
              p.x = zmin;
              p.y = (i + 0.5) / pixels_x * (xmax - xmin) + xmin;
              p.z = (j + 0.5) / pixels_y * (ymax - ymin) + ymin;
              pend.x = zmax;
              pend.y = p.y;
              pend.z = p.z;
            }
          else if(xaxis == 2 && yaxis == 1)
            {
              p.x = zmin;
              p.y = (j + 0.5) / pixels_y * (ymax - ymin) + ymin;
              p.z = (i + 0.5) / pixels_x * (xmax - xmin) + xmin;
              pend.x = zmax;
              pend.y = p.y;
              pend.z = p.z;
            }
          else
            terminate("invalid combination of axes");
#ifndef OPTIMIZE_MEMORY_USAGE
          set_integers_for_pointer(&p);
          set_integers_for_pointer(&pend);
#endif
          if(ttstart < 0)
            terminate("ttstart < 0");

          if(DT[ttstart].t[0] < 0)
            {
              printf("this start tetra (ttstart=%d) was deleted\n", ttstart);
              terminate("deleted start tetrahedron");
            }

          // Find first tetra that contains point pp
          tt0 = get_tetra(T, &p, &moves, ttstart, &ret, &ret);

          // Update starting tetras in both row and column, if applicable
          if(j == 0)
            ttrow = tt0;

          ttstart = tt0;


          tt = tt0;

          pstart = p;           /* this is the starting point */
          pold = p;
          previous_tetra = -1;

          count_moves = 0;

          /*
             mpi_printf("pixel=%d|%d\n", i, j);
           */

          // Now walk through the tetras along the ray. ppnew is set
          // to the exit point from this/entry point to next tetra
          while((ret = image_get_next_tetra(T, tt, &pstart, &pend, &next_tetra, &pnew, &previous_tetra)))
            {
              //    printf("moves=%d tt=%d\n", count_moves, tt);
#ifndef OPTIMIZE_MEMORY_USAGE
              set_integers_for_pointer(&pold);
              set_integers_for_pointer(&pnew);
#endif
              if(count_moves > MAX_COUNT_MOVES - 10)
                {
                  printf("i/j=(%d|%d)  tetra=%d next=%d   x=%g y=%g z=%g\n", i, j, (int) (tt), (int) (next_tetra), pnew.x, pnew.y, pnew.z);
                }

              count_moves++;

              if(count_moves > MAX_COUNT_MOVES)
                terminate("count_moves > MAX_COUNT_MOVES");

              // calculate the contribution to this pixel for this tetrahedron
              calc_picture_contribution(T, tt, &pold, &pnew, &sigma, &sigmatemp, &sigmaweight, weight_flag, gradients_flag
#if defined(COOLING) && !defined(GRACKLE)
                                        , &sigmaszy
#endif
#ifdef GFM_STELLAR_EVOLUTION
                                        , &sigmametal
#endif
#ifdef GFM_AGN_RADIATION
                                        , &sigmaagnbol
#endif
#ifdef TRACER_MC
                                        , &sigmatracernum, &sigmatrweight
#endif
#ifdef CHEM_IMAGE
                                        , &sigmadust, &sigmah2, &sigmahp, &sigmaco
#endif
                );

              /*
                 printf("old=(%g|%g|%g)   new=(%g|%g|%g)  dens=%g   dens[i * pixels_y + j]=%g \n",
                 pold.x, pold.y, pold.z,
                 pnew.x, pnew.y, pnew.z, sigma, dens[i * pixels_y + j]);
               */

              dens[i * pixels_y + j] += sigma;
              temp[i * pixels_y + j] += sigmatemp;
              weight[i * pixels_y + j] += sigmaweight;
#if defined(COOLING) && !defined(GRACKLE)
              szy[i * pixels_y + j] += sigmaszy * BOLTZMANN * THOMPSON / (ELECTRONMASS * CLIGHT * CLIGHT);
#endif
#ifdef GFM_STELLAR_EVOLUTION
              metal[i * pixels_y + j] += sigmametal;
#endif
#ifdef GFM_AGN_RADIATION
              agnbol[i * pixels_y + j] += sigmaagnbol;
#endif
#ifdef TRACER_MC
              tracernum[i * pixels_y + j] += sigmatracernum;
              trweight[i * pixels_y + j] += sigmatrweight;
#endif
#ifdef CHEM_IMAGE
              dust[i * pixels_y + j] += sigmadust;
              xH2[i * pixels_y + j] += sigmah2;
              xHP[i * pixels_y + j] += sigmahp;
              xCO[i * pixels_y + j] += sigmaco;
#endif
              // update "old" values
              previous_tetra = tt;
              tt = next_tetra;

              pold = pnew;
            }

          calc_picture_contribution(T, tt, &pold, &pnew, &sigma, &sigmatemp, &sigmaweight, weight_flag, gradients_flag
#if defined(COOLING) && !defined(GRACKLE)
                                    , &sigmaszy
#endif
#ifdef GFM_STELLAR_EVOLUTION
                                    , &sigmametal
#endif
#ifdef GFM_AGN_RADIATION
                                    , &sigmaagnbol
#endif
#ifdef TRACER_MC
                                    , &sigmatracernum, &sigmatrweight
#endif
#ifdef CHEM_IMAGE
                                    , &sigmadust, &sigmah2, &sigmahp, &sigmaco
#endif
            );


          dens[i * pixels_y + j] += sigma;
          temp[i * pixels_y + j] += sigmatemp;
          weight[i * pixels_y + j] += sigmaweight;
#if defined(COOLING) && !defined(GRACKLE)
          szy[i * pixels_y + j] += sigmaszy * BOLTZMANN * THOMPSON / (ELECTRONMASS * CLIGHT * CLIGHT);
#endif
#ifdef GFM_STELLAR_EVOLUTION
          metal[i * pixels_y + j] += sigmametal;
#endif
#ifdef GFM_AGN_RADIATION
          agnbol[i * pixels_y + j] += sigmaagnbol;
#endif
#ifdef TRACER_MC
          tracernum[i * pixels_y + j] += sigmatracernum;
          trweight[i * pixels_y + j] += sigmatrweight;
#endif
#ifdef CHEM_IMAGE
          dust[i * pixels_y + j] += sigmadust;
          xH2[i * pixels_y + j] += sigmah2;
          xHP[i * pixels_y + j] += sigmahp;
          xCO[i * pixels_y + j] += sigmaco;
#endif
        }
    }


  MPI_Reduce(dens, denssum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(temp, tempsum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(weight, weightsum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#if defined(COOLING) && !defined(GRACKLE)
  MPI_Reduce(szy, szysum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
#ifdef GFM_STELLAR_EVOLUTION
  MPI_Reduce(metal, metalsum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
#ifdef GFM_AGN_RADIATION
  MPI_Reduce(agnbol, agnbolsum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
#ifdef TRACER_MC
  MPI_Reduce(tracernum, tracernumsum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(trweight, trweightsum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
#ifdef CHEM_IMAGE
  MPI_Reduce(dust, dustsum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(xH2, h2sum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(xHP, hpsum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(xCO, cosum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#endif

  if(ThisTask == 0)
    {
      for(i = 0, wp = weightsum, tp = tempsum, dp = denssum
#if defined(COOLING) && !defined(GRACKLE)
          , sp = szysum
#endif
#ifdef GFM_STELLAR_EVOLUTION
          , mp = metalsum
#endif
#ifdef GFM_AGN_RADIATION
          , ap = agnbolsum
#endif
#ifdef TRACER_MC
          , trp = tracernumsum, trwp = trweightsum
#endif
#ifdef CHEM_IMAGE
          , dst = dustsum, dh2 = h2sum, dhp = hpsum, dco = cosum
#endif
          ; i < pixels_x; i++)

        for(j = 0; j < pixels_y; j++, dp++, tp++, wp++
#if defined(COOLING) && !defined(GRACKLE)
            , sp++
#endif
#ifdef GFM_STELLAR_EVOLUTION
            , mp++
#endif
#ifdef GFM_AGN_RADIATION
            , ap++
#endif
#ifdef TRACER_MC
            , trp++, trwp++
#endif
#ifdef CHEM_IMAGE
            , dst++, dh2++, dhp++, dco++
#endif
          )
          {
            if(*dp > 0)
              {
                *tp /= *dp;
              }

#ifdef CHEM_IMAGE
            if(*dp > 0)
              {
                *dst /= *dp;
                *dh2 /= *dp;
                *dhp /= *dp;
                *dco /= *dp;
              }
#endif

            if(*wp > 0)
              *dp /= *wp;

#ifdef GFM_STELLAR_EVOLUTION
            if(*dp > 0)
              *mp /= *dp;
#endif
#ifdef GFM_AGN_RADIATION
            if(*dp > 0)
              *ap /= *dp;
#endif
#ifdef TRACER_MC
            if(*trwp > 0)
              *trp /= *trwp;
#endif

            if(weight_flag > 0)
              {
                *dp *= (zmax - zmin);
#ifdef TRACER_MC
                *trp *= (zmax - zmin);
#endif
              }
          }

      my_fwrite(denssum, sizeof(float), pixels_x * pixels_y, fd);
      my_fwrite(tempsum, sizeof(float), pixels_x * pixels_y, fd);
#if defined(COOLING) && !defined(GRACKLE)
      my_fwrite(szysum, sizeof(float), pixels_x * pixels_y, fd);
#endif
#ifdef GFM_STELLAR_EVOLUTION
      my_fwrite(metalsum, sizeof(float), pixels_x * pixels_y, fd);
#endif
#ifdef GFM_AGN_RADIATION
      my_fwrite(agnbolsum, sizeof(float), pixels_x * pixels_y, fd);
#endif
#ifdef TRACER_MC
      my_fwrite(tracernumsum, sizeof(float), pixels_x * pixels_y, fd);
#endif
#ifdef CHEM_IMAGE
      my_fwrite(dustsum, sizeof(float), pixels_x * pixels_y, fd);
      my_fwrite(h2sum, sizeof(float), pixels_x * pixels_y, fd);
      my_fwrite(hpsum, sizeof(float), pixels_x * pixels_y, fd);
      my_fwrite(cosum, sizeof(float), pixels_x * pixels_y, fd);
#endif
      fclose(fd);
    }

#ifdef CHEM_IMAGE
  myfree(cosum);
  myfree(xCO);
  myfree(hpsum);
  myfree(xHP);
  myfree(h2sum);
  myfree(xH2);
  myfree(dustsum);
  myfree(dust);
#endif
#ifdef TRACER_MC
  myfree(trweightsum);
  myfree(trweight);
  myfree(tracernumsum);
  myfree(tracernum);
#endif
#ifdef GFM_AGN_RADIATION
  myfree(agnbolsum);
  myfree(agnbol);
#endif
#ifdef GFM_STELLAR_EVOLUTION
  myfree(metalsum);
  myfree(metal);
#endif
#if defined(COOLING) && !defined(GRACKLE)
  myfree(szysum);
  myfree(szy);
#endif
  myfree(weightsum);
  myfree(weight);
  myfree(tempsum);
  myfree(temp);
  myfree(denssum);
  myfree(dens);

  myfree(T->DTF);

  CPU_Step[CPU_MAKEIMAGES] += measure_time();
}
#endif

/* this is the code for simple tetrahedra interpolation */
#ifdef SIMPLE_TETRAHEDRA_PROJECTION
void calc_picture_contribution(tessellation * T, tetra * t, point * p0, point * p1, double *sigma, double *sigmatemp, double *sigmaweight, int weight_flag, int gradients_flag
#if defined(COOLING) && !defined(GRACKLE)
                               , double *sigmaszy
#endif
  )
{
  double rho = 0;
  int k, li;

  for(k = 0; k < 4; k++)
    {
      if(t->p[k]->task == ThisTask)
        {
          li = t->p[k]->index;

          if(li >= NumGas)
            li -= NumGas;

          rho += SphP[li].Density;
        }
    }

  rho /= 4;

  double dx = p0->x - p1->x;
  double dy = p0->y - p1->y;
  double dz = p0->z - p1->z;

  double r = sqrt(dx * dx + dy * dy + dz * dz);

  if(weight_flag > 0)
    {
      *sigma = pow(rho, weight_flag) * r;
      *sigmaweight = pow(rho, weight_flag - 1) * r;
    }
  else
    {
      *sigma = rho * r;
      *sigmaweight = 0;
    }

  *sigmatemp = 0;
#if defined(COOLING) && !defined(GRACKLE)
  *sigmaszy = 0;
#endif
}
#else

/** These 6 permutations of 0,1,2,3 define the intersection tests
    below. The first two numbers are the 6 edges. The final two are
    just the two remaining points that must be tested against. */
const char pairs[6][4] = {
  {0, 1, 2, 3},
  {0, 2, 1, 3},
  {0, 3, 1, 2},
  {1, 2, 0, 3},
  {1, 3, 0, 2},
  {2, 3, 0, 1}
};




/** Checks that the specified position is on the inside of all the
    face planes that make up the specified cell. */
double assert_contains(tessellation * T, int cell, MyDouble p0[3])
{
#ifdef VORONOI_DYNAMIC_UPDATE
  point *DP = T->DP;

  MyDouble cell_p[3];
  cell_p[0] = P[cell].Pos[0];
  cell_p[1] = P[cell].Pos[1];
  cell_p[2] = P[cell].Pos[2];

  // if cell center is across the boundary, wrap it
  periodic_wrap_point(cell_p, p0);

  MyDouble nb_p[3];
  double m[3];
  double c[3];
  double q[3];

  int edge = SphP[cell].first_connection;
  int last_edge = SphP[cell].last_connection;

  double mindist = MAX_DOUBLE_NUMBER;

  while(1)
    {
      const int neighbor = DC[edge].dp_index;
      myassert((DC[edge].task != ThisTask) || (DC[edge].index != cell));

      nb_p[0] = DP[neighbor].x;
      nb_p[1] = DP[neighbor].y;
      nb_p[2] = DP[neighbor].z;
      // if neighbor is across the boundary, wrap it
      periodic_wrap_point(nb_p, p0);


      int i;
      for(i = 0; i < 3; ++i)
        {
          // m is the edge midpoint, which is a point on the face 
          m[i] = 0.5 * (nb_p[i] + cell_p[i]);
          // c is the vector from entry point p0 to the point on the face
          c[i] = m[i] - p0[i];
          /* q is the edge vector to the neighboring cell, which is a
             normal vector of the face plane, pointing outward. */
          q[i] = nb_p[i] - cell_p[i];
        }

      double cdotq = c[0] * q[0] + c[1] * q[1] + c[2] * q[2];
      // distance from point to face, counting POSITIVE INSIDE. (If
      // q and c point in the same direction, point is inside and
      // cdotq is positive.)
      double dist_to_face = cdotq / sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);
      //myassert(dist_to_face > -1e-6);
      if(dist_to_face < mindist)
        mindist = dist_to_face;

      if(edge == last_edge)
        break;

      edge = DC[edge].next;
    }
  return mindist;
#else
  myassert(0);
  return 0;
#endif
}


/* Calculates the intersections between a ray and the internal Voronoi
   faces in a Delaunay tetrahedron. The intersections are returned in
   list (which must be an array of 8 elements), with nlist indicating
   how many intersections there are. (There are always at least 2, the
   start and end points.) */
void calc_delaunay_intersections(tessellation * T, int tt, point * p0, point * p1, intersection_list * list, int *nlist)
{
  tetra *DT = T->DT;
  point *DP = T->DP;
  tetra_center *DTC = T->DTC;
  tetra *t = &DT[tt];
  tetra_center *tc = &DTC[tt];

  int i, j, k, l;
  double qx, qy, qz;
  double dx, dy, dz;
  double cx, cy, cz;
  double px, py, pz;
  double rx, ry, rz;
  double s, dist2, dist2max;

  //  printf("\n");

  *nlist = 2;

  // initialize list[0] to be the starting point p0. indA and indB are
  // both set to the index (0-4) of the closest point. s=0.
  list[0].s = 0;
  list[0].p.x = p0->x;
  list[0].p.y = p0->y;
  list[0].p.z = p0->z;

  for(i = 0, dist2max = 1.0e30; i < 4; i++)
    {
      dx = p0->x - DP[t->p[i]].x;
      dy = p0->y - DP[t->p[i]].y;
      dz = p0->z - DP[t->p[i]].z;

      dist2 = dx * dx + dy * dy + dz * dz;
      if(dist2 < dist2max)
        {

          dist2max = dist2;
          list[0].indA = list[0].indB = i;
        }
    }

  // initialize list[1] to be the end point p1. indA and indB are
  // both set to the index (0-4) of the closest point.
  list[1].s = 1.0;
  list[1].p.x = p1->x;
  list[1].p.y = p1->y;
  list[1].p.z = p1->z;

  for(i = 0, dist2max = 1.0e30; i < 4; i++)
    {
      dx = p1->x - DP[t->p[i]].x;
      dy = p1->y - DP[t->p[i]].y;
      dz = p1->z - DP[t->p[i]].z;

      dist2 = dx * dx + dy * dy + dz * dz;
      if(dist2 < dist2max)
        {
          dist2max = dist2;
          list[1].indA = list[1].indB = i;
        }
    }

  // c is the vector from the tetra center to tetra point 0
  // XXX CLOBBERED BELOW!
  cx = tc->cx - DP[t->p[0]].x;
  cy = tc->cy - DP[t->p[0]].y;
  cz = tc->cz - DP[t->p[0]].z;

  // d is the vector from entry to exit point (the full line segment)
  dx = p1->x - p0->x;
  dy = p1->y - p0->y;
  dz = p1->z - p0->z;

  // c is the vector from entry point p0 to tetra center 
  cx = tc->cx - p0->x;
  cy = tc->cy - p0->y;
  cz = tc->cz - p0->z;


  /*

     printf("\n%g %g %g %g\n",
     (t->cx - p0->x) * (t->cx - p0->x) + (t->cy - p0->y) * (t->cy - p0->y) + (t->cz - p0->z) * (t->cz - p0->z),
     (t->cx - p1->x) * (t->cx - p1->x) + (t->cy - p1->y) * (t->cy - p1->y) + (t->cz - p1->z) * (t->cz - p1->z),
     (t->cx - p2->x) * (t->cx - p2->x) + (t->cy - p2->y) * (t->cy - p2->y) + (t->cz - p2->z) * (t->cz - p2->z),
     (t->cx - p3->x) * (t->cx - p3->x) + (t->cy - p3->y) * (t->cy - p3->y) + (t->cz - p3->z) * (t->cz - p3->z));
   */


  for(i = 0; i < 6; i++)
    {
      // pairs is a [6][4] const array defined above that defines 6 different
      // permutations of the tetrahedron points 0,1,2,3.
      j = pairs[i][0];
      k = pairs[i][1];

      point *dpj = &DP[t->p[j]];
      point *dpk = &DP[t->p[k]];

      // q is the tetra edge vector from pk to pj.
      qx = dpj->x - dpk->x;
      qy = dpj->y - dpk->y;
      qz = dpj->z - dpk->z;

      // s = c.q / d.q
      // This is the standard formula for the intersection between a
      // ray and a plane. s is the point where the ray p0+s*d
      // intersects the plane which is perpendicular to q and goes
      // through c, i.e. the Voronoi face corresponding to the j-k
      // edge. 
      s = (cx * qx + cy * qy + cz * qz) / (dx * qx + dy * qy + dz * qz);

      if(s > 0 && s < 1)
        {
          /* The ray intersects the j-k Voronoi face inside the
             tetrahedron, which means we'll have at least one
             split. But we need to check whether it intersects another
             face before. */

          // set p to the intersection point
          px = p0->x + s * dx;
          py = p0->y + s * dy;
          pz = p0->z + s * dz;

          // r is the vector from pj to p
          rx = px - dpj->x;
          ry = py - dpj->y;
          rz = pz - dpj->z;

          /* Now we check whether p is really within voronoi cell l,
             by comparing the distance from p to l. Since Voronoi cell
             j is defined by the points closer to j than any other
             point. By construction p is on the j-k face, so is
             equally distant from j and k. We thus only need to check
             the distance from p to the two other points in the
             tetrahedron. */
          l = pairs[i][2];
          point *dpl = &DP[t->p[l]];

          // rr is the vector from pl to p
          double rrx, rry, rrz;
          rrx = px - dpl->x;
          rry = py - dpl->y;
          rrz = pz - dpl->z;


          //      if((rx*rx + ry*ry + rz*rz) <rad2)

          if((rx * rx + ry * ry + rz * rz) < (rrx * rrx + rry * rry + rrz * rrz))
            {
              /* r is shorter than rr, the intersection point p is
                 closer to pj than pl, which means p can not be in pl's
                 voronoi cell. */

              // Now, in the same way, we check whether p is really
              // within voronoi cell m.
              int m = pairs[i][3];
              point *dpm = &DP[t->p[m]];

              rrx = px - dpm->x;
              rry = py - dpm->y;
              rrz = pz - dpm->z;

              if((rx * rx + ry * ry + rz * rz) < (rrx * rrx + rry * rry + rrz * rrz))
                {
                  /* p is indeed on the j-k voronoi face. We have
                     found a new cut-point. Perform an insertion sort
                     on the list. */

                  for(l = 0; l < *nlist - 1; l++)
                    {
                      // find the place where to insert the cut point
                      if(list[l].s < s && s < list[l + 1].s)
                        {
                          // make space in the list and add the new element
                          memmove(&list[l + 2], &list[l + 1], (*nlist - (l + 1)) * sizeof(intersection_list));

                          list[l + 1].s = s;
                          list[l + 1].p.x = px;
                          list[l + 1].p.y = py;
                          list[l + 1].p.z = pz;
                          list[l + 1].indA = j;
                          list[l + 1].indB = k;
                          (*nlist)++;
                          break;
                        }
                    }
                }
            }

        }
    }
}


/** This little function returns the Voronoi cell that the line
    segment [i, i+1] in the intersection_list is internal to,
    expressed as the point in the tetra */
int decode_intersection_list(int i, intersection_list * list)
{
  if(list[i].indA == list[i + 1].indA)
    return list[i].indA;
  else if(list[i].indA == list[i + 1].indB)
    return list[i].indA;
  else if(list[i].indB == list[i + 1].indA)
    return list[i].indB;
  else if(list[i].indB == list[i + 1].indB)
    return list[i].indB;
  else
    return -1;
}

/** Calculates the surface density for tetrahedron tt for the line
    segment [pp0, pp1] (which is known to be inside this tetra). The
    strategy appears to be to split the line segment into segments
    that are within the 4 Voronoi cells that exist within the
    tetrahedron. For each of these segments, the density field is
    described by the density at the point (ie, Voronoi center) and the
    gradient. It is thus trivial to integrate the column density
    within that segment. Most of the logic goes into finding the
    splitting points.  */
void calc_picture_contribution(tessellation * T, int tt, point * p0, point * p1, double *sigma, double *sigmatemp, double *sigmaweight, int weight_flag, int gradients_flag
#if defined(COOLING) && !defined(GRACKLE)
                               , double *sigmaszy
#endif
#ifdef GFM_STELLAR_EVOLUTION
                               , double *sigmametal
#endif
#ifdef GFM_AGN_RADIATION
                               , double *sigmaagnbol
#endif
#ifdef TRACER_MC
                               , double *sigmatracernum, double *sigmatrweight
#endif
#ifdef CHEM_IMAGE
                               , double *sigmadust, double *sigmah2, double *sigmahp, double *sigmaco
#endif
  )
{
  tetra *DT = T->DT;
  point *DP = T->DP;

  int i, j, li;
  double dx, dy, dz;
  double mx, my, mz;

  tetra *t = &DT[tt];

  intersection_list list[8];
  int nlist;

  double rho, proj_rhotemp_sum = 0;

  calc_delaunay_intersections(T, tt, p0, p1, list, &nlist);

  double proj_sum = 0;
  double proj_weight_sum = 0;
#if defined(COOLING) && !defined(GRACKLE)
  double proj_szy_sum = 0;
#endif
#ifdef GFM_STELLAR_EVOLUTION
  double proj_rhometal_sum = 0;
#endif
#ifdef GFM_AGN_RADIATION
  double proj_rhoagnbol_sum = 0;
#endif
#ifdef TRACER_MC
  double proj_rhotracernum_sum = 0;
  double proj_trweight_sum = 0;
#endif
#ifdef CHEM_IMAGE
  double proj_dust_sum = 0;
  double proj_h2_sum = 0;
  double proj_hp_sum = 0;
  double proj_co_sum = 0;
#endif


  /* now that all cuts have been established, we can integrate the
     column density for the segments */

  for(i = 0; i < nlist - 1; i++)
    {
      //      printf("(%d|%d)  s0=%g s1=%g  (%d|%d)  (%d|%d)\n", i, nlist, list[i].s, list[i+1].s, list[i].indA, list[i].indB,  list[i+1].indA, list[i+1].indB);

      // figure out which voronoi cell the segment is internal to
      j = decode_intersection_list(i, list);

      // not sure what could cause this to fail...
      if(j >= 0)
        {
          // check that this voronoi cell indeed is on this processor
          if(DP[t->p[j]].task == ThisTask)
            {
              li = DP[t->p[j]].index;

              // check that this point actually has a hydro quantity
              // (what are failure scenarios?)
              if(li >= 0 && li < NumGas)
                {
                  /* now determine the interpolated density value at the
                     middle of the two points */

                  // Is the SphP[i].Center different from DP[i] ?? Yes.

                  mx = 0.5 * (list[i].p.x + list[i + 1].p.x) - SphP[li].Center[0];
                  my = 0.5 * (list[i].p.y + list[i + 1].p.y) - SphP[li].Center[1];
                  mz = 0.5 * (list[i].p.z + list[i + 1].p.z) - SphP[li].Center[2];

                  if(gradients_flag == 1)
                    rho = SphP[li].Density + SphP[li].Grad.drho[0] * mx + SphP[li].Grad.drho[1] * my + SphP[li].Grad.drho[2] * mz;
                  else if(gradients_flag == 0)
                    rho = SphP[li].Density;
                  else
                    terminate("gradients_flag != 0 && gradients_flag != 1");

                  // get the column length (this seems suboptimal, we already have s)
                  // r=(list[i+1].s-list[i].s)*l_seg;
                  dx = list[i].p.x - list[i + 1].p.x;
                  dy = list[i].p.y - list[i + 1].p.y;
                  dz = list[i].p.z - list[i + 1].p.z;

                  double r = sqrt(dx * dx + dy * dy + dz * dz);

                  if(weight_flag > 0)
                    {
                      double pow_rho = pow(rho, weight_flag), pow_rho_sum = pow(rho, weight_flag - 1);
                      proj_sum += pow_rho * r;
#ifdef VORONOI_PROJ_TEMP
                      double meanweight = 4.0 / (1. + 3. * HYDROGEN_MASSFRAC) * PROTONMASS;
#if defined(COOLING) && !defined(GRACKLE)
                      meanweight = 4.0 / (1. + 3. * HYDROGEN_MASSFRAC + 4. * HYDROGEN_MASSFRAC * SphP[li].Ne) * PROTONMASS;
#endif
                      double temp_in_K = GAMMA_MINUS1 * SphP[li].Utherm / BOLTZMANN * All.UnitEnergy_in_cgs / All.UnitMass_in_g * meanweight;
#ifdef SGCHEM
                      double yn, en, yntot;

                      yn = rho * All.UnitDensity_in_cgs / ((1.0 + 4.0 * ABHE) * PROTONMASS);
                      en = SphP[li].Utherm * rho * All.UnitEnergy_in_cgs / pow(All.UnitLength_in_cm, 3);
#if CHEMISTRYNETWORK == 1
                      yntot = (1.0 + ABHE - SphP[li].TracAbund[IH2] + SphP[li].TracAbund[IHP] + SphP[li].TracAbund[IHEP] + SphP[li].TracAbund[IHEPP]) * yn;
#else
                      yntot = (1.0 + ABHE - SphP[li].TracAbund[IH2] + SphP[li].TracAbund[IHP]) * yn;
#endif
                      temp_in_K = (GAMMA - 1.0) * en / (yntot * BOLTZMANN);
#endif /* SGCHEM */

                      proj_rhotemp_sum += temp_in_K * pow_rho * r;
#if defined(COOLING) && !defined(GRACKLE)
                      proj_szy_sum += temp_in_K * SphP[li].Ne * (rho * All.UnitDensity_in_cgs * HYDROGEN_MASSFRAC / PROTONMASS) * (r * All.UnitLength_in_cm);   // in pysical units
#endif
#else
                      proj_rhotemp_sum += SphP[li].Utherm * pow_rho * r;
#endif

#ifdef GFM_STELLAR_EVOLUTION
                      proj_rhometal_sum += SphP[li].Metallicity * pow_rho * r;
#endif
#ifdef GFM_AGN_RADIATION
                      proj_rhoagnbol_sum += SphP[li].AGNBolIntensity * pow_rho * r;
#endif
#ifdef TRACER_MC
                      proj_rhotracernum_sum += pow(get_number_of_tracers(li) * All.ReferenceTracerMCMass / P[li].Mass, 3.0) * pow_rho * r;
                      proj_trweight_sum += pow(get_number_of_tracers(li) * All.ReferenceTracerMCMass / P[li].Mass, 2.0) * pow_rho_sum * r;
#endif
#ifdef CHEM_IMAGE
                      proj_dust_sum += SphP[li].DustTemp * pow_rho * r;
                      proj_h2_sum += SphP[li].TracAbund[IH2] * pow_rho * r;
                      proj_hp_sum += SphP[li].TracAbund[IHP] * pow_rho * r;
                      proj_co_sum += SphP[li].TracAbund[ICO] * pow_rho * r;
#endif

                      proj_weight_sum += pow_rho_sum * r;
                    }
                  else
                    {
                      proj_sum += rho * r;
#ifdef VORONOI_PROJ_TEMP
                      double meanweight = 4.0 / (1. + 3. * HYDROGEN_MASSFRAC) * PROTONMASS;
#if defined(COOLING) && !defined(GRACKLE)
                      meanweight = 4.0 / (1. + 3. * HYDROGEN_MASSFRAC + 4. * HYDROGEN_MASSFRAC * SphP[li].Ne) * PROTONMASS;
#endif
                      double temp_in_K = GAMMA_MINUS1 * SphP[li].Utherm / BOLTZMANN * All.UnitEnergy_in_cgs / All.UnitMass_in_g * meanweight;
#ifdef SGCHEM
                      double yn, en, yntot;

                      yn = rho * All.UnitDensity_in_cgs / ((1.0 + 4.0 * ABHE) * PROTONMASS);
                      en = SphP[li].Utherm * rho * All.UnitEnergy_in_cgs / pow(All.UnitLength_in_cm, 3);
#if CHEMISTRYNETWORK == 1
                      yntot = (1.0 + ABHE - SphP[li].TracAbund[IH2] + SphP[li].TracAbund[IHP] + SphP[li].TracAbund[IHEP] + SphP[li].TracAbund[IHEPP]) * yn;
#else
                      yntot = (1.0 + ABHE - SphP[li].TracAbund[IH2] + SphP[li].TracAbund[IHP]) * yn;
#endif
                      temp_in_K = (GAMMA - 1.0) * en / (yntot * BOLTZMANN);
#endif /* SGCHEM */

                      proj_rhotemp_sum += temp_in_K * rho * r;
#if defined(COOLING) && !defined(GRACKLE)
                      proj_szy_sum += temp_in_K * SphP[li].Ne * (rho * All.UnitDensity_in_cgs * HYDROGEN_MASSFRAC / PROTONMASS) * (r * All.UnitLength_in_cm);
#endif
#else
                      proj_rhotemp_sum += SphP[li].Utherm * rho * r;
#endif


#ifdef GFM_STELLAR_EVOLUTION
                      proj_rhometal_sum += SphP[li].Metallicity * rho * r;
#endif
#ifdef GFM_AGN_RADIATION
                      proj_rhoagnbol_sum += SphP[li].AGNBolIntensity * rho * r;
#endif
#ifdef TRACER_MC
                      proj_rhotracernum_sum += rho * r * get_number_of_tracers(li) * All.ReferenceTracerMCMass / P[li].Mass;
#endif
#ifdef CHEM_IMAGE
                      proj_dust_sum += SphP[li].DustTemp * rho * r;
                      proj_h2_sum += SphP[li].TracAbund[IH2] * rho * r;
                      proj_hp_sum += SphP[li].TracAbund[IHP] * rho * r;
                      proj_co_sum += SphP[li].TracAbund[ICO] * rho * r;
#endif

                    }
                }
            }
        }
    }

  *sigma = proj_sum;
  *sigmatemp = proj_rhotemp_sum;
  *sigmaweight = proj_weight_sum;
#if defined(COOLING) && !defined(GRACKLE)
  *sigmaszy = proj_szy_sum;
#endif
#ifdef GFM_STELLAR_EVOLUTION
  *sigmametal = proj_rhometal_sum;
#endif
#ifdef GFM_AGN_RADIATION
  *sigmaagnbol = proj_rhoagnbol_sum;
#endif
#ifdef TRACER_MC
  *sigmatracernum = proj_rhotracernum_sum;
  *sigmatrweight = proj_trweight_sum;
#endif
#ifdef CHEM_IMAGE
  *sigmadust = proj_dust_sum;
  *sigmah2 = proj_h2_sum;
  *sigmahp = proj_hp_sum;
  *sigmaco = proj_co_sum;
#endif
}

#endif



void extract_position_from_axes(int i, int j, int xaxis, int yaxis, int zaxis, int pixels_x, int pixels_y, double xmin, double xmax, double ymin, double ymax, double zval, double *p)
{
  // The logic below sets the point p depending on which axes
  // we have chosen.
  if(xaxis == 0 && yaxis == 1)
    {
      p[0] = (i + 0.5) / pixels_x * (xmax - xmin) + xmin;
      p[1] = (j + 0.5) / pixels_y * (ymax - ymin) + ymin;
      p[2] = zval;
    }
  else if(xaxis == 1 && yaxis == 0)
    {
      p[0] = (j + 0.5) / pixels_y * (ymax - ymin) + ymin;
      p[1] = (i + 0.5) / pixels_x * (xmax - xmin) + xmin;
      p[2] = zval;
    }
  else if(xaxis == 0 && yaxis == 2)
    {
      p[0] = (i + 0.5) / pixels_x * (xmax - xmin) + xmin;
      p[1] = zval;
      p[2] = (j + 0.5) / pixels_y * (ymax - ymin) + ymin;
    }
  else if(xaxis == 2 && yaxis == 0)
    {
      p[0] = (j + 0.5) / pixels_y * (ymax - ymin) + ymin;
      p[1] = zval;
      p[2] = (i + 0.5) / pixels_x * (xmax - xmin) + xmin;
    }
  else if(xaxis == 1 && yaxis == 2)
    {
      p[0] = zval;
      p[1] = (i + 0.5) / pixels_x * (xmax - xmin) + xmin;
      p[2] = (j + 0.5) / pixels_y * (ymax - ymin) + ymin;
    }
  else if(xaxis == 2 && yaxis == 1)
    {
      p[0] = zval;
      p[1] = (j + 0.5) / pixels_y * (ymax - ymin) + ymin;
      p[2] = (i + 0.5) / pixels_x * (xmax - xmin) + xmin;
    }
  else
    terminate("invalid combination of axes");
}


/** Generates an image of the density field in a slice given by the
    specified coordinate ranges.
*/
void make_3d_voronoi_slice_image(int num, int gradients_flag, int pixels_x, int pixels_y, int xaxis, int yaxis, int zaxis, double xmin, double xmax, double ymin, double ymax, double zval)
{
  CPU_Step[CPU_MISC] += measure_time();

  float *density = 0, *temp = 0, *met = 0, *velocity = 0, *Bfield = 0, *vorticity = 0, *photon_density = 0, *chem_elements = 0, *density_trmc = 0;
  FILE *fd = 0, *fdtemp = 0, *fdmet = 0, *fdvel = 0, *fdmag = 0, *fdvort = 0, *fdphot = 0, *fdchem = 0, *fdtr = 0;
  char buf[1000];
#ifdef CHEM_IMAGE
  float *dust, *xH2, *xHP, *xCO;
  FILE *fddust = 0, *fdh2 = 0, *fdhp = 0, *fdco = 0;
#endif

  if(gradients_flag == 1)
    sprintf(buf, "slice");
  else if(gradients_flag == 0)
    sprintf(buf, "slice_nograds");
  else
    terminate("gradients_flag != 1 && gradients_flag != 0");

  mpi_printf("we start to generate the slice image... gradients_flag=%d\n", gradients_flag);

  open_image_files(buf, num, &fd, &fdtemp, &fdmet, &fdvel,
#ifdef MHD
                   &fdmag,
#else
                   0,
#endif
#ifdef OUTPUT_VORTICITY
                   &fdvort,
#else
                   0,
#endif
#ifdef RT_ADVECT
                   &fdphot,
#else
                   0,
#endif
#ifdef GFM_STELLAR_EVOLUTION
                   &fdchem,
#else
                   0,
#endif
#ifdef TRACER_MC
                   &fdtr,
#else
                   0,
#endif
#ifdef CHEM_IMAGE
                   &fddust, &fdh2, &fdhp, &fdco
#else
                   0, 0, 0, 0
#endif
    );

  write_image_header(fd, pixels_x, pixels_y, 0);
  write_image_header(fdtemp, pixels_x, pixels_y, 0);
  write_image_header(fdmet, pixels_x, pixels_y, 0);
  write_image_header(fdvel, pixels_x, pixels_y, 0);
  write_image_header(fdmag, pixels_x, pixels_y, 0);
  write_image_header(fdvort, pixels_x, pixels_y, 0);
  write_image_header(fdphot, pixels_x, pixels_y, 0);
  write_image_header(fdchem, pixels_x, pixels_y, 0);
  write_image_header(fdtr, pixels_x, pixels_y, 0);
#ifdef CHEM_IMAGE
  write_image_header(fddust, pixels_x, pixels_y, 0);
  write_image_header(fdh2, pixels_x, pixels_y, 0);
  write_image_header(fdhp, pixels_x, pixels_y, 0);
  write_image_header(fdco, pixels_x, pixels_y, 0);
#endif

  density = mymalloc("density", pixels_x * pixels_y * sizeof(float));
  temp = mymalloc("temp", pixels_x * pixels_y * sizeof(float));
  met = mymalloc("met", pixels_x * pixels_y * sizeof(float));
  velocity = mymalloc("velocity", 3 * pixels_x * pixels_y * sizeof(float));
  Bfield = mymalloc("Bfield", 3 * pixels_x * pixels_y * sizeof(float));
  vorticity = mymalloc("vorticity", 3 * pixels_x * pixels_y * sizeof(float));
  photon_density = mymalloc("photon_density", RT_N_DIR * pixels_x * pixels_y * sizeof(float));
  chem_elements = mymalloc("chem_elements", GFM_N_CHEM_ELEMENTS * pixels_x * pixels_y * sizeof(float));
  density_trmc = mymalloc("density_trmc", pixels_x * pixels_y * sizeof(float));
#ifdef CHEM_IMAGE
  dust = mymalloc("dust", pixels_x * pixels_y * sizeof(float));
  xH2 = mymalloc("xH2", pixels_x * pixels_y * sizeof(float));
  xHP = mymalloc("xHP", pixels_x * pixels_y * sizeof(float));
  xCO = mymalloc("xCO", pixels_x * pixels_y * sizeof(float));
#endif

  MaxNray = pixels_x * pixels_y;
  Ray = mymalloc_movable(&Ray, "Ray", MaxNray * sizeof(ray_data));

  /* For simplicity, we use the "Ray" struct from the projection code. */

  setup_rays(pixels_x, pixels_y, xaxis, yaxis, zaxis, xmin, xmax, ymin, ymax, zval, zval);

  fill_slice(Ray, Nray, gradients_flag, pixels_x, pixels_y, density, temp, met, velocity, Bfield, vorticity, photon_density, chem_elements, density_trmc
#ifdef CHEM_IMAGE
             , dust, xH2, xHP, xCO
#endif
    );

  my_fwrite(density, sizeof(float), pixels_x * pixels_y, fd);
  my_fwrite(temp, sizeof(float), pixels_x * pixels_y, fdtemp);
  my_fwrite(met, sizeof(float), pixels_x * pixels_y, fdmet);
  my_fwrite(velocity, sizeof(float), 3 * pixels_x * pixels_y, fdvel);
  my_fwrite(Bfield, sizeof(float), 3 * pixels_x * pixels_y, fdmag);
  my_fwrite(vorticity, sizeof(float), 3 * pixels_x * pixels_y, fdvort);
  my_fwrite(photon_density, sizeof(float), RT_N_DIR * pixels_x * pixels_y, fdphot);
  my_fwrite(chem_elements, sizeof(float), GFM_N_CHEM_ELEMENTS * pixels_x * pixels_y, fdchem);
  my_fwrite(density_trmc, sizeof(float), pixels_x * pixels_y, fdtr);
#ifdef CHEM_IMAGE
  my_fwrite(dust, sizeof(float), pixels_x * pixels_y, fddust);
  my_fwrite(xH2, sizeof(float), pixels_x * pixels_y, fdh2);
  my_fwrite(xHP, sizeof(float), pixels_x * pixels_y, fdhp);
  my_fwrite(xCO, sizeof(float), pixels_x * pixels_y, fdco);

  if(fddust)
    fclose(fddust);
  if(fdh2)
    fclose(fdh2);
  if(fdhp)
    fclose(fdhp);
  if(fdco)
    fclose(fdco);
#endif

  if(fd)
    fclose(fd);
  if(fdtemp)
    fclose(fdtemp);
  if(fdmet)
    fclose(fdmet);
  if(fdvel)
    fclose(fdvel);
  if(fdmag)
    fclose(fdmag);
  if(fdvort)
    fclose(fdvort);
  if(fdphot)
    fclose(fdphot);
  if(fdchem)
    fclose(fdchem);
  if(fdtr)
    fclose(fdtr);

  myfree(Ray);
#ifdef CHEM_IMAGE
  myfree(xCO);
  myfree(xHP);
  myfree(xH2);
  myfree(dust);
#endif
  myfree(density_trmc);
  myfree(chem_elements);
  myfree(photon_density);
  myfree(vorticity);
  myfree(Bfield);
  myfree(velocity);
  myfree(met);
  myfree(temp);
  myfree(density);

  CPU_Step[CPU_MAKEIMAGES] += measure_time();
}


FILE *FdFaces;

/*! Writes a file listing all the Voronoi faces that intersect the
  slice plane as a list of line segments (x0,y0,x1,y1). Note that all
  intersections across the entire box are written, regardless of the
  x/y limits made on the image slice. Also note that this only works
  correctly on one processor because no attempt is made to gather the
  tessellations from a decomposed domain.  */
void make_3d_voronoi_listfaces(tessellation * T, int num, int xaxis, int yaxis, int zaxis, double zval)
{
  CPU_Step[CPU_MISC] += measure_time();

  int i, nr, bit;
  char buf[1000], msg[1000];

  tetra_center *tmpDTC = T->DTC;
  T->DTC = mymalloc_movable(&T->DTC, "DTC", T->MaxNdt * sizeof(tetra_center));
  T->DTF = mymalloc_movable(&T->DTF, "DTF", T->MaxNdt * sizeof(char));
  for(i = 0; i < T->Ndt; i++)
    T->DTF[i] = 0;
  compute_circumcircles(T);

  Edge_visited = mymalloc("Edge_visited", T->Ndt * sizeof(unsigned char));

  for(i = 0; i < T->Ndt; i++)
    Edge_visited[i] = 0;

  sprintf(buf, "%s/faces_list_%03d.txt", All.OutputDir, num);

  if(ThisTask == 0)
    {
      if(!(FdFaces = fopen(buf, "w")))
        {
          sprintf(msg, "can't open file `%s' for writing snapshot.\n", buf);
          terminate(msg);
        }
    }


  for(i = 0; i < T->Ndt; i++)
    {
      if(T->DT[i].t[0] < 0)     /* deleted? */
        continue;

      bit = 1;
      nr = 0;

      // Loop over all edges in tetra i.
      while(Edge_visited[i] != EDGE_ALL)
        {
          if((Edge_visited[i] & bit) == 0)
            make_3d_voronoi_listfaces_check_for_cut(T, i, nr, xaxis, yaxis, zaxis, zval);

          bit <<= 1;
          nr++;
        }
    }

  if(ThisTask == 0)
    fclose(FdFaces);

  myfree(Edge_visited);

  myfree(T->DTF);
  myfree(T->DTC);
  T->DTC = tmpDTC;

  CPU_Step[CPU_MAKEIMAGES] += measure_time();
}


/*! Checks whether the triangle defined by the three points in pp
  intersects the slice plane, and if so writes it to the output
  file.  */
void check_for_cut(double pp[3][3], int xaxis, int yaxis, int zaxis, double zval)
{
  int i, ii, count;
  double x[3], y[3], w;

  for(i = 0, count = 0; i < 3; i++)
    {
      ii = i + 1;
      if(ii > 2)
        ii = 0;

      // The edge only intersects the slice plane if the points are on
      // opposite sides.
      if((pp[i][zaxis] - zval) * (pp[ii][zaxis] - zval) < 0)
        {
          // determine the x and y of the intersection point
          w = (zval - pp[i][zaxis]) / (pp[ii][zaxis] - pp[i][zaxis]);
          x[count] = pp[i][xaxis] + w * (pp[ii][xaxis] - pp[i][xaxis]);
          y[count] = pp[i][yaxis] + w * (pp[ii][yaxis] - pp[i][yaxis]);
          count++;
        }
    }

  if(count == 2 && ThisTask == 0)
    fprintf(FdFaces, "%g %g %g %g\n", x[0], y[0], x[1], y[1]);
}



/*! Checks whether the face corresponding to edge nr in tetrahedron tt
    intersects the slice plane. */
void make_3d_voronoi_listfaces_check_for_cut(tessellation * T, int tt, int nr, int xaxis, int yaxis, int zaxis, double zval)
{
  point *DP = T->DP;
  tetra *DT = T->DT;
  tetra_center *DTC = T->DTC;
  int i, j, k, l, m, ii, jj, kk, ll, count, nr_next, flag, nn;
  tetra *prev, *next;
  tetra_center *prevc, *nextc;
  double pp[3][3];

  tetra *t = &DT[tt];
  tetra_center *tc = &DTC[tt];

  i = edge_start[nr];
  j = edge_end[nr];
  k = edge_opposite[nr];
  l = edge_nexttetra[nr];

  Edge_visited[tt] |= (1 << nr);

  if(T->Nvf + 1 >= T->MaxNvf)
    terminate("Nvf + 1 >= MaxNvf");


  pp[0][0] = tc->cx;
  pp[0][1] = tc->cy;
  pp[0][2] = tc->cz;

  count = 0;

  flag = 0;

  if(DP[t->p[i]].task == ThisTask && DP[t->p[i]].index >= 0 && DP[t->p[i]].index < NumGas)
    flag = 1;

  if(DP[t->p[j]].task == ThisTask && DP[t->p[j]].index >= 0 && DP[t->p[j]].index < NumGas)
    flag = 1;

  prev = t;
  prevc = tc;
  do
    {
      nn = prev->t[l];
      next = &DT[nn];
      nextc = &DTC[nn];

      if(prev != t && next != t)
        {
          pp[1][0] = prevc->cx;
          pp[1][1] = prevc->cy;
          pp[1][2] = prevc->cz;

          pp[2][0] = nextc->cx;
          pp[2][1] = nextc->cy;
          pp[2][2] = nextc->cz;

          if(flag)
            check_for_cut(pp, xaxis, yaxis, zaxis, zval);
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
        terminate("ll < 0 || ii < 0 || jj < 0");

      kk = 6 - (ll + ii + jj);


      /* need to determine the edge number to be able to flag it */

      for(nr_next = 0; nr_next < 6; nr_next++)
        if((edge_start[nr_next] == ii && edge_end[nr_next] == jj) || (edge_start[nr_next] == jj && edge_end[nr_next] == ii))
          {
            if((Edge_visited[nn] & (1 << nr_next)) && next != t)
              terminate("can't be");

            Edge_visited[nn] |= (1 << nr_next);
            break;
          }

      prevc = nextc;
      prev = next;
      i = ii;
      l = ll;
      j = jj;
      k = kk;

      count++;

      if(count > 1000)
        terminate("count > 1000");
    }
  while(next != t);

}

#endif


#ifdef TWODIMS


void make_2d_voronoi_image(int num, int pixels_x, int pixels_y)
{
  CPU_Step[CPU_MISC] += measure_time();

  char buf[1000], msg[1000];
  float *dens, *denssum, *dp;
  FILE *fd = 0;
  point *p;
  int pp;
  int tt0, ttstart, ttrow;
  tetra *t0;
  double l_dx, l_dy;
  int i, j, k, kmin, li, moves, ret, no, task;
  double r2, r2min, rho_L;
  peanokey key;

  if(Mesh.Ndp >= Mesh.MaxNdp)
    {
      terminate("Ndp >= MaxNdp");
    }

  pp = Mesh.Ndp;
  p = &Mesh.DP[pp];


  sprintf(buf, "%s/density_field_%03d", All.OutputDir, num);

  if(ThisTask == 0)
    {
      if(!(fd = fopen(buf, "w")))
        {
          sprintf(msg, "can't open file `%s' for writing snapshot.\n", buf);
          terminate(msg);
        }

      my_fwrite(&pixels_x, sizeof(int), 1, fd);
      my_fwrite(&pixels_y, sizeof(int), 1, fd);
    }

  dens = mymalloc("dens", pixels_x * pixels_y * sizeof(float));
  denssum = mymalloc("denssum", pixels_x * pixels_y * sizeof(float));

  for(i = 0, dp = dens; i < pixels_x; i++)
    for(j = 0; j < pixels_y; j++)
      *dp++ = 0;

  ttrow = 0;

  point *DP = Mesh.DP;

  for(i = 0; i < pixels_x; i++)
    {
      ttstart = ttrow;

      for(j = 0; j < pixels_y; j++)
        {
          p->x = (i + 0.5) / pixels_x * boxSize_X;
          p->y = (j + 0.5) / pixels_y * boxSize_Y;
          p->z = 0;

          key = peano_hilbert_key((int) ((p->x - DomainCorner[0]) * DomainFac), (int) ((p->y - DomainCorner[1]) * DomainFac), (int) ((p->z - DomainCorner[2]) * DomainFac), BITS_PER_DIMENSION);

          no = 0;
          while(TopNodes[no].Daughter >= 0)
            no = TopNodes[no].Daughter + (key - TopNodes[no].StartKey) / (TopNodes[no].Size >> 3);

          no = TopNodes[no].Leaf;
          task = DomainTask[no];

          if(task == ThisTask)
            {
#ifndef OPTIMIZE_MEMORY_USAGE
              set_integers_for_point(&Mesh, pp);
#endif
              tt0 = get_triangle(&Mesh, pp, &moves, &ret, ttstart);
              t0 = &Mesh.DT[tt0];

              for(k = 0, kmin = 0, r2min = 1.0e30; k < 3; k++)
                {
                  r2 = (p->x - DP[t0->p[k]].x) * (p->x - DP[t0->p[k]].x) + (p->y - DP[t0->p[k]].y) * (p->y - DP[t0->p[k]].y);
                  if(r2 < r2min)
                    {
                      r2min = r2;
                      kmin = k;
                    }
                }

              li = DP[t0->p[kmin]].index;

              if(li >= NumGas)
                li -= NumGas;

              if(DP[t0->p[kmin]].task == ThisTask)
                {
                  l_dx = p->x - SphP[li].Center[0];
                  l_dy = p->y - SphP[li].Center[1];
                }
              else
                {
                  l_dx = p->x - PrimExch[li].Center[0];
                  l_dy = p->y - PrimExch[li].Center[1];
                }
#ifdef PERIODIC
#if !defined(REFLECTIVE_X)
              if(l_dx < -boxHalf_X)
                l_dx += boxSize_X;
              if(l_dx > boxHalf_X)
                l_dx -= boxSize_X;
#endif
#if !defined(REFLECTIVE_Y)
              if(l_dy < -boxHalf_Y)
                l_dy += boxSize_Y;
              if(l_dy > boxHalf_Y)
                l_dy -= boxSize_Y;
#endif
#endif
              if(DP[t0->p[kmin]].task == ThisTask)
                {
                  rho_L = SphP[li].Density + SphP[li].Grad.drho[0] * l_dx + SphP[li].Grad.drho[1] * l_dy;
                }
              else
                {
                  rho_L = PrimExch[li].Density + GradExch[li].drho[0] * l_dx + GradExch[li].drho[1] * l_dy;
                }

              dens[i * pixels_y + j] = rho_L;

              ttstart = tt0;

              if(j == 0)
                ttrow = tt0;
            }
        }
    }


  MPI_Reduce(dens, denssum, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      my_fwrite(denssum, sizeof(float), pixels_x * pixels_y, fd);
      fclose(fd);
    }

  myfree(denssum);
  myfree(dens);

  CPU_Step[CPU_MAKEIMAGES] += measure_time();
}


void make_2d_voronoi_image_zoomed(tessellation * T, int num, int pixels_x, int pixels_y, double xmin, double xmax, double ymin, double ymax)
{
  CPU_Step[CPU_MISC] += measure_time();

  float *density = 0, *temp = 0, *velocity = 0, *vorticity = 0;
  FILE *fd = 0, *fdtemp = 0, *fdvel = 0, *fdvort = 0;

  open_image_files("field", num, &fd, &fdtemp, 0, &fdvel, 0,
#ifdef OUTPUT_VORTICITY
                   &fdvort,
#else
                   0,
#endif
                   0, 0, 0, 0, 0, 0, 0);

  write_image_header(fd, pixels_x, pixels_y, 0);
  write_image_header(fdtemp, pixels_x, pixels_y, 0);
  write_image_header(fdvel, pixels_x, pixels_y, 0);
  write_image_header(fdvort, pixels_x, pixels_y, 0);

  density = mymalloc("density", pixels_x * pixels_y * sizeof(float));
  temp = mymalloc("temp", pixels_x * pixels_y * sizeof(float));
  velocity = mymalloc("velocity", 3 * pixels_x * pixels_y * sizeof(float));
  vorticity = mymalloc("vorticity", 3 * pixels_x * pixels_y * sizeof(float));

  MaxNray = pixels_x * pixels_y;
  Ray = mymalloc_movable(&Ray, "Ray", MaxNray * sizeof(ray_data));

  /* For simplicity, we use the "Ray" struct from the projection
     code. */

  setup_rays(pixels_x, pixels_y, 0, 1, 2, xmin, xmax, ymin, ymax, 0.0, 0.0);

  fill_slice(Ray, Nray, 1, pixels_x, pixels_y, density, temp, 0, velocity, 0, vorticity, 0, 0, 0);

  my_fwrite(density, sizeof(float), pixels_x * pixels_y, fd);
  my_fwrite(temp, sizeof(float), pixels_x * pixels_y, fdtemp);
  my_fwrite(velocity, sizeof(float), 3 * pixels_x * pixels_y, fdvel);
  my_fwrite(vorticity, sizeof(float), 3 * pixels_x * pixels_y, fdvort);

  if(fd)
    fclose(fd);
  if(fdtemp)
    fclose(fdtemp);
  if(fdvel)
    fclose(fdvel);
  if(fdvort)
    fclose(fdvort);



  myfree(Ray);
  myfree(vorticity);
  myfree(velocity);
  myfree(temp);
  myfree(density);

  CPU_Step[CPU_MAKEIMAGES] += measure_time();
}


#endif


#if !defined(ONEDIMS) && !defined(TWODIMS)

void make_3d_voronoi_grid(int num, int pixels_x, int pixels_y, int pixels_z, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
{
  CPU_Step[CPU_MISC] += measure_time();

  float *density = 0, *temp = 0, *met = 0, *velocity = 0, *Bfield = 0, *vorticity = 0, *photon_density = 0, *chem_elements = 0, *density_trmc = 0;
  FILE *fd = 0, *fdvel = 0, *fdtemp = 0, *fdmet = 0, *fdmag = 0, *fdvort = 0, *fdphot = 0, *fdchem = 0, *fdtr = 0;
#ifdef CHEM_IMAGE
  float *dust = 0, *xH2 = 0, *xHP = 0, *xCO = 0;
  FILE *fddust = 0, *fdh2 = 0, *fdhp = 0, *fdco = 0;
#endif
  int i, j, gradients_flag;
  double xstep;

#ifdef VORONOI_NOGRADS
  gradients_flag = 0;
#else
  gradients_flag = 1;
#endif

  open_image_files("grid", num, &fd, &fdtemp, &fdmet, &fdvel,
#ifdef MHD
                   &fdmag,
#else
                   0,
#endif
#ifdef OUTPUT_VORTICITY
                   &fdvort,
#else
                   0,
#endif
#ifdef RT_ADVECT
                   &fdphot,
#else
                   0,
#endif
#ifdef GFM_STELLAR_EVOLUTION
                   &fdchem,
#else
                   0,
#endif
#ifdef TRACER_MC
                   &fdtr,
#else
                   0,
#endif
#ifdef CHEM_IMAGE
                   &fddust, &fdh2, &fdhp, &fdco
#else
                   0, 0, 0, 0
#endif
    );

  write_image_header(fd, pixels_x, pixels_y, pixels_z);
  write_image_header(fdtemp, pixels_x, pixels_y, pixels_z);
  write_image_header(fdmet, pixels_x, pixels_y, pixels_z);
  write_image_header(fdvel, pixels_x, pixels_y, pixels_z);
  write_image_header(fdmag, pixels_x, pixels_y, pixels_z);
  write_image_header(fdvort, pixels_x, pixels_y, pixels_z);
  write_image_header(fdphot, pixels_x, pixels_y, pixels_z);
  write_image_header(fdchem, pixels_x, pixels_y, pixels_z);
  write_image_header(fdtr, pixels_x, pixels_y, pixels_z);

#ifdef CHEM_IMAGE
  write_image_header(fddust, pixels_x, pixels_y, pixels_z);
  write_image_header(fdh2, pixels_x, pixels_y, pixels_z);
  write_image_header(fdhp, pixels_x, pixels_y, pixels_z);
  write_image_header(fdco, pixels_x, pixels_y, pixels_z);
#endif

  MaxNray = pixels_y * pixels_z;
  Ray = mymalloc_movable(&Ray, "Ray", MaxNray * sizeof(ray_data));

  xstep = (xmax - xmin) / pixels_x;

  // First create rays at the starting positions in the xmin plane,
  // with length appropriate to take them to the next x-slice
  setup_rays(pixels_y, pixels_z, 1, 2, 0, ymin, ymax, zmin, zmax, xmin, xmin + xstep);

  density = mymalloc("density", pixels_y * pixels_z * sizeof(float));
  temp = mymalloc("temp", pixels_y * pixels_z * sizeof(float));
  met = mymalloc("met", pixels_y * pixels_z * sizeof(float));
  velocity = mymalloc("velocity", 3 * pixels_y * pixels_z * sizeof(float));
  Bfield = mymalloc("Bfield", 3 * pixels_y * pixels_z * sizeof(float));
  vorticity = mymalloc("vorticity", 3 * pixels_y * pixels_z * sizeof(float));
  photon_density = mymalloc("photon_density", RT_N_DIR * pixels_y * pixels_z * sizeof(float));
  chem_elements = mymalloc("chem_elements", GFM_N_CHEM_ELEMENTS * pixels_y * pixels_z * sizeof(float));
  density_trmc = mymalloc("density_trmc", pixels_y * pixels_z * sizeof(float));

#ifdef CHEM_IMAGE
  dust = mymalloc("dust", pixels_y * pixels_z * sizeof(float));
  xH2 = mymalloc("H2", pixels_y * pixels_z * sizeof(float));
  xHP = mymalloc("HP", pixels_y * pixels_z * sizeof(float));
  xCO = mymalloc("CO", pixels_y * pixels_z * sizeof(float));
#endif

  mpi_printf("Extracting grid info...\n");

  // loop over x-slices
  for(i = 0; i < pixels_x; i++)
    {

      // reset ray propagation lengths
      for(j = 0; j < Nray; ++j)
        {
          Ray[j].len = 0;
          Ray[j].target_len = xstep;
        }

      // now we simply use fill_slice for the rays
      fill_slice(Ray, Nray, gradients_flag, pixels_y, pixels_z, density, temp, met, velocity, Bfield, vorticity, photon_density, chem_elements, density_trmc
#ifdef CHEM_IMAGE
                 , dust, xH2, xHP, xCO
#endif
        );

#ifdef VORONOI_DYNAMIC_UPDATE
      // step rays to the next x-position
      int left_this_task, rays_left;
      do
        {
          left_this_task = advance_rays_for_one_cell(0, 0, 1);
          exchange_rays();

          MPI_Allreduce(&left_this_task, &rays_left, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        }
      while(rays_left);
#else
      /* without being able to trace the rays, we must redo the mesh
         search from scratch */
      setup_rays(pixels_y, pixels_z, 1, 2, 0, ymin, ymax, zmin, zmax, xmin + (i + 1) * xstep, xmax);
#endif

      // write the slice
      if(ThisTask == 0)
        {
          mpi_printf("Writing slice %d/%d      \r", i, pixels_x);
          int b, j, k;
          for(j = 0; j < pixels_y; j++)
            for(k = 0; k < pixels_z; k++)
              for(b = 0; b < 3; b++)
                {
                  myassert(gsl_finite(velocity[3 * (j * pixels_z + k) + b]));
                }

          // these don't write if the fd is null
          my_fwrite(density, sizeof(float), pixels_y * pixels_z, fd);
          my_fwrite(temp, sizeof(float), pixels_y * pixels_z, fdtemp);
          my_fwrite(met, sizeof(float), pixels_y * pixels_z, fdmet);
          my_fwrite(velocity, sizeof(float), 3 * pixels_y * pixels_z, fdvel);
          my_fwrite(Bfield, sizeof(float), 3 * pixels_y * pixels_z, fdmag);
          my_fwrite(vorticity, sizeof(float), 3 * pixels_y * pixels_z, fdvort);
          my_fwrite(photon_density, sizeof(float), RT_N_DIR * pixels_y * pixels_z, fdvort);
          my_fwrite(chem_elements, sizeof(float), GFM_N_CHEM_ELEMENTS * pixels_y * pixels_z, fdchem);
          my_fwrite(density_trmc, sizeof(float), pixels_y * pixels_z, fdtr);
#ifdef CHEM_IMAGE
          my_fwrite(dust, sizeof(float), pixels_y * pixels_z, fddust);
          my_fwrite(xH2, sizeof(float), pixels_y * pixels_z, fdh2);
          my_fwrite(xHP, sizeof(float), pixels_y * pixels_z, fdhp);
          my_fwrite(xCO, sizeof(float), pixels_y * pixels_z, fdco);
#endif

        }

    }

#ifdef IMAGE_FOOTERS
  write_image_footer(fd, xmin, xmax, ymin, ymax, zmin, zmax);
  write_image_footer(fdtemp, xmin, xmax, ymin, ymax, zmin, zmax);
  write_image_footer(fdmet, xmin, xmax, ymin, ymax, zmin, zmax);
  write_image_footer(fdvel, xmin, xmax, ymin, ymax, zmin, zmax);
#ifdef CHEM_IMAGE
  write_image_footer(fddust, xmin, xmax, ymin, ymax, zmin, zmax);
  write_image_footer(fdh2, xmin, xmax, ymin, ymax, zmin, zmax);
  write_image_footer(fdhp, xmin, xmax, ymin, ymax, zmin, zmax);
  write_image_footer(fdco, xmin, xmax, ymin, ymax, zmin, zmax);
#endif
#endif

  if(fd)
    fclose(fd);
  if(fdtemp)
    fclose(fdtemp);
  if(fdmet)
    fclose(fdmet);
  if(fdvel)
    fclose(fdvel);
  if(fdmag)
    fclose(fdmag);
  if(fdvort)
    fclose(fdvort);
  if(fdphot)
    fclose(fdphot);
  if(fdchem)
    fclose(fdchem);
  if(fdtr)
    fclose(fdtr);
#ifdef CHEM_IMAGE
  if(fddust)
    fclose(fddust);
  if(fdh2)
    fclose(fdh2);
  if(fdhp)
    fclose(fdhp);
  if(fdco)
    fclose(fdco);
#endif

#ifdef CHEM_IMAGE
  myfree(xCO);
  myfree(xHP);
  myfree(xH2);
  myfree(dust);
#endif
  myfree(density_trmc);
  myfree(chem_elements);
  myfree(photon_density);
  myfree(vorticity);
  myfree(Bfield);
  myfree(velocity);
  myfree(met);
  myfree(temp);
  myfree(density);
  myfree(Ray);

  CPU_Step[CPU_MAKEIMAGES] += measure_time();
}
#endif // #if !defined(ONEDIMS) && !defined(TWODIMS)



/** Fills an image slice with information for the rays. Is used by
    both the slice image and "grid of slices" code. The data is
    returned in the output arrays, which should be allocated with
    pixels_x*pixels_y entries (or 3x that in the case of the vector
    quantities).

*/
void fill_slice(ray_data * Ray, int Nray, int gradients_flag,
                int pixels_x, int pixels_y,
                float *density, float *temperature, float *metallicity, float *velocity, float *Bfield, float *vorticity, float *photon_density, float *chem_elements, float *density_trmc
#ifdef CHEM_IMAGE
                , float *dust, float *xH2, float *xHP, float *xCO
#endif
  )
{
  float *dens, *temp, *vel, *mag, *vort, *met, *phot, *chem_elem, *denstrmc;
  int i, sph_idx, pix;
#if defined(RT_ADVECT) || defined(GFM_STELLAR_EVOLUTION)
  int j;
#endif
#ifdef CHEM_IMAGE
  float *dst, *H2, *HP, *CO;
#endif
  double l_dx, l_dy, l_dz;

  dens = mymalloc("dens", pixels_x * pixels_y * sizeof(float));
  temp = mymalloc("temp", pixels_x * pixels_y * sizeof(float));
  met = mymalloc("met", pixels_x * pixels_y * sizeof(float));
  vel = mymalloc("vel", 3 * pixels_x * pixels_y * sizeof(float));
  mag = mymalloc("mag", 3 * pixels_x * pixels_y * sizeof(float));
  vort = mymalloc("vort", 3 * pixels_x * pixels_y * sizeof(float));
  phot = mymalloc("phot", RT_N_DIR * pixels_x * pixels_y * sizeof(float));
  chem_elem = mymalloc("chem_elem", GFM_N_CHEM_ELEMENTS * pixels_x * pixels_y * sizeof(float));
  denstrmc = mymalloc("denstrmc", pixels_x * pixels_y * sizeof(float));
#ifdef CHEM_IMAGE
  dst = mymalloc("dust", pixels_x * pixels_y * sizeof(float));
  H2 = mymalloc("xH2", pixels_x * pixels_y * sizeof(float));
  HP = mymalloc("xHP", pixels_x * pixels_y * sizeof(float));
  CO = mymalloc("xCO", pixels_x * pixels_y * sizeof(float));
#endif

  memset(dens, 0, pixels_x * pixels_y * sizeof(float));
  memset(temp, 0, pixels_x * pixels_y * sizeof(float));
  memset(met, 0, pixels_x * pixels_y * sizeof(float));
  memset(vel, 0, 3 * pixels_x * pixels_y * sizeof(float));
  memset(mag, 0, 3 * pixels_x * pixels_y * sizeof(float));
  memset(vort, 0, 3 * pixels_x * pixels_y * sizeof(float));
  memset(phot, 0, RT_N_DIR * pixels_x * pixels_y * sizeof(float));
  memset(chem_elem, 0, GFM_N_CHEM_ELEMENTS * pixels_x * pixels_y * sizeof(float));
  memset(denstrmc, 0, pixels_x * pixels_y * sizeof(float));
#ifdef CHEM_IMAGE
  memset(dst, 0, pixels_x * pixels_y * sizeof(float));
  memset(H2, 0, pixels_x * pixels_y * sizeof(float));
  memset(HP, 0, pixels_x * pixels_y * sizeof(float));
  memset(CO, 0, pixels_x * pixels_y * sizeof(float));
#endif

  // Now we simply loop over the pixels *we own*
  for(i = 0; i < Nray; ++i)
    {
      sph_idx = Ray[i].index;
      pix = Ray[i].pixel;

      l_dx = Ray[i].pos[0] - SphP[sph_idx].Center[0];
      l_dy = Ray[i].pos[1] - SphP[sph_idx].Center[1];
      l_dz = Ray[i].pos[2] - SphP[sph_idx].Center[2];

#ifdef PERIODIC
#if !defined(REFLECTIVE_X)
      if(l_dx < -boxHalf_X)
        l_dx += boxSize_X;
      if(l_dx > boxHalf_X)
        l_dx -= boxSize_X;
      if(l_dz > boxHalf_Z)
        l_dz -= boxSize_Z;
#endif
#if !defined(REFLECTIVE_Y)
      if(l_dy < -boxHalf_Y)
        l_dy += boxSize_Y;
      if(l_dy > boxHalf_Y)
        l_dy -= boxSize_Y;
      if(l_dz > boxHalf_Z)
        l_dz -= boxSize_Z;
#endif
#endif

      if(gradients_flag == 1)
        dens[pix] = SphP[sph_idx].Density + SphP[sph_idx].Grad.drho[0] * l_dx + SphP[sph_idx].Grad.drho[1] * l_dy + SphP[sph_idx].Grad.drho[2] * l_dz;
      else if(gradients_flag == 0)
        dens[pix] = SphP[sph_idx].Density;
      else
        terminate("gradients_flag != 1 && gradients_flag != 0");

      if(temperature)
        {
#ifdef VORONOI_PROJ_TEMP
          double meanweight = 4.0 / (1. + 3. * HYDROGEN_MASSFRAC) * PROTONMASS;

#if defined(COOLING) && !defined(GRACKLE)
          meanweight = 4.0 / (1. + 3. * HYDROGEN_MASSFRAC + 4. * HYDROGEN_MASSFRAC * SphP[pix].Ne) * PROTONMASS;
#endif
          double temp_in_K = GAMMA_MINUS1 * SphP[sph_idx].Utherm / BOLTZMANN * All.UnitEnergy_in_cgs / All.UnitMass_in_g * meanweight;

#ifdef SGCHEM
          double yn, en, yntot;

          yn = SphP[sph_idx].Density * All.UnitDensity_in_cgs / ((1.0 + 4.0 * ABHE) * PROTONMASS);
          en = SphP[sph_idx].Utherm * SphP[sph_idx].Density * All.UnitEnergy_in_cgs / pow(All.UnitLength_in_cm, 3);
#if CHEMISTRYNETWORK == 1
    yntot = (1.0 + ABHE - SphP[sph_idx].TracAbund[IH2] + SphP[sph_idx].TracAbund[IHP] + SphP[sph_idx].TracAbund[IHEP] + SphP[sph_idx].TracAbund[IHEPP]) * yn;
#else
    yntot = (1.0 + ABHE - SphP[sph_idx].TracAbund[IH2] + SphP[sph_idx].TracAbund[IHP]) * yn;
#endif
	  temp_in_K = (GAMMA - 1.0) * en / (yntot * BOLTZMANN);
#endif /* SGCHEM */

          temp[pix] = temp_in_K;
#else
          temp[pix] = SphP[sph_idx].Utherm;
#endif
        }

#ifdef METALS
      if(metallicity)
        met[pix] = SphP[sph_idx].Metallicity;
#endif

      if(velocity)
        {
          vel[3 * pix + 0] = P[sph_idx].Vel[0] + SphP[sph_idx].Grad.dvel[0][0] * l_dx + SphP[sph_idx].Grad.dvel[0][1] * l_dy + SphP[sph_idx].Grad.dvel[0][2] * l_dz;
          vel[3 * pix + 1] = P[sph_idx].Vel[1] + SphP[sph_idx].Grad.dvel[1][0] * l_dx + SphP[sph_idx].Grad.dvel[1][1] * l_dy + SphP[sph_idx].Grad.dvel[1][2] * l_dz;
          vel[3 * pix + 2] = P[sph_idx].Vel[2] + SphP[sph_idx].Grad.dvel[2][0] * l_dx + SphP[sph_idx].Grad.dvel[2][1] * l_dy + SphP[sph_idx].Grad.dvel[2][2] * l_dz;
        }

#ifdef MHD
      if(Bfield)
        {
          double mag_x = SphP[sph_idx].B[0] + SphP[sph_idx].Grad.dB[0][0] * l_dx + SphP[sph_idx].Grad.dB[0][1] * l_dy + SphP[sph_idx].Grad.dB[0][2] * l_dz;
          double mag_y = SphP[sph_idx].B[1] + SphP[sph_idx].Grad.dB[1][0] * l_dx + SphP[sph_idx].Grad.dB[1][1] * l_dy + SphP[sph_idx].Grad.dB[1][2] * l_dz;
          double mag_z = SphP[sph_idx].B[2] + SphP[sph_idx].Grad.dB[2][0] * l_dx + SphP[sph_idx].Grad.dB[2][1] * l_dy + SphP[sph_idx].Grad.dB[2][2] * l_dz;
          mag[3 * pix + 0] = mag_x;
          mag[3 * pix + 1] = mag_y;
          mag[3 * pix + 2] = mag_z;
        }
#endif

      vort[3 * pix + 0] = SphP[sph_idx].Grad.dvel[2][1] - SphP[sph_idx].Grad.dvel[1][2];
      vort[3 * pix + 1] = SphP[sph_idx].Grad.dvel[0][2] - SphP[sph_idx].Grad.dvel[2][0];
      vort[3 * pix + 2] = SphP[sph_idx].Grad.dvel[1][0] - SphP[sph_idx].Grad.dvel[0][1];

#ifdef RT_ADVECT
      if(photon_density)
        {
          for(j = 0; j < RT_N_DIR; ++j)
            phot[RT_N_DIR * pix + j] = SphP[sph_idx].DensPhot[j];
        }
#endif

#ifdef GFM_STELLAR_EVOLUTION
      if(chem_elements)
        {
          for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
            chem_elem[GFM_N_CHEM_ELEMENTS * pix + j] = SphP[sph_idx].MetalsFraction[j];
        }
#endif

#ifdef TRACER_MC
      if(denstrmc)
        denstrmc[pix] = dens[pix] * get_number_of_tracers(sph_idx) * All.ReferenceTracerMCMass / P[sph_idx].Mass;
#endif

#ifdef CHEM_IMAGE
      dst[pix] += SphP[sph_idx].DustTemp;
      H2[pix] += SphP[sph_idx].TracAbund[IH2];
      HP[pix] += SphP[sph_idx].TracAbund[IHP];
      CO[pix] += SphP[sph_idx].TracAbund[ICO];
#endif

    }

  // Now assemble the densities. Since all points except the ones that
  // were set by the requisite tasks are zero, we just sum the contribution
  MPI_Reduce(dens, density, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(temperature)
    MPI_Reduce(temp, temperature, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(metallicity)
    MPI_Reduce(met, metallicity, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(velocity)
    MPI_Reduce(vel, velocity, 3 * pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(Bfield)
    MPI_Reduce(mag, Bfield, 3 * pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(vorticity)
    MPI_Reduce(vort, vorticity, 3 * pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(photon_density)
    MPI_Reduce(phot, photon_density, RT_N_DIR * pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(chem_elements)
    MPI_Reduce(chem_elem, chem_elements, GFM_N_CHEM_ELEMENTS * pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(density_trmc)
    MPI_Reduce(denstrmc, density_trmc, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#ifdef CHEM_IMAGE
  if(dust)
    MPI_Reduce(dst, dust, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(xH2)
    MPI_Reduce(H2, xH2, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(xHP)
    MPI_Reduce(HP, xHP, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(dust)
    MPI_Reduce(CO, xCO, pixels_x * pixels_y, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#endif

#ifdef CHEM_IMAGE
  myfree(CO);
  myfree(HP);
  myfree(H2);
  myfree(dst);
#endif
  myfree(denstrmc);
  myfree(chem_elem);
  myfree(phot);
  myfree(vort);
  myfree(mag);
  myfree(vel);
  myfree(met);
  myfree(temp);
  myfree(dens);
}



/** If pos is more than a half-box away from ref, wrap it so that they
    are in the same octant. */
void periodic_wrap_point(MyDouble pos[3], MyDouble ref[3])
{
  MyDouble dx, dy, dz;

  dx = pos[0] - ref[0];
  dy = pos[1] - ref[1];
  dz = pos[2] - ref[2];
#ifdef PERIODIC
  if(dx > boxHalf_X)
    pos[0] -= boxSize_X;
  if(dx < -boxHalf_X)
    pos[0] += boxSize_X;
  if(dy > boxHalf_Y)
    pos[1] -= boxSize_Y;
  if(dy < -boxHalf_Y)
    pos[1] += boxSize_Y;
  if(dz > boxHalf_Z)
    pos[2] -= boxSize_Z;
  if(dz < -boxHalf_Z)
    pos[2] += boxSize_Z;
#endif
}




/** Calculates the intersections between a ray at position p0 and
    direction dir and the Voronoi faces across all the connections to
    the cell and returns the element in the DC array representing the
    face the cell will cross. Cell is the index of the primary mesh
    cell the ray is in. Previous is the cell that we came from so, if
    set, it is ignored in the intersection test (making sure we don't
    detect a face we have just entered through due to numerical
    truncation.) The distance to the face is returned in length.

    Note that unlike the funciton in Sunrise that this was back-ported
    from, it DOES use periodic boundary conditions. If an edge to a
    point on the opposite side of the box is encountered, it will be
    wrapped to get the correct intersections.
*/
int find_next_voronoi_cell(tessellation * T, int cell, MyDouble p0[3], double dir[3], int previous, double *length)
{
#ifdef VORONOI_DYNAMIC_UPDATE
  point *DP = T->DP;
  //myassert(DP[cell].index >= 0);
  //myassert(DP[cell].index < T->Ndp);

  MyDouble cell_p[3];
  cell_p[0] = P[cell].Pos[0];
  cell_p[1] = P[cell].Pos[1];
  cell_p[2] = P[cell].Pos[2];

  // if mesh point is across the boundary, wrap it
  periodic_wrap_point(cell_p, p0);

  MyDouble nb_p[3];
  double m[3];
  double c[3];
  double q[3];
  double s;

  // init next to -1, which means it didn't cross the face
  int next = -1;
  // initialize length to huge
  *length = HUGE_VAL;

  int edge = SphP[cell].first_connection;
  int last_edge = SphP[cell].last_connection;
  int iter = 0;

  while(1)
    {
      ++iter;
      const int neighbor = DC[edge].dp_index;

#if (defined(REFLECTIVE_X) && defined(REFLECTIVE_Y) && defined(REFLECTIVE_Z))

      if(DC[edge].image_flags == 1)     //only ignore the edge we entered through if the edge does not cross the box
        {
          // ignore the edge we entered through
          if((DC[edge].index == previous) && (DC[edge].task == ThisTask))
            {
              if(edge == last_edge)
                break;
              edge = DC[edge].next;
              continue;
            }
        }
#else
      myassert((DC[edge].task != ThisTask) || (DC[edge].index != cell));

      // ignore the edge we entered through
      if((DC[edge].index == previous) && (DC[edge].task == ThisTask))
        {

          if(edge == last_edge)
            break;
          edge = DC[edge].next;
          continue;
        }
#endif

      nb_p[0] = DP[neighbor].x;
      nb_p[1] = DP[neighbor].y;
      nb_p[2] = DP[neighbor].z;

#if !(defined(REFLECTIVE_X) && defined(REFLECTIVE_Y) && defined(REFLECTIVE_Z))
      // if neighbor is across the boundary, wrap it
      periodic_wrap_point(nb_p, p0);
#endif

      int i;
      for(i = 0; i < 3; ++i)
        {
          // m is the edge midpoint, which is a point on the face plane
          m[i] = 0.5 * (nb_p[i] + cell_p[i]);
          // c is the vector from point p0 to m.
          c[i] = m[i] - p0[i];
          /* q is the edge vector to the neighboring cell, which is a
             normal vector of the plane */
          q[i] = nb_p[i] - cell_p[i];
        }

      // sanity check: by construction we know that the point is inside
      // the cell, which means that c.q>0. If this is not true, it's
      // because some numerical error has put the point on the other
      // side of the face. In this case we short-circuit the process and
      // say it is on the face, the propagation distance is zero (if
      // it's heading towards the face, ie d.q>0). We must also guard
      // against the case where the ray is perfectly on the face, in
      // which case we will get NaN. In that case we simply ignore the
      // face.
      double cdotq = c[0] * q[0] + c[1] * q[1] + c[2] * q[2];
      double ddotq = dir[0] * q[0] + dir[1] * q[1] + dir[2] * q[2];

      if(cdotq > 0)
        {
          // s = c.nq / d.q
          // This is the standard formula for the intersection between a
          // ray and a plane. s is the point where the ray p0+s*dir
          // intersects the plane which is perpendicular to q and goes
          // through c, i.e. the Voronoi face corresponding to the j-k
          // edge. 
          s = cdotq / ddotq;
        }
      else
        {
          // point on wrong (outside) side of face. Something's up.
          // If this is due to numerical truncation, distance to face
          // should be small. If it's large, it's likely we aren't in
          // the cell we think we are in.

          // double dist_to_face = cdotq / sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);
          //myassert(dist_to_face > -1e-6);

          if(ddotq > 0)
            {
              // it's heading away from the cell, so it must have been
              // supposed to intersect this face
              // set s=0.
              s = 0;
            }
          else
            // it's heading into the cell, so it must have come in on
            // this face. ignore the face
            s = HUGE_VAL;
        }

      if(s >= 0 && s < *length)
        {
          /* The ray intersects the Voronoi face closer to the
             starting point than any previous points. This is our
             current candidate exit face. */
          next = edge;
          *length = s;
        }

      if(edge == last_edge)
        break;

      myassert(edge != DC[edge].next);
      edge = DC[edge].next;
    }

  // set length to the physical length instead of the fractional length
  *length *= sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
  return next;
#else
  myassert(0);
  return 0;
#endif
}

/*! This function works exactly as find_next_voronoi_cell, but in addition it also skips the previous cell
 *  when the ray changes the domain.
 *
 *  \param previous The index of the cell on the original task, where the cell is not imported.
 *  \param task_of_previous The task of the cell 'previous' / the original task.
 */

int find_next_voronoi_cell2(tessellation * T, int cell, MyDouble p0[3], double dir[3], int previous, int task_of_previous, double *length)
{
#ifdef VORONOI_DYNAMIC_UPDATE
  point *DP = T->DP;
  //myassert(DP[cell].index >= 0);
  //myassert(DP[cell].index < T->Ndp);

  MyDouble cell_p[3];
  cell_p[0] = P[cell].Pos[0];
  cell_p[1] = P[cell].Pos[1];
  cell_p[2] = P[cell].Pos[2];

  // if mesh point is across the boundary, wrap it
  periodic_wrap_point(cell_p, p0);

  MyDouble nb_p[3];
  double m[3];
  double c[3];
  double q[3];
  double s;

  // init next to -1, which means it didn't cross the face
  int next = -1;
  // initialize length to huge
  *length = HUGE_VAL;

  int edge = SphP[cell].first_connection;
  int last_edge = SphP[cell].last_connection;
  int iter = 0;

  int ignored = 0;

  while(1)
    {
      ++iter;
      const int neighbor = DC[edge].dp_index;

#if (defined(REFLECTIVE_X) && defined(REFLECTIVE_Y) && defined(REFLECTIVE_Z))

      if(DC[edge].image_flags == 1)     //only ignore the edge we entered through if the edge does not cross the box
        {
          //ignore the edge we entered through
          if((DC[edge].index == previous) && (DC[edge].task == task_of_previous))
            {

              ignored = 1;

              if(edge == last_edge)
                break;
              edge = DC[edge].next;
              continue;
            }
        }
#else


      myassert((DC[edge].task != ThisTask) || (DC[edge].index != cell));


      //ignore the edge we entered through
      if((DC[edge].index == previous) && (DC[edge].task == task_of_previous))
        {
          ignored = 1;

          if(edge == last_edge)
            break;
          edge = DC[edge].next;
          continue;
        }


#endif

      nb_p[0] = DP[neighbor].x;
      nb_p[1] = DP[neighbor].y;
      nb_p[2] = DP[neighbor].z;

#if !(defined(REFLECTIVE_X) && defined(REFLECTIVE_Y) && defined(REFLECTIVE_Z))
      // if neighbor is across the boundary, wrap it
      periodic_wrap_point(nb_p, p0);
#endif

      int i;
      for(i = 0; i < 3; ++i)
        {
          // m is the edge midpoint, which is a point on the face plane
          m[i] = 0.5 * (nb_p[i] + cell_p[i]);
          // c is the vector from point p0 to m.
          c[i] = m[i] - p0[i];
          /* q is the edge vector to the neighboring cell, which is a
             normal vector of the plane */
          q[i] = nb_p[i] - cell_p[i];
        }

      // sanity check: by construction we know that the point is inside
      // the cell, which means that c.q>0. If this is not true, it's
      // because some numerical error has put the point on the other
      // side of the face. In this case we short-circuit the process and
      // say it is on the face, the propagation distance is zero (if
      // it's heading towards the face, ie d.q>0). We must also guard
      // against the case where the ray is perfectly on the face, in
      // which case we will get NaN. In that case we simply ignore the
      // face.
      double cdotq = c[0] * q[0] + c[1] * q[1] + c[2] * q[2];
      double ddotq = dir[0] * q[0] + dir[1] * q[1] + dir[2] * q[2];


      if(cdotq > 0)
        {
          // s = c.nq / d.q
          // This is the standard formula for the intersection between a
          // ray and a plane. s is the point where the ray p0+s*dir
          // intersects the plane which is perpendicular to q and goes
          // through c, i.e. the Voronoi face corresponding to the j-k
          // edge.
          s = cdotq / ddotq;
        }
      else
        {
          // point on wrong (outside) side of face. Something's up.
          // If this is due to numerical truncation, distance to face
          // should be small. If it's large, it's likely we aren't in
          // the cell we think we are in.

          // double dist_to_face = cdotq / sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);
          //myassert(dist_to_face > -1e-6);

          if(ddotq > 0)
            {
              // it's heading away from the cell, so it must have been
              // supposed to intersect this face
              // set s=0.
              s = 0;
            }
          else
            // it's heading into the cell, so it must have come in on
            // this face. ignore the face
            s = HUGE_VAL;
        }

      if(s >= 0 && s < *length)
        {
          /* The ray intersects the Voronoi face closer to the
             starting point than any previous points. This is our
             current candidate exit face. */
          next = edge;
          *length = s;
        }

      if(edge == last_edge)
        break;

      myassert(edge != DC[edge].next);
      edge = DC[edge].next;
    }

#if defined(DEBUG) && !(defined(REFLECTIVE_X) && defined(REFLECTIVE_Y) && defined(REFLECTIVE_Z))
  //assert that the previous cell gets ignored unless the ray just started
  if(!ignored && (!(SphP[cell].Center[0] == p0[0] && SphP[cell].Center[1] == p0[1] && SphP[cell].Center[2] == p0[2])))
    {
      myassert(0);
    }
#endif


  //printf("Tested %d voronoi faces\n", i);

  // set length to the physical length instead of the fractional length
  *length *= sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
  return next;
#else
  myassert(0);
  return 0;
#endif
}


#endif /* VORONOI */
