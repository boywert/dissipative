/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/voronoi_check.c
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

#include "allvars.h"
#include "proto.h"
#include "voronoi.h"

#ifdef VORONOI

#if !defined(TWODIMS) && !defined(ONEDIMS)      /* three-dimensional test code */

int points_compare(const void *a, const void *b);


void check_for_min_distance(tessellation * T)
{
  point *DP = T->DP;
  int i, j;
  double r2, r2min;
  char msg[200];

  for(i = 0, r2min = 1.0e30; i < T->Ndp; i++)
    {
      printf("i=%d\n", i);

      for(j = i + 1; j < T->Ndp; j++)
        {
          r2 = (DP[i].x - DP[j].x) * (DP[i].x - DP[j].x) + (DP[i].y - DP[j].y) * (DP[i].y - DP[j].y) + (DP[i].z - DP[j].z) * (DP[i].z - DP[j].z);
          if(r2 < r2min)
            r2min = r2;

          if(r2min == 0)
            {
              sprintf(msg, "i=%d j=%d equal.  DP[i].index=%d DP[j].index=%d\n", i, j, DP[i].index, DP[j].index);
            terminate(msg)}
        }
    }

  printf("min distance=%g\n", sqrt(r2min));
}


void check_links(tessellation * T)
{
  tetra *DT = T->DT;
  int i, j, s, c, flag = 0;
  int pl[3], pr[3];
  char msg[200];


  for(i = 0; i < T->Ndt; i++)
    {
      if(DT[i].t[0] < 0)        /* deleted ? */
        continue;

      for(j = 0; j < 4; j++)
        {
          if(DT[DT[i].t[j]].t[DT[i].s[j]] != i)
            {
              printf("LINK for tetra=%d j=%d DT[i].s[j]=%d incorrect %d\n", i, j, DT[i].s[j], (int) (DT[DT[i].t[j]].t[DT[i].s[j]]));
            }
        }

      for(j = 0; j < 4; j++)
        {
          for(s = 0, c = 0; s < 4; s++)
            if(s != j)
              pl[c++] = DT[i].p[s];

          for(s = 0, c = 0; s < 4; s++)
            if(s != DT[i].s[j])
              pr[c++] = DT[DT[i].t[j]].p[s];

          /* sort the points */

          mysort(&pl[0], 3, sizeof(int), points_compare);
          mysort(&pr[0], 3, sizeof(int), points_compare);

          for(s = 0; s < 3; s++)
            {
              if(pl[s] != pr[s])
                {
                  sprintf(msg, "LINK for i=%d j=%d incorrect. points of triangles don't match up s=%d\n", i, j, s);
                  flag = 1;
                }
            }

          if(flag)
            terminate(msg);
        }
    }

  printf("links ok\n");
}


void check_orientations(tessellation * T)
{
  tetra *DT = T->DT;
  point *DP = T->DP;
  int i, ivol;
  double vol, volmin = 1.0e30;
  char msg[200];

  for(i = 0; i < T->Ndt; i++)
    {
      tetra *t = &DT[i];

      point *p0 = &DP[t->p[0]];
      point *p1 = &DP[t->p[1]];
      point *p2 = &DP[t->p[2]];
      point *p3 = &DP[t->p[3]];

      if(t->t[0] < 0)           /* deleted ? */
        continue;

      if(isInfinity(p0) || isInfinity(p1) || isInfinity(p2) || isInfinity(p3))
        continue;

      vol = calculate_tetra_volume(p0, p1, p2, p3);
      ivol = Orient3d_Exact(p0, p1, p2, p3);

      if(ivol <= 0)
        {
          sprintf(msg, "Tetra %d is NEGATIVE (%d %d %d %d) oriented or FLAT: ivol=%d vol=%g\n", i, (int) (t->p[0]), (int) (t->p[1]), (int) (t->p[2]), (int) (t->p[3]), ivol, vol);
          terminate(msg);
        }

      if(vol < volmin)
        volmin = vol;
    }

  printf("orientations ok, volmin=%g\n", volmin);
}





void check_tetras(tessellation * T, int npoints)
{
  tetra *DT = T->DT;
  point *DP = T->DP;
  int i, j, res, res_exact;
  char msg[200];

  for(i = 0; i < T->Ndt; i++)
    {
      if((i % 100) == 0)
        printf("check tetra i=%d/%d\n", i, T->Ndt);

      tetra *t = &DT[i];

      point *p0 = &DP[t->p[0]];
      point *p1 = &DP[t->p[1]];
      point *p2 = &DP[t->p[2]];
      point *p3 = &DP[t->p[3]];

      if(t->t[0] < 0)           /* deleted ? */
        continue;

      if(isInfinity(p0) || isInfinity(p1) || isInfinity(p2) || isInfinity(p3))
        continue;

      if(test_tetra_orientation(p0, p1, p2, p3) > 0)
        {
          /*
             printf("Tetra %d is positive oriented\n", i);
           */
        }
      else
        {
          sprintf(msg, "Tetra %d is NEGATIVE oriented\n", i);
          terminate(msg);
        }


      for(j = 0; j < npoints; j++)
        {
          if(t->p[0] != j)
            if(t->p[1] != j)
              if(t->p[2] != j)
                if(t->p[3] != j)
                  {
                    res = InSphere_Errorbound(p0, p1, p2, p3, &DP[j]);

                    if(res >= 0)
                      {
                        res_exact = InSphere_Exact(p0, p1, p2, p3, &DP[j]);

                        if(res_exact > 0)
                          {
                            sprintf(msg, "ERROR tetra=%d: point=%d  in tetra with edges=%d|%d|%d|%d   res=%d|%d\n",
                                    i, j, (int) (t->p[0]), (int) (t->p[1]), (int) (t->p[2]), (int) (t->p[3]), res, res_exact);
                            terminate(msg);
                          }
                      }
                  }
        }
    }

  printf("Tetrahedra OK\n");
}


int points_compare(const void *a, const void *b)
{
  if(*((int *) a) < *((int *) b))
    return -1;

  if(*((int *) a) > *((int *) b))
    return +1;

  return 0;
}

#endif



#ifdef TWODIMS                  /* two-dimensional test code */

void check_triangles(tessellation * T, int npoints)
{
  int i, j, res, res_exact;
  char msg[200];

  tetra *DT = T->DT;

  for(i = 0; i < T->Ndt; i++)
    {
      if(DT[i].p[0] == DPinfinity)
        continue;
      if(DT[i].p[1] == DPinfinity)
        continue;
      if(DT[i].p[2] == DPinfinity)
        continue;

      if(Orient2d_Exact(T, DT[i].p[0], DT[i].p[1], DT[i].p[2]) != 1)
        {
          sprintf(msg, "Triangle %d is NEGATIVE oriented or FLAT\n", i);
          terminate(msg);
        }

      for(j = 0; j < npoints; j++)
        {
          if(DT[i].p[0] != j)
            if(DT[i].p[1] != j)
              if(DT[i].p[2] != j)
                {
                  res = InCircle_Quick(T, DT[i].p[0], DT[i].p[1], DT[i].p[2], j);

                  if(res > 0)
                    {
                      res_exact = InCircle_Exact(T, DT[i].p[0], DT[i].p[1], DT[i].p[2], j);

                      if(res_exact > 0)
                        {
                          sprintf(msg, "ERROR: point=%d lies in triangle=%d with edges=%d|%d|%d   res=%d|%d\n", j, i, (int) (DT[i].p[0]), (int) (DT[i].p[1]), (int) (DT[i].p[2]), res, res_exact);
                          terminate(msg);
                        }
                    }
                }
        }
    }

  printf("triangles ok\n");
}

void check_orientations(tessellation * T)
{
  int i, ivol;
  double vol, volmin = 1.0e30;
  char msg[200];

  tetra *DT = T->DT;

  for(i = 0; i < T->Ndt; i++)
    {
      if(DT[i].p[0] == DPinfinity)
        continue;
      if(DT[i].p[1] == DPinfinity)
        continue;
      if(DT[i].p[2] == DPinfinity)
        continue;

      vol = test_triangle_orientation(T, DT[i].p[0], DT[i].p[1], DT[i].p[2]);
      ivol = Orient2d_Exact(T, DT[i].p[0], DT[i].p[1], DT[i].p[2]);

      if(ivol <= 0)
        {
          double vol2 = Orient2d_Quick(T, DT[i].p[0], DT[i].p[1], DT[i].p[2]);

          sprintf(msg, "Triangle %d is NEGATIVE (%d %d %d) oriented or FLAT: ivol=%d vol=%g|%g\n", i, (int) (DT[i].p[0]), (int) (DT[i].p[1]), (int) (DT[i].p[2]), ivol, vol, vol2);
          terminate(msg);
        }

      if(vol < volmin)
        volmin = vol;
    }

  printf("orientations ok, volmin=%g\n", volmin);
}



void check_links(tessellation * T)
{
  int i, j;
  char msg[200];

  tetra *DT = T->DT;

  for(i = 0; i < T->Ndt; i++)
    {
      for(j = 0; j < 3; j++)
        {
          if(DT[DT[i].t[j]].t[DT[i].s[j]] != i)
            {
              sprintf(msg, "LINK for i=%d j=%d  incorrect\n", i, j);
              terminate(msg);
            }
        }
    }
}

#endif

#endif /* VORONOI */
