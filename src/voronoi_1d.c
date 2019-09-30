/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/voronoi_1d.c
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
#include <gmp.h>

#include "allvars.h"
#include "proto.h"
#include "voronoi.h"


#if defined (VORONOI) && defined (ONEDIMS) && !defined (ONEDIMS_SPHERICAL)      /* will only be compiled in 1D case */

void write_voronoi_mesh(tessellation * T, char *fname, int writeTask, int lastTask)
{
  terminate("You stupid idiot should not do that!");
}

void initialize_and_create_first_tetra(tessellation * T)
{
  char msg[200];

  if(NTask > 1)
    {
      mpi_printf("1D code works only for 1 CPU\n");
      endrun();
    }

  T->MaxNdp = NumGas + 4;
  T->MaxNdt = 4 + T->MaxNdp * 2;
  T->MaxNvf = T->MaxNdt;

  if(NumGas == 0)
    {
      sprintf(msg, "NumGas=%d on Task=%d, but need at least one particle!\n", NumGas, ThisTask);
      terminate(msg);
    }

  T->Ndp = 0;
  T->Nvf = 0;
  T->Ndt = 0;


  T->VF = mymalloc_movable(&T->VF, "VF", T->MaxNvf * sizeof(face));

  T->DP = mymalloc_movable(&T->DP, "DP", (T->MaxNdp + 5) * sizeof(point));
  T->DP += 5;

  T->DT = mymalloc_movable(&T->DT, "DT", T->MaxNdt * sizeof(tetra));
}

void compute_circumcircles(tessellation * T)
{
}

void set_integers_for_point(tessellation * T, int pp)
{
}


int insert_point(tessellation * T, int pp, int ttstart) /* returns a triangle that (currently) contains the point p */
{
  return 0;
}

int voronoi_ghost_search(tessellation * T)
{
  return voronoi_ghost_search_alternative(T);
}

int count_undecided_tetras(tessellation * T)
{
  return 0;
}

int voronoi_ghost_search_alternative(tessellation * T)
{
  double xl, xr;
  int index_l, index_r;


#if defined(REFLECTIVE_X)
  xl = -P[0].Pos[0];
  index_l = 0;

  xr = boxSize_X + (boxSize_X - P[NumGas - 1].Pos[0]);
  index_r = NumGas - 1;
#else
  xl = P[NumGas - 1].Pos[0] - boxSize_X;
  index_l = NumGas - 1;

  xr = P[0].Pos[0] + boxSize_X;
  index_r = 0;
#endif

  point *DP = T->DP;

  DP[-1].x = xl;
  DP[-1].y = 0;
  DP[-1].z = 0;
  DP[-1].task = ThisTask;
  DP[-1].ID = P[index_l].ID;
  DP[-1].index = index_l + NumGas;      /* this is a mirrored local point */
#if defined(REFLECTIVE_X)
  DP[-1].image_flags = REFL_X_FLAGS;
#if (REFLECTIVE_X == 2)
  DP[-1].image_flags |= OUTFLOW_X;
#endif
#endif
  DP[NumGas].x = xr;
  DP[NumGas].y = 0;
  DP[NumGas].z = 0;
  DP[NumGas].task = ThisTask;
  DP[NumGas].ID = P[index_r].ID;
  DP[NumGas].index = index_r + NumGas;  /* this is a mirrored local point */
#if defined(REFLECTIVE_X)
  DP[NumGas].image_flags = REFL_X_FLAGS;
#if (REFLECTIVE_X == 2)
  DP[NumGas].image_flags |= OUTFLOW_X;
#endif
#endif
  return 0;
}

void compute_voronoi_faces_and_volumes(void)
{
  int i;

  tessellation *T = &Mesh;

  T->Nvf = 0;
  point *DP = T->DP;
  face *VF = T->VF;

  for(i = -1; i < NumGas; i++)
    {
      VF[T->Nvf].p1 = i;
      VF[T->Nvf].p2 = i + 1;

      VF[T->Nvf].cx = 0.5 * (DP[i].x + DP[i + 1].x);

      VF[T->Nvf].cy = 0;
      VF[T->Nvf].cz = 0;
      VF[T->Nvf].area = 1;

      T->Nvf++;
    }

  for(i = 0; i < NumGas; i++)
    {
      SphP[i].Volume = VF[i + 1].cx - VF[i].cx;
      SphP[i].Center[0] = 0.5 * (VF[i + 1].cx + VF[i].cx);
      SphP[i].Center[1] = 0;
      SphP[i].Center[2] = 0;

      SphP[i].SurfaceArea = 2.;

      /*
         printf("i=%d x=%g volume=%g l=%g r=%g\n", i, P[i].Pos[0], SphP[i].Volume, VF[i].cx, VF[i + 1].cx);
       */
    }
}







static struct voronoi_1D_data
{
  double x;
  int index;
}
 *mp;

static int *Id;

void voronoi_1D_order(void)
{
  int i;

  mpi_printf("begin 1D order...\n");

  if(NumGas)
    {
      mp = (struct voronoi_1D_data *) mymalloc("mp", sizeof(struct voronoi_1D_data) * NumGas);
      Id = (int *) mymalloc("Id", sizeof(int) * NumGas);

      for(i = 0; i < NumGas; i++)
        {
          mp[i].index = i;
          mp[i].x = P[i].Pos[0];
        }

      mysort(mp, NumGas, sizeof(struct voronoi_1D_data), voronoi_1D_compare_key);


      for(i = 0; i < NumGas; i++)
        Id[mp[i].index] = i;

      voronoi_1D_reorder_gas();

      myfree(Id);
      myfree(mp);
    }

  mpi_printf("1D order done.\n");
}


int voronoi_1D_compare_key(const void *a, const void *b)
{
  if(((struct voronoi_1D_data *) a)->x < (((struct voronoi_1D_data *) b)->x))
    return -1;

  if(((struct voronoi_1D_data *) a)->x > (((struct voronoi_1D_data *) b)->x))
    return +1;

  return 0;
}

void voronoi_1D_reorder_gas(void)
{
  int i;
  struct particle_data Psave, Psource;
  struct sph_particle_data SphPsave, SphPsource;
  int idsource, idsave, dest;

  for(i = 0; i < NumGas; i++)
    {
      if(Id[i] != i)
        {
          Psource = P[i];
          SphPsource = SphP[i];

          idsource = Id[i];
          dest = Id[i];

          do
            {
              Psave = P[dest];
              SphPsave = SphP[dest];
              idsave = Id[dest];

              P[dest] = Psource;
              SphP[dest] = SphPsource;
              Id[dest] = idsource;

              if(dest == i)
                break;

              Psource = Psave;
              SphPsource = SphPsave;
              idsource = idsave;

              dest = idsource;
            }
          while(1);
        }
    }
}




#endif
