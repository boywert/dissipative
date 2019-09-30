/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/dmpic/project.c
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
#include <sys/stat.h>
#include <sys/types.h>
#include <gsl/gsl_math.h>
#include <inttypes.h>

#include "../allvars.h"
#include "../proto.h"

/*! \file project.c
 *  \brief code for making dark matter pictures on the fly
 */

#ifdef DMPIC

#include "dmpic.h"
#include "../domain.h"

struct picdata                  /* variables that influence picture geometry */
{
  double TimeTarget;
  double CenterX;
  double CenterY;
  double CenterZ;
  double AvecX, AvecY, AvecZ;
  double BvecX, BvecY, BvecZ;
  double Depth;
  double EyeAngle;
  double LengthX;
  double LengthY;
  double CvecX, CvecY, CvecZ;
}
 *Pic;


static float *LocMassMap, *LocDensMap, *LocVelDispMap, *LocVelZMap;
static float *MassMap, *DensMap, *VelDispMap, *VelZMap;


void dmpic_make_image(void)
{
  int i;

  CPU_Step[CPU_MISC] += measure_time();

  mpi_printf("DMPIC: Begin to compute local densities...  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

  ngb_treefree();
  domain_free();

  domain_Decomposition();

  ngb_treeallocate();
  ngb_treebuild(NumGas);


  /* this structure will hold auxiliary information for each particle, needed only during group finding */
  PS = (struct subfind_data *) mymalloc_movable(&PS, "PS", All.MaxPart * sizeof(struct subfind_data));
  memset(PS, 0, NumPart * sizeof(struct subfind_data));

  /* First, we save the original location of the particles, in order to be able to revert to this layout later on */
  for(i = 0; i < NumPart; i++)
    {
      PS[i].OriginTask = ThisTask;
      PS[i].OriginIndex = i;
    }

  force_treeallocate(NumPart, All.MaxPart);
  force_treebuild(NumPart, 0, 0);

  if(All.PicCurrentNr == 0)
    subfind_density_hsml_guess();
  else
    for(i = 0; i < NumPart; i++)
      PS[i].Hsml = P[i].DM_Hsml;

  /* calculate velocity dispersion etc. */

  subfind_density(0);

  for(i = 0; i < NumPart; i++)
    {
      P[i].DM_Rho = PS[i].Density;
      P[i].DM_Hsml = PS[i].Hsml;
      P[i].DM_VelDisp = PS[i].SubfindVelDisp;
    }

  myfree(Father);
  myfree(Nextnode);
#ifdef BLACK_HOLES
  myfree(Tree_AuxBH_Points);
#endif
  myfree(Tree_Points);
  force_treefree();

  myfree(PS);

  CPU_Step[CPU_FOF] += measure_time();


  /* now create the actual image */
  dmpic_project_and_smooth(All.PicCurrentNr++);
}






integertime dmpic_get_next_pictime(void)
{
  integertime ti = 0;

  if(All.PicCurrentNr < All.PicCount)
    {
      if(All.ComovingIntegrationOn)
        ti = (integertime) (log(Pic[All.PicCurrentNr].TimeTarget / All.TimeBegin) / All.Timebase_interval);
      else
        ti = (integertime) ((Pic[All.PicCurrentNr].TimeTarget - All.TimeBegin) / All.Timebase_interval);
    }

  return ti;
}


void dmpic_init(void)
{
  FILE *fd;
  int rep, nr;
  char buf[2048];

  for(rep = 0; rep < 2; rep++)
    {
      if(rep == 1)
        Pic = mymalloc("Pic", All.PicCount * sizeof(struct picdata));

      All.PicCount = 0;

      if(!(fd = fopen(All.PicList, "r")))
        terminate("can't read picture list in file '%s'\n", All.PicList);

      while(1)
        {
          if(fgets(buf, 2000, fd) != buf)
            break;

          double pictime, cx, cy, cz, ax, ay, az, bx, by, bz, depth, angle;

          int count = sscanf(buf, " %d %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg ",
                             &nr, &pictime, &cx, &cy, &cz, &ax, &ay, &az, &bx, &by, &bz, &depth, &angle);

          if(count == 13)
            {
              if(rep == 1)
                {
                  Pic[All.PicCount].TimeTarget = pictime;
                  Pic[All.PicCount].CenterX = cx;
                  Pic[All.PicCount].CenterY = cy;
                  Pic[All.PicCount].CenterZ = cz;
                  Pic[All.PicCount].AvecX = ax;
                  Pic[All.PicCount].AvecY = ay;
                  Pic[All.PicCount].AvecZ = az;
                  Pic[All.PicCount].BvecX = bx;
                  Pic[All.PicCount].BvecY = by;
                  Pic[All.PicCount].BvecZ = bz;
                  Pic[All.PicCount].Depth = depth;
                  Pic[All.PicCount].EyeAngle = angle;
                }

              All.PicCount++;
            }
        }

      fclose(fd);
    }

  All.PicCurrentNr = 0;

  mpi_printf("\nDMPIC: found %d image descriptions.\n", All.PicCount);
}


void dmpic_geometry(int nr)
{
  double norm;

  norm = sqrt(Pic[nr].AvecX * Pic[nr].AvecX + Pic[nr].AvecY * Pic[nr].AvecY + Pic[nr].AvecZ * Pic[nr].AvecZ);

  Pic[nr].AvecX /= norm;
  Pic[nr].AvecY /= norm;
  Pic[nr].AvecZ /= norm;

  Pic[nr].CvecX = Pic[nr].AvecY * Pic[nr].BvecZ - Pic[nr].AvecZ * Pic[nr].BvecY;
  Pic[nr].CvecY = Pic[nr].AvecZ * Pic[nr].BvecX - Pic[nr].AvecX * Pic[nr].BvecZ;
  Pic[nr].CvecZ = Pic[nr].AvecX * Pic[nr].BvecY - Pic[nr].AvecY * Pic[nr].BvecX;

  norm = sqrt(Pic[nr].CvecX * Pic[nr].CvecX + Pic[nr].CvecY * Pic[nr].CvecY + Pic[nr].CvecZ * Pic[nr].CvecZ);

  Pic[nr].CvecX /= norm;
  Pic[nr].CvecY /= norm;
  Pic[nr].CvecZ /= norm;

  Pic[nr].BvecX = Pic[nr].CvecY * Pic[nr].AvecZ - Pic[nr].CvecZ * Pic[nr].AvecY;
  Pic[nr].BvecY = Pic[nr].CvecZ * Pic[nr].AvecX - Pic[nr].CvecX * Pic[nr].AvecZ;
  Pic[nr].BvecZ = Pic[nr].CvecX * Pic[nr].AvecY - Pic[nr].CvecY * Pic[nr].AvecX;

  Pic[nr].LengthX = 2 * Pic[nr].Depth * tan(0.5 * Pic[nr].EyeAngle * M_PI / 180.0);
  Pic[nr].LengthY = Pic[nr].LengthX * ((double) All.PixelsY) / All.PixelsX;

  mpi_printf("\nDMPIC: Basis vectors: nr=%d\n", nr);
  mpi_printf("DMPIC: A = ( %10.6f | %10.6f | %10.6f )\n", Pic[nr].AvecX, Pic[nr].AvecY, Pic[nr].AvecZ);
  mpi_printf("DMPIC: B = ( %10.6f | %10.6f | %10.6f )\n", Pic[nr].BvecX, Pic[nr].BvecY, Pic[nr].BvecZ);
  mpi_printf("DMPIC: C = ( %10.6f | %10.6f | %10.6f )\n", Pic[nr].CvecX, Pic[nr].CvecY, Pic[nr].CvecZ);
  mpi_printf("\n");
  mpi_printf("DMPIC: LengthX=%g  LengthY=%g\n", Pic[nr].LengthX, Pic[nr].LengthY);
}


void dmpic_project_and_smooth(int nr)
{
  int i, j, n, flag, count = 0;
  double r, r2, h, h2, hinv, wk, u;
  double sum, hmin, x, y, z, xx, yy, zz, xxx, yyy, hmax, xt, yt, zt, xp, yp;
  double vz, pixelsizeX, pixelsizeY, mass;
  int dx, dy, nx, ny, ii, jj, kk;


  /* allocate local picture data */
  LocMassMap = mymalloc_clear("LocMassMap", All.PixelsX * All.PixelsY * sizeof(float));
  LocDensMap = mymalloc_clear("LocDensMap", All.PixelsX * All.PixelsY * sizeof(float));
  LocVelDispMap = mymalloc_clear("LocVelDispMap", All.PixelsX * All.PixelsY * sizeof(float));
  LocVelZMap = mymalloc_clear("LocVelZMap", All.PixelsX * All.PixelsY * sizeof(float));

  dmpic_geometry(nr);

  pixelsizeX = Pic[nr].LengthX / All.PixelsX;
  pixelsizeY = Pic[nr].LengthY / All.PixelsY;

  if(pixelsizeX < pixelsizeY)
    hmin = 1.001 * pixelsizeX / 2;
  else
    hmin = 1.001 * pixelsizeY / 2;

  hmax = 100.0 * pixelsizeX;



  double ascale = Pic[nr].TimeTarget;
  double box = All.BoxSize * ascale;
  double zoffset, zmin, zmax, ma;
  int boxc;

  if(box < 4 * Pic[nr].Depth)
    {
      zoffset = 2 * Pic[nr].Depth;
      zmin = 2 * Pic[nr].Depth - 0.5 * box;
      zmax = 2 * Pic[nr].Depth + 0.5 * box;

      if(Pic[nr].LengthX > Pic[nr].LengthY)
        ma = Pic[nr].LengthX;
      else
        ma = Pic[nr].LengthY;

      ma *= (2 * Pic[nr].Depth + 0.5 * box) / (2 * Pic[nr].Depth);

      boxc = ceil((ma - 0.5 * box) / box);
    }
  else if(box < 8 * Pic[nr].Depth)
    {
      zoffset = 4 * Pic[nr].Depth - 0.5 * box;
      zmin = 0;
      zmax = 4 * Pic[nr].Depth;

      if(Pic[nr].LengthX > Pic[nr].LengthY)
        ma = Pic[nr].LengthX;
      else
        ma = Pic[nr].LengthY;

      ma *= 2;

      if(4 * Pic[nr].Depth > ma)
        ma = Pic[nr].Depth;

      boxc = ceil((ma - 0.5 * box) / box);
    }
  else
    {
      zoffset = 4 * Pic[nr].Depth - 0.5 * box;
      zmin = 0;
      zmax = 4 * Pic[nr].Depth;

      if(Pic[nr].LengthX > Pic[nr].LengthY)
        ma = Pic[nr].LengthX;
      else
        ma = Pic[nr].LengthY;

      ma *= 2;

      if(4 * Pic[nr].Depth > ma)
        ma = Pic[nr].Depth;

      boxc = ceil((ma - 0.5 * box) / box);
    }

  mpi_printf("DMPIC: Depth=%g  boxc=%d\n", Pic[nr].Depth, boxc);


  for(n = 0; n < NumPart; n++)
    {
      double xtmp, ytmp, ztmp;

      vz = Pic[nr].CvecX * P[n].Vel[0] + Pic[nr].CvecY * P[n].Vel[1] + Pic[nr].CvecZ * P[n].Vel[2];

      x = GRAVITY_NEAREST_X(P[n].Pos[0] - Pic[nr].CenterX);
      y = GRAVITY_NEAREST_Y(P[n].Pos[1] - Pic[nr].CenterY);
      z = GRAVITY_NEAREST_Z(P[n].Pos[2] - Pic[nr].CenterZ);

      x *= ascale;
      y *= ascale;
      z *= ascale;

      for(ii = -boxc; ii <= boxc; ii++)
        for(jj = -boxc; jj <= boxc; jj++)
          for(kk = -boxc; kk <= boxc; kk++)
            {
              flag = 0;

              xx = x + ii * box;
              yy = y + jj * box;
              zz = z + kk * box;

              zt = Pic[nr].CvecX * xx + Pic[nr].CvecY * yy + Pic[nr].CvecZ * zz;

              zt += zoffset;


              if(zt > zmin && zt < zmax)
                {
                  xt = Pic[nr].AvecX * xx + Pic[nr].AvecY * yy + Pic[nr].AvecZ * zz;
                  yt = Pic[nr].BvecX * xx + Pic[nr].BvecY * yy + Pic[nr].BvecZ * zz;

                  xt *= Pic[nr].Depth / zt;
                  yt *= Pic[nr].Depth / zt;

                  h = ascale * P[n].DM_Hsml * Pic[nr].Depth / zt;

                  if(h < hmin)
                    h = hmin;
                  if(h > hmax)
                    h = hmax;

                  if(fabs(xt) - h > 0.5 * Pic[nr].LengthX)
                    continue;

                  if(fabs(yt) - h > 0.5 * Pic[nr].LengthY)
                    continue;

                  double dd;

                  if(zt < 0.1 * Pic[nr].Depth)
                    dd = 0.1 * Pic[nr].Depth;
                  else
                    dd = zt;

                  mass = P[n].Mass / (dd * dd) * erfc(0.5 * zt / Pic[nr].Depth) * erfc(0.1 * Pic[nr].Depth / zt);

                  h2 = h * h;
                  hinv = 1.0 / h;
                  nx = h / pixelsizeX + 1;
                  ny = h / pixelsizeY + 1;

                  /* xp,yp central pixel of region covered by the particle on the mesh */

                  xp = (floor(xt / pixelsizeX) + 0.5) * pixelsizeX;
                  yp = (floor(yt / pixelsizeY) + 0.5) * pixelsizeY;

                  /* determine kernel normalizaton */

                  sum = 0;

                  for(dx = -nx; dx <= nx; dx++)
                    for(dy = -ny; dy <= ny; dy++)
                      {
                        xx = xp + dx * pixelsizeX - xt;
                        yy = yp + dy * pixelsizeY - yt;
                        r2 = xx * xx + yy * yy;

                        if(r2 < h2)
                          {
                            r = sqrt(r2);
                            u = r * hinv;

                            if(u < 0.5)
                              wk = (2.546479089470 + 15.278874536822 * (u - 1) * u * u);
                            else
                              wk = 5.092958178941 * (1.0 - u) * (1.0 - u) * (1.0 - u);

                            sum += wk;
                          }
                      }

                  if(sum < 1.0e-10)
                    continue;

                  for(dx = -nx; dx <= nx; dx++)
                    for(dy = -ny; dy <= ny; dy++)
                      {
                        xxx = xp + dx * pixelsizeX + 0.5 * Pic[nr].LengthX;
                        yyy = yp + dy * pixelsizeY + 0.5 * Pic[nr].LengthY;

                        if(xxx >= 0 && xx < Pic[nr].LengthX && yyy >= 0 && yyy < Pic[nr].LengthY)
                          {
                            i = xxx / pixelsizeX;
                            j = yyy / pixelsizeY;

                            if(i >= 0 && i < All.PixelsX)
                              if(j >= 0 && j < All.PixelsY)
                                {
                                  xx = xp + dx * pixelsizeX - xt;
                                  yy = yp + dy * pixelsizeY - yt;
                                  r2 = xx * xx + yy * yy;

                                  if(r2 < h2)
                                    {
                                      r = sqrt(r2);
                                      u = r * hinv;

                                      if(u < 0.5)
                                        wk = (2.546479089470 + 15.278874536822 * (u - 1) * u * u);
                                      else
                                        wk = 5.092958178941 * (1.0 - u) * (1.0 - u) * (1.0 - u);

                                      double fac = wk / sum;
                                      LocMassMap[i * All.PixelsY + j] += mass * fac;
                                      LocDensMap[i * All.PixelsY + j] += mass * P[n].DM_Rho * fac;
                                      LocVelDispMap[i * All.PixelsY + j] += mass * P[n].DM_Rho * P[n].DM_VelDisp * fac;
                                      LocVelZMap[i * All.PixelsY + j] += mass * P[n].DM_Rho * vz * fac;

                                      flag = 1;
                                    }
                                }
                          }
                      }
                }

              if(flag)
                count++;
            }
    }


  /* allocate global picture data */
  MassMap = mymalloc_clear("MassMap", All.PixelsX * All.PixelsY * sizeof(float));
  DensMap = mymalloc_clear("DensMap", All.PixelsX * All.PixelsY * sizeof(float));
  VelDispMap = mymalloc_clear("VelDispMap", All.PixelsX * All.PixelsY * sizeof(float));
  VelZMap = mymalloc_clear("VelZMap", All.PixelsX * All.PixelsY * sizeof(float));

  /* consolidate across processors */
  MPI_Reduce(LocMassMap, MassMap, All.PixelsX * All.PixelsY, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(LocDensMap, DensMap, All.PixelsX * All.PixelsY, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(LocVelDispMap, VelDispMap, All.PixelsX * All.PixelsY, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(LocVelZMap, VelZMap, All.PixelsX * All.PixelsY, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);


  dmpic_save_picture(nr);

  myfree(VelZMap);
  myfree(VelDispMap);
  myfree(DensMap);
  myfree(MassMap);


  myfree(LocVelZMap);
  myfree(LocVelDispMap);
  myfree(LocDensMap);
  myfree(LocMassMap);


  long long countall;
  sumup_large_ints(1, &count, &countall);
  mpi_printf("DMPIC: %ld particles contributed to picture\n", countall);
}


void dmpic_save_minmax(float *map, FILE * fd)
{
  int i;
  float max, min;
  min = max = map[0];


  for(i = 0; i < All.PixelsX * All.PixelsY; i++)
    {
      if(map[i] > max)
        max = map[i];

      if(map[i] < min)
        min = map[i];
    }

  fwrite(&min, sizeof(float), 1, fd);
  fwrite(&max, sizeof(float), 1, fd);
}


void dmpic_save_picture(int nr)
{
  FILE *fd;
  char buf[1024];

  if(ThisTask == 0)
    {
      mkdir(All.PicDataDir, 02755);
      sprintf(buf, "%s/%s_%04d.dat", All.PicDataDir, All.PicDataName, nr);

      if(!(fd = fopen(buf, "w")))
        terminate("can't open file `%s'\n", buf);

      fwrite(&All.PixelsX, sizeof(int), 1, fd);
      fwrite(&All.PixelsY, sizeof(int), 1, fd);

      dmpic_save_minmax(MassMap, fd);
      dmpic_save_minmax(VelDispMap, fd);
      dmpic_save_minmax(DensMap, fd);
      dmpic_save_minmax(VelZMap, fd);

      fwrite(MassMap, All.PixelsX * All.PixelsY, sizeof(float), fd);
      fwrite(VelDispMap, All.PixelsX * All.PixelsY, sizeof(float), fd);
      fwrite(DensMap, All.PixelsX * All.PixelsY, sizeof(float), fd);
      fwrite(VelZMap, All.PixelsX * All.PixelsY, sizeof(float), fd);
      fclose(fd);

      mpi_printf("DMPIC: frame-file '%s' written.\n", buf);
    }
}









#endif
