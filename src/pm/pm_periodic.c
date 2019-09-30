/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        pm/pm_periodic.c
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

#include "../allvars.h"
#include "../proto.h"

/*! \file  pm_periodic.c
 *  \brief Routines for periodic PM-force computation
 *
 *
 * These routines support two different strategies for doing the particle data exchange to assemble the density field
 * and to read out the forces and potentials:
 *
 * The default scheme sends the particle positions to the target slabs, and bins them there. This works usually well for
 * homogenuously loaded boxes, but can be problematic for zoom-in runs. In the latter case,  PM_ZOOM_OPTIMIZED can be
 * activated, where the data is binned on the originating processor followed by assembly of the binned density field.
 *
 * In addition, the routines can be either used with a slab-based FFT (as is traditionally done in FFTW), or with a
 * column-based FFT. The latter requires more communication and is hence usually slower than the slab-based one.
 * But if the number of MPI ranks exceeds the number of cells per dimension, then the column-based one can still scale
 * and offers a balanced memory consumption, whereas this is not the case for the slab-based approach. To select the
 * column-based FFT, the switch FFT_COLUMN_BASED can be activated.
 *
 * The switches PM_ZOOM_OPTIMIZED and FFT_COLUMN_BASED may also be combined, such that there are 4 main modes of how the
 * PM routines may operate.
 *
 * It is also possible to use non-cubical boxes, by means of setting one or several of the LONG_X, LONG_Y, and LONG_Z
 * options in the config file. The values need to be integers, and then BoxSize is stretched by that factor in the
 * corresponding dimension.
 *
 * Finally, one may also use the TreePM routine for simulations where gravity is perdiodic only in two spatial dimensions.
 * The non-periodic dimension is selected via the GRAVITY_TALLBOX flag. Also in this case, arbitrarily stretched boxes can
 * be used, and one can use PM_ZOOM_OPTIMIZED and/or FFT_COLUMN_BASED if desired.
 *
 * Much of the code is multi-threaded, so there should be some speed-up if OpenMP is used with NUM_THREADS > 1, but the
 * benefit may be limited because the data transfer steps (which weigh in quite heavily) are not accelerated by this.
 *
 * If eight times the particle load per processor exceeds 2^31 ~ 2 billion, one should activate NUMPART_PER_TASK_LARGE.
 * The code will check this condition and terminate if this is violated, so there should hopefully be no severe risk
 * to accidentally forget this.
 */

#if defined(PMGRID) && defined(PERIODIC)

#define GRIDX (PMGRID * STRETCHX * DBX + DBX_EXTRA)
#define GRIDY (PMGRID * STRETCHY * DBY + DBY_EXTRA)
#define GRIDZ (PMGRID * STRETCHZ * DBZ + DBZ_EXTRA)

#define  GRIDz (GRIDZ/2 + 1)
#define  GRID2 (2*GRIDz)

#if (GRIDX > 1024) || (GRIDY > 1024) || (GRIDZ > 1024)
typedef long long large_array_offset;   /* use a larger data type in this case so that we can always address all cells of the 3D grid with a single index */
#else
typedef unsigned int large_array_offset;
#endif

#ifdef NUMPART_PER_TASK_LARGE
typedef long long large_numpart_type;   /* if there is a risk that the local particle number times 8 overflows a 32-bit integer, this data type should be used */
#else
typedef int large_numpart_type;
#endif

/* short-cut macros for accessing different 3D arrays */
#define FI(x,y,z)  (((large_array_offset) GRID2) * (GRIDY * (x) + (y)) + (z))
#define FC(c,z)    (((large_array_offset) GRID2) * ((c) - myplan.base_firstcol) + (z))
#ifndef FFT_COLUMN_BASED
#define NI(x,y,z)  (((large_array_offset) GRIDZ) * ((y) + (x) * myplan.nslab_y) + (z))
#endif

/* variables for power spectrum estimation */
#ifndef BINS_PS
#define BINS_PS  2000           /* number of bins for power spectrum computation */
#endif
#ifndef POWERSPEC_FOLDFAC
#define POWERSPEC_FOLDFAC 16.    /* folding factor to obtain an estimate of the power spectrum on very small scales */
#endif

static void pmforce_measure_powerspec(int flag, int *typeflag);
static void pmforce_do_powerspec(int *typeflag);
static double pmperiodic_tallbox_long_range_potential(double x, double y, double z);
static void pmforce_setup_tallbox_kernel(void);

static fft_plan myplan;         /*!< In this structure, various bookkeeping variables for the distributed FFTs are stored */

/*! \var maxfftsize
 *  \brief maximum size of the local fft grid among all tasks
 */
static size_t maxfftsize;

/*! \var rhogrid
 *  \brief This array hold the local part of the density field and
 *  after the FFTs the local part of the potential
 *
 *  \var forcegrid
 *  \brief This array will contain the force field
 *
 *  \var workspace
 *  \brief Workspace array used during the FFTs
 */
static fft_real *rhogrid, *forcegrid, *workspace;

/*! \brief Array containing the FFT of #rhogrid
 *
 *  This pointer points to the same array as #rhogrid,
 *  because in-place FFTs are used.
 */
static fft_complex *fft_of_rhogrid;

#if defined(GRAVITY_TALLBOX)
static fft_real *kernel;        /*!< If the tallbox option is used, the code will construct and store the k-space Greens function by FFTing it from real space */
static fft_complex *fft_of_kernel;
#endif

/* Variable for power spectrum calculation */
static double power_spec_totmass, power_spec_totmass2;
static long long power_spec_totnumpart;

/*! \brief This routine generates the FFTW-plans to carry out the FFTs later on.
 *
 *  Some auxiliary variables for bookkeeping are also initialized.
 */
void pm_init_periodic(void)
{
#ifdef LONG_X
  if(LONG_X != (int) (LONG_X))
    terminate("LONG_X must be an integer if used with PMGRID");
#endif

#ifdef LONG_Y
  if(LONG_Y != (int) (LONG_Y))
    terminate("LONG_Y must be an integer if used with PMGRID");
#endif

#ifdef LONG_Z
  if(LONG_Z != (int) (LONG_Z))
    terminate("LONG_Z must be an integer if used with PMGRID");
#endif

  All.Asmth[0] = ASMTH * All.BoxSize / PMGRID;
  All.Rcut[0] = RCUT * All.Asmth[0];

  /* Set up the FFTW-3 plan files. */
  int ndimx[1] = { GRIDX };     /* dimension of the 1D transforms */
  int ndimy[1] = { GRIDY };     /* dimension of the 1D transforms */
  int ndimz[1] = { GRIDZ };     /* dimension of the 1D transforms */

  int max_GRID2 = 2 * (imax(imax(GRIDX, GRIDY), GRIDZ) / 2 + 1);

  /* temporarily allocate some arrays to make sure that out-of-place plans are created */
  rhogrid = (fft_real *) mymalloc("rhogrid", max_GRID2 * sizeof(fft_real));
  forcegrid = (fft_real *) mymalloc("forcegrid", max_GRID2 * sizeof(fft_real));

#ifdef DOUBLEPRECISION_FFTW
  int alignflag = 0;
#else
  /* for single precision, the start of our FFT columns is presently only guaranteed to be 8-byte aligned */
  int alignflag = FFTW_UNALIGNED;
#endif

  myplan.forward_plan_zdir = FFTW(plan_many_dft_r2c) (1, ndimz, 1, rhogrid, 0, 1, GRID2, (fft_complex *) forcegrid, 0, 1, GRIDz, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

#ifndef FFT_COLUMN_BASED
  int stride = GRIDz;
#else
  int stride = 1;
#endif

  myplan.forward_plan_ydir =
    FFTW(plan_many_dft) (1, ndimy, 1, (fft_complex *) rhogrid, 0, stride, GRIDz * GRIDY, (fft_complex *) forcegrid, 0, stride, GRIDz * GRIDY, FFTW_FORWARD,
                         FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  myplan.forward_plan_xdir =
    FFTW(plan_many_dft) (1, ndimx, 1, (fft_complex *) rhogrid, 0, stride, GRIDz * GRIDX, (fft_complex *) forcegrid, 0, stride, GRIDz * GRIDX, FFTW_FORWARD,
                         FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  myplan.backward_plan_xdir =
    FFTW(plan_many_dft) (1, ndimx, 1, (fft_complex *) rhogrid, 0, stride, GRIDz * GRIDX, (fft_complex *) forcegrid, 0, stride, GRIDz * GRIDX, FFTW_BACKWARD,
                         FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  myplan.backward_plan_ydir =
    FFTW(plan_many_dft) (1, ndimy, 1, (fft_complex *) rhogrid, 0, stride, GRIDz * GRIDY, (fft_complex *) forcegrid, 0, stride, GRIDz * GRIDY, FFTW_BACKWARD,
                         FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  myplan.backward_plan_zdir = FFTW(plan_many_dft_c2r) (1, ndimz, 1, (fft_complex *) rhogrid, 0, 1, GRIDz, forcegrid, 0, 1, GRID2, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  myfree(forcegrid);
  myfree(rhogrid);

#ifndef FFT_COLUMN_BASED

  my_slab_based_fft_init(&myplan, GRIDX, GRIDY, GRIDZ);

  maxfftsize = imax(myplan.largest_x_slab * GRIDY, myplan.largest_y_slab * GRIDX) * ((size_t) GRID2);

#else

  my_column_based_fft_init(&myplan, GRIDX, GRIDY, GRIDZ);

  maxfftsize = myplan.max_datasize;

#endif

#if defined(GRAVITY_TALLBOX)
  kernel = (fft_real *) mymalloc("kernel", maxfftsize * sizeof(fft_real));
  fft_of_kernel = (fft_complex *) kernel;

  pmforce_setup_tallbox_kernel();
#endif
}

/* Below, the two functions
 *
 *           pmforce_ ...... _prepare_density()
 * and
 *           pmforce_ ...... _readout_forces_or_potential(int dim)
 *
 * are defined in two different versions, one that works better for uniform
 * simulations, the other for zoom-in runs. Only one of the two sets is used,
 * depending on the setting of PM_ZOOM_OPTIMIZED.
 */



#ifdef PM_ZOOM_OPTIMIZED

static void mysort_pmperiodic(void *b, size_t n, size_t s, int (*cmp) (const void *, const void *));
static int pm_periodic_compare_sortindex(const void *a, const void *b);

/*! \brief This structure links the particles to the mesh cells, to which they contribute their mass
 *
 * Each particle will have eight items of this structure in the #part array.
 * For each of the eight mesh cells the CIC assignment will contribute,
 * one item of this struct exists.
 */
static struct part_slab_data
{
  large_array_offset globalindex;       /*!< index in the global density mesh */
  large_numpart_type partindex; /*!< contains the local particle index shifted by 2^3, the first three bits encode to which part of the CIC assignment this item belongs to */
  large_array_offset localindex;        /*!< index to a local copy of the corresponding mesh cell of the global density array (used during local mass and force assignment) */
} *part;                        /*!< array of part_slab_data linking the local particles to their mesh cells */

static size_t *localfield_sendcount, *localfield_first, *localfield_offset, *localfield_recvcount;
static large_array_offset *localfield_globalindex, *import_globalindex;
static fft_real *localfield_data, *import_data;


void pmforce_zoom_optimized_prepare_density(int mode, int *typelist)
{
  large_numpart_type i;
  int level, recvTask;
  MPI_Status status;

  double to_slab_fac = PMGRID / All.BoxSize;    /* note: This is the same as GRIDX / (All.BoxSize * LONG_X), and similarly for each dimension */

  if(mode == 2)
    to_slab_fac *= POWERSPEC_FOLDFAC;
  if(mode == 3)
    to_slab_fac *= POWERSPEC_FOLDFAC * POWERSPEC_FOLDFAC;

  part = (struct part_slab_data *) mymalloc("part", 8 * (NumPart * sizeof(struct part_slab_data)));
  large_numpart_type *part_sortindex = (large_numpart_type *) mymalloc("part_sortindex", 8 * (NumPart * sizeof(large_numpart_type)));

  /* determine the cells each particle accesses */
#pragma omp parallel for
  for(i = 0; i < NumPart; i++)
    {
      MyDouble *pos;

#ifdef CELL_CENTER_GRAVITY
      MyDouble posw[3], xtmp, ytmp, ztmp;
      if(P[i].Type == 0)
        {
          posw[0] = WRAP_X(SphP[i].Center[0]);
          posw[1] = WRAP_Y(SphP[i].Center[1]);
          posw[2] = WRAP_Z(SphP[i].Center[2]);

          pos = posw;
        }
      else
#endif
        pos = P[i].Pos;

      int slab_x = (int) (to_slab_fac * pos[0]);
      int slab_y = (int) (to_slab_fac * pos[1]);
      int slab_z = (int) (to_slab_fac * pos[2]);

      if(mode >= 2)
        {
          slab_x %= GRIDX;
          slab_y %= GRIDY;
          slab_z %= GRIDZ;
        }
      else
        {
          if(slab_x >= GRIDX)
            slab_x -= GRIDX;
          if(slab_y >= GRIDY)
            slab_y -= GRIDY;
          if(slab_z >= GRIDZ)
            slab_z -= GRIDZ;
        }

      large_numpart_type index_on_grid = ((large_numpart_type) i) << 3;

      for(int xx = 0; xx < 2; xx++)
        for(int yy = 0; yy < 2; yy++)
          for(int zz = 0; zz < 2; zz++)
            {
              int slab_xx = slab_x + xx;
              int slab_yy = slab_y + yy;
              int slab_zz = slab_z + zz;

              if(slab_xx >= GRIDX)
                slab_xx -= GRIDX;
              if(slab_yy >= GRIDY)
                slab_yy -= GRIDY;
              if(slab_zz >= GRIDZ)
                slab_zz -= GRIDZ;

              large_array_offset offset = FI(slab_xx, slab_yy, slab_zz);

              part[index_on_grid].partindex = (i << 3) + (xx << 2) + (yy << 1) + zz;
              part[index_on_grid].globalindex = offset;
              part_sortindex[index_on_grid] = index_on_grid;
              index_on_grid++;
            }
    }

  /* note: num_on_grid will be  8 times larger than the particle number, but num_field_points will generally be much smaller */

  large_array_offset num_field_points;
  large_numpart_type num_on_grid = ((large_numpart_type) NumPart) << 3;

  /* bring the part-field into the order of the accessed cells. This allows the removal of duplicates */
  mysort_pmperiodic(part_sortindex, num_on_grid, sizeof(large_numpart_type), pm_periodic_compare_sortindex);

  if(num_on_grid > 0)
    num_field_points = 1;
  else
    num_field_points = 0;

  /* determine the number of unique field points */
#pragma omp parallel for reduction(+ : num_field_points)
  for(i = 1; i < num_on_grid; i++)
    {
      if(part[part_sortindex[i]].globalindex != part[part_sortindex[i - 1]].globalindex)
        num_field_points++;
    }

  /* allocate the local field */
  localfield_globalindex = (large_array_offset *) mymalloc_movable(&localfield_globalindex, "localfield_globalindex", num_field_points * sizeof(large_array_offset));
  localfield_data = (fft_real *) mymalloc_movable(&localfield_data, "localfield_data", num_field_points * sizeof(fft_real));
  localfield_first = (size_t *) mymalloc_movable(&localfield_first, "localfield_first", NTask * sizeof(size_t));
  localfield_sendcount = (size_t *) mymalloc_movable(&localfield_sendcount, "localfield_sendcount", NTask * sizeof(size_t));
  localfield_offset = (size_t *) mymalloc_movable(&localfield_offset, "localfield_offset", NTask * sizeof(size_t));
  localfield_recvcount = (size_t *) mymalloc_movable(&localfield_recvcount, "localfield_recvcount", NTask * sizeof(size_t));

  for(i = 0; i < NTask; i++)
    {
      localfield_first[i] = 0;
      localfield_sendcount[i] = 0;
    }

  /* establish the cross link between the part[ ]-array and the local list of
   * mesh points. Also, count on which CPU the needed field points are stored.
   */
  for(i = 0, num_field_points = 0; i < num_on_grid; i++)
    {
      if(i > 0)
        if(part[part_sortindex[i]].globalindex != part[part_sortindex[i - 1]].globalindex)
          num_field_points++;

      part[part_sortindex[i]].localindex = num_field_points;

      if(i > 0)
        if(part[part_sortindex[i]].globalindex == part[part_sortindex[i - 1]].globalindex)
          continue;

      localfield_globalindex[num_field_points] = part[part_sortindex[i]].globalindex;

#ifndef FFT_COLUMN_BASED
      int slab = part[part_sortindex[i]].globalindex / (GRIDY * GRID2);
      int task = myplan.slab_to_task[slab];
#else
      int task, column = part[part_sortindex[i]].globalindex / (GRID2);

      if(column < myplan.pivotcol)
        task = column / myplan.avg;
      else
        task = (column - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;
#endif

      if(localfield_sendcount[task] == 0)
        localfield_first[task] = num_field_points;

      localfield_sendcount[task]++;
    }
  num_field_points++;

  for(i = 1, localfield_offset[0] = 0; i < NTask; i++)
    localfield_offset[i] = localfield_offset[i - 1] + localfield_sendcount[i - 1];

  myfree_movable(part_sortindex);
  part_sortindex = NULL;

  /* now bin the local particle data onto the mesh list */
#pragma omp parallel for
  for(i = 0; i < num_field_points; i++)
    localfield_data[i] = 0;

  for(i = 0; i < num_on_grid; i += 8)
    {
      int pindex = (part[i].partindex >> 3);

      MyDouble *pos;
#ifdef CELL_CENTER_GRAVITY
      MyDouble posw[3], xtmp, ytmp, ztmp;
      if(P[pindex].Type == 0)
        {
          posw[0] = WRAP_X(SphP[pindex].Center[0]);
          posw[1] = WRAP_Y(SphP[pindex].Center[1]);
          posw[2] = WRAP_Z(SphP[pindex].Center[2]);

          pos = posw;
        }
      else
#endif
        pos = P[pindex].Pos;

      int slab_x = (int) (to_slab_fac * pos[0]);
      int slab_y = (int) (to_slab_fac * pos[1]);
      int slab_z = (int) (to_slab_fac * pos[2]);

      double dx = to_slab_fac * pos[0] - slab_x;
      double dy = to_slab_fac * pos[1] - slab_y;
      double dz = to_slab_fac * pos[2] - slab_z;

      double weight = P[pindex].Mass;

      if(mode)                  /* only for power spectrum calculation */
        if(typelist[P[pindex].Type] == 0)
          continue;

      localfield_data[part[i + 0].localindex] += weight * (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
      localfield_data[part[i + 1].localindex] += weight * (1.0 - dx) * (1.0 - dy) * dz;
      localfield_data[part[i + 2].localindex] += weight * (1.0 - dx) * dy * (1.0 - dz);
      localfield_data[part[i + 3].localindex] += weight * (1.0 - dx) * dy * dz;
      localfield_data[part[i + 4].localindex] += weight * (dx) * (1.0 - dy) * (1.0 - dz);
      localfield_data[part[i + 5].localindex] += weight * (dx) * (1.0 - dy) * dz;
      localfield_data[part[i + 6].localindex] += weight * (dx) * dy * (1.0 - dz);
      localfield_data[part[i + 7].localindex] += weight * (dx) * dy * dz;
    }

  rhogrid = (fft_real *) mymalloc("rhogrid", maxfftsize * sizeof(fft_real));

  /* clear local FFT-mesh density field */
  large_array_offset ii;
#pragma omp parallel for
  for(ii = 0; ii < maxfftsize; ii++)
    rhogrid[ii] = 0;


  /* exchange data and add contributions to the local mesh-path */
  MPI_Alltoall(localfield_sendcount, sizeof(size_t), MPI_BYTE, localfield_recvcount, sizeof(size_t), MPI_BYTE, MPI_COMM_WORLD);

  for(level = 0; level < (1 << PTask); level++) /* note: for level=0, target is the same task */
    {
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
        {
          if(level > 0)
            {
              import_data = (fft_real *) mymalloc("import_data", localfield_recvcount[recvTask] * sizeof(fft_real));
              import_globalindex = (large_array_offset *) mymalloc("import_globalindex", localfield_recvcount[recvTask] * sizeof(large_array_offset));

              if(localfield_sendcount[recvTask] > 0 || localfield_recvcount[recvTask] > 0)
                {
                  myMPI_Sendrecv(localfield_data + localfield_offset[recvTask],
                                 localfield_sendcount[recvTask] * sizeof(fft_real), MPI_BYTE, recvTask, TAG_NONPERIOD_A, import_data, localfield_recvcount[recvTask] * sizeof(fft_real), MPI_BYTE,
                                 recvTask, TAG_NONPERIOD_A, MPI_COMM_WORLD, &status);

                  myMPI_Sendrecv(localfield_globalindex + localfield_offset[recvTask],
                                 localfield_sendcount[recvTask] * sizeof(large_array_offset), MPI_BYTE, recvTask, TAG_NONPERIOD_B,
                                 import_globalindex, localfield_recvcount[recvTask] * sizeof(large_array_offset), MPI_BYTE, recvTask, TAG_NONPERIOD_B, MPI_COMM_WORLD, &status);
                }
            }
          else
            {
              import_data = localfield_data + localfield_offset[ThisTask];
              import_globalindex = localfield_globalindex + localfield_offset[ThisTask];
            }

          /* note: here every element in rhogrid is only accessed once, so there should be no race condition */
#pragma omp parallel for
          for(i = 0; i < localfield_recvcount[recvTask]; i++)
            {
              /* determine offset in local FFT slab */
#ifndef FFT_COLUMN_BASED
              large_array_offset offset = import_globalindex[i] - myplan.first_slab_x_of_task[ThisTask] * GRIDY * ((large_array_offset) GRID2);
#else
              large_array_offset offset = import_globalindex[i] - myplan.base_firstcol * ((large_array_offset) GRID2);
#endif
              rhogrid[offset] += import_data[i];
            }

          if(level > 0)
            {
              myfree(import_globalindex);
              myfree(import_data);
            }
        }
    }
#ifdef MODGRAV_EFF_MASS
  modgrav_add_amr_nodes_to_periodic_grid(rhogrid, to_slab_fac, maxfftsize, myplan);
#endif
}


/* Function to read out the force component corresponding to spatial dimension 'dim'.
 * If dim is negative, potential values are read out and assigned to particles.
 */
void pmforce_zoom_optimized_readout_forces_or_potential(int dim)
{
#ifdef EVALPOTENTIAL
#ifdef GRAVITY_TALLBOX
  double fac = All.G / (((double) GRIDX) * GRIDY * GRIDZ);      /* to get potential  */
#else
  double fac = 4 * M_PI * All.G / (pow(All.BoxSize, 3) * STRETCHX * STRETCHY * STRETCHZ);       /* to get potential  */
#endif
#endif

  large_numpart_type i;
  int level, recvTask;
  MPI_Status status;

  fft_real *grid;

  if(dim < 0)
    grid = rhogrid;
  else
    grid = forcegrid;

  double to_slab_fac = PMGRID / All.BoxSize;

  for(level = 0; level < (1 << PTask); level++) /* note: for level=0, target is the same task */
    {
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
        {
          if(level > 0)
            {
              import_data = (fft_real *) mymalloc("import_data", localfield_recvcount[recvTask] * sizeof(fft_real));
              import_globalindex = (large_array_offset *) mymalloc("import_globalindex", localfield_recvcount[recvTask] * sizeof(large_array_offset));

              if(localfield_sendcount[recvTask] > 0 || localfield_recvcount[recvTask] > 0)
                {
                  myMPI_Sendrecv(localfield_globalindex + localfield_offset[recvTask],
                                 localfield_sendcount[recvTask] * sizeof(large_array_offset), MPI_BYTE, recvTask, TAG_NONPERIOD_C,
                                 import_globalindex, localfield_recvcount[recvTask] * sizeof(large_array_offset), MPI_BYTE, recvTask, TAG_NONPERIOD_C, MPI_COMM_WORLD, &status);
                }
            }
          else
            {
              import_data = localfield_data + localfield_offset[ThisTask];
              import_globalindex = localfield_globalindex + localfield_offset[ThisTask];
            }

          for(i = 0; i < localfield_recvcount[recvTask]; i++)
            {
#ifndef FFT_COLUMN_BASED
              large_array_offset offset = import_globalindex[i] - myplan.first_slab_x_of_task[ThisTask] * GRIDY * ((large_array_offset) GRID2);
#else
              large_array_offset offset = import_globalindex[i] - myplan.base_firstcol * ((large_array_offset) GRID2);
#endif
              import_data[i] = grid[offset];
            }

          if(level > 0)
            {
              myMPI_Sendrecv(import_data, localfield_recvcount[recvTask] * sizeof(fft_real), MPI_BYTE, recvTask,
                             TAG_NONPERIOD_A, localfield_data + localfield_offset[recvTask], localfield_sendcount[recvTask] * sizeof(fft_real), MPI_BYTE, recvTask, TAG_NONPERIOD_A, MPI_COMM_WORLD,
                             &status);

              myfree(import_globalindex);
              myfree(import_data);
            }
        }
    }

  /* read out the froce/potential values, which all have been assembled in localfield_data */
#pragma omp parallel for
  for(i = 0; i < NumPart; i++)
    {
      large_numpart_type j = (i << 3);

      MyDouble *pos;

#ifdef CELL_CENTER_GRAVITY
      MyDouble posw[3], xtmp, ytmp, ztmp;
      if(P[i].Type == 0)
        {
          posw[0] = WRAP_X(SphP[i].Center[0]);
          posw[1] = WRAP_Y(SphP[i].Center[1]);
          posw[2] = WRAP_Z(SphP[i].Center[2]);

          pos = posw;
        }
      else
#endif
        pos = P[i].Pos;

      int slab_x = (int) (to_slab_fac * pos[0]);
      double dx = to_slab_fac * pos[0] - slab_x;

      int slab_y = (int) (to_slab_fac * pos[1]);
      double dy = to_slab_fac * pos[1] - slab_y;

      int slab_z = (int) (to_slab_fac * pos[2]);
      double dz = to_slab_fac * pos[2] - slab_z;

      double value =
        +localfield_data[part[j + 0].localindex] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz)
        + localfield_data[part[j + 1].localindex] * (1.0 - dx) * (1.0 - dy) * dz
        + localfield_data[part[j + 2].localindex] * (1.0 - dx) * dy * (1.0 - dz)
        + localfield_data[part[j + 3].localindex] * (1.0 - dx) * dy * dz
        + localfield_data[part[j + 4].localindex] * (dx) * (1.0 - dy) * (1.0 - dz)
        + localfield_data[part[j + 5].localindex] * (dx) * (1.0 - dy) * dz + localfield_data[part[j + 6].localindex] * (dx) * dy * (1.0 - dz) +
        localfield_data[part[j + 7].localindex] * (dx) * dy * dz;

      if(dim < 0)
        {
#ifdef EVALPOTENTIAL
          P[i].PM_Potential += value * fac;
#endif
        }
      else
        P[i].GravPM[dim] += value;
    }
}


#else


/*
 *  Here come the routines for a different communication algorithm that is better suited for a homogenuously loaded boxes.
 */
static struct partbuf
{
  MyFloat Mass;
  MyFloat Pos[3];
} *partin, *partout;

static size_t nimport, nexport;

static size_t *Sndpm_count, *Sndpm_offset;
static size_t *Rcvpm_count, *Rcvpm_offset;

static void pmforce_uniform_optimized_prepare_density(int mode)
{
  int i, j;

  double to_slab_fac = PMGRID / All.BoxSize;

  if(mode == 2)
    to_slab_fac *= POWERSPEC_FOLDFAC;
  if(mode == 3)
    to_slab_fac *= POWERSPEC_FOLDFAC * POWERSPEC_FOLDFAC;

  /* We here enlarge NTask such that each thread gets his own cache line for send_count/send_offset.
   * This should hopefully prevent a performance penalty from 'false sharing' for these variables
   */
  int multiNtask = roundup_to_multiple_of_cacheline_size(NTask * sizeof(size_t)) / sizeof(size_t);

  Sndpm_count = (size_t *) mymalloc("Sndpm_count", MaxThreads * multiNtask * sizeof(size_t));
  Sndpm_offset = (size_t *) mymalloc("Sndpm_offset", MaxThreads * multiNtask * sizeof(size_t));
  Rcvpm_count = (size_t *) mymalloc("Rcvpm_count", NTask * sizeof(size_t));
  Rcvpm_offset = (size_t *) mymalloc("Rcvpm_offset", NTask * sizeof(size_t));

  /* determine the slabs/columns each particles accesses */
#pragma omp parallel private(j)
  {
    size_t *send_count = Sndpm_count + get_thread_num() * multiNtask;

    /* each threads needs to do theloop to clear its send_count[] array */
    for(j = 0; j < NTask; j++)
      send_count[j] = 0;

#pragma omp for schedule(static) private(i)
    for(i = 0; i < NumPart; i++)
      {
        MyDouble *pos;

#ifdef CELL_CENTER_GRAVITY
        MyDouble posw[3], xtmp, ytmp, ztmp;
        if(P[i].Type == 0)
          {
            posw[0] = WRAP_X(SphP[i].Center[0]);
            posw[1] = WRAP_Y(SphP[i].Center[1]);
            posw[2] = WRAP_Z(SphP[i].Center[2]);

            pos = posw;
          }
        else
#endif
          pos = P[i].Pos;

        int slab_x = (int) (to_slab_fac * pos[0]);
        int slab_xx = slab_x + 1;

        if(mode >= 2)
          {
            slab_x %= GRIDX;
            slab_xx %= GRIDX;
          }
        else
          {
            if(slab_x >= GRIDX)
              slab_x -= GRIDX;

            if(slab_xx >= GRIDX)
              slab_xx -= GRIDX;
          }

#ifndef FFT_COLUMN_BASED
        int task0 = myplan.slab_to_task[slab_x];
        int task1 = myplan.slab_to_task[slab_xx];

        send_count[task0]++;
        if(task0 != task1)
          send_count[task1]++;
#else
        int slab_y = (int) (to_slab_fac * pos[1]);
        int slab_yy = slab_y + 1;

        if(mode >= 2)
          {
            slab_y %= GRIDY;
            slab_yy %= GRIDY;
          }
        else
          {
            if(slab_y >= GRIDY)
              slab_y -= GRIDY;

            if(slab_yy >= GRIDY)
              slab_yy -= GRIDY;
          }

        int column0 = slab_x * GRIDY + slab_y;
        int column1 = slab_x * GRIDY + slab_yy;
        int column2 = slab_xx * GRIDY + slab_y;
        int column3 = slab_xx * GRIDY + slab_yy;

        int task0, task1, task2, task3;

        if(column0 < myplan.pivotcol)
          task0 = column0 / myplan.avg;
        else
          task0 = (column0 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

        if(column1 < myplan.pivotcol)
          task1 = column1 / myplan.avg;
        else
          task1 = (column1 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

        if(column2 < myplan.pivotcol)
          task2 = column2 / myplan.avg;
        else
          task2 = (column2 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

        if(column3 < myplan.pivotcol)
          task3 = column3 / myplan.avg;
        else
          task3 = (column3 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

        send_count[task0]++;
        if(task1 != task0)
          send_count[task1]++;
        if(task2 != task1 && task2 != task0)
          send_count[task2]++;
        if(task3 != task0 && task3 != task1 && task3 != task2)
          send_count[task3]++;
#endif
      }
  }


  /* collect thread-specific offset table and collect the results from the other threads */
  for(i = 0, Sndpm_offset[0] = 0; i < NTask; i++)
    for(j = 0; j < MaxThreads; j++)
      {
        int ind_prev, ind = j * multiNtask + i;
        if(ind > 0)
          {
            if(j == 0)
              ind_prev = (MaxThreads - 1) * multiNtask + i - 1;
            else
              ind_prev = ind - multiNtask;

            Sndpm_offset[ind] = Sndpm_offset[ind_prev] + Sndpm_count[ind_prev];
          }
      }

  for(j = 1; j < MaxThreads; j++)
    for(i = 0; i < NTask; i++)
      Sndpm_count[i] += Sndpm_count[i + j * multiNtask];

  MPI_Alltoall(Sndpm_count, sizeof(size_t), MPI_BYTE, Rcvpm_count, sizeof(size_t), MPI_BYTE, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, nexport = 0, Rcvpm_offset[0] = 0, Sndpm_offset[0] = 0; j < NTask; j++)
    {
      nexport += Sndpm_count[j];
      nimport += Rcvpm_count[j];

      if(j > 0)
        {
          Sndpm_offset[j] = Sndpm_offset[j - 1] + Sndpm_count[j - 1];
          Rcvpm_offset[j] = Rcvpm_offset[j - 1] + Rcvpm_count[j - 1];
        }
    }

  /* allocate import and export buffer */
  partin = (struct partbuf *) mymalloc("partin", nimport * sizeof(struct partbuf));
  partout = (struct partbuf *) mymalloc("partout", nexport * sizeof(struct partbuf));


#pragma omp parallel private(j)
  {
    size_t *send_count = Sndpm_count + get_thread_num() * multiNtask;
    size_t *send_offset = Sndpm_offset + get_thread_num() * multiNtask;

    for(j = 0; j < NTask; j++)
      send_count[j] = 0;

    /* fill export buffer */
#pragma omp for schedule(static)
    for(i = 0; i < NumPart; i++)
      {
        MyDouble *pos;

#ifdef CELL_CENTER_GRAVITY
        MyDouble posw[3], xtmp, ytmp, ztmp;
        if(P[i].Type == 0)
          {
            posw[0] = WRAP_X(SphP[i].Center[0]);
            posw[1] = WRAP_Y(SphP[i].Center[1]);
            posw[2] = WRAP_Z(SphP[i].Center[2]);

            pos = posw;
          }
        else
#endif
          pos = P[i].Pos;

        int slab_x = (int) (to_slab_fac * pos[0]);
        int slab_xx = slab_x + 1;

        if(mode >= 2)
          {
            slab_x %= GRIDX;
            slab_xx %= GRIDX;
          }
        else
          {
            if(slab_x >= GRIDX)
              slab_x -= GRIDX;

            if(slab_xx >= GRIDX)
              slab_xx -= GRIDX;
          }

#ifndef FFT_COLUMN_BASED
        int task0 = myplan.slab_to_task[slab_x];
        int task1 = myplan.slab_to_task[slab_xx];

        size_t ind0 = send_offset[task0] + send_count[task0]++;
        partout[ind0].Mass = P[i].Mass;
        for(j = 0; j < 3; j++)
          partout[ind0].Pos[j] = pos[j];

        if(task0 != task1)
          {
            size_t ind1 = send_offset[task1] + send_count[task1]++;
            partout[ind1].Mass = P[i].Mass;
            for(j = 0; j < 3; j++)
              partout[ind1].Pos[j] = pos[j];
          }
#else
        int slab_y = (int) (to_slab_fac * pos[1]);
        int slab_yy = slab_y + 1;

        if(mode >= 2)
          {
            slab_y %= GRIDY;
            slab_yy %= GRIDY;
          }
        else
          {
            if(slab_y >= GRIDY)
              slab_y -= GRIDY;

            if(slab_yy >= GRIDY)
              slab_yy -= GRIDY;
          }

        int column0 = slab_x * GRIDY + slab_y;
        int column1 = slab_x * GRIDY + slab_yy;
        int column2 = slab_xx * GRIDY + slab_y;
        int column3 = slab_xx * GRIDY + slab_yy;

        int task0, task1, task2, task3;

        if(column0 < myplan.pivotcol)
          task0 = column0 / myplan.avg;
        else
          task0 = (column0 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

        if(column1 < myplan.pivotcol)
          task1 = column1 / myplan.avg;
        else
          task1 = (column1 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

        if(column2 < myplan.pivotcol)
          task2 = column2 / myplan.avg;
        else
          task2 = (column2 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

        if(column3 < myplan.pivotcol)
          task3 = column3 / myplan.avg;
        else
          task3 = (column3 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

        size_t ind0 = send_offset[task0] + send_count[task0]++;
        partout[ind0].Mass = P[i].Mass;
        for(j = 0; j < 3; j++)
          partout[ind0].Pos[j] = pos[j];

        if(task1 != task0)
          {
            size_t ind1 = send_offset[task1] + send_count[task1]++;
            partout[ind1].Mass = P[i].Mass;
            for(j = 0; j < 3; j++)
              partout[ind1].Pos[j] = pos[j];
          }
        if(task2 != task1 && task2 != task0)
          {
            size_t ind2 = send_offset[task2] + send_count[task2]++;
            partout[ind2].Mass = P[i].Mass;
            for(j = 0; j < 3; j++)
              partout[ind2].Pos[j] = pos[j];
          }
        if(task3 != task0 && task3 != task1 && task3 != task2)
          {
            size_t ind3 = send_offset[task3] + send_count[task3]++;
            partout[ind3].Mass = P[i].Mass;
            for(j = 0; j < 3; j++)
              partout[ind3].Pos[j] = pos[j];
          }
#endif
      }
  }

  /* collect the send_count[] results from the other threads */
  for(j = 1; j < MaxThreads; j++)
    for(i = 0; i < NTask; i++)
      Sndpm_count[i] += Sndpm_count[i + j * multiNtask];

  int flag_big = 0, flag_big_all;
  for(i = 0; i < NTask; i++)
    if(Sndpm_count[i] * sizeof(struct partbuf) > MPI_MESSAGE_SIZELIMIT_IN_BYTES)
      flag_big = 1;

  /* produce a flag if any of the send sizes is above our transfer limit, in this case we will
   * transfer the data in chunks.
   */
  MPI_Allreduce(&flag_big, &flag_big_all, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  /* exchange particle data */
  myMPI_Alltoallv(partout, Sndpm_count, Sndpm_offset, partin, Rcvpm_count, Rcvpm_offset, sizeof(struct partbuf), flag_big_all, MPI_COMM_WORLD);

  myfree(partout);

  /* allocate density field */
  rhogrid = (fft_real *) mymalloc("rhogrid", maxfftsize * sizeof(fft_real));

  /* clear local FFT-mesh density field */
  large_array_offset ii;
#pragma omp parallel for
  for(ii = 0; ii < maxfftsize; ii++)
    rhogrid[ii] = 0;

#ifndef FFT_COLUMN_BASED
  /* bin particle data onto mesh, in multi-threaded fashion */
#pragma omp parallel private(i)
  {
    int tid = get_thread_num();

    int first_y, count_y;
    subdivide_evenly(GRIDY, MaxThreads, tid, &first_y, &count_y);
    int last_y = first_y + count_y - 1;

    for(i = 0; i < nimport; i++)
      {
        int slab_y = (int) (to_slab_fac * partin[i].Pos[1]);
        int slab_yy = slab_y + 1;
        double dy = to_slab_fac * partin[i].Pos[1] - slab_y;

        if(mode >= 2)
          {
            slab_y %= GRIDY;
            slab_yy %= GRIDY;
          }
        else
          {
            if(slab_y >= GRIDY)
              slab_y -= GRIDY;

            if(slab_yy >= GRIDY)
              slab_yy -= GRIDY;
          }

        int flag_slab_y, flag_slab_yy;

        if(slab_y >= first_y && slab_y <= last_y)
          flag_slab_y = 1;
        else
          flag_slab_y = 0;

        if(slab_yy >= first_y && slab_yy <= last_y)
          flag_slab_yy = 1;
        else
          flag_slab_yy = 0;

        if(flag_slab_y || flag_slab_yy)
          {
            double mass = partin[i].Mass;

            int slab_x = (int) (to_slab_fac * partin[i].Pos[0]);
            int slab_z = (int) (to_slab_fac * partin[i].Pos[2]);
            int slab_xx = slab_x + 1;
            int slab_zz = slab_z + 1;

            double dx = to_slab_fac * partin[i].Pos[0] - slab_x;
            double dz = to_slab_fac * partin[i].Pos[2] - slab_z;

            if(mode >= 2)
              {
                slab_x %= GRIDX;
                slab_z %= GRIDZ;
                slab_xx %= GRIDX;
                slab_zz %= GRIDZ;
              }
            else
              {
                if(slab_x >= GRIDX)
                  slab_x -= GRIDX;
                if(slab_z >= GRIDZ)
                  slab_z -= GRIDZ;

                if(slab_xx >= GRIDX)
                  slab_xx -= GRIDX;
                if(slab_zz >= GRIDZ)
                  slab_zz -= GRIDZ;
              }

            int flag_slab_x, flag_slab_xx;

            if(myplan.slab_to_task[slab_x] == ThisTask)
              {
                slab_x -= myplan.first_slab_x_of_task[ThisTask];
                flag_slab_x = 1;
              }
            else
              flag_slab_x = 0;

            if(myplan.slab_to_task[slab_xx] == ThisTask)
              {
                slab_xx -= myplan.first_slab_x_of_task[ThisTask];
                flag_slab_xx = 1;
              }
            else
              flag_slab_xx = 0;

            if(flag_slab_x)
              {
                if(flag_slab_y)
                  {
                    rhogrid[FI(slab_x, slab_y, slab_z)] += (mass * (1.0 - dx) * (1.0 - dy) * (1.0 - dz));
                    rhogrid[FI(slab_x, slab_y, slab_zz)] += (mass * (1.0 - dx) * (1.0 - dy) * (dz));
                  }

                if(flag_slab_yy)
                  {
                    rhogrid[FI(slab_x, slab_yy, slab_z)] += (mass * (1.0 - dx) * (dy) * (1.0 - dz));
                    rhogrid[FI(slab_x, slab_yy, slab_zz)] += (mass * (1.0 - dx) * (dy) * (dz));
                  }
              }

            if(flag_slab_xx)
              {
                if(flag_slab_y)
                  {
                    rhogrid[FI(slab_xx, slab_y, slab_z)] += (mass * (dx) * (1.0 - dy) * (1.0 - dz));
                    rhogrid[FI(slab_xx, slab_y, slab_zz)] += (mass * (dx) * (1.0 - dy) * (dz));
                  }

                if(flag_slab_yy)
                  {
                    rhogrid[FI(slab_xx, slab_yy, slab_z)] += (mass * (dx) * (dy) * (1.0 - dz));
                    rhogrid[FI(slab_xx, slab_yy, slab_zz)] += (mass * (dx) * (dy) * (dz));
                  }
              }
          }
      }
  }

#else

  struct data_cols
  {
    int col0, col1, col2, col3;
    double dx, dy;
  } *aux;

  aux = mymalloc("aux", nimport * sizeof(struct data_cols));

#pragma omp parallel for
  for(i = 0; i < nimport; i++)
    {
      int slab_x = (int) (to_slab_fac * partin[i].Pos[0]);
      int slab_xx = slab_x + 1;

      int slab_y = (int) (to_slab_fac * partin[i].Pos[1]);
      int slab_yy = slab_y + 1;

      aux[i].dx = to_slab_fac * partin[i].Pos[0] - slab_x;
      aux[i].dy = to_slab_fac * partin[i].Pos[1] - slab_y;

      if(mode >= 2)
        {
          slab_x %= GRIDX;
          slab_xx %= GRIDX;
          slab_y %= GRIDY;
          slab_yy %= GRIDY;
        }
      else
        {
          if(slab_x >= GRIDX)
            slab_x -= GRIDX;
          if(slab_xx >= GRIDX)
            slab_xx -= GRIDX;

          if(slab_y >= GRIDY)
            slab_y -= GRIDY;
          if(slab_yy >= GRIDY)
            slab_yy -= GRIDY;
        }

      aux[i].col0 = slab_x * GRIDY + slab_y;
      aux[i].col1 = slab_x * GRIDY + slab_yy;
      aux[i].col2 = slab_xx * GRIDY + slab_y;
      aux[i].col3 = slab_xx * GRIDY + slab_yy;
    }

#pragma omp parallel private(i)
  {
    int tid = get_thread_num();

    int first_col, last_col, count_col;
    subdivide_evenly(myplan.base_ncol, MaxThreads, tid, &first_col, &count_col);
    last_col = first_col + count_col - 1;
    first_col += myplan.base_firstcol;
    last_col += myplan.base_firstcol;

    for(i = 0; i < nimport; i++)
      {
        int flag0, flag1, flag2, flag3;
        int col0 = aux[i].col0;
        int col1 = aux[i].col1;
        int col2 = aux[i].col2;
        int col3 = aux[i].col3;

        if(col0 >= first_col && col0 <= last_col)
          flag0 = 1;
        else
          flag0 = 0;

        if(col1 >= first_col && col1 <= last_col)
          flag1 = 1;
        else
          flag1 = 0;

        if(col2 >= first_col && col2 <= last_col)
          flag2 = 1;
        else
          flag2 = 0;

        if(col3 >= first_col && col3 <= last_col)
          flag3 = 1;
        else
          flag3 = 0;

        if(flag0 || flag1 || flag2 || flag3)
          {
            double mass = partin[i].Mass;

            double dx = aux[i].dx;
            double dy = aux[i].dy;

            int slab_z = (int) (to_slab_fac * partin[i].Pos[2]);
            int slab_zz = slab_z + 1;

            double dz = to_slab_fac * partin[i].Pos[2] - slab_z;

            if(mode >= 2)
              {
                slab_z %= GRIDZ;
                slab_zz %= GRIDZ;
              }
            else
              {
                if(slab_z >= GRIDZ)
                  slab_z -= GRIDZ;

                if(slab_zz >= GRIDZ)
                  slab_zz -= GRIDZ;
              }

            if(flag0)
              {
                rhogrid[FC(col0, slab_z)] += (mass * (1.0 - dx) * (1.0 - dy) * (1.0 - dz));
                rhogrid[FC(col0, slab_zz)] += (mass * (1.0 - dx) * (1.0 - dy) * (dz));
              }

            if(flag1)
              {
                rhogrid[FC(col1, slab_z)] += (mass * (1.0 - dx) * (dy) * (1.0 - dz));
                rhogrid[FC(col1, slab_zz)] += (mass * (1.0 - dx) * (dy) * (dz));
              }

            if(flag2)
              {
                rhogrid[FC(col2, slab_z)] += (mass * (dx) * (1.0 - dy) * (1.0 - dz));
                rhogrid[FC(col2, slab_zz)] += (mass * (dx) * (1.0 - dy) * (dz));
              }

            if(flag3)
              {
                rhogrid[FC(col3, slab_z)] += (mass * (dx) * (dy) * (1.0 - dz));
                rhogrid[FC(col3, slab_zz)] += (mass * (dx) * (dy) * (dz));
              }
          }
      }
  }

  myfree(aux);

#endif
}


/* If dim<0, this function reads out the potential, otherwise Cartesian force components.
 */
static void pmforce_uniform_optimized_readout_forces_or_potential(int dim)
{
#ifdef EVALPOTENTIAL
#ifdef GRAVITY_TALLBOX
  double fac = All.G / (((double) GRIDX) * GRIDY * GRIDZ);      /* to get potential  */
#else
  double fac = 4 * M_PI * All.G / (pow(All.BoxSize, 3) * STRETCHX * STRETCHY * STRETCHZ);       /* to get potential  */
#endif
#endif

  double to_slab_fac = PMGRID / All.BoxSize;

  double *flistin = (double *) mymalloc("flistin", nimport * sizeof(double));
  double *flistout = (double *) mymalloc("flistout", nexport * sizeof(double));

  fft_real *grid;

  if(dim < 0)
    grid = rhogrid;
  else
    grid = forcegrid;

  size_t i;
#pragma omp parallel for
  for(i = 0; i < nimport; i++)
    {
      flistin[i] = 0;

      int slab_x = (int) (to_slab_fac * partin[i].Pos[0]);
      int slab_y = (int) (to_slab_fac * partin[i].Pos[1]);
      int slab_z = (int) (to_slab_fac * partin[i].Pos[2]);

      double dx = to_slab_fac * partin[i].Pos[0] - slab_x;
      double dy = to_slab_fac * partin[i].Pos[1] - slab_y;
      double dz = to_slab_fac * partin[i].Pos[2] - slab_z;

      if(slab_x >= GRIDX)
        slab_x -= GRIDX;
      if(slab_y >= GRIDY)
        slab_y -= GRIDY;
      if(slab_z >= GRIDZ)
        slab_z -= GRIDZ;

      int slab_xx = slab_x + 1;
      int slab_yy = slab_y + 1;
      int slab_zz = slab_z + 1;

      if(slab_xx >= GRIDX)
        slab_xx -= GRIDX;
      if(slab_yy >= GRIDY)
        slab_yy -= GRIDY;
      if(slab_zz >= GRIDZ)
        slab_zz -= GRIDZ;


#ifndef FFT_COLUMN_BASED
      if(myplan.slab_to_task[slab_x] == ThisTask)
        {
          slab_x -= myplan.first_slab_x_of_task[ThisTask];

          flistin[i] += grid[FI(slab_x, slab_y, slab_z)] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz)
            + grid[FI(slab_x, slab_y, slab_zz)] * (1.0 - dx) * (1.0 - dy) * (dz) + grid[FI(slab_x, slab_yy, slab_z)] * (1.0 - dx) * (dy) * (1.0 - dz) + grid[FI(slab_x, slab_yy, slab_zz)] * (1.0 -
                                                                                                                                                                                              dx) *
            (dy) * (dz);
        }

      if(myplan.slab_to_task[slab_xx] == ThisTask)
        {
          slab_xx -= myplan.first_slab_x_of_task[ThisTask];

          flistin[i] += grid[FI(slab_xx, slab_y, slab_z)] * (dx) * (1.0 - dy) * (1.0 - dz)
            + grid[FI(slab_xx, slab_y, slab_zz)] * (dx) * (1.0 - dy) * (dz) + grid[FI(slab_xx, slab_yy, slab_z)] * (dx) * (dy) * (1.0 - dz) + grid[FI(slab_xx, slab_yy, slab_zz)] * (dx) * (dy) * (dz);
        }
#else
      int column0 = slab_x * GRIDY + slab_y;
      int column1 = slab_x * GRIDY + slab_yy;
      int column2 = slab_xx * GRIDY + slab_y;
      int column3 = slab_xx * GRIDY + slab_yy;

      if(column0 >= myplan.base_firstcol && column0 <= myplan.base_lastcol)
        {
          flistin[i] += grid[FC(column0, slab_z)] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz) + grid[FC(column0, slab_zz)] * (1.0 - dx) * (1.0 - dy) * (dz);
        }
      if(column1 >= myplan.base_firstcol && column1 <= myplan.base_lastcol)
        {
          flistin[i] += grid[FC(column1, slab_z)] * (1.0 - dx) * (dy) * (1.0 - dz) + grid[FC(column1, slab_zz)] * (1.0 - dx) * (dy) * (dz);
        }

      if(column2 >= myplan.base_firstcol && column2 <= myplan.base_lastcol)
        {
          flistin[i] += grid[FC(column2, slab_z)] * (dx) * (1.0 - dy) * (1.0 - dz) + grid[FC(column2, slab_zz)] * (dx) * (1.0 - dy) * (dz);
        }

      if(column3 >= myplan.base_firstcol && column3 <= myplan.base_lastcol)
        {
          flistin[i] += grid[FC(column3, slab_z)] * (dx) * (dy) * (1.0 - dz) + grid[FC(column3, slab_zz)] * (dx) * (dy) * (dz);
        }
#endif
    }

  /* exchange the potential component data */
  int flag_big = 0, flag_big_all;
  for(i = 0; i < NTask; i++)
    if(Sndpm_count[i] * sizeof(double) > MPI_MESSAGE_SIZELIMIT_IN_BYTES)
      flag_big = 1;

  /* produce a flag if any of the send sizes is above our transfer limit, in this case we will
   * transfer the data in chunks.
   */
  MPI_Allreduce(&flag_big, &flag_big_all, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  /* exchange  data */
  myMPI_Alltoallv(flistin, Rcvpm_count, Rcvpm_offset, flistout, Sndpm_count, Sndpm_offset, sizeof(double), flag_big_all, MPI_COMM_WORLD);


  /* now assign them to the correct particles */
  int multiNtask = roundup_to_multiple_of_cacheline_size(NTask * sizeof(size_t)) / sizeof(size_t);

#pragma omp parallel            /* note: this is *not* a parallel for, instead each threads needs to do the whole loop */
  {
    size_t *send_count = Sndpm_count + get_thread_num() * multiNtask;
    size_t *send_offset = Sndpm_offset + get_thread_num() * multiNtask;

    int j;
    for(j = 0; j < NTask; j++)
      send_count[j] = 0;

    int i;
#pragma omp for schedule(static)
    for(i = 0; i < NumPart; i++)
      {
        MyDouble *pos;

#ifdef CELL_CENTER_GRAVITY
        MyDouble posw[3], xtmp, ytmp, ztmp;
        if(P[i].Type == 0)
          {
            posw[0] = WRAP_X(SphP[i].Center[0]);
            posw[1] = WRAP_Y(SphP[i].Center[1]);
            posw[2] = WRAP_Z(SphP[i].Center[2]);

            pos = posw;
          }
        else
#endif
          pos = P[i].Pos;

        int slab_x = (int) (to_slab_fac * pos[0]);
        int slab_xx = slab_x + 1;

        if(slab_x >= GRIDX)
          slab_x -= GRIDX;

        if(slab_xx >= GRIDX)
          slab_xx -= GRIDX;

#ifndef FFT_COLUMN_BASED
        int task0 = myplan.slab_to_task[slab_x];
        int task1 = myplan.slab_to_task[slab_xx];

        double value = flistout[send_offset[task0] + send_count[task0]++];

        if(task0 != task1)
          value += flistout[send_offset[task1] + send_count[task1]++];
#else
        int slab_y = (int) (to_slab_fac * pos[1]);
        int slab_yy = slab_y + 1;

        if(slab_y >= GRIDY)
          slab_y -= GRIDY;

        if(slab_yy >= GRIDY)
          slab_yy -= GRIDY;

        int column0 = slab_x * GRIDY + slab_y;
        int column1 = slab_x * GRIDY + slab_yy;
        int column2 = slab_xx * GRIDY + slab_y;
        int column3 = slab_xx * GRIDY + slab_yy;

        int task0, task1, task2, task3;

        if(column0 < myplan.pivotcol)
          task0 = column0 / myplan.avg;
        else
          task0 = (column0 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

        if(column1 < myplan.pivotcol)
          task1 = column1 / myplan.avg;
        else
          task1 = (column1 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

        if(column2 < myplan.pivotcol)
          task2 = column2 / myplan.avg;
        else
          task2 = (column2 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

        if(column3 < myplan.pivotcol)
          task3 = column3 / myplan.avg;
        else
          task3 = (column3 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

        double value = flistout[send_offset[task0] + send_count[task0]++];

        if(task1 != task0)
          value += flistout[send_offset[task1] + send_count[task1]++];

        if(task2 != task1 && task2 != task0)
          value += flistout[send_offset[task2] + send_count[task2]++];

        if(task3 != task0 && task3 != task1 && task3 != task2)
          value += flistout[send_offset[task3] + send_count[task3]++];
#endif
        if(dim < 0)
          {
#ifdef EVALPOTENTIAL
            P[i].PM_Potential += value * fac;
#endif
          }
        else
          P[i].GravPM[dim] += value;
      }
  }

  int j;
  /* restore total Sndpm_count */
  for(j = 1; j < MaxThreads; j++)
    for(i = 0; i < NTask; i++)
      Sndpm_count[i] += Sndpm_count[i + j * multiNtask];

  myfree(flistout);
  myfree(flistin);
}
#endif

#ifdef POWERSPECTRUM_IN_POSTPROCESSING
#ifndef HAVE_HDF5
#error POWERSPECTRUM_IN_POSTPROCESSING requires HAVE_HDF5
#endif

static void pmforce_uniform_optimized_prepare_density_powerspectrum(int mode, int *typelist)
{
  mpi_printf( "Starting to compute density field, foldings: %d.\n", mode-1 );
  MPI_Barrier( MPI_COMM_WORLD );
  
  /* We here enlarge NTask such that each thread gets his own cache line for send_count/send_offset.
   * This should hopefully prevent a performance penalty from 'false sharing' for these variables
   */
  int multiNtask = roundup_to_multiple_of_cacheline_size(NTask * sizeof(size_t)) / sizeof(size_t);

  Sndpm_count = (size_t *) mymalloc("Sndpm_count", MaxThreads * multiNtask * sizeof(size_t));
  Sndpm_offset = (size_t *) mymalloc("Sndpm_offset", MaxThreads * multiNtask * sizeof(size_t));
  Rcvpm_count = (size_t *) mymalloc("Rcvpm_count", NTask * sizeof(size_t));
  Rcvpm_offset = (size_t *) mymalloc("Rcvpm_offset", NTask * sizeof(size_t));
  
  /* allocate density field */
  rhogrid = (fft_real *) mymalloc("rhogrid", maxfftsize * sizeof(fft_real));

  /* clear local FFT-mesh density field */
#pragma omp parallel for
  for(large_array_offset ii = 0; ii < maxfftsize; ii++)
    rhogrid[ii] = 0;
  
  power_spec_totmass = 0;
  power_spec_totmass2 = 0;
  
  /* read snapshot and add particles by chunks, filename is stored in All.InputFileName */
  if(All.ICFormat != 3 && All.ICFormat != 1)
    {
      mpi_printf("ICFormat=%d not supported for POWERSPECTRUM_IN_POSTPROCESSING.\n", All.ICFormat);
      endrun();
    }
  
  int num_files = find_files(All.InputFileName);
  
  int nrounds = ((num_files-1) / NTask) + 1;
  
#ifndef POWERSPECTRUM_IN_POSTPROCESSING_ICS
  hid_t hdf5_file = 0;
#endif
  
  char fname[1000];
  
  /* count particles we have to read */
  int npart[NTYPES];
  long long npartall[NTYPES];
  double masses[NTYPES];
  
  for(int type=0; type<NTYPES; type++)
    npart[type] = 0;
  
  for(int r=0; r < nrounds; r++)
    {
      int ThisFile = ThisTask + r * NTask;
      if(ThisFile < num_files)
        {
#ifdef POWERSPECTRUM_IN_POSTPROCESSING_ICS
          if(num_files > 1)
            sprintf(fname, "%s.%d", All.InputFileName, ThisFile);
          else
            sprintf(fname, "%s", All.InputFileName);
          
          FILE *fd = 0;
          if(!(fd = fopen(fname, "r")))
            {
              char buf[500];
              sprintf(buf, "can't open file `%s' for reading initial conditions.\n", fname);
              terminate(buf);
            }
          
          int blksize;
          my_fread( &blksize, sizeof(int), 1, fd );
          read_header_attributes(fd);
          
          for(int type = 0; type < NTYPES; type++) {
            npart[type] += header.npart[type];
            masses[type] = header.mass[type];
          }
          
          All.Time = All.TimeBegin;
          
          fclose(fd);
#else
          if(num_files > 1)
            sprintf(fname, "%s.%d.hdf5", All.InputFileName, ThisFile);
          else
            sprintf(fname, "%s.hdf5", All.InputFileName);
      
          hdf5_file = my_H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
          if(hdf5_file < 0)
            terminate("cannot read initial conditions file %s", fname);
      
          /* get total number of particles we have to read */
          hid_t hdf5_headergrp, hdf5_attribute;
          hssize_t scalar_attr_dim = 1;
          hssize_t vector_attr_dim = NTYPES;
      
          hdf5_headergrp = my_H5Gopen(hdf5_file, "/Header");
          
          int npart_load[NTYPES];
          hdf5_attribute = my_H5Aopen_name(hdf5_headergrp, "NumPart_ThisFile");
          my_H5Aread(hdf5_attribute, H5T_NATIVE_INT, npart_load, "NumPart_ThisFile", vector_attr_dim);
          my_H5Aclose(hdf5_attribute, "NumPart_ThisFile");
          for(int type=0; type<NTYPES; type++)
            npart[type] += npart_load[type];
      
          hdf5_attribute = my_H5Aopen_name(hdf5_headergrp, "MassTable");
          my_H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, masses, "MassTable", vector_attr_dim);
          my_H5Aclose(hdf5_attribute, "MassTable");
      
          hdf5_attribute = my_H5Aopen_name(hdf5_headergrp, "Time");
          my_H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.Time, "Time", scalar_attr_dim);
          my_H5Aclose(hdf5_attribute, "Time");
      
          my_H5Gclose(hdf5_headergrp, "/Header");
          my_H5Fclose(hdf5_file, fname);
#endif
        }
    }
     
  for(int type = 0; type < NTYPES; type++)
    if(npart[type] > 0 && typelist[type] == 0)
      npart[type] = 0;
#ifdef TRACER_MC
  npart[3] = 0; // stores the shitty tracers
#endif

  sumup_large_ints(NTYPES, npart, npartall);
  MPI_Bcast(masses, NTYPES, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  MPI_Bcast(&All.Time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  set_softenings();
  
  if (ThisTask == 0)
    {
      printf( "Time: %g\n", All.Time );
      for(int type=0; type<NTYPES; type++)
        {
          if(typelist[type])
            printf( "Type %d: npart=%lld, mass=%g\n", type, npartall[type], masses[type] );
        }
    }
  
  MPI_Barrier( MPI_COMM_WORLD );
  
  power_spec_totnumpart = 0;
  for(int type=0; type<NTYPES; type++)
    power_spec_totnumpart += npartall[type];
  
  double mtot = 0;
  double m2tot = 0;

#ifdef POWERSPECTRUM_IN_POSTPROCESSING_ICS
  FILE *fd = 0;
#ifdef NTYPES_ICS
  int headersize = sizeof(header_ICs);
#else
  int headersize = sizeof(header);
#endif
#endif

  /* read masses and positions, and add the particles to density grid */
  for(int r=0; r < nrounds; r++)
    {
#ifdef POWERSPECTRUM_IN_POSTPROCESSING_ICS
      int npartcum[NTYPES+1];
      for(int type=0; type<=NTYPES; type++)
        npartcum[type] = 0;
#endif
      
      int ThisFile = ThisTask + r * NTask;
      if(ThisFile < num_files)
        {
#ifdef POWERSPECTRUM_IN_POSTPROCESSING_ICS
          if(num_files > 1)
            sprintf(fname, "%s.%d", All.InputFileName, ThisFile);
          else
            sprintf(fname, "%s", All.InputFileName);
          
          fd = fopen(fname, "r");
          
          {int blksize; my_fread( &blksize, sizeof(int), 1, fd );}
          read_header_attributes(fd);
          {int blksize; my_fread( &blksize, sizeof(int), 1, fd );}
          
          for(int type = 0; type < NTYPES; type++) {
            npart[type] = header.npart[type];
            npartcum[type+1] = npartcum[type] + npart[type];
          }
#else
          if(num_files > 1)
            sprintf(fname, "%s.%d.hdf5", All.InputFileName, ThisFile);
          else
            sprintf(fname, "%s.hdf5", All.InputFileName);
      
          hdf5_file = my_H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
          if(hdf5_file < 0)
            terminate("cannot read initial conditions file %s", fname);
          
          /* get total number of particles we have to read */
          hid_t hdf5_headergrp, hdf5_attribute;
          hssize_t vector_attr_dim = NTYPES;
      
          hdf5_headergrp = my_H5Gopen(hdf5_file, "/Header");
          hdf5_attribute = my_H5Aopen_name(hdf5_headergrp, "NumPart_ThisFile");
          my_H5Aread(hdf5_attribute, H5T_NATIVE_INT, npart, "NumPart_ThisFile", vector_attr_dim);
          my_H5Aclose(hdf5_attribute, "NumPart_ThisFile");
          my_H5Gclose(hdf5_headergrp, "/Header");
#endif
          for(int type = 0; type < NTYPES; type++)
            if(npart[type] > 0 && typelist[type] == 0)
              npart[type] = 0;
          
#ifdef TRACER_MC
          npart[3] = 0; // stores the shitty tracers
#endif
        }
      else
        {
          for(int type=0; type<NTYPES; type++)
            npart[type] = 0;
        }
      
      if(nrounds > 1)
        mpi_printf( "Round %d of %d, reading files %d to %d.\n", r+1, nrounds, r*NTask, imin((r+1)*NTask-1, num_files) );
      
      for(int type=0; type < NTYPES; type++)
        {
          if(typelist[type] == 0)
            continue;
          
          mpi_printf( "Loading type=%d.\n", type );
          MPI_Barrier( MPI_COMM_WORLD );
      
          int totcount = 0;
      
          int done = 0;
          int iter = 0;
          while(!done)
            {
              double *pmass, *ppos;
              pmass = ppos = 0;
          
              int count = 0;
          
              if(ThisFile < num_files && npart[type] > totcount)
                {
                  count = dmin( npart[type]-totcount, 1000000 );
              
                  pmass = mymalloc("pmass", count * sizeof(double) );
                  ppos  = mymalloc("ppos",  count * sizeof(double) * 3 );
              
#ifdef POWERSPECTRUM_IN_POSTPROCESSING_ICS
#ifdef INPUT_IN_DOUBLEPRECISION
                  fseek( fd, 4 + headersize + 4 + 4 + npartcum[type] * sizeof(double) * 3 + totcount * sizeof(double) * 3, SEEK_SET );
                  my_fread( ppos, count * sizeof(double) * 3, 1, fd);
                  
                  if(masses[type] == 0)
                    {
                      fseek( fd, 4 + headersize + 4 + 4 + npartcum[NTYPES] * sizeof(double) * 6 + npartcum[NTYPES] * sizeof(int) * 6 + totcount * sizeof(double), SEEK_SET );
                      my_fread( pmass, count * sizeof(float), 1, fd);
                    }
                  else
                    {
                      if(iter == 0) mpi_printf("Using header mass value of %g for type %d.\n", masses[type], type );
                      for(large_array_offset p=0; p<count; p++)
                        pmass[p] = masses[type];
                    }
#else
                  fseek( fd, 4 + headersize + 4 + 4 + npartcum[type] * sizeof(float) * 3 + totcount * sizeof(float) * 3, SEEK_SET );
                  float *temp = mymalloc("temp", count * sizeof(float) * 3);
                  my_fread( temp, sizeof(float), count * 3, fd);
                  
                  for(large_array_offset p=0; p<count*3; p++)
                    ppos[p] = temp[p];
                  myfree( temp );
                  
                  if(masses[type] == 0)
                    {
                      fseek( fd, 4 + headersize + 4 + 4 + npartcum[NTYPES] * sizeof(float) * 6 + npartcum[NTYPES] * sizeof(int) * 6 + totcount * sizeof(float), SEEK_SET );                      
                      temp = mymalloc("temp", count * sizeof(float));
                      my_fread( temp, count * sizeof(float), 1, fd);

                      for(large_array_offset p=0; p<count; p++)
                        pmass[p] = temp[p];
                      myfree( temp );
                    }
                  else
                    {
                      if(iter == 0) mpi_printf("Using header mass value of %g for type %d.\n", masses[type], type );
                      for(large_array_offset p=0; p<count; p++)
                        pmass[p] = masses[type];
                    }
#endif

                  for(large_array_offset p = 0; p < count; p++)
                    {
                      mtot += pmass[p];
                      m2tot += pmass[p] * pmass[p];
                    }
#else
                  char buf[1000];
                  sprintf(buf, "/PartType%d", type);
                  hid_t hdf5_grp = my_H5Gopen(hdf5_file, buf);
                  
                  hid_t hdf5_dataset_masses;
                  if(type == 5)
                    {
                      hdf5_dataset_masses = my_H5Dopen_if_existing( hdf5_grp, "BH_Mass" );
              
                      if(hdf5_dataset_masses < 0)
                        {
                          if(iter == 0) mpi_printf( "Did not find Dataset BH_Mass for type=5, using Masses instead.\n" );
                          hdf5_dataset_masses = my_H5Dopen_if_existing( hdf5_grp, "Masses" );
                        }
                      else
                        {
                          if(iter == 0) mpi_printf( "Using Dataset BH_Mass for type=5.\n" );
                        }
                    }
                  else
                    {
                      hdf5_dataset_masses = my_H5Dopen_if_existing( hdf5_grp, "Masses" );
                    }
              
                  if(hdf5_dataset_masses < 0)
                    {
                      if(iter == 0) mpi_printf("No Masses found for type %d in file %d, using header value of %g instead.\n", type, ThisFile, masses[type] );
                      for(large_array_offset p=0; p<count; p++)
                        pmass[p] = masses[type];
                    }
                  else
                    {
                      hid_t hdf5_space_masses = H5Dget_space( hdf5_dataset_masses );
                      if( H5Sget_simple_extent_ndims( hdf5_space_masses ) != 1 )
                        terminate( "Masses should be a one-dimensional array." )
          
                      hsize_t dims_masses;
                      H5Sget_simple_extent_dims( hdf5_space_masses, &dims_masses, NULL );
                      if(dims_masses != npart[type])
                        terminate( "Unexpected length of mass array." );
          
                      hid_t type_masses = H5Dget_type( hdf5_dataset_masses );
                      hid_t type_masses_mem = H5Tget_native_type( type_masses, H5T_DIR_DESCEND );
                  
                      hsize_t data_offset = totcount;
                      hsize_t data_count = count;
                      H5Sselect_hyperslab(hdf5_space_masses, H5S_SELECT_SET, &data_offset, NULL, &data_count, NULL);
                  
                      hsize_t memory_dim = count;
                      hid_t memspace = H5Screate_simple(1, &memory_dim, NULL);
                  
                      if(H5Tequal(type_masses_mem,H5T_NATIVE_DOUBLE))
                        {
                          //my_H5Dread( hdf5_dataset_masses, H5T_NATIVE_DOUBLE, hdf5_space_masses, hdf5_space_masses, H5P_DEFAULT, pmass, "Masses" );
                          my_H5Dread( hdf5_dataset_masses, H5T_NATIVE_DOUBLE, memspace, hdf5_space_masses, H5P_DEFAULT, pmass, "Masses" );
                        }
                      else if(H5Tequal(type_masses_mem,H5T_NATIVE_FLOAT))
                        {
                          float *temp = mymalloc("temp", count * sizeof(float) );
                          //my_H5Dread( hdf5_dataset_masses, H5T_NATIVE_FLOAT, hdf5_space_masses, hdf5_space_masses, H5P_DEFAULT, temp, "Masses" );
                          my_H5Dread( hdf5_dataset_masses, H5T_NATIVE_FLOAT, memspace, hdf5_space_masses, H5P_DEFAULT, temp, "Masses" );
                          for(large_array_offset p=0; p<count; p++)
                            pmass[p] = temp[p];
                          myfree( temp );
                        }
                      else
                        terminate( "Masses is neither of type float nor double." );
          
                      H5Tclose( type_masses );
                      H5Tclose( type_masses_mem );
                      my_H5Dclose(hdf5_dataset_masses, "Masses");
                    }
          
                  for(large_array_offset p = 0; p < count; p++)
                    {
                      mtot += pmass[p];
                      m2tot += pmass[p] * pmass[p];
                    }
          
                  // Check if center of mass exists!
                  hid_t hdf5_dataset_coordinates;
          
                  int wrapping = 0;
                  if(type == 0)
                    {
                      hdf5_dataset_coordinates = my_H5Dopen_if_existing( hdf5_grp, "CenterOfMass" );
              
                      if(hdf5_dataset_coordinates < 0)
                        {
                          if(iter == 0) mpi_printf( "Did not find Dataset CenterOfMass for type=0, using Coordinates instead.\n" );
                          hdf5_dataset_coordinates = my_H5Dopen_if_existing( hdf5_grp, "Coordinates" );
                        }
                      else
                        {
                          if(iter == 0) mpi_printf( "Using Dataset CenterOfMass for type=0.\n" );
                          wrapping = 1;
                        }
                    }
                  else
                    {
                      hdf5_dataset_coordinates = my_H5Dopen_if_existing( hdf5_grp, "Coordinates" );
                    }
          
                  if(hdf5_dataset_coordinates < 0)
                    terminate( "Dataset Coordinates not found." );
          
                  hid_t hdf5_space_coordinates = H5Dget_space( hdf5_dataset_coordinates );
                  if( H5Sget_simple_extent_ndims( hdf5_space_coordinates ) != 2 )
                    terminate( "Coordinates should be a two-dimensional array." )
          
                  hsize_t dims_coordinates[2];
                  H5Sget_simple_extent_dims( hdf5_space_coordinates, dims_coordinates, NULL );
                  if(dims_coordinates[0] != npart[type] || dims_coordinates[1] != 3)
                    terminate( "Unexpected length of coordinates array." );
          
                  hid_t type_coordinates = H5Dget_type( hdf5_dataset_coordinates );
                  hid_t type_coordinates_mem = H5Tget_native_type( type_coordinates, H5T_DIR_DESCEND );
              
                  hsize_t data_offset[2];
                  data_offset[0] = totcount;
                  data_offset[1] = 0;
                  hsize_t data_count[2];
                  data_count[0] = count;
                  data_count[1] = 3;
                  H5Sselect_hyperslab(hdf5_space_coordinates, H5S_SELECT_SET, data_offset, NULL, data_count, NULL);
              
                  hsize_t memory_dim[2];
                  memory_dim[0] = count;
                  memory_dim[1] = 3;
                  hid_t memspace = H5Screate_simple(2, memory_dim, NULL);
          
                  if(H5Tequal(type_coordinates_mem,H5T_NATIVE_DOUBLE))
                    {
                      //my_H5Dread( hdf5_dataset_coordinates, H5T_NATIVE_DOUBLE, hdf5_space_coordinates, hdf5_space_coordinates, H5P_DEFAULT, ppos, "Coordinates" );
                      my_H5Dread( hdf5_dataset_coordinates, H5T_NATIVE_DOUBLE, memspace, hdf5_space_coordinates, H5P_DEFAULT, ppos, "Coordinates" );
                    }
                  else if(H5Tequal(type_coordinates_mem,H5T_NATIVE_FLOAT))
                    {
                      float *temp = mymalloc("temp", count * sizeof(float) * 3);
                      //my_H5Dread( hdf5_dataset_coordinates, H5T_NATIVE_FLOAT, hdf5_space_coordinates, hdf5_space_coordinates, H5P_DEFAULT, temp, "Coordinates" );
                      my_H5Dread( hdf5_dataset_coordinates, H5T_NATIVE_FLOAT, memspace, hdf5_space_coordinates, H5P_DEFAULT, temp, "Coordinates" );
                      for(large_array_offset p=0; p<count*3; p++)
                        ppos[p] = temp[p];
                      myfree( temp );
                    }
                  else
                    terminate( "Coordinates is neither of type float nor double." );
          
                  H5Tclose( type_coordinates );
                  H5Tclose( type_coordinates_mem );
                  my_H5Dclose(hdf5_dataset_coordinates, "Coordinates");
          
                  if(wrapping)
                    {
                      double xtmp, ytmp, ztmp;
                      for(large_array_offset p = 0; p < count; p++)
                        {
                          ppos[p*3+0] = WRAP_X(ppos[p*3+0]);
                          ppos[p*3+1] = WRAP_Y(ppos[p*3+1]);
                          ppos[p*3+2] = WRAP_Z(ppos[p*3+2]);
                        }
                    }
            
                  my_H5Gclose(hdf5_grp, buf);
#endif
                }
      
              large_array_offset i, j;

              double to_slab_fac = PMGRID / All.BoxSize;

              if(mode == 2)
                to_slab_fac *= POWERSPEC_FOLDFAC;
              if(mode == 3)
                to_slab_fac *= POWERSPEC_FOLDFAC * POWERSPEC_FOLDFAC;

              /* determine the slabs/columns each particles accesses */
#pragma omp parallel private(j)
              {
                size_t *send_count = Sndpm_count + get_thread_num() * multiNtask;

                /* each threads needs to do theloop to clear its send_count[] array */
                for(j = 0; j < NTask; j++)
                  send_count[j] = 0;

#pragma omp for schedule(static) private(i)
                for(i = 0; i < count; i++)
                  {
                    MyDouble *pos = &(ppos[i*3]);

                    int slab_x = (int) (to_slab_fac * pos[0]);
                    int slab_xx = slab_x + 1;

                    if(mode >= 2)
                      {
                        slab_x %= GRIDX;
                        slab_xx %= GRIDX;
                      }
                    else
                      {
                        if(slab_x >= GRIDX)
                          slab_x -= GRIDX;

                        if(slab_xx >= GRIDX)
                          slab_xx -= GRIDX;
                      }

#ifndef FFT_COLUMN_BASED
                    int task0 = myplan.slab_to_task[slab_x];
                    int task1 = myplan.slab_to_task[slab_xx];

                    send_count[task0]++;
                    if(task0 != task1)
                      send_count[task1]++;
#else
                    int slab_y = (int) (to_slab_fac * pos[1]);
                    int slab_yy = slab_y + 1;

                    if(mode >= 2)
                      {
                        slab_y %= GRIDY;
                        slab_yy %= GRIDY;
                      }
                    else
                      {
                        if(slab_y >= GRIDY)
                          slab_y -= GRIDY;

                        if(slab_yy >= GRIDY)
                          slab_yy -= GRIDY;
                      }

                    int column0 = slab_x * GRIDY + slab_y;
                    int column1 = slab_x * GRIDY + slab_yy;
                    int column2 = slab_xx * GRIDY + slab_y;
                    int column3 = slab_xx * GRIDY + slab_yy;

                    int task0, task1, task2, task3;

                    if(column0 < myplan.pivotcol)
                      task0 = column0 / myplan.avg;
                    else
                      task0 = (column0 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

                    if(column1 < myplan.pivotcol)
                      task1 = column1 / myplan.avg;
                    else
                      task1 = (column1 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

                    if(column2 < myplan.pivotcol)
                      task2 = column2 / myplan.avg;
                    else
                      task2 = (column2 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

                    if(column3 < myplan.pivotcol)
                      task3 = column3 / myplan.avg;
                    else
                      task3 = (column3 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

                    send_count[task0]++;
                    if(task1 != task0)
                      send_count[task1]++;
                    if(task2 != task1 && task2 != task0)
                      send_count[task2]++;
                    if(task3 != task0 && task3 != task1 && task3 != task2)
                      send_count[task3]++;
#endif
                  }
              }


              /* collect thread-specific offset table and collect the results from the other threads */
              for(i = 0, Sndpm_offset[0] = 0; i < NTask; i++)
                for(j = 0; j < MaxThreads; j++)
                  {
                    large_array_offset ind_prev, ind = j * multiNtask + i;
                    if(ind > 0)
                      {
                        if(j == 0)
                          ind_prev = (MaxThreads - 1) * multiNtask + i - 1;
                        else
                          ind_prev = ind - multiNtask;

                        Sndpm_offset[ind] = Sndpm_offset[ind_prev] + Sndpm_count[ind_prev];
                      }
                  }

              for(j = 1; j < MaxThreads; j++)
                for(i = 0; i < NTask; i++)
                  Sndpm_count[i] += Sndpm_count[i + j * multiNtask];

              MPI_Alltoall(Sndpm_count, sizeof(size_t), MPI_BYTE, Rcvpm_count, sizeof(size_t), MPI_BYTE, MPI_COMM_WORLD);

              for(j = 0, nimport = 0, nexport = 0, Rcvpm_offset[0] = 0, Sndpm_offset[0] = 0; j < NTask; j++)
                {
                  nexport += Sndpm_count[j];
                  nimport += Rcvpm_count[j];

                  if(j > 0)
                    {
                      Sndpm_offset[j] = Sndpm_offset[j - 1] + Sndpm_count[j - 1];
                      Rcvpm_offset[j] = Rcvpm_offset[j - 1] + Rcvpm_count[j - 1];
                    }
                }

              /* allocate import and export buffer */
              partin = (struct partbuf *) mymalloc("partin", nimport * sizeof(struct partbuf));
              partout = (struct partbuf *) mymalloc("partout", nexport * sizeof(struct partbuf));

#pragma omp parallel private(j)
              {
                size_t *send_count = Sndpm_count + get_thread_num() * multiNtask;
                size_t *send_offset = Sndpm_offset + get_thread_num() * multiNtask;

                for(j = 0; j < NTask; j++)
                  send_count[j] = 0;

                /* fill export buffer */
#pragma omp for schedule(static)
                for(i = 0; i < count; i++)
                  {
                    MyDouble *pos = &ppos[i*3];

                    int slab_x = (int) (to_slab_fac * pos[0]);
                    int slab_xx = slab_x + 1;

                    if(mode >= 2)
                      {
                        slab_x %= GRIDX;
                        slab_xx %= GRIDX;
                      }
                    else
                      {
                        if(slab_x >= GRIDX)
                          slab_x -= GRIDX;

                        if(slab_xx >= GRIDX)
                          slab_xx -= GRIDX;
                      }

#ifndef FFT_COLUMN_BASED
                    int task0 = myplan.slab_to_task[slab_x];
                    int task1 = myplan.slab_to_task[slab_xx];

                    size_t ind0 = send_offset[task0] + send_count[task0]++;
                    partout[ind0].Mass = pmass[i];
                    for(j = 0; j < 3; j++)
                      partout[ind0].Pos[j] = pos[j];

                    if(task0 != task1)
                      {
                        size_t ind1 = send_offset[task1] + send_count[task1]++;
                        partout[ind1].Mass = pmass[i];
                        for(j = 0; j < 3; j++)
                          partout[ind1].Pos[j] = pos[j];
                      }
#else
                    int slab_y = (int) (to_slab_fac * pos[1]);
                    int slab_yy = slab_y + 1;

                    if(mode >= 2)
                      {
                        slab_y %= GRIDY;
                        slab_yy %= GRIDY;
                      }
                    else
                      {
                        if(slab_y >= GRIDY)
                          slab_y -= GRIDY;

                        if(slab_yy >= GRIDY)
                          slab_yy -= GRIDY;
                      }

                    int column0 = slab_x * GRIDY + slab_y;
                    int column1 = slab_x * GRIDY + slab_yy;
                    int column2 = slab_xx * GRIDY + slab_y;
                    int column3 = slab_xx * GRIDY + slab_yy;

                    int task0, task1, task2, task3;

                    if(column0 < myplan.pivotcol)
                      task0 = column0 / myplan.avg;
                    else
                      task0 = (column0 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

                    if(column1 < myplan.pivotcol)
                      task1 = column1 / myplan.avg;
                    else
                      task1 = (column1 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

                    if(column2 < myplan.pivotcol)
                      task2 = column2 / myplan.avg;
                    else
                      task2 = (column2 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

                    if(column3 < myplan.pivotcol)
                      task3 = column3 / myplan.avg;
                    else
                      task3 = (column3 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

                    size_t ind0 = send_offset[task0] + send_count[task0]++;
                    partout[ind0].Mass = pmass[i];
                    for(j = 0; j < 3; j++)
                      partout[ind0].Pos[j] = pos[j];

                    if(task1 != task0)
                      {
                        size_t ind1 = send_offset[task1] + send_count[task1]++;
                        partout[ind1].Mass = pmass[i];
                        for(j = 0; j < 3; j++)
                          partout[ind1].Pos[j] = pos[j];
                      }
                    if(task2 != task1 && task2 != task0)
                      {
                        size_t ind2 = send_offset[task2] + send_count[task2]++;
                        partout[ind2].Mass = pmass[i];
                        for(j = 0; j < 3; j++)
                          partout[ind2].Pos[j] = pos[j];
                      }
                    if(task3 != task0 && task3 != task1 && task3 != task2)
                      {
                        size_t ind3 = send_offset[task3] + send_count[task3]++;
                        partout[ind3].Mass = pmass[i];
                        for(j = 0; j < 3; j++)
                          partout[ind3].Pos[j] = pos[j];
                      }
#endif
                  }
              }

              /* collect the send_count[] results from the other threads */
              for(j = 1; j < MaxThreads; j++)
                for(i = 0; i < NTask; i++)
                  Sndpm_count[i] += Sndpm_count[i + j * multiNtask];

              int flag_big = 0, flag_big_all;
              for(i = 0; i < NTask; i++)
                if(Sndpm_count[i] * sizeof(struct partbuf) > MPI_MESSAGE_SIZELIMIT_IN_BYTES)
                  flag_big = 1;

              /* produce a flag if any of the send sizes is above our transfer limit, in this case we will
               * transfer the data in chunks.
               */
              MPI_Allreduce(&flag_big, &flag_big_all, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

              /* exchange particle data */
              myMPI_Alltoallv(partout, Sndpm_count, Sndpm_offset, partin, Rcvpm_count, Rcvpm_offset, sizeof(struct partbuf), flag_big_all, MPI_COMM_WORLD);

              myfree(partout);
      
#ifndef FFT_COLUMN_BASED
              /* bin particle data onto mesh, in multi-threaded fashion */
#pragma omp parallel private(i)
              {
                int tid = get_thread_num();

                int first_y, count_y;
                subdivide_evenly(GRIDY, MaxThreads, tid, &first_y, &count_y);
                int last_y = first_y + count_y - 1;

                for(i = 0; i < nimport; i++)
                  {
                    int slab_y = (int) (to_slab_fac * partin[i].Pos[1]);
                    int slab_yy = slab_y + 1;
                    double dy = to_slab_fac * partin[i].Pos[1] - slab_y;

                    if(mode >= 2)
                      {
                        slab_y %= GRIDY;
                        slab_yy %= GRIDY;
                      }
                    else
                      {
                        if(slab_y >= GRIDY)
                          slab_y -= GRIDY;

                        if(slab_yy >= GRIDY)
                          slab_yy -= GRIDY;
                      }

                    int flag_slab_y, flag_slab_yy;

                    if(slab_y >= first_y && slab_y <= last_y)
                      flag_slab_y = 1;
                    else
                      flag_slab_y = 0;

                    if(slab_yy >= first_y && slab_yy <= last_y)
                      flag_slab_yy = 1;
                    else
                      flag_slab_yy = 0;

                    if(flag_slab_y || flag_slab_yy)
                      {
                        double mass = partin[i].Mass;

                        int slab_x = (int) (to_slab_fac * partin[i].Pos[0]);
                        int slab_z = (int) (to_slab_fac * partin[i].Pos[2]);
                        int slab_xx = slab_x + 1;
                        int slab_zz = slab_z + 1;

                        double dx = to_slab_fac * partin[i].Pos[0] - slab_x;
                        double dz = to_slab_fac * partin[i].Pos[2] - slab_z;

                        if(mode >= 2)
                          {
                            slab_x %= GRIDX;
                            slab_z %= GRIDZ;
                            slab_xx %= GRIDX;
                            slab_zz %= GRIDZ;
                          }
                        else
                          {
                            if(slab_x >= GRIDX)
                              slab_x -= GRIDX;
                            if(slab_z >= GRIDZ)
                              slab_z -= GRIDZ;

                            if(slab_xx >= GRIDX)
                              slab_xx -= GRIDX;
                            if(slab_zz >= GRIDZ)
                              slab_zz -= GRIDZ;
                          }

                        int flag_slab_x, flag_slab_xx;

                        if(myplan.slab_to_task[slab_x] == ThisTask)
                          {
                            slab_x -= myplan.first_slab_x_of_task[ThisTask];
                            flag_slab_x = 1;
                          }
                        else
                          flag_slab_x = 0;

                        if(myplan.slab_to_task[slab_xx] == ThisTask)
                          {
                            slab_xx -= myplan.first_slab_x_of_task[ThisTask];
                            flag_slab_xx = 1;
                          }
                        else
                          flag_slab_xx = 0;

                        if(flag_slab_x)
                          {
                            if(flag_slab_y)
                              {
                                rhogrid[FI(slab_x, slab_y, slab_z)] += (mass * (1.0 - dx) * (1.0 - dy) * (1.0 - dz));
                                rhogrid[FI(slab_x, slab_y, slab_zz)] += (mass * (1.0 - dx) * (1.0 - dy) * (dz));
                              }

                            if(flag_slab_yy)
                              {
                                rhogrid[FI(slab_x, slab_yy, slab_z)] += (mass * (1.0 - dx) * (dy) * (1.0 - dz));
                                rhogrid[FI(slab_x, slab_yy, slab_zz)] += (mass * (1.0 - dx) * (dy) * (dz));
                              }
                          }

                        if(flag_slab_xx)
                          {
                            if(flag_slab_y)
                              {
                                rhogrid[FI(slab_xx, slab_y, slab_z)] += (mass * (dx) * (1.0 - dy) * (1.0 - dz));
                                rhogrid[FI(slab_xx, slab_y, slab_zz)] += (mass * (dx) * (1.0 - dy) * (dz));
                              }

                            if(flag_slab_yy)
                              {
                                rhogrid[FI(slab_xx, slab_yy, slab_z)] += (mass * (dx) * (dy) * (1.0 - dz));
                                rhogrid[FI(slab_xx, slab_yy, slab_zz)] += (mass * (dx) * (dy) * (dz));
                              }
                          }
                      }
                  }
              }

#else

              struct data_cols
              {
                int col0, col1, col2, col3;
                double dx, dy;
              } *aux;

              aux = mymalloc("aux", nimport * sizeof(struct data_cols));

#pragma omp parallel for
              for(i = 0; i < nimport; i++)
                {
                  int slab_x = (int) (to_slab_fac * partin[i].Pos[0]);
                  int slab_xx = slab_x + 1;

                  int slab_y = (int) (to_slab_fac * partin[i].Pos[1]);
                  int slab_yy = slab_y + 1;

                  aux[i].dx = to_slab_fac * partin[i].Pos[0] - slab_x;
                  aux[i].dy = to_slab_fac * partin[i].Pos[1] - slab_y;

                  if(mode >= 2)
                    {
                      slab_x %= GRIDX;
                      slab_xx %= GRIDX;
                      slab_y %= GRIDY;
                      slab_yy %= GRIDY;
                    }
                  else
                    {
                      if(slab_x >= GRIDX)
                        slab_x -= GRIDX;
                      if(slab_xx >= GRIDX)
                        slab_xx -= GRIDX;

                      if(slab_y >= GRIDY)
                        slab_y -= GRIDY;
                      if(slab_yy >= GRIDY)
                        slab_yy -= GRIDY;
                    }

                  aux[i].col0 = slab_x * GRIDY + slab_y;
                  aux[i].col1 = slab_x * GRIDY + slab_yy;
                  aux[i].col2 = slab_xx * GRIDY + slab_y;
                  aux[i].col3 = slab_xx * GRIDY + slab_yy;
                }

#pragma omp parallel private(i)
              {
                int tid = get_thread_num();

                int first_col, last_col, count_col;
                subdivide_evenly(myplan.base_ncol, MaxThreads, tid, &first_col, &count_col);
                last_col = first_col + count_col - 1;
                first_col += myplan.base_firstcol;
                last_col += myplan.base_firstcol;

                for(i = 0; i < nimport; i++)
                  {
                    int flag0, flag1, flag2, flag3;
                    int col0 = aux[i].col0;
                    int col1 = aux[i].col1;
                    int col2 = aux[i].col2;
                    int col3 = aux[i].col3;

                    if(col0 >= first_col && col0 <= last_col)
                      flag0 = 1;
                    else
                      flag0 = 0;

                    if(col1 >= first_col && col1 <= last_col)
                      flag1 = 1;
                    else
                      flag1 = 0;

                    if(col2 >= first_col && col2 <= last_col)
                      flag2 = 1;
                    else
                      flag2 = 0;

                    if(col3 >= first_col && col3 <= last_col)
                      flag3 = 1;
                    else
                      flag3 = 0;

                    if(flag0 || flag1 || flag2 || flag3)
                      {
                        double mass = partin[i].Mass;

                        double dx = aux[i].dx;
                        double dy = aux[i].dy;

                        int slab_z = (int) (to_slab_fac * partin[i].Pos[2]);
                        int slab_zz = slab_z + 1;

                        double dz = to_slab_fac * partin[i].Pos[2] - slab_z;

                        if(mode >= 2)
                          {
                            slab_z %= GRIDZ;
                            slab_zz %= GRIDZ;
                          }
                        else
                          {
                            if(slab_z >= GRIDZ)
                              slab_z -= GRIDZ;

                            if(slab_zz >= GRIDZ)
                              slab_zz -= GRIDZ;
                          }

                        if(flag0)
                          {
                            rhogrid[FC(col0, slab_z)] += (mass * (1.0 - dx) * (1.0 - dy) * (1.0 - dz));
                            rhogrid[FC(col0, slab_zz)] += (mass * (1.0 - dx) * (1.0 - dy) * (dz));
                          }

                        if(flag1)
                          {
                            rhogrid[FC(col1, slab_z)] += (mass * (1.0 - dx) * (dy) * (1.0 - dz));
                            rhogrid[FC(col1, slab_zz)] += (mass * (1.0 - dx) * (dy) * (dz));
                          }

                        if(flag2)
                          {
                            rhogrid[FC(col2, slab_z)] += (mass * (dx) * (1.0 - dy) * (1.0 - dz));
                            rhogrid[FC(col2, slab_zz)] += (mass * (dx) * (1.0 - dy) * (dz));
                          }

                        if(flag3)
                          {
                            rhogrid[FC(col3, slab_z)] += (mass * (dx) * (dy) * (1.0 - dz));
                            rhogrid[FC(col3, slab_zz)] += (mass * (dx) * (dy) * (dz));
                          }
                      }
                  }
              }

              myfree(aux);

#endif
              myfree(partin);
      
              if(ppos)
                myfree( ppos );
              if(pmass)
                myfree( pmass );
          
              /* update done */
              totcount += count;
          
              int localdone = (totcount == npart[type]);          
              MPI_Allreduce(&localdone, &done, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
          
              long long totcountall;
              sumup_large_ints(1, &totcount, &totcountall);
              if(!done || iter > 0) mpi_printf( "iter=%d, %lld particles done.\n", iter, totcountall );
              iter++;
            }
        }
    
      if(ThisFile < num_files)
        {
#ifdef POWERSPECTRUM_IN_POSTPROCESSING_ICS
          fclose(fd);
#else
          my_H5Fclose(hdf5_file, fname);
#endif
        }
    }
  
  MPI_Allreduce(&mtot, &power_spec_totmass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&m2tot, &power_spec_totmass2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  double rho_min = 1e200;
  double rho_max =    0.;
  
  for(large_array_offset i=0; i<maxfftsize; i++)
    {
      if(rhogrid[i] > rho_max)
        rho_max = rhogrid[i];
      if(rhogrid[i] < rho_min)
        rho_min = rhogrid[i];
      
      if(rhogrid[i] < 0)
        printf( "Task %d: cell i=%lld rho=%g\n.", ThisTask, (long long)i, rhogrid[i] );
    }
  
  double rho_min_all, rho_max_all;
  MPI_Allreduce(&rho_min, &rho_min_all, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&rho_max, &rho_max_all, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  
  double massMesh = 0;
#pragma omp parallel for
    for(large_array_offset ii = 0; ii < maxfftsize; ii++)
      massMesh += rhogrid[ii];
  
  double massMeshTot;
  MPI_Allreduce(&massMesh, &massMeshTot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  mpi_printf( "Total mass: %g rho min/max %g/%g, total mass on mesh: %g.\n", power_spec_totmass, rho_min_all, rho_max_all, massMeshTot );
}
#endif

/*! Calculates the long-range periodic force given the particle positions
 *  using the PM method.  The force is Gaussian filtered with Asmth, given in
 *  mesh-cell units. We carry out a CIC charge assignment, and compute the
 *  potential by fast Fourier transform methods. The potential is finite-differenced
 *  using a 4-point finite differencing formula, and the forces are
 *  interpolated tri-linearly to the particle positions. The CIC kernel is
 *  deconvolved.
 *
 *  For mode=0, normal force calculation, mode=1, only density field construction
 *  for a power spectrum calculation. In the later case, typelist flags the particle
 *  types that should be included in the density field.
 */
void pmforce_periodic(int mode, int *typelist)
{
  int x, y, z, xx, yy, zz;

  double tstart = second();

  if(mode == 0)
    mpi_printf("PM-PERIODIC: Starting periodic PM calculation.  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

#ifndef NUMPART_PER_TASK_LARGE
  if((((long long) NumPart) << 3) >= (((long long) 1) << 31))
    terminate("We are dealing with a too large particle number per MPI rank - enabling NUMPART_PER_TASK_LARGE might help.");
#endif

  double asmth2 = All.Asmth[0] * All.Asmth[0];
  double d = All.BoxSize / PMGRID;
  double dhalf = 0.5 * d;

#ifdef GRAVITY_TALLBOX
  double fac = All.G / (((double) GRIDX) * GRIDY * GRIDZ);      /* to get potential  */
#else
  double fac = 4 * M_PI * All.G / (pow(All.BoxSize, 3) * STRETCHX * STRETCHY * STRETCHZ);       /* to get potential  */
#endif

  fac *= 1 / (2 * d);           /* for finite differencing */

#ifdef POWERSPECTRUM_IN_POSTPROCESSING
  pmforce_uniform_optimized_prepare_density_powerspectrum(mode, typelist);
#else
#ifdef PM_ZOOM_OPTIMIZED
  pmforce_zoom_optimized_prepare_density(mode, typelist);
#else
  pmforce_uniform_optimized_prepare_density(mode);
#endif
#endif

  /* allocate the memory to hold the FFT fields */

  forcegrid = (fft_real *) mymalloc("forcegrid", maxfftsize * sizeof(fft_real));

  workspace = forcegrid;

#ifndef FFT_COLUMN_BASED
  fft_of_rhogrid = (fft_complex *) & rhogrid[0];
#else
  fft_of_rhogrid = (fftw_complex *) & workspace[0];
#endif

  /* Do the FFT of the density field */
#ifndef FFT_COLUMN_BASED
  my_slab_based_fft(&myplan, &rhogrid[0], &workspace[0], 1);
#else
  my_column_based_fft(&myplan, rhogrid, workspace, 1);  /* result is in workspace, not in rhogrid ! */
#endif


  if(mode != 0)
    {
      pmforce_measure_powerspec(mode - 1, typelist);
    }
  else
    {
      /* multiply with Green's function in order to obtain the potential (or forces for spectral diffencing) */

      double kfacx = 2.0 * M_PI / (GRIDX * d);
      double kfacy = 2.0 * M_PI / (GRIDY * d);
      double kfacz = 2.0 * M_PI / (GRIDZ * d);

#ifdef FFT_COLUMN_BASED
#pragma omp parallel for private(x, y, z, xx, yy, zz)
      for(large_array_offset ip = 0; ip < myplan.second_transposed_ncells; ip++)
        {
          large_array_offset ipcell = ip + ((large_array_offset) myplan.second_transposed_firstcol) * GRIDX;
          y = ipcell / (GRIDX * GRIDz);
          int yr = ipcell % (GRIDX * GRIDz);
          z = yr / GRIDX;
          x = yr % GRIDX;
#else
#pragma omp parallel for private(y, z)
      for(x = 0; x < GRIDX; x++)
        for(y = myplan.slabstart_y; y < myplan.slabstart_y + myplan.nslab_y; y++)
          for(z = 0; z < GRIDz; z++)
            {
#endif
              if(x >= (GRIDX / 2))
                xx = x - GRIDX;
              else
                xx = x;
              if(y >= (GRIDY / 2))
                yy = y - GRIDY;
              else
                yy = y;
              if(z >= (GRIDZ / 2))
                zz = z - GRIDZ;
              else
                zz = z;

              double kx = kfacx * xx;
              double ky = kfacy * yy;
              double kz = kfacz * zz;

              double k2 = kx * kx + ky * ky + kz * kz;

              if(k2 > 0)
                {
                  double smth = -exp(-k2 * asmth2) / k2;

                  /* do deconvolution */

                  double fx = 1, fy = 1, fz = 1;

                  if(xx != 0)
                    {
                      fx = kx * dhalf;
                      fx = sin(fx) / fx;
                    }
                  if(yy != 0)
                    {
                      fy = ky * dhalf;
                      fy = sin(fy) / fy;
                    }
                  if(zz != 0)
                    {
                      fz = kz * dhalf;
                      fz = sin(fz) / fz;
                    }

                  double ff = 1 / (fx * fy * fz);
                  double deconv = ff * ff * ff * ff;

                  smth *= deconv;       /* deconvolution */

#ifndef FFT_COLUMN_BASED
                  large_array_offset ip = ((large_array_offset) GRIDz) * (GRIDX * (y - myplan.slabstart_y) + x) + z;
#endif

#ifdef GRAVITY_TALLBOX
                  double re = fft_of_rhogrid[ip][0] * fft_of_kernel[ip][0] - fft_of_rhogrid[ip][1] * fft_of_kernel[ip][1];
                  double im = fft_of_rhogrid[ip][0] * fft_of_kernel[ip][1] + fft_of_rhogrid[ip][1] * fft_of_kernel[ip][0];

                  fft_of_rhogrid[ip][0] = re * deconv * exp(-k2 * asmth2);
                  fft_of_rhogrid[ip][1] = im * deconv * exp(-k2 * asmth2);
#else
                  fft_of_rhogrid[ip][0] *= smth;
                  fft_of_rhogrid[ip][1] *= smth;
#endif
                }
            }

#ifdef FFT_COLUMN_BASED
      if(myplan.second_transposed_firstcol == 0)
        fft_of_rhogrid[0][0] = fft_of_rhogrid[0][1] = 0.0;
#else
      if(myplan.slabstart_y == 0)
        fft_of_rhogrid[0][0] = fft_of_rhogrid[0][1] = 0.0;
#endif

      /* Do the inverse FFT to get the potential/forces */

#ifndef FFT_COLUMN_BASED
      my_slab_based_fft(&myplan, &rhogrid[0], &workspace[0], -1);
#else
      my_column_based_fft(&myplan, workspace, rhogrid, -1);
#endif

      /* Now rhogrid holds the potential/forces */


#ifdef EVALPOTENTIAL
#ifdef PM_ZOOM_OPTIMIZED
      pmforce_zoom_optimized_readout_forces_or_potential(-1);
#else
      pmforce_uniform_optimized_readout_forces_or_potential(-1);
#endif
#endif

      /* get the force components by finite differencing of the potential for each dimension,
       * and send the results back to the right CPUs
       */
      for(int dim = 2; dim >= 0; dim--) /* Calculate each component of the force. */
        {
          /* we do the x component last, because for differencing the potential in the x-direction, we need to construct the transpose */

#ifndef FFT_COLUMN_BASED
          if(dim == 0)
            {
              my_slab_transposeA(&myplan, rhogrid, forcegrid);  /* compute the transpose of the potential field for finite differencing */
              /* note: for the x-direction, we difference the transposed field */

#pragma omp parallel for private(y, z)
              for(x = 0; x < GRIDX; x++)
                for(y = 0; y < myplan.nslab_y; y++)
                  for(z = 0; z < GRIDZ; z++)
                    {
                      int xrr = x + 2, xll = x - 2, xr = x + 1, xl = x - 1;
                      if(xr >= GRIDX)
                        xr -= GRIDX;
                      if(xrr >= GRIDX)
                        xrr -= GRIDX;
                      if(xl < 0)
                        xl += GRIDX;
                      if(xll < 0)
                        xll += GRIDX;

                      forcegrid[NI(x, y, z)] = fac * ((4.0 / 3) * (rhogrid[NI(xl, y, z)] - rhogrid[NI(xr, y, z)]) - (1.0 / 6) * (rhogrid[NI(xll, y, z)] - rhogrid[NI(xrr, y, z)]));
                    }

              my_slab_transposeB(&myplan, forcegrid, rhogrid);  /* reverse the transpose from above */
            }
          else
            {
#pragma omp parallel for private(x, z)
              for(y = 0; y < GRIDY; y++)
                for(x = 0; x < myplan.nslab_x; x++)
                  for(z = 0; z < GRIDZ; z++)
                    {
                      if(dim == 1)
                        {
                          int yr = y + 1, yl = y - 1, yrr = y + 2, yll = y - 2;
                          if(yr >= GRIDY)
                            yr -= GRIDY;
                          if(yrr >= GRIDY)
                            yrr -= GRIDY;
                          if(yl < 0)
                            yl += GRIDY;
                          if(yll < 0)
                            yll += GRIDY;

                          forcegrid[FI(x, y, z)] = fac * ((4.0 / 3) * (rhogrid[FI(x, yl, z)] - rhogrid[FI(x, yr, z)]) - (1.0 / 6) * (rhogrid[FI(x, yll, z)] - rhogrid[FI(x, yrr, z)]));
                        }
                      else if(dim == 2)
                        {
                          int zr = z + 1, zl = z - 1, zrr = z + 2, zll = z - 2;
                          if(zr >= GRIDZ)
                            zr -= GRIDZ;
                          if(zrr >= GRIDZ)
                            zrr -= GRIDZ;
                          if(zl < 0)
                            zl += GRIDZ;
                          if(zll < 0)
                            zll += GRIDZ;

                          forcegrid[FI(x, y, z)] = fac * ((4.0 / 3) * (rhogrid[FI(x, y, zl)] - rhogrid[FI(x, y, zr)]) - (1.0 / 6) * (rhogrid[FI(x, y, zll)] - rhogrid[FI(x, y, zrr)]));
                        }
                    }
            }

#else

          if(dim == 2)
            {
#pragma omp parallel for
              for(large_array_offset i = 0; i < myplan.base_ncol; i++)
                {
                  fft_real *forcep = &forcegrid[GRID2 * i];
                  fft_real *potp = &rhogrid[GRID2 * i];

                  for(int z = 0; z < GRIDZ; z++)
                    {
                      int zr = z + 1;
                      int zl = z - 1;
                      int zrr = z + 2;
                      int zll = z - 2;

                      if(zr >= GRIDZ)
                        zr -= GRIDZ;
                      if(zrr >= GRIDZ)
                        zrr -= GRIDZ;
                      if(zl < 0)
                        zl += GRIDZ;
                      if(zll < 0)
                        zll += GRIDZ;

                      forcep[z] = fac * ((4.0 / 3) * (potp[zl] - potp[zr]) - (1.0 / 6) * (potp[zll] - potp[zrr]));
                    }
                }
            }
          else if(dim == 1)
            {
              fft_real *scratch = mymalloc("scratch", myplan.fftsize * sizeof(fft_real));       /* need a third field as scratch space */
              memcpy(scratch, rhogrid, myplan.fftsize * sizeof(fft_real));

              my_fft_swap23(&myplan, scratch, forcegrid);

#pragma omp parallel for private(forcep, potp)
              for(large_array_offset i = 0; i < myplan.ncol_XZ; i++)
                {
                  fft_real *forcep = &scratch[GRIDY * i];
                  fft_real *potp = &forcegrid[GRIDY * i];

                  for(int y = 0; y < GRIDY; y++)
                    {
                      int yr = y + 1;
                      int yl = y - 1;
                      int yrr = y + 2;
                      int yll = y - 2;

                      if(yr >= GRIDY)
                        yr -= GRIDY;
                      if(yrr >= GRIDY)
                        yrr -= GRIDY;
                      if(yl < 0)
                        yl += GRIDY;
                      if(yll < 0)
                        yll += GRIDY;

                      forcep[y] = fac * ((4.0 / 3) * (potp[yl] - potp[yr]) - (1.0 / 6) * (potp[yll] - potp[yrr]));
                    }
                }

              my_fft_swap23back(&myplan, scratch, forcegrid);
              myfree(scratch);
            }
          else if(dim == 0)
            {
              fft_real *scratch = mymalloc("scratch", myplan.fftsize * sizeof(fft_real));       /* need a third field as scratch space */
              memcpy(scratch, rhogrid, myplan.fftsize * sizeof(fft_real));

              my_fft_swap13(&myplan, scratch, forcegrid);

#pragma omp parallel for private(forcep, potp)
              for(large_array_offset i = 0; i < myplan.ncol_YZ; i++)
                {
                  fft_real *forcep = &scratch[GRIDX * i];
                  fft_real *potp = &forcegrid[GRIDX * i];

                  for(int x = 0; x < GRIDX; x++)
                    {
                      int xr = x + 1;
                      int xl = x - 1;
                      int xrr = x + 2;
                      int xll = x - 2;

                      if(xr >= GRIDX)
                        xr -= GRIDX;
                      if(xrr >= GRIDX)
                        xrr -= GRIDX;
                      if(xl < 0)
                        xl += GRIDX;
                      if(xll < 0)
                        xll += GRIDX;

                      forcep[x] = fac * ((4.0 / 3) * (potp[xl] - potp[xr]) - (1.0 / 6) * (potp[xll] - potp[xrr]));
                    }
                }

              my_fft_swap13back(&myplan, scratch, forcegrid);
              myfree(scratch);
            }
#endif

#ifdef PM_ZOOM_OPTIMIZED
          pmforce_zoom_optimized_readout_forces_or_potential(dim);
#else
          pmforce_uniform_optimized_readout_forces_or_potential(dim);
#endif
        }

    }

  /* free stuff */

  myfree(forcegrid);
  myfree(rhogrid);

#ifdef PM_ZOOM_OPTIMIZED
  myfree(localfield_recvcount);
  myfree(localfield_offset);
  myfree(localfield_sendcount);
  myfree(localfield_first);
  myfree(localfield_data);
  myfree(localfield_globalindex);
  myfree(part);
#else
#ifndef POWERSPECTRUM_IN_POSTPROCESSING
  myfree(partin);
#endif
  myfree(Rcvpm_offset);
  myfree(Rcvpm_count);
  myfree(Sndpm_offset);
  myfree(Sndpm_count);
#endif

  double tend = second();

  if(mode == 0)
    mpi_printf("PM-PERIODIC: done.  (took %g seconds)\n", timediff(tstart, tend));
}


#ifdef GRAVITY_TALLBOX

/*! This function sets-up the Greens function for calculating the tall-box potential
 *  in real space, with suitable zero padding in the direction of the tall box.
 */
static void pmforce_setup_tallbox_kernel(void)
{
  int i, j, k, ii, jj, kk;
  double xx, yy, zz;

  double d = All.BoxSize / PMGRID;

  mpi_printf("PM-PERIODIC: Setting up tallbox kernel (GRIDX=%d, GRIDY=%d, GRIDZ=%d)\n", GRIDX, GRIDY, GRIDZ);

  /* now set up kernel and its Fourier transform */

  for(i = 0; i < maxfftsize; i++)       /* clear local field */
    kernel[i] = 0;

#ifndef FFT_COLUMN_BASED
  for(i = myplan.slabstart_x; i < (myplan.slabstart_x + myplan.nslab_x); i++)
    for(j = 0; j < GRIDY; j++)
      {
#else
  for(int c = myplan.base_firstcol; c < (myplan.base_firstcol + myplan.base_ncol); c++)
    {
      i = c / GRIDY;
      j = c % GRIDY;
#endif

      for(k = 0; k < GRIDZ; k++)
        {
          if(i >= (GRIDX / 2))
            ii = i - GRIDX;
          else
            ii = i;
          if(j >= (GRIDY / 2))
            jj = j - GRIDY;
          else
            jj = j;
          if(k >= (GRIDZ / 2))
            kk = k - GRIDZ;
          else
            kk = k;

          xx = ii * d;
          yy = jj * d;
          zz = kk * d;

          double pot = pmperiodic_tallbox_long_range_potential(xx, yy, zz);

#ifndef FFT_COLUMN_BASED
          size_t ip = FI(i - myplan.slabstart_x, j, k);
#else
          size_t ip = FC(c, k);
#endif
          kernel[ip] = pot / All.BoxSize;
        }

#ifndef FFT_COLUMN_BASED
    }
#else
    }
#endif

  /* Do the FFT of the kernel */

  fft_real *workspc = (fft_real *) mymalloc("workspc", maxfftsize * sizeof(fft_real));

#ifndef FFT_COLUMN_BASED
  my_slab_based_fft(&myplan, kernel, workspc, 1);
#else
  my_column_based_fft(&myplan, kernel, workspc, 1);     /* result is in workspace, not in kernel */
  memcpy(kernel, workspc, maxfftsize * sizeof(fft_real));
#endif

  myfree(workspc);

  mpi_printf("PM-PERIODIC: Done setting up tallbox kernel\n");
}

static double pmperiodic_tallbox_long_range_potential(double x, double y, double z)
{
  x /= All.BoxSize;
  y /= All.BoxSize;
  z /= All.BoxSize;

  double r = sqrt(x * x + y * y + z * z);

  if(r == 0)
    return 0;

  double xx, yy, zz;
  switch (GRAVITY_TALLBOX)
    {
    case 0:
      xx = y;
      yy = z;
      zz = x;
      break;
    case 1:
      xx = x;
      yy = z;
      zz = y;
      break;
    case 2:
      xx = x;
      yy = y;
      zz = z;
      break;
    }
  x = xx;
  y = yy;
  z = zz;

  /* the third dimension, z, is now the non-periodic one */

  double lmin = imin(BOXX, BOXY);

  double alpha = 2.0 / lmin;

  double sum1 = 0.0;

  const int nmax = 4;

  for(int nx = -nmax; nx <= nmax; nx++)
    for(int ny = -nmax; ny <= nmax; ny++)
      {
        double dx = x - nx * BOXX;
        double dy = y - ny * BOXY;
        double r = sqrt(dx * dx + dy * dy + z * z);
        if(r > 0)
          sum1 += erfc(alpha * r) / r;
      }

  int nxmax = (int) (2 * alpha * (BOXX / lmin) + 0.5);
  int nymax = (int) (2 * alpha * (BOXY / lmin) + 0.5);

  double alpha2 = alpha * alpha;

  double sum2 = 0.0;

  for(int nx = -nxmax; nx <= nxmax; nx++)
    for(int ny = -nymax; ny <= nymax; ny++)
      {
        if(nx != 0 || ny != 0)
          {
            double kx = (2.0 * M_PI / BOXX) * nx;
            double ky = (2.0 * M_PI / BOXY) * ny;
            double k2 = kx * kx + ky * ky;
            double k = sqrt(k2);

            sum2 += cos(kx * x + ky * y) * (exp(k * z) * erfc(k / (2 * alpha) + alpha * z) + exp(-k * z) * erfc(k / (2 * alpha) - alpha * z)) / k;
          }
      }

  sum2 *= M_PI / (BOXX * BOXY);

  double psi = 2.0 * alpha / sqrt(M_PI) + (2 * sqrt(M_PI) / (BOXX * BOXY) * (exp(-alpha2 * z * z) / alpha + sqrt(M_PI) * z * erf(alpha * z))) - (sum1 + sum2);

  return psi;
}
#endif


#ifdef PM_ZOOM_OPTIMIZED

/*! \brief Sort function for #part array indices.
 *
 * Sorts the indices into the #part array by the global index of the corresponding #part_slab_data struct.
 *
 * \param a index to be compared
 * \param b index to be compared
 * \return sort result
 */
static int pm_periodic_compare_sortindex(const void *a, const void *b)
{
  if(part[*(int *) a].globalindex < part[*(int *) b].globalindex)
    return -1;

  if(part[*(int *) a].globalindex > part[*(int *) b].globalindex)
    return +1;

  return 0;
}

/*! \brief Implements the sorting function for mysort_pmperiodic()
 *
 * The index array is sorted using a merge sort algorithm.
 *
 * \param b index array to sort
 * \param n number of elements to sort
 * \param t temporary buffer array
 */
static void msort_pmperiodic_with_tmp(large_numpart_type * b, size_t n, large_numpart_type * t)
{
  large_numpart_type *tmp;
  large_numpart_type *b1, *b2;
  size_t n1, n2;

  if(n <= 1)
    return;

  n1 = n / 2;
  n2 = n - n1;
  b1 = b;
  b2 = b + n1;

  msort_pmperiodic_with_tmp(b1, n1, t);
  msort_pmperiodic_with_tmp(b2, n2, t);

  tmp = t;

  while(n1 > 0 && n2 > 0)
    {
      if(part[*b1].globalindex <= part[*b2].globalindex)
        {
          --n1;
          *tmp++ = *b1++;
        }
      else
        {
          --n2;
          *tmp++ = *b2++;
        }
    }

  if(n1 > 0)
    memcpy(tmp, b1, n1 * sizeof(large_numpart_type));

  memcpy(b, t, (n - n2) * sizeof(large_numpart_type));
}

/*! \brief Sort the index array b of n entries using the sort kernel
 * cmp.
 *
 * The parameter s is set to sizeof(int). The index array b
 * is sorted according to the globalindex field of the referenced item in the #part array
 *
 * \param b the index array to sort
 * \param n number of entries in array b
 * \param size of each entry (must be sizeof(int))
 * \param cmp comparator function
 */
static void mysort_pmperiodic(void *b, size_t n, size_t s, int (*cmp) (const void *, const void *))
{
  const size_t size = n * s;

  large_numpart_type *tmp = (large_numpart_type *) mymalloc("tmp", size);

  msort_pmperiodic_with_tmp((large_numpart_type *) b, n, tmp);

  myfree(tmp);
}
#endif


/*----------------------------------------------------------------------------------------------------*/
/*           Here comes code for the power-sepctrum computation                                       */
/*----------------------------------------------------------------------------------------------------*/


static char power_spec_fname[MAXLEN_PATH];

void calculate_power_spectra(int num, long long *ntot_type_all)
{
  int i, typeflag[NTYPES];

  power_spec_totnumpart = 0;

  for(i = 0; i < NTYPES; i++)
    {
      typeflag[i] = 1;
      power_spec_totnumpart += ntot_type_all[i];
    }

  sprintf(power_spec_fname, "%s/powerspec_%03d.txt", All.OutputDir, num);

  pmforce_do_powerspec(typeflag);       /* calculate power spectrum for all particle types */

  if(ntot_type_all)
    for(i = 0; i < NTYPES; i++)
      {
        if(ntot_type_all[i] > 0)
          {
            int j;
            for(j = 0; j < NTYPES; j++)
              typeflag[j] = 0;

            typeflag[i] = 1;
            power_spec_totnumpart = ntot_type_all[i];

            sprintf(power_spec_fname, "%s/powerspec_type%d_%03d.txt", All.OutputDir, i, num);
            
            pmforce_do_powerspec(typeflag);     /* calculate power spectrum for type i */
          }
      }
}


void pmforce_do_powerspec(int *typeflag)
{
  mpi_printf("POWERSPEC: Begin power spectrum. (typeflag=[%d|%d|%d|%d|%d|%d])\n", typeflag[0], typeflag[1], typeflag[2], typeflag[3], typeflag[4], typeflag[5]);
  
  MPI_Barrier( MPI_COMM_WORLD );

  double tstart = second();

  pmforce_periodic(1, typeflag);        /* calculate regular power spectrum for selected particle types */

  pmforce_periodic(2, typeflag);        /* calculate folded power spectrum for selected particle types  */
  
  pmforce_periodic(3, typeflag);        /* calculate double folded power spectrum for selected particle types  */

  double tend = second();

  mpi_printf("POWERSPEC: End power spectrum. took %g seconds\n", timediff(tstart, tend));
}



void pmforce_measure_powerspec(int flag, int *typeflag)
{
  long long CountModes[BINS_PS];
  double SumPowerUncorrected[BINS_PS];       /* without binning correction (as for shot noise) */
  double PowerUncorrected[BINS_PS];  /* without binning correction */
  double DeltaUncorrected[BINS_PS];  /* without binning correction */
  double ShotLimit[BINS_PS];
  double Kbin[BINS_PS];

#ifndef POWERSPECTRUM_IN_POSTPROCESSING
  double mass = 0;
  double mass2 = 0;
  for(int i = 0; i < NumPart; i++)
    if(typeflag[P[i].Type])
      {
        double m = P[i].Mass;
#ifdef BLACK_HOLES
        if(P[i].Type == 5)
          m = BPP(i).BH_Mass;
#endif
        mass += m;
        mass2 += m * m;
      }

  MPI_Allreduce(&mass, &power_spec_totmass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&mass2, &power_spec_totmass2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  double d = All.BoxSize / PMGRID;
  double dhalf = 0.5 * d;

  double fac = 1.0 / power_spec_totmass;

  double K0 = 2 * M_PI / All.BoxSize;   /* minimum k */
  //double K1 = K0 * All.BoxSize / get_default_softening_of_particletype(1);      /* maximum k */
  double K1 = K0 * (POWERSPEC_FOLDFAC * POWERSPEC_FOLDFAC * PMGRID / 2);      /* maximum k */
  double binfac = BINS_PS / (log(K1) - log(K0));
  
  mpi_printf( "K0=%g, K1=%g, binfac=%g, totmass=%g\n", K0, K1, binfac, power_spec_totmass );

  double kfacx = 2.0 * M_PI / (STRETCHX * All.BoxSize);
  double kfacy = 2.0 * M_PI / (STRETCHY * All.BoxSize);
  double kfacz = 2.0 * M_PI / (STRETCHZ * All.BoxSize);

  for(int i = 0; i < BINS_PS; i++)
    {
      SumPowerUncorrected[i] = 0;
      CountModes[i] = 0;
    }
  
#ifdef FFT_COLUMN_BASED
  for(large_array_offset ip = 0; ip < myplan.second_transposed_ncells; ip++)
    {
      large_array_offset ipcell = ip + ((large_array_offset) myplan.second_transposed_firstcol) * GRIDX;
      int y = ipcell / (GRIDX * GRIDz);
      int yr = ipcell % (GRIDX * GRIDz);
      int z = yr / GRIDX;
      int x = yr % GRIDX;
#else
  for(int y = myplan.slabstart_y; y < myplan.slabstart_y + myplan.nslab_y; y++)
    for(int x = 0; x < GRIDX; x++)
      for(int z = 0; z < GRIDz; z++)
        {
#endif
          int count_double;

          if(z >= 1 && z < (GRIDZ + 1) / 2)     /* these modes need to be counted twice due to the storage scheme for the FFT of a real field */
            count_double = 1;
          else
            count_double = 0;
          
          int xx, yy, zz;

          if(x >= (GRIDX / 2))
            xx = x - GRIDX;
          else
            xx = x;

          if(y >= (GRIDY / 2))
            yy = y - GRIDY;
          else
            yy = y;

          if(z >= (GRIDZ / 2))
            zz = z - GRIDZ;
          else
            zz = z;

          double kx = kfacx * xx;
          double ky = kfacy * yy;
          double kz = kfacz * zz;

          double k2 = kx * kx + ky * ky + kz * kz;

          if(k2 > 0)
            {
              /* do deconvolution */

              double fx = 1, fy = 1, fz = 1;
              if(xx != 0)
                {
                  fx = kx * dhalf;
                  fx = sin(fx) / fx;
                }
              if(yy != 0)
                {
                  fy = ky * dhalf;
                  fy = sin(fy) / fy;
                }
              if(zz != 0)
                {
                  fz = kz * dhalf;
                  fz = sin(fz) / fz;
                }
              double ff = 1 / (fx * fy * fz);
              double smth = ff * ff * ff * ff;
              /*
               * Note: The Fourier-transform of the density field (rho_hat) must be multiplied with ff^2
               * in order to do the de-convolution. Thats why po = rho_hat^2 gains a factor of ff^4. (Christian Arnold)
               */

              /* end deconvolution */

#ifndef FFT_COLUMN_BASED
              large_array_offset ip = ((large_array_offset) GRIDz) * (GRIDX * (y - myplan.slabstart_y) + x) + z;
#endif

              double po = (fft_of_rhogrid[ip][0] * fft_of_rhogrid[ip][0] + fft_of_rhogrid[ip][1] * fft_of_rhogrid[ip][1]);

              po *= fac * fac * smth;

              double k = sqrt(k2);

              if(flag == 1)
                k *= POWERSPEC_FOLDFAC;
              if(flag == 2)
                k *= POWERSPEC_FOLDFAC * POWERSPEC_FOLDFAC;

              if(k >= K0 && k < K1)
                {
                  int bin = log(k / K0) * binfac;

                  SumPowerUncorrected[bin] += po;
                  CountModes[bin] += 1;

                  if(count_double)
                    {
                      SumPowerUncorrected[bin] += po;
                      CountModes[bin] += 1;
                    }
                }
            }
        }

  /* Now compute the power spectrum */

  long long *countbuf = (long long int *) mymalloc("countbuf", NTask * BINS_PS * sizeof(long long));
  double *powerbuf = (double *) mymalloc("powerbuf", NTask * BINS_PS * sizeof(double));

  MPI_Allgather(CountModes, BINS_PS * sizeof(long long), MPI_BYTE, countbuf, BINS_PS * sizeof(long long), MPI_BYTE, MPI_COMM_WORLD);

  for(int i = 0; i < BINS_PS; i++)
    {
      CountModes[i] = 0;
      for(int n = 0; n < NTask; n++)
        CountModes[i] += countbuf[n * BINS_PS + i];
    }

  MPI_Allgather(SumPowerUncorrected, BINS_PS * sizeof(double), MPI_BYTE, powerbuf, BINS_PS * sizeof(double), MPI_BYTE, MPI_COMM_WORLD);

  for(int i = 0; i < BINS_PS; i++)
    {
      SumPowerUncorrected[i] = 0;
      for(int n = 0; n < NTask; n++)
        SumPowerUncorrected[i] += powerbuf[n * BINS_PS + i];
    }

  myfree(powerbuf);
  myfree(countbuf);

  for(int i = 0; i < BINS_PS; i++)
    {
      Kbin[i] = exp((i + 0.5) / binfac + log(K0));

      if(CountModes[i] > 0)
        PowerUncorrected[i] = SumPowerUncorrected[i] / CountModes[i];
      else
        PowerUncorrected[i] = 0;

      DeltaUncorrected[i] = 4 * M_PI * pow(Kbin[i], 3) / pow(2 * M_PI / All.BoxSize, 3) * PowerUncorrected[i];

      ShotLimit[i] = 4 * M_PI * pow(Kbin[i], 3) / pow(2 * M_PI / All.BoxSize, 3) * (power_spec_totmass2 / (power_spec_totmass*power_spec_totmass));
    }

  if(ThisTask == 0)
    {
      FILE *fd;
      
      if(flag == 0)
        {
          if(!(fd = fopen(power_spec_fname, "w"))) /* store the unfolded spectrum */
            terminate("can't open file `%s`\n", power_spec_fname);
        }
      else if(flag == 1 || flag == 2)
        {
          if(!(fd = fopen(power_spec_fname, "a"))) /* append the file, store the folded spectrum */
            terminate("can't open file `%s`\n", power_spec_fname);
        }
      else
        terminate("Something wrong.\n");

      fprintf(fd, "%g\n", All.Time);
      fprintf(fd, "%d\n", (int) BINS_PS);
      fprintf(fd, "%g\n", power_spec_totmass);
      fprintf(fd, "%15lld\n", (long long) (power_spec_totnumpart));

      for(int i = 0; i < BINS_PS; i++)
        fprintf(fd, "%g %g %g %g %g\n", Kbin[i], DeltaUncorrected[i], PowerUncorrected[i], (double) CountModes[i], ShotLimit[i]);
      
      fclose(fd);
    }
}

#endif
