/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/dtypes.h
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

#ifndef DTYPES_H
#define DTYPES_H



#ifndef FFTW
#define CONCAT(prefix, name) prefix ## name
#ifdef DOUBLEPRECISION_FFTW
#define FFTW(x) CONCAT(fftw_, x)
#else
#define FFTW(x) CONCAT(fftwf_, x)
#endif
#endif



#ifndef SPECIAL_BOUNDARY
#ifndef LONGIDS
typedef unsigned int MyIDType;
#define MPI_MYIDTYPE MPI_UNSIGNED
#else
typedef unsigned long long MyIDType;
#define MPI_MYIDTYPE MPI_UNSIGNED_LONG_LONG
#endif
#else
#ifndef LONGIDS
typedef int MyIDType;
#define MPI_MYIDTYPE MPI_INT
#else
typedef long long MyIDType;
#define MPI_MYIDTYPE MPI_LONG_LONG_INT
#endif
#endif

#ifndef DOUBLEPRECISION         /* default is single-precision */
typedef float MySingle;
typedef float MyFloat;
typedef float MyDouble;
#define MPI_MYFLOAT MPI_FLOAT
#define MPI_MYDOUBLE MPI_FLOAT
#else
#if (DOUBLEPRECISION == 2)      /* mixed precision */
typedef float MySingle;
typedef float MyFloat;
typedef double MyDouble;
#define MPI_MYFLOAT MPI_FLOAT
#define MPI_MYDOUBLE MPI_DOUBLE
#else
#if (DOUBLEPRECISION == 3)      /* mixed precision, fewer single precision variables */
typedef float  MySingle;
typedef double MyFloat;
typedef double MyDouble;
#define MPI_MYFLOAT MPI_FLOAT
#define MPI_MYDOUBLE MPI_DOUBLE
#else
/* everything double-precision */
typedef double MySingle;
typedef double MyFloat;
typedef double MyDouble;
#define MPI_MYFLOAT MPI_DOUBLE
#define MPI_MYDOUBLE MPI_DOUBLE
#endif
#endif
#endif

#ifdef OUTPUT_IN_DOUBLEPRECISION
typedef double MyOutputFloat;
#else
typedef float MyOutputFloat;
#endif

#ifdef INPUT_IN_DOUBLEPRECISION
typedef double MyInputFloat;
#else
typedef float MyInputFloat;
#endif


#ifndef NGB_TREE_DOUBLEPRECISION
typedef float MyNgbTreeFloat;
#define MAX_NGBRANGE_NUMBER MAX_FLOAT_NUMBER
#else
typedef double MyNgbTreeFloat;
#define MAX_NGBRANGE_NUMBER MAX_DOUBLE_NUMBER
#endif


#if defined(PMGRID) || defined(POWERSPEC_GRID)

#include <fftw3.h>

#ifdef DOUBLEPRECISION_FFTW
typedef double fft_real;
typedef fftw_complex fft_complex;
#else
typedef float fft_real;
typedef fftwf_complex fft_complex;
#endif
typedef ptrdiff_t fft_ptrdiff_t;


typedef struct
{
  int NgridX, NgridY, NgridZ;
  int Ngridz, Ngrid2;

  FFTW(plan) forward_plan_zdir;
  FFTW(plan) forward_plan_xdir;
  FFTW(plan) forward_plan_ydir;

  FFTW(plan) backward_plan_zdir;
  FFTW(plan) backward_plan_ydir;
  FFTW(plan) backward_plan_xdir;

#ifndef FFT_COLUMN_BASED

  int *slab_to_task;            /*!< Maps a slab index to the task responsible for the slab */
  int *slabs_x_per_task;
  int *first_slab_x_of_task;      /*!< Array containing the index of the first slab of each task */
  int *slabs_y_per_task;        /*!< Array containing the number of slabs each task is responsible for */
  int *first_slab_y_of_task;    /*!< Array containing the index of the first slab of each task */

  int nslab_x, slabstart_x, nslab_y, slabstart_y;
  int largest_x_slab;             /*!< size of the largest slab in x direction */
  int largest_y_slab;             /*!< size of the largest slab in y direction */

#else
  size_t max_datasize;
  size_t fftsize;


  int base_firstcol, base_ncol, base_lastcol;
  int transposed_firstcol, transposed_ncol;
  int second_transposed_firstcol, second_transposed_ncol;
  size_t second_transposed_ncells;

  int firstcol_XZ, ncol_XZ;
  int firstcol_YZ, ncol_YZ;

  int pivotcol;                 /* to go from column number to task */
  int avg;
  int tasklastsection;

  size_t *offsets_send_A;
  size_t *offsets_recv_A;
  size_t *offsets_send_B;
  size_t *offsets_recv_B;
  size_t *offsets_send_C;
  size_t *offsets_recv_C;
  size_t *offsets_send_D;
  size_t *offsets_recv_D;
  size_t *offsets_send_13;
  size_t *offsets_recv_13;
  size_t *offsets_send_23;
  size_t *offsets_recv_23;
  size_t *offsets_send_13back;
  size_t *offsets_recv_13back;
  size_t *offsets_send_23back;
  size_t *offsets_recv_23back;

  size_t *count_send_A;
  size_t *count_recv_A;
  size_t *count_send_B;
  size_t *count_recv_B;
  size_t *count_send_C;
  size_t *count_recv_C;
  size_t *count_send_D;
  size_t *count_recv_D;
  size_t *count_send_13;
  size_t *count_recv_13;
  size_t *count_send_23;
  size_t *count_recv_23;
  size_t *count_send_13back;
  size_t *count_recv_13back;
  size_t *count_send_23back;
  size_t *count_recv_23back;
#endif
}
fft_plan;

#endif

#endif
