/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/cuda_util.h
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

#ifndef CUDA_UTIL_H
#define CUDA_UTIL_H

#include <cuda.h>
#include <cuda_runtime.h>

#define CudaSafeCall( err )     __cudaSafeCall( err, __FILE__, __LINE__ )
#define CudaCheckError()        __cudaCheckError( __FILE__, __LINE__ )

#define CuDrvSafeCall( err )     __cuDrvSafeCall( err, __FILE__, __LINE__ )

void __cudaCheckError(const char *file, const int line);
void __cudaSafeCall(cudaError_t err, const char *file, const int line);
void __cuDrvSafeCall(CUresult err, const char *file, const int line);


#endif
