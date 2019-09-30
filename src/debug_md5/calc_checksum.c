/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/debug_md5/calc_checksum.c
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

#include "../allvars.h"
#include "../proto.h"


#include "Md5.h"

void calc_memory_checksum(void *base, size_t bytes)
{
  MD5_CTX sum;
  union
  {
    unsigned char digest[16];
    int val[4];
  } u, uglob;


  MD5Init(&sum);
  MD5UpdateLong(&sum, base, bytes);
  MD5Final(&sum);

  int i;

  for(i = 0; i < 16; i++)
    u.digest[i] = sum.digest[i];

  MPI_Allreduce(u.val, uglob.val, 4, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      printf("Step=%d  MD5=", All.NumCurrentTiStep);
      for(i = 0; i < 16; i++)
        printf("%02x", uglob.digest[i]);
      printf("\n");
    }
}


#ifdef RESTART_DEBUG

void log_restart_debug(void)
{
  MD5_CTX sum;
  union
  {
    unsigned char digest[16];
    int val[4];
  } u, uglob_P, uglob_SphP, uglob_StarP, uglob_BHP;
  int i;

  MD5Init(&sum);
  MD5UpdateLong(&sum, (void *) P, NumPart * sizeof(struct particle_data));
  MD5Final(&sum);

  for(i = 0; i < 16; i++)
    u.digest[i] = sum.digest[i];

  MPI_Allreduce(u.val, uglob_P.val, 4, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  MD5Init(&sum);
  MD5UpdateLong(&sum, (void *) SphP, NumGas * sizeof(struct sph_particle_data));
  MD5Final(&sum);

  for(i = 0; i < 16; i++)
    u.digest[i] = sum.digest[i];

  MPI_Allreduce(u.val, uglob_SphP.val, 4, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

#ifdef USE_SFR
  MD5Init(&sum);
  MD5UpdateLong(&sum, (void *) SphP, N_star * sizeof(struct star_particle_data));
  MD5Final(&sum);

  for(i = 0; i < 16; i++)
    u.digest[i] = sum.digest[i];

  MPI_Allreduce(u.val, uglob_StarP.val, 4, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

#ifdef BLACK_HOLES
  MD5Init(&sum);
  MD5UpdateLong(&sum, (void *) BHP, NumBHs * sizeof(struct bh_particle_data));
  MD5Final(&sum);

  for(i = 0; i < 16; i++)
    u.digest[i] = sum.digest[i];

  MPI_Allreduce(u.val, uglob_BHP.val, 4, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

  if(ThisTask == 0)
    {
      fprintf(FdRestartTest, "\n");
      fprintf(FdRestartTest, "Step=%8d  P[]        ", All.NumCurrentTiStep);
      for(i = 0; i < 16; i++)
        fprintf(FdRestartTest, "%02x", uglob_P.digest[i]);
      fprintf(FdRestartTest, "\n");
      fprintf(FdRestartTest, "               SphP[]     ");
      for(i = 0; i < 16; i++)
        fprintf(FdRestartTest, "%02x", uglob_SphP.digest[i]);
      fprintf(FdRestartTest, "\n");
#ifdef USE_SFR
      fprintf(FdRestartTest, "               StarP[]    ");
      for(i = 0; i < 16; i++)
        fprintf(FdRestartTest, "%02x", uglob_StarP.digest[i]);
      fprintf(FdRestartTest, "\n");
#endif
#ifdef BLACK_HOLES
      fprintf(FdRestartTest, "               BHP[]      ");
      for(i = 0; i < 16; i++)
        fprintf(FdRestartTest, "%02x", uglob_BHP.digest[i]);
      fprintf(FdRestartTest, "\n");
#endif
      fflush(FdRestartTest);
    }
}
#endif

