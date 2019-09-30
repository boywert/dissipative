/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/voronoi_exchange.c
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

#include "../allvars.h"
#include "../proto.h"
#include "../voronoi.h"
#include "RT.h"
#include "RT_proto.h"

#ifdef VORONOI


void exchange_primitive_variables_RT(void)
{
  if(All.TotNumGas == 0)
    return;

  //  TIMER_START(CPU_MESH_EXCHANGE);

  int listp;
  struct rt_primexch *tmpPrimExch;
  int i, j, p, task, off;
  int ngrp, recvTask, place;

  tmpPrimExch = (struct rt_primexch *) mymalloc("rt_tmpPrimExch", Mesh_nexport * sizeof(struct rt_primexch));

  /* prepare data for export */
  for(j = 0; j < NTask; j++)
    Mesh_Send_count[j] = 0;

  for(i = 0; i < NumGasInMesh; i++)
    {
      p = List_InMesh[i];

      listp = List_P[p].firstexport;
      while(listp >= 0)
        {
          if((task = ListExports[listp].origin) != ThisTask)
            {
              place = ListExports[listp].index;
              off = Mesh_Send_offset[task] + Mesh_Send_count[task]++;

	      for(int num3=0;num3<MRT_BINS;num3++)
		{
		  tmpPrimExch[off].OldCons_DensPhot[num3] = SphP[place].OldCons_DensPhot[num3] ;
		  if(TimeBinSynchronized[P[place].TimeBinHydro])
		    {
		      tmpPrimExch[off].DensPhot[num3] = SphP[place].DensPhot[num3] ;
		      tmpPrimExch[off].modFN[num3] = SphP[place].modFN[num3] ;
		      for(int num1=0;num1<3;num1++)
			{
			  tmpPrimExch[off].RT_F[num3][num1] = SphP[place].RT_F[num3][num1] ;
			  tmpPrimExch[off].FN[num3][num1] = SphP[place].FN[num3][num1] ;
			  for(int num2=0;num2<3;num2++)
			    tmpPrimExch[off].PT[num3][num1][num2] = SphP[place].PT[num3][num1][num2] ;
			}
		    }
		  else
		    {
		      tmpPrimExch[off].DensPhot[num3] = SphP[place].Cons_DensPhot[num3]/SphP[place].Volume ;
		      double modFN ;

		      modFN = sqrt(SphP[place].Cons_RT_F[num3][0]*SphP[place].Cons_RT_F[num3][0] + SphP[place].Cons_RT_F[num3][1]*SphP[place].Cons_RT_F[num3][1] + SphP[place].Cons_RT_F[num3][2]*SphP[place].Cons_RT_F[num3][2])/SphP[place].Cons_DensPhot[num3] ;

                      tmpPrimExch[off].modFN[num3] = modFN ;
                      for(int num1=0;num1<3;num1++)
                        {
			  double FN ;
			  FN = SphP[i].Cons_RT_F[num3][num1] / SphP[i].Cons_DensPhot[num3] ;
                          tmpPrimExch[off].RT_F[num3][num1] =  SphP[place].Cons_RT_F[num3][num1]/SphP[place].Volume ;
                          tmpPrimExch[off].FN[num3][num1] = FN ;
                          for(int num2=0;num2<3;num2++)
                            tmpPrimExch[off].PT[num3][num1][num2] = SphP[place].PT[num3][num1][num2] ;
			}
		    }
		}
            }
          listp = ListExports[listp].nextexport;
        }
    }


  /* exchange data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Mesh_Send_count[recvTask] > 0 || Mesh_Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&tmpPrimExch[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask]
                           * sizeof(struct rt_primexch), MPI_BYTE, recvTask, TAG_DENS_A,
                           &RTPrimExch[Mesh_Recv_offset[recvTask]], Mesh_Recv_count[recvTask] * sizeof(struct rt_primexch), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  myfree(tmpPrimExch);

  //  TIMER_STOP(CPU_MESH_EXCHANGE);
}

void exchange_primitive_variables_and_gradients_RT(void)
{
  if(All.TotNumGas == 0)
    return;

  //TIMER_START(CPU_MESH_EXCHANGE);

  int listp;
  struct rt_grad_data *tmpGradExch;
  struct rt_primexch *tmpPrimExch;

  int i, j, p, task, off;
  int ngrp, recvTask, place;

  tmpPrimExch = (struct rt_primexch *) mymalloc("rt_tmpPrimExch", Mesh_nexport * sizeof(struct rt_primexch));
  tmpGradExch = (struct rt_grad_data *) mymalloc("rt_tmpGradExch", Mesh_nexport * sizeof(struct rt_grad_data));


  /* prepare data for export */
  for(j = 0; j < NTask; j++)
    Mesh_Send_count[j] = 0;

  for(i = 0; i < NumGasInMesh; i++)
    {
      p = List_InMesh[i];

      /* in case previous steps already lowered the Mass, update OldMass to yield together with metallicity vector conservative estimate of metal mass of each species contained in cell */
#ifdef MRT_LSF_GRADIENTS
      for(int num1=0; num1<MRT_BINS; num1++)
	{
	  if(SphP[p].Cons_DensPhot[num1] < SphP[p].OldCons_DensPhot[num1])
	    SphP[p].OldCons_DensPhot[num1] = SphP[p].Cons_DensPhot[num1] ;
	}
#endif


      listp = List_P[p].firstexport;
      while(listp >= 0)
        {
          if((task = ListExports[listp].origin) != ThisTask)
            {
              place = ListExports[listp].index;
              off = Mesh_Send_offset[task] + Mesh_Send_count[task]++;

              int num1, num2, num3 ;
	      for(num3=0;num3<MRT_BINS;num3++)
		{
		  tmpPrimExch[off].OldCons_DensPhot[num3] = SphP[place].OldCons_DensPhot[num3] ;
		  if(TimeBinSynchronized[P[place].TimeBinHydro])
		    {
		      tmpPrimExch[off].DensPhot[num3] = SphP[place].DensPhot[num3] ;
		      tmpPrimExch[off].modFN[num3] = SphP[place].modFN[num3] ;
		      for(int num1=0;num1<3;num1++)
			{
			  tmpPrimExch[off].RT_F[num3][num1] = SphP[place].RT_F[num3][num1] ;
			  tmpPrimExch[off].FN[num3][num1] = SphP[place].FN[num3][num1] ;
			  for(int num2=0;num2<3;num2++)
			    tmpPrimExch[off].PT[num3][num1][num2] = SphP[place].PT[num3][num1][num2] ;
			}
		    }
		  else
		    {
		      tmpPrimExch[off].DensPhot[num3] = SphP[place].Cons_DensPhot[num3]/SphP[place].Volume ;
		      double modFN ;

		      modFN = sqrt(SphP[place].Cons_RT_F[num3][0]*SphP[place].Cons_RT_F[num3][0] + SphP[place].Cons_RT_F[num3][1]*SphP[place].Cons_RT_F[num3][1] + SphP[place].Cons_RT_F[num3][2]*SphP[place].Cons_RT_F[num3][2])/SphP[place].Cons_DensPhot[num3] ;

                      tmpPrimExch[off].modFN[num3] = modFN ;
                      for(int num1=0;num1<3;num1++)
                        {
			  double FN ;
			  FN = SphP[i].Cons_RT_F[num3][num1] / SphP[i].Cons_DensPhot[num3] ;
                          tmpPrimExch[off].RT_F[num3][num1] =  SphP[place].Cons_RT_F[num3][num1]/SphP[place].Volume ;
                          tmpPrimExch[off].FN[num3][num1] = FN ;
                          for(int num2=0;num2<3;num2++)
                            tmpPrimExch[off].PT[num3][num1][num2] = SphP[place].PT[num3][num1][num2] ;
			}
		    }
		}
              tmpGradExch[off] = SphP[place].RTGrad;
            }
          listp = ListExports[listp].nextexport;
        }
    }


  /* exchange data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Mesh_Send_count[recvTask] > 0 || Mesh_Recv_count[recvTask] > 0)
            {
              /* exchange the data */
              MPI_Sendrecv(&tmpPrimExch[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask]
                           * sizeof(struct rt_primexch), MPI_BYTE, recvTask, TAG_DENS_A,
                           &RTPrimExch[Mesh_Recv_offset[recvTask]], Mesh_Recv_count[recvTask] * sizeof(struct rt_primexch), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

              MPI_Sendrecv(&tmpGradExch[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask]
                           * sizeof(struct rt_grad_data), MPI_BYTE, recvTask, TAG_HYDRO_A,
                           &RTGradExch[Mesh_Recv_offset[recvTask]], Mesh_Recv_count[recvTask] * sizeof(struct rt_grad_data), MPI_BYTE, recvTask, TAG_HYDRO_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            }
        }
    }

  myfree(tmpGradExch);
  myfree(tmpPrimExch);

  //  TIMER_STOP(CPU_MESH_EXCHANGE);

  /* note: because the sequence is the same as before, we don't have to do the sorts again */
}



#endif /* VORONOI */
