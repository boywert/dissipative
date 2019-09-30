/*! * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/anisotropic_RT/RT_semi_HYPRE.c
 * \date        MM/YYYY
 * \author      Rahul Kannan
 * \brief        
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 */


/*! \file RT_semi_HYPRE.c
 *  \brief main driver for an moment based RT with the VET formalism
 *
 *  
 */



#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#include "../allvars.h"
#include "../proto.h"
#include "../domain.h"
#include "../voronoi.h"



#ifdef MRT_SETUP_SPECIAL_BOUNDARIES

void mrt_setup(void) 
{
  int i ;
  for(i=0;i<NumGas;i++)
    {

      /*      if(P[i].ID >= 40000000 && P[i].ID <= 49999999)
	{
	  SphP[i].RT_F[0][0] = 1e-57 * pow(All.UnitLength_in_cm, 2) * All.UnitTime_in_s ;
	  SphP[i].DensPhot[0] = 1.01 * SphP[i].RT_F[0][0]/c_internal_units ;
	  SphP[i].RT_F[0][1] = 0.0 ;
	  SphP[i].RT_F[0][2] = 0.0 ;

	  SphP[i].Cons_DensPhot[0] = SphP[i].DensPhot[0]*SphP[i].Volume ;
	  SphP[i].Cons_RT_F[0][0] = SphP[i].RT_F[0][0]*SphP[i].Volume ;
	  SphP[i].Cons_RT_F[0][1] = SphP[i].RT_F[0][1]*SphP[i].Volume ;
	  SphP[i].Cons_RT_F[0][2] = SphP[i].RT_F[0][2]*SphP[i].Volume ;
	  }*/

      if((P[i].ID >= 60000000 && P[i].ID <= 69999999))
	{
	  SphP[i].RT_F[0][0] = MINDENSPHOT ;
	  SphP[i].DensPhot[0] = MINDENSPHOT ;	  
	  SphP[i].RT_F[0][1] = MINDENSPHOT ;
	  SphP[i].RT_F[0][2] = MINDENSPHOT ;

	  SphP[i].Cons_DensPhot[0] = SphP[i].DensPhot[0]*SphP[i].Volume ;
	  SphP[i].Cons_RT_F[0][0] = SphP[i].RT_F[0][0]*SphP[i].Volume ;
	  SphP[i].Cons_RT_F[0][1] = SphP[i].RT_F[0][1]*SphP[i].Volume ;
	  SphP[i].Cons_RT_F[0][2] = SphP[i].RT_F[0][2]*SphP[i].Volume ;

	  SphP[i].Density = 1e-13*0.0028555628 ;
	  SphP[i].Utherm = 0.43574726*1e-3 ;
	}



      /*      if((P[i].ID >= 40000000 && P[i].ID <= 49999999))
	{
	  for(int num1=0;num1<MRT_BINS;num1++)
	    {
	      SphP[i].DensPhot[num1] = 1.03e4 * pow(All.UnitLength_in_cm, 2) * All.UnitTime_in_s / All.UnitEnergy_in_cgs / c_internal_units; 
	      // }SphP[i].Cons_DensPhot[num1]
	      
	      SphP[i].Cons_DensPhot[num1] = SphP[i].DensPhot[num1] * SphP[i].Volume ;
	      if(isnan(SphP[i].Cons_DensPhot[num1]))
		terminate("Something went wrong\n") ;
	      
	      for(int j=0;j<3;j++)
		{
		  if(j==1)
		    SphP[i].RT_F[num1][j] = 0.99 * c_internal_units * SphP[i].DensPhot[num1] ;
		  else
		    SphP[i].RT_F[num1][j] = MINDENSPHOT;
		  
		  SphP[i].Cons_RT_F[num1][j] = SphP[i].RT_F[num1][j] * SphP[i].Volume ;
		}
	      
	    }
	}
	  /*	  SphP[i].RT_F[0][0] = MINDENSPHOT ;
	  SphP[i].DensPhot[0] = MINDENSPHOT ;	  
	  SphP[i].RT_F[0][1] = MINDENSPHOT ;
	  SphP[i].RT_F[0][2] = MINDENSPHOT ;

	  SphP[i].Cons_DensPhot[0] = SphP[i].DensPhot[0]*SphP[i].Volume ;
	  SphP[i].Cons_RT_F[0][0] = SphP[i].RT_F[0][0]*SphP[i].Volume ;
	  SphP[i].Cons_RT_F[0][1] = SphP[i].RT_F[0][1]*SphP[i].Volume ;
	  SphP[i].Cons_RT_F[0][2] = SphP[i].RT_F[0][2]*SphP[i].Volume ;

	  SphP[i].Density = 1e-13*0.0028555628 ;
	  SphP[i].Utherm = 0.0010152911 ;
	  }*/


      /*  if(P[i].ID >= 40000000 && P[i].ID <= 59999999)
	{
	  //double F = 1.03e4 * pow(All.UnitLength_in_cm, 2) * All.UnitTime_in_s / All.UnitEnergy_in_cgs ;
	  //SphP[i].Cons_DensPhot[0] += dt * F * (All.BoxSize/258.0) ;
	  // SphP[i].Cons_RT_F[0][0] = MINDENSPHOT ;
	  //SphP[i].Cons_DensPhot[0] = MINDENSPHOT ;

	  /*	  double F = 1.03e4 * pow(All.UnitLength_in_cm, 2) * All.UnitTime_in_s / All.UnitEnergy_in_cgs ;


	  SphP[i].RT_F[0][1] = 0.99 * F;
	  SphP[i].Cons_RT_F[0][1] = SphP[i].RT_F[0][1] * SphP[i].Volume;*/

      /*	  SphP[i].RT_F[0][0] = MINDENSPHOT ;
	  SphP[i].DensPhot[0] = MINDENSPHOT ;	  
	  SphP[i].RT_F[0][1] = MINDENSPHOT ;
	  SphP[i].RT_F[0][2] = MINDENSPHOT ;

	  SphP[i].Cons_DensPhot[0] = SphP[i].DensPhot[0]*SphP[i].Volume ;
	  SphP[i].Cons_RT_F[0][0] = SphP[i].RT_F[0][0]*SphP[i].Volume ;
	  SphP[i].Cons_RT_F[0][1] = SphP[i].RT_F[0][1]*SphP[i].Volume ;
	  SphP[i].Cons_RT_F[0][2] = SphP[i].RT_F[0][2]*SphP[i].Volume ;


	  
	  SphP[i].Density = 10.0 * 0.0028555628 ;
	  SphP[i].Utherm = 1.0152911 ;
	  }*/
    //	}
  //    }
    }
}

/*static struct inexch
{
  MyIDType ID ;
  double DensPhot[1] ;
  double RT_F[1][3] ;
}
  *InExch ;

void mrt_setup(void) 
{
  InExch = (struct inexch *) mymalloc("InExch", Mesh_nimport * sizeof(struct inexch));

  exchange_vector() ;

  mpi_printf("RT: Setting up special boundary conditions\n") ;
  point *DP = Mesh.DP;
  face *VF = Mesh.VF;
  tetra *DT = Mesh.DT;

  double n[3], mdls ;
  int dp, vf, particle ;

  double Fstar = 1.03e4 * pow(All.UnitLength_in_cm, 2) * All.UnitTime_in_s / All.UnitEnergy_in_cgs ; 

  MyIDType id_K, id_L ;
  double E_K, E_L, Fx_K, Fy_K, Fz_K, Fx_L, Fy_L, Fz_L ;

  int num1 = IR_BINS ;

  for(int i = 0; i < NumGas; i++)
    { 
      if(P[i].ID >= 40000000 && P[i].ID <= 49999999)
	{

	  SphP[i].Density = 1e-5*0.0028555628 ;                                                                                           
          SphP[i].Utherm = 1.0152911 ; 

	  id_K = P[i].ID ;
	  int num_fl = 0 ;
	  double totnum = 0.0 ;
	  E_K = SphP[i].DensPhot[num1] ;
	  Fx_K = SphP[i].RT_F[num1][0] ;
	  Fy_K = SphP[i].RT_F[num1][1] ;
	  Fz_K = SphP[i].RT_F[num1][2] ;


	  int q = SphP[i].first_connection;
	  while(q >= 0)
	    {
	      dp = DC[q].dp_index;
	      vf = DC[q].vf_index;
	      particle = DP[dp].index;
	      
	      
	      if(particle < 0) 
		{
		  q = DC[q].next;
		  continue;
		}
	      
	      if(DP[dp].task == ThisTask)
		{
		  if(particle >= NumGas)
		    particle -= NumGas;

		  id_L = P[particle].ID ;
		  E_L = SphP[particle].DensPhot[num1] ;
		  Fx_L = SphP[particle].RT_F[num1][0] ;
		  Fy_L = SphP[particle].RT_F[num1][1] ;
		  Fz_L = SphP[particle].RT_F[num1][2] ;
		}
	      else
		{
		  id_L = InExch[particle].ID ;
                  E_L = InExch[particle].DensPhot[num1] ;
                  Fx_L = InExch[particle].RT_F[num1][0] ;
                  Fy_L = InExch[particle].RT_F[num1][1] ;
                  Fz_L = InExch[particle].RT_F[num1][2] ;
		}


	      if(id_L == id_K+10000000)
		{
		  //		  printf("Entered Here\n") ;
		  SphP[i].DensPhot[0] = (Fstar - Fy_L + c_internal_units*E_L)/c_internal_units ;
		  SphP[i].RT_F[0][1] = 0.99 * SphP[i].DensPhot[0] *c_internal_units ;
		  SphP[i].RT_F[0][0] = 0.0 ;
		  SphP[i].RT_F[0][2] = 0.0 ;
		}


	      if(q == SphP[i].last_connection)
		break;
	      
	      q = DC[q].next;
	    }
	}
      if(P[i].ID >= 30000000 && P[i].ID <= 39999999)
	{
	  SphP[i].RT_F[0][0] = MINDENSPHOT ;
	  SphP[i].DensPhot[0] = MINDENSPHOT ;	  
	  SphP[i].RT_F[0][1] = MINDENSPHOT ;
	  SphP[i].RT_F[0][2] = MINDENSPHOT ;

	  SphP[i].Cons_DensPhot[0] = SphP[i].DensPhot[0]*SphP[i].Volume ;
	  SphP[i].Cons_RT_F[0][0] = SphP[i].RT_F[0][0]*SphP[i].Volume ;
	  SphP[i].Cons_RT_F[0][1] = SphP[i].RT_F[0][1]*SphP[i].Volume ;
	  SphP[i].Cons_RT_F[0][2] = SphP[i].RT_F[0][2]*SphP[i].Volume ;

	  SphP[i].Density = 1e-13*0.0028555628 ;
	  SphP[i].Utherm = 0.0010152911 ;
	}
    }
  myfree(InExch) ;
  //  terminate("\n") ;
}

void exchange_vector()
{
  int listp;
  int j, p, task, off;
  int ngrp, recvTask, place;    // there was earlier the following defined: sendTask 


  struct inexch *tmpInExch;
  tmpInExch = (struct inexch *) mymalloc("tmpInExch", Mesh_nexport * sizeof(struct inexch));

  //   prepare data for export 
  for(j = 0; j < NTask; j++)
    Mesh_Send_count[j] = 0;

  for(p = 0; p < NumGas; p++)
    {
      if(P[p].Type == 0)
        {
          listp = List_P[p].firstexport;

          while(listp >= 0)
            {
              if((task = ListExports[listp].origin) != ThisTask)
                {
                  place = ListExports[listp].index;
                  off = Mesh_Send_offset[task] + Mesh_Send_count[task]++;

                  tmpInExch[off].ID = P[place].ID ;
                  tmpInExch[off].DensPhot[0] = SphP[place].DensPhot[0] ;
                  tmpInExch[off].RT_F[0][0] = SphP[place].RT_F[0][0] ;
                  tmpInExch[off].RT_F[0][1] = SphP[place].RT_F[0][1] ;
                  tmpInExch[off].RT_F[0][2] = SphP[place].RT_F[0][2] ; ;
                }
              listp = ListExports[listp].nextexport;
            }
        }
    }

  // exchange data 
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Mesh_Send_count[recvTask] > 0 || Mesh_Recv_count[recvTask] > 0)
            {
              // get the particles 
              MPI_Sendrecv(&tmpInExch[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask] *
                           sizeof(struct inexch), MPI_BYTE, recvTask, TAG_DENS_A, &InExch[Mesh_Recv_offset[recvTask]],
                           Mesh_Recv_count[recvTask] * sizeof(struct inexch), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  myfree(tmpInExch);
}
*/

#endif
