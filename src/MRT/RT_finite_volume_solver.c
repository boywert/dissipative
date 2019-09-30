/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/MRT/RT_finite_volume_solver.c
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../allvars.h"
#include "../proto.h"
#include "../voronoi.h"
#include "RT.h"
#include "RT_proto.h"

#ifdef MRT

static struct flux_list_data
{
  int task, index;
 
  double dDensPhot[MRT_BINS] ;
  double dRT_F[MRT_BINS][3] ;  
}
 *FluxList;

static int Nflux, MaxNflux;


struct primexch *PrimExch;
struct grad_data *GradExch;

struct rt_primexch *RTPrimExch ;
struct rt_grad_data *RTGradExch ;

/*!< state on a face determined by riemann solver */
struct state_face state_face;

/*!< flux through a face */
struct fluxes fluxes;
struct geometry geom;


void compute_interface_fluxes_RT(tessellation * T)
{

  mpi_printf("RT: computing interface fluxes....\n") ;
  int num1, num2, num3 ;

 
  int i, j;
  double count = 0, count_reduced = 0, tot_count, tot_count_reduced;
  double face_dt, hubble_a, atime;


  MaxNflux = T->Indi.AllocFacNflux;
  Nflux = 0;
  FluxList = mymalloc_movable(&FluxList, "FluxList", MaxNflux * sizeof(struct flux_list_data));

  face *VF = T->VF;
  point *DP = T->DP;


  for(i = 0; i < T->Nvf; i++)
    {
      struct state state_L, state_center_L, delta_time_L, delta_space_L;
      struct state state_R, state_center_R, delta_time_R, delta_space_R;

      face_dt = 0;              /* the default is that this face is not active */

      /* calculate normal vectors */
      if(face_get_normals(T, i, &geom))
        continue;

      /* get the values of the states at the center of the cells */
      if(face_get_state_RT(T, VF[i].p1, i, &state_center_L))
        continue;

      if(face_get_state_RT(T, VF[i].p2, i, &state_center_R))
        continue;

      /* only treat faces where one of the two sides is active */
      if(!TimeBinSynchronized[state_center_L.timeBin] && !TimeBinSynchronized[state_center_R.timeBin])
        continue;

      /* clarify whether the face should be done by this task (it may be present also on another task) */
      if(face_check_responsibility_of_this_task(T, VF[i].p1, VF[i].p2, &state_center_L, &state_center_R))
        continue;

      /* calculate timestep of the face */
      face_dt = face_timestep_RT(&state_center_L, &state_center_R, &hubble_a, &atime); 

      //face_dt /= ((double) All.RTNumSubCycles) ;

      if(!(face_dt > 0))
        continue;

      /* now estimate the velocity of the midpoint of the face based on the velocities of the generators of the mesh. */
      double vel_face[3];

      if(All.ComovingIntegrationOn)
        for(j = 0; j < 3; j++)
          {
            state_center_L.velVertex[j] /= atime;       /* convert vertex motion to peculiar velocity */
            state_center_R.velVertex[j] /= atime;
          }

      /* rough motion of mid-point of edge */
      vel_face[0] = 0.5 * (state_center_L.velVertex[0] + state_center_R.velVertex[0]);
      vel_face[1] = 0.5 * (state_center_L.velVertex[1] + state_center_R.velVertex[1]);
      vel_face[2] = 0.5 * (state_center_L.velVertex[2] + state_center_R.velVertex[2]);

      double cx, cy, cz, facv;

      cx = VF[i].cx - 0.5 * (DP[VF[i].p2].x + DP[VF[i].p1].x);
      cy = VF[i].cy - 0.5 * (DP[VF[i].p2].y + DP[VF[i].p1].y);
      cz = VF[i].cz - 0.5 * (DP[VF[i].p2].z + DP[VF[i].p1].z);

      facv = (cx * (state_center_L.velVertex[0] - state_center_R.velVertex[0]) +
              cy * (state_center_L.velVertex[1] - state_center_R.velVertex[1]) + cz * (state_center_L.velVertex[2] - state_center_R.velVertex[2])) / geom.nn;

      /* put in a limiter for highly distorted cells */
      double cc = sqrt(cx * cx + cy * cy + cz * cz);
      if(cc > 0.9 * geom.nn)
        facv *= (0.9 * geom.nn) / cc;

      vel_face[0] += facv * geom.nx;
      vel_face[1] += facv * geom.ny;
      vel_face[2] += facv * geom.nz;

#if defined(VORONOI_STATIC_MESH) || defined(AMR)
      vel_face[0] = 0;
      vel_face[1] = 0;
      vel_face[2] = 0;
#endif

      double vel_face_turned[3];
      
      /* turn the face velocity */
      vel_face_turned[0] = vel_face[0] * geom.nx + vel_face[1] * geom.ny + vel_face[2] * geom.nz;
      vel_face_turned[1] = vel_face[0] * geom.mx + vel_face[1] * geom.my + vel_face[2] * geom.mz;
      vel_face_turned[2] = vel_face[0] * geom.px + vel_face[1] * geom.py + vel_face[2] * geom.pz;

      state_convert_to_local_frame_RT(&state_center_L, vel_face, hubble_a, atime);
      state_convert_to_local_frame_RT(&state_center_R, vel_face, hubble_a, atime);



      /* copy center state to state at interface, then add extrapolation terms */
      state_L = state_center_L;
      state_R = state_center_R;


#ifdef MRT_TIME_EXTRAPOLATION
      calculate_div_F(&delta_time_L, &state_center_L) ;
      calculate_div_F(&delta_time_R, &state_center_R) ;
#endif
     

      face_do_time_extrapolation_RT(&delta_time_L, &state_center_L, atime);
      face_do_time_extrapolation_RT(&delta_time_R, &state_center_R, atime);

      face_do_spatial_extrapolation_RT(&delta_space_L, &state_center_L, &state_center_R);
      face_do_spatial_extrapolation_RT(&delta_space_R, &state_center_R, &state_center_L);

      face_add_extrapolations_RT(&state_L, &delta_time_L, &delta_space_L);
      face_add_extrapolations_RT(&state_R, &delta_time_R, &delta_space_R);

      
#ifdef MRT
#if defined(BOUNDARY_REFL_FLUIDSIDE_MINID) && defined(BOUNDARY_REFL_FLUIDSIDE_MAXID) && defined(BOUNDARY_REFL_SOLIDSIDE_MINID) && defined(BOUNDARY_REFL_SOLIDSIDE_MAXID)
      MyIDType idL = state_L.ID;
      MyIDType idR = state_R.ID;

      int nn1 ;
      double F = 1e-57 * pow(All.UnitLength_in_cm,2) * All.UnitTime_in_s ;

      if(idR >= BOUNDARY_REFL_FLUIDSIDE_MINID && idR < BOUNDARY_REFL_FLUIDSIDE_MAXID && idL >= BOUNDARY_REFL_SOLIDSIDE_MINID && idL < BOUNDARY_REFL_SOLIDSIDE_MAXID)
	{
	  for(nn1=0;nn1<MRT_BINS;nn1++)
	    {
	      state_L.DensPhot[nn1] = 1.01*(F - state_R.RT_F[nn1][0] + c_internal_units*state_R.DensPhot[nn1])/c_internal_units ;
	      state_L.RT_F[nn1][0] = F ;
	      state_L.RT_F[nn1][1] = 0.0;
	      state_L.RT_F[nn1][2] = 0.0;
	      state_L.OldCons_DensPhot[nn1] = 1.01*(F - state_R.RT_F[nn1][0] + c_internal_units*state_R.DensPhot[nn1])/c_internal_units ;
	    }
	}
      else if(idL >= BOUNDARY_REFL_FLUIDSIDE_MINID && idL < BOUNDARY_REFL_FLUIDSIDE_MAXID && idR >= BOUNDARY_REFL_SOLIDSIDE_MINID && idR < BOUNDARY_REFL_SOLIDSIDE_MAXID)
	{
	  for(nn1=0;nn1<MRT_BINS;nn1++)
	    {
	      state_R.DensPhot[nn1] = 1.01*(F - state_L.RT_F[nn1][0] + c_internal_units*state_L.DensPhot[nn1])/c_internal_units ;
	      state_R.RT_F[nn1][0] = F ;
	      state_R.RT_F[nn1][1] = 0.0;
	      state_R.RT_F[nn1][2] = 0.0;
	      state_L.OldCons_DensPhot[nn1] = 1.01*(F - state_R.RT_F[nn1][0] + c_internal_units*state_R.DensPhot[nn1])/c_internal_units ;
	    }
	}
#endif
#endif
      
      /* turn the velocities to get velx perpendicular and vely and velz in the plane of the face */
      face_turn_velocities_RT(&state_L, &geom);
      face_turn_velocities_RT(&state_R, &geom);



      ///#ifndef MESHRELAX

#ifdef MRT_RIEMANN_ROSUNOV
      double press2 = godunov_flux_3d_rosunov_RT(&state_L, &state_R, &state_face, &fluxes);
#else
      double press2 = godunov_flux_3d_HLLE_RT(&state_L, &state_R, &state_face, &fluxes);
#endif

      //#endif /* end of MESHRELAX */


      /* turn the velocity field back */
      //      face_turnback_velocities_RT(&state_face, &geom);

      /* add the face velocity again */
      state_face.velx += vel_face[0];
      state_face.vely += vel_face[1];
      state_face.velz += vel_face[2];

      //#ifndef MESHRELAX

#ifndef MRT_COMOVING
      flux_convert_to_lab_frame_RT(&state_L, &state_R, vel_face_turned, &fluxes);
#endif
      face_turn_momentum_flux_RT(&fluxes, &geom);


      face_limit_fluxes_RT(&state_L, &state_R, &state_center_L, &state_center_R, &fluxes, face_dt, &count, &count_reduced);

      /*#else /* MESHRELAX */

      /* just solve the advection equation instead of Riemann problem */

      /*solve_advection(&state_L, &state_R, &state_face, &geom, vel_face);
      face_clear_fluxes(&fluxes);
      face_add_fluxes_advection(&state_face, &fluxes, &geom, vel_face);
      face_set_scalar_states_and_fluxes(&state_L, &state_R, &state_face, &fluxes);

      #endif*/


      /* now apply the flux to update the conserved states of the cells */


      if(face_dt > 0)           /* selects active faces */
        {
          int k, p, q;
          double dir;
          double fac = face_dt * VF[i].area;

#ifndef MUSCL_HANCOCK
#ifndef RUNGE_KUTTA_FULL_UPDATE   
          fac *= 0.5;
#else
          fac *= 1.;      
#endif          
#endif

          for(k = 0; k < 2; k++)
            {
              if(k == 0)
                {
                  q = VF[i].p1;
                  p = DP[q].index;
                  dir = -fac;
#if defined(TRACER_MC) || defined(MHD_CT) || defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
                  qother = VF[i].p2;
#endif
                }
              else
                {
                  q = VF[i].p2;
                  p = DP[q].index;
                  dir = +fac;
#if defined(TRACER_MC) || defined(MHD_CT)
                  qother = VF[i].p1;
#endif
                }

	      // This piece of code is designed to ensure only outflow at inflow/outflow boundaries in problems where this is desired
#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
	      if(k == 0)
		{
		  if(DP[q].index >= NumGas || DP[qother].index >= NumGas) /* ghost point */
		    {
#if defined(REFLECTIVE_Z) && defined(EXTERNALSHEARBOX)
		      /* check if we are at z boundary */
		      if(((DP[q].image_flags & OUTFLOW_Z) && (DP[q].image_flags & REFL_Z_FLAGS)) ||
			 ((DP[qother].image_flags & OUTFLOW_Z) && (DP[qother].image_flags & REFL_Z_FLAGS)))
			{
			  double vec;
			  if(REFLECTIVE_Z == 2){
			    if(DP[q].index >= NumGas) /* ghost point */
			      // p is a ghost cell!  If mass is being removed from it, this is inflow
			      if(dir*fluxes.mass < 0.)
				break;
			    
			    if(DP[qother].index >= NumGas) /* ghost point */
			      //pother is a ghost cell!  If p is gaining mass then we have inflow!
			      if(dir*fluxes.mass > 0.)
				break;
			  }
			}
#endif
		    }

		}
#endif

	      

              if(DP[q].task == ThisTask)
                {
                  if(DP[q].index >= NumGas)     /* this is a local ghost point */
                    {
                      if(DP[VF[i].p1].ID == DP[VF[i].p2].ID)    /* this may happen for reflective points */
                        continue;
                      p -= NumGas;
                    }

		  for(num1=0;num1<MRT_BINS;num1++)
		    {
		      SphP[p].Cons_DensPhot[num1] += dir * fluxes.DensPhot[num1] ;
		      SphP[p].Cons_RT_F[num1][0] += dir * fluxes.RT_F[num1][0] ;
		      SphP[p].Cons_RT_F[num1][1] += dir * fluxes.RT_F[num1][1] ;
		      SphP[p].Cons_RT_F[num1][2] += dir * fluxes.RT_F[num1][2] ;
		    }

                }
              else
                {
                  /* here we have a foreign ghost point */
                  if(DP[q].originalindex < 0)
                    terminate("should not happen");

                  if(Nflux >= MaxNflux)
                    {
                      T->Indi.AllocFacNflux *= ALLOC_INCREASE_FACTOR;
                      MaxNflux = T->Indi.AllocFacNflux;

		      FluxList = myrealloc_movable(FluxList, MaxNflux * sizeof(struct flux_list_data));

                      if(Nflux >= MaxNflux)
                        terminate("Nflux >= MaxNflux");
                    }

                  FluxList[Nflux].task = DP[q].task;
                  FluxList[Nflux].index = DP[q].originalindex;

		  for(num1=0;num1<MRT_BINS;num1++)
		    {
		      FluxList[Nflux].dDensPhot[num1] = dir * fluxes.DensPhot[num1] ;
		      FluxList[Nflux].dRT_F[num1][0] = dir * fluxes.RT_F[num1][0] ;
		      FluxList[Nflux].dRT_F[num1][1] = dir * fluxes.RT_F[num1][1] ;
		      FluxList[Nflux].dRT_F[num1][2] = dir * fluxes.RT_F[num1][2] ;
		    }
                  Nflux++;
                }
            }
        }
    }
  /* end of big loop over all faces */

  //  TIMER_STOPSTART(CPU_FLUXES, CPU_FLUXES_COMM);

  /* now exchange the flux-list and apply it when needed */
  apply_flux_list_RT();

  //TIMER_STOPSTART(CPU_FLUXES_COMM, CPU_FLUXES);

  myfree(FluxList);

  double in[2] = {count, count_reduced},  out[2];
  MPI_Reduce(in, out, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(ThisTask == 0)
    {
      tot_count = out[0];
      tot_count_reduced = out[1];

      printf("RT FLUX: exchanged fluxes over %g faces, with %g reduced (fraction %g), cumulative fraction %g\n",
             tot_count, tot_count_reduced, tot_count_reduced / (tot_count + 1.0e-30), All.TotCountReducedFluxes / (All.TotCountFluxes + 1.0e-30));
      All.TotCountReducedFluxes += tot_count_reduced;
      All.TotCountFluxes += tot_count;
    }

  //  TIMER_STOP(CPU_RT_FLUXES);
}


int face_get_state_RT(tessellation * T, int p, int i, struct state *st)
{
  int particle;

  double aBegin;

  point *DP = T->DP;
  face *VF = T->VF;

  particle = DP[p].index;

  if(particle < 0)
    return -1;

  if(particle >= NumGas && DP[p].task == ThisTask)
    particle -= NumGas;

  /* interpolation vector for the left state */
  if(DP[p].task == ThisTask)
    {
      st->dx = VF[i].cx - SphP[particle].Center[0];
      st->dy = VF[i].cy - SphP[particle].Center[1];
      st->dz = VF[i].cz - SphP[particle].Center[2];
    }
  else
    {
      st->dx = VF[i].cx - PrimExch[particle].Center[0];
      st->dy = VF[i].cy - PrimExch[particle].Center[1];
      st->dz = VF[i].cz - PrimExch[particle].Center[2];
    }

  /* correct for periodicity */
#ifdef PERIODIC
#if !defined(REFLECTIVE_X) && !defined(ONEDIMS_SPHERICAL)
  if(st->dx < -boxHalf_X)
    st->dx += boxSize_X;
  if(st->dx > boxHalf_X)
    st->dx -= boxSize_X;
#endif
#if !defined(REFLECTIVE_Y)
  if(st->dy < -boxHalf_Y)
    st->dy += boxSize_Y;
  if(st->dy > boxHalf_Y)
    st->dy -= boxSize_Y;
#endif
#if !defined(REFLECTIVE_Z)
  if(st->dz < -boxHalf_Z)
    st->dz += boxSize_Z;
  if(st->dz > boxHalf_Z)
    st->dz -= boxSize_Z;
#endif
#endif


  if(DP[p].task == ThisTask)
    {
      st->velGas[0] = P[particle].Vel[0];
      st->velGas[1] = P[particle].Vel[1];
      st->velGas[2] = P[particle].Vel[2];

      st->velVertex[0] = SphP[particle].VelVertex[0];
      st->velVertex[1] = SphP[particle].VelVertex[1];
      st->velVertex[2] = SphP[particle].VelVertex[2];

      st->rho = SphP[particle].Density ;
      for(int num1=0;num1<MRT_BINS;num1++)
	{ 
	  if(TimeBinSynchronized[P[particle].TimeBinHydro])
	    {
	      st->DensPhot[num1] = SphP[particle].DensPhot[num1] ;
	      st->OldCons_DensPhot[num1] = SphP[particle].OldCons_DensPhot[num1] ;
	      st->RT_F[num1][0] = SphP[particle].RT_F[num1][0] ;
	      st->RT_F[num1][1] = SphP[particle].RT_F[num1][1] ;
	      st->RT_F[num1][2] = SphP[particle].RT_F[num1][2] ;
	      st->modFN[num1] = SphP[particle].modFN[num1] ;
	    }
	  else
	    {
	      st->DensPhot[num1] = SphP[particle].Cons_DensPhot[num1]/SphP[particle].Volume ;

	      double modFN = sqrt(SphP[particle].Cons_RT_F[num1][0]*SphP[particle].Cons_RT_F[num1][0] + SphP[particle].Cons_RT_F[num1][1]*SphP[particle].Cons_RT_F[num1][1] + SphP[particle].Cons_RT_F[num1][2]*SphP[particle].Cons_RT_F[num1][2])/SphP[particle].Cons_DensPhot[num1] ;

	      st->modFN[num1] = modFN ;

	      st->OldCons_DensPhot[num1] = SphP[particle].OldCons_DensPhot[num1] ;
	      st->RT_F[num1][0] = SphP[particle].Cons_RT_F[num1][0]/SphP[particle].Volume ;
	      st->RT_F[num1][1] = SphP[particle].Cons_RT_F[num1][1]/SphP[particle].Volume ;
	      st->RT_F[num1][2] = SphP[particle].Cons_RT_F[num1][2]/SphP[particle].Volume ;
	    }
	}
      

      st->grad = &SphP[particle].Grad;
      st->rtgrad = &SphP[particle].RTGrad;

      st->timeBin = P[particle].TimeBinHydro;

      st->volume = SphP[particle].Volume;


      aBegin = SphP[particle].TimeLastPrimUpdate;

      st->surfacearea = SphP[particle].SurfaceArea;
      st->activearea = SphP[particle].ActiveArea;
      st->csnd = get_sound_speed(particle);
      st->ID = P[particle].ID;
    }
  else
    {
      st->velGas[0] = PrimExch[particle].VelGas[0];
      st->velGas[1] = PrimExch[particle].VelGas[1];
      st->velGas[2] = PrimExch[particle].VelGas[2];

      st->velVertex[0] = PrimExch[particle].VelVertex[0];
      st->velVertex[1] = PrimExch[particle].VelVertex[1];
      st->velVertex[2] = PrimExch[particle].VelVertex[2];

      st->rho = PrimExch[particle].Density ;

      for(int num1=0;num1<MRT_BINS;num1++)
	{
	  st->DensPhot[num1] = RTPrimExch[particle].DensPhot[num1] ;
	  st->OldCons_DensPhot[num1] = RTPrimExch[particle].OldCons_DensPhot[num1] ;
	  st->RT_F[num1][0] = RTPrimExch[particle].RT_F[num1][0] ;
	  st->RT_F[num1][1] = RTPrimExch[particle].RT_F[num1][1] ;
	  st->RT_F[num1][2] = RTPrimExch[particle].RT_F[num1][2] ;
	  st->modFN[num1] = RTPrimExch[particle].modFN[num1] ;
	}

      
      st->grad = &GradExch[particle];
      st->rtgrad = &RTGradExch[particle];

      st->timeBin = PrimExch[particle].TimeBinHydro;    /* This is the hydro timestep */

      st->volume = PrimExch[particle].Volume;

      aBegin = PrimExch[particle].TimeLastPrimUpdate;

      st->surfacearea = PrimExch[particle].SurfaceArea;
      st->activearea = PrimExch[particle].ActiveArea;
      st->csnd = PrimExch[particle].Csnd;
      st->ID = DP[p].ID;
    }

  st->dtExtrapolation = All.Time - aBegin;


  /* check for reflecting or outflowing boundaries */
  face_boundary_check_vertex(T, p, &st->velVertex[0], &st->velVertex[1], &st->velVertex[2]);



  return 0;
}


double face_timestep_RT(struct state *state_L, struct state *state_R, double *hubble_a, double *atime)
{
  //  integertime ti_begin_L, ti_begin_R;
  short int timeBin;
  double face_dt;

  double frac = 1.0 / ((double) (All.RTNumSubCycles)) ;


  /* take the minimum of the two */
  timeBin = state_L->timeBin;
  if(timeBin > state_R->timeBin)
    timeBin = state_R->timeBin;

  /* determine most recent start of the time bins */
  //ti_begin_L = ((All.Ti_Current >> state_L->timeBin) << state_L->timeBin) + (((integertime) 1) << timeBin) ;
  //ti_begin_R = ((All.Ti_Current >> state_R->timeBin) << state_R->timeBin) + (((integertime) 1) << timeBin) ;

  /* compute the half-step prediction times */

  if(All.ComovingIntegrationOn)
    {
      /* calculate scale factor at middle of timestep */
      *atime = All.TimeBegin * exp((All.Ti_Current + (((integertime) 1) << (timeBin - 1))) * All.Timebase_interval);
      *hubble_a = hubble_function(*atime);
    }
  else
    *atime = *hubble_a = 1.0;

  /* set the actual time-step for the face */
  face_dt = (((integertime) 1) << timeBin) * All.Timebase_interval;

  if(All.ComovingIntegrationOn)
    {
      /* converts to delta_t */
      state_L->dt_half /= *hubble_a;
      state_R->dt_half /= *hubble_a;
      face_dt /= *hubble_a;

      face_dt /= *atime;        /* we need dt/a, the (1/a) takes care of the gradient in the cosmological euler equations */

      state_L->dtExtrapolation /= *hubble_a;
      state_L->dtExtrapolation /= *atime;
      state_R->dtExtrapolation /= *hubble_a;
      state_R->dtExtrapolation /= *atime;
    }

  state_L->dt_half = frac*0.5*face_dt ;
  state_R->dt_half = frac*0.5*face_dt ; 

  return frac*face_dt;
}



void state_convert_to_local_frame_RT(struct state *st, double *vel_face, double hubble_a, double atime)
{
  if(All.ComovingIntegrationOn)
    {
      st->velGas[0] /= atime;   /* convert to peculiar velocity */
      st->velGas[1] /= atime;
      st->velGas[2] /= atime;
    }

  st->velx = st->velGas[0] - vel_face[0];
  st->vely = st->velGas[1] - vel_face[1];
  st->velz = st->velGas[2] - vel_face[2];

  if(All.ComovingIntegrationOn)
    {
      st->velx -= atime * hubble_a * st->dx;    /* need to get the physical velocity relative to the face */
      st->vely -= atime * hubble_a * st->dy;
      st->velz -= atime * hubble_a * st->dz;
    }
}


void face_do_time_extrapolation_RT(struct state *delta, struct state *st, double atime)
{
  /* st is the state at the center of the cell */


#if defined (MESHRELAX) || defined (DISABLE_TIME_EXTRAPOLATION)
  /* do not time extrapolation */
  (void) st;
  (void) atime;
  memset(delta, 0, sizeof(struct state));
  return;
#endif

  struct grad_data *grad = st->grad;
  struct rt_grad_data *rtgrad = st->rtgrad; 

#ifndef MUSCL_HANCOCK           /* we use the default Runge Kutta time integration scheme */
  double dt_half = st->dtExtrapolation;
#else
  double dt_half = st->dt_half;
#endif

  if(All.ComovingIntegrationOn)
    dt_half /= atime;

  delta->velx = -dt_half * (1.0 / st->rho * grad->dpress[0] + st->velx * grad->dvel[0][0] + st->vely * grad->dvel[0][1] + st->velz * grad->dvel[0][2]);

  delta->vely = -dt_half * (1.0 / st->rho * grad->dpress[1] + st->velx * grad->dvel[1][0] + st->vely * grad->dvel[1][1] + st->velz * grad->dvel[1][2]);
  
  delta->velz = -dt_half * (1.0 / st->rho * grad->dpress[2] + st->velx * grad->dvel[2][0] + st->vely * grad->dvel[2][1] + st->velz * grad->dvel[2][2]);


 
#ifdef MRT_TIME_EXTRAPOLATION
  //terminate("Entered here\n") ;
  int num1 ;
  for(num1=0;num1<MRT_BINS;num1++)
    { 
      delta->DensPhot[num1] = -dt_half * delta->divF[num1] ;
      delta->RT_F[num1][0] = delta->RT_F[num1][1] = delta->RT_F[num1][2] = 0.0 ;
      /*      delta->RT_F[num1][0] = -dt_half * c_internal_units*c_internal_units*(rtgrad->dPT[num1][0][0][0] + rtgrad->dPT[num1][1][0][1] + rtgrad->dPT[num1][2][0][2]) ;
      delta->RT_F[num1][1] = -dt_half * c_internal_units*c_internal_units*(rtgrad->dPT[num1][0][1][0] + rtgrad->dPT[num1][1][1][1] + rtgrad->dPT[num1][2][1][2]) ;
      delta->RT_F[num1][2] = -dt_half * c_internal_units*c_internal_units*(rtgrad->dPT[num1][0][2][0] + rtgrad->dPT[num1][1][2][1] + rtgrad->dPT[num1][2][2][2]) ;*/
	delta->modFN[num1] = 0.0 ;
    }
#else/*Time extrapolation induces too much noise into the solution. Set to zero for now*/
  int num1 ;
  for(num1=0;num1<MRT_BINS;num1++)
    {
      delta->DensPhot[num1] = 0.0 ;
      delta->RT_F[num1][0] = 0.0 ;
      delta->RT_F[num1][1] = 0.0 ;
      delta->RT_F[num1][2] = 0.0 ;
      delta->modFN[num1] = 0.0 ;
    }
#endif
}

void face_do_spatial_extrapolation_RT(struct state *delta, struct state *st, struct state *st_other)
{
#ifdef DISABLE_SPATIAL_RECONSTRUCTION
  memset(delta, 0, sizeof(struct state));
  return;
#endif


  struct grad_data *grad = st->grad;
  struct rt_grad_data *rtgrad = st->rtgrad;

  double dx[3];
  dx[0] = st->dx;
  dx[1] = st->dy;
  dx[2] = st->dz;
 
  double r[3];
  r[0] = -st_other->dx + st->dx;
  r[1] = -st_other->dy + st->dy;
  r[2] = -st_other->dz + st->dz;

  face_do_spatial_extrapolation_single_quantity(&delta->velx, st->velx, st_other->velx, grad->dvel[0], dx, r);
  face_do_spatial_extrapolation_single_quantity(&delta->vely, st->vely, st_other->vely, grad->dvel[1], dx, r);
  face_do_spatial_extrapolation_single_quantity(&delta->velz, st->velz, st_other->velz, grad->dvel[2], dx, r);
  
  int num1 ;
  for(num1=0;num1<MRT_BINS;num1++)
    {
      face_do_spatial_extrapolation_single_quantity(&delta->DensPhot[num1], st->DensPhot[num1], st_other->DensPhot[num1], rtgrad->dDensPhot[num1], dx, r) ;
      face_do_spatial_extrapolation_single_quantity(&delta->modFN[num1], st->modFN[num1], st_other->modFN[num1], rtgrad->dmodFN[num1], dx, r) ;
      face_do_spatial_extrapolation_single_quantity(&delta->RT_F[num1][0], st->RT_F[num1][0]/st->DensPhot[num1], st_other->RT_F[num1][0]/st_other->DensPhot[num1], rtgrad->dFN[num1][0], dx, r) ;
      face_do_spatial_extrapolation_single_quantity(&delta->RT_F[num1][1], st->RT_F[num1][1]/st->DensPhot[num1], st_other->RT_F[num1][1]/st_other->DensPhot[num1], rtgrad->dFN[num1][1], dx, r) ;
      face_do_spatial_extrapolation_single_quantity(&delta->RT_F[num1][2], st->RT_F[num1][2]/st->DensPhot[num1], st_other->RT_F[num1][2]/st_other->DensPhot[num1], rtgrad->dFN[num1][2], dx, r) ;
    }
}

void face_add_extrapolations_RT(struct state *st_face, struct state *delta_time, struct state *delta_space)
{
  if(st_face->rho <= 0)
    return;

#if !defined(MESHRELAX) && !defined(DISABLE_TIME_EXTRAPOLATION) 
  delta_time->flag = 1 ;
  face_add_extrapolation_RT(st_face, delta_time);
#endif

#if !defined(DISABLE_SPATIAL_EXTRAPOLATION)
  delta_space->flag = 0 ;
  face_add_extrapolation_RT(st_face, delta_space);
#endif

}

 void face_add_extrapolation_RT(struct state *st_face, struct state *delta)
{

  st_face->velx += delta->velx;
  st_face->vely += delta->vely;
  st_face->velz += delta->velz;

  int num1 ;
  double Ncell, Nface, Fx_face, Fy_face, Fz_face, FN_face ;
  if(delta->flag)
    {
      //printf("Entered Here time\n") ;
      for(num1=0;num1<MRT_BINS;num1++)
	{
#ifdef MRT_IR_PHOTON_TRAPPING
	  if(num1<UV_BINS || num1>=(UV_BINS+IR_BINS))
	    {
#endif

	  //	  Ncell = st_face->DensPhot[num1] ;
	  Nface = st_face->DensPhot[num1] + delta->DensPhot[num1] ;

	  if(Nface<MINDENSPHOT)
	    Nface = MINDENSPHOT ;

	  Fx_face = st_face->RT_F[num1][0] + delta->RT_F[num1][0] ;
	  Fy_face = st_face->RT_F[num1][1] + delta->RT_F[num1][1] ;
	  Fz_face = st_face->RT_F[num1][2] + delta->RT_F[num1][2] ;


	   if(isnan(st_face->RT_F[num1][0]))
	     terminate("TIME, Initial F is a nan %d \t %g \t%g \n", num1, st_face->RT_F[num1][0], delta->RT_F[num1][0]) ;


	  double oldmodFN = sqrt(st_face->RT_F[num1][0]*st_face->RT_F[num1][0] + st_face->RT_F[num1][1]*st_face->RT_F[num1][1] + st_face->RT_F[num1][2]*st_face->RT_F[num1][2])/st_face->DensPhot[num1] ;
	  double newmodFN = sqrt(Fx_face*Fx_face + Fy_face*Fy_face +Fz_face*Fz_face)/Nface ;

	  if(oldmodFN == 0.0)
	    oldmodFN = 1.0 ;

	  if(newmodFN == 0.0)
	    newmodFN = 1.0 ;

	  st_face->DensPhot[num1] = Nface ;
	  st_face->RT_F[num1][0] = Fx_face * oldmodFN / newmodFN ;
	  st_face->RT_F[num1][1] = Fy_face * oldmodFN / newmodFN ;
	  st_face->RT_F[num1][2] = Fz_face * oldmodFN / newmodFN ;     


	   if(isnan(st_face->RT_F[num1][0]))
	     terminate("TIME, Initial F is a nan %g %g %g \t%g \t %g \t %g \t %g \n", st_face->RT_F[num1][0], st_face->RT_F[num1][1], st_face->RT_F[num1][2], Fx_face, oldmodFN, newmodFN, st_face->DensPhot[num1]) ;

#ifdef MRT_IR_PHOTON_TRAPPING
	}
#endif

	}
    }
  else
    {      
      //printf("Entered Here space\n") ;
      for(num1=0;num1<MRT_BINS;num1++)
	{
#ifdef MRT_IR_PHOTON_TRAPPING
	  if(num1<UV_BINS || num1>=(UV_BINS+IR_BINS))
	    {
#endif

	  Ncell = st_face->DensPhot[num1] ;
	  Nface = st_face->DensPhot[num1] + delta->DensPhot[num1] ;
	  
	  if(Nface<MINDENSPHOT)
	    Nface = MINDENSPHOT ;
	  
	  Fx_face = st_face->RT_F[num1][0] ;//*Nface/Ncell + delta->RT_F[num1][0]*Nface ;
	  Fy_face = st_face->RT_F[num1][1] ;//*Nface/Ncell + delta->RT_F[num1][1]*Nface ;
	  Fz_face = st_face->RT_F[num1][2] ;//*Nface/Ncell + delta->RT_F[num1][2]*Nface ;
	  FN_face = st_face->modFN[num1] + delta->modFN[num1] ;
	  
	   if(isnan(st_face->RT_F[num1][0]))
	     terminate("SPACE - Initial F is a nan %d \t %g \t%g \t%g\n", num1, Nface, Ncell, st_face->RT_F[num1][0]) ;
	  
	  
	   if(isnan(Fx_face))
	     terminate(" Fx is a nan %g \t %g \t %g||||| Nface = %g \t Ncell = %g \n", Fx_face, st_face->RT_F[num1][0], delta->RT_F[num1][0], Nface, Ncell) ;
      
	  
	  if(FN_face > c_internal_units)
	    FN_face = 0.99 * c_internal_units ;
	  
	  double FN_ulim ;
	  
	  if(Nface != 0.0)
	    FN_ulim = sqrt(Fx_face*Fx_face + Fy_face*Fy_face + Fz_face*Fz_face)/Nface ;
	  else
	    FN_ulim = 1.0 ;
	  
	  if(FN_ulim == 0.0)
	    FN_ulim = 1.0 ;
	  
	  st_face->DensPhot[num1] = Nface ;
	  st_face->RT_F[num1][0] = Fx_face * FN_face / FN_ulim ;
	  st_face->RT_F[num1][1] = Fy_face * FN_face / FN_ulim ;
	  st_face->RT_F[num1][2] = Fz_face * FN_face / FN_ulim ;
#ifdef MRT_IR_PHOTON_TRAPPING
	}
#endif
	}
    }
}



void face_turn_velocities_RT(struct state *st, struct geometry *geom)
{
  double velx, vely, velz;

  velx = st->velx;
  vely = st->vely;
  velz = st->velz;

  st->velx = velx * geom->nx + vely * geom->ny + velz * geom->nz;
  st->vely = velx * geom->mx + vely * geom->my + velz * geom->mz;
  st->velz = velx * geom->px + vely * geom->py + velz * geom->pz;

  /*Turn Velocities*/

  double Fx, Fy, Fz ;  

  int j, k, num1 ;

  for(num1=0;num1<MRT_BINS;num1++)
    {
      Fx = st->RT_F[num1][0] ;
      Fy = st->RT_F[num1][1] ;
      Fz = st->RT_F[num1][2] ;
  
      st->RT_F[num1][0] = Fx * geom->nx + Fy * geom->ny + Fz * geom->nz;
      st->RT_F[num1][1] = Fx * geom->mx + Fy * geom->my + Fz * geom->mz;
      st->RT_F[num1][2] = Fx * geom->px + Fy * geom->py + Fz * geom->pz;

      double modF, n[3], chi, f ;
      modF = sqrt(st->RT_F[num1][0]*st->RT_F[num1][0] + st->RT_F[num1][1]*st->RT_F[num1][1] + st->RT_F[num1][2]*st->RT_F[num1][2]) ;
      
      
      if(st->DensPhot[num1] == 0.0)
	f = modF/c_internal_units ;
      else
	f = modF/c_internal_units/st->DensPhot[num1] ;
      
      
      if(f>1.0000000001)
	terminate("%g\t%g\t%g\t%g\t%g\t%g\t%g\nAgain!!!! Before = %g %g %g\n", f, st->RT_F[num1][0], st->RT_F[num1][1], st->RT_F[num1][2], modF, st->DensPhot[num1], c_internal_units, Fx, Fy, Fz) ;
      
      chi = (3.0 + 4.0*f*f)/(5.0 + 2.0*sqrt(4.0 - 3.0*f*f)) ;
      
      
      
      if(modF == 0.0)
	modF = 1.0 ;
      
      
      for(j=0;j<3;j++)
	n[j] = st->RT_F[num1][j]/modF ;
      
      /*Calculate the Pressure tensor in turned co-ordinate system*/  
      
      for(j=0;j<3;j++)
	{
	  for(k=0;k<3;k++)
	    {
	      if(j==k)
		st->PT[num1][j][k] =  ((1.0 - chi)/2.0) + ((3.0*chi - 1.0)/2.0)*n[j]*n[k] ;
	      else
		st->PT[num1][j][k] = ((3.0*chi - 1.0)/2.0)*n[j]*n[k] ;
	      
	      if(isnan(st->PT[num1][j][k]))
		{ 
		  terminate("\nWHAT!!!!\n %g\t%g\t%g\t%g\n%g\t%g\t%g\t%g\t%g\t%g\n", st->DensPhot[num1], modF, st->RT_F[num1][0], st->RT_F[num1][1], st->RT_F[num1][2], f, chi, n[0], n[1], n[2]);
		}
	      
	      st->PT[num1][j][k] *= st->DensPhot[num1] ;
	    }
	}
    }

}





#ifndef MRT_COMOVING
void flux_convert_to_lab_frame_RT(struct state *st_L, struct state *st_R, double *vel_face, struct fluxes *flux)
{  
  double N, Fx, Fy, Fz ;
  int num1 ;

  for(num1=0;num1<MRT_BINS;num1++)
    {
      N = 0.5 * (st_L->DensPhot[num1] + st_R->DensPhot[num1]) ;
      
      Fx = 0.5 * (st_L->RT_F[num1][0] + st_R->RT_F[num1][0]) ;
      Fy = 0.5 * (st_L->RT_F[num1][1] + st_R->RT_F[num1][1]) ;
      Fz = 0.5 * (st_L->RT_F[num1][2] + st_R->RT_F[num1][2]) ;
      
      flux->DensPhot[num1] -= N * vel_face[0] ;
      
      flux->RT_F[num1][0] -= Fx * vel_face[0] ;
      flux->RT_F[num1][1] -= Fy * vel_face[0] ;
      flux->RT_F[num1][2] -= Fz * vel_face[0] ;
    }
}
#endif

void face_turn_momentum_flux_RT(struct fluxes *flux, struct geometry *geom)
{

  int num1 ;

  for(num1=0;num1<MRT_BINS;num1++)
    {
      double Fx = flux->RT_F[num1][0] ;
      double Fy = flux->RT_F[num1][1] ;
      double Fz = flux->RT_F[num1][2] ;
      
      flux->RT_F[num1][0] = Fx*geom->nx + Fy*geom->mx + Fz*geom->px ;
      flux->RT_F[num1][1] = Fx*geom->ny + Fy*geom->my + Fz*geom->py ;
      flux->RT_F[num1][2] = Fx*geom->nz + Fy*geom->mz + Fz*geom->pz ;     
    }
}




void face_limit_fluxes_RT(struct state *st_L, struct state *st_R, struct state *st_center_L, struct state *st_center_R, struct fluxes *flux, double dt, double *count, double *count_reduced)
{
  *count = *count + 1.0;

  /* choose upwind mass to determine a stability bound on the maximum allowed mass exchange,
     (we do this to prevent negative masses under all circumstances) */

  double upwind_activearea, reduc_fac;

  double upwind_cons_densphot;

  integertime upwind_timebin, downstream_timebin;

  for(int num1=0; num1<MRT_BINS; num1++)
    {
      if(flux->DensPhot[num1] > 0)
	{
	  upwind_cons_densphot = st_L->OldCons_DensPhot[num1];
	  upwind_activearea = st_L->activearea;
	  upwind_timebin = st_L->timeBin;
	  downstream_timebin = st_R->timeBin;
	}
      else
	{
	  upwind_cons_densphot = st_R->OldCons_DensPhot[num1];
	  upwind_activearea = st_R->activearea;
	  upwind_timebin = st_R->timeBin;
	  downstream_timebin = st_L->timeBin;
	}

      if(upwind_timebin > downstream_timebin)
	dt *= pow(2, upwind_timebin - downstream_timebin);
      
      if(fabs(flux->DensPhot[num1] * dt * upwind_activearea) > 0.9 * upwind_cons_densphot)
	{
	  /*#if defined(BOUNDARY_REFL_FLUIDSIDE_MINID) && defined(BOUNDARY_REFL_FLUIDSIDE_MAXID) && defined(BOUNDARY_REFL_SOLIDSIDE_MINID) && defined(BOUNDARY_REFL_SOLIDSIDE_MAXID)
	  MyIDType idL = st_L->ID;
	  MyIDType idR = st_R->ID;

	  if((idR >= BOUNDARY_REFL_FLUIDSIDE_MINID && idR < BOUNDARY_REFL_FLUIDSIDE_MAXID && idL >= BOUNDARY_REFL_SOLIDSIDE_MINID && idL < BOUNDARY_REFL_SOLIDSIDE_MAXID)  || (idL >= BOUNDARY_REFL_FLUIDSIDE_MINID && idL < BOUNDARY_REFL_FLUIDSIDE_MAXID && idR >= BOUNDARY_REFL_SOLIDSIDE_MINID && idR < BOUNDARY_REFL_SOLIDSIDE_MAXID))
	    {
	      reduc_fac = 1.0 ;
	    }
	  else
	    {
	      reduc_fac = 0.9 * upwind_cons_densphot / fabs(flux->DensPhot[num1] * dt * upwind_activearea);
	      }*/
	  //#else
	  reduc_fac = 0.9 * upwind_cons_densphot / fabs(flux->DensPhot[num1] * dt * upwind_activearea);
	  //#endif
	  
	  *count_reduced = *count_reduced + 1.0;

	  flux->DensPhot[num1] *= reduc_fac;
	  flux->RT_F[num1][0] *= reduc_fac;
	  flux->RT_F[num1][1] *= reduc_fac;
	  flux->RT_F[num1][2] *= reduc_fac;
	}
    }
}


void apply_flux_list_RT(void)
{
  int i, j, p, nimport, ngrp, recvTask;
#if defined(MAXSCALARS) || defined(TGCHEM)
  int k;
#endif

  /* now exchange the flux-list and apply it when needed */

  mysort(FluxList, Nflux, sizeof(struct flux_list_data), flux_list_data_compare);

  for(j = 0; j < NTask; j++)
    Send_count[j] = 0;

  for(i = 0; i < Nflux; i++)
    Send_count[FluxList[i].task]++;

  if(Send_count[ThisTask] > 0)
    terminate("Send_count[ThisTask]");

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  struct flux_list_data *FluxListGet = (struct flux_list_data *) mymalloc("FluxListGet", nimport * sizeof(struct flux_list_data));

  /* exchange particle data */
  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&FluxList[Send_offset[recvTask]],
                           Send_count[recvTask] * sizeof(struct flux_list_data), MPI_BYTE,
                           recvTask, TAG_DENS_A,
                           &FluxListGet[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct flux_list_data), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  /* apply the fluxes */

  for(i = 0; i < nimport; i++)
    {
      p = FluxListGet[i].index;



#if defined(BOUNDARY_INFLOWOUTFLOW_MINID) && defined(BOUNDARY_INFLOWOUTFLOW_MAXID)
      if(P[p].ID >= BOUNDARY_INFLOWOUTFLOW_MINID && P[p].ID < BOUNDARY_INFLOWOUTFLOW_MAXID)
        continue;
#endif

#if defined(BOUNDARY_REFL_SOLIDSIDE_MINID) && defined(BOUNDARY_REFL_SOLIDSIDE_MAXID)
      if(P[p].ID >= BOUNDARY_REFL_SOLIDSIDE_MINID && P[p].ID < BOUNDARY_REFL_SOLIDSIDE_MAXID)
        continue;
#endif


      for(int num1=0; num1<MRT_BINS; num1++)
	{
	  SphP[p].Cons_DensPhot[num1] += FluxListGet[i].dDensPhot[num1] ;
	  SphP[p].Cons_RT_F[num1][0] += FluxListGet[i].dRT_F[num1][0] ;
	  SphP[p].Cons_RT_F[num1][1] += FluxListGet[i].dRT_F[num1][1] ;
	  SphP[p].Cons_RT_F[num1][2] += FluxListGet[i].dRT_F[num1][2] ;
	}
    }
  myfree(FluxListGet);
}

#endif
