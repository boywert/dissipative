/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/finite_volume_solver.c
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
#ifdef MCMA
#include <assert.h>
#include "SGChem/f2c.h"
#include "SGChem/sgchem_def.h"
#endif

#include "allvars.h"
#include "proto.h"
#include "voronoi.h"

#include "helm_eos.h"

#ifndef DG

static struct flux_list_data
{
  int task, index;
  double dM, dP[3];
#ifdef MHD
  double dB[3];
#ifdef MHD_CT
  double dA[3];
#endif
#ifdef MHD_THERMAL_ENERGY_SWITCH
  double dEtherm;
#endif
#endif
  

#ifdef MRT
  double dnHI, dnHII, dnHeI, dnHeII, dnHeIII, dne ;
#endif
  
#ifndef ISOTHERM_EQS
  double dEnergy;
#endif
#ifdef USE_ENTROPY_FOR_COLD_FLOWS
  double dEntropy;
#endif
#ifdef TRACER_FIELD
  double dConservedTracer;
#endif
#ifdef MAXSCALARS
  double dConservedScalars[MAXSCALARS];
#endif
#ifdef TRACER_MC
  int pother_task;
  int pother_index;
  MyIDType pother_ID;
#endif
#ifdef OUTPUT_CELL_SPIN
  double dSpin[3];
  double dCenterOffsetMass[3];
#endif
#ifdef TGCHEM
  double Abund[TGCHEM_NUM_ABUNDANCES];
#endif
#ifdef COSMIC_RAYS
  double dCR_Energy;
#endif
#ifdef FLD
  double dFLD;
#endif
#ifdef GFM_CHEMTAGS
  double dMassMetalsChemTags[GFM_N_CHEM_TAGS];
#endif
}
 *FluxList;

static int Nflux, MaxNflux;

struct primexch *PrimExch;
struct grad_data *GradExch;
#ifdef TVD_SLOPE_LIMITER
struct grad_data *GradExchUl;
#endif
#ifdef RT_ADVECT
struct rt_primexch *RTPrimExch;
struct rt_grad_data *RTGradExch;
#endif
#ifdef SECOND_DERIVATIVES
struct hessian_data *HessianExch;
#endif
#ifdef TRACER_MC
static MyFloat *CellTracerMasses;
#endif
#ifdef TGCHEM
static double *AbundFluxes;
#endif

/*!< state on a face determined by riemann solver */
struct state_face state_face;

/*!< flux through a face */
struct fluxes fluxes;
#if defined(VISCOSITY) || defined(THERMAL_CONDUCTION) || defined(TRACER_DIFFUSION)
struct fluxes diffusionfluxes;
#endif

struct geometry geom;


#ifdef ACTIVE_CELL_SPIN
void face_calc_spin_velocities(struct state *st);
void face_add_spin_velocities(struct state *state_L, struct state *state_R);
#endif

#ifdef ONEDIMS_SPHERICAL
void apply_spherical_source_terms();
#endif

static void face_add_extrapolation_with_check(struct state *st_face, struct state *delta, struct fvs_stat *stat);
static void fvs_initialize_statistics(struct fvs_stat *stat);
static void fvs_evaluate_statistics(struct fvs_stat *stat);

#ifdef VORONOI_BACKUP_RESTORE_FACE_AREAS
void backup_face_areas(tessellation * T);
void restore_face_areas(tessellation * T);
#endif

void compute_interface_fluxes(tessellation * T)
{

 
#ifdef NOHYDRO
  return;
#endif
  TIMER_START(CPU_FLUXES);

#if defined(DO_NOT_MOVE_GAS) || defined(OTVET_NOTMOVEGAS)
  TIMER_STOP(CPU_FLUXES);
  return;                       /* do noting in this case */
#endif

#ifdef COSMIC_RAYS_EXTRA_DIAGNOSTICS
  double InitialCREnergy = 0;
  for(int i=0; i<NumGas; i++)
    if(P[i].Mass != 0 && P[i].ID != 0 && P[i].Type == 0)
      InitialCREnergy += SphP[i].CR_Energy;
#endif

  int i, j;
  double count = 0, count_reduced = 0, tot_count, tot_count_reduced;
  double face_dt, hubble_a, atime;
  struct fvs_stat stat;
#ifdef MHD
  double sqrtatime;
#endif

#ifdef GODUNOV_STATS
  FILE *fdstats;
  char buf[1000];

  sprintf(buf, "%s/godunov_stats_%d.txt", All.OutputDir, ThisTask);
  if(!(fdstats = fopen(buf, "w")))
    terminate("error in opening file '%s'", buf);
#endif

#ifdef VORONOI_BACKUP_RESTORE_FACE_AREAS
  backup_face_areas(T);
#endif

  fvs_initialize_statistics(&stat);

  MaxNflux = T->Indi.AllocFacNflux;
  Nflux = 0;
  FluxList = mymalloc_movable(&FluxList, "FluxList", MaxNflux * sizeof(struct flux_list_data));

  face *VF = T->VF;
  point *DP = T->DP;

#ifdef TRACER_MC
  start_MC_tracer(N_tracer);

  CellTracerMasses = (MyFloat *) mymalloc("CellTracerMasses", NumGas * sizeof(MyFloat));
#ifdef AMR
  for(i = 0; i < NumGas; i++)
    CellTracerMasses[i] = P[i].Mass;
#else
  for(i = 0; i < NumGasInMesh; i++)
    {
      int p = List_InMesh[i];
      CellTracerMasses[p] = P[p].Mass;
    }
#endif
#endif

#ifdef TGCHEM
  AbundFluxes = (double *) mymalloc("AbundFluxes", NumGas * TGCHEM_NUM_ABUNDANCES * sizeof(double));

  memset(AbundFluxes, 0, NumGas * TGCHEM_NUM_ABUNDANCES * sizeof(double));
#endif

  for(i = 0; i < T->Nvf; i++)
    {
      struct state state_L, state_center_L, delta_time_L, delta_space_L;
      struct state state_R, state_center_R, delta_time_R, delta_space_R;

#if defined(SECOND_DERIVATIVES) && defined(RECONSTRUCT_GRADIENTS)
      struct grad_data delta_grad_time_L, delta_grad_space_L;
      struct grad_data delta_grad_time_R, delta_grad_space_R;
#endif

      face_dt = 0;              /* the default is that this face is not active */

      /* calculate normal vectors */
      if(face_get_normals(T, i, &geom))
        continue;

      /* get the values of the states at the center of the cells */
      if(face_get_state(T, VF[i].p1, i, &state_center_L))
        continue;

      if(face_get_state(T, VF[i].p2, i, &state_center_R))
        continue;

      /* only treat faces where one of the two sides is active */
      if(!TimeBinSynchronized[state_center_L.timeBin] && !TimeBinSynchronized[state_center_R.timeBin])
        continue;

      /* clarify whether the face should be done by this task (it may be present also on another task) */
      if(face_check_responsibility_of_this_task(T, VF[i].p1, VF[i].p2, &state_center_L, &state_center_R))
        continue;

      /* calculate timestep of the face */
      face_dt = face_timestep(&state_center_L, &state_center_R, &hubble_a, &atime);
#ifdef MHD
      sqrtatime = sqrt(atime);
#endif

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

#if defined(RIEMANN_HLLC) || defined(RIEMANN_ROSUNOV) || defined(RIEMANN_HLLD) || defined(RIEMANN_HLL) || defined(SPECIAL_RELATIVITY) || defined(GENERAL_RELATIVITY) 
      double vel_face_turned[3];
      /* for these riemann solvers, the riemann problem is not solved in the restframe of
         the face, instead the mesh motion is accounted for via an advection step */

      /* turn the face velocity */
      vel_face_turned[0] = vel_face[0] * geom.nx + vel_face[1] * geom.ny + vel_face[2] * geom.nz;
      vel_face_turned[1] = vel_face[0] * geom.mx + vel_face[1] * geom.my + vel_face[2] * geom.mz;
      vel_face_turned[2] = vel_face[0] * geom.px + vel_face[1] * geom.py + vel_face[2] * geom.pz;
#endif


#if defined(VISCOSITY) && defined(SECOND_DERIVATIVES)
      /*modify the predictor step with a viscous kick (source term) */
      double local_dynvisc = 0;
      double local_bulkvisc = 0;
      face_extrapolate_viscous_kick(&state_center_L, local_dynvisc, local_bulkvisc);
      face_extrapolate_viscous_kick(&state_center_R, local_dynvisc, local_bulkvisc);
#endif
#if defined(SECOND_DERIVATIVES) && defined(RECONSTRUCT_GRADIENTS)
      face_time_advance_gradients(&delta_grad_time_L, &state_center_L);
      face_time_advance_gradients(&delta_grad_time_R, &state_center_R);

      face_space_extrapolate_gradients(&delta_grad_space_L, &state_center_L);
      face_space_extrapolate_gradients(&delta_grad_space_R, &state_center_R);
#endif

      state_convert_to_local_frame(&state_center_L, vel_face, hubble_a, atime);
      state_convert_to_local_frame(&state_center_R, vel_face, hubble_a, atime);



      /* copy center state to state at interface, then add extrapolation terms */
      state_L = state_center_L;
      state_R = state_center_R;


     

      face_do_time_extrapolation(&delta_time_L, &state_center_L, atime);
      face_do_time_extrapolation(&delta_time_R, &state_center_R, atime);

      face_do_spatial_extrapolation(&delta_space_L, &state_center_L, &state_center_R);
      face_do_spatial_extrapolation(&delta_space_R, &state_center_R, &state_center_L);

      face_add_extrapolations(&state_L, &delta_time_L, &delta_space_L, &stat);
      face_add_extrapolations(&state_R, &delta_time_R, &delta_space_R, &stat);


#ifdef JEANS_RIEMANN_PRESSURE_LIMIT
      double jeans_pressure, celldim, ncells, ulimit;
      ncells = JEANS_RIEMANN_PRESSURE_LIMIT;
      celldim = 2.0 * pow(All.MinVolume * 3.0 / 4.0 / M_PI / 2.0, 1.0 / 3.0);

      jeans_pressure = ncells * ncells * All.G * celldim * celldim * state_L.rho * state_L.rho / M_PI / GAMMA;
      ulimit = jeans_pressure / (GAMMA_MINUS1) / state_L.rho;
      state_L.press = fmax(state_L.press, jeans_pressure);

      jeans_pressure = ncells * ncells * All.G * celldim * celldim * state_R.rho * state_R.rho / M_PI / GAMMA;
      ulimit = jeans_pressure / (GAMMA_MINUS1) / state_R.rho;
      state_R.press = fmax(state_R.press, jeans_pressure);
#endif


#ifdef ACTIVE_CELL_SPIN
      face_calc_spin_velocities(&state_L);
      face_calc_spin_velocities(&state_R);
#endif

#if defined(SECOND_DERIVATIVES) && defined(RECONSTRUCT_GRADIENTS)
      face_add_gradient_extrapolations(&state_L, &delta_grad_time_L, &delta_grad_space_L);
      face_add_gradient_extrapolations(&state_R, &delta_grad_time_R, &delta_grad_space_R);
#endif

#ifdef MHD
      if(All.ComovingIntegrationOn)
        {
          state_L.Bx /= sqrtatime;
          state_L.By /= sqrtatime;
          state_L.Bz /= sqrtatime;

          state_R.Bx /= sqrtatime;
          state_R.By /= sqrtatime;
          state_R.Bz /= sqrtatime;
        }
#endif

#ifdef COSMIC_RAYS
      if(All.ComovingIntegrationOn)
        {
          state_L.crPressure /= atime;
          state_R.crPressure /= atime;
        }
#endif

#ifndef MESHRELAX
#ifndef ISOTHERM_EQS
      /* check for crazy values */
      if(state_L.press < 0 || state_R.press < 0 || state_L.rho < 0 || state_R.rho < 0)
        {
          printf("i=%d press_L=%g press_R=%g rho_L=%g rho_R=%g\n", i, state_L.press, state_R.press, state_L.rho, state_R.rho);
          printf("area=%g lx=%g ly=%g   rx=%g ry=%g\n", VF[i].area, state_L.dx, state_L.dy, state_R.dx, state_R.dy);
          terminate("found crazy values");
        }
#else
      if(state_L.press < 0 || state_R.press < 0 || state_L.rho < 0 || state_R.rho < 0)
        {
          printf("i=%d rho_L=%g rho_R=%g\n", i, state_L.rho, state_R.rho);
          printf("area=%g lx=%g ly=%g   rx=%g ry=%g\n", VF[i].area, state_L.dx, state_L.dy, state_R.dx, state_R.dy);
          terminate("found crazy values");
        }
#endif
#endif

      /* mirror velocity in case of reflecting boundaries */
      face_boundary_check(&T->DP[VF[i].p1], &state_L.velx, &state_L.vely, &state_L.velz);
      face_boundary_check(&T->DP[VF[i].p2], &state_R.velx, &state_R.vely, &state_R.velz);

#ifdef MHD
      /* mirror magnetic field in case of reflecting boundaries */
      face_boundary_check(&T->DP[VF[i].p1], &state_L.Bx, &state_L.By, &state_L.Bz);
      face_boundary_check(&T->DP[VF[i].p2], &state_R.Bx, &state_R.By, &state_R.Bz);
#endif

      /* turn the velocities to get velx perpendicular and vely and velz in the plane of the face */
      face_turn_velocities(&state_L, &geom);
      face_turn_velocities(&state_R, &geom);

#ifdef ACTIVE_CELL_SPIN
      face_add_spin_velocities(&state_L, &state_R);
#endif

#ifdef SPECIAL_BOUNDARY
      /*
         if((state_R.ID == -1) && (state_L.ID == -2))
         get_boundary_cell_state(&state_R, &state_L);

         if((state_L.ID == -1) && (state_R.ID == -2))
         get_boundary_cell_state(&state_L, &state_R);

       */
      int orient;
      if((state_R.ID == -1) && (state_L.ID == -2))
        {
          if(state_R.velx < 0)
            orient = 1;
          else
            orient = -1;
          get_boundary_cell_state(&state_R, &state_L, orient);
        }
      if((state_L.ID == -1) && (state_R.ID == -2))
        {
          if(state_L.velx > 0)
            orient = 1;
          else
            orient = -1;
          get_boundary_cell_state(&state_L, &state_R, orient);
        }

      if((state_L.ID <= -3) || (state_R.ID <= -3))
        continue;

      if((state_L.ID == -2) && (state_R.ID == -2))
        continue;

      if((state_L.ID == -2) && (state_R.ID > 0))
        terminate("Not allowed: solid surface cells in contact with normal fluid cells: ID_L=%d ID_R=%d", state_L.ID, state_R.ID);
      if((state_R.ID == -2) && (state_L.ID > 0))
        terminate("Not allowed: solid surface cells in contact with normal fluid cells: ID_R=%d ID_L=%d", state_R.ID, state_L.ID);
#endif


#if defined(BOUNDARY_REFL_FLUIDSIDE_MINID) && defined(BOUNDARY_REFL_FLUIDSIDE_MAXID) && defined(BOUNDARY_REFL_SOLIDSIDE_MINID) && defined(BOUNDARY_REFL_SOLIDSIDE_MAXID)
      //#ifndef MRT
      MyIDType idL = state_L.ID;
      MyIDType idR = state_R.ID;
      //#else
      //idL = state_L.ID;
      //idR = state_R.ID;
      //#endif
      if(idR >= BOUNDARY_REFL_FLUIDSIDE_MINID && idR < BOUNDARY_REFL_FLUIDSIDE_MAXID && idL >= BOUNDARY_REFL_SOLIDSIDE_MINID && idL < BOUNDARY_REFL_SOLIDSIDE_MAXID)
        {
          state_L.velx = -state_R.velx;
          state_L.vely = state_R.vely;
          state_L.velz = state_R.velz;

          state_L.rho = state_R.rho;
          state_L.press = state_R.press;
        }
      else if(idL >= BOUNDARY_REFL_FLUIDSIDE_MINID && idL < BOUNDARY_REFL_FLUIDSIDE_MAXID && idR >= BOUNDARY_REFL_SOLIDSIDE_MINID && idR < BOUNDARY_REFL_SOLIDSIDE_MAXID)
        {
          state_R.velx = -state_L.velx;
          state_R.vely = state_L.vely;
          state_R.velz = state_L.velz;

          state_R.rho = state_L.rho;
          state_R.press = state_L.press;

	  /*   #ifdef MRT
	       state_R.RT_F[nn1][0] = -state_L.RT_F[nn1][0];
	       state_R.RT_F[nn1][1] = state_L.RT_F[nn1][1];
	       state_R.RT_F[nn1][2] = state_L.RT_F[nn1][2];
	       
	       state_R.DensPhot[nn1] = state_L.DensPhot[nn1];
	       state_R.modFN[nn1] = state_L.modFN[nn1] ;
	       #endif*/
	}
#endif


#ifdef COFFEE_PROBLEM
      MyIDType idL = state_L.ID;
      MyIDType idR = state_R.ID;

      if(idR >= 10000000 && idR < 20000000 && idL >= 20000000)
        {
          state_L.velx = -state_R.velx;
          state_L.vely = state_R.vely;
          state_L.velz = state_R.velz;

          state_L.rho = state_R.rho;
          state_L.press = state_R.press;
          state_L.tracer = state_R.tracer;
        }

      if(idL >= 10000000 && idL < 20000000 && idR >= 20000000)
        {
          state_R.velx = -state_L.velx;
          state_R.vely = state_L.vely;
          state_R.velz = state_L.velz;

          state_R.rho = state_L.rho;
          state_R.press = state_L.press;
          state_R.tracer = state_L.tracer;
        }
#endif

#ifndef MESHRELAX

      /* call Riemann solver */

      double press;
#ifdef SPECIAL_RELATIVITY
      press = godunov_flux_3d_hlle_special_relativity(&state_L, &state_R, vel_face_turned, &state_face, &fluxes);
#else
#ifdef GENERAL_RELATIVITY
      press = godunov_flux_3d_hlle_general_relativity(&state_L, &state_R, vel_face_turned, &state_face, &fluxes, &geom, VF[i].cx, VF[i].cy, VF[i].cz);
#else        
#ifdef RIEMANN_GAMMA
      press = godunov_flux_3d_gamma(&state_L, &state_R, &state_face);
#else
#ifdef RIEMANN_HLLC
      press = godunov_flux_3d_hllc(&state_L, &state_R, &state_face, &fluxes);
#else
#ifdef RIEMANN_ROSUNOV
      press = godunov_flux_3d_rosunov(&state_L, &state_R, &state_face, &fluxes);
#else
#ifdef RIEMANN_HLLD
      press = godunov_flux_3d_hlld(&state_L, &state_R, vel_face_turned, &state_face, &fluxes);
#else
#ifdef RIEMANN_HLL
      press = godunov_flux_3d_hll(&state_L, &state_R, &state_face, &fluxes);
#else
      press = godunov_flux_3d(&state_L, &state_R, &state_face); /* exact ideal gas solver */
#endif
#endif
#endif
#endif
#endif
#endif
#endif


      if(press < 0)
        terminate("press < 0: ID_L: %d, ID_R: %d", VF[i].p1, VF[i].p2);

#ifdef TGCHEM
      state_face.gamma = (state_L.gamma + state_R.gamma) / 2;
#endif

#if defined(ENTROPY_MACH_THRESHOLD)  || (TRACER_MC_MACHMAX) || (TRACER_PART_MACHMAX)
#if defined(RIEMANN_HLLD)
      state_L.mach = state_R.mach = 1;  /* no mach numbers yet for RIEMANN_HLLD scheme: FIX ME for entropy scheme */
#else
      get_mach_numbers(&state_L, &state_R, press);
#endif
#endif

#ifdef GODUNOV_STATS
      get_mach_numbers(&state_L, &state_R, press);
      if(st_L.rho > 1.0e-6 && st_R.rho > 1.0e-6)
        fprintf(fdstats, "%g %g %g   %g %g %g  %g %g %g  %g %g %g\n",
                state_L.rho, state_L.velx, state_L.press, state_L.rho, state_L.velx, state_L.press, state_face.rho, state_face.velx, state_face.press, state_L.mach, state_R.mach, VF[i].area);
#endif /* GODUNOV_STATS */

#endif /* end of MESHRELAX */

#ifdef MHD_CT
      double state_face_velx = state_face.velx;
#endif

#if defined(COSMIC_RAYS) || defined(MHD_THERMAL_ENERGY_SWITCH)
      double face_velx = state_face.velx + vel_face_turned[0];
#endif

      /* turn the velocity field back */
      face_turnback_velocities(&state_face, &geom);

      /* add the face velocity again */
      state_face.velx += vel_face[0];
      state_face.vely += vel_face[1];
      state_face.velz += vel_face[2];

#ifndef MESHRELAX

#if defined(RIEMANN_HLLC) || defined(RIEMANN_ROSUNOV) || defined(RIEMANN_HLLD) || defined(RIEMANN_HLL) || defined(SPECIAL_RELATIVITY) || defined(GENERAL_RELATIVITY) /* for non-exact Riemann solver, fluxes are already computed in the local frame, so convert to lab frame and turn momentum fluxes to the lab orientation  */
#if !defined(SPECIAL_RELATIVITY) && !defined(GENERAL_RELATIVITY)
      flux_convert_to_lab_frame(&state_L, &state_R, vel_face_turned, &fluxes);
#endif
      face_turn_momentum_flux(&fluxes, &geom);

#else

      /* calculate fluxes for exact Riemann problem */
      /* compute net flux with dot-product of outward normal and area of face */
      /* multiplication with area and time-step comes later */

      face_get_fluxes(&state_L, &state_R, &state_face, &fluxes, &geom, vel_face);

#endif


      /* set the face states and fluxes of those quantities that are passively advected */
      face_set_scalar_states_and_fluxes(&state_L, &state_R, &state_face, &fluxes);


#if defined(VISCOSITY) || defined(THERMAL_CONDUCTION) || defined(TRACER_DIFFUSION)
      face_get_gradients(&state_L, &state_R, &state_face, &fluxes);
      face_clear_fluxes(&diffusionfluxes);
#ifdef VISCOSITY
      double dynvisc, bulkvisc;
      dynvisc = local_get_dynvisc_coefficient(VF[i].cx, VF[i].cy, VF[i].cz, state_face.rho, state_face.press);
      bulkvisc = local_get_bulkvisc_coefficient(VF[i].cx, VF[i].cy, VF[i].cz, state_face.rho, state_face.press);
      face_get_viscous_fluxes(&state_face, &diffusionfluxes, &geom, dynvisc, bulkvisc);
#endif
#ifdef THERMAL_CONDUCTION
      face_get_conduction_fluxes(&state_face, &diffusionfluxes, &geom);
#endif
#ifdef TRACER_DIFFUSION
      face_get_scalar_diffusion_fluxes(&state_face, &diffusionfluxes, &geom);
#endif
#endif

#ifdef NON_IDEAL_MHD
#ifdef MHD_POWELL
#ifdef AMBIPOLAR_DIFFUSION
      face_add_ambipolar_fluxes(&state_center_L, &state_center_R, &delta_time_L, &delta_time_R, &geom, &fluxes, atime, sqrtatime);
#endif
#ifdef OHMIC_DIFFUSION
      face_add_ohmic_fluxes(&state_center_L, &state_center_R, &delta_time_L, &delta_time_R, &geom, &fluxes, atime, sqrtatime);
#endif
#endif
#ifdef MHD_CT
#ifdef OHMIC_DIFFUSION
      face_add_ohmic_heating(&state_center_L, &state_center_R, &delta_time_L, &delta_time_R, &geom, &fluxes, atime, sqrtatime);
#endif
#endif
#endif


      face_limit_fluxes(&state_L, &state_R, &state_center_L, &state_center_R, &fluxes, face_dt, &count, &count_reduced);

#if defined(BOUNDARY_REFL_ACTS_AS_SOURCE) && defined(BOUNDARY_REFL_FLUIDSIDE_MINID) && defined(BOUNDARY_REFL_FLUIDSIDE_MAXID) && defined(BOUNDARY_REFL_SOLIDSIDE_MINID) && defined(BOUNDARY_REFL_SOLIDSIDE_MAXID)
      {
        MyIDType idL = state_L.ID;
        MyIDType idR = state_R.ID;

        if(idR >= BOUNDARY_REFL_FLUIDSIDE_MINID && idR < BOUNDARY_REFL_FLUIDSIDE_MAXID && idL >= BOUNDARY_REFL_SOLIDSIDE_MINID && idL < BOUNDARY_REFL_SOLIDSIDE_MAXID)
          {
            double fac = (state_center_L.velx) * geom.nx + (state_center_L.vely) * geom.ny + (state_center_L.velz) * geom.nz;
            double dm_flux = state_center_L.rho * fac;

            fluxes.mass += dm_flux;
            fluxes.momentum[0] += state_center_L.velx * dm_flux;
            fluxes.momentum[1] += state_center_L.vely * dm_flux;
            fluxes.momentum[2] += state_center_L.velz * dm_flux;

            fluxes.energy += 0.5 * dm_flux * ((state_center_L.velx * state_center_L.velx +
                                               state_center_L.vely * state_center_L.vely + state_center_L.velz * state_center_L.velz)) + state_center_L.press / (GAMMA_MINUS1) * fac;
          }
        else if(idL >= BOUNDARY_REFL_FLUIDSIDE_MINID && idL < BOUNDARY_REFL_FLUIDSIDE_MAXID && idR >= BOUNDARY_REFL_SOLIDSIDE_MINID && idR < BOUNDARY_REFL_SOLIDSIDE_MAXID)
          {
            double fac = (state_center_R.velx) * geom.nx + (state_center_R.vely) * geom.ny + (state_center_R.velz) * geom.nz;
            double dm_flux = state_center_R.rho * fac;

            fluxes.mass += dm_flux;
            fluxes.momentum[0] += state_center_R.velx * dm_flux;
            fluxes.momentum[1] += state_center_R.vely * dm_flux;
            fluxes.momentum[2] += state_center_R.velz * dm_flux;

            fluxes.energy += 0.5 * dm_flux * ((state_center_R.velx * state_center_R.velx +
                                               state_center_R.vely * state_center_R.vely + state_center_R.velz * state_center_R.velz)) + state_center_R.press / (GAMMA_MINUS1) * fac;
          }
      }
#endif

      /* put in cosmological factors */
      if(All.ComovingIntegrationOn)
        {
          fluxes.momentum[0] *= atime;
          fluxes.momentum[1] *= atime;
          fluxes.momentum[2] *= atime;
          fluxes.energy *= atime * atime;
#ifdef MHD
          fluxes.B[0] *= sqrtatime;
          fluxes.B[1] *= sqrtatime;
          fluxes.B[2] *= sqrtatime;
#ifdef MHD_POWELL
          state_face.Bx *= sqrtatime;
#endif
#endif
#ifdef COSMIC_RAYS
          face_velx *= atime;
#endif
        }

#else /* MESHRELAX */

      /* just solve the advection equation instead of Riemann problem */

      solve_advection(&state_L, &state_R, &state_face, &geom, vel_face);
      face_clear_fluxes(&fluxes);
      face_add_fluxes_advection(&state_face, &fluxes, &geom, vel_face);
      face_set_scalar_states_and_fluxes(&state_L, &state_R, &state_face, &fluxes);

#endif

#ifndef ISOTHERM_EQS
      if(!gsl_finite(fluxes.energy))
        {
          printf("i=%d eFlux-Bummer: %g %g %g\n", i, fluxes.energy, state_face.press, state_face.rho);
          printf("rho_L=%g velx_L=%g vely_L=%g velz_L=%g press_L=%g\n", state_L.rho, state_L.velx, state_L.vely, state_L.velz, state_L.press);
          printf("rho_R=%g velx_R=%g vely_R=%g velz_R=%g press_R=%g\n", state_R.rho, state_R.velx, state_R.vely, state_R.velz, state_R.press);
#ifdef VARIABLE_GAMMA
          printf("gammaE_L=%g gammaE_R=%g, gammaE_face=%g\n", state_L.gammaE, state_R.gammaE, state_face.gammaE);
          print_state_face_info(&state_face);
#endif
          print_particle_info(i);
          terminate("infinity encountered");
        }
#endif


#ifdef MRT_DO_NOT_MOVE_GAS
      fluxes.mass = 0.0 ;
      fluxes.momentum[0] = 0.0 ;
      fluxes.momentum[1] = 0.0 ;
      fluxes.momentum[2] = 0.0 ;
      fluxes.energy = 0.0 ;
#endif

      
      
      /* now apply the flux to update the conserved states of the cells */


      if(face_dt > 0)           /* selects active faces */
        {
          int k, p, q;
          double dir;
          double fac = face_dt * VF[i].area;
#if defined(MAXSCALARS) || defined(TGCHEM)
          int m;
#endif

#ifndef MUSCL_HANCOCK
#ifndef RUNGE_KUTTA_FULL_UPDATE   
          fac *= 0.5;
#else
          fac *= 1.;      
#endif          
#endif

#if defined(MHD_POWELL) && !defined(MHD_POWELL_SPLIT)
          struct state *state_center, *delta_time;
#endif
          for(k = 0; k < 2; k++)
            {
#if defined(TRACER_MC) || defined(MHD_CT) || defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
              int qother;
#endif
              if(k == 0)
                {
                  q = VF[i].p1;
                  p = DP[q].index;
                  dir = -fac;
#if defined(TRACER_MC) || defined(MHD_CT) || defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
                  qother = VF[i].p2;
#endif
#if defined(MHD_POWELL) && !defined(MHD_POWELL_SPLIT)
                  state_center = &state_center_L;
                  delta_time = &delta_time_L;
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
#if defined(MHD_POWELL) && !defined(MHD_POWELL_SPLIT)
                  state_center = &state_center_R;
                  delta_time = &delta_time_R;
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

#ifdef MHD_CT
              get_A_fluxes(k, q, qother, vel_face_turned[0], state_face_velx, &state_L, &state_R, &fluxes);

#ifdef NON_IDEAL_MHD
#ifdef OHMIC_DIFFUSION
              face_add_ohmic_A_fluxes(q, qother, &state_center_L, &state_center_R, &fluxes, atime, sqrtatime);
#endif
#endif

#endif

#ifdef TRACER_MC
              int pother_task = DP[qother].task;
              int pother_index = DP[qother].index;
              MyIDType pother_ID = DP[qother].ID;
              if(pother_task == ThisTask)
                {
                  if(pother_index >= NumGas)
                    pother_index -= NumGas;
                }
              else
                pother_index = DP[qother].originalindex;
#endif

              if(DP[q].task == ThisTask)
                {
                  if(DP[q].index >= NumGas)     /* this is a local ghost point */
                    {
                      if(DP[VF[i].p1].ID == DP[VF[i].p2].ID)    /* this may happen for reflective points */
                        continue;
                      p -= NumGas;
                    }

                  /* note: this will be executed if P[p] is a local point, independent of active or not */
#if defined(ENTROPY_MACH_THRESHOLD) || (TRACER_MC_MACHMAX) || (TRACER_PART_MACHMAX)
                  if(SphP[p].MaxMach < state_L.mach)
                    SphP[p].MaxMach = state_L.mach;
                  if(SphP[p].MaxMach < state_R.mach)
                    SphP[p].MaxMach = state_R.mach;
#endif

#ifdef SPECIAL_BOUNDARY
                  if((P[p].ID == -2) || (P[p].ID <= -3))
                    continue;
#endif

#if defined(BOUNDARY_INFLOWOUTFLOW_MINID) && defined(BOUNDARY_INFLOWOUTFLOW_MAXID)
                  if(P[p].ID >= BOUNDARY_INFLOWOUTFLOW_MINID && P[p].ID < BOUNDARY_INFLOWOUTFLOW_MAXID)
                    continue;
#endif

#if defined(BOUNDARY_REFL_SOLIDSIDE_MINID) && defined(BOUNDARY_REFL_SOLIDSIDE_MAXID)
                  if(P[p].ID >= BOUNDARY_REFL_SOLIDSIDE_MINID && P[p].ID < BOUNDARY_REFL_SOLIDSIDE_MAXID)
                    continue;
#endif

#ifdef COFFEE_PROBLEM
                  if(P[p].ID >= 20000000)
                    continue;
#endif

                  P[p].Mass += dir * fluxes.mass;
                  SphP[p].Momentum[0] += dir * fluxes.momentum[0];
                  SphP[p].Momentum[1] += dir * fluxes.momentum[1];
                  SphP[p].Momentum[2] += dir * fluxes.momentum[2];                 

#ifdef MRT
		  SphP[p].nHI += dir * fluxes.nHI ;
		  SphP[p].nHII += dir * fluxes.nHII ;
		  SphP[p].ne += dir * fluxes.ne ;
#ifdef MRT_INCLUDE_HE
		  SphP[p].nHeI += dir * fluxes.nHeI ;
		  SphP[p].nHeII += dir * fluxes.nHeII ;
		  SphP[p].nHeIII += dir * fluxes.nHeIII ;
#endif
#endif
	

#if defined(TRACER_MC)
                  if(dir * fluxes.mass < 0)
                    {
                      if(pother_index > All.MaxPart)
                        terminate("pother_index > All.MaxPart");

                      consider_moving_tracers(p, pother_task, pother_index, pother_ID, -dir * fluxes.mass / CellTracerMasses[p]);
                      CellTracerMasses[p] += dir * fluxes.mass;
                    }
#endif

#ifdef MHD
                  SphP[p].BConserved[0] += dir * fluxes.B[0];
                  SphP[p].BConserved[1] += dir * fluxes.B[1];
                  SphP[p].BConserved[2] += dir * fluxes.B[2];
#ifdef MHD_CT
                  SphP[p].AConserved[0] += dir * fluxes.A[0];
                  SphP[p].AConserved[1] += dir * fluxes.A[1];
                  SphP[p].AConserved[2] += dir * fluxes.A[2];
#endif
#if defined(MHD_POWELL) && !defined(MHD_POWELL_SPLIT)
                  double Velx = state_center->velx + delta_time->velx + vel_face[0];
                  double Vely = state_center->vely + delta_time->vely + vel_face[1];
                  double Velz = state_center->velz + delta_time->velz + vel_face[2];

                  if(All.ComovingIntegrationOn)
                    {
                      Velx += atime * hubble_a * state_center->dx;
                      Vely += atime * hubble_a * state_center->dy;
                      Velz += atime * hubble_a * state_center->dz;
                    }

                  double Bx = state_center->Bx + delta_time->Bx;
                  double By = state_center->By + delta_time->By;
                  double Bz = state_center->Bz + delta_time->Bz;

                  SphP[p].BConserved[0] += dir * Velx * state_face.Bx;
                  SphP[p].BConserved[1] += dir * Vely * state_face.Bx;
                  SphP[p].BConserved[2] += dir * Velz * state_face.Bx;

                  SphP[p].Momentum[0] += dir * Bx * state_face.Bx;
                  SphP[p].Momentum[1] += dir * By * state_face.Bx;
                  SphP[p].Momentum[2] += dir * Bz * state_face.Bx;

                  SphP[p].Energy += dir * (Bx * Velx + By * Vely + Bz * Velz) * state_face.Bx;

                  {
                    double dMomX = dir * Bx * state_face.Bx;
                    double dMomY = dir * By * state_face.Bx;
                    double dMomZ = dir * Bz * state_face.Bx;

                    All.Powell_Momentum[0] += dMomX;
                    All.Powell_Momentum[1] += dMomY;
                    All.Powell_Momentum[2] += dMomZ;

                    double dx = SphP[p].Center[0] - 0.5 * All.BoxSize;
                    double dy = SphP[p].Center[1] - 0.5 * All.BoxSize;
                    double dz = SphP[p].Center[2] - 0.5 * All.BoxSize;

                    All.Powell_Angular_Momentum[0] += dy * dMomZ - dz * dMomY;
                    All.Powell_Angular_Momentum[1] += dz * dMomX - dx * dMomZ;
                    All.Powell_Angular_Momentum[2] += dx * dMomY - dy * dMomX;
                    All.Powell_Energy += dir * (Bx * Velx + By * Vely + Bz * Velz) * state_face.Bx;
                  }
#endif
#ifdef MHD_THERMAL_ENERGY_SWITCH
                  double ptot = SphP[p].Pressure + 0.5 * (SphP[p].B[0] * SphP[p].B[0] + SphP[p].B[1] * SphP[p].B[1] + SphP[p].B[2] * SphP[p].B[2]);
#ifdef COSMIC_RAYS
                  ptot += SphP[p].CR_Pressure;
#endif
                  SphP[p].Etherm += dir * ptot * face_velx;
#endif
#endif

#ifdef COSMIC_RAYS
                  SphP[p].Energy -= dir * SphP[p].CR_Pressure * face_velx * All.cf_atime;
                  SphP[p].CR_Energy += dir * SphP[p].CR_Pressure * face_velx;
#endif

#ifdef FLD
                  SphP[p].n_gamma += dir * fluxes.dFLD;
#endif

#ifdef MAXSCALARS
                  for(m = 0; m < N_Scalar; m++)
                    {
//                  assert(*(MyFloat *) (((char *) (&SphP[p])) + scalar_elements[m].offset_mass) + dir * fluxes.scalars[m]>=0);
                      *(MyFloat *) (((char *) (&SphP[p])) + scalar_elements[m].offset_mass) += dir * fluxes.scalars[m];
                    }
#endif

#ifdef GFM_CHEMTAGS
                  for(m = 0; m < GFM_N_CHEM_TAGS; m++)
                    SphP[p].MassMetalsChemTags[m] += dir * fluxes.chemtags[m];
#endif

#if !defined(ISOTHERM_EQS) 
                 
                  SphP[p].Energy += dir * fluxes.energy;      
#endif

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
                  SphP[p].Entropy += dir * fluxes.entropy;
#endif

#ifdef TRACER_FIELD
                  SphP[p].ConservedTracer += dir * fluxes.tracer;
#endif

#ifdef OUTPUT_CELL_SPIN
                  double dx, dy, dz;
                  {
                    double xtmp, ytmp, ztmp;
                    dx = NEAREST_X(geom.cx - SphP[p].Center[0] - SphP[p].CenterOffset[0]);
                    dy = NEAREST_Y(geom.cy - SphP[p].Center[1] - SphP[p].CenterOffset[1]);
                    dz = NEAREST_Z(geom.cz - SphP[p].Center[2] - SphP[p].CenterOffset[2]);
                  }

                  SphP[p].Spin[0] += dir * (dy * fluxes.momentum[2] - dz * fluxes.momentum[1]);
                  SphP[p].Spin[1] += dir * (dz * fluxes.momentum[0] - dx * fluxes.momentum[2]);
                  SphP[p].Spin[2] += dir * (dx * fluxes.momentum[1] - dy * fluxes.momentum[0]);

                  SphP[p].CenterOffsetMass[0] += dir * (dx * fluxes.mass);
                  SphP[p].CenterOffsetMass[1] += dir * (dy * fluxes.mass);
                  SphP[p].CenterOffsetMass[2] += dir * (dz * fluxes.mass);
#endif

#if defined(VISCOSITY) || defined(THERMAL_CONDUCTION) || defined(TRACER_DIFFUSION)
                  SphP[p].Momentum[0] -= dir * diffusionfluxes.momentum[0];
                  SphP[p].Momentum[1] -= dir * diffusionfluxes.momentum[1];
                  SphP[p].Momentum[2] -= dir * diffusionfluxes.momentum[2];
                  SphP[p].Energy -= dir * diffusionfluxes.energy;
#if defined(TRACER_FIELD) && defined(TRACER_DIFFUSION)
                  SphP[p].ConservedTracer -= dir * diffusionfluxes.tracer;
#endif
#endif

#ifdef TGCHEM
                  for(m = 0; m < TGCHEM_NUM_ABUNDANCES; m++)
                    AbundFluxes[TGCHEM_NUM_ABUNDANCES * p + m] += dir * fluxes.pcabund[m];
#endif

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
#ifdef VERBOSE
                      printf("Task=%d: increase memory allocation, MaxNflux=%d Indi.AllocFacNflux=%g\n", ThisTask, MaxNflux, T->Indi.AllocFacNflux);
#endif
		      FluxList = myrealloc_movable(FluxList, MaxNflux * sizeof(struct flux_list_data));

                      if(Nflux >= MaxNflux)
                        terminate("Nflux >= MaxNflux");
                    }

                  FluxList[Nflux].task = DP[q].task;
                  FluxList[Nflux].index = DP[q].originalindex;

                  FluxList[Nflux].dM = dir * fluxes.mass;

#ifdef TRACER_MC
                  FluxList[Nflux].pother_task = pother_task;
                  FluxList[Nflux].pother_index = pother_index;
                  FluxList[Nflux].pother_ID = pother_ID;
#endif

                  FluxList[Nflux].dP[0] = dir * fluxes.momentum[0];
                  FluxList[Nflux].dP[1] = dir * fluxes.momentum[1];
                  FluxList[Nflux].dP[2] = dir * fluxes.momentum[2];

#ifdef MRT
		  FluxList[Nflux].dnHI = dir * fluxes.nHI ;
		  FluxList[Nflux].dnHII = dir * fluxes.nHII ;
		  FluxList[Nflux].dne = dir * fluxes.ne ;
#ifdef MRT_INCLUDE_HE
		  FluxList[Nflux].dnHeI = dir * fluxes.nHeI ;
		  FluxList[Nflux].dnHeII = dir * fluxes.nHeII ;
		  FluxList[Nflux].dnHeIII = dir * fluxes.nHeIII ;
#endif
#endif
	
                  
#if !defined(ISOTHERM_EQS) 

                  FluxList[Nflux].dEnergy = dir * fluxes.energy;
#endif

#ifdef MHD
                  FluxList[Nflux].dB[0] = dir * fluxes.B[0];
                  FluxList[Nflux].dB[1] = dir * fluxes.B[1];
                  FluxList[Nflux].dB[2] = dir * fluxes.B[2];
#ifdef MHD_CT
                  FluxList[Nflux].dA[0] = dir * fluxes.A[0];
                  FluxList[Nflux].dA[1] = dir * fluxes.A[1];
                  FluxList[Nflux].dA[2] = dir * fluxes.A[2];
#endif
#if defined(MHD_POWELL) && !defined(MHD_POWELL_SPLIT)
                  double Velx = state_center->velx + delta_time->velx + vel_face[0];
                  double Vely = state_center->vely + delta_time->vely + vel_face[1];
                  double Velz = state_center->velz + delta_time->velz + vel_face[2];

                  if(All.ComovingIntegrationOn)
                    {
                      Velx += atime * hubble_a * state_center->dx;
                      Vely += atime * hubble_a * state_center->dy;
                      Velz += atime * hubble_a * state_center->dz;
                    }

                  double Bx = state_center->Bx + delta_time->Bx;
                  double By = state_center->By + delta_time->By;
                  double Bz = state_center->Bz + delta_time->Bz;

                  FluxList[Nflux].dB[0] += dir * Velx * state_face.Bx;
                  FluxList[Nflux].dB[1] += dir * Vely * state_face.Bx;
                  FluxList[Nflux].dB[2] += dir * Velz * state_face.Bx;

                  FluxList[Nflux].dP[0] += dir * Bx * state_face.Bx;
                  FluxList[Nflux].dP[1] += dir * By * state_face.Bx;
                  FluxList[Nflux].dP[2] += dir * Bz * state_face.Bx;
#ifndef ISOTHERM_EQS
                  FluxList[Nflux].dEnergy += dir * (Bx * Velx + By * Vely + Bz * Velz) * state_face.Bx;
#endif

                  {
                    double dMomX = dir * Bx * state_face.Bx;
                    double dMomY = dir * By * state_face.Bx;
                    double dMomZ = dir * Bz * state_face.Bx;

                    All.Powell_Momentum[0] += dMomX;
                    All.Powell_Momentum[1] += dMomY;
                    All.Powell_Momentum[2] += dMomZ;

                    double dx = PrimExch[p].Center[0] - 0.5 * All.BoxSize;
                    double dy = PrimExch[p].Center[1] - 0.5 * All.BoxSize;
                    double dz = PrimExch[p].Center[2] - 0.5 * All.BoxSize;

                    All.Powell_Angular_Momentum[0] += dy * dMomZ - dz * dMomY;
                    All.Powell_Angular_Momentum[1] += dz * dMomX - dx * dMomZ;
                    All.Powell_Angular_Momentum[2] += dx * dMomY - dy * dMomX;
                    All.Powell_Energy += dir * (Bx * Velx + By * Vely + Bz * Velz) * state_face.Bx;
                  }
#endif
#ifdef MHD_THERMAL_ENERGY_SWITCH
                  double ptot = PrimExch[p].Pressure + 0.5 * (PrimExch[p].B[0] * PrimExch[p].B[0] + PrimExch[p].B[1] * PrimExch[p].B[1] + PrimExch[p].B[2] * PrimExch[p].B[2]);
#ifdef COSMIC_RAYS
                  ptot += PrimExch[p].CR_Pressure;
#endif
                  FluxList[Nflux].dEtherm = dir * ptot * face_velx;
#endif
#endif

#ifdef COSMIC_RAYS
                  FluxList[Nflux].dEnergy -= dir * PrimExch[p].CR_Pressure * face_velx * All.cf_atime;
                  FluxList[Nflux].dCR_Energy = dir * PrimExch[p].CR_Pressure * face_velx;
#endif

#ifdef MAXSCALARS
                  for(m = 0; m < N_Scalar; m++)
                    FluxList[Nflux].dConservedScalars[m] = dir * fluxes.scalars[m];
#endif

#ifdef GFM_CHEMTAGS
                  for(m = 0; m < GFM_N_CHEM_TAGS; m++)
                    FluxList[Nflux].dMassMetalsChemTags[m] = dir * fluxes.chemtags[m];
#endif

#ifdef FLD
                  FluxList[Nflux].dFLD = dir * fluxes.dFLD;
#endif

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
                  FluxList[Nflux].dEntropy = dir * fluxes.entropy;
#endif
#ifdef TRACER_FIELD
                  FluxList[Nflux].dConservedTracer = dir * fluxes.tracer;
#endif
#ifdef OUTPUT_CELL_SPIN
                  double dx, dy, dz;
                  {
                    double xtmp, ytmp, ztmp;
                    dx = NEAREST_X(geom.cx - PrimExch[p].Center[0] - PrimExch[p].CenterOffset[0]);
                    dy = NEAREST_Y(geom.cy - PrimExch[p].Center[1] - PrimExch[p].CenterOffset[1]);
                    dz = NEAREST_Z(geom.cz - PrimExch[p].Center[2] - PrimExch[p].CenterOffset[2]);
                  }

                  FluxList[Nflux].dSpin[0] = dir * (dy * fluxes.momentum[2] - dz * fluxes.momentum[1]);
                  FluxList[Nflux].dSpin[1] = dir * (dz * fluxes.momentum[0] - dx * fluxes.momentum[2]);
                  FluxList[Nflux].dSpin[2] = dir * (dx * fluxes.momentum[1] - dy * fluxes.momentum[0]);

                  FluxList[Nflux].dCenterOffsetMass[0] = dir * (dx * fluxes.mass);
                  FluxList[Nflux].dCenterOffsetMass[1] = dir * (dy * fluxes.mass);
                  FluxList[Nflux].dCenterOffsetMass[2] = dir * (dz * fluxes.mass);
#endif
#if defined(VISCOSITY) || defined(THERMAL_CONDUCTION) || defined(TRACER_DIFFUSION)
                  FluxList[Nflux].dP[0] -= dir * diffusionfluxes.momentum[0];
                  FluxList[Nflux].dP[1] -= dir * diffusionfluxes.momentum[1];
                  FluxList[Nflux].dP[2] -= dir * diffusionfluxes.momentum[2];
#ifndef ISOTHERM_EQS
                  FluxList[Nflux].dEnergy -= dir * diffusionfluxes.energy;
#endif
#if defined(TRACER_FIELD) && defined(TRACER_DIFFUSION)
                  FluxList[Nflux].dConservedTracer -= dir * diffusionfluxes.tracer;
#endif
#endif
#ifdef TGCHEM
                  for(m = 0; m < TGCHEM_NUM_ABUNDANCES; m++)
                    FluxList[Nflux].Abund[m] = dir * fluxes.pcabund[m];
#endif
                  Nflux++;
                }
            }
        }
    }
  /* end of big loop over all faces */

  TIMER_STOPSTART(CPU_FLUXES, CPU_FLUXES_COMM);

  /* now exchange the flux-list and apply it when needed */
  apply_flux_list();

  TIMER_STOPSTART(CPU_FLUXES_COMM, CPU_FLUXES);

#ifdef TRACER_MC
  myfree(CellTracerMasses);
  finish_MC_tracer();
#endif

#ifdef TGCHEM
  myfree(AbundFluxes);
#endif

  myfree(FluxList);

  double in[2] = {count, count_reduced},  out[2];
  MPI_Reduce(in, out, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(ThisTask == 0)
    {
      tot_count = out[0];
      tot_count_reduced = out[1];

      printf("FLUX: exchanged fluxes over %g faces, with %g reduced (fraction %g), cumulative fraction %g\n",
             tot_count, tot_count_reduced, tot_count_reduced / (tot_count + 1.0e-30), All.TotCountReducedFluxes / (All.TotCountFluxes + 1.0e-30));
      All.TotCountReducedFluxes += tot_count_reduced;
      All.TotCountFluxes += tot_count;
    }

  fvs_evaluate_statistics(&stat);



#ifdef MESHRELAX
  for(i = 0; i < NumGas; i++)
    {
      if(P[i].Mass < 0)
        {
          terminate("negative mass reached for cell=%d mass=%g", P[i].ID, P[i].Mass);

          P[i].Mass = 0;
          SphP[i].Energy = 0;
          SphP[i].Momentum[0] = 0;
          SphP[i].Momentum[1] = 0;
          SphP[i].Momentum[2] = 0;
#ifdef USE_ENTROPY_FOR_COLD_FLOWS
          SphP[i].Entropy = 0;
#endif
        }
    }
#endif

#ifdef GODUNOV_STATS
  endrun();
#endif

#ifdef AXISYMMETRY
  apply_axisymmetric_source_terms(T);
#endif

#ifdef ONEDIMS_SPHERICAL
  apply_spherical_source_terms();
#endif
      
#ifdef GENERAL_RELATIVITY
  apply_general_relativity_source_terms();
#endif      
      
#ifdef TGCHEM
  for(i = 0; i < NumGas; i++)
    for(j = 0; j < TGCHEM_NUM_ABUNDANCES; j++)
      if(AbundFluxes[TGCHEM_NUM_ABUNDANCES * i + j])
        {
          SphP[i].Abund[j] = (SphP[i].Abund[j] * SphP[i].OldMass + AbundFluxes[TGCHEM_NUM_ABUNDANCES * i + j]) / P[i].Mass;

          SphP[i].Abund[j] = dmin(TGCD.AbMax[j], dmax(0, SphP[i].Abund[j]));
        }
#endif

#if defined(MHD_POWELL) && defined(VERBOSE) && !defined(MHD_POWELL_SPLIT)
  double Powell_Momentum[3];
  double Powell_Angular_Momentum[3];
  double Powell_Energy;

  MPI_Reduce(All.Powell_Momentum, Powell_Momentum, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(All.Powell_Angular_Momentum, Powell_Angular_Momentum, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&All.Powell_Energy, &Powell_Energy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    printf("MHD_POWELL: Total ST contribution: Mom=%g,%g,%g   AngMom=%g,%g,%g   Energy=%g\n",
           Powell_Momentum[0], Powell_Momentum[1], Powell_Momentum[2], Powell_Angular_Momentum[0], Powell_Angular_Momentum[1], Powell_Angular_Momentum[2], Powell_Energy);
#endif

#ifdef GFM_DUST_CAP
  gfm_cap_dust();
#endif

#ifdef VORONOI_BACKUP_RESTORE_FACE_AREAS
  restore_face_areas(T);
#endif

#ifdef COSMIC_RAYS_EXTRA_DIAGNOSTICS
  double FinalCREnergy = 0;
  for(int i=0; i<NumGas; i++)
    if(P[i].Mass != 0 && P[i].ID != 0 && P[i].Type == 0)
      FinalCREnergy += SphP[i].CR_Energy;
  
  double dCREnergy = FinalCREnergy - InitialCREnergy;
  All.TotalCREnergyChangeAdiabatic += dCREnergy;
#endif

  TIMER_STOP(CPU_FLUXES);
}



#ifdef VORONOI_BACKUP_RESTORE_FACE_AREAS
void backup_face_areas(tessellation * T)
{
  for(int i = 0; i < T->Nvf; i++)
    T->VF[i].area_backup = T->VF[i].area;
}

void restore_face_areas(tessellation * T)
{
  for(int i = 0; i < T->Nvf; i++)
    T->VF[i].area = T->VF[i].area_backup;
}
#endif


int face_get_state(tessellation * T, int p, int i, struct state *st)
{
  int particle;
#if defined(MAXSCALARS) || defined(TGCHEM)
  int j;
#endif
  double aBegin;
#ifdef MHD_CT
  double aBBegin;
#endif

  point *DP = T->DP;
  face *VF = T->VF;

#ifdef AXISYMMETRY
  if(p == -100)
    {
      particle = i;
      p = 0;
    }
  else
#endif
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

#ifdef ACTIVE_CELL_SPIN
  {
    double xtmp, ytmp, ztmp;
    if(DP[p].task == ThisTask)
      {
        st->dxL = NEAREST_X(VF[i].cx - SphP[particle].Center[0] - SphP[particle].CenterOffset[0]);
        st->dyL = NEAREST_Y(VF[i].cy - SphP[particle].Center[1] - SphP[particle].CenterOffset[1]);
        st->dzL = NEAREST_Z(VF[i].cz - SphP[particle].Center[2] - SphP[particle].CenterOffset[2]);
      }
    else
      {
        st->dxL = NEAREST_X(VF[i].cx - PrimExch[particle].Center[0] - PrimExch[particle].CenterOffset[0]);
        st->dyL = NEAREST_Y(VF[i].cy - PrimExch[particle].Center[1] - PrimExch[particle].CenterOffset[1]);
        st->dzL = NEAREST_Z(VF[i].cz - PrimExch[particle].Center[2] - PrimExch[particle].CenterOffset[2]);
      }
  }
#endif

#ifdef ONEDIMS_SPHERICAL
  if(DP[p].task == ThisTask)
    st->radius = SphP[particle].Center[0];
  else
    st->radius = PrimExch[particle].Center[0];
#endif

  if(DP[p].task == ThisTask)
    {
      st->velGas[0] = P[particle].Vel[0];
      st->velGas[1] = P[particle].Vel[1];
      st->velGas[2] = P[particle].Vel[2];

      st->velVertex[0] = SphP[particle].VelVertex[0];
      st->velVertex[1] = SphP[particle].VelVertex[1];
      st->velVertex[2] = SphP[particle].VelVertex[2];

      st->rho = SphP[particle].Density;
#if defined(SPECIAL_RELATIVITY) || defined(GENERAL_RELATIVITY)
      st->utherm = SphP[particle].Utherm;
#endif
      st->press = SphP[particle].Pressure;

#ifdef COSMIC_RAYS
      st->crPressure = SphP[particle].CR_Pressure;
#endif
      
#ifdef MRT
      st->nHI = SphP[particle].HI ;
      st->nHII =  SphP[particle].HII ;
      st->ne = SphP[particle].Ne ;
#ifdef MRT_INCLUDE_HE
      st->nHeI = SphP[particle].HeI ; 
      st->nHeII = SphP[particle].HeII ;
      st->nHeIII = SphP[particle].HeIII ;
#endif
#endif

      
#ifdef USE_ENTROPY_FOR_COLD_FLOWS
      st->A = SphP[particle].A;
#endif

#ifdef TVD_SLOPE_LIMITER
      st->gradul = &SphP[particle].GradUl;
#endif
      st->grad = &SphP[particle].Grad;

      st->timeBin = P[particle].TimeBinHydro;

      st->volume = SphP[particle].Volume;

#ifdef TGCHEM
      for(j = 0; j < TGCHEM_NUM_ABUNDANCES; j++)
        st->pcabund[j] = SphP[particle].Abund[j];

      st->gamma = SphP[particle].Gamma;
#endif

#ifdef MHD
      st->Bx = SphP[particle].B[0];
      st->By = SphP[particle].B[1];
      st->Bz = SphP[particle].B[2];
#ifdef MHD_CT
      st->Ax = SphP[particle].A[0];
      st->Ay = SphP[particle].A[1];
      st->Az = SphP[particle].A[2];
#endif
#ifdef MHD_POWELL
      st->divB = SphP[particle].DivB;
#endif
#ifdef NON_IDEAL_MHD
#if defined(OHMIC_DIFFUSION) || defined(AMBIPOLAR_DIFFUSION) || defined(IMPLICIT_OHMIC_DIFFUSION)
      st->CurlB[0] = SphP[particle].CurlB[0];
      st->CurlB[1] = SphP[particle].CurlB[1];
      st->CurlB[2] = SphP[particle].CurlB[2];
#endif
#if defined(OHMIC_DIFFUSION)
      double xtmp, ytmp, ztmp;
      st->dcp[0] = NEAREST_X(DP[p].x - SphP[particle].Center[0]);
      st->dcp[1] = NEAREST_Y(DP[p].y - SphP[particle].Center[1]);
      st->dcp[2] = NEAREST_Z(DP[p].z - SphP[particle].Center[2]);
#endif
#endif
#endif

#ifdef VARIABLE_GAMMA
      st->gammaE = SphP[particle].GammaE;
      st->gammaC = SphP[particle].GammaC;
#endif

#ifdef MAXSCALARS
      for(j = 0; j < N_Scalar; j++)
        st->scalars[j] = *(MyFloat *) (((char *) (&SphP[particle])) + scalar_elements[j].offset);
#endif

#ifdef GFM_CHEMTAGS
      for(j = 0; j < GFM_N_CHEM_TAGS; j++)
        st->chemtagsfraction[j] = SphP[particle].MassMetalsChemTagsFraction[j];
#endif

#ifdef FLD
      st->n_gamma = SphP[particle].n_gamma;
      st->R2 = SphP[particle].R2;
#endif

#ifdef TRACER_FIELD
      st->tracer = SphP[particle].Tracer;
#endif

#ifdef SECOND_DERIVATIVES
      st->hessian = &SphP[particle].Hessian;
#endif

#ifdef LOCALLY_ISOTHERM_DISK
      double dx, dy, r;
      dx = P[particle].Pos[0] - boxHalf_X;
      dy = P[particle].Pos[1] - boxHalf_Y;

      r = sqrt(dx * dx + dy * dy);

      st->localSoundSpeed = get_isotherm_disk_sound_speed(particle);
      st->inDiskFlag = get_isotherm_disk_flag(particle);
#endif

#ifdef ACTIVE_CELL_SPIN
      st->wx = SphP[particle].Omega[0];
      st->wy = SphP[particle].Omega[1];
      st->wz = SphP[particle].Omega[2];
#endif

      aBegin = SphP[particle].TimeLastPrimUpdate;

#ifdef MHD_CT
      aBBegin = SphP[particle].TimeLastBUpdate;
#endif

      st->oldmass = SphP[particle].OldMass;
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

      st->rho = PrimExch[particle].Density;

      st->press = PrimExch[particle].Pressure;
#if defined(SPECIAL_RELATIVITY) || defined(GENERAL_RELATIVITY)
      st->utherm = PrimExch[particle].Utherm;
#endif

#ifdef COSMIC_RAYS
      st->crPressure = PrimExch[particle].CR_Pressure;
#endif
      
#ifdef MRT
      st->nHI = PrimExch[particle].HI ;
      st->nHII =  PrimExch[particle].HII ;
      st->ne = PrimExch[particle].Ne ;
#ifdef MRT_INCLUDE_HE
      st->nHeI = PrimExch[particle].HeI ; 
      st->nHeII = PrimExch[particle].HeII ;
      st->nHeIII = PrimExch[particle].HeIII ;
#endif
#endif
      
#ifdef USE_ENTROPY_FOR_COLD_FLOWS
      st->A = PrimExch[particle].A;
#endif

#ifdef TVD_SLOPE_LIMITER
      st->gradul = &GradExchUl[particle];
#endif
      st->grad = &GradExch[particle];

      st->timeBin = PrimExch[particle].TimeBinHydro;    /* This is the hydro timestep */

      st->volume = PrimExch[particle].Volume;

#ifdef TGCHEM
      for(j = 0; j < TGCHEM_NUM_ABUNDANCES; j++)
        st->pcabund[j] = PrimExch[particle].Abund[j];

      st->gamma = PrimExch[particle].Gamma;
#endif

#ifdef MHD
      st->Bx = PrimExch[particle].B[0];
      st->By = PrimExch[particle].B[1];
      st->Bz = PrimExch[particle].B[2];
#ifdef MHD_CT
      st->Ax = PrimExch[particle].A[0];
      st->Ay = PrimExch[particle].A[1];
      st->Az = PrimExch[particle].A[2];
#endif
#ifdef MHD_POWELL
      st->divB = PrimExch[particle].DivB;
#endif
#ifdef NON_IDEAL_MHD
#if defined(OHMIC_DIFFUSION) || defined(AMBIPOLAR_DIFFUSION) || defined(IMPLICIT_OHMIC_DIFFUSION)
      st->CurlB[0] = PrimExch[particle].CurlB[0];
      st->CurlB[1] = PrimExch[particle].CurlB[1];
      st->CurlB[2] = PrimExch[particle].CurlB[2];
#endif
#if defined(OHMIC_DIFFUSION)
      double xtmp, ytmp, ztmp;
      st->dcp[0] = NEAREST_X(DP[p].x - PrimExch[particle].Center[0]);
      st->dcp[1] = NEAREST_Y(DP[p].y - PrimExch[particle].Center[1]);
      st->dcp[2] = NEAREST_Z(DP[p].z - PrimExch[particle].Center[2]);
#endif
#endif
#endif

#ifdef VARIABLE_GAMMA
      st->gammaE = PrimExch[particle].GammaE;
      st->gammaC = PrimExch[particle].GammaC;
#endif

#ifdef MAXSCALARS
      for(j = 0; j < N_Scalar; j++)
        st->scalars[j] = PrimExch[particle].Scalars[j];
#endif

#ifdef GFM_CHEMTAGS
      for(j = 0; j < GFM_N_CHEM_TAGS; j++)
        st->chemtagsfraction[j] = PrimExch[particle].MassMetalsChemTagsFraction[j];
#endif

#ifdef TRACER_FIELD
      st->tracer = PrimExch[particle].Tracer;
#endif

#ifdef SECOND_DERIVATIVES
      st->hessian = &HessianExch[particle];
#endif

#ifdef FLD
      st->n_gamma = PrimExch[particle].n_gamma;
      st->R2 = PrimExch[particle].R2;
#endif

#ifdef LOCALLY_ISOTHERM_DISK
      double dx, dy, r;
      dx = PrimExch[particle].Center[0] - boxHalf_X;
      dy = PrimExch[particle].Center[1] - boxHalf_Y;

      r = sqrt(dx * dx + dy * dy);

      st->localSoundSpeed = All.AspectRatio * sqrt(All.G * All.CentralMass / r);
      if(r > All.inner_radius && r < All.outer_radius)
        st->inDiskFlag = 1;
      else
        st->inDiskFlag = 0;
#endif

#ifdef ACTIVE_CELL_SPIN
      st->wx = PrimExch[particle].Omega[0];
      st->wy = PrimExch[particle].Omega[1];
      st->wz = PrimExch[particle].Omega[2];
#endif

      aBegin = PrimExch[particle].TimeLastPrimUpdate;

#ifdef MHD_CT
      aBBegin = PrimExch[particle].TimeLastBUpdate;
#endif

      st->oldmass = PrimExch[particle].OldMass;
      st->surfacearea = PrimExch[particle].SurfaceArea;
      st->activearea = PrimExch[particle].ActiveArea;
      st->csnd = PrimExch[particle].Csnd;
      st->ID = DP[p].ID;
    }

  st->dtExtrapolation = All.Time - aBegin;

#ifdef MHD_CT
  st->dtBExtrapolation = All.Time - aBBegin;
#endif

  /* check for reflecting or outflowing boundaries */
  face_boundary_check_vertex(T, p, &st->velVertex[0], &st->velVertex[1], &st->velVertex[2]);


#ifdef FLD_TEST_BOUNDARY
#ifdef REFLECTIVE_Y
  if(DP[p].task == ThisTask)
    {
      if((T->DP[p].image_flags & REFL_Y_FLAGS))
        {
          if(P[particle].ID >= FLD_UPPER_BOUNDARY_MINID && P[particle].ID < FLD_UPPER_BOUNDARY_MAXID)
            {
              st->velGas[0] = 0.;
              st->velGas[1] = 0.;
              st->velGas[2] = 0.;

              st->rho = All.fld_density;
              st->press = GAMMA_MINUS1 * All.fld_density * All.fld_u;
            }
        }
    }
#endif
#endif

  return 0;
}

void face_boundary_check_vertex(tessellation * T, int p, MyFloat * velx, MyFloat * vely, MyFloat * velz)
{
  /* check for reflecting or outflowing boundaries */
#if defined(REFLECTIVE_X)
  if((T->DP[p].image_flags & REFL_X_FLAGS))     // && !(DP[p].image_flags & OUTFLOW_X))
    *velx *= -1;
#endif
#if defined(REFLECTIVE_Y)
  if((T->DP[p].image_flags & REFL_Y_FLAGS))     // && !(DP[p].image_flags & OUTFLOW_Y))
    *vely *= -1;
#endif
#if defined(REFLECTIVE_Z)
  if((T->DP[p].image_flags & REFL_Z_FLAGS))     // && !(DP[p].image_flags & OUTFLOW_Z))
    *velz *= -1;
#endif

#ifdef ONEDIMS_SPHERICAL
  if(p == -1)
    *velx *= -1;
#endif
}


void face_boundary_check(point * p, double *velx, double *vely, double *velz)
{
  /* check for reflecting or outflowing boundaries */
#if defined(REFLECTIVE_X)
  if((p->image_flags & REFL_X_FLAGS) && !(p->image_flags & OUTFLOW_X))
    *velx *= -1;
#endif
#if defined(REFLECTIVE_Y)
  if((p->image_flags & REFL_Y_FLAGS) && !(p->image_flags & OUTFLOW_Y))
    *vely *= -1;
#endif
#if defined(REFLECTIVE_Z)
  if((p->image_flags & REFL_Z_FLAGS) && !(p->image_flags & OUTFLOW_Z))
    *velz *= -1;
#endif

#ifdef ONEDIMS_SPHERICAL
  if(p == &Mesh.DP[-1])
    *velx *= -1;
#endif
}


int face_check_responsibility_of_this_task(tessellation * T, int p1, int p2, struct state *st_L, struct state *st_R)
{
  int low_p, high_p;
  struct state *low_state, *high_state;

  point *DP = T->DP;

  if(DP[p1].ID < DP[p2].ID)
    {
      low_p = p1;
      high_p = p2;
      low_state = st_L;
      high_state = st_R;
    }
  else if(DP[p1].ID > DP[p2].ID)
    {
      low_p = p2;
      high_p = p1;
      low_state = st_R;
      high_state = st_L;
    }
  else
    {
      /* equality of the IDs should only occur for reflective boundaries */
      if(DP[p1].task == ThisTask && DP[p1].index < NumGas)
        {
          low_p = p1;
          high_p = p2;
          low_state = st_L;
          high_state = st_R;
        }
      else
        {
          low_p = p2;
          high_p = p1;
          low_state = st_R;
          high_state = st_L;
        }
    }

  if(TimeBinSynchronized[low_state->timeBin])   /* the one with the lower ID is active */
    {
      /* we need to check whether the one with the lower ID is a local particle */
      if(DP[low_p].task == ThisTask && DP[low_p].index < NumGas)
        return 0;
    }
  else if(TimeBinSynchronized[high_state->timeBin])     /* only the side with the higher ID is active */
    {
      /* we need to check whether we hold the one with the higher ID, if yes, we'll do it */
      if(DP[high_p].task == ThisTask && DP[high_p].index < NumGas)
        return 0;
    }

  return -1;                    /* we can skip this face on the local task */
}


double face_timestep(struct state *state_L, struct state *state_R, double *hubble_a, double *atime)
{
  integertime ti_begin_L, ti_begin_R;
  short int timeBin;
  double face_dt;

  /* determine most recent start of the time bins */
  ti_begin_L = (All.Ti_Current >> state_L->timeBin) << state_L->timeBin;
  ti_begin_R = (All.Ti_Current >> state_R->timeBin) << state_R->timeBin;

  /* take the minimum of the two */
  timeBin = state_L->timeBin;
  if(timeBin > state_R->timeBin)
    timeBin = state_R->timeBin;

  /* compute the half-step prediction times */
  state_L->dt_half = (All.Ti_Current + (((integertime) 1) << (timeBin - 1)) - ti_begin_L) * All.Timebase_interval;
  state_R->dt_half = (All.Ti_Current + (((integertime) 1) << (timeBin - 1)) - ti_begin_R) * All.Timebase_interval;

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
#ifdef MHD_CT
      state_L->dtBExtrapolation /= *hubble_a;
      state_L->dtBExtrapolation /= *atime;
      state_R->dtBExtrapolation /= *hubble_a;
      state_R->dtBExtrapolation /= *atime;
#endif
    }

  return face_dt;
}

void state_convert_to_local_frame(struct state *st, double *vel_face, double hubble_a, double atime)
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

#ifdef ACTIVE_CELL_SPIN
void face_calc_spin_velocities(struct state *st)
{
  double fac = 1.0;
  double dvel = sqrt(st->dvx * st->dvx + st->dvy * st->dvy + st->dvz * st->dvz);
  if(dvel > st->csnd)
    fac = st->csnd / dvel;

  //double dist = (st->dx-st->dxL)*(st->dx-st->dxL) + (st->dy-st->dyL)*(st->dy-st->dyL) + (st->dz-st->dzL)*(st->dz-st->dzL);
  //if (dist > 0.8 * get_cell_radius_from_volume( st->volume ))
  //  fac = 0.;

  st->dvx = fac * (st->wy * st->dzL - st->wz * st->dyL);
  st->dvy = fac * (st->wz * st->dxL - st->wx * st->dzL);
  st->dvz = fac * (st->wx * st->dyL - st->wy * st->dxL);
}

void face_add_spin_velocities(struct state *state_L, struct state *state_R)
{
  state_L->velx += state_L->dvx;
  state_R->velx += state_R->dvx;
  state_L->vely += state_L->dvy;
  state_R->vely += state_R->dvy;
  state_L->velz += state_L->dvz;
  state_R->velz += state_R->dvz;
/*
  if (state_L->velx > state_R->velx) {
    if (state_L->velx + state_L->dvx < state_R->velx + state_R->dvx) {
      double av = 0.5 * (state_L->velx + state_R->velx);
      state_L->velx = av;
      state_R->velx = av;
    } else {
      state_L->velx += state_L->dvx;
      state_R->velx += state_R->dvx;
    }
  } if (state_L->velx < state_R->velx) {
    if (state_L->velx + state_L->dvx > state_R->velx + state_R->dvx) {
      double av = 0.5 * (state_L->velx + state_R->velx);
      state_L->velx = av;
      state_R->velx = av;
    } else {
      state_L->velx += state_L->dvx;
      state_R->velx += state_R->dvx;
    }
  }

  if (state_L->vely > state_R->vely) {
    if (state_L->vely + state_L->dvy < state_R->vely + state_R->dvy) {
        double av = 0.5 * (state_L->vely + state_R->vely);
        state_L->vely = av;
        state_R->vely = av;
      } else {
        state_L->vely += state_L->dvy;
        state_R->vely += state_R->dvy;
      }
  } if (state_L->vely < state_R->vely) {
    if (state_L->vely + state_L->dvy > state_R->vely + state_R->dvy) {
      double av = 0.5 * (state_L->vely + state_R->vely);
      state_L->vely = av;
      state_R->vely = av;
    } else {
      state_L->vely += state_L->dvy;
      state_R->vely += state_R->dvy;
    }
  }

  if (state_L->velz > state_R->velz) {
    if (state_L->velz + state_L->dvz < state_R->velz + state_R->dvz) {
      double av = 0.5 * (state_L->velz + state_R->velz);
      state_L->velz = av;
      state_R->velz = av;
    } else {
      state_L->velz += state_L->dvz;
      state_R->velz += state_R->dvz;
    }
  } if (state_L->velz < state_R->velz) {
    if (state_L->velz + state_L->dvz > state_R->velz + state_R->dvz) {
      double av = 0.5 * (state_L->velz + state_R->velz);
      state_L->velz = av;
      state_R->velz = av;
    } else {
      state_L->velz += state_L->dvz;
      state_R->velz += state_R->dvz;
    }
  }
*/
}
#endif

void face_do_time_extrapolation(struct state *delta, struct state *st, double atime)
{
  /* st is the state at the center of the cell */

  /* the code still allows for emtpy cells but we are going to divide by rho, so ... */
  if(st->rho <= 0)
    return;

#if defined (MESHRELAX) || defined (DISABLE_TIME_EXTRAPOLATION)
  /* do not time extrapolation */
  (void) st;
  (void) atime;
  memset(delta, 0, sizeof(struct state));
  return;
#endif

  struct grad_data *grad = st->grad;

#ifndef MUSCL_HANCOCK           /* we use the default Runge Kutta time integration scheme */
  double dt_half = st->dtExtrapolation;
#ifdef MHD_CT
  double dt_Bhalf = st->dtBExtrapolation;
#endif
#else
  double dt_half = st->dt_half;
#ifdef MHD_CT
  double dt_Bhalf = st->dt_half;
#endif
#endif

  if(All.ComovingIntegrationOn)
    dt_half /= atime;

  delta->rho = -dt_half * (st->velx * grad->drho[0] + st->rho * grad->dvel[0][0] + st->vely * grad->drho[1] + st->rho * grad->dvel[1][1] + st->velz * grad->drho[2] + st->rho * grad->dvel[2][2]);

#if defined(SPECIAL_RELATIVITY) || defined(GENERAL_RELATIVITY)
  delta->utherm = -dt_half * (st->velx * grad->dutherm[0] + st->utherm * grad->dvel[0][0] +
                              st->vely * grad->dutherm[1] + st->utherm * grad->dvel[1][1] + st->velz * grad->dutherm[2] + st->utherm * grad->dvel[2][2]);
#endif

  delta->velx = -dt_half * (1.0 / st->rho * grad->dpress[0] + st->velx * grad->dvel[0][0] + st->vely * grad->dvel[0][1] + st->velz * grad->dvel[0][2]);

  delta->vely = -dt_half * (1.0 / st->rho * grad->dpress[1] + st->velx * grad->dvel[1][0] + st->vely * grad->dvel[1][1] + st->velz * grad->dvel[1][2]);

  delta->velz = -dt_half * (1.0 / st->rho * grad->dpress[2] + st->velx * grad->dvel[2][0] + st->vely * grad->dvel[2][1] + st->velz * grad->dvel[2][2]);

#ifdef TGCHEM
  delta->press = -dt_half * (st->gamma * st->press * (grad->dvel[0][0] + grad->dvel[1][1] + grad->dvel[2][2]) + st->velx * grad->dpress[0] + st->vely * grad->dpress[1] + st->velz * grad->dpress[2]);
#else
#ifndef VARIABLE_GAMMA
  delta->press = -dt_half * (GAMMA * st->press * (grad->dvel[0][0] + grad->dvel[1][1] + grad->dvel[2][2]) + st->velx * grad->dpress[0] + st->vely * grad->dpress[1] + st->velz * grad->dpress[2]);
#else
  delta->press = -dt_half * (st->gammaE * st->press * (grad->dvel[0][0] + grad->dvel[1][1] + grad->dvel[2][2]) + st->velx * grad->dpress[0] + st->vely * grad->dpress[1] + st->velz * grad->dpress[2]);
#endif
#endif

#ifdef ONEDIMS_SPHERICAL
  delta->velx += dt_half * 2. * st->press / (st->rho * st->radius);
#endif

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
  delta->A = -dt_half * (st->velx * grad->dA[0] + st->vely * grad->dA[1] + st->velz * grad->dA[2]);
#endif

#ifdef TGCHEM
  delta->gamma = -dt_half * (st->velx * grad->dgamma[0] + st->vely * grad->dgamma[1] + st->velz * grad->dgamma[2]);
#endif

#ifdef VARIABLE_GAMMA
  delta->gammaE = -dt_half * (st->velx * grad->dgammaE[0] + st->vely * grad->dgammaE[1] + st->velz * grad->dgammaE[2]);
  delta->gammaC = -dt_half * (st->velx * grad->dgammaC[0] + st->vely * grad->dgammaC[1] + st->velz * grad->dgammaC[2]);
#endif

#ifdef MHD
  delta->velx += -dt_half * (1.0 / st->rho * (st->By * grad->dB[1][0] + st->Bz * grad->dB[2][0] - st->By * grad->dB[0][1] - st->Bz * grad->dB[0][2]) / atime);

  delta->vely += -dt_half * (1.0 / st->rho * (st->Bx * grad->dB[0][1] + st->Bz * grad->dB[2][1] - st->Bx * grad->dB[1][0] - st->Bz * grad->dB[1][2]) / atime);

  delta->velz += -dt_half * (1.0 / st->rho * (st->Bx * grad->dB[0][2] + st->By * grad->dB[1][2] - st->Bx * grad->dB[2][0] - st->By * grad->dB[2][1]) / atime);

  delta->Bx =
    -dt_half * (-st->velx * grad->dB[1][1] - grad->dvel[0][1] * st->By + st->vely * grad->dB[0][1] +
                grad->dvel[1][1] * st->Bx + st->velz * grad->dB[0][2] + grad->dvel[2][2] * st->Bx - st->velx * grad->dB[2][2] - grad->dvel[0][2] * st->Bz);

  delta->By =
    -dt_half * (+st->velx * grad->dB[1][0] + grad->dvel[0][0] * st->By - st->vely * grad->dB[0][0] -
                grad->dvel[1][0] * st->Bx - st->vely * grad->dB[2][2] - grad->dvel[1][2] * st->Bz + st->velz * grad->dB[1][2] + grad->dvel[2][2] * st->By);

  delta->Bz =
    -dt_half * (-st->velz * grad->dB[0][0] - grad->dvel[2][0] * st->Bx + st->velx * grad->dB[2][0] +
                grad->dvel[0][0] * st->Bz + st->vely * grad->dB[2][1] + grad->dvel[1][1] * st->Bz - st->velz * grad->dB[1][1] - grad->dvel[2][1] * st->By);
#ifdef MHD_CT
  delta->Ax = 0;                //-dt_Bhalf * (-(st->vely * st->Bz - st->velz * st->By)); // atime is there to cancel the factor for spatial gradients
  delta->Ay = 0;                //-dt_Bhalf * ( (st->velx * st->Bz - st->velz * st->Bx));
  delta->Az = 0;                //-dt_Bhalf * (-(st->velx * st->By - st->vely * st->Bx));
#endif
#ifdef MHD_POWELL
  /*
     This seems to hurt more than it helps, possibly the divB estimate is not good enough...
     delta->velx += -dt_half * st->divB * st->Bx / st->rho;
     delta->vely += -dt_half * st->divB * st->By / st->rho;
     delta->velz += -dt_half * st->divB * st->Bz / st->rho;

     delta->Bx += -dt_half * st->divB * st->velx;
     delta->By += -dt_half * st->divB * st->vely;
     delta->Bz += -dt_half * st->divB * st->velz;
   */
  /* energy extrapolation does not change pressure... */
#endif
#endif

#if defined(MAXSCALARS)
  int k;
  for(k = 0; k < N_Scalar; k++)
    {
#if defined(GFM_DISCARD_ENRICHMENT_GRADIENTS)
      delta->scalars[k] = 0;
#else
      delta->scalars[k] = -dt_half * (st->velx * grad->dscalars[k][0] + st->vely * grad->dscalars[k][1] + st->velz * grad->dscalars[k][2]);
#endif
    }
#endif

#ifdef TRACER_FIELD
  delta->tracer = -dt_half * (st->velx * grad->dtracer[0] + st->vely * grad->dtracer[1] + st->velz * grad->dtracer[2]);
#endif

#ifdef COSMIC_RAYS
  delta->crPressure = -dt_half * (st->crPressure * (grad->dvel[0][0] + grad->dvel[1][1] + grad->dvel[2][2]) +
                                  st->velx * grad->dcrPressure[0] + st->vely * grad->dcrPressure[1] + st->velz * grad->dcrPressure[2]);
#endif

#ifdef LOCALLY_ISOTHERM_DISK
  delta->localSoundSpeed = 0.0;
#endif
  
}

void face_do_spatial_extrapolation(struct state *delta, struct state *st, struct state *st_other)
{
#ifdef DISABLE_SPATIAL_RECONSTRUCTION
  memset(delta, 0, sizeof(struct state));
  return;
#endif

#ifdef NO_RECONSTRUCTION_AT_STRONG_SHOCKS
  if(dmax(st->press, st_other->press) > 100. * dmin(st->press, st_other->press))
    {
      memset(delta, 0, sizeof(struct state));
      return;
    }
#endif

#ifdef TVD_SLOPE_LIMITER
  struct grad_data *grad = st->gradul;
#else
  struct grad_data *grad = st->grad;
#endif

  double dx[3];
  dx[0] = st->dx;
  dx[1] = st->dy;
  dx[2] = st->dz;
 
  double r[3];
  r[0] = -st_other->dx + st->dx;
  r[1] = -st_other->dy + st->dy;
  r[2] = -st_other->dz + st->dz;
  
  face_do_spatial_extrapolation_single_quantity(&delta->rho, st->rho, st_other->rho, grad->drho, dx, r);

  face_do_spatial_extrapolation_single_quantity(&delta->velx, st->velx, st_other->velx, grad->dvel[0], dx, r);
  face_do_spatial_extrapolation_single_quantity(&delta->vely, st->vely, st_other->vely, grad->dvel[1], dx, r);
  face_do_spatial_extrapolation_single_quantity(&delta->velz, st->velz, st_other->velz, grad->dvel[2], dx, r);

  face_do_spatial_extrapolation_single_quantity(&delta->press, st->press, st_other->press, grad->dpress, dx, r);

#if defined(SPECIAL_RELATIVITY) || defined(GENERAL_RELATIVITY)
  face_do_spatial_extrapolation_single_quantity(&delta->utherm, st->utherm, st_other->utherm, grad->dutherm, dx, r);
#endif

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
  face_do_spatial_extrapolation_single_quantity(&delta->A, st->A, st_other->A, grad->dA, dx, r);
#endif

#ifdef TGCHEM
  face_do_spatial_extrapolation_single_quantity(&delta->gamma, st->gamma, st_other->gamma, grad->dgamma, dx, r);
#endif

#ifdef VARIABLE_GAMMA
  face_do_spatial_extrapolation_single_quantity(&delta->gammaE, st->gammaE, st_other->gammaE, grad->dgammaE, dx, r);
  face_do_spatial_extrapolation_single_quantity(&delta->gammaC, st->gammaC, st_other->gammaC, grad->dgammaC, dx, r);
#endif

#ifdef MHD
  face_do_spatial_extrapolation_single_quantity(&delta->Bx, st->Bx, st_other->Bx, grad->dB[0], dx, r);
  face_do_spatial_extrapolation_single_quantity(&delta->By, st->By, st_other->By, grad->dB[1], dx, r);
  face_do_spatial_extrapolation_single_quantity(&delta->Bz, st->Bz, st_other->Bz, grad->dB[2], dx, r);
#endif

#ifdef MHD_CT
  face_do_spatial_extrapolation_single_quantity(&delta->Ax, st->Ax, st_other->Ax, grad->dA[0], dx, r);
  face_do_spatial_extrapolation_single_quantity(&delta->Ay, st->Ay, st_other->Ay, grad->dA[1], dx, r);
  face_do_spatial_extrapolation_single_quantity(&delta->Az, st->Az, st_other->Az, grad->dA[2], dx, r);
#endif
  
  
#ifdef MAXSCALARS
  int k;
  for(k = 0; k < N_Scalar; k++)
    {
#if defined(GFM_DISCARD_ENRICHMENT_GRADIENTS)
      delta->scalars[k] = 0;
#else
      face_do_spatial_extrapolation_single_quantity(&delta->scalars[k], st->scalars[k], st_other->scalars[k], grad->dscalars[k], dx, r);
#endif
    }
#endif

#ifdef TRACER_FIELD
  face_do_spatial_extrapolation_single_quantity(&delta->tracer, st->tracer, st_other->tracer, grad->dtracer, dx, r);
#endif

#ifdef LOCALLY_ISOTHERM_DISK
  /*In this case the sound speed is a externally-imposed function of space */
  double csnd_grad[3];
  int j;
  for(j = 0; j < 3; j++)
    csnd_grad[j] = 0.5 / st->localSoundSpeed / st->rho * (grad->dpress[j] - st->localSoundSpeed * st->localSoundSpeed * grad->drho[j]);

  face_do_spatial_extrapolation_single_quantity(&delta->localSoundSpeed, st->localSoundSpeed, st_other->localSoundSpeed, csnd_grad, dx, r);
#endif

#ifdef COSMIC_RAYS
  face_do_spatial_extrapolation_single_quantity(&delta->crPressure, st->crPressure, st_other->crPressure, grad->dcrPressure, dx, r);
#endif
}

void face_do_spatial_extrapolation_single_quantity(double *delta, double st, double st_other, MySingle * grad, double *dx, double *r)
{
#ifndef TVD_SLOPE_LIMITER
  (void) st;
  (void) st_other;
  (void) r;
  *delta = grad[0] * dx[0] + grad[1] * dx[1] + grad[2] * dx[2];
#else
  (void) dx;

  if(st_other == st)
    {
      *delta = 0;
      return;
    }
    
// for GR version try double rf = (grad[0] * r[0] + grad[1] * r[1] + grad[2] * r[2]) / (st_other - st); -> possibly more noisy
  double rf = 2. * (grad[0] * r[0] + grad[1] * r[1] + grad[2] * r[2]) / (st_other - st) - 1.0;

#if !defined(TVD_SLOPE_LIMITER_VANLEER) && !defined(TVD_SLOPE_LIMITER_SUPERBEE) && !defined(TVD_SLOPE_LIMITER_ALBADA) && !defined(TVD_SLOPE_LIMITER_MINBEE) && !defined(TVD_SLOPE_LIMITER_MINMOD) && !defined(TVD_SLOPE_LIMITER_MC)
#error "Choose a slope limiter"
#endif

  double psi;

#ifdef TVD_SLOPE_LIMITER_MINMOD
  if(rf <= 0.)
    {
      psi = 0.;
    }
  else if (rf <= 1.)
    {
      psi = rf;
    }
  else
    {
      psi = 1.;
    }
#endif

#ifdef TVD_SLOPE_LIMITER_VANLEER
  if(rf <= 0.)
    {
      psi = 0.;
    }
  else
    {
      double psi_R = 2. / (1. + rf);
      psi = dmin(2. * rf / (1. + rf), psi_R);
    }
#endif

#ifdef TVD_SLOPE_LIMITER_MC
  if(rf <= 0.)
    {
      psi = 0.;
    }
  else
    {
      double t1 = 2.0 * rf;
      double t2 = 0.5*(1.0+rf);
      psi = dmin(dmin(t1,t2),2);
    }
/*  if(rf <= 0.)   incorrect???
    {
      psi = 0.;
    }
  else if (3.*rf <= 1.)
    {
      psi = 2. * rf;
    }
  else
    {
      double psi_R = 2. / (1. + rf);
      psi = dmin( dmin( 0.5*(1. + rf), psi_R ), 2. );
    }*/   
#endif

#ifdef TVD_SLOPE_LIMITER_SUPERBEE
  if(rf <= 0.)
    {
      psi = 0.;
    }
  else if(rf <= 0.5)
    {
      psi = 2. * rf;
    }
  else if(rf <= 1.)
    {
      psi = 1.;
    }
  else
    {
      double psi_R = 2. / (1. + rf);
      psi = dmin(dmin(rf, psi_R), 2.);
    }
#endif

#ifdef TVD_SLOPE_LIMITER_ALBADA
  if(rf <= 0.)
    {
      psi = 0.;
    }
  else
    {
      double psi_R = 2. / (1. + rf);
      psi = dmin((rf * (1. + rf)) / (1. + rf * rf), psi_R);
    }
#endif

#ifdef TVD_SLOPE_LIMITER_MINBEE
  if(rf <= 0.)
    {
      psi = 0.;
    }
  else if(rf <= 1.)
    {
      psi = rf;
    }
  else
    {
      double psi_R = 2. / (1. + rf);
      psi = dmin(1., psi_R);
    }
#endif

  /* limiter is applied */
  *delta = 0.5 * psi * (st_other - st);
#endif
}

void face_add_extrapolations(struct state *st_face, struct state *delta_time, struct state *delta_space, struct fvs_stat *stat)
{
  stat->count_disable_extrapolation += 1;

#if defined(BOUNDARY_INFLOWOUTFLOW_MINID) && defined(BOUNDARY_INFLOWOUTFLOW_MAXID)
  if(st_face->ID >= BOUNDARY_INFLOWOUTFLOW_MINID && st_face->ID < BOUNDARY_INFLOWOUTFLOW_MAXID)
    return;
#endif

#if defined(BOUNDARY_REFL_SOLIDSIDE_MINID) && defined(BOUNDARY_REFL_SOLIDSIDE_MAXID)
  if(st_face->ID >= BOUNDARY_REFL_SOLIDSIDE_MINID && st_face->ID < BOUNDARY_REFL_SOLIDSIDE_MAXID)
    return;
#endif

  if(st_face->rho <= 0)
    return;

  if(st_face->rho + delta_time->rho + delta_space->rho < 0 || st_face->press + delta_time->press + delta_space->press < 0)
    return;

#if defined (USE_ENTROPY_FOR_COLD_FLOWS)
  if(st_face->A + delta_time->A + delta_space->A < 0)
    return;
#endif

#if defined(SPECIAL_RELATIVITY) || defined(GENERAL_RELATIVITY)
  if(st_face->utherm + delta_time->utherm + delta_space->utherm < 0)
    return;
#endif

  stat->count_disable_extrapolation -= 1;

#if !defined(MESHRELAX) && !defined(DISABLE_TIME_EXTRAPOLATION) 
#ifdef MRT
  delta_time->flag = 1 ;
 #endif
  face_add_extrapolation(st_face, delta_time, stat);
#endif

#if !defined(DISABLE_SPATIAL_EXTRAPOLATION)
#ifdef MRT
  delta_space->flag = 0 ;
#endif
  face_add_extrapolation(st_face, delta_space, stat);
#endif

}

 void face_add_extrapolation(struct state *st_face, struct state *delta, struct fvs_stat *stat)
{
  st_face->rho += delta->rho;
  st_face->velx += delta->velx;
  st_face->vely += delta->vely;
  st_face->velz += delta->velz;
  st_face->press += delta->press;

#ifdef COSMIC_RAYS
  if((st_face->crPressure + delta->crPressure) >= 0.0)
    st_face->crPressure += delta->crPressure;
  else
    stat->count_CR_limiter += 1;
#endif

#if defined(SPECIAL_RELATIVITY) || defined(GENERAL_RELATIVITY)
  st_face->utherm += delta->utherm;
#endif

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
  st_face->A += delta->A;
#endif

#ifdef TGCHEM
  st_face->gamma += delta->gamma;
#endif

#ifdef VARIABLE_GAMMA
  st_face->gammaE += delta->gammaE;
  st_face->gammaC += delta->gammaC;
#endif

#ifdef MHD
#ifndef ONEDIMS
  /* in one dimension, Bx has to be constant! */
  st_face->Bx += delta->Bx;
#endif
  st_face->By += delta->By;
  st_face->Bz += delta->Bz;

#ifdef MHD_CT
  st_face->Ax += delta->Ax;
  st_face->Ay += delta->Ay;
  st_face->Az += delta->Az;
#endif

#endif

#ifdef MAXSCALARS
  int k;
  for(k = 0; k < N_Scalar; k++)
    st_face->scalars[k] += delta->scalars[k];
#endif

#ifdef TRACER_FIELD
  st_face->tracer += delta->tracer;
#endif

#ifdef LOCALLY_ISOTHERM_DISK
  st_face->localSoundSpeed += delta->localSoundSpeed;
#endif

}

void face_add_extrapolation_with_check(struct state *st_face, struct state *delta, struct fvs_stat *stat)
{
  stat->count_disable_extrapolation += 1;

#if defined(BOUNDARY_INFLOWOUTFLOW_MINID) && defined(BOUNDARY_INFLOWOUTFLOW_MAXID)
  if(st_face->ID >= BOUNDARY_INFLOWOUTFLOW_MINID && st_face->ID < BOUNDARY_INFLOWOUTFLOW_MAXID)
    return;
#endif

#if defined(BOUNDARY_REFL_SOLIDSIDE_MINID) && defined(BOUNDARY_REFL_SOLIDSIDE_MAXID)
  if(st_face->ID >= BOUNDARY_REFL_SOLIDSIDE_MINID && st_face->ID < BOUNDARY_REFL_SOLIDSIDE_MAXID)
    return;
#endif

  if(st_face->rho <= 0)
    return;

  if(st_face->rho + delta->rho < 0 || st_face->press + delta->press < 0)
    return;

  stat->count_disable_extrapolation -= 1;

  face_add_extrapolation(st_face, delta, stat);
}


void face_turn_velocities(struct state *st, struct geometry *geom)
{
  double velx, vely, velz;

  velx = st->velx;
  vely = st->vely;
  velz = st->velz;

  st->velx = velx * geom->nx + vely * geom->ny + velz * geom->nz;
  st->vely = velx * geom->mx + vely * geom->my + velz * geom->mz;
  st->velz = velx * geom->px + vely * geom->py + velz * geom->pz;

#ifdef ACTIVE_CELL_SPIN
  double dvx, dvy, dvz;

  dvx = st->dvx;
  dvy = st->dvy;
  dvz = st->dvz;

  st->dvx = dvx * geom->nx + dvy * geom->ny + dvz * geom->nz;
  st->dvy = dvx * geom->mx + dvy * geom->my + dvz * geom->mz;
  st->dvz = dvx * geom->px + dvy * geom->py + dvz * geom->pz;
#endif

#ifdef MHD
  double Bx, By, Bz;

  Bx = st->Bx;
  By = st->By;
  Bz = st->Bz;

  st->Bx = Bx * geom->nx + By * geom->ny + Bz * geom->nz;
  st->By = Bx * geom->mx + By * geom->my + Bz * geom->mz;
  st->Bz = Bx * geom->px + By * geom->py + Bz * geom->pz;
#endif


}



void solve_advection(struct state *st_L, struct state *st_R, struct state_face *st_face, struct geometry *geom, double *vel_face)
{
  double ev = vel_face[0] * geom->nx + vel_face[1] * geom->ny + vel_face[2] * geom->nz;

  if(ev < 0)
    {
      st_face->rho = st_L->rho;
      st_face->velx = st_L->velx;
      st_face->vely = st_L->vely;
      st_face->velz = st_L->velz;
      st_face->press = st_L->press;
#if defined(SPECIAL_RELATIVITY) || defined(GENERAL_RELATIVITY)
      st_face->utherm = st_L->utherm;
#endif
    }
  else
    {
      st_face->rho = st_R->rho;
      st_face->velx = st_R->velx;
      st_face->vely = st_R->vely;
      st_face->velz = st_R->velz;
      st_face->press = st_R->press;
#if defined(SPECIAL_RELATIVITY) || defined(GENERAL_RELATIVITY)
      st_face->utherm = st_R->utherm;
#endif
    }
}



void face_turnback_velocities(struct state_face *st_face, struct geometry *geom)
{
  double velx, vely, velz;

  velx = st_face->velx;
  vely = st_face->vely;
  velz = st_face->velz;

  st_face->velx = velx * geom->nx + vely * geom->mx + velz * geom->px;
  st_face->vely = velx * geom->ny + vely * geom->my + velz * geom->py;
  st_face->velz = velx * geom->nz + vely * geom->mz + velz * geom->pz;
}




void face_set_scalar_states_and_fluxes(struct state *st_L, struct state *st_R, struct state_face *st_face, struct fluxes *flux)
{
#if defined(MAXSCALARS) || defined(TGCHEM)
  int i;
#endif


#ifdef MRT
  double rs = 0.0 ;
  double ls = 0.0 ;
  if(flux->mass > 0)
    ls = 1.0 ;
  else
    rs = 1.0;

  flux->nHI = (ls*st_L->nHI + rs*st_R->nHI) * flux->mass ;
  flux->nHII = (ls*st_L->nHII + rs*st_R->nHII) * flux->mass ;
  flux->ne = (ls*st_L->ne + rs*st_R->ne) * flux->mass ;
#ifdef MRT_INCLUDE_HE
  flux->nHeI = (ls*st_L->nHeI + rs*st_R->nHeI) * flux->mass ;
  flux->nHeII = (ls*st_L->nHeII + rs*st_R->nHeII) * flux->mass ;
  flux->nHeIII = (ls*st_L->nHeIII + rs*st_R->nHeIII) * flux->mass ;
#endif
#endif

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
  /* choose upwind entropy for the entropy flux, i.e. we suppress dissipation */
  if(flux->mass > 0)
    st_face->A = st_L->A;
  else
    st_face->A = st_R->A;
  
  flux->entropy = flux->mass * log(st_face->A);
#endif


#ifdef MAXSCALARS

  double normfac, normifac;

  if(flux->mass > 0)
    st_face->scalars = st_L->scalars;
  else
    st_face->scalars = st_R->scalars;

  /* Normalize species here (other than those handled by SGCHEM) */
  normfac = 0;

  for(i = 0; i < N_Scalar; i++)
    {
      flux->scalars[i] = st_face->scalars[i] * flux->mass;

      if(scalar_elements[i].type == SCALAR_TYPE_SPECIES)
        normfac += st_face->scalars[i];
    }

  if(normfac != 0)
    {
      normifac = 1.0 / normfac;

      for(i = 0; i < N_Scalar; i++)
        if(scalar_elements[i].type == SCALAR_TYPE_SPECIES || scalar_elements[i].type == SCALAR_TYPE_NORMALIZE)
          flux->scalars[i] *= normifac;
    }

#ifdef SGCHEM
  double normfac_H,  normifac_H;
  double normfac_He, normifac_He;
  double normfac_C,  normifac_C;
  double normfac_D,  normifac_D;
  double normfac_M,  normifac_M;
  double normifac_HD;
  double normfac_Celem, normfac_Oelem, normfac_Melem;
  long sgflux_count, new_sgflux_count;
#if CHEMISTRYNETWORK == 15
  double element_flux[3];
  double species_flux[NSPEC_CMA];
#endif

  normfac_H  = normfac_He  = normfac_C  = normfac_D  = normfac_M   = normfac_Celem = normfac_Oelem = normfac_Melem = 0;
  normifac_H = normifac_He = normifac_C = normifac_D = normifac_HD = normifac_M = 1.0;
  sgflux_count = new_sgflux_count = 0;

  for(i = 0; i < N_Scalar; i++)
    {
      /* Handle the species in the SGCHEM set; details depend on choice of chemistry network */
#if CHEMISTRYNETWORK == 15
      switch (scalar_elements[i].type) {
	/* Helium and M can be dealt with using standard CMA */
      case SCALAR_TYPE_HELIUM:
        normfac_He += st_face->scalars[i];
        break;
      case SCALAR_TYPE_METAL:
        normfac_M  += st_face->scalars[i];
        break;
      case SCALAR_TYPE_MCMA:
	/* Other species need modified CMA */
        species_flux[sgflux_count] = st_face->scalars[i];
        sgflux_count++;
        break;
      default:
        break;
      }
      /* For debugging: should only have 10 fluxes at this point for network 15 */
      //assert(sgflux_count == NSPEC_CMA);
#else /* CHEMISTRYNETWORK == 15 */
      if(scalar_elements[i].type == SCALAR_TYPE_H2)
        normfac_H += 2.0 * st_face->scalars[i];
#if CHEMISTRYNETWORK == 1
      if(scalar_elements[i].type == SCALAR_TYPE_HYDROGEN || scalar_elements[i].type == SCALAR_TYPE_HD)
        normfac_H += st_face->scalars[i];

      if(scalar_elements[i].type == SCALAR_TYPE_DEUTERIUM || scalar_elements[i].type == SCALAR_TYPE_HD)
        normfac_D += st_face->scalars[i];

      if(scalar_elements[i].type == SCALAR_TYPE_HELIUM)
        normfac_He += st_face->scalars[i];
#endif
#if CHEMISTRYNETWORK == 4
      if(scalar_elements[i].type == SCALAR_TYPE_HYDROGEN)
        normfac_H += st_face->scalars[i];
#endif
#if CHEMISTRYNETWORK == 5
      if(scalar_elements[i].type == SCALAR_TYPE_HYDROGEN)
        normfac_H += st_face->scalars[i];
      if(scalar_elements[i].type == SCALAR_TYPE_CARBON)
        normfac_C += st_face->scalars[i];
#endif
#endif /* CHEMISTRYNETWORK == 15 */

#ifdef SGCHEM_VARIABLE_Z
      if(scalar_elements[i].type == SCALAR_TYPE_CELEM)
        normfac_Celem += st_face->scalars[i];

      if(scalar_elements[i].type == SCALAR_TYPE_OELEM)
        normfac_Oelem += st_face->scalars[i];

      if(scalar_elements[i].type == SCALAR_TYPE_MELEM)
        normfac_Melem += st_face->scalars[i];
#endif
    }

  /* For debugging: should only have 10 fluxes at this point for network 15 */
  //assert(sgflux_count == NSPEC_CMA);
#if CHEMISTRYNETWORK == 15
#ifndef MCMA
#error "SGCHEM network 15 requires MCMA"
#endif
  /* Correct the fluxes for the MCMA species */
  element_flux[0] = 1.0;
#ifdef SGCHEM_VARIABLE_Z
  element_flux[1] = normfac_Celem;
  element_flux[2] = normfac_Oelem;
#else
  element_flux[1] = All.CarbAbund;
  element_flux[2] = All.OxyAbund;
#endif

  CMA_CORRECT(species_flux, element_flux, &sgflux_count);

  if (normfac != 0)
    normifac = 1.0 / normfac;
  if (normfac_He != 0)
    normifac_He = ABHE / normfac_He;
  if (normfac_M != 0)
#ifdef SGCHEM_VARIABLE_Z
    normifac_M = normfac_Melem / normfac_M;
#else
    normifac_M = All.MAbund / normfac_M;
#endif

  new_sgflux_count = 0;
  for(i = 0; i < N_Scalar; i++) {
    switch (scalar_elements[i].type) {
    case SCALAR_TYPE_HELIUM:
      flux->scalars[i] *= normifac_He;
      break;
    case SCALAR_TYPE_METAL:
      flux->scalars[i] *= normifac_M;
      break;
    case SCALAR_TYPE_MCMA:
      flux->scalars[i] = species_flux[new_sgflux_count] * flux->mass;
      new_sgflux_count++;
      break;
    default:
      break;
    }
  }
  /* For debugging: check that we're counting correctly */
  assert(sgflux_count == new_sgflux_count);
#else /* CHEMISTRYNETWORK == 15 */
  if(normfac_H != 0)
    normifac_H = 1.0 / normfac_H;
  if(normfac_C != 0)
#ifdef SGCHEM_VARIABLE_Z
    normifac_C = normfac_Celem / normfac_C;
#else
    normifac_C = All.CarbAbund / normfac_C;
#endif

#if CHEMISTRYNETWORK == 1
  if (normfac_He != 0)
    normifac_He = ABHE / normfac_He;
  if (normfac_D != 0)
    normifac_D = All.DeutAbund / normfac_D;
  /* Kluge for HD: use either scaling for H or for D, depending on
   * which is most different from zero. In the future, we will switch 
   * to MCMA, which avoids this problem
   */
  if (normifac_H == 1.0) {
    normifac_HD = normifac_D;
  }
  else if (normifac_H > 1.0) {
    if (normifac_D > 1.0)
      normifac_HD = (normifac_H > normifac_D) ? normifac_H : normifac_D;
    if (normifac_D < 1.0)
      normifac_HD = (normifac_H > 1.0/normifac_D) ? normifac_H : normifac_D;
  }
  else if (normifac_H < 1.0) {
    if (normifac_D > 1.0)
      normifac_HD = (1.0/normifac_H > normifac_D) ? normifac_H : normifac_D;
    if (normifac_D < 1.0)
      normifac_HD = (1.0/normifac_H > 1.0/normifac_D) ? normifac_H : normifac_D;
  }
#endif

  for(i = 0; i < N_Scalar; i++)
    {
      if(scalar_elements[i].type == SCALAR_TYPE_HYDROGEN || scalar_elements[i].type == SCALAR_TYPE_H2)
        flux->scalars[i] *= normifac_H;
#if CHEMISTRYNETWORK == 1
      if(scalar_elements[i].type == SCALAR_TYPE_DEUTERIUM)
	flux->scalars[i] *= normifac_D;
      if(scalar_elements[i].type == SCALAR_TYPE_HELIUM)
	flux->scalars[i] *= normifac_He;
      if(scalar_elements[i].type == SCALAR_TYPE_HD)
	flux->scalars[i] *= normifac_HD;
#endif
#if CHEMISTRYNETWORK == 5
      if(scalar_elements[i].type == SCALAR_TYPE_CARBON)
        flux->scalars[i] *= normifac_C;
#endif
    }
#endif /* CHEMISTRYNETWORK == 15 */
#endif /* SGCHEM */

#endif


#ifdef TRACER_FIELD
  if(flux->mass > 0)
    st_face->tracer = st_L->tracer;
  else
    st_face->tracer = st_R->tracer;

  flux->tracer = st_face->tracer * flux->mass;
#endif

#ifdef TGCHEM
  for(i = 0; i < TGCHEM_NUM_ABUNDANCES; i++)
    {
      if(flux->mass > 0)
        st_face->pcabund[i] = st_L->pcabund[i];
      else
        st_face->pcabund[i] = st_R->pcabund[i];

      flux->pcabund[i] = st_face->pcabund[i] * flux->mass;
    }
#endif


}

#if defined(RIEMANN_HLLC) || defined(RIEMANN_ROSUNOV) || defined(RIEMANN_HLLD) || defined(RIEMANN_HLL) || defined(SPECIAL_RELATIVITY) || defined(GENERAL_RELATIVITY)
void flux_convert_to_lab_frame(struct state *st_L, struct state *st_R, double *vel_face, struct fluxes *flux)
{
#if (defined(SPECIAL_RELATIVITY) || defined(GENERAL_RELATIVITY)) && !defined(VORONOI_STATIC_MESH) 
#error "Fixme, before using RELATIVITY with moving mesh."
#endif
    
  double momx = flux->momentum[0];
  double momy = flux->momentum[1];
  double momz = flux->momentum[2];

  flux->momentum[0] += vel_face[0] * flux->mass;
  flux->momentum[1] += vel_face[1] * flux->mass;
  flux->momentum[2] += vel_face[2] * flux->mass;

  flux->energy += momx * vel_face[0] + momy * vel_face[1] + momz * vel_face[2] + 0.5 * flux->mass * (vel_face[0] * vel_face[0] + vel_face[1] * vel_face[1] + vel_face[2] * vel_face[2]);

#ifdef MHD
  double Bx;
  Bx = 0.5 * (st_L->Bx + st_R->Bx);

  flux->B[0] -= vel_face[0] * Bx;
  flux->B[1] -= vel_face[1] * Bx;
  flux->B[2] -= vel_face[2] * Bx;
#endif
}
#endif


void face_turn_momentum_flux(struct fluxes *flux, struct geometry *geom)
{
  /* flux->momentum vector needs to be turned in case the HLLC or Rosunov Riemann solvers are used */

  double momx = flux->momentum[0];
  double momy = flux->momentum[1];
  double momz = flux->momentum[2];

  flux->momentum[0] = momx * geom->nx + momy * geom->mx + momz * geom->px;
  flux->momentum[1] = momx * geom->ny + momy * geom->my + momz * geom->py;
  flux->momentum[2] = momx * geom->nz + momy * geom->mz + momz * geom->pz;


#ifdef MHD
  double Bx = flux->B[0];
  double By = flux->B[1];
  double Bz = flux->B[2];

  flux->B[0] = Bx * geom->nx + By * geom->mx + Bz * geom->px;
  flux->B[1] = Bx * geom->ny + By * geom->my + Bz * geom->py;
  flux->B[2] = Bx * geom->nz + By * geom->mz + Bz * geom->pz;
#endif
}

void face_get_fluxes(struct state *st_L, struct state *st_R, struct state_face *st_face, struct fluxes *flux, struct geometry *geom, double *vel_face)
{
  double fac;

  /* calculate fluxes for ordinary Riemann solver */

  fac = (st_face->velx - vel_face[0]) * geom->nx + (st_face->vely - vel_face[1]) * geom->ny + (st_face->velz - vel_face[2]) * geom->nz;

  flux->mass = st_face->rho * fac;

  flux->momentum[0] = (st_face->rho * st_face->velx * fac + st_face->press * geom->nx);
  flux->momentum[1] = (st_face->rho * st_face->vely * fac + st_face->press * geom->ny);
  flux->momentum[2] = (st_face->rho * st_face->velz * fac + st_face->press * geom->nz);


#ifndef ISOTHERM_EQS
#ifdef TGCHEM
  flux->energy = (0.5 * st_face->rho * (st_face->velx * st_face->velx +
                                        st_face->vely * st_face->vely +
                                        st_face->velz * st_face->velz) +
                  st_face->press / (st_face->gamma - 1)) * fac + st_face->press * (st_face->velx * geom->nx + st_face->vely * geom->ny + st_face->velz * geom->nz);
#else
#ifdef VARIABLE_GAMMA
  flux->energy = (0.5 * st_face->rho * (st_face->velx * st_face->velx +
                                        st_face->vely * st_face->vely +
                                        st_face->velz * st_face->velz) +
                  st_face->press / (st_face->gammaE - 1.0)) * fac + st_face->press * (st_face->velx * geom->nx + st_face->vely * geom->ny + st_face->velz * geom->nz);
#else
  flux->energy = (0.5 * st_face->rho * (st_face->velx * st_face->velx +
                                        st_face->vely * st_face->vely +
                                        st_face->velz * st_face->velz) +
                  st_face->press / GAMMA_MINUS1) * fac + st_face->press * (st_face->velx * geom->nx + st_face->vely * geom->ny + st_face->velz * geom->nz);
#endif
#endif
#endif

#if defined (VS_TURB) || defined (AB_TURB)
  flux->energy += st_face->rho * fac * get_turb_pot(geom->cx, geom->cy, geom->cz);
#endif

#ifdef FLD
  if(flux->mass > 0)
    {
      st_face->n_gamma = st_L->n_gamma;
      st_face->R2 = st_L->R2;
    }
  else
    {
      st_face->n_gamma = st_R->n_gamma;
      st_face->R2 = st_R->R2;
    }

  flux->dFLD = (3. - st_face->R2) / 2. * fac * st_face->n_gamma;
#endif
}


void face_limit_fluxes(struct state *st_L, struct state *st_R, struct state *st_center_L, struct state *st_center_R, struct fluxes *flux, double dt, double *count, double *count_reduced)
{
  *count = *count + 1.0;

  /* choose upwind mass to determine a stability bound on the maximum allowed mass exchange,
     (we do this to prevent negative masses under all circumstances) */

  double upwind_mass, upwind_activearea, reduc_fac;
  integertime upwind_timebin, downstream_timebin;
#ifdef GFM_NORMALIZED_METAL_ADVECTION
  double *upwind_scalars;    /* cell centered upwind scalars */
#ifdef GFM_CHEMTAGS
  double *upwind_chemtagsfraction;
#endif
#endif

  if(flux->mass > 0)
    {
      upwind_mass = st_L->oldmass;
      upwind_activearea = st_L->activearea;
      upwind_timebin = st_L->timeBin;
      downstream_timebin = st_R->timeBin;
#ifdef GFM_NORMALIZED_METAL_ADVECTION
      upwind_scalars = st_center_L->scalars;
#ifdef GFM_CHEMTAGS
      upwind_chemtagsfraction = st_center_L->chemtagsfraction;
#endif
#endif
    }
  else
    {
      upwind_mass = st_R->oldmass;
      upwind_activearea = st_R->activearea;
      upwind_timebin = st_R->timeBin;
      downstream_timebin = st_L->timeBin;
#ifdef GFM_NORMALIZED_METAL_ADVECTION
      upwind_scalars = st_center_R->scalars;
#ifdef GFM_CHEMTAGS
      upwind_chemtagsfraction = st_center_R->chemtagsfraction;
#endif
#endif
    }

  if(upwind_timebin > downstream_timebin)
    dt *= pow(2, upwind_timebin - downstream_timebin);

  if(fabs(flux->mass * dt * upwind_activearea) > 0.9 * upwind_mass)
    {
      reduc_fac = 0.9 * upwind_mass / fabs(flux->mass * dt * upwind_activearea);

      *count_reduced = *count_reduced + 1.0;

      flux->mass *= reduc_fac;
      flux->energy *= reduc_fac;
      flux->momentum[0] *= reduc_fac;
      flux->momentum[1] *= reduc_fac;
      flux->momentum[2] *= reduc_fac;

      /* remark: do not reduce the magnetic field flux, as it is not coupled to the mass flux */
#ifdef USE_ENTROPY_FOR_COLD_FLOWS
      flux->entropy *= reduc_fac;
#endif
#ifdef MAXSCALARS
      for(int i = 0; i < N_Scalar; i++)
        flux->scalars[i] *= reduc_fac;
#endif
#ifdef TRACER_FIELD
      flux->tracer *= reduc_fac;
#endif

#ifdef TGCHEM
      for(int i = 0; i < TGCHEM_NUM_ABUNDANCES; i++)
        flux->pcabund[i] *= reduc_fac;
#endif
    }

  
  
#ifdef GFM_NORMALIZED_METAL_ADVECTION
  int flag = 0;

  for(int i = 0; i < N_Scalar; i++)
    {
      if(scalar_elements[i].type == SCALAR_TYPE_SPECIES)
        {
          if(flux->scalars[i] * flux->mass < 0) /* they have opposite sign */
            flag = 1;

          /* there is a danger that we drain all the metals from the cell */
          if(fabs(flux->scalars[i] * dt * upwind_activearea) > 0.9 * upwind_mass * fabs(upwind_scalars[i]))
            flag = 1;
        }
    }

  if(flag) /* in this case, revert to 1st order donor cell advection, and make sure that the abundance vector is normalized */
    {
      /* calculate normalization factor */
      double normfac = 0;

      for(int i = 0; i < N_Scalar; i++)
        {
          if(scalar_elements[i].type == SCALAR_TYPE_SPECIES)
            normfac += upwind_scalars[i];
        }

      if(normfac != 0)
        {
          double normifac = 1.0 / normfac;

          for(int i = 0; i < N_Scalar; i++)
            {
              if(scalar_elements[i].type == SCALAR_TYPE_SPECIES)
                flux->scalars[i] = normifac * upwind_scalars[i] * flux->mass;
            }
        }
    }

#ifdef GFM_CHEMTAGS
  static int first_entry_of_metals = -1;
  
  if(first_entry_of_metals < 0)
    for(int i = 0; i < N_Scalar; i++)
      if(scalar_elements[i].type == SCALAR_TYPE_SPECIES)
        {
          first_entry_of_metals = i;
          break;
        }

  double flux_metals = 0;
  for(int j = 2; j < GFM_N_CHEM_ELEMENTS; j++)
    flux_metals += flux->scalars[first_entry_of_metals + j];    /* total metal mass flux */

  double sum_frac = 
    upwind_chemtagsfraction[GFM_SNIA_CHEMTAG] + upwind_chemtagsfraction[GFM_SNII_CHEMTAG] + upwind_chemtagsfraction[GFM_AGB_CHEMTAG]; 

  if(sum_frac > 0)
    {
      flux->chemtags[GFM_SNIA_CHEMTAG] = flux_metals * upwind_chemtagsfraction[GFM_SNIA_CHEMTAG] / sum_frac;
      flux->chemtags[GFM_SNII_CHEMTAG] = flux_metals * upwind_chemtagsfraction[GFM_SNII_CHEMTAG] / sum_frac;
      flux->chemtags[GFM_AGB_CHEMTAG] =  flux_metals * upwind_chemtagsfraction[GFM_AGB_CHEMTAG]  / sum_frac;
    }
  else
    {
      flux->chemtags[GFM_SNIA_CHEMTAG] = 0;
      flux->chemtags[GFM_SNII_CHEMTAG] = 0;
      flux->chemtags[GFM_AGB_CHEMTAG] = 0;
    }

#if defined(GFM_SPLITFE)
  double flux_iron = flux->scalars[first_entry_of_metals + element_index_Iron];

  double sumironfrac = 0;

#ifdef GFM_SPLITFE_ADDINAGB
  sumironfrac += upwind_chemtagsfraction[GFM_FESNIA_CHEMTAG];
  sumironfrac += upwind_chemtagsfraction[GFM_FESNII_CHEMTAG];
#else
  sumironfrac = upwind_scalars[first_entry_of_metals + element_index_Iron]; /* this includes the AGB iron */
#endif

  if(sumironfrac > 0)
    {
      flux->chemtags[GFM_FESNIA_CHEMTAG] = flux_iron * upwind_chemtagsfraction[GFM_FESNIA_CHEMTAG] / sumironfrac;
      flux->chemtags[GFM_FESNII_CHEMTAG] = flux_iron * upwind_chemtagsfraction[GFM_FESNII_CHEMTAG] / sumironfrac;
    }
  else
    {
      flux->chemtags[GFM_FESNIA_CHEMTAG] = 0;
      flux->chemtags[GFM_FESNII_CHEMTAG] = 0;
    }
#endif

#ifdef GFM_RPROCESS
  if(sum_frac > 0)
    flux->chemtags[GFM_NSNS_CHEMTAG] = flux_metals * upwind_chemtagsfraction[GFM_NSNS_CHEMTAG] / sum_frac;
  else
    flux->chemtags[GFM_NSNS_CHEMTAG] = 0;
#endif



#endif  // end of GFM_CHEMTAGS
#endif  // end of GFM_NORMALIZED_METAL_ADVECTION
}


void face_clear_fluxes(struct fluxes *flux)
{
  flux->mass = 0;
  flux->momentum[0] = 0;
  flux->momentum[1] = 0;
  flux->momentum[2] = 0;
  flux->energy = 0;
#ifdef MHD
  flux->B[0] = 0;
  flux->B[1] = 0;
  flux->B[2] = 0;
#endif
#ifdef MHD_CT
  flux->A[0] = 0;
  flux->A[1] = 0;
  flux->A[2] = 0;
#endif
#ifdef MRT_LSF_GRADIENTS
  flux->nHI = 0 ;
  flux->nHII = 0 ;
  flux->nHeI = 0 ;
  flux->ne = 0 ;
  flux->nHeII = 0 ;
  flux->nHeIII = 0 ;
#endif
}

void face_add_fluxes_advection(struct state_face *st_face, struct fluxes *flux, struct geometry *geom, double *vel_face)
{
  double fac = -vel_face[0] * geom->nx - vel_face[1] * geom->ny - vel_face[2] * geom->nz;

  flux->mass += st_face->rho * fac;

  flux->momentum[0] += st_face->rho * st_face->velx * fac;
  flux->momentum[1] += st_face->rho * st_face->vely * fac;
  flux->momentum[2] += st_face->rho * st_face->velz * fac;

#ifdef TGCHEM
  flux->energy += 0.5 * st_face->rho * fac * (st_face->velx * st_face->velx + st_face->vely * st_face->vely + st_face->velz * st_face->velz) + st_face->press / (st_face->gamma - 1) * fac;
#else
  flux->energy += 0.5 * st_face->rho * fac * (st_face->velx * st_face->velx + st_face->vely * st_face->vely + st_face->velz * st_face->velz) + st_face->press / GAMMA_MINUS1 * fac;
#endif
}

int flux_list_data_compare(const void *a, const void *b)
{
  if(((struct flux_list_data *) a)->task < (((struct flux_list_data *) b)->task))
    return -1;

  if(((struct flux_list_data *) a)->task > (((struct flux_list_data *) b)->task))
    return +1;

  return 0;
}



void apply_flux_list(void)
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


#ifdef SPECIAL_BOUNDARY
      if((P[p].ID == -2) || (P[p].ID <= -3))
        continue;
#endif

#if defined(BOUNDARY_INFLOWOUTFLOW_MINID) && defined(BOUNDARY_INFLOWOUTFLOW_MAXID)
      if(P[p].ID >= BOUNDARY_INFLOWOUTFLOW_MINID && P[p].ID < BOUNDARY_INFLOWOUTFLOW_MAXID)
        continue;
#endif

#if defined(BOUNDARY_REFL_SOLIDSIDE_MINID) && defined(BOUNDARY_REFL_SOLIDSIDE_MAXID)
      if(P[p].ID >= BOUNDARY_REFL_SOLIDSIDE_MINID && P[p].ID < BOUNDARY_REFL_SOLIDSIDE_MAXID)
        continue;
#endif


#ifdef COFFEE_PROBLEM
      if(P[p].ID >= 20000000)
        continue;
#endif

      P[p].Mass += FluxListGet[i].dM;

#if defined(TRACER_MC)
      if(FluxListGet[i].dM < 0)
        {
          if(FluxListGet[i].pother_index >= All.MaxPart)
            terminate("FluxListGet[i].pother_index >= All.MaxPart");

          consider_moving_tracers(p, FluxListGet[i].pother_task, FluxListGet[i].pother_index, FluxListGet[i].pother_ID, -FluxListGet[i].dM / CellTracerMasses[p]);
          CellTracerMasses[p] += FluxListGet[i].dM;
        }
#endif

      SphP[p].Momentum[0] += FluxListGet[i].dP[0];
      SphP[p].Momentum[1] += FluxListGet[i].dP[1];
      SphP[p].Momentum[2] += FluxListGet[i].dP[2];

#ifdef MRT

      SphP[p].nHI += FluxListGet[i].dnHI ; 
      SphP[p].nHII += FluxListGet[i].dnHII ; 
      SphP[p].ne += FluxListGet[i].dne ;
#ifdef MRT_INCLUDE_HE
      SphP[p].nHeI += FluxListGet[i].dnHeI ; 
      SphP[p].nHeII += FluxListGet[i].dnHeII ; 
      SphP[p].nHeIII += FluxListGet[i].dnHeIII ;
#endif
#endif

#ifdef MHD
      SphP[p].BConserved[0] += FluxListGet[i].dB[0];
      SphP[p].BConserved[1] += FluxListGet[i].dB[1];
      SphP[p].BConserved[2] += FluxListGet[i].dB[2];
#ifdef MHD_CT
      SphP[p].AConserved[0] += FluxListGet[i].dA[0];
      SphP[p].AConserved[1] += FluxListGet[i].dA[1];
      SphP[p].AConserved[2] += FluxListGet[i].dA[2];
#endif
#ifdef MHD_THERMAL_ENERGY_SWITCH
      SphP[p].Etherm += FluxListGet[i].dEtherm;
#endif
#endif

#ifdef MAXSCALARS
      for(k = 0; k < N_Scalar; k++)
        *(MyFloat *) (((char *) (&SphP[p])) + scalar_elements[k].offset_mass) += FluxListGet[i].dConservedScalars[k];
#endif

#ifdef GFM_CHEMTAGS
      for(int m = 0; m < GFM_N_CHEM_TAGS; m++)
        SphP[p].MassMetalsChemTags[m] += FluxListGet[i].dMassMetalsChemTags[m];
#endif

#ifndef ISOTHERM_EQS
      SphP[p].Energy += FluxListGet[i].dEnergy;
#endif

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
      SphP[p].Entropy += FluxListGet[i].dEntropy;
#endif

#ifdef TRACER_FIELD
      SphP[p].ConservedTracer += FluxListGet[i].dConservedTracer;
#endif

#ifdef OUTPUT_CELL_SPIN
      SphP[p].Spin[0] += FluxListGet[i].dSpin[0];
      SphP[p].Spin[1] += FluxListGet[i].dSpin[1];
      SphP[p].Spin[2] += FluxListGet[i].dSpin[2];

      SphP[p].CenterOffsetMass[0] += FluxListGet[i].dCenterOffsetMass[0];
      SphP[p].CenterOffsetMass[1] += FluxListGet[i].dCenterOffsetMass[1];
      SphP[p].CenterOffsetMass[2] += FluxListGet[i].dCenterOffsetMass[2];
#endif

#ifdef TGCHEM
      for(k = 0; k < TGCHEM_NUM_ABUNDANCES; k++)
        AbundFluxes[TGCHEM_NUM_ABUNDANCES * p + k] += FluxListGet[i].Abund[k];
#endif

#ifdef COSMIC_RAYS
      SphP[p].CR_Energy += FluxListGet[i].dCR_Energy;
#endif
#ifdef FLD
      SphP[p].n_gamma += FluxListGet[i].dFLD;
#endif
    }


  myfree(FluxListGet);
}

void fvs_initialize_statistics(struct fvs_stat *stat)
{
  stat->count_disable_extrapolation = 0;
#ifdef COSMIC_RAYS
  stat->count_CR_limiter = 0;
#endif
}

void fvs_evaluate_statistics(struct fvs_stat *stat)
{
#ifdef VERBOSE
  int count_disable_extrapolation = 0;
  MPI_Reduce(&stat->count_disable_extrapolation, &count_disable_extrapolation, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  mpi_printf("FLUX: Disabled extrapolation for %d interfaces.\n", count_disable_extrapolation);
#endif

#ifdef COSMIC_RAYS
  int count_CR_limiter = 0;
  MPI_Reduce(&stat->count_CR_limiter, &count_CR_limiter, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  mpi_printf("COSMIC_RAYS: limited CR pressure extrapolation for %d interfaces.\n", count_CR_limiter);
#endif
}

#ifdef AXISYMMETRY
void apply_axisymmetric_source_terms(tessellation * T)
{
  int idx, i;
  double dt, atime = 1;
  struct state state, delta_time, delta_space;

  if(All.ComovingIntegrationOn)
    terminate("axisymmetric option not implemented for cosmological integrations");

  memset(&delta_space, 0, sizeof(delta_space)); /* make sure that the spatial part is zero */
  memset(&delta_time, 0, sizeof(delta_time));

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(T->DP[0].task != ThisTask)
        terminate("DP[0].task != ThisTask");

      face_get_state(T, -100, i, &state);

      state.velx = state.velGas[0];
      state.vely = state.velGas[1];
      state.velz = state.velGas[2];

      dt = (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;

      state.dt_half = 0.5 * dt;

      face_do_time_extrapolation(&delta_time, &state, atime);

      face_add_extrapolations(&state, &delta_time, &delta_space);

      /* we have now extrapolated the state of the cell half a timestep to the future. 
         Let's now add the axisymmetric source terms with this */

      /* radial momentum */
      SphP[i].Momentum[0] += dt * SphP[i].Volume / P[i].Pos[0] * (state.rho * state.velz * state.velz + state.press);

      /* azimuthal momentum */
      SphP[i].Momentum[2] += dt * SphP[i].Volume / P[i].Pos[0] * (-state.rho * state.velx * state.velz);
    }
}
#endif


#ifdef ONEDIMS_SPHERICAL
void apply_spherical_source_terms()
{
  int idx, i;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      double Pressure = SphP[i].Pressure;
      double dt_Extrapolation = All.Time - SphP[i].TimeLastPrimUpdate;
      struct grad_data *grad = &SphP[i].Grad;
#ifndef VARIABLE_GAMMA
      Pressure +=
        -dt_Extrapolation * (GAMMA * Pressure * (grad->dvel[0][0] + grad->dvel[1][1] + grad->dvel[2][2]) + P[i].Vel[0] * grad->dpress[0] + P[i].Vel[1] * grad->dpress[1] +
                             P[i].Vel[2] * grad->dpress[2]);
#else
      Pressure +=
        -dt_Extrapolation * (SphP[i].GammaE * Pressure * (grad->dvel[0][0] + grad->dvel[1][1] + grad->dvel[2][2]) + P[i].Vel[0] * grad->dpress[0] + P[i].Vel[1] * grad->dpress[1] +
                             P[i].Vel[2] * grad->dpress[2]);
#endif

      double dt = 0.5 * (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;
      SphP[i].Momentum[0] += dt * Pressure * (Mesh.VF[i + 1].area - Mesh.VF[i].area);
    }
}
#endif

#ifdef OUTPUT_CELL_SPIN
void add_spin_source_term_from_grid_movement()
{
  int i, k;
  for(i = 0; i < NumGas; i++)
    if(P[i].Type == 0 && P[i].ID != 0 && P[i].Mass > 0.)
      {
        double dt = (P[i].TimeBin ? (((integertime) 1) << P[i].TimeBin) : 0) * All.Timebase_interval;

        {
          double xtmp, ytmp, ztmp;
          SphP[i].CenterOffset[0] = NEAREST_X(-SphP[i].Center[0] + SphP[i].CenterOld[0] + SphP[i].Momentum[0] * dt / P[i].Mass + SphP[i].CenterOffsetMass[0] / P[i].Mass);
          SphP[i].CenterOffset[1] = NEAREST_Y(-SphP[i].Center[1] + SphP[i].CenterOld[1] + SphP[i].Momentum[1] * dt / P[i].Mass + SphP[i].CenterOffsetMass[1] / P[i].Mass);
          SphP[i].CenterOffset[2] = NEAREST_Z(-SphP[i].Center[2] + SphP[i].CenterOld[2] + SphP[i].Momentum[2] * dt / P[i].Mass + SphP[i].CenterOffsetMass[2] / P[i].Mass);
        }

        for(k = 0; k < 3; k++)
          SphP[i].CenterOffsetMass[k] = 0;

        double dx, dy, dz;
        {
          double xtmp, ytmp, ztmp;
          dx = NEAREST_X(SphP[i].Center[0] + SphP[i].CenterOffset[0] - SphP[i].CenterOld[0]);
          dy = NEAREST_Y(SphP[i].Center[1] + SphP[i].CenterOffset[1] - SphP[i].CenterOld[1]);
          dz = NEAREST_Z(SphP[i].Center[2] + SphP[i].CenterOffset[2] - SphP[i].CenterOld[2]);
        }

        SphP[i].Spin[0] -= dy * SphP[i].Momentum[2] - dz * SphP[i].Momentum[1];
        SphP[i].Spin[1] -= dz * SphP[i].Momentum[0] - dx * SphP[i].Momentum[2];
        SphP[i].Spin[2] -= dx * SphP[i].Momentum[1] - dy * SphP[i].Momentum[0];

        SphP[i].CenterOld[0] = SphP[i].Center[0] + SphP[i].CenterOffset[0];
        SphP[i].CenterOld[1] = SphP[i].Center[1] + SphP[i].CenterOffset[1];
        SphP[i].CenterOld[2] = SphP[i].Center[2] + SphP[i].CenterOffset[2];
      }
}
#endif
#endif
