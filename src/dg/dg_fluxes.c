/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dg/dg_fluxes.c
 * \date        10/2014
 * \author		Kevin Schaal
 * \brief		Calculation of the fluxes accross the interfaces
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

//3d done

#include "../allvars.h"
#include "../proto.h"

#ifdef DG

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dg_core_inline.h"


static struct flux_list_data
{
  int task, index;
  double dM[NOF_OUTER_QUAD_POINTS];     //flux for every quad point
  double dP[NOF_OUTER_QUAD_POINTS][3];
  double dEnergy[NOF_OUTER_QUAD_POINTS];
  enum E_interface interface;
}
 *FluxList;

static int Nflux, MaxNflux;

struct primexch *PrimExch;
struct grad_data *GradExch;

/*!< state on a face determined by riemann solver */
struct state_face state_face;

/*!< flux through a face */
struct fluxes fluxes;

struct geometry geom;

#ifdef DISCONTINUITY_DETECTION
static const double dd_epsilon = 0;
#endif

/*!
 * Calculate the change of the weights due to fluxes accross the interfaces
 * a2 = a2 - c * dt * R_outer(a1)
 */
void subtract_R_outer(tessellation * T, double c, CBV(a1), CBV(a2))
{

  TIMER_START(CPU_DG_OUTER);

  mpi_printf("DG: Calculating outer integral\n");

  int i, j, quadp;
  double count = 0, tot_count;
  double face_dt, hubble_a, atime;

#ifdef TWODIMS
  const double norm_fac = 0.5;
#else
  const double norm_fac = 0.25;
#endif

#ifdef DISCONTINUITY_DETECTION
  double density_temp[2];

  for(i = 0; i < NumGas; i++)
    {
      SphP[i].Inflow_boundaries = 0;
    }
#endif


  MaxNflux = T->Indi.AllocFacNflux;
  Nflux = 0;
  FluxList = mymalloc_movable(&FluxList, "FluxList", MaxNflux * sizeof(struct flux_list_data));

  face *VF = T->VF;
  point *DP = T->DP;

  int face_dir;
  int interface_for_state;
  int opposing_interface_for_state;

  for(i = 0; i < T->Nvf; i++)   //loop over interfaces
    {

      struct state state_L;
      struct state state_R;


      face_dt = 0;              /* the default is that this face is not active */

      /* calculate normal vectors */
      if(face_get_normals(T, i, &geom))
        continue;

      if(fabs(geom.nx) > 0.5)
        {
          face_dir = 1;

          if(geom.cx - DP[VF[i].p1].x > 0)
            {
              interface_for_state = RIGHT;
              opposing_interface_for_state = LEFT;
            }
          else
            {
              interface_for_state = LEFT;
              opposing_interface_for_state = RIGHT;
            }
        }
      else if(fabs(geom.ny) > 0.5)
        {
          face_dir = 2;

          if(geom.cy - DP[VF[i].p1].y > 0)
            {
              interface_for_state = FRONT;
              opposing_interface_for_state = BACK;
            }
          else
            {
              interface_for_state = BACK;
              opposing_interface_for_state = FRONT;
            }
        }
      else
        {
          face_dir = 3;

          if(geom.cz - DP[VF[i].p1].z > 0)
            {
              interface_for_state = TOP;
              opposing_interface_for_state = BOTTOM;
            }
          else
            {
              interface_for_state = BOTTOM;
              opposing_interface_for_state = TOP;
            }


        }

      //projection matrices
      double (*P_X_L)[NOF_BASE_FUNCTIONS];
      double (*P_X_R)[NOF_BASE_FUNCTIONS];

      if(DP[VF[i].p1].level == DP[VF[i].p2].level)
        {
          P_X_L = 0;
          P_X_R = 0;
        }
      else if(DP[VF[i].p1].level > DP[VF[i].p2].level)
        {
          assign_projection_matrix(VF[i].p1, VF[i].p2, &P_X_R);
          P_X_L = 0;
        }
      else                      //level p2 > level p1
        {
          assign_projection_matrix(VF[i].p2, VF[i].p1, &P_X_L);
          P_X_R = 0;
        }

      for(quadp = 0; quadp < Nof_outer_quad_points; quadp++)    //loop over quadrature points
        {

          /* get the values of the states at the quadrature point */
          if(dg_face_get_state(T, a1, interface_for_state, quadp, VF[i].p1, i, &state_L, P_X_L))
            continue;

          if(dg_face_get_state(T, a1, opposing_interface_for_state, quadp, VF[i].p2, i, &state_R, P_X_R))
            continue;

          /* only treat faces where one of the two sides is active */
          if(!TimeBinSynchronized[state_L.timeBin] && !TimeBinSynchronized[state_R.timeBin])
            continue;

          /* clarify whether the face should be done by this task (it may be present also on another task) */
          if(face_check_responsibility_of_this_task(T, VF[i].p1, VF[i].p2, &state_L, &state_R))
            continue;

          /* calculate timestep of the face */
          face_dt = face_timestep(&state_L, &state_R, &hubble_a, &atime);       //minimum timestep of the two states

          if(!(face_dt > 0))
            continue;

          /* now estimate the velocity of the midpoint of the face based on the velocities of the generators of the mesh. */
          double vel_face[3];

          if(All.ComovingIntegrationOn)
            for(j = 0; j < 3; j++)
              {
                state_L.velVertex[j] /= atime;  /* convert vertex motion to peculiar velocity */
                state_R.velVertex[j] /= atime;
              }

          vel_face[0] = 0;
          vel_face[1] = 0;
          vel_face[2] = 0;

          state_convert_to_local_frame(&state_L, vel_face, hubble_a, atime);
          state_convert_to_local_frame(&state_R, vel_face, hubble_a, atime);

          /* check for positivity */
          if(state_L.press < 0 || state_R.press < 0 || state_L.rho < 0 || state_R.rho < 0)
            {
              printf("i=%d press_L=%.16e press_R=%.16e rho_L=%.16e rho_R=%.16e\n", i, state_L.press, state_R.press, state_L.rho, state_R.rho);
              printf("area=%g lx=%g ly=%g   rx=%g ry=%g\n", VF[i].area, state_L.dx, state_L.dy, state_R.dx, state_R.dy);
              terminate("found negative values");
            }

          /* mirror velocity in case of reflecting boundaries */
          face_boundary_check(&T->DP[VF[i].p1], &state_L.velx, &state_L.vely, &state_L.velz);
          face_boundary_check(&T->DP[VF[i].p2], &state_R.velx, &state_R.vely, &state_R.velz);

          /* turn the velocities to get velx perpendicular and vely and velz in the plane of the face */
          face_turn_velocities(&state_L, &geom);
          face_turn_velocities(&state_R, &geom);

          /* call Riemann solver */

          double press;
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
          double vel_face_turned[3];
          // Set vel face to 0, since we're in AMR (see finite_volume_solver.c)
          vel_face_turned[0] = 0.0;
          vel_face_turned[1] = 0.0;
          vel_face_turned[2] = 0.0;
          press = godunov_flux_3d_hlld(&state_L, &state_R, vel_face_turned, &state_face, &fluxes);
#else
#ifdef RIEMANN_HLL
          press = godunov_flux_3d_hll(&state_L, &state_R, &state_face, &fluxes);
#else
          press = godunov_flux_3d(&state_L, &state_R, &state_face);     /* exact ideal gas solver */
#endif
#endif
#endif
#endif
#endif

          count++;

          //printf("Left state: rho=%e, p=%e, vx=%e, vy=%e, vz=%e\n", state_L.rho, state_L.press, state_L.velx, state_L.vely, state_L.velz);
          //printf("Right state: rho=%e, p=%e, vx=%e, vy=%e, vz=%e\n", state_R.rho, state_R.press, state_R.velx, state_R.vely, state_R.velz);

#ifdef RIEMANN_HLLC
          if(press < Epsilon_p)
            {
              //use HLL flux
              press = godunov_flux_3d_hll(&state_L, &state_R, &state_face, &fluxes);
            }
#endif

          if(press < 0)
            {
              printf("ID_L: %d, ID_R: %d, pressure: %e\n", VF[i].p1, VF[i].p2, press);
              printf("Interface: (%g, %g, %g), n=(%g, %g, %g)\n", VF[i].cx, VF[i].cy, VF[i].cz, geom.nx, geom.ny, geom.nz);
              printf("Left state: rho=%e, p=%e, vx=%e, vy=%e, vz=%e\n", state_L.rho, state_L.press, state_L.velx, state_L.vely, state_L.velz);
              printf("Right state: rho=%e, p=%e, vx=%e, vy=%e, vz=%e\n", state_R.rho, state_R.press, state_R.velx, state_R.vely, state_R.velz);
              terminate("press < 0");
            }

          /* turn the velocity field back */
          face_turnback_velocities(&state_face, &geom);

#if defined(RIEMANN_HLLC) || defined(RIEMANN_ROSUNOV) || defined(RIEMANN_HLLD) || defined(RIEMANN_HLL)  /* for non-exact Riemann solver, fluxes are already computed in the local frame, so convert to lab frame and turn momentum fluxes to the lab orientation  */
          face_turn_momentum_flux(&fluxes, &geom);
#else

          /* calculate fluxes for exact Riemann problem */
          /* compute net flux with dot-product of outward normal and area of face */
          /* multiplication with area and time-step comes later */

          face_get_fluxes(&state_L, &state_R, &state_face, &fluxes, &geom, vel_face);
#endif

          /* put in cosmological factors */
          if(All.ComovingIntegrationOn)
            {
              fluxes.momentum[0] *= atime;
              fluxes.momentum[1] *= atime;
              fluxes.momentum[2] *= atime;
              fluxes.energy *= atime * atime;

            }


          if(!gsl_finite(fluxes.energy))
            {
              printf("i=%d eFlux-Bummer: %g %g %g\n", i, fluxes.energy, state_face.press, state_face.rho);
              printf("rho_L=%g velx_L=%g vely_L=%g velz_L=%g press_L=%g\n", state_L.rho, state_L.velx, state_L.vely, state_L.velz, state_L.press);
              printf("rho_R=%g velx_R=%g vely_R=%g velz_R=%g press_R=%g\n", state_R.rho, state_R.velx, state_R.vely, state_R.velz, state_R.press);

              terminate("infinity encountered");
            }


          /* now apply the flux to update the conserved states of the cells */


          if(face_dt > 0)       /* selects active faces */
            {
              int k, p, q, l;
              double dir;
              enum E_interface interface;

              for(k = 0; k < 2; k++)
                {
                  if(k == 0)
                    {
                      q = VF[i].p1;
                      p = DP[q].index;
                      assert(p == q);

                      dir = 1;
                    }
                  else
                    {
                      q = VF[i].p2;
                      p = DP[q].index;
                      dir = -1;
                    }

                  //determine position of the edge

                  if(face_dir == 1)
                    {
                      if(geom.cx - DP[q].x > 0)
                        {
                          interface = RIGHT;
                        }
                      else
                        {
                          interface = LEFT;
                        }
                    }
                  else if(face_dir == 2)
                    {
                      if(geom.cy - DP[q].y > 0)
                        {
                          interface = FRONT;
                        }
                      else
                        {
                          interface = BACK;
                        }
                    }
                  else
                    {
                      assert(face_dir == 3);

                      if(geom.cz - DP[q].z > 0)
                        {
                          interface = TOP;
                        }
                      else
                        {
                          interface = BOTTOM;
                        }
                    }

                  if(DP[q].task == ThisTask)
                    {
                      if(DP[q].index >= NumGas) /* this is a local ghost point */
                        {
                          if(DP[VF[i].p1].ID == DP[VF[i].p2].ID)        /* this may happen for reflective points */
                            continue;
                          p -= NumGas;
                        }

                      /* note: this will be executed if P[p] is a local point, independent of active or not */

#ifdef DG_DEBUG
                      double dt = (P[p].TimeBinHydro ? (((integertime) 1) << P[p].TimeBinHydro) : 0) * All.Timebase_interval;
                      assert(dt == face_dt);
#endif


#ifdef DISCONTINUITY_DETECTION
                      if(quadp == 0)
                        {
                          density_temp[k] = a2[p][0][0];
                        }
#endif

                      for(l = 0; l < Nof_base_functions; l++)
                        {
                          //a2 = a2 - c * dt * R_outer(a1)
                          a2[p][l][0] =
                            a2[p][l][0] - c * face_dt * norm_fac * VF[i].area / SphP[p].Volume * dir * fluxes.mass * GET_Outer_base_values(interface, quadp,
                                                                                                                                           l) * GET_Outer_quad_points_weights(interface, quadp);
                          a2[p][l][1] =
                            a2[p][l][1] - c * face_dt * norm_fac * VF[i].area / SphP[p].Volume * dir * fluxes.momentum[0] * GET_Outer_base_values(interface, quadp,
                                                                                                                                                  l) * GET_Outer_quad_points_weights(interface, quadp);
                          a2[p][l][2] =
                            a2[p][l][2] - c * face_dt * norm_fac * VF[i].area / SphP[p].Volume * dir * fluxes.momentum[1] * GET_Outer_base_values(interface, quadp,
                                                                                                                                                  l) * GET_Outer_quad_points_weights(interface, quadp);
                          a2[p][l][3] =
                            a2[p][l][3] - c * face_dt * norm_fac * VF[i].area / SphP[p].Volume * dir * fluxes.momentum[2] * GET_Outer_base_values(interface, quadp,
                                                                                                                                                  l) * GET_Outer_quad_points_weights(interface, quadp);
                          a2[p][l][4] =
                            a2[p][l][4] - c * face_dt * norm_fac * VF[i].area / SphP[p].Volume * dir * fluxes.energy * GET_Outer_base_values(interface, quadp,
                                                                                                                                             l) * GET_Outer_quad_points_weights(interface, quadp);
                        }

#ifdef DISCONTINUITY_DETECTION

                      if(quadp == Nof_outer_quad_points - 1 && a2[p][0][0] > density_temp[k] + dd_epsilon)      // inflow boundary
                        {
                          SphP[p].Inflow_boundaries += (1 << interface);
                        }
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
                          FluxList = myrealloc_movable(FluxList, MaxNflux * sizeof(struct flux_list_data));

                          if(Nflux >= MaxNflux)
                            terminate("Nflux >= MaxNflux");
                        }

                      FluxList[Nflux].task = DP[q].task;
                      FluxList[Nflux].index = DP[q].originalindex;

                      FluxList[Nflux].dM[quadp] = c * face_dt * norm_fac * VF[i].area * dir * fluxes.mass;
                      FluxList[Nflux].dP[quadp][0] = c * face_dt * norm_fac * VF[i].area * dir * fluxes.momentum[0];
                      FluxList[Nflux].dP[quadp][1] = c * face_dt * norm_fac * VF[i].area * dir * fluxes.momentum[1];
                      FluxList[Nflux].dP[quadp][2] = c * face_dt * norm_fac * VF[i].area * dir * fluxes.momentum[2];
                      FluxList[Nflux].dEnergy[quadp] = c * face_dt * norm_fac * VF[i].area * dir * fluxes.energy;


                      FluxList[Nflux].interface = interface;

                      if(quadp == Nof_outer_quad_points - 1)
                        {
                          Nflux++;
                        }

                    }
                }               //end k=0/k=1 (left/right state)
            }
        }                       //end loop over quadrature points



    }                           //end loop over all faces



  TIMER_START(CPU_DG_FLIMBAL);

  MPI_Barrier(MPI_COMM_WORLD);

  TIMER_STOP(CPU_DG_FLIMBAL);

  /* now exchange the flux-list and apply it when needed */
  dg_apply_flux_list(a2);

  myfree(FluxList);

  MPI_Reduce(&count, &tot_count, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(ThisTask == 0)
    {
      printf("FLUX: exchanged fluxes over %g faces\n", tot_count);
      All.TotCountFluxes += tot_count;
    }

  TIMER_STOP(CPU_DG_OUTER);
}


int dg_face_get_state(tessellation * T, CBV(a), int interface, int q, int p, int i, struct state *st, double (*P_X)[NOF_BASE_FUNCTIONS])
{
  int particle;

  point *DP = T->DP;
  face *VF = T->VF;

  double w[5];                  //state at quadrature point
  double pv[5];                 //primitive variables at quadrature point

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
#if !defined(REFLECTIVE_X)
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
      calc_state_at_outer_quad_point(a[particle], interface, q, w, P_X);
      w_to_primvars(w, pv);

      if(pv[PRES] < 0 || pv[RHO] < 0)
        {
          /*
             double pos_abs[3];
             double cell_center[3];

             cell_center[0]=SphP[particle].Center[0];
             cell_center[1]=SphP[particle].Center[1];
             cell_center[2]=SphP[particle].Center[2];

             abs_pos_outer_quad_points(interface, q, amr_length[Mesh.DP[particle].level], cell_center, pos_abs);

             double w2[5];
             double pv2[5];
             calc_state_at_pos(a[particle], GET_Outer_quad_points_x(interface, q), GET_Outer_quad_points_y(interface, q), GET_Outer_quad_points_z(interface, q), w2);
             w_to_primvars(w2,pv2);

             printf("SECOND METHOD: pressure=%e, density=%e\n", pv2[PRES], pv2[RHO]);
             printf("cell center: (%f, %f, %f)\n", cell_center[0], cell_center[1], cell_center[2]);

           */

          printf("particle index: %d, level=%d, interface: %d, quadpoint: %d, p1:(%f, %f, %f), p2:(%f, %f, %f), pressure=%e, density=%e\n",
                 particle, Mesh.DP[particle].level, interface, q, Mesh.DP[VF[i].p1].x, Mesh.DP[VF[i].p1].y, Mesh.DP[VF[i].p1].z, Mesh.DP[VF[i].p2].x, Mesh.DP[VF[i].p2].y, Mesh.DP[VF[i].p2].z,
                 pv[PRES], pv[RHO]);
          terminate("Bad values at quadrature point!\n");
        }


      st->velGas[0] = pv[VX];
      st->velGas[1] = pv[VY];
      st->velGas[2] = pv[VZ];

      st->velVertex[0] = SphP[particle].VelVertex[0];
      st->velVertex[1] = SphP[particle].VelVertex[1];
      st->velVertex[2] = SphP[particle].VelVertex[2];

      st->rho = pv[RHO];
      st->press = pv[PRES];

      st->grad = &SphP[particle].Grad;

      st->timeBin = P[particle].TimeBinHydro;

      st->volume = SphP[particle].Volume;       //todocheckwhyneeded, also other quantities

      st->oldmass = SphP[particle].OldMass;
      st->surfacearea = SphP[particle].SurfaceArea;
      st->activearea = SphP[particle].ActiveArea;
      st->csnd = get_sound_speed(particle);
      st->ID = P[particle].ID;
    }
  else
    {
      calc_state_at_outer_quad_point(Aexch[particle], interface, q, w, P_X);
      w_to_primvars(w, pv);

      if(pv[PRES] < 0 || pv[RHO] < 0)
        {
          printf("particle index: %d, interface: %d, quadpoint: %d, p1:(%f, %f, %f), p2:(%f, %f, %f), pressure=%f, density=%f\n",
                 particle, interface, q, Mesh.DP[VF[i].p1].x, Mesh.DP[VF[i].p1].y, Mesh.DP[VF[i].p1].z, Mesh.DP[VF[i].p2].x, Mesh.DP[VF[i].p2].y, Mesh.DP[VF[i].p2].z, pv[PRES], pv[RHO]);
          terminate("Bad values at quadrature point!\n");
        }

      st->velGas[0] = pv[VX];
      st->velGas[1] = pv[VY];
      st->velGas[2] = pv[VZ];

      st->velVertex[0] = PrimExch[particle].VelVertex[0];
      st->velVertex[1] = PrimExch[particle].VelVertex[1];
      st->velVertex[2] = PrimExch[particle].VelVertex[2];

      st->rho = pv[RHO];
      st->press = pv[PRES];

      st->grad = &GradExch[particle];

      st->timeBin = PrimExch[particle].TimeBinHydro;

      st->volume = PrimExch[particle].Volume;

      st->oldmass = PrimExch[particle].OldMass;
      st->surfacearea = PrimExch[particle].SurfaceArea;
      st->activearea = PrimExch[particle].ActiveArea;
      st->csnd = PrimExch[particle].Csnd;
      st->ID = DP[p].ID;
    }

  /* check for reflecting or outflowing boundaries */
  face_boundary_check_vertex(T, p, &st->velVertex[0], &st->velVertex[1], &st->velVertex[2]);

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
  short int timeBin;
  double face_dt;

  /* take the minimum of the two */
  timeBin = state_L->timeBin;
  if(timeBin > state_R->timeBin)
    timeBin = state_R->timeBin;

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
      st->vely -= atime * hubble_a * st->dy;    //todowhyneeded?
      st->velz -= atime * hubble_a * st->dz;
    }
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
    }
  else
    {
      st_face->rho = st_R->rho;
      st_face->velx = st_R->velx;
      st_face->vely = st_R->vely;
      st_face->velz = st_R->velz;
      st_face->press = st_R->press;
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

void face_turn_momentum_flux(struct fluxes *flux, struct geometry *geom)
{
  /* flux->momentum vector needs to be turned in case the HLLC or Rosunov Riemann solvers are used */

  double momx = flux->momentum[0];
  double momy = flux->momentum[1];
  double momz = flux->momentum[2];

  flux->momentum[0] = momx * geom->nx + momy * geom->mx + momz * geom->px;
  flux->momentum[1] = momx * geom->ny + momy * geom->my + momz * geom->py;
  flux->momentum[2] = momx * geom->nz + momy * geom->mz + momz * geom->pz;
}

void face_get_fluxes(struct state *st_L, struct state *st_R, struct state_face *st_face, struct fluxes *flux, struct geometry *geom, double *vel_face)
{                               //todo st_L, st_R not needed!
  double fac;

  /* calculate fluxes for ordinary Riemann solver */

  fac = (st_face->velx - vel_face[0]) * geom->nx + (st_face->vely - vel_face[1]) * geom->ny + (st_face->velz - vel_face[2]) * geom->nz;

  flux->mass = st_face->rho * fac;

  flux->momentum[0] = (st_face->rho * st_face->velx * fac + st_face->press * geom->nx);
  flux->momentum[1] = (st_face->rho * st_face->vely * fac + st_face->press * geom->ny);
  flux->momentum[2] = (st_face->rho * st_face->velz * fac + st_face->press * geom->nz);


  flux->energy = (0.5 * st_face->rho * (st_face->velx * st_face->velx +
                                        st_face->vely * st_face->vely +
                                        st_face->velz * st_face->velz) +
                  st_face->press / GAMMA_MINUS1) * fac + st_face->press * (st_face->velx * geom->nx + st_face->vely * geom->ny + st_face->velz * geom->nz);
}




void face_clear_fluxes(struct fluxes *flux)
{
  flux->mass = 0;
  flux->momentum[0] = 0;
  flux->momentum[1] = 0;
  flux->momentum[2] = 0;
  flux->energy = 0;
}

void face_add_fluxes_advection(struct state_face *st_face, struct fluxes *flux, struct geometry *geom, double *vel_face)
{
  double fac = -vel_face[0] * geom->nx - vel_face[1] * geom->ny - vel_face[2] * geom->nz;

  flux->mass += st_face->rho * fac;

  flux->momentum[0] += st_face->rho * st_face->velx * fac;
  flux->momentum[1] += st_face->rho * st_face->vely * fac;
  flux->momentum[2] += st_face->rho * st_face->velz * fac;

  flux->energy += 0.5 * st_face->rho * fac * (st_face->velx * st_face->velx + st_face->vely * st_face->vely + st_face->velz * st_face->velz) + st_face->press / GAMMA_MINUS1 * fac;
}


int flux_list_data_compare(const void *a, const void *b)
{
  if(((struct flux_list_data *) a)->task < (((struct flux_list_data *) b)->task))
    return -1;

  if(((struct flux_list_data *) a)->task > (((struct flux_list_data *) b)->task))
    return +1;

  return 0;
}



void dg_apply_flux_list(CBV(a))
{
  int i, j, p, nimport, ngrp, recvTask;

#ifdef DISCONTINUITY_DETECTION
  double density_temp;
#endif

  /* now exchange the flux-list and apply it when needed */

  mysort(FluxList, Nflux, sizeof(struct flux_list_data), flux_list_data_compare);

  TIMER_START(CPU_DG_FLEXCH);

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

  TIMER_STOP(CPU_DG_FLEXCH);

  /* apply the fluxes */

  int l, quadp;

  for(i = 0; i < nimport; i++)
    {
      p = FluxListGet[i].index;

#ifdef DISCONTINUITY_DETECTION
      density_temp = a[p][0][0];
#endif

      for(quadp = 0; quadp < Nof_outer_quad_points; quadp++)
        {
          for(l = 0; l < Nof_base_functions; l++)
            {
              a[p][l][0] =
                a[p][l][0] - FluxListGet[i].dM[quadp] / SphP[p].Volume * GET_Outer_base_values(FluxListGet[i].interface, quadp, l) * GET_Outer_quad_points_weights(FluxListGet[i].interface, quadp);
              a[p][l][1] =
                a[p][l][1] - FluxListGet[i].dP[quadp][0] / SphP[p].Volume * GET_Outer_base_values(FluxListGet[i].interface, quadp, l) * GET_Outer_quad_points_weights(FluxListGet[i].interface, quadp);
              a[p][l][2] =
                a[p][l][2] - FluxListGet[i].dP[quadp][1] / SphP[p].Volume * GET_Outer_base_values(FluxListGet[i].interface, quadp, l) * GET_Outer_quad_points_weights(FluxListGet[i].interface, quadp);
              a[p][l][3] =
                a[p][l][3] - FluxListGet[i].dP[quadp][2] / SphP[p].Volume * GET_Outer_base_values(FluxListGet[i].interface, quadp, l) * GET_Outer_quad_points_weights(FluxListGet[i].interface, quadp);
              a[p][l][4] =
                a[p][l][4] - FluxListGet[i].dEnergy[quadp] / SphP[p].Volume * GET_Outer_base_values(FluxListGet[i].interface, quadp, l) * GET_Outer_quad_points_weights(FluxListGet[i].interface,
                                                                                                                                                                        quadp);
            }
        }

#ifdef DISCONTINUITY_DETECTION
      if(a[p][0][0] > density_temp + dd_epsilon)        // inflow boundary
        {
          SphP[p].Inflow_boundaries += (1 << FluxListGet[i].interface);
        }
#endif

    }


  myfree(FluxListGet);
}

#endif //DG
