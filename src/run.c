/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/run.c
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
#include <unistd.h>
#include <ctype.h>

#include "allvars.h"
#include "proto.h"
#include "voronoi.h"
#include "domain.h"

#ifdef MRT
#include "MRT/RT_proto.h"
#include "MRT/RT.h"
#endif

static void do_second_order_source_terms_first_half(void);
static void do_second_order_source_terms_second_half(void);

/* TO DO:
 * - modify restart code to still support continuation of run beyond a new MaxTime from restart files.
 * - consolidate global variables into ones that need to be kept across a restart, and ones that don't.
 * - introduce separate timesteps for gravity and hydrodynamics
 */

/*! \file run.c
 *  \brief Contains the main simulation loop
 */

/*! \brief contains the main simulation loop that iterates over
 * single timesteps.
 *
 * The loop terminates when the cpu-time limit is
 * reached, when a `stop' file is found in the output directory, or
 * when the simulation ends because we arrived at TimeMax.
 *
 * If the simulation is started from initial conditions, a domain decomposition performed,
 * the gravitational forces are computed and the Voronoi mesh is constructed.
 *
 * The main loop is structured as follow:
 *  - find new timesteps: find_timesteps()
 *  - first gravitational half kick: do_gravity_step_first_half()
 *  - gradients are calculated: calculate_gradients()
 *  - vertex velocities are assigned: set_vertex_velocities()
 *  - computation of the hydro flux: compute_interface_fluxes()
 *  - (de)refinement of hydro cells: do_derefinements_and_refinements()
 *  - drifting particles to next sync point: find_next_sync_point()
 *  (Afterwards the timebins are updated, so different particles might now be active then before)
 *  - (if needed) a new domain decomposition: domain_Decomposition()
 *  - computation of gravitational forces: within do_gravity_step_second_half()
 *  - construction of the Voronoi mesh: create_mesh()
 *  - update of primitive variables: update_primitive_variables()
 *  - second gravitational half kick: do_gravity_step_second_half()
 *
 */
void run(void)
{
  CPU_Step[CPU_MISC] += measure_time();

#ifdef SIDM
  reset_timestep_counter();
#endif

  if(RestartFlag != 1)          /* if we have restarted from restart files, no need to do the setup sequence */
    {
      mark_active_timebins();

      output_log_messages();

      set_non_standard_physics_for_current_time();

      ngb_treefree();
      domain_free();
      domain_Decomposition();   /* do domain decomposition if needed */

      ngb_treeallocate();
      ngb_treebuild(NumGas);

      calculate_non_standard_physics_prior_mesh_construction();

      create_mesh();

      mesh_setup_exchange();

#ifdef MRT

      if(TimeBinSynchronized[All.HighestOccupiedTimeBin]) 
	add_source_fluxes() ;
#if defined(MRT_CHEMISTRY_PS2011) || defined(MRT_CHEMISTRY_PS2009)
      mrt_update_chemistry();  
#endif  
#ifdef MRT_IR_ONLY_CHEMISTRY
      mrt_IR_chemistry() ;
#endif

#endif


      update_primitive_variables();

      //#ifdef MRT_SETUP_SPECIAL_BOUNDARIES
	// mrt_setup() ;
      //#endif


      calculate_non_standard_physics_end_of_step();

      exchange_primitive_variables();

      calculate_gradients();

      set_vertex_velocities();  /* determine the speed of the mesh-generating vertices */

      ngb_update_velocities();  /* update the neighbor tree with the new vertex and cell velocities */


#ifdef RT_ADVECT
      rt_voronoi_exchange_primitive_variables();
      rt_calculate_green_gauss_gradients();
      rt_voronoi_exchange_primitive_variables_and_gradients();
#endif

#ifdef SECOND_DERIVATIVES
      exchange_primitive_variables_and_gradients();
      calculate_green_gauss_hessian(&Mesh);
#endif

      do_second_order_source_terms_second_half();

      do_gravity_step_second_half();
    }

#if defined(VORONOI_STATIC_MESH) || defined(AMR_STATIC_MESH)
  if(RestartFlag == 1)
    {
      int n_hydro_backup = TimeBinsHydro.NActiveParticles;
      int *time_bin_hydro = (int*)malloc( NumGas * sizeof(int) );
      int *hydro_particles = (int*)malloc( n_hydro_backup * sizeof(int) );
      for(int j=0; j<TimeBinsHydro.NActiveParticles; j++)  
       hydro_particles[j] = TimeBinsHydro.ActiveParticleList[j];
         
      for(int j=0; j<NumGas; j++)
       {
         time_bin_hydro[j] = P[j].TimeBinHydro;
         P[j].TimeBinHydro = All.HighestActiveTimeBin;
         TimeBinsHydro.ActiveParticleList[j] = j;         
       }
      TimeBinsHydro.NActiveParticles = NumGas;

      create_mesh();
      mesh_setup_exchange();

      for(int j=0; j<NumGas; j++)
        P[j].TimeBinHydro = time_bin_hydro[j];

      TimeBinsHydro.NActiveParticles = n_hydro_backup;
      for(int j=0; j<TimeBinsHydro.NActiveParticles; j++)
        TimeBinsHydro.ActiveParticleList[j] = hydro_particles[j];

      free( time_bin_hydro );
      free( hydro_particles );
    }
#endif

  while(1)                      /* main loop */
    {
      if(RestartFlag != 1)      /* if we are starting from restart files, skip in the first iteration the parts until the restart files were written  */
        {

#ifdef TRACER_MC_CHECKS
          check_tracer_lists();
#endif
          compute_statistics();

          flush_everything();

          create_snapshot_if_desired();

          if(All.Ti_Current >= TIMEBASE)        /* we reached the final time */
            {
              mpi_printf("\nFinal time=%g reached. Simulation ends.\n", All.TimeMax);

              if(All.Ti_lastoutput != All.Ti_Current)   /* make a snapshot at the final time in case none has produced at this time */
                produce_dump(); /* this will be overwritten if All.TimeMax is increased and the run is continued */

              break;
            }

          find_timesteps_without_gravity();     /* find-timesteps */

#ifdef SIDM
          reset_timestep_counter();
#endif

#ifdef NUCLEAR_NETWORK_DETONATE_CORE
          detonate_core();
#endif

#ifndef OTVET_NOGRAVITY
          find_gravity_timesteps_and_do_gravity_step_first_half();      /* gravity half-step for hydrodynamics */
          /* kicks collisionless particles by half a step */
#endif

#if (defined(SELFGRAVITY) || defined(EXTERNALGRAVITY) || defined(EXACT_GRAVITY_FOR_PARTICLE_TYPE) || defined(DG_EXTERNAL_ACCELERATION)) && !defined(MESHRELAX)
          update_timesteps_from_gravity();
#endif

#ifdef MHD_CT
          do_mhd_ct_source_terms();
#endif

          do_second_order_source_terms_first_half();

          exchange_primitive_variables();

          /* let's reconstruct gradients for every cell using Green-Gauss gradient estimation */
          calculate_gradients();

          /* determine the speed of the mesh-generating vertices */
          set_vertex_velocities();

          /* update the neighbor tree with the new vertex and cell velocities */
          ngb_update_velocities();


          exchange_primitive_variables_and_gradients();

#ifdef TRACER_PARTICLE
          tracer_particle_assign_cell_properties_and_timestep();
#endif

          /*and the second derivatives if necessary after Voronoi exchange */

#ifdef RT_ADVECT
          rt_voronoi_exchange_primitive_variables();
          rt_calculate_green_gauss_gradients();
          rt_voronoi_exchange_primitive_variables_and_gradients();
#endif

#ifdef SECOND_DERIVATIVES
          calculate_green_gauss_hessian(&Mesh);
          voronoi_exchange_hessians();
#endif

#ifdef RUNGE_KUTTA_FULL_UPDATE
          rk_save_conservative_variables();
#endif


#ifdef DG
          /* compute a discontinuous Galerkin step and update the conserved variables */
          dg_compute_step();
#else
	  //#ifdef MRT_SETUP_SPECIAL_BOUNDARIES
	    //mrt_setup() ;
	  //#endif
          /* compute intercell flux with Riemann solver and update the cells with the fluxes */
          compute_interface_fluxes(&Mesh);
#ifdef MRT
	  mrt_run();
#endif
#endif

#ifdef BECDM
          do_becdm_kinetic_drift();
#endif

#ifdef OPTIMIZE_MESH_MEMORY_FOR_REFINEMENT
#ifndef VORONOI_STATIC_MESH
          free_mesh_structures_not_needed_for_derefinement_refinement();
#endif
#endif

#ifdef REFINEMENT
          do_derefinements_and_refinements();
#endif

          write_cpu_log();      /* output some CPU usage log-info (accounts for everything needed up to completion of the current sync-point) */

          find_next_sync_point();       /* find next synchronization time */

          make_list_of_active_particles();

          output_log_messages();        /* write some info to log-files */

#if !defined(VORONOI_STATIC_MESH) && !defined(AMR)
#ifdef OPTIMIZE_MESH_MEMORY_FOR_REFINEMENT
          free_all_remaining_mesh_structures();
#else
          free_mesh();
#endif
#endif

          /* Check whether we should write a restart file.
           * Note that at this place we do not need to store the mesh, not the gravity tree.
           */
          if(check_for_interruption_of_run())
            return;
        }
      else
        RestartFlag = 0;

      set_non_standard_physics_for_current_time();

#if defined(VORONOI_STATIC_MESH) && !defined(VORONOI_STATIC_MESH_DO_DOMAIN_DECOMPOSITION)       /* may only be used if there is no gravity */
#ifndef MRT_SOURCES
      if(NumPart > NumGas)
        do_box_wrapping();
#endif
#ifdef TRACER_MC
      domain_resize_storage_tracer(0);
#endif
#else

#ifdef AMR
      if(amr_check_domain_decomposition())
#else
      if(All.HighestActiveTimeBin >= All.SmallestTimeBinWithDomainDecomposition)        /* only do this for sufficiently large steps */
#endif
        {
#ifdef VORONOI_STATIC_MESH
          free_mesh();
#endif

          ngb_treefree();
          domain_free();

          drift_all_particles();

          domain_Decomposition();       /* do new domain decomposition, will also make a new chained-list of synchronized particles */

          ngb_treeallocate();
          ngb_treebuild(NumGas);

#if defined(VORONOI_STATIC_MESH) || defined(AMR_STATIC_MESH)
          create_mesh();
          mesh_setup_exchange();
#endif
        }
#ifdef AMR
      {
        amr_update_nodes();
      }
#endif

#endif



#ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE
      special_particle_update_list();
#endif

#if (defined(CIRCUMSTELLAR_IRRADIATION) || defined(ALPHA_VISCOSITY) || defined(CIRCUMSTELLAR_REFINEMENTS)) && !defined (EXTERNALGRAVITY)
      source_particle_update_list();
#endif


      calculate_non_standard_physics_prior_mesh_construction();

#if !defined(VORONOI_STATIC_MESH) && !defined(AMR_STATIC_MESH)
      create_mesh();
      mesh_setup_exchange();
#endif

#if !defined(MUSCL_HANCOCK) && !defined(DG)
#ifdef RUNGE_KUTTA_FULL_UPDATE
#ifdef GENERAL_RELATIVITY
      update_primitive_variables();
      do_second_order_source_terms_second_half(); // new instead of doing it below, required for unsplit that src terms are know before flux
#endif
      rk_finish_step();
#else
      exchange_primitive_variables_and_gradients();
#ifdef MRT_SETUP_SPECIAL_BOUNDARIES
      mrt_setup() ;
#endif

      compute_interface_fluxes(&Mesh);
#ifdef MRT
      mrt_run();
#endif

#endif
#endif

#ifdef OUTPUT_CELL_SPIN
      add_spin_source_term_from_grid_movement();
#endif

#ifdef MRT

      if(TimeBinSynchronized[All.HighestOccupiedTimeBin]) 
	add_source_fluxes() ;
#if defined(MRT_CHEMISTRY_PS2011) || defined(MRT_CHEMISTRY_PS2009)
      mrt_update_chemistry();  
#endif  
#ifdef MRT_IR_ONLY_CHEMISTRY
      mrt_IR_chemistry() ;
#endif

#endif

      
      update_primitive_variables();     /* these effectively closes off the hydro step */
      /* the masses and positions are updated, let's get new forces and potentials */
      //#ifdef MRT_SETUP_SPECIAL_BOUNDARIES
      //mrt_setup() ;
      //#endif


#ifdef MHD_CT
      do_mhd_ct_source_terms();
      correct_ctr_b();
#endif

#ifndef GENERAL_RELATIVITY
      do_second_order_source_terms_second_half();
#endif

#ifndef OTVET_NOGRAVITY
      do_gravity_step_second_half();    /* this closes off the gravity half-step */
#endif

      /* do any extra physics, Strang-split (update both primitive and conserved variables as needed ) */
      calculate_non_standard_physics_end_of_step();

#ifdef FLD
      fld();
#endif
    }

  restart(0);                   /* write a restart file at final time - can be used to continue simulation beyond final time */

  write_cpu_log();              /* output final cpu measurements */
}

void do_second_order_source_terms_first_half(void)
{
#ifdef INSPIRAL
  do_inspiral_source_terms_first_half();
#endif
#ifdef RELAXOBJECT_BINARY
  do_binary_source_terms_first_half();
#endif
#ifdef MHD
  do_mhd_source_terms_first_half();
#endif
#ifdef DUST_LIVE
  do_drag_step_first_half();
#endif
#ifdef COSMIC_RAYS_DIFFUSION_EXPLICIT
  do_cr_diffusion();
#endif
#ifdef GENERAL_RELATIVITY
  do_general_relativity_source_terms();
#endif
#ifdef MRT_RADIATION_PRESSURE
  do_radiation_pressure_source_terms() ;
#endif
#ifdef MRT_COMOVING
  do_comoving_frame_source_terms() ;
#endif
}

void do_second_order_source_terms_second_half(void)
{
#ifdef GENERAL_RELATIVITY
  do_general_relativity_source_terms();
#endif
#ifdef COSMIC_RAYS_DIFFUSION_EXPLICIT
  do_cr_diffusion();
#endif
#ifdef DUST_LIVE
  do_drag_step_second_half();
#endif
#ifdef MHD
  do_mhd_source_terms_second_half();
#endif
#ifdef RELAXOBJECT_BINARY
  do_binary_source_terms_second_half();
#endif
#ifdef INSPIRAL
  do_inspiral_source_terms_second_half();
#endif
#ifdef MRT_RADIATION_PRESSURE
  do_radiation_pressure_source_terms() ;
#endif
#ifdef MRT_COMOVING
  do_comoving_frame_source_terms() ;
#endif
}

/*! \brief calls extra modules after drift operator
 *
 * This routine is called after the particles are drifted
 * to the next syncpoint, but before a new domain decomposition
 * is performed.
 */
void set_non_standard_physics_for_current_time(void)
{
#if defined(COOLING) && defined(GFM_COOLING_METAL)
  read_cooling_tables_current_time();
#endif
#if defined(WINDTUNNEL) && defined(WINDTUNNEL_EXTERNAL_SOURCE)
  interpolate_from_wind_table(All.Time, &All.InjectionDensity, &All.InjectionVelocity);
#endif

#if defined(COOLING) && !defined(SIMPLE_COOLING) && !defined(GRACKLE)
  IonizeParams();               /* set UV background for the current time */
#endif

//#ifdef SINK_PARTICLES
//  sink_particles();
//#endif

#ifdef NOH_PROBLEM
  set_special_noh_boundary_conditions();
#endif
}

/*! \brief calls extra modules after the gravitational force is recomputed
 *
 */
void calculate_non_standard_physics_with_valid_gravity_tree(void)
{
  /*  *** NOTICE *** if HIERARCHICAL_GRAVITY is adopted, this function is carried out once per synchronization time,
   *  with in general only a partial tree that does not necessarily contain all particles. The latter is the case only for steps
   *  where the highest timesteps are active ("full timesteps").
   */

#ifdef BLACK_HOLES
  blackhole_find_neighboring_holes_and_potmin();
  blackhole_do_mergers();
#ifdef REPOSITION_ON_POTMIN
  blackhole_reposition();
#endif
#ifdef BH_NEW_CENTERING
  if(TimeBinSynchronized[All.HighestOccupiedTimeBin])  /* do this only on global steps */
    blackhole_centering();
#endif
#ifdef BH_FRICTION
  blackhole_friction_apply();
#endif
#ifdef BH_HARMONIC_OSCILLATOR_FORCE
  blackhole_harmonic_force();
#endif
#endif

#ifdef ISM_LOCAL_RADIATION_PRESSURE
  ism_find_clump_centers();     /* loop over all (star forming) gas particles, set field SphP[i].ClumpCenter to clump center location                   */
  ism_find_clump_radii();
  ism_find_clump_interior_properties();
  ism_kick_particles();
#endif

/* Will try to place it elsewhere */
/*
#ifdef ISM_HII_HEATING
  ism_HII_heating();
#endif
*/

#ifdef OTVET
  do_otvet();
#endif
}

void calculate_non_standard_physics_with_valid_gravity_tree_always(void)
{
#ifdef SIDM
  double Ekin_before_sidm_parts, Ekin_after_sidm_parts;
  double Ekin_before_all_parts, Ekin_after_all_parts;
  Ekin_before_sidm_parts = sidm_calc_kinetic_energy_sidm_parts();
  Ekin_before_all_parts = sidm_calc_kinetic_energy_all_parts();
  sidm_DoScatter();
  Ekin_after_sidm_parts = sidm_calc_kinetic_energy_sidm_parts();
  Ekin_after_all_parts = sidm_calc_kinetic_energy_all_parts();
  All.sidm_EnergyInjectedCheck_sidm_parts += Ekin_after_sidm_parts - Ekin_before_sidm_parts;
  All.sidm_EnergyInjectedCheck_all_parts += Ekin_after_all_parts - Ekin_before_all_parts;

  mpi_printf("SIDM: injected energy (sidm particles) in step (double-check)  = %g (relative=%g)\n", Ekin_after_sidm_parts - Ekin_before_sidm_parts, fabs(Ekin_after_sidm_parts - Ekin_before_sidm_parts)/Ekin_after_sidm_parts);
  mpi_printf("SIDM: injected energy (sidm particles) in total (double-check) = %g (relative=%g)\n", All.sidm_EnergyInjectedCheck_sidm_parts, All.sidm_EnergyInjectedCheck_sidm_parts/Ekin_after_sidm_parts);
  mpi_printf("SIDM: injected energy (all particles) in step (double-check)  = %g  (relative=%g)\n", Ekin_after_all_parts - Ekin_before_all_parts, fabs(Ekin_after_all_parts - Ekin_before_all_parts)/Ekin_after_all_parts);
  mpi_printf("SIDM: injected energy (all particles) in total (double-check) = %g  (relative=%g)\n", All.sidm_EnergyInjectedCheck_all_parts, All.sidm_EnergyInjectedCheck_all_parts/Ekin_after_all_parts);
#endif

#if defined(DUST_LIVE) && defined(DL_GRAIN_BINS)
#if defined(DL_SNE_DESTRUCTION) || defined(DL_SHATTERING) || defined(DL_COAGULATION)
  do_shattering_coagulation();
#endif
#endif
}

/*! \brief calls extra modules before the Voronoi mesh is build
 *
 */
void calculate_non_standard_physics_prior_mesh_construction(void)
{
#ifdef BLACK_HOLES
  /* note: we estimate the gas density around the BH before the gravity tree is constructed such
   * that the later gravity tree construction can export current values for the BH_Hsml/BH_U, etc.
   * fields to the Tree_Points array, if needed. These are then later used in (gravity) tree-searches
   * for the local potential minimum and for other BHs nearby that are candidates for BH-mergers.
   */
  blackhole_density();
#ifdef BH_BIPOLAR_FEEDBACK
  blackhole_bipolar();
#endif
#ifdef BONDI_DISK_VORTICITY
  blackhole_disk_vorticity();
#endif
#endif

#if defined(FOF) && (defined(BLACK_HOLES) || defined(GFM_WINDS_VARIABLE) || defined(GFM_BIPOLAR_WINDS) || defined(GFM_WINDS_LOCAL))
  /* this will find new black hole seed halos and/or assign host halo masses for the variable wind model */
  if(All.Time >= All.TimeNextOnTheFlyFoF && TimeOfLastDomainConstruction == All.Time)
    {
#if (defined(GFM_WINDS_VARIABLE) && (GFM_WINDS_VARIABLE==1)) || defined(GFM_WINDS_LOCAL)
      fof_fof(-2);
#else
      fof_fof(-1);
#endif

      if(All.ComovingIntegrationOn)
        All.TimeNextOnTheFlyFoF *= All.TimeBetOnTheFlyFoF;
      else
        All.TimeNextOnTheFlyFoF += All.TimeBetOnTheFlyFoF;
    }
#endif



#ifdef GFM_STELLAR_EVOLUTION
  start_enrichment();
  find_cells_to_enrich();
  evolve_active_stars();
  do_chemical_enrichment();
#if defined(DUST_LIVE) && defined(DL_PRODUCTION)
  create_dust_particles();
#ifdef DL_REFINEMENT
  refine_dust_particles();
#endif
#endif
#ifdef GFM_DUST
  dust_growth_and_destruction();
#endif
#if defined(FM_RADIATION_FEEDBACK)
  find_radiation_feedback_cells();
  do_radiation_stellar_feedback();
#endif
#ifdef ISM_HII_HEATING
  ism_HII_heating();
#endif
#ifdef FM_EARLY_STAR_FEEDBACK
  do_early_stellar_feedback();
#endif

#if defined(GFM_STELLAR_FEEDBACK) || defined(GFM_WINDS_LOCAL) || defined(FM_STAR_FEEDBACK)
#if defined(FM_STAR_FEEDBACK) && defined(DELAYED_COOLING)
  find_feedback_cells();
#endif
  do_stellar_feedback();
#if defined(FM_STAR_FEEDBACK) && defined(OUTPUT_STELLAR_FEEDBACK)
  output_stellar_feedback_statistics();
#endif
#endif

  /* logs only for highest time bin */
  if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
    output_stellar_evolution_statistics();

  end_enrichment();
#endif

#if defined(GFM_WINDS) || defined(GFM_WINDS_LOCAL)
  do_winds();
#endif

#ifdef GFM_CHECKS
  check_AuxDataID_references();
#endif

#ifdef BLACK_HOLES
  blackhole_accretion();
#endif

#if defined(DUST_LIVE) && defined(DL_GRAIN_BINS)
#ifdef DL_SNE_DESTRUCTION
  update_sn_rates();
#endif
#if defined(DL_GROWTH) || defined(DL_SPUTTERING)
  do_dust_transfer();
#endif
#ifdef DL_SNE_DESTRUCTION
  do_dust_transfer_sne();
#endif
#endif

#ifdef GFM_WINDS_LOCAL
  create_winds_local();
#endif

#if defined(COOLING) && defined(USE_SFR) && !defined(LOCAL_FEEDBACK)
  sfr_create_star_particles();
#endif

#if defined(LOCAL_FEEDBACK) && defined(USE_SFR)
#if defined(EXTERNALSHEARBOX_KSRATE_RANDOM) || defined(EXTERNALSHEARBOX_MIXED_INJECTION)
  if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
    if(All.NumCurrentTiStep > 0)
      create_dummy_particles();   /* dummy particles are place holders for SN event positions; they have to be created, do their feedback, and be destroyed in the same timestep.  Currently feedback only happens at global sync points, so dummy particles can only be created at global sync points */
  
#endif
#ifndef EXTERNALSHEARBOX_KSRATE_RANDOM
  create_star_particles();      /* call this even if star particles are not activated to print information to log file */
#endif
#endif

#ifdef COSMIC_RAYS
#ifdef COSMIC_RAYS_SN_INJECTION
  deposit_cosmic_rays_from_supernovae();
#endif
#endif

#ifdef GFM_WINDS_STRIPPING
  start_stripping();            /* stellar_evolution_util.c */

  find_cells_to_strip();        /* stellar_density.c */

  strip_active_winds();         /* stellar_evolution_main.c */

  do_chemical_stripping();      /* stellar_evolution_main.c */

  end_stripping();              /* stellar_evolution_util.c */
#endif

#if defined(CIRCUMSTELLAR) && defined(CIRCUMSTELLAR_SINKS)
  circumstellar_swallow_gas();
#endif

#ifdef SINKS
  sinks();
#endif

#ifdef SINK_PARTICLES
  sink_particles();
#endif

#ifdef DVR_RENDER
#if (DVR_RENDER==1)
  dvr_render_main();
#endif
#endif

#ifdef GFM_DUST
#ifdef GFM_DUST_CAP
  gfm_check_dust(1);
#endif
#endif

#ifdef ACCRETE_ONTO_CENTRAL_POTENTIAL
  accrete_onto_central_potential();
#endif

#ifdef PERTURB_VELOCITIES
  perturb_velocities();
#endif

}


/*! \brief calls extra modules at the end of the run loop
 *
 * The second gravitational half kick is already applied to the
 * particles and the voronoi mesh is updated.
 *
 */
void calculate_non_standard_physics_end_of_step(void)
{
#ifdef SHOCK_FINDER_ON_THE_FLY
  shock_finder_on_the_fly();
#endif

#if defined(NUCLEAR_NETWORK) && !defined(NUCLEAR_NETWORK_TIMESTEP_LIMITER)
  network_main( NULL );
#endif

#if defined (VS_TURB) || defined (AB_TURB)
  update_primitive_variables();
  reset_turb_temp();
#if defined(POWERSPEC_GRID)
  if(All.Time >= All.TimeNextTurbSpectrum)
    {
      powerspec_turb(All.FileNumberTurbSpectrum, 0);
#ifdef TRACER_MC
      powerspec_turb(All.FileNumberTurbSpectrum, TRACER_MC);
#endif
      All.FileNumberTurbSpectrum++;
      All.TimeNextTurbSpectrum += All.TimeBetTurbSpectrum;
    }
#endif
#endif

#ifdef COSMIC_RAYS_STREAMING
  cosmic_rays_do_streaming();
#endif

#if defined(COSMIC_RAYS_DIFFUSION) && !defined(COSMIC_RAYS_DIFFUSION_EXPLICIT)
  do_cr_diffusion();
#endif

#ifdef TURBULENT_METALDIFFUSION
  //if(All.metaldiff_Ti_endstep == All.Ti_Current)
  if(1)
    {
      turbulent_metal_mixing();
      /*if(All.ComovingIntegrationOn)
        {
          if(All.Time>0.125)
	    turbulent_metal_mixing();
	  else
	    mpi_printf("TURBULENT_METALDIFFUSION: Nothing to do on this timestep: All.metaldiff_Ti_endstep, All.Ti_Current=%d,%d\n",All.metaldiff_Ti_endstep,All.Ti_Current);
	    }*/
    }
  else
    mpi_printf("TURBULENT_METALDIFFUSION: Nothing to do on this timestep: All.metaldiff_Ti_endstep, All.Ti_Current=%d,%d\n",All.metaldiff_Ti_endstep,All.Ti_Current);
#endif

#ifdef MONOTONE_CONDUCTION
  mpi_printf("CONDUCTION: EndStep = %llu\tCurrentStep = %llu\n", (long long) All.conduction_Ti_endstep, (long long) All.Ti_Current) ;
  if(All.conduction_Ti_endstep == All.Ti_Current)
    {
      if(All.ComovingIntegrationOn)
	{
	  if(All.Time>0.125)
	      monotone_conduction();
	  else
	    mpi_printf("CONDUCTION: Nothing to do, Very High Redhsift\n", All.Time);
	}
      else
	monotone_conduction();
    }
  else
    mpi_printf("CONDUCTION: Nothing to do on this timestep, begin=%llu, end=%llu, current=%llu\n", (long long) All.conduction_Ti_begstep, (long long) All.conduction_Ti_endstep, (long long) All.Ti_Current);
#endif

#ifdef IMPLICIT_OHMIC_DIFFUSION
  mpi_printf("NON-IDEAL MHD: Ohm diffusion, EndStep = %llu\tCurrentStep = %llu\n", (long long)All.ohmdiffusion_Ti_endstep, (long long)All.Ti_Current);
  if(All.ohmdiffusion_Ti_endstep == All.Ti_Current && All.Ti_Current > 0)
    ohmic_diffusion();
  else
    mpi_printf("NON-IDEAL MHD: Ohm diffusion, nothing to do on this timestep, begin=%llu, end=%llu, current=%llu\n", (long long)All.ohmdiffusion_Ti_begstep, (long long)All.ohmdiffusion_Ti_endstep, (long long)All.Ti_Current);
#endif

#if defined(BLACK_HOLES) && defined(GFM_AGN_RADIATION)
  if(All.HighestActiveTimeBin >= All.SmallestTimeBinWithAGNRad) /* only do this for sufficiently large steps */
    do_agn_radiation();
  CellsWithAGNRadiation = 0;
#endif


#ifdef ATOMIC_DM
  ADM_cooling();
#endif

#ifdef COOLING
#ifdef USE_SFR
#ifdef RADCOOL
  set_radcool_units_for_current_time();
#endif
#ifdef LOCAL_FEEDBACK
  cooling_only();
#else
  cooling_and_starformation();
#endif
#else
#ifdef RADCOOL
  set_radcool_units_for_current_time();
#endif
  cooling_only();
#endif

#if defined(BLACK_HOLES) && defined(GFM_AGN_RADIATION)
  agn_radiation_info();
#endif
#endif

#if defined(BLACK_HOLES) && defined(BH_ADIOS_WIND)
  blackhole_update_wind_affected_cells();
#endif

#ifdef COSMIC_RAYS
#ifdef COSMIC_RAYS_COOLING
  do_cosmic_ray_cooling();
#endif
#ifdef COSMIC_RAYS_ALFVEN_COOLING
  do_cosmic_ray_Alfven_cooling();
#endif
#endif

#if defined(TRACER_MC) && defined(TRACER_MC_NUM_FLUID_QUANTITIES)
  record_tracer_parent_fluid_properties();
#endif

#ifdef TGSET
  tgset_get_nh_max();
#endif

#ifdef HEALRAY
  healray();
#endif

#ifdef TGCHEM
  tgchem();
#endif

#ifdef SGCHEM
#ifndef ADVECTION_ONLY
  evolve_chemistry();
#endif
#endif

#ifdef BAROTROPIC
  mpi_printf("BAROTROPIC: Applying the Barotropic EOS \n");
  apply_barotropic_eos();
#endif

#if defined(LOCAL_FEEDBACK)
  if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
    {                           /* only do this for full steps */
      inject_feedback();
#ifdef EXTERNALSHEARBOX_KSRATE_UPDATE_PARAM
      update_shearbox_param();
#endif
    
#ifndef EXTERNALSHEARBOX_KSRATE_RANDOM
  compute_sfr();
#endif
    }
#endif



#ifdef RELAXOBJECT
  relaxobject();
#endif

#ifdef RT_ADVECT
  if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)    /* only do this for full steps */
    rt_advect_main();
#endif

#ifdef PREHEATING
  impose_preheating();
#endif

#ifdef GFM_PREENRICH
  if((All.PreEnrichTime >= 0) && (All.Time >= All.PreEnrichTime))
    {
      /* do the pre-enrichment */
      gfm_preenrich_gas();

      /* flag as done */
      All.PreEnrichTime = -1;
    }
#endif

#if defined(BLACK_HOLES) && defined(USE_SFR)
  blackhole_energy_log_info();
#endif

#ifdef SNE_FEEDBACK
  mpi_printf("Calling SNE routines\n");
  sne_feedback();
#endif
}


/*! \brief checks whether the run must interrupted
 *
 * The run is interrupted either if the stop file is present or,
 * if 85% of the CPU time are up. This routine also handles the
 * regular writing of restart files. The restart file is also
 * written if the restart file is present
 *
 * \return 1 if the run has to be interrupted, 0 otherwise
 */
int check_for_interruption_of_run(void)
{

  /* Check whether we need to interrupt the run */
  int stopflag = 0;
  if(ThisTask == 0)
    {
      FILE *fd;
      char stopfname[MAXLEN_PATH];

      sprintf(stopfname, "%sstop", All.OutputDir);
      if((fd = fopen(stopfname, "r")))  /* Is the stop-file present? If yes, interrupt the run. */
        {
          fclose(fd);
          printf("stop-file detected. stopping.\n");
          stopflag = 1;
          unlink(stopfname);
        }

      sprintf(stopfname, "%srestart", All.OutputDir);
      if((fd = fopen(stopfname, "r")))  /* Is the restart-file present? If yes, write a user-requested restart file. */
        {
          fclose(fd);
          printf("restart-file detected. writing restart files.\n");
          stopflag = 3;
          unlink(stopfname);
        }

      if(CPUThisRun > 0.85 * All.TimeLimitCPU)  /* are we running out of CPU-time ? If yes, interrupt run. */
        {
          printf("reaching time-limit. stopping.\n");
          stopflag = 2;
        }
    }

  MPI_Bcast(&stopflag, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if(stopflag)
    {
      restart(0);               /* write restart file */

      MPI_Barrier(MPI_COMM_WORLD);

      if(stopflag == 3)
        return 0;

      if(stopflag == 2 && ThisTask == 0)
        {
          FILE *fd;
          char contfname[MAXLEN_PATH];
          sprintf(contfname, "%scont", All.OutputDir);
          if((fd = fopen(contfname, "w")))
            fclose(fd);

          if(All.ResubmitOn)
            execute_resubmit_command();
        }
      return 1;
    }

  /* is it time to write a regular restart-file? (for security) */
  if(ThisTask == 0)
    {
      if((CPUThisRun - All.TimeLastRestartFile) >= All.CpuTimeBetRestartFile)
        {
          All.TimeLastRestartFile = CPUThisRun;
          stopflag = 3;
        }
      else
        stopflag = 0;
    }

  MPI_Bcast(&stopflag, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if(stopflag == 3)
    {
      restart(0);               /* write an occasional restart file */
      stopflag = 0;
    }
  return 0;
}






/*! \brief Returns the next output time that is equal or larger than
 *  ti_curr
 *
 *  \param ti_curr current simulation time
 */
integertime find_next_outputtime(integertime ti_curr)
{
  int i, iter = 0;
  integertime ti, ti_next;
  double next, time;

  DumpFlagNextSnap = 1;
  ti_next = -1;

  if(All.OutputListOn)
    {
      for(i = 0; i < All.OutputListLength; i++)
        {
          time = All.OutputListTimes[i];

          if(time >= All.TimeBegin && time <= All.TimeMax)
            {
              if(All.ComovingIntegrationOn)
                ti = (integertime) (log(time / All.TimeBegin) / All.Timebase_interval);
              else
                ti = (integertime) ((time - All.TimeBegin) / All.Timebase_interval);

#ifdef PROCESS_TIMES_OF_OUTPUTLIST
              /* first, determine maximum output interval based on All.MaxSizeTimestep */
              integertime timax = (integertime) (All.MaxSizeTimestep / All.Timebase_interval);

              /* make it a power 2 subdivision */
              integertime ti_min = TIMEBASE;
              while(ti_min > timax)
                ti_min >>= 1;
              timax = ti_min;

              double multiplier = ti / ((double) timax);

              /* now round this to the nearest multiple of timax */
              ti = ((integertime) (multiplier + 0.5)) * timax;
#endif
              if(ti >= ti_curr)
                {
                  if(ti_next == -1)
                    {
                      ti_next = ti;
                      DumpFlagNextSnap = All.OutputListFlag[i];
                    }

                  if(ti_next > ti)
                    {
                      ti_next = ti;
                      DumpFlagNextSnap = All.OutputListFlag[i];
                    }
                }
            }
        }
    }
  else
    {
#ifdef TGSET
      time = All.TimeBetSnapshot * All.TimeBegin;
#else
      if(All.ComovingIntegrationOn)
        {
          if(All.TimeBetSnapshot <= 1.0)
            terminate("TimeBetSnapshot > 1.0 required for your simulation.\n");
        }
      else
        {
          if(All.TimeBetSnapshot <= 0.0)
            terminate("TimeBetSnapshot > 0.0 required for your simulation.\n");
        }

      time = All.TimeOfFirstSnapshot;
#endif
      iter = 0;

      while(time < All.TimeBegin)
        {
          if(All.ComovingIntegrationOn)
            time *= All.TimeBetSnapshot;
          else
            time += All.TimeBetSnapshot;

          iter++;

          if(iter > 1000000)
            terminate("Can't determine next output time.\n");
        }

      while(time <= All.TimeMax)
        {
          if(All.ComovingIntegrationOn)
            ti = (integertime) (log(time / All.TimeBegin) / All.Timebase_interval);
          else
            ti = (integertime) ((time - All.TimeBegin) / All.Timebase_interval);

          if(ti >= ti_curr)
            {
              ti_next = ti;
              break;
            }

          if(All.ComovingIntegrationOn)
            time *= All.TimeBetSnapshot;
          else
            time += All.TimeBetSnapshot;

          iter++;

          if(iter > 1000000)
            terminate("Can't determine next output time.\n");
        }
    }

  if(ti_next == -1)
    {
      ti_next = 2 * TIMEBASE;   /* this will prevent any further output */

      mpi_printf("\nRUN: There is no valid time for a further snapshot file.\n");
    }
  else
    {
      if(All.ComovingIntegrationOn)
        next = All.TimeBegin * exp(ti_next * All.Timebase_interval);
      else
        next = All.TimeBegin + ti_next * All.Timebase_interval;

#ifdef TIMESTEP_OUTPUT_LIMIT
      mpi_printf("\nRUN: Limiting timestep to %g to fulfill output frequency", 0.1 * (next - All.Time));
      All.TimestepOutputLimit = 0.1 * (next - All.Time);
#endif

      mpi_printf("\nRUN: Setting next time for snapshot file to Time_next= %g  (DumpFlag=%d)\n\n", next, DumpFlagNextSnap);
    }

  return ti_next;
}



/*! \brief executes the resubmit command
 *
 */
void execute_resubmit_command(void)
{
  char buf[1000];
  sprintf(buf, "%s", All.ResubmitCommand);
#ifndef NOCALLSOFSYSTEM
  system(buf);
#endif
}
