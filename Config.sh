#!/bin/bash            # this line only there to enable syntax highlighting in this file

##################################################
#  Enable/Disable compile-time options as needed #
##################################################


#--------------------------------------- Basic operation mode of code
NTYPES=6                       # number of particle types
PERIODIC
#TWODIMS
#AXISYMMETRY                    # This is for axisymmetry in cylindrical coordinates (requires TWODIMS and a stationary mesh)
#ONEDIMS
#LONG_X=10.0
#LONG_Y=2.0
#LONG_Z=10.0
#REFLECTIVE_X=1 #=2            # if set to 2, the boundary is inflow/outflow
#REFLECTIVE_Y=1 #=2
#REFLECTIVE_Z=1 #=2

#COOLING
#UVB_SELF_SHIELDING            # gas is self-shielded from the cosmic background based on its density
#USE_SFR
#QUICK_LYALPHA                # turns dense and cold gas to stars immediately
#QUICK_LYALPHA_LATETIMEONLY   # cooling and star formation only after a certain time
#SFR_KEEP_CELLS
#GAMMA=1.4
#ISOTHERM_EQS
#USE_ENTROPY_FOR_COLD_FLOWS
#ENTROPY_MACH_THRESHOLD=1.1
#PREHEATING
#NOHYDRO

#----------------------------------------MPI/Threading Hybrid
#NUM_THREADS=4                           # use OpenMP, with the given number of threads per MPI task
IMPOSE_PINNING
#IMPOSE_PINNING_OVERRIDE_MODE
#GENERIC_ASYNC                           # enables asynchronous communication scheme

#--------------------------------------- Mesh Type
#AMR
VORONOI

#----------------------------------------- SR/GR
#SPECIAL_RELATIVITY
#SPECIAL_RELATIVITY_HLLC
#SR_HLLC_ZERO_COMPVEL
#GENERAL_RELATIVITY
#METRIC_TYPE=1
#ATMOSPHERE_GENERAL_RELATIVITY
#ADIABATIC_GENERAL_RELATIVITY=0

#----------------------------------------- MHD
#MHD
#MHD_CT
#MHD_CT_IC
#MHD_CT_PERTURB_POSITIONS
#MHD_DIVBCLEANING
#MHD_POWELL
#MHD_POWELL_LIMIT_TIMESTEP
#MHD_POWELL_SPLIT
#MHD_SEEDFIELD
#MHD_SEEDPSPEC
#MHD_THERMAL_ENERGY_SWITCH

#----------------------------------------- NON-IDEAL MHD
#NON_IDEAL_MHD
#OHMIC_DIFFUSION
#AMBIPOLAR_DIFFUSION
#IMPLICIT_OHMIC_DIFFUSION
#OHM_CRANK_NICHOLSON
#ONLY_OHMIC_DIFFUSION
#ONLY_AMBIPOLAR_DIFFUSION
#NON_IDEAL_MHD_EXPLICIT_LIMIT_TIMESTEP

#----------------------------------------- COSMIC RAYS
#DIFFUSION
#COSMIC_RAYS
#COSMIC_RAYS_EXTRA_DIAGNOSTICS
#COSMIC_RAYS_COOLING
#COSMIC_RAYS_ALFVEN_COOLING
#COSMIC_RAYS_STREAMING
#COSMIC_RAYS_STREAMING_GLOBAL_CHI
#COSMIC_RAYS_STREAMING_EXPLICIT
#COSMIC_RAYS_DIFFUSION
#COSMIC_RAYS_DIFFUSION_CONSTANT_TIMESTEP
#COSMIC_RAYS_DIFFUSION_GLOBAL_TIMESTEP
#COSMIC_RAYS_DIFFUSION_EXPLICIT
#COSMIC_RAYS_DIFFUSION_EXPLICIT_LIMITER
#COSMIC_RAYS_DIFFUSION_FULL_NORMAL_GRADIENT
#COSMIC_RAYS_DIFFUSION_ALWAYS_USE_PRECONDITIONER
#COSMIC_RAYS_DIFFUSION_ANISOTROPIC
#COSMIC_RAYS_DIFFUSION_BOUNDARY_X
#COSMIC_RAYS_DIFFUSION_BOUNDARY_Y
#COSMIC_RAYS_DIFFUSION_BOUNDARY_Z
#COSMIC_RAYS_DIFFUSION_LIMITER
#COSMIC_RAYS_DIFFUSION_OLD
#COSMIC_RAYS_SN_INJECTION
#COSMIC_RAYS_SHOCK_ACCELERATION
#COSMIC_RAYS_IN_ICS
#COSMIC_RAYS_MAGNETIC_OBLIQUITY
#OUTPUT_CR_PRESSURE_GRADIENT

#--------------------------------------- Riemann solver
#VARIABLE_GAMMA
#RIEMANN_HLL
#RIEMANN_HLLC
#RIEMANN_ROSUNOV
#RIEMANN_HLLD
#RIEMANN_GAMMA

#AMR_CONNECTIONS
#AMR_GRADIENTS
#AMR_REDUCE_DOMAIN_DECOMPOISTION

#--------------------------------------- Reconstruction
#TVD_SLOPE_LIMITER
#TVD_SLOPE_LIMITER_VANLEER
#TVD_SLOPE_LIMITER_SUPERBEE
#TVD_SLOPE_LIMITER_ALBADA
#TVD_SLOPE_LIMITER_MINBEE
#TVD_SLOPE_LIMITER_MINMOD
#TVD_SLOPE_LIMITER_MC
#GRADIENT_LIMITER_DUFFELL
#DISABLE_TIME_EXTRAPOLATION              # use only when you know exactly what you are doing. activating this option will make your results wrong but can tell you about the behaviour of your code
#DISABLE_SPATIAL_EXTRAPOLATION           # use only when you know exactly what you are doing. activating this option will make your results wrong but can tell you about the behaviour of your code
#NO_SCALAR_GRADIENTS                     # disables time and spatial extrapolation for passive scalar fields
#GRADIENTS_GREEN_GAUSS                   # original (now depreciated) gradient estimate, reduced hydro scheme to first order

#--------------------------------------- Mesh motion and regularization
#VORONOI_STATIC_MESH
#VORONOI_STATIC_MESH_DO_DOMAIN_DECOMPOSITION  # for VORONOI_STATIC_MESH force domain decomposition if there exist non-gas particles
REGULARIZE_MESH_CM_DRIFT
REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED
REGULARIZE_MESH_FACE_ANGLE
#REGULARIZE_MESH_LLOYD
#OUTPUT_MESH_FACE_ANGLE
#STICKY_POINTS_ON_REFLECTIVE_SURFACE     # if reflective boundaries are used, allows points to move only tangentially at boundary

#--------------------------------------- Time integration options
#FORCE_EQUAL_TIMESTEPS    # this chooses a variable but global timestep
TREE_BASED_TIMESTEPS     # non-local timestep criterion (take 'signal speed' into account)
#DECOUPLE_TIMESTEPS       # allows different timebins for gravity and hydro. use only WITHOUT FORCE_EQUAL_TIMESTEPS
#MUSCL_HANCOCK           # original (now depreciated) time integration scheme, only first order
#RUNGE_KUTTA_FULL_UPDATE
#PM_TIMESTEP_BASED_ON_TYPES=2+4      # select particle types that should be considered in setting the PM timestep
#NO_PMFORCE_IN_SHORT_RANGE_TIMESTEP  # if this is on, PM force is not included in short-range timestep criterion

#--------------------------------------- Direct Volume Raycast Rendering
#DVR_RENDER=0                             # direct volumetric raycasting, 0->stand-alone, 1->on-the-fly
#DVR_RENDER_SMOOTH                        # smooth field
#DVR_RENDER_ORTHOGONAL                    # orthogonal projection
#DVR_NUM_FIELDS=3                         # can be used to set output to a subset of fields, defaults to 13 (all fields), not allowed to be <3
#DVR_STAY_IN_BOX                          # do not ray-trace beyond simulation volume

#--------------------------------------- Image generation
#VORONOI_MESHOUTPUT                      # 2D and 3D mesh output
#VORONOI_IMAGES_FOREACHSNAPSHOT
#VORONOI_FREQUENT_IMAGES                 # creates images with frequency 'TimeBetweenImages' given in parameterfile, independent of snapshots
#VORONOI_FIELD_DUMP_PIXELS_X=1536
#VORONOI_FIELD_DUMP_PIXELS_Y=150
#VORONOI_VELOCITY_FIELD_2D
#VORONOI_FIELD_COMPENSATE_VX=4.0
#VORONOI_FIELD_COMPENSATE_VY=0
#VORONOI_NEW_IMAGE
#VORONOI_PROJ_TEMP                       #project T instead of u
#VORONOI_PROJ                            # do projection along any predefined direction
#VORONOI_MULTIPLE_PROJECTIONS            # do face-on and edge-on projections by swapping y and z axes
#VORONOI_NOGRADS                         # add an additional set of images where density gradients are not taken into account
#IMAGE_FOOTERS
#VORONOI_PROJ_TAU                        # Projection using opacities, can be used to compute photospheres for stellar calculations -> also use VORONOI_PROJ and OPACITIES
#VORONOI_PROJ_SUBSTEPS=10
#PROJ_WRITE_RAYS
#OPACITIES                               # include opacities for stellar material, tables have to be supplied

#--------------------------------------- Refinement and derefinement
REFINEMENT_SPLIT_CELLS
REFINEMENT_MERGE_CELLS
#REFINEMENT_MERGE_PAIRS
#REFINEMENT_VOLUME_LIMIT
#REFINEMENT_HIGH_RES_GAS
#REFINEMENT_AROUND_BH=0                    # spatial refinement scheme near BHs (0: default, 1: ignore cell shape constraints and always refine)
#DEREFINE_ONLY_DENSE_GAS
#NODEREFINE_BACKGROUND_GRID
#DEREFINE_GENTLY
#OPTIMIZE_MESH_MEMORY_FOR_REFINEMENT       # deletes the mesh structures not needed for refinement/derefinemet to lower the peak memory consumption
#REFINEMENT_AROUND_DM                      # refine around DM particles according to their softening length (useful for binary systems)
#GMC_REFINEMENT
#JEANS_DEREFINEMENT_DENSITY_THRESHOLD
#NO_TARGET_MASS_CONDITION
#DISC_REFINE_ONLY
#REFINE_ONLY_WITH_TRACER
#ROTATING_HIGHRES_REGION
#RAMP_REFINE
#SNE_RAMP_REFINE

#--------------------------------------- Mesh-relaxing or mesh-adding (this will not carry out a simulation)
#MESHRELAX                     # this keeps the mass constant and only regularizes the mesh
#MESHRELAX_DENSITY_IN_INPUT
#ADDBACKGROUNDGRID=16
#AMR_REMAP

#--------------------------------------- Gravity treatment
SELFGRAVITY                   # switch on for self-gravity
#HIERARCHICAL_GRAVITY         # use hierarchical splitting of the time integration of the gravity
#CELL_CENTER_GRAVITY          # uses geometric centers to calculate gravity of cells, only possible with HIERARCHICAL_GRAVITY
#NO_GAS_SELFGRAVITY            # switch off gas self-gravity in tree
#GRAVITY_NOT_PERIODIC          # if gravity is not to be treated periodically
#GRAVITY_TALLBOX               # special switch for making treating gravity in z-extended box, with x/y periodic, and z nonperiodic. LONG_Z may be used but must be an integer.
#ALLOW_DIRECT_SUMMATION
#DIRECT_SUMMATION_THRESHOLD=1000
#EXACT_GRAVITY_FOR_PARTICLE_TYPE=4 #N-squared fashion gravity for a small number of particles of the given type
#NO_SELFGRAVITY_TYPE=1         # exclude particle type from self-gravity (can be used with exact gravity)
#NO_GRAVITY_TYPE=1             # disable computation of gravity on particle type
#EXACT_GRAVITY_REACTION        # include reaction to other particle types when using exact gravity
#EXTERNALGRAVITY               # switch on for external potential
#EXTERNALGY=0.0
#EXTERNALDISKPOTENTIAL
#EXTERNALSHEARBOX
#EXTERNALSHEARBOX_KSRATE_RANDOM
#EXTERNALSHEARBOX_KSRATE_UPDATE_PARAM
#ENFORCE_JEANS_STABILITY_OF_CELLS_EEOS
ENFORCE_JEANS_STABILITY_OF_CELLS    # this imposes an adaptive floor for the temperature
#EVALPOTENTIAL                 # computes gravitational potential
#EXTERNALSHEETY
#COMPUTE_POTENTIAL_ENERGY
#RANDOMIZE_DOMAINCENTER
#ACCRETE_ONTO_CENTRAL_POTENTIAL # Allow mass to be accreted onto the central potential (needs CENTRAL_MASS_POTENTIAL)


#--------------------------------------- Gravity softening
#NSOFTTYPES=4                  # Number of different softening values to which particle types can be mapped.
#MULTIPLE_NODE_SOFTENING       # If a tree node is to be used which is softened, this is done with the softenings of its different mass components
#INDIVIDUAL_GRAVITY_SOFTENING=2+4  # bitmask with particle types where the softenig type should be chosen with that of parttype 1 as a reference type
#ADAPTIVE_HYDRO_SOFTENING
#NSOFTTYPES_HYDRO=64           # this is only relevant for ADAPTIVE_HYDRO_SOFTENING can can be set to override default value of 64


#--------------------------------------- TreePM Options
PMGRID=512
#ASMTH=1.25
#RCUT=6.0

#PLACEHIGHRESREGION=2
#ENLARGEREGION=1.1
#GRIDBOOST=2
#ONLY_PM

#FFT_COLUMN_BASED
#PM_ZOOM_OPTIMIZED

#--------------------------------------- Things that are always recommended
#AUTO_SWAP_ENDIAN_READIC                # Enables automatic ENDIAN swapping for reading ICs
CHUNKING                 # will calculated the gravity force in interleaved blocks. This can reduce imbalances in case multiple iterations due to insufficient buffer size need to be done


#---------------------------------------- Single/Double Precision
DOUBLEPRECISION=1
DOUBLEPRECISION_FFTW
#OUTPUT_IN_DOUBLEPRECISION                # snapshot files will be written in double precision
#INPUT_IN_DOUBLEPRECISION                 # initial conditions are in double precision
#OUTPUT_COORDINATES_IN_DOUBLEPRECISION    # will always output coordinates in double precision
#NGB_TREE_DOUBLEPRECISION                 # if this is enabled, double precision is used for the neighbor node extension


#---------------------------------------- On the fly FOF groupfinder
#FOF                                # enable FoF output
#FOF_PRIMARY_LINK_TYPES=2           # 2^type for the primary dark matter type
#FOF_SECONDARY_LINK_TYPES=1+16+32   # 2^type for the types linked to nearest primaries
#FOF_SECONDARY_LINK_TARGET_TYPES=   # should normally be set to a list of all dark matter types (in zoom runs), if not set defaults to FOF_PRIMARY_LINK_TYPES
#FOF_GROUP_MIN_LEN=32               # default is 32
#FOF_LINKLENGTH=0.16                # Linkinglength for FoF (default=0.2)
#FOF_FUZZ_SORT_BY_NEAREST_GROUP=0   # sort fuzz particles by nearest group and generate offset table in catalog (=1 writes nearest group number to snapshot)
#FOF_STOREIDS                       # store IDs in group/subfind catalogue, do not order particles in snapshot files by group order
#USE_AREPO_FOF_WITH_GADGET_FIX      # Needed in order to run FOF with Arepo on Gadget snapshot files, if gas is present and should be linked to the FOFs
#ADD_GROUP_PROPERTIES               # This can be used to calculate additional properties for an already existing group catalogue. These are then added as additional columns to the HDF5 group catalogues.
#ADD_MAGNETIC_GROUP_PROPERTIES

#---------------------------------------- Subfind
#SUBFIND                            # enables substructure finder
#SAVE_HSML_IN_SNAPSHOT              # stores hsml, density, and velocity dispersion values in the snapshot files

#SUBFIND_MEASURE_H2MASS             # special measuremenat option for mass in molecular hydrogen
#SUBFIND_CALC_MORE                  # calculates also the velocity dispersion in the local density estimate (this is automatically enabled by several other options, e.g. SAVE_HSML_IN_SNAPSHOT)
#SUBFIND_EXTENDED_PROPERTIES        # adds calculation of further quantities related to angular momentum in different components

#--------------------------------------- SFR/feedback model
#METALS
#MIN_METALLICITY_ON_STARTUP
#STELLARAGE

#SOFTEREQS
#MODIFIED_EOS
#SLOW_RELAX_TO_EOS
#STEEPER_SFR_FOR_STARBURST

#-------------------------------------- AGN stuff
#BLACK_HOLES               # enables Black-Holes (master switch)
#BH_THERMALFEEDBACK        # quasar-mode: couple a fraction of the BH luminosity into surrounding
#BH_THERMALFEEDBACK_ACC    # quasar-mode: bursty quasar-mode, accumulate thermal energy
#BH_NF_RADIO               # radio-mode model based on Nulsen & Fabian theory
#DRAINGAS=1                # non-stochastic smooth accretion (1: on, 2: on + cell rho, 3: on + gas drained from all cells within hsml)
#BH_EXACT_INTEGRATION      # integrates analytically mass accretion
#BH_BONDI_DEFAULT          # default Bondi prescription
#BH_BONDI_DENSITY          # Bondi -> density dependent
#BH_BONDI_CAPTURE          # Bondi -> capture mass dependent
#BH_BONDI_DISK_VORTICITY   # Bondi -> vorticity dependent
#BH_DO_NOT_PREVENT_MERGERS # When this is enabled, BHs can merge irrespective of their relative velocity
#BH_USE_GASVEL_IN_BONDI    # only when this is enabled, the surrounding gas velocity is used in addition to the sounds speed in the Bondi rate
#BH_USE_ALFVEN_SPEED_IN_BONDI  # when this is enabled the alfven speed is added to the gas sound speed in the Bondi rate and the total gas pressure around the BH includes the magnetic contribution when compared the the reference pressure in BH_PRESSURE_CRITERION (requires MHD)
#MASSIVE_SEEDS             # BH seeds assigned large dynamical mass, such that ideally no repositioning is needed anymore
#MASSIVE_SEEDS_MERGER
#BH_NEW_CENTERING          # an alternative to the BH_FRICTION and REPOSITION_ON_POTMIN switches
#BH_DRAG                   # Drag on black-holes due to accretion: current implementation simply double-accounts for the accretion momentum transfer, and should therefore not be used
#REPOSITION_ON_POTMIN      # repositions hole on potential minimum (requires EVALPOTENTIAL)
#BH_PRESSURE_CRITERION
#BH_RELATIVE_NGB_DEVIATION # Maximum NGB number deviation calculated relative to total number of neighbours
#OUTPUT_BLACK_HOLE_TIMESTEP #outputs the 3 time-steps for BH particles
#BH_FRICTION				# Estimates the local DM density around BH and applies a friction force to the relative velocity, meant as a replacement for REPOSITION_ON_POTMIN
#BH_FRICTION_AGGRESSIVE
#BH_HARMONIC_OSCILLATOR_FORCE

#-------------------------------------- AGN spin evolution and recoil merger kicks
#BH_RECOIL_KICK				# Includes the remnant recoil of a BH merger
#BH_SPIN_EVOLUTION			# When this is enabled, spin evolution of black holes is computed
#BH_SPIN_MODEL=0			# Spin model to be used: 0-Prolonged spin model. 1-Chaotic spin model. 2-Mass dependend model. 3-Self-gravity dependend model.

#-------------------------------------- other AGN stuff
#UNIFIED_FEEDBACK        # activates BH_THERMALFEEDBACK at high Mdot and BH_BUBBLES FEEDBACK al low Mdot (-->OBSOLETE: replaced by BH_NEW_RADIO)
#BH_BUBBLES              # calculate bubble energy directly from the black hole accretion rate (-->OBSOLETE: replaced by BH_NEW_RADIO)
#BH_MAGNETIC_BUBBLES     # inject part of the  bubble energy as magnetic energy
#BH_MAGNETIC_DIPOLAR_BUBBLES #inject part of the bubble energy as magnetic energy, field arranged as a dipole with random orientation
#BH_ADIOS_WIND
#BH_ADIOS_DENS_DEP_EFFICIANCY  # makes the radiative efficiency density dependend
#BH_ADIOS_WIND_WITH_QUASARTHRESHOLD  # use a threshold value ("qusarthrehold") of bondi-rate over Eddington rate to decide about quasar mode vs. adios wind
#BH_ADIOS_WIND_WITH_VARIABLE_QUASARTHRESHOLD  # scales the threshold with black hole mass (with a factor (M_BH/M_ref)^2, where M_ref = 10^8 Msun)
#BH_ADIOS_WIND_DIRECTIONAL  # puts in momentum preferentially along a random direction
#BH_ADIOS_RANDOMIZED        # inputs momentum along alternating random directions
#BH_ADIOS_ONLY_ABOVE_MINIMUM_DENSITY   # disable ADIOS wind if density around blackhole drops below a certain fraction of the star formation density
#BH_CONTINOUS_MODE_SWITCH # calculates fraction of thermal and mechanical feedback energy depending on eddington factor and mass (continously in both quantities)

#-------------------------------------- Black Hole Refinement and Bipolar Options
#REFINEMENT_AROUND_BH_FIXED
#SUPPRESS_SF_IN_REFINEMENT_REGION
#BH_BIPOLAR_FEEDBACK

#---------------------------------------- Passive Tracers
#TRACER_FIELD                        # passive scalar field which is advected in proportion to fluid mass fluxes

#TRACER_MC=3                         # Monte Carlo tracer particles: master switch (value specifies output parttype)
#GENERATE_TRACER_MC_IN_ICS           # add a fixed number (given in the parameter file) of MC tracers to each gas cell in ICs
#TRACER_MC_NUM_FLUID_QUANTITIES=13   # number of fluid quantities to be stored for MC tracers - must match the number in TRACER_MC_STORE_WHAT
#TRACER_MC_STORE_WHAT=1+2+4          # bit mask for quantities to store (see allvars.h for bitmask)
#TRACER_NO_RESET_EACH_SNAP           # do not set tracked fluid quantities to zero after writing each snapshot
#TRACER_MC_CHECKS                    # carries out frequent consistency checks

#TRACER_PARTICLE=2                   # Velocity Field tracer particles: master switch (value specified parttype)
#GENERATE_TRACER_PARTICLE_IN_ICS     # add tracer particles at positions of cell vertices in ICs
#TRACER_PART_NUM_FLUID_QUANTITIES=8  # number of fluid quantities to be stored for velocity tracers - must match the value given to TRACER_PART_STORE_WHAT
#TRACER_PART_STORE_WHAT=1+2+4        # bit mask for quantities to store (see allvars.h for bitmask)

#TRACER_TRAJECTORY
#TRACER_TRAJECTORY_GENERATE

#OUTPUT_MCTRNUM                      # write number of MC tracers in each cell to the output snapshots

#-------------------------------------------- Things for special behaviour
#READ_DM_AS_GAS
#NO_ID_UNIQUE_CHECK
#RUNNING_SAFETY_FILE            # if file './running' exists, do not start the run
#LOAD_TYPES=1+2+4+16+32
#READ_COORDINATES_IN_DOUBLE
#IDS_OFFSET=1                   # offset for gas particles if created from DM
#TILE_ICS
#COMBINETYPES                   # reads in the IC file types 4+5 as type 3 (useful for doing gas runs of Aquarius ICs)
#MULTIPLE_RESTARTS
#TOLERATE_WRITE_ERROR
#OPTIMIZE_MEMORY_USAGE          # optimize for memory, not for speed. Note: this is dangerous for high dynamic range simulations with mixed precision, since some position variables are singles instead of doubles
#SUBBOX_SNAPSHOTS
#PROCESS_TIMES_OF_OUTPUTLIST
#EXTENDED_GHOST_SEARCH          # This extends the ghost search to the full 3x3 domain instead of the principal domain
#DOUBLE_STENCIL                 # this will ensure that the boundary region of the local mesh is deep enough to have a valid double stencil for all local cells
#TETRA_INDEX_IN_FACE            # adds an index to each entry of VF[] and DC[] to one of the tetrahedra that share this edge
VORONOI_DYNAMIC_UPDATE          # keeps track of mesh connectivity, which speeds up mesh construction
#COFFEE_PROBLEM
#NOH_PROBLEM
#SHIFT_BY_HALF_BOX
#DISABLE_VELOCITY_CSND_SLOPE_LIMITING
NO_MPI_IN_PLACE
NO_ISEND_IRECV_IN_DOMAIN
FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
#USE_MPIALLTOALLV_IN_DOMAINDECOMP
#MPI_HYPERCUBE_ALLGATHERV       # some MPI-libraries may use quite a bit of internal storage for MPI_Allgatherv. This uses hypercubes instead as a work-around
#MPISENDRECV_CHECKSUM
#NOTREERND
#ENLARGE_DYNAMIC_RANGE_IN_TIME  # This extends the dynamic range of the integer timeline from 32 to 64 bit
#NOSTOP_WHEN_BELOW_MINTIMESTEP
#TIMESTEP_OUTPUT_LIMIT          # Limit timesteps to write snaps on time for output lists with huge range
#DO_NOT_CREATE_STAR_PARTICLES
#DMPIC                          # enable special image code for dark matter simulations
#ALLOWEXTRAPARAMS
#RADIATIVE_RATES                # used in non-equilibrium chemistry model
#FIX_SPH_PARTICLES_AT_IDENTICAL_COORDINATES  # this can be used to load SPH ICs that contain identical particle coordinates
#VEL_POWERSPEC                  # compiles in a code module that allows via restart-flag 7 the calculation of a gas velocity power spectrum of a snapshot
#VEL_POWERSPEC_BOX
#ADJ_BOX_POWERSPEC              # compiles in a code module that allows via restart-flag 7 the calculation of gas power spectra of a snapshot with an adjustable box (user defined center and size)
#DISABLE_OPTIMIZE_DOMAIN_MAPPING
#RECOMPUTE_POTENTIAL_IN_SNAPSHOT   # needed for postprocess option 18 that can be used to calculate potential values for a snapshot
#ACTIVATE_MINIMUM_OPENING_ANGLE   # this does not open tree nodes under the relative opening criterion any more if their opening angle has dropped below a minimum angle
#USE_DIRECT_IO_FOR_RESTARTS     # Try to use O_DIRECT for low-level read/write operations of restart files to circumvent the linux kernel page caching
#PERTURB_VELOCITIES             # continuously perturb velocities when running simulation

#CUDA                       # enables CUDA support in Arepo
#CUDA_INSTRUMENT            # This enables instrumentation support for the nvidia profiler
#USE_DSDE                   # try to use a dynamic sparse data exchange paradigm to get rid off sparse MPI_Alltoall patterns on large partitions
#USE_NBC_FOR_IBARRIER       # use the NBC library to implement non-blocking collectives (only relevant when USE_DSDE is used)

#HUGEPAGES                  # use huge pages for memory allocation, through hugetlbfs library
#DETAILEDTIMINGS            # creates individual timings entries for primary/secondary kernels to diagnose work-load balancing

#PERFORMANCE_TEST_SPARSE_MPI_ALLTOALL
#BITS_PER_DIMENSION=42      # Peano-Hilbert order
#OVERRIDE_PEANOGRID_WARNING


#--------------------------------------- Output/Input options
#UPDATE_GRADIENTS_FOR_OUTPUT
#REDUCE_FLUSH
#OUTPUT_REFBHCOUNTER
#OUTPUT_EVERY_STEP
#GODUNOV_STATS
#OUTPUT_CPU_CSV
#OUTPUT_TASK
#OUTPUT_TIMEBIN_HYDRO
#OUTPUT_PRESSURE_GRADIENT
#OUTPUT_DENSITY_GRADIENT
#OUTPUT_VELOCITY_GRADIENT
#OUTPUT_BFIELD_GRADIENT
#OUTPUT_VERTEX_VELOCITY
#OUTPUT_VERTEX_VELOCITY_DIVERGENCE
#OUTPUT_VOLUME
#OUTPUT_CENTER_OF_MASS
#OUTPUT_SURFACE_AREA
#OUTPUT_PRESSURE
#OUTPUTPOTENTIAL
#OUTPUTACCELERATION
#OUTPUTTIMESTEP
#OUTPUT_SOFTENINGS            # output particle softenings
#OUTPUTGRAVINTERACTIONS       # output gravitatational interactions (from the tree) of particles
HAVE_HDF5                     # needed when HDF5 I/O support is desired
#HDF5_FILTERS                  # activate snapshot compression and checksum for HDF5 output
#OUTPUT_XDMF                   #writes an .xmf file for each snapshot, which can be read by visit (with the hdf5 snapshot)
#OUTPUTCOOLRATE                # outputs cooling rate, and conduction rate if enabled
#OUTPUT_DIVVEL                 # output  velocity divergence
#OUTPUT_CURLVEL                 # output  velocity curl
#OUTPUT_COOLHEAT               # output actual energy loss/gain in cooling/heating routine
#OUTPUT_VORTICITY
#OUTPUT_CELL_SPIN
#MEASURE_DISSIPATION_RATE      # measures and outputs dissipation rate. Note: requires USE_ENTROPY_FOR_COLD_FLOWS, even though it will then always use the thermal energy update
#OUTPUT_MACHNUM                # output maximum mach number of a cell
#OUTPUT_TASK
#OUTPUT_ENTROPY
#OUTPUT_CSND

#--------------------------------------- Testing and Debugging options
DEBUG                         # enables core-dumps
#DEBUG_ENABLE_FPU_EXCEPTIONS   # tries to enable FPU exceptions
#RESTART_DEBUG
#VERBOSE                       # reports readjustments of buffer sizes
HOST_MEMORY_REPORTING         # reports after start-up the available system memory by analyzing /proc/meminfo
#VTUNE_INSTRUMENT
#FORCETEST=0.001               # calculates for given fraction of particles direct summation forces to check accuracy of tree force
#FORCETEST_TESTFORCELAW=1      # this enables a special test to measure the effective force law of the code, can be set to 1 or 2

#--------------------------------------- Static Disk Potential
#DISK_POTENTIAL
#DISK_MASS_M0=1.0
#DISK_SCALE_R0=1.0

#--------------------------------------- Static NFW Potential
#STATICNFW
#NFW_C=12
#NFW_M200=100.0
#NFW_Eps=0.01
#NFW_DARKFRACTION=0.87

#--------------------------------------- Static Isothermal Sphere Potential
#STATICISO
#ISO_M200=100.0
#ISO_R200=160.0
#ISO_Eps=0.1
#ISO_FRACTION=0.9

#--------------------------------------- Static Hernquist Potential
#STATICHQ
#HQ_M200=186.015773
#HQ_C=10.0
#HQ_DARKFRACTION=0.9

#--------------------------------------- Growing Disk Potential
#GROWING_DISK_POTENTIAL

#--------------------------------------- Dark energy
#DARKENERGY # Enables Dark Energy
#TIMEDEPDE  # read w(z) from a DE file
#RESCALEVINI # rescale v_ini in read_ic / read_ic_cluster
#EXTERNALHUBBLE # reads the hubble function from the DE file
#TIMEDEPGRAV # resacles H and G according to DE model
#DARKENERGY_DEBUG # enable writing of drift/kick table

#--------------------------------------- Glass making/ 2nd-order initial conditions / Initial conditions options
#SECOND_ORDER_ICS
#LONGIDS
#OFFSET_FOR_NON_CONTIGUOUS_IDS
#GENERATE_GAS_IN_ICS
#SPLIT_PARTICLE_TYPE=4+8
#NTYPES_ICS=6 # number of particle types in ICs, if not NTYPES (only works for 6, and non-HDF5 ICs!)

#-------------------------------------- Simple turbulence test
#VS_TURB
#POWERSPEC_GRID=128

#AB_TURB

#READ_LEGACY_ICS

#--------------------------------------- Degenerate Equation of State
#EOS_DEGENERATE
#EOS_COULOMB_CORRECTIONS
#EOS_COULOMB_CORRECTIONS_SMOOTH
#EOS_NSPECIES=3
#RELAXOBJECT
#RELAXOBJECT_COOLING
#RELAXOBJECT_COOLING2
#RELAXOBJECT_BINARY
#RELAX_RUNTIME
#INSPIRAL
#PASSIVE_SCALARS=3

#-------------------------------------- Nuclear Network
#NUCLEAR_NETWORK
#NETWORK_NSE
#NETWORK_PARDISO
#NETWORK_SCREENING
#REACLIB1
#NUCLEAR_NETWORK_DETONATE
#NUCLEAR_NETWORK_DETONATE_CORE
#NUCLEAR_NETWORK_DETONATE_POSITION
#NUCLEAR_NETWORK_TIMESTEP_LIMITER
#NUCLEAR_NETWORK_USE_SHOCKFINDER
#NUCLEAR_NETWORK_LIMIT_COMPOSITION_CHANGE

#--------------------------------------- OPAL Equation of State
#EOS_OPAL
#EOS_NSPECIES=3               # first species is X(H), which is relevant for the EOS

#-------------------------------------- Radiative transfer options
#RT_ENABLE                  #RT master switch
#RT_COOLING_PHOTOHEATING
#RT_ADVECT                   # enable advection of radiation field
#RT_CGMETHOD                   # enables CG method solution of the RT advection
#RT_SLOWLIGHT                # enable slow light approximation
#RT_N_DIR=2                 # track this number of locally brightest sources (one of them is diffuse field)
#RT_COMBINE_N_DIR_IN_OUTPUT  # writes only a single summed photon/photon-density field into output files
#RT_ALLOW_ABSORBING_CELLS    # if this is set, all cells with ID>=1000000000 will absorb radiation
#RT_SPREAD_SOURCE
#RT_STELLAR_SOURCES
#RT_HEALPIX_NSIDE=1          # if this is set, a discretization of the solid angle is used instead of brightest source selection
#RT_INCLUDE_HE
#SOURCE_PERIODIC

#DO_NOT_MOVE_GAS
#HYDROGEN_ONLY

#-------------------------------------- Calculate Hessian Matrix
#SECOND_DERIVATIVES
#SLOPE_LIMIT_HESSIANS
#RECONSTRUCT_GRADIENTS

#-------------------------------------- Navier-Stokes Terms
#GLOBAL_VISCOSITY 	   #needs dynamic and bulk coefficients
#USE_KINEMATIC_VISCOSITY    #needs only one input parameter
#ALPHA_VISCOSITY=2            #for accretion disks
#LOCAL_VISCOSITY=1          #=1 Sutherland viscosity/ =2 Spitzer viscosity
#THERMAL_CONDUCTION
#TRACER_DIFFUSION           #requires TRACER_FIELD switched on

#-------------------------------------- Circumstellar Disks
#CIRCUMSTELLAR              #Master switch
#CIRCUMSTELLAR_WBOUNDARIES
#CIRCUMSTELLAR_IRRADIATION
#CIRCUMSTELLAR_SINKS
#CIRCUMSTELLAR_PLANET_GROWTH     #Requires BLACK_HOLES turned on
#GRAVITY_FROM_STARS_PLANETS_ONLY #Requires EXTERNALGRAVITY turned on
#CENTRAL_MASS_POTENTIAL     #Point-mass potential
#BINARY_POTENTIAL     #Fixed star-planet circular orbit
#LOCALLY_ISOTHERM_DISK      #Isothermal Equation of state at each radii.

#-------------------------------------- Special Boundaries within domain
#SPECIAL_BOUNDARY          #Main Switch
#COAXIAL_BOUNDARIES              #e.g. Couette flow-type boundaries

#-------------------------------------- Windtunnel
#WINDTUNNEL
#WINDTUNNEL_COORD=0                    # sets the coordinate in which the wind blows (0,1,2 for x,y,z)
#WINDTUNNEL_EXTERNAL_SOURCE
#WINDTUNNEL_FIXVARIABLESININJECTIONREGION # enables a region with fixed properties
#WINDTUNNEL_REFINEMENT_VOLUME_LIMIT #Volume refinement option for windtunnel setup. REFINEMENT_VOLUME_LIMIT should also be enabled.


#--------------------------------------- Boundaries with optional inflow/outflow
#BOUNDARY_INFLOWOUTFLOW_MINID=10000000   # defines the ID range describing inflow/outflow nozzle of wind-tunnel
#BOUNDARY_INFLOWOUTFLOW_MAXID=20000000
#BOUNDARY_REFL_FLUIDSIDE_MINID=30000000  # defines the ID ranges describing a reflective boundary
#BOUNDARY_REFL_FLUIDSIDE_MAXID=30000000  # defines the ID ranges describing a reflective boundary
#BOUNDARY_REFL_SOLIDSIDE_MINID=40000000
#BOUNDARY_REFL_SOLIDSIDE_MAXID=40000000
#BOUNDARY_REFL_ACTS_AS_SOURCE            # makes the boundary act as a source (using the inner values)
#BOUNDARY_STICKY_MINID=50000000          # this-ID range specifies cells that will not be moved, and neighbors of these cells will only do mesh regularization motions
#BOUNDARY_STICKY_MAXID=60000000
#STICKYFLAGS
#OUTPUT_STICKYFLAGS

#-------------------------------------- GFM - Galaxy Formation Module
#GFM                                    #master switch
#GFM_STELLAR_EVOLUTION=0                #stellar evolution: 0->default, 1->no mass loss (beta value changes + MassMetallicity & MassMetals inconsistent internally with cell dynamical mass) 2->call only test routine
#GFM_CONST_IMF=1                        #0 for Chabrier (default), 1 for a pure power-law (requires parameter IMFslope, e.g. -2.35 for Salpeter)
#GFM_VARIABLE_IMF=0                     #0 for a pure power-law that depends on DM-veldisp
#GFM_PREENRICH                          #pre enrich gas at given redshift
#GFM_EXACT_NUMNGB                       #use direct neighbor count instead of kernel weighted neighbor count
#GFM_WINDS                              #decoupled ISM winds
#GFM_WINDS_VARIABLE=0                   #decoupled ISM winds: 0->scale winds with halo mass, requires FoF, 1->sigma winds
#GFM_WINDS_VARIABLE_HUBBLE              #add an additional H(z)^(-1/3) factor to the wind scaling, such that it scales with halo mass not halo velocity dispersion
#GFM_WINDS_HUBBLESCALING                #scale the wind energy fraction with the Hubble rate, limit the maximum to 1
#GFM_WINDS_MASSSCALING                  #scale the wind energy mass loading with halo mass (equivalent to scaling the wind energy fraction with halo virial radius)
#GFM_WIND_ENERGY_METAL_DEPENDENCE       #this can be used to decrease the wind energy for high metallicity (mimicking higher cooling losses)
#GFM_WIND_ENERGY_METAL_DEPENDENCE_TANH  #this selects an alternative functional form for the transition, requires GFM_WIND_ENERGY_METAL_DEPENDENCE
#GFM_WINDS_STRIPPING                    #wind metal stripping
#GFM_WINDS_THERMAL                      #not only give the wind kinetic energy but also thermal energy
#GFM_WINDS_THERMAL_NEWDEF               #with this switch, the thermal energy is specified as a fraction of the total energy
#GFM_BIPOLAR_WINDS=1                    #decoupled ISM winds: bipolar winds: 0->default, 1->relative to motion of FOF group, 3->parallel to spin of star-forming gas in halo
#GFM_WINDS_LOCAL                        #energy-driven decoupled local sigma winds
#GFM_STELLAR_FEEDBACK                   #local SNIa and AGB energy and momentum feedback
#GFM_PRIMORDIAL_RATES                   #updated coefficients for primordial chemistry and cooling
#GFM_COOLING_METAL                      #metal line cooling
#GFM_UVB_CORRECTIONS                    #reionization energy corrections
#GFM_AGN_RADIATION                      #cooling suppression/heating due to AGN radiation field (proximity effect)
#GFM_STELLAR_PHOTOMETRICS               #calculate stellar magnitudes for different filters based on GALAXEV/BC03
#GFM_OUTPUT_MASK=1+2+4+8+16+32+64+128   #which fields to output (see io_fields.c)
#GFM_CHECKS                             #this checks the consistency of the AuxDataID/PID indices of stars and black holes every timestep
#GFM_DISCARD_ENRICHMENT_GRADIENTS       #this disables the gradient extrapolation of the passively advected metallicity scalar variables
#GFM_NORMALIZED_METAL_ADVECTION         #this introduces an additional pseudo element for all untracked metals and normalizes the extrapolated abundance vectors to unity
#GFM_OUTPUT_BIRTH_POS                   #output BirthPos and BirthVel for all star particles
#GFM_CHEMTAGS                           #see documentation/modules_GFM_chemtags
#GFM_WINDS_SAVE_PARTTYPE=2              #save wind particles as separate particle type instead of mixed with 4 (stars)
#GFM_DISCRETE_ENRICHMENT                #allow stars to enrich nearby gas from stellar evolution only above some delta mass fraction threshold
#GFM_SPLITFE                            #see documentation/modules_GFM_chemtags
#GFM_SPLITFE_ADDINAGB                   #add in the AGB iron half-half on the two iron SNIa/SNII tags such that the sum of them should be equal to the total iron
#GFM_RPROCESS                           #see documentation/modules_GFM_chemtags, must have GFM_SPLITFE toggled as well
#GFM_LAMBDA                             #output all cooling rates

#-------------------------------------- Dust physics
#GFM_DUST                               #formation and evolution of dust, requires GFM_STELLAR_EVOLUTION
#GFM_DUST_DESTMODE=0                    #dust destruction mode: 0->default (uses supernova rate), 1->constant destruction timescale
#GFM_DUST_SPUTTERING=1                  #sputtering of dust grains by gas-phase metals: 0->first principles calculation, 1->using empirical timescale
#GFM_DUST_COOLING                       #high temperature dust cooling
#GFM_DUST_MRN                           #MRN grain size distribution; otherwise single grain size with size in mu specified in parameter file
#GFM_DUST_CAP                           #cap negative dust masses

#-------------------------------------- Live dust physics
#DUST_LIVE=3                            #turns on live dust particles, value specifies output parttype; parttype 3 is recommended
#DL_STOPPING_TIME_CORRECTION            #include higher-order corrections to stopping timescale for supersonic flow, makes analytic tests more difficult
#DL_DRAG_SEMI_IMPLICIT                  #make use of drag analytic solution for velocity updates, instead of requiring explicit drag timesteps
#DL_GRAIN_BINS=10                       #track grain size distribution information for dust particles, using the specified number of bins; requires cooling, star formation
#DL_GRAIN_BINS_PIECEWISE_LINEAR         #allow grain size distribution bins to be piecewise linear, not just piecewise constant
#DL_GROWTH                              #enable growth of dust mass by accumulating gas-phase metals
#DL_SPUTTERING                          #loss of dust mass due to thermal sputtering
#DL_SNE_DESTRUCTION                     #loss of dust mass due to supernova shocks
#DL_SHATTERING                          #grain size shattering
#DL_COAGULATION                         #grain size coagulation
#DL_SHATTERING_DETAILED_INTEGRALS       #do not use piecewise constant approximation for integrals used in shattering and coagulation calculations to compute grain collision cross sections, but use slower piecewise linear integrals
#DL_PRODUCTION                          #creation of dust particles from star particles
#DL_REFINEMENT                          #refinement of large dust particles

#---------------------------------- Modified Gas Cooling
#RADCOOL                                #Include the effects of local radiation fields on gas cooling rates, requires GFM, if PMGRID not defined then uncomment GRAVITY_NOT_PERIODIC, currently not compatible with PLACEHIGHRESREGION
#RADCOOL_HOTHALO                        #Include radiation field from HOT GAS, works only with RADCOOL option
#RADCOOL_HOTHALO_METAL_BOOST            #Include an additional boost factor to account for the additional luminosity emitted in emission lines

#-------------------------------------- FM - Star formation and feedback module
#FM_SFR                                #turns on star formation (needs USE_SFR)
#FM_STAR_FEEDBACK                      #turns on stellar feedback
#FM_STAR_FEEDBACK_KICK_TYPE=0          #direction of the velocity kick: 0->random, 1->radial, 2->bipolar
#NON_STOCHASTIC_MOMENTUM_FEEDBACK      #enables non-probabilistic momentum feedback
#INJECT_INTO_SINGLE_CELL               #feedback energy is injected only over the closest neighbour of a star particle (requires NON_STOCHASTIC_MOMENTUM_FEEDBACK)
#DIRECT_MOMENTUM_INJECTION_FEEDBACK    #SN momentum is injected non-stochastically but determined separately from kinetic energy
#OUTPUT_SF_PROBABILITY                 #enables output of the probability of transforming gas cell into star particle
#TEST_SFR                              #only calls the SF initialization and saves Kennicutt law (and the gas effective EOS if available)
#USE_POLYTROPIC_EQSTATE                #imposes a minimum temperature to star forming gas (through a polytropic equation of state)
#DELAYED_COOLING                       #turns on delayed cooling model (Stinson et al. 2006)
#SHUTOFFTIME_UPDATE=0                  #update of the cooling shutoff time: 0->max(actual,computed); 1->sum(actual,computed)
#DELAYED_COOLING_TURB                  #turns on delayed cooling model with turbulent enegy advection (Teyssier et al. 2012)
#INSTANTANEOUS_DEPOSITION              #inject the SN energy after a minimum stellar age in at most FeedbackInjectionEvents times
#EXPLICIT_COOLING                      #switch to a 2-nd order explicit method for cooling if (u^{n+1} - u_{n}) < tol * u^{n}
#COMPUTE_SFR_FROM_H2                   #links the SFR to the H2 gas fraction
#TEST_COOLING_METAL                    #call only cooling test routine (save cooling function with metal cooling for solar metallicity)
#OUTPUT_STELLAR_FEEDBACK               #outputs SNII number, feedback energy and mass released for stellar particles and log files for feedback (requires GFM_STELLAR_EVOLUTION)
#OUTPUT_MOLECULAR_FRACTION             #outputs the H2 gas fraction (requires COMPUTE_SFR_FROM_H2 switched on)
#OUTPUT_OPTICAL_DEPTH                  #outputs the gas optical depth (requires COMPUTE_SFR_FROM_H2 switched on)
#RADPRESS_OPT_THIN                     #adds radiative pressure in optically thin approximation. If GFM active only young stars are considered. Needs OTVET
#RADPRESS_OPT_THIN_LUMPERMASS          #source emits at a rate proportional to mass (IonizingLumPerSolarMass in parameterfile). Otherwise constant given by IonizingLumPerSolarMass
#RADPRESS_OPT_THICK                    #adds radiation pressure feedback using radiative transfer. Needs OTVET active
#FM_RADIATION_FEEDBACK                 #inputs momentum to gas particles within stromgren radius, keep cells at 10^4 K and prevents star formation in them.
#FM_RADIATION_FEEDBACK_DEBUG           #extra output fields for FM_RADIATION_FEEDBACK
#FM_EARLY_STAR_FEEDBACK                #inputs momentum to neighbor gas particles according to the luminosity emmited by stars
#FM_EARLY_STAR_FEEDBACK_KICK_TYPE=1    #direction of the early-feedback velocity kick: 0->random, 1->radial, 2->bipolar
#OUTPUT_EARLY_STELLAR_FEEDBACK         #enable output of logging info (such as total momentum injected) for early stellar feedback
#FM_MASS_WEIGHT_SN                     #feedback energy weighted by mass instead of volume
#FM_VAR_SN_EFF	                       #SN efficiency scales with neighboring gas metallicity
#FM_MOLEC_COOLING                      #approx extra molecular cooling contribution addition based on fit to GRACKLE cooling curves
#FM_SN_COOLING_RADIUS_BOOST            #returns momentum and energy to the ISM accounting for an unresolved energy conserving early ST blast wave phase
#FM_STOCHASTIC_HII_PHOTOIONIZATION     #Photoionization is carried out stochastically based on mass in stromgren sphere

#-------------------------------------- Conduction
#MONOTONE_CONDUCTION                   # Monotonicity Preserving anisotropic diffusion solver for thermal conduction
#CONDUCTION_ISOTROPIC                  # Isotropic conduction
#CONDUCTION_ANISOTROPIC                # Anisotropic Conduction
#CONDUCTION_CONSTANT                   # Set Conduction coefficient constant
#CONDUCTION_SATURATION                 # Saturation of Conduction coefficient at low densities
#IMPLICIT_TI                           # Implicit time integration , backwards euler scheme -no limitation of time step
#SEMI_IMPLICIT_TI                      # Semi Implicit Time Integration, stable upto ncfl=4
#RESTRICT_KAPPA                        # Set a maximum value for diffusivity, done in order to avoid very small timesteps
#MULTIPLE_TIMESTEP                     # A approximated implicit time integration scheme which works on multiple time steps
#NON_LINEAR_SLOPE_LIMITERS             # ADIITIONAL SLOPE LIMITERS TO LIMIT NUMERICAL DIFFUSION


#-------------------------------------- MRT
#MRT                           # Moment based RT - currently support M1 closure
#MRT_COMOVING                  # Solve the RT equations in the comoving reference frame
#PHOTDENS_IN_ICS               # Read and write photon density in ICs   
#MRT_TIME_EXTRAPOLATION        # Include the Runge Kutta time extrapolation - needed for second order convergence
#MRT_COOLING_HEATING           # Cooling and Heating 
#MRT_RADIATION_PRESSURE        # Include radition pressure as soucre term for the momentum conservation equation
#MRT_INCLUDE_HE                # Include Helium in heating and cooling
#MRT_LSF_GRADIENTS             # Use the least square fit gradient estimates
#MRT_RIEMANN_ROSUNOV           # Use the Rosunov (GLF) riemann solver
#MRT_RIEMANN_HLLE              # Use the Harten-Lax-van Leer flux function
#MRT_MULTI_FREQUENCY           # Multi Frequency radiative transfer
#MRT_CHEMISTRY_PS2009          # Petkova & Springel 2009 Chemistry
#MRT_CHEMISTRY_PS2011          # Petkova & Springel 2011 Chemistry
#MRT_NOCOLLISION_IONIZATION    # No Collisional ionisation
#MRT_SLOWLIGHT                 # Reduce speed of light
#MRT_DO_NOT_MOVE_GAS           # Do not move gas in MRT move
#MRT_CONSTANT_KAPPA            # Constant Kappa
#MRT_IR                        # Include IR radiative transfer - ala Rosdahl+15
#MRT_IR_ONLY_CHEMISTRY         # Only do absorption terms and not the cooling terms - assume there is no absorption of energy, only the flux is absorbed
#MRT_IR_LTE                    # Do gas heating and cooling according to single fluid gas-dust-IR radiation LTE assumption
#MRT_IR_LTE_SEMI_IMPLICIT      # Follows  Rosdahl+15 iterative approach. Can lead to large number of interative steps (sometimes infinite)
#MRT_IR_LTE_GSL                # Use the GSL provided Implicit Bulirsch-Stoer method of Bader and Deuflhard to solve the energy equation (Prefered Method) 
#MRT_IR_PHOTON_TRAPPING        # Trap unresolved IR photons - not compatible with gradient extrapolations (hardcoded - no need to turn off time and space extrapolations)
#MRT_IR_GRAIN_KAPPA            # Use grain opacities calculated from Draine & Lee 1984, Laor & Draine 1993, requires local dust-to-gas ratio from GFM_DUST
#MRT_NO_UV                     # Do not include UV RT (mainly for testing purposes)
#MRT_SETUP_SPECIAL_BOUNDARIES  # Setup special boundary conditions (mainly for testing purposes)
#MRT_SOURCES                   # Enable source treatment for GFM stellar particles and black holes
#MRT_STARS                     # Include ionizing photons from GFM stellar particles
#MRT_BH                        # Include photons from black hole particles

#-------------------------------------- OTVET IMPLEMENTATION
#OTVET                                 #Master switch
#OTVET_CHEMISTRY_PS2009                #uses Petkova&Springel 2009 original chemical network, with CGmethod solving transport+absorption
#OTVET_CHEMISTRY_PS2011                #uses Petkova&Springel 2011 chemical network, with CGmethod solving only transport
#OTVET_NOGRAVITY                       #builds the tree but does not apply the kicks. Needs improvement to not account optionally for self-gravity of the gas, instead of all gravity like now
#OTVET_NOTMOVEGAS                      #Does not move gas according to flow
#OTVET_OUTPUT_ET                       #outputs the Eddington Tensor for all cells
#EDDINGTON_TENSOR_STARS                #activate stars as sources of radiation
#OTVET_MODIFY_EDDINGTON_TENSOR         #fully anisotropic Eddington Tensor
#OTVET_FLUXLIMITER                     #activate flux limited diffusion, expresion from Petkova & Springel (2009)
#OTVET_CHANGEFLUXLIMITER               #change flux limiter formula to Levermore and Pomraning (1981)
#OTVET_FIXTIMESTEP                     #fixed timestep, activate together with FORCE_EQUAL_TIMESTEPS
#OTVET_COOLING_HEATING                 #activate cooling and heating following Petkova & Springel (2009)
#OTVET_MULTI_FREQUENCY                 #relaxes the assumption of monochromatic emission
#OTVET_INCLUDE_HE                      #follow also Helium
#OTVET_SCATTER_SOURCE                  #distributes luminosity in an SPH-way in otvet_Ngb_source neighbour gas cells
#OTVET_OUTPUT_SOURCEHSML               #if OTVET_SCATTER_SOURCE is active, this outputs the HSML and density of sources
#OTVET_CHECK_PHOTONCOUNT               #checks explicity photon conservation during OTVET transport. Only makes sense if OTVET_CHEMISTRY_PS2011 isactive.
#OTVET_NOCOLLISION_IONIZATION          #switch off collisional ionization
#OTVET_SILENT                          # --inactive--
#OTVET_MODIFY_EDDINGTON_TENSOR         # --inactive--


#-------------------------------------- TG's switches for primordial simulations
#TGSET                                 #some custom settings
#TGCHEM                                #primordial chemistry and cooling network (Greif 2014)
#TGCHEM_TEST                           #primordial chemistry and cooling network test (Greif 2014)
#HEALRAY                               #adaptive ray-tracing (Greif 2014)
#SINKS                                 #sink particles (under construction)


#-------------------------------------- SIDM - Self-Interacting DM
SIDM=2                                #activate and set types
SIDM_MAXWELLIAN			       #Maxwellian cross section
#SIDM_CONST_CROSS                      #constant cross section
#SIDM_STATES=2                         #number of DM states (for inelastic models)
#SIDM_REACTIONS=5                      #number of scatter reactions (for inelasitc models)
#SIDM_NO_SCATTER                       #DEBUG: no scattering at all
#SIDM_NO_TIMESTEP                      #DEBUG: do not change timestep
#SIDM_NO_KINEMATICS                    #DEBUG: do not change particle velocities, but still run through full scattering process
#SIDM_NO_NGB_SEL                       #DEBUG: take closest particle to scatter with; will select particles in wrong state for multiple states; scatter state check is then turned off
#SIDM_NO_MASSCHANGE                    #DEBUG: no mass change during inelastic scattering
#SIDM_NO_ENERGYCHANGE                  #DEBUG: no energy change during inelastic scattering   NOTE: to get fully inelastic behaviour turn on SIDM_NO_MASSCHANGE and SIDM_NO_ENERGYCHANGE

#-------------------------------------- ISM - new detailed ISM/stellar feedback model
#ISM                                   #master switch (needs at minimum also GFM and GFM_STELLAR_EVOLUTION)
#ISM_LOCAL_RADIATION_PRESSURE          #local stellar radiation pressure from star forming clumps (GMC)
#ISM_LONG_RANGE_RADIATION_PRESSURE     #tbd
#ISM_HII_PHOTO_HEATING                 #tbd
#ISM_H2_SFR                            #tbd
#ISM_OUTPUT_FIELDS=1+2+4


#-------------------------------------- On-the-fly shock finder

#SHOCK_FINDER_BEFORE_OUTPUT             #Use this flag if you want to run the shock finder before a snapshot dump, no additional flags or parameters needed.
#SHOCK_FINDER_ON_THE_FLY                #Run the shock finder at every local timestep, no additional flags or parameters needed.


#--------------------------------------- Post-processing shock finder, please read the instructions in shock_finder.h
#SHOCK_FINDER_POST_PROCESSING           #post-processing shock finder
#SHOCK_FINDER_AREPO			            #standard operating mode
#UNLIMITED_GRADIENTS		            #standard option
#ZONE_JUMP_P			                #standard option
#ZONE_JUMP_T                            #standard option
#SHOCK_DIR_GRAD_T                       #standard option
#SHOCK_JUMP_T                           #standard option
#SURFACE_SPHERE_APPROX                  #use this for 2d sims
#SURFACE_ANGLE_APPROX                   #use this for 3d sims
#RESET_WRONG_JUMPS                      #standard option
#RESET_WRONG_RHO_JUMPS                  #standard option
#RESET_WRONG_P_JUMPS                    #standard option
#OUTPUT_GRAVITY_FRACTION                #standard option
#SKIP_BORDER                            #for non-periodic boundaries of the snapshot/subbox


#--------------------------------------- atomic dark matter (in fluid approximation)
#ATOMIC_DM                              #master switch


#--------------------------------------- Binary stellar systems
#BINARYLOG
#SPECIAL_SOFTENINGS
#REDUCE_SOFTENINGS
#ID_RGCORE=1000000000
#DMLOWESTTIMEBIN
#DMFIXED

#-------------------------------------- FLD
#FLD

#FLD_CONES
#FLD_NCONES=12

#FLD_CONST_KAPPA
#FLD_MARSHAK

#FLD_HYPRE
#FLD_HYPRE_IJ1
#FLD_HYPRE_IJ2
#HYPRE_PCG

#FLD_MG
#FLD_MG_GS

#FLD_ANISOTROPIC_CIRCULAR
#FLD_NO_TEMP_UPDATE
#FLD_SILENT


#FLD_TEST_BOUNDARY
#FLD_UPPER_BOUNDARY_MINID=1
#FLD_UPPER_BOUNDARY_MAXID=1000000
#FLD_LOWER_BOUNDARY_MINID=5000000
#FLD_LOWER_BOUNDARY_MAXID=6000000



#-------------------------------------- Calculate Quantities In Post Process
#CALCULATE_QUANTITIES_IN_POSTPROCESS
#POWERSPECTRUM_IN_POSTPROCESSING
#POWERSPECTRUM_IN_POSTPROCESSING_ICS


#-------------------------------------- Discontinuous Galerkin (DG)
#DG                                     #master switch
#DG_SET_IC_FROM_AVERAGES                #loads ordinary non-DG ICs and sets initial weights according to density, velocity and internal energy
#DG_TEST_PROBLEM                        #initial conditions are created from file src/dg/test_problems.c
#DG_VERBOSE                            #additional output for debugging
#DG_DEBUG                              #run in debug mode, additional checks are active
#DEGREE_K=1                             #spatial degree of the scheme
#CALC_QUADRATURE_DATA                   #calculate the quadrature data instead of reading it from a table
#RIEMANN_HLLC                           #use the hllc riemann solver instead of the normal one


#MINMOD_B                              #use the slope limiter as a total variaton bounded limiter
#DISCONTINUITY_DETECTION               #limit only when a discontinuity is found
#OUTPUT_DG_DISCONTINUITIES
#OUTPUT_DG_INFLOW_BOUNDARIES
#ANGLE_BOUND                           #use the angle bound methdod
#CHARACTERISTIC_LIMITER                #limit the characteristic variables instead of the conserved variables
#CONSERVED_LIMITER                     #limit the conserved variables
#POSITIVITY_LIMITER                     #keep the cell average values positive
#FIX_MEAN_VALUES                        #reset negative values

#-------------------------------------  Spiral Potential as used by Dobbs
#SPIRAL

#------------------------------------   Supernova Energy or Momentum cons.
#SNE_FEEDBACK
#CLUSTERED_SNE
#INJECT_TRACER_INTO_SN

#-----------------Deprecated
#CONDUCTION

#------------------------------------- Grackle (DO NOT USE, UNTESTED, SEE README)
#GRACKLE                               #master switch
#GRACKLE_H2                            #Turn on H2 cooling and chemistry
#GRACKLE_D                             #Turn on Deuterium cooling and chemistry
#GRACKLE_TAB                           #Run in tabulated mode
#GRACKLE_ABUNDANCE_IN_ICS              #Use abundances in ICs instead of converging on startup
#GRACKLE_VERBOSE

#------------------------------------- Modified gravity solver (Private to Ewald Puchwein, Volker Springel and Christian Arnold)
#MODGRAV                               #master switch
#MODGRAV_EFF_MASS                      #use the effective mass scheme to obtain the forces (currently the only method implemented)
#MODGRAV_INTERPOLATE_PHI
#BAROTROPIC
#PRIMCHEM

# SINK partices                                                                                                                                                    
#DUMP_SINK_PARTICLE_INFO
#SINK_PARTICLES
#SINK_PARTICLES_REFINEMENT_LIMIT

#TURBULENT_METALDIFFUSION
#TURBULENT_METALDIFFUSION_EXPLICIT

#------------------------------------- external Galaxy potential
#GALPOT                               

#------------------------------------- SGChem chemistry module
#SGCHEM_VARIABLE_Z                     #Allow metallicity and dust-to-gas ratio to vary between different cells
#SGCHEM_VARIABLE_ISRF                  #Allow interstellar radiation field strength to vary spatially
#SGCHEM_VARIABLE_CRION                 #Allow cosmic ray ionization rate to vary spatially
