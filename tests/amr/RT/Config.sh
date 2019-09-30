#!/bin/bash            # this line only there to enable syntax highlighting in this file

##################################################
#  Enable/Disable compile-time options as needed #
##################################################

#--------------------------------------- Basic operation mode of code
NTYPES=6                       # number of particle types
PERIODIC
TWODIMS
#AXISYMMETRY                    # This is for axisymmetry in cylindrical coordinates (requires TWODIMS and a stationary mesh)
#ONEDIMS
#LONG_X=10.0
LONG_Y=4.0
#LONG_Z=10.0
#REFLECTIVE_X=1 #=2            # if set to 2, the boundary is inflow/outflow */
REFLECTIVE_Y=1 #=2
#REFLECTIVE_Z=1 #=2
#MHD
#MHD_DIVBCLEANING
#MHD_POWELL
#MHD_SEEDFIELD
#COOLING
#UVB_SELF_SHIELDING            # gas is self-shielded from the cosmic background based on its density
#USE_SFR
#SFR_KEEP_CELLS
GAMMA=1.4
#ISOTHERM_EQS
#USE_ENTROPY_FOR_COLD_FLOWS
#ENTROPY_MACH_THRESHOLD=1.1
#PREHEATING

#----------------------------------------MPI/Threading Hybrid
#NUM_THREADS=4                           # use OpenMP, with the given number of threads per MPI task
#THREAD_COSTS_EXACT                      # measure gravity cost exactly for all threads instead of extrapolating from thread 0
#THREAD_COSTS_IMPROVE_MEM_AFFINITY       # special memory allocation for thread gravity cost measurement - only relevant if more threads than cores per socket are used 
#IMPOSE_PINNING
#IMPOSE_PINNING_OVERRIDE_MODE

#--------------------------------------- Mesh Type
AMR
#VORONOI

AMR_GRADIENTS
#AMR_CONNECTIONS

#--------------------------------------- Riemann solver
#VARIABLE_GAMMA
#RIEMANN_HLL
#RIEMANN_HLLC
#RIEMANN_ROSUNOV
#RIEMANN_HLLD
#RIEMANN_GAMMA

#--------------------------------------- Reconstruction
#TVD_SLOPE_LIMITER
#TVD_SLOPE_LIMITER_VANLEER
#TVD_SLOPE_LIMITER_SUPERBEE
#TVD_SLOPE_LIMITER_ALBADA
#TVD_SLOPE_LIMITER_MINBEE
#DISABLE_TIME_EXTRAPOLATION              # use only when you know exactly what you are doing. activating this option will make your results wrong but can tell you about the behaviour of your code
#DISABLE_SPATIAL_EXTRAPOLATION           # use only when you know exactly what you are doing. activating this option will make your results wrong but can tell you about the behaviour of your code
#LIMITER_ENFORCE_TVD                     # prevents the spatial reconstruction from introducing any new extrema
#NO_SCALAR_GRADIENTS                     # disables time and spatial extrapolation for passive scalar fields
#GRADIENTS_LEAST_SQUARE_FIT



#--------------------------------------- Mesh motion and regularization
#VORONOI_STATIC_MESH
#VORONOI_STATIC_MESH_DO_DOMAIN_DECOMPOSITION  # if run with VORONOI_STATIC_MESH and there exist non-gas particles, domain decomposition can be forced with this. however, whether it is a good choice depends on the setup, as the mesh construction after the domain decomposition is slow for VORONOI_STATIC_MESH
#REGULARIZE_MESH_CM_DRIFT
#REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED
#REGULARIZE_MESH_FACE_ANGLE
#REGULARIZE_MESH_OPTIMAL
#OUTPUT_MESH_FACE_ANGLE
#CALCULATE_VERTEX_VELOCITY_DIVERGENCE
#STICKY_POINTS_ON_REFLECTIVE_SURFACE     # if reflective boundaries are used, allows points to move only tangentially at boundary

#--------------------------------------- Time integration options
FORCE_EQUAL_TIMESTEPS    # this chooses a variable but global timestep
#TREE_BASED_TIMESTEPS     # non-local timestep criterion (take 'signal speed' into account)


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

#--------------------------------------- Mesh-relaxing or mesh-adding (this will not carry out a simulation)
#MESHRELAX                     # this keeps the mass constant and only regularizes the mesh
#MESHRELAX_DENSITY_IN_INPUT
#ADDBACKGROUNDGRID=16
#AMR_REMAP

#--------------------------------------- Gravity treatment
#SELFGRAVITY                   # switch on for self-gravity     
NO_GAS_SELFGRAVITY            # switch off gas self-gravity in tree 
#GRAVITY_NOT_PERIODIC          # if gravity is not to be treated periodically
#EXACT_GRAVITY_FOR_PARTICLE_TYPE=4 #N-squared fashion gravity for a small number of particles of the given type
EXTERNALGRAVITY               # switch on for external potential
EXTERNALGY=-0.1
#EXTERNALDISKPOTENTIAL
#EXTERNALSHEARBOX
#FIXED_GRAVITATIONAL_SOFTENINGS_FOR_CELLS
#ENFORCE_JEANS_STABILITY_OF_CELLS_EEOS
#ENFORCE_JEANS_STABILITY_OF_CELLS    # this imposes an adaptive floor for the temperature
#EVALPOTENTIAL                 # computes gravitational potential
#EXTERNALSHEETY
#COMPUTE_POTENTIAL_ENERGY
#ALLOW_NODE_SOFTENING
#LOCALIZED_SOFTENINGS

#--------------------------------------- TreePM Options
#PMGRID=512
#ASMTH=1.25
#RCUT=6.0

#PLACEHIGHRESREGION=2
#ENLARGEREGION=1.1
#GRIDBOOST=2
#ONLY_PM
#PMPERIODIC_LOWMEM_THREADED       # This replaces the standard pm_periodic.c with a version that uses less peak memory and is slightly faster
                                  # but this routine only works well for homogeneously sampled boxes. It can also be used with threads, in this case both OPENMP NUM_THREADS should be set
#NEW_FFT


#--------------------------------------- Things that are always recommended
#AUTO_SWAP_ENDIAN_READIC        # Enables automatic ENDIAN swapping for reading ICs
#PEANOHILBERT_EXTEND_DYNAMIC_RANGE
CHUNKING                 # will calculated the gravity force in interleaved blocks. This can reduce imbalances in case multiple iterations due to insufficient buffer size need to be done
                          

#---------------------------------------- Single/Double Precision
DOUBLEPRECISION=1
DOUBLEPRECISION_FFTW
OUTPUT_IN_DOUBLEPRECISION # snapshot files will be written in double precision
INPUT_IN_DOUBLEPRECISION  # initial conditions are in double precision

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

#---------------------------------------- Subfind
#SUBFIND                            # enables substructure finder
#SAVE_HSML_IN_SNAPSHOT              # this will store hsml and density values in the snapshot files 

#SUBFIND_MEASURE_H2MASS             # special measuremenat option for mass in molecular hydrogen
  
#--------------------------------------- SFR/feedback model
#METALS
#MIN_METALLICITY_ON_STARTUP
#STELLARAGE

#SOFTEREQS
#SLOW_RELAX_TO_EOS

#-------------------------------------- AGN stuff
#BLACK_HOLES               # enables Black-Holes (master switch)
#BH_THERMALFEEDBACK        # quasar-mode: couple a fraction of the BH luminosity into surrounding
#BH_THERMALFEEDBACK_ACC    # quasar-mode: bursty quasar-mode, accumulate thermal energy 
#BH_NF_RADIO               # radio-mode model based on Nulsen & Fabian theory
#DRAINGAS=1                # non-stochastic smooth accretion (1: on, 2: on + cell rho)
#BH_EXACT_INTEGRATION      # integrates analytically mass accretion
#BH_BONDI_DEFAULT          # default Bondi prescription
#BH_BONDI_DENSITY          # Bondi -> density dependent
#BH_BONDI_CAPTURE          # Bondi -> capture mass dependent
#BH_BONDI_DISK_VORTICITY   # Bondi -> vorticity dependent
#BH_DO_NOT_PREVENT_MERGERS # When this is enabled, BHs can merge irrespective of their relative velocity 
#BH_USE_GASVEL_IN_BONDI    # only when this is enabled, the surrounding gas velocity is used in addition to the sounds speed in the Bondi rate
#MASSIVE_SEEDS             # BH seeds assigned large dynamical mass, such that ideally no repositioning is needed anymore
#MASSIVE_SEEDS_MERGER
#BH_DRAG                   # Drag on black-holes due to accretion: current implementation simply double-accounts for the accretion momentum transfer, and should therefore not be used
#REPOSITION_ON_POTMIN      # repositions hole on potential minimum (requires EVALPOTENTIAL)
#REPOSITION_ON_POTMIN_VEL  # set new BH velocity
#BH_PRESSURE_CRITERION
#BH_RELATIVE_NGB_DEVIATION # Maximum NGB number deviation calculated relative to total number of neighbours
#OUTPUT_BLACK_HOLE_TIMESTEP #outputs the 3 time-steps for BH particles

#-------------------------------------- old/unused AGN stuff
#UNIFIED_FEEDBACK        # activates BH_THERMALFEEDBACK at high Mdot and BH_BUBBLES FEEDBACK al low Mdot (-->OBSOLETE: replaced by BH_NEW_RADIO)
#BH_BUBBLES              # calculate bubble energy directly from the black hole accretion rate (-->OBSOLETE: replaced by BH_NEW_RADIO)
#SWALLOWGAS              # Enables stochastic accretion of gas particles consistent with growth rate of hole (-->OBSOLETE: replaced by DRAINGAS)

#---------------------------------------- Passive Tracers
#TRACER_FIELD
#TRACER_PARTICLE=2                 # advect massless tracer particles of type TRACER_PARTICLE with velocity field
#TRACER_MC=1 #=2                   # Monte Carlo tracer particles (=1 to enable, >=2 to output as that partType)
#GENERATE_TRACER_PARTICLE_IN_ICS   # add tracer particles at positions of cell vertices in ICs
#GENERATE_TRACER_MC_IN_ICS         # add a fixed number (given in the parameter file) of MC tracers to each gas cell in ICs

#TRACER_PART_NUM_FLUID_QUANTITIES=8        # number of fluid quantities to be stored for velocity tracers - must match the value given to TRACER_PART_STORE_WHAT
#TRACER_PART_STORE_WHAT=1+2+4+8+16+32+64+128   # bit mask for quantities to store (see allvars.h for bitmask)

#TRACER_MC_NUM_FLUID_QUANTITIES=13                  # number of fluid quantities to be stored for MC tracers - must match the value gien to TRACER_MC_STORE_WHAT
#TRACER_MC_STORE_WHAT=1+2+4+8+16+32+64+128+256+512+1024+2048+4096     # bit mask for quantities to store (see allvars.h for bitmask)

#TRACER_MC_SKIPLOAD=3                   # skip reading this particle type when reading initial conditions from a snapshot file
#FOF_DISABLE_SNAP_REWRITE               # do not rewrite a snapshot file when RestartFlag==3


#-------------------------------------------- Things for special behaviour
#READ_DM_AS_GAS
#NO_ID_UNIQUE_CHECK
#RUNNING_SAFETY_FILE                # if file './running' exists, do not start the run
#LOAD_TYPES=1+2+4+16+32
#READ_COORDINATES_IN_DOUBLE 
#IDS_OFFSET=1       #offset for gas particles if created from DM
#TILE_ICS
#COMBINETYPES        # reads in the IC file types 4+5 as type 3 (useful for doing gas runs of Aquarius ICs)
#USE_RANDOM_GENERATOR
#MULTIPLE_RESTARTS
#TOLERATE_WRITE_ERROR
#OPTIMIZE_MEMORY_USAGE                   #optimize for memory, not for speed. Note: this is dangerous for high dynamic range simulations with mixed precision, since some position variables are singles instead of doubles
#SUBBOX_SNAPSHOTS
#INDIVIDUAL_GRAVITY_SOFTENING=4+8   # sets (via bitmask) the particle types for which the softening will be determined according to their mass instead of according to the parameter file
#PROCESS_TIMES_OF_OUTPUTLIST
#EXTENDED_GHOST_SEARCH           # This extends the ghost search to the full 3x3 domain instead of the principal domain
#ALTERNATIVE_GHOST_SEARCH        # This switches on the "old" routines that find the ghost neighbours

#DOUBLE_STENCIL                 # this will ensure that the boundary region of the local mesh is deep enough to have a valid double stencil for all local cells

#TETRA_INDEX_IN_FACE            # adds an index to each entry of VF[] and DC[] to one of the tetrahedra that share this edge 

#VORONOI_DYNAMIC_UPDATE          # keeps track of mesh connectivity, which speeds up mesh construction
#COFFEE_PROBLEM
#NOH_PROBLEM
#SHIFT_BY_HALF_BOX
#DISABLE_VELOCITY_CSND_SLOPE_LIMITING
#NO_FANCY_MPI_CONSTRUCTS
#NO_MPI_IN_PLACE
#NO_ISEND_IRECV_IN_DOMAIN
#FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
#USE_MPIALLTOALLV_IN_DOMAINDECOMP
#MPISENDRECV_SIZELIMIT=100
#MPI_HYPERCUBE_ALLGATHERV     # some MPI-libraries may use quite a bit of internal storage for MPI_Allgatherv. This uses hypercubes instead as a work-around
#MPISENDRECV_CHECKSUM
#NOTREERND
#ENLARGE_DYNAMIC_RANGE_IN_TIME   # This extends the dynamic range of the integer timeline from 32 to 64 bit
#NOSTOP_WHEN_BELOW_MINTIMESTEP
#NOTYPEPREFIX_FFTW
#DO_NOT_CREATE_STAR_PARTICLES
#DMPIC                          # enable special image code for dark matter simulations   
#ALLOWEXTRAPARAMS
#RADIATIVE_RATES               # used in non-equilibrium chemistry model
#FIX_SPH_PARTICLES_AT_IDENTICAL_COORDINATES  # this can be used to load SPH ICs that contain identical particle coordinates
#VEL_POWERSPEC                 # compiles in a code module that allows via restart-flag 7 the calculation of a gas velocity power spectrum of a snapshot
#ADJ_BOX_POWERSPEC         # compiles in a code module that allows via restart-flag 7 the calculation of gas power spectra of a snapshot with an adjustable box (user defined center and size)
#ALTERNATIVE_FFT
#DISABLE_OPTIMIZE_DOMAIN_MAPPING

#CUDA                       # enables CUDA support in Arepo
#CUDA_INSTRUMENT            # This enables instrumentation support for the nvidia profiler
#GPU_TREE                   # uses the GPU to calculate the tree grav interactions
#GPU_TREE_CALC_CPU          # does the computation on the CPU instead (for debugging)
#GPU_TREE_VERBOSE           # output more informations

#GPU_PM                     # uses the GPU for the FFTs of the PM force

#--------------------------------------- Output/Input options
#UPDATE_GRADIENTS_FOR_OUTPUT
#REDUCE_FLUSH
#OUTPUT_REFBHCOUNTER             
#OUTPUT_EVERY_STEP
#GODUNOV_STATS
#OUTPUT_CPU_CSV
OUTPUT_PRESSURE_GRADIENT
OUTPUT_DENSITY_GRADIENT
OUTPUT_VELOCITY_GRADIENT
#OUTPUT_VERTEX_VELOCITY
#OUTPUT_VERTEX_VELOCITY_DIVERGENCE  # requires CALCULATE_VERTEX_VELOCITY_DIVERGENCE
#OUTPUT_CENTER_OF_MASS
OUTPUT_SURFACE_AREA
OUTPUT_PRESSURE
#OUTPUTPOTENTIAL
#OUTPUTACCELERATION
#OUTPUTTIMESTEP
#OUTPUT_SOFTENINGS            # output particle softenings
#OUTPUTGRAVINTERACTIONS       # output gravitatational interactions (from the tree) of particles
HAVE_HDF5                     # needed when HDF5 I/O support is desired
#PARAMS_IN_SNAP                # add the compiler flags and parameter file values to every snapshot file (requires HAVE_HDF5)
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

#--------------------------------------- Testing and Debugging options
DEBUG                         # enables core-dumps 
#DEBUG_ENABLE_FPU_EXCEPTIONS   # tries to enable FPU exceptions
#CHECKSUM_DEBUG
#RESTART_DEBUG
#VERBOSE                       # reports readjustments of buffer sizes
#HOST_MEMORY_REPORTING         # reports after start-up the available system memory by analyzing /proc/meminfo
#FORCETEST=0.001               # calculates for given fraction of particles direct summation forces to check accuracy of tree force

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

#-------------------------------------- Simple turbulence test
#VS_TURB
#POWERSPEC_GRID=128
#GAMMA=1.01

#AB_TURB

#--------------------------------------- Planets/Materials (Robert)
#MATERIALS                   # master switch for all extension for planetary physics 
#NUMBER_OF_MATERIALS=1
#READ_LEGACY_ICS

#--------------------------------------- Degenerate Equation of State
#EOS_DEGENERATE
#EOS_COULOMB_CORRECTIONS
#EOS_NSPECIES=3
#RELAXOBJECT
#PASSIVE_SCALARS=3

#-------------------------------------- Nuclear Network
#NUCLEAR_NETWORK
#NETWORK_NSE
#NETWORK_PARDISO
#NETWORK_SCREENING
#REACLIB1

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


#-------------------------------------- GFM - Galaxy Formation Module
#GFM                                    #master switch
#GFM_STELLAR_EVOLUTION=0                #stellar evolution: 0->default, 1->no mass loss (beta value changes + MassMetallicity & MassMetals inconsistent internally with cell dynamical mass) 2->call only test routine 
#GFM_CONST_IMF=1                        #0 for Chabrier (default), 1 for a pure power-law (requires parameter IMFslope, e.g. -2.35 for Salpeter)
#GFM_VARIABLE_IMF=0                     #0 for a pure power-law that depends on DM-veldisp
#GFM_PREENRICH                          #pre enrich gas at given redshift
#GFM_WINDS                              #decoupled ISM winds 
#GFM_WINDS_VARIABLE=0                   #decoupled ISM winds: 0->scale winds with halo mass, requires FoF, 1->sigma winds
#GFM_WINDS_HUBBLESCALING                #scale the wind energy fraction with the Hubble rate, limit the maximum to 1
#GFM_WINDS_MASSSCALING                  #scale the wind energy mass loading with halo mass (equivalent to scaling the wind energy fraction with halo virial radius)
#GFM_WINDS_STRIPPING                    #wind metal stripping
#GFM_WINDS_THERMAL                      #not only give the wind kinetic energy but also thermal energy
#GFM_BIPOLAR_WINDS=1                    #decoupled ISM winds: bipolar winds: 0->default, 1->relative to motion of FOF group, 3->parallel to spin of star-forming gas in halo
#GFM_WINDS_LOCAL                        #energy-driven decoupled local sigma winds
#GFM_STELLAR_FEEDBACK                   #local SNIa and AGB energy and momentum feedback
#GFM_PRIMORDIAL_RATES                   #updated coefficients for primordial chemistry and cooling
#GFM_COOLING_METAL                      #metal line cooling
#GFM_UVB_CORRECTIONS                    #reionization energy corrections               
#GFM_AGN_RADIATION                      #cooling suppression/heating due to AGN radiation field (proximity effect)
#GFM_STELLAR_PHOTOMETRICS               #calculate stellar magnitudes for different filters based on GALAXEV/BC03 
#GFM_OUTPUT_MASK=1+2+4+8+16+32+64+128   #which fields to output (search GFM_OUTPUT_MASK in io.c to see which fields the bits encode)
#GFM_DUST                               #production of dust from stellar ejecta, requires GFM_STELLAR_EVOLUTION

#-------------------------------------- FM - Star formation and feedback module
#FM_SFR                                #turns on star formation (needs USE_SFR)
#FM_STAR_FEEDBACK                      #turns on stellar feedback
#FM_STAR_FEEDBACK_KICK_TYPE=0          #direction of the velocity kick: 0->random, 1->radial, 2->bipolar
#NON_STOCHASTIC_MOMENTUM_FEEDBACK      #enables non-probabilistic momentum feedback
#INJECT_INTO_SINGLE_CELL               #feedback energy is injected only over the closest neighbour of a star particle (requires NON_STOCHASTIC_MOMENTUM_FEEDBACK)
#OUTPUT_SF_PROBABILITY                 #enables output of the probability of transforming gas cell into star particle
#TEST_SFR                              #only calls the SF initialization and saves Kennicutt law (and the gas effective EOS if available)
#USE_POLYTROPIC_EQSTATE                #imposes a minimum temperature to star forming gas (through a polytropic equation of state)
#DELAYED_COOLING                       #turns on delayed cooling model (Stinson et al. 2006)
#SHUTOFFTIME_UPDATE=0                  #update of the cooling shutoff time: 0->max(actual,computed); 1->sum(actual,computed)
#DELAYED_COOLING_TURB                  #turns on delayed cooling model with turbulent enegy advection (Teyssier et al. 2012)
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

#-------------------------------------- Conduction
#CONDUCTION
#CONDUCTION_CONSTANT
#CONDUCTION_SATURATION

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
#HEALRAY                               #adaptive ray-tracing (Greif 2014)
#SINKS                                 #sink particles (under construction)


#-------------------------------------- SIDM - Self-Interacting DM
#SIDM=2                                #activate and set types 
#SIDM_CONST_CROSS                      #constant transfer cross section
#SIDM_ADM                              #Francis-Yan ADM

#-------------------------------------- ISM - new detailed ISM/stellar feedback model 
#ISM                                   #master switch (needs at minimum also GFM and GFM_STELLAR_EVOLUTION)
#ISM_LOCAL_RADIATION_PRESSURE          #local stellar radiation pressure from star forming clumps (GMC)
#ISM_LONG_RANGE_RADIATION_PRESSURE     #tbd
#ISM_HII_PHOTO_HEATING                 #tbd
#ISM_H2_SFR                            #tbd
#ISM_OUTPUT_FIELDS=1+2+4                     


#--------------------------------------- Post-processing shock finder, please read the instructions in shock_finder.h
#SHOCK_FINDER				#master switch
#SHOCK_FINDER_AREPO			#standard operating mode 
#UNLIMITED_GRADIENTS			
#ZONE_JUMP_P				
#ZONE_JUMP_T
#SHOCK_DIR_GRAD_T
#SHOCK_JUMP_T
#SURFACE_ANGLE_APPROX
#RESET_WRONG_JUMPS
#SKIP_BORDER				#for non-periodic boundaries of the snapshot/subbox

                                 
#--------------------------------------- atomic dark matter (in fluid approximation)
#ATOMIC_DM                              #master switch 

#--------------------------------------- Star index in binary stellar systems
#STAR_INDEX
