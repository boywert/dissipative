
Config.sh Options
=================

A number of features of AREPO are controlled with **compile-time options** in the makefile rather 
than by the parameter file. This has been done in order to allow the generation of highly optimised 
binaries by the compiler, even when the underlying source code allows for many different ways to 
run the code. The ``Template-Config.sh`` file contains a list of all available compile-time options, 
with most of them commented out by default. Below, a brief guide to these options is given.

.. note::

  This page documents ``Template-Config.sh`` options which are not considered part of any 
  non-standard module. Documentation for module specific compile-time options can be found 
  under their respective pages.

.. contents::
  :local:
  


Basic operating modes
---------------------

NTYPES=6
  The number of distinct particle types used. No longer hard-coded to six.

PERIODIC
  Has to be always enabled for legacy reasons. Use REFLECTIVE_X/Y/Z to set
  non-periodic boundary conditions for the simulation domain.

REFLECTIVE_X/Y/Z=1
  Sets non-periodic boundary conditions for the X/Y/Z boundary. A value of 1
  makes a boundary a reflective boundary. A value of 2 sets inflow/outflow
  boundaries.

LONG_X/Y/Z=2.
  These options can be used to distort the periodic simulation cube along the given direction
  with the given factor into a parallelepiped of arbitrary aspect ratio. The box size in the
  given direction increases from the value in the parameterfile by the factor given (e.g. if
  Boxsize is set to 100. and LONG_X=4. is set the simulation domain extends from 0 to 400 along 
  X and from 0 to 100 along Y and Z.)

TWODIMS
  This effectively switches off one dimension, i.e. the code follows only 2d hydrodynamics in the 
  xy-plane. This only works with ``NOGRAVITY``, and if all coordinates of the third axis are exactly equal. 
  Can be useful for idealized tests.

ONEDIMS
  This switches to 1D. The only valid spatial coordinate is the first dimension. All values of the
  other two coordinates have to be zero.

ONEDIMS_SPHERICAL
  This can be used in addition to ONEDIMS to obtain 1D spherical coordinates. The first dimension is
  used as the radial coordinate.

Gravity treatment
-----------------

GRAVITY_NOT_PERIODIC
  Controls whether gravity is periodic or not independently from ``PERIODIC``. Further, ``PMGRID`` 
  can be used in non-periodic mode as well for the gravity.

RANDOMIZE_DOMAINCENTER
  For non-periodic gravity simulations: decreases time-correlated gravity force errors and hence 
  minimize the center-of-mass drift of isolated objects.


.. _gravity-details:

Tree gravity & softening
^^^^^^^^^^^^^^^^^^^^^^^^

In the default configuration, the code uses a small table of possible gravitational softening 
lengths, which are specified in the parameterfile through the ``SofteningComovingTypeX`` and 
``SofteningMaxPhysTypeX`` options, where X is an integer that gives the "softening type". Each 
particle type is mapped to one of these softening types through the ``SofteningTypeOfPartTypeY`` 
parameters, where ``Y`` gives the particle type. The number of particle types and the number of 
softening types do not necessarily have to be equal (although the default for both of them is 6). 
Several particle types can be mapped to the same softening if desired.

ADAPTIVE_HYDRO_SOFTENING
  When this is enabled, the gravitational softening lengths of hdyro cells are varied according to 
  their radius. To this end, the radius of a cell is multiplied by the parameter ``GasSoftFactor``. 
  Then, the closest softening from a logarithmicaly spaced table of possible softenings is adopted 
  for the cell. The minimum softening in the table is specified by the parameter 
  ``MinimumComovingHydroSoftening``, and the larger ones are spaced a factor 
  ``AdaptiveHydroSofteningSpacing`` apart. The resulting minimum and maximum softening values are 
  reported in the stdout log file. *Generally recommended for hydro simulations.*

MULTIPLE_NODE_SOFTENING
  If the tree walk wants to use a 'softened node' (i.e. where the maximum gravitational softening 
  of some particles in the node is larger than the node distance and larger than the target 
  particle's softening), the node is opened by default (because there could be mass components 
  with a still smaller softening hidden in the node). This can cause a subtantial performance 
  penalty in some cases. By setting this option, this can be avoided. The code will then be 
  allowed to use softened nodes, but it does that by evaluating the node-particle interaction for 
  each mass component with different softening type separately (but by neglecting possible shifts 
  in their centers of masses). This also requires that each tree node computes and stores a vector 
  with these different masses. It is therefore desirable to not make the table of softening types 
  excessively large. This option can be combined with adaptive hydro softening. In this case, 
  particle type 0 needs to be mapped to softening type 0 in the parameterfile, and no other particle 
  type may be mapped to softening type 0 (the code will issue an error message if one doesn't obey 
  to this). *Generally recommended for hydro simulations.*

INDIVIDUAL_GRAVITY_SOFTENING=2+4
  The code can also be asked to set the softening types of some of the particle types automatically 
  based on particle mass. The particle types to which this is applied are set by this compile time 
  option through a bitmask encoding the types. The code by default assumes that the softening of 
  particle type 1 should be the reference. To this end, the code determines the average mass of 
  type 1 particles, and the types selected through this option then compute a desired softening 
  length by scaling the type-1 softening with the cube root of the mass ratio. Then, the softening 
  type that is closest to this desired softening is assigned to the particle (*choosing only from those 
  softening values explicitly input as a SofteningComovingTypeX parameter*). This option is primarily 
  useful for zoon simulations, where one may for example lump all boundary dark matter particles 
  together into type 2 or 3, but yet provide a set of softening types over which they are automatically 
  distributed according to their mass. If both ``ADAPTIVE_HYDRO_SOFTENING`` and 
  ``MULTIPLE_NODE_SOFTENING`` are set, the softening types considered for assignment exclude softening 
  type 0. Note: particles that accrete matter (BHs or sinks) get their softening updated if needed.

NSOFTTYPES=6
  This can be changed to modify the number of available softening types. These must be explicitly 
  input as SofteningComovingTypeX parameters, and so the value of ``NSOFTTYPES`` must match the number 
  of these entries in the parameter file.

NSOFTTYPES_HYDRO=64
  This is only relevant if adaptive hydro softening is enabled and can be set to override the default 
  value of 64 for the length of the logarithmically spaced softening table. The sum of ``NSOFTTYPES`` 
  and ``NSOFTTYPES_HYDRO`` may not exceed 254 (this is checked). *Default value generally recommended.*


Particle-mesh (PM) gravity
^^^^^^^^^^^^^^^^^^^^^^^^^^

PMGRID=512
  Dimension of particle-mesh grid covering the domain.
  This enables the TreePM method, i.e. the long-range force is computed with a PM-algorithm, 
  and the short range force with the tree. The parameter has to be set to the size of the mesh that
  should be used, e.g. 64, 128, 256, etc. The mesh dimensions need not necessarily be a power of two, 
  but the FFT is fastest for such a choice.  Note: If the simulation is not in a periodic box, then 
  a FFT method for vacuum boundaries is employed, using a mesh with dimension twice that specified 
  by ``PMGRID``. *Recommended value:* for periodic boxes, twice the cube root of the number of dark 
  matter particles, i.e. 512 for a :math:`256^3` simulation.

ASMTH=1.25
  This factor expressed the adopted force split scale in the TreePM approach in units of the grid 
  cell size. Setting this value overrides the default value of 1.25, in mesh-cells, which defines the 
  long-range/short-range force split.

RCUT=6.0
  This determines the maximum radius, in units of the force split scale, out to which the tree 
  calculation in TreePM mode considers tree nodes. If a tree node is more distant, the corresponding 
  branch is discarded. The default value is 4.5, given in mesh-cells.

PLACEHIGHRESREGION=2
  If this option is set (will only work together with ``PMGRID``), then the long range force is 
  computed in two stages: One Fourier-grid is used to cover the whole simulation volume, allowing the 
  computation of the large-scale force.  A second Fourier mesh is placed on the region occupied by 
  "high-resolution" particles, allowing the computation of an intermediate-scale force. Finally, the
  force on very small scales is computed by the tree. This procedure can be useful for 
  "zoom-simulations", where the majority of particles (the high-res particles) are occupying only a 
  small fraction of the volume. To activate this option, the parameter needs to be set to an integer 
  that encodes the particle types that make up the high-res particles in the form of a bit mask. For 
  example, if types 0, 1, and 4 are the high-res particles, then the parameter should be set to
  ``PLACEHIGHRESREGION=1+2+16``, i.e. to the sum :math:`2^0+2^1+2^4`. The spatial region covered by 
  the high-res grid is determined automatically from the initial conditions. The region is recalculated 
  if one of the selected particles is falling outside of the high-resolution region. Note: If a periodic 
  box is used, the high-res zone is not allowed to intersect the box boundaries.

ENLARGEREGION=1.1
  This is only relevant when ``PLACEHIGHRESREGION`` is activated. The size of the high resolution 
  box will be automatically determined as the minimum size required to contain the selected 
  particle type(s), in a "shrink-wrap" fashion. This region is be expanded on the fly, as needed. 
  However, in order to prevent a situation where this size needs to be enlarged frquently, such as 
  when the particle set is (slowly) expanding, the minimum size is multiplied by the factor 
  ``ENLARGEREGION`` (if defined). Then even if the set is expanding, this will only rarely 
  trigger a recalculation of the high resolution mesh geometry, which is in general also associated 
  with a change of the force split scale. *Recommended value:* between 1.1 and 1.4. 

GRIDBOOST=1
  Normally, if ``PLACEHIGHRESREGION`` is enabled, the code will try to offer an effective grid size 
  for the high-resolution patch that is equivalent to ``PMGRID``. Because zero-padding has to be used 
  for the high-res inset, this gives a total mesh twice as large, which corresponds to ``GRIDBOOST=2``. 
  This value can here be increased by hand, to e.g. 4 or 8, to increase the resolution of the high-res 
  PM grid. The total mesh size used for the high-resolution FFTs is given by ``GRIDBOOST*PMGRID``.

ONLY_PM
  This option disables the short-range gravitational force calculation by the tree, i.e. the 
  gravitational force only consists of the PM-force. *Not recommended.*

PM_ZOOM_OPTIMIZED
  This option enables a different communication algorithm in the PM calculations (the one of Gadget3) 
  which works well independent of the data layout, in particular it can cope well with highly clustered 
  particle distributions that occupy only a small subset of the total simulated volume. However, 
  this method is a bit slower than the default approach (used when the option is disabled), which 
  is best matched for homogenously sampled periodic boxes. *Recommended only for large zoom runs.*

FFT_COLUMN_BASED
  When this is enabled, the FFT calculations are not parallelized in terms of a slab-decomposition 
  but rather through a column based approach. This scales to larger number of MPI ranks but is slower 
  in absolute terms as twice as many transpose operations neeed to be performed. It is hence only 
  worthwhile to use this option for very large number of MPI ranks that exceed the 1D mesh dimension.
  *Recommended only for very large runs (MPI task number exceeding 4096 or 8192?).*

CHUNKING
  This will calculate the gravity forces in interleaved blocks and is a pure optimization option 
  that can help to improve imbalances in the communication patterns in case multiple iterations 
  in a force calculation are needed due to an insufficient buffer size. *Always recommended.*



Voronoi Mesh
------------

Mesh motion & regularization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

VORONOI_STATIC_MESH
  Assumes the mesh to be static, i.e. to not change with time. The vertex velocities of all mesh-generating
  points is set to zero and domain decomposition is disabled.

VORONOI_STATIC_MESH_DO_DOMAIN_DECOMPOSITION
  Enables domain decomposition together with ``VORONOI_STATIC_MESH`` (which is otherwise then disabled), 
  in case non-gas particle types exist and the use of domain decompotions is desired. Note that on one 
  hand it may be advantageous in case the non-gas particles mix well or cluster strongly, but on the 
  other hand the mesh construction that follows the domain decomposition is slow for a static mesh, so 
  whether or not using this new flag is overall advantageous depends on the problem.


Mesh refinement & derefinement
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

REFINEMENT_HIGH_RES_GAS
  Limits the dynamical (de-)refinements of cells to cells which are either already present in the ICs 
  or are created with ``GENERATE_GAS_IN_ICS`` from type 1 particles. This adds an additional integer 
  quantity ``AllowRefinement`` to PartType0 in the snapshots indicating if a gas cell is allowed to be 
  refined and if it is, how often this cell has already been split: if 0, no splitting allowed. If odd 
  (starting at 1), the cell was already present in the ICs. If even (starting at 2), the cell was 
  generated from a type 1 particle. For values of 3 or more, ``floor((AllowRefinement-1)/2.0)`` gives the 
  number of times the cell was split. Note: the interaction of SofteningComovingOfPartType[3] with 
  ``REFINEMENT_HIGH_RES_GAS`` which occurs in ``convert_cell_into_star()``!

NODEREFINE_BACKGROUND_GRID
  The background grid will be prevented from derefining, when refinement is used. In practice, when 
  enabled this option will compute the mean cell volume during initialization (it also appears to 
  require this as a ``MeanVolume`` parameter which should be removed). Derefinement is then 
  disallowed during the run for all cells with ``Volume > 0.1 * MeanVolume``.

DEREFINE_GENTLY
  Prevent derefinement in regions of strong local changes of pressure or density. Cells will only be 
  refined if they temperature of surrounding cells is sufficiently similar (within a factor 
  GENTLE_DEREFINE_FACTOR, which has currently a default value of 1.2).

OPTIMIZE_MESH_MEMORY_FOR_REFINEMENT
  If activated some grid structures not needed for mesh refinement/derefinement are freed before the 
  function do_derefinements_and_refinements is called. The remaining mesh structures are freed after this 
  step as usual.



Hydrodynamics
-------------

Riemann solver
^^^^^^^^^^^^^^


Slope Limiting and Reconstruction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TVD_SLOPE_LIMITER
  Apply a componentwise total variation dimishing slope limiter.
  Use one of the ``TVD_SLOPE_LIMITER_*`` options to select the limiter.

TVD_SLOPE_LIMITER_ALBADA
  Chooses the van Albada limiter. Use together with option TVD_SLOPE_LIMITER.

TVD_SLOPE_LIMITER_MINBEE
  Chooses the minbee limiter. Use together with option TVD_SLOPE_LIMITER.

TVD_SLOPE_LIMITER_SUPERBEE
  Chooses the superbee limiter. Use together with option TVD_SLOPE_LIMITER.

TVD_SLOPE_LIMITER_VANLEER
  Chooses the van Leer limiter. Use together with option TVD_SLOPE_LIMITER.

TVD_SLOPE_LIMITER_MINMOD
  todo

TVD_SLOPE_LIMITER_MC
  todo

DISABLE_TIME_EXTRAPOLATION
  todo

DISABLE_SPATIAL_EXTRAPOLATION
  todo

NO_SCALAR_GRADIENTS
  todo

GRADIENTS_GREEN_GAUSS
  todo

GRADIENT_LIMITER_DUFFELL
  todo



Cooling & star formation
------------------------

MODIFIED_EOS
  The effective equation of state is modified close to the SF density threshold to make it harder.

SFR_KEEP_CELLS
  Forces derefinement of cells that formed stars?

SLOW_RELAX_TO_EOS
  Temperature of star-forming gas is relaxed *slowly* to the eEOS *only if* SLOW_RELAX_TO_EOS is defined.

STEEPER_SFR_FOR_STARBURST
  Changes the scaling of the star formation rate with density to :math:`\rm{SFR} \propto \rho^{\alpha}`
  for densities greater than the turnover point for the eEOS (where run-away sets in). Here, 
  :math:`\alpha` corresponds to the ``StarburstPowerLawIndex`` parameter which must then be set. 
  See :ref:`sfr-params` for documentation.

Output options
--------------

HAVE_HDF5
  If this is set, the code will be compiled with support for input and output in the HDF5 format. 
  You need to have the HDF5 libraries and headers installed on your computer for this option to work. 
  The HDF5 format can then be selected as format "3" in Arepo's parameterfile. *Always recommended.*

REDUCE_FLUSH
  If enabled files and stdout are only flushed after a certain time defined in the parameter file 
  (standard behaviour: everything is flashed most times something is written to it).

EVALPOTENTIAL
  When this option is set, the code will compute the gravitational potential energy each time a 
  global statistics is computed. This can be useful for testing global energy conservation.

TOLERATE_WRITE_ERROR
  TODO.
  Tries to store files for which things go wrong in an alternative output directory, as specified by 
  the ``AlternativeOutputDir`` parameter.

SUBBOX_SNAPSHOTS
  Instructs the code to write additional "subbox"-type snapshots in addition to normal snapshots. 
  Subbox snapshots are sub-volumes, and contain the full particle information within those sub-volumes. 
  They are meant to be written very frequently, for high time frequency analysis or visualization. 
  An arbitrary number of subboxes can be requested at once. Configuration is through a number of 
  required parameters - see :ref:`subbox-params` for documentation.


Optional additional outputs
^^^^^^^^^^^^^^^^^^^^^^^^^^^

OUTPUTPOTENTIAL
  This will force the code to compute gravitational potentials for all particles each time a 
  snapshot file is generated. These values are then included in the snapshot files. Note that the 
  computation of the values of the potential costs additional time.

OUTPUTACCELERATION
  This will include the physical acceleration of each particle in snapshot files.

OUTPUT_REFBHCOUNTER
  Write out refinement statistics of cells close to black holes. (todo: Config dependencies)

OUTPUT_CELL_SPIN
  Output finite-volume tracking of cell internal angular momentum/spin in order to check for global
  conservation of angular momentum.




Special Behaviour
-----------------

Generally options in this category are recommended only for specialized simulations.
They affect the physics, physical setup, or type of simulation.

ISOTHERM_EQS
  This special option makes the gas behave like an isothermal gas with equation of state 
  :math:`P = c_s^2 \rho`. The sound-speed :math:`c_s` is set by the thermal energy per unit mass 
  in the intial conditions, i.e. :math:`c_s^2=u`. If the value for :math:`u` is zero, then the 
  initial gas temperature in the parameter file is used to define the sound speed according to 
  :math:`c_s^2= k T/m_p`, where :math:`m_p` is the proton mass.

LONGIDS
  If this is set, the code assumes that particle-IDs are stored as 64-bit long integers. This is 
  only really needed if you want to go beyond ~2 billion particles.

RUNNING_SAFETY_FILE
  If used, run will not start in case a file './running' exists. This may be useful if the same job 
  is submitted several times (to deal with long queue wait times), to prevent the job from actually 
  starting several times. The file is automatically written at the beginning of the run (must be 
  manually deleted by the user after the run has ended).



Single- or double-precision
---------------------------

DOUBLEPRECISION=1
  If this option is disabled, the code uses single precision variables for storing most data 
  associated with particles and cells. If it is activated (or set to 1), then all particle/cell 
  data is stored in double precision internally. Alternatively, setting ``DOUBLEPRECISION=2`` 
  enables a mixed mode where some of the more crucial data is kept in double precision internally, 
  whereas the bulk is stored in single precision. Note that output files are nevertheless written 
  by converting the values that are saved to single precision. Finally, the new option 
  ``DOUBLEPRECISION=3`` uses a MySingle datatype, which for this setting becomes float, otherwise 
  it is treated like MyFloat.

DOUBLEPRECISION_FFTW
  This enables the use of double precision variables for storing the fields that are processed 
  with the Fourier transforms. If this is set, the code will also then use the double-precision 
  version of FTTW for the PM gravity calculations.

OUTPUT_IN_DOUBLEPRECISION
  When enabled, the code stores snapshot files (as well as group catalogues) in double precision. 
  This applies to all floating-point values, which then take 8 bytes each instead of 4.

INPUT_IN_DOUBLEPRECISION
  When this is enabled, initial conditions files or snapshot files used as code input are read as 
  double precision variables. (Normally, an error message should be produced if the input data is 
  instead stored in single precision.) *This should be deprecitated in the future for HDF5 ICs 
  from which the precision of all input datasets is determined automatically.*

AUTO_SWAP_ENDIAN_READIC
  Enables automatic swapping of the endianness when reading ICs in the file formats 1 and 2. 
  For file format 3, HDF5 is supporting this functionality as well. *This should be deprecitated 
  in the future for HDF5 ICs.*



Technical Code Behavior
-----------------------

These options are largely numerical, related to performance, optimizations, etc.

NOSTOP_WHEN_BELOW_MINTIMESTEP
  If this is activated, the code will not terminate when the timestep falls below the value of 
  ``MinSizeTimestep`` specified in the parameterfile. This is useful for runs where one wants to 
  enforce a constant timestep for all particles. This can be done by activating this option, and 
  by setting ``MinSizeTimestep`` and ``MaxSizeTimestep`` to an equal value.

NOTREERND
  If ``NOTREERND`` is not set, the tree construction for very close particle pairs may move a into 
  a randomized subnode among the 8 nodes, on a scale much smaller than the softening length. In this 
  case, we now enlarge the size of the enclosing tree nodes such that they are guaranteed to enclose 
  all the listed daughter particles.

MPI_HYPERCUBE_ALLGATHERV
  Implements a workaround for some MPI-libraries that use quite a lot of internal storage
  for MPI_Allgatherv.

OPTIMIZE_MEMORY_USAGE
  Optimize for memory, not for speed, related to the Voronoi mesh. Reduces sizeof(point) but requires 
  more bit shifts / operations for integer mapping. Note: this is dangerous for high dynamic range 
  simulations with mixed precision, since some position variables are singles instead of doubles.


NUM_THREADS
  todo

IMPOSE_PINNING
  todo

IMPOSE_PINNING_OVERRIDE_MODE
  todo

GENERIC_ASYNC
  todo


CUDA
  Activates support for CUDA. Currently no part of the code is optimized using cuda.

CUDA_INSTRUMENT
  Enables support for the cuda profiler. Sections marked with TIMER_START/STOP will be shown in 
  the timeline in nvvp.

VTUNE_INSTRUMENT
  Enables support for intel vtune. Sections marked with TIMER_START/STOP will be shown in the 
  timeline.



On the fly and post-processed analysis
--------------------------------------

Friends of friends
^^^^^^^^^^^^^^^^^^

FOF
  Master switch to enable the friends-of-friends group finder code. 
  This will then usually be applied automatically before snapshot files are written 
  (unless disabled selectively for certain output dumps).

FOF_PRIMARY_LINK_TYPES=2
  This option selects the particle types that are processed by the friends-of-friends linking 
  algorithm. A default linking length of 0.2 is assumed for this particle type unless specified 
  otherwise.

FOF_SECONDARY_LINK_TYPES=1+16+32
  With this option, FOF groups can be augmented by particles/cells of other particle types that 
  they "enclose". To this end, for each particle among the types selected by the bit mask specifed 
  with ``FOF_SECONDARY_LINK_TYPES``, the nearest among ``FOF_PRIMARY_LINK_TYPES`` is found and then 
  the particle is attached to whatever group this particle is in.

FOF_SECONDARY_LINK_TARGET_TYPES=1 
  *The motivation why this exists is unclear to me at present...*
  A new option to make the secondary linking work better in zoom runs (after the FOF groups have been 
  found, the tree is newly constructed for all the secondary link targets). This should normally be set 
  to all dark matter particle types. If not set, it defaults to ``FOF_PRIMARY_LINK_TYPES``, which 
  reproduces the old behaviour.

FOF_GROUP_MIN_LEN=32
  This sets the minimum total number of particles/cells (primary plus secondary) that a group must 
  have before it is stored in the group catalogue. The default value for this is 32.

FOF_LINKLENGTH=0.2
  This can be set to override the default value of the FOF linking length in units of the mean 
  particle spacing of the primary link type.

FOF_FUZZ_SORT_BY_NEAREST_GROUP=0
  DESCRIPTION NEEDED.

FOF_STOREIDS
  Normally, the snapshots produced with a FOF group catalogue are stored in group order, such that 
  the particle set making up a group can be inferred as a contiguous block of particles in the 
  snapsot file, making it redundant to separately store the IDs of the particles making up a group 
  in the group catalogue. By activating this option, one can nevertheless force to create the 
  corresponding lists of IDs as part of the group catalogue output.

ADD_GROUP_PROPERTIES
  This can be used to calculate additional properties for an already existing group catalogue. 
  The snapshot is loaded, the substructure calculations are made, and the new properties are then 
  added as additional datasets to the HDF5 group catalog files. The run then exits.

USE_AREPO_FOF_WITH_GADGET_FIX
  Recommened when applying Arepo's FOF group finder in post-processing mode on prexisting Gadget 
  snapshot files that contain gas. Then the SPH smoothing lengths stored in the Gadget file can be 
  used as initial guesses for the secondary linking stage.

Subfind
^^^^^^^

SUBFIND
  When enabled, this automatically runs the Subfind analysis of all FOF groups after they have been 
  found. This snapshot files are brought into subhalo order within each group.

SAVE_HSML_IN_SNAPSHOT
  When activated, this will store the hsml-values used for estimating total matter density around 
  every point and the corresonding densities in the snapshot files associated with a run of Subfind.

SUBFIND_MEASURE_H2MASS
  This is a special measuremenat option for the mass in molecular hydrogen found in each subhalo.

SUBFIND_CALC_MORE
  todo.
  Additional calculations are carried out in the Subfind algorithm, which may be expensive. 
  (i) The velocity dispersion in the local density estimate.
  (ii) The DM density around every particle is stored in the snapshot if this is set 
  together with ``SAVE_HSML_IN_SNAPSHOT``.

SUBFIND_EXTENDED_PROPERTIES
  todo.
  Additional calculations are carried out, which may be expensive.
  (i) Further quantities related to the angular momentum in different components.
  (ii) The kinetic, thermal and potential binding energies for SO halos.


Visualization
^^^^^^^^^^^^^

VORONOI_MESHOUTPUT
  todo.

VORONOI_IMAGES_FOREACHSNAPSHOT
  todo.

VORONOI_FREQUENT_IMAGES
  Can be used to create images with frequency ``TimeBetweenImages`` set in the parameter file, 
  independent of snapshot files. Note that this is only guaranteed to work correctly with a constant 
  and equal image spacing if TimeBetweenImages >= MaxSizeTimestep.

VORONOI_FIELD_DUMP_PIXELS_X=1024
  todo.

VORONOI_FIELD_DUMP_PIXELS_Y=1024
  todo.

VORONOI_VELOCITY_FIELD_2D
  todo.

VORONOI_FIELD_COMPENSATE_VX=0.0
  todo.

VORONOI_FIELD_COMPENSATE_VY=0.0
  todo.
  
VORONOI_NEW_IMAGE
  todo.

VORONOI_PROJ_TEMP
  todo.

VORONOI_PROJ
  todo.

VORONOI_MULTIPLE_PROJECTIONS
  todo.

VORONOI_NOGRADS
  todo.

DVR_RENDER
  Simplistic direct volumetric front-to-back alpha-composite RGBA orthogonal/perspective raycast 
  rendering using recursive F2B compositing equation (which allows for early ray termination and ray 
  profile mappings like first-hit, MIP).

DVR_RENDER_SMOOTH
  Interpolates scalar field from connected cells for smoother field.

DVR_RENDER_ORTHOGONAL
  todo.

DVR_NUM_FIELDS
  todo.

DVR_STAY_IN_BOX
  todo.
