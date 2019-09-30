
Parameter File Options
======================

A number of features of AREPO are controlled with **run-time options**. The parameterfile for 
AREPO is a simple text file, consisting of pairs of tags and values. For each parameter, a separate 
line needs to be specified, first listing the name (tag) of the parameter, and then the assigned 
value, separated by whitespace. It is allowed to add further text behind the assigned parameter 
value. The order of the parameters is arbitrary, but each one needs to occur exactly one time, 
otherwise an error message will be produced. Empty lines, or lines beginning with a %-sign, are 
ignored and treated as comments.

Here parameters which are typically required by default (e.g. without any additional compile-time 
code options enabled) are described:

.. note::

  All options on this page are required to be present in the parameter file for all runs, 
  unless specified otherwise. Documentation for module-specific parameters can be found under 
  their respective pages.

.. contents::
  :local:

  

Relevant files and file formats
-------------------------------

InitCondFile
  The filename of the initial conditions file. Can be a relative or absolute path. The ICs 
  can be distributed in one or more files, as with snapshots. If using ICs with multiple files, 
  only the basename without the trailing ".n" should be specified here. Similarly, the ".hdf5" 
  extension can be omitted, and the code will automatically append it for ``ICFormat=3``.
  If a restart from a snapshot with the "2" option is desired, one needs to specify the snapshot 
  file here.

OutputDir
  Pathname of the output directory of the code. Can be a relative or absolute path. Must exist.

AlternativeOutputDir
  **(Optional)**. Only required if ``TOLERATE_WRITE_ERROR`` is set.
  TODO.

SnapshotFileBase
  Basename of snapshot files produced by the code.

OutputListFilename
  File with a list of the desired output times. Can be specified with a relative or absolute path. 
  See :doc:`core_output_list` for more information.

ICFormat
  The file format of the initial conditions. Currently, three different formats are supported, 
  selected by one of the choices "1", "2", or "3". Format "1" is the traditional fortran-style 
  unformatted format familiar from GADGET-1. Format "2" is a variant of this format, where each
  block of data is preceeded by a 4-character block-identifier. Finally, format "3" selects the 
  HDF-5 format (recommended).

SnapFormat
  Similar to ``ICFormat``, this parameter selects the file-format of snapshot dumps produced by 
  the code. *Recommended value: 3*.



CPU and memory limits
---------------------

MaxMemSize
  The memory allocate per MPI task, in megabytes. A contiguous memory arena of this total size is 
  allocated at startup, and then partitioned internally within AREPO for memory allocation and 
  deallocation requests. Can generally be set to ~95% of the total available, e.g. (memory per node 
  / number of MPI tasks per node), to leave room for OS tasks and MPI buffers. This value can be 
  changed on a restart to increase the amount of memory available to each task.

TimeLimitCPU
  CPU-time limit for the present submission of the code. If 85 percent of this time have been 
  reached at the end of a timestep, the code terminates itself and produces restart files. The 
  extra 15% is used to guarantee that there is enough time to safely finish the current time 
  step and write the restart files. This CPU time refers to the wall-lock time on a single 
  processor only.

CpuTimeBetRestartFile
  The value specfied here gives the time in seconds the code will run before it writes regularly 
  produced restart files. This can be useful to protect against unexpected interruptions (for 
  example due to a hardware problem) of a simulation, particularly if it is run for a long time. 
  It is then possible to resume a simulation from the last restart file, reducing the potential 
  loss to the elapsed CPU-time since this was produced.

ResubmitOn
  If set to "1", the code will try to resubmit itself to the queuing system when an interruption 
  of the run due to the CPU-time limit occurs. The resubmission itself is done by executing the 
  program/script given with ``ResubmitCommand``.

ResubmitCommand
  The name of a script file or program that is executed for automatic resubmission of the job to 
  the queuing system. Note that the file given here needs to be executable.



Basic code options
------------------

.. note::
  
  Some of the four parameters here are redundant in the sense that each requires a specific Config.sh 
  option to be enabled. For instance, if ``COOLING`` is set, then ``CoolingOn`` must equal 1, and vice 
  versa. This is to provide rough assurance that you are running a given parameter file with the 
  correct executable.

ComovingIntegrationOn
  If set to "0", the code assumes plain Newtonian physics, with time, positions, velocities, and 
  masses measured in the internal system of units.
  If set to "1", the code assumes that a cosmological integration in comoving coordinates should be  
  carried out, assuming an expanding universe described by the 'Cosmological parameters' below. In 
  a cosmological integration, the time variable is the scale factor.

PeriodicBoundariesOn
  If set to "1", periodic boundary conditions are assumed, with a cubical box-size of side-length
  ``BoxSize``. Particle coordinates are expected to be in the range ``[0,BoxSize[``. Can only be 
  set to zero if ``GRAVITY_NOT_PERIODIC`` is set. Note: refers to gravity only!

CoolingOn
  If set to "1", gas looses energy through a (optically-thin) radiative cooling model at each timestep.
  Can only be set to zero if ``COOLING`` is not set.

StarformationOn
  If set to "1", gas can (stochastically) convert into collisionless star particles based on a 
  star formation model.
  Can only be set to zero if ``USE_SFR`` is not set.



Cosmological parameters
-----------------------

Omega0
  Gives the total matter density, i.e. :math:`\Omega_m`, in units of the critical density at z=0 
  for cosmological simulations. Ignored if ``ComovingIntegrationOn=0``.

OmegaLambda
  Gives the vacuum energy density, i.e. :math:`\Omega_\Lambda` (cosmological constant) at z=0 for 
  cosmological simulations. Ignored if ``ComovingIntegrationOn=0``.

OmegaBaryon
  Gives the baryon density, i.e. :math:`\Omega_b`, in units of the critical densty at z=0 for 
  cosmological simulations. Ignored if ``ComovingIntegrationOn=0``.

HubbleParam
  This gives the Hubble constant ("little h") at z=0 in units of 100 km/sec/Mpc.  Note that this 
  parameter has been basically absorbed into the definition of the internal code units, such that 
  for gravitational dynamics and adiabatic gas dynamics the actual value assigned for ``HubbleParam`` 
  is not used by the code. Only used in cosmological simulations when conversions to physical 
  cgs units are required (e.g. for radiative cooling physics).
  Ignored if ``ComovingIntegrationOn=0``.

BoxSize
  The boxsize for the simulation, in units of the ``UnitLength_in_cm`` parameter. All particles 
  and gas cells in the ICs must have Coordinates within the range ``[0,BoxSize]`` in each dimension.



Output control
--------------

OutputListOn
  If set to "1", the code tries to read a list of desired output times from the file given in
  ``OutputListFilename`` (see :doc:`core_output_list` for documentation). Otherwise, output times 
  are generated equally spaced from the values assigned for ``TimeOfFirstSnapshot`` and 
  ``TimeBetSnapshot``.

TimeOfFirstSnapshot
  The time of the first desired snapshot file in case a file with output times is not specified. 
  For cosmological simulations, the value given here is the scale factor of the first desired output.

TimeBetSnapshot
  The time interval between two subsequent snapshot files in case a file with output times is not 
  specified. For cosmological simulations, this is a multiplicative factor applied to the time of 
  the last snapshot, such that the snapshots will have a constant logarithmic spacing in the scale 
  factor. Otherwise, the parameter is an additive constant that gives the linear spacing between 
  snapshot times.

TimeBetStatistics
  The code can be asked to measure the total kinetic, thermal, and potential energy in regular 
  intervals, and to write the results to the file given in ``EnergyFile``. The time interval between 
  two such measurements is given by this parameter, in an analogous way as with ``TimeBetSnapshot``. 
  Note that the compile time option ``EVALPOTENTIAL`` needs to be activated to obtain a measurement 
  of the gravitational potential energy.

NumFilesPerSnapshot
  The number of separate files requested for each snapshot dump. Each file of the snapshot will hold 
  the data of one or several processors, up to all of them. ``NumFilesPerSnapshot`` must hence lie 
  between 1 and the number of processors used. Distributing a snapshot onto several files can be done 
  in parallel and may lead to much better I/O performance, depending on the hardware configuration. 
  It can also help to avoid problems due to big files (>2GB) for large simulations. Note that initial
  conditions may also be distributed into several files, the number of which is automatically recognised 
  by the code and does not have to be equal to ``NumFilesPerSnapshot`` (it may also be larger than the 
  number of processors).

NumFilesWrittenInParallel
  The number of files the code may read or write simultaneously when writing or reading 
  snapshot/restart files. If the value of this parameter is larger than the number of 
  processors, it is capped by that.

FlushCpuTimeDiff
  **(Optional)**. Only if ``REDUCE_FLUSH`` is set.
  TODO.



Time integration
----------------

TimeBegin
  This sets the starting time of a simulation when the code is started from initial conditions. 
  For cosmological integrations, the value specified here is taken as the initial scale factor.

TimeMax
  This sets the final time for the simulation. The code normally tries to run until this time is 
  reached. For cosmological integrations, the value given here is the final scale factor.

TypeOfTimestepCriterion
  This parameter can in principle be used to select different kinds of timestep criteria for 
  gravitational dynamics. However, AREPO presently only supports the standard criterion "0".
  *Required value: 0*.

ErrTolIntAccuracy
  This dimensionless parameter controls the accuracy of the timestep criterion selected by
  ``TypeOfTimestepCriterion``. It is the variable :math:`\eta`, where the cosmological timestep 
  for collisionless particles scales as :math:`\Delta t \propto \eta^{1/2}`.
  *Typical value: 0.012*.

CourantFac
  This sets the value of the Courant parameter used in the determination of the hydrodynamical 
  timestep of gas cells. The hydrodynamical timestep is this value times the standard Courant 
  condition calculated for each cell.
  *Typical values: 0.1 - 0.4*.

MaxSizeTimestep
  This gives the maximum timestep a particle may take. This should be set to a sensible value in 
  order to protect against too large timesteps for particles with very small acceleration. For 
  cosmological simulations, the parameter given here is the maximum allowed step in the logarithm 
  of the expansion factor.
  For comoving runs, this is in units of :math:`\Delta \rm{ln} a`.
  Note that the definition of ``MaxSizeTimestep`` has **changed** compared to Arepo-1.1 for 
  cosmological simulations (what does this mean?).

MinSizeTimestep
  If a particle requests a timestep smaller than the value specified here, the code will normally 
  terminate with a warning message. If compiled with the ``NOSTOP_WHEN_BELOW_MINTIMESTEP`` option, 
  the code will instead force the timesteps to be at least as large as ``MinSizeTimestep``.


.. _sfr-params:

Star formation model
--------------------

.. note::

  These fundamental parameters of the `Springel & Hernquist (2003)`_ star formation model are 
  required if ``USE_SFR`` is set.

CritPhysDensity
  The critical physical density above which star formation may take place (in :math:`\rm{cm}^{-3}`).
  Used instead of ``CritOverDensity`` for non-comoving runs.

MaxSfrTimescale
  This is the variable :math:`t_0^\star`, the star-formation timescale at the threshold density, 
  such that the local star-formation timescale is then calculated as 
  :math:`t_\star(\rho) = t_0^\star (\rho / \rho_{\rm th})^{-1/2}`.
  Given in internal time units, or Gyr?
  *Illustris value: 2.27*.

CritOverDensity
  The critical (over-)density above which star formation may take place, where the threshold 
  density is then :math:`\rho_{\rm th} = \rm{CritOverDensity} \times ( 3 \Omega_b H^2 / (8 \pi G) )` 
  (redshift independent). Used in place of a critical physical density for comoving integrations.
  *Illustris value: 57.7* (which then comes out to the usual :math:`n_{\rm sfr} = 0.13 \rm{cm}^{-3}`).

TempSupernova
  The "supernova temperature" :math:`T_{\rm SN}` of the hot intercloud medium, in Kelvin.
  *Illustris value: 5.73e7*.

TempClouds
  The "cold cloud temperature" :math:`T_{\rm c}`, in Kelvin.
  *Illustris value: 1000.0*.

FactorSN
  The variable :math:`\beta`, giving the mass fraction of massive stars (:math:`> 8 M_\odot`) formed 
  for each initial population of stars. This is thus determined by the stellar IMF. Therefore, if 
  ``GFM_STELLAR_EVOLUTION`` is set this parameter is not allowed and is instead calculated internally 
  based on the e.g. Chabrier IMF.

FactorEVP
  The variable :math:`A`, giving the efficiency of the cloud evaporation process.
  *Illustris value: 573.0*.

TemperatureThresh
  Star formation is prevented for cells which are hotter than the eEOS and hotter than the 
  TemperatureThresh parameter (in Kelvin). If this parameter is very large (e.g. 1e20), then 
  nothing is changed compared to the base model. If this parameter is small (e.g. 0, 1e4, or 1e5) 
  then star-formation will be prevented in hot halo gas.
  *Illustris value: 0*.

FactorForSofterEQS
  **(Optional)**. Required only if ``SOFTEREQS`` is set.
  The interpolate weight factor :math:`q` (as in `Vogelsberger+ (2013)`_), where the temperature of gas 
  above the star-formation threshold is then set as 
  :math:`T = q \times T_{\rm eEOS} + (1-q) \times T_{\rm iso}`. Here 'eEOS' is the original model 
  and 'iso' is an isothermal constant given by the ``TempForSofterEQS`` parameter below.
  *Illustris value: 0.3*.

TempForSofterEQS
  **(Optional)**. Required only if ``SOFTEREQS`` is set.
  The variable :math:`T_{\rm iso}` described above for the (cold) isothermal temperature, in Kelvin, 
  with which to soften the effective equation of state model.
  *Illustris value: 1e4*.

StarburstPowerLawIndex
  **(Optional)**. Required only if ``STEEPER_SFR_FOR_STARBURST`` is set.
  The variable :math:`\alpha` which sets a new scaling of the star formation rate with density above 
  the turnover point for the eEOS, of :math:`\rm{SFR} \propto \rho^{\alpha}` (since 
  :math:`t_\star \propto \rho^{-\alpha}`). A value of 
  :math:`\alpha=0.5` therefore recovers the default behavior, while for example 
  :math:`\alpha=1.0` results in higher star formation rates at high densities.

.. _Vogelsberger+ (2013): http://adsabs.harvard.edu/abs/2013MNRAS.436.3031V
.. _Springel & Hernquist (2003): http://adsabs.harvard.edu/abs/2003MNRAS.339..289S


Treatment of empty space, temperature limits, cooling
-----------------------------------------------------

InitGasTemp
  This sets the initial gas temperature (assuming either a mean molecular weight corresponding to 
  full ionization or full neutrality, depending on whether the temperature is above or below 10^4 K) 
  in Kelvin when initial conditions are read. However, the gas temperature is only set to a certain
  temperature if ``InitGasTemp>0``, and if the temperature of the gas particles in the initial 
  conditions file is zero, otherwise the initial gas temperature is left at the value stored in the 
  IC file.

MinGasTemp
  A minimum temperature floor imposed by the code. This may be set to zero, but it may be desirable 
  to prevent the gas from becoming too cold, e.g. for resolution reasons or because of lower limits 
  in the implemented cooling function. (This value is converted by the code to a minimum thermal 
  energy per unit mass assuming the mean molecular weight of neutral gas).

MinimumDensityOnStartUp
  TODO.

LimitUBelowThisDensity
  TODO.

LimitUBelowCertainDensityToThisValue
  TODO.

MinEgySpec
  TODO.

TreecoolFile
  **(Optional)**. Only required if ``COOLING`` is set.
  TODO.

SelfShieldingFile
  **(Optional)**. Only required if ``UVB_SELF_SHIELDING`` is set.
  TODO.



Tree gravity
------------

TypeOfOpeningCriterion
  This selects the type of cell-opening criterion used in the tree walks. A value of "0" results 
  in standard Barnes & Hut, while "1" selects the relative opening criterion.
  *Required value: 1* (only implemented option in AREPO).

ErrTolTheta
  This gives the maximum opening angle :math:`\theta` if the BH criterion is used for the tree walk. 
  If the relative opening criterion is used instead, a first force estimate is computed using the BH 
  algorithm, which is then recomputed with the relative opening criterion.
  *Typical value: 0.7*.

ErrTolForceAcc
  The accuracy parameter for the relative opening criterion for the tree walk.
  *Typical value: 0.0025*.

MultipleDomains
  TODO.
  Pure numerical parameter related to computational efficiency, run output should be invariant to its value.
  *Typical values: ? (highly dependent on run type)*.

TopNodeFactor
  TODO.
  Pure numerical parameter related to computational efficiency, run output should be invariant to its value.
  *Typical values: ? (highly dependent on run type)*.

ActivePartFracForNewDomainDecomp
  TOOD.
  Pure numerical parameter related to computational efficiency, run output should be invariant to its value.
  *Typical values: 0.05 - 0.2 (highly dependent on run type)*.

MaxRMSDisplacementFac
  **(Optional)**. Only required with ``LEGACY_DISPLACEMENT_CONSTRAINT`` enabled. 
  This parameter is an additional timestep limitation criterion for the long-range integration in 
  case the TreePM algorithm is used. It limits the long-range timestep such that the rms-displacement 
  of particles per step is at most ``MaxRMSDisplacementFac`` times the mean particle separation, or 
  the mesh-scale of the short-range/long-range force split, whichever is smaller. 
  *Typical value: 0.125*.



Initial density estimate
------------------------

DesNumNgb
  This sets the desired number of nearest neighbors for an initial density/size estimate for gas 
  cells during code startup.
  *Typical value: 32 or 64*.

MaxNumNgbDeviation
  This sets the allowed variation of the number of neighbours around the target value ``DesNumNgb``.
  As with all deviation parameters, this is purely for computational efficiency.
  *Typical value: 0 or 1*.



System of units
---------------

.. unit_params_start

UnitLength_in_cm
  This sets the internal length unit in **cm/h**, where H_0 = 100 h km/sec/Mpc. 
  For example, the standard choice of ``3.085678e21`` sets the length unit to :math:`1.0 \, \rm{kpc/h}`.

UnitMass_in_g
  This sets the internal mass unit in **g/h**, where H_0 = 100 h km/sec/Mpc.  
  For example, the standard choice of ``1.989e43`` sets the mass unit to :math:`10^{10} M_\odot/\rm{h}`.

UnitVelocity_in_cm_per_s
  This sets the internal velocity unit in **cm/sec**. 
  For example, the standard choice of ``1e5`` sets the velocity unit to :math:`1.0 \, \rm{km/sec}`. 
  Note that the specification of ``UnitLength_in_cm``, ``UnitMass_in_g``, and
  ``UnitVelocity_in_cm_per_s`` also determines the internal unit of time.

GravityConstantInternal
  The numerical value of the gravitational constant G in internal units depends on the system of 
  units you choose. For example, for the choices above, G=43007.1 in internal units.  For 
  ``GravityConstantInternal=0``, the code calculates the value corresponding to the physical value
  of G automatically.  However, you might want to set G yourself.  For example, by specifying:
  ``GravityConstantInternal=1``, 
  ``UnitLength_in_cm=1``, 
  ``UnitMass_in_g=1``, and
  ``UnitVelocity_in_cm_per_s=1``, 
  one obtains a *natural* system of units. Note that the code will nevertheless try to use the 
  *correct* value of the Hubble constant in this case, so you should not set 
  ``GravityConstantInternal`` in cosmological integrations.

.. unit_params_end



Gravitational softening lengths
-------------------------------

SofteningComovingType%d
  A Plummer equivalent gravitational softening length, to be referenced by one or more specific 
  particle types. For cosmological simulations in comoving coordinates, this is interpreted as a 
  comoving softening length. Code length units.

SofteningMaxPhysType%d
  When comoving integration is used, this parameter gives the maximum physical gravitational 
  softening length corresponding to ``SofteningComovingTypeX`` (referenced by one or more specific 
  particle types depending on the entries of ``SofteningTypeOfPartTypeN``). Depening on the relative 
  settings of the *Comoving* and *MaxPhys* softenings, the code will hence switch from a softening 
  constant in comoving units to one constant in physical units.
  For example, if the *MaxPhys* value is exactly half the *Comoving* value, then particles using 
  this softening type will have comoving softening until :math:`z=1` and fixed physical softenings 
  after that point in time. Code length units.

SofteningTypeOfPartType%d
  For each particle type in the simulation which is involved gravitational calculations, it must be 
  assigned to a "softening type", a 0-based integer index corresponding to one of the above 
  ``SofteningComovingTypeX``/``SofteningMaxPhysTypeX`` entry pairs.

GasSoftFactor
  The gravitational softening length of a gas cell is this value times the cellsize, which is 
  calculated as the radius of the volume-equilvalent-sphere, i.e. 
  :math:`r_{\rm cell} = (3 V_{\rm cell} / 4\pi)^{1/3}`.
  *Typical value: 2.5*.

MinimumComovingHydroSoftening
  **(Optional)**. Only required if ``ADAPTIVE_HYDRO_SOFTENING`` is set 
  (see :ref:`gravity-details` for more information). 
  If this treatment for gas softenings is based used, a discrete spectrum of possible softening 
  lengths :math:`\epsilon_i` for gas cells is created at startup. It contains 
  ``i=0, ..., NSOFTTYPES_HYDRO-1`` entries, controlled by this 'minimum' parameter and the following 
  'spacing' parameter, such that :math:`\epsilon_i = \rm{min} \times \rm{spacing}^i`. Code length units.
  *Typical values: a reasonable fraction of SofteningComovingType0, e.g. 0.1 - 1.0 times this value*.

AdaptiveHydroSofteningSpacing
  **(Optional)**. Only required if ``ADAPTIVE_HYDRO_SOFTENING`` is set.
  The logarithmic spacing for the adaptive gas softenings table, as described above. 
  Must be larger than one.
  *Typical values: 1.01 - 1.4*.



Mesh regularization & refinement
--------------------------------

.. note::

  The following options for mesh regularization are generally required if ``VORONOI`` is set, unless 
  ``VORONOI_STATIC_MESH`` is also enabled.

CellMaxAngleFactor
  **(Optional, generally set)**. Required only if ``REGULARIZE_MESH_FACE_ANGLE`` is enabled, which is the 
  recommended choice.
  TODO.

CellShapingFactor
  **(Optional, generally not set)**. Required only if ``REGULARIZE_MESH_FACE_ANGLE`` is not enabled.
  TODO.

CellShapingSpeed
  TODO.

.. note::

  The following options for mesh refinement are required if ``REFINEMENT`` is enabled (generally true).

ReferenceGasPartMass
  For comoving runs, it can either be given a non-zero value, in which case the other *MassFactor 
  parameters take its value as a basis, or it can be given 0, in which case the code calculates the 
  mean cell mass, to which the *MassFactor parameters will refer. For non-comoving runs, it must be 
  given non-zero, otherwise the run will exit. *Recommended value: 0 (for periodic cosmo boxes)*.
  If ``REFINEMENT_HIGH_RES_GAS`` is enabled, then: if ReferenceGasPartMass==0 in the parameter file, then 
  all gas present in the ICs will be allowed to be (de-)refined (and the code calculates the reference
  mass as the mean mass of those cells for which (de-)refinement is allowed), and if that is not desired,
  then ReferenceGasPartMass should be set to the correct value, in which case only gas with initial 
  mass<1.2*ReferenceGasPartMass will be allowed to be (de-)refined.

TargetGasMassFactor
  The target gas cell mass, where (de-)refinement is triggered if a given cell deviates by more than a 
  factor of 2.0 above or below this value. Multiplicative factor with respect to the mean cell mass. 
  *Recommended value: 1.0*.

RefinementCriterion
  TODO.

DerefinementCriterion
  TODO.



Subfind
-------

DesLinkNgb
  **(Optional)**. Only required if ``SUBFIND`` is enabled at compile-time.
  The (integer) minimum number of particles/cells, of all types, for Subfind groups. If a Subfind group 
  is identified with fewer than this number of total particles/cells, it is discarded. Note that 
  this means many small friends-of-friends groups (with a nomimal minimum number of 32 member particles) 
  may frequently have no sufficiently large Subfind groups, and so will have ``GroupFirstSub==0`` 
  indicating that that FoF has no central subhalo in addition to no satellite subhalos.
  *Typical value: 20*.

ErrTolThetaSubfind
  **(Optional)**. Only required if ``SUBFIND`` is enabled at compile-time.
  This has the same meaning as the ``ErrTolTheta`` parameter, i.e. the tree opening angle used to 
  control the accuracy of the gravity calculation, for uses within the Subfind algorithm.
  *Typical value: 0.7*.



.. _subbox-params:

Subboxes
--------

.. note::

  Options for subbox output are required only if ``SUBBOX_SNAPSHOTS`` is set.

SubboxCoordinatesPath
  Absolute path (including filename) to a text file containing position information for requested 
  subboxes. Each line of this file adds one subbox, and should include six floating point numbers 
  separated by spaces. These correspond to ``X_min X_max Y_min Y_max Z_min Z_max``, i.e. the 
  minimum and maximum coordinates along each axes, in code length units, which should be included 
  in that subbox. Periodic boundaries are not handled for subboxes, so these cannot cross a box 
  boundary.

SubboxMinTime
  The simulation time (scale factor for cosmological runs) at which subbox outputs should begin.

SubboxMaxTime
  The simulation time (scale factor for cosmological runs) at which subbox outputs should stop.

SubboxSyncModulo
  Subbox snapshots are typically written once for each global (largest) timestep of the run. The 
  total number of such timesteps depends sensitively on the simulation setup and particularly 
  resolution (for example, it was ~4000 for Illustris-1, ~2260 for Illustris-2, and ~1430 for 
  Illustris-3). This number may be undesirably large, in which case setting this parameter to an 
  integer value greater :math:`M>=1` means that only every :math:`M^{\rm th}` global timestep will 
  output a subbox snapshot. For example, the total number will be reduced by half if this parameter 
  is set to 2, and reduced to a quarter if it is set to 4.
  *Typical value: 1*.

SubboxNumFilesPerSnapshot
  In analogy to the ``NumFilesPerSnapshot`` parameter for full snapshots, the number of file chunks 
  to use for subbox snapshots (may be different).

SubbxNumFilesWrittenInParallel
  In analogy to the ``NumFilesWrittenInParallel`` parameter for full snapshot (typo here sic), the 
  number of MPI tasks which should simultaneously write to disk.



Visualization
-------------

TimeBetweenImages
  **(Optional)**. Only required if ``VORONOI_FREQUENT_IMAGES`` is enabled at compile-time.
  Sets the desired time between on-the-fly image generation.
