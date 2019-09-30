
SHOCK_FINDER
============

Shock finder implementation. The Mach number field provides valuable insights into the gas dynamics 
and can help in analysis of many simulations. The general method as well as interpretations of common 
shock patterns are described in the papers listed below.

This module generates several output fields. In the simplest configurations these are just two:

=================   ===============    ==================================    ==============================================================================================
HDF5 Dataset        Physical Factor    Units                                 Description
=================   ===============    ==================================    ==============================================================================================
Machnumber          :math:`1`          None                                  The Mach number of the gas cell, zero if no shock is present.
EnergyDissipation   :math:`1/a`        UnitMass/UnitLength*UnitVelocity^3    The dissipated energy (amount of kinetic energy irreversibly transformed into thermal energy).
=================   ===============    ==================================    ==============================================================================================

The units of EnergyDissipation correspond to Energy / Time. 
When running in post-processing mode on existing snapshots, additional information is also saved:

===================  =================  ====================================  ==============================================================================================
HDF5 Dataset         Physical Factor    Units                                 Description
===================  =================  ====================================  ==============================================================================================
Coordinates          :math:`a/h`        UnitLength                            Same as in snapshots.
Volume               :math:`a^3 / h^3`  UnitLength^3                          Same as in snapshots.
Velocities           :math:`\sqrt{a}`   UnitVelocity                          Same as in snapshots.
Density              :math:`h^2 / a^3`  UnitMass/UnitLength^3                 Same as in snapshots.
InternalEnergy       :math:`1`          UnitVelocity^2                        Same as in snapshots.
Temperature          :math:`1`          Kelvin                                Calculated in update_primitive_variables.c
Surface              :math:`a^2 / h^2`  UnitLength^2                          The surface area of the shock surface inside the cell.
Gen.Int.EnergyFlux   :math:`h^2 / a^3`  UnitMass/UnitLength^3*UnitVelocity^3  Generated internal energy flux due to shocks, multiply this by the surface area.
Shock direction      :math:`1`          none                                  Points towards the pre-shock region.
PreShockSoundSpeed   :math:`1`          UnitVelocity                          The sound speed in the pre-shock region.
PreShockDensity      :math:`h^2 / a^3`  UnitMass/UnitLength^3                 The pre-shock density.
PreShockTemperature  :math:`1`          Kelvin                                The pre-shock temperature.
===================  =================  ====================================  ==============================================================================================

Here, multiply by the given physical factor to remove all little h and scale factors, to arrive at physical units. 

Note: For (visual) projections of the Mach number field make sure to weight the Mach number with the 
energy dissipation, otherwise you underestimate the Mach number in your projections, since in most of 
the cells no shock is present and the Mach number as well as the energy dissipation is zero.


Usage
-----

The shock finder can be run in several modes, which are described here.

(1) On-the-fly, Before Snapshot Output
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you want to run the shock finder once before each snapshot dump and save the results into that snapshot.
Then, compile only with the additional flag ``SHOCK_FINDER_BEFORE_OUTPUT``. No additional parameters or flags 
are needed, and the shock finder will run with the default configuration.

Two additional fields ``Machnumber`` and ``EnergyDissipation`` will be generated for all gas cells (described above).

(2) On-the-fly, Continuous
^^^^^^^^^^^^^^^^^^^^^^^^^^

If you want to run the shock finder on-the-fly (during every local time step of the entire simulation).
Then, compile only with the additional flag ``SHOCK_FINDER_ON_THE_FLY``. No additional parameters or flags 
are needed, and the shock finder will run with the default configuration.

Two additional fields ``Machnumber`` and ``EnergyDissipation`` will be generated for all gas cells (described above).

(3) Post-processing Snapshots
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you want to run the post-processing shock finder and write standard format snapshot files which include 
the additional shock finder fields. Then, compile only with the additional flag ``SHOCK_FINDER_BEFORE_OUTPUT``. 
No additional parameters or flags are needed. The snapshot file base of the new file is "shocks", e.g. the 
new file is called "shocks_000.hdf5". Start with::

  mpirun -np <NumberOfProcessors> ./Arepo <ParameterFile> 15 <SnapNum> [optionally: <SubboxNum>]

Two additional fields ``Machnumber`` and ``EnergyDissipation`` will be generated for all gas cells (described above).

(4) Post-processing Snapshots, Detailed
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the post-processing operating mode, run the shock finder with the following flags in Config.sh::

  SHOCK_FINDER_POST_PROCESSING
  SHOCK_FINDER_AREPO
  UNLIMITED_GRADIENTS
  ZONE_JUMP_P
  ZONE_JUMP_T
  SHOCK_DIR_GRAD_T
  SHOCK_JUMP_T
  SURFACE_ANGLE_APPROX (for 2D: SURFACE_SPHERE_APPROX)
  RESET_WRONG_JUMPS
  RESET_WRONG_RHO_JUMPS

In this case, the usage of the post-processing shock finder is also with restartflag 15::

  mpirun -np <NumberOfProcessors> ./Arepo <ParameterFile> 15 <SnapNum> [optionally: <SubboxNum>]

Many additional properties of the detected shocks are saved, as described above.


Additional Parameters
---------------------

For running in post-processing modes the following parameters may be set:

* ``RayMemFac`` Controls the amount of memory for the rays.
  *Default value: 1.0*.

* ``MachMin`` The minimum Machnumber inside a cell for being in the shock zone.
  *Default value: 1.3*.

* ``RayStepsMax`` The maximum distance from the shock surface before the Mach number gets calculated 
  in number of cells.
  *Default value: 5*.

* ``OutputDirShockFinder`` The directory of the output of the shock finder.
  *Default value: output_shock_finder*.

* ``NumFilesPerOutput`` The number of output files written in parallel, the number of tasks has to be 
  a multiple of this number.
  *Default value: 1*.

* ``DistToBorder`` Only if ``SKIP_BORDER`` is used, there will be no shock zone closer to the border. 
  Make sure the skipped area consists of at least one cell layer.
  *Default value: varies*.

* ``SubboxCoordinatesPath`` Only if ``SUBBOX_SNAPSHOTS`` is used and a subbox snapshot is being processed.
  Then this is the only relevant subbox parameter for the sohck finder. Note: Include the coordinates of 
  ALL subboxes into this file.
  *Default value: subbox.txt*.


Additional Config.sh Options
----------------------------

For running in post-processing modes the following options are relevant/recommended and may be set:

FIX_SPH_PARTICLES_AT_IDENTICAL_COORDINATES
  Set this flag if there is a problem with the mesh construction. 
  This may happen if we read in a Snapshot file with single precision data.

NO_ID_UNIQUE_CHECK
  Set this flag if there is a problem in test_id_uniqueness. 
  This happens if a processor gets zero particles while reading in a snapshot.

SHOCK_JUMP_T_FLOOR=10000
  Set this flag (in Kelvin) in order to account for post-processing reionization. 
  This affects the fields ``PreShockTemperature``, ``PreShockSoundSpeed``, ``Machnumber`` and 
  ``Gen.Int.EnergyFlux``.

SKIP_BORDER
  Set this flag if the snapshot/subbox has reflective boundaries. 
  Rays are stopped if they have a certain distance to such a boundary.

SDUBBOX_SNAPSHOTS
  Set this flag if reading in a subbox snapshot.


Authors
-------

  * Kevin Schaal (kevin.schaal@h-its.org) (http://www.kmschaal.de/)


Usage Policy and Citation
-------------------------

This module is open and freely available for use in any project.

Papers to cite:

  * `Schaal et al. 2015, MNRAS, 446, 399 <http://adsabs.harvard.edu/abs/2015MNRAS.446.3992S>`_
  * `Schaal et al. 2016, MNRAS, 461, 4441 <http://adsabs.harvard.edu/abs/2016MNRAS.461.4441S>`_
