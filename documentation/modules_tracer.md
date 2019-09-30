
TRACER_*
========

Implementation of tracer particles that allow one to follow the flow of baryonic mass between 
gas cells (as well as stars/BHs if relevant) in a quasi-Lagrangian sense. Otherwise, because 
gas cells exchange mass each timestep in the finite volume scheme of AREPO, it is fundamentally 
impossible to determine the Lagrangian history (i.e. thermal or dynamical history) of any gas 
parcel.

There are three different tracer particle modules implemented, and one passive scalar field: 

1. ``TRACER_MC`` are "Monte Carlo tracers" than belong to parents (and so exist only as children 
to those cells/particles), without phase space coordinates (no positions or velocities). They are 
exchanged between parents in proportion to mass fluxes.

2. ``TRACER_PARTICLE`` are "Velocity field tracers" which have positions and velocities and are 
advected through space according to the local (linear gradient interpolated) gas velocity field.

3. ``TRACER_TRAJECTORY`` is a customized/extended version of velocity field tracers, intended 
primarily for recording tracer properties through time for e.g. nucleosynthesis postprocessing. 
It is not generalized.

4. ``TRACER_FIELD`` enables a passive scalar field exchanged in direct proportion to mass flux.
Useful for idealized problems.

-----


TRACER_MC
---------

Usage
^^^^^

Tracer "particles" which are probabilistically transfered between gas cells and star/BH particles 
based on the mass flux through each face (gas<->gas), the mass ratio of the spawned star (gas->star), 
the return of mass due to stellar mass loss (star->gas), and the accreted mass for sink type 
particles (gas->BH).

.. note::

  For on-the-fly group finding, although MC tracers are included in `GroupLenType` and 
  `SubhaloLenType`, they are not output in group order. However this is not generally a drawback, 
  since following tracers across snapshots anyways requires a global snapshot search in most cases.

To run a simulation with MC tracers use the ``TRACER_MC=N`` configuration flag, where ``N`` is an 
empty particle type that the tracers will output as in snapshots. 

Currently, the only way to introduce MC tracers is uniformly in all gas cells during the read of the
initial condition file. To do this, set ``TracerMCPerCell`` in the parameter file to the starting 
number of tracers per cell, and set the ``GENERATE_TRACER_MC_IN_ICS`` config flag. The note for
velocity field tracers about SPLIT_PARTICLE_TYPE applies to Monte Carlo tracers as well.

Additional Parameters
^^^^^^^^^^^^^^^^^^^^^

* ``TracerMCPerCell`` the number of Monte Carlo tracers to generate for each parent gas cell after 
  reading the initial conditions file. Only applies if ``GENERATE_TRACER_MC_IN_ICS`` is enabled. 
  Recommended value: from 1 to 10 (integer).


Additional Config.sh Options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TRACER_MC=3
  Master switch for MC tracers. Value specifies particle type to use for output.

GENERATE_TRACER_MC_IN_ICS
  If enabled, one or more tracers are created within each gas cell during the read of the initial 
  conditions (RestartFlag == 0 only). The number of tracers desired is set by the ``TracerMCPerCell`` 
  parameter. Note that using this option when gas cells have different initial masses, as in the 
  case of multi-mass zoom simulations, violates the typical implicit assumption that all Monte Carlo 
  tracers represent fluid elements of equal mass. For the typical (old) scheme, this requires only 
  additional thought in post-processing. For the directional (new) scheme, this violates a key 
  assumption of the scheme and will modify the evolution of tracers, with more thought required.

TRACER_MC_NUM_FLUID_QUANTITIES
  The number of fluid quantities to store for each tracer. Must match the number of fields added 
  to the ``TRACER_MC_STORE_WHAT`` bitmask below.

TRACER_MC_STORE_WHAT
  A bitmask specifying which fluid quantities to store for each tracer. For example, ``1+8`` stores 
  the maximum temperature and maximum density encountered. These values are generally recorded from 
  the parent gas cell for each tracer, although some (e.g. "wind counter") are calculated based on 
  more complex metrics. Maximum-type fields (e.g. "Tmax" or "Tmax Time") are reset to zero immediately 
  after writing a snapshot *unless* the ``TRACER_NO_RESET_EACH_SNAP`` option is set. See ``allvars.h`` 
  for bitmask options. See `public Illustris documentation 
  <http://www.illustris-project.org/data/docs/specifications/#tracerquant>`_ for some further field 
  descriptions.

TRACER_NO_RESET_EACH_SNAP
  For maximum-type fluid quantity fields (e.g. "Tmax" or "Tmax Time"), disable zero'ing immediately 
  after writing a snapshot. Such quantities are then global over the entire simulation time up until 
  that snapshot has been written. This makes post-processing easier, at the cost of loss of information.

TRACER_MC_CHECKS
  Enable frequent (and somewhat expensive) consistency checks on the internal tracer linked list 
  structures, terminating if any problems are encountered. Also increases memory usage.

Authors
^^^^^^^

* Shy Genel (shygenelastro@gmail.com)
* Dylan Nelson (dnelson@mpa-garching.mpg.de)
* Mark Vogelsberger (mvogelsb@mit.edu)


Usage Policy and Citation
^^^^^^^^^^^^^^^^^^^^^^^^^

There is no requirement to contact the authors before using these modules,
or for co-authorship on publications. We are happy, however, to provide advice or 
actively collaborate with interested users.
The only exception to the above is questions broadly related to
"cosmological gas accretion from the IGM into galaxies". On science related to this,
please contact the developers prior to use of the tracer particles.
  
Papers to cite:

* `Genel et al. 2013 MNRAS, 435, 1426 <http://adsabs.harvard.edu/abs/2013MNRAS.435.1426G>`_
* Genel et al. 2016 (in prep)


-----


TRACER_PARTICLE
---------------

Usage
^^^^^

Tracer particles which are passively advected through the domain by interpolating the velocity 
field of the underlying gas. These particles live in P and in timebins, but have no effect on 
the simulation.

To run a simulation with "velocity field" tracers use the ``TRACER_PARTICLE=N`` configuration flag, 
where N is an empty particle type (1-NTYPES) that the tracers will use and output as in snapshots.

Tracers can be introduced directly into the initial conditions by creating them with the target 
particle type. In this case, only position is important, velocity is set by the code, and mass 
should not be included. If the ICs have non-uniform gas masses and the goal is to sample the gas 
density field with the tracers, their distribution must be constructed manually in the ICs. 
Likewise to increase the ratio of tracers to gas cells above one.

Alternatively, the ``GENERATE_TRACER_PARTICLE_IN_ICS`` config flag can be used to place one tracer 
exactly at the position of each gas point during the read of the initial conditions. Useful to match 
the density field of any complicated gas geometry in the ICs when the masses are uniform.

.. note::

    ``GENERATE_TRACER_PARTICLE_IN_ICS`` generates the tracer particles directly from gas particles, 
    regardless of whether those have been read from the ICs or generated using ``GENERATE_GAS_IN_ICS``.

.. note::

    FoF can be run either on the fly or post-processed with velocity field tracers, which are 
    included explicitly in FoF groups (nearest DM particle attachment), but not in subfind groups.


Additional Parameters
^^^^^^^^^^^^^^^^^^^^^

* ``MinimumTracerHsml`` is a minimum distance (in code units) used when a tracer searches for its 
  nearest cell. This should be a rough estimate of the minimum cell size.


Additional Config.sh Options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TRACER_PARTICLE=2
  Master switch for Velocity Field tracer particles. Value specifies particle type to use.

GENERATE_TRACER_PARTICLE_IN_ICS
  If enabled, a single tracer particle is placed at the position of each gas cell position during 
  the read of the initial conditions (RestartFlag == 0 only).

TRACER_PART_NUM_FLUID_QUANTITIES
  The number of fluid quantities to store for each tracer. Must match the number of fields added 
  to the ``TRACER_PART_STORE_WHAT`` bitmask below.

TRACER_PART_STORE_WHAT
  A bitmask specifying which fluid quantities to store for each tracer. For example, ``1+8`` stores 
  the maximum temperature and maximum density encountered. These values are recorded from the parent 
  (closest) gas cell for each tracer. They are reset to zero immediately after writing a snapshot 
  *unless* the ``TRACER_NO_RESET_EACH_SNAP`` option is set. See ``allvars.h`` for bitmask options.


Authors
^^^^^^^

* Shy Genel (shygenelastro@gmail.com)
* Dylan Nelson (dnelson@mpa-garching.mpg.de)
* Mark Vogelsberger (mvogelsb@mit.edu)


Usage Policy and Citation
^^^^^^^^^^^^^^^^^^^^^^^^^

There is no requirement to contact the authors before using this module,
or for co-authorship on publications. We are happy, however, to provide advice or 
actively collaborate with interested users.
The only exception to the above is questions broadly related to
"cosmological gas accretion from the IGM into galaxies". On science related to this,
please contact the developers prior to use of the tracer particles.
  
Papers to cite:

* `Genel et al. 2013 MNRAS, 435, 1426 <http://adsabs.harvard.edu/abs/2013MNRAS.435.1426G>`_
* Genel et al. 2016 (in prep)


-----


TRACER_TRAJECTORY
-----------------

Usage
^^^^^

Enable the ``TRACER_TRAJECTORY`` option and use either the auto-generation or the init file to 
specify the initial tracer positions. Output destination file and frequencies must be specified 
by the two parameters, respectively. Fields to record and save can be modified by hand in 
``tracer_trajectory.c``.


Additional Parameters
^^^^^^^^^^^^^^^^^^^^^

* ``NumberOfTracersToGenerate`` is required only if ``TRACER_TRAJECTORY_GENERATE`` is set. In this 
  case, this number of new velocity tracers (per gas cell) are generated at startup and randomly 
  placed within each cell.

* ``TracerInitFile`` is required if ``TRACER_TRAJECTORY_GENERATE`` is not set. In this case, previously 
  saved tracer positions are read from this file. The format is simple binary: the total number of 
  tracers (4 byte int) as a header, followed by xyz coordinate triplets (8 byte doubles each).

* ``TracerOutputFile`` specifies a file name to write the (simple binary) output into. Seems not 
  entirely used ("tracer.dat" still hard-coded).

* ``TracerOutputConfFile`` specifies an external text file which controls the timing of tracer outputs. 
  The first line is an integer, the total number of tracer output steps. Subsequent lines, one per step, 
  each contain two spaced-separated floats, the first being an output time, the second being the 
  output time period. For each output time in the file, tracers are written until that time with the 
  cadence given by the output time period.


Additional Config.sh Options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TRACER_TRAJECTORY
  Master switch to enable saving of tracer "trajectories" (parent gas cell properties) in time.
  Requires that ``TRACER_PARTICLE`` be also enabled.

TRACER_TRAJECTORY_GENERATE
  If specified, velocity tracers are generated at startup within each gas cell. The number desired is 
  given by the ``NumberOfTracersToGenerate`` parameter. If not specified, tracers are assumed to have 
  been already generated in the past and saved, in which case they are read from the ``TracerInitFile``.

TRACER_TRAJECTORY_EXTENDED_OUTPUT
  Normally six values (x,y,z,rho,temp,u) are saved for each tracer, plus its unique ID. If this 
  option is enabled, four additional values (v_x,v_y,v_z,dedt) as well as ``EOS_NSPECIES`` abundances 
  are also saved.


Authors
^^^^^^^

* Ruediger Pakmor (ruediger.pakmor@h-its.org)


Usage Policy and Citation
^^^^^^^^^^^^^^^^^^^^^^^^^

TODO (Ruediger).

Papers to cite:

* todo.


-----


TRACER_FIELD
------------

Usage
^^^^^

Enable the ``TRACER_FIELD`` option and specify the value of a ``TracerField`` for each gas cell 
in the initial conditions. For example, to tag a specific region of space (e.g. the lower half 
of the box in the 2D coffee problem, or the inner slab of a KH instability problem) this could be 
set to unity in some cells and zero elsewhere. In subsequent snapshot outputs, the amount of 
``TracerField`` in a given gas cell gives an indication of the amount of mass from the tagged 
region which has reached (or remains in) this cell.


Additional Parameters
^^^^^^^^^^^^^^^^^^^^^

* None.


Additional Config.sh Options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TRACER_FIELD
  Master switch to enable tracking of a passive scalar field. If enabled, a ``TracerField`` field 
  should be included for each gas cell in the initial conditions. The "amount" of this field, 
  which is then fluxed between cells in direct proportion to mass, is then output in snapshot 
  files for each gas cell.


Authors
^^^^^^^

* Volker Springel (volker.springel@h-its.org)


Usage Policy and Citation
^^^^^^^^^^^^^^^^^^^^^^^^^

Part of original codebase. There is no requirement to contact the author before using 
this module, or for co-authorship on publications.

Papers to cite:

  * `Springel et al. 2010, MNRAS, 401, 791 <http://adsabs.harvard.edu/abs/2010MNRAS.401..791S>`_


------


Reference
---------

**Note: This is an experimental section which pulls doxygen info from the source code. Testing.**

.. doxygenfile:: tracer_mc.c

.. doxygenfile:: tracer_particle.c

.. doxygenfile:: tracer_trajectory.c

----------

Testing:

.. .. doxygenstruct:: tessellation
.. .. doxygenfunction:: consider_moving_tracers
