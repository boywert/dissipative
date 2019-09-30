
BLACK_HOLES
===========

The module for black holes includes methods for BH seeding/merging, movement/positioning, 
accretion, and feedback.

The majority of these aspects were originally developed for the "Illustris model".
Subsequent additions (particularly to the feedback model) were developed for the 
"Auriga model" as well as (particularly to the low-state feedback model) for the 
"IllustrisTNG model".

.. _Sijacki+ 2007: http://adsabs.harvard.edu/abs/2007MNRAS.380..877S
.. _Vogelsberger+ 2013: http://adsabs.harvard.edu/abs/2013MNRAS.436.3031V
.. _Weinberger+ 2016: http://tbd


Usage
-----

TODO


Additional Parameters
---------------------

* ``BlackHoleAccretionFactor`` todo.
* ``BlackHoleFeedbackFactor`` todo.
* ``BlackHoleEddingtonFactor`` todo.
* ``SeedBlackHoleMass`` todo.
* ``MinFoFMassForNewSeed`` todo.
* ``DesNumNgbBlackHole`` todo.
* ``BlackHoleMaxAccretionRadius`` todo.
* ``BlackHoleRadiativeEfficiency`` todo.

If the ADIOS low-state model is enabled:

* ``QuasarThreshold`` todo.
* ``RadioFeedbackFactor`` todo.
* ``RadioFeedbackReiorientationFactor`` todo.
* ``RadioFeedbackMinDensityFactor`` todo.

If the bubble low-state (quasar-mode) is enabled:

* todo

If the new centering is enabled:

* ``BlackHoleCenteringMassMultiplier`` todo.

If the AGN radiation is enabled:

* ``TreecoolFileAGN`` todo.
* ``SelfShieldingDensity`` todo.
* ``ObscurationFactor`` todo.
* ``ObscurationSlope`` todo.


Additional Config.sh Options
----------------------------

BLACK_HOLES
  todo


Seeding and Accretion
^^^^^^^^^^^^^^^^^^^^^

DRAINGAS=1
  todo

BH_EXACT_INTEGRATION
  todo

BH_BONDI_DEFAULT
  todo

BH_BONDI_DENSITY
  todo

BH_BONDI_CAPTURE
  todo

BH_BONDI_DISK_VORTICITY
  todo

BH_DO_NOT_PREVENT_MERGERS
  todo

BH_USE_GASVEL_IN_BONDI
  todo

BH_USE_ALFVEN_SPEED_IN_BONDI
  Take into account the alfven speed of the plasma in the bh accretion rate calculation (requires ``MHD``).

MASSIVE_SEEDS
  Once BH seeds are generated the total mass of a predefined number of gas neighbors is transfered to 
  the BH seed dynamical mass. Ideally with this option no repositioning is needed anymore, while the 
  mass conservation is ensured.

MASSIVE_SEEDS_MERGER
  todo


Positioning
^^^^^^^^^^^

BH_NEW_CENTERING
  todo

BH_FRICTION
  todo

BH_FRICTION_AGGRESSIVE
  todo

BH_HARMONIC_OSCILLATOR_FORCE
  todo

BH_DRAG
  todo

REPOSITION_ON_POTMIN
  todo


Other
^^^^^

BH_PRESSURE_CRITERION
  todo

BH_RELATIVE_NGB_DEVIATION
  todo

OUTPUT_BLACK_HOLE_TIMESTEP
  todo

REFINEMENT_AROUND_BH_FIXED
  todo

SUPPRESS_SF_IN_REFINEMENT_REGION
  todo

BH_BIPOLAR_FEEDBACK
  todo


Feedback
^^^^^^^^

.. warning::

  Status of certain aspects of the bubble-based models with current code versions needed.

BH_THERMALFEEDBACK
  todo

BH_THERMALFEEDBACK_ACC
  todo

BH_NF_RADIO
  todo

UNIFIED_FEEDBACK
  todo (obsolete?)

BH_BUBBLES
  todo (obsolete?)

BH_MAGNETIC_BUBBLES
  todo

BH_MAGNETIC_DIPOLAR_BUBBLES
  todo

BH_ADIOS_WIND
  kinetic black hole feedback in the low accretion state. Without BH_ADIOS_WIND_DIRECTIONAL or 
  BH_ADIOS_RANDOMIZED, this will kick the gas cells within a given radius around the black hole
  radially outward in a momentum conserving way. Do NOT use this together with BH_BUBBLES, 
  BH_MAGNETIC_BUBBLES, BH_MAGNETIC_DIPOLAR_BUBBLES, UNIFIED_FEEDBACK, BH_NF_RADIO.
  Problems: The model is known to produce large, low density regions around black holes, in particular
  in cosmological simulations with low resolution.

BH_ADIOS_DENS_DEP_EFFICIANCY
  todo

BH_ADIOS_WIND_WITH_QUASARTHRESHOLD
  todo

BH_ADIOS_WIND_WITH_VARIABLE_QUASARTHRESHOLD
  todo

BH_ADIOS_WIND_DIRECTIONAL
  requires BH_ADIOS_WIND. This introduces a preferred direction for the wind, randomly chosen at every
  injection event, in which the wind is more prominent. Do NOT use together with BH_ADIOS_RANDOMIZED.

BH_ADIOS_RANDOMIZED
  requires BH_ADIOS_WIND. This changes the ADIOS_WIND model and the kinetic kicks are in one direction.
  The injection only happens once the available energy exceeds some fraction of the dark matter velocity
  dispersion times enclosed gas mass.
  This preserves momentum only in a time-average sense, but behaves nicer in the injection regions, 
  without any (known) numerical effects.
  Problems: This model can create quite fast velocity kicks, which lead to very small timesteps and a
  large timebin-hierarchy.
  Model-paper: Weinberger et al. (2016)

BH_ADIOS_ONLY_ABOVE_MINIMUM_DENSITY
  todo

BH_CONTINOUS_MODE_SWITCH
  todo


Authors
-------

* Volker Springel
* Mark Vogelsberger
* Debora Sijacki
* Rainer Weinberger
* and others


Usage Policy and Citation
-------------------------

Please contact the authors before using this code for a new project. 
Co-authorship on papers may be requested.

To use aspects of BLACK_HOLES for a new simulation project, including un-modified versions 
of the "Illustris", "Auriga", or "IllustrisTNG" models un-modified, please 
**contact together (at a minimum): Mark, Volker, Lars, and Annalisa.**

Papers to cite:

* Depends on application.
