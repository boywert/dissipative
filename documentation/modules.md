
Code Modules List
=================

The AREPO code contains a number of physics modules that take the
form of extensions of the code for specific science applications. 
Examples include MHD, non-equilibrium chemistry, tracer particles,
etc.

Usually, these modules have been developed by a single person or a
small group of people, and are proprietary to them. If other people
want to use these code parts for their own project, they have to
check with the corresponding authors first. Depending on the policy
they have defined, such a use may (or may not) be possible and may 
involve a request for co-authorship in publications. 

In any case, please check this list for a first info about the
usage policies of different, clearly defined modules. In case of
doubt, please feel free to discuss with the corresponding module
authors.


-----


Solvers and Equation Systems
----------------------------

AMR
^^^
Use cartesian-based adaptive mesh refinement instead of the usual Voronoi approach for the 
discretization of the hydrodynamics.

* Documentation: :doc:`modules_amr`


MHD
^^^
Ideal MHD implementation. This implementation uses cell-centered magnetic fields
and employs the Powell/Dedner scheme for divergence cleaning.

* Documentation: :doc:`modules_mhd`


MHD_CT
^^^^^^
Constrained transport Ideal MHD implementation. More accurate and robust (and experimental) 
than the Powell/Dedner scheme.

* Documentation: :doc:`modules_mhd_ct`


GLOBAL_VISCOSITY
^^^^^^^^^^^^^^^^
Solve Navier-Stokes instead of Euler equations.

* Documentation: TODO.


SPECIAL_RELATIVITY
^^^^^^^^^^^^^^^^^^
Special relativistic hydro solver for Arepo.

* Documentation: :doc:`modules_special_relativity`


DG
^^
Discontinuous galerkin method for hydrodynamics instead of usual linear gradient approach.

* Documentation: TODO.


RT
^^
Radiative transfer module (cone method).

* Documentation: TODO.


OTVET
^^^^^
Radiative transfer module (optically thin variable eddington tensor method).

* Documentation: TODO.


FLD
^^^
Radiative transfer module (flux limited diffusion method).

* Documentation: TODO.


CONDUCTION
^^^^^^^^^^
Thermal conduction (multiple methods?).

* Documentation: TODO.


-----


Physical Models (Galaxy)
------------------------

GFM
^^^
The 'galaxy formation model' with several components originally developed for the Illustris 
simulation. Includes galactic-scale winds, stellar evolution and enrichment.

* Documentation: :doc:`modules_gfm`


GFM_CHEMTAGS
^^^^^^^^^^^^
Chemical tagging to follow in more detail the origin of mass/metals from different processes.

* Documentation: :doc:`modules_gfm_chemtags`


GFM_DUST
^^^^^^^^
Production and tracking of dust molecules.

* Documentation: :doc:`modules_gfm_dust`


BLACK_HOLES
^^^^^^^^^^^
Treatment of supermassive blackholes, including their seeding, accretion, and feedback. 

* Documentation: :doc:`modules_black_holes`


GRACKLE
^^^^^^^
Replaces the standard cooling implementation with functions from the Grackle library.

* Documentation: :doc:`modules_grackle`


COSMIC_RAYS
^^^^^^^^^^^
Cosmic Ray implementation. Work in progress.

* Documentation: :doc:`modules_cosmic_rays`


FM_SFR
^^^^^^
Stellar feedback and ISM physics module of FM, intended for 'resolved' ISM resolution runs.

* Documentation: TODO.


SIDM
^^^^
Self-interacting dark matter model.

* Documentation: TODO.


-----


Physical Models (Other)
-----------------------

SGCHEM
^^^^^^
A set of routines for modelling gas-phase chemistry and cooling.
Currently implements a very simple H, C, O network.

* Documentation: :doc:`modules_sgchem`


CIRCUMSTELLAR
^^^^^^^^^^^^^
Models for circumstellar disks.

* Documentation: TODO.


EOS_OPAL
^^^^^^^^
OPAL equation of state for stellar applications.

* Documentation: :doc:`modules_eos_opal`


-----


Analysis
--------

TRACER_*
^^^^^^^^
Implementation of tracer particles that allow to follow the flow of gas (as well as stars/BHs 
if relevant) in a 'Lagrangian' way. Also the passive scalar field.

* Documentation: :doc:`modules_tracer`


SHOCK_FINDER 
^^^^^^^^^^^^
Hydrodynamic shock finder implementation. 

* Documentation: :doc:`modules_shock_finder`

