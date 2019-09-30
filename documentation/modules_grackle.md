
GRACKLE
=======

.. warning::

  This module is not yet fully tested, do not use!!


This module replaces the standard cooling implementation with functions from the Grackle library
(see https://grackle.readthedocs.org). Grackle provides:

* Two options for primordial chemistry:

  #. non-equilibrium primordial chemistry network. This is the default behaviour when GRACKLE is
     uncommented in Config.sh. The network can be switched on in three stages. By default atomic
     species are enabled (HI, HII, HeI, HeII, HeIII and e). GRACKLE_H2 enables molecular hydrogen
     cooling and chemistry (HM, H2I, H2II). GRACKLE_D enables deuterium (DI, DII, HDI). The fractional
     abundances of these species are tracked and advected.

  #. tabulated H and He cooling rates calculated from Cloudy, activated with GRACKLE_TAB

* tabulated metal cooling rates calculated from Cloudy

* optional photo-heating and photo-ionization from a UV-background


Usage
-----

Grackle must be installed separately (see https://grackle.readthedocs.org) and then linked to Arepo as an 
external library. This will require modifying the Arepo Makefile to correctly find and link Grackle. Set
GRACKLE_INCL and GRACKLE_LIBS appropriately in the SYSTYPE section for your machine e.g.::

  GRACKLE_INCL = -I/path/to/grackle/installation/include
  GRACKLE_LIBS = -L/path/to/grackle/installation/lib -lgrackle

making sure that the installation location is on your LD_LIBRARY_PATH.

When the code is run (with non-equilibrium chemistry), primordial abundances of the relevant species are
assumed and the correct ionization states based on internal energy provided in the ICs are obtained by 
iterating on cells until converged. This behaviour can be overidden with ``GRACKLE_ABUNDANCE_IN_ICS``,
whereupon Arepo will look for the species abundances in the ICs.

When enabled, Grackle completely replaces the standard cooling implementation; this means that any calls
to functions found in cooling.c will fail, with the exception of cooling_only() which has a new definition
in grackle.c. Calling this function evolves the chemistry of the cells, carries out cooling and updates the
pressure of the cells. To reiterate, because GRACKLE replaces the standard cooling implementation it should 
be assumed that any extra modules beyond 'vanilla' Arepo that rely in some way on the cooling routines 
(e.g. star formation) will not be compatible with GRACKLE until they are appropriately modified; however, 
this should not be particularly difficult.

When Grackle is enabled, cell temperatures and fractional abundances of relevant species are outputed in 
snapshots.


Additional Parameters
---------------------

Certain parameters must be provided in the parameterfile:

* ``GrackleOn`` Master switch, 1 for on, 0 for off
* ``GrackleRadiativeCooling`` Use radiative cooling
* ``GrackleMetalCooling`` Use metal cooling
* ``GrackleUVB`` UV background, 0 for off, 1 for Faucher-Giguere et al. (2009), 2 for Haardt & Madau (2012)
* ``GrackleDataFile`` path to tabulated Cloudy datafile
* ``GrackleInitialMetallicity`` If metals are not being tracked in the simulation, 
  this value gives adefault metallicity to be used by Grackle


Additional Config.sh Options
----------------------------

GRACKLE_H2
  see above

GRACKLE_D
  see above


Authors
-------

* Matthew C. Smith (m.c.smith@ast.cam.ac.uk)


Usage Policy and Citation
-------------------------

If Grackle is used, the authors of the library request the following recognition:

If you used the Grackle library in your work, please cite it as “the Grackle chemistry and cooling library 
(The Enzo Collaboration et al. 2014; Kim, J. et al. 2014).” Also, please add a footnote to 
https://grackle.readthedocs.org/.

* The Enzo Collaboration, Bryan, G. L., Norman, M. L., et al. 2014, ApJS, 211, 19
* Kim, J.-h., Abel, T., Agertz, O., et al. 2014, ApJS, 210, 14
