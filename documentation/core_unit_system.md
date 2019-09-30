
Unit System
===========

These parameters control the units used for a simulation:

.. include:: core_param_options.md
  :start-after: .. unit_params_start
  :end-before: .. unit_params_end


Unit metadata in output
-----------------------

Some information about the physical unit can be stored in the hdf5 output.
The unit is stored as a bunch of hdf5 attributes on each dataset which has
unit information available.


The length_scaling, mass_scaling, and velocity_scaling attributes describe
how the quantity has to change if the length, mass and velocity units are
changed, i.e. if the quantity has to be multiplied with 0.01^n when you go
from cm to meter, length_scaling would be set to n. a_scaling and
h_scaling are just the exponents of the scaling with the cosmological a
and h factors. The attribute to_cgs gives the factor to convert the output
to cgs units (excluding a and h factors).


.. note::

  Almost every unit can be described with these 6 numbers. A bit more tricky
  is the temperature unit in Kelvin. Either, the unit system needs to be
  extend with another scaling factor, or one adopts a system of units with
  the Boltzmann constant :math:`k_B` equal to 1 and set the to_cgs factor for Kelvin
  quantities to the real value of :math:`k_B`.
  Also to be addressed: magnitudes.


To enable unit output, the function init_unit must be called after
init_field in ``io_fields.c``. Examples are::


  init_units(IO_SFR, 0.0, 0.0, -1., 1.0, 1.0, SOLAR_MASS/SEC_PER_YEAR);
  /* Msun/yr */

  init_units(IO_GFM_AGN_RADIATION, 0.0, 0.0, -3.0, 1.0, 3.0, 1.0);
  /* bolometric intensity in physical cgs units of erg/s/cm^2 */


Some python code that can handle the output for unit conversions and guessing 
the usual name can be found in inspector gadget (unit.py). This part of the code 
does not depend on any other inspector gadget related code and should be fairly 
easy to portable to other i/o libraries. The unit class represent a unit and 
handles conversion/guessing a common name for a set of unit parameters. The 
Quantity class is derived from numpy.ndarray and is a numpy array with unit 
information. This class tries to carry on unit information through most 
mathematical operations. The idea is to not alter the outcome of any operation 
and just augment the result with some unit information. I.e. if you add up 1 
meter with 1 centimetre you end up with 2 and a unit of none instead of 1.1 m. 
The file can be found at:

* https://bitbucket.org/abauer/inspector_gadget/src/cc69deae3a5ffffa8eee916d40ebe5f7b949bc82/gadget/units.py?at=master&fileviewer=file-view-default

