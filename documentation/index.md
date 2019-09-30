
AREPO
=====

AREPO is a massively parallel code for hydrodynamical cosmological simulations. It is a 
flexible code that can be applied to a variety of different types of simulations, offering 
a number of sophisticated simulation algorithms. An account of the numerical algorithms 
employed by the code is given in the original code paper, subsequent publications, and this 
documentation.

AREPO was written by Volker Springel (volker.springel@h-its.org) with further development 
by many authors. Many parts of the code are in active development.


.. toctree::
   :maxdepth: 2
   :caption: General Usage

   getting_started
   core_config_options
   core_param_options

   core_unit_system
   core_output_list


Beyond the simplest test problems one or more specialized code options will typically be used. 
These are broadly termed 'modules' and are activated through compile-time flags. Each may have 
additional compile-time options and/or parameter file values.

.. toctree::
   :caption: Special Behavior
   :glob:
   :maxdepth: 1

   modules
   modules_*


.. toctree::
   :maxdepth: 2
   :caption: Development

   code_changelog
   code_style_guide
   code_comm_structure


To contribute to this documentation, update either in-code comments or supplementary files in 
the ``documentation/`` directory, then commit. A new version will be built and placed here 
automatically after a few minutes. :ref:`genindex`.
