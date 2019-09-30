
Getting Started
===============

Prerequisites
-------------

AREPO needs the following non-standard libraries for compilation:

* **MPI** - the Message Passing Interface (version 1.0 or higher). Many vendor supplied versions 
  exist, in addition to excellent open source implementations, e.g.
  `OpenMPI <http://www.open-mpi.org/>`_ or 
  `MPICH <http://www-unix.mcs.anl.gov/mpi/mpich/>`_.

* **GSL** - the `GNU scientific library <http://www.gnu.org/software/gsl>`_  is a required 
  open-source package. AREPO needs this library for a few simple cosmological integrations at 
  start-up, and for random number generation.

* **GMP** - the `GNU multiple precision arithmetic library <http://www.gnu.org/software/gmp>`_. 
  is a required open-source package. AREPO needs this library for big number integer arithmetic 
  in order to calculate exacy geometric predicates if needed.

* **FFTW** - the `Fastest Fourier Transform in the West <http://www.fftw.org>`_ is a required 
  open-source library, needed for simulations that use the TreePM algorithm. Note that the
  MPI-capable version **3.x** of FFTW is required, and that FFTW needs to be explicitly compiled 
  with parallel support enabled.  This can be achieved by passing the option ``--enable-mpi`` to 
  the configure script. You should also add ``--enable-type-prefix`` to obtain the libraries in 
  both a single and double precision version.

* **HDF5** - the `Hierarchical Data Format <http://hdf.ncsa.uiuc.edu/HDF5>`_. This library has 
  been developed by NCSA, is open-source, and optional. AREPO can be compiled without this library, 
  but then input/output in HDF5 format is not supported. The old binary file formats are no longer 
  updated, so HDF5 is highly recommended.

Note that if any of the above libraries is not installed in standard locations on your system, the 
``Makefile`` provided with the code may need slight adjustments. Similarly, compiler options,
particularly with respect to optimisations, may need adjustment to the
C-compiler that is used. 

The provided ``Makefile`` is compatible with GNU-make, so typing ``make`` will then build the 
executable assuming the default input configuration filename of ``Config.sh`` and the default 
output executable name ``Arepo``. You customize these with the ``CONFIG=`` and ``EXEC=`` options 
to make, for example::

  make CONFIG=Config_setup2.sh EXEC=Arepo_setup2


Configuration
-------------

There are two important ways in which the AREPO code is configured and controlled:

1. The ``Config.sh`` file contains a substantial number of compile-time options that need to be set 
appropriately for the type of simulation that is to be carried out. The code needs to be
recompiled whenever one of these options is changed.

* See :doc:`core_config_options`.

2. There is also a text parameterfile e.g. ``param.txt`` that is passed as an argument at run-time 
to the code. This file specifies a number of variables through simple keyword-value pairs.

* See :doc:`core_param_options`.


Compiling
---------

To build the code, do the following:

#. Copy the file ``Template-Config.sh`` to ``Config.sh``.
#. Edit ``Config.sh`` as needed for your application.
#. Specify the target system using the ``SYSTYPE`` option (see below).
#. Run ``make`` to produce the executable.

The ``Template-Config.sh`` contains a list of all available compile-time options, with most 
of them commented out by default. Before compiling for a run, this template should be copied to a
new file called ``Config.sh`` for customziation. To activate a certain feature, the corresponding 
parameter can then be commented in, and given the desired value, where appropriate. 

.. warning::

  Whenever one of the compile-time options described below is modified, a full 
  recompilation of the code may be necessary. To guarantee that this is done when a simple ``make`` is
  specified, all source files have been specified dependent on the Config.sh options. Alternatively, 
  one can also issue the command ``make clean``, which will erase all object files, followed by ``make``.

This technique has the disadvantage that different simulations may require different 
binaries of AREPO. If several simulations are run concurrently, there is hence the danger that a
simulation is started/resumed with the *wrong* binary. Note that while AREPO checks the plausibility 
of some of the most important code options, this is not done for all of them. To minimise the risk of 
using the wrong executable for a simulation, it is recommended to produce a separate executable for 
each simulation that is run. For example, a good strategy is to make a copy of the whole code together 
with its makefile in the output directory of each simulation run, and then to use this copy to compile 
the code and to run the simulation.

The ``SYSTYPE`` can be set in two ways, with an environment variable or with a file in the 
build directory. If present, the file has priority over the environment variable.

* To set with an environment variable e.g. ``export SYSTYPE=Magny`` either at the command 
  line or permanently in your ``.bashrc`` (or equivalent) file.

* Or, to set with a file:

  #. copy the ``Template-Makefile.systype`` file to ``Makefile.systype``
  #. uncomment your system in ``Makefile.systype``

.. note::

  The ``Config.sh`` file should not be checked in to the repository. During development, 
  new compile-time options should be added to the ``Template-Config.sh`` file only. Usually, they 
  should be added there in the disabled/default version.


Running
-------

In order to start the simulation code, a "parameterfile" needs to be specified. An additional 
optional numerical parameter can be used to signal whether a continuation from a set of restart 
files, or from a snapshot file, is desired. A typical command to start the code looks like
the following::

  mpiexec -np 8 ./Arepo parameterfile [restartflag]

This would start the code using 8 processors, assuming that the parallel environment uses the 
``mpiexec`` command to start MPI applications. Depending on the operating system, other commands 
may be required for this task, e.g. ``mpirun``. Note that the code can in principle be started 
using an arbitrary number of processors, but the communication algorithms will be most efficient for
powers of 2. It is also possible to use a single processor only, in which case the code behaves 
like a serial code.

The optional ``restartflag`` takes a numeric value, which invokes a partial mode of code operation. 
The primary start-up options are just 0, 1, or 2. 

* 0: (the default, if omitted) requests the code start from initial conditions.
* 1: signals a continuation from restart files.
* 2: can be used to restart from a snapshot file produced by the code.

Higher values for the start-up options will execute the special postprocessing options built into 
the code. These are currently:

* 3: Execution of FOF/SUBFIND on a snapshot dump.
* 4: Making of an image slice.
* 5: Making of a projected image.
* 6: Conversion of snapshots from one file format to another.
* 7: Calculate a velocity power spectrum for the gas cells.
* 8: Make a grid/orthographic projection (using ray-tracing through the Voronoi mesh).
* 9: Make a projection along an arbitrary axis.
* 10: Make a perspective camera projection.
* 11: Calculate power spectra of various quantities for TRACER_PARTICLEs.
* 12: Calculate two-point correlation function (for a given particle type).
* 13: Calculate power spectrum (for a given particle type).
* 14: Write out the Voronoi mesh.
* 15: Run the post-processing shock finder.
* 16: Write out a two-dimensional slice of the Voronoi mesh.
* 17: Write out snapshot dump with measured gradients.
* 18: Recalculate gravitational potential values.

More information about the additional parameters required for these options can be obtained by 
starting AREPO without any options.

For **examples of actual simulation setups** see the ``tests/`` directory.
