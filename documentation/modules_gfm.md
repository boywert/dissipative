
GFM
===

The 'galaxy formation module' includes aspects for gas (metal) cooling and heating, stellar 
evolution and enrichment, and galactic-scale (decoupled) winds from star formation feedback.

The majority of these aspects were originally developed for the "Illustris model", and are 
described in `Vogelsberger+ 2013`_ with observational comparisons in `Torrey+ 2014`_. 
Subsequent additions, particularly to the winds, were developed for the "IllustrisTNG model", 
and are not yet described in published work.

.. _Vogelsberger+ 2013: http://adsabs.harvard.edu/abs/2013MNRAS.436.3031V
.. _Torrey+ 2014: http://adsabs.harvard.edu/abs/2014MNRAS.438.1985T


Usage
-----

The ``GFM`` master switch enables the module, and acts on top of ``USE_SFR`` and/or ``COOLING``. 
As the number of configuration options is large, not all possible permutations and combinations 
have been tested. However, it is theoretically possible to use only specific subsets of the 
funcionality described herein.

In particular, ``GFM_STELLAR_EVOLUTION`` is the master switch for 
stellar evolution with mass and metal return to nearby gas cells, which can optionally be enabled 
in a passive only mode (no mass loss). Similarly, ``GFM_WINDS`` is the master switch for the 
kinetic wind scheme, and ``GFM_METAL_COOLING`` enables the use of pre-tabulted Cloudy tables 
for metal line cooling. These three options have many optional sub-options, and many required 
parameters.

An exact set of configuration options and corresponding parameters can be used to enable the 
"Illustris model" or "IllustrisTNG model" (which also depend on e.g. BH-related options). See 
the `Usage Policy and Citation`_ section.


Additional Parameters
---------------------

For output, if ``GFM_STELLAR_PHOTOMETRICS`` is enabled, then:

* ``PhotometricsTablePath`` is a filesystem path to a *directory* which contains a file named 
  ``stellar_photometrics.hdf5`` which contains stellar population luminosities in several bands as a 
  function of the population age and initial metallicity. See below for the exact format of this 
  file and script used to produce it for Illustris/TNG runs (BC03 models).

Enrichment:
^^^^^^^^^^^

SNIa power-law delay-time-distribution (DTD) is the default, which implies the meanings for the 
parameters ``SNIa_Rate_Norm`` and ``SNIa_Rate_TAU`` as given below.

* ``IMF_MinMass_Msun`` is the lower mass limit for the normalization of the IMF (e.g. Chabrier 2003).
  In units of solar masses.
  *Illustris/TNG value: 0.1*.

* ``IMF_MaxMass_Msun`` is the upper mass limit for the normalization of the IMF (e.g. Chabrier 2003).
  In units of solar masses.
  *Illustris/TNG value: 100*.

* ``AGB_MassTransferOn`` if set to "1", mass and metals from the mass loss and chemical enrichment of stars 
  in the AGB phase are returned to surrounding gas cells, representing the local ISM. Note here that 
  "AGB winds" only include mass transfer (no energy or momentum transfer to the gas). If set to "0" this 
  physical process is disabled.

* ``SNIa_MassTransferOn`` if set to "1", mass and metals from the mass loss and chemical enrichment of 
  supernovae Type-Ia (SNIa) are returned to surrounding gas cells. The Type Ia rates are based on a power-law 
  delay time distribution (DTD). If set to "0" this process is disabled.

* ``SNII_MassTransferOn`` if set to "1", mass and metals from the mass loss and chemical enrichment of 
  supernovae Type-II (SNII, core-collapse) are returned to surrounding gas cells. If set to "0" this process 
  is disabled.

* ``SNII_MinMass_Msun`` is the lower mass limit for stars ending in SNII. Note that the number of SNII events 
  are calculated at every time step by integrating the chosen IMF between the minimum and maximum SNII mass 
  limits.
  *Illustris value: 6.0, IllustrisTNG value: 8.0*.

* ``SNII_MaxMass_Msun`` is the upper mass limit for stars ending in SNII.
  *Illustris/TNG value: 100*.

* ``SNIa_Rate_TAU`` is the minimum stellar age for SNIa going off, in Gyr. This is the offset time between 
  the birth of a single-age stellar population and the first expected SNIa event. Note: 40 Myr corresponds 
  to the main-sequence life time of 8 solar mass stars, which is the adopted upper mass limit for 
  progenitors of SNIa events.
  *Illustris/TNG value: 0.04*.

* ``SNIa_Rate_Norm`` is the integrated number of SNIa per solar mass formed, i.e. the normalization of the 
  DTD for SNIa, a power-law in time for time greater than ``SNIa_Rate_TAU``.
  *Illustris/TNG value: 0.0013 (adopted from Maoz+ 2012)*.

* ``YieldTablePath`` is a filesystem path to a *directory* which contains several yield tables used by the 
  enrichment routines. By changing the files in this directory, the physical yield models can be updated. 
  Four files are required: ``AGB.hdf5``, ``SNIa.hdf5``, ``SNII.hdf5``, and ``Lifetimes.hdf5``. 
  Note that **there are multiple generations of the yield files (be careful which you have)**. See below 
  for details.

* ``DesNumNgbEnrichment`` the destination number of neighbors for stellar enrichment (integer value), assuming
  that the typical neighbor mass is ``ReferenceGasPartMass``. This parameter is contains
  **explicit resolution dependence**, e.g. if increasing one "resolution step" corresponding to a factor of 
  8 in particle mass, to enrich the same gas mass (and therefore spatial volume), this value should also 
  be increased by a factor of 8. However, this leads to very large neighbor searches for high resolution 
  simulations, given that a *reasonable minimum must be enforced* to avoid strange enrichment patterns.
  *Typical values: 20 to 256*.

* ``MaxNumNgbDeviationEnrichment`` maximum deviation (e.g. integer plus/minus allowed) in the enrichment 
  neighbor search. As with all neighbor search deviation parameters, it is primarily meant to increase 
  computational efficiency.
  *Typical values: 0, 1, or between sqrt(DesNumNgbEnrichment) and DesNumNgbEnrichment/8*.

* ``IMFslope`` is required only if ``GFM_CONST_IMF==1`` is set. In this case, it gives the slope 
  of a pure power law IMF.


Initial Gas Metallicity:
^^^^^^^^^^^^^^^^^^^^^^^^

* ``PreEnrichTime`` is required only if ``GFM_PREENRICH`` is enabled.
  The simulation time (scale factor for comoving runs) at which point to set the total metallicity 
  and individual species metals for all gas cells. At this point, any existing values are overwritten. 
  The desired values are specified by the following parameter.

* ``PreEnrichAbundanceFile`` is required only if ``GFM_PREENRICH`` is enabled.
  The full path (including file name) to a text file specifying the initial abundances to use for 
  gas pre-enrichment. The format of the file is: one element per line. For each line, the element 
  name should be followed by its pre-enrichment mass fraction. Example contents of such a file::

    Hydrogen         0.76
    Helium           0.24
    Carbon           1.0e-10
    Nitrogen         1.0e-10
    Oxygen           1.0e-10
    Neon             1.0e-10
    Magnesium        1.0e-10
    Silicon          1.0e-10
    Iron             1.0e-10
    OtherMetals      1.0e-10

* ``CoolingTablePath`` is required only if ``GFM_COOLING_METAL`` is enabled.
  The full path (including file name) to the pre-tabulated cooling table in HDF5 format.
  Note that **there are multiple versions of the cooling tables (be careful which you have)**. See below 
  for details.

* ``MinMetalTemp`` is required only if ``GFM_COOLING_METAL`` is enabled.
  The minimum gas temperature (in Kelvin) below which metal line cooling will be turned 
  off. For unrestricted metal line cooling set MinMetalTemp=0. To avoid low temperature fine-structure 
  cooling set MinMetalTemp=1e4K.


Winds:
^^^^^^

* ``WindEnergyIn1e51erg`` gives the energy of each SNII in units of 1e51 erg. This parameter is required if 
  **GFM_STELLAR_EVOLUTION & (GFM_WINDS | GFM_WINDS_LOCAL)**, otherwise the alternative parameter 
  ``WindEnergyFraction`` is present. They have the same physical meaning in principle, but in the former case 
  the wind is decoupled from the SF sub-grid model. This means that the wind energy is determined only by the 
  IMF (how many SNII) and the WindEnergyIn1e51erg parameter. To reproduce the Illustris settings (which had 
  WindEnergyFraction=3.0 and a specific choice of SF sub-grid parameters), one should now set 
  WindEnergyIn1e51erg=1.0946745 instead. This is the variable :math:`\rm{egy}_w` (`Vogelsberger+ 2013`_) 
  in units of :math:`\rm{egy}_w^0`.

* ``WindEnergyFraction`` see above.

* ``VariableWindVelFactor`` is required only if either of ``GFM_WINDS_LOCAL`` or ``GFM_WINDS_VARIABLE`` 
  is enabled. It is the variable :math:`\kappa_w` (`Vogelsberger+ 2013`_), the coefficient of the primary 
  wind velocity scaling. The wind velocity (i.e. initial wind particle launch speed) is then calculated as 
  :math:`v_w = \kappa_w \sigma`. If ``GFM_WINDS_VARIABLE=1`` (decoupled sigma winds) or ``GFM_WINDS_LOCAL``, 
  then :math:`\sigma = \sigma_{\rm DM}^{\rm 3D} / \sqrt{3} = \sigma_{\rm DM}^{\rm 1D}` the local 
  one-dimensional dark matter velocity dispersion, estimated by Subfind at the position of this gas cell. 
  If ``GFM_WINDS_VARAIBLE=0`` (decoupled halo-mass based sigma winds) then 
  :math:`\sigma = v_{\rm vir} / \sqrt{2}` where the virial velocity is computed from the FoF halo mass.
  *Illustris value: 3.7*.

* ``VariableWindSpecMomentum`` is the variable :math:`\rm{mom}_w` (`Vogelsberger+ 2013`_) in units of km/s. 
  If this parameter is zero, winds are "energy-driven", otherwise they can be "momentum-driven".
  *Illustris value: 0*.

* ``WindFreeTravelMaxTimeFactor`` is the maximum free travel time (in units of the Hubble time at the 
  current simulation redshift) for wind particles under the "hydrodynamic decoupling" scheme. Only applies 
  if they haven't yet reached their minimum density.
  *Illustris value: 0.025*.

* ``WindFreeTravelDensFac`` is a density below which wind particles will recouple under the 
  "hydrodynamic decoupling" scheme. This is in units of the star formation threshold density. Only applies 
  if they haven't yet reached their maximum travel time.
  *Illustris value: 0.05*.

* ``TimeBetOnTheFlyFoF`` is a multiplicative factor setting the relative spacing of the on-the-fly friends 
  of friends halo finder. This is used for setting e.g. wind parameters in halos, and so has a physical 
  impact on the simulation. The FoF algorithm is run at a scale factor equal to the current scale factor 
  times this parameter value.
  *Illustris value: 1.03*.

* ``MinWindVel`` is a minimum wind velocity to enforce globally, for winds being kicked in halos of all 
  masses. It is in physical code units (km/s).
  *Illustris value: 0*.

For setting the thermal content and metallicity of winds:

* ``ThermalWindFraction`` is required only if ``GFM_WINDS_THERMAL_NEWDEF`` is enabled. In this case, 
  it is the variable :math:`\alpha` such that the thermal energy of wind particles at launch is calculated 
  as :math:`U_{\rm therm} = \alpha / (1-\alpha) * v_w^2 / 2` where :math:`v_w` is the wind launch velocity 
  as described for ``VariableWindVelFactor``. Otherwise, :math:`U_{\rm therm}=0`.

* ``ThermalWindFactor`` is required only if ``GFM_WINDS_THERMAL`` is enabled. It is an alternative 
  definition of the above. In this case, it is the variable :math:`\alpha` such that the thermal energy 
  of wind particles at launch is calculated as :math:`U_{\rm therm} = \alpha * (3/2) * \sigma^2` where 
  :math:`\sigma` is as described for ``VariableWindVelFactor``. Otherwise, :math:`U_{\rm therm}=0`.

* ``WindDumpFactor`` is required only if ``GFM_WINDS_STRIPPING`` is enabled. This parameter is equal to
  :math:`(1.0 - \gamma_w)` (`Vogelsberger+ 2013`_) where the metallicity of a newly launched wind particle 
  is then set as :math:`Z_w = \gamma_w Z_{\rm ISM}`.
  *Illustris value: 0.6 (such that the wind metallicity is 40% of the original ISM metallicity)*.

For wind models with additional complexity (scalings), the following parameters may apply:

* ``WindEnergyReductionFactor`` is required only if ``GFM_WIND_ENERGY_METAL_DEPENDENCE`` is enabled.
  TODO.

* ``WindEnergyReductionMetallicity`` is required only if ``GFM_WIND_ENERGY_METAL_DEPENDENCE`` is enabled.
  TODO.

* ``WindEnergyReductionExponent`` is required only if ``GFM_WIND_ENERGY_METAL_DEPENDENCE`` is enabled.
  TODO.

* ``WindSuppressionRedshift`` is required only if ``GFM_WINDS_HUBBLESCALING`` is enabled.
  It sets the redshift :math:`z_s` below which the wind energy is suppressed. Then the wind energy, 
  before going into the mass-loading calculation, is multiplied by an efficiency factor :math:`\epsilon` 
  (capped at 1.0) which is calculated as :math:`\epsilon = h(a) / h(1 / (1+z_s))` where 
  :math:`h(x) = \Omega_0 / x^3 + (1 - \Omega_0 - \Omega_L) / x^2` is the usual hubble function.

* ``VariableWindMassScale`` is required only if ``GFM_WINDS_MASSSCALING`` is enabled.
  TODO.



Additional Config.sh Options
----------------------------

GFM
  Master switch required to enable stellar enrichment *and/or* galactic-scale wind routines, with any 
  of their associated sub-options. Enabled alone, it should have no physical impact on the run (is 
  this still true?), with some important numerical modifications. Namely, the creation of the 
  auxiliary ``StarP`` array of structures, which then stores additional fields for each star particle 
  which need not be contained inside the primary ``P``.

Stellar Evolution and Enrichment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

GFM_STELLAR_EVOLUTION=0
  todo

GFM_CONST_IMF=1
  Sets the desired model for the stellar initial mass function (IMF).
  ``GFM_CONST_IMF=0`` implements a Chabrier IMF (and is default, even if not used explicitly).
  ``GFM_CONST_IMF=1`` implements a pure power-law IMF, and requires a parameter ``IMFslope`` 
  (e.g., -2.35 will then give a Salpeter IMF).

GFM_VARIABLE_IMF=0
  As an alternative to a constant (in time and space) IMF, ``GFM_VARIABLE_IMF=0`` implements a 
  power-law IMF that depends on the DM velocity dispersion around the star-forming cell. In 
  particular, the slope is given by :math:`\alpha = -(2.3 * \log(\sigma/200) + 2.13)` following 
  `Spiniello+ 2014 <http://adsabs.harvard.edu/abs/2014MNRAS.438.1483S>`_ where 
  :math:`\sigma = \sigma_{\rm DM}^{\rm 1D}`, the local one-dimensional dark matter velocity 
  dispersion. This is currently the only option for a variable IMF.

GFM_PREENRICH
  Request that the total metallicity and individual metal abundances be set manually for all gas 
  cells in the simulation, intended to be at high redshift. The time this occurs is set by the 
  ``PreEnrichTime`` parameter, and the values used by the ``PreEnrichAbundanceFile`` contents.

GFM_EXACT_NUMNGB
  Require the exact number of neighbor cells for stellar enrichment and/or local feedback.
  TODO: why is this not degenerate with ``MaxNumNgbDeviationEnrichment=0``?


Galactic-Scale Winds
^^^^^^^^^^^^^^^^^^^^

GFM_WINDS
  Master switch to enable the kinetic (decoupled) wind generation scheme. In particular, all 
  star-forming gas cells have, in addition to a probability for turning into a star particle, 
  a probability for being turned into a wind-phase particle and kicked out of the galaxy. 
  The velocity and thermal energy of new wind particles is set using scalings as determined 
  by further options. The probability is given by the ratio of the desired wind mass to the 
  gas cell mass, where the desired wind mass is given by the stellar mass expected to form 
  multiplied by the mass loading :math:`\eta`, which is itself set using scalings as 
  determined by further options. This proceeds largely as in the original 
  `SH03 model <http://adsabs.harvard.edu/abs/2003MNRAS.339..289S>`_ description (2nd half of 
  paper).

GFM_WINDS_VARIABLE=0
  Determines the primary scaling of the wind launch velocity. If ``GFM_WINDS_VARIABLE=1`` then 
  decoupled sigma winds are used, whereas if ``GFM_WINDS_VARAIBLE=0`` then decoupled halo-mass 
  based sigma winds are used. See the ``VariableWindVelFactor`` parameter for details.

GFM_WINDS_VARIABLE_HUBBLE
  If enabled, the wind launch velocity scales with the Hubble factor. In particular, the 
  :math:`\sigma` as described under the ``VariableWindVelFactor`` parameter is multiplied by 
  a factor of :math:`(H_0/H(z))^{(1/3)}` with the effect of lowering the wind velocity at high 
  redshift, and leaving the :math:`z=0` velocity unchanged.

GFM_WINDS_HUBBLESCALING
  If enabled, the wind energy fraction scales with the Hubble factor.
  Requires the parameter ``WindSuppressionRedshift`` (see for more details).

GFM_WINDS_MASSSCALING
  todo

GFM_WIND_ENERGY_METAL_DEPENDENCE
  todo

GFM_WIND_ENERGY_METAL_DEPENDENCE_TANH
  todo

GFM_WINDS_STRIPPING
  todo
  Requires the parameter ``WindDumpFactor``.

GFM_WINDS_THERMAL
  Allows specification of the thermal component of the newly launched wind particles in terms of 
  the ``ThermalWindFactor`` parameter, which is otherwise zero.

GFM_WINDS_THERMAL_NEWDEF
  Allows specification of the thermal component of the wind as an energy fraction, using the 
  ``ThermalWindFraction`` parameter. The thermal energy for wind particles at launch is otherwise zero.

GFM_BIPOLAR_WINDS=1
  Selects a model for wind directionality. Three options: (i) using "local" (relative to group motion) 
  v-cross-a wind direction (GFM_BIPOLAR_WINDS=1, requires FOF). (ii) using absolute (relative to box) 
  v-cross-a wind direction (GFM_BIPOLAR_WINDS=0). (iii) using isotropic wind direction (default, i.e. 
  without GFM_BIPOLAR_WINDS). All three options can be used independently of GFM_WINDS_VARIABLE.
  Newer: GFM_BIPOLAR_WINDS=3 option launches the wind parallel to the spin of the star-forming gas in 
  the halo.


Cooling and Heating
^^^^^^^^^^^^^^^^^^^

GFM_PRIMORDIAL_RATES
  todo

GFM_COOLING_METAL
  todo

GFM_UVB_CORRECTIONS
  todo

GFM_AGN_RADIATION
  todo.
  A simple ADAF radiative efficiency scaling for advection dominated BH accretion flows can be used 
  to calculate the bolometric AGN luminosity. 
  Can account for the B-band obscuration fraction.
  The maximum size of the AGN sphere is limited to a multiple (currently 2.0) of the virial radius.
  The UV/X background is combined with AGN radiation field in a 4D table such that no switch between 
  different radiation fields is required. Primordial cooling uses linearity of photo-heating and 
  photo-ionisation rates to adjust cooling to correct radiation background which can be superposition 
  of UV/X plus AGN radiation field.
  The bolometric intensity is capped for a cell such that the ionization parameter is at its maximum 
  reliable value: GFM_MAX_IONIZATION_PARAMETER, instead of turning off the AGN radiation altogether 
  for cells with such high ionization parameters.


Additional Output
^^^^^^^^^^^^^^^^^

GFM_STELLAR_PHOTOMETRICS
  todo

GFM_OUTPUT_MASK=1+2+4+...
  todo

GFM_OUTPUT_BIRTH_POS
  todo


Other
^^^^^

GFM_WINDS_LOCAL
  todo.
  Non-local energy driven sigma winds, which are created around stars and not from the ISM.

GFM_STELLAR_FEEDBACK
  todo

GFM_CHECKS
  todo

GFM_DISCARD_ENRICHMENT_GRADIENTS
  If enabled, the advection of **all passive scalars** (independent of type) reverts to the 
  first order advection scheme, ignoring the gradients.

GFM_NORMALIZED_METAL_ADVECTION
  Adds a fiducial "OtherMetals" chemical element (the 10th entry of the ``GFM_Metals`` output vector), 
  so that the chemical abundance vector at interfaces can always be normalized to unity. The total 
  metallicity is also normalized. If this option is enabled, the vector of metals is treated as 
  ``SCALAR_TYPE_SPECIES`` instead of the default ``SCALAR_TYPE_PASSIVE``, meaning that the pure 
  advection of these scalers is modified such that their sum is one. Furthermore, in cases where the 
  FV fluxes of a scalar and mass have opposite sign, or when the FV flux of the scalar is a dominant 
  fraction of the FV mass flux, then the fluxes are reverted to the first order donor cell advection 
  scheme.

GFM_WINDS_SAVE_PARTTYPE=2
  Saves wind particles as a separate particle type when writing out snapshots, instead of being mixed 
  in with stars (PartType4) as by default. The particle type number is given by the definition value 
  (e.g. 2 in this case). **Not fully tested and should be verified before use.**

GFM_DISCRETE_ENRICHMENT
  Restricts the frequency of stellar enrichment events. In particular, allow stars to enrich nearby gas 
  only when (i) the fractional mass loss of the star particle is above some threshold fraction, *or* 
  (ii) the age of the star particle is below some value. If this option is not enabled, star particles 
  enrich their surroundings every timestep they are active, which can computationally expensive, and also 
  implies the frequency of enrichment scales (gets higher) with increasing (better) numerical resolution. 
  Using this option attempts to address those two issues. The threshold values are currently hard-coded 
  in ``GFM/stellar_evolution_main.c`` and are not parameters. They are :math:`\Delta M = 0.0001` and 
  :math:`t_{\rm age} = 0.1 \rm{Gyr}`, and should be changed as needed for a given simulation.


Other GFM-Dependent Modules
---------------------------

See :doc:`modules_gfm_chemtags` for information on::

  GFM_CHEMTAGS
  GFM_SPLITFE
  GFM_SPLITFE_ADDINAGB
  GFM_RPROCESS

See :doc:`modules_gfm_dust` for information on::

  GFM_DUST
  GFM_DUST_DESTMODE=0
  GFM_DUST_SPUTTERING=1


Authors
-------

* Mark Vogelsberger
* Debora Sijacki
* Paul Torrey
* Shy Genel
* Volker Springel
* and others


Usage Policy and Citation
-------------------------

Please contact the authors before using this code for a new project. 
Co-authorship on papers may be requested.

To use aspects of GFM for a new simulation project, including un-modified versions of the 
"Illustris" or "IllustrisTNG" models un-modified, please 
**contact together (at a minimum): Mark, Shy, Volker, Lars, and Annalisa.**

Papers to cite:

* Depends on application.


External Data Files
-------------------

Several components of `GFM` require external binary/HDF5 data files. In particular, photometrics, 
yields, and cooling. For some, multiple versions exist and have been used. Many physical assumptions 
go into each file, and their creation can be documented here.

* `stellar_photometrics.hdf5` (`md5sum: f4bcd628b35036f346b4e47f4997d55e`). Only one version exists. 
  Can be re-created with the following script::

    def makeStellarPhotometricsHDF5_BC03():
        """ Create stellar_photometrics.hdf5 file using BC03 models, as used for Illustris and IllustrisTNG runs.
        Bands: UBVK (Buser U,B3,V,IR K filter + Palomar200 IR detectors + atmosphere.57) in Vega, griz (sdss) in AB
        Requires: http://www.bruzual.org/bc03/Original_version_2003/bc03.models.padova_1994_chabrier_imf.tar.gz
        Produces: 87f665fe5cdac109b229973a2b48f848  stellar_photometrics.hdf5
        Original: f4bcd628b35036f346b4e47f4997d55e  stellar_photometrics.hdf5
          (all datasets between the two satisfy np.allclose(rtol=1e-8,atol=8e-4))
        """
        import numpy as np
        import h5py
        import glob

        filenames1 = sorted(glob.glob("bc2003_hr_m*_chab_ssp.1color")) # m22-m72
        filenames2 = sorted(glob.glob("bc2003_hr_m*_chab_ssp.1ABmag")) # m22-m72

        # linear metallicities (mass_metals/mass_total), not in solar!
        Zvals = [0.0001, 0.0004, 0.004, 0.008, 0.02, 0.05] 
        bandNames = ['U','B','V','K','g','r','i','z'] 

        nAgeVals = 220
        assert len(Zvals) == len(filenames1) == len(filenames2)

        # allocate
        ages = np.zeros(nAgeVals)
        mags = {}
        for bandName in bandNames:
            mags[bandName] = np.zeros( [len(Zvals),nAgeVals] )

        # load BC03 model files
        for i in range(len(Zvals)):
            data1 = np.loadtxt(filenames1[i])
            data2 = np.loadtxt(filenames2[i])

            # verify expected number of rows/ages, and that we process the correct metallicity files
            assert data1.shape[0] == data2.shape[0] == nAgeVals
            with open(filenames1[i], 'r') as f:
                assert "Z=%g" % Zvals[i] in f.read()
            with open(filenames2[i], 'r') as f:
                assert "Z=%g" % Zvals[i] in f.read()

            ages = data1[:,0] - 9.0 # log yr -> log Gyr, same in all files
            mags['U'][i,:] = data1[:,2]
            mags['B'][i,:] = data1[:,3]
            mags['V'][i,:] = data1[:,4]   
            mags['K'][i,:] = data1[:,5]

            mags['g'][i,:] = data2[:,2]
            mags['r'][i,:] = data2[:,2] - data2[:,4]
            mags['i'][i,:] = data2[:,2] - data2[:,5]
            mags['z'][i,:] = data2[:,2] - data2[:,6]

        # write output
        with h5py.File('stellar_photometrics.hdf5', 'w') as f:

            f["N_LogMetallicity"]    = np.array([len(Zvals)], dtype='int32')
            f["N_LogAgeInGyr"]       = np.array([nAgeVals], dtype='int32')
            f["LogMetallicity_bins"] = np.log10(np.array(Zvals, dtype='float64'))
            f["LogAgeInGyr_bins"]    = ages

            for bandName in bandNames:
                f["Magnitude_"+bandName] = mags[bandName]

The exact format of the `Yields/*` files and the actual scripts and original data sources used to make them should be 
better documented, this is a work in progress. 

* `Yields/AGB.hdf5` may be (or symlinked to) any of the following versions:

  * `AGB_Karakas.hdf5` (`md5sum: 924b37d1700d422fc724c440fdf08134`) details needed (used for original Illustris runs).
  * `AGB_Margio.hdf5` (`md5sum: fb8cba2c5965be76e8cc5744d65894aa`)details needed.

* `Yields/Lifetimes.hdf5` (`md5sum: ba13e3b1f0cb32041c3beb7ef9021575`) details needed.

* `Yields/SNIa.hdf5` may be (or symlinked to) any of the following versions:

  * `SNIa_Illustris.hdf5` (`md5sum: d1f792eb624db366757359b8dcd10e4f`) details needed (used for original Illustris runs).
  * `SNIa_nomoto97_w7.hdf5` (`md5sum: 78075f760b224d18f8b4a34ba1bed81c`) details needed.
  * `SNIa_nomoto97_w7_withSplitFe_withEu.hdf5` (`md5sum: 633d980a3b0f103469adbe9c05313b5b`) details needed.
  * `SNIa_nomoto97_wdd2.hdf5` (`md5sum: 6a97a3e9956e76d93c4aa9edb422d3ae`) details needed.

* `Yields/SNII.hdf5` may be (or symlinked to) any of the following versions:

  * `SNII_Illustris.hdf5` (`md5sum: c9b036bbadb2f5705532568041f51bbe`) details needed (used for original Illustris runs).
  * `SNII_imfweighted.hdf5` (`md5sum: a8ca0d7f3e100a2d8065616a660b3e7f`) details needed.
  * `SNII_kobayashi_kobHe_m8to120_concat.hdf5` (`md5sum: 987f9fdaac2f747d5c6aeacb71d7ad05`) details needed.
  * `SNII_kobayashi_kobHe_m8to120_imfrenorm.hdf5` (`md5sum: 1ce3987aeff507a620e3fc313ebbd4d4`) details needed.
  * `SNII_kobayashi_kobHe_m8to120_imfrenorm_withSplitFe_withEu.hdf5` (`md5sum: adaf84ce81f115052738ddfc440595b2`) details needed.
  * `SNII_sn13.hdf5` (`md5sum: 6841a7fabf8bb56ae229306a074ee6f8`) details needed.

The exact format of the `Cooling/*` files and the actual scripts and original data sources (CLOUDY) used to make them should be 
better documented, this is a work in progress. 

* `Cooling/cooling_metal_AGN_Compton.hdf5` (`md5sum: 95455b258cd0deff8536a9c7f55e41b0`) details needed.
* `Cooling/cooling_metal_AGN_Compton_self_shielding_Rahmati12.hdf5` (`md5sum: db36eeada6bb23ec7622f07d4975d54e`) details needed (used for original Illustris runs).
* `Cooling/cooling_metal_UVB.hdf5` (`md5sum: 216689a0a8411bb8357b5a300b66aef5`) details needed.
* `Cooling/cooling_metal_UVB_self_shielding_Rahmati12.hdf5` (`md5sum: 181744af14bbc9778fb0b7ade9f72687`) details needed.
