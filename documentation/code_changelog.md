
Change Log
==========

Meant to be a much coarser representation of code changes over time than the svn/git repository. In 
particular, substantial additions and/or changes to code algorithms can be noted. Changes affecting 
default physical behavior or output should be included, as should changes affecting future code 
development. Module-specific improvements or fixes are not here.

Mar 16, 2016
  * Refactored I/O, migrate default fields to **new I/O scheme**, add support for physical units in 
    output. (abauer)

Feb 16, 2016
  * Removed ``USE_RANDOM_NUMBERS`` option (this is now the default), using an indexed look-up table 
    is gone. (vspringel)

Dec 28, 2015
  * Added asynchronous communication pattern ``GENERIC_ASYNC``. (vspringel)

Nov 25, 2015
  * **Changed the default time integration scheme to Runge Kutta and the default gradient estimate to 
    least squares fit gradients.** The switches ``RUNGE_KUTTA`` and ``GRADIENTS_LEAST_SQUARE_FIT`` are 
    now obsolete. Instead, there are two new switches that have to be enabled explicitely to get the 
    old behaviour: ``MUSCL_HANCOCK`` and ``GRADIENTS_GREEN_GAUSS``. Note, however, that using either 
    of those switches will reduced the order of the hydro scheme to first order. So it is strongly 
    recommended to only use them for testing purposes. (rpakmor)

Sep 16, 2015
  * Migration from SVN to git/bitbucket.

June 16, 2015
  * **Revision of the generic communication pattern for the parallel tree walks.** The primary and secondary 
    walks are now moved to separate static functions in each file that uses the routines. This allows the 
    communication routines to more flexibly handle memory management, and leads to a further simplification 
    of using the routines, which is now effectively a one-liner. The specification of ``BufferSize`` is not 
    necessary any more. Instead, the code will automatically make a close to optimal use of the still 
    available memory. Also, out of memory situations mid-way into the communication ("can't allocate 
    DataGet issue") are safely avoided. The new communication scheme is now also implemented for the 
    tree-based gravity, simplifying the associated code. (vspringel)

Apr 1, 2015
  * Change to new communication structure, fixes to most code parts necessary. (vspringel)

Mar 23, 2015
  * Moved to new determination of long-range PM timestep (vspringel)
  * The old scheme, with the parameter ``MaxRMSDisplacementFac``, can still be used if activated 
    with ``LEGACY_DISPLACEMENT_CONSTRAINT``.

Feb 22, 2015
  * Moved to FFTW3 and a home-grown parallel FFT: only the serial routines from FFTW3 are used. (vspringel)
  * For zoom-simulations, ``PM_ZOOM_OPTIMIZED`` should be enabled, which corresponds to the old communication 
    strategy. For periodically sampled boxes, one can disable this option for an alternative communication 
    scheme, which is then a bit more efficient. The PM_ZOOM_OPRIMIZED routines should essentially be about 
    equally fast as the old ones.
  * There is also a new option ``FFT_COLUMN_BASED``, which removes the scaling barrier of a slab-decomposed 
    parallel FFT. However, the intrinsic speed of this parallel FFT is slower because more transpose 
    operations have to be done. One a very large number of MPI ranks is used, it will however be more 
    efficient to use this one, also for reasons of memory balance.
  * All these new PM routines are hybrid parallelized with OpenMP.

Feb 20, 2015
  * Direct summation gravity for small particle numbers. (vspringel)

Feb 16, 2015
  * This large commit implements: (vspringel/rpakmor)

    * a new organization of the timestepping for gravity and hydro, so that **hydro timesteps can be 
      decoupled and nested inside gravity timesteps** if desired.
    * there is a new Runge-Kutta integration for the fluxes that can get away with one mesh 
      construction per step
    * related to this, there is also an alternative least square gradient estimator which keeps 
      better accuracy for distorted meshes
    * there is also a new **hierarchical timestepping scheme for gravitational forces**
    * gravitational softening is now organized via softening types (there can be a different number of 
      those than particle types). Tree nodes can now be softened, in which case their mass components 
      with different softenings are evaluated separately.

Nov 19, 2013
  * Fixed a bug that affected various SUBFIND group properties:
    gas and stellar masses, metallicities, sizes (radii),
    stellar photometry, all BH properties (mass and Mdot). (sgenel)

Aug 1, 2013
  * Releasing the restriction to 6 particle types using the new switch NTYPES (sgenel)

Jan 7, 2013
  * Reorganized cooling_and_starformation() and get_starformation_rate() such that: (sgenel)

    #. non-SF (too hot) gas above the density threshold will cool also with metals, not just with 
       primordial cooling
    #. SphP.Ne of such gas will not be based on the model hot-phase temperature, because the 
       multi-phase model doesn't apply to such gas.

  * Now all cells (except for dense cold cells) first cool with cool_cell(), and only then the decision re: 
    star-formation is made based on density and temperature.

Nov 29, 2012
  * Electron abundance (IO_NE) in snapshot files is changed for star-forming gas: it is now calculated based 
    on the effective temperature (Utherm), not the hot-phase temperature of the multi-phase medium anymore,
    which also makes it consistent with the neutral abundance IO_NH. (sgenel)

Sep 18, 2012
  * For the ``ENFORCE_JEANS_STABILITY_OF_CELL`` option: the internal energy of the cell is updated 
    consistently with the new value of the pressure floor. (fmarinacci)

Aug 23, 2012
  * Significant change in the way the export of particles (in particular with respect to their tree 
    node indices) is handled during the gravitational tree walks. (vspringel)
  * Instead of allocating a block of slots for such tree nodes as a fixed-length integer array (given 
    by NODELISTLENGTH), the node indices are now stored in a separate buffer such that the number of 
    nodes associated with one particle export is fully flexible. This avoids wasting buffer space and 
    communicating unneeded data, as well as avoiding multiple exports of the same particle (if more 
    nodes need to be accessed than possible for the old fixed-length buffer).
  * Eventually, all other tree walks (there are quite a few!) should be changed to the new scheme, which 
    is simpler and more efficient.

Aug 15, 2012
  * Removed PartAllocFactor (and the other AllocFactors for BHs, Stars, Tracers). The storage needed 
    will be sized automatically and adjusts to what the simulation requires (vspringel)

Aug 1, 2012
  * Large code cleanup, e.g. src/ folder created. (abauer)

July 28, 2012
  * Removed some extremely subtle difference between collective and serial subfind code, such that the 
    agreement between the results is now essentially perfect. In rare cases, it can 
    however nevertheless happen that the gravitational potential calculated by the two codes can have 
    very small differences due to numerical round-off. In even much rarer cases, this can lead to a 
    difference in the final particle order, or even in the substructure length. (vspringel)

July 2, 2012
  * Changed the logic of ``ENFORCE_JEANS_STABILITY_OF_CELLS`` and ``SOFTEREQS``:
    ENFORCE_JEANS_STABILITY_OF_CELLS is now enforced as the last step in set_pressure_of_cell(),
    directly on the pressure, instead of affecting the internal energy, i.e. a direct "pressure floor".
    It also operates now in conjunction with various equations of state, and generates a warning
    in case it is compiled together with any such equation of state, to make sure the user is aware.
    SOFTEREQS is now enforced in sfr_eEOS.c as part of the sub-grid model, affecting the internal energy,
    instead of just being applied on the pressure. (sgenel)

July 2, 2012
  * Bugfixes in collective subfind, assume broken until now. (vspringel)

June 29, 2012
  * Replaced the determination of the top-level tree with a version that scales much better. (vspringel)

June 21, 2012
  * Added the variable P[i].SofteningType. This variable is used only if ``REFINEMENT_HIGH_RES_GAS`` is on 
    and ``INDIVIDUAL_GRAVITY_SOFTENING`` is not defined. It is usually equal to P[i].Type unless a star particle 
    forms from a cell that does not belong to the high-res region. In the previous implementation this star 
    particle was stored as a type = 3 particle to give it a large gravitational softening but this caused problems 
    to the GFM module that assumes that all star particles are of type = 4. Now all stars are stored as type = 4 
    and the gravitational softening is set according to P[i].SofteningType in the function 
    get_softening_of_particle(). (fmarinacci)

June 2, 2012
  * Bugfix Sourceterms: Since the big rewrite of the trees, the list of active particles is already updated in 
    find_next_sync_point_and_dritf(). Therefore, all source terms at the end of the timestep have to be applied 
    BEFORE the call to find_next_sync, which went wrong for do_mhd_source_terms_divb (MHD) and
    set_non_standard_physics_for_current_time (including COOLING: UV Background, WINDTUNNEL & NOH_PROBLEM). (rpakmor)

May 29, 2012
  * Prepared a first alternative top-level tree construction, which should
    be faster on large core counts because it avoids collective operations
    for disjoint branches of the tree. To this end, the intracommunicator is
    split hierachically and intercommuncation is used as well. (vspringel)
  * New top-level tree construction now default as it scales better to high numbers of MPI tasks.

May 25, 2012
  * Changed the radius at which the jeans temperature floor is computed to the softening lenght of the 
    gas cell. (fmarinacci)
  * Added option ``ENFORCE_JEANS_STABILITY_OF_CELLS_EEOS`` to ensure Jeans stability for two-phase gas. 
    (mvogelsberger)

May 10, 2012
  * Modified the optimizing of the work-load balancing such that memory-load balance is also taken 
    into account. The latter should never get worse than about ~2 even in extreme situations. (vspringel)

May 4, 2012
  * Some **large architectural changes have been made**, in particular: (vspringel)

    * a **separate neighbour search tree** for the gas has been introduced (most code parts need fixes)
    * this tree is dynamically updated every step, using shrink wrap size information (TreeDomainUpdateFrequency parameter removed)
    * the gravity tree is reconstructed as optimal tree every timestep, and is balanced on a a per timebin basis
    * the domain decomposition needs not be repeated before a tree construction any more
    * improved balancing algorithm in the domain decomposition seeks balance of several time levels simultaneously
    * the peak memory consumption is reduced because mesh and gravity tree do not need to be stored simultaneously
    * also, file size of restart files reduced because they do not contain the mesh and the gravity tree any more
    * the neighbor finding tree walk has been accelerated
    * first code fragements try out SSE instructions available on Intel/AMD chips

Mar 15, 2012
    * Added ``TOLERATE_WRITE_ERROR``, which only gives warning with node information if write fails. (mvogelsberger)
    * Added user-requested restartfile writing: "touch restart" in outputdir -> write restartfile and continue run. 
    * Added option ``MULTIPLE_RESTARTS`` for incremental restart files.

Mar 10, 2012
  * Added option SUBBOX_SNAPSHOTS to write out subbox of simulation volume with higher output frequency. (mvogelsberger)

Jan 1, 2012
  * Added TREECOOL file according to Faucher-Giguere, most recent December version which gives WMAP-7 
    compatible reionization. (mvogelsberger)

Dec 29, 2011
  * Modified FOF and SUBFIND very substantially, in order to run more efficiently. **Snapshots will now be stored 
    automatically in (sub)group order** in case FOF or SUBFIND is enabled, respectively. Estimated DM densities can 
    be included in the right order in the snapshot files if desired. Working by ~Mar 2012. (vspringel)

Dec 4, 2011
  * Added MPI/Threading hybrid code for gravity loops and related treewalks in NUM_THREADS (mvogelsberger/vspringel)

Apr 27, 2011
  * Changed loading of IC files so that they can contain any field appropriate for the particle type (e.g., metallicity 
    for gas and stars, ages for stars) and it will be read into the arrays. Any fields not present are initialized to 
    zero, along with an info message. This only works for HDF5 files. Also added capability to only load certain particle 
    types, which is used when making projections to avoid useless building of tree for non-gas particles. (pjonsson)

Apr 10, 2011
  * Added code to impose MPI-task binding on a predefined set of cores. (vspringel)

Feb 20, 2011
  * Introduced functions for different refinement and derefinement criteria, which are selected by the new integer parameters 
    ``RefinementCriterion`` and ``DerefinementCriterion``. Added ``INDIVIDUAL_GRAVITY_SOFTENING`` option. (vspringel)

Jan 7, 2011
  * Femoved the ``REGULARIZE_MESH_WITH_INVERSE_ZELDOVICH`` option, and the code that goes with it (was anyway broken at 
    the moment). This feature never worked satisfactorily well in all situations, and was also expensive to calculate. 
    Hence it appears better to focus on refinement/derefinement instead of trying to stear the mesh motion in a highly 
    non-trivial way. (vspringel)

Jan 3, 2011
  * Added TVD slope limiter (option: ``TVD_SLOPE_LIMITER``). This will hopefully
    improve results and we might think about enabling it by default in the future. (rpakmor)

Jan 2, 2011
  * Modified opening critertion for nodes slightly: they are now always opened if the maximal softening for any node 
    particle/cell is larger than the target particle's softening, and if the distance is smaller than this softening 
    (i.e. the node would be softened). In this way, potential large forces still hidden in the node due to small 
    softening will still be seen. *Note: a similar thing is also done in gadget*. (vspringel)

Dec 16, 2010
  * Eliminated drand48() calls in favor of the gsl random number generator. The generated/applied random numbers should 
    now be the same independent of whether or not an interruption with a restart file happened. (vspringel)

Dec 1, 2010
  * Subfind integration started. (vspringel)

Nov 12, 2010
  * Modified the build-strategy of the code: the makefile options are now set in a separated file "Config.sh", which 
    should be created from a template file "Template-Config.sh" (only the latter is under revision control). (vspringel)

Oct 3, 2010
  * "_new" image projection routines. (vspringel)

Sep 11, 2010
  * Modified the treatmeant of refinement and derefinement of cells. The code now calculates the involved volume fractions 
    of the cells for the redistribution of the conserved quantities exactly instead of doing this approximately. This is 
    done by a separate construction of the relevant voronoi cells, such that this approach is still quick. In order to 
    facilitate this, the voronoi_2d/3d functions have been generalized to do their work on a mesh that is discribed by the 
    tessellation structure, which is passed to the functions as a pointer. This structure in turn contains the base pointers 
    of the various tables that describe the meshes. (vspringel)

Aug 25, 2010
  * Modified the slope limiting of the velocity gradient. The maximum allowed change between centre and sufrace of a cell is 
    now 1.0*csnd instead of 0.5*csnd, and the limiting occurs for the vx, vy, vz gradients independently. This reduces the 
    numerical diffusivity that may be introduced by the slope limiting. (vspringel)

July 15, 2010
  * Modified entropy treatment slightly. Instead of M * A, the entropy variable is now taken as M * ln(A). For both 
    formulations, there is conservation law, because DA/Dt = 0 and also DlnA/Dt = 0. However, the formulation Entropy 
    = M lnA is perferred, because this should give a more physical result for adding together gas of different specific 
    entropy. This is because MlnA is directly proportional to the thermodynamic entropy associated with this gas element. 
    (vspringel)

Jul 14, 2010
  * Changed default cooling of gas above two phase medium density threshold. Now the gas is only cooled on the slow relaxation 
    time scale if the ``SLOW_RELAX_TO_EOS`` is set, i.e. the default behaviour was changed! (mvogelsberger)

Jun 20, 2010
  * New ``REGULARIZE_MESH_FACE_ANGLE`` option for MaxFaceAngle regularization scheme. (vspringel)

Jun 18, 2010
  * Re-enabled additional condition on refinement eligibility that restricts the refinement to cells that are sufficiently 
    non-distorted. At present, the condition d < 3.0 * All.CellShapingFactor * cellrad is demanded. Without this, it could 
    occasionally happen that a cell that is split is right away split again, which could then cause inaccuracies in
    the hydro. (vspringel)

Jun 12, 2010
  * Voronoi dynamic update scheme. (vspringel)

May 31, 2010
  * Modified how the radiative cooling is integrated into the timestepping. The new scheme places the cooling after the timestep 
    has fully completed so that it can already work with the updated primitive variables. This prevents that small errors in the 
    temperature evolution can occur if the cooling is weak and the entropy scheme is used. However, the spawning of new stars 
    still needs to be done directly before the mesh construction, it has hence been decoupled from the function that does the 
    radiative cooling.

Mar 19, 2010
  * Implemented a new way to treat the fluxes for faces that are present on more than one CPU. Instead of calculating those twice 
    and applying the result only in a one-sided fashion, they are only calculated once, normally on the CPU that holds the side 
    with the lower ID. In this way, strict conservation of more safely insured because the fluxes for the two independent 
    calculations may differ by round-off error.
  * Also, put in a conservative limiter for the maximum allowed flux over a face to prevent that cells with negative mass may 
    appear under conditions of too poor timestepping. The fraction of faces that receive such a reduction of the flux is monitored.

Feb 20, 2010
  * Major change in treatment of tree updates and in the coupling of gravity to the hydro in the cosmological case. (vspringel)

Dec 23, 2009
  * Added update of UV background. (vspringel)

May 26, 2009
  * AREPO code base split from GADGET-3. (vspringel).