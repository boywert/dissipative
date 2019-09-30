
Explicit Output List
====================

If the ``OutputListOn`` parameter is set to 1, then the user can explicitly specify times 
for the simulation to output. These times are input as a separate text file, whose name 
is specified in the ``OutputListFilename`` parameter. Each line requests one code output, 
and can be specified with one or more numbers per line. The first value (mandatory) 
specifies the requested simulation time for the output (scale factor in the case of 
comoving integrations). The second value (optional), is used as the DumpFlag, described 
below.


DumpFlag
--------

This short documents summarizes the functioning of the DumpFlag, i.e. of the second column input (optional) which can be specified in file given as OutputListFilename.
The flag can assume a few values, implying a different combination of output files:

* 0 - no snapshots and no halo catalogs
* 1 - both snapshots and halo catalogs (default)
* 2 - only snapshots 
* 3 - mini snapshots and halo catalogs
* 4 - only halo catalogs

Note: in this implementation, the naming scheme is kept fixed, i.e. for example full and mini snapshots have the same name roots (``snap_00*/snap_00*.*.hdf5``) but simply differ in size.

Note: in all cases, the counter for the dump keeps running (see example below)  

Note: as of 2016/03/21, option ``DumpFlag == 4`` might conflict with the dump of the SUBBOXES: to be checked and fixed.

Example
^^^^^^^

The following entries in OutputListFile::

  0.0079  1
  0.0081  1
  0.0083  0
  0.0085  2
  0.0088  3
  0.0090  3
  0.0093  4
  0.0095  4
  0.0098  4
  0.01    1

produce the following dumps

* snapdir_000
* snapdir_001
* snapdir_003
* snapdir_004
* snapdir_005
* snapdir_009

* groups_000
* groups_001
* groups_004
* groups_005
* groups_006
* groups_007
* groups_008
* groups_009

