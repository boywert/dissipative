
AMR
===

Provides an AMR mesh instead of the Voronoi mesh.


Usage
-----

TODO


Additional Parameters
---------------------

* ``Any`` todo.


Additional Config.sh Options
----------------------------

AMR
  This switch changes the mesh to an adaptive mesh refinement (AMR) mesh

AMR_CONNECTIONS
  Compute the Mesh.DC connection information in case of an AMR mesh, like for a Voronoi mesh. 
  For each cell Mesh.DC contains a linked list over all neighboring cells.

AMR_GRADIENTS
  Use the AMR gradient estimator instead of the least square fitting method. 
  This gradient estimator computes gradients using finit differences from neighboring cells.

AMR_REDUCE_DOMAIN_DECOMPOISTION
  Only do a domain decomposition when a snapshot is written. This is especially usefull for static meshes.

AMR_REMAP
  Remap Voronoi/SPH ICs on an AMR mesh. This option is currently broken/partiually implemented.


Authors
-------

  * Andreas Bauer (andreas.bauer@h-its.org)
  * Ruediger Pakmor (ruediger.pakmor@h-its.org)


Usage Policy and Citation
-------------------------

Please contact the author before using this code for a new project. Co-authorship on
a first paper may be required.

Papers to cite:

  * No papers yet.
