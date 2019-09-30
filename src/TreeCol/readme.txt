Module: TreeCol

 Description:

  Obtains column density maps of the sky as seen by each cell in the
  simulation volume. Maps are obtained during the walk of the gravitational
  tree - hence the name. Works by recording the properties of the tree nodes
  (i.e position in sky, angular area and column), as they are encountered
  during the tree walk, storing them on a map for later use. (using a
  HEALPix sphere). In principle, the code can make a map of anything that is
  stored in the tree. Currently makes maps of total column, column of H2
  and column of CO (the latter two in conjunction with the SGChem routine).

 Authors:
  Simon Glover (glover@uni-heidelberg.de),
  Paul Clark (paul.clark@astro.cf.ac.uk),
  Rowan Smith (rowanjsmith.astro@googlemail.com),
  Tilman Hartwig (hartwig@iap.fr)

 Acknowledgements:

  The authors thank Anna Schauer for her help with testing and debugging the
  TREE_RAD_VEL option.

 Usage policy:

 Please contact the authors before using this code, so that we can
 avoid any unnecessary project overlap. Co-authorship on publications
 produced with the code is not required, although we are happy to
 actively collaborate with users of the code that are interested in
 doing so.

 Papers to cite:

   We request that users always cite the main code paper for the algorithm:
        Clark, Glover & Klessen 2012, MNRAS, 420, 745

   In addition, we request that the following papers be cited where
   appropriate (i.e. depending on what the code is being used for). NOTE:
   This list will grow as new features are added (by us and other authors):

   1) For first use with present-day chemistry & cloud heating and cooling
        Glover & Clark, 2012, MNRAS, 421, 9

   2) For first use with primordial chemistry and/or with TREE_RAD_VEL option
        Hartwig et al, 2015, ApJ, 799, 114; Hartwig et al, 2015, MNRAS, 452, 1233
