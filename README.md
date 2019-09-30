
AREPO
=====

AREPO is a massively parallel code for hydrodynamical cosmological
simulations. It is a flexible code that can be applied to a variety of
different types of simulations, offering a number of sophisticated
simulation algorithms.

A full account of the numerical algorithms employed by the code is given
in the original code paper and subsequent code development papers. 
Instructions for both usage and development of the code are given in the 
included code documentation.


Documentation
=============

The documentation is built automatically from source after each commit and is available here:

* [Arepo Online Documentation](http://www.illustris-project.org/w/arepo_docs/)

It covers both user and development oriented information. 
A [PDF Version](http://www.illustris-project.org/w/arepo_docs/Arepo.pdf) also exists.

To build a local copy of the documentation:

```
  doxygen
  pip install sphinx breathe
  sphinx-build -b html documentation/ doxygen/sphinx_html/

  sphinx-build -b latex documentation/ doxygen/sphinx_latex/
  (cd doxygen/sphinx_latex && pdflatex Arepo.tex)
```


Developers
==========

See [DEVELOPERS](DEVELOPERS).

# dissipative
