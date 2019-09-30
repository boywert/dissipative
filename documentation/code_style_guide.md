
Coding Guidelines for AREPO
===========================

Code extensions
---------------

* Non-standard code extensions should always be written such that they
  can be switched off if not needed, and have no side effects on
  existing code. Normally this means that they have to be enclosed in
  conditional compilation precompiler statements (``#ifdef``), especially
  if variables in the global structures of the code need to be
  allocated for the extension. However, if the extension's execution
  can also be controlled at run-time by simple variables, then
  consider introducing a parameterfile variable to control the
  extension. In general, the number of Makefile symbols (Config.sh 
  options) to control conditional compilation should be kept to a minimum.

* Do not place any substantial piece of code belonging to your
  extension into existing functions of the code. Write your own
  functions for the code extension, and only place a function call (if
  needed bracketed by an ``#ifdef``) into the appropriate place of the
  primary code. Also, place your extension functions into separate
  source files.


General code-style principles
-----------------------------

* Code formatting: Try to be consistent with the code formatting of
  the main code, which is more or less GNU-style, with a few small
  differences. You can run the indent-command with the options::
  
    indent -gnu -npsl -npcs -nbs -nsaf -nsai -nsaw -nprs -bap -pmt -l110 *.c

  on your source file(s) to make the indention consistent.

* Name functions all in lower case as a "command" that is descriptive
  of what the function does. Words are separated by underscores, e.g.
  ``calculate_normal_vector_for_triangle(...)``.
  For all functions, input arguments should come first, output arguments
  last.

* Global variables (which use should be kept to a minimum) start with
  an upper case character. They are nouns, with words separated by
  mixed lower/upper case characters, and contain no underscores. Use
  abbreviations that are clear, e.g. ``NumForceCalculations``.

* Local variables start with lowercase, and should have descriptive
  names too, except for simple loop iterators and the like. Try to
  narrow the scope of a variable as much as possible, for example by
  declaring them inside a block where they are needed. Declaration and
  initialization should be combined in one command where possible, for 
  example::

    int n = get_particle_count();

  instead of::

    int n;
    n = get_particle_count();

* Avoid repition of code, write a function instead. Break up long
  functions into smaller more managable pieces.

* Preprocessor macros that have arguments should be avoided whenever
  possible. If needed, their names should be fully capitalized.

* Magic numbers (including numerical constants) in the code should be
  avoid and instead be replaced by a symbolic constant declared with
  #ifdef in a header file and made all uppercase.

* Address all warnings emitted by the compiler when compiler with
  ``-Wall``. Unless there are good reasons, the code should compile
  without any warning.


Code comments
-------------

Include consistent commenting in your code. The meaning of all
global variables should be commented where they are introduced,
ideally with doxygen syntax, for example::

    int MyGlobalCount;   /*!< counts the number of timesteps */

Function should be preceeded by a brief explanation of what the
function does, including warnings/instruction, if any, about how the
function may be used. For example::

    /*!  
     * Insert the point P[i] into the current tessellation. Start
     * the search at Delaunay triangle DT[t], and return the index
     * of the last accessed triangle.
     */
    int insert_point(int i, int t)  
    {
      ...
    }

You do not need to state here *how* the function achieves what it
does, but this can be stated if appropriated in comments in the
function body. There, avoid superfluous comments that just reflect
what's obvious from the code anyway, like::

    do_domain_decomposition();   /* call domain decomposition */

Instead, focus on comments that help one to quickly understand/check
what the code tries to do at an algoritmic level. If complicated
formulae are implemented, try to include in a comment a reference to
the equation that is implemented. Avoid c++ style comments in C-code.
