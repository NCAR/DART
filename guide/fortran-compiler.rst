##################
Fortran90 compiler
##################

The DART software is written in standard Fortran 90, with no
compiler-specific extensions. It has been compiled and run with several
versions of each of the following:

- `GNU Fortran Compiler (known as "gfortran") <http://gcc.gnu.org/fortran>`_ (free)
- `Intel Fortran Compiler for Linux and OSX <http://software.intel.com/en-us/intel-composer-xe>`_
- `IBM XL Fortran Compiler <http://www-01.ibm.com/software/awdtools/fortran/>`_
- `Portland Group Fortran Compiler <http://www.pgroup.com/>`_
- `Lahey Fortran Compiler <http://www.lahey.com/>`_
- `NAG Fortran compiler <https://www.nag.com/nag-compiler>`_
- `PathScale Fortran compiler <https://en.wikipedia.org/wiki/PathScale>`_

Since recompiling the code is a necessity to experiment with different
models, there are no DART binaries to distribute. 

If you are unfamiliar with Fortran and/or wonder why we would choose this language, 
see the :doc:`Why Fortran? <./dart-design-philosophy>` discussion for more information.

