##################
Fortran compiler
##################

DART is written in Fortran. We recommend using a modern Fortran compiler
that supports Fortran 2008 or later. The code has been tested with several
compilers.

- `GNU Fortran Compiler <http://gcc.gnu.org/fortran>`__ known as "gfortran" (free)
- `Intel Fortran Compiler for Linux <https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html#gs.mnuvoc>`_
- `Nvidia Fortran Compiler  <https://developer.nvidia.com/hpc-sdk>`_
- `NAG Fortran compiler <https://nag.com/fortran-compiler/>`_
- `Cray Compiler Environment <https://cpe.ext.hpe.com/docs/24.11/cce/index.html>`__ (available on Cray systems)

Since recompiling the code is a necessity to experiment with different
models, there are no DART binaries to distribute. 

If you are unfamiliar with Fortran and/or wonder why we would choose this language, 
see the :doc:`Why Fortran? <./dart-design-philosophy>` discussion for more information.

