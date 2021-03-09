module fft_mod
==============

Overview
--------

Performs simultaneous fast Fourier transforms (FFTs) between real grid space and complex Fourier space.

.. container::

   This routine computes multiple 1-dimensional FFTs and inverse FFTs. There are 2d and 3d versions between type real
   grid point space and type complex Fourier space. There are single (32-bit) and full (64-bit) versions.
   On Cray and SGI systems, vendor-specific scientific library routines are used, otherwise a user may choose a NAG
   library version or stand-alone version using Temperton's FFT.

| 

Other modules used
------------------

.. container::

   ::

      platform_mod
           fms_mod
         fft99_mod

Public interface
----------------

.. container::

   ::

      use fft_mod [, only:  fft_grid_to_fourier,
                            fft_fourier_to_grid,
                            fft_init,
                            fft_end ]

   fft_grid_to_fourier:
      Given multiple sequences of real data values, this routine computes the complex Fourier transform for all
      sequences.
   fft_fourier_to_grid:
      Given multiple sequences of Fourier space transforms, this routine computes the inverse transform and returns the
      real data values for all sequences.
   fft_init:
      This routine must be called to initialize the size of a single transform and setup trigonometric constants.
   fft_end:
      This routine is called to unset the transform size and deallocate memory.

| 

Public data
-----------

.. container::

   None.

Public routines
---------------

a. .. rubric:: Fft_grid_to_fourier
      :name: fft_grid_to_fourier

   ::

      fourier = fft_grid_to_fourier ( grid )

   **DESCRIPTION**
      Given multiple sequences of real data values, this routine computes the complex Fourier transform for all
      sequences.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``grid``                                                  | Multiple sequence of real data values. The first          |
      |                                                           | dimension must be n+1 (where n is the size of a single    |
      |                                                           | sequence).                                                |
      |                                                           | [real(R4_KIND), dimension(:,:)]                           |
      |                                                           | [real(R8_KIND), dimension(:,:)]                           |
      |                                                           | [real(R4_KIND), dimension(:,:,:)]                         |
      |                                                           | [real(R8_KIND), dimension(:,:,:)]                         |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``fourier``                                               | Multiple sequences of transformed data in complex Fourier |
      |                                                           | space. The first dimension must equal n/2+1 (where n is   |
      |                                                           | the size of a single sequence). The remaining dimensions  |
      |                                                           | must be the same size as the input argument "grid".       |
      |                                                           | [complex(R4_KIND), dimension(lenc,size(grid,2))]          |
      |                                                           | [complex(R8_KIND), dimension(lenc,size(grid,2))]          |
      |                                                           | [complex(R4_KIND),                                        |
      |                                                           | dimension(lenc,size(grid,2),size(grid,3))]                |
      |                                                           | [complex(R8_KIND),                                        |
      |                                                           | dimension(lenc,size(grid,2),size(grid,3))]                |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **NOTE**
      The complex Fourier components are passed in the following format.

      ::

                 fourier (1)     = cmplx ( a(0), b(0) )
                 fourier (2)     = cmplx ( a(1), b(1) )
                     :              :
                     :              :
                 fourier (n/2+1) = cmplx ( a(n/2), b(n/2) )

      where n = length of each real transform

b. .. rubric:: Fft_fourier_to_grid
      :name: fft_fourier_to_grid

   ::

      grid = fft_fourier_to_grid ( fourier )

   **DESCRIPTION**
      Given multiple sequences of Fourier space transforms, this routine computes the inverse transform and returns the
      real data values for all sequences.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``fourier``                                               | Multiple sequence complex Fourier space transforms. The   |
      |                                                           | first dimension must equal n/2+1 (where n is the size of  |
      |                                                           | a single real data sequence).                             |
      |                                                           | [real(R4_KIND), dimension(:,:)]                           |
      |                                                           | [real(R8_KIND), dimension(:,:)]                           |
      |                                                           | [real(R4_KIND), dimension(:,:,:)]                         |
      |                                                           | [real(R8_KIND), dimension(:,:,:)]                         |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``grid``                                                  | Multiple sequence of real data values. The first          |
      |                                                           | dimension must be n+1 (where n is the size of a single    |
      |                                                           | sequence). The remaining dimensions must be the same size |
      |                                                           | as the input argument "fourier".                          |
      |                                                           | [complex(R4_KIND), dimension(leng1,size(fourier,2))]      |
      |                                                           | [complex(R8_KIND), dimension(leng1,size(fourier,2))]      |
      |                                                           | [complex(R4_KIND),                                        |
      |                                                           | dimension(leng1,size(fourier,2),size(fourier,3))]         |
      |                                                           | [complex(R8_KIND),                                        |
      |                                                           | dimension(leng1,size(fourier,2),size(fourier,3))]         |
      +-----------------------------------------------------------+-----------------------------------------------------------+

c. .. rubric:: Fft_init
      :name: fft_init

   ::

      call fft_init ( n )

   **DESCRIPTION**
      This routine must be called once to initialize the size of a single transform. To change the size of the transform
      the routine fft_exit must be called before re-initialing with fft_init.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``n``                                                     | The number of real values in a single sequence of data.   |
      |                                                           | The resulting transformed data will have n/2+1 pairs of   |
      |                                                           | complex values.                                           |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

d. .. rubric:: Fft_end
      :name: fft_end

   ::

      call fft_end 

   **DESCRIPTION**
      This routine is called to unset the transform size and deallocate memory. It can not be called unless fft_init has
      already been called. There are no arguments.

Data sets
---------

.. container::

   None.

Error messages
--------------

.. container::

   **Error in fft_grid_to_fourier**
      fft_init must be called
      The initialization routine fft_init must be called before routines fft_grid_to_fourier.
   **Error in fft_grid_to_fourier**
      size of first dimension of input data is wrong
      The real grid point field must have a first dimension equal to n+1 (where n is the size of each real transform).
      This message occurs when using the SGI/Cray fft.
   **Error in fft_grid_to_fourier**
      length of input data too small
      The real grid point field must have a first dimension equal to n (where n is the size of each real transform).
      This message occurs when using the NAG or Temperton fft.
   **Error in fft_grid_to_fourier**
      float kind not supported for nag fft
      32-bit real data is not supported when using the NAG fft. You may try modifying this part of the code by
      uncommenting the calls to the NAG library or less consider using the Temperton fft.
   **Error in fft_fourier_to_grid**
      fft_init must be called
      The initialization routine fft_init must be called before routines fft_fourier_to_grid.
   **Error in fft_fourier_to_grid**
      size of first dimension of input data is wrong
      The complex Fourier field must have a first dimension equal to n/2+1 (where n is the size of each real transform).
      This message occurs when using the SGI/Cray fft.
   **Error in fft_fourier_to_grid**
      length of input data too small
      The complex Fourier field must have a first dimension greater than or equal to n/2+1 (where n is the size of each
      real transform). This message occurs when using the NAG or Temperton fft.
   **Error in fft_fourier_to_grid**
      float kind not supported for nag fft
      float kind not supported for nag fft 32-bit real data is not supported when using the NAG fft. You may try
      modifying this part of the code by uncommenting the calls to the NAG library or less consider using the Temperton
      fft.
   **FATAL in fft_init**
      attempted to reinitialize fft
      You must call fft_exit before calling fft_init for a second time.
   **Error in fft_end**
      attempt to un-initialize fft that has not been initialized
      You can not call fft_end unless fft_init has been called.

References
----------

.. container::

   #. For the SGI/Cray version refer to the manual pages for DZFFTM, ZDFFTM, SCFFTM, and CSFFTM.
   #. For the NAG version refer to the NAG documentation for routines C06FPF, C06FQF, and C06GQF.

| 

Compiler specifics
------------------

.. container::

   None.

| 

Precompiler options
-------------------

.. container::

   -D **NAGFFT**
      -D NAGFFT On non-Cray/SGI machines, set to use the NAG library FFT routines. Otherwise the Temperton FFT is used
      by default.
   -D **test_fft**
      Provides source code for a simple test program. The program generates several sequences of real data. This data is
      transformed to Fourier space and back to real data, then compared to the original real data.

| 

Loader options
--------------

.. container::

   On SGI machines the scientific library needs to be loaded by linking with:

   ::

              -lscs

   If using the NAG library, the following loader options (or something similar) may be necessary:

   ::

              -L/usr/local/lib -lnag

Test PROGRAM
------------

.. container::

   None.

| 

Notes
-----

.. container::

   The routines are overloaded for 2d and 3d versions. The 2d versions copy data into 3d arrays then calls the 3d
   interface.
   On SGI/Cray machines:
   There are single (32-bit) and full (64-bit) versions. For Cray machines the single precision version does not apply.
   On non-SGI/CRAY machines:
   The NAG library option uses the "full" precision NAG routines (C06FPF,C06FQF,C06GQF). Users may have to specify a
   64-bit real compiler option (e.g., -r8).
   The stand-alone Temperton FFT option works for the real precision specified at compile time. If you compiled with
   single (32-bit) real precision then FFT's cannot be computed at full (64-bit) precision.

| 
