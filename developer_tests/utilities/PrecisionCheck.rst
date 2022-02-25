PROGRAM ``PrecisionCheck``
==========================

Overview
--------

| This is a self-contained program to explore the interaction between the compiler options to 'autopromote' variables
  from one precision to another and the intrinsic F90 mechanism for getting consistent behavior without relying on
  autopromotion - namely, the ``SELECT_INT_KIND()`` and ``SELECT_REAL_KIND()`` functions. The most portable code
  explicity types the variables to avoid relying on compiler flags. The core DART code abides by these rules; some
  pieces that are derived from dynamical models may have original code fragments.
| All that is required is to compile the single file and run the resulting executable. There are no required libraries -
  any F90 compiler should have no trouble with this program. There is no input of any kind.
| You are encouraged to view the source code. It's pretty obvious what is being tested.

Examples
--------

The following examples have differences from the default configuration highlighted in boldface. You are strongly
encouraged to test your compiler and its autopromotion options. The Absoft compiler actually does what I consider to be
reasonable and logical (as long as you know that "-dp" means **d**\ emote **p**\ recision). Many other compilers are
surprising.

PowerPC chipset : Absoft Pro Fortran 9.0
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. container:: unix

   ::

      [~/DART/utilities] % f90 PrecisionCheck.f90
      [~/DART/utilities] % ./a.out
       
       This explores the use of the intrinsic SELECTED_[REAL,INT]_KIND() functions
       and the interplay with the compiler options. You are encouraged to use the
       "autopromotion" flags on your compiler and compare the results.
       
      ----------------------------------------------
       "integer"
       DIGITS    =   31
       HUGE      =   2147483647
       KIND      =   4
      ----------------------------------------------
       "integer(i4)" i4 = SELECTED_INT_KIND(8)
       DIGITS    =   31
       HUGE      =   2147483647
       KIND      =   4
      ----------------------------------------------
       "integer(i8)" i8 = SELECTED_INT_KIND(13)
       DIGITS    =   63
       HUGE      =   9223372036854775807
       KIND      =   8
      ----------------------------------------------
       "real"
       DIGITS    =   24
       EPSILON   =   1.192093E-07
       HUGE      =   3.402823E+38
       KIND      =   4
       PRECISION =   6
      ----------------------------------------------
       "real(r4)" r4 = SELECTED_REAL_KIND(6,30)
       DIGITS    =   24
       EPSILON   =   1.192093E-07
       HUGE      =   3.402823E+38
       KIND      =   4
       PRECISION =   6
      ----------------------------------------------
       "real(r8)" r8 = SELECTED_REAL_KIND(13)
       DIGITS    =   53
       EPSILON   =   2.220446049250313E-016
       HUGE      =   1.797693134862315E+308
       KIND      =   8
       PRECISION =   15
      ----------------------------------------------
       "double precision"
       DIGITS    =   53
       EPSILON   =   2.220446049250313E-016
       HUGE      =   1.797693134862315E+308
       KIND      =   8
       PRECISION =   15

PowerPC chipset : Absoft Pro Fortran 9.0 : "-dp"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. container:: unix

   ::

      [~/DART/utilities] % f90 -dp PrecisionCheck.f90
      [~/DART/utilities] % ./a.out
       
       This explores the use of the intrinsic SELECTED_[REAL,INT]_KIND() functions
       and the interplay with the compiler options. You are encouraged to use the
       "autopromotion" flags on your compiler and compare the results.
       
      ----------------------------------------------
       "integer"
       DIGITS    =   31
       HUGE      =   2147483647
       KIND      =   4
      ----------------------------------------------
       "integer(i4)" i4 = SELECTED_INT_KIND(8)
       DIGITS    =   31
       HUGE      =   2147483647
       KIND      =   4
      ----------------------------------------------
       "integer(i8)" i8 = SELECTED_INT_KIND(13)
       DIGITS    =   63
       HUGE      =   9223372036854775807
       KIND      =   8
      ----------------------------------------------
       "real"
       DIGITS    =   24
       EPSILON   =   1.192093E-07
       HUGE      =   3.402823E+38
       KIND      =   4
       PRECISION =   6
      ----------------------------------------------
       "real(r4)" r4 = SELECTED_REAL_KIND(6,30)
       DIGITS    =   24
       EPSILON   =   1.192093E-07
       HUGE      =   3.402823E+38
       KIND      =   4
       PRECISION =   6
      ----------------------------------------------
       "real(r8)" r8 = SELECTED_REAL_KIND(13)
       DIGITS    =   53
       EPSILON   =   2.220446049250313E-016
       HUGE      =   1.797693134862315E+308
       KIND      =   8
       PRECISION =   15
      ----------------------------------------------
       "double precision"
       DIGITS    =   24
       EPSILON   =   1.192093E-07
       HUGE      =   3.402823E+38
       KIND      =   4
       PRECISION =   6

PowerPC chipset : Absoft Pro Fortran 9.0 : "-n113"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. container:: unix

   ::

      [~/DART/utilities] % f90 -N113 PrecisionCheck.f90
      [~/DART/utilities] % ./a.out
       
       This explores the use of the intrinsic SELECTED_[REAL,INT]_KIND() functions
       and the interplay with the compiler options. You are encouraged to use the
       "autopromotion" flags on your compiler and compare the results.
       
      ----------------------------------------------
       "integer"
       DIGITS    =   31
       HUGE      =   2147483647
       KIND      =   4
      ----------------------------------------------
       "integer(i4)" i4 = SELECTED_INT_KIND(8)
       DIGITS    =   31
       HUGE      =   2147483647
       KIND      =   4
      ----------------------------------------------
       "integer(i8)" i8 = SELECTED_INT_KIND(13)
       DIGITS    =   63
       HUGE      =   9223372036854775807
       KIND      =   8
      ----------------------------------------------
       "real"
       DIGITS    =   53
       EPSILON   =   2.220446049250313E-016
       HUGE      =   1.797693134862315E+308
       KIND      =   8
       PRECISION =   15
      ----------------------------------------------
       "real(r4)" r4 = SELECTED_REAL_KIND(6,30)
       DIGITS    =   24
       EPSILON   =   1.192093E-07
       HUGE      =   3.402823E+38
       KIND      =   4
       PRECISION =   6
      ----------------------------------------------
       "real(r8)" r8 = SELECTED_REAL_KIND(13)
       DIGITS    =   53
       EPSILON   =   2.220446049250313E-016
       HUGE      =   1.797693134862315E+308
       KIND      =   8
       PRECISION =   15
      ----------------------------------------------
       "double precision"
       DIGITS    =   53
       EPSILON   =   2.220446049250313E-016
       HUGE      =   1.797693134862315E+308
       KIND      =   8
       PRECISION =   15
