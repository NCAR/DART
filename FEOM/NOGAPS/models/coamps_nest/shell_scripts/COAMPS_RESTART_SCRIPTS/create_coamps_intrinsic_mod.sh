#!/bin/bash
#
# This code may (or may not) be part of the COAMPS distribution,
# So it is not protected by the DART copyright agreement.
#
# DART $Id$
#
########################################################################
#
# SCRIPT:     create_coamps_intrinsic_mod.sh
# AUTHOR:     T. R. Whitcomb
#             Naval Research Laboratory
#
# Create a "coamps_intrinsic_mod" module that patches and then
# assembles all of the subroutines that I borrowed directly from
# COAMPS utility package

# This program expects the path to the COAMPS utility package as
# a command line argument

COAMPS_UTIL=$1


# directories for COAMPS utility programs
DIAG=diagnostic

INTRINSIC_SUBS="ij2ll 
                ll2ij 
                s2pint 
                seaprs 
                z2zint 
                utom 
                vtom
                filt9"

MODULE_FILE='coamps_intrinsic_mod.f90'

# Process Fortran routines
for routine in ${INTRINSIC_SUBS}
do
  f_file=${routine}.f
  f90_file=${routine}.f90
  cp -f ${COAMPS_UTIL}/${DIAG}/${f_file}  ./${f90_file}

  # Switch line continuations from a single character in column 6
  # to an end-of-line ampersand
  sed -i -e "N;s/\n\(\s*\)[0-9]\s/\&\n\1  /;P;D;" ${f90_file}

  # Switch to the DART-style real
  sed -i -e "s/\(\s*\)real\(\s*\)/\1real(kind=r8)\2/g" ${f90_file}

  # Some compilers require END SUBROUTINE instead of just END
  sed -i -e "s/\(\s*\)end\s*$/\1end subroutine ${routine}/g" ${f90_file}

  # Switch comment style
  sed -i -e "s/^c/!/g" ${f90_file}

  # Force double-precision logs
  sed -i -e "s/alog/dlog/g" ${f90_file}

  # We'll put the contents of constant.h in our file
  sed -i -e "/include/d" ${f90_file}
done

# Process header file
header=constant.h
orig_header=${COAMPS_UTIL}/${DIAG}/${header}
sed -n -e '/=/ s/^\(\s*\)\(.*\)$/\1real(kind=r8), parameter :: \2/gp' < ${orig_header} > ${header} 
sed -i -e 's/\(pi\s*=\).*/\1 3.141592741012573/' ${header}

# Write module header
cat <<END_HEADER > ${MODULE_FILE}
!------------------------------
! MODULE:       coamps_intrinsic_mod
! AUTHOR:       T. R. Whitcomb
!               Naval Research Laboratory
! MODIFIED:     07 July 2008
! DART VERSION: Jamaica
!
! Collect the routines taken directly from the COAMPS utility
! package that have been translated to Fortran 90 and patched
! to use real(kind=r8) instead of whatever was there before. 
! This module is automatically generated.
!------------------------------
module coamps_intrinsic_mod
  use types_mod, only : r8

  implicit none

  !------------------------------
  ! BEGIN PUBLIC INTERFACE
  !------------------------------

  ! All routines are accessible
  public
  
  !------------------------------
  ! END PUBLIC INTERFACE
  !------------------------------

  !------------------------------
  ! BEGIN EXTERNAL INTERFACES
  !------------------------------
  ! [none]
  !------------------------------
  ! END EXTERNAL INTERFACES
  !------------------------------
END_HEADER

# Write constants information
cat <<END_HEADER >> ${MODULE_FILE}
  !------------------------------
  ! BEGIN TYPES AND CONSTANTS 
  !------------------------------
END_HEADER

cat ${header} >> ${MODULE_FILE}
rm ${header}

cat <<END_FOOTER >> ${MODULE_FILE}
  !------------------------------
  ! END TYPES AND CONSTANTS 
  !------------------------------

  !------------------------------
  ! BEGIN MODULE VARIABLES
  !------------------------------
  ! [none]
  !------------------------------
  ! END MODULE VARIABLES
  !------------------------------

contains
END_FOOTER

# Now insert all the routines
for intrinsic_routine in ${INTRINSIC_SUBS}
do
  intrinsic_filename=${intrinsic_routine}.f90
  cat ${intrinsic_filename} >> ${MODULE_FILE}
  rm ${intrinsic_filename}
done

# Write the module footer
cat <<END_FOOTER >> ${MODULE_FILE}
end module coamps_intrinsic_mod
END_FOOTER


# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

