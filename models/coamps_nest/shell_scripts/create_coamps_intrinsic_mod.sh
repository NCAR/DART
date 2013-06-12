#!/bin/bash
#
# This code may (or may not) be part of the COAMPS distribution,
# So it is not protected by the DART copyright agreement.
#
# DART $Id$
#
########################################################################
#
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
DIAG=diagnostic_cdf/src

INTRINSIC_SUBS="ij2ll 
                ll2ij 
                s2pint 
				sfcp
                z2zint 
                utom 
                vtom
                uvg2uv
                pstd
                tstd
                deriv_x
                deriv_y
                vor
                filt9"

MODULE_FILE='coamps_intrinsic_mod.f90'

# Process Fortran routines
for routine in ${INTRINSIC_SUBS}
do
  f_file=${routine}.f
  f90_file=${routine}.f90

  if [ -e ${COAMPS_UTIL}/${DIAG}/${f90_file} ]; then 
    cp -f ${COAMPS_UTIL}/${DIAG}/${f90_file}  ./${f90_file}
  else
    cp -f ${COAMPS_UTIL}/${DIAG}/${f_file}  ./${f90_file}
  fi

  # Switch line continuations from a single character in column 6
  # to an end-of-line ampersand
  sed -i -e "N;s/\n\(\s*\)[0-9]\s/\&\n\1  /;P;D;" ${f90_file}

  # Switch to the DART-style real
  sed -i -e "s/\(\s*\)real\(\s*\)/\1real(kind=r8)\2/g" ${f90_file}

  # Some compilers require END SUBROUTINE instead of just END
  sed -i -e "s/\(\s*\)end\s*$/\1end subroutine ${routine}/g" ${f90_file}

  # Remove whitespace before 'subroutine' and 'end subroutine'
  sed -i -e "s/^\s*\(subroutine.*$\)/\1/g" ${f90_file}
  sed -i -e "s/^\s*\(end subroutine.*$\)/\1/g" ${f90_file}

  # Switch comment style
  sed -i -e "s/^c/!/gi" ${f90_file}

  # Force double-precision logs
  sed -i -e "s/alog/dlog/g" ${f90_file}
  sed -i -e "s/amod/dmod/g" ${f90_file}

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
cat <<END_SPACE >> ${MODULE_FILE}

END_SPACE
  rm ${intrinsic_filename}
done

# Write routines that define standard exner function and theta
# In coamps/diagnostic these are done all together in astate.f
# but we need individual functions.
cat <<END_FOOTER >> ${MODULE_FILE}
  ! define_mean_exner
  ! ---------------------
  ! Defines the mean exner based on the standard 
  ! atmosphere on sigma surfaces in a 3D domain.
  !  PARAMETERS
  !   OUT  exner   Mean exner function 
  !   IN   mn      Horizontal dimension
  !   IN   kk      Vertical levels
  !   IN   ztop    Domain top
  !   IN   zsfc    Topography height
  !   IN   sigma   Sigma levels
  subroutine define_mean_exner(exner, mn, kk, ztop, zsfc, sigma)
    real(kind=r8), dimension(:,:), intent(out) :: exner
    integer,                       intent(in)  :: mn
    integer,                       intent(in)  :: kk
    real(kind=r8),                 intent(in)  :: ztop
    real(kind=r8), dimension(:),   intent(in)  :: zsfc
    real(kind=r8), dimension(:),   intent(in)  :: sigma

    real(kind=r8) :: zout(mn)
    real(kind=r8) :: pres(mn)
    integer       :: k

    do k=1,kk
      zout(:)=( sigma(k)*(ztop-zsfc(:))/ztop + zsfc(:) )*g
      call pstd(mn, zout, pres)
      exner(:,k)=( (pres(:)*100.d0 )/p00)**rocp
    enddo
  end subroutine define_mean_exner

  ! define_mean_theta
  ! ---------------------
  ! Defines the mean theta based on the standard 
  ! atmosphere on sigma surfaces in a 3D domain.
  !  PARAMETERS
  !   OUT  theta   Mean theta function 
  !   IN   mn      Horizontal dimension
  !   IN   kk      Vertical levels
  !   IN   ztop    Domain top
  !   IN   zsfc    Topography height
  !   IN   sigma   Sigma levels
  subroutine define_mean_theta(theta, mn, kk, ztop, zsfc, sigma)
    real(kind=r8), dimension(:,:), intent(out) :: theta
    integer,                       intent(in)  :: mn
    integer,                       intent(in)  :: kk
    real(kind=r8),                 intent(in)  :: ztop
    real(kind=r8), dimension(:),   intent(in)  :: zsfc
    real(kind=r8), dimension(:),   intent(in)  :: sigma

    real(kind=r8) :: zout(mn)
    real(kind=r8) :: pres(mn)
    integer       :: k

    do k=1,kk
      zout(:)=( sigma(k)*(ztop-zsfc(:))/ztop + zsfc(:) )*g
      call pstd(mn, zout, pres)
      call tstd(mn, pres, zout)
      theta(:,k)=zout(:)*(p00/(pres(:)*100.d0))**rocp
    enddo
  end subroutine define_mean_theta

  ! compute_altimeter
  ! ---------------------
  ! Computes the altimeter of a pressure value at
  ! a gemotric height
  !  PARAMETERS
  !   IN   pres    pressure value  (Pa)
  !   IN   zhgt    pressure height (m)
  function compute_altimeter(pres, zhgt) result(altimeter)
    real(kind=r8), intent(in) :: pres
    real(kind=r8), intent(in) :: zhgt
    real(kind=r8)             :: altimeter

    real(kind=r8), parameter :: k1 = 0.190284_r8
    real(kind=r8), parameter :: k2 = 8.4228807E-5_r8

    real(kind=r8), parameter  :: hPa_to_Pa= 100.0_r8

    altimeter = hPa_to_Pa * ((pres/hPa_to_Pa - 0.3_r8) ** k1 + k2 * zhgt) ** (1.0_r8 / k1)

    return
  end function compute_altimeter
END_FOOTER

# Write the module footer
cat <<END_FOOTER >> ${MODULE_FILE}

end module coamps_intrinsic_mod
END_FOOTER

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

