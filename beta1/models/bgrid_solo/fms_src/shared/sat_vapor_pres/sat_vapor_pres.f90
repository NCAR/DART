!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module sat_vapor_pres_mod

!-----------------------------------------------------------------------
!
!                 saturation vapor pressure lookup
!
!      routines for computing the saturation vapor pressure (es),
!      the derivation of es with respect to temperature, and
!      initialization of the look-up table.
!
!-----------------------------------------------------------------------
!
!                               usage
!                               -----
!
!                       call lookup_es  (temp,es)
!
!                       call lookup_des (temp,des)
!
!    arguments
!    ---------
!      temp    intent in       temperature in degrees kelvin
!      es      intent out      saturation vapor pressure in Pascals
!      des     intent out      derivative of saturation vapor pressure
!                              with respect to temperature (Pascals/degree)
!
!-----------------------------------------------------------------------

! <CONTACT EMAIL="bw@gfdl.noaa.gov">
!   Bruce Wyman
! </CONTACT>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

! <OVERVIEW>
!   Routines for determining the saturation vapor pressure (<TT>ES</TT>) 
!   and the derivative of <TT>ES</TT> with respect to temperature.
! </OVERVIEW>

! <DESCRIPTION>
!   This module contains routines for determining the saturation vapor
!   pressure (<TT>ES</TT>) from lookup tables constructed using equations given
!   in the Smithsonian tables.  The <TT>ES</TT> lookup tables are valid between
!   -160C and +100C (approx 113K to 373K).

!   The values of <TT>ES</TT> are computed over ice from -160C to -20C,
!   over water from 0C to 100C, and a blended value (over water and ice)
!   from -20C to 0C.

!   This version was written for non-vector machines.
!   See the <LINK SRC="#NOTES">notes</LINK> section for details on vectorization.

! </DESCRIPTION>

! <PUBLIC>
!   Description summarizing public interface.
! </PUBLIC>

 use  constants_mod, only:  TFREEZE
 use        fms_mod, only:  write_version_number,   &
                            error_mesg, FATAL

implicit none
private

 public :: lookup_es, lookup_des, sat_vapor_pres_init
 public :: compute_es
 public :: escomp, descomp ! for backward compatibility
                           ! use lookup_es, lookup_des instead

!-----------------------------------------------------------------------
 
! <INTERFACE NAME="lookup_es">

!   <OVERVIEW>
!     For the given temperatures, returns the saturation vapor pressures.
!   </OVERVIEW>
!   <DESCRIPTION>
!     For the given temperatures these routines return the
!     saturation vapor pressure (esat). The return values are derived from
!     lookup tables (see notes below).
!   </DESCRIPTION>
!   <TEMPLATE>
!     call lookup_es( temp, esat )
!   </TEMPLATE>
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(scalar),(:),(:,:),(:,:,:)">
!     Temperature in degrees Kelvin.
!   </IN>
!   <OUT NAME="esat" UNITS="pascal" TYPE="real" DIM="(scalar),(:),(:,:),(:,:,:)">
!     Saturation vapor pressure in pascals.
!             May be a scalar, 1d, 2d, or 3d array.
!             Must have the same order and size as temp.
!   </OUT>
!   <ERROR MSG="table overflow, nbad=##" STATUS="FATAL">
!     Temperature(s) provided to the saturation vapor pressure lookup
!          are outside the valid range of the lookup table (-160 to 100 deg C).
!          This may be due to a numerical instability in the model.
!          Information should have been printed to standard output to help
!          determine where the instability may have occurred.
!          If the lookup table needs a larger temperature range,
!          then parameters in the module header must be modified.
!   </ERROR> *

 interface lookup_es
   module procedure lookup_es_0d, lookup_es_1d, lookup_es_2d, lookup_es_3d
 end interface
! for backward compatibility (to be removed soon)
 interface escomp
   module procedure lookup_es_0d, lookup_es_1d, lookup_es_2d, lookup_es_3d
 end interface
! </INTERFACE>
!-----------------------------------------------------------------------
! <INTERFACE NAME="lookup_des">

!   <OVERVIEW>
!     For the given temperatures, returns the derivative of saturation vapor pressure
!     with respect to temperature.
!   </OVERVIEW>
!   <DESCRIPTION>
!     For the given temperatures these routines return the derivative of esat w.r.t.
!     temperature (desat). The return values are derived from
!     lookup tables (see notes below).
!   </DESCRIPTION>
!   <TEMPLATE>
!     call lookup_des( temp, desat )
!   </TEMPLATE>
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(scalar),(:),(:,:),(:,:,:)">
!     Temperature in degrees Kelvin.
!   </IN>
!   <OUT NAME="desat" UNITS="pascal" TYPE="real" DIM="(scalar),(:),(:,:),(:,:,:)">
!     Derivative of saturation vapor pressure w.r.t. temperature
!                 in pascals/degree. May be a scalar, 1d, 2d, or 3d array.
!                 Must have the same order and size as temp.
!   </OUT>

 interface lookup_des
   module procedure lookup_des_0d, lookup_des_1d, lookup_des_2d, lookup_des_3d
 end interface
! </INTERFACE>
! for backward compatibility (to be removed soon)
 interface descomp
   module procedure lookup_des_0d, lookup_des_1d, lookup_des_2d, lookup_des_3d
 end interface
!-----------------------------------------------------------------------
! <INTERFACE NAME="compute_es">

!   <OVERVIEW>
!     For the given temperatures, computes the saturation vapor pressures. 
!   </OVERVIEW>
!   <DESCRIPTION>
!     Computes saturation vapor pressure for the given temperature using
!     the equations given in the Smithsonian Meteorological Tables.
!     Between -20C and 0C a blended value over ice and water is returned.
!   </DESCRIPTION>
!   <TEMPLATE>
!     es = compute_es ( temp )
!   </TEMPLATE>
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(scalar),(:),(:,:),(:,:,:)">
!     Temperature in degrees Kelvin.
!   </IN>
!   <OUT NAME="es" UNITS="pascal" TYPE="real" DIM="(scalar),(:),(:,:),(:,:,:)">
!     Saturation vapor pressure in pascals.
!             May be a scalar, 1d, 2d, or 3d array.
!             Must have the same order and size as temp.
!   </OUT>

 interface compute_es
   module procedure compute_es_0d, compute_es_1d, compute_es_2d, compute_es_3d
 end interface
! </INTERFACE>
!-----------------------------------------------------------------------
 interface temp_check
   module procedure temp_check_0d, temp_check_1d, temp_check_2d, temp_check_3d
 end interface
!-----------------------------------------------------------------------
!  cvs version and tag name

character(len=128) :: version = '$Revision$'
character(len=128) :: tagname = '$Id$'

!-----------------------------------------------------------------------
!  parameters for table size and resolution

 integer, parameter :: tcmin = -160  ! minimum temperature (degC) in lookup table
 integer, parameter :: tcmax =  100  ! maximum temperature (degC) in lookup table
 integer, parameter :: esres =  10   ! table resolution (increments per degree)
 integer, parameter :: nsize = (tcmax-tcmin)*esres+1    !  lookup table size
 integer, parameter :: nlim  = nsize-1

 real    :: tmin, tmax          !  lookup table limits in degK
 real    :: dtres, dtinv, teps

 real ::   TABLE(nsize)    !  sat vapor pres (es)
 real ::  DTABLE(nsize)    !  first derivative of es
 real :: D2TABLE(nsize)    ! second derivative of es

 logical :: module_is_initialized = .FALSE.

!-----------------------------------------------------------------------

contains

!#######################################################################

! <SUBROUTINE NAME="lookup_es_0d" INTERFACE="lookup_es">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(scalar)"></IN>
!   <OUT NAME="esat" UNITS="pascal" TYPE="real" DIM="(scalar)"></OUT>
! </SUBROUTINE>
 subroutine lookup_es_0d ( temp, esat )

 real, intent(in)  :: temp
 real, intent(out) :: esat

 real    :: tmp, del
 integer :: ind
!-----------------------------------------------

   if (.not.module_is_initialized) call sat_vapor_pres_init

   tmp = temp-tmin
   ind = int(dtinv*(tmp+teps))
   del = tmp-dtres*real(ind)
   esat = TABLE(ind+1) + del*(DTABLE(ind+1) + del*D2TABLE(ind+1))
!!!esat = TABLE(ind+1) + del*DTABLE(ind+1)

     if (ind < 0 .or. ind > nlim) call temp_check ( 1, temp )

!-----------------------------------------------

 end subroutine lookup_es_0d

!#######################################################################

! <SUBROUTINE NAME="lookup_es_1d" INTERFACE="lookup_es">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:)"></IN>
!   <OUT NAME="esat" UNITS="pascal" TYPE="real" DIM="(:)"></OUT>
! </SUBROUTINE>
 subroutine lookup_es_1d ( temp, esat )

 real, intent(in)  :: temp(:)
 real, intent(out) :: esat(:)

 real    :: tmp, del
 integer :: ind, i, n
!-----------------------------------------------

   if (.not.module_is_initialized) call sat_vapor_pres_init

   n = 0
   do i = 1, size(temp,1)
     tmp = temp(i)-tmin
     ind = int(dtinv*(tmp+teps))
     del = tmp-dtres*real(ind)
     esat(i) = TABLE(ind+1) + del*(DTABLE(ind+1) + del*D2TABLE(ind+1))
!!!!!esat(i) = TABLE(ind+1) + del*DTABLE(ind+1)
     if (ind < 0 .or. ind > nlim) n = n+1
   enddo

   if ( n > 0 ) call temp_check ( n, temp )

!-----------------------------------------------

 end subroutine lookup_es_1d

!#######################################################################

! <SUBROUTINE NAME="lookup_es_2d" INTERFACE="lookup_es">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:,:)"></IN>
!   <OUT NAME="esat" UNITS="pascal" TYPE="real" DIM="(:,:)"></OUT>
! </SUBROUTINE>
 subroutine lookup_es_2d ( temp, esat )

 real, intent(in)  :: temp(:,:)
 real, intent(out) :: esat(:,:)

 real    :: tmp, del
 integer :: ind, i, j, n
!-----------------------------------------------

   if (.not.module_is_initialized) call sat_vapor_pres_init

   n = 0
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j)-tmin
     ind = int(dtinv*(tmp+teps))
     del = tmp-dtres*real(ind)
     esat(i,j) = TABLE(ind+1) + del*(DTABLE(ind+1) + del*D2TABLE(ind+1))
!!!!!esat(i,j) = TABLE(ind+1) + del*DTABLE(ind+1)
     if (ind < 0 .or. ind > nlim) n = n+1
   enddo
   enddo

   if ( n > 0 ) call temp_check_2d ( n, temp )

!-----------------------------------------------

 end subroutine lookup_es_2d

!#######################################################################

! <SUBROUTINE NAME="lookup_es_3d" INTERFACE="lookup_es">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:,:,:)"></IN>
!   <OUT NAME="esat" UNITS="pascal" TYPE="real" DIM="(:,:,:)"></OUT>
! </SUBROUTINE>
 subroutine lookup_es_3d ( temp, esat )

 real, intent(in)  :: temp(:,:,:)
 real, intent(out) :: esat(:,:,:)

 real    :: tmp, del
 integer :: ind, i, j, k, n
!-----------------------------------------------

   if (.not.module_is_initialized) call sat_vapor_pres_init

   n = 0
   do k = 1, size(temp,3)
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j,k)-tmin
     ind = int(dtinv*(tmp+teps))
     del = tmp-dtres*real(ind)
     esat(i,j,k) = TABLE(ind+1) + del*(DTABLE(ind+1) + del*D2TABLE(ind+1))
!!!!!esat(i,j,k) = TABLE(ind+1) + del*DTABLE(ind+1)
     if (ind < 0 .or. ind > nlim) n = n+1
   enddo
   enddo
   enddo

   if ( n > 0 ) call temp_check ( n, temp )

!-----------------------------------------------

 end subroutine lookup_es_3d

!#######################################################################
!  routines for computing derivative of es
!#######################################################################

! <SUBROUTINE NAME="lookup_des_0d" INTERFACE="lookup_des">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(scalar)"></IN>
!   <OUT NAME="desat" UNITS="pascal" TYPE="real" DIM="(scalar)"></OUT>
! </SUBROUTINE>
 subroutine lookup_des_0d ( temp, desat )

 real, intent(in)  :: temp
 real, intent(out) :: desat

 real    :: tmp, del
 integer :: ind
!-----------------------------------------------

   if (.not.module_is_initialized) call sat_vapor_pres_init

   tmp = temp-tmin
   ind = int(dtinv*(tmp+teps))
   del = tmp-dtres*real(ind)
   desat = DTABLE(ind+1) + 2.*del*D2TABLE(ind+1)

   if (ind < 0 .or. ind > nlim) call temp_check ( 1, temp )

!-----------------------------------------------

 end subroutine lookup_des_0d

!#######################################################################

! <SUBROUTINE NAME="lookup_des_1d" INTERFACE="lookup_des">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:)"></IN>
!   <OUT NAME="desat" UNITS="pascal" TYPE="real" DIM="(:)"></OUT>
! </SUBROUTINE>
 subroutine lookup_des_1d ( temp, desat )

 real, intent(in)  :: temp (:)
 real, intent(out) :: desat(:)

 real    :: tmp, del
 integer :: ind, i, n
!-----------------------------------------------

   if (.not.module_is_initialized) call sat_vapor_pres_init

   n = 0
   do i = 1, size(temp,1)
     tmp = temp(i)-tmin
     ind = int(dtinv*(tmp+teps))
     del = tmp-dtres*real(ind)
     desat(i) = DTABLE(ind+1) + 2.*del*D2TABLE(ind+1)
     if (ind < 0 .or. ind > nlim) n = n+1
   enddo

   if ( n > 0 ) call temp_check ( n, temp )

!-----------------------------------------------

 end subroutine lookup_des_1d

!#######################################################################

! <SUBROUTINE NAME="lookup_des_2d" INTERFACE="lookup_des">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:,:)"></IN>
!   <OUT NAME="desat" UNITS="pascal" TYPE="real" DIM="(:,:)"></OUT>
! </SUBROUTINE>
 subroutine lookup_des_2d ( temp, desat )

 real, intent(in)  :: temp (:,:)
 real, intent(out) :: desat(:,:)

 real    :: tmp, del
 integer :: ind, i, j, n
!-----------------------------------------------
   
   if (.not.module_is_initialized) call sat_vapor_pres_init
   
   n = 0
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j)-tmin
     ind = int(dtinv*(tmp+teps))
     del = tmp-dtres*real(ind)
     desat(i,j) = DTABLE(ind+1) + 2.*del*D2TABLE(ind+1)
     if (ind < 0 .or. ind > nlim) n = n+1
   enddo
   enddo

   if ( n > 0 ) call temp_check ( n, temp )

!-----------------------------------------------

 end subroutine lookup_des_2d

!#######################################################################
! <SUBROUTINE NAME="lookup_des_3d" INTERFACE="lookup_des">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:,:,:)"></IN>
!   <OUT NAME="desat" UNITS="pascal" TYPE="real" DIM="(:,:,:)"></OUT>
! </SUBROUTINE>
 subroutine lookup_des_3d ( temp, desat )

 real, intent(in)  :: temp (:,:,:)
 real, intent(out) :: desat(:,:,:)

 real    :: tmp, del
 integer :: ind, i, j, k, n
!-----------------------------------------------

   if (.not.module_is_initialized) call sat_vapor_pres_init

   n = 0
   do k = 1, size(temp,3)
   do j = 1, size(temp,2)
   do i = 1, size(temp,1)
     tmp = temp(i,j,k)-tmin
     ind = int(dtinv*(tmp+teps))
     del = tmp-dtres*real(ind)
     desat(i,j,k) = DTABLE(ind+1) + 2.*del*D2TABLE(ind+1)
     if (ind < 0 .or. ind > nlim) n = n+1
   enddo
   enddo
   enddo

   if ( n > 0 ) call temp_check ( n, temp )
   
!-----------------------------------------------

 end subroutine lookup_des_3d

!#######################################################################
!#######################################################################

! <SUBROUTINE NAME="sat_vapor_pres_init">

!   <OVERVIEW>
!     Initializes the lookup tables for saturation vapor pressure. 
!   </OVERVIEW>
!   <DESCRIPTION>
!     Initializes the lookup tables for saturation vapor pressure.
!     This routine will be called automatically the first time
!     <B>lookup_es</B> or <B>lookup_des</B> is called,
!     the user does not need to call this routine.
!     There are no arguments.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call sat_vapor_pres_init
!   </TEMPLATE>

! </SUBROUTINE>
 subroutine sat_vapor_pres_init

!  =================================================================
!  +                                                               +
!  +             construction of the es table                      +
!  +                                                               +
!  + this table is constructed from es equations from the          +
!  + smithsonian tables.  the es input is computed from values     +
!  + (in one-tenth of a degree increments) of es over ice          +
!  + from -153c to 0c and values of es over water from 0c to 102c. +
!  + output table contains these data interleaved with their       +
!  + derivatives with respect to temperature except between -20c   +
!  + and 0c where blended (over water and over ice) es values and  +
!  + derivatives are calculated.                                   +
!  +   note: all es computation is done in pascals                 +
!  =================================================================

 real    :: tem(3), es(3), hdtinv
 integer :: i, n

! increment used to generate derivative table
  real, parameter :: tinrc = .01           
  real, parameter :: tfact = 1./(2.*tinrc)

! return silently if this routine has already been called
      if (module_is_initialized) return

! write version number to log file
      call write_version_number (version, tagname)

! global variables
      tmin = real(tcmin)+TFREEZE   ! minimum valid temp in table
      tmax = real(tcmax)+TFREEZE   ! maximum valid temp in table
      dtinv = real(esres)
      dtres = 1./dtinv
      teps = 1./real(2*esres)
! local variables
      hdtinv = dtinv*0.5

! compute es tables from tcmin to tcmax
! estimate es derivative with small +/- difference

      do i = 1, nsize
         tem(1) = tmin + dtres*real(i-1)
         tem(2) = tem(1)-tinrc
         tem(3) = tem(1)+tinrc
         es = compute_es (tem)
          TABLE(i) = es(1)
         DTABLE(i) = (es(3)-es(2))*tfact
      enddo

! compute one-half second derivative using centered differences
! differencing des values in the table

      do i = 2, nsize-1
         D2TABLE(i) = 0.25*dtinv*(DTABLE(i+1)-DTABLE(i-1))
      enddo
    ! one-sided derivatives at boundaries
         D2TABLE(1)     = 0.50*dtinv*(DTABLE(2)    -DTABLE(1))
         D2TABLE(nsize) = 0.50*dtinv*(DTABLE(nsize)-DTABLE(nsize-1))

    module_is_initialized = .true.

 end subroutine sat_vapor_pres_init

!#######################################################################
!#######################################################################
!-------------------------------------------------------------------
!                Computation of the es values
!
!   Saturation vapor pressure (es) values are computed from
!   equations in the Smithsonian meteorological tables page 350.
!   For temperatures < 0C, sat vapor pres is computed over ice.
!   For temperatures > -20C, sat vapor pres is computed over water.
!   Between -20C and 0C the returned value is blended (over water
!   and over ice).  All sat vapor pres values are returned in pascals.
!
!   Reference:  Smithsonian meteorological tables, page 350.
!-------------------------------------------------------------------

! <FUNCTION NAME="compute_es_1d" INTERFACE="compute_es">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:)"></IN>
!   <OUT NAME="es" UNITS="pascal" TYPE="real" DIM="(:)"></OUT>
! </FUNCTION>
 function compute_es_1d (tem) result (es)
 real, intent(in) :: tem(:)
 real :: es(size(tem))

 real, parameter :: TBASW = TFREEZE+100.
 real, parameter :: TBASI = TFREEZE
 real, parameter :: ESBASW = 101324.60
 real, parameter :: ESBASI =    610.71

 real    :: x, esice, esh2o
 integer :: i

   do i = 1, size(tem)

!  compute es over ice 

     if (tem(i) < TBASI) then
         x = -9.09718*(TBASI/tem(i)-1.0) - 3.56654*log10(TBASI/tem(i)) &
             +0.876793*(1.0-tem(i)/TBASI) + log10(ESBASI)
         esice =10.**(x)
     else
         esice = 0.
     endif

!  compute es over water greater than -20 c.
!  values over 100 c may not be valid
!  see smithsonian meteorological tables page 350.

     if (tem(i) > -20.+TBASI) then
         x = -7.90298*(TBASW/tem(i)-1) + 5.02808*log10(TBASW/tem(i)) &
             -1.3816e-07*(10**((1-tem(i)/TBASW)*11.344)-1)        &
             +8.1328e-03*(10**((TBASW/tem(i)-1)*(-3.49149))-1)    &
             +log10(ESBASW)
         esh2o = 10.**(x)
     else
         esh2o = 0.
     endif

!  derive blended es over ice and supercooled water between -20c and 0c

     if (tem(i) <= -20.+TBASI) then
         es(i) = esice
     else if (tem(i) >= TBASI) then
         es(i) = esh2o
     else
         es(i) = 0.05*((TBASI-tem(i))*esice + (tem(i)-TBASI+20.)*esh2o)
     endif

   enddo
   
 end function compute_es_1d

!--------------------------------------------------------

! <FUNCTION NAME="compute_es_0d" INTERFACE="compute_es">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(scalar)"></IN>
!   <OUT NAME="es" UNITS="pascal" TYPE="real" DIM="(scalar)"></OUT>
! </FUNCTION>
 function compute_es_0d (tem) result (es)
 real, intent(in) :: tem
 real :: es
 real, dimension(1) :: tem1, es1

   tem1(1) = tem
   es1 = compute_es_1d (tem1)
   es = es1(1)

 end function compute_es_0d

!--------------------------

! <FUNCTION NAME="compute_es_2d" INTERFACE="compute_es">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:,:)"></IN>
!   <OUT NAME="es" UNITS="pascal" TYPE="real" DIM="(:,:)"></OUT>
! </FUNCTION>
 function compute_es_2d (tem) result (es)
 real, intent(in) :: tem(:,:)
 real, dimension(size(tem,1),size(tem,2)) :: es
 integer :: j

    do j = 1, size(tem,2)
      es(:,j) = compute_es_1d (tem(:,j))
    enddo

 end function compute_es_2d

!--------------------------
! <FUNCTION NAME="compute_es_3d" INTERFACE="compute_es">
!   <IN NAME="temp" UNIT="degrees Kelvin" TYPE="real" DIM="(:,:,:)"></IN>
!   <OUT NAME="es" UNITS="pascal" TYPE="real" DIM="(:,:,:)"></OUT>
! </FUNCTION>
 function compute_es_3d (tem) result (es)
 real, intent(in) :: tem(:,:,:)
 real, dimension(size(tem,1),size(tem,2),size(tem,3)) :: es
 integer :: j, k

    do k = 1, size(tem,3)
    do j = 1, size(tem,2)
      es(:,j,k) = compute_es_1d (tem(:,j,k))
    enddo
    enddo

 end function compute_es_3d

!#######################################################################
!#######################################################################

 subroutine error_handler ( n )
 integer, intent(in) :: n
 character(len=28) :: mesg

   write (mesg,'(a21,i7)') 'table overflow, nbad=', n

   call error_mesg ('sat_vapor_pres_mod', mesg, FATAL)

!  print *, 'ERROR: ' // mesg
!  stop 111

 end subroutine error_handler

!#######################################################################

 function check_1d ( temp ) result ( nbad )
 real   , intent(in)  :: temp(:)
 integer :: nbad, ind, i
 real    :: tmp

   nbad = 0
   do i = 1, size(temp,1)
     ind = int(dtinv*(temp(i)-tmin+teps))
     if (ind < 0 .or. ind > nlim) nbad = nbad+1
   enddo

 end function check_1d

!------------------------------------------------

 function check_2d ( temp ) result ( nbad )
 real   , intent(in)  :: temp(:,:)
 integer :: nbad
 integer :: j, ind
 real    :: tmp
    nbad = 0
    do j = 1, size(temp,2)
      nbad = nbad + check_1d ( temp(:,j) )
    enddo
 end function check_2d

!#######################################################################

 subroutine temp_check_0d ( nbad, temp )
 integer, intent(in) :: nbad
 real   , intent(in) :: temp

   call error_handler (nbad)

 end subroutine temp_check_0d

!--------------------------------------------------------------

 subroutine temp_check_1d ( nbad, temp )
 integer, intent(in) :: nbad
 real   , intent(in) :: temp(:)
 integer :: i

   print *, 'Bad temperatures (dimension 1): ', (check_1d(temp(i:i)),i=1,size(temp,1))
   call error_handler (nbad)

 end subroutine temp_check_1d

!--------------------------------------------------------------

 subroutine temp_check_2d ( nbad, temp )
 integer, intent(in) :: nbad
 real   , intent(in) :: temp(:,:)
 integer :: i, j

   print *, 'Bad temperatures (dimension 1): ', (check_1d(temp(i,:)),i=1,size(temp,1))
   print *, 'Bad temperatures (dimension 2): ', (check_1d(temp(:,j)),j=1,size(temp,2))
   call error_handler (nbad)

 end subroutine temp_check_2d

!--------------------------------------------------------------

 subroutine temp_check_3d ( nbad, temp )
 integer, intent(in) :: nbad
 real, intent(in)  :: temp(:,:,:)
 integer :: i, j, k

   print *, 'Bad temperatures (dimension 1): ', (check_2d(temp(i,:,:)),i=1,size(temp,1))
   print *, 'Bad temperatures (dimension 2): ', (check_2d(temp(:,j,:)),j=1,size(temp,2))
   print *, 'Bad temperatures (dimension 3): ', (check_2d(temp(:,:,k)),k=1,size(temp,3))
   call error_handler (nbad)

 end subroutine temp_check_3d

!#######################################################################

end module sat_vapor_pres_mod

! <INFO>

!   <REFERENCE>            
!     Smithsonian Meteorological Tables Page 350.
!   </REFERENCE>

!   <BUG>                  
!     No error checking is done to make sure that the size of the
!     input and output fields match.
!   </BUG>

!   <NOTE>
!     1. <B>Vectorization</B><BR/>
!        To create a vector version the lookup routines need to be modified.
!    The local variables: tmp, del, ind, should be changed to arrays
!    with the same size and order as input array temp.
!
!     2. <B>Construction of the <TT>ES</TT> tables</B><BR/>
!         The tables are constructed using the saturation vapor pressure (<TT>ES</TT>)
!    equations in the Smithsonian tables. The tables are valid between
!    -160C to +100C with increments at 1/10 degree. Between -160C and -20C
!    values of <TT>ES</TT> over ice are used, between 0C and 100C values of<TT> ES</TT>
!    over water are used, between -20C and 0C blended values of <TT>ES</TT>
!    (over water and over ice) are used.
!
!    There are three tables constructed: <TT>ES</TT>, first derivative 
!       (<TT>ES'</TT>), and
!    second derivative (<TT>ES</TT>'').  The ES table is constructed directly from
!    the equations in the Smithsonian tables. The <TT>ES</TT>' table is constructed
!    by bracketing temperature values at +/- 0.01 degrees. The <TT>ES</TT>'' table
!    is estimated by using centered differencing of the <TT>ES</TT>' table.
!
!     3. <B>Determination of <TT>es</TT> and <TT>es'</TT> from lookup tables</B><BR/>
!         Values of the saturation vapor pressure (<TT>es</TT>) and the 
!    derivative (<TT>es'</TT>) are determined at temperature (T) from the lookup 
!    tables (<TT>ES</TT>, <TT>ES'</TT>, <TT>ES''</TT>)
!    using the following formula.
!<PRE>
!    es (T) = ES(t) + ES'(t) * dt + 0.5 * ES''(t) * dt**2
!    es'(T) = ES'(t) + ES''(t) * dt
!
!    where     t = lookup table temperature closest to T
!             dt = T - t
!</PRE>
!
!     4. Internal (private) parameters<BR/>
!       These parameters can be modified to increase/decrease the size/range
!    of the lookup tables.
!<PRE>
!!    tcmin   The minimum temperature (in deg C) in the lookup tables.
!!              [integer, default: tcmin = -160]
!!
!!    tcmax   The maximum temperature (in deg C) in the lookup tables.
!!              [integer, default: tcmin = +100]
!!</PRE>
!!   </NOTE>
!
!!   <TESTPROGRAM NAME="test_sat_vapor_pres">
!<PRE>
!use sat_vapor_pres_mod
!implicit none
!
!integer, parameter :: ipts=500, jpts=100, kpts=50, nloop=1
!real, dimension(ipts,jpts,kpts) :: t,es,esn,des,desn
!integer :: n
!
!! generate temperatures between 120K and 340K
!  call random_number (t)
!  t = 130. + t * 200.
!
!! initialize the tables (optional)
!  call sat_vapor_pres_init
!
!! compute actual es and "almost" actual des
!   es = compute_es  (t)
!  des = compute_des (t)
!
!do n = 1, nloop
!! es and des
!  call lookup_es  (t, esn)
!  call lookup_des (t,desn)
!enddo
!
!! terminate, print deviation from actual
!  print *, 'size=',ipts,jpts,kpts,nloop
!  print *, 'err es  = ', sum((esn-es)**2)
!  print *, 'err des = ', sum((desn-des)**2)
!
!contains
!
!!----------------------------------
!! routine to estimate derivative
!
! function compute_des (tem) result (des)
! real, intent(in) :: tem(:,:,:)
! real, dimension(size(tem,1),size(tem,2),size(tem,3)) :: des,esp,esm
! real, parameter :: tdel = .01
!    esp = compute_es (tem+tdel)
!    esm = compute_es (tem-tdel)
!    des = (esp-esm)/(2*tdel)
! end function compute_des
!!----------------------------------
!
!end program test_sat_vapor_pres
!</PRE>
!   </TESTPROGRAM>
! </INFO>
