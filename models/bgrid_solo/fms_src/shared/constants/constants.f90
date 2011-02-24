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

module constants_mod

! <CONTACT EMAIL="bw@gfdl.noaa.gov">
!   Bruce Wyman
! </CONTACT>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

! <OVERVIEW>
!    Defines useful constants for Earth in mks units.
! </OVERVIEW>

! <DESCRIPTION>
!   Constants are defined as real(r8) parameters.They are accessed through
!   the "use" statement. While the local name of constant may be changed,
!   their values can not be redefined.
! </DESCRIPTION>

!-----------------------------------------------------
!
!  Defines useful constants for Earth in mks units.
!
!-----------------------------------------------------

use types_mod, only : r8
use fms_mod, only: write_version_number

implicit none
private

character(len=128) :: version='$Revision$'
character(len=128) :: tagname='$Id$'
logical :: do_log = .true.
logical :: module_is_initialized = .FALSE.
!-----------------------------------------------------------------------

!------------ physical constants ---------------
! <DATA NAME="RADIUS" UNITS="meters" TYPE="real(r8)" DEFAULT="6376.e3">
!   radius of the earth
! </DATA>
! <DATA NAME="OMEGA" UNITS="1/sec" TYPE="real(r8)" DEFAULT="7.292e-5">
!   rotation rate of the planet (earth)
! </DATA>
! <DATA NAME="GRAV" UNITS="m/s2" TYPE="real(r8)" DEFAULT="9.80">
!   acceleration due to gravity
! </DATA>
! <DATA NAME="RDGAS" UNITS="J/kg/deg" TYPE="real(r8)" DEFAULT="287.04">
!   gas constant for dry air
! </DATA>
! <DATA NAME=" KAPPA" UNITS="" TYPE="real(r8)" DEFAULT="2./7.">
!   RDGAS / CP
! </DATA>
! <DATA NAME="CP" UNITS="J/kg/deg" TYPE="real(r8)" DEFAULT="RDGAS/KAPPA">
!   specific heat capacity of dry air at constant pressure
! </DATA>

real(r8), public, parameter :: RADIUS = 6376.e3   
real(r8), public, parameter :: OMEGA  = 7.292e-5 
real(r8), public, parameter :: GRAV   = 9.80    
real(r8), public, parameter :: RDGAS  = 287.04 
real(r8), public, parameter :: KAPPA  = 2./7.  
real(r8), public, parameter :: cp     = RDGAS/KAPPA 

!------------ water vapor constants ---------------
! <DATA NAME="RVGAS" UNITS="J/Kg/deg" TYPE="real(r8)" DEFAULT="461.50">
!   gas constant for water vapor
! </DATA>
! <DATA NAME="DENS_H2O" UNITS="Kg/m3" TYPE="real(r8)" DEFAULT="1000.">
!   density of liquid water
! </DATA>
! <DATA NAME="HLV" UNITS="J/Kg" TYPE="real(r8)" DEFAULT="2.500e6">
!   latent heat of evaporation
! </DATA>
! <DATA NAME="HLF" UNITS="J/kg" TYPE="real(r8)" DEFAULT="3.34e5">
!   latent heat of fusion
! </DATA>
! <DATA NAME="HLS" UNITS="J/Kg" TYPE="real(r8)" DEFAULT="2.834e6">
!   latent heat of sublimation
! </DATA>
! <DATA NAME="TFREEZE" UNITS="deg K" TYPE="real(r8)" DEFAULT="273.16">
!   temp where fresh water freezes
! </DATA>

real(r8), public, parameter :: RVGAS = 461.50 
real(r8), public, parameter :: DENS_H2O = 1000. 
real(r8), public, parameter :: HLV = 2.500e6   
real(r8), public, parameter :: HLF = 3.34e5   
real(r8), public, parameter :: HLS = 2.834e6 
real(r8), public, parameter :: TFREEZE = 273.16    

!------------ miscellaneous constants ---------------
! <DATA NAME="STEFAN" UNITS="(W/m2/deg4" TYPE="real(r8)" DEFAULT="5.6734e-8">
!   Stefan-Boltzmann constant
! </DATA>
! <DATA NAME="VONKARM"  TYPE="real(r8)" DEFAULT="0.40">
!   Von Karman constant
! </DATA>
! <DATA NAME="PI" TYPE="real(r8)" DEFAULT="3.14159265358979323846">
!   is it enough?
! </DATA>

real(r8), public, parameter :: STEFAN  =  5.6734e-8 
real(r8), public, parameter :: VONKARM =  0.40     
real(r8), public, parameter :: PI      =  3.14159265358979323846 
!-----------------------------------------------------------------------

public constants_init

contains

! <SUBROUTINE NAME="constants_init">

!   <OVERVIEW>
!     A optional initialization routine. The only purpose of this routine 
!     is to write the version and tag name information to the log file.   
!   </OVERVIEW>
!   <DESCRIPTION>
!     The only purpose of this routine is to write the version and
!     tag name information to the log file.  This routine does not
!     have to be called.  If it is called more than once or called
!     from other than the root PE it will return silently.
!     There are no arguments.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call constants_init
!   </TEMPLATE>

subroutine constants_init

  if (module_is_initialized) return
  module_is_initialized = .TRUE.

  if (.not.do_log) return
  call write_version_number (version,tagname)
  do_log = .false.

end subroutine constants_init
! </SUBROUTINE>

!-----------------------------------------------------------------------

end module constants_mod

! <INFO>

!   <FUTURE>               
!   1.  Renaming of constants.
!   </FUTURE>               
!   <FUTURE>               
!   2.  Additional constants.
!   </FUTURE>
!   <NOTE>
!      <B>NOTES ON USAGE:</B>
!
!    All constants have been declared as type REAL, PARAMETER.
!
!    The value a constant can not be changed in a users program.
!    New constants can be defined in terms of values from the
!    constants module using a parameter statement.<br><br>
!
!    The name given to a particular constant may be changed.<br><br>
!
!    Constants can be used on the right side on an assignment statement
!   (their value can not be reassigned). 
!
!<B>EXAMPLES:</B>
!<PRE>
!     use constants_mod, only:  TFREEZE, grav_new =&gt; GRAV
!     real(r8), parameter :: grav_inv = 1.0 / grav_new
!     tempc(:,:,:) = tempk(:,:,:) - TFREEZE
!     geopotential(:,:) = height(:,:) * grav_new
!</PRE>

!   </NOTE>

! </INFO>

