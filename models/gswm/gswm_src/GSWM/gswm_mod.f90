module gswm_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                       !!
!!                   GNU General Public License                          !!
!!                                                                       !!
!! This file is part of the Flexible Modeling System (FMS).              !!
!!                                                                       !!
!! FMS is free software; you can redistribute it and/or modify           !!
!! it and are expected to follow the terms of the GNU General Public     !!
!! License as published by the Free Software Foundation.                 !!
!!                                                                       !!
!! FMS is distributed in the hope that it will be useful,                !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of        !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         !!
!! GNU General Public License for more details.                          !!
!!                                                                       !!
!! You should have received a copy of the GNU General Public License     !!
!! along with FMS; if not, write to:                                     !!
!!          Free Software Foundation, Inc.                               !!
!!          59 Temple Place, Suite 330                                   !!
!!          Boston, MA  02111-1307  USA                                  !!
!! or see:                                                               !!
!!          http://www.gnu.org/licenses/gpl.txt                          !!
!!                                                                       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-----------------------------------------------------------------------
!
!     Global Scale Wave Model module routines
!
!-----------------------------------------------------------------------
!---------------- m o d u l e   i n f o r m a t i o n ------------------

use time_manager_mod, only: time_type, set_time, get_time,  &
                            operator(+), operator (<), operator (>), &
                            operator (/=), operator (/), operator (*)

use location_mod,     only: location_type, get_location, set_location, &
                            get_dist, vert_is_level, &
                            LocationDims, LocationName, LocationLName

use types_mod

!-----------------------------------------------------------------------

implicit none
private

public  get_model_size

!-----------------------------------------------------------------------
! let CVS fill strings ... DO NOT EDIT ...

character(len=128) :: version = "$Id$"
character(len=128) :: tag = "$Name$"

character(len=128) :: &
   source = "$Source$", &
   revision = "$Revision$", &
   revdate  = "$Date$"


!-----------------------------------------------------------------------
!---- namelist (saved in file input.nml) ----
!

   integer, dimension(2) :: spectral_coefficients = (/0.0, 0.0 /)

   namelist /gswm_nml/ spectral_coefficients

!-----------------------------------------------------------------------
! Additional stuff currently in main_nml 

   integer, dimension(4) :: current_time = (/ 0, 0, 0, 0 /)
   logical :: override = .false.
   integer :: days=0, hours=0, minutes=0, seconds=0
   integer :: dt_atmos = 0
   real    :: noise_sd = 0.0
   integer :: dt_bias = -1
   logical :: output_state_vector = .false.  ! output prognostic variables

   namelist /main_nml/ current_time, override, dt_atmos, &
                          days, hours, minutes, seconds, noise_sd, &
                          dt_bias, output_state_vector 

!-----------------------------------------------------------------------
! More stuff from atmos_solo driver
! ----- model time -----

   type (time_type) :: Time, Time_init, Time_end, Time_step_atmos

! ----- coupled model initial date -----

   integer :: date_init(6)

!-----------------------------------------------------------------------
! Public definition of variable types
integer, parameter :: TYPE_PS = 0, TYPE_T = 1, TYPE_U = 2, TYPE_V = 3, TYPE_TRACER = 4


!-----------------------------------------------------------------------
!---- private data ----

type (gswm_type) :: Dynam    ! for example

integer                            :: model_size
real,    dimension(:,:,:), pointer :: omega
integer, dimension(4)              :: atmos_axes
integer                            :: num_levels
integer                            :: ntracers 
real, dimension(:), pointer        :: v_lons, v_lats, t_lons, t_lats



! this is from setatmos.f
! common/interp1/zpt(37,101,3),zpu(37,101,3),zpr(37,101,3),                 
!    +zpxh(37,101,3),do2(37,101,3) 
integer, parameter :: numlat = 37, numlon = 101, numvert = 3
real(r8) ::  zpt(numlat,numlon,numvert),  zpu(numlat,numlon,numvert), &
             zpr(numlat,numlon,numvert), zpxh(numlat,numlon,numvert), &
             do2(numlat,numlon,numvert) 



CONTAINS
!#######################################################################







FUNCTION PSI(XI,ETA,DEC)
! This was rewritten from psi.f
    real(r8), INTENT(IN) :: XI
    real(r8), INTENT(IN) :: ETA
    real(r8), INTENT(IN) :: DEC
    real(r8)             :: PSI

    PSI=XI+DEC*COS(ETA)

    RETURN

END FUNCTION psi




!#######################################################################
END MODULE gswm_mod
