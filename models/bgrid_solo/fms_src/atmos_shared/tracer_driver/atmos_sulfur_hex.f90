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

module atmos_sulfur_hex_mod
! <CONTACT EMAIL="jbg@gfdl.noaa.gov">
!   Jeff Greenblatt
! </CONTACT>

! <REVIEWER EMAIL="wfc@gfdl.noaa.gov">
!   William Cooke
! </REVIEWER>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

! <OVERVIEW>
!     This code allows the implementation of sulfur hexafluoride
!     tracer in the FMS framework.
! </OVERVIEW>

! <DESCRIPTION>

! </DESCRIPTION>

! <DATASET NAME="Sulfur hexaflouride emissions">
!
! Monthly.emissions contains the estimated global emission rate of SF6 in
! Gg/yr for 62 months between December 1988 and January 1994, inclusive.
! These are based on the annual estimates of Levin and Hesshaimer
! (submitted), and have been linearly interpolated to monthly values. The
! last half of 1993 has been extrapolated using the trend for the previous 12
! months. 
!
!   The dataset can be obtained from the contact person above.
! </DATASET>  
! <INFO>

!   <REFERENCE>
!Levin, I. and V. Hessahimer: Refining of atmospheric
! transport model entries by the globally observed passive tracer
! distributions of 85Krypton and Sulfur Hexafluoride (SF6). Submitted to the
! Journal of Geophysical Research.
! </REFERENCE>
!</INFO>

use types_mod, only : r8
use              fms_mod, only : file_exist,           &
!                                 open_file,            &
                                 stdlog,               &
                                 close_file,           &
                                 write_version_number
use     time_manager_mod, only : time_type,            &
                                 set_date,             &
                                 operator( > ),        &
                                 operator( < ),        &
                                 operator( >= )
!use     diag_manager_mod, only : send_data,            &
!                                 register_diag_field,  &
!                                 register_static_field
use   tracer_manager_mod, only : get_tracer_index
use    field_manager_mod, only : MODEL_ATMOS
use atmos_tracer_utilities_mod, only : interp_emiss
use        constants_mod, only : grav, PI

implicit none
private
!-----------------------------------------------------------------------
!----- interfaces -------

public  atmos_sf6_sourcesink, atmos_sulfur_hex_init, atmos_sulfur_hex_end

!-----------------------------------------------------------------------
!----------- namelist -------------------
!-----------------------------------------------------------------------
!
!  When initializing additional tracers, the user needs to make the
!  following changes.
!
!  Add an integer variable below for each additional tracer. This should
!  be initialized to zero. 
!
!-----------------------------------------------------------------------

! tracer number for radon
integer :: nsf6     =0


!--- identification numbers for  diagnostic fields and axes ----

integer :: id_emiss

!--- Arrays to help calculate tracer sources/sinks ---
real(r8), allocatable, dimension(:,:) :: sf6_grid

integer, parameter :: NUM_SF6_RATE = 62 !number of entries in file 'monthly.emissions'
type sf6_rate_type
  type(time_type) :: Time
  real(r8) :: rate
end type sf6_rate_type
type(sf6_rate_type), dimension(NUM_SF6_RATE) :: sf6_rate

character(len=6), parameter :: module_name = 'tracer'

logical :: module_is_initialized=.FALSE.
logical :: used

!---- version number -----
character(len=128) :: version = '$Revision$'
character(len=128) :: tagname = '$Id$'
!-----------------------------------------------------------------------

contains

!#######################################################################
!<SUBROUTINE NAME="atmos_sf6_sourcesink">
!<OVERVIEW>
! A routine to calculate the sources and sinks of sulfur hexafluoride.
!</OVERVIEW>
!<DESCRIPTION>
! A routine to calculate the sources and sinks of sulfur hexafluoride.
!</DESCRIPTION>
!<TEMPLATE>
!call atmos_sf6_sourcesink (lon, lat, land, pwt, sf6, sf6_dt, 
!        Time, is, ie, js, je, kbot)
!
!</TEMPLATE>
!   <IN NAME="lon" TYPE="real(r8)" DIM="(:,:)">
!     Longitude of the centre of the model gridcells.
!   </IN>
!   <IN NAME="lat" TYPE="real(r8)" DIM="(:,:)">
!     Latitude of the centre of the model gridcells.
!   </IN>
!   <IN NAME="land" TYPE="real(r8)" DIM="(:,:)">
!     Land/sea mask.
!   </IN>
!   <IN NAME="pwt" TYPE="real(r8)" DIM="(:,:,:)">
!     The pressure weighting array. = dP/grav
!   </IN>
!   <IN NAME="sf6" TYPE="real(r8)" DIM="(:,:,:)">
!     The array of the sulfur hexafluoride mixing ratio.
!   </IN>
!   <IN NAME="Time" TYPE="type(time_type)">
!     Model time.
!   </IN>
!   <IN NAME="is, ie, js, je" TYPE="integer">
!     Local domain boundaries.
!   </IN>
!   <IN NAME="kbot" TYPE="integer, optional" DIM="(:,:)">
!     Integer array describing which model layer intercepts the surface.
!   </IN>

!   <OUT NAME="sf6_dt" TYPE="real(r8)" DIM="(:,:,:)">
!     The array of the tendency of the sulfur hexafluoride mixing ratio.
!   </OUT>
!
subroutine atmos_sf6_sourcesink (lon, lat, land, pwt, sf6, sf6_dt, &
        Time, is, ie, js, je, kbot)
!-----------------------------------------------------------------------
   real(r8), intent(in),  dimension(:,:)   :: lon, lat
   real(r8), intent(in),  dimension(:,:)   :: land
   real(r8), intent(in),  dimension(:,:,:) :: pwt, sf6
   real(r8), intent(out), dimension(:,:,:) :: sf6_dt
     type(time_type), intent(in) :: Time     
integer, intent(in),  dimension(:,:), optional :: kbot
integer, intent(in)                    :: is, ie, js, je
!-----------------------------------------------------------------------
   real(r8), dimension(size(sf6,1),size(sf6,2),size(sf6,3)) ::  &
         source, sink
logical, dimension(size(sf6,1),size(sf6,2)) ::  maskeq,masknh
integer :: i,j,kb,kd, id,jd
real(r8) :: rate ! sf6 interpolated emission rate
!-----------------------------------------------------------------------

      id=size(sf6,1); jd=size(sf6,2); kd=size(sf6,3)

source=0.0
! Interpolate SF6 global emission rate from sf6_rate (time dependent). For now
! just use first or last entry if time falls outside the bounds of the table.

      if (Time < sf6_rate(1)%Time) then
        rate=0. !previously sf6_rate(1)%Rate
      else
        if (Time > sf6_rate(size(sf6_rate))%Time) then
          rate=sf6_rate(size(sf6_rate))%Rate !just keep fixed past end of array
        else
          do i=1,size(sf6_rate)-1 !This can be optimized with efficient search
            if (Time >= sf6_rate(i)%Time .and. Time < sf6_rate(i+1)%Time) then
              rate=sf6_rate(i)%Rate
              exit
            endif
          enddo
        endif
      endif

      if (present(kbot)) then
          do j=1,jd
          do i=1,id
             kb=kbot(i,j)
             source(i,j,kb)=sf6_grid(i,j+js-1)*rate/pwt(i,j,kb)
          enddo
          enddo
      else
          do j=1,jd
            source(:,j,kd)=sf6_grid(:,j+js-1)*rate/pwt(:,j,kd)
          enddo
      endif
      
      sink=0.0

      sf6_dt=source+sink

end subroutine atmos_sf6_sourcesink
!</SUBROUTINE>

!#######################################################################
!<SUBROUTINE NAME="atmos_sulfur_hex_init">
!<OVERVIEW>
! The constructor routine for the sulfur hexafluoride module.
!</OVERVIEW>
!<DESCRIPTION>
! A routine to initialize the sulfur hexafluoride module.
!</DESCRIPTION>
!<TEMPLATE>
!call atmos_sulfur_hex_init (lonb, latb, r, axes, Time, mask)
!</TEMPLATE>
!   <IN NAME="lonb" TYPE="real(r8)" DIM="(:)">
!     The longitudes for the local domain.
!   </IN>
!   <IN NAME="latb" TYPE="real(r8)" DIM="(:)">
!     The latitudes for the local domain.
!   </IN>
!   <INOUT NAME="r" TYPE="real(r8)" DIM="(:,:,:,:)">
!     Tracer fields dimensioned as (nlon,nlat,nlev,ntrace). 
!   </INOUT>
!   <IN NAME="mask" TYPE="real(r8), optional" DIM="(:,:,:)">
!      optional mask (0. or 1.) that designates which grid points
!           are above (=1.) or below (=0.) the ground dimensioned as
!           (nlon,nlat,nlev).
!   </IN>
!   <IN NAME="Time" TYPE="type(time_type)">
!     Model time.
!   </IN>
!   <IN NAME="axes" TYPE="integer" DIM="(4)">
!     The axes relating to the tracer array dimensioned as
!      (nlon, nlat, nlev, ntime)
!   </IN>

 subroutine atmos_sulfur_hex_init (lonb, latb, r, axes, Time, mask)

!-----------------------------------------------------------------------
!
!   r    = tracer fields dimensioned as (nlon,nlat,nlev,ntrace)
!   mask = optional mask (0. or 1.) that designates which grid points
!          are above (=1.) or below (=0.) the ground dimensioned as
!          (nlon,nlat,nlev).
!
!-----------------------------------------------------------------------
real(r8),            intent(in),    dimension(:)               :: lonb, latb
real(r8),            intent(inout), dimension(:,:,:,:)         :: r
real(r8),            intent(in),    dimension(:,:,:), optional :: mask
type(time_type), intent(in)                                :: Time
integer        , intent(in)                                :: axes(4)

logical :: flag
integer :: n
!-----------------------------------------------------------------------
!
!  When initializing additional tracers, the user needs to make changes 
!  to two namelists. 
!
!  In what is core_namelist in the run_script (atmos_*_input.nml)
!  ntrace and ntprog need to be changed in &atmosphere_nml so that 
!  the number is the numbers of tracers in this module plus 
!  one (specific humidity I believe)
!
!  In what is phys_namelist in the run_script ( atmos_param_*_input.nml)
!  the namelist &tracer_driver_nml needs to be extended with numbers 
!  corresponding to each tracer. Thes numbers should be positive and 
!  non-zero.
!
!-----------------------------------------------------------------------
      integer  log_unit,unit,io,index,ntr,nt

      if (module_is_initialized) return

!---- write namelist ------------------

      call write_version_number (version, tagname)

      n = get_tracer_index(MODEL_ATMOS,'sf6')
      if (n>0) then
        nsf6=n
        if (nsf6 > 0 ) write (*,30) 'SF6',nsf6
        if (nsf6 > 0 ) write (stdlog(),30) 'SF6',nsf6
      endif

  30        format (A,' was initialized as tracer number ',i2)
      !Read in emission files
      

!     id_emiss = register_static_field ( 'tracers',                    &
!                     'sf6emiss', axes(1:2),       &
!                     'sulfhexemiss', 'g/m2/s')

   allocate (sf6_grid(size(lonb)-1,size(latb)-1))


      call sf6_init(Time)

      module_is_initialized = .TRUE.
!-----------------------------------------------------------------------
 end subroutine atmos_sulfur_hex_init
!</SUBROUTINE>



!######################################################################

subroutine sf6_init(Time)
type(time_type), intent(in) :: Time
!-------------------------------------------------
!-------------------------------------------------
      integer      :: i,j,unit !,imon,irec,n,io
      real(r8)         :: dtr,deg_90, deg_180, gxdeg, gydeg
      real(r8)         :: GEIA(720, 360)
      integer, parameter :: k6=selected_int_kind(6) ! find kind sufficient for 6 digit integer precision
      integer(kind=k6) :: t ! temporary time variable of kind k6
      integer      :: y,m,d ! calendar vars
      real(r8),pointer         :: data_out1(:,:) ! temporary sf6_grid output
!-------------------------------------------------
      real(r8) :: MW_air=28.9644 ! molecular wt. of air (gm/mole)
      real(r8) :: MW_sf6=86.0 ! molecular wt. of sf6 (gm/mole) PLEASE CHECK!
      logical :: used, opened

      dtr=PI/180.
      deg_90= -90.*dtr; deg_180= -180.*dtr ! -90 and -180 degrees are the southwest boundaries of the emission field you are reading in.

! Read in GEIA SF6 emission distribution grid and determine sizes:
!
! GEIA grid arrangement (taken from snf's togasp.f code documentation):
!     i  = 1 is centered at 179.75w
!     i increases eastward
!     j  = 1 is centered at 89.75s
!     j increases northward
!
!       88.5s - | - - - | - - - | - - - | - -
!               | (1,3) | (2,3) | (3,3) |
!       89s   - | - - - | - - - | - - - | - -
!               | (1,2) | (2,2) | (3,2) |
!       89.5s - | - - - | - - - | - - - | - -
!               | (1,1) | (2,1) | (3,1) |
!       90s -   | - - - | - - - | - - - |
!              180w   179.5w  179w    178.5w

do unit = 30,100
INQUIRE(unit=unit, opened= opened)
if (.NOT. opened) exit
enddo
      open(unit,file='distribution.grid', form='formatted',action='read')
      do j = 1, 360 !rearrange input array so begins at 0 E
        read(unit,'(5e16.8)') (GEIA(I,J), I=361,720)
        read(unit,'(5e16.8)') (GEIA(I,J), I=  1,360)
      end do
      call close_file(unit) 

      gxdeg=360./size(GEIA,1)*dtr
      gydeg=180./size(GEIA,2)*dtr

      call interp_emiss (GEIA, deg_180, deg_90, gxdeg, gydeg, sf6_grid)

! Note: must do scaling of global integral to 3.1828e-2 (for kg/m2/s) or 2750.
! (for kg/m2/day), according to S.-M. Fan.

! Scale to same units as used for radon (see radon_sourcesink):

      sf6_grid=sf6_grid * grav * MW_air / MW_sf6

! Now read in emission rate table. Comments from README file supplied by
! Song-Miao Fan:
! "Monthly.emissions contains the estimated global emission rate of SF6 in
! Gg/yr for 62 months between December 1988 and January 1994, inclusive.
! These are based on the annual estimates of Levin and Hesshaimer
! (submitted), and have been linearly interpolated to monthly values. The
! last half of 1993 has been extrapolated using the trend for the previous 12
! months. (Ref: Levin, I. and V. Hessahimer: Refining of atmospheric
! transport model entries by the globally observed passive tracer
! distributions of 85Krypton and Sulfur Hexafluoride (SF6). Submitted to the
! Journal of Geophysical Research)."
!
! Units are Gg/yr. When multiplied by sf6_grid cell, units will be kg/m2/s.
!
! Note time is in integer YYMMDD format, which must be converted to Time_type
! for storage in the sf6_rate array.
!
do unit = 30,100
INQUIRE(unit=unit, opened= opened)
if (.NOT. opened) exit
enddo
      open(unit,file='monthly.emissions', form='formatted',action='read')
      do j = 1, size(sf6_rate)
        read(unit,'(i6,2x,f7.5)') t, sf6_rate(j)%rate
! convert YYMMDD into components:
        y=int(t/10000)
        m=int((t-y*10000)/100)
        d=mod(t,100)
! shift start year to 1981 (start date 1981.12.15):
        y=y-88+1981 ! y was 2-digit year
! now convert to time_type format and store:
        sf6_rate(j)%Time=set_date(y, m, d)
      end do
      call close_file(unit) 

!         if (id_emiss > 0 ) &
!         used = send_data ( id_emiss, sf6_grid, Time )
         
      end subroutine sf6_init

!######################################################################
!<SUBROUTINE NAME="sulfur_hex_end">
!<OVERVIEW>
!  The destructor routine for the sulfur hexafluoride module.
!</OVERVIEW>
! <DESCRIPTION>
! This subroutine is the exit routine for the sulfur hexafluoride module.
! </DESCRIPTION>
!<TEMPLATE>
! call atmos_sulfur_hex_end
!</TEMPLATE>

 subroutine atmos_sulfur_hex_end
 
      module_is_initialized = .FALSE.

 end subroutine atmos_sulfur_hex_end
!</SUBROUTINE>


end module atmos_sulfur_hex_mod



