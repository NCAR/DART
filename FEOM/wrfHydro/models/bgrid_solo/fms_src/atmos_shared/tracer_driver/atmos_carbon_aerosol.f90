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
module atmos_carbon_aerosol_mod
! <CONTACT EMAIL="wfc@gfdl.noaa.gov">
!   William Cooke
! </CONTACT>

! <REVIEWER EMAIL="lwh@gfdl.noaa.gov">
!   Larry Horowitz
! </REVIEWER>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

! <OVERVIEW>
!     This code allows the implementation of black and organic carbon 
!     tracers in the FMS framework.
! </OVERVIEW>

! <DESCRIPTION>
!   This module presents the method of Cooke et al. (1999, 2002) 

!   In its present implementation the black and organic carbon tracers 
!   are from the combustion of fossil fuel.

!   While the code here should provide insights into the carbonaceous 
!   aerosol cycle it is provided here more as an example of how to implement 
!   a tracer module in the FMS infrastructure. The parameters of the model 
!   should be checked and set to values corresponding to previous works if 
!   a user wishes to try to reproduce those works.
! </DESCRIPTION>

! <DATASET NAME="Black carbon emissions">
!   The black carbon emission dataset is that derived in Cooke et al. (1999)
!   The dataset can be obtained from the contact person above.
! </DATASET>  

! <DATASET NAME="Organic carbon emissions">
!   The organic carbon emission dataset is that derived in Cooke et al. (1999)
!   The dataset can be obtained from the contact person above.
! </DATASET>  

! <INFO>

!   <REFERENCE>
!Cooke, W. F. and J. J. N. Wilson,  A global black carbon aerosol model,
!J. Geophys. Res., 101, 19395-19409, 1996.
! </REFERENCE>
!   <REFERENCE>
!Cooke, W. F., C. Liousse, H. Cachier and J. Feichter, 
!Construction of a 1 x 1 fossil fuel emission dataset for carbonaceous
!aerosol and implementation and radiative impact in the ECHAM-4 model, 
!J. Geophys. Res., 104, 22137-22162, 1999 </REFERENCE>
!   <REFERENCE>
! Cooke, W.F., V. Ramaswamy and P. Kasibathla, 
! A GCM study of the global carbonaceous aerosol distribution.
!J. Geophys. Res., 107, accepted, 2002
! </REFERENCE>
! </INFO>
use types_mod, only : r8
use              fms_mod, only : file_exist,           &
                                 close_file,           &
                                 stdlog,               &
                                 write_version_number
use     time_manager_mod, only : time_type
!use     diag_manager_mod, only : send_data,            &
!                                 register_diag_field,  &
!                                 register_static_field
use   tracer_manager_mod, only : get_tracer_index, &
                                 set_tracer_atts
use    field_manager_mod, only : MODEL_ATMOS
use atmos_tracer_utilities_mod, only : interp_emiss
use        constants_mod, only : PI

implicit none
private
!-----------------------------------------------------------------------
!----- interfaces -------

public  atmos_blackc_sourcesink,   &
        atmos_organic_sourcesink,  &
        atmos_carbon_aerosol_init, &
        atmos_carbon_aerosol_end

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
!  Add id_tracername for each additional tracer. These are used in
!  initializing and outputting the tracer fields.
!
!-----------------------------------------------------------------------

! tracer number for radon
integer :: nbcphobic=0
integer :: nbcphilic=0
integer :: nocphobic=0
integer :: nocphilic=0

!--- identification numbers for  diagnostic fields and axes ----

integer :: id_emissoc, id_emissbc 

!--- Arrays to help calculate tracer sources/sinks ---
real(r8), allocatable, dimension(:,:) :: bcsource,ocsource

character(len=6), parameter :: module_name = 'tracer'

logical :: module_is_initialized = .FALSE.
logical :: used

!---- version number -----
character(len=128) :: version = '$Revision$'
character(len=128) :: tagname = '$Id$'
!-----------------------------------------------------------------------

contains

!#######################################################################

!<SUBROUTINE NAME ="atmos_blackc_sourcesink">
!<OVERVIEW>
!  A subroutine to calculate the source and sinks of black carbon aerosol.
!</OVERVIEW>
!
!<DESCRIPTION>
! 
! This routine calculates the source and sink terms for black carbon.
! Simply put, the hydrophobic aerosol has sources from emissions and 
! sinks from dry deposition and transformation into hydrophilic aerosol.
! The hydrophilic aerosol also has emission sources and has sinks of wet 
! and dry deposition.
!
! The following schematic shows how the black carbon scheme 
! is implemented. 

!<PRE>
! +------------+  Trans-   +------------+
! |  Hydro-    | formation |  Hydro-    |
! |  phobic    |           |  philic    |
! |  black     |---------->|  black     |
! |  carbon    |           |  carbon    |
! |            |           |            |
! +------------+           +------------+
!    ^      |                ^    |   |
!    |      |                |    |   |
!    |      =                |    =   =
!  Source  Dry            Source Dry Wet
!          Dep.                  Dep Dep
!  
!</PRE>

! The transformation time used here is 1 day, which corresponds to an
! e-folding time of 1.44 days. This can be varied as necessary.

!</DESCRIPTION>
!<TEMPLATE>
!call atmos_blackc_sourcesink (lon, lat, land, pwt, &
!                         black_cphob, black_cphob_dt,  &
!                         black_cphil, black_cphil_dt,  &
!                         Time, is, ie, js, je, kbot)
!</TEMPLATE>
!   <IN NAME="lon" TYPE="real(r8)" DIM="(:,:)">
!     Longitude of the centre of the model gridcells
!   </IN>
!   <IN NAME="lat" TYPE="real(r8)" DIM="(:,:)">
!     Latitude of the centre of the model gridcells
!   </IN>
!   <IN NAME="land" TYPE="real(r8)" DIM="(:,:)">
!     Land/sea mask.
!   </IN>
!   <IN NAME="pwt" TYPE="real(r8)" DIM="(:,:,:)">
!     The pressure weighting array. = dP/grav
!   </IN>
!   <IN NAME="black_cphob" TYPE="real(r8)" DIM="(:,:,:)">
!     The array of the hydrophobic black carbon aerosol mixing ratio
!   </IN>
!   <IN NAME="black_cphil" TYPE="real(r8)" DIM="(:,:,:)">
!     The array of the hydrophilic black carbon aerosol mixing ratio
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

!   <OUT NAME="black_cphob_dt" TYPE="real(r8)" DIM="(:,:,:)">
!     The array of the tendency of the hydrophobic black carbon aerosol mixing ratio.
!   </OUT>
!   <OUT NAME="black_cphil_dt" TYPE="real(r8)" DIM="(:,:,:)">
!     The array of the tendency of the hydrophilic black carbon aerosol mixing ratio.
!   </OUT>

 subroutine atmos_blackc_sourcesink (lon, lat, land, pwt, &
                               black_cphob, black_cphob_dt,  &
                               black_cphil, black_cphil_dt,  &
                               Time, is, ie, js, je, kbot)

!-----------------------------------------------------------------------
   real(r8), intent(in),  dimension(:,:)   :: lon, lat
   real(r8), intent(in),  dimension(:,:)   :: land
   real(r8), intent(in),  dimension(:,:,:) :: pwt, black_cphob,black_cphil
   real(r8), intent(out), dimension(:,:,:) :: black_cphob_dt,black_cphil_dt
type(time_type), intent(in)            :: Time
integer, intent(in)                    :: is, ie, js, je
integer, intent(in),  dimension(:,:), optional :: kbot
!-----------------------------------------------------------------------
   real(r8), dimension(size(black_cphob,1),size(black_cphob,2),size(black_cphob,3)) ::  &
         sourcephob, sinkphob, sourcephil, sinkphil
   real(r8)  dtr
integer  i,j,kb,id,jd,kd,lat1
!-----------------------------------------------------------------------

      id=size(black_cphob,1); jd=size(black_cphob,2); kd=size(black_cphob,3)

      dtr= PI/180.

!----------- compute black carbon source ------------

      sourcephob = 0.0
      sourcephil = 0.0

          do j=1,jd
           sourcephob(:,j,kd)=0.8*bcsource(:,j+js-1)/pwt(:,j,kd)
           sourcephil(:,j,kd)=0.2*bcsource(:,j+js-1)/pwt(:,j,kd) +&
                              8.038e-6*black_cphob(:,j,kd)
          enddo


!------- compute black carbon phobic sink --------------
!
!  BCphob has a half-life time of 1.0days 
!   (corresponds to an e-folding time of 1.44 days)
!
!  sink = 1./(86400.*1.44) = 8.023e-6
!

    where (black_cphob(:,:,:) >= 0.0)
       sinkphob(:,:,:) = -8.038e-6*black_cphob(:,:,:)
    elsewhere
       sinkphob(:,:,:) = 0.0
    endwhere

       sinkphil(:,:,:) = 0.0

!------- tendency ------------------

      black_cphob_dt=sourcephob+sinkphob
      black_cphil_dt=sourcephil+sinkphil
      

!-----------------------------------------------------------------------

 end subroutine atmos_blackc_sourcesink
!</SUBROUTINE >

!#######################################################################

!<SUBROUTINE NAME ="atmos_organic_sourcesink">
!<OVERVIEW>
!  A subroutine to calculate the source and sinks of organic carbon aerosol.
!</OVERVIEW>
!<DESCRIPTION>

! This routine calculates the source and sink terms for organic carbon.
! Simply put, the hydrophobic aerosol has sources from emissions and 
! sinks from dry deposition and transformation into hydrophilic aerosol.
! The hydrophilic aerosol also has emission sources and has sinks of wet 
! and dry deposition.
!
! The following schematic shows how the organic carbon scheme 
! is implemented. 

!<PRE>
! +------------+  Trans-   +------------+
! |  Hydro-    | formation |  Hydro-    |
! |  phobic    |           |  philic    |
! |  organic   |---------->|  organic   |
! |  carbon    |           |  carbon    |
! |            |           |            |
! +------------+           +------------+
!    ^      |                ^    |   |
!    |      |                |    |   |
!    |      =                |    =   =
!  Source  Dry            Source Dry Wet
!          Dep.                  Dep Dep
!</PRE>
!
! The transformation time used here is 2 days, which corresponds to an
! e-folding time of 2.88 days. This can be varied as necessary.
!  
!</DESCRIPTION>
!<TEMPLATE>
!call atmos_organic_sourcesink (lon, lat, land, pwt, organic_carbon, organic_carbon_dt,  &
!                              Time, is, ie, js, je, kbot)
!</TEMPLATE>
!   <IN NAME="lon" TYPE="real(r8)" DIM="(:,:)">
!     Longitude of the centre of the model gridcells
!   </IN>
!   <IN NAME="lat" TYPE="real(r8)" DIM="(:,:)">
!     Latitude of the centre of the model gridcells
!   </IN>
!   <IN NAME="land" TYPE="real(r8)" DIM="(:,:)">
!     Land/sea mask.
!   </IN>
!   <IN NAME="pwt" TYPE="real(r8)" DIM="(:,:,:)">
!     The pressure weighting array. = dP/grav
!   </IN>
!   <IN NAME="organic_carbon" TYPE="real(r8)" DIM="(:,:,:)">
!     The array of the organic carbon aerosol mixing ratio
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

!   <OUT NAME="organic_carbon_dt" TYPE="real(r8)" DIM="(:,:,:)">
!     The array of the tendency of the organic carbon aerosol mixing ratio.
!   </OUT>

 subroutine atmos_organic_sourcesink (lon, lat, land, pwt, organic_carbon, organic_carbon_dt,  &
                              Time, is, ie, js, je, kbot)

!-----------------------------------------------------------------------
   real(r8), intent(in),  dimension(:,:)   :: lon, lat
   real(r8), intent(in),  dimension(:,:)   :: land
   real(r8), intent(in),  dimension(:,:,:) :: pwt, organic_carbon
   real(r8), intent(out), dimension(:,:,:) :: organic_carbon_dt
     type(time_type), intent(in) :: Time
integer, intent(in)                    :: is, ie, js, je 
integer, intent(in),  dimension(:,:), optional :: kbot
!-----------------------------------------------------------------------
   real(r8), dimension(size(organic_carbon,1),size(organic_carbon,2),size(organic_carbon,3)) ::  &
         source, sink
   real(r8)  dtr
integer  i,j,kb,id,jd,kd,lat1
!-----------------------------------------------------------------------

      id=size(organic_carbon,1); jd=size(organic_carbon,2); kd=size(organic_carbon,3)

      dtr=PI/180.

!----------- compute organic carbon source ------------

      source = 0.0

      if (present(kbot)) then
          do j=1,jd
          do i=1,id
             kb=kbot(i,j)
             source(i,j,kb)=ocsource(i,j+js-1)/pwt(i,j,kb)
          enddo
          enddo
      else
          do j=1,je-js+1
           source(:,j,kd)= ocsource(:,j+js-1)/pwt(:,j,kd)
          enddo
      endif


!------- compute organic carbon sink --------------
!
!  OCphob has a half-life time of 2.0days 
!   (corresponds to an e-folding time of 2.88 days)
!
!  sink = 1./(86400.*2.88) = 4.019e-6
!

    where (organic_carbon(:,:,:) >= 0.0)
         sink(:,:,:) = -4.019e-6*organic_carbon(:,:,:)
    elsewhere
       sink(:,:,:) = 0.0
    endwhere

!------- tendency ------------------

      organic_carbon_dt=source+sink
      

!-----------------------------------------------------------------------

 end subroutine atmos_organic_sourcesink
!</SUBROUTINE>


!#######################################################################

!<SUBROUTINE NAME ="atmos_carbon_aerosol_init">

!<OVERVIEW>
! Subroutine to initialize the carbon aerosol module.
!</OVERVIEW>
!<DESCRIPTION>
! This subroutine querys the tracer manager to find the indices for the 
! various carbonaceous aerosol tracers. It also registers the emission 
! fields for diagnostic purposes.
!  
!</DESCRIPTION>
!<TEMPLATE>
!call atmos_carbon_aerosol_init (lonb, latb, r, axes, Time, mask)
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

 subroutine atmos_carbon_aerosol_init (lonb, latb, r, axes, Time, mask)

!-----------------------------------------------------------------------
!
!   r    = tracer fields dimensioned as (nlon,nlat,nlev,ntrace)
!   mask = optional mask (0. or 1.) that designates which grid points
!          are above (=1.) or below (=0.) the ground dimensioned as
!          (nlon,nlat,nlev).
!
!-----------------------------------------------------------------------
real(r8), dimension(:),    intent(in) :: lonb, latb
real(r8),            intent(inout), dimension(:,:,:,:) :: r
integer        , intent(in)                        :: axes(4)
type(time_type), intent(in)                        :: Time
real(r8),            intent(in),    dimension(:,:,:), optional :: mask

integer :: n

   if (module_is_initialized) return

!----- set initial value of carbon ------------

   n = get_tracer_index(MODEL_ATMOS,'bcphob')
   if (n>0) then
      nbcphobic = n
      call set_tracer_atts(MODEL_ATMOS,'bcphob','hphobic_bc','g/g')
      if (nbcphobic > 0 ) write (*,30) 'Hydrophobic BC',nbcphobic
      if (nbcphobic > 0 ) write (stdlog(),30) 'Hydrophobic BC',nbcphobic
   endif

   n = get_tracer_index(MODEL_ATMOS,'bcphil')
   if (n>0) then
      nbcphilic=n
      call set_tracer_atts(MODEL_ATMOS,'bcphil','hphilic_bc','g/g')
      if (nbcphilic > 0 ) write (*,30) 'Hydrophilic BC',nbcphilic
      if (nbcphilic > 0 ) write (stdlog(),30) 'Hydrophilic BC',nbcphilic
   endif

   n = get_tracer_index(MODEL_ATMOS,'ocphob')
   if (n>0) then
      nocphobic=n
      call set_tracer_atts(MODEL_ATMOS,'ocphob','hphobic_oc','g/g')
      if (nocphobic > 0 ) write (*,30) 'Hydrophobic OC',nocphobic
      if (nocphobic > 0 ) write (stdlog(),30) 'Hydrophobic OC',nocphobic
   endif

   n = get_tracer_index(MODEL_ATMOS,'ocphil')
   if (n>0) then
      nocphilic=n
      call set_tracer_atts(MODEL_ATMOS,'ocphil','hphilic_oc','g/g')
      if (nocphilic > 0 ) write (*,30) 'Hydrophilic OC',nocphilic
      if (nocphilic > 0 ) write (stdlog(),30) 'Hydrophilic OC',nocphilic
   endif

  30        format (A,' was initialized as tracer number ',i2)
      !Read in emission files
!
!   id_emissbc = register_static_field ( 'tracers',                    &
!                     'bcemiss', axes(1:2),       &
!                     'bcemiss', 'g/m2/s')
!   id_emissoc = register_static_field ( 'tracers',                    &
!                     'ocemiss', axes(1:2),       &
!                     'ocemiss', 'g/m2/s')
!
   allocate (bcsource(size(lonb)-1,size(latb)-1))
   allocate (ocsource(size(lonb)-1,size(latb)-1))
   call tracer_input(lonb, latb, Time)

   call write_version_number (version, tagname)
   module_is_initialized = .TRUE.


!-----------------------------------------------------------------------

end subroutine atmos_carbon_aerosol_init
!</SUBROUTINE>


!<SUBROUTINE NAME ="atmos_carbon_aerosol_end">
!<OVERVIEW>
!  The destructor routine for the carbon aerosol module.
!</OVERVIEW>
! <DESCRIPTION>
! This subroutine writes the version name to logfile and exits. 
! </DESCRIPTION>
!<TEMPLATE>
!call atmos_carbon_aerosol_end
!</TEMPLATE>
subroutine atmos_carbon_aerosol_end

   module_is_initialized = .FALSE.

 end subroutine atmos_carbon_aerosol_end
!</SUBROUTINE>


!#######################################################################
 subroutine tracer_input(lonb, latb, Time)
real(r8), dimension(:),    intent(in) :: lonb, latb
type(time_type),intent(in) :: Time

integer      :: i, j, unit, io
real(r8)         :: emiss
real(r8)         :: dtr, deg_90, deg_180, deg3p6, deg3!, modxdeg, modydeg
real(r8)         :: ZCARBONSEASON(12)
real(r8)         :: bcsource1(100,60)
logical :: opened
!
! This is the Rotty seaonality for fossil fuel emissions of sulfate.
!
      DATA ZCARBONSEASON/1.146,1.139,1.081,0.995,0.916,0.920,0.910, &
                         0.907,0.934,0.962,1.019,1.072/
!
      dtr= PI/180.
      deg_90= -90.*dtr; deg_180= -180.*dtr 
      ! -90 and -180 degrees are the southwest boundaries of the 
      ! emission  field you are reading in.
      deg3p6 = 3.6*dtr; deg3 = 3.*dtr;
      ! 3.6 degrees longitude and 3 degree latitude is the resolution 
      ! of the r30 emission data that I used in SKYHI.

! initialise the BC phobic and philic
! read in the emission sources here. 

do unit = 30,100
INQUIRE(unit=unit, opened= opened)
if (.NOT. opened) exit
enddo
         open (unit,file='INPUT/r30.bc.ann', form='formatted', action='read')
        do io= 1,6000
           read  (unit, FMT=1968, end=11) i,j,emiss
           bcsource1(i,j)=emiss
        enddo
  11    call close_file (unit)
1968  FORMAT(2I3,e11.4)
1969  FORMAT(2I3,f11.3)
! Interpolate the R30 emission field to the resolution of the model.
        call interp_emiss ( bcsource1, 0.0_r8, deg_90, deg3p6, deg3, &
                     bcsource)
                    

 write(*,*) 'Reading OC emissions'
!Now let's do the OC 
!
         bcsource1 = 0.0E+00
do unit = 30,100
INQUIRE(unit=unit, opened= opened)
if (.NOT. opened) exit
enddo
         open (unit,file='INPUT/r30.oc.ann', form='formatted', action='read')
         do io= 1,6000
            read  (unit, FMT=1968, end=13) i,j,emiss
            bcsource1(i,j)=emiss
         enddo
  13     call close_file (unit)

! Interpolate the R30 emission field to the resolution of the model.
         call interp_emiss ( bcsource1, 0.0_r8, deg_90, deg3p6, deg3, &
                              ocsource)
          
! Send the emission data to the diag_manager for output.
!         if (id_emissbc > 0 ) &
!           used = send_data ( id_emissbc, bcsource, Time )
!         if (id_emissoc > 0 ) &
!           used = send_data ( id_emissoc, ocsource, Time )

end subroutine tracer_input

end module atmos_carbon_aerosol_mod



