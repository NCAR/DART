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

module atmos_tracer_utilities_mod
! <CONTACT EMAIL="wfc@gfdl.noaa.gov">
!   William Cooke
! </CONTACT>

! <REVIEWER EMAIL="bw@gfdl.noaa.gov">
!   Bruce Wyman
! </REVIEWER>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

! <OVERVIEW>
!     This code provides some utility routines for atmospheric tracers in the FMS framework.
! </OVERVIEW>
! <DESCRIPTION>
!    This module gives utility routines which can be used to provide 
!    consistent removal mechanisms for atmospheric tracers. 
!
!    In particular it provides schemes for wet and dry deposiiton that 
!    can be easily utilized.
!
! </DESCRIPTION>

use types_mod, only : r8
use            fms_mod, only : lowercase, &
                               write_version_number, &
                               stdlog
use   time_manager_mod, only : time_type
!use   diag_manager_mod, only : send_data, &
!                               register_diag_field
use tracer_manager_mod, only : query_method, &
                               get_tracer_names, &
                               get_number_tracers
use  field_manager_mod, only : MODEL_ATMOS, parse
use      constants_mod, only : grav, rdgas, PI
use   horiz_interp_mod, only : horiz_interp

implicit none
private
!-----------------------------------------------------------------------
!----- interfaces -------

public  wet_deposition,    &
        dry_deposition,    &
        interp_emiss,      &
        atmos_tracer_utilities_end, &
        atmos_tracer_utilities_init


!---- version number -----
logical :: module_is_initialized = .FALSE.

character(len=128) :: version = '$Revision$'
character(len=128) :: tagname = '$Id$'

character(len=7), parameter :: mod_name = 'tracers'
!-----------------------------------------------------------------------
!--- identification numbers for  diagnostic fields and axes ----
integer, parameter :: max_tracers = 30
integer :: id_tracer_ddep(max_tracers), id_tracer_wdep(max_tracers)
character(len=32),  dimension(max_tracers) :: tracer_names     = ' '
character(len=32),  dimension(max_tracers) :: tracer_units     = ' '
character(len=128), dimension(max_tracers) :: tracer_longnames = ' '
character(len=32),  dimension(max_tracers) :: tracer_wdep_names     = ' '
character(len=32),  dimension(max_tracers) :: tracer_wdep_units     = ' '
character(len=128), dimension(max_tracers) :: tracer_wdep_longnames = ' '
character(len=32),  dimension(max_tracers) :: tracer_ddep_names     = ' '
character(len=32),  dimension(max_tracers) :: tracer_ddep_units     = ' '
character(len=128), dimension(max_tracers) :: tracer_ddep_longnames = ' '

real(r8), allocatable :: blon_out(:), blat_out(:)

contains

!
! ######################################################################
!
!<SUBROUTINE NAME="atmos_tracer_utilities_init">
!<OVERVIEW>
! This is a routine to create and register the dry and wet deposition 
! fields of the tracers.
!</OVERVIEW>
!<DESCRIPTION>
!  This routine creates diagnostic names for dry and wet deposition fields of the tracers.
!  It takes the tracer name and appends "ddep" for the dry deposition field and "wdep" for 
!  the wet deposition field. This names can then be entered in the diag_table for 
!  diagnostic output of the tracer dry and wet deposition. The module name associated with
!  these fields in "tracers". The units of the deposition fields are assumed to be kg/m2/s.
!</DESCRIPTION>
!<TEMPLATE>
! call atmos_tracer_utilities_init(lonb,latb, mass_axes, Time)
!</TEMPLATE>
!   <IN NAME="lonb" TYPE="real(r8)" DIM="(:)">
!     The longitudes for the local domain.
!   </IN>
!   <IN NAME="latb" TYPE="real(r8)" DIM="(:)">
!     The latitudes for the local domain.
!   </IN>
!   <IN NAME="mass_axes" TYPE="integer" DIM="(3)">
!     The axes relating to the tracer array.
!   </IN>
!   <IN NAME="Time" TYPE="type(time_type)">
!     Model time.
!   </IN>

subroutine atmos_tracer_utilities_init(lonb, latb, mass_axes, Time)

! Routine to initialize the tracer identification numbers. 
! This registers the 2D fields for the wet and dry deposition.
real(r8), dimension(:),    intent(in) :: lonb, latb
integer, dimension(3), intent(in) :: mass_axes
type(time_type),       intent(in) :: Time

integer :: ntrace
!
integer :: n, unit
character(len=128) :: name

! Make local copies of the local domain dimensions for use 
! in interp_emiss.
      allocate ( blon_out(size(lonb)))
      allocate ( blat_out(size(latb)))
!      allocate ( data_out(size(lonb)-1, size(latb)-1))
      blon_out = lonb
      blat_out = latb
      
      do n = 1, max_tracers
         write ( tracer_names(n),     100 ) n
         write ( tracer_longnames(n), 102 ) n
         tracer_units(n) = 'none'
      enddo
  100 format ('tr',i2.2)
  102 format ('tracer ',i2.2)

call get_number_tracers(MODEL_ATMOS, num_tracers= ntrace)
   do n = 1, ntrace
!--- set tracer tendency names where tracer names have changed ---

call get_tracer_names(MODEL_ATMOS,n,tracer_names(n),tracer_longnames(n),tracer_units(n))
      write (name,100) n
      if (trim(tracer_names(n)) /= name) then
          tracer_ddep_names(n) = trim(tracer_names(n)) // 'ddep'
          tracer_wdep_names(n) = trim(tracer_names(n)) // 'wdep'
      endif
      write (name,102) n
      if (trim(tracer_longnames(n)) /= name) then
          tracer_wdep_longnames(n) = &
                  trim(tracer_longnames(n)) // ' wet deposition for tracers'
          tracer_ddep_longnames(n) = &
                  trim(tracer_longnames(n)) // ' dry deposition for tracers'
      endif

!     id_tracer_ddep(n) = register_diag_field ( mod_name,               &
!                     trim(tracer_ddep_names(n)), mass_axes(1:2), Time, &
!                     trim(tracer_ddep_longnames(n)),                   &
!                     'kg/m2/s', missing_value=-999.     )
!     id_tracer_wdep(n) = register_diag_field ( mod_name,               &
!                     trim(tracer_wdep_names(n)), mass_axes(1:2), Time, &
!                     trim(tracer_wdep_longnames(n)),                   &
!                     'kg/m2/s', missing_value=-999.    )
   enddo

 
      call write_version_number (version, tagname)

         call write_namelist_values (stdlog(),ntrace)

      module_is_initialized = .TRUE.

end subroutine atmos_tracer_utilities_init
!</SUBROUTINE>
!
!#######################################################################
!
subroutine write_namelist_values (unit, ntrace)
    integer, intent(in) :: unit, ntrace
    integer :: n

    write (unit,10)
    do n = 1, ntrace
       write (unit,11) trim(tracer_wdep_names(n)),     &
                       trim(tracer_wdep_longnames(n)), &
                       trim(tracer_wdep_units(n))
       write (unit,11) trim(tracer_ddep_names(n)),     &
                       trim(tracer_ddep_longnames(n)), &
                       trim(tracer_ddep_units(n))
    enddo

 10 format (' &TRACER_DIAGNOSTICS_NML', &
          /,'    TRACER:  names  longnames  (units)')
 11 format (a16,2x,a,2x,'(',a,')')

 end subroutine write_namelist_values

!
!#######################################################################
!
!<SUBROUTINE NAME = "dry_deposition">
subroutine dry_deposition( n, is, js, u, v, T, pwt, pfull, &
                           u_star, landmask, dsinku, tracer, Time)
!
!<OVERVIEW>
! Routine to calculate the fraction of tracer to be removed by dry 
! deposition.
!</OVERVIEW>
!<DESCRIPTION>
! There are two types of dry deposition coded.
!
! 1) Wind driven derived dry deposition velocity.
!
! 2) Fixed dry deposition velocity.
! 
! The theory behind the wind driven dry deposition velocity calculation
! assumes that the deposition can be modeled as a parallel resistance type 
! problem.
!
!  Total resistance to HNO3-type dry deposition, 
!<PRE>       R = Ra + Rb
!  resisa = aerodynamic resistance
!  resisb = surface resistance (laminar layer + uptake)
!         = 5/u*  [s/cm]        for neutral stability
!      Vd = 1/R
!</PRE>
! For the fixed dry deposition velocity, there is no change in the 
! deposition velocity but the variation of the depth of the surface 
! layer implies that there is variation in the amount deposited.
!
! To utilize this section of code add one of the following lines as 
! a method for the tracer of interest in the field table.
!<PRE>
! "dry_deposition","wind_driven","surfr=XXX"
!     where XXX is the total resistance defined above.
!
! "dry_deposition","fixed","land=XXX, sea=YYY"
!     where XXX is the dry deposition velocity (m/s) over land
!       and YYY is the dry deposition velocity (m/s) over sea.
!</PRE>
!</DESCRIPTION>
!<TEMPLATE>
! call dry_deposition( n, is, js, u, v, T, pwt, pfull, 
!                           u_star, landmask, dsinku, tracer, Time)
!</TEMPLATE>
!
!  <IN NAME="n" TYPE="integer">
!    The tracer number.
!  </IN>
!  <IN NAME="is, js" TYPE="integer">
!    Start indices for array (computational indices).
!  </IN>
!  <IN NAME="u" TYPE="real(r8)" DIM="(:,:)">
!    U wind field.
!  </IN>
!  <IN NAME="v" TYPE="real(r8)" DIM="(:,:)">
!    V wind field.
!  </IN>
!  <IN NAME="T" TYPE="real(r8)" DIM="(:,:)">
!    Temperature.
!  </IN>
!  <IN NAME="pwt" TYPE="real(r8)" DIM="(:,:)">
!     Pressure differential of half levels.
!  </IN>
!  <IN NAME="pfull" TYPE="real(r8)" DIM="(:,:)">
!     Full pressure levels.
!  </IN>
!  <IN NAME="u_star" TYPE="real(r8)" DIM="(:,:)">
!     Friction velocity.
!  </IN>
!  <IN NAME="landmask" TYPE="logical">
!     Land - sea mask.
!  </IN>
!
!  <OUT NAME="dsinku" TYPE="real(r8)" DIM="(:,:)">
!    The amount of tracer in the surface layer which is dry deposited per second.
!  </OUT>
!
integer, intent(in)                 :: n, is, js
real(r8), intent(in), dimension(:,:)    :: u, v, T, pwt, pfull, u_star, tracer
logical, intent(in), dimension(:,:) :: landmask
type(time_type), intent(in)         :: Time
real(r8), intent(out), dimension(:,:)   :: dsinku

real(r8),dimension(size(u,1),size(u,2)) :: hwindv,frictv,resisa,xxfm,dz
integer :: i,j, flagsr
real(r8)    :: land_dry_dep_vel, sea_dry_dep_vel, surfr
logical :: used, flag
character(len=80) :: name,control,scheme

! Default zero
dsinku = 0.0
flag = query_method ('dry_deposition',MODEL_ATMOS,n,name,control)

if (.not. flag) return

! delta z = dp/(rho * grav)
! delta z = RT/g*dp/p    pwt = dp/g
dz(:,:) = pwt(:,:)*rdgas*T(:,:)/pfull(:,:)

call get_drydep_param(name,control,scheme,land_dry_dep_vel,sea_dry_dep_vel)
if(lowercase(scheme)=='wind_driven') then
! Calculate horizontal wind velocity and aerodynamic resistance:
!   where xxfm=(u*/u) is drag coefficient, Ra=u/(u*^2), 
!   and  u*=sqrt(momentum flux)  is friction velocity.
!
!****  Compute dry sinks (loss frequency, need modification when 
!****    different vdep values are to be used for species)
        flagsr=parse(control,'surfr',surfr)
        if(flagsr == 0) surfr=500.
        hwindv=sqrt(u**2+v**2)
        frictv=u_star
        resisa=hwindv/(u_star*u_star)
        where (frictv .lt. 0.1) frictv=0.1
        dsinku=(1./(surfr/frictv + resisa))/dz
!
else if (lowercase(scheme)=='fixed') then
!
! For the moment let's try to calculate the delta-z of the bottom 
! layer and using a simple dry deposition velocity times the 
! timestep, idt, calculate the fraction of the lowest layer which 
! deposits.
       where (landmask(:,:))
! dry dep value over the land surface divided by the height of the box.
         dsinku(:,:) = land_dry_dep_vel /dz(:,:)
      elsewhere
! dry dep value over the sea surface divided by the height of the box.
         dsinku(:,:) = sea_dry_dep_vel/dz(:,:)
      endwhere
endif

dsinku(:,:) = MAX(dsinku(:,:), 0.0_r8)
where(tracer>0)
  dsinku=dsinku*tracer
elsewhere
  dsinku=0.0
endwhere

! Now save the dry deposition to the diagnostic manager
! delta z = dp/(rho * grav)
! delta z *rho  = dp/g
! tracer(kgtracer/kgair) * dz(m)* rho(kgair/m3) = kgtracer/m2
! so rho drops out of the equation
         if (id_tracer_ddep(n) > 0 ) then
!         used = send_data ( id_tracer_ddep(n), dsinku*pwt, Time, &
!         is_in =is,js_in=js)
endif
end subroutine dry_deposition
!</SUBROUTINE>
!
!#######################################################################
!
!<SUBROUTINE NAME = "wet_deposition">
!<TEMPLATE>
!CALL wet_deposition(n, T, pfull, phalf, rain, snow, qdt, tracer, tracer_dt, Time, cloud_param, is, js)
!</TEMPLATE>
subroutine wet_deposition(n, T, pfull, phalf, rain, snow, qdt, tracer, tracer_dt, Time, cloud_param, is, js)
!
!<OVERVIEW>
! Routine to calculate the fraction of tracer removed by wet deposition
!</OVERVIEW>
!
!<IN NAME="n" TYPE="integer">
!   Tracer number
!</IN>
!<IN NAME="is, js" TYPE="integer">
!   start indices for array (computational indices)
!</IN>
!<IN NAME="T" TYPE="real(r8)" DIM="(:,:,:)">
!   Temperature
!</IN>
!<IN NAME="pfull" TYPE="real(r8)" DIM="(:,:,:)">
!   Full level pressure field
!</IN>
!<IN NAME="phalf" TYPE="real(r8)" DIM="(:,:,:)">
!   Half level pressure field
!</IN>
!<IN NAME="rain" TYPE="real(r8)" DIM="(:,:)">
!   Precipitation in the form of rain
!</IN>
!<IN NAME="snow" TYPE="real(r8)" DIM="(:,:)">
!   Precipitation in the form of snow
!</IN>
!<IN NAME="qdt" TYPE="real(r8)" DIM="(:,:,:)">
!   The tendency of the specific humidity due to the cloud parametrization
!</IN>
!<IN NAME="tracer" TYPE="real(r8)" DIM="(:,:,:)">
!   The tracer field 
!</IN>
!<IN NAME="Time" TYPE="type(time_type)">
!   The time structure for submitting wet deposition as a diagnostic
!</IN>
!<IN NAME="cloud_param" TYPE="character">
!   Is this a convective (convect) or large scale (lscale) cloud parametrization?
!</IN>
!  <OUT NAME="tracer_dt" TYPE="real(r8)" DIM="(:,:,:)">
!  The tendency of the tracer field due to wet deposition.
! </OUT>
!<DESCRIPTION>
! Schemes allowed here are 
!
! 1) Deposition removed in the same fractional amount as the modeled precipitation rate is to 
!    a standardized precipitation rate.
!    Basically this scheme assumes that a fractional area of the gridbox is affected by 
!    precipitation and that this precipitation rate is due to a cloud of standardized cloud 
!    liquid water content. Removal is constant throughout the column where precipitation is occuring.
!
! 2) Removal according to Henry's Law. This law states that the ratio of the concentation in 
!    cloud water and the partial pressure in the interstitial air is a constant. In this 
!    instance, the units for Henry's constant are kg/L/Pa (normally it is M/L/Pa)
!    Parameters for a large number of species can be found at
!    http://www.mpch-mainz.mpg.de/~sander/res/henry.html

! To utilize this section of code add one of the following lines as 
! a method for the tracer of interest in the field table.
!<PRE>
! "wet_deposition","henry","henry=XXX, dependence=YYY"
!     where XXX is the Henry's constant for the tracer in question
!       and YYY is the temperature dependence of the Henry's Law constant.
!
! "wet_deposition","fraction","lslwc=XXX, convlwc=YYY"
!     where XXX is the liquid water content of a standard large scale cloud
!       and YYY is the liquid water content of a standard convective cloud.
!</PRE>

!</DESCRIPTION>
!
integer, intent(in)                 :: n, is, js
real(r8), intent(in), dimension(:,:,:)  :: T, pfull,phalf, qdt, tracer
real(r8), intent(in), dimension(:,:)    :: rain, snow
character(len=*),intent(in)         :: cloud_param
type (time_type)  , intent(in)      :: Time
real(r8), intent(out), dimension(:,:,:) :: tracer_dt
!
real(r8), dimension(size(T,1),size(T,2),size(pfull,3))   :: wsinku
real(r8), dimension(size(T,1),size(T,2)) :: Htemp, dz, washout,scav_factor, sum_wdep
integer, dimension(size(T,1),size(T,2)) :: ktopcd, kendcd
integer :: i,j,k,kd, flaglw
real(r8)    :: Henry_constant, Henry_variable, inv298p15, clwc, wash, premin, prenow, hwtop
!real(r8), dimension(size(rain,1),size(rain,2)) :: prenow,hwtop
logical :: used,flag
character(len=80) :: name,control,scheme
tracer_dt = 0.0E+00
ktopcd = 0
kendcd = 0

flag = query_method ('wet_deposition',MODEL_ATMOS,n,name,control)
if(.not. flag) return
call get_wetdep_param(name,control,scheme,Henry_constant,Henry_variable)
if(lowercase(scheme)=='henry') then
! Henry_constant = [X](aq) / Px(g) 
! where [X](aq) is the concentration of tracer X in precipitation
!       Px(g) is the partial pressure of the tracer in the air
! [X](aq) = Mixing ratio (MR) in cloud / qdt
! Px(g)   = MR (non cloud) * Pfull 
!
! [X](aq)/Px = MR(incloud)/qdt /(Pfull MR non cloud) = H
! => MR(in cloud) = H * qdt * Pfull* MR(non cloud) 
! MR (total) = MR(incloud) + MR(noncloud)
!            = MR(noncloud) * ( 1 + H*Pfull*qdt)
! MR(incloud) = H*Pfull*qdt * MR(total)/(1+H*Pfull*qdt)
! Fraction removed = MR(incloud)/MR(total) =
!  H*Pfull*qdt/(1+H*Pfull*qdt)
!

 if(Henry_constant > 0 ) then
  inv298p15 = 1/298.15
  kd = size(T,3)
  do k = 1, kd
  ! Calculate the temperature dependent part of Henry's constant
  ! exp( k *(1/T - 1/298.15))
   Htemp(:,:) = exp(Henry_variable*(1/T(:,:,k)-inv298p15))
   tracer_dt(:,:,k) = 0.0
   scav_factor(:,:) = 0.0
   where (qdt(:,:,k) < 0.0)
   !qdt is -ve so need to multiply by -1.0
    scav_factor(:,:) = -1.0*Henry_constant*Htemp*pfull(:,:,k)*qdt(:,:,k)
    tracer_dt(:,:,k) = scav_factor(:,:)/(1+scav_factor(:,:))
   endwhere
  enddo
 endif 
endif



if(lowercase(scheme)=='fraction') then
tracer_dt = 0.0
!-----------------------------------------------------------------------
!
!     Compute areal(r8) fractions experiencing wet deposition:
!
!     Set minimum precipitation rate below which no wet removal
!     occurs to 0.01 cm/day ie 1.16e-6 mm/sec (kg/m2/s)
         premin=1.16e-6
!
!     Large scale cloud liquid water content (kg/m3)
!     and below cloud washout efficiency (cm-1):
            flaglw =parse(control,'lslwc',clwc)
            if (flaglw == 0 ) clwc=0.5e-3
            wash=1.0  
!
!     When convective adjustment occurs, use convective cloud liquid water content:
!
            if(trim(cloud_param) .eq. 'convect') then
              flaglw = parse(control,'convlwc',clwc)
              if (flaglw == 0) clwc=2.0e-3
              wash=0.3 
            end if
!
      do j=1,size(rain,2)
        do i=1,size(rain,1)
          tracer_dt(i,j,:)=0.0
          washout(i,j)=0.0
          prenow = rain(i,j) + snow(i,j)
          if(prenow .gt. premin) then      
!
! Assume that the top of the cloud is where the highest model level 
! specific humidity is reduced. And the the bottom of the cloud is the
! lowest model level where specific humidity is reduced.
!
            ktopcd(i,j) = 0
            do k = size(t,3),1,-1
             if (qdt(i,j,k) < 0.0 ) ktopcd(i,j) = k
            enddo
            kendcd(i,j) = 0
            do k = 1,size(t,3)
             if (qdt(i,j,k) < 0.0 ) kendcd(i,j) = k
            enddo
!
!     Thickness of precipitating cloud deck:
!
            if(ktopcd(i,j).gt.1) then
            hwtop = 0.0
            do k=ktopcd(i,j),kendcd(i,j)
             hwtop=hwtop+(phalf(i,j,k+1)-phalf(i,j,k))*rdgas*T(i,j,k)/grav/pfull(i,j,k)
            enddo
            do k=ktopcd(i,j),kendcd(i,j)
!     Areal(r8) fraction affected by precip clouds (max = 0.5):
             tracer_dt(i,j,k)=prenow/(clwc*hwtop)
            end do  
            endif

            washout(i,j)=prenow*wash
          end if
        end do
      end do
endif

! Now multiply by the tracer mixing ratio to get the actual tendency.
where (tracer_dt(:,:,:) .gt. 0.5) tracer_dt(:,:,:)=0.5
tracer_dt(:,:,:) = MAX(tracer_dt(:,:,:), 0.0_r8)
where(tracer>0)
tracer_dt=tracer_dt*tracer
elsewhere
tracer_dt=0.0
endwhere

sum_wdep=0.0
do k=1,size(tracer_dt,3)
! delta z = dp/(rho * grav)
! delta z = RT/g*dp/p
! tracer(kgtracer/kgair) * dz(m)* rho(kgair/m3) = kgtracer/m2
! so rho drops out of the equation
sum_wdep=sum_wdep + tracer_dt(:,:,k)*(phalf(:,:,k+1)-phalf(:,:,k))/grav
enddo

         if (id_tracer_wdep(n) > 0 ) then
!        used = send_data ( id_tracer_wdep(n), sum_wdep, Time , &
!         is_in =is,js_in=js)
endif
end subroutine wet_deposition
!</SUBROUTINE>
!
!#######################################################################
!
subroutine get_drydep_param(text_in_scheme,text_in_param,scheme,land_dry_dep_vel,sea_dry_dep_vel)
!
! Subroutine to initialiize the parameters for the dry deposition scheme.
! If the dry dep scheme is 'fixed' then the dry_deposition velocity value
! has to be set.
! If the dry dep scheme is 'wind_driven' then the dry_deposition
! velocity value will be calculated. So set to a dummy value of 0.0
! INTENT IN
!  text_in_scheme   : The text that has been parsed from tracer table as 
!                     the dry deposition scheme to be used.
!  text_in_param    : The parameters that are associated with the dry 
!                     deposition scheme.
! INTENT OUT
!  scheme           : The scheme that is being used.
!  land_dry_dep_vel : Dry deposition velocity over the land
!  sea_dry_dep_vel  : Dry deposition velocity over the sea
!
character(len=*), intent(in)    :: text_in_scheme, text_in_param
character(len=*), intent(out)   :: scheme
real(r8), intent(out)               :: land_dry_dep_vel, sea_dry_dep_vel

integer :: m,m1,n,lentext, flag
character(len=32) :: dummy

!Default
scheme                  = 'None'
land_dry_dep_vel=0.0
sea_dry_dep_vel=0.0

if(lowercase(trim(text_in_scheme(1:4))).eq.'wind') then
scheme                  = 'Wind_driven'
land_dry_dep_vel=0.0
sea_dry_dep_vel=0.0
endif

if(lowercase(trim(text_in_scheme(1:5))).eq.'fixed') then
scheme                 = 'fixed'
flag=parse(text_in_param,'land',land_dry_dep_vel)
flag=parse(text_in_param,'sea', sea_dry_dep_vel)
endif

end subroutine get_drydep_param
!
!#######################################################################
!
subroutine get_wetdep_param(text_in_scheme,text_in_param,scheme,henry_constant,henry_temp)
! 
! Routine to initialize the parameters for the wet deposition scheme.
! INTENT IN
!  text_in_scheme : Text read from the tracer table which provides information on which 
!                   wet deposition scheme to use.
!  text_in_param  : Parameters associated with the wet deposition scheme. These will be 
!                   parsed in this routine.
! INTENT OUT 
!  scheme         : Wet deposition scheme to use. 
!                   Choices are None, Fraction and Henry
!  henry_constant : Henry's constant for the tracer (see wet_deposition for explanation of Henry's Law)
!  henry_temp     : The temperature dependence of the Henry's Law constant.
!
!
character(len=*), intent(in)    :: text_in_scheme, text_in_param
character(len=*), intent(out)   :: scheme
real(r8), intent(out)               :: henry_constant, henry_temp

integer :: m,m1,n,lentext, flag
character(len=32) :: dummy

!Default
scheme                  = 'None'
henry_constant=0.0
henry_temp=0.0

if(trim(lowercase(text_in_scheme(1:8))).eq.'fraction') then
scheme                 = 'Fraction'
henry_constant=0.0
henry_temp=0.0
endif

if(trim(lowercase(text_in_scheme(1:5))).eq.'henry') then
scheme                 = 'Henry'
flag=parse(text_in_param,'henry',     henry_constant)
flag=parse(text_in_param,'dependence',henry_temp    )
endif

end subroutine get_wetdep_param
!
!#######################################################################
!
!<SUBROUTINE NAME="interp_emiss">
subroutine interp_emiss(global_source, start_lon, start_lat, &
                        lon_resol, lat_resol, data_out)
!
!<OVERVIEW>
! A routine to interpolate emission fields of arbitrary resolution onto the 
! resolution of the model.
!</OVERVIEW>
!<DESCRIPTION>
! Routine to interpolate emission fields (or any 2D field) to the model 
! resolution. The local section of the global field is returned to the 
! local processor.
!</DESCRIPTION>
! 
!<TEMPLATE>
! call interp_emiss(global_source, start_lon, start_lat, &
!                        lon_resol, lat_resol, data_out)
!</TEMPLATE>
! INTENT IN
!<IN NAME="global_source" TYPE="real(r8)" DIM="(:,:)">
!  Global emission field.
!</IN>
!<IN NAME="start_lon" TYPE="real(r8)">
!  Longitude of starting point of emission field 
!  (in radians). This is the westernmost boundary of the 
!  global field.
!</IN>
!<IN NAME="start_lat" TYPE="real(r8)">
!  Latitude of starting point of emission field
!  (in radians). This is the southern boundary of the 
!  global field.
!</IN>
!<IN NAME="lon_resol" TYPE="real(r8)">
!  Longitudinal resolution of the emission data (in radians).
!</IN>
!<IN NAME="lat_resol" TYPE="real(r8)">
!  Latitudinal resolution of the emission data (in radians).
!</IN>
! 
! INTENT OUT
!<OUT NAME="data_out" TYPE="real(r8)" DIM="(:,:)">
!  Interpolated emission field on the local PE. 
!</OUT>

real(r8), intent(in)  :: global_source(:,:)
real(r8), intent(in)  :: start_lon,start_lat,lon_resol,lat_resol
real(r8), intent(out) :: data_out(:,:)

real(r8) :: modydeg,modxdeg, tpi
integer :: i, j, nlon_in, nlat_in
real(r8) :: blon_in(size(global_source,1)+1)
real(r8) :: blat_in(size(global_source,2)+1)
! Set up the global surface boundary condition longitude-latitude boundary values

   tpi = 2. *PI
   nlon_in = size(global_source,1)
   nlat_in = size(global_source,2)
! For some reason the input longitude needs to be incremented by 180 degrees.
   do i = 1, nlon_in+1
      blon_in(i) = start_lon + float(i-1)*lon_resol + PI
   enddo
      if (abs(blon_in(nlon_in+1)-blon_in(1)-tpi) < epsilon(blon_in)) &
              blon_in(nlon_in+1)=blon_in(1)+tpi

   do j = 2, nlat_in
      blat_in(j) = start_lat + float(j-1)*lat_resol
   enddo
      blat_in(1)         = -0.5*PI
      blat_in(nlat_in+1) =  0.5*PI

! Now interpolate the global data to the model resolution
   call horiz_interp (global_source, blon_in, blat_in,    &
                        blon_out, blat_out, data_out )


end subroutine interp_emiss
!</SUBROUTINE>
!
!######################################################################
!<SUBROUTINE NAME="tracer_utilities_end">
!<OVERVIEW>
!  The destructor routine for the tracer utilities module.
!</OVERVIEW>
! <DESCRIPTION>
! This subroutine writes the version name to logfile and exits. 
! </DESCRIPTION>

subroutine atmos_tracer_utilities_end
 

   deallocate(blon_out, blat_out)
   module_is_initialized = .FALSE.

 end subroutine atmos_tracer_utilities_end
!</SUBROUTINE>

! ######################################################################
!

end module atmos_tracer_utilities_mod



