! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id$

program dart_to_noah

!----------------------------------------------------------------------
! purpose: redistribute changes to the diagnostic snow water equivalent
!          variable into the prognostic variables needed by noah.
!
! method: Read DART output and the noah restart file.
!         Calculate the increment and redistribute the implied changes
!         so that the noah restart file is consistent with the new
!         snow water equivalent.
!----------------------------------------------------------------------

use        types_mod, only : r8

use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             logfileunit, open_file, close_file, &
                             error_handler, E_MSG

use time_manager_mod, only : time_type, print_time, print_date, operator(-), get_time

use        model_mod, only : static_init_model, read_model_time

use  netcdf_utilities_mod, only : nc_check, nc_synchronize_file, &
                                  nc_get_variable, &
                                  nc_put_variable, &
                                  nc_get_variable_size, &
                                  nc_get_variable_num_dimensions, &
                                  nc_get_global_attribute, &
                                  nc_open_file_readonly, &
                                  nc_open_file_readwrite, &
                                  nc_close_file

use netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"
character(len=*), parameter :: routine  = "dart_to_noah"

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

character(len=256) :: dart_analysis_file = 'dart_posterior.nc'
character(len=256) :: noah_restart_file = 'noah_restart.nc'

namelist /dart_to_noah_nml/ dart_analysis_file, &
                           noah_restart_file

!----------------------------------------------------------------------

integer         :: iunit, io,i,j,k
type(time_type) :: dart_time, noah_time

integer :: ncid_dart, ncid_noah
integer :: numdims
integer :: dimlens(NF90_MAX_VAR_DIMS)

real(r8), allocatable :: innov_swe(:,:)
real(r8), allocatable :: innov(:,:)

real(r8), allocatable :: dart_sneqv(:,:)
real(r8), allocatable :: noah_sneqv(:,:)
real(r8), allocatable :: wt_swe(:,:)

real(r8), allocatable :: prior_snowh(:,:)
real(r8), allocatable :: poste_snowh(:,:)

real(r8), allocatable :: prior_snice(:,:,:)
real(r8), allocatable :: poste_snice(:,:,:)
real(r8), allocatable :: wt_ice(:,:)
real(r8), allocatable :: innov_ice(:,:)
real(r8), allocatable :: wt_ice_l(:,:,:)

real(r8), allocatable :: prior_snliq(:,:,:)
real(r8), allocatable :: poste_snliq(:,:,:)
real(r8), allocatable :: wt_liq(:,:)
real(r8), allocatable :: innov_liq(:,:)
real(r8), allocatable :: wt_liq_l(:,:,:)

real(r8), allocatable :: snow_density(:,:)
real(r8), allocatable :: sum_sn(:,:)
real(r8), allocatable :: sum_ice(:,:)
real(r8), allocatable :: sum_liq(:,:)

real(r8) :: missing_value

character(len=*), parameter :: SNEQV = 'SNEQV'
character(len=*), parameter :: SNOWH = 'SNOWH'
character(len=*), parameter :: SNICE = 'SNICE'
character(len=*), parameter :: SNLIQ = 'SNLIQ'

character(len=512) :: string1, string2, string3

!----------------------------------------------------------------------

call initialize_utilities(progname=routine)

!----------------------------------------------------------------------
! Call model_mod:static_init_model() which reads the noah namelists
! to set grid sizes, etc.
!>@todo if we want to manually check the date strings instead of
!> using static_init_model() and read_model_time(), we can greatly 
!> reduce the number of source code files needed.  
!----------------------------------------------------------------------

call static_init_model()

! Read the namelist to get the input filename. 

call find_namelist_in_file("input.nml", "dart_to_noah_nml", iunit)
read(iunit, nml = dart_to_noah_nml, iostat = io)
call check_namelist_read(iunit, io, "dart_to_noah_nml")

write(string1,*)'converting'
write(string2,*)'   DART file "'//trim(dart_analysis_file)//'"'
write(string3,*)'to noah file "'//trim(noah_restart_file)//'"'
call error_handler(E_MSG,routine,string1,text2=string2,text3=string3)

!----------------------------------------------------------------------
! Reads the valid time, the state, and the target time.
!----------------------------------------------------------------------

dart_time = read_model_time(dart_analysis_file)
noah_time = read_model_time(noah_restart_file)
!>@todo make sure the times match!

call print_date( dart_time,'dart_to_noah:noah model date')
call print_time( dart_time,'dart_to_noah:DART model time')
call print_date( dart_time,'dart_to_noah:noah model date',logfileunit)
call print_time( dart_time,'dart_to_noah:DART model time',logfileunit)

!----------------------------------------------------------------------
! update the current noah state vector
! SWE is diagnostic  and marked 'NOUPDATE' in the model_nml
! SNICE, SNLIQ, SNOWH are prognostic but are not in the model_nml AT ALL.
!----------------------------------------------------------------------

ncid_dart = nc_open_file_readonly(dart_analysis_file, routine)
ncid_noah = nc_open_file_readwrite(noah_restart_file, routine)


call nc_get_global_attribute(ncid_noah,'missing_value',missing_value,'dart_to_noah')

!>@todo check the number of dimensions to make sure they are what we expect

call nc_get_variable_num_dimensions(ncid_dart, SNEQV, numdims, routine)
call nc_get_variable_size(ncid_dart, SNEQV, dimlens(1:numdims), routine)

allocate(dart_sneqv(dimlens(1),dimlens(2)))
allocate(noah_sneqv(dimlens(1),dimlens(2)))
allocate(innov_swe(dimlens(1),dimlens(2)))
allocate(innov(dimlens(1),dimlens(2)))
allocate(wt_swe(dimlens(1),dimlens(2)))

call nc_get_variable(ncid_dart, SNEQV, dart_sneqv, routine)
call nc_get_variable(ncid_noah, SNEQV, noah_sneqv, routine)

innov_swe = dart_sneqv - noah_sneqv

call nc_close_file(ncid_dart,routine)
!----------------------------------------------------------------------
!>@todo check the number of dimensions to make sure they are what we expect
call nc_get_variable_num_dimensions(ncid_noah, SNOWH, numdims, routine)
call nc_get_variable_size(ncid_noah, SNOWH, dimlens(1:numdims), routine)
allocate(prior_snowh(dimlens(1),dimlens(2)))
allocate(poste_snowh(dimlens(1),dimlens(2)))
allocate(snow_density(dimlens(1),dimlens(2)))
allocate(sum_sn(dimlens(1),dimlens(2)))
call nc_get_variable(ncid_noah, SNOWH, prior_snowh, routine)

call nc_get_variable_num_dimensions(ncid_noah, SNICE, numdims, routine)
call nc_get_variable_size(ncid_noah, SNICE, dimlens(1:numdims), routine)
allocate(prior_snice(dimlens(1),dimlens(2),dimlens(3)))
allocate(poste_snice(dimlens(1),dimlens(2),dimlens(3)))
allocate(wt_ice(dimlens(1),dimlens(3)))
allocate(sum_ice(dimlens(1),dimlens(3)))
allocate(innov_ice(dimlens(1),dimlens(3)))
allocate(wt_ice_l(dimlens(1),dimlens(2),dimlens(3)))
call nc_get_variable(ncid_noah, SNICE, prior_snice, routine)

call nc_get_variable_num_dimensions(ncid_noah, SNLIQ, numdims, routine)
call nc_get_variable_size(ncid_noah, SNLIQ, dimlens(1:numdims), routine)
allocate(prior_snliq(dimlens(1),dimlens(2),dimlens(3)))
allocate(poste_snliq(dimlens(1),dimlens(2),dimlens(3)))
allocate(wt_liq(dimlens(1),dimlens(3)))
allocate(sum_liq(dimlens(1),dimlens(3)))
allocate(innov_liq(dimlens(1),dimlens(3)))
allocate(wt_liq_l(dimlens(1),dimlens(2),dimlens(3)))
call nc_get_variable(ncid_noah, SNLIQ, prior_snliq, routine)

!----------------------------------------------------------------------
! use the innov and the prior_* to come up with a poste_*
! assuming that the snow_density won't change
do i = 1,dimlens(1)
   do j = 1, dimlens(3)
      sum_sn(i,j) = 0.0_r8
      sum_ice(i,j)= 0.0_r8
      sum_liq(i,j)= 0.0_r8
      do k = 1, dimlens(2)
         if(prior_snice(i,k,j)>=0.0_r8 .and. prior_snliq(i,k,j)>=0.0_r8) then
            sum_sn(i,j) = sum_sn(i,j) + prior_snliq(i,k,j) + prior_snice(i,k,j)
            sum_ice(i,j)= sum_ice(i,j)+ prior_snice(i,k,j)
            sum_liq(i,j)= sum_liq(i,j)+ prior_snliq(i,k,j)
         endif
      enddo
      if(noah_sneqv(i,j)>0.0_r8) then
         wt_swe(i,j) = sum_sn(i,j)/noah_sneqv(i,j)
         innov(i,j) = innov_swe(i,j)*wt_swe(i,j)
      endif
      if( wt_swe(i,j)>0.0_r8) then
         wt_liq(i,j) = sum_liq(i,j)/sum_sn(i,j)
         wt_ice(i,j) = 1.0_r8-wt_liq(i,j) 
         innov_liq(i,j) = innov(i,j)*wt_liq(i,j)
         innov_ice(i,j) = innov(i,j)*wt_ice(i,j)
      endif
      do k = 1,dimlens(2)
         if(prior_snice(i,k,j)>=0.0_r8 .and. sum_ice(i,j)>0.0_r8 ) then
           wt_ice_l(i,k,j) = prior_snice(i,k,j)/sum_ice(i,j)
         endif
         if(prior_snliq(i,k,j)>=0.0_r8 .and. sum_liq(i,j)>0.0_r8 ) then
           wt_liq_l(i,k,j) = prior_snliq(i,k,j)/sum_liq(i,j)
         endif
      enddo   
         
   enddo
enddo

do i = 1,dimlens(1)
   do j = 1, dimlens(3)
      do k = 1, dimlens(2)
         if(prior_snice(i,k,j)>=0.0_r8) then
            poste_snice(i,k,j) = prior_snice(i,k,j)+innov_ice(i,j)*wt_ice_l(i,k,j)
         else
            poste_snice(i,k,j) = missing_value
         endif
          if(prior_snliq(i,k,j)>=0.0_r8) then
            poste_snliq(i,k,j) = prior_snliq(i,k,j)+innov_liq(i,j)*wt_liq_l(i,k,j)
         else
            poste_snliq(i,k,j) = missing_value
         endif
         if(poste_snice(i,k,j)<0.0_r8 .and. poste_snice(i,k,j)/=missing_value) then
            poste_snice(i,k,j) = prior_snice(i,k,j)
         endif
         if(poste_snliq(i,k,j)<0.0_r8 .and. poste_snliq(i,k,j)/=missing_value) then
            poste_snliq(i,k,j) = prior_snliq(i,k,j)
         endif
      enddo
    enddo
enddo            

!call calculate_posterior_snowh(innov,prior_snowh,poste_snowh)
Do i = 1,dimlens(1)
   Do j = 1,dimlens(3)
       if(prior_snowh(i,j)>0.0_r8 .and. sum_sn(i,j)>0.0_r8) then
          snow_density(i,j) = sum_sn(i,j)/prior_snowh(i,j)
          poste_snowh(i,j) = innov(i,j)/snow_density(i,j) + prior_snowh(i,j)
       else if (sum_sn(i,j)==0.0_r8) then
          poste_snowh(i,j) = prior_snowh(i,j)
       else
          poste_snowh(i,j) = missing_value
       endif
       if(poste_snowh(i,j)<0.0_r8 .and. poste_snowh(i,j)/=missing_value) then
          poste_snowh(i,j) = prior_snowh(i,j)
       endif
   enddo
enddo
!test---------------------------------------------------------
!print*, prior_snowh(40,60),poste_snowh(40,60)
!print*, noah_sneqv(40,60),dart_sneqv(40,60)
!print*, prior_snice(40,3,60),poste_snice(40,3,60)
!print*, prior_snliq(40,3,60),poste_snliq(40,3,60)
!-------------------------------------------------------------
call nc_put_variable(ncid_noah, SNOWH, poste_snowh, routine)
call nc_put_variable(ncid_noah, SNEQV, dart_sneqv, routine)
call nc_put_variable(ncid_noah, SNICE, poste_snice, routine)
call nc_put_variable(ncid_noah, SNLIQ, poste_snliq, routine)

call nc_close_file(ncid_noah,routine)

deallocate(prior_snowh, poste_snowh)
deallocate(prior_snice, poste_snice)
deallocate(prior_snliq, poste_snliq)
deallocate(dart_sneqv, noah_sneqv)
deallocate(wt_swe,wt_ice,wt_liq,wt_liq_l,wt_ice_l,innov,innov_swe,innov_ice,innov_liq)
deallocate(sum_sn,sum_ice,sum_liq,snow_density)

!----------------------------------------------------------------------

call finalize_utilities(routine)

!contains
!subroutine calculate_posterior_snowh(innovation, prior, posterior)
!real(r8), intent(in)  :: innovation(:,:)
!real(r8), intent(in)  :: prior(:,:)
!real(r8), intent(out) :: posterior(:,:)
!end subroutine calculate_posterior_snowh


end program dart_to_noah

