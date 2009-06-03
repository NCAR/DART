! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
PROGRAM dart_to_wrf

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

use        types_mod, only : r8, missing_r8
use time_manager_mod, only : time_type, write_time, read_time, get_date, set_date, operator(-), &
                             get_time, print_time, set_calendar_type, GREGORIAN, julian_day
use    utilities_mod, only : get_unit, file_exist, open_file, close_file, &
                             error_handler, E_ERR, E_MSG, initialize_utilities, &
                             register_module, logfileunit, nmlfileunit, timestamp, &
                             find_namelist_in_file, check_namelist_read, &
                             nc_check, do_nml_file, do_nml_term
use  assim_model_mod, only : open_restart_read, open_restart_write, aread_state_restart, &
                             awrite_state_restart
use model_mod, only        : max_state_variables,  &
                             num_state_table_columns, read_wrf_dimensions, &
                             num_bounds_table_columns, &
                             get_number_of_wrf_variables,  &
                             get_variable_bounds,          &
                             set_variable_bound_defaults,          &
                             get_variable_size_from_file, wrf_dom, &
                             fill_default_state_table, &
                             trans_1Dto3D, trans_1Dto2D,&
                             set_wrf_date

use                          netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"


!-----------------------------------------------------------------------
! Model namelist parameters with default values.
!-----------------------------------------------------------------------

logical :: output_state_vector  = .false.  ! state vs. prognostic format
logical :: default_state_variables = .true.   ! use default state list?
character(len=129) :: wrf_state_variables(num_state_table_columns,max_state_variables) = 'NULL'
character(len=129) :: wrf_state_bounds(num_bounds_table_columns,max_state_variables) = 'NULL'
integer :: num_moist_vars       = 3
integer :: num_domains          = 1
integer :: calendar_type        = GREGORIAN
integer :: assimilation_period_seconds = 21600
logical :: surf_obs             = .true.
logical :: soil_data            = .true.
logical :: h_diab               = .false.
logical :: allow_obs_below_vol  = .false.
character(len = 72) :: adv_mod_command = './wrf.exe'
real (kind=r8) :: center_search_half_length = 500000.0_r8
integer :: center_spline_grid_scale = 10
integer :: vert_localization_coord  =  3  ! 1,2,3 == level,pressure,height
! candidates for including in the WRF netcdf files:
logical :: polar = .false.         ! wrap over the poles
logical :: periodic_x = .false.    ! wrap in longitude or x
logical :: periodic_y = .false.    ! used for single column model, wrap in y
!JPH -- single column model flag 
logical :: scm        = .false.    ! using the single column model


namelist /model_nml/ output_state_vector, num_moist_vars, &
                     num_domains, calendar_type, surf_obs, soil_data, h_diab, &
                     default_state_variables, wrf_state_variables, &
                     wrf_state_bounds, &
                     adv_mod_command, assimilation_period_seconds, &
                     allow_obs_below_vol, vert_localization_coord, &
                     center_search_half_length, center_spline_grid_scale, &
                     polar, periodic_x, periodic_y, scm


!-------------------------------------------------------------

type(wrf_dom) :: wrf

real(r8), pointer :: dart(:)
real(r8), pointer :: wrf_var_3d(:,:,:), wrf_var_2d(:,:)
real(r8)          :: lb, ub
type(time_type)   :: dart_time(2)
integer           :: number_dart_values, ndays, &
                     year, month, day, hour, minute, second
integer           :: ind, dart_ind, my_index
character(len=19) :: timestring
character(len=1)  :: idom

logical, parameter :: debug = .true.
integer            :: ncid(50), io, var_id, id, iunit, dart_unit
integer            :: var_element_list(max_state_variables)

write(*,*) 'DART to WRF' !(.false./F)?'

call initialize_utilities('dart_to_wrf')
call register_module(source, revision, revdate)

! Begin by reading the namelist input
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

call set_calendar_type(calendar_type)

call error_handler(E_MSG,'dart_to_wrf', &
   'Converting a dart_state_vector to a WRF netcdf file', &
   source, revision, revdate)

allocate(wrf%dom(num_domains))

! get default state variable table if asked
if ( default_state_variables ) then
  wrf_state_variables = 'NULL'
  call fill_default_state_table(wrf_state_variables)
endif

if ( debug ) then
  print*,'WRF state vector table'
  print*,'default_state_variables = ',default_state_variables
  print*,wrf_state_variables
endif


! open wrf data netCDF file 'wrfinput_d0x'
! we get sizes of the WRF geometry and resolution

! approach is to do everything in a loop, without any regard to the actual
! order.  just use the namelist to read and unroll one at a time.

! big loop over domains, just like in static_init
number_dart_values = 0
WRFDomains : do id = 1,num_domains

   write(idom,'(I1)') id

   if(file_exist('wrfinput_d0'//idom)) then

      call nc_check( nf90_open('wrfinput_d0'//idom, NF90_WRITE, ncid(id)), &
                      'dart_to_wrf','open wrfinput_d0'//idom )

   else

      call error_handler(E_ERR,'dart_to_wrf', &
           'Please put wrfinput_d0'//idom//' in the work directory.', source, revision,revdate)

   endif

! read WRF dimensions
   call read_wrf_dimensions(ncid(id),wrf%dom(id)%bt, wrf%dom(id)%bts, &
                                 wrf%dom(id)%sn, wrf%dom(id)%sns, &
                                 wrf%dom(id)%we, wrf%dom(id)%wes, &
                                 wrf%dom(id)%sls)


! get the number of wrf variables wanted in this domain's state
   wrf%dom(id)%number_of_wrf_variables = get_number_of_wrf_variables(id,wrf_state_variables,var_element_list)

  if ( debug ) then
    print*,'Domain ',id,' number of wrf variables is ',wrf%dom(id)%number_of_wrf_variables
  endif

! allocate and store the table locations of the variables valid on this domain
   allocate(wrf%dom(id)%var_index_list(wrf%dom(id)%number_of_wrf_variables))
   wrf%dom(id)%var_index_list = var_element_list(1:wrf%dom(id)%number_of_wrf_variables)

! allocate bounds lists and instructions
   allocate(wrf%dom(id)%lower_bound(wrf%dom(id)%number_of_wrf_variables))
   allocate(wrf%dom(id)%upper_bound(wrf%dom(id)%number_of_wrf_variables))
   allocate(wrf%dom(id)%clamp_or_fail(wrf%dom(id)%number_of_wrf_variables))
   call set_variable_bound_defaults(wrf%dom(id)%number_of_wrf_variables, &
                                    wrf%dom(id)%lower_bound, &
                                    wrf%dom(id)%upper_bound, &
                                    wrf%dom(id)%clamp_or_fail)

! assign default bounds and link the bounds to the correct variable locations
   do ind = 1,wrf%dom(id)%number_of_wrf_variables

      ! actual location in state variable table
      my_index =  wrf%dom(id)%var_index_list(ind)

      call get_variable_bounds(wrf_state_bounds, &
                               wrf_state_variables(1,my_index), &
                               wrf%dom(id)%lower_bound(ind), &
                               wrf%dom(id)%upper_bound(ind), &
                               wrf%dom(id)%clamp_or_fail(ind))

      if ( debug ) then
         write(*,*) 'Bounds for variable ',  &
                     trim(wrf_state_variables(1,my_index)), &
                     ' are ',wrf%dom(id)%lower_bound(ind), &
                     wrf%dom(id)%upper_bound(ind), &
                     wrf%dom(id)%clamp_or_fail(ind)
      endif

   enddo


! allocate var size
  allocate(wrf%dom(id)%var_size(3,wrf%dom(id)%number_of_wrf_variables))

! allocate stagger
  allocate(wrf%dom(id)%stagger(wrf%dom(id)%number_of_wrf_variables))

! add to the dart state vector size
   do ind = 1,wrf%dom(id)%number_of_wrf_variables

      ! actual location in state variable table
      my_index =  wrf%dom(id)%var_index_list(ind)

      ! get stagger and variable size
      call get_variable_size_from_file(ncid(id),id,  &
                                       wrf_state_variables(1,my_index), &
                                       wrf%dom(id)%bt, wrf%dom(id)%bts, &
                                       wrf%dom(id)%sn, wrf%dom(id)%sns, &
                                       wrf%dom(id)%we, wrf%dom(id)%wes, &
                                       wrf%dom(id)%stagger(ind),        &
                                       wrf%dom(id)%var_size(:,ind))

      if ( debug ) then
         print*,'variable size ',trim(wrf_state_variables(1,my_index)),' ',wrf%dom(id)%var_size(:,ind)
      endif

      number_dart_values = number_dart_values &
                         + (wrf%dom(id)%var_size(1,ind) &
                           * wrf%dom(id)%var_size(2,ind) &
                           * wrf%dom(id)%var_size(3,ind)) 
   enddo

enddo WRFDomains

if ( debug ) then
   print*,'dart vector size ',number_dart_values
endif

! allocate dart state vector
allocate(dart(number_dart_values))

! open dart file
dart_unit = open_restart_read("dart_wrf_vector")

! read dart vector
call aread_state_restart(dart_time(2), dart, dart_unit, dart_time(1))

! record wrf.info
iunit = get_unit()
open(unit = iunit, file = 'wrf.info')
call write_time(iunit, dart_time(1))
call write_time(iunit, dart_time(2))
call get_date(dart_time(2), year, month, day, hour, minute, second)
write (iunit,FMT='(I4,5I3.2)') year, month, day, hour, minute, second

write (iunit,*) num_domains
write (iunit,*) adv_mod_command
close(iunit)

! set new times
call set_wrf_date(timestring, year, month, day, hour, minute, second)
ndays = julian_day(year, month, day)

! loop through domains again and pull each variable from the state
dart_ind = 1
WRFDomains2 : do id = 1,num_domains

   do ind = 1,wrf%dom(id)%number_of_wrf_variables

      ! actual location in state variable table
      my_index = wrf%dom(id)%var_index_list(ind)

      if ( debug ) then
         write(*,*) 'Rolling up variable ',trim(wrf_state_variables(1,my_index))
      endif

      ! get stagger and variable size
      call nc_check( nf90_inq_varid(ncid(id),wrf_state_variables(1,my_index), &
                     var_id), 'dart_to_wrf', &
                     'inq_var_id '//wrf_state_variables(1,my_index))

      if (  wrf%dom(id)%var_size(3,ind) == 1 ) then

         if ( debug ) then
            write(*,*) trim(wrf_state_variables(1,my_index)),' is 2D'
            write(*,*) 'size  ',wrf%dom(id)%var_size(:,ind)
         endif

         allocate(wrf_var_2d(wrf%dom(id)%var_size(1,ind),wrf%dom(id)%var_size(2,ind)))

         wrf_var_2d = 0.0_r8

         call trans_1Dto2D( dart(dart_ind:),wrf_var_2d, &
                       wrf%dom(id)%var_size(1,ind),wrf%dom(id)%var_size(2,ind))

         call nc_check( nf90_put_var(ncid(id), var_id, wrf_var_2d), &
                     'dart_to_wrf','put_var '//wrf_state_variables(1,my_index) )


         deallocate(wrf_var_2d)

      else

         if ( debug ) then
            write(*,*) trim(wrf_state_variables(1,my_index)),' is 3D'
            write(*,*) 'size  ',wrf%dom(id)%var_size(:,ind)
         endif

         allocate(wrf_var_3d(wrf%dom(id)%var_size(1,ind),wrf%dom(id)%var_size(2,ind),wrf%dom(id)%var_size(3,ind)))

         wrf_var_3d = 0.0_r8

         call trans_1Dto3D( dart(dart_ind:),wrf_var_3d, &
                      wrf%dom(id)%var_size(1,ind), &
                      wrf%dom(id)%var_size(2,ind),wrf%dom(id)%var_size(3,ind))
 
         ! check for bounds
         lb = wrf%dom(id)%lower_bound(my_index)
         ub = wrf%dom(id)%upper_bound(my_index)

         ! check bounds and fail if requested
         if ( trim(wrf%dom(id)%clamp_or_fail(my_index)) == 'FAIL' ) then

            if (minval(wrf_var_3d) < lb ) then
               call error_handler(E_ERR,'dart_to_wrf', &
               'Variable '//trim(wrf_state_variables(1,my_index))// &
               ' failed lower bounds check.', source, revision,revdate)
            endif

            if (maxval(wrf_var_3d) > ub ) then
               call error_handler(E_ERR,'dart_to_wrf', &
               'Variable '//trim(wrf_state_variables(1,my_index))// &
               ' failed upper bounds check.', source, revision,revdate)
            endif

         endif ! bounds check failure request

         ! apply bounds if asked
         if ( lb /= missing_r8 ) then

            if ( debug ) then
              write(*,*) 'Setting lower bound ',lb,' on ', &
                          trim(wrf_state_variables(1,my_index))
            endif

            wrf_var_3d = max(lb,wrf_var_3d)
  
         endif

         if ( ub /= missing_r8 ) then

            if ( debug ) then
              write(*,*) 'Setting upper bound ',ub,' on ', &
                          trim(wrf_state_variables(1,my_index))
            endif

            wrf_var_3d = min(ub,wrf_var_3d)
  
         endif

!         if ( trim(wrf_state_variables(1,my_index)) == 'QVAPOR' .or.    &
!              trim(wrf_state_variables(1,my_index)) == 'QRAIN' .or.     &
!              trim(wrf_state_variables(1,my_index)) == 'QCLOUD' ) then
!
!            if ( debug ) then
!              write(*,*) 'Setting 0 lower bound on ', &
!                          trim(wrf_state_variables(1,my_index))
!            endif

!            wrf_var_3d = max(0.0_r8,wrf_var_3d)

!         endif

         call nc_check( nf90_put_var(ncid(id), var_id, wrf_var_3d), &
                     'dart_to_wrf','put_var '//wrf_state_variables(1,my_index) )

         deallocate(wrf_var_3d)

      endif 

      dart_ind = dart_ind &
                 + (wrf%dom(id)%var_size(1,ind) &
                  * wrf%dom(id)%var_size(2,ind) &
                  * wrf%dom(id)%var_size(3,ind))
 
   enddo

   ! output times too
   call nc_check( nf90_inq_varid(ncid(id), "Times", var_id), 'dart_to_wrf', &
               'inq_varid Times' )
   call nc_check( nf90_put_var(ncid(id), var_id, timestring), 'dart_to_wrf', &
               'put_var Times' )
   call nc_check( nf90_put_att(ncid(id), nf90_global, "START_DATE", timestring), &
                'dart_to_wrf','put_att START_DATE' )
   call nc_check( nf90_put_att(ncid(id), nf90_global, "JULYR", year), &
                'dart_to_wrf','put_att JULYR' )
   call nc_check( nf90_put_att(ncid(id), nf90_global, "JULDAY", ndays), &
                'dart_to_wrf','put_att JULDY' ) 

   ! at this point we are done with this domain
   call nc_check ( nf90_sync(ncid(id)),'dart_to_wrf','sync wrfinput' )
   call nc_check ( nf90_close(ncid(id)),'dart_to_wrf','close wrfinput' )

enddo WRFDomains2

deallocate(wrf%dom)
deallocate(dart)

write(logfileunit,*)'FINISHED dart_to_wrf.'
write(logfileunit,*)

call timestamp(source,revision,revdate,'end') ! That closes the log file, too.

end PROGRAM dart_to_wrf

