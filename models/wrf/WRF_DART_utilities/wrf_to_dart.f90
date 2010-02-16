! DART software - Copyright © 2004 - 2010 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
 
PROGRAM wrf_to_dart

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

use        types_mod, only : r8
use time_manager_mod, only : time_type, read_time, set_date, print_time
use    utilities_mod, only : get_unit, file_exist, open_file, close_file, &
                             error_handler, E_ERR, E_MSG, initialize_utilities, &
                             register_module, logfileunit, nmlfileunit, &
                             nc_check, finalize_utilities
use  assim_model_mod, only : static_init_assim_model, open_restart_write, &
                             awrite_state_restart, close_restart
use        model_mod, only : max_state_variables, num_state_table_columns, &
                             get_wrf_state_variables, get_number_domains, &
                             get_model_size, get_wrf_static_data, trans_3Dto1D, &
                             trans_2Dto1D, get_wrf_date, wrf_static_data_for_dart

use                          netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

logical, parameter :: debug = .false.

character(len=129) :: wrf_state_variables(num_state_table_columns,max_state_variables)

type(wrf_static_data_for_dart) :: wrf

real(r8), pointer :: dart(:)
real(r8), allocatable :: wrf_var_3d(:,:,:), wrf_var_2d(:,:)
type(time_type)   :: dart_time(2)
integer           :: number_dart_values, num_domains, &
                     year, month, day, hour, minute, second
integer           :: ndims, idims(2), dimids(2)
integer           :: i, ivtype, ind, dart_ind, my_index
character(len=80) :: varname
character(len=19) :: timestring
character(len=2)  :: idom

integer            :: ncid(50), var_id, id, iunit, dart_unit

write(*,*) 'WRF to DART'

call initialize_utilities('wrf_to_dart')
call register_module(source, revision, revdate)

call static_init_assim_model()
call get_wrf_state_variables(wrf_state_variables)

num_domains        = get_number_domains()
number_dart_values = get_model_size()

call error_handler(E_MSG,'wrf_to_dart', &
   'Converting a WRF netcdf file to a dart state vector', &
   source, revision, revdate)

if ( debug ) then
   print*,'dart vector size ',number_dart_values
endif

! allocate dart state vector
allocate(dart(number_dart_values))

! loop through domains again and append each variable to the state
dart_ind = 1
WRFDomains2 : do id = 1,num_domains

   wrf = get_wrf_static_data(id)
   write(idom,'(i2.2)') id

   call nc_check( nf90_open('wrfinput_d' // idom, NF90_NOWRITE, ncid(id)), &
                  'wrf_to_dart', 'open wrfinput_d' // idom )

   do ind = 1,wrf%number_of_wrf_variables

      ! actual location in state variable table
      my_index = wrf%var_index_list(ind)

      ! get stagger and variable size
      call nc_check( nf90_inq_varid(ncid(id),wrf_state_variables(1,my_index), &
                     var_id), 'wrf_to_dart', &
                     'inq_var_id '//wrf_state_variables(1,my_index))

      if (  wrf%var_size(3,ind) == 1 ) then

         allocate(wrf_var_2d(wrf%var_size(1,ind),wrf%var_size(2,ind)))

         call nc_check( nf90_get_var(ncid(id), var_id, wrf_var_2d), &
                     'wrf_to_dart','get_var '//wrf_state_variables(1,my_index) )

         call trans_2Dto1D( dart(dart_ind:),wrf_var_2d, &
                            wrf%var_size(1,ind),wrf%var_size(2,ind))

         deallocate(wrf_var_2d)

      else

         allocate(wrf_var_3d(wrf%var_size(1,ind),wrf%var_size(2,ind),wrf%var_size(3,ind)))

         call nc_check( nf90_get_var(ncid(id), var_id, wrf_var_3d), &
                     'wrf_to_dart','get_var '//wrf_state_variables(1,my_index) )

         call trans_3Dto1D( dart(dart_ind:),wrf_var_3d, wrf%var_size(1,ind), & 
                            wrf%var_size(2,ind),wrf%var_size(3,ind))

         deallocate(wrf_var_3d)

      endif 

      dart_ind = dart_ind + wrf%var_size(1,ind) * wrf%var_size(2,ind) * wrf%var_size(3,ind)
 
   enddo

enddo WRFDomains2

!---
!  output

if(debug) write(*,*) ' state output '
call nc_check( nf90_inq_varid(ncid(1), "Times", var_id), 'wrf_to_dart', &
               'inq_varid Times' )
call nc_check( nf90_inquire_variable(ncid(1), var_id, varname, xtype=ivtype, &
               ndims=ndims, dimids=dimids), 'wrf_to_dart', &
               'inquire_variable Times' )
do i=1,ndims
   call nc_check( nf90_inquire_dimension(ncid(1), dimids(i), &
                   len=idims(i)),'wrf_to_dart','inquire_dimensions Times' )
   if(debug) write(*,*) ' dimension ',i,idims(i)
enddo

call nc_check( nf90_get_var(ncid(1), var_id, timestring, &
               start = (/ 1, idims(2) /)), 'wrf_to_dart','get_var Times' )

call get_wrf_date(timestring, year, month, day, hour, minute, second)
dart_time(1) = set_date(year, month, day, hour, minute, second)

call print_time(dart_time(1),str='Time from wrfinput_d0x:')

iunit = get_unit()
if(file_exist('wrf.info')) then
   open(unit = iunit, file = 'wrf.info')
   dart_time(1) = read_time(iunit)
   close(iunit)
endif

call print_time(dart_time(1),str='Time written to dart vector file:')

!  open DART data file

dart_unit = open_restart_write("dart_wrf_vector")

call awrite_state_restart(dart_time(1), dart, dart_unit)

if(debug) write(*,*) ' returned from state output '

call close_restart(dart_unit)

do id=1,num_domains
   call nc_check ( nf90_sync(ncid(id)),'wrf_to_dart','sync wrfinput' )
   call nc_check ( nf90_close(ncid(id)),'wrf_to_dart','close wrfinput' )
enddo

deallocate(dart)

write(logfileunit,*)'FINISHED wrf_to_dart.'
write(logfileunit,*)

call finalize_utilities('wrf_to_dart')  ! closes log file.
 
end PROGRAM wrf_to_dart
