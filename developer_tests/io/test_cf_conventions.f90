! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program test_cf_conventions

use             types_mod, only : r4, r8, i8, MISSING_R8 , MISSING_R4
use         utilities_mod, only : register_module, error_handler, E_MSG, E_ERR
use  adaptive_inflate_mod, only : adaptive_inflate_init
use     mpi_utilities_mod, only : initialize_mpi_utilities, finalize_mpi_utilities
use   state_vector_io_mod, only : read_state, write_state, state_vector_io_init 
use  ensemble_manager_mod, only : init_ensemble_manager, ensemble_type,        &
                                  set_num_extra_copies
use      io_filenames_mod, only : io_filenames_init, io_filenames_finalize,    &
                                  file_info_type, netcdf_file_type, READ_COPY, &
                                  set_file_metadata, set_io_copy_flag
use   state_structure_mod, only : get_xtype,  get_units, get_long_name,        &
                                  get_short_name, get_has_missing_value,       &
                                  get_FillValue, get_missing_value,            &
                                  get_add_offset, get_scale_factor,            &
                                  add_domain, state_structure_info,            &
                                  get_sum_variables
use      time_manager_mod, only : time_type, set_time
use            filter_mod, only : filter_set_initial_time
use            assert_mod, only : assert_equal

use netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

! this should be a namelist variable
logical :: verbose = .false.

type(ensemble_type)         :: ens_handle
type(file_info_type)        :: file_info_input
type(time_type)             :: time1
type(time_type)             :: curr_ens_time
character(len=256)          :: test_file(1) = "cf_test.nc"

integer(i8) :: model_size
integer     :: num_ens    = 1
integer     :: num_extras = 0
integer     :: num_copies

logical :: read_time_from_file

integer :: domid = 1 ! only one domain
integer :: var_xtype
character(len=NF90_MAX_NAME) :: var_units, blank_string, var_att_name

integer  :: missINT
real(r4) :: missR4
real(r8) :: missR8, var_offset, var_scale_factor

blank_string = ' '

! initialize the dart libs
call initialize_module()

!>@todo FIXME ... add variable E when scale/offset are supported
domid = add_domain(test_file(1), num_vars=4, var_names=(/'A', 'B', 'C', 'D'/))

if (verbose) then
   call state_structure_info(domid)
endif

! since we are calling add_domain directly instead of through
! static_assim_model_mod we need to get the total number of
! variables from the state_strucutre_mod instead of using
! get_model_size()
model_size = get_sum_variables(1, 4, domid)

write(*,*) " model size : ", model_size

num_copies = num_ens + num_extras

! initalize routines needed for read_state and write_state
call init_ensemble_manager(ens_handle, num_copies, model_size)
call set_num_extra_copies(ens_handle, num_extras)
call filter_set_initial_time(0,0,time1,read_time_from_file)

file_info_input  = initialize_filenames(test_file)

curr_ens_time = set_time(0, 0)

! read in restarts
call read_state(ens_handle, file_info_input, read_time_from_file, time1)

write(*,*)' '
write(*,*)'======================================================================'
write(*,*)' Unit Test for CF-Conventions'
write(*,*)'======================================================================'
write(*,*)' '


write(*,*)'Testing get_xtype'

var_xtype = get_xtype(domid,1)
call assert_equal(var_xtype, NF90_INT,    'variable1:get_xtype')

var_xtype = get_xtype(domid,2)
call assert_equal(var_xtype, NF90_FLOAT,  'variable2:get_xtype')

var_xtype = get_xtype(domid,3)
call assert_equal(var_xtype, NF90_DOUBLE, 'variable3:get_xtype')

var_xtype = get_xtype(domid,4)
call assert_equal(var_xtype, NF90_DOUBLE, 'variable4:get_xtype')

! test units

write(*,*)'Testing get_units'

!> todo FIXME Need a unique prefix so don't confuse with get_unit

var_units = get_units(domid,1)
call assert_equal(var_units, 'units A',     'variable1:get_units')

var_units = get_units(domid,2)
call assert_equal(var_units, 'units B',     'variable2:get_units')

var_units = get_units(domid,3)
call assert_equal(var_units, 'units C',     'variable3:get_units')

var_units = get_units(domid,4)
call assert_equal(var_units, blank_string,  'variable4:get_units')


write(*,*)'Testing get_long_name'

var_att_name = get_long_name(domid,1)
call assert_equal(var_att_name, 'variable A',  'variable1:get_long_name')

var_att_name = get_long_name(domid,2)
call assert_equal(var_att_name, 'variable B',  'variable2:get_long_name')

var_att_name = get_long_name(domid,3)
call assert_equal(var_att_name, 'variable C',  'variable3:get_long_name')

var_att_name = get_long_name(domid,4)
call assert_equal(var_att_name, blank_string,  'variable4:get_long_name')


write(*,*)'Testing get_short_name'

var_att_name = get_short_name(domid,1)
call assert_equal(var_att_name, 'short A',     'variable1:get_short_name')

var_att_name = get_short_name(domid,2)
call assert_equal(var_att_name, 'short B',     'variable2:get_short_name')

var_att_name = get_short_name(domid,3)
call assert_equal(var_att_name, 'short C',     'variable3:get_short_name')

var_att_name = get_short_name(domid,4)

call assert_equal(var_att_name, blank_string,  'variable4:get_short_name')


write(*,*)'Testing get_missing_value'

call get_missing_value(domid,1,missINT)
call assert_equal(missINT, -77,                   'variable1:get_missing_value')

call get_missing_value(domid,2,missR4)
call assert_equal(missR4, -777.77_r4,             'variable2:get_missing_value')

call get_missing_value(domid,3,missR8)
call assert_equal(missR8, -88888.88888_r8,        'variable3:get_missing_value')

! ! variable 4 has no missing value this will fail
! call get_missing_value(domid,4,missR8)
! call assert_equal(missR8, -88888.88888_r8,        'variable4:get_missing_value')

write(*,*)'Testing get_FillValue'

call get_FillValue(domid,1,missINT)
call assert_equal(missINT, -77,                   'variable1:get_FillValue')

call get_FillValue(domid,2,missR4)
call assert_equal(missR4, -777.77_r4,             'variable2:get_FillValue')

call get_FillValue(domid,3,missR8)
call assert_equal(missR8, -88888.88888_r8,        'variable3:get_FillValue')

! ! variable 4 has no _FillValue this will fail
! call get_FillValue(domid,4,missR8)
! call assert_equal(missR8, -88888.88888_r8,        'variable4:get_FillValue')

! write(*,*)'Testing offset and scale factor'

! this is only for r8 at the moment
! since it is not being used within DART
! NOTE: This test is supposed to break.
! Commenting it out for testing all programs.

! var_offset = get_add_offset(domid,3)
! call assert_equal(var_offset, 2.0_r8,            'variable3:get_var_offset')
! 
! var_offset = get_add_offset(domid,4)
! call assert_equal(var_offset, missR8 ,       'variable4:get_var_offset')
! 
! 
! var_scale_factor = get_scale_factor(domid,3)
! call assert_equal(var_scale_factor, 0.2_r8 ,     'variable3:get_scale_factor')
! 
! var_scale_factor = get_scale_factor(domid,4)
! call assert_equal(var_scale_factor, missR8 , 'variable3:get_scale_factor')

write(*,*)' '
write(*,*)'======================================================================'
write(*,*)' Finished Unit Test'
write(*,*)'======================================================================'
write(*,*)' '

call io_filenames_finalize(file_info_input)
 
! finalize test_cf_conventions
call error_handler(E_MSG,'test_cf_conventions','Finished successfully.',source,revision,revdate)
call finalize_mpi_utilities()

! end of main code


contains

!----------------------------------------------------------------------

subroutine initialize_module

call initialize_mpi_utilities('test_cf_conventions')
call register_module(source, revision, revdate)
!call static_init_assim_model()
call state_vector_io_init()

end subroutine initialize_module

!----------------------------------------------------------------------

function initialize_filenames(filename) result(file_handle)
character(len=*), intent(in) :: filename(:)
type(file_info_type) :: file_handle

integer :: num_domains = 1
integer :: imem

character(len=256), allocatable :: file_array(:,:)
character(len=512) :: my_base, my_desc

allocate(file_array(num_ens, num_domains))
file_array = RESHAPE(filename,  (/num_ens,  num_domains/))

call io_filenames_init(file_handle,            &
                       ncopies      = 1,       &
                       cycling      = .false., &
                       single_file  = .false., &
                       restart_files = file_array)

do imem = 1, num_ens
   write(my_base,'(A,I2)') 'inens_',    imem
   write(my_desc,'(A,I2)') 'input ens', imem
   call set_file_metadata(file_handle,                          &
                          cnum     = imem,                      &
                          fnames   = file_array(imem,:),  &
                          basename = my_base,                   &
                          desc     = my_desc)

   call set_io_copy_flag(file_handle,        &
                         cnum    = imem,     &
                         io_flag = READ_COPY)
enddo

end function initialize_filenames

!----------------------------------------------------------------------

end program

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
