! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> start of test of setting and getting variable sizes from the state struct.

program test_sizes

use           types_mod, only : r8, i8, MISSING_R8

use       utilities_mod, only : initialize_utilities, finalize_utilities, &
                                error_handler, E_ERR, E_MSG

use model_mod,           only : static_init_model

use          assert_mod, only : assert_equal

use state_structure_mod, only : add_domain,                 &
                                get_domain_size,            &
                                get_num_domains,            &
                                get_variable_size,          &
                                get_variable_name,          &
                                get_num_variables,          &
                                get_num_dims,               &
                                get_dim_lengths,            &
                                get_dim_length,             &
                                add_dimension_to_variable,  &
                                finished_adding_domain,     &
                                state_structure_info,       &
                                get_kind_string,            &
                                get_kind_index,             &
                                get_varid_from_kind,        &
                                get_varids_from_kind,       &
                                get_dim_name,               &
                                get_io_num_dims,            &
                                get_io_dim_ids,             &
                                get_io_dim_lengths,         &
                                get_io_num_unique_dims,     &
                                get_io_unique_dim_name,     &
                                get_io_unique_dim_length,   &
                                add_time_unlimited,         &
                                get_unlimited_dimid,        &
                                set_var_id,                 &
                                get_io_clamping_maxval,     &
                                get_io_clamping_minval,     &
                                do_io_clamping,             &
                                do_io_update,               &
                                get_index_start,            &
                                get_index_end,              &
                                get_sum_variables,          &
                                get_sum_variables_below,    &
                                get_model_variable_indices, &
                                get_dart_vector_index,      &
                                get_num_varids_from_kind

! diagnostic files 
use state_structure_mod, only : create_diagnostic_structure, &
                                end_diagnostic_structure
 

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

integer, allocatable :: int_array(:)
integer  :: ds, did1, did2, did3
integer  :: num_domains, var_size
real(r8) :: realnum
logical  :: truefalse
integer  :: iloc, jloc, kloc, var_id, dom_id, kind_index 
integer(i8) :: bigint, index_in

character(len=512) :: kind_string

character(len=512) :: string1, string2, string3

! main code here
 
! initialize the dart libs
call initialize_utilities('test_sizes')

write(*,*)
write(*,*)'======================================================================'
write(*,*)'Testing the add_domain_blank method.'
write(*,*)'======================================================================'
write(*,*)

did1 = add_domain(100_i8)
ds   = get_domain_size(did1)
call assert_equal(ds, 100, 'domain1:get_domain_size')

did2 = add_domain(10_i8)
ds   = get_domain_size(did2)
call assert_equal(ds, 10, 'domain2:get_domain_size')

! ds   = get_domain_size(-7) ! Correctly generates fatal error
! ds   = get_variable_size(did1,2) ! Correctly generates fatal error

call state_structure_info(did1)
call state_structure_info(did2)

num_domains = get_num_domains()
call assert_equal(num_domains, 2, 'get_num_domains')

var_size = get_variable_size(did1,1)
call assert_equal(var_size, 100, 'domain1:variable1:get_variable_size')

var_size = get_variable_size(did2,1)
call assert_equal(var_size, 10, 'domain2:variable1:get_variable_size')

string3 = get_variable_name(did1,1)
call assert_equal(string3, 'state', 'domain1:variable1:get_variable_name')

ds = get_num_variables(did1)
call assert_equal(ds, 1, 'domain1:get_num_variables')

ds = get_num_dims(did2,1)
call assert_equal(ds, 1, 'domain2:get_num_dims')

ds = get_dim_length(did2,1,1)
call assert_equal(ds, 10, 'domain2:variable1:dimension1:get_dim_length')

allocate(int_array(ds))
int_array = get_dim_lengths(did1,1)
call assert_equal(int_array(1), 100, 'domain1:variable1:get_dim_lengths')

int_array = get_dim_lengths(did2,1)
call assert_equal(int_array(1), 10, 'domain2:variable1:get_dim_lengths')

string3 =  get_kind_string(did1,1)
call assert_equal(string3, 'KIND_RAW_STATE_VARIABLE','domain1:variable1:get_kind_string')

ds = get_kind_index(did1,1)
call assert_equal(ds, 0, 'domain1:variable1:get_kind_index')

int_array(1) = get_varid_from_kind(did1,ds)
call assert_equal(int_array(1), 1, 'domain1:get_varid_from_kind')

! try to find a variable that does not exist
int_array(1) = get_varid_from_kind(did1,did2)
call assert_equal(int_array(1), -1, 'domain1:domain2:get_varid_from_kind')

call get_varids_from_kind(did1, ds, int_array)
call assert_equal(int_array(1), 1, 'domain1:get_varids_from_kind')

string3 = get_dim_name(did1, 1, 1)
call assert_equal(string3, 'model_size', 'domain1:variable1:dimension1:get_dim_name')

ds = get_io_num_dims(did1, 1)
call assert_equal(ds, 1, 'domain1:variable1:get_io_num_dims')

int_array = get_io_dim_ids(did1, 1)
call assert_equal(int_array, (/ 1 /), 'domain1:variable1:get_io_dim_ids')

int_array = get_io_dim_lengths(did1, 1)
call assert_equal(int_array, (/ 100 /), 'domain1:variable1:get_io_dim_lengths')

ds = get_io_num_unique_dims(did1)
call assert_equal(ds, 1, 'domain1:get_io_num_unique_dims')

string3 = get_io_unique_dim_name(did1,1)
call assert_equal(string3, 'model_size', 'domain1:variable1:get_io_unique_dim_name')

ds = get_io_unique_dim_length(did1,1)
call assert_equal(ds, 100, 'domain1:variable1:get_io_unique_dim_length')

! Routine interface only at this time ...
! ds = add_time_unlimited(did1,1)
! call assert_equal(ds, 100, 'domain1:variable1:add_time_unlimited')

ds = get_unlimited_dimid(did1)   ! no unlimited dimension in this case
call assert_equal(ds, -1, 'domain1:get_unlimited_dimid')

!>@ no way to test this
! call set_var_id(did1,1,400)

realnum = get_io_clamping_maxval(did1,1)
call assert_equal(realnum, MISSING_R8, 'domain1:variable1:get_io_clamping_maxval')

realnum = get_io_clamping_minval(did1,1)
call assert_equal(realnum, MISSING_R8, 'domain1:variable1:get_io_clamping_minval')

truefalse = do_io_clamping(did1,1)
call assert_equal(truefalse, .false., 'domain1:variable1:do_io_clamping')

truefalse = do_io_update(did1,1)
call assert_equal(truefalse, .true., 'domain1:variable1:do_io_update')

! get_index_start_from_[varname,varid]

ds = get_index_start(did1,'state')
call assert_equal(ds, 1, 'domain1:variable1:get_index_start_from_varname')

ds = get_index_start(did1,1)
call assert_equal(ds, 1, 'domain1:variable1:get_index_start_from_varid')

! get_index_end_from_[varname,varid]

ds = get_index_end(did1,'state')
call assert_equal(ds, 100, 'domain1:variable1:get_index_end_from_varname')

ds = get_index_end(did1,1)
call assert_equal(ds, 100, 'domain1:variable1:get_index_end_from_varid')

bigint = get_sum_variables(1,1,did1) ! return the number of elements in variables
call assert_equal(bigint, 100_i8, 'domain1:get_sum_variables')

bigint = get_sum_variables_below(1,did1) ! return the number of elements in variables
call assert_equal(bigint, 0_i8, 'domain1:get_sum_variables_below')

index_in = 50_i8
call get_model_variable_indices(index_in, iloc, jloc, kloc, var_id, &
                                dom_id, kind_index, kind_string)
call assert_equal(iloc,      50, 'domain1:      iloc:get_model_variable_indices')
call assert_equal(jloc,       1, 'domain1:      jloc:get_model_variable_indices')
call assert_equal(kloc,       1, 'domain1:      kloc:get_model_variable_indices')
call assert_equal(var_id,     1, 'domain1:    var_id:get_model_variable_indices')
call assert_equal(dom_id,     1, 'domain1:    dom_id:get_model_variable_indices')
call assert_equal(kind_index, 0, 'domain1:kind_index:get_model_variable_indices')
call assert_equal(kind_string, 'KIND_RAW_STATE_VARIABLE', 'domain1:kind_string:get_model_variable_indices')

index_in = get_dart_vector_index(iloc, jloc, kloc, dom_id, var_id)
call assert_equal(index_in, 50_i8, 'domain1:get_dart_vector_index')

ds = get_num_varids_from_kind(did1, kind_index)
call assert_equal(ds, 1, 'domain1:get_num_varids_from_kind')

write(*,*)
write(*,*)'======================================================================'
write(*,*)'Testing the add_domain_from_file method.'
write(*,*)'======================================================================'
write(*,*)
! netcdf simple {
! dimensions:
!        level = 3 ;
!        dump_trucks = 5 ;
!        brown_trout = 5 ;
!        lat = 4 ;
!        lon = 5 ;
!        palm_trees = 5 ;
!        time = UNLIMITED ; // (1 currently)
! variables:
!         int A(level) ;
!                 A:units = "meters" ;
!                 A:long_name = "variable A" ;
!                 A:short_name = "short A" ;
!                 A:missing_value = 22 ;
!                 A:scale_factor = 0.2 ;
!                 A:add_offset = 2 ;
!         float B(level) ;
!                 B:units = "meters/second" ;
!                 B:long_name = "variable B" ;
!                 B:short_name = "short B" ;
!                 B:missing_value = 111.11 ;
!                 B:scale_factor = 0.1 ;
!                 B:add_offset = 1 ;
!         double C(level) ;
!                 C:units = "meters/kg" ;
!                 C:long_name = "variable C" ;
!                 C:short_name = "short C" ;
!                 C:missing_value = 111.11 ;
!                 C:scale_factor = 0.1 ;
!                 C:add_offset = 1 ;
!         double temp(time, lon, lat, level) ;
!                 temp:units = "palm trees" ;
!                 temp:long_name = "ambient spectacular temperature from some really great planet and season" ;
!                 temp:short_name = "temperature" ;
!         float time(time) ;
!                 time:units = "hours" ;
! 
! // global attributes:
!                 :title = "simple_file" ;

call static_init_model()
num_domains = get_num_domains()
call assert_equal(num_domains, 3, 'get_num_domains')

did3 = 3; ! The 'blank' test added 2 domains ... this is the third.

call state_structure_info(did3)

write(*,*)'Early exit ... Friday COB.'
stop ! TJH LEFT OFF HERE

var_size = get_variable_size(did3,1)
call assert_equal(var_size, 100, 'domain1:variable1:get_variable_size')

var_size = get_variable_size(did3,1)
call assert_equal(var_size, 10, 'domain2:variable1:get_variable_size')

string3 = get_variable_name(did3,1)
call assert_equal(string3, 'state', 'domain1:variable1:get_variable_name')

ds = get_num_variables(did3)
call assert_equal(ds, 1, 'domain1:get_num_variables')

ds = get_num_dims(did3,1)
call assert_equal(ds, 1, 'domain2:get_num_dims')

ds = get_dim_length(did3,1,1)
call assert_equal(ds, 10, 'domain2:variable1:dimension1:get_dim_length')

allocate(int_array(ds))
int_array = get_dim_lengths(did3,1)
call assert_equal(int_array(1), 100, 'domain1:variable1:get_dim_lengths')

int_array = get_dim_lengths(did3,1)
call assert_equal(int_array(1), 10, 'domain2:variable1:get_dim_lengths')

string3 =  get_kind_string(did3,1)
call assert_equal(string3, 'KIND_RAW_STATE_VARIABLE','domain1:variable1:get_kind_string')

ds = get_kind_index(did3,1)
call assert_equal(ds, 0, 'domain1:variable1:get_kind_index')

int_array(1) = get_varid_from_kind(did3,ds)
call assert_equal(int_array(1), 1, 'domain1:get_varid_from_kind')

! try to find a variable that does not exist
int_array(1) = get_varid_from_kind(did3,did3)
call assert_equal(int_array(1), -1, 'domain1:domain2:get_varid_from_kind')

call get_varids_from_kind(did3, ds, int_array)
call assert_equal(int_array(1), 1, 'domain1:get_varids_from_kind')

string3 = get_dim_name(did3, 1, 1)
call assert_equal(string3, 'model_size', 'domain1:variable1:dimension1:get_dim_name')

ds = get_io_num_dims(did3, 1)
call assert_equal(ds, 1, 'domain1:variable1:get_io_num_dims')

int_array = get_io_dim_ids(did3, 1)
call assert_equal(int_array, (/ 1 /), 'domain1:variable1:get_io_dim_ids')

int_array = get_io_dim_lengths(did3, 1)
call assert_equal(int_array, (/ 100 /), 'domain1:variable1:get_io_dim_lengths')

ds = get_io_num_unique_dims(did3)
call assert_equal(ds, 1, 'domain1:get_io_num_unique_dims')

string3 = get_io_unique_dim_name(did3,1)
call assert_equal(string3, 'model_size', 'domain1:variable1:get_io_unique_dim_name')

ds = get_io_unique_dim_length(did3,1)
call assert_equal(ds, 100, 'domain1:variable1:get_io_unique_dim_length')

! Routine interface only at this time ...
! ds = add_time_unlimited(did3,1)
! call assert_equal(ds, 100, 'domain1:variable1:add_time_unlimited')

ds = get_unlimited_dimid(did3)   ! no unlimited dimension in this case
call assert_equal(ds, -1, 'domain1:get_unlimited_dimid')

!>@ no way to test this
! call set_var_id(did3,1,400)

realnum = get_io_clamping_maxval(did3,1)
call assert_equal(realnum, MISSING_R8, 'domain1:variable1:get_io_clamping_maxval')

realnum = get_io_clamping_minval(did3,1)
call assert_equal(realnum, MISSING_R8, 'domain1:variable1:get_io_clamping_minval')

truefalse = do_io_clamping(did3,1)
call assert_equal(truefalse, .false., 'domain1:variable1:do_io_clamping')

truefalse = do_io_update(did3,1)
call assert_equal(truefalse, .true., 'domain1:variable1:do_io_update')

! get_index_start_from_[varname,varid]

ds = get_index_start(did3,'state')
call assert_equal(ds, 1, 'domain1:variable1:get_index_start_from_varname')

ds = get_index_start(did3,1)
call assert_equal(ds, 1, 'domain1:variable1:get_index_start_from_varid')

! get_index_end_from_[varname,varid]

ds = get_index_end(did3,'state')
call assert_equal(ds, 100, 'domain1:variable1:get_index_end_from_varname')

ds = get_index_end(did3,1)
call assert_equal(ds, 100, 'domain1:variable1:get_index_end_from_varid')

bigint = get_sum_variables(1,1,did3) ! return the number of elements in variables
call assert_equal(bigint, 100_i8, 'domain1:get_sum_variables')

bigint = get_sum_variables_below(1,did3) ! return the number of elements in variables
call assert_equal(bigint, 0_i8, 'domain1:get_sum_variables_below')

index_in = 50_i8
call get_model_variable_indices(index_in, iloc, jloc, kloc, var_id, &
                                dom_id, kind_index, kind_string)
call assert_equal(iloc,      50, 'domain1:      iloc:get_model_variable_indices')
call assert_equal(jloc,       1, 'domain1:      jloc:get_model_variable_indices')
call assert_equal(kloc,       1, 'domain1:      kloc:get_model_variable_indices')
call assert_equal(var_id,     1, 'domain1:    var_id:get_model_variable_indices')
call assert_equal(dom_id,     1, 'domain1:    dom_id:get_model_variable_indices')
call assert_equal(kind_index, 0, 'domain1:kind_index:get_model_variable_indices')
call assert_equal(kind_string, 'KIND_RAW_STATE_VARIABLE', 'domain1:kind_string:get_model_variable_indices')

index_in = get_dart_vector_index(iloc, jloc, kloc, dom_id, var_id)
call assert_equal(index_in, 50_i8, 'domain1:get_dart_vector_index')

ds = get_num_varids_from_kind(did3, kind_index)
call assert_equal(ds, 1, 'domain1:get_num_varids_from_kind')

write(*,*)
write(*,*)'======================================================================'
write(*,*)'Testing the add_domain_from_spec method.'
write(*,*)'======================================================================'
write(*,*)

! call add_dimension_to_variable(did1, 1, 'lon', 4)
! ds = get_num_dims(did1,1)
! call assert_equal(ds, 2, 'domain1:variable1:get_num_dims')

! ! Should do nothing for this (add_domain_blank) method.
! call finished_adding_domain(did1)

call finalize_utilities('test_sizes')

end program

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
