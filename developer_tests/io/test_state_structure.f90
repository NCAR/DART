! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program test_state_structure

use           types_mod, only : r8, i8, i4, MISSING_R8

use       utilities_mod, only : register_module, initialize_utilities, &
                                finalize_utilities, error_handler, &
                                E_ERR, E_MSG

use        obs_kind_mod, only : get_name_for_quantity,  &
                                QTY_U_WIND_COMPONENT,   &
                                QTY_V_WIND_COMPONENT,   &
                                QTY_SURFACE_PRESSURE,   &
                                QTY_STATE_VARIABLE,     &
                                QTY_TEMPERATURE,        &
                                QTY_SALINITY

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

use netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

! global variables
integer, parameter :: VAR1 = 1
integer, parameter :: VAR2 = 2
integer, parameter :: VAR3 = 3
integer, parameter :: VAR4 = 4
integer, parameter :: VAR5 = 5 ! variable does not exist

integer, parameter :: DIM1 = 1
integer, parameter :: DIM2 = 2
integer, parameter :: DIM3 = 3
integer, parameter :: DIM4 = 4

integer, parameter :: b1sz = 100
integer, parameter :: b2sz =  10

! number of variables for adding domain from file
integer, parameter :: f1nv = 4
integer, parameter :: f2nv = 3

! number of variables for adding domain from spec
integer, parameter :: spnv = 4

! spec dimension sizes
integer, parameter :: nsd1 = 2
integer, parameter :: nsd2 = 3
integer, parameter :: nsd3 = 4
integer, parameter :: nsd4 = 5

! spec dimension sizes
integer, parameter :: nsdid1 = 1
integer, parameter :: nsdid2 = 2
integer, parameter :: nsdid3 = 3
integer, parameter :: nsdid4 = 4

! spec variable sizes
integer, parameter :: spv1sz = nsd1
integer, parameter :: spv2sz = nsd1*nsd2
integer, parameter :: spv3sz = nsd1*nsd2*nsd3
integer, parameter :: spv4sz = nsd1*nsd2*nsd3*nsd4

! file dimension sizes
integer, parameter :: f1d1 = 3
integer, parameter :: f1d2 = 4
integer, parameter :: f1d3 = 5
integer, parameter :: f1d4 = 1

integer, parameter :: f2d1 = 4
integer, parameter :: f2d2 = 5
integer, parameter :: f2d3 = 6
integer, parameter :: f2d4 = 1

! file dimension ids
integer, parameter :: f1did1 = 1
integer, parameter :: f1did2 = 4
integer, parameter :: f1did3 = 5
integer, parameter :: f1did4 = 7

integer, parameter :: f2did1 = 1
integer, parameter :: f2did2 = 2
integer, parameter :: f2did3 = 3
integer, parameter :: f2did4 = 4

integer, parameter :: f1v1sz = f1d1
integer, parameter :: f1v2sz = f1d1
integer, parameter :: f1v3sz = f1d1
integer, parameter :: f1v4sz = f1d1*f1d2*f1d3

integer, parameter :: f2v1sz = f2d1
integer, parameter :: f2v2sz = f2d1
integer, parameter :: f2v3sz = f2d1*f2d2*f2d3

! domain sizes
integer(i8), parameter :: dsize1 = b1sz 
integer(i8), parameter :: dsize2 = b2sz 
integer(i8), parameter :: dsize3 = f1v1sz+f1v2sz+f1v3sz+f1v4sz
integer(i8), parameter :: dsize4 = f2v1sz+f2v2sz+f2v3sz
integer(i8), parameter :: dsize5 = spv1sz+spv2sz+spv3sz+spv4sz

! model sizes after adding a domain
integer(i8), parameter :: msize1 = dsize1
integer(i8), parameter :: msize2 = dsize1+dsize2
integer(i8), parameter :: msize3 = dsize1+dsize2+dsize3
integer(i8), parameter :: msize4 = dsize1+dsize2+dsize3+dsize4
integer(i8), parameter :: msize5 = dsize1+dsize2+dsize3+dsize4+dsize5

integer :: did1, did2, did3, did4, did5

integer, parameter :: ndoms = 5

integer(i8), parameter, dimension(ndoms+1) :: msize = (/ 0_i8, msize1, msize2, msize3, msize4, msize5 /)

! simple1.nc metadata
character(len=NF90_MAX_NAME), dimension(f1nv)   ::   var_names1
integer,                      dimension(f1nv)   ::   kind_list1
real(r8),                     dimension(f1nv,2) ::  clamp_vals1
logical,                      dimension(f1nv)   :: update_list1

! simple2.nc metadata
character(len=NF90_MAX_NAME), dimension(f2nv)   ::   var_names2
integer,                      dimension(f2nv)   ::   kind_list2
real(r8),                     dimension(f2nv,2) ::  clamp_vals2
logical,                      dimension(f2nv)   :: update_list2

! domain information for add_domain_from_spec
character(len=NF90_MAX_NAME), dimension(spnv)   ::   var_names3
integer,                      dimension(spnv)   ::   kind_list3
real(r8),                     dimension(spnv,2) ::  clamp_vals3
logical,                      dimension(spnv)   :: update_list3

character(len=32)            :: kindString
character(len=NF90_MAX_NAME) :: varName
integer(i8)                  :: int8Val
integer                      :: intVal
real(r8)                     :: realVal
logical                      :: trueFalse

integer(i8) :: state_indx
integer     :: dom

character(len=512) :: string1

! namelist variables
logical :: debug = .false.

! namelist items we are going to create/overwrite
namelist /test_state_structure_nml/  debug

! main code here
 
! initialize the dart libs
call initialize_module()

call error_handler(E_ERR,'test_state_structure ',&
                   'Has not been tested yet with new naming conventions.',source,revision,revdate)

call initialize_domains()

did1 = add_domain(int(b1sz,i8))
did2 = add_domain(int(b2sz,i8))
did3 = add_domain('simple1.nc', f1nv,   var_names1, &
                  kind_list1, clamp_vals1, update_list1)
did4 = add_domain('simple2.nc', f2nv,   var_names2, &
                  kind_list2, clamp_vals2, update_list2)
did5 = add_domain(spnv, var_names3, kind_list3, clamp_vals3, update_list3)

! add dimensions for domains from spec
call fill_domain_structure_for_spec(did5)

if ( debug ) then
   write(*,*)
   write(*,*)'=============================================================='
   write(*,*)'State structure information.'
   write(*,*)'=============================================================='
   write(*,*)
   
   call state_structure_info(did1)
   call state_structure_info(did2)
   call state_structure_info(did3)
   call state_structure_info(did4)
   call state_structure_info(did5)
endif 

write(*,*)
write(*,*)'=============================================================='
write(*,*)'Testing domain sizes.'
write(*,*)'=============================================================='
write(*,*)

call assert_equal(get_num_domains(), ndoms, 'get_num_domains')

write(*,*)' ... Testing Domain 1'
call test_sizes_domain(did1, exp_dom_size     = b1sz, &
                             exp_num_dom_vars = 1, & ! only one var for blank
                             exp_unlim_dimid  = -1) ! no unlim dim

write(*,*)' ... Testing Domain 2'
call test_sizes_domain(did2, exp_dom_size     = b2sz, &
                             exp_num_dom_vars = 1, &
                             exp_unlim_dimid  = -1)

write(*,*)' ... Testing Domain 3'
call test_sizes_domain(did3, exp_dom_size     = f1v1sz+f1v2sz+f1v3sz+f1v4sz, &
                             exp_num_dom_vars = 4, &
                             exp_unlim_dimid  = 7)

write(*,*)' ... Testing Domain 4'
call test_sizes_domain(did4, exp_dom_size     = f2v1sz+f2v2sz+f2v3sz, &
                             exp_num_dom_vars = 3, &
                             exp_unlim_dimid  = 4)

write(*,*)' ... Testing Domain 5'
call test_sizes_domain(did5, exp_dom_size     = spv1sz+spv2sz+spv3sz+spv4sz, &
                             exp_num_dom_vars = 4, &
                             exp_unlim_dimid  = -1)

write(*,*)
write(*,*)'=============================================================='
write(*,*)'Testing variable information.'
write(*,*)'=============================================================='
write(*,*)

write(*,*)' ... Testing Domain 1'
call test_variable_info(did1, VAR1, &
                        exp_var_size    = b1sz, &
                        exp_num_dims    = 1, &
                        exp_kind_indx   =  QTY_STATE_VARIABLE, &
                        exp_kind_string = 'QTY_STATE_VARIABLE',&
                        exp_var_name    = 'state')

write(*,*)' ... Testing Domain 2'
call test_variable_info(did2, VAR1, &
                        exp_var_size    = b2sz, &
                        exp_num_dims    = 1, &
                        exp_kind_indx   =  QTY_STATE_VARIABLE, &
                        exp_kind_string = 'QTY_STATE_VARIABLE',&
                        exp_var_name    = 'state')

write(*,*)' ... Testing Domain 3 Var 1'
call test_variable_info(did3, VAR1, &
                        exp_var_size    = f1d1, &
                        exp_num_dims    = 1, &
                        exp_kind_indx   =  QTY_STATE_VARIABLE, &
                        exp_kind_string = 'QTY_STATE_VARIABLE',&
                        exp_var_name    = 'A')

write(*,*)' ... Testing Domain 3 Var 4'
call test_variable_info(did3, VAR4, &
                        exp_var_size    = f1d1*f1d2*f1d3, &
                        exp_num_dims    = 3, &
                        exp_kind_indx   =  QTY_TEMPERATURE, &
                        exp_kind_string = 'QTY_TEMPERATURE',&
                        exp_var_name    = 'temp')

write(*,*)' ... Testing Domain 4 Var 1'
call test_variable_info(did4, VAR1, &
                        exp_var_size    = f2d1, &
                        exp_num_dims    = 1, &
                        exp_kind_indx   =  QTY_STATE_VARIABLE, &
                        exp_kind_string = 'QTY_STATE_VARIABLE',&
                        exp_var_name    = 'B')

call test_variable_info(did4, VAR3, &
                        exp_var_size    = f2d1*f2d2*f2d3, &
                        exp_num_dims    = 3, &
                        exp_kind_indx   =  QTY_TEMPERATURE, &
                        exp_kind_string = 'QTY_TEMPERATURE',&
                        exp_var_name    = 'temp')

write(*,*)' ... Testing Domain 5 Var 1'
call test_variable_info(did5, VAR1, &
                        exp_var_size    = spv1sz, &
                        exp_num_dims    = 1, &
                        exp_kind_indx   = QTY_SURFACE_PRESSURE, &
                        exp_kind_string = 'QTY_SURFACE_PRESSURE',&
                        exp_var_name    = 'PS')

write(*,*)' ... Testing Domain 5 Var 4'
call test_variable_info(did5, VAR4, &
                        exp_var_size    = spv4sz, &
                        exp_num_dims    = 4, &
                        exp_kind_indx   =  QTY_V_WIND_COMPONENT, &
                        exp_kind_string = 'QTY_V_WIND_COMPONENT',&
                        exp_var_name    = 'V ')

write(*,*)
write(*,*)'=============================================================='
write(*,*)'Testing dimension information.'
write(*,*)'=============================================================='
write(*,*)

write(*,*)' ... Testing Domain 1 Var 1'
call test_dimension_info(did1, VAR1, & 
                         dim_ids            = (/DIM1/), &
                         exp_dim_names      = (/'domain_size'/), &
                         exp_dim_lengths    = (/b1sz/), &
                         io_dim_ids         = (/DIM1/), &
                         exp_io_numdims     = 1, &
                         exp_io_dim_names   = (/'domain_size'/), &
                         exp_io_dim_lengths = (/ b1sz/), &
                         exp_io_dim_ids     = (/ 1/), &
                         exp_io_unique_dim_ids = (/1/), &
                         exp_io_unique_numdims = 1, &
                         exp_io_unique_dim_length = (/ b1sz/), &  
                         exp_io_unique_dim_names = (/'domain_size'/) )

write(*,*)' ... Testing Domain 2 Var 1'
call test_dimension_info(did2, VAR1, & 
                         dim_ids            = (/DIM1/), &
                         exp_dim_names      = (/'domain_size'/), &
                         exp_dim_lengths    = (/b2sz/), &
                         io_dim_ids         = (/DIM1/), &
                         exp_io_numdims     = 1, &
                         exp_io_dim_names   = (/'domain_size'/), &
                         exp_io_dim_lengths = (/ b2sz/), &
                         exp_io_dim_ids     = (/ 1/), &
                         exp_io_unique_dim_ids = (/1/), &
                         exp_io_unique_numdims = 1, &
                         exp_io_unique_dim_length = (/ b2sz/), &  
                         exp_io_unique_dim_names = (/'domain_size'/) )

write(*,*)' ... Testing Domain 3 Var 1'
call test_dimension_info(did3, VAR1, & 
                         dim_ids            = (/DIM1/), &
                         exp_dim_names      = (/'level'/), &
                         exp_dim_lengths    = (/ f1d1/), &
                         io_dim_ids         = (/DIM1/), &
                         exp_io_numdims     = 1, &
                         exp_io_dim_names   = (/'level', 'lat  ', 'lon  ', 'time '/), &
                         exp_io_dim_lengths = (/ f1d1/), &
                         exp_io_dim_ids     = (/DIM1/), &
                         exp_io_unique_dim_ids = (/DIM1, DIM2, DIM3, DIM4/), &
                         exp_io_unique_numdims = 4, &
                         exp_io_unique_dim_length = (/ f1d1, f1d2, f1d3, f1d4/), & 
                         exp_io_unique_dim_names = (/'level', 'lat  ', 'lon  ','time '/) )

write(*,*)' ... Testing Domain 3 Var 4'
call test_dimension_info(did3, VAR4, & 
                         dim_ids            = (/DIM1, DIM2, DIM3/), &
                         exp_dim_names      = (/'level', 'lat  ', 'lon  '/), &
                         exp_dim_lengths    = (/ f1d1, f1d2, f1d3/), &
                         io_dim_ids         = (/DIM1, DIM2, DIM3, DIM4/), &
                         exp_io_numdims     = 4, &
                         exp_io_dim_names   = (/'level', 'lat  ', 'lon  ','time '/), &
                         exp_io_dim_lengths = (/ f1d1, f1d2, f1d3, f1d4/), &
                         exp_io_dim_ids     = (/f1did1, f1did2, f1did3, f1did4/), &
                         exp_io_unique_dim_ids = (/DIM1, DIM2, DIM3, DIM4/), &
                         exp_io_unique_numdims = 4, &
                         exp_io_unique_dim_length = (/ f1d1, f1d2, f1d3, f1d4/), &  
                         exp_io_unique_dim_names = (/'level', 'lat  ', 'lon  ','time '/) )

write(*,*)' ... Testing Domain 4 Var 2'
call test_dimension_info(did4, VAR2, & 
                         dim_ids            = (/DIM1/), &
                         exp_dim_names      = (/'level'/), &
                         exp_dim_lengths    = (/ f2d1/), &
                         io_dim_ids         = (/DIM1/), &
                         exp_io_numdims     = 1, &
                         exp_io_dim_names   = (/'level'/), &
                         exp_io_dim_lengths = (/ f2d1/), &
                         exp_io_dim_ids     = (/ f2did1/), &
                         exp_io_unique_dim_ids = (/DIM1, DIM2, DIM3, DIM4/), &
                         exp_io_unique_numdims = 4, &
                         exp_io_unique_dim_length = (/ f2d1, f2d2, f2d3, f2d4/), &  
                         exp_io_unique_dim_names = (/'level', 'lat  ', 'lon  ','time '/) )

write(*,*)' ... Testing Domain 4 Var 3'
call test_dimension_info(did4, VAR3, & 
                         dim_ids            = (/DIM1, DIM2, DIM3/), &
                         exp_dim_names      = (/'level', 'lat  ', 'lon  '/), &
                         exp_dim_lengths    = (/ f2d1, f2d2, f2d3/), &
                         io_dim_ids         = (/DIM1, DIM2, DIM3, DIM4/), &
                         exp_io_numdims     = 4, &
                         exp_io_dim_names   = (/'level', 'lat  ', 'lon  ','time '/), &
                         exp_io_dim_lengths = (/ f2d1, f2d2, f2d3, f2d4/), &
                         exp_io_dim_ids     = (/ f2did1, f2did2, f2did3, f2did4/), &
                         exp_io_unique_dim_ids = (/DIM1, DIM2, DIM3, DIM4/), &
                         exp_io_unique_numdims = 4, &
                         exp_io_unique_dim_length = (/ f2d1, f2d2, f2d3, f2d4/), &  
                         exp_io_unique_dim_names = (/'level', 'lat  ', 'lon  ','time '/) )

write(*,*)' ... Testing Domain 5 Var 2'
call test_dimension_info(did5, VAR2, & 
                         dim_ids            = (/DIM1, DIM2/), &
                         exp_dim_names      = (/'dim1', 'dim2'/), &
                         exp_dim_lengths    = (/ nsd1, nsd2/), &
                         io_dim_ids         = (/DIM1, DIM2/), &
                         exp_io_numdims     = 2, &
                         exp_io_dim_names   = (/'dim1', 'dim2'/), &
                         exp_io_dim_lengths = (/ nsd1, nsd2/), &
                         exp_io_dim_ids     = (/ nsdid1, nsdid2/), &
                         exp_io_unique_dim_ids = (/1,2,3,4,5,6,7,8,9,10/), &
                         exp_io_unique_numdims = 10, &
                         exp_io_unique_dim_length = (/ nsd1, &
                                                       nsd1, nsd2, &
                                                       nsd1, nsd2, nsd3, &
                                                       nsd1, nsd2, nsd3,nsd4/), &  
                         exp_io_unique_dim_names = (/'dim1', &
                                                     'dim1', 'dim2', &
                                                     'dim1', 'dim2', 'dim3', &
                                                     'dim1', 'dim2', 'dim3', 'dim4'/) )
write(*,*)' ... Testing Domain 5 Var 4'
call test_dimension_info(did5, VAR4, & 
                         dim_ids            = (/DIM1, DIM2, DIM3, DIM4/), &
                         exp_dim_names      = (/'dim1', 'dim2', 'dim3', 'dim4'/), &
                         exp_dim_lengths    = (/ nsd1, nsd2, nsd3, nsd4/), &
                         io_dim_ids         = (/DIM1, DIM2, DIM3, DIM4/), &
                         exp_io_numdims     = 4, &
                         exp_io_dim_names   = (/'dim1', 'dim2', 'dim3', 'dim4'/), &
                         exp_io_dim_lengths = (/ nsd1, nsd2, nsd3, nsd4/), &
                         exp_io_dim_ids     = (/ nsdid1, nsdid2, nsdid3, nsdid4/), &
                         exp_io_unique_dim_ids = (/1,2,3,4,5,6,7,8,9,10/), &
                         exp_io_unique_numdims = 10, &
                         exp_io_unique_dim_length = (/ nsd1, &
                                                       nsd1, nsd2, &
                                                       nsd1, nsd2, nsd3, &
                                                       nsd1, nsd2, nsd3,nsd4/), &  
                         exp_io_unique_dim_names = (/'dim1', &
                                                     'dim1', 'dim2', &
                                                     'dim1', 'dim2', 'dim3', &
                                                     'dim1', 'dim2', 'dim3', 'dim4'/) )

write(*,*)
write(*,*)'=============================================================='
write(*,*)'Testing update and clamping.'
write(*,*)'=============================================================='
write(*,*)

write(*,*)' ... Testing Domain 1 Var 1'
call test_update_and_clamping(did1, VAR1, &
                              exp_clamp_min = MISSING_R8, &
                              exp_clamp_max = MISSING_R8, &
                              exp_clamping  = .false., &
                              exp_update    = .true.)

write(*,*)' ... Testing Domain 2 Var 1'
call test_update_and_clamping(did2, VAR1, &
                              exp_clamp_min = MISSING_R8, &
                              exp_clamp_max = MISSING_R8, &
                              exp_clamping  = .false., &
                              exp_update    = .true.)

write(*,*)' ... Testing Domain 3 Var 1'
call test_update_and_clamping(did3, VAR1, &
                              exp_clamp_min = clamp_vals1(VAR1,1), &
                              exp_clamp_max = clamp_vals1(VAR1,2), &
                              exp_clamping  = .true., &
                              exp_update    = update_list1(VAR1))

write(*,*)' ... Testing Domain 3 Var 4'
call test_update_and_clamping(did3, VAR4, &
                              exp_clamp_min = clamp_vals1(VAR4,1), &
                              exp_clamp_max = clamp_vals1(VAR4,2), &
                              exp_clamping  = .true., &
                              exp_update    = update_list1(VAR4))

write(*,*)' ... Testing Domain 4 Var 1'
call test_update_and_clamping(did4, VAR1, &
                              exp_clamp_min = clamp_vals2(VAR1,1), &
                              exp_clamp_max = clamp_vals2(VAR1,2), &
                              exp_clamping  = .false., &
                              exp_update    = update_list1(VAR1))

write(*,*)' ... Testing Domain 4 Var 3'
call test_update_and_clamping(did4, VAR3, &
                              exp_clamp_min = clamp_vals2(VAR3,1), &
                              exp_clamp_max = clamp_vals2(VAR3,2), &
                              exp_clamping  = .true., &
                              exp_update    = update_list1(VAR3))

write(*,*)' ... Testing Domain 5 Var 1'
call test_update_and_clamping(did5, VAR1, &
                              exp_clamp_min = clamp_vals3(VAR1,1), &
                              exp_clamp_max = clamp_vals3(VAR1,2), &
                              exp_clamping  = .false., &
                              exp_update    = update_list3(VAR1))

write(*,*)' ... Testing Domain 5 Var 3'
call test_update_and_clamping(did5, VAR3, &
                              exp_clamp_min = clamp_vals3(VAR3,1), &
                              exp_clamp_max = clamp_vals3(VAR3,2), &
                              exp_clamping  = .true., &
                              exp_update    = update_list3(VAR3))

write(*,*)
write(*,*)'=============================================================='
write(*,*)'Testing start and end indices'
write(*,*)'=============================================================='
write(*,*)

write(*,*)' ... Testing Domain 1 Var 1'
call test_indices(did1, VAR1, &
                  var_name  = 'state', &
                  exp_start = int(1,i8) , &
                  exp_end   = int(b1sz,i8))

write(*,*)' ... Testing Domain 2 Var 1'
call test_indices(did2, VAR1, &
                  var_name  = 'state', &
                  exp_start = msize1+1 , &
                  exp_end   = msize1+b2sz)

write(*,*)' ... Testing Domain 3 Var 1'
call test_indices(did3, VAR1, &
                  var_name  = 'A', &
                  exp_start = msize2+1 , &
                  exp_end   = msize2+f1v1sz )

write(*,*)' ... Testing Domain 3 Var 4'
call test_indices(did3, VAR4, &
                  var_name  = 'temp', &
                  exp_start = msize2+1+f1v1sz+f1v2sz+f1v3sz, &
                  exp_end   = msize2+  f1v1sz+f1v2sz+f1v3sz+f1v4sz )

write(*,*)' ... Testing Domain 4 Var 1'
call test_indices(did4, VAR1, &
                  var_name  = 'B', &
                  exp_start = msize3+1 , &
                  exp_end   = msize3+f2v1sz )

write(*,*)' ... Testing Domain 4 Var 3'
call test_indices(did4, VAR3, &
                  var_name  = 'temp', &
                  exp_start = msize3+1+f2v1sz+f2v2sz, &
                  exp_end   = msize3 + f2v1sz+f2v2sz+f2v3sz )

write(*,*)' ... Testing Domain 5 Var 2'
call test_indices(did5, VAR2, &
                  var_name  = 'T', &
                  exp_start = msize4+1+spv1sz, &
                  exp_end   = msize4 + spv1sz +spv2sz )

write(*,*)' ... Testing Domain 5 Var 4'
call test_indices(did5, VAR4, &
                  var_name  = 'V', &
                  exp_start = msize4+1+spv1sz+spv2sz+spv3sz, &
                  exp_end   = msize4 + spv1sz+spv2sz+spv3sz+spv4sz )

write(*,*)
write(*,*)'=============================================================='
write(*,*)'Testing sum variables below'
write(*,*)'=============================================================='
write(*,*)

write(*,*)' ... Testing Domain 1 Var 1'
call test_sum_variables_below(did1, VAR1, VAR1, &
                              exp_sum_range  = int(b1sz,i8), & 
                              exp_sum_below1 = int(0,i8), &
                              exp_sum_below2 = int(0,i8))

write(*,*)' ... Testing Domain 2 Var 1'
call test_sum_variables_below(did2, VAR1, VAR1, &
                              exp_sum_range  = int(b2sz,i8), & 
                              exp_sum_below1 = int(b1sz,i8), &
                              exp_sum_below2 = int(b1sz,i8))

write(*,*)' ... Testing Domain 3 Var 2 to Var 3'
call test_sum_variables_below(did3, VAR2, VAR3, &
                              exp_sum_range  = int(f1v2sz+f1v3sz,i8), & 
                              exp_sum_below1 = int(msize2+f1v1sz,i8), &
                              exp_sum_below2 = int(msize2+f1v1sz+f1v2sz,i8))

write(*,*)' ... Testing Domain 4 Var 1 to Var 3'
call test_sum_variables_below(did4, VAR1, VAR3, &
                              exp_sum_range  = int(f2v1sz+f2v2sz+f2v3sz,i8), & 
                              exp_sum_below1 = int(msize3,i8), &
                              exp_sum_below2 = int(msize3+f2v1sz+f2v2sz,i8))

write(*,*)' ... Testing Domain 5 Var 1 to Var 3'
call test_sum_variables_below(did5, VAR1, VAR3, &
                              exp_sum_range  = int(spv1sz+spv2sz+spv3sz,i8), & 
                              exp_sum_below1 = int(msize4,i8), &
                              exp_sum_below2 = int(msize4+spv1sz+spv2sz,i8))

write(*,*)
write(*,*)'=============================================================='
write(*,*)'Testing state indicies'
write(*,*)'=============================================================='
write(*,*)

! DOMAIN1
state_indx = 10_i8
dom = find_domain(state_indx)
write(*,'(A,I5,A,I1)') '  ... Testing state_index ', state_indx, ' Domain ', dom
call test_state_indicies(state_indx, &
                         exp_iloc = 10, &
                         exp_jloc =  1, &
                         exp_kloc =  1, &
                         exp_varid = 1, &
                         exp_domid = dom, &
                         exp_kind_index  =  QTY_STATE_VARIABLE, &
                         exp_kind_string = 'QTY_STATE_VARIABLE')

! DOMAIN2
state_indx = 109_i8
dom = find_domain(state_indx)
write(*,'(A,I5,A,I1)') '  ... Testing state_index ', state_indx, ' Domain ', dom
call test_state_indicies(state_indx, &
                         exp_iloc =  9, &
                         exp_jloc =  1, &
                         exp_kloc =  1, &
                         exp_varid = 1, &
                         exp_domid = dom, &
                         exp_kind_index  =  QTY_STATE_VARIABLE, &
                         exp_kind_string = 'QTY_STATE_VARIABLE')

! DOMAIN3
state_indx = 112_i8
dom = find_domain(state_indx)
write(*,'(A,I5,A,I1)') '  ... Testing state_index ', state_indx, ' Domain ', dom
call test_state_indicies(state_indx, &
                         exp_iloc =  2, &
                         exp_jloc =  1, &
                         exp_kloc =  1, &
                         exp_varid = 1, &
                         exp_domid = dom, &
                         exp_kind_index  =  QTY_STATE_VARIABLE, &
                         exp_kind_string = 'QTY_STATE_VARIABLE')

! DOMAIN3
state_indx = 115_i8
dom = find_domain(state_indx)
write(*,'(A,I5,A,I1)') '  ... Testing state_index ', state_indx, ' Domain ', dom
call test_state_indicies(state_indx, &
                         exp_iloc =  2, &
                         exp_jloc =  1, &
                         exp_kloc =  1, &
                         exp_varid = 2, &
                         exp_domid = dom, &
                         exp_kind_index  =  QTY_STATE_VARIABLE, &
                         exp_kind_string = 'QTY_STATE_VARIABLE')

! DOMAIN3
state_indx = 179_i8
dom = find_domain(state_indx)
write(*,'(A,I5,A,I1)') '  ... Testing state_index ', state_indx, ' Domain ', dom
call test_state_indicies(state_indx, &
                         exp_iloc =  3, &
                         exp_jloc =  4, &
                         exp_kloc =  5, &
                         exp_varid = 4, &
                         exp_domid = dom, &
                         exp_kind_index  =  QTY_TEMPERATURE, &
                         exp_kind_string = 'QTY_TEMPERATURE')

! DOMAIN3
state_indx = 178_i8
dom = find_domain(state_indx)
write(*,'(A,I5,A,I1)') '  ... Testing state_index ', state_indx, ' Domain ', dom
call test_state_indicies(state_indx, &
                         exp_iloc =  2, &
                         exp_jloc =  4, &
                         exp_kloc =  5, &
                         exp_varid = 4, &
                         exp_domid = dom, &
                         exp_kind_index  =  QTY_TEMPERATURE, &
                         exp_kind_string = 'QTY_TEMPERATURE')

! DOMAIN4
state_indx = 179_i8
dom = find_domain(state_indx)
write(*,'(A,I5,A,I1)') '  ... Testing state_index ', state_indx, ' Domain ', dom
call test_state_indicies(state_indx, &
                         exp_iloc =  3, &
                         exp_jloc =  4, &
                         exp_kloc =  5, &
                         exp_varid = 4, &
                         exp_domid = dom, &
                         exp_kind_index  =  QTY_TEMPERATURE, &
                         exp_kind_string = 'QTY_TEMPERATURE')

! DOMAIN4
state_indx = 187_i8
dom = find_domain(state_indx)
write(*,'(A,I5,A,I1)') '  ... Testing state_index ', state_indx, ' Domain ', dom
call test_state_indicies(state_indx, &
                         exp_iloc =  4, &
                         exp_jloc =  1, &
                         exp_kloc =  1, &
                         exp_varid = 2, &
                         exp_domid = dom, &
                         exp_kind_index  =  QTY_SALINITY, &
                         exp_kind_string = 'QTY_SALINITY')

! DOMAIN4
state_indx = 194_i8
dom = find_domain(state_indx)
write(*,'(A,I5,A,I1)') '  ... Testing state_index ', state_indx, ' Domain ', dom
call test_state_indicies(state_indx, &
                         exp_iloc =  3, &
                         exp_jloc =  2, &
                         exp_kloc =  1, &
                         exp_varid = 3, &
                         exp_domid = dom, &
                         exp_kind_index  =  QTY_TEMPERATURE, &
                         exp_kind_string = 'QTY_TEMPERATURE')

! DOMAIN4
state_indx = 307_i8
dom = find_domain(state_indx)
write(*,'(A,I5,A,I1)') '  ... Testing state_index ', state_indx, ' Domain ', dom
call test_state_indicies(state_indx, &
                         exp_iloc =  4, &
                         exp_jloc =  5, &
                         exp_kloc =  6, &
                         exp_varid = 3, &
                         exp_domid = dom, &
                         exp_kind_index  =  QTY_TEMPERATURE, &
                         exp_kind_string = 'QTY_TEMPERATURE')

! DOMAIN5
state_indx = 308_i8
dom = find_domain(state_indx)
write(*,'(A,I5,A,I1)') '  ... Testing state_index ', state_indx, ' Domain ', dom
call test_state_indicies(state_indx, &
                         exp_iloc =  1, &
                         exp_jloc =  1, &
                         exp_kloc =  1, &
                         exp_varid = 1, &
                         exp_domid = dom, &
                         exp_kind_index  =  QTY_SURFACE_PRESSURE, &
                         exp_kind_string = 'QTY_SURFACE_PRESSURE')

! DOMAIN5
state_indx = 318_i8
dom = find_domain(state_indx)
write(*,'(A,I5,A,I1)') '  ... Testing state_index ', state_indx, ' Domain ', dom
call test_state_indicies(state_indx, &
                         exp_iloc =  1, &
                         exp_jloc =  2, &
                         exp_kloc =  1, &
                         exp_varid = 3, &
                         exp_domid = dom, &
                         exp_kind_index  =  QTY_U_WIND_COMPONENT, &
                         exp_kind_string = 'QTY_U_WIND_COMPONENT')

write(*,*)
write(*,*)'=============================================================='
write(*,*)'Testing varids from kind'
write(*,*)'=============================================================='
write(*,*)

write(*,'(2A)')'  ... Testing Domain 1 Kind ', trim(get_name_for_quantity(QTY_STATE_VARIABLE))
call test_varids_from_kind(did1, QTY_STATE_VARIABLE, &
                           exp_num_varids_from_kind = 1, &
                           exp_varids_from_kind = (/1/) )

write(*,'(2A)')'  ... Testing Domain 1 Kind ', trim(get_name_for_quantity(QTY_TEMPERATURE))
call test_varids_from_kind(did1, QTY_TEMPERATURE, &
                           exp_num_varids_from_kind = 0, &
                           exp_varids_from_kind = (/-1/) )

write(*,'(2A)')'  ... Testing Domain 2 Kind ', trim(get_name_for_quantity(QTY_STATE_VARIABLE))
call test_varids_from_kind(did2, QTY_STATE_VARIABLE, &
                           exp_num_varids_from_kind = 1, &
                           exp_varids_from_kind = (/1/) )

write(*,'(2A)')'  ... Testing Domain 2 Kind ', trim(get_name_for_quantity(QTY_TEMPERATURE))
call test_varids_from_kind(did2, QTY_TEMPERATURE, &
                           exp_num_varids_from_kind = 0, &
                           exp_varids_from_kind = (/-1/) )

write(*,'(2A)')'  ... Testing Domain 3 Kind ', trim(get_name_for_quantity(QTY_STATE_VARIABLE))
call test_varids_from_kind(did3, QTY_STATE_VARIABLE, &
                           exp_num_varids_from_kind = 3, &
                           exp_varids_from_kind = (/1, 2, 3/) )

write(*,'(2A)')'  ... Testing Domain 3 Kind ', trim(get_name_for_quantity(QTY_TEMPERATURE))
call test_varids_from_kind(did3, QTY_TEMPERATURE, &
                           exp_num_varids_from_kind = 1, &
                           exp_varids_from_kind = (/4/) )

write(*,'(2A)')'  ... Testing Domain 4 Kind ', trim(get_name_for_quantity(QTY_STATE_VARIABLE))
call test_varids_from_kind(did4, QTY_STATE_VARIABLE, &
                           exp_num_varids_from_kind = 1, &
                           exp_varids_from_kind = (/1/) )

write(*,'(2A)')'  ... Testing Domain 4 Kind ', trim(get_name_for_quantity(QTY_SALINITY))
call test_varids_from_kind(did4, QTY_SALINITY, &
                           exp_num_varids_from_kind = 1, &
                           exp_varids_from_kind = (/2/) )

write(*,'(2A)')'  ... Testing Domain 4 Kind ', trim(get_name_for_quantity(QTY_TEMPERATURE))
call test_varids_from_kind(did4, QTY_TEMPERATURE, &
                           exp_num_varids_from_kind = 1, &
                           exp_varids_from_kind = (/3/) )

write(*,'(2A)')'  ... Testing Domain 5 Kind ', trim(get_name_for_quantity(QTY_SURFACE_PRESSURE))
call test_varids_from_kind(did5, QTY_SURFACE_PRESSURE, &
                           exp_num_varids_from_kind = 1, &
                           exp_varids_from_kind = (/1/) )

write(*,'(2A)')'  ... Testing Domain 5 Kind ', trim(get_name_for_quantity(QTY_TEMPERATURE))
call test_varids_from_kind(did5, QTY_TEMPERATURE, &
                           exp_num_varids_from_kind = 1, &
                           exp_varids_from_kind = (/2/) )

write(*,'(2A)')'  ... Testing Domain 5 Kind ', trim(get_name_for_quantity(QTY_U_WIND_COMPONENT))
call test_varids_from_kind(did5, QTY_U_WIND_COMPONENT, &
                           exp_num_varids_from_kind = 1, &
                           exp_varids_from_kind = (/3/) )

write(*,'(2A)')'  ... Testing Domain 5 Kind ', trim(get_name_for_quantity(QTY_V_WIND_COMPONENT))
call test_varids_from_kind(did5, QTY_V_WIND_COMPONENT, &
                           exp_num_varids_from_kind = 1, &
                           exp_varids_from_kind = (/4/) )

print*, ' '
! finalize test_state_structure
call error_handler(E_MSG,'test_state_structure',&
                   'Finished successfully.',source,revision,revdate)


call finalize_utilities()

! end of main code

contains

!----------------------------------------------------------------------

subroutine initialize_domains()

   !!!! DOMAIN 1 and DOMAIN 2 !!!!

   ! ! info for blank domain
   ! b1sz = 100_i8
   ! b2sz = 10_i8
   
   !!!! DOMAIN 3 !!!!

   ! info for domain from file1
   var_names1 = (/'A   ', 'B   ', 'C   ', 'temp'/)
   kind_list1 = (/QTY_STATE_VARIABLE, QTY_STATE_VARIABLE,  &
                  QTY_STATE_VARIABLE, QTY_TEMPERATURE/)
   ! variable A
   clamp_vals1(1,1) = 0.0_r8
   clamp_vals1(1,2) = MISSING_R8
   ! variable B
   clamp_vals1(2,1) = MISSING_R8
   clamp_vals1(2,2) = MISSING_R8
   ! variable C
   clamp_vals1(3,1) = MISSING_R8
   clamp_vals1(3,2) = 10.0_r8
   ! variable temp
   clamp_vals1(4,1) =  1.0_r8
   clamp_vals1(4,2) = 20.0_r8
   
   update_list1 = (/.true., .true., .false., .false./)
   
   !!!! DOMAIN 4 !!!!

   ! info for domain from file2
   var_names2 = (/'B   ', 'C   ', 'temp'/)
   kind_list2 = (/QTY_STATE_VARIABLE, QTY_SALINITY,  &
                  QTY_TEMPERATURE/)
   ! variable B
   clamp_vals2(1,1) = MISSING_R8
   clamp_vals2(1,2) = MISSING_R8
   ! variable C
   clamp_vals2(2,1) = MISSING_R8
   clamp_vals2(2,2) = 10.0_r8
   ! variable temp
   clamp_vals2(3,1) =  1.0_r8
   clamp_vals2(3,2) = 20.0_r8
   
   update_list2 = (/.true., .true., .false./)

   !!!! DOMAIN 5 !!!!

   ! infor for domain from spec
   var_names3 = (/'PS', 'T ', 'U ', 'V '/)
   kind_list3 = (/QTY_SURFACE_PRESSURE, QTY_TEMPERATURE, &
                  QTY_U_WIND_COMPONENT, QTY_V_WIND_COMPONENT/)
   ! variable PS
   clamp_vals3(1,1) = MISSING_R8
   clamp_vals3(1,2) = MISSING_R8
   ! variable T
   clamp_vals3(2,1) = MISSING_R8
   clamp_vals3(2,2) = 10.0_r8
   ! variable U
   clamp_vals3(3,1) = MISSING_R8
   clamp_vals3(3,2) = 10.0_r8
   ! variable V
   clamp_vals3(4,1) = MISSING_R8
   clamp_vals3(4,2) = 10.0_r8

   update_list3 = (/.false., .false., .false., .false./)

end subroutine initialize_domains

!----------------------------------------------------------------------

subroutine fill_domain_structure_for_spec(dom_id)
   integer,          intent(in) :: dom_id
   
   integer :: i
   do i = 1, get_num_variables(dom_id)
      if ( i >= 1) call add_dimension_to_variable(dom_id, i, 'dim1', nsd1)
      if ( i >= 2) call add_dimension_to_variable(dom_id, i, 'dim2', nsd2)
      if ( i >= 3) call add_dimension_to_variable(dom_id, i, 'dim3', nsd3)
      if ( i >= 4) call add_dimension_to_variable(dom_id, i, 'dim4', nsd4)
   enddo 

   call finished_adding_domain(dom_id)

end subroutine fill_domain_structure_for_spec

!----------------------------------------------------------------------

subroutine test_sizes_domain(dom_id, exp_dom_size, exp_num_dom_vars, exp_unlim_dimid)
   integer, intent(in) :: dom_id
   integer, intent(in) :: exp_dom_size
   integer, intent(in) :: exp_num_dom_vars
   integer, intent(in) :: exp_unlim_dimid


   string1 = write_message_dom(dom_id, 'get_domain_size')
   intVal = get_domain_size(dom_id)
   call assert_equal(intVal, exp_dom_size, string1)

   string1 = write_message_dom(dom_id, 'get_num_variables')
   intVal = get_num_variables(dom_id)
   call assert_equal(intVal, exp_num_dom_vars, string1)

   string1 = write_message_dom(dom_id, 'get_unlimited_dimid')
   intVal = get_unlimited_dimid(dom_id)
   call assert_equal(intVal, exp_unlim_dimid, string1)

end subroutine test_sizes_domain

!----------------------------------------------------------------------

subroutine test_variable_info(dom_id, var_id, exp_var_size, exp_num_dims, exp_kind_indx, exp_kind_string, exp_var_name)
   integer, intent(in) :: dom_id
   integer, intent(in) :: var_id
   integer, intent(in) :: exp_var_size
   integer, intent(in) :: exp_num_dims
   integer, intent(in) :: exp_kind_indx
   character(len=*), intent(in) :: exp_kind_string
   character(len=*), intent(in) :: exp_var_name
   
   string1 = write_message_var(dom_id, var_id, 'get_variable_size ')
   intVal = get_variable_size(dom_id,var_id)
   call assert_equal(intVal, exp_var_size, string1)
   
   string1 = write_message_var(dom_id, var_id, 'get_variable_name ')
   varName = get_variable_name(dom_id,var_id)
   call assert_equal(varName, exp_var_name, string1)
   
   string1 = write_message_var(dom_id, var_id, 'get_num_dims ')
   intVal = get_num_dims(dom_id,var_id)
   call assert_equal(intVal, exp_num_dims, string1)

   string1 = write_message_var(dom_id, var_id, 'get_kind_string ')
   kindString = get_kind_string(dom_id,var_id)
   call assert_equal(kindString, exp_kind_string, string1)

   string1 = write_message_var(dom_id, var_id, 'get_kind_index ')
   intVal = get_kind_index(dom_id,var_id)
   call assert_equal(intVal, exp_kind_indx, string1)

end subroutine test_variable_info

!----------------------------------------------------------------------

subroutine test_dimension_info(dom_id, var_id, dim_ids, exp_io_numdims, exp_dim_names, exp_dim_lengths,  io_dim_ids, exp_io_dim_names, exp_io_dim_lengths, exp_io_dim_ids, exp_io_unique_dim_ids, exp_io_unique_numdims, exp_io_unique_dim_length, exp_io_unique_dim_names)
   integer,          intent(in) :: dom_id
   integer,          intent(in) :: var_id
   integer,          intent(in) :: dim_ids(:)
   character(len=*), intent(in) :: exp_dim_names(:)
   integer,          intent(in) :: exp_dim_lengths(:)
   integer,          intent(in) :: io_dim_ids(:)
   integer,          intent(in) :: exp_io_numdims
   character(len=*), intent(in) :: exp_io_dim_names(:)
   integer,          intent(in) :: exp_io_dim_lengths(:)
   integer,          intent(in) :: exp_io_dim_ids(:)
   integer,          intent(in) :: exp_io_unique_dim_ids(:)
   integer,          intent(in) :: exp_io_unique_numdims
   integer,          intent(in) :: exp_io_unique_dim_length(:)
   character(len=*), intent(in) :: exp_io_unique_dim_names(:)

   integer, allocatable :: int_array(:)
   integer :: i, num_dims

   ! user interface for dimensions

   num_dims = get_num_dims(dom_id,var_id)
   do i = 1,num_dims
      string1 = write_message_dim(dom_id, var_id, dim_ids(i), 'get_dim_name ')
      varName = get_dim_name(dom_id,var_id,dim_ids(i))
      call assert_equal(varName, exp_dim_names(i), string1)

      string1 = write_message_dim(dom_id, var_id, dim_ids(i), 'get_dim_length ')
      intVal = get_dim_length(dom_id,var_id,dim_ids(i))
      call assert_equal(intVal, exp_dim_lengths(i), string1)
   
   enddo
   
   allocate( int_array(num_dims) )
   
   string1 = write_message_var(dom_id, var_id, 'get_dim_lengths ')
   int_array = get_dim_lengths(dom_id,var_id)
   call assert_equal(int_array, exp_dim_lengths, string1)

   deallocate( int_array )


   string1 = write_message_var(dom_id, var_id, 'get_io_num_dims')
   num_dims = get_io_num_dims(dom_id, var_id)
   call assert_equal(num_dims, exp_io_numdims, string1)

   allocate( int_array(num_dims) )

   string1 = write_message_var(dom_id, var_id, 'get_io_dim_lengths ')
   int_array = get_io_dim_lengths(dom_id,var_id)
   call assert_equal(int_array, exp_io_dim_lengths, string1)

   deallocate( int_array )

   ! IO interface for dimensions
   num_dims = get_io_num_unique_dims(dom_id)
   do i = 1,num_dims
      string1 = write_message_dim(dom_id, var_id, exp_io_unique_dim_ids(i), 'get_io_unique_dim_name ')
      varName = get_io_unique_dim_name(dom_id,exp_io_unique_dim_ids(i))
      call assert_equal(varName, exp_io_unique_dim_names(i), string1)

      string1 = write_message_dim(dom_id, var_id, exp_io_unique_dim_ids(i), 'get_io_unique_dim_lengths ')
      intVal = get_io_unique_dim_length(dom_id,exp_io_unique_dim_ids(i))
      call assert_equal(intVal, exp_io_unique_dim_length(i), string1)
   
   enddo
   
   allocate( int_array(num_dims) )

   string1 = write_message_var(dom_id, var_id, 'get_io_dim_ids ')
   int_array = get_io_dim_ids(dom_id, var_id)
   call assert_equal(int_array, exp_io_dim_ids, string1)

   string1 = write_message_var(dom_id, var_id, 'get_io_num_unique_dims')
   intVal = get_io_num_unique_dims(dom_id)
   call assert_equal(intVal, exp_io_unique_numdims, string1)

   deallocate( int_array )

end subroutine test_dimension_info

!----------------------------------------------------------------------

subroutine test_update_and_clamping(dom_id, var_id, exp_clamp_min, exp_clamp_max, exp_clamping, exp_update)
   integer,  intent(in) :: dom_id
   integer,  intent(in) :: var_id
   real(r8), intent(in) :: exp_clamp_min
   real(r8), intent(in) :: exp_clamp_max
   logical,  intent(in) :: exp_clamping
   logical,  intent(in) :: exp_update

   string1 = write_message_var(dom_id, var_id, 'get_io_clamping_minval')
   realVal = get_io_clamping_minval(dom_id,var_id)
   call assert_equal(realVal, exp_clamp_min, string1)
   
   string1 = write_message_var(dom_id, var_id, 'get_io_clamping_maxval')
   realVal = get_io_clamping_maxval(dom_id,var_id)
   call assert_equal(realVal, exp_clamp_max, string1)
   
   string1 = write_message_var(dom_id, var_id, 'do_io_clamping')
   trueFalse = do_io_clamping(dom_id,var_id)
   call assert_equal(trueFalse, exp_clamping, string1)
   
   string1 = write_message_var(dom_id, var_id, 'do_io_update')
   trueFalse = do_io_update(dom_id,var_id)
   call assert_equal(trueFalse, exp_update, string1)

end subroutine test_update_and_clamping

!----------------------------------------------------------------------

subroutine test_indices(dom_id, var_id, var_name, exp_start, exp_end)
   integer,          intent(in) :: dom_id
   integer,          intent(in) :: var_id
   character(len=*), intent(in) :: var_name
   integer(i8),      intent(in) :: exp_start
   integer(i8),      intent(in) :: exp_end

   string1 = write_message_var(dom_id, var_id, 'get_index_start_from_varname')
   int8val = get_index_start(dom_id,var_name)
   call assert_equal(int8val, exp_start, string1)
   string1 = write_message_var(dom_id, var_id, 'get_index_start_from_varid')
   int8val = get_index_start(dom_id,var_id)
   call assert_equal(int8val, exp_start, string1)


   string1 = write_message_var(dom_id, var_id, 'get_index_end_from_varname')
   int8val = get_index_end(dom_id,var_name)
   call assert_equal(int8val, exp_end, string1)

   string1 = write_message_var(dom_id, var_id, 'get_index_end_from_varid')
   int8val = get_index_end(dom_id,var_id)
   call assert_equal(int8val, exp_end, string1)

end subroutine test_indices

!----------------------------------------------------------------------

subroutine test_sum_variables_below(dom_id, var_id1, var_id2, exp_sum_range, exp_sum_below1, exp_sum_below2)
   integer,     intent(in) :: dom_id
   integer,     intent(in) :: var_id1
   integer,     intent(in) :: var_id2
   integer(i8), intent(in) :: exp_sum_range
   integer(i8), intent(in) :: exp_sum_below1
   integer(i8), intent(in) :: exp_sum_below2

   string1 = write_message_2var(dom_id, var_id1, var_id2, 'get_sum_variables')
   int8Val = get_sum_variables(var_id1,var_id2,dom_id)
   call assert_equal(int8Val, exp_sum_range, string1)
   
   string1 = write_message_var(dom_id, var_id1, 'get_sum_variables_below')
   int8Val = get_sum_variables_below(var_id1,dom_id)
   call assert_equal(int8Val, exp_sum_below1, string1)

   string1 = write_message_var(dom_id, var_id2, 'get_sum_variables_below')
   int8Val = get_sum_variables_below(var_id2,dom_id)
   call assert_equal(int8Val, exp_sum_below2, string1)

end subroutine test_sum_variables_below

!----------------------------------------------------------------------

subroutine test_state_indicies(index_in, exp_iloc, exp_jloc, exp_kloc, exp_varid, exp_domid, exp_kind_index, exp_kind_string)
   integer(i8),      intent(in) :: index_in
   integer,          intent(in) :: exp_iloc
   integer,          intent(in) :: exp_jloc
   integer,          intent(in) :: exp_kloc
   integer,          intent(in) :: exp_varid
   integer,          intent(in) :: exp_domid
   integer,          intent(in) :: exp_kind_index
   character(len=*), intent(in) :: exp_kind_string
   
   integer :: iloc, jloc, kloc, var_id, dom_id, kind_index
   character(len=32) :: kind_string

   call get_model_variable_indices(index_in, iloc, jloc, kloc, var_id, &
                                   dom_id, kind_index, kind_string)

   string1 = write_message_indx(index_in, 'iloc',       'get_model_variable_indices')
   call assert_equal(iloc,       exp_iloc,         string1)
   string1 = write_message_indx(index_in, 'jloc',       'get_model_variable_indices')
   call assert_equal(jloc,       exp_jloc,         string1)
   string1 = write_message_indx(index_in, 'kloc',       'get_model_variable_indices')
   call assert_equal(kloc,       exp_kloc,         string1)
   string1 = write_message_indx(index_in, 'var_id',     'get_model_variable_indices')
   call assert_equal(var_id,     exp_varid,        string1)
   string1 = write_message_indx(index_in, 'dom_id',     'get_model_variable_indices')
   call assert_equal(dom_id,     exp_domid,        string1)
   string1 = write_message_indx(index_in, 'kind_index', 'get_model_variable_indices')
   call assert_equal(kind_index, exp_kind_index,   string1)
   string1 = write_message_indx(index_in, 'kind_string', 'get_model_variable_indices')
   call assert_equal(kind_string,exp_kind_string , string1)
   
   string1 = write_message_indx(index_in, 'index_out',    'get_dart_vector_index')
   int8Val = get_dart_vector_index(iloc, jloc, kloc, dom_id, var_id)
   call assert_equal(int8Val, index_in, string1)

end subroutine test_state_indicies

!----------------------------------------------------------------------

subroutine test_varids_from_kind(dom_id, kind_index, exp_num_varids_from_kind, exp_varids_from_kind)
   integer, intent(in) :: dom_id
   integer, intent(in) :: kind_index
   integer, intent(in) :: exp_num_varids_from_kind
   integer, intent(in) :: exp_varids_from_kind(:)
   
   integer :: size1, nVarIds
   integer, allocatable :: int_array(:)

   string1 = write_message_kind(dom_id, kind_index, 'get_num_varids_from_kind')
   nVarIds = get_num_varids_from_kind(dom_id, kind_index)
   call assert_equal(nVarIds, exp_num_varids_from_kind , string1)
   
   size1 = size(exp_varids_from_kind,1)
   allocate( int_array(size1) )

   string1 = write_message_kind(dom_id, kind_index, 'get_num_varids_from_kind')
   call get_varids_from_kind(dom_id, kind_index, int_array)
   call assert_equal(int_array, exp_varids_from_kind, string1)
   
   if ( nVarIds == 1 ) then ! only works for case when one to one maping of kinds and variables
      string1 = write_message_kind(dom_id, kind_index, 'get_varid_from_kind')
      intVal  = get_varid_from_kind(dom_id, kind_index)
      call assert_equal(intVal, exp_varids_from_kind(1), string1)
   endif

   deallocate( int_array )

end subroutine test_varids_from_kind

!----------------------------------------------------------------------

function find_domain(state_indx) result(dom)
   integer(i8), intent(in) :: state_indx
   integer :: dom

   integer :: idom
   
   dom = -1
   do idom = 1, get_num_domains()
      if ( state_indx >= msize(idom) .and. &
           state_indx <= msize(idom+1) ) then
         dom = idom
         return
      endif
   enddo
         
end function find_domain

!----------------------------------------------------------------------

function write_message_dom(dom_id, routine) result(msg)
   integer,          intent(in) :: dom_id
   character(len=*), intent(in) :: routine
   
   character(len=512) :: msg

   write(msg, '(''dom'',I1,'':'', A)') dom_id, trim(routine) 

end function write_message_dom

!----------------------------------------------------------------------

function write_message_2var(dom_id, var_id1, var_id2, routine) result(msg)
   integer,          intent(in) :: dom_id
   integer,          intent(in) :: var_id1
   integer,          intent(in) :: var_id2
   character(len=*), intent(in) :: routine
   
   character(len=512) :: msg

   write(msg, '(''dom'',I1,'':var'',I0.2,'':var'',I0.2,'':'',A)') dom_id, var_id1, var_id2, trim(routine)

end function write_message_2var

!----------------------------------------------------------------------

function write_message_var(dom_id, var_id, routine) result(msg)
   integer,          intent(in) :: dom_id
   integer,          intent(in) :: var_id
   character(len=*), intent(in) :: routine
   
   character(len=512) :: msg

   write(msg, '(''dom'',I1,'':var'',I0.2,'':'',A)') dom_id, var_id, trim(routine)

end function write_message_var

!----------------------------------------------------------------------

function write_message_kind(dom_id, kind_index, routine) result(msg)
   integer,          intent(in) :: dom_id
   integer,          intent(in) :: kind_index
   character(len=*), intent(in) :: routine
   
   character(len=512) :: msg
   character(len=32)  :: kind_string

   kind_string = get_name_for_quantity(kind_index)

   write(msg, '(''dom'',I1,'':kind:'',A,'':'',A)') dom_id, kind_string, trim(routine)

end function write_message_kind

!----------------------------------------------------------------------

function write_message_dim(dom_id, var_id, dim_id, routine) result(msg)
   integer,          intent(in) :: dom_id
   integer,          intent(in) :: var_id
   integer,          intent(in) :: dim_id
   character(len=*), intent(in) :: routine
   
   character(len=512) :: msg

   write(msg, '(''dom'',I1,'':var'',I0.2'':dim'',I0.2,'':'',A)') dom_id, var_id, dim_id, trim(routine) 

end function write_message_dim

!----------------------------------------------------------------------

function write_message_indx(index_in, test, routine) result(msg)
   integer(I8),      intent(in) :: index_in
   character(len=*), intent(in) :: test
   character(len=*), intent(in) :: routine
   
   character(len=512) :: msg

   write(msg, '(I4,A,'':'',A)') index_in, trim(test), trim(routine) 

end function write_message_indx

!----------------------------------------------------------------------

subroutine initialize_module()

  call initialize_utilities('test_state_structure')
  call register_module(source, revision, revdate)

end subroutine initialize_module


end program

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
