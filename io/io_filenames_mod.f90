! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module io_filenames_mod

!> \defgroup io_filenames io_filenames
!> Aim is to store the io filenames
!>  * Restarts
!>  * Diagnostics
!>  * Inflation files
!>
!> Any module can set the filenames here, then state_vector_io_mod
!> can read from this module to get the filenames.
!> Maybe this is a bit lazy, but I just want a way for the 
!> different modules to set the filenames
!>
!> Diagnostic files could have different netcdf variable ids
!> @{

use utilities_mod,        only : do_nml_file, nmlfileunit, do_nml_term, &
                                 check_namelist_read, find_namelist_in_file, &
                                 file_exist, E_ERR, error_handler, nc_check
use model_mod,            only : construct_file_name_in
use state_structure_mod,  only : get_num_domains, get_dim_length, get_dim_name, &
                                 get_num_dims, get_num_variables, get_variable_name
use ensemble_manager_mod, only : is_single_restart_file_in
use ensemble_manager_mod, only : is_single_restart_file_in
use netcdf

implicit none

private

! These should probably be set and get functions rather than 
! direct access

public :: io_filenames_init, restart_files_in, restart_files_out, &
          set_filenames

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

! How do people name there restart files?
! What about domains?

! public arrays of filenames. Do we need arrays for restarts AND extras?
character(len=2048), allocatable :: restart_files_in(:,:), restart_files_out(:,:,:)

! Namelist options
character(len=512) :: restart_in_stub  = 'input'
character(len=512) :: restart_out_stub = 'output'
logical            :: overwrite_input = .false.  ! sets the output file = input file


! Should probably get num_domains, num_restarts from elsewhere. In here for now
namelist / io_filenames_nml / restart_in_stub, restart_out_stub, overwrite_input

contains

!----------------------------------
!> read namelist and set up filename arrays
subroutine io_filenames_init()

integer :: iunit, io
!call register_module(source, revision, revdate)

! Read the namelist entry
call find_namelist_in_file("input.nml", "io_filenames_nml", iunit)
read(iunit, nml = io_filenames_nml, iostat = io)
call check_namelist_read(iunit, io, "io_filenames_nml")

! Write the namelist values to the log file
if (do_nml_file()) write(nmlfileunit, nml=io_filenames_nml)
if (do_nml_term()) write(     *     , nml=io_filenames_nml)

end subroutine io_filenames_init

!----------------------------------

subroutine set_filenames(ens_size, inflation_in, inflation_out)

integer,              intent(in) :: ens_size
character(len = 129), intent(in) :: inflation_in(2), inflation_out(2)
character(len = 4)   :: extension, dom_str

integer :: num_domains
integer :: dom, i ! loop variables
integer :: num_files

num_files = ens_size + 10 !> @toto do you worry about spare copies?

num_domains = get_num_domains()

allocate(restart_files_in(num_files, num_domains))
allocate(restart_files_out(num_files, num_domains, 2)) ! for prior and posterior filenames

do dom = 1, num_domains

   if (is_single_restart_file_in()) then ! reading first restart for now

      restart_files_in(:, dom) = construct_file_name_in(restart_in_stub, dom, 1)

      do i = 1, ens_size  ! restarts
         write(extension, '(i4.4)') i
         write(dom_str, '(i1.1)') dom
         restart_files_out(i, dom, 1) = 'prior_member_d0' // trim(dom_str) // '.' // extension
         if (overwrite_input) then
            restart_files_out(i, dom, 2) = restart_files_in(i, dom)
         else
            restart_files_out(i, dom, 2) = construct_file_name_out(restart_out_stub, dom, i)
         endif
      enddo

   else

      do i = 1, ens_size  ! restarts
         restart_files_in(i, dom)  = construct_file_name_in(restart_in_stub, dom, i)
         write(extension, '(i4.4)') i
         write(dom_str, '(i1.1)') dom
         restart_files_out(i, dom, 1) = 'prior_member_d0' // trim(dom_str) // '.' // extension
         if (overwrite_input) then
            restart_files_out(i, dom, 2) = restart_files_in(i, dom)
         else
            restart_files_out(i, dom, 2) = construct_file_name_out(restart_out_stub, dom, i)
         endif
      enddo

   endif

enddo

! input extras
do dom = 1, num_domains
   ! mean -never used
   write(restart_files_in(ens_size + 1, dom), '(A, i2.2, A)') 'mean_d', dom
   ! sd -never used
   write(restart_files_in(ens_size + 2, dom), '(A, i2.2, A)') 'sd_d',   dom
   ! prior inf copy
   write(restart_files_in(ens_size + 3, dom), '(A, A, i2.2, A)') trim(inflation_in(1)), '_mean_d', dom
   ! prior inf sd copy
   write(restart_files_in(ens_size + 4, dom), '(A, A, i2.2, A)') trim(inflation_in(1)), '_sd_d', dom
   ! post inf copy
   write(restart_files_in(ens_size + 5, dom), '(A, A, i2.2, A)') trim(inflation_in(2)), '_mean_d', dom
   ! post inf sd copy
   write(restart_files_in(ens_size + 6, dom), '(A, A, i2.2, A)') trim(inflation_in(2)), '_sd_d', dom
enddo

! output extras
do dom = 1, num_domains
   ! Prior
   ! mean
   write(restart_files_out(ens_size + 1, dom, 1), '(A, i2.2, A)') 'PriorDiag_mean_d', dom, '.nc'
   ! sd
   write(restart_files_out(ens_size + 2, dom, 1), '(A, i2.2, A)') 'PriorDiag_sd_d', dom, '.nc'
   ! prior inf copy
   write(restart_files_out(ens_size + 3, dom, 1), '(A, i2.2, A)') 'PriorDiag_inf_mean_d', dom, '.nc'
   ! prior inf sd copy
   write(restart_files_out(ens_size + 4, dom, 1), '(A, i2.2, A)') 'PriorDiag_inf_sd_d', dom, '.nc'
   ! post inf copy - not used
   write(restart_files_out(ens_size + 5, dom, 1), '(A, A, i2.2, A)') trim(inflation_out(2)), '_mean_d', dom, '.nc'
   ! post inf sd copy - not used
   write(restart_files_out(ens_size + 6, dom, 1), '(A, A, i2.2, A)') trim(inflation_out(2)), '_sd_d', dom, '.nc'

   ! Posterior
   ! mean
   write(restart_files_out(ens_size + 1, dom, 2), '(A, i2.2, A)') 'mean_d', dom, '.nc'
   ! sd
   write(restart_files_out(ens_size + 2, dom, 2), '(A, i2.2, A)') 'sd_d', dom, '.nc'
   ! prior inf copy
   write(restart_files_out(ens_size + 3, dom, 2), '(A, A, i2.2, A)') trim(inflation_out(1)), '_mean_d', dom
   ! prior inf sd copy
   write(restart_files_out(ens_size + 4, dom, 2), '(A, A, i2.2, A)') trim(inflation_out(1)), '_sd_d', dom
   ! post inf copy
   write(restart_files_out(ens_size + 5, dom, 2), '(A, A, i2.2, A)') trim(inflation_out(2)), '_mean_d', dom
   ! post inf sd copy
   write(restart_files_out(ens_size + 6, dom, 2), '(A, A, i2.2, A)') trim(inflation_out(2)), '_sd_d', dom

   ! Storage for copies that would have gone in the Prior_diag.nc if we were to write it
   write(restart_files_out(ens_size + 7, dom, 2), '(A, i2.2, A)') restart_files_out(ens_size + 1, dom, 1)
   write(restart_files_out(ens_size + 8, dom, 2), '(A, i2.2, A)') restart_files_out(ens_size + 2, dom, 1)
   write(restart_files_out(ens_size + 9, dom, 2), '(A, i2.2, A)') restart_files_out(ens_size + 3, dom, 1)
   write(restart_files_out(ens_size + 10, dom, 2), '(A, i2.2, A)') restart_files_out(ens_size + 4, dom, 1)


enddo


! check that the netcdf files match the variables for this domain
! to prevent overwritting unwanted files.
do i = 1, ens_size+6
   do dom = 1, num_domains
      ! check the prior files
      if(file_exist(restart_files_out(i,dom,1))) &
         call check_correct_variables(restart_files_out(i,dom,1),dom)

      ! check the posterior files
      if(file_exist(restart_files_out(i,dom,2))) &
         call check_correct_variables(restart_files_out(i,dom,2),dom)
   enddo
enddo

end subroutine set_filenames

!-------------------------------------------------------
!> Check that the netcdf file matches the variables
!> for this domain
!> Do you want to overload this to take a filename or 
!> netcdf file id?
!> Do we need an nc_check warning rather than error out?
!> This checks that an existing output netcdf file contains:
!>     - each variable (matched by name)
!>     - correct dimensions for each variable (matched by name and size)
subroutine check_correct_variables(netcdf_filename, dom)

character(len=*), intent(in) :: netcdf_filename
integer, intent(in) :: dom

integer :: ncfile ! netcdf file id
integer :: i ! loop index variable
integer :: j ! loop index dimension
integer :: ret ! nc_check return value

integer :: var_id ! variable id
integer :: ndims ! number of dimensions
integer, dimension(NF90_MAX_VAR_DIMS) :: dimids ! dimension ids for a variable
character(len=NF90_MAX_NAME), dimension(NF90_MAX_VAR_DIMS) :: name ! dimension names for a variables
integer, dimension(NF90_MAX_VAR_DIMS) :: length
integer :: xtype ! do we care about this? Yes.
character(len=512) :: msgstring ! message handler

ret = nf90_open(netcdf_filename, NF90_NOWRITE, ncfile)
call nc_check(ret, 'check_correct_variables opening ', netcdf_filename)

do i = 1, get_num_variables(dom)

   ! get variable id from necfile
   ret = nf90_inq_varid(ncfile, get_variable_name(dom,i), var_id)
   write(msgstring,*) 'no match for variable ',  trim(get_variable_name(dom,i)), &
                      ' in ', trim(netcdf_filename)
   call nc_check(ret, 'check_correct_variables', msgstring)

   ! get dimension information from ncfile
   ret = nf90_inquire_variable(ncfile, var_id, ndims=ndims, dimids=dimids, xtype=xtype)
   call nc_check(ret, 'check_correct_variables', 'nf90_inquire_variable')

   ! check number of dimensions are the same - should you worry about the unlimited dimension?
   if (ndims /= get_num_dims(dom,i)) then
      write(msgstring,*) 'ndims ', get_num_dims(dom,i), ' in state does not', &
                         ' match ndims ', ndims, ' in ', trim(netcdf_filename)
      call error_handler(E_ERR, 'check_correct_variables', msgstring)
   endif

   ! check if the dimensions are what we expect. The dimensions should be same size same order.
   do j = 1, get_num_dims(dom,i)

      ! get dimension names and lengths from ncfile
      ret = nf90_inquire_dimension(ncfile, dimids(j), name=name(j), len=length(j))
      call nc_check(ret, 'check_correct_variables', 'nf90_inquire_dimension')

      ! check that the dimension names are the same
      if (get_dim_name(dom,i,j) /= name(j)) then
         write(msgstring,*) 'dim name', trim(get_dim_name(dom,i,j)), ' in state does', &
                            ' not match dim name', name(j), ' in ', trim(netcdf_filename)
         call error_handler(E_ERR, 'check_correct_variables', msgstring)
      endif

      ! check that the dimension lengths are the same
      if (get_dim_length(dom,i,j) /= length(j)) then
         write(msgstring,*) 'dimension ', name(j), "'s length ", &
                            get_dim_length(dom,i,j), ' in state does not match', &
                            ' dimension length ', length(j), ' in ', trim(netcdf_filename)
         call error_handler(E_ERR, 'check_correct_variables', msgstring)
      endif

   enddo

enddo

ret = nf90_close(ncfile)
call nc_check(ret, 'check_correct_variables closing', netcdf_filename)


end subroutine check_correct_variables

!--------------------------------------------------------------------
!> construct restart file name for writing
function construct_file_name_out(stub, domain, copy)

character(len=512), intent(in) :: stub
integer,            intent(in) :: domain
integer,            intent(in) :: copy
character(len=1024)            :: construct_file_name_out

write(construct_file_name_out, '(A,  A, i2.2, A, i2.2)') TRIM(stub), '_d', domain, '.', copy

end function construct_file_name_out

!----------------------------------
!> Destroy module storage
subroutine end_io_filenames()

deallocate(restart_files_in, restart_files_out)

end subroutine end_io_filenames

!----------------------------------
!> @}
end module io_filenames_mod