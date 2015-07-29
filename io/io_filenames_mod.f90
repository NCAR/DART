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
                                 file_exist, E_ERR, E_MSG, error_handler, nc_check, &
                                 open_file, find_textfile_dims
use model_mod,            only : construct_file_name_in
use state_structure_mod,  only : get_num_domains, get_dim_length, get_dim_name, &
                                 get_num_dims, get_num_variables, get_variable_name
use ensemble_manager_mod, only : is_single_restart_file_in

use netcdf

implicit none

private

! These should probably be set and get functions rather than 
! direct access
public :: io_filenames_init, set_filenames, get_input_file, get_output_file

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

! How do people name there restart files?
! What about domains?

character(len=2048), allocatable :: restart_files_in(:,:), restart_files_out(:,:,:)

!-------------------------------------------------------------------------------
! Namelist options
!-------------------------------------------------------------------------------
character(len=512) :: restart_in_stub  = 'input'
character(len=512) :: restart_out_stub = 'output'

logical            :: overwrite_input  = .false. ! sets output file = input file
logical            :: domain_extension = .false. ! add _d0X to filenames

logical            :: rpointer          = .false. ! define a list of restart files
character(len=512) :: rpointer_file(10) = 'null'  ! list of restarts
! JH should we enforce the pointer file name?
!-------------------------------------------------------------------------------

! Should probably get num_domains, num_restarts from elsewhere. In here for now
namelist / io_filenames_nml / restart_in_stub, restart_out_stub, &
overwrite_input, domain_extension, rpointer, rpointer_file

contains

!-------------------------------------------------------------------------------
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

!-------------------------------------------------------------------------------
subroutine set_filenames(ens_size, inflation_in, inflation_out)

integer,              intent(in) :: ens_size
character(len = *),   intent(in) :: inflation_in(2), inflation_out(2)

character(len = 32)   :: ext
character(len = 32)   :: dom_str = ''

integer :: num_domains
integer :: dom, i
integer :: idom, icopy ! loop variables
integer :: num_files

integer :: iunit, ios
integer :: nlines
character(len=512) :: fname
character(len=256) :: msgstring

num_files   = ens_size + 10 !> @toto do you worry about spare copies?
num_domains = get_num_domains()

allocate(restart_files_in(num_files , num_domains))
allocate(restart_files_out(num_files, num_domains, 2)) ! for prior and posterior filenames

do idom = 1, num_domains

   ! optional domain string
   if (num_domains > 1 .or. domain_extension) then
      write(dom_str, '(A, i2.2)') '_d', idom 
   else
      write(dom_str, '(A)') ''
   endif

   ! Read filenames from rpointer_file
   if(rpointer) then
      
      write(msgstring,*) "reading restarts from ", trim(rpointer_file(idom))
      call error_handler(E_MSG,'set_filenames', &
                         msgstring, source, revision, revdate)

      if ( .not. file_exist(rpointer_file(idom)) ) then
         msgstring = 'io_filenames_mod:rpointer '//trim(rpointer_file(idom))//&
                     ' not found'
         call error_handler(E_ERR,'set_filenames', &
                msgstring, source, revision, revdate)
      endif
      
      ! Check the dimensions of the pointer file
      call find_textfile_dims(trim(rpointer_file(idom)), nlines)
      if( nlines < ens_size) then
         write(msgstring,*) 'io_filenames_mod: expecting ',ens_size, &
                            'files in ', trim(rpointer_file(idom)),  &
                            'and only found ', nlines
         call error_handler(E_ERR,'set_filenames', &
                            msgstring, source, revision, revdate)
      endif 

      ! Read filenames in
      iunit = open_file(trim(rpointer_file(idom)),action = 'read')
        
      do icopy = 1, ens_size
         read(iunit,'(A)',iostat=ios) restart_files_in(icopy, idom)
      enddo
   else ! Construct restarts

      if (is_single_restart_file_in()) then ! reading first restart for now
         !> @todo should we not append a copy number to single file?
         restart_files_in(:, idom) = construct_file_name_in(restart_in_stub, idom, 1)
      else
         do icopy = 1, ens_size  ! restarts
            restart_files_in(icopy, idom) = construct_file_name_in(restart_in_stub, idom,icopy)
         enddo
      endif

   endif ! rpointer

   ! Construct the output files
   do icopy = 1, ens_size  ! output restarts
      
      write(ext, '(2A, i4.4, A)') trim(dom_str), ".", icopy, ".nc"
      write(restart_files_out(icopy, idom, 1),'(2A)') 'prior_member', ext

      if (overwrite_input) then
         restart_files_out(icopy, idom, 2) = restart_files_in(icopy, idom)
      else
         restart_files_out(icopy, idom, 2) = construct_file_name_out(restart_out_stub, icopy, idom)
      endif
   enddo

   ! input extras
   ! mean -never used
   write(restart_files_in(ens_size+1,idom), '(A)') 'xxx_mean' 
   ! sd -never used
   write(restart_files_in(ens_size+2,idom), '(A)') 'xxx_sd'  
   ! prior inf copy
   write(restart_files_in(ens_size+3,idom), '(3A)') trim(inflation_in(1)), '_mean'
   ! prior inf sd copy
   write(restart_files_in(ens_size+4,idom), '(3A)') trim(inflation_in(1)), '_sd'  
   ! post inf copy
   write(restart_files_in(ens_size+5,idom), '(3A)') trim(inflation_in(2)), '_mean'
   ! post inf sd copy
   write(restart_files_in(ens_size+6,idom), '(3A)') trim(inflation_in(2)), '_sd'  

   ! output extras
   ! Prior
   ! mean
   write(restart_files_out(ens_size+1,idom,1),'(A)') 'PriorDiag_mean'
   ! sd
   write(restart_files_out(ens_size+2,idom,1),'(A)') 'PriorDiag_sd' 
   ! prior inf copy (should be the same as trim(inflation_in(1)), '_mean', should we write this out?) 
   ! JH turn on adaptive inflation and test if files are the same.
   write(restart_files_out(ens_size+3,idom,1),'(A)') 'PriorDiag_inf_mean'
   ! prior inf sd copy (should be the same as trim(inflation_in(1)), '_sd', should we write this out?)
   ! JH turn on adaptive inflation and test if files are the same.
   write(restart_files_out(ens_size+4,idom,1),'(A)') 'PriorDiag_inf_sd' 
   ! post inf copy - not used
   write(restart_files_out(ens_size+5,idom,1),'(2A)') trim(inflation_out(2)), 'xxx_mean'
   ! post inf sd copy - not used
   write(restart_files_out(ens_size+6,idom,1),'(2A)') trim(inflation_out(2)), 'xxx_sd_d'

   ! Posterior
   ! mean
   write(restart_files_out(ens_size+1,idom,2),'(A)') 'mean'
   ! sd
   write(restart_files_out(ens_size+2,idom,2),'(A)') 'sd' 
   ! prior inf copy
   write(restart_files_out(ens_size+3,idom,2),'(2A)') trim(inflation_out(1)), '_mean'
   ! prior inf sd copy
   write(restart_files_out(ens_size+4,idom,2),'(2A)') trim(inflation_out(1)), '_sd'
   ! post inf copy
   write(restart_files_out(ens_size+5,idom,2),'(2A)') trim(inflation_out(2)), '_mean'
   ! post inf sd copy
   write(restart_files_out(ens_size+6,idom,2),'(2A)') trim(inflation_out(2)), '_sd'
   
   ! Add extension to extra files
   write(ext, '(2A)') trim(dom_str), ".nc"
   do icopy = ens_size + 1, ens_size + 6
      write(restart_files_in(icopy,idom), '(2A)') &
            trim(restart_files_in(icopy,idom)), ext

      write(restart_files_out(icopy,idom,1), '(2A)') &
            trim(restart_files_out(icopy,idom,1)), ext

      write(restart_files_out(icopy,idom,2), '(2A)') &
            trim(restart_files_out(icopy,idom,2)), ext
   enddo

   ! Storage for copies that would have gone in the Prior_diag.nc if we 
   ! were to write it
   restart_files_out(ens_size+7 ,idom,2) =  restart_files_out(ens_size+1,idom, 1)
   restart_files_out(ens_size+8 ,idom,2) =  restart_files_out(ens_size+2,idom, 1)
   restart_files_out(ens_size+9 ,idom,2) =  restart_files_out(ens_size+3,idom, 1)
   restart_files_out(ens_size+10,idom,2) =  restart_files_out(ens_size+4,idom, 1)

enddo ! domain loop

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
integer, dimension(NF90_MAX_VAR_DIMS) :: unique_dim_length
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
function construct_file_name_out(stub, copy, domain)

character(len=512), intent(in) :: stub
integer,            intent(in) :: copy
integer,            intent(in) :: domain
character(len=1024)            :: construct_file_name_out

character(len = 32)   :: ext = ''

if (get_num_domains() > 1 .or. domain_extension)then 
   write(ext, '(A, i2.2)') '_d', domain
endif

write(construct_file_name_out, '( 2A,".",i4.4,".nc")') TRIM(stub), trim(ext), copy

end function construct_file_name_out

!----------------------------------
!> Return the appropriate input file for copy and domain
function get_input_file(copy, domain)

integer, intent(in) :: copy
integer, intent(in) :: domain

character(len=2048) :: get_input_file

get_input_file = restart_files_in(copy, domain)

end function get_input_file

!----------------------------------
!> Return the appropriate output file for copy and domain
function get_output_file(copy, domain, isprior)

integer, intent(in) :: copy
integer, intent(in) :: domain
logical, intent(in) :: isprior

character(len=2048) :: get_output_file

if(isprior) then
   get_output_file = restart_files_out(copy, domain, 1)
else
   get_output_file = restart_files_out(copy, domain, 2)
endif

end function get_output_file

!----------------------------------
!> Destroy module storage
subroutine end_io_filenames()

deallocate(restart_files_in, restart_files_out)

end subroutine end_io_filenames

!----------------------------------
!> @}
end module io_filenames_mod
