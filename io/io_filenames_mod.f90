! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module io_filenames_mod

!> \defgroup io_filenames_mod io_filenames_mod
!> Filenames for state vector IO
!> Aim is to store the io filenames for state read/writes
!>  * Restarts
!>  * Diagnostics
!>  * Inflation files
!>
!> Usage:
!> A file_info_type is created with a call to file_info_type = io_filenames_init()
!>
!> The file_info_type is public and contains:
!>   * file_options_type (private)
!>   * 3 restart_names_types (private):
!>        - restart_files_in
!>        - restart_files_out_prior (see state_space_diag_mod.f90 for where this is used.)
!>        - restart_files_out
!>
!> The file_options_type contains all the read/write options that typically
!> come from the filter or perfect_model_obs namelist. These options are passed from 
!> the calling routine to io_file_namesames_init()
!> 
!> The restart_names_types contain a 2D array of filenames (num files, num_domains).
!>
!> The file_info_type is passed to the state IO routines: read_state, write_state, 
!> and filter_state_space_diagnostics (diagnostic file)
!> The internals of the file_info_type are accessed through the accessor functions
!> listed below. assert_file_info_initialized() and assert_restart_names_initialized()
!> should be used to check that the file_info_type has been initialized before
!> attempting any IO.
!>
!> Diagnostic files could have different netcdf variable ids
!> @{

use types_mod,            only : r4, r8, MISSING_R8
use utilities_mod,        only : file_exist, E_ERR, E_MSG, error_handler, &
                                 nc_check, open_file, find_textfile_dims
use model_mod,            only : construct_file_name_in
use state_structure_mod,  only : get_num_domains, get_dim_length, get_dim_name, &
                                 get_io_num_dims, get_num_variables, get_variable_name, &
                                 get_units, get_long_name, get_short_name, get_missing_value, &
                                 get_FillValue, get_xtype, get_add_offset, get_scale_factor, &
                                 get_has_missing_value
use ensemble_manager_mod, only : ensemble_type

use copies_on_off_mod,    only : ENS_MEAN_COPY, ENS_SD_COPY, &
                                 PRIOR_INF_COPY, PRIOR_INF_SD_COPY, &
                                 POST_INF_COPY, POST_INF_SD_COPY, &
                                 SPARE_PRIOR_MEAN, SPARE_PRIOR_SPREAD, &
                                 SPARE_PRIOR_INF_MEAN, SPARE_PRIOR_INF_SPREAD, &
                                 SPARE_POST_INF_MEAN, SPARE_POST_INF_SPREAD, &
                                 query_copy_present

use netcdf

implicit none

private

! File_info_type initialization and assertions.
public :: io_filenames_init, &
          end_io_filenames, &
          file_info_type, &
          assert_file_info_initialized, &
          restart_names_type, &
          assert_restart_names_initialized
! Accessor functions:
public :: get_input_file, &
          get_output_file, &
          get_file_description, &
          get_read_from_netcdf, &
          get_write_to_netcdf, &
          get_output_restart, &
          get_output_mean, &
          get_restart_out_base, &
          get_restart_in_base, &
          get_read_from_single_file, &
          get_write_to_single_file

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

! Holds all the input options to io_filenames_init
type file_options_type
   private
   logical            :: initialized       = .false.
   logical            :: netcdf_read       = .true. ! Netcdf or DART format
   logical            :: netcdf_write      = .true. ! Netcdf or DART format
   character(len=512) :: restart_in_base   = 'input'
   character(len=512) :: restart_out_base  = 'output'

   logical            :: overwrite_input   = .false. ! sets output file = input file
   logical            :: domain_extension  = .false. ! add _d0X to filenames

   logical            :: rpointer          = .false. ! define a list of restart files
   character(len=512) :: rpointer_file(10) = 'null'  ! list of restarts
   character(len=512) :: inflation_in(2)   = 'null'  ! prior, post
   character(len=512) :: inflation_out(2)  = 'null' ! prior, post

   logical :: output_restart = .false.  ! Should these be true or false?
   logical :: output_mean = .false.

   logical :: single_restart_file_in  = .false. ! all copies read from 1 file
   logical :: single_restart_file_out = .false. ! all copies written to 1 file
end type

! Holds an array of restart file names to be used with an ensemble handle
type restart_names_type
   private
   logical                          :: initialized       = .false.
   character(len=256), allocatable  :: filenames(:,:)  ! num_files x num_domains
   character(len=512),  allocatable :: file_description(:,:)  !  information about file
end type


! Fileninfo type
! Composed of four types, - these types are private
!  * One type containing the file options
!  * Three restart_names_types: in, prior_out, posterior_out.
! File info type is public so the restart names can be passed to read and write routines
! using, for example, file_info_type%restart_files_out_prior
type file_info_type
   type(file_options_type) :: options

   type(restart_names_type) :: restart_files_in
   type(restart_names_type) :: restart_files_out_prior
   type(restart_names_type) :: restart_files_out

end type

character(len=512) :: msgstring ! message handler

contains

!-------------------------------------------------------------------------------
! Accessor functions for file_info type
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Test whether file_info_type has been initialized
!> Error out if not, giving the name of the calling routine.
subroutine assert_file_info_initialized(file_info, routine_name)

type(file_info_type), intent(in) :: file_info
character(len=*),     intent(in) :: routine_name

if ( file_info%options%initialized .eqv. .false.) then
   call error_handler(E_ERR, routine_name, ':: io_filenames_init must be used to initialize file_info_type')
endif

end subroutine assert_file_info_initialized

!-------------------------------------------------------------------------------
!> Test whether file_info_type has been initialized for routines that only
!> have access the %restart_files(in/out/prior)
!> Error out if not, giving the name of the calling routine.
subroutine assert_restart_names_initialized(restart_names, routine_name)

type(restart_names_type), intent(in) :: restart_names
character(len=*),     intent(in) :: routine_name

if ( restart_names%initialized .eqv. .false.) then
   call error_handler(E_ERR, routine_name, ':: io_filenames_init must be used to initialize file_info_type')
endif

end subroutine assert_restart_names_initialized

!-------------------------------------------------------------------------------
! Accessor functions for file_info_type
!-------------------------------------------------------------------------------

function get_read_from_single_file(file_info)

type(file_info_type), intent(in) :: file_info
logical :: get_read_from_single_file

get_read_from_single_file = file_info%options%single_restart_file_in

end function get_read_from_single_file

!-------------------------------------------------------------------------------

function get_write_to_single_file(file_info)

type(file_info_type), intent(in) :: file_info
logical :: get_write_to_single_file

get_write_to_single_file = file_info%options%single_restart_file_out

end function get_write_to_single_file

!-------------------------------------------------------------------------------

function get_read_from_netcdf(file_info)

type(file_info_type), intent(in) :: file_info
logical :: get_read_from_netcdf

get_read_from_netcdf = file_info%options%netcdf_read

end function get_read_from_netcdf

!-------------------------------------------------------------------------------
function get_write_to_netcdf(file_info)

type(file_info_type), intent(in) :: file_info
logical :: get_write_to_netcdf

get_write_to_netcdf = file_info%options%netcdf_write

end function get_write_to_netcdf

!-------------------------------------------------------------------------------
function get_output_mean(file_info)

type(file_info_type), intent(in) :: file_info
logical :: get_output_mean

get_output_mean = file_info%options%output_mean

end function get_output_mean

!-------------------------------------------------------------------------------
function get_output_restart(file_info)

type(file_info_type), intent(in) :: file_info
logical :: get_output_restart

get_output_restart = file_info%options%output_restart

end function get_output_restart

!-------------------------------------------------------------------------------
function get_restart_out_base(file_info)

type(file_info_type), intent(in) :: file_info
character(len=512) :: get_restart_out_base

get_restart_out_base = file_info%options%restart_out_base

end function get_restart_out_base

!-------------------------------------------------------------------------------
function get_restart_in_base(file_info)

type(file_info_type), intent(in) :: file_info
character(len=512) :: get_restart_in_base

get_restart_in_base = file_info%options%restart_in_base

end function get_restart_in_base
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Initialize file_info type
!> The filenames come from above - this is so filter_nml or perfect_nml
!> can hold the filenames.
!> Previously the netcdf filenames were in io_filenames_mod - the user would
!> have to change this namelist between running perfect_model_obs and
!> filter.
function io_filenames_init(ens_handle, single_restart_file_in,  single_restart_file_out, &
             restart_in_base, restart_out_base, output_restart, netcdf_read, netcdf_write, &
             output_restart_mean, domain_extension, rpointer, rpointer_file, overwrite_input, &
             inflation_in, inflation_out) result(file_info)

type(ensemble_type), intent(in) :: ens_handle
logical,             intent(in) :: single_restart_file_in ! all copies read from one file
logical,             intent(in) :: single_restart_file_out ! all copies written to one file
character(len=*),    intent(in) :: restart_in_base
character(len=*),    intent(in) :: restart_out_base
logical,             intent(in) :: output_restart
logical,             intent(in) :: netcdf_read  ! Netcdf or DART format
logical,             intent(in) :: netcdf_write ! Netcdf of DART format
! Optional arguments - apply to filter only
logical,            optional, intent(in) :: output_restart_mean
logical,            optional, intent(in) :: domain_extension ! add _d0X to filenames
logical,            optional, intent(in) :: rpointer          ! define a list of restart files
character(len=*),   optional, intent(in) :: rpointer_file(10) ! list of restarts
logical,            optional, intent(in) :: overwrite_input  ! sets output file = input file
character(len = *), optional, intent(in) :: inflation_in(2), inflation_out(2)
type(file_info_type) :: file_info


! load up filenames structure
file_info%options%initialized = .true.

file_info%options%single_restart_file_in = single_restart_file_in
file_info%options%single_restart_file_out = single_restart_file_out
file_info%options%restart_in_base = restart_in_base
file_info%options%restart_out_base = restart_out_base
file_info%options%output_restart = output_restart
file_info%options%netcdf_read  = netcdf_read
file_info%options%netcdf_write = netcdf_write

if(present(output_restart_mean)) file_info%options%output_mean = output_restart_mean
if(present(domain_extension))    file_info%options%domain_extension = domain_extension
if(present(rpointer))            file_info%options%rpointer = rpointer
if(present(rpointer))            file_info%options%rpointer_file = rpointer_file
if(present(overwrite_input))     file_info%options%overwrite_input = overwrite_input
if(present(inflation_in))        file_info%options%inflation_in = inflation_in
if(present(inflation_out))       file_info%options%inflation_out = inflation_out

!> @todo Need to be clearer about single_restart_file in and
!> what to do with the filenames.
!> Perfect_model_obs is single_restart_file_in - the name list 
!> restart_in_file should be the whole name of the file
!> filter can have single_restart_file_in - all members in one file, 
!> but how does this interact with perturb_from_single_instance?
!> Should we have a file_info%options%single_file_name?
if (netcdf_read) then
   call error_handler(E_MSG, 'io_filenames_init', 'cannot read multiple files from single netcdf file')
endif

! loads up filenames and checks any existing files are the correct shape
call set_filenames(ens_handle, file_info)

file_info%restart_files_in%initialized        = .true.
file_info%restart_files_out_prior%initialized = .true.
file_info%restart_files_out%initialized       = .true.

end function io_filenames_init

!-------------------------------------------------------------------------------
!> Loads up filenames and checks any existing files are the correct shape
!> Only the owners of the ensemble members check existing restart files.
subroutine set_filenames(ens_handle, file_info)

type(ensemble_type),  intent(in) :: ens_handle
type(file_info_type),  intent(inout) :: file_info

character(len = 32)   :: dom_str = ''

integer :: num_domains
integer :: dom, i
integer :: idom, icopy ! loop variables
integer :: num_files
integer :: ens_size ! number of actual copies

integer :: iunit, ios
integer :: nlines
character(len=256) :: file_string
integer :: copy

ens_size = ens_handle%num_copies - ens_handle%num_extras
num_files   = ens_handle%num_copies
num_domains = get_num_domains()

allocate(file_info%restart_files_in%filenames(num_files , num_domains))
allocate(file_info%restart_files_out_prior%filenames(num_files , num_domains)) ! prior
allocate(file_info%restart_files_out%filenames(num_files , num_domains)) ! posterior

allocate(file_info%restart_files_out_prior%file_description(num_files , num_domains)) ! prior
allocate(file_info%restart_files_out%file_description(num_files , num_domains)) ! posterior

file_info%restart_files_in%filenames = 'null'
file_info%restart_files_out_prior%filenames = 'null'
file_info%restart_files_out%filenames = 'null'

do idom = 1, num_domains

   ! optional domain string
   if (num_domains > 1 .or. file_info%options%domain_extension) then
      write(dom_str, '(A, i2.2)') '_d', idom 
   else
      write(dom_str, '(A)') ''
   endif

   ! Read filenames from rpointer_file
   if(file_info%options%rpointer) then
      
      write(msgstring,*) "reading restarts from ", trim(file_info%options%rpointer_file(idom))
      call error_handler(E_MSG,'set_filenames', &
                         msgstring, source, revision, revdate)

      if ( .not. file_exist(file_info%options%rpointer_file(idom)) ) then
         msgstring = 'io_filenames_mod:rpointer '//trim(file_info%options%rpointer_file(idom))//&
                     ' not found'
         call error_handler(E_ERR,'set_filenames', &
                msgstring, source, revision, revdate)
      endif
      
      ! Check the dimensions of the pointer file
      call find_textfile_dims(trim(file_info%options%rpointer_file(idom)), nlines)
      if( nlines < ens_size) then
         write(msgstring,*) 'io_filenames_mod: expecting ',ens_size, &
                            'files in ', trim(file_info%options%rpointer_file(idom)),  &
                            'and only found ', nlines
         call error_handler(E_ERR,'set_filenames', &
                            msgstring, source, revision, revdate)
      endif 

      ! Read filenames in
      iunit = open_file(trim(file_info%options%rpointer_file(idom)),action = 'read')
        
      do icopy = 1, ens_size
         read(iunit,'(A)',iostat=ios) file_info%restart_files_in%filenames(icopy, idom)
      enddo
   else ! Construct restarts

      if (file_info%options%single_restart_file_in) then ! reading first restart for now
         !> @todo should we not append a copy number to single file?
         file_info%restart_files_in%filenames(:, idom) = construct_file_name_in(file_info%options%restart_in_base, idom, 1)
      else
         do icopy = 1, ens_size  ! restarts
            file_info%restart_files_in%filenames(icopy, idom) = construct_file_name_in(file_info%options%restart_in_base, idom,icopy)
         enddo
      endif

   endif ! rpointer

   ! Construct the output files:
   do icopy = 1, ens_size  ! output restarts
      
      ! prior member file names and descriptions
      write(file_string, '(''prior_member'',A, I4.4,''.nc'')') trim(dom_str), icopy
      file_info%restart_files_out_prior%filenames(icopy, idom) = file_string

      write(file_string,'(A,I2.2,A)') 'dart prior member ', icopy, trim(dom_str)
      file_info%restart_files_out_prior%file_description(icopy,idom) = file_string

      !  restart file names and descriptions
      if (file_info%options%overwrite_input) then
         file_info%restart_files_out%filenames(icopy, idom) = file_info%restart_files_in%filenames(icopy, idom)
      else
         file_info%restart_files_out%filenames(icopy, idom) = construct_file_name_out(file_info, icopy, idom)
      endif

      write(file_string,'(A,I2.2,A)') 'dart output member ', icopy, trim(dom_str)
      file_info%restart_files_out%file_description(icopy,idom) = file_string
   enddo


   ! ---- input extras -----------------
   ! mean -never used
   ! sd -never used
   ! prior inf copy
   if (query_copy_present(PRIOR_INF_COPY)) &
      write(file_info%restart_files_in%filenames(PRIOR_INF_COPY,idom),    '(2A)') trim(file_info%options%inflation_in(1)), '_mean'
   ! prior inf sd copy
   if (query_copy_present(PRIOR_INF_SD_COPY)) &
      write(file_info%restart_files_in%filenames(PRIOR_INF_SD_COPY,idom), '(2A)') trim(file_info%options%inflation_in(1)), '_sd'
   ! post inf copy
   if (query_copy_present(POST_INF_COPY)) &
      write(file_info%restart_files_in%filenames(POST_INF_COPY,idom),     '(2A)') trim(file_info%options%inflation_in(2)), '_mean'
   ! post inf sd copy
   if (query_copy_present(POST_INF_SD_COPY))  &
      write(file_info%restart_files_in%filenames(POST_INF_SD_COPY,idom),  '(2A)') trim(file_info%options%inflation_in(2)), '_sd'
   !-----------------------------------


   ! ---- output extras ---------------
   ! PRIOR
   ! mean
   if (query_copy_present(ENS_MEAN_COPY)) &
      call write_output_file_info(file_info%restart_files_out_prior, ENS_MEAN_COPY, idom, &
             'PriorDiag_mean','dart prior ensemble mean', dom_str)

   ! sd
   if (query_copy_present(ENS_SD_COPY)) &
      call write_output_file_info(file_info%restart_files_out_prior, ENS_SD_COPY, idom, &
             'PriorDiag_sd', 'dart prior ensemble sd', dom_str)

   ! prior inf copy (should be the same as trim(inflation_in(1)), '_mean', should we write this out?)
   if (query_copy_present(PRIOR_INF_COPY)) &
      call write_output_file_info(file_info%restart_files_out_prior, PRIOR_INF_COPY, idom, &
             'PriorDiag_inf_mean', 'dart prior inflation mean', dom_str)

   ! prior inf sd copy (should be the same as trim(inflation_in(1)), '_sd', should we write this out?)
   if (query_copy_present(PRIOR_INF_SD_COPY)) &
      call write_output_file_info(file_info%restart_files_out_prior, PRIOR_INF_SD_COPY, idom, &
             'PriorDiag_inf_sd', 'dart prior inflation sd', dom_str)

   ! post inf copy - not used
   ! post inf sd copy - not used

   ! POSTERIOR
   ! mean
   if (query_copy_present(ENS_MEAN_COPY)) &
      call write_output_file_info(file_info%restart_files_out, ENS_MEAN_COPY, idom, &
             'mean', 'dart posterior ensemble mean', dom_str)

   ! sd
   if (query_copy_present(ENS_SD_COPY)) &
      call write_output_file_info(file_info%restart_files_out, ENS_SD_COPY, idom, &
             'sd', 'dart posterior ensemble sd', dom_str)

   ! potentially updated prior inflation mean ( analysis )
   if (query_copy_present(PRIOR_INF_COPY)) then
      write(file_string,*) trim(file_info%options%inflation_out(1))//'_mean'
      call write_output_file_info(file_info%restart_files_out, PRIOR_INF_COPY, idom, &
             file_string, 'dart prior inflation mean', dom_str)
   endif

   ! potentially updated prior inflation sd ( analysis )
   if (query_copy_present(PRIOR_INF_SD_COPY)) then
      write(file_string,*) trim(file_info%options%inflation_out(1))//'_sd'
      call write_output_file_info(file_info%restart_files_out, PRIOR_INF_SD_COPY, idom, &
             file_string, 'dart prior inflation sd', dom_str)
   endif

   ! potentially updated posterior inflation mean ( analysis )
   if (query_copy_present(POST_INF_COPY)) then
      write(file_string,*) trim(file_info%options%inflation_out(2))//'_mean'
      call write_output_file_info(file_info%restart_files_out, POST_INF_COPY, idom, &
             file_string, 'dart posterior inflation mean', dom_str)
   endif

   ! potentially updated posterior inflation sd ( analysis )
   if (query_copy_present(POST_INF_SD_COPY)) then
      write(file_string,*) trim(file_info%options%inflation_out(2))//'_sd'
      call write_output_file_info(file_info%restart_files_out, POST_INF_SD_COPY, idom, &
             file_string, 'dart posterior inflation sd', dom_str)
   endif

   ! Filename for copies that would have gone in the Posterior_diag.nc if we
   ! were to write it
   if (query_copy_present(SPARE_POST_INF_MEAN)) &
      call write_output_file_info(file_info%restart_files_out, SPARE_POST_INF_MEAN, idom, &
             'PosteriorDiag_inf_mean', 'dart posterior inflation mean', dom_str)

   if (query_copy_present(SPARE_POST_INF_SPREAD)) &
      call write_output_file_info(file_info%restart_files_out, SPARE_POST_INF_SPREAD, idom, &
             'PosteriorDiag_inf_sd', 'dart posterior inflation sd', dom_str)

   !-----------------------------------

   ! Filenames for copies that would have gone in the Prior_diag.nc if we were to write it, and
   ! we saved writing these copies until the end of filter.
   ! Assuming that there will always be an ens_mean_copy and ens_sd_copy.
   if ( query_copy_present(SPARE_PRIOR_MEAN)) &
      file_info%restart_files_out%filenames(SPARE_PRIOR_MEAN ,idom) =  file_info%restart_files_out_prior%filenames(ENS_MEAN_COPY,idom)

   if (query_copy_present(SPARE_PRIOR_SPREAD)) &
      file_info%restart_files_out%filenames(SPARE_PRIOR_SPREAD ,idom) =  file_info%restart_files_out_prior%filenames(ENS_SD_COPY,idom)

   ! Testing for inflation copies (these may not be exist if there is no inflation)
   if (query_copy_present(SPARE_PRIOR_INF_MEAN) .and. query_copy_present(PRIOR_INF_COPY)) &
      file_info%restart_files_out%filenames(SPARE_PRIOR_INF_MEAN ,idom) =  file_info%restart_files_out_prior%filenames(PRIOR_INF_COPY,idom)

   if (query_copy_present(SPARE_PRIOR_INF_SPREAD) .and. query_copy_present(PRIOR_INF_SD_COPY)) &
      file_info%restart_files_out%filenames(SPARE_PRIOR_INF_SPREAD,idom) =  file_info%restart_files_out_prior%filenames(PRIOR_INF_SD_COPY,idom)


enddo ! domain loop

! check that the netcdf files match the variables for this domain
! to prevent overwritting unwanted files.
do i = 1, ens_handle%my_num_copies ! just have owners check
   copy = ens_handle%my_copies(i)
   do dom = 1, num_domains
      ! check the prior files
      if(file_exist(file_info%restart_files_out_prior%filenames(copy,dom))) &
         call check_correct_variables(file_info%restart_files_out_prior%filenames(copy,dom),dom)

      ! check the posterior files
      if(file_exist(file_info%restart_files_out%filenames(copy,dom))) &
         call check_correct_variables(file_info%restart_files_out%filenames(copy,dom),dom)
   enddo
enddo

end subroutine set_filenames

!-------------------------------------------------------
!> Write output file name and description
subroutine write_output_file_info(restart_type, copy_number, domain, stub, desc, domain_string)      
type(restart_names_type), intent(inout) :: restart_type
integer,                  intent(in)    :: copy_number
integer,                  intent(in)    :: domain
character(len=*),         intent(in)    :: stub
character(len=*),         intent(in)    :: desc
character(len=*),         intent(in)    :: domain_string
   
character(len=256) :: string1

write(string1,'(2A,''.nc'')') trim(stub), trim(domain_string)
restart_type%filenames(copy_number,domain) = adjustl(trim(string1))

write(string1,'(2A)') trim(desc), trim(domain_string)
restart_type%file_description(copy_number,domain) = adjustl(trim(string1))

end subroutine write_output_file_info

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
   if (ndims /= get_io_num_dims(dom,i)) then
      write(msgstring,*) 'ndims ', get_io_num_dims(dom,i), ' in state does not', &
                         ' match ndims ', ndims, ' in ', trim(netcdf_filename)
      call error_handler(E_ERR, 'check_correct_variables', msgstring)
   endif
   
   ! check that the attributes are the same as the state structure
   call check_attributes(ncFile, netcdf_filename, var_id, dom, i)

   ! check if the dimensions are what we expect. The dimensions should be same size same order.
   do j = 1, get_io_num_dims(dom,i)

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
!> check that cf-convention attributes are consistent across restarts
subroutine check_attributes(ncFile, filename, ncVarId, domid, varid)

integer,          intent(in) :: ncFile
character(len=*), intent(in) :: filename
integer,          intent(in) :: ncVarID
integer,          intent(in) :: domid
integer,          intent(in) :: varid

integer  :: spvalINT
real(r4) :: spvalR4
real(r8) :: spvalR8

call check_attributes_name(ncFile, filename, ncVarId, 'units'     , get_units     (domid, varid) )
call check_attributes_name(ncFile, filename, ncVarId, 'long_name' , get_long_name (domid, varid) )
call check_attributes_name(ncFile, filename, ncVarId, 'short_name', get_short_name(domid, varid) )

if ( get_has_missing_value(domid, varid) ) then
   select case (get_xtype(domid,varid))
      case (NF90_INT)

         call get_FillValue(domid, varid, spvalINT)
         call check_attribute_value_int(ncFile, filename, ncVarID, '_FillValue', spvalINT)
         call get_missing_value(domid, varid, spvalINT)
         call check_attribute_value_int(ncFile, filename, ncVarID, 'missing_value', spvalINT)
      case (NF90_FLOAT)

         call get_FillValue(domid, varid, spvalR4)
         call check_attribute_value_r4(ncFile, filename, ncVarID, '_FillValue', spvalR4)
         call get_missing_value(domid, varid, spvalR4)
         call check_attribute_value_r4(ncFile, filename, ncVarID, 'missing_value', spvalR4)
      case (NF90_DOUBLE)

         call get_FillValue(domid, varid, spvalR8)
         call check_attribute_value_r8(ncFile, filename, ncVarID, '_FillValue', spvalR8)
         call get_missing_value(domid, varid, spvalR8)
         call check_attribute_value_r8(ncFile, filename, ncVarID, 'missing_value', spvalR8)
      case default
         call error_handler(E_ERR, 'check_attributes', 'unknown xtype')
   end select
endif
         
!@>todo FIXME : for now we are only storing r8 offset and scale since DART is not using them
call check_attribute_value_r8(ncFile, filename, ncVarID, 'add_offset'  , get_add_offset(domid,varid))
call check_attribute_value_r8(ncFile, filename, ncVarID, 'scale_factor', get_scale_factor(domid,varid))

end subroutine check_attributes

!--------------------------------------------------------------------
!> check integer values are the same
subroutine check_attribute_value_int(ncFile, filename, ncVarID, att_string, spvalINT)

integer,          intent(in) :: ncFile
character(len=*), intent(in) :: filename
integer,          intent(in) :: ncVarID
character(len=*), intent(in) :: att_string
integer,          intent(in) :: spvalINT

integer :: ret_spvalINT

if ( nf90_get_att(ncFile, ncVarID, att_string, ret_spvalINT) == NF90_NOERR ) then
   if (spvalINT /= ret_spvalINT) then
      write(msgstring,*) ' variable attribute, ', trim(att_string), ' in state', spvalINT, &
                         ' does not match ', trim(att_string), ' ', ret_spvalINT, ' in ', trim(filename)
      call error_handler(E_ERR, 'check_attributes', msgstring)
   endif
endif

end subroutine check_attribute_value_int

!--------------------------------------------------------------------
!> check r4 values are the same
subroutine check_attribute_value_r4(ncFile, filename, ncVarID, att_string, spvalR4)

integer,          intent(in) :: ncFile
character(len=*), intent(in) :: filename
integer,          intent(in) :: ncVarID
character(len=*), intent(in) :: att_string
real(r4),         intent(in) :: spvalR4

real(r4) :: ret_spvalR4

if ( nf90_get_att(ncFile, ncVarID, att_string, ret_spvalR4) == NF90_NOERR ) then
   if (spvalR4 /= ret_spvalR4) then
      write(msgstring,*) ' variable attribute, ', trim(att_string), ' in state', spvalR4, &
                         ' does not match ', trim(att_string), ' ', ret_spvalR4, ' in ', trim(filename)
      call error_handler(E_ERR, 'check_attribute_value_r4', msgstring)
   endif
endif

end subroutine check_attribute_value_r4

!--------------------------------------------------------------------
!> check r8 values are the same
subroutine check_attribute_value_r8(ncFile, filename, ncVarID, att_string, spvalR8)

integer,          intent(in) :: ncFile
character(len=*), intent(in) :: filename
integer,          intent(in) :: ncVarID
character(len=*), intent(in) :: att_string
real(r8),         intent(in) :: spvalR8

real(r8) :: ret_spvalR8

if ( nf90_get_att(ncFile, ncVarID, att_string, ret_spvalR8) == NF90_NOERR ) then
   if (spvalR8 /= ret_spvalR8) then
      write(msgstring,*) ' variable attribute, ', trim(att_string), ' in state', spvalR8, &
                         ' does not match ', trim(att_string), ' ', ret_spvalR8, ' in ', trim(filename)
      call error_handler(E_ERR, 'check_attribute_value_r8', msgstring)
   endif
endif

end subroutine check_attribute_value_r8

!--------------------------------------------------------------------
!> check attribute name is consistent across restarts
subroutine check_attributes_name(ncFile, filename, ncVarId, att_string, comp_string)

integer,                    intent(in) :: ncFile
character(len=*),           intent(in) :: filename
integer,                    intent(in) :: ncVarID
character(len=*),           intent(in) :: att_string
character(len=*),           intent(in) :: comp_string

character(len=NF90_MAX_NAME) :: att_name

if ( nf90_get_att(ncFile, ncVarID, att_string, att_name) == NF90_NOERR ) then
   ! inflation files will all have unitless attributes while restarts may
   ! have some measurement real units
   if (comp_string /= att_name .and. trim(att_name) /= 'unitless') then
      write(msgstring,*) ' variable attribute ,', trim(att_string), ' in state : ', trim(comp_string), &
                         ', does not match ', trim(att_name), ' in ', trim(filename)
      call error_handler(E_ERR, 'check_attributes_name', msgstring)
   end if
endif

end subroutine check_attributes_name

!--------------------------------------------------------------------
!> construct restart file name for writing
function construct_file_name_out(file_info, copy, domain)

type(file_info_type), intent(in) :: file_info
integer,             intent(in) :: copy
integer,             intent(in) :: domain
character(len=256) :: construct_file_name_out

character(len = 32)   :: ext = ''

if (get_num_domains() > 1 .or. file_info%options%domain_extension)then
   write(ext, '(A, i2.2)') '_d', domain
endif

write(construct_file_name_out, '( 2A,".",i4.4,".nc")') TRIM(file_info%options%restart_out_base), trim(ext), copy

end function construct_file_name_out

!----------------------------------
!> Return the appropriate input file for copy and domain
function get_input_file(name_handle, copy, domain)

type(restart_names_type), intent(in) :: name_handle
integer,             intent(in) :: copy
integer,             intent(in) :: domain

character(len=256) :: get_input_file

get_input_file = name_handle%filenames(copy, domain)

end function get_input_file

!----------------------------------
!> Return the appropriate output file for copy and domain
function get_output_file(name_handle, copy, domain)

type(restart_names_type), intent(in) :: name_handle
integer,             intent(in) :: copy
integer,             intent(in) :: domain

character(len=256) :: get_output_file

get_output_file = name_handle%filenames(copy, domain)

end function get_output_file

!----------------------------------
!> Return whether the file is an input, output or prior files
function get_file_description(name_handle, copy, domain)

type(restart_names_type), intent(in) :: name_handle
integer,                  intent(in) :: copy
integer,                  intent(in) :: domain

character(len=512) :: get_file_description

get_file_description= name_handle%file_description(copy, domain)

end function get_file_description

!----------------------------------
!> Destroy module storage
!>@ todo should be called somewhere

subroutine end_io_filenames(file_info)

type(file_info_type), intent(inout) :: file_info

call error_handler(E_ERR,'end_io_filenames','test this routine', &
          source, revision, revdate)

deallocate(file_info%restart_files_in%filenames)
deallocate(file_info%restart_files_out_prior%filenames)
deallocate(file_info%restart_files_out%filenames)

file_info%options%initialized = .false.
file_info%restart_files_in%initialized = .false.
file_info%restart_files_out_prior%initialized = .false.
file_info%restart_files_out%initialized = .false.

end subroutine end_io_filenames

!----------------------------------
end module io_filenames_mod
!> @}
