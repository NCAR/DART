! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

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
!>   * IO file names
!>   * IO metadata such as units, copy number, description, ...
!>
!> The restart_names_types contain a 2D array of filenames (num files, num_domains).
!>
!> The file_info_type is passed to the state IO routines: read_state, write_state,
!> and filter_state_space_diagnostics (diagnostic file)
!> 
!> The internals of the file_info_type are accessed through the accessor functions
!> listed below. assert_file_info_initialized() and assert_restart_names_initialized()
!> should be used to check that the file_info_type has been initialized before
!> attempting any IO.
!>
!> Diagnostic files could have different netcdf variable ids
!> @{

use types_mod,            only : r4, r8, MISSING_R8, MAX_NUM_DOMS, digits12
use utilities_mod,        only : file_exist, E_ERR, E_MSG, E_WARN, error_handler,&
                                 open_file, close_file, find_textfile_dims, &
                                 do_output
use time_manager_mod,     only : time_type
use mpi_utilities_mod,    only : my_task_id
use state_structure_mod,  only : get_num_domains, &
                                 get_dim_length, &
                                 get_dim_name, &
                                 get_num_dims, &
                                 get_io_num_dims, &
                                 get_num_variables, &
                                 get_variable_name, &
                                 get_units, &
                                 get_long_name, &
                                 get_short_name, &
                                 get_missing_value, &
                                 get_FillValue, &
                                 get_xtype, &
                                 get_add_offset, &
                                 get_scale_factor, &
                                 get_has_missing_value, &
                                 get_has_FillValue, &
                                 do_io_update, &
                                 get_io_dim_length
use ensemble_manager_mod, only : ensemble_type
use netcdf_utilities_mod, only : nc_check

use netcdf

implicit none

interface set_io_copy_flag
   module procedure set_io_copy_flag_range
   module procedure set_io_copy_flag_single
end interface

interface set_file_metadata
   module procedure set_explicit_file_metadata
   module procedure set_stage_file_metadata
end interface

private

! File_info_type initialization and assertions.
public :: io_filenames_init, &
          io_filenames_finalize, &
          file_info_type, &
          netcdf_file_type, &
          stage_metadata_type, &
          set_file_metadata, &
          set_member_file_metadata, &
          set_io_copy_flag, &
          assert_file_info_initialized, &
          assert_restart_names_initialized, &
          file_info_dump, &
          combine_file_info, &
          check_file_info_variable_shape

! Accessor functions:
public :: get_restart_filename, &
          get_single_file, &
          get_cycling, &
          get_file_description, &
          get_copy_name, &
          get_stage_metadata, &
          single_file_initialized, &
          inherit_copy_units, &
          copy_is_clamped, &
          force_copy_back, &
          noutput_state_variables

! Query functions:
public :: query_read_copy, &
          query_write_copy, &
          query_copy_present

! Parameters
public :: READ_COPY, &
          WRITE_COPY, &
          READ_WRITE_COPY, &
          NO_IO, &
          COPY_NOT_PRESENT

character(len=*), parameter :: source = 'io_filenames_mod.f90'

! IO Parameters
integer, parameter :: NO_IO            = -1
integer, parameter :: READ_COPY        =  1
integer, parameter :: WRITE_COPY       =  2
integer, parameter :: READ_WRITE_COPY  =  3

! Used to test if a copy is not in use, e.g. the spare copies may not be in use.
integer, parameter :: COPY_NOT_PRESENT = -1

!>@todo FIXME need to make these derived types components private again, and
!>            These types should be moved into a separate NetCDF module.

!-------------------------------------------------------------------------------
!> output (netcdf) file descriptor (diagnostic file handle)
!> basically, we want to keep a local mirror of the unlimited dimension
!> coordinate variable (i.e. time) because dynamically querying it
!> causes unacceptable performance degradation over "long" integrations.
!> The logical model_mod_will_write_state_variables determines whether the model is going
!> to define and write the state variables in the diagnostic files. If .false. dart will
!> define and write the state variables.
type netcdf_file_type
   !private
   integer :: ncid                       ! the "unit" -- sorta
   integer :: Ntimes                     ! the current working length
   integer :: NtimesMAX                  ! the allocated length.
   real(digits12),  pointer :: rtimes(:) ! times -- as a 64bit (at least) real
   type(time_type), pointer :: times(:)  ! times -- as the models use
   character(len=80)        :: fname     ! filename ...
   ! The following only applies to single file
   logical :: model_mod_will_write_state_variables = .false.
   integer :: diag_id = -1  ! to access state_structure
end type netcdf_file_type

! Holds an array of restart file names to be used with an ensemble handle
!>@todo FIXME 'file_description' -> 'copy_description' -> 'file metadata'?
type stage_metadata_type
   !private
   logical                         :: initialized       = .false.
   integer                         :: noutput_ens = 0   ! num_output_state_members
   integer                         :: num_copies  = 0
   logical,            allocatable :: clamp_vars(:)     ! num_copies
   logical,            allocatable :: inherit_units(:)  ! num_copies
   logical,            allocatable :: force_copy_back(:)! num_copies
   integer,            allocatable :: io_flag(:)        ! read = 1, write = 2, read/write = 3
   integer,            allocatable :: my_copy_number(:) ! num_copies
   character(len=256), allocatable :: copy_name(:)      ! num_copies
   character(len=256), allocatable :: long_name(:)      ! num_copies
   character(len=256), allocatable :: filenames(:,:)    ! num_copies x num_domains
   character(len=512), allocatable :: file_description(:,:)  !  information about file
   type(netcdf_file_type)          :: ncFileID
end type

! Fileinfo type
! Composed of four types, - these types are private
!  * One type containing the file options
!  * Three restart_names_types: in, prior_out, posterior_out.
! File info type is public so the restart names can be passed to read and write routines
! using, for example, file_info_type%restart_files_out_prior

type file_info_type
   !private
   logical            :: initialized                = .false.
   logical            :: singlefile_initialized     = .false.
   logical            :: check_output_compatibility = .false.
   logical            :: cycling                    = .false. ! model will cycle
   logical            :: single_file                = .false. ! all copies read from 1 file
   character(len=32)  :: root_name                  = 'null'

   type(stage_metadata_type) :: stage_metadata

end type

character(len=512) :: msgstring

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

if ( .not. file_info%initialized ) then
   call error_handler(E_ERR, routine_name, &
           ':: io_filenames_init must be used to initialize file_info_type', source)
endif

end subroutine assert_file_info_initialized

!-------------------------------------------------------------------------------
!> Test whether file_info_type has been initialized for routines that only
!> have access the %restart_files(in/out/prior)
!> Error out if not, giving the name of the calling routine.


subroutine assert_restart_names_initialized(restart_names, routine_name)

type(stage_metadata_type), intent(in) :: restart_names
character(len=*),         intent(in) :: routine_name

if ( .not. restart_names%initialized ) then
   call error_handler(E_ERR, routine_name, &
           ':: io_filenames_init must be used to initialize file_info_type', source)
endif

end subroutine assert_restart_names_initialized


!-------------------------------------------------------------------------------
!> Accessor function for file_info_type


function get_cycling(file_info)

type(file_info_type), intent(in) :: file_info
logical :: get_cycling

get_cycling = file_info%single_file

end function get_cycling


!-------------------------------------------------------------------------------
!> Accessor function for file_info_type


function get_single_file(file_info)

type(file_info_type), intent(in) :: file_info
logical :: get_single_file

get_single_file = file_info%single_file

end function get_single_file


!-------------------------------------------------------------------------------
!> Initialize file_info type


subroutine io_filenames_init(file_info, ncopies, cycling, single_file, &
                             restart_files, root_name, check_output_compatibility)

type(file_info_type),       intent(out):: file_info          !< structure with expanded list of filenames
integer,                    intent(in) :: ncopies            !< number of ensemble copies
logical,                    intent(in) :: cycling            !< model will cycle
logical,                    intent(in) :: single_file        !< all copies read from one file
character(len=*), optional, intent(in) :: restart_files(:,:) !< list of restarts one for each domain
character(len=*), optional, intent(in) :: root_name          !< base if restart_files not given
logical,          optional, intent(in) :: check_output_compatibility !< ensure netCDF variables exist in output BEFORE spending a ton of core hours

integer :: ndomains, esize

file_info%single_file = single_file
file_info%cycling     = cycling

ndomains = get_num_domains()

allocate(file_info%stage_metadata%force_copy_back( ncopies))
allocate(file_info%stage_metadata%clamp_vars(      ncopies))
allocate(file_info%stage_metadata%inherit_units(   ncopies))
allocate(file_info%stage_metadata%io_flag(         ncopies))
allocate(file_info%stage_metadata%my_copy_number(  ncopies))
allocate(file_info%stage_metadata%copy_name(       ncopies))
allocate(file_info%stage_metadata%long_name(       ncopies))
allocate(file_info%stage_metadata%filenames(       ncopies , ndomains))
allocate(file_info%stage_metadata%file_description(ncopies , ndomains))

file_info%stage_metadata%force_copy_back  = .false.
file_info%stage_metadata%clamp_vars       = .false.
file_info%stage_metadata%inherit_units    = .true.
file_info%stage_metadata%io_flag          = NO_IO
file_info%stage_metadata%my_copy_number   = -1
file_info%stage_metadata%noutput_ens      = 0
file_info%stage_metadata%num_copies       = ncopies
file_info%stage_metadata%copy_name        = 'copy_name_not_set'
file_info%stage_metadata%long_name        = 'long_name_not_set'
file_info%stage_metadata%filenames        = 'null'
file_info%stage_metadata%file_description = 'null'

!>@todo FIXME JPH : Should these be required interfaces?
if(present(restart_files)) then
   ! restart files are expected to be (ens_size x dom_size) arrays
   esize = SIZE(restart_files,1)
   file_info%stage_metadata%filenames(1:esize,:) = restart_files(:,:)
endif

if(present(root_name))                  file_info%root_name                     = root_name
if(present(check_output_compatibility)) file_info%check_output_compatibility    = check_output_compatibility
file_info%initialized = .true.

end subroutine io_filenames_init


!-------------------------------------------------------------------------------
!> Check that the netcdf files variables have the correct shape
!> to prevevent overwriting unwanted files


subroutine check_file_info_variable_shape(file_info, ens_handle)

type(file_info_type),  intent(inout) :: file_info
type(ensemble_type),   intent(in)    :: ens_handle

integer :: num_domains, idom, icopy, my_copy

num_domains = get_num_domains()

! check that the netcdf files match the variables for this domain
! to prevent overwriting unwanted files.
do icopy = 1, ens_handle%my_num_copies ! just have owners check
   my_copy = ens_handle%my_copies(icopy)
   do idom = 1, num_domains
      if(file_exist(file_info%stage_metadata%filenames(my_copy,idom))) &
         call check_correct_variables(file_info%stage_metadata%filenames(my_copy,idom),idom)
   enddo
enddo

end subroutine check_file_info_variable_shape


!-------------------------------------------------------
!> read file list names


subroutine set_member_file_metadata(file_info, ens_size, my_copy_start)

type(file_info_type), intent(inout) :: file_info
integer,              intent(in)    :: ens_size
integer,              intent(in)    :: my_copy_start

character(len=256) :: fname, desc
character(len=128) :: stage_name, basename
integer :: icopy, idom
integer :: offset

offset = my_copy_start - 1

if (my_copy_start <= 0)    return

!>@todo FIXME : Is it sufficient to check if the first file exists? JPH
! construct filenames if they do not exist, using the root name as a base
if (file_info%stage_metadata%filenames(1,1) == 'null' .or. &
    file_info%stage_metadata%filenames(1,1) == '') then

   if (file_info%root_name == 'null') then
      write(msgstring,*) 'Unable to construct file names.', &
                         ' No stage name or file list given'
      call error_handler(E_ERR,'set_member_file_metadata', msgstring, source)
   endif

   ! Construct file names
   call error_handler(E_MSG,'set_member_file_metadata', &
             'no file list given for stage "'//trim(file_info%root_name)// &
             '" so using default names', source)

   stage_name = file_info%root_name
   if (file_info%single_file) then
      if (get_num_domains() > 1) then
         write(msgstring,*) 'single file input is currently only supported for 1 domain models'
         call error_handler(E_ERR,'set_member_file_metadata', msgstring, source)
      endif
      write(fname,'(2A)') trim(stage_name),'.nc'
      write(desc, '(A)') 'ensemble member single file'
      call set_explicit_file_metadata(file_info, 1, (/fname/), stage_name, desc)
   else
      do icopy = 1, ens_size
         write(basename,'(A,I4.4)')  'member_', icopy
         write(desc,'(A,I4)') 'ensemble member ', icopy
         call set_file_metadata(file_info, icopy, stage_name, basename, desc, offset)
      enddo
   endif

else ! ensemble member files have been defined, just set metadata

   do idom = 1, get_num_domains()
      do icopy = 1, ens_size

         write(desc,'(2A,I4)') trim(file_info%root_name), ' ensemble member ', icopy
         file_info%stage_metadata%file_description(offset+icopy,idom) = trim(desc)
         file_info%stage_metadata%my_copy_number(  offset+icopy)      = offset + icopy
         file_info%stage_metadata%copy_name(       offset+icopy)      = trim(desc)

      enddo
   enddo
endif

file_info%stage_metadata%initialized = .true.

end subroutine set_member_file_metadata


!-------------------------------------------------------
!> Explicitly construct a file name from a list of files

subroutine set_explicit_file_metadata(file_info, cnum, fnames, basename, desc)

type(file_info_type), intent(inout) :: file_info
integer,              intent(in)    :: cnum
character(len=*),     intent(in)    :: fnames(:)
character(len=*),     intent(in)    :: basename
character(len=*),     intent(in)    :: desc

character(len=256) :: string1
integer :: idom

if (cnum <= 0) return

! Array of files (possibly multiple domains)
file_info%stage_metadata%my_copy_number(cnum) = cnum
file_info%stage_metadata%copy_name(     cnum) = trim(basename)

! check that the number of domains matches the size of file names provided
if (get_num_domains() /= size(fnames(:),1)) then
   write(msgstring,'(A,I2,A,I2)') 'The number of domains, ', get_num_domains(), &
                                  ', does not match the number of filenames, ', size(fnames(:),1)
   call error_handler(E_ERR, 'set_explicit_file_metadata', msgstring, source)
endif

do idom = 1, get_num_domains()
    if (get_num_domains() > 1) then
        write(string1,'(2A,I4)') trim(desc), 'for domain ', idom
    else
        write(string1,*) trim(desc)
    endif
    file_info%stage_metadata%filenames(cnum,idom)        = trim(fnames(idom))
    file_info%stage_metadata%file_description(cnum,idom) = trim(string1)
enddo

file_info%stage_metadata%initialized = .true.

end subroutine set_explicit_file_metadata


!-------------------------------------------------------
!> Write file name and description

subroutine set_stage_file_metadata(file_info, cnum, stage, basename, desc, offset)

! EXAMPLE: to construct a file name of 'input_mean_d01.nc'
!   cnum is ensemble copy number and is not used in the file name
!   stage_name is 'input'
!   basename is 'mean'
!   desc is 'prior ensemble mean'

type(file_info_type), intent(inout) :: file_info
integer,              intent(in)    :: cnum
character(len=*),     intent(in)    :: stage
character(len=*),     intent(in)    :: basename
character(len=*),     intent(in)    :: desc
integer,              intent(in), optional :: offset

character(len=256) :: string1
character(len=32)  :: stage_name
character(len=32)  :: dom_str
integer :: idom, my_copy

if (cnum <= 0) return

if (present(offset)) then
   my_copy = cnum + offset
else
   my_copy = cnum
endif

if (trim(file_info%root_name) == 'null') then
   stage_name = stage
else
   stage_name = trim(file_info%root_name)
endif

if (my_copy <= 0) return

file_info%stage_metadata%my_copy_number(my_copy) = my_copy
file_info%stage_metadata%copy_name(     my_copy) = trim(basename)

if (get_num_domains() > 1) then
   do idom = 1, get_num_domains()
     write(dom_str, '(A, i2.2)') '_d', idom
     write(string1,'(2A,''.nc'')') trim(stage_name)//'_'//trim(basename), trim(dom_str)
     file_info%stage_metadata%filenames(my_copy,idom) = trim(string1)

     write(string1,'(A,1x,A,'' for domain '',i4)') trim(stage_name), trim(desc), idom
     file_info%stage_metadata%file_description(my_copy,idom) = trim(string1)
  enddo
else
  write(string1,'(A,''.nc'')') trim(stage_name)//'_'//trim(basename)
  file_info%stage_metadata%filenames(my_copy,1) = trim(string1)

  write(string1,'(A,1x,A)') trim(stage_name), trim(desc)
  file_info%stage_metadata%file_description(my_copy,1) = trim(string1)
endif

file_info%stage_metadata%initialized = .true.

end subroutine set_stage_file_metadata


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

ret = nf90_open(netcdf_filename, NF90_NOWRITE, ncfile)
call nc_check(ret, 'check_correct_variables opening ', netcdf_filename)

do i = 1, get_num_variables(dom)

   if ( .not. do_io_update(dom,i) ) cycle

   ! get variable id from necfile
   ret = nf90_inq_varid(ncfile, get_variable_name(dom,i), var_id)
   write(msgstring,*) 'no match for variable ',  trim(get_variable_name(dom,i)), &
                      ' in ', trim(netcdf_filename)
   call nc_check(ret, 'check_correct_variables', msgstring)

   ! get dimension information from ncfile
   ret = nf90_inquire_variable(ncfile, var_id, ndims=ndims, dimids=dimids)
   call nc_check(ret, 'check_correct_variables', 'nf90_inquire_variable')

   ! check number of dimensions are the same.  we are comparing
   ! the full number of dimensions including time and member if present.
   if (ndims /= get_io_num_dims(dom,i)) then
      write(msgstring,*) 'ndims ', get_io_num_dims(dom,i), ' in state does not', &
                         ' match ndims ', ndims, ' in ', trim(netcdf_filename)
      call error_handler(E_ERR, 'check_correct_variables', msgstring, source)
   endif

   ! check that the attributes are the same as the state structure
   call check_attributes(ncFile, netcdf_filename, var_id, dom, i)

   ! check if the dimensions are what we expect. The dimensions should be same size same order.
   ! in this case check only the dims for the variables, excluding time and member.
   do j = 1, get_num_dims(dom,i)

      ! get dimension names and lengths from ncfile
      ret = nf90_inquire_dimension(ncfile, dimids(j), name=name(j), len=length(j))
      call nc_check(ret, 'check_correct_variables', 'nf90_inquire_dimension')

      ! check that the dimension names are the same
      if (get_dim_name(dom,i,j) /= name(j)) then
         write(msgstring,*) 'dim name', trim(get_dim_name(dom,i,j)), ' in state does', &
                            ' not match dim name', name(j), ' in ', trim(netcdf_filename)
         call error_handler(E_ERR, 'check_correct_variables', msgstring, source)
      endif

      ! check that the dimension lengths are the same
      if (get_io_dim_length(dom,i,j) /= length(j)) then
         write(msgstring,*) 'dimension ', trim(name(j)), "'s length ", &
                            get_io_dim_length(dom,i,j), ' in state does not match', &
                            ' dimension length ', length(j), ' in ', trim(netcdf_filename)
         call error_handler(E_ERR, 'check_correct_variables', msgstring, source)
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
         call get_missing_value(domid, varid, spvalINT)
         call check_attribute_value_int(ncFile, filename, ncVarID, 'missing_value', spvalINT)

      case (NF90_FLOAT)
         call get_missing_value(domid, varid, spvalR4)
         call check_attribute_value_r4(ncFile, filename, ncVarID, 'missing_value', spvalR4)

      case (NF90_DOUBLE)
         call get_missing_value(domid, varid, spvalR8)
         call check_attribute_value_r8(ncFile, filename, ncVarID, 'missing_value', spvalR8)

      case default
         !>@todo FIXME report the variable with the unsupported xtype
         call error_handler(E_ERR, 'check_attributes', 'unknown xtype', source)
   end select
endif

if ( get_has_FillValue(   domid, varid) ) then
   select case (get_xtype(domid, varid))
      case (NF90_INT)
         call get_FillValue(domid, varid, spvalINT)
         call check_attribute_value_int(ncFile, filename, ncVarID, '_FillValue', spvalINT)

      case (NF90_FLOAT)
         call get_FillValue(domid, varid, spvalR4)
         call check_attribute_value_r4(ncFile, filename, ncVarID, '_FillValue', spvalR4)

      case (NF90_DOUBLE)
         call get_FillValue(domid, varid, spvalR8)
         call check_attribute_value_r8(ncFile, filename, ncVarID, '_FillValue', spvalR8)

      case default
         !>@todo FIXME report the variable with the unsupported xtype
         call error_handler(E_ERR, 'check_attributes', 'unknown xtype', source)
   end select
endif

!>@todo FIXME : for now we are only storing r8 offset and scale since DART is not using them
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
      call error_handler(E_ERR, 'check_attributes', msgstring, source)
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
   if (ret_spvalR4 /= ret_spvalR4) then 
      return
   endif
   if (spvalR4 /= ret_spvalR4) then
      write(msgstring,*) ' variable attribute, ', trim(att_string), ' in state', spvalR4, &
                         ' does not match ', trim(att_string), ' ', ret_spvalR4, ' in ', trim(filename)
      call error_handler(E_ERR, 'check_attribute_value_r4', msgstring, source)
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
   if (ret_spvalR8 /= ret_spvalR8) then
      return
   endif
   if (spvalR8 /= ret_spvalR8) then
      write(msgstring,*) ' variable attribute, ', trim(att_string), ' in state', spvalR8, &
                         ' does not match ', trim(att_string), ' ', ret_spvalR8, ' in ', trim(filename)
      call error_handler(E_ERR, 'check_attribute_value_r8', msgstring, source)
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
      call error_handler(E_ERR, 'check_attributes_name', msgstring, source)
   end if
endif

end subroutine check_attributes_name


!--------------------------------------------------------------------
!> construct restart file name for reading


function construct_file_names(file_info, ens_size, copy, domain)

type(file_info_type), intent(in) :: file_info
integer,             intent(in) :: ens_size
integer,             intent(in) :: copy
integer,             intent(in) :: domain
character(len=256) :: construct_file_names

character(len=32)  :: dom_str

dom_str = ''

if (get_num_domains() > 1) write(dom_str, '(A, i2.2)') '_d', domain

if (copy <= ens_size) then
   write(construct_file_names, '(A, ''_member_'', I4.4, A, ''.nc'')') trim(file_info%root_name), copy, trim(dom_str)
else
   !>@todo what to do with copy .... should be string to describe inflation/mean/etc.
!  write(construct_file_names, '( 2A,".",i4.4,".nc")') trim(file_info%root_name), trim(dom_str), copy
   write(construct_file_names, '(''uninitialized'')')
endif

end function construct_file_names


!----------------------------------
!> Return the appropriate input file for copy and domain


function get_restart_filename(name_handle, copy, domain)

type(stage_metadata_type), intent(in) :: name_handle
integer,             intent(in) :: copy
integer,             intent(in) :: domain

character(len=256) :: get_restart_filename

get_restart_filename = name_handle%filenames(copy, domain)

end function get_restart_filename


!----------------------------------
!> Return whether the file is an input, output or prior files


function get_file_description(name_handle, copy, domain)

type(stage_metadata_type), intent(in) :: name_handle
integer,                  intent(in) :: copy
integer,                  intent(in) :: domain

character(len=512) :: get_file_description

get_file_description= name_handle%file_description(copy, domain)

end function get_file_description


!----------------------------------
!> Return stage metadata from file handle


function get_stage_metadata(file_info)

type(file_info_type), intent(in) :: file_info

type(stage_metadata_type) :: get_stage_metadata

get_stage_metadata = file_info%stage_metadata

end function get_stage_metadata


!----------------------------------
!> Destroy module storage
!>@todo FIXME should be called somewhere


subroutine io_filenames_finalize(file_info)

type(file_info_type), intent(inout) :: file_info

deallocate(file_info%stage_metadata%filenames)

file_info%initialized               = .false.
file_info%stage_metadata%initialized = .false.

end subroutine io_filenames_finalize


!----------------------------------
!> routine to summarize the contents of the file_info_type
!>

subroutine file_info_dump(file_info,context)

type(file_info_type),       intent(in) :: file_info
character(len=*), optional, intent(in) :: context

integer :: i,j

! we only want task 0 to do the print statements
if (.not. do_output()) return

if (present(context)) then
   write(*,*) trim(context)
   write(*,*)'file_info%initialized                ', file_info%initialized
   write(*,*)'file_info%check_output_compatibility ', file_info%check_output_compatibility
   write(*,*)'file_info%single_file                ', file_info%single_file
   write(*,*)'file_info%root_name                  ', file_info%root_name
   write(*,*)'file_info%stage_metadata%initialized ', file_info%stage_metadata%initialized
endif

do i = 1,size(file_info%stage_metadata%filenames,1)
   do j = 1,size(file_info%stage_metadata%filenames,2)
      if (file_info%stage_metadata%filenames(i,j) /= 'null') then
         write(*,'(A,2I4,2A)')   'file_info%stage_metadata%filenames(       ',i,j,' ) ', &
                           trim(file_info%stage_metadata%filenames(           i,j))
         write(*,'(A,2I4,2A)')   'file_info%stage_metadata%file_description(',i,j,' ) ', &
                           trim(file_info%stage_metadata%file_description(    i,j))
         write(*,'(A,2I4,A,I4)') 'file_info%stage_metadata%my_copy_number(  ',i,j,' ) ', &
                                file_info%stage_metadata%my_copy_number(      i)
         write(*,'(A, I4,A,I4)') 'file_info%stage_metadata%io_flag(         ',i,'     ) ', &
                                file_info%stage_metadata%io_flag(             i)
         write(*,'(A, I4,2A)')   'file_info%stage_metadata%copy_name(       ',i,'     ) ', &
                                  file_info%stage_metadata%copy_name(i)
      endif
   enddo
enddo

end subroutine file_info_dump

!----------------------------------
!> Combine multiple file_info_type into a single file_info


function combine_file_info(file_info_array) result(file_info_out)

type(file_info_type), intent(in) :: file_info_array(:)
type(file_info_type) :: file_info_out

integer :: i, j, k, num_domains, num_files

num_domains = get_num_domains()

num_files = size(file_info_array(1)%stage_metadata%filenames(:,:),1)

call io_filenames_init(file_info_out, num_files, .false., .false., root_name='combine_files')

do i = 1, size(file_info_array(:),1)
   do j = 1, num_files
      do k = 1, num_domains
         if (trim(file_info_array(i)%stage_metadata%filenames(j, k)) /= trim('null'))then
            file_info_out%stage_metadata%filenames(                   j, k) = &
                                file_info_array(i)%stage_metadata%filenames(j, k)
            file_info_out%stage_metadata%file_description(            j, k) = &
                         file_info_array(i)%stage_metadata%file_description(j, k)
            file_info_out%stage_metadata%io_flag(                     j   ) = &
                         file_info_array(i)%stage_metadata%io_flag(j)
            file_info_out%stage_metadata%clamp_vars(                     j   ) = &
                         file_info_array(i)%stage_metadata%clamp_vars(j)
            file_info_out%stage_metadata%copy_name(                      j   ) = &
                         file_info_array(i)%stage_metadata%copy_name(j)
         endif
      enddo
   enddo
enddo

file_info_out%stage_metadata%initialized = .true.
file_info_out%initialized = .true.

end function combine_file_info

!-------------------------------------------------------
!> Set whether a copy should be read/written for a range
!> of copies c1->c2.  Optional argument to set if the
!> copy could have units and if the variables should be
!> clamped.  If this information is available it grabs
!> it from the state structure and stores it in files
!> created from scratch.


subroutine set_io_copy_flag_range(file_info, c1, c2, io_flag, num_output_ens, &
                                  inherit_units, clamp_vars, force_copy_back)

type(file_info_type),      intent(inout) :: file_info      !< stage name handle
integer,                   intent(in)    :: c1             !< start copy to read
integer,                   intent(in)    :: c2             !< end copy to read
integer,                   intent(in)    :: io_flag        !< read = 1, write = 2, read/write = 3
integer, optional,         intent(in)    :: num_output_ens !< number of output ensembles to diag files
logical, optional,         intent(in)    :: inherit_units  !< inherit units from state structure
logical, optional,         intent(in)    :: clamp_vars
logical, optional,         intent(in)    :: force_copy_back

integer :: i

if (c1 <=0 .or. c2 <=0) return

do i = c1, c2
  file_info%stage_metadata%io_flag(i) = io_flag

  if(present(inherit_units)  )    &
     file_info%stage_metadata%inherit_units( i) = inherit_units

  if(present(clamp_vars) )    &
     file_info%stage_metadata%clamp_vars(i) = clamp_vars

  if(present(force_copy_back) ) &
     file_info%stage_metadata%force_copy_back(i) = force_copy_back

  if(present(num_output_ens)) &
     file_info%stage_metadata%noutput_ens = num_output_ens
enddo

end subroutine set_io_copy_flag_range

!-------------------------------------------------------
!> Set whether a copy should be read/written for a single
!> copy cnum.  Optional argument to set if the
!> copy could have units and if the variables should be
!> clamped.  If this information is available it grabs
!> it from the state structure and stores it in files
!> created from scratch.


subroutine set_io_copy_flag_single(file_info, cnum, io_flag, inherit_units, &
                                   clamp_vars, force_copy_back)

type(file_info_type),      intent(inout) :: file_info     !< stage name handle
integer,                   intent(in)    :: cnum          !< start copy to read
integer,                   intent(in)    :: io_flag       !< read = 1, write = 2, read/write = 3
logical, optional,         intent(in)    :: inherit_units !< inherit units from state structure
logical, optional,         intent(in)    :: clamp_vars
logical, optional,         intent(in)    :: force_copy_back

if (cnum <=0) return

file_info%stage_metadata%io_flag(cnum)   = io_flag

if(present(inherit_units)  ) &
   file_info%stage_metadata%inherit_units(   cnum ) = inherit_units
if(present(clamp_vars) ) &
   file_info%stage_metadata%clamp_vars(      cnum ) = clamp_vars
if(present(force_copy_back) ) &
   file_info%stage_metadata%force_copy_back( cnum ) = force_copy_back

end subroutine set_io_copy_flag_single


!----------------------------------
!> Determine whether a copy could potentially have units.
!> Copies such as Spread and Inflation never have units.


function inherit_copy_units(name_handle, copy)
type(stage_metadata_type), intent(in) :: name_handle
integer,                   intent(in) :: copy
logical :: inherit_copy_units

if (copy <= 0) return

inherit_copy_units = name_handle%inherit_units(copy)

end function inherit_copy_units

!----------------------------------
!> Determine if a variable should be clamped or not. If
!> clapming information is stored in the state structure
!> clamping will be applied accordingly


function copy_is_clamped(name_handle, copy)
type(stage_metadata_type), intent(in) :: name_handle
integer,                   intent(in)  :: copy
logical :: copy_is_clamped

if (copy <= 0) return

copy_is_clamped = name_handle%clamp_vars(copy)

end function copy_is_clamped

!----------------------------------
!>

function force_copy_back(name_handle, copy)
type(stage_metadata_type), intent(in) :: name_handle
integer,                   intent(in)  :: copy
logical :: force_copy_back

if (copy <= 0) return

force_copy_back = name_handle%force_copy_back(copy)

end function force_copy_back

!-------------------------------------------------------
!> returns true/false depending on whether you should read this copy


function query_read_copy(name_handle, c)

type(stage_metadata_type), intent(in) :: name_handle !< stage name handle
integer,                   intent(in) :: c           !< copy number
logical :: query_read_copy

query_read_copy = .false.

if (.not. assert_valid_copy(name_handle, c)) return

if (name_handle%io_flag(c) == READ_COPY .or. &
    name_handle%io_flag(c) == READ_WRITE_COPY) then
   query_read_copy = .true.
endif

end function query_read_copy


!-------------------------------------------------------
!> returns true/false depending on whether you should write this copy


function query_write_copy(name_handle, c)

type(stage_metadata_type), intent(in) :: name_handle !< stage name handle
integer,                   intent(in) :: c           !< copy number
logical :: query_write_copy

query_write_copy = .false.

if (.not. assert_valid_copy(name_handle, c)) return

if (name_handle%io_flag(c) == WRITE_COPY .or. &
    name_handle%io_flag(c) == READ_WRITE_COPY) then
   query_write_copy = .true.
endif

end function query_write_copy


!------------------------------------------------------------------
!> Test whether a copy is part of the ensemble


function query_copy_present(copy)

integer, intent(in) :: copy
logical :: query_copy_present

if (copy == COPY_NOT_PRESENT .or. copy <= 0 ) then
   query_copy_present = .false.
else
   query_copy_present = .true.
endif

end function


!------------------------------------------------------------------
!> set netcdf file type


function get_copy_name(file_handle, cnum)

type(file_info_type), intent(in) :: file_handle
integer,              intent(in) :: cnum
character(len=256) :: get_copy_name

get_copy_name = file_handle%stage_metadata%copy_name(cnum)

end function get_copy_name


!#! !------------------------------------------------------------------
!#! !> set netcdf file type
!#!
!#!
!#! subroutine nc_set_netcdf_info(file_handle, ncFileInfo)
!#!
!#! type(file_info_type),   intent(inout) :: file_handle
!#! type(netcdf_file_type), intent(in)    :: ncFileInfo
!#!
!#! file_handle%stage_metadata%ncFileID = ncFileInfo
!#!
!#! end subroutine nc_set_netcdf_info
!#!
!#!
!#!
!#!
!#! !------------------------------------------------------------------
!#! !> get diagnostic id
!#!
!#!
!#! subroutine nc_set_diag_id(ncFileInfo, diagnostic_id)
!#!
!#! type(netcdf_file_type), intent(inout) :: ncFileInfo
!#! integer,                intent(in)    :: diagnostic_id
!#!
!#! ncFileInfo%diag_id = diagnostic_id
!#!
!#! end subroutine nc_set_diag_id
!#!
!#!
!#! !------------------------------------------------------------------
!#! !> get diagnostic id
!#!
!#!
!#! function nc_get_diag_id(ncFileInfo)
!#!
!#! type(netcdf_file_type), intent(in)    :: ncFileInfo
!#! integer :: nc_get_diag_id
!#!
!#! nc_get_diag_id = ncFileInfo%diag_id
!#!
!#! end function nc_get_diag_id
!#!
!#!
!#! !------------------------------------------------------------------
!#! !> set netcdf file name
!#!
!#!
!#! subroutine nc_set_fname(ncFileInfo, filename)
!#!
!#! type(netcdf_file_type), intent(inout) :: ncFileInfo
!#! character(len=*) :: filename
!#!
!#! ncFileInfo%fname = filename
!#!
!#! end subroutine nc_set_fname
!#!
!#!
!#! !------------------------------------------------------------------
!#! !> get netcdf file name
!#!
!#!
!#! function nc_get_fname(ncFileInfo)
!#!
!#! type(netcdf_file_type), intent(in)    :: ncFileInfo
!#! character(len=80) :: nc_get_fname
!#!
!#! nc_get_fname = ncFileInfo%fname
!#!
!#! end function nc_get_fname
!#!
!#!
!#! !------------------------------------------------------------------
!#! !> get netcdf file id
!#!
!#!
!#! subroutine nc_set_ncid(ncFileInfo, my_ncid)
!#!
!#! type(netcdf_file_type), intent(inout) :: ncFileInfo
!#! integer,                intent(in)    :: my_ncid
!#!
!#! ncFileInfo%ncid = my_ncid
!#!
!#! end subroutine nc_set_ncid
!#!
!#!
!#! !------------------------------------------------------------------
!#! !> get netcdf file id
!#!
!#!
!#! function nc_get_ncid(ncFileInfo)
!#!
!#! type(netcdf_file_type), intent(in)    :: ncFileInfo
!#! integer :: nc_get_ncid
!#!
!#! nc_get_ncid = ncFileInfo%ncid
!#!
!#! end function nc_get_ncid


!------------------------------------------------------------------
!> return whether the init routines have been called?
!> or whether it's been marked as using a single file?
!>@todo fixme: document what this routine does.

function single_file_initialized(file_handle) result(is_initialized)

type(file_info_type), intent(inout) :: file_handle
logical :: is_initialized

is_initialized = file_handle%singlefile_initialized

end function single_file_initialized


!----------------------------------
!> Number of output ensemble members
!>@todo need to change the name to represent that this relates to the number of
!>      copies to write.
function noutput_state_variables(file_handle)
type(file_info_type), intent(in) :: file_handle !< file information handle
integer :: noutput_state_variables

noutput_state_variables = file_handle%stage_metadata%noutput_ens

end function noutput_state_variables


!------------------------------------------------------------------
!> Assert if a copy is in a valid range

function assert_valid_copy(name_handle, copy) result(valid_copy)
type(stage_metadata_type), intent(in) :: name_handle !< stage name handle
integer,                   intent(in) :: copy        !< ensemble copy number
logical :: valid_copy

integer :: size_fnames

size_fnames = size(name_handle%io_flag(:),1)

valid_copy = .false.
if(copy <= size_fnames .and. copy > 0) valid_copy = .true.

end function assert_valid_copy

!----------------------------------
end module io_filenames_mod
!> @}

