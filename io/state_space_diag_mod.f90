! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!>@todo FIXME need to rename to single_file_output_mod or something more
!>            representative
module state_space_diag_mod

!> \defgroup state_space_diag_mod state_space_diag_mod
!> Diagnostic file creation, writing, and closing.
!>
!> Usage:
!>   call init_singlefile_output(...)
!>   call write_singlefile(...)
!>   call finalize_singlefile_output(...)
!>
!>   these calls must all be collective calls. This allows tasks other than zero
!>   to write diagnostic information and also allows parallel IO within state_space_diag_mod.
!>
!> note 'all time steps' means all time steps on the output_interval
!> set in filter.
!> different options for diagnostic files:
!> 1. single file - all copies, all time steps in one file.
!>      model_mod::nc_write_model_atts is called which then returns
!>      a flag whether DART should define and write the state variables.
!>      this allows the model to still write whatever it wants to 
!>      the diagnostic file, and still have DART write the state if needed
!>      using the following routines:
!>         * write_model_attributes
!>         * write_model_variables
!>
!> 2. one copy per file, but all timesteps
!>      not sure if there is any desire for this.
!>      routines to write:
!>         * init_diag_one_copy_per_file
!>         * dart_nc_write_model_atts_one_copy_per_file
!>         * dart_nc_write_model_vars_one_copy_per_file
!>         * finalize_diag_output_one_copy_per_file
!>
!> 3. one copy per file, one timestep
!>      this is for large models
!>      here IO time is a concern, e.g. 0.1 degree POP where 60% of
!>      the run time of filter is in transposing and reading/writing
!>      restart files.
!>          * This is only for ONE timestep runs. - filter must have the output_interval = 1
!>
!>  a large amount of code in this module was moved from assim_model_mod and smoother_mod.
!>  some routines are only used by the program rms_diag.f90. It is believed that this program
!>  is not in use. There has been some discusion on whether to deprecate assim_model_type
!>  also.
!>
!> @{

use        types_mod,     only : r8, i8, digits12, metadatalength
use time_manager_mod,     only : time_type, get_time, get_calendar_type, &
                                 THIRTY_DAY_MONTHS, JULIAN, GREGORIAN, NOLEAP, &
                                 operator(<), operator(>), operator(+), &
                                 operator(-), operator(/), operator(*), &
                                 operator(==), operator(/=)
use ensemble_manager_mod, only : ensemble_type, map_task_to_pe, get_copy, &
                                 all_copies_to_all_vars, all_vars_to_all_copies, &
                                 get_allow_transpose, allocate_vars
use assim_model_mod,      only : assim_model_type, get_model_size
use model_mod,            only : nc_write_model_vars, nc_write_model_atts
use mpi_utilities_mod,    only : my_task_id, broadcast_flag, task_count
use utilities_mod,        only : error_handler, E_MSG, E_ERR, E_DBG, E_WARN, &
                                 file_to_text, find_textfile_dims, nc_check, &
                                 register_module
use io_filenames_mod,     only : file_info_type, stage_metadata_type, &
                                 noutput_state_variables, get_stage_metadata, &
                                 get_copy_name, netcdf_file_type, file_info_dump, &
                                 get_restart_filename, READ_COPY, WRITE_COPY
use state_structure_mod,  only : get_num_domains, get_io_num_dims, &
                                 get_dim_name, get_num_dims, get_dim_length, &
                                 get_io_num_unique_dims, get_io_unique_dim_name, &
                                 get_io_unique_dim_length, get_dim_lengths, &
                                 get_num_variables, get_variable_name, &
                                 get_variable_size, set_var_id, &
                                 get_index_start, get_index_end,  &
                                 create_diagnostic_structure, &
                                 end_diagnostic_structure
                                 
use model_mod,            only : read_model_time
use netcdf
use typesizes ! Part of netcdf?

implicit none
private

public :: init_singlefile_output, &          
          finalize_singlefile_output, &
          get_netcdf_file_type, &
          read_singlefile, &
          write_singlefile

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$url: https://svn-dares-dart.cgd.ucar.edu/DART/branches/rma_single_file/io/state_space_diag_mod.f90 $"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

logical :: module_initialized = .false.

! global storage for error/message string output
character(len=512)  :: msgstring

!-------------------------------------------------------------------------------

contains


!-------------------------------------------------------------------------------
! Routines that were in assim_model_mod.
!-------------------------------------------------------------------------------
!> Creates a diagnostic file(s).
!> Calls the model for any model specific attributes to be written
!> Leaves the diagnostic file open and passes out a handle: ncFileID (netcdf_file_type)

subroutine init_singlefile_output(ens_handle, file_handle)

type(ensemble_type),  intent(inout) :: ens_handle
type(file_info_type), intent(inout) :: file_handle

! Typical sequence:
! NF90_OPEN             ! create netCDF dataset: enter define mode
!    NF90_def_dim       ! define dimenstions: from name and length
!    NF90_def_var       ! define variables: from name, type, and dims
!    NF90_put_att       ! assign attribute values
! NF90_ENDDEF           ! end definitions: leave define mode
!    NF90_put_var       ! provide values for variable
! NF90_CLOSE            ! close: save updated netCDF dataset
!
! Time is a funny beast ... 
! Many packages decode the time:units attribute to convert the offset to a calendar
! date/time format. Using an offset simplifies many operations, but is not the
! way we like to see stuff plotted. The "approved" calendars are:
! gregorian or standard 
!      Mixed Gregorian/Julian calendar as defined by Udunits. This is the default. 
!  noleap   Modern calendar without leap years, i.e., all years are 365 days long. 
!  360_day  All years are 360 days divided into 30 day months. 
!  julian   Julian calendar. 
!  none     No calendar. 
!
! location is another one ...
!

! Variables for outputing input.nml variable
integer :: metadata_length, nlines, linelen, createmode
character(len=129),            allocatable :: textblock(:)
character(len=metadatalength), allocatable :: state_meta(:)

! Netcdf variables
type(netcdf_file_type) :: ncFileID
integer ::   MemberDimID
integer ::     TimeDimID,      TimeVarID
integer :: MetaDataDimID,  MetadataVarID
integer ::   nlinesDimID,   linelenDimID, nmlVarID

! For single file output, dart will always have control
! of the file formatting
logical :: local_model_mod_will_write_state_variables = .false.

! Local variables including counters and storing names
character(len=256) :: fname, copyname
integer :: icopy, ivar, ret, ens_size, num_output_ens

if (my_task_id() == 0) then
   if(.not. byteSizesOK()) then
       call error_handler(E_ERR,'init_singlefile_output', &
      'Compiler does not support required kinds of variables.',source,revision,revdate) 
   end if

   num_output_ens = file_handle%stage_metadata%noutput_ens
   if (num_output_ens <= 0) num_output_ens = 1

   allocate(state_meta(num_output_ens))
   ! Set up the metadata for the output file
   do ivar = 1, num_output_ens
      write(state_meta(ivar), '(a15, 1x, i6)') 'ensemble member', ivar
      write(file_handle%stage_metadata%copy_name(ivar),'(a,i4.4)') 'ens_mem_',ivar
   end do
   
   metadata_length = LEN(state_meta(1))
   
   ! NetCDF large file support
   createmode = NF90_64BIT_OFFSET
   
   ! Create the file
   fname          = get_restart_filename(file_handle%stage_metadata, 1, 1)
   ncFileID%fname = fname
   call nc_check(nf90_create(path=trim(fname), cmode=createmode, ncid=ncFileID%ncid), &
                 'init_singlefile_output', 'create '//trim(fname))

   ncFileID%ncid = ncFileID%ncid
   
   write(msgstring,*) trim(fname), ' is ncFileID ',ncFileID%ncid
   call error_handler(E_MSG,'init_singlefile_output',msgstring,source,revision,revdate)
   
   ! Define the dimensions
   call nc_check(nf90_def_dim(ncid=ncFileID%ncid, &
                 name="metadatalength", len=metadata_length, dimid=MetaDataDimID), &
                 'init_singlefile_output', 'def_dim metadatalength '//trim(fname))
   
   call nc_check(nf90_def_dim(ncid=ncFileID%ncid, &
                 name="member", len=num_output_ens, dimid=MemberDimID), &
                 'init_singlefile_output', 'def_dim member '//trim(fname))
   
   call nc_check(nf90_def_dim(ncid=ncFileID%ncid, name="time", &
                 len=nf90_unlimited, dimid=TimeDimID), &
                 'init_singlefile_output', 'def_dim time '//trim(fname))

   !----------------------------------------------------------------------------
   ! Find dimensions of namelist file ... will save it as a variable.
   !----------------------------------------------------------------------------
   
   ! All DART programs require input.nml, so it is unlikely this can fail, but
   ! still check and in this case, error out if not found.
   call find_textfile_dims("input.nml", nlines, linelen)
   if (nlines <= 0 .or. linelen <= 0) then
      call error_handler(E_MSG,'init_singlefile_output', &
                         'cannot open/read input.nml to save in diagnostic file', &
                         source,revision,revdate)
   endif
   
   allocate(textblock(nlines))
   textblock = ''
   
   call nc_check(nf90_def_dim(ncid=ncFileID%ncid, &
                 name="NMLlinelen", len=LEN(textblock(1)), dimid=linelenDimID), &
                 'init_singlefile_output', 'def_dim NMLlinelen '//trim(fname))
   
   call nc_check(nf90_def_dim(ncid=ncFileID%ncid, &
                 name="NMLnlines", len=nlines, dimid=nlinesDimID), &
                 'init_singlefile_output', 'def_dim NMLnlines '//trim(fname))
   
   !----------------------------------------------------------------------------
   ! Write Global Attributes 
   !----------------------------------------------------------------------------
   
   call nc_check(nf90_put_att(ncFileID%ncid, NF90_GLOBAL, "title", fname), &
                 'init_singlefile_output', 'put_att title '//trim(fname))
   call nc_check(nf90_put_att(ncFileID%ncid, NF90_GLOBAL, "assim_model_source", source ), &
                 'init_singlefile_output', 'put_att assim_model_source '//trim(fname))
   call nc_check(nf90_put_att(ncFileID%ncid, NF90_GLOBAL, "assim_model_revision", revision ), &
                 'init_singlefile_output', 'put_att assim_model_revision '//trim(fname))
   call nc_check(nf90_put_att(ncFileID%ncid, NF90_GLOBAL, "assim_model_revdate", revdate ), &
                 'init_singlefile_output', 'put_att assim_model_revdate '//trim(fname))
   
   !    Metadata for each Copy
   call nc_check(nf90_def_var(ncid=ncFileID%ncid, name="MemberMetadata", xtype=nf90_char,    &
                 dimids=(/ MetaDataDimID, MemberDimID /),  varid=metadataVarID), &
                 'init_singlefile_output', 'def_var MemberMetadata')
   call nc_check(nf90_put_att(ncFileID%ncid, metadataVarID, "long_name",       &
                 "Metadata for each copy/member"), 'init_singlefile_output', 'put_att long_name')
   
   !    Input namelist 
   call nc_check(nf90_def_var(ncid=ncFileID%ncid,name="inputnml", xtype=nf90_char,    &
                 dimids = (/ linelenDimID, nlinesDimID /),  varid=nmlVarID), &
                 'init_singlefile_output', 'def_var inputnml')
   call nc_check(nf90_put_att(ncFileID%ncid, nmlVarID, "long_name",       &
                 "input.nml contents"), 'init_singlefile_output', 'put_att input.nml')
   
   !    Time -- the unlimited dimension
   call nc_check(nf90_def_var(ncFileID%ncid, name="time", xtype=nf90_double, dimids=TimeDimID, &
                 varid =TimeVarID), 'init_singlefile_output', 'def_var time' )
   ret = nc_write_calendar_atts(ncFileID, TimeVarID)     ! comes from time_manager_mod
   if ( ret /= 0 ) then
      write(msgstring, *)'nc_write_calendar_atts  bombed with error ', ret
      call error_handler(E_MSG,'init_singlefile_output',msgstring,source,revision,revdate)
   endif
   
   ! Create the time "mirror" with a static length. There is another routine
   ! to increase it if need be. For now, just pick something.
   ncFileID%Ntimes    = 0
   ncFileID%NtimesMAX = 1000
   allocate(ncFileID%rtimes(ncFileID%NtimesMAX), ncFileID%times(ncFileID%NtimesMAX) )
   
   !----------------------------------------------------------------------------
   ! Leave define mode so we can fill
   !----------------------------------------------------------------------------
   
   call nc_check(nf90_enddef(ncFileID%ncid), 'init_singlefile_output', 'enddef '//trim(fname))

   call nc_check(nf90_sync(ncFileID%ncid), 'init_singlefile_output', 'sync '//trim(fname))               
   !----------------------------------------------------------------------------
   ! Define the model-specific components
   !----------------------------------------------------------------------------
   ret =  nc_write_model_atts( ncFileID%ncid, local_model_mod_will_write_state_variables)
   if ( ret /= 0 ) then
      write(msgstring, *)'nc_write_model_atts  bombed with error ', ret
      call error_handler(E_MSG,'init_diag_output',msgstring,source,revision,revdate)
   endif

   if ( .not. local_model_mod_will_write_state_variables ) then
      call write_model_attributes(ncFileID, MemberDimID, TimeDimID)
   endif

   !----------------------------------------------------------------------------
   ! Create variables and attributes.
   ! The locations are part of the model (some models have multiple grids).
   ! They are written by model_mod:nc_write_model_atts
   !----------------------------------------------------------------------------

   ens_size = ens_handle%num_copies - ens_handle%num_extras
   do icopy = ens_size+1, ens_handle%num_copies 
      copyname = trim(get_copy_name(file_handle,icopy))
      if ( file_handle%stage_metadata%io_flag(icopy) == WRITE_COPY ) then
         call write_extra_attributes(ncFileID, TimeDimID, copyname)
      endif
   enddo 
   

   
   !----------------------------------------------------------------------------
   ! Fill the coordinate variables.
   ! Write the input namelist as a netCDF variable.
   ! The time variable is filled as time progresses.
   !----------------------------------------------------------------------------
   
   call nc_check(nf90_put_var(ncFileID%ncid, metadataVarID, state_meta ), &
                 'init_singlefile_output', 'put_var MetaDataVarID')
    
   call file_to_text("input.nml", textblock)
   
   call nc_check(nf90_put_var(ncFileID%ncid, nmlVarID, textblock ), &
                 'init_singlefile_output', 'put_var nmlVarID')
   
   deallocate(textblock)
   
   !----------------------------------------------------------------------------
   ! sync to disk, but leave open
   !----------------------------------------------------------------------------
   
   call nc_check(nf90_sync(ncFileID%ncid), 'init_singlefile_output', 'sync '//trim(fname))               
   !----------------------------------------------------------------------------
   ! sync again, but still leave open
   !----------------------------------------------------------------------------
   
   call nc_check(nf90_sync(ncFileID%ncid), 'init_singlefile_output', 'sync '//trim(fname))               
   call set_netcdf_file_type(file_handle, ncFileID)
 
endif

! Broadcast the value of model_mod_will_write_state_variables to every task
! This keeps track of whether the model_mod or dart code will write state_variables.
call broadcast_flag(local_model_mod_will_write_state_variables, 0)
ncFileID%model_mod_will_write_state_variables = local_model_mod_will_write_state_variables

file_handle%singlefile_initialized = .true.

end subroutine init_singlefile_output


!-------------------------------------------------------------------------------
!>

function finalize_singlefile_output(ncFileID) result(ierr)

type(netcdf_file_type), intent(inout) :: ncFileID
integer             :: ierr

ierr = 0

if (my_task_id()==0) then
   ierr = NF90_close(ncFileID%ncid)
   if(associated(ncFileID%rtimes)) deallocate(ncFileID%rtimes, ncFileID%times )
endif

call end_diagnostic_structure()

ncFileID%fname     = "notinuse"
ncFileID%ncid      = -1
ncFileID%Ntimes    = -1
ncFileID%NtimesMax = -1

end function finalize_singlefile_output


!-------------------------------------------------------
!> read a single netcdf file containing all of the members
!> and possibly inflation information
subroutine read_singlefile(state_ens_handle, file_info, use_time_from_file, time)

type(ensemble_type),  intent(inout) :: state_ens_handle
type(file_info_type), intent(in)    :: file_info
logical,              intent(in)    :: use_time_from_file
type(time_type),      intent(inout) :: time

type(stage_metadata_type) :: name_handle

! NetCDF IO variables
integer :: my_ncid, varid, MemDimID, TimeDimID, ret, icopy, ivar, jdim, istart, iend
integer :: ens_size, mem_size, time_size, var_size, domain
integer :: num_output_ens, ncount, ndims
real(r8), allocatable :: var_block(:)
character(len=NF90_MAX_NAME) :: fname, dimname, varname, copyname, extraname
integer, dimension(NF90_MAX_VAR_DIMS) :: dim_lengths
integer, dimension(NF90_MAX_VAR_DIMS) :: start_point

! check whether file_info handle is initialized
if (task_count() > 1) then
   call error_handler(E_ERR,'read_singlefile', &
   'current version of the system does not support MPI with single input/output files', &
    source, revision, revdate, text2='compile without MPI, or use multi-file i/o')
endif

call allocate_vars(state_ens_handle)

!>@todo FIXME aren't we reading in the initial data here?
!> so what's in copies that needs to be moved to vars?
!call all_copies_to_all_vars(state_ens_handle)

! do this once
name_handle = get_stage_metadata(file_info)

ens_size = state_ens_handle%num_copies - state_ens_handle%num_extras

domain = 1 !>@todo : only a single domain for single file read supported. need
           !>        to consider case for multiple domains.
fname = get_restart_filename(file_info%stage_metadata, 1, domain)

! debug only: 
!call file_info_dump(file_info,'read_singlefile : ')

! read time from input file if time not set in namelist
!>@todo Check time consistency across files? This is assuming they are consistent.
!debug: print*, 'single file : use_time_from_file = ', use_time_from_file, trim(fname)
if(use_time_from_file) then
   time = read_model_time(fname)
endif

state_ens_handle%time = time

ret = nf90_open(fname, NF90_NOWRITE, my_ncid)
call nc_check(ret, 'read_singlefile: nf90_open', fname)

ret = nf90_inq_dimid(my_ncid, "time", TimeDimID)
call nc_check(ret, 'read_singlefile', 'inq_varid time : '//trim(fname))

ret = nf90_inquire_dimension(my_ncid, TimeDimID, len=time_size) 
call nc_check(ret, 'read_singlefile', 'inquire_dimension time '//trim(fname))

istart = 1
do ivar = 1, get_num_variables(domain)
   var_size = get_variable_size(domain, ivar)
   allocate( var_block(var_size) )
   
   varname = get_variable_name(domain, ivar)
   ret      = nf90_inq_varid(my_ncid, varname, varid)
   call nc_check(ret, 'read_singlefile', 'inq_varid '//trim(varname)//' : '//trim(fname))

   ! if member dimension is defined then we are reading from a single file
   !>@todo there is probably a cleaner way of doing this
   ret = nf90_inq_dimid(my_ncid, "member", MemDimID)
   if( ret == 0 ) then
      ret = nf90_inquire_dimension(my_ncid, MemDimID, len=mem_size) 
      call nc_check(ret, 'read_singlefile', 'inq_varid member : '//trim(fname))
      num_output_ens = ens_size
      if (mem_size < ens_size) num_output_ens = mem_size
   endif

   ncount = 0
   ndims = get_num_dims(domain, ivar)
   do jdim = 1, ndims
      dimname = get_dim_name(domain, ivar, jdim)
      if ( trim(dimname) == 'time' .or. trim(dimname) == 'member') cycle

      ncount = ncount+1
      dim_lengths(ncount) = get_dim_length(domain, ivar, jdim)
   enddo

   dim_lengths(ncount + 1) = 1 ! member
   dim_lengths(ncount + 2) = 1 ! time

   iend   = istart + var_size - 1
   do icopy = 1, num_output_ens
      start_point(1:ncount) = 1
      start_point(ncount + 1) = icopy  ! member
      start_point(ncount + 2) = time_size ! time

      ret = nf90_get_var(my_ncid, varid, var_block, count=dim_lengths(1:ncount+2), &
                         start=start_point(1:ncount+2))
      call nc_check(ret, 'read_singlefile', 'get_var '//trim(varname)//' : '//trim(fname))
      state_ens_handle%copies(icopy, istart:iend) = var_block
 

   enddo

   istart = iend + 1

   deallocate( var_block )
enddo

! Read extra opies into ensemble handle, this could include
!   {variable}_{mean,sd}
!   {variable}_priorinf_{mean,sd}
!   {variable}_postinf_{mean,sd}

do ivar = 1, get_num_variables(domain) ! assuming one domain for single files
   var_size = get_variable_size(domain, ivar)
   varname  = get_variable_name(1,ivar) 
   allocate( var_block(var_size) )

   do icopy = ens_size+1, state_ens_handle%num_copies 
      copyname = get_copy_name(file_info,icopy)

      if ( file_info%stage_metadata%io_flag(icopy) == READ_COPY) then

         ncount = 0
         ndims = get_num_dims(domain, ivar)
         do jdim = 1, ndims
            dimname = get_dim_name(domain, ivar, jdim)
            if ( trim(dimname) == 'time' .or. trim(dimname) == 'member') cycle

            ncount = ncount+1
            dim_lengths(ncount) = get_dim_length(domain, ivar, jdim)
         enddo

         dim_lengths(ncount + 1) = 1 ! time

         start_point(1:ncount)    = 1
         start_point(ncount + 1) = time_size ! time

         write(extraname,'(a,"_",a)') trim(varname), trim(copyname)
         ret = nf90_inq_varid(my_ncid, extraname, varid)
         call nc_check(ret, 'read_singlefile', 'inq_varid '//trim(extraname))

         ret = nf90_get_var(my_ncid, varid, var_block, count=dim_lengths(1:ncount+1), &
                            start=start_point(1:ncount+1))
         call nc_check(ret, 'read_singlefile', 'get_var '//trim(varname)//' : '//trim(fname))
         state_ens_handle%copies(icopy, istart:iend) = var_block

         istart = 1
         iend   = var_size
         state_ens_handle%copies(icopy, istart:iend) = var_block
      endif
   enddo

   deallocate( var_block )
enddo 

ret = nf90_close(my_ncid)
call nc_check(ret, 'read_singlefile: nf90_close', fname)

! and now distribute the var data to the copies
! call all_vars_to_all_copies(state_ens_handle)

end subroutine read_singlefile


!-----------------------------------------------------------
!> several options for diagnostic files:
!>  * model_mod writes the diagnostic file (all copies, all output time steps in one file)
!>  * dart writes the diagnostic file (all copies, all output time steps in one file)
!>  * dart writes the diagnostic files (one copy, all timesteps per file)
!>  * single timestep runs:
!>       - If num_output_ens is 0, the prior ENS_MEAN_COPY, ENS_SD_COPY, INF_COPY, INF_SD_COPY
!>         are saved until the end.  This is for large models to reduce IO (transpose time).

subroutine write_singlefile(ens_handle, file_info)

type(ensemble_type),         intent(inout) :: ens_handle
type(file_info_type),        intent(inout) :: file_info


! Local variables
integer                :: copyindex, my_ncid
integer                :: num_output_ens, ens_size
integer(i8)            :: model_size
character(len=128)     :: copyname
type(time_type)        :: curr_ens_time
type(netcdf_file_type) :: ncFileID
real(r8), allocatable  :: temp_ens(:)

! assumes that mean and spread have already been computed
! make sure vars is up-to-date
call allocate_vars(ens_handle)
call all_copies_to_all_vars(ens_handle)

model_size = get_model_size()

! task 0 needs some space
if (my_task_id() == 0) then
   allocate(temp_ens(model_size))
else
   allocate(temp_ens(1))
endif

curr_ens_time = ens_handle%current_time
call get_netcdf_file_type(file_info, ncFileID) 

my_ncid        = ncFileID%ncid
num_output_ens = noutput_state_variables(file_info)

! Output ensemble members
do copyindex = 1, num_output_ens
   call get_copy(map_task_to_pe(ens_handle, 0), ens_handle, copyindex, temp_ens)
   call write_model_variables(ncFileID,  temp_ens, copyindex, curr_ens_time)
enddo

! Output Extras
ens_size = ens_handle%num_copies - ens_handle%num_extras

do copyindex = ens_size+1, ens_handle%num_copies
   if ( file_info%stage_metadata%io_flag(copyindex) == WRITE_COPY ) then
      copyname = get_copy_name(file_info, copyindex)
      call get_copy(map_task_to_pe(ens_handle, 0), ens_handle, copyindex, temp_ens)
      if(my_task_id() == 0) then
         call write_extra_variable(ncFileID, temp_ens, copyname, curr_ens_time)
      endif
   endif
enddo

call set_netcdf_file_type(file_info, ncFileID)

deallocate(temp_ens)

!>@todo FIXME we haven't changed anything in the vars array
!> so there's no need to transpose back.  save time by skipping this.
!call all_vars_to_all_copies(ens_handle)

end subroutine write_singlefile


!-----------------------------------------------------------
!> write out a single variable appending time stamp curr_ens_time

subroutine write_extra_variable(ncFileID, model_state, cname, curr_ens_time)

type(netcdf_file_type), intent(inout) :: ncFileID
real(r8),               intent(in)    :: model_state(:)
character(len=*),       intent(in)    :: cname
type(time_type),        intent(in)    :: curr_ens_time

integer, dimension(NF90_MAX_VAR_DIMS) :: dim_lengths
integer, dimension(NF90_MAX_VAR_DIMS) :: start_point
integer :: istart, iend
integer :: ivar, jdim
integer :: ndims, ncount
integer :: ret ! netcdf return code
integer :: var_id ! netcdf variable id
integer :: domain
integer :: my_ncid
integer :: timeindex
integer :: is1, id1
character(len=NF90_MAX_NAME) :: dimname, varname, extraname

my_ncid = ncFileID%ncid
timeindex = nc_get_tindex(ncFileID, curr_ens_time)

if ( timeindex < 0 ) then
   call get_time(curr_ens_time,is1,id1)
   write(msgstring,*)'model time (d,s)',id1,is1,' not in ',ncFileID%fname
   write(msgstring,'(''model time (d,s) ('',i8,i5,'') is index '',i6, '' in ncFileID '',i10)') &
          id1,is1,timeindex,ncFileID%ncid
   call error_handler(E_ERR,'write_singlefile', msgstring, source, revision, revdate)
endif

call get_time(curr_ens_time,is1,id1)
write(msgstring,'(''model time (d,s) ('',i8,i5,'') is index '',i6, '' in ncFileID '',i10)') &
       id1,is1,timeindex,ncFileID%ncid
call error_handler(E_DBG,'write_singlefile', msgstring, source, revision, revdate)

!#! print*, 'write_extra_variable ', trim(varname), ' timeind [',timeindex,'], size state = ', size(model_state)

domain    = 1 !>@todo ONLY ONE DOMAIN FOR SINGLE FILE OUTPUT

do ivar = 1, get_num_variables(domain)

   istart   = get_index_start(   domain, ivar )
   iend     = get_index_end(     domain, ivar )
   ndims    = get_num_dims(      domain, ivar )
   varname  = get_variable_name( domain, ivar )

   !#! print*, 'number of dimensions ', ndims
   ncount = 0
   do jdim = 1, ndims
      dimname = get_dim_name(domain, ivar, jdim)
      if ( trim(dimname) == 'time' .or. trim(dimname) == 'member') cycle

      ncount = ncount+1
      dim_lengths(ncount) = get_dim_length(domain, ivar, jdim)
   enddo

   !#! dim_lengths(ncount + 1) = 1 ! member
   dim_lengths(ncount+1) = 1 ! time

   start_point(1:ncount) = 1
   !#! start_point(ncount + 1) = memindex  ! member
   start_point(ncount + 1) = timeindex ! time

   write(extraname,'(a,"_",a)') trim(varname), trim(cname)
   !#! write(*        ,'(a,"_",a)') trim(varname), trim(cname)

   ret = nf90_inq_varid(my_ncid, extraname, var_id)
   call nc_check(ret, 'write_extra_variable', 'inq_varid '//trim(extraname))

   ret = nf90_put_var(my_ncid, var_id, model_state(istart:iend), &
                count=dim_lengths(1:ncount+2), start=start_point(1:ncount+2))
   call nc_check(ret, 'write_extra_variable', 'put_var '//trim(extraname))

enddo

end subroutine write_extra_variable


!-------------------------------------------------------------------------------
! DART writing diagnostic files
!-------------------------------------------------------------------------------
!> Write state to the last time slice of a file
!> This routine is called from write_singlefile if the model_mod
!> nc_write_model_vars has returned 
!> model_mod_will_write_state_varaibles = .false.

subroutine write_model_variables(ncFileID, model_state, memindex, curr_ens_time)

type(netcdf_file_type), intent(inout) :: ncFileID
real(r8),               intent(in) :: model_state(:)
integer,                intent(in) :: memindex
type(time_type),        intent(in) :: curr_ens_time

integer, dimension(NF90_MAX_VAR_DIMS) :: dim_lengths
integer, dimension(NF90_MAX_VAR_DIMS) :: start_point
integer :: istart, iend
integer :: ivar, jdim
integer :: ndims, ncount
integer :: ret ! netcdf return code
integer :: var_id ! netcdf variable id
integer :: domain
integer :: my_ncid
integer :: timeindex
integer :: is1, id1
character(len=NF90_MAX_NAME) :: dimname, varname

! may not be needed
ncFileID%diag_id = create_diagnostic_structure()
my_ncid          = ncFileID%ncid

timeindex = nc_get_tindex(ncFileID, curr_ens_time)

if ( timeindex < 0 ) then
   call get_time(curr_ens_time,is1,id1)
   write(msgstring,*)'model time (d,s)',id1,is1,' not in ',ncFileID%fname
   write(msgstring,'(''model time (d,s) ('',i8,i5,'') is index '',i6, '' in ncFileID '',i10)') &
          id1,is1,timeindex,ncFileID%ncid
   call error_handler(E_ERR,'write_model_variables', msgstring, source, revision, revdate)
endif

call get_time(curr_ens_time,is1,id1)
write(msgstring,'(''model time (d,s) ('',i8,i5,'') is index '',i6, '' in ncFileID '',i10)') &
       id1,is1,timeindex,ncFileID%ncid
call error_handler(E_DBG,'write_model_variables', msgstring, source, revision, revdate)


domain = 1 !>@todo ONLY ONE DOMAIN FOR SINGLE FILE OUTPUT

do ivar = 1, get_num_variables(domain)

   istart  = get_index_start(  domain, ivar)
   iend    = get_index_end(    domain, ivar)
   ndims   = get_num_dims(     domain, ivar)
   varname = get_variable_name(domain, ivar)

   ncount = 0
   do jdim = 1, ndims
      dimname = get_dim_name(domain, ivar, jdim)
      if ( trim(dimname) == 'time' .or. trim(dimname) == 'member') cycle

      ncount = ncount+1
      dim_lengths(ncount) = get_dim_length(domain, ivar, jdim)

   enddo

   dim_lengths(ncount + 1) = 1 ! member
   dim_lengths(ncount + 2) = 1 ! time

   start_point(1:ncount) = 1
   start_point(ncount + 1) = memindex  ! member
   start_point(ncount + 2) = timeindex ! time

   ret = nf90_inq_varid(my_ncid, varname, var_id)
   call nc_check(ret, 'write_model_variables', 'inq_varid '//trim(varname))

   ret = nf90_put_var(my_ncid, var_id, model_state(istart:iend), &
                count=dim_lengths(1:ncount+2), start=start_point(1:ncount+2))
   call nc_check(ret, 'write_model_variables', 'put_var '//trim(varname))

enddo

end subroutine write_model_variables


!-----------------------------------------------------------
!>

subroutine write_extra_attributes(ncFileID, time_dimId, cname)

type(netcdf_file_type), intent(inout) :: ncFileID
integer,          intent(in) :: time_dimId
character(len=*), intent(in) :: cname

integer :: ivar, jdim, ndims, countdims, domain
integer :: ret, my_ncid, new_varid, my_dimid, my_xtype
character(len=NF90_MAX_VAR_DIMS) :: varname, dimname, extraname
integer ::  model_dimids(NF90_MAX_VAR_DIMS)

!--------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file 
!--------------------------------------------------------------------

call nc_check(nf90_sync(ncFileID%ncid),  'write_extra_attributes', 'nf90_sync') ! Ensure netCDF file is current
call nc_check(nf90_Redef(ncFileID%ncid), 'write_extra_attributes', 'nf90_Redef')

my_ncid  = ncFileID%ncid

if(my_task_id()==0) then

   if (r8 == digits12) then
      my_xtype = nf90_double
   else
      my_xtype = nf90_real
   endif

   domain = 1 !>@todo ONLY ONE DOMAIN FOR SINGLE FILES
   ! Define dimensions for state
   do ivar = 1, get_num_variables(domain)
      varname = get_variable_name(domain, ivar)
      ndims   = get_num_dims(domain, ivar)
   
      countdims = 0
      do jdim = 1, ndims
         dimname = get_dim_name(domain, ivar, jdim)
         if ( trim(dimname) == 'time' .or. trim(dimname) == 'member') cycle
   
         countdims = countdims+1
         ret = nf90_inq_dimid(my_ncid, dimname, dimid=my_dimid)
         call nc_check(ret, 'write_extra_attributes', 'inq_dimid '//trim(dimname))
         model_dimids(countdims) = my_dimid
         !#! print*, ivar, jdim, ndims, trim(varname), ' ', trim(dimname), my_dimid
      enddo
      model_dimids(countdims+1) = time_dimId
   
      write(extraname,'(a,"_",a)') trim(varname), trim(cname)
      ret = nf90_def_var(my_ncid, name   = extraname,& 
                               xtype  = my_xtype, &
                               dimids = model_dimids(1:countdims+1), &
                               varid  = new_varid)
      call nc_check(ret, 'write_extra_attributes', 'defining variable '//trim(extraname))
   enddo
endif

! Leave define mode so we can fill
call nc_check(nf90_enddef(ncfileID%ncid), 'write_extra_attributes', 'nf90_enddef')

end subroutine write_extra_attributes


!-------------------------------------------------------------------------------
!> Called if the model_mod is NOT going to write the state variables to
!> the diagnostic file.
!> This routine defines the state variables in the ncFileID
!> Time and copy dimensions are already defined by init_singlefile_output
!> State variables are defined:
!>   variable(dim1, dim2, ...,  copy, time)
!> If there are multiple domains the variables and dimensions are
!> given the suffix _d0*, where * is the domain number.

subroutine write_model_attributes(ncFileID, member_dimID, time_dimId)

type(netcdf_file_type), intent(inout) :: ncFileID
integer,                intent(in)    :: member_dimID
integer,                intent(in)    :: time_dimId

integer :: domain, my_ncid ! local variable
integer :: ret ! netcdf return code
integer :: dummy
integer :: dimids(NF90_MAX_VAR_DIMS)
integer :: ndims
integer :: my_xtype ! precision for netcdf variable
integer :: ivar, jdim ! loop variables
integer :: new_varid
character(len=metadatalength) :: dimname, varname

call nc_check(nf90_sync(ncFileID%ncid), 'write_model_attributes', 'nf90_sync') ! Ensure netCDF file is current
call nc_check(nf90_Redef(ncFileID%ncid), 'write_model_attributes', 'nf90_Redef')

ncFileID%diag_id = create_diagnostic_structure()

my_ncid  = ncFileID%ncid

if(my_task_id()==0) then

   if (r8 == digits12) then
      my_xtype = nf90_double
   else
      my_xtype = nf90_real
   endif

   domain = 1 !>@todo ONLY ONE DOMAIN FOR SINGLE FILES
   ! Define dimensions for state
   do ivar = 1, get_num_variables(domain)

      ndims = get_num_dims(domain, ivar)

      do jdim = 1, ndims
         dimname = get_dim_name(domain, ivar, jdim)
         ret = nf90_def_dim(my_ncid, dimname, get_dim_length(domain, ivar, jdim), dummy)
         !>@todo if we already have a unique names we can take this test out
         if(ret /= NF90_NOERR .and. ret /= NF90_ENAMEINUSE) then
            call nc_check(ret, 'write_model_attributes', 'defining dimensions')
         endif
      enddo

      ! Define variables
      ! query the dimension ids
      do jdim = 1, ndims
         dimname = get_dim_name(domain, ivar, jdim)
         ret = nf90_inq_dimid(my_ncid, dimname, dimids(jdim))
         call nc_check(ret, 'write_model_attributes', 'querying dimensions')
      enddo

      dimids(ndims + 1) = member_dimID
      dimids(ndims + 2) = time_dimId

      varname = get_variable_name(domain, ivar)
      ret = nf90_def_var(my_ncid, name   = varname, &
                               xtype  = my_xtype, &
                               dimids = dimids(1:ndims+2), &
                               varid  = new_varid)
      call nc_check(ret, 'write_model_attributes', 'defining variable '//trim(varname))
      call set_var_id(domain, ivar, new_varid)

   enddo

! Leave define mode
ret = nf90_enddef(my_ncid)
call nc_check(ret, 'write_model_attributes', 'enddef')

endif

end subroutine write_model_attributes


!-------------------------------------------------------------------------------
!> routine write callendar attributes

function nc_write_calendar_atts(ncFileID, TimeVarID) result(ierr)

type(netcdf_file_type), intent(in) :: ncFileID
integer,                intent(in) :: TimeVarID
integer                            :: ierr

integer :: my_ncid

ierr = 0

my_ncid = ncFileID%ncid

call nc_check(nf90_put_att(my_ncid, TimeVarID, "long_name", "time"), &
              'nc_write_calendar_atts', 'put_att long_name '//trim(ncFileID%fname))
call nc_check(nf90_put_att(my_ncid, TimeVarID, "axis", "T"), &
              'nc_write_calendar_atts', 'put_att axis '//trim(ncFileID%fname))
call nc_check(nf90_put_att(my_ncid, TimeVarID, "cartesian_axis", "T"), &
              'nc_write_calendar_atts', 'put_att cartesian_axis '//trim(ncFileID%fname))

select case( get_calendar_type() )
case(THIRTY_DAY_MONTHS)
!  call get_date_thirty(time, year, month, day, hour, minute, second)
case(GREGORIAN)
   call nc_check(nf90_put_att(my_ncid, TimeVarID, "calendar", "gregorian" ), &
              'nc_write_calendar_atts', 'put_att calendar '//trim(ncFileID%fname))
   call nc_check(nf90_put_att(my_ncid, TimeVarID, "units", "days since 1601-01-01 00:00:00"), &
              'nc_write_calendar_atts', 'put_att units '//trim(ncFileID%fname))
case(JULIAN)
   call nc_check(nf90_put_att(my_ncid, TimeVarID, "calendar", "julian" ), &
              'nc_write_calendar_atts', 'put_att calendar '//trim(ncFileID%fname))
case(NOLEAP)
   call nc_check(nf90_put_att(my_ncid, TimeVarID, "calendar", "no_leap" ), &
              'nc_write_calendar_atts', 'put_att calendar '//trim(ncFileID%fname))
case default
   call nc_check(nf90_put_att(my_ncid, TimeVarID, "calendar", "no calendar" ), &
              'nc_write_calendar_atts', 'put_att calendar '//trim(ncFileID%fname))
   call nc_check(nf90_put_att(my_ncid, TimeVarID, "month_lengths", &
              (/31,28,31,30,31,30,31,31,30,31,30,31/)), &
              'nc_write_calendar_atts', 'put_att month_lengths '//trim(ncFileID%fname))
   call nc_check(nf90_put_att(my_ncid, TimeVarID, "units", &
              'days since 0000-01-01 00:00:00'), &
              'nc_write_calendar_atts', 'put_att units '//trim(ncFileID%fname))
end select

end function nc_write_calendar_atts


!-------------------------------------------------------------------------------
!> We need to compare the time of the current assim_model to the
!> netcdf time coordinate variable (the unlimited dimension).
!> If they are the same, no problem ...
!> If it is earlier, we need to find the right index and insert ...
!> If it is the "future", we need to add another one ...
!> If it is in the past but does not match any we have, we're in trouble.
!> The new length of the "time" variable is returned.
!>
!> A "times" array has been added to mirror the times that are stored
!> in the netcdf time coordinate variable. While somewhat unpleasant, it
!> is SUBSTANTIALLY faster than reading the netcdf time variable at every
!> turn -- which caused a geometric or exponential increase in overall 
!> netcdf I/O. (i.e. this was really bad)
!>
!> The time mirror is maintained as a time_type, so the comparison with
!> the state time uses the operators for the time_type. The netCDF file,
!> however, has time units of a different convention. The times are
!> converted only when appending to the time coordinate variable.    


function nc_get_tindex(ncFileID, statetime) result(timeindex)

type(netcdf_file_type), intent(inout) :: ncFileID
type(time_type), intent(in) :: statetime
integer                     :: timeindex

integer  :: nDimensions, nVariables, nAttributes, unlimitedDimID, TimeVarID
integer  :: xtype, ndims, nAtts, nTlen
integer  :: secs, days, my_ncid, i

character(len=NF90_MAX_NAME)          :: varname
integer, dimension(NF90_MAX_VAR_DIMS) :: dimids

timeindex = -1  ! assume bad things are going to happen

my_ncid = ncFileID%ncid

! Make sure we're looking at the most current version of the netCDF file.
! Get the length of the (unlimited) Time Dimension 
! If there is no length -- simply append a time to the dimension and return ...
! Else   get the existing times ["days since ..."] and convert to time_type 
!        if the statetime < earliest netcdf time ... we're in trouble
!        if the statetime does not match any netcdf time ... we're in trouble
!        if the statetime > last netcdf time ... append a time ... 

call nc_check(NF90_Sync(my_ncid), 'nc_get_tindex', 'sync '//trim(ncFileID%fname))    
call nc_check(NF90_Inquire(my_ncid, nDimensions, nVariables, nAttributes, unlimitedDimID), &
              'nc_get_tindex', 'inquire '//trim(ncFileID%fname))
call nc_check(NF90_Inq_Varid(my_ncid, "time", TimeVarID), &
              'nc_get_tindex', 'inq_varid time '//trim(ncFileID%fname))
call nc_check(NF90_Inquire_Variable(my_ncid, TimeVarID, varname, xtype, ndims, dimids, nAtts), &
              'nc_get_tindex', 'inquire_variable time '//trim(ncFileID%fname))
call nc_check(NF90_Inquire_Dimension(my_ncid, unlimitedDimID, varname, nTlen), &
              'nc_get_tindex', 'inquire_dimension unlimited '//trim(ncFileID%fname))
! Sanity check all cases first.

if ( ndims /= 1 ) then
   write(msgstring,*)'"time" expected to be rank-1' 
   call error_handler(E_WARN,'nc_get_tindex',msgstring,source,revision,revdate)
   timeindex = timeindex -   1
endif
if ( dimids(1) /= unlimitedDimID ) then
   write(msgstring,*)'"time" must be the unlimited dimension'
   call error_handler(E_WARN,'nc_get_tindex',msgstring,source,revision,revdate)
   timeindex = timeindex -  10
endif
if ( timeindex < -1 ) then
   write(msgstring,*)'trouble deep ... can go no farther. Stopping.'
   call error_handler(E_ERR,'nc_get_tindex',msgstring,source,revision,revdate)
endif

! convert statetime to time base of "days since ..."
call get_time(statetime, secs, days)

if (ncFileID%Ntimes < 1) then          ! First attempt at writing a state ...

   write(msgstring,*)'current unlimited  dimension length',nTlen, &
                     'for ncFileID ',trim(ncFileID%fname)
   call error_handler(E_DBG,'nc_get_tindex',msgstring,source,revision,revdate)
   write(msgstring,*)'current time array dimension length',ncFileID%Ntimes
   call error_handler(E_DBG,'nc_get_tindex',msgstring,source,revision,revdate)

   nTlen = nc_append_time(ncFileID, statetime)

   write(msgstring,*)'Initial time array dimension length',ncFileID%Ntimes
   call error_handler(E_DBG,'nc_get_tindex',msgstring,source,revision,revdate)

endif


TimeLoop : do i = 1,ncFileId%Ntimes

   if ( statetime == ncFileID%times(i) ) then
      timeindex = i
      exit TimeLoop
   endif

enddo TimeLoop



if ( timeindex <= 0 ) then   ! There was no match. Either the model
                             ! time precedes the earliest file time - or - 
                             ! model time is somewhere in the middle  - or - 
                             ! model time needs to be appended.

   if (statetime < ncFileID%times(1) ) then

      call error_handler(E_DBG,'nc_get_tindex', &
              'Model time precedes earliest netCDF time.', source,revision,revdate)

      write(msgstring,*)'          model time (days, seconds) ',days,secs
      call error_handler(E_DBG,'nc_get_tindex',msgstring,source,revision,revdate)

      call get_time(ncFileID%times(1),secs,days)
      write(msgstring,*)'earliest netCDF time (days, seconds) ',days,secs
      call error_handler(E_DBG,'nc_get_tindex',msgstring,source,revision,revdate)

      call error_handler(E_ERR,'nc_get_tindex', &
              'Model time precedes earliest netCDF time.', source,revision,revdate)
      timeindex = -2

   else if ( statetime < ncFileID%times(ncFileID%Ntimes) ) then  

      ! It is somewhere in the middle without actually matching an existing time.
      ! This is very bad.

      write(msgstring,*)'model time does not match any netCDF time.'
      call error_handler(E_DBG,'nc_get_tindex',msgstring,source,revision,revdate)
      write(msgstring,*)'model time (days, seconds) is ',days,secs
      call error_handler(E_DBG,'nc_get_tindex',msgstring,source,revision,revdate)

      BadLoop : do i = 1,ncFileId%Ntimes   ! just find times to print before exiting

         if ( ncFileId%times(i) > statetime ) then
            call get_time(ncFileID%times(i-1),secs,days)
            write(msgstring,*)'preceding netCDF time (days, seconds) ',days,secs
            call error_handler(E_DBG,'nc_get_tindex',msgstring,source,revision,revdate)

            call get_time(ncFileID%times(i),secs,days)
            write(msgstring,*)'subsequent netCDF time (days, seconds) ',days,secs
            call error_handler(E_ERR,'nc_get_tindex',msgstring,source,revision,revdate)
            timeindex = -3
            exit BadLoop
         endif

      enddo BadLoop

   else ! we must need to append ... 

      timeindex = nc_append_time(ncFileID, statetime)

      write(msgstring,'(''appending model time (d,s) ('',i8,i5,'') as index '',i6, '' in ncFileID '',i10)') &
          days,secs,timeindex,my_ncid
      call error_handler(E_DBG,'nc_get_tindex',msgstring,source,revision,revdate)

   endif
   
endif

end function nc_get_tindex

!-------------------------------------------------------------------------------
!> The current time is appended to the "time" coordinate variable.
!> The new length of the "time" variable is returned.
!>
!> This REQUIRES that "time" is a coordinate variable AND it is the
!> unlimited dimension. If not ... bad things happen.

function nc_append_time(ncFileID, time) result(lngth)

type(netcdf_file_type), intent(inout) :: ncFileID
type(time_type), intent(in) :: time
integer                     :: lngth

integer  :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer  :: TimeVarID
integer  :: secs, days, my_ncid
real(digits12) :: realtime         ! gets promoted to nf90_double ...

character(len=NF90_MAX_NAME)          :: varname
integer                               :: xtype, ndims, nAtts
integer, dimension(NF90_MAX_VAR_DIMS) :: dimids

type(time_type), allocatable, dimension(:) :: temptime   ! only to reallocate mirror
real(digits12),  allocatable, dimension(:) :: tempRtime  ! only to reallocate mirror

lngth = -1 ! assume a bad termination

my_ncid = ncFileID%ncid

call nc_check(NF90_Inquire(my_ncid, nDimensions, nVariables, nAttributes, unlimitedDimID), &
              'nc_append_time', 'inquire '//ncFileID%fname)
call nc_check(NF90_Inq_Varid(my_ncid, "time", TimeVarID), 'nc_append_time', 'inq_varid time')
call nc_check(NF90_Inquire_Variable(my_ncid, TimeVarID, varname, xtype, ndims, dimids, nAtts), &
             'nc_append_time', 'inquire_variable time')

if ( ndims /= 1 ) call error_handler(E_ERR,'nc_append_time', &
           '"time" expected to be rank-1',source,revision,revdate)

if ( dimids(1) /= unlimitedDimID ) call error_handler(E_ERR,'nc_append_time', &
           'unlimited dimension expected to be slowest-moving',source,revision,revdate)

! make sure the mirror and the netcdf file are in sync
call nc_check(NF90_Inquire_Dimension(my_ncid, unlimitedDimID, varname, lngth ), &
           'nc_append_time', 'inquire_dimension unlimited')

if (lngth /= ncFileId%Ntimes) then
   write(msgstring,*)'netCDF file has length ',lngth,' /= mirror has length of ',ncFileId%Ntimes
   call error_handler(E_ERR,'nc_append_time', &
           'time mirror and netcdf file time dimension out-of-sync', &
           source,revision,revdate,text2=msgstring)
endif

! make sure the time mirror can handle another entry.
if ( lngth == ncFileID%NtimesMAX ) then   

   write(msgstring,*)'doubling mirror length of ',lngth,' of ',ncFileID%fname
   call error_handler(E_DBG,'nc_append_time',msgstring,source,revision,revdate)

   allocate(temptime(ncFileID%NtimesMAX), tempRtime(ncFileID%NtimesMAX)) 
   temptime  = ncFileID%times            ! preserve
   tempRtime = ncFileID%rtimes           ! preserve

   deallocate(ncFileID%times, ncFileID%rtimes)

   ncFileID%NtimesMAX = 2 * ncFileID%NtimesMAX  ! double length of exising arrays

   allocate(ncFileID%times(ncFileID%NtimesMAX), ncFileID%rtimes(ncFileID%NtimesMAX) )

   ncFileID%times(1:lngth)  = temptime    ! reinstate
   ncFileID%rtimes(1:lngth) = tempRtime   ! reinstate

   deallocate(temptime, tempRtime)

endif

call get_time(time, secs, days)         ! get time components to append
realtime = days + secs/86400.0_digits12 ! time base is "days since ..."
lngth           = lngth + 1             ! index of new time 
ncFileID%Ntimes = lngth                 ! new working length of time mirror

call nc_check(nf90_put_var(my_ncid, TimeVarID, realtime, start=(/ lngth /) ), &
           'nc_append_time', 'put_var time')

ncFileID%times( lngth) = time
ncFileID%rtimes(lngth) = realtime

write(msgstring,*)'ncFileID (',my_ncid,') : ',trim(adjustl(varname)), &
         ' (should be "time") has length ',lngth, ' appending t= ',realtime
call error_handler(E_DBG,'nc_append_time',msgstring,source,revision,revdate)

end function nc_append_time


!------------------------------------------------------------------
!> set handle netcdf_file_type

subroutine set_netcdf_file_type(file_handle, ncFileInfo)

type(file_info_type),   intent(inout) :: file_handle
type(netcdf_file_type), intent(in)    :: ncFileInfo

file_handle%stage_metadata%ncFileID = ncFileInfo

end subroutine set_netcdf_file_type


!------------------------------------------------------------------
!> return handle netcdf_file_type

subroutine get_netcdf_file_type(file_handle, ncFileInfo)

type(file_info_type),   intent(in)    :: file_handle
type(netcdf_file_type), intent(inout) :: ncFileInfo

ncFileInfo = file_handle%stage_metadata%ncFileID

end subroutine get_netcdf_file_type

!------------------------------------------------------------------
!------------------------------------------------------------------

!> @}
end module state_space_diag_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
