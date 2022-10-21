! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module direct_netcdf_mod

!> \defgroup direct_netcdf_mod direct_netcdf_mod
!> @{ 
!>
!> Netcdf IO for a domain.
!> The idea is to be generic rather than having a converter for each module.
!> Also aim to go to the filesystem only once, e.g.
!> \verbatim
!>      wrf => dart => wrf
!> \endverbatim
!> rather than
!> \verbatim
!>      wrf => wrf_to_dart => dart => dart_to_wrf => wrf
!> \endverbatim
!>
!> Every task needs the dimensions of each state variable to calculate
!> what it is going to recieve in the IO transpose (this is calculated in add_domain )
!> \par Aim of the limited transpose:
!>
!>To limit how much of the state vector you can read at once.
!> You can limit the transpose by memory using <code>buffer_state_io</code>.
!>
!>What you (potentially) gain from this:
!>
!>* Don't have to have the whole state vector.
!>* Don't have to use a parallel IO library.
!>
!>If buffer_state_io is false you have the regular transpose + IO, except:
!>  1. You are reading directly from a netcdf file, not a dart state vector file.
!>  2. You only transpose the copies that are being written/read.
!>
!> This code has multiple places where round-robin layout of state onto task is assumed.
!> These are in the section labelled: \n
!> \verbatim
!>    !--------------------------------------------------------
!>    ! Routines that are making the assumption that the ensemble
!>    ! distribution is round-robin (distribution type 1)
!>    !--------------------------------------------------------
!> \endverbatim
!>
!> read_transpose() is the read routine.
!> transpose_write() is the write routine.
!>
!> Note dart_index is an inout variable to read_transpose() and transpose_write().
!> This is making the assumption that the calling code is using dart_index in the following way:
!>  * dart_index going in to the subroutines is where the domain starts in the state vector.
!>  * dart_index coming out of the subroutines is where the domain ends.
!>  * The domain is contiguous in the state_vector.

use types_mod,            only : r4, r8, i4, i8, MISSING_R8, MISSING_R4, MISSING_I, &
                                 digits12, metadatalength

use options_mod,          only : get_missing_ok_status

use ensemble_manager_mod, only : ensemble_type, map_pe_to_task, map_task_to_pe, &
                                 all_copies_to_all_vars, all_vars_to_all_copies, &
                                 get_copy_owner_index, get_copy, get_ensemble_time, &
                                 allocate_vars

use time_manager_mod,     only : time_type, get_time, get_calendar_type, &
                                 THIRTY_DAY_MONTHS, JULIAN, GREGORIAN, NOLEAP, &
                                 operator(<), operator(>), operator(+), &
                                 operator(-), operator(/), operator(*), &
                                 operator(==), operator(/=)

use utilities_mod,        only : error_handler, file_to_text, &
                                 find_textfile_dims, file_exist, &
                                 E_MSG, E_ALLMSG, E_ERR, E_DBG, E_WARN

use netcdf_utilities_mod, only : nc_check

use mpi_utilities_mod,    only : task_count, send_to, receive_from, my_task_id, &
                                 broadcast_flag

use state_structure_mod,  only : get_num_variables, get_sum_variables,  &
                                 get_sum_variables_below, get_dim_length, &
                                 get_variable_name, get_io_clamping_maxval,   &
                                 get_io_clamping_minval, do_io_clamping,         &
                                 get_io_num_dims, get_dim_lengths,            &
                                 get_variable_size, get_io_num_unique_dims,   &
                                 get_io_unique_dim_name, get_dim_name,        &
                                 get_io_unique_dim_length, get_scale_factor, &
                                 set_var_id, get_domain_size, do_io_update, &
                                 get_units, get_long_name, get_short_name, &
                                 get_has_missing_value, get_FillValue, &
                                 get_missing_value, get_add_offset, get_xtype, &
                                 get_has_FillValue, &
                                 get_index_start, get_index_end , get_num_dims, &
                                 create_diagnostic_structure, &
                                 end_diagnostic_structure, &
                                 has_unlimited_dim

use io_filenames_mod,     only : get_restart_filename, inherit_copy_units, &
                                 stage_metadata_type, get_file_description, &
                                 get_copy_name, copy_is_clamped, query_read_copy, &
                                 query_write_copy, force_copy_back, file_info_type, &
                                 netcdf_file_type, READ_COPY, WRITE_COPY, &
                                 noutput_state_variables

use assim_model_mod,      only : get_model_size, read_model_time, write_model_time

!>@todo FIXME : should move to assim_model_mod.f90
use model_mod,            only : nc_write_model_atts

use typesizes

use netcdf

implicit none
private

public :: read_transpose,            &
          transpose_write,           &
          initialize_single_file_io, &          
          finalize_single_file_io,   &
          read_single_file,          &
          write_single_file,         &
          write_augmented_state,     &
          read_variables,            &
          nc_get_num_times

character(len=*), parameter :: source = 'direct_netcdf_mod.f90'

! only a single MPI Task reads and writes reads the state variable,
! when using single_file_{input,output} = .true.
integer, parameter :: SINGLE_IO_TASK_ID = 0 

! module global variables
integer :: ret 
character(len=512) :: msgstring, msgstring2

!>@todo FIXME:
!> this should be in a namelist somewhere and passed in from
!> a higher level routine.
!>
!> we had this in the ROMS branch but no one can remember why.
!> we think it's because we copied a template file into place
!> for filter to overwrite and it might not have had the right
!> analysis time in the file.  it's true that filter doesn't
!> change the time in the file, but when we're doing direct
!> updates of a netcdf file it's also true that we rarely
!> update the file we read from; in case of an error you've
!> destroyed your input needed to rerun the job.
!> (but we don't remember for sure if this was the reason.)

logical :: overwrite_time_in_output_file = .false.

contains

!=================================================
! PUBLIC ROUTINES
!=================================================

!=================================================
! mutiple file IO
!=================================================

!-------------------------------------------------
!> Read and transpose variables in cyclic-cyclic distribution (round robbin)

subroutine read_transpose(state_ens_handle, name_handle, domain, dart_index, read_single_vars)

type(ensemble_type),       intent(inout) :: state_ens_handle !< Ensemble handle
type(stage_metadata_type), intent(in)    :: name_handle      !< Name handle
integer,                   intent(in)    :: domain           !< Which domain to read
integer(i8),               intent(inout) :: dart_index       !< This is for multiple domains
logical,                   intent(in)    :: read_single_vars !< Read one variable at a time

if (task_count() == 1) then
   call read_transpose_single_task(state_ens_handle, name_handle, domain, dart_index)
else
   call read_transpose_multi_task( state_ens_handle, name_handle, domain, dart_index, read_single_vars)
endif

end subroutine read_transpose


!-------------------------------------------------
!> This code transposes and writes out the state vector copies
!> This can either be done on a single task or with multiple tasks

subroutine transpose_write(state_ens_handle, name_handle, domain, &
                           dart_index, write_single_vars, write_single_precision)

type(ensemble_type),       intent(inout) :: state_ens_handle
type(stage_metadata_type), intent(in)    :: name_handle
integer,                   intent(in)    :: domain
integer(i8),               intent(inout) :: dart_index
logical,                   intent(in)    :: write_single_vars !< write one variable at a time
logical,                   intent(in)    :: write_single_precision

if (task_count() == 1) then
   call transpose_write_single_task(state_ens_handle, name_handle, domain, &
                                    dart_index, write_single_precision)
else
   call transpose_write_multi_task(state_ens_handle, name_handle, domain, &
                                   dart_index, write_single_vars, write_single_precision)
endif

end subroutine transpose_write

!=================================================
! single file IO
!=================================================

!-------------------------------------------------
!> Creates a template file for single file io
!> Calls the model for any model specific attributes to be written

subroutine initialize_single_file_io(ens_handle, file_handle)

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
!  no_calendar No calendar. 
!
! location is another one ...

!@todo FIXME : need to have this work for multiple domains

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
character(len=256) :: fname
integer :: icopy, ivar, ret, ens_size, num_output_ens, domain

if (my_task_id() == 0) then
   if(.not. byteSizesOK()) then
       call error_handler(E_ERR,'initialize_single_file_io', &
      'Compiler does not support required kinds of variables.',source) 
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
   call nc_check(nf90_create(path=fname, cmode=createmode, ncid=ncFileID%ncid), &
                 'initialize_single_file_io', 'create '//trim(fname))

   ncFileID%ncid = ncFileID%ncid
   
   ! Define the dimensions
   call nc_check(nf90_def_dim(ncid=ncFileID%ncid, &
                 name="metadatalength", len=metadata_length, dimid=MetaDataDimID), &
                 'initialize_single_file_io', 'def_dim metadatalength '//trim(fname))
   
   call nc_check(nf90_def_dim(ncid=ncFileID%ncid, &
                 name="member", len=num_output_ens, dimid=MemberDimID), &
                 'initialize_single_file_io', 'def_dim member '//trim(fname))
   
   call nc_check(nf90_def_dim(ncid=ncFileID%ncid, name="time", &
                 len=nf90_unlimited, dimid=TimeDimID), &
                 'initialize_single_file_io', 'def_dim time '//trim(fname))

   !----------------------------------------------------------------------------
   ! Find dimensions of namelist file ... will save it as a variable.
   !----------------------------------------------------------------------------
   
   ! All DART programs require input.nml, so it is unlikely this can fail, but
   ! still check and in this case, error out if not found.
   call find_textfile_dims("input.nml", nlines, linelen)
   if (nlines <= 0 .or. linelen <= 0) then
      call error_handler(E_MSG,'initialize_single_file_io', &
              'cannot open/read input.nml to save in diagnostic file', source)
   endif
   
   allocate(textblock(nlines))
   textblock = ''
   
   call nc_check(nf90_def_dim(ncid=ncFileID%ncid, &
                 name="NMLlinelen", len=LEN(textblock(1)), dimid=linelenDimID), &
                 'initialize_single_file_io', 'def_dim NMLlinelen '//trim(fname))
   
   call nc_check(nf90_def_dim(ncid=ncFileID%ncid, &
                 name="NMLnlines", len=nlines, dimid=nlinesDimID), &
                 'initialize_single_file_io', 'def_dim NMLnlines '//trim(fname))
   
   !----------------------------------------------------------------------------
   ! Write Global Attributes 
   !----------------------------------------------------------------------------
   
   call nc_check(nf90_put_att(ncFileID%ncid, NF90_GLOBAL, "title", fname), &
                 'initialize_single_file_io', 'put_att title '//trim(fname))
   call nc_check(nf90_put_att(ncFileID%ncid, NF90_GLOBAL, "assim_model_source", source ), &
                 'initialize_single_file_io', 'put_att assim_model_source '//trim(fname))
   
   !    Metadata for each Copy
   call nc_check(nf90_def_var(ncid=ncFileID%ncid, name="MemberMetadata", xtype=nf90_char,    &
                 dimids=(/ MetaDataDimID, MemberDimID /),  varid=metadataVarID), &
                 'initialize_single_file_io', 'def_var MemberMetadata')
   call nc_check(nf90_put_att(ncFileID%ncid, metadataVarID, "long_name",       &
                 "Metadata for each copy/member"), 'initialize_single_file_io', 'put_att long_name')
   
   !    Input namelist 
   call nc_check(nf90_def_var(ncid=ncFileID%ncid,name="inputnml", xtype=nf90_char,    &
                 dimids = (/ linelenDimID, nlinesDimID /),  varid=nmlVarID), &
                 'initialize_single_file_io', 'def_var inputnml')
   call nc_check(nf90_put_att(ncFileID%ncid, nmlVarID, "long_name",       &
                 "input.nml contents"), 'initialize_single_file_io', 'put_att input.nml')
   
   !    Time -- the unlimited dimension
   call nc_check(nf90_def_var(ncFileID%ncid, name="time", xtype=nf90_double, dimids=TimeDimID, &
                 varid =TimeVarID), 'initialize_single_file_io', 'def_var time' )
   ret = nc_write_calendar_atts(ncFileID, TimeVarID)     ! comes from time_manager_mod
   if ( ret /= 0 ) then
      write(msgstring, *)'nc_write_calendar_atts  bombed with error ', ret
      call error_handler(E_MSG,'initialize_single_file_io',msgstring,source)
   endif
   
   ! Create the time "mirror" with a static length. There is another routine
   ! to increase it if need be. For now, just pick something.
   ncFileID%Ntimes    = 0
   ncFileID%NtimesMAX = 1000
   allocate(ncFileID%rtimes(ncFileID%NtimesMAX), ncFileID%times(ncFileID%NtimesMAX) )
   
   !----------------------------------------------------------------------------
   ! Leave define mode so we can fill
   !----------------------------------------------------------------------------
   
   call nc_check(nf90_enddef(ncFileID%ncid), 'initialize_single_file_io', 'enddef '//trim(fname))

   call nc_check(nf90_sync(ncFileID%ncid), 'initialize_single_file_io', 'sync '//trim(fname))               
   !----------------------------------------------------------------------------
   ! Define the model-specific components
   !----------------------------------------------------------------------------
   ! for single file io, we are assuming a single domain, no lower order models
   ! have multiple domains.

   domain = 1
   call nc_write_model_atts( ncFileID%ncid, domain)

   if ( .not. local_model_mod_will_write_state_variables ) then
      call write_model_attributes(ncFileID, MemberDimID, TimeDimID, file_handle)
   endif

   !----------------------------------------------------------------------------
   ! Create variables and attributes.
   ! The locations are part of the model (some models have multiple grids).
   ! They are written by model_mod:nc_write_model_atts
   !----------------------------------------------------------------------------

   ens_size = ens_handle%num_copies - ens_handle%num_extras
   do icopy = ens_size+1, ens_handle%num_copies 
      if ( file_handle%stage_metadata%io_flag(icopy) == WRITE_COPY ) then
         call write_extra_attributes( ncFileID, TimeDimID, file_handle, icopy )
      endif
   enddo 
   

   
   !----------------------------------------------------------------------------
   ! Fill the coordinate variables.
   ! Write the input namelist as a netCDF variable.
   ! The time variable is filled as time progresses.
   !----------------------------------------------------------------------------
   
   call nc_check(nf90_put_var(ncFileID%ncid, metadataVarID, state_meta ), &
                 'initialize_single_file_io', 'put_var MetaDataVarID')
    
   call file_to_text("input.nml", textblock)
   
   call nc_check(nf90_put_var(ncFileID%ncid, nmlVarID, textblock ), &
                 'initialize_single_file_io', 'put_var nmlVarID')
   
   deallocate(textblock)
   
   !----------------------------------------------------------------------------
   ! sync to disk, but leave open
   !----------------------------------------------------------------------------
   
   call nc_check(nf90_sync(ncFileID%ncid), 'initialize_single_file_io', 'sync '//trim(fname))               
   !----------------------------------------------------------------------------
   ! sync again, but still leave open
   !----------------------------------------------------------------------------
   
   call nc_check(nf90_sync(ncFileID%ncid), 'initialize_single_file_io', 'sync '//trim(fname))               
   file_handle%stage_metadata%ncFileID = ncFileID
 
endif

! Broadcast the value of model_mod_will_write_state_variables to every task
! This keeps track of whether the model_mod or dart code will write state_variables.
call broadcast_flag(local_model_mod_will_write_state_variables, 0)
ncFileID%model_mod_will_write_state_variables = local_model_mod_will_write_state_variables

file_handle%singlefile_initialized = .true.

end subroutine initialize_single_file_io


!-------------------------------------------------------------------------------
!>

subroutine finalize_single_file_io(file_handle)

type(file_info_type), intent(in) :: file_handle

type(netcdf_file_type) :: ncFileID
integer :: ierr

ncFileID = file_handle%stage_metadata%ncFileID

if (my_task_id()==0) then
   ierr = nf90_close(ncFileID%ncid)
   call nc_check(ierr, 'finalize_single_file_io: nf90_close', ncFileID%fname)
   if(associated(ncFileID%rtimes)) deallocate(ncFileID%rtimes, ncFileID%times )
endif

ncFileID%fname     = "notinuse"
ncFileID%ncid      = -1
ncFileID%Ntimes    = -1
ncFileID%NtimesMax = -1

call end_diagnostic_structure()

end subroutine finalize_single_file_io


!-------------------------------------------------------
!> read a single netcdf file containing all of the members
!> and possibly inflation information

subroutine read_single_file(state_ens_handle, file_handle, use_time_from_file, mtime, pert_from_single_copy)

type(ensemble_type),  intent(inout) :: state_ens_handle      !! ensemble handle to store data
type(file_info_type), intent(in)    :: file_handle           !! file handle for file names
logical,              intent(in)    :: use_time_from_file    !! read time from file
type(time_type),      intent(inout) :: mtime                 !! external time
logical, optional,    intent(in)    :: pert_from_single_copy !! reading single file and perturbing

! NetCDF IO variables
integer                               :: my_ncid, varid, TimeDimID, MemDimID, ret 
character     (len=NF90_MAX_NAME    ) :: fname, varname, copyname
integer, dimension(NF90_MAX_VAR_DIMS) :: dim_lengths
integer, dimension(NF90_MAX_VAR_DIMS) :: dim_start_point

integer :: icopy, ivar, domain  
integer :: ens_size, extra_size, time_size, ndims
integer :: my_pe, recv_pe, recv_start, recv_end, start_rank
integer :: send_start, send_end
logical :: do_perturb, is_sender, is_receiver, is_extra_copy
integer(i8) :: var_size, elm_count, start_pos, end_pos, start_point

real(r8), allocatable :: var_block(:)

do_perturb = .false.
if ( present(pert_from_single_copy) ) then
   do_perturb = pert_from_single_copy
endif

! if perturbing we only need to read the first ensemble member
if ( do_perturb ) ens_size = 1

! grab ensemble size and number of extra copies
ens_size   = state_ens_handle%num_copies - state_ens_handle%num_extras
extra_size = state_ens_handle%num_extras

!>@todo : only a single domain for single file read supported. need           
!>        to consider case for multiple domains.
domain = 1 

fname = get_restart_filename(file_handle%stage_metadata, 1, domain)

!>@todo Check time consistency across files? This is assuming they are consistent.
! read time from input file if time not set in namelist
if( use_time_from_file ) mtime = read_model_time(fname)

state_ens_handle%time = mtime

ret = nf90_open(fname, NF90_NOWRITE, my_ncid)
call nc_check(ret, 'read_single_file: nf90_open', fname)

ret = nf90_inq_dimid(my_ncid, "member", MemDimID)
if (ret /= NF90_NOERR) then
   call error_handler(E_ERR,'direct_netcdf_mod:', &
         'If using single_file_in/single_file_out = .true. ', &
          source, text2='you must have a member dimension in your input/output file.')
endif

ret = nf90_inq_dimid(my_ncid, "time", TimeDimID)
call nc_check(ret, 'read_single_file', 'inq_varid time : '//trim(fname))

ret = nf90_inquire_dimension(my_ncid, TimeDimID, len=time_size) 
call nc_check(ret, 'read_single_file', 'inquire_dimension time '//trim(fname))

! mpi task variables
my_pe      = my_task_id()
is_sender  = (my_pe == SINGLE_IO_TASK_ID)

! recv_* and send_* are PE's that variables are sent and received
call get_pe_loops(ens_size, recv_start, recv_end, send_start, send_end)

call check_singlefile_member_info(my_ncid, fname, ens_size, do_perturb)

COPY_LOOP: do icopy = 1, ens_size+extra_size

   ! only SINGLE_IO_TASK_ID reads and distributes data 
   ! extra copies into ensemble handle, this could include
   !   {variable}_{mean,sd}
   !   {variable}_priorinf_{mean,sd}
   !   {variable}_postinf_{mean,sd}

   if ( file_handle%stage_metadata%io_flag(icopy) /= READ_COPY ) cycle

   is_extra_copy  = (icopy >  ens_size)

   ! check that copy infomation is valid

   ! starting position in the copies array
   start_pos  = 1

   VAR_LOOP: do ivar = 1, get_num_variables(domain)

      var_size = get_variable_size(domain, ivar)
      varname  = get_variable_name(domain, ivar)
   
      ! work out netcdf dimension id's and lenghts
      call get_dimension_info(icopy, domain, ivar, time_size, is_extra_copy, &
                              ndims, dim_start_point, dim_lengths)

      ! read variable from SINGLE_IO_TASK_ID
      if ( is_sender ) then

         allocate(var_block(var_size))

         ! append the copyname to the variable (ex. ps_mean)
         if ( is_extra_copy ) then
            copyname = get_copy_name(file_handle,icopy)
            write(varname,'(a,"_",a)') trim(varname), trim(copyname)
         endif

         ! get the variable id from the copy we are reading
         ret = nf90_inq_varid(my_ncid, varname, varid)
         call nc_check(ret, 'read_single_file', 'inq_varid '//trim(varname)//' : '//trim(fname))

         ! if it is an ensemble member we are just reading one member from one variable 
         ret = nf90_get_var(my_ncid, varid, var_block,  &
                            count=dim_lengths(1:ndims), &
                            start=dim_start_point(1:ndims))
         call nc_check(ret, 'read_single_file', 'get_var '//trim(varname)//' : '//trim(fname))

      endif
      
      !>@todo FIXME : just reading and distributing one variable at a time.  
      !>              Should probably do this in a block?
      !start_rank = get_start_rank(start_var, domain)

      ! this is where we left off writing variables in the case of multiple tasks
      start_rank = get_start_rank(ivar, domain)
      
      RECEIVING_PE_LOOP: do recv_pe = 0, task_count()-1

         ! work out the start in var_block on the receiving pe.
         ! elm_count is the sum of variables that have been read.
         elm_count = num_elements_on_pe(recv_pe, start_rank, var_size) 
         end_pos   = start_pos + elm_count - 1 ! ending position in the copies array

         ! work out the start in var_block corresponding to the receiving pe
         start_point = find_start_point(recv_pe, start_rank) 
         is_receiver = (my_pe == recv_pe)

         if ( is_receiver ) then 
            if ( is_sender ) then ! just copy directly into copies

               ! non-contiguous array. directly copying into the copies array
               ! since no send required.
               state_ens_handle%copies(icopy, start_pos:end_pos ) = &
               var_block(start_point:elm_count*task_count():task_count())

            else ! post receive

               call wait_to_receive(state_ens_handle, SINGLE_IO_TASK_ID, icopy, &
                                    start_pos, end_pos)
            endif ! is_sender

            ! update starting point in the state vector
            start_pos = start_pos + elm_count

         elseif ( is_sender ) then ! send variables to receivers
            
            call send_to_waiting_task(state_ens_handle, recv_pe, start_point, &
                                      elm_count, var_size, var_block)

         endif ! is_receiver

      enddo RECEIVING_PE_LOOP

      if( is_sender ) then
         deallocate(var_block)
      endif

   enddo VAR_LOOP
enddo COPY_LOOP

ret = nf90_close(my_ncid)
call nc_check(ret, 'read_single_file: nf90_close', fname)

end subroutine read_single_file


!-----------------------------------------------------------
!> write all variable to a single file including all member,
!> and optionally inflation, mean, sd

subroutine write_single_file(ens_handle, file_handle)

type(ensemble_type),         intent(inout) :: ens_handle
type(file_info_type),        intent(inout) :: file_handle


! Local variables
integer                :: member_index, copy_index, my_ncid
integer                :: num_output_ens, ens_size
integer(i8)            :: model_size
character(len=128)     :: copyname
type(time_type)        :: curr_ens_time
type(netcdf_file_type) :: ncFileID
real(r8), allocatable  :: temp_ens(:)
integer                :: timeindex
integer                :: is1, id1

! Assumes that mean and spread have already been computed
! make sure vars is up-to-date
call allocate_vars(ens_handle)
call all_copies_to_all_vars(ens_handle)

model_size = get_model_size()

! SINGLE_IO_TASK_ID needs some space
if (my_task_id() == SINGLE_IO_TASK_ID) then
   allocate(temp_ens(model_size))
else
   allocate(temp_ens(1))
endif

 ! SINGLE_IO_TASK_ID writes out all files
if (my_task_id() == SINGLE_IO_TASK_ID) then
  
   curr_ens_time = ens_handle%current_time
   ncFileID = file_handle%stage_metadata%ncFileID
   
   timeindex = nc_get_tindex(ncFileID, curr_ens_time)
   
   if ( timeindex < 0 ) then
      call get_time(curr_ens_time,is1,id1)
      write(msgstring,*)'model time (d,s)',id1,is1,' not in ',ncFileID%fname
      write(msgstring,'(''model time (d,s) ('',i8,i5,'') is index '',i6, '' in ncFileID '',i10)') &
             id1,is1,timeindex,ncFileID%ncid
      call error_handler(E_ERR,'write_model_variables', msgstring, source)
   endif
   
   call get_time(curr_ens_time,is1,id1)
   write(msgstring,'(''model time (d,s) ('',i8,i5,'') is index '',i6, '' in ncFileID '',i10)') &
          id1,is1,timeindex,ncFileID%ncid
   call error_handler(E_DBG,'write_model_variables', msgstring, source)
   
   my_ncid = ncFileID%ncid
endif 

num_output_ens = noutput_state_variables(file_handle)

! Output ensemble members
do member_index = 1, num_output_ens
   call get_copy(map_task_to_pe(ens_handle, 0), ens_handle, member_index, temp_ens)
   if(my_task_id() == SINGLE_IO_TASK_ID) &
      call write_model_variables(ncFileID,  temp_ens, member_index, curr_ens_time, timeindex)
enddo

! Output Extras
ens_size = ens_handle%num_copies - ens_handle%num_extras

do copy_index = ens_size+1, ens_handle%num_copies
   if ( file_handle%stage_metadata%io_flag(copy_index) == WRITE_COPY ) then
      copyname = get_copy_name(file_handle, copy_index)
      call get_copy(map_task_to_pe(ens_handle, 0), ens_handle, copy_index, temp_ens)
      if(my_task_id() == SINGLE_IO_TASK_ID) &
         call write_extra_variables(ncFileID, temp_ens, copyname, curr_ens_time, timeindex)
   endif
enddo
   
if (my_task_id() == SINGLE_IO_TASK_ID) then ! SINGLE_IO_TASK_ID writes out all files
   !>@ todo FIXME : This sync is not necessary but ensures that all of the variable
   !>  information gets written out if filter happens to crash in the middle of a
   !> run. This may slow down the code for longer runs.
   call nc_check(nf90_sync(my_ncid), 'write_single_file', 'nf90_sync')
   
   file_handle%stage_metadata%ncFileID = ncFileID
endif 

deallocate(temp_ens)

end subroutine write_single_file


!-----------------------------------------------------------
!> insert the mean and sd into the input file if the user requests
!> this assumes that the files has already been closed and the
!> variables read in.

subroutine write_augmented_state(ens_handle, file_handle)

type(ensemble_type),  intent(inout) :: ens_handle
type(file_info_type), intent(inout) :: file_handle

! Local variables
integer                :: copy_index, ens_size, domain
integer                :: TimeDimID, my_ncid, ret
integer(i8)            :: model_size
type(time_type)        :: curr_ens_time
type(netcdf_file_type) :: ncFileID
real(r8), allocatable  :: temp_ens(:)
character(len=256) :: fname, copyname

! assumes that mean and spread have already been computed
! make sure vars is up-to-date
call allocate_vars(ens_handle)
call all_copies_to_all_vars(ens_handle)

model_size = get_model_size()

! SINGLE_IO_TASK_ID needs some space
if (my_task_id() == SINGLE_IO_TASK_ID) then
   allocate(temp_ens(model_size))
else
   allocate(temp_ens(1))
endif

curr_ens_time = ens_handle%current_time

!>@todo : only a single domain for single file read supported.
!>        need to consider case for multiple domains.

domain = 1

if (my_task_id() == SINGLE_IO_TASK_ID) then
   fname = get_restart_filename(file_handle%stage_metadata, 1, domain)
   
   ret = nf90_open(fname, NF90_WRITE, my_ncid)
   call nc_check(ret, 'write_augmented_state: nf90_open', fname)
   
   ret = nf90_inq_dimid(my_ncid, "time", TimeDimID)
   call nc_check(ret, 'write_augmented_state', 'inq_varid time : '//trim(fname))
   
   ncFileID%ncid  = my_ncid
   ncFileID%fname = fname
endif

! Output Mean and SD

ens_size = ens_handle%num_copies - ens_handle%num_extras

! all tasks must execute this loop in order to get_copy from all tasks
do copy_index = ens_size+1, ens_handle%num_copies
   if ( file_handle%stage_metadata%io_flag(copy_index) == WRITE_COPY ) then
      call get_copy(map_task_to_pe(ens_handle, 0), ens_handle, copy_index, temp_ens)
      if(my_task_id() == SINGLE_IO_TASK_ID) then
         copyname = get_copy_name(file_handle, copy_index)
         call write_extra_attributes( ncFileID, TimeDimID, file_handle, copy_index)
         call write_extra_variables(  ncFileID, temp_ens,  copyname, curr_ens_time, 1)
      endif
   endif
enddo

if (my_task_id() == SINGLE_IO_TASK_ID) then
   ret = nf90_close(my_ncid)
   call nc_check(ret, 'write_augmented_state: nf90_close', fname)
endif

deallocate(temp_ens)

end subroutine write_augmented_state


!-------------------------------------------------------------------------------
!> Read in variables from start_var to end_var
!> Read the latest time slice in the file
subroutine read_variables(ncfile_in, var_block, start_var, end_var, domain)

integer,  intent(in)    :: ncfile_in
real(r8), intent(inout) :: var_block(:)
integer,  intent(in)    :: start_var
integer,  intent(in)    :: end_var
integer,  intent(in)    :: domain

integer :: i
integer(i8) :: istart, iend
integer(i8) :: var_size
integer :: num_dims
integer, allocatable :: counts(:)
integer, allocatable :: slice_start(:) ! slice of variable
integer :: ret, var_id, unlim_dimID
logical :: missing_possible

missing_possible = get_missing_ok_status()

istart = 1

do i = start_var, end_var

   var_size = get_variable_size(domain, i)
   iend = istart + var_size - 1

   ! number of dimensions and length the variable
   num_dims = get_io_num_dims(domain, i)
   allocate(counts(num_dims))
   allocate(slice_start(num_dims))
   counts(:) = 1

   slice_start(:) = 1 ! default to read all dimensions start at 1

   if (has_unlimited_dim(domain)) then

      counts(num_dims) = 1 ! one slice of unlimited dimesion
      counts(1:num_dims-1) = get_dim_lengths(domain, i) ! the state
      
      ! read latest time slice - hack to get started with tiegcm
      ! not sure if it will always be the last time slice
      ret = nf90_inquire(ncfile_in, unlimitedDimID=unlim_dimID)
      call nc_check(ret, 'read_variables: nf90_inquire', 'unlimitedDimID')

      if (unlim_dimID /= -1) then ! unlimited dimension exists
         ret = nf90_inquire_dimension(ncfile_in, unlim_dimID, len=slice_start(num_dims))
         call nc_check(ret, 'read_variables: nf90_inquire_dimension', 'unlimitedDim length')
         if (slice_start(num_dims) == 0) slice_start(num_dims) = 1 ! newly created file
      else  ! file does not have an unlimited dimension because it was created by DART
          slice_start(num_dims) = 1
      endif

   else

      counts(1:get_num_dims(domain,i)) = get_dim_lengths(domain, i) ! the state
   endif

   ret = nf90_inq_varid(ncfile_in, get_variable_name(domain, i), var_id)
   call nc_check(ret, 'read_variables: nf90_inq_varid',trim(get_variable_name(domain,i)) )

   ret = nf90_get_var(ncfile_in, var_id, var_block(istart:iend), count=counts, start=slice_start)
   call nc_check(ret, 'read_variables: nf90_get_var',trim(get_variable_name(domain,i)) )

   if (missing_possible) call set_dart_missing_value(var_block(istart:iend), domain, i)

   istart = istart + var_size

   deallocate(counts, slice_start)

enddo

end subroutine read_variables


!=================================================
! HELPER ROUTINES
!=================================================

!-------------------------------------------------
!> Single processor version of read_transpose.  Reads ens_size whole vectors from
!> netcdf files and fills up a row of %copies for each file.

subroutine read_transpose_single_task(state_ens_handle, name_handle, domain, dart_index)

type(ensemble_type),      intent(inout) :: state_ens_handle
type(stage_metadata_type), intent(in)   :: name_handle
integer,                   intent(in)   :: domain
integer(i8),              intent(inout) :: dart_index !< This is for multiple domains

real(r8), allocatable :: vector(:)

integer :: ncfile !< netcdf input file identifier
character(len=256) :: netcdf_filename

integer(i8) :: block_size, istart, iend
integer :: copy , start_var

istart     = dart_index ! position in state_ens_handle%vars
block_size = 0

! need to read into a tempory array, then fill up copies
allocate(vector(get_domain_size(domain)))

COPIES: do copy = 1, state_ens_handle%my_num_copies

   start_var = 1 ! read first variable first

   ! open netcdf file
   if (query_read_copy(name_handle, copy)) then
      netcdf_filename = get_restart_filename(name_handle, copy, domain)
      ret = nf90_open(netcdf_filename, NF90_NOWRITE, ncfile)
      call nc_check(ret, 'read_transpose_single_task: opening', netcdf_filename)
   endif

   block_size = get_domain_size(domain)

   iend = istart + block_size -1

   if (query_read_copy(name_handle, copy)) then
      call read_variables(ncfile, vector, 1, get_num_variables(domain), domain)
      ! close netcdf file
      ret = nf90_close(ncfile)
      call nc_check(ret, 'read_transpose_single_task: closing', netcdf_filename)
      state_ens_handle%copies(copy, istart:iend) = vector

   endif

enddo COPIES

! update starting point
istart = istart + block_size

dart_index = istart

deallocate(vector)

end subroutine read_transpose_single_task


!-------------------------------------------------
!> Single processor version of transpose write.  Takes copies array one row
!> at a time and writes copy to a netcdf file.

subroutine transpose_write_single_task(state_ens_handle, name_handle, domain, &
                     dart_index, write_single_precision)

type(ensemble_type),       intent(inout) :: state_ens_handle
type(stage_metadata_type), intent(in)    :: name_handle
integer,                   intent(in)    :: domain
integer(i8),               intent(inout) :: dart_index
logical,                   intent(in)    :: write_single_precision

! netcdf variables
integer :: ncfile_out
character(len=256) :: netcdf_filename_out

real(r8), allocatable :: vector(:)

integer(i8) :: block_size , istart, iend
integer :: copy , start_var, end_var
integer :: time_owner, time_owner_index
logical :: clamp_vars, force_copy
type(time_type) :: dart_time

! need to read into a tempory array to fill with one copies
allocate(vector(get_domain_size(domain)))

istart = dart_index ! position in state_ens_handle%vars
block_size = 0

! need to read into a temporary array, then fill up copies

COPIES: do copy = 1, state_ens_handle%my_num_copies

   start_var = 1 ! read first variable first

   ! open netcdf file
   if (query_write_copy(name_handle, copy)) then
      netcdf_filename_out = get_restart_filename(name_handle, copy, domain)

      if(file_exist(netcdf_filename_out)) then
         ret = nf90_open(netcdf_filename_out, NF90_WRITE, ncfile_out)
         call nc_check(ret, 'transpose_write: opening', trim(netcdf_filename_out))
         call nc_write_global_att_clamping(ncfile_out, copy, domain)

         if (overwrite_time_in_output_file) then
            call get_copy_owner_index(state_ens_handle, copy, time_owner, time_owner_index)
            call get_ensemble_time(state_ens_handle, time_owner_index, dart_time)
            call write_model_time(ncfile_out, dart_time)
         endif

      else ! create and open file
         !>@todo This is grabbing the time assuming the ensemble is var complete.
         !> Should we instead have all copies time in the ensemble handle?
         call get_copy_owner_index(state_ens_handle, copy, time_owner, time_owner_index)
         call get_ensemble_time(state_ens_handle, time_owner_index, dart_time)
         ncfile_out = create_and_open_state_output(name_handle, domain, copy, &
                                                   dart_time, write_single_precision)
         !>@todo if multiple domains exist in the same file, only the variables
         !>      from the first domain are created by create_and_open_state_output()
         !>      and since the file exists, the variables for the additional domains
         !>      never get defined in the netCDF file.
      endif
   endif

   block_size = get_domain_size(domain)
   iend       = istart + block_size -1

   if (query_write_copy(name_handle, copy)) then

      vector = state_ens_handle%copies(copy, istart:iend)

      ! for a single task the end var will always be the last element.
      ! do not need to limit memory since the entire state is all on
      ! a single processor.
      end_var = get_num_variables(domain)

      ! actual copy, may need clamping
      clamp_vars = copy_is_clamped(name_handle, copy)
      force_copy = force_copy_back(name_handle, copy)

      call write_variables(ncfile_out, vector, start_var, end_var, &
                           domain, clamp_vars, force_copy)

      ! close netcdf file
      ret = nf90_close(ncfile_out)
      call nc_check(ret, 'transpose_write closing', netcdf_filename_out)
   endif

enddo COPIES

! update starting point
istart = istart + block_size

dart_index = istart

deallocate(vector)

end subroutine transpose_write_single_task


!-------------------------------------------------
! Multiple tasks
!-------------------------------------------------
!> Read in variables from model restart file and transpose so that every processor
!> has all copies of a subset of state variables (fill state_ens_handle%copies)
!> Read and transpose data according to the memory limit imposed by
!> read_var_by_var.


subroutine read_transpose_multi_task(state_ens_handle, name_handle, domain, &
                dart_index, read_var_by_var)

type(ensemble_type),       intent(inout) :: state_ens_handle
type(stage_metadata_type), intent(in)    :: name_handle
integer,                   intent(in)    :: domain
integer(i8),               intent(inout) :: dart_index !< This is for multiple domains
logical,                   intent(in)    :: read_var_by_var !< Read one variable at a time 

integer :: start_var, end_var   !< start/end variables in a read block
integer :: start_rank           !< starting rank containg variable of interest
integer :: recv_start, recv_end !< start/end variables for receives
integer :: send_start, send_end !< start/end variables for sends
integer :: recv_pe, sending_pe  !< PEs sending and receiving data
integer(i8) :: elm_count        !< number of elements to send
integer(i8) :: block_size       !< number of state elements in a block
integer(i8) :: istart, iend     !< position in state vector copies array
integer :: ens_size             !< ensemble size
integer :: my_pe                !< task or pe?
integer :: ensemble_member      !< the ensmeble_member you are receiving.
integer :: my_copy              !< which copy a pe is reading, from 1 to ens_handle%num_copies
integer :: c                    !< copies_read loop index
integer(i8) :: start_point
integer :: copies_read
integer :: num_state_variables
integer :: dummy_loop
logical :: is_reader            !< pe is a reader or not
real(r8), allocatable :: var_block(:) !< for reading in variables

! netcdf variables
integer :: ncfile !< netcdf input file identifier
character(len=256) :: netcdf_filename !< different for each task

ens_size = state_ens_handle%num_copies ! have the extras, incase you need to read inflation restarts

my_pe = state_ens_handle%my_pe
num_state_variables = get_num_variables(domain)

! need to calculate RECEIVING_PE_LOOP start:end, group size, sending_pe start:end for each group.
call get_pe_loops(ens_size, recv_start, recv_end, send_start, send_end)

if (my_pe < ens_size) then
   is_reader = .true.
else
   is_reader = .false.
endif

copies_read = 0

COPIES: do c = 1, ens_size
   if (copies_read >= ens_size) exit

   ! what to do if a variable is larger than the memory limit?
   start_var = 1 ! read first variable first from the var_block
   istart    = dart_index ! position in state_ens_handle%copies

   my_copy = copies_read + my_pe + 1

   ! open netcdf file
   ! You have already opened this once to read the variable info. Should you just leave it open
   ! on the readers?
   if (is_reader) then

      if (query_read_copy(name_handle, my_copy)) then
         netcdf_filename = get_restart_filename(name_handle, my_copy, domain)
         ret = nf90_open(netcdf_filename, NF90_NOWRITE, ncfile)
         call nc_check(ret, 'read_transpose opening', netcdf_filename)
      endif

   endif

   ! Read the state variables into a buffer to distribute.
   ! If possible, reading all the variables into a single buffer 
   ! (i.e. a large block_size) is preferable to reading variables
   ! into multiple buffers. Huge DART states may require using
   ! the same buffer multiple times.

   VARIABLE_LOOP: do dummy_loop = 1, num_state_variables

      if (start_var > num_state_variables) exit VARIABLE_LOOP

      ! calculate how many variables will be read into one buffer
      if (read_var_by_var) then 
         end_var = start_var
      else
         end_var = num_state_variables
      endif

      block_size = get_sum_variables(start_var, end_var, domain)

      if (is_reader) then
         if (query_read_copy(name_handle, my_copy)) then

            allocate(var_block(block_size))
            call read_variables(ncfile, var_block, start_var, end_var, domain)

         endif
      endif

      start_rank = get_start_rank(start_var, domain)

      ! loop through and post recieves
      RECEIVING_PE_LOOP: do recv_pe = recv_start, recv_end

         ! work out elm_count on the receiving pe
         elm_count = num_elements_on_pe(recv_pe, start_rank, block_size)
         iend = istart + elm_count -1

         ! work out the start in var_block corresponding to the receiving pe
         start_point = find_start_point(recv_pe, start_rank)

         if (my_pe == recv_pe) then ! get ready to recieve from each reader

            ensemble_member =  copies_read + 1

            RECEIVE_FROM_EACH: do sending_pe = send_start, send_end

               if (query_read_copy(name_handle, sending_pe + copies_read + 1)) then

                  if(sending_pe == recv_pe) then ! just copy
                     state_ens_handle%copies(ensemble_member, istart:iend ) = &
                     var_block(start_point:elm_count*task_count():task_count())
                  else ! post receive
                     call wait_to_receive(state_ens_handle, sending_pe, ensemble_member, istart, iend)
                  endif

               endif

               ensemble_member = ensemble_member + 1

            enddo RECEIVE_FROM_EACH

            ! update starting point

            istart = istart + elm_count

         elseif (is_reader) then ! sending

            if (query_read_copy(name_handle, my_copy)) then

               call send_to_waiting_task(state_ens_handle, recv_pe, start_point, &
                                         elm_count, block_size, var_block)

            endif

         endif

      enddo RECEIVING_PE_LOOP

      start_var = end_var + 1

      if (is_reader) then
         if (query_read_copy(name_handle, my_copy)) deallocate(var_block)
      endif

   enddo VARIABLE_LOOP

   ! keep track of how many copies have been read.
   copies_read = copies_read + task_count()

   ! close netcdf file
   if (is_reader) then
      if (query_read_copy(name_handle, my_copy)) then
         ret = nf90_close(ncfile)
         call nc_check(ret, 'read_transpose closing', netcdf_filename)
      endif
   endif

enddo COPIES

dart_index = istart

end subroutine read_transpose_multi_task


!-------------------------------------------------
!> Transpose from state_ens_handle%copies to the writers according to
!> the memory limit imposed by write_var_by_var.
!>
!> This is assuming round-robin layout of state on procesors (distribution type 1
!> in the ensemble handle).

subroutine transpose_write_multi_task(state_ens_handle, name_handle, domain, &
                dart_index, write_var_by_var, write_single_precision)

type(ensemble_type),       intent(inout) :: state_ens_handle
type(stage_metadata_type), intent(in)    :: name_handle
integer,                   intent(in)    :: domain
integer(i8),               intent(inout) :: dart_index
logical,                   intent(in)    :: write_var_by_var !< Write a single variable, one at a time
logical,                   intent(in)    :: write_single_precision

integer(i8) :: i
integer :: start_var, end_var !< start/end variables in a read block
integer :: my_pe !< task or pe?
integer :: recv_pe, sending_pe
real(r8), allocatable :: var_block(:) !< for reading in variables
integer(i8) :: block_size !< number of variables in a block
integer(i8) :: elm_count !< number of elements to send
integer(i8) :: istart!< position in state_ens_handle%copies
integer(i8) :: iend
integer :: ens_size !< ensemble size
integer :: start_rank
integer :: recv_start, recv_end
integer :: send_start, send_end
integer :: ensemble_member
integer :: dummy_loop
integer :: my_copy !< which copy a pe is reading, starting from 1 to num_copies
integer :: c !< copies_read loop index
integer :: copies_written

integer :: ncfile_out !< netcdf output file handle
character(len=256) :: netcdf_filename_out !< different for each task

integer :: num_state_variables

! single file
type(time_type) :: dart_time
integer :: time_owner, time_owner_index

logical :: is_writer, clamp_vars, force_copy

ens_size = state_ens_handle%num_copies ! have the extras incase you want to read inflation restarts
my_pe = state_ens_handle%my_pe
num_state_variables = get_num_variables(domain)

! need to calculate RECEIVING_PE_LOOP start:end, group size, sending_pe start:end for each group.
! Flipped send and recv compared to read_transpose.
call get_pe_loops(ens_size, send_start, send_end, recv_start, recv_end)
if (my_pe < ens_size) then  ! I am a writer
   is_writer = .true.
else
   is_writer = .false.
endif

copies_written = 0

COPIES : do c = 1, ens_size
   if (copies_written >= ens_size) exit

   start_var = 1 ! collect first variable first
   istart = dart_index ! position in state_ens_handle%copies

   my_copy = copies_written + my_pe + 1

   ! writers open netcdf output file. This is a copy of the input file
   if (is_writer) then
      if ( query_write_copy(name_handle, my_copy)) then
         netcdf_filename_out = get_restart_filename(name_handle, my_copy, domain)

         if(file_exist(netcdf_filename_out)) then
            ret = nf90_open(netcdf_filename_out, NF90_WRITE, ncfile_out)
            call nc_check(ret, 'transpose_write: opening', trim(netcdf_filename_out))
            call nc_write_global_att_clamping(ncfile_out, my_copy, domain)

            if (overwrite_time_in_output_file) then
               call get_copy_owner_index(state_ens_handle, my_copy, time_owner, time_owner_index)
               call get_ensemble_time(state_ens_handle, time_owner_index, dart_time)
               call write_model_time(ncfile_out, dart_time)
            endif

         else ! create and open output file

            !>@todo This is grabbing the time assuming the ensemble is var complete.
            !> Should we instead have all copies time in the ensemble handle?
            call get_copy_owner_index(state_ens_handle, my_copy, time_owner, time_owner_index)
            call get_ensemble_time(state_ens_handle, time_owner_index, dart_time)

            ncfile_out = create_and_open_state_output(name_handle, domain, my_copy, &
                            dart_time, write_single_precision)
         endif
      endif

   endif

   VARIABLE_LOOP: do dummy_loop = 1, num_state_variables
      if (start_var > num_state_variables) exit ! instead of using do while loop

      ! calculate how many variables will be sent to writer
      if (write_var_by_var) then 
         end_var = start_var
      else
         end_var = num_state_variables
      endif

      block_size = get_sum_variables(start_var, end_var, domain)

      if (is_writer) then
         if (query_write_copy(name_handle, my_copy)) then
            allocate(var_block(block_size))
         endif
      endif

      start_rank =  get_start_rank(start_var, domain)

      SENDING_PE_LOOP: do sending_pe = send_start, send_end

         ! work out elm_count on the sending pe
         elm_count = num_elements_on_pe(sending_pe, start_rank, block_size)
         iend = istart + elm_count -1

         ! work out the start in var_block corresponding to the sending_pe
         i = find_start_point(sending_pe, start_rank)

         if (my_pe /= sending_pe ) then ! post recieves

            if (query_write_copy(name_handle, my_copy)) then
               call recv_variables_to_write(state_ens_handle, sending_pe, i, elm_count, block_size, var_block)
            endif

         else ! send to the collector

            ensemble_member =  copies_written + 1

            do recv_pe = recv_start, recv_end ! no if statement because everyone sends

               if (query_write_copy(name_handle, recv_pe + copies_written + 1)) then

                  if ( recv_pe /= my_pe ) then
                     call send_variables_to_write(state_ens_handle, recv_pe, ensemble_member, istart, iend)
                  else ! if sender = receiver just copy
                     var_block(i:elm_count*task_count():task_count()) = state_ens_handle%copies(ensemble_member, istart:iend)
                  endif

               endif

               ensemble_member = ensemble_member + 1

            enddo

            ! update starting point
            istart = istart + elm_count

         endif

      enddo SENDING_PE_LOOP

      if (is_writer) then ! I am a writer
         if ( query_write_copy(name_handle, my_copy)) then
            !var_block = MISSING_R8  ! if you want to create a file for bitwise testing
            clamp_vars = copy_is_clamped(name_handle, my_copy)
            force_copy = force_copy_back(name_handle, my_copy)
            call write_variables(ncfile_out, var_block, start_var, end_var, &
                                 domain, clamp_vars, force_copy)
            deallocate(var_block)
         endif
      endif

      start_var = end_var + 1

   enddo VARIABLE_LOOP

   ! keep track of how many copies have been written
   copies_written = copies_written + task_count()

   ! close netcdf file
   if (is_writer) then
      if (query_write_copy(name_handle, my_copy)) then
         ret = nf90_close(ncfile_out)
         call nc_check(ret, 'transpose_write', 'closing')
      endif
   endif

enddo COPIES

dart_index = istart

end subroutine transpose_write_multi_task


!-------------------------------------------------------------------------------
!> Check a variable for out of bounds and clamp or fail if needed.
!> If the variable has clamping limits, this routine returns .TRUE.
!> If the variable is unbounded, this routine returns .FALSE.
!> The return value is not an indication of whether or not the values have
!> actually been modified.
!-------------------------------------------------------------------------------

subroutine clamp_variable(dom_id, var_index, variable)

integer,     intent(in) :: dom_id      ! domain id
integer,     intent(in) :: var_index   ! variable index
real(r8), intent(inout) :: variable(:) ! variable

real(r8) :: minclamp, maxclamp, my_minmax(2)
character(len=NF90_MAX_NAME) :: varname ! for informational log messages
logical  :: allow_missing ! used in CLM for state variables

! if neither bound is set, return early
minclamp = get_io_clamping_minval(dom_id, var_index)
maxclamp = get_io_clamping_maxval(dom_id, var_index)

if (minclamp == missing_r8 .and. maxclamp == missing_r8) return

! if we get here, either the min, max or both have a clamping value.
  
!>@todo this is what the code needs to be for CLM and any other
! model that allows missing values in the state.  right now that
! is defined in assim_tools_mod but i don't think we can use it
! because of circular module dependencies.  it should be defined
! maybe in filter?  and set into some low level module (like types
! or constants or options_mod so anyone can query it).
!
! if we allow missing values in the state (which jeff has never
! liked because it makes the statistics funny), then these next
! two lines need to be:
allow_missing = get_missing_ok_status()

if (allow_missing) then
   my_minmax(1) = minval(variable, mask=(variable /= missing_r8))
   my_minmax(2) = maxval(variable, mask=(variable /= missing_r8))
else
   ! get the min/max for this variable before we start
   my_minmax(1) = minval(variable)
   my_minmax(2) = maxval(variable)
endif
     
varname = get_variable_name(dom_id, var_index)

! is lower bound set?
if ( minclamp /= missing_r8 ) then ! missing_r8 is flag for no clamping
   if ( my_minmax(1) < minclamp ) then
      !>@todo again, if we're allowing missing in state, this has to be masked:
       if (allow_missing) then
          where(variable /= missing_r8) variable = max(minclamp, variable)
       else
          variable = max(minclamp, variable)
       endif
   
! TJH TOO VERBOSE      write(msgstring, *) trim(varname)// ' lower bound ', minclamp, ' min value ', my_minmax(1)
! TJH TOO VERBOSE      call error_handler(E_ALLMSG, 'clamp_variable', msgstring, &
! TJH TOO VERBOSE                         source)
   endif
endif ! min range set

! is upper bound set?
if ( maxclamp /= missing_r8 ) then ! missing_r8 is flag for no clamping
   if ( my_minmax(2) > maxclamp ) then
      !>@todo again, if we're allowing missing in state, this has to be masked:
      if (allow_missing) then
         where(variable /= missing_r8) variable = min(maxclamp, variable)
      else
         variable = min(maxclamp, variable)
      endif

! TJH TOO VERBOSE      write(msgstring, *) trim(varname)// ' upper bound ', maxclamp, ' max value ', my_minmax(2)
! TJH TOO VERBOSE      call error_handler(E_ALLMSG, 'clamp_variable', msgstring, &
! TJH TOO VERBOSE                         source)
   endif

endif ! max range set

end subroutine clamp_variable



!-------------------------------------------------------------------------------
!> Write variables from start_var to end_var no clamping
!-------------------------------------------------------------------------------

subroutine write_variables(ncid, var_block, start_var, end_var, domain, &
                           do_file_clamping, force_copy)

integer,  intent(in)    :: ncid
real(r8), intent(inout) :: var_block(:)
integer,  intent(in)    :: start_var
integer,  intent(in)    :: end_var
integer,  intent(in)    :: domain
logical,  intent(in)    :: do_file_clamping
logical,  intent(in)    :: force_copy

integer(i8) :: istart, iend, var_size
integer :: i, ret, var_id, unlim_dimID
integer :: num_dims
integer, allocatable :: counts(:)
integer, allocatable :: slice_start(:) ! slice of variable

logical :: missing_possible

missing_possible = get_missing_ok_status()

!>@todo reduce output in log file?
! clamp_variable() currently prints out a line per variable per ensemble member.
! this results in a lot of output in the log file.  we may want to enable or
! disable the clamping output with a namelist or some other mechanism.  it would
! be nice to print a single value per variable across all ensemble members, but
! at this point we've gathered the variables for a single ensemble onto different
! tasks, so only N tasks (where N = number of ensemble members) have the information.
! we'd need a selective gather or a loop of send_to() calls to get the info into
! a single task for writing.

istart = 1
do i = start_var, end_var

   var_size = get_variable_size(domain, i)
   iend = istart + var_size - 1
  
   ! Some diagnostic variables do not need to be  updated. 
   ! This information is stored in the state structure and
   ! set by the model.
   if ( do_io_update(domain, i) .or. force_copy ) then
      ! diagnostic files do not get clamped but restart may be clamped
      if ( do_io_clamping(domain, i) .and. do_file_clamping) then
         call clamp_variable(domain, i, var_block(istart:iend))
      endif
     
      ! number of dimensions and length of each variable
      num_dims = get_io_num_dims(domain, i)
      allocate(counts(num_dims))
      allocate(slice_start(num_dims))
      slice_start(:) = 1 ! default to read all dimensions starting at 1
      counts(:) = 1

      if (has_unlimited_dim(domain)) then

         counts(num_dims) = 1 ! one slice of unlimited dimesion
         counts(1:get_num_dims(domain, i)) = get_dim_lengths(domain, i)

         ! write the latest time slice - HK hack to get started with tiegcm
         ! not sure if it will always be the last time slice
         ret = nf90_inquire(ncid, unlimitedDimID=unlim_dimID)
         call nc_check(ret, 'write_variables: nf90_inquire', 'unlimitedDimID')
         if (unlim_dimID /= -1) then ! unlimited dimension exists
            ret = nf90_inquire_dimension(ncid, unlim_dimID, len=slice_start(num_dims))
            call nc_check(ret, 'write_variables: nf90_inquire dimension', 'unlimitedDim length')
            if (slice_start(num_dims) == 0) slice_start(num_dims) = 1 ! newly created file
         else  ! file does not have an unlimited dimension because it was created by DART
            slice_start(num_dims) = 1
         endif

      else

         counts(1:get_num_dims(domain, i)) = get_dim_lengths(domain, i)
      endif


!>@todo FIXME, the first variable in the second domain is not found when using coamps_nest.
      ret = nf90_inq_varid(ncid, trim(get_variable_name(domain, i)), var_id)
      call nc_check(ret, 'write_variables:', 'nf90_inq_varid "'//trim(get_variable_name(domain,i))//'"')

      if (missing_possible) call set_model_missing_value(var_block(istart:iend), domain, i)

      ret = nf90_put_var(ncid, var_id, var_block(istart:iend), count=counts, start=slice_start)
      call nc_check(ret, 'write_variables:', 'nf90_put_var "'//trim(get_variable_name(domain,i))//'"')

      deallocate(counts, slice_start)
   endif

   istart = istart + var_size

enddo

end subroutine write_variables

!-------------------------------------------------------------------------------
!> Write variables from start_var to end_var for actual ensemble members
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Create the output files
!>
!> A 'blank' domain is one variable called state, with dimension = model size.
!> It is used when the model has not supplied any netcdf info but
!>     direct_netcdf_write = .true.
!> The file is intentionally left OPEN.
!-------------------------------------------------------------------------------

function create_and_open_state_output(name_handle, dom_id, copy_number, &
                dart_time, single_precision_output) result(ncfile_out)

type(stage_metadata_type), intent(in) :: name_handle
integer,                   intent(in) :: dom_id
integer,                   intent(in) :: copy_number
type(time_type),           intent(in) :: dart_time
logical,                   intent(in) :: single_precision_output
integer :: ncfile_out

character(len=*), parameter :: routine = 'create_and_open_state_output'

integer :: ret
integer :: create_mode
integer :: i, j
integer :: new_dimid
integer :: new_varid
integer :: ndims
integer :: xtype
integer :: dimids(NF90_MAX_VAR_DIMS)

character(len=256) :: filename

filename = get_restart_filename(name_handle, copy_number, dom_id)

write(msgstring,*) 'Creating output file ', trim(filename)
call error_handler(E_ALLMSG, routine, msgstring)

! What file options do you want?
create_mode = ior(NF90_CLOBBER, NF90_64BIT_OFFSET)
ret = nf90_create(filename, create_mode, ncfile_out)
call nc_check(ret, routine, 'nf90_create "'//trim(filename)//'"')

ret = nf90_enddef(ncfile_out)
call nc_check(ret, routine, 'end define mode')

! write grid information
call nc_write_model_atts(ncfile_out, dom_id)
call nc_check(nf90_Redef(ncfile_out), routine, 'redef ')

! filename discription
call nc_write_file_information(ncfile_out, filename, &
          get_file_description(name_handle, copy_number, dom_id))

! revision information
call nc_write_revision_info(ncfile_out)

! clamping information
call nc_write_global_att_clamping(ncfile_out, copy_number, dom_id, &
          from_scratch=.true.)

! define dimensions, loop around unique dimensions
do i = 1, get_io_num_unique_dims(dom_id)
   if ( trim(get_io_unique_dim_name(dom_id, i)) == 'time' ) then
      ret = nf90_def_dim(ncfile_out, 'time', NF90_UNLIMITED, new_dimid)
   else
      ret = nf90_def_dim(ncfile_out, get_io_unique_dim_name(dom_id, i), &
                       get_io_unique_dim_length(dom_id, i), new_dimid)
   endif
   !>@todo if we already have a unique names we can take this test out
   !HK no you can not, because nc_write_model_atts can create dimemsions.
   if(ret /= NF90_NOERR .and. ret /= NF90_ENAMEINUSE) then
      call nc_check(ret, routine, &
              'defining dimensions'//trim(get_io_unique_dim_name(dom_id, i)))
   endif
enddo

! define variables
do i = 1, get_num_variables(dom_id) ! loop around state variables
   if ( do_io_update(dom_id, i) .or. &
        force_copy_back(name_handle, copy_number) ) then

      ! double or single precision?
      ndims = get_io_num_dims(dom_id, i)
  
      if (single_precision_output) then
         xtype = NF90_REAL
      else ! write output that is the precision of filter
         xtype = get_xtype(dom_id, i)
         if (r8 == r4 .and. xtype == NF90_DOUBLE) xtype = NF90_REAL
      endif
  
      ! query the dimension ids
      do j = 1, ndims
         ret = nf90_inq_dimid(ncfile_out, get_dim_name(dom_id, i, j), dimids(j))
         call nc_check(ret, routine, 'querying dimensions')
      enddo
  
      ! define variable name and attributes
      write(msgstring,*) '"'//trim(get_variable_name(dom_id, i))//'"'
      ret = nf90_def_var(ncfile_out, trim(get_variable_name(dom_id, i)), &
                         xtype=xtype, dimids=dimids(1:ndims), varid=new_varid)
      call nc_check(ret, routine, 'defining variable '//trim(msgstring))
  
      call set_var_id(dom_id, i, new_varid)
  
      call write_variable_attributes(filename, ncfile_out, new_varid, dom_id, i, &
                                     name_handle, copy_number)
   endif

enddo

ret = nf90_enddef(ncfile_out)
call nc_check(ret, routine, 'nf90_enddef end define mode')

call write_model_time(ncfile_out, dart_time)

end function create_and_open_state_output


!-------------------------------------------------
!> Write model attributes if they exist. This is only
!> used for files that are created by DART.

subroutine write_variable_attributes(filename, ncFileID, ncVarID, domid, &
                                     varid, name_handle, copy_number)

character(len=*),          intent(in) :: filename
integer,                   intent(in) :: ncFileID
integer,                   intent(in) :: ncVarID
integer,                   intent(in) :: domid
integer,                   intent(in) :: varid
type(stage_metadata_type), intent(in) :: name_handle
integer,                   intent(in) :: copy_number

! long name attribute
if ( get_long_name(domid, varid) /= ' ' ) then
  call nc_check(nf90_put_att(ncFileID,ncVarID,'long_name',get_long_name(domid, varid)),&
                'write_variable_attributes','long_name in : '//trim(filename))
endif

! short name attribute
if ( get_short_name(domid, varid) /= ' ' ) then
  call nc_check(nf90_put_att(ncFileID,ncVarID,'short_name',get_short_name(domid, varid)),&
                'write_variable_attributes','short_name in : '//trim(filename))
endif

! check to see if template file has missing value attributes
if ( get_has_missing_value(domid, varid) .or. &
     get_has_FillValue(    domid, varid) ) then
   select case ( get_xtype(domid, varid) )
      case ( NF90_INT )
         call nc_write_missing_value_int(ncFileID, filename, ncVarID, domid, varid)
      case ( NF90_FLOAT )
         call nc_write_missing_value_r4 (ncFileID, filename, ncVarID, domid, varid)
      case ( NF90_DOUBLE )
         call nc_write_missing_value_r8 (ncFileID, filename, ncVarID, domid, varid)
   end select
endif

!>@todo FIXME: also need to have different routines for different different types
!>             of numbers for add_offset and scale_factor.  Keeping it simple for
!>             now since they are not being used.
if (get_scale_factor(domid, varid) /= MISSING_R8) then
  call nc_check(nf90_put_att(ncFileID,ncVarID,'scale_factor',get_scale_factor(domid, varid)),&
                'write_variable_attributes','scale_factor '//trim(filename))
endif

if (get_add_offset(domid, varid) /= MISSING_R8) then
  call nc_check(nf90_put_att(ncFileID,ncVarID,'add_offset',get_add_offset(domid, varid)),&
                'write_variable_attributes','add_offset '//trim(filename))
endif

!! ONLY for name_handle specific attributes.

! attributes that contain a file description, including whether it was clamped or not
if ( copy_is_clamped(name_handle, copy_number) ) then
   call  nc_write_variable_att_clamping(ncFileID, filename, ncVarID, domid, varid)
   write(msgstring,'(2A)') trim(get_file_description(name_handle, copy_number, domid)), '[clamped]'
   call nc_check(nf90_put_att(ncFileID,NF90_GLOBAL,'DART_note',msgstring),&
                    'write_variable_attributes','note in : '//trim(filename))
else
   write(msgstring,'( A)') trim(get_file_description(name_handle, copy_number, domid))
   call nc_check(nf90_put_att(ncFileID,NF90_GLOBAL,'DART_note',msgstring),&
                    'write_variable_attributes','note in : '//trim(filename))
endif

!>@todo FIXME: should this be independent of the name_handle?  right now it is
!>             just used to put unitless for the sd and inflation values
!>             attributes for variables without units such as inflation and sd
if( inherit_copy_units(name_handle, copy_number) ) then
   if (get_units(domid, varid) /= ' ') then
      call nc_check(nf90_put_att(ncFileID,ncVarID,'units',get_units(domid, varid)),&
                   'write_variable_attributes','units in :'//trim(filename))
   endif
else 
   call nc_check(nf90_put_att(ncFileID,ncVarID,'units','unitless'),&
                'write_variable_attributes','units in :'//trim(filename))
endif

!>@todo FIXME: put clamping values with min_val, max_val, valid_range, or whatever the proper CF-range

end subroutine write_variable_attributes


!-------------------------------------------------
!> Write global clamping attributes to files that already exist
!>@todo FIXME ? use derived type for copy number and get long name


subroutine nc_write_file_information(ncFileID, filename, description)

integer,          intent(in) :: ncFileID
character(len=*), intent(in) :: filename
character(len=*), intent(in) :: description

call nc_check(nf90_put_att(ncFileID,NF90_GLOBAL,'DART_file_information',description),&
             'nc_write_file_information','file_information'//trim(filename))

end subroutine nc_write_file_information


!-------------------------------------------------
!> Write clamping to variable attributes to files created from scratch

subroutine nc_write_variable_att_clamping(ncFileID, filename, ncVarID, domid, varid)

integer,          intent(in) :: ncFileID
character(len=*), intent(in) :: filename
integer,          intent(in) :: ncVarID
integer,          intent(in) :: domid
integer,          intent(in) :: varid

real(r8) :: clamp_val

if ( do_io_clamping(domid, varid) ) then
   clamp_val = get_io_clamping_maxval(domid, varid)

   ! max clamping attribute
   if ( clamp_val /= MISSING_R8 ) then
      call nc_check(nf90_put_att(ncFileID,ncVarID,'DART_clamp_max',clamp_val),&
                   'nc_write_variable_att_clamping','DART_clamp_max'//trim(filename))
   endif

   ! min clamping attribute
   clamp_val = get_io_clamping_minval(domid, varid)
   if ( clamp_val /= MISSING_R8 ) then
      call nc_check(nf90_put_att(ncFileID,ncVarID,'DART_clamp_min',clamp_val),&
                   'nc_write_variable_att_clamping','DART_clamp_min'//trim(filename))
   endif

endif

end subroutine nc_write_variable_att_clamping


!-------------------------------------------------
!> Write global clamping attributes for variables that have clamping
!>@todo FIXME check to see if this has a performance impact.

subroutine nc_write_global_att_clamping(ncFileID, copy, domid, from_scratch)

integer, intent(in)  :: ncFileID
integer, intent(in)  :: copy
integer, intent(in)  :: domid
logical, intent(in), optional :: from_scratch

integer  :: ivar
real(r8) :: clamp_val
character(len=NF90_MAX_NAME) :: clamp_max, clamp_min, att_name
logical :: need_netcdf_def_mode

need_netcdf_def_mode = .true.

if (present(from_scratch)) then
   if (from_scratch) need_netcdf_def_mode = .false.
endif

if (need_netcdf_def_mode) call nc_check(nf90_Redef(ncFileID),'nc_write_global_att_clamping',   'redef ')

do ivar = 1,get_num_variables(domid)
   if ( do_io_clamping(domid, ivar) ) then

     write(clamp_min,*)  'NA'
     write(clamp_max,*)  'NA'
    
     clamp_val = get_io_clamping_maxval(domid, ivar)
     if ( clamp_val /= MISSING_R8 ) write(clamp_max,*)  clamp_val
    
     clamp_val = get_io_clamping_minval(domid, ivar)
     if ( clamp_val /= MISSING_R8 ) write(clamp_min,*)  clamp_val
    
     write(msgstring,'(''min_val = '',A15,'' , max val = '',A15)') trim(clamp_min), trim(clamp_max)
     write(att_name,'(2A)')  'DART_clamp_', trim(get_variable_name(domid, ivar))
     call nc_check(nf90_put_att(ncFileID,NF90_GLOBAL,att_name, msgstring), &
                'nc_write_global_att_clamping','DART_clamping_range')
   endif
enddo

if (need_netcdf_def_mode) call nc_check(nf90_enddef(ncFileID), 'nc_write_global_att_clamping', 'end define mode')

end subroutine nc_write_global_att_clamping


!-------------------------------------------------
!> Write revision information
!>@todo this should only be done for _new_ files that DART creates - performance issue
!>@todo this will change when we move to GIT


subroutine nc_write_revision_info(ncFileID)

integer,          intent(in) :: ncFileID

! we are going to need these to record the creation date in the netCDF file.
! This is entirely optional, but nice.

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'DART_creation_date', str1), &
           'nc_write_revision_info', 'creation put ')
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'DART_source', source), &
           'nc_write_revision_info', 'source put ')

end subroutine nc_write_revision_info


!-------------------------------------------------
!> Write model integer missing_value/_FillValue attributes if they exist


subroutine nc_write_missing_value_int(ncFileID, filename, ncVarID, domid, varid)

integer,          intent(in) :: ncFileID
character(len=*), intent(in) :: filename
integer,          intent(in) :: ncVarID
integer,          intent(in) :: domid
integer,          intent(in) :: varid

integer :: missingValINT, spvalINT

call get_missing_value(domid, varid, missingValINT)
if (missingValINT /= MISSING_I) then
   call nc_check(nf90_put_att(ncFileID,ncVarID,'missing_value',missingValINT), &
                 'nc_write_missing_value_int','missing_value '//trim(filename))
endif
call get_fillValue(domid, varid, spvalINT)
if (spvalINT /= MISSING_I) then
   call nc_check(nf90_put_att(ncFileID,ncVarID,'_FillValue',spvalINT), &
                 'nc_write_missing_value_int','_FillValue'//trim(filename))
endif

end subroutine nc_write_missing_value_int


!-------------------------------------------------
!> Write model r4 missing_value/_FillValue attributes if they exist


subroutine nc_write_missing_value_r4(ncFileID, filename, ncVarID, domid, varid)

integer,          intent(in) :: ncFileID
character(len=*), intent(in) :: filename
integer,          intent(in) :: ncVarID
integer,          intent(in) :: domid
integer,          intent(in) :: varid

real(r4) :: missingValR4, spvalR4

call get_missing_value(domid, varid, missingValR4)
if (missingValR4 /= MISSING_R4) then
   call nc_check(nf90_put_att(ncFileID,ncVarID,'missing_value',missingValR4), &
                 'nc_write_missing_value_r4','missing_value '//trim(filename))
endif
call get_fillValue(domid, varid, spvalR4)
if (spValR4 /= MISSING_R4) then
   call nc_check(nf90_put_att(ncFileID,ncVarID,'_FillValue',spvalR4), &
                 'nc_write_missing_value_r4','_FillValue'//trim(filename))
endif

end subroutine nc_write_missing_value_r4


!-------------------------------------------------
!> Write model r8 missing_value/_FillValue attributes if they exist


subroutine nc_write_missing_value_r8(ncFileID, filename, ncVarID, domid, varid)

integer,          intent(in) :: ncFileID
character(len=*), intent(in) :: filename
integer,          intent(in) :: ncVarID
integer,          intent(in) :: domid
integer,          intent(in) :: varid

real(digits12) :: missingValR8, spvalR8

call get_missing_value(domid, varid, missingValR8)
if (missingValR8 /= MISSING_R8) then
   call nc_check(nf90_put_att(ncFileID,ncVarID,'missing_value',missingValR8), &
                 'nc_write_missing_value_r8','missing_value '//trim(filename))
endif
call get_fillValue(domid, varid, spvalR8)
if (spvalR8 /= MISSING_R8) then
   call nc_check(nf90_put_att(ncFileID,ncVarID,'_FillValue',spvalR8), &
                 'nc_write_missing_value_r8','_FillValue'//trim(filename))
endif

end subroutine nc_write_missing_value_r8


!-------------------------------------------------
!> Find pes for loop indices
!> if you have M tasks and N ensemble members, only the first
!> N + extra copy tasks will send, everyone will receive.


subroutine get_pe_loops(ens_size, recv_start, recv_end, send_start, send_end)

integer, intent(in)    :: ens_size
integer, intent(out)   :: recv_start !< for RECIEVING_PE_LOOP
integer, intent(out)   :: recv_end   !< for RECIEVING_PE_LOOP
integer, intent(out)   :: send_start !< for RECEIVE_FROM_EACH_LOOP
integer, intent(out)   :: send_end   !< for RECEIVE_FROM_EACH_LOOP

recv_start = 0
recv_end   = task_count() -1

send_start  = recv_start

if (ens_size > task_count()) then
   send_end = send_start + task_count() -1
else
   send_end = send_start + ens_size -1
endif

end subroutine get_pe_loops


!--------------------------------------------------------
!--------------------------------------------------------
! Routines that are making the assumption that the ensemble
! distribution is round-robin (distribution type 1)
!------------------------------------------------------
!--------------------------------------------------------
!> Send elements of variables to a processor. This routine must be called
!> with a corresponding 'wait_to_receive'.
!> The data on the sender are non-contiguous with a stride of task_count. The
!> start is different depending on which pe is the recv and which variables
!> are being sent (these are calculated in the calling routine).
!> The data on the receiver are contiguous (it is the %copies array)


subroutine send_to_waiting_task(state_ens_handle, recv_pe, start, elm_count, block_size, variable_block)

type(ensemble_type), intent(in) :: state_ens_handle
integer,             intent(in) :: recv_pe ! receiving pe
integer(i8),         intent(in) :: start ! start in variable block on sender.
integer(i8),         intent(in) :: elm_count ! how many elements
integer(i8),         intent(in) :: block_size ! size of info on sender - the receiver only
                                              ! gets part of this.
real(r8),            intent(in) :: variable_block(block_size) ! variable info

real(r8), allocatable :: buffer(:) ! for making send array contiguous

if (state_ens_handle%distribution_type == 1) then

   !>@todo MPI vector data type or packing could/should be used here?
   
   ! need a contiguous array to send variables with MPI
   allocate(buffer(elm_count))
   buffer = variable_block(start:elm_count*task_count():task_count())
   call send_to(map_pe_to_task(state_ens_handle, recv_pe), buffer)
   deallocate(buffer)

else

   call error_handler(E_ERR, 'send_to_waiting_task', 'distributions other than 1 not supported')

endif

end subroutine send_to_waiting_task


!--------------------------------------------------------
!> Send the data from a pe to a writer.
!> Note this may be 1 variable or many.  Start is the start index in %copies
!> on the sending pe. Finish is the last index in %copies to send.
!> If all variables are transposed at once,
!> start = 1,
!> finish = ens_handle%my_num_vars  (on sending pe)
!> This routine must be called with a corresponding 'recv_variables_to_write.'
!> The data on the sender is non-contiguous since it is a ROW of %copies.


subroutine send_variables_to_write(state_ens_handle, recv_pe, &
                ensemble_member, start, finish)

type(ensemble_type), intent(in) :: state_ens_handle
integer,             intent(in) :: recv_pe ! receiving pe
integer,             intent(in) :: ensemble_member
integer(i8),         intent(in) :: start  ! start in copies array on sender.
integer(i8),         intent(in) :: finish ! end in copies array on sender

real(r8), allocatable :: buffer(:) ! for making send array contiguous

if (state_ens_handle%distribution_type == 1) then

   ! MPI vector data type or packing should be used here.
   allocate(buffer(finish - start + 1))
   buffer = state_ens_handle%copies(ensemble_member, start:finish)
   call send_to(map_pe_to_task(state_ens_handle, recv_pe), buffer)
   deallocate(buffer)

else

   call error_handler(E_ERR, 'send_variables_to_write', 'distributions other than 1 not supported')

endif

end subroutine send_variables_to_write


!--------------------------------------------------------
!> Receive data from a reader. Start and finish are the local indicies
!> in the %copies array for the data being received.
!> If all variables are transposed at once,
!> start = 1,
!> finish = ens_handle%my_num_vars  (on receiveing pe)
!> This routine must be called with a corresponding 'send_to_waiting_task'.
!> The data on the sender is non-contiguous since it is a ROW of %copies.


subroutine wait_to_receive(state_ens_handle, recv_pe, &
                ensemble_member, start, finish)

type(ensemble_type), intent(inout) :: state_ens_handle
integer,             intent(in)    :: recv_pe ! receiving pe
integer,             intent(in)    :: ensemble_member
integer(i8),         intent(in)    :: start  ! start in copies array on sender.
integer(i8),         intent(in)    :: finish ! end in copies array on sender

real(r8), allocatable :: buffer(:) ! for making send array contiguous

if (state_ens_handle%distribution_type == 1) then

   ! MPI vector data type or packing should be used here.
   allocate(buffer(finish - start + 1))
   call receive_from(map_pe_to_task(state_ens_handle, recv_pe), buffer)
   state_ens_handle%copies(ensemble_member, start:finish) = buffer
   deallocate(buffer)

else

   call error_handler(E_ERR, 'wait_to_receive', 'distributions other than 1 not supported')

endif

end subroutine wait_to_receive


!--------------------------------------------------------
!> Receive data to write/collect. The data is put non-contiguously into
!> variable_block.  Variable_block is the block of data writen to a
!> netcdf file. It may be 1 or more variables.
!> start and elm_count depend on the sending pe. These are calculated in the
!> calling code. The stride is task_count() since we are assuming a round robin
!> distribution of state vector onto processors.
!> This routine must be called with a corresponding 'send_variables_to_write'


subroutine recv_variables_to_write(state_ens_handle, sending_pe, start, &
                elm_count, block_size, variable_block)

type(ensemble_type), intent(in)    :: state_ens_handle
integer,             intent(in)    :: sending_pe !! sending_pe
integer(i8),         intent(in)    :: start      !! start in vars array on receiver.
integer(i8),         intent(in)    :: elm_count  !! how many elements
integer(i8),         intent(in)    :: block_size !! size of info on sender - the receiver only gets part of this.
real(r8),            intent(inout) :: variable_block(block_size) !! variable info

real(r8), allocatable :: buffer(:) !! for making send array contiguous

if (state_ens_handle%distribution_type == 1) then

   ! MPI vector data type or packing should be used here.
   allocate(buffer(elm_count))
   call receive_from(map_pe_to_task(state_ens_handle, sending_pe), &
                           buffer)
   variable_block(start:elm_count*task_count():task_count()) = buffer
   deallocate(buffer)

else

   call error_handler(E_ERR, 'recv_variables_to_write', 'distributions other than 1 not supported')

endif

end subroutine recv_variables_to_write


!-------------------------------------------------------------------------------
!> Write state to the last time slice of a file when using "single file io".

subroutine write_model_variables(ncFileID, model_state, memindex, curr_ens_time, timeindex)

type(netcdf_file_type), intent(inout) :: ncFileID
real(r8),               intent(in) :: model_state(:)
integer,                intent(in) :: memindex
type(time_type),        intent(in) :: curr_ens_time
integer,                intent(in) :: timeindex

integer, dimension(NF90_MAX_VAR_DIMS) :: dim_lengths
integer, dimension(NF90_MAX_VAR_DIMS) :: start_point
integer(i8) :: istart, iend
integer :: ivar, jdim
integer :: ndims
integer :: ret ! netcdf return code
integer :: var_id ! netcdf variable id
integer :: domain
integer :: my_ncid
character(len=NF90_MAX_NAME) :: dimname, varname

! may not be needed
ncFileID%diag_id = create_diagnostic_structure()
my_ncid          = ncFileID%ncid

!>@todo ONLY ONE DOMAIN FOR SINGLE FILE OUTPUT
domain = 1

do ivar = 1, get_num_variables(domain)

   istart  = get_index_start(  domain, ivar)
   iend    = get_index_end(    domain, ivar)
   ndims   = get_io_num_dims(  domain, ivar)
   varname = get_variable_name(domain, ivar)

   dim_lengths(:) = 1
   start_point(:) = 1

   do jdim = 1, ndims
      dimname = get_dim_name(domain, ivar, jdim)
      if (dimname == 'time') then
         start_point(jdim) = timeindex
      else if (dimname == 'member') then
         start_point(jdim) = memindex
      else
         dim_lengths(jdim) = get_dim_length(domain, ivar, jdim)
      endif
   enddo

   ret = nf90_inq_varid(my_ncid, varname, var_id)
   call nc_check(ret, 'write_model_variables', 'inq_varid '//trim(varname))

   ret = nf90_put_var(my_ncid, var_id, model_state(istart:iend), &
                count=dim_lengths(1:ndims), start=start_point(1:ndims))
   call nc_check(ret, 'write_model_variables', 'put_var '//trim(varname))

enddo

end subroutine write_model_variables


!-------------------------------------------------------------------------------
!> Write model attributes to for each variable when using "single file io" 

subroutine write_model_attributes(ncFileID, member_dimID, time_dimId, file_handle)

type(netcdf_file_type), intent(inout) :: ncFileID
integer,                intent(in)    :: member_dimID
integer,                intent(in)    :: time_dimId
type(file_info_type),   intent(in)    :: file_handle

integer :: domain, my_ncid ! local variable
integer :: ret ! netcdf return code
integer :: copynumber, dummy
integer :: dimids(NF90_MAX_VAR_DIMS)
integer :: ndims
integer :: my_xtype ! precision for netcdf variable
integer :: ivar, jdim ! loop variables
integer :: new_varid
character(len=256) :: fname
character(len=NF90_MAX_NAME) :: dimname, varname

call nc_check(nf90_Redef(ncFileID%ncid), 'write_model_attributes', 'nf90_Redef')

ncFileID%diag_id = create_diagnostic_structure()

my_ncid    = ncFileID%ncid
fname      = ncFileID%fname
copynumber = 1 ! assuming it is an ensemble memeber

if(my_task_id()==SINGLE_IO_TASK_ID) then

   !>@todo we should always write double, unless the single output namelist
   !>flag is set.  this gives you no way to write single precision if you've
   !>compiled with normal r8 - FIXME
   if (r8 == digits12) then
      my_xtype = nf90_double
   else
      my_xtype = nf90_real
   endif

   !>@todo ONLY ONE DOMAIN FOR SINGLE FILES
   domain = 1

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
      call write_variable_attributes(fname, my_ncid, new_varid, domain, ivar, &
                                     file_handle%stage_metadata, copynumber)

   enddo

! Leave define mode
ret = nf90_enddef(my_ncid)
call nc_check(ret, 'write_model_attributes', 'enddef')

endif

end subroutine write_model_attributes


!-----------------------------------------------------------
!> Write out extra copies (if specified in the filter_nml) in the
!> last time slot.  This is only used for the when writting single
!> file diagnostc/restart files

subroutine write_extra_variables(ncFileID, model_state, copyname, curr_ens_time, timeindex)

type(netcdf_file_type), intent(inout) :: ncFileID
real(r8),               intent(in)    :: model_state(:)
character(len=*),       intent(in)    :: copyname
type(time_type),        intent(in)    :: curr_ens_time
integer,                intent(in)    :: timeindex

integer, dimension(NF90_MAX_VAR_DIMS) :: dim_lengths
integer, dimension(NF90_MAX_VAR_DIMS) :: start_point
integer :: istart, iend
integer :: ivar, jdim
integer :: ndims
integer :: ret ! netcdf return code
integer :: var_id ! netcdf variable id
integer :: domain, dcount
integer :: my_ncid
character(len=NF90_MAX_NAME) :: dimname, varname, extraname

my_ncid = ncFileID%ncid

!>@todo ONLY ONE DOMAIN FOR SINGLE FILE OUTPUT
domain    = 1

do ivar = 1, get_num_variables(domain)

   istart   = get_index_start(   domain, ivar )
   iend     = get_index_end(     domain, ivar )
   ndims    = get_io_num_dims(   domain, ivar )
   varname  = get_variable_name( domain, ivar )

   ! set the defaults and then change any needed below
   dim_lengths(:) = 1
   start_point(:) = 1
   dcount = 0

   do jdim = 1, ndims
      dimname = get_dim_name(domain, ivar, jdim)
      if (dimname == 'time') then
         dcount = dcount + 1
         start_point(dcount) = timeindex
      else if (dimname == 'member') then
         continue   ! extra vars have no member dimension
      else 
         dcount = dcount + 1
         dim_lengths(dcount) = get_dim_length(domain, ivar, jdim)
      endif
   enddo

   write(extraname,'(a,"_",a)') trim(varname), trim(copyname)
   ret = nf90_inq_varid(my_ncid, extraname, var_id)
   call nc_check(ret, 'write_extra_variables', 'inq_varid '//trim(extraname))

   ret = nf90_put_var(my_ncid, var_id, model_state(istart:iend), &
                count=dim_lengths(1:dcount), start=start_point(1:dcount))
   call nc_check(ret, 'write_extra_variables', 'put_var '//trim(extraname))

enddo

end subroutine write_extra_variables

!-----------------------------------------------------------
!> Write out extra copies (if specified in the filter_nml) attribute in the
!> last time slot.  This is only  used when writting single
!> file diagnostc/restart files.

subroutine write_extra_attributes(ncFileID, time_dimId, file_handle, copy_index)

type(netcdf_file_type), intent(inout) :: ncFileID
integer,                intent(in)    :: time_dimId
type(file_info_type),   intent(in)    :: file_handle
integer,                intent(in)    :: copy_index

integer :: ivar, jdim, ndims, domain
integer :: ret, my_ncid, new_varid, my_dimid, my_xtype
character(len=NF90_MAX_NAME) :: varname, dimname, extraname
character(len=256)           :: fname, copyname
integer ::  model_dimids(NF90_MAX_VAR_DIMS)

!--------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file 
!--------------------------------------------------------------------

call nc_check(nf90_Redef(ncFileID%ncid), 'write_extra_attributes', 'nf90_Redef')

my_ncid  = ncFileID%ncid
fname    = ncFileID%fname
copyname = get_copy_name(file_handle,copy_index)

if(my_task_id()==SINGLE_IO_TASK_ID) then

   !>@todo we should always write double, unless the single output namelist
   !>flag is set.  this gives you no way to write single precision if you've
   !>compiled with normal r8 - FIXME
   if (r8 == digits12) then
      my_xtype = nf90_double
   else
      my_xtype = nf90_real
   endif

   !>@todo ONLY ONE DOMAIN FOR SINGLE FILES
   domain = 1

   ! Define dimensions for state
   do ivar = 1, get_num_variables(domain)
      varname = get_variable_name(domain, ivar)
      ndims   = get_num_dims(domain, ivar)
   
      do jdim = 1, ndims
         dimname = get_dim_name(domain, ivar, jdim)
         ret = nf90_inq_dimid(my_ncid, dimname, dimid=my_dimid)
         call nc_check(ret, 'write_extra_attributes', 'inq_dimid '//trim(dimname))
         model_dimids(jdim) = my_dimid
      enddo
   
      model_dimids(ndims+1) = time_dimId
   
      write(extraname,'(a,"_",a)') trim(varname), trim(copyname)
      ret = nf90_inq_varid(my_ncid, extraname, new_varid)
      if(ret /= NF90_NOERR .and. ret /= NF90_ENAMEINUSE) then
         ret = nf90_def_var(my_ncid, name   = extraname,& 
                                  xtype  = my_xtype, &
                                  dimids = model_dimids(1:ndims+1), &
                                  varid  = new_varid)
         call nc_check(ret, 'write_extra_attributes', 'defining variable '//trim(extraname))
         call write_variable_attributes(fname, my_ncid, new_varid, domain, ivar, &
                                        file_handle%stage_metadata, copy_index)
      endif
   enddo
endif

! Leave define mode so we can fill
call nc_check(nf90_enddef(ncfileID%ncid), 'write_extra_attributes', 'nf90_enddef')

end subroutine write_extra_attributes


!------------------------------------------------------------------
!> helper function to return the netcdf dimension id's and lengths

subroutine get_dimension_info(copy, dom_id, var_id, timestep, is_extra, numdims, start, lengths)

integer, intent(in)  :: copy 
integer, intent(in)  :: dom_id 
integer, intent(in)  :: var_id 
integer, intent(in)  :: timestep
logical, intent(in)  :: is_extra
integer, intent(out) :: numdims
integer, intent(out) :: start  (NF90_MAX_VAR_DIMS)
integer, intent(out) :: lengths(NF90_MAX_VAR_DIMS)

logical :: is_member
integer :: member_id, dcount, jdim
character(len=NF90_MAX_NAME) :: dimname

! set defaults for start, count, then modify in the loop below
start  (:) = 1
lengths(:) = 1
dcount     = 0
member_id  = -1
is_member  = .false.

numdims = get_io_num_dims(dom_id, var_id)
do jdim = 1, numdims
   dimname = get_dim_name(dom_id, var_id, jdim)

   if ( dimname == 'time' ) then

      ! we agreed the expected behavior was to read the last time in file
      if ( is_extra ) then 
         ! extra copies have time but it might not be the jdim dimension
         dcount        = dcount + 1
         start(dcount) = timestep
      else
         start(jdim)   = timestep 
      endif

   else if ( dimname == 'member') then

      if( is_extra ) cycle ! extra copies do not have a member dimension

      member_id = jdim
      is_member = .true.

   else

      if ( is_extra ) then
         ! extra copies spatial dimensions might not be the same as jdim
         dcount = dcount + 1

         lengths(dcount) = get_dim_length(dom_id, var_id, jdim)
         start  (dcount) = 1 
      else
         lengths(jdim) = get_dim_length(dom_id, var_id, jdim)
      endif

   endif
   
   if ( is_extra ) numdims = dcount

enddo

! if the copy is an ensemble member we need set the netcdf variable start
! dimension.
if ( is_member ) start(member_id) = copy  ! member number

end subroutine get_dimension_info

!------------------------------------------------------------------
!> we're in read_single_file - there better be multiple members or
!> we can't read in an ensemble from this file. the only exception is
!> if we are reading in a single member to perturb into an ensemble.
!> otherwise it's an error.
   
subroutine check_singlefile_member_info(ncid, fname, ens_size, do_pert)
integer,                      intent(in) :: ncid
character(len=NF90_MAX_NAME), intent(in) :: fname
integer,                      intent(in) :: ens_size
logical,                      intent(in) :: do_pert

integer :: ret, MemDimID, member_size

ret = nf90_inq_dimid(ncid, "member", MemDimID)
if( ret == 0 ) then ! has member id

   ret = nf90_inquire_dimension(ncid, MemDimID, len=member_size) 
   call nc_check(ret, 'check_singlefile_member_info', 'inq_varid member : '//trim(fname))

   ! are there enough members to start from this file?  if you're perturbing a single member
   ! to generate an ensemble it's ok to have 1.  otherwise you have to have at least
   ! 'ens_size' members (more is ok).
   if (do_pert) then
      if (member_size < 1) then
         write(msgstring,  *) 'input file contains ', member_size, ' ensemble members; '
         write(msgstring2, *) 'requires at least 1 member to perturb'
         call error_handler(E_ERR, 'check_singlefile_member_info: ', msgstring, source)
      endif

   else if (member_size < ens_size) then
      write(msgstring,  *) 'input file only contains ', member_size, ' ensemble members; '
      write(msgstring2, *) 'requested ensemble size is ', ens_size
      call error_handler(E_ERR, 'check_singlefile_member_info: ', msgstring, source)

   endif

else
   ! if you don't have a member dimension, it's only ok if you're reading in a single
   ! array and not an ensemble.  
   if ( .not. do_pert ) then
      write(msgstring,  *) 'input file does not contain a "member" dimension; '
      write(msgstring2, *) 'cannot read in an ensemble of values'
      call error_handler(E_ERR, 'check_singlefile_member_info: ', msgstring, source)
   endif

endif

end subroutine check_singlefile_member_info


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
type(time_type),        intent(in) :: statetime
integer                            :: timeindex

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
   call error_handler(E_WARN,'nc_get_tindex',msgstring,source)
   timeindex = timeindex -   1
endif
if ( dimids(1) /= unlimitedDimID ) then
   write(msgstring,*)'"time" must be the unlimited dimension'
   call error_handler(E_WARN,'nc_get_tindex',msgstring,source)
   timeindex = timeindex -  10
endif
if ( timeindex < -1 ) then
   write(msgstring,*)'trouble deep ... can go no farther. Stopping.'
   call error_handler(E_ERR,'nc_get_tindex',msgstring,source)
endif

! convert statetime to time base of "days since ..."
call get_time(statetime, secs, days)

if (ncFileID%Ntimes < 1) then          ! First attempt at writing a state ...

   write(msgstring,*)'current unlimited  dimension length',nTlen, &
                     'for ncFileID ',trim(ncFileID%fname)
   call error_handler(E_DBG,'nc_get_tindex',msgstring,source)
   write(msgstring,*)'current time array dimension length',ncFileID%Ntimes
   call error_handler(E_DBG,'nc_get_tindex',msgstring,source)

   nTlen = nc_append_time(ncFileID, statetime)

   write(msgstring,*)'Initial time array dimension length',ncFileID%Ntimes
   call error_handler(E_DBG,'nc_get_tindex',msgstring,source)

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
              'Model time precedes earliest netCDF time.', source)

      write(msgstring,*)'          model time (days, seconds) ',days,secs
      call error_handler(E_DBG,'nc_get_tindex',msgstring,source)

      call get_time(ncFileID%times(1),secs,days)
      write(msgstring,*)'earliest netCDF time (days, seconds) ',days,secs
      call error_handler(E_DBG,'nc_get_tindex',msgstring,source)

      call error_handler(E_ERR,'nc_get_tindex', &
              'Model time precedes earliest netCDF time.', source)
      timeindex = -2

   else if ( statetime < ncFileID%times(ncFileID%Ntimes) ) then  

      ! It is somewhere in the middle without actually matching an existing time.
      ! This is very bad.

      write(msgstring,*)'model time does not match any netCDF time.'
      call error_handler(E_DBG,'nc_get_tindex',msgstring,source)
      write(msgstring,*)'model time (days, seconds) is ',days,secs
      call error_handler(E_DBG,'nc_get_tindex',msgstring,source)

      BadLoop : do i = 1,ncFileId%Ntimes   ! just find times to print before exiting

         if ( ncFileId%times(i) > statetime ) then
            call get_time(ncFileID%times(i-1),secs,days)
            write(msgstring,*)'preceding netCDF time (days, seconds) ',days,secs
            call error_handler(E_DBG,'nc_get_tindex',msgstring,source)

            call get_time(ncFileID%times(i),secs,days)
            write(msgstring,*)'subsequent netCDF time (days, seconds) ',days,secs
            call error_handler(E_ERR,'nc_get_tindex',msgstring,source)
            timeindex = -3
            exit BadLoop
         endif

      enddo BadLoop

   else ! we must need to append ... 

      timeindex = nc_append_time(ncFileID, statetime)

      write(msgstring,'(''appending model time (d,s) ('',i8,i5,'') as index '',i6, '' in ncFileID '',i10)') &
          days,secs,timeindex,my_ncid
      call error_handler(E_DBG,'nc_get_tindex',msgstring,source)

   endif
   
endif

end function nc_get_tindex

!-------------------------------------------------------------------------------
! Helper routine for nc_get_tindex
!-------------------------------------------------------------------------------
!>
!> The current time is appended to the "time" coordinate variable.
!> The new length of the "time" variable is returned.
!>
!> This REQUIRES that "time" is a coordinate variable AND it is the
!> unlimited dimension. If not ... bad things happen.

function nc_append_time(ncFileID, dart_time) result(lngth)

type(netcdf_file_type), intent(inout) :: ncFileID
type(time_type),        intent(in)    :: dart_time
integer                               :: lngth

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
           '"time" expected to be rank-1',source)

if ( dimids(1) /= unlimitedDimID ) call error_handler(E_ERR,'nc_append_time', &
           'unlimited dimension expected to be slowest-moving',source)

! make sure the mirror and the netcdf file are in sync
call nc_check(NF90_Inquire_Dimension(my_ncid, unlimitedDimID, varname, lngth ), &
           'nc_append_time', 'inquire_dimension unlimited')

if (lngth /= ncFileId%Ntimes) then
   write(msgstring,*)'netCDF file has length ',lngth,' /= mirror has length of ',ncFileId%Ntimes
   call error_handler(E_ERR,'nc_append_time', &
           'time mirror and netcdf file time dimension out-of-sync', &
           source,text2=msgstring)
endif

! make sure the time mirror can handle another entry.
if ( lngth == ncFileID%NtimesMAX ) then   

   write(msgstring,*)'doubling mirror length of ',lngth,' of ',ncFileID%fname
   call error_handler(E_DBG,'nc_append_time',msgstring,source)

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

call get_time(dart_time, secs, days)    ! get time components to append
realtime = days + secs/86400.0_digits12 ! time base is "days since ..."
lngth           = lngth + 1             ! index of new time 
ncFileID%Ntimes = lngth                 ! new working length of time mirror

call nc_check(nf90_put_var(my_ncid, TimeVarID, realtime, start=(/ lngth /) ), &
           'nc_append_time', 'put_var time')

ncFileID%times( lngth) = dart_time
ncFileID%rtimes(lngth) = realtime

write(msgstring,*)'ncFileID (',my_ncid,') : "',trim(varname), &
         '" (should be "time") has length ',lngth, ' appending t= ',realtime
call error_handler(E_DBG,'nc_append_time',msgstring,source)

end function nc_append_time

!-------------------------------------------------------------------------------
!> routine write calendar attributes
!>@todo FIXME this duplicates code in dart_time_io_mod.f90 - it should
!>be one place or the other.

function nc_write_calendar_atts(ncFileID, TimeVarID) result(ierr)

type(netcdf_file_type), intent(in) :: ncFileID
integer,                intent(in) :: TimeVarID
integer                            :: ierr

integer :: my_ncid

ierr = 0

my_ncid = ncFileID%ncid

call nc_check(nf90_put_att(my_ncid, TimeVarID, "long_name", "time"), &
              'nc_write_calendar_atts', 'put_att long_name '//trim(ncFileID%fname))

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
   call nc_check(nf90_put_att(my_ncid, TimeVarID, "calendar", "none" ), &
              'nc_write_calendar_atts', 'put_att calendar '//trim(ncFileID%fname))

   ! ncview (actually, probably udunits2) crashes or errors out or
   ! displays misleading plot axes if you use 'days since ...' as the units.
   ! if you simply use 'days' it works much better.

   call nc_check(nf90_put_att(my_ncid, TimeVarID, "units", "days"), &
               'nc_write_calendar_atts', 'put_att units '//trim(ncFileID%fname))
end select

end function nc_write_calendar_atts

!--------------------------------------------------------------------

function nc_get_num_times(fname) result(time_size) 

character(len=256), intent(in) :: fname 
integer :: time_size

integer :: my_ncid, TimeDimID, ret

time_size = -1 ! assume bad size

ret = nf90_open(fname, NF90_NOWRITE, my_ncid)
call nc_check(ret, 'nc_get_num_times: nf90_open', fname)

ret = nf90_inq_dimid(my_ncid, "time", TimeDimID)
if (ret == NF90_NOERR) then
   ret = nf90_inquire_dimension(my_ncid, TimeDimID, len=time_size) 
   call nc_check(ret, 'nc_get_num_times', 'inquire_dimension time '//trim(fname))
endif

ret = nf90_close(my_ncid)
call nc_check(ret, 'nc_get_num_times: nf90_close', fname)

end function nc_get_num_times


!=================================================
! Routines for distributing data round-robbin
!=================================================

!--------------------------------------------------------
!> Calculate number of elements going to the receiving pe (read_transpose)
!> Or being sent from the sending pe (transpose_write) for a given
!> start_rank and block_size.
!> block_size is the number of elements in a block of variables. There may
!> be 1 variable or all variables depending on buffer_state_io.
!> start_rank is the pe that owns the 1st element of the 1st variable in
!> the variable_block.
!>@todo FIXME ? This should go in ensemble manager.

function num_elements_on_pe(pe, start_rank, block_size) result(elm_count)

integer, intent(in) :: pe
integer, intent(in) :: start_rank
integer(i8), intent(in) :: block_size

integer(i8) :: elm_count, remainder

elm_count = block_size/task_count()
remainder = mod(block_size, task_count())

! mop up leftovers CHECK THESE.
! How many elements a pe gets depends on its rank relative to start_rank.
if ( (start_rank <= pe) .and. (pe) < (start_rank + remainder)) elm_count = elm_count + 1
if ( pe < (start_rank + remainder - task_count() )) elm_count = elm_count + 1

end function num_elements_on_pe


!--------------------------------------------------------
!> Give the rank of the processor that owns the start of a variable


function get_start_rank(variable, domain)

integer, intent(in) :: variable
integer, intent(in) :: domain

integer :: get_start_rank

get_start_rank = mod(get_sum_variables_below(variable, domain), task_count())

end function get_start_rank


!-------------------------------------------------------
!> Find i, the start point in var_block for a given recv_pe
!> This is assuming round robin. - will have to query the
!> ensemble manager to find this for different disrtibutions

function find_start_point(recv_pe, start_rank)

integer, intent(in)  :: recv_pe !< the receiver
integer, intent(in)  :: start_rank !< the pe that owns the 1st element of the var_block
integer              :: find_start_point

if (start_rank < recv_pe) then
   find_start_point = recv_pe - start_rank + 1
elseif(start_rank > recv_pe) then
   find_start_point = recv_pe + task_count() - start_rank + 1
else ! recv_pe = start_rank
   find_start_point = 1
endif

end function find_start_point


!--------------------------------------------------------
!> replace the netCDF missing_value or _FillValue with
!> the DART missing value.

subroutine set_dart_missing_value(array, domain, variable)

real(r8), intent(inout) :: array(:)
integer,  intent(in)    :: domain
integer,  intent(in)    :: variable

integer        :: model_missing_valueINT
real(r4)       :: model_missing_valueR4
real(digits12) :: model_missing_valueR8

! check to see if variable has missing value attributes
if ( get_has_missing_value(domain, variable) ) then

   select case ( get_xtype(domain, variable) )
      case ( NF90_INT )
         call get_missing_value(domain, variable, model_missing_valueINT)
         where(array == model_missing_valueINT) array = MISSING_R8
      case ( NF90_FLOAT )
         call get_missing_value(domain, variable, model_missing_valueR4)
         where(array == model_missing_valueR4) array = MISSING_R8
      case ( NF90_DOUBLE )
         call get_missing_value(domain, variable, model_missing_valueR8)
         where(array == model_missing_valueR8) array = MISSING_R8
   end select

endif

if ( get_has_FillValue(domain, variable) ) then

   select case ( get_xtype(domain, variable) )
      case ( NF90_INT )
         call get_FillValue(domain, variable, model_missing_valueINT)
         where(array == model_missing_valueINT) array = MISSING_R8
      case ( NF90_FLOAT )
         call get_FillValue(domain, variable, model_missing_valueR4)
         where(array == model_missing_valueR4) array = MISSING_R8
      case ( NF90_DOUBLE )
         call get_FillValue(domain, variable, model_missing_valueR8)
         where(array == model_missing_valueR8) array = MISSING_R8
   end select

endif

end subroutine set_dart_missing_value

!--------------------------------------------------------
!> replace the DART missing value code with the 
!> original netCDF missing_value (or _FillValue) value.

subroutine set_model_missing_value(array, domain, variable)

real(r8), intent(inout) :: array(:)
integer,  intent(in)    :: domain
integer,  intent(in)    :: variable

integer        :: model_missing_valueINT
real(r4)       :: model_missing_valueR4
real(digits12) :: model_missing_valueR8

! check to see if variable has missing value attributes
if ( get_has_missing_value(domain, variable) ) then

   select case ( get_xtype(domain, variable) )
      case ( NF90_INT )
         call get_missing_value(domain, variable, model_missing_valueINT)
         where(array == MISSING_R8) array = model_missing_valueINT
      case ( NF90_FLOAT )
         call get_missing_value(domain, variable, model_missing_valueR4)
         where(array == MISSING_R8) array = model_missing_valueR4
      case ( NF90_DOUBLE )
         call get_missing_value(domain, variable, model_missing_valueR8)
         where(array == MISSING_R8) array = model_missing_valueR8
   end select

endif

if ( get_has_FillValue(domain, variable) ) then

   select case ( get_xtype(domain, variable) )
      case ( NF90_INT )
         call get_FillValue(domain, variable, model_missing_valueINT)
         where(array == MISSING_R8)   array = model_missing_valueINT
      case ( NF90_FLOAT )
         call get_FillValue(domain, variable, model_missing_valueR4)
         where(array == MISSING_R8)   array = model_missing_valueR4
      case ( NF90_DOUBLE )
         call get_FillValue(domain, variable, model_missing_valueR8)
         where(array == MISSING_R8)   array = model_missing_valueR8
   end select

endif

end subroutine set_model_missing_value

!--------------------------------------------------------
!--------------------------------------------------------

!> @}
end module direct_netcdf_mod

