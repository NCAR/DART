! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> Simple program that opens an obs_seq file and loops over the observations
!> and copies them to a new output file. This is intended to be a
!> template for programs that want to alter existing obs in some simple way.
!>
!> NOTE:
!> Observations in a sequence file are processed in ascending time order.
!> For efficiency at creation time, sequence files support observations in a
!> different physical order than temporal order.  This program will reorder
!> the observations to be in temporal order, so the output observations
!> may have a different observation index number than in the input file.

program cam_dart_obs_preprocessor

use        types_mod, only : r4, r8, missing_r8, metadatalength

use    utilities_mod, only : initialize_utilities,            &
                             find_namelist_in_file, check_namelist_read,       &
                             error_handler, E_ERR, E_WARN, E_MSG, nmlfileunit, &
                             do_nml_file, do_nml_term,                         &
                             open_file, close_file, finalize_utilities

use     location_mod, only : location_type, get_location, set_location,        &
                             write_location, LocationDims, operator(/=),       &
                             query_location

use      obs_def_mod, only : obs_def_type, get_obs_def_type_of_obs,            &
                             get_obs_def_location, read_obs_def,               &
                             set_obs_def_time, get_obs_def_time

use     obs_kind_mod, only : max_defined_types_of_obs,                         &
                             get_name_for_type_of_obs,                         &
                             get_index_for_type_of_obs,                        &
                             read_type_of_obs_table

use        model_mod, only : static_init_model, end_model, obs_too_high

use time_manager_mod, only : set_calendar_type

use obs_sequence_mod, only : obs_sequence_type, obs_type, write_obs_seq,       &
                             init_obs, assignment(=),                          &
                             init_obs_sequence, static_init_obs_sequence,      &
                             read_obs_seq_header, read_obs_seq, get_num_obs,   &
                             get_first_obs, get_last_obs, get_next_obs,        &
                             insert_obs_in_seq, get_num_copies, get_num_qc,    &
                             get_copy_meta_data, get_qc_meta_data,             &
                             set_copy_meta_data, set_qc_meta_data,             &
                             destroy_obs, destroy_obs_sequence,                &
                             delete_seq_head, delete_seq_tail,                 &
                             get_num_key_range, get_obs_key, get_qc,           &
                             copy_partial_obs, get_next_obs_from_key,          &
                             get_obs_def, set_obs_def,                         &
                             print_obs_seq_summary, validate_obs_seq_time

implicit none

character(len=*), parameter :: source = 'cam_dart_obs_preprocessor.f90'

type(obs_sequence_type) :: seq_in, seq_out
type(obs_type)          :: obs_in, next_obs_in
type(obs_type)          :: obs_out, prev_obs_out
type(obs_def_type)      :: this_obs_def
logical                 :: is_this_last
integer                 :: size_seq_in, size_seq_out
integer                 :: num_copies_in, num_qc_in
integer                 :: num_inserted, iunit, io, j
integer                 :: max_num_obs, file_id
character(len=129)      :: read_format
logical                 :: pre_I_format
character(len=512)      :: msgstring, msgstring1, msgstring2

character(len=metadatalength) :: meta_data


type(location_type) :: obs_loc
real(r8)            :: obs_loc_array(LocationDims)
real(r8)            :: vert_value
integer             :: which_vert 
integer             :: vert_status

!----------------------------------------------------------------
! Namelist input with default values

character(len=256)   :: filename_in  = ''
character(len=256)   :: filename_out = ''
logical              :: verbose      = .false.
character(len=32)    :: calendar     = 'Gregorian'
integer              :: print_every  = 5000

namelist /cam_dart_obs_preprocessor_nml/ &
         filename_in, filename_out, &
         verbose, calendar, print_every

!----------------------------------------------------------------
! Start of the program:
!
! Process each input observation sequence file in turn, optionally
! selecting observations to insert into the output sequence file.
!----------------------------------------------------------------

call setup()

! Read the namelist entry
call find_namelist_in_file("input.nml", "cam_dart_obs_preprocessor_nml", iunit)
read(iunit, nml = cam_dart_obs_preprocessor_nml, iostat = io)
call check_namelist_read(iunit, io, "cam_dart_obs_preprocessor_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=cam_dart_obs_preprocessor_nml)
if (do_nml_term()) write(     *     , nml=cam_dart_obs_preprocessor_nml)

! the default is a gregorian calendar.  if you are using a different type
! set it in the namelist.  this only controls how it prints out the first
! and last timestamps in the obs_seq files.
call set_calendar_type(calendar)

! if you add anything to the namelist, you can process it here.

! end of namelist processing and setup


! single pass algorithm (unlike other obs tools).

call read_obs_seq_header(filename_in, num_copies_in, num_qc_in, &
   size_seq_in, max_num_obs, file_id, read_format, pre_I_format, &
   close_the_file = .true.)

if (max_num_obs == 0) then
   write(msgstring,*) 'No obs in input sequence file ', trim(filename_in)
   call error_handler(E_ERR,'main',msgstring,source)
endif

write(msgstring, *) 'Starting to process input sequence file: '
write(msgstring1,*)  trim(filename_in)
call error_handler(E_MSG,'main',msgstring,source,text2=msgstring1)

call read_obs_seq(filename_in, 0, 0, 0, seq_in)

! sanity check - ensure the linked list times are in increasing time order
call validate_obs_seq_time(seq_in, filename_in)

! output is same size (or less) than input, generally.
! if this program is going to dup obs, account for it here.
size_seq_out = max_num_obs

! blank line, start of actually creating output file
call error_handler(E_MSG,' ',' ')

! Initialize individual observation variables
call init_obs(     obs_in,  num_copies_in, num_qc_in)
call init_obs(next_obs_in,  num_copies_in, num_qc_in)
call init_obs(     obs_out, num_copies_in, num_qc_in)
call init_obs(prev_obs_out, num_copies_in, num_qc_in)
   
! create the output sequence here 
call init_obs_sequence(seq_out, num_copies_in, num_qc_in, size_seq_out)
do j=1, num_copies_in
   meta_data = get_copy_meta_data(seq_in, j)
   call set_copy_meta_data(seq_out, j, meta_data)
enddo
do j=1, num_qc_in
   meta_data = get_qc_meta_data(seq_in, j)
   call set_qc_meta_data(seq_out, j, meta_data)
enddo

!-------------------------------------------------------------------------------
! Start to insert obs from sequence_in into sequence_out
!
! NOTE: insert_obs_in_seq CHANGES the obs passed in.
!       Must pass a copy of incoming obs to insert_obs_in_seq.
!
!-------------------------------------------------------------------------------

num_inserted = 0

if ( get_first_obs(seq_in, obs_in) )  then

   is_this_last = .false.
   next_obs_in = obs_in

   ObsLoop : do while ( .not. is_this_last )

      obs_in = next_obs_in

      ! obs_out will be modified when it is inserted in the output sequence
      ! so we have to make a copy of obs_in before modifiying it.
      obs_out = obs_in


      ! ----------- MODIFY HERE ----------------
      ! modify obs_out to match what you need

      ! If you need something in the obs_def type, here's how you do that.
      ! This just queries the location and then makes a call to unpackage
      ! the location into an array.

      ! call get_obs_def(obs_out, this_obs_def)
      ! obs_loc       = get_obs_def_location(this_obs_def)
      ! obs_loc_array = get_location(obs_loc)
 
      ! If you change something in the observation definition 
      ! you have to update the actual observation.

      ! call set_obs_def(obs_out, this_obs_def)

      ! if you do NOT want to insert this observation in the output
      ! file, set a condition that skips the next block.  you *must*
      ! still call get_next_obs() before looping.

      ! ----------- MODIFY HERE ----------------

      call get_obs_def(obs_out, this_obs_def)
      obs_loc       = get_obs_def_location(this_obs_def)
      obs_loc_array = get_location(obs_loc)

      vert_value = obs_loc_array(3)
      which_vert = int(query_location(obs_loc, "which_vert"))

      call obs_too_high(vert_value, which_vert, vert_status)
      if (verbose) print *, '   vert_value = ', vert_value, ' which_vert = ', which_vert, ' vert_status = ', vert_status

      !If vert_status is good (0), then insert the obs into the new obs_seq file
      if (vert_status == 0) then

         ! Since the stride through the observation sequence file is always
         ! guaranteed to be in temporally-ascending order, we can use the
         ! 'previous' observation as the starting point to search for the
         ! correct insertion point.  This speeds up the insert code a lot.
   
         if (num_inserted > 0) then
            call insert_obs_in_seq(seq_out, obs_out, prev_obs_out)
         else
            call insert_obs_in_seq(seq_out, obs_out)
         endif
   
         prev_obs_out = obs_out  ! update position in seq for next insert
         num_inserted = num_inserted + 1
   
         if (print_every > 0) then
            if (mod(num_inserted,print_every) == 0) then
               print*, 'inserted number ',num_inserted,' of ',size_seq_out
            endif
         endif

      endif

      ! ---- THIS MUST BE EXECUTED BEFORE THE NEXT LOOP ----

      call get_next_obs(seq_in, obs_in, next_obs_in, is_this_last)

   enddo ObsLoop

else
   write(msgstring, *)'No first observation in ',trim(filename_in)
   call error_handler(E_WARN,'main', msgstring,source)
endif

write(msgstring ,*) '---------  Obs seqs '
write(msgstring1,*) 'Number of obs input sequence    :         ', size_seq_in
write(msgstring2,*) 'Number of obs copied to output  :         ', num_inserted
call error_handler(E_MSG,'main',msgstring,source, text2=msgstring1, text3=msgstring2)
call error_handler(E_MSG,'main','---------------------------------------------------------',source)
 
if (num_inserted == 0) then
   write(msgstring ,*) 'No obs will be written to ', trim(filename_out)
   call error_handler(E_WARN,'main',msgstring,source)
else
   write(msgstring, *) 'Starting to write output sequence file ', trim(filename_out)
   call error_handler(E_MSG,'main',msgstring,source)

   call print_obs_seq_summary(seq_out, filename_out)
   call write_obs_seq(seq_out, filename_out)
endif

! clean up

call destroy_obs_sequence(seq_in)
call destroy_obs_sequence(seq_out)
call destroy_obs(     obs_in )
call destroy_obs(next_obs_in )
call destroy_obs(     obs_out)
!call destroy_obs(prev_obs_out)  ! copy of something already deleted

call shutdown()

!---------------------------------------------------------------------
! end of main program.
!---------------------------------------------------------------------


contains


!---------------------------------------------------------------------
subroutine setup()

! Initialize modules used that require it
call initialize_utilities(source)
call static_init_obs_sequence()

! read in cam-specific stuff
call static_init_model()


end subroutine setup


!---------------------------------------------------------------------
subroutine shutdown()

call end_model()
call finalize_utilities(source)

end subroutine shutdown

!---------------------------------------------------------------------

end program cam_dart_obs_preprocessor

