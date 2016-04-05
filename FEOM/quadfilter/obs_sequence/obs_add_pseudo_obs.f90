! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> program that opens an obs_seq file and loops over the obs
!> and copies them to a new output file, adding a new obs type
!> QUAD_FILTER_SQUARED_ERROR after each one.  intended for use
!> only with the 'quad_filter' option (see filter.html for more info).


program obs_add_pseudo_obs

use        types_mod, only : r8, MISSING_R8, metadatalength
use    utilities_mod, only : register_module, initialize_utilities,            &
                             find_namelist_in_file, check_namelist_read,       &
                             error_handler, E_ERR, E_MSG, nmlfileunit,         &
                             do_nml_file, do_nml_term, get_next_filename,      &
                             open_file, close_file, finalize_utilities
use      obs_def_mod, only : obs_def_type, set_obs_def_error_variance,         &
                             set_obs_def_kind
use     obs_kind_mod, only : QUAD_FILTER_SQUARED_ERROR
use obs_sequence_mod, only : obs_sequence_type, obs_type, write_obs_seq,       &
                             init_obs, assignment(=), get_obs_def,             &
                             init_obs_sequence, static_init_obs_sequence,      &
                             read_obs_seq_header, read_obs_seq, get_num_obs,   &
                             get_first_obs, get_next_obs,                      &
                             insert_obs_in_seq, get_num_copies, get_num_qc,    &
                             get_copy_meta_data, get_qc_meta_data,             &
                             set_copy_meta_data, set_qc_meta_data,             &
                             destroy_obs, destroy_obs_sequence,                &
                             set_obs_def, set_obs_values

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

type(obs_sequence_type) :: seq_in, seq_out
type(obs_type)          :: obs_in, next_obs_in
type(obs_type)          :: obs_out, prev_obs_out
logical                 :: is_this_last
integer                 :: size_seq_in, size_seq_out
integer                 :: num_copies_in, num_qc_in
integer                 :: num_inserted, iunit, io, i
integer                 :: max_num_obs, file_id
character(len=129)      :: read_format
logical                 :: pre_I_format
character(len=512)      :: msgstring, msgstring1
type(obs_def_type)      :: this_obs_def

character(len=metadatalength) :: meta_data

! could go into namelist if you wanted more control of output
integer, parameter      :: print_every = 5000

integer  :: observation_copy_index
real(r8) :: value(1)


!----------------------------------------------------------------
! Namelist input with default values

character(len=256)   :: filename_in  = 'obs_seq.out'
character(len=256)   :: filename_out = 'obs_seq_with_pseudo_obs.out'


namelist /obs_add_pseudo_obs_nml/ &
         filename_in, filename_out

!----------------------------------------------------------------
! Start of the program:
!
! Process each input observation sequence file in turn, optionally
! selecting observations to insert into the output sequence file.
!----------------------------------------------------------------

call initialize()

! Read the namelist entry
call find_namelist_in_file("input.nml", "obs_add_pseudo_obs_nml", iunit)
read(iunit, nml = obs_add_pseudo_obs_nml, iostat = io)
call check_namelist_read(iunit, io, "obs_add_pseudo_obs_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=obs_add_pseudo_obs_nml)
if (do_nml_term()) write(     *     , nml=obs_add_pseudo_obs_nml)


! end of namelist processing and setup


! read in existing obs sequence and compute output values

call read_obs_seq_header(filename_in, num_copies_in, num_qc_in, &
   size_seq_in, max_num_obs, file_id, read_format, pre_I_format, &
   close_the_file = .true.)

if (max_num_obs == 0) then
   write(msgstring,*) 'No obs in input sequence file ', trim(filename_in)
   call error_handler(E_ERR,'obs_add_pseudo_obs',msgstring)
endif

write(msgstring, *) 'Starting to process input sequence file: '
write(msgstring1,*)  trim(filename_in)
call error_handler(E_MSG,'obs_add_pseudo_obs',msgstring, &
                   text2=msgstring1)

call read_obs_seq(filename_in, 0, 0, 0, seq_in)

! output is exactly twice the input size - each obs is duplicated.
size_seq_out = size_seq_in * 2

! blank line, start of actually creating output file
call error_handler(E_MSG,' ',' ')

! Initialize individual observation variables
call init_obs(     obs_in,  num_copies_in, num_qc_in)
call init_obs(next_obs_in,  num_copies_in, num_qc_in)
call init_obs(     obs_out, num_copies_in, num_qc_in)
call init_obs(prev_obs_out, num_copies_in, num_qc_in)
   
! create the output sequence here 
call init_obs_sequence(seq_out, num_copies_in, num_qc_in, size_seq_out)
do i=1, num_copies_in
   meta_data = get_copy_meta_data(seq_in, i)
   call set_copy_meta_data(seq_out, i, meta_data)
enddo
do i=1, num_qc_in
   meta_data = get_qc_meta_data(seq_in, i)
   call set_qc_meta_data(seq_out, i, meta_data)
enddo

! figure out which copy is the actual observation value
observation_copy_index = get_copy_index_for_obs(seq_in)

!-------------------------------------------------------------
! Start to insert obs from sequence_in into sequence_out
!
! NOTE: insert_obs_in_seq CHANGES the obs passed in.
!       Must pass a copy of incoming obs to insert_obs_in_seq.
!--------------------------------------------------------------

! this must be called because it initializes next_obs_in for the loop below
if (.not. get_first_obs(seq_in, next_obs_in) )  then
   write(msgstring, *)'no first observation in ',trim(filename_in)
   call error_handler(E_ERR,'obs_add_pseudo_obs', msgstring, source,revision,revdate)
endif

num_inserted = 0
is_this_last = .false.

ObsLoop : do while ( .not. is_this_last )

   obs_in = next_obs_in

   ! obs_out will be modified when it is inserted in the output sequence
   ! so we have to make a copy of obs_in before modifiying it.
   obs_out = obs_in

   ! Since the stride through the observation sequence file is always
   ! guaranteed to be in temporally-ascending order, we can use the
   ! 'previous' observation as the starting point to search for the
   ! correct insertion point.  This speeds up the insert code a lot.

   ! insert original obs, unchanged, to output 
   if (num_inserted > 0) then
      call insert_obs_in_seq(seq_out, obs_out, prev_obs_out)
   else
      call insert_obs_in_seq(seq_out, obs_out)
   endif

   prev_obs_out = obs_out  ! update position in seq for next insert
   num_inserted = num_inserted + 1

   ! ----------- QUAD OBS MODIFICATIONS HERE ----------------
   ! change obs value to MISSING_R8, err MISSING_R8
   ! change obs type to QUAD_FILTER_SQUARED_ERROR

   call get_obs_def(obs_out, this_obs_def)

   call set_obs_def_kind(this_obs_def, QUAD_FILTER_SQUARED_ERROR)
   call set_obs_def_error_variance(this_obs_def, MISSING_R8)

   call set_obs_def(obs_out, this_obs_def)

   value(1) = MISSING_R8
   call set_obs_values(obs_out, value, observation_copy_index)

   ! ----------- QUAD OBS MODIFICATIONS HERE ----------------

   ! insert the companion quad obs here
   call insert_obs_in_seq(seq_out, obs_out, prev_obs_out)

   prev_obs_out = obs_out  ! update position in seq for next insert
   num_inserted = num_inserted + 1

   if (print_every > 0) then
      if (mod(num_inserted, print_every) == 0) then
         print*, 'inserted number ',num_inserted,' of ',size_seq_out
      endif
   endif

   call get_next_obs(seq_in, obs_in, next_obs_in, is_this_last)

enddo ObsLoop

write(msgstring, *) 'Number of obs in the input  seq file :', size_seq_in
call error_handler(E_MSG, 'obs_add_pseudo_obs: ', msgstring)

write(msgstring, *) 'Number of obs in the output seq file :', size_seq_out
call error_handler(E_MSG, 'obs_add_pseudo_obs: ', msgstring)

call error_handler(E_MSG,' ',' ')

call write_obs_seq(seq_out, filename_out)

! clean up

call destroy_obs_sequence(seq_in)
call destroy_obs_sequence(seq_out)
call destroy_obs(     obs_in )
call destroy_obs(next_obs_in )
call destroy_obs(     obs_out)
!call destroy_obs(prev_obs_out)  ! copy of something already deleted

call finalize()

!---------------------------------------------------------------------
! end of main program.
!---------------------------------------------------------------------


contains
!---------------------------------------------------------------------

subroutine initialize()

! Initialize modules used that require it
call initialize_utilities('obs_add_pseudo_obs')
call register_module(source,revision,revdate)
call static_init_obs_sequence()

end subroutine initialize

!---------------------------------------------------------------------

subroutine finalize()

call finalize_utilities('obs_add_pseudo_obs')

end subroutine finalize

!---------------------------------------------------------------------

function get_copy_index_for_obs(seq)

type(obs_sequence_type), intent(in) :: seq
integer                             :: get_copy_index_for_obs

integer :: i

! Determine which data copy in sequence has actual obs value

do i = 1, get_num_copies(seq)
   get_copy_index_for_obs= i
   ! Need to look for 'observation'
   if(index(get_copy_meta_data(seq, i), 'observation') > 0) return
end do

! Not returning sooner means 'observation' not found in metadata; die
call error_handler(E_ERR,'get_copy_index_for_obs', &
   'Did not find observation copy with metadata "observation"', &
      source, revision, revdate)

end function get_copy_index_for_obs

!---------------------------------------------------------------------

end program obs_add_pseudo_obs

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
