! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program sample_grid

! Generate an identity observation based on every Nth point in the
! state vector.  The locations and kinds are queried, and then
! the observations are generated based on the kinds listed in the
! namelist and the corresponding specific type of observation
! that corresponds to that kind.  The output can be empty observations
! ready to go into perfect_model_obs, or exact observations
! based on the values in the state vector, or perturbed observations
! based on the expected error.  

! FIXME: i don't have a good way to specify the error yet, but
! i can call the ncep or ecmwf routines based on the kind and
! location to get their errors.  for now i've got a namelist
! that specifies a percentage of the obs value.

use types_mod,         only : r8

use  time_manager_mod, only : time_type, set_calendar_type, set_date, &
                              operator(>=), increment_time, get_time, &
                              operator(-), GREGORIAN, operator(+),    &
                              set_time_missing, set_time,             &
                              operator(/=), print_time, print_date
 
use utilities_mod,     only : register_module, do_output,                &
                              error_handler, nmlfileunit, E_MSG, E_ERR,  &
                              find_namelist_in_file,                     &
                              check_namelist_read, logfileunit,          &
                              do_nml_file, do_nml_term, open_file,       &
                              close_file, initialize_utilities,          &
                              finalize_utilities

use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, &
                              init_obs_sequence, get_num_obs, &
                              set_copy_meta_data, set_qc_meta_data

use  location_mod,     only : location_type

! everything here
use  obs_kind_mod

use assim_model_mod,   only : static_init_assim_model, get_model_size,   &
                              open_restart_read, open_restart_write,     &
                              awrite_state_restart, aread_state_restart, &
                              close_restart, get_state_meta_data

use obs_utilities_mod, only : add_obs_to_seq

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

integer                 :: iunit, model_size, ocount, kindindex, typeindex
integer                 :: i, io, num_kinds, state_kind, num_types
integer                 :: num_copies, num_qc, otype
integer, allocatable    :: map_list(:)
integer, parameter      :: max_list_len = 500
character(len = 128)    :: msgstring, msgstring1
real(r8), allocatable   :: state(:)
real(r8)                :: oval, err, qc
logical                 :: done, first_obs
logical, allocatable    :: usekind(:)
type(location_type)     :: loc
type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: obs_time, prev_time
type(time_type)         :: restart_time, restart_advance_time

!----------------------------------------------------------------
! These variables are namelist-controllable.
!
character(len = 128) :: input_file_name   = "filter_ics"
character(len = 128) :: output_file_name  = "obs_seq.out"
integer              :: observation_days  = -1
integer              :: observation_secs  = -1
integer              :: use_every_Nth     = 1
character(len = paramname_length) :: use_only_kinds(max_list_len) = ''
character(len = paramname_length) :: type_kind_map(2, max_list_len) = ''
logical              :: add_data          = .true.
logical              :: add_noise         = .true.
integer              :: expected_error_option   = 1
real(r8)             :: expected_error_fraction = 0.10
logical              :: debug             = .false.

! expect error options 1=percentage of obs value, 2=ncep, 3=ecmwf

! FIXME: add_data and add_noise are always true.
! add code to support no data or data with no error

namelist /sample_grid_nml/  &
   input_file_name,   &
   output_file_name,  &
   observation_days,  &
   observation_secs,  &
   use_every_Nth,     &
   use_only_kinds,    &
   type_kind_map,     &
   add_data,          &
   add_noise,         &
   expected_error_option,   &
   expected_error_fraction, &
   debug


! this could be added to the namelist if needed.
! normally the files are always restart files with a single
! timestep and not a model advance file with two times.

logical              :: input_is_model_advance_file  = .false.


!----------------------------------------------------------------
! program start
!----------------------------------------------------------------

call initialize_utilities('sample_grid')


call register_module(source,revision,revdate)


! Read the namelist entry and print it
call find_namelist_in_file("input.nml", "sample_grid_nml", iunit)
read(iunit, nml = sample_grid_nml, iostat = io)
call check_namelist_read(iunit, io, "sample_grid_nml")

if (do_nml_file()) write(nmlfileunit, nml=sample_grid_nml)
if (do_nml_term()) write(     *     , nml=sample_grid_nml)

! FIXME:
! do any error checks here
! 1) obs days/time both >= 0 or < 0
! 2) use_every_Nth >= 1

! time setup
call set_calendar_type(GREGORIAN)

! Initialize the model so we can get the size.
call static_init_assim_model()
model_size = get_model_size()

write(msgstring, *) 'Model size/restart data length =', model_size
call error_handler(E_MSG,'',msgstring)

! make space for the state data
allocate(state(model_size))

! are we generating observations from particular kinds, or the entire vector?
num_kinds = get_num_raw_obs_kinds()
allocate(usekind(num_kinds))
usekind = .true.

if (use_only_kinds(1) /= '') then
   usekind = .false.

   done = .false.
   KindList:do i=1, max_list_len
      if (use_only_kinds(i) == '') then
         done = .true.
         exit KindList
      endif
      kindindex = get_raw_obs_kind_index(use_only_kinds(i))   
      if (kindindex < 0) then
         write(msgstring, *) 'unrecognized KIND string: '//trim(use_only_kinds(i))
         call error_handler(E_ERR,'sample_grid', msgstring, &
                            source,revision,revdate)
      endif
      usekind(kindindex) = .true.

   enddo KindList

   if (.not. done) then
      write(msgstring, *) 'cannot have more than ', max_list_len, ' kinds'
      call error_handler(E_ERR,'sample_grid', msgstring, &
                         source,revision,revdate)
   endif

   write(msgstring, *) 'Generating observations based only on items in state vector items of kind:'
   call error_handler(E_MSG,'',msgstring)
   do i=1, num_kinds
      if (usekind(i)) then
         write(msgstring, *) '   ', trim(get_raw_obs_kind_name(i))
         call error_handler(E_MSG,'',msgstring)
      endif
   enddo

endif

! for every kind, we have to have a more specific obs type to
! create an observation that maps to that kind.  the list must
! be here.
if (type_kind_map(1, 1) == '') then
   write(msgstring, *) 'must specify specific types for all generic kinds being used'
   call error_handler(E_ERR,'sample_grid', msgstring, &
                            source,revision,revdate)
endif

num_types = get_num_obs_kinds()

allocate(map_list(num_kinds))

map_list(:) = -1

done = .false.
TypeList:do i=1, max_list_len
   if (type_kind_map(1, i) == '') then
      done = .true.
      exit TypeList
   endif
if (debug) print *, 'next type/kind map: ', trim(type_kind_map(1, i)), ' ', trim(type_kind_map(2, i))
   typeindex = get_obs_kind_index(type_kind_map(1, i))   
   if (typeindex < 0) then
      write(msgstring, *) 'unrecognized TYPE string: '//trim(type_kind_map(1, i))
      call error_handler(E_ERR,'sample_grid', msgstring, &
                         source,revision,revdate)
   endif
   kindindex = get_raw_obs_kind_index(type_kind_map(2, i))   
   if (kindindex < 0) then
      write(msgstring, *) 'unrecognized KIND string: '//trim(type_kind_map(2, i))
      call error_handler(E_ERR,'sample_grid', msgstring, &
                         source,revision,revdate)
   endif

   if (typeindex > num_types) then ! shouldn't happen
      write(msgstring, *) 'bad TYPE index: ', typeindex
      call error_handler(E_ERR,'sample_grid', msgstring, &
                         source,revision,revdate)
   endif
   if (kindindex > num_kinds) then ! shouldn't happen either
      write(msgstring, *) 'bad KIND index: ', kindindex
      call error_handler(E_ERR,'sample_grid', msgstring, &
                         source,revision,revdate)
   endif

   ! how to get from a generic kind to a specific type
   ! note that generating observations this way you can only
   ! create one specific observation type for any given generic
   ! state vector kind.

   map_list(kindindex) = typeindex
if (debug) print *, 'setting kindindex ', kindindex, ' to type ', typeindex
  

enddo TypeList

if (.not. done) then
   write(msgstring, *) 'cannot have more than ', max_list_len, ' kinds'
   call error_handler(E_ERR,'sample_grid', msgstring, &
                      source,revision,revdate)
endif

write(msgstring, *) 'Each generic kind in the state vector will generate a specific obs type of:'
call error_handler(E_MSG,'',msgstring)
do i=1, num_kinds
   if (usekind(i)) then
      write(msgstring, *) '   ', trim(get_raw_obs_kind_name(i)), ' -> ', trim(get_obs_kind_name(map_list(i)))
      call error_handler(E_MSG,'',msgstring)
   endif
enddo



restart_time         = set_time_missing()
restart_advance_time = set_time_missing()

! open the state vector restart file
iunit = open_restart_read(trim(input_file_name))

! Read in the advance time if present
if (input_is_model_advance_file) then
   call aread_state_restart(restart_time, state, iunit, restart_advance_time)
else
   call aread_state_restart(restart_time, state, iunit)
endif
call close_restart(iunit)


! if days/secs from namelist are both -1, then use the time
! from the restart file as the time for the observations.  
! if they are both >= 0, override the time with the namelist value.
! it's an error for one to be >= 0 and one to be < 0.
! (already checked for right after namelist read in)
if (observation_days < 0 .and. observation_secs < 0) then
   obs_time = restart_time
else
   obs_time = set_time(observation_secs, observation_days)
endif

if (add_data) then
   if (add_noise) then
      ! each obs has 2 values - truth and observation
      num_copies = 2
   else
      ! each obs has 1 value - unadjusted observation
      num_copies = 1
   endif
   num_qc     = 1
else
   ! just the location, time, kind, and expected error. no data.
   num_copies = 0
   num_qc     = 0
endif

! call the initialization code, and initialize two empty observation types
call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)
first_obs = .true.

! create a new, empty obs_seq file.
call init_obs_sequence(obs_seq, num_copies, num_qc, model_size)

! set up the metadata if needed in this case
if (add_data) then
   call set_copy_meta_data(obs_seq, 1, 'observations')
   call set_qc_meta_data(obs_seq, 1, 'Data QC')
   if (add_noise) then
      call set_copy_meta_data(obs_seq, 2, 'truth')
   endif
endif

! Set the DART data quality control.   0 is good data. 
! increasingly larger QC values are more questionable quality data.
qc = 0.0_r8

! main loop - for every item in the state vector, generate an observation
! at that location, of that type.  FIXME: do we want this to be an empty
! obs_seq.in file, or we can add gaussian here and make both obs and truth...

if (debug) print *, 'starting loop for ', model_size, ' items'
stateloop: do i=1, model_size
   call get_state_meta_data(-1 * i, loc, state_kind)
   if (state_kind < 1 .or. state_kind > num_kinds) then
      write(msgstring, *) 'bad KIND from get_state_meta_data, ', state_kind, ' for index ', i 
      write(msgstring1, *) 'must be between 1 and ', num_kinds
      call error_handler(E_ERR,'sample_grid', msgstring, &
                         source,revision,revdate, text2=msgstring1)

   endif
   if (.not. usekind(state_kind)) cycle stateloop

   ! FIXME: how do use_every and kinds interact?
   ! how do we order these cycles?

   if (modulo(i, use_every_Nth) /= 1) cycle stateloop

if (debug) print *, 'ready to create obs ', i

   oval = state(i)
   err = compute_error(state_kind, loc, oval)
   otype = map_list(state_kind)

   ! create an obs for it:
   call create_3d_obs_custom(add_data, loc, oval, otype, err, obs_time, qc, obs)
   call add_obs_to_seq(obs_seq, obs, obs_time, prev_obs, prev_time, first_obs)
   
end do stateloop

! if we added any obs to the sequence, write it out to a file now.
ocount = get_num_obs(obs_seq)
if ( ocount > 0 ) then
   if (debug) print *, 'writing obs_seq, obs_count = ', ocount
   call write_obs_seq(obs_seq, output_file_name)
endif

write(msgstring, *) 'used ', ocount, ' of ', model_size, ' items in the state vector'
call error_handler(E_MSG,'sample_grid', msgstring)

! clean up state
deallocate(state)

call finalize_utilities() 

!----------------------------------------------------------------
!----------------------------------------------------------------

contains

function compute_error(otype, oloc, oval)
 integer,             intent(in) :: otype
 type(location_type), intent(in) :: oloc
 real(r8),            intent(in) :: oval

 real(r8) :: compute_error

! FIXME:  options here could use the ecmwf or ncep error tables

! select based on expected_error_option

! for now, only option 1 works
compute_error = expected_error_fraction * oval

end function compute_error

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   create_3d_obs - subroutine that is used to create an observation
!                   type from observation data.  
!
!       NOTE: assumes the code is using the threed_sphere locations module, 
!             that the observation has a single data value and a single
!             qc value, and that this obs type has no additional required
!             data (e.g. gps and radar obs need additional data per obs)
!
!    add_data - logical to say whether there is data/qc or no
!    oloc  - threed_sphere location of observation
!    oval  - observation value
!    okind - observation kind
!    oerr  - observation error
!    otime - dart time_type of observation
!    oqc   - quality control value
!    obs   - output observation type
!
!     created Oct. 2007 Ryan Torn, NCAR/MMM
!     adapted for more generic use 11 Mar 2010, nancy collins, ncar/image
!     unadapted for specific use 3 sept 2013, nancy collins, ncar/image
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine create_3d_obs_custom(add_data, oloc, oval, okind, oerr, otime, oqc, obs)
use        types_mod, only : r8
use obs_def_mod,      only : obs_def_type, set_obs_def_time, set_obs_def_kind, &
                             set_obs_def_error_variance, set_obs_def_location
use obs_sequence_mod, only : obs_type, set_obs_values, set_qc, set_obs_def
use time_manager_mod, only : time_type, generate_seed
use     location_mod, only : location_type
use   random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian

 logical,             intent(in)    :: add_data
 type(location_type), intent(in)    :: oloc
 real(r8),            intent(in)    :: oval, oerr, oqc
 integer,             intent(in)    :: okind
 type(time_type),     intent(in)    :: otime
 type(obs_type),      intent(inout) :: obs

real(r8)           :: obs_val(2), qc_val(1)
type(obs_def_type) :: obs_def
integer            :: seed

type(random_seq_type), save  :: random_seq
logical,               save  :: first_time = .true.

if (first_time) then
   ! Initialize a repeatable random sequence for perturbations
   seed = generate_seed(otime)
   call init_random_seq(random_seq,seed)
   first_time = .false.
endif

call set_obs_def_location(obs_def, oloc)
call set_obs_def_kind(obs_def, okind)
call set_obs_def_time(obs_def, otime)
call set_obs_def_error_variance(obs_def, oerr * oerr)
call set_obs_def(obs, obs_def)

if (add_data) then
   if (add_noise) then
      obs_val(1) = random_gaussian(random_seq, oval, oerr)
      obs_val(2) = oval
      call set_obs_values(obs, obs_val)
   else
      obs_val(1) = oval
      call set_obs_values(obs, obs_val(1:1), 1)
   endif

   qc_val(1)  = oqc
   call set_qc(obs, qc_val)
endif

end subroutine create_3d_obs_custom



end program sample_grid

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
