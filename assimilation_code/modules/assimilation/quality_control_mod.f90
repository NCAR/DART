! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!------------------------------------------------------------------------------
!> quality_control_mod.f90
!>
!> This module contains routines related to quality control.
!>
!------------------------------------------------------------------------------
module quality_control_mod


use     types_mod,    only : r8

use utilities_mod,    only : error_handler, E_ERR, E_MSG, &
                             do_output, do_nml_file, do_nml_term, nmlfileunit, &
                             find_namelist_in_file, check_namelist_read

use location_mod,     only : location_type

use obs_sequence_mod, only : obs_sequence_type, init_obs, get_obs_from_key, &
                             get_obs_def, obs_type

use obs_def_mod,      only : get_obs_def_type_of_obs, obs_def_type

!------------------------------------------------------------------------------

implicit none

private

public :: initialize_qc, input_qc_ok, get_dart_qc, check_outlier_threshold, &
          good_dart_qc, set_input_qc, &
          DARTQC_ASSIM_GOOD_FOP, DARTQC_EVAL_GOOD_FOP, &
          DARTQC_ASSIM_FAILED_POST_FOP, DARTQC_EVAL_FAILED_POST_FOP, &
          DARTQC_FAILED_FOP, DARTQC_NOT_IN_NAMELIST, &
          DARTQC_BAD_INCOMING_QC, DARTQC_FAILED_OUTLIER_TEST, &
          DARTQC_FAILED_VERT_CONVERT

character(len=*), parameter :: source = 'quality_control_mod.f90'

! Dart quality control variables
integer, parameter :: DARTQC_ASSIM_GOOD_FOP        = 0
integer, parameter :: DARTQC_EVAL_GOOD_FOP         = 1
integer, parameter :: DARTQC_ASSIM_FAILED_POST_FOP = 2
integer, parameter :: DARTQC_EVAL_FAILED_POST_FOP  = 3
integer, parameter :: DARTQC_FAILED_FOP            = 4
integer, parameter :: DARTQC_NOT_IN_NAMELIST       = 5
integer, parameter :: DARTQC_BAD_INCOMING_QC       = 6
integer, parameter :: DARTQC_FAILED_OUTLIER_TEST   = 7
integer, parameter :: DARTQC_FAILED_VERT_CONVERT   = 8
!integer, parameter :: DARTQC_OUTSIDE_DOMAIN        = 9  ! we have no way (yet) for the model_mod to signal this

!------------------------------------------------------------------------------
! namelist parameters
!>@todo ... in the (upcoming) html ... mention that outlier_threshold = 3.0 is a good value

real(r8) :: input_qc_threshold = 3.0_r8  ! values larger than input_qc_threshold will be rejected
real(r8) :: outlier_threshold  = -1.0_r8

logical  :: enable_special_outlier_code = .false. ! user defined outlier threshold code

namelist / quality_control_nml / input_qc_threshold, outlier_threshold,  &
   enable_special_outlier_code
!------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------
!> Subroutine initialize_qc
!>
!> read the quality control namelist
!------------------------------------------------------------------------------
subroutine initialize_qc()

character(len=512) :: msgstring

integer :: iunit, io

! Read the namelist entry
call find_namelist_in_file("input.nml", "quality_control_nml", iunit)
read(iunit, nml = quality_control_nml, iostat = io)
call check_namelist_read(iunit, io, "quality_control_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=quality_control_nml)
if (do_nml_term()) write(     *     , nml=quality_control_nml)

if(do_output()) then
   write(msgstring, '(A,I4)') 'Will reject obs with Data QC larger than ',nint(input_qc_threshold)
   call error_handler(E_MSG,'quality_control_mod:', msgstring)

   ! Let user know what settings the outlier threshold has, and data qc.
   if (outlier_threshold <= 0.0_r8) then
      write(msgstring, '(A)') 'No observation outlier threshold rejection will be done'
   else
      write(msgstring, '(A,F12.6,A)') 'Will reject obs values more than', outlier_threshold, ' sigma from mean'
   endif
   call error_handler(E_MSG,'quality_control_mod:', msgstring)

   ! if doing something special with outlier threshold, say so
   if (enable_special_outlier_code) then
      call error_handler(E_MSG,'quality_control_mod:', &
              'special outlier threshold handling enabled', source)
   endif

endif

end subroutine initialize_qc

!------------------------------------------------------------------------------
!> Function set_input_qc
!>
!> if observation is not geing evaluated or assimilated, skip it and return
!> 1000 in qc field so the observation is not used again 
!------------------------------------------------------------------------------
function set_input_qc(input_qc, assimilate_this_ob, evaluate_this_ob)

integer,  intent(in) :: input_qc   !! result from forward operator calculation 
logical,  intent(in) :: assimilate_this_ob  !! true if obs was assimilated
logical,  intent(in) :: evaluate_this_ob  !! true if obs was evaluated
integer              :: set_input_qc  !!  return value

if(input_qc == 0 .and. (assimilate_this_ob .or. evaluate_this_ob)) then
   set_input_qc = 0.0_r8
else
   set_input_qc = 1000.0_r8 + input_qc   ! give some clue about why
endif

end function set_input_qc

!------------------------------------------------------------------------------
!> Function input_qc_ok
!>
!> Do checks on input_qc value with namelist control
!>
!> NOTE: this code has changed since the original version.
!> qc values equal to the threshold are kept now; only qc
!> values LARGER THAN the threshold are rejected.  
!> e.g. to keep only obs with a data qc of 0, set the 
!> threshold to 0.
!>
!------------------------------------------------------------------------------
function input_qc_ok(input_qc, qc_to_use)

real(r8), intent(in)  :: input_qc   !! result from forward operator calculation 
integer,  intent(out) :: qc_to_use  !! resulting DART QC
logical               :: input_qc_ok !! true if input_qc is good

qc_to_use = DARTQC_ASSIM_GOOD_FOP

if(nint(input_qc) <= nint(input_qc_threshold)) then
   input_qc_ok = .true.
else
   input_qc_ok = .false.
   qc_to_use = DARTQC_BAD_INCOMING_QC 
endif

end function input_qc_ok

!------------------------------------------------------------------------------
!> Subroutine get_dart_qc
!>
!> Consolidate forward operator qc
!> find the min and max istatus values across all ensemble members.  These are
!> either set by dart code, or returned by the model-specific model_interpolate()
!> routine, or by forward operator code in obs_def_xxx_mod files.
!>
!------------------------------------------------------------------------------
subroutine get_dart_qc(istatus, ens_size, assimilate_this_ob, evaluate_this_ob, isprior, dart_qc)

integer, intent(in)    :: ens_size !!  number of ensemble members
integer, intent(in)    :: istatus(ens_size) !!   0=ok, >0 is error, <0 reserved for system use
logical, intent(in)    :: assimilate_this_ob !!   true if obs was assimilated
logical, intent(in)    :: evaluate_this_ob !!  true if obs was evaluated
logical, intent(in)    :: isprior !!  true for prior pass; false for posterior
integer, intent(inout) :: dart_qc !!   resulting dart qc

logical :: failed_fop, inconsistent

! If we didn't call the forward operators, return now.
if (dart_qc == DARTQC_BAD_INCOMING_QC) return

failed_fop   = (maxval(istatus(:)) > 0)
inconsistent = (minval(istatus(:)) /= 0)

! Now do a case statement to figure out what the qc result should be
! For prior, have to test for a bunch of stuff
! FIXME: note that this case statement doesn't cover every possibility;
! if there's an error in the code and a minus value gets into the forward
! operator istatus without being caught, it will fail all cases below.
! add another line for 'internal inconsistency' to be safe.

if(isprior) then

   ! only do this for the priors.  the posterior failures
   ! don't change whether the prior worked or not.

   if((.not. assimilate_this_ob) .and. (.not. evaluate_this_ob)) then
      dart_qc = DARTQC_NOT_IN_NAMELIST
   else if(failed_fop) then            ! At least one forward operator failed
      dart_qc = DARTQC_FAILED_FOP
   else if(evaluate_this_ob) then          ! Observation to be evaluated only
      dart_qc = DARTQC_EVAL_GOOD_FOP
   else if(.not. inconsistent) then           ! All clear, assimilate this ob
      dart_qc = DARTQC_ASSIM_GOOD_FOP
   
   else   ! 'should not happen'
      dart_qc = -99  ! inconsistent istatus codes
      call error_handler(E_ERR, 'get_dart_qc: ', &
              'internal error, should not happen: bad else in qc case', source)
   endif

else ! posterior

   ! For failed posterior, only update qc if prior was successful
   if(failed_fop) then
      ! Both the following 2 tests and assignments were on single executable lines,
      ! but one compiler (gfortran) was confused by this, so they were put in
      ! if/endif blocks.
      if(dart_qc == DARTQC_ASSIM_GOOD_FOP) then
         dart_qc = DARTQC_ASSIM_FAILED_POST_FOP

      else if(dart_qc == DARTQC_EVAL_GOOD_FOP) then
         dart_qc = DARTQC_EVAL_FAILED_POST_FOP
    
      endif
   endif

endif

end subroutine get_dart_qc

!------------------------------------------------------------------------------
!> Subroutine check_outlier_threshold
!>
!> Check on the outlier threshold quality control.  If existing DART QC
!> is good and the obs is outside the threshold, set the QC to indicate that.
!------------------------------------------------------------------------------
subroutine check_outlier_threshold(obs_prior_mean, obs_prior_var, obs_val, obs_err_var, &
                                   obs_seq, this_obs_key, dart_qc)

real(r8),                intent(in)    :: obs_prior_mean !!  prior observation mean 
real(r8),                intent(in)    :: obs_prior_var !!  prior observation variance
real(r8),                intent(in)    :: obs_val !!  observation value
real(r8),                intent(in)    :: obs_err_var !!  observation error variance
type(obs_sequence_type), intent(in)    :: obs_seq !!  observation sequence
integer,                 intent(in)    :: this_obs_key !!  index for this observation
integer,                 intent(inout) :: dart_qc !!  possibly modified DART QC

real(r8) :: error, diff_sd, ratio

logical  :: failed

! FIX ME: could be a loop if we consider groups

! Check on the outlier threshold quality control:
! BUG FIX: this was using the incoming qc threshold before.
! it should only be doing the outlier test on dart qc values of 0 or 1.
! if there is already a different qc code set, leave it alone.
! only if it is still successful (assim or eval, 0 or 1), then check
! for failing outlier test.

! If bad QC flag exit here.
if ( .not. good_dart_qc(dart_qc)) return

error   = obs_prior_mean - obs_val
diff_sd = sqrt(obs_prior_var + obs_err_var)

if (diff_sd /= 0.0_r8) then
   ratio = abs(error / diff_sd)
else
   ratio = 0.0_r8
endif

! if special handling requested, pass in the outlier ratio for this obs,
! the default outlier threshold value, and enough info to extract the specific 
! obs type for this obs. the function should return .true. if this is an 
! outlier, .false. if it is ok.
if (enable_special_outlier_code) then
   failed = failed_outlier(ratio, this_obs_key, obs_seq)
else 
   if ( outlier_threshold < 0 ) return ! don't modify dart_qc if no outlier check
   failed = (ratio > outlier_threshold)
endif

if (failed) then
   dart_qc = DARTQC_FAILED_OUTLIER_TEST
endif

end subroutine check_outlier_threshold

!------------------------------------------------------------------------------
!> Function failed_outlier
!>
!> Intended to be easy to modify if needed to use other criteria
!> when evaluating the outlier threshold.  Set the namelist to
!> turn this call on (enable_special_outlier_code) and then add your
!> own tests to this section.  Return 'true' if the tests fail,
!> return 'false' if the obs is ok.
!------------------------------------------------------------------------------
function failed_outlier(ratio, this_obs_key, seq)

! return true if the observation value is too far away from the ensemble mean
! and should be rejected and not assimilated.

real(r8),                intent(in) :: ratio !!   observation error ratio
integer,                 intent(in) :: this_obs_key !!  index for this observation
type(obs_sequence_type), intent(in) :: seq !! the observation sequence
logical                             :: failed_outlier !! true if outlier test fails

! the default test is:  if (ratio > outlier_threshold) failed_outlier = .true.
! but you can add code here to do different tests for different observation 
! types.  this function is only called if this namelist item is set to true:
!  enable_special_outlier_code = .true.
! in the &quality_control_nml namelist.  it is intended to be customized by the user.

integer :: this_obs_type

type(obs_def_type)   :: obs_def
type(obs_type), save :: observation

logical :: first_time = .true.

! make sure there's space to hold the observation.   this is a memory
! leak in that we never release this space, but we only allocate it
! one time so it doesn't grow.  if you are going to access the values
! of the observation or the qc values, you must change the 0, 0 below
! to match what's in your obs_seq file.  the example code below uses
! only things in the obs_def derived type and so doesn't need space
! allocated for either copies or qcs.
if (first_time) then
   call init_obs(observation, 0, 0)
   first_time = .false.
endif

! if you want to do something different based on the observation specific type:
call get_obs_from_key(seq, this_obs_key, observation)
call get_obs_def(observation, obs_def)
this_obs_type = get_obs_def_type_of_obs(obs_def)

! note that this example uses the specific type (e.g. RADIOSONDE_TEMPERATURE)
! to make decisions.  you have the observation so any other part (e.g. the
! time, the value, the error) is available to you as well.

select case (this_obs_type)

! example of specifying a different threshold value for one obs type:
!   case (RADIOSONDE_TEMPERATURE)
!      if (ratio > some_other_value) then
!         failed_outlier = .true.
!      else
!         failed_outlier = .false.
!      endif

! accept all values of this observation type no matter how far
! from the ensemble mean:
!   case (AIRCRAFT_U_WIND_COMPONENT, AIRCRAFT_V_WIND_COMPONENT)
!      failed_outlier = .false.

   case default
      if (ratio > outlier_threshold) then
         failed_outlier = .true.
      else
         failed_outlier = .false.
      endif
      ! outlier_threshold < 0 means don't do outlier check.
      if ( outlier_threshold < 0 ) failed_outlier = .false.

end select

end function failed_outlier

!------------------------------------------------------------------------------
!> Function good_dart_qc
!>
!> Check incoming dart QC to see if the obs was processed ok.
!>
!------------------------------------------------------------------------------
function good_dart_qc(qc_value)

integer, intent(in) :: qc_value  !! incoming DART QC (not data qc)
logical             :: good_dart_qc  !! whether this is a good value or not


if ( (qc_value == DARTQC_ASSIM_GOOD_FOP) .or. &
     (qc_value == DARTQC_EVAL_GOOD_FOP)  .or. &
     (qc_value == DARTQC_FAILED_OUTLIER_TEST)) then
   good_dart_qc = .true.
else
   good_dart_qc = .false.
endif

end function good_dart_qc

!------------------------------------------------------------------------------
end module quality_control_mod

