! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> Implements strongly coupled localization for oned locations

module strongly_coupled_localization_mod

use      types_mod, only : r8, i8
use  utilities_mod, only : error_handler, E_ERR, nmlfileunit, find_namelist_in_file, &
                           check_namelist_read, do_nml_file, do_nml_term
use   location_mod, only : get_close_type, location_type, get_location, write_location
use   ensemble_manager_mod, only : ensemble_type

implicit none
private

public get_close_state_strongly_coupled

character(len=*), parameter :: source = 'oned/strongly_coupled_localization_mod.f90'

logical, save         :: module_initialized = .false.
integer, save         :: s_model = -1, o_model = -1

character(len = 512) :: errstring

!-----------------------------------------------------------------
! Namelist with default values
! strongly_coupled -> logical, default false, 
!     true if strongly coupled localization should be done
! state_model -> string specifying the low-order assimilating model
! obs_model   -> string specifying the low-order model that generated obs
! Note that these can be the same model, just two different instances

! Turn on strongly coupled get_close computation
logical :: strongly_coupled = .false.
! state_model and obs_model can be used if strongly_coupled is true
! Currently supported values are Lorenz-63 and Lorenz-96
! Model that is currently doing the DA
character(len=256) :: state_model = 'none'
! Model that generated the observations
character(len=256) :: obs_model = 'none'

namelist /strongly_coupled_localization_nml/ strongly_coupled, state_model, obs_model

contains

!-----------------------------------------------------------------
subroutine initialize_module()

integer :: iunit, io

if (module_initialized) return

module_initialized = .true.

! Read the namelist
call find_namelist_in_file("input.nml", "strongly_coupled_localization_nml", iunit)
read(iunit, nml = strongly_coupled_localization_nml, iostat = io)
call check_namelist_read(iunit, io, "strongly_coupled_localization_nml")

! Write the namelist values to the log file
if(do_nml_file()) write(nmlfileunit, nml=strongly_coupled_localization_nml)
if(do_nml_term()) write(     *     , nml=strongly_coupled_localization_nml)

! Only cases supported are Lorenz_96 with Lorenz_63
! Specify one for obs and one for state

if(index('Lorenz-96', trim(state_model)) > 0) then
   s_model = 96
else if(index('Lorenz-63', trim(state_model)) > 0) then
   s_model = 63
else
   write(errstring,*) 'state_model in strongly_coupled_localization_nml &
      must be Lorenz-96 or Lorenz-63'
   call error_handler(E_ERR, 'get_close_state_strongly_coupled', errstring, source)
endif

if(index('Lorenz-96', trim(obs_model)) > 0) then
   o_model = 96
else if(index('Lorenz-63', trim(obs_model)) > 0) then
   o_model = 63
else
   write(errstring,*) 'obs_model in strongly_coupled_localization_nml &
      must be Lorenz-96 or Lorenz-63'
   call error_handler(E_ERR, 'get_close_state_strongly_coupled', errstring, source)
endif

end subroutine initialize_module

!-----------------------------------------------------------------

subroutine get_close_state_strongly_coupled(gc, base_loc, base_type, &
   locs, loc_qtys, loc_indx, num_close, close_ind, dist, ensemble_handle)


type(get_close_type),          intent(in)  :: gc
type(location_type),           intent(in)  :: base_loc, locs(:)
integer,                       intent(in)  :: base_type, loc_qtys(:)
integer(i8),                   intent(in)  :: loc_indx(:)
integer,                       intent(inout) :: num_close, close_ind(:)
real(r8),            optional, intent(inout) :: dist(:)
type(ensemble_type), optional, intent(in)  :: ensemble_handle

! The extra distance for crossing the model boundary
real(r8), parameter :: same_extra_distance = 0.1
real(r8), parameter :: different_extra_distance = 0.2

integer :: s_model, o_model

if(.not. module_initialized) call initialize_module

! If strongly coupled is false, just return
if(.not. strongly_coupled) return

! Simple demonstration of low order model strongly coupled DA
! Make distance 'across the boundary' an extra addition
! Could update num_close and close_ind for additional efficiency

if(s_model == o_model) then
   ! If they are both the same model (but different instances) increase distance
   dist(1:num_close) = dist(1:num_close) + same_extra_distance
else
   ! If they are different models increase distance
   dist(1:num_close) = dist(1:num_close) + different_extra_distance
endif

end subroutine get_close_state_strongly_coupled

!--------------------------------------------------------------------------

end module strongly_coupled_localization_mod

