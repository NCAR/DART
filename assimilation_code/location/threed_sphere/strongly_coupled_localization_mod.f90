! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> Implements strongly coupled localization for threed_sphere locations

module strongly_coupled_localization_mod

use      types_mod, only : r8, i8
use  utilities_mod, only : error_handler, E_ERR, nmlfileunit, find_namelist_in_file, &
                           check_namelist_read, do_nml_file, do_nml_term
use   location_mod, only : get_close_type, location_type, get_location, write_location
use   ensemble_manager_mod, only : ensemble_type

implicit none
private

public get_close_state_strongly_coupled

character(len=*), parameter :: source = 'threed_sphere/strongly_coupled_localization_mod.f90'

logical, save         :: module_initialized = .false.
integer, save         :: s_model = -1, o_model = -1

! Define an integer string for the supported Earth system components
integer, parameter :: ATMOSPHERE = 1, LAND = 2

character(len = 512) :: errstring, msgstring, msgstring1, msgstring2

!-----------------------------------------------------------------
! Namelist with default values
! strongly_coupled -> logical, default false, 
!     true if strongly coupled localization should be done
! state_model -> string specifying Earth system component for the assimilating model
! obs_model   -> string specifying Earth system component of model that generated obs
! Currently supported values for state_model and obs_model are Atmosphere and Land

! Turn on strongly coupled get_close computation
logical :: strongly_coupled = .false.
! state_model and obs_model can be used if strongly_coupled is true
! Currently supported options are 'Atmosphere' and 'Land'
! Model that is currently doing the DA
character(len=256) :: state_model = 'none'
! Model that generated the observations
character(len=256) :: obs_model = 'none'

namelist /strongly_coupled_localization_nml/ strongly_coupled, state_model, obs_model


contains

!-----------------------------------------------------------------
subroutine initialize_module()

integer :: iunit, io, i, k, typecount, type_index

if (module_initialized) return

module_initialized = .true.

! Read the namelist
call find_namelist_in_file("input.nml", "strongly_coupled_localization_nml", iunit)
read(iunit, nml = strongly_coupled_localization_nml, iostat = io)
call check_namelist_read(iunit, io, "strongly_coupled_localization_nml")

! Write the namelist values to the log file
if(do_nml_file()) write(nmlfileunit, nml=strongly_coupled_localization_nml)
if(do_nml_term()) write(     *     , nml=strongly_coupled_localization_nml)

! Only cases currently supported are Atmosphere and Land
! Specify one for obs and one for state     

if(index('Atmosphere', trim(state_model)) > 0) then
   s_model = ATMOSPHERE
else if(index('Lorenz-63', trim(state_model)) > 0) then
   s_model = LAND
else
   write(errstring,*) 'state_model in strongly_coupled_localization_nml must be Atmosphere or LAND'
   call error_handler(E_ERR, 'get_close_state_strongly_coupled', errstring, source)
endif

if(index('Atmosphere', trim(obs_model)) > 0) then
   o_model = ATMOSPHERE
else if(index('Land', trim(obs_model)) > 0) then
   o_model = LAND
else
   write(errstring,*) 'obs_model in strongly_coupled_localization_nml must be Atmosphere or Land'
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

! Template for applying localization for strongly coupled 
real(r8):: obs_loc(3)
real(r8):: state_loc(3), maxdist
integer:: i, bt

if(.not. module_initialized) call initialize_module

write(*, *) 'in get_close_state_strongly_coupled', strongly_coupled

! If strongly coupled is false, just return
if(.not. strongly_coupled) return

! Figure out the appropriate base_type classification to get the maximum localization 
! distance for this observation
if (base_type < 0) then
   if (gc%nt > 1) then
      write(msgstring,  '(A)') 'no support for identity obs if per-obs-type cutoffs are specified'
      write(msgstring1, '(A)') 'contact dart support if you have a need for this combination'
      call error_handler(E_ERR, 'get_close_state_strongly_coupled', msgstring, source, text2=msgstring1)
   endif
   bt = 1
else
   ! map from type index to gtt index
   if (base_type < 1 .or. base_type > size(gc%type_to_cutoff_map)) then
      write(msgstring,'(A,I8)')'base_type out of range, is ', base_type
      write(msgstring1,'(A,2I8)')'must be between ', 1, size(gc%type_to_cutoff_map)
      call write_location (0, base_loc, charstring=msgstring2)
      call error_handler(E_ERR, 'get_close_state_strongly_coupled', msgstring, source, &
                         text2=msgstring1, text3=msgstring2)
   endif
   bt = gc%type_to_cutoff_map(base_type)
   if (bt < 1 .or. bt > gc%nt) then
      write(msgstring,'(A,I8)')'mapped type index out of range, is ', bt
      write(msgstring1,'(A,2I8)')'must be between ', 1, gc%nt
      write(msgstring2, '(A)')'internal error, should not happen.  Contact DART Support'
      call error_handler(E_ERR, 'get_close_state_strongly_coupled', msgstring, source, &
                         text2=msgstring1, text3=msgstring2)
   endif
endif

! Local variable for what the maxdist is in this particular case.
maxdist = gc%gtt(bt)%maxdist
write(*, *) 'maxdist ', maxdist
write(*, *) 'base_type ', bt

! Look at the obs location
obs_loc = get_location(base_loc)
write(*, *) 'obs_loc ', obs_loc

! Have looked at doing localization of CLM obs as function of CAM model level or
! normalized scale height. To get normalized scale height, in cam model_nml
! set vertical_localization_coord to SCALEHEIGHT and 
! no_normalization_of_scale_heights to .false. 
! This means the scale height of ps is 0 and values decrease as one moves up from the
! ps value.
! JLA CURRENTLY HARDCODING TO FORCE THE CONVERSION IN ASSIM_TOOLS
! What to do for missing_r8 for vertical location value?


! Loop through state variables that are horizontally close and add in some vertical localizattion
do i = 1, num_close
   state_loc = get_location(locs(close_ind(i)))
write(*, *) 'stateloc', state_loc(1), state_loc(2), state_loc(3)

   ! Add on some additional distances 
   ! Add on 20% of the maxdist just for being across the model boundary
   dist(i) = dist(i) + 0.2_r8 * maxdist

   if(.false.) then

      !------------------------------------------------------------------------------
      ! This block is an example for vertical localization as function of model level
      ! Surface variables have missing_r8 for vertical location for now
      ! No additonal cost for surface
      if(state_loc(3) >= 0.0_r8) then
         ! Add on additional distance that is function of model level
         ! Even more cheating by knowing the model has 32 levels and level 32 is at the surface
         dist(i) = dist(i) + maxdist * (32 - state_loc(3)) / 10.0_r8
      endif
      !------------------------------------------------------------------------------

   else

      !------------------------------------------------------------------------------
      ! This block is an example for vertical localization as function of 
      ! normalized scale height
      ! Surface variables have missing_r8 for vertical location for now
      ! No additonal cost for surface
      if(state_loc(3) >= 0.0_r8) then
         ! Add on additional distance that is function of normalized scale height
         dist(i) = dist(i) + maxdist * state_loc(3) / 4.0_r8
      endif
      !------------------------------------------------------------------------------
   endif

enddo

end subroutine get_close_state_strongly_coupled

!--------------------------------------------------------------------------

end module strongly_coupled_localization_mod
