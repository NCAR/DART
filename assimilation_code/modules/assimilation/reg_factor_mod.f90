! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module reg_factor_mod

use     types_mod, only : r8, i8
use utilities_mod, only : get_unit, open_file, error_handler, &
                          E_ERR, nmlfileunit, find_namelist_in_file,    &
                          check_namelist_read, do_nml_file, do_nml_term

use time_manager_mod, only : time_type, get_time

implicit none
private

public :: comp_reg_factor

character(len=*), parameter :: source = 'reg_factor_mod.f90'

!============================================================================

logical :: namelist_initialized = .false.

!---- namelist with default values

integer :: select_regression = 1
! Value 1 selects default: Compute using sampling theory for any ensemble size
! Value 2 selects L96 file format: Works for archived 40 observation L96 files
! Value 3 selects bgrid archive default: Reads in file from bgrid experiments
character(len = 129) :: input_reg_file = "time_mean_reg"
character(len = 129) :: reg_diagnostics_file = "reg_diagnostics"
logical              :: save_reg_diagnostics = .false.

namelist / reg_factor_nml / select_regression, input_reg_file, &
                            save_reg_diagnostics, reg_diagnostics_file

!============================================================================


! Flags for loading startup
logical :: first_call = .true.
! Unit for output diagnostics
integer :: diag_unit
! Size of regression input files
integer :: num_obs, model_size

! Global storage for time mean regression factors from file
real(r8), allocatable :: time_mean_reg(:, :)

! Global storage for bgrid mean regression factor file
real(r8), allocatable :: obs_state_reg(:)


CONTAINS


function comp_reg_factor(num_groups, regress, obs_time, &
   obs_index, state_index, obs_state_ind, obs_state_max)

! Computes factor by which to multiply regression coefficients
! for a given distribution of sample regressions OR computes
! factor for a single sample of a regression using some other
! methodology (for instance time mean from previous runs). Could
! also implement the standard distance dependence method, too.

integer,         intent(in) :: num_groups
real(r8),        intent(in) :: regress(num_groups)
type(time_type), intent(in) :: obs_time
integer,         intent(in) :: obs_index
integer(i8),     intent(in) :: state_index
integer,         intent(in), optional :: obs_state_ind, obs_state_max

real(r8) :: comp_reg_factor

real(r8) :: sum_reg2, sum_reg_reg

integer :: i, j, ii, jj, iunit, io, secs, days

!--------------------------------------------------------
! Initialize namelist if not already done
if(.not. namelist_initialized) then


   namelist_initialized = .true.

   ! Read the namelist entry
   call find_namelist_in_file("input.nml", "reg_factor_nml", iunit)
   read(iunit, nml = reg_factor_nml, iostat = io)
   call check_namelist_read(iunit, io, "reg_factor_nml")

   ! Record the namelist values used for the run ...
   if (do_nml_file()) write(nmlfileunit, nml=reg_factor_nml)
   if (do_nml_term()) write(     *     , nml=reg_factor_nml)

   ! See if diagnostic output is requested, if so, open file
   if(save_reg_diagnostics) then
      diag_unit = open_file(reg_diagnostics_file, action = 'write')
   endif

endif
!---------------------------------------------------------

!_____________________________________________________________________
if(select_regression == 1) then

! Get regression directly from sampling theory
! If only one group, don't know what else to do
   if(num_groups == 1) then
      comp_reg_factor = 1.0_r8
   else

      sum_reg_reg = 0.0_r8
      sum_reg2 = sum(regress * regress)
      do i = 1, num_groups
         do j = i + 1, num_groups
            sum_reg_reg = sum_reg_reg + regress(i) * regress(j)
         end do                                               
      end do
      if (sum_reg2 /= 0.0_r8) then
         comp_reg_factor = 2.0_r8 * sum_reg_reg / (sum_reg2 * (num_groups - 1))
      else
         comp_reg_factor = 0.0_r8
      endif

      if(comp_reg_factor < 0.0_r8) comp_reg_factor = 0.0_r8

      ! Write out diagnostic information
      if(save_reg_diagnostics) then
       
! DATA REDUCTION FOR WORKSHOP PURPOSES
         if(obs_index <= 4 .and. state_index > 0) then

         call get_time(obs_time, secs, days)
         write(diag_unit, 22) days, secs, obs_index, state_index, comp_reg_factor
         22 format(4(i7, 1x), e14.4)
         endif
      endif

   endif

!___________________________________________________________________

else if(select_regression == 2) then

! Table lookup version for time mean, temporary implementation
! This only works for a model with a time invariant observation set
   if(first_call) then
      first_call = .false.
! Read in the regression statistics file
      iunit = get_unit()
      open(unit = iunit, file = input_reg_file)
      read(iunit, *) num_obs, model_size
      allocate(time_mean_reg(num_obs, model_size))
      do j = 1, num_obs
         do i = 1, model_size
            read(iunit, *) jj, ii, time_mean_reg(j, i)
         end do
      end do
      close(iunit)
   endif

   comp_reg_factor = time_mean_reg(obs_index, state_index)

   if(comp_reg_factor < 0.0_r8) comp_reg_factor = 0.0_r8

!_____________________________________________________________________

else if(select_regression == 3) then

   if(first_call) then
      first_call = .false.
      iunit = get_unit()
      open(unit = iunit, file = 'obs_state_reg_file')
      allocate(obs_state_reg(obs_state_max))
      do i = 1, obs_state_max
         read(iunit, 11) obs_state_reg(i)
         close(unit = iunit)
11    format(f5.3)
      end do 
   end if

   comp_reg_factor = obs_state_reg(obs_state_ind)

!_____________________________________________________________________

else
   call error_handler(E_ERR,'comp_reg_factor', &
      'Illegal value for namelist parameter select_regression',source)
endif

end function comp_reg_factor

end module reg_factor_mod

