! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

! BEGIN DART PREPROCESS KIND LIST
! END DART PREPROCESS KIND LIST

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
! use column_calculation_mod, only :simulate_column_ob, &
!                                   vert_interp_weights
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS SIMULATE_COLUMN_OB
!         call simulate_column_ob(nlayers_obs, nlayers_model, avgkernel_obs, &
!                                  prsi_obs, prsi_model, profile_model, hofx, islog)
! END DART PREPROCESS SIMULATE_COLUMN_OB

! BEGIN DART PREPROCESS VERT_INTERP_WEIGHTS
!         call vert_interp_weights(nlev, obl, vec, wi, wf)
! END DART PREPROCESS VERT_INTERP_WEIGHTS

! BEGIN DART PREPROCESS MODULE CODE
module column_calculation_mod

    ! Generic module for column calculations
    ! This module provides the subroutine to simulate column observations
    ! and the subroutine to calculate vertical interpolation weights.
    ! It is used in the obs_def_profile_mopitt module to simulate MOPITT CO column observations.
    ! It is also used in other modules that require similar calculations TropOMI, TEMPO, etc.
    ! This module is designed to be independent of the specific observation type,
    ! allowing for reuse in different contexts where column calculations are needed.
    ! This module is orignially form JEDI UFO and is now part of the DART system.
    ! Jerome Barre, 2025-07-28

    use typeSizes
    use utilities_mod, only : error_handler, E_MSG
    implicit none
    private
    
    public :: simulate_column_ob, vert_interp_weights
    
    contains

    subroutine simulate_column_ob(nlayers_obs, nlayers_model, avgkernel_obs, &
        prsi_obs, prsi_model, profile_model, hofx, islog)

        integer, intent(in   ) :: nlayers_obs, nlayers_model
        real(r8), intent(in   ), dimension(nlayers_obs) :: avgkernel_obs
        real(r8), intent(in   ), dimension(nlayers_obs+1) :: prsi_obs
        real(r8), intent(in   ), dimension(nlayers_model+1) :: prsi_model
        real(r8), intent(in   ), dimension(nlayers_model) :: profile_model
        real(r8), intent(  out) :: hofx
        real(r8) :: wf_a, wf_b, avgkernel, grav, M_dryair
        real(r8), dimension(nlayers_obs) :: profile_obslayers
        real(r8), dimension(nlayers_obs+1) :: pobs
        real(r8), dimension(nlayers_model+1) :: pmod
        integer, parameter :: max_string=800
        character(len=max_string) :: err_msg
        integer :: k, j, wi_a, wi_b
        logical :: islog
        character(len=max_string) :: msgstring

        hofx = 0.0_r8
        profile_obslayers = 0.0_r8
        grav = 9.80665_r8
        M_dryair = 0.0289645_r8 ! kg/mol, dry air molar mass
        avgkernel = 1.0_r8
        do k=1,nlayers_obs
            ! get obs layer bound model indexes and weights for staggered
            ! obs and geoval levels
            call vert_interp_weights(nlayers_model+1, pobs(k), pmod, wi_a, wf_a)
            call vert_interp_weights(nlayers_model+1, pobs(k+1), pmod, wi_b, wf_b)

            !check if pmod is monotonic and decreasing
            if ((pmod(wi_a+1) < pmod(wi_a)) .or. (pmod(wi_b+1) < pmod(wi_b))) then
                write(msgstring, *) "Error: inverted pressure coordinate in geovals, &
                &convention: top->bottom, decreasing pressures"
                call error_handler(E_MSG,'simulate_column_ob',msgstring,source,revision,revdate)
            end if

            ! when multiple mopdel levels are in a obs layer
            if ( wi_a < wi_b ) then
                profile_obslayers(k) = profile_obslayers(k) + profile_model(wi_a) * &
                (pmod(wi_a+1)-pmod(wi_a)) * wf_a / (M_dryair*grav)
                do j=wi_a+1,wi_b-1
                    profile_obslayers(k) = profile_obslayers(k) + profile_model(j) * &
                    (pmod(j+1)-pmod(j)) / (M_dryair*grav)
                enddo
                profile_obslayers(k) = profile_obslayers(k) + profile_model(wi_b) * &
                (pmod(wi_b+1)-pmod(wi_b)) * (1.0_r8-wf_b) / (M_dryair*grav)

            ! when multiple obs layers are in a model level
            else if ( wi_a == wi_b ) then
                profile_obslayers(k) = profile_obslayers(k) + profile_model(wi_a) * &
                (pmod(wi_a+1)-pmod(wi_a)) * (wf_a-wf_b) / (M_dryair*grav)

            ! if pressures coordinates are inverted return exception
            else if ( wi_a > wi_b ) then
                write(msgstring, *) "Error: inverted pressure coordinate in obs, &
                &convention: top->bottom, decreasing pressures"
                call error_handler(E_MSG,'simulate_column_ob',msgstring,source,revision,revdate)
            end if

            if (islog) then
                hofx = hofx + (avgkernel * log10(profile_obslayers(k)))
            else if (.not. islog) then
                ! if not log, then just multiply by the avg kernel
                hofx = hofx + (avgkernel * profile_obslayers(k))
            else
                write(msgstring, *) "Error: islog must be .true. or .false."
                call error_handler(E_MSG,'simulate_column_ob',msgstring,source,revision,revdate)
            end if
        end do

        if (islog) then
            hofx = 10.0**hofx
        endif

    end subroutine simulate_column_ob

    subroutine vert_interp_weights(nlev,obl,vec,wi,wf)

        implicit none
        integer,         intent(in ) :: nlev       !Number of model levels
        real(r8), intent(in ) :: obl        !Observation location
        real(r8), intent(in ) :: vec(nlev)  !Structured vector of grid points
        integer,         intent(out) :: wi         !Index for interpolation
        real(r8), intent(out) :: wf         !Weight for interpolation

        integer         :: k

        if (vec(1) < vec(nlev)) then !Pressure increases with index
            if (obl < vec(1)) then
                wi = 1
                wf = 1.0
            elseif (obl > vec(nlev)) then
                wi = nlev - 1
                wf = 0.0
            else
                do k = 1,nlev-1
                    if (obl >= vec(k) .and. obl <= vec(k+1)) then
                        wi = k
                    endif
                enddo
                wf = (vec(wi+1) - obl)/(vec(wi+1) - vec(wi))
            endif
        else !Pressure decreases with index
            if (obl > vec(1)) then
                wi = 1
                wf = 1.0
            elseif (obl < vec(nlev)) then
                wi = nlev - 1
                wf = 0.0
            else
                do k = 1,nlev-1
                    if (obl >= vec(k+1) .and. obl <= vec(k)) then
                        wi = k
                    endif
                enddo
                wf = (vec(wi+1) - obl)/(vec(wi+1) - vec(wi))
            endif
        endif

    end subroutine vert_interp_weights

end module column_calculation_mod
! END DART PREPROCESS MODULE CODE

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$