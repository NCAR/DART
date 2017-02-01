! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program write_sqg_restart

!============================================================
! read in the model state from a netCDF file and write out 
! to a DART readable restart file
!============================================================

   use        types_mod, only : r8
   use    utilities_mod, only : initialize_utilities, finalize_utilities
   use time_manager_mod, only : time_type, set_time
   use  assim_model_mod, only : open_restart_write, awrite_state_restart
   use        model_mod, only : get_model_size, sqg_to_dart, &
                                model_static, get_model_static_data
   use     spectral_mod
   use          sqg_mod, only : init, xy_to_sp, sp_to_xy, write_diag

   implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

   integer :: model_unit, model_size, model_days, model_seconds

   real(r8), allocatable, dimension(:)     :: model_state
   real,     allocatable, dimension(:,:)   :: thxyB, thxyT
   real,     allocatable, dimension(:,:,:) :: theta

   character(len=129) :: infile, operation

   type(time_type)    :: model_time
   type(model_static) :: sqg_static

   call initialize_utilities('write_sqg_restart')

   ! get filename to read from user:
   print*, 'Enter netCDF filename to read [eg. sqgRestart.nc]: '
   read*, infile
   ! add base-state:
   print*, 'add / remove base-state [eg. add/remove/none]: '
   read*, operation
   ! get model time from user:
   print*, 'Enter model time in DAYS, SECONDS [eg. 0 0]: '
   read*, model_days, model_seconds

   model_unit = open_restart_write('sqgRestart')
   model_time = set_time(model_seconds, model_days)
   model_size = get_model_size()
   sqg_static = get_model_static_data()

   ! allocate space for all variables:
   allocate( model_state(model_size) )
   allocate( theta(2*kmax,2*lmax,2)  )
   allocate( thxyB(2*kmax,2*lmax) )
   allocate( thxyT(2*kmax,2*lmax) )

   ! read in the theta fields from the netCDF file:
   call init(infile,thxyB,thxyT)

   ! optionally, add/remove jet:
   if ( trim(operation) .eq. 'add' .or. trim(operation) .eq. 'remove' ) then
      call toggle_base_state(thxyB,thxyT,sqg_static,operation)
   endif

   ! convert theta into DART state-vector:
   theta(:,:,1) = thxyB
   theta(:,:,2) = thxyT
   call sqg_to_dart(theta,model_state)

   ! write to the file:
   call awrite_state_restart(model_time, model_state, model_unit)

   ! deallocate allocated space:
   deallocate( model_state )
   deallocate( theta       )
   deallocate( thxyB ) ; deallocate( thxyT )

   call finalize_utilities('write_sqg_restart')

   stop

contains

subroutine toggle_base_state(thxyB,thxyT,sqg_static,operation)

   implicit none

   real, dimension(2*kmax,2*lmax), intent(inout) :: thxyB, thxyT
   type(model_static),             intent(in)    :: sqg_static
   character(len=*),               intent(in)    :: operation

   complex,  dimension(2*kmax,2*lmax)      :: thspB, thspT
   integer                                 :: j
   
   if ( trim(operation) .eq. 'add' ) then
      print*,'ADDING base state'
   elseif( trim(operation) .eq. 'remove' ) then
      print*,'REMOVING base state'
   endif

   ! map into spectral space at the same resolution:
   call xy_to_sp(cmplx(thxyB,0.),thspB,2*kmax,2*lmax,kmax,lmax)
   call xy_to_sp(cmplx(thxyT,0.),thspT,2*kmax,2*lmax,kmax,lmax)

   if ( trim(operation) .eq. 'add' ) then
      thspB = thspB + sqg_static%thbB
      thspT = thspT + sqg_static%thbT
   elseif( trim(operation) .eq. 'remove' ) then
      thspB = thspB - sqg_static%thbB
      thspT = thspT - sqg_static%thbT
   endif

   ! map into grid-point space space at the same resolution:
   call sp_to_xy(thspB,thxyB,kmax,lmax,2*kmax,2*lmax)
   call sp_to_xy(thspT,thxyT,kmax,lmax,2*kmax,2*lmax)

   do j = 1, 2*lmax
      if ( trim(operation) .eq. 'add' ) then
         thxyB(:,j) = thxyB(:,j) - sqg_static%lam * real(j-1) * YL/real(2*lmax) 
         thxyT(:,j) = thxyT(:,j) - sqg_static%lam * real(j-1) * YL/real(2*lmax)
      elseif( trim(operation) .eq. 'remove' ) then
         thxyB(:,j) = thxyB(:,j) + sqg_static%lam * real(j-1) * YL/real(2*lmax) 
         thxyT(:,j) = thxyT(:,j) + sqg_static%lam * real(j-1) * YL/real(2*lmax)
      endif
   enddo

   return
end subroutine toggle_base_state

end program write_sqg_restart

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
