! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

module ice_postprocessing_mod

use types_mod, only : r8
use  utilities_mod, only : E_MSG, error_handler
use  netcdf_utilities_mod, only : nc_get_variable_size, &
                                  nc_get_variable, nc_open_file_readonly, &
                                  nc_close_file
use netcdf

implicit none
private

character(len=256), parameter :: source = 'ice_postprocessing'

! routines available from this module
public :: area_simple_squeeze, volume_simple_squeeze, &
          cice_rebalancing, read_cice_state_variable

! -----------------------------------------------------------------------------
contains
! -----------------------------------------------------------------------------

subroutine read_cice_state_variable(varname, var_array, filename)

   ! read a 3D variable from a CICE restart file. This code is based on
   ! a similar routine in the CAM-FV model module, 'read_cam_phis_array',
   ! but needed to be generalized to function for any user-supplied 3D 
   ! variable. It replaces the previous 'get_3d_variable' routine in the 
   ! 'dart_to_cice' module.

   character(len=*),      intent(in)  :: filename
   character(len=*),      intent(in)  :: varname
   character(len=*),      parameter   :: routine = 'read_cice_state_variable'
   real(r8), allocatable, intent(out) :: var_array(:,:,:)

   integer :: ncid, nsize(3)     ! ni, nj, ncat

   ncid = nc_open_file_readonly(filename, routine)

   call nc_get_variable_size(ncid, varname, nsize, routine)
   allocate(var_array(nsize(1), nsize(2), nsize(3)))

   call nc_get_variable(ncid, varname, var_array, routine)

   call nc_close_file(ncid, routine)

end subroutine read_cice_state_variable

! ----------------------------------------------------------------------------- 

 subroutine area_simple_squeeze(qice001, qice002,     &
                                qice003, qice004,     &
                                qice005, qice006,     &
                                qice007, qice008,     &
                                sice001, sice002,     &
                                sice003, sice004,     &
                                sice005, sice006,     &
                                sice007, sice008,     &
                                qsno001, qsno002,     &
                                qsno003, aicen,       &
                                vicen, vsnon,         &
                                aicen_original,       &
                                vicen_original,       &
                                vsnon_original,       &
                                Tsfcn,                &
                                Ncat, Nx, Ny)
 
    real(r8), intent(inout), dimension(Nx,Ny,Ncat) :: &
                                        aicen, vicen, vsnon,  &
                                        qice001, qice002, qice003, qice004, &
                                        qice005, qice006, qice007, qice008, &
                                        sice001, sice002, sice003, sice004, &
                                        sice005, sice006, sice007, sice008, &
                                        qsno001, qsno002, qsno003,   &
                                        Tsfcn
    real(r8), intent(in), dimension(Nx,Ny,Ncat) :: &
                                        aicen_original,    &
                                        vicen_original,    &
                                        vsnon_original  
    integer, intent(in) :: Ncat, Nx, Ny
 
    real(r8), dimension(Nx,Ny,Ncat) :: hicen_original, hsnon_original, aicen_temp
    real(r8), dimension(Nx,Ny) :: aice, aice_temp
    real(r8), dimension(Ncat) :: hcat_midpoint
    real(r8), dimension(0:Ncat) :: hin_max
    real(r8) :: squeeze, cc1, cc2, x1, Si0new, Ti, qsno_hold, qi0new
    real(r8), parameter :: Tsmelt = 0._r8,        &
                           cc3 = 3._r8,           &
                           c1 = 1._r8,            &
                           phi_init = 0.75_r8,    &
                           dSin0_frazil = 3.0_r8, &
                           sss = 34.7_r8
 
    integer :: i, j, n ! loop indices

    ! calculate dependent variables 
    cc1 = cc3/real(Ncat,kind=r8)
    cc2 = 15.0_r8*cc1
    Si0new = sss - dSin0_frazil              
     
    ! calculate bounds and midpoints of thickness distribution
    hin_max(0) = 0.0_r8
    do n = 1, Ncat
       x1 = real(n-1,kind=r8) / real(Ncat,kind=r8)
       hin_max(n) = hin_max(n-1) &
                     + cc1 + cc2*(c1 + tanh(cc3*(x1-c1)))
       hcat_midpoint(n)=0.5_r8*(hin_max(n-1)+hin_max(n))
    enddo
    
    ! Begin process 
    sice001  = max(0.0_r8, sice001)  ! salinities must be non-negative
    sice002  = max(0.0_r8, sice002)  ! salinities must be non-negative
    sice003  = max(0.0_r8, sice003)  ! salinities must be non-negative
    sice004  = max(0.0_r8, sice004)  ! salinities must be non-negative
    sice005  = max(0.0_r8, sice005)  ! salinities must be non-negative
    sice006  = max(0.0_r8, sice006)  ! salinities must be non-negative
    sice007  = max(0.0_r8, sice007)  ! salinities must be non-negative
    sice008  = max(0.0_r8, sice008)  ! salinities must be non-negative
    qice001  = min(0.0_r8, qice001)  ! enthalpies (ice) must be non-positive
    qice002  = min(0.0_r8, qice002)  ! enthalpies (ice) must be non-positive
    qice003  = min(0.0_r8, qice003)  ! enthalpies (ice) must be non-positive
    qice004  = min(0.0_r8, qice004)  ! enthalpies (ice) must be non-positive
    qice005  = min(0.0_r8, qice005)  ! enthalpies (ice) must be non-positive
    qice006  = min(0.0_r8, qice006)  ! enthalpies (ice) must be non-positive
    qice007  = min(0.0_r8, qice007)  ! enthalpies (ice) must be non-positive
    qice008  = min(0.0_r8, qice008)  ! enthalpies (ice) must be non-positive
    qsno001  = min(0.0_r8, qsno001)  ! enthalphies (snow) must be non-positive
    qsno002  = min(0.0_r8, qsno002)  ! enthalphies (snow) must be non-positive
    qsno003  = min(0.0_r8, qsno003)  ! enthalphies (snow) must be non-positive
    aicen = min(1.0_r8,aicen)  ! concentrations must not exceed 1 
    Tsfcn = min(Tsmelt,Tsfcn)  ! ice/snow surface must not exceed melting
 
    ! calculate aice, which might be negative or >1 at this point
    aice = aicen(:,:,1)
    do n = 2, Ncat  
       aice = aice+aicen(:,:,n)
    enddo
    
    ! set negative aicen to zero
    aicen = max(0.0_r8,aicen)   ! concentrations must be non-negative
    vicen = max(0.0_r8,vicen)   ! same for volumes (ice)
    vsnon = max(0.0_r8,vsnon)   ! same for volumes (snow)
 
    ! reclaculate aice, now it should be non-negative\
    aice_temp = aicen(:,:,1)
    do n = 2, Ncat
       aice_temp = aice_temp + aicen(:,:,n)
    enddo
  
    ! if aice <0, then set every category to 0
    do j = 1, Ny
       do i = 1, Nx
          if (aice(i,j)<0._r8) then
             aicen(i,j,:) = 0._r8
          endif
       enddo
    enddo
 
    ! shift negative concentration values 
    do n=1, Ncat
       do j=1, Ny
          do i=1, Nx
             if (aice_temp(i,j) > 0._r8 .and. aice(i,j)>0._r8) then
                aicen(i,j,n) = aicen(i,j,n) - (aice_temp(i,j)-aice(i,j))*aicen(i,j,n)/aice_temp(i,j)
             endif  
          enddo
       enddo
    enddo
 
    ! now squeeze aicen 
     do j = 1, Ny
       do i = 1, Nx
          if (aice(i,j) > 1.0_r8) then
             squeeze        = 1.0_r8 / aice(i,j)
             aicen(i,j,:)   = aicen(i,j,:)*squeeze
          endif
       enddo
    enddo
 
    ! update vsnon and vicen using conserved category thickness values
    aicen_temp = aicen_original
    where(aicen_temp==0) aicen_temp = -999
 
    do n=1,Ncat
       do j=1,Ny
          do i=1,Nx
             hicen_original(i,j,n) = vicen_original(i,j,n)/aicen_temp(i,j,n)
             hsnon_original(i,j,n) = vsnon_original(i,j,n)/aicen_temp(i,j,n)
          end do
       end do
    end do
 
    where(hicen_original < 0)  hicen_original = 0.0_r8
    where(hsnon_original < 0)  hsnon_original = 0.0_r8
 
    vicen  = aicen*hicen_original
    vsnon  = aicen*hsnon_original
 
   ! consider special cases
    do n = 1, Ncat
       do j = 1, Ny
          do i = 1, Nx
          ! If there is no ice post-adjustment... 
             if (aicen(i,j,n)==0._r8 ) then
                vicen(i,j,n)   = 0._r8
                vsnon(i,j,n) = 0._r8
                sice001(i,j,n) = 0._r8
                sice002(i,j,n) = 0._r8
                sice003(i,j,n) = 0._r8
                sice004(i,j,n) = 0._r8
                sice005(i,j,n) = 0._r8
                sice006(i,j,n) = 0._r8
                sice007(i,j,n) = 0._r8
                sice008(i,j,n) = 0._r8
                qice001(i,j,n) = 0._r8
                qice002(i,j,n) = 0._r8
                qice003(i,j,n) = 0._r8
                qice004(i,j,n) = 0._r8
                qice005(i,j,n) = 0._r8
                qice006(i,j,n) = 0._r8
                qice007(i,j,n) = 0._r8
                qice008(i,j,n) = 0._r8
                qsno001(i,j,n) = 0._r8
                qsno002(i,j,n) = 0._r8
                qsno003(i,j,n) = 0._r8
                Tsfcn(i,j,n)   = -1.836_r8
             ! If the adjustment introduced new ice.. 
             else if (aicen(i,j,n)>0._r8 .and. aicen_original(i,j,n)==0._r8) then
                ! allow no snow volume or enthalpy
                vsnon(i,j,n) = 0._r8
                qsno001(i,j,n) = 0._r8
                qsno002(i,j,n) = 0._r8
                qsno003(i,j,n) = 0._r8
 
                ! require ice volume for thickness = category boundary midpoint
                vicen(i,j,n) =  aicen(i,j,n) * hcat_midpoint(n)
 
                ! salinity of mushy ice, see add_new_ice in ice_therm_itd.F90
                Si0new = sss - dSin0_frazil ! given our choice of sss
                sice001(i,j,n) = Si0new
                sice002(i,j,n) = Si0new
                sice003(i,j,n) = Si0new
                sice004(i,j,n) = Si0new
                sice005(i,j,n) = Si0new
                sice006(i,j,n) = Si0new
                sice007(i,j,n) = Si0new
                sice008(i,j,n) = Si0new
 
                ! temperature and enthalpy
                Ti          = min(liquidus_temperature_mush(Si0new/phi_init), -0.1_r8)
                qi0new      = enthalpy_mush(Ti, Si0new)
                qice001(i,j,n) = qi0new
                qice002(i,j,n) = qi0new
                qice003(i,j,n) = qi0new
                qice004(i,j,n) = qi0new
                qice005(i,j,n) = qi0new
                qice006(i,j,n) = qi0new
                qice007(i,j,n) = qi0new
                qice008(i,j,n) = qi0new
                Tsfcn(i,j,n)  = Ti
             endif
          enddo
       enddo
    enddo
 
 end subroutine area_simple_squeeze
! -----------------------------------------------------------------------------
 
 subroutine volume_simple_squeeze(qice001, qice002,     &
                                  qice003, qice004,     &
                                  qice005, qice006,     &
                                  qice007, qice008,     &
                                  sice001, sice002,     &
                                  sice003, sice004,     &
                                  sice005, sice006,     &
                                  sice007, sice008,     &
                                  qsno001, qsno002,     &
                                  qsno003, aicen,       &
                                  vicen, vsnon,         &
                                  aicen_original,      &
                                  vicen_original,      &
                                  vsnon_original,      &
                                  Tsfcn,               &
                                  Ncat, Nx, Ny)
 
    real(r8), intent(inout), dimension(Nx,Ny,Ncat) :: &
                                        aicen, vicen, vsnon,  &
                          qice001, qice002, qice003, qice004, &
                          qice005, qice006, qice007, qice008, &
                          sice001, sice002, sice003, sice004, &
                          sice005, sice006, sice007, sice008, &
                                   qsno001, qsno002, qsno003, &
                                                        Tsfcn
    real(r8), intent(in), dimension(Nx,Ny,Ncat) :: &
                                           aicen_original, &
                                           vicen_original, &
                                           vsnon_original   
    integer, intent(in) :: Ncat, Nx, Ny
 
    real(r8), dimension(Nx,Ny,Ncat) :: hicen_original, hsnon_original, aicen_temp
    real(r8), dimension(Nx,Ny) :: aice, vice, vice_temp
    real(r8), dimension(0:Ncat) :: hin_max
    real(r8), dimension(Ncat) :: hcat_midpoint
    real(r8) :: squeeze, cc1, cc2, x1, Si0new, Ti, qsno_hold, qi0new
    real(r8), parameter :: Tsmelt = 0._r8,        &
                           cc3 = 3._r8,           &
                           c1 = 1._r8,            &
                           phi_init = 0.75_r8,    &
                           dSin0_frazil = 3.0_r8, &
                           sss = 34.7_r8
 
    integer :: i, j, n ! loop indices

   ! calculate dependent variables 
     cc1 = cc3/real(Ncat,kind=r8)
     cc2 = 15.0_r8*cc1
     Si0new = sss - dSin0_frazil              
 
    ! calculate bounds and midpoints of thickness distribution
    hin_max(0) = 0.0_r8
    do n = 1, Ncat
       x1 = real(n-1,kind=r8) / real(Ncat,kind=r8)
       hin_max(n) = hin_max(n-1) &
                     + cc1 + cc2*(c1 + tanh(cc3*(x1-c1)))
       hcat_midpoint(n)=0.5_r8*(hin_max(n-1)+hin_max(n))
    enddo

    ! Begin process 
    sice001  = max(0.0_r8, sice001)  ! salinities must be non-negative
    sice002  = max(0.0_r8, sice002)  ! salinities must be non-negative
    sice003  = max(0.0_r8, sice003)  ! salinities must be non-negative
    sice004  = max(0.0_r8, sice004)  ! salinities must be non-negative
    sice005  = max(0.0_r8, sice005)  ! salinities must be non-negative
    sice006  = max(0.0_r8, sice006)  ! salinities must be non-negative
    sice007  = max(0.0_r8, sice007)  ! salinities must be non-negative
    sice008  = max(0.0_r8, sice008)  ! salinities must be non-negative
    qice001  = min(0.0_r8, qice001)  ! enthalpies (ice) must be non-positive
    qice002  = min(0.0_r8, qice002)  ! enthalpies (ice) must be non-positive
    qice003  = min(0.0_r8, qice003)  ! enthalpies (ice) must be non-positive
    qice004  = min(0.0_r8, qice004)  ! enthalpies (ice) must be non-positive
    qice005  = min(0.0_r8, qice005)  ! enthalpies (ice) must be non-positive
    qice006  = min(0.0_r8, qice006)  ! enthalpies (ice) must be non-positive
    qice007  = min(0.0_r8, qice007)  ! enthalpies (ice) must be non-positive
    qice008  = min(0.0_r8, qice008)  ! enthalpies (ice) must be non-positive
    qsno001  = min(0.0_r8, qsno001)  ! enthalphies (snow) must be non-positive
    qsno002  = min(0.0_r8, qsno002)  ! enthalphies (snow) must be non-positive
    qsno003  = min(0.0_r8, qsno003)  ! enthalphies (snow) must be non-positive    ! aicen = min(1.0_r8,aicen)  ! concentrations must not exceed 1 
    Tsfcn = min(Tsmelt,Tsfcn)  ! ice/snow surface must not exceed melting
 
    ! calculate aice, which might be negative or >1 at this point
    vice = vicen(:,:,1)
    do n = 2, Ncat  
       vice = vice+vicen(:,:,n)
    enddo
 
    ! set negative aicen to zero
    vicen = max(0.0_r8,vicen)   ! same for volumes (ice)
 
    ! reclaculate aice, now it should be non-negative
    vice_temp = vicen(:,:,1)
    do n = 2, Ncat
       vice_temp = vice_temp + vicen(:,:,n)
    enddo
 
    ! if vice <0, then set every category to 0
    do j = 1, Ny
       do i = 1, Nx
          if (vice(i,j)<0._r8) then
             vicen(i,j,:) = 0._r8
          endif
       enddo
    enddo
 
    ! shift negative post-adjustment volume values 
    do n=1, Ncat
       do j=1, Ny
          do i=1, Nx
             if (vice_temp(i,j) > 0._r8 .and. vice(i,j)>0._r8) then
                vicen(i,j,n) = vicen(i,j,n) - (vice_temp(i,j)-vice(i,j))*vicen(i,j,n)/vice_temp(i,j)
             endif  
          enddo
       enddo
    enddo
 
    ! calculate orignal caterogy thickness values
    aicen_temp = aicen_original
    where(aicen_temp==0) aicen_temp = -999
 
    do n=1,Ncat
       do j=1,Ny
          do i=1,Nx
             hicen_original(i,j,n) = vicen_original(i,j,n)/aicen_temp(i,j,n)
             hsnon_original(i,j,n) = vsnon_original(i,j,n)/aicen_temp(i,j,n)
          end do
       end do
    end do
 
    where(hicen_original < 0)  hicen_original = 0.0_r8
    where(hsnon_original < 0)  hsnon_original = 0.0_r8
   
    ! calculate the area implied by original category thickness and updated volume
    where(hicen_original /= 0.0_r8)
        aicen = vicen/hicen_original
    else where
        aicen = 0.0_r8
    end where
    
    aice = aicen(:,:,1)
    do n = 2, Ncat  
       aice = aice+aicen(:,:,n)
    enddo
 
    ! now squeeze aicen implied by original category thickness and updated volume
    do j = 1, Ny
       do i = 1, Nx
          if (aice(i,j) > 1.0_r8) then
             squeeze = 1.0_r8/aice(i,j)
             aicen(i,j,:) = aicen(i,j,:)*squeeze
          endif
       enddo
    enddo
 
    ! recalculate volume and snow volume with squeezed vicen
    vicen = aicen*hicen_original
    vsnon = aicen*hsnon_original
 
    ! consider special cases
    do n = 1, Ncat
       do j = 1, Ny
          do i = 1, Nx
          ! If there is no ice post-adjustment... 
             if (vicen(i,j,n)==0._r8 ) then
                aicen(i,j,n)   = 0._r8
                vsnon(i,j,n)   = 0._r8
                sice001(i,j,n) = 0._r8
                sice002(i,j,n) = 0._r8
                sice003(i,j,n) = 0._r8
                sice004(i,j,n) = 0._r8
                sice005(i,j,n) = 0._r8
                sice006(i,j,n) = 0._r8
                sice007(i,j,n) = 0._r8
                sice008(i,j,n) = 0._r8
                qice001(i,j,n) = 0._r8
                qice002(i,j,n) = 0._r8
                qice003(i,j,n) = 0._r8
                qice004(i,j,n) = 0._r8
                qice005(i,j,n) = 0._r8
                qice006(i,j,n) = 0._r8
                qice007(i,j,n) = 0._r8
                qice008(i,j,n) = 0._r8
                qsno001(i,j,n) = 0._r8
                qsno002(i,j,n) = 0._r8
                qsno003(i,j,n) = 0._r8
                Tsfcn(i,j,n)   = -1.836_r8
             ! If the adjustment introduced new ice.. 
             else if (aicen(i,j,n)>0._r8 .and. aicen_original(i,j,n)==0._r8) then
                ! allow no snow volume or enthalpy
                vsnon(i,j,n) = 0._r8
                qsno001(i,j,n) = 0._r8
                qsno002(i,j,n) = 0._r8
                qsno003(i,j,n) = 0._r8
 
                ! require ice volume for thickness = category boundary midpoint
                aicen(i,j,n) =  vicen(i,j,n)/hcat_midpoint(n)
 
                ! salinity of mushy ice, see add_new_ice in ice_therm_itd.F90
                Si0new = sss - dSin0_frazil ! given our choice of sss
                sice001(i,j,n) = Si0new
                sice002(i,j,n) = Si0new
                sice003(i,j,n) = Si0new
                sice004(i,j,n) = Si0new
                sice005(i,j,n) = Si0new
                sice006(i,j,n) = Si0new
                sice007(i,j,n) = Si0new
                sice008(i,j,n) = Si0new
 
                ! temperature and enthalpy
                Ti          = min(liquidus_temperature_mush(Si0new/phi_init), -0.1_r8)
                qi0new      = enthalpy_mush(Ti, Si0new)
                qice001(i,j,n) = qi0new
                qice002(i,j,n) = qi0new
                qice003(i,j,n) = qi0new
                qice004(i,j,n) = qi0new
                qice005(i,j,n) = qi0new
                qice006(i,j,n) = qi0new
                qice007(i,j,n) = qi0new
                qice008(i,j,n) = qi0new
                Tsfcn(i,j,n)  = Ti
             endif
          enddo
       enddo
    enddo
 
 end subroutine volume_simple_squeeze
! -----------------------------------------------------------------------------
 
 subroutine cice_rebalancing(qice001, qice002,     &
                             qice003, qice004,     &
                             qice005, qice006,     &
                             qice007, qice008,     &
                             sice001, sice002,     &
                             sice003, sice004,     &
                             sice005, sice006,     &
                             sice007, sice008,     &
                             qsno001, qsno002,     &
                             qsno003, aicen,       &
                             vicen, vsnon,         &
                             aicen_original,       &
                             vicen_original,       &
                             vsnon_original,       &
                             Tsfcn,                &
                             Ncat, Nx, Ny)
 
    real(r8), intent(inout), dimension(Nx,Ny,Ncat) :: &
                                        aicen, vicen, vsnon,  &
                          qice001, qice002, qice003, qice004, &
                          qice005, qice006, qice007, qice008, &
                          sice001, sice002, sice003, sice004, &
                          sice005, sice006, sice007, sice008, &
                                   qsno001, qsno002, qsno003, &
                                                       Tsfcn
    real(r8), intent(in), dimension(Nx,Ny,Ncat) :: &
                                           aicen_original, &
                                           vicen_original, &
                                           vsnon_original   
    integer, intent(in) :: Ncat, Nx, Ny
 
    real(r8), dimension(Nx,Ny) :: aice, vice, vsno, aice_temp, vice_temp, vsno_temp 
    real(r8), dimension(0:Ncat) :: hin_max
    real(r8), dimension(Ncat) :: hcat_midpoint
    real(r8) :: squeeze, cc1, cc2, x1, Si0new, Ti, qsno_hold, qi0new, hicen
    real(r8), parameter :: Tsmelt = 0._r8,        &
                           cc3 = 3._r8,           &
                           c1 = 1._r8,            &
                           phi_init = 0.75_r8,    &
                           dSin0_frazil = 3.0_r8, &
                           sss = 34.7_r8
 
    integer :: i, j, n ! loop indices

    ! calculate dependent variables 
    cc1 = cc3/real(Ncat,kind=r8)
    cc2 = 15.0_r8*cc1
    Si0new = sss - dSin0_frazil              
 
    ! calculate bounds and midpoints of thickness distribution
    hin_max(0) = 0.0_r8
    do n = 1, Ncat
       x1 = real(n-1,kind=r8) / real(Ncat,kind=r8)
       hin_max(n) = hin_max(n-1) &
                   + cc1 + cc2*(c1 + tanh(cc3*(x1-c1)))
       hcat_midpoint(n)=0.5_r8*(hin_max(n-1)+hin_max(n))
    enddo

    ! Begin process 
    call error_handler(E_MSG, source,  'beginning cice postprocessing process...')
    sice001  = max(0.0_r8, sice001)  ! salinities must be non-negative
    sice002  = max(0.0_r8, sice002)  ! salinities must be non-negative
    sice003  = max(0.0_r8, sice003)  ! salinities must be non-negative
    sice004  = max(0.0_r8, sice004)  ! salinities must be non-negative
    sice005  = max(0.0_r8, sice005)  ! salinities must be non-negative
    sice006  = max(0.0_r8, sice006)  ! salinities must be non-negative
    sice007  = max(0.0_r8, sice007)  ! salinities must be non-negative
    sice008  = max(0.0_r8, sice008)  ! salinities must be non-negative
    qice001  = min(0.0_r8, qice001)  ! enthalpies (ice) must be non-positive
    qice002  = min(0.0_r8, qice002)  ! enthalpies (ice) must be non-positive
    qice003  = min(0.0_r8, qice003)  ! enthalpies (ice) must be non-positive
    qice004  = min(0.0_r8, qice004)  ! enthalpies (ice) must be non-positive
    qice005  = min(0.0_r8, qice005)  ! enthalpies (ice) must be non-positive
    qice006  = min(0.0_r8, qice006)  ! enthalpies (ice) must be non-positive
    qice007  = min(0.0_r8, qice007)  ! enthalpies (ice) must be non-positive
    qice008  = min(0.0_r8, qice008)  ! enthalpies (ice) must be non-positive
    qsno001  = min(0.0_r8, qsno001)  ! enthalphies (snow) must be non-positive
    qsno002  = min(0.0_r8, qsno002)  ! enthalphies (snow) must be non-positive
    qsno003  = min(0.0_r8, qsno003)  ! enthalphies (snow) must be non-positive
    Tsfcn = min(Tsmelt,Tsfcn)  ! ice/snow surface must not exceed melting
    aicen = min(1.0_r8,aicen)  ! concentrations must not exceed 1 
    
    ! calculate aggregates for post-adjustment category variables  
    call error_handler(E_MSG, source,  'calculating aggregates...')
    aice = aicen(:,:,1)
    vice = vicen(:,:,1)
    vsno = vsnon(:,:,1)
    do n = 2, Ncat  
       aice = aice+aicen(:,:,n)
       vice = vice+vicen(:,:,n)
       vsno = vsno+vsnon(:,:,n)
    enddo
 
    ! impose bounds on categories
    call error_handler(E_MSG, source,  'imposing bounds on categories...')
    aicen = max(0.0_r8,aicen) ! concentration must be non-negative
    vicen = max(0.0_r8,vicen) ! volumes (ice) must be non-negative
    vsnon = max(0.0_r8,vsnon) ! volumes (snow) must be non-negative
 
    ! re-calculate aggregates once bounds are enforced
    call error_handler(E_MSG, source,  'recalculating aggregates...')
    aice_temp = aicen(:,:,1)
    vice_temp = vicen(:,:,1)
    vsno_temp = vsnon(:,:,1)
    do n = 2, Ncat  
       aice_temp = aice_temp+aicen(:,:,n)
       vice_temp = vice_temp+vicen(:,:,n)
       vsno_temp = vsno_temp+vsnon(:,:,n)
    enddo

    call error_handler(E_MSG, source,  'begin squeezing...')
    do j = 1, Ny
       do i = 1, Nx
       ! if the post-adjustment concentartion was 0 or less than 0, remove all ice 
          if (aice(i,j) <= 0.0_r8) then
             aicen(i,j,:) = 0.0_r8
             vicen(i,j,:) = 0.0_r8
             vsnon(i,j,:) = 0.0_r8
          else if (aice(i,j) > 0.0_r8) then
             do n=1,Ncat
                if (aice_temp(i,j) > 0.0_r8 .and. aice(i,j) > 0.0_r8) then
                   aicen(i,j,n) = aicen(i,j,n) - (aice_temp(i,j) - aice(i,j))*aicen(i,j,n)/aice_temp(i,j) 
                endif
                if (vice_temp(i,j) > 0.0_r8 .and. vice(i,j) > 0.0_r8) then
                   vicen(i,j,n) = vicen(i,j,n) - (vice_temp(i,j) - vice(i,j))*vicen(i,j,n)/vice_temp(i,j)
                endif
                if (vsno_temp(i,j) > 0.0_r8 .and. vsno(i,j) > 0.0_r8) then
                   vsnon(i,j,n) = vsnon(i,j,n) - (vsno_temp(i,j) - vsno(i,j))*vsnon(i,j,n)/vsno_temp(i,j)
                endif
                if (aicen(i,j,n) > 0.0_r8) then
                   hicen = vicen(i,j,n)/aicen(i,j,n)
                   if (n == Ncat) then
                      if (hicen < hin_max(n-1)) then
                          aicen(i,j,n) = vicen(i,j,n)/hin_max(n-1)
                      endif
                   else
                      if (hicen > hin_max(n) .or. hicen < hin_max(n-1)) then
                         aicen(i,j,n) = vicen(i,j,n)/hcat_midpoint(n)
                      endif
                   endif
                else
                   vicen(i,j,n) = 0.0_r8
                   vsnon(i,j,n) = 0.0_r8
                endif
             enddo
 
             ! recalculate the aggregate area
             aice = aicen(:,:,1)
             do n = 2, Ncat  
                aice = aice+aicen(:,:,n)
             enddo
 
             ! If the post-adjustment concentration is greater than 1, squeeze it down
             if (aice(i,j) > 1.0_r8) then
                squeeze = 1.0_r8/aice(i,j)
                aicen(i,j,:) = aicen(i,j,:)*squeeze
             endif        
          endif
 
          !! if ice exists in both the post-adjustment and post-bounds variables, 
          !! shift the post-adjustment negative values of each category variable 
          !do n=1,Ncat
          !   if (aice_temp(i,j) > 0.0_r8 .and. aice(i,j) > 0.0_r8) then
          !      aicen(i,j,n) = aicen(i,j,n) - (aice_temp(i,j) - aice(i,j))*aicen(i,j,n)/aice_temp(i,j)
          !   endif
          !   if (vice_temp(i,j) > 0.0_r8 .and. vice(i,j) > 0.0_r8) then
          !      vicen(i,j,n) = vicen(i,j,n) - (vice_temp(i,j) - vice(i,j))*vicen(i,j,n)/vice_temp(i,j)
          !   endif
          !   if (vsno_temp(i,j) > 0.0_r8 .and. vsno(i,j) > 0.0_r8) then
          !      vsnon(i,j,n) = vsnon(i,j,n) - (vsno_temp(i,j) - vsno(i,j))*vsnon(i,j,n)/vsno_temp(i,j)
          !   endif
          !enddo
 
          !! If the post-adjustment concentration is greater than 1, squeeze it down
          !if (aice(i,j) > 1.0_r8) then
          !   squeeze = 1.0_r8/aice(i,j)
          !   aicen(i,j,:) = aicen(i,j,:)*squeeze
          !endif
 
          ! Adjust the volume, snow, salinities and enthalphies to be consistent with the squeezed concentrations
          do n=1, Ncat
             ! if the adjustment and the original category both have ice in them... 
             if (aicen(i,j,n) > 0.0_r8 .and. aicen_original(i,j,n) > 0.0_r8) then
                ! calculate the volume corresponding to the area and midpoint thickness, if there's no volume
                if (vicen(i,j,n) == 0.0_r8) vicen(i,j,n) = aicen(i,j,n)*hcat_midpoint(n)
                ! calculate the enthalphy required to accomodate any new snow in the category
                if (vsnon(i,j,n) > 0.0_r8 .and. vsnon_original(i,j,n) == 0.0_r8) then
                   Ti = min(liquidus_temperature_mush(Si0new/phi_init), -0.1_r8)
                   qsno_hold = snow_enthaply(Ti)
                   qsno001(i,j,n) = qsno_hold
                   qsno002(i,j,n) = qsno_hold
                   qsno003(i,j,n) = qsno_hold
                endif
                ! if the adjustment doesn't have ice but the original does...
             else if (aicen(i,j,n) == 0.0_r8 .and. aicen_original(i,j,n) > 0.0_r8) then
                vicen(i,j,n) = 0.0_r8
                sice001(i,j,n) = 0._r8
                sice002(i,j,n) = 0._r8
                sice003(i,j,n) = 0._r8
                sice004(i,j,n) = 0._r8
                sice005(i,j,n) = 0._r8
                sice006(i,j,n) = 0._r8
                sice007(i,j,n) = 0._r8
                sice008(i,j,n) = 0._r8
                qice001(i,j,n) = 0._r8
                qice002(i,j,n) = 0._r8
                qice003(i,j,n) = 0._r8
                qice004(i,j,n) = 0._r8
                qice005(i,j,n) = 0._r8
                qice006(i,j,n) = 0._r8
                qice007(i,j,n) = 0._r8
                qice008(i,j,n) = 0._r8
                qsno001(i,j,n) = 0._r8
                qsno002(i,j,n) = 0._r8
                qsno003(i,j,n) = 0._r8
                vsnon(i,j,n) = 0.0_r8
                Tsfcn(i,j,n) = -1.836_r8
             ! if the adjustment has ice but the original doesn't... 
             else if (aicen(i,j,n)>0.0_r8 .and. aicen_original(i,j,n) == 0.0_r8) then
                if (vicen(i,j,n) == 0.0_r8) vicen(i,j,n) =  aicen(i,j,n) * hcat_midpoint(n)
                sice001(i,j,n) = Si0new
                sice002(i,j,n) = Si0new
                sice003(i,j,n) = Si0new
                sice004(i,j,n) = Si0new
                sice005(i,j,n) = Si0new
                sice006(i,j,n) = Si0new
                sice007(i,j,n) = Si0new
                sice008(i,j,n) = Si0new
                Ti = min(liquidus_temperature_mush(Si0new/phi_init), -0.1_r8)
                qi0new = enthalpy_mush(Ti, Si0new)
                qice001(i,j,n) = qi0new
                qice002(i,j,n) = qi0new
                qice003(i,j,n) = qi0new
                qice004(i,j,n) = qi0new
                qice005(i,j,n) = qi0new
                qice006(i,j,n) = qi0new
                qice007(i,j,n) = qi0new
                qice008(i,j,n) = qi0new
 
                if (vsnon(i,j,n) == 0.0_r8 .and. vsnon_original(i,j,n) > 0.0_r8) then
                   qsno001(i,j,n) = 0._r8
                   qsno002(i,j,n) = 0._r8
                   qsno003(i,j,n) = 0._r8
                else if (vsnon(i,j,n) > 0.0_r8 .and. vsnon_original(i,j,n) == 0.0_r8) then
                   qsno_hold = snow_enthaply(Ti)
                   qsno001(i,j,n) = qsno_hold
                   qsno002(i,j,n) = qsno_hold
                   qsno003(i,j,n) = qsno_hold
                endif
                Tsfcn(i,j,n) = Ti
             ! If neither the adjustment nor the original category have ice in them... 
             else if (aicen(i,j,n) == 0.0_r8) then
                vicen(i,j,n) = 0.0_r8
                vsnon(i,j,n) = 0.0_r8
             endif
          enddo
       enddo
    enddo
 
 call error_handler(E_MSG, source,  'finishing squeezing and reinitalizing associated variables...')
 
 end subroutine cice_rebalancing
 !------------------------------------------------------------------------

!------------------------------------------------------------------
! FUNCTIONS       
!------------------------------------------------------------------
function enthalpy_mush(zTin, zSin) result(zqin)

   ! enthalpy of mush from mush temperature and bulk salinity

   real(r8), intent(in) :: &
        zTin, & ! ice layer temperature (C)
        zSin    ! ice layer bulk salinity (ppt)

   real(r8) :: &
        zqin    ! ice layer enthalpy (J m-3) 

   real(r8) :: &
        phi     ! ice liquid fraction 

! from shr_const_mod.F90
   real(r8),parameter :: SHR_CONST_CPSW  = 3.996e3_R8   ! specific heat of sea water ~ J/kg/K
   real(R8),parameter :: SHR_CONST_CPICE = 2.11727e3_R8 ! specific heat of fresh ice ~ J/kg/K
   real(R8),parameter :: SHR_CONST_RHOSW = 1.026e3_R8   ! density of sea water ~ kg/m^3
   real(R8),parameter :: SHR_CONST_RHOICE= 0.917e3_R8   ! density of ice        ~ kg/m^3
   real(R8),parameter :: SHR_CONST_LATICE= 3.337e5_R8   ! latent heat of fusion ~ J/kg


! from cice/src/drivers/cesm/ice_constants.F90
   real(r8) :: cp_ocn, cp_ice, rhoi, rhow, Lfresh

   cp_ice    = SHR_CONST_CPICE  ! specific heat of fresh ice (J/kg/K)
   cp_ocn    = SHR_CONST_CPSW   ! specific heat of ocn    (J/kg/K)
   rhoi      = SHR_CONST_RHOICE ! density of ice (kg/m^3)
   rhow      = SHR_CONST_RHOSW  ! density of seawater (kg/m^3)
   Lfresh    = SHR_CONST_LATICE ! latent heat of melting of fresh ice (J/kg)

   phi = liquid_fraction(zTin, zSin)

   zqin = phi * (cp_ocn * rhow - cp_ice * rhoi) * zTin + &
          rhoi * cp_ice * zTin - (1._r8 - phi) * rhoi * Lfresh

 end function enthalpy_mush

 function liquid_fraction(zTin, zSin) result(phi)

   ! liquid fraction of mush from mush temperature and bulk salinity

   real(r8), intent(in) :: &
        zTin, & ! ice layer temperature (C)
        zSin    ! ice layer bulk salinity (ppt)

   real(r8) :: &
        phi , & ! liquid fraction
        Sbr     ! brine salinity (ppt)

   real (r8), parameter :: puny = 1.0e-11_r8 ! cice/src/drivers/cesm/ice_constants.F90

   Sbr = max(liquidus_brine_salinity_mush(zTin),puny)
   phi = zSin / max(Sbr, zSin)

 end function liquid_fraction

 function snow_enthaply(Ti) result(qsno)
   real(r8), intent(in) :: Ti

   real(r8),parameter :: rhos = 330.0_r8, &
                       Lfresh = 2.835e6_r8 - 2.501e6_r8, &
                       cp_ice = 2106._r8
   real(r8) :: qsno

   qsno = -rhos*(Lfresh - cp_ice*min(0.0_r8,Ti))
end function snow_enthaply

function liquidus_brine_salinity_mush(zTin) result(Sbr)

   ! liquidus relation: equilibrium brine salinity as function of temperature
   ! based on empirical data from Assur (1958)

   real(r8), intent(in) :: &
        zTin         ! ice layer temperature (C)

   real(r8) :: &
        Sbr          ! ice brine salinity (ppt)

   real(r8) :: &
        t_high   , & ! mask for high temperature liquidus region
        lsubzero     ! mask for sub-zero temperatures

   !constant numbers from ice_constants.F90
   real(r8), parameter :: &
        c1      = 1.0_r8 , &
        c1000   = 1000_r8

   ! liquidus relation - higher temperature region
   real(r8), parameter :: &
        az1_liq = -18.48_r8 ,&
        bz1_liq =   0.0_r8

   ! liquidus relation - lower temperature region
   real(r8), parameter :: &
        az2_liq = -10.3085_r8,  &
        bz2_liq =  62.4_r8

   ! liquidus break
   real(r8), parameter :: &
        Tb_liq = -7.6362968855167352_r8
        
   ! basic liquidus relation constants
   real(r8), parameter :: &
        az1p_liq = az1_liq / c1000, &
        bz1p_liq = bz1_liq / c1000, &
        az2p_liq = az2_liq / c1000, &
        bz2p_liq = bz2_liq / c1000

   ! temperature to brine salinity
   real(r8), parameter :: &
      J1_liq = bz1_liq / az1_liq         , &
      K1_liq = c1 / c1000                , &
      L1_liq = (c1 + bz1p_liq) / az1_liq , &
      J2_liq = bz2_liq  / az2_liq        , &
      K2_liq = c1 / c1000                , &
      L2_liq = (c1 + bz2p_liq) / az2_liq

   t_high   = merge(1._r8, 0._r8, (zTin > Tb_liq))
   lsubzero = merge(1._r8, 0._r8, (zTin <= 1._r8))

   Sbr = ((zTin + J1_liq) / (K1_liq * zTin + L1_liq)) * t_high + &
         ((zTin + J2_liq) / (K2_liq * zTin + L2_liq)) * (1._r8 - t_high)

   Sbr = Sbr * lsubzero

end function liquidus_brine_salinity_mush

function liquidus_temperature_mush(Sbr) result(zTin)

   ! liquidus relation: equilibrium temperature as function of brine salinity
   ! based on empirical data from Assur (1958)

   real(r8), intent(in) :: &
        Sbr    ! ice brine salinity (ppt)

   real(r8) :: &
        zTin   ! ice layer temperature (C)

   real(r8) :: &
        t_high ! mask for high temperature liquidus region

   ! liquidus break
   real(r8), parameter :: &
      Sb_liq =  123.66702800276086_r8    ! salinity of liquidus break

   ! constant numbers from ice_constants.F90
   real(r8), parameter :: &
        c1      = 1.0_r8 , &
        c1000   = 1000_r8

   ! liquidus relation - higher temperature region
   real(r8), parameter :: &
        az1_liq = -18.48_r8 ,&
        bz1_liq =   0.0_r8

   ! liquidus relation - lower temperature region
   real(r8), parameter :: &
        az2_liq = -10.3085_r8,  &
        bz2_liq =  62.4_r8

   ! basic liquidus relation constants
   real(r8), parameter :: &
        az1p_liq = az1_liq / c1000, &
        bz1p_liq = bz1_liq / c1000, &
        az2p_liq = az2_liq / c1000, &
        bz2p_liq = bz2_liq / c1000

 ! brine salinity to temperature
   real(r8), parameter :: &
      M1_liq = az1_liq            , &
      N1_liq = -az1p_liq          , &
      O1_liq = -bz1_liq / az1_liq , &
      M2_liq = az2_liq            , &
      N2_liq = -az2p_liq          , &
      O2_liq = -bz2_liq / az2_liq

   t_high = merge(1._r8, 0._r8, (Sbr <= Sb_liq))

   zTin = ((Sbr / (M1_liq + N1_liq * Sbr)) + O1_liq) * t_high + &
         ((Sbr / (M2_liq + N2_liq * Sbr)) + O2_liq) * (1._r8 - t_high)

end function liquidus_temperature_mush


end module ice_postprocessing_mod
