! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program dart_to_cice

!----------------------------------------------------------------------
! purpose: implement a 'partition function' to modify the icepack state
!          to be consistent with the states from assimilation
!
! method: Read in restart (restart with prior) and out restart (restart
!         with posterior) written by DART after filter.
!
! author: M Wieringa (2023) based on C Bitz (2016) and C Riedel (2022)
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             file_exist, error_handler, E_ERR, E_MSG, to_upper
use  netcdf_utilities_mod, only : nc_check
use netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"

!------------------------------------------------------------------
! SET NAMELIST AND ALLOCATE VARIABLES    
!------------------------------------------------------------------
character(len=256) :: dart_to_cice_input_file = 'post_filter_restart.nc'
character(len=256) :: original_cice_restart_file = 'pre_filter_restart.nc'
character(len=256) :: postprocessed_output_file = 'postprocessed_restart.nc'
character(len=128) :: balance_method = 'simple_squeeze'
character(len=128) :: postprocess = 'cice'
character(len=15)  :: r_snw_name  = 'r_snw'
integer :: gridpt_oi = 3

namelist /dart_to_cice_nml/ dart_to_cice_input_file,    &
                            original_cice_restart_file, &
                            postprocessed_output_file,  &
                            balance_method,             &
                            postprocess,                &
                            r_snw_name,                 &
                            gridpt_oi

! general variable iniatlization
character(len=512) :: string1, string2, msgstring
character(len=15)  :: varname
character(len=128) :: method
character(len=3)   :: nchar

integer :: iunit, io, ncid, dimid, l, n, VarID, Ncat, Nx 
integer, parameter :: Nilyr = 8   ! number of layers in ice, hardwired
integer, parameter :: Nslyr = 3   ! number of layers in snow, hardwired

real(r8), allocatable :: aicen_original(:), vicen_original(:), vsnon_original(:)
real(r8), allocatable :: aicen(:), vicen(:), vsnon(:), Tsfcn(:)
real(r8), allocatable :: qice(:,:)
real(r8), allocatable :: sice(:,:)
real(r8), allocatable :: qsno(:,:)

!------------------------------------------------------------------
! INIALIZE AND PERFORM CHECKS ON FILES                
!------------------------------------------------------------------
call initialize_utilities(progname='dart_to_cice')

call find_namelist_in_file("input.nml", "dart_to_cice_nml", iunit)
read(iunit, nml = dart_to_cice_nml, iostat = io)
call check_namelist_read(iunit, io, "dart_to_cice_nml")

method = balance_method
call to_upper(method)

write(string1,*) 'converting DART output file "'// &
                 &trim(dart_to_cice_input_file)//'" to one CICE will like'
write(string2,*) 'using the "'//trim(balance_method)//'" method.'
call error_handler(E_MSG,'dart_to_cice',string1,text2=string2)

if ( .not. file_exist(dart_to_cice_input_file) ) then
   write(string1,*) 'cannot open "', trim(dart_to_cice_input_file),'" for updating.'
   call error_handler(E_ERR,'dart_to_cice:filename not found ',trim(dart_to_cice_input_file))
endif

if ( .not. file_exist(original_cice_restart_file) ) then
   write(string1,*) 'cannot open "', trim(original_cice_restart_file),'" for reading.'
   call error_handler(E_ERR,'dart_to_cice:filename not found ',trim(original_cice_restart_file))
endif

!------------------------------------------------------------------
! READ VARIABLES FROM RESTART FILES               
!------------------------------------------------------------------
! Read the pre-assim variables
call nc_check( nf90_open(trim(original_cice_restart_file), NF90_NOWRITE, ncid), &
                  'dart_to_cice', 'open "'//trim(original_cice_restart_file)//'"')

call nc_check(nf90_inq_dimid(ncid, "ncat", dimid), &
                  'dart_to_cice', 'inquire ncat dimid from "'//trim(original_cice_restart_file)//'"')
call nc_check(nf90_inquire_dimension(ncid, dimid, len=Ncat), &
                  'dart_to_cice', 'inquire ncat from "'//trim(original_cice_restart_file)//'"')
call nc_check(nf90_inq_dimid(ncid,"ni",dimid), &
                   'dart_to_cice', 'inquire ni dimid from "'//trim(original_cice_restart_file)//'"')
call nc_check(nf90_inquire_dimension(ncid,dimid,len=Nx),&
                   'dart_to_cice', 'inquire ni from "'//trim(original_cice_restart_file)//'"')


allocate(aicen(Ncat), vicen(Ncat), vsnon(Ncat), Tsfcn(Ncat), aicen_original(Ncat), vicen_original(Ncat), vsnon_original(Ncat))
allocate(qice(Nilyr, Ncat), sice(Nilyr, Ncat), qsno(Nslyr, Ncat))

call get_variable(ncid, 'aicen', aicen_original, original_cice_restart_file, gridpt_oi, Ncat)
call get_variable(ncid, 'vicen', vicen_original, original_cice_restart_file, gridpt_oi, Ncat)
call get_variable(ncid, 'vsnon', vsnon_original, original_cice_restart_file, gridpt_oi, Ncat)
call get_variable(ncid, 'Tsfcn', Tsfcn, original_cice_restart_file, gridpt_oi, Ncat)
do l=1, Nilyr
  write(nchar,'(i3.3)') l
  call get_variable(ncid, 'qice'//trim(nchar), qice(l,:), original_cice_restart_file, gridpt_oi, Ncat)
  call get_variable(ncid, 'sice'//trim(nchar), sice(l,:), original_cice_restart_file, gridpt_oi, Ncat)
enddo
do l=1, Nslyr
  write(nchar,'(i3.3)') l
  call get_variable(ncid, 'qsno'//trim(nchar), qsno(l,:), original_cice_restart_file, gridpt_oi, Ncat)
enddo

call nc_check(nf90_close(ncid),'dart_to_cice', 'close '//trim(original_cice_restart_file))

! Read the post-assim variables 
call nc_check( nf90_open(trim(dart_to_cice_input_file), NF90_NOWRITE, ncid), &
                  'dart_to_cice', 'open "'//trim(dart_to_cice_input_file)//'"')

call get_variable(ncid, 'aicen', aicen, dart_to_cice_input_file, gridpt_oi, Ncat)
call get_variable(ncid, 'vicen', vicen, dart_to_cice_input_file, gridpt_oi, Ncat)
call get_variable(ncid, 'vsnon', vsnon, dart_to_cice_input_file, gridpt_oi, Ncat)

call nc_check(nf90_close(ncid),'dart_to_cice', 'close '//trim(dart_to_cice_input_file))

! write(*,*) 'aicen dimensions are ', shape(aicen)
!------------------------------------------------------------------
! PERFORM POSTPROCESSING 
!------------------------------------------------------------------
if (postprocess == 'cice') then
    write(*,*) 'calling cice postprocessing...'
    call cice_rebalancing(qice, sice, qsno,     &
                          aicen, vicen, vsnon,  &
                          aicen_original,       &
                          vicen_original,       &
                          vsnon_original,       &
                          Tsfcn,               &
                          Ncat, Nilyr, Nslyr)
    write(*,*) 'cice postprocessing function completed...'

  else if (postprocess == 'aice') then 
    write(*,*) 'calling aice postprocessing...'
    call area_simple_squeeze(qice, sice, qsno,     &
                             aicen, vicen, vsnon,  &
                             aicen_original,       &
                             vicen_original,       &
                             vsnon_original,       &
                             Tsfcn,                &
                             Ncat, Nilyr, Nslyr)
    write(*,*) 'aice postprocessing function completed...'
  else if (postprocess == 'vice') then
    write(*,*) 'calling vice postprocessing...'
    call volume_simple_squeeze(qice, sice, qsno,     &
                               aicen, vicen, vsnon,  &
                               aicen_original,       &
                               vicen_original,       &
                               vsnon_original,       &
                               Tsfcn,                &
                               Ncat, Nilyr, Nslyr)
    write(*,*) 'vice postprocessing function completed...'
  else
    write(*,*) 'No valid postprocessing method called. No adjustments will be made.'
end if

!------------------------------------------------------------------
! WRITE VARIABLES TO RESTART FILE
!------------------------------------------------------------------
call nc_check( nf90_open(trim(postprocessed_output_file), NF90_WRITE, ncid), &
                  'dart_to_cice', 'open "'//trim(postprocessed_output_file)//'"')

varname='aicen'
io = nf90_inq_varid(ncid, trim(varname), VarID)
call nc_check(io, 'dart_to_cice', 'inq_varid '//trim(varname)//' '//trim(postprocessed_output_file))
io = nf90_put_var(ncid, VarID, aicen, start=(/gridpt_oi,1/), count=(/1,Ncat/))
call nc_check(io, 'dart_to_cice', 'put_var '//trim(varname)//' '//trim(postprocessed_output_file))

varname='vicen'
io = nf90_inq_varid(ncid, trim(varname), VarID)
call nc_check(io, 'dart_to_cice', 'inq_varid '//trim(varname)//' '//trim(postprocessed_output_file))
io = nf90_put_var(ncid, VarID, vicen, start=(/gridpt_oi,1/), count=(/1,Ncat/))
call nc_check(io, 'dart_to_cice', 'put_var '//trim(varname)//' '//trim(postprocessed_output_file))

varname='vsnon'
io = nf90_inq_varid(ncid, trim(varname), VarID)
call nc_check(io, 'dart_to_cice', 'inq_varid '//trim(varname)//' '//trim(postprocessed_output_file))
io = nf90_put_var(ncid, VarID, vsnon, start=(/gridpt_oi,1/), count=(/1,Ncat/))
call nc_check(io, 'dart_to_cice', 'put_var '//trim(varname)//' '//trim(postprocessed_output_file))

varname='Tsfcn'
io = nf90_inq_varid(ncid, trim(varname), VarID)
call nc_check(io, 'dart_to_cice', 'inq_varid '//trim(varname)//' '//trim(postprocessed_output_file))
io = nf90_put_var(ncid, VarID, Tsfcn, start=(/gridpt_oi,1/), count=(/1,Ncat/))
call nc_check(io, 'dart_to_cice', 'put_var '//trim(varname)//' '//trim(postprocessed_output_file))

do l=1, Nilyr
  write(nchar,'(i3.3)') l
  varname='qice'//trim(nchar)
  io = nf90_inq_varid(ncid, trim(varname), VarID)
  call nc_check(io, 'dart_to_cice', 'inq_varid '//trim(varname)//' '//trim(postprocessed_output_file))
  io = nf90_put_var(ncid, VarID, qice(l,:), start=(/gridpt_oi,1/), count=(/1,Ncat/))
  call nc_check(io, 'dart_to_cice', 'put_var '//trim(varname)//' '//trim(postprocessed_output_file))
 
  varname='sice'//trim(nchar)
  io = nf90_inq_varid(ncid, trim(varname), VarID)
  call nc_check(io, 'dart_to_cice', 'inq_varid '//trim(varname)//' '//trim(postprocessed_output_file))
  io = nf90_put_var(ncid, VarID, sice(l,:), start=(/gridpt_oi,1/), count=(/1,Ncat/))
  call nc_check(io, 'dart_to_cice', 'put_var '//trim(varname)//' '//trim(postprocessed_output_file))
enddo

do l=1, Nslyr
  write(nchar,'(i3.3)') l
  varname='qsno'//trim(nchar)
  io = nf90_inq_varid(ncid, trim(varname), VarID)
  call nc_check(io, 'dart_to_cice', 'inq_varid '//trim(varname)//' '//trim(postprocessed_output_file))
  io = nf90_put_var(ncid, VarID, qsno(l,:), start=(/gridpt_oi,1/), count=(/1,Ncat/))
  call nc_check(io, 'dart_to_cice', 'put_var '//trim(varname)//' '//trim(postprocessed_output_file))
enddo

call nc_check(nf90_close(ncid),'dart_to_cice', 'close '//trim(postprocessed_output_file))

!------------------------------------------------------------------
! DEALLOCATE AND FINALIZE                
!------------------------------------------------------------------
deallocate(aicen, vicen, vsnon, Tsfcn, aicen_original, vicen_original, vsnon_original)
deallocate(qice, sice, qsno )

call finalize_utilities('dart_to_cice')


contains

!------------------------------------------------------------------
! SUBROUTINES            
!------------------------------------------------------------------
subroutine get_variable(ncid, varname, var, filename, space_index, Ncat)
  integer,               intent(in)  :: ncid, Ncat
  character(len=*),      intent(in)  :: varname
  real(r8),              intent(out) :: var(Ncat)
  character(len=*),      intent(in)  :: filename
  integer,               intent(in)  :: space_index
  
  integer :: VarID, ndims, dimIDs
  real(r8) :: holder(4,ncat)
  
  write(6,*) 'Getting data for ',trim(varname)
  
  io = nf90_inq_varid(ncid, trim(varname), VarID)
  call nc_check(io, 'dart_to_cice', 'inq_varid '//trim(msgstring))
  
  call nc_check(nf90_get_var(ncid, VarID, holder), 'dart_to_cice', &
           'get_var '//trim(msgstring))
  
  
  var(:) = holder(space_index,:)
  
end subroutine get_variable

subroutine area_simple_squeeze(qice, sice, qsno,    &
                               aicen, vicen, vsnon, &
                               aicen_original,      &
                               vicen_original,      &
                               vsnon_original,      &
                               Tsfcn,               &
                               Ncat, Nilyr, Nslyr)
  
  real(r8), intent(inout) :: aicen(Ncat), vicen(Ncat), vsnon(Ncat),  &
                             qice(Nilyr,Ncat), sice(Nilyr,Ncat), qsno(Nslyr,Ncat),     &
                             Tsfcn(Ncat)
  real(r8), intent(in) :: aicen_original(Ncat), vicen_original(Ncat), vsnon_original(Ncat)     
  integer, intent(in) :: Ncat, Nilyr, Nslyr

  real(r8) :: aice, aice_temp
  real(r8) :: hin_max(0:Ncat)
  real(r8) :: hcat_midpoint(Ncat), hicen_original(Ncat), hsnon_original(Ncat)
  real(r8) :: aicen_temp(Ncat)
  real(r8) :: squeeze, cc1, cc2, x1, Si0new, Ti, qsno_hold, qi0new 
  real(r8), parameter :: Tsmelt = 0._r8,        &
                         cc3 = 3._r8,           &
                         c1 = 1._r8,            &
                         phi_init = 0.75_r8,    &
                         dSin0_frazil = 3.0_r8, &
                         sss = 34.7_r8

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
    sice  = max(0.0_r8, sice)  ! salinities must be non-negative
    qice  = min(0.0_r8, qice)  ! enthalpies (ice) must be non-positive
    qsno  = min(0.0_r8, qsno)  ! enthalphies (snow) must be non-positive
    aicen = min(1.0_r8,aicen)  ! concentrations must not exceed 1 
    Tsfcn = min(Tsmelt,Tsfcn)  ! ice/snow surface must not exceed melting

  ! calculate aice, which might be negative or >1 at this point
    aice = sum(aicen)
   
  ! set negative aicen to zero
    aicen = max(0.0_r8,aicen)   ! concentrations must be non-negative
    vicen = max(0.0_r8,vicen)   ! same for volumes (ice)
    vsnon = max(0.0_r8,vsnon)   ! same for volumes (snow)

  ! reclaculate aice, now it should be non-negative
    aice_temp = sum(aicen)
 
  ! if aice <0, then set every category to 0
    if (aice < 0._r8) then
      aicen(:) = 0._r8
    end if

  ! shift negative concentration values 
    do n=1, Ncat
      if (aice_temp > 0._r8 .and. aice > 0._r8) then
        aicen(n) = aicen(n) - (aice_temp-aice)*aicen(n)/aice_temp
      endif  
    enddo
 
  ! now squeeze aicen 
    if (aice > 1.0_r8) then
      squeeze  = 1.0_r8/aice
      aicen(:) = aicen(:)*squeeze
    endif

  ! update vsnon and vicen using conserved category thickness values
    aicen_temp = aicen_original
    where(aicen_temp==0) aicen_temp = -999
 
    do n=1,Ncat
      hicen_original(n) = vicen_original(n)/aicen_temp(n)
      hsnon_original(n) = vsnon_original(n)/aicen_temp(n)
    end do
 
    where(hicen_original < 0)  hicen_original = 0.0_r8
    where(hsnon_original < 0)  hsnon_original = 0.0_r8

    vicen  = aicen*hicen_original
    vsnon  = aicen*hsnon_original

  ! consider special cases 
  do n = 1, Ncat
  ! If there is no ice post-adjustment... 
    if (aicen(n) == 0._r8) then
      vicen(n)  = 0._r8
      vsnon(n)  = 0._r8
      sice(:,n) = 0._r8
      qice(:,n) = 0._r8
      qsno(:,n) = 0._r8
      Tsfcn(n)  = -1.836_r8
  ! If the adjustment introduced new ice.. 
    else if (aicen(n) > 0._r8 .and. aicen_original(n) == 0._r8) then
      ! allow no snow volume or enthalpy
      vsnon(n) = 0._r8
      qsno(:,n) = 0._r8
 
      ! require ice volume for thickness = category boundary midpoint
      vicen(n) = aicen(n) * hcat_midpoint(n)
 
      ! salinity of mushy ice, see add_new_ice in ice_therm_itd.F90
      Si0new    = sss - dSin0_frazil ! given our choice of sss
      sice(:,n) = Si0new

      ! temperature and enthalphy 
      Ti        = min(liquidus_temperature_mush(Si0new/phi_init), -0.1_r8)
      qi0new    = enthalpy_mush(Ti, Si0new)
      qice(:,n) = qi0new
      Tsfcn(n) = Ti
    endif
  enddo

end subroutine area_simple_squeeze

subroutine volume_simple_squeeze(qice, sice, qsno,    &
                                 aicen, vicen, vsnon, &
                                 aicen_original,      &
                                 vicen_original,      &
                                 vsnon_original,      &
                                 Tsfcn,               &
                                 Ncat, Nilyr, Nslyr)
  
  real(r8), intent(inout) :: aicen(Ncat), vicen(Ncat), vsnon(Ncat),  &
                             qice(Nilyr,Ncat), sice(Nilyr,Ncat), qsno(Nslyr,Ncat),     &
                             Tsfcn(Ncat)
  real(r8), intent(in) :: aicen_original(Ncat), vicen_original(Ncat), vsnon_original(Ncat)     
  integer, intent(in) :: Ncat, Nilyr, Nslyr

  real(r8) :: aice, vice, vice_temp
  real(r8) :: hin_max(0:Ncat)
  real(r8) :: hcat_midpoint(Ncat), hicen_original(Ncat), hsnon_original(Ncat)
  real(r8) :: aicen_temp(Ncat)
  real(r8) :: squeeze, cc1, cc2, x1, Si0new, Ti, qsno_hold, qi0new 
  real(r8), parameter :: Tsmelt = 0._r8,        &
                         cc3 = 3._r8,           &
                         c1 = 1._r8,            &
                         phi_init = 0.75_r8,    &
                         dSin0_frazil = 3.0_r8, &
                         sss = 34.7_r8

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
    sice  = max(0.0_r8, sice)  ! salinities must be non-negative
    qice  = min(0.0_r8, qice)  ! enthalpies (ice) must be non-positive
    qsno  = min(0.0_r8, qsno)  ! enthalphies (snow) must be non-positive
    ! aicen = min(1.0_r8,aicen)  ! concentrations must not exceed 1 
    Tsfcn = min(Tsmelt,Tsfcn)  ! ice/snow surface must not exceed melting

  ! calculate vice, which might be negative or >1 at this point
    vice = sum(vicen)
   
  ! set negative values to zero
    vicen = max(0.0_r8,vicen)   ! same for volumes (ice)

  ! reclaculate vice, now it should be non-negative
    vice_temp = sum(vicen)
 
  ! if vice <0, then set every category to 0
    if (vice < 0._r8) then
      vicen(:) = 0._r8
    end if

  ! shift negative post-adjustment volume values 
    do n=1, Ncat
      if (vice_temp > 0._r8 .and. vice > 0._r8) then
        vicen(n) = vicen(n) - (vice_temp-vice)*vicen(n)/vice_temp
      endif  
    enddo
  
  ! calculate original category thickness values 
    aicen_temp = aicen_original
    where(aicen_temp==0) aicen_temp = -999
 
    do n=1,Ncat
      hicen_original(n) = vicen_original(n)/aicen_temp(n)
      hsnon_original(n) = vsnon_original(n)/aicen_temp(n)
    end do
 
    where(hicen_original < 0)  hicen_original = 0.0_r8
    where(hsnon_original < 0)  hsnon_original = 0.0_r8

  ! calculate the area implied by original category thickness and updated volume
    aicen = vicen/hicen_original
    aice  = sum(aicen)   

  ! now squeeze the aicen via the aice implied by original category thickness 
    if (aice > 1.0_r8) then
      squeeze  = 1.0_r8/aice
      aicen(:) = aicen(:) * squeeze
    endif

  ! recalculate volume and snow volume with squeezed vicen
    vicen = aicen*hicen_original
    vsnon = aicen*hsnon_original

  ! consider special cases 
  do n = 1, Ncat
  ! If there is no ice post-adjustment... 
    if (vicen(n) == 0._r8) then
      aicen(n)  = 0._r8
      vsnon(n)  = 0._r8
      sice(:,n) = 0._r8
      qice(:,n) = 0._r8
      qsno(:,n) = 0._r8
      Tsfcn(n)  = -1.836_r8
  ! If the adjustment introduced new ice.. 
    else if (vicen(n) > 0._r8 .and. vicen_original(n) == 0._r8) then
      ! allow no snow volume or enthalpy
      vsnon(n) = 0._r8
      qsno(:,n) = 0._r8
 
      ! require ice volume for thickness = category boundary midpoint
      aicen(n) = vicen(n)/hcat_midpoint(n)
 
      ! salinity of mushy ice, see add_new_ice in ice_therm_itd.F90
      Si0new    = sss - dSin0_frazil ! given our choice of sss
      sice(:,n) = Si0new

      ! temperature and enthalphy 
      Ti        = min(liquidus_temperature_mush(Si0new/phi_init), -0.1_r8)
      qi0new    = enthalpy_mush(Ti, Si0new)
      qice(:,n) = qi0new
      Tsfcn(n) = Ti
    endif
  enddo

end subroutine volume_simple_squeeze

subroutine cice_rebalancing(qice, sice, qsno,     &
                            aicen, vicen, vsnon,  &
                            aicen_original,       &
                            vicen_original,       &
                            vsnon_original,       &
                            Tsfcn,                &
                            Ncat, Nilyr, Nslyr)

  real(r8), intent(inout) :: aicen(Ncat), vicen(Ncat), vsnon(Ncat),  &
                             qice(Nilyr,Ncat), sice(Nilyr,Ncat), qsno(Nslyr,Ncat),     &
                             Tsfcn(Ncat)
  real(r8), intent(in) :: aicen_original(Ncat), vicen_original(Ncat), vsnon_original(Ncat)     
  integer, intent(in) :: Ncat, Nilyr, Nslyr

  real(r8) :: aice, vice, vsno, aice_temp, vice_temp, vsno_temp
  real(r8) :: hin_max(0:Ncat)
  real(r8) :: hcat_midpoint(Ncat)
  real(r8) :: squeeze, cc1, cc2, x1, Si0new, Ti, qsno_hold, qi0new 
  real(r8), parameter :: Tsmelt = 0._r8,        &
                         cc3 = 3._r8,           &
                         c1 = 1._r8,            &
                         phi_init = 0.75_r8,    &
                         dSin0_frazil = 3.0_r8, &
                         sss = 34.7_r8
  
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

  ! Begin process with the variables post-adjustment 
    sice = max(0.0_r8,sice)   ! salinities must be non-negative
    qice = min(0.0_r8,qice)   ! enthalpy (ice) must be non-positive
    qsno = min(0.0_r8,qsno)   ! enthalphy (snow) must be non-positive
    Tsfcn = min(Tsmelt,Tsfcn) ! surface temperature must be less than the melting point 
    aicen = min(1.0_r8,aicen) ! concentration must be less than 1

  ! calculate aggregates for post-adjustment category variables
    aice = sum(aicen)
    vice = sum(vicen)
    vsno = sum(vsnon)
  
  ! impose bounds on categories
    aicen = max(0.0_r8,aicen) ! concentration must be non-negative
    vicen = max(0.0_r8,vicen) ! volumes (ice) must be non-negative
    vsnon = max(0.0_r8,vsnon) ! volumes (snow) must be non-negative 

  ! re-calculate aggregates once bounds are enforced
    aice_temp = sum(aicen)
    vice_temp = sum(vicen)
    vsno_temp = sum(vsnon)

  ! If the post-adjustment concentration was 0 or less than 0, remove all ice  
    if (aice <= 0.0_r8) then
      aicen(:) = 0.0_r8
      vicen(:) = 0.0_r8
      vsnon(:) = 0.0_r8
    endif

  ! If ice exists in both the post-adjustment and post-bounds variables, 
  ! shift the post-adjustment negative values of each category variable 
    do n=1,Ncat
      if (aice_temp > 0._r8 .and. aice > 0._r8) then
         aicen(n) = aicen(n) - (aice_temp-aice)*aicen(n)/aice_temp
     endif
      if (vice_temp > 0._r8 .and. vice > 0._r8) then
       vicen(n) = vicen(n) - (vice_temp-vice)*vicen(n)/vice_temp
      endif
      if (vsno_temp > 0._r8 .and. vsno > 0._r8) then
       vsnon(n) = vsnon(n) - (vsno_temp-vsno)*vsnon(n)/vsno_temp
      endif
    enddo
   
  ! If the post-adjustment concentration is greater than 1, squeeze it down
    if (aice > 1.0_r8) then
      squeeze = 1.0_r8/aice
      aicen(:) = aicen(:)*squeeze
    endif

  ! Adjust the volume, snow, salinities and enthalphies to be consistent with the squeezed concentrations
    do n=1,Ncat
    ! if the adjustment and the original category both have ice in them... 
      if (aicen(n) > 0.0_r8 .and. aicen_original(n) > 0.0_r8) then
        ! calculate the volume corresponding to the area and midpoint thickness, if there's no volume
        if (vicen(n) == 0.0_r8) vicen(n) = aicen(n)*hcat_midpoint(n)
        ! calculate the enthalphy required to accomodate any new snow in the category
        if (vsnon(n) > 0.0_r8 .and. vsnon_original(n) == 0.0_r8) then
          Ti = min(liquidus_temperature_mush(Si0new/phi_init), -0.1_r8)
          qsno_hold = snow_enthaply(Ti)
          qsno(:,n) = qsno_hold
        endif
    ! if the adjustment doesn't have ice but the original does...
      else if (aicen(n) == 0.0_r8 .and. aicen_original(n) > 0.0_r8) then
        vicen(n) = 0.0_r8
        qice(:,n) = 0.0_r8
        sice(:,n) = 0.0_r8
        qsno(:,n) = 0.0_r8
        vsnon(n) = 0.0_r8
        Tsfcn(n) = -1.836_r8
    ! if the adjustment has ice but the original doesn't... 
      else if (aicen(n)>0.0_r8 .and. aicen_original(n) == 0.0_r8) then
        if (vicen(n) == 0.0_r8) vicen(n) =  aicen(n) * hcat_midpoint(n)
        sice(:,n) = Si0new
        Ti = min(liquidus_temperature_mush(Si0new/phi_init), -0.1_r8)
        qi0new = enthalpy_mush(Ti, Si0new)
        qice(:,n) = qi0new

        if (vsnon(n) == 0.0_r8 .and. vsnon_original(n) > 0.0_r8) then
          qsno(:,n) = 0.0_r8
        else if (vsnon(n) > 0.0_r8 .and. vsnon_original(n) == 0.0_r8) then
          qsno_hold = snow_enthaply(Ti)
          qsno(:,n) = qsno_hold
        endif
        Tsfcn(n) = Ti
    ! If neither the adjustment nor the original category have ice in them... 
      else if (aicen(n) == 0.0_r8) then
        vicen(n) = 0.0_r8
        vsnon(n) = 0.0_r8
      endif
    enddo
  
end subroutine cice_rebalancing

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

!------------------------------------------------------------------
! END             
!------------------------------------------------------------------
end program dart_to_cice

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
