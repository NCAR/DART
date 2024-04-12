! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program dart_to_cice

!----------------------------------------------------------------------
! purpose: implement a 'partition function' to modify the cice state
!          to be consistent with the states from assimilation
!
! method: Read in restart (restart with prior) and out restart (restart
!         with posterior) written by DART after filter.
!
! author: C Bitz June 2016
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             file_exist, error_handler, E_ERR, E_MSG, to_upper
use  netcdf_utilities_mod, only : nc_check
use netcdf


implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = &
   "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"

!------------------------------------------------------------------

character(len=256) :: dart_to_cice_input_file = 'dart_restart.nc'
character(len=256) :: original_cice_input_file = 'cice_restart.nc'
character(len=256) :: previous_cice_input_file = 'pre_restart.nc'
character(len=128) :: balance_method = 'simple_squeeze'
character(len=15)  :: r_snw_name  = 'r_snw'
integer :: gridpt_oi = 3

namelist /dart_to_cice_nml/ dart_to_cice_input_file, &
                            original_cice_input_file, &
                            previous_cice_input_file, &
                            balance_method, &
                            r_snw_name, &
                            gridpt_oi

character(len=512) :: string1, string2, msgstring
character(len=15)  :: varname
character(len=128) :: method

integer :: Nx
integer :: Ncat   ! number of categories in ice-thickness dist
integer, parameter :: Nilyr = 8   ! number of layers in ice, hardwired
integer, parameter :: Nslyr = 3   ! number of layers in snow, hardwired

real(r8), allocatable :: aicen_original(:)
real(r8), allocatable :: vicen_original(:)
real(r8), allocatable :: vsnon_original(:)
!real(r8), allocatable :: aice_original(:,:)
!real(r8), allocatable :: hicen_original(:)
!real(r8), allocatable :: hsnon_original(:)
logical :: sst_present = .true.
logical :: sst_org_present = .true.

real(r8) :: sst,sst_original
real(r8), allocatable :: aicen(:)
real(r8), allocatable :: vicen(:)
real(r8), allocatable :: vsnon(:)
real(r8), allocatable :: Tsfcn(:)
real(r8), allocatable :: qice(:,:)
real(r8), allocatable :: sice(:,:)
real(r8), allocatable :: qsno(:,:)

character (len=3) :: nchar
integer :: iunit,io,ncid,dimid,l,n,VarID
real(r8) :: aice,aice_temp
real(r8) :: vice,vice_temp
real(r8) :: vsno,vsno_temp
real(r8), parameter :: Tsmelt = 0._r8
real(r8), parameter :: c1  = 1.0_r8
real(r8), parameter :: &
     phi_init = 0.75_r8, &
     dSin0_frazil = 3.0_r8
real(r8), parameter :: sss = 34.7_r8
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(r8) :: squeeze,cc1,cc2,cc3,x1,Si0new,Ti,qsno_hold,qi0new
real(r8), allocatable ::  hin_max(:)
real(r8), allocatable ::  hcat_midpoint(:)

call initialize_utilities(progname='dart_to_cice')

call find_namelist_in_file("input.nml", "dart_to_cice_nml", iunit)
read(iunit, nml = dart_to_cice_nml, iostat = io)
call check_namelist_read(iunit, io, "dart_to_cice_nml")

method = balance_method
call to_upper(method)

! check on namelist stuff, and whether files exist
write(string1,*) 'converting DART output file "'// &
                 &trim(dart_to_cice_input_file)//'" to one CICE will like'
write(string2,*) 'using the "'//trim(balance_method)//'" method.'
call error_handler(E_MSG,'dart_to_cice',string1,text2=string2)

if ( .not. file_exist(dart_to_cice_input_file) ) then
   write(string1,*) 'cannot open "', trim(dart_to_cice_input_file),'" for updating.'
   call error_handler(E_ERR,'dart_to_cice:filename not found ',trim(dart_to_cice_input_file))
endif

if ( .not. file_exist(original_cice_input_file) ) then
   write(string1,*) 'cannot open "', trim(original_cice_input_file),'" for reading.'
   call error_handler(E_ERR,'dart_to_cice:filename not found ',trim(original_cice_input_file))
endif


call nc_check( nf90_open(trim(original_cice_input_file), NF90_NOWRITE, ncid), &
                  'dart_to_cice', 'open "'//trim(original_cice_input_file)//'"')

call nc_check(nf90_inq_dimid(ncid,"ncat",dimid), &
                  'dart_to_cice', 'inquire ncat dimid from "'//trim(original_cice_input_file)//'"')
call nc_check(nf90_inquire_dimension(ncid,dimid,len=Ncat), &
                   'dart_to_cice', 'inquire ncat from "'//trim(original_cice_input_file)//'"')
call nc_check(nf90_inq_dimid(ncid,"ni",dimid), &
                   'dart_to_cice', 'inquire ni dimid from "'//trim(original_cice_input_file)//'"')
call nc_check(nf90_inquire_dimension(ncid,dimid,len=Nx),&
                   'dart_to_cice', 'inquire ni from "'//trim(original_cice_input_file)//'"')

allocate(aicen_original(NCAT),vicen_original(NCAT),vsnon_original(NCAT),Tsfcn(NCAT),qice(Nilyr,NCAT),sice(Nilyr,NCAT),qsno(Nslyr,NCAT))
call get_variable(ncid,'aicen',aicen_original,original_cice_input_file,gridpt_oi,Ncat)
call get_variable(ncid,'vicen',vicen_original,original_cice_input_file,gridpt_oi,Ncat)
call get_variable(ncid,'vsnon',vsnon_original,original_cice_input_file,gridpt_oi,Ncat)
call get_variable(ncid,'Tsfcn',Tsfcn,dart_to_cice_input_file,gridpt_oi,Ncat)
call get_variable1d(ncid,'sst',sst_original,dart_to_cice_input_file,gridpt_oi,sst_org_present)
do l=1, Nilyr
  write(nchar,'(i3.3)') l
  call get_variable(ncid,'qice'//trim(nchar),qice(l,:),dart_to_cice_input_file,gridpt_oi,Ncat)
  call get_variable(ncid,'sice'//trim(nchar),sice(l,:),dart_to_cice_input_file,gridpt_oi,Ncat)
enddo
do l=1, Nslyr
  write(nchar,'(i3.3)') l
  call get_variable(ncid,'qsno'//trim(nchar),qsno(l,:),dart_to_cice_input_file,gridpt_oi,Ncat)
enddo
call nc_check(nf90_close(ncid),'dart_to_cice', 'close '//trim(original_cice_input_file))
!!!!!!!!!
call nc_check( nf90_open(trim(dart_to_cice_input_file), NF90_NOWRITE, ncid), &
                  'dart_to_cice', 'open "'//trim(dart_to_cice_input_file)//'"')
allocate(aicen(NCAT),vicen(NCAT),vsnon(NCAT))
call get_variable(ncid,'aicen',aicen,dart_to_cice_input_file,gridpt_oi,Ncat)
call get_variable(ncid,'vicen',vicen,dart_to_cice_input_file,gridpt_oi,Ncat)
call get_variable(ncid,'vsnon',vsnon,dart_to_cice_input_file,gridpt_oi,Ncat)
call get_variable1d(ncid,'sst',sst,dart_to_cice_input_file,gridpt_oi,sst_present)
call nc_check(nf90_close(ncid),'dart_to_cice', 'close '//trim(dart_to_cice_input_file))
!!!!!!!!!!!!!!!!!!!!!!!!!
qice = min(0.0_r8,qice)
sice = max(0.0_r8,sice)
qsno = min(0.0_r8,qsno)
aicen = min(1.0_r8,aicen)
Tsfcn = min(Tsmelt,Tsfcn)
!!!!!!
aice = sum(aicen)
vice = sum(vicen)
vsno = sum(vsnon)
!!!!!!
aicen = max(0.0_r8,aicen)
vicen = max(0.0_r8,vicen)
vsnon = max(0.0_r8,vsnon)
!!!!!
aice_temp = sum(aicen)
vice_temp = sum(vicen)
vsno_temp = sum(vsnon)
!!!!!
if (aice<0.0_r8) then
  aicen(:) = 0.0_r8
  vicen(:) = 0.0_r8
  vsnon(:) = 0.0_r8
endif
!!!!!
do n=1,NCAT
  if (aice_temp > 0._r8 .and. aice>0._r8) then
      aicen(n) = aicen(n) - (aice_temp-aice)*aicen(n)/aice_temp
    endif
    if (vice_temp > 0._r8 .and. vice>0._r8) then
      vicen(n) = vicen(n) - (vice_temp-vice)*vicen(n)/vice_temp
    endif
    if (vsno_temp > 0._r8 .and. vsno > 0._r8) then
      vsnon(n) = vsnon(n) - (vsno_temp-vsno)*vsnon(n)/vsno_temp
    endif
enddo
!!!!
if (aice>1.0_r8) then
  squeeze = 1.0_r8/aice
  aicen(:) = aicen(:)*squeeze
endif
!!!!!!
if (sst_present) then
  if (aice == 0.0_r8) sst = 0.0_r8
endif
where(aicen==-999) aicen = 0.0_r8
!!!!!!
cc1 = 3._r8/real(Ncat,kind=r8)
cc2 = 15.0_r8*cc1
cc3 = 3._r8
allocate( hin_max(0:Ncat) )
allocate( hcat_midpoint(Ncat) )
hin_max(0) = 0._r8
do n = 1, NCAT
  x1 = real(n-1,kind=r8) / real(Ncat,kind=r8)
  hin_max(n) = hin_max(n-1) &
             + cc1 + cc2*(c1 + tanh(cc3*(x1-c1)))
  hcat_midpoint(n)=0.5_r8*(hin_max(n-1)+hin_max(n))
enddo
!!!!!!!
do n=1,NCAT
  if (aicen(n) > 0.0_r8 .and. aicen_original(n) > 0.0_r8) then
    if (vicen(n) == 0.0_r8) then
      vicen(n) = aicen(n)*hcat_midpoint(n)
    endif
  endif
  if (aicen(n) == 0.0_r8 .and. aicen_original(n) > 0.0_r8) then
    vicen(n) = 0.0_r8
    qice(:,n) = 0.0_r8
    sice(:,n) = 0.0_r8
    qsno(:,n) = 0.0_r8
    vsnon(n) = 0.0_r8
    Tsfcn(n) = -1.8_r8
  else if (aicen(n)>0.0_r8 .and. aicen_original(n) == 0.0_r8) then
    if (vicen(n) == 0.0_r8) vicen(n) =  aicen(n) * hcat_midpoint(n)
    Si0new = sss - dSin0_frazil
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
  endif
  if (aicen(n) == 0.0_r8) then
    vicen(n) = 0.0_r8
    vsnon(n) = 0.0_r8
  endif
enddo
!!!!!!!!
call nc_check( nf90_open(trim(original_cice_input_file), NF90_WRITE, ncid), &
                  'dart_to_cice', 'open "'//trim(original_cice_input_file)//'"')
varname='aicen'
io = nf90_inq_varid(ncid, trim(varname), VarID)
call nc_check(io, 'dart_to_cice', &
                 'inq_varid '//trim(varname)//' '//trim(original_cice_input_file))
io = nf90_put_var(ncid, VarID, aicen,start=(/gridpt_oi,1/),count=(/1,NCAT/))
call nc_check(io, 'dart_to_cice', &
                 'put_var '//trim(varname)//' '//trim(original_cice_input_file))
!!!!
varname='vicen'
io = nf90_inq_varid(ncid, trim(varname), VarID)
call nc_check(io, 'dart_to_cice', &
                 'inq_varid '//trim(varname)//' '//trim(original_cice_input_file))
io = nf90_put_var(ncid, VarID, vicen,start=(/gridpt_oi,1/),count=(/1,NCAT/))
call nc_check(io, 'dart_to_cice', &
                 'put_var '//trim(varname)//' '//trim(original_cice_input_file))
!!!!
varname='vsnon'
io = nf90_inq_varid(ncid, trim(varname), VarID)
call nc_check(io, 'dart_to_cice', &
                 'inq_varid '//trim(varname)//' '//trim(original_cice_input_file))
io = nf90_put_var(ncid, VarID, vsnon,start=(/gridpt_oi,1/),count=(/1,NCAT/))
call nc_check(io, 'dart_to_cice', &
                 'put_var '//trim(varname)//' '//trim(original_cice_input_file))
!!!!
varname='Tsfcn'
io = nf90_inq_varid(ncid, trim(varname), VarID)
call nc_check(io, 'dart_to_cice', &
                 'inq_varid '//trim(varname)//' '//trim(original_cice_input_file))
io = nf90_put_var(ncid, VarID, Tsfcn,start=(/gridpt_oi,1/),count=(/1,NCAT/))
call nc_check(io, 'dart_to_cice', &
                 'put_var '//trim(varname)//' '//trim(original_cice_input_file))
!!!!!
if (sst_present) then
  varname='sst'
  io = nf90_inq_varid(ncid, trim(varname), VarID)
  call nc_check(io, 'dart_to_cice', &
                 'inq_varid '//trim(varname)//' '//trim(original_cice_input_file))
  io = nf90_put_var(ncid, VarID, sst,start=(/gridpt_oi/))!,count=(/1/))
  call nc_check(io, 'dart_to_cice', &
                 'put_var '//trim(varname)//' '//trim(original_cice_input_file))
endif
!!!!!
do l=1, Nilyr
  write(nchar,'(i3.3)') l
  varname='qice'//trim(nchar)
  io = nf90_inq_varid(ncid, trim(varname), VarID)
  call nc_check(io, 'dart_to_cice', &
                 'inq_varid '//trim(varname)//' '//trim(original_cice_input_file))
  io = nf90_put_var(ncid, VarID, qice(l,:),start=(/gridpt_oi,1/),count=(/1,NCAT/))
  call nc_check(io, 'dart_to_cice', &
                 'put_var '//trim(varname)//' '//trim(original_cice_input_file))
  !!!!!!!!!!
  varname='sice'//trim(nchar)
  io = nf90_inq_varid(ncid, trim(varname), VarID)
  call nc_check(io, 'dart_to_cice', &
                 'inq_varid '//trim(varname)//' '//trim(original_cice_input_file))
  io = nf90_put_var(ncid, VarID, sice(l,:),start=(/gridpt_oi,1/),count=(/1,NCAT/))
  call nc_check(io, 'dart_to_cice', &
                 'put_var '//trim(varname)//' '//trim(original_cice_input_file))
enddo
!!!!
do l=1, Nslyr
  write(nchar,'(i3.3)') l
  varname='qsno'//trim(nchar)
  io = nf90_inq_varid(ncid, trim(varname), VarID)
  call nc_check(io, 'dart_to_cice', &
                 'inq_varid '//trim(varname)//' '//trim(original_cice_input_file))
  io = nf90_put_var(ncid, VarID, qsno(l,:),start=(/gridpt_oi,1/),count=(/1,NCAT/))
  call nc_check(io, 'dart_to_cice', &
                 'put_var '//trim(varname)//' '//trim(original_cice_input_file))
enddo

call nc_check(nf90_close(ncid),'dart_to_cice', 'close '//trim(original_cice_input_file))

deallocate( aicen, vicen, vsnon, Tsfcn)
deallocate( qice, sice, qsno )


call finalize_utilities('dart_to_cice')


contains

subroutine get_variable(ncid,varname,var,filename,space_index,ncat)
integer,               intent(in)  :: ncid,ncat
character(len=*),      intent(in)  :: varname
real(r8),              intent(out) :: var(ncat)
character(len=*),      intent(in)  :: filename
integer,               intent(in)  :: space_index

integer :: VarID, ndims, dimIDs
real(r8) :: holder(4,ncat)

write(6,*) 'Getting data for ',trim(varname)

io = nf90_inq_varid(ncid, trim(varname), VarID)
call nc_check(io, 'dart_to_cice', 'inq_varid '//trim(msgstring))

call nc_check(nf90_get_var(ncid, VarID, holder), 'dart_to_cice', &
         'get_var '//trim(msgstring))


var(:) = holder(gridpt_oi,:)

end subroutine get_variable
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_variable1d(ncid,varname,var,filename,space_index,var_present)
integer,               intent(in)  :: ncid
character(len=*),      intent(in)  :: varname
real(r8),              intent(out) :: var
character(len=*),      intent(in)  :: filename
integer,               intent(in)  :: space_index
logical,               intent(inout) :: var_present

integer :: VarID, ndims, dimIDs
real(r8) :: holder(4)

write(6,*) 'Getting data for ',trim(varname)

io = nf90_inq_varid(ncid, trim(varname), VarID)
if(io /= nf90_NoErr) then
  write(6,*) "No netcdf ID for ",trim(varname)
  var_present = .false.
  return
endif
call nc_check(io, 'dart_to_cice', 'inq_varid '//trim(msgstring))

call nc_check(nf90_get_var(ncid, VarID, holder), 'dart_to_cice', &
         'get_var '//trim(msgstring))


var = holder(gridpt_oi)

end subroutine get_variable1d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function snow_enthaply(Ti) result(qsno)
    real(r8), intent(in) :: Ti

    real(r8),parameter :: rhos = 330.0_r8, &
                        Lfresh = 2.835e6_r8 - 2.501e6_r8, &
                        cp_ice = 2106._r8
    real(r8) :: qsno

    qsno = -rhos*(Lfresh - cp_ice*min(0.0_r8,Ti))
  end function snow_enthaply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end program dart_to_cice

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
