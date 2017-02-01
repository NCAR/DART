! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

MODULE wrf_data_module

use     types_mod, only : r8
use utilities_mod, only : register_module, error_handler, E_ERR, E_MSG

use netcdf

implicit none
private

public :: wrf_data, wrf_bdy_data, wrf_open_and_alloc, wrfbdy_open_and_alloc, &
          wrf_dealloc, wrfbdy_dealloc, wrf_io, wrfbdy_io, set_wrf_date, get_wrf_date

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

TYPE wrf_data

   integer :: ncid                                    ! netcdf id for file
   integer :: sls_id, sls, bt_id, bt, sn_id, sn, we_id, we
   integer :: u_id, v_id, w_id, ph_id, phb_id, t_id, mu_id, mub_id, &
              tsk_id, qv_id, qc_id, qr_id, qi_id, qs_id, qg_id, qnice_id, &
              u10_id, v10_id, t2_id, th2_id, q2_id, ps_id, &
              tslb_id, smois_id, sh2o_id, hdiab_id

   integer :: n_moist
   logical :: surf_obs
   logical :: soil_data
   logical :: h_diab


!---
!  arrays for data

   real(r8), pointer :: u(:,:,:)
   real(r8), pointer :: v(:,:,:)
   real(r8), pointer :: w(:,:,:)
   real(r8), pointer :: ph(:,:,:)
   real(r8), pointer :: phb(:,:,:)
   real(r8), pointer :: t(:,:,:)
   real(r8), pointer :: tsk(:,:)
   real(r8), pointer :: qv(:,:,:)
   real(r8), pointer :: qc(:,:,:)
   real(r8), pointer :: qr(:,:,:)
   real(r8), pointer :: qi(:,:,:)
   real(r8), pointer :: qs(:,:,:)
   real(r8), pointer :: qg(:,:,:)
   real(r8), pointer :: qnice(:,:,:)
   real(r8), pointer :: mu(:,:)
   real(r8), pointer :: mub(:,:)
   real(r8), pointer :: u10(:,:)
   real(r8), pointer :: v10(:,:)
   real(r8), pointer :: t2(:,:)
   real(r8), pointer :: th2(:,:)
   real(r8), pointer :: q2(:,:)
   real(r8), pointer :: ps(:,:)
   real(r8), pointer :: tslb(:,:,:)
   real(r8), pointer :: smois(:,:,:)
   real(r8), pointer :: sh2o(:,:,:)
   real(r8), pointer :: hdiab(:,:,:)

end type

TYPE wrf_bdy_data

   integer :: ncid                                       ! netcdf id for file
   integer :: time_id, time, bdywdth_id, bdywdth,      &
              bt_id, bt, sn_id, sn, we_id, we  
   integer :: uxs_id , uxe_id , uys_id , uye_id ,      &
              utxs_id, utxe_id, utys_id, utye_id,      &
              vxs_id , vxe_id , vys_id , vye_id ,      &
              vtxs_id, vtxe_id, vtys_id, vtye_id,      &
              wxs_id , wxe_id , wys_id , wye_id ,      &
              wtxs_id, wtxe_id, wtys_id, wtye_id,      &
              phxs_id , phxe_id , phys_id , phye_id ,  &
              phtxs_id, phtxe_id, phtys_id, phtye_id,  &
              txs_id , txe_id , tys_id , tye_id ,      &
              ttxs_id, ttxe_id, ttys_id, ttye_id,      &
              muxs_id , muxe_id , muys_id , muye_id ,  &
              mutxs_id, mutxe_id, mutys_id, mutye_id,  &
              qvxs_id , qvxe_id , qvys_id , qvye_id ,  &
              qvtxs_id, qvtxe_id, qvtys_id, qvtye_id,  &
              qcxs_id , qcxe_id , qcys_id , qcye_id ,  &
              qctxs_id, qctxe_id, qctys_id, qctye_id,  &
              qrxs_id , qrxe_id , qrys_id , qrye_id ,  &
              qrtxs_id, qrtxe_id, qrtys_id, qrtye_id,  &
              qixs_id , qixe_id , qiys_id , qiye_id ,  &
              qitxs_id, qitxe_id, qitys_id, qitye_id,  &
              qsxs_id , qsxe_id , qsys_id , qsye_id ,  &
              qstxs_id, qstxe_id, qstys_id, qstye_id,  &
              qgxs_id , qgxe_id , qgys_id , qgye_id ,  &
              qgtxs_id, qgtxe_id, qgtys_id, qgtye_id,  &
              qnicexs_id , qnicexe_id , qniceys_id , qniceye_id ,  &
              qnicetxs_id, qnicetxe_id, qnicetys_id, qnicetye_id

   integer :: n_moist

   real(r8), pointer :: uxs(:,:,:,:) , uxe(:,:,:,:) , uys(:,:,:,:) , uye(:,:,:,:)
   real(r8), pointer :: utxs(:,:,:,:), utxe(:,:,:,:), utys(:,:,:,:), utye(:,:,:,:)
   real(r8), pointer :: vxs(:,:,:,:) , vxe(:,:,:,:) , vys(:,:,:,:) , vye(:,:,:,:)
   real(r8), pointer :: vtxs(:,:,:,:), vtxe(:,:,:,:), vtys(:,:,:,:), vtye(:,:,:,:)
   real(r8), pointer :: wxs(:,:,:,:) , wxe(:,:,:,:) , wys(:,:,:,:) , wye(:,:,:,:)
   real(r8), pointer :: wtxs(:,:,:,:), wtxe(:,:,:,:), wtys(:,:,:,:), wtye(:,:,:,:)
   real(r8), pointer :: phxs(:,:,:,:) , phxe(:,:,:,:) , phys(:,:,:,:) , phye(:,:,:,:)
   real(r8), pointer :: phtxs(:,:,:,:), phtxe(:,:,:,:), phtys(:,:,:,:),phtye(:,:,:,:)
   real(r8), pointer :: txs(:,:,:,:) , txe(:,:,:,:) , tys(:,:,:,:) , tye(:,:,:,:)
   real(r8), pointer :: ttxs(:,:,:,:), ttxe(:,:,:,:), ttys(:,:,:,:), ttye(:,:,:,:)
   real(r8), pointer :: muxs(:,:,:) , muxe(:,:,:) , muys(:,:,:) , muye(:,:,:)
   real(r8), pointer :: mutxs(:,:,:), mutxe(:,:,:), mutys(:,:,:),mutye(:,:,:)
   real(r8), pointer :: qvxs(:,:,:,:) , qvxe(:,:,:,:) , qvys(:,:,:,:) , qvye(:,:,:,:)
   real(r8), pointer :: qvtxs(:,:,:,:), qvtxe(:,:,:,:), qvtys(:,:,:,:),qvtye(:,:,:,:)
   real(r8), pointer :: qcxs(:,:,:,:) , qcxe(:,:,:,:) , qcys(:,:,:,:) , qcye(:,:,:,:)
   real(r8), pointer :: qctxs(:,:,:,:), qctxe(:,:,:,:), qctys(:,:,:,:),qctye(:,:,:,:)
   real(r8), pointer :: qrxs(:,:,:,:) , qrxe(:,:,:,:) , qrys(:,:,:,:) , qrye(:,:,:,:)
   real(r8), pointer :: qrtxs(:,:,:,:), qrtxe(:,:,:,:), qrtys(:,:,:,:),qrtye(:,:,:,:)
   real(r8), pointer :: qixs(:,:,:,:) , qixe(:,:,:,:) , qiys(:,:,:,:) , qiye(:,:,:,:)
   real(r8), pointer :: qitxs(:,:,:,:), qitxe(:,:,:,:), qitys(:,:,:,:),qitye(:,:,:,:)
   real(r8), pointer :: qsxs(:,:,:,:) , qsxe(:,:,:,:) , qsys(:,:,:,:) , qsye(:,:,:,:)
   real(r8), pointer :: qstxs(:,:,:,:), qstxe(:,:,:,:), qstys(:,:,:,:),qstye(:,:,:,:)
   real(r8), pointer :: qgxs(:,:,:,:) , qgxe(:,:,:,:) , qgys(:,:,:,:) , qgye(:,:,:,:)
   real(r8), pointer :: qgtxs(:,:,:,:), qgtxe(:,:,:,:), qgtys(:,:,:,:),qgtye(:,:,:,:)
   real(r8), pointer :: qnicexs(:,:,:,:) , qnicexe(:,:,:,:) , qniceys(:,:,:,:) , qniceye(:,:,:,:)
   real(r8), pointer :: qnicetxs(:,:,:,:), qnicetxe(:,:,:,:), qnicetys(:,:,:,:),qnicetye(:,:,:,:)

end type

logical, save :: module_initialized = .false.

contains

subroutine initialize_module

  call register_module(source, revision, revdate)
  module_initialized = .true.

end subroutine initialize_module


!**********************************************************************

subroutine wrf_open_and_alloc( wrf, file_name, mode, debug )

implicit none

type(wrf_data)     :: wrf
character (len=*)  :: file_name   ! filename from which dimensions, 
                                  ! variable id's are read
integer            :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer            :: mode, istatus
logical            :: debug

character (len=80) :: name

if ( .not. module_initialized ) call initialize_module

if(debug) write(6,*) 'Opening ',file_name
call check ( nf90_open(file_name, mode, wrf%ncid) )
if(debug) write(6,*) ' wrf%ncid is ',wrf%ncid

! get wrf grid dimensions

call check ( nf90_inq_dimid(wrf%ncid, "bottom_top", wrf%bt_id) )
call check ( nf90_inquire_dimension(wrf%ncid, wrf%bt_id, name, wrf%bt) )

call check ( nf90_inq_dimid(wrf%ncid, "south_north", wrf%sn_id) )
call check ( nf90_inquire_dimension(wrf%ncid, wrf%sn_id, name, wrf%sn) )

call check ( nf90_inq_dimid(wrf%ncid, "west_east", wrf%we_id) )
call check ( nf90_inquire_dimension(wrf%ncid, wrf%we_id, name, wrf%we) )

call check ( nf90_inq_dimid(wrf%ncid, "soil_layers_stag", wrf%sls_id) )
call check ( nf90_inquire_dimension(wrf%ncid, wrf%sls_id, name, wrf%sls) )

if(debug) write(6,*) ' dimensions bt, sn, we are ',wrf%bt,wrf%sn,wrf%we

!---
! get wrf variable ids and allocate space for wrf variables

call check ( nf90_inq_varid(wrf%ncid, "U", wrf%u_id))
if(debug) write(6,*) ' u_id = ',wrf%u_id
allocate(wrf%u(wrf%we+1,wrf%sn,wrf%bt))

call check ( nf90_inq_varid(wrf%ncid, "V", wrf%v_id))
if(debug) write(6,*) ' v_id = ',wrf%v_id
allocate(wrf%v(wrf%we,wrf%sn+1,wrf%bt))

call check ( nf90_inq_varid(wrf%ncid, "W", wrf%w_id))
if(debug) write(6,*) ' w_id = ',wrf%w_id
allocate(wrf%w(wrf%we,wrf%sn,wrf%bt+1))

call check ( nf90_inq_varid(wrf%ncid, "PH", wrf%ph_id))
if(debug) write(6,*) ' ph_id = ',wrf%ph_id
allocate(wrf%ph(wrf%we,wrf%sn,wrf%bt+1))

call check ( nf90_inq_varid(wrf%ncid, "PHB", wrf%phb_id))
if(debug) write(6,*) ' phb_id = ',wrf%phb_id
allocate(wrf%phb(wrf%we,wrf%sn,wrf%bt+1))

call check ( nf90_inq_varid(wrf%ncid, "T", wrf%t_id))
if(debug) write(6,*) ' t_id = ',wrf%t_id
allocate(wrf%t(wrf%we,wrf%sn,wrf%bt))

call check ( nf90_inq_varid(wrf%ncid, "MU", wrf%mu_id))
if(debug) write(6,*) ' mu_id = ',wrf%mu_id
allocate(wrf%mu(wrf%we,wrf%sn))

call check ( nf90_inq_varid(wrf%ncid, "MUB", wrf%mub_id))
if(debug) write(6,*) ' mub_id = ',wrf%mub_id
allocate(wrf%mub(wrf%we,wrf%sn))

if(wrf%n_moist > 0) then
   call check ( nf90_inq_varid(wrf%ncid, "QVAPOR", wrf%qv_id))
   allocate(wrf%qv(wrf%we,wrf%sn,wrf%bt))
endif

if(wrf%n_moist > 1) then
   call check ( nf90_inq_varid(wrf%ncid, "QCLOUD", wrf%qc_id))
   allocate(wrf%qc(wrf%we,wrf%sn,wrf%bt))
endif

if(wrf%n_moist > 2) then
   call check ( nf90_inq_varid(wrf%ncid, "QRAIN", wrf%qr_id))
   allocate(wrf%qr(wrf%we,wrf%sn,wrf%bt))
endif

if(wrf%n_moist > 3) then
   call check ( nf90_inq_varid(wrf%ncid, "QICE", wrf%qi_id))
   allocate(wrf%qi(wrf%we,wrf%sn,wrf%bt))
endif

if(wrf%n_moist > 4) then
   call check ( nf90_inq_varid(wrf%ncid, "QSNOW", wrf%qs_id))
   allocate(wrf%qs(wrf%we,wrf%sn,wrf%bt))
endif

if(wrf%n_moist > 5) then
   call check ( nf90_inq_varid(wrf%ncid, "QGRAUP", wrf%qg_id))
   allocate(wrf%qg(wrf%we,wrf%sn,wrf%bt))
endif

if(wrf%n_moist > 6) then
   call check ( nf90_inq_varid(wrf%ncid, "QNICE", wrf%qnice_id))
   allocate(wrf%qnice(wrf%we,wrf%sn,wrf%bt))
endif

if(wrf%n_moist > 7) then
   write(6,*) 'n_moist = ',wrf%n_moist
   call error_handler(E_ERR,'wrf_open_and_alloc', &
         'n_moist is too large.', source, revision, revdate)
endif

if( wrf%surf_obs ) then

   call check ( nf90_inq_varid(wrf%ncid, "U10", wrf%u10_id))
   allocate(wrf%u10(wrf%we,wrf%sn))

   call check ( nf90_inq_varid(wrf%ncid, "V10", wrf%v10_id))
   allocate(wrf%v10(wrf%we,wrf%sn))

   call check ( nf90_inq_varid(wrf%ncid, "T2", wrf%t2_id))
   allocate(wrf%t2(wrf%we,wrf%sn))

   call check ( nf90_inq_varid(wrf%ncid, "TH2", wrf%th2_id))
   allocate(wrf%th2(wrf%we,wrf%sn))

   call check ( nf90_inq_varid(wrf%ncid, "Q2", wrf%q2_id))
   allocate(wrf%q2(wrf%we,wrf%sn))

   allocate(wrf%ps(wrf%we,wrf%sn))
   istatus = nf90_inq_varid(wrf%ncid, "PSFC", wrf%ps_id)
   if(istatus /= nf90_noerr) then
      call error_handler(E_MSG,'PSFC', &
           trim(nf90_strerror(istatus)), source, revision, revdate)
      if(mode == NF90_WRITE) then
         call error_handler(E_MSG,'wrf_open_and_alloc', &
              'creates PSFC', source, revision, revdate)
         call check(nf90_Inquire(wrf%ncid, nDimensions, nVariables, nAttributes, unlimitedDimID))
         call check(nf90_Redef(wrf%ncid))
         call check(nf90_def_var(wrf%ncid, name="PSFC", xtype=nf90_real, &
              dimids= (/ wrf%we_id,  wrf%sn_id, unlimitedDimID/), varid=wrf%ps_id) )
         call check(nf90_enddef(wrf%ncid))
         wrf%ps(:,:) = 0.0_r8
         call check( nf90_put_var(wrf%ncid, wrf%ps_id, wrf%ps, start = (/ 1, 1, 1 /)))
      endif
   endif

endif

if( wrf%soil_data ) then

   call check ( nf90_inq_varid(wrf%ncid, "TSLB", wrf%tslb_id))
   if(debug) write(6,*) ' tslb_id = ',wrf%tslb_id
   allocate(wrf%tslb(wrf%we,wrf%sn,wrf%sls))

   call check ( nf90_inq_varid(wrf%ncid, "SMOIS", wrf%smois_id))
   if(debug) write(6,*) ' smois_id = ',wrf%smois_id
   allocate(wrf%smois(wrf%we,wrf%sn,wrf%sls))

   call check ( nf90_inq_varid(wrf%ncid, "SH2O", wrf%sh2o_id))
   if(debug) write(6,*) ' sh2o_id = ',wrf%sh2o_id
   allocate(wrf%sh2o(wrf%we,wrf%sn,wrf%sls))

   call check ( nf90_inq_varid(wrf%ncid, "TSK", wrf%tsk_id))
   if(debug) write(6,*) ' tsk_id = ',wrf%tsk_id
   allocate(wrf%tsk(wrf%we,wrf%sn))

endif

if( wrf%h_diab ) then

   allocate(wrf%hdiab(wrf%we,wrf%sn,wrf%bt))
   istatus = nf90_inq_varid(wrf%ncid, "H_DIABATIC", wrf%hdiab_id)
   if(istatus /= nf90_noerr) then
      call error_handler(E_MSG,'H_DIABATIC', &
           trim(nf90_strerror(istatus)), source, revision, revdate)
      if(mode == NF90_WRITE) then
         call error_handler(E_MSG,'wrf_open_and_alloc', &
              'creates H_DIABATIC', source, revision, revdate)
         call check(nf90_Inquire(wrf%ncid, nDimensions, nVariables, nAttributes, unlimitedDimID))
         call check(nf90_Redef(wrf%ncid))
         call check(nf90_def_var(wrf%ncid, name="H_DIABATIC", xtype=nf90_real, &
              dimids= (/ wrf%we_id,  wrf%sn_id, wrf%bt_id, unlimitedDimID/), varid=wrf%hdiab_id) )
         call check(nf90_put_att(wrf%ncid, wrf%hdiab_id, "FieldType", 104))
         call check(nf90_put_att(wrf%ncid, wrf%hdiab_id, "MemoryOrder", "XYZ"))
         call check(nf90_put_att(wrf%ncid, wrf%hdiab_id, "description", &
              "PREVIOUS TIMESTEP CONDENSATIONAL HEATING"))
         call check(nf90_put_att(wrf%ncid, wrf%hdiab_id, "units", ""))
         call check(nf90_put_att(wrf%ncid, wrf%hdiab_id, "stagger", ""))
         call check(nf90_enddef(wrf%ncid))
         wrf%hdiab(:,:,:) = 0.0_r8
         call check( nf90_put_var(wrf%ncid, wrf%hdiab_id, wrf%hdiab, start = (/ 1, 1, 1, 1 /)))
      endif
   endif

endif

contains

  ! Internal subroutine - checks error status after each netcdf, prints 
  !                       text message each time an error code is returned. 
  subroutine check(istatus)
    integer, intent ( in) :: istatus 
    if(istatus /= nf90_noerr) call error_handler(E_ERR,'wrf_open_and_alloc', &
         trim(nf90_strerror(istatus)), source, revision, revdate)
  end subroutine check

end subroutine wrf_open_and_alloc

!**********************************************************************

subroutine wrf_dealloc( wrf )

implicit none

type(wrf_data) :: wrf

if ( .not. module_initialized ) call initialize_module

deallocate(wrf%u)
deallocate(wrf%v)
deallocate(wrf%w)
deallocate(wrf%ph)
deallocate(wrf%phb)
deallocate(wrf%t)
deallocate(wrf%mu)
deallocate(wrf%mub)
if(wrf%n_moist > 0) then
   deallocate(wrf%qv)
endif

if(wrf%n_moist > 1) then
   deallocate(wrf%qc)
endif

if(wrf%n_moist > 2) then
   deallocate(wrf%qr)
endif

if(wrf%n_moist > 3) then
   deallocate(wrf%qi)
endif

if(wrf%n_moist > 4) then
   deallocate(wrf%qs)
endif

if(wrf%n_moist > 5) then
   deallocate(wrf%qg)
endif

if(wrf%n_moist > 6) then
   deallocate(wrf%qnice)
endif

if(wrf%n_moist > 7) then
   write(6,*) 'n_moist = ',wrf%n_moist
   call error_handler(E_ERR,'wrf_dealloc', &
         'n_moist is too large.', source, revision, revdate)
endif

if( wrf%surf_obs ) then
   deallocate(wrf%u10)
   deallocate(wrf%v10)
   deallocate(wrf%t2)
   deallocate(wrf%th2)
   deallocate(wrf%q2)
   deallocate(wrf%ps)
endif

if( wrf%soil_data ) then
  deallocate(wrf%tslb)
  deallocate(wrf%smois)
  deallocate(wrf%sh2o)
  deallocate(wrf%tsk)
endif

if( wrf%h_diab ) then
   deallocate(wrf%hdiab)
endif

end subroutine wrf_dealloc

!---------------------------------------------------------------

subroutine wrfbdy_open_and_alloc( wrfbdy, file_name, mode, debug )

implicit none

type(wrf_bdy_data) :: wrfbdy
character (len=*)  :: file_name   ! filename from which dimensions, 
                                  ! variable id's are read
integer            :: mode
logical            :: debug

character (len=80) :: name

if ( .not. module_initialized ) call initialize_module

if (debug) write(6,*) '   in wrfbdy_open_and_alloc, file =', file_name

call check ( nf90_open(file_name, mode, wrfbdy%ncid) ) 
if(debug) write(6,*) ' wrfbdy%ncid is ',wrfbdy%ncid

! get wrfbdy dimensions

call check ( nf90_inq_dimid(wrfbdy%ncid, "Time", wrfbdy%time_id) )
call check ( nf90_inquire_dimension(wrfbdy%ncid, wrfbdy%time_id, name, wrfbdy%time) )

call check ( nf90_inq_dimid(wrfbdy%ncid, "bdy_width", wrfbdy%bdywdth_id) )
call check ( nf90_inquire_dimension(wrfbdy%ncid, wrfbdy%bdywdth_id, name, wrfbdy%bdywdth) )

call check ( nf90_inq_dimid(wrfbdy%ncid, "bottom_top", wrfbdy%bt_id) )
call check ( nf90_inquire_dimension(wrfbdy%ncid, wrfbdy%bt_id, name, wrfbdy%bt) )

call check ( nf90_inq_dimid(wrfbdy%ncid, "south_north", wrfbdy%sn_id) )
call check ( nf90_inquire_dimension(wrfbdy%ncid, wrfbdy%sn_id, name, wrfbdy%sn) )

call check ( nf90_inq_dimid(wrfbdy%ncid, "west_east", wrfbdy%we_id) )
call check ( nf90_inquire_dimension(wrfbdy%ncid, wrfbdy%we_id, name, wrfbdy%we) )

if(debug) write(6,*) ' dimensions bt, sn, we are ',wrfbdy%bt,wrfbdy%sn,wrfbdy%we

!---
! get wrfbdy variable ids and allocate space for wrfbdy variables

  !-- u on bdy
call check ( nf90_inq_varid(wrfbdy%ncid, "U_BXS", wrfbdy%uxs_id) )
if(debug) write(6,*) ' uxs_id = ',wrfbdy%uxs_id
allocate( wrfbdy%uxs( wrfbdy%sn, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

call check ( nf90_inq_varid(wrfbdy%ncid, "U_BXE", wrfbdy%uxe_id) )
if(debug) write(6,*) ' uxe_id = ',wrfbdy%uxe_id
allocate( wrfbdy%uxe( wrfbdy%sn, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

call check ( nf90_inq_varid(wrfbdy%ncid, "U_BYS", wrfbdy%uys_id) )
if(debug) write(6,*) ' uys_id = ',wrfbdy%uys_id
allocate( wrfbdy%uys( wrfbdy%we+1, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

call check ( nf90_inq_varid(wrfbdy%ncid, "U_BYE", wrfbdy%uye_id) )
if(debug) write(6,*) ' uye_id = ',wrfbdy%uye_id
allocate( wrfbdy%uye( wrfbdy%we+1, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

  !-- u tendency on bdy
call check ( nf90_inq_varid(wrfbdy%ncid, "U_BTXS", wrfbdy%utxs_id) )
if(debug) write(6,*) ' utxs_id = ',wrfbdy%utxs_id
allocate( wrfbdy%utxs( wrfbdy%sn, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

call check ( nf90_inq_varid(wrfbdy%ncid, "U_BTXE", wrfbdy%utxe_id) )
if(debug) write(6,*) ' utxe_id = ',wrfbdy%utxe_id
allocate( wrfbdy%utxe( wrfbdy%sn, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

call check ( nf90_inq_varid(wrfbdy%ncid, "U_BTYS", wrfbdy%utys_id) )
if(debug) write(6,*) ' utys_id = ',wrfbdy%utys_id
allocate( wrfbdy%utys( wrfbdy%we+1, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

call check ( nf90_inq_varid(wrfbdy%ncid, "U_BTYE", wrfbdy%utye_id) )
if(debug) write(6,*) ' utye_id = ',wrfbdy%utye_id
allocate( wrfbdy%utye( wrfbdy%we+1, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

  !-- v on bdy
call check ( nf90_inq_varid(wrfbdy%ncid, "V_BXS", wrfbdy%vxs_id) )
if(debug) write(6,*) ' vxs_id = ',wrfbdy%vxs_id
allocate( wrfbdy%vxs( wrfbdy%sn+1, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

call check ( nf90_inq_varid(wrfbdy%ncid, "V_BXE", wrfbdy%vxe_id) )
if(debug) write(6,*) ' vxe_id = ',wrfbdy%vxe_id
allocate( wrfbdy%vxe( wrfbdy%sn+1, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

call check ( nf90_inq_varid(wrfbdy%ncid, "V_BYS", wrfbdy%vys_id) )
if(debug) write(6,*) ' vys_id = ',wrfbdy%vys_id
allocate( wrfbdy%vys( wrfbdy%we, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

call check ( nf90_inq_varid(wrfbdy%ncid, "V_BYE", wrfbdy%vye_id) )
if(debug) write(6,*) ' vye_id = ',wrfbdy%vye_id
allocate( wrfbdy%vye( wrfbdy%we, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

  !-- v tendency on bdy
call check ( nf90_inq_varid(wrfbdy%ncid, "V_BTXS", wrfbdy%vtxs_id) )
if(debug) write(6,*) ' vtxs_id = ',wrfbdy%vtxs_id
allocate( wrfbdy%vtxs( wrfbdy%sn+1, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

call check ( nf90_inq_varid(wrfbdy%ncid, "V_BTXE", wrfbdy%vtxe_id) )
if(debug) write(6,*) ' vtxe_id = ',wrfbdy%vtxe_id
allocate( wrfbdy%vtxe( wrfbdy%sn+1, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

call check ( nf90_inq_varid(wrfbdy%ncid, "V_BTYS", wrfbdy%vtys_id) )
if(debug) write(6,*) ' vtys_id = ',wrfbdy%vtys_id
allocate( wrfbdy%vtys( wrfbdy%we, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

call check ( nf90_inq_varid(wrfbdy%ncid, "V_BTYE", wrfbdy%vtye_id) )
if(debug) write(6,*) ' vtye_id = ',wrfbdy%vtye_id
allocate( wrfbdy%vtye( wrfbdy%we, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

  !-- w on bdy
call check ( nf90_inq_varid(wrfbdy%ncid, "W_BXS", wrfbdy%wxs_id) )
allocate( wrfbdy%wxs( wrfbdy%sn, wrfbdy%bt+1, wrfbdy%bdywdth, wrfbdy%time ) )

call check ( nf90_inq_varid(wrfbdy%ncid, "W_BXE", wrfbdy%wxe_id) )
allocate( wrfbdy%wxe( wrfbdy%sn, wrfbdy%bt+1, wrfbdy%bdywdth, wrfbdy%time ) )

call check ( nf90_inq_varid(wrfbdy%ncid, "W_BYS", wrfbdy%wys_id) )
allocate( wrfbdy%wys( wrfbdy%we, wrfbdy%bt+1, wrfbdy%bdywdth, wrfbdy%time ) )

call check ( nf90_inq_varid(wrfbdy%ncid, "W_BYE", wrfbdy%wye_id) )
allocate( wrfbdy%wye( wrfbdy%we, wrfbdy%bt+1, wrfbdy%bdywdth, wrfbdy%time ) )

  !-- w tendency on bdy
call check ( nf90_inq_varid(wrfbdy%ncid, "W_BTXS", wrfbdy%wtxs_id) )
allocate( wrfbdy%wtxs( wrfbdy%sn, wrfbdy%bt+1, wrfbdy%bdywdth, wrfbdy%time ) )

call check ( nf90_inq_varid(wrfbdy%ncid, "W_BTXE", wrfbdy%wtxe_id) )
allocate( wrfbdy%wtxe( wrfbdy%sn, wrfbdy%bt+1, wrfbdy%bdywdth, wrfbdy%time ) )

call check ( nf90_inq_varid(wrfbdy%ncid, "W_BTYS", wrfbdy%wtys_id) )
allocate( wrfbdy%wtys( wrfbdy%we, wrfbdy%bt+1, wrfbdy%bdywdth, wrfbdy%time ) )

call check ( nf90_inq_varid(wrfbdy%ncid, "W_BTYE", wrfbdy%wtye_id) )
allocate( wrfbdy%wtye( wrfbdy%we, wrfbdy%bt+1, wrfbdy%bdywdth, wrfbdy%time ) )

  !-- height on bdy
call check ( nf90_inq_varid(wrfbdy%ncid, "PH_BXS", wrfbdy%phxs_id) )
if(debug) write(6,*) ' phxs_id = ',wrfbdy%phxs_id
allocate( wrfbdy%phxs( wrfbdy%sn, wrfbdy%bt+1, wrfbdy%bdywdth, wrfbdy%time ) )

call check ( nf90_inq_varid(wrfbdy%ncid, "PH_BXE", wrfbdy%phxe_id) )
if(debug) write(6,*) ' phxe_id = ',wrfbdy%phxe_id
allocate( wrfbdy%phxe( wrfbdy%sn, wrfbdy%bt+1, wrfbdy%bdywdth, wrfbdy%time ) )

call check ( nf90_inq_varid(wrfbdy%ncid, "PH_BYS", wrfbdy%phys_id) )
if(debug) write(6,*) ' phys_id = ',wrfbdy%phys_id
allocate( wrfbdy%phys( wrfbdy%we, wrfbdy%bt+1, wrfbdy%bdywdth, wrfbdy%time ) )

call check ( nf90_inq_varid(wrfbdy%ncid, "PH_BYE", wrfbdy%phye_id) )
if(debug) write(6,*) ' phye_id = ',wrfbdy%phye_id
allocate( wrfbdy%phye( wrfbdy%we, wrfbdy%bt+1, wrfbdy%bdywdth, wrfbdy%time ) )

  !-- height tendency on bdy
call check ( nf90_inq_varid(wrfbdy%ncid, "PH_BTXS", wrfbdy%phtxs_id) )
if(debug) write(6,*) ' phtxs_id = ',wrfbdy%phtxs_id
allocate( wrfbdy%phtxs( wrfbdy%sn, wrfbdy%bt+1, wrfbdy%bdywdth, wrfbdy%time ) )

call check ( nf90_inq_varid(wrfbdy%ncid, "PH_BTXE", wrfbdy%phtxe_id) )
if(debug) write(6,*) ' phtxe_id = ',wrfbdy%phtxe_id
allocate( wrfbdy%phtxe( wrfbdy%sn, wrfbdy%bt+1, wrfbdy%bdywdth, wrfbdy%time ) )

call check ( nf90_inq_varid(wrfbdy%ncid, "PH_BTYS", wrfbdy%phtys_id) )
if(debug) write(6,*) ' phtys_id = ',wrfbdy%phtys_id
allocate( wrfbdy%phtys( wrfbdy%we, wrfbdy%bt+1, wrfbdy%bdywdth, wrfbdy%time ) )

call check ( nf90_inq_varid(wrfbdy%ncid, "PH_BTYE", wrfbdy%phtye_id) )
if(debug) write(6,*) ' phtye_id = ',wrfbdy%phtye_id
allocate( wrfbdy%phtye( wrfbdy%we, wrfbdy%bt+1, wrfbdy%bdywdth, wrfbdy%time ) )

  !-- t on bdy
call check ( nf90_inq_varid(wrfbdy%ncid, "T_BXS", wrfbdy%txs_id) )
if(debug) write(6,*) ' txs_id = ',wrfbdy%txs_id
allocate( wrfbdy%txs( wrfbdy%sn, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

call check ( nf90_inq_varid(wrfbdy%ncid, "T_BXE", wrfbdy%txe_id) )
if(debug) write(6,*) ' txe_id = ',wrfbdy%txe_id
allocate( wrfbdy%txe( wrfbdy%sn, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

call check ( nf90_inq_varid(wrfbdy%ncid, "T_BYS", wrfbdy%tys_id) )
if(debug) write(6,*) ' tys_id = ',wrfbdy%tys_id
allocate( wrfbdy%tys( wrfbdy%we, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

call check ( nf90_inq_varid(wrfbdy%ncid, "T_BYE", wrfbdy%tye_id) )
if(debug) write(6,*) ' tye_id = ',wrfbdy%tye_id
allocate( wrfbdy%tye( wrfbdy%we, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

  !-- t tendency on bdy
call check ( nf90_inq_varid(wrfbdy%ncid, "T_BTXS", wrfbdy%ttxs_id) )
if(debug) write(6,*) ' ttxs_id = ',wrfbdy%ttxs_id
allocate( wrfbdy%ttxs( wrfbdy%sn, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

call check ( nf90_inq_varid(wrfbdy%ncid, "T_BTXE", wrfbdy%ttxe_id) )
if(debug) write(6,*) ' ttxe_id = ',wrfbdy%ttxe_id
allocate( wrfbdy%ttxe( wrfbdy%sn, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

call check ( nf90_inq_varid(wrfbdy%ncid, "T_BTYS", wrfbdy%ttys_id) )
if(debug) write(6,*) ' ttys_id = ',wrfbdy%ttys_id
allocate( wrfbdy%ttys( wrfbdy%we, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

call check ( nf90_inq_varid(wrfbdy%ncid, "T_BTYE", wrfbdy%ttye_id) )
if(debug) write(6,*) ' ttye_id = ',wrfbdy%ttye_id
allocate( wrfbdy%ttye( wrfbdy%we, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

  !-- mu on bdy
call check ( nf90_inq_varid(wrfbdy%ncid, "MU_BXS", wrfbdy%muxs_id) )
if(debug) write(6,*) ' muxs_id = ',wrfbdy%muxs_id
allocate( wrfbdy%muxs( wrfbdy%sn, wrfbdy%bdywdth, wrfbdy%time ) )

call check ( nf90_inq_varid(wrfbdy%ncid, "MU_BXE", wrfbdy%muxe_id) )
if(debug) write(6,*) ' muxe_id = ',wrfbdy%muxe_id
allocate( wrfbdy%muxe( wrfbdy%sn, wrfbdy%bdywdth, wrfbdy%time ) )

call check ( nf90_inq_varid(wrfbdy%ncid, "MU_BYS", wrfbdy%muys_id) )
if(debug) write(6,*) ' muys_id = ',wrfbdy%muys_id
allocate( wrfbdy%muys( wrfbdy%we, wrfbdy%bdywdth, wrfbdy%time ) )

call check ( nf90_inq_varid(wrfbdy%ncid, "MU_BYE", wrfbdy%muye_id) )
if(debug) write(6,*) ' muye_id = ',wrfbdy%muye_id
allocate( wrfbdy%muye( wrfbdy%we, wrfbdy%bdywdth, wrfbdy%time ) )

  !-- mu tendency on bdy
call check ( nf90_inq_varid(wrfbdy%ncid, "MU_BTXS", wrfbdy%mutxs_id) )
if(debug) write(6,*) ' mutxs_id = ',wrfbdy%mutxs_id
allocate( wrfbdy%mutxs( wrfbdy%sn, wrfbdy%bdywdth, wrfbdy%time ) )

call check ( nf90_inq_varid(wrfbdy%ncid, "MU_BTXE", wrfbdy%mutxe_id) )
if(debug) write(6,*) ' mutxe_id = ',wrfbdy%mutxe_id
allocate( wrfbdy%mutxe( wrfbdy%sn, wrfbdy%bdywdth, wrfbdy%time ) )

call check ( nf90_inq_varid(wrfbdy%ncid, "MU_BTYS", wrfbdy%mutys_id) )
if(debug) write(6,*) ' mutys_id = ',wrfbdy%mutys_id
allocate( wrfbdy%mutys( wrfbdy%we, wrfbdy%bdywdth, wrfbdy%time ) )

call check ( nf90_inq_varid(wrfbdy%ncid, "MU_BTYE", wrfbdy%mutye_id) )
if(debug) write(6,*) ' mutye_id = ',wrfbdy%mutye_id
allocate( wrfbdy%mutye( wrfbdy%we, wrfbdy%bdywdth, wrfbdy%time ) )

if(wrfbdy%n_moist > 0) then
   !-- qv on bdy
   call check ( nf90_inq_varid(wrfbdy%ncid, "QVAPOR_BXS", wrfbdy%qvxs_id) )
   allocate( wrfbdy%qvxs( wrfbdy%sn, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QVAPOR_BXE", wrfbdy%qvxe_id) )
   allocate( wrfbdy%qvxe( wrfbdy%sn, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QVAPOR_BYS", wrfbdy%qvys_id) )
   allocate( wrfbdy%qvys( wrfbdy%we, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QVAPOR_BYE", wrfbdy%qvye_id) )
   allocate( wrfbdy%qvye( wrfbdy%we, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   !-- qv tendency on bdy
   call check ( nf90_inq_varid(wrfbdy%ncid, "QVAPOR_BTXS", wrfbdy%qvtxs_id) )
   allocate( wrfbdy%qvtxs( wrfbdy%sn, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QVAPOR_BTXE", wrfbdy%qvtxe_id) )
   allocate( wrfbdy%qvtxe( wrfbdy%sn, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QVAPOR_BTYS", wrfbdy%qvtys_id) )
   allocate( wrfbdy%qvtys( wrfbdy%we, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QVAPOR_BTYE", wrfbdy%qvtye_id) )
   allocate( wrfbdy%qvtye( wrfbdy%we, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )
endif

if(wrfbdy%n_moist > 1) then
   !-- qc on bdy
   call check ( nf90_inq_varid(wrfbdy%ncid, "QCLOUD_BXS", wrfbdy%qcxs_id) )
   allocate( wrfbdy%qcxs( wrfbdy%sn, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QCLOUD_BXE", wrfbdy%qcxe_id) )
   allocate( wrfbdy%qcxe( wrfbdy%sn, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QCLOUD_BYS", wrfbdy%qcys_id) )
   allocate( wrfbdy%qcys( wrfbdy%we, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QCLOUD_BYE", wrfbdy%qcye_id) )
   allocate( wrfbdy%qcye( wrfbdy%we, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   !-- qc tendency on bdy
   call check ( nf90_inq_varid(wrfbdy%ncid, "QCLOUD_BTXS", wrfbdy%qctxs_id) )
   allocate( wrfbdy%qctxs( wrfbdy%sn, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QCLOUD_BTXE", wrfbdy%qctxe_id) )
   allocate( wrfbdy%qctxe( wrfbdy%sn, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QCLOUD_BTYS", wrfbdy%qctys_id) )
   allocate( wrfbdy%qctys( wrfbdy%we, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QCLOUD_BTYE", wrfbdy%qctye_id) )
   allocate( wrfbdy%qctye( wrfbdy%we, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )
endif

if(wrfbdy%n_moist > 2) then
   !-- qr on bdy
   call check ( nf90_inq_varid(wrfbdy%ncid, "QRAIN_BXS", wrfbdy%qrxs_id) )
   allocate( wrfbdy%qrxs( wrfbdy%sn, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QRAIN_BXE", wrfbdy%qrxe_id) )
   allocate( wrfbdy%qrxe( wrfbdy%sn, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QRAIN_BYS", wrfbdy%qrys_id) )
   allocate( wrfbdy%qrys( wrfbdy%we, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QRAIN_BYE", wrfbdy%qrye_id) )
   allocate( wrfbdy%qrye( wrfbdy%we, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   !-- qr tendency on bdy
   call check ( nf90_inq_varid(wrfbdy%ncid, "QRAIN_BTXS", wrfbdy%qrtxs_id) )
   allocate( wrfbdy%qrtxs( wrfbdy%sn, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QRAIN_BTXE", wrfbdy%qrtxe_id) )
   allocate( wrfbdy%qrtxe( wrfbdy%sn, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QRAIN_BTYS", wrfbdy%qrtys_id) )
   allocate( wrfbdy%qrtys( wrfbdy%we, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QRAIN_BTYE", wrfbdy%qrtye_id) )
   allocate( wrfbdy%qrtye( wrfbdy%we, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )
endif

if(wrfbdy%n_moist > 3) then
   !-- qi on bdy
   call check ( nf90_inq_varid(wrfbdy%ncid, "QICE_BXS", wrfbdy%qixs_id) )
   allocate( wrfbdy%qixs( wrfbdy%sn, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QICE_BXE", wrfbdy%qixe_id) )
   allocate( wrfbdy%qixe( wrfbdy%sn, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QICE_BYS", wrfbdy%qiys_id) )
   allocate( wrfbdy%qiys( wrfbdy%we, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QICE_BYE", wrfbdy%qiye_id) )
   allocate( wrfbdy%qiye( wrfbdy%we, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   !-- qi tendency on bdy
   call check ( nf90_inq_varid(wrfbdy%ncid, "QICE_BTXS", wrfbdy%qitxs_id) )
   allocate( wrfbdy%qitxs( wrfbdy%sn, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QICE_BTXE", wrfbdy%qitxe_id) )
   allocate( wrfbdy%qitxe( wrfbdy%sn, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QICE_BTYS", wrfbdy%qitys_id) )
   allocate( wrfbdy%qitys( wrfbdy%we, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QICE_BTYE", wrfbdy%qitye_id) )
   allocate( wrfbdy%qitye( wrfbdy%we, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )
endif

if(wrfbdy%n_moist > 4) then
   !-- qs on bdy
   call check ( nf90_inq_varid(wrfbdy%ncid, "QSNOW_BXS", wrfbdy%qsxs_id) )
   allocate( wrfbdy%qsxs( wrfbdy%sn, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QSNOW_BXE", wrfbdy%qsxe_id) )
   allocate( wrfbdy%qsxe( wrfbdy%sn, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QSNOW_BYS", wrfbdy%qsys_id) )
   allocate( wrfbdy%qsys( wrfbdy%we, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QSNOW_BYE", wrfbdy%qsye_id) )
   allocate( wrfbdy%qsye( wrfbdy%we, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   !-- qs tendency on bdy
   call check ( nf90_inq_varid(wrfbdy%ncid, "QSNOW_BTXS", wrfbdy%qstxs_id) )
   allocate( wrfbdy%qstxs( wrfbdy%sn, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QSNOW_BTXE", wrfbdy%qstxe_id) )
   allocate( wrfbdy%qstxe( wrfbdy%sn, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QSNOW_BTYS", wrfbdy%qstys_id) )
   allocate( wrfbdy%qstys( wrfbdy%we, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QSNOW_BTYE", wrfbdy%qstye_id) )
   allocate( wrfbdy%qstye( wrfbdy%we, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )
endif

if(wrfbdy%n_moist > 5) then
   !-- qg on bdy
   call check ( nf90_inq_varid(wrfbdy%ncid, "QGRAUP_BXS", wrfbdy%qgxs_id) )
   allocate( wrfbdy%qgxs( wrfbdy%sn, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QGRAUP_BXE", wrfbdy%qgxe_id) )
   allocate( wrfbdy%qgxe( wrfbdy%sn, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QGRAUP_BYS", wrfbdy%qgys_id) )
   allocate( wrfbdy%qgys( wrfbdy%we, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QGRAUP_BYE", wrfbdy%qgye_id) )
   allocate( wrfbdy%qgye( wrfbdy%we, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   !-- qg tendency on bdy
   call check ( nf90_inq_varid(wrfbdy%ncid, "QGRAUP_BTXS", wrfbdy%qgtxs_id) )
   allocate( wrfbdy%qgtxs( wrfbdy%sn, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QGRAUP_BTXE", wrfbdy%qgtxe_id) )
   allocate( wrfbdy%qgtxe( wrfbdy%sn, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QGRAUP_BTYS", wrfbdy%qgtys_id) )
   allocate( wrfbdy%qgtys( wrfbdy%we, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QGRAUP_BTYE", wrfbdy%qgtye_id) )
   allocate( wrfbdy%qgtye( wrfbdy%we, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )
endif

if(wrfbdy%n_moist > 6) then
   !-- qnice on bdy
   call check ( nf90_inq_varid(wrfbdy%ncid, "QNICE_BXS", wrfbdy%qnicexs_id) )
   allocate( wrfbdy%qnicexs( wrfbdy%sn, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QNICE_BXE", wrfbdy%qnicexe_id) )
   allocate( wrfbdy%qnicexe( wrfbdy%sn, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QNICE_BYS", wrfbdy%qniceys_id) )
   allocate( wrfbdy%qniceys( wrfbdy%we, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QNICE_BYE", wrfbdy%qniceye_id) )
   allocate( wrfbdy%qniceye( wrfbdy%we, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   !-- qnice tendency on bdy
   call check ( nf90_inq_varid(wrfbdy%ncid, "QNICE_BTXS", wrfbdy%qnicetxs_id) )
   allocate( wrfbdy%qnicetxs( wrfbdy%sn, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QNICE_BTXE", wrfbdy%qnicetxe_id) )
   allocate( wrfbdy%qnicetxe( wrfbdy%sn, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QNICE_BTYS", wrfbdy%qnicetys_id) )
   allocate( wrfbdy%qnicetys( wrfbdy%we, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )

   call check ( nf90_inq_varid(wrfbdy%ncid, "QNICE_BTYE", wrfbdy%qnicetye_id) )
   allocate( wrfbdy%qnicetye( wrfbdy%we, wrfbdy%bt, wrfbdy%bdywdth, wrfbdy%time ) )
endif

if(wrfbdy%n_moist > 7) then
   write(6,*) 'n_moist = ',wrfbdy%n_moist
   call error_handler(E_ERR,'wrfbdy_open_and_alloc', &
         'n_moist is too large.', source, revision, revdate)
endif

contains

  ! Internal subroutine - checks error status after each netcdf, prints 
  !                       text message each time an error code is returned. 
  subroutine check(istatus)
    integer, intent ( in) :: istatus 
    if(istatus /= nf90_noerr) call error_handler(E_ERR,'wrfbdy_open_and_alloc', &
         trim(nf90_strerror(istatus)), source, revision, revdate)
  end subroutine check

end subroutine wrfbdy_open_and_alloc

!---------------------------------------------------------------

subroutine wrfbdy_dealloc( wrfbdy )

implicit none

type(wrf_bdy_data) :: wrfbdy

if ( .not. module_initialized ) call initialize_module

  !-- u on bdy
deallocate( wrfbdy%uxs )
deallocate( wrfbdy%uxe )
deallocate( wrfbdy%uys )
deallocate( wrfbdy%uye )

  !-- u tendency on bdy
deallocate( wrfbdy%utxs )
deallocate( wrfbdy%utxe )
deallocate( wrfbdy%utys )
deallocate( wrfbdy%utye )

  !-- v on bdy
deallocate( wrfbdy%vxs )
deallocate( wrfbdy%vxe )
deallocate( wrfbdy%vys )
deallocate( wrfbdy%vye )

  !-- v tendency on bdy
deallocate( wrfbdy%vtxs )
deallocate( wrfbdy%vtxe )
deallocate( wrfbdy%vtys )
deallocate( wrfbdy%vtye )

  !-- w on bdy
deallocate( wrfbdy%wxs )
deallocate( wrfbdy%wxe )
deallocate( wrfbdy%wys )
deallocate( wrfbdy%wye )

  !-- w tendency on bdy
deallocate( wrfbdy%wtxs )
deallocate( wrfbdy%wtxe )
deallocate( wrfbdy%wtys )
deallocate( wrfbdy%wtye )

  !-- height on bdy
deallocate( wrfbdy%phxs )
deallocate( wrfbdy%phxe )
deallocate( wrfbdy%phys )
deallocate( wrfbdy%phye )

  !-- height tendency on bdy
deallocate( wrfbdy%phtxs )
deallocate( wrfbdy%phtxe )
deallocate( wrfbdy%phtys )
deallocate( wrfbdy%phtye )

  !-- t on bdy
deallocate( wrfbdy%txs )
deallocate( wrfbdy%txe )
deallocate( wrfbdy%tys )
deallocate( wrfbdy%tye )

  !-- t tendency on bdy
deallocate( wrfbdy%ttxs )
deallocate( wrfbdy%ttxe )
deallocate( wrfbdy%ttys )
deallocate( wrfbdy%ttye )

  !-- mu on bdy
deallocate( wrfbdy%muxs )
deallocate( wrfbdy%muxe )
deallocate( wrfbdy%muys )
deallocate( wrfbdy%muye )

  !-- mu tendency on bdy
deallocate( wrfbdy%mutxs )
deallocate( wrfbdy%mutxe )
deallocate( wrfbdy%mutys )
deallocate( wrfbdy%mutye )

if(wrfbdy%n_moist > 0) then
   !-- qv on bdy
   deallocate( wrfbdy%qvxs )
   deallocate( wrfbdy%qvxe )
   deallocate( wrfbdy%qvys )
   deallocate( wrfbdy%qvye )

   !-- qv tendency on bdy
   deallocate( wrfbdy%qvtxs )
   deallocate( wrfbdy%qvtxe )
   deallocate( wrfbdy%qvtys )
   deallocate( wrfbdy%qvtye )
endif

if(wrfbdy%n_moist > 1) then
   !-- qc on bdy
   deallocate( wrfbdy%qcxs )
   deallocate( wrfbdy%qcxe )
   deallocate( wrfbdy%qcys )
   deallocate( wrfbdy%qcye )

   !-- qc tendency on bdy
   deallocate( wrfbdy%qctxs )
   deallocate( wrfbdy%qctxe )
   deallocate( wrfbdy%qctys )
   deallocate( wrfbdy%qctye )
endif

if(wrfbdy%n_moist > 2) then
   !-- qr on bdy
   deallocate( wrfbdy%qrxs )
   deallocate( wrfbdy%qrxe )
   deallocate( wrfbdy%qrys )
   deallocate( wrfbdy%qrye )

   !-- qr tendency on bdy
   deallocate( wrfbdy%qrtxs )
   deallocate( wrfbdy%qrtxe )
   deallocate( wrfbdy%qrtys )
   deallocate( wrfbdy%qrtye )
endif

if(wrfbdy%n_moist > 3) then
   !-- qi on bdy
   deallocate( wrfbdy%qixs )
   deallocate( wrfbdy%qixe )
   deallocate( wrfbdy%qiys )
   deallocate( wrfbdy%qiye )

   !-- qi tendency on bdy
   deallocate( wrfbdy%qitxs )
   deallocate( wrfbdy%qitxe )
   deallocate( wrfbdy%qitys )
   deallocate( wrfbdy%qitye )
endif

if(wrfbdy%n_moist > 4) then
   !-- qs on bdy
   deallocate( wrfbdy%qsxs )
   deallocate( wrfbdy%qsxe )
   deallocate( wrfbdy%qsys )
   deallocate( wrfbdy%qsye )

   !-- qs tendency on bdy
   deallocate( wrfbdy%qstxs )
   deallocate( wrfbdy%qstxe )
   deallocate( wrfbdy%qstys )
   deallocate( wrfbdy%qstye )
endif

if(wrfbdy%n_moist > 5) then
   !-- qg on bdy
   deallocate( wrfbdy%qgxs )
   deallocate( wrfbdy%qgxe )
   deallocate( wrfbdy%qgys )
   deallocate( wrfbdy%qgye )

   !-- qg tendency on bdy
   deallocate( wrfbdy%qgtxs )
   deallocate( wrfbdy%qgtxe )
   deallocate( wrfbdy%qgtys )
   deallocate( wrfbdy%qgtye )
endif

if(wrfbdy%n_moist > 6) then
   !-- qnice on bdy
   deallocate( wrfbdy%qnicexs )
   deallocate( wrfbdy%qnicexe )
   deallocate( wrfbdy%qniceys )
   deallocate( wrfbdy%qniceye )

   !-- qnice tendency on bdy
   deallocate( wrfbdy%qnicetxs )
   deallocate( wrfbdy%qnicetxe )
   deallocate( wrfbdy%qnicetys )
   deallocate( wrfbdy%qnicetye )
endif

if(wrfbdy%n_moist > 7) then
   write(6,*) 'n_moist = ',wrfbdy%n_moist
   call error_handler(E_ERR,'wrfbdy_dealloc', &
         'n_moist is too large.', source, revision, revdate)
endif

end subroutine wrfbdy_dealloc

!******************************************************************************

subroutine wrf_io( wrf, in_or_out, debug )

implicit none

type(wrf_data)    :: wrf
character (len=6) :: in_or_out
logical           :: debug

integer :: k, ndims, lngth, dimids(5), istatus

!----------------------------------------------------------------------

if ( .not. module_initialized ) call initialize_module

if(debug) then

   if (in_or_out  == "OUTPUT") then
      call error_handler(E_MSG,'wrf_io','Writing to the WRF restart netCDF file.', &
           source,revision,revdate)
   else
      call error_handler(E_MSG,'wrf_io','Reading the WRF restart netCDF file.', &
           source,revision,revdate)
   endif

endif

!----------------------------------------------------------------------
! Reading or Writing the variables. ignoring count, stride, map ...
!----------------------------------------------------------------------

call check( nf90_inquire_variable(wrf%ncid, wrf%u_id, ndims=ndims, dimids=dimids) )
call check( nf90_inquire_dimension(wrf%ncid, dimids(ndims), len=lngth) )

if(debug) write(6,*) 'len= ',lngth,' n_moist = ',wrf%n_moist

if (in_or_out  == "OUTPUT") then
   call check( nf90_put_var(wrf%ncid, wrf%u_id,    wrf%u,    start = (/ 1, 1, 1, 1 /)))
   call check( nf90_put_var(wrf%ncid, wrf%v_id,    wrf%v,    start = (/ 1, 1, 1, 1 /)))
   call check( nf90_put_var(wrf%ncid, wrf%w_id,    wrf%w,    start = (/ 1, 1, 1, 1 /)))
   call check( nf90_put_var(wrf%ncid, wrf%ph_id,   wrf%ph,   start = (/ 1, 1, 1, 1 /)))
   call check( nf90_put_var(wrf%ncid, wrf%t_id,    wrf%t,    start = (/ 1, 1, 1, 1 /)))
   call check( nf90_put_var(wrf%ncid, wrf%mu_id,   wrf%mu,   start = (/ 1, 1, 1 /)))
   if(wrf%n_moist > 0) then
      call check( nf90_put_var(wrf%ncid, wrf%qv_id, wrf%qv, start = (/ 1, 1, 1, 1 /)))
   endif
   if(wrf%n_moist > 1) then
      call check( nf90_put_var(wrf%ncid, wrf%qc_id, wrf%qc, start = (/ 1, 1, 1, 1 /)))
   endif
   if(wrf%n_moist > 2) then
      call check( nf90_put_var(wrf%ncid, wrf%qr_id, wrf%qr, start = (/ 1, 1, 1, 1 /)))
   endif
   if(wrf%n_moist > 3) then
      call check( nf90_put_var(wrf%ncid, wrf%qi_id, wrf%qi, start = (/ 1, 1, 1, 1 /)))
   endif
   if(wrf%n_moist > 4) then
      call check( nf90_put_var(wrf%ncid, wrf%qs_id, wrf%qs, start = (/ 1, 1, 1, 1 /)))
   endif
   if(wrf%n_moist > 5) then
      call check( nf90_put_var(wrf%ncid, wrf%qg_id, wrf%qg, start = (/ 1, 1, 1, 1 /)))
   endif
   if(wrf%n_moist > 6) then
      call check( nf90_put_var(wrf%ncid, wrf%qnice_id, wrf%qnice, start = (/ 1, 1, 1, 1 /)))
   endif
   if(wrf%n_moist > 7) then
      write(6,*) 'n_moist = ',wrf%n_moist
      call error_handler(E_ERR,'wrf_io', &
           'n_moist is too large.', source, revision, revdate)
   endif
   if( wrf%surf_obs ) then
      call check( nf90_put_var(wrf%ncid, wrf%u10_id, wrf%u10, start = (/ 1, 1, 1 /)))
      call check( nf90_put_var(wrf%ncid, wrf%v10_id, wrf%v10, start = (/ 1, 1, 1 /)))
      call check( nf90_put_var(wrf%ncid, wrf%t2_id, wrf%t2, start = (/ 1, 1, 1 /)))
      call check( nf90_put_var(wrf%ncid, wrf%th2_id, wrf%th2, start = (/ 1, 1, 1 /)))
      call check( nf90_put_var(wrf%ncid, wrf%q2_id, wrf%q2, start = (/ 1, 1, 1 /)))
      istatus = nf90_put_var(wrf%ncid, wrf%ps_id, wrf%ps, start = (/ 1, 1, 1 /))
      if(istatus /= nf90_noerr) then
         call error_handler(E_MSG,'PSFC', &
              trim(nf90_strerror(istatus)), source, revision, revdate)
         call error_handler(E_MSG,'wrf_io', &
              'creates PSFC', source, revision, revdate)
         call check(nf90_Redef(wrf%ncid))
         call check(nf90_def_var(wrf%ncid, name="PSFC", xtype=nf90_real, &
              dimids= (/ wrf%we_id,  wrf%sn_id/), varid=wrf%ps_id) )
         call check(nf90_enddef(wrf%ncid))
         call check( nf90_put_var(wrf%ncid, wrf%ps_id, wrf%ps, start = (/ 1, 1, 1 /)))
      endif
   endif
   if( wrf%soil_data ) then   
      call check( nf90_put_var(wrf%ncid, wrf%tslb_id, wrf%tslb, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrf%ncid, wrf%smois_id, wrf%smois, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrf%ncid, wrf%sh2o_id, wrf%sh2o, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrf%ncid, wrf%tsk_id,  wrf%tsk,  start = (/ 1, 1, 1 /)))
   endif
   if( wrf%h_diab ) then
      istatus = nf90_put_var(wrf%ncid, wrf%hdiab_id, wrf%hdiab, start = (/ 1, 1, 1, 1 /))
      if(istatus /= nf90_noerr) then
         call error_handler(E_ERR,'H_DIABATIC', &
              trim(nf90_strerror(istatus)), source, revision, revdate)
      endif
   endif
else
   call check( nf90_get_var(wrf%ncid, wrf%u_id,    wrf%u,    start = (/ 1, 1, 1, lngth /)))
   call check( nf90_get_var(wrf%ncid, wrf%v_id,    wrf%v,    start = (/ 1, 1, 1, lngth /)))
   call check( nf90_get_var(wrf%ncid, wrf%w_id,    wrf%w,    start = (/ 1, 1, 1, lngth /)))
   call check( nf90_get_var(wrf%ncid, wrf%ph_id,   wrf%ph,   start = (/ 1, 1, 1, lngth /)))
   call check( nf90_get_var(wrf%ncid, wrf%phb_id,  wrf%phb,  start = (/ 1, 1, 1, lngth /)))
   call check( nf90_get_var(wrf%ncid, wrf%t_id,    wrf%t,    start = (/ 1, 1, 1, lngth /)))
   call check( nf90_get_var(wrf%ncid, wrf%mu_id,   wrf%mu,   start = (/ 1, 1, lngth /)))
   call check( nf90_get_var(wrf%ncid, wrf%mub_id,  wrf%mub,  start = (/ 1, 1, lngth /)))
   if(wrf%n_moist > 0) then
      call check( nf90_get_var(wrf%ncid, wrf%qv_id, wrf%qv, start = (/ 1, 1, 1, lngth /)))
   endif
   if(wrf%n_moist > 1) then
      call check( nf90_get_var(wrf%ncid, wrf%qc_id, wrf%qc, start = (/ 1, 1, 1, lngth /)))
   endif
   if(wrf%n_moist > 2) then
      call check( nf90_get_var(wrf%ncid, wrf%qr_id, wrf%qr, start = (/ 1, 1, 1, lngth /)))
   endif
   if(wrf%n_moist > 3) then
      call check( nf90_get_var(wrf%ncid, wrf%qi_id, wrf%qi, start = (/ 1, 1, 1, lngth /)))
   endif
   if(wrf%n_moist > 4) then
      call check( nf90_get_var(wrf%ncid, wrf%qs_id, wrf%qs, start = (/ 1, 1, 1, lngth /)))
   endif
   if(wrf%n_moist > 5) then
      call check( nf90_get_var(wrf%ncid, wrf%qg_id, wrf%qg, start = (/ 1, 1, 1, lngth /)))
   endif
   if(wrf%n_moist > 6) then
      call check( nf90_get_var(wrf%ncid, wrf%qnice_id, wrf%qnice, start = (/ 1, 1, 1, lngth /)))
   endif
   if(wrf%n_moist > 7) then
      write(6,*) 'n_moist = ',wrf%n_moist
      call error_handler(E_ERR,'wrf_io', &
           'n_moist is too large.', source, revision, revdate)
   endif
   if( wrf%surf_obs ) then
      call check( nf90_get_var(wrf%ncid, wrf%u10_id,  wrf%u10,  start = (/ 1, 1, lngth /)))
      call check( nf90_get_var(wrf%ncid, wrf%v10_id,  wrf%v10,  start = (/ 1, 1, lngth /)))
      call check( nf90_get_var(wrf%ncid, wrf%t2_id,  wrf%t2,  start = (/ 1, 1, lngth /)))
      call check( nf90_get_var(wrf%ncid, wrf%th2_id,  wrf%th2,  start = (/ 1, 1, lngth /)))
      call check( nf90_get_var(wrf%ncid, wrf%q2_id,  wrf%q2,  start = (/ 1, 1, lngth /)))
      istatus = nf90_inq_varid(wrf%ncid, "PSFC", wrf%ps_id)
      if(istatus == nf90_noerr) then
         call check( nf90_get_var(wrf%ncid, wrf%ps_id,  wrf%ps,  start = (/ 1, 1, lngth /)))
      else
         call error_handler(E_MSG,'PSFC', &
              trim(nf90_strerror(istatus)), source, revision, revdate)
         call error_handler(E_MSG,'wrf_io', &
              'sets PSFC to zero', source, revision, revdate)
         wrf%ps(:,:) = 0.0_r8
      endif
   endif
   if( wrf%soil_data ) then
      call check( nf90_get_var(wrf%ncid, wrf%tslb_id, wrf%tslb, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrf%ncid, wrf%smois_id, wrf%smois, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrf%ncid, wrf%sh2o_id, wrf%sh2o, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrf%ncid, wrf%tsk_id,  wrf%tsk,  start = (/ 1, 1, lngth /)))
   endif
   if( wrf%h_diab ) then
      istatus = nf90_get_var(wrf%ncid, wrf%hdiab_id,  wrf%hdiab,  start = (/ 1, 1, 1, lngth /))
      if(istatus /= nf90_noerr) then
         call error_handler(E_MSG,'H_DIABATIC', &
              trim(nf90_strerror(istatus)), source, revision, revdate)
         call error_handler(E_MSG,'wrf_io', &
              'sets H_DIABATIC to zero', source, revision, revdate)
         wrf%hdiab(:,:,:) = 0.0_r8
      endif
   endif
endif

if(debug) then

   do k=1,wrf%bt
      write(6,*) ' k, corner vals for u '
      write(6,*) k, wrf%u(1,1,k),wrf%u(wrf%we+1,1,k),  &
           wrf%u(1,wrf%sn,k),wrf%u(wrf%we+1,wrf%sn,k)
   enddo

   write(6,*) ' '

   do k=1,wrf%bt
      write(6,*) ' k, corner vals for v '
      write(6,*) k, wrf%v(1,1,k),wrf%v(wrf%we,1,k),  &
           wrf%v(1,wrf%sn+1,k),wrf%v(wrf%we,wrf%sn+1,k)
   enddo

   write(6,*) ' '

   write(6,*) ' corner vals for mu '
   write(6,*) wrf%mu(1,1),wrf%mu(wrf%we,1),  &
        wrf%mu(1,wrf%sn),wrf%mu(wrf%we,wrf%sn)

   write(6,*) ' '

   write(6,*) ' corner vals for mub '
   write(6,*) wrf%mub(1,1),wrf%mub(wrf%we,1),  &
        wrf%mub(1,wrf%sn),wrf%mub(wrf%we,wrf%sn)

endif

contains

  ! Internal subroutine - checks error status after each netcdf, prints 
  !                       text message each time an error code is returned. 
  subroutine check(istatus)
    integer, intent ( in) :: istatus 
    if(istatus /= nf90_noerr) call error_handler(E_ERR,'wrf_io', &
         trim(nf90_strerror(istatus)), source, revision, revdate) 
  end subroutine check

end subroutine wrf_io

!---------------------------------------------------------------

subroutine wrfbdy_io( wrfbdy, in_or_out, debug )

! Reads or writes wrfbdy; includes all variables (bdy fields
! and tendencies) from file specified in wrfbdy%ncid.

implicit none

type(wrf_bdy_data) :: wrfbdy
character (len=6)  :: in_or_out
logical            :: debug

integer :: ndims, lngth, dimids(5)

!----------------------------------------------------------------------

if ( .not. module_initialized ) call initialize_module

if(debug) then

   if (in_or_out  == "OUTPUT") then
      call error_handler(E_MSG,'wrfbdy_io','Writing to the WRF restart netCDF file.', &
           source,revision,revdate)
   else
      call error_handler(E_MSG,'wrfbdy_io','Reading the WRF restart netCDF file.', &
           source,revision,revdate)
   endif

endif

call check( nf90_inquire_variable(wrfbdy%ncid, wrfbdy%uxs_id, ndims=ndims, dimids=dimids) )
call check( nf90_inquire_dimension(wrfbdy%ncid, dimids(ndims), len=lngth) )

if (in_or_out  == "OUTPUT") then
   !-- u on boundary
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%uxs_id, wrfbdy%uxs, start = (/ 1, 1, 1, 1 /)))
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%uxe_id, wrfbdy%uxe, start = (/ 1, 1, 1, 1 /)))
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%uys_id, wrfbdy%uys, start = (/ 1, 1, 1, 1 /)))
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%uye_id, wrfbdy%uye, start = (/ 1, 1, 1, 1 /)))
   !-- u tendencies on boundary
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%utxs_id, wrfbdy%utxs, start = (/ 1, 1, 1, 1 /)))
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%utxe_id, wrfbdy%utxe, start = (/ 1, 1, 1, 1 /)))
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%utys_id, wrfbdy%utys, start = (/ 1, 1, 1, 1 /)))
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%utye_id, wrfbdy%utye, start = (/ 1, 1, 1, 1 /)))
   !-- v on boundary
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%vxs_id, wrfbdy%vxs, start = (/ 1, 1, 1, 1 /)))
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%vxe_id, wrfbdy%vxe, start = (/ 1, 1, 1, 1 /)))
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%vys_id, wrfbdy%vys, start = (/ 1, 1, 1, 1 /)))
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%vye_id, wrfbdy%vye, start = (/ 1, 1, 1, 1 /)))
   !-- v tendencies on boundary
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%vtxs_id, wrfbdy%vtxs, start = (/ 1, 1, 1, 1 /)))
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%vtxe_id, wrfbdy%vtxe, start = (/ 1, 1, 1, 1 /)))
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%vtys_id, wrfbdy%vtys, start = (/ 1, 1, 1, 1 /)))
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%vtye_id, wrfbdy%vtye, start = (/ 1, 1, 1, 1 /)))
   !-- w on boundary
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%wxs_id, wrfbdy%wxs, start = (/ 1, 1, 1, 1 /)))
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%wxe_id, wrfbdy%wxe, start = (/ 1, 1, 1, 1 /)))
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%wys_id, wrfbdy%wys, start = (/ 1, 1, 1, 1 /)))
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%wye_id, wrfbdy%wye, start = (/ 1, 1, 1, 1 /)))
   !-- w tendencies on boundary
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%wtxs_id, wrfbdy%wtxs, start = (/ 1, 1, 1, 1 /)))
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%wtxe_id, wrfbdy%wtxe, start = (/ 1, 1, 1, 1 /)))
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%wtys_id, wrfbdy%wtys, start = (/ 1, 1, 1, 1 /)))
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%wtye_id, wrfbdy%wtye, start = (/ 1, 1, 1, 1 /)))
   !-- height on boundary
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%phxs_id, wrfbdy%phxs, start = (/ 1, 1, 1, 1 /)))
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%phxe_id, wrfbdy%phxe, start = (/ 1, 1, 1, 1 /)))
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%phys_id, wrfbdy%phys, start = (/ 1, 1, 1, 1 /)))
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%phye_id, wrfbdy%phye, start = (/ 1, 1, 1, 1 /)))
   !-- height tendencies on boundary
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%phtxs_id, wrfbdy%phtxs, start = (/ 1, 1, 1, 1 /)))
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%phtxe_id, wrfbdy%phtxe, start = (/ 1, 1, 1, 1 /)))
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%phtys_id, wrfbdy%phtys, start = (/ 1, 1, 1, 1 /)))
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%phtye_id, wrfbdy%phtye, start = (/ 1, 1, 1, 1 /)))
   !-- t on boundary
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%txs_id, wrfbdy%txs, start = (/ 1, 1, 1, 1 /)))
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%txe_id, wrfbdy%txe, start = (/ 1, 1, 1, 1 /)))
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%tys_id, wrfbdy%tys, start = (/ 1, 1, 1, 1 /)))
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%tye_id, wrfbdy%tye, start = (/ 1, 1, 1, 1 /)))
   !-- t tendencies on boundary
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%ttxs_id, wrfbdy%ttxs, start = (/ 1, 1, 1, 1 /)))
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%ttxe_id, wrfbdy%ttxe, start = (/ 1, 1, 1, 1 /)))
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%ttys_id, wrfbdy%ttys, start = (/ 1, 1, 1, 1 /)))
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%ttye_id, wrfbdy%ttye, start = (/ 1, 1, 1, 1 /)))
   !-- mu on boundary
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%muxs_id, wrfbdy%muxs, start = (/ 1, 1, 1 /)))
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%muxe_id, wrfbdy%muxe, start = (/ 1, 1, 1 /)))
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%muys_id, wrfbdy%muys, start = (/ 1, 1, 1 /)))
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%muye_id, wrfbdy%muye, start = (/ 1, 1, 1 /)))
   !-- mu tendencies on boundary
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%mutxs_id, wrfbdy%mutxs, start = (/ 1, 1, 1 /)))
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%mutxe_id, wrfbdy%mutxe, start = (/ 1, 1, 1 /)))
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%mutys_id, wrfbdy%mutys, start = (/ 1, 1, 1 /)))
   call check( nf90_put_var(wrfbdy%ncid, wrfbdy%mutye_id, wrfbdy%mutye, start = (/ 1, 1, 1 /)))
   if(wrfbdy%n_moist > 0) then
      !-- qv on boundary
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qvxs_id, wrfbdy%qvxs, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qvxe_id, wrfbdy%qvxe, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qvys_id, wrfbdy%qvys, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qvye_id, wrfbdy%qvye, start = (/ 1, 1, 1, 1 /)))
      !-- qv tendencies on boundary
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qvtxs_id, wrfbdy%qvtxs, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qvtxe_id, wrfbdy%qvtxe, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qvtys_id, wrfbdy%qvtys, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qvtye_id, wrfbdy%qvtye, start = (/ 1, 1, 1, 1 /)))
   endif
   if(wrfbdy%n_moist > 1) then
      !-- qc on boundary
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qcxs_id, wrfbdy%qcxs, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qcxe_id, wrfbdy%qcxe, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qcys_id, wrfbdy%qcys, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qcye_id, wrfbdy%qcye, start = (/ 1, 1, 1, 1 /)))
      !-- qc tendencies on boundary
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qctxs_id, wrfbdy%qctxs, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qctxe_id, wrfbdy%qctxe, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qctys_id, wrfbdy%qctys, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qctye_id, wrfbdy%qctye, start = (/ 1, 1, 1, 1 /)))
   endif
   if(wrfbdy%n_moist > 2) then
      !-- qr on boundary
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qrxs_id, wrfbdy%qrxs, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qrxe_id, wrfbdy%qrxe, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qrys_id, wrfbdy%qrys, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qrye_id, wrfbdy%qrye, start = (/ 1, 1, 1, 1 /)))
      !-- qr tendencies on boundary
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qrtxs_id, wrfbdy%qrtxs, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qrtxe_id, wrfbdy%qrtxe, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qrtys_id, wrfbdy%qrtys, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qrtye_id, wrfbdy%qrtye, start = (/ 1, 1, 1, 1 /)))
   endif
   if(wrfbdy%n_moist > 3) then
      !-- qi on boundary
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qixs_id, wrfbdy%qixs, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qixe_id, wrfbdy%qixe, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qiys_id, wrfbdy%qiys, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qiye_id, wrfbdy%qiye, start = (/ 1, 1, 1, 1 /)))
      !-- qi tendencies on boundary
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qitxs_id, wrfbdy%qitxs, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qitxe_id, wrfbdy%qitxe, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qitys_id, wrfbdy%qitys, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qitye_id, wrfbdy%qitye, start = (/ 1, 1, 1, 1 /)))
   endif
   if(wrfbdy%n_moist > 4) then
      !-- qs on boundary
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qsxs_id, wrfbdy%qsxs, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qsxe_id, wrfbdy%qsxe, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qsys_id, wrfbdy%qsys, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qsye_id, wrfbdy%qsye, start = (/ 1, 1, 1, 1 /)))
      !-- qs tendencies on boundary
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qstxs_id, wrfbdy%qstxs, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qstxe_id, wrfbdy%qstxe, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qstys_id, wrfbdy%qstys, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qstye_id, wrfbdy%qstye, start = (/ 1, 1, 1, 1 /)))
   endif
   if(wrfbdy%n_moist > 5) then
      !-- qg on boundary
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qgxs_id, wrfbdy%qgxs, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qgxe_id, wrfbdy%qgxe, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qgys_id, wrfbdy%qgys, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qgye_id, wrfbdy%qgye, start = (/ 1, 1, 1, 1 /)))
      !-- qg tendencies on boundary
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qgtxs_id, wrfbdy%qgtxs, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qgtxe_id, wrfbdy%qgtxe, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qgtys_id, wrfbdy%qgtys, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qgtye_id, wrfbdy%qgtye, start = (/ 1, 1, 1, 1 /)))
   endif
   if(wrfbdy%n_moist > 6) then
      !-- qnice on boundary
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qnicexs_id, wrfbdy%qnicexs, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qnicexe_id, wrfbdy%qnicexe, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qniceys_id, wrfbdy%qniceys, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qniceye_id, wrfbdy%qniceye, start = (/ 1, 1, 1, 1 /)))
      !-- qnice tendencies on boundary
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qnicetxs_id, wrfbdy%qnicetxs, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qnicetxe_id, wrfbdy%qnicetxe, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qnicetys_id, wrfbdy%qnicetys, start = (/ 1, 1, 1, 1 /)))
      call check( nf90_put_var(wrfbdy%ncid, wrfbdy%qnicetye_id, wrfbdy%qnicetye, start = (/ 1, 1, 1, 1 /)))
   endif
   if(wrfbdy%n_moist > 7) then
      write(6,*) 'n_moist = ',wrfbdy%n_moist
      call error_handler(E_ERR,'wrfbdy_io', &
           'n_moist is too large.', source, revision, revdate)
   endif
else
   !-- u on boundary
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%uxs_id, wrfbdy%uxs, start = (/ 1, 1, 1, lngth /)))
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%uxe_id, wrfbdy%uxe, start = (/ 1, 1, 1, lngth /)))
   if(debug) write(6,*) ' calling netcdf read for uys '
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%uys_id, wrfbdy%uys, start = (/ 1, 1, 1, lngth /)))
   if(debug) write(6,*) ' calling netcdf read for uye '
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%uye_id, wrfbdy%uye, start = (/ 1, 1, 1, lngth /)))
   !-- u tendencies on boundary
   if(debug) write(6,*) ' calling netcdf read for utxs '
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%utxs_id, wrfbdy%utxs, start = (/ 1, 1, 1, lngth /)))
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%utxe_id, wrfbdy%utxe, start = (/ 1, 1, 1, lngth /)))
   if(debug) write(6,*) ' calling netcdf read for utys '
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%utys_id, wrfbdy%utys, start = (/ 1, 1, 1, lngth /)))
   if(debug) write(6,*) ' calling netcdf read for utye '
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%utye_id, wrfbdy%utye, start = (/ 1, 1, 1, lngth /)))
   !-- v on boundary
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%vxs_id, wrfbdy%vxs, start = (/ 1, 1, 1, lngth /)))
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%vxe_id, wrfbdy%vxe, start = (/ 1, 1, 1, lngth /)))
   if(debug) write(6,*) ' calling netcdf read for vys '
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%vys_id, wrfbdy%vys, start = (/ 1, 1, 1, lngth /)))
   if(debug) write(6,*) ' calling netcdf read for vye '
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%vye_id, wrfbdy%vye, start = (/ 1, 1, 1, lngth /)))
   !-- v tendencies on boundary
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%vtxs_id, wrfbdy%vtxs, start = (/ 1, 1, 1, lngth /)))
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%vtxe_id, wrfbdy%vtxe, start = (/ 1, 1, 1, lngth /)))
   if(debug) write(6,*) ' calling netcdf read for vtys '
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%vtys_id, wrfbdy%vtys, start = (/ 1, 1, 1, lngth /)))
   if(debug) write(6,*) ' calling netcdf read for vtye '
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%vtye_id, wrfbdy%vtye, start = (/ 1, 1, 1, lngth /)))
   !-- w on boundary
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%wxs_id, wrfbdy%wxs, start = (/ 1, 1, 1, lngth /)))
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%wxe_id, wrfbdy%wxe, start = (/ 1, 1, 1, lngth /)))
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%wys_id, wrfbdy%wys, start = (/ 1, 1, 1, lngth /)))
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%wye_id, wrfbdy%wye, start = (/ 1, 1, 1, lngth /)))
   !-- w tendencies on boundary
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%wtxs_id, wrfbdy%wtxs, start = (/ 1, 1, 1, lngth /)))
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%wtxe_id, wrfbdy%wtxe, start = (/ 1, 1, 1, lngth /)))
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%wtys_id, wrfbdy%wtys, start = (/ 1, 1, 1, lngth /)))
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%wtye_id, wrfbdy%wtye, start = (/ 1, 1, 1, lngth /)))
   !-- height on boundary
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%phxs_id, wrfbdy%phxs, start = (/ 1, 1, 1, lngth /)))
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%phxe_id, wrfbdy%phxe, start = (/ 1, 1, 1, lngth /)))
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%phys_id, wrfbdy%phys, start = (/ 1, 1, 1, lngth /)))
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%phye_id, wrfbdy%phye, start = (/ 1, 1, 1, lngth /)))
   !-- height tendencies on boundary
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%phtxs_id, wrfbdy%phtxs, start = (/ 1, 1, 1, lngth /)))
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%phtxe_id, wrfbdy%phtxe, start = (/ 1, 1, 1, lngth /)))
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%phtys_id, wrfbdy%phtys, start = (/ 1, 1, 1, lngth /)))
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%phtye_id, wrfbdy%phtye, start = (/ 1, 1, 1, lngth /)))
   !-- t on boundary
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%txs_id, wrfbdy%txs, start = (/ 1, 1, 1, lngth /)))
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%txe_id, wrfbdy%txe, start = (/ 1, 1, 1, lngth /)))
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%tys_id, wrfbdy%tys, start = (/ 1, 1, 1, lngth /)))
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%tye_id, wrfbdy%tye, start = (/ 1, 1, 1, lngth /)))
   !-- t tendencies on boundary
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%ttxs_id, wrfbdy%ttxs, start = (/ 1, 1, 1, lngth /)))
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%ttxe_id, wrfbdy%ttxe, start = (/ 1, 1, 1, lngth /)))
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%ttys_id, wrfbdy%ttys, start = (/ 1, 1, 1, lngth /)))
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%ttye_id, wrfbdy%ttye, start = (/ 1, 1, 1, lngth /)))
   !-- mu on boundary
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%muxs_id, wrfbdy%muxs, start = (/ 1, 1, lngth /)))
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%muxe_id, wrfbdy%muxe, start = (/ 1, 1, lngth /)))
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%muys_id, wrfbdy%muys, start = (/ 1, 1, lngth /)))
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%muye_id, wrfbdy%muye, start = (/ 1, 1, lngth /)))
   !-- mu tendencies on boundary
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%mutxs_id, wrfbdy%mutxs, start = (/ 1, 1, lngth /)))
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%mutxe_id, wrfbdy%mutxe, start = (/ 1, 1, lngth /)))
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%mutys_id, wrfbdy%mutys, start = (/ 1, 1, lngth /)))
   call check( nf90_get_var(wrfbdy%ncid, wrfbdy%mutye_id, wrfbdy%mutye, start = (/ 1, 1, lngth /)))
   if(wrfbdy%n_moist > 0) then
      !-- qv on boundary
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qvxs_id, wrfbdy%qvxs, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qvxe_id, wrfbdy%qvxe, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qvys_id, wrfbdy%qvys, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qvye_id, wrfbdy%qvye, start = (/ 1, 1, 1, lngth /)))
      !-- qv tendencies on boundary
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qvtxs_id, wrfbdy%qvtxs, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qvtxe_id, wrfbdy%qvtxe, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qvtys_id, wrfbdy%qvtys, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qvtye_id, wrfbdy%qvtye, start = (/ 1, 1, 1, lngth /)))
   endif
   if(wrfbdy%n_moist > 1) then
      !-- qc on boundary
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qcxs_id, wrfbdy%qcxs, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qcxe_id, wrfbdy%qcxe, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qcys_id, wrfbdy%qcys, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qcye_id, wrfbdy%qcye, start = (/ 1, 1, 1, lngth /)))
      !-- qc tendencies on boundary
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qctxs_id, wrfbdy%qctxs, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qctxe_id, wrfbdy%qctxe, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qctys_id, wrfbdy%qctys, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qctye_id, wrfbdy%qctye, start = (/ 1, 1, 1, lngth /)))
   endif
   if(wrfbdy%n_moist > 2) then
      !-- qr on boundary
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qrxs_id, wrfbdy%qrxs, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qrxe_id, wrfbdy%qrxe, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qrys_id, wrfbdy%qrys, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qrye_id, wrfbdy%qrye, start = (/ 1, 1, 1, lngth /)))
      !-- qr tendencies on boundary
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qrtxs_id, wrfbdy%qrtxs, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qrtxe_id, wrfbdy%qrtxe, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qrtys_id, wrfbdy%qrtys, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qrtye_id, wrfbdy%qrtye, start = (/ 1, 1, 1, lngth /)))
   endif
   if(wrfbdy%n_moist > 3) then
      !-- qi on boundary
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qixs_id, wrfbdy%qixs, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qixe_id, wrfbdy%qixe, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qiys_id, wrfbdy%qiys, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qiye_id, wrfbdy%qiye, start = (/ 1, 1, 1, lngth /)))
      !-- qi tendencies on boundary
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qitxs_id, wrfbdy%qitxs, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qitxe_id, wrfbdy%qitxe, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qitys_id, wrfbdy%qitys, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qitye_id, wrfbdy%qitye, start = (/ 1, 1, 1, lngth /)))
   endif
   if(wrfbdy%n_moist > 4) then
      !-- qs on boundary
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qsxs_id, wrfbdy%qsxs, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qsxe_id, wrfbdy%qsxe, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qsys_id, wrfbdy%qsys, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qsye_id, wrfbdy%qsye, start = (/ 1, 1, 1, lngth /)))
      !-- qs tendencies on boundary
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qstxs_id, wrfbdy%qstxs, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qstxe_id, wrfbdy%qstxe, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qstys_id, wrfbdy%qstys, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qstye_id, wrfbdy%qstye, start = (/ 1, 1, 1, lngth /)))
   endif
   if(wrfbdy%n_moist > 5) then
      !-- qg on boundary
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qgxs_id, wrfbdy%qgxs, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qgxe_id, wrfbdy%qgxe, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qgys_id, wrfbdy%qgys, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qgye_id, wrfbdy%qgye, start = (/ 1, 1, 1, lngth /)))
      !-- qg tendencies on boundary
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qgtxs_id, wrfbdy%qgtxs, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qgtxe_id, wrfbdy%qgtxe, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qgtys_id, wrfbdy%qgtys, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qgtye_id, wrfbdy%qgtye, start = (/ 1, 1, 1, lngth /)))
   endif
   if(wrfbdy%n_moist > 6) then
      !-- qnice on boundary
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qnicexs_id, wrfbdy%qnicexs, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qnicexe_id, wrfbdy%qnicexe, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qniceys_id, wrfbdy%qniceys, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qniceye_id, wrfbdy%qniceye, start = (/ 1, 1, 1, lngth /)))
      !-- qnice tendencies on boundary
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qnicetxs_id, wrfbdy%qnicetxs, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qnicetxe_id, wrfbdy%qnicetxe, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qnicetys_id, wrfbdy%qnicetys, start = (/ 1, 1, 1, lngth /)))
      call check( nf90_get_var(wrfbdy%ncid, wrfbdy%qnicetye_id, wrfbdy%qnicetye, start = (/ 1, 1, 1, lngth /)))
   endif
   if(wrfbdy%n_moist > 7) then
      write(6,*) 'n_moist = ',wrfbdy%n_moist
      call error_handler(E_ERR,'wrfbdy_io', &
           'n_moist is too large.', source, revision, revdate)
   endif
endif

contains

  ! Internal subroutine - checks error status after each netcdf, prints 
  !                       text message each time an error code is returned. 
  subroutine check(istatus)
    integer, intent ( in) :: istatus 
    if(istatus /= nf90_noerr) call error_handler(E_ERR,'wrfbdy_io', &
         trim(nf90_strerror(istatus)), source, revision, revdate) 
  end subroutine check

end subroutine wrfbdy_io

!#######################################################

subroutine set_wrf_date (timestring, year, month, day, hour, minute, second)

implicit none

integer,           intent(in)  :: year, month, day, hour, minute, second
character(len=19), intent(out) :: timestring

character(len=4)  :: ch_year
character(len=2)  :: ch_month, ch_day, ch_hour, ch_minute, ch_second

if ( .not. module_initialized ) call initialize_module

write(ch_year,'(i4)') year
write(ch_month,'(i2)') month
if (ch_month(1:1) == " ") ch_month(1:1) = "0"
write(ch_day,'(i2)') day
if (ch_day(1:1) == " ") ch_day(1:1) = "0"
write(ch_hour,'(i2)') hour
if (ch_hour(1:1) == " ") ch_hour(1:1) = "0"
write(ch_minute,'(i2)') minute
if (ch_minute(1:1) == " ") ch_minute(1:1) = "0"
write(ch_second,'(i2)') second
if (ch_second(1:1) == " ") ch_second(1:1) = "0"

timestring(1:4)   = ch_year
timestring(5:5)   = "-"
timestring(6:7)   = ch_month
timestring(8:8)   = "-"
timestring(9:10)  = ch_day
timestring(11:11) = "_"
timestring(12:13) = ch_hour
timestring(14:14) = ":"
timestring(15:16) = ch_minute
timestring(17:17) = ":"
timestring(18:19) = ch_second

end subroutine set_wrf_date

!#######################################################

subroutine get_wrf_date (tstring, year, month, day, hour, minute, second)

implicit none
!--------------------------------------------------------
! Returns integers taken from wrf%timestring
! It is assumed that the tstring char array is as YYYY-MM-DD_hh:mm:ss

integer,           intent(out) :: year, month, day, hour, minute, second
character(len=19), intent(in)  :: tstring

if ( .not. module_initialized ) call initialize_module

read(tstring(1:4),'(i4)') year
read(tstring(6:7),'(i2)') month
read(tstring(9:10),'(i2)') day
read(tstring(12:13),'(i2)') hour
read(tstring(15:16),'(i2)') minute
read(tstring(18:19),'(i2)') second

end subroutine get_wrf_date

END MODULE wrf_data_module

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
