! DART software - Copyright 2004 - 2011 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

#include <misc.h>
#include <preproc.h>

module inicFileMod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!----------------------------------------------------------------------- 
! Purpose: 
! read and writes initial data netCDF history files
!
! Method: 
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------
! $Id$
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varder
  use clm_varmap  , only : begpatch, endpatch, numland, numpatch, &
                           landvec,  patchvec, begland, endland
  use clm_varpar  , only : nlevsno, nlevsoi, nlevlak, maxpatch, rtmlon, rtmlat
  use clm_varcon  , only : spval
  use fileutils   , only : getfil
#if (defined SPMD)
  use spmdMod     , only : masterproc, npes, compute_mpigs_patch, compute_mpigs_land
  use mpishorthand, only : mpir8, mpiint, mpilog, mpicom
#else
  use spmdMod     , only : masterproc
#endif
#if (defined RTM)
  use RtmMod      , only : volr
#endif
  implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

! netcdf data

  integer, private  :: ncid                !netCDF dataset id
  integer, private  :: dimid               !netCDF dimension id 
  integer, private  :: varid               !netCDF variable id

! input dataset dimensions

  integer, private  :: numland_dim         !value for [numland] from dataset
  integer, private  :: maxpatch_dim        !value for [maxpatch] from dataset
  integer, private  :: nlevsoi_dim         !value for [nlevsoi] from dataset
  integer, private  :: nlevsno_dim         !value for [nlevsno] from dataset
  integer, private  :: nlevtot_dim         !number of total (snow+soil) levels from dataset  
  integer, private  :: rtmlon_dim          !number of RTM longitudes
  integer, private  :: rtmlat_dim          !number of RTM latitudes

! methods

  public  :: do_inicwrite 
  private :: patch_to_land 
  private :: land_to_patch
  private :: set_init_filename

  INTERFACE patch_to_land
     MODULE procedure patch_to_land_1d_int
     MODULE procedure patch_to_land_1d_real
     MODULE procedure patch_to_land_2d_real
  END INTERFACE

  INTERFACE land_to_patch
     MODULE procedure land_to_patch_1d_int
     MODULE procedure land_to_patch_1d_real
     MODULE procedure land_to_patch_2d_real
  END INTERFACE

  SAVE

!=======================================================================
CONTAINS
!=======================================================================

  subroutine inicrd ()

!----------------------------------------------------------------------- 
! 
! Purpose: 
! read initial data from netCDF instantaneous initial data history file 
! for variables:
!   snlsno, dzsno, zsno, zisno, h2ocan, h2osno, snowdp, snowage, 
!   h2osoi_liq, h2osoi_ice, t_veg, t_grnd, t_soisno, t_lake
!
! Method: 
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------

    use shr_kind_mod, only: r8 => shr_kind_r8
    use clm_varctl, only : finidat
    implicit none
    include 'netcdf.inc'

! ------------------------ local variables -----------------------------
    integer :: i,j,k,l,m,n        !loop indices
    character(len=256) :: locfn   !local file name
    integer :: ndim               !input dimension       
    integer :: ret                !netcdf return code
#if (defined SPMD)
    integer :: numrecvv(0:npes-1) !vector of items to be received  
    integer :: displsv(0:npes-1)  !displacement vector
    integer :: numsend            !number of items to be sent
    integer :: ier                !mpi return code
#endif
    integer , allocatable :: ibuf1dl(:,:)
    integer , allocatable :: ibuf1dp(:)
    real(r8), allocatable :: rbuf1dl(:,:)
    real(r8), allocatable :: rbuf1dp(:)
    real(r8), allocatable :: rbuf2dl(:,:,:)
    real(r8), allocatable :: rbuf2dp(:,:)
! --------------------------------------------------------------------

! Open netCDF data file and read data

    if (masterproc) then

       call getfil (finidat, locfn, 0)
       call wrap_open (locfn, nf_nowrite, ncid)

! check input dimensions

       call wrap_inq_dimid (ncid, 'numland', dimid)
       call wrap_inq_dimlen (ncid, dimid, numland_dim)
       if (numland_dim /= numland) then
          write (6,*) 'INICRD error: numland values disagree'
          write (6,*) 'finidat numland = ',numland_dim,' model numland = ',numland
          call endrun
       end if

       call wrap_inq_dimid (ncid, 'maxpatch', dimid)
       call wrap_inq_dimlen (ncid, dimid, maxpatch_dim)
       if (maxpatch_dim /= maxpatch) then
          write (6,*) 'INICRD error: maxpatch values disagree'
          write (6,*) 'finidat maxpatch = ',maxpatch_dim,' model maxpatch = ',maxpatch
          call endrun
       end if

       call wrap_inq_dimid (ncid, 'nlevsno', dimid)
       call wrap_inq_dimlen (ncid, dimid, nlevsno_dim)
       if (nlevsno_dim /= nlevsno) then
          write (6,*) 'INICRD error: nlevsno values disagree'
          write (6,*) 'finidat levsno = ',nlevsno_dim,' model nlevsno = ',nlevsno
          call endrun
       end if

       call wrap_inq_dimid (ncid, 'nlevsoi', dimid)
       call wrap_inq_dimlen (ncid, dimid, nlevsoi_dim)
       if (nlevsoi_dim /= nlevsoi) then
          write (6,*) 'INICRD error: nlevsoi values disagree'
          write (6,*) 'finidat nlevsoi = ',nlevsoi_dim,' model nlevsoi = ',nlevsoi
          call endrun
       end if

       call wrap_inq_dimid (ncid, 'nlevtot', dimid)
       call wrap_inq_dimlen (ncid, dimid, nlevtot_dim)
       if (nlevtot_dim /= nlevsoi+nlevsno) then
          write (6,*) 'INICRD error: nlevtot values disagree'
          write (6,*) 'finidat nlevtot = ',nlevtot_dim,' model nlevtot = ',nlevsno+nlevsoi
          call endrun 
       end if

#if (defined RTM)
       ret = nf_inq_dimid (ncid, 'rtmlon', dimid)
       if (ret == NF_NOERR) then
          call wrap_inq_dimlen (ncid, dimid, rtmlon_dim)
          if (rtmlon_dim /= rtmlon) then
             write (6,*) 'INICRD error: rtmlon values disagree'
             write (6,*) 'finidat rtmlon = ',rtmlon_dim,' model rtmlon = ',rtmlon
             call endrun
          end if
       endif

       ret = nf_inq_dimid (ncid, 'rtmlat', dimid)
       if (ret == NF_NOERR) then
          call wrap_inq_dimlen (ncid, dimid, rtmlat_dim)
          if (rtmlat_dim /= rtmlat) then
             write (6,*) 'INICRD error: rtmlat values disagree'
             write (6,*) 'finidat rtmlat = ',rtmlat_dim,' model rtmlat = ',rtmlat
             call endrun
          end if
       endif
#endif

    endif ! if-masterproc block

! Obtain data - for the snow interfaces, are only examing the snow 
! interfaces above zi=0 which is why zisno and zsno have the same 
! level dimension below

    allocate (rbuf1dl(numland,maxpatch))
    allocate (ibuf1dl(numland,maxpatch))
    allocate (rbuf1dp(begpatch:endpatch))
    allocate (ibuf1dp(begpatch:endpatch))
    
    ! Read in zisno
    ! NOTE: zi(0) is set to 0 in routine iniTimeConst
    allocate (rbuf2dp(-nlevsno+0:-1,begpatch:endpatch))
    allocate (rbuf2dl(numland,maxpatch,-nlevsno+0:-1)) 
    if (masterproc) then
       call wrap_inq_varid (ncid, 'ZISNO_INI', varid)
       call wrap_get_var_realx(ncid, varid, rbuf2dl)
    endif
    call land_to_patch (rbuf2dl, rbuf2dp, nlevsno) 
    do k = begpatch,endpatch
       clm(k)%zi(-nlevsno+0:-1) = rbuf2dp(-nlevsno+0:-1,k)
    end do
    deallocate (rbuf2dl)
    deallocate (rbuf2dp)
    
    ! Read in zsno
    allocate (rbuf2dl(numland,maxpatch,-nlevsno+1:0))
    allocate (rbuf2dp(-nlevsno+1: 0,begpatch:endpatch))
    if (masterproc) then
       call wrap_inq_varid (ncid, 'ZSNO_INI', varid)
       call wrap_get_var_realx(ncid, varid, rbuf2dl)
    endif
    call land_to_patch (rbuf2dl, rbuf2dp, nlevsno)
    do k = begpatch,endpatch
       clm(k)%z (-nlevsno+1:0) = rbuf2dp(-nlevsno+1:0,k)
    end do
    deallocate (rbuf2dl)
    deallocate (rbuf2dp)
    
    ! Read in dzsno
    allocate (rbuf2dl(numland,maxpatch,-nlevsno+1:0))
    allocate (rbuf2dp(-nlevsno+1: 0,begpatch:endpatch))
    if (masterproc) then
       call wrap_inq_varid (ncid, 'DZSNO_INI', varid)
       call wrap_get_var_realx(ncid, varid, rbuf2dl)
    endif
    call land_to_patch (rbuf2dl, rbuf2dp, nlevsno)
    do k = begpatch,endpatch
       clm(k)%dz(-nlevsno+1:0) = rbuf2dp(-nlevsno+1:0,k)
    end do
    deallocate (rbuf2dl)
    deallocate (rbuf2dp)

    ! Read in h2osoi_liq
    allocate (rbuf2dl(numland,maxpatch,-nlevsno+1:nlevsoi)) 
    allocate (rbuf2dp(-nlevsno+1:nlevsoi,begpatch:endpatch))
    if (masterproc) then
       call wrap_inq_varid (ncid, 'H2OSOI_LIQ_INI', varid)
       call wrap_get_var_realx(ncid, varid, rbuf2dl)
    endif
    call land_to_patch (rbuf2dl, rbuf2dp, nlevsno+nlevsoi)
    do k = begpatch,endpatch
       clm(k)%h2osoi_liq(-nlevsno+1:nlevsoi) = rbuf2dp(-nlevsno+1:nlevsoi,k)
    end do
    deallocate (rbuf2dl)
    deallocate (rbuf2dp)
    
    ! Read in h2osoi_ice
    allocate (rbuf2dl(numland,maxpatch,-nlevsno+1:nlevsoi)) 
    allocate (rbuf2dp(-nlevsno+1:nlevsoi,begpatch:endpatch))
    if (masterproc) then
       call wrap_inq_varid (ncid, 'H2OSOI_ICE_INI', varid)
       call wrap_get_var_realx(ncid, varid, rbuf2dl)
    endif
    call land_to_patch (rbuf2dl, rbuf2dp, nlevsno+nlevsoi)
    do k = begpatch,endpatch
       clm(k)%h2osoi_ice(-nlevsno+1:nlevsoi) = rbuf2dp(-nlevsno+1:nlevsoi,k)
    end do
    deallocate (rbuf2dl)
    deallocate (rbuf2dp)
    
    ! Read in t_soisno
    allocate (rbuf2dl(numland,maxpatch,-nlevsno+1:nlevsoi)) 
    allocate (rbuf2dp(-nlevsno+1:nlevsoi,begpatch:endpatch))
    if (masterproc) then
       call wrap_inq_varid (ncid, 'T_SOISNO_INI', varid)
       call wrap_get_var_realx(ncid, varid, rbuf2dl)
    endif
    call land_to_patch (rbuf2dl, rbuf2dp, nlevsno+nlevsoi)
    do k = begpatch,endpatch
       clm(k)%t_soisno(-nlevsno+1:nlevsoi) = rbuf2dp(-nlevsno+1:nlevsoi,k)
    end do
    deallocate (rbuf2dl)
    deallocate (rbuf2dp)
    
    ! Read in t_lake
    allocate(rbuf2dl(numland,maxpatch,1:nlevlak))
    allocate(rbuf2dp(1:nlevlak,begpatch:endpatch))
    if (masterproc) then
       call wrap_inq_varid (ncid, 'T_LAKE_INI', varid)
       call wrap_get_var_realx(ncid, varid, rbuf2dl)
    endif
    call land_to_patch (rbuf2dl, rbuf2dp, nlevlak)
    do k = begpatch,endpatch
       clm(k)%t_lake(1:nlevlak) = rbuf2dp(1:nlevlak,k)
    end do
    deallocate (rbuf2dl)
    deallocate (rbuf2dp)
    
    ! Read in t_veg
    if (masterproc) then
       call wrap_inq_varid (ncid, 'T_VEG_INI', varid)
       call wrap_get_var_realx (ncid, varid, rbuf1dl)
    endif
    call land_to_patch (rbuf1dl, rbuf1dp)
    do k = begpatch,endpatch
       clm(k)%t_veg  = rbuf1dp(k)
    end do
    
    ! Read in t_grnd
    if (masterproc) then
       call wrap_inq_varid (ncid, 'T_GRND_INI', varid)
       call wrap_get_var_realx (ncid, varid, rbuf1dl)
    endif
    call land_to_patch (rbuf1dl, rbuf1dp)
    do k = begpatch,endpatch
       clm(k)%t_grnd = rbuf1dp(k)
    end do
    
    ! Read in h2ocan
    if (masterproc) then
       call wrap_inq_varid (ncid, 'H2OCAN_INI', varid)
       call wrap_get_var_realx(ncid, varid, rbuf1dl)
    endif
    call land_to_patch (rbuf1dl, rbuf1dp)
    do k = begpatch,endpatch
       clm(k)%h2ocan = rbuf1dp(k)
    end do

    ! Read in h2osno
    if (masterproc) then
       call wrap_inq_varid (ncid, 'H2OSNO_INI', varid)
       call wrap_get_var_realx(ncid, varid, rbuf1dl)
    endif
    call land_to_patch (rbuf1dl, rbuf1dp)
    do k = begpatch,endpatch
       clm(k)%h2osno = rbuf1dp(k)
    end do
    
    ! Read in snowdp
    if (masterproc) then
       call wrap_inq_varid (ncid, 'SNOWDP_INI', varid)
       call wrap_get_var_realx(ncid, varid, rbuf1dl)
    endif
    call land_to_patch (rbuf1dl, rbuf1dp)
    do k = begpatch,endpatch
       clm(k)%snowdp = rbuf1dp(k)
    end do
    
    ! Read in snowage
    if (masterproc) then
       call wrap_inq_varid (ncid, 'SNOWAGE_INI', varid)
       call wrap_get_var_realx(ncid, varid, rbuf1dl)
    endif
    call land_to_patch (rbuf1dl, rbuf1dp)
    do k = begpatch,endpatch
       clm(k)%snowage= rbuf1dp(k)
    end do
    
    ! Read in snlsno
    if (masterproc) then
       call wrap_inq_varid (ncid, 'SNLSNO_INI', varid)
       call wrap_get_var_int(ncid, varid, ibuf1dl)
    endif
    call land_to_patch (ibuf1dl, ibuf1dp)
    do k = begpatch,endpatch
       clm(k)%snl = ibuf1dp(k)
    end do

#if (defined RTM)

    if (masterproc) then
       ret = nf_inq_varid (ncid, 'RTMVOLR', varid)
       if (ret == NF_NOERR) then
          write(6,*)'INICFILE: reading in river volr'
          call wrap_get_var_realx(ncid, varid, volr)
       endif
    endif

#endif

    deallocate (ibuf1dl)
    deallocate (rbuf1dl)
    deallocate (ibuf1dp)
    deallocate (rbuf1dp)

    return
  end subroutine inicrd

!=======================================================================

  subroutine inicwrt ()

!----------------------------------------------------------------------- 
! Purpose: 
! write initial data to netCDF history file
!
! Method: 
! 
! Author: Mariana Vertenstein
!-----------------------------------------------------------------------

    use shr_kind_mod, only: r8 => shr_kind_r8
    use clm_varctl, only : caseid, ctitle, version, fsurdat
    use time_manager, only : get_nstep
    use fileutils, only : set_filename, putfil
    use clm_varctl, only : archive_dir, mss_wpass, mss_irt
    implicit none
    include 'netcdf.inc'

! ------------------------ local variables ---------------------------
    integer :: i,j,k,l,m,n                     !loop indices
    integer :: dim1_id(1)                      !netCDF dimension id for 1-d variables   
    integer :: dim2_id(2)                      !netCDF dimension id for 2-d variables
    integer :: dim3_id(3)                      !netCDF dimension id for 3-d variables
    integer :: omode                           !netCDF dummy variable
    integer :: status                          !netCDF error status
    character(len=256) :: loc_fn               !local 
    character(len=256) :: rem_dir              !remote (archive) directory
    character(len=256) :: rem_fn               !remote (archive) filename
    character(len=256) :: str                  !global attribute string 
    character(len=256) :: name                 !name of attribute
    character(len=256) :: unit                 !units of attribute
    character(len=256) :: mode                 !field mode (aver, inst, max, min, etc)
    character(len=  8) :: curdate              !current date
    character(len=  8) :: curtime              !current time 
    integer :: snlsno_id                       !netCDF variable id
    integer :: dzsno_id                        !netCDF variable id 
    integer :: zsno_id                         !netCDF variable id
    integer :: zisno_id                        !netCDF variable id
    integer :: h2osoi_liq_id                   !netCDF variable id 
    integer :: h2osoi_ice_id                   !netCDF variable id
    integer :: t_soisno_id                     !netCDF variable id  
    integer :: t_lake_id                       !netCDF variable id  
    integer :: t_veg_id                        !netCDF variable id
    integer :: t_grnd_id                       !netCDF variable id
    integer :: h2ocan_id                       !netCDF variable id
    integer :: h2osno_id                       !netCDF variable id  
    integer :: snowdp_id                       !netCDF variable id
    integer :: snowage_id                      !netCDF variable id    
#if (defined RTM)
    integer :: volr_id                         !netCDF variable id
#endif
    integer , allocatable :: ibuf1dl(:,:)
    integer , allocatable :: ibuf1dp(:)
    real(r8), allocatable :: rbuf1dl(:,:)
    real(r8), allocatable :: rbuf1dp(:)
    real(r8), allocatable :: rbuf2dl(:,:,:)
    real(r8), allocatable :: rbuf2dp(:,:)
! --------------------------------------------------------------------

! --------------------------------------------------------------------
! create initial conditions file for writing
! --------------------------------------------------------------------

    if (masterproc) then

       loc_fn = set_init_filename()
       write(6,*)
       write(6,*)'(INICFILEMOD): Writing clm2 initial conditions dataset at ',&
            trim(loc_fn), 'at nstep = ',get_nstep()
       write(6,*)

! create new netCDF file (in defined mode)

       call wrap_create (trim(loc_fn), nf_clobber, ncid)

! set fill mode to "no fill" to optimize performance

       status = nf_set_fill (ncid, nf_nofill, omode)
       if (status /= nf_noerr) then
          write (6,*) ' netCDF error = ',nf_strerror(status)
          call endrun
       end if

! define dimensions 

       call wrap_def_dim (ncid, 'numland'  , numland        ,numland_dim)
       call wrap_def_dim (ncid, 'maxpatch' , maxpatch       ,maxpatch_dim)
       call wrap_def_dim (ncid, 'nlevsno'  , nlevsno        ,nlevsno_dim)
       call wrap_def_dim (ncid, 'nlevsoi'  , nlevsoi        ,nlevsoi_dim)
       call wrap_def_dim (ncid, 'nlevtot'  , nlevsno+nlevsoi,nlevtot_dim)
#if (defined RTM)
       call wrap_def_dim (ncid, 'rtmlon'   , rtmlon         ,rtmlon_dim)
       call wrap_def_dim (ncid, 'rtmlat'   , rtmlat         ,rtmlat_dim)
#endif

! define global attributes

       str = 'CF1.0'
       call wrap_put_att_text (ncid, NF_GLOBAL, 'conventions', trim(str))
       
       call datetime(curdate, curtime)
       str = 'created on ' // curdate // ' ' // curtime
       call wrap_put_att_text (ncid, NF_GLOBAL,'history', trim(str))
       
       call getenv ('LOGNAME', str)
       call wrap_put_att_text (ncid, NF_GLOBAL, 'logname',trim(str))
       
       call getenv ('HOST', str)
       call wrap_put_att_text (ncid, NF_GLOBAL, 'host', trim(str))
       
       str = 'Community Land Model: CLM2'
       call wrap_put_att_text (ncid, NF_GLOBAL, 'source', trim(str))
       
       str = '$Name$' 
       call wrap_put_att_text (ncid, NF_GLOBAL, 'version', trim(str))
       
       str = '$Id$'
       call wrap_put_att_text (ncid, NF_GLOBAL, 'revision_id', trim(str))
       
       str = ctitle 
       call wrap_put_att_text (ncid,nf_global,'case_title',trim(str))

       str = caseid
       call wrap_put_att_text (ncid,nf_global,'case_id',trim(str))

! define current date

       mode = 'instantaneous initial conditions'

       name = 'current date as 8 digit integer (YYYYMMDD)'
       call wrap_def_var (ncid, 'mcdate', nf_int, 0, 0, varid)
       call wrap_put_att_text (ncid, varid, 'long_name',name)
       call wrap_put_att_text (ncid, varid, 'mode'     ,mode)

       name = 'current seconds of current date'
       unit = 's'
       call wrap_def_var (ncid, 'mcsec',  nf_int, 0, 0, varid)
       call wrap_put_att_text (ncid, varid, 'long_name',name)
       call wrap_put_att_text (ncid, varid, 'units'    ,unit)
       call wrap_put_att_text (ncid, varid, 'mode'     ,mode)

! define single-level fields (numland x maxpatch)

       dim2_id(1) = numland_dim; dim2_id(2) = maxpatch_dim
       name = 'vegetation temperature (T_VEG)'
       unit = 'K'
       call wrap_def_var (ncid, 'T_VEG_INI', nf_double, 2, dim2_id, t_veg_id)
       call wrap_put_att_text (ncid, t_veg_id, 'long_name',name)
       call wrap_put_att_text (ncid, t_veg_id, 'units'    ,unit)
       call wrap_put_att_text (ncid, t_veg_id, 'mode'     ,mode)

       dim2_id(1) = numland_dim; dim2_id(2) = maxpatch_dim
       name = 'ground temperature (T_GRND)'
       unit = 'K'
       call wrap_def_var (ncid, 'T_GRND_INI', nf_double, 2, dim2_id, t_grnd_id)
       call wrap_put_att_text (ncid, t_grnd_id, 'long_name',name)
       call wrap_put_att_text (ncid, t_grnd_id, 'units'    ,unit)
       call wrap_put_att_text (ncid, t_grnd_id, 'mode'     ,mode)

       dim2_id(1) = numland_dim; dim2_id(2) = maxpatch_dim
       name = 'canopy water (H2OCAN)'
       unit = 'kg/m2'
       call wrap_def_var (ncid, 'H2OCAN_INI', nf_double, 2, dim2_id, h2ocan_id)
       call wrap_put_att_text (ncid, h2ocan_id, 'long_name',name)
       call wrap_put_att_text (ncid, h2ocan_id, 'units'    ,unit)
       call wrap_put_att_text (ncid, h2ocan_id, 'mode'     ,mode)

       dim2_id(1) = numland_dim; dim2_id(2) = maxpatch_dim
       name = 'snow water (H2OSNO)'
       unit = 'kg/m2'
       call wrap_def_var (ncid, 'H2OSNO_INI', nf_double, 2, dim2_id, h2osno_id)
       call wrap_put_att_text (ncid, h2osno_id, 'long_name',name)
       call wrap_put_att_text (ncid, h2osno_id, 'units'    ,unit)
       call wrap_put_att_text (ncid, h2osno_id, 'mode'     ,mode)

       dim2_id(1) = numland_dim; dim2_id(2) = maxpatch_dim
       name = 'snow depth (SNOWDP)'
       unit = 'm'
       call wrap_def_var (ncid, 'SNOWDP_INI', nf_double, 2, dim2_id, snowdp_id)
       call wrap_put_att_text (ncid, snowdp_id, 'long_name',name)
       call wrap_put_att_text (ncid, snowdp_id, 'units'    ,unit)
       call wrap_put_att_text (ncid, snowdp_id, 'mode'     ,mode)

       dim2_id(1) = numland_dim; dim2_id(2) = maxpatch_dim
       name = 'snow age (SNOWAGE)'
       call wrap_def_var (ncid, 'SNOWAGE_INI', nf_double, 2, dim2_id, snowage_id)
       call wrap_put_att_text (ncid, snowage_id, 'long_name',name)
       call wrap_put_att_text (ncid, snowage_id, 'mode'     ,mode)

       dim2_id(1) = numland_dim; dim2_id(2) = maxpatch_dim
       name = 'number of snow layers (SNLSNO)'
       call wrap_def_var (ncid, 'SNLSNO_INI', nf_int, 2, dim2_id, snlsno_id)
       call wrap_put_att_text (ncid, snlsno_id, 'long_name',name)
       call wrap_put_att_text (ncid, snlsno_id, 'mode'     ,mode)

! define multi-level fields (numland x maxpatch x numlev)

       dim3_id(1) = numland_dim; dim3_id(2) = maxpatch_dim; dim3_id(3) = nlevtot_dim
       name = 'soil-snow temperature'
       unit = 'K'
       call wrap_def_var (ncid, 'T_SOISNO_INI', nf_double, 3, dim3_id, t_soisno_id)
       call wrap_put_att_text (ncid, t_soisno_id, 'long_name',name)
       call wrap_put_att_text (ncid, t_soisno_id, 'units'    ,unit)
       call wrap_put_att_text (ncid, t_soisno_id, 'mode'     ,mode)

       dim3_id(1) = numland_dim; dim3_id(2) = maxpatch_dim; dim3_id(3) = nlevsoi_dim
       name = 'lake temperature'
       unit = 'K'
       call wrap_def_var (ncid, 'T_LAKE_INI', nf_double, 3, dim3_id, t_lake_id)
       call wrap_put_att_text (ncid, t_lake_id, 'long_name',name)
       call wrap_put_att_text (ncid, t_lake_id, 'units'    ,unit)
       call wrap_put_att_text (ncid, t_lake_id, 'mode'     ,mode)

       dim3_id(1) = numland_dim; dim3_id(2) = maxpatch_dim; dim3_id(3) = nlevtot_dim
       name = 'liquid water'
       unit = 'kg/m2'
       call wrap_def_var (ncid, 'H2OSOI_LIQ_INI', nf_double, 3, dim3_id, h2osoi_liq_id)
       call wrap_put_att_text (ncid, h2osoi_liq_id, 'long_name',name)
       call wrap_put_att_text (ncid, h2osoi_liq_id, 'units'    ,unit)
       call wrap_put_att_text (ncid, h2osoi_liq_id, 'mode'     ,mode)

       dim3_id(1) = numland_dim; dim3_id(2) = maxpatch_dim; dim3_id(3) = nlevtot_dim
       name = 'ice lens'
       unit = 'kg/m2'                                                       
       call wrap_def_var (ncid, 'H2OSOI_ICE_INI', nf_double, 3, dim3_id, h2osoi_ice_id)
       call wrap_put_att_text (ncid, h2osoi_ice_id, 'long_name',name)
       call wrap_put_att_text (ncid, h2osoi_ice_id, 'units'    ,unit)
       call wrap_put_att_text (ncid, h2osoi_ice_id, 'mode'     ,mode)

       dim3_id(1) = numland_dim; dim3_id(2) = maxpatch_dim; dim3_id(3) = nlevsno_dim
       name = 'snow layer depth'
       unit = 'm'
       call wrap_def_var (ncid, 'ZSNO_INI', nf_double, 3, dim3_id, zsno_id)
       call wrap_put_att_text (ncid, zsno_id, 'long_name',name)
       call wrap_put_att_text (ncid, zsno_id, 'units'    ,unit)
       call wrap_put_att_text (ncid, zsno_id, 'mode'     ,mode)

       dim3_id(1) = numland_dim; dim3_id(2) = maxpatch_dim; dim3_id(3) = nlevsno_dim
       name = 'snow layer thickness'
       unit = 'm'
       call wrap_def_var (ncid, 'DZSNO_INI', nf_double, 3, dim3_id, dzsno_id)
       call wrap_put_att_text (ncid, dzsno_id, 'long_name',name)
       call wrap_put_att_text (ncid, dzsno_id, 'units'    ,unit)
       call wrap_put_att_text (ncid, dzsno_id, 'mode'     ,mode)

       dim3_id(1) = numland_dim; dim3_id(2) = maxpatch_dim; dim3_id(3) = nlevsno_dim
       name = 'snow interface depth'
       unit = 'm'
       call wrap_def_var (ncid, 'ZISNO_INI', nf_double, 3, dim3_id, zisno_id)
       call wrap_put_att_text (ncid, zisno_id, 'long_name',name)
       call wrap_put_att_text (ncid, zisno_id, 'units'    ,unit)
       call wrap_put_att_text (ncid, zisno_id, 'mode'     ,mode)

#if (defined RTM)
       dim2_id(1) = rtmlon_dim ; dim2_id(2) = rtmlat_dim

       name = 'water volumn in cell (volr)'
       unit = 'm3'
       call wrap_def_var (ncid, 'RTMVOLR', nf_double, 2, dim2_id, volr_id)
       call wrap_put_att_text (ncid, volr_id, 'long_name',name)
       call wrap_put_att_text (ncid, volr_id, 'units'    ,unit)
       call wrap_put_att_text (ncid, volr_id, 'mode'     ,mode)
#endif

! finish creating netCDF file (end define mode)

       status = nf_enddef(ncid)

    endif  !end of if-masterproc block

! --------------------------------------------------------------------
! Write single-level variables: [numland x maxpatch] and 
! multi-level variables: [numland x maxpatch x lev]
! NOTE: for non-spmd begpatch=1 and endpatch=numpatch
! --------------------------------------------------------------------

! Convert clm derived type components to patch vectors
! Convert initial data from subgrid patch form to land point form

    allocate (ibuf1dl(numland,maxpatch))
    allocate (rbuf1dl(numland,maxpatch))
    allocate (ibuf1dp(begpatch:endpatch))
    allocate (rbuf1dp(begpatch:endpatch))

    ! Write out zisno
    allocate (rbuf2dp(-nlevsno+0:-1,begpatch:endpatch))
    allocate (rbuf2dl(numland,maxpatch,-nlevsno+0:-1)) 
    do k = begpatch,endpatch
       rbuf2dp(-nlevsno+0:-1,k) = clm(k)%zi(-nlevsno+0:-1)  
    end do
    call patch_to_land (rbuf2dp, rbuf2dl, nlevsno)
    if (masterproc) call wrap_put_var_realx (ncid, zisno_id, rbuf2dl)
    deallocate (rbuf2dl)
    deallocate (rbuf2dp)

    ! Write out zsno
    allocate (rbuf2dl(numland,maxpatch,-nlevsno+1:0))
    allocate (rbuf2dp(-nlevsno+1: 0,begpatch:endpatch))
    do k = begpatch,endpatch
       rbuf2dp(-nlevsno+1:0,k) = clm(k)%z(-nlevsno+1:0)  
    end do
    call patch_to_land (rbuf2dp, rbuf2dl, nlevsno)
    if (masterproc) call wrap_put_var_realx (ncid, zsno_id, rbuf2dl)
    deallocate (rbuf2dl)
    deallocate (rbuf2dp)

    ! Write out dzsno
    allocate (rbuf2dl(numland,maxpatch,-nlevsno+1:0))
    allocate (rbuf2dp(-nlevsno+1: 0,begpatch:endpatch))
    do k = begpatch,endpatch
       rbuf2dp(-nlevsno+1: 0,k) = clm(k)%dz(-nlevsno+1: 0)  
    end do
    call patch_to_land (rbuf2dp, rbuf2dl, nlevsno)
    if (masterproc) call wrap_put_var_realx (ncid, dzsno_id, rbuf2dl)
    deallocate (rbuf2dl)
    deallocate (rbuf2dp)

    ! Write out h2osoi_liq
    allocate (rbuf2dl(numland,maxpatch,-nlevsno+1:nlevsoi)) 
    allocate (rbuf2dp(-nlevsno+1:nlevsoi,begpatch:endpatch))
    do k = begpatch,endpatch
       rbuf2dp(-nlevsno+1:nlevsoi,k) = clm(k)%h2osoi_liq(-nlevsno+1:nlevsoi) 
    end do
    call patch_to_land (rbuf2dp, rbuf2dl, nlevsno+nlevsoi)
    if (masterproc) call wrap_put_var_realx (ncid, h2osoi_liq_id, rbuf2dl)
    deallocate (rbuf2dl)
    deallocate (rbuf2dp)

    ! Write out h2osoi_ice
    allocate (rbuf2dl(numland,maxpatch,-nlevsno+1:nlevsoi)) 
    allocate (rbuf2dp(-nlevsno+1:nlevsoi,begpatch:endpatch))
    do k = begpatch,endpatch
       rbuf2dp(-nlevsno+1:nlevsoi,k) = clm(k)%h2osoi_ice(-nlevsno+1:nlevsoi) 
    end do
    call patch_to_land (rbuf2dp, rbuf2dl, nlevsno+nlevsoi)
    if (masterproc) call wrap_put_var_realx (ncid, h2osoi_ice_id, rbuf2dl)
    deallocate (rbuf2dl)
    deallocate (rbuf2dp)

    ! Write out t_soisno
    allocate (rbuf2dl(numland,maxpatch,-nlevsno+1:nlevsoi)) 
    allocate (rbuf2dp(-nlevsno+1:nlevsoi,begpatch:endpatch))
    do k = begpatch,endpatch
       rbuf2dp(-nlevsno+1:nlevsoi,k) = clm(k)%t_soisno(-nlevsno+1:nlevsoi) 
    end do
    call patch_to_land (rbuf2dp, rbuf2dl, nlevsno+nlevsoi)
    if (masterproc) call wrap_put_var_realx (ncid, t_soisno_id, rbuf2dl)
    deallocate (rbuf2dl)
    deallocate (rbuf2dp)

    ! Write out t_lake
    allocate(rbuf2dl(numland,maxpatch,1:nlevlak))
    allocate(rbuf2dp(1:nlevlak,begpatch:endpatch))
    do k = begpatch,endpatch
       rbuf2dp(1:nlevlak,k) = clm(k)%t_lake(1:nlevlak)
    end do
    call patch_to_land (rbuf2dp, rbuf2dl, nlevlak)
    if (masterproc) call wrap_put_var_realx (ncid, t_lake_id, rbuf2dl)
    deallocate (rbuf2dl)
    deallocate (rbuf2dp)

    ! Write out t_veg
    do k = begpatch,endpatch
       rbuf1dp(k) = clm(k)%t_veg   
    end do
    call patch_to_land (rbuf1dp, rbuf1dl)
    if (masterproc) call wrap_put_var_realx (ncid, t_veg_id, rbuf1dl)

    ! Write out t_grnd
    do k = begpatch,endpatch
       rbuf1dp(k) = clm(k)%t_grnd  
    end do
    call patch_to_land (rbuf1dp, rbuf1dl)
    if (masterproc) call wrap_put_var_realx (ncid, t_grnd_id, rbuf1dl)

    ! Write out h2ocan
    do k = begpatch,endpatch
       rbuf1dp(k) = clm(k)%h2ocan  
    end do
    call patch_to_land (rbuf1dp, rbuf1dl)
    if (masterproc) call wrap_put_var_realx (ncid, h2ocan_id, rbuf1dl)

    ! Write out h2osno
    do k = begpatch,endpatch
       rbuf1dp(k) = clm(k)%h2osno  
    end do
    call patch_to_land (rbuf1dp, rbuf1dl)
    if (masterproc) call wrap_put_var_realx (ncid, h2osno_id, rbuf1dl)

    ! Write out snowdp
    do k = begpatch,endpatch
       rbuf1dp(k) = clm(k)%snowdp
    end do
    call patch_to_land (rbuf1dp, rbuf1dl)
    if (masterproc) call wrap_put_var_realx (ncid, snowdp_id, rbuf1dl)

    ! Write out snowage
    do k = begpatch,endpatch
       rbuf1dp(k)= clm(k)%snowage
    end do
    call patch_to_land (rbuf1dp, rbuf1dl)
    if (masterproc) call wrap_put_var_realx (ncid, snowage_id, rbuf1dl)

    ! Write out snlsno
    do k = begpatch,endpatch
       ibuf1dp(k) = clm(k)%snl     
    end do
    call patch_to_land (ibuf1dp, ibuf1dl)
    if (masterproc) call wrap_put_var_int (ncid, snlsno_id, ibuf1dl)

#if (defined RTM)
    ! Write out volr
    if (masterproc) call wrap_put_var_realx (ncid, volr_id, volr)
#endif

    deallocate (ibuf1dl)
    deallocate (rbuf1dl)
    deallocate (ibuf1dp)
    deallocate (rbuf1dp)

! archive initial conditions dataset (Mass Store currently)

    if (masterproc) then
       call wrap_close (ncid)
       if (mss_irt > 0) then 
          rem_dir = trim(archive_dir) // '/init/'
          rem_fn = set_filename(rem_dir, loc_fn)
! kdr orig;          call putfil (loc_fn, rem_fn, mss_wpass, mss_irt, .true.)
          call putfil (loc_fn, rem_fn, mss_wpass, mss_irt, .false.)
       endif
    endif

    return
  end subroutine inicwrt

!=======================================================================
! BEGIN GENERIC PROCEDURE DEFINITIONS
!=======================================================================

  logical function do_inicwrite()

! kdr added last function to this list
    use time_manager, only : get_curr_date, get_prev_date, is_last_step

    use clm_varctl, only : hist_crtinic

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Determine if initial dataset is to be written at this time step
!
! Method: 
!
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------

! ------------------------ local variables ------------------------
    integer :: yr         !nstep year (0 -> ...)
    integer :: yrm1       !nstep-1 year (0 -> ...)
    integer :: daym1      !nstep-1 day (1 -> 31)
    integer :: day        !nstep day (1 -> 31)
    integer :: mon        !nstep month (1 -> 12)
    integer :: monm1      !nstep-1 month (1 -> 12)
    integer :: mcsec      !nstep time of day [seconds] 
    integer :: mcsecm1    !nstep-1 time of day [seconds]
! -----------------------------------------------------------------

    ! Set calendar for current time step and previous time step

    call get_curr_date (yr, mon, day, mcsec) 
    call get_prev_date (yrm1, monm1, daym1, mcsecm1)

    ! Determine if time to write out initial dataset

    do_inicwrite = .false.
    if (hist_crtinic /= 'NONE') then
       if (hist_crtinic == 'MONTHLY') then
          if (mon /= monm1 .and. monm1 /= -1) do_inicwrite = .true.
       else if (hist_crtinic == 'YEARLY') then
          if (monm1 == 12 .and. mon == 1)  do_inicwrite = .true.
! kdr
       else if (hist_crtinic == 'ENDOFRUN') then
          if (is_last_step())  do_inicwrite = .true.
       endif
    endif

  end function do_inicwrite

!=======================================================================

  character(len=256) function set_init_filename ()

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Determine initial dataset filenames
! 
! Method: 
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------

    use clm_varctl  , only : caseid
    use time_manager, only : get_curr_date

! ------------------------ local variables ------------------------
    character(len=256) :: cdate       !date char string
    integer :: day                    !day (1 -> 31)
    integer :: mon                    !month (1 -> 12)
    integer :: yr                     !year (0 -> ...)
    integer :: sec                    !seconds into current day
! -----------------------------------------------------------------

    call get_curr_date (yr, mon, day, sec) 
    write(cdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr,mon,day,sec
    set_init_filename = "./"//trim(caseid)//".clm2.i."//trim(cdate)//".nc"

  end function set_init_filename

!=======================================================================

!----------------------------------------------------------------------- 
! 
! Purpose: 
! [numland] x [maxpatch] array from 1d subgrid patches
!
! Method: 
! Map a subgrid input vector [fldin] of length [numpatch] to a 2-d
! [numland] x [maxpatch] output array [fldout]. Not all land points have
! [maxpatch] subgrid patches. Many have less. [numpatch] is some number <=
! [numland]*[maxpatch], i.e., is valid subgrid patches only. This routine
! converts a field from its [numpatch] representation to a [numland] x 
! [maxpatch] representation, setting values for non-valid subgrid patches 
! to that of the first valid subgrid patch using mapping from clm_map
! o numland  = number of land grid cells
! o maxpatch = maximum number of subgrid patches per grid cell
! o numpatch = actual number of subgrid patches (<= numland*maxpatch)
! 
!-----------------------------------------------------------------------

  subroutine patch_to_land_1d_int (fldin, fldout)
! ------------------------ arguments ---------------------------------
    integer, intent(in)  :: fldin(begpatch:endpatch)          
    integer, intent(out) :: fldout(numland,maxpatch) 
! --------------------------------------------------------------------

! ------------------------ local variables ----------------------
    integer l,m,k                   !indices
#if (defined SPMD)
    integer :: ier                   !MPI error status
    integer :: ibuf1d(numpatch)      !MPI temporary buffer 
    integer :: numrecvv(0:npes-1)    !vector of items to be received  
    integer :: displsv(0:npes-1)     !displacement vector
    integer :: numsend               !number of items to be sent
! ---------------------------------------------------------------
    call compute_mpigs_patch(1, numsend, numrecvv, displsv)
    call mpi_gatherv (fldin(begpatch), numsend , mpiint,  &
         ibuf1d, numrecvv, displsv, mpiint, 0, mpicom, ier)
    if (masterproc) then
       do m = 1, maxpatch           !subgrid patches for each land point
          do l = 1, numland         !land point index for [lsmlon] x [lsmlat] grid
             k = landvec%patch(l,m) !subgrid patch index: [1] to [numpatch]
             fldout(l,m) = ibuf1d(k)
          end do
       end do
    endif
#else 
    do m = 1, maxpatch              !subgrid patches for each land point
       do l = 1, numland            !land point index for [lsmlon] x [lsmlat] grid
          k = landvec%patch(l,m)    !subgrid patch index: [1] to [numpatch]
          fldout(l,m) = fldin(k)
       end do
    end do
#endif
    return
  end subroutine patch_to_land_1d_int

!=======================================================================

  subroutine patch_to_land_1d_real (fldin, fldout)
! ------------------------ arguments ---------------------------------
    real(r8), intent(in)  :: fldin(begpatch:endpatch)          
    real(r8), intent(out) :: fldout(numland,maxpatch) 
! --------------------------------------------------------------------
! ------------------------ local variables ----------------------
    integer l,m,k                   !indices
#if (defined SPMD)
    integer  :: ier                 !MPI error status
    real(r8) :: buf1d(numpatch)     !MPI temporary buffer 
    integer  :: numrecvv(0:npes-1)  !vector of items to be received  
    integer  :: displsv(0:npes-1)   !displacement vector
    integer  :: numsend             !number of items to be sent
! ---------------------------------------------------------------
    call compute_mpigs_patch(1, numsend, numrecvv, displsv)
    call mpi_gatherv (fldin(begpatch), numsend , mpir8, &
         buf1d, numrecvv, displsv, mpir8  , 0, mpicom, ier)
    if (masterproc) then
       do m = 1, maxpatch           !subgrid patches for each land point
          do l = 1, numland         !land point index for [lsmlon] x [lsmlat] grid
             k = landvec%patch(l,m) !subgrid patch index: [1] to [numpatch]
             fldout(l,m) = buf1d(k)
          end do
       end do
    endif
#else 
    do m = 1, maxpatch               !subgrid patches for each land point
       do l = 1, numland             !land point index for [lsmlon] x [lsmlat] grid
          k = landvec%patch(l,m)     !subgrid patch index: [1] to [numpatch]
          fldout(l,m) = fldin(k)
       end do
    end do
#endif
    return
  end subroutine patch_to_land_1d_real

!=======================================================================

  subroutine patch_to_land_2d_real (fldin, fldout, nlev)
! ------------------------ arguments ---------------------------------
    integer , intent(in)  :: nlev
    real(r8), intent(in)  :: fldin(nlev,begpatch:endpatch)
    real(r8), intent(out) :: fldout(numland,maxpatch,nlev) 
! --------------------------------------------------------------------
! ------------------------ local variables ----------------------
    integer l,m,k,n                    !indices
#if (defined SPMD)                    
    integer  :: ier                   !MPI error status
    real(r8) :: buf2d(nlev,numpatch)  !MPI temporary buffer
    integer  :: numrecvv(0:npes-1)    !vector of items to be received  
    integer  :: displsv(0:npes-1)     !displacement vector
    integer  :: numsend               !number of items to be sent
! ---------------------------------------------------------------
    call compute_mpigs_patch(nlev, numsend, numrecvv, displsv)
    call mpi_gatherv (fldin(1,begpatch), numsend , mpir8, &
         buf2d, numrecvv, displsv, mpir8, 0, mpicom, ier)
    if (masterproc) then
       do m = 1, maxpatch           !subgrid patches for each land point
          do l = 1, numland         !land point index for [lon] x [lat] grid
             k = landvec%patch(l,m) !subgrid patch index: [1] to [numpatch]
             do n = 1,nlev          !level index 
                fldout(l,m,n) = buf2d(n,k)
             end do
          end do
       end do
    endif
#else
    do m = 1,maxpatch               !subgrid patches for each land point
       do l = 1,numland             !land point index for [lon] x [lat] grid
          k = landvec%patch(l,m)    !subgrid patch index: [1] to [numpatch]
          do n = 1,nlev             !level index 
             fldout(l,m,n) = fldin(n,k)
          end do
       end do
    end do
#endif
    return
  end subroutine patch_to_land_2d_real

!=======================================================================

  subroutine land_to_patch_1d_int (fldin, fldout)
! ------------------------ arguments ---------------------------------
    integer, intent(in)  :: fldin(numland,maxpatch)
    integer, intent(out) :: fldout(begpatch:endpatch)
! --------------------------------------------------------------------
! ------------------------ local variables ----------------------
    integer l,m,k                   !indices
#if (defined SPMD)
    integer :: ier                  !MPI error status
    integer :: ibuf1d(numpatch)     !MPI temporary buffer
    integer :: numsendv(0:npes-1)   !vector of items to be sent
    integer :: displsv(0:npes-1)    !displacement vector
    integer :: numrecv              !number of items to be received
! ---------------------------------------------------------------
    if (masterproc) then
       do m = 1, maxpatch              !subgrid patches for each land point
          do l = 1, numland            !land point index 
             if (landvec%wtxy(l,m) > 0.) then
                k = landvec%patch(l,m) !subgrid patch index
                ibuf1d(k) = fldin(l,m)
             end if
          end do
       end do
    endif
    call compute_mpigs_patch(1, numrecv, numsendv, displsv)
    call mpi_scatterv (ibuf1d, numsendv, displsv, mpiint, &
         fldout(begpatch), numrecv , mpiint, 0, mpicom, ier)
#else
    do m = 1, maxpatch                !subgrid patches for each land point
       do l = 1, numland              !land point index 
          if (landvec%wtxy(l,m) > 0.) then
             k = landvec%patch(l,m)   !subgrid patch index
             fldout(k) = fldin(l,m)
          endif
       end do
    end do
#endif
    return
  end subroutine land_to_patch_1d_int

!=======================================================================

  subroutine land_to_patch_1d_real (fldin, fldout)
! ------------------------ arguments ---------------------------------
    real(r8), intent(in)  :: fldin(numland,maxpatch)
    real(r8), intent(out) :: fldout(begpatch:endpatch)
! --------------------------------------------------------------------
! ------------------------ local variables ----------------------
    integer l,m,k  !indices
#if (defined SPMD)
    integer  :: ier                  !MPI error status
    real(r8) :: buf1d(numpatch)      !MPI temporary buffer
    integer  :: numsendv(0:npes-1)   !vector of items to be sent
    integer  :: displsv(0:npes-1)    !displacement vector
    integer  :: numrecv              !number of items to be received
! ---------------------------------------------------------------
    if (masterproc) then
       do m = 1, maxpatch             !subgrid patches for each land point
          do l = 1, numland           !land point index 
             if (landvec%wtxy(l,m) > 0.) then
                k = landvec%patch(l,m) !subgrid patch index
                buf1d(k) = fldin(l,m)
             end if
          end do
       end do
    endif
    call compute_mpigs_patch(1, numrecv, numsendv, displsv)
    call mpi_scatterv (buf1d, numsendv, displsv, mpir8, &
         fldout(begpatch), numrecv , mpir8  , 0, mpicom, ier)
#else
    do m = 1, maxpatch              !subgrid patches for each land point
       do l = 1, numland            !land point index 
          if (landvec%wtxy(l,m) > 0.) then
             k = landvec%patch(l,m) !subgrid patch index
             fldout(k) = fldin(l,m)
          endif
       end do
    end do
#endif
    return
  end subroutine land_to_patch_1d_real

!=======================================================================

  subroutine land_to_patch_2d_real (fldin, fldout, nlev)
! ------------------------ arguments ---------------------------------
    integer , intent(in)  :: nlev
    real(r8), intent(in)  :: fldin (numland,maxpatch,nlev) 
    real(r8), intent(out) :: fldout(nlev,begpatch:endpatch)
! --------------------------------------------------------------------
! ------------------------ local variables ----------------------
    integer l,m,k,n                   !indices
#if (defined SPMD)                   
    integer  :: ier                   !MPI error status
    real(r8) :: buf2d(nlev,numpatch)  !MPI temporary buffer
    integer  :: numsendv(0:npes-1)    !vector of items to be sent
    integer  :: displsv(0:npes-1)     !displacement vector
    integer  :: numrecv               !number of items to be received
! ---------------------------------------------------------------
    if (masterproc) then
       do m = 1, maxpatch             !subgrid patches for each land point
          do l = 1, numland           !land point index 
             if (landvec%wtxy(l,m) > 0.) then
                k = landvec%patch(l,m) !subgrid patch index
                do n = 1,nlev
                   buf2d(n,k) = fldin(l,m,n)
                end do
             end if
          end do
       end do
    endif
    call compute_mpigs_patch(nlev, numrecv, numsendv, displsv)
    call mpi_scatterv (buf2d, numsendv, displsv, mpir8, &
         fldout(1,begpatch), numrecv , mpir8  , 0, mpicom, ier)
#else
    do m = 1, maxpatch               !subgrid patches for each land point
       do l = 1, numland             !land point index 
          k = landvec%patch(l,m)  !subgrid patch index
          if (landvec%wtxy(l,m) > 0.) then
             do n = 1,nlev
                fldout(n,k) = fldin(l,m,n)
             end do
          endif
       end do
    end do
#endif
    return
  end subroutine land_to_patch_2d_real

!=======================================================================
! END GENERIC PROCEDURE DEFINITIONS
!=======================================================================

end module inicFileMod







