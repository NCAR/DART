! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research

#include <misc.h>
#include <preproc.h>

module inicFileMod

! <next few lines under version control, do not edit>
! $URL: http://subversion.ucar.edu/DAReS/DART/trunk/models/cam/model_mod.f90 $
! $Id: model_mod.f90 2721 2007-03-27 00:08:01Z thoar $
! $Revision: 2721 $
! $Date: 2007-03-26 18:08:01 -0600 (Mon, 26 Mar 2007) $

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: inicFileMod
!
! !DESCRIPTION:
!
! Read and writes CLM initial data netCDF files
!
! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use clm_varpar     , only : nlevsno, nlevsoi, nlevlak, rtmlon, rtmlat
  use spmdMod        , only : masterproc, iam
  use shr_sys_mod    , only : shr_sys_flush
  use ncdio
  use abortutils     , only : endrun

! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: do_inicwrite     ! true=> write initial data file
  public :: inicfile         ! creat initial dataset
  public :: inicfields       ! read/write initial dataset fields
#ifdef COUP_CAM
  public :: inicperp         ! perpetual read
#endif
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 1/04
!
!EOP
!
! !PRIVATE METHODS:
  private :: set_init_filename  ! set initial filename
  character(len=256), private :: loc_finidat  ! local finidat name
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: inicfile
!
! !INTERFACE:
  subroutine inicfile(flag)
!
! !DESCRIPTION:
! Write instantaneous initial data to netCDF initial data file
!
! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
    use clmtype
    use clm_varctl  , only : finidat, caseid, ctitle, version, fsurdat, archive_dir, mss_wpass, mss_irt
    use time_manager, only : get_nstep, get_curr_date
    use fileutils   , only : set_filename, putfil, getfil
    use ncdio       , only : check_ret, check_dim
    use decompMod   , only : get_proc_global
#ifdef RTM
    use RunoffMod   , only : get_proc_rof_global
#endif
    use shr_sys_mod , only : shr_sys_getenv
!
! !ARGUMENTS:
    implicit none
    character(len=*) :: flag    ! flag to specify read or write
!
! !CALLED FROM:
! subroutine driver
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    character(len=256) :: loc_fn   ! local
    character(len=256) :: rem_dir  ! remote (archive) directory
    character(len=256) :: rem_fn   ! remote (archive) filename
    character(len=256) :: str      ! global attribute string
    character(len=  8) :: curdate  ! current date
    character(len=  8) :: curtime  ! current time
    integer :: ncid                ! netCDF dataset id
    integer :: dimid               ! netCDF dimension id
    integer :: varid               ! netCDF variable id
    integer :: ier                 ! error code
    integer :: omode               ! netCDF dummy variable
    integer :: numg                ! total number of gridcells across all processors
    integer :: numl                ! total number of landunits across all processors
    integer :: numc                ! total number of columns across all processors
    integer :: nump                ! total number of pfts across all processors
    integer :: nrof_lnd            ! total number of land runoff points across all procs
    integer :: nrof_ocn            ! total number of ocean runoff points across all procs
    character(len=32) :: subname='inicfile' ! subroutine name
!------------------------------------------------------------------------

    ! Get relevant sizes

    call get_proc_global(numg, numl, numc, nump)
#ifdef RTM
    call get_proc_rof_global(nrof_lnd, nrof_ocn)
#endif

    ! Read/Write initial file

    if (flag == 'write') then

       if (masterproc) then

          ! Create new netCDF file (in define mode) and set fill mode
          ! to "no fill" to optimize performance

          loc_fn = set_init_filename()
          write(6,*)
          write(6,*)'INICFILE: Writing clm initial conditions dataset at ',&
               trim(loc_fn), ' at nstep = ',get_nstep()
          write(6,*)
          call check_ret(nf_create(trim(loc_fn), nf_clobber, ncid), subname)
          call check_ret(nf_set_fill(ncid, nf_nofill, omode), subname)

          ! Define dimensions

          call check_ret( nf_def_dim(ncid, 'gridcell', numg           , dimid), subname )
          call check_ret( nf_def_dim(ncid, 'landunit', numl           , dimid), subname )
          call check_ret( nf_def_dim(ncid, 'column'  , numc           , dimid), subname )
          call check_ret( nf_def_dim(ncid, 'pft'     , nump           , dimid), subname )

          call check_ret( nf_def_dim(ncid, 'levsoi'  , nlevsoi        , dimid), subname )
          call check_ret( nf_def_dim(ncid, 'levlak'  , nlevlak        , dimid), subname )
          call check_ret( nf_def_dim(ncid, 'levsno'  , nlevsno        , dimid), subname )
          call check_ret( nf_def_dim(ncid, 'levtot'  , nlevsno+nlevsoi, dimid), subname )
#ifdef RTM
          call check_ret( nf_def_dim(ncid, 'ocnrof'  , nrof_ocn       , dimid), subname )
          call check_ret( nf_def_dim(ncid, 'lndrof'  , nrof_lnd       , dimid), subname )
          call check_ret( nf_def_dim(ncid, 'rtmlon'  , rtmlon         , dimid), subname )
          call check_ret( nf_def_dim(ncid, 'rtmlat'  , rtmlat         , dimid), subname )
#endif

          ! Define global attributes

          str = 'CF1.0'
          call check_ret(nf_put_att_text(ncid, NF_GLOBAL, 'conventions', len_trim(str), trim(str)), subname)

          call getdatetime(curdate, curtime)
          str = 'created on ' // curdate // ' ' // curtime
          call check_ret(nf_put_att_text(ncid, NF_GLOBAL, 'history', len_trim(str), trim(str)), subname)

          call shr_sys_getenv ('LOGNAME', str, ier)
          call check_ret(nf_put_att_text(ncid, NF_GLOBAL, 'logname', len_trim(str), trim(str)), subname)

          call shr_sys_getenv ('HOST', str, ier)
          call check_ret(nf_put_att_text(ncid, NF_GLOBAL, 'host', len_trim(str), trim(str)), subname)

          str = 'Community Land Model: CLM3'
          call check_ret(nf_put_att_text(ncid, NF_GLOBAL, 'source', len_trim(str), trim(str)), subname)
          
          str = '$Name: cam3_1_brnchT_release01 $' 
          call check_ret(nf_put_att_text(ncid, NF_GLOBAL, 'version', len_trim(str), trim(str)), subname)
          
          str = '$Id: inicFileMod.F90,v 1.10.10.18 2005/03/10 21:02:14 rosinski Exp $'
          call check_ret(nf_put_att_text(ncid, NF_GLOBAL, 'revision_id', len_trim(str), trim(str)), subname)

          str = ctitle
          call check_ret(nf_put_att_text(ncid, NF_GLOBAL, 'case_title', len_trim(str), trim(str)), subname)

          str = caseid
          call check_ret(nf_put_att_text(ncid, NF_GLOBAL, 'case_id', len_trim(str), trim(str)), subname)

          str = fsurdat
          call check_ret(nf_put_att_text(ncid, NF_GLOBAL, 'surface_dataset', len_trim(str), trim(str)), subname)

          ! Define variables

          call inicfields(flag='define', ncid=ncid)

          ! Leave define model

          call check_ret(nf_enddef(ncid), subname)

       end if

       call inicfields(flag='write' , ncid=ncid)

       ! Archive initial conditions dataset (Mass Store currently)

       if (masterproc) then
          call check_ret(nf_close(ncid), subname)
          if (mss_irt > 0) then
             rem_dir = trim(archive_dir) // '/init/'
             rem_fn = set_filename(rem_dir, loc_fn)
! kdr orig;             call putfil (loc_fn, rem_fn, mss_wpass, mss_irt, .true.)
          call putfil (loc_fn, rem_fn, mss_wpass, mss_irt, .false.)
          end if
       end if

    else if (flag == 'read') then

       if (masterproc) then
          call getfil(finidat, loc_finidat, 0)
          call check_ret(nf_open(loc_finidat, nf_nowrite, ncid), subname)
#if ( !defined SCAM )
          call check_dim(ncid, 'gridcell', numg)
          call check_dim(ncid, 'landunit', numl)
          call check_dim(ncid, 'column'  , numc)
          call check_dim(ncid, 'pft'     , nump)
          call check_dim(ncid, 'levsno'  , nlevsno)
          call check_dim(ncid, 'levsoi'  , nlevsoi)
          call check_dim(ncid, 'levlak'  , nlevlak) 
#endif
       end if

       call inicfields(flag='read', ncid=ncid)

       if (masterproc) then
          call check_ret(nf_close(ncid), subname)
       end if

    end if

  end subroutine inicfile

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: inicfields
!
! !INTERFACE:
  subroutine inicfields(flag, ncid)
!
! !DESCRIPTION:
! Read/Write initial data from/to netCDF instantaneous initial data file
!
! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
    use clmtype
    use time_manager, only : get_curr_date
    use decompMod   , only : get_proc_bounds, get_proc_global, map_dc2sn, &
                             get_sn_land1d, get_sn_cols1d, get_sn_pfts1d
    use clm_varcon  , only : sb, denice, denh2o
#ifdef DGVM
    use DGVMMod     , only : resetWeightsDGVM, gatherWeightsDGVM
#endif
#ifdef RTM
    use RunoffMod   , only : get_proc_rof_global
    use RtmMod      , only : volr, Rtmfluxini
#endif
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: flag  ! flag to determine if define, write or read data
    integer, intent(in) :: ncid           ! netCDF dataset id
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: g,l,c,p,j,i         ! indices
    integer :: yr                  ! current year (0 -> ...)
    integer :: mon                 ! current month (1 -> 12)
    integer :: day                 ! current day (1 -> 31)
    integer :: mcsec               ! seconds of current date
    integer :: mcdate              ! current date
    logical :: readvar             ! determine if variable is on initial file
    integer :: begp, endp          ! per-proc beginning and ending pft indices
    integer :: begc, endc          ! per-proc beginning and ending column indices
    integer :: begl, endl          ! per-proc beginning and ending landunit indices
    integer :: begg, endg          ! per-proc gridcell ending gridcell indices
    integer :: numg                ! total number of gridcells across all processors
    integer :: numl                ! total number of landunits across all processors
    integer :: numc                ! total number of columns across all processors
    integer :: nump                ! total number of pfts across all processors
    integer :: nrof_lnd            ! total number of land runoff points across all procs
    integer :: nrof_ocn            ! total number of ocean runoff points across all procs
    integer :: ier                 ! error status
    integer, pointer :: iltemp(:)  ! temporary
    integer, pointer :: ictemp(:)  ! temporary
    integer, pointer :: iptemp(:)  ! temporary
    integer, pointer :: ictype(:)  ! temporary
    integer, pointer :: iptype(:)  ! temporary
    type(gridcell_type), pointer :: gptr  ! pointer to gridcell derived subtype
    type(landunit_type), pointer :: lptr  ! pointer to landunit derived subtype
    type(column_type)  , pointer :: cptr  ! pointer to column derived subtype
    type(pft_type)     , pointer :: pptr  ! pointer to pft derived subtype
!------------------------------------------------------------------------

    ! Set pointers into derived type

    gptr => clm3%g
    lptr => clm3%g%l
    cptr => clm3%g%l%c
    pptr => clm3%g%l%c%p

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
    call get_proc_global(numg, numl, numc, nump)
#ifdef RTM
    call get_proc_rof_global(nrof_lnd, nrof_ocn)
#endif

    allocate(iltemp(numl), ictemp(numc), iptemp(nump), ictype(numc), iptype(nump), stat=ier)
    if (ier /= 0) then
       write(6,*)'allocation error from inicfields'; call endrun()
    end if

    ! Determine necessary indices
    ! Perform I/O

    if (flag == 'define') then

       ! Define time info

       call ncd_defvar(ncid=ncid, varname='mcdate', xtype=nf_int, &
            long_name='current date as 8 digit integer (YYYYMMDD)')

       call ncd_defvar(ncid=ncid, varname='mcsec', xtype=nf_int,  &
            long_name='current seconds of current date', units='s')

       ! Define gridcell info

       call ncd_defvar(ncid=ncid, varname='grid1d_lon', xtype=nf_double,  &
            dim1name='gridcell', long_name='gridcell longitude', units='degrees_east')

       call ncd_defvar(ncid=ncid, varname='grid1d_lat', xtype=nf_double,  &
            dim1name='gridcell', long_name='gridcell latitude', units='degrees_north')

       call ncd_defvar(ncid=ncid, varname='grid1d_ixy', xtype=nf_int,  &
            dim1name='gridcell', long_name='2d longitude index of corresponding gridcell')

       call ncd_defvar(ncid=ncid, varname='grid1d_jxy', xtype=nf_int,  &
            dim1name='gridcell', long_name='2d latitude index of corresponding gridcell')

       ! Define landunit info

       call ncd_defvar(ncid=ncid, varname='land1d_lon', xtype=nf_double,  &
            dim1name='landunit', long_name='landunit longitude', units='degrees_east')

       call ncd_defvar(ncid=ncid, varname='land1d_lat', xtype=nf_double,  &
            dim1name='landunit', long_name='landunit latitude', units='degrees_north')

       call ncd_defvar(ncid=ncid, varname='land1d_ixy', xtype=nf_int,  &
            dim1name='landunit', long_name='2d longitude index of corresponding landunit')

       call ncd_defvar(ncid=ncid, varname='land1d_jxy', xtype=nf_int,  &
            dim1name='landunit', long_name='2d latitude index of corresponding landunit')

       call ncd_defvar(ncid=ncid, varname='land1d_gi', xtype=nf_int,  &
            dim1name='landunit', long_name='1d grid index of corresponding landunit')

       call ncd_defvar(ncid=ncid, varname='land1d_wtxy', xtype=nf_double,  &
            dim1name='landunit', long_name='landunit weight relative to corresponding gridcell')

       call ncd_defvar(ncid=ncid, varname='land1d_ityplun', xtype=nf_int,  &
            dim1name='landunit', long_name='landunit type (vegetated,urban,lake,wetland or glacier)')

       ! Define column info

       call ncd_defvar(ncid=ncid, varname='cols1d_lon', xtype=nf_double,  &
            dim1name='column', long_name='column longitude', units='degrees_east')

       call ncd_defvar(ncid=ncid, varname='cols1d_lat', xtype=nf_double,  &
            dim1name='column', long_name='column latitude', units='degrees_north')

       call ncd_defvar(ncid=ncid, varname='cols1d_ixy', xtype=nf_int,   &
            dim1name='column', long_name='2d longitude index of corresponding column')

       call ncd_defvar(ncid=ncid, varname='cols1d_jxy', xtype=nf_int,   &
            dim1name='column', long_name='2d latitude index of corresponding column')

       call ncd_defvar(ncid=ncid, varname='cols1d_gi', xtype=nf_int,   &
            dim1name='column', long_name='1d grid index of corresponding column')

       call ncd_defvar(ncid=ncid, varname='cols1d_li', xtype=nf_int,   &
            dim1name='column', long_name='1d landunit index of corresponding column')

       call ncd_defvar(ncid=ncid, varname='cols1d_wtxy', xtype=nf_double,   &
            dim1name='column', long_name='column weight relative to corresponding gridcell')

       call ncd_defvar(ncid=ncid, varname='cols1d_wtlnd', xtype=nf_double,   &
            dim1name='column', long_name='column weight relative to corresponding landunit')

       call ncd_defvar(ncid=ncid, varname='cols1d_ityplun', xtype=nf_int,   &
            dim1name='column', long_name='column landunit type (vegetated,urban,lake,wetland or glacier)')

       ! Deine pft info

       call ncd_defvar(ncid=ncid, varname='pfts1d_lon', xtype=nf_double,  &
            dim1name='pft', long_name='pft longitude', units='degrees_east')

       call ncd_defvar(ncid=ncid, varname='pfts1d_lat', xtype=nf_double,  &
            dim1name='pft', long_name='pft latitude', units='degrees_north')

       call ncd_defvar(ncid=ncid, varname='pfts1d_ixy', xtype=nf_int,  &
            dim1name='pft', long_name='2d longitude index of corresponding pft')

       call ncd_defvar(ncid=ncid, varname='pfts1d_jxy', xtype=nf_int,  &
            dim1name='pft', long_name='2d latitude index of corresponding pft')

       call ncd_defvar(ncid=ncid, varname='pfts1d_gi', xtype=nf_int,  &
            dim1name='pft', long_name='1d grid index of corresponding pft')

       call ncd_defvar(ncid=ncid, varname='pfts1d_li', xtype=nf_int,  &
            dim1name='pft', long_name='1d landunit index of corresponding pft')

       call ncd_defvar(ncid=ncid, varname='pfts1d_ci', xtype=nf_int,  &
            dim1name='pft', long_name='1d column index of corresponding pft')

       call ncd_defvar(ncid=ncid, varname='pfts1d_wtxy', xtype=nf_double,  &
            dim1name='pft', long_name='pft weight relative to corresponding gridcell')

       call ncd_defvar(ncid=ncid, varname='pfts1d_wtlnd', xtype=nf_double,  &
            dim1name='pft', long_name='pft weight relative to corresponding landunit')

       call ncd_defvar(ncid=ncid, varname='pfts1d_wtcol', xtype=nf_double,  &
            dim1name='pft', long_name='pft weight relative to corresponding column')

       call ncd_defvar(ncid=ncid, varname='pfts1d_itypveg', xtype=nf_int,  &
            dim1name='pft', long_name='pft vegetation type')

       call ncd_defvar(ncid=ncid, varname='pfts1d_ityplun', xtype=nf_int,  &
            dim1name='pft', long_name='pft landunit type (vegetated,urban,lake,wetland or glacier)')

       ! Define initial file variables

       call ncd_defvar(ncid=ncid, varname='T_VEG', xtype=nf_double,  &
            dim1name='pft', long_name='vegetation temperature', units='K')

       call ncd_defvar(ncid=ncid, varname='T_GRND', xtype=nf_double,  &
            dim1name='column', long_name='ground temperature', units='K')

       call ncd_defvar(ncid=ncid, varname='EFLX_LWRAD_OUT', xtype=nf_double,  &
            dim1name='pft', long_name='emitted infrared (longwave) radiation', units='watt/m^2')

       call ncd_defvar(ncid=ncid, varname='H2OCAN', xtype=nf_double,  &
            dim1name='pft', long_name='canopy water', units='kg/m2')

       call ncd_defvar(ncid=ncid, varname='H2OSNO', xtype=nf_double,  &
            dim1name='column', long_name='snow water', units='kg/m2')

       call ncd_defvar(ncid=ncid, varname='SNOWDP', xtype=nf_double,  &
            dim1name='column', long_name='snow depth', units='m')

       call ncd_defvar(ncid=ncid, varname='SNOWAGE', xtype=nf_double,  &
            dim1name='column', long_name='snow age', units='unitless')

       call ncd_defvar(ncid=ncid, varname='SNLSNO', xtype=nf_int,  &
            dim1name='column', long_name='number of snow layers', units='unitless')

       call ncd_defvar(ncid=ncid, varname='T_SOISNO', xtype=nf_double,   &
            dim1name='column', dim2name='levtot', long_name='soil-snow temperature', units='K')

       call ncd_defvar(ncid=ncid, varname='T_LAKE', xtype=nf_double,  &
            dim1name='column', dim2name='levlak', long_name='lake temperature', units='K')

       call ncd_defvar(ncid=ncid, varname='H2OSOI_LIQ', xtype=nf_double,  &
            dim1name='column', dim2name='levtot', long_name='liquid water', units='kg/m2')

       call ncd_defvar(ncid=ncid, varname='H2OSOI_ICE', xtype=nf_double,   &
            dim1name='column', dim2name='levtot', long_name='ice lens', units='kg/m2')

       call ncd_defvar(ncid=ncid, varname='ZSNO', xtype=nf_double,  &
            dim1name='column', dim2name='levsno',long_name='snow layer depth', units='m')

       call ncd_defvar(ncid=ncid, varname='DZSNO', xtype=nf_double,  &
            dim1name='column', dim2name='levsno', long_name='snow layer thickness', units='m')

       call ncd_defvar(ncid=ncid, varname='ZISNO', xtype=nf_double,  &
            dim1name='column', dim2name='levsno', long_name='snow interface depth', units='m')

#ifdef RTM
       call ncd_defvar(ncid=ncid, varname='RTMVOLR', xtype=nf_double,  &
            dim1name='rtmlon', dim2name='rtmlat', &
            long_name='water volumn in cell (volr)', units='m3')
#endif

#ifdef DGVM

       call ncd_defvar(ncid=ncid, varname='ITYPVEG', xtype=nf_int,  &
            dim1name='pft', long_name='plant functional type')

       call ncd_defvar (ncid=ncid, varname='FPCGRID', xtype=nf_double,  &
            dim1name='pft', long_name='plant functional type cover - fraction of vegetated area')

       call ncd_defvar (ncid=ncid, varname='LAI_IND', xtype=nf_double,  &
            dim1name='pft', long_name='LAI per individual', units='m2/m2')

       call ncd_defvar(ncid=ncid, varname='CROWNAREA', xtype=nf_double,  &
            dim1name='pft', long_name='crownarea', units='m2')

       call ncd_defvar(ncid=ncid, varname='LITTERAG', xtype=nf_double,  &
            dim1name='pft', long_name='above ground litter', units='grams C/m2 gridcell area')

       call ncd_defvar(ncid=ncid, varname='LITTERBG', xtype=nf_double,  &
            dim1name='pft', long_name='below ground litter', units='grams C/m2 gridcell area')

       call ncd_defvar(ncid=ncid, varname='CPOOL_FAST', xtype=nf_double,  &
            dim1name='pft', long_name='fast carbon pool', units='grams C/m2 gridcell area')

       call ncd_defvar(ncid=ncid, varname='CPOOL_SLOW', xtype=nf_double,  &
            dim1name='pft', long_name='slow carbon pool', units='grams C/m2 gridcell area')

       call ncd_defvar(ncid=ncid, varname='PRESENT', xtype=nf_int,  &
            dim1name='pft', long_name='if pft is present is patch, 1=>yes,=>no)')

       call ncd_defvar(ncid=ncid, varname='NIND', xtype=nf_double,  &
            dim1name='pft', long_name='number of individuals',units= 'number of individuals/m2 gridcell area')

       call ncd_defvar(ncid=ncid, varname='LM_IND', xtype=nf_double,  &
            dim1name='pft',long_name='individual leaf mass',units='grams C')

       call ncd_defvar(ncid=ncid, varname='SM_IND', xtype=nf_double,  &
            dim1name='pft',long_name='individual stem mass', units='grams C')

       call ncd_defvar(ncid=ncid, varname='HM_IND', xtype=nf_double,  &
            dim1name='pft',long_name='individual heartwood mass',units='grams C')

       call ncd_defvar(ncid=ncid, varname='RM_IND', xtype=nf_double,  &
            dim1name='pft',long_name='individual root mass',units='grams C')

       call ncd_defvar(ncid=ncid, varname='HTOP', xtype=nf_double,  &
            dim1name='pft',long_name='canopy top',units='m')

#endif

    else if (flag == 'write') then

       ! Write output data (first write current date and seconds of current date)

       if (masterproc) then
          call get_curr_date (yr, mon, day, mcsec)
          mcdate = yr*10000 + mon*100 + day

          call ncd_ioglobal(varname='mcdate', data=mcdate, ncid=ncid, flag=flag, readvar=readvar)
          call ncd_ioglobal(varname='mcsec' , data=mcsec , ncid=ncid, flag=flag, readvar=readvar)
       end if

       ! Write gridcell info

       call ncd_iolocal(varname='grid1d_lon', data=gptr%londeg, datadim1='gridcell', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='grid1d_lat', data=gptr%latdeg, datadim1='gridcell', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='grid1d_ixy', data=gptr%ixy   , datadim1='gridcell', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='grid1d_jxy', data=gptr%jxy   , datadim1='gridcell', ncid=ncid, flag=flag)

       ! Write landunit info

       call ncd_iolocal(varname='land1d_lon'    , data=lptr%londeg , datadim1='landunit', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='land1d_lat'    , data=lptr%latdeg , datadim1='landunit', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='land1d_ixy'    , data=lptr%ixy    , datadim1='landunit', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='land1d_jxy'    , data=lptr%jxy    , datadim1='landunit', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='land1d_wtxy'   , data=lptr%wtgcell, datadim1='landunit', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='land1d_ityplun', data=lptr%itype  , datadim1='landunit', ncid=ncid, flag=flag)
       if (masterproc) then
          call get_sn_land1d(iltemp, type1d='gridcell')
          call ncd_ioglobal(varname='land1d_gi', data=iltemp, ncid=ncid, flag=flag)
       end if

       ! Write column info

       call ncd_iolocal(varname='cols1d_lon'  , data=cptr%londeg , datadim1='column', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='cols1d_lat'  , data=cptr%latdeg , datadim1='column', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='cols1d_ixy'  , data=cptr%ixy    , datadim1='column', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='cols1d_jxy'  , data=cptr%jxy    , datadim1='column', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='cols1d_wtxy' , data=cptr%wtgcell, datadim1='column', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='cols1d_wtlnd', data=cptr%wtlunit, datadim1='column', ncid=ncid, flag=flag)
       if (masterproc) then
          call get_sn_cols1d(ictemp, type1d='gridcell')
          call ncd_ioglobal(varname='cols1d_gi', data=ictemp, ncid=ncid, flag=flag)
          call get_sn_cols1d(ictemp, type1d='landunit')
          call ncd_ioglobal(varname='cols1d_li', data=ictemp, ncid=ncid, flag=flag)
          do c = 1,numc
             l = cptr%landunit(c)
             ictype(c) = lptr%itype(l)
          end do
          call map_dc2sn(ictype, ictemp, type1d=namec)
          call ncd_ioglobal(varname='cols1d_ityplun', data=ictemp, ncid=ncid, flag=flag)
       end if

       ! Write pft info

       call ncd_iolocal(varname='pfts1d_lon'    , data=pptr%londeg , datadim1='pft', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='pfts1d_lat'    , data=pptr%latdeg , datadim1='pft', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='pfts1d_ixy'    , data=pptr%ixy    , datadim1='pft', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='pfts1d_jxy'    , data=pptr%jxy    , datadim1='pft', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='pfts1d_wtxy'   , data=pptr%wtgcell, datadim1='pft', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='pfts1d_wtlnd'  , data=pptr%wtlunit, datadim1='pft', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='pfts1d_wtcol'  , data=pptr%wtcol  , datadim1='pft', ncid=ncid, flag=flag)
       call ncd_iolocal(varname='pfts1d_itypveg',data=pptr%itype   , datadim1='pft', ncid=ncid, flag=flag)
       if (masterproc) then
          call get_sn_pfts1d(iptemp, type1d='gridcell')
          call ncd_ioglobal(varname='pfts1d_gi', data=iptemp, ncid=ncid, flag=flag)
          call get_sn_pfts1d(iptemp, type1d='landunit')
          call ncd_ioglobal(varname='pfts1d_li', data=iptemp, ncid=ncid, flag=flag)
          call get_sn_pfts1d(iptemp, type1d='column')
          call ncd_ioglobal(varname='pfts1d_ci', data=iptemp, ncid=ncid, flag=flag)
          do p = 1,nump
             l = pptr%landunit(p)
             iptype(p) = lptr%itype(l)
          end do
          call map_dc2sn(iptype, iptemp, type1d=namep)
          call ncd_ioglobal(varname='pfts1d_ityplun', data=iptemp, ncid=ncid, flag=flag)
       end if

    end if ! end of if-write flag block

    if (flag == 'write' .or. flag == 'read') then

       ! For the snow interfaces, are only examing the snow interfaces
       ! above zi=0 which is why zisno and zsno have the same level dimension below
       ! (Note - for zisno, zi(0) is set to 0 in routine iniTimeConst)

       call ncd_iolocal(varname='T_VEG', data=pptr%pes%t_veg, datadim1='pft', &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) call endrun()

       call ncd_iolocal(varname='T_GRND', data=cptr%ces%t_grnd, datadim1='column', &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) call endrun()

       call ncd_iolocal(varname='EFLX_LWRAD_OUT', data=pptr%pef%eflx_lwrad_out, datadim1='pft', &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag == 'read' .and. .not. readvar) then
          write(6,*)'EFLX_LWRAD_OUT will be created at startup using t_grnd'
!dir$ concurrent
!cdir nodep
          do p = begp,endp
             c = pptr%column(p)
             pptr%pef%eflx_lwrad_out(p) = sb*(cptr%ces%t_grnd(c))**4
          end do
       end if

       call ncd_iolocal(varname='H2OCAN', data=pptr%pws%h2ocan, datadim1='pft', &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) call endrun()

       call ncd_iolocal(varname='H2OSNO', data=cptr%cws%h2osno, datadim1='column', &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) call endrun()

       call ncd_iolocal(varname='SNOWDP', data=cptr%cps%snowdp, datadim1='column', &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) call endrun()

       call ncd_iolocal(varname='SNOWAGE', data=cptr%cps%snowage, datadim1='column', &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) call endrun()

       call ncd_iolocal(varname='SNLSNO', data=cptr%cps%snl, datadim1='column', &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) call endrun()

       call ncd_iolocal(varname='T_SOISNO', data=cptr%ces%t_soisno, datadim1='column', datadim2='levtot', &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) call endrun()

       call ncd_iolocal(varname='T_LAKE', data=cptr%ces%t_lake, datadim1='column', datadim2='levlak', &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) call endrun()

       call ncd_iolocal(varname='H2OSOI_LIQ', data=cptr%cws%h2osoi_liq, datadim1='column', datadim2='levtot', &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) call endrun()

       call ncd_iolocal(varname='H2OSOI_ICE', data=cptr%cws%h2osoi_ice, datadim1='column', datadim2='levtot', &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) call endrun()

       call ncd_iolocal(varname='ZSNO', data=cptr%cps%z, datadim1='column', datadim2='levsno', &
            lowerb2=-nlevsno+1, upperb2=0, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) call endrun()

       call ncd_iolocal(varname='DZSNO', data=cptr%cps%dz, datadim1='column', datadim2='levsno', &
            lowerb2=-nlevsno+1, upperb2=0, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) call endrun()

       call ncd_iolocal(varname='ZISNO', data=cptr%cps%zi, datadim1='column', datadim2='levsno', &
            lowerb2=-nlevsno, upperb2=-1, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) call endrun()

#ifdef RTM

       ! If read field - also determine initial fluxout from volr

       call ncd_ioglobal(varname='RTMVOLR', data=volr, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) volr(:,:) = 0.
       if (flag=='read') call Rtmfluxini()

#endif


#ifdef DGVM

       call ncd_iolocal(varname='ITYPVEG', data=pptr%itype, &
            datadim1='pft', ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) call endrun()

       call ncd_iolocal (varname='FPCGRID', data=pptr%pdgvs%fpcgrid, &
            datadim1='pft', ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) call endrun()

       call ncd_iolocal (varname='LAI_IND', data=pptr%pdgvs%lai_ind, &
            datadim1='pft', ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) call endrun()

       call ncd_iolocal(varname='CROWNAREA', data=pptr%pdgvs%crownarea, &
            datadim1='pft', ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) call endrun()

       call ncd_iolocal(varname='LITTERAG', data=pptr%pdgvs%litterag, &
            datadim1='pft', ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) call endrun()

       call ncd_iolocal(varname='LITTERBG', data=pptr%pdgvs%litterbg, &
            datadim1='pft', ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) call endrun()

       call ncd_iolocal(varname='CPOOL_FAST', data=pptr%pdgvs%cpool_fast, &
            datadim1='pft', ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) call endrun()

       call ncd_iolocal(varname='CPOOL_SLOW', data=pptr%pdgvs%cpool_slow, &
            datadim1='pft', ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) call endrun()

!dir$ concurrent
!cdir nodep
       do p = begp,endp
          iptemp(p) = 0
          if (pptr%pdgvs%present(p)) iptemp(p) = 1
       end do
       call ncd_iolocal(varname='PRESENT', data=iptemp, &
            datadim1='pft', ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) call endrun()
!dir$ concurrent
!cdir nodep
       do p = begp,endp
          pptr%pdgvs%present(p) = .false.
          if (iptemp(p) == 1) pptr%pdgvs%present(p) = .true.
       end do

       call ncd_iolocal(varname='NIND', data=pptr%pdgvs%nind, &
            datadim1='pft', ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) call endrun()

       call ncd_iolocal(varname='LM_IND', data=pptr%pdgvs%lm_ind, &
            datadim1='pft', ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) call endrun()

       call ncd_iolocal(varname='SM_IND', data=pptr%pdgvs%sm_ind, &
            datadim1='pft', ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) call endrun()

       call ncd_iolocal(varname='HM_IND', data=pptr%pdgvs%hm_ind, &
            datadim1='pft', ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) call endrun()

       call ncd_iolocal(varname='RM_IND', data=pptr%pdgvs%rm_ind, &
            datadim1='pft', ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) call endrun()

       call ncd_iolocal(varname='HTOP', data=pptr%pps%htop, &
            datadim1='pft', ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) call endrun()

#endif

    end if   ! end of if-read or if-write block

    ! Upon read only - calculate the following

    if (flag == 'read') then

       ! Obtain column average of h2ocan (driver.F90 needs this upon input for an initial run)
       ! The following should not be vectorized

       do c = begc,endc
          cptr%cws%pws_a%h2ocan(c) = 0.
       end do
       do p = begp,endp
          c = pptr%column(p)
          cptr%cws%pws_a%h2ocan(c) = cptr%cws%pws_a%h2ocan(c) + pptr%pws%h2ocan(p) * pptr%wtcol(p)
       end do

       ! Determine volumetric soil water

       do j = 1,nlevsoi
!dir$ concurrent
!cdir nodep
          do c = begc,endc
             l = cptr%landunit(c)
             if (.not. lptr%lakpoi(l)) then
                cptr%cws%h2osoi_vol(c,j) = cptr%cws%h2osoi_liq(c,j)/(cptr%cps%dz(c,j)*denh2o) &
                                          +cptr%cws%h2osoi_ice(c,j)/(cptr%cps%dz(c,j)*denice)
             end if
          end do
       end do

#ifdef DGVM
       ! Determine new subgrid weights and areas (obtained from new value of fpcgrid read in above)

       call resetWeightsDGVM(begg, endg, begc, endc, begp, endp)
#ifdef SPMD
       call gatherWeightsDGVM()
#endif
#endif

    end if

    deallocate(iltemp, ictemp, iptemp, ictype, iptype)

  end subroutine inicfields

#ifdef COUP_CAM

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: inicperp
!
! !INTERFACE:
  subroutine inicperp
!
! !DESCRIPTION:
! Read perpetual initial data fields
!
! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
    use clmtype
    use clm_varctl  , only : finidat
    use clm_varcon  , only : zlnd, denice, denh2o
    use decompMod   , only : get_proc_bounds, get_proc_global
    use fileutils   , only : getfil
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ncid                           ! netCDF dataset id
    integer :: j,c,l                          ! indices
    integer :: begp, endp                     ! per-proc beginning and ending pft indices
    integer :: begc, endc                     ! per-proc beginning and ending column indices
    integer :: begl, endl                     ! per-proc beginning and ending landunit indices
    integer :: begg, endg                     ! per-proc gridcell ending gridcell indices
    integer :: numg                           ! total number of gridcells across all processors
    integer :: numl                           ! total number of landunits across all processors
    integer :: numc                           ! total number of columns across all processors
    integer :: nump                           ! total number of pfts across all processors
    logical :: readvar                        ! determine if variable is on initial file
    type(gridcell_type), pointer :: gptr      ! pointer to gridcell derived subtype
    type(landunit_type), pointer :: lptr      ! pointer to landunit derived subtype
    type(column_type)  , pointer :: cptr      ! pointer to column derived subtype
    type(pft_type)     , pointer :: pptr      ! pointer to pft derived subtype
    logical, save :: opened_finidat = .false. ! true => finidat was opened for read
    character(len=32) :: subname='inicperp'   ! subroutine name
!------------------------------------------------------------------------

    ! Set pointers into derived type

    gptr => clm3%g
    lptr => clm3%g%l
    cptr => clm3%g%l%c
    pptr => clm3%g%l%c%p

    ! Determine necessary processor subgrid bounds

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
    call get_proc_global(numg, numl, numc, nump)

    ! Read namelist

    if (masterproc) then
       if (.not. opened_finidat) then
          call getfil(finidat, loc_finidat, 0)
          call check_ret(nf_open(loc_finidat, nf_nowrite, ncid), subname)
	  write(6,*)'INICPERP: opened netcdf file ',loc_finidat
	  call shr_sys_flush(6)
#if ( !defined SCAM )
          call check_dim(ncid, 'gridcell', numg)
          call check_dim(ncid, 'landunit', numl)
          call check_dim(ncid, 'column'  , numc)
          call check_dim(ncid, 'pft'     , nump)
          call check_dim(ncid, 'levsno'  , nlevsno)
          call check_dim(ncid, 'levsoi'  , nlevsoi)
          call check_dim(ncid, 'levlak'  , nlevlak) 
#endif
          opened_finidat = .true.
       else
          call check_ret(nf_open(loc_finidat, nf_nowrite, ncid), subname)
       end if
    end if

    call ncd_iolocal(varname='ZSNO', data=cptr%cps%z, datadim1='column', datadim2='levsno', &
         lowerb2=-nlevsno+1, upperb2=0, ncid=ncid, flag='read', readvar=readvar)
    if (.not. readvar) call endrun()

    call ncd_iolocal(varname='DZSNO', data=cptr%cps%dz, datadim1='column', datadim2='levsno', &
         lowerb2=-nlevsno+1, upperb2=0, ncid=ncid, flag='read', readvar=readvar)
    if (.not. readvar) call endrun()

    call ncd_iolocal(varname='ZISNO', data=cptr%cps%zi, datadim1='column', datadim2='levsno', &
         lowerb2=-nlevsno, upperb2=-1, ncid=ncid, flag='read', readvar=readvar)
    if (.not. readvar) call endrun()

    call ncd_iolocal(varname='H2OSNO', data=cptr%cws%h2osno, datadim1='column', &
         ncid=ncid, flag='read', readvar=readvar)
    if (.not. readvar) call endrun()

    call ncd_iolocal(varname='SNOWDP', data=cptr%cps%snowdp, datadim1='column', &
         ncid=ncid, flag='read', readvar=readvar)
    if (.not. readvar) call endrun()

    call ncd_iolocal(varname='SNOWAGE', data=cptr%cps%snowage, datadim1='column', &
         ncid=ncid, flag='read', readvar=readvar)
    if (.not. readvar) call endrun()

    call ncd_iolocal(varname='SNLSNO', data=cptr%cps%snl, datadim1='column', &
         ncid=ncid, flag='read', readvar=readvar)
    if (.not. readvar) call endrun()

    call ncd_iolocal(varname='H2OSOI_LIQ', data=cptr%cws%h2osoi_liq, datadim1='column', datadim2='levtot', &
         ncid=ncid, flag='read', readvar=readvar)
    if (.not. readvar) call endrun()

    call ncd_iolocal(varname='H2OSOI_ICE', data=cptr%cws%h2osoi_ice, datadim1='column', datadim2='levtot', &
         ncid=ncid, flag='read', readvar=readvar)
    if (.not. readvar) call endrun()

    if (masterproc) then
       call check_ret(nf_close(ncid), subname)
    end if

    ! Determine volumetric soil water

    do j = 1,nlevsoi
!dir$ concurrent
!cdir nodep
       do c = begc,endc
          l = cptr%landunit(c)
          if (.not. lptr%lakpoi(l)) then
             cptr%cws%h2osoi_vol(c,j) = cptr%cws%h2osoi_liq(c,j)/(cptr%cps%dz(c,j)*denh2o) &
                                       +cptr%cws%h2osoi_ice(c,j)/(cptr%cps%dz(c,j)*denice)
          end if
       end do
    end do

!dir$ concurrent
!cdir nodep
     do c = begc,endc
        cptr%cps%frac_sno(c) = cptr%cps%snowdp(c) / (10.*zlnd + cptr%cps%snowdp(c))
     end do

  end subroutine inicperp

#endif

!=======================================================================
! BEGIN GENERIC PROCEDURE DEFINITIONS
!=======================================================================

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_inicwrite
!
! !INTERFACE:
  logical function do_inicwrite()
!
! !DESCRIPTION:
! Determine if initial dataset is to be written at this time step
!
! !USES:
! kdr added last function to this list
    use time_manager, only : get_curr_date, get_prev_date, get_step_size, is_last_step
    use clm_varctl  , only : hist_crtinic
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine driver
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein

!EOP
!
! !LOCAL VARIABLES:
    integer :: yr         !nstep year (0 -> ...)
    integer :: yrm1       !nstep-1 year (0 -> ...)
    integer :: daym1      !nstep-1 day (1 -> 31)
    integer :: day        !nstep day (1 -> 31)
    integer :: mon        !nstep month (1 -> 12)
    integer :: monm1      !nstep-1 month (1 -> 12)
    integer :: mcsec      !nstep time of day [seconds]
    integer :: mcsecm1    !nstep-1 time of day [seconds]
    integer :: mcsecp1    !nstep+1 time of day [seconds]
    integer :: dayp1      !nstep+1 day (1 -> 31)
    integer :: monp1      !nstep+1 month (1 -> 12)
    integer :: yrp1       !nstep+1 year (0 -> ...)
    integer :: dtime      !timestep size [seconds]
!-----------------------------------------------------------------------

    ! Set calendar for current, previous, and next time steps

    dtime = get_step_size()
    call get_curr_date (yr  , mon  , day  , mcsec  )
    call get_prev_date (yrm1, monm1, daym1, mcsecm1)
    call get_curr_date (yrp1, monp1, dayp1, mcsecp1, offset=dtime)

    ! Determine if time to write out initial dataset

    do_inicwrite = .false.
    if (hist_crtinic /= 'NONE') then
       if      (hist_crtinic == '6-HOURLY') then
          if (mod(mcsecp1,21600) == 0) do_inicwrite = .true.
       elseif  (hist_crtinic == 'DAILY') then
          if (day /= dayp1)  do_inicwrite = .true.
       else if (hist_crtinic == 'MONTHLY') then
          if (mon /= monp1)  do_inicwrite = .true.
       else if (hist_crtinic == 'YEARLY') then
          if (mon == 12 .and. monp1 == 1)  do_inicwrite = .true.
! kdr
       else if (hist_crtinic == 'ENDOFRUN') then
          if (is_last_step())  do_inicwrite = .true.
       endif
    endif

  end function do_inicwrite

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_init_filename
!
! !INTERFACE:
  character(len=256) function set_init_filename ()
!
! !DESCRIPTION:
! Determine initial dataset filenames
!
! !USES:
    use clm_varctl  , only : caseid
    use time_manager, only : get_curr_date, get_step_size
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine inicwrt in this module
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    character(len=256) :: cdate       !date char string
    integer :: day                    !day (1 -> 31)
    integer :: mon                    !month (1 -> 12)
    integer :: yr                     !year (0 -> ...)
    integer :: sec                    !seconds into current day
    integer :: dtime                  !timestep size [seconds]
!-----------------------------------------------------------------------

    dtime = get_step_size()
    call get_curr_date (yr, mon, day, sec, offset=dtime)
    write(cdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr,mon,day,sec
    set_init_filename = "./"//trim(caseid)//".clm2.i."//trim(cdate)//".nc"

  end function set_init_filename

end module inicFileMod
