!--------------------------------------------------------------------
! This program that converts the following emissions data into 
!    into WRF input data data files. The potential fields are:
!       1) biogenic emissions data
!--------------------------------------------------------------------

program map_megan2_emissions

   use area_mapper, only : xlong => lon, xlat => lat
   use bio_types
   use misc_definitions_module, only : PROJ_PS, PROJ_LATLON, PROJ_LATLON, PROJ_CASSINI

   implicit none

   INTEGER, parameter :: lower = 0
   INTEGER, parameter :: upper = 1
   INTEGER :: icnt
   INTEGER :: ids, ide, jds, jde
   INTEGER :: i, j, n
   INTEGER :: ii, jj
   integer :: ncid
   integer :: ngatts
   integer :: cell
   integer :: map_proj
   integer :: ierr, astat
   integer :: dimid, varid
   integer :: nlon_megan, nlat_megan
   integer :: domain, domains
   integer :: mnth, mnth_s, mnth_e, megan_month
   integer :: start_lai_mnth = 1
   integer :: end_lai_mnth   = 12
   integer :: xndx_megan(2)
   integer :: yndx_megan(2)
   integer, allocatable :: ix(:,:,:)                        ! index used by interpolation
   integer, allocatable :: jy(:,:,:)                        ! index used by interpolation

   real    :: missing_value
   real    :: scale_factor
   real    :: wrk_sum
   real    :: ds1, ds2
   real    :: xl, xu
   real    :: yl, yu, dy
   real    :: wrf_lon_min
   real    :: wrf_lon_max
   real    :: wrf_lat_min
   real    :: wrf_lat_max
   real    :: cen_lon
   real    :: cen_lat
   real    :: stand_lon
   real    :: truelat1
   real    :: truelat2
   real    :: loninc
   real    :: latinc
   real    :: knowni
   real    :: knownj
   real    :: dx
   real(8), allocatable :: xedge_megan(:)
   real(8), allocatable :: yedge_megan(:)
   real, allocatable :: megan_lons(:)
   real, allocatable :: megan_lats(:)
   real, allocatable :: msebio_isop(:,:)
   real, allocatable :: pftp_bt(:,:)
   real, allocatable :: pftp_nt(:,:)
   real, allocatable :: pftp_sb(:,:)
   real, allocatable :: pftp_hb(:,:)
   real, allocatable :: mlai(:,:,:)
   real, allocatable :: mtsa(:,:,:)
   real, allocatable :: mswdown(:,:,:)
   real, allocatable :: tmp3(:,:,:)
   real, allocatable :: ax(:,:,:)                        ! weight coef. all domain
   real, allocatable :: by(:,:,:)                        ! weight coef. all domain

   character(len=19)   :: proj_name(0:3) = (/ 'LATLON             ', 'LAMBERT            ', &
                                              'POLAR STEREOGRAPHIC', 'MERCATOR           ' /)
   CHARACTER (LEN=132) :: varname
   CHARACTER (LEN=132) :: filespec
   CHARACTER (LEN=80)  :: message
   CHARACTER (LEN=80)  :: attribute
   CHARACTER (LEN=80)  :: units_attribute
   CHARACTER (LEN=80)  :: description_attribute
   CHARACTER (LEN=80)  :: stagger_attribute
   CHARACTER (LEN=80)  :: coor_attribute
   CHARACTER (LEN=80)  :: memord_attribute
   CHARACTER (LEN=80)  :: inpname
   CHARACTER (LEN=80)  :: outpname
   CHARACTER (LEN=80)  :: wrf_dir
   CHARACTER (LEN=80)  :: megan_dir
   CHARACTER (LEN=19)  :: Times(1)
   CHARACTER (LEN=3)   :: num
   CHARACTER (LEN=3)   :: char_mnth(12)

   logical :: has_area_map
   logical :: new_grid

   type glb_att
     integer :: len
     integer :: type
     character(len=132)  :: name
     integer(1), pointer :: attr_byte(:)
     integer(2), pointer :: attr_short(:)
     integer, pointer    :: attr_int(:)
     real, pointer       :: attr_real(:)
     real(8), pointer    :: attr_dbl(:)
     character(len=256)  :: attr_char
   end type

   type(glb_att), allocatable :: attrs(:)

   namelist /control/ domains, start_lai_mnth, end_lai_mnth, &
                      wrf_dir, megan_dir

!---------------------------------------------------------------------
!	... include files
!---------------------------------------------------------------------
   include 'netcdf.inc'

   char_mnth(:) = (/ 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', &
                     'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec' /)
   wrf_dir   = '.'
   megan_dir = '.'

!-----------------------------------------------------------------
!     read control variables
!-----------------------------------------------------------------
   read(*,nml=control,iostat=ierr)
   if( ierr /= 0 ) then
     write(*,*) 'convert_emissions: failed to read namelist; error = ',ierr
     stop 'bio_emiss abort'
   endif
!-----------------------------------------------------------------
!     check namelist inputs
!-----------------------------------------------------------------
   if( domains < 1 ) then
     write(*,*) 'convert_emissions: domains must be >= 1'
     stop 'bio_emiss abort'
   endif
   if( start_lai_mnth < 1 .or. start_lai_mnth > 12 ) then
     write(*,*) 'convert_emissions: start month must be in set [1,12]'
     stop 'bio_emiss abort'
   endif
   if( end_lai_mnth < 1 .or. end_lai_mnth > 12 ) then
     write(*,*) 'convert_emissions: end month must be in set [1,12]'
     stop 'bio_emiss abort'
   endif
   if( end_lai_mnth < start_lai_mnth ) then
     write(*,*) 'convert_emissions: end month must >= start_month'
      stop 'bio_emiss abort'
   endif

!-----------------------------------------------------------------
!     loop over domains
!-----------------------------------------------------------------
domain_loop : &
   do domain = 1,domains

      write(*,*) ' '
      write(*,*) '========================================='
      write(*,*) 'Domain = ',domain
      write(*,*) '========================================='
      write(*,*) ' '

!-----------------------------------------------------------
!     ... read wrfinput file
!-----------------------------------------------------------
      call wrf_file

!-----------------------------------------------------------
!     ... read and interpolate megan datasets
!-----------------------------------------------------------
      mnth_s        = 1
      mnth_e        = 12
      scale_factor  = 1.
      inpname       = 'TAS.nc'
      varname       = 'TAS_AVE'
      missing_value = -32768.
      CALL  megan2_bioemiss
      inpname       = 'DSW.nc'
      varname       = 'DSW_AVE'
      CALL  megan2_bioemiss

      mnth_e        = 1
      scale_factor  = 1.e-3
      missing_value = -32768.
      do megan_month = start_lai_mnth,end_lai_mnth
         write(num,'(i3)') 100+megan_month
         inpname       = 'laiv2003' // num(2:3) // '_30sec.nc'
         varname       = 'LAI_for_' // char_mnth(megan_month) // '_2003_(m2_per_m2)'
         CALL  megan2_bioemiss
      end do

      scale_factor  = 1./68.
      inpname       = 'isoall200021_30sec.nc'
      varname       = 'isoprene'
      CALL  megan2_bioemiss
      scale_factor  = 1.
      inpname       = 'btr200121_30sec.nc'
      varname       = 'Broadleaf_tree_cover_fraction_for_year_2001_(m2_per_m2)'
      missing_value = -128.
      CALL  megan2_bioemiss
      inpname       = 'ntr200121_30sec.nc'
      varname       = 'Needleleaf_tree_cover_fraction_for_year_2001_(m2_per_m2)'
      CALL  megan2_bioemiss
      inpname       = 'hrb200121_30sec.nc'
      varname       = 'Herbaceous_vegetation_cover_fraction_for_year_2001_(m2_per_m2)'
      CALL  megan2_bioemiss
      inpname       = 'shr200121_30sec.nc'
      varname       = 'Shrub_cover_fraction_for_year_2001_(m2_per_m2)'
      CALL  megan2_bioemiss

!-----------------------------------------------------------
!     ... write wrfbiochemi_d<nn> file
!-----------------------------------------------------------
      write(*,*) 'map_megan2_emissions: Before write_bioemiss'
      CALL  write_bioemiss
      write(*,*) 'map_megan2_emissions: After write_bioemiss'
      write(*,*) ' '
      write(*,'(''map_megan2_emissions: '',i3,'' cached grid(s)'')') grid_cnt

!-----------------------------------------------------------
!     cleanup domain variables
!-----------------------------------------------------------
      deallocate( mlai, msebio_isop, &
                  pftp_bt, pftp_nt, pftp_hb, pftp_sb, &
                  mtsa, mswdown )
      if( allocated( xlong ) ) then
        deallocate( xlong )
      endif
      if( allocated( xlat ) ) then
        deallocate( xlat )
      endif
      do n = 1,grid_cnt
        if( associated( grid_specs(n)%lon ) ) then
          deallocate( grid_specs(n)%lon )
        endif
        if( associated( grid_specs(n)%lat ) ) then
          deallocate( grid_specs(n)%lat )
        endif
        if( associated( grid_specs(n)%model_area_type ) ) then
          do j = 1,jde
            do i = 1,ide
              if( associated( grid_specs(n)%model_area_type(i,j)%dcell_lon_ndx ) ) then
                deallocate( grid_specs(n)%model_area_type(i,j)%dcell_lon_ndx )
              endif
              if( associated( grid_specs(n)%model_area_type(i,j)%dcell_lat_ndx ) ) then
                deallocate( grid_specs(n)%model_area_type(i,j)%dcell_lat_ndx )
              endif
              if( associated( grid_specs(n)%model_area_type(i,j)%wght ) ) then
                deallocate( grid_specs(n)%model_area_type(i,j)%wght )
              endif
            end do
          end do
          deallocate( grid_specs(n)%model_area_type )
        endif
      end do
!-----------------------------------------------------------------------
!   ... deallocate global attributes
!-----------------------------------------------------------------------
      if( ngatts > 0 ) then
        call dealloc_bioemiss_glb_atts
        deallocate( attrs )
      endif
      grid_cnt = 0
   end do domain_loop

   CONTAINS

   subroutine wrf_file
!---------------------------------------------------------------------
!   read wrf file
!---------------------------------------------------------------------

   use area_mapper, only : proj_init
   use constants_module, only : rad_per_deg, earth_radius_m

!---------------------------------------------------------------------
!   local variables
!---------------------------------------------------------------------
   real, allocatable :: wrk(:,:)

   write(num,'(i3)') 100+domain
   inpname = 'wrfinput_d' // num(2:3)
   filespec = trim( wrf_dir ) // '/' // trim( inpname )
!---------------------------------------------------------------------
!   open wrf input file
!---------------------------------------------------------------------
   message = 'wrf_file: Failed to open ' // trim(inpname)
   call handle_ncerr( nf_open( trim(filespec), nf_noclobber, ncid ), message )       
!---------------------------------------------------------------------
!   get wrf dimesions
!---------------------------------------------------------------------
   message = 'Failed to get lon dimension id'
   call handle_ncerr( nf_inq_dimid( ncid, 'west_east', dimid ), message )
   message = 'Failed to get lon dimension'
   call handle_ncerr( nf_inq_dimlen( ncid, dimid, ide ), message )
   message = 'Failed to get lat dimension id'
   call handle_ncerr( nf_inq_dimid( ncid, 'south_north', dimid ), message )
   message = 'Failed to get lat dimension'
   call handle_ncerr( nf_inq_dimlen( ncid, dimid, jde ), message )
!---------------------------------------------------------------------
!   get wrf map projection variables
!---------------------------------------------------------------------
   message = 'Failed to get map_proj'
   call handle_ncerr( nf_get_att_int( ncid, nf_global, 'MAP_PROJ', map_proj ), message )
   write(*,*) ' '
   write(*,*) 'wrf_file: MAP_PROJ is ',trim(proj_name(map_proj))
   write(*,*) ' '
   if( map_proj == PROJ_LATLON .or. map_proj == PROJ_CASSINI ) then
     message = 'Failed to get XLONG variable id'
     call handle_ncerr( nf_inq_varid( ncid, 'XLONG', varid ), message )
     allocate( wrk(ide,jde),stat=astat )
     if( astat /= 0 ) then
       write(*,*) 'wrf_file: failed to allocate wrk variable; error = ',astat
       stop 'Alloc error'
     endif
     message = 'Failed to read XLONG variable'
     call handle_ncerr( nf_get_var_real( ncid, varid, wrk ), message )
     cen_lon = wrk(1,1)
     loninc = wrk(2,1) - cen_lon
     message = 'Failed to get XLAT variable id'
     call handle_ncerr( nf_inq_varid( ncid, 'XLAT', varid ), message )
     message = 'Failed to read XLAT variable'
     call handle_ncerr( nf_get_var_real( ncid, varid, wrk ), message )
     cen_lat = wrk(1,1)
     latinc = wrk(1,2) - cen_lat
     knowni = 1.
     knownj = 1.
     deallocate( wrk )
     map_proj = PROJ_LATLON
     dx = earth_radius_m * latinc * rad_per_deg
   else
     message = 'wrf_file: Failed to get cen_lon'
     call handle_ncerr( nf_get_att_real( ncid, nf_global, 'CEN_LON', cen_lon ), message )
     write(*,*) 'wrf_file: CEN_LON = ',cen_lon
     message = 'wrf_file: Failed to get cen_lat'
     call handle_ncerr( nf_get_att_real( ncid, nf_global, 'CEN_LAT', cen_lat ), message )
     write(*,*) 'wrf_file: CEN_LAT = ',cen_lat
     message = 'wrf_file: Failed to get stand_lon'
     call handle_ncerr( nf_get_att_real( ncid, nf_global, 'STAND_LON', stand_lon ), message )
     write(*,*) 'wrf_file: STAND_LON = ',stand_lon
     message = 'wrf_file: Failed to get truelat1'
     call handle_ncerr( nf_get_att_real( ncid, nf_global, 'TRUELAT1', truelat1 ), message )
     write(*,*) 'wrf_file: TRUELAT1 = ',truelat1
     message = 'wrf_file: Failed to get truelat2'
     call handle_ncerr( nf_get_att_real( ncid, nf_global, 'TRUELAT2', truelat2 ), message )
     write(*,*) 'wrf_file: TRUELAT2 = ',truelat2
     message = 'wrf_file: Failed to get dx'
     call handle_ncerr( nf_get_att_real( ncid, nf_global, 'DX', dx ), message )
     write(*,*) 'wrf_file: DX = ',dx
   endif

!---------------------------------------------------------------------
!   initialize map projection
!---------------------------------------------------------------------
   call proj_init( map_proj, cen_lon, cen_lat, truelat1, truelat2, &
                   stand_lon, loninc, latinc, knowni, knownj, &
                   dx, ide, jde )

   ids = 1
   jds = 1

   message = 'Failed to get Times id'
   call handle_ncerr( nf_inq_varid( ncid, 'Times', varid ), message )
   message = 'Failed to read Times'
   call handle_ncerr( nf_get_var_text( ncid, varid, Times ), message )

   write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
   write(*,*) 'wrf_file: time = ',trim(Times(1))
   write(*,*) 'wrf_file: grid dimensions'
   write(*,*) 'wrf_file: ids,ide,jds,jde'
   write(*,'(4i6)') ids,ide,jds,jde
   write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
!---------------------------------------------------------------------
!   get wrfinput_<domain> global attributes
!---------------------------------------------------------------------
   call get_bioemiss_glb_atts
!---------------------------------------------------------------------
!   close wrf file
!---------------------------------------------------------------------
   message = 'wrf_file: Failed to close ' // trim(inpname)
   call handle_ncerr( nf_close( ncid ), message )       
!---------------------------------------------------------------------
!   allocate final bioemission variables
!---------------------------------------------------------------------
   allocate( mlai(ide,jde,12),stat=astat ) 
   if( astat /= 0 ) then
     write(*,*) 'wrf_file: failed to allocate mlai; error = ',astat
     stop 'allocate failed'
   endif
   allocate( msebio_isop(ide,jde),stat=astat ) 
   if( astat /= 0 ) then
     write(*,*) 'wrf_file: failed to allocate msebio_isop; error = ',astat
     stop 'allocate failed'
   endif
   allocate( pftp_bt(ide,jde),stat=astat ) 
   if( astat /= 0 ) then
     write(*,*) 'wrf_file: failed to allocate pftp_bt; error = ',astat
     stop 'allocate failed'
   endif
   allocate( pftp_nt(ide,jde),stat=astat ) 
   if( astat /= 0 ) then
     write(*,*) 'wrf_file: failed to allocate pftp_nt; error = ',astat
     stop 'allocate failed'
   endif
   allocate( pftp_sb(ide,jde),stat=astat ) 
   if( astat /= 0 ) then
     write(*,*) 'wrf_file: failed to allocate pftp_sb; error = ',astat
     stop 'allocate failed'
   endif
   allocate( pftp_hb(ide,jde),stat=astat ) 
   if( astat /= 0 ) then
     write(*,*) 'wrf_file: failed to allocate pftp_hb; error = ',astat
     stop 'allocate failed'
   endif
   allocate( mtsa(ide,jde,12),stat=astat ) 
   if( astat /= 0 ) then
     write(*,*) 'wrf_file: failed to allocate mtsa; error = ',astat
     stop 'allocate failed'
   endif
   allocate( mswdown(ide,jde,12),stat=astat ) 
   if( astat /= 0 ) then
     write(*,*) 'wrf_file: failed to allocate mswdown; error = ',astat
     stop 'allocate failed'
   endif

   end subroutine wrf_file

   subroutine megan2_bioemiss
!---------------------------------------------------------------------
!   map megan dataset to wrf grid
!---------------------------------------------------------------------

   use area_mapper, only : area_interp
   use constants_module, only : rad_per_deg, earth_radius_m
   use bio_types

!---------------------------------------------------------------------
!   local variables
!---------------------------------------------------------------------
    integer :: il, iu, jl, ju
    integer :: n
    integer :: status
    real    :: wrf_lon, wrf_lat
    real    :: data_dx
    real, allocatable :: wrk_data(:,:)
    logical :: debug = .false.

    write(*,*) ' '
    write(*,*) 'Reading megan2 bio emiss file ' // trim(inpname)
!---------------------------------------------------------------------
!   open megan dataset file
!---------------------------------------------------------------------
    message = 'megan2_bioemiss: Failed to open ' // trim(inpname)
    filespec = trim( megan_dir ) // '/' // trim(inpname)
    call handle_ncerr( nf_open( trim(filespec), nf_noclobber, ncid ), message )       
!---------------------------------------------------------------------
!   get megan dataset dimesions
!---------------------------------------------------------------------
    message = 'megan2_bioemiss: Failed to get lon dimension id'
    call handle_ncerr( nf_inq_dimid( ncid, 'lon', dimid ), message )
    message = 'megan2_bioemiss: Failed to get lon dimension'
    call handle_ncerr( nf_inq_dimlen( ncid, dimid, nlon_megan ), message )
    message = 'megan2_bioemiss: Failed to get lat dimension id'
    call handle_ncerr( nf_inq_dimid( ncid, 'lat', dimid ), message )
    message = 'megan2_bioemiss: Failed to get lat dimension'
    call handle_ncerr( nf_inq_dimlen( ncid, dimid, nlat_megan ), message )
    write(*,*) 'megan2_bioemiss:  nlon_megan, nlat_megan = ',nlon_megan,nlat_megan

!---------------------------------------------------------------------
!   allocate working variable
!---------------------------------------------------------------------
    if( allocated( megan_lons ) ) then
      deallocate( megan_lons )
    endif
    allocate( megan_lons(nlon_megan),stat=astat )
    if( astat /= 0 ) then
      write(*,*) 'megan2_bioemiss: Failed to allocate megan_lons; error = ',astat
      stop 'allocate failed'
    endif
!---------------------------------------------------------------------
!   read megan longitude variable
!---------------------------------------------------------------------
    message = 'megan2_bioemiss: Failed to get lon variable id'
    call handle_ncerr( nf_inq_varid( ncid, 'lon', varid ), message )
    message = 'megan2_bioemiss: Failed to read lon variable'
    call handle_ncerr( nf_get_var_real( ncid, varid, megan_lons ), message )

    if( allocated( megan_lats ) ) then
      deallocate( megan_lats )
    endif
    allocate( megan_lats(nlat_megan),stat=ierr )
    if( ierr /= 0 ) then
      write(*,*) 'megan2_bioemiss: Failed to allocate megan_lats; error = ',ierr
      stop 'allocate failed'
    endif
!---------------------------------------------------------------------
!   read megan latitude variable
!---------------------------------------------------------------------
    message = 'megan2_bioemiss: Failed to get lat variable id'
    call handle_ncerr( nf_inq_varid( ncid, 'lat', varid ), message )
    message = 'megan2_bioemiss: Failed to read lat variable'
    call handle_ncerr( nf_get_var_real( ncid, varid, megan_lats ), message )

!---------------------------------------------------------------------
!   determine interpolation type; bilinear or area conserving
!---------------------------------------------------------------------
    data_dx = earth_radius_m * (megan_lats(2) - megan_lats(1)) * rad_per_deg
    has_area_map = data_dx < dx
    write(*,*) 'megan2_bioemiss: data_dx,dx,has_area_map = ',data_dx,dx,has_area_map

!-------------------------------------------------------------
!   check for match against prior datasets
!-------------------------------------------------------------
   if( grid_cnt >= grid_max ) then
     write(*,*) 'megan2_bioemiss: reached grid cache max: ',grid_max
     stop
   endif
   grid_ndx = 0
   new_grid = .true.
   do n = 1,grid_cnt
     if( grid_specs(n)%nlons /= nlon_megan .or. grid_specs(n)%nlats /= nlat_megan ) then
       cycle
     endif
     if( any( grid_specs(n)%lon(:) /= megan_lons(:) ) ) then
       cycle
     endif
     if( any( grid_specs(n)%lat(:) /= megan_lats(:) ) ) then
       cycle
     endif
     grid_ndx = n
     new_grid = .false.
     exit
   end do
!-------------------------------------------------------------
!   new data grid to cache
!-------------------------------------------------------------
   if( new_grid ) then
     grid_cnt = grid_cnt + 1
     grid_specs(grid_cnt)%nlons = nlon_megan
     grid_specs(grid_cnt)%nlats = nlat_megan
     grid_specs(grid_cnt)%has_area_map = has_area_map
     allocate( grid_specs(grid_cnt)%lon(nlon_megan),stat=ierr )
     if( ierr /= 0 ) then
       write(*,*) 'megan2_bioemiss: Failed to allocate megan_lats; error = ',ierr
       stop 'allocate failed'
     endif
     allocate( grid_specs(grid_cnt)%lat(nlat_megan),stat=ierr )
     if( ierr /= 0 ) then
       write(*,*) 'megan2_bioemiss: Failed to allocate megan_lats; error = ',ierr
       stop 'allocate failed'
     endif
     grid_specs(grid_cnt)%lon(:) = megan_lons(:)
     grid_specs(grid_cnt)%lat(:) = megan_lats(:)
     if( has_area_map ) then
       allocate( grid_specs(grid_cnt)%model_area_type(ide,jde),stat=astat )
       if( astat /= 0 ) then
         write(*,*) 'proj_init; failed to allocate model_area_type: error = ',astat
         stop
       endif
       grid_specs(grid_cnt)%model_area_type(:,:)%has_data = .false.
       grid_specs(grid_cnt)%model_area_type(:,:)%active_dcell_cnt = 0
       grid_specs(grid_cnt)%model_area_type(:,:)%total_dcell_cnt  = 0
       grid_specs(grid_cnt)%model_area_type(:,:)%interior_dcell_cnt = 0
       grid_specs(grid_cnt)%model_area_type(:,:)%partial_dcell_cnt  = 0
     endif
     grid_ndx = grid_cnt
     write(*,*) 'megan2_bioemiss: file ' // trim(inpname),' has a new grid'
   endif

is_area_map : &
    if( has_area_map ) then
!---------------------------------------------------------------------
!   form megan longitude edges
!---------------------------------------------------------------------
      allocate( xedge_megan(nlon_megan+1),stat=astat )
      if( astat /= 0 ) then
        write(*,*) 'megan2_bioemiss: Failed to allocate xedge_megan; error = ',astat
        stop 'allocate error'
      endif
      xedge_megan(2:nlon_megan) = .5_8*(megan_lons(1:nlon_megan-1) + megan_lons(2:nlon_megan))
      xedge_megan(1)            = megan_lons(1) - .5_8*(megan_lons(2) - megan_lons(1))
      xedge_megan(nlon_megan+1) = megan_lons(nlon_megan) + .5_8*(megan_lons(nlon_megan) - megan_lons(nlon_megan-1))
      write(*,'(''megan2_bioemiss: xcen_megan(1,2)  = '',1p,2g22.15)') megan_lons(1:2)
      write(*,'(''megan2_bioemiss: xedge_megan(1,2) = '',1p,2g22.15)') xedge_megan(1:2)
      write(*,'(''megan2_bioemiss: dx = '',1pg22.15)') int( 1./(megan_lons(2) - megan_lons(1)) )
!---------------------------------------------------------------------
!   form megan latitude edges
!---------------------------------------------------------------------
      allocate( yedge_megan(nlat_megan+1),stat=ierr )
      if( ierr /= 0 ) then
        write(*,*) 'megan2_bioemiss: Failed to allocate yedge_megan; error = ',ierr
        stop 'allocate error'
      endif

      yedge_megan(2:nlat_megan) = .5_8*(megan_lats(1:nlat_megan-1) + megan_lats(2:nlat_megan))
      yedge_megan(1)            = megan_lats(1) - .5_8*(megan_lats(2) - megan_lats(1))
      yedge_megan(nlat_megan+1) = megan_lats(nlat_megan) + .5_8*(megan_lats(nlat_megan) - megan_lats(nlat_megan-1))

      write(*,'(''megan2_bioemiss: nlon_megan,nlat_megan = '',i6,1x,i6)') nlon_megan,nlat_megan
      write(*,'(''megan2_bioemiss: ycen_megan  = '',1p,2g22.15)') megan_lats(nlat_megan-1:nlat_megan)
      write(*,'(''megan2_bioemiss: yedge_megan = '',1p,2g22.15)') yedge_megan(nlat_megan:nlat_megan+1)

      if( allocated( wrk_data ) ) then
        deallocate( wrk_data )
      endif
      allocate( wrk_data(ide,jde),stat=status )
      if( status /= 0 ) then
        write(*,*) 'megan2_bioemiss: allocate for wrk_data failed; error = ',ierr
        stop 'allocation error'
      endif
!---------------------------------------------------------------------
!   area conserving interpolation
!---------------------------------------------------------------------
      call area_interp( xedge_megan, yedge_megan, nlon_megan, nlat_megan, int(missing_value,2), &
                        wrk_data, ncid, varname, grid_ndx, new_grid )
      if( varname(:3) == 'LAI' ) then
        mlai(:,:,megan_month) = scale_factor * wrk_data(:,:)
      elseif( trim(varname) == 'isoprene' ) then
        msebio_isop(:,:) = scale_factor * wrk_data(:,:)
      elseif( varname(:14) == 'Broadleaf_tree' ) then
        pftp_bt(:,:) = scale_factor * wrk_data(:,:)
      elseif( varname(:15) == 'Needleleaf_tree' ) then
        pftp_nt(:,:) = scale_factor * wrk_data(:,:)
      elseif( varname(:5) == 'Shrub' ) then
        pftp_sb(:,:) = scale_factor * wrk_data(:,:)
      elseif( varname(:10) == 'Herbaceous' ) then
        pftp_hb(:,:) = scale_factor * wrk_data(:,:)
      endif
      deallocate( wrk_data )
    else is_area_map
!---------------------------------------------------------------------
!   form megan coordinate limits
!---------------------------------------------------------------------
      wrf_lon_min = minval( xlong(ids:ide,jds:jde) )
      wrf_lon_max = maxval( xlong(ids:ide,jds:jde) )
      wrf_lat_min = minval( xlat(ids:ide,jds:jde) )
      wrf_lat_max = maxval( xlat(ids:ide,jds:jde) )
      write(*,*) ' '
      write(*,'('' megan2_bioemiss: model lon limits = '',1p,2g25.16)') wrf_lon_min,wrf_lon_max
      write(*,'('' megan2_bioemiss: model lat limits = '',1p,2g25.16)') wrf_lat_min,wrf_lat_max
      write(*,*) ' '
      write(*,'('' megan2_bioemiss: megan lon limits = '',1p,2g25.16)') megan_lons(1),megan_lons(nlon_megan)
      write(*,'('' megan2_bioemiss: megan lat limits = '',1p,2g25.16)') megan_lats(1),megan_lats(nlat_megan)
      write(*,*) ' '
!---------------------------------------------------------------
!     allocate memory space to store interpolation coef.
!---------------------------------------------------------------
      ierr = 0
      allocate( ax(ids:ide,jds:jde,0:1), &
                by(ids:ide,jds:jde,0:1), stat=status )
      ierr = ierr + status
      allocate( ix(ids:ide,jds:jde,0:1), &
                jy(ids:ide,jds:jde,0:1), stat=status )
      ierr = ierr + status
      if( ierr /= 0 ) then
        write(*,*) 'allocate for ax ... jy failed; error = ',ierr
      endif
!---------------------------------------------------------------------
!   set bilinear interp variables
!---------------------------------------------------------------------
lat_loop : &
      do j = jds,jde
        do i = ids,ide
!---------------------------------------------------------------------
!   longitudes
!---------------------------------------------------------------------
          wrf_lon = xlong(i,j)
          if( wrf_lon >= megan_lons(nlon_megan) .or. &
              wrf_lon < megan_lons(1) ) then
            ix(i,j,lower) = nlon_megan
          else
            do n = 2,nlon_megan
              if( wrf_lon < megan_lons(n) ) then
                ix(i,j,lower) = min( nlon_megan-1, max(n-1,1) )
                exit
              endif
            end do
          endif
          ix(i,j,upper) = mod( ix(i,j,lower),nlon_megan ) + 1
          data_dx = megan_lons(ix(i,j,upper)) - megan_lons(ix(i,j,lower))
          if( data_dx < 0. ) then
            data_dx = 360. + data_dx
          endif
          ds1 = wrf_lon - megan_lons(ix(i,j,lower))
          if( ds1 < 0. ) then
            ds1 = 360. + ds1
          endif
          ax(i,j,lower) = ds1/data_dx
          ax(i,j,upper) = 1.0 - ax(i,j,lower)
!---------------------------------------------------------------------
!   latitudes
!---------------------------------------------------------------------
          wrf_lat = xlat(i,j)
          if( wrf_lat < megan_lats(1) ) then
            jy(i,j,0:1) = -1
            by(i,j,0:1) = 0.
          elseif( wrf_lat > megan_lats(nlat_megan) ) then
            jy(i,j,0:1) = -2
            by(i,j,0:1) = 0.
          else
            do n = 1,nlat_megan
              if( wrf_lat < megan_lats(n) ) then
                exit
              endif
            end do
            jy(i,j,lower) = min( nlat_megan-1, max(n-1,1) )
            jy(i,j,upper) = jy(i,j,lower) + 1
            by(i,j,lower) = (wrf_lat - megan_lats(jy(i,j,lower))) &
                             /(megan_lats(jy(i,j,upper)) - megan_lats(jy(i,j,lower)))
            by(i,j,upper) = 1.0 - by(i,j,lower)
          endif
        end do
      end do lat_loop
!---------------------------------------------------------------------
!   form dataset index limits
!---------------------------------------------------------------------
      write(*,*) 'megan2_bioemiss: count of points <,> data min,max lat = ',count(ix(:,:,0) == nlon_megan )
      xndx_megan(1) = minval( ix(:,:,lower) )
      xndx_megan(2) = maxval( ix(:,:,upper) )
      write(*,*) 'xndx_megan = ',xndx_megan(:)
      write(*,*) 'megan2_bioemiss: count of points < data min lat = ',count(jy(:,:,0) == -1)
      write(*,*) 'megan2_bioemiss: count of points > data max lat = ',count(jy(:,:,0) == -2)
      yndx_megan(1) = minval( jy(:,:,lower),mask=jy(:,:,lower)>0 )
      yndx_megan(2) = maxval( jy(:,:,upper),mask=jy(:,:,upper)>0 )
      write(*,*) 'yndx_megan = ',yndx_megan(:)

      if( debug ) then
      write(*,*) ' '
      write(*,*) 'megan2_bioemiss: bilinear interp diagnostics'
      write(*,*) 'megan2_bioemiss: ix'
      write(*,*) ix(ids,jds,:)
      write(*,*) 'megan2_bioemiss: ax'
      write(*,*) ax(ids,jds,:)
      write(*,*) 'megan2_bioemiss: megan lons'
      write(*,*) megan_lons(ix(ids,jds,0)),megan_lons(ix(ids,jds,1))
      write(*,*) 'megan2_bioemiss: wrf lon = ',xlong(ids,jds)
      write(*,*) 'megan2_bioemiss: jy'
      write(*,*) jy(ids,jds,:)
      write(*,*) 'megan2_bioemiss: by'
      write(*,*) by(ids,jds,:)
      write(*,*) 'megan2_bioemiss: megan lats'
      write(*,*) megan_lats(jy(ids,jds,0)),megan_lats(jy(ids,jds,1))
      write(*,*) 'megan2_bioemiss: wrf lat = ',xlat(ids,jds)
      write(*,*) ' '
      do j = jds,jde
        do i = ids,ide
          if( ix(i,j,lower) == nlon_megan ) then
      write(*,*) 'megan2_bioemiss: bilinear interp diagnostics'
      write(*,*) 'megan2_bioemiss: ix'
      write(*,*) ix(i,j,:)
      write(*,*) 'megan2_bioemiss: ax'
      write(*,*) ax(i,j,:)
      write(*,*) 'megan2_bioemiss: megan lons'
      write(*,*) megan_lons(ix(i,j,0)),megan_lons(ix(i,j,1))
      write(*,*) 'megan2_bioemiss: wrf lon = ',xlong(i,j)
      write(*,*) 'megan2_bioemiss: jy'
      write(*,*) jy(i,j,:)
      write(*,*) 'megan2_bioemiss: by'
      write(*,*) by(i,j,:)
      write(*,*) 'megan2_bioemiss: megan lats'
      write(*,*) megan_lats(jy(i,j,0)),megan_lats(jy(i,j,1))
      write(*,*) 'megan2_bioemiss: wrf lat = ',xlat(i,j)
            stop 'diagnostics'
          endif
        end do
      end do
      stop 'diagnostics'
      endif

!---------------------------------------------------------------------
!   allocate and read dataset variable
!---------------------------------------------------------------------
       if( allocated( tmp3 ) ) then
          deallocate( tmp3 )
       endif
       allocate( tmp3(xndx_megan(1):xndx_megan(2),yndx_megan(1):yndx_megan(2),mnth_s:mnth_e),stat=ierr )
                 
       if( ierr /= 0 ) then
         write(message,*) 'Failed to allocate tmp3 for lai; error = ',ierr
         stop 'bio_emiss abort'
       endif
       message = 'Failed to get variable id'
       call handle_ncerr( nf_inq_varid( ncid, trim(varname), varid ), message )
       message = 'Failed to read variable'
       if( mnth_s == mnth_e ) then
          call handle_ncerr( nf_get_vara_real( ncid, varid, &
                                               (/ xndx_megan(1),yndx_megan(1) /), &
                                               (/ xndx_megan(2)-xndx_megan(1)+1,yndx_megan(2)-yndx_megan(1)+1 /), &
                                               tmp3 ), message )
       else
          call handle_ncerr( nf_get_vara_real( ncid, varid, &
                                               (/ xndx_megan(1),yndx_megan(1),mnth_s /), &
                                               (/ xndx_megan(2)-xndx_megan(1)+1, &
                                                  yndx_megan(2)-yndx_megan(1)+1,mnth_e-mnth_s+1 /), &
                                               tmp3 ), message )
       endif

       write(*,*)  'dataset size = ',size(tmp3)
       write(*,*)  'dataset min,max values = ',minval(tmp3(:,:,:)),maxval(tmp3(:,:,:))
       write(*,*)  'dataset missing value count = ',count(tmp3(:,:,:) == missing_value )
       write(*,*)  '% valid data = ',100.* real( count(tmp3(:,:,:) /= missing_value ) ) /real(size(tmp3))

!---------------------------------------------------------------------
!   replace missing values with zero
!---------------------------------------------------------------------
       where( tmp3(:,:,:) == missing_value )
          tmp3(:,:,:) = 0.
       endwhere
!---------------------------------------------------------------------
!   set wrf bioemission variable
!---------------------------------------------------------------------
mnth_loop : &
       do mnth = mnth_s,mnth_e
         do j = jds,jde
           do i = ids,ide
             jl = jy(i,j,0)
             if( jl > 0 ) then
               il = ix(i,j,0)
               iu = ix(i,j,1)
               ju = jy(i,j,1)
               wrk_sum  = tmp3(il,jl,mnth)*ax(i,j,upper)*by(i,j,upper) &
                        + tmp3(il,ju,mnth)*ax(i,j,upper)*by(i,j,lower) &
                        + tmp3(iu,jl,mnth)*ax(i,j,lower)*by(i,j,upper) &
                        + tmp3(iu,ju,mnth)*ax(i,j,lower)*by(i,j,lower)
             else
               wrk_sum  = 0.
             endif
             if( varname(:3) == 'LAI' ) then
               mlai(i,j,megan_month) = scale_factor * wrk_sum
             elseif( trim(varname) == 'isoprene' ) then
               msebio_isop(i,j) = scale_factor * wrk_sum
             elseif( varname(:14) == 'Broadleaf_tree' ) then
               pftp_bt(i,j) = scale_factor * wrk_sum
             elseif( varname(:15) == 'Needleleaf_tree' ) then
               pftp_nt(i,j) = scale_factor * wrk_sum
             elseif( varname(:10) == 'Herbaceous' ) then
               pftp_hb(i,j) = scale_factor * wrk_sum
             elseif( varname(:5) == 'Shrub' ) then
               pftp_sb(i,j) = scale_factor * wrk_sum
             elseif( varname(:7) == 'TAS_AVE' ) then
               mtsa(i,j,mnth) = scale_factor * wrk_sum
             elseif( varname(:7) == 'DSW_AVE' ) then
               mswdown(i,j,mnth) = scale_factor * wrk_sum
             endif
           end do
         end do
       end do mnth_loop
    endif is_area_map

!---------------------------------------------------------------------
!   exit
!---------------------------------------------------------------------
    if( has_area_map ) then
      deallocate( xedge_megan, yedge_megan )
    else
       deallocate( ix, jy, ax, by, tmp3 )
    endif

    write(*,*) ' Finished megan2 bio emiss dataset ',trim(inpname)
!---------------------------------------------------------------------
!   close megan dataset file
!---------------------------------------------------------------------
    message = 'Failed to close ' // trim(inpname)
    call handle_ncerr( nf_close( ncid ), message )       

    end subroutine megan2_bioemiss

    subroutine write_bioemiss
!---------------------------------------------------------------------
!	... write the netcdf bio emission file
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!	... local variables
!---------------------------------------------------------------------
      integer :: lon_id
      integer :: lat_id
      integer :: time_id
      integer :: months_id
      integer :: string_id
      integer :: dims(4)
      integer :: start_ndx(4)
      integer :: length(4)
      character(len=132) :: message, text
      character(len=10)  :: ctime
      character(len=8)   :: cdate
      character(len=10)  :: t_string(2)

      write(num,'(i3)') 100+domain
      outpname = 'wrfbiochemi_d' // num(2:3)
!-----------------------------------------------------------------------
!     	... create netcdf bio emission file and enter define mode
!-----------------------------------------------------------------------
      message = 'write_bioemiss: Failed to create ' // trim( outpname )
      call handle_ncerr( nf_create( trim( outpname ), nf_clobber, ncid ), message )

!-----------------------------------------------------------------------
!     	... define the dimensions
!-----------------------------------------------------------------------
      call handle_ncerr( nf_def_dim( ncid, 'west_east', ide, lon_id ), &
                         'write_bioemiss: Failed to define longitude dimension' )
      call handle_ncerr( nf_def_dim( ncid, 'south_north', jde, lat_id ), &
                         'write_bioemiss: Failed to define latitude dimension' )
      call handle_ncerr( nf_def_dim( ncid, 'months_per_year_stag', 12, months_id ), &
                         'write_bioemiss: Failed to define months_per_year_stag dimension' )
      call handle_ncerr( nf_def_dim( ncid, 'DateStrLen', 19, string_id ), &
                         'write_bioemiss: Failed to define DateStrLen dimension' )
      call handle_ncerr( nf_def_dim( ncid, 'Time', nf_unlimited, time_id ), &
                         'write_bioemiss: Failed to create Time dimension' )

!-----------------------------------------------------------------------
!     	... define the variables
!-----------------------------------------------------------------------
      dims(1:2) = (/ string_id, time_id /)
      call handle_ncerr( nf_def_var( ncid, 'Times', nf_char, 2, dims(1:2), varid ), &
                         'write_bioemiss: Failed to define Times variable' )
      dims(1:3) = (/ lon_id, lat_id, time_id /)
      call handle_ncerr( nf_def_var( ncid, 'XLONG', nf_float, 3, dims(1:3), varid ), &
                         'write_bioemiss: Failed to define XLONG variable' )
      call handle_ncerr( nf_def_var( ncid, 'XLAT', nf_float, 3, dims(1:3), varid ), &
                         'write_bioemiss: Failed to define XLAT variable' )
      call handle_ncerr( nf_def_var( ncid, 'MSEBIO_ISOP', nf_float, 3, dims(1:3), varid ), &
                         'write_bioemiss: Failed to define MSEBIO_ISOP variable' )
      call handle_ncerr( nf_def_var( ncid, 'PFTP_BT', nf_float, 3, dims(1:3), varid ), &
                         'write_bioemiss: Failed to define PFTP_BT variable' )
      call handle_ncerr( nf_def_var( ncid, 'PFTP_NT', nf_float, 3, dims(1:3), varid ), &
                         'write_bioemiss: Failed to define PFTP_NT variable' )
      call handle_ncerr( nf_def_var( ncid, 'PFTP_SB', nf_float, 3, dims(1:3), varid ), &
                         'write_bioemiss: Failed to define PFTP_SB variable' )
      call handle_ncerr( nf_def_var( ncid, 'PFTP_HB', nf_float, 3, dims(1:3), varid ), &
                         'write_bioemiss: Failed to define PFTP_HB variable' )
      dims(:) = (/ lon_id, lat_id, months_id, time_id /)
      call handle_ncerr( nf_def_var( ncid, 'MLAI', nf_float, 4, dims, varid ), &
                         'write_bioemiss: Failed to define MLAI variable' )
      call handle_ncerr( nf_def_var( ncid, 'MTSA', nf_float, 4, dims, varid ), &
                         'write_bioemiss: Failed to define MTSA variable' )
      call handle_ncerr( nf_def_var( ncid, 'MSWDOWN', nf_float, 4, dims, varid ), &
                         'write_bioemiss: Failed to define MSWDOWN variable' )

!-----------------------------------------------------------------------
!   ... define variable attributes
!-----------------------------------------------------------------------
      varname = 'XLONG'
      units_attribute       = 'degree east'
      description_attribute = 'LONGITUDE, WEST IS NEGATIVE'
      stagger_attribute     = ''
      memord_attribute      = 'XY '
      coor_attribute        = ' '
      CALL write_bioemiss_attributes

      varname = 'XLAT'
      units_attribute       = 'degree north'
      description_attribute = 'LATITUDE, SOUTH IS NEGATIVE'
      stagger_attribute     = ''
      memord_attribute      = 'XY '
      coor_attribute        = ' '
      CALL write_bioemiss_attributes

      varname = 'MSEBIO_ISOP'
      units_attribute       = 'mol km^-2 hr^-1'
      description_attribute = 'isoprene emission factor'
      stagger_attribute     = ''
      memord_attribute      = 'XY '
      coor_attribute        = 'XLONG XLAT'
      CALL write_bioemiss_attributes

      varname = 'PFTP_BT'
      description_attribute = 'MEGAN2 Broadleaf PFT % coverage'
      units_attribute       = '%'
      CALL write_bioemiss_attributes

      varname = 'PFTP_NT'
      description_attribute = 'MEGAN2 Needleleaf PFT % coverage'
      CALL write_bioemiss_attributes

      varname = 'PFTP_SB'
      description_attribute = 'MEGAN2 Shrub and Bush PFT % coverage'
      CALL write_bioemiss_attributes

      varname = 'PFTP_HB'
      description_attribute = 'MEGAN2 Herbs PFT % coverage'
      CALL write_bioemiss_attributes

      varname = 'MLAI'
      units_attribute   = ''
      stagger_attribute = 'Z'
      memord_attribute  = 'XYZ'
      description_attribute = 'Monthly Leaf area index for MEGAN2'
      CALL write_bioemiss_attributes

      varname = 'MTSA'
      units_attribute   = 'K'
      description_attribute = 'Monthly surface air temperature'
      CALL write_bioemiss_attributes

      varname = 'MSWDOWN'
      units_attribute   = 'W/m2'
      description_attribute = 'Monthly SWdown'
      CALL write_bioemiss_attributes

!-----------------------------------------------------------------------
!   ... define global attributes
!-----------------------------------------------------------------------
      message = 'global_attributes: Failed to write title'
      text    = 'MEGAN2 bio emissions'
      call handle_ncerr( nf_put_att_text( ncid, nf_global, 'Title', len_trim(text), trim(text) ), message )
      message = 'global_attributes: Failed to write History'
      call date_and_time( cdate, ctime )
      t_string(1) = cdate(1:4) // '-' // cdate(5:6) // '-' // cdate(7:8)
      t_string(2) = ctime(1:2) // ':' // ctime(3:4)
      text    = 'Created on ' // trim(t_string(1)) // ' at ' // trim(t_string(2))
      call handle_ncerr( nf_put_att_text( ncid, nf_global, 'History', len_trim(text), trim(text) ), message )
      message = 'global_attributes: Author to write Files'
      text    = 'megan_bio_emiss'
      call handle_ncerr( nf_put_att_text( ncid, nf_global, 'Author', len_trim(text), trim(text) ), message )

      if( ngatts > 0 ) then
        call set_bioemiss_glb_atts
      endif

!-----------------------------------------------------------------------
!     	... leave define mode
!-----------------------------------------------------------------------
      call handle_ncerr( nf_enddef( ncid ), 'write_bioemiss: Failed to leave define mode' )

!-----------------------------------------------------------------------
!     	... write the variables
!-----------------------------------------------------------------------
      start_ndx(1:2) = (/ 1,1 /)
      length(1:2)    = (/ 19, 1 /)
      call handle_ncerr( nf_inq_varid( ncid, 'Times', varid ), &
                         'write_bioemiss: Failed to get Times variable id' )
      call handle_ncerr( nf_put_vara_text( ncid, varid, start_ndx(:2), length(:2), Times ), &
                         'write_bioemiss: Failed to write Times variable' )
      start_ndx(1:3) = (/ 1,1,1 /)
      length(1:3)    = (/ ide, jde, 1 /)
      call handle_ncerr( nf_inq_varid( ncid, 'XLONG', varid ), &
                         'write_bioemiss: Failed to get xlong variable id' )
      call handle_ncerr( nf_put_vara_real( ncid, varid, start_ndx(:3), length(:3), xlong ), &
                         'write_bioemiss: Failed to write xlong variable' )
      call handle_ncerr( nf_inq_varid( ncid, 'XLAT', varid ), &
                         'write_bioemiss: Failed to get xlat variable id' )
      call handle_ncerr( nf_put_vara_real( ncid, varid, start_ndx(:3), length(:3), xlat ), &
                         'write_bioemiss: Failed to write xlat variable' )
      call handle_ncerr( nf_inq_varid( ncid, 'MSEBIO_ISOP', varid ), &
                         'write_bioemiss: Failed to get msebio_isop variable id' )
      call handle_ncerr( nf_put_vara_real( ncid, varid, start_ndx(:3), length(:3), msebio_isop ), &
                         'write_bioemiss: Failed to write msebio_isop variable' )
      call handle_ncerr( nf_inq_varid( ncid, 'PFTP_BT', varid ), &
                         'write_bioemiss: Failed to get pftp_bt variable id' )
      call handle_ncerr( nf_put_vara_real( ncid, varid, start_ndx(:3), length(:3), pftp_bt ), &
                         'write_bioemiss: Failed to write pftp_bt variable' )
      call handle_ncerr( nf_inq_varid( ncid, 'PFTP_NT', varid ), &
                         'write_bioemiss: Failed to get pftp_nt variable id' )
      call handle_ncerr( nf_put_vara_real( ncid, varid, start_ndx(:3), length(:3), pftp_nt ), &
                         'write_bioemiss: Failed to write pftp_nt variable' )
      call handle_ncerr( nf_inq_varid( ncid, 'PFTP_HB', varid ), &
                         'write_bioemiss: Failed to get pftp_hb variable id' )
      call handle_ncerr( nf_put_vara_real( ncid, varid, start_ndx(:3), length(:3), pftp_hb ), &
                         'write_bioemiss: Failed to write pftp_hb variable' )
      call handle_ncerr( nf_inq_varid( ncid, 'PFTP_SB', varid ), &
                         'write_bioemiss: Failed to get pftp_sb variable id' )
      call handle_ncerr( nf_put_vara_real( ncid, varid, start_ndx(:3), length(:3), pftp_sb ), &
                         'write_bioemiss: Failed to write pftp_sb variable' )

      start_ndx(:) = (/ 1,1,1,1 /)
      length(:)    = (/ ide, jde, 1, 1 /)
      do mnth = 1,12
         start_ndx(3) = mnth
         call handle_ncerr( nf_inq_varid( ncid, 'MLAI', varid ), &
                            'write_bioemiss: Failed to get mlai variable id' )
         call handle_ncerr( nf_put_vara_real( ncid, varid, start_ndx, length, mlai(:,:,mnth) ), &
                            'write_bioemiss: Failed to write mlai variable' )
         call handle_ncerr( nf_inq_varid( ncid, 'MTSA', varid ), &
                            'write_bioemiss: Failed to get mtsa variable id' )
         call handle_ncerr( nf_put_vara_real( ncid, varid, start_ndx, length, mtsa(:,:,mnth) ), &
                            'write_bioemiss: Failed to write mtsa variable' )
         call handle_ncerr( nf_inq_varid( ncid, 'MSWDOWN', varid ), &
                            'write_bioemiss: Failed to get mswdown variable id' )
         call handle_ncerr( nf_put_vara_real( ncid, varid, start_ndx, length, mswdown(:,:,mnth) ), &
                            'write_bioemiss: Failed to write mswdown variable' )
      end do
!---------------------------------------------------------------------
!   close wrf file
!---------------------------------------------------------------------
       message = 'Failed to close ' // trim(outpname)
       call handle_ncerr( nf_close( ncid ), message )       

   end subroutine write_bioemiss

   subroutine write_bioemiss_attributes
!---------------------------------------------------------------------
!   write common variable attributes
!---------------------------------------------------------------------

      message = 'write_bioemiss: Failed to get ' // trim(varname) // ' variable id'
      call handle_ncerr( nf_inq_varid( ncid, trim(varname), varid ), message )
      message = 'write_bioemiss: Failed to create ' // trim(varname) // ' attribute'
      call handle_ncerr( nf_put_att_text( ncid, varid, 'MemoryOrder', 3, memord_attribute ), message )
      call handle_ncerr( nf_put_att_text( ncid, varid, 'description', &
                                          len_trim(description_attribute), trim(description_attribute) ), message )
      call handle_ncerr( nf_put_att_text( ncid, varid, 'units', &
                                          len_trim(units_attribute), trim(units_attribute) ), message )
      call handle_ncerr( nf_put_att_text( ncid, varid, 'stagger', &
                                          len_trim(stagger_attribute), trim(stagger_attribute) ), message )
      if( coor_attribute /= ' ' ) then
         call handle_ncerr( nf_put_att_text( ncid, varid, 'coordinates', &
                                             len_trim(coor_attribute), trim(coor_attribute) ), message )
      endif
      ii = 104
      call handle_ncerr( nf_put_att_int( ncid, varid, 'FieldType', nf_int, 1, ii ), message )

   end subroutine write_bioemiss_attributes

   subroutine get_bioemiss_glb_atts
!---------------------------------------------------------------------
!   read the global attributes
!---------------------------------------------------------------------

   implicit none

!---------------------------------------------------------------------
!   local variables
!---------------------------------------------------------------------
   integer              :: m
   integer              :: astat
   integer              :: attr_len
   character(len=132)   :: message
   character(len=132)   :: attr_name

!---------------------------------------------------------------------
!   get global attr count
!---------------------------------------------------------------------
   message = 'glb_attr: Failed to get glb attr count'
   call handle_ncerr( nf_inq_natts( ncid, ngatts ), message )       
!---------------------------------------------------------------------
!   allocate variables
!---------------------------------------------------------------------
   if( ngatts > 0 ) then
     allocate( attrs(ngatts),stat=astat )
     if( astat /= 0 ) then
       write(*,*) 'glb_attr: failed to allocate type glb_att'
       stop 'Alloc err'
     endif
     attrs(:)%name = ' '
   endif
!---------------------------------------------------------------------
!   loop over glb attributes
!---------------------------------------------------------------------
glb_attr_loop : &
   do m = 1,ngatts
     write(message,*) 'glb_attr: Failed to get glb attr # ',m,' name'
     call handle_ncerr( nf_inq_attname( ncid, nf_global, m, attr_name ), message )       
     attrs(m)%name = attr_name
     write(message,*) 'glb_attr: Failed to get glb attr # ',m,' type,len'
     call handle_ncerr( nf_inq_att( ncid, nf_global, trim(attr_name), attrs(m)%type, attr_len ), message )       
     attrs(m)%len = attr_len
     message = 'glb_attr: Failed to get ' // trim(attr_name)
     select case( attrs(m)%type )
       case( nf_byte )
         allocate( attrs(m)%attr_byte(attr_len),stat=astat )
         if( astat /= 0 ) then
           write(*,*) 'glb_attr: failed to allocate attr_byte'
           stop 'Alloc err'
         endif
         call handle_ncerr( nf_get_att_int1( ncid, nf_global, trim(attr_name), attrs(m)%attr_byte ), message )       
       case( nf_char )
         attrs(m)%attr_char = ' '
         call handle_ncerr( nf_get_att_text( ncid, nf_global, trim(attr_name), attrs(m)%attr_char ), message )       
       case( nf_short )
         allocate( attrs(m)%attr_short(attr_len),stat=astat )
         if( astat /= 0 ) then
           write(*,*) 'glb_attr: failed to allocate attr_short'
           stop 'Alloc err'
         endif
         call handle_ncerr( nf_get_att_int2( ncid, nf_global, trim(attr_name), attrs(m)%attr_short ), message )       
       case( nf_int )
         allocate( attrs(m)%attr_int(attr_len),stat=astat )
         if( astat /= 0 ) then
           write(*,*) 'glb_attr: failed to allocate attr_int'
           stop 'Alloc err'
         endif
         call handle_ncerr( nf_get_att_int( ncid, nf_global, trim(attr_name), attrs(m)%attr_int ), message )       
       case( nf_float )
         allocate( attrs(m)%attr_real(attr_len),stat=astat )
         if( astat /= 0 ) then
           write(*,*) 'glb_attr: failed to allocate attr_real'
           stop 'Alloc err'
         endif
         call handle_ncerr( nf_get_att_real( ncid, nf_global, trim(attr_name), attrs(m)%attr_real ), message )       
       case( nf_double )
         allocate( attrs(m)%attr_dbl(attr_len),stat=astat )
         if( astat /= 0 ) then
           write(*,*) 'glb_attr: failed to allocate attr_dbl'
           stop 'Alloc err'
         endif
         call handle_ncerr( nf_get_att_double( ncid, nf_global, trim(attr_name), attrs(m)%attr_dbl ), message )       
     end select
   end do glb_attr_loop

   end subroutine get_bioemiss_glb_atts

   subroutine set_bioemiss_glb_atts
!---------------------------------------------------------------------
!   set the global attributes
!---------------------------------------------------------------------

   implicit none

!---------------------------------------------------------------------
!   local variables
!---------------------------------------------------------------------
   integer              :: m
   integer              :: attr_len
   integer              :: attr_xtype
   integer              :: slen
   character(len=132)   :: message
   character(len=132)   :: attr_name

!---------------------------------------------------------------------
!   loop over glb attributes
!---------------------------------------------------------------------
glb_attr_loop : &
   do m = 1,ngatts
     attr_name = trim(attrs(m)%name)
     slen      = len_trim(attr_name)
     write(message,*) 'set_glb_att: Failed to define glb att ',trim(attr_name)
     attr_len   = attrs(m)%len
     attr_xtype = attrs(m)%type
     select case( attrs(m)%type )
       case( nf_byte )
         call handle_ncerr( nf_put_att_int1( ncid, nf_global, attr_name(:slen), attr_xtype, attr_len, attrs(m)%attr_byte ), message )       
       case( nf_char )
         call handle_ncerr( nf_put_att_text( ncid, nf_global, attr_name(:slen), attr_len, attrs(m)%attr_char ), message )       
       case( nf_short )
         call handle_ncerr( nf_put_att_int2( ncid, nf_global, attr_name(:slen), attr_xtype, attr_len, attrs(m)%attr_short ), message )       
       case( nf_int )
         call handle_ncerr( nf_put_att_int( ncid, nf_global, attr_name(:slen), attr_xtype, attr_len, attrs(m)%attr_int ), message )       
       case( nf_float )
         call handle_ncerr( nf_put_att_real( ncid, nf_global, attr_name(:slen), attr_xtype, attr_len, attrs(m)%attr_real ), message )       
       case( nf_double )
         call handle_ncerr( nf_put_att_double( ncid, nf_global, attr_name(:slen), attr_xtype, attr_len, attrs(m)%attr_dbl ), message )       
     end select
   end do glb_attr_loop

   end subroutine set_bioemiss_glb_atts

   subroutine dealloc_bioemiss_glb_atts
!---------------------------------------------------------------------
!   deallocate variables
!---------------------------------------------------------------------

   integer :: m

   do m = 1,ngatts
     select case( attrs(m)%type )
       case( nf_byte )
         deallocate( attrs(m)%attr_byte )
       case( nf_short )
         deallocate( attrs(m)%attr_short )
       case( nf_int )
         deallocate( attrs(m)%attr_int )
       case( nf_float )
         deallocate( attrs(m)%attr_real )
       case( nf_double )
         deallocate( attrs(m)%attr_dbl )
     end select
   end do

   end subroutine dealloc_bioemiss_glb_atts

   subroutine handle_ncerr( ret, mes )
!---------------------------------------------------------------------
!	... netcdf error handling routine
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!	... dummy arguments
!---------------------------------------------------------------------
   integer, intent(in) :: ret
   character(len=*), intent(in) :: mes

   if( ret /= nf_noerr ) then
      write(*,*) nf_strerror( ret )
      stop 'netcdf error'
   endif

   end subroutine handle_ncerr

end program map_megan2_emissions
