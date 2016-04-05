
      module module_mozart_lib

      implicit none

!---------------------------------------------------------------
!     include files
!---------------------------------------------------------------
      include 'netcdf.inc'

!---------------------------------------------------------------
!     public procedures
!---------------------------------------------------------------
      public  :: init_mozart_lib
      public  :: bc_interpolate4d
      public  :: ic_interpolate4d
      public  :: exit_mozart_lib
      public  :: get_moz_time_ndx

!---------------------------------------------------------------
!     private procedures
!---------------------------------------------------------------
      private :: handle_error
      private :: lotim

!---------------------------------------------------------------
!     public variables
!---------------------------------------------------------------
      integer, public  :: nlon, nlat, nlev                    ! data dimension length
      integer, public  :: year, month, day
      integer, public  :: hour, minute, second
      integer, public  :: ntime_m
      real, public, allocatable :: moz_times(:)
      character(len=9) :: moz_var_suffix = '_VMR_inst'

!---------------------------------------------------------------
!     variables for reading netCDF data
!---------------------------------------------------------------
      integer, private :: ncid = 0                           ! netCDF file ID
      integer, private :: nstt(4), ncnt(4)                   ! start index and counter array
      integer, private :: start_index = 1                    ! monthly mean data
      integer, private :: ntime                              ! data dimension length
      integer, private :: nscalar, nstring                   ! data dimension length
      integer, private, allocatable :: moz_date(:)
      integer, private, allocatable :: moz_datesec(:) 

!---------------------------------------------------------------
!     variables used by interpolation
!---------------------------------------------------------------
      real, private, parameter      :: ps0 = 1.e5

      integer, private, allocatable :: ix(:,:,:)
      integer, private, allocatable :: jy(:,:,:)                        ! index used by interpolation
      real, private, allocatable    :: ax(:,:,:)                        ! weight coef. all domain
      real, private, allocatable    :: by(:,:,:)                        ! weight coef. all domain
      real, private, allocatable    :: ps_moz(:,:)
      real, private, allocatable    :: ps_mozi(:,:)
      real, private, allocatable    :: hyam(:)
      real, private, allocatable    :: hybm(:)
      character(len=16) :: dbg_species = ' '

      type t_type
        integer :: ndx
        integer :: ncid_lo, ncid_hi
        integer :: lo_moz_ndx, hi_moz_ndx
        integer :: lo_buf_ndx, hi_buf_ndx
        integer :: gap_date, gap_secs
        real    :: dels
        real, pointer :: conc(:,:,:,:)
        real, pointer :: ps(:,:,:)
        logical :: in_gap
        logical :: t_interp
      end type t_type

      type(t_type), private :: time_type

      contains

      subroutine init_mozart_lib( moz_dir, moz_fn, x2d, y2d, wrf_date, &
                                  wrf_datesec, nx, ny, nspec )
!---------------------------------------------------------------
!     initialize netcdf file
!---------------------------------------------------------------

      use mo_calendar, only : diffdat
      use utils,       only : wrf2mz_map

!---------------------------------------------------------------
!     dummy arguments
!---------------------------------------------------------------
      integer, intent(in)      :: wrf_date
      integer, intent(in)      :: wrf_datesec
      integer, intent(in)      :: nx, ny
      integer, intent(in)      :: nspec
      character(*), intent(in) :: moz_dir
      character(*), intent(inout) :: moz_fn
      real, intent(in)         :: x2d(nx,ny)
      real, intent(in)         :: y2d(nx,ny)

!---------------------------------------------------------------
!     local variables
!---------------------------------------------------------------
      integer :: status, tstat
      integer :: dimid, varid
      integer :: i, j, k, n 
      integer :: vid
      real    :: model_x, model_y, model_dx
      real, allocatable  :: x1d(:), y1d(:)
      character(len=128) :: filenm
      character(len=32)  :: spcnam
      logical            :: found
      logical            :: monotonic_moz_x

      filenm = trim( moz_dir ) // adjustl( moz_fn )
!---------------------------------------------------------------
!     open the initial mozart netCDF file
!---------------------------------------------------------------
      status = nf_open( trim(filenm), nf_nowrite, ncid )
      if( status /= nf_noerr ) then
         write(*,*) 'failed to open ',trim(filenm)
         call handle_error( status )
      end if
      write(*,*) 'init_mozart_lib: opened ',trim(filenm)

!---------------------------------------------------------------
!     read times
!---------------------------------------------------------------
      call read_moz_times( ncid )

      time_type%ndx = lotim( wrf_date, wrf_datesec, moz_date, moz_datesec, ntime_m )
      if( time_type%ndx == 0 ) then
        status = nf_close( ncid )
        call next_flnm( moz_fn, .false. )
        filenm = trim( moz_dir ) // adjustl( moz_fn )
!---------------------------------------------------------------
!     open the initial mozart netCDF file
!---------------------------------------------------------------
        status = nf_open( trim(filenm), nf_nowrite, ncid )
        if( status /= nf_noerr ) then
          write(*,*) 'failed to open ',trim(filenm)
          call handle_error( status )
        end if
        write(*,*) 'init_mozart_lib: opened ',trim(filenm)
        call read_moz_times( ncid )
        time_type%ndx = lotim( wrf_date, wrf_datesec, moz_date, moz_datesec, ntime_m )
        if( time_type%ndx == 0 ) then
          write(*,*) 'init_mozart_lib: time ',wrf_date,' ',wrf_datesec
          write(*,*) '                 is before '
          write(*,*) '                 ',moz_date(1),' ',moz_datesec(1)
          write(*,*) '                 the first time in mozart file ',trim(moz_fn)
          stop
        end if
      end if
!---------------------------------------------------------------
!     check mozart variables
!---------------------------------------------------------------
      call chk_moz_vars( nspec, filenm )

!---------------------------------------------------------------
!     get the longitude dimension length of the MOZART data
!---------------------------------------------------------------
      status = nf_inq_dimid( ncid, 'lon', dimid )
      if( status /= nf_noerr )  call handle_error( status )

      status = nf_inq_dimlen( ncid, dimid, nlon )
      if( status /= nf_noerr )  call handle_error( status )

!---------------------------------------------------------------
!     get the latitude dimension length of the MOZART data
!---------------------------------------------------------------
      status = nf_inq_dimid( ncid, 'lat', dimid )
      if( status /= nf_noerr )  call handle_error( status )

      status = nf_inq_dimlen( ncid, dimid, nlat )
      if( status /= nf_noerr )  call handle_error( status )

!---------------------------------------------------------------
!     get the vertical dimension length of the MOZART data
!---------------------------------------------------------------
      status = nf_inq_dimid( ncid, 'lev', dimid )
      if( status /= nf_noerr )  call handle_error( status )

      status = nf_inq_dimlen( ncid, dimid, nlev )
      if( status /= nf_noerr )  call handle_error( status )

!---------------------------------------------------------------
!     read vertical coordinates
!---------------------------------------------------------------
      allocate( hyam(nlev), hybm(nlev), stat=status )
      if( status /= 0 ) then
         write(*,*) 'failed to allocate hyam, hybm; error = ',status
         stop
      end if

      status = nf_inq_varid( ncid, 'hyam', varid )
      if( status /= nf_noerr )  call handle_error( status )

      status = nf_get_var_real( ncid, varid, hyam )
      if( status /= nf_noerr )  call handle_error( status )

      status = nf_inq_varid( ncid, 'hybm', varid )
      if( status /= nf_noerr )  call handle_error( status )

      status = nf_get_var_real( ncid, varid, hybm )
      if( status /= nf_noerr )  call handle_error( status )

!---------------------------------------------------------------
!     read longitudes
!---------------------------------------------------------------
      allocate( x1d(nlon), y1d(nlat), stat=status )
      if( status /= 0 ) then
         write(*,*) 'failed to allocate x1d, y1d; error = ',status
         stop
      end if

      status = nf_inq_varid( ncid, 'lon', varid )
      if( status /= nf_noerr )  call handle_error( status )

      status = nf_get_var_real( ncid, varid, x1d )
      if( status /= nf_noerr )  call handle_error( status )

!---------------------------------------------------------------
!     check mozart longitudes for monotonicity
!---------------------------------------------------------------
      monotonic_moz_x = x1d(nlon) > x1d(1)
      if( .not. monotonic_moz_x ) then
        do i = 2,nlon
          if( x1d(i) < x1d(i-1) ) then
            x1d(i:nlon) = x1d(i:nlon) + 360.
            exit
          endif
        end do
      endif

!---------------------------------------------------------------
!     read latitudes
!---------------------------------------------------------------
      status = nf_inq_varid( ncid, 'lat', varid )
      if( status /= nf_noerr )  call handle_error( status )

      status = nf_get_var_real( ncid, varid, y1d )
      if( status /= nf_noerr )  call handle_error( status )

!---------------------------------------------------------------
!     allocate memory space to store interpolation coef.
!---------------------------------------------------------------
      tstat = 0
      allocate( ax(nx,ny,0:1), by(nx,ny,0:1), stat=status )
      tstat = tstat + status
      allocate( ix(nx,ny,0:1), jy(nx,ny,0:1), stat=status )
      tstat = tstat + status
      if( tstat /= 0 ) then
         write(*,*) 'allocate for ax ... jy failed; error = ',tstat
      end if

!---------------------------------------------------------------
!     horizontal interpolation coefs.
!---------------------------------------------------------------
!---------------------------------------------------------------
! all domain
!---------------------------------------------------------------
      do j = 1, ny 
         do i = 1, nx
            model_x = 360. + x2d(i,j)
            if( monotonic_moz_x ) then
              model_x = mod( model_x,360. )
            endif
            if( model_x >= x1d(nlon) ) then
               ix(i,j,0) = nlon
            else if( model_x < x1d(1) ) then
               ix(i,j,0) = nlon
            else
               do n = 1, nlon
                  if( model_x < x1d(n) ) then
                     ix(i,j,0) = min( nlon-1, max(n-1,1) )
                     exit
                  end if
               end do
            end if
            ix(i,j,1) = mod( ix(i,j,0),nlon ) + 1
            model_dx = x1d(ix(i,j,1)) - x1d(ix(i,j,0))
            if( model_dx < 0. ) then
               model_dx = 360. + model_dx
            end if
            ax(i,j,0) = (model_x - x1d(ix(i,j,0)))/model_dx
            ax(i,j,1) = 1.0 - ax(i,j,0)
            model_y = y2d(i,j)
            do n = 1, nlat
               if( model_y < y1d(n) ) then
                 exit
               end if
            end do
            jy(i,j,0) = min( nlat-1, max(n-1,1) )
            jy(i,j,1) = jy(i,j,0) + 1
            by(i,j,0) = (model_y - y1d(jy(i,j,0)))/(y1d(jy(i,j,1)) - y1d(jy(i,j,0)))
            by(i,j,1) = 1.0 - by(i,j,0)
         end do
      end do

!---------------------------------------------------------------
!     release memory
!---------------------------------------------------------------
      deallocate( x1d, y1d )

!---------------------------------------------------------------
!     allocate memory space for reading
!---------------------------------------------------------------
      allocate( time_type%conc(nlon,nlat,nlev,2), stat=status )
      if( status /= 0 ) then
         write(*,*) 'failed to allocate  time_type%conc; error = ',status
         stop
      end if
!---------------------------------------------------------------
!     allocate mozart surface pressure arrays
!---------------------------------------------------------------
      allocate( ps_moz(nlon,nlat), ps_mozi(nx,ny), &
                time_type%ps(nlon,nlat,2), stat=status )
      if( status /= 0 ) then
         write(*,*) 'failed to allocate ps_moz,ps_mozi; error = ',status
         stop
      end if

!---------------------------------------------------------------
!     setup time_type
!---------------------------------------------------------------
      time_type%lo_moz_ndx = 0
      time_type%hi_moz_ndx = 0
      time_type%lo_buf_ndx = 1
      time_type%hi_buf_ndx = 2
      time_type%ncid_lo    = 0
      time_type%ncid_hi    = 0
      time_type%in_gap     = .false.
      time_type%t_interp   = .false.

      if( time_type%ndx > 0 ) then
        time_type%lo_moz_ndx = time_type%ndx
        time_type%ncid_lo    = ncid
        time_type%ncid_hi    = ncid
      endif

      write(*,*) 'finished init_mozart_lib'

      end subroutine init_mozart_lib

      subroutine read_moz_times( ncid )
!---------------------------------------------------------------
!     read times from current mozart input file
!---------------------------------------------------------------

!---------------------------------------------------------------
!     dummy arguments
!---------------------------------------------------------------
      integer, intent(in) :: ncid

!---------------------------------------------------------------
!     local variables
!---------------------------------------------------------------
      integer :: varid
      integer :: status

      status = nf_inq_dimid( ncid, 'time', varid )
      if( status /= nf_noerr )  call handle_error( status )

      status = nf_inq_dimlen( ncid, varid, ntime_m )
      if( status /= nf_noerr ) call handle_error( status )

      if( allocated( moz_date ) ) then
        deallocate( moz_date, moz_datesec )
      end if
      allocate( moz_date(ntime_m), moz_datesec(ntime_m), stat=status )
      if( status /= 0 ) then
        write(*,*) 'failed to allocate date, datesec; error = ',status
        stop
      end if

      status = nf_inq_varid( ncid, 'date', varid )
      if( status /= nf_noerr )  call handle_error( status )

      status = nf_get_var_int( ncid, varid, moz_date )
      if( status /= nf_noerr )  call handle_error( status )

      status = nf_inq_varid( ncid, 'datesec', varid )
      if( status /= nf_noerr )  call handle_error( status )

      status = nf_get_var_int( ncid, varid, moz_datesec )
      if( status /= nf_noerr )  call handle_error( status )

      end subroutine read_moz_times

      subroutine get_moz_time_ndx( moz_dir, moz_fn, wrf_date, wrf_datesec, nspec )
!---------------------------------------------------------------
!     get mozart time index for wrf_date, wrf_datesec
!---------------------------------------------------------------

      use mo_calendar, only : diffdat

!---------------------------------------------------------------
!     dummy arguments
!---------------------------------------------------------------
      integer, intent(in)  :: wrf_date
      integer, intent(in)  :: wrf_datesec
      integer, intent(in)  :: nspec
      character(len=*), intent(in)    :: moz_dir
      character(len=*), intent(inout) :: moz_fn

!---------------------------------------------------------------
!     local variables
!---------------------------------------------------------------
      integer :: i, n
      integer :: status
      character(len=128) :: filenm
      logical :: found

      write(*,*) ' '
      write(*,*) 'get_moz_time_ndx; moz_dir,moz_fn = ',trim(moz_dir),trim(moz_fn)
      write(*,*) 'get_moz_time_ndx; wrf_date,wrf_datesec,ntime_m = ',wrf_date,wrf_datesec,ntime_m

      n = lotim( wrf_date, wrf_datesec, moz_date, moz_datesec, ntime_m )
      if( n > 0 ) then
!---------------------------------------------------------------
!     wrf time in present mozart dataset
!---------------------------------------------------------------
        write(*,*) 'get_moz_time_ndx; moz_tndx = ',n
        time_type%lo_moz_ndx = n
        time_type%in_gap     = .false.
        if( time_type%t_interp ) then
          time_type%hi_moz_ndx = n + 1
        endif
        if( time_type%ncid_hi /= time_type%ncid_lo ) then
          status = nf_close( time_type%ncid_lo )
          if( status /= 0 ) then
            filenm = trim( moz_dir ) // adjustl( moz_fn )
            write(*,*) 'get_moz_time_ndx: failed to close ',trim(filenm),' ; error = ',status
            stop
          end if
          time_type%ncid_lo = time_type%ncid_hi
        endif
        time_type%ncid_hi    = time_type%ncid_lo
      else if( n < 0 ) then
!---------------------------------------------------------------
!     wrf time after present mozart dataset
!---------------------------------------------------------------
         time_type%ncid_lo    = ncid
         time_type%lo_moz_ndx = ntime_m
         call next_flnm( moz_fn, .true. )
         filenm = trim( moz_dir ) // adjustl( moz_fn )
!---------------------------------------------------------------
!     open the input netCDF file
!---------------------------------------------------------------
         status = nf_open( trim(filenm), nf_nowrite, ncid )
         if( status /= nf_noerr ) call handle_error( status )
         write(*,*) 'get_moz_time_ndx: opened ',trim(filenm)
!---------------------------------------------------------------
!     check mozart variables
!---------------------------------------------------------------
         call chk_moz_vars( nspec, filenm )
         time_type%gap_date = moz_date(ntime_m)
         time_type%gap_secs = moz_datesec(ntime_m)
         call read_moz_times( ncid )
         n = lotim( wrf_date, wrf_datesec, moz_date, moz_datesec, ntime_m )
         time_type%ndx = n
         if( n > 0 ) then
           write(*,*) 'get_moz_time_ndx; moz_tndx = ',n
           status = nf_close( time_type%ncid_lo )
           if( status /= 0 ) then
             filenm = trim( moz_dir ) // adjustl( moz_fn )
             write(*,*) 'get_moz_time_ndx: failed to close ',trim(filenm),' ; error = ',status
             stop
           end if
           time_type%in_gap     = .false.
           time_type%ncid_lo    = ncid
           time_type%ncid_hi    = ncid
           time_type%lo_moz_ndx = n
         else if( n == 0 ) then
           time_type%in_gap     = .true.
           time_type%hi_moz_ndx = 1
           time_type%ncid_hi = ncid
           time_type%dels = time_type%dels/diffdat( time_type%gap_date, time_type%gap_secs, moz_date(1), moz_datesec(1) )
           time_type%t_interp = .true.
         else
           write(*,*) 'get_moz_time_ndx: failed to find ',wrf_date,' : ',wrf_datesec
           write(*,*) '                  in file ',trim(filenm)
           stop
         end if
      else if( time_type%in_gap ) then
        time_type%dels = time_type%dels/diffdat( time_type%gap_date, time_type%gap_secs, moz_date(1), moz_datesec(1) )
      end if

      end subroutine get_moz_time_ndx

      subroutine next_flnm( filenm, incr )
!---------------------------------------------------------------
!     increment mozart filename
!---------------------------------------------------------------

!---------------------------------------------------------------
!     dummy arguments
!---------------------------------------------------------------
      character(len=*), intent(inout) :: filenm
      logical, intent(in)             :: incr

!---------------------------------------------------------------
!     local variables
!---------------------------------------------------------------
      integer :: il, iu, ios
      integer :: file_number
      character(len=5) :: numa
      logical :: found

      found = .false.
      iu    = scan( trim(filenm), '0123456789', back=.true. )
      do il = iu-1,1,-1
         if( index( '0123456789', filenm(il:il) ) == 0 ) then
            found = .true.
            exit
         end if
      end do
      if( .not. found ) then
         write(*,*) 'next_filenm: mozart file ',trim(filenm),' is not proper format'
         stop
      end if

      il = il + 1
      write(*,*) ' '
      write(*,*) 'next_flnm; trying to increment file ',trim(filenm)
      write(*,*) 'next_flnm; il, iu = ',il,iu
      read(filenm(il:iu),*,iostat=ios) file_number
      if( ios /= 0 ) then
         write(*,*) 'next_filenm: failed to read ',filenm(il:iu),' ; error = ',ios
         stop
      end if
      write(*,*) 'next_flnm; file_number = ',file_number
      if( incr ) then
        file_number = file_number + 1
      else
        file_number = file_number - 1
      endif
      write(numa,'(i5)') file_number+10000
      filenm(il:iu) = numa(2:5)

      write(*,*) 'next_flnm; new file = ',trim(filenm)

      end subroutine next_flnm

      subroutine read_mozart_ps( bndx, tndx, ncid )
!---------------------------------------------------------------
!     read mozart surface pressure
!---------------------------------------------------------------

!---------------------------------------------------------------
!     dummy arguments
!---------------------------------------------------------------
      integer, intent(in) :: bndx              ! buffer index
      integer, intent(in) :: tndx              ! time index
      integer, intent(in) :: ncid              ! netcdf file index

!---------------------------------------------------------------
!     local variables
!---------------------------------------------------------------
      integer :: status
      integer :: varid

      status = nf_inq_varid( ncid, 'PS', varid )
      if( status /= nf_noerr )  then
        write(*,*) 'mozart_ps: failed to get PS id'
        call handle_error( status )
      end if

      nstt(1:3) = (/ 1, 1, tndx /)
      ncnt(1:3) = (/ nlon, nlat, 1 /)
      status = nf_get_vara_real( ncid, varid, nstt(1:3), ncnt(1:3), time_type%ps(:,:,bndx) )
      if( status /= nf_noerr ) then
        write(*,*) 'mozart_ps: failed to read PS'
        call handle_error( status )
      end if

      end subroutine read_mozart_ps

      subroutine tinterp_mozart_ps( dels )
!---------------------------------------------------------------
!     time interpolation of mozart surface pressure
!---------------------------------------------------------------

!---------------------------------------------------------------
!     dummy arguments
!---------------------------------------------------------------
      real, intent(in) :: dels                ! linear interp factor

!---------------------------------------------------------------
!     local variables
!---------------------------------------------------------------
      integer :: j
      integer :: lo_ndx, hi_ndx
      real    :: delsm1

      if( time_type%t_interp ) then
        delsm1 = 1. - dels
!---------------------------------------------------------------
!     interpolate mozart surface pressure to wrf grid
!---------------------------------------------------------------
        do j = 1,nlat 
          ps_moz(:,j) = time_type%ps(:,j,1) * dels &
                      + time_type%ps(:,j,2) * delsm1
        end do
      else
        do j = 1,nlat 
          ps_moz(:,j) = time_type%ps(:,j,1)
        end do
      endif

      end subroutine tinterp_mozart_ps

      subroutine hinterp_mozart_ps( nx, ny )
!---------------------------------------------------------------
!     horizontal interpolation of mozart surface pressure
!---------------------------------------------------------------

!---------------------------------------------------------------
!     dummy arguments
!---------------------------------------------------------------
      integer, intent(in) :: nx                ! wrf x dimension
      integer, intent(in) :: ny                ! wrf y dimension

!---------------------------------------------------------------
!     local variables
!---------------------------------------------------------------
      integer :: i, j

!---------------------------------------------------------------
!     interpolate mozart surface pressure to wrf grid
!---------------------------------------------------------------
      do j = 1,ny 
        do i = 1,nx 
          ps_mozi(i,j) = ps_moz(ix(i,j,0),jy(i,j,0))*ax(i,j,1)*by(i,j,1) &
                       + ps_moz(ix(i,j,0),jy(i,j,1))*ax(i,j,1)*by(i,j,0) &
                       + ps_moz(ix(i,j,1),jy(i,j,0))*ax(i,j,0)*by(i,j,1) &
                       + ps_moz(ix(i,j,1),jy(i,j,1))*ax(i,j,0)*by(i,j,0)
        end do
      end do

      end subroutine hinterp_mozart_ps

      subroutine bc_interpolate4d( ndx, wrfxs, wrfxe, wrfys, wrfye, &
                                   ps_wrf, znu, p_top, nx, ny, &
                                   nz, nw )
!---------------------------------------------------------------
!     interpolate four-dimensional field ... 
!---------------------------------------------------------------

      use utils, only : wrf2mz_map

!---------------------------------------------------------------
!     input arguments
!---------------------------------------------------------------
      integer, intent(in)      :: ndx
      integer, intent(in)      :: nx, ny, nz, nw          ! dimensions
      real, intent(in)         :: p_top                   ! wrf top reference pressure (Pa)
      real, intent(in)         :: ps_wrf(nx,ny)           ! wrf surface pressure (Pa)
      real, intent(in)         :: znu(nz)                 ! sigma coordinates
      real, intent(out)        :: wrfxs(ny,nz,nw)         ! wrfchem vmr(ppm)
      real, intent(out)        :: wrfxe(ny,nz,nw)         ! wrfchem vmr(ppm)
      real, intent(out)        :: wrfys(nx,nz,nw)         ! wrfchem vmr(ppm)
      real, intent(out)        :: wrfye(nx,nz,nw)         ! wrfchem vmr(ppm)

!---------------------------------------------------------------
!     local variables
!---------------------------------------------------------------
      integer :: status, varid
      real, allocatable :: wrk(:)
      real, allocatable :: wrk1(:)
      real, allocatable :: p_moz(:)
      real, allocatable :: p_wrf(:)
      real    :: mozval(nlon,nlat,nlev)
      character(len=20) :: mozspn

      integer :: i, j, k, n

      if( wrf2mz_map(ndx)%moz_cnt > 0 ) then
        allocate( wrk(nlev), wrk1(nz), p_moz(nlev), p_wrf(nz), stat=status )
        if( status /= 0 ) then
          write(*,*) 'bc_interpolated4d: failed to allocate wrk ... p_wrf; error = ',status
          stop
        end if
!---------------------------------------------------------------
!     read mozart pressure
!---------------------------------------------------------------
        if( ndx == 1 ) then
          call read_mozart_ps( time_type%lo_buf_ndx, time_type%lo_moz_ndx, time_type%ncid_lo )
          if( time_type%t_interp ) then
            call read_mozart_ps( time_type%hi_buf_ndx, time_type%hi_moz_ndx, time_type%ncid_hi )
          endif
          call tinterp_mozart_ps( time_type%dels )
          call hinterp_mozart_ps( nx, ny )
        endif
!---------------------------------------------------------------
!     read mozart species
!---------------------------------------------------------------
        call read_moz_species( ndx, time_type%lo_buf_ndx, time_type%lo_moz_ndx, time_type%ncid_lo )
        if( time_type%t_interp ) then
          call read_moz_species( ndx, time_type%hi_buf_ndx, time_type%hi_moz_ndx, time_type%ncid_hi )
        endif
        call tinterp_moz_species( time_type%dels, mozval )
!---------------------------------------------------------------
!     horizontally interpolate species at boundaries
!---------------------------------------------------------------
!     west
!---------------------------------------------------------------
        do j = 1,ny
          wrk(:) = mozval(ix(1,j,0),jy(1,j,0),:)*ax(1,j,1)*by(1,j,1) &
                 + mozval(ix(1,j,0),jy(1,j,1),:)*ax(1,j,1)*by(1,j,0) &
                 + mozval(ix(1,j,1),jy(1,j,0),:)*ax(1,j,0)*by(1,j,1) &
                 + mozval(ix(1,j,1),jy(1,j,1),:)*ax(1,j,0)*by(1,j,0)
          p_moz(:) = ps_mozi(1,j)*hybm(:) + ps0*hyam(:)
          p_wrf(:) = ps_wrf(1,j)*znu(nz:1:-1) + (1. - znu(nz:1:-1))*p_top
          call vinterp( p_moz, p_wrf, wrk, wrk1, nz, nlev )
          do n = 1,nw
            wrfxs(j,:,n) = wrk1(nz:1:-1)
          end do
        end do
!---------------------------------------------------------------
!     east
!---------------------------------------------------------------
        do j = 1, ny
          wrk(:) = mozval(ix(nx,j,0),jy(nx,j,0),:)*ax(nx,j,1)*by(nx,j,1) &
                 + mozval(ix(nx,j,0),jy(nx,j,1),:)*ax(nx,j,1)*by(nx,j,0) &
                 + mozval(ix(nx,j,1),jy(nx,j,0),:)*ax(nx,j,0)*by(nx,j,1) &
                 + mozval(ix(nx,j,1),jy(nx,j,1),:)*ax(nx,j,0)*by(nx,j,0)
          p_moz(:) = ps_mozi(nx,j)*hybm(:) + ps0*hyam(:)
          p_wrf(:) = ps_wrf(nx,j)*znu(nz:1:-1) + (1. - znu(nz:1:-1))*p_top
          call vinterp( p_moz, p_wrf, wrk, wrk1, nz, nlev )
          do n = 1,nw
            wrfxe(j,:,n) = wrk1(nz:1:-1)
          end do
        end do
!---------------------------------------------------------------
!     north
!---------------------------------------------------------------
        do i = 1,nx
          wrk(:) = mozval(ix(i,ny,0),jy(i,ny,0),:)*ax(i,ny,1)*by(i,ny,1) &
                 + mozval(ix(i,ny,0),jy(i,ny,1),:)*ax(i,ny,1)*by(i,ny,0) &
                 + mozval(ix(i,ny,1),jy(i,ny,0),:)*ax(i,ny,0)*by(i,ny,1) &
                 + mozval(ix(i,ny,1),jy(i,ny,1),:)*ax(i,ny,0)*by(i,ny,0)
          p_moz(:) = ps_mozi(i,ny)*hybm(:) + ps0*hyam(:)
          p_wrf(:) = ps_wrf(i,ny)*znu(nz:1:-1) + (1. - znu(nz:1:-1))*p_top
          call vinterp( p_moz, p_wrf, wrk, wrk1, nz, nlev )
          do n = 1,nw
            wrfye(i,:,n) = wrk1(nz:1:-1)
          end do
        end do
!---------------------------------------------------------------
!     south
!---------------------------------------------------------------
        do i = 1,nx
          wrk(:) = mozval(ix(i,1,0),jy(i,1,0),:)*ax(i,1,1)*by(i,1,1) &
                 + mozval(ix(i,1,0),jy(i,1,1),:)*ax(i,1,1)*by(i,1,0) &
                 + mozval(ix(i,1,1),jy(i,1,0),:)*ax(i,1,0)*by(i,1,1) &
                 + mozval(ix(i,1,1),jy(i,1,1),:)*ax(i,1,0)*by(i,1,0)
          p_moz(:) = ps_mozi(i,1)*hybm(:) + ps0*hyam(:)
          p_wrf(:) = ps_wrf(i,1)*znu(nz:1:-1) + (1. - znu(nz:1:-1))*p_top
          call vinterp( p_moz, p_wrf, wrk, wrk1, nz, nlev )
          do n = 1,nw
            wrfys(i,:,n) = wrk1(nz:1:-1)
          end do
        end do
        deallocate( wrk, wrk1, p_moz, p_wrf )
      else
        do n = 1,nw
          wrfxs(:,:,n) = wrf2mz_map(ndx)%wrf_conc
          wrfxe(:,:,n) = wrf2mz_map(ndx)%wrf_conc
          wrfys(:,:,n) = wrf2mz_map(ndx)%wrf_conc
          wrfye(:,:,n) = wrf2mz_map(ndx)%wrf_conc
        end do
      end if

      if( trim(wrf2mz_map(ndx)%wrf_name) == trim(dbg_species) ) then
         write(*,*) ' '
         write(*,*) 'mozart ',trim(dbg_species),' east bc values'
         write(*,'(10(1p,g15.7))') wrfxe(:,1,1)
      end if

      end subroutine bc_interpolate4d

      subroutine ic_interpolate4d( ndx, conc, ps_wrf, znu, p_top, &
                                   nx, ny, nz )
!---------------------------------------------------------------
!     interpolate four-dimensional field ... 
!---------------------------------------------------------------

      use utils, only : wrf2mz_map

!---------------------------------------------------------------
!     input arguments
!---------------------------------------------------------------
      integer, intent(in)      :: ndx
      integer, intent(in)      :: nx, ny, nz              ! dimensions
      real, intent(in)         :: p_top                   ! wrf top reference pressure (Pa)
      real, intent(in)         :: ps_wrf(nx,ny)           ! wrf surface pressure (Pa)
      real, intent(in)         :: znu(nz)                 ! sigma coordinates
      real, intent(out)        :: conc(nx,ny,nz)          ! wrfchem vmr(ppm)

!---------------------------------------------------------------
!     local variables
!---------------------------------------------------------------
      integer :: i, j, k, n
      integer :: status, varid
      real    :: mozval(nlon,nlat,nlev)
      real, allocatable :: wrk(:)
      real, allocatable :: wrk1(:)
      real, allocatable :: p_moz(:)
      real, allocatable :: p_wrf(:)
      character(len=20) :: mozspn

      if( wrf2mz_map(ndx)%moz_cnt > 0 ) then
        allocate( wrk(nlev), wrk1(nz), p_moz(nlev), p_wrf(nz), stat=status )
        if( status /= 0 ) then
          write(*,*) 'failed to allocate wrk ... p_wrf; error = ',status
          stop
        end if
!---------------------------------------------------------------
!     read mozart pressure
!---------------------------------------------------------------
        if( ndx == 1 ) then
          call read_mozart_ps( time_type%lo_buf_ndx, time_type%lo_moz_ndx, time_type%ncid_lo )
          if( time_type%t_interp ) then
            call read_mozart_ps( time_type%hi_buf_ndx, time_type%hi_moz_ndx, time_type%ncid_hi )
          endif
          call tinterp_mozart_ps( time_type%dels )
          call hinterp_mozart_ps( nx, ny )
        endif
!---------------------------------------------------------------
!     read mozart species
!---------------------------------------------------------------
        call read_moz_species( ndx, time_type%lo_buf_ndx, time_type%lo_moz_ndx, time_type%ncid_lo )
        if( time_type%t_interp ) then
          call read_moz_species( ndx, time_type%hi_buf_ndx, time_type%hi_moz_ndx, time_type%ncid_hi )
        endif
        call tinterp_moz_species( time_type%dels, mozval )
!---------------------------------------------------------------
!     horizontally interpolate species
!---------------------------------------------------------------
        do j = 1,ny
          do i = 1,nx
            wrk(:)  = mozval(ix(i,j,0),jy(i,j,0),:)*ax(i,j,1)*by(i,j,1) &
                    + mozval(ix(i,j,0),jy(i,j,1),:)*ax(i,j,1)*by(i,j,0) &
                    + mozval(ix(i,j,1),jy(i,j,0),:)*ax(i,j,0)*by(i,j,1) &
                    + mozval(ix(i,j,1),jy(i,j,1),:)*ax(i,j,0)*by(i,j,0)
            p_moz(:) = ps_mozi(i,j)*hybm(:) + ps0*hyam(:)
            p_wrf(:) = ps_wrf(i,j)*znu(nz:1:-1) + (1. - znu(nz:1:-1))*p_top
            call vinterp( p_moz, p_wrf, wrk, wrk1, nz, nlev )
            conc(i,j,:) = wrk1(nz:1:-1)
            if( trim(wrf2mz_map(ndx)%wrf_name) == trim(dbg_species) ) then
              if( i == 91 .and. j == 93 ) then
                write(*,*) 'interpolation points'
                write(*,*) 'ix = ',ix(i,j,:)
                write(*,*) 'jy = ',jy(i,j,:)
                write(*,*) 'interpolation wghts'
                write(*,*) 'ax = ',ax(i,j,:)
                write(*,*) 'by = ',by(i,j,:)
                write(*,*) 'mozart ',trim(dbg_species),' values @ i,j =',ix(i,j,0),jy(i,j,0)
                write(*,'(10(1p,g15.7))') mozval(ix(i,j,0),jy(i,j,0),:)
                write(*,*) 'horiz interpolated mozart ',trim(dbg_species),' values'
                write(*,'(10(1p,g15.7))') wrk
                write(*,*) 'vert  interpolated mozart ',trim(dbg_species),' values'
                write(*,'(10(1p,g15.7))') wrk1
              end if
            end if
          end do
        end do
        if( trim(wrf2mz_map(ndx)%wrf_name) == trim(dbg_species) ) then
          write(*,*) ' '
          write(*,*) 'mozart ',trim(dbg_species),' ic values'
          write(*,'(10(1p,g15.7))') conc(nx,:,1)
        end if
        deallocate( wrk, wrk1, p_moz, p_wrf )
      else
        conc(:,:,:) = wrf2mz_map(ndx)%wrf_conc
      end if

      end subroutine ic_interpolate4d

      subroutine read_moz_species( sndx, bndx, tndx, ncid )
!---------------------------------------------------------------
!     read mozart species
!---------------------------------------------------------------

      use utils, only : wrf2mz_map

!---------------------------------------------------------------
!     dummy arguments
!---------------------------------------------------------------
      integer, intent(in) :: sndx              ! species index
      integer, intent(in) :: bndx              ! buffer index
      integer, intent(in) :: tndx              ! time index
      integer, intent(in) :: ncid              ! netcdf file index

!---------------------------------------------------------------
!     local variables
!---------------------------------------------------------------
      integer :: n
      integer :: nspec
      integer :: status
      integer :: varid
      real    :: accum(nlon,nlat,nlev)
      real    :: tmpval(nlon,nlat,nlev)
      character(len=20) :: mozspn

      nspec = wrf2mz_map(sndx)%moz_cnt
      if( nspec > 0 ) then
        nstt(1:4) = (/ 1, 1, 1, tndx /)
        ncnt(1:4) = (/ nlon, nlat, nlev, 1 /)
        accum(:,:,:) = 0.
        do n = 1,nspec
          if( wrf2mz_map(sndx)%moz_ext(n) ) then
            mozspn = trim(wrf2mz_map(sndx)%moz_names(n)) // moz_var_suffix
          else
            mozspn = trim(wrf2mz_map(sndx)%moz_names(n))
          end if
          status = nf_inq_varid( ncid, mozspn, varid )
          if( status /= nf_noerr )  then
            write(*,*) 'read_moz_species: failed to get id of ',mozspn
            call handle_error( status )
          end if
          status = nf_get_vara_real( ncid, varid, nstt, ncnt, tmpval )
          if( status /= nf_noerr ) then
            write(*,*) 'read_moz_species: failed to read ',mozspn
           call handle_error( status )
          end if
          if( wrf2mz_map(sndx)%moz_wght(n) == 1. ) then
            accum(:,:,:) = accum(:,:,:) + tmpval(:,:,:)
          else
            accum(:,:,:) = accum(:,:,:) + wrf2mz_map(sndx)%moz_wght(n)*tmpval(:,:,:)
          end if
        end do
        accum(:,:,:) = accum(:,:,:)*wrf2mz_map(sndx)%wrf_wght
      end if

      time_type%conc(:,:,:,bndx) = accum(:,:,:)

      end subroutine read_moz_species

      subroutine tinterp_moz_species( dels, mozval )
!---------------------------------------------------------------
!     time interpolation of mozart species concentration
!---------------------------------------------------------------

!---------------------------------------------------------------
!     dummy arguments
!---------------------------------------------------------------
      real, intent(in)  :: dels                ! linear interp factor
      real, intent(out) :: mozval(nlon,nlat,nlev)

!---------------------------------------------------------------
!     local variables
!---------------------------------------------------------------
      integer :: j, k
      integer :: lo_ndx, hi_ndx
      real    :: delsm1

      lo_ndx = time_type%lo_buf_ndx
      if( time_type%t_interp ) then
        hi_ndx = time_type%hi_buf_ndx
        delsm1 = 1. - dels
        do k = 1,nlev 
          do j = 1,nlat 
            mozval(:,j,k) = time_type%conc(:,j,k,lo_ndx) * dels &
                          + time_type%conc(:,j,k,hi_ndx) * delsm1
          end do
        end do
      else
        do k = 1,nlev 
          do j = 1,nlat 
            mozval(:,j,k) = time_type%conc(:,j,k,lo_ndx)
          end do
        end do
      endif

      end subroutine tinterp_moz_species

      subroutine handle_error( status )
!---------------------------------------------------------------
!     handle errors produced by calling netCDF functions
!---------------------------------------------------------------

!---------------------------------------------------------------
!     dummy arguments :
!---------------------------------------------------------------
      integer, intent(in) :: status

!---------------------------------------------------------------
!     print the error information from processing NETcdf file
!---------------------------------------------------------------
      write(*,*) nf_strerror( status )

!---------------------------------------------------------------
!     exit from the bconLib
!---------------------------------------------------------------
      call exit_mozart_lib( flag=1 )

      end subroutine handle_error

      subroutine exit_mozart_lib( flag )
!---------------------------------------------------------------
!     exit from module_mozart_lib 
!---------------------------------------------------------------

!---------------------------------------------------------------
!     dummy arguments
!---------------------------------------------------------------
      integer, optional, intent(in) :: flag

!---------------------------------------------------------------
!     local variables
!---------------------------------------------------------------
      integer :: status

!---------------------------------------------------------------
!     deallocate
!---------------------------------------------------------------
      if( allocated(ix) ) deallocate( ix )
      if( allocated(jy) ) deallocate( jy )
      if( allocated(ax) ) deallocate( ax )
      if( allocated(by) ) deallocate( by )
      if( associated(time_type%conc) ) deallocate( time_type%conc )
      if( associated(time_type%ps) ) deallocate( time_type%ps )

!-----------------------------------------------------------------------
!     close netCDF file
!-----------------------------------------------------------------------
      if( time_type%ncid_lo /= 0 ) then
        status = nf_close( time_type%ncid_lo )
        if( time_type%ncid_hi /= 0 .and. &
            (time_type%ncid_hi /= time_type%ncid_lo) ) then
          status = nf_close( time_type%ncid_hi )
        end if
      end if

!-----------------------------------------------------------------------
!     output information
!-----------------------------------------------------------------------
      if( present(flag) ) then
        select case( flag )
          case( 1 ); write(*,*) 'fail to process netCDF file...'
          case default; write(*,*) 'unknown error(s) occurred ...'
        endselect
        stop ' in module_mozart_lib ...'
      else
        write(*,*) 'successfully exited from module_mozart_lib ...'
      end if

      end subroutine exit_mozart_lib

      subroutine vinterp( p_moz, p_wrf, src, interp, nz, nlev )
!-----------------------------------------------------------------------
!   	... vertically interpolate input data
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!   	... dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) :: nz, nlev
      real, intent(in)    :: p_moz(nlev)
      real, intent(in)    :: p_wrf(nz)
      real, intent(in)    :: src(nlev)
      real, intent(out)   :: interp(nz)

!-----------------------------------------------------------------------
!   	... local variables
!-----------------------------------------------------------------------
      integer :: i
      integer :: k, kl, ku
      real    :: delp, pinterp

level_loop : &
      do k = 1,nz
         pinterp = p_wrf(k)
         if( pinterp <= p_moz(1) ) then
            interp(k) = src(1)
         else if( pinterp > p_moz(nlev) ) then
            interp(k) = src(nlev)
         else
            do ku = 2,nlev
               if( pinterp <= p_moz(ku) ) then
                  kl = ku - 1
                  delp = log( pinterp/p_moz(kl) ) &
                         / log( p_moz(ku)/p_moz(kl) )
                  interp(k) = src(kl) + delp * (src(ku) - src(kl))
                  exit
               end if
            end do
         end if
      end do level_loop

      end subroutine vinterp

      integer function lotim( cdate, csec, date, datesec, ntim )
!-----------------------------------------------------------------------
! 	... return the index of the time sample that is the lower
!           bound of the interval that contains the input date.  if
!           (cdate,csec) is earlier than the first time sample then 0 is
!           returned.  if (cdate,csec) is later than the last time sample then
!           -index is returned.  if (cdate,csec) is equal to the date of a
!           dynamics time sample then that index is returned.
!-----------------------------------------------------------------------

      use mo_calendar,  only : diffdat

      implicit none

!-----------------------------------------------------------------------
! 	... dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) :: cdate    ! date in yyyymmdd
      integer, intent(in) :: csec     ! seconds relative to date
      integer, intent(in) :: ntim     ! number times
      integer, intent(in) :: date(ntim)
      integer, intent(in) :: datesec(ntim)

!-----------------------------------------------------------------------
! 	... local variables
!-----------------------------------------------------------------------
      integer :: n
      real    :: dtime

!-----------------------------------------------------------------------
!     	... find latest date that is earlier than or equal to (date,sec)
!-----------------------------------------------------------------------
      do n = 1,ntim
        dtime = diffdat( cdate, csec, date(n), datesec(n) )
        if( dtime > 0. ) then
          lotim = n - 1
          if( lotim > 0 ) then
            time_type%dels = dtime/diffdat( date(lotim), datesec(lotim), date(n), datesec(n) )
            time_type%t_interp = time_type%dels /= 1.
          else
            time_type%dels = dtime
          endif
          exit
        endif
      end do

      if( n > ntim ) then
        if( dtime == 0. ) then
          lotim = ntim
          time_type%dels = 1.
          time_type%t_interp = .false.
        else
          lotim = -ntim
          time_type%dels = dtime
        endif
      endif

      end function lotim

      subroutine chk_moz_vars( nspec, filenm )
!---------------------------------------------------------------
!     check wrf to mozart variable mapping
!---------------------------------------------------------------

      use utils,       only : wrf2mz_map

!---------------------------------------------------------------
!     dummy arguments
!---------------------------------------------------------------
      integer, intent(in) :: nspec
      character(len=*), intent(in) :: filenm

!---------------------------------------------------------------
!     local variables
!---------------------------------------------------------------
      integer :: i, n
      integer :: vid
      integer :: status
      character(len=32) :: spcnam

      do n = 1,nspec
         write(*,*) 'checking wrf variable ',trim(wrf2mz_map(n)%wrf_name)
         do i = 1,wrf2mz_map(n)%moz_cnt
            spcnam = ' '
            if( wrf2mz_map(n)%moz_ext(i) ) then
               spcnam = trim(wrf2mz_map(n)%moz_names(i)) // trim(moz_var_suffix)
            else
               spcnam = trim(wrf2mz_map(n)%moz_names(i))
            end if
            write(*,*) 'len ',trim(wrf2mz_map(n)%moz_names(i)),' = ',len_trim(wrf2mz_map(n)%moz_names(i))
            status = nf_inq_varid( ncid, trim(spcnam), vid )
            if( status /= nf_noerr ) then
               write(*,*) 'chk_moz_vars: could not find ',spcnam,' in ',trim(filenm)
               call handle_error( status )
            end if
         end do
      end do

      end subroutine chk_moz_vars

      end module module_mozart_lib
