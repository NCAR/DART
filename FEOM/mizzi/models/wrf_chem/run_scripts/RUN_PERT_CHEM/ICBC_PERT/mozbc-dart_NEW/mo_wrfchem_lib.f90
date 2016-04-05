
      module module_wrfchem_lib

      implicit none

!------------------------------------------------------------------
!     include files
!------------------------------------------------------------------
      include 'netcdf.inc'

!------------------------------------------------------------------
!     public  procedures
!------------------------------------------------------------------
      public  :: init_wrfchem_lib, exit_wrfchem_lib
      public  :: wrfchem_readscalar, wrfchem_read2d,  wrfchem_read3d
      public  :: wrfchem_bc_write4d
      public  :: wrfchem_ic_write4d
      public  :: wrf_ps

!------------------------------------------------------------------
!     private procedures
!------------------------------------------------------------------
      private :: handle_error

!------------------------------------------------------------------
!     public variables
!------------------------------------------------------------------
      integer, public  :: nx, ny, nz, nw, ntime      ! data dimension length
      integer, private :: nscalar, nstring           ! data dimension length
      integer, public  :: year_start, month_start, day_start
      integer, public  :: hour_start, minute_start, second_start
      integer, public  :: year_end, month_end, day_end
      integer, public  :: hour_end, minute_end, second_end
      integer, public  :: dt                        ! temporal resolution
      integer, public  :: domain = 1                ! domain number
      integer, public  :: map_proj
      integer, target, private :: ncids(2)
      integer, allocatable :: wrf_date(:)
      integer, allocatable :: wrf_datesec(:)
      real, public     :: cen_lon
      real, public     :: cen_lat
      real, public     :: stand_lon
      real, public     :: truelat1
      real, public     :: truelat2
      real, public     :: dx
      real, allocatable :: wrf_times(:)
      character(len=32) :: met_file_prefix = 'met_em'
      character(len=8)  :: met_file_suffix = '.nc'
      character(len=16) :: surf_press_name = 'PSFC'
      character(len=3)  :: domain_name = 'd01'
      character(len=1)  :: met_file_separator = '.'

!------------------------------------------------------------------
!     variables for reading netCDF data
!------------------------------------------------------------------
      integer, pointer, private :: ncid_bc                     ! netCDF file ID for bndy file
      integer, pointer, private :: ncid_ic                     ! netCDF file ID for ic file
      integer, private :: nstt(4), ncnt(4)                     ! start index and counter array
      character(len=5), private :: bc_suffix(4,2)
      character(len=3), private :: bc_memord(4)

      contains

      subroutine init_wrfchem_lib( wrfchem_bdy_fn, wrfchem_input_fn, wrf_dir, nspec, &
                                   def_missing_var, do_bc, do_ic )
!------------------------------------------------------------------
!     initialize netcdf file
!------------------------------------------------------------------

      use utils, only : wrf2mz_map

!------------------------------------------------------------------
!     dummy arguments
!------------------------------------------------------------------
      integer, intent(in) :: nspec
      logical, intent(in) :: do_bc
      logical, intent(in) :: do_ic
      logical, intent(in) :: def_missing_var
      character(len=*), intent(in) ::  wrfchem_bdy_fn, wrfchem_input_fn
      character(*), intent(in) ::  wrf_dir

!------------------------------------------------------------------
!     local variables
!------------------------------------------------------------------
      integer, parameter :: ic = 1
      integer, parameter :: bc = 2
      character(len=19)   :: proj_name(0:3) = (/ 'LATLON             ', 'LAMBERT            ', &
                                                 'POLAR STEREOGRAPHIC', 'MERCATOR           ' /)

      integer :: status
      integer :: ncid
      integer :: vid
      integer :: open_option
      integer :: att_type
      integer :: n, ns, i
      integer :: l, m
      integer :: file
      integer :: xdimid(2),  ydimid(2), zdimid(2), wdimid(2), tdimid(2)
      integer :: dims(4)
      character(len=64) :: attribute
      character(len=64) :: filespec(2)
      character(len=32) :: varname
      character(len=32) :: dtstring
      character(len=19) :: tstring
      character(len=3)  :: memorder
      logical :: is_gas
      logical :: missing(nspec)

      ncid_ic => ncids(ic)
      ncid_bc => ncids(bc)
      ncids(:) = 0
!------------------------------------------------------------------
!     open the wrfinput or wrfbdy NETCDF file
!------------------------------------------------------------------
file_loop1 : &
      do file = ic,bc
        open_option = nf_write
        if( file == ic .and. do_ic ) then
          filespec(file) = wrfchem_input_fn
        else if( file == bc .and. do_bc ) then
          filespec(file) = wrfchem_bdy_fn
        else
          cycle
        end if
        status = nf_open( filespec(file), open_option, ncids(file) )
        if( status /= nf_noerr ) then
          write(*,*) 'Failed to open ',trim(filespec(file)),' with open_option = ',open_option
          call handle_error( status )
        else
          write(*,*) 'Opened ',trim(filespec(file))
        end if
!------------------------------------------------------------------
!     get times from wrfbdy file and convert to date,sec
!------------------------------------------------------------------
        call read_wrf_times( file )
        write(*,*) ' '
        write(*,*) 'wrf times in ',trim(filespec(file))
        do n = 1,ntime
          write(*,*) n,wrf_date(n),wrf_datesec(n)
        end do
      end do file_loop1

!------------------------------------------------------------------
!     check wrf ps file for existence
!------------------------------------------------------------------
      if( .not. check_wrfps( wrf_dir, wrf_date(1), wrf_datesec(1) ) ) then
        stop
      end if
      if( ntime > 1 ) then
        if( .not. check_wrfps( wrf_dir, wrf_date(ntime), wrf_datesec(ntime) ) ) then
           stop
        end if
      end if

file_loop2 : &
      do file = ic,bc
        ncid = ncids(file)
        if( ncid > 0 ) then
!------------------------------------------------------------------
!     get spacial dimension lengths from wrfinput file
!------------------------------------------------------------------
          status = nf_inq_dimid( ncid, 'west_east', xdimid(file) )
          if( status /= nf_noerr ) call handle_error( status )

          status = nf_inq_dimlen( ncid, xdimid(file), nx )
          if( status /= nf_noerr ) call handle_error( status )

          status = nf_inq_dimid( ncid, 'south_north', ydimid(file) )
          if( status /= nf_noerr ) call handle_error( status )

          status = nf_inq_dimlen( ncid, ydimid(file), ny )
          if( status /= nf_noerr ) call handle_error( status )

          status = nf_inq_dimid( ncid, 'bottom_top', zdimid(file) )
          if( status /= nf_noerr ) call handle_error( status )

          status = nf_inq_dimlen( ncid, zdimid(file), nz )
          if( status /= nf_noerr ) call handle_error( status )

          status = nf_inq_dimid( ncid, 'Time', tdimid(file) )
          if( status /= nf_noerr ) call handle_error( status )

          if( file == bc ) then
            status = nf_inq_dimid( ncid, 'bdy_width', wdimid(file) )
            if( status /= nf_noerr ) call handle_error( status )

            status = nf_inq_dimlen( ncid, wdimid(file), nw )
            if( status /= nf_noerr ) call handle_error( status )

            bc_suffix(:,1) = (/ '_BXS ', '_BXE ', '_BYS ', '_BYE ' /)
            bc_suffix(:,2) = (/ '_BTXS', '_BTXE', '_BTYS', '_BTYE' /)
            bc_memord(:)   = (/ 'XSZ', 'XEZ', 'YSZ', 'YEZ' /)
          end if
        else
          cycle file_loop2
        end if

        missing(:nspec) = .false.
!------------------------------------------------------------------
!     check wrfinput file for wrf names
!------------------------------------------------------------------
        do n = 1,nspec
          if( file == ic ) then
            status = nf_inq_varid( ncid, trim(wrf2mz_map(n)%wrf_name), vid )
            if( status /= nf_noerr ) then
              missing(n) = .true.
            end if
          elseif( file == bc ) then
            do l = 1,2
              do m = 1,4
                varname = trim( wrf2mz_map(n)%wrf_name ) // trim( bc_suffix(m,l) )
                status = nf_inq_varid( ncid, trim(varname), vid )
                if( status /= nf_noerr ) then
                  missing(n) = .true.
                endif
              end do
            end do
          endif
        end do

        if( .not. def_missing_var .and. any( missing(:) ) ) then
          write(*,*) ' '
          write(*,*) 'init_wrfchem: could not find the following variables in ',trim(filespec(file))
          write(*,*) ' '
          do n = 1,nspec
            if( missing(n) ) then
              write(*,*) trim(wrf2mz_map(n)%wrf_name)
            endif
          end do
          write(*,*) ' '
          call handle_error( status )
        end if

any_missing : &
        if( any( missing(:) ) ) then
          status = nf_redef( ncid )
          if( status /= nf_noerr ) then
            write(*,*) 'init_wrfchem: could not put ',trim(filespec(file)),' in redefine mode'
            call handle_error( status )
          end if
          if( file == ic ) then
            dims(:) = (/ xdimid(file), ydimid(file), zdimid(file), tdimid(file) /)
          else
            dims(:) = (/ xdimid(file), zdimid(file), wdimid(file), tdimid(file) /)
          endif
          do n = 1,nspec
            if( missing(n) ) then
              varname = trim(wrf2mz_map(n)%wrf_name)
              is_gas = wrf2mz_map(n)%wrf_wght == 1.e6
              if( file == ic ) then
                write(*,*) 'inserting ',trim(varname)
                call insert_var( file, varname, 'XYZ', dims, is_gas )
              else
                do l = 1,2
                  do m = 1,4
                    if( m <= 2 ) then
                      dims(1) = ydimid(file)
                    else
                      dims(1) = xdimid(file)
                    endif
                    varname = trim(wrf2mz_map(n)%wrf_name) // trim(bc_suffix(m,l))
                    write(*,*) 'inserting ',trim(varname)
                    call insert_var( file, trim(varname), bc_memord(m), dims, is_gas )
                  end do
                end do
              endif
            endif
          end do
          status = nf_enddef( ncid )
          if( status /= nf_noerr ) then
            write(*,*) 'init_wrfchem: could not take ',trim(filespec(file)),' out of redefine mode'
            call handle_error( status )
          end if
        end if any_missing
      end do file_loop2

!------------------------------------------------------------------
!     get map projection variables
!------------------------------------------------------------------
      if( ncid_ic > 0 ) then
        ncid = ncid_ic
      else
        ncid = ncid_bc
      endif
      status = nf_get_att_int( ncid, nf_global, 'MAP_PROJ', map_proj )
      if( status /= nf_noerr ) then
         write(*,*) 'init_wrfchem: could not get MAP_PROJ'
         call handle_error( status )
      else
         write(*,*) 'init_wrfchem: MAP_PROJ = ',trim(proj_name(map_proj))
      end if
      status = nf_get_att_real( ncid, nf_global, 'CEN_LON', cen_lon )
      if( status /= nf_noerr ) then
         write(*,*) 'init_wrfchem: could not get CEN_LON'
         call handle_error( status )
      else
         write(*,*) 'init_wrfchem: CEN_LON = ',cen_lon
      end if
      status = nf_get_att_real( ncid, nf_global, 'CEN_LAT', cen_lat )
      if( status /= nf_noerr ) then
         write(*,*) 'init_wrfchem: could not get CEN_LAT'
         call handle_error( status )
      else
         write(*,*) 'init_wrfchem: CEN_LAT = ',cen_lat
      end if
      status = nf_get_att_real( ncid, nf_global, 'STAND_LON', stand_lon )
      if( status /= nf_noerr ) then
         write(*,*) 'init_wrfchem: could not get STAND_LON'
         call handle_error( status )
      else
         write(*,*) 'init_wrfchem: STAND_LON = ',stand_lon
      end if
      status = nf_get_att_real( ncid, nf_global, 'TRUELAT1', truelat1 )
      if( status /= nf_noerr ) then
         write(*,*) 'init_wrfchem: could not get TRUELAT1'
         call handle_error( status )
      else
         write(*,*) 'init_wrfchem: TRUELAT1 = ',truelat1
      end if
      status = nf_get_att_real( ncid, nf_global, 'TRUELAT2', truelat2 )
      if( status /= nf_noerr ) then
         write(*,*) 'init_wrfchem: could not get TRUELAT2'
         call handle_error( status )
      else
         write(*,*) 'init_wrfchem: TRUELAT2 = ',truelat2
      end if
      status = nf_get_att_real( ncid, nf_global, 'DX', dx )
      if( status /= nf_noerr ) then
         write(*,*) 'init_wrfchem: could not get DX'
         call handle_error( status )
      else
         write(*,*) 'init_wrfchem: DX = ',dx
      end if

      write(*,*) 'finished init_wrfchem_lib'

      end subroutine init_wrfchem_lib

      subroutine read_wrf_times( file )
!------------------------------------------------------------------
!     read wrf times in wrfinput file
!------------------------------------------------------------------

      use utils,       only : wrf2mz_time
      use mo_calendar, only : addsec2dat, diffdat

!------------------------------------------------------------------
!     dummy arguments
!------------------------------------------------------------------
      integer, intent(in) :: file

!------------------------------------------------------------------
!     local variables
!------------------------------------------------------------------
      integer, parameter :: bc = 2
      integer :: i
      integer :: ncid
      integer :: sdimid
      integer :: tdimid
      integer :: nstring
      integer :: delta_secs
      integer :: vid
      integer :: status
      character(len=19) :: dtstring

      ncid = ncids(file)
!------------------------------------------------------------------
!     get the time dimension length
!------------------------------------------------------------------
      status = nf_inq_dimid( ncid, 'Time', tdimid )
      if( status /= nf_noerr )  call handle_error( status )

      status = nf_inq_dimlen( ncid, tdimid, ntime )
      if( status /= nf_noerr ) call handle_error( status )

      if( allocated( wrf_date ) ) then
         deallocate( wrf_date )
      end if
      if( allocated( wrf_datesec ) ) then
         deallocate( wrf_datesec )
      end if
      allocate( wrf_date(ntime+1), wrf_datesec(ntime+1), stat=status )
      if( status /= 0 ) then
         write(*,*) 'read_wrf_times: failed to allocate wrf_date, wrf_datesec; error = ',status
         stop
      end if

!------------------------------------------------------------------
!     get the string dimension length
!------------------------------------------------------------------
      status = nf_inq_dimid( ncid, 'DateStrLen', sdimid )
      if( status /= nf_noerr )  call handle_error( status )

      status = nf_inq_dimlen( ncid, sdimid, nstring )
      if( status /= nf_noerr ) call handle_error( status )

!------------------------------------------------------------------
!     read times
!------------------------------------------------------------------
      status = nf_inq_varid( ncid, 'Times', vid )
      if( status /= nf_noerr )  call handle_error( status )

      ncnt(1:2) = (/ nstring, 1 /)
      nstt(1)   = 1
      do i = 1,ntime
         nstt(2) = i
         status = nf_get_vara_text( ncid, vid, nstt(1:2), ncnt(1:2), dtstring )
         if( status /= nf_noerr ) then
            call handle_error( status )
         end if
         call wrf2mz_time( dtstring, wrf_date(i), wrf_datesec(i) )
      end do

      if( file == bc ) then
         ntime = ntime + 1
         delta_secs = int( 86400. * diffdat( wrf_date(1), wrf_datesec(1), wrf_date(2), wrf_datesec(2) ) )
         write(*,*) 'read_wrf_times: delta_secs = ',delta_secs
         wrf_date(ntime)    = wrf_date(ntime-1)
         wrf_datesec(ntime) = wrf_datesec(ntime-1)
         call addsec2dat( delta_secs, wrf_date(ntime), wrf_datesec(ntime) )
      end if

      end subroutine read_wrf_times

      logical function check_wrfps( wrf_dir, date, secs )
!------------------------------------------------------------------
!     check met_em file for existence
!------------------------------------------------------------------

      use utils, only : mz2wrf_time

!------------------------------------------------------------------
!     dummy arguments
!------------------------------------------------------------------
      integer, intent(in) :: date
      integer, intent(in) :: secs
      character(len=*), intent(in) :: wrf_dir

!------------------------------------------------------------------
!     local variables
!------------------------------------------------------------------
      logical :: lexist
      character(len=132) :: filenm
      character(len=19)  :: tstring

      call mz2wrf_time( tstring, date, secs )
      if( len_trim( met_file_suffix ) > 0 ) then
         filenm = trim( wrf_dir ) // trim( met_file_prefix) // trim(met_file_separator) // domain_name // trim(met_file_separator) // tstring // trim(met_file_suffix)
      else
         filenm = trim( wrf_dir ) // trim( met_file_prefix) // trim(met_file_separator) // domain_name // trim(met_file_separator) // tstring
      end if
      inquire( exist=lexist, file = trim( filenm ) )
      check_wrfps = lexist
      if( .not. lexist ) then
         write(*,*) 'check_wrfps: file ',trim(filenm)
         write(*,*) '             does not exist'
      end if

      end function check_wrfps

      subroutine wrfchem_readscalar( vname, val, filespec )
!------------------------------------------------------------------
!     read scalar data 
!------------------------------------------------------------------

!------------------------------------------------------------------
!     dummy arguments
!------------------------------------------------------------------
      character(*), intent(in) :: vname
      character(*), optional, intent(in) :: filespec
      real, intent(out)        :: val

!------------------------------------------------------------------
!     local variables
!------------------------------------------------------------------
      integer, parameter :: ic = 1
      integer :: status
      integer :: varid
      integer :: ncid

      if( present( filespec ) ) then
        status = nf_open( trim(filespec), nf_nowrite, ncid )
        if( status /= nf_noerr ) then
          write(*,*) 'wrfchem_readscalar: failed to open ',trim(filespec)
          call handle_error( status )
        else
          write(*,*) 'wrfchem_readscalar: Opened ',trim(filespec)
        end if
      else
        ncid = ncids(ic)
      endif
!------------------------------------------------------------------
!     get variable ID, if no such a variable, return
!------------------------------------------------------------------
      status = nf_inq_varid( ncid, vname, varid )
      if( status /= nf_noerr ) call handle_error( status )

!------------------------------------------------------------------
!     read value
!------------------------------------------------------------------
      status = nf_get_var_real( ncid, varid, val )
      if( status /= nf_noerr ) call handle_error( status )

      if( present( filespec ) ) then
        status = nf_close( ncid )
      endif

      end subroutine wrfchem_readscalar

      subroutine wrfchem_read2d( vname, val, filespec )
!------------------------------------------------------------------
!     read one-dimensional data
!------------------------------------------------------------------

!------------------------------------------------------------------
!     dummy arguments
!------------------------------------------------------------------
      character(*), intent(in) :: vname
      character(*), optional, intent(in) :: filespec
      real, intent(out)        :: val(nz)

!------------------------------------------------------------------
!     local variables
!------------------------------------------------------------------
      integer, parameter :: ic = 1
      integer :: status
      integer :: varid
      integer :: ncid

      if( present( filespec ) ) then
        status = nf_open( trim(filespec), nf_nowrite, ncid )
        if( status /= nf_noerr ) then
          write(*,*) 'wrfchem_read2d: failed to open ',trim(filespec)
          call handle_error( status )
        else
          write(*,*) 'wrfchem_read2d: Opened ',trim(filespec)
        end if
      else
        ncid = ncids(ic)
      endif
!------------------------------------------------------------------
!     get variable ID, if no such a variable, return
!------------------------------------------------------------------
      status = nf_inq_varid( ncid, vname, varid )
      if( status /= nf_noerr ) call handle_error( status )

!------------------------------------------------------------------
!     read variable
!------------------------------------------------------------------
      nstt(1:2) = (/ 1, 1 /)
      ncnt(1:2) = (/ nz, 1 /)
      status = nf_get_vara_real( ncid, varid, nstt(1:2), ncnt(1:2), val(:) )
      if( status /= nf_noerr ) call handle_error( status )

      if( present( filespec ) ) then
        status = nf_close( ncid )
      endif

      end subroutine wrfchem_read2d

      subroutine wrfchem_read3d( vname, val, filespec )
!------------------------------------------------------------------
!     read three-dimensional data
!------------------------------------------------------------------

!------------------------------------------------------------------
!     dummy arguments
!------------------------------------------------------------------
      character(len=*), intent(in) :: vname
      character(len=*), optional, intent(in) :: filespec
      real, intent(out)        :: val(nx,ny)

!------------------------------------------------------------------
!     local variables
!------------------------------------------------------------------
      integer, parameter :: ic = 1
      integer :: status
      integer :: varid
      integer :: ncid

      if( present( filespec ) ) then
        status = nf_open( trim(filespec), nf_nowrite, ncid )
        if( status /= nf_noerr ) then
          write(*,*) 'wrfchem_read2d: failed to open ',trim(filespec)
          call handle_error( status )
        else
          write(*,*) 'wrfchem_read2d: Opened ',trim(filespec)
        end if
      else
        ncid = ncids(ic)
      endif
!------------------------------------------------------------------
!     get variable ID, if no such a variable, return
!------------------------------------------------------------------
      status = nf_inq_varid( ncid, vname, varid )
      if( status /= nf_noerr ) call handle_error( status )

!------------------------------------------------------------------
!     read value
!------------------------------------------------------------------
      nstt(1:3) = (/ 1, 1, 1 /)
      ncnt(1:3) = (/ nx, ny, 1 /)
      status = nf_get_vara_real( ncid, varid, nstt(1:3), ncnt(1:3), val(:,:) )
      if( status /= nf_noerr ) call handle_error( status )

      if( present( filespec ) ) then
        status = nf_close( ncid )
      endif

      end subroutine wrfchem_read3d

      subroutine wrfchem_bc_write4d( ndx, valxs, valxe, valys, valye, &
                                     vatxs, vatxe, vatys, vatye, it, &
                                     write_conc, write_tend )
!------------------------------------------------------------------
!     write four-dimensional data
!------------------------------------------------------------------
      
      use utils, only : wrf2mz_map

!------------------------------------------------------------------
!     dummy arguments
!------------------------------------------------------------------
      integer, intent(in)      :: ndx
      integer, intent(in)      :: it
      real, intent(in)         :: valxs(ny,nz,nw)
      real, intent(in)         :: valxe(ny,nz,nw)
      real, intent(in)         :: valys(nx,nz,nw)
      real, intent(in)         :: valye(nx,nz,nw)
      real, intent(in)         :: vatxs(ny,nz,nw)
      real, intent(in)         :: vatxe(ny,nz,nw)
      real, intent(in)         :: vatys(nx,nz,nw)
      real, intent(in)         :: vatye(nx,nz,nw)
      logical, intent(in)      :: write_conc
      logical, intent(in)      :: write_tend

!------------------------------------------------------------------
!     local variables
!------------------------------------------------------------------
      integer :: status
      integer :: varid
      integer :: i, m
      character(20) :: bcname

      nstt(1:3) = (/ 1, 1, 1 /)
      ncnt(3:4) = (/ nw, 1 /)
is_concentration : &
      if( write_conc ) then
         nstt(4) = it
         do m = 1,4
            bcname = trim(wrf2mz_map(ndx)%wrf_name) // trim(bc_suffix(m,1))
!------------------------------------------------------------------
!     get variable id
!------------------------------------------------------------------
            status = nf_inq_varid( ncid_bc, trim(bcname), varid )
            if( status /= nf_noerr )  then
               write(*,*) 'write4d_bc: failed to ',trim(bcname),' id'
               call handle_error( status )
            end if
!------------------------------------------------------------------
!     set start and count arrays
!------------------------------------------------------------------
            if( m <= 2 ) then
               ncnt(1:2) = (/ ny, nz /)
            else
               ncnt(1:2) = (/ nx, nz /)
            end if
!------------------------------------------------------------------
!     write BC values
!------------------------------------------------------------------
            select case( m )
! west
            case( 1 )
               status = nf_put_vara_real( ncid_bc, varid, nstt, ncnt, valxs )
! east       
            case( 2 )
               status = nf_put_vara_real( ncid_bc, varid, nstt, ncnt, valxe )
! north
            case( 3 )
               status = nf_put_vara_real( ncid_bc, varid, nstt, ncnt, valys )
! south
            case( 4 )
               status = nf_put_vara_real( ncid_bc, varid, nstt, ncnt, valye )
            end select
            if( status /= nf_noerr ) then
               write(*,*) 'write4d_bc: failed to write ',trim(bcname)
               call handle_error( status )     
            end if
         end do
      end if is_concentration
is_trend : &
      if( write_tend ) then
         nstt(4) = it - 1
         do m = 1,4
            bcname = trim(wrf2mz_map(ndx)%wrf_name) // trim(bc_suffix(m,2))
!------------------------------------------------------------------
!     get variable id
!------------------------------------------------------------------
            status = nf_inq_varid( ncid_bc, trim(bcname), varid )
            if( status /= nf_noerr )  then
               write(*,*) 'write4d_bc: failed to ',trim(bcname),' id'
               call handle_error( status )
            end if
!------------------------------------------------------------------
!     set start and count arrays
!------------------------------------------------------------------
            if( m <= 2 ) then
               ncnt(1:2) = (/ ny, nz /)
            else
               ncnt(1:2) = (/ nx, nz /)
            end if
!------------------------------------------------------------------
! Tendencies 
!------------------------------------------------------------------
            select case( m )
            case( 1 )           ! west
               status = nf_put_vara_real( ncid_bc, varid, nstt, ncnt, vatxs )
            case( 2 )           ! east
               status = nf_put_vara_real( ncid_bc, varid, nstt, ncnt, vatxe )
            case( 3 )           ! north
               status = nf_put_vara_real( ncid_bc, varid, nstt, ncnt, vatys )
            case( 4 )           ! south
               status = nf_put_vara_real( ncid_bc, varid, nstt, ncnt, vatye )
            end select
            if( status /= nf_noerr ) then
               call handle_error( status )
               write(*,*) 'write4d_bc: failed to write ',trim(bcname)
            end if
         end do
      end if is_trend

      end subroutine wrfchem_bc_write4d 

      subroutine wrfchem_ic_write4d( ndx, conc, it )
!------------------------------------------------------------------
!     write four-dimensional data
!------------------------------------------------------------------

      use utils, only : wrf2mz_map

!------------------------------------------------------------------
!     dummy arguments
!------------------------------------------------------------------
      integer, intent(in)      :: ndx
      integer, intent(in)      :: it
      real, intent(in)         :: conc(nx,ny,nz)

!------------------------------------------------------------------
!     local variables
!------------------------------------------------------------------
      integer :: status
      integer :: varid

      nstt(:) = (/ 1, 1, 1, 1 /)
      ncnt(:) = (/ nx, ny, nz, 1 /)
!------------------------------------------------------------------
!     get variable id
!------------------------------------------------------------------
      status = nf_inq_varid( ncid_ic, trim(wrf2mz_map(ndx)%wrf_name), varid )
      if( status /= nf_noerr )  then
         write(*,*) 'ic_write4d; failed to get ',trim(wrf2mz_map(ndx)%wrf_name),' id'
         call handle_error( status )
      end if
!------------------------------------------------------------------
!     write IC values
!------------------------------------------------------------------
      status = nf_put_vara_real( ncid_ic, varid, nstt, ncnt, conc )
      if( status /= nf_noerr ) then
         write(*,*) 'ic_write4d; failed to write ',trim(wrf2mz_map(ndx)%wrf_name)
         call handle_error( status )     
      end if

      end subroutine wrfchem_ic_write4d 

      subroutine handle_error( status )
!------------------------------------------------------------------
!     handle errors produced by calling netCDF functions
!------------------------------------------------------------------

!------------------------------------------------------------------
!     dummy arguments :
!------------------------------------------------------------------
      integer, intent(in) :: status

      print*, nf_strerror( status )

!------------------------------------------------------------------
!     exit from the bconLib
!------------------------------------------------------------------
      call exit_wrfchem_lib( flag=1 )

      end subroutine handle_error

      subroutine exit_wrfchem_lib( flag )
!------------------------------------------------------------------
!     exit from bconLib
!------------------------------------------------------------------

!------------------------------------------------------------------
!     dummy arguments
!------------------------------------------------------------------
      integer, optional, intent(in) :: flag

!------------------------------------------------------------------
!     local variables
!------------------------------------------------------------------
      integer :: status

!------------------------------------------------------------------
!     close netCDF files
!------------------------------------------------------------------
      if( ncid_bc /= 0 ) then
         status = nf_close( ncid_bc )
      end if
      if( ncid_ic /= 0 ) then
         status = nf_close( ncid_ic )
      end if

!------------------------------------------------------------------
!     output information
!------------------------------------------------------------------
      if( present(flag) ) then
        select case( flag )
          case( 1 ); print*, 'fail to process netCDF file...'
          case( 2 ); print*, 'no such a species to save ...' 
          case default; print*, 'unknown error(s) occurred ...'
        end select
        stop ' in module_wrfchem_lib ...'
      else
        write(*,*) 'successfully exited from module_wrfchem_lib ...'
      end if 

      end subroutine exit_wrfchem_lib

      subroutine wrf_ps( wrf_dir, date, datesec, ps_wrf )
!---------------------------------------------------------------
!     read the wrf surface pressure
!---------------------------------------------------------------

      use utils, only : mz2wrf_time

!------------------------------------------------------------------
!     dummy arguments
!------------------------------------------------------------------
      integer, intent(in) :: date
      integer, intent(in) :: datesec
      real, intent(out)   :: ps_wrf(nx,ny)
      character(len=*), intent(in) :: wrf_dir

!------------------------------------------------------------------
!     local variables
!------------------------------------------------------------------
      integer :: status
      integer :: varid
      integer :: ncid_bc
      character(len=132) :: filenm
      character(len=19)  :: tstring

      call mz2wrf_time( tstring, date, datesec )
      if( len_trim( met_file_suffix ) > 0 ) then
         filenm = trim( wrf_dir ) // trim( met_file_prefix) // trim(met_file_separator) // domain_name // trim(met_file_separator) // tstring // trim(met_file_suffix)
      else
         filenm = trim( wrf_dir ) // trim( met_file_prefix) // trim(met_file_separator) // domain_name // trim(met_file_separator) // tstring
      end if

      status = nf_open( trim(filenm), nf_write, ncid_bc )
      if( status /= nf_noerr ) then
         write(*,*) 'wrf_ps: failed to open ',trim(filenm)
         call handle_error( status )
      else
         write(*,*) 'wrf_ps: Opened ',trim(filenm)
      end if

      status = nf_inq_varid( ncid_bc, trim(surf_press_name), varid )
      if( status /= nf_noerr )  then
         write(*,*) 'wrf_ps: failed to find PSFC'
         call handle_error( status )
      end if

      nstt(1:3) = (/ 1, 1, 1 /)
      ncnt(1:3) = (/ nx, ny, 1 /)
      status = nf_get_vara_real( ncid_bc, varid, nstt(1:3), ncnt(1:3), ps_wrf )
      if( status /= nf_noerr ) then
         write(*,*) 'wrf_ps: failed to read ',trim(surf_press_name)
         call handle_error( status )
      end if

      status = nf_close( ncid_bc )

      end subroutine wrf_ps

      subroutine insert_var( file, varname, MemOrder, dims, is_gas )
!------------------------------------------------------------------
!     put new variable in netcdf file
!------------------------------------------------------------------

!------------------------------------------------------------------
!     dummy arguments
!------------------------------------------------------------------
      integer, intent(in) :: file
      integer, intent(in) :: dims(:)
      character(len=*), intent(in) :: varname
      character(len=*), intent(in) :: MemOrder
      logical, intent(in) :: is_gas

!------------------------------------------------------------------
!     local variables
!------------------------------------------------------------------
      integer, parameter :: ic = 1
      integer :: ncid
      integer :: status
      integer :: vid
      character(len=64) :: attribute

      ncid = ncids(file)
      status = nf_def_var( ncid, trim(varname), nf_float, 4, dims, vid )
      if( status /= nf_noerr ) then
        write(*,*) 'insert_var: Failed to define ' // trim(varname) // ' variable'
        call handle_error( status )
      endif

      status = nf_put_att_int( ncid, vid, 'FieldType', nf_int, 1, 104 )
      if( status /= nf_noerr ) then
        write(*,*) 'insert_var: Failed to create ' // trim(varname) // ' int attribute'
        call handle_error( status )
      endif

      status = nf_put_att_text( ncid, vid, 'MemoryOrder', len_trim(MemOrder), trim(MemOrder) )
      if( status /= nf_noerr ) then
        write(*,*) 'insert_var: Failed to create ' // trim(varname) // ' MemoryOrder attribute'
        call handle_error( status )
      endif

      attribute = trim(varname) // ' concentration'
      status = nf_put_att_text( ncid, vid, 'description', len_trim(attribute), trim(attribute) )
      if( status /= nf_noerr ) then
        write(*,*) 'insert_var: Failed to create ' // trim(varname) // ' concentration attribute'
        call handle_error( status )
      endif

      if( is_gas ) then
        attribute = 'ppmv'
      else
        attribute = 'ug/kg-dryair'
      endif
      status = nf_put_att_text( ncid, vid, 'units', len_trim(attribute), trim(attribute) )
      if( status /= nf_noerr ) then
        write(*,*) 'insert_var: Failed to create ' // trim(varname) // ' units attribute'
        call handle_error( status )
      endif

      attribute = ' '
      status = nf_put_att_text( ncid, vid, 'stagger', len_trim(attribute), trim(attribute) )
      if( status /= nf_noerr ) then
        write(*,*) 'insert_var: Failed to create ' // trim(varname) // ' stagger attribute'
        call handle_error( status )
      endif

      if( file == ic ) then
        attribute = 'XLONG XLAT'
        status = nf_put_att_text( ncid, vid, 'coordinates', len_trim(attribute), trim(attribute) )
        if( status /= nf_noerr ) then
          write(*,*) 'insert_var: Failed to create ' // trim(varname) // ' coordinates attribute'
          call handle_error( status )
        endif
      endif

      end subroutine insert_var

      end module module_wrfchem_lib
