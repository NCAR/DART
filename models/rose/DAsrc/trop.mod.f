c--------------------------------------------------------------------
      module trop
c--------------------------------------------------------------------

      implicit none

      private
      public :: trop_open, trop_close, trop_retrieve

      integer :: ncid

      integer :: tlen            ! size of time dimension
      integer :: zlen            ! size of altitude dimension
      integer :: ylen            ! size of latitude dimension
      integer :: xlen            ! size of longitude dimension

      real, allocatable :: trop_tim(:)
      real, allocatable :: trop_lev(:)
      real, allocatable :: trop_lat(:)
      real, allocatable :: trop_lon(:)
      real, allocatable :: trop_var(:,:,:,:)
      integer, allocatable :: mtime(:,:)
     
      contains

c--------------------------------------------------------------------
      subroutine trop_open( tfile )
c--------------------------------------------------------------------

      implicit none 
      include 'netcdf.inc'

      character(len=*) :: tfile

      integer :: ncerr, i

      integer :: tid, zid, yid, xid, mtimeid
      integer, dimension(2) ::  start, count

      print *, 'trop_open()'

c... open netcdf file and inquire about the size of the model dimensions

      ncerr = nf_open( tfile, NF_NOWRITE, ncid )
      print *, 'opening with '//nf_strerror(ncerr)

      ncerr = nf_inq_dimid( ncid, 'time', tid)
      ncerr = nf_inq_dimid( ncid, 'lev', zid)
      ncerr = nf_inq_dimid( ncid, 'lat', yid)
      ncerr = nf_inq_dimid( ncid, 'lon', xid)
      ncerr = nf_inq_dimid( ncid, 'mtime', mtimeid)

      ncerr = nf_inq_dimlen( ncid, tid, tlen)
      ncerr = nf_inq_dimlen( ncid, zid, zlen)
      ncerr = nf_inq_dimlen( ncid, yid, ylen)
      ncerr = nf_inq_dimlen( ncid, xid, xlen)
      print *, 'trop array lengths', tlen, zlen, ylen, xlen

c... allocate memory for the time, alt, lat, and constituent arrays

      allocate( trop_tim(tlen) )
      allocate( trop_lev(zlen) )
      allocate( trop_lat(ylen) )
      allocate( trop_lon(xlen) )
      allocate( trop_var(zlen,xlen,ylen,tlen) )
      allocate( mtime(2,tlen) )

      ncerr = nf_inq_varid( ncid, 'time', tid )
      ncerr = nf_inq_varid( ncid, 'lev', zid )
      ncerr = nf_inq_varid( ncid, 'lat', yid )
      ncerr = nf_inq_varid( ncid, 'lon', xid )
      ncerr = nf_inq_varid( ncid, 'mtime', mtimeid )

c... copy the time, alt, and lat arrays from the netcdf file

      ncerr = nf_get_vara_real( ncid, tid, 1, tlen, trop_tim)
      ncerr = nf_get_vara_real( ncid, zid, 1, zlen, trop_lev)
      ncerr = nf_get_vara_real( ncid, yid, 1, ylen, trop_lat)
      ncerr = nf_get_vara_real( ncid, xid, 1, xlen, trop_lon)
      start(1) = 1
      start(2) = 1
      count(1) = 2
      count(2) = tlen 
      ncerr = nf_get_vara_int( ncid, mtimeid, start, count, mtime ) 

      print 40,trop_tim
 40   format('time',8f8.0)

      print 45,trop_lev
 45   format('lev ',8f8.1)

      print 50,trop_lat
 50   format('lat ',8f8.1)

      print 55,trop_lon
 55   format('lon ',8f8.1)

      print 60,(mtime(1,i),mtime(2,i),i=1,tlen)
 60   format('mtime ',10i6)

      end subroutine trop_open 

c--------------------------------------------------------------------
      subroutine trop_var_get( varname )
c--------------------------------------------------------------------

      implicit none 
      include 'netcdf.inc'

      character(len=*) :: varname 
      logical :: found 

      integer :: ncerr, varid
      integer, dimension(4) ::  start, count

      ncerr = nf_inq_varid(ncid, varname, varid)

      if ( ncerr .eq. 0 ) then 
        found = .TRUE.
      else
        found = .FALSE.
        print *, varname//nf_strerror(ncerr)
        stop
      endif
      
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
      count(1) = zlen 
      count(2) = xlen 
      count(3) = ylen 
      count(4) = tlen 

      ncerr = nf_get_vara_real(ncid, varid, start, count, trop_var ) 
c      print *, nf_strerror(ncerr)

      end subroutine trop_var_get


c--------------------------------------------------------------------
      subroutine trop_retrieve(varname,t_date,t_var)
c--------------------------------------------------------------------

      implicit none 

c... parameters
      character(len=*) :: varname 
      integer, dimension(2,52) :: t_date
      real, dimension(7,32,36,52) :: t_var

      call trop_var_get( varname )
      t_date = mtime
      t_var = trop_var

      end subroutine trop_retrieve

c--------------------------------------------------------------------
      subroutine trop_close()
c--------------------------------------------------------------------
c    close the netcdf file and free up memory allocated for lat, alt,
c    time, and constituent arrays 
c--------------------------------------------------------------------

      implicit none 
      include 'netcdf.inc'

      integer :: ncerr

      print *, 'trop_close()'
      ncerr = nf_close( ncid )

      deallocate( trop_tim )
      deallocate( trop_lev )
      deallocate( trop_lat )
      deallocate( trop_lon )
      deallocate( trop_var )

      end subroutine trop_close

c--------------------------------------------------------------------
      end module trop
c--------------------------------------------------------------------
