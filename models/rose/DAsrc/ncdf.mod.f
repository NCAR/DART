!----------------------------------------------------------------------
      module ncdf
!----------------------------------------------------------------------

      use params

      implicit none

      private
      public :: nc_init, nc_update, check_err, nc_add_global

      integer, parameter :: mtime_len = 3

      contains

!----------------------------------------------------------------------
      subroutine nc_init( nc_name, var_names, var_units, var_ids, 
     $                    var_cnt, ncid, nstp_id, time_id, mtime_id )
!----------------------------------------------------------------------

      implicit none

      include 'netcdf.inc'

      character(len=*), intent(in) :: nc_name

      integer, intent(in) :: var_cnt       ! number of variables 
      character (LEN=*), intent(in), dimension(var_cnt) :: var_names
      character (LEN=*), intent(in), dimension(var_cnt) :: var_units

      integer, intent(out) :: ncid         ! id of new file
      integer, intent(out) :: nstp_id      ! variable ids
      integer, intent(out) :: time_id
      integer, intent(out) :: mtime_id
      integer, intent(out), dimension(var_cnt) ::  var_ids
      
!... local variables
      integer  iret                        ! error status return

      integer  lon_dim                     ! dimension ids
      integer  lat_dim
      integer  lev_dim
      integer  time_dim
      integer  mtime_dim
      
      integer  lon_id                      ! variable ids
      integer  lat_id
      integer  z_id
      
      integer, parameter :: mtime_rank = 2 ! number of dimensions for
      integer, parameter :: var_rank = 4   ! each variable
      
      integer  mtime_dims(mtime_rank)      ! variable shapes
      integer  dims(var_rank)
      
      integer :: i
      real, dimension(nx) :: lon
      real, dimension(ny) :: lat
      real, dimension(nz) :: z

!... enter define mode

      print *, 'creating ncdf: ', nc_name
      iret = nf_create( nc_name, NF_CLOBBER, ncid)
      call check_err(iret)
      
!... define dimensions

      iret = nf_def_dim(ncid, 'lon', nx, lon_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'lat', ny, lat_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'lev', nz, lev_dim)
      call check_err(iret)

      iret = nf_def_dim(ncid, 'time', NF_UNLIMITED, time_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'mtime', mtime_len, mtime_dim)
      call check_err(iret)
      
!... define variables

      iret = nf_def_var(ncid, 'nstep', NF_INT, 1, time_dim, nstp_id)
      call check_err(iret)

      iret = nf_def_var(ncid, 'lon', NF_REAL, 1, lon_dim, lon_id)
      call check_err(iret)
      
      iret = nf_def_var(ncid, 'lat', NF_REAL, 1, lat_dim, lat_id)
      call check_err(iret)

      iret = nf_def_var(ncid, 'z', NF_REAL, 1, lev_dim, z_id)
      call check_err(iret)

      iret = nf_def_var(ncid, 'time', NF_REAL, 1, time_dim, time_id)
      call check_err(iret)

      mtime_dims(2) = time_dim
      mtime_dims(1) = mtime_dim

      iret = nf_def_var(ncid, 'mtime', NF_INT, mtime_rank, mtime_dims, 
     $                  mtime_id)
      call check_err(iret)
      
      dims(1) = lev_dim
      dims(2) = lon_dim
      dims(3) = lat_dim
      dims(4) = time_dim

      do i=1, var_cnt
        print *, trim(var_names(i))
        iret = nf_def_var(ncid, trim(var_names(i)), NF_REAL,
     $           var_rank, dims, var_ids(i))

        iret = nf_put_att_text(ncid,  var_ids(i), 'units', 
     $	        len(trim(var_units(i))), trim(var_units(i)) )

        call check_err(iret)
      enddo

!... assign attributes

      iret = nf_put_att_text(ncid, nstp_id, 'long_name', 15,
     $                       'model time step')
      iret = nf_put_att_text(ncid, lon_id, 'long_name', 35,
     $                       'geographic longitude (-west, +east)')
      iret = nf_put_att_text(ncid, lat_id, 'long_name', 36,
     $                       'geographic latitude (-south, +north)')
      iret = nf_put_att_text(ncid, z_id, 'long_name', 20,
     $                       'vertical height (km)')
      iret = nf_put_att_text(ncid, time_id, 'long_name', 28,
     $                       'elapsed model time (minutes)')
      iret = nf_put_att_text(ncid, mtime_id, 'long_name', 30,
     $                       'model times (year, doy, utsec)')

      iret = nf_put_att_text(ncid, lon_id, 'units', 7, 'degrees')
      iret = nf_put_att_text(ncid, lat_id, 'units', 7, 'degrees')
      iret = nf_put_att_text(ncid, z_id,   'units', 2, 'km')
      iret = nf_put_att_text(ncid, time_id,'units', 7, 'minutes')

!... leave define mode

      iret = nf_enddef(ncid)
      call check_err(iret)
      
!... store lon, lat, z

      do i=0,nx-1
        lon(i+1) = i*11.25
      enddo
      do i=0,ny-1
        lat(i+1) = i*5-87.5
      enddo
      do i=0,nz-1
        z(i+1) = i*2.5+17.5
      enddo

      iret = nf_put_var_real(ncid, lon_id, lon)
      iret = nf_put_var_real(ncid, lat_id, lat)
      iret = nf_put_var_real(ncid, z_id, z)

      iret = nf_sync(ncid) 
      call check_err(iret)

      end subroutine nc_init 

!----------------------------------------------------------------------
      subroutine nc_add_global( ncid )
!----------------------------------------------------------------------
!... sets global attributes for the netcdf output files
!----------------------------------------------------------------------

      use params, only : year0, day0, ut0, ntime, nseg
      
      implicit none

      include 'netcdf.inc'

      integer, intent(in) :: ncid          ! id of new file

!... local variables

      integer :: iret                      ! error status return
      integer :: intval(1)                 ! global attribute value

!... enter define mode

      iret = nf_redef(ncid)
      call check_err(iret)

!... add global variables

      iret = nf_put_att_text( ncid, NF_GLOBAL, 'model_version', 8, 
     $                       'rose 2.0')

      intval(1) = year0 
      iret = nf_put_att_int(ncid, NCGLOBAL, 'year0', NF_INT, 1, intval)
      
      intval(1) = day0 
      iret = nf_put_att_int(ncid, NCGLOBAL, 'day0',  NF_INT, 1, intval)
      
      intval(1) = ut0 
      iret = nf_put_att_int(ncid, NCGLOBAL, 'ut0',   NF_INT, 1, intval)
      
      intval(1) = ntime 
      iret = nf_put_att_int(ncid, NCGLOBAL, 'ntime', NF_INT, 1, intval)

      intval(1) = nseg           ! segment number for run (0=initial run)
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nseg',  NF_INT, 1, intval)

!... leave define mode

      iret = nf_enddef(ncid)
      call check_err(iret)

      end subroutine nc_add_global

!----------------------------------------------------------------------
      subroutine nc_update( ncid, frame, nstp_id, mtime_id, time_id,    
     $                      mtime, var_cnt, var_ids, vars )
!----------------------------------------------------------------------
!... updates netcdf output files with current model variables
!----------------------------------------------------------------------
       
      use dynam, only : nstep, ntime

      implicit none

      include 'netcdf.inc'

      integer, intent(in) :: ncid       ! id of netcdf file
      integer, intent(in) :: frame      ! unlimeted dimension index for output
      integer, intent(in) :: nstp_id    ! model step variable id
      integer, intent(in) :: time_id    ! elapsed time id
      integer, intent(in) :: mtime_id   ! id of model times (year, doy, utsec)
      integer, intent(in) :: var_cnt    ! number of variables to be written
      integer, intent(in) :: mtime(mtime_len)

      integer, dimension(var_cnt), intent(in) ::  var_ids
      real, dimension(nz, nx, ny, var_cnt), intent(in) :: vars 


!... local variables 

      integer :: iret                  ! error status return
      integer :: i
      real    :: t_elap
       
      integer, dimension(mtime_len,1) :: mtime_dummy

!... rank (number of dimensions) for each variable

      integer, parameter :: mtime_rank = 2
      integer, parameter :: var_rank = 4

!... starts and counts for array sections of record variables

      integer :: mtime_start(mtime_rank) = (/  1, 0 /)
      integer :: mtime_count(mtime_rank) = (/  mtime_len, 1 /)
      integer :: start(var_rank)         = (/  1,  1,  1, 0 /)
      integer :: count(var_rank)         = (/ nz, nx, ny, 1 /)
       
      print *, 'nc_update, frame:', frame
     
!... store nstep time and mtime

      iret = nf_put_vara_int(ncid, nstp_id, frame, 1, nstep)
      call check_err(iret)

      t_elap = 60.0 * real(nstep) / real(ntime)
      iret = nf_put_vara_real(ncid, time_id, frame, 1, t_elap )
      call check_err(iret)

      mtime_start(2) = frame
      mtime_dummy(1,1) = mtime(1) 
      mtime_dummy(2,1) = mtime(2)
      mtime_dummy(3,1) = mtime(3)
      iret = nf_put_vara_int(ncid, mtime_id, mtime_start,
     $                       mtime_count, mtime_dummy )
      call check_err(iret)

!... save dynamical fields

      start(4) = frame

      do i=1, var_cnt
        iret = nf_put_vara_real(ncid, var_ids(i), start, count,
     1                          vars(:,:,:,i) ) 
        call check_err(iret)
      enddo

c... force output queue to be flushed and file to be updated
      iret = nf_sync(ncid) 
      call check_err(iret)

      end subroutine nc_update 

c----------------------------------------------------------------------
      subroutine check_err(iret)
c----------------------------------------------------------------------

      implicit none

      include 'netcdf.inc'

      integer iret

      if (iret .ne. NF_NOERR) then
        print *, nf_strerror(iret)
        stop
      endif

      end subroutine check_err

      end module ncdf
