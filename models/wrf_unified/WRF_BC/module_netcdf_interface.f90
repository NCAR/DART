! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

MODULE module_netcdf_interface

  use        types_mod, only : r8
  use    utilities_mod, only : error_handler, E_ERR

  use netcdf

  implicit none
  private

public :: get_dims_cdf,      &
        get_gl_att_real_cdf, &
        put_gl_att_real_cdf, &
        get_var_3d_real_cdf, &
        get_var_2d_real_cdf, &
        put_var_3d_real_cdf, &
        put_var_2d_real_cdf, &
        get_times_cdf,       &
        put_time_cdf,        & 
        variable_exist

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

character(len=128) :: errstring

CONTAINS

!--------------------------------------------------------------------

subroutine get_dims_cdf( fname, var, idims, ndims, debug )

  implicit none

  character (len=80),     intent(in) :: fname
  character (len=*),      intent(in) :: var
  logical,                intent(in) :: debug
  integer, intent(out), dimension(4) :: idims
  integer,               intent(out) :: ndims

  integer  :: cdfid, id_data, i
  integer  :: dimids(4)

  if(debug) write(6,*) ' open netcdf file ', trim(fname)

  call check( nf90_open(fname, NF90_NOWRITE, cdfid) )

  call check( nf90_inq_varid(cdfid, var, id_data) )

  call check( nf90_Inquire_Variable(cdfid, id_data, ndims=ndims, dimids=dimids) )

  do i=1,ndims
     call check( nf90_inquire_dimension(cdfid, dimids(i), len=idims(i)) )
  enddo

  call check( nf90_close(cdfid) )

contains

  ! Internal subroutine - checks error status after each netcdf, prints 
  !                       text message each time an error code is returned. 
  subroutine check(istatus)
    integer, intent (in) :: istatus

    if(istatus /= nf90_noerr) call error_handler(E_ERR, 'get_dims_cdf', &
       trim(nf90_strerror(istatus)), source, revision, revdate)

  end subroutine check

end subroutine get_dims_cdf

!-------------------------------------------------------------------------------

subroutine get_gl_att_real_cdf( fname, att_name, value, debug )

  implicit none

  character (len=80), intent(in) :: fname
  character (len=*),  intent(in) :: att_name
  logical,            intent(in) :: debug
  real(r8),          intent(out) :: value

  integer :: cdfid

  if(debug) write(6,*) ' open netcdf file ', trim(fname)

  call check( nf90_open(fname, NF90_NOWRITE, cdfid) )

  call check( nf90_get_att(cdfid, nf90_global, att_name, value) )

  if(debug) write(6,*) ' global attribute ',att_name,' is ',value

  call check( nf90_close(cdfid) )

contains

  ! Internal subroutine - checks error status after each netcdf, prints 
  !                       text message each time an error code is returned. 
  subroutine check(istatus)
    integer, intent (in) :: istatus

    if(istatus /= nf90_noerr) call error_handler(E_ERR, 'get_gl_att_real_cdf', &
       trim(nf90_strerror(istatus)), source, revision, revdate)

  end subroutine check

end subroutine get_gl_att_real_cdf

!-------------------------------------------------------------------------------

subroutine put_gl_att_real_cdf( fname, att_name, value, debug )

  implicit none

  character (len=80), intent(in) :: fname
  character (len=*),  intent(in) :: att_name
  logical,            intent(in) :: debug
  real(r8),           intent(in) :: value

  integer :: cdfid

  if(debug) write(6,*) ' open netcdf file ', trim(fname)

  call check( nf90_open(fname, NF90_WRITE, cdfid) )

  call check( nf90_redef(cdfid) )
  call check( nf90_put_att(cdfid, nf90_global, att_name, value) )

  if(debug) write(6,*) ' global attribute ',att_name,' is ',value

  call check( nf90_close(cdfid) )

contains

  ! Internal subroutine - checks error status after each netcdf, prints 
  !                       text message each time an error code is returned. 
  subroutine check(istatus)
    integer, intent (in) :: istatus

    if(istatus /= nf90_noerr) call error_handler(E_ERR, 'put_gl_att_real_cdf', &
       trim(nf90_strerror(istatus)), source, revision, revdate)

  end subroutine check


end subroutine put_gl_att_real_cdf


!--------------------------------------------------------------------

subroutine get_var_3d_real_cdf( fname, var, data, &
     i1, i2, i3, time, debug )

  implicit none

  integer,                    intent(in) :: i1, i2, i3, time
  character (len=80),         intent(in) :: fname
  logical,                    intent(in) :: debug
  character (len=*),          intent(in) :: var
  real(r8), dimension(i1,i2,i3), intent(out) :: data

  integer            :: cdfid, id_data
  character (len=80) :: varnam
  integer            :: ndims, idims(4), dimids(4)
  integer            :: i, ivtype

  if(debug) write(6,*) ' open netcdf file ', trim(fname)

  call check( nf90_open(fname, NF90_NOWRITE, cdfid) )

  call check( nf90_inq_varid(cdfid, var, id_data) )

  call check( nf90_Inquire_Variable(cdfid, id_data, name=varnam, xtype=ivtype, ndims=ndims, dimids=dimids) )

  if(debug) then
     write(6,*) ' number of dims for ',var,' ',ndims
     write(unit=*, fmt='(a,i6)') ' ivtype=', ivtype
     write(unit=*, fmt='(a, a)') ' varnam=', trim(varnam)
     write(unit=*, fmt='(a,i6)') ' kind(data)=', kind(data)
  endif

  do i=1,ndims
     call check( nf90_inquire_dimension(cdfid, dimids(i), len=idims(i)) )
     if(debug) write(6,*) ' dimension ',i,idims(i)
  enddo

!  check the dimensions

  if(  (i1 /= idims(1)) .or.  &
       (i2 /= idims(2)) .or.  &
       (i3 /= idims(3)) .or.  &
       (time > idims(4))     )  then

     write(6,*) ' error in get_var_3d_real read, dimension problem '
     write(6,*) i1, idims(1)
     write(6,*) i2, idims(2)
     write(6,*) i3, idims(3)
     write(6,*) time, idims(4)
     stop

  end if

!  get the data
  
  call check( nf90_get_var(cdfid, id_data, data, start = (/ 1, 1, 1, time /)) )

  if(debug) write(unit=*, fmt='(a,e24.12)') ' Sample data=', data(1,1,1)

  call check( nf90_close(cdfid) )

contains

  ! Internal subroutine - checks error status after each netcdf, prints 
  !                       text message each time an error code is returned. 
  subroutine check(istatus)
    integer, intent (in) :: istatus

    if(istatus /= nf90_noerr) call error_handler(E_ERR, 'get_var_3d_real_cdf', &
       trim(nf90_strerror(istatus)), source, revision, revdate)

  end subroutine check

end subroutine get_var_3d_real_cdf

!--------------------------------------------------------------------

subroutine get_var_2d_real_cdf( fname, var, data, &
     i1, i2, time, debug )

  implicit none

  integer,                 intent(in) :: i1, i2, time
  character (len=80),      intent(in) :: fname
  logical,                 intent(in) :: debug
  character (len=*),       intent(in) :: var
  real(r8), dimension(i1,i2), intent(out) :: data

  integer            :: cdfid, id_data
  character (len=80) :: varnam
  integer            :: ndims, idims(4), dimids(4)
  integer            :: i, ivtype

  if(debug) write(6,*) ' open netcdf file ', trim(fname)

  call check( nf90_open(fname, NF90_NOWRITE, cdfid) )

  call check( nf90_inq_varid(cdfid, var, id_data) )

  call check( nf90_Inquire_Variable(cdfid, id_data, name=varnam, xtype=ivtype, ndims=ndims, dimids=dimids) )

  if(debug) then
     write(6,*) ' number of dims for ',var,' ',ndims
     write(unit=*, fmt='(a,i6)') ' ivtype=', ivtype
     write(unit=*, fmt='(a, a)') ' varnam=', trim(varnam)
  endif

  do i=1,ndims
     call check( nf90_inquire_dimension(cdfid, dimids(i), len=idims(i)) )
     if(debug) write(6,*) ' dimension ',i,idims(i)
  enddo

!  check the dimensions

  if(  (i1 /= idims(1)) .or.  &
       (i2 /= idims(2)) .or.  &
       (time > idims(3))     )  then

     write(6,*) ' error in get_var_2d_real read, dimension problem '
     write(6,*) i1, idims(1)
     write(6,*) i2, idims(2)
     write(6,*) time, idims(4)
     stop

  end if

!  get the data

  call check( nf90_get_var(cdfid, id_data, data, start = (/ 1, 1, time /)) )

  if(debug) write(unit=*, fmt='(a,e24.12)') ' Sample data=', data(1,1)

  call check( nf90_close(cdfid) )

contains

  ! Internal subroutine - checks error status after each netcdf, prints 
  !                       text message each time an error code is returned. 
  subroutine check(istatus)
    integer, intent (in) :: istatus

    if(istatus /= nf90_noerr) call error_handler(E_ERR, 'get_var_2d_real_cdf', &
       trim(nf90_strerror(istatus)), source, revision, revdate)

  end subroutine check

end subroutine get_var_2d_real_cdf

!--------------------------------------------------------------------

subroutine put_var_3d_real_cdf( fname, var, data, &
     i1, i2, i3, time, debug )

  implicit none

  integer,                       intent(in) :: i1, i2, i3, time
  character (len=80),            intent(in) :: fname
  logical,                       intent(in) :: debug
  character (len=*),             intent(in) :: var
  real(r8), dimension(i1,i2,i3), intent(in) :: data

  integer            :: cdfid, id_data
  character (len=80) :: varnam
  integer            :: ndims, idims(4), dimids(4)
  integer            :: i, ivtype

  if(debug) write(6,*) ' open netcdf file ', trim(fname)

  call check( nf90_open(fname, NF90_WRITE, cdfid) )

  call check( nf90_inq_varid(cdfid, var, id_data) )

  call check( nf90_Inquire_Variable(cdfid, id_data, name=varnam, xtype=ivtype, ndims=ndims, dimids=dimids) )

  if(debug) write(6,*) ' number of dims for ',var,' ',ndims

  do i=1,ndims
     call check( nf90_inquire_dimension(cdfid, dimids(i), len=idims(i)) )
     if(debug) write(6,*) ' dimension ',i,idims(i)
  enddo

!  check the dimensions

  if(  (i1 /= idims(1)) .or.  &
       (i2 /= idims(2)) .or.  &
       (i3 /= idims(3)) .or.  &
       (time > idims(4))     )  then

     write(6,*) ' error in put_var_3d_real read, dimension problem '
     write(6,*) i1, idims(1)
     write(6,*) i2, idims(2)
     write(6,*) i3, idims(3)
     write(6,*) time, idims(4)
     stop

  end if

!  write the data
  
  call check( nf90_put_var(cdfid, id_data, data, start = (/ 1, 1, 1, time /)) )

  call check( nf90_close(cdfid) )

contains

  ! Internal subroutine - checks error status after each netcdf, prints 
  !                       text message each time an error code is returned. 
  subroutine check(istatus)
    integer, intent (in) :: istatus

    if(istatus /= nf90_noerr) call error_handler(E_ERR, 'put_var_3d_real_cdf', &
       trim(nf90_strerror(istatus)), source, revision, revdate)

  end subroutine check

end subroutine put_var_3d_real_cdf

!--------------------------------------------------------------------

subroutine put_var_2d_real_cdf( fname, var, data, &
     i1, i2, time, debug )

  implicit none

  integer,                    intent(in) :: i1, i2, time
  character (len=80),         intent(in) :: fname
  logical,                    intent(in) :: debug
  character (len=*),          intent(in) :: var
  real(r8), dimension(i1,i2), intent(in) :: data

  integer            :: cdfid, id_data
  character (len=80) :: varnam
  integer            :: ndims, idims(3), dimids(3)
  integer            :: i, ivtype

  if(debug) write(6,*) ' open netcdf file ', trim(fname)

  call check( nf90_open(fname, NF90_WRITE, cdfid) )

  call check( nf90_inq_varid(cdfid, var, id_data) )

  call check( nf90_Inquire_Variable(cdfid, id_data, name=varnam, xtype=ivtype, ndims=ndims, dimids=dimids) )

  if(debug) write(6,*) ' number of dims for ',var,' ',ndims

  do i=1,ndims
     call check( nf90_inquire_dimension(cdfid, dimids(i), len=idims(i)) )
     if(debug) write(6,*) ' dimension ',i,idims(i)
  enddo

!  check the dimensions

  if(  (i1 /= idims(1)) .or.  &
       (i2 /= idims(2)) .or.  &
       (time > idims(3))     )  then

     write(6,*) ' error in put_var_2d_real read, dimension problem '
     write(6,*) i1, idims(1)
     write(6,*) i2, idims(2)
     write(6,*) time, idims(3)
     stop

  end if

!  write the data
  
  call check( nf90_put_var(cdfid, id_data, data, start = (/ 1, 1, time /)) )

  call check( nf90_close(cdfid) )

contains

  ! Internal subroutine - checks error status after each netcdf, prints 
  !                       text message each time an error code is returned. 
  subroutine check(istatus)
    integer, intent (in) :: istatus

    if(istatus /= nf90_noerr) call error_handler(E_ERR, 'put_var_2d_real_cdf', &
       trim(nf90_strerror(istatus)), source, revision, revdate)

  end subroutine check

end subroutine put_var_2d_real_cdf

!--------------------------------------------------------------------

subroutine get_times_cdf( fname, time_name, times, n_times, max_times, debug )

  implicit none

  integer,            intent(in)  :: max_times
  integer,            intent(out) :: n_times
  character (len=80), intent(in)  :: fname, time_name
  character (len=19), intent(out) :: times(max_times)
  logical,            intent(in)  :: debug

  integer            :: cdfid, id_time
  character (len=80) :: varnam, time1
  integer            :: ndims, idims(max_times)
  integer            :: istart(max_times),iend(max_times), dimids(max_times)
  integer            :: i, ivtype

  if(debug) write(6,*) ' open netcdf file ', trim(fname)

  call check( nf90_open(fname, NF90_NOWRITE, cdfid) )

  call check( nf90_inq_varid(cdfid, time_name, id_time) )

  call check( nf90_Inquire_Variable(cdfid, id_time, name=varnam, xtype=ivtype, ndims=ndims, dimids=dimids) )

  do i=1,ndims
     call check( nf90_inquire_dimension(cdfid, dimids(i), len=idims(i)) )
     if(debug) write(6,*) ' dimension ',i,idims(i)
  enddo

!  get the times

  n_times = idims(2)
  if (n_times > max_times) then
    write(errstring, '(2(A,I6))') 'number of times in file ', n_times, &
                                  '  is larger than allocated space ', max_times
    call error_handler(E_ERR, 'get_times_cdf', errstring, source, revision, revdate, &
                       text2='increase max_times in [pert,update]_wrf_bc.f90 and recompile')
  endif

  do i=1,idims(2)
    istart(1) = 1
    iend(1) = idims(1)
    istart(2) = i
    iend(2) = 1

    call check( nf90_get_var(cdfid, id_time, times(i), start = (/ 1, i /), &
                             count = (/idims(1), 1/) ) )
!    rcode = NF_GET_VARA_TEXT  ( cdfid, id_time,  &
!                                istart, iend,    &
!                                times(i)          )
    time1 = times(i)

    if(debug) write(6,*) trim(fname), time1(1:19)
  enddo

  call check( nf90_close(cdfid) )
contains

  ! Internal subroutine - checks error status after each netcdf, prints 
  !                       text message each time an error code is returned. 
  subroutine check(istatus)
    integer, intent (in) :: istatus

    if(istatus /= nf90_noerr) call error_handler(E_ERR, 'get_times_cdf', &
       trim(nf90_strerror(istatus)), source, revision, revdate)

  end subroutine check

end subroutine get_times_cdf

!--------------------------------------------------------------------

subroutine put_time_cdf( fname, time_name, char, itime, debug )

  implicit none

  integer,            intent(in)  :: itime
  character (len=80), intent(in)  :: fname, time_name
  character (len=19), intent(in)  :: char
  logical,            intent(in)  :: debug

  integer            :: cdfid, id_time

  if(debug) write(6,*) ' open netcdf file ', trim(fname)

  call check( nf90_open(fname, NF90_WRITE, cdfid) )

  call check( nf90_inq_varid(cdfid, time_name, id_time) )

  call check( nf90_put_var(cdfid, id_time, char, start = (/ 1, itime /)) )

  call check( nf90_close(cdfid) )

contains

  ! Internal subroutine - checks error status after each netcdf, prints 
  !                       text message each time an error code is returned. 
  subroutine check(istatus)
    integer, intent (in) :: istatus

    if(istatus /= nf90_noerr) call error_handler(E_ERR, 'put_time_cdf', &
       trim(nf90_strerror(istatus)), source, revision, revdate)

  end subroutine check

end subroutine put_time_cdf

function variable_exist( fname, var_name )

  implicit none

  character(len=*), intent(in) :: fname, var_name

  integer :: cdfid, var_id, rcode
  logical :: variable_exist

  variable_exist = .false.
  call check( nf90_open(fname, nf90_nowrite, cdfid) )
  rcode = nf90_inq_varid(cdfid, var_name, var_id)
  call check( nf90_close(cdfid) )
  if ( rcode == nf90_noerr )  variable_exist = .true.
  return 
 
contains

  ! Internal subroutine - checks error status after each netcdf, prints
  !                       text message each time an error code is returned.
  subroutine check(istatus)
    integer, intent (in) :: istatus

    if(istatus /= nf90_noerr) call error_handler(E_ERR, 'variable_exist', &
       trim(nf90_strerror(istatus)), source, revision, revdate)

  end subroutine check

end function variable_exist

END MODULE module_netcdf_interface

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
