! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
module module_netcdf_interface

!
! <next four lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
!

  use        types_mod, only : r8
  use    utilities_mod, only : error_handler, E_ERR
  use netcdf

  implicit none
  private

  public get_dims_cdf,          &
         get_gl_att_real_cdf,   &
         get_var_3d_real_cdf,   &
         get_var_2d_real_cdf,   &
         put_var_3d_real_cdf


! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

CONTAINS

!--------------------------------------------------------------------

  subroutine get_dims_cdf( file, var, dims, ndims, debug )
        
    implicit none

    include 'netcdf.inc'

    character (len=80),     intent(in) :: file
    character (len=*),      intent(in) :: var
    logical,                intent(in) :: debug
    integer, dimension(4), intent(out) :: dims
    integer,               intent(out) :: ndims

    integer cdfid, rcode, id_time
    character (len=80) :: varnam
    integer :: natts, dimids(10)
    integer :: i, ivtype

    cdfid = ncopn(file, NCNOWRIT, rcode )
  
    if( rcode == 0) then
       if(debug) write(6,*) ' open netcdf file ', trim(file)
    else
       write(6,*) ' error opening netcdf file ', trim(file)
       stop
    end if

    id_time = ncvid( cdfid, var, rcode )

    rcode = nf_inq_var( cdfid, id_time, varnam, ivtype, ndims, dimids, natts )
    if(debug) then
       write(6,*) ' number of dims for ',var,' ',ndims
    endif
    do i=1,ndims
       rcode = nf_inq_dimlen( cdfid, dimids(i), dims(i) )
       if(debug) write(6,*) ' dimension ',i,dims(i)
    enddo

    call ncclos(cdfid,rcode)

  end subroutine get_dims_cdf

!-------------------------------------------------------------------------------

  subroutine get_gl_att_real_cdf( file, att_name, value, debug )
        
    implicit none

    include 'netcdf.inc'

    character (len=80), intent(in) :: file
    character (len=*),  intent(in) :: att_name
    logical,            intent(in) :: debug
    real(r8),          intent(out) :: value

    integer :: cdfid, rcode

    cdfid = ncopn(file, NCNOWRIT, rcode )

    if( rcode == 0) then
       if(debug) write(6,*) ' open netcdf file ', trim(file)
    else
       write(6,*) ' error opening netcdf file ', trim(file)
       stop
    end if

    call check( nf90_get_att(cdfid, nf90_global, att_name, value) )

    call ncclos(cdfid,rcode)

    if(debug) write(6,*) ' global attribute ',att_name,' is ',value

contains

  ! Internal subroutine - checks error status after each netcdf, prints 
  !                       text message each time an error code is returned. 
  subroutine check(istatus)
    integer, intent ( in) :: istatus

    if(istatus /= nf90_noerr) call error_handler(E_ERR, 'get_gl_att_real_cdf', &
       trim(nf90_strerror(istatus)), source, revision, revdate)

  end subroutine check

  end subroutine get_gl_att_real_cdf

!--------------------------------------------------------------------

  subroutine get_var_3d_real_cdf( file, var, data, &
                                  i1, i2, i3, time, debug )
        
    implicit none

    include 'netcdf.inc'

    integer,                        intent(in) :: i1, i2, i3, time
    character (len=80),             intent(in) :: file
    logical,                        intent(in) :: debug
    character (len=*),              intent(in) :: var
    real(r8), dimension(i1,i2,i3), intent(out) :: data

    integer cdfid, rcode, id_data
    character (len=80) :: varnam
    integer :: ndims, natts, idims(10), dimids(10)
    integer :: i, ivtype

    cdfid = ncopn(file, NCNOWRIT, rcode )

    if( rcode == 0) then
       if(debug) write(6,*) ' open netcdf file ', trim(file)
    else
       write(6,*) ' error openiing netcdf file ', trim(file)
       stop
    end if

    id_data = ncvid( cdfid, var, rcode )

    rcode = nf_inq_var( cdfid, id_data, varnam, ivtype, ndims, dimids, natts )

    if(debug) then
       write(6,*) ' number of dims for ',var,' ',ndims
       write(unit=*, fmt='(a,i6)') ' ivtype=', ivtype
       write(unit=*, fmt='(a, a)') ' varnam=', trim(varnam)
       write(unit=*, fmt='(a,i6)') ' kind(data)=', kind(data)
    endif

    do i=1,ndims
       rcode = nf_inq_dimlen( cdfid, dimids(i), idims(i) )
       if(debug) write(6,*) ' dimension ',i,idims(i)
    enddo

    !  check the dimensions

    if( (i1 /= idims(1)) .or.  &
         (i2 /= idims(2)) .or.  &
         (i3 /= idims(3)) .or.  &
         (time > idims(4))     )  then

       write(6,*) ' error in 3d_var_real read, dimension problem '
       write(6,*) i1, idims(1)
       write(6,*) i2, idims(2)
       write(6,*) i3, idims(3)
       write(6,*) time, idims(4)
       write(6,*) ' error stop 1'
       stop

    end if

    !  get the data

    call check( nf90_get_var(cdfid,id_data,data) )

    if(debug) then
       write(unit=*, fmt='(a,e24.12)') ' Sample data=', data(1,1,1)
    endif

    call ncclos(cdfid,rcode)

contains

  ! Internal subroutine - checks error status after each netcdf, prints 
  !                       text message each time an error code is returned. 
  subroutine check(istatus)
    integer, intent ( in) :: istatus

    if(istatus /= nf90_noerr) call error_handler(E_ERR, 'get_var_3d_real_cdf', &
       trim(nf90_strerror(istatus)), source, revision, revdate)

  end subroutine check

  end subroutine get_var_3d_real_cdf

!--------------------------------------------------------------------

  subroutine get_var_2d_real_cdf( file, var, data, &
                                  i1, i2, time, debug )
        
    implicit none

    include 'netcdf.inc'

    integer,                     intent(in) :: i1, i2, time
    character (len=80),          intent(in) :: file
    logical,                     intent(in) :: debug
    character (len=*),           intent(in) :: var
    real(r8), dimension(i1,i2), intent(out) :: data

    integer cdfid, rcode, id_data
    character (len=80) :: varnam
    integer :: ndims, natts, idims(10), dimids(10)
    integer :: i, ivtype

    cdfid = ncopn(file, NCNOWRIT, rcode )

    if( rcode == 0) then
       if(debug) write(6,*) ' open netcdf file ', trim(file)
    else
       write(6,*) ' error opening netcdf file ', trim(file)
       stop
    end if

    id_data = ncvid( cdfid, var, rcode )

    rcode = nf_inq_var( cdfid, id_data, varnam, ivtype, ndims, dimids, natts )
    if(debug) then
       write(6,*) ' number of dims for ',var,' ',ndims
    endif
    do i=1,ndims
       rcode = nf_inq_dimlen( cdfid, dimids(i), idims(i) )
       if(debug) then
          write(6,*) ' dimension ',i,idims(i)
          write(unit=*, fmt='(a,i6)') ' ivtype=', ivtype
          write(unit=*, fmt='(a, a)') ' varnam= ', trim(varnam)
       endif
    enddo

    !  check the dimensions

    if( (i1 /= idims(1)) .or.  &
         (i2 /= idims(2)) .or.  &
         (time > idims(3))     )  then

       write(6,*) ' error in 2d_var_real read, dimension problem '
       write(6,*) i1, idims(1)
       write(6,*) i2, idims(2)
       write(6,*) time, idims(4)
       write(6,*) ' error stop 2'
       stop

    end if

    !  get the data

    call check( nf90_get_var(cdfid,id_data,data) )

    if(debug) then
       write(unit=*, fmt='(a,e24.12)') ' Sample data=', data(1,1)
    endif

    call ncclos(cdfid,rcode)

contains

  ! Internal subroutine - checks error status after each netcdf, prints 
  !                       text message each time an error code is returned. 
  subroutine check(istatus)
    integer, intent ( in) :: istatus

    if(istatus /= nf90_noerr) call error_handler(E_ERR, 'get_var_2d_real_cdf', &
       trim(nf90_strerror(istatus)), source, revision, revdate)

  end subroutine check

  end subroutine get_var_2d_real_cdf

!--------------------------------------------------------------------

  subroutine put_var_3d_real_cdf( file, var, data, &
                                  i1, i2, i3, time, debug )
        
    implicit none

    include 'netcdf.inc'

    integer,                       intent(in) :: i1, i2, i3, time
    character (len=80),            intent(in) :: file
    logical,                       intent(in) :: debug
    character (len=*),             intent(in) :: var
    real(r8), dimension(i1,i2,i3), intent(in) :: data

    integer cdfid, rcode, id_data
    character (len=80) :: varnam
    integer :: ndims, natts, idims(10), dimids(10)
    integer :: i, ivtype

    cdfid = ncopn(file, NCWRITE, rcode )

    if( rcode == 0) then
       if(debug) write(6,*) ' open netcdf file ', trim(file)
    else
       write(6,*) ' error openiing netcdf file ', trim(file)
       stop
    end if

    id_data = ncvid( cdfid, var, rcode )

    rcode = nf_inq_var( cdfid, id_data, varnam, ivtype, ndims, dimids, natts )
    if(debug) then
       write(6,*) ' number of dims for ',var,' ',ndims
    endif
    do i=1,ndims
       rcode = nf_inq_dimlen( cdfid, dimids(i), idims(i) )
       if(debug) write(6,*) ' dimension ',i,idims(i)
    enddo

    !  check the dimensions

    if(  (i1 /= idims(1)) .or.  &
         (i2 /= idims(2)) .or.  &
         (i3 /= idims(3)) .or.  &
         (time > idims(4))     )  then

       write(6,*) ' error in 3d_var_real read, dimension problem '
       write(6,*) i1, idims(1)
       write(6,*) i2, idims(2)
       write(6,*) i3, idims(3)
       write(6,*) time, idims(4)
       write(6,*) ' error stop 3'
       stop

    end if

    !  write the data

    call check( nf90_put_var(cdfid,id_data,data) )

    call ncclos(cdfid,rcode)

contains

  ! Internal subroutine - checks error status after each netcdf, prints 
  !                       text message each time an error code is returned. 
  subroutine check(istatus)
    integer, intent ( in) :: istatus

    if(istatus /= nf90_noerr) call error_handler(E_ERR, 'put_var_3d_real_cdf', &
       trim(nf90_strerror(istatus)), source, revision, revdate)

  end subroutine check

  end subroutine put_var_3d_real_cdf

!--------------------------------------------------------------------

END MODULE module_netcdf_interface
