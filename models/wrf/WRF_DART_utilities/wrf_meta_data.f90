! program to read in wrf metadata and static data for DART

!  $Source$
!  $Revision$
!  $Date$

implicit none

! let CVS fill strings ... DO NOT EDIT ...
character(len=128) :: &
   source   = "$Source$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

TYPE wrf_static_data_for_dart

   integer :: bt, sn, we

   real :: p_top, dx, dy, dt
   real, dimension(:), pointer :: dn, dnw
   integer :: n_moist
   real, dimension(:,:), pointer :: mub, latitude, longitude
   real, dimension(:,:,:), pointer :: phb

   integer :: model_size, number_of_wrf_variables
   integer, dimension(:), pointer :: var_type
   integer, dimension(:,:), pointer :: var_index
   integer, dimension(:,:), pointer :: var_size

end type

integer, parameter :: TYPE_MU = 0, TYPE_T = 1, TYPE_U = 2, TYPE_V = 3, &
                      TYPE_W = 4,  TYPE_GZ = 5,                        &
                      TYPE_QV = 6,  TYPE_QC = 7,  TYPE_QR = 8,         &
                      TYPE_QI = 9,  TYPE_QS = 10, TYPE_QG = 11

include 'netcdf.inc'

character (len = 80)      :: path
     integer              :: mode
     integer              :: ncid, bt_id, we_id, sn_id

integer :: status
character (len=80) :: name
logical :: debug
integer :: var_id, ind, i
integer, dimension(5) :: count, start, stride, map

type(wrf_static_data_for_dart) :: wrf

integer :: n_values

!----------

debug = .true.
mode = 0
status = nf_open('wrfinput', mode, ncid) 
if (status /= nf_noerr) call handle_err(1,status)
if(debug) write(6,*) ' ncid is ',ncid

! get wrf grid dimensions

status = nf_inq_dimid(ncid, "bottom_top", bt_id)
if (status /= nf_noerr) call handle_err(2,status)
status = nf_inq_dim(ncid, bt_id, name, wrf%bt)
if (status /= nf_noerr) call handle_err(3,status)

status = nf_inq_dimid(ncid, "south_north", sn_id)
if (status /= nf_noerr) call handle_err(2,status)
status = nf_inq_dim(ncid, sn_id, name, wrf%sn)
if (status /= nf_noerr) call handle_err(3,status)

status = nf_inq_dimid(ncid, "west_east", we_id)
if (status /= nf_noerr) call handle_err(2,status)
status = nf_inq_dim(ncid, we_id, name, wrf%we)
if (status /= nf_noerr) call handle_err(3,status)

write(6,*) ' dimensions bt, sn, we are ',wrf%bt,wrf%sn,wrf%we

! get meta data adn static data we need

count  = 1
start  = 1
stride = 1
map    = 1

call netcdf_read_write_var( "P_TOP",ncid, var_id, wrf%p_top,        &
                            start, count, stride, map, 'INPUT ', debug  )
if(debug) write(6,*) ' p_top is ',wrf%p_top

status = nf_get_att_real(ncid, nf_global, "DX", wrf%dx)
if (status /= nf_noerr) call handle_err(4,status)
if(debug) write(6,*) ' dx is ',wrf%dx

status = nf_get_att_real(ncid, nf_global, "DY", wrf%dy)
if (status /= nf_noerr) call handle_err(4,status)
if(debug) write(6,*) ' dy is ',wrf%dy

status = nf_get_att_real(ncid, nf_global, "DT", wrf%dt)
if (status /= nf_noerr) call handle_err(4,status)
if(debug) write(6,*) ' dt is ',wrf%dt

wrf%n_moist = 3  ! hardwired for now, need a way to find this

!  get 1D (z) static data defining grid levels

count(1)  = wrf%bt
allocate(wrf%dn(1:wrf%bt))
call netcdf_read_write_var( "DN",ncid, var_id, wrf%dn,        &
                            start, count, stride, map, 'INPUT ', debug  )
if(debug) write(6,*) ' dn ',wrf%dn

count(1)  = wrf%bt
allocate(wrf%dnw(1:wrf%bt))
call netcdf_read_write_var( "DNW",ncid, var_id, wrf%dnw,        &
                            start, count, stride, map, 'INPUT ', debug  )
if(debug) write(6,*) ' dnw is ',wrf%dnw

!  get 2D (x,y) base state for mu, latitude, longitude

count(1)  = wrf%we
count(2)  = wrf%sn
allocate(wrf%mub(1:wrf%we,1:wrf%sn))
call netcdf_read_write_var( "MUB",ncid, var_id, wrf%mub,        &
                            start, count, stride, map, 'INPUT ', debug  )
if(debug) then
    write(6,*) ' corners of mub '
    write(6,*) wrf%mub(1,1),wrf%mub(wrf%we,1),  &
               wrf%mub(1,wrf%sn),wrf%mub(wrf%we,wrf%sn)
end if

allocate(wrf%longitude(1:wrf%we,1:wrf%sn))
call netcdf_read_write_var( "XLONG",ncid, var_id, wrf%longitude,        &
                            start, count, stride, map, 'INPUT ', debug  )
allocate(wrf%latitude(1:wrf%we,1:wrf%sn))
call netcdf_read_write_var( "XLAT",ncid, var_id, wrf%latitude,        &
                            start, count, stride, map, 'INPUT ', debug  )

! get 3D base state geopotential

count(1)  = wrf%we
count(2)  = wrf%sn
count(3)  = wrf%bt+1
allocate(wrf%phb(1:wrf%we,1:wrf%sn,1:wrf%bt+1))
call netcdf_read_write_var( "PHB",ncid, var_id, wrf%phb,        &
                            start, count, stride, map, 'INPUT ', debug  )
if(debug) then
    write(6,*) ' corners of phb '
    write(6,*) wrf%phb(1,1,1),wrf%phb(wrf%we,1,1),  &
               wrf%phb(1,wrf%sn,1),wrf%phb(wrf%we,wrf%sn,1)
    write(6,*) wrf%phb(1,1,wrf%bt+1),wrf%phb(wrf%we,1,wrf%bt+1),  &
               wrf%phb(1,wrf%sn,wrf%bt+1),wrf%phb(wrf%we,wrf%sn,wrf%bt+1)
end if

! close data file, we have all we need

status = nf_close(ncid)
if (status /= nf_noerr) call handle_err(4,status)

!  build the map into the 1D DART vector for WRF data

  wrf%number_of_wrf_variables = 6 + wrf%n_moist
  allocate(wrf%var_type(12))  ! map for var_type, we'll always carry all
  wrf%var_type(1)  = TYPE_U   !  possible types
  wrf%var_type(2)  = TYPE_V
  wrf%var_type(3)  = TYPE_W
  wrf%var_type(4)  = TYPE_GZ
  wrf%var_type(5)  = TYPE_T
  wrf%var_type(6)  = TYPE_MU
  wrf%var_type(7)  = TYPE_QV
  wrf%var_type(8)  = TYPE_QC
  wrf%var_type(9)  = TYPE_QR
  wrf%var_type(10) = TYPE_QI
  wrf%var_type(11) = TYPE_QS
  wrf%var_type(12) = TYPE_QG

  allocate(wrf%var_index(2,6 + wrf%n_moist)) ! indices into 1D array
  allocate(wrf%var_size(3,6 + wrf%n_moist)) ! dimension of variables
  
  ind = 1                         ! *** u field ***
  wrf%var_size(1,ind) = wrf%we + 1  
  wrf%var_size(2,ind) = wrf%sn
  wrf%var_size(3,ind) = wrf%bt
  wrf%var_index(1,ind) = 1
  wrf%var_index(2,ind) = wrf%var_index(1,ind) - 1 +  &
                       wrf%var_size(1,ind)*wrf%var_size(2,ind)*wrf%var_size(3,ind)

  ind = ind + 1                   ! *** v field ***
  wrf%var_size(1,ind) = wrf%we
  wrf%var_size(2,ind) = wrf%sn + 1
  wrf%var_size(3,ind) = wrf%bt
  wrf%var_index(1,ind) = wrf%var_index(2,ind-1) + 1
  wrf%var_index(2,ind) = wrf%var_index(1,ind) - 1 +   &
                       wrf%var_size(1,ind)*wrf%var_size(2,ind)*wrf%var_size(3,ind)

  ind = ind + 1                   ! *** w field ***
  wrf%var_size(1,ind) = wrf%we
  wrf%var_size(2,ind) = wrf%sn
  wrf%var_size(3,ind) = wrf%bt + 1
  wrf%var_index(1,ind) = wrf%var_index(2,ind-1) + 1
  wrf%var_index(2,ind) = wrf%var_index(1,ind) - 1 +    &
                       wrf%var_size(1,ind)*wrf%var_size(2,ind)*wrf%var_size(3,ind)

  ind = ind + 1                   ! *** geopotential field ***
  wrf%var_size(1,ind) = wrf%we
  wrf%var_size(2,ind) = wrf%sn
  wrf%var_size(3,ind) = wrf%bt + 1
  wrf%var_index(1,ind) = wrf%var_index(2,ind-1) + 1
  wrf%var_index(2,ind) = wrf%var_index(1,ind) - 1 +   &
                       wrf%var_size(1,ind)*wrf%var_size(2,ind)*wrf%var_size(3,ind)

  ind = ind + 1                   ! *** theta field ***
  wrf%var_size(1,ind) = wrf%we
  wrf%var_size(2,ind) = wrf%sn
  wrf%var_size(3,ind) = wrf%bt
  wrf%var_index(1,ind) = wrf%var_index(2,ind-1) + 1
  wrf%var_index(2,ind) = wrf%var_index(1,ind) - 1 +   &
                       wrf%var_size(1,ind)*wrf%var_size(2,ind)*wrf%var_size(3,ind)

  ind = ind + 1                   ! *** mu field ***
  wrf%var_size(1,ind) = wrf%we
  wrf%var_size(2,ind) = wrf%sn
  wrf%var_size(3,ind) = 1
  wrf%var_index(1,ind) = wrf%var_index(2,ind-1) + 1
  wrf%var_index(2,ind) = wrf%var_index(1,ind) - 1 +   &
                       wrf%var_size(1,ind)*wrf%var_size(2,ind)*wrf%var_size(3,ind)

  do i = 1, wrf%n_moist
    ind = ind + 1                   ! *** moisture field ***
    wrf%var_size(1,ind) = wrf%we
    wrf%var_size(2,ind) = wrf%sn
    wrf%var_size(3,ind) = wrf%bt
    wrf%var_index(1,ind) = wrf%var_index(2,ind-1) + 1
    wrf%var_index(2,ind) = wrf%var_index(1,ind) - 1 +   &
                         wrf%var_size(1,ind)*wrf%var_size(2,ind)*wrf%var_size(3,ind)
  enddo

  wrf%model_size = wrf%var_index(2,wrf%number_of_wrf_variables)

! we're finished here

end

!****************************1***********************************

subroutine handle_err(ifn,status)
implicit none
integer :: ifn, status
character (len=80) :: nf_strerror

write(6,*) ' error for netcdf function ',ifn
write(6,*) ' status code = ',status
write(6,'(a80)') nf_strerror(status)

stop

end subroutine handle_err

!**********************************************************************

subroutine netcdf_read_write_var( variable, ncid, var_id, var,          &
                                  start, count, stride, map, in_or_out, debug )

implicit none
include 'netcdf.inc'
integer :: ncid, var_id
real, dimension(*) :: var
character (len=6) :: in_or_out
integer, dimension(5) :: start, count, stride, map
integer :: status
character (len=*) :: variable
logical :: debug


write(6,*) ' var for io is ',variable
status = nf_inq_varid(ncid, variable, var_id)
if(debug) write(6,*) variable, ' id = ',var_id
if (status /= nf_noerr) call handle_err(4,status)

if( in_or_out(1:5) == "INPUT" ) then
  write(6,*) ' call netcdf read ', ncid, var_id
  status = nf_get_var_real(ncid, var_id, var, start, count, stride, map)
  write(6,*) ' returned netcdf read '
else if( in_or_out(1:6) == "OUTPUT" ) then
  status = nf_put_var_real(ncid, var_id, var, start, count, stride, map)
else
  write(6,*) ' unknown IO function for var_id ',var_id, in_or_out
  stop
end if
if (status /= nf_noerr) call handle_err(100,status)

end subroutine netcdf_read_write_var

!**********************************************************************
