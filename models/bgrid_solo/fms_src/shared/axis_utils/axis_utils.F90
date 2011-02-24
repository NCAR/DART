!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module axis_utils_mod
  !
  !<CONTACT EMAIL="mjh@gfdl.noaa.gov">M.J. Harrison</CONTACT>
  !
  !<REVIEWER EMAIL="bw@gfdl.noaa.gov">Bruce Wyman</REVIEWER>
  !

  !</OVERVIEW>
  ! A set of utilities for manipulating axes and extracting axis
  ! attributes
  !</OVERVIEW>

  !<DESCRIPTION>
  !
  ! subroutine get_axis_cart(axis,cart) : Returns X,Y,Z or T cartesian attribute
  ! subroutine get_axis_bounds(axis,axis_bound,axes) : Return axis_bound either from an array of 
  !                                                    available axes, or defined based on axis mid-points
  ! function get_axis_modulo : Returns true if axis has the modulo attribute
  ! function get_axis_fold   : Returns is axis is folded at a boundary (non-standard meta-data)
  ! function lon_in_range    : Returns lon_strt <= longitude <= lon_strt+360
  ! subroutine tranlon       : Returns monotonic array of longitudes s.t., lon_strt <= lon(:) <= lon_strt+360.
  ! subroutine nearest_index : Return index of nearest point along axis
  !
  !</DESCRIPTION>
  !

  use mpp_io_mod, only : default_axis, default_att, atttype, axistype, mpp_get_atts, &
     mpp_modify_meta, mpp_get_axis_data
  use mpp_mod, only : mpp_error, FATAL, stdout
  use fms_mod, only : lowercase, string_array_index  

  implicit none

# include <netcdf.inc>

  public get_axis_cart, get_axis_bounds, get_axis_modulo, get_axis_fold, lon_in_range, tranlon, frac_index, nearest_index, interp_1d

  private

  integer, parameter :: maxatts = 100
  real, parameter    :: epsln= 1.e-10
  real, parameter    :: fp5 = 0.5, f360 = 360.0
  character(len=256) :: version = '$Revision$'
  character(len=256) :: tagname = '$Id$'

  interface interp_1d
     module procedure interp_1d_1d
     module procedure interp_1d_2d
     module procedure interp_1d_3d
  end interface

contains


  subroutine get_axis_cart(axis, cart)      

    type(axistype), intent(in) :: axis
    character(len=1), intent(out) :: cart
    character(len=1) :: axis_cart
    character(len=16), dimension(2) :: lon_names, lat_names
    character(len=16), dimension(3) :: z_names
    character(len=16), dimension(2) :: t_names
    character(len=16), dimension(2) :: lon_units, lat_units
    character(len=8) , dimension(4) :: z_units
    character(len=3) , dimension(4) :: t_units
    character(len=32) :: name
    integer :: i,j

    lon_names = (/'lon','x  '/)
    lat_names = (/'lat','y  '/)
    z_names = (/'depth ','height','z     '/)
    t_names = (/'time','t   '/)
    lon_units = (/'degrees_e   ', 'degrees_east'/)
    lat_units = (/'degrees_n    ', 'degrees_north'/)
    z_units = (/'cm ','m  ','pa ','hpa'/)
    t_units = (/'sec', 'min','hou','day'/)

    call mpp_get_atts(axis,cartesian=axis_cart)
    cart = 'N'

    if (axis_cart == 'x' ) cart = 'X'
    if (axis_cart == 'y' ) cart = 'Y'
    if (axis_cart == 'z' ) cart = 'Z'
    if (axis_cart == 't' ) cart = 'T'

    if (cart /= 'X' .and. cart /= 'Y' .and. cart /= 'Z' .and. cart /= 'T') then
       call mpp_get_atts(axis,name=name)
       name = lowercase(name)
       do i=1,size(lon_names)
          if (lowercase(name(1:3)) == trim(lon_names(i))) cart = 'X'
       enddo
       do i=1,size(lat_names)
          if (name(1:3) == trim(lat_names(i))) cart = 'Y'
       enddo
       do i=1,size(z_names)
          if (name == trim(z_names(i))) cart = 'Z'
       enddo
       do i=1,size(t_names)
          if (name(1:3) == t_names(i)) cart = 'T'
       enddo
    end if

    if (cart /= 'X' .and. cart /= 'Y' .and. cart /= 'Z' .and. cart /= 'T') then
       call mpp_get_atts(axis,units=name)
       name = lowercase(name)
       do i=1,size(lon_units)
          if (trim(name) == trim(lon_units(i))) cart = 'X'
       enddo
       do i=1,size(lat_units)
          if (trim(name) == trim(lat_units(i))) cart = 'Y'
       enddo
       do i=1,size(z_units)
          if (trim(name) == trim(z_units(i))) cart = 'Z'
       enddo
       do i=1,size(t_units)
          if (name(1:3) == trim(t_units(i))) cart = 'T'
       enddo
    end if

    return

  end subroutine get_axis_cart


  subroutine get_axis_bounds(axis,axis_bound,axes)

    type(axistype), intent(in) :: axis
    type(axistype), intent(out) :: axis_bound
    type(axistype), intent(in), dimension(:) :: axes

    type(atttype), dimension(:), allocatable :: att
    real, dimension(:), allocatable :: data, tmp

    integer :: i, lngth
    character(len=128) :: bounds_name = 'none', name, units
    character(len=256) :: longname
    character(len=1) :: cartesian

    axis_bound = default_axis
    allocate(att(maxatts))
    att = default_att
    call mpp_get_atts(axis,atts=att)

    do i=1,maxatts
       if (att(i)%type == NF_CHAR) then
          !            if (str_contains(att(i)%name,'bounds') .or. str_contains(att(i)%name,'edge')) then
          if (string_array_index('bounds',(/att(i)%name/)) .or. string_array_index('edge',(/att(i)%name/))) then             
             bounds_name = att(i)%catt(1:att(i)%len)
          endif
       endif
    enddo

    if (trim(bounds_name) /= 'none') then
       do i=1,size(axes)
          call mpp_get_atts(axes(i),name=name)
          if (trim(lowercase(name)) == trim(lowercase(bounds_name))) then
             axis_bound = axes(i)
          endif
       enddo
       call mpp_get_atts(axis_bound,lngth=lngth)
       if (lngth < 1) call mpp_error(FATAL,'error locating boundary axis')
    else
       call mpp_get_atts(axis,name=name,units=units,longname=longname,&
            cartesian=cartesian,lngth=lngth)
       name = trim(name)//'_bounds'
       longname = trim(longname)//' bounds'
       allocate(tmp(lngth))
       call mpp_get_axis_data(axis,tmp)
       allocate(data(lngth+1))
       do i=2,lngth
          data(i)= tmp(i-1)+fp5*(tmp(i)-tmp(i-1))
       enddo
       data(1)= tmp(1)- fp5*(tmp(2)-tmp(1))
       if (abs(data(1)) < epsln) data(1) = 0.0
       data(lngth+1)= tmp(lngth)+ fp5*(tmp(lngth)-tmp(lngth-1))         
       if (data(1) == 0.0) then
          if (abs(data(lngth+1)-360.) > epsln) data(lngth+1)=360.0
       endif
       call mpp_modify_meta(axis_bound,name=name,units=units,longname=&
            longname,cartesian=cartesian,data=data)
       deallocate(tmp)
       deallocate(data)
    endif

    return
  end subroutine get_axis_bounds

  function get_axis_modulo(axis)

    type(axistype) :: axis
    logical :: get_axis_modulo
    integer :: natt, i
    type(atttype), dimension(:), allocatable :: atts


    call mpp_get_atts(axis,natts=natt)
    allocate(atts(natt))
    call mpp_get_atts(axis,atts=atts)

    get_axis_modulo=.false.
    do i = 1,natt
       if (trim(lowercase(atts(i)%name)) == 'modulo') get_axis_modulo = .true.
    enddo

    deallocate(atts)

    return
  end function get_axis_modulo

  function get_axis_fold(axis)

    type(axistype) :: axis
    logical :: get_axis_fold
    integer :: natt, i
    type(atttype), dimension(:), allocatable :: atts


    call mpp_get_atts(axis,natts=natt)
    allocate(atts(natt))
    call mpp_get_atts(axis,atts=atts)

    get_axis_fold=.false.
    do i = 1,natt
       if (atts(i)%catt == 'fold_top') get_axis_fold = .true.
    enddo

    deallocate(atts)

    return
  end function get_axis_fold

  function lon_in_range(lon, l_strt)
    real :: lon, l_strt, lon_in_range, l_end

    lon_in_range = lon
    l_end = l_strt+360.

    if (abs(lon_in_range - l_strt) < 1.e-4) then
       lon_in_range = l_strt
       return
    endif

    if (abs(lon_in_range - l_end) < 1.e-4) then
       lon_in_range = l_strt
       return
    endif

    do
       if (lon_in_range < l_strt) then          
          lon_in_range = lon_in_range +  f360;
       else if (lon_in_range  >  l_end) then
          lon_in_range  = lon_in_range - f360;
       else
          exit
       end if
    end do

  end function lon_in_range

  subroutine tranlon(lon, lon_start, istrt)

    ! returns array of longitudes s.t.  lon_strt <= lon < lon_strt+360.
    ! also, the first istrt-1 entries are moved to the end of the array
    !
    ! e.g.
    !        lon =      0 1 2 3 4 5  ...  358 359; lon_strt = 3 ==>
    !        tranlon =  3 4 5 6 7 8  ...  359 360 361 362; istrt = 4

    real, intent(inout), dimension(:) :: lon
    real, intent(in) :: lon_start
    integer, intent(out) :: istrt


    integer :: lngth, i
    real :: lon_strt, tmp(size(lon)-1)

    lngth = size(lon)

    do i=1,lngth
       lon(i) = lon_in_range(lon(i),lon_start)
    enddo

    istrt=0
    do i=1,lngth-1
       if (lon(i+1) < lon(i)) then
          istrt=i+1 
          exit
       endif
    enddo

    if (istrt>1) then ! grid is not monotonic
       if (abs(lon(lngth)-lon(1)) < epsln) then 
          tmp = cshift(lon(1:lngth-1),istrt-1)
          lon(1:lngth-1) = tmp
          lon(lngth) = lon(1)
       else
          lon = cshift(lon,istrt-1)
       endif
       lon_strt = lon(1)
       do i=2,lngth+1
          lon(i) = lon_in_range(lon(i),lon_strt)
          lon_strt = lon(i)
       enddo
    endif

    return
  end subroutine tranlon

  function frac_index (value, array)
    !=======================================================================
    !
    !     nearest_index = index of nearest data point within "array" corresponding to
    !            "value".
    !
    !     inputs:
    !
    !     value  = arbitrary data...same units as elements in "array"
    !     array  = array of data points  (must be monotonically increasing)
    !
    !     output:
    !
    !     nearest_index =  index of nearest data point to "value"
    !             if "value" is outside the domain of "array" then nearest_index = 1
    !             or "ia" depending on whether array(1) or array(ia) is
    !             closest to "value"
    !
    !             note: if "array" is dimensioned array(0:ia) in the calling
    !                   program, then the returned index should be reduced
    !                   by one to account for the zero base.
    !
    !     example:
    !
    !     let model depths be defined by the following:
    !     parameter (km=5)
    !     dimension z(km)
    !     data z /5.0, 10.0, 50.0, 100.0, 250.0/
    !
    !     k1 = nearest_index (12.5, z, km)
    !     k2 = nearest_index (0.0, z, km)
    !
    !     k1 would be set to 2, and k2 would be set to 1 so that
    !     z(k1) would be the nearest data point to 12.5 and z(k2) would
    !     be the nearest data point to 0.0
    !
    !=======================================================================

    integer :: ia, i, ii
    real :: value, frac_index
    real, dimension(:) :: array
    logical keep_going
    
    ia = size(array)

    do i=2,ia
       if (array(i) < array(i-1)) then
          write (stdout(),*) '=> Error: "frac_index" array must be monotonically increasing when searching for nearest value to ',&
                              value
          write (stdout(),*) '          array(i) < array(i-1) for i=',i 
          write (stdout(),*) '          array(i) for i=1..ia follows:'
          do ii=1,ia
             write (stdout(),*) 'i=',ii, ' array(i)=',array(ii)
          enddo
          call mpp_error(FATAL,' "frac_index" array must be monotonically increasing.')
       endif
    enddo
    if (value < array(1) .or. value > array(ia)) then
       if (value < array(1))  frac_index = 1.
       if (value > array(ia)) frac_index = float(ia)
    else
       i=1
       keep_going = .true.
       do while (i <= ia .and. keep_going)
          i = i+1
          if (value <= array(i)) then
             frac_index = float(ia) + (value-array(i-1))/(array(i)-array(i-1)) 
             keep_going = .false.
          endif
       enddo
    endif
  end function frac_index

  function nearest_index (value, array)
    !=======================================================================
    !
    !     nearest_index = index of nearest data point within "array" corresponding to
    !            "value".
    !
    !     inputs:
    !
    !     value  = arbitrary data...same units as elements in "array"
    !     array  = array of data points  (must be monotonically increasing)
    !     ia     = dimension of "array"
    !
    !     output:
    !
    !     nearest_index =  index of nearest data point to "value"
    !             if "value" is outside the domain of "array" then nearest_index = 1
    !             or "ia" depending on whether array(1) or array(ia) is
    !             closest to "value"
    !
    !             note: if "array" is dimensioned array(0:ia) in the calling
    !                   program, then the returned index should be reduced
    !                   by one to account for the zero base.
    !
    !     example:
    !
    !     let model depths be defined by the following:
    !     parameter (km=5)
    !     dimension z(km)
    !     data z /5.0, 10.0, 50.0, 100.0, 250.0/
    !
    !     k1 = nearest_index (12.5, z, km)
    !     k2 = nearest_index (0.0, z, km)
    !
    !     k1 would be set to 2, and k2 would be set to 1 so that
    !     z(k1) would be the nearest data point to 12.5 and z(k2) would
    !     be the nearest data point to 0.0
    !
    !=======================================================================

    integer :: nearest_index, ia, i, ii
    real :: value
    real, dimension(:) :: array
    logical keep_going

    ia = size(array)

    do i=2,ia
       if (array(i) < array(i-1)) then
          write (stdout(),*) '=> Error: "nearest_index" array must be monotonically increasing when searching for nearest value to ',value
          write (stdout(),*) '          array(i) < array(i-1) for i=',i 
          write (stdout(),*) '          array(i) for i=1..ia follows:'
          do ii=1,ia
             write (stdout(),*) 'i=',ii, ' array(i)=',array(ii)
          enddo
          call mpp_error(FATAL,' "nearest_index" array must be monotonically increasing.')
       endif
    enddo
    if (value < array(1) .or. value > array(ia)) then
       if (value < array(1))  nearest_index = 1
       if (value > array(ia)) nearest_index = ia
    else
       i=1
       keep_going = .true.
       do while (i <= ia .and. keep_going)
          i = i+1
          if (value <= array(i)) then
             nearest_index = i
             if (array(i)-value > value-array(i-1)) nearest_index = i-1
             keep_going = .false.
          endif
       enddo
    endif
  end function nearest_index

  subroutine interp_1d_1d(grid_1,grid_2,data1,data2)  

    real, dimension(:), intent(in) :: grid_1, data1
    real, dimension(:), intent(inout) :: grid_2, data2

    integer :: n1, n2, i, n, ext
    real :: w

    n1 = size(grid_1)
    n2 = size(grid_2)


    do i=2,n1
       if (grid_1(i) <= grid_1(i-1)) call mpp_error(FATAL, 'grid1 not monotonic')
    enddo

    do i=2,n2
       if (grid_2(i) <= grid_2(i-1)) call mpp_error(FATAL, 'grid2 not monotonic')
    enddo

    if (grid_1(1) > grid_2(1) ) call mpp_error(FATAL, 'grid2 lies outside grid1')
    if (grid_1(n1) < grid_2(n2) ) call mpp_error(FATAL, 'grid2 lies outside grid1')

    do i=1,n2
       n = nearest_index(grid_2(i),grid_1)

       if (grid_1(n) < grid_2(i)) then
          w = (grid_2(i)-grid_1(n))/(grid_1(n+1)-grid_1(n))
          data2(i) = (1.-w)*data1(n) + w*data1(n+1)
       else
          w = (grid_2(i)-grid_1(n-1))/(grid_1(n)-grid_1(n-1))
          data2(i) = (1.-w)*data1(n-1) + w*data1(n)        
       endif
    enddo


    return

  end subroutine interp_1d_1d

  subroutine interp_1d_2d(grid_1,grid_2,data1,data2)  

    real, dimension(:,:), intent(in) :: grid_1, data1
    real, dimension(:,:), intent(inout) :: grid_2, data2

    integer :: n1, n2, i, n
    real :: w

    n1 = size(grid_1,1)
    n2 = size(grid_2,1)

    if (n1 /= n2) call mpp_error(FATAL,'grid size mismatch')

    do n=1,n1
       call interp_1d_1d(grid_1(n,:),grid_2(n,:),data1(n,:),data2(n,:))
    enddo



    return

  end subroutine interp_1d_2d

  subroutine interp_1d_3d(grid_1,grid_2,data1,data2)  

    real, dimension(:,:,:), intent(in) :: grid_1, data1
    real, dimension(:,:,:), intent(inout) :: grid_2, data2

    integer :: n1, n2, m1, m2, i, n, m
    real :: w

    n1 = size(grid_1,1)
    n2 = size(grid_2,1)
    m1 = size(grid_1,2)
    m2 = size(grid_2,2)


    if (n1 /= n2 .or. m1 /= m2) call mpp_error(FATAL,'grid size mismatch')

    do m=1,m1
       do n=1,n1
          call interp_1d_1d(grid_1(n,m,:),grid_2(n,m,:),data1(n,m,:),data2(n,m,:))
       enddo
    enddo



    return

  end subroutine interp_1d_3d


end module axis_utils_mod

