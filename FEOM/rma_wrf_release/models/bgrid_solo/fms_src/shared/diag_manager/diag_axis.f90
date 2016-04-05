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

module diag_axis_mod

!-----------------------------------------------------------------------

use mpp_domains_mod, only: domain1d, domain2d, mpp_get_compute_domain, &
                           mpp_get_domain_components, null_domain1d, &
                           null_domain2d, operator(/=)
use         fms_mod, only: error_mesg, write_version_number

implicit none
private

public  diag_axis_init, get_diag_axis, get_domain1d, get_domain2d, &
        get_axis_length, get_axis_global_length

integer, parameter :: FATAL = 2

! --- maximum number of independent axes ---

 integer, parameter :: max_axes = 100

! --- character string lengths (function of netcdf) ---

 integer, parameter :: mxch  = 128
 integer, parameter :: mxchl = 256

! --- counter of number of axes defined ---

 integer :: num_def_axes = 0

! --- defined type ---

 type diag_axis_type
      character(len=mxch)  :: name
      character(len=mxchl) :: units, long_name
      character(len=1)     :: cart_name
      real, pointer        :: data(:)
      integer              :: length, direction, edges, set
      type(domain1d) :: Domain
      type(domain2d) :: Domain2
 end type diag_axis_type

!-----------------------------------------------------------------------
!
!  name        = short name for axis
!  units       = units for axis
!  long_name   = long name for axis
!  cart_name   = cartesian axis ('x','y','z', or 't')
!  length      = number of coordinates for axis
!  direction   = if +1, data are in a positive direction (default:+1)
!                if -1, data are in a negative direction
!  data        = array of coordinate values for this axis
!
!-----------------------------------------------------------------------

! ---- storage for axis set names ----

  integer, parameter :: max_num_axis_sets = 25
  integer            :: num_axis_sets = 0
  character(len=128) :: Axis_sets(max_num_axis_sets)

! ---- global storage for all defined axes ----

 type (diag_axis_type), dimension(max_axes) :: Axes

 logical :: module_is_initialized = .FALSE.

 character(len=128) :: &
   version='$Revision$'
 character(len=128) :: &
   tagname='$Id$'

contains

!#######################################################################

 function diag_axis_init (name, data, units, cart_name, long_name,     &
                          direction, set_name, edges, Domain, Domain2) &
                  result (index)

! increment axis counter and fill in axes     

!-----------------------------------------------------------------------
!  name                = short name for axis
!  data                = array of coordinate values for this axis
!  units               = units for axis
!  cart_name           = cartesian axis ("x",'y','z','t')
!  direction(optional) = if +1, data are in a up   direction
!                      = if -1, data are in a down direction
!                      = if  0, neither up or down (default)
!  long_name(optional) = long name for axis (default: name)
!  edges    (optional) = axis id for the previously defined "edges axis"
!-----------------------------------------------------------------------

character(len=*), intent(in) :: name
real            , intent(in) :: data(:)
character(len=*), intent(in) :: units
character(len=*), intent(in) :: cart_name

character(len=*), intent(in), optional :: long_name, set_name
integer         , intent(in), optional :: direction, edges
type(domain1d)  , intent(in), optional :: Domain
type(domain2d)  , intent(in), optional :: Domain2
type(domain1d) :: domain_x, domain_y
integer ::  index, ierr, axlen

integer ::  i, set

 if ( .not.module_is_initialized ) then
      call write_( version, tagname )
 endif

!---- is there an axis set? ----

 if ( present(set_name) ) then
     set = get_axis_set_num (set_name)
!    ---- add new set name ----
       if (set == 0) then
           num_axis_sets = num_axis_sets + 1
           if (num_axis_sets > max_num_axis_sets)   &
                              call error_mesg('axis_mod',  &
                              'exceeded max_num_axis_sets', FATAL)
           set = num_axis_sets
           Axis_sets(set) = set_name
       endif
 else
     set = 0
 endif


!---- see if axis already exists ----

 do i = 1, num_def_axes
   if (name == Axes(i)%name .and. set == Axes(i)%set)     &
                        call error_mesg ('diag_axis_mod', &
                        'axis_name and axis_set already exists', FATAL)
 enddo

!---- register axis ----

 num_def_axes = num_def_axes + 1
 if (num_def_axes > max_axes) call error_mesg ('diag_axis_mod',  &
                                               'too many axes', FATAL)
 index = num_def_axes



!---- check and then save cart_name name ----

 if (cart_name == 'x' .or.cart_name == 'X' ) then
      Axes(index)%cart_name = 'X'
 else if (cart_name == 'y' .or.cart_name == 'Y' ) then
      Axes(index)%cart_name = 'Y'
 else if (cart_name == 'z' .or.cart_name == 'Z' ) then
      Axes(index)%cart_name = 'Z'
 else if (cart_name == 't' .or.cart_name == 'T' ) then
      Axes(index)%cart_name = 'T'
 else
      call error_mesg('axis_mod', 'Invalid cart_name name.', FATAL)
 endif

!---- allocate storage for coordinate values of axis ----
 axlen = size(data)
 if ( Axes(index)%cart_name == 'T' ) axlen = 0
 allocate ( Axes(index)%data(1:axlen) )

 Axes(index)%name   = name
 Axes(index)%data   = data(1:axlen)
 Axes(index)%units  = units

 Axes(index)%length = axlen
 Axes(index)%set    = set


!---- check and then save cart_name name ----

 if (cart_name == 'x' .or.cart_name == 'X' ) then
      Axes(index)%cart_name = 'X'
 else if (cart_name == 'y' .or.cart_name == 'Y' ) then
      Axes(index)%cart_name = 'Y'
 else if (cart_name == 'z' .or.cart_name == 'Z' ) then
      Axes(index)%cart_name = 'Z'
 else if (cart_name == 't' .or.cart_name == 'T' ) then
      Axes(index)%cart_name = 'T'
 else
      call error_mesg('axis_mod', 'Invalid cart_name name.', FATAL)
 endif

!---- long name (optional) default use short name ----

if (present(long_name))then
   Axes(index)%long_name = long_name
else
   Axes(index)%long_name = name
endif

!---- axis direction (-1, 0, or +1) ----

 if (present(direction))then
   if(abs(direction) > 1) call error_mesg('axis_mod',     &
                                 'direction not +/- 1',FATAL)
   Axes(index)%direction = direction
 else
   Axes(index)%direction = 0
 endif

!---- domain2d type ----

 if ( present(Domain2) ) then
    Axes(index)%Domain2 = Domain2
    call mpp_get_domain_components(Domain2, domain_x, domain_y)
    if ( Axes(index)%cart_name == 'X' ) Axes(index)%Domain = domain_x
    if ( Axes(index)%cart_name == 'Y' ) Axes(index)%Domain = domain_y
 else
    Axes(index)%Domain2 = null_domain2d 
!---- domain1d type ----

    if ( present(Domain)) then
       Axes(index)%Domain = Domain
    else
       Axes(index)%Domain = null_domain1d
    endif
 endif





!---- have axis edges been defined ? ----

     Axes(index)%edges = 0
 if (present(edges))then
   if ( edges > 0 .and. edges < num_def_axes ) then
       ierr=0
       if ( Axes(edges)%cart_name /= Axes(index)%cart_name) ierr=1
       if ( Axes(edges)%length    /= Axes(index)%length+1 ) ierr=ierr+2
       if ( Axes(edges)%set       /= Axes(index)%set      ) ierr=ierr+4
       if ( ierr > 0 )   call error_mesg ('axis_mod',     &
                        'Edges axis does not match axis', FATAL)
       Axes(index)%edges = edges
   endif
 endif

 module_is_initialized = .TRUE.

!-----------------------------------------------------------------------

end function diag_axis_init

!#######################################################################

 subroutine get_diag_axis (id, name, units, long_name, cart_name, &
                           direction, edges, Domain, data)

! Return information about the axis with index id

!-----------------------------------------------------------------------
!  id         =  axis number
!  name       = short name for axis
!  units      = units for axis
!  long_name  = long name for axis
!  cart_name  = cartesian axis ("x",'y','z', or 't')
!  direction  = if +1, data are in a positive direction (default:+1)
!               if -1, data are in a negative direction
!  data       = array of coordinate values for this axis
!-----------------------------------------------------------------------

 character(len=*), intent(out) :: name, units, long_name, cart_name
 integer, intent(in) :: id
 type(domain1d), intent(out) :: Domain
 integer, intent(out) :: direction, edges
 real, intent(out) :: data(:)

 if (id < 1 .or. id > num_def_axes) call error_mesg('diag_axis_mod',   &
                                         'illegal value for id', FATAL)

 name      = Axes(id)%name
 units     = Axes(id)%units
 long_name = Axes(id)%long_name
 cart_name = Axes(id)%cart_name
 direction = Axes(id)%direction
 edges     = Axes(id)%edges
 Domain    = Axes(id)%Domain

 if (Axes(id)%length > size(data)) call error_mesg ('diag_axis_mod', &
                                     'array data is too small', FATAL)
 data(1:Axes(id)%length) = Axes(id)%data

!-----------------------------------------------------------------------

 end subroutine get_diag_axis

!#######################################################################

 function get_axis_length (id) result (length)
   integer, intent(in) :: id
   integer             :: length
   
    if ( Axes(id)%Domain /= null_domain1d ) then
        call mpp_get_compute_domain(Axes(id)%Domain,size=length)
    else
        length = Axes(id) % length
    endif

 end function get_axis_length

!#######################################################################

 function get_axis_global_length (id) result (length)
   integer, intent(in) :: id
   integer             :: length

      length = Axes(id) % length

 end function get_axis_global_length

!#######################################################################

 function get_domain1d (id) result (Domain1)
   integer, intent(in) :: id
   type(domain1d) :: Domain1

   
   if (Axes(id)%Domain .NE. NULL_DOMAIN1D) then
      Domain1 = Axes(id)%Domain
   else
      Domain1 = NULL_DOMAIN1D
   endif

 end function get_domain1d

!#######################################################################

 function get_domain2d (ids) result (Domain2)
   integer, intent(in) :: ids(:)
   type(domain2d) :: Domain2

   integer :: i, id, flag


   if ( size(ids) < 1 .or. size(ids) > 4 ) call error_mesg  &
                          ('get_domain2d in diag_axis_mod', &
                           'input argument has incorrect size', FATAL)

   Domain2 = null_domain2d

   flag = 0
   do i = 1, size(ids)
      id = ids(i)

      if ( Axes(id)%cart_name == 'X' .or.  &
           Axes(id)%cart_name == 'Y' ) flag = flag + 1

!     --- both x/y axes found ---
      if ( flag == 2 ) then
         if (Axes(id)%Domain2 .NE. NULL_DOMAIN2D) then
            Domain2 = Axes(id)%Domain2
         endif
         exit
      endif

   enddo

 end function get_domain2d

!#######################################################################

 function get_axis_num (axis_name, set_name) result (num)

! Returns index into axis table corresponding to a given axis name

 character(len=*), intent(in)           :: axis_name
 character(len=*), intent(in), optional ::  set_name
 integer :: num, set, n

 if (present(set_name)) then
     set = get_axis_set_num (trim(set_name))
 else
     set = 0
 endif

 num = 0

 do n = 1, num_def_axes
   if ( axis_name == Axes(n)%name .and. Axes(n)%set == set ) then
     num = n
     return
   endif
 enddo

 end function get_axis_num

!#######################################################################

 function get_axis_set_num (set_name) result (num)

! Returns index in axis set table corresponding to a given axis set name

 character(len=*), intent(in) :: set_name
 integer                      :: num, naxis, iset

    num = 0
    do iset = 1, num_axis_sets
      if (set_name == Axis_sets(iset))then
        num = iset
        return
      endif
    enddo

 end function get_axis_set_num

!#######################################################################

end module diag_axis_mod
