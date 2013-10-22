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

module diag_output_mod

!-----------------------------------------------------------------------

use       mpp_io_mod, only: axistype, fieldtype, mpp_io_init,     &
                            mpp_open,  mpp_write_meta, mpp_write, &
                            mpp_flush, mpp_close, mpp_get_id,     &
                            MPP_APPEND, MPP_WRONLY, MPP_OVERWR,   &
                            MPP_NETCDF, MPP_MULTI

use  mpp_domains_mod, only: domain1d, domain2d, mpp_define_domains, mpp_get_pelist,&
                            mpp_get_global_domain, mpp_get_compute_domains, &
                            null_domain1d, null_domain2d, operator(/=), mpp_get_layout
use mpp_mod, only : mpp_npes, mpp_pe
use    diag_axis_mod, only: diag_axis_init, get_diag_axis,           &
                            get_axis_length, get_axis_global_length, &
                            get_domain1d, get_domain2d

use time_manager_mod, only: get_calendar_type, THIRTY_DAY_MONTHS, &
                            JULIAN, GREGORIAN, NOLEAP

use          fms_mod, only: error_mesg, mpp_pe, write_version_number

use platform_mod, only: r8_kind

implicit none
private

public :: diag_output_init, write_axis_meta_data, write_field_meta_data, &
          done_meta_data, diag_field_out, diag_output_end, diag_flush, &
          diag_fieldtype, NETCDF

!-----------------------------------------------------------------------
!------------------------- interfaces ----------------------------------
!
!   diag_output_init:   initializes output file and time axis
!
!   write_axis_meta_data:  writes meta data for axes
!                           called for each field in a file
!
!   write_field_meta_data:  writes meta data for fields and returns a
!                             fieldtype, called for each field in a file
!
!   done_meta_data:   writes axis data, called once per file after all
!                       "write_meta_data" calls and before the first
!                       "diag_field_out" call
!
!   diag_field_out:   writes field data to an output file
!
!   diag_flush:       called periodically to flush buffer and insure
!                       that data is not lost if execution fails
!
!   diag_output_end:  called once to write all data at the end of
!                       a program
!
!------------------------ data type ------------------------------------
!
!   fieldtype:   made public from mpp_io_mod,
!                  needed by several interfaces in diag_output_mod
!
!------------------------ data -----------------------------------------
!
!   NETCDF:    valid format types, needed in "diag_output_init"
!
!-----------------------------------------------------------------------

integer, parameter :: NETCDF = 1
integer, parameter :: FATAL = 2

integer, parameter :: mxch  = 128
integer, parameter :: mxchl = 256

integer :: current_file_unit = -1

integer, dimension(2,2) :: max_range
data max_range / -32767, 32767, &
                   -127,   127 /

integer, dimension(2) :: missval = (/ -32768, -128 /)

integer, parameter :: max_axis_num = 20
integer :: num_axis_in_file = 0
integer                :: axis_in_file   (max_axis_num)
logical                :: time_axis_flag (max_axis_num)
logical                :: edge_axis_flag (max_axis_num)
type(axistype)         :: Axis_types     (max_axis_num)

!-----------------------------------------------------------------------

type diag_fieldtype
   type(fieldtype)         :: Field
   type(domain2d)          :: Domain
   real                    :: miss, miss_pack
   logical                 :: miss_present, miss_pack_present
end type

!-----------------------------------------------------------------------

logical :: module_is_initialized = .FALSE.

character(len=128), private :: version= &
  '$Revision$'
character(len=128), private :: tagname= &
  '$Id$'

contains

!#######################################################################

subroutine diag_output_init ( file_name, format, file_title,  &
                            time_name, time_units,          &
                            file_unit, time_id )

 character(len=*), intent(in)  :: file_name, file_title,  &
                                  time_name, time_units
 integer         , intent(in)  :: format
 integer         , intent(out) :: file_unit, time_id

!-----------------------------------------------------------------------
!
!        Registers the time axis and opens the output file
!
! INPUT: file_name  = output file name (character, max len=128)
!        format     = file format (integer, only type NETCDF allowed)
!        file_title = descriptive title for the file (char, max len=128)
!        time_name  = name of the time axis (character, max len=128)
!                       note: this name will also be used for the
!                             long name
!        time_units = units for the time axis (character, max len=128)
!                       note: this string may contain date information,
!                             e.g., days since 1979-01-01 00:00:00.0
!
! OUTPUT: file_unit = file unit assign to the output file, needed
!                     for subsequent calls to diag_output_mod (integer)
!         time_id   = axis id assigned to the time axis (integer)
!
!-----------------------------------------------------------------------

 real(KIND=r8_kind), dimension(1) :: tdata
 integer :: form, threading, fileset

!-----------------------------------------------------------------------
!---- initialize mpp_io ----

   if ( .not.module_is_initialized ) then
        call mpp_io_init ()
        module_is_initialized = .TRUE.
   endif
   call write_version_number( version, tagname )

!---- register time axis (use file_name as axis set_name) ----

   time_id = diag_axis_init ( time_name, tdata, time_units, 'T',  &
                              time_name, set_name=file_name )

!---- set up output file ----

   select case (format)
     case (NETCDF)
         form      = MPP_NETCDF
         threading = MPP_MULTI
         fileset   = MPP_MULTI
     case default
         call error_mesg ('diag_output_init', 'invalid format', FATAL)
   end select

!---- open output file (return file_unit id) -----

!!!call mpp_open ( file_unit, file_name, action=MPP_WRONLY,        &
   call mpp_open ( file_unit, file_name, action=MPP_OVERWR,        &
                   form=form, threading=threading, fileset=fileset )

!---- write global attributes ----

   if ( file_title(1:1) /= ' ' ) then
       call mpp_write_meta ( file_unit, 'title', cval=trim(file_title))
   endif

!-----------------------------------------------------------------------

end subroutine diag_output_init

!#######################################################################

subroutine write_axis_meta_data ( file_unit, axes )

 integer         ,  intent(in)  :: file_unit, axes(:)

!-----------------------------------------------------------------------
!
! INPUT: file_name  = output file name (character, max len=128)
!        axes       = array of axis id's (including the time axis)
!                       (integer, dimension(:))
!
!-----------------------------------------------------------------------

 type(domain1d)          :: Domain
 type(domain1d)          :: Edge_Domain

   character(len=mxch)  :: axis_name, axis_units
   character(len=mxchl) :: axis_long_name
   character(len=1)     :: axis_cart_name
   integer              :: axis_direction, axis_edges
   real, allocatable    :: axis_data(:)
   integer, allocatable    :: axis_extent(:), pelist(:)

   integer :: calendar, id_axis, id_time_axis
   integer :: i, indx, num, length, edges_index
   integer :: gbegin, gend, gsize, ndivs
   logical :: use_range

!-----------------------------------------------------------------------
!---- save the current file_unit ----

 if ( num_axis_in_file == 0 ) current_file_unit = file_unit

!-----------------------------------------------------------------------
!---- dummy checks ----

 num = size(axes)
 if ( num < 1 ) call error_mesg ( 'write_axis_meta_data', &
                                  'number of axes < 1', FATAL)

 if ( file_unit /= current_file_unit ) call error_mesg  &
           ( 'write_axis_meta_data',  &
             'writing meta data out-of-order to different files', FATAL)

!-----------------------------------------------------------------------
!---- check all axes ----
!---- write axis meta data for new axes ----

 do i = 1, num

!-----------------------------------------------------------------------

   id_axis = axes(i)
   indx = get_axis_index ( id_axis )

!---- skip axes already written -----

   if ( indx > 0 ) cycle

!---- create new axistype (then point to) -----

   num_axis_in_file = num_axis_in_file + 1
   axis_in_file(num_axis_in_file) = id_axis
   edge_axis_flag(num_axis_in_file) = .false.

   length = get_axis_global_length ( id_axis )
   allocate ( axis_data(length) )

   call get_diag_axis ( id_axis, axis_name, axis_units,     &
                        axis_long_name, axis_cart_name,     &
                        axis_direction, axis_edges, Domain, &
                        axis_data )

 if ( Domain .ne. null_domain1d ) then
   if (length > 0) then
      call mpp_write_meta ( file_unit, Axis_types(num_axis_in_file),  &
                            axis_name,      axis_units,       &
                            axis_long_name, axis_cart_name,   &
                            axis_direction, Domain, axis_data )
   else
      call mpp_write_meta ( file_unit, Axis_types(num_axis_in_file),  &
                            axis_name,      axis_units,       &
                            axis_long_name, axis_cart_name,   &
                            axis_direction, Domain            )
   endif
 else
   if (length > 0) then
      call mpp_write_meta ( file_unit, Axis_types(num_axis_in_file),  &
                            axis_name,      axis_units,       &
                            axis_long_name, axis_cart_name,   &
                            axis_direction,    data=axis_data )
   else
      call mpp_write_meta ( file_unit, Axis_types(num_axis_in_file),  &
                            axis_name,      axis_units,       &
                            axis_long_name, axis_cart_name,   &
                            axis_direction                    )
   endif
 endif

!---- write additional attribute (calendar_type) for time axis ----
!---- NOTE: calendar attribute is compliant with CF convention 
!---- http://www.cgd.ucar.edu/cms/eaton/netcdf/CF-current.htm#cal

   if ( axis_cart_name == 'T' ) then
      time_axis_flag (num_axis_in_file) = .true.
      id_time_axis = mpp_get_id( Axis_types(num_axis_in_file) )
      calendar = get_calendar_type ( )
      select case (calendar)
         case (THIRTY_DAY_MONTHS)
            call mpp_write_meta ( file_unit, id_time_axis, &
                                 'calendar_type', cval='THIRTY_DAY_MONTHS')
            call mpp_write_meta ( file_unit, id_time_axis, &
                                 'calendar', cval='360')
         case (JULIAN)
            call mpp_write_meta ( file_unit, id_time_axis, &
                                 'calendar_type', cval='JULIAN')
            call mpp_write_meta ( file_unit, id_time_axis, &
                                 'calendar', cval='JULIAN')
         case (GREGORIAN)
            call mpp_write_meta ( file_unit, id_time_axis, &
                                 'calendar_type', cval='GREGORIAN')
            call mpp_write_meta ( file_unit, id_time_axis, &
                                 'calendar', cval='GREGORIAN')
         case (NOLEAP)
            call mpp_write_meta ( file_unit, id_time_axis, &
                                 'calendar_type', cval='NOLEAP')
            call mpp_write_meta ( file_unit, id_time_axis, &
                                 'calendar', cval='NOLEAP')
      end select
   else
      time_axis_flag (num_axis_in_file) = .false.
   endif

   deallocate (axis_data)

!-----------------------------------------------------------------------
!  ------------- write axis containing edge information ---------------

!  --- this axis has no edges -----
   if ( axis_edges <= 0 ) cycle

!  --- was this axis edge previously defined? ---
   id_axis = axis_edges
   edges_index = get_axis_index ( id_axis )
   if ( edges_index > 0 ) cycle

!  ---- get data for axis edges ----

   length = get_axis_global_length ( id_axis )
   allocate ( axis_data(length) )

   call get_diag_axis ( id_axis, axis_name, axis_units,     &
                        axis_long_name, axis_cart_name,     &
                        axis_direction, axis_edges, Domain, &
                        axis_data )

!  ---- write edges attribute to original axis ----

   call mpp_write_meta ( file_unit, mpp_get_id(Axis_types(num_axis_in_file)), &
                         'edges', cval=axis_name )

!  ---- add edges index to axis list ----
!  ---- assume this is not a time axis ----

   num_axis_in_file = num_axis_in_file + 1
   axis_in_file(num_axis_in_file) = id_axis
   edge_axis_flag(num_axis_in_file) = .true.
   time_axis_flag (num_axis_in_file) = .false.

!  ---- write edges axis to file ----

 if ( Domain /= null_domain1d ) then
! assume domain decomposition is irregular and loop through all prev and next
! domain pointers extracting domain extents.  Assume all pes are used in
! decomposition
   call mpp_get_global_domain(Domain,begin=gbegin,end=gend,size=gsize)
   call mpp_get_layout(Domain, ndivs)
   if (ndivs .EQ. 1) then
      call mpp_write_meta ( file_unit, Axis_types(num_axis_in_file),  &
                         axis_name,      axis_units,       &
                         axis_long_name, axis_cart_name,   &
                         axis_direction,  data=axis_data )
   else
      if (ALLOCATED(axis_extent)) deallocate(axis_extent)
      allocate(axis_extent(0:ndivs-1))
      call mpp_get_compute_domains(Domain,size=axis_extent(0:ndivs-1))
      gend=gend+1
      axis_extent(ndivs-1)= axis_extent(ndivs-1)+1
      if (ALLOCATED(pelist)) deallocate(pelist)      
      allocate(pelist(0:ndivs-1))
      call mpp_get_pelist(Domain,pelist)
      call mpp_define_domains((/gbegin,gend/),ndivs,Edge_Domain,&
                           pelist=pelist(0:ndivs-1), &
                           extent=axis_extent(0:ndivs-1))
      call mpp_write_meta ( file_unit, Axis_types(num_axis_in_file),  &
                         axis_name,      axis_units,       &
                         axis_long_name, axis_cart_name,   &
                         axis_direction, Edge_Domain,  data=axis_data )
    endif
 else

   call mpp_write_meta ( file_unit, Axis_types(num_axis_in_file),  &
                         axis_name,      axis_units,       &
                         axis_long_name, axis_cart_name,   &
                         axis_direction, data=axis_data )

 endif


   deallocate (axis_data)

!-----------------------------------------------------------------------

 enddo

!-----------------------------------------------------------------------

end subroutine write_axis_meta_data

!#######################################################################

function write_field_meta_data ( file_unit, name, axes, units,      &
                          long_name, range, pack, mval, avg_name )  &
                 result ( Field )

 integer         ,  intent(in)  :: file_unit, axes(:)
 character(len=*),  intent(in)  :: name, units, long_name
 real,    optional, intent(in)  :: range(2), mval
 integer, optional, intent(in)  :: pack
 character(len=*), optional, intent(in)  :: avg_name

 type(diag_fieldtype)          :: Field

!-----------------------------------------------------------------------
!
! INPUT: file_name  = output file name (character, max len=128)
!        name       = field name (character, max len=128)
!        axes       = array of axis id's (including the time axis)
!                       (integer, dimension(4))
!        units      = field units (character, max len=128)
!        long_name  = field long_name (character, max len=256)
!
! OPTIONAL INPUT:
!        range      = valid range (min, max), if min > max the range
!                        will be ignored (real, dimension(2))
!        pack       = packing flag, only valid when range specified
!                      (1=64 bit, 2=32 bit, 4=16 bit, 8=8 bit)
!                      (integer)
!        mval       = missing value, must be within valid range (real)
!        avg_name   = name of variable containing time averaging info
!                      (character, max len=128)
!
! RETURNS:  Field = diag_fieldtype, will be needed for
!                     subsequent calls to diag_output_mod
!
!-----------------------------------------------------------------------

   real    :: rscale, add
   integer :: i, indx, num, ipack, np
   logical :: use_range

   integer :: axis_indices(4)

!-----------------------------------------------------------------------
!---- dummy checks ----

 num = size(axes)
 if ( num < 1 ) call error_mesg ( 'write_meta_data', &
                                  'number of axes < 1', FATAL)
 if ( num > 4 ) call error_mesg ( 'write_meta_data', &
                                  'number of axes > 4', FATAL)

 if ( file_unit /= current_file_unit ) call error_mesg  &
           ( 'write_meta_data',  &
             'writing meta data out-of-order to different files', FATAL)


!---- check all axes for this field ----
!---- set up indexing to axistypes ----

 do i = 1, num

   indx = get_axis_index ( axes(i) )

!---- point to existing axistype -----

   if ( indx > 0 ) then
       axis_indices(i) = indx
   else
       call error_mesg ('write_field_meta_data',   &
                 'axis data not written for field '//trim(name), FATAL)
   endif

 enddo

!-----------------------------------------------------------------------
!--------------------- write field meta data ---------------------------

!       ---- select packing? ----
!            (packing option only valid with range option)

   ipack = 2
   if (present(pack))  ipack = pack

!---- check range ----

   use_range = .false.
   add   = 0.0
   rscale = 1.0

   if ( present(range) ) then
        if ( range(2) > range(1) ) then
             use_range = .true.
!            ---- set packing parameters ----
             if ( ipack > 2 ) then
                np = ipack/4
                add   = 0.5*(range(1)+range(2))
                rscale = (range(2)-range(1))/ &
                        real(max_range(2,np)-max_range(1,np))
             endif
        endif
   endif

!       ---- select packing? ----

   if ( present(mval) ) then
      Field%miss = mval
      Field%miss_present = .true.
      if (ipack > 2 ) then
        np = ipack/4
        Field%miss_pack = real(missval(np))*rscale+add
        Field%miss_pack_present = .true.
      else
        Field%miss_pack = mval
        Field%miss_pack_present = .false.
      endif
   else
      Field%miss_present      = .false.
      Field%miss_pack_present = .false.
   endif

!------ write meta data and return fieldtype -------

   if ( use_range ) then

     if ( Field%miss_present ) then

         call mpp_write_meta ( file_unit, Field%Field,            &
                               Axis_types(axis_indices(1:num)),   &
                               name, units, long_name,            &
                               range(1), range(2),                &
                               missing=Field%miss_pack,           &
                               scale=rscale, add=add, pack=ipack   )
     else

         call mpp_write_meta ( file_unit, Field%Field,            &
                               Axis_types(axis_indices(1:num)),   &
                               name, units,  long_name,           &
                               range(1), range(2),                &
                               scale=rscale, add=add, pack=ipack   )
     endif

   else

     if ( Field%miss_present ) then

         call mpp_write_meta ( file_unit, Field%Field,            &
                               Axis_types(axis_indices(1:num)),   &
                               name, units, long_name,            &
                               missing=Field%miss_pack,           &
                               pack=ipack                         )
     else

         call mpp_write_meta ( file_unit, Field%Field,            &
                               Axis_types(axis_indices(1:num)),   &
                               name, units, long_name,            &
                               pack=ipack                         )
     endif

   endif

!---- write additional attribute for time averaging -----

   if ( present(avg_name) ) then
      if ( avg_name(1:1) /= ' ' ) then
          call mpp_write_meta ( file_unit, mpp_get_id(Field%Field),       &
                                'time_avg_info',                 &
                                cval=trim(avg_name) // '_T1,' // &
                                     trim(avg_name) // '_T2,' // &
                                     trim(avg_name) // '_DT'  )
      endif
   endif

!---- get axis domain ----


   Field%Domain = get_domain2d ( axes )

!-----------------------------------------------------------------------

end function write_field_meta_data

!#######################################################################

subroutine done_meta_data (file_unit)

  integer,  intent(in)  :: file_unit
  
  type(domain1d)  :: Domain
  type(domain1d)  :: Edge_Domain
  integer, allocatable    :: axis_extent(:)
  integer :: i, gbegin, gend, gsize

!---- write data for all non-time axes ----

   do i = 1, num_axis_in_file
       if (time_axis_flag(i)) cycle
       call mpp_write ( file_unit, Axis_types(i) )
   enddo

   num_axis_in_file = 0

end subroutine done_meta_data

!#######################################################################

subroutine diag_field_out ( file_unit, Field, data, time )

integer             , intent(in)    :: file_unit
type(diag_fieldtype), intent(inout)    :: Field
real                , intent(inout) :: data(:,:,:)
real(KIND=r8_kind), optional , intent(in)    :: time

!---- replace original missing value with (un)packed missing value ----
!print *, 'PE,name,miss_pack_present=',mpp_pe(), &
!  trim(Field%Field%name),Field%miss_pack_present

   if ( Field%miss_pack_present ) then
!print *, 'PE,name,miss,miss_pack',mpp_pe(), &
!  trim(Field%Field%name),Field%miss,Field%miss_pack
        where ( data == Field%miss ) data = Field%miss_pack
   endif

!---- output data ----

  if ( Field%Domain /= null_domain2d ) then
      call mpp_write (file_unit, Field%Field, Field%Domain, data, time)
  else
      call mpp_write (file_unit, Field%Field, data, time)
  endif

end subroutine diag_field_out

!#######################################################################

subroutine diag_flush (file_unit)

integer        , intent(in) :: file_unit

call mpp_flush (file_unit)

end subroutine diag_flush

!#######################################################################

subroutine diag_output_end (file_unit)

integer        , intent(in) :: file_unit

if (.not. module_is_initialized ) return

call mpp_close (file_unit)

module_is_initialized = .FALSE.

end subroutine diag_output_end

!#######################################################################

function get_axis_index ( num ) result ( indx )

  integer, intent(in) :: num
  integer             :: indx
  integer             :: i

!---- get the array index for this axis type ----
!---- set up pointers to axistypes ----
!---- write axis meta data for new axes ----

   indx = 0
   do i = 1, num_axis_in_file
      if ( num == axis_in_file(i) ) then
          indx = i
          exit
      endif
   enddo

end function get_axis_index

!#######################################################################

!subroutine modify_domain_index ( Domain, add )
!
!  type(domain1d), pointer :: Domain
!  integer, intent(in) :: add
!
!   if (.not.associated(Domain)) return
!
!   if ( Domain%Compute%end_index == Domain%Global%end_index ) then
!        Domain%Compute%end_index = Domain%Compute%end_index + add
!        Domain%Compute%size      = Domain%Compute%size      + add
!   endif
!
!        Domain%Global %end_index = Domain%Global %end_index + add
!        Domain%Global %size      = Domain%Global %size      + add
!
!end subroutine modify_domain_index

!#######################################################################

end module diag_output_mod

