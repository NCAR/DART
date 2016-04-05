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
module diag_manager_mod

! <CONTACT EMAIL="mh@gfdl.noaa.gov">
!   Matt Harrison
! </CONTACT>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

! <OVERVIEW>
!   <TT>diag_manager_mod</TT> is a set of simple calls for parallel diagnostics on
!   distributed systems. It is geared toward the writing of data in netCDF format.
! </OVERVIEW>

! <DESCRIPTION>
!   <TT>diag_manager_mod</TT> provides a convenient set of interfaces for
!   writing data to disk.  It is built upon the parallel I/O interface
!   <LINK SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/shared/mpp/mpp_io.html"><TT>mpp_io</TT></LINK>.  A single
!   group of calls to the <TT>diag_manager_mod</TT> interfaces provides data to disk
!   at any number of sampling and/or averaging intervals specified at run-time.
!   Run-time specification of diagnostics are input through the diagnostics table,
!   which is described in the
!   <LINK SRC="diag_table_tk.html">diag_table_tk</LINK> documentation.
!   
!   <B>Features of <TT>diag_manager_mod</TT> include:</B>
!   Simple, minimal API.<BR/>
!   Run-time choice of diagnostics.<BR/>
!   Self-describing files: comprehensive header information
!   (metadata) in the file itself.<BR/>
!   Strong parallel write performance.<BR/>
!   Integrated netCDF capability: <LINK
!   SRC="http://www.unidata.ucar.edu/packages/netcdf/">netCDF</LINK> is a
!   data format widely used in the climate/weather modeling
!   community. netCDF is considered the principal medium of data storage
!   for <TT>diag_manager_mod</TT>. Raw unformatted
!   fortran I/O capability is also available.<BR/>
!   Requires off-line post-processing: a tool for this purpose,
!   <TT>mppnccombine</TT>, is available. GFDL users may use
!   <TT>~hnv/pub/mppnccombine</TT>. Outside users may obtain the
!   source <LINK
!   SRC="ftp://ftp.gfdl.gov/perm/hnv/mpp/mppnccombine.c">here</LINK>.  It
!   can be compiled on any C compiler and linked with the netCDF
!   library. The program is free and is covered by the <LINK
!   SRC="ftp://ftp.gfdl.gov/perm/hnv/mpp/LICENSE">GPL license</LINK>.
! </DESCRIPTION>
use time_manager_mod, only: get_time, set_time, get_date, set_date,    &
                            increment_date, operator(-), operator(>=), &
                            operator(>), operator(<), operator(==),    &
                            time_type, increment_time, month_name,     &
                            get_calendar_type, NO_CALENDAR

use       mpp_io_mod, only: mpp_open, MPP_RDONLY

use          fms_mod, only: error_mesg, FATAL, WARNING, NOTE,          &
                            close_file, stdlog, write_version_number,  &
                            file_exist, mpp_pe                         

use    diag_axis_mod, only: diag_axis_init, get_axis_length

use  diag_output_mod, only: diag_fieldtype, diag_output_init, write_axis_meta_data,  &
                            write_field_meta_data, done_meta_data, &
                            diag_field_out, diag_output_end,       &
                            diag_flush 

implicit none
private

!public  diag_manager_init, send_data, diag_manager_end,  &
!        register_diag_field, register_static_field, &
!        diag_axis_init, get_base_time, get_base_date, need_data

! Specify storage limits for fixed size tables used for pointers, etc.
integer, parameter :: max_fields_per_file = 150
integer, parameter :: max_out_per_in_field = 10
integer, parameter :: max_files = 10
integer ::            num_files = 0
integer, parameter :: max_input_fields = 200
integer ::            num_input_fields = 0
integer, parameter :: max_output_fields = 200
integer ::            num_output_fields = 0

! Global data for all files
type (time_type) :: base_time
integer :: base_year, base_month, base_day, base_hour, base_minute, base_second
character(len = 256) :: global_descriptor

type file_type
   character(len=128) :: name
   integer :: output_freq
   integer :: output_units
   integer :: format
   integer :: time_units
   character(len=128) :: long_name
   integer :: fields(max_fields_per_file)
   integer :: num_fields
   integer :: file_unit
   integer :: time_axis_id
   type (time_type) :: last_flush
   type(diag_fieldtype) :: f_avg_start, f_avg_end, f_avg_nitems
end type file_type

type input_field_type
   character(len=128) :: module_name, field_name, long_name, units
   integer :: axes(3)
   integer :: num_axes
   logical :: missing_value_present, range_present
   real :: missing_value, range(2)
   integer :: output_fields(max_out_per_in_field)
   integer :: num_output_fields
   integer :: size(3), total_elements
   logical :: static, register
end type input_field_type

type output_field_type
   integer :: input_field
   integer :: output_file
   character(len=128) :: output_name
   logical :: time_average
   integer :: pack
   real, pointer :: buffer(:, :, :)
   type(time_type) :: last_output, next_output, next_next_output
   integer :: total_elements, num_elements
   logical :: static
!   type(diag_fieldtype), pointer :: f_type
   type(diag_fieldtype) :: f_type
   integer :: axes(3), num_axes
end type output_field_type


type (file_type) :: files(max_files)
type (input_field_type) :: input_fields(max_input_fields)
type (output_field_type) :: output_fields(max_output_fields)

logical :: first_send_data_call = .true.
logical :: module_is_initialized = .false.

integer, parameter :: EVERY_TIME =  0
integer, parameter :: END_OF_RUN = -1

integer, parameter :: DIAG_SECONDS = 1, DIAG_MINUTES = 2, DIAG_HOURS = 3
integer, parameter :: DIAG_DAYS = 4, DIAG_MONTHS = 5, DIAG_YEARS = 6
character (len=10) :: time_unit_list(6) = (/'seconds   ', 'minutes   ', &
   'hours     ', 'days      ', 'months    ', 'years     '/)

character (len = 7) :: avg_name = 'average'

! version number of this module
  character(len=128) :: version = '$Revision$'
  character(len=128) :: tagname = '$Id$'

! <INTERFACE NAME="send_data">

!   <OVERVIEW>
!     Send data over to output fields. 
!   </OVERVIEW>
!   <DESCRIPTION>
!     send_data is overloaded for 1 to 3-d arrays.
!     diag_field_id corresponds to the id returned from a previous call to
!     register_diag_field.  The field array is restricted to the computational
!     range of the array. Optional argument is_in can be used to update
!     sub-arrays of the entire field.  Additionally, an optional logical or real
!     mask can be used to apply missing values to the array.  For the real
!     mask, the mask is applied if the mask value is less than 0.5.  The
!     weight array is currently not implemented.
!   </DESCRIPTION>
!   <IN NAME="diag_field_id" TYPE="integer"  > </IN>
!   <IN NAME="field" TYPE="real" DIM="(:,:,:)" > </IN>
!   <IN NAME="time" TYPE="time_type"  > </IN>

interface send_data
   module procedure send_data_0d
   module procedure send_data_1d
   module procedure send_data_2d
   module procedure send_data_3d
end interface
! </INTERFACE>

contains

!-------------------------------------------------------------------------

! <FUNCTION NAME="register_diag_field">

!   <OVERVIEW>
!     Register Diagnostic Field.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Return field index for subsequent calls to <LINK SRC="#send_data"> send_data </LINK>
!   </DESCRIPTION>
!   <TEMPLATE>
!     register_diag_field(module_name, field_name, axes, init_time, &
!     long_name, units, missing_value, range)
!   </TEMPLATE>

!   <IN NAME="module_name" TYPE="character(len=*)"> </IN>
!   <IN NAME="field_name" TYPE="character(len=*)"> </IN>
!   <IN NAME="axes" TYPE="integer" DIM="(:)"> </IN>
!   <IN NAME="init_time" TYPE="time_type"> </IN>
!   <IN NAME="long_name" TYPE="character(len=*)"> </IN>
!   <IN NAME="units" TYPE="character(len=*)"> </IN>
!   <IN NAME="missing_value" TYPE="real"> </IN>
!   <IN NAME="range" TYPE="real" DIM="(2)"> </IN>

function register_diag_field(module_name, field_name, axes, init_time, &
   long_name, units, missing_value, range)

! Indicates the calling modules intent to supply data for this field.

integer ::  register_diag_field
character(len=*), intent(in) :: module_name, field_name
integer, intent(in) :: axes(:)
type(time_type), intent(in) :: init_time
character(len=*), optional, intent(in) :: long_name, units
real, optional, intent(in) :: missing_value, range(2)

integer :: field, num_axes, i_size, j_size, k_size, j, ind, file_num, freq
integer :: total_elements, output_units

! Call register static, then set static back to false
register_diag_field = register_static_field(module_name, field_name, axes, &
   long_name, units, missing_value, range)

if(register_diag_field >0) then
   input_fields(register_diag_field)%static = .false.

   field = register_diag_field


! Need set output fields to not static, too
! Need to loop through all output_fields associated and allocate their buffers
! Would have to consider horizontal operations at this point
!   num_output_fields = input_fields(field)%num_output_fields
   do j = 1, input_fields(field)%num_output_fields
!! do j = 1, num_output_fields
      ind = input_fields(field)%output_fields(j)
      output_fields(ind)%static = .false.
! Set up times in output_fields
      output_fields(ind)%last_output = init_time
! Get output frequency from for the appropriate output file
      file_num = output_fields(ind)%output_file
      freq = files(file_num)%output_freq
      output_units = files(file_num)%output_units
      output_fields(ind)%next_output = &
         diag_time_inc(init_time, freq, output_units)
      output_fields(ind)%next_next_output = &
         diag_time_inc(output_fields(ind)%next_output, freq, output_units)
   end do

endif

end function register_diag_field
! </FUNCTION>

!-------------------------------------------------------------------------
! <FUNCTION NAME="register_static_field">

!   <OVERVIEW>
!     Register Static Field.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Return field index for subsequent call to send_data.
!   </DESCRIPTION>
!   <TEMPLATE>
!     register_static_field(module_name, field_name, axes, &
!     long_name, units, missing_value, range, require)
!   </TEMPLATE>

!   <IN NAME="module_name" TYPE="character(len=*)"> </IN>
!   <IN NAME="field_name" TYPE="character(len=*)"> </IN>
!   <IN NAME="axes" TYPE="integer" DIM="(:)"> </IN>
!   <IN NAME="long_name" TYPE="character(len=*)"> </IN>
!   <IN NAME="units" TYPE="character(len=*)"> </IN>
!   <IN NAME="missing_value" TYPE="real"> </IN>
!   <IN NAME="range" TYPE="real" DIM="(2)"> </IN>

function register_static_field(module_name, field_name, axes, &
   long_name, units, missing_value, range, require)

integer register_static_field
character(len=*), intent(in) :: module_name, field_name
integer, intent(in) :: axes(:)
character(len=*), optional, intent(in) :: long_name, units
real, optional, intent(in) :: missing_value, range(2)
logical, optional, intent(in) :: require  ! require static field to be in every file
                                          ! e.g. 2-d axes
integer :: field, num_axes, j, out_num, file_num, freq, siz(3)
character(len=128) :: error_string

if (.not.module_is_initialized) call error_mesg ('diag_manager_mod',  &
                       'module has not been initialized', FATAL)

register_static_field = find_input_field(module_name, field_name)

if (PRESENT(require)) then
   if (require) then
      call init_input_field(module_name, field_name)
      register_static_field = find_input_field(module_name, field_name)
      do j=1, num_files
! need to think about what to do if the axes are not present, e.g. file only
! contains data slices
         call init_output_field(module_name, field_name,field_name, &
                               files(j)%name,.false.,2)
      enddo
   endif
endif

! Negative index returned if this field is not used in table
if(register_static_field < 0) return

! Store information for this input field into input field table
field = register_static_field

! Set static to true, if called by register_diag_field this is flipped to false
input_fields(field)%static = .true.
! Set flag that this field was registered
input_fields(field)%register = .true.

! Store the axis info
num_axes = size(axes)
input_fields(field)%axes(1:num_axes) = axes
input_fields(field)%num_axes = num_axes

! Need to check for present, otherwise defaults
if(present(long_name)) then
   input_fields(field)%long_name = trim(long_name)
else
   input_fields(field)%long_name = input_fields(field)%field_name
endif
if(present(units)) then
   input_fields(field)%units = trim(units)
else
   input_fields(field)%units = 'none'
endif
if(present(missing_value)) then
! ??? What about passing missing value through to the output_fields???
   input_fields(field)%missing_value = missing_value
   input_fields(field)%missing_value_present = .true.
else
   input_fields(field)%missing_value_present = .false.
endif
if(present(range)) then
   input_fields(field)%range = range
   input_fields(field)%range_present = .true.
else
   input_fields(field)%range = (/ 1., 0. /)
   input_fields(field)%range_present = .false.
endif

! Next need to compute the size of the domain for this set of axes
! Comes from axis handler; assume result is i_size, j_size, k_size
! for three-d problem
siz = 1
do j = 1, num_axes
   siz(j) = get_axis_length(axes(j))
end do

! Default length for axes is 1
do j = 1, 3
   input_fields(field)%size(j) = siz(j)
end do
input_fields(field)%total_elements = siz(1) * siz(2) * siz(3)

! Do not allow zero size variable through ?
if (input_fields(field)%total_elements <= 0) then
     write (error_string,11) trim(field_name)
  11 format ('field ',a,' has zero size')
     call error_mesg ('register_static_field in diag_manager_mod', &
                      trim(error_string), WARNING)
endif

! Need to loop through all output_fields associated and allocate their buffers
! Would have to consider horizontal operations at this point
!num_output_fields = input_fields(field)%num_output_fields
do j = 1, input_fields(field)%num_output_fields
!do j = 1, num_output_fields
   out_num = input_fields(field)%output_fields(j)
   allocate(output_fields(out_num)%buffer(siz(1), siz(2), siz(3)))
   output_fields(out_num)%buffer = 0.0
   output_fields(out_num)%total_elements = input_fields(field)%total_elements
   output_fields(out_num)%num_elements = 0
! Reset to false in register_field if this is not static
   output_fields(out_num)%static = .true.
! For now, axes are copied straight across to output fields
   output_fields(out_num)%axes = input_fields(field)%axes
   output_fields(out_num)%num_axes = input_fields(field)%num_axes
end do


end function register_static_field
! </FUNCTION>

!-------------------------------------------------------------------------
! <FUNCTION NAME="send_data_0d" INTERFACE="send_data">
!   <IN NAME="diag_field_id" TYPE="integer"  > </IN>
!   <IN NAME="field" TYPE="real"  > </IN>
!   <IN NAME="time" TYPE="time_type"  > </IN>
! </FUNCTION>

function send_data_0d(diag_field_id, field, time, mask, rmask, weight)

logical :: send_data_0d
integer, intent(in) :: diag_field_id
real, intent(in) :: field
type (time_type), intent(in) :: time
logical, optional :: mask
real, optional :: rmask, weight

real :: field_out(1, 1, 1)
real :: weight_out(1, 1, 1)
logical :: mask_out(1, 1, 1)


! First copy the data to a three d array with last element 1
field_out(1, 1, 1) = field

! Default values for mask and weight
mask_out = .true.
weight = 1.0

if(present(mask)) mask_out(1, 1, 1) = mask
if(present(rmask)) then
    if (rmask < 0.5) mask_out(1, 1, 1) = .false.
endif
if(present(weight)) weight_out(1, 1, 1) = weight

send_data_0d = send_data_3d(diag_field_id, field_out, time, 1, 1, 1, mask_out, weight_out)

end function send_data_0d

!-------------------------------------------------------------------------
! <FUNCTION NAME="send_data_1d" INTERFACE="send_data">
!   <IN NAME="diag_field_id" TYPE="integer"  > </IN>
!   <IN NAME="field" TYPE="real" DIM="(:)" > </IN>
!   <IN NAME="time" TYPE="time_type"  > </IN>
! </FUNCTION>

function send_data_1d(diag_field_id, field, time, is_in, mask, rmask, weight)

logical :: send_data_1d
integer, intent(in) :: diag_field_id
real, intent(in) :: field(:)
type (time_type), intent(in) :: time
integer, optional :: is_in
logical, optional :: mask(:)
real, optional :: rmask(:), weight(:)

real :: field_out(size(field), 1, 1)
real :: weight_out(size(field), 1, 1)
logical :: mask_out(size(field), 1, 1)


! First copy the data to a three d array with last element 1
field_out(:, 1, 1) = field

! Default values for mask and weight
mask_out = .true.
weight_out = 1.0

if(present(mask)) mask_out(:, 1, 1) = mask
if(present(rmask)) where (rmask < 0.5) mask_out(:, 1, 1) = .false.
if(present(weight)) weight_out(:, 1, 1) = weight

send_data_1d = send_data_3d(diag_field_id, field_out, time, is_in, 1, 1, mask_out, weight_out)

end function send_data_1d


!-------------------------------------------------------------------------
! <FUNCTION NAME="send_data_2d" INTERFACE="send_data">
!   <IN NAME="diag_field_id" TYPE="integer"  > </IN>
!   <IN NAME="field" TYPE="real" DIM="(:,:)" > </IN>
!   <IN NAME="time" TYPE="time_type"  > </IN>
! </FUNCTION>

function send_data_2d(diag_field_id, field, time, is_in, js_in, &
    mask, rmask, weight)

logical :: send_data_2d
integer, intent(in) :: diag_field_id
real, intent(in) :: field(:, :)
type (time_type), intent(in) :: time
integer, optional :: is_in, js_in
logical, optional :: mask(:, :)
real, optional :: rmask(:, :), weight(:, :)

real :: field_out(size(field, 1), size(field, 2), 1)
real :: weight_out(size(field, 1), size(field, 2), 1)
logical :: mask_out(size(field, 1), size(field, 2), 1)


! First copy the data to a three d array with last element 1
field_out(:, :, 1) = field

! Default values for mask and weight
mask_out = .true.
weight_out = 1.0

if(present(mask)) mask_out(:, :, 1) = mask
if(present(rmask)) where (rmask < 0.5) mask_out(:, :, 1) = .false.
if(present(weight)) weight_out(:, :, 1) = weight

send_data_2d = send_data_3d(diag_field_id, field_out, time, is_in, js_in, 1, mask_out, weight_out)

end function send_data_2d

!-------------------------------------------------------------------------
! <FUNCTION NAME="send_data_3d" INTERFACE="send_data">
!   <IN NAME="diag_field_id" TYPE="integer"  > </IN>
!   <IN NAME="field" TYPE="real" DIM="(:,:,:)" > </IN>
!   <IN NAME="time" TYPE="time_type"  > </IN>
! </FUNCTION>

function send_data_3d(diag_field_id, field, time, is_in, js_in, ks_in, &
             mask, rmask, weight)

logical :: send_data_3d
integer, intent(in) :: diag_field_id
real, intent(in) :: field(:, :, :)
type (time_type), intent(in) :: time
integer, optional :: is_in, js_in, ks_in
logical, optional :: mask(:, :, :)
real, optional :: rmask(:, :, :)
real, optional :: weight(:, :, :)

logical :: average
integer :: i, out_num, file_num, num_elements, total_elements, num, n1, n2, n3
integer :: sub_size, freq, units, is, js, ks, kount
character(len=128) :: error_string

if (.not.module_is_initialized) call error_mesg ('send_data in diag_manager_mod', &
                              'module has not been initialized', FATAL)

! If is, js, or ks not present default them to 1
is = 1
if(present(is_in)) is = is_in
js = 1
if(present(js_in)) js = js_in
ks = 1
if(present(ks_in)) ks = ks_in

! Whole bunch of things must be done on first call
! DO NEED TO DO THIS ONLY ON FIRST NON-STATIC FIELD
! Originally, distinguished static fields here; don't need to?
if(first_send_data_call .and. .not. input_fields(diag_field_id)%static) then
   call output_setup()
   first_send_data_call = .false.
! Also set the time last flushed for all files to this time
   do i = 1, num_files
      files(i)%last_flush = time
   end do
endif

! If diag_field_id is < 0 it means that this field is unused, return
if(diag_field_id < 0) then
   send_data_3d = .false.
   return
else
   send_data_3d = .true.
endif


! IF NON-STATIC FIELD COMES THROUGH AFTER OUTPUT_SETUP NEED TO WRITE IT?
! COULD ALSO WRITE STATIC FIELDS AT END FOR NOW???

! Loop through each output field that depends on this input field
do i = 1, input_fields(diag_field_id)%num_output_fields

! Get index to an output field
   out_num = input_fields(diag_field_id)%output_fields(i)

! Get index to output file for this field
   file_num = output_fields(out_num)%output_file

! Output frequency for this file is
   freq = files(file_num)%output_freq
   units = files(file_num)%output_units

! Is this output field being time averaged?
   average = output_fields(out_num)%time_average

! Initialize output time for fields output every time step
! this will guarantee that first time here field will not be written
   if (freq == EVERY_TIME) then
     if (output_fields(out_num)%next_output == output_fields(out_num)%last_output) &
         output_fields(out_num)%next_output = time
   endif

! Is it time to output for this field; CAREFUL ABOUT > vs >= HERE
!  if(time > output_fields(out_num)%next_output .and. &
!     .not.output_fields(out_num)%static) then
   if( .not.output_fields(out_num)%static .and. freq /= END_OF_RUN .and. &
        time > output_fields(out_num)%next_output ) then
! A non-static field that has skipped a time level is an error
!del  if(time >= output_fields(out_num)%next_next_output .and. .not. &
!del     output_fields(out_num)%static) then
      if(time >  output_fields(out_num)%next_next_output .and. freq > 0) then
         call error_mesg('send_data', &
            'Attempt to output time level without input data in interval', WARNING)
      end if
! Static fields need to have a dummy time here
!     if (output_fields(out_num)%static) then
!        output_fields(out_num)%next_output = base_time
!        output_fields(out_num)%last_output = base_time
!     end if

! Check to see that total number of elements stored is a multiple of the
! size of the grid; if not, die
! There are clever ways to do checksums on elements here which should
! be considered.
      num_elements = output_fields(out_num)%num_elements
! error if no data is in the buffer
      if ( num_elements == 0 ) then
         write (error_string,'(a,"/",a)')  &
                trim(input_fields(diag_field_id)%module_name), &
                trim(output_fields(out_num)%output_name)
         call error_mesg ('send_data in diag_manager_mod', &
                      'module/output_field '//trim(error_string)//&
                    &', attempt to write buffer with no data', FATAL)
      endif

      total_elements = output_fields(out_num)%total_elements
      if(num_elements / total_elements * total_elements /= num_elements) then
         call error_mesg('send_data', &
            'Attempt to write incompletely filled output buffer', WARNING)
      endif

! If average get size: Average intervals are last_output, next_output
      if(average) then
         num = num_elements / total_elements
!!       if (mpp_pe() == 0) write(*, *) 'number of elements in average is ', num
         if (input_fields(diag_field_id)%missing_value_present) then
            where (output_fields(out_num)%buffer /=            &
                    input_fields(diag_field_id)%missing_value) &
            output_fields(out_num)%buffer = output_fields(out_num)%buffer / num
         else
            output_fields(out_num)%buffer = output_fields(out_num)%buffer / num
         endif
      endif

! Output field
         call diag_data_out(file_num, out_num, &
            output_fields(out_num)%buffer, output_fields(out_num)%next_output)

! Take care of cleaning up the time counters and the storeage size
      output_fields(out_num)%last_output = output_fields(out_num)%next_output
    if (freq == EVERY_TIME) then
      output_fields(out_num)%next_output = time
    else
      output_fields(out_num)%next_output = output_fields(out_num)%next_next_output
      output_fields(out_num)%next_next_output = &
         diag_time_inc(output_fields(out_num)%next_next_output, freq, units)
    endif
      output_fields(out_num)%num_elements = 0
      output_fields(out_num)%buffer = 0.0

   end if

! Finished output of previously buffered data, now deal with buffering new data
! Get number of elements and extent for this subdomain
   n1 = size(field, 1); n2 = size(field, 2); n3 = size(field, 3)
   sub_size = n1 * n2 * n3

! Take care of submitted field data
   if(average) then
      output_fields(out_num)%buffer(is:is+n1-1,js: js+n2-1, ks:ks+n3-1) = &
         output_fields(out_num)%buffer(is:is+n1-1,js: js+n2-1, ks:ks+n3-1) + &
         field
   else
       output_fields(out_num)%buffer(is:is+n1-1,js: js+n2-1, ks:ks+n3-1) = field
   endif

! Increment number of sub_field elements seen
   output_fields(out_num)%num_elements = &
      output_fields(out_num)%num_elements + sub_size

! If mask is present, check to see if a missing value is defined
   if ((present(mask).or.present(rmask)) .and. &
      .not. input_fields(diag_field_id)%missing_value_present) then
      kount = 0
      if (present( mask)) kount = kount + count(.not.mask)
      if (present(rmask)) kount = kount + count(rmask < 0.5)
      if (kount > 0) call error_mesg('send_data in diag_manager_mod', &
           'Mask will be ignored since missing values were not specified', &
            NOTE)
   end if

! If mask (or rmask) and missing value present, then insert missing value
      if (present(mask) .and. input_fields(diag_field_id)%missing_value_present) then
          where (.not.mask) output_fields(out_num)%buffer(is:is+n1-1, js:js+n2-1, ks:ks+n3-1) = &
                             input_fields(diag_field_id)%missing_value
      endif
      if (present(rmask) .and. input_fields(diag_field_id)%missing_value_present) then
          where (rmask < 0.5) output_fields(out_num)%buffer(is:is+n1-1, js:js+n2-1, ks:ks+n3-1) = &
                               input_fields(diag_field_id)%missing_value
      endif

end do

end function send_data_3d

!-------------------------------------------------------------------------
! <SUBROUTINE NAME="diag_manager_end">

!   <OVERVIEW>
!     Exit Diagnostics Manager.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Flushes diagnostic buffers where necessary. Close diagnostics files.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call diag_manager_end (time)
!   </TEMPLATE>
!   <IN NAME="TIME" TYPE="time_type"></IN>

subroutine diag_manager_end (time)

! PROBABLY JUST WRITE STATIC FIELDS HERE FOR NOW???

type(time_type), intent(in) :: time

integer :: i, file_num, freq, num_elements, total_elements, num, j, input_num
integer :: out_num
logical :: average
character(len=128) :: message

if(.not.module_is_initialized) return

do i = 1, num_output_fields

! only output non-static fields here
   input_num = output_fields(i)%input_field
   if (input_fields(input_num)%static) cycle

! skip fields that were not registered
   if (.not.input_fields(input_num)%register) cycle

! Get index to output file for this field
   file_num = output_fields(i)%output_file

! Output frequency for this file is
   freq = files(file_num)%output_freq

! Is it time to output for this field; CAREFUL ABOUT >= vs > HERE
! For end should be >= because no more data is coming 
!  if(time >= output_fields(i)%next_output) then
! Assume buffer is full for every time step output?
   if(time >= output_fields(i)%next_output .or. freq == END_OF_RUN) then

      if(time >= output_fields(i)%next_next_output .and. freq > 0) then
         write (message,'(a,"/",a)') trim(input_fields(input_num)%module_name), &
                                     trim(output_fields(i)%output_name)
         call error_mesg('diag_manager_end', 'module/output_field ' //  &
                        &trim(message)//', attempt to output time level &
                        &without input data in interval', WARNING)
      end if

! Check to see that total number of elements stored is a multiple of the
! size of the grid; if not, die
      num_elements = output_fields(i)%num_elements
! error if no data is in the buffer
      if ( num_elements == 0 ) then
         write (message,'(a,"/",a)')  &
                trim(input_fields(input_num)%module_name), &
                trim(output_fields(i)%output_name)
         call error_mesg ('diag_manager_end in diag_manager_mod', &
                      'module/output_field '//trim(message)//&
                    &', attempt to write buffer with no data', FATAL)
      endif

      total_elements = output_fields(i)%total_elements
      if(num_elements / total_elements * total_elements /= num_elements) then
         call error_mesg('diag_manager_end', &
            'Attempt to write incompletely filled output buffer', WARNING)
      endif

! If average, get size: Need time average intervals here
      if(output_fields(i)%time_average) then
         num = num_elements / total_elements
         if (input_fields(input_num)%missing_value_present) then
            where (output_fields(i)%buffer /=        &
                    input_fields(input_num)%missing_value) &
            output_fields(i)%buffer = output_fields(i)%buffer / num
         else
            output_fields(i)%buffer = output_fields(i)%buffer / num
         endif
      endif

! Output field
      if (freq == END_OF_RUN) output_fields(i)%next_output = time
      call diag_data_out(file_num, i, &
            output_fields(i)%buffer, output_fields(i)%next_output)

   end if
 
end do

! If there are no non-static fields then need to setup output data
if(first_send_data_call) then
   call output_setup()
   first_send_data_call = .false.
endif

! Output static fields and signal end of output for each file
do i = 1, num_files
! Loop to look for static fields
   do j = 1, files(i)%num_fields
!bug  input_num = output_fields(j)%input_field
!bug  if(input_fields(input_num)%static) &
!bug     call diag_data_out(i, j, output_fields(j)%buffer, time)
      out_num = files(i)%fields(j)
      input_num = output_fields(out_num)%input_field

! skip fields that were not registered
      if (.not.input_fields(input_num)%register) cycle

! only output static fields here
      if (.not.output_fields(out_num)%static) cycle

! error if no data is in the buffer
      if ( output_fields(out_num)%num_elements == 0 ) then
         write (message,'(a,"/",a)')  &
                trim(input_fields(input_num)%module_name), &
                trim(output_fields(out_num)%output_name)
         call error_mesg ('diag_manager_end in diag_manager_mod', &
                      'module/output_field '//trim(message)//&
                    &', attempt to write buffer with no data', FATAL)
      endif

      call diag_data_out(i, out_num, output_fields(out_num)%buffer, time)
   end do

! Close up this file   
   call diag_output_end(files(i)%file_unit)
end do

module_is_initialized = .FALSE.

end subroutine diag_manager_end
! </SUBROUTINE>

!-------------------------------------------------------------------------

subroutine init_file(name, output_freq, output_units, format, time_units, long_name)

character(len=*), intent(in) :: name, long_name
integer, intent(in) :: output_freq, output_units, format, time_units

! Get a number for this file
num_files = num_files + 1
if(num_files > max_files) then
   call error_mesg('init_file', 'max_files for diagnostics exceeded', FATAL)
endif

files(num_files)%name = trim(name)
files(num_files)%output_freq = output_freq
files(num_files)%output_units = output_units
files(num_files)%format = format
files(num_files)%time_units = time_units
files(num_files)%long_name = trim(long_name)
files(num_files)%num_fields = 0
! This value should be updated with first call to send_data for non-static field
files(num_files)%last_flush = set_time(0, 0)

end subroutine init_file

!--------------------------------------------------------------------------

subroutine init_input_field(module_name, field_name)

character(len=*), intent(in) :: module_name, field_name

! Get a number for this input_field if not already set up
if(find_input_field(module_name, field_name) < 0) then
   num_input_fields = num_input_fields + 1
   if(num_input_fields > max_input_fields) then
      call error_mesg('init_input_field', 'max_input_fields exceeded', FATAL)
   end if
else
! If this is already initialized don't need to do anything
   return
end if

input_fields(num_input_fields)%module_name = trim(module_name)
input_fields(num_input_fields)%field_name = trim(field_name)
input_fields(num_input_fields)%num_output_fields = 0

! Set flag that this field has not been registered
input_fields(num_input_fields)%register = .false.

end subroutine init_input_field

!---------------------------------------------------------------------------

subroutine init_output_field(module_name, field_name, output_name, output_file,&
   time_average, pack)

character(len=*), intent(in) :: module_name, field_name, output_name, output_file
logical, intent(in) :: time_average
integer, intent(in) :: pack

integer :: out_num, in_num, file_num, num_fields

! Get a number for this output field
num_output_fields = num_output_fields + 1
if(num_output_fields > max_output_fields) then
   call error_mesg('init_output_field', 'max_output_fields exceeded', FATAL)
endif
out_num = num_output_fields

write(*,*)'TJH init_output_field :<',module_name,'>'
write(*,*)'TJH init_output_field :<',field_name,'>'

! First, find the index to the associated input field
in_num = find_input_field(module_name, field_name)

write(*,*)'TJH init_output_field : in_num = ',in_num

if(in_num < 0) then
   call error_mesg('init_output_field', &
      'module_name / field_name pair not registered as an input field', FATAL)
endif

! Add this output field into the list for this input field
input_fields(in_num)%num_output_fields = &
   input_fields(in_num)%num_output_fields + 1
if(input_fields(in_num)%num_output_fields > max_out_per_in_field) then
   call error_mesg('init_output_field', 'max_out_per_in_field exceeded', FATAL)
endif
input_fields(in_num)%output_fields(input_fields(in_num)%num_output_fields) &
   = out_num

! Also put pointer to input field in this output field
output_fields(out_num)%input_field = in_num



! Next, find the number for the corresponding file
file_num = find_file(output_file)
if(file_num < 0) then
   call error_mesg('init_output_field', 'file has not been registered', FATAL)
end if

! Insert this field into list for this file
files(file_num)%num_fields = files(file_num)%num_fields + 1
if(files(file_num)%num_fields > max_fields_per_file) then
   call error_mesg('init_output_field', 'max_fields_per_file exceeded', FATAL)
endif
num_fields = files(file_num)%num_fields
files(file_num)%fields(num_fields) = out_num

! Set the file for this output field
output_fields(out_num)%output_file = file_num

! Enter the other data for this output field
output_fields(out_num)%output_name = trim(output_name)
output_fields(out_num)%pack = pack

! cannot time average fields output every time
if (files(file_num)%output_freq == EVERY_TIME) then
  output_fields(out_num)%time_average = .false.
else
  output_fields(out_num)%time_average = time_average
endif


end subroutine init_output_field

!-------------------------------------------------------------------------

! <SUBROUTINE NAME="diag_manager_init">

!   <OVERVIEW>
!     Initialize Diagnostics Manager.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Open and read diag_table. Select fields and files for diagnostic output.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call diag_manager_init()
!   </TEMPLATE>

subroutine diag_manager_init()

type tableB_type
   character(len=128) :: module_name,field_name,output_name,name
   character(len=50) :: time_sampling
   logical :: time_avg
   character(len=50) :: other_ops
   integer :: pack
end type tableB_type

type tableA_type
   character(len=128) :: name
   integer :: output_freq
   character(len=10) :: output_freq_units
   integer :: format
   character(len=10) :: time_units
   character(len=128) :: long_name
end type tableA_type

character(len=256) :: record
character(len=9)   :: amonth

integer :: iunit,n,m,num_fields,time_units, output_freq_units, nfiles,nfields
integer :: j, log_unit, name_len, nrecs

type(tableB_type) :: textB
type(tableA_type) :: textA

!  If the module was already initialized do nothing
if (module_is_initialized) return

if (.not.file_exist('diag_table') ) &
call error_mesg('diag_manager_init','file diag_table nonexistent',FATAL)

call mpp_open(iunit, 'diag_table', action=MPP_RDONLY)

! Read in the global file labeling string
read(iunit, *, end = 99, err=99) global_descriptor

! Read in the base date
read(iunit, *, end = 99, err = 99) base_year, base_month, base_day, &
   base_hour, base_minute, base_second

! Set up the time type for base time
if (get_calendar_type() /= NO_CALENDAR) then
   base_time = set_date(base_year, base_month, base_day, base_hour, &
                        base_minute, base_second)
   amonth = month_name(base_month)
else
! No calendar - ignore year and month
   base_time = set_time(base_hour*3600+base_minute*60+base_second, base_day)
   base_year  = 0
   base_month = 0
   amonth = 'day'
end if

nrecs=0
nfiles=0
do while (nfiles <= max_files)
   read(iunit,'(a)',end=86,err=85) record
   nrecs=nrecs+1
   if (record(1:1) == '#') cycle
   read(record,*,err=85,end=85) textA
   ! test file format to make sure its OK
   if (textA%format .gt. 2 .or. textA%format .lt. 1) cycle
   nfiles=nfiles+1
   time_units = 0
   output_freq_units = 0
   do j = 1, size(time_unit_list)
      if(textA%time_units == time_unit_list(j)) time_units = j
      if(textA%output_freq_units == time_unit_list(j)) output_freq_units = j
   end do
   if(time_units == 0) &
        call error_mesg('diag_manager_init','invalid time units',FATAL)
   if(output_freq_units == 0) & 
        call error_mesg('diag_manager_init','invalid output frequency units',FATAL)
   ! remove trailing .nc extension from file name 
   name_len = len_trim(textA%name)
   if (textA%name(name_len-2:name_len) == '.nc') textA%name = textA%name(1:name_len-3)
   ! assign values to file_types
   call init_file(textA%name,textA%output_freq, output_freq_units, &
        textA%format, time_units,textA%long_name)
85 continue
enddo
call error_mesg('diag_manager_init','too many files in table', FATAL)
86 continue



rewind(iunit)
!if (nfiles .lt. 1) call error_mesg('diag_manager_init','error reading file records',FATAL)

nfields=0;nrecs=0
do while (nfields <= max_output_fields)
   read(iunit,'(a)',end=94,err=93) record
   nrecs=nrecs+1
   if (record(1:1) == '#') cycle
   read(record,*,end=93,err=93) textB
   if (textB%pack .gt. 8 .or. textB%pack .lt. 1) cycle
   nfields=nfields+1
   !   assign values to field_types
   call init_input_field(textB%module_name,textB%field_name)
   !   remove trailing .nc extension
   name_len= len_trim(textB%name)
   if (textB%name(name_len-2:name_len) == '.nc') &
        textB%name = textB%name(1:name_len-3)

   write(*,*)'TJH diag_manager_init:textB%module_name<',trim(textB%module_name),'>'
   write(*,*)'TJH diag_manager_init:textB%field_name <',trim(textB%field_name),'>'
   write(*,*)'TJH diag_manager_init:textB%output_name<',trim(textB%output_name),'>'
   write(*,*)'TJH diag_manager_init:textB%name       <',trim(textB%name),'>' 

   call init_output_field(textB%module_name,textB%field_name,textB%output_name,&
        textB%name,textB%time_avg,textB%pack)
93 continue
enddo
call error_mesg('diag_manager_init','too many fields in table', FATAL)
94 continue




!if (nfields .lt. 1) call error_mesg('diag_manager_init','error reading field records',FATAL)

call close_file(iunit)

! version number to logfile
  call write_version_number (version, tagname)

  log_unit = stdlog()
  if ( mpp_pe() == 0 ) then
       write (log_unit,95) base_year, trim(amonth), base_day, &
                           base_hour, base_minute, base_second
  endif
  call close_file (log_unit)
95 format ('base date used = ',i4,1x,a,2i3,2(':',i2.2),' gmt')


module_is_initialized = .true.

return

99 continue

call error_mesg('diag_manager_init','error reading table',FATAL)

! AT THIS POINT NEED TO PUT IN THE GLOBAL BASE DATE AND DESCRIPTOR TOO


end subroutine diag_manager_init
! </SUBROUTINE>

!-------------------------------------------------------------------------

function find_file(name)

integer :: find_file
character(len=*), intent(in) :: name

integer :: i

find_file = -1
do i = 1, num_files
   if(files(i)%name == name) then
      find_file = i
      return
   end if
end do

end function find_file

!-------------------------------------------------------------------------

function find_input_field(module_name, field_name)

integer find_input_field
character(len=*), intent(in) :: module_name, field_name

integer :: i

find_input_field = -1
do i = 1, num_input_fields
   if(trim(input_fields(i)%module_name) == trim(module_name) .and. &
      lcase(trim(input_fields(i)%field_name)) == &
      lcase(trim(field_name))) then 
      find_input_field = i
      return
   endif
end do

end function find_input_field

!-------------------------------------------------------------------------

subroutine output_setup()

! WARNING: Assumes that all data structures are fully initialized

integer i, j, field_num, n_fields, axes(4), input_field_num, num_axes
character(len=128) time_units, avg
logical :: file_time_avg
integer :: time_axis_id(1)


! Set up for output, focused on netcdf

! First, get a file_unit and a time axis identifier for each file
do i = 1, num_files
   file_time_avg = .false.

! Skip this file if no fields are to be output
   if (files(i)%num_fields == 0) cycle

! Here is where time_units string must be set up; time since base date
!  write(*, *) 'file ', i, trim(global_descriptor)
!  write(*, *) 'file ', i, 'time units label is '
   write(time_units, 11) trim(time_unit_list(files(i)%time_units)), base_year, &
      base_month, base_day, base_hour, base_minute, base_second
!  write(*, *) time_units
 11 format(a, ' since ', i4.4, '-', i2.2, '-', i2.2, ' ', i2.2, ':', i2.2, ':', i2.2) 
   call diag_output_init(files(i)%name, files(i)%format, global_descriptor, &
      files(i)%long_name, time_units, files(i)%file_unit, files(i)%time_axis_id)

! Loop through all fields with this file to output axes
   do j = 1, files(i)%num_fields
      field_num = files(i)%fields(j)
      input_field_num = output_fields(field_num)%input_field
      if (.not.input_fields(input_field_num)%register) cycle
! Put the time axis in the axis field
      num_axes = output_fields(field_num)%num_axes
      axes(1:num_axes) = output_fields(field_num)%axes(1:num_axes)
      axes(num_axes + 1) = files(i)%time_axis_id
      call write_axis_meta_data(files(i)%file_unit, axes(1:num_axes + 1))
   end do

! Now output metadata for each of the fields associated with this file
   do j = 1, files(i)%num_fields
      field_num = files(i)%fields(j)
      input_field_num = output_fields(field_num)%input_field
      if (.not.input_fields(input_field_num)%register) cycle
! Put the time axis in the axis field
      num_axes = output_fields(field_num)%num_axes
      axes(1:num_axes) = output_fields(field_num)%axes(1:num_axes)
      if (.not. output_fields(field_num)%static) then
        num_axes=num_axes+1
        axes(num_axes) = files(i)%time_axis_id
      endif
      if(output_fields(field_num)%time_average) then
         avg = avg_name
         file_time_avg = .true.
      else
         avg = ""
      end if
      if(input_fields(input_field_num)%missing_value_present) then
         output_fields(field_num)%f_type = write_field_meta_data(files(i)%file_unit, &
            output_fields(field_num)%output_name, axes(1:num_axes), &
            input_fields(input_field_num)%units, &
            input_fields(input_field_num)%long_name, &
            input_fields(input_field_num)%range, output_fields(field_num)%pack,&
            input_fields(input_field_num)%missing_value, avg_name = avg)

! NEED TO TAKE CARE OF TIME AVERAGING INFO TOO BOTH CASES
      else
         output_fields(field_num)%f_type = write_field_meta_data(files(i)%file_unit, &
            output_fields(field_num)%output_name, axes(1:num_axes), &
            input_fields(input_field_num)%units, &
            input_fields(input_field_num)%long_name, &
            input_fields(input_field_num)%range, output_fields(field_num)%pack,&
               avg_name = avg)
      endif
   end do

! If any of the fields in the file are time averaged, need to output the axes
! Use double precision since time axis is double precision
   if(file_time_avg) then
      time_axis_id(1) = files(i)%time_axis_id
      files(i)%f_avg_start = write_field_meta_data(files(i)%file_unit, &
         avg_name // '_T1', time_axis_id, time_units, &
         "Start time for average period", pack=1)
      files(i)%f_avg_end = write_field_meta_data(files(i)%file_unit, &
         avg_name // '_T2', time_axis_id, time_units, &
         "End time for average period", pack=1)
      files(i)%f_avg_nitems = write_field_meta_data(files(i)%file_unit, &
         avg_name // '_DT', time_axis_id,     &
         trim(time_unit_list(files(i)%time_units)), &
         "Length of average period", pack=1)
   end if

! Let lower levels know that all meta data has been sent
   call done_meta_data(files(i)%file_unit)

end do

! ALSO NEED TO SEND DATA FOR ANY STATIC FIELD THAT HAS ALREADY BEEN WRITTEN

end subroutine output_setup

!-------------------------------------------------------------------------

subroutine diag_data_out(file, field, dat, time)

integer, intent(in) :: file, field
real, intent(inout) :: dat(:, :, :)
type(time_type), intent(in) :: time
integer :: i, num
real :: dif, time_data(1, 1, 1), dt_time(1, 1, 1), start_dif

! Check to see if this field was registered (may not need this here)
num = output_fields(field)%input_field
if (.not.input_fields(num)%register) then
  if (mpp_pe() == 0) then
    call error_mesg ('diag_data_out in diag_manager_mod', &
                     'module '//trim(input_fields(num)%module_name)//&
                    &', field '//trim(input_fields(num)%field_name)//&
                    &' not registered', NOTE)
  endif
  return
endif

dif = get_date_dif(time, base_time, files(file)%time_units)

call diag_field_out(files(file)%file_unit,output_fields(field)%f_type, dat, dif)

! *** inserted this line because start_dif < 0 for static fields ***
if (output_fields(field)%static) return

start_dif = get_date_dif(output_fields(field)%last_output, base_time, &
   files(file)%time_units)

! Need to write average axes out; only do this when first averaged axis for
! this file is written for now. Might want to reconsider when to do this???
! Is this field the first averaged one in its file?
do i = 1, files(file)%num_fields
   num = files(file)%fields(i)
   if(output_fields(num)%time_average .and. .not.output_fields(num)%static .and. &
      input_fields(output_fields(num)%input_field)%register) then
      if(num == field) then
! Output the axes if this is first time-averaged field
         time_data(1, 1, 1) = start_dif
         call diag_field_out(files(file)%file_unit, files(file)%f_avg_start, &
            time_data, dif)
         time_data(1, 1, 1) = dif
         call diag_field_out(files(file)%file_unit, files(file)%f_avg_end, &
            time_data, dif)
! Compute the number of items in this averaging period
!        dt_time(1, 1, 1) = output_fields(field)%num_elements / &
!           output_fields(field)%total_elements
! Compute the length of the average
         dt_time(1, 1, 1) = dif - start_dif
         call diag_field_out(files(file)%file_unit, files(file)%f_avg_nitems, &
            dt_time, dif)
      endif
      goto 10
   end if
end do

! If write time is greater than last_flush for this file, flush it
10 if(time > files(file)%last_flush) then
   call diag_flush(files(file)%file_unit)
   files(file)%last_flush = time
endif

end subroutine diag_data_out

!-------------------------------------------------------------------------

function get_date_dif(t2, t1, units)

real :: get_date_dif
type(time_type), intent(in) :: t2, t1
integer, intent(in) :: units

integer :: year, month, day, hour, minute, second, dif_seconds, dif_days
type(time_type) :: dif_time

! Compute time axis label value
if(t2 < t1)   call error_mesg('get_date_dif', &
                't2 is less than t1', FATAL)

dif_time = t2 - t1
!del call get_date(dif_time, year, month, day, hour, minute, second) ! not used
call get_time(dif_time, dif_seconds, dif_days)

if(units == DIAG_SECONDS) then
   get_date_dif = dif_seconds + 86400 * dif_days
else if(units == DIAG_MINUTES) then
   get_date_dif = 1440 * dif_days + dif_seconds / 60.
else if(units == DIAG_HOURS) then
   get_date_dif = 24 * dif_days + dif_seconds / 3600.
else if(units == DIAG_DAYS) then
   get_date_dif = dif_days + dif_seconds / 86400.
else if(units == DIAG_MONTHS) then
   call error_mesg('diag_data_out', 'months not supported as output units', FATAL)
else if(units == DIAG_YEARS) then
call error_mesg('diag_data_out', 'years not supported as output units', FATAL)
else
   call error_mesg('diag_data_out', 'illegal time units', FATAL)
end if

end function get_date_dif

!-------------------------------------------------------------------------

subroutine print_file(i)

integer, intent(in) :: i
integer :: j

write(*, *)'----------------------------------------------------'
write(*, *) 'contents of file entry ', i
write(*, *) 'name'
write(*, *) trim(files(i)%name)
write(*, *) 'output_freq ', files(i)%output_freq
write(*, *) 'format ', files(i)%format
write(*, *) 'time units '
write(*, *) time_unit_list(files(i)%time_units)
write(*, *) 'long_name '
write(*, *) trim(files(i)%long_name)
write(*, *) 'num fields ', files(i)%num_fields
do j = 1, files(i)%num_fields
   write(*, *) 'field ', j, ' is ', files(i)%fields(j)
end do
write(*, *)'----------------------------------------------------'

end subroutine print_file

!------------------------------------------------------------------------

subroutine print_input_field(i)

integer, intent(in) :: i
integer :: j

write(*, *)'----------------------------------------------------'
write(*, *) 'contents of field entry ', i
write(*, *) 'module and field name'
write(*, *) trim(input_fields(i)%module_name)
write(*, *) trim(input_fields(i)%field_name)
write(*, *) 'long_name '
write(*, *) trim(input_fields(i)%long_name)
write(*, *) 'units '
write(*, *) trim(input_fields(i)%units)
write(*, *) 'axes ', input_fields(i)%axes(1), input_fields(i)%axes(2), &
   input_fields(i)%axes(3)
write(*, *) 'num_axes ', input_fields(i)%num_axes
write(*, *) 'missing value and range present ', &
   input_fields(i)%missing_value_present, input_fields(i)%range_present
write(*, *) 'missing value and range ', input_fields(i)%missing_value, &
   input_fields(i)%range(1), input_fields(i)%range(2)
write(*, *) 'sizes and total ', input_fields(i)%size(1), &
   input_fields(i)%size(2), input_fields(i)%size(3), &
   input_fields(i)%total_elements
write(*, *) 'static ', input_fields(i)%static
write(*, *) 'num_output_fields ', input_fields(i)%num_output_fields
do j = 1, input_fields(i)%num_output_fields
   write(*, *) 'output field ', j, ' is ', input_fields(i)%output_fields(j)
end do
write(*, *)'----------------------------------------------------'

end subroutine print_input_field

!-------------------------------------------------------------------------

subroutine print_output_field(i)

integer, intent(in) :: i
integer :: seconds, days

write(*, *)'----------------------------------------------------'
write(*, *) 'contents of output field ', i
write(*, *) 'input field and output file ', output_fields(i)%input_field, &
   output_fields(i)%output_file
write(*, *) 'output name ', trim(output_fields(i)%output_name)
write(*, *) 'time average ', output_fields(i)%time_average
write(*, *) 'pack ', output_fields(i)%pack
write(*, *) 'total_elements and num_elements ', &
   output_fields(i)%total_elements, output_fields(i)%num_elements
write(*, *) 'static ', output_fields(i)%static
write(*, *) 'num_axes ', output_fields(i)%num_axes
write(*, *) 'axes ', output_fields(i)%axes
call get_time(output_fields(i)%last_output, seconds, days)
write(*, *) 'last output ', seconds, ' seconds ', days, ' days '
call get_time(output_fields(i)%next_output, seconds, days)
write(*, *) 'next output ', seconds, ' seconds ', days, ' days '
call get_time(output_fields(i)%next_next_output, seconds, days)
write(*, *) 'next_next_output ', seconds, ' seconds ', days, ' days '
write(*, *)'----------------------------------------------------'

end subroutine print_output_field

!-------------------------------------------------------------------------

function diag_time_inc(time, output_freq, output_units)

type (time_type) :: diag_time_inc
type (time_type), intent(in) :: time
integer, intent(in) :: output_freq, output_units

! special values for output frequency are -1 for output at end of run
! and 0 for every timestep.  Need to check for these here?
! Return zero time increment, hopefully this value is never used

if (output_freq == END_OF_RUN .or. output_freq == EVERY_TIME) then
    diag_time_inc = time
    return
endif

! Make sure calendar was not set after initialization
if (base_year == 0 .and. get_calendar_type() /= NO_CALENDAR) &
  call error_mesg('diag_time_inc',  &
        'calendar_type was set after diag_manager_init', FATAL)
   

if(output_units == DIAG_SECONDS) then
   if (get_calendar_type() == NO_CALENDAR) then
      diag_time_inc = increment_time(time, output_freq, 0)
   else
      diag_time_inc = increment_date(time, 0, 0, 0, 0, 0, output_freq)
   endif
else if(output_units == DIAG_MINUTES) then
   if (get_calendar_type() == NO_CALENDAR) then
      diag_time_inc = increment_time(time, output_freq*60, 0)
   else
      diag_time_inc = increment_date(time, 0, 0, 0, 0, output_freq, 0)
   endif
else if(output_units == DIAG_HOURS) then
   if (get_calendar_type() == NO_CALENDAR) then
      diag_time_inc = increment_time(time, output_freq*3600, 0)
   else
      diag_time_inc = increment_date(time, 0, 0, 0, output_freq, 0, 0)
   endif
else if(output_units == DIAG_DAYS) then
   if (get_calendar_type() == NO_CALENDAR) then
      diag_time_inc = increment_time(time, 0, output_freq)
   else
      diag_time_inc = increment_date(time, 0, 0, output_freq, 0, 0, 0)
   endif
else if(output_units == DIAG_MONTHS) then
   if (get_calendar_type() == NO_CALENDAR) then
      call error_mesg('diag_time_inc', &
         'output units of months not allowed with no calendar', FATAL)
   else
      diag_time_inc = increment_date(time, 0, output_freq, 0, 0, 0, 0)
   endif
else if(output_units == DIAG_YEARS) then
   if (get_calendar_type() == NO_CALENDAR) then
      call error_mesg('diag_time_inc', &
         'output units of years not allowed with no calendar', FATAL)
   else
      diag_time_inc = increment_date(time, output_freq, 0, 0, 0, 0, 0)
   endif
else 
   call error_mesg('diag_time_inc','illegal output units',FATAL)
endif

end function diag_time_inc

!-------------------------------------------------------------------------

! <FUNCTION NAME="get_base_time">

!   <OVERVIEW>
!     Return base time for diagnostics. 
!   </OVERVIEW>
!   <DESCRIPTION>
!     Return base time for diagnostics (note: base time must be >= model time).
!   </DESCRIPTION>
!   <TEMPLATE>
!     call get_base_time ()
!   </TEMPLATE>

 function get_base_time ()
 type(time_type) :: get_base_time

   if (.not.module_is_initialized) call error_mesg (  &
                        'get_base_time in diag_manager_mod', &
                        'module has not been initialized', FATAL)

   get_base_time = base_time

 end function get_base_time
! </FUNCTION>

!-------------------------------------------------------------------------

! <SUBROUTINE NAME="get_base_date">

!   <OVERVIEW>
!     Return base date for diagnostics.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Return date information for diagnostic reference time.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call get_base_date (year, month, day, hour, minute, second)
!   </TEMPLATE>

 subroutine get_base_date (year, month, day, hour, minute, second)
   integer, intent(out) :: year, month, day, hour, minute, second

   if (.not.module_is_initialized) call error_mesg (  &
                        'get_base_date in diag_manager_mod', &
                        'module has not been initialized', FATAL)

   year   = base_year
   month  = base_month
   day    = base_day
   hour   = base_hour
   minute = base_minute
   second = base_second

 end subroutine get_base_date
! </SUBROUTINE>

!-------------------------------------------------------------------------

 function lcase (cs) 
!
    integer, parameter :: co=iachar('a')-iachar('A') ! case offset
    
    character(len=*), intent(in) :: cs ! character string 
    character(len=len(cs)) :: lcase 
    character :: ca(len(cs)) ! character array
    
    ca=transfer(cs,"x",len(cs)) 
    where (ca >= "A" .and. ca <= "Z") ca=achar(iachar(ca)+co) 
    lcase=transfer(ca,cs) 
    return
    
 end function lcase 

!-------------------------------------------------------------------------

! <FUNCTION NAME="need_data">

!   <OVERVIEW>
!     Determine whether data is needed for the current model time step.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Determine whether data is needed for the current model time step.
!     Since diagnostic data are buffered, the "next" model time is passed
!     instead of the current model time. This call can be used to minimize
!     overhead for complicated diagnostics.
!   </DESCRIPTION>
!   <TEMPLATE>
!     need_data(diag_field_id,next_model_time)
!   </TEMPLATE>

!   <IN NAME="inext_model_time" TYPE="time_type"  >
!     next_model_time = current model time + model time_step
!   </IN>
!   <IN NAME="diag_field_id" TYPE="integer"  > </IN>

function need_data(diag_field_id,next_model_time)
!
! next_model_time = current model time + model time_step
!
type (time_type), intent(in) :: next_model_time
integer, intent(in) :: diag_field_id
logical :: need_data
integer :: i, out_num
 
! loop through output fields

need_data=.false.

! If diag_field_id is < 0 it means that this field is unused, return
if (diag_field_id < 0 ) return

do i= 1,input_fields(diag_field_id)%num_output_fields
! Get index to an output field
   out_num = input_fields(diag_field_id)%output_fields(i)
   if (.not.output_fields(out_num)%static) then
      if (next_model_time > output_fields(out_num)%next_output) need_data=.true.
! Is this output field being time averaged?
! assume average data based on every timestep
! needs to be changed when different forms of averaging are implemented 
      if (output_fields(out_num)%time_average) need_data=.true. 
   endif
enddo

return

end function need_data
! </FUNCTION>

!-------------------------------------------------------------------------
 
end module diag_manager_mod

! <INFO>

!   <COMPILER NAME="COMPILING AND LINKING SOURCE">
!     Any module or program unit using <TT>diag_manager_mod</TT> must contain the line

!   <PRE>
!   use diag_manager_mod
!   </PRE>

!   If netCDF output is desired, the cpp flag <TT>-Duse_netCDF</TT>
!   must be turned on. The loader step requires an explicit link to the
!   netCDF library (typically something like <TT>-L/usr/local/lib
!   -lnetcdf</TT>, depending on the path to the netCDF library).
!   <LINK SRC="http://www.unidata.ucar.edu/packages/netcdf/guidef">netCDF
!   release 3 for fortran</LINK> is required.
!   </COMPILER>
!   <PRECOMP FLAG="PORTABILITY"> 
!     <TT>diag_manager_mod</TT> uses standard f90.
!   </PRECOMP>
!   <LOADER FLAG="ACQUIRING SOURCE">
!     GFDL users can checkout diag_manager_mod using the cvs command
!     <TT>setenv CVSROOT '/home/fms/cvs';cvs co diag_manager</TT>.  
!   </LOADER>

! </INFO>
