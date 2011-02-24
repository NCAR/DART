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
module field_manager_mod

! <CONTACT EMAIL="wfc@gfdl.noaa.gov">
! William Cooke
! </CONTACT>

! <CONTACT EMAIL="mh2@gfdl.noaa.gov">
! Matthew Harrison
! </CONTACT>

! <HISTORY
!  SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/shared/field_manager/field_manager.F90"/>

! <OVERVIEW>
! The field manager reads entries from a field table and stores this information along with the type 
! of field it belongs to. This allows the component models to query the field manager to see if 
! non-default methods of operation are desired. In essence the field table is a powerful type of namelist.
! Default values can be provided for all the fields through a namelist, individual fields can be modified 
! through the field table however.
!</OVERVIEW>

! <DESCRIPTION>
! 
! The field table consists of entries in the following format.
!
! The first line of an entry should consist of three quoted strings.
! The first quoted string will tell the field manager what type of field it is.
! At present the supported types of fields are 
! "tracer" for tracers, 
! "xland_mix" for cross-land mixing, and,
! "checkerboard" for checkerboard null mode.
! 
! The second quoted string will tell the field manager which model the field is 
! being applied to.
! The supported types at present are
! "atmos_mod" for the atmosphere model,
! "ocean_mod" for the ocean model,
! "land_mod" for the land model, and,
! "ice_mod" for the ice model.
! 
! The third quoted string should be a unique name that can be used as a query.
!
! The second and following lines of each entry are called methods in this context.
! Methods can be developed within any module and these modules can query the field manager to
! find any methods that are supplied in the field table.
!
! These lines can consist of two or three quoted strings. The first string will be an identifier 
! that the querying module will ask for. The second string will be a name that the querying module 
! can use to set up values for the module. The third string, if present, can supply parameters to 
! the calling module that can be parsed and used to further modify values.
! 
! An entry is ended with a backslash (/) as the final character in a row.
! 
! Comments can be inserted in the field table by having a # as the first character in the line.
!
! An example of a field table entry could be
! <PRE>
!"tracer","atmos_mod","sphum"/
!"tracer","atmos_mod","sf6"
!"longname","sulf_hex"
!"advection_scheme_horiz","2nd_order"
!"Profile_type","Fixed","surface_value = 0.0E+00"/
! </PRE>
! 
! In this example we have two field entries. 
! 
! The first is a simple declaration of a tracer called "sphum". 
!
! The second is for a tracer called "sf6". Methods that are being applied to this tracer include
! initiating the long name of the tracer to be "sulf_hex", changing the horizontal advection scheme 
! to be second order, and initiating the tracer with a profile with fixed values, in this example all zero.
!
! </DESCRIPTION>

use types_mod, only : r8
use    fms_mod, only : lowercase, stdlog,   &
                       file_exist,  &
                       write_version_number, NOTE, FATAL, error_mesg
use utilities_mod, only : open_file

implicit none
private

! <DATA NAME="NUM_MODELS" TYPE="integer, parameter" DEFAULT="5">
!   Number of models.
! </DATA>

! <DATA NAME="module_is_initialized" TYPE="logical" DEFAULT=".false.">
!   Field manager is initialized.
! </DATA>

! <DATA NAME="MODEL_ATMOS" TYPE="integer, parameter" DEFAULT="1">
!   Atmospheric model.
! </DATA>! 

! <DATA NAME="MODEL_OCEAN" TYPE="integer, parameter" DEFAULT="2">
!   Ocean model.
! </DATA>

! <DATA NAME="MODEL_LAND" TYPE="integer, parameter" DEFAULT="3">
!   Land model.
! </DATA>

! <DATA NAME="MODEL_ICE" TYPE="integer, parameter" DEFAULT="4">
!   Ice model.
! </DATA>

character(len=128) :: version = '$Revision$'
character(len=128) :: tagname = '$Id$'

integer, private :: num_fields = 0
integer, parameter, public :: NUM_MODELS = 5
integer, parameter         :: MAX_FIELDS = 150
integer, parameter         :: MAX_FIELD_METHODS = 100
logical, public :: module_is_initialized  = .false.
integer, parameter, public  :: MODEL_ATMOS = 1, MODEL_OCEAN = 2, MODEL_LAND = 3, MODEL_ICE = 4,&
                               MODEL_DEFAULT = 5

! <TYPE NAME="method_type">
! <DESCRIPTION>
! This method_type is a way to allow a component module to alter the parameters it needs for various tracers.
! In essence this is a way to modify a namelist. A namelist can supply default parameters for all tracers. This 
! method will allow the user to modify these default parameters for an individual tracer. An example could be that 
! the user wishes to use second order advection on a tracer and also use fourth order advection on a second tracer 
! within the same model run. The default advection could be second order and the field table would then indicate 
! that the second tracer requires fourth order advection. This would be parsed by the advection routine.
!
! </DESCRIPTION>

type, public :: method_type

  ! <DATA NAME="method_type :: method_type" TYPE="character" DIM="(128)">
  !   This string represents a tag that a module using this method can key on.
  !   Typically this should contain some reference to the module that is calling it.
  ! </DATA> ! 

  ! <DATA NAME="method_type :: method_name" TYPE="character" DIM="(128)">
  !   This is the name of a method which the module can parse and use to assign 
  !   different default values to a field method.
  ! </DATA> 

  ! <DATA NAME="method_type :: method_control" TYPE="character" DIM="(256)">
  !   This is the string containing parameters that the module can use as values 
  !   for a field method. These should override default values within the module.
  ! </DATA>

  character(len=128) :: method_type
  character(len=128) :: method_name
  character(len=256) :: method_control
end type
! </TYPE> NAME="method_type"

! <TYPE NAME="method_type_short">
! <DESCRIPTION>
!   This method_type is the same as method_type except that the method_control string is not present.
!   This is used when you wish to change to a scheme within a module but do not need to pass 
!   parameters.
! </DESCRIPTION>

type, public :: method_type_short

  ! <DATA NAME="method_type_short :: method_type" TYPE="character" DIM="(128)">
  !   see method_type :: method_type above.
  ! </DATA> ! 

  ! <DATA NAME="method_type_short :: method_name" TYPE="character" DIM="(128)">
  !   see method_type :: method_name above.
  ! </DATA> 

  character(len=128) :: method_type
  character(len=128) :: method_name
end type
! </TYPE> NAME="method_type_short"


type, private :: field_type
  character(len=32) :: field_type, field_name
  integer :: model, num_methods
  type(method_type) :: methods(MAX_FIELD_METHODS)
end type  

type, private :: field_names_type
  character(len=32) :: fld_type, mod_name, fld_name
!  integer :: num_methods
end  type field_names_type


type(field_type), private :: fields(MAX_FIELDS)
type(method_type), public :: default_method

public field_manager_init, &
       field_manager_end,  &
       get_field_info,     &
       get_field_method,   &
       find_field_index,   &
       get_field_methods,  &
       parse

interface parse
   module procedure parse_real,    parse_reals,    &
                    parse_integer, parse_integers, &
                    parse_string,  parse_strings
end interface
contains

! <SUBROUTINE NAME="field_manager_init">
!   <OVERVIEW>
!     Routine to initialize the field manager.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This routine reads from a file containing formatted strings. 
!     These formatted strings contain information on which schemes are needed within
!     various modules. The field manager does not initialize any of those schemes 
!     however. It simply holds the information and is queried by the appropriate 
!     module.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call field_manager_init(nfields, table_name)
!   </TEMPLATE>

subroutine field_manager_init(nfields, table_name)

! <OUT NAME="nfields" TYPE="integer">
!   The number of fields.
! </OUT>

integer, intent(out) :: nfields

! <IN NAME="table_name" TYPE="character, optional"
!     DIM="(len=128)" DEFAULT="field_table">
!   The name of the field table. The default name is field_table.
! </IN>

character(len=128), intent(in), optional :: table_name

character(len=128) :: tbl_name
character(len=32) :: model_name
character(len=1024) :: record
integer :: n, iunit, num_methods, m, model, icount, l
logical :: flag_method

type(field_names_type) :: text_names
type(method_type) :: text_method
type(method_type_short) :: text_method_short

if (module_is_initialized) then
   nfields = num_fields
   return
endif

!!!call mpp_io_init()

if (.not.PRESENT(table_name)) then
   tbl_name = 'field_table'
else
   tbl_name = trim(table_name)
endif

if (.not. file_exist(trim(tbl_name))) then
!   <ERROR MSG="No field table available, so no fields are being registered." STATUS="NOTE">
!      The field table does not exist.
!   </ERROR>
  call error_mesg('field_manager_init', 'No field table available, so no fields are being registered.', NOTE)
nfields = 0
return
endif


!!!call mpp_open(iunit,file=trim(tbl_name), form=MPP_ASCII, action=MPP_RDONLY)
iunit = open_file(trim(tbl_name), form = 'formatted', action = 'READ') 

do while (.TRUE.)
   read(iunit,'(a)',end=89,err=99) record
   if (record(1:1) == "#" ) cycle
   read(record,*,end=79,err=79) text_names
   select case (text_names%mod_name)
   case ('atmos_mod')
      model = MODEL_ATMOS
   case ('ocean_mod')
      model = MODEL_OCEAN
   case ('land_mod')
      model = MODEL_LAND
   case ('ice_mod')
      model = MODEL_ICE
   case default
      model = MODEL_DEFAULT
   end select
   if (find_field_index(model, text_names%fld_name) < 0) then
      num_fields = num_fields + 1


!     <ERROR MSG="max fields exceeded" STATUS="FATAL">
!       Maximum number of fields for this module has been exceeded.
!     </ERROR>
      if (num_fields > MAX_FIELDS) call error_mesg('field_manager_init', 'max fields exceeded', FATAL)
      fields(num_fields)%model = model
      fields(num_fields)%field_name = lowercase(trim(text_names%fld_name))
      fields(num_fields)%field_type = lowercase(trim(text_names%fld_type))
      fields(num_fields)%num_methods = 0

! Check to see that the first line is not the only line
      if ( record(LEN_TRIM(record):LEN_TRIM(record)) == "/") cycle

      flag_method = .TRUE.
      m = 1
      do while (flag_method)
         read(iunit,'(a)',end=99,err=99) record
! If the last character in the line is / then this is the end of the field methods
         if ( record(LEN_TRIM(record):LEN_TRIM(record)) == "/") then
            flag_method = .FALSE.
            if (LEN_TRIM(record) == 1) cycle
         endif
! If the first character in the line is # then it is treated as a comment
         if (record(1:1) == "#" ) cycle

         icount = 0
         do l= 1, LEN_TRIM(record)
            if (record(l:l) == '"' ) then
               icount = icount + 1
            endif
         enddo     
!     <ERROR MSG="Too many fields in tracer entry." STATUS="FATAL">
!       There are more that 3 fields in the tracer entry. This is probably due
!       to separating the parameters entry into multiple strings. 
!       The entry should look like <BR/>       
!       "Type","Name","Control1=XXX,Control2=YYY" <BR/>
!        and not like<BR/>
!       "Type","Name","Control1=XXX","Control2=YYY"
!     </ERROR>
      if (icount > 6 ) call error_mesg('field_manager_init', 'Too many fields in tracer entry.'//trim(record), FATAL)
         if (icount == 6 ) then
         read(record,*,end=99,err=99) text_method
            fields(num_fields)%methods(m)%method_type = lowercase(trim(text_method%method_type))
            fields(num_fields)%methods(m)%method_name = lowercase(trim(text_method%method_name))
            fields(num_fields)%methods(m)%method_control = lowercase(trim(text_method%method_control))
         else
! If there is no control string then the last string can be omitted and there are only 4 '"' in the record.
            read(record,*,end=99,err=99) text_method_short
            fields(num_fields)%methods(m)%method_type = lowercase(trim(text_method_short%method_type))
            fields(num_fields)%methods(m)%method_name = lowercase(trim(text_method_short%method_name))
            fields(num_fields)%methods(m)%method_control = " "
         endif 
         fields(num_fields)%num_methods = fields(num_fields)%num_methods + 1
!     <ERROR MSG="Maximum number of methods for field exceeded" STATUS="FATAL">
!       Maximum number of methods allowed for entries in the field table has been exceeded.
!     </ERROR>
      if (fields(num_fields)%num_methods > MAX_FIELD_METHODS) &
         call error_mesg('field_manager_init', 'Maximum number of methods for field exceeded', FATAL)
         m = m + 1
      enddo
   else

!     <ERROR MSG="field with identical name and model name duplicate found, skipping" STATUS="NOTE">
!       The name of the field and the model name are identical. Skipping that field.
!     </ERROR>
         call error_mesg('field_manager_init', 'field with identical name and model name duplicate found, skipping', NOTE)
      flag_method = .TRUE.
      do while (flag_method)
         read(iunit,end=99,err=99) record
         if ( record(LEN_TRIM(record):LEN_TRIM(record)) == "/") then
            flag_method = .FALSE.
         endif
      enddo
   endif
!   read(iunit,'(a)',end=99,err=99) record(1:1)
!   if (record(1:1) /= ' ') backspace(iunit)
79 continue
enddo
         
89 continue

  call write_version_number (version, tagname)
module_is_initialized = .true.
nfields = num_fields

default_method%method_type = 'none'
default_method%method_name = 'none'
default_method%method_control = 'none'

close(iunit)
return

99 continue

!     <ERROR MSG="error reading field table" STATUS="FATAL">
!       There is an error in reading the field table.
!     </ERROR>
call error_mesg('field_manager_init', 'error reading field table', FATAL)

end subroutine field_manager_init
! </SUBROUTINE>

! <SUBROUTINE NAME="field_manager_end">
!   <OVERVIEW>
!     Destructor for field manager.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This subroutine writes to the logfile that the user is exiting field_manager and 
!     changes the initialized flag to false.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call field_manager_end
!   </TEMPLATE>
subroutine field_manager_end

integer :: log_unit

log_unit = stdlog()
   write (log_unit,'(/,(a))') 'Exiting field_manager, have a nice day ...'

module_is_initialized = .false.

end subroutine field_manager_end
! </SUBROUTINE>



! <FUNCTION NAME="find_field_index">
!   <OVERVIEW>
!     Function to return the index of the field.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This function when passed a model number and a field name will 
!     return the index of the field within the field manager. This index 
!     can be used to access other information from the field manager.
!   </DESCRIPTION>
!   <TEMPLATE>
!     value=find_field_index( model, field_name )
!   </TEMPLATE>

function find_field_index(model, field_name)
! 
!   <IN NAME="model" TYPE="integer">
!     The number indicating which model is used.
!   </IN>

integer, intent(in) :: model
character(len=*) :: field_name
integer :: find_field_index, i

find_field_index = -1

do i=1,num_fields
   if (fields(i)%model == model .and. fields(i)%field_name == lowercase(field_name)) then
      find_field_index = i
      return
   endif
enddo

return

end function find_field_index
! </FUNCTION>


! <SUBROUTINE NAME="get_field_info">
!   <OVERVIEW>
!     This routine allows access to field information given an index.
!   </OVERVIEW>
!   <DESCRIPTION>
!     When passed an index, this routine will return the type of field, 
!     the name of the field, the model which the field is associated and 
!     the number of methods associated with the field.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call get_field_info( n,fld_type,fld_name,model,num_methods )
!   </TEMPLATE>
subroutine get_field_info(n,fld_type,fld_name,model,num_methods)
!
!   <IN NAME="n" TYPE="integer">
!     The field index.
!   </IN>
integer,          intent(in)  :: n

!   <OUT NAME="fld_type" TYPE="character" DIM="(*)">
!     The field type.
!   </OUT>

!   <OUT NAME="fld_name" TYPE="character" DIM="(*)">
!     The name of the field.
!   </OUT>

!   <OUT NAME="model" TYPE="integer">
!     The number indicating which model is used.
!   </OUT>

!   <OUT NAME="num_methods" TYPE="integer">
!     The number of methods.
!   </OUT>
character (len=*),intent(out) :: fld_type, fld_name
integer, intent(out) :: model, num_methods

!   <ERROR MSG="invalid field index" STATUS="FATAL">
!     The field index is invalid because it is less than 1 or greater than the 
!     number of fields.
!   </ERROR>
if (n < 1 .or. n > num_fields) call error_mesg('get_field_info field_manager','invalid field index', FATAL)

fld_type = fields(n)%field_type
fld_name=fields(n)%field_name
model=fields(n)%model
num_methods = fields(n)%num_methods

end subroutine get_field_info
! </SUBROUTINE>


! <SUBROUTINE NAME="get_field_method">
!   <OVERVIEW>
!     A routine to get a specified method.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This routine, when passed a field index and a method index will 
!     return the method text associated with the field(n) method(m).
!   </DESCRIPTION>
!   <TEMPLATE>
!     call get_field_method( n,m,method )
!   </TEMPLATE>
subroutine get_field_method(n,m,method)
!
!   <IN NAME="n" TYPE="integer">
!     The field index.
!   </IN>
!   <IN NAME="m" TYPE="integer">
!     The method index.
!   </IN>
integer,          intent(in)  :: n, m

type(method_type) ,intent(out) :: method

!   <ERROR MSG="invalid field index" STATUS="FATAL">
!     The field index is invalid because it is less than 1 or greater than the 
!     number of fields.
!   </ERROR>
if (n < 1 .or. n > num_fields) call error_mesg('get_field_method field_manager','invalid field index', FATAL)

!   <ERROR MSG="invalid method index" STATUS="FATAL">
!     The method index is invalid because it is less than 1 or greater than 
!     the number of methods.
!   </ERROR>
if (m < 1 .or. m > fields(n)%num_methods) call error_mesg('field_manager','invalid method index', FATAL)

method = fields(n)%methods(m)

end subroutine get_field_method
! </SUBROUTINE>


! <SUBROUTINE NAME="get_field_methods">
!   <OVERVIEW>
!     A routine to obtain all the methods associated with a field.
!   </OVERVIEW>
!   <DESCRIPTION>
!     When passed a field index, this routine will return the text 
!     associated with all the methods attached to the field.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call get_field_methods( n,methods )
!   </TEMPLATE>
subroutine get_field_methods(n,methods)
!
!   <IN NAME="n" TYPE="integer">
!     The field index.
!   </IN>
integer,          intent(in)  :: n

type(method_type),intent(out) :: methods(:)

!   <ERROR MSG="invalid field index" STATUS="FATAL">
!     The field index is invalid because it is less than 1 or greater than the 
!     number of fields.
!   </ERROR>
if (n < 1 .or. n > num_fields) call error_mesg('field_manager','invalid field index', FATAL)

!   <ERROR MSG="method array too small" STATUS="FATAL">
!     The method array is smaller than the number of methods.
!   </ERROR>
if (size(methods) <  fields(n)%num_methods) call error_mesg('field manager','method array too small', FATAL)

methods(1:fields(n)%num_methods) = fields(n)%methods(1:fields(n)%num_methods)

end subroutine get_field_methods
! </SUBROUTINE>

  
! <FUNCTION NAME="parse">
!   <OVERVIEW>
!     A function to parse an integer or an array of integers, 
!     a real or an array of reals, a string or an array of strings.
!   </OVERVIEW>
!   <DESCRIPTION>
!  Parse is an integer function that decodes values from a text string.
!  The text string has the form: "label=list" where "label" is an
!  arbitrary user defined label describing the values being decoded,
!  and "list" is a list of one or more values separated by commas.
!  The values may be integer, real, or character.
!  Parse returns the number of values decoded.
!   </DESCRIPTION>
!   <TEMPLATE>
!     number = parse(text, label, value)
!   </TEMPLATE>


function parse_reals ( text, label, values ) result (parse)
!
!   <IN NAME="text" TYPE="character(len=*)">
!     The text string from which the values will be parsed.
!   </IN>
!   <IN NAME="label" TYPE="character(len=*)">
!     A label which describes the values being decoded. 
!   </IN>
!   <OUT NAME="value" TYPE="integer, real, character(len=*)">
!     The value or values that have been decoded.
!   </OUT>
!   <OUT NAME="parse" TYPE="integer">
!     The number of values that have been decoded. This allows 
!     a user to define a large array and fill it partially with 
!     values from a list. This should be the size of the value array.
!   </OUT>
character(len=*), intent(in)  :: text, label
real(r8),             intent(out) :: values(:)
include 'parse.inc'
end function parse_reals
! </FUNCTION>

function parse_integers ( text, label, values ) result (parse)
character(len=*), intent(in)  :: text, label
integer,          intent(out) :: values(:)
include 'parse.inc'
end function parse_integers

function parse_strings ( text, label, values ) result (parse)
character(len=*), intent(in)  :: text, label
character(len=*), intent(out) :: values(:)
include 'parse.inc'
end function parse_strings

!---- scalar overloads -----

function parse_real ( text, label, value ) result (parse)
character(len=*), intent(in)  :: text, label
real(r8),             intent(out) :: value
integer :: parse
real(r8) :: values(1)
   parse = parse_reals ( text, label, values )
   if (parse > 0) value = values(1)
end function parse_real

function parse_integer ( text, label, value ) result (parse)
character(len=*), intent(in)  :: text, label
integer,          intent(out) :: value
integer :: parse
integer :: values(1)
   parse = parse_integers ( text, label, values )
   if (parse > 0) value = values(1)
end function parse_integer

function parse_string ( text, label, value ) result (parse)
character(len=*), intent(in)  :: text, label
character(len=*), intent(out) :: value
integer :: parse
character(len=len(value)) :: values(1)
   parse = parse_strings ( text, label, values )
   if (parse > 0) value = values(1)
end function parse_string


end module field_manager_mod

#ifdef test_field_manager

program test

use field_manager_mod
use mpp_mod, only : mpp_exit

implicit none

integer :: i, j, nfields, num_methods, model
character(len=128) :: field_type, field_name, str
real(r8) :: param
logical :: flag
type(method_type), dimension(20) :: methods

call field_manager_init(nfields)


do i=1,nfields
   call get_field_info(i,field_type, field_name, model, num_methods)
   call get_field_methods(i,methods)
      write(*,'(/,a)') '------------------------------------------'
      write(*,'(/,a,a)') 'field_type = ', trim(field_type)
      write(*,'(a,a)') 'field_name = ', trim(field_name)
      write(*,'(a,i4)') 'model      = ', model
      write(*,'(a,i4,/)') 'num_methods= ',num_methods
      do j=1,num_methods
         write(*,'(a,i4)') 'method ',j
         write(*,'(a,a)') '   method_type=    ',trim(methods(j)%method_type)
         write(*,'(a,a)') '   method_name=    ',trim(methods(j)%method_name)
         write(*,'(a,a)') '   method_control= ',trim(methods(j)%method_control)
         if (len_trim(methods(j)%method_control) > 0) then
             str= parse_string(methods(j)%method_control,'file',flag)
             if (flag) then
                 write(*,*) 'filename = ', trim(str)
             endif
             param = parse_real(methods(j)%method_control,'param1',flag)
             if (flag) then
                 write(*,'(a,f7.3)') 'param1 = ', param
             endif             
         endif
      enddo
enddo

call field_manager_end

call mpp_exit

end program test

#endif
