program transform_names

! Test the functions that will go into model_mod for use by aether_to_dart
! and dart_to_aether to convert Aether field names to CF compliant DART names.
!

! use netcdf
! use typesizes
use types_mod, only : MISSING_I

! Why not? character (len=NF90_MAX_NAME) :: aether_name, dart_name
character (len=256) :: aether_name, dart_name

aether_name = ''
read '(A)', aether_name

dart_name = aeth_name_to_dart(aether_name)
print*, trim(dart_name), '||end'

contains
!-----------------------------------------------------------------------------
! Translate an Aether field name (not CF-compliant) into a form filter likes.
! E.g.  'Perp.\ Ion\ Velocity\ \(Meridional\)\ \(O+\)', ->
!       'Opos_Perp_Ion_Vel_Merid' 
function aeth_name_to_dart(varname)

! character(len=NF90_MAX_NAME), intent(in) :: varname
character(len=256), intent(in) :: varname

! character(len=NF90_MAX_NAME) :: aeth
character(len=256) :: aeth
character(len=128) :: aeth_name_to_dart
character(len=32)  :: parts(8), var_root
integer :: char_num, first, i_parts, aeth_len, end_str

aeth = trim(varname)
aeth_len = len_trim(varname)
parts = ''

! Look for the last ' '.  The characters after that are the species.
! If there's no ' ', the whole string is the species.
char_num = 0
char_num = scan(trim(aeth),' ',back=.true.)
print*,'species blank at ',char_num
var_root = aeth(char_num+1:aeth_len)
print*,'species var_root = ', var_root
end_str = char_num
  
! purge_chars removes unwanted [()\] 
! Remove blanks from front and end.
parts(1) = purge_chars( trim(var_root),')(\' )
print*,'parts(1) = ',parts(1)

! Tranform remaining pieces of varname into DART versions.
char_num = MISSING_I
first = 1
i_parts = 2
do
   ! This returns the position of the first blank *within the substring* passed in.
   char_num = scan(aeth(first:end_str),' ',back=.false.)
   print*,'char_num, aeth substring = ',char_num, aeth(first:end_str),'||end'
   if (char_num > 0 .and. first < aeth_len) then
      parts(i_parts) = purge_chars(aeth(first:first+char_num-1), '.)(\' )

      first   = first + char_num 
      print*,'parts(i_parts), first, aeth_len = ' ,parts(i_parts), first , aeth_len
      i_parts = i_parts + 1
   else
      exit
   endif
enddo

! Construct the DART field name from the parts
aeth_name_to_dart = trim(parts(1))
i_parts = 2
do
if (trim(parts(i_parts)) /= '') then
   aeth_name_to_dart = trim(aeth_name_to_dart)//'_'//trim(parts(i_parts))
   print*,'i_parts, aeth_name_to_dart = ' ,i_parts, aeth_name_to_dart 
   i_parts = i_parts + 1
else
   exit
endif
enddo

end function aeth_name_to_dart

!-----------------------------------------------------------------
! Replace undesirable characters with better.

function purge_chars(ugly_string, chars)

character (len=*), intent(in) :: ugly_string, chars
character (len=32) :: purge_chars, temp_str

integer :: char_num, end_str, pm_num

! Trim is not needed here
purge_chars = ugly_string
char_num = MISSING_I
do 
   ! Returns 0 if chars are not found
   char_num = scan(trim(purge_chars),chars)
   ! Need to change it to a char that won't be found by scan in the next iteration,
   ! and can be easily removed.
   print*,'purge_chars: purge_chars, char_num = ',trim(purge_chars),' ', char_num
   if (char_num > 0) then
      purge_chars(char_num:char_num) = ' '
   else
      exit
   endif
enddo

! Replace + and - with pos and neg.  Assume there's only 1.
temp_str = trim(adjustl(purge_chars))
end_str = len_trim(temp_str)
pm_num = scan(trim(temp_str),'+-',back=.false.)
if (pm_num == 0) then
   purge_chars = trim(temp_str)
else
   if (temp_str(pm_num:pm_num) == '+') then
      purge_chars = temp_str(1:pm_num-1)//'pos'
   else if (temp_str(pm_num:pm_num) == '-') then
      purge_chars = temp_str(1:pm_num-1)//'neg'
   endif
   if (pm_num+1 <= end_str) &
      purge_chars = trim(purge_chars)//temp_str(pm_num+1:end_str)
endif


end function purge_chars

end program

