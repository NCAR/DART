! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program html_function

implicit none

! Interactive creation of html for function documentation for DART
integer :: num, num_args, i, max_len, len, j
character(len = 299) :: c1, c2, function_name, function_type, spaces
character(len = 299), allocatable :: arg_name(:), arg_type(:)
logical :: ispointer

write(22, *) '<!--============= DESCRIPTION OF A FUNCTION ========================-->'
write(22, *) ''

spaces = '                                                                      '

write(*, *) 'input the number of arguments'
read(*, *) num_args
allocate(arg_name(0:num_args), arg_type(0:num_args))

write(*, *) 'input the function name'
read(*, *) arg_name(0)
write(22, *) '<A NAME="' // trim(arg_name(0)) // '"></A>'

write(22, *) '<br />'
write(22, *) '<div class=routine>'

! Loop to read in the argument names
do i = 1, num_args
   write(*, *) 'input argument name ', i
   read(*, *) arg_name(i)
end do

c1 = ''
do i = 1, num_args
   c1 = trim(c1) // trim(arg_name(i))
   if(i /= num_args) c1 = trim(c1) // ', '
end do

! Write out the command line
write(22, *) '<em class=call> function ', trim(arg_name(0)), '(', trim(c1), ') </em>'
write(22, *) '<pre>'

max_len = 0
do i = 0, num_args
   write(*, *) 'input the type for ', trim(arg_name(i))
   read(*, *) arg_type(i)
   if (i > 0) then
      write(*, *) 'input 1 if allocatable'
      read(*, *) j
      if(j == 1) arg_type(i) = trim(arg_type(i)) // ', allocatable'
      write(*, *) 'input 1 if optional'
      read(*, *) j
      if(j == 1) arg_type(i) = trim(arg_type(i)) // ', optional'
   endif
   write(*, *) 'input 1 if pointer'
   read(*, *) j
   ispointer = (j == 1)
   if(ispointer) arg_type(i) = trim(arg_type(i)) // ', pointer'
   write(*, *) 'input dimension (0, 1, 2, 3)'
   read(*, *) j
   if(j == 1) then
      arg_type(i) = trim(arg_type(i)) // ', dimension(:)'
   else if (j == 2) then
      arg_type(i) = trim(arg_type(i)) // ', dimension(:, :)'
   else if (j == 3) then
      arg_type(i) = trim(arg_type(i)) // ', dimension(:, :, :)'
   endif
   if (i > 0 .or. ispointer) then   ! cannot give intent on pointer args
      write(*, *) '1=intent(in); 2=intent(out); 3=intent(inout)'
      read(*, *) j
      if(j == 1) then
         arg_type(i) = trim(arg_type(i)) // ', intent(in)'
      elseif(j == 2) then
         arg_type(i) = trim(arg_type(i)) // ', intent(out)'
      elseif(j == 3) then
         arg_type(i) = trim(arg_type(i)) // ', intent(inout)'
      endif
   endif

   len = len_trim(arg_type(i))
   if(len > max_len) max_len = len
end do

do i = 0, num_args
   c1 = trim(arg_type(i)) // spaces
   write(22, *) c1(1:max_len) // ' :: <em class=code>' // trim(arg_name(i)) // '</em>'
end do

write(22, *) '</pre>'
write(22, *) '</div>'
write(22, *) ' '
write(22, *) '<div class=indent1>'
write(22, *) '<!-- Description -->'

write(22, *) ' '
write(22, *) '<P>'
write(22, *) 'input function description here '
write(22, *) '</P>'
write(22, *) ' '

! Write table header
write(22, *) '<TABLE width=100% border=0 summary="argument description" cellpadding=3>'
write(22, *) '<TBODY valign=top>'


! Write the table of argument descriptions
do i = 0, num_args
   write(*, *) 'input description of ', trim(arg_name(i))
   read(*, 31) c1
   31 format(a299)
   write(22, *) '<TR><TD><em class=code> ' // trim(arg_name(i)) // &
      '  </em></TD>' 
   if (i == 0) then
      write(22, *) '    <TD> Function return value. ' // trim(c1) // '</TD>'
   else
      write(22, *) '    <TD>' // trim(c1) // '</TD>'
   endif
   write(22, *) '</TR>'
end do


!End of table
write(22, *) '</TBODY>'
write(22, *) '</TABLE>'
write(22, *) '</div>'
write(22, *) '<br />'
write(22, *) ''
write(22, *) '<!--================================================================-->'

write(*,*) 'output file is fort.22 -- move or rename before running program again'

end program html_function

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
