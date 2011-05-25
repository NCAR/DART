! DART software - Copyright 2004 - 2011 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program html_namelist

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

implicit none

! Interactive creation of html for function documentation for DART
integer :: num, num_args, i, max_len, len, j
character(len = 299) :: c1, c2, function_name, function_type, spaces, default_value
character(len = 299), allocatable :: arg_name(:), arg_type(:)

write(22, *) '<!--============== DESCRIPTION OF A NAMELIST ========================-->'

spaces = '                                                                      '

! Put in the standard header stuff
write(22, *) '<A NAME="Namelist"></A>'
write(22, *) '<BR><HR><BR>'
write(22, *) '<P>We adhere to the F90 standard of starting a namelist with an ampersand'
write(22, *) "'&' and terminating with a slash '/'."
write(22, *) "<div class=namelist><pre>"

! Now write the namelist sequence
write(*, *) 'input the number of namelist entries'
read(*, *) num_args
allocate(arg_name(0:num_args), arg_type(0:num_args))
write(*, *) 'input the namelist name'
read(*, *) arg_name(0)

! Loop to read in the argument names
do i = 1, num_args
   write(*, *) 'input namelist entry name ', i
   read(*, *) arg_name(i)
end do

c1 = ''
do i = 1, num_args
   c1 = trim(c1) // trim(arg_name(i))
   if(i /= num_args) c1 = trim(c1) // ', '
end do

! Write out the command line if possible
write(22, *) '<em class=call>namelist / ' , trim(arg_name(0)), ' / </em>'
write(22, *) c1
write(22, *) '</pre></div>'

write(22, *) '<H3 class=indent1>Discussion</H3>'
write(22, *) '<P>'
write(22, *) 'put in general discussion ????????'
write(22, *) '</em>'
write(22, *) '</P>'

write(22, *) '<P>This namelist is read in a file called <em class=file>input.nml</em>'
write(22, *) '</P>'


! Prepare the table header
write(22, *) '<TABLE border=0 cellpadding=3 width=100%>'
write(22, *) '<TR><TH align=left>Contents    </TH>'
write(22, *) '    <TH align=left>Type        </TH>'
write(22, *) '    <TH align=left>Description </TH></TR>'

do i = 1, num_args
   write(*, *) 'input the type for ', trim(arg_name(i))
   read(*, *) arg_type(i)
end do
do i = 1, num_args
   write(*, *) 'input default value for ', trim(arg_name(i))
   read(*, *) default_value
   write(*, *) 'input description for ', trim(arg_name(i))
   read(*, 31) c1
   31 format(a299)
   c1 = trim(c1) // ' Default: ' // trim(default_value)
   
   write(22, *) '<TR><!--contents--><TD valign=top>' // trim(arg_name(i)) // '</TD>'
   write(22, *) '    <!--  type  --><TD>' // trim(arg_type(i)) // '</TD>'
   write(22, *) '<!--descript--><TD>' // trim(c1) // '</TD></TR>'
end do


! Write the end of the table
write(22, *) '</TABLE>'

write(22, *) '<P>'
write(22, *) '<!--================================================================-->'

end program html_namelist
