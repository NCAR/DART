! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program html_namelist

implicit none

! Interactive creation of html for namelist documentation for DART
integer :: num, num_args, i, max_len, len, j
character(len = 299) :: c1, c2, function_name, function_type, spaces
character(len = 299), allocatable :: arg_name(:), arg_type(:)

write(22, *) '<!--============== DESCRIPTION OF A NAMELIST =======================-->'
write(22, *) ''

spaces = '                                                                      '

! Put in the standard header stuff
write(22, *) '<A NAME="Namelist"></A>'
write(22, *) '<div class="top">[<a href="#">top</a>]</div><hr />'
write(22, *) '<H2>NAMELIST</H2>'
write(22, *) '<P>'
write(22, *) 'This namelist is read from the file <em class=file>input.nml</em>.'
write(22, *) 'Namelists start with an ampersand'
write(22, *) "'&amp;' and terminate with a slash '/'."
write(22, *) "Character strings that contain a '/' must be"
write(22, *) 'enclosed in quotes to prevent them from'
write(22, *) 'prematurely terminating the namelist.'
write(22, *) '</P>'
write(22, *) "<div class=namelist>"
write(22, *) "<pre>"

! Now write the namelist sequence
write(*, *) 'input the number of namelist entries'
read(*, *) num_args
allocate(arg_name(0:num_args), arg_type(0:num_args))
write(*, *) 'input the namelist name'
read(*, *) arg_name(0)

! Write out the namelist start
write(22, *) '&amp;' , trim(arg_name(0))

! Loop to read in the argument names
do i = 1, num_args
   write(*, *) 'input namelist entry ', i, ' as name = value '
   read(*, 31) c1
   read(c1, *) arg_name(i)
   write(22, *) '   ' // trim(c1) // ', '
end do

write(22, *) '/'
write(22, *) '</pre></div>'

write(22, *) '<P>'
write(22, *) 'put in any general namelist discussion information'
write(22, *) '</P>'
write(22, *) '<br /><br />'

! Prepare the table header
write(22, *) '<div>'
write(22, *) '<TABLE border=0 cellpadding=3 width=100% summary="namelist description">'
write(22, *) '<THEAD align=left>'
write(22, *) '<TR><TH> Item </TH>'
write(22, *) '    <TH> Type </TH>'
write(22, *) '    <TH> Description </TH>'
write(22, *) '</TR>'
write(22, *) ''
write(22, *) '<TBODY valign=top>'
do i = 1, num_args
   write(*, *) 'input the type for ', trim(arg_name(i))
   read(*, *) arg_type(i)
end do
do i = 1, num_args
   write(*, *) 'input description for ', trim(arg_name(i))
   read(*, 31) c1
   31 format(a299)
   
   write(22, *) '<TR><TD> ' // trim(arg_name(i)) // ' </TD>'
   write(22, *) '    <TD> ' // trim(arg_type(i)) // ' </TD>'
   write(22, *) '    <TD> ' // trim(c1) // ' </TD>'
   write(22, *) '</TR>'
end do


! Write the end of the table
write(22, *) '</TBODY>'
write(22, *) '</TABLE>'
write(22, *) '</div>'

write(22, *) '<br /><br />'

write(22, *) '<!--================================================================-->'

write(*,*) 'output file is fort.22 -- move or rename before running program again'

end program html_namelist

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
