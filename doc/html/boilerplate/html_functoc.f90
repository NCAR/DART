! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without chiteme, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program html_functoc

implicit none

! Interactive creation of html for interface list documentation for DART
integer :: num, num_items, i, max_len, len, j
character(len = 299) :: c1, c2, module_name, item_name, spaces

write(22, *) '<!--============== DESCRIPTION OF PUBLIC INTERFACES ================-->'
write(22, *) ''

spaces = '                                                                      '

! get the name of the module
write(*, *) 'input the name of this module'
read(*, *) module_name

! Put in the standard header stuff
write(22, *) '<A NAME="Interface"></A>'
write(22, *) '<div class="top">[<a href="#">top</a>]</div><hr />'
write(22, *) '<H2><PUBLIC INTERFACES</H2>'
write(22, *) ''
write(22, *) "<TABLE summary='public interfaces'>"
write(22, *) '<TR><TD><em class=call>use ' // trim(module_name), ', only : </em></TD>'

write(*, *) 'input the names of the subroutines and functions in this module'
write(*, *) 'input "end" when done'
readloop: do i=1, 10000
   write(*, *) 'next routine name ("end" to finish)'
   read(*, *) item_name
   if (trim(item_name) == "end") exit readloop
   if (i==1) then
      write(22, *) '                   <TD><A HREF="#'//trim(item_name)//'"> '//trim(item_name)//' </A></TD></TR>'
   else
      write(22, *) '<TR><TD>&nbsp;</TD><TD><A HREF="#'//trim(item_name)//'"> '//trim(item_name)//' </A></TD></TR>'
   endif
end do readloop

write(22, *) '</TABLE>'

write(22, *) ''
write(22, *) '<P>'
write(22, *) '   A note about documentation style.'
write(22, *) '   Optional arguments are enclosed in brackets'
write(22, *) '   <em class=optionalcode>[like this]</em>.'
write(22, *) '</P>'
write(22, *) ''

write(22, *) '<!--================================================================-->'

write(*,*) 'output file is fort.22 -- move or rename before running program again'

end program html_functoc

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
