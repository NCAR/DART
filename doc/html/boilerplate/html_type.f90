! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without chiteme, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program html_type

implicit none

! Interactive creation of html for derived type documentation for DART
integer :: num, num_items, i, max_len, len, j
character(len = 299) :: c1, c2, type_name, spaces
character(len = 299), allocatable :: item_name(:), item_type(:)

write(22, *) '<!--============== DESCRIPTION OF A DERIVED TYPE ===================-->'
write(22, *) ''

spaces = '                                                                      '

! get the count and name of the type
write(*, *) 'input the number of items in this type'
read(*, *) num_items
allocate(item_name(0:num_items), item_type(0:num_items))
write(*, *) 'input the derived type name'
read(*, *) item_name(0)
write(*, *) 'enter 1 if contents are private'
read(*, *) j

! Put in the standard header stuff
write(22, *) '<A NAME="' // trim(item_name(0)) // '"></A>'
write(22, *) '<br />'
write(22, *) '<div class=type><pre>'
write(22, *) '<em class=call>type ' // trim(item_name(0)) // '</em>'
if (j==1) write(22, *) '   private'

! Loop to read in the item names
do i = 1, num_items
   write(*, *) 'input name for item ', i
   read(*, *) item_name(i)
   write(*, *) 'input the type for ', trim(item_name(i))
   read(*, *) item_type(i)
   write(22, *) '   ' // trim(item_type(i)) // ' :: ' // trim(item_name(i))
end do

write(22, *) 'end type ' // trim(item_name(0))
write(22, *) '</pre>'
write(22, *) '</div>'

write(22, *) ''
write(22, *) '<div class=indent1>'
write(22, *) '<!-- Description -->'
write(22, *) ''
write(22, *) '<P>'
write(22, *) 'put in any general type discussion information'
write(22, *) '</P>'

! Prepare the table header
write(22, *) '<TABLE border=0 cellpadding=3 width=100% summary="derived type description">'
write(22, *) '<THEAD align=left>'
write(22, *) '<TR><TH> Component </TH>'
write(22, *) '    <TH> Description </TH>'
write(22, *) '</TR>'
write(22, *) '</THEAD>'
write(22, *) ''
write(22, *) '<TBODY valign=top>'
do i = 1, num_items
   write(*, *) 'input description for ', trim(item_name(i))
   read(*, 31) c1
   31 format(a299)
   
   write(22, *) '<TR><TD> ' // trim(item_name(i)) // ' </TD>'
   write(22, *) '    <TD> ' // trim(c1) // ' </TD>'
   write(22, *) '</TR>'
end do


! Write the end of the table
write(22, *) '</TBODY>'
write(22, *) '</TABLE>'
write(22, *) ''
write(22, *) '</div>'
write(22, *) '<br />'
write(22, *) ''

write(22, *) '<!--================================================================-->'

write(*,*) 'output file is fort.22 -- move or rename before running program again'

end program html_type

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
