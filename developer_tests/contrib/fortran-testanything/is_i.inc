! Template parameter: wp (working precision)
! Template free identifiers: testline, tests
subroutine is(got, expected, msg)
   integer(kind=wp), intent(in) :: got, expected
   character(len=*), intent(in), optional :: msg
   character(len=:), allocatable :: testmsg, idmsg
   character(len=120) gotmsg, expectedmsg
   logical good

   if (present(msg)) then
      allocate(character(len=len_trim(msg)+20) :: testmsg, idmsg)
      write (unit=idmsg, fmt='(A,A,A)') 'Failed test: "', trim(msg), '"'
      testmsg = trim(msg)
   else
      allocate(character(len=30) :: testmsg, idmsg)
      write (unit=idmsg, fmt='(A,I0)') 'Failed test no. ', tests + 1
      testmsg = ""
   end if
   write (unit=gotmsg,      fmt='(A,I0)') '     got: ', got
   write (unit=expectedmsg, fmt='(A,I0)') 'expected: ', expected

   good = got == expected
   call testline(good, testmsg, idmsg, gotmsg, expectedmsg)
end
