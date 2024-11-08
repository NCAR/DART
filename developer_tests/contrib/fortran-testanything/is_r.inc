! Template parameter: wp (working precision)
! Template free identifiers: testline, tests
subroutine isabs(got, expected, eps, msg)
   real(kind=wp), intent(in) :: got, expected
   character(len=*), intent(in), optional :: msg
   real(kind=wp), intent(in), optional :: eps
   character(len=:), allocatable :: testmsg, idmsg
   character(len=120) gotmsg, expectedmsg
   real(kind=wp) tolerance
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
   write (unit=gotmsg,      fmt='(A,G0)') '     got: ', got
   write (unit=expectedmsg, fmt='(A,G0)') 'expected: ', expected

   if (present(eps)) then
      tolerance = eps
   else
      tolerance = epsilon(got)
   end if
   ! eps = 0.5e-10_wp
   ! Absolute accuracy within the 10 least significant digits
   good = abs(got - expected) < tolerance
   call testline(good, testmsg, idmsg, gotmsg, expectedmsg)
end

subroutine isrel(got, expected, eps, msg)
   real(kind=wp), intent(in) :: got, expected
   character(len=*), intent(in), optional :: msg
   real(kind=wp), intent(in), optional :: eps
   real(kind=wp) tolerance

   ! eps = (abs(a) + abs(b)) * 0.5e-10_wp
   ! Relative accuracy within the 10 most significant digits
   tolerance = (abs(got) + abs(expected))
   if (present(eps)) then
      tolerance = tolerance * eps
   else
      tolerance = tolerance * epsilon(got)
   end if
   call isabs(got, expected, tolerance, msg)
end

subroutine isnear(got, expected, eps, msg)
   real(kind=wp), intent(in) :: got, expected
   character(len=*), intent(in), optional :: msg
   real(kind=wp), intent(in), optional :: eps
   character(len=:), allocatable :: testmsg, idmsg
   character(len=120) gotmsg, expectedmsg
   real(kind=wp) tolerance
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
   write (unit=gotmsg,      fmt='(A,G0)') '     got: ', got
   write (unit=expectedmsg, fmt='(A,G0)') 'expected: ', expected

   if (present(eps)) then
      tolerance = eps
   else
      tolerance = epsilon(got) ! minimun eps for which 1 + eps /= 1
   end if
   ! Relative accuracy around 1.0_wp
   ! Semantics of isnear means using <=, and not <, c.f. epsilon(got)
   good = abs(got / expected - 1.0_wp) <= tolerance
   call testline(good, testmsg, idmsg, gotmsg, expectedmsg)
end

