! Copyright 2015 Dennis Decker Jensen
! See <http://testanything.org> and <https://metacpan.org/pod/Test::More>
! Tectonics: gfortran -g -Wall -Wextra -std=f2008ts -c test.f08

module test_base
   use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
   implicit none

   ! Kept as variables instead of aliases,
   ! so that test output or diagnostic output can be redirected
   integer :: test_unit = output_unit, diag_unit = error_unit

   integer :: tests = 0, todos = 0
   character(len=120) :: todomsg = ""

   interface todo
      module procedure todo_i, todo_s, todo_s_i, todo
   end interface

contains

   subroutine diag(msg)
      character(len=*), intent(in) :: msg
      write (diag_unit, '("# ",A)') trim(msg) ! only trailing spaces
   end

   subroutine note(msg)
      character(len=*), intent(in) :: msg
      write (test_unit, '("# ",A)') trim(msg)
   end

   subroutine testline(ok, msg, idmsg, gotmsg, expectedmsg)
      logical, intent(in) :: ok
      character(len=*), intent(in) :: msg, idmsg, gotmsg, expectedmsg

      tests = tests + 1
      if (.not. ok) call out("not ")
      write (test_unit, '("ok ",I0)', advance="NO") tests

      if (msg /= "" .or. todos > 0) call out(" - ")

      if (msg /= "") call out(trim(msg))

      if (todos > 0) then
         todos = todos - 1
         if (msg /= "") call out(" ")
         call out("# TODO")
         if (todomsg .ne. "") then
            call out(": ")
            call out(trim(todomsg))
         end if
      end if
      if (todos == 0) todomsg = ""

      write (test_unit, *) ""

      if (.not. ok) then
         ! 3 spaces prepended = 4 spaces indentation after # on diag
         if (idmsg /= "") call diag("   " // idmsg)
         if (gotmsg /= "") call diag("   " // gotmsg)
         if (expectedmsg /= "") call diag("   " // expectedmsg)
      end if
   contains
      subroutine out(str)
         character(len=*), intent(in) :: str
         write (test_unit, '(A)', advance="NO") str
      end
   end subroutine testline

   subroutine ok(condition, msg)
      logical, intent(in) :: condition
      character(len=*), intent(in), optional :: msg
      if (present(msg)) then
         call testline(condition, msg, "", "", "")
      else
         call testline(condition,  "", "", "", "")
      end if
   end

   subroutine pass(msg)
      character(len=*), intent(in), optional :: msg
      call ok(.true., msg)
   end

   subroutine fail(msg)
      character(len=*), intent(in), optional :: msg
      call ok(.false., msg)
   end

   subroutine todo_s_i(msg, howmany)
      character(len=*), intent(in) :: msg
      integer, intent(in) :: howmany
      todomsg = msg
      todos = howmany
   end

   subroutine todo
      call todo_s_i("", 1)
   end 

   subroutine todo_s(msg)
      character(len=*), intent(in) :: msg
      call todo_s_i(msg, 1)
   end

   subroutine todo_i(howmany)
      integer, intent(in) :: howmany
      call todo_s_i("", howmany)
   end 

end module test_base

module test_planning
   use test_base, only: test_unit, tests
   implicit none

   integer, private :: planned = 0 

contains

   subroutine bail_out(msg)
      character(len=*), intent(in), optional :: msg
      if (present(msg)) then
         write (test_unit, '("Bail out! ",A)') msg
      else
         write (test_unit, '("Bail out!")')
      end if
      stop
   end

   subroutine plan(tests)
      integer, intent(in) :: tests

      select case (tests)
      case (:-1)
         call bail_out("A plan with a negative number of tests")
      case (0)
         write (test_unit, '("1..0")')
         stop ! The same as skip_all without a given reason
      case (1:)
         if (planned > 0) &
            & call bail_out("More than one plan in test output")
         planned = tests
         write (test_unit, '("1..",I0)') planned
      end select
   end

   subroutine done_testing(howmany)
      integer, intent(in), optional :: howmany

      ! Put plan at the end of test output
      if (present(howmany)) then
         call plan(howmany)
      else
         if (planned == 0) call plan(tests)
         ! else - We already have a plan
      end if
   end

   subroutine skip_all(msg)
      character(len=*), intent(in), optional :: msg
      if (present(msg)) then
         write (test_unit, '("1..0 # Skipped: ",A)') msg
      else
         write (test_unit, '("1..0 # Skipped all")')
      end if
      stop
   end

end module test_planning

! Template instances of integer kinds for "is"

module is_i8_mod
   use, intrinsic :: iso_fortran_env, only: wp => int8
   use, non_intrinsic :: test_base, only: testline, tests
contains
   include "is_i.inc"
end

module is_i16_mod
   use, intrinsic :: iso_fortran_env, only: wp => int16
   use, non_intrinsic :: test_base, only: testline, tests
contains
   include "is_i.inc"
end

module is_i32_mod
   use, intrinsic :: iso_fortran_env, only: wp => int32
   use, non_intrinsic :: test_base, only: testline, tests
contains
   include "is_i.inc"
end

module is_i64_mod
   use, intrinsic :: iso_fortran_env, only: wp => int64
   use, non_intrinsic :: test_base, only: testline, tests
contains
   include "is_i.inc"
end

module is_i
   use is_i8_mod,  only: is_i8  => is
   use is_i16_mod, only: is_i16 => is
   use is_i32_mod, only: is_i32 => is
   use is_i64_mod, only: is_i64 => is
   interface is
      module procedure is_i8, is_i16, is_i32, is_i64
   end interface
end

! Template instances of real kinds for "is"

module is_r32_mod
   use, intrinsic :: iso_fortran_env, only: wp => real32
   use, non_intrinsic :: test_base, only: testline, tests
contains
   include "is_r.inc"
end

module is_r64_mod
   use, intrinsic :: iso_fortran_env, only: wp => real64
   use, non_intrinsic :: test_base, only: testline, tests
contains
   include "is_r.inc"
end

module is_r
   use is_r32_mod, only: isrel_r32 => isrel, isabs_r32 => isabs, &
                           & isnear_r32 => isnear
   use is_r64_mod, only: isrel_r64 => isrel, isabs_r64 => isabs, &
                           & isnear_r64 => isnear
   interface isrel
      module procedure isrel_r32, isrel_r64
   end interface

   interface isabs
      module procedure isabs_r32, isabs_r64
   end interface

   interface isnear
      module procedure isnear_r32, isnear_r64
   end interface
end

module test_more
   use test_base, only: testline, tests, test_unit
   use test_planning, only: bail_out ! for negative skips
   use is_i, only: is, is_i8, is_i16, is_i32, is_i64
   use is_r, only: isabs, isrel, isnear, &
                     & isabs_r32, isrel_r32, isnear_r32, &
                     & isabs_r64, isrel_r64, isnear_r64

   ! Complex numbers cannot be compared, hence no is_c module

   implicit none

   interface skip
      module procedure skip_i, skip_s, skip_s_i, skip
   end interface

   interface is
      module procedure is_s, is_l
   end interface

contains

   subroutine skip_s_i(msg, howmany)
      character(len=*), intent(in) :: msg
      integer, intent(in) :: howmany
      character(len=120) skipmsg
      integer i

      if (howmany <= 0) then
         call bail_out("Skipped non-positive number of tests")
      end if

      if (msg == "") then
         skipmsg = "# SKIP"
      else
         skipmsg = "# SKIP: " // trim(msg)
      end if

      do i = 1, howmany
         tests = tests + 1
         write (test_unit, '("ok ",I0," ",A)') tests, trim(skipmsg)
      end do
   end

   subroutine skip
      call skip_s_i("", 1)
   end 

   subroutine skip_s(msg)
      character(len=*), intent(in) :: msg
      call skip_s_i(msg, 1)
   end

   subroutine skip_i(howmany)
      integer, intent(in) :: howmany
      call skip_s_i("", howmany)
   end 

   ! Duplicates of is_i routines in file is_i.inc and ditto is_r
   ! They are not factored any further, because it is easier
   ! to see all the output together rather than in separate routines

   subroutine is_s(got, expected, msg)
      character(len=*), intent(in) :: got
      character(len=*), intent(in) :: expected
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
      write (unit=gotmsg,      fmt='(A,A,A)') '     got: "', got, '"'
      write (unit=expectedmsg, fmt='(A,A,A)') 'expected: "', expected, '"'

      good = got == expected
      call testline(good, testmsg, idmsg, gotmsg, expectedmsg)
   end

   subroutine is_l(got, expected, msg)
      logical, intent(in) :: got, expected
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
      write (unit=gotmsg,      fmt='(A,L1)') '     got: ', got
      write (unit=expectedmsg, fmt='(A,L1)') 'expected: ', expected

      good = got .eqv. expected
      call testline(good, testmsg, idmsg, gotmsg, expectedmsg)
   end

end module test_more

module test
   use test_base, only: test_unit, diag_unit, &
                     & ok, diag, note, pass, fail, todo
   use test_planning, only: plan, done_testing, skip_all, bail_out
   use test_more, only: is, isabs, isrel, isnear, skip
end module test

