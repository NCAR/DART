program test_optional_namelist

use utilities_mod, only : initialize_utilities, finalize_utilities,&
                          find_namelist_in_file, check_namelist_read, &
                          do_nml_file, do_nml_term, nmlfileunit, &
                          NAMELIST_NOT_PRESENT

use test ! fortran-testanything

implicit none

integer :: iunit, io

integer :: beans = 1
integer :: good  = 2
integer :: heart = 3

namelist /namelist_nml/ beans, good, heart

call initialize_utilities('test_optional_namelist')

! Read the namelist entry
call find_namelist_in_file("input.nml", "namelist_nml", iunit, optional_nml=.true.)
print *, "Namelist file unit: ", iunit
if (iunit /= NAMELIST_NOT_PRESENT) then
    read(iunit, nml = namelist_nml, iostat = io)
    call check_namelist_read(iunit, io, "namelist_nml")
endif

if (do_nml_file()) write(nmlfileunit, nml=namelist_nml)
if (do_nml_term()) write(     *     , nml=namelist_nml)

call ok(beans == 1, "Beans should be 1")
call ok(good == 2, "Good should be 2")
call ok(heart == 3, "Heart should be 3")

call finalize_utilities()

end program test_optional_namelist
