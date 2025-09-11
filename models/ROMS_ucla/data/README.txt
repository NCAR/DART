At one point, the grid information, vertical information, and variable size information all came from
different files. Hernan then extended ROMS such that the restart file contained all that information.

If you get an errors from "routine: get_grid", you are most likely using an old ROMS file.

Files capable of testing can be found at http://www.image.ucar.edu/pub/DART/ROMS/dart_roms_test_data.tar.gz 

For people with access to NCAR's supercomputer, they exist on spinning disk at
/glade/p/image/DART_test_cases/roms/wc13
