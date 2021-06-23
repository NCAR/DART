#########
Obs Error
#########


This directory is where to add modules that compute/set the
observational errors for different types of real-world observations.

For the 2 existing files, the data source is:

ECMWF errors: 
http://www.ecmwf.int/research/ifsdocs/CY25r1/Observations/Observations-03-3.html

NCEP errors: 
a 2005 version of the GFS observation error tables.

(Note that the return values from these modules should be 
the ERROR STANDARD DEVIATION.  In the obs_seq files, the
value stored with each observation will be the variance.)

Each center uses different errors, and these separate files
make it easy to collect these values in one place, and switch
them in and out depending on the needs of the user who is
creating new obs_seq files for DART.

Anyone who wants to contribute another error module is
more than welcome to add files here.

IMPORTANT:
Each file should have the same module name; i.e. the source
file names will differ but the module name inside the file
must be the same across all modules in this directory.

All the subroutines must also have the same names and
calling sequence. They must return appropriate values 
for each observation type that is required.  
If errors for a new observation type is added, it should
be added to all the files in this directory.

This way the user can change between error values by editing
the filename in the path_names_xxx files and recompiling
without changing the code.

Thanks to Ryan Torn for the idea and initial contributions.

