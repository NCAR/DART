# Shell script to work with asynchronous filter integration
# This script needs to be piped to the filter program with the
# filter namelist async variable set to .true.
# Make sure that file filter_control is cleared out by a higher level
# script.

ls filter_control > .async_garb
if($status != 0) then
   echo FATAL ERROR IN ASYNC FILTER: filter_control file not found
   exit(4)
endif

# First line of filter_control should have number of model states to be integrated
set num = `head -1 filter_control`
echo number of lines expected is $num


# Would be nice to check to make sure total filter_control
# file length is consistent with this number of model sets
#wc -l filter_control


# Create a directory for each member to run in for namelists
for i = 1, 10
   echo i
endfor

