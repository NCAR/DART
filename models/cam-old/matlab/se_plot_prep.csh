#!/bin/tcsh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download


if ($#argv == 0) then
   echo 'Usage: se_plot_prep.csh name_root "copies"'
   exit
endif

set name = $1
set copies = ($2)
echo "copies = $copies"

if (! -d $name) then
   mkdir ${name}
endif

cd ${name}

foreach c ($copies)
   ncks -F -d copy,${c} -o ${name}.c${c}.nc ../${name}.nc
end

exit 0


