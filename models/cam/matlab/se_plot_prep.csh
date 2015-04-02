#!/bin/tcsh

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

exit
