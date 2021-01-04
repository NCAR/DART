#!/bin/csh 
# 
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

cd ../work

# specify the directry where the original/processed data is located
set BASEOBSDIR = /scratch/04270/l_jing/DART_NoahMP/observation/GRACE

set subobsdir = 'dart_seq_grace'

foreach YEAR  ( 2003 2004 2005 2006 2007 2008 2009 2010 ) 
   foreach FILE  ($BASEOBSDIR/nc_grace/*${YEAR}* ) 

      echo "Converting $FILE ..."
      
      echo $FILE >! file_list.txt
   
      # run convert_daily_grace to begin the processing
      ./convert_daily_grace || exit 1
 
      mv -v obs_seq.* ${BASEOBSDIR}/${subobsdir}
      
   end
end

exit 0


