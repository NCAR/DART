#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# The special values(!) which may be found on the CLM initial/restart files.
# CLM 3.6.71 set spvals = ('1e+36' '-1e+36')
# ICE 4.0
#  -1.836 is the value of Tsfcn (and eicen? no; embedded in larger numbers) 
#         over land for member 1, but NOT 2 or 3 (or ....?)
set spvals = ('\-1\.836' '1e+30' )

if ($#argv == 0) then
   echo "Usage: edit find_FillValue_vars.csh to provide spvalS"
   echo "       find_FillValue_vars.csh file.nc"
   exit
else
   set file = $1
   if (-e $file:r.FillValue_list) then
      echo $file:r.FillValue_list exists: move and retry
      exit
   endif
endif

# Generate a list of variables from the input file.
set head = $file:r.head
ncdump -h $file | grep double  >! $head
ncdump -h $file | grep ' int ' >> $head

touch $file:r.FillValue_list

set num_vars = `wc -l $head`
echo num_vars = $num_vars
foreach spv ($spvals)
   # Make a list of vars which have it.
   set spvars = ()
   set n = 1
   # Check each variable for this spval.  
   while ($n <= $num_vars[1])
      head -$n $head | tail -1 >! varstring
      # Replace the ( with a ' ' so that the variable name becomes a separate word.
      set string = `sed -e "s#(# #g" varstring`
      set var = $string[2]
      ncks -v $var $file | grep "$spv" >! spvals
      set num_spvals = `wc -l spvals`
      if ($num_spvals[1] > 0) then
         set spvars = ($spvars $var)
      endif
   
      @ n++
   end
   echo spval = $spv
   echo spvars = $spvars
   
   # convert the word list into the string which will be written out 
   # and read in by analyses2initial.csh and used to add _FillValue 
   # to the necessary variables in CLM initial files.
   echo "#spvars = $#spvars"
   set n = 0
   set l = 1
   @ num_lines = ($#spvars / 10) + 1
   while ($l <= $num_lines)
      @ nend = $l * 10
      if ($l == $num_lines) @ nend = $#spvars 
      @ n++
      set var_list = $spvars[$n]
      while ($n < $nend)
         @ n++
         set var_list = "$var_list|$spvars[$n]"
      end
      echo $spv $var_list >> $file:r.FillValue_list
      @ l++
   end
end
   
rm varstring spvals

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

