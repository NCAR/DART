#!/bin/tcsh

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$

# Single use script to back up $project data to $campaign.
# Submit from casper:$CASEROOT

setenv CASEROOT $cwd
set CASE           = $CASEROOT:t
# set local_arch     = `./xmlquery DOUT_S_ROOT --value`
# set components     = (lnd  atm ice  rof)
# set models         = (clm2 cam cice mosart)
set components     = (esp cpl lnd  rof)
# set components     = (cpl)

set project    = /glade/p/nsc/ncis0006/Reanalyses
set campaign   = /gpfs/csfs1/cisl/dares/Reanalyses
set year  = 2011

cd ${project}/${CASE}

#-----------------------------------------------------
# echo "Forcing starts at "`date`
# cd cpl/hist
#    
# mkdir ../Uncompressed
# set n = 1
# while ($n <= 80)
#    set nn = `printf %04d $n`
#    mkdir ../Uncompressed/${nn}
#    mv ${nn}/*.nc ../Uncompressed/${nn}
#    @ n++
# end
# 
# # It appears that leaving / at the end of the source dir 
# # allows '' to be appended to the destination directory,
# # which is what I want for the $MEMB directories.
# set yr_mo = `printf %4d-%02d ${year} 3`
# ${CASEROOT}/mv_to_campaign.csh $CASE ${yr_mo} ${project}/${CASE}/cpl/hist/ \
#                                              ${campaign}/${CASE}/cpl/hist
# cd ../..

#-----------------------------------------------------
# cd esp/hist
# echo " "
# echo "Location for obs space is `pwd`"
# 
# set month = 1
# while ($month <= 3) 
#    set yr_mo = `printf %4d-%02d ${year} ${month}`
#    
#    if (! -d $yr_mo) mkdir $yr_mo
#    if (-f Diags_NTrS_${yr_mo}.tgz) then
#       mv Diags_NTrS_${yr_mo}.tgz                 $yr_mo
#       mv ${CASE}.cam_obs_seq_final.${yr_mo}.tgz  $yr_mo
#    endif
#    ${CASEROOT}/mv_to_campaign.csh $CASE ${yr_mo} \
#              ${project}/${CASE}/esp/hist/$yr_mo \
#             ${campaign}/${CASE}/esp/hist
#    @ month++
# end
#    
# cd ../..
#    
   
set m = 1
while ($m <= $#components)
   ls $components[$m]/hist/0001/*.h[a0]* >& /dev/null
   if ($status == 0) then
      cd $components[$m]/hist
      echo " "
      echo "Location for history is `pwd`"

      # I wanted to move compressed files to campaign storage,
      # but leave uncompressed files here.  So the compression
      # created new copies instead of replacing the uncompressed.
      # I need to separate them to allow mv_to_campaign.csh to work
      # on the directory.
      ls 0001/*.gz
      if ($status == 0) then
         mkdir ../Uncompressed
         set n = 1
         while ($n <= 80)
            set nn = `printf %04d $n`
            mkdir ../Uncompressed/${nn}
            mv ${nn}/*.nc ../Uncompressed/${nn}
            @ n++
         end
      endif
   else if ($components[$m] == esp) then
      # Don't skip this component.
      # The obs space files (should) have been moved into YYYY-MM,
      # which will be handled just like the member directories in other components.
      cd $components[$m]/hist
      echo " "
      echo "Location for history is `pwd`"
   else
      echo "------------ "
      echo "Skipping $components[$m]/hist"
      echo "------------ "
      @ m++
      continue
   endif

   ${CASEROOT}/mv_to_campaign.csh $CASE ${year} \
      ${project}/${CASE}/$components[$m]/hist/ \
     ${campaign}/${CASE}/$components[$m]/hist
   cd ../..

   @ m++
end

exit

# <next few lines under version control, do not edit>
# $URL$
# $Id$
# $Revision$
# $Date$
