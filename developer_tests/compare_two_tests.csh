#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# Script to compare the output files from running test_dart.csh on two
# code bases (in preparation to merge, for example). The files under version
# control are not considered, and the executables themselves are not considered.
#
# If the file is a netCDF file, the compare_states utility is run and will only
# create output if there is a difference in one of the numeric fields.
# If the file is a text file, 'xxdiff' is used to compare the files.
# Otherwise, the file must be a binary file and 'cmp' is used.
#
# The list of files to check is created from the BRANCH1 directory.
# Nuisance files (dart_log.*, input*default) and compiled executables
# are removed from each directory to facilitate meaningful comparison.

# FIXME ... this script should not be run in either the BRANCH1 or the BRANCH2
# because the compare_states executable will be removed by the findexec step.

set nonomatch
set SNAME = $0

set BRANCH1 = /glade/work/thoar/DART/clean_rma_trunk
set BRANCH2 = /glade/work/thoar/DART/candidate_rma_trunk
set orgdir = `pwd`

cp $BRANCH1/assimilation_code/programs/compare_states/work/compare_states .
cp $BRANCH1/assimilation_code/programs/compare_states/work/input.nml      .

# We want to remove the newly-compiled executables.
# long-winded alias to find binary executables and not mkmf_ for example.
alias findexec "find . -type f -print | grep -v '\.svn' | xargs file | grep -i executable | grep -v 'text executable' | sed -e 's/: .*//'"

foreach CANDIDATE ( $BRANCH1 $BRANCH2 )
   cd $CANDIDATE
   echo "Removing nuisance files   from $CANDIDATE"
   \rm `find . -name dart_log.nml`
   \rm `find . -name dart_log.out`
   \rm `find . -name input*default`
   echo "Removing DART executables from $CANDIDATE"
   \rm `findexec`
   cd -
end

# Want to compare the unversioned output files from running test_dart.csh
# The unversioned files have a '?' in column 1, and we need to remove the
# question mark and the blank spaces before the filename. We want FileList
# to contain the filenames relative to the base directory.

cd ${BRANCH1}
set FileList = `svn status | grep '^?' | sed -e "s#. *##"`
cd ${orgdir}

foreach FileName ( $FileList )    # loop through all the files

   # some versions of 'file' identify *.nml files as Sendmail 

   set OrgFile = ${BRANCH1}/${FileName}
   set NewFile = ${BRANCH2}/${FileName}

   set asciicheck = `(file ${OrgFile} | grep ASCII)`
   set isascii = $status
   set textcheck = `(file ${OrgFile} | grep " text")`
   set istext = $status
   set textcheck = `(file ${OrgFile} | grep " FORTRAN")`
   set isfortranish = $status
   set mailcheck = `(file ${OrgFile} | grep " Sendmail")`
   set isnamelist = $status

   if ($isascii == 0 || $istext == 0 || $isfortranish == 0 || $isnamelist == 0) then
      set isascii = 0
   else
      set isascii = 1
   endif

   switch (  $OrgFile:e )
   case 'nc'
      set netcdf = 1
      breaksw
   default:
      set netcdf = 0
      breaksw
   endsw

   if ( -f ${NewFile} && $isascii == 0 ) then ;# ASCII file exists

      set numdiffs = `cmp ${OrgFile} ${NewFile} | wc`
      if ( $numdiffs[1] > 0 ) then
         echo "${OrgFile} and ${NewFile} differ ... comparing ..."
         xxdiff ${OrgFile} ${NewFile}   
      else
         echo "${OrgFile} and ${NewFile} are identical ... nothing to do."
      endif

   else if ( $netcdf > 0 ) then ;# netCDF file, use compare_states to summarize
      echo ${OrgFile} ${NewFile} | ./compare_states

   else if ( -f ${NewFile} ) then ;# binary file exists, compare

      cmp ${OrgFile} ${NewFile}
      set DIFFSTAT = $status
      if ( $DIFFSTAT != 0 ) then
         echo "${OrgFile} and ${NewFile} differ ... non_ASCII"
      endif

   else
      echo "${NewFile} does not exist"
   endif

end

exit 0  

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
