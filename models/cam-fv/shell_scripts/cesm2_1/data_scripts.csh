#! /bin/tcsh

# This script should be created by setup_advanced,
# probably as a 'here' document.
# It will be run by each script that needs any of this data.
# !  Be careful about where it is run; don't let it overwrite arguments
#    which are passed in.


The values here will be provided by setup_advanced and xmlquery.
setenv data_CASEROOT     ${CASEROOT}
setenv data_CASE         ${CASEROOT}:t
setenv data_DART_src     ~/DART/reanalysis_git
setenv data_CESM_python  ${CIMEROOT}/scripts/lib/CIME
setenv data_scratch      /glade/scratch/${USER}/${CASE}
setenv data_proj_space   /glade/p/nsc/ncis0006/Reanalyses/${CASE}
setenv data_campaign     /gpfs/csfs1/cisl/dares/Reanalyses
setenv data_NINST        ${NINST}


Do we need or want all 4 date parts in here?
   As a single string or 4 parts?
   I think just the year and month, which are used in directory names
     and a few other ways (obs_diag).
   More specific dates will be provided as arguments to the scripts.
! Carefully escape the special characters.

set CONTINUE_RUN = `./xlmquery CONTINUE_RUN --value`
if ($CONTINUE_RUN == FALSE) then
   set START_DATE = `./xlmquery CONTINUE_RUN --value`
   set parts = `echo $START_DATE | sed -e "s#-# #"`
   setenv data_year $parts[1]
   setenv data_month $parts[2]

else if ($CONTINUE_RUN == FALSE) then
   # Get date from an rpointer file
   set FILE = `head -n 1 ${data_scratch}/run/rpointer.atm_0001`
   set FILE = $FILE:r
   set ATM_DATE_EXT = $FILE:e
   set ATM_DATE     = `echo $FILE:e | sed -e "s#-# #g"`
   setenv data_year   `echo $ATM_DATE[1] | bc`
   setenv data_month  `echo $ATM_DATE[2] | bc`

else
   echo "env_run.xml: CONTINUE_RUN must be FALSE or TRUE (case sensitive)"
   exit

endif
