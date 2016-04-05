#!/bin/csh -f


setenv case                 cesm_6hr_popcam
setenv mach         yellowstone
setenv cesmroot     /glade/p/cesm/cseg/collections/$cesmtag
setenv caseroot     /glade/p/work/${USER}/cases/${case}
setenv exeroot      /glade/scratch/${USER}/${case}/bld
setenv rundir       /glade/scratch/${USER}/${case}/run
setenv archdir      /glade/scratch/${USER}/archive/${case}
setenv dartroot     /glade/u/home/${USER}/DART_svn/DART


set nonomatch       # suppress "rm" warnings if wildcard does not match anything

# The FORCE options are not optional.
# The VERBOSE options are useful for debugging though
# some systems don't like the -v option to any of the following
switch ("`hostname`")
   case be*:
      # NCAR "bluefire"
      set   MOVE = '/usr/local/bin/mv -fv'
      set   COPY = '/usr/local/bin/cp -fv --preserve=timestamps'
      set   LINK = '/usr/local/bin/ln -fvs'
      set REMOVE = '/usr/local/bin/rm -fr'

   breaksw
   default:
      # NERSC "hopper", NWSC "yellowstone"
      set   MOVE = '/bin/mv -fv'
      set   COPY = '/bin/cp -fv --preserve=timestamps'
      set   LINK = '/bin/ln -fvs'
      set REMOVE = '/bin/rm -fr'

   breaksw
endsw


cd ${caseroot}

./Buildconf/cam.user_nl.csh 
./Buildconf/pop2.user_nl.csh 

@ inst = 1
while ($inst <= $num_instances)

   # following the CESM strategy for 'inst_string'
   set inst_string = `printf _%04d $inst`

   # ===========================================================================
   set fname = "user_nl_cam${inst_string}"
   # ===========================================================================
   # For a HOP TEST ... empty_htapes = .false.
   # For a HOP TEST ... use a default fincl1
   # inithist == 'ENDOFRUN' ensures that CAM writes an initial file every time it stops.
   # fincl1,nhtfrq,mfilt all control the history file containing a REQUIRED PHIS field.
   #AK changed this
   echo " inithist      = 'ENDOFRUN'"                     >> ${fname}
   echo " ncdata        = 'cam_initial${inst_string}.nc'" >> ${fname}
   echo " empty_htapes  = .true. "                        >> ${fname}
   echo "fincl2        = 'UBOT','U200','U850','VBOT','V200','V850','TBOT','T200','T850','PRECT'"  >> ${fname}
   echo "fincl3        = 'PHIS','SST'" >> ${fname}
   echo " nhtfrq        = 0,-24,-6"                     >> ${fname}
   echo " mfilt         = 1,1,1"                             >> ${fname}
   echo " avgflag_pertape = 'A','A','I' "                         >> ${fname}   
   # ===========================================================================
##   set fname = "user_nl_clm${inst_string}"
   # ===========================================================================

   # Customize the land namelists
   # The filename is built using the REFCASE/REFDATE/REFTOD information.
   #
   # This is the time to consider how DART and CESM will interact.  If you intend
   # on assimilating flux tower observations (nominally at 30min intervals),
   # then it is required to create a .h1. file with the instantaneous flux
   # variables every 30 minutes. Despite being in a namelist, these values
   # HAVE NO EFFECT once CONTINUE_RUN = TRUE so now is the time to set these.
   #
   # DART's forward observation operators for these fluxes just reads them
   # from the .h1. file rather than trying to create them from the subset of
   # CLM variables that are available in the DART state vector.
   #
   # For a HOP TEST ... hist_empty_htapes = .false.
   # For a HOP TEST ... use a default hist_fincl1
   #
   # FIXME ... add documentation for configuring CLM history files

##   @ thirtymin = $assim_n * 2
#AK changed this becuase with clm bgc off there NEP is not a variable.
##   echo "hist_empty_htapes = .true."                 >> $fname
##   echo "hist_fincl1 = 'TSA'"                        >> $fname
###   echo "hist_fincl2 = 'NEP','FSH','EFLX_LH_TOT_R'"  >> $fname
##   echo "hist_nhtfrq = -$assim_n "                 >> $fname
##   echo "hist_mfilt  = 1 "                 >> $fname
##   echo "hist_avgflag_pertape = 'A' "             >> $fname

   # ===========================================================================
   set fname = "user_nl_pop2${inst_string}"
   # ===========================================================================

   # POP Namelists
   # init_ts_suboption = 'data_assim'   for non bit-for-bit restarting (assimilation mode)
   # init_ts_suboption = 'rest'         --> default behavior
   #
   # README:
   # Configuring the contents of the history files for POP is best explained in
   # the section marked "POP2: TAVG Settings" in the cesm1_1_1 pop2 namelist documentation
   # http://www.cesm.ucar.edu/models/cesm1.1/cesm/doc/modelnl/nl_pop2.html
   #
   # and the CESM-specific documentation for the tavg output variables in the pop2
   # online documentation:
   # http://www.cesm.ucar.edu/models/cesm1.1/pop2/doc/users/node78.html
   #
   # In CESM1_1_1 keep the values for tavg_file_freq_opt and tavg_freq_opt identical.
   # pop2/trunk_tags/cesm_pop_2_1_20130412  explains the issue.
   #
   # DEFAULT values for these are:
   # tavg_file_freq_opt = 'nmonth' 'nmonth' 'once'
   # tavg_freq_opt      = 'nmonth' 'nday'   'once'
   # The  first entry indicates we get a monthly average once a month.
   # The second entry indicates we get a monthly average as it is being created.
   # The  third entry indicates  we get a daily timeslice
   #
   # Default copies of SourceMods/src.pop2/ocn.*.tavg.csh files are provided in the
   # DART_SourceMods_cesm1_1_1.tar bundle.

   echo "init_ts_suboption  = 'data_assim'"             >> $fname
   echo "tavg_file_freq_opt = 'nmonth' 'never' 'never' " >> $fname
   echo "tavg_freq_opt = 'nmonth' 'never' 'never' " >>$fname
   # ===========================================================================
##   set fname = "user_nl_cice${inst_string}"
   # ===========================================================================
   # CICE Namelists

##   echo "ice_ic = '${run_refcase}.cice${inst_string}.r.${run_refdate}-${run_reftod}.nc'" >> $fname

   # ===========================================================================
##   set fname = "user_nl_rtm${inst_string}"
   # ===========================================================================
   # RIVER RUNOFF CAN START FROM AN OLD CLM RESTART FILE
   # you can specify the RTM filename here and override the settings from
   # RUN_REFCASE/RUN_REFDATE/RUN_REFTOD (something you cannot do with CLM).

##   echo "finidat_rtm = '${run_refcase}.rtm${inst_string}.r.${run_refdate}-${run_reftod}.nc'" >> $fname

   @ inst ++
end

./preview_namelists
