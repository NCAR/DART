#!/bin/tcsh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next three lines automatically updated by CVS, do not edit>
# $Id$
# $Source$
# $Name$
#
# Filter server waits for requests from DART clients to either advance an
# ensemble of model states or to assimilate a set of observations for a group
# of regions. Filter_server waits for one of three files to be created
# and takes the following actions:
# 1. go_end_filter: terminate filter_server.
# 2. go_advance_model: advance ensemble of model states.
# 3. go_assim_region: assimilate a set of regions.
#
#### LSF options for BSUB
### -J      job name
### -o      output listing filename 
### -q      queue
### Queue charging    cheapest ..... most expensive   (at least on lightning)
### Queue name (standby, economy, [regular,debug], premium)
#
#BSUB -J filter_server
#BSUB -o filter_server.%J.log
#BSUB -q debug
#BSUB -n 2

# If they exist, log the values of the batch environment vars.
if ($?LSB_JOBNAME)     then
   echo "LSB_JOBNAME     is $LSB_JOBNAME"
endif
if ($?LSB_JOBFILENAME) then
   echo "LSB_JOBFILENAME is $LSB_JOBFILENAME"
endif
if ($?LSB_MCPU_HOSTS)  then
   echo "LSB_MCPU_HOSTS  is $LSB_MCPU_HOSTS"
endif
if ($?LS_SUBCWD )      then
   echo "LS_SUBCWD       is $LS_SUBCWD"
endif
if ($?LSB_HOSTS)       then
   echo "LSB_HOSTS       is $LSB_HOSTS"
endif
if ($?LSB_EXECHOSTS)   then
   echo "LSB_EXECHOSTS   is $LSB_EXECHOSTS"
endif

# Determine number of processors -- one of three ways.
# 1) Batch jobs set a variable LSB_HOSTS
# 2) Interactive jobs can have a NPROCS environment variable defined.
# 3) Interactive jobs default to 1 (one).
#
# list of hosts/machines is in $PROCNAMES

if ($?LSB_HOSTS) then
   set NPROCS = `echo $LSB_HOSTS | wc -w`
   set PROCNAMES = "${LSB_HOSTS}"
else if ($?NPROCS) then
   set PROCNAMES = $host
   set iproc = 2
   while($iproc <= $NPROCS)
      set PROCNAMES = "$PROCNAMES $host"
      @ iproc ++
   end
else
   set NPROCS = 1
   set PROCNAMES = $host
endif

# The working directory is set one of two ways,
# 1) Batch jobs set $LS_SUBCWD that is the directory FROM WHICH 
#    you submitted the batch job.
# 2) the current working directory

if ( $?LS_SUBCWD ) then
   set DARTWORKDIR = $LS_SUBCWD
   cd $DARTWORKDIR
else
   setenv DARTWORKDIR `pwd`
endif

# Initialize a log of filter_server's actions
rm -f filter_server.log

### Output to confirm job characteristics
echo "filter_server Running on host "`hostname`   > filter_server.log
echo "Initialized at "`date`                     >> filter_server.log
echo "DARTWORKDIR is $DARTWORKDIR"               >> filter_server.log
echo "This job has allocated $NPROCS tasks."     >> filter_server.log
echo "They will run on the following nodes: "    >> filter_server.log
echo $PROCNAMES                                  >> filter_server.log


# Hang around forever for now and wait for go_advance file to appear
# The program 'filter' creates the file 'go_advance_file'.
while(1 == 1)

   #------------------------------------------------------------------------------------ 
   # If go_end_filter exists then stop this process
   #------------------------------------------------------------------------------------ 
   if( -e go_end_filter ) then
      echo "terminating normally at " `date` >> filter_server.log
      rm -f go_end_filter
      exit
   endif
 
   #------------------------------------------------------------------------------------ 
   # Check to see if the go_advance_model file exists
   #------------------------------------------------------------------------------------ 
   if( -e go_advance_model ) then

      # First line of filter_control should have number of model states to be integrated
      set nensmbl = `head -1 filter_control`
      echo "advancing $nensmbl members at " `date`>> filter_server.log
     
      # figure # batches of runs to do, from # ensemble members and # processors
      @ nbatch = $nensmbl / $NPROCS
      if ($nensmbl % $NPROCS != 0 ) @ nbatch++
      echo "$nbatch batches will be executed" >> filter_server.log

      # Create a directory for each member to run in for namelists
      set element = 0
      set batch = 1
      while($batch <= $nbatch)
         foreach proc ( $PROCNAMES )
            @ element++
            if ($element > $nensmbl) goto all_elements_done

            # set up a remote command and arguments -- easier to echo 
            set RCMD = $DARTWORKDIR/advance_model.csh
            set RCMDARGS = "$DARTWORKDIR $element ${USER}_tempdir${element}"
            echo "$proc csh $RCMD $RCMDARGS &"  >> filter_server.log
            ssh   $proc csh $RCMD $RCMDARGS &

         end
         echo "waiting to finish batch $batch at" `date`  >> filter_server.log
         wait

         @ batch++
      end
      all_elements_done:

      # Need to have all backgrounds completed before continuing
      wait

      # finished with advance_model so remove the go_advance_model file
      echo "Completed this advance at " `date` >> filter_server.log
      echo "---------" >> filter_server.log
      rm -f go_advance_model

   endif

   #------------------------------------------------------------------------------------ 
   # Check to see if the go_assim_regions file exists.
   # If it does, we can assimilate the regions.
   #------------------------------------------------------------------------------------ 
   if( -e go_assim_regions ) then

      # First line of assim_region_control should have number of regions to be assimilated
      set nregions = `head -1 assim_region_control`
      echo "assimilating $nregions regions at " `date`>> filter_server.log

      # figure # batches of runs to do, from # ensemble members and # processors
      @ nbatch = $nregions / $NPROCS
      if ($nregions % $NPROCS != 0 ) @ nbatch++
      echo "$nbatch batches will be executed" >> filter_server.log
      #
      # Create a directory for each member to run in for namelists
      set element = 0
      set batch = 1
      while($batch <= $nbatch)
         foreach proc ( $PROCNAMES )
            @ element++
            if ($element > $nregions) goto all_regions_done

            # set up a remote command and arguments -- easier to echo 
            set RCMD = $DARTWORKDIR/assim_region.csh
            set RCMDARGS = "$DARTWORKDIR $element ${USER}_tempdir${element}"
            echo "$proc csh $RCMD $RCMDARGS &"  >> filter_server.log
            ssh   $proc csh $RCMD $RCMDARGS &

         end
         echo "waiting to finish batch $batch" `date`  >> filter_server.log
         wait

         @ batch++
      end
      all_regions_done:

      # Need to have all backgrounds completed before continuing
      wait

      # signal to async_filter.csh to continue
      echo "Completed this assimilation at " `date` >> filter_server.log
      echo "---------" >> filter_server.log
      rm -f go_assim_regions

   endif

   #------------------------------------------------------------------------------------ 
   # No files found, wait and check again
   #------------------------------------------------------------------------------------ 
   sleep 1

end
