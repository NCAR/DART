#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# $Id$

From chriss@ucar.edu Fri Dec 20 15:11:17 2002
Date: Thu, 19 Dec 2002 17:05:58 -0700 (MST)
From: Chris Snyder <chriss@ucar.edu>
To: Jeffrey Anderson <jla@ucar.edu>, Dale Barker <dmbarker@ucar.edu>
Subject: more scripts

Hi: Attached is the "remote" script that
 runs a forecast under rsh on some other
 machine.  --CS

------------------------------------------------
#csh -f
#
#set verbose

echo ' input arguments ' $argv

set LocalMachine    = $argv[1] ; shift
set RemoteMachine   = $argv[1] ; shift
set RemoteDir       = $argv[1] ; shift
set NumberOfMembers = $argv[1] ; shift
set NeededFiles     = ($argv[*])  

set Rundir = /mmmtmp/chriss/$LocalMachine
mkdir $Rundir

cd $Rundir
rm -f $LocalMachine.test

# construct list of members
set ensmems = () 
set iloop = 1 
while ( $iloop <= $NumberOfMembers )  
  set iens = `pad $iloop`
  set ensmems = ( $ensmems $iens ) ; @ iloop++
end

# loop through possible jobs

foreach iens ($ensmems[*])
    # try to start the job, if we can delete the file in "Need_to_run"
    # directory on the remote machine, then we can run it

  rsh $RemoteMachine "\rm $RemoteDir/Need_to_run/Member.$iens" > & $LocalMachine.test

    # if $LocalMachine.test is empty, then run the member, else go to the next

  if( -z $LocalMachine.test ) then  # run the member

    cd $Rundir
    mkdir Member.$iens
    cd Member.$iens

      echo " starting member " $iens " at " `date` \
           " on " $LocalMachine > $LocalMachine.status
      rcp  $LocalMachine.status \
           $RemoteMachine":"$RemoteDir/Running/Member.$iens.$LocalMachine

      # copy needed files from remote

        foreach file ($NeededFiles[*])
          rcp $RemoteMachine":"$RemoteDir/Member.$iens/$file $file
        end

      # run the job

        $NeededFiles[1] > &out.member.$iens  # note that the executable
                                             # must be the first file
                                             # in the NeededFiles list
        if($status != 0) then
          touch abort.member.$iens # check for abort
        endif

      # remove the original files

        rm -f $NeededFiles[*]
        set Files = `ls`  # remaining file, copy back to remote
        foreach file ($Files[*])
          echo " copying file " $file " to " $RemoteDir/Member.$iens/$file
          rcp $file $RemoteMachine":"$RemoteDir/Member.$iens/$file
        end

        echo " finished member " $iens " at " `date` \
             " on " $LocalMachine >> $LocalMachine.status

        rsh $RemoteMachine \
            "rm -f $RemoteMachine":"$RemoteDir/Running/Member.$iens.$LocalMachine"
        rcp  $LocalMachine.status \
             $RemoteMachine":"$RemoteDir/Finished/Member.$iens.$LocalMachine
        rcp  $LocalMachine.status \
             $RemoteMachine":"$RemoteDir/Member.$iens/Member.$iens.timing

      # clean

        rm -f * ; cd ../ ; rmdir Member.$iens
        rm -f $LocalMachine.test

  else

        rm -f $LocalMachine.test

  endif

end

cd $Rundir
cd ../
\rm -r $Rundir

exit

