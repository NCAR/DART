#!/bin/csh

set DARTDIR=/glade/p/work/hendric/DART/johnny_rma_trunk

cd $DARTDIR


# Copy the original template to change back later
set MKMF_TEMP = "$DARTDIR/mkmf/mkmf.template"
cp $MKMF_TEMP $MKMF_TEMP.original

# Set todays date
set tds=`date +"%m.%d.%y"`

# Add build directory if it does not exist
if (! -d ${DARTDIR}/buildlogs) then
   mkdir ${DARTDIR}/buildlogs
endif 

# Make pass fail log 
set   PASS_FAIL="${DARTDIR}/buildlogs/PASS_FAIL.$tds.log"
echo "PASS FAIL LOG $tds \n" > $PASS_FAIL
echo "The top-level DART directory (DARTDIR) is $DARTDIR\n" | tee -a $PASS_FAIL
   
foreach COMPILER ( intel.linux gfortran pgi.linux) 
   # intel.linux      
   # pgi.cray
   # intel.osx        
   # pgi.linux
   # absoft.osx       
   # lahey.linux      
   # pgi.osx
   # g95              
   # original         
   # sgi.altix
   # gfortran         
   # pathscale.linux  
   # xlf.aix
   
   module purge

   switch ( $COMPILER )
   case intel.linux:
      module load intel
      breaksw
   case gfortran:
      module load gnu
      breaksw
   case pgi.linux:
      module load pgi
      breaksw
   default:
      echo "$COMPILER not supported"
      exit 0
      breaksw 
   endsw

   module load netcdf
   module load ncarcompilers
   module load ncarbinlibs
   module load ncarenv

   cp $MKMF_TEMP.${COMPILER} $MKMF_TEMP

   #----------------------------------------------------------------------
   # Compile 'filter' for a wide range of models.
   #----------------------------------------------------------------------
   set MPIFLAG = -mpi
   
   foreach MODEL ( POP wrf cam mpas_atm ) 
       
       # change directory to model work
       cd ${DARTDIR}/models/${MODEL}/work
   
       # create build log with time stamp
       set BUILDLOG="${DARTDIR}/buildlogs/${MODEL}.$tds.log"
       set BUILDERR="${DARTDIR}/buildlogs/${MODEL}.$tds.err"
   
       echo "${MODEL} BUILD LOG    $tds" > $BUILDLOG
       echo "${MODEL} BUILD ERRORS $tds" > $BUILDERR
   
       # header
       echo "Compiling $MODEL with $COMPILER\n" | tee -a $PASS_FAIL $BUILDLOG $BUILDERR
      
       # build and run preprocess 
       rm -f *.mod *.o preprocess Makefile
       (csh mkmf_preprocess >>& $BUILDLOG)  || (echo "preprocess MAKE FAILED" && break)
         
       echo "--------------------------------------" | tee -a $BUILDLOG $BUILDERR
       echo "building ${MODEL} preprocess"           | tee -a $BUILDLOG $BUILDERR
       echo "--------------------------------------" | tee -a $BUILDLOG $BUILDERR

       (make >>& $BUILDLOG)

       if( -e preprocess) then

          # assume the test is going to fail
          set pf="FAILED"
          # run preprocess
          (./preprocess >& /dev/null) && set pf="PASSED"
   
          printf "  %-25s %10s\n" "./preprocess" "$pf" | tee -a $PASS_FAIL
          
          if( $pf == "PASSED") then

             # make all programs in work directory
             foreach TARGET ( mkmf_* )
             
                set PROG = `echo $TARGET | sed -e 's#mkmf_##'`
             
                switch ( $TARGET )
                case mkmf_preprocess:
                   breaksw
                default:
                   rm -f ${PROG} *.o *.mod Makefile
                   (csh $TARGET -mpi >& $BUILDLOG)  || (echo "${PROG} MAKE FAILED" && breaksw)
   
                   echo "--------------------------------------" | tee -a $BUILDLOG $BUILDERR
                   echo "building ${MODEL} ${PROG}             " | tee -a $BUILDLOG $BUILDERR
                   echo "--------------------------------------" | tee -a $BUILDLOG $BUILDERR

                   (make >>& $BUILDLOG)
                   if( -e ${PROG}) then
                     printf "  %-25s %10s\n" "${PROG}" "PASSED" | tee -a $PASS_FAIL
                   else
                     printf "  %-25s %10s\n" "${PROG}" "FAILED" | tee -a $PASS_FAIL
                   endif
   
                   breaksw
                endsw
             end # programs
          endif # ./preprocess passed
       else
          echo "make preprocess FAILED             " >> $PASS_FAIL
          echo "aborting all other tests for $MODEL" >> $PASS_FAIL
       endif # make preprocess passed
   
   end # model
end # compiler

# copy the original template file back
echo "\nswitching back your original mkmf.template"

cp $MKMF_TEMP.original $MKMF_TEMP

