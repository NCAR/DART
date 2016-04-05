#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
#     This script requires Job Array syntax - available via LSF and PBS and ... 
#     Fundamentally - N identical jobs get submitted,
#     each with a unique value of 'mem_id'
#
# Function:
#
#     This file is a batch script that runs on coral to convert an
#     ensemble of NOGAPS output to a DART experiment directory called 'HOME'
#
#     The individual NOGAPS ensemble members are converted and staged into
#     HOME/mbr??? directories. The converters are MPI-aware, so this script
#     actually stages things and then fires off batch jobs to do the 
#     conversions.
#
# Usage notes:
#
#     LSF: the -J directive indirectly specifies the ensemble size.
#     PBS: the -t directive indirectly specifies the ensemble size.
#     All remaining settings come from config.csh
#
#=============================================================================
#
# This block of directives constitutes the preamble for the LSF queuing system 
#
# the normal way to submit to the queue is:    bsub < run_nogapsIC_to_dart.csh
#
# an explanation of the most common directives follows:
# -J    Job name
# -q    queue
# -W    hr:mn
# -n    number of processors (really)
# -o    STDOUT filename  ( %J == unique job id, %I == job array index )
# -e    STDERR filename
#
#=============================================================================
#
#BSUB -J mbr[1-10]
#BSUB -q economy
#BSUB -W 0:10
#BSUB -n 8
#BSUB -o convert.%J.%I.log
#BSUB -e convert.%J.%I.err
#
#=============================================================================
#
# This block of directives constitutes the preamble for the PBS queuing system 
#
# the normal way to submit to the queue is:    qsub run_nogapsIC_to_dart.csh
#
# an explanation of the most common directives follows:
# -N     Job name
# -t     Job array indices (specifies ensemble size in this case)
# -r n   Declare job non-rerunable
# -e <arg>  filename for standard error 
# -o <arg>  filename for standard out 
# -q <arg>   Queue name (small, medium, long, verylong)
# -l nodes=xx:ppn=2   requests two processors on the node.
#                     (ppn == Processors Per Node)
#
#=============================================================================
#
#PBS -N mbr
#PBS -t 1-20
#PBS -r n
#PBS -e convert.err
#PBS -o convert.log
#PBS -q medium
#PBS -l nodes=2:ppn=2
#
#=============================================================================

# Check for the existence of some variables to glean
# which queuing mechanism we are working with.
# Set 'queue-independent' variables for use for the remainder 
# of the script.

if ($?LSB_QUEUE) then

   #-------------------------------------------------------------------
   # This is used by LSF
   #-------------------------------------------------------------------

   setenv ORIGINALDIR $LS_SUBCWD 
   setenv JOBNAME     $LSB_OUTPUTFILE:ar
   setenv JOBID       $LSB_JOBID
   setenv MYQUEUE     $LSB_QUEUE
   setenv MYHOST      $LSB_SUB_HOST
   setenv mem_id      $LSB_JOBINDEX
   setenv f_mbr       $LSB_JOBINDEX_END

else if ($?PBS_QUEUE) then

   #-------------------------------------------------------------------
   # This is used by PBS - f_mbr cannot be set with PBS Job Array ...
   #-------------------------------------------------------------------

   setenv ORIGINALDIR $PBS_O_WORKDIR
   setenv JOBNAME     $PBS_JOBNAME
   setenv JOBID       $PBS_JOBID:ar
   setenv MYQUEUE     $PBS_QUEUE
   setenv MYHOST      $PBS_O_HOST
   setenv mem_id      $PBS_ARRAYID
   setenv f_mbr       xx

else

   #-------------------------------------------------------------------
   # You can run this interactively to check syntax, file motion, etc.
   # These are all just make-believe.
   #-------------------------------------------------------------------

   setenv ORIGINALDIR `pwd` 
   setenv JOBNAME     mbr001
   setenv JOBID       $$
   setenv MYQUEUE     Interactive
   setenv MYHOST      $host
   setenv mem_id      3
   setenv f_mbr       20

endif

#----------------------------------------------------------------------
# Echo job attributes ... nothing more.
#----------------------------------------------------------------------

echo
echo "${JOBNAME} ($JOBID) submitted   from $ORIGINALDIR"
echo "${JOBNAME} ($JOBID) submitted   from $MYHOST"
echo "${JOBNAME} ($JOBID) running in queue $MYQUEUE"
echo "${JOBNAME} ($JOBID) running       on $host"
echo "${JOBNAME} ($JOBID) member $mem_id of $f_mbr started at "`date`
echo

#----------------------------------------------------------------------
# Get the constants from the config.csh file that was present in the
# original submission directory so they can be kept consistent
# across all subsequent scripts.
#----------------------------------------------------------------------

source ${ORIGINALDIR}/config.csh

if ( $archive != "newton" ) then
    echo "ERROR You must enter newton as the archive machine"
    echo "ERROR ${ORIGINALDIR}/config.csh has it as $archive"
    exit 1
endif

if (! $?experiment_dir ) then
    echo "ERROR There is no viable experiment_dir variable"
    echo "ERROR This comes from ${ORIGINALDIR}/config.csh"
    exit 1
endif

echo "Staging the experiment in $experiment_dir"
echo "assimilation date is $dtg"

#########################################################
#
# In the experiment directory; create one subdirectory 
# for each ensemble member and move into it.
#
#########################################################

set index = `printf "%03d" $mem_id`
set dirname = "mbr$index"

echo "Creating member directory ${experiment_dir}/${dirname}"

mkdir -p ${experiment_dir}/${dirname}
mkdir -p ${experiment_dir}/${dirname}/bin
mkdir -p ${experiment_dir}/${dirname}/log
mkdir -p ${experiment_dir}/${dirname}/specfiles
mkdir -p ${experiment_dir}/${dirname}/outp

cd ${experiment_dir}/${dirname}

#########################################################
echo "putting the climo files into local storage"
#########################################################

# TJH climo is readonly ... could just use it in-situ ...
# In the batch array world -- all members are trying to do this.
mkdir -p ${climo}
ln -sfv ${climo_dir}/climo${resolution}/* ${climo}/

cp ${ocards_files}/idisis.th.car idisis                 || exit 1
cp ${ocards_files}/isisocd.jun05jul05.32.scaled isisocd || exit 1
cp ${climo}/noggeom${resolution}l${n_levels}_full.txt \
                                  specfiles/noggeom.txt || exit 1

#
#  stage the specfiles from a tar file to local storage
#
# ( note:  initial-time info assumed to already exist
#          on scratch from running of ET routine )
#

#FIXME: this is where the IC and non-IC scripts differ;
#   this one is set to use the perturb files for testing, but
#   really it should be picking up the analysis files from the
#   last model advance inside a single run, or else from where 
#   the initial conditions from a previous job are.  (the outer
#   run script should try to make these look as similar as possible).

echo "untarring specfiles for ensemble member ${mem_id}"
tar -xvf ${perturb_dir}/${dtg}.${dirname}.ic.specfiles.tar 

echo "Moving specfiles for ${dirname}"

mv ${dtg}/spect${resolution}/${dirname}/shist000000  specfiles/shist000000
mv ${dtg}/spect${resolution}/${dirname}/c3grid000000 specfiles/c3grid000000
mv ${dtg}/spect${resolution}/${dirname}/c3land000000 specfiles/c3land000000
mv ${dtg}/spect${resolution}/${dirname}/trpfil       specfiles/trpfil

#########################################################
# It seems that things had evolved to the point that
# HOME == TMP == temp and all of them were in common use.
# Since they all are (believe it or not) '.' ... there is
# no need to make things more complicated.  All of these are not needed.
#
# The data directory
# (location of read-only data like initial fields and climatology)
#
# set HOME = `pwd`
# set DATA = `pwd`/specfiles
# set GFLD = `pwd`/outp
# set flds = `pwd`/specfiles
#
# The working scratch directory.
# (where output data is written during the model run)
#
# set TMP = `pwd`
#
# X create a logical link to ${HOME}.  This is sometimes needed
# X if you use very long pathnames that exceed the size of the 
# X character variables that hold them in the Fortran code.
# X ln -fs ${HOME} temp
#   just use '.' 
#
#########################################################

cat <<EOF3 >! rdifil
 &rdilst
 iproc  = $iproc,
 jproc  = $jproc,
 jsplit = $jsplit,
 lgtrdy = f,
 lnmode = t,
 loutp  = t,
 lfcst  = t,
 lmpi2  = f,
 al2al  = t/
EOF3

cat <<EOF4 >! filist
 &namfil
 ocards  = 'isisocd',
 idfile  = 'idisis',
 ifilin  = './specfiles',
 ifilout = './outp',
 sstdir  = './f/f$1',
 icedir  = './f/f$1',
 snodir  = './f/f$1',
 pstdir0 = './',
 pstdir1 = './specfiles',
 hstdir  = './specfiles',
 clmdir  = '${climo}'/
EOF4

cat <<EOF5 >! namlsts
 &modlst
 dt      = 600.,
 dtrad   =  12.0,
 taui    =   0.0,
 taue    = $endtau,
 tauh    =   3.0,
 vistsh  = .0000354, .000707, 2.828,
 visvor  = .0000354, .000707, 2.828,
 visdiv  = .0000354, .000707, 2.828,
 dtgfnoc = '${dtg}',
 nmodev  = 3,
 jskip   = 2/
 &keylst
 lgeosfilt=t/
EOF5

echo "right before running ${NOGAPS2DART_exec_name:t}"

#---------------------------------------------------------
# convert the NOGAPS file to a DART start vector
#
# input.nml:&nogaps_to_dart_nml has the following settings: 
#
# nogaps_restart_filename    = 'specfiles/shist000000',
# nogaps_data_time_filename  = 'CRDATE.dat'
# nogaps_to_dart_output_file = 'dart_new_vector',
#---------------------------------------------------------

ln -sf specfiles/noggeom.txt .
cp ${ORIGINALDIR}/input.nml  . || exit 1
cp ${ORIGINALDIR}/config.csh . || exit 1

# CRDATE.dat is the model state time - 
# read by NOGAPS codes and the DART converters.

echo $dtg >! CRDATE.dat

${MPI} ${NOGAPS2DART_exec_name}

if ( $status != 0 ) then
   echo "ERROR ${NOGAPS2DART_exec_name:t} failed."
   exit 1
endif

#---------------------------------------------------------
# DART expects the states to be called filter_ic.NNNN 
# in the experiment directory. This is specified by 
# input.nml:restart_in_file_name = "filter_ic",
#---------------------------------------------------------

set icname   = `printf "filter_ic.%04d" $mem_id`

mv dart_new_vector ../${icname}

echo
echo "${JOBNAME} ($JOBID) job $mem_id of $f_mbr finished at "`date`

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

# original collaboration information - do not edit
# $orgURL: https://svn2.assembla.com/svn/ngdart/shell_scripts/run_nogapsIC_to_dart.csh $
# $orgId: run_nogapsIC_to_dart.csh 111 2010-06-09 21:55:44Z thoar $
# $orgRevision: 111 $
# $orgDate: 2010-06-09 15:55:44 -0600 (Wed, 09 Jun 2010) $
