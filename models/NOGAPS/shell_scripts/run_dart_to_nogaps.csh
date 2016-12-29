#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
#############################################################################
#
# Function:
#     This file is a batch script that runs on coral to convert a
#     set of DART restart files to a set of NOGAPS ...
#


# Set some default values.

set experiment_name       = ""

# base directories

set scratch_dir       = "/home/coral/nancy/nogaps/work"
set top_dir           = "/home/coral/nancy/nogaps"


# archive directories  

set archive              = "newton"
set archive_ens_dir      = "/home/coral/nancy/nogaps/outp"
set archive_climo_dir    = "/home/coral/hansenj"
set archive_analysis_dir = "/home/coral/hansenj/analyses_T47L42"


# NOGAPS executable

set resolution        = 47
set NOGAPS_exec_dir   = "~/subversion/virgin/models/NOGAPS/work"
set NOGAPS_exec_name  = dart_to_nogaps

# ending lead time and wall clock 

set endtau = '168'
set clock_time = '0:30'


# various MPI parameters

set iproc             = 1
set jproc             = 1
# jproc was 16
@ nproc = $iproc * $jproc
set jsplit=1
#set jsplit=8
set MPI='mpirun.lsf '




while ($#argv)

    switch ($argv[1])

    case -h*:
         echo 'run_dart_to_nogaps [-a archive] [-s scratch_dir] experiment_name'
         exit 1
    breaksw
             
    case -a*:
         shift
         if ($#argv < 1) then
         echo 'Usage: specify the archive machine -a'
         exit 1
         else
             set archive = $argv[1]
         endif
    breaksw

    case -p*:
         shift
         if ($#argv < 1) then
         echo 'Usage: specify the archive ensemble directory -p'
         exit 1
         else
             set archive_ens_dir = $argv[1]
         endif
    breaksw

    case -r*:
         shift
         if ($#argv < 1) then
         echo 'Usage: specify the resolution -r'
         exit 1
         else
             set resolution = $argv[1]
         endif
    breaksw

    case -iproc:
         shift
         if ($#argv < 1) then
         echo 'Usage: specify the number of processors in the i direction with -iproc'
         exit 1
         else
             set iproc = $argv[1]
         endif
    breaksw

    case -date:
         shift
         if ($#argv < 1) then
         echo 'Usage: specify date -date'
         exit 1
         else
             set dtg = $argv[1]
         endif
    breaksw

    case -jproc:
         shift
         if ($#argv < 1) then
         echo 'Usage: specify the number of processors in the j direction with -jproc'
         exit 1
         else
             set iproc = $argv[1]
         endif
    breaksw

    case -c*:
         shift
         if ($#argv < 1) then
         echo 'Usage: specify the archive climo directory -c'
         exit 1
         else
             set archive_climo_dir = $argv[1]
         endif
    breaksw

    case -s*:
         shift
         if ($#argv < 1) then
         echo 'Usage: specify the scratch directory with -a'
         exit 1
         else
             set scratch_dir = $argv[1]
         endif
    breaksw

    default:
         set experiment_name = $argv[1]
    breaksw

    endsw
    shift
end        

#########################################################
#
# Tests
# 1. See if the experiment directory exists.
# 2. See how many subdirectories exist for the experiment directory.
# 3. See if any initial conditions exist in the subdirectories.
#
#########################################################


if ($experiment_name == "") then
    echo 'You must enter an experiment name'
    exit 1
endif

if (! -d $scratch_dir) then
    echo 'You must enter a valid scratch base directory'
    exit 1
endif

echo "archive is $archive"

if ( $archive != "newton" ) then
    echo 'You must enter newton as the archive machine'
    exit 1
endif

cd $scratch_dir

if ($status == 1) then
    echo "Error: unable to change to the scratch directory $scratch_dir"
    exit 1
endif

echo "experiment_name is $experiment_name"
echo date is $dtg
set dtgp6 = `echo $dtg\ +6 | /home/coral/hansenj/bin/newdtg`
echo dtgp6 is $dtgp6 

# The top level directory is the experiment name

mkdir $experiment_name

#########################################################
#
#  retrieve the climo files and put them on scratch.
#
#########################################################

rm -rf ${scratch_dir}/climo$resolution
mkdir  ${scratch_dir}/climo$resolution
cp ${archive_climo_dir}/climo$resolution/* ${scratch_dir}climo$resolution/

#########################################################
#
#  retrieve the analysis files and put them on scratch.
#
#########################################################

cp ${archive_analysis_dir}/anal.${dtg}.tar $scratch_dir
tar -xvf ${scratch_dir}anal.${dtg}.tar
rm -f    ${scratch_dir}anal.${dtg}.tar
echo "just finished untar"
ls -l


mkdir    ${scratch_dir}/${experiment_name}/initial_shist

cp ${scratch_dir}/${dtg}/spect47/shist000000 ${scratch_dir}/${experiment_name}/initial_shist/



cd ${experiment_name}

mkdir data
mkdir queue

#########################################################
#
# Create one subdirectory for each ensemble member
# Determine the ensemble members by assuming the names
# are shist00n0000 where n is the ensemble number.
#
#########################################################

#--------------------------------------------------------
#  begin loop over ensemble members
#--------------------------------------------------------

set file = ( shist000000 )

# @ mem_id = $i_mbr  FIXME
@ mem_id = 1
@  i_mbr = 1
@  f_mbr = 1

while ( $mem_id <= $f_mbr ) 

    @ new_mem_id = ( $mem_id + 17 )

    set index = `echo ${file[${mem_id}]} | cut -c6-8`
    set dirname = "mbr$index"
    echo "Creating member directory ${dirname}"

# Make the ensemble member directories in the scratch directory.

    mkdir ${dirname}
    mkdir ${dirname}/bin
    mkdir ${dirname}/log
    mkdir ${dirname}/specfiles
    mkdir ${dirname}/outp

# Copy the input files to the specfiles directory.
# ( note:  initial-time info assumed to already exist
#          on scratch from running of ET routine )

    echo "ready to copy c3 and shist files, cwd is "`pwd`
    echo "dirname is ${dirname}"
    cp ${scratch_dir}/${experiment_name}/initial_shist/${file[${mem_id}]} ${dirname}/specfiles/shist000000
    cp ${scratch_dir}/${dtg}/spect47/c3grid000000 ${dirname}/specfiles/c3grid000000
    cp ${scratch_dir}/${dtg}/spect47/c3land000000 ${dirname}/specfiles/c3land000000
    cp ${scratch_dir}/${dtg}/spect47/trpfil       ${dirname}/specfiles/trpfil
    cp ${scratch_dir}climo47/noggeom.txt          ${dirname}/specfiles/noggeom.txt

    echo "cp ${scratch_dir}/${dtg}/spect47/c3grid000000 ${dirname}/specfiles/c3grid000000"
    echo "cp ${scratch_dir}/${dtg}/spect47/c3land000000 ${dirname}/specfiles/c3land000000"

    cp ${NOGAPS_exec_dir}/${NOGAPS_exec_name} ${dirname}/bin/
    cp ${NOGAPS_exec_dir}/dart_old_vector ${dirname}    # FIXME!
    set NOGAPS_exec_mpi = ${scratch_dir}/${experiment_name}/${dirname}/bin/${NOGAPS_exec_name}

    # the input.nml file with the dart namelist
    cp ${NOGAPS_exec_dir}/input.nml ${dirname}

#################################################    
#
# Now generate a submit script for the given ensemble member
#

    if (-d ${dirname}) then
         echo member_name is ${dirname}
    endif

        cat <<EOFINIT > "${dirname}/bin/ufcst$dirname"
#!/bin/csh 
#
#########################################################
#
#    This file is submitted to the Load Sharing Facility
#    schedular on the IBM SP3 at NAVO.
#
#########################################################

# BSUB -J m${index}
# BSUB -q standby
# BSUB -W $clock_time
# BSUB -n $nproc
# BSUB -o /home/coral/nancy/nogaps/log/ufcst${dirname}.out
# BSUB -e /home/coral/nancy/nogaps/log/ufcst${dirname}.err
# BSUB -m 'cr0139en'

########################################
#
#  Set some environment variables that will be global to all of the
#  ensemble member runs.
#
########################################

# The home directory
# (location of standard output and standard error files
#  created by the script submitted to the LoadLeveler)
#

  set HOME     = $scratch_dir$experiment_name/${dirname}

#
# terrain and climotological files
#
  set climo    = ${scratch_dir}climo$resolution
#
#
# The data directory
# (location of read-only data like initial fields and
#  climatology)
#

  set DATA=\$HOME/specfiles
  set GFLD=\$HOME/outp

#
# The working scratch directory.
# (where output data is written during the model run)
#

  set TMP=\$HOME

#
# The initial fields directory
#

  set flds=\$HOME/specfiles

#
#
########################################

set echo
echo $dtg >> \$TMP/CRDATE.dat

# create a logical link to \$TMP.  This is sometimes needed
# if you use very long pathnames that exceed the size of the 
# character variables that hold them in the Fortran code.

ln -fs \$TMP temp

######################################## 

cd \$HOME

cat <<EOF3 > rdifil
 &rdilst
 iproc = $iproc,
 jproc = $jproc,
 jsplit= $jsplit,
 lgtrdy=f,
 lnmode=t,
 loutp= t,
 lfcst= t,
 lmpi2=f,
 al2al=t/
EOF3

cat <<EOF4 > filist
 &namfil
 ocards='isisocd',
 idfile='idisis',
 ifilin='\$TMP/specfiles/',
 ifilout='\$TMP/outp/',
 sstdir='temp/f/f$1',
 icedir='temp/f/f$1',
 snodir='temp/f/f$1',
 pstdir0='\$TMP/',
 pstdir1='\$TMP/specfiles/',
 hstdir='\$TMP/specfiles/',
 clmdir='\$climo/'/
EOF4

cat <<EOF5 > namlsts
 &modlst
 dt=600.,
 dtrad= 12.0,
 taui=   0.0,
 taue= $endtau,
 tauh=   3.0,
 kdiff1= 2,
 kdiff2= 14,
 vistsh=.0000354,.000707,2.828,
 visvor=.0000354,.000707,2.828,
 visdiv=.0000354,.000707,2.828,
 dtgfnoc='$dtg',
 nmodev=3,
 jskip=2/
 &keylst
 lgeosfilt=t/
EOF5

cp /home/coral/hansenj/ocards_files/idisis.th.car idisis
cp /home/coral/hansenj/ocards_files/isisocd.jun05jul05.32.scaled isisocd

#---------------------------------------------------------
# run the model
#---------------------------------------------------------

$MPI $NOGAPS_exec_mpi

if ( \$status != 0 ) then
       echo status = \$status
       echo "nogaps forecast failed "
       exit
endif

# the current time and the advance_to time are now
# in the dart_date.time file here.

########################################
#
# cleanup processing
#
########################################

#--------------------------------------------------------------------------
# Remove this job from the queue directory
#--------------------------------------------------------------------------

cd ${scratch_dir}/${experiment_name}/queue
rm -rf ${dirname}

cd ${scratch_dir}/${experiment_name}/${dirname}/specfiles 

#--------------------------------------------------------------------------
# Rename the spectral history files to reflect the
# given ensemble member's designation 
#--------------------------------------------------------------------------

set mbrindex = $index

foreach histfile (shist*)
  set histindex = \`echo \$histfile | cut -c9-11\`
  echo mv \$histfile ${scratch_dir}/${experiment_name}/shist\${mbrindex}\${histindex}
end

foreach histfile (c3grid*)
  set histindex = \`echo \$histfile | cut -c10-12\`
  echo mv \$histfile ${scratch_dir}/${experiment_name}/c3grid\${mbrindex}\${histindex}
end

foreach histfile (c3land*)
  set histindex = \`echo \$histfile | cut -c10-12\`
  echo mv \$histfile $scratch_dir$experiment_name/c3land\$mbrindex\$histindex
end

echo mv trpfil $scratch_dir$experiment_name/trpfil\$mbrindex
echo mv hccnmpi $scratch_dir$experiment_name/hccnmpi\$mbrindex
cp noggeom.txt $scratch_dir$experiment_name/

cd $scratch_dir$experiment_name

#--------------------------------------------------------------------------
# Move the spectral history files and other output to the archive machine
#--------------------------------------------------------------------------

endif

EOFINIT

chmod 777 "${dirname}/bin/ufcst${dirname}"

if (! -e ${scratch_dir}${experiment_name}/${dirname}/bin/submit_script) then

    cat <<EOF >> ${scratch_dir}${experiment_name}/${dirname}/bin/submit_script 
#!/bin/csh 
# This script submits all the jobs to the batch system.
#
# Check and see how many jobs are already queued.
set count = `bjobs -u \$USER | fgrep \$USER | wc -l`
echo "There are \$count jobs already in the queue"
EOF
chmod +x ${scratch_dir}${experiment_name}/${dirname}/bin/submit_script 

endif

cat <<EOF >> ${scratch_dir}${experiment_name}/${dirname}/bin/submit_script 

bsub < ${scratch_dir}${experiment_name}/${dirname}/bin/ufcst${dirname} 

EOF

# Run the batch script to submit the jobs.

   cd ${scratch_dir}${experiment_name}/${dirname}
   bin/submit_script
   cd ${scratch_dir}${experiment_name}

#-------------------------------------------------------
# End of the while loop on the ensemble members 
#-------------------------------------------------------

   @ mem_id += 1

end

echo 'Done...'

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

# original collaboration information - do not edit
# $orgURL: https://svn2.assembla.com/svn/ngdart/shell_scripts/advance_model.csh $
# $orgId: run_dart_to_nogaps.csh 99 2010-06-03 21:08:12Z thoar $
# $orgRevision: 111 $
# $orgDate: 2010-06-09 15:55:44 -0600 (Wed, 09 Jun 2010) $

