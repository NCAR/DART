#!/bin/tcsh

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$
#
#PBS  -N diags_rean
#PBS  -A NCIS0006
#PBS  -q share
# Resources I want:
#    select=#nodes
#    ncpus=#CPUs/node
#    mpiprocs=#MPI_tasks/node
#PBS  -l select=1:ncpus=1:mpiprocs=1
#PBS  -l walltime=04:00:00
# Send email after a(bort) or e(nd)
#PBS  -m ae
#PBS  -M raeder@ucar.edu
# Send standard output and error to this file.
# It's helpful to use the $casename here.
#PBS  -o diags_rean.eo
#PBS  -j oe 
#--------------------------------------------

# obs_diag can now take a namelist argument that specifies a file
# containing a list of obs_seq filenames.
if ($?LS_SUBCWD) then
   cd $LS_SUBCWD
else if ($?PBS_O_WORKDIR) then
   cd $PBS_O_WORKDIR
else
#    if ($#argv == 0) then
#       echo 'Two scripts are required to generate observation space diagnostics '  
#       echo 'from obs_seq.[dates].final files; this one and ~/Scripts/matlab_cesm.csh.'  
#       echo '1) Go to the $EXEDIR directory, which contains the run directory.'  
#       echo '2) If needed, make a directory, e.g. "Obs_seq_final", to hold the obs_seq.final files.'  
#       echo '3) Link (or copy: >hsi  >cd hsi_dir  >get obs*)   the obs_seq.final files into Obs_seq_final.'  
#       echo '4) In EXEDIR run this script, giving it the directory where'  
#       echo '   the diagnostics will be created and Obs_seq_final (no ../):'  
#       echo '    > > >  ~/Scripts/diags_cesm.csh Diag_[details]_[dates] Obs_seq_final  < < < '
#       echo '5) Go to the diag_[dates] directory. '  
#       echo '6) Execute ~/Scripts/matlab_cesm.csh'  
#       exit
#    endif
endif

# DART source directory on this machine
# These things should be gathered from env_*.xml files in CASEROOT.
set DART = ~/DART/reanalysis
set cam = 'cam-fv'

# Use big endian obs_diag for output from IBM
# set endian = '_big_endian'
set endian = ' '
set year = 2013
# set diag_dir = 	Diags_NTrS_${year}-${mm}
set mo   = 4

set mm = `printf %02d $mo`
set yymm = ${year}-${mm}
set proj_dir = /glade/p/nsc/ncis0006/Reanalyses/f.e21.FHIST_BGC.f09_025.CAM6assim.011/esp/hist/$yymm

set diag_dir = 	Diags_N.PaAmAt_${year}-${mm}
echo "diag_dir = $diag_dir"
echo "proj_dir = $proj_dir"


if (! -d $diag_dir) then
   mkdir $diag_dir
   cd $diag_dir
   echo "In $diag_dir"
else
   echo "$diag_dir exists; choose another name"
   exit 20
endif

pwd

# ls -1 does not work; unusable formatting.
ls ../*obs_seq_final*${year}-${mm}*[^rz] >! obs.list 
if ($status != 0) then
   echo "Making obs.list failed.  Exiting"
   exit
endif


if (-e ../input.nml) then
   echo "using ./input.nml" 
   cp ../input.nml .
else
   echo "can't find an input.nml file"
   exit 30
endif

if ($mo == 12) then
   @ year_last = $year + 1
   @ mo_last = 1
else
   @ year_last = $year
   @ mo_last = $mo + 1
endif

ex input.nml<< ex_end
/obs_diag_nml/
/obs_sequence_name/
s;= '.*';= "";
/obs_sequence_list/
s;= '.*';= "./obs.list";
/first_bin_center/
s;=  $year, 1;=  $year,$mo;
/last_bin_center/
s;=  $year, 2;=  $year_last,$mo_last;
wq
ex_end

if ($?LS_SUBCWD) then
else if ($?PBS_O_WORKDIR) then
else
   vi input.nml
endif

echo "Running ${DART}/models/${cam}/work/obs_diag${endian}"
${DART}/models/${cam}/work/obs_diag${endian} >&! obs_diag.out 
set ostat = $status
if ($ostat != 0) then
   echo "ERROR: obs_diag failed.  Exiting"
   exit 40
endif

# Create the obs_seq_tar file name
# Extract the case name from the first file name in the obs.list.
# OR; set CASE = `./xmlquery CASE --value`
# ${CASE}.dart.e.cam_obs_seq_final.YYYY-MM-DD-SSSSS
set obs_seq = `head -n 1 obs.list` 
# Extract the file name from the path (:t), then strip off 
# everything after the $CASENAME (:r), and finally add the pieces
# needed for the tar file name.
set obs_seq_tar = $obs_seq:t:r:r:r:r.cam_obs_seq_final.${yymm}.tgz 

# cd ../${obs_dir}
if (! -d ${proj_dir}) then
   mkdir ${proj_dir}
else
   echo "ERROR: $proj_dir exists."
   echo "       Refusing to overwrite"
   exit
endif

tar -c -z -f ${proj_dir}/$obs_seq_tar ../*obs_seq*${yymm}*
if ($status != 0) then
   echo "ERROR: tar of obs_seq_finals failed.  Exiting"
   exit 50
endif

cd ..
if (-f $obs_seq_tar) then
   echo "obs_seq tar file was about to be removed.  Exiting"
   exit 60
endif
rm *obs_seq*${year}-${mm}*

exit

# <next few lines under version control, do not edit>
# $URL$
# $Id$
# $Revision$
# $Date$
