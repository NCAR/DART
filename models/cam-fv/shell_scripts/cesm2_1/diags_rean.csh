#!/bin/tcsh

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$
#
#PBS  -N diags_rean
#PBS  -A YOUR_ACCOUNT
#PBS  -q share
# Resources I want:
#    select=#nodes
#    ncpus=#CPUs/node
#    mpiprocs=#MPI_tasks/node
#PBS  -l select=1:ncpus=1:mpiprocs=1
#PBS  -l walltime=04:00:00
# Send email after a(bort) or e(nd)
#PBS  -m ae
#PBS  -M YOUR_EMAIL
# Send standard output and error to this file.
# It's helpful to use the $casename here.
#PBS  -o diags_rean.eo
#PBS  -j oe 
#PBS  -k eod 
# #PBS  -W depend=afterok:3047237.chadmin1.ib0.cheyenne.ucar.edu

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

# Comments to separate multiple job output that's sent to the same diags_rean.eo file.
echo "==========================================="
echo "Running obs_diag and tarring obs_seq.finals on `date`"

if (! -f data_scripts.csh) then
   echo "ERROR: no data_scripts.csh.  Submit from the CASEROOT directory"
   exit
endif

# Get CASE environment variables from the central variables file.
source ./data_scripts.csh
if ($status != 0) exit 57
# Convert to numbers
# @ num_month = `echo $data_month | bc`
# @ num_year  = `echo $data_year  | bc`
echo "data_month = $data_month"
echo "data_year  = $data_year"
echo "data_proj_space = ${data_proj_space}"
echo "data_DART_src   = ${data_DART_src}"
echo "data_DOUT_S_ROOT   = ${data_DOUT_S_ROOT}"
echo "data_CASEROOT   = ${data_CASEROOT}"
echo "data_CASE       = ${data_CASE}"

# Use big endian obs_diag for output from IBM
# set endian = '_big_endian'
set endian = ' '

set yymm = `printf %4d-%02d $data_year $data_month`

# ? Give this as an argument
set diag_dir = 	${data_DOUT_S_ROOT}/esp/hist/Diags_NTrS_${yymm}
set proj_dir = ${data_proj_space}/${data_CASE}/esp/hist/${yymm}
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
ls ../*obs_seq_final*${yymm}*[^rz] >! obs.list 
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

if ($data_month == 12) then
   @ year_last = $data_year + 1
   @ mo_last = 1
else
   @ year_last = $data_year
   @ mo_last = $data_month + 1
endif

ex input.nml<< ex_end
/obs_diag_nml/
/obs_sequence_name/
s;= '.*';= "";
/obs_sequence_list/
s;= '.*';= "./obs.list";
/first_bin_center/
s;=  BOGUS_YEAR, 1;=  $data_year,$data_month;
/last_bin_center/
s;=  BOGUS_YEAR, 2;=  $year_last,$mo_last;
wq
ex_end

if ($?LS_SUBCWD) then
else if ($?PBS_O_WORKDIR) then
else
   vi input.nml
endif

echo "Running ${data_DART_src}/models/cam-fv/work/obs_diag${endian}"
${data_DART_src}/models/cam-fv/work/obs_diag${endian} >&! obs_diag.out 
set ostat = $status
if ($ostat != 0) then
   echo "ERROR: obs_diag failed.  Exiting"
   exit 40
endif

# cd ../${obs_dir}
if (! -d ${proj_dir}) then
   mkdir ${proj_dir}
else
   echo "ERROR: $proj_dir exists."
   echo "       Refusing to overwrite"
   exit
endif

# Archive all the obs_seq_finals for this month in a tar file.
set obs_seq_tar = ${data_CASE}.cam_obs_seq_final.${yymm}.tgz 
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
# rm *obs_seq*${yymm}*

exit

