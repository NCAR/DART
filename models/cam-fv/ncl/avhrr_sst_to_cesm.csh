#!/bin/csh -x
# #!/bin/csh -fvx    -f interferes with getting the module alias?

#PBS  -N avhrr_sst_to_cesm
#PBS  -A NCIS0006
#PBS  -q share
# Resources I want:
#    select=#nodes
#    ncpus=#CPUs/node
#    mpiprocs=#MPI_tasks/node
#PBS  -l select=1:ncpus=1:mpiprocs=1
#PBS  -l walltime=03:00:00
# Send email after a(bort) or e(nd)
#PBS  -m ae
#PBS  -M raeder@ucar.edu
# Send standard output and error to this file.
# It's helpful to use the $casename here.
#PBS  -o avhrr_sst_to_cesm.eo
#PBS  -j oe 
#--------------------------------------------

# Set the dates and submit this as an interactive job using
# > /glade/u/home/raeder/Scripts/qcmd_share_1cpu.csh ./avhrr_sst_to_cesm.csh 
# or if I need to log off soon
# > remote ./qcmd_avhrr.csh

# User set variables.
set product = avhrr

# These files have sst and sea ice in them.
set dir_in  = /gpfs/fs1/collections/rda/data/ds277.7/${product}/
# 2017/avhrr-only-v2.2017
set dir_out = /glade/scratch/raeder/SST/${product}

set dom_dir  = /gpfs/fs1/work/raeder/Models/CAM_init/SST
set dom_file = domain.ocn.d025.120821.nc

if ($#argv == 1) then
   set years = ($1)
else
   set years = ( 2018 2019 )
endif

cd $PBS_O_WORKDIR

echo "Starting at `date`"

# If not in my environment already:
module load gnu
module load ncl

set echo verbose
#----------------------------------------------------------------------------------------------------------------
if !(-e ${dir_out} ) then
   mkdir ${dir_out}
endif
cd ${dir_out}
#----------------------------------------------------------------------------------------------------------------
#untar
foreach yr ( $years )
#----------------------------------------------------------------------------------------------------------------
#    set tar_files = `ls -1 ${dir_in}/${product}_nc_v2_${yr}??.tar`
#    foreach tar_file ( ${tar_files} )
#       tar -xvf ${tar_file} 
   if (! -d $yr) then
      mkdir $yr
      cd $yr
      cp ${dir_in}/${yr}/${product}*${yr}* .
# Don't add days from bracketing years!  
#    Use multiple files in the stream file mechanism 
#    to handle the last day of each year.
#       @ yrp = $yr + 1
#       # * because some directories have compressed files, others don't
#       cp ${dir_in}/${yrp}/${product}*${yrp}0101.nc* .
#       @ yrm = $yr - 1
#       cp ${dir_in}/${yrm}/${product}*${yrm}1231.nc* .
      gunzip *.gz
   else
      cd $yr
# Temporary, for recovering from mistaken assumption about existence of directory.
      cp ${dir_in}/${yr}/${product}*${yr}* .
#       @ yrp = $yr + 1
#       cp ${dir_in}/${yrp}/${product}*${yrp}0101.nc* .
#       @ yrm = $yr - 1
#       cp ${dir_in}/${yrm}/${product}*${yrm}1231.nc* .
      gunzip *.gz
# end Temp
   endif
#----------------------------------------------------------------------------------------------------------------
#  Output file for first stage.
#  Eventual name will be something like
#  old; avhrr-only-v2.20100101_cat_20101231_filled_c130829.times
#   set file_out = ${product}-only-v2.${yr}0101_cat_${yr}1231_gregorian.nc
#  new(padding before and after this year; avhrr-only-v2.2010_filled_c130829.times
   set file_out = ${product}-only-v2.${yr}_gregorian.nc
   echo ${file_out}

# Want to use ncrcat.

# Change 'time' to a record (unlimited) dimension.
# -v sst,ice  = subset the variables
# (which it should have been to start with).
   foreach f (`ls -1 avhrr*.nc`)
      if (! -f temp_$f) then
         # Do first, to make zlev(zlev) a scalar var 
         # which will be removed by ncks.
         ncwa -O -a zlev -o nozlev.nc $f
         ncks --mk_rec_dmn time  -v sst,ice nozlev.nc temp_$f
      endif
   end
   rm nozlev.nc

#  Concatenate all of the daily files into a yearly file.
#  -O = overwrite existing output file
#  -h = don't propagate the history (I may want to)
#  Could make time the record dimension in the ncecat call; -u time
#    No, that caused an error.
#    Let ncecat make the file with the default record dimension 'record',
#    then fix it all.
   if (! -f ../temp_${file_out}) then
      ncrcat -O temp_* ../temp_${file_out} 
   else
      echo "`pwd`/../temp_${file_out} exists. Exiting"
      exit 5
   endif 

   cd ..

# Add vars from 1/4 degree domain file
# 	nj = 720 ;
# 	ni = 1440 ;
# 	int mask(nj, ni) ;
# 	double area(nj, ni) ;
   if (! -f ${dom_file}) then
      ncrename -d ni,lon -d nj,lat ${dom_dir}/${dom_file} ${dom_file} 
   endif
# -C; don't copy over the ancillary variables of mask, area (xc,xv,yc,yv)
# -A; append onto temp_${file_out}
   ncks -A -C -v mask,area ${dom_file} temp_${file_out}   
   
# Change variable precision from short to float
   ncap2 -s 'sst=float(sst);ice=float(ice);area=float(area)' temp_${file_out} ${file_out} 
   if ($status == 0) rm temp_${file_out}

# rename variables: sst -> SST_cpl, ice -> ice_cov
   ncrename -v sst,SST_cpl -v ice,ice_cov ${file_out} 

# Add 'calendar=gregorian attribute' to the time variable
   ncatted -a calendar,time,c,c,gregorian ${file_out} 

# The previous section seems to be very fast (2019-10-23)
#----------------------------------------------------------------------------------------------------------
# This section seems to be take ~45 minutes/year.

#  Need to add mask, area (using ./fill*?)
#  Actually, maybe I don't.  
#    2011 failed to do this, and shows no evidence in the file history that it was done,
#    but it was usable.

#    setenv CYEAR $yr
#    ncl /gpfs/fs1/work/raeder/Models/CAM_init/SST/Make_yearly_from_daily/fill_sst_forcing.ncl

#    cd ..
   
   echo "Finished $yr at `date`"

end # next year
