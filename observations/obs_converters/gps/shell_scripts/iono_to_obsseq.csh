#!/bin/csh -v
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#

# makes multiple obs_seq files, each containing observations from
# all the ionospheric observations (ionPrf) for a day
#
# downloads the daily tar file from the CDAAC (set your user name
# and password below), untars it into 100s of individual netcdf files,
# one profile per file.
#
# runs the converter to create a single obs_seq.ion.YYYYMMDD file
# with all the obs for that day.
#


# set start and end YYYYMMDD here. 
set startd = 20130311
set   endd = 20130312

# set these to download from the CDAAC repo
set cdaac_uname = your_user_name
set cdaac_pword = your_password

# web address of the data repo
set repo_base = http://cdaac-www.cosmic.ucar.edu/cdaac/rest/tarservice/data

# set once and should be able to leave as-is
set input_dir    = ../raw
set output_dir   = ../daily
set advance_exec = ./advance_time
set nml_template = ./input.nml.template
set COPY         = "/bin/cp -f"
set MOVE         = "/bin/mv -f"


# should be ok from here down

set cycle_interval = 1d 

# loop from start to end time.
set dtg = $startd
while ( $dtg <= $endd )
  
  set julian=(`echo $dtg 0 -j | $advance_exec`)
  set jul_year = $julian[1]
  set jul_doy  = `printf %03d $julian[2]`

  # cut off hours 
  set ymd = `echo $dtg | cut -c 1-8`

  # get year, month, date for use below
  set y = `echo $dtg | cut -c 1-4`
  set m = `echo $dtg | cut -c 5-6`
  set d = `echo $dtg | cut -c 7-8`

  # set the format of the input and output files, tars, dirs, etc
  set input_doydir =  $input_dir/cosmic2013/ionPrf/${jul_year}.{$jul_doy}
  set input_tar    =             ${jul_year}.{$jul_doy}.tar

  set output_file  = \'$output_dir/obs_seq.ion.${y}-${m}-${d}\'

  # download from CDAAC server
  wget --http-user=$cdaac_uname --http-passwd=$cdaac_pword $repo_base/cosmic2013/ionPrf/${jul_year}.${jul_doy} -O $input_dir/$input_tar

  # untar
  cd $input_dir
  tar -xvf $input_tar
  cd -

  # list all files in the directory into a fixed filename
  ls $input_doydir/*_nc >! file_list.txt    

  echo input will come from directory $input_doydir
  echo output will go into file $output_file

  # set up the conversion parameters in the namelist

  sed -e "s;^.*output_file *=.*;output_file=\ ${output_file};" \
      $nml_template > input.nml      

  # do the conversion here
  ./convert_cosmic_ionosphere

  # advance to next time
  set dtg = `echo $dtg ${cycle_interval} | $advance_exec`

end

echo 'Finished'

exit 0

