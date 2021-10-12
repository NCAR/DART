#!/bin/bash

# L3 mapped
# daily, 9km/4km
# Modis Aqua / VIIRS chlorophyll-a OCI algorithm

work_dir=/Users/gharamti/Documents/DART_git/MITgcm-NBLING-DART/observations/obs_converters/ocean_color/work
data_dir=/Users/gharamti/Documents/DART_git/MITgcm-NBLING-DART/observations/obs_converters/ocean_color/data
mkdir -p $data_dir 

date1=2018-10-01
date2=2018-10-02
resol=4km #4km, 9km
obsname=chlor_a

obs_seq_tag=obs_seq_chl

minlat=10.0
maxlat=30.0
minlon=30.0
maxlon=50.0

crop_region=true

# MODIS/AQUA [2002]: NASA satellite orbitting earth from south to north; MODIS: Moderate Resolution Imaging Spectroradiometer
# VIIRS [2011]: MODIS successor and multi-disciplinary intrument; VIIRS: Visible and Infrared Imager/Radiometer Suite
instrument=VIIRS 

# webpage credentials
username=XXXXXXXX
password=XXXXXXXX
web_page=https://oceandata.sci.gsfc.nasa.gov/api/file_search

touch ~/.netrc

echo "machine urs.earthdata.nasa.gov login ${username} password ${password}" > ~/.netrc
chmod  0600 ~/.netrc

touch ~/.urs_cookies

echo -e ""
echo -e "Data type:                  : Ocean Color; $obsname"
echo -e "Date range                  : '$date1' -- '$date2'"
echo -e "Data will be downloaded from: $web_page"
echo -e "Ocean Color Instrument      : $instrument"
echo -e "Resolution                  : $resol \n"
echo -e "---------------------------------------------------------------------------------\n"

cd $data_dir

# Download data: 2 different instruments; MODIS AQUA & VIIRS
if [[ "$instrument" != "AQUA" && "$instrument" != "VIIRS" ]]; then 
   echo "Instrument: $instrument is not supported."
   exit
fi

inst=`echo $instrument | tr '[:upper:]' '[:lower:]'`

for file in `curl -d "sensor=$inst&sdate=$date1&edate=$date2&dtype=L3m&addurl=1&results_as_file=1&search=*DAY*$obsname*$resol*" $web_page | grep getfile`; do
    curl -L -O -n -c ~/.urs_cookies -b ~/.urs_cookies $file
    files+=(`echo $file | rev | cut -d'/' -f1-1 | rev`)
done 

# Modifications to the netcdf file 
echo -e "\nPreparing the netcdf files for conversion .."

cp ${work_dir}/input.nml ${data_dir}
ln -sf ${work_dir}/advance_time .

for k in `seq 1 ${#files[@]}`; do
    
    ncfile=${files[`expr $k - 1`]}
 
    year=`echo $ncfile | cut -c 2-5`
    days=`echo $ncfile | cut -c 6-8`
    pryr=`expr $year - 1`

    greg_days=(`echo ${pryr}1231 ${days}d -g | ./advance_time`)

    if $crop_region; then 
       # Red Sea region
       ncks -O -d lat,$minlat,$maxlat -d lon,$minlon,$maxlon -v chlor_a $ncfile tmp.nc
       mv tmp.nc $ncfile
    fi

    # Add the time dimension
    ncap2 -O -s 'defdim("time",1)' $ncfile out.nc

    # Now, add the time variable
    ncap2 -O -s "time[time]=${greg_days[0]}.0" out.nc out2.nc

    # Append the rest of the time attributes
    ncap2 -O -s 'time@longname="Time";time@units="days since 1601-01-01 00:00:00";time@calendar="gregorian"' out2.nc $ncfile
done


# Run the converter
echo -e "\nRunning the ocean color converter .."

ln -sf ${work_dir}/convert_sat_chl .
sed -i -e "s/.*debug.*/   debug           = .false./" input.nml

for k in `seq 1 ${#files[@]}`; do
    
    ncfile=${files[`expr $k - 1`]}

    echo -e "\nRaw obs NETCDF file: $ncfile"

    # Update the input/output filename
    sed -i' ' -e "s/.*file_in.*/   file_in         = '$ncfile'/" input.nml
    sed -i' ' -e "s/.*file_out.*/   file_out        = '$obs_seq_tag'/" input.nml

    # Run the converter
    ./convert_sat_chl
done

# cleaning
rm input.nml*
rm advance_time
rm convert_sat_chl
rm dart_log* out*

echo -e "Done.\n"
