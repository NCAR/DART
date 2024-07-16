#!/bin/bash
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

function usage {
  echo "Usage: $0 [converter_name] [data_format]"
  echo "Valid data formats: netcdf, hdf, hdf5, csv, txt"
}

converter_name=$1
data_format=$2

if [[ $# -le 1 ]]; then
  usage
  exit 0
elif [[ -d "$converter_name" ]]; then
  echo "$converter_name already exists! Please try a different converter name."
  exit 1
fi

mkdir -p "$converter_name/data"
mkdir -p "$converter_name/shell_scripts"
mkdir -p "$converter_name/work"

cp "template/new_converter.rst" "$converter_name/${converter_name}_to_obs.rst"
cp "template/new_converter.rst" "$converter_name/readme.rst"

sed -i -e "s;template_converter;$converter_name;" "$converter_name/${converter_name}_to_obs.rst" 2>/dev/null
sed -i -e "s;template_converter;$converter_name;" "$converter_name/readme.rst" 2>/dev/null

cp "template/work/quickbuild.sh" "$converter_name/work/quickbuild.sh"
chmod +x "$converter_name/work/quickbuild.sh"
sed -i -e "s;template_converter;$converter_name;" -e "s;template_location;threed_sphere;" "$converter_name/work/quickbuild.sh" 2>/dev/null

case $data_format in
  netcdf)
    cp "template/threed_sphere_netcdf_converter_mod.f90" "$converter_name/${converter_name}_to_obs.f90"
    ;;
  hdf)
    cp "template/threed_sphere_hdf_converter_mod.f90" "$converter_name/${converter_name}_to_obs.f90"
    ;;
  hdf5)
    cp "template/threed_sphere_hdf5_converter_mod.f90" "$converter_name/${converter_name}_to_obs.f90"
    ;;
  csv)
    cp "template/threed_sphere_csv_converter_mod.f90" "$converter_name/${converter_name}_to_obs.f90"
    ;;
  txt)
    cp "template/text/text_to_obs.f90" "$converter_name/${converter_name}_to_obs.f90"
    cp "template/text/text_to_obs.rst" "$converter_name/${converter_name}_to_obs.rst"
    cp -r "template/text/data" "$converter_name"
    cp -r "template/text/shell_scripts" "$converter_name"
    cp -r "template/text/work" "$converter_name"
    ;;
  *)
    echo "$data_format is not a valid data format! Please use netcdf, hdf, hdf5, csv, or txt."
    exit 1
    ;;
esac

cp "template/work/threed_input.nml" "$converter_name/work/input.nml"

echo -e "$converter_name created successfully!\nTo get started,\n  cd $converter_name/work"
exit 0