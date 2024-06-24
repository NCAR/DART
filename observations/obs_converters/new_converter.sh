#!/bin/bash
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

function usage {
  echo "Usage: $0 [converter_name] [location_mod] [data_format]"
}

converter_name=$1
location_mod=$2
data_format=$3

if [[ $# -le 2 ]]; then
  usage
  exit 0
elif [[ -d "$converter_name" ]]; then
  echo "$converter_name already exists! Please try a different converter name."
  exit 1
elif [[ ! -d "../../assimilation_code/location/$location_mod" ]]; then
  echo "$location_mod is not a proper location module! Please try a different location name."
  exit 1
fi

mkdir -p "$converter_name/work"
cp "template/new_converter.rst" "$converter_name/${converter_name}_to_obs.rst"
cp "template/new_converter.rst" "$converter_name/readme.rst"

sed -i -e "s;template_converter;$converter_name;" "$converter_name/${converter_name}_to_obs.rst" 2>/dev/null
sed -i -e "s;template_converter;$converter_name;" "$converter_name/readme.rst" 2>/dev/null

cp "template/work/quickbuild.sh" "$converter_name/work/quickbuild.sh"
chmod +x "$converter_name/work/quickbuild.sh"
sed -i -e "s;template_converter;$converter_name;" -e "s;template_location;$location_mod;" "$converter_name/work/quickbuild.sh" 2>/dev/null

case $location_mod in
  threed | threed_cartesian | threed_sphere)
    case $data_format in
      netcdf)
        cp "template/threed_converter_mod.f90" "$converter_name/${converter_name}_to_obs.f90"
        ;;
      *)
        echo "$data_format is not a valid data format! Please use a valid data format."
        exit 1
        ;;
    esac
    cp "template/work/threed_input.nml" "$converter_name/work/input.nml"
    ;;
  oned)
    case $data_format in
      netcdf)
        cp "template/oned_converter_mod.f90" "$converter_name/${converter_name}_to_obs.f90"
        ;;
      *)
        echo "$data_format is not a valid data format! Please use a valid data format."
        exit 1
        ;;
    esac
    cp "template/work/oned_input.nml" "$converter_name/work/input.nml"
    ;;
  *)
    echo "$location_mod is not a valid location module! Please use a valid location module."
    exit 1
    ;;
esac

echo -e "$converter_name created successfully!\nTo get started,\n  cd $converter_name/work"
exit 0