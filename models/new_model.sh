#!/bin/bash
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

function usage { echo "Usage: $0 [model_name] [location_choice]"; }

model_name=$1
location_choice=$2

if [[ ($# -le 1) ]]; then
  usage
  exit 0
elif [[ -d "$model_name" ]]; then
  echo "A directory called $model_name already exists! Please try a different model name."
  exit 1
elif [[ ! -d "../assimilation_code/location/$location_choice" ]]; then
  echo "$location_choice is not a valid location name! Please try a different location name."
  echo "Available location names are located in ../assimilation_code/location"
  exit 1
fi

mkdir -p "$model_name/work"
cp "template/new_model.rst" "$model_name/readme.rst"
sed -i '' -e "s;template_model;$model_name;" "$model_name/readme.rst" 2>/dev/null
cp "template/work/quickbuild.sh" "$model_name/work/quickbuild.sh"
sed -i '' -e  "s;template_model;$model_name;" -e "s;template_location;$location_mod;" "$model_name/work/quickbuild.sh" 2>/dev/null

case $location_choice in
  threed | threed_cartesian | threed_sphere) # threed does not compile, but is a good framework
    cp "template/threed_model_mod.f90" "$model_name/model_mod.f90"
    cp "template/work/threed_input.nml" "$model_name/work/input.nml"
    ;;
  oned)
    cp "template/oned_model_mod.f90" "$model_name/model_mod.f90"
    cp "template/work/oned_input.nml" "$model_name/work/input.nml"
    ;;
  twod)
    cp "template/twod_model_mod.f90" "$model_name/model_mod.f90"
    cp "template/work/twod_input.nml" "$model_name/work/input.nml"
    ;;
  *) # default is a threed_sphere template. As other location templates added, they can be added to case
    cp "template/threed_model_mod.f90" "$model_name/model_mod.f90"
    cp "template/work/threed_input.nml" "$model_name/work/input.nml"
    ;;
esac

echo -e "$model_name created successfully!\nTo get started,\n\tcd $model_name/work"
exit 0
