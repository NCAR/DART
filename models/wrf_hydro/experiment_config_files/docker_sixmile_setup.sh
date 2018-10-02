wrf_hydro_local=/Users/${USER}/WRF_Hydro/wrf_hydro_nwm_public
dart_local=/Users/${USER}/WRF_Hydro/wrf_hydro_dart
domain_local=/Users/${USER}/Downloads/sixmile
wrf_hydro_py_local=/Users/${USER}/WRF_Hydro/wrf_hydro_py

docker run -it \
       -v ${wrf_hydro_local}:/home/docker/wrf_hydro_nwm_public \
       -v ${dart_local}:/home/docker/wrf_hydro_dart \
       -v ${domain_local}:/home/docker/domain \
       -v ${wrf_hydro_py_local}:/home/docker/wrf_hydro_py \
       wrfhydro/dev:conda


## Inside docker 
