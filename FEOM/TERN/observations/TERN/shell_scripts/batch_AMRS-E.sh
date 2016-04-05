#!/bin/bash
#PBS -l ncpus=1
#PBS -l walltime=06:00:00
#PBS -l mem=500mb
#PBS -P xa5
#PBS -q express

#module load gdal python netcdf
module use /projects/xa5/modules
module load pythonlib/netCDF4/1.0.4
module load gdal

python /g/data1/xa5/pyeMAST/convert2netCDF.py
