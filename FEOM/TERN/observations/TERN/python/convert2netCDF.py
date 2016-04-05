"""
Created on Thu Jan 16 13:11:10 2014

@author: Daisy Duursma
Purpose of code is to take ascii data and convert it to NetCDF-CF
Built in data checks to make sure:
    lat and longs are as expected
    no data value is -9999
    ****to do: 
        check to see if data is uniform and resample if cells are not in desired locations
        read in standardized metadata and add data to netCDF header
        
    

"""

import numpy as np
import gdal
import netCDF4
from osgeo import gdalconst



def split_asciigrid(f_name):
    """
    Split an Arc ASCII grid file into the seperate elements. 
        -lats
        -longs
        -data
    Check that data is as expected and nodata = -9999
    """
    #open ascii file
    ds =  gdal.Open(f_name, gdalconst.GA_ReadOnly)
    band = ds.GetRasterBand(1)
    #get array of data values
    a = ds.ReadAsArray()
    #make sure no data is  -9999
    nodata = band.GetNoDataValue()
    #if nodata does not = -9999 replace nodata value with -9999
    if nodata!=-9999:
        a[a==nodata] = -9999
    #get length of lat and long    
    nlat,nlon = np.shape(a)
    #make tuples of lat and long
    b = ds.GetGeoTransform()
    lon = np.arange(nlon)*b[1]+b[0]
    lat = np.arange(nlat)*b[5]+b[3]
    #return errors if shape of data is not correct. 
    cols = ds.RasterXSize
    rows = ds.RasterYSize
    if len(a)!=rows:
        raise ValueError("Array contains data but not of length "+str(rows))
    if len(a[1])!=cols:
        raise ValueError("Array contains data but not of length "+str(cols))
    if len(a[1])!=len(lon):
        raise ValueError("Longitude vector contains data but not of length "+str(cols))
    if len(a)!=len(lat):
        raise ValueError("Latitude vector contains data but not of length "+str(cols))
    #return items needed
    return lon,lat,nlat,nlon,a,b,nodata
    




def make_netCDF(lond,latd,nlat,nlon,a,b,nodata,f_tm,out_file):
    """
    make a new netCDF file and write the lat,long, and data to file
    """
    #make empty netCDF to serve as container for variables, dimensions and attributes
    nco = netCDF4.Dataset(out_file,mode = 'w', clobber = True, format = 'NETCDF4_CLASSIC')
    
    # defines the sizes of all variables in terms of dimensions
    nco.createDimension('time',1)
    nco.createDimension('lon',nlon)
    nco.createDimension('lat',nlat)
    #create the variable names, datatype and dimensions, tuple containing the name from the dimensions
    lon = nco.createVariable('lon','f4',('lon',))
    lat = nco.createVariable('lat','f4',('lat',))
    time = nco.createVariable('time','i2',('time',))
    layer = nco.createVariable(f_stnd_nm, 'f4',  ('time','lat', 'lon',),fill_value=nodata)
    # add the variable attributes
    lon.units = 'degrees'
    lon.standard_name = 'longitude'
    lon.long_name = 'longitude'
    lon._CoordinateAxisType = 'Lon'
    lat.units = 'degrees'
    lat.standard_name = 'latitude'
    lat.long_name = 'latitude'
    lat._CoordinateAxisType = 'Lat'
    time.unit = "%s of %s" %(f_tm_unit,f_tm)
    time.standard_name = 'time'
    time.long_name = 'time'
    time._CoordinateAxisType = 'Time'
    time.calander = 'gregorian'
    layer.units = f_unit
    #layer._FillValue = -9999.
    layer.long_name = f_lng_nm
    layer.standard_name  = f_stnd_nm
    # add the global attributes
    nco.Conventions = 'CF-1.6'
    nco.spatial_coverage = "Australia"
    nco.long_name = f_lng_nm
    nco.geospatial_lat_min = min(lat)
    nco.geospatial_lat_max = max(lat)
    nco.geospatial_lat_units = 'degrees'
    nco.geospatial_lat_resolution = b[1]
    nco.geospatial_lon_min = min(lon)
    nco.geospatial_lon_max = max(lon)
    nco.geospatial_lon_units = 'degrees'
    nco.geospatial_lon_resolution = b[1]
    nco.crs_name = 'Lon/Lat Coords in WGS84'
    nco.crs_longitude_of_prime_meridian = 0.0
    nco.crs_semi_major_axis = 6378137.0
    nco.crs_inverse_flattening = 298.257223563
    #write lon,lat, and data
    lon[:]=lond
    lat[:]=latd
    layer[:] = a
    #close netcdf
    nco.close()
    
    
###########################################



f_tm_unit = 'day'
#unit of data, this will come out of the metadata.csv
f_unit = 'vol/vol'
#long name of data, to come out of metadata.csv
f_lng_nm = 'AMSR-E Soil moisture (VUA-NASA_LPRM)'
#standard name of data, to come out of metadata.csv
f_stnd_nm = 'AMSR-E_SM'
# Raijin input directory
iDir ="/g/data1/xa5/fenner/pyeMAST"
# Raijin output directory
oDir ="/g/data1/xa5/fenner/pyeMAST"
#
# import libraries used only for MAIN
import os
# change the working directory
os.chdir(iDir)
# loop through all files in the working directory
for r,d,f in os.walk('.'):
    # loop through the files one by one
    for files in f:
        # work only on the .flt files i.e. ignore the header files
        if files.endswith('.flt')==1:
                # generate a tidy output file name
                out_file = '%s/%s.nc' %(oDir,str.split(files[0:20])[0])
                # extract the time from the exitsitng file name, set it var.                
                f_tm = str.split(files[12:20])[0]
                # print for debug
                  #print 'processing ... %s ...  f_tm=' % (out_file,f_tm)
                # extract data from original file - one by one                 
                lond,latd,nlat,nlon,a,b,nodata = split_asciigrid(files)
                # make he new netcdf file
                make_netCDF(lond,latd,nlat,nlon,a,b,nodata,f_tm,out_file)
