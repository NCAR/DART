import rioxarray as rxr
import xarray as xr
import geopandas as gpd
import numpy as np
import sys 

hdf_file = sys.argv[1]
year = sys.argv[2]
month = sys.argv[3]
day = sys.argv[4]

pre_file = rxr.open_rasterio(hdf_file)

# translate to geopandas to extract lat/lon from the hdf file's EASE-Grid CRS
mid_file = gpd.GeoDataFrame(data = pre_file.Ice_Surface_Temperature_NP[0], 
                            geometry=gpd.points_from_xy(pre_file.x, pre_file.y),
                            crs='EPSG:3408')

laty = mid_file.to_crs('EPSG:4326').geometry.y
lonx = mid_file.to_crs('EPSG:4326').geometry.x

laty = laty.where((laty <= 90) & (laty >= 0))
lonx = lonx.where((lonx <= 180) & (lonx >= -180))

lon, lat = np.meshgrid(lonx, laty)

pre_file['tsfc'] = pre_file.Ice_Surface_Temperature_NP[0]/100 - 273.15
pre_file['lat'] = (('y','x'), lat)
pre_file['lon'] = (('y','x'), lon)

post_file = pre_file.drop_vars({'Sea_Ice_by_Reflectance_NP', 
                                'Ice_Surface_Temperature_NP',
                                'Sea_Ice_by_Reflectance_SP',
                                'Ice_Surface_Temperature_SP'
                                })

post_file = post_file.drop_dims('band')
post_file = post_file.drop_vars('spatial_ref')

post_file = post_file.where(post_file.tsfc > -100, drop=True)
post_file = post_file.fillna(-800)

# post_file.to_netcdf('MOD29E1D_'+year+'_'+month+'_'+day+'_preprocessed.nc')
post_file.to_netcdf('test.nc')
