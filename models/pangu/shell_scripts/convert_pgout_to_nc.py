import sys
import xarray as xr
import numpy as np
from netCDF4 import Dataset,num2date,date2num
from datetime import datetime

landmask = np.load('/glade/work/chennuo/data/global_land_data/landmask.npy')
terrain = np.load('/glade/work/chennuo/data/global_land_data/terrain.npy')

f_date = sys.argv[1]
ninst = sys.argv[2]

data_dir = sys.argv[3]
output_sfc_path = data_dir + '/output_surface'+ninst+'.forecast-'+f_date+'.npy'
output_upper_path = data_dir + '/output_upper'+ninst+'.forecast-'+f_date+'.npy'

npy_sl = np.load(output_sfc_path)
npy_pl = np.load(output_upper_path)

pg_msl = npy_sl[0,:]
pg_u10 = npy_sl[1,:]
pg_v10 = npy_sl[2,:]
pg_t2m = npy_sl[3,:]

pg_z = npy_pl[0,:]
pg_q = npy_pl[1,:]
pg_t = npy_pl[2,:]
pg_u = npy_pl[3,:]
pg_v = npy_pl[4,:]

# Create a new NetCDF file
work_dir = './'
savename_nc = work_dir + 'pginput'+ninst+'.nc'
nctest = Dataset(savename_nc, "w", format="NETCDF4")

# Dataset Attributes
nctest.DLON = 0.25
nctest.DLAT = 0.25
nctest.DT = 21600
# https://github.com/NCAR/DART/blob/82f7e914400aae86486bd715d852a284cb656540/models/wrf/module_map_utils.f90#L71
nctest.MAP_PROJ_CHAR = "LATLON"
nctest.MAP_PROJ = 0
nctest.PERIODIC_X = 0
nctest.POLAR = 0

# Create dimensions
time = nctest.createDimension("time", None)  # unlimited dimension
south_north = nctest.createDimension("south_north", 721)
west_east = nctest.createDimension("west_east", 1440)
bottom_top = nctest.createDimension("bottom_top", 13)

# Create coordinate variables
Times = nctest.createVariable("time", "i", ("time",), fill_value=0)
cam_date = nctest.createVariable("date", "i", ("time",), fill_value=-999)
cam_datesec = nctest.createVariable("datesec", "f8", ("time",), fill_value=-999)

XLAT = nctest.createVariable("XLAT", "f4", ("south_north", "west_east"), fill_value=-900)
XLONG = nctest.createVariable("XLONG", "f4", ("south_north", "west_east"), fill_value=-900)
P_TOP = nctest.createVariable("P_TOP", "f4", ("time",), fill_value=0)

# Create 3D data variable
U = nctest.createVariable(
    "U",
    "f4",
    (
        "time",
        "bottom_top",
        "south_north",
        "west_east",
    ),
	fill_value=-999
)
V = nctest.createVariable(
    "V",
    "f4",
    (
        "time",
        "bottom_top",
        "south_north",
        "west_east",
    ),
	fill_value=-999
)
T = nctest.createVariable(
    "T",
    "f4",
    (
        "time",
        "bottom_top",
        "south_north",
        "west_east",
    ),
	fill_value=-999
)
Q = nctest.createVariable(
    "Q",
    "f4",
    (
        "time",
        "bottom_top",
        "south_north",
        "west_east",
    ),
	fill_value=-999
)

# Z = nctest.createVariable("Z","f4",("time","bottom_top","south_north", "west_east",)) # geopotential
PH = nctest.createVariable(
    "PH",
    "f4",
    (
        "time",
        "bottom_top",
        "south_north",
        "west_east",
    ),
	fill_value=-999
)  # geopotential
PHB = nctest.createVariable(
    "PHB",
    "f4",
    (
        "time",
        "bottom_top",
        "south_north",
        "west_east",
    ),
	fill_value=-999
)  # geopotential

# 2D variables
MUB = nctest.createVariable(
    "MUB",
    "f4",
    (
        "time",
        "south_north",
        "west_east",
    ),
	fill_value=-999
)
MU = nctest.createVariable(
    "MU",
    "f4",
    (
        "time",
        "south_north",
        "west_east",
    ),
	fill_value=-999
)
U10 = nctest.createVariable(
    "U10",
    "f4",
    (
        "time",
        "south_north",
        "west_east",
    ),
	fill_value=-999
)
V10 = nctest.createVariable(
    "V10",
    "f4",
    (
        "time",
        "south_north",
        "west_east",
    ),
	fill_value=-999
)
PSFC = nctest.createVariable(
    "PSFC",
    "f4",
    (
        "time",
        "south_north",
        "west_east",
    ),
	fill_value=-999
)
T2 = nctest.createVariable(
    "T2",
    "f4",
    (
        "time",
        "south_north",
        "west_east",
    ),
	fill_value=-999
)
XLAND = nctest.createVariable(
    "XLAND",
    "f4",
    (
        "south_north",
        "west_east",
    ),
)  # landmask
HGT = nctest.createVariable(
    "HGT",
    "f4",
    (
        "south_north",
        "west_east",
    ),
	fill_value=-999
)  # terrain height

# Assign data to coordinate variables

timevalue = datetime.strptime(f_date, '%Y-%m-%d-%H')
time_unit_out= "days since 1600-01-01 00:00:00"
Times.calendar = "gregorian"
Times[:] = date2num(timevalue,time_unit_out,calendar=Times.calendar)
Times.setncattr('unit',time_unit_out)
#Times.units = "days since 2024-01-20 06:00:00"
Times.long_name = "time"

tyear = timevalue.year
tday = timevalue.day
tmonth = timevalue.month
tdate = int(f'{tyear}{str(tmonth).zfill(2)}{str(tday).zfill(2)}')

cam_date[:] = np.array([tdate])
cam_date.long_name = "current date (YYYYMMDD)"

thour = timevalue.hour
thour*3600
cam_datesec[:] = np.array([thour*3600])
cam_datesec.long_name = "current seconds of current date"


P_TOP[:] = np.array([5000])
P_TOP.units = "Pa"
P_TOP.long_name = "PRESSURE TOP OF THE MODEL"
P_TOP.MemoryOrder = ""
P_TOP.FieldType = 104
P_TOP.stagger = ""

lat = np.linspace(-90, 90, 721)
long = np.linspace(0, 359.75, 1440)
long = (long + 180) % 360 - 180
XLONG[:], XLAT[:] = np.meshgrid(long, lat)
# stagger on T-grid is empty ""
XLONG.FieldType = 104
XLONG.MemoryOrder = "XY"
XLONG.long_name = "LONGITUDE, WEST IS NEGATIVE"
XLONG.stagger = ""
XLONG.units = "degrees_east"
XLONG.description = "LONGITUDE, WEST IS NEGATIVE"

XLAT.FieldType = 104
XLAT.MemoryOrder = "XY"
XLAT.long_name = "LATITUDE, SOUTH IS NEGATIVE"
XLAT.stagger = ""
XLAT.units = "degrees_north"
XLAT.description = "LATITUDE, SOUTH IS NEGATIVE"

# Assign data to data variable
U[0,:] = pg_u[:, ::-1, :]
U.FieldType = 104
U.MemoryOrder = "XYZ"
U.long_name = "x-wind component"
U.stagger = ""
U.units = "m s-1"
U.description = "x-wind component"

V[0,:] = pg_v[:, ::-1, :]
V.FieldType = 104
V.MemoryOrder = "XYZ"
V.long_name = "y-wind component"
V.stagger = ""
V.units = "m s-1"
V.description = "y-wind component"

T[0,:] = pg_t[:, ::-1, :]
T.FieldType = 104
T.MemoryOrder = "XYZ"
T.long_name = "Temperature"
T.stagger = ""
T.units = "K"
T.description = "air_temperature"

Q[0,:] = pg_q[:, ::-1, :]
Q.FieldType = 104
Q.MemoryOrder = "XYZ"
Q.long_name = "Specific_humidity"
Q.stagger = ""
Q.units = "kg kg-1"
Q.description = "specific_humidity"


ph_base = np.mean(pg_z[:, ::-1, :], axis=(-1, -2))
ph_base = np.broadcast_to(ph_base[None, :, None, None], (1, 13, 721, 1440))
PHB[0,:] = ph_base
PHB.FieldType = 104
PHB.MemoryOrder = "XYZ"
PHB.long_name = "base-state geopotential"
PHB.stagger = ""
PHB.units = "m2 s-2"
PHB.description = "base-state geopotential"

ph_pert = pg_z[:, ::-1, :] - ph_base
PH[0,:] = ph_pert
PH.FieldType = 104
PH.MemoryOrder = "XYZ"
PH.long_name = "perturbation geopotential"
PH.stagger = ""
PH.units = "m2 s-2"
PH.description = "perturbation geopotential"

PSFC[0,:] = pg_msl[::-1, :]
PSFC.FieldType = 104
PSFC.MemoryOrder = "XYZ"
PSFC.long_name = "Mean sea level pressure"
PSFC.stagger = ""
PSFC.units = "Pa"
PSFC.description = "air_pressure_at_mean_sea_level"

mu_base = np.nanmean(pg_msl[::-1, :], axis=(-1, -2))
mu_base = np.ones((721, 1440))*mu_base
MUB[0,:] = mu_base
MUB.FieldType = 104
MUB.MemoryOrder = "XY"
MUB.long_name = "base state dry air mass in column"
MUB.stagger = ""
MUB.units = "Pa"
MUB.description = "base state dry air mass in column"

MU[0,:] = pg_msl[::-1, :] - mu_base
MU.FieldType = 104
MU.MemoryOrder = "XY"
MU.long_name = "Perturbation dry column air mass"
MU.stagger = ""
MU.units = "Pa"
MU.description = "perturbation dry air mass in column"

U10[0,:] = pg_u10[::-1, :]
U10.FieldType = 104
U10.MemoryOrder = "XY"
U10.long_name = "10 meter U wind component"
U10.stagger = ""
U10.units = "m s-1"
U10.description = "U at 10 M"

V10[0,:] = pg_v10[::-1, :]
V10.FieldType = 104
V10.MemoryOrder = "XY"
V10.long_name = "10 meter V wind component"
V10.stagger = ""
V10.units = "m s-1"
V10.description = "V at 10 M"

T2[0,:] = pg_t2m[::-1, :].data
T2.FieldType = 104
T2.MemoryOrder = "XY"
T2.long_name = "2 meter temperature"
T2.stagger = ""
T2.units = "K"
T2.description = "TEMP at 2 M"

XLAND[:] = landmask
XLAND.FieldType = 104
XLAND.MemoryOrder = "XY"
XLAND.long_name = "land mask"
XLAND.stagger = ""
XLAND.units = ""
XLAND.description = "TLAND MASK (1 FOR LAND, 2 FOR WATER)"

HGT[:] = terrain
HGT.FieldType = 104
HGT.MemoryOrder = "XY"
HGT.long_name = "Terrain Height"
HGT.stagger = ""
HGT.units = ""
HGT.description = "Terrain Height"

# Close the NetCDF file
nctest.close()
