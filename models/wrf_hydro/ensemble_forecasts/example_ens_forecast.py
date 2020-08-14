import datetime
import pathlib
import pickle
import wrfhydropy

# Number of cores available (18 is using interactive, shared queue example)
# and using one core per model advance.
n_cores = 18 

da_exp_dir = pathlib.Path(
    '/glade/scratch/jamesmcc/wrfhydro_dart/flo_cut/runs/bucket_multiphys_six_routelink_params')

# Ensemble
ens_pkl = da_exp_dir / 'WrfHydroEnsembleSim.pkl'
ens = pickle.load(ens_pkl.open('rb'))
ens.restore_members(ens_pkl.parent)
ens.ncores = 1

# "refresh" / reuse / repurpose the ensemble
del ens._compose_dir
# There are some relative paths that need edited.
ens.jobs = []
for mm in ens.members:
    mm.jobs = []
    mm.base_hrldas_namelist['noahlsm_offline']['indir'] = './FORCING_AnA_channel-only'


# Ensemble forecast
first_time = datetime.datetime(2018, 9, 16, 0)
init_times = [first_time + datetime.timedelta(hours=h) for h in range(24)]

mem_dirs = sorted(da_exp_dir.glob('member_*'))
restart_dirs = [mem_dirs for ii in range(len(init_times))]

forc_dirs = sorted(da_exp_dir.glob('member_*/FORCING_perturbed'))
forcing_dirs = [forc_dirs for ii in range(len(init_times))]

fcst = wrfhydropy.CycleSimulation(
    init_times=init_times,
    restart_dirs=restart_dirs,
    forcing_dirs=forcing_dirs,
    ncores=n_cores
)

# Add ens
fcst.add(ens)

# Add Job
model_start_time = init_times[0]
model_end_time = init_times[0] + datetime.timedelta(hours=18)
#exe_cmd = 'mpirun -np {0} ./wrf_hydro.exe' # for a scheduler
exe_cmd = 'mpirun -np 1 ./wrf_hydro.exe'
job_name = 'ensfcst'
job = wrfhydropy.Job(
    exe_cmd=exe_cmd,
    job_id=job_name,
    restart=True,
    model_start_time= model_start_time,
    model_end_time=model_end_time
)
fcst.add(job)

# Write to disk
fcst_dir = pathlib.Path('/glade/scratch/jamesmcc/ens_fcst_example/')
if not fcst_dir.exists():
    fcst_dir.mkdir()
    os.chdir(str(fcst_dir))
    fcst.compose()

#fcst.pickle('WrfHydroCycle.pkl')

fcst.run(n_concurrent=n_cores)

# Collect!

# This will eventually be integrated with the forecast object.
from wrfhydropy.core.collection import open_whp_dataset
hydro_rst_files = sorted(fcst_dir.glob('*/*/HYDRO_RST.*DOMAIN1'))
len(hydro_rst_files) == 24*19*80

# Remove the time zero restart file, it is a full restart file and the variables
# are more extenzive. Could alternatively add to the hydro_rst_vars_drop.
first_file = 'HYDRO_RST.2018-08-01_00:00_DOMAIN1'
hydro_rst_files = [ff for ff in hydro_rst_files if first_file not in ff.name]
len(hydro_rst_files) == 24*18*80

# Just leave qlink1
hydro_rst_vars_drop = [
    'hlink',
    'qlakei',
    'qlakeo',
    'qlink2',
    'resht',
    'z_gwsubbas'
]
    
hydro_rst_ds = open_whp_dataset(
    paths=hydro_rst_files,
    drop_variables=hydro_rst_vars_drop,
    n_cores=n_cores
)

hydro_rst_ds

# Out[23]: 
# <xarray.Dataset>
# Dimensions:         (lead_time: 18, links: 1642, member: 80, reference_time: 24)
# Coordinates:
#   * reference_time  (reference_time) datetime64[ns] 2018-09-16 ... 2018-09-16T23:00:00
#   * member          (member) int64 0 1 10 11 12 13 14 15 ... 75 76 77 78 79 8 9
#   * lead_time       (lead_time) timedelta64[ns] 01:00:00 02:00:00 ... 18:00:00
# Dimensions without coordinates: links
# Data variables:
#     qlink1          (lead_time, member, reference_time, links) float32 0.28081888 ... 127.21303
#     valid_time      (lead_time, reference_time) datetime64[ns] 2018-09-16T01:00:00 ... 2018-09-17T17:00:00

hydro_rst_ds.to_netcdf(fcst_dir / 'HYDRO_RST_all.nc')

# (dart) jamesmcc@cheyenne6[1119]:/glade/scratch/jamesmcc/ens_fcst_example> ncdump -h HYDRO_RST_all.nc 
# netcdf HYDRO_RST_all {
# dimensions:
# 	reference_time = 24 ;
# 	member = 80 ;
# 	lead_time = 18 ;
# 	links = 1642 ;
# variables:
# 	int64 reference_time(reference_time) ;
# 		reference_time:units = "hours since 2018-09-16 00:00:00" ;
# 		reference_time:calendar = "proleptic_gregorian" ;
# 	int64 member(member) ;
# 	int64 lead_time(lead_time) ;
# 		lead_time:units = "hours" ;
# 	float qlink1(lead_time, member, reference_time, links) ;
# 	int64 valid_time(lead_time, reference_time) ;
# 		valid_time:units = "hours since 2018-09-16 01:00:00" ;
# 		valid_time:calendar = "proleptic_gregorian" ;
