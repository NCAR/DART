import numpy as np
import lorenz as l96
import netCDF4 as nc
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("reference", type=str, help="Reference File to Begin Creating Ensemble")
parser.add_argument("destination", type=str, help="Filepath to save ensemble")
parser.add_argument("-n", "--n_members", type=int, default=10)
parser.add_argument("-i", "--init_style", type=int, help="0 = Use Default to Initialise, 1 = Use Reference to Initialise", default=0)
parser.add_argument("-p", "--perturb_style", type=int, help="0 = Perturb Starting Location, 1 = Perturb All Locations", default=1)
args = parser.parse_args()

if args.init_style == 0:
    default_cycles = 100
    init_x = np.zeros(40)
    init_x[0] = 1
    reference = l96.vanilla_RK4(x=init_x, F=8, size=40, cycles=default_cycles, step=0.05)[1][-1]
elif args.init_style == 1:
    my_file = nc.Dataset(args.reference, "r")
    reference = my_file["state"][:]
    my_file.close()
else:
    print("Unrecognised Initialisation Style")
    raise SystemExit

RNG = np.random.default_rng()
ref_file = nc.Dataset(args.reference, "r")
for i_mem in range(args.n_members):
    new_filepath = os.path.join(args.destination, "l96_mem{:05}.nc".format(i_mem + 1))
    mem_file = nc.Dataset(new_filepath, "w")
    mem_file.setncatts(ref_file.__dict__)
    for dim_name, dimension in ref_file.dimensions.items():
        mem_file.createDimension(dim_name, (len(dimension) if not dimension.isunlimited() else None))
    for var_name, variable in ref_file.variables.items():
        mem_file.createVariable(var_name, variable.datatype, variable.dimensions)
        mem_file[var_name].setncatts(ref_file[var_name].__dict__)
        mem_file[var_name][:] = ref_file[var_name][:]

    if args.perturb_style == 0:
        mem_file["state"][:] = reference
        mem_file["state"][0] += RNG.normal(0, 1)
    elif args.perturb_style == 1:
        state_length = len(reference)
        mem_file["state"][:] = reference + np.array([RNG.normal(0,1) for i_perturb in range(state_length)])
    else:
        print("Unrecognised Perturb Style")
        raise SystemExit
    mem_file.close()
ref_file.close()

