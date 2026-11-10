"""
preprocess_ocean_obs.py  --  Filter an obs_seq file by bathymetric depth and
                             optionally replace observation error with a depth-
                             dependent model.

When --obs-type is given, only that type is depth-filtered and error-updated;
all other types are kept unchanged. Without --obs-type, all obs are processed.

Usage examples
--------------
# All obs, depth filter + depth-dependent error (defaults):
python preprocess_ocean_obs.py obs_seq.all obs_seq.all_trim \
    --roms-file roms_restart.nc

# All obs, depth filter only (keep original error variance):
python preprocess_ocean_obs.py obs_seq.all obs_seq.all_trim \
    --roms-file roms_restart.nc --no-depth-error

# One obs type, depth filter + error update, custom depth cut-off
# (other types in the file are kept unchanged):
python preprocess_ocean_obs.py obs_seq.all obs_seq.deep \
    --roms-file roms_restart.nc \
    --obs-type SATELLITE_SSH \
    --min-depth 500

# All obs, error update only (no depth filtering):
python preprocess_ocean_obs.py obs_seq.all obs_seq.err_updated \
    --roms-file roms_restart.nc --no-depth-filter

"""

# Module imports
from scipy.spatial import cKDTree

import argparse
import sys
import pandas                                as pd
import numpy                                 as np
import xarray                                as xr
import pydartdiags.obs_sequence.obs_sequence as obsq


# CLI
def parse_args():
    p = argparse.ArgumentParser(
        description="Trim an obs_seq file by ROMS bathymetric depth.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required positional arguments
    p.add_argument("obs_in",  help="Input obs_seq file")
    p.add_argument("obs_out", help="Output obs_seq file")
    p.add_argument("--roms-file", required=True,
                   help="ROMS NetCDF file with h, mask_rho, lon_rho, lat_rho")
    p.add_argument("--obs-type", default=None,
                   help="Process only this obs type; all others pass through unchanged")

    # Depth filter
    filt = p.add_argument_group("depth filter")
    filt.add_argument("--no-depth-filter", dest="use_depth_filter",
                      action="store_false", default=True,
                      help="Keep obs at all depths")
    filt.add_argument("--min-depth", type=float, default=200.0,
                      help="Remove obs shallower than this [m]")

    # Depth-dependent error
    err = p.add_argument_group("depth-dependent error")
    err.add_argument("--no-depth-error", dest="use_depth_error",
                     action="store_false", default=True,
                     help="Keep original obs error variance")
    err.add_argument("--sigma-min", type=float, default=0.04,
                     help="Deep-ocean obs error std-dev [obs units]")
    err.add_argument("--sigma-max", type=float, default=0.08,
                     help="Shallow-water obs error std-dev [obs units]")
    err.add_argument("--h0", type=float, default=500.0,
                     help="Depth scale for error transition [m]")

    args = p.parse_args()

    if not args.use_depth_filter and not args.use_depth_error:
        p.error("--no-depth-filter and --no-depth-error together leave nothing to do.")

    return args


# Depth lookup: vectorised via KDTree
def lookup_depths(lon_obs, lat_obs, glon, glat, bath, mask):
    """
    Return the bathymetric depth at the nearest wet grid point for every
    observation. This is an approximate depth estimate because observations 
    are not necessarily collocated with the model grid. 
    """
    wet    = mask > 0.5   # 0: land, 1: wet
    tree   = cKDTree(np.column_stack([glon[wet], glat[wet]]))
    _, idx = tree.query(np.column_stack([lon_obs, lat_obs]))
    return bath[wet][idx]


# Depth-dependent error model
def depth_dependent_error_std(depth, sigma_min, sigma_max, h0):
    """Exponential depth-to-error mapping; returns std-dev (not variance)."""
    return sigma_min + (sigma_max - sigma_min) * np.exp(-depth / h0)


# Main
def main():
    args = parse_args()

    # Read the input obs_seq and ROMS grid
    print(f"\nReading obs_seq:  {args.obs_in}")
    obs_seq = obsq.ObsSequence(args.obs_in)

    print(f"Reading ROMS grid: {args.roms_file}")
    with xr.open_dataset(args.roms_file) as roms:    
    	bath = roms["h"].values
    	mask = roms["mask_rho"].values
    	glon = np.mod(roms["lon_rho"].values, 360.0)  # DART-style
    	glat = roms["lat_rho"].values

    # Which observation type(s) to process?
    df = obs_seq.df.copy()

    # First, find out what obs are available in file: 
    available = sorted(df["type"].unique())
    print(f"\nObs types in file: {available}")

    if args.obs_type:
        # User selected obs type
        if args.obs_type not in available:
            print(f"ERROR: '{args.obs_type}' not found in file.\n"
                  f"Available: {available}")
            sys.exit(1)
        obs_df   = df[df["type"] == args.obs_type].copy()
        obs_rest = df[df["type"] != args.obs_type].copy()
        label    = args.obs_type
    else:
        obs_df   = df.copy()
        obs_rest = None
        label    = "all"

    Nobs_in = len(obs_df)
    print(f"\nInput obs ({label}): {Nobs_in}")

    if Nobs_in == 0:   # empty?
        print("No observations to process — exiting.")
        sys.exit(0)

    # Figure out the options for depth filtering and 
    # the error update. 
    if args.use_depth_filter or args.use_depth_error:
        # Depth lookup needed for filter and/or error update

        print("Looking up bathymetric depths ...")
        lon_obs = obs_df["longitude"].values     
        lat_obs = obs_df["latitude"].values
        depth   = lookup_depths(lon_obs, lat_obs, glon, glat, bath, mask)
    else:
        depth = None

    if args.use_depth_filter and depth is not None:
        keep       = np.isfinite(depth) & (depth > args.min_depth)
        obs_kept   = obs_df.loc[keep].copy()
        depth_kept = depth[keep]
    else:
        obs_kept   = obs_df.copy()
        depth_kept = depth  

    if args.use_depth_error and depth_kept is not None:
        sigma = depth_dependent_error_std(
            depth_kept, args.sigma_min, args.sigma_max, args.h0)
        obs_kept["obs_err_var"] = sigma ** 2

    # Merge back and re-number:
    if obs_rest is not None:
        obs_final = pd.concat([obs_rest, obs_kept], ignore_index=True)
        obs_final = obs_final.sort_values(["days", "seconds", "obs_num"]). \
                    reset_index(drop=True)
    else:
        obs_final = obs_kept.reset_index(drop=True)

    obs_final["obs_num"] = np.arange(1, len(obs_final) + 1)

    # Write the output obs_seq
    obs_seq.df = obs_final

    obs_seq.update_attributes_from_df()
    obs_seq._update_linked_list(obs_final)  # private, no public API for this?
    obs_seq.write_obs_seq(args.obs_out)

    # Print summary
    print(f"\nSummary ({label})")
    
    if args.use_depth_filter:
        print(f"  Kept:    {len(obs_kept)}  (depth > {args.min_depth:.0f} m)")
        print(f"  Removed: {Nobs_in - len(obs_kept)}")
        if len(obs_kept):
            print(f"  Depth range kept: {depth_kept.min():.1f} – {depth_kept.max():.1f} m")
    else:
        print(f"  Depth filter: off | all {Nobs_in} obs kept")

    if args.use_depth_error:
        print(f"  Error model: sigma = {args.sigma_min} + "
              f"({args.sigma_max} - {args.sigma_min}) * exp(-h / {args.h0})")
    else:
        print(f"  Error model: Original obs error variance kept")

    if obs_rest is not None:
        print(f"\n  Other obs types kept unchanged: {len(obs_rest)}")
    
    print(f"\n  Wrote: {args.obs_out}")

if __name__ == "__main__":
    main()
