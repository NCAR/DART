"""
Plot output from test_quad_irreg_interp.

Reads the six text files written by the Fortran program and produces
two side-by-side panels:
  Left  – data grid: scatter of the 9x5 corner points, coloured by data value
  Right – sample grid: pcolormesh of the 210x150 interpolated values

Run from the work/ directory where the output files were created:
    python ../plot_irreg_interp.py
or supply a path:
    python plot_irreg_interp.py --dir /path/to/work
"""

import argparse
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# ---------------------------------------------------------------------------
# Grid dimensions – must match the Fortran program constants
# ---------------------------------------------------------------------------
NX, NY   = 9, 5          # data grid
NRX, NRY = 210, 150      # sampling grid


def load_1d(path, n):
    """Read n values, one per line, return 1-D array."""
    data = np.loadtxt(path)
    if data.size != n:
        raise ValueError(f"{path}: expected {n} values, got {data.size}")
    return data


def load_2d(path, ni, nj):
    """
    Read ni*nj values written by the Fortran loops
        do i=1,ni; do j=1,nj; write(*,*) arr(i,j); enddo; enddo
    Returns array shaped (ni, nj).
    """
    data = np.loadtxt(path)
    if data.size != ni * nj:
        raise ValueError(f"{path}: expected {ni*nj} values, got {data.size}")
    return data.reshape(ni, nj)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--dir', default='.', help='Directory containing the output text files (default: .)')
    parser.add_argument('--out', default=None, help='Save figure to this file instead of displaying it')
    args = parser.parse_args()

    d = args.dir

    def path(fname):
        p = os.path.join(d, fname)
        if not os.path.exists(p):
            sys.exit(f"File not found: {p}\n"
                     f"Run test_quad_irreg_interp first, then re-run this script from the work/ directory.")
        return p

    # ------------------------------------------------------------------
    # Load data grid  (9 x 5)
    # ------------------------------------------------------------------
    data_lons = load_2d(path('data_lons_2d_irreg_test.txt'), NX, NY)
    data_lats = load_2d(path('data_lats_2d_irreg_test.txt'), NX, NY)
    data_vals = load_2d(path('data_data_2d_irreg_test.txt'), NX, NY)

    # ------------------------------------------------------------------
    # Load sample grid  (210 x 150)
    # ------------------------------------------------------------------
    sample_lons = load_1d(path('sample_lons_1d_irreg_test.txt'), NRX)
    sample_lats = load_1d(path('sample_lats_1d_irreg_test.txt'), NRY)
    sample_vals = load_2d(path('sample_data_2d_irreg_test.txt'), NRX, NRY)

    # sample_vals is (nrx, nry) = (lon-index, lat-index).
    # pcolormesh expects (ny, nx) so transpose before plotting.
    LON, LAT = np.meshgrid(sample_lons, sample_lats)   # both (nry, nrx)
    sample_plot = sample_vals.T            # shape (nry, nrx)

    # ------------------------------------------------------------------
    # Colour scale: shared between both panels, excluding MISSING values
    # DART MISSING_R8 is a large negative number; mask values below -1e30
    # ------------------------------------------------------------------
    MISSING_THRESHOLD = -1.0e29

    data_masked   = np.where(data_vals   < MISSING_THRESHOLD, np.nan, data_vals)
    sample_masked = np.where(sample_plot < MISSING_THRESHOLD, np.nan, sample_plot)

    valid = np.concatenate([data_masked[np.isfinite(data_masked)],
                            sample_masked[np.isfinite(sample_masked)]])
    vmin, vmax = valid.min(), valid.max()

    cmap = plt.get_cmap('RdYlBu_r')
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    # ------------------------------------------------------------------
    # Figure
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(14, 6),
                             constrained_layout=True)

    # ---------- Left panel: data grid (scatter) ----------
    ax = axes[0]
    lons_flat = data_lons.ravel()
    lats_flat = data_lats.ravel()
    vals_flat = data_masked.ravel()

    sc = ax.scatter(lons_flat, lats_flat,
                    c=vals_flat, cmap=cmap, norm=norm,
                    s=80, edgecolors='k', linewidths=0.5, zorder=3)

    # Draw quad edges (connect adjacent grid points)
    for i in range(NX):
        ax.plot(data_lons[i, :], data_lats[i, :], 'k-', lw=0.4, alpha=0.4)
    for j in range(NY):
        ax.plot(data_lons[:, j], data_lats[:, j], 'k-', lw=0.4, alpha=0.4)

    ax.set_xlabel('Longitude (°)')
    ax.set_ylabel('Latitude (°)')
    ax.set_title(f'Data grid ({NX}×{NY} corners)')
    fig.colorbar(sc, ax=ax, label='Data value', shrink=0.85)

    # ---------- Right panel: interpolated values at sample points (scatter) ----------
    ax = axes[1]
    sc2 = ax.scatter(LON.ravel(), LAT.ravel(),
                     c=sample_masked.ravel(),
                     cmap=cmap, norm=norm,
                     s=2, linewidths=0)

    # Overlay the data grid corners for reference
    ax.scatter(data_lons.ravel(), data_lats.ravel(),
               c='k', s=20, zorder=3, label='Data corners')
    ax.legend(loc='upper right', fontsize=8)

    ax.set_xlabel('Longitude (°)')
    ax.set_ylabel('Latitude (°)')
    ax.set_title(f'Interpolated values at sample points ({NRX}×{NRY})')
    fig.colorbar(sc2, ax=ax, label='Interpolated value', shrink=0.85)

    fig.suptitle('test_quad_irreg_interp output', fontsize=13)

    if args.out:
        fig.savefig(args.out, dpi=150)
        print(f'Saved to {args.out}')
    else:
        plt.show()


if __name__ == '__main__':
    main()
