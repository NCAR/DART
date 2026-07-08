import argparse
import datetime
import os
import pathlib
import sys

# ###################################################################
def get_ens_dotfile_end_datetime(run_dir):
    """Use the the .model_end_time files to get the current ensemble time."""
    run_dir = pathlib.Path(run_dir)
    mem_dirs = sorted(run_dir.glob("member_*"))

    def read_dot_file(file):
        with open(file) as f:
            content = f.readline()
        return datetime.datetime.strptime(content, '%Y-%m-%d %H:%M:%S')

    end_times = [read_dot_file(mm / '.model_end_time') for mm in mem_dirs]
    if not all([end_times[0] == ee for ee in end_times]):
        raise ValueError("Not all ensemble members at the same time (HYDRO_RST files).")

    return end_times[0]

# ####################################################################
def get_ensemble_time(
    run_dir: str,
    with_previous: int=0
):
    current_time = get_ens_dotfile_end_datetime(run_dir)

    if with_previous != 0:
        previous_time = current_time - datetime.timedelta(hours=with_previous)
        return previous_time, current_time
    else:
        return current_time


if __name__ == "main":

    parser = argparse.ArgumentParser(
        description='Get the time of a simulation ensemble'
    )

    parser.add_argument(
        '--run_dir',
        required=False,
        metavar='/abs/path/to/run/directory',
        help='Path to the experiment directory. (Default is director where this script or a ' +
             'symlink (not resolved) to it lives).',
        default= os.path.dirname(os.path.abspath(__file__))
    )

    parser.add_argument(
        '--with_previous',
        required=False,
        metavar='delta_time_hours',
        help='returns a tuple"previous|current"',
        default='0'
    )

    args = parser.parse_args()
    run_dir = pathlib.PosixPath(args.run_dir)
    with_previous = int(args.with_previous)

    if with_previous == 1:
        previous_time, current_time = get_ensemble_time(args.run_dir, args.with_previous)
        print(
            previous_time.strftime('%Y-%m-%d_%H:%M') +
            '|' +
            current_time.strftime('%Y-%m-%d_%H:%M')
        )
    else:
        current_time = get_ensemble_time(args.run_dir, args.with_previous)
        print(current_time.strftime('%Y-%m-%d_%H:%M'))
    
    sys.exit(0)
