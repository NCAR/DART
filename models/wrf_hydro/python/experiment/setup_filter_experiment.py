import argparse
from hydrodartpy import setup_experiment
import sys

if __name__ == "__main__":

    # Arguments to this script.
    # python setup_experiment.py --help   
    the_desc = "Setup a WRF-Hydro-DART experiment."
    the_epilog = """
    Additional Examples: 
    regular: python setup_experiment.py experiment_config_files/test_exp.yaml
    debug:  ipython --pdb -c "%run setup_experiment.py experiment_config_files/test_exp.yaml"
    """
    parser = argparse.ArgumentParser(
        description=the_desc,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=the_epilog
    )

    # Single positional argument is config_file
    parser.add_argument(
        'config_file',
        metavar='config_file.yaml',
        type=str,
        nargs=1,
        help='The YAML experiment configuration file of arbitrary name.'
    )
    args = parser.parse_args()

    return_code = setup_experiment(config_file=args.config_file[0])
    
    sys.exit(return_code)
