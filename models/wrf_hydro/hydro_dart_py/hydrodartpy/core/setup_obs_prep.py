import pickle
import warnings
import wrfhydropy
from .setup_usgs_daily import setup_usgs_daily


def setup_obs_prep(config, config_file=None):

    obs_prep_config = config['observation_preparation']

    ## THIS SCTIPT WILL EVENTUALLY HANDLE MORE OBS PREP

    if 'USGS_daily' in obs_prep_config.keys():
        if obs_prep_config['USGS_daily']['prepare']:
            success = setup_usgs_daily(config, config_file=config_file)

    if 'success' not in locals():
        warnings.warn("Observation prepartion was not requested.")
        success = 1
    elif success != 0:
        warnings.warn("Observation prepartion was not succesful.")

    return success


if __name__ == "__main__":
    pass
