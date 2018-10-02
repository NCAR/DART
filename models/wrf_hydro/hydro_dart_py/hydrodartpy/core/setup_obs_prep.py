import pickle
import wrfhydropy

from .setup_usgs_daily import setup_usgs_daily


def setup_obs_prep(config):

    obs_prep_config = config['observation_preparation']

    ## THIS SCTIPT WILL EVENTUALLY HANDLE MORE OBS PREP

    if 'USGS_daily' in obs_prep_config.keys():
        if obs_prep_config['USGS_daily']['prepare']:
            success = setup_usgs_daily(config)


if __name__ == "__main__":
    pass
