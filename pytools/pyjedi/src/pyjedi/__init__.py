import logging as log
import pandas as pd
import netCDF4 as nc4
from datetime import datetime, timedelta
import numpy as np
import re

def ioda2df(iodaFile, mask_obsvalue=None, sort_by=None, derived=False):
    """Creates pandas dataframe from netCDF.

    Note:

        ObsBias = obsvalue - corvalue,  OR
        ObsBias = initial_obsvalue - corvalue
        corvalue = @ObsValue - @ObsBias = bgvalue + fg_depar, if @ObsBias is valid  # noqa
        corvalue = @ObsValue, otherwise

    Args:
        mask_obsvalue (None, optional): Mask ObsValue _FillValue and NaN
            Set this to something like netCDF4.default_fillvals[dtype]
            where dtype can be 'f4' or 'f8' for floating point arrays
        sort_by (list of str, optional): Sort dataframe by columns
        derived: Set to true if using @DerivedObsValue

    Returns:
        pandas.DataFrame: Dataframe of 1D arrays in the NetCDF with
        exceptions for DateTime and Station ID which are 2D arrays.

    """
    df = pd.DataFrame()
    with nc4.Dataset(iodaFile) as ds:
        groups = [ds.groups[k] for k in ds.groups.keys()]
        dsets = []
        for i, j in enumerate(groups):
            dsets += [ds[j.name + '/' + s] for s in groups[i].variables]

        # With epoch time representation need offset (window start)
        assert('dateTime' in ds.groups['MetaData'].variables)
        epochdatetime = ds.groups['MetaData'].variables['dateTime'].units.split()[-1]
        dt_epoch = datetime.strptime(epochdatetime, '%Y-%m-%dT%H:%M:%SZ')

        loc = 0
        for var in dsets:
            fullVarName = var.group().name + "/" + var.name
            if (var.dimensions[0] != 'Location'):
                print('WARNING: Can only process variables with Location as the first dimension')
                print('WARNING:     skipping variable: ', fullVarName)
                continue

            if var.group().name != 'VarMetaData':
                if var.ndim == 2:
                    if var.name in ("dateTime", "stationIdentifier"):
                        outarr = ["".join(i) for i in var[:].astype(str)]
                        df.insert(loc, fullVarname, outarr)
                        loc += 1
                    # multi-channel BT
                    if var.name in ("brightnessTemperature", "emissivity"):
                        if 'Channel' not in var.dimensions:
                            log.warning("Expected channel dim not found")
                            continue

                        # 2D BT are usually of (Location, Channel) shape
                        # but the order may not be guaranteed so
                        # force nchan as first dim before looping over
                        # individual channels
                        channelIndex = var.dimensions.index('Channel')
                        outarr = np.moveaxis(var[:], channelIndex, 0)
                        for ci in range(0, channelIndex):
                            df.insert(loc, fullVarName, outarr[ci, :])
                            loc += 1

                if var.ndim == 1:
                    df.insert(loc, fullVarName, var[:])
                    loc += 1

    # - v1. mask ObsValue >= mask_obsvalue e.g., nc4.default_fillvals['f4']
    #   v2.               <= -3.3687953E38, numpy.finfo('float32').min * 0.99
    if mask_obsvalue is not None:
        if derived:
            r = re.compile('DerivedObsValue/.*')  # v2
        else:
            r = re.compile('ObsValue/.*')  # v2

        # note: filter returns an iterator in Python3 so we use
        # list(filter(...)) as a fail-safe approach
        obscols = list(filter(r.match, df.columns))
        # alternatively, use list comprehension
        # obscols = [col for col in df.columns if '@ObsValue' in col]

        mask = df[obscols] <= mask_obsvalue
        df[obscols] = df[obscols].where(~mask, np.nan)

    # - sort dataframe by sequence number
    if sort_by is not None:
        df.sort_values(by=sort_by, kind='mergesort', inplace=True)

    return df
