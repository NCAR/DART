import logging as log
import pandas as pd
import netCDF4 as nc4
import numpy as np
import re
import pydartdiags.obs_sequence.obs_sequence as obsq

def _iodaDtime2obsqDtime(epoch, iodaDtime):
    """Convert IODA datetime into obs_seq datetime

    This function will convert the IODA epoch style dateTime format to
    the obs_seq days, seconds, time format.

    Args:
        1. epoch - IODA style epoch (ISO 8601 string format) datetime
        2. iodaDtime - integer vector holding the offset in seconds from
                       the epoch reference datetime

    Returns:
        Tuple containing the obs_seq (seconds, days, time) values.

    """

    # Steps
    #    1. Convert the IODA datetime to numpy datetimes
    #    2. Generate the days, seconds and time values from the numpy datetimes
    #       Use python datetime objects to get the days since 1601-01-01
    #       via the Gregorian calendar.

    # Create the numpy datetimes
    epochDT = np.datetime64(epoch.rstrip('Z'))
    numLocs = len(iodaDtime)
    dTime = np.full(numLocs, epochDT, dtype='datetime64[s]')
    offset = iodaDtime.astype('timedelta64[s]')
    dTime = dTime + offset

    # Create the time vector
    time = np.char.replace(dTime.astype('str'), 'T', ' ')

    # Create the days
    gregDTref = np.datetime64('1601-01-01')
    days = (np.array(time, dtype='datetime64[D]') - gregDTref).astype(np.int32)

    # Create the seconds
    seconds = np.array(time, dtype='datetime64[s]').astype(np.int32) % 86400
    
    return (seconds, days, time)

def _ioda2iodaDF(iodaFile, mask_obsvalue=None, sort_by=None, derived=False):
    """Creates pandas dataframe from a IODA observation data file.

    Args:
        1. iodaFile - IODA observation data file
        mask_obsvalue (None, optional): Mask ObsValue _FillValue and NaN
            Set this to something like netCDF4.default_fillvals[dtype]
            where dtype can be 'f4' or 'f8' for floating point arrays
        sort_by (list of str, optional): Sort dataframe by columns
        derived: Set to true if using @DerivedObsValue

    Returns:
        pandas.DataFrame: IODA formatted dataframe

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

    return (df, epochdatetime)

def _buildObsSeqFromObsqDF(obsqDF, maxNumObs=None, nonQcNames=None, qcNames=None):
    """Construct a pyDARTdiags ObsSequence ojbect from an obs_seq format pandas dataframe

    Args:
        1. obsqDF - input obs_seq formatted dataframe
        2. maxNumObs - maximum number of obs for the file header
        3. nonQcNames - default setting for the ObsSequence object "non_qc_copie_names" attribute
        4. qcNames - default setting for the ObsSequence object "qc_copie_names" attribute

        The default for the maxNumObs (-1) tells this function to use the length of the
        input dataframe for the ObsSequence object's "max_num_obs" entry in its header attribute

    Return:
        This function returns an ObsSequence object based on the input dataframe

    """

    # set up defaults for function arguments
    if maxNumObs is None:
        maxNumObs = -1
    if nonQcNames is None:
        nonQcNames = [ 'observation' ]
    if qcNames is None:
        qcNames = [ 'Data_QC' ]
 
    # Steps:
    #   1. Construct an empty ObsSequence object
    #   2. assign the input dataframe to objects data frame (shallow copy)
    #   3. update the object attributes from the input dataframe
    #   4. fix the "copie" related ojbect attributes
    #   5. generate the header attribute from the ojbect attributes
    obsSeq = obsq.ObsSequence(file=None)
    obsSeq.df = obsqDF
    obsSeq.update_attributes_from_df()

    # An actual file header was never read, so the obsSeq methods have no
    # way of figuring out which columns are non-qc vs qc. They need hints
    # which are given by setting the non_qc_copie_names and qc_copie_names
    # attributes.
    obsSeq.non_qc_copie_names = nonQcNames
    obsSeq.qc_copie_names = qcNames
    obsSeq.n_non_qc = len(obsSeq.non_qc_copie_names)
    obsSeq.n_qc = len(obsSeq.qc_copie_names)

    if maxNumObs < 0:
        obsSeq.create_header(len(obsSeq.df))
    else:
        obsSeq.create_header(maxNumObs)

    return obsSeq

def _iodaDF2obsqDF(iodaDF, epochDT, iodaVarName, iodaVarType, vertCoordName, vertCoordUnits):
    """Reformat a pandas dataframe built with ioda2iodaDF into a dataframe suitable
       for buildObsSeqFromObsqDF.

     Args:
         1. iodaDF - pandas dataframe in the ioda layout
         2. epochDT - IODA dateTime epoch value in ISO 8601 string format

     Return:
         This function returns a pandas dataframe in the obs_seq layout. This dataframe
         is suitable as input to the buildObsSeqFromObsqDF function.

    """

    # The ioda dataframe layout has these features:
    #   1. Each separate variable is in a separate column
    #   2. The datetime variable (MetaData/dateTime) has
    #     a. A reference datetime in a variable attribute named "units"
    #     b. An offset in seconds from that reference in the variable values (int64)
    #   3. There are additional metadata variable (eg MetaData/staionIdentification)
    #      beyond lat, lon, datetime, and vertical coordinate

    # convert the datetime from the iodaDF into the obs_seq (seconds, days, time) format
    iodaDT = np.array(iodaDF['MetaData/dateTime'], dtype=np.int64)
    (obsqSec, obsqDays, obsqTime) = _iodaDtime2obsqDtime(epochDT, iodaDT)
    
    # For now, follow the radiosonde example which should be close or good for
    # other conventional obs types. Put in radiance obs types later.
    obsValName = 'ObsValue/' + iodaVarName
    obsErrorName = 'ObsError/' + iodaVarName
    obsQcName = 'PreQC/' + iodaVarName
    numLocs = len(iodaDT)
    emptyList = [ ]

    obsqDF = pd.DataFrame()
    obsqDF.insert(0, 'obs_num', np.arange(1, (numLocs + 1), 1, np.int32))
    obsqDF.insert(1, 'observation', np.array(iodaDF[obsValName], dtype=np.float32))
    obsqDF.insert(2, 'Data_QC', np.array(iodaDF[obsQcName], dtype=np.float32))
    obsqDF.insert(3, 'linked_list', np.full(numLocs, '1, 1, -1', dtype=object))

    obsqDF.insert(4, 'longitude', np.array(iodaDF['MetaData/longitude'], dtype=np.float32))
    obsqDF.insert(5, 'latitude', np.array(iodaDF['MetaData/latitude'], dtype=np.float32))
    obsqDF.insert(6, 'vertical', np.array(iodaDF[vertCoordName], dtype=np.float32))
    obsqDF.insert(7, 'vert_unit', np.full(numLocs, vertCoordUnits, dtype=object))

    obsqDF.insert(8, 'type', np.full(numLocs, iodaVarType, dtype=object))
    obsqDF.insert(9, 'metadata', np.full(numLocs, '', dtype=object))
    obsqDF.insert(10, 'external_FO', np.full(numLocs, '', dtype=object))

    obsqDF.insert(11, 'seconds', obsqSec)
    obsqDF.insert(12, 'days', obsqDays)
    obsqDF.insert(13, 'time', obsqTime)

    obsqDF.insert(14, 'obs_err_var', np.array(iodaDF[obsErrorName], dtype=np.float32))

    # Remove rows with missing obs data, or missing location, timestamp data
    obsqDF.dropna(subset=['observation', 'longitude', 'latitude', 'vertical', 'seconds', 'days'], inplace=True)

    return obsqDF
