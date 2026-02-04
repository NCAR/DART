import logging as log
import pandas as pd
import netCDF4 as nc4
import numpy as np
import re
import yaml
import pydartdiags.obs_sequence.obs_sequence as obsq

log.basicConfig(level=log.INFO)
log = log.getLogger(__name__)
missingVal = -888888.0

def _expandChanNums(chanNumsStr):
    """Expand channel numbers YAML spec into a list of channel numbers

    Args: chanNumsStr (str): Comma-separated string of channel numbers with optional ranges

    The optional ranges is of the form (<start>-<end>) where both <start> and <end> are
    inclusive. Eg, 12-15 expands to [12, 13, 14, 15].

    """
    chanNums = []
    for part in chanNumsStr.split(','):
        part = part.strip()
        if '-' in part:
            start, end = part.split('-')
            chanNums.extend(range(int(start), int(end) + 1))
        else:
            chanNums.append(int(part))
    return chanNums

def _parseYamlConfig(configFile):
    """Parse the YAML configuration file

    Args:
        configFile (str): Path to the YAML configuration file

    Returns:
        Tuple containing the processed "observation variables" configuration, and
        the processed "observation category" configuration

    The configuration file should contain the following sections:
      - observation variables:
        Required
        - List of IODA variable names and types (and optional channels) to convert
      - observation category:
        Required
        - name
          Name of the category. For now, accept: "radiance" or "conventional"
          Some possibilities for the future could be "radar", "gpsro", etc.
        The remaining contents of this section are dependent on the name of the category
        - vertical coordinate:
          Required for radiance and conventional obs types
          name and units of the vertical coordinate variable for conventional
          units and data value of the vertical coordinate variable for radiance
        - channel numbers:
          Required for radiance obs types
          list of numbers which can include ranges of numbers (eg, 12-15 for 12, 13, 14, 15)
        - metadata:
          Required for radiance obst types
          sensor key and rttov sensor db (path to the db file)

    """

    with open(configFile, 'r') as file:
        config = yaml.safe_load(file)
    iodaVarsConfig = config['ioda to obsq converter']['observation variables']
    obsCategoryConfig = config['ioda to obsq converter']['observation category']

    # check for proper values of the obs category name
    if (obsCategoryConfig['name'] not in ['radiance', 'conventional']):
        raise ValueError("Invalid observation category name: " + obsCategoryConfig['name'])
    obsCategory = obsCategoryConfig['name']

    # check for proper obs category configuration
    if (obsCategory == 'conventional'):
        # Conventional obs type.
        if ('vertical coordinate' not in obsCategoryConfig):
            raise ValueError("Must specify 'vertical coordinate' for conventional observations.")
        else:
            # need name and units for vertical coordinate
            if ('name' not in obsCategoryConfig['vertical coordinate']):
                raise ValueError("Must specify 'name' for vertical coordinate in conventional observations.")
            if ('units' not in obsCategoryConfig['vertical coordinate']):
                raise ValueError("Must specify 'units' for vertical coordinate in conventional observations.")
    elif (obsCategory == 'radiance'):
        # Radiance obs type.
        if ('vertical coordinate' not in obsCategoryConfig):
            raise ValueError("Must specify 'vertical coordinate' for radiance observations.")
        else:
            # need units and data value for vertical coordinate
            if ('units' not in obsCategoryConfig['vertical coordinate']):
                raise ValueError("Must specify 'units' for vertical coordinate in radiance observations.")
            if ('data value' not in obsCategoryConfig['vertical coordinate']):
                raise ValueError("Must specify 'data value' for vertical coordinate in radiance observations.")
        if ('channel numbers' not in obsCategoryConfig):
            raise ValueError("Must specify 'channel numbers' for radiance observations.")
        if ('metadata' not in obsCategoryConfig):
            raise ValueError("Must specify 'metadata' for radiance observations.")
        else:
            # need sensor key and rttov sensor db for metadata
            if ('sensor key' not in obsCategoryConfig['metadata']):
                raise ValueError("Must specify 'sensor key' for metadata in radiance observations.")
            if ('rttov sensor db' not in obsCategoryConfig['metadata']):
                raise ValueError("Must specify 'rttov sensor db' for metadata in radiance observations.")
            if ('sat az variable' not in obsCategoryConfig['metadata']):
                raise ValueError("Must specify 'sat az variable' for metadata in radiance observations.")
            if ('sat ze variable' not in obsCategoryConfig['metadata']):
                raise ValueError("Must specify 'sat ze variable' for metadata in radiance observations.")

        # Expand the channel numbers into the variable names by creating variable
        # names with the given name and each channel number appended to the end of the given name.
        # This is how the separate channels are stored in the ioda dataframe (see _ioda2iodaDF below).
        channelNumbers = _expandChanNums(obsCategoryConfig['channel numbers'])
        newIodaVarsConfig = []
        for iodaVarConfig in iodaVarsConfig:
            varName = iodaVarConfig['name']
            varType = iodaVarConfig['type']
            for chanNum in channelNumbers:
                newConfig = {}
                newConfig['name'] = f"{varName}_{chanNum}"
                newConfig['type'] = varType
                newIodaVarsConfig.append(newConfig)
        iodaVarsConfig = newIodaVarsConfig

        # Replace the channel numbers string in the category config with the expanded channel numbers
        obsCategoryConfig['channel numbers'] = channelNumbers

        # Build the sensor db dictionary, and extract the sensor key entry
        # Do this here to be able to look up the sensor entry in one shot.
        sensorDB = _buildSensorDict(obsCategoryConfig)
        sensorKey = obsCategoryConfig['metadata']['sensor key']
        sensorDbEntry = sensorDB.get(sensorKey, None)
        if (sensorDbEntry == None):
            raise ValueError("No entry in the sensor DB for sensor key: " + sensorKey)
        obsCategoryConfig['metadata']['sensor db entry'] = sensorDbEntry

    else:
        # Currently, can only use conventional or radiance
        raise ValueError("Invalid observation category name: " + obsCategory, \
                         ", must use one of 'radiance' or 'conventional'")

    return iodaVarsConfig, obsCategoryConfig

def _buildSensorDict(obsCategoryConfig):
    """Build the sensor database dictionary from the observation category config.

    Args:
        obsCategoryConfig (dict): Observation category configuration

    Returns:
        dict: Sensor database dictionary
    """

    sensorDB = {}
    sensorDF = pd.read_csv(obsCategoryConfig['metadata']['rttov sensor db'], header=None)
    sensorDF.columns = ['key', 'platform_id', 'sat_id', 'sensor_id', 'spectral_band', 'rttov_coeff_file']
    sensorDB = sensorDF.set_index('key').to_dict(orient='index')

    return sensorDB

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

def _ioda2iodaDF(iodaFile, obsCategoryConfig, mask_obsvalue=None, sort_by=None, derived=False):
    """Creates pandas dataframe from a IODA observation data file.

    Args:
        1. iodaFile - IODA observation data file
        2. obsCategoryConfig (dict): Observation category configuration
        3. mask_obsvalue (None, optional): Mask ObsValue _FillValue and NaN
             Set this to something like netCDF4.default_fillvals[dtype]
             where dtype can be 'f4' or 'f8' for floating point arrays
        4. sort_by (list of str, optional): Sort dataframe by columns
        5 derived: Set to true if using @DerivedObsValue

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

        # Grab the channel numbers if available
        channelNums = None
        if ('Channel' in ds.variables):
            channelNums = ds.variables['Channel'][:]

        # Grab the observation category config channel numbers
        obsChannels = obsCategoryConfig.get('channel numbers', [])

        loc = 0
        for var in dsets:
            fullVarName = var.group().name + "/" + var.name
            if (var.dimensions[0] != 'Location'):
                log.warning('Can only process variables with Location as the first dimension')
                log.warning('    Skipping variable: %s', fullVarName)
                continue

            if var.group().name != 'VarMetaData':
                if var.ndim == 2:
                    if var.name in ("dateTime", "stationIdentifier"):
                        outarr = ["".join(i) for i in var[:].astype(str)]
                        df.insert(loc, fullVarname, outarr)
                        loc += 1
                    # multi-channel BT
                    if var.name in ("brightnessTemperature", "emissivity"):
                        if (var.dimensions[1] != 'Channel'):
                            log.warning("Can only process variables with Channel as the second dimension")
                            log.warning("    Skipping variable: %s", fullVarName)
                            continue

                        # If obsChannels is empty, skip the variable
                        if not obsChannels:
                            log.warning("No observation channels specified")
                            log.warning("    Please add desired channels to the YAML configuration file.")
                            log.warning("    Skipping variable: %s", fullVarName)
                            continue

                        # Walk through the observation channels and form separate columns
                        # for each channel
                        for obsChannel in obsChannels:
                            # Find the index of the channel in the Channel dimension
                            if (obsChannel in channelNums):
                                channelIndex = np.where(channelNums == obsChannel)[0][0]
                            else:
                                log.warning("Channel number %s not found in Channel dimension", obsChannel)
                                log.warning("    Skipping channel number: %s", obsChannel)
                                continue

                            # Insert the column for the specific channel
                            fullVarNameChannel = f"{fullVarName}_{obsChannel}"
                            if fullVarNameChannel in df.columns:
                                log.warning("Column %s already exists")
                                log.warning("    Skipping column: %s", fullVarNameChannel)
                                continue

                            # Insert the data for the specific channel
                            outarr = var[:, channelIndex]
                            df.insert(loc, fullVarNameChannel, outarr)
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

def _buildRadianceMetadata(iodaDF, obsCategoryConfig, chanNum):
    """Build metadata for radiance observations.

    Args:
        1. iodaDF - pandas dataframe in the ioda layout
        2. obsCategoryConfig - Observation category configuration
        3. chanNum - Channel number

    Return:
        A list of metadata entries for the radiance observations.
    """

    # Grab the sensor database info
    sensorDbEntry = obsCategoryConfig['metadata']['sensor db entry']
    spectralBand = sensorDbEntry['spectral_band']
    platformId = sensorDbEntry['platform_id']
    satId = sensorDbEntry['sat_id']
    sensorId = sensorDbEntry['sensor_id']

    # Grab the instrument metadata variable names
    satAzVar = obsCategoryConfig['metadata']['sat az variable']
    satZeVar = obsCategoryConfig['metadata']['sat ze variable']

    # Form the line in the metadata containing the instrument info
    instrumentInfo = f"    {platformId} {satId} {sensorId} {chanNum}"

    # For now we only handle infrared and microwave instruments
    radMetaData = []
    if (spectralBand == 'ir'):
        # Format contains five strings:
        #     'visir'
        #     '    sat_az sat_ze sun_az sun_ze'
        #     '    {instrumentInfo}'
        #     '    specularity'
        #     '    oldkey'
        specularityInfo = f"    {missingVal}"
        for i, (satAz, satZe) in enumerate(zip(iodaDF[satAzVar], iodaDF[satZeVar])):
            orientationInfo = f"    {satAz} {satZe} {missingVal} {missingVal}"
            oldKeyInfo = f"    {i + 1}"
            radMetaData.append([ 'visr', orientationInfo, instrumentInfo, specularityInfo, oldKeyInfo ])
    elif (spectralBand == 'mw'):
        # Format contains six strings:
        #     'mw'
        #     '    sat_az sat_ze
        #     '    {instrumentInfo}'
        #     '    magfield cosbk'
        #     '    fastem_p1 fastem_p2 fastem_p3 fastem_p4 fastem_p5'
        #     '    oldkey'
        magfieldInfo = f"    {missingVal} {missingVal}"
        fastemInfo = f"    {missingVal} {missingVal} {missingVal} {missingVal} {missingVal}"
        for i, (satAz, satZe) in enumerate(zip(iodaDF[satAzVar], iodaDF[satZeVar])):
            orientationInfo = f"    {satAz} {satZe}"
            oldKeyInfo = f"    {i + 1}"
            radMetaData.append([ 'mw', orientationInfo, instrumentInfo, magfieldInfo, fastemInfo, oldKeyInfo ])
    else:
        raise ValueError("Unrecognized spectral band: " + spectralBand)

    return radMetaData

def _iodaDF2obsqDF(iodaDF, epochDT, iodaVarName, iodaVarType, obsCategoryConfig):
    """Reformat a pandas dataframe built with ioda2iodaDF into a dataframe suitable
       for buildObsSeqFromObsqDF.

     Args:
         1. iodaDF - pandas dataframe in the ioda layout
         2. epochDT - IODA dateTime epoch value in ISO 8601 string format
         3. iodaVarName - IODA dataframe column name
         4. iodaVarType - IODA dataframe column type
         5. obsCategoryConfig - Observation category configuration

     Return:
         This function returns a pandas dataframe in the obs_seq layout. This dataframe
         is suitable as input to the buildObsSeqFromObsqDF function.

    """

    # Set up obs category specific variables. At this point we only recognize
    # two categories:
    #     conventional
    #     radiance
    obsCategory = obsCategoryConfig['name']
    vertCoordName = None
    vertCoordUnits = None
    vertCoordValue = None
    if (obsCategory == 'conventional'):
        vertCoordName = obsCategoryConfig['vertical coordinate']['name']
        vertCoordUnits = obsCategoryConfig['vertical coordinate']['units']
    elif (obsCategory == 'radiance'):
        vertCoordUnits = obsCategoryConfig['vertical coordinate']['units']
        vertCoordValue = obsCategoryConfig['vertical coordinate']['data value']

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
    if obsCategory == 'conventional':
        obsqDF.insert(6, 'vertical', np.array(iodaDF[vertCoordName], dtype=np.float32))
        obsqDF.insert(7, 'vert_unit', np.full(numLocs, vertCoordUnits, dtype=object))
    elif obsCategory == 'radiance':
        obsqDF.insert(6, 'vertical', np.full(numLocs, vertCoordValue, dtype=np.float32))
        obsqDF.insert(7, 'vert_unit', np.full(numLocs, vertCoordUnits, dtype=object))

    obsqDF.insert(8, 'type', np.full(numLocs, iodaVarType, dtype=object))
    if (obsCategory == 'conventional'):
        obsqDF.insert(9, 'metadata', np.full(numLocs, '', dtype=object))
    elif (obsCategory == 'radiance'):
        chanNum = int(iodaVarName.split('_')[-1])
        radMetaData = _buildRadianceMetadata(iodaDF, obsCategoryConfig, chanNum)
        obsqDF.insert(9, 'metadata', radMetaData)
    obsqDF.insert(10, 'external_FO', np.full(numLocs, '', dtype=object))

    obsqDF.insert(11, 'seconds', obsqSec)
    obsqDF.insert(12, 'days', obsqDays)
    obsqDF.insert(13, 'time', obsqTime)

    obsqDF.insert(14, 'obs_err_var', np.array(iodaDF[obsErrorName], dtype=np.float32))

    # Remove rows with missing obs data, or missing location, timestamp data
    obsqDF.dropna(subset=['observation', 'longitude', 'latitude', 'vertical', 'seconds', 'days'], inplace=True)

    return obsqDF
