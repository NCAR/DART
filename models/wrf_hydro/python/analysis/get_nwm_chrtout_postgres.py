import pandas as pd
from sqlalchemy import create_engine

def get_nwm_chrtout(nwm_version, forecast_cycle, site_nos_list):
    """Get NWM streamflow from databse by version and site.
    
    Keyword Arguments:
    nwm_version -- one of 'v_1_2_1'
    forecast_cycle -- one of 'analysis_assim', 'short_range', 'medium_range', long_range'
    site_nos_list -- list of multiple gage identifiers from the RouteLink file, typically a USGS identifier.
    """

    # TODO: more sanity checking on arguments

    #Correct version number for v1_2 vs v1_2_1.
    if nwm_version == 'v1_2':
        nwm_version = 'v1_2_1'

    # Construct engine and query strings
    engineString = 'postgresql://nwmload:nwmIngest@hydro-c1-web/' + nwm_version # database credentials
    tableName = forecast_cycle + '_channel'

    sites_list=list(map(str.strip,site_nos_list))
    sites_str = '\',\''.join(sites_list)


    channelQuery = "select * from " + tableName + \
                  " RIGHT JOIN gage_link ON " + \
                  tableName + ".feature_id=gage_link.feature_id where gage_link.site_no IN ('" + sites_str + "')"

    #channelQuery = "select * from " + tableName + \
    #              " RIGHT JOIN gage_link WHERE gage_link.site_no IN ('" + sites_str + "')"

    # Create an engine
    engine = create_engine(engineString)

    # Get the data
    channel_data = pd.read_sql_query(channelQuery, con=engine)
    channel_data = channel_data.rename(index=str,columns={'streamflow':'streamflow_mod'})
    
    return channel_data

