import geopandas as gp
import matplotlib as mp
mp.use('TkAgg')  
import matplotlib.pyplot as plt

def plot_channel_network(shape_file, value_frame=None, column_name=None, **plot_args):
    """
    Plot the channel network with an optionally supplied value dataframe.
    
    Keyword Arguments:
    shapefile -- A path/file to use for the domain.
    value_frame -- A pandas dataframe with feautre_id column named "ID"  and a value column
                   to be joined to the shapefile.

    Example: 

    """

    network = gp.GeoDataFrame.from_file(shape_file)

    if(value_frame is not None):
        network = gp.pd.merge(network, value_frame, on=('ID'), how='inner')
        if(column_name is None):
            if value_frame.shape[1] != 2:
                print('More than 2 columns in value_frame and column_name not set.')
                return None
            column_name = list(set(value_frame.columns.values).difference({"ID"}))[0]
        network.plot(column=column_name, **plot_args)
    else: 
        network.plot(**plot_args)

    return plt
