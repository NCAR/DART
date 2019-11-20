import geopandas as gp
#import matplotlib as mp
#mp.use('TkAgg')  
#import matplotlib.pyplot as plt
#import shapely
import xarray as xa

os.chdir('python')
from plot_channel_network import *

## The shape file of the domain
network_file = '/Users/jamesmcc/Downloads/sixmile_NY_shapefiles/'

## The data to join to the shape file by feature_id (NHD comid) but named as 'ID'.
rl_file = '/Users/jamesmcc/Downloads/sixmile_test_domain/DOMAIN/RouteLink.nc'
rl_frame = xa.open_dataset(rl_file).to_dataframe()
rl_frame = rl_frame.rename(index=str, columns={ 'link' : 'ID' })
rl_subframe = rl_frame[ ['ID', 'order']]
rl_subframe3 = rl_frame[ ['ID', 'order', 'BtmWdth']]

rl_gage = rl_frame

rl_gage=rl_frame.loc[rl_frame.gages != b'               ', ['gages','lat','lon']]
rl_gage['gages']=rl_gage['gages'].str.decode("utf-8")
rl_gage = rl_gage.reset_index()

#######################################################
flowlines_file = '/Users/jamesmcc/Downloads/sixmile_NY_shapefiles/flowlines.shp'
flowlines = gp.GeoDataFrame.from_file(flowlines_file)
flow = flowlines.query('ID=="*"')
from plotnine import *

(ggplot() +
 geom_map(flowlines, fill=None))


#######################################################

+ geom_text(
     political,
     aes('geometry.centroid.x', 'geometry.centroid.y', label='fmt_labels(name, ClaimedBy)'),
     size=8,
     fontweight='bold'
 )
 + geom_text(
     cities,
     aes('geometry.centroid.x', 'geometry.centroid.y', label='name'),
     size=8,
     ha='left',
     nudge_x=.20
 )
 + labs(title="The Political Territories of Westeros")
 + scale_fill_brewer(type='qual', palette=8)
 + scale_x_continuous(expand=(0, 0, 0, 1))
 + scale_y_continuous(expand=(0, 1, 0, 0))
 + scale_size_continuous(range=(0.4, 1))
 + theme_void()
 + theme(figure_size=(8, 12), panel_background=element_rect(fill=water_color))
)







#######################################################
plot = plot_channel_network( network_file, value_frame=rl_subframe3, column_name='BtmWdth', \
                             legend=True, cmap='viridis' )
plot.scatter(rl_gage['lon'], \
             rl_gage['lat'], \
             label=rl_gage['gages'], \
             color=range(rl_gage.shape[0]), \ 
             alpha=0.9)
plot.legend()
plot.show()


plot = plot_channel_network( network_file, rl_subframe, \
                             legend=True, cmap='gist_earth',
                             title='Sixmile, near Ithaca NY')
ax=plot.gca()
ax.set_facecolor('lightgrey')
plot.show()

plot = plot_channel_network( network_file )
plot.show()


