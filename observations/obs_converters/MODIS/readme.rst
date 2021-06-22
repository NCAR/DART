DART observations and MODIS products.
=====================================

There are many MODIS products, in many formats. This document will list all of the data products and formats that have
DART programs to convert them to observation sequence files.

Programs
--------

+-------------------------+-------------------------------------------------------------------------------------------+
| :doc:`./MOD15A2_to_obs` | Converts `MODIS Land Product Subsets <http://daac.ornl.gov/MODIS/modis.shtml>`__ Leaf     |
|                         | Area Index (**LAI**) and Fraction of Photosynthetically Active Radiation (**FPAR**) 8 day |
|                         | composite `[MOD15A2] <https://lpdaac.usgs.gov/products/modis_products_table/mod15a2>`__   |
+-------------------------+-------------------------------------------------------------------------------------------+

Plans
-----

#. Support MOD15A2 'Global Tool' records.
#. The work that remains is to get the IGBP landcover code for the site and incorporate that into the observation
   metadata. I *almost* have everything I need. Once that happens, the forward observation operator can be made to be
   much more accurate by only using model landunits that have the right landcover class.
#. Support more products. Put in a request to help me prioritize.
