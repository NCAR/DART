#!/bin/bash

# 2/12/23
# Extracting site level grid cell from CAM reanalysis
# This site level extract will be fed to 'CAM4_reanalysis_biascorrection.m'
# bias corrected for site level met forcing

# CAM4 Reanalysis

CAM4_grid_path="/glade/collections/rda/data/ds199.1/"
CAM4_site_path="/glade/work/bmraczka/CAM4_NR1/"

# NR1 location (40.03, -105.55) or (40.03, 254.45)
# CAM6 grid is 1.25x0.95  or 288 longitude grids and 192 latitude grids
# This lat/lon corresponds with doma_lon= 204 (255 degrees east) ; doma_lat= 138  ; (40.0524 degrees)

# NR1 location (40.03, -105.55) or (40.03, 254.45)
# CAM6 grid is 1.95x2.5  or 288 longitude grids and 96 latitude grids
# This lat/lon corresponds with doma_lon= 102 (255 degrees east) ; doma_lat= 69  ; (40.7368 degrees)


for YEAR in {1997..2010..1}

do

echo "  "
echo "entering YEAR loop where value of YEAR is:"
echo ${YEAR}
echo "  "
   for NINST in {0001..0080..1}  # Will create 0001, 0002, ... 0080
   do
echo "  "
echo "entering NINST loop where value of YEAR is:"
echo ${NINST}
echo "  "  
 
   #Create output directory if necessary
   OUTDIR=${CAM4_site_path}${NINST}
   if [[ ! -d ${OUTDIR} ]] ; then
        mkdir ${OUTDIR}
   fi
  
   
   ncks -d doma_nx,102,102 -d doma_ny,69,69 -d a2x6h_nx,102,102 -d a2x6h_ny,69,69 \
         ${CAM4_grid_path}/CAM_DATM.cpl_${NINST}.ha2x1dx6h.${YEAR}.nc \
         ${CAM4_site_path}${NINST}/CAM4_NR1.cpl_${NINST}.ha2x1dx6h.${YEAR}.nc
 
   echo "  "
   echo "Completed CAM4 NR1 extract"
   echo "  "

   done
done


exit 0


