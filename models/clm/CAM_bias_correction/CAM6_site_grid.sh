#!/bin/bash

# Extracts site level grid cell from CAM6 reanalysis

# CAM6 Reanalysis

CAM6_grid_path="/glade/collections/rda/data/ds345.0/cpl_unzipped/"
# Site level  grid cell from CAM6 reanalysis
CAM6_site_path="/glade/work/bmraczka/CAM6_NR1/"

# The US-NR1 flux tower site is used as an example 
# NR1 location (40.03, -105.55) or (40.03, 254.45)
# CAM6 grid is 1.25x0.95  or 288 longitude grid points and 192 latitude grid points
# This lat/lon corresponds with doma_lon= 204 (255 degrees east) ; doma_lat= 138  ; (40.0524 degrees)

for YEAR in {2011..2020..1}
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
      OUTDIR=${CAM6_site_path}${NINST}
      if [[ ! -d ${OUTDIR} ]] ; then
           mkdir ${OUTDIR}
      fi


      # SOLAR
      ncks -d doma_nx,204,204 -d doma_ny,138,138 -d a2x1hi_nx,204,204 -d a2x1hi_ny,138,138 \
            ${CAM6_grid_path}${NINST}/f.e21.FHIST_BGC.f09_025.CAM6assim.011.cpl_${NINST}.ha2x1hi.${YEAR}.nc \
            ${CAM6_site_path}${NINST}/CAM6_NR1.cpl_${NINST}.ha2x1hi.${YEAR}.nc
      echo "  "
      echo "Finished CAM6 SOLAR, 1/3 completed"
      echo "  "

      # NON-SOLAR                                       
      ncks -d doma_nx,204,204 -d doma_ny,138,138 -d a2x3h_nx,204,204 -d a2x3h_ny,138,138 \
            ${CAM6_grid_path}${NINST}/f.e21.FHIST_BGC.f09_025.CAM6assim.011.cpl_${NINST}.ha2x3h.${YEAR}.nc \
            ${CAM6_site_path}${NINST}/CAM6_NR1.cpl_${NINST}.ha2x3h.${YEAR}.nc

      echo "  "
      echo "Finished CAM6 NON-SOLAR, 2/3 completed"
      echo "  "
      # 1 hr STATE
      ncks -d doma_nx,204,204 -d doma_ny,138,138 -d a2x1h_nx,204,204 -d a2x1h_ny,138,138 \
            ${CAM6_grid_path}${NINST}/f.e21.FHIST_BGC.f09_025.CAM6assim.011.cpl_${NINST}.ha2x1h.${YEAR}.nc \
            ${CAM6_site_path}${NINST}/CAM6_NR1.cpl_${NINST}.ha2x1h.${YEAR}.nc
      echo "  "
      echo "Finished CAM6 1 hr STATE, 3/3 completed"
      echo "  "

   done
done


exit 0


