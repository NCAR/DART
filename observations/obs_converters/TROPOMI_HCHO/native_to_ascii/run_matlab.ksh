#!/bin/ksh -aeux
. /etc/profile.d/lmod.sh
module load matlab
matlab -nosplash -nodesktop -r 'tropomi_no2_extract(2019,7,12,3,0,0,2019,7,13,3,0,0,0.,360.,-90.,90.)'
exit

