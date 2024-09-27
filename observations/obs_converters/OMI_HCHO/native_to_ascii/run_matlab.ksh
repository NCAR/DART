#!/bin/ksh -aeux
. /etc/profile.d/lmod.sh
module load matlab
matlab -nosplash -nodesktop -r 'omi_no2_extract(2014,7,13,21,0,0,2014,7,14,3,0,0,0.,360.,-90.,90.)'
exit

