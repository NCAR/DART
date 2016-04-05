#!/bin/bash

source ~/ncScripts/ncFilters.sh ## will have to eventually distribte these nco "filters".

ncVarDiff sh2ox Prior_Diag.nc validation/Prior_Diag.nc | grep sh2ox
ncVarDiff sh2ox Posterior_Diag.nc validation/Posterior_Diag.nc | grep sh2ox

exit 0
