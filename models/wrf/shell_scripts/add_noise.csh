#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

# echo calling add_noise
# ../add_noise.csh $wrfsecs $wrfdays $state_copy $ensemble_member $temp_dir $CENTRALDIR
# The following are passed from the advance model script:
set wrfsecs = $1
set wrfdays = $2
set state_copy = $3
set ensemble_member = $4
set temp_dir = $5
set CENTRALDIR = $6
#
# notes - do we need to do a 'touch' here for files that don't yet exist?
# The remaining lines in this section control the additive noise (Dowell et al. 2008).
#
# if run_add_pert is 1, then grid_refl_obs and add_pert_where_high_refl will be run
# at intervals of time_between_perturbations (in seconds)
# other definitions that don't need to be passed:

set MOVE = '/bin/mv -f'
set  REMOVE = '/bin/rm -rf'

set run_add_pert = 0
set time_between_perturbations = 300

set refl_min = 25.0                 # minimum reflectivity threshold (dBZ) for where noise is added

set lh = 4000.0                     # horizontal length scale (m) for perturbations
set lv = 2000.0                     # vertical length scale (m) for perturbations

set transition_time_days = 145917   # gregorian day when transition occurs from first perturbation magnitude to second
set transition_time_secs = 83400    # seconds in day when transition occurs from first perturbation magnitude to second

set u_noise_1 = 1.0                  # std. dev. of u noise (m/s), before smoothing, during first time period
set v_noise_1 = 1.0                  # std. dev. of v noise (m/s), before smoothing, during first time period
set w_noise_1 = 0.0                  # std. dev. of w noise (m/s), before smoothing, during first time period
set t_noise_1 = 1.0                  # std. dev. of temperature noise (K), before smoothing, during first time period
set td_noise_1 = 1.0                 # std. dev. of dewpoint noise (K), before smoothing, during first time period
set qv_noise_1 = 0.0                 # std. dev. of water vapor noise (g/kg), before smoothing, during first time period

set u_noise_2 = 0.25                 # std. dev. of u noise (m/s), before smoothing, during second time period
set v_noise_2 = 0.25                 # std. dev. of v noise (m/s), before smoothing, during second time period
set w_noise_2 = 0.0                  # std. dev. of w noise (m/s), before smoothing, during second time period
set t_noise_2 = 0.25                 # std. dev. of temperature noise (K), before smoothing, during second time period
set td_noise_2 = 0.25                # std. dev. of dewpoint noise (K), before smoothing, during second time period
set qv_noise_2 = 0.0                 # std. dev. of water vapor noise (g/kg), before smoothing, during second time period

#
#****************************** END OF USER SPECIFICATIONS *********************************
###################################################################################################
#
# This code section is for the additive noise (Dowell et al. 2008).
# Run grid_refl_obs and add_pert_where_high_refl, if required.
#
#
# GSR passed information needed below: wrfsecs, wrfdays
   set need_to_run_grid_refl_obs = 0
   set perturbations_needed = 0

   if ($run_add_pert == 1) then
     if (-e ../last_perturbation_time_${ensemble_member}.txt) then
       set secday = `head -1 ../last_perturbation_time_${ensemble_member}.txt`
       set pertsecs = $secday[1]
       set pertdays = $secday[2]
       ${REMOVE} ../last_perturbation_time_${ensemble_member}.txt
     else
       set pertsecs = $wrfsecs
       set pertdays = $wrfdays
       @ pertsecs = $pertsecs - $time_between_perturbations
     endif

     @ time_since_last_perturbation = 86400 * ($wrfdays - $pertdays) + ($wrfsecs - $pertsecs)

     echo "wrfdays wrfsecs = ${wrfdays} ${wrfsecs}"
     echo "time_since_last_perturbation = ${time_since_last_perturbation}"
     echo "time_between_perturbations = ${time_between_perturbations}"

     if ($time_since_last_perturbation >= $time_between_perturbations ) then
       set perturbations_needed = 1
       if ($state_copy == 1) then
         set need_to_run_grid_refl_obs = 1
       endif
     endif

   endif

   if ( ($need_to_run_grid_refl_obs == 1) && ($ensemble_member == 1) ) then

     cd ..

     set days_end = $wrfdays
     set seconds_end = $wrfsecs
     set days_begin = $days_end
     @ seconds_begin = $seconds_end - $time_between_perturbations

     echo "days_end seconds_end days_begin seconds_begin = ${days_end} ${seconds_end} ${days_begin} ${seconds_begin}"

     while ($seconds_begin < 0)
       @ days_begin = $days_begin - 1
       @ seconds_begin = $seconds_begin + 86400
     end     

     echo "member ${ensemble_member} calling grid_refl_obs"
     grid_refl_obs obs_seq.out $refl_min $days_begin $seconds_begin $days_end $seconds_end wrfinput_d01
     
     ${MOVE} refl_obs.txt refl_obs_${wrfsecs}_${wrfdays}.txt
     cat > finished.grid_refl_obs_${wrfsecs}_${wrfdays} <<EOF
grid_refl_obs has finished
EOF
     
     cd $temp_dir

   endif

   if ($perturbations_needed == 1) then

     while (! -e ../finished.grid_refl_obs_${wrfsecs}_${wrfdays})
       echo "member ${ensemble_member} waiting for ../finished.grid_refl_obs_${wrfsecs}_${wrfdays} to appear"
       sleep 2
     end   

     @ transition_relative_time = 86400 * ($wrfdays - $transition_time_days) + ($wrfsecs - $transition_time_secs)

     if ($transition_relative_time <= 0) then
       echo "(1) member ${ensemble_member} calling add_pert_where_high_refl"
       ../add_pert_where_high_refl ../refl_obs_${wrfsecs}_${wrfdays}.txt wrfinput_d01 $lh $lv $u_noise_1 $v_noise_1 $w_noise_1 $t_noise_1 $td_noise_1 $qv_noise_1 $ensemble_member $wrfdays $wrfsecs
     else
       echo "(2) member ${ensemble_member} calling add_pert_where_high_refl"
       ../add_pert_where_high_refl ../refl_obs_${wrfsecs}_${wrfdays}.txt wrfinput_d01 $lh $lv $u_noise_2 $v_noise_2 $w_noise_2 $t_noise_2 $td_noise_2 $qv_noise_2 $ensemble_member $wrfdays $wrfsecs
     endif


     cd ${CENTRALDIR}
     cat > last_perturbation_time_${ensemble_member}.txt <<EOF
     ${wrfsecs} ${wrfdays}
EOF
     cd $temp_dir

   endif
#
#
# end of additive noise section
#
###################################################################################################
echo 'add_noise finished'

exit $status

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

