#!/bin/bash
#
# This code may (or may not) be part of the COAMPS distribution,
# So it is not protected by the DART copyright agreement.
#
# DART $Id$
#
########################################################################
#
# Helper functions for main script
#
########################################################################
# Filtered directory listings to output file.
function set_field_list
{
  usage="Usage: `basename $0` outfile dtg targetDir [awk_cmd]"
  outfile=${1} ; dtg_loc=${2} ; targDir=${3} ; awk_cmd=${4}
  if [ -f  ${outfile} ] ; then rm -f ${outfile}; fi
  echo ${targDir} >  ${outfile}
  echo ${dtg_loc} >> ${outfile}
  if [ ! -d ${targDir} ] ; then return; fi
  if [ -n "${awk_cmd}" ]; then
    ls -1 ${targDir} | grep `echo ${dtg_loc}` | ${awk_cmd} >> ${outfile}
  else
    ls -1 ${targDir} | grep `echo ${dtg_loc}` >> ${outfile}
  fi
  return
}
function check_return
{
   if [ ${1} != 0 ]; then
    echo "error in ${2}"
    exit ${1}
   fi
}
function add_trailing_slash
{
	if [ ! -n "$1" ]; then return; fi
	if [ ${1: (-1)} == \/ ]; then echo ${1}; else echo ${1}/; fi
}

usage="Usage: `basename $0` -d dtg -i icycle -m member -n cycle_num -c config_file -h help"

# Parse the options and grab the file name
while getopts ":d:i:m:n:c:h:" option; do
  case $option in
    d  ) dtg=$OPTARG;;
	i  ) CYCLE_HR=$OPTARG;;
	m  ) member=$OPTARG;;
	n  ) CYCLE_NO=$OPTARG;;
	c  ) PATH_CONFIG=$OPTARG;;
	h  ) echo $usage
         exit 0;;
	\? ) echo $usage
	     exit 65;;
	*  ) echo $usage
         exit 65;;
  esac
done

PERL='perl -i -p -e'
SED='sed -i -e'

########################################################################
# Set default parameter values
########################################################################
: ${PATH_CONFIG:=CONFIG/FILE/GOES/HERE}
. ${PATH_CONFIG}
: ${CYCLE_HR:=$icycle}
: ${dtg:=$cdtg_beg}
: ${member:=mean}
: ${CYCLE_NO=2}
: ${coamps_nml=namelist.coa}

########################################################################
#   Path locations needed for this script
########################################################################
WORKDIR=`pwd`
NAVDAS_BIN=${NAVDAS_HOME}/bin
NAVDAS_AWK=${NAVDAS_HOME}/src/scripts/awk
NAVDAS_RUN=${WORKDIR}/obs_preproc
NAVDAS_LOG=${NAVDAS_RUN}/log

if [ ! -d ${NAVDAS_RUN} ] ; then mkdir ${NAVDAS_RUN} ; fi
if [ ! -d ${NAVDAS_LOG} ] ; then mkdir ${NAVDAS_LOG} ; fi
cd ${NAVDAS_RUN}

########################################################################
#  setup date-time-groups
########################################################################
dtgbk=`echo ${dtg}\ -${CYCLE_HR} | $NAVDAS_BIN/newdtg_nav.exe`
check_return $? newdtg_nav.exe

dtgbk2=`echo ${dtgbk}\ -${CYCLE_HR} | $NAVDAS_BIN/newdtg_nav.exe`
check_return $? newdtg_nav.exe

CRDATE=${dtg}
CDTG_AN=${dtg}

cp ${DART_BASE}/namelist ./${coamps_nml}
NESTS_AN=`awk -f ${NAVDAS_AWK}/grid_num.awk ${coamps_nml} | tr -d ','`

########################################################################
#   Path locations needed for navdas
########################################################################
LANDSEA=${NAVDAS_DATABASE}/etc/landseatable
INNOVATIONS=${WORKDIR}/innovations
OBSERVATIONS=`printf ${OBS_PATH} ${dtg}`
NOGAPSPATH=${NOGAPS_PATH}

if [[ "$member" =~ ^[0-9]+$ ]]; then 
  BACKGROUND=`printf ${COAMPS_DATA} ${member}`
else
  BACKGROUND=`printf ${COAMPS_MEAN_DATA} ${member}`
fi
BACKGROUND2=`printf ${NOGAPS_DETERMIN} ${dtg}`
if [ ${CYCLE_NO} -ge 1 ]; then
 BACKGROUND3=`printf ${NOGAPS_DETERMIN} ${dtgbk}`
else
 BACKGROUND3=`printf ${NOGAPS_DETERMIN} ${dtg}`
fi
TFILE_OBS=$OBSERVATIONS
SATW_OBS=$OBSERVATIONS
TOVS_OBS=$OBSERVATIONS
TC_BOGUS_OBS=$OBSERVATIONS

########################################################################
#  Update namelist variables from generic values   
########################################################################
echo "Updating COAMPS namelist for observation preproc"
  ${PERL} "s:^(\s*cdtg\s*=\s*).\d+.(,).*:\1'${dtg}'\2:" ${coamps_nml}
  ${PERL} "s:^(\s*dsngff\s*=\s*\').*(\',):\1${BACKGROUND2}\2:" ${coamps_nml}
  ${PERL} "s:^(\s*dsnrff\s*=\s*\').*(\',):\1${BACKGROUND}\2:" ${coamps_nml}
  ${SED} "s/^\(\s*icycle\s*=\s*\).*,\s*/\1${CYCLE_HR},/" ${coamps_nml}

########################################################################
#  make directories if they do not exist   
########################################################################
  if [ ! -d ${INNOVATIONS} ]           ; then mkdir ${INNOVATIONS} ; fi
  if [ ! -d ${INNOVATIONS}/acft ]      ; then mkdir ${INNOVATIONS}/acft ; fi
  if [ ! -d ${INNOVATIONS}/satw ]      ; then mkdir ${INNOVATIONS}/satw ; fi
  if [ ! -d ${INNOVATIONS}/grid ]      ; then mkdir ${INNOVATIONS}/grid ; fi
  if [ ! -d ${BACKGROUND} ]            ; then mkdir ${BACKGROUND} ; fi
  if [ ! -d ${BACKGROUND}/background ] ; then mkdir ${BACKGROUND}/background ; fi
  if [ ! -d ${BACKGROUND}/nogaps ]     ; then mkdir ${BACKGROUND}/nogaps ; fi

###################################################
##    list of observation variables              ##
###################################################
#    c_typ_rej  == variable types  
# 
#  c_var_rej(1) ='Z',
#  c_var_rej(2) ='T',
#  c_var_rej(3) ='u',
#  c_var_rej(4) ='v',
#  c_var_rej(5) ='prh',
#  c_var_rej(6) ='o3',
#  c_var_rej(7) ='dd',
#  c_var_rej(8) ='ff',
#  c_var_rej(9) ='thk',
#  c_var_rej(10)='P_msl',
#  c_var_rej(11)='P_stn',
#  c_var_rej(12)='T_dpd',
#  c_var_rej(13)='bT',
#  c_var_rej(14)='TPPW',
#  c_var_rej(15)='q',
 
###################################################
#   list observation variable parameters to	 ##
#   be excluded from this analysis      	 ##
###################################################

#    definition of rejected observations 
#    order is not important in definition
#    of c_var_rej....examples are above  

  cat << eof > ${INNOVATIONS}/var_rej.${dtg}
 &var_rej

  c_var_rej(1)='',
  c_var_rej(2)='',
  c_var_rej(3)='',
  c_var_rej(4)='',
  c_var_rej(5)='',
  c_var_rej(6)='',
  c_var_rej(7)='',
  c_var_rej(8) ='',
  c_var_rej(9) ='',
  c_var_rej(10)='',
  c_var_rej(11)='',
  c_var_rej(12)='',
  c_var_rej(13)='',
  c_var_rej(14)='',
  c_var_rej(15)='',

 /
eof

###################################################
##    list instrument types                      ##
###################################################
 
  cat << eof > ${INNOVATIONS}/full_list.${dtg}
 &typ_rej
 c_typ_rej(1) ='sfc land'	 ,
 c_typ_rej(2) ='sfc ship'	 ,
 c_typ_rej(3) ='man-airep'	 ,
 c_typ_rej(4) ='man-Yairep'	 ,
 c_typ_rej(5) ='airep'		 ,
 c_typ_rej(6) ='airep_asc'	 ,
 c_typ_rej(7) ='airep_des'	 ,
 c_typ_rej(8) ='airep_lvl'	 ,
 c_typ_rej(9) ='airep_msg'	 ,
 c_typ_rej(10)='amdar'		 ,
 c_typ_rej(11)='amdar_asc'	 ,
 c_typ_rej(12)='amdar_des'	 ,
 c_typ_rej(13)='amdar_lvl'	 ,
 c_typ_rej(14)='acars'		 ,
 c_typ_rej(15)='acars_asc'	 ,
 c_typ_rej(16)='acars_des'	 ,
 c_typ_rej(17)='acars_lvl'	 ,
 c_typ_rej(18)='mdcrs'		 ,
 c_typ_rej(19)='mdcrs_asc'	 ,
 c_typ_rej(20)='mdcrs_des'	 ,
 c_typ_rej(21)='mdcrs_lvl'	 ,
 c_typ_rej(22)='cld wnds1'	 ,
 c_typ_rej(23)='cld wnds2'	 ,
 c_typ_rej(24)='METEO-5A'	 ,
 c_typ_rej(25)='METEO-7A'	 ,
 c_typ_rej(26)='GOES-10'	 ,
 c_typ_rej(27)='GOES-8'		 ,
 c_typ_rej(28)='GMSAFW'		 ,
 c_typ_rej(29)='GMSC'		 ,
 c_typ_rej(30)='ssmi ff1'	 ,
 c_typ_rej(31)='ssmi ff2'	 ,
 c_typ_rej(32)='scat winds'	 ,
 c_typ_rej(33)='Aus synth'	 ,
 c_typ_rej(34)='raob'		 ,
 c_typ_rej(35)='pibal'		 ,
 c_typ_rej(36)='analytic'	 ,
 c_typ_rej(37)='amsub qprf'	 ,
 c_typ_rej(38)='tovs T'		 ,
 c_typ_rej(39)='ssmi prh'	 ,
 c_typ_rej(40)='TC synth'	 ,
 c_typ_rej(41)='atovs bT'	 ,
 c_typ_rej(42)='rtovs bT'	 ,
 c_typ_rej(43)='ssmt1 bT'	 ,
 c_typ_rej(44)='ssmt2 bT'	 ,
 c_typ_rej(45)='ssmi TPPW'	 ,
 c_typ_rej(46)='WSat TPPW'       ,
 c_typ_rej(50)='METEO-5C'	 ,
 c_typ_rej(51)='METEO-7C'	 ,
 c_typ_rej(52)='GOES-W'		 ,
 c_typ_rej(53)='GOES-E'		 ,
 c_typ_rej(54)='GOES-11'	 ,
 c_typ_rej(55)='GOES-12'	 ,
 c_typ_rej(56)='GMSX'		 ,
 c_typ_rej(60)='raob_sig'		 ,
 /
eof
 
###############################################################
##    list instruments to be excluded from this analysis     ##
###############################################################

# c_typ_rej(1) ='sfc land'	 ,
# c_typ_rej(2) ='sfc ship'	 ,
  cat << eof > ${INNOVATIONS}/typ_rej.${dtg}
 &typ_rej
 c_typ_rej(3) ='man-airep'       ,
 c_typ_rej(4) ='man-Yairep'      ,
 c_typ_rej(7) ='airep_des'       ,
 c_typ_rej(12)='amdar_des'       ,
 c_typ_rej(16)='acars_des'       ,
 c_typ_rej(20)='mdcrs_des'       ,
 c_typ_rej(33)='Aus synth'       ,
 c_typ_rej(37)='amsub qprf'	 ,
 c_typ_rej(39)='ssmi prh'	 ,
 c_typ_rej(43)='ssmt1 bT'	 ,
 c_typ_rej(44)='ssmt2 bT'	 ,
 /
eof
###########################################################
##    list of changes to navdas parameter switches       ##
###########################################################

if [ -f  navdas_user_${dtg} ] ; then rm -f navdas_user_${dtg}; fi
cat << eof > navdas_user_${dtg}
 &navdas_user
 /
eof
cp -f navdas_user_${dtg} ${NAVDAS_LOG}/navdas_user_${dtg}

########################################################################
#  export environtmental variables that navdas needs (Grumble, grumble).
#  Be sure to add trailing slashes to paths.
########################################################################
export dtg
export CRDATE
export CDTG_AN
export CYCLE_HR
export CYCLE_NO
export NESTS_AN
export RUNTYPE='unclassified'

export OBSERVATIONS=`add_trailing_slash $OBSERVATIONS`
export TFILE_OBS=`add_trailing_slash $TFILE_OBS`
export TC_BOGUS_OBS=`add_trailing_slash $TC_BOGUS_OBS`
export SATW_OBS=`add_trailing_slash $SATW_OBS`
export TOVS_OBS=`add_trailing_slash $TOVS_OBS`

export LANDSEA=`add_trailing_slash ${LANDSEA}`
export INNOVATIONS=`add_trailing_slash $INNOVATIONS`
export BACKGROUND=`add_trailing_slash ${BACKGROUND}`
export BACKGROUND2=`add_trailing_slash ${BACKGROUND2}`
export BACKGROUND3=`add_trailing_slash ${BACKGROUND3}`
export NOGAPSPATH=`add_trailing_slash ${NOGAPSPATH}`
export RUNOUTPUT=`add_trailing_slash ${NAVDAS_LOG}`

########################################################################
#  setup coamps grid for navdas 
########################################################################

#  cleanup and files that may exists from a previous execution
   for i in a b `seq 0 ${NESTS_AN}`; do
     if [ -f  parm_${i}_${dtg} ] ; then rm -f parm_${i}_${dtg}; fi
     if [ -f  grid_${i}_${dtg} ] ; then rm -f grid_${i}_${dtg}; fi
   done
   if [ -f  parm_${dtg} ]           ; then rm -f parm_$dtg; fi
   if [ -f  coamps_grid_${dtg}.nl ] ; then rm -f coamps_grid_${dtg}.nl; fi

#  filter out namelist variables unknown to navdas 
   if [ -f  namelist ] ; then rm -f namelist; fi
   $NAVDAS_AWK/namelist_grid_nav.awk  < ${coamps_nml} >  namelist
   $NAVDAS_AWK/namelist_atmos_nav.awk < ${coamps_nml} >> namelist
   $NAVDAS_AWK/namelist_dset_nav.awk  < ${coamps_nml} >> namelist
   $NAVDAS_AWK/namelist_coam_nav.awk  < ${coamps_nml} >> namelist
   cp ./namelist ./navdas_namelist_${dtg}

   ${NAVDAS_BIN}/setup_coamps.exe > ${NAVDAS_LOG}/setup_coamps_${NESTS_AN}.$dtg
   check_return $? setup_coamps.exe
#
########################################################################
#  Coamps background file names  
########################################################################
   set_field_list bk_fld_list.${dtg}    ${dtg}    ${BACKGROUND} ${NAVDAS_AWK}/find_bck_flds.awk
   set_field_list bk_fld_list.${dtgbk}  ${dtgbk}  ${BACKGROUND} ${NAVDAS_AWK}/find_fcst_coflds.awk
   set_field_list bk_fld_list.${dtgbk2} ${dtgbk2} ${BACKGROUND} ${NAVDAS_AWK}/find_fcst_coflds.awk
   coamps_ftype=`sed -n '3p' bk_fld_list.${dtgbk} | wc -c`

  if [ $CYCLE_NO = "1" ] ; then 
     echo "BACKGROUND3 set to ${dtg} analysis: $BACKGROUND3"
     echo "BACKGROUND4 set to ${dtg} analysis: $BACKGROUND4"
  else

# check if coamps background file list empty
  numlines=` wc -l < bk_fld_list.${dtgbk} `
  if [ $numlines -lt 69 ]; then
     echo "coamps background list bk_fld_list.${dtgbk} is empty"
     numlines=` wc -l < bk_fld_list.${dtgbk2} `
     if [ $numlines -lt 69 ]; then
        echo "coamps background list bk_fld_list.${dtgbk2} is empty"
        if [ $CYCLE_NO -ge 2 ]; then
          export CYCLE_NO=1
          export BACKGROUND3=`add_trailing_slash $BACKGROUND2`
          export BACKGROUND4=`add_trailing_slash $BACKGROUND2`
          echo "reset CYCLE_NO to $CYCLE_NO"
          echo "reset BACKGROUND3 to ${dtg} analysis: $BACKGROUND3"
          echo "reset BACKGROUND4 to ${dtg} analysis: $BACKGROUND4"
        else
          echo "BACKGROUND3 set to ${dtg} analysis: $BACKGROUND3"
          echo "BACKGROUND4 set to ${dtg} analysis: $BACKGROUND4"
        fi
     else
        echo "coamps 6 hr background fields missing ..... reset for 12 hr !!! "
        export BACKGROUND3=`add_trailing_slash $BACKGROUND4`
        export BACKGROUND4=`add_trailing_slash $BACKGROUND4`
        echo "reset BACKGROUND3 to ${dtgbk2} forecast: $BACKGROUND3"
     fi
  fi

  fi

########################################################################
#  boundary conditions  file names
########################################################################
 
   set_field_list bk_fld_list2.${dtg} ${dtg} $BACKGROUND2

   if [ ${CYCLE_NO} -ge 2 ]; then dtg_tmp=${dtgbk}; else dtg_tmp=${dtg}; fi
   set_field_list bk_fld_list3.${dtg} ${dtg_tmp} $BACKGROUND3

   if [ ${CYCLE_NO} -ge 2 ]; then dtg_tmp=${dtgbk2}; else dtg_tmp=${dtg}; fi
   set_field_list bk_fld_list4.${dtg} ${dtg_tmp} $BACKGROUND4
   nogaps_ftype=`sed -n '3p' bk_fld_list3.${dtg} | wc -c`

#  check if nogaps background file list empty
   numlines=` wc -l < bk_fld_list3.${dtg} `
   if [ $numlines -lt 69 ]; then
      echo "nogaps $BACKGROUND3 background fields missing....."
      numlines=` wc -l < bk_fld_list4.${dtg} `
      if [ $numlines -lt 69 ]; then
        echo "nogaps $BACKGROUND4 background fields missing....."
        export CYCLE_NO=1
        export BACKGROUND3=`add_trailing_slash $BACKGROUND2`
        export BACKGROUND4=`add_trailing_slash $BACKGROUND2`
        echo "reset CYCLE_NO to $CYCLE_NO"
        echo "reset BACKGROUND3 to ${dtg} analysis: $BACKGROUND3"
        echo "reset BACKGROUND4 to ${dtg} analysis: $BACKGROUND3"
        cp -f bk_fld_list2.${dtg} bk_fld_list3.${dtg}
        cp -f bk_fld_list2.${dtg} bk_fld_list4.${dtg}

        numlines=` wc -l < bk_fld_list2.${dtg} `
        nogaps_ftype=`sed -n '3p' bk_fld_list2.${dtg} | wc -c`
        if [ $numlines -lt 69 ]; then
          echo " .... missing $BACKGROUND2  analysis"
          echo " .... abort run !!!"
          exit 99
        fi
      else
        echo "reset BACKGROUND3 to $BACKGROUND4"
        BACKGROUND3=$BACKGROUND4
        echo "reset BACKGROUND3 to $BACKGROUND3"
        cp -f bk_fld_list4.${dtg} bk_fld_list3.${dtg}
      fi
   fi

#  check if nogaps background file lists exits for previous dtg
   if [ ! -f bk_fld_list2.${dtgbk} ]; then 
      if [ $BACKGROUND3 != $BACKGROUND2 ] ; then
         set_field_list bk_fld_list2.${dtgbk} ${dtgbk} $BACKGROUND3
      fi
   fi

   if [ $CYCLE_NO = "1" ] ; then 
     if [ ! -f  restart_${dtg} ] ; then
cat << eof > restart_$dtg
 &restartnl
      iupd = 0,
 &end
eof
     fi
   fi

########################################################################
#  Start processes obs
########################################################################

   coamps_grid_nml=coamps_grid_${dtg}.nl
########################################################################
#  Reset the filename format based on what is actually available
########################################################################
   if [ $nogaps_ftype -lt 64 ]; then IDBMS_NGPS=4; else IDBMS_NGPS=2; fi
   ${SED} "s/^\(\s*IDBMS\s*=\s*\).*,\s*/\1${IDBMS_NGPS},/" ${coamps_grid_nml}

########################################################################
#  list boundary conditions
########################################################################

   ${NAVDAS_BIN}/bclist_coamps.exe < ${coamps_grid_nml} > ${NAVDAS_LOG}/bclist_coamps.${dtg}
   check_return $? bclist_coamps.exe

########################################################################
#  list available fields useable for background, then decide which to use   
########################################################################

   ${NAVDAS_BIN}/bklist_coamps.exe < ${coamps_grid_nml} > ${NAVDAS_LOG}/bklist_coamps.${dtg}
   check_return $? bklist_coamps.exe

########################################################################
#  generate grid 0 nogaps background fields  
########################################################################

   ${NAVDAS_BIN}/bkinit_coamps.exe < ${coamps_grid_nml} > ${NAVDAS_LOG}/bkinit_coamps.${dtg}
   check_return $? bkinit_coamps.exe

########################################################################
#  process conventional observations   
########################################################################

   echo 'unclassified' | ${NAVDAS_BIN}/innov_coamps.exe >  ${NAVDAS_LOG}/innov_coamps_${NESTS_AN}.$dtg
   check_return $? innov_coamps.exe

########################################################################
#  process satellite winds  
########################################################################


   echo 'wvwnd' | ${NAVDAS_BIN}/satwnd_coamps_new.exe > ${NAVDAS_LOG}/satwnd_coamps_wvwnd_${NESTS_AN}.${dtg}
   check_return $? "satwnd_coamps_new.exe wvwnd"

   echo 'scatt' | ${NAVDAS_BIN}/satwnd_coamps_new.exe > ${NAVDAS_LOG}/satwnd_coamps_scatt_${NESTS_AN}.${dtg}
   check_return $? "satwnd_coamps_new.exe scatt"

   echo 'Qscat' | ${NAVDAS_BIN}/satwnd_coamps_new.exe > ${NAVDAS_LOG}/satwnd_coamps_Qscat_${NESTS_AN}.${dtg}
   check_return $? "satwnd_coamps_new.exe Qscat"

   echo 'Ascat' | ${NAVDAS_BIN}/satwnd_coamps_new.exe > ${NAVDAS_LOG}/satwnd_coamps_Ascat_${NESTS_AN}.${dtg}
   check_return $? "satwnd_coamps_new.exe Ascat"

   echo 'wnsat' | ${NAVDAS_BIN}/satwnd_coamps_new.exe > ${NAVDAS_LOG}/satwnd_coamps_wnsat_${NESTS_AN}.${dtg}
   check_return $? "satwnd_coamps_new.exe wnsat"

   echo 'wnspw' | ${NAVDAS_BIN}/satwnd_coamps_new.exe > ${NAVDAS_LOG}/satwnd_coamps_wnspw_${NESTS_AN}.${dtg}
   check_return $? "satwnd_coamps_new.exe wnspw"

   echo 'ssmiv' | ${NAVDAS_BIN}/satwnd_coamps_new.exe > ${NAVDAS_LOG}/satwnd_coamps_ssmiv_${NESTS_AN}.${dtg}
   check_return $? "satwnd_coamps_new.exe ssmiv"

   echo 'modis' | ${NAVDAS_BIN}/satwnd_coamps_new.exe > ${NAVDAS_LOG}/satwnd_coamps_modis_${NESTS_AN}.${dtg}
   check_return $? "satwnd_coamps_new.exe modis"

   echo 'avhrr' | ${NAVDAS_BIN}/satwnd_coamps_new.exe > ${NAVDAS_LOG}/satwnd_coamps_avhrr_${NESTS_AN}.${dtg}
   check_return $? "satwnd_coamps_new.exe avhrr"
 
########################################################################
#  process tovs retrievals   
########################################################################

   ${NAVDAS_BIN}/t_atovs_coamps.exe > ${NAVDAS_LOG}/t_atovs_coamps_${NESTS_AN}.$dtg
   check_return $? t_atovs_coamps.exe

########################################################################
#  process tovs radiances   
########################################################################

#   ${NAVDAS_BIN}/coamps_rttovs_old.exe < ${coamps_grid_nml} > ${NAVDAS_LOG}/coamps_rtovs_old_${NESTS_AN}.$dtg
#   check_return $? coamps_rttovs.exe

########################################################################
#  process tppw humidity retrievals   
########################################################################

   echo 'both ' | ${NAVDAS_BIN}/totppw_coamps.exe > ${NAVDAS_LOG}/totppw_coamps_${NESTS_AN}.$dtg
   check_return $? totppw_coamps.exe

########################################################################
#  sort the observations to be used in analysis   
########################################################################

   echo 'unclassified' | ${NAVDAS_BIN}/sort_coamps.exe > ${NAVDAS_LOG}/sort_coamps_${NESTS_AN}.$dtg
   check_return $? sort_coamps.exe

#######################################################################
#   normal exit
#######################################################################

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

