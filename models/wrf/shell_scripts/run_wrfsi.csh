#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

# Purpose: Run wrfsi V2.1
#
# set echo

#Install wrfsi.

 set    START_DATE      = $1
 set    FORECAST_LENGTH = $2
 set    DATA_INTERVAL   = $3
 setenv	DATAROOT          $4
 setenv INSTALLROOT       $5
 setenv WEST_EAST_GRIDS	  $6
 setenv SOUTH_NORTH_GRIDS $7
 setenv VERTICAL_GRIDS	  $8
 setenv GRID_DISTANCE	  $9
 setenv MY_MOAD_KNOWN_LAT $10
 setenv MY_MOAD_KNOWN_LON $11
 setenv MY_MOAD_STAND_LONS $12
 setenv MY_NUM_DOMAINS    $13
 setenv LLI2              $14
 setenv LLJ2              $15
 setenv URI2              $16
 setenv URJ2              $17
 setenv DATASOURCE        $18
 setenv CASE_NAME         $19
 setenv GRIB_DATA         $20
 setenv RUN_GRID_GEN      $21
 setenv GRIB_DATA_DIR     $22

 set MY_NUM_ACTIVE_SUBNESTS = `expr ${MY_NUM_DOMAINS} \- 1`

 setenv SOURCE_ROOT   $INSTALLROOT
 setenv GEOG_DATAROOT /usr/local/wrfsi/SI_GEOG
 setenv GRIB_DATADIR  ${DATAROOT}/${GRIB_DATA}

 set OS = `uname`

 if((${OS} == 'AIX') || (${OS} == 'Linux')) then
    setenv PATH_TO_PERL /usr/bin
 else
    setenv PATH_TO_PERL /usr/local/bin
 endif

# setenv NETCDF	/usr/local/netcdf
 setenv NETCDF	/ocotillo/users/caya/netcdf

 set BGN_CCYYMMDDHH = ${START_DATE}
 set END_CCYYMMDDHH = `advance_cymdh ${START_DATE} ${FORECAST_LENGTH}`

 set INTERVAL_IN_HOUR = ${DATA_INTERVAL}

 setenv BGN_CCYYMMDDHH ${BGN_CCYYMMDDHH}
 setenv END_CCYYMMDDHH ${END_CCYYMMDDHH}
 
 setenv	MY_PROJECTION		lambert
 setenv MY_MOAD_STAND_LATS	"30.0, 60.0"
 set    VERTIVAL_LEVELS	=	( "1.000, 0.990, 0.978, 0.964, 0.946," \
                                  "0.922, 0.894, 0.860, 0.817, 0.766," \
                                  "0.707, 0.644, 0.576, 0.507, 0.444," \
                                  "0.380, 0.324, 0.273, 0.228, 0.188," \
                                  "0.152, 0.121, 0.093, 0.069, 0.048," \
                                  "0.029, 0.014, 0.000" )

#--------------------------------------------------------------------

 setenv	INSTALLROOT	${SOURCE_ROOT}
 setenv	MOAD_DATAROOT	${DATAROOT}/data
 setenv	EXT_DATAROOT	${DATAROOT}/data
 setenv TEMPLATES       ${SOURCE_ROOT}/templates/${CASE_NAME}
 setenv MPICH		/usr/local/mpich

#--------------------------------------------------------------------
 set user_desc = $user

 set date_string_ccyymmddhh = ${BGN_CCYYMMDDHH}
 set s_year  = `echo ${date_string_ccyymmddhh} | cut -c1-4`
 set s_month = `echo ${date_string_ccyymmddhh} | cut -c5-6`
 set s_day   = `echo ${date_string_ccyymmddhh} | cut -c7-8`
 set s_hour  = `echo ${date_string_ccyymmddhh} | cut -c9-10`

 set date_string_ccyymmddhh = ${END_CCYYMMDDHH}
 set e_year  = `echo ${date_string_ccyymmddhh} | cut -c1-4`
 set e_month = `echo ${date_string_ccyymmddhh} | cut -c5-6`
 set e_day   = `echo ${date_string_ccyymmddhh} | cut -c7-8`
 set e_hour  = `echo ${date_string_ccyymmddhh} | cut -c9-10`

 @ INTERVAL_IN_SECOND = 3600 * ${INTERVAL_IN_HOUR}
 
#--------------------------------------------------------------------

install_wrfsi:

 set need_install = no 

#If any of the following executable is missing, we need to install.

 foreach f ( \
	grib_prep.exe		\
	gridgen_model.exe	\
	hinterp.exe		\
	siscan			\
	vinterp.exe		\
	)

    if(! -f $INSTALLROOT/bin/$f ) then
       set need_install = yes
    endif
 end

 if ( ${need_install} == "yes" ) then
    set cwd = `pwd`
    cd ${INSTALLROOT}
    ${PATH_TO_PERL}/perl install_wrfsi.pl
    cd ${cwd}
 endif

#--------------------------------------------------------------------

#get_avn_data:

# get_avn.csh
 
#--------------------------------------------------------------------

 mkdir -p ${DATAROOT}/templates
 cd ${DATAROOT}/templates

# mkdir -p ${TEMPLATES}
# cd ${TEMPLATES}

 rm -f script.sed

cat > script.sed << EOF
 /SIMULATION_NAME/c\
 SIMULATION_NAME = '${CASE_NAME}',
 /USER_DESC/c\
 USER_DESC = '${user}'
 /START_YEAR/c\
 START_YEAR = ${s_year},
 /START_MONTH/c\
 START_MONTH = ${s_month},
 /START_DAY/c\
 START_DAY = ${s_day},
 /START_HOUR/c\
 START_HOUR = ${s_hour},
 /START_MINUTE/c\
 START_MINUTE = 00,
 /START_SECOND/c\
 START_SECOND = 00,
 /END_YEAR/c\
 END_YEAR = ${e_year},
 /END_MONTH/c\
 END_MONTH = ${e_month},
 /END_DAY/c\
 END_DAY = ${e_day},
 /END_HOUR/c\
 END_HOUR = ${e_hour},
 /END_MINUTE/c\
 END_MINUTE = 00,
 /END_SECOND/c\
 END_SECOND = 00,
 /INTERVAL/c\
 INTERVAL = ${INTERVAL_IN_SECOND}
 /NUM_DOMAINS/c\
 NUM_DOMAINS = ${MY_NUM_DOMAINS},
 /XDIM/c\
 XDIM = ${WEST_EAST_GRIDS},
 /YDIM/c\
 YDIM = ${SOUTH_NORTH_GRIDS},
 /DOMAIN_ORIGIN_LLI/c\
 DOMAIN_ORIGIN_LLI = 1,${LLI2},
 /DOMAIN_ORIGIN_LLJ/c\
 DOMAIN_ORIGIN_LLJ = 1,${LLJ2},
 /DOMAIN_ORIGIN_URI/c\
 DOMAIN_ORIGIN_URI = ${WEST_EAST_GRIDS},${URI2},
 /DOMAIN_ORIGIN_URJ/c\
 DOMAIN_ORIGIN_URJ = ${SOUTH_NORTH_GRIDS},${URJ2},
 /MAP_PROJ_NAME/c\
 MAP_PROJ_NAME = '${MY_PROJECTION}',
 /MOAD_KNOWN_LAT/c\
 MOAD_KNOWN_LAT = ${MY_MOAD_KNOWN_LAT},
 /MOAD_KNOWN_LON/c\
 MOAD_KNOWN_LON = ${MY_MOAD_KNOWN_LON},
 /MOAD_STAND_LATS/c\
 MOAD_STAND_LATS = ${MY_MOAD_STAND_LATS},
 /MOAD_STAND_LONS/c\
 MOAD_STAND_LONS = ${MY_MOAD_STAND_LONS},
 /MOAD_DELTA_X/c\
 MOAD_DELTA_X = ${GRID_DISTANCE},
 /MOAD_DELTA_Y/c\
 MOAD_DELTA_Y = ${GRID_DISTANCE},
 /TOPO_30S/c\
 TOPO_30S = '${GEOG_DATAROOT}/topo_30s',
 /LANDUSE_30S/c\
 LANDUSE_30S = '${GEOG_DATAROOT}/landuse_30s',
 /SOILTYPE_TOP_30S/c\
 SOILTYPE_TOP_30S = '${GEOG_DATAROOT}/soiltype_top_30s',
 /SOILTYPE_BOT_30S/c\
 SOILTYPE_BOT_30S = '${GEOG_DATAROOT}/soiltype_bot_30s',
 /GREENFRAC/c\
 GREENFRAC = '${GEOG_DATAROOT}/greenfrac',
 /SOILTEMP_1DEG/c\
 SOILTEMP_1DEG = '${GEOG_DATAROOT}/soiltemp_1deg',
 /ALBEDO_NCEP/c\
 ALBEDO_NCEP = '${GEOG_DATAROOT}/albedo_ncep',
 /MAXSNOWALB/c\
 MAXSNOWALB = '${GEOG_DATAROOT}/maxsnowalb',
 /ISLOPE/c\
 ISLOPE = '${GEOG_DATAROOT}/islope'
 /NUM_ACTIVE_SUBNESTS/c\
 NUM_ACTIVE_SUBNESTS = ${MY_NUM_ACTIVE_SUBNESTS},
 /INIT_ROOT/c\
 INIT_ROOT = '${DATASOURCE}',
 /LBC_ROOT/c\
 LBC_ROOT = '${DATASOURCE}',
 /LEVELS/c\
 LEVELS = ${VERTIVAL_LEVELS}
 /ANALPATH/c\
 ANALPATH = '${EXT_DATAROOT}/extprd',
 /LBCPATH/c\
 LBCPATH = '${EXT_DATAROOT}/extprd',
 /LSMPATH/c\
 LSMPATH = '${EXT_DATAROOT}/extprd',
 /CONSTANTS_PATH/c\
 CONSTANTS_PATH = '${EXT_DATAROOT}/extprd'
EOF

 sed -f script.sed \
    ${SOURCE_ROOT}/templates/default/wrfsi.nl > wrfsi.nl

#--------------------------------------------------------------------

grid_gen:

 if ( $RUN_GRID_GEN ) then
    mkdir -p ${TEMPLATES}
    cp -pv ${DATAROOT}/templates/wrfsi.nl ${TEMPLATES}/wrfsi.nl

    cd ${INSTALLROOT}/etc
    ./window_domain_rt.pl -w wrfsi -t ${TEMPLATES} -c \
           -s ${SOURCE_ROOT}

    cp -r ${MOAD_DATAROOT} ${INSTALLROOT}/domains/${CASE_NAME}
 else
    cp -rf ${INSTALLROOT}/domains/${CASE_NAME} ${MOAD_DATAROOT}
 endif

# THE SCRIPT window_domain_rt.pl GENERATE THE DIRECTORY STRUCTURE OF
# $MOAD_DATAROOT. SINCE WE CHOSE $EXT_DATAROOT == $MOAD_DATAROOT,
# ${EXT_DATAROOT}/extprd HAS TO BE CREATED AFTER window_domain_rt.pl IS RUN.

 mkdir -p ${EXT_DATAROOT}/extprd

#--------------------------------------------------------------------

grip_prep:

# link to appropriate GRIB DATA
 set END_DATE = `advance_cymdh $START_DATE $FORECAST_LENGTH`  # End time of forecast
 mkdir -p $GRIB_DATADIR
 cd $GRIB_DATADIR

 set DATE = $START_DATE

 while ( $DATE <= $END_DATE )

     set MM = `echo $DATE | cut -c5-6`
     set DD = `echo $DATE | cut -c7-8`
     set HH = `echo $DATE | cut -c9-10`

   if ($GRIB_DATA == AVN) then
     set YY = `echo $DATE | cut -c3-4`
     set GRIB_FILE = fnl_${YY}${MM}${DD}_${HH}_00
   else
     set YY = `echo $DATE | cut -c1-4`
     set GRIB_FILE = ${YY}${MM}${DD}${HH}.AWIP3D00.tm00
   endif

   if ( -e $GRIB_DATA_DIR/$GRIB_FILE ) then
     ln -s $GRIB_DATA_DIR/$GRIB_FILE $GRIB_FILE
   else
     echo FAILED TO FIND $GRIB_DATA_DIR/$GRIB_FILE
     exit 1
   endif

   set DATE=`advance_cymdh ${DATE} ${INTERVAL_IN_HOUR}`

   end

 cd ${MOAD_DATAROOT}/static

 cp -pv $SOURCE_ROOT/extdata/static/Vtable.$DATASOURCE .

cat >! ./grib_prep.nl << EOF
&filetimespec
 START_YEAR   = ${s_year},
 START_MONTH  = ${s_month},
 START_DAY    = ${s_day},
 START_HOUR   = ${s_hour},
 START_MINUTE = 00,
 START_SECOND = 00,
 END_YEAR     = ${e_year},
 END_MONTH    = ${e_month},
 END_DAY      = ${e_day},
 END_HOUR     = ${e_hour},
 END_MINUTE   = 00,
 END_SECOND   = 00,
 INTERVAL     = ${INTERVAL_IN_SECOND}
/
&gpinput_defs
 SRCNAME = '$DATASOURCE', 'NNRP', 'NNRPSFC', 'AWIP', 'RUCH', 'NNRP', 'NNRPSFC', 'SST'
 SRCVTAB = '$DATASOURCE', 'NNRP', 'NNRPSFC', 'AWIP', 'RUCH', 'NNRPSFC', 'NNRPSFC', 'SST'
 SRCPATH = '${GRIB_DATADIR}',
           '${GRIB_DATADIR}/sfc',
           '${GRIB_DATADIR}',
           '/path/to/RUCH',
           '/path/to/nnrp/grib',
           '/path/to/nnrp/sfc/grib',
           '/public/data/grids/ncep/sst/grib'
 SRCCYCLE = 6, 6, 6, 6, 12, 12, 24
 SRCDELAY = 3, 4, 4, 3, 0, 0, 36
/
EOF

 cd ${MOAD_DATAROOT}

wrfsi:

 $INSTALLROOT/etc/wrfsi.pl -d ${MOAD_DATAROOT} \
	${BGN_CCYYMMDDHH} ${FORECAST_LENGTH} ${DATASOURCE}

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

