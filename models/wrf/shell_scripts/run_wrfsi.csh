#!/bin/csh -f
#----------------------------------------------------------------------
# Purpose: Run wrfsi
#-----------------------------------------------------------------------
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

### setenv	MACHINE		"ibmncar"

 setenv SOURCE_ROOT   $INSTALLROOT
 setenv GEOG_DATAROOT /usr/local/wrfsi/SI_GEOG
 setenv GRIB_DATADIR  ${DATAROOT}/${GRIB_DATA}

 set OS = `uname`

 if((${OS} == 'AIX') || (${OS} == 'Linux')) then
    setenv PATH_TO_PERL /usr/bin
 else
    setenv PATH_TO_PERL /usr/local/bin
 endif

 if(${OS} == 'AIX') then
    if ( ! -d /home/blackforest/$user/bin/netcdf_links ) then
     mkdir -p /home/blackforest/$user/bin/netcdf_links
     set cwd = `pwd`
     cd /home/blackforest/$user/bin/netcdf_links
     ln -sf /usr/local/apps/netcdf-3.5.1/bin bin
     ln -sf /usr/local/apps/netcdf-3.5.1/man man
     ln -sf /usr/local/apps/netcdf-3.5.1/include include
     ln -sf /usr/local/apps/netcdf-3.5.1/lib32/r4i4 lib
     cd ${cwd}
    endif

    setenv NETCDF	/home/blackforest/$user/bin/netcdf_links
# else
#    setenv NETCDF	/usr/local/netcdf
    setenv NETCDF	/ocotillo/users/caya/netcdf
 endif

 set BGN_CCYYMMDDHH = ${START_DATE}
 set END_CCYYMMDDHH = `advance_cymdh ${START_DATE} ${FORECAST_LENGTH}`

 set INTERVAL_IN_HOUR = ${DATA_INTERVAL}

 setenv BGN_CCYYMMDDHH ${BGN_CCYYMMDDHH}
 setenv END_CCYYMMDDHH ${END_CCYYMMDDHH}
 setenv INTERVAL_IN_HOUR ${INTERVAL_IN_HOUR}

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

 if ( ! -d ${DATAROOT} ) mkdir -p ${DATAROOT}
 cd ${DATAROOT}

 if(! -d templates ) mkdir -p templates
 cd templates

# if(! -d ${TEMPLATES} ) mkdir -p ${TEMPLATES}
# cd ${TEMPLATES}

 rm -f script.sed

cat > script.sed << EOF
 s/MY_SIMULATION_NAME/${CASE_NAME}/g
 s/MY_USER_DESC/${user}/g
 s/MY_START_YEAR/${s_year}/g
 s/MY_START_MONTH/${s_month}/g
 s/MY_START_DAY/${s_day}/g
 s/MY_START_HOUR/${s_hour}/g
 s/MY_START_MINUTE/00/g
 s/MY_START_SECOND/00/g
 s/MY_END_YEAR/${e_year}/g
 s/MY_END_MONTH/${e_month}/g
 s/MY_END_DAY/${e_day}/g
 s/MY_END_HOUR/${e_hour}/g
 s/MY_END_MINUTE/00/g
 s/MY_END_SECOND/00/g
 s/MY_INTERVAL/${INTERVAL_IN_SECOND}/g
 s/MY_NUM_DOMAINS/${MY_NUM_DOMAINS}/g
 s/MY_PROJECTION/${MY_PROJECTION}/g
 s/LLI2/${LLI2}/g
 s/LLJ2/${LLJ2}/g
 s/URI2/${URI2}/g
 s/URJ2/${URJ2}/g
 s/MY_INIT_DATA/${DATASOURCE}/g
 s/MY_MOAD_KNOWN_LAT/${MY_MOAD_KNOWN_LAT}/g
 s/MY_MOAD_KNOWN_LON/${MY_MOAD_KNOWN_LON}/g
 s/MY_MOAD_STAND_LATS/${MY_MOAD_STAND_LATS}/g
 s/MY_MOAD_STAND_LONS/${MY_MOAD_STAND_LONS}/g
 s/GRID_DISTANCE/${GRID_DISTANCE}/g
 s/WEST_EAST_GRIDS/${WEST_EAST_GRIDS}/g
 s/SOUTH_NORTH_GRIDS/${SOUTH_NORTH_GRIDS}/g
 s/VERTIVAL_LEVELS/${VERTIVAL_LEVELS}/g
 s/MY_NUM_ACTIVE_SUBNESTS/${MY_NUM_ACTIVE_SUBNESTS}/g
EOF

 sed -f script.sed \
    ${SOURCE_ROOT}/templates/default/wrfsi.nl.template > sample

 m4 -DMY_GEOG_ROOT=${GEOG_DATAROOT} \
    -DMY_EXT_DATAROOT=${EXT_DATAROOT} \
    sample > wrfsi.nl

 cp wrfsi.nl ${DATAROOT}/wrfsi.nl

#--------------------------------------------------------------------

grid_gen:

 if(! -d ${DATAROOT}/data/cdl ) mkdir -p ${MOAD_DATAROOT}/cdl
 if(! -d ${DATAROOT}/static ) mkdir -p ${MOAD_DATAROOT}/static
 if ( ! -d ${MOAD_DATAROOT}/extprd ) mkdir -p ${MOAD_DATAROOT}/extprd
 if ( ! -d ${MOAD_DATAROOT}/log    ) mkdir -p ${MOAD_DATAROOT}/log
 if ( ! -d ${MOAD_DATAROOT}/siprd  ) mkdir -p ${MOAD_DATAROOT}/siprd

 if ( $RUN_GRID_GEN ) then
    cp ${DATAROOT}/wrfsi.nl ${MOAD_DATAROOT}/cdl/wrfsi.cdl
    cp ${DATAROOT}/wrfsi.nl ${MOAD_DATAROOT}/static/.

    $INSTALLROOT/etc/window_domain_rt.pl -w wrfsi \
           -t ${TEMPLATES} -c \
           -s ${SOURCE_ROOT}

    $SOURCE_ROOT/bin/gridgen_model.exe ${MOAD_DATAROOT}
    cp ${MOAD_DATAROOT}/static/static.wrfsi.d01 \
       ${MOAD_DATAROOT}/static/static.wrfsi
    cp -r ${MOAD_DATAROOT}/static ${INSTALLROOT}/domains/${CASE_NAME}
 else
    cp -rf ${INSTALLROOT}/domains/${CASE_NAME}/static ${MOAD_DATAROOT}
    cp -f ${DATAROOT}/wrfsi.nl ${MOAD_DATAROOT}/cdl/wrfsi.cdl
    cp -f ${DATAROOT}/wrfsi.nl ${MOAD_DATAROOT}/static/.
 endif
 
#--------------------------------------------------------------------

grip_prep:

# link to appropriate GRIB DATA
 set END_DATE = `advance_cymdh $START_DATE $FORECAST_LENGTH`  # End time of forecast
 if ( ! -e $GRIB_DATADIR ) mkdir $GRIB_DATADIR
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

 cp $SOURCE_ROOT/extdata/static/Vtable.$DATASOURCE .

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

 $INSTALLROOT/etc/grib_prep.pl -l ${FORECAST_LENGTH} \
                               -s ${START_DATE} \
                               -t ${DATA_INTERVAL} \
                                  ${DATASOURCE}

 cd ${MOAD_DATAROOT}

 setenv EXT_DATAROOT ${MOAD_DATAROOT}

 $INSTALLROOT/etc/wrfprep.pl -f ${FORECAST_LENGTH} -s ${START_DATE}

wrfsi:

 $INSTALLROOT/etc/wrfsi.pl -d ${MOAD_DATAROOT} \
	${BGN_CCYYMMDDHH} ${FORECAST_LENGTH} ${DATASOURCE}

 exit ( 0 )

#mv ${EXT_DATAROOT}/siprd/wrf_real_input_em* /mmm/users/mslee/SIOUT_100KM/.
