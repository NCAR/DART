#! /bin/csh -f
#
# This code may (or may not) be part of the CESM distribution,
# So it is not protected by the DART copyright agreement.
#
# DART $Id$

set exedir = $RUNDIR; cd $exedir

#------------------------------------------------------------------------------
# specify input data files
#------------------------------------------------------------------------------
# If the user changes any input datasets - be sure they have unique filenames.
# Do not duplicate existing input file names.
# Note that streams namelist input has the form
#      streams = 'stream1.txt year_align year_first year_last ',
#                'stream2.txt year_align year_first year_last ',
#                ...
#                'streamN.txt year_align year_first year_last '
# where
# streamN.txt is the stream description file containing input stream details
# year_first  is the first year of data that will be used
# year_last   is the last  year of data that will be used
# year_align  is the model year that will be aligned with data for year_first
#------------------------------------------------------------------------------

set DOMAINPATH = $DIN_LOC_ROOT/atm/datm7/domain.clm
set DOMAINFILE = $DOMAINPATH/domain.lnd.fv1.9x2.5_USGS.090106.nc

echo DOMAINFILE = $DOMAINFILE >! $CASEBUILD/datm.input_data_list

set atm_inst_counter = 1

while ($atm_inst_counter <= $NINST_ATM)
    if ($NINST_ATM > 1) then
       set atm_inst_string = `printf _%04d $atm_inst_counter`
    else
       set atm_inst_string = ''
    endif
    set atm_in_filename = datm_atm_in$atm_inst_string

set FFN        = "unused "

cat >! $atm_in_filename << EOF1
 &shr_strdata_nml
   dataMode       = 'CPLHIST'
   domainFile     = '$DOMAINFILE'
   streams        = 'CPLHIST.3Hrly.f09in.stream$atm_inst_string.Solar.txt  2000 2000 2000 ',
                    'CPLHIST.3Hrly.f09in.stream$atm_inst_string.Precip.txt 2000 2000 2000 ',
                    'CPLHIST.3Hrly.f09in.stream$atm_inst_string.Other.txt  2000 2000 2000 ',
                    'presaero.stream.txt 1 1 1'
   vectors        = 'u:v'
   mapmask        = 'nomask',
                    'nomask',
                    'nomask',
		    'nomask'
   tintalgo       = 'coszen',
                    'nearest',
                    'linear',
		    'linear'
  /
EOF1

cat >! CPLHIST.3Hrly.f09in.stream$atm_inst_string.Solar.txt << EOF1
<streamstemplate>
      <general_comment>
         streams template for datm in CCSM4
      </general_comment>
<stream>
      <comment>
         Stream description file for CPL history 3-hourly Solar data at 0.9x1.25 resolution for BGC spinup
      </comment>
      <dataSource>
         CPL
      </dataSource>
      <domainInfo>
         <variableNames>
            time    time
            doma_lon      lon
            doma_lat      lat
            doma_area    area
            doma_mask    mask
         </variableNames>
         <filePath>
            /glade/proj3/DART/CAM_DATM/2000
         </filePath>
         <fileNames>
            CAM_DATM.cpl$atm_inst_string.ha2x1dx6h.2000.nc
         </fileNames>
      </domainInfo>
      <fieldInfo>
         <variableNames>
            a2x6h_Faxa_swndr     swndr
            a2x6h_Faxa_swvdr     swvdr
            a2x6h_Faxa_swndf     swndf
            a2x6h_Faxa_swvdf     swvdf
         </variableNames>
         <filePath>
            /glade/proj3/DART/CAM_DATM/2000
         </filePath>
         <tInterpAlgo>
            coszen
         </tInterpAlgo>
         <offset>
            -21600
         </offset>
         <fileNames>
            CAM_DATM.cpl$atm_inst_string.ha2x1dx6h.2000.nc
         </fileNames>
      </fieldInfo>
      <!-- Information on the program that created this file -->
      <build_streams_documentation>
         This CCSM stream text file was created by build_streams using the command line:
               /blhome/afox/Fei/cesm1_1_beta01/scripts/ccsm_utils/Tools/build_streams -t datm.template.streams.xml -p /glade/data01/CMIP5/CCSM/csm/%c/cpl/hist -c enstest_0727 -b 2000 -e 2000 -s CPLHIST3Hrly.Solar
         For more information on build_streams:
             /blhome/yfzhang/cesm1_1_beta01/scripts/ccsm_utils/Tools/build_streams -help
      </build_streams_documentation>
</stream>
</streamstemplate>
EOF1

cat >! CPLHIST.3Hrly.f09in.stream$atm_inst_string.Precip.txt << EOF1
<streamstemplate>
      <general_comment>
         streams template for datm in CCSM4
      </general_comment>
<stream>
      <comment>
         Stream description file for CPL history 3-hourly Precip data at 0.9x1.25 resolution for BGC spinup
      </comment>
      <dataSource>
         CPL
      </dataSource>
      <domainInfo>
         <variableNames>
            time        time
            doma_lon          lon
            doma_lat          lat
            doma_area        area
            doma_mask        mask
         </variableNames>
         <filePath>
            /glade/proj3/DART/CAM_DATM/2000
         </filePath>
         <fileNames>
            CAM_DATM.cpl$atm_inst_string.ha2x1dx6h.2000.nc
         </fileNames>
      </domainInfo>
      <fieldInfo>
         <variableNames>
            a2x6h_Faxa_rainc     rainc
            a2x6h_Faxa_rainl     rainl
            a2x6h_Faxa_snowc     snowc
            a2x6h_Faxa_snowl     snowl
         </variableNames>
         <filePath>
            /glade/proj3/DART/CAM_DATM/2000
         </filePath>
         <tInterpAlgo>
            nearest
         </tInterpAlgo>
         <offset>
            -10800
         </offset>
         <fileNames>
            CAM_DATM.cpl$atm_inst_string.ha2x1dx6h.2000.nc
         </fileNames>
      </fieldInfo>
      <!-- Information on the program that created this file -->
      <build_streams_documentation>
         This CCSM stream text file was created by build_streams using the command line:
               /blhome/yfzhang/cesm1_1_beta01/scripts/ccsm_utils/Tools/build_streams -t datm.template.streams.xml -p /glade/data01/CMIP5/CCSM/csm/%c/cpl/hist -c enstest_0727 -b 2000 -e 2000 -s CPLHIST3Hrly.Precip
         For more information on build_streams:
             /blhome/yfzhang/cesm1_1_beta01/scripts/ccsm_utils/Tools/build_streams -help
      </build_streams_documentation>
</stream>
</streamstemplate>
EOF1

cat >! CPLHIST.3Hrly.f09in.stream$atm_inst_string.Other.txt << EOF1
<streamstemplate>
      <general_comment>
         streams template for datm in CCSM4
      </general_comment>
<stream>
      <comment>
         Stream description file for CPL history 3-hourly non-Precipitation and non-Solar data at 0.9x1.25 resolution for BGC spinup
      </comment>
      <dataSource>
         CPL
      </dataSource>
      <domainInfo>
         <variableNames>
            time    time
            doma_lon      lon
            doma_lat      lat
            doma_area    area
            doma_mask    mask
         </variableNames>
         <filePath>
            /glade/proj3/DART/CAM_DATM/2000
         </filePath>
         <fileNames>
            CAM_DATM.cpl$atm_inst_string.ha2x1dx6h.2000.nc
         </fileNames>
      </domainInfo>
      <fieldInfo>
         <variableNames>
            a2x6h_Sa_z           z
            a2x6h_Sa_u           u
            a2x6h_Sa_v           v
            a2x6h_Sa_tbot        tbot
            a2x6h_Sa_ptem        ptem
            a2x6h_Sa_shum        shum
            a2x6h_Sa_pbot        pbot
            a2x6h_Faxa_lwdn      lwdn
            a2x6h_Sa_dens        dens
            a2x6h_Sa_pslv        pslv
         </variableNames>
         <filePath>
            /glade/proj3/DART/CAM_DATM/2000
         </filePath>
         <tInterpAlgo>
            linear
         </tInterpAlgo>
         <offset>
            -10800
         </offset>
         <fileNames>
            CAM_DATM.cpl$atm_inst_string.ha2x1dx6h.2000.nc
         </fileNames>
      </fieldInfo>
      <!-- Information on the program that created this file -->
      <build_streams_documentation>
         This CCSM stream text file was created by build_streams using the command line:
               /blhome/yfzhang/cesm1_1_beta01/scripts/ccsm_utils/Tools/build_streams -t datm.template.streams.xml -p /glade/data01/CMIP5/CCSM/csm/%c/cpl/hist -c enstest_0727 -b 2000 -e 2000 -s CPLHIST3Hrly.nonSolarNonPrecip
         For more information on build_streams:
             /blhome/yfzhang/cesm1_1_beta01/scripts/ccsm_utils/Tools/build_streams -help
      </build_streams_documentation>
</stream>
</streamstemplate>
EOF1


#$CASETOOLS/listfilesin_streams -input_data_list -t CPLHIST.3Hrly.f09in.stream.Solar.txt >> $CASEBUILD/datm.input_data_list
#$CASETOOLS/listfilesin_streams -input_data_list -t CPLHIST.3Hrly.f09in.stream.Precip.txt >> $CASEBUILD/datm.input_data_list
#$CASETOOLS/listfilesin_streams -input_data_list -t CPLHIST.3Hrly.f09in.stream.Other.txt >> $CASEBUILD/datm.input_data_list

 @ atm_inst_counter = $atm_inst_counter + 1
end

# ---
# End of the loop over the atmospheric instance counter
# ---

cat >! presaero.stream.txt << EOF1
<streamstemplate>
      <general_comment>
         streams template for datm in CCSM4
      </general_comment>
<stream>
      <comment>
         Stream description file for aerosol deposition
      </comment>
      <dataSource>
         presaero
      </dataSource>
      <domainInfo>
         <variableNames>
            time    time
            lon      lon
            lat      lat
            area    area
            mask    mask
         </variableNames>
         <filePath>
            $DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/aero
         </filePath>
         <fileNames>
            aerosoldep_monthly_2000_mean_1.9x2.5_c090421.nc
         </fileNames>
      </domainInfo>
      <fieldInfo>
         <variableNames>
            BCDEPWET   bcphiwet
            BCPHODRY   bcphodry
            BCPHIDRY   bcphidry
            OCDEPWET   ocphiwet
            OCPHIDRY   ocphidry
            OCPHODRY   ocphodry
            DSTX01WD   dstwet1
            DSTX01DD   dstdry1
            DSTX02WD   dstwet2
            DSTX02DD   dstdry2
            DSTX03WD   dstwet3
            DSTX03DD   dstdry3
            DSTX04WD   dstwet4
            DSTX04DD   dstdry4
         </variableNames>
         <filePath>
            $DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/aero
         </filePath>
         <offset>
            0
         </offset>
         <fileNames>
            aerosoldep_monthly_2000_mean_1.9x2.5_c090421.nc
         </fileNames>
      </fieldInfo>
      <!-- Information on the program that created this file -->
      <build_streams_documentation>
         This CCSM stream text file was created by build_streams using the command line:
               /blhome/yfzhang/cesm1_1_beta01/scripts/ccsm_utils/Tools/build_streams -t datm.template.streams.xml -s presaero -b 1 -e 1 -p $DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/aero -c aerosoldep_monthly_2000_mean_1.9x2.5_c090421.nc -dp $DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/aero -do aerosoldep_monthly_2000_mean_1.9x2.5_c090421.nc
         For more information on build_streams:
             /blhome/yfzhang/cesm1_1_beta01/scripts/ccsm_utils/Tools/build_streams -help
      </build_streams_documentation>
</stream>
</streamstemplate>
EOF1


$CASETOOLS/listfilesin_streams -input_data_list -t presaero.stream.txt >> $CASEBUILD/datm.input_data_list


set base_filename = "datm_in"
set atm_inst_counter = 1
while ($atm_inst_counter <= $NINST_ATM)
    if ($NINST_ATM > 1) then
       set atm_inst_string = `printf _%04d $atm_inst_counter`
    else
       set atm_inst_string = ""
    endif
    set in_filename = ${base_filename}${atm_inst_string}

cat >! ${in_filename} << EOF1
  &datm_nml
    atm_in = 'datm_atm_in$atm_inst_string'
    decomp = '1d'
    factorFn = '$FFN'
    iradsw   = 1
    presaero = .true.
  /
EOF1

    @ atm_inst_counter = $atm_inst_counter + 1
end

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

