#!/bin/sh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

#short-term archive script - move model output out of run directory
#to free disc space for next invocation of model run
#must be executed from run directory

# This must exist as ${CASEROOT}/Tools/st_archive.sh 
# ./xmlchange -file env_run.xml -id DOUT_S_SAVE_INT_REST_FILES -val 'TRUE'
#
# The original has been substantially modified to restage the CAM INITIAL files
# for successive coldstarts - as opposed to restage restart files. As of
# Nov 10, 2011 ... DART is using CAM initial files. This will change SOON.

#function dispose:
#moves output files to specified area of st archive and will
#process interim files along the way: 
#arg1 => interim files flag
#arg2 => destination
#remaining args => actual files to be processed 
dispose() {
    if [ "$1" == "ifiles_y" ] && [ "$DOUT_S_SAVE_INT_REST_FILES" != "TRUE" ]; then
	shift
	shift
	rm $*       2> /dev/null
    else
	shift
	dest=$1
	mkdir -p $dest
	shift
	mv $* $dest 2> /dev/null
    fi
}

#function get_inst_suffix:
#Gets a string corresponding to the current instance index of a given
#component, with a leading '.'; this can be appended to a file name.
#In the general case, the returned string is something like ".###",
#but if there is only a single instance of the given component, then
#the returned string is empty, because in that case the file names
#don't have any instance number (or the associated extra '.').
#
#arg1 => instance index for a given component
#arg2 => number of instances of this component
#
#result is returned in $inst_suffix
get_inst_suffix() {
#   echo "get_inst_suffixa $1 $2 ${inst_suffix}"
    if [ $2 -eq 1 ]; then   # only one instance of this component
	inst_suffix=""
    else                    # multiple instances of this component
        inst_suffix=`printf _%04d $1`
    fi
#   echo "get_inst_suffixb $1 $2 ${inst_suffix}"
}

echo ""
echo "st_archive.sh: start of short-term archiving"

#validate required env var settings
if [ -z "$DOUT_S_ROOT" ]; then
    echo "st_archive.sh: error, environment variable DOUT_S_ROOT is required "
    echo "               for root location of short-term archive"
    echo "st_archive.sh: exiting"
    exit 1
fi

sta=${DOUT_S_ROOT}/.sta-$$-`date +%Y%m%d%H%M%S%N`
mkdir -p ${sta} 2> /dev/null

if [ $? -ne 0 ]; then
    echo "st_archive.sh: error, unable to create short-term archive directory"
    echo "st_archive.sh: exiting"
    exit 1
fi
mv ${DOUT_S_ROOT}/* ${sta}

if [ -z "$DOUT_S_SAVE_INT_REST_FILES" ]; then
    echo "st_archive.sh: warning, environment variable DOUT_S_SAVE_INT_REST_FILES is not "
    echo "               set - using "FALSE" as default for saving interim restart files"
    export DOUT_S_SAVE_INT_REST_FILES=FALSE
fi

if [ "$DOUT_S_SAVE_INT_REST_FILES" == "FALSE" ]; then
    echo "st_archive.sh: restart files from end of run will be saved, "
    echo "               interim restart files will be deleted"
fi

## Check component instance counts

if [ -z "$NINST_ATM" ]; then
    echo "st_archive.sh: warning, NINST_ATM not set -- using 1 instance"
    export NINST_ATM=1
fi

if [ -z "$NINST_LND" ]; then
    echo "st_archive.sh: warning, NINST_LND not set -- using 1 instance"
    export NINST_LND=1
fi

if [ -z "$NINST_ICE" ]; then
    echo "st_archive.sh: warning, NINST_ICE not set -- using 1 instance"
    export NINST_ICE=1
fi

if [ -z "$NINST_OCN" ]; then
    echo "st_archive.sh: warning, NINST_OCN not set -- using 1 instance"
    export NINST_OCN=1
fi

if [ -z "$NINST_GLC" ]; then
    echo "st_archive.sh: warning, NINST_GLC not set -- using 1 instance"
    export NINST_GLC=1
fi

#create directory for restart files
set ${CASE}.cpl.r.*
cplfile=`ls -rt $* 2> /dev/null | tail -1`
dname=`echo $cplfile | sed "s/\.nc//; s/^.*\.r\.//;"`
if [ -d ${sta}/rest/${dname} ]; then
    rm -rf ${sta}/rest/${dname}
fi
mkdir -p ${sta}/rest/${dname}
if [ $? -ne 0 ]; then
    echo "st_archive.sh: error, unable to create rest directory"
    echo "st_archive.sh: exiting"
    exit 1
fi

#populate temp directory with pointer files
set rpointer.*
if [ $# -le 0 ]; then
    echo "st_archive.sh: error, script should be invoked from run directory..."
    echo "               expecting restart pointer files of the form 'rpointer.<component>'"
    echo "               but did not find any: $*"
    echo "st_archive.sh: exiting"
    exit 1
fi
# TJH mv $* ${sta}/rest/${dname}

set rpointer.*;                                                                                                       dispose ifiles_n ${sta}/cpl/logs $*
set cpl.log.*;                                                                                                        dispose ifiles_n ${sta}/cpl/logs $*
set ccsm*.log.*;                                                                                                      dispose ifiles_n ${sta}/cpl/logs $*
set ${CASE}.cpl.r.*;                                                                                                  dispose ifiles_y ${sta}/cpl/rest $*
set ${CASE}.cpl.h* ;                                                                                                  dispose ifiles_n ${sta}/cpl/hist $*


#possible tweaking - remove assimilate_dir.* directories? anything else?
set dart_log.*;                                                                                                       dispose ifiles_n ${sta}/dart/logs $*
set True_State.*.nc;                                                                                                  dispose ifiles_n ${sta}/dart/hist $*
set Prior_Diag.*.nc;                                                                                                  dispose ifiles_n ${sta}/dart/hist $*
set Posterior_Diag.*.nc;                                                                                              dispose ifiles_n ${sta}/dart/hist $*
set obs_seq.*.out;                                                                                                    dispose ifiles_n ${sta}/dart/hist $*
set obs_seq.*.final;                                                                                                  dispose ifiles_n ${sta}/dart/hist $*
set pr*inflate_restart*;  latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null;    dispose ifiles_n ${sta}/dart/rest $*
set po*inflate_restart*;  latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null;    dispose ifiles_n ${sta}/dart/rest $*


IDX=1
while [ $IDX -le $NINST_ATM ]
do
    get_inst_suffix $IDX $NINST_ATM
    set atm${inst_suffix}.log.*;                                                                                                   dispose ifiles_n ${sta}/atm/logs $*
    set ${CASE}.cam${inst_suffix}.r.*;                                                                                             dispose ifiles_y ${sta}/atm/rest $*
    set ${CASE}.cam${inst_suffix}.rs.*;                                                                                            dispose ifiles_y ${sta}/atm/rest $*
    set ${CASE}.cam${inst_suffix}.ra.*;                                                                                            dispose ifiles_y ${sta}/atm/rest $*
    set ${CASE}.cam${inst_suffix}.rh0.*;                                                                                           dispose ifiles_y ${sta}/atm/rest $*
    set ${CASE}.cam${inst_suffix}.rh1.*;                                                                                           dispose ifiles_y ${sta}/atm/rest $*
    set ${CASE}.cam${inst_suffix}.rh2.*;                                                                                           dispose ifiles_y ${sta}/atm/rest $*
    set ${CASE}.cam${inst_suffix}.rh3.*;                                                                                           dispose ifiles_y ${sta}/atm/rest $*
    set ${CASE}.cam${inst_suffix}.rh4.*;                                                                                           dispose ifiles_y ${sta}/atm/rest $*
    set ${CASE}.cam${inst_suffix}.rh5.*;                                                                                           dispose ifiles_y ${sta}/atm/rest $*
    set ${CASE}.cam${inst_suffix}.h0.*;                                                                                            dispose ifiles_n ${sta}/atm/hist $*
    set ${CASE}.cam${inst_suffix}.h1.*;                                                                                            dispose ifiles_n ${sta}/atm/hist $*
    set ${CASE}.cam${inst_suffix}.h2.*;                                                                                            dispose ifiles_n ${sta}/atm/hist $*
    set ${CASE}.cam${inst_suffix}.h3.*;                                                                                            dispose ifiles_n ${sta}/atm/hist $*
    set ${CASE}.cam${inst_suffix}.h4.*;                                                                                            dispose ifiles_n ${sta}/atm/hist $*
    set ${CASE}.cam${inst_suffix}.h5.*;                                                                                            dispose ifiles_n ${sta}/atm/hist $*
    set ${CASE}.cam${inst_suffix}.hs.*;                                                                                            dispose ifiles_n ${sta}/atm/hist $*
    set ${CASE}.cam${inst_suffix}.i.*;    latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_n ${sta}/atm/init $*
    set ${CASE}.camice${inst_suffix}.r.*;                                                                                          dispose ifiles_y ${sta}/ice/rest $*
    set ${CASE}.camdom${inst_suffix}.r.*;                                                                                          dispose ifiles_y ${sta}/ocn/rest $*
    set ${CASE}.camsom${inst_suffix}.r.*;                                                                                          dispose ifiles_y ${sta}/ocn/rest $*
    set ${CASE}.datm${inst_suffix}.r.* ;                                                                                           dispose ifiles_y ${sta}/atm/rest $*
    set ${CASE}.datm${inst_suffix}.rs* ;                                                                                           dispose ifiles_y ${sta}/atm/rest $*
    set ${CASE}.datm${inst_suffix}.h.* ;                                                                                           dispose ifiles_n ${sta}/atm/hist $*
    IDX=`expr $IDX + 1`
done


IDX=1
while [ $IDX -le $NINST_LND ]
do
    get_inst_suffix $IDX $NINST_LND
    set lnd${inst_suffix}.log.*;                                                                                                   dispose ifiles_n ${sta}/lnd/logs $*
    set ${CASE}.clm?${inst_suffix}.r.*;   latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_n ${sta}/lnd/rest $*
    set ${CASE}.clm?${inst_suffix}.rh0.*;                                                                                          dispose ifiles_y ${sta}/lnd/rest $*
    set ${CASE}.clm?${inst_suffix}.rh1.*;                                                                                          dispose ifiles_y ${sta}/lnd/rest $*
    set ${CASE}.clm?${inst_suffix}.rh2.*;                                                                                          dispose ifiles_y ${sta}/lnd/rest $*
    set ${CASE}.clm?${inst_suffix}.rh3.*;                                                                                          dispose ifiles_y ${sta}/lnd/rest $*
    set ${CASE}.clm?${inst_suffix}.rh4.*;                                                                                          dispose ifiles_y ${sta}/lnd/rest $*
    set ${CASE}.clm?${inst_suffix}.rh5.*;                                                                                          dispose ifiles_y ${sta}/lnd/rest $*
    set ${CASE}.clm?${inst_suffix}.h0.*;                                                                                           dispose ifiles_n ${sta}/lnd/hist $*
    set ${CASE}.clm?${inst_suffix}.h1.*;                                                                                           dispose ifiles_n ${sta}/lnd/hist $*
    set ${CASE}.clm?${inst_suffix}.h2.*;                                                                                           dispose ifiles_n ${sta}/lnd/hist $*
    set ${CASE}.clm?${inst_suffix}.h3.*;                                                                                           dispose ifiles_n ${sta}/lnd/hist $*
    set ${CASE}.clm?${inst_suffix}.h4.*;                                                                                           dispose ifiles_n ${sta}/lnd/hist $*
    set ${CASE}.clm?${inst_suffix}.h5.*;                                                                                           dispose ifiles_n ${sta}/lnd/hist $*
    set ${CASE}.clm?${inst_suffix}.hv.*;                                                                                           dispose ifiles_n ${sta}/lnd/hist $*
    set ${CASE}.clm?${inst_suffix}.i.*;                                                                                            dispose ifiles_y ${sta}/lnd/init $*
    set ${CASE}.dlnd${inst_suffix}.r.* ;                                                                                           dispose ifiles_y ${sta}/lnd/rest $*
    set ${CASE}.dlnd${inst_suffix}.rs* ;                                                                                           dispose ifiles_y ${sta}/lnd/rest $*
    set ${CASE}.dlnd${inst_suffix}.h.* ;                                                                                           dispose ifiles_n ${sta}/lnd/hist $*
    IDX=`expr $IDX + 1`
done


IDX=1
while [ $IDX -le $NINST_ICE ]
do
    get_inst_suffix $IDX $NINST_ICE
    set ice${inst_suffix}.log.*;                                                                                                      dispose ifiles_n ${sta}/ice/logs $*
    set ${CASE}.cice${inst_suffix}.r.*;      latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_n ${sta}/ice/rest $*
    set ${CASE}.cice${inst_suffix}.r.[0-9]*;                                                                                          dispose ifiles_y ${sta}/ice/rest $*
    set ${CASE}.cice${inst_suffix}.r.volpn*;                                                                                          dispose ifiles_y ${sta}/ice/rest $*
    set ${CASE}.cice${inst_suffix}.r.dEdd*;                                                                                           dispose ifiles_y ${sta}/ice/rest $*
    set ${CASE}.cice${inst_suffix}.r.age*;                                                                                            dispose ifiles_y ${sta}/ice/rest $*
    set ${CASE}.cice${inst_suffix}.r.aero*;                                                                                           dispose ifiles_y ${sta}/ice/rest $*
    set ${CASE}.cice${inst_suffix}.h*;                                                                                                dispose ifiles_n ${sta}/ice/hist $*
    set ${CASE}.cice${inst_suffix}.i.*;                                                                                               dispose ifiles_y ${sta}/ice/init $*
    set ${CASE}.dice${inst_suffix}.r.* ;                                                                                              dispose ifiles_y ${sta}/ice/rest $*
    set ${CASE}.dice${inst_suffix}.rs* ;                                                                                              dispose ifiles_y ${sta}/ice/rest $*
    set ${CASE}.dice${inst_suffix}.h.* ;                                                                                              dispose ifiles_n ${sta}/ice/hist $*
    IDX=`expr $IDX + 1`
done


IDX=1
while [ $IDX -le $NINST_OCN ]
do
    get_inst_suffix $IDX $NINST_OCN
    set ocn${inst_suffix}.log.*;                                                                                                               dispose ifiles_n ${sta}/ocn/logs $*
    set ${CASE}.pop${inst_suffix}.r.*.hdr;                                                                                                     dispose ifiles_y ${sta}/ocn/rest $*
    set ${CASE}.pop${inst_suffix}.r.*0000;                                                                                                     dispose ifiles_y ${sta}/ocn/rest $*
    set ${CASE}.pop${inst_suffix}.r.*0000.nc;                                                                                                  dispose ifiles_y ${sta}/ocn/rest $*
    set ${CASE}.pop${inst_suffix}.rh.ecosys.*.hdr;                                                                                             dispose ifiles_y ${sta}/ocn/rest $*
    set ${CASE}.pop${inst_suffix}.rh.ecosys.*0000;                                                                                             dispose ifiles_y ${sta}/ocn/rest $*
    set ${CASE}.pop${inst_suffix}.rh.ecosys.*0000.nc;                                                                                          dispose ifiles_y ${sta}/ocn/rest $*
    set ${CASE}.pop${inst_suffix}.rh.*.hdr;                                                                                                    dispose ifiles_y ${sta}/ocn/rest $*
    set ${CASE}.pop${inst_suffix}.rh.*0000;                                                                                                    dispose ifiles_y ${sta}/ocn/rest $*
    set ${CASE}.pop${inst_suffix}.rh.*0000.nc;                                                                                                 dispose ifiles_y ${sta}/ocn/rest $*
    set ${CASE}.pop${inst_suffix}.ro.*;                                                                                                        dispose ifiles_y ${sta}/ocn/rest $*
    set ${CASE}.pop${inst_suffix}.d?*;                                                                                                         dispose ifiles_n ${sta}/ocn/hist $*
    set ${CASE}.pop${inst_suffix}.h*;                                                                                                          dispose ifiles_n ${sta}/ocn/hist $*
    set ${CASE}.docn${inst_suffix}.r.* ;                                                                                                       dispose ifiles_y ${sta}/ocn/rest $*
    set ${CASE}.docn${inst_suffix}.rs* ;                                                                                                       dispose ifiles_y ${sta}/ocn/rest $*
    set ${CASE}.docn${inst_suffix}.h.* ;                                                                                                       dispose ifiles_n ${sta}/ocn/hist $*
    IDX=`expr $IDX + 1`
done


IDX=1
while [ $IDX -le $NINST_GLC ]
do
    get_inst_suffix $IDX $NINST_GLC
    set glc${inst_suffix}.log.*;                                                                                                      dispose ifiles_n ${sta}/glc/logs $*
    set ${CASE}.cism${inst_suffix}.r.[0-9]*;                                                                                          dispose ifiles_y ${sta}/glc/rest $*
    set ${CASE}.cism${inst_suffix}.r.volpn*;                                                                                          dispose ifiles_y ${sta}/glc/rest $*
    set ${CASE}.cism${inst_suffix}.r.dEdd*;                                                                                           dispose ifiles_y ${sta}/glc/rest $*
    set ${CASE}.cism${inst_suffix}.r.age*;                                                                                            dispose ifiles_y ${sta}/glc/rest $*
    set ${CASE}.cism${inst_suffix}.r.aero*;                                                                                           dispose ifiles_y ${sta}/glc/rest $*
    set ${CASE}.cism${inst_suffix}.h*;                                                                                                dispose ifiles_n ${sta}/glc/hist $*
    set ${CASE}.cism${inst_suffix}.i.*;                                                                                               dispose ifiles_y ${sta}/glc/init $*
    IDX=`expr $IDX + 1`
done


cp ${sta}/rest/${dname}/* .

mv ${sta}/* ${DOUT_S_ROOT}
rm -fr ${sta}

echo "st_archive.sh: short-term archiving completed successfully"

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

