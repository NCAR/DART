#! /bin/sh

dartdir=$1
cosmodir=$2
filein=$3
fileout=$4
is=$5
ie=$6
ndin=$7
ndout=$8

rundir=/home/jkeller/DART_COSMO/models/cosmo/run
progdir=/home/jkeller/DART_COSMO/models/cosmo/work

dir=`pwd`
mdir=${cosmodir}/m`printf %${ndout}.${ndout}i 0`
if [ ! -s $mdir ]; then
    mkdir $mdir
fi
if [ ! -s $mdir/${fileout} ]; then
    cp ${cosmodir}/${fileout} ${mdir}
fi
laf=${mdir}/${fileout}

for ifile in `seq ${is} ${ie}`; do
    mdir=${cosmodir}/m`printf %${ndout}.${ndout}i $ifile`
    if [ ! -s $mdir ]; then
	mkdir $mdir
    fi
    in=${dartdir}/${filein}`printf %${ndin}.${ndin}i $ifile`
    out=${mdir}/${fileout}
    rm -rf ${out}
    rm -rf ${rundir}/script1.sed
    echo "s#OLDCOSMOFILE#${laf}#g"  > ${rundir}/script1.sed
    echo "s#NEWCOSMOFILE#${out}#g"  >> ${rundir}/script1.sed
    echo "s#DARTFILE#${in}#g"       >> ${rundir}/script1.sed

    sed -f ${rundir}/script1.sed ${rundir}/input.nml.template > ${rundir}/input.nml
    cd ${rundir}
    ${progdir}/dart_to_cosmo
    cd ${dir}
done
