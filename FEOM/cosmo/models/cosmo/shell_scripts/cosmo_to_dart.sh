#! /bin/sh

cosmodir=$1
dartdir=$2
filein=$3
fileout=$4
is=$5
ie=$6
ndin=$7
ndout=$8

rundir=/home/jkeller/DART_COSMO/models/cosmo/run
progdir=/home/jkeller/DART_COSMO/models/cosmo/work

dir=`pwd`

for ifile in `seq ${is} ${ie}`; do
    mdir=${cosmodir}/m`printf %${ndin}.${ndin}i $ifile`

    in=${mdir}/${filein}
    out=${dartdir}/${fileout}`printf %${ndout}.${ndout}i $ifile`

    rm ${rundir}/script1.sed
    echo "s#OLDCOSMOFILE#${in}#g"   > ${rundir}/script1.sed
    echo "s#NEWCOSMOFILE#dummy#g"  >> ${rundir}/script1.sed
    echo "s#DARTFILE#${out}#g"     >> ${rundir}/script1.sed

    sed -f ${rundir}/script1.sed ${rundir}/input.nml.template > ${rundir}/input.nml
    cat ${rundir}/input.nml
    exit
    cd ${rundir}
    ${progdir}/cosmo_to_dart
    cd ${dir}
done
