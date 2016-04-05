#! /bin/sh

cosmodir=$1
sdate=$2
nhours=$3
is=$4
ie=$5
nd=$6
nnodes=$7

nprocs=$(($nnodes * 4))

maxjobs=3

rundir=/home/jkeller/DART_COSMO/models/cosmo/run
progdir=/home/jkeller/DART_COSMO/models/cosmo/work

bdir=${cosmodir}
execdir=/daten03/cosmonrw/cosmo/cosmo_4.11uc_orig

username=`whoami`

for irun in `seq ${is} ${ie}`; do
    rfile=${rundir}/finished.`printf %${nd}.${nd}i $irun`

    mdir=${cosmodir}/m`printf %${nd}.${nd}i $irun`
    cd ${mdir}

    rm -rf ${rundir}/script1.sed
    echo "s#STARTDATE#${sdate}#g"   > ${rundir}/script1.sed
    echo "s#NHOURS#${nhours}#g"     >> ${rundir}/script1.sed
    echo "s#MEMBERDIR#${mdir}#g"    >> ${rundir}/script1.sed
    echo "s#BOUNDARYDIR#${bdir}#g"  >> ${rundir}/script1.sed
    echo "s#NNODES#${nnodes}#g"     >> ${rundir}/script1.sed
    echo "s#NPROCS#${nprocs}#g"     >> ${rundir}/script1.sed
    echo "s#EXECDIR#${execdir}#g"   >> ${rundir}/script1.sed
    echo "s#READYFILE#${rfile}#g"   >> ${rundir}/script1.sed
    echo "s#USERNAME#${username}#g"    >> ${rundir}/script1.sed

    sed -f ${rundir}/script1.sed ${rundir}/cosmo_job.template > ${mdir}/cosmo_job
    echo "cat 'sleep 5'" > test_job
    echo "cat 'finished' \> ${rfile}" >> test_job

    rm ${rfile}
    cd ${dir}
done

irun=$is

while [ ! -s ${rfile} ]; do
    while [ $irun -le $ie ]; do
        njobs=`qstat | grep ${username} | wc -l`
	if [ $njobs -lt $maxjobs ]; then
	    mdir=${cosmodir}/m`printf %${nd}.${nd}i $irun`
	    cd ${mdir}
	    qsub cosmo_job
	    cd ${dir}
	    irun=$(($irun+1))
	else
	    sleep 10
	fi
    done
    exit
done
