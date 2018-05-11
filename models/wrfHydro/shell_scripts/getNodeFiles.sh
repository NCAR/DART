#! /bin/bash

currentDir=$1

cd $currentDir || exit 1
pwd
nodeDir=`head -1 nodeDir.control` || exit 1
echo $nodeDir

for uu in `tail -n+2 nodeDir.control | sort | uniq`
do
    echo $uu 
    #ssh $uu ls ${nodeDir}/OUTPUT/
    ssh $uu mv ${nodeDir}/OUTPUT/\*  ${currentDir}/OUTPUT/.
    ssh $uu rm -rf ${nodeDir}
done 


exit 0
