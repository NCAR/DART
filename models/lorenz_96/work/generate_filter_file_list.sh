#!/bin/bash

input_projectname=$1
output_projectname=$2
n_members=$3
input_filelist=lorenz_input_file_list.txt
output_filelist=lorenz_output_file_list.txt

echo "Input Ensemble from : ${input_projectname}"
echo "Output Ensemble to : ${output_projectname}"

if [ -f ${input_filelist} ]; then
    rm ${input_filelist}
fi
if [ -f ${output_filelist} ]; then
    rm ${output_filelist}
fi

for i_mem in $(seq -f %05g 1 1 ${n_members}); do
    echo ${input_projectname}/l96_mem${i_mem}.nc  >> ${input_filelist}
#    echo ${output_projectname}/l96_mem${i_mem}.nc  >> ${output_filelist}  # presently unsurported
done
echo ${output_projectname}/l96.nc >> ${output_filelist}

if [ ! -d ${output_projectname} ]; then
	mkdir ${output_projectname}
fi
