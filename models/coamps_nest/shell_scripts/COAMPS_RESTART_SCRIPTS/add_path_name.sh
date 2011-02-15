#!/bin/bash
#
# SCRIPT:   add_path_name.sh
# AUTHOR:   T. R. Whitcomb
#           Naval Research Laboratory
#
# Adds a file name to all entries with path_names at the beginning -
# use this for adding new files into the DART mkmf framework.

for pathfile in `ls -1 path_names*`
do
  echo "$1" >> $pathfile
done
