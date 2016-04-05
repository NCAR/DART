#!/bin/csh


# run a single instance of prepbufr.  
# requires several arguments: 
# 1. thedirectory to cd into first.
# 2. year
# 3. month
# 4. day


mkdir $1
cd $1
cp -f ../prepbufr.csh .
cp -f ../input.nml .
./prepbufr.csh $2 $3 $4 $4
rm input.nml 
rm prepbufr.csh
cd ..

exit 0
