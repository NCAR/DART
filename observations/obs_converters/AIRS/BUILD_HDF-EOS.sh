#!/bin/sh
# 
# updated 4 Dec 2020

echo 
echo  'These converters require either the HDF-EOS or the HDF-EOS5 libraries.'
echo  'These libraries are, in general, not compatible with each other.'
echo  'There is a compatibility library that "provides uniform access to HDF-EOS2'
echo  'and 5 files though one set of API calls." which sounds great.'
echo 
echo  'The HDF-EOS5 libraries are installed on the supercomputers, and are'
echo  'available via MacPorts (hdfeos5). The HDF-EOS libraries are older and'
echo  'are much less available. Consequently, I have used the HDF-EOS5 interfaces'
echo  'where possible.'
echo 
echo  'If the he5_hdfeos libraries are installed on your system, you are in luck.'
echo  'On our system, it has been useful to define variables like:'
echo 
echo  'setenv("NCAR_INC_HDFEOS5",    "/glade/u/apps/ch/opt/hdf-eos5/5.1.16/intel/19.0.5/include")'
echo  'setenv("NCAR_LDFLAGS_HDFEOS5","/glade/u/apps/ch/opt/hdf-eos5/5.1.16/intel/19.0.5/lib")'
echo  'setenv("NCAR_LIBS_HDFEOS5","-Wl,-Bstatic -lGctp -lhe5_hdfeos -lsz -lz -Wl,-Bdynamic")'
echo  'which we then use in when compiling convert_airs_L2'
echo 
echo  'If you need to build the HDF-EOS and/or the HDF-EOS5 libraries, you may '
echo  'try to follow the steps outlined in this script. They will need to be '
echo  'modified for your system.'
echo 
echo  'You will have to edit this script, first, by removing the early exit ...'
echo

exit

# ------------------------------------------------------------------------------
##
## The NASA Earthdata Data Access Services portal serves as the download site:
## https://wiki.earthdata.nasa.gov/display/DAS/Toolkit+Downloads
##
## The following packages were downloaded:
##
##    zlib-1.2.11.tar.gz
##    jpegsrc.v9b.tar.gz
##    hdf-4.2.13.tar.gz
##    HDF-EOS2.20v1.00.tar.Z
##    HDF-EOS2.20v1.00_TestDriver.tar.Z
##    szip-2.1.1.tar.gz
##    hdf5-1.8.19.tar.gz
##    HDF-EOS5-5.1.16.tar.Z
##    HDF-EOS5-5.1.16_TESTDRIVERS.tar.Z
##
## The documentation files were downloaded:
##
##    HDF-EOS_REF.pdf
##    HDF-EOS_UG.pdf
##    HDF-EOS5_REF.pdf
##    HDF-EOS5_UG.pdf
##
## Some other useful websites for HDF and HDF-related products are:
## https://portal.hdfgroup.org/display/support/Downloads
## https://hdfeos.org/software/library.php#HDF-EOS2
## https://opensource.gsfc.nasa.gov/projects/HDF-EOS2/index.php

# Change this to 'true' to uncompress the packages. You only need to uncompress them
# once, but you may need to run this script several times.

if ( `false` ); then
  
  for i in zlib-1.2.11.tar.gz \
           jpegsrc.v9b.tar.gz \
           hdf-4.2.13.tar.gz \
           szip-2.1.1.tar.gz \
           hdf5-1.8.19.tar.gz  
  do
    tar -zxvf $i
  done

  uncompress HDF-EOS2.20v1.00.tar.Z
  uncompress HDF-EOS2.20v1.00_TestDriver.tar.Z
  uncompress HDF-EOS5.1.16.tar.Z
  uncompress HDF-EOS5.1.16_TESTDRIVERS.tar.Z

  tar -xvf HDF-EOS2.20v1.00.tar
  tar -xvf HDF-EOS5-5.1.16.tar

fi

# ------------------------------------------------------------------------------
# start with smaller libs, work up to HDF-EOS.
# ------------------------------------------------------------------------------

# set the installation location of the final libraries
H4_PREFIX=/glade/work/${USER}/local/hdf-eos
H5_PREFIX=/glade/work/${USER}/local/hdf-eos5

# make the target install dirs
mkdir -p ${H4_PREFIX}/{lib,bin,include,man,man/man1,share} 
mkdir -p ${H5_PREFIX}/{lib,bin,include,man,man/man1,share} 

# record the build script and environment
echo                         >  ${H4_PREFIX}/build_environment_log.txt
echo 'the build script'      >> ${H4_PREFIX}/build_environment_log.txt
cat    $0                    >> ${H4_PREFIX}/build_environment_log.txt
echo                         >> ${H4_PREFIX}/build_environment_log.txt
echo '=====================' >> ${H4_PREFIX}/build_environment_log.txt
echo 'the build environment' >> ${H4_PREFIX}/build_environment_log.txt
echo                         >> ${H4_PREFIX}/build_environment_log.txt
env | sort                   >> ${H4_PREFIX}/build_environment_log.txt

# start with smaller libs, work up to HDF-EOS.

echo ''
echo '======================================================================'
if [ -f ${H4_PREFIX}/lib/libz.a ]; then
   echo 'zlib already exists - no need to build.'
else

   export CC='icc'
   export CFLAGS='-fPIC'
   export FFLAGS='-fPIC'

   echo 'building zlib at '`date`
   cd zlib-1.2.11 || exit 1
   ./configure --prefix=${H4_PREFIX} || exit 1
   make clean     || exit 1
   make           || exit 1
   make test      || exit 1
   make install   || exit 1
   cd ..
fi


echo ''
echo '======================================================================'
if [ -f ${H4_PREFIX}/lib/libsz.a ]; then
   echo 'szip already exists - no need to build.'
else

   export CC='icc'
   export CFLAGS='-fPIC'
   export FFLAGS='-fPIC'

   echo 'building szip at '`date`
   cd szip-2.1.1 || exit 1
   ./configure --prefix=${H4_PREFIX} || exit 1
   make clean     || exit 1
   make           || exit 1
   make test      || exit 1
   make install   || exit 1
   cd ..
fi

echo ''
echo '======================================================================'
# This is peculiar - on Cheyenne:
# If I build with --libdir=H4_PREFIX, subsequent linking works.
# If I build with --libdir=H4_PREFIX/lib, subsequent linking FAILS with an 
# undefined reference to 'rpl_malloc'.
if [ -f ${H4_PREFIX}/lib/libjpeg.a ]; then
   echo 'jpeg already exists - no need to build.'
else
   echo 'buiding jpeg at '`date`
   cd jpeg-9b        || exit 2
   ./configure CC='icc -Df2cFortran' CFLAGS='-fPIC' FFLAGS='-fPIC' \
               --prefix=${H4_PREFIX} || exit 2
   make clean        || exit 2
   make              || exit 2
   make test         || exit 2
   make install      || exit 2
   cd ..
   cd ${H4_PREFIX}
   \ln -s lib/libjpeg* .
   cd -
fi
 
echo ''
echo '======================================================================'
if [ -f ${H4_PREFIX}/lib/libmfhdf.a ]; then
   echo 'hdf4 already exists - no need to build.'
else
   echo 'building hdf4 at '`date`
   # (apparently there is no 'make test')
   
   cd hdf-4.2.13 || exit 3
   ./configure CC='icc -Df2cFortran' CFLAGS='-fPIC' FFLAGS='-fPIC' \
               --prefix=${H4_PREFIX} \
               --disable-netcdf \
               --with-zlib=${H4_PREFIX} \
               --with-jpeg=${H4_PREFIX} || exit 3
   make clean    || exit 3
   make          || exit 3
   make install  || exit 3
   cd ..
fi
 
echo ''
echo '======================================================================'
if [ -f ${H4_PREFIX}/lib/libhdfeos.a ]; then
   echo 'hdf-eos already exists - no need to build.'
else
   echo 'building HDF-EOS2.20v1.00 at '`date`
   echo 'after expanding the .tar.gz file, the source is in "hdfeos"'
   cd hdfeos    || exit 4
   # (the CC options are crucial to provide Fortran interoperability)
   ./configure CC='icc -Df2cFortran' CFLAGS='-fPIC' FFLAGS='-fPIC' \
               --prefix=${H4_PREFIX} \
               --enable-install-include \
               --with-zlib=${H4_PREFIX} \
               --with-jpeg=${H4_PREFIX} \
               --with-hdf=${H4_PREFIX} || exit 4
   make clean   || exit 4
   make         || exit 4
   make install || exit 4
   cd ..
fi

#-------------------------------------------------------------------------------
# HDF-EOS5 record the build script and environment
#-------------------------------------------------------------------------------

echo                         >  ${H5_PREFIX}/build_environment_log.txt
echo 'the build script'      >> ${H5_PREFIX}/build_environment_log.txt
cat    $0                    >> ${H5_PREFIX}/build_environment_log.txt
echo                         >> ${H5_PREFIX}/build_environment_log.txt
echo '=====================' >> ${H5_PREFIX}/build_environment_log.txt
echo 'the build environment' >> ${H5_PREFIX}/build_environment_log.txt
echo                         >> ${H5_PREFIX}/build_environment_log.txt
env | sort                   >> ${H5_PREFIX}/build_environment_log.txt

echo '======================================================================'
if [ -f ${H5_PREFIX}/lib/libhdf5.a ]; then
   echo 'hdf5 already exists - no need to build.'
else
   echo 'building hdf5 at '`date`
   
   cd hdf5-1.8.19 || exit 3
   ./configure CC='icc -Df2cFortran' CFLAGS='-fPIC' FFLAGS='-fPIC' \
               --prefix=${H5_PREFIX} \
               --enable-fortran \
               --enable-fortran2003 \
               --enable-production \
               --with-zlib=${H4_PREFIX} || exit 3
   make clean          || exit 3
   make                || exit 3
   make check          || exit 3
   make install        || exit 3
   make check-install  || exit 3
   cd ..
fi

echo ''
echo '======================================================================'
if [ -f ${H5_PREFIX}/lib/libhe5_hdfeos.a ]; then
   echo 'hdf-eos5 already exists - no need to build.'
else
   echo 'building HDF-EOS5.1.16 at '`date`
   echo 'after expanding the .tar.Z file, the source is in "hdfeos5"'
   cd hdfeos5    || exit 4
   # (the CC options are crucial to provide Fortran interoperability)
   ./configure CC='icc -Df2cFortran' CFLAGS='-fPIC' FFLAGS='-fPIC' \
               --prefix=${H5_PREFIX} \
               --enable-install-include \
               --with-zlib=${H4_PREFIX} \
               --with-hdf5=${H5_PREFIX} || exit 4
   make clean   || exit 4
   make         || exit 4
   make check   || exit 4
   make install || exit 4
   cd ..
fi

exit 0
