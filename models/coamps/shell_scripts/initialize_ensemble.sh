#!/bin/bash
#
# This code may (or may not) be part of the COAMPS distribution,
# So it is not protected by the DART copyright agreement.
#
# DART $Id$
#
# AUTHOR:   T. R. Whitcomb
#           Naval Research Laboratory
#
# Creates and populates the ensemble member subdirectories for a
# COAMPS ensemble - each subdirectory contains everything needed
# to run the model.
#
# No command line arguments are required.
######

COAMPS_DIR=/coamps/data/directory/goes/here

# How many members?
NUM_MEMS=10

for dir in `seq 1 ${NUM_MEMS}`
do
  # Generate the name of the directory currently holding COAMPS data.  
  # Information will be copied *FROM* this directory into MEMBER_DIR.
  # 
  # The printf statement allows the substitution of the ensemble member
  # number into COAMPS_DIR if needed, otherwise it leaves it as is.
  # It goes without saying that you should not use any character
  # string that printf will recognize in your directory names.
  MEMBER_DATA_DIR=`printf ${COAMPS_DIR} ${dir}`
  
  # Recast the ensemble member number to a directory name - just a 
  # fixed number of digits is the default
  MEMBER_DIR=`printf "%05d" ${dir}`
 
  # These are tasks that only have to be done once
  if [ ! -e ${MEMBER_DIR} ]
      then
      # create the ensemble directory if it's not there
      # This could be done with just the second command, but I
      # think this more explicit version is clearer
      echo "Creating directory for ensemble member ${MEMBER_DIR}..."
      mkdir -p ${MEMBER_DIR}
      mkdir -p ${MEMBER_DIR}/data
      mkdir -p ${MEMBER_DIR}/data/backward
      mkdir -p ${MEMBER_DIR}/data/forward

      # Copy in the restart file, boundary conditions, and 
      # COAMPS domain information files - only copy the boundary
      # conditions over if we need them (i.e. not peridic)
      # Modify this to search for one boundary file in particular
      # since if the wildcard expands, this throws an error
      if [ -e ${MEMBER_DATA_DIR}/bdwwnd_sig_*1a*00000000_bndyfld ]
        then
        cp -f ${MEMBER_DATA_DIR}/*bndyfld ${MEMBER_DIR}/data
      fi
      cp -f ${MEMBER_DATA_DIR}/datahd*  ${MEMBER_DIR}/data
      cp -f ${MEMBER_DATA_DIR}/terrht*  ${MEMBER_DIR}/data

      # We also need all the surface fields at initial time.  This
      # includes things like land/sea flags, surface roughness, etc.
      cp -f ${MEMBER_DATA_DIR}/*sfc*00000000*  ${MEMBER_DIR}/data
  fi

  # Copy our template namelist into each ensemble member's
  # directory and replace so it finds its own data
  echo "Copying COAMPS namelist for ensemble member ${MEMBER_DIR}"
  cp -f namelist ${MEMBER_DIR}/namelist
  sed -i.bak "s/ENSDIR/${MEMBER_DIR}/" ${MEMBER_DIR}/namelist
  rm ${MEMBER_DIR}/*.bak

  # Copy 
  echo "Copying pristine restart file for ensemble member ${MEMBER_DIR}"
  rsync -v ${MEMBER_DATA_DIR}/restart* ${MEMBER_DIR}/data

  # Copy in the restart.vars file that identifies the members of the
  # DART state vector.  We also need the DART input.nml file here.
  ln -sf `pwd`/restart.vars ${MEMBER_DIR}/data/restart.vars
  ln -sf `pwd`/convert.vars ${MEMBER_DIR}/data/convert.vars
  ln -sf `pwd`/input.nml ${MEMBER_DIR}/data/input.nml

done

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

