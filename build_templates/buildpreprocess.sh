#!/usr/bin/env bash

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

#-------------------------
# Generate build_templates/version_mod.f90 with the version from git describe.
# Arguments:
#  none
# Globals:
#  DART - root of DART
#-------------------------
function dartversion() {

local git_version out tmp

git_version=""
if command -v git >/dev/null 2>&1 && [ -e "$DART/.git" ]; then
   git_version=$(cd "$DART" && git describe --tags --dirty --always 2>/dev/null) || git_version=""
fi

# The version becomes a Fortran character literal.
# Restrict it to characters that cannot break out of the quotes and cap length to 64 characters
git_version=$(printf '%s' "$git_version" | LC_ALL=C tr -c 'A-Za-z0-9._+-' '_' | cut -c1-64)
git_version="${git_version:-unknown}"

out="$DART/build_templates/version_mod.f90"
tmp="$out.$$"

cat > "$tmp" << EOF
! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> Module to provide version information for DART programs
!>
!> Generated at build time by dartversion() in
!> \$DART/build_templates/buildpreprocess.sh - do not edit, do not commit.

module version_mod

implicit none
private

character(len=*), parameter :: dart_version_string = &
   '$git_version'

public :: get_dart_version

contains

! Return DART version in string
function get_dart_version() result(version)

character(len=len(dart_version_string)) :: version

version = dart_version_string

end function get_dart_version

end module version_mod
EOF

# leave the file alone if the version has not changed, so make does not
# rebuild everything that uses version_mod on every build
if cmp -s "$tmp" "$out" ; then
  \rm -f -- "$tmp"
else
  mv -f -- "$tmp" "$out"
fi

}

#-------------------------
# Build and run preprocess
# Arguements:
#  none
# Globals:
#  DART - root of DART
#-------------------------
function buildpreprocess() {

 local pp_dir=$DART/assimilation_code/programs/preprocess
 dartversion
 # run preprocess if it is in the current directory
 if [ -f preprocess ]; then
   ./preprocess
   return
 fi

# link to preprocess if it is already built, run
if [ -f $pp_dir/preprocess ]; then
   ln -s $pp_dir/preprocess .
   ./preprocess 
   return
fi

 # build preprocess, link, run
 cd $pp_dir
 $DART/build_templates/mkmf -x -p $pp_dir/preprocess \
      -a $DART $pp_dir/path_names_preprocess
 cd -
 ln -s $pp_dir/preprocess .
 ./preprocess
}

#-------------------------
# clean up *.mod *.o for preprocess
#-------------------------
function cleanpreprocess() {

 local pp_dir=$DART/assimilation_code/programs/preprocess
 cd $pp_dir
 \rm -f -- *.o *.mod Makefile preprocess
 cd -

}

