! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> This program assists in constructing a table which can be read
!> by filter at run-time to disable or alter how the assimilation
!> of different types of observations impact the state vector values
!> based on their kind.  This tool allows users to group related
!> collections of observation types and state vector kinds by name
!> and then express the relationship of the named groups to each
!> other in a concise way.
!>
!> At run time, filter can read the output file from this tool
!> and use it to control the impact at assimilation time based
!> on the relationships specified.
!>
!> The first version of this tool requires the last numeric
!> column to be 0.0, but future extensions may allow different
!> values to be used.
!>
!> All the listed observation types and state vector kinds
!> must be known by the system.  If they are not, look at the
!> &preprocess_nml :: input_items namelist which specifies
!> which obs_def_xxx_mod.f90 files are included, which is
!> where observation types are defined.  kinds are defined
!> in the obs_kinds/DEFAULT_obs_kinds_mod.F90 file and
!> are static.  (note you must add new kinds in 2 places 
!> if you do alter this file.)
!>



! program to read an ascii file with directions for which state and observation
! KINDS should impact which other state and observation KINDS.

! the format of the ascii input file is:
!
! # rest of line is comment after hash mark
! GROUP groupname1
!  KIND_xxx  KIND_xxx  KIND_xxx
!  KIND_xxx
! END GROUP
!
! GROUP groupname2
!  KIND_xxx  KIND_xxx  KIND_xxx
!  KIND_xxx
! END GROUP
!
! GROUP groupnameM
! ALL EXCEPT KIND_xxx KIND_xxx
! KIND_xxx
! END GROUP
! # to choose all kinds except a select few
!
! # FIXME: this is not supported yet - should it be?  yes.
! GROUP groupnameN
! ALL EXCEPT groupnameY
! END GROUP
!
! also ALLTYPES, ALLKINDS, as well as ALL
!
! IMPACT
!  KIND_xxx    KIND_xxx     0.0
!  KIND_xxx    groupname1   0.0
!  groupname1  KIND_xxx     0.0
!  groupname1  groupname1   0.0
! END IMPACT

! # this also is not supported yet, should it be?
! GROUP groupnameX
! # not only kinds, but other groups, recursively?
!  different_groupname  # but what about loops?  even possible?
! END GROUP

! alternative output is triplets of 'kind1 kind2 value'
! which gets put into the 2d array at run time?  is that better?

! maybe prototype the application of these values inside assim_tools
! first before deciding what the output should be.  inside the
! filter run i think it should be a 2d array for speed, but is that
! the most flexible for input to filter?

! or finally, the input to this tool could be the ascii text file,
! and the 2d table could be computed and filled in the initialization
! code (if not too slow) and there is no separate tool needed?

! i started to add types, but they can only appear in column 1 of the
! triplet line.  rethink this.

! the output of this tool is a ascii file containing lines:
!  KIND1_string   KIND2_string    0.0
!
! the 2d table building will now be done at run time in assim_tools
! use would be if (impact_table(kind1, kind2) > 0.0) then ok else 
! change the increments. this could be done in assim_tools in a generic way.

program obs_impact_tool

use      types_mod, only : r8
use  utilities_mod, only : register_module, initialize_utilities, finalize_utilities
use obs_impact_mod, only : create_impact_table

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"
character(len=128), parameter :: id  = "$Id$"


real(r8), allocatable :: table(:,:)

! initialization and setup

call initialize_utilities('obs_impact_tool')
call register_module(id)

! build and output impact_table
call create_impact_table()

! clean up
call finalize_utilities('obs_impact_tool')

end program

