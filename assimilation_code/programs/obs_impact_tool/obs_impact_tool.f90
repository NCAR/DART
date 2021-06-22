! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> This program users to construct a table that is read by filter at 
!> run-time to localize the impact of sets of observation types on sets of 
!> state vector quantities. 

!> The standard DART algorithms compute increments for an observation and then
!> compute corresponding increments for each model state variable due to that
!> observation. To do this, DART computes a sample regression coefficient using
!> the prior ensemble distributions of a state variable and the observation. The
!> increments for each member of the observation are multiplied by this
!> regression coefficient and then added to the corresponding prior ensemble
!> member for the state variable. However, in many cases, it is appropriate to
!> reduce the impact of an observation on a state variable; this is called
!> localization. The standard DART algorithms allow users to specify a
!> localization that is a function of the horizontal (and optionally vertical)
!> distance between the observation and the state variable. The localization is
!> a value between 0 and 1 and multiplies the regression coefficient when
!> updating state ensemble members.
!>
!> Sometimes, it may be desirable to do an additional localization that is a
!> function of the type of observation and the state vector quantity. This
!> program allows users to construct a table that is read by filter at run-time
!> to localize the impact of sets of observation types on sets of state
!> vectorquantities. Users can create named sets of observation types and sets
!> of state vector quantities and specify a localization for the impact of the
!> specified observation types on the state vector quantities.
!>
!> An example would be to create a subset of observations of tracer
!> concentration for a variety of tracers, and a subset of dynamic state
!> variable quantities like temperatures and wind components. It has been common
!> to set this localization value to 0 so that tracer observations have no
!> impact on dynamic state quantities, however, the tool allows values between 0
!> and 1 to be specified.
!> 
!> All the listed observation types and state vector quantities
!> must be known by the system.  If they are not, look at the
!> &preprocess_nml :: input_items namelist which specifies
!> which obs_def_xxx_mod.f90 files are included, which is
!> where observation types are defined.  Quantities are defined
!> in the assimilation_code/modules/observations/DEFAULT_obs_kinds_mod.F90 file.
!> (Note you must add new quantities in 2 places 
!> if you do alter this file.)




! the format of the ascii input file is:
!
! # rest of line is comment after hash mark
! GROUP groupname1
!  QTY_xxx  QTY_xxx  QTY_xxx
!  QTY_xxx
! END GROUP
!
! GROUP groupname2
!  QTY_xxx  QTY_xxx  QTY_xxx
!  QTY_xxx
! END GROUP
!
! GROUP groupnameM
!  ALL EXCEPT QTY_xxx QTY_xxx
!  QTY_xxx
! END GROUP
!
! # to choose all quantities except a select few
! GROUP groupnameN
!  ALL EXCEPT groupnameY
! END GROUP
!
! also ALLTYPES, ALLQTYS, as well as ALL
!
! IMPACT
!  QTY_xxx     QTY_xxx      0.0
!  QTY_xxx     groupname1   0.0
!  groupname1  QTY_xxx      0.0
!  groupname1  groupname2   0.0
! END IMPACT

! GROUP groupnameX
!  different_groupname  # no circular dependencies allowed
! END GROUP

! the output of this tool is a ascii file containing lines:
!  QTY1_string   QTY2_string    0.0
!

program obs_impact_tool

use      types_mod, only : r8
use  utilities_mod, only : initialize_utilities, finalize_utilities, &
                           find_namelist_in_file, check_namelist_read, E_MSG,         &
                           do_nml_file, do_nml_term, nmlfileunit, error_handler
use obs_impact_mod, only : create_impact_table

integer :: funit, ios

! namelist: input/output names, values, etc
character(len=512) :: input_filename  = ''
character(len=512) :: output_filename = ''
logical :: debug = .false.  ! .true. for more output

! namelist
namelist /obs_impact_tool_nml/  &
   input_filename,  &
   output_filename, &
   debug


! initialization and setup

call initialize_utilities('obs_impact_tool')

call find_namelist_in_file("input.nml", "obs_impact_tool_nml", funit)
read(funit, nml = obs_impact_tool_nml, iostat = ios)
call check_namelist_read(funit, ios, "obs_impact_tool_nml")

if (do_nml_file()) write(nmlfileunit, nml=obs_impact_tool_nml)
if (do_nml_term()) write(     *     , nml=obs_impact_tool_nml)

if (debug) call error_handler(E_MSG, 'obs_impact_tool', ' debug on')


! build and output impact_table
call create_impact_table(input_filename, output_filename, debug) 

! clean up
call finalize_utilities('obs_impact_tool')

end program

