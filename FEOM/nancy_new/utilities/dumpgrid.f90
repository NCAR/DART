! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program dumpgrid

! Program to output the locations and kinds for any state vector to standard output.
! Uses the required 'get_state_meta_data()' function in the model_mod.f90
! and so will work with any model in DART.

use utilities_mod,     only : register_module, error_handler, E_MSG, E_ERR,  &
                              open_file, close_file, initialize_utilities,   &
                              finalize_utilities

use  location_mod,     only : location_type, write_location

use  obs_kind_mod,     only : get_raw_obs_kind_name, paramname_length

use assim_model_mod,   only : static_init_assim_model, get_model_size,   &
                              get_state_meta_data


implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"
character(len=128), parameter :: id  = "$Id$"

integer               :: model_size, i
integer               :: itemkind
type(location_type)   :: itemloc
character(len=128)    :: locstring
character(len=paramname_length) :: pstring


!----------------------------------------------------------------
! program start
!----------------------------------------------------------------

 
call initialize_utilities('dumpgrid')
!call register_module(id)


! Initialize the model so we can get the size and metadata.
call static_init_assim_model()
model_size = get_model_size()

write(*, *) 'Model size =', model_size

! get the locations and kinds for each item in the state vector.

do i=1, model_size
   call get_state_meta_data(i, itemloc, itemkind)
   call write_location(0, itemloc, charstring=locstring)
   pstring = get_raw_obs_kind_name(itemkind)

   write(*, *) i, trim(locstring)//', '//pstring
enddo


call finalize_utilities() 

end program dumpgrid

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
