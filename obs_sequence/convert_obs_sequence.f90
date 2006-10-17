! Data Assimilation Research Testbed -- DART
! Copyright 2004-2006, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program convert_obs_sequence

! <next five lines automatically updated by CVS, do not edit>                             
! $Source$         
! $Revision$                                                                       
! $Date$                                                            
! $Author$
! $Name$

use    utilities_mod, only : timestamp, register_module, initialize_utilities 
use obs_sequence_mod, only : read_obs_seq, write_obs_seq, &
                             static_init_obs_sequence, obs_sequence_type
use  assim_model_mod, only : static_init_assim_model

implicit none

! CVS Generated file description for error handling, do not edit                          
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

type(obs_sequence_type) :: seq
character(len = 129)    :: file_name_in, file_name_out

! Record the current time, date, etc. to the logfile
call initialize_utilities('convert_obs_sequence')
call register_module(source,revision,revdate)

! Initialize the assim_model module, need this to get model
! state meta data for locations of identity observations
call static_init_assim_model()

! Initialize the obs_sequence module
call static_init_obs_sequence()

! Read the existing sequence in
write(*, *) 'Input filename for sequence ( obs_seq.out or obs_seq.final are common )'
read(*, *) file_name_in
call read_obs_seq(file_name_in, 0, 0, 0, seq)


! Write the sequence to a file
write(*, *) 'Output filename for sequence ( must be different than input name )'
read(*, *) file_name_out
call write_obs_seq(seq, file_name_out)

! Clean up
call timestamp(string1=source,string2=revision,string3=revdate,pos='end')

end program convert_obs_sequence
