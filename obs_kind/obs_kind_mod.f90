! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module obs_kind_mod

! <next five lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
! $Name$ 

! This module is designed to provide general information about observation 
! types. It is not clear at present whether this is going to be viable or
! not given the requirement of underlying models and domains. For now, kind
! simply is an index that identifies what kind of observation this is 
! for instance temperature, pressure, satellite radiance, etc. This may be
! a good place to store other parameters associated with more complicated
! observation types.

! Revised obs_kind treats identity observations (including the index) as 
! just another kind of observation. Identity observations of an extended 
! state vector are just represented by indices that go beyond the state
! vector size into the extended state vector.
! Initial revision 15 April, 2004

use utilities_mod, only : register_module, error_handler, E_ERR

implicit none
private

public :: KIND_U, KIND_V, KIND_PS, KIND_T, KIND_QV, KIND_P, KIND_W, &
          KIND_QR, KIND_RHO, KIND_U10, KIND_V10, KIND_T2, KIND_Q2, &
          KIND_QG, KIND_QS, &
          KIND_VR, KIND_REF, KIND_TD, KIND_TD2

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

! ADD A LONG TABLE OF DEFINED BUFR INDICES, ETC.
! Definition of observation kind types:

integer, parameter :: KIND_U   =   1,  & ! zonal wind component
                      KIND_V   =   2,  & ! meridional wind component
                      KIND_PS  =   3,  & ! Surface pressure
                      KIND_T   =   4,  & ! Temperature
                      KIND_QV  =   5,  & ! Specific humidity (mixing ratio)
                      KIND_P   =   6,  & ! Pressure
                      KIND_W   =   7,  & ! Vertical velocity
                      KIND_QR  =   8,  & ! Rainwater mixing ratio
                      KIND_TD  =  10,  & ! Dew point temperature
                      KIND_RHO =  11,  & ! Density
                      KIND_VR  = 100,  & ! Doppler radar radial velocity
                      KIND_REF = 101,  & ! Radar reflectivity
                      KIND_U10 =  14,  & ! zonal wind component at 10 m AGL
                      KIND_V10 =  15,  & ! meridional wind component at 10 m AGL
                      KIND_T2  =  16,  & ! Temperature at 2 m AGL
                      KIND_Q2  =  17,  & ! Specific humidity (mixing ratio) at 2 m AGL
                      KIND_QG  =  20,  & ! Graupel mixing ratio
                      KIND_QS  =  21,  & ! Snow mixing ratio
                      KIND_TD2 = 204     ! Dew point temperature at 2 m AGL

    ! NOTICE   at one point   KIND_U10 = 200
    ! NOTICE   at one point   KIND_V10 = 201
    ! NOTICE   at one point   KIND_T2  = 202
    ! NOTICE   at one point   KIND_Q2  = 203

logical, save :: module_initialized = .false.

contains


!----------------------------------------------------------------------------

subroutine initialize_module

   call register_module(source, revision, revdate)
   module_initialized = .true.

end subroutine initialize_module


end module obs_kind_mod
