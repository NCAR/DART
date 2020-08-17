! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module dart_gitm_mod

! This is the interface between the GITM modules and DART.
! To reduce the possibility of scoping issues, all the
! unrestricted GITM modules are confined to this module.

use ModConstants
use ModSizeGitm
use ModPlanet 

use typesizes
use netcdf

use    utilities_mod, only : error_handler, E_ERR, E_WARN, E_MSG

implicit none
private

! these routines must be public and you cannot change
! the arguments - they will be called *from* the DART code.
public :: get_nLatsPerBlock, &
          get_nLonsPerBlock, &
          get_nAltsPerBlock, &
          get_nSpecies,      &
          get_nSpeciesTotal, &
          get_nIons,         &
          get_nSpeciesAll,   &
          decode_gitm_indices

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

character(len=256) :: string1, string2

contains

!===================================================================
! All the public interfaces ... nothing more.
!===================================================================

! @todo FIXME - should this now get the sizes from the netcdf file
! and not include GITM code?  (i think yes.)

integer function get_nLatsPerBlock()
   get_nLatsPerBlock = nLats
end function get_nLatsPerBlock

integer function get_nLonsPerBlock()
   get_nLonsPerBlock = nLons
end function get_nLonsPerBlock

integer function get_nAltsPerBlock()
   get_nAltsPerBlock = nAlts
end function get_nAltsPerBlock

integer function get_nSpecies()
   get_nSpecies = nSpecies   ! From ModPlanet, hopefully
end function get_nSpecies

integer function get_nSpeciesTotal()
   get_nSpeciesTotal = nSpeciesTotal   ! From ModPlanet, hopefully
end function get_nSpeciesTotal

integer function get_nIons()
   get_nIons = nIons   ! From ModPlanet, hopefully
end function get_nIons

integer function get_nSpeciesAll()
   get_nSpeciesAll = nSpeciesAll   ! From ModPlanet, hopefully
end function get_nSpeciesAll


subroutine decode_gitm_indices( varname, gitm_varname, gitm_dim, gitm_index, &
                                long_name, units)
! The rosetta stone relating the user input 'strings' to integer indices. 
!
! progvar%varname      = varname
! progvar%long_name    = long_name
! progvar%units        = units
! progvar%gitm_varname = gitm_varname
! progvar%gitm_dim     = gitm_dim
! progvar%gitm_index   = gitm_index

character(len=*),              intent(in)  :: varname
character(len=*),              intent(out) :: gitm_varname
integer,                       intent(out) :: gitm_dim, gitm_index
character(len=NF90_MAX_NAME),  intent(out) :: long_name
character(len=NF90_MAX_NAME),  intent(out) :: units



   long_name    = 'something real'
   units        = 'furlongs/fortnight'

   select case (trim(varname))

   ! The first hunk of these all come from the NDensityS variable, defined to be:
   ! do iSpecies=1,nSpeciesTotal
   !    write(iRestartUnit_) NDensityS(:,:,:,iSpecies,iBlock)
   ! enddo

   case ('iO_3P_NDensityS') 
      gitm_varname = 'NDensityS'
      gitm_dim     = 4
      gitm_index   = iO_3P_
      long_name    = 'density of O3P molecules'
      units        = 'mol/m3'

   case ('iO2_NDensityS')
      gitm_varname = 'NDensityS'
      gitm_dim     = 4
      gitm_index = iO2_
      long_name    = 'density of O2 molecules'
      units        = 'mol/m3'

   case ('iN2_NDensityS')
      gitm_varname = 'NDensityS'
      gitm_dim     = 4
      gitm_index = iN2_
      long_name    = 'density of N2 molecules'
      units        = 'mol/m3'

   case ('iN_4S_NDensityS')
      gitm_varname = 'NDensityS'
      gitm_dim     = 4
      gitm_index = iN_4S_
      long_name    = 'density of N4S molecules'
      units        = 'mol/m3'

   case ('iNO_NDensityS')
      gitm_varname = 'NDensityS'
      gitm_dim     = 4
      gitm_index = iNO_
      long_name    = 'density of NO molecules'
      units        = 'mol/m3'

   case ('iN_2D_NDensityS')
      gitm_varname = 'NDensityS'
      gitm_dim     = 4
      gitm_index = iN_2D_
      long_name    = 'density of N2D molecules'
      units        = 'mol/m3'

   case ('iN_2P_NDensityS')
      gitm_varname = 'NDensityS'
      gitm_dim     = 4
      gitm_index = iN_2P_
      long_name    = 'density of N2P molecules'
      units        = 'mol/m3'

   case ('iH_NDensityS')
      gitm_varname = 'NDensityS'
      gitm_dim     = 4
      gitm_index = iH_
      long_name    = 'density of H molecules'
      units        = 'mol/m3'

   case ('iHe_NDensityS')
      gitm_varname = 'NDensityS'
      gitm_dim     = 4
      gitm_index = iHe_
      long_name    = 'density of He molecules'
      units        = 'mol/m3'

   case ('iCO2_NDensityS')
      gitm_varname = 'NDensityS'
      gitm_dim     = 4
      gitm_index = iCO2_
      long_name    = 'density of CO2 molecules'
      units        = 'mol/m3'

   case ('iO_1D_NDensityS')
      gitm_varname = 'NDensityS'
      gitm_dim     = 4
      gitm_index = iO_1D_
      long_name    = 'density of O1D molecules'
      units        = 'mol/m3'

   ! The next hunk of these all pertain to the IDensityS variable:
   ! do iSpecies=1,nIons
   !    write(iRestartUnit_) IDensityS(:,:,:,iSpecies,iBlock)
   ! enddo

   case ('iO_4SP_IDensityS')
      gitm_varname = 'IDensityS'
      gitm_dim     = 4
      gitm_index   = iO_4SP_
      long_name    = 'density of O4SP ions'
      units        = 'mol/m3'

   case ('iO2P_IDensityS')
      gitm_varname = 'IDensityS'
      gitm_dim     = 4
      gitm_index   = iO2P_
      long_name    = 'density of O2P ions'
      units        = 'mol/m3'

   case ('iN2P_IDensityS')
      gitm_varname = 'IDensityS'
      gitm_dim     = 4
      gitm_index   = iN2P_
      long_name    = 'density of N2P ions'
      units        = 'mol/m3'

   case ('iNP_IDensityS')
      gitm_varname = 'IDensityS'
      gitm_dim     = 4
      gitm_index   = iNP_
      long_name    = 'density of NP ions'
      units        = 'mol/m3'

   case ('iNOP_IDensityS')
      gitm_varname = 'IDensityS'
      gitm_dim     = 4
      gitm_index   = iNOP_
      long_name    = 'density of NOP ions'
      units        = 'mol/m3'

   case ('iO_2DP_IDensityS')
      gitm_varname = 'IDensityS'
      gitm_dim     = 4
      gitm_index   = iO_2DP_
      long_name    = 'density of O2DP ions'
      units        = 'mol/m3'

   case ('iO_2PP_IDensityS')
      gitm_varname = 'IDensityS'
      gitm_dim     = 4
      gitm_index   = iO_2PP_
      long_name    = 'density of O2PP ions'
      units        = 'mol/m3'

   case ('iHP_IDensityS')
      gitm_varname = 'IDensityS'
      gitm_dim     = 4
      gitm_index   = iHP_
      long_name    = 'density of HP ions'
      units        = 'mol/m3'

   case ('iHeP_IDensityS')
      gitm_varname = 'IDensityS'
      gitm_dim     = 4
      gitm_index   = iHeP_
      long_name    = 'density of HeP ions'
      units        = 'mol/m3'

   case ('ie_IDensityS')
      gitm_varname = 'IDensityS'
      gitm_dim     = 4
      gitm_index   = ie_
      long_name    = 'density of the electrons'
      units        = 'mol/m3'

   case ('Temperature') ! write(iRestartUnit_) Temperature(:,:,:,iBlock)*TempUnit(:,:,:)
      gitm_varname = 'Temperature'
      gitm_dim     = -1
      gitm_index   = -1
      long_name    = 'temperature (quantity tied to the square of velocity of the particles)'
      units        = 'Kelvin'

   case ('ITemperature') ! write(iRestartUnit_) ITemperature(:,:,:,iBlock)
      gitm_varname = 'ITemperature'
      gitm_dim     = -1
      gitm_index   = -1
      long_name    = 'ion temperature (quantity tied to the square of velocity of the ions)'
      units        = 'Kelvin'

   case ('eTemperature') ! write(iRestartUnit_) eTemperature(:,:,:,iBlock)
      gitm_varname = 'eTemperature'
      gitm_dim     = -1
      gitm_index   = -1
      long_name    = 'electron temperature (quantity tied to the square of velocity of the electrons)'
      units        = 'Kelvin'

   case ('U_Velocity_component') ! write(iRestartUnit_) Velocity(:,:,:,iBlock)
      gitm_varname = 'Velocity'
      gitm_dim     = 4
      gitm_index   = 1
      long_name    = 'the U-component of the velocity of the particles' 
      units        = 'm/s'

   case ('V_Velocity_component') ! write(iRestartUnit_) Velocity(:,:,:,iBlock)
      gitm_varname = 'Velocity'
      gitm_dim     = 4
      gitm_index   = 2
      long_name    = 'the V-component of the velocity of the particles' 
      units        = 'm/s'

   case ('W_Velocity_component') ! write(iRestartUnit_) Velocity(:,:,:,iBlock)
      gitm_varname = 'Velocity'
      gitm_dim     = 4
      gitm_index   = 3
      long_name    = 'the W-component of the velocity of the particles' 
      units        = 'm/s'

   case ('U_IVelocity_component') ! write(iRestartUnit_) Velocity(:,:,:,iBlock)
      gitm_varname = 'IVelocity'
      gitm_dim     = 4
      gitm_index   = 1
      long_name    = 'the U-component of the velocity of the ions' 
      units        = 'm/s'

   case ('V_IVelocity_component') ! write(iRestartUnit_) Velocity(:,:,:,iBlock)
      gitm_varname = 'IVelocity'
      gitm_dim     = 4
      gitm_index   = 2
      long_name    = 'the V-component of the velocity of the ions' 
      units        = 'm/s'

   case ('W_IVelocity_component') ! write(iRestartUnit_) IVelocity(:,:,:,iBlock)
      gitm_varname = 'IVelocity'
      gitm_dim     = 4
      gitm_index   = 3
      long_name    = 'the W-component of the velocity of the ions' 
      units        = 'm/s'

   case ('iO_3P_VerticalVelocity')
      gitm_varname = 'VerticalVelocity'
      gitm_dim     = 4
      gitm_index   = iO_3P_
      long_name    = 'the vertical velocity of the O3P molecule' 
      units        = 'm/s'

   case ('iO2_VerticalVelocity')
      gitm_varname = 'VerticalVelocity'
      gitm_dim     = 4
      gitm_index   = iO2_
      long_name    = 'the vertical velocity of the O2 molecule' 
      units        = 'm/s'

   case ('iN2_VerticalVelocity')
      gitm_varname = 'VerticalVelocity'
      gitm_dim     = 4
      gitm_index   = iN2_
      long_name    = 'the vertical velocity of the N2 molecule' 
      units        = 'm/s'

   case ('iN_4S_VerticalVelocity')
      gitm_varname = 'VerticalVelocity'
      gitm_dim     = 4
      gitm_index   = iN_4S_
      long_name    = 'the vertical velocity of the N4S molecule' 
      units        = 'm/s'

   case ('iNO_VerticalVelocity')
      gitm_varname = 'VerticalVelocity'
      gitm_dim     = 4
      gitm_index   = iNO_
      long_name    = 'the vertical velocity of the NO molecule' 
      units        = 'm/s'

   case ('iHe_VerticalVelocity')
      gitm_varname = 'VerticalVelocity'
      gitm_dim     = 4
      gitm_index   = iHE_
      long_name    = 'the vertical velocity of the He molecule' 
      units        = 'm/s'

   case ('TEC')
      gitm_varname = 'TEC'
      gitm_dim     = -1
      gitm_index   = -1
      long_name    = 'Vertically integrated total electron content'
      units        = '10^16 electron/m^2'

   case ('f107') ! write(iRestartUnit_) f107_est !Alex !Does DART assume that anything that has gitm_dim = -1 is 3D? 
      gitm_varname = 'f107'
      gitm_dim     = -1
      gitm_index   = -1
      long_name    = 'f107 solar flux index'
      units        = '1 Solar Flux Unit 10^-22 Wa m^-2 Hz^-1'

   case ('Rho')
      gitm_varname = 'Rho'
      gitm_dim     = -1
      gitm_index   = -1
      long_name    = 'mass density'
      units        = 'kg/m3'

   case default

      write(string1,*)'unknown GITM variable '//trim(varname)
      call error_handler(E_ERR,'define_var_dims',string1,source,revision,revdate)

   end select


end subroutine decode_gitm_indices




!===================================================================
! End of dart_gitm_mod
!===================================================================
end module dart_gitm_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
