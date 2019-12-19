!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModConstants

!-------------------------------------------------------------------
! Useful physical constants
!-------------------------------------------------------------------

  real, parameter :: Avogadros_Number       = 6.02e23               ! mol./mole
  real, parameter :: AMU                    = 1.6726e-27            ! Kg
  real, parameter :: Mass_Ion               = 1.6726e-27            ! Kg
  real, parameter :: Mass_Proton            = Mass_Ion              ! Kg
  real, parameter :: Mass_Electron          = 9.1094e-31            ! Kg
  real, parameter :: Boltzmanns_Constant    = 1.3807e-23            ! J/K
  real, parameter :: Planck_Constant        = 6.6261e-34            ! Js
  real, parameter :: Element_Charge         = 1.6022e-19            ! C (J/eV)
  real, parameter :: Speed_Light            = 2.9979e8              ! m/s
  real, parameter :: Univ_Gas_Constant      = Avogadros_Number*   &
                                              Boltzmanns_Constant   ! J/(moleK)
  real, parameter :: RGAS                   = Univ_Gas_Constant*   &
                                              1.0E+07               !erg/(moleK)
  real, parameter :: pi                     = 3.141592653589793
  real, parameter :: twopi = 2*pi

  real, parameter :: P00                    = 1.0e5                   ! Pa
  real, parameter :: T00                    = 273.0                   ! K

  real, parameter :: SpecificHeatVolume     = 3./2.
  real, parameter :: SpecificHeatPressure   = SpecificHeatVolume + 1.0
  real, parameter :: Gamma_const                  = &
       SpecificHeatPressure/SpecificHeatVolume

end module ModConstants
