! DART software - Copyright UCAR. This open source software is provided
! by ucar, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/dares/dart/dart_download
! 
! $Id$
!----------------------------------------------------------------
!>
!> this routines supports unit conversions between chemistry values
!> in the cam state and the observed quantities
!>
!----------------------------------------------------------------

module chem_tables_mod


use types_mod
use utilities_mod
use obs_kind_mod

implicit none
private

public :: init_chem_tables,     &
          finalize_chem_tables, &
          get_molar_mass,       &
          get_volume_mixing_ratio

!public :: chem_convert_factors, molar_mass_dry_air

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = &
   "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"


type chem_convert
 character(len=31) :: netcdf_varname
 real(r8)          :: convert_factor
 character(len=31) :: quantity_name
end type


! table of chemistry conversion factors from mmr to vmr and back
type(chem_convert), allocatable :: chem_conv_table(:) 

real(r8), parameter :: molar_mass_dry_air = 28.9644_r8

integer :: num_qtys
character(len=256) :: string1, string2
logical :: module_initialized = .false.

contains 

!--------------------------------------------------------------------
!> call once to allocate and initialize the table

subroutine init_chem_tables()

! to add more quantities, add entries here following the pattern.
! there can be only one entry per 'QTY_xxx'

integer :: i

if (module_initialized) return

num_qtys = get_num_quantities()

allocate(chem_conv_table(0:num_qtys))

! initialize all entries to an empty name and a factor of 1.0
do i=0, num_qtys
   call set_entry('', 1.0_r8, i)
enddo

! and now add entries for real items
call add_entry('H',     1.0074_r8, 'QTY_ATOMIC_H_MIXING_RATIO')
call add_entry('O',    15.9994_r8, 'QTY_ATOMIC_OXYGEN_MIXING_RATIO')
call add_entry('O2',   31.9988_r8, 'QTY_MOLEC_OXYGEN_MIXING_RATIO')
call add_entry('O3',   47.9982_r8, 'QTY_O3') 
call add_entry('N2',   28.0135_r8, 'QTY_NITROGEN')


!%!  'N2O                            ',       44.01288,     'QTY_N2O                         ', &
!%!  'NO                             ',       30.00614,     'QTY_NO                          ', &
!%!  'NO2                            ',       46.00554,     'QTY_NO2                         ', &
!%!  'NO3                            ',       62.00494,     'QTY_NO3                         ', &
!%!  'HNO3                           ',       63.01234,     'QTY_HNO3                        ', &
!%!  'HO2NO2                         ',       79.01174,     'QTY_HO2NO2                      ', &
!%!  'N2O5                           ',      108.01048,     'QTY_N2O5                        ', &
!%!  'H2                             ',         2.0148,     'QTY_H2                          ', &
!%!  'OH                             ',        17.0068,     'QTY_OH                          ', &
!%!  'HO2                            ',        33.0062,     'QTY_HO2                         ', &
!%!  'H2O2                           ',        34.0136,     'QTY_H2O2                        ', &
!%!  'CH4                            ',        16.0406,     'QTY_CH4                         ', &
!%!  'CO                             ',        28.0104,     'QTY_CO                          ', &
!%!  'CH3O2                          ',         47.032,     'QTY_CH3O2                       ', &
!%!  'CH3OOH                         ',        48.0394,     'QTY_CH3OOH                      ', &
!%!  'CH2O                           ',        30.0252,     'QTY_CH2O                        ', &
!%!  'CH3OH                          ',          32.04,     'QTY_CH3OH                       ', &
!%!  'C2H5OH                         ',        46.0658,     'QTY_C2H5OH                      ', &
!%!  'C2H4                           ',        28.0516,     'QTY_C2H4                        ', &
!%!  'EO                             ',        61.0578,     'QTY_EO                          ', &
!%!  'EO2                            ',        77.0572,     'QTY_EO2                         ', &
!%!  'CH3COOH                        ',        60.0504,     'QTY_CH3COOH                     ', &
!%!  'GLYALD                         ',        60.0504,     'QTY_GLYALD                      ', &
!%!  'C2H6                           ',        30.0664,     'QTY_C2H6                        ', &
!%!  'C2H5O2                         ',        61.0578,     'QTY_C2H5O2                      ', &
!%!  'C2H5OOH                        ',        62.0652,     'QTY_C2H5OOH                     ', &
!%!  'CH3CHO                         ',         44.051,     'QTY_CH3CHO                      ', &
!%!  'CH3CO3                         ',        75.0424,     'QTY_CH3CO3                      ', &
!%!  'CH3COOOH                       ',        76.0498,     'QTY_CH3COOOH                    ', &
!%!  'C3H6                           ',        42.0774,     'QTY_C3H6                        ', &
!%!  'C3H8                           ',        44.0922,     'QTY_C3H8                        ', &
!%!  'C3H7O2                         ',        75.0836,     'QTY_C3H7O2                      ', &
!%!  'C3H7OOH                        ',         76.091,     'QTY_C3H7OOH                     ', &
!%!  'PO2                            ',         91.083,     'QTY_PO2                         ', &
!%!  'POOH                           ',        92.0904,     'QTY_POOH                        ', &
!%!  'CH3COCH3                       ',        58.0768,     'QTY_CH3COCH3                    ', &
!%!  'RO2                            ',        89.0682,     'QTY_RO2                         ', &
!%!  'ROOH                           ',        90.0756,     'QTY_ROOH                        ', &
!%!  'BIGENE                         ',        56.1032,     'QTY_BIGENE                      ', &
!%!  'ENEO2                          ',       105.1088,     'QTY_ENEO2                       ', &
!%!  'MEK                            ',        72.1026,     'QTY_MEK                         ', &
!%!  'MEKO2                          ',        103.094,     'QTY_MEKO2                       ', &
!%!  'MEKOOH                         ',       104.1014,     'QTY_MEKOOH                      ', &
!%!  'BIGALK                         ',        72.1438,     'QTY_BIGALK                      ', &
!%!  'ALKO2                          ',       103.1352,     'QTY_ALKO2                       ', &
!%!  'ALKOOH                         ',       104.1426,     'QTY_ALKOOH                      ', &
!%!  'ISOP                           ',        68.1142,     'QTY_ISOP                        ', &
!%!  'ISOPO2                         ',       117.1198,     'QTY_ISOPO2                      ', &
!%!  'ISOPOOH                        ',       118.1272,     'QTY_ISOPOOH                     ', &
!%!  'MVK                            ',        70.0878,     'QTY_MVK                         ', &
!%!  'MACR                           ',        70.0878,     'QTY_MACR                        ', &
!%!  'MACRO2                         ',       119.0934,     'QTY_MACRO2                      ', &
!%!  'MACROOH                        ',       120.1008,     'QTY_MACROOH                     ', &
!%!  'MCO3                           ',       101.0792,     'QTY_MCO3                        ', &
!%!  'HYDRALD                        ',        100.113,     'QTY_HYDRALD                     ', &
!%!  'HYAC                           ',        74.0762,     'QTY_HYAC                        ', &
!%!  'CH3COCHO                       ',        72.0614,     'QTY_CH3COCHO                    ', &
!%!  'XO2                            ',       149.1186,     'QTY_XO2                         ', &
!%!  'XOOH                           ',        150.126,     'QTY_XOOH                        ', &
!%!  'C10H16                         ',       136.2284,     'QTY_C10H16                      ', &
!%!  'TERPO2                         ',        185.234,     'QTY_TERPO2                      ', &
!%!  'TERPOOH                        ',       186.2414,     'QTY_TERPOOH                     ', &
!%!  'TOLUENE                        ',        92.1362,     'QTY_TOLUENE                     ', &
!%!  'CRESOL                         ',       108.1356,     'QTY_CRESOL                      ', &
!%!  'TOLO2                          ',       173.1406,     'QTY_TOLO2                       ', &
!%!  'TOLOOH                         ',        174.148,     'QTY_TOLOOH                      ', &
!%!  'XOH                            ',       190.1474,     'QTY_XOH                         ', &
!%!  'BIGALD                         ',        98.0982,     'QTY_BIGALD                      ', &
!%!  'GLYOXAL                        ',        58.0356,     'QTY_GLYOXAL                     ', &
!%!  'PAN                            ',      121.04794,     'QTY_PAN                         ', &
!%!  'ONIT                           ',      119.07434,     'QTY_ONIT                        ', &
!%!  'MPAN                           ',      147.08474,     'QTY_MPAN                        ', &
!%!  'ISOPNO3                        ',      162.11794,     'QTY_ISOPNO3                     ', &
!%!  'ONITR                          ',      147.12594,     'QTY_ONITR                       ', &
!%!  'SOA                            ',        144.132,     'QTY_SOA                         ', &
!%!  'SO2                            ',        64.0648,     'QTY_SO2                         ', &
!%!  'DMS                            ',        62.1324,     'QTY_DMS                         ', &
!%!  'NH3                            ',       17.02894,     'QTY_NH3                         ', &
!%!  'NH4                            ',       18.03634,     'QTY_NH4                         ', &
!%!  'NH4NO3                         ',       80.04128,     'QTY_NH4NO3                      ', &
!%!  'Rn                             ',          222.0,     'QTY_Rn                          ', &
!%!  'Pb                             ',          207.2,     'QTY_Pb                          ', &
!%!  'HCN                            ',       27.02514,     'QTY_HCN                         ', &
!%!  'CH3CN                          ',       41.05094,     'QTY_CH3CN                       ', &
!%!  'C2H2                           ',        26.0368,     'QTY_C2H2                        ', &
!%!  'HCOOH                          ',        46.0246,     'QTY_HCOOH                       ', &
!%!  'HOCH2OO                        ',        63.0314,     'QTY_HOCH2OO                     ', &
!%!  'H2SO4                          ',        98.0784,     'QTY_H2SO4                       ', &
!%!  'SOAG                           ',         12.011,     'QTY_SOAG                        ', &
!%!  'so4_a1                         ',      115.10734,     'QTY_so4_a1                      ', &
!%!  'pom_a1                         ',         12.011,     'QTY_pom_a1                      ', &
!%!  'soa_a1                         ',         12.011,     'QTY_soa_a1                      ', &
!%!  'bc_a1                          ',         12.011,     'QTY_bc_a1                       ', &
!%!  'dst_a1                         ',     135.064039,     'QTY_dst_a1                      ', &
!%!  'ncl_a1                         ',      58.442468,     'QTY_ncl_a1                      ', &
!%!  'num_a1                         ',         1.0074,     'QTY_num_a1                      ', &
!%!  'so4_a2                         ',      115.10734,     'QTY_so4_a2                      ', &
!%!  'soa_a2                         ',         12.011,     'QTY_soa_a2                      ', &
!%!  'ncl_a2                         ',      58.442468,     'QTY_ncl_a2                      ', &
!%!  'num_a2                         ',         1.0074,     'QTY_num_a2                      ', &
!%!  'dst_a3                         ',     135.064039,     'QTY_dst_a3                      ', &
!%!  'ncl_a3                         ',      58.442468,     'QTY_ncl_a3                      ', &
!%!  'so4_a3                         ',      115.10734,     'QTY_so4_a3                      ', &
!%!  'num_a3                         ',         1.0074,     'QTY_num_a3                      ', &
!%!  'CO01                           ',        28.0104,     'QTY_CO01                        ', &
!%!  'CO02                           ',        28.0104,     'QTY_CO02                        ', &
!%!  'CO03                           ',        28.0104,     'QTY_CO03                        ', &
!%!  'CO04                           ',        28.0104,     'QTY_CO04                        ', &
!%!  'CO05                           ',        28.0104,     'QTY_CO05                        ', &
!%!  'CO06                           ',        28.0104,     'QTY_CO06                        ', &
!%!  'CO07                           ',        28.0104,     'QTY_CO07                        ', &
!%!  'CO08                           ',        28.0104,     'QTY_CO08                        ', &
!%!  'CO09                           ',        28.0104,     'QTY_CO09                        ', &
!%!  'CB1                            ',         12.011,     'QTY_CB1                         ', &
!%!  'CB2                            ',         12.011,     'QTY_CB2                         ', &
!%!  'OC1                            ',         12.011,     'QTY_OC1                         ', &
!%!  'OC2                            ',         12.011,     'QTY_OC2                         ', &
!%!  'CB101                          ',         12.011,     'QTY_CB101                       ', &
!%!  'CB201                          ',         12.011,     'QTY_CB201                       ', &
!%!  'OC101                          ',         12.011,     'QTY_OC101                       ', &
!%!  'OC201                          ',         12.011,     'QTY_OC201                       ', &
!%!  'CB102                          ',         12.011,     'QTY_CB102                       ', &
!%!  'CB202                          ',         12.011,     'QTY_CB202                       ', &
!%!  'OC102                          ',         12.011,     'QTY_OC102                       ', &
!%!  'OC202                          ',         12.011,     'QTY_OC202                       '    


end subroutine init_chem_tables

!--------------------------------------------------------------------
!> call at end to free the table space

subroutine finalize_chem_tables

deallocate(chem_conv_table)

end subroutine finalize_chem_tables

!--------------------------------------------------------------------
!>

subroutine add_entry(netcdf_varname, convert_factor, quantity_name)

character(len=*), intent(in) :: netcdf_varname
real(r8),         intent(in) :: convert_factor
character(len=*), intent(in) :: quantity_name

integer :: qty_index

! get qty indx, error if not found
qty_index = get_index_for_quantity(quantity_name)
if (qty_index < 0) then
   write(string1,'(3A)') 'quantity string "', trim(quantity_name), &
                         '" not found in known quantities list'
   write(string2, *) 'check obs_kind_mod.f90 for valid quantities; defined by preprocess namelist'
   call error_handler(E_ERR, 'chem_convert_factor', string1, &
                      source, revision, revdate, text2=string2)
endif

! build type, add to array at indx

chem_conv_table(qty_index) = chem_convert(netcdf_varname, convert_factor, quantity_name)

end subroutine add_entry

!--------------------------------------------------------------------
!>

subroutine set_entry(netcdf_varname, convert_factor, quantity_index)

character(len=*), intent(in) :: netcdf_varname
real(r8),         intent(in) :: convert_factor
integer,          intent(in) :: quantity_index

character(len=32) :: quantity_name

quantity_name = get_name_for_quantity(quantity_index)
if (quantity_name == '') then
   write(string1,'(3A)') 'quantity index "', quantity_index, &
                         '" not found in known quantities list'
   write(string2, *) 'check obs_kind_mod.f90 for valid quantities; defined by preprocess namelist'
   call error_handler(E_ERR, 'chem_convert_factor', string1, &
                      source, revision, revdate, text2=string2)
endif

chem_conv_table(quantity_index) = chem_convert(netcdf_varname, convert_factor, quantity_name)

end subroutine set_entry

!--------------------------------------------------------------------
!>
function get_molar_mass(qty)

integer, intent(in) :: qty
real(r8)            :: get_molar_mass

if (qty < 0 .or. qty > num_qtys) then
   write(string1,'(A,I6,A,I6)') 'quantity number ', qty, &
                                   ' must be between 0 and ', num_qtys
   call error_handler(E_ERR, 'get_molar_mass', string1, &
                      source, revision, revdate)
endif

get_molar_mass = chem_conv_table(qty)%convert_factor

end function get_molar_mass

!--------------------------------------------------------------------
!>
function get_volume_mixing_ratio(qty)

integer, intent(in) :: qty
real(r8)            :: get_volume_mixing_ratio

if (qty < 0 .or. qty > num_qtys) then
   write(string1,'(A,I6,A,I6)') 'quantity number ', qty, &
                                   ' must be between 0 and ', num_qtys
   call error_handler(E_ERR, 'get_volume_mixing_ratio', string1, &
                      source, revision, revdate)
endif
if (chem_conv_table(qty)%convert_factor /= 1.0_r8) then
   get_volume_mixing_ratio = molar_mass_dry_air / chem_conv_table(qty)%convert_factor
else
   get_volume_mixing_ratio  = 1.0_r8
endif

end function get_volume_mixing_ratio

!--------------------------------------------------------------------

end module chem_tables_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$


