! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
MODULE wrf_data_module

! <next four lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$ 

use        types_mod, only : r8

implicit none
private

public wrf_data

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

TYPE wrf_data

   integer :: ncid  ! netcdf id for file
   integer :: bt_id, bt, sn_id, sn, we_id, we
   integer :: u_id, v_id, w_id, ph_id, phb_id, t_id,   &
                        mu_id, mub_id,                           &
                        qv_id, qc_id, qr_id, qi_id, qs_id, qg_id 
   integer :: ptop_id
   logical :: ice_micro


!---
!  arrays for data

   real(r8), pointer :: u(:,:,:)
   real(r8), pointer :: v(:,:,:)
   real(r8), pointer :: w(:,:,:)
   real(r8), pointer :: ph(:,:,:)
   real(r8), pointer :: phb(:,:,:)
   real(r8), pointer :: t(:,:,:)
   real(r8), pointer :: qv(:,:,:)
   real(r8), pointer :: qc(:,:,:)
   real(r8), pointer :: qr(:,:,:)
   real(r8), pointer :: qi(:,:,:)
   real(r8), pointer :: qs(:,:,:)
   real(r8), pointer :: qg(:,:,:)
   real(r8), pointer :: mu(:,:)
   real(r8), pointer :: mub(:,:)

end type
END MODULE wrf_data_module
