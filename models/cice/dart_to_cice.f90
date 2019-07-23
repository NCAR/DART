! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program dart_to_cice

!----------------------------------------------------------------------
! purpose: implement a 'partition function' to modify the cice state 
!          to be consistent with the states from assimilation
!
! method: Read in restart (restart with prior) and out restart (restart 
!         with posterior) written by DART after filter. 
!
! author: C Bitz June 2016
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             file_exist, error_handler, E_ERR, E_MSG, to_upper
use  netcdf_utilities_mod, only : nc_check
use netcdf 

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = &
   "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"

!------------------------------------------------------------------

character(len=256) :: dart_to_cice_input_file = 'dart_restart.nc'
character(len=256) :: original_cice_input_file = 'cice_restart.nc'
character(len=256) :: previous_cice_input_file = 'pre_restart.nc'
character(len=128) :: balance_method = 'simple_squeeze'
character(len=15)  :: r_snw_name  = 'r_snw'

namelist /dart_to_cice_nml/ dart_to_cice_input_file, & 
                            original_cice_input_file, &
                            previous_cice_input_file, &
                            balance_method, &
                            r_snw_name

character(len=512) :: string1, string2, msgstring
character(len=15)  :: varname
character(len=128) :: method

integer :: Nx, Ny
integer :: Ncat   ! number of categories in ice-thickness dist
integer, parameter :: Nilyr = 8   ! number of layers in ice, hardwired
integer, parameter :: Nslyr = 8   ! number of layers in snow, hardwired

real(r8), allocatable :: aicen_original(:,:,:)
real(r8), allocatable :: vicen_original(:,:,:)
real(r8), allocatable :: vsnon_original(:,:,:)
real(r8), allocatable :: aice_original(:,:)
!real(r8), allocatable :: vice_original(:,:)
!real(r8), allocatable :: vsno_original(:,:)
real(r8), allocatable :: hicen_original(:,:,:)
real(r8), allocatable :: hsnon_original(:,:,:)

real(r8), allocatable :: aicen_pre(:,:,:)
real(r8), allocatable :: vicen_pre(:,:,:)
real(r8), allocatable :: vsnon_pre(:,:,:)
real(r8), allocatable :: aice_pre(:,:)
!real(r8), allocatable :: vice_pre(:,:)
!real(r8), allocatable :: vsno_pre(:,:)

real(r8), allocatable :: aicen(:,:,:)
real(r8), allocatable :: vicen(:,:,:)
real(r8), allocatable :: vsnon(:,:,:)
real(r8), allocatable :: Tsfcn(:,:,:)
real(r8), allocatable :: sice001(:,:,:)
real(r8), allocatable :: sice002(:,:,:)
real(r8), allocatable :: sice003(:,:,:)
real(r8), allocatable :: sice004(:,:,:)
real(r8), allocatable :: sice005(:,:,:)
real(r8), allocatable :: sice006(:,:,:)
real(r8), allocatable :: sice007(:,:,:)
real(r8), allocatable :: sice008(:,:,:)
real(r8), allocatable :: qice001(:,:,:)
real(r8), allocatable :: qice002(:,:,:)
real(r8), allocatable :: qice003(:,:,:)
real(r8), allocatable :: qice004(:,:,:)
real(r8), allocatable :: qice005(:,:,:)
real(r8), allocatable :: qice006(:,:,:)
real(r8), allocatable :: qice007(:,:,:)
real(r8), allocatable :: qice008(:,:,:)
real(r8), allocatable :: qsno001(:,:,:)
real(r8), allocatable :: qsno002(:,:,:)
real(r8), allocatable :: qsno003(:,:,:)
real(r8), allocatable :: aice(:,:)
!real(r8), allocatable :: vice(:,:)
!real(r8), allocatable :: vsno(:,:)

!Parameters
real(r8), allocatable :: r_snw(:,:)

!Temporary variables
real(r8), allocatable :: aice_temp(:,:)
real(r8), allocatable :: increment_aice(:,:)
!real(r8), allocatable :: increment_vice(:,:)
!real(r8), allocatable :: increment_vsno(:,:)

real(r8), allocatable :: tendency_aice(:,:)
!real(r8), allocatable :: tendency_vice(:,:)
!real(r8), allocatable :: tendency_vsno(:,:)
real(r8), allocatable :: tendency_aicen(:,:,:)
!real(r8), allocatable :: tendency_vicen(:,:,:)
!real(r8), allocatable :: tendency_vsnon(:,:,:)

real(r8) :: R, weight_aicen     !, weight_vicen, weight_vsnon

integer  :: i, j, n
! integer :: k
integer  :: VarID, ncid,ncid2, iunit, io, ndims
real(r8) :: squeeze

real(r8), parameter :: &        !from ice_shortwave.F90
     rsnw_max   = 1.6_r8, &   
     rsnw_min   = -2.0_r8

real(r8), parameter :: &      ! from ice_therm_vertical.F90
     phi_init = 0.75_r8, &    ! initial liquid fraction of frazil ice 
     dSin0_frazil = 3.0_r8    ! bulk salinity reduction of newly formed frazil 

real(r8), parameter :: &      ! from pop constants.F90
     sss = 34.7_r8  !  use ocn_ref_salinity since we dont know truth

real(r8), parameter :: &
     c1  = 1.0_r8

real(r8), parameter :: Tsmelt = 0._r8 ! from ice_constant

real(r8) :: cc1, cc2, cc3, x1, Si0new, qi0new, Ti
real(r8), allocatable ::  hin_max(:)
real(r8), allocatable ::  hcat_midpoint(:)

!----------------------------------------------------------------------

call initialize_utilities(progname='dart_to_cice')

call find_namelist_in_file("input.nml", "dart_to_cice_nml", iunit)
read(iunit, nml = dart_to_cice_nml, iostat = io)
call check_namelist_read(iunit, io, "dart_to_cice_nml")

method = balance_method
call to_upper(method)

! check on namelist stuff, and whether files exist
write(string1,*) 'converting DART output file "'// &
                 &trim(dart_to_cice_input_file)//'" to one CICE will like'
write(string2,*) 'using the "'//trim(balance_method)//'" method.'
call error_handler(E_MSG,'dart_to_cice',string1,text2=string2)

if ( .not. file_exist(dart_to_cice_input_file) ) then
   write(string1,*) 'cannot open "', trim(dart_to_cice_input_file),'" for updating.'
   call error_handler(E_ERR,'dart_to_cice:filename not found ',trim(dart_to_cice_input_file))
endif

if ( .not. file_exist(original_cice_input_file) ) then
   write(string1,*) 'cannot open "', trim(original_cice_input_file),'" for reading.'
   call error_handler(E_ERR,'dart_to_cice:filename not found ',trim(original_cice_input_file))
endif

! open original restart file with read only
call nc_check( nf90_open(trim(original_cice_input_file), NF90_NOWRITE, ncid), &
                  'dart_to_cice', 'open "'//trim(original_cice_input_file)//'"')

! get the original ice concentration, ice volume and snow volume (FYI it is allocated in routine)
call get_3d_variable(ncid, 'aicen', aicen_original, original_cice_input_file)
call get_3d_variable(ncid, 'vicen', vicen_original, original_cice_input_file)
call get_3d_variable(ncid, 'vsnon', vsnon_original, original_cice_input_file)

call nc_check(nf90_close(ncid),'dart_to_cice', 'close '//trim(original_cice_input_file))

! open posterior restart file with read and write 
call nc_check( nf90_open(trim(dart_to_cice_input_file), NF90_WRITE, ncid), &
                  'dart_to_cice', 'open "'//trim(dart_to_cice_input_file)//'"')

! get the key restart variables (FYI allocated in the routine)
call get_3d_variable(ncid, 'aicen', aicen, dart_to_cice_input_file)
call get_3d_variable(ncid, 'vicen', vicen, dart_to_cice_input_file)
call get_3d_variable(ncid, 'vsnon', vsnon, dart_to_cice_input_file)
call get_3d_variable(ncid, 'Tsfcn', Tsfcn, dart_to_cice_input_file)
call get_3d_variable(ncid, 'sice001', sice001, dart_to_cice_input_file)
call get_3d_variable(ncid, 'sice002', sice002, dart_to_cice_input_file)
call get_3d_variable(ncid, 'sice003', sice003, dart_to_cice_input_file)
call get_3d_variable(ncid, 'sice004', sice004, dart_to_cice_input_file)
call get_3d_variable(ncid, 'sice005', sice005, dart_to_cice_input_file)
call get_3d_variable(ncid, 'sice006', sice006, dart_to_cice_input_file)
call get_3d_variable(ncid, 'sice007', sice007, dart_to_cice_input_file)
call get_3d_variable(ncid, 'sice008', sice008, dart_to_cice_input_file)
call get_3d_variable(ncid, 'qice001', qice001, dart_to_cice_input_file)
call get_3d_variable(ncid, 'qice002', qice002, dart_to_cice_input_file)
call get_3d_variable(ncid, 'qice003', qice003, dart_to_cice_input_file)
call get_3d_variable(ncid, 'qice004', qice004, dart_to_cice_input_file)
call get_3d_variable(ncid, 'qice005', qice005, dart_to_cice_input_file)
call get_3d_variable(ncid, 'qice006', qice006, dart_to_cice_input_file)
call get_3d_variable(ncid, 'qice007', qice007, dart_to_cice_input_file)
call get_3d_variable(ncid, 'qice008', qice008, dart_to_cice_input_file)
call get_3d_variable(ncid, 'qsno001', qsno001, dart_to_cice_input_file)
call get_3d_variable(ncid, 'qsno002', qsno002, dart_to_cice_input_file)
call get_3d_variable(ncid, 'qsno003', qsno003, dart_to_cice_input_file)

! get the parameter variables in the restart
call get_2d_variable(ncid, r_snw_name, r_snw, dart_to_cice_input_file)

Nx   = size(aicen,1)
Ny   = size(aicen,2)
Ncat = size(aicen,3)

allocate(aice(Nx,Ny))
allocate(aice_temp(Nx,Ny))
allocate(hicen_original(Nx,Ny,Ncat))
allocate(hsnon_original(Nx,Ny,Ncat))

SELECT CASE (method)

  CASE ('SIMPLE_SQUEEZE')

      allocate(aice_original(Nx,Ny))
      do n = 2, Ncat
        aice_original  = aice_original + aicen_original(:,:,n)
      end do 

  CASE ('TENDENCY_WEIGHT')  
      !This method redistributes the total increment into each catetory    
      !based on each category's prior tendency weight
      !prior states:     aicen_original(Nx,Ny,Ncat), aice_original(Nx,Ny),
      !                  aicen_pre(Nx,Ny,Ncat),      aice_pre(Nx,Ny),  !previous day
      !                  vicen_original(Nx,Ny,Ncat), vice_original(Nx,Ny),
      !                  vicen_pre(Nx,Ny,Ncat),      vice_pre(Nx,Ny),  !previous day
      !                  vsnon_original(Nx,Ny,Ncat), vsno_original(Nx,Ny),
      !                  vsnon_pre(Nx,Ny,Ncat),      vsno_pre(Nx,Ny)   !previous day
      !posterior states: aicen(Nx,Ny,Ncat),aice(Nx,Ny),
      !                  vicen(Nx,Ny,Ncat),vice(Nx,Ny),
      !                  vsnon(Nx,Ny,Ncat),vsno(Nx,Ny)
      !weights:          weight_aicen,weight_vicen
      !                  weight_vsnon      
      !total increment : increment_aice(Nx,Ny), increment_vice(Nx,Ny),
      !                  increment_vsno(Nx,Ny)

      allocate (aicen_pre(Nx,Ny,Ncat))
      allocate (vicen_pre(Nx,Ny,Ncat))
      allocate (vsnon_pre(Nx,Ny,Ncat))

      allocate (tendency_aicen(Nx,Ny,Ncat))
      !allocate (tendency_vicen(Nx,Ny,Ncat))
      !allocate (tendency_vsnon(Nx,Ny,Ncat))
   
      allocate (tendency_aice(Nx,Ny))
      !allocate (tendency_vice(Nx,Ny))
      !allocate (tendency_vsno(Nx,Ny))

      allocate (aice_original(Nx,Ny))
      allocate (aice_pre(Nx,Ny))
      allocate (increment_aice(Nx,Ny))
       
      !allocate (vice(Nx,Ny))
      !allocate (vice_original(Nx,Ny))
      !allocate (vice_pre(Nx,Ny))
      !allocate (increment_vice(Nx,Ny))

      !allocate (vsno(Nx,Ny))
      !allocate (vsno_original(Nx,Ny))
      !allocate (vsno_pre(Nx,Ny))
      !allocate (increment_vsno(Nx,Ny))

      !Open the restart file from the previous day (beginning of the current
      !day's forecast)
      if ( .not. file_exist(previous_cice_input_file)) then
         write(string1,*)'cannot open "',trim(previous_cice_input_file),'" for updating.'
         call error_handler(E_ERR,'dart_to_cice',string1)
      endif

      call nc_check(nf90_open(trim(previous_cice_input_file),NF90_NOWRITE,ncid2), &
           'dart_to_cice', 'open '//trim(previous_cice_input_file))     
     
      call get_3d_variable(ncid2,'aicen',aicen_pre,previous_cice_input_file)
      call get_3d_variable(ncid2,'vicen',vicen_pre,previous_cice_input_file)     
      call get_3d_variable(ncid2,'vsnon',vsnon_pre,previous_cice_input_file)

      call nc_check(nf90_close(ncid2),'dart_to_cice','close'//trim(previous_cice_input_file))

      !calculate the weights
      aice            = aicen(:,:,1)
      aice_original   = aicen_original(:,:,1)
      aice_pre        = aicen_pre(:,:,1)
   
      !====================================================================
      !seaice/snow thickness is preserved for each category
      !hence vicen and vsnon are calculated based on their relationship
      !with aicen
      !===================================================================
      !vice            = vicen(:,:,1)
      !vice_original   = vicen_original(:,:,1)
      !vice_pre        = vicen_pre(:,:,1)

      !vsno            = vsnon(:,:,1)
      !vsno_original   = vsnon_original(:,:,1)
      !vsno_pre        = vsnon_pre(:,:,1)

      do n = 2, Ncat
         aice         = aice          + aicen(:,:,n)
         aice_original= aice_original + aicen_original(:,:,n)
         aice_pre     = aice_pre      + aicen_pre(:,:,n)
      
         !vice         = vice          + vicen(:,:,n)
         !vice_original= vice_original + vicen_original(:,:,n)
         !vice_pre     = vice_pre      + vicen_pre(:,:,n)

         !vsno         = vsno          + vsnon(:,:,n)
         !vsno_original= vsno_original + vsnon_original(:,:,n)
         !vsno_pre     = vsno_pre      + vsnon_pre(:,:,n)
      end do

      increment_aice   = aice - aice_original
      !increment_vice   = vice - vice_original
      !increment_vsno   = vsno - vsno_original

      tendency_aice    = aice_original - aice_pre
      !tendency_vice    = vice_original - vice_pre
      !tendency_vsno    = vsno_original - vsno_pre

      tendency_aicen   = aicen_original - aicen_pre
      !tendency_vicen   = vicen_original - vicen_pre
      !tendency_vsnon   = vsnon_original - vsnon_pre

      do n = 1, Ncat
         do j= 1, Ny
            do i = 1, Nx

                if (abs(increment_aice(i,j))>0._r8) then

                   R = abs(tendency_aice(i,j)/increment_aice(i,j))

                   if (R > 0.5_r8) then

                       weight_aicen = tendency_aicen(i,j,n)/tendency_aice(i,j)
                       aicen(i,j,n) = increment_aice(i,j)*weight_aicen + &
                                  aicen_original(i,j,n)    

                    else if (aice_original(i,j)>0.0_r8) then
              
                       weight_aicen = aicen_original(i,j,n)/aice_original(i,j)
                       aicen(i,j,n) = increment_aice(i,j)*weight_aicen + &
                                  aicen_original(i,j,n)

                    endif

                endif

              !  if (abs(increment_vice(i,j))>0._r8) then
                   
              !     R = abs(tendency_vice(i,j)/increment_vice(i,j))
                   
              !     if (R > 0.5_r8) then
                       
              !         weight_vicen = tendency_vicen(i,j,n)/tendency_vice(i,j)
              !         vicen(i,j,n) = increment_vice(i,j)*weight_vicen + &
              !                    vicen_original(i,j,n)
              
              !    else if (vice_original(i,j)>0.001_r8) then
                       
              !         weight_vicen = vicen_original(i,j,n)/vice_original(i,j)
              !         vicen(i,j,n) = increment_vice(i,j)*weight_vicen + &
              !                    vicen_original(i,j,n)

              !     endif
               
              !  endif

              ! if (abs(increment_vsno(i,j))>0._r8) then
                  
              !     R = abs(tendency_vsno(i,j)/increment_vsno(i,j))
                   
              !     if (R > 0.5_r8) then
                       
              !         weight_vsnon = tendency_vsnon(i,j,n)/tendency_vsno(i,j)
              !         vsnon(i,j,n) = increment_vsno(i,j)*weight_vsnon + &
              !                    vsnon_original(i,j,n)

              !     else if (vsno_original(i,j)>0.0_r8) then
                       
              !         weight_vsnon = vsnon_original(i,j,n)/vsno_original(i,j)
              !         vsnon(i,j,n) = increment_vsno(i,j)*weight_vsnon + &
              !                    vsnon_original(i,j,n)

              !     endif
                
              !  endif

            end do
         end do
      end do                       
 
   CASE DEFAULT

      write(string1,*)'input.nml:dart_to_cice_nml:balance_method "'//trim(balance_method)//'" unsupported.'
      write(string2,*)'valid values are "simple_squeeze", "tendency_weight", or "prior weight"'
      call error_handler(E_ERR,'dart_to_cice',string1, source, revision, revdate, text2=string2)

   CASE ('PRIOR_WEIGHT')  
      !Fei
      !This method redistributes the total increment into each category given the
      !prior areal weight 
      !prior states    : aicen_original(Nx,Ny,Ncat), aice_original(Nx,Ny),
      !                  vicen_original(Nx,Ny,Ncat), vice_original(Nx,Ny),
      !                  vsnon_original(Nx,Ny,Ncat), vsno_original(Nx,Ny)
      !posterior states: aicen(Nx,Ny,Ncat),aice(Nx,Ny),
      !                  vicen(Nx,Ny,Ncat),vice(Nx,Ny),
      !                  vsnon(Nx,Ny,Ncat),vsno(Nx,Ny)
      !weights         : weight_aicen(Nx,Ny,Ncat),weight_vicen(Nx,Ny,Ncat), 
      !                  weight_vsnon(Nx,Ny,Ncat)
      !total increment : increment_aice(Nx,Ny), increment_vice(Nx,Ny),
      !                  increment_vsno(Nx,Ny)
     
      allocate (aice_original(Nx,Ny))
      allocate (increment_aice(Nx,Ny))

    !  allocate (vice(Nx,Ny))
    !  allocate (vice_original(Nx,Ny))
    !  allocate (increment_vice(Nx,Ny))

    !  allocate (vsno(Nx,Ny))
    !  allocate (vsno_original(Nx,Ny))
    !  allocate (increment_vsno(Nx,Ny))

      aice          = aicen(:,:,1)
      aice_original = aicen_original(:,:,1)
      
    !  vice          = vicen(:,:,1)
    !  vice_original = vicen_original(:,:,1)

    !  vsno          = vsnon(:,:,1)
    !  vsno_original = vsnon_original(:,:,1)

      do n = 2, Ncat
        aice           = aice          + aicen(:,:,n)
        aice_original  = aice_original + aicen_original(:,:,n)    

     !   vice           = vice          + vicen(:,:,n)
     !   vice_original  = vice_original + vicen_original(:,:,n)

      !  vsno           = vsno          + vsnon(:,:,n)
      !  vsno_original  = vsno_original + vsnon_original(:,:,n) 
      end do
      
      increment_aice   = aice - aice_original
      !increment_vice   = vice - vice_original
      !increment_vsno   = vsno - vsno_original


     do n = 1, Ncat
        do j = 1, Ny
          do i = 1, Nx
             if ( aice_original(i,j)>0 ) then 
                weight_aicen         = aicen_original(i,j,n)/aice_original(i,j)
                aicen(i,j,n)         = increment_aice(i,j)*weight_aicen + &
                                       aicen_original(i,j,n)
             end if
       !      if ( vice_original(i,j)>0 ) then
       !         weight_vicen         = vicen_original(i,j,n)/vice_original(i,j)
       !         vicen(i,j,n)         = increment_vice(i,j)*weight_vicen + &
       !                                vicen_original(i,j,n)
       !      end if
       !      if ( vsno_original(i,j)>0 ) then
       !         weight_vsnon         = vsnon_original(i,j,n)/vsno_original(i,j)
       !         vsnon(i,j,n)         = increment_vsno(i,j)*weight_vsnon + &
       !                                vsnon_original(i,j,n)
       !      end if
          end do
        end do
     end do

      
END SELECT

     !Fei
     ! the SIMPLE SQUEEZE codes are generic, so I moved the codes toward the end.

     !===============================================================
     !now calculate vsnon = aicen*hsnon_original, and vicen =
     !aicen*hicen_original
     !===============================================================
     
   !  aicen   = max(0.0_r8,aicen)   ! concentrations must be non-negative
     vicen   = max(0.0_r8,vicen)   ! same for volumes 
     vsnon   = max(0.0_r8,vsnon)
     sice001 = max(0.0_r8,sice001) ! same for salinities 
     sice002 = max(0.0_r8,sice002)
     sice003 = max(0.0_r8,sice003)
     sice004 = max(0.0_r8,sice004)
     sice005 = max(0.0_r8,sice005)
     sice006 = max(0.0_r8,sice006)
     sice007 = max(0.0_r8,sice007)
     sice008 = max(0.0_r8,sice008)
     qice001 = min(0.0_r8,qice001) ! enthalpies must be non-positive
     qice002 = min(0.0_r8,qice002)
     qice003 = min(0.0_r8,qice003)
     qice004 = min(0.0_r8,qice004)
     qice005 = min(0.0_r8,qice005)
     qice006 = min(0.0_r8,qice006)
     qice007 = min(0.0_r8,qice007)
     qice008 = min(0.0_r8,qice008)
     qsno001 = min(0.0_r8,qsno001)
     qsno002 = min(0.0_r8,qsno002)
     qsno003 = min(0.0_r8,qsno003)
     aicen   = min(1.0_r8,aicen)    ! concentrations must not exceed 1 
     Tsfcn   = min(Tsmelt,Tsfcn)    ! ice/sno surface must not exceed melting

     ! post-process the parameters
     r_snw   = min(rsnw_max,r_snw)
     r_snw   = max(rsnw_min,r_snw)

     ! calculate aice, which might be negative or >1 at this point
     aice = aicen(:,:,1)
     do n = 2, Ncat  
       aice = aice+aicen(:,:,n)
     enddo

     ! set negative aicen to zero
     aicen   = max(0.0_r8,aicen)

     ! reclaculate aice, now it should be non-negative
     aice_temp = aicen(:,:,1)
     do n = 2, Ncat
        aice_temp = aice_temp + aicen(:,:,n)
     enddo

     ! if aice <0, then set every category to 0
     do j = 1, Ny
        do i = 1, Nx
           if (aice(i,j)<0._r8) then
              aicen(i,j,:) = 0._r8
           endif
        enddo
     enddo
     
     ! aice_temp - aice are the magnitudes of the negative values
     ! this block moves the negative values around
     !DEBUG
     !i = 89
     !j = 372
     !n = 2
     !write(*,*)'aice:',aice(i,j),'aice_temp:',aice_temp(i,j)
     !write(*,*)'aicen:',aicen(i,j,:)
     !write(*,*)'aice_temp-aice',aice_temp(i,j)-aice(i,j)
     !write(*,*)'aicen/aice_temp',aicen(i,j,n)/aice_temp(i,j)
     !write(*,*)'aicen - delta*ratio:',aicen(i,j,n) -(aice_temp(i,j)-aice(i,j))*aicen(i,j,n)/aice_temp(i,j)

     do n=1, Ncat
       do j=1, Ny
          do i=1, Nx
             if (aice_temp(i,j) > 0._r8 .and. aice(i,j)>0._r8) then
                aicen(i,j,n) = aicen(i,j,n) - (aice_temp(i,j)-aice(i,j))*aicen(i,j,n)/aice_temp(i,j)
             endif  
          enddo
       enddo
     enddo

     !i = 89
     !j = 372
     !n = 2
     !write(*,*)'after moving negative values:',aicen(i,j,n)
     ! now squeeze aicen 

      do j = 1, Ny
        do i = 1, Nx
           if (aice(i,j) > 1.0_r8) then
              squeeze        = 1.0_r8 / aice(i,j)
              aicen(i,j,:)   = aicen(i,j,:)*squeeze
           endif
        enddo
     enddo

     !update vsnon and vicen
     where(aicen_original==0) aicen_original = -999

     do n=1,Ncat
        do j=1,Ny
           do i=1,Nx
              hicen_original(i,j,n) = vicen_original(i,j,n)/aicen_original(i,j,n)
              hsnon_original(i,j,n) = vsnon_original(i,j,n)/aicen_original(i,j,n)
           end do
        end do
     end do

     where(hicen_original<0)  hicen_original = 0.0_r8
     where(hsnon_original<0)  hsnon_original = 0.0_r8
     where(aicen_original==-999) aicen_original=0.0_r8
         vicen  = aicen*hicen_original
         vsnon  = aicen*hsnon_original

     cc1 = 3._r8/real(Ncat,kind=r8)
     cc2 = 15.0_r8*cc1
     cc3 = 3._r8
     allocate( hin_max(0:Ncat) )
     allocate( hcat_midpoint(Ncat) )
     hin_max(0) = 0._r8

     do n = 1, ncat
        x1 = real(n-1,kind=r8) / real(Ncat,kind=r8)
        hin_max(n) = hin_max(n-1) &
             + cc1 + cc2*(c1 + tanh(cc3*(x1-c1)))
        hcat_midpoint(n)=0.5_r8*(hin_max(n-1)+hin_max(n))
     enddo

    do n = 1, Ncat
        do j = 1, Ny
           do i = 1, Nx
            
              if (aicen(i,j,n)==0._r8 ) then
                 vicen(i,j,n)   = 0._r8
                 sice001(i,j,n) = 0._r8
                 sice002(i,j,n) = 0._r8
                 sice003(i,j,n) = 0._r8
                 sice004(i,j,n) = 0._r8
                 sice005(i,j,n) = 0._r8
                 sice006(i,j,n) = 0._r8
                 sice007(i,j,n) = 0._r8
                 sice008(i,j,n) = 0._r8

                 qice001(i,j,n) = 0._r8
                 qice002(i,j,n) = 0._r8
                 qice003(i,j,n) = 0._r8
                 qice004(i,j,n) = 0._r8
                 qice005(i,j,n) = 0._r8
                 qice006(i,j,n) = 0._r8
                 qice007(i,j,n) = 0._r8
                 qice008(i,j,n) = 0._r8

                 vsnon(i,j,n) = 0._r8
                 qsno001(i,j,n) = 0._r8
                 qsno002(i,j,n) = 0._r8
                 qsno003(i,j,n) = 0._r8

                 Tsfcn(i,j,n)   = -1.836_r8
               endif

               if (aicen(i,j,n)>0._r8 .and. aicen_original(i,j,n)==0._r8) then

               ! allow no snow volume or enthalpy
                 vsnon(i,j,n) = 0._r8
                 qsno001(i,j,n) = 0._r8
                 qsno002(i,j,n) = 0._r8
                 qsno003(i,j,n) = 0._r8

               ! require ice volume for thickness = category boundary midpoint
                 vicen(i,j,n) =  aicen(i,j,n) * hcat_midpoint(n)

               ! salinity of mushy ice, see add_new_ice in ice_therm_itd.F90
                 Si0new = sss - dSin0_frazil ! given our choice of sss
                 sice001(i,j,n) = Si0new
                 sice002(i,j,n) = Si0new
                 sice003(i,j,n) = Si0new
                 sice004(i,j,n) = Si0new
                 sice005(i,j,n) = Si0new
                 sice006(i,j,n) = Si0new
                 sice007(i,j,n) = Si0new
                 sice008(i,j,n) = Si0new

                 Ti = min(liquidus_temperature_mush(Si0new/phi_init), -0.1_r8)
                 qi0new = enthalpy_mush(Ti, Si0new)
                 qice001(i,j,n) = qi0new
                 qice002(i,j,n) = qi0new
                 qice003(i,j,n) = qi0new
                 qice004(i,j,n) = qi0new
                 qice005(i,j,n) = qi0new
                 qice006(i,j,n) = qi0new
                 qice007(i,j,n) = qi0new
                 qice008(i,j,n) = qi0new

                 Tsfcn(i,j,n) = Ti
              endif
              
            enddo
        enddo
     enddo

   !for testing make something to fix
!  aicen(10,10,1)=1.1
!  write(*,*) (aicen(10,10,k), k=1,5)
   
   varname='aicen'
   io = nf90_inq_varid(ncid, trim(varname), VarID)
   call nc_check(io, 'dart_to_cice', &
                 'inq_varid '//trim(varname)//' '//trim(dart_to_cice_input_file))
   io = nf90_put_var(ncid, VarID, aicen)
   call nc_check(io, 'dart_to_cice', &
                 'put_var '//trim(varname)//' '//trim(dart_to_cice_input_file))
   
   varname='vicen'
   io = nf90_inq_varid(ncid, trim(varname), VarID)
   call nc_check(io, 'dart_to_cice', &
                 'inq_varid '//trim(varname)//' '//trim(dart_to_cice_input_file))
   io = nf90_put_var(ncid, VarID, vicen)
   call nc_check(io, 'dart_to_cice', &
                 'put_var '//trim(varname)//' '//trim(dart_to_cice_input_file))
   
   varname='vsnon'
   io = nf90_inq_varid(ncid, trim(varname), VarID)
   call nc_check(io, 'dart_to_cice', &
                 'inq_varid '//trim(varname)//' '//trim(dart_to_cice_input_file))
   io = nf90_put_var(ncid, VarID, vsnon)
   call nc_check(io, 'dart_to_cice', &
                 'put_var '//trim(varname)//' '//trim(dart_to_cice_input_file))

   varname='Tsfcn'
   io = nf90_inq_varid(ncid, trim(varname), VarID)
   call nc_check(io, 'dart_to_cice', &
                 'inq_varid '//trim(varname)//' '//trim(dart_to_cice_input_file))
   io = nf90_put_var(ncid, VarID, Tsfcn)
   call nc_check(io, 'dart_to_cice', &
                 'put_var '//trim(varname)//' '//trim(dart_to_cice_input_file))

   varname='sice001'
   io = nf90_inq_varid(ncid, trim(varname), VarID)
   call nc_check(io, 'dart_to_cice', &
                 'inq_varid '//trim(varname)//' '//trim(dart_to_cice_input_file))
   io = nf90_put_var(ncid, VarID, sice001)
   call nc_check(io, 'dart_to_cice', &
                 'put_var '//trim(varname)//' '//trim(dart_to_cice_input_file))

   varname='sice002'
   io = nf90_inq_varid(ncid, trim(varname), VarID)
   call nc_check(io, 'dart_to_cice', &
                 'inq_varid '//trim(varname)//' '//trim(dart_to_cice_input_file))
   io = nf90_put_var(ncid, VarID, sice002)
   call nc_check(io, 'dart_to_cice', &
                 'put_var '//trim(varname)//' '//trim(dart_to_cice_input_file))

   varname='sice003'
   io = nf90_inq_varid(ncid, trim(varname), VarID)
   call nc_check(io, 'dart_to_cice', &
                 'inq_varid '//trim(varname)//' '//trim(dart_to_cice_input_file))
   io = nf90_put_var(ncid, VarID, sice003)
   call nc_check(io, 'dart_to_cice', &
                 'put_var '//trim(varname)//' '//trim(dart_to_cice_input_file))

   varname='sice004'
   io = nf90_inq_varid(ncid, trim(varname), VarID)
   call nc_check(io, 'dart_to_cice', &
                 'inq_varid '//trim(varname)//' '//trim(dart_to_cice_input_file))
   io = nf90_put_var(ncid, VarID, sice004)
   call nc_check(io, 'dart_to_cice', &
                 'put_var '//trim(varname)//' '//trim(dart_to_cice_input_file))

   varname='sice005'
   io = nf90_inq_varid(ncid, trim(varname), VarID)
   call nc_check(io, 'dart_to_cice', &
                 'inq_varid '//trim(varname)//' '//trim(dart_to_cice_input_file))
   io = nf90_put_var(ncid, VarID, sice005)
   call nc_check(io, 'dart_to_cice', &
                 'put_var '//trim(varname)//' '//trim(dart_to_cice_input_file))

   varname='sice006'
   io = nf90_inq_varid(ncid, trim(varname), VarID)
   call nc_check(io, 'dart_to_cice', &
                 'inq_varid '//trim(varname)//' '//trim(dart_to_cice_input_file))
   io = nf90_put_var(ncid, VarID, sice006)
   call nc_check(io, 'dart_to_cice', &
                 'put_var '//trim(varname)//' '//trim(dart_to_cice_input_file))

   varname='sice007'
   io = nf90_inq_varid(ncid, trim(varname), VarID)
   call nc_check(io, 'dart_to_cice', &
                 'inq_varid '//trim(varname)//' '//trim(dart_to_cice_input_file))
   io = nf90_put_var(ncid, VarID, sice007)
   call nc_check(io, 'dart_to_cice', &
                 'put_var '//trim(varname)//' '//trim(dart_to_cice_input_file))

   varname='sice008'
   io = nf90_inq_varid(ncid, trim(varname), VarID)
   call nc_check(io, 'dart_to_cice', &
                 'inq_varid '//trim(varname)//' '//trim(dart_to_cice_input_file))
   io = nf90_put_var(ncid, VarID, sice008)
   call nc_check(io, 'dart_to_cice', &
                 'put_var '//trim(varname)//' '//trim(dart_to_cice_input_file))

   varname='qice001'
   io = nf90_inq_varid(ncid, trim(varname), VarID)
   call nc_check(io, 'dart_to_cice', &
                 'inq_varid '//trim(varname)//' '//trim(dart_to_cice_input_file))
   io = nf90_put_var(ncid, VarID, qice001)
   call nc_check(io, 'dart_to_cice', &
                 'put_var '//trim(varname)//' '//trim(dart_to_cice_input_file))

   varname='qice002'
   io = nf90_inq_varid(ncid, trim(varname), VarID)
   call nc_check(io, 'dart_to_cice', &
                 'inq_varid '//trim(varname)//' '//trim(dart_to_cice_input_file))
   io = nf90_put_var(ncid, VarID, qice002)
   call nc_check(io, 'dart_to_cice', &
                 'put_var '//trim(varname)//' '//trim(dart_to_cice_input_file))

   varname='qice003'
   io = nf90_inq_varid(ncid, trim(varname), VarID)
   call nc_check(io, 'dart_to_cice', &
                 'inq_varid '//trim(varname)//' '//trim(dart_to_cice_input_file))
   io = nf90_put_var(ncid, VarID, qice003)
   call nc_check(io, 'dart_to_cice', &
                 'put_var '//trim(varname)//' '//trim(dart_to_cice_input_file))

   varname='qice004'
   io = nf90_inq_varid(ncid, trim(varname), VarID)
   call nc_check(io, 'dart_to_cice', &
                 'inq_varid '//trim(varname)//' '//trim(dart_to_cice_input_file))
   io = nf90_put_var(ncid, VarID, qice004)
   call nc_check(io, 'dart_to_cice', &
                 'put_var '//trim(varname)//' '//trim(dart_to_cice_input_file))

   varname='qice005'
   io = nf90_inq_varid(ncid, trim(varname), VarID)
   call nc_check(io, 'dart_to_cice', &
                 'inq_varid '//trim(varname)//' '//trim(dart_to_cice_input_file))
   io = nf90_put_var(ncid, VarID, qice005)
   call nc_check(io, 'dart_to_cice', &
                 'put_var '//trim(varname)//' '//trim(dart_to_cice_input_file))

   varname='qice006'
   io = nf90_inq_varid(ncid, trim(varname), VarID)
   call nc_check(io, 'dart_to_cice', &
                 'inq_varid '//trim(varname)//' '//trim(dart_to_cice_input_file))
   io = nf90_put_var(ncid, VarID, qice006)
   call nc_check(io, 'dart_to_cice', &
                 'put_var '//trim(varname)//' '//trim(dart_to_cice_input_file))

   varname='qice007'
   io = nf90_inq_varid(ncid, trim(varname), VarID)
   call nc_check(io, 'dart_to_cice', &
                 'inq_varid '//trim(varname)//' '//trim(dart_to_cice_input_file))
   io = nf90_put_var(ncid, VarID, qice007)
   call nc_check(io, 'dart_to_cice', &
                 'put_var '//trim(varname)//' '//trim(dart_to_cice_input_file))

   varname='qice008'
   io = nf90_inq_varid(ncid, trim(varname), VarID)
   call nc_check(io, 'dart_to_cice', &
                 'inq_varid '//trim(varname)//' '//trim(dart_to_cice_input_file))
   io = nf90_put_var(ncid, VarID, qice008)
   call nc_check(io, 'dart_to_cice', &
                 'put_var '//trim(varname)//' '//trim(dart_to_cice_input_file))

   varname='qsno001'
   io = nf90_inq_varid(ncid, trim(varname), VarID)
   call nc_check(io, 'dart_to_cice', &
                 'inq_varid '//trim(varname)//' '//trim(dart_to_cice_input_file))
   io = nf90_put_var(ncid, VarID, qsno001)
   call nc_check(io, 'dart_to_cice', &
                 'put_var '//trim(varname)//' '//trim(dart_to_cice_input_file))

   varname='qsno002'
   io = nf90_inq_varid(ncid, trim(varname), VarID)
   call nc_check(io, 'dart_to_cice', &
                 'inq_varid '//trim(varname)//' '//trim(dart_to_cice_input_file))
   io = nf90_put_var(ncid, VarID, qsno002)
   call nc_check(io, 'dart_to_cice', &
                 'put_var '//trim(varname)//' '//trim(dart_to_cice_input_file))

   varname='qsno003'
   io = nf90_inq_varid(ncid, trim(varname), VarID)
   call nc_check(io, 'dart_to_cice', &
                 'inq_varid '//trim(varname)//' '//trim(dart_to_cice_input_file))
   io = nf90_put_var(ncid, VarID, qsno003)
   call nc_check(io, 'dart_to_cice', &
                 'put_var '//trim(varname)//' '//trim(dart_to_cice_input_file))

   varname=r_snw_name
   io = nf90_inq_varid(ncid, trim(varname), VarID)
   call nc_check(io, 'dart_to_cice', &
                 'inq_varid '//trim(varname)//' '//trim(dart_to_cice_input_file))
   io = nf90_put_var(ncid, VarID, r_snw)
   call nc_check(io, 'dart_to_cice', &
                 'put_var '//trim(varname)//' '//trim(dart_to_cice_input_file))

call nc_check(nf90_close(ncid),'dart_to_cice', 'close '//trim(dart_to_cice_input_file))

deallocate(aicen_original, hin_max, hcat_midpoint)
deallocate( aicen, vicen, vsnon, Tsfcn, aice )
deallocate( sice001, sice002, sice003, sice004, sice005, sice006, sice007, sice008 )
deallocate( qice001, qice002, qice003, qice004, qice005, qice006, qice007, qice008 )
deallocate( qsno001, qsno002, qsno003 )
deallocate(r_snw)

call finalize_utilities('dart_to_cice')

!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------

subroutine get_3d_variable(ncid, varname, var, filename)

integer,               intent(in)  :: ncid
character(len=*),      intent(in)  :: varname
real(r8), allocatable, intent(out) :: var(:,:,:)
character(len=*),      intent(in)  :: filename

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, dimLengths
character(len=NF90_MAX_NAME)          :: dimName

write(msgstring,*) trim(varname)//' '//trim(filename)

io = nf90_inq_varid(ncid, trim(varname), VarID)
call nc_check(io, 'dart_to_cice', 'inq_varid '//trim(msgstring))

io = nf90_inquire_variable(ncid, VarID, dimids=dimIDs, ndims=ndims)
call nc_check(io, 'dart_to_cice', 'inquire_variable '//trim(msgstring))

if (ndims /= 3) then
   write(string2,*) 'expected 3 dimension, got ', ndims
   call error_handler(E_ERR,'dart_to_cice',msgstring,text2=string2)
endif

dimLengths = 1
DimensionLoop : do i = 1,ndims

   write(string1,'(''inquire dimension'',i2,A)') i,trim(msgstring)
   io = nf90_inquire_dimension(ncid, dimIDs(i), name=dimname, len=dimLengths(i))
   call nc_check(io, 'dart_to_cice', string1)

enddo DimensionLoop

allocate( var(dimLengths(1), dimLengths(2), dimLengths(3)) )

call nc_check(nf90_get_var(ncid, VarID, var), 'dart_to_cice', &
         'get_var '//trim(msgstring))

end subroutine get_3d_variable

!==============================================================

subroutine get_2d_variable(ncid, varname, var, filename)

integer,               intent(in)  :: ncid
character(len=*),      intent(in)  :: varname
real(r8), allocatable, intent(out) :: var(:,:)
character(len=*),      intent(in)  :: filename

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, dimLengths
character(len=NF90_MAX_NAME)          :: dimName

write(msgstring,*) trim(varname)//' '//trim(filename)

io = nf90_inq_varid(ncid, trim(varname), VarID)
call nc_check(io, 'dart_to_cice', 'inq_varid '//trim(msgstring))

io = nf90_inquire_variable(ncid, VarID, dimids=dimIDs, ndims=ndims)
call nc_check(io, 'dart_to_cice', 'inquire_variable '//trim(msgstring))

if (ndims /= 2) then
   write(string2,*) 'expected 2 dimension, got ', ndims
   call error_handler(E_ERR,'dart_to_cice',msgstring,text2=string2)
endif

dimLengths = 1
DimensionLoop : do i = 1,ndims

   write(string1,'(''inquire dimension'',i2,A)') i,trim(msgstring)
   io = nf90_inquire_dimension(ncid, dimIDs(i), name=dimname, len=dimLengths(i))
   call nc_check(io, 'dart_to_cice', string1)

enddo DimensionLoop

allocate( var(dimLengths(1), dimLengths(2)) )

call nc_check(nf90_get_var(ncid, VarID, var), 'dart_to_cice', &
         'get_var '//trim(msgstring))

end subroutine get_2d_variable


!=======================================================================
! Mushy Layer Formulation - Assur (1958) liquidus
! functions from cice/src/source/ice_therm_mushy.F90 by Adrian Turner
!=======================================================================

function liquidus_brine_salinity_mush(zTin) result(Sbr)

    ! liquidus relation: equilibrium brine salinity as function of temperature
    ! based on empirical data from Assur (1958)

    real(r8), intent(in) :: &
         zTin         ! ice layer temperature (C)

    real(r8) :: &
         Sbr          ! ice brine salinity (ppt)

    real(r8) :: &
         t_high   , & ! mask for high temperature liquidus region
         lsubzero     ! mask for sub-zero temperatures

    !constant numbers from ice_constants.F90
    real(r8), parameter :: &
         c1      = 1.0_r8 , &
         c1000   = 1000_r8

    ! liquidus relation - higher temperature region
    real(r8), parameter :: &
         az1_liq = -18.48_r8 ,&
         bz1_liq =   0.0_r8

    ! liquidus relation - lower temperature region
    real(r8), parameter :: &
         az2_liq = -10.3085_r8,  &
         bz2_liq =  62.4_r8

    ! liquidus break
    real(r8), parameter :: &
         Tb_liq = -7.6362968855167352_r8
         
    ! basic liquidus relation constants
    real(r8), parameter :: &
         az1p_liq = az1_liq / c1000, &
         bz1p_liq = bz1_liq / c1000, &
         az2p_liq = az2_liq / c1000, &
         bz2p_liq = bz2_liq / c1000

    ! temperature to brine salinity
    real(r8), parameter :: &
       J1_liq = bz1_liq / az1_liq         , &
       K1_liq = c1 / c1000                , &
       L1_liq = (c1 + bz1p_liq) / az1_liq , &
       J2_liq = bz2_liq  / az2_liq        , &
       K2_liq = c1 / c1000                , &
       L2_liq = (c1 + bz2p_liq) / az2_liq

    t_high   = merge(1._r8, 0._r8, (zTin > Tb_liq))
    lsubzero = merge(1._r8, 0._r8, (zTin <= 1._r8))

    Sbr = ((zTin + J1_liq) / (K1_liq * zTin + L1_liq)) * t_high + &
          ((zTin + J2_liq) / (K2_liq * zTin + L2_liq)) * (1._r8 - t_high)

    Sbr = Sbr * lsubzero

  end function liquidus_brine_salinity_mush

!=======================================================================

  function liquidus_temperature_mush(Sbr) result(zTin)

    ! liquidus relation: equilibrium temperature as function of brine salinity
    ! based on empirical data from Assur (1958)

    real(r8), intent(in) :: &
         Sbr    ! ice brine salinity (ppt)

    real(r8) :: &
         zTin   ! ice layer temperature (C)

    real(r8) :: &
         t_high ! mask for high temperature liquidus region

    ! liquidus break
    real(r8), parameter :: &
       Sb_liq =  123.66702800276086_r8    ! salinity of liquidus break

    ! constant numbers from ice_constants.F90
    real(r8), parameter :: &
         c1      = 1.0_r8 , &
         c1000   = 1000_r8

    ! liquidus relation - higher temperature region
    real(r8), parameter :: &
         az1_liq = -18.48_r8 ,&
         bz1_liq =   0.0_r8

    ! liquidus relation - lower temperature region
    real(r8), parameter :: &
         az2_liq = -10.3085_r8,  &
         bz2_liq =  62.4_r8

    ! basic liquidus relation constants
    real(r8), parameter :: &
         az1p_liq = az1_liq / c1000, &
         bz1p_liq = bz1_liq / c1000, &
         az2p_liq = az2_liq / c1000, &
         bz2p_liq = bz2_liq / c1000

  ! brine salinity to temperature
    real(r8), parameter :: &
       M1_liq = az1_liq            , &
       N1_liq = -az1p_liq          , &
       O1_liq = -bz1_liq / az1_liq , &
       M2_liq = az2_liq            , &
       N2_liq = -az2p_liq          , &
       O2_liq = -bz2_liq / az2_liq

    t_high = merge(1._r8, 0._r8, (Sbr <= Sb_liq))

    zTin = ((Sbr / (M1_liq + N1_liq * Sbr)) + O1_liq) * t_high + &
          ((Sbr / (M2_liq + N2_liq * Sbr)) + O2_liq) * (1._r8 - t_high)

  end function liquidus_temperature_mush

!=======================================================================

  function enthalpy_mush(zTin, zSin) result(zqin)

    ! enthalpy of mush from mush temperature and bulk salinity

    real(r8), intent(in) :: &
         zTin, & ! ice layer temperature (C)
         zSin    ! ice layer bulk salinity (ppt)

    real(r8) :: &
         zqin    ! ice layer enthalpy (J m-3) 

    real(r8) :: &
         phi     ! ice liquid fraction 

! from shr_const_mod.F90
    real(r8),parameter :: SHR_CONST_CPSW  = 3.996e3_R8   ! specific heat of sea water ~ J/kg/K
    real(R8),parameter :: SHR_CONST_CPICE = 2.11727e3_R8 ! specific heat of fresh ice ~ J/kg/K
    real(R8),parameter :: SHR_CONST_RHOSW = 1.026e3_R8   ! density of sea water ~ kg/m^3
    real(R8),parameter :: SHR_CONST_RHOICE= 0.917e3_R8   ! density of ice        ~ kg/m^3
    real(R8),parameter :: SHR_CONST_LATICE= 3.337e5_R8   ! latent heat of fusion ~ J/kg


! from cice/src/drivers/cesm/ice_constants.F90
    real(r8) :: cp_ocn, cp_ice, rhoi, rhow, Lfresh

    cp_ice    = SHR_CONST_CPICE  ! specific heat of fresh ice (J/kg/K)
    cp_ocn    = SHR_CONST_CPSW   ! specific heat of ocn    (J/kg/K)
    rhoi      = SHR_CONST_RHOICE ! density of ice (kg/m^3)
    rhow      = SHR_CONST_RHOSW  ! density of seawater (kg/m^3)
    Lfresh    = SHR_CONST_LATICE ! latent heat of melting of fresh ice (J/kg)

    phi = liquid_fraction(zTin, zSin)
 
    zqin = phi * (cp_ocn * rhow - cp_ice * rhoi) * zTin + &
           rhoi * cp_ice * zTin - (1._r8 - phi) * rhoi * Lfresh

  end function enthalpy_mush

  !=======================================================================

  function liquid_fraction(zTin, zSin) result(phi)

    ! liquid fraction of mush from mush temperature and bulk salinity

    real(r8), intent(in) :: &
         zTin, & ! ice layer temperature (C)
         zSin    ! ice layer bulk salinity (ppt)

    real(r8) :: &
         phi , & ! liquid fraction
         Sbr     ! brine salinity (ppt)

    real (r8), parameter :: puny = 1.0e-11_r8 ! cice/src/drivers/cesm/ice_constants.F90

    Sbr = max(liquidus_brine_salinity_mush(zTin),puny)
    phi = zSin / max(Sbr, zSin)

  end function liquid_fraction

!=======================================================================

end program dart_to_cice

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
