! Data Assimilation Research Testbed -- DART
! Copyright 2004-2009, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
program clm_ens_avg

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!----------------------------------------------------------------------
! purpose: Average an ensemble of CLM initial/restart NetCDF files
!          for use with the ensemble mean CAM files to be used as ICs for forecasts.
!
! method: Read new ensemble file of CLM, 'initial', snow and water fields.
!         Shift snow layers so that top layers line up.  Average them. Shift them back down.
!         Preserve special values in some columns.
!         (Other fields will be averaged using NCO ncea).  Write the averaged snow and water 
!         fields into the pre-existing (from NCO ncea) ensemble averaged file.
!         Not generically coded to handle arbitrary lists of fields in arbitrary order.
!         Little abstraction/object-oriented programming.
!         But does distinguish CLM3 from 3.6 by reading :source attribute from initial files.
!         Updated for CAM3.6; mss_XXX,flx_XXX, etc. fields
!            Note that at least 1 dimension name changed (levsoi -> levgrnd), 
!            and that the the dimension names in 3.5 which were not used (levsoi, ...?), 
!            hence not passed to the new file, may be different in 3.6.
!            Flx and mss fields don't have special values anywhere, but snw_rds should be 
!            averaged only over members with snow in the relevant layer(s).
!          ? flx_ are dimensioned (levsno+1,column), the extra being the first ground layer.
!               NetCDF will allow just (levsno,column) of them to be read, 
!               and write those back out again.  (The ground layer will be averaged by ncea.)
!
! author: Kevin Raeder 4/22/09
! upgrade to CLM 3.6 Kevin Raeder 8/20/09
!         based on model_mod:read_cam_init 
!
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : get_unit, file_exist, nc_check, &
                             initialize_utilities, finalize_utilities
use netcdf
use typeSizes

! eventually use        model_mod, only : read_clm_init
! use  assim_model_mod, only : &
! will want eventually for CLM; write it now?     get_model_size , &

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

character (len = 128) :: file_in = 'snow_water_ens.nc', file_out = 'clm_ens_avg.nc' &
                        ,glob_att_name = 'source', glob_att_val, model_ver

! lo_top is the index of the top layer of snow when all the layers on "on the ground"
! hi_bot is the index of the bottom layer of snow when all the layers have been "raised up"
!        to allign all the top layers of snow into index 1.
! zero_  are the complements to lo_top and hi_bot, used for zeroing the vacated layers.
integer :: file_unit, ncfldid, ncfileid, e,c,k, glob_att_len
integer :: num_col,num_ens, num_levsno, t_soisno_num, lo_top, hi_bot, zero_bot, zero_top

logical                :: do_output = .false. &
                         ,lcam4 = .true.

! Temporary allocatable storage to read in a native format for cam state.

! 2D fields (col,lvl) with ens dimension
real(r8), allocatable, dimension(:,:,:) :: &
          dzsnoE,zsnoE,zisnoE,t_soisnoE,h2osoi_iceE,h2osoi_liqE &
         ,snw_rdsE,qflx_snofrz_lyrE &
         ,mss_bcphoE,mss_bcphiE,mss_ocphoE,mss_ocphiE &
         ,mss_dst1E,mss_dst2E,mss_dst3E,mss_dst4E &
         ,flx_absdvE,flx_absdnE,flx_absivE,flx_absinE
! 2D fields 
real(r8), allocatable, dimension(:,:)   :: &
          dzsno,zsno,zisno,t_soisno,h2osoi_ice,h2osoi_liq &
         ,snw_rds,qflx_snofrz_lyr &
         ,mss_bcpho,mss_bcphi,mss_ocpho,mss_ocphi &
         ,mss_dst1,mss_dst2,mss_dst3,mss_dst4 &
         ,flx_absdv,flx_absdn,flx_absiv,flx_absin
! 1D fields (col) with ens dimension
real(r8), allocatable, dimension(:,:)   :: snowdpE,h2osnoE
! 1D fields (col) 
real(r8), allocatable, dimension(:)     :: snowdp,h2osno
! 1D field with ens dimension for storing # of snow layers
integer,  allocatable, dimension(:,:)   :: snlsnoE
! 1D field 
integer,  allocatable, dimension(:)     :: snlsno
real(r8) :: snowdp_sum, ICE_sum,  LIQ_sum, DZ_sum, SNOWmax, TSOI_sum, SNODIFF, &
            snw_rds_sum,qflx_snofrz_lyr_sum, &
            mss_bcpho_sum,mss_bcphi_sum,mss_ocpho_sum,mss_ocphi_sum, &
            mss_dst1_sum,mss_dst2_sum,mss_dst3_sum,mss_dst4_sum, &
            flx_absdv_sum,flx_absdn_sum,flx_absiv_sum,flx_absin_sum

print*,'==========================================================================='
print*,'clm_ens_avg.f90:'

call initialize_utilities('clm_ens_avg')

! Static init assim model calls static_init_model which now (merge/MPI) calls read_cam_init)
! call static_init_assim_model()

! Read in the array dimensions from the ens file
call nc_check(nf90_open(path = trim(file_in), mode = nf90_write, ncid = ncfileid), &
      'clm_ens_avg', 'opening '//trim(file_in))

! Get model version from the CLM initial file global attribute :source
!     :source = "Community Land Model: CLM3" ;  or CLM3.6
call nc_check(nf90_inquire_attribute(ncfileid, nf90_global, trim(glob_att_name), len = glob_att_len), &
              'clm_ens_avg','inquire; get CAM version')
call nc_check(nf90_get_att(ncfileid, nf90_global, glob_att_name,  glob_att_val), &
              'clm_ens_avg','get_att; get CAM version')
! assim_model example:
! call nc_check(nf90_put_att(ncFileID%ncid, NF90_GLOBAL, "assim_model_source", source ), &

model_ver = glob_att_val(22:glob_att_len)
print*,'glob_att_name, Model version are ',trim(glob_att_name),' ', trim(model_ver)

call read_clm_init_size(ncfileid, num_ens, num_levsno, num_col)

! ens fields
allocate (dzsnoE     (num_levsno,num_col,num_ens), zsnoE           (num_levsno,num_col,num_ens) &
         ,zisnoE     (num_levsno,num_col,num_ens), t_soisnoE       (num_levsno,num_col,num_ens) &
         ,h2osoi_iceE(num_levsno,num_col,num_ens), h2osoi_liqE     (num_levsno,num_col,num_ens) &
         ,snw_rdsE   (num_levsno,num_col,num_ens), qflx_snofrz_lyrE(num_levsno,num_col,num_ens) &
         ,mss_bcphoE (num_levsno,num_col,num_ens), mss_bcphiE      (num_levsno,num_col,num_ens) &
         ,mss_ocphoE (num_levsno,num_col,num_ens), mss_ocphiE      (num_levsno,num_col,num_ens) &
         ,mss_dst1E  (num_levsno,num_col,num_ens), mss_dst2E       (num_levsno,num_col,num_ens) &
         ,mss_dst3E  (num_levsno,num_col,num_ens), mss_dst4E       (num_levsno,num_col,num_ens) &
         ,flx_absdvE (num_levsno,num_col,num_ens), flx_absdnE      (num_levsno,num_col,num_ens) &
         ,flx_absivE (num_levsno,num_col,num_ens), flx_absinE      (num_levsno,num_col,num_ens) )

allocate (snowdpE(num_col,num_ens), h2osnoE(num_col,num_ens), snlsnoE(num_col,num_ens))
! ens avg fields
allocate (dzsno     (num_levsno,num_col), zsno            (num_levsno,num_col) &
         ,zisno     (num_levsno,num_col), t_soisno        (num_levsno,num_col) &
         ,h2osoi_ice(num_levsno,num_col), h2osoi_liq      (num_levsno,num_col) &
         ,snw_rds   (num_levsno,num_col), qflx_snofrz_lyr (num_levsno,num_col) &
         ,mss_bcpho (num_levsno,num_col), mss_bcphi       (num_levsno,num_col) &
         ,mss_ocpho (num_levsno,num_col), mss_ocphi       (num_levsno,num_col) &
         ,mss_dst1  (num_levsno,num_col), mss_dst2        (num_levsno,num_col) &
         ,mss_dst3  (num_levsno,num_col), mss_dst4        (num_levsno,num_col) &
         ,flx_absdv (num_levsno,num_col), flx_absdn       (num_levsno,num_col) &
         ,flx_absiv (num_levsno,num_col), flx_absin       (num_levsno,num_col) )

allocate (snowdp( num_col), h2osno( num_col), snlsno( num_col))

! Initialize output fields to 0.
dzsno           = 0.0_r8
zsno            = 0.0_r8
zisno           = 0.0_r8
t_soisno        = 0.0_r8
h2osoi_ice      = 0.0_r8
h2osoi_liq      = 0.0_r8
snowdp          = 0.0_r8
h2osno          = 0.0_r8
snlsno          = 0
snw_rds         = 0.0_r8
qflx_snofrz_lyr = 0.0_r8
mss_bcpho       = 0.0_r8
mss_bcphi       = 0.0_r8
mss_ocpho       = 0.0_r8
mss_ocphi       = 0.0_r8
mss_dst1        = 0.0_r8
mss_dst2        = 0.0_r8
mss_dst3        = 0.0_r8
mss_dst4        = 0.0_r8
flx_absdv       = 0.0_r8
flx_absdn       = 0.0_r8
flx_absiv       = 0.0_r8
flx_absin       = 0.0_r8


! Read the ensembles of fields from the file of concatenated clminput_#.nc files
!                                      (member,levels,field)
! [?generalize program with:   call read_clm_init(ens(1,1,fld),field_name(fld), num_ens, levsno, fld)]

! Multi-level fields
call read_real(ncfileid, 'DZSNO',           num_levsno,num_col,num_ens, dzsnoE )
call read_real(ncfileid, 'T_SOISNO',        num_levsno,num_col,num_ens, t_soisnoE )
call read_real(ncfileid, 'H2OSOI_ICE',      num_levsno,num_col,num_ens, h2osoi_iceE )
call read_real(ncfileid, 'H2OSOI_LIQ',      num_levsno,num_col,num_ens, h2osoi_liqE )
if (trim(model_ver) == 'CLM3.6') then
   call read_real(ncfileid, 'snw_rds',         num_levsno,num_col,num_ens, snw_rdsE )
   call read_real(ncfileid, 'qflx_snofrz_lyr', num_levsno,num_col,num_ens, qflx_snofrz_lyrE )
   call read_real(ncfileid, 'mss_bcpho',       num_levsno,num_col,num_ens, mss_bcphoE )
   call read_real(ncfileid, 'mss_bcphi',       num_levsno,num_col,num_ens, mss_bcphiE )
   call read_real(ncfileid, 'mss_ocpho',       num_levsno,num_col,num_ens, mss_ocphoE )
   call read_real(ncfileid, 'mss_ocphi',       num_levsno,num_col,num_ens, mss_ocphiE )
   call read_real(ncfileid, 'mss_dst1',        num_levsno,num_col,num_ens, mss_dst1E )
   call read_real(ncfileid, 'mss_dst2',        num_levsno,num_col,num_ens, mss_dst2E )
   call read_real(ncfileid, 'mss_dst3',        num_levsno,num_col,num_ens, mss_dst3E )
   call read_real(ncfileid, 'mss_dst4',        num_levsno,num_col,num_ens, mss_dst4E )
   call read_real(ncfileid, 'flx_absdv',       num_levsno,num_col,num_ens, flx_absdvE )
   call read_real(ncfileid, 'flx_absdn',       num_levsno,num_col,num_ens, flx_absdnE )
   call read_real(ncfileid, 'flx_absiv',       num_levsno,num_col,num_ens, flx_absivE )
   call read_real(ncfileid, 'flx_absin',       num_levsno,num_col,num_ens, flx_absinE )
endif 
   
! Single level fields
call read_2Dreal(ncfileid, 'SNOWDP', num_col,num_ens, snowdpE)
call read_2Dreal(ncfileid, 'H2OSNO', num_col,num_ens, h2osnoE)

! Integer single level fields
call read_2Dint(ncfileid, 'SNLSNO', num_col,num_ens, snlsnoE)

call nc_check(nf90_close(ncfileid), 'clm_ens_avg', 'close input file')

! find maximum number of snow layers amoung the ensemble members for each column.
snlsno = minval(snlsnoE,2)

! Shift snowy columns up to top levels, so that all top layers are in index 1.
do c=1,num_col
   ! average the snow depth for additional filtering of calculations
! Led to wrong values?   snlsno(c) = minval(snlsnoE(c,:))
!   if (c == 1) write(*,'(A,5I4)') 'snlsno, snlsnoE = ',snlsno(c),(snlsnoE (c,e),e=1,num_ens)
!   write(*,'(A,(5I5))') 'col,snlsno, snlsnoE = ',c,snlsno(c),(snlsnoE (c,e),e=1,num_ens)
   SNOWmax = maxval(snowdpE(c,:)) 
   if (SNOWmax > 0.0_r8) then
      ! only need to test first member because these spvals come from a mask, not physical variables
      if (h2osoi_liqE(1,c,1) < 1.e+35) then
         snowdp_sum = 0.0_r8 
         do e=1,num_ens 
            snowdp_sum = snowdp_sum + snowdpE(c,e) 
!            if (c == 1) print*,'snowdpE, snowdp_sum = ',snowdpE(c,e), snowdp_sum 
!           lo_top, hi_bot, ... refer to positions in array as found on the NetCDF file.
!           Those positions are counted in the opposite direction from how CLM numbers them.
            lo_top = num_levsno + snlsnoE(c,e) + 1
            if (lo_top > 1 .and. snlsnoE(c,e) < 0) then
               hi_bot = -snlsnoE(c,e) 
!               write(*,'(A,5I4)') 'snlsnoE,lo_top,hi_bot = ',snlsnoE (c,e),lo_top,hi_bot
!               write(*,'(A,5F8.5)') 'dzsnoE before shift up = ',(dzsnoE(k,c,e),k=1,num_levsno)
               dzsnoE          (1:hi_bot,c,e) = dzsnoE          (lo_top:num_levsno,c,e) 
               h2osoi_liqE     (1:hi_bot,c,e) = h2osoi_liqE     (lo_top:num_levsno,c,e) 
               h2osoi_iceE     (1:hi_bot,c,e) = h2osoi_iceE     (lo_top:num_levsno,c,e) 
               t_soisnoE       (1:hi_bot,c,e) = t_soisnoE       (lo_top:num_levsno,c,e) 
               if (trim(model_ver) == 'CLM3.6') then
                  snw_rdsE        (1:hi_bot,c,e) = snw_rdsE        (lo_top:num_levsno,c,e)          
                  qflx_snofrz_lyrE(1:hi_bot,c,e) = qflx_snofrz_lyrE(lo_top:num_levsno,c,e)
                  mss_bcphoE      (1:hi_bot,c,e) = mss_bcphoE      (lo_top:num_levsno,c,e)        
                  mss_bcphiE      (1:hi_bot,c,e) = mss_bcphiE      (lo_top:num_levsno,c,e)        
                  mss_ocphoE      (1:hi_bot,c,e) = mss_ocphoE      (lo_top:num_levsno,c,e)        
                  mss_ocphiE      (1:hi_bot,c,e) = mss_ocphiE      (lo_top:num_levsno,c,e)        
                  mss_dst1E       (1:hi_bot,c,e) = mss_dst1E       (lo_top:num_levsno,c,e)        
                  mss_dst2E       (1:hi_bot,c,e) = mss_dst2E       (lo_top:num_levsno,c,e)        
                  mss_dst3E       (1:hi_bot,c,e) = mss_dst3E       (lo_top:num_levsno,c,e)        
                  mss_dst4E       (1:hi_bot,c,e) = mss_dst4E       (lo_top:num_levsno,c,e)        
                  flx_absdvE      (1:hi_bot,c,e) = flx_absdvE      (lo_top:num_levsno,c,e)        
                  flx_absdnE      (1:hi_bot,c,e) = flx_absdnE      (lo_top:num_levsno,c,e)        
                  flx_absivE      (1:hi_bot,c,e) = flx_absivE      (lo_top:num_levsno,c,e)        
                  flx_absinE      (1:hi_bot,c,e) = flx_absinE      (lo_top:num_levsno,c,e)        
               endif 

!               write(*,'(A,5F8.5)') ' dzsnoE after shift up = ',(dzsnoE(k,c,e),k=1,num_levsno)
!              0 the rest of the vacated layers
               zero_top = hi_bot +1 
               dzsnoE          (zero_top:num_levsno,c,e) = 0.0_r8 
               h2osoi_liqE     (zero_top:num_levsno,c,e) = 0.0_r8 
               h2osoi_iceE     (zero_top:num_levsno,c,e) = 0.0_r8 
               t_soisnoE       (zero_top:num_levsno,c,e) = 0.0_r8 
               if (trim(model_ver) == 'CLM3.6') then
                  snw_rdsE        (zero_top:num_levsno,c,e) = 0.0_r8
                  qflx_snofrz_lyrE(zero_top:num_levsno,c,e) = 0.0_r8
                  mss_bcphoE      (zero_top:num_levsno,c,e) = 0.0_r8
                  mss_bcphiE      (zero_top:num_levsno,c,e) = 0.0_r8
                  mss_ocphoE      (zero_top:num_levsno,c,e) = 0.0_r8
                  mss_ocphiE      (zero_top:num_levsno,c,e) = 0.0_r8
                  mss_dst1E       (zero_top:num_levsno,c,e) = 0.0_r8
                  mss_dst2E       (zero_top:num_levsno,c,e) = 0.0_r8
                  mss_dst3E       (zero_top:num_levsno,c,e) = 0.0_r8
                  mss_dst4E       (zero_top:num_levsno,c,e) = 0.0_r8
                  flx_absdvE      (zero_top:num_levsno,c,e) = 0.0_r8
                  flx_absdnE      (zero_top:num_levsno,c,e) = 0.0_r8
                  flx_absivE      (zero_top:num_levsno,c,e) = 0.0_r8
                  flx_absinE      (zero_top:num_levsno,c,e) = 0.0_r8
               endif 
!               write(*,'(A,5F8.5)') '    dzsnoE after 0s up = ',(dzsnoE(k,c,e),k=1,num_levsno)
            endif ! full column exclusion
         enddo ! end of ensemble loop
         snowdp(c) = snowdp_sum / num_ens 
!         if (c == 1) print*,'snowdp = ',snowdp(c)
      else 
         snowdp(c) = 0.0_r8 
      endif ! end of test of whether it is spval point
   endif ! end of test of snow presence
enddo ! end of column loop

! Average the re-aligned layers over the ensemble, or the non-0 ens members for T. 
do c=1,num_col  
if (snowdp(c) >= 0.01_r8 ) then
   ! snowdp_sum is just for testing whether the averaging leads to the right accumulation.
   snowdp_sum = 0.0_r8 
   do k=1,num_levsno
      DZ_sum   = 0.0_r8 
      ICE_sum  = 0.0_r8 
      LIQ_sum  = 0.0_r8 
      TSOI_sum = 0.0_r8 
      t_soisno_num = 0 
      if (trim(model_ver) == 'CLM3.6') then
         snw_rds_sum = 0.0_r8
         qflx_snofrz_lyr_sum = 0.0_r8
         mss_bcphi_sum       = 0.0_r8
         mss_bcpho_sum       = 0.0_r8
         mss_ocphi_sum       = 0.0_r8
         mss_ocpho_sum       = 0.0_r8
         mss_dst1_sum        = 0.0_r8
         mss_dst2_sum        = 0.0_r8
         mss_dst3_sum        = 0.0_r8
         mss_dst4_sum        = 0.0_r8
         flx_absdv_sum       = 0.0_r8
         flx_absdn_sum       = 0.0_r8
         flx_absiv_sum       = 0.0_r8
         flx_absin_sum       = 0.0_r8
      endif 
      do e=1,num_ens  
         if (dzsnoE(k,c,e) > 0.0_r8) then
            DZ_sum               = DZ_sum              + dzsnoE          (k,c,e) 
            LIQ_sum              = LIQ_sum             + h2osoi_liqE     (k,c,e) 
            ICE_sum              = ICE_sum             + h2osoi_iceE     (k,c,e) 
            TSOI_sum             = TSOI_sum            + t_soisnoE       (k,c,e) 
            if (trim(model_ver) == 'CLM3.6') then
               snw_rds_sum          = snw_rds_sum         + snw_rdsE        (k,c,e)   
               qflx_snofrz_lyr_sum  = qflx_snofrz_lyr_sum + qflx_snofrz_lyrE(k,c,e)   
               mss_bcphi_sum        = mss_bcphi_sum       + mss_bcphiE      (k,c,e)   
               mss_bcpho_sum        = mss_bcpho_sum       + mss_bcphoE      (k,c,e)   
               mss_ocphi_sum        = mss_ocphi_sum       + mss_ocphiE      (k,c,e)   
               mss_ocpho_sum        = mss_ocpho_sum       + mss_ocphoE      (k,c,e)   
               mss_dst1_sum         = mss_dst1_sum        + mss_dst1E       (k,c,e)   
               mss_dst2_sum         = mss_dst2_sum        + mss_dst2E       (k,c,e)   
               mss_dst3_sum         = mss_dst3_sum        + mss_dst3E       (k,c,e)   
               mss_dst4_sum         = mss_dst4_sum        + mss_dst4E       (k,c,e)   
               flx_absdv_sum        = flx_absdv_sum       + flx_absdvE      (k,c,e)   
               flx_absdn_sum        = flx_absdn_sum       + flx_absdnE      (k,c,e)   
               flx_absiv_sum        = flx_absiv_sum       + flx_absivE      (k,c,e)   
               flx_absin_sum        = flx_absin_sum       + flx_absinE      (k,c,e)   
            endif 
            ! T  and snw_rds averaged only over snowy *members* in this layer.
            t_soisno_num = t_soisno_num + 1 
 
         endif
!         write(*,'(A,5F8.5)') ' DZ_sum    dzsnoE = ',DZ_sum,dzsnoE(k,c,e)
      enddo

      dzsno          (k,c) = DZ_sum              / num_ens 
      h2osoi_liq     (k,c) = LIQ_sum             / num_ens 
      h2osoi_ice     (k,c) = ICE_sum             / num_ens 
      if (trim(model_ver) == 'CLM3.6') then
         qflx_snofrz_lyr(k,c) = qflx_snofrz_lyr_sum / num_ens 
         mss_bcphi      (k,c) = mss_bcphi_sum       / num_ens
         mss_bcpho      (k,c) = mss_bcpho_sum       / num_ens
         mss_ocphi      (k,c) = mss_ocphi_sum       / num_ens
         mss_ocpho      (k,c) = mss_ocpho_sum       / num_ens
         mss_dst1       (k,c) = mss_dst1_sum        / num_ens
         mss_dst2       (k,c) = mss_dst2_sum        / num_ens
         mss_dst3       (k,c) = mss_dst3_sum        / num_ens
         mss_dst4       (k,c) = mss_dst4_sum        / num_ens
         flx_absdv      (k,c) = flx_absdv_sum       / num_ens
         flx_absdn      (k,c) = flx_absdn_sum       / num_ens
         flx_absiv      (k,c) = flx_absiv_sum       / num_ens
         flx_absin      (k,c) = flx_absin_sum       / num_ens
      endif 
      if (t_soisno_num > 0) then
         t_soisno(k,c) = TSOI_sum / t_soisno_num 
         if (trim(model_ver) == 'CLM3.6') snw_rds (k,c) = snw_rds_sum / t_soisno_num
      endif

!     Derived fields
      h2osno(c) = h2osno(c) + h2osoi_liq(k,c) + h2osoi_ice(k,c) 
      snowdp_sum = snowdp_sum + dzsno(k,c) 
!      write(*,'(A,5F8.5)') 'snowdp_sum    dzsno = ',snowdp_sum,dzsno(k,c)

   enddo  ! end of k loop

!  Check whether averaged snow total equals sum of average snow layers.
   SNODIFF = snowdp_sum - snowdp(c) 
   if (SNODIFF < 0.0_r8) SNODIFF = -SNODIFF 
   if (SNODIFF > 0.01_r8) then
      print*,"Column sum of dzsno != snowdp ", c, snowdp_sum, snowdp(c)
   endif

else if (snowdp(c) > 0.0_r8 .and. snowdp(c) < 0.01_r8 ) then
!   print*,'low snow column set to 0; ',c
   snlsno(c) = 0 
   snowdp(c) = 0.0_r8 
   dzsno     (:,c) = 0.0_r8 
   t_soisno       (1:num_levsno,c) = 0.0_r8 
   h2osoi_liq     (1:num_levsno,c) = 0.0_r8 
   h2osoi_ice     (1:num_levsno,c) = 0.0_r8 
   if (trim(model_ver) == 'CLM3.6') then
      snw_rds        (1:num_levsno,c) = 0.0_r8
      qflx_snofrz_lyr(1:num_levsno,c) = 0.0_r8
      mss_bcpho      (1:num_levsno,c) = 0.0_r8
      mss_bcphi      (1:num_levsno,c) = 0.0_r8
      mss_ocpho      (1:num_levsno,c) = 0.0_r8
      mss_ocphi      (1:num_levsno,c) = 0.0_r8
      mss_dst1       (1:num_levsno,c) = 0.0_r8
      mss_dst2       (1:num_levsno,c) = 0.0_r8
      mss_dst3       (1:num_levsno,c) = 0.0_r8
      mss_dst4       (1:num_levsno,c) = 0.0_r8
      flx_absdv      (1:num_levsno,c) = 0.0_r8
      flx_absdn      (1:num_levsno,c) = 0.0_r8
      flx_absiv      (1:num_levsno,c) = 0.0_r8
      flx_absin      (1:num_levsno,c) = 0.0_r8
   endif 

else if (snowdp(c) == 0.0_r8) then
   ! Need this copy because some columns have 1.e+36
   ! Assumes that snowdp=0 will find all the sp val columns,
   ! and that all the ens members have the same sp vals
   ! or 0s in the snow layers. 
   ! No copy needed for mss_XXX, flx_XXX because they never have spvals
   t_soisno  (1:num_levsno,c) = t_soisnoE  (1:num_levsno,c,1) 
   h2osoi_liq(1:num_levsno,c) = h2osoi_liqE(1:num_levsno,c,1) 
   h2osoi_ice(1:num_levsno,c) = h2osoi_iceE(1:num_levsno,c,1) 
endif ! end of snowdp if
enddo ! end of column loop

! Re-align the averaged fields so that bottom layers line up at the bottom.
do c=1,num_col 
if (snowdp(c) >= 0.01_r8 ) then        ! Only for columns with snow
   hi_bot = - snlsno(c) 
   lo_top = num_levsno + snlsno(c) +1
   if (hi_bot < num_levsno) then     ! Only for columns that are not full
      h2osoi_liq     (lo_top:num_levsno,c) = h2osoi_liq     (1:hi_bot,c) 
      h2osoi_ice     (lo_top:num_levsno,c) = h2osoi_ice     (1:hi_bot,c) 
      t_soisno       (lo_top:num_levsno,c) = t_soisno       (1:hi_bot,c) 
      dzsno          (lo_top:num_levsno,c) = dzsno          (1:hi_bot,c) 
      if (trim(model_ver) == 'CLM3.6') then
         snw_rds        (lo_top:num_levsno,c) = snw_rds        (1:hi_bot,c)
         qflx_snofrz_lyr(lo_top:num_levsno,c) = qflx_snofrz_lyr(1:hi_bot,c)
         mss_bcpho      (lo_top:num_levsno,c) = mss_bcpho      (1:hi_bot,c)
         mss_bcphi      (lo_top:num_levsno,c) = mss_bcphi      (1:hi_bot,c)
         mss_ocpho      (lo_top:num_levsno,c) = mss_ocpho      (1:hi_bot,c)
         mss_ocphi      (lo_top:num_levsno,c) = mss_ocphi      (1:hi_bot,c)
         mss_dst1       (lo_top:num_levsno,c) = mss_dst1       (1:hi_bot,c)
         mss_dst2       (lo_top:num_levsno,c) = mss_dst2       (1:hi_bot,c)
         mss_dst3       (lo_top:num_levsno,c) = mss_dst3       (1:hi_bot,c)
         mss_dst4       (lo_top:num_levsno,c) = mss_dst4       (1:hi_bot,c)
         flx_absdv      (lo_top:num_levsno,c) = flx_absdv      (1:hi_bot,c)
         flx_absdn      (lo_top:num_levsno,c) = flx_absdn      (1:hi_bot,c)
         flx_absiv      (lo_top:num_levsno,c) = flx_absiv      (1:hi_bot,c)
         flx_absin      (lo_top:num_levsno,c) = flx_absin      (1:hi_bot,c)
      endif 

      ! 0 the rest of the vacated layers
      zero_bot = lo_top -1 
      dzsno(1:zero_bot,c) = 0.0_r8 
      zsno (1:zero_bot,c) = 0.0_r8 
      zisno(1:zero_bot,c) = 0.0_r8 
      if (h2osoi_liq(1,c) < 1.e+35) then
         ! Some columns of these fields are assigned special
         ! values in assigment loop (above) 
         h2osoi_liq     (1:zero_bot,c) = 0.0_r8 
         h2osoi_ice     (1:zero_bot,c) = 0.0_r8 
         t_soisno       (1:zero_bot,c) = 0.0_r8 
         if (trim(model_ver) == 'CLM3.6') then
            snw_rds        (1:zero_bot,c) = 0.0_r8
            qflx_snofrz_lyr(1:zero_bot,c) = 0.0_r8
            mss_bcpho      (1:zero_bot,c) = 0.0_r8
            mss_bcphi      (1:zero_bot,c) = 0.0_r8
            mss_ocpho      (1:zero_bot,c) = 0.0_r8
            mss_ocphi      (1:zero_bot,c) = 0.0_r8
            mss_dst1       (1:zero_bot,c) = 0.0_r8
            mss_dst2       (1:zero_bot,c) = 0.0_r8
            mss_dst3       (1:zero_bot,c) = 0.0_r8
            mss_dst4       (1:zero_bot,c) = 0.0_r8
            flx_absdv      (1:zero_bot,c) = 0.0_r8
            flx_absdn      (1:zero_bot,c) = 0.0_r8
            flx_absiv      (1:zero_bot,c) = 0.0_r8
            flx_absin      (1:zero_bot,c) = 0.0_r8
         endif 
      endif
   endif ! full column exclusion
   ! Calculate new snow layer depths and interfaces.
   !  zisno are the tops of the layers, 5(num_levsno) is the bottom layer.
   !  Distances are measured upward, but called < 0,  from the bottom(0.)
   !  Don't forget that snlsno are also < 0
   zsno  (num_levsno,c) = -(dzsno(num_levsno,c) * 0.5_r8) 
   zisno (num_levsno,c) = - dzsno(num_levsno,c) 
   do k = num_levsno-1, num_levsno+snlsno(c)+1, -1
       zisno(k,c) = zisno(k+1,c) - dzsno(k,c) 
       zsno (k,c) = zisno(k,c)  + (dzsno(k,c) * 0.5_r8) 
   enddo ! end of layer loop
endif ! end of snowdp if
enddo ! end of column loop

! Open the pre-existing, ensemble averaged CLM initial file, which has bad snow and
! water fields in it, which will be replaced here.
! Get channel for output.
call nc_check(nf90_open(path = trim(file_out), mode = nf90_write, ncid = ncfileid), &
           'clm_ens_avg', 'opening '//trim(file_out))

! Multi-level fields
call write_real(ncfileid, 'DZSNO',           num_levsno,num_col, dzsno )
call write_real(ncfileid, 'T_SOISNO',        num_levsno,num_col, t_soisno )
call write_real(ncfileid, 'H2OSOI_ICE',      num_levsno,num_col, h2osoi_ice )
call write_real(ncfileid, 'H2OSOI_LIQ',      num_levsno,num_col, h2osoi_liq )
if (trim(model_ver) == 'CLM3.6') then
   call write_real(ncfileid, 'snw_rds',         num_levsno,num_col, snw_rds )
   call write_real(ncfileid, 'qflx_snofrz_lyr', num_levsno,num_col, qflx_snofrz_lyr )
   call write_real(ncfileid, 'mss_bcpho',       num_levsno,num_col, mss_bcpho )
   call write_real(ncfileid, 'mss_bcphi',       num_levsno,num_col, mss_bcphi )
   call write_real(ncfileid, 'mss_ocpho',       num_levsno,num_col, mss_ocpho )
   call write_real(ncfileid, 'mss_ocphi',       num_levsno,num_col, mss_ocphi )
   call write_real(ncfileid, 'mss_dst1',        num_levsno,num_col, mss_dst1 )
   call write_real(ncfileid, 'mss_dst2',        num_levsno,num_col, mss_dst2 )
   call write_real(ncfileid, 'mss_dst3',        num_levsno,num_col, mss_dst3 )
   call write_real(ncfileid, 'mss_dst4',        num_levsno,num_col, mss_dst4 )
   call write_real(ncfileid, 'flx_absdv',       num_levsno,num_col, flx_absdv )
   call write_real(ncfileid, 'flx_absdn',       num_levsno,num_col, flx_absdn )
   call write_real(ncfileid, 'flx_absiv',       num_levsno,num_col, flx_absiv )
   call write_real(ncfileid, 'flx_absin',       num_levsno,num_col, flx_absin )
endif 

! Single level fields
call write_1Dreal(ncfileid, 'SNOWDP', num_col, snowdp)
call write_1Dreal(ncfileid, 'H2OSNO', num_col, h2osno)

! Single level integer fields
call write_1Dint(ncfileid, 'SNLSNO', num_col, snlsno)

call nc_check(nf90_close(ncfileid), 'clm_ens_avg', 'close output file')

deallocate (dzsnoE     , zsnoE            &
           ,zisnoE     , t_soisnoE        &
           ,h2osoi_iceE, h2osoi_liqE      &
           ,snw_rdsE   , qflx_snofrz_lyrE &
           ,mss_bcphoE , mss_bcphiE       &
           ,mss_ocphoE , mss_ocphiE       &
           ,mss_dst1E  , mss_dst2E        &
           ,mss_dst3E  , mss_dst4E        &
           ,flx_absdvE , flx_absdnE       &
           ,flx_absivE , flx_absinE       )
deallocate (snowdpE, h2osnoE, snlsnoE)
deallocate (dzsno     , zsno             &
           ,zisno     , t_soisno         &
           ,h2osoi_ice, h2osoi_liq       &
           ,snw_rds   , qflx_snofrz_lyr  &
           ,mss_bcpho , mss_bcphi        &
           ,mss_ocpho , mss_ocphi        &
           ,mss_dst1  , mss_dst2         &
           ,mss_dst3  , mss_dst4         &
           ,flx_absdv , flx_absdn        &
           ,flx_absiv , flx_absin        )
deallocate (snowdp, h2osno, snlsno)

call finalize_utilities()

contains


subroutine read_clm_init_size(ncfileid, num_ens, num_levsno, num_col)
!=======================================================================
! subroutine read_clm_init_size(ncfileid)

!
! Gets the number, names, and sizes of field dimensions from a CAM init netcdf file
! in file_name (regardless of dynamical core).
! Called by static_init_model (only).

integer,  intent(in)  :: ncfileid
integer,  intent(out) :: num_ens, num_levsno, num_col

integer :: i,j, num_dims
integer, allocatable, dimension(:) :: dim_sizes
character(len=NF90_MAX_NAME), allocatable  :: dim_names(:)

!------------------------------------
! learn how many dimensions are defined in this file.
call nc_check(nf90_inquire(ncfileid, num_dims), 'read_clm_init_size', 'inquire num_dims')

! where to deallocate?
allocate (dim_names(num_dims), dim_sizes(num_dims))

! Cycle through dimids until there aren't any more
! Dimension ids are sequential integers on the NetCDF file.
do i = 1,num_dims
   call nc_check(nf90_inquire_dimension(ncfileid, i, dim_names(i), dim_sizes(i)), &
                 'read_clm_init_size', 'inquire for '//trim(dim_names(i)))
   if (dim_names(i) == 'ensemble') num_ens    = dim_sizes(i)
   if (dim_names(i) == 'levsno')   num_levsno = dim_sizes(i)
   if (dim_names(i) == 'column')   num_col    = dim_sizes(i)
   write(*,*) 'Dims info = ',i, trim(dim_names(i)), dim_sizes(i)
end do


end subroutine read_clm_init_size

subroutine read_real(ncfileid,cfield,dim1,dim2,dim3,field)
!=======================================================================
! subroutine read_real(ncfileid,cfield,dim1,dim2,dim3,field)

character(len=*),           intent(in) :: cfield
integer,                    intent(in)  :: ncfileid,dim1,dim2,dim3
real(r8), dimension(:,:,:), intent(out) :: field
! 

integer :: ncfldid

call nc_check(nf90_inq_varid(ncfileid, trim(cfield), ncfldid), 'read_real', 'inq_varid '//trim(cfield))
PRINT*,'read ',trim(cfield),'  using id ',ncfldid
call nc_check(nf90_get_var  (ncfileid, ncfldid, field ,start=(/1,1,1/) &
             ,count=(/dim1,dim2,dim3/)), 'read_real', 'get_var '//trim(cfield))

return
end subroutine read_real

subroutine read_2Dreal(ncfileid,cfield,dim1,dim2,field)
!=======================================================================
! subroutine read_2Dreal(ncfileid,cfield,dim1,dim2,field)

character(len=*),         intent(in)  :: cfield
integer,                  intent(in)  :: ncfileid,dim1,dim2
real(r8), dimension(:,:), intent(out) :: field
! 

integer :: ncfldid

call nc_check(nf90_inq_varid(ncfileid, trim(cfield), ncfldid), 'read_2Dreal', 'inq_varid '//trim(cfield))
PRINT*,'read ',trim(cfield),'  using id ',ncfldid
call nc_check(nf90_get_var  (ncfileid, ncfldid, field ,start=(/1,1/) &
             ,count=(/dim1,dim2/)), 'read_2Dreal', 'get_var '//trim(cfield))

return
end subroutine read_2Dreal

subroutine read_2Dint(ncfileid,cfield,dim1,dim2,field)
!=======================================================================
! subroutine read_2Dint(ncfileid,cfield,dim1,dim2,field)

character(len=*), intent(in) :: cfield
integer, intent(in)  :: ncfileid,dim1,dim2
integer, dimension(:,:),   intent(out) :: field
! 

integer :: ncfldid

call nc_check(nf90_inq_varid(ncfileid, trim(cfield), ncfldid), 'read_2Dint', 'inq_varid '//trim(cfield))
PRINT*,'read ',trim(cfield),'  using id ',ncfldid
call nc_check(nf90_get_var  (ncfileid, ncfldid, field ,start=(/1,1/) &
             ,count=(/dim1,dim2/)), 'read_2Dint', 'get_var '//trim(cfield))

return
end subroutine read_2Dint

subroutine write_real(ncfileid,cfield,dim1,dim2,field)
!=======================================================================
! subroutine write_real(ncfileid,cfield,dim1,dim2,field)

character(len=*),          intent(in) :: cfield
integer,                   intent(in) :: ncfileid,dim1,dim2
real(r8),  dimension(:,:), intent(in) :: field


integer :: ncfldid

call nc_check(nf90_inq_varid(ncfileid, trim(cfield), ncfldid), 'write_real', 'inq_varid '//trim(cfield))
PRINT*,'wrote ',trim(cfield),'  using id ',ncfldid

if (dim2 == 1) then
   call nc_check(nf90_put_var  (ncfileid, ncfldid, field ,start=(/1/) &
                ,count=(/dim1/)), 'write_real', 'put_var '//trim(cfield))
else
   call nc_check(nf90_put_var  (ncfileid, ncfldid, field ,start=(/1,1/) &
                ,count=(/dim1,dim2/)), 'write_real', 'put_var '//trim(cfield))
endif

return
end subroutine write_real

subroutine write_1Dreal(ncfileid,cfield,dim1,field)
!=======================================================================
! subroutine write_1Dreal(ncfileid,cfield,dim1,field)

character(len=*),        intent(in) :: cfield
integer,                 intent(in) :: ncfileid,dim1
real(r8),  dimension(:), intent(in) :: field


integer :: ncfldid

call nc_check(nf90_inq_varid(ncfileid, trim(cfield), ncfldid), 'write_1Dfield', 'inq_varid '//trim(cfield))
PRINT*,'wrote ',trim(cfield),'  using id ',ncfldid

call nc_check(nf90_put_var  (ncfileid, ncfldid, field ,start=(/1/) &
             ,count=(/dim1/)), 'write_1Dreal', 'put_var '//trim(cfield))

return
end subroutine write_1Dreal

subroutine write_1Dint(ncfileid,cfield,dim1,field)
!=======================================================================
! subroutine write_1Dint(ncfileid,cfield,dim1,field)

character(len=*),       intent(in) :: cfield
integer,                intent(in) :: ncfileid,dim1
integer,  dimension(:), intent(in) :: field


integer :: ncfldid

call nc_check(nf90_inq_varid(ncfileid, trim(cfield), ncfldid), 'write_1Dint', 'inq_varid '//trim(cfield))
PRINT*,'wrote ',trim(cfield),'  using id ',ncfldid

call nc_check(nf90_put_var  (ncfileid, ncfldid, field ,start=(/1/) &
             ,count=(/dim1/)), 'write_1Dint', 'put_var '//trim(cfield))

return
end subroutine write_1Dint

end program clm_ens_avg
