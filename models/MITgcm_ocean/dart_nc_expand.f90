program nc_reduce

use netcdf_utilities_mod, only : nc_get_variable, nc_define_dimension, nc_define_real_variable, &
                                 nc_put_variable, nc_check, nc_open_file_readonly, &
                                 nc_open_file_readwrite, nc_close_file, nc_create_file, &
                                 nc_get_variable_size

use types_mod,            only : r4

use utilities_mod,        only : initialize_utilities, finalize_utilities

use netcdf

implicit none

integer            :: ncid, ret, new_ncid, ncid_comp
character(len=NF90_MAX_NAME) :: new_name


integer, parameter :: ndim_3d = 3
integer, parameter :: ndim_2d = 2
integer, parameter :: hgrid = 500
integer, parameter :: vgrid = 50

real(r4), allocatable  :: psal(:,:,:), ptmp(:,:,:), uvel(:,:,:), vvel(:,:,:)
real(r4), allocatable  :: no3(:,:,:), po4(:,:,:), o2(:,:,:), phy(:,:,:), alk(:,:,:)
real(r4), allocatable  :: dic(:,:,:), dop(:,:,:), don(:,:,:), fet(:,:,:)
real(r4), allocatable  :: eta(:,:), chl(:,:)
real(r4), allocatable :: psal_f(:), ptmp_f(:), uvel_f(:), vvel_f(:)
real(r4), allocatable :: no3_f(:), po4_f(:), o2_f(:), phy_f(:), alk_f(:)
real(r4), allocatable :: dic_f(:), dop_f(:), don_f(:), fet_f(:)
real(r4), allocatable :: eta_f(:), chl_f(:)

! Dimensions
real(r4)            :: xg(hgrid), xc(hgrid), yg(hgrid), yc(hgrid)
real(r4)            :: zc(vgrid)
integer            :: i,j,k   ! loop counter
integer            :: ct_3d, ct_2d, dimarr_3d_ct, dimarr_2d_ct
integer            :: psalsize(ndim_3d), ptmpsize(ndim_3d), uvelsize(ndim_3d)
integer            :: vvelsize(ndim_3d), no3size(ndim_3d), po4size(ndim_3d)
integer            :: o2size(ndim_3d), physize(ndim_3d), alksize(ndim_3d)
integer            :: dicsize(ndim_3d), dopsize(ndim_3d), donsize(ndim_3d), fetsize(ndim_3d)
integer            :: etasize(ndim_2d), chlsize(ndim_2d)


call initialize_utilities('dart_nc_expand')

ncid = nc_open_file_readonly('mem01.nc')

call nc_get_variable(ncid, 'XC', xc)
call nc_get_variable(ncid, 'XG', xg)
call nc_get_variable(ncid, 'YC', yc)
call nc_get_variable(ncid, 'YG', yg)
call nc_get_variable(ncid, 'ZC', zc)


! Get the size, allocate arrays, and assign values.
call nc_get_variable_size(ncid, 'PSAL', psalsize)
call nc_get_variable_size(ncid, 'PTMP', ptmpsize)
call nc_get_variable_size(ncid, 'UVEL', uvelsize)
call nc_get_variable_size(ncid, 'VVEL', vvelsize)
call nc_get_variable_size(ncid, 'NO3', no3size)
call nc_get_variable_size(ncid, 'PO4', po4size)
call nc_get_variable_size(ncid, 'O2', o2size)
call nc_get_variable_size(ncid, 'PHY', physize)
call nc_get_variable_size(ncid, 'ALK', alksize)
call nc_get_variable_size(ncid, 'DIC', dicsize)
call nc_get_variable_size(ncid, 'DOP', dopsize)
call nc_get_variable_size(ncid, 'DON', donsize)
call nc_get_variable_size(ncid, 'FET', fetsize)
call nc_get_variable_size(ncid, 'ETA', etasize)
call nc_get_variable_size(ncid, 'CHL', chlsize)

allocate(psal(psalsize(1), psalsize(2), psalsize(3)))
call nc_get_variable(ncid, 'PSAL', psal)

allocate(ptmp(ptmpsize(1), ptmpsize(2), ptmpsize(3)))
call nc_get_variable(ncid, 'PTMP', ptmp)

allocate(uvel(uvelsize(1), uvelsize(2), uvelsize(3)))
call nc_get_variable(ncid, 'UVEL', uvel)

allocate(vvel(vvelsize(1), vvelsize(2), vvelsize(3)))
call nc_get_variable(ncid, 'VVEL', vvel)

allocate(no3(no3size(1), no3size(2), no3size(3)))
call nc_get_variable(ncid, 'NO3', no3)

allocate(po4(po4size(1), po4size(2), po4size(3)))
call nc_get_variable(ncid, 'PO4', po4)

allocate(o2(o2size(1), o2size(2), o2size(3)))
call nc_get_variable(ncid, 'O2', o2)

allocate(phy(physize(1), physize(2), physize(3)))
call nc_get_variable(ncid, 'PHY', phy)

allocate(alk(alksize(1), alksize(2), alksize(3)))
call nc_get_variable(ncid, 'ALK', alk)

allocate(dic(dicsize(1), dicsize(2), dicsize(3)))
call nc_get_variable(ncid, 'DIC', dic)

allocate(dop(dopsize(1), dopsize(2), dopsize(3)))
call nc_get_variable(ncid, 'DOP', dop)

allocate(don(donsize(1), donsize(2), donsize(3)))
call nc_get_variable(ncid, 'DON', don)

allocate(fet(fetsize(1), fetsize(2), fetsize(3)))
call nc_get_variable(ncid, 'FET', fet)

allocate(eta(etasize(1), etasize(2)))
call nc_get_variable(ncid, 'ETA', eta)

allocate(chl(chlsize(1), chlsize(2)))
call nc_get_variable(ncid, 'CHL', chl)

! counts are from the compressed file
ncid_comp = nc_open_file_readonly('output_mem01.nc')
call nc_get_variable_size(ncid_comp, 'psal_f', ct_3d)
call nc_get_variable_size(ncid_comp, 'chl_f', ct_2d)


allocate(psal_f(ct_3d))
allocate(ptmp_f(ct_3d))
allocate(uvel_f(ct_3d))
allocate(vvel_f(ct_3d))
allocate(no3_f(ct_3d))
allocate(po4_f(ct_3d))
allocate(o2_f(ct_3d))
allocate(phy_f(ct_3d))
allocate(alk_f(ct_3d))
allocate(dic_f(ct_3d))
allocate(dop_f(ct_3d))
allocate(don_f(ct_3d))
allocate(fet_f(ct_3d))
allocate(chl_f(ct_2d))
allocate(eta_f(ct_2d))


dimarr_3d_ct = 1
dimarr_2d_ct = 1

do k=1,psalsize(3)
   do i=1,psalsize(1)
      do j=1,psalsize(2)
         if (psal(i,j,k) /= -999.) then          
             psal(i,j,k) = psal_f(dimarr_3d_ct) 
             ptmp(i,j,k) = ptmp_f(dimarr_3d_ct) 
             uvel(i,j,k) = uvel_f(dimarr_3d_ct) 
             vvel(i,j,k) = vvel_f(dimarr_3d_ct) 
             no3(i,j,k)  = no3_f(dimarr_3d_ct)  
             po4(i,j,k)  = po4_f(dimarr_3d_ct)  
             o2(i,j,k)   = o2_f(dimarr_3d_ct)   
             phy(i,j,k)  = phy_f(dimarr_3d_ct)  
             alk(i,j,k)  = alk_f(dimarr_3d_ct)  
             dic(i,j,k)  = dic_f(dimarr_3d_ct)  
             dop(i,j,k)  = dop_f(dimarr_3d_ct)  
             don(i,j,k)  = don_f(dimarr_3d_ct)  
             fet(i,j,k)  = fet_f(dimarr_3d_ct)  
            dimarr_3d_ct = dimarr_3d_ct + 1
         endif
      enddo
   enddo
enddo

do i=1,chlsize(1)
   do j=1,chlsize(2)
      if (chl(i,j) /= -999.) then
         
         eta(i,j) = eta_f(dimarr_2d_ct)
         chl(i,j) = chl_f(dimarr_2d_ct)
         
         dimarr_2d_ct = dimarr_2d_ct + 1
      endif
   enddo
enddo


! Start creating the new netcdf and define the new 1-d dimension. 
new_name = 'unsquished_mem01.nc'
new_ncid = nc_create_file(new_name, 'unsquished file')
call nc_define_dimension(new_ncid, 'XG', hgrid)
call nc_define_dimension(new_ncid, 'XC', hgrid)
call nc_define_dimension(new_ncid, 'YG', hgrid)
call nc_define_dimension(new_ncid, 'YC', hgrid)
call nc_define_dimension(new_ncid, 'ZC', vgrid)

! Put all the (new) variables in
call nc_define_real_variable(new_ncid, 'PSAL', (/'XC','YC','ZC'/))
call nc_define_real_variable(new_ncid, 'PTMP', (/'XC','YC','ZC'/))
call nc_define_real_variable(new_ncid, 'UVEL', (/'XC','YC','ZC'/))
call nc_define_real_variable(new_ncid, 'VVEL', (/'XG','YC','ZC'/))
call nc_define_real_variable(new_ncid, 'ETA',  (/'XC','YC'/))
call nc_define_real_variable(new_ncid, 'NO3',  (/'XC','YC','ZC'/))
call nc_define_real_variable(new_ncid, 'PO4',  (/'XC','YC','ZC'/))
call nc_define_real_variable(new_ncid, 'O2',   (/'XC','YC','ZC'/))
call nc_define_real_variable(new_ncid, 'PHY',  (/'XC','YC','ZC'/))
call nc_define_real_variable(new_ncid, 'ALK',  (/'XC','YC','ZC'/))
call nc_define_real_variable(new_ncid, 'DIC',  (/'XC','YC','ZC'/))
call nc_define_real_variable(new_ncid, 'DOP',  (/'XC','YC','ZC'/))
call nc_define_real_variable(new_ncid, 'DON',  (/'XC','YC','ZC'/))
call nc_define_real_variable(new_ncid, 'FET',  (/'XC','YC','ZC'/))
call nc_define_real_variable(new_ncid, 'CHL',  (/'XC','YC'/))

! Close the file
call nc_close_file(new_ncid)

! Write the information 
new_ncid = nc_open_file_readwrite(new_name)
call nc_put_variable(new_ncid, 'PSAL', psal)
call nc_put_variable(new_ncid, 'PTMP', ptmp)
call nc_put_variable(new_ncid, 'UVEL', uvel)
call nc_put_variable(new_ncid, 'VVEL', vvel)
call nc_put_variable(new_ncid, 'ETA', eta)
call nc_put_variable(new_ncid, 'NO3', no3)
call nc_put_variable(new_ncid, 'PO4', po4)
call nc_put_variable(new_ncid, 'O2', o2)
call nc_put_variable(new_ncid, 'PHY', phy)
call nc_put_variable(new_ncid, 'ALK', alk)
call nc_put_variable(new_ncid, 'DIC', dic)
call nc_put_variable(new_ncid, 'DOP', dop)
call nc_put_variable(new_ncid, 'DON', don)
call nc_put_variable(new_ncid, 'FET', fet)
call nc_put_variable(new_ncid, 'CHL', chl)

call nc_put_variable(new_ncid, 'XC', xc)
call nc_put_variable(new_ncid, 'XG', xg)
call nc_put_variable(new_ncid, 'YC', yc)
call nc_put_variable(new_ncid, 'YG', yg)
call nc_put_variable(new_ncid, 'ZC', zc)

call nc_close_file(new_ncid)

call finalize_utilities('dart_nc_reduce')

end program nc_reduce
