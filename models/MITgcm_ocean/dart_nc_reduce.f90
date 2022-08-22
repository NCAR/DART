program nc_reduce

use netcdf_utilities_mod, only : nc_get_variable, nc_define_dimension, nc_define_real_variable, &
                                 nc_put_variable, nc_check, nc_open_file_readonly, &
                                 nc_open_file_readwrite, nc_close_file, nc_create_file, &
                                 nc_get_variable_size, nc_define_double_variable

use types_mod,            only : r4, r8

use utilities_mod,        only : initialize_utilities, finalize_utilities

use netcdf

implicit none

integer            :: ncid, ret, new_ncid
character(len=NF90_MAX_NAME) :: new_name


integer, parameter :: ndim_3d=3
integer, parameter :: ndim_2d=2
real(r4), allocatable  :: psal(:,:,:), ptmp(:,:,:), uvel(:,:,:), vvel(:,:,:)
real(r4), allocatable  :: no3(:,:,:), po4(:,:,:), o2(:,:,:), phy(:,:,:), alk(:,:,:)
real(r4), allocatable  :: dic(:,:,:), dop(:,:,:), don(:,:,:), fet(:,:,:)
real(r4), allocatable  :: eta(:,:), chl(:,:)
real(r4), allocatable :: psal_f(:), ptmp_f(:), uvel_f(:), vvel_f(:)
real(r4), allocatable :: no3_f(:), po4_f(:), o2_f(:), phy_f(:), alk_f(:)
real(r4), allocatable :: dic_f(:), dop_f(:), don_f(:), fet_f(:)
real(r4), allocatable :: eta_f(:), chl_f(:)

! Dimensions
!real(r4)            :: xg(2000), xc(2000), yg(2000), yc(2000)
real(r4)            :: xg(500), xc(500), yg(500), yc(500)
real(r8)            :: zc(50)
integer            :: i,j,k   ! loop counter
integer            :: ct_3d, ct_2d, dimarr_3d_ct, dimarr_2d_ct
integer            :: psalsize(ndim_3d), ptmpsize(ndim_3d), uvelsize(ndim_3d)
integer            :: vvelsize(ndim_3d), no3size(ndim_3d), po4size(ndim_3d)
integer            :: o2size(ndim_3d), physize(ndim_3d), alksize(ndim_3d)
integer            :: dicsize(ndim_3d), dopsize(ndim_3d), donsize(ndim_3d), fetsize(ndim_3d)
integer            :: etasize(ndim_2d), chlsize(ndim_2d)
real(r4), allocatable :: dimarr_3d(:,:)
real(r4), allocatable :: dimarr_2d(:,:)
integer, allocatable :: dimind_3d(:,:)
integer, allocatable :: dimind_2d(:,:)


call initialize_utilities('dart_nc_reduce')

ncid = nc_open_file_readonly('mem01.nc')

call nc_get_variable(ncid, 'XC', xc)
call nc_get_variable(ncid, 'XG', xg)
call nc_get_variable(ncid, 'YC', yc)
call nc_get_variable(ncid, 'YG', yg)
call nc_get_variable(ncid, 'ZC', zc)

write(*,*) 'xc'
write(*,*) xc(3)

write(*,*) 'xg'
write(*,*) xg(3)

write(*,*) 'yc'
write(*,*) yc(3)

write(*,*) 'yg'
write(*,*) yg(3)


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

! ul = size(pack(psal, psal /= -999.0))
! write(*,*) psalsize
! write(*,*) o2size
! write(*,*) etasize

ct_3d = 0
ct_2d = 0
! 
! 
do i=1,psalsize(1)
   do j=1,psalsize(2)
      if (chl(i,j) /= -999.) then
      ct_2d = ct_2d + 1 
      endif
      do k=1,psalsize(3)
         if (psal(i,j,k) /= -999.) then
         ct_3d = ct_3d + 1
         endif
      enddo
   enddo
enddo

allocate(dimarr_3d(ct_3d, 3))
allocate(dimarr_2d(ct_2d, 2))
allocate(dimind_3d(ct_3d, 3))
allocate(dimind_2d(ct_2d, 2))

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

! > EL change 06/23: make the depth the outer loop for this. This will make sure the 2d components 
! > are the first terms of the 3d components. 
do k=1,psalsize(3)
   do i=1,psalsize(1)
      do j=1,psalsize(2)
         if (psal(i,j,k) /= -999.) then
            dimarr_3d(dimarr_3d_ct, 1) = xc(i)
            dimarr_3d(dimarr_3d_ct, 2) = yc(j)
            dimarr_3d(dimarr_3d_ct, 3) = zc(k)
            dimind_3d(dimarr_3d_ct, 1) = i
            dimind_3d(dimarr_3d_ct, 2) = j
            dimind_3d(dimarr_3d_ct, 3) = k 
            
            psal_f(dimarr_3d_ct) = psal(i,j,k)
            ptmp_f(dimarr_3d_ct) = ptmp(i,j,k)
            uvel_f(dimarr_3d_ct) = uvel(i,j,k)
            vvel_f(dimarr_3d_ct) = vvel(i,j,k)
            no3_f(dimarr_3d_ct) = no3(i,j,k)
            po4_f(dimarr_3d_ct) = po4(i,j,k)
            o2_f(dimarr_3d_ct) = o2(i,j,k)
            phy_f(dimarr_3d_ct) = phy(i,j,k)
            alk_f(dimarr_3d_ct) = alk(i,j,k)
            dic_f(dimarr_3d_ct) = dic(i,j,k)
            dop_f(dimarr_3d_ct) = dop(i,j,k)
            don_f(dimarr_3d_ct) = don(i,j,k)
            fet_f(dimarr_3d_ct) = fet(i,j,k)
            dimarr_3d_ct = dimarr_3d_ct + 1
         endif
      enddo
   enddo
enddo

do i=1,chlsize(1)
   do j=1,chlsize(2)
      if (chl(i,j) /= -999.) then
         dimarr_2d(dimarr_2d_ct, 1) = xc(i)
         dimarr_2d(dimarr_2d_ct, 2) = yc(j)
         
         dimind_2d(dimarr_2d_ct, 1) = i
         dimind_2d(dimarr_2d_ct, 2) = j 
         eta_f(dimarr_2d_ct) = eta(i,j)
         chl_f(dimarr_2d_ct) = chl(i,j)
         
         dimarr_2d_ct = dimarr_2d_ct + 1
      endif
   enddo
enddo


! Start creating the new netcdf and define the new 1-d dimension. 
new_name = 'output_mem01.nc'
new_ncid = nc_create_file(new_name, 'squished file')
print*, 'ct_3d', ct_3d, 'ct_2d', ct_2d
call nc_define_dimension(new_ncid, 'useful_info_3d', ct_3d)
call nc_define_dimension(new_ncid, 'useful_info_2d', ct_2d)

! Put all the (new) variables in
call nc_define_real_variable(new_ncid, 'PSAL', 'useful_info_3d')
call nc_define_real_variable(new_ncid, 'PTMP', 'useful_info_3d')
call nc_define_real_variable(new_ncid, 'UVEL', 'useful_info_3d')
call nc_define_real_variable(new_ncid, 'VVEL', 'useful_info_3d')
call nc_define_real_variable(new_ncid, 'ETA', 'useful_info_2d')
call nc_define_real_variable(new_ncid, 'NO3', 'useful_info_3d')
call nc_define_real_variable(new_ncid, 'PO4', 'useful_info_3d')
call nc_define_real_variable(new_ncid, 'O2', 'useful_info_3d')
call nc_define_real_variable(new_ncid, 'PHY', 'useful_info_3d')
call nc_define_real_variable(new_ncid, 'ALK', 'useful_info_3d')
call nc_define_real_variable(new_ncid, 'DIC', 'useful_info_3d')
call nc_define_real_variable(new_ncid, 'DOP', 'useful_info_3d')
call nc_define_real_variable(new_ncid, 'DON', 'useful_info_3d')
call nc_define_real_variable(new_ncid, 'FET', 'useful_info_3d')
call nc_define_real_variable(new_ncid, 'CHL', 'useful_info_2d')
call nc_define_real_variable(new_ncid, 'XC_3D', 'useful_info_3d')
call nc_define_real_variable(new_ncid, 'XC_2D', 'useful_info_2d')
call nc_define_real_variable(new_ncid, 'YC_3D', 'useful_info_3d')
call nc_define_real_variable(new_ncid, 'YC_2D', 'useful_info_2d')
call nc_define_double_variable(new_ncid, 'ZC_3D', 'useful_info_3d')

! Close the file
call nc_close_file(new_ncid)

! Write the information 
new_ncid = nc_open_file_readwrite(new_name)
call nc_put_variable(new_ncid, 'PSAL', psal_f)
call nc_put_variable(new_ncid, 'PTMP', ptmp_f)
call nc_put_variable(new_ncid, 'UVEL', uvel_f)
call nc_put_variable(new_ncid, 'VVEL', vvel_f)
call nc_put_variable(new_ncid, 'ETA', eta_f)
call nc_put_variable(new_ncid, 'NO3', no3_f)
call nc_put_variable(new_ncid, 'PO4', po4_f)
call nc_put_variable(new_ncid, 'O2', o2_f)
call nc_put_variable(new_ncid, 'PHY', phy_f) 
call nc_put_variable(new_ncid, 'ALK', alk_f)
call nc_put_variable(new_ncid, 'DIC', dic_f)
call nc_put_variable(new_ncid, 'DOP', dop_f)
call nc_put_variable(new_ncid, 'DON', don_f)
call nc_put_variable(new_ncid, 'FET', fet_f)
call nc_put_variable(new_ncid, 'CHL', chl_f)
call nc_put_variable(new_ncid, 'XC_3D', dimarr_3d(:, 1))
call nc_put_variable(new_ncid, 'YC_3D', dimarr_3d(:, 2))
call nc_put_variable(new_ncid, 'ZC_3D', dimarr_3d(:, 3))
call nc_put_variable(new_ncid, 'XC_2D', dimarr_2d(:, 1))
call nc_put_variable(new_ncid, 'YC_2D', dimarr_2d(:, 2))


call nc_close_file(new_ncid)

call finalize_utilities('dart_nc_reduce')

end program nc_reduce
