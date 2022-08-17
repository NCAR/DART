module netcdf_test

use netcdf
contains

function nc_open_file_readonly(filename, context)

character(len=*), intent(in)  :: filename
character(len=*), intent(in), optional :: context
integer                       :: nc_open_file_readonly

character(len=*), parameter :: routine = 'nc_open_file_readonly'
integer :: ret, ncid

ret = nf90_open(filename, NF90_NOWRITE, ncid)
nc_open_file_readonly = ncid

end function nc_open_file_readonly


subroutine nc_define_var_int_Nd(ncid, varname, dimnames, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: dimnames(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_define_var_int_Nd'
integer :: i, ret, ndims, varid, dimids(NF90_MAX_VAR_DIMS)

ndims = size(dimnames)

do i=1, ndims
   ret = nf90_inq_dimid(ncid, dimnames(i), dimids(i))
enddo

ret = nf90_def_var(ncid, varname, nf90_int, dimids(1:ndims), varid=varid)

end subroutine nc_define_var_int_Nd


function nc_create_file(filename, context)

character(len=*), intent(in)  :: filename
character(len=*), intent(in), optional :: context
integer                       :: nc_create_file

character(len=*), parameter :: routine = 'nc_create_file'
integer :: ret, ncid, oldmode

ret = nf90_create(filename, NF90_CLOBBER, ncid)
nc_create_file = ncid

! faster if we don't fill the vars first with 'fill' value.
! this works if we are planning to write all vars.  (which we are.)

ret = nf90_set_fill(ncid, NF90_NOFILL, oldmode)

end function nc_create_file


subroutine nc_get_variable_size_Nd(ncid, varname, varsize, context, filename)      

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
integer,          intent(out) :: varsize(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_variable_size_Nd'
integer :: ret, i, ndims, varid, dimids(NF90_MAX_VAR_DIMS)


ret = nf90_inq_varid(ncid, varname, varid)

ret = nf90_inquire_variable(ncid, varid, dimids=dimids, ndims=ndims)

! if (ndims > size(varsize)) &
!    call nc_check(NF90_EDIMSIZE, routine, 'variable '//trim(varname)//' return varsize array too small', &
!                  context, filename, ncid)
! 
! ! in case the var is larger than ndims, set unused dims to -1
! varsize(:) = -1
do i=1, ndims
   ret = nf90_inquire_dimension(ncid, dimids(i), len=varsize(i))
enddo

end subroutine nc_get_variable_size_Nd


subroutine nc_get_double_4d(ncid, varname, varvals, context, filename, &
   nc_start, nc_count, nc_stride, nc_map)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
real(8),   intent(out) :: varvals(:,:,:,:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename
integer,          intent(in), optional :: nc_start(:)
integer,          intent(in), optional :: nc_count(:)
integer,          intent(in), optional :: nc_stride(:)
integer,          intent(in), optional :: nc_map(:)

character(len=*), parameter :: routine = 'nc_get_double_4d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
ret = nf90_get_var(ncid, varid, varvals, nc_start, nc_count, nc_stride, nc_map)

end subroutine nc_get_double_4d


subroutine nc_get_real_1d(ncid, varname, varvals, context, filename, &
   nc_start, nc_count, nc_stride, nc_map)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
real(4),   intent(out) :: varvals(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename
integer,          intent(in), optional :: nc_start(:)
integer,          intent(in), optional :: nc_count(:)
integer,          intent(in), optional :: nc_stride(:)
integer,          intent(in), optional :: nc_map(:)

character(len=*), parameter :: routine = 'nc_get_real_1d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
ret = nf90_get_var(ncid, varid, varvals, nc_start, nc_count, nc_stride, nc_map)

end subroutine nc_get_real_1d


subroutine nc_get_variable_size_1d(ncid, varname, varsize, context, filename)      

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
integer,          intent(out) :: varsize
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_variable_size_1d'
integer :: ret, ndims, varid, dimids(NF90_MAX_VAR_DIMS)

ret = nf90_inq_varid(ncid, varname, varid)
ret = nf90_inquire_variable(ncid, varid, dimids=dimids, ndims=ndims)
ret = nf90_inquire_dimension(ncid, dimids(1), len=varsize)

end subroutine nc_get_variable_size_1d


subroutine nc_put_real_1d(ncid, varname, varvals, context, filename, &
   nc_start, nc_count, nc_stride, nc_map)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
real(4),         intent(in) :: varvals(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename
integer,          intent(in), optional :: nc_start(:)
integer,          intent(in), optional :: nc_count(:)
integer,          intent(in), optional :: nc_stride(:)
integer,          intent(in), optional :: nc_map(:)

character(len=*), parameter :: routine = 'nc_put_real_1d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
ret = nf90_put_var(ncid, varid, varvals, nc_start, nc_count, nc_stride, nc_map)

end subroutine nc_put_real_1d


subroutine nc_define_dimension(ncid, dimname, dimlen, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: dimname
integer,          intent(in) :: dimlen
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_define_dimension'
integer :: ret, dimid

ret = nf90_def_dim(ncid, dimname, dimlen, dimid)

end subroutine nc_define_dimension

!--------------------------------------------------------------------

subroutine nc_define_unlimited_dimension(ncid, dimname, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: dimname
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_define_unlimited_dimension'
integer :: ret, dimid

ret = nf90_def_dim(ncid, dimname, NF90_UNLIMITED, dimid)

end subroutine nc_define_unlimited_dimension


function nc_open_file_readwrite(filename, context)

character(len=*), intent(in)  :: filename
character(len=*), intent(in), optional :: context
integer                       :: nc_open_file_readwrite

character(len=*), parameter :: routine = 'nc_open_file_readwrite'
integer :: ret, ncid, oldmode

ret = nf90_open(filename, NF90_WRITE, ncid)
nc_open_file_readwrite = ncid

! faster if we don't fill the vars first with 'fill' value.
! this works if we are planning to write all vars.  (which we are.)

ret = nf90_set_fill(ncid, NF90_NOFILL, oldmode)

end function nc_open_file_readwrite


subroutine nc_close_file(ncid, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_close_file'
integer :: ret

ret = nf90_close(ncid)
end subroutine nc_close_file


subroutine nc_define_var_real_1d(ncid, varname, dimname, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: dimname
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_define_var_real_1d'
integer :: ret, dimid, varid

ret = nf90_inq_dimid(ncid, dimname, dimid)
ret = nf90_def_var(ncid, varname, nf90_real, dimid, varid)

end subroutine nc_define_var_real_1d


subroutine nc_get_real_3d(ncid, varname, varvals, context, filename, &
   nc_start, nc_count, nc_stride, nc_map)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
real(4),         intent(out) :: varvals(:,:,:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename
integer,          intent(in), optional :: nc_start(:)
integer,          intent(in), optional :: nc_count(:)
integer,          intent(in), optional :: nc_stride(:)
integer,          intent(in), optional :: nc_map(:)

character(len=*), parameter :: routine = 'nc_get_real_3d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
ret = nf90_get_var(ncid, varid, varvals, nc_start, nc_count, nc_stride, nc_map)

end subroutine nc_get_real_3d


subroutine nc_get_real_2d(ncid, varname, varvals, context, filename, &
   nc_start, nc_count, nc_stride, nc_map)
integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
real(4),         intent(out) :: varvals(:,:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename
integer,          intent(in), optional :: nc_start(:)
integer,          intent(in), optional :: nc_count(:)
integer,          intent(in), optional :: nc_stride(:)
integer,          intent(in), optional :: nc_map(:)

character(len=*), parameter :: routine = 'nc_get_real_2d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
ret = nf90_get_var(ncid, varid, varvals, nc_start, nc_count, nc_stride, nc_map)

end subroutine nc_get_real_2d


subroutine nc_get_double_1d(ncid, varname, varvals, context, filename, &
   nc_start, nc_count, nc_stride, nc_map)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
real(8),   intent(out) :: varvals(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename
integer,          intent(in), optional :: nc_start(:)
integer,          intent(in), optional :: nc_count(:)
integer,          intent(in), optional :: nc_stride(:)
integer,          intent(in), optional :: nc_map(:)

character(len=*), parameter :: routine = 'nc_get_double_1d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)

ret = nf90_get_var(ncid, varid, varvals, nc_start, nc_count, nc_stride, nc_map)

end subroutine nc_get_double_1d


subroutine nc_put_double_1d(ncid, varname, varvals, context, filename, &
   nc_start, nc_count, nc_stride, nc_map)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
real(8),   intent(in) :: varvals(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename
integer,          intent(in), optional :: nc_start(:)
integer,          intent(in), optional :: nc_count(:)
integer,          intent(in), optional :: nc_stride(:)
integer,          intent(in), optional :: nc_map(:)

character(len=*), parameter :: routine = 'nc_put_double_1d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)

ret = nf90_put_var(ncid, varid, varvals, nc_start, nc_count, nc_stride, nc_map)

end subroutine nc_put_double_1d


subroutine nc_define_var_double_1d(ncid, varname, dimname, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: dimname
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_define_var_double_1d'
integer :: ret, dimid, varid

ret = nf90_inq_dimid(ncid, dimname, dimid)

ret = nf90_def_var(ncid, varname, nf90_double, dimid, varid)

end subroutine nc_define_var_double_1d



end module netcdf_test


program nc_reduce

use netcdf_test

implicit none
integer            :: ncid, status, new_ncid
character(len=NF90_MAX_NAME) :: varname, new_name
integer, parameter :: ndim_3d=3
integer, parameter :: ndim_2d=2
real(4), allocatable  :: psal(:,:,:), ptmp(:,:,:), uvel(:,:,:), vvel(:,:,:)
real(4), allocatable  :: no3(:,:,:), po4(:,:,:), o2(:,:,:), phy(:,:,:), alk(:,:,:)
real(4), allocatable  :: dic(:,:,:), dop(:,:,:), don(:,:,:), fet(:,:,:)
real(4), allocatable  :: eta(:,:), chl(:,:)
real(4), allocatable :: psal_f(:), ptmp_f(:), uvel_f(:), vvel_f(:)
real(4), allocatable :: no3_f(:), po4_f(:), o2_f(:), phy_f(:), alk_f(:)
real(4), allocatable :: dic_f(:), dop_f(:), don_f(:), fet_f(:)
real(4), allocatable :: eta_f(:), chl_f(:)

! Dimensions
!real(4)            :: xg(2000), xc(2000), yg(2000), yc(2000)
real(4)            :: xg(500), xc(500), yg(500), yc(500)
real(4)            :: zc(50)
logical            :: fill_var
integer            :: ul
integer            :: i,j,k   ! loop counter
integer            :: ct_3d, ct_2d, dimarr_3d_ct, dimarr_2d_ct
integer            :: psalsize(ndim_3d), ptmpsize(ndim_3d), uvelsize(ndim_3d)
integer            :: vvelsize(ndim_3d), no3size(ndim_3d), po4size(ndim_3d)
integer            :: o2size(ndim_3d), physize(ndim_3d), alksize(ndim_3d)
integer            :: dicsize(ndim_3d), dopsize(ndim_3d), donsize(ndim_3d), fetsize(ndim_3d)
integer            :: etasize(ndim_2d), chlsize(ndim_2d)
real(4), allocatable :: dimarr_3d(:,:)
real(4), allocatable :: dimarr_2d(:,:)
integer, allocatable :: dimind_3d(:,:)
integer, allocatable :: dimind_2d(:,:)

! The non_nan values in the variable
integer :: non_nan

ncid = nc_open_file_readonly('mem01.nc')

call nc_get_real_1d(ncid, 'XC', xc)
call nc_get_real_1d(ncid, 'XG', xg)
call nc_get_real_1d(ncid, 'YC', yc)
call nc_get_real_1d(ncid, 'YG', yg)
call nc_get_real_1d(ncid, 'ZC', zc)

write(*,*) 'xc'
write(*,*) xc(3)

write(*,*) 'xg'
write(*,*) xg(3)

write(*,*) 'yc'
write(*,*) yc(3)

write(*,*) 'yg'
write(*,*) yg(3)


! Get the size, allocate arrays, and assign values.
call nc_get_variable_size_Nd(ncid, 'PSAL', psalsize)
call nc_get_variable_size_Nd(ncid, 'PTMP', ptmpsize)
call nc_get_variable_size_Nd(ncid, 'UVEL', uvelsize)
call nc_get_variable_size_Nd(ncid, 'VVEL', vvelsize)
call nc_get_variable_size_Nd(ncid, 'NO3', no3size)
call nc_get_variable_size_Nd(ncid, 'PO4', po4size)
call nc_get_variable_size_Nd(ncid, 'O2', o2size)
call nc_get_variable_size_Nd(ncid, 'PHY', physize)
call nc_get_variable_size_Nd(ncid, 'ALK', alksize)
call nc_get_variable_size_Nd(ncid, 'DIC', dicsize)
call nc_get_variable_size_Nd(ncid, 'DOP', dopsize)
call nc_get_variable_size_Nd(ncid, 'DON', donsize)
call nc_get_variable_size_Nd(ncid, 'FET', fetsize)
call nc_get_variable_size_Nd(ncid, 'ETA', etasize)
call nc_get_variable_size_Nd(ncid, 'CHL', chlsize)

allocate(psal(psalsize(1), psalsize(2), psalsize(3)))
call nc_get_real_3d(ncid, 'PSAL', psal)

allocate(ptmp(ptmpsize(1), ptmpsize(2), ptmpsize(3)))
call nc_get_real_3d(ncid, 'PTMP', ptmp)

allocate(uvel(uvelsize(1), uvelsize(2), uvelsize(3)))
call nc_get_real_3d(ncid, 'UVEL', uvel)

allocate(vvel(vvelsize(1), vvelsize(2), vvelsize(3)))
call nc_get_real_3d(ncid, 'VVEL', vvel)

allocate(no3(no3size(1), no3size(2), no3size(3)))
call nc_get_real_3d(ncid, 'NO3', no3)

allocate(po4(po4size(1), po4size(2), po4size(3)))
call nc_get_real_3d(ncid, 'PO4', po4)

allocate(o2(o2size(1), o2size(2), o2size(3)))
call nc_get_real_3d(ncid, 'O2', o2)

allocate(phy(physize(1), physize(2), physize(3)))
call nc_get_real_3d(ncid, 'PHY', phy)

allocate(alk(alksize(1), alksize(2), alksize(3)))
call nc_get_real_3d(ncid, 'ALK', alk)

allocate(dic(dicsize(1), dicsize(2), dicsize(3)))
call nc_get_real_3d(ncid, 'DIC', dic)

allocate(dop(dopsize(1), dopsize(2), dopsize(3)))
call nc_get_real_3d(ncid, 'DOP', dop)

allocate(don(donsize(1), donsize(2), donsize(3)))
call nc_get_real_3d(ncid, 'DON', don)

allocate(fet(fetsize(1), fetsize(2), fetsize(3)))
call nc_get_real_3d(ncid, 'FET', fet)

allocate(eta(etasize(1), etasize(2)))
call nc_get_real_2d(ncid, 'ETA', eta)

allocate(chl(chlsize(1), chlsize(2)))
call nc_get_real_2d(ncid, 'CHL', chl)

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

write(*,*) '3d_values'
write(*,*) no3_f(154311)
write(*,*) dimarr_3d(154311, :)
write(*,*) dimind_3d(154311, :)
write(*,*) '2d_values'
write(*,*) chl_f(154311)
write(*,*) dimarr_2d(154311,:)
write(*,*) dimind_2d(154311, :)
! 
write(*,*) 'original values'
! write(*,*) no3(254,1214,1)
write(*,*) chl(781,1205)

write(*,*) '1-d values'
write(*,*) 

! Start creating the new netcdf and define the new 1-d dimension. 
new_name = 'output_mem01.nc'
status = nf90_create(new_name, NF90_CLOBBER, new_ncid)
call nc_define_dimension(new_ncid, 'useful_info_3d', ct_3d)
call nc_define_dimension(new_ncid, 'useful_info_2d', ct_2d)

! Put all the (new) variables in
call nc_define_var_real_1d(new_ncid, 'PSAL', 'useful_info_3d')
call nc_define_var_real_1d(new_ncid, 'PTMP', 'useful_info_3d')
call nc_define_var_real_1d(new_ncid, 'UVEL', 'useful_info_3d')
call nc_define_var_real_1d(new_ncid, 'VVEL', 'useful_info_3d')
call nc_define_var_real_1d(new_ncid, 'ETA', 'useful_info_2d')
call nc_define_var_real_1d(new_ncid, 'NO3', 'useful_info_3d')
call nc_define_var_real_1d(new_ncid, 'PO4', 'useful_info_3d')
call nc_define_var_real_1d(new_ncid, 'O2', 'useful_info_3d')
call nc_define_var_real_1d(new_ncid, 'PHY', 'useful_info_3d')
call nc_define_var_real_1d(new_ncid, 'ALK', 'useful_info_3d')
call nc_define_var_real_1d(new_ncid, 'DIC', 'useful_info_3d')
call nc_define_var_real_1d(new_ncid, 'DOP', 'useful_info_3d')
call nc_define_var_real_1d(new_ncid, 'DON', 'useful_info_3d')
call nc_define_var_real_1d(new_ncid, 'FET', 'useful_info_3d')
call nc_define_var_real_1d(new_ncid, 'CHL', 'useful_info_2d')
call nc_define_var_real_1d(new_ncid, 'XC_3D', 'useful_info_3d')
call nc_define_var_real_1d(new_ncid, 'XC_2D', 'useful_info_2d')
call nc_define_var_real_1d(new_ncid, 'YC_3D', 'useful_info_3d')
call nc_define_var_real_1d(new_ncid, 'YC_2D', 'useful_info_2d')
call nc_define_var_real_1d(new_ncid, 'ZC_3D', 'useful_info_3d')

! Close the file
call nc_close_file(new_ncid)

! Write the information 
status = nc_open_file_readwrite(new_name)
call nc_put_real_1d(new_ncid, 'PSAL', psal_f)
call nc_put_real_1d(new_ncid, 'PTMP', ptmp_f)
call nc_put_real_1d(new_ncid, 'UVEL', uvel_f)
call nc_put_real_1d(new_ncid, 'VVEL', vvel_f)
call nc_put_real_1d(new_ncid, 'ETA', eta_f)
call nc_put_real_1d(new_ncid, 'NO3', no3_f)
call nc_put_real_1d(new_ncid, 'PO4', po4_f)
call nc_put_real_1d(new_ncid, 'O2', o2_f)
call nc_put_real_1d(new_ncid, 'PHY', phy_f) 
call nc_put_real_1d(new_ncid, 'ALK', alk_f)
call nc_put_real_1d(new_ncid, 'DIC', dic_f)
call nc_put_real_1d(new_ncid, 'DOP', dop_f)
call nc_put_real_1d(new_ncid, 'DON', don_f)
call nc_put_real_1d(new_ncid, 'FET', fet_f)
call nc_put_real_1d(new_ncid, 'CHL', chl_f)
call nc_put_real_1d(new_ncid, 'XC_3D', dimarr_3d(:, 1))
call nc_put_real_1d(new_ncid, 'YC_3D', dimarr_3d(:, 2))
call nc_put_real_1d(new_ncid, 'ZC_3D', dimarr_3d(:, 3))
call nc_put_real_1d(new_ncid, 'XC_2D', dimarr_2d(:, 1))
call nc_put_real_1d(new_ncid, 'YC_2D', dimarr_2d(:, 2))


call nc_close_file(new_ncid)

! Start writing the results: 


end program nc_reduce
