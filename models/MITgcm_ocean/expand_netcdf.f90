! Uncompress a netcdf fil
program expand_netcdf

use netcdf_utilities_mod, only: nc_open_file_readonly, nc_get_dimension_size, &
                                nc_define_dimension, nc_create_file, &
                                nc_get_variable, nc_close_file, nc_put_variable, &
                                nc_define_real_variable, nc_end_define_mode, &
                                nc_add_attribute_to_variable

use types_mod,        only : r4, MISSING_R4

use utilities_mod,    only : initialize_utilities, finalize_utilities
                             
use netcdf

implicit none

integer :: ncid, ncid_comp, dimid(1), dimlen, ret
integer :: Nx,Ny,Nz
integer :: nvars ! total number of variables in compressed file
integer :: id, n, c! loop variables
integer :: i,j,k, ncomp3d, ncomp2d
character(len=NF90_MAX_NAME) :: varname
real(r4), allocatable :: vals3d(:,:,:), vals2d(:,:), vals_comp(:)
integer, allocatable :: Xcomp_ind(:), Ycomp_ind(:), Zcomp_ind(:)

call initialize_utilities('expand_netcdf')

ncid_comp = nc_open_file_readonly('compressed.nc')
ncid = nc_create_file('expanded.nc')

! get the Nx,Ny,Nz
Nx = nc_get_dimension_size(ncid_comp, 'XC')
Ny = nc_get_dimension_size(ncid_comp, 'YC')
Nz = nc_get_dimension_size(ncid_comp, 'ZC')

! define Nx,Ny,Nz in the expanded file
call nc_define_dimension(ncid, 'X', Nx)
call nc_define_dimension(ncid, 'Y', Ny)
call nc_define_dimension(ncid, 'Z', Nz)

! get the compressed size
ncomp2d = nc_get_dimension_size(ncid_comp, 'comp2d')
ncomp3d = nc_get_dimension_size(ncid_comp, 'comp3d')

allocate(vals_comp(ncomp3d))
allocate(vals2d(Nx,Ny), vals3d(Nx,Ny,Nz))

! read in
allocate(Xcomp_ind(ncomp3d), Ycomp_ind(ncomp3d), Zcomp_ind(ncomp3d))
call nc_get_variable(ncid_comp, 'Ycomp_ind', Ycomp_ind)
call nc_get_variable(ncid_comp, 'Xcomp_ind', Xcomp_ind)
call nc_get_variable(ncid_comp, 'Zcomp_ind', Zcomp_ind)


! get the number of variables
ret = nf90_inquire(ncid_comp, nVariables=nvars)

! define variables
do id = 1, nvars
   ret = nf90_inquire_variable(ncid_comp, id, varname, dimids=dimid)

   !  is a it a compressed state variable?
   if (var_of_interest(varname)) then

      ! inquire dimention length (2d or 3d)
      ret = nf90_inquire_dimension(ncid_comp, dimid(1), len=dimlen)
      
      ! define expanded variable
      if (dimlen == ncomp3d) then
         call nc_define_real_variable(ncid, varname, (/'X','Y','Z'/))
      else
         call nc_define_real_variable(ncid, varname, (/'X','Y'/))
      endif

      call nc_add_attribute_to_variable(ncid, varname, 'missing_value', MISSING_R4)

   endif
enddo

call nc_end_define_mode(ncid)

! write variables
do id = 1, nvars
   ret = nf90_inquire_variable(ncid_comp, id, varname, dimids=dimid)

   !  is a it a compressed state variable?
   if (var_of_interest(varname)) then
   
      ! inquire dimention length (2d or 3d)
      ret = nf90_inquire_dimension(ncid_comp, dimid(1), len=dimlen)
   
      ! read in compressed variable
      if (dimlen == ncomp3d) then
         call nc_get_variable(ncid_comp, varname, vals_comp)
         vals3d = MISSING_R4
      else
         call nc_get_variable(ncid_comp, varname, vals_comp(1:ncomp2d))
         vals2d = MISSING_R4
      endif

      ! expand
      c = 1
      do n = 1, ncomp3d
         i = Xcomp_ind(n)
         j = Ycomp_ind(n)
         k = Zcomp_ind(n)
         if (k == 1 .and. dimlen == ncomp2d) then
            vals2d(i,j) = vals_comp(c)
            c = c + 1
         else
            vals3d(i,j,k) = vals_comp(n)
         endif
      enddo
 
      ! write expanded variable
      if (dimlen == ncomp3d) then
         call nc_put_variable(ncid, varname, vals3d)
      else
         call nc_put_variable(ncid, varname, vals2d)
      endif
   
   endif
enddo

call nc_close_file(ncid_comp)
call nc_close_file(ncid)

call finalize_utilities('expand_netcdf')

contains

   ! logical to ignore compression variables
   function var_of_interest(varname)
   character(len=*), intent(in) :: varname
   logical :: var_of_interest
   
   select case (varname)
      case ('XGcomp', 'XCcomp', 'YGcomp', 'YCcomp', 'ZCcomp', 'Xcomp_ind', 'Ycomp_ind', 'Zcomp_ind')
         var_of_interest = .false.
      case ('XC', 'YC', 'ZC', 'XG', 'YG')
         var_of_interest = .false.
      case default
         var_of_interest = .true.
    end select

   end function var_of_interest

end program expand_netcdf
