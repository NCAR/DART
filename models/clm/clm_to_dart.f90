! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program clm_to_dart

!-------------------------------------------------------------------------------
! purpose: Replace 'bogus' values in unused snow layers with _FillValue
!
! USAGE:  The clm filename is read from the clm_in namelist
!         <edit clm_to_dart_output_file in input.nml:clm_to_dart_nml>
!         clm_to_dart
!
! author: Tim Hoar 12 July 2011
!         and again 2 June 2021
!-------------------------------------------------------------------------------

use        types_mod, only : r8

use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             error_handler, E_MSG

use netcdf_utilities_mod, only : nc_check, &
                                 nc_open_file_readwrite, &
                                 nc_close_file, &
                                 nc_synchronize_file, &
                                 nc_begin_define_mode, &
                                 nc_end_define_mode, &
                                 nc_get_variable, &
                                 nc_put_variable, &
                                 nc_get_variable_size, &
                                 nc_get_variable_dimension_names, &
                                 nc_get_attribute_from_variable, &
                                 nc_get_dimension_size

use time_manager_mod, only : time_type, print_time, print_date

use netcdf

implicit none

character(len=*), parameter :: source = 'clm_to_dart.f90'

!-------------------------------------------------------------------------------
! namelist variables 
!-------------------------------------------------------------------------------

character(len=256) :: clm_restart_file = 'clm_restart.nc'
integer            :: verbose = 0

namelist /clm_to_dart_nml/ clm_restart_file, verbose

!-------------------------------------------------------------------------------
! global storage
!-------------------------------------------------------------------------------

integer :: iunit, ncid, ncolumn, xtype
integer :: io, nvariables, ivar, ndims
integer :: i, j, nlevsno, numsnowlevels
integer :: dimids(NF90_MAX_VAR_DIMS)
integer :: dimlen(NF90_MAX_VAR_DIMS)
character(len=NF90_MAX_NAME) :: dimnames(NF90_MAX_VAR_DIMS)
character(len=NF90_MAX_NAME) :: varname

integer,  allocatable :: SNLSNO(:)     ! "negative number of snow layers"
real(r8), allocatable :: variable(:,:)
real(r8)              :: FillValue

character(len=512) :: string1

!===============================================================================

call initialize_utilities(progname='clm_to_dart')

! Read the namelist to get the filename.

call find_namelist_in_file("input.nml", "clm_to_dart_nml", iunit)
read(iunit, nml = clm_to_dart_nml, iostat = io)
call check_namelist_read(iunit, io, "clm_to_dart_nml") ! closes, too.

ncid = nc_open_file_readwrite(clm_restart_file, 'open for bogus snow value replacement')

! If SNLSNO does not exist, there is nothing in the file that has bogus
! values in the snow layers. Issue warning about wrong file (i.e. not a restart)?

call nc_get_variable_size(ncid, 'SNLSNO', ncolumn)
allocate(SNLSNO(ncolumn))
call nc_get_variable(ncid,'SNLSNO',SNLSNO)
if (verbose > 1) write(*,*)'minval of SNLSNO is ',minval(SNLSNO)

nlevsno = nc_get_dimension_size(ncid,'levsno')

!>@todo  SNLSNO has a _FillValue of -9999 but there should never be any of these 

! Get the number of variables in the netCDF file.
! We will query each of them to see if they use the 'levsno' or 'levtot'
! dimension. If so, we have to replace the bogus values - being careful not
! to clobber any 'trace' snow amounts.

io = nf90_inquire(ncid, nvariables=nvariables)
call nc_check(io, source, 'determining number of variables in file')

if (verbose > 0) write(*,*)'There are ',nvariables,' variables in the file.'

VARIABLES : do ivar = 1,nvariables
   write(string1,*)'inquire variable number ',ivar
   io = nf90_inquire_variable(ncid, ivar, varname, xtype, ndims, dimids)
   call nc_check(io, source, string1)

!  One stop shopping ... 
!  call nc_get_variable_info(ncid, varname, xtype, ndims, dimlen, dimnames)

   if (xtype /= NF90_DOUBLE) cycle VARIABLES
   if (ndims /= 2)           cycle VARIABLES

   call nc_get_variable_dimension_names(ncid, varname, dimnames)

   if (trim(dimnames(2)) /= 'column') cycle VARIABLES

   if ( verbose > 1 ) &
      write(*,*)trim(string1),' varname ',trim(varname), &
             ' dimensions are ', trim(dimnames(1)), ' ', trim(dimnames(2))

   ! For 2D variables, the levels are always the first dimension

   select case (dimnames(1))
      case ( 'levtot', 'levsno' )
         if (verbose > 0) then
            write(string1,*)'variable # ',ivar,' is "',trim(varname),'" - replacing bogus'
            call error_handler(E_MSG,'clm_to_dart',string1)
         endif

         call nc_get_variable_size(ncid, varname, dimlen)
         allocate( variable( dimlen(1),dimlen(2) ) )
         call nc_get_variable(ncid, varname, variable)
         call nc_get_attribute_from_variable(ncid,varname,'_FillValue',FillValue)

         ! replace the bogus values for layers we KNOW to be unused
         ! The SNLSNO has the negative number of snow layers, so ...

         do j = 1, dimlen(2)  ! loop over columns
            numsnowlevels = abs(SNLSNO(j))
            do i = 1, nlevsno - numsnowlevels  ! loop over layers
               variable(i,j) = FillValue
            enddo
         enddo

         ! CLM uses some indeterminate values instead of the _FillValue code 

         if (varname == 'T_SOISNO') then
            where(variable < 1.0_r8) variable = FillValue
            !do j = 1,nj  ! T_SOISNO has missing data in lake columns
            !  if (cols1d_ityplun(j) == icol_deep_lake) then
            !  !  write(*,*)'Found a lake column resetting the following:'
            !  !  write(*,*)variable(:,j)
            !     variable(:,j) = FillValue
            !  endif
            !enddo
         endif

         if ((varname == 'H2OSOI_LIQ')  .or. &
             (varname == 'H2OSOI_ICE')) then
            where(variable < 0.0_r8) variable = FillValue
            !do j = 1,nj  ! missing data in lake columns
            !  if (cols1d_ityplun(j) == icol_deep_lake) then
            !     variable(:,j) = FillValue
            !  endif
            !enddo
         endif

         ! update the netCDF file - in place

         call nc_put_variable(ncid,varname,variable,'replacing bogus values')
         deallocate(variable)

      case default
         continue
   end select

enddo VARIABLES

call nc_close_file(ncid)

deallocate(SNLSNO)

call finalize_utilities('clm_to_dart')

end program clm_to_dart

