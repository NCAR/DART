! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program clm_to_dart

!-------------------------------------------------------------------------------
! purpose: Update CLM output files to be suitable as input to filter
!                 Current functionality:
!                 - replaces random values in empty snow layers with a consistent _FillValue
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
                             error_handler, E_MSG, E_ERR

use netcdf_utilities_mod, only : nc_check, &
                                 nc_open_file_readwrite, &
                                 nc_close_file, &
                                 nc_synchronize_file, &
                                 nc_begin_define_mode, &
                                 nc_end_define_mode, &
                                 nc_variable_exists, &
                                 nc_get_variable, &
                                 nc_put_variable, &
                                 nc_get_variable_size, &
                                 nc_get_variable_dimension_names, &
                                 nc_get_attribute_from_variable, &
                                 nc_get_dimension_size

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

real(r8), allocatable :: SNOW_DEPTH(:) ! "snow depth"
real(r8), allocatable :: H2OSNO(:)     ! "snow water" (in column - includes traces)
real(r8), allocatable :: frac_sno(:)   ! "fraction of ground covered by snow (0 to 1)"
integer,  allocatable :: SNLSNO(:)     ! "negative number of snow layers"

real(r8), allocatable :: variable(:,:)
real(r8)              :: FillValue
real(r8)              :: missingValue

character(len=512) :: string1, string2

!===============================================================================

call initialize_utilities(progname='clm_to_dart')

! Read the namelist to get the filename.

call find_namelist_in_file("input.nml", "clm_to_dart_nml", iunit)
read(iunit, nml = clm_to_dart_nml, iostat = io)
call check_namelist_read(iunit, io, "clm_to_dart_nml") ! closes, too.

ncid = nc_open_file_readwrite(clm_restart_file, 'open for unused snow layer value replacement')

call get_snow_metadata()

nlevsno = nc_get_dimension_size(ncid,'levsno') ! The number of snow layers

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

   ! immediately skip the variables that cannot contain snow layers 
   if (xtype /= NF90_DOUBLE) cycle VARIABLES
   if (ndims /= 2)           cycle VARIABLES

   ! Now that we are guaranteed to have a 2D variable - skip anything
   ! that is not dimensioned (*levels*,column)
   call nc_get_variable_dimension_names(ncid, varname, dimnames)
   if (trim(dimnames(2)) /= 'column') cycle VARIABLES

   if ( verbose > 1 ) &
      write(*,*)trim(string1),' varname ',trim(varname), &
             ' dimensions are ', trim(dimnames(1)), ' ', trim(dimnames(2))

   ! For 2D variables, the levels are always the first dimension

   select case (dimnames(1))
      case ( 'levtot', 'levsno')
         if (verbose > 0) then
            write(string1,*)'variable # ',ivar,' is "',trim(varname), &
                            '" - replacing indeterminate values'
            call error_handler(E_MSG,'clm_to_dart',string1)
         endif

         call nc_get_variable_size(ncid, varname, dimlen)
         allocate( variable( dimlen(1),dimlen(2) ) )
         call nc_get_variable(ncid, varname, variable)
         call nc_get_attribute_from_variable(ncid,varname,'_FillValue',FillValue)
         call nc_get_attribute_from_variable(ncid,varname,'missing_value',missingValue)

         ! Replace the bogus values for layers we KNOW to be unused.
         ! The SNLSNO has the negative number of snow layers, so ...

         do j = 1, dimlen(2)  ! loop over columns
            numsnowlevels = abs(SNLSNO(j))

            if (verbose > 2) then
               ! Debug block to check what happens in related variables
               ! when there is no snow in the layer closest to the ground.
               if (numsnowlevels == 0 .and. &
                   variable(nlevsno,j) /= FillValue .and. &
                   variable(nlevsno,j) > 0.0_r8) then
                  write(*,*)trim(varname), nlevsno, j, &
                      variable(nlevsno,j), H2OSNO(j), frac_sno(j), SNOW_DEPTH(j)
               endif
            endif

            
                ! Prevent unused snow layers from being updated by the assimilation.
                ! Unused snow layers have indeterminate values. The indeterminate values are replaced with 
                ! FillValue which prevents filter from updating the unused snow layers during the assimilation.
                do i = 1, nlevsno - numsnowlevels  ! loop over unused layers

                    
                     ! trace amounts of snow are in the level closest to the ground
                     ! frac_sno(j) seems to be a reliable indicator of a trace of snow
                     if (frac_sno(j) > 0.0_r8 .and. i == nlevsno) cycle
                       
                     variable(i,j) = FillValue
            
                enddo

         enddo

         ! update the netCDF file - in place

         call nc_put_variable(ncid,varname,variable,'replacing bogus values')
         deallocate(variable)

      case default
         continue
   end select

enddo VARIABLES

call nc_close_file(ncid)

deallocate(SNLSNO, H2OSNO, frac_sno, SNOW_DEPTH)

call finalize_utilities('clm_to_dart')

!===============================================================================
contains
!===============================================================================


!-------------------------------------------------------------------------------
!> The SNLSNO variable contains the information about how many snow layers are
!> being used in each column.
!> H2OSNO is used to detect the columns that have trace amounts of snow.

subroutine get_snow_metadata()

! If SNLSNO does not exist, there is nothing in 
! the file that has bogus values in the snow layers.

if ( .not. nc_variable_exists(ncid,'SNLSNO') ) then
   write(string1,*)'"'//trim(clm_restart_file)//'" has no "SNLSNO" variable.'
   string2 = 'CLM restart files have "SNLSNO", clm_to_dart only works with restart files.'
   call error_handler(E_ERR,'get_snow_metadata',string1,source,text2=string2)
endif

! All variables are dimensioned the same size
call nc_get_variable_size(ncid, 'SNLSNO', ncolumn)
allocate(SNLSNO(ncolumn),H2OSNO(ncolumn),frac_sno(ncolumn),SNOW_DEPTH(ncolumn))

call nc_get_variable(ncid,'SNLSNO',SNLSNO)
call nc_get_variable(ncid,'frac_sno',frac_sno)
call nc_get_variable(ncid,'SNOW_DEPTH',SNOW_DEPTH)

! ctsm5.1.dev043 replaced the H2OSNO variable with H2OSNO_NO_LAYERS
if ( nc_variable_exists(ncid,'H2OSNO')) then
   call nc_get_variable(ncid,'H2OSNO',H2OSNO)

else if (nc_variable_exists(ncid,'H2OSNO_NO_LAYERS')) then
   call nc_get_variable(ncid,'H2OSNO_NO_LAYERS',H2OSNO)

else
   string1 = 'neither "H2OSNO" nor "H2OSNO_NO_LAYERS" exists in file'
   string2 = 'one or the other is needed'
   call error_handler(E_ERR,'get_snow_metadata',string1,source,text2=string2)
endif

if (verbose > 1) write(*,*)'minval of SNLSNO is ',minval(SNLSNO)

!SNLSNO has a _FillValue of -9999 but there should never be any of these 

call nc_get_attribute_from_variable(ncid,'SNLSNO','_FillValue',FillValue)
if (any(SNLSNO == FillValue)) then
   write(string1,*)'SNLSNO has at least one _FillValue ... unexpected, unable to proceed.'
   call error_handler(E_ERR,'get_snow_metadata',string1,source)
endif

end subroutine get_snow_metadata


end program clm_to_dart
