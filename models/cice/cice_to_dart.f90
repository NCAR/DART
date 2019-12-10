! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id$

!> This program reads in ensemble CICE parameters (2D) and 
!> attach them to the CICE restart files

program cice_to_dart

use        types_mod, only : r8,metadatalength
use    utilities_mod, only : error_handler, E_MSG, &
                             find_namelist_in_file, check_namelist_read, &
                             nmlfileunit,do_nml_file, do_nml_term, &
                             initialize_utilities, finalize_utilities
use netcdf_utilities_mod, only : nc_check
use netcdf

!version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = &
"$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"

! Contents specified in the input.nml: &cice_parameter_nml namelist
integer, parameter :: max_parameters = 10
integer, parameter :: num_parameter_columns = 1
character(len=NF90_MAX_NAME) :: parameter_table(max_parameters,num_parameter_columns)

character(len=256) :: cice_restart_input_file = 'cice_restart.nc'
character(len=256) :: parameter_input_file    = 'parameter_prior.nc'
character(len=metadatalength) :: cice_parameters(max_parameters) = ''
 
namelist /cice_parameter_nml/      &
         cice_restart_input_file,  &
         parameter_input_file,     &
         cice_parameters

integer :: iunit, io, ios
integer :: ipar   ! loop index for parameters
integer :: npar   ! number of parameters

integer :: varid, ncid,ncid2, dimid1, dimid2
real(r8), allocatable :: par(:,:,:)

! message string
character(len=512) :: string1

call initialize_utilities('cice_to_dart')

! read the namelist
call find_namelist_in_file('input.nml','cice_parameter_nml',iunit)
read(iunit, nml = cice_parameter_nml, iostat = io)
call check_namelist_read(iunit,io, 'cice_parameter_nml')

! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=cice_parameter_nml)
if (do_nml_term()) write(     *     , nml=cice_parameter_nml)

! verify the cice_parameters namelist was filled in correctly
! returns cice_parameters which has parameter names and number of parameters
call verify_parameters(cice_parameters,npar,parameter_table)
! open the parameter input file
call nc_check(nf90_open(trim(parameter_input_file),NF90_nowrite,ncid), &
              'cice_to_dart','open '//trim(parameter_input_file))
! get dims
call nc_check(nf90_inq_dimid(ncid,'ni',dimid1), &
              'cice_to_dart','inq dimid ni')
call nc_check(nf90_inquire_dimension(ncid,dimid1,len=ni), &
              'cice_to_dart','inq dimlen ni')
call nc_check(nf90_inq_dimid(ncid,'nj',dimid2), &
              'cice_to_dart','inq dimid nj')
call nc_check(nf90_inquire_dimension(ncid,dimid2,len=nj), &
              'cice_to_dart','inq dimlen nj')

allocate (par(npar,ni,nj))


do ipar =1, npar

   call nc_check(nf90_inq_varid(ncid,cice_parameters(ipar),varid), &
                 'cice_to_dart','inquire parameter '//trim(cice_parameters(ipar)))
   call nc_check(nf90_get_var(ncid,varid,par(ipar,:,:)),&
                 'cice_to_dart','get parameter '//trim(cice_parameters(ipar)))
end do

call nc_check(nf90_close(ncid),'cice_to_dart','close'//trim(parameter_input_file))

call nc_check(nf90_open(trim(cice_restart_input_file),NF90_WRITE,ncid2), &
                'cice_to_dart','open '//trim(cice_restart_input_file))

call nc_check(nf90_inq_dimid(ncid2,'ni',dimid1), &
             'cice_to_dart', 'inq dimid ni')
call nc_check(nf90_inq_dimid(ncid2,'nj',dimid2), &
             'cice_to_dart', 'inq dimid nj')


print *,par(1,1,1)
do ipar=1, npar

   ios =  nf90_inq_varid(ncid2,cice_parameters(ipar),varid)
   if (ios/=nf90_noerr) then !variable does not exist
   call nc_check(nf90_redef(ncid2), &
             'cice_to_dart','redef')
   call nc_check(nf90_def_var(ncid2,cice_parameters(ipar),nf90_double,dimids=(/dimid1,dimid2/),varid=varid), &
                'cice to dart', 'def parameter '//trim(cice_parameters(ipar)))
   !print *, cice_parameters(ipar)
   call nc_check(nf90_enddef(ncid2), &
             'cice_to_dart','redef')
   !print *, dimid1,dimid2,size(par(ipar,1,:)),size(par(ipar,:,1))
   end if
   call nc_check(nf90_put_var(ncid2,varid,par(ipar,:,:)),&
                'cice to dart', 'put parameter '//trim(cice_parameters(ipar)))
   
end do

call nc_check(nf90_close(ncid2),'cice to dart','close'//trim(cice_restart_input_file))

call finalize_utilities()

contains

subroutine verify_parameters(parameters, ngood,table)
character (len=*), intent(inout) :: parameters(:)
integer,           intent(out)   :: ngood
character (len=*), intent(out)   :: table(:,:)

integer:: nrows, i
character(len=NF90_MAX_NAME) :: parname

nrows = size(table,1)

ngood = 0

if ( parameters(1) == '') then ! no parameters found in the namelist
   string1 = 'please specify parameters you want to estimate'
   call error_handler(E_MSG,'verify_parameters', string1, source, revision, revdate)
endif

MyLoop : do i = 1, nrows
           parname = trim(parameters(i))
           if(parameters(i) == ' ') exit MyLoop
           ngood = ngood + 1
         end do MyLoop
print *,'ngood:',ngood
end subroutine verify_parameters

end program cice_to_dart

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
