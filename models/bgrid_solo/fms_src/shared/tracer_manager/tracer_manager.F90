!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module tracer_manager_mod
! <CONTACT EMAIL="wfc@gfdl.noaa.gov">
!   William Cooke
! </CONTACT>

! <REVIEWER EMAIL="mjh@gfdl.noaa.gov">
!   Matt Harrison
! </REVIEWER>

! <REVIEWER EMAIL="bw@gfdl.noaa.gov">
!   Bruce Wyman
! </REVIEWER>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

! <OVERVIEW>
!   Code to manage the simple addition of tracers to the FMS code.
!     This code keeps track of the numbers and names of tracers included
!     in a tracer table.
! </OVERVIEW>

! <DESCRIPTION>
!     This code is a grouping of calls which will allow the simple
!     introduction of tracers into the FMS framework. It is designed to
!     allow users of a variety of component models interact easily with
!     the dynamical core of the model. 
!     
!     In calling the tracer manager routines the user must provide a
!     parameter identifying the model that the user is working with. This
!     parameter is defined within field_manager as MODEL_X 
!     where X is one of [ATMOS, OCEAN, LAND, ICE].
!
!     In many of these calls the argument list includes model and tracer_index. These 
!     are the parameter corresponding to the component model and the tracer_index N is 
!     the Nth tracer within the component model. Therefore a call with MODEL_ATMOS and 5 
!     is different from a call with MODEL_OCEAN and 5.
!
! </DESCRIPTION>


!----------------------------------------------------------------------

use types_mod, only : r8
use           fms_mod, only : lowercase,          &
                              write_version_number, error_mesg, &
                              FATAL, WARNING, NOTE, stdlog
use field_manager_mod, only : field_manager_init, &
                              get_field_info,     &
                              get_field_methods,  &
                              MODEL_ATMOS,        &
                              MODEL_LAND,         &
                              MODEL_OCEAN,        &
                              MODEL_ICE,          &
                              NUM_MODELS,         &
                              method_type,        &
                              default_method,     &
                              parse

implicit none
private

!-----------------------------------------------------------------------

public  tracer_manager_end,        &
        check_if_prognostic,       &
        assign_tracer_field,       &
        add_members_to_family,     &
        split_family_into_members, &
        find_family_members,       &
        get_family_name,           &
        get_tracer_indices,        &
        get_tracer_index,          &
        get_tracer_names,          &
        get_tracer_field,          &
        get_tracer_tendency,       &
        get_tracer_tlevels,        &
        query_combined,            &
        query_method,              &
        set_tracer_profile,        &
        register_tracers,          &
        set_tracer_atts,           &
        get_number_tracers

!-----------------------------------------------------------------------

integer            :: num_tracer_fields = 0
integer, parameter :: MAX_TRACER_FIELDS = 30
integer, parameter :: MAX_TRACER_METHOD = 20

integer ::  total_tracers(NUM_MODELS)
integer ::   prog_tracers(NUM_MODELS)
integer ::   diag_tracers(NUM_MODELS)
integer :: family_tracers(NUM_MODELS)

type, private ::  tracer_type
   character(len=32)        :: tracer_name, tracer_units
   character(len=128)       :: tracer_longname, tracer_family
   integer                  :: num_methods, model
   type(method_type)        :: methods(MAX_TRACER_METHOD)
   logical                  :: is_prognostic, is_family, is_combined
   real(r8), pointer, dimension(:,:,:,:) :: field_tlevels
   real(r8), pointer, dimension(:,:,:) :: field, field_tendency, weight
end type tracer_type

type, private ::  tracer_name_type
   character (len=32)    :: model_name, tracer_name, tracer_units
   character (len=128)   :: tracer_longname, tracer_family
end type tracer_name_type

type(tracer_type)  :: tracers(MAX_TRACER_FIELDS)

character(len=128) :: version = '$Revision$'
character(len=128) :: tagname = '$Id$'
logical            :: module_is_initialized = .false.

logical            :: verbose_local
integer            :: TRACER_ARRAY(NUM_MODELS,MAX_TRACER_FIELDS)

contains

!
!#######################################################################
!
! <SUBROUTINE NAME="tracer_manager_init">
!   <OVERVIEW>
!     Routine to initialize the tracer manager
!   </OVERVIEW>
!   <DESCRIPTION>
!     This routine writes the version and tagname to the logfile and 
!     sets the module initialization flag.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call tracer_manager_init
!   </TEMPLATE>
subroutine tracer_manager_init
! version number to logfile

  call write_version_number (version, tagname)
  module_is_initialized = .TRUE.


end subroutine tracer_manager_init
! </SUBROUTINE>
!
!#######################################################################
!
! <SUBROUTINE NAME="register_tracers">

!   <OVERVIEW>
!      A routine to register the tracers included in a component model.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This routine returns the total number of valid tracers, the number of
! prognostic and diagnostic tracers and the number of families of
! tracers.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call register_tracers(model, num_tracers,num_prog,num_diag,num_family)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter to identify which model is being used.
!   </IN>
!   <OUT NAME="num_tracers" TYPE="integer">
!    The total number of valid tracers within the component model.
!   </OUT>
!   <OUT NAME="num_prog" TYPE="integer">
!     The number of prognostic tracers within the component model.
!   </OUT>
!   <OUT NAME="num_diag" TYPE="integer">
!     The number of diagnostic tracers within the component model.
!   </OUT>
!   <OUT NAME="num_family" TYPE="integer">
!     The number of family tracers within the component model.
!   </OUT>

subroutine register_tracers(model, num_tracers,num_prog,num_diag,num_family)
!
! read tracer table and store tracer information associated with "model"
! in "tracers" array. 

integer,  intent(in) :: model ! model being used
integer, intent(out) :: num_tracers, num_prog, num_diag, num_family
character(len=1024)  :: record
character(len=80)    :: warnmesg

character(len=32) :: model_name, name_type, type, name
character(len=128) :: units, longname, family
integer :: iunit,n,m, ntf, mdl, num_tracer_methods, nfields, swop
integer :: j, log_unit, num_in_family, num_methods, ns, num_tracer_comp_model
logical :: flag,is_family_member(MAX_TRACER_FIELDS), flag_type
type(method_type), dimension(20) :: methods

num_tracers = 0; num_prog = 0; num_diag = 0; num_family = 0

call field_manager_init(nfields)
call tracer_manager_init

!   <ERROR MSG="invalid model type" STATUS="FATAL">
!     The index for the model type is invalid.
!   </ERROR>
if (model .ne. MODEL_ATMOS .and. model .ne. MODEL_LAND .and. &
    model .ne. MODEL_OCEAN .and. model .ne. MODEL_ICE) &
    call error_mesg('register_tracers', 'invalid model type', FATAL)


!   <ERROR MSG="No tracers are available to be registered." STATUS="NOTE">
!      No tracers are available to be registered. This means that the field
!      table does not exist or is empty.
!   </ERROR>
if (nfields == 0 ) then
  call error_mesg('register_tracers', 'No tracers are available to be registered.', NOTE)
return
endif

! search through field entries for model tracers
num_tracer_comp_model = 0

do n=1,nfields
   call get_field_info(n,type,name,mdl,num_methods)
   if (mdl == model .and. type == 'tracer') then
      if (get_tracer_index(model, name) < 0) then      
         num_tracer_fields = num_tracer_fields + 1
         num_tracer_comp_model = num_tracer_comp_model + 1
         TRACER_ARRAY(model,num_tracer_comp_model)  = num_tracer_fields
!   <ERROR MSG="MAX_TRACER_FIELDS exceeded" STATUS="FATAL">
!     The maximum number of tracer fields has been exceeded.
!   </ERROR>
         if(num_tracer_fields > MAX_TRACER_FIELDS) &
            call error_mesg('register_tracers', 'MAX_TRACER_FIELDS exceeded', FATAL)
         tracers(num_tracer_fields)%model          = model
         tracers(num_tracer_fields)%tracer_name    = name
         tracers(num_tracer_fields)%tracer_units   = 'none'
         tracers(num_tracer_fields)%tracer_longname = tracers(num_tracer_fields)%tracer_name
         tracers(num_tracer_fields)%tracer_family   = 'orphan'
         num_tracer_methods     = 0
         methods = default_method ! initialize methods array
         call get_field_methods(n,methods)
         do j=1,num_methods
            select case (methods(j)%method_type) 
            case ('units')
               tracers(num_tracer_fields)%tracer_units   = methods(j)%method_name
            case ('longname')
               tracers(num_tracer_fields)%tracer_longname = methods(j)%method_name
            case ('family')
               tracers(num_tracer_fields)%tracer_family = methods(j)%method_name
            case default
               num_tracer_methods = num_tracer_methods+1
               tracers(num_tracer_fields)%methods(num_tracer_methods) = methods(j)
            end select
         enddo
         tracers(num_tracer_fields)%num_methods = num_tracer_methods
         flag_type = query_method ('tracer_type',model,num_tracer_comp_model,name_type)
         if (flag_type .and. name_type == 'diagnostic') then
            tracers(num_tracer_fields)%is_prognostic = .false.
         else   
            tracers(num_tracer_fields)%is_prognostic = .true.
         endif   
         tracers(num_tracer_fields)%is_family = .false.
         tracers(num_tracer_fields)%is_combined = .false.
         if (tracers(num_tracer_fields)%is_prognostic) then
            num_prog = num_prog+1
         else
            num_diag = num_diag+1
         endif
      endif
   endif
enddo

!Now recycle through the tracers and get the family names

ntf = num_tracer_fields
do n=1,ntf
   if (tracers(n)%tracer_family /= 'orphan') then
      call find_family_members(tracers(n)%model,tracers(n)%tracer_family,is_family_member)
      num_in_family = 0
      do m=1,ntf
         if (is_family_member(m)) num_in_family = num_in_family+1
      end do
      if (num_in_family < 2) then ! do not set up a family tracer
         write(warnmesg,911) num_in_family,trim(tracers(n)%tracer_family)
911      format('register_tracers : There is only ',i2,' tracers for tracer family ', (a),'. Making an orphan.')
!   <ERROR MSG="There is only 1 tracer for tracer family X. Making an orphan." STATUS="NOTE">
!     A tracer has been given a family name but that family has only this member. Therefore it should be an orphan.
!   </ERROR>
         call error_mesg('tracer_manager', warnmesg, NOTE)
         tracers(n)%tracer_family = 'orphan'
         cycle
      else
         m = get_tracer_index(tracers(n)%model,tracers(n)%tracer_family)
         if (m < 0) then 
         num_tracer_fields = num_tracer_fields + 1
         num_tracer_comp_model = num_tracer_comp_model + 1
         TRACER_ARRAY(model,num_tracer_comp_model)  = num_tracer_fields
            num_family = num_family+1
!   <ERROR MSG="MAX_TRACER_FIELDS needs to be increased" STATUS="FATL">
!     The number of tracer fields has exceeded the maximum allowed. 
!     The parameter MAX_TRACER_FIELDS needs to be increased.
!   </ERROR>
            if (num_tracer_fields .gt. MAX_TRACER_FIELDS) &
                call error_mesg('register_tracers', 'MAX_TRACER_FIELDS needs to be increased', FATAL)
            write(*,*) 'defining new tracer family: ',trim(tracers(n)%tracer_family)
            tracers(num_tracer_fields)%is_family = .true.
            tracers(num_tracer_fields)%tracer_name = trim(tracers(n)%tracer_family)
            tracers(num_tracer_fields)%model = tracers(n)%model
            tracers(num_tracer_fields)%tracer_units = tracers(n)%tracer_units
            tracers(num_tracer_fields)%tracer_longname = tracers(num_tracer_fields)%tracer_name
            tracers(num_tracer_fields)%tracer_family = 'orphan'
            num_methods = tracers(n)%num_methods
            tracers(num_tracer_fields)%num_methods = num_methods
            do ns=1,num_methods
               tracers(num_tracer_fields)%methods(ns) = tracers(n)%methods(ns)
            end do
            call check_family_parameters(model,tracers(n)%tracer_name) ! make sure family parameters are same for tracers within family
         endif
      endif
   endif
enddo

num_tracers = num_prog + num_diag + num_family
! Make the number of tracers available publicly.
total_tracers(model)  = num_tracers
prog_tracers(model)   = num_prog
diag_tracers(model)   = num_diag
family_tracers(model) = num_family
! Now sort through the tracer fields and sort them so that the 
! prognostic tracers are first. This should include the family tracers.

do n=1, num_tracers
  if (.not.check_if_prognostic(model,n) .and. n.le.(num_prog+num_family)) then 
  ! This is a diagnostic tracer so find a prognostic tracer to swop with
    do m = n, num_tracers
       if (check_if_prognostic(model,m) .and. .not.check_if_prognostic(model,n)) then
           swop = TRACER_ARRAY(model,n)
           TRACER_ARRAY(model,n) = TRACER_ARRAY(model,m)
           TRACER_ARRAY(model,m) = swop
           cycle
       endif
    enddo
  endif
enddo


do n=1, num_tracer_fields
   if(TRACER_ARRAY(model,n)> 0 ) &
      call print_tracer_info(TRACER_ARRAY(model,n))
enddo


log_unit = stdlog()
!   write (log_unit,'(/,80("="),/(a))') trim(version), trim(tagname)
    select case (model)
      case (MODEL_ATMOS)
         model_name = "atmospheric"
      case (MODEL_OCEAN)
         model_name = "oceanic"
      case (MODEL_ICE)
         model_name = "ice"
      case (MODEL_LAND)
         model_name = "land"
      case default
    end select

   write (log_unit,15) trim(model_name),total_tracers(model)

15 format ('Number of tracers in field table for ',A,' model = ',i4)


end subroutine register_tracers
!</SUBROUTINE>


! <SUBROUTINE NAME="get_number_tracers">
!   <OVERVIEW>
!      A routine to return the number of tracers included in a component model.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This routine returns the total number of valid tracers, the number of
! prognostic and diagnostic tracers and the number of families of
! tracers.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call get_number_tracers(model, num_tracers,num_prog,num_diag,num_family)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter to identify which model is being used.
!   </IN>
!   <OUT NAME="num_tracers" TYPE="integer, optional">
!    The total number of valid tracers within the component model.
!   </OUT>
!   <OUT NAME="num_prog" TYPE="integer, optional">
!     The number of prognostic tracers within the component model.
!   </OUT>
!   <OUT NAME="num_diag" TYPE="integer, optional">
!     The number of diagnostic tracers within the component model.
!   </OUT>
!   <OUT NAME="num_family" TYPE="integer, optional">
!     The number of family tracers within the component model.
!   </OUT>
subroutine get_number_tracers(model, num_tracers , num_prog, num_diag, num_family)

integer,  intent(in) :: model
integer, intent(out), optional :: num_tracers, num_prog, num_diag, num_family

!   <ERROR MSG="Model number is invalid." STATUS="FATAL">
!     The index of the component model is invalid.
!   </ERROR>
if (model .ne. MODEL_ATMOS .and. model .ne. MODEL_LAND .and. &
    model .ne. MODEL_OCEAN .and. model .ne. MODEL_ICE)  &
    call error_mesg("get_number_tracers", "Model number is invalid.", FATAL)

if (present(num_tracers)) num_tracers = total_tracers(model)
if (present(num_prog))    num_prog    = prog_tracers(model)
if (present(num_diag))    num_diag    = diag_tracers(model)
if (present(num_family))  num_family  = family_tracers(model)

end subroutine get_number_tracers
!</SUBROUTINE>


! <SUBROUTINE NAME="get_tracer_indices">

!   <OVERVIEW>
!     Routine to return the component model tracer indices as defined within
!     the tracer manager.
!   </OVERVIEW>
!   <DESCRIPTION>
!     If several models are being used or redundant
! tracers have been written to the tracer_table, then the indices in
! the component model and the tracer manager may not have a one to one
! correspondence. Therefore the component model needs to know what index
! to pass to calls to tracer_manager routines in order that the correct
! tracer information be accessed.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call get_tracer_indices(model, ind, prog_ind, diag_ind, fam_ind)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter to identify which model is being used.
!   </IN>
!   <OUT NAME="ind" TYPE="integer, optional" DIM="(:)" >
! An array containing the tracer manager defined indices for
!             all the tracers within the component model.
!   </OUT>
!   <OUT NAME="prog_ind" TYPE="integer, optional" DIM="(:)" >
! An array containing the tracer manager defined indices for
!             the prognostic tracers within the component model.
!   </OUT>
!   <OUT NAME="diag_ind" TYPE="integer, optional" DIM="(:)" >
! An array containing the tracer manager defined indices for
!             the diagnostic tracers within the component model.
!   </OUT>
!   <OUT NAME="fam_ind" TYPE="integer, optional" DIM="(:)" >
! An array containing the tracer manager defined indices for
!             the family tracers within the component model.
!   </OUT>
subroutine get_tracer_indices(model, ind, prog_ind, diag_ind, fam_ind)

integer, intent(in) :: model
integer, intent(out), dimension(:), optional :: ind, prog_ind, diag_ind, fam_ind

integer :: i, j, nf, np, nd, n

nf=0;nd=0;np=0;n=0

do i = 1, MAX_TRACER_FIELDS
j = TRACER_ARRAY(model,i)
 if ( j .gt. 0) then
   if ( model == tracers(j)%model) then
      if (PRESENT(ind)) then
         n=n+1
!   <ERROR MSG="index array size too small in get_tracer_indices" STATUS="Fatal">
!     The global index array is too small and cannot contain all the tracer numbers.
!   </ERROR>
         if (n > size(ind)) &
         call error_mesg('get_tracer_indices', 'index array size too small in get_tracer_indices', FATAL)
         ind(n) = i
      endif

      if (tracers(j)%is_family.and.PRESENT(fam_ind)) then
         nf=nf+1
!   <ERROR MSG="family array size too small in get_tracer_indices" STATUS="FATAL">
!     The family index array is too small and cannot contain all the tracer numbers.
!   </ERROR>
         if (nf > size(fam_ind)) &
         call error_mesg('get_tracer_indices', 'family array size too small in get_tracer_indices', FATAL)
         fam_ind(nf) = i
         cycle
      endif

      if (tracers(j)%is_prognostic.and.PRESENT(prog_ind)) then
         np=np+1
!   <ERROR MSG="prognostic array size too small in get_tracer_indices" STATUS="FATAL">
!     The prognostic index array is too small and cannot contain all the tracer numbers.
!   </ERROR>
         if (np > size(prog_ind)) &
         call error_mesg('get_tracer_indices', 'prognostic array size too small in get_tracer_indices', FATAL)
         prog_ind(np) = i
      else if (.not.tracers(j)%is_prognostic .and. .not. tracers(j)%is_family &
             .and.PRESENT(diag_ind)) then
         nd = nd+1
!   <ERROR MSG="diagnostic array size too small in get_tracer_indices" STATUS="FATAL">
!     The diagnostic index array is too small and cannot contain all the tracer numbers.
!   </ERROR>
         if (nd > size(diag_ind)) &
         call error_mesg('get_tracer_indices', 'diagnostic array size too small in get_tracer_indices', FATAL)
         diag_ind(nd) = i
      endif
   endif
 endif
enddo

return
end subroutine get_tracer_indices
!</SUBROUTINE>

!<FUNCTION NAME= "get_tracer_index">
!   <OVERVIEW>
!     Function which returns the number assigned to the tracer name.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This is a function which returns the index, as implied within the component model.
!   </DESCRIPTION>
!   <TEMPLATE>
!     value=get_tracer_index(model, name, indices, verbose)
!   </TEMPLATE>
!   <IN NAME="model" TYPE="integer">
!     A parameter to identify which model is being used.
!   </IN>
!   <IN NAME="name" TYPE="character">
!     The name of the tracer (as assigned in the field table).
!   </IN>
!   <IN NAME="indices" TYPE="integer, optional" DIM="(:)">
!     An array of the component model indices. This array can be found by 
!            calling get_tracer_indices.
!   </IN>
!   <IN NAME="verbose" TYPE="logical, optional">
!     A flag to allow the message saying that a tracer with this name has not 
!     been found. This should only be used for debugging purposes.
!   </IN>
!   <OUT NAME="get_tracer_index" TYPE="integer">
!     The index of the tracer named "name". 
!     If indices is passed then the result is the array index which 
! corresponds to tracer named "name".
!   </OUT>
function get_tracer_index(model, name, indices, verbose)

integer, intent(in)                         :: model
character(len=*), intent(in)                :: name
integer, intent(in), dimension(:), optional :: indices
logical, intent(in), optional               :: verbose
integer :: get_tracer_index

integer :: i

get_tracer_index = -1

if (PRESENT(indices)) then
    do i = 1, size(indices)
       if (model == tracers(indices(i))%model .and. lowercase(trim(name)) == trim(tracers(indices(i))%tracer_name)) then
           get_tracer_index = i
           exit
       endif
    enddo
else
    do i=1, num_tracer_fields
       if (lowercase(trim(name)) == trim(tracers(TRACER_ARRAY(model,i))%tracer_name)) then
           get_tracer_index = i
           exit
       endif
    enddo
end if



verbose_local=.FALSE.
if (present(verbose)) verbose_local=verbose

if (verbose_local) then
! <ERROR MSG="tracer with this name not found: X" STATUS="NOTE">
  if (get_tracer_index == -1 ) call error_mesg('get_tracer_index', 'tracer with this name not found: '//trim(name), NOTE)
! </ERROR>
endif
   
return

end function get_tracer_index
!</FUNCTION>

! <SUBROUTINE NAME="assign_tracer_field" >
!   <OVERVIEW>
!     Routine to point the appropriate field within the tracer_type to the 
! appropriate field within the component model.
!   </OVERVIEW>
!   <DESCRIPTION>
!     The generality provided here is that one can point the three
! dimensional tracer field at either a two time level scheme [data and
! tendency] or a three time level scheme [data_tlevels]. The tracer manager
! points the appropriate tracer_type field at the data supplied from the 
! component model.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call assign_tracer_field(model,index, data, data_tlevels, tendency)
!   </TEMPLATE>
!   <IN NAME="model" TYPE="integer">
!     A parameter representing the component model in use.
!   </IN>
!   <IN NAME="index" TYPE="integer">
!     The tracer number that you wish to assign a tracer
!                  field for.
!   </IN>
!   <IN NAME="data" TYPE="real(r8), target, optional" DIM="(:,:,:)" >
!     The 3D field that is associated with the present time 
!                  step in the component model.
!   </IN>
!   <IN NAME="tendency" TYPE="real(r8), target, optional" DIM="(:,:,:)" >
!     The 3D field that is associated with the tendency time
!                  step in the component model.
!   </IN>
!   <IN NAME="data_tlevels" TYPE="real(r8), target, optional" DIM="(:,:,:,:)" >
!     The 4D field that is associated with the tracer field 
!                  in the component model.
!   </IN>
subroutine assign_tracer_field(model,index, data, data_tlevels, tendency)

integer, intent(in)                        :: model, index
real(r8), intent(in), dimension(:,:,:), target, optional   :: data, tendency
real(r8), intent(in), dimension(:,:,:,:), target, optional :: data_tlevels

integer :: check

!   <ERROR MSG="invalid index" STATUS="FATAL">
!     The index that has been passed to this routine is invalid.
!   </ERROR>
if (index < 0 .or. index > num_tracer_fields) call error_mesg('assign_tracer_field', 'invalid index', FATAL)

if (PRESENT(data)) tracers(TRACER_ARRAY(model,index))%field => data
if (PRESENT(data_tlevels)) tracers(TRACER_ARRAY(model,index))%field_tlevels => data_tlevels
if (PRESENT(tendency)) tracers(TRACER_ARRAY(model,index))%field_tendency => tendency


check =0
if (PRESENT(data)) check = check + 1
if (PRESENT(data_tlevels)) check = check + 1
if (PRESENT(tendency)) check = check + 1

!   <ERROR MSG="At least one of data, data_tlevels or tendency must be passed in here." STATUS="FATAL">
!     At least one of data, data_tlevels or tendency must be passed to assign_tracer_field
!     Otherwise there is not much point in calling this routine.
!   </ERROR>
if (check == 0) call error_mesg('assign_tracer_field', &
    'At least one of data, data_tlevels or tendency must be passed in here.', FATAL)

return
end subroutine assign_tracer_field
!</SUBROUTINE>

!
!#######################################################################
!

! <SUBROUTINE NAME="tracer_manager_end" >
!   <OVERVIEW>
!     Routine to write to the log file that the tracer manager is ending.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Routine to write to the log file that the tracer manager is ending.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call tracer_manager_end
!   </TEMPLATE>

subroutine tracer_manager_end

integer :: log_unit

log_unit = stdlog()
   write (log_unit,'(/,(a))') 'Exiting tracer_manager, have a nice day ...'

module_is_initialized = .FALSE.

end subroutine tracer_manager_end
!</SUBROUTINE>


!
!#######################################################################
!
subroutine print_tracer_info(i)
!
! Routine to print out the components of the tracer.
! This is useful for informational purposes.
! Used in register_tracers.
!
! Arguments:
! INTENT IN
!  i            : index of the tracer that is being printed.
!
integer, intent(in) :: i
integer :: j,log_unit

character(len=80) :: name,control
logical   :: flag


log_unit = stdlog()
write(log_unit, *)'----------------------------------------------------'
write(log_unit, *) 'Contents of tracer entry ', i
write(log_unit, *) 'Model type and field name'
write(log_unit, *) 'Model                : ', tracers(i)%model
write(log_unit, *) 'Field name           : ', trim(tracers(i)%tracer_name)
write(log_unit, *) 'Tracer family        : ', trim(tracers(i)%tracer_family)
write(log_unit, *) 'Tracer units         : ', trim(tracers(i)%tracer_units)
write(log_unit, *) 'Tracer longname      : ', trim(tracers(i)%tracer_longname)
write(log_unit, *) 'Tracer is_prognostic : ', tracers(i)%is_prognostic
write(log_unit, *) 'Tracer is_family     : ', tracers(i)%is_family
write(log_unit, *) 'Tracer is_combined   : ', tracers(i)%is_combined

if (associated(tracers(i)%field)) then
   write(log_unit, '(a,3i5)') 'Size of tracer field : ', size(tracers(i)%field,1), &
                                                         size(tracers(i)%field,2), &
                                                         size(tracers(i)%field,3)
else
   write(log_unit,*) 'Tracer field         : not associated'
endif

if (associated(tracers(i)%field_tendency)) then
   write(log_unit, '(a,3i5)') 'Size of tracer tendency: ', size(tracers(i)%field_tendency,1), &
                                                           size(tracers(i)%field_tendency,2), &
                                                           size(tracers(i)%field_tendency,3)
else
   write(log_unit,*) 'Tracer tendency      : not associated'
endif

if (associated(tracers(i)%field_tlevels)) then
   write(log_unit, '(a,3i5)') 'Size of tracer tlevels : ', size(tracers(i)%field_tlevels,1), &
                                                           size(tracers(i)%field_tlevels,2), &
                                                           size(tracers(i)%field_tlevels,3), &
                                                           size(tracers(i)%field_tlevels,4)
else
   write(log_unit,*) 'Tracer tlevels       : not associated'
endif



do j=1,tracers(i)%num_methods
   write(log_unit, '(a,i4,a,a)') ' Method ',j,' type     : ',trim(tracers(i)%methods(j)%method_type)
   write(log_unit, '(a,i4,a,a)') ' Method ',j,' name     : ',trim(tracers(i)%methods(j)%method_name)
   write(log_unit, '(a,i4,a,a)') ' Method ',j,' control  : ',trim(tracers(i)%methods(j)%method_control)
enddo





write(log_unit, *)'----------------------------------------------------'

900 FORMAT(A,2(1x,E14.6))   ! are these even used
901 FORMAT(E14.6,1x,E14.6)  ! are these even used


end subroutine print_tracer_info

!
!#######################################################################
!
!<FUNCTION NAME= "get_tracer_field">
!   <OVERVIEW>
!     A function to retrieve the present timestep data.
!   </OVERVIEW>
!   <DESCRIPTION>
!      Function to point to the 3D field associated with a tracer.
!   </DESCRIPTION>
!   <TEMPLATE>
!     array=get_tracer_field(model, tracer_index)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter representing the component model in use.
!   </IN>
!   <IN NAME="tracer_index" TYPE="integer">
!     The tracer number within the component model.
!   </IN>
!   <OUT NAME="data"  TYPE="real(r8), pointer" DIM="(:,:,:)">
!     The tracer field is returned in this array.
!   </OUT>
function get_tracer_field(model, tracer_index) result (data)

integer              :: model, tracer_index
real(r8), pointer        :: data(:,:,:)

integer :: n

!   <ERROR MSG="invalid index" STATUS="FATAL">
!     The index that has been passed to this routine is invalid.
!          Check the index that is being passed corresponds to a valid
!          tracer name.
!   </ERROR>
if (tracer_index < 1 .or. tracer_index > num_tracer_fields) &
    call error_mesg('get_tracer_field', 'invalid index ', FATAL)
!Convert local model index to tracer_manager index
!   <ERROR MSG="invalid index" STATUS="FATAL">
!     The index that has been passed to this routine is invalid.
!          Check the index that is being passed corresponds to a valid
!          tracer name.
!   </ERROR>
if (TRACER_ARRAY(model,tracer_index) < 1 .or. TRACER_ARRAY(model,tracer_index) > num_tracer_fields) &
    call error_mesg('get_tracer_field', 'invalid index', FATAL)
!   <ERROR MSG="tracer field array not allocated" STATUS="FATAL">
!      The tracer array has not been allocated. This means that a
!          call to assign_tracer_field is absent in the code.
!   </ERROR>
if (.not. associated(tracers(TRACER_ARRAY(model,tracer_index))%field)) &
    call error_mesg('get_tracer_field', 'tracer field array not allocated', FATAL)
data =>  tracers(TRACER_ARRAY(model,tracer_index))%field

end function get_tracer_field
!</FUNCTION>
!
!#######################################################################
!
!<FUNCTION NAME= "get_tracer_tlevels">
!   <OVERVIEW>
!     A function to retrieve the three time levels  data.
!   </OVERVIEW>
!   <DESCRIPTION>
!      Function to point to the 4D field associated with a tracer.
!   </DESCRIPTION>
!   <TEMPLATE>
!     array=get_tracer_tlevels(model, tracer_index)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter representing the component model in use.
!   </IN>
!   <IN NAME="tracer_index" TYPE="integer">
!     The tracer number within the component model.
!   </IN>
!   <OUT NAME="data"  TYPE="real(r8), pointer" DIM="(:,:,:,:)">
!     The tracer field is returned in this array.
!   </OUT>
function get_tracer_tlevels(model, tracer_index) result (data)

integer  :: model, tracer_index
real(r8), pointer        :: data(:,:,:,:)

integer :: n

!   <ERROR MSG="invalid index" STATUS="FATAL">
!     The index that has been passed to this routine is invalid.
!          Check the index that is being passed corresponds to a valid
!          tracer name.
!   </ERROR>
if (tracer_index < 1 .or. tracer_index > num_tracer_fields) &
    call error_mesg('get_tracer_tlevels', 'invalid index', FATAL)
!Convert local model index to tracer_manager index
!   <ERROR MSG="invalid index" STATUS="FATAL">
!     The index that has been passed to this routine is invalid.
!          Check the index that is being passed corresponds to a valid
!          tracer name.
!   </ERROR>
if (TRACER_ARRAY(model,tracer_index) < 1 .or. TRACER_ARRAY(model,tracer_index) > num_tracer_fields) &
    call error_mesg('get_tracer_tlevels', 'invalid index', FATAL)
!   <ERROR MSG="tracer field array not allocated" STATUS="FATAL">
!      The tracer array has not been allocated. This means that a
!          call to assign_tracer_field is absent in the code.
!   </ERROR>
if (.not. associated(tracers(TRACER_ARRAY(model,tracer_index))%field_tlevels)) &
    call error_mesg('get_tracer_tlevels', 'tracer field array not allocated', FATAL)
data =>  tracers(TRACER_ARRAY(model,tracer_index))%field_tlevels

end function get_tracer_tlevels
!</FUNCTION>
!
!#######################################################################
!
!<FUNCTION NAME= "get_tracer_tendency">
!   <OVERVIEW>
!     A function to retrieve the tendency data.
!   </OVERVIEW>
!   <DESCRIPTION>
!      Function to point to the 3D field associated with a tracer.
!   </DESCRIPTION>
!   <TEMPLATE>
!     array=get_tracer_tendency(model, tracer_index)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter representing the component model in use.
!   </IN>
!   <IN NAME="tracer_index" TYPE="integer">
!     The tracer number within the component model.
!   </IN>
!   <OUT NAME="data"  TYPE="real(r8), pointer" DIM="(:,:,:)">
!     The tracer tendency field is returned in this array.
!   </OUT>
function get_tracer_tendency(model, tracer_index) result (data)

integer, intent(in)  :: model
integer :: tracer_index
real(r8), pointer        :: data(:,:,:)

integer :: n

!   <ERROR MSG="invalid index" STATUS="FATAL">
!     The index that has been passed to this routine is invalid.
!          Check the index that is being passed corresponds to a valid
!          tracer name.
!   </ERROR>
if (tracer_index < 1 .or. tracer_index > num_tracer_fields) &
    call error_mesg('get_tracer_tendency', 'invalid index', FATAL)
!Convert local model index to tracer_manager index
!   <ERROR MSG="invalid index" STATUS="FATAL">
!     The index that has been passed to this routine is invalid.
!          Check the index that is being passed corresponds to a valid
!          tracer name.
!   </ERROR>
if (TRACER_ARRAY(model,tracer_index) < 1 .or. TRACER_ARRAY(model,tracer_index) > num_tracer_fields) &
    call error_mesg('get_tracer_tendency', 'invalid index', FATAL)
!   <ERROR MSG="tracer tendency field array not allocated" STATUS="FATAL">
!      The tracer array has not been allocated. This means that a
!          call to assign_tracer_field is absent in the code.
!   </ERROR>
if (.not. associated(tracers(TRACER_ARRAY(model,tracer_index))%field_tendency)) &
    call error_mesg('get_tracer_tendency', 'tracer tendency field array not allocated', FATAL)
data =>  tracers(TRACER_ARRAY(model,tracer_index))%field_tendency

end function get_tracer_tendency
!</FUNCTION>

!#######################################################################
!
! <SUBROUTINE NAME="get_tracer_names" >
!   <OVERVIEW>
!     Routine to find the names associated with a tracer number.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This routine can return the name, long name and units associated
!     with a tracer.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call get_tracer_names(model,n,name,longname, units)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter representing the component model in use.
!   </IN>
!   <IN NAME="n" TYPE="integer">
!     Tracer number.
!   </IN>
!   <OUT NAME="name" TYPE="character" >
!     Field name associated with tracer number.
!   </OUT>
!   <OUT NAME="longname" TYPE="character, optional" >
!     The long name associated with tracer number.
!   </OUT>
!   <OUT NAME="units" TYPE="character, optional" >
!     The units associated with tracer number.
!   </OUT>

subroutine get_tracer_names(model,n,name,longname, units)

integer,          intent(in)  :: model, n
character (len=*),intent(out) :: name
character (len=*), intent(out), optional :: longname, units

if (n < 1 .or. n > num_tracer_fields) &
    call error_mesg('get_tracer_names', 'invalid local tracer index for '//trim(name), FATAL)
!Convert local model index to tracer_manager index
if (TRACER_ARRAY(model,n) < 1 .or. TRACER_ARRAY(model,n) > num_tracer_fields) &
    call error_mesg('get_tracer_names', 'invalid tracer index for '//trim(name), FATAL)

name=tracers(TRACER_ARRAY(model,n))%tracer_name
if (PRESENT(longname)) longname =tracers(TRACER_ARRAY(model,n))%tracer_longname
if (PRESENT(units)) units =tracers(TRACER_ARRAY(model,n))%tracer_units

end subroutine get_tracer_names
!</SUBROUTINE>
!
!#######################################################################
!
! <SUBROUTINE NAME="get_family_name" >
!   <OVERVIEW>
!     Routine to return the family name for tracer n.
!   </OVERVIEW>
!   <DESCRIPTION>
!     You may wish to use this routine to retrieve the name of the family
!     that a tracer belongs to.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call get_family_name(model,n,name)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter representing the component model in use.
!   </IN>
!   <IN NAME="n" TYPE="integer">
!     Tracer number that you want the family name for.
!   </IN>
!   <OUT NAME="name" TYPE="character" >
!     The family name.
!   </OUT>
subroutine get_family_name(model,n,name)

integer,          intent(in)  :: model, n
character (len=*),intent(out) :: name

if (n < 1 .or. n > num_tracer_fields) call error_mesg('get_family_name', 'invalid tracer index', FATAL)
!Convert local model index to tracer_manager index
name=tracers(TRACER_ARRAY(model,n))%tracer_family

end subroutine get_family_name
!</SUBROUTINE>
!
!#######################################################################
!
!<FUNCTION NAME= "check_if_prognostic">
!   <OVERVIEW>
!    Function to see if a tracer is prognostic or diagnostic.
!   </OVERVIEW>
!   <DESCRIPTION>
!    All tracers are assumed to be prognostic when read in from the field_table
!    However a tracer can be changed to a diagnostic tracer by adding the line
!    "tracer_type","diagnostic"
!    to the tracer description in field_table.
!   </DESCRIPTION>
!   <TEMPLATE>
!     logical =check_if_prognostic(model, n)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter representing the component model in use.
!   </IN>
!   <IN NAME="n" TYPE="integer">
!     Tracer number that you want the family name for.
!   </IN>
!   <OUT NAME="check_if_prognostic" TYPE="logical">
!     A logical flag set TRUE if the tracer is 
!                        prognostic.
!   </OUT>
function check_if_prognostic(model, n)

integer, intent(in) :: model, n
logical             :: check_if_prognostic

if (n < 1 .or. n > num_tracer_fields) call error_mesg('check_if_prognostic', 'invalid tracer index', FATAL)
!Convert local model index to tracer_manager index

check_if_prognostic = tracers(TRACER_ARRAY(model,n))%is_prognostic

end function check_if_prognostic
!</FUNCTION>
!
!#######################################################################
!
! <SUBROUTINE NAME="find_family_members" >
!   <OVERVIEW>
!     Subroutine to find which tracers are members of family family_name.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Subroutine to find which tracers are members of family family_name.
! This will return a logical array where the array positions 
! corresponding to the tracer numbers for family members are set .TRUE.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call find_family_members(model, family_name,is_family_member)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter representing the component model in use.
!   </IN>
!   <IN NAME="family_name" TYPE="character">
!      The family name of the members one is seeking.
!   </IN>
!   <OUT NAME="is_family_member" TYPE="logical" DIM ="(:)">
!      A logical array where the tracer number is used as 
!                     the index to signify which tracer is part of the family.
!                     i.e. If tracers 1, 3, and 7 are part of the same family
!                     then is_family_member(1), is_family_member(3), and 
!                     is_family_member(7) are set TRUE.
!   </OUT>
subroutine find_family_members(model, family_name,is_family_member)

integer, intent(in)         :: model
character(len=*),intent(in) :: family_name
logical , intent(out)       :: is_family_member(:)
integer ::n

is_family_member = .false.

if (size(is_family_member) < num_tracer_fields) &
    call error_mesg('find_family_members', 'array too short', FATAL)

do n=1,num_tracer_fields
   if(trim(tracers(TRACER_ARRAY(model,n))%tracer_family) == trim(family_name) .and. &
      tracers(TRACER_ARRAY(model,n))%model == model ) then
      is_family_member(n)= .true.
   endif
enddo

end subroutine find_family_members
!</SUBROUTINE>

!
!#######################################################################
!
subroutine check_family_parameters(model, family_name)
!
! Subroutine to check that the advection and diffusion schemes and 
! controls of each tracer in a family is the same.
! INTENT IN
!  model       : The model that you are calling this subroutine from.
!  family_name : The name that has been given to the family of tracers of interest
!
character(len=*),intent(in) :: family_name
integer, intent(in)         :: model
logical            :: is_family_member(num_tracer_fields), flag
integer            :: nfam,n
character(len=128) :: name,advectvert, advecthoriz, diffusvert, diffushoriz

nfam = get_tracer_index(model, family_name)
flag=.true.

if (nfam > 0) then 
   if( query_method ('advection_scheme_vert',model,nfam,name)) advectvert = name
   if( query_method ('advection_scheme_horiz',model,nfam,name)) advecthoriz = name
   if( query_method ('diffusion_scheme_vert',model,nfam,name)) diffusvert = name
   if( query_method ('diffusion_scheme_horiz',model,nfam,name)) diffushoriz = name
   call find_family_members(model,family_name,is_family_member)
   do n=1,num_tracer_fields
      if(is_family_member(n)) then
         if( query_method ('advection_scheme_vert',model,n,name))  then
            if (name .ne. advectvert)  flag=.FALSE.
         endif
         if( query_method ('advection_scheme_horiz',model,n,name)) then
            if (name .ne. advecthoriz)  flag=.FALSE.
         endif
         if( query_method ('diffusion_scheme_vert',model,n,name))  then
            if (name .ne. diffusvert)  flag=.FALSE.
         endif
         if( query_method ('diffusion_scheme_horiz',model,n,name)) then
            if (name .ne. diffushoriz)  flag=.FALSE.
         endif
      endif
      if(.not. flag) call error_mesg('check_family_parameters', &
           'Family members do not have the same parameters for advection and diffusion', FATAL)
   enddo
endif

end subroutine check_family_parameters
!
!#######################################################################
!
! <SUBROUTINE NAME="add_members_to_family" >
!   <OVERVIEW>
!     Routine to sum up the members of a family of tracers so that they may
! be advected and diffused as one tracer.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Routine to sum up the members of a family of tracers so that they may
! be advected and diffused as one tracer. This should only be used in
! conjunction with split_family_into_members and should be placed before
! the advection scheme is called. 
!   </DESCRIPTION>
!   <TEMPLATE>
!     call add_members_to_family(model,family_name, cur, prev, next)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter representing the component model in use.
!   </IN>
!   <IN NAME="n" TYPE="integer">
!     Tracer number.
!   </IN>
!   <IN NAME="cur" TYPE="integer, optional">
!     Array index for the current time step. This is only of use 
!                with a three timestep model.
!   </IN>
!   <IN NAME="prev" TYPE="integer, optional">
!     Array index for the previous time step. This is only of use 
!                with a three timestep model.
!   </IN>
!   <IN NAME="next" TYPE="integer, optional">
!     Array index for the next time step. This is only of use 
!                with a three timestep model.
!   </IN>

!   <NOTE>
! This should be used with extreme caution. 
! Unless the family member distributions are similar to each other spatially, 
! advection as one tracer and subsequent splitting will result in a different
! result to advecting each tracer separately. The user should understand the 
! possible repercussions of this before using it.
!</NOTE>
subroutine add_members_to_family(model,family_name, cur, prev, next)

integer, intent(in)                     :: model
character(len=*), intent(in)            :: family_name
integer, intent(in), optional           :: cur, prev, next

!
logical :: is_family_member(num_tracer_fields)
integer :: nfam,n, min_indx, t0, tm1, tp1, siz(num_tracer_fields,3)

nfam=0
nfam = get_tracer_index(model,family_name)
!This is the family number within whichever model
!is calling this routine. So need to convert it to 
! the tracer_manager tracer number.
nfam = TRACER_ARRAY(model,nfam)
! If the tracer is not a family it doesn't need to be summed.
if (.not.tracers(nfam)%is_family) return


if (nfam > 0) then
   call find_family_members(model, family_name,is_family_member)    
   if (associated(tracers(nfam)%field).and.associated(tracers(nfam)%field_tendency)) then 
       tracers(nfam)%field = 0.0 ! initialize
       do n=1,num_tracer_fields
          if (is_family_member(n)) then
              if (.not.associated(tracers(n)%field)) &
                  call error_mesg('add_members_to_family', 'current tracer field not associated', FATAL)
              siz(n,1) = size(tracers(n)%field,1)
              siz(n,2) = size(tracers(n)%field,2)
              siz(n,3) = size(tracers(n)%field,3)
              tracers(nfam)%field=tracers(nfam)%field+tracers(n)%field 
          endif
       enddo
! Now divide the members by the family total
! The members should then add up to 1.
       do n=1,num_tracer_fields
          if (is_family_member(n)) then
             if (.not.associated(tracers(n)%weight)) allocate(tracers(n)%weight(siz(n,1),siz(n,2),siz(n,3)))
              where(tracers(nfam)%field /= 0.0)
                  tracers(n)%weight = tracers(n)%field/tracers(nfam)%field
              elsewhere
                  tracers(n)%weight = 0.0
              end where
              tracers(n)%is_combined = .true.
          endif
       enddo
    else if (associated(tracers(nfam)%field_tlevels)) then
       if (.not.PRESENT(cur).or..not.PRESENT(prev).or..not.PRESENT(next)) then
           call error_mesg('add_members_to_family', &
               'need to specify time level indices to add family members', FATAL)
       endif
       ! pointer array indices start at 1
       min_indx = min(cur,prev,next)
       t0 = cur-min_indx+1
       tm1 = prev-min_indx+1
       tp1 = next-min_indx+1
       tracers(nfam)%field_tlevels(:,:,:,:)=0.0 ! initialize
       do n=1,num_tracer_fields
          if (is_family_member(n)) then
             if (.not.associated(tracers(n)%field_tlevels)) &
                 call error_mesg('add_members_to_family', 'tracer time levels not associated', FATAL)
             tracers(nfam)%field_tlevels(:,:,:,t0)=tracers(nfam)%field_tlevels(:,:,:,t0)+&
                  tracers(n)%field_tlevels(:,:,:,t0)
             tracers(nfam)%field_tlevels(:,:,:,tm1)=tracers(nfam)%field_tlevels(:,:,:,tm1)+&
                  tracers(n)%field_tlevels(:,:,:,tm1)
             siz(n,1) = size(tracers(n)%field_tlevels,1)
             siz(n,2) = size(tracers(n)%field_tlevels,2)
             siz(n,3) = size(tracers(n)%field_tlevels,3)
          endif
       enddo
! Now divide the members by the family total
! The members should then add up to 1.
       do n=1,num_tracer_fields
          if (is_family_member(n)) then
              if (.not.associated(tracers(n)%weight)) allocate(tracers(n)%weight(siz(n,1),siz(n,2),siz(n,3)))
              where(tracers(nfam)%field_tlevels(:,:,:,t0) /= 0.0)  
                  tracers(n)%weight(:,:,:) = tracers(n)%field_tlevels(:,:,:,t0)/&
                                                tracers(nfam)%field_tlevels(:,:,:,t0)
              elsewhere
                  tracers(n)%weight = 0.0
              end where
              tracers(n)%is_combined = .true.
          endif
       enddo
    else
       call error_mesg('add_members_to_family', 'need to associate tracer pointers to use families', FATAL)
    endif
 endif

 return

end subroutine add_members_to_family
!</SUBROUTINE>

!
!#######################################################################
!
! Subroutine that sets the present value of the member of a tracer 
! family according to the fraction of the family that it was in the 
! previous step.
!
! <SUBROUTINE NAME="split_family_into_members" >
!   <OVERVIEW>
!      Subroutine that sets the present value of the member of a tracer 
! family according to the fraction of the family that it was in the 
! previous step.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Subroutine that sets the present value of the member of a tracer 
! family according to the fraction of the family that it was in the 
! previous step.
!
! This splits the transported family into the constituent members. This
! should only be used in conjunction with <I>add_members_to_family</I> and
! should be placed after the advection scheme is called.
!
!   </DESCRIPTION>
!   <TEMPLATE>
!     call split_family_into_members(model,family_name,cur,prev,next)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter representing the component model in use.
!   </IN>
!   <IN NAME="family_name" TYPE="character">
!     The name of the family of tracers that you would 
!                like to split up.
!   </IN>
!   <IN NAME="cur" TYPE="integer, optional">
!     Array index for the current time step. This is only of use 
!                with a three timestep model.
!   </IN>
!   <IN NAME="prev" TYPE="integer, optional">
!     Array index for the previous time step. This is only of use 
!                with a three timestep model.
!   </IN>
!   <IN NAME="next" TYPE="integer, optional">
!     Array index for the next time step. This is only of use 
!                with a three timestep model.
!   </IN>

!   <NOTE>
! This should be used with extreme caution. 
! Unless the family member distributions are similar to each other spatially, 
! advection as one tracer and subsequent splitting will result in a different
! result to advecting each tracer separately. The user should understand the 
! possible repercussions of this before using it.
!</NOTE>
subroutine split_family_into_members(model,family_name,cur,prev,next)

integer, intent(in) :: model
character(len=*),intent(in) :: family_name
integer, intent(in), optional :: cur, prev, next

logical :: is_family_member(num_tracer_fields)
integer :: nfam,n, min_indx, t0, tm1, tp1

! need to make sure tracers have been combined already

nfam=0
nfam = get_tracer_index(model,family_name)
! This is the family tracer number within whichever model
! is calling this routine. So need to convert it to 
! the tracer_manager tracer number.
nfam = TRACER_ARRAY(model,nfam)

! If the family name is not the same as the tracer name then the
! tracer is not a family tracer (It may be a member of a family)
if (tracers(nfam)%tracer_name .ne. tracers(nfam)%tracer_family) return

write(*,*) 'split',model, trim(family_name),nfam,TRACER_ARRAY(model,nfam)

if (nfam > 0) then 
   call find_family_members(model,family_name,is_family_member)
   if (associated(tracers(nfam)%field).and.associated(tracers(nfam)%field_tendency)) then 
      do n=1,num_tracer_fields
         if (is_family_member(n)) then           
            if (.not. tracers(n)%is_combined) &
                call error_mesg('split_family_into_members', &
                   'call to split family into members when fields are not combined', FATAL)
            tracers(n)%field = tracers(n)%weight*tracers(nfam)%field 
            tracers(n)%is_combined = .false.
         endif
      enddo
   else if (associated(tracers(nfam)%field_tlevels)) then
      if (.not.PRESENT(cur).or..not.PRESENT(prev).or..not.PRESENT(next)) then
         call error_mesg('split_family_into_members', &
                         'need to specify time level indices to split family members', FATAL)
      endif
      ! pointer array indices start at 1
      min_indx = min(cur,prev,next)
      t0 = cur-min_indx+1
      tm1 = prev-min_indx+1
      tp1 = next-min_indx+1
      do n=1,num_tracer_fields
         if (is_family_member(n)) then
            if (.not. tracers(n)%is_combined) call error_mesg('split_family_into_members', &
                         'call to split family into members when fields are not combined', FATAL)
            tracers(n)%field_tlevels(:,:,:,tp1) = tracers(nfam)%field_tlevels(:,:,:,tp1)*tracers(n)%weight(:,:,:)
            tracers(n)%is_combined = .false.
         endif
      enddo
   else
      call error_mesg('split_family_into_members', 'need to associate tracer pointers to use families', FATAL)
   endif
endif


end subroutine split_family_into_members
!</SUBROUTINE>

!
!#######################################################################
!
! <SUBROUTINE NAME="set_tracer_profile" >
!   <OVERVIEW>
!     Subroutine to set the tracer field to the wanted profile.
!   </OVERVIEW>
!   <DESCRIPTION>
!     If the profile type is 'fixed' then the tracer field values are set 
! equal to the surface value.
! If the profile type is 'profile' then the top/bottom of model and
! surface values are read and an exponential profile is calculated,
! with the profile being dependent on the number of levels in the
! component model. This should be called from the part of the dynamical
! core where tracer restarts are called in the event that a tracer
! restart file does not exist.
!
!  This can be activated by adding a method to the field_table
! e.g.
!  "profile_type","fixed","surface_value = 1e-12"
!  would return values of surf_value = 1e-12 and a multiplier of 1.0
!  One can use these to initialize the entire field with a value of 1e-12.
!
!  "profile_type","profile","surface_value = 1e-12, top_value = 1e-15"
!   In a 15 layer model this would return values of surf_value = 1e-12 and 
!   multiplier = 0.6309573 i.e 1e-15 = 1e-12*(0.6309573^15)
!   In this case the model should be MODEL_ATMOS as you have a "top" value.
!
!   If you wish to initialize the ocean model, one can use bottom_value instead
!   of top_value.

!   </DESCRIPTION>
!   <TEMPLATE>
!     call set_tracer_profile(model, n, surf_value, multiplier)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter representing the component model in use.
!   </IN>
!   <IN NAME="n" TYPE="integer">
!     Tracer number.
!   </IN>
!   <OUT NAME="surf_value" TYPE="real(r8)">
!     The surface value that will be initialized for the tracer
!   </OUT>
!   <OUT NAME="multiplier" TYPE="real(r8)">
!     The vertical multiplier for the tracer
!                   Level(k-1) = multiplier * Level(k)
!   </OUT>
subroutine set_tracer_profile(model, n, surf_value, multiplier)

integer,  intent(in)  :: model, n
real(r8),     intent(out) :: surf_value, multiplier

integer :: numlevels,m
real(r8) :: top_value, bottom_value
character(len=80) :: scheme, control,profile_type
integer :: flag

!default values
profile_type  = 'Fixed'
surf_value = 0.0E+00
top_value  = surf_value
bottom_value = surf_value
multiplier = 1.0

if ( query_method ( 'profile_type',model,n,scheme,control)) then
!Change the tracer_number to the tracer_manager version

  if(lowercase(trim(scheme(1:5))).eq.'fixed') then
    profile_type                   = 'Fixed'
    flag =parse(control,'surface_value',surf_value)
    multiplier = 1.0
  endif

  if(lowercase(trim(scheme(1:7))).eq.'profile') then
    profile_type                   = 'Profile'
    flag=parse(control,'surface_value',surf_value)
    write(*,*) 'surf ',surf_value, control,flag
    if (surf_value .eq. 0.0) &
      call error_mesg('set_tracer_profile', 'Cannot have a zero surface value for an exponential profile', FATAL)
    select case (tracers(TRACER_ARRAY(model,n))%model)
      case (MODEL_ATMOS)
        flag=parse(control,'top_value',top_value)
        if(flag == 0) &
           call error_mesg('set_tracer_profile', 'Parameter top_value needs to be defined for the tracer profile.', NOTE)
      case (MODEL_OCEAN)
        flag =parse(control,'bottom_value',bottom_value)
        if(flag == 0) &
           call error_mesg('set_tracer_profile', 'Parameter bottom_value needs to be defined for the tracer profile.', NOTE)
      case default
    end select

! If profile type is profile then set the surface value to the input
! value and calculate the vertical multiplier.
! 
! Assume an exponential decay/increase from the surface to the top level
!  C = C0 exp ( -multiplier* level_number)
!  => multiplier = exp [ ln(Ctop/Csurf)/number_of_levels]
!
    if (associated(tracers(TRACER_ARRAY(model,n))%field)) numlevels = size(tracers(TRACER_ARRAY(model,n))%field,3) -1
    if (associated(tracers(TRACER_ARRAY(model,n))%field_tlevels)) &
                 numlevels = size(tracers(TRACER_ARRAY(model,n))%field_tlevels,3)-1
    select case (tracers(TRACER_ARRAY(model,n))%model)
      case (MODEL_ATMOS)
        multiplier = exp( log (top_value/surf_value) /numlevels)
      case (MODEL_OCEAN)
        multiplier = exp( log (bottom_value/surf_value) /numlevels)
      case default
    end select
  endif !scheme.eq.profile

  write(*,700) 'Tracer ',trim(tracers(TRACER_ARRAY(model,n))%tracer_name),    &
                            ' initialized with surface value of ',surf_value, &
                            ' and vertical multiplier of ',multiplier
  700 FORMAT (3A,E12.6,A,F10.6)

endif ! end of query scheme

end subroutine set_tracer_profile
!</SUBROUTINE>

!
!#######################################################################
!
! <FUNCTION NAME="query_method" >
!   <OVERVIEW>
!     A function to query the "methods" associated with each tracer.
!   </OVERVIEW>
!   <DESCRIPTION>
!     A function to query the "methods" associated with each tracer. The
!  "methods" are the parameters of the component model that can be
!  adjusted by user by placing formatted strings, associated with a
!  particular tracer, within the field table.
!  These methods can control the advection, wet deposition, dry
!  deposition or initial profile of the tracer in question. Any
!  parametrization can use this function as long as a routine for parsing
!  the name and control strings are provided by that routine.
!   </DESCRIPTION>
!   <TEMPLATE>
!     logical =query_method  (method_type, model, n, name, control)
!   </TEMPLATE>

!   <IN NAME="method_type" TYPE="character">
!     The method that is being requested.
!   </IN>
!   <IN NAME="model" TYPE="integer">
!     A parameter representing the component model in use.
!   </IN>
!   <IN NAME="n" TYPE="integer">
!     Tracer number that you want the family name for.
!   </IN>
!   <OUT NAME="name" TYPE="character">
!     A string containing the modified name to be used with
!     method_type. i.e. "2nd_order" might be the default for 
!     advection. One could use "4th_order" here to modify 
!     that behaviour.
!   </OUT>
!   <OUT NAME="control" TYPE="character, optional">
!     A string containing the modified parameters that are 
!     associated with the method_type and name.
!   </OUT>
!   <OUT NAME="query_method" TYPE="logical">
!      A flag to show whether method_type exists with regard to
!      tracer n. If method_type is not present then one must
!      have default values.
!   </OUT>

!<NOTE>
!  At present the tracer manager module allows the initialization of a tracer
!  profile if a restart does not exist for that tracer. 
!  Options for this routine are as follows
!
!  Tracer profile setup
!  ==================================================================
!  |method_type  |method_name  |method_control                      |
!  ==================================================================
!  |profile_type |fixed        |surface_value = X                   |
!  |profile_type |profile      |surface_value = X, top_value = Y    |(atmosphere)
!  |profile_type |profile      |surface_value = X, bottom_value = Y |(ocean)
!  ==================================================================
!
!</NOTE>
function query_method  (method_type, model, n, name, control)
!
!  A function to query the schemes associated with each tracer. 
!  
!  INTENT IN
!   method_type  : The method that is being requested.
!   model        : The model that you are calling this function from.
!   n            : The tracer number.
!  INTENT OUT
!   name         : A string containing the modified name to be used with
!                  method_type. i.e. "2nd_order" might be the default for 
!                  advection. One could use "4th_order" here to modify 
!                  that behaviour.
!   control      : A string containing the modified parameters that are 
!                  associated with the method_type and name.
!   query_method : A flag to show whether method_type exists with regard 
!                  to tracer n. If method_type is not present then one
!                  must have default values.

character(len=*), intent(in)            :: method_type
integer         , intent(in)            :: model, n
character(len=*), intent(out)           :: name
character(len=*), intent(out), optional :: control
logical                                 :: query_method

integer :: m, n1

!Convert the local model tracer number to the tracer_manager version.

name=" "
query_method = .false.
n1 = TRACER_ARRAY(model,n)
do m = 1,tracers(n1)%num_methods
   if (tracers(n1)%methods(m)%method_type .eq. lowercase(method_type) ) then
      name         = tracers(n1)%methods(m)%method_name
      if (PRESENT(control)) control = tracers(n1)%methods(m)%method_control
      query_method = .true.
   endif
enddo

end function query_method
!</FUNCTION>

! <FUNCTION NAME="query_combined" >
!   <OVERVIEW>
!     A function to query whether families of tracers have been combined already.
!   </OVERVIEW>
!   <DESCRIPTION>
!     A function to query whether families of tracers have been combined already.
!  This function should only be used in conjunction with add_members_to_family 
!  and split_family_into_members.
!   </DESCRIPTION>
!   <TEMPLATE>
!     logical =query_combined (model, index)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter representing the component model in use.
!   </IN>
!   <IN NAME="index" TYPE="integer">
!     Tracer number.
!   </IN>
!   <OUT NAME="query_combined" TYPE="logical" >
!     A flag to show whether the tracer family has been combined.
!   </OUT>
function query_combined (model, index) 

integer,intent(in) :: model, index
logical :: query_combined

if (index < 0 .or. index > num_tracer_fields) call error_mesg('query_combined', 'invalid tracer index', FATAL)

query_combined = .false.
if (tracers(TRACER_ARRAY(model,index))%is_combined) query_combined = .true.

return
end function query_combined
!</FUNCTION>

!<SUBROUTINE NAME="set_tracer_atts">
!   <OVERVIEW>
!     A subroutine to allow the user set the tracer longname and units from the 
!     tracer initialization routine.
!   </OVERVIEW>
!   <DESCRIPTION>
!     A function to allow the user set the tracer longname and units from the 
!     tracer initialization routine. It seems sensible that the user who is 
!     coding the tracer code will know what units they are working in and it 
!     is probably safer to set the value in the tracer code rather than in 
!     the field table.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call set_tracer_atts(model, name, longname, units)
!   </TEMPLATE>

!   <IN NAME="model" TYPE="integer">
!     A parameter representing the component model in use.
!   </IN>
!   <IN NAME="name" TYPE="character">
!     Tracer name.
!   </IN>
!   <OUT NAME="longname" TYPE="character, optional">
!     A string describing the longname of the tracer for output to NetCDF files
!   </OUT>
!   <OUT NAME="units" TYPE="character, optional">
!     A string describing the units of the tracer for output to NetCDF files
!   </OUT>
!   <OUT NAME="set_tracer_atts" TYPE="character, optional">
!     A flag to show that 
!   </OUT>
subroutine set_tracer_atts(model, name, longname, units)

integer, intent(in)                    :: model
character(len=*), intent(in)           :: name
character(len=*), intent(in), optional :: longname, units

integer :: n

n = get_tracer_index(model,name)

    tracers(TRACER_ARRAY(model,n))%tracer_units   = units
    tracers(TRACER_ARRAY(model,n))%tracer_longname = longname

end subroutine set_tracer_atts
!</SUBROUTINE>

end module tracer_manager_mod

#ifdef test_tracer_manager

program test

use field_manager_mod, only : MODEL_ATMOS, MODEL_OCEAN, field_manager_init
use tracer_manager_mod
use fms_mod, only : error_mesg, FATAL
use mpp_domains_mod, only : domain2d, mpp_define_domains, cyclic_global_domain, &
                            mpp_get_data_domain, mpp_get_compute_domain, mpp_update_domains, &
                            mpp_get_global_domain
use diag_manager_mod, only : register_diag_field, send_data, diag_axis_init, diag_manager_init, diag_manager_end
use time_manager_mod

implicit none

integer, parameter :: nx = 100, ny = 100, nz = 10, ntsteps = 10, delt = 1, two_delt = 2*delt
real(r8), parameter :: advect_speed = .23
integer :: num_tracer_prog_atmos, num_tracer_diag_atmos, num_tracer_fam_atmos, num_tracers_atmos
integer :: num_tracer_prog_ocean, num_tracer_diag_ocean, num_tracer_fam_ocean, num_tracers_ocean
integer :: i, axes(3)
integer, allocatable, dimension(:) :: atmos_prog_ind, atmos_diag_ind, atmos_fam_ind, atmos_prog_tracer_diagnostic_id
integer, allocatable, dimension(:) :: atmos_diag_tracer_diagnostic_id
integer, allocatable, dimension(:) :: ocean_prog_ind, ocean_diag_ind, ocean_fam_ind

real(r8), allocatable, dimension(:,:,:,:,:), target :: atmos_tracers, ocean_tracers, atmos_fam_tracers
real(r8), allocatable, dimension(:,:,:,:), target :: atmos_diag_tracers, ocean_diag_tracers
real(r8), allocatable, dimension(:,:,:), target :: atmos_tendency, ocean_tendency  
real(r8), allocatable, dimension(:,:,:,:) :: atmos_vel, ocean_vel

real(r8) ::  missing_value = -1.e10
real(r8), dimension(max(nx,ny,nz)) :: data
logical :: result, flag

type(time_type) :: model_time
type(domain2d) :: domain ! just using a single domain type for all tracers here since
                         ! allocating storage in a single array.  This is inefficient in the
                         ! case where different tracers use different order advection schemes
                         ! need to consider extending mpp_update_domains to only update a portion
                         ! of the halo OR using a derived type for storing tracer arrays with 
                         ! a domain2d component (this seems unworkable because diag_axis_init would
                         ! need different calls for different domains ...)
   
character(len=128) :: name, longname, units, scheme_name, control, family_name
integer :: halo, ndivs, isc, iec, jsc, jec, isd, ied, jsd, jed, tau, taum1, taup1, tmp, n, m, order, mm
integer :: co2_index, color_index, temp_index,k,nfields
real(r8) :: surf_value,multiplier

call set_calendar_type(JULIAN)
call diag_manager_init

taum1=1
tau=2
taup1=3

model_time = set_date(1980,1,1,0,0,0)

call register_tracers(MODEL_ATMOS, num_tracers_atmos, num_tracer_prog_atmos, num_tracer_diag_atmos, &
                     num_tracer_fam_atmos)

call register_tracers(MODEL_OCEAN, num_tracers_ocean, num_tracer_prog_ocean, num_tracer_diag_ocean, &
                     num_tracer_fam_ocean)


call set_tracer_atts(MODEL_ATMOS,"rh","relative_humidity","percent")

do i=1,max(nx,ny,nz)
   data(i) = float(i)
enddo

ndivs = 1


if (num_tracer_prog_atmos > 0) then
   allocate(atmos_prog_ind(num_tracer_prog_atmos))
   allocate(atmos_prog_tracer_diagnostic_id(num_tracer_prog_atmos))
   call get_tracer_indices(MODEL_ATMOS,prog_ind=atmos_prog_ind)
   co2_index = get_tracer_index(MODEL_ATMOS,'co2',atmos_prog_ind)
   temp_index = get_tracer_index(MODEL_ATMOS,'temp',atmos_prog_ind)
   if (co2_index /= -1 )) write(*,'(a)') 'co2 exists in tracer table '
   halo=1 ! default halo size for 2nd order advection
   do n=1,num_tracer_prog_atmos
      flag =  query_method ('advection_scheme_horiz',MODEL_ATMOS,atmos_prog_ind(n),name)
      if (trim(name) == '4th_order') halo = 2
   enddo
   call mpp_define_domains( (/1,nx,1,ny/), (/1,ndivs/), domain, xflags = CYCLIC_GLOBAL_DOMAIN, &
                            yflags = CYCLIC_GLOBAL_DOMAIN, xhalo = halo, yhalo = halo)
   axes(1) = diag_axis_init('x',data(1:nx),units='meters',cart_name='x',domain2=domain)
   axes(2) = diag_axis_init('y',data(1:ny),units='meters',cart_name='y',domain2=domain)
   axes(3) = diag_axis_init('z',data(1:nz),units='mb',cart_name='z')
   call mpp_get_data_domain(domain,isd,ied,jsd,jed)
   call mpp_get_compute_domain(domain,isc,iec,jsc,jec)
   allocate(atmos_tracers(isd:ied,jsd:jed,nz,num_tracer_prog_atmos,1:3))
   allocate(atmos_tendency(isd:ied,jsd:jed,nz))
   allocate(atmos_vel(isd:ied,jsd:jed,nz,2))
   atmos_vel=advect_speed
   do i=1, num_tracer_prog_atmos
      write(*,'(a,i3,a)') 'Assigning tracer index ',atmos_prog_ind(i),' field to target array' 
      call get_tracer_names(MODEL_ATMOS,atmos_prog_ind(i),name,longname,units)
      atmos_prog_tracer_diagnostic_id(i) =  register_diag_field('tracers', trim(name),axes(1:3),model_time,&
                                trim(longname), trim(units), missing_value)      
      call assign_tracer_field(MODEL_ATMOS, atmos_prog_ind(i), data_tlevels=atmos_tracers(:,:,:,i,:))
   enddo
endif

if (num_tracer_prog_ocean > 0) then
   allocate(ocean_prog_ind(num_tracer_prog_ocean))
!   allocate(ocean_prog_tracer_diagnostic_id(num_tracer_prog_ocean))
   call get_tracer_indices(MODEL_OCEAN,prog_ind=ocean_prog_ind)
endif

if (num_tracer_diag_atmos > 0) then
   allocate(atmos_diag_ind(num_tracer_diag_atmos))
   allocate(atmos_diag_tracer_diagnostic_id(num_tracer_diag_atmos))
   allocate(atmos_diag_tracers(isd:ied,jsd:jed,nz,num_tracer_diag_atmos))
   call get_tracer_indices(MODEL_ATMOS,diag_ind=atmos_diag_ind)
   color_index = get_tracer_index(MODEL_ATMOS,'color',atmos_diag_ind)
   do m=1,num_tracer_diag_atmos
      call get_tracer_names(MODEL_ATMOS,atmos_diag_ind(m),name,longname,units)
      atmos_diag_tracer_diagnostic_id(m) = register_diag_field('tracers',trim(name),axes(1:3),model_time,&
                         trim(longname),trim(units),missing_value)
   enddo
endif

if (num_tracer_fam_atmos > 0) then
   allocate(atmos_fam_ind(num_tracer_fam_atmos))
   allocate(atmos_fam_tracers(isd:ied,jsd:jed,nz,num_tracer_fam_atmos,1:3))
   call get_tracer_indices(MODEL_ATMOS,fam_ind=atmos_fam_ind)
   do i=1,num_tracer_fam_atmos
      call assign_tracer_field(MODEL_ATMOS, atmos_fam_ind(i), data_tlevels=atmos_fam_tracers(:,:,:,i,:))
   enddo
endif

! atmos initialization

do m=1,num_tracer_prog_atmos
   call init_tracer(atmos_tracers(isc:iec,jsc:jec,:,m,tau),float(m), domain)
if(query_method ( 'profile_type',MODEL_ATMOS,atmos_prog_ind(m),longname,units)) then
       call get_tracer_names(MODEL_ATMOS,atmos_prog_ind(m),name,longname,units)
       call set_tracer_profile(MODEL_ATMOS,atmos_prog_ind(m),surf_value,multiplier)
write(*,'(2E12.4,i)') surf_value,multiplier,m
       atmos_tracers(:,:,:,m,tau) = surf_value

!  select case (size(atmos_tracers,5))
!  case (0)
! The atmosphere has k=1 at the top and k=kmax at the surface
       do k = size(atmos_tracers,3)-1,1,-1
         atmos_tracers(:,:,k,m,tau) = atmos_tracers(:,:,k+1,m,tau)*multiplier
       enddo
!  case (1)
!! The ocean has k=1 at the surface and k=kmax at the bottom
!       do k = 2,size(atmos_tracers,3)
!         atmos_tracers(:,:,k,m,tau) = atmos_tracers(:,:,k-1,m,tau)*multiplier
!       enddo
!  case default
!  end select
endif
!   atmos_tracers(isc:iec,jsc:jec,:,m,tau) = 1.0
   call mpp_update_domains(atmos_tracers(:,:,:,m,tau), domain)
   atmos_tracers(:,:,:,m,taum1) = atmos_tracers(:,:,:,m,tau)
enddo


! atmos time loop

!atmos_tendency = 0.0

do n=1,ntsteps
! combine family members
   do m=1,num_tracer_fam_atmos
! update pointers for family tracers
      call get_tracer_names(MODEL_ATMOS, atmos_fam_ind(m),name)
      call add_members_to_family(MODEL_ATMOS,name,cur=tau,prev=taum1,next=taup1)
      order = 2 ! Default advection scheme
      flag = query_method ('advection_scheme_horiz',MODEL_ATMOS,atmos_fam_ind(m),scheme_name)
      if(flag) then ! There is an advection scheme in the table so use that instead.
      select case (trim(scheme_name))
         case ('2nd_order') 
            order = 2
         case ('4th_order')
            order = 4
         case default
            call error_mesg('tracer manager', 'invalid tracer advection scheme', FATAL)
      end select
      endif
      call lateral_advection(order, atmos_tendency(:,:,:),tracer_now = atmos_fam_tracers(:,:,:,m,tau),&
                              vel_now = atmos_vel(:,:,:,:))
      atmos_fam_tracers(isc:iec,jsc:jec,:,m,taup1) = atmos_fam_tracers(isc:iec,jsc:jec,:,m,taum1) - &
                                             atmos_tendency(isc:iec,jsc:jec,:)*two_delt
      call mpp_update_domains(atmos_fam_tracers(:,:,:,m,taup1), domain)      
      call split_family_into_members(MODEL_ATMOS,name,cur=tau,prev=taum1,next=taup1)
      do mm = 1, num_tracer_prog_atmos
         call get_family_name(MODEL_ATMOS,atmos_prog_ind(mm),family_name)
         if (trim(name) == trim(family_name) .and. atmos_prog_tracer_diagnostic_id(mm) > 0) &
              result =  send_data(atmos_prog_tracer_diagnostic_id(mm),&
                                 atmos_tracers(isc:iec,jsc:jec,:,mm,tau),model_time)
      enddo
   enddo
   do m=1,num_tracer_prog_atmos
      call get_family_name(MODEL_ATMOS,atmos_prog_ind(m),family_name)
      if (family_name == 'orphan') then ! this is an orphan tracer so advect without splitting
         order=2! Default advection scheme
         flag = query_method ('advection_scheme_horiz',MODEL_ATMOS,atmos_prog_ind(m),scheme_name,control)
         if(flag) then ! There is an advection scheme in the table so use that instead.
         select case (trim(scheme_name))
         case ('2nd_order') 
            order = 2
         case ('4th_order')
            order = 4
         case default
            call error_mesg('tracer_manager', 'invalid tracer advection scheme', FATAL)
         end select
         endif
         call lateral_advection(order, atmos_tendency(:,:,:),tracer_now = atmos_tracers(:,:,:,m,tau),&
              vel_now = atmos_vel(:,:,:,:))
         atmos_tracers(isc:iec,jsc:jec,:,m,taup1) = atmos_tracers(isc:iec,jsc:jec,:,m,taum1) - &
              atmos_tendency(isc:iec,jsc:jec,:)*two_delt
         call mpp_update_domains(atmos_tracers(:,:,:,m,taup1), domain)
         if (atmos_prog_tracer_diagnostic_id(m) > 0) result =  send_data(atmos_prog_tracer_diagnostic_id(m),&
                                       atmos_tracers(isc:iec,jsc:jec,:,m,tau),model_time)
      endif
   enddo
   do m=1,num_tracer_diag_atmos
      if (m.eq.color_index) call calculate_color(atmos_tracers(:,:,:,temp_index,tau),atmos_diag_tracers(:,:,:,m))
      if (atmos_diag_tracer_diagnostic_id(m) > 0) &
           result = send_data(atmos_diag_tracer_diagnostic_id(m),&
          atmos_diag_tracers(isc:iec,jsc:jec,:,m),model_time)
   enddo
   tmp=taup1
   taup1=taum1
   taum1 = tau
   tau = tmp
   if (n.eq. ntsteps) call diag_manager_end(model_time)
   model_time = model_time + set_time(delt,0)      
enddo




contains 

  subroutine lateral_advection(order, tendency, tracer_now, vel_now, tracer_prev, vel_prev, tracer_next, vel_next)

    integer, intent(in) :: order
    real(r8), intent(inout), dimension(:,:,:) :: tendency, tracer_now
    real(r8), intent(in), dimension(:,:,:,:) :: vel_now
    real(r8), intent(in), dimension(:,:,:), optional :: tracer_prev, tracer_next
    real(r8), intent(in), dimension(:,:,:,:), optional :: vel_prev, vel_next

    integer :: ni, nj, nk, i, j, k
    real(r8) :: a,b
    real(r8), allocatable, dimension(:) :: fe, fn, fw, fs

    a=7.0/12.0;b=-1.0/12.0
    ni=size(tendency,1);nj=size(tendency,2);nk=size(tendency,3)

    allocate(fe(nk), fw(nk), fn(nk), fs(nk))

    select case (order)
    case (2)
       do j=2,nj-1
          do i=2,ni-1
             fe = 0.5*(tracer_now(i,j,:)+tracer_now(i+1,j,:))*vel_now(i+1,j,:,1)
             fn = 0.5*(tracer_now(i,j,:)+tracer_now(i,j+1,:))*vel_now(i,j+1,:,2)
             fw = 0.5*(tracer_now(i,j,:)+tracer_now(i-1,j,:))*vel_now(i,j,:,1) 
             fs = 0.5*(tracer_now(i,j,:)+tracer_now(i,j-1,:))*vel_now(i,j,:,2)
             tendency(i,j,:) = fe - fw + fn - fs
          enddo
       enddo
    case (4)
       do j=3,nj-2
          do i=3,ni-2
             fe = (a*(tracer_now(i,j,:)+tracer_now(i+1,j,:))+b*(tracer_now(i+2,j,:)+tracer_now(i-1,j,:)))*vel_now(i+1,j,:,1)
             fn = (a*(tracer_now(i,j,:)+tracer_now(i,j+1,:))+b*(tracer_now(i,j+2,:)+tracer_now(i,j-1,:)))*vel_now(i,j+1,:,2)
             fw = (a*(tracer_now(i,j,:)+tracer_now(i-1,j,:))+b*(tracer_now(i+1,j,:)+tracer_now(i-2,j,:)))*vel_now(i,j,:,1)
             fs = (a*(tracer_now(i,j,:)+tracer_now(i,j-1,:))+b*(tracer_now(i,j+1,:)+tracer_now(i,j-2,:)))*vel_now(i,j,:,2)
             tendency(i,j,:) = fe - fw + fn - fs
          enddo
       enddo
    case default
       call error_mesg('Tracer manager','invalid advection scheme', FATAL)
    end select

    deallocate(fe, fw, fn, fs)

    return
  end subroutine lateral_advection

  subroutine init_tracer(tracer, param, domain)

    real(r8), intent(out), dimension(:,:,:) :: tracer
    real(r8), intent(in)                    :: param
    type(domain2d), intent(in)          :: domain

    integer :: ni, nj, nk, i, j, isg, ieg, jsg, jeg, isc, iec, jsc, jec
    real(r8) :: xscale, yscale
    
    call mpp_get_global_domain(domain, xbegin=isg, xend=ieg, ybegin=jsg, yend=jeg, xsize=ni, ysize=nj)
    call mpp_get_compute_domain(domain,isc, iec, jsc, jec)

    xscale = ni/10
    yscale = nj/10

    do j=1,size(tracer,2)
       do i=1,size(tracer,1)
          tracer(i,j,:) = param*exp(-1*(isc+i-1-ni/2)**2/xscale**2)*exp(-1*(jsc+j-1-nj/2)**2/yscale**2)
       enddo
    enddo

    return
  end subroutine init_tracer

  subroutine calculate_color(temp,color)

    real(r8), dimension(:,:,:), intent(in) :: temp
    real(r8), dimension(:,:,:), intent(out) :: color

    color(:,:,:) = 0.5*temp(:,:,:)

    return

  end subroutine calculate_color

end program test
#endif


