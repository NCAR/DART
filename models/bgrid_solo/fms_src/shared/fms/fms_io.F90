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
module fms_io_mod

  use mpp_io_mod, only : mpp_open, mpp_close, mpp_io_init, mpp_io_exit, &
                         MPP_NETCDF, MPP_ASCII, MPP_MULTI, MPP_SINGLE, &
                         mpp_read, mpp_write, mpp_write_meta, &
                         mpp_get_info, mpp_get_atts, MPP_IEEE32, &
                         MPP_OVERWR, fieldtype, axistype, atttype, &
                         MPP_RDONLY, MPP_NATIVE, MPP_DELETE, &
                         default_field, default_axis, default_att, &
                         mpp_get_fields, MPP_SEQUENTIAL, MPP_DIRECT
  use mpp_domains_mod, only : domain2d, domain1d, mpp_get_domain_components, &
                              mpp_get_compute_domain, mpp_get_data_domain, &
                              mpp_get_global_domain, NULL_DOMAIN1D, &
                              NULL_DOMAIN2D, mpp_global_field, operator( == )
                             
  use mpp_mod, only : mpp_error, FATAL, NOTE, mpp_pe, mpp_root_pe, mpp_npes, &
                      stdlog, stdout, mpp_broadcast, ALL_PES, &
                      mpp_chksum
  implicit none
  private

  integer, parameter, private :: max_files=50
  integer, parameter, private :: max_fields=50
  integer, parameter, private :: max_axes=20
  integer, parameter, private :: max_atts=20

  type buff_type
     real, dimension(:,:,:), pointer :: buffer
  end type buff_type

  type file_type
     integer :: unit ! mpp_io unit for netcdf file
     integer :: ndim, nvar, natt, max_ntime
     logical :: domain_present(max_fields)
     character(len=128) :: filename
     integer :: siz(max_fields,4)  ! X/Y/Z/T extent of fields (data domain size for 
     ! distributed writes;global size for reads)
     integer :: gsiz(max_fields,4) ! global X/Y/Z/T extent of fields
     integer :: unit_tmpfile(max_fields)
     character(len=128) :: fieldname(max_fields)
     type(buff_type), dimension(:), pointer :: field_buffer 
     type(fieldtype), dimension(max_fields) :: fields
     type(axistype),  dimension(max_axes) :: axes   ! spatial axes
     type(atttype),  dimension(max_atts) :: atts
     type(domain2D) :: domain(max_fields)
  end type file_type

  interface read_data
     module procedure read_data_3d_new
     module procedure read_data_2d_new
     module procedure read_data_1d_new
     module procedure read_data_scalar_new
     module procedure read_data_2d, read_ldata_2d, read_idata_2d
     module procedure read_data_3d, read_data_4d
     module procedure read_cdata_2d,read_cdata_3d,read_cdata_4d
  end interface

  interface write_data
     module procedure write_data_3d_new
     module procedure write_data_2d_new
     module procedure write_data_1d_new
     !   module procedure write_data_scalar_new
     module procedure write_data_2d, write_ldata_2d, write_idata_2d
     module procedure write_data_3d, write_data_4d
     module procedure write_cdata_2d,write_cdata_3d,write_cdata_4d
  end interface


  integer, private :: num_files_r=0 ! number of currently opened files for reading
  integer, private :: num_files_w=0 ! number of currently opened files for writing

  character(len=32), private :: threading_read, fileset_read, threading_write,&
       fileset_write, format ! global i/o settings

  integer, private :: thread_r, thread_w, fset_r, fset_w, form

  logical, private :: module_is_initialized = .FALSE.
  logical, private :: read_all_pe = .TRUE.
  character(len=64) :: iospec_ieee32 = '-N ieee_32'
  
  !------ private data, pointer to current 2d domain ------
  ! entrained from fms_mod.  This will be deprecated in the future.

  type(domain2D), pointer, private :: Current_domain

  integer, private :: is,ie,js,je      ! compute domain
  integer, private :: isd,ied,jsd,jed  ! data domain
  integer, private :: isg,ieg,jsg,jeg  ! global domain


  type(file_type), dimension(max_files), private :: files_read
  type(file_type), dimension(max_files), private :: files_write
  type(file_type) :: default_file

  public :: read_data, write_data, fms_io_init, fms_io_exit, field_size
  public :: open_namelist_file, open_restart_file, open_ieee32_file, close_file 
  public :: set_domain, nullify_domain, get_domain_decomp

  private :: lookup_field_w, lookup_axis, unique_axes
!  ---- version number -----

  character(len=128) :: version = '$Id$'
  character(len=128) :: tagname = '$Name$'

contains

  subroutine fms_io_init()
    ! assign default values for default_file
    IMPLICIT NONE

    integer  :: i,j, unit, io_status
    logical :: file_exist

    namelist /fms_io_nml/ threading_read, fileset_read, threading_write,&
         fileset_write, format, read_all_pe, iospec_ieee32

    call mpp_io_init()
    if (module_is_initialized) return

    threading_read='multi';fileset_read='single';format='netcdf'
    threading_write='multi';fileset_write='multi'

    call mpp_open(unit, 'input.nml',form=MPP_ASCII,action=MPP_RDONLY)
    read(unit,fms_io_nml,iostat=io_status)
    write(stdlog(), fms_io_nml)
    if (io_status > 0) then
       call mpp_error(FATAL,'=>Error reading input.nml')
    endif
    call mpp_close (unit)

    ! take namelist options if present

    select case (fileset_read) 
    case ('multi')
       fset_r = MPP_MULTI
    case ('single')
       fset_r = MPP_SINGLE
    case default
       call mpp_error(FATAL,'invalid fileset option')
    end select

    select case (threading_read) 
    case ('multi')
       thread_r = MPP_MULTI
    case ('single')
       thread_r = MPP_SINGLE
    case default
       call mpp_error(FATAL,'invalid threading option')
    end select

    ! take namelist options if present

    select case (fileset_write) 
    case ('multi')
       fset_w = MPP_MULTI
    case ('single')
       fset_w = MPP_SINGLE
    case default
       call mpp_error(FATAL,'invalid fileset option')
    end select

    select case (threading_write) 
    case ('multi')
       thread_w = MPP_MULTI
    case ('single')
       thread_w = MPP_SINGLE
    case default
       call mpp_error(FATAL,'invalid threading option')
    end select

    select case(format)
    case ('netcdf')
       form=MPP_NETCDF
    case default
       call mpp_error(FATAL,'only NetCDF I/O currently supported with fms_io')
    end select
    default_file%unit = -1
    default_file%ndim = 0  
    default_file%nvar = 0
    default_file%max_ntime = -1
    default_file%natt = 0
    default_file%siz(:,:) = 0
    default_file%gsiz(:,:) = 0
    default_file%unit_tmpfile(:) = -1
    default_file%filename = 'none'
    default_file%fieldname = 'none' 
    default_file%domain_present(:)=.false.
    default_file%domain(:)=NULL_DOMAIN2D
    default_file%fields(:) = default_field
    default_file%axes(:) = default_axis
    default_file%atts(:) = default_att      

    do i=1, max_files
       files_write(i) = default_file
       files_read(i)  = default_file
    enddo

    !---- initialize module domain2d pointer ----

    nullify (Current_domain)

    module_is_initialized = .TRUE.
    write (stdlog(),'(/,80("="),/(a))') trim(version), trim(tagname)

  end subroutine fms_io_init

  subroutine fms_io_exit()
    ! <OVERVIEW>
    ! This routine is called after all fields have been written to temporary files
    ! The netcdf files are created here
    ! </OVERVIEW>
    IMPLICIT NONE
    integer, parameter :: max_axis_size=10000
    integer :: i,j,k,unit,unit2,index_field
    integer :: num_axes, t_axis_id, x_axis_id, y_axis_id, z_axis_id
    integer, dimension(4) :: size_field, global_size ! x/y/z/t
    integer :: siz_x_axes(max_axes), siz_y_axes(max_axes), siz_z_axes(max_axes), max_t_size
    integer :: x_axes(max_axes), y_axes(max_axes), z_axes(max_axes)
    integer :: num_x_axes, num_y_axes, num_z_axes
    type(domain1d) :: domain_x(max_fields), domain_y(max_fields), x_domains(max_axes), y_domains(max_axes)
    real, dimension(max_axis_size) :: axisdata
    real :: tlev

    character (len=128) :: axisname,filename, fieldname,temp_name

    type(domain2D) :: domain


    do i=1,max_axis_size
       axisdata(i) = i
    enddo

    ! each field has an associated domain type (may be undefined).
    ! each file only needs to write unique axes (i.e. if 2 fields share an identical axis, then only write the axis once)
    ! unique axes are defined by the global size and domain decomposition (i.e. can support identical axis sizes with
    ! different domain decomposition)

    do i = 1, num_files_w

       ! determine maximum number of time levels for this file
       files_write(i)%max_ntime = 0
       do j=1,files_write(i)%nvar
          files_write(i)%max_ntime = max(files_write(i)%max_ntime,files_write(i)%siz(j,4))
       enddo

       allocate(files_write(i)%field_buffer(files_write(i)%nvar))

       filename = files_write(i)%filename
       call mpp_open(unit,trim(filename),action=MPP_OVERWR,form=form,threading=thread_w,&
            fileset=fset_w)

       do j = 1, max_fields
          if (files_write(i)%domain_present(j)) then
             call mpp_get_domain_components(files_write(i)%domain(j), domain_x(j), domain_y(j))
          else
             domain_x(j) = NULL_DOMAIN1D
             domain_y(j) = NULL_DOMAIN1D
          endif
       enddo

       siz_x_axes = -1;siz_y_axes=-1;siz_z_axes=-1
       x_axes = -1; y_axes = -1; z_axes=-1
       x_domains(:) = NULL_DOMAIN1D; y_domains(:) = NULL_DOMAIN1D

       x_axes = unique_axes(files_write(i)%gsiz(:,1),domain_x(:))
       do j=1,max_axes
          if (x_axes(j) > 0) then
             siz_x_axes(j) = files_write(i)%gsiz(x_axes(j),1) ! global array sizes
             x_domains(j) = domain_x(x_axes(j))
          endif
       end do
       y_axes = unique_axes(files_write(i)%gsiz(:,2),domain_y(:))
       do j=1,max_axes
          if (y_axes(j) > 0) then
             siz_y_axes(j) = files_write(i)%gsiz(y_axes(j),2) ! global array sizes
             y_domains(j) = domain_y(y_axes(j))
          endif
       end do
       z_axes = unique_axes(files_write(i)%gsiz(:,3))
       do j=1,max_axes
          if (z_axes(j) > 0) then
             siz_z_axes(j) = files_write(i)%gsiz(z_axes(j),3) ! global array sizes
          endif
       end do

       num_axes=0
       num_x_axes=0
       j=1
       do while (x_axes(j) > 0) 
          if (j < 10) then
             write(axisname,'(a,i1)') 'xaxis_',j
          else
             write(axisname,'(a,i2)') 'xaxis_',j
          endif
          num_axes=num_axes+1
          if (files_write(i)%domain_present(x_axes(j))) then
             call mpp_write_meta(unit,files_write(i)%axes(num_axes),axisname,'none',axisname, &
                  data=axisdata(1:siz_x_axes(j)),domain=domain_x(x_axes(j)),cartesian='X')
          else
!             if((mpp_pe()==mpp_root_pe().and.thread_w==MPP_SINGLE).or.(thread_w==MPP_MULTI)) &
                  call mpp_write_meta(unit,files_write(i)%axes(num_axes),axisname,'none',axisname, &
                  data=axisdata(1:siz_x_axes(j)),cartesian='X')
          endif
          j=j+1
          num_x_axes=num_x_axes+1
          if (num_axes > max_axes) call mpp_error(FATAL,'# axes exceeded max_axes in fms_io')
       enddo

       num_y_axes=0
       j=1
       do while (y_axes(j) > 0) 
          if (j < 10) then
             write(axisname,'(a,i1)') 'yaxis_',j
          else
             write(axisname,'(a,i2)') 'yaxis_',j
          endif
          num_axes=num_axes+1
          if (files_write(i)%domain_present(y_axes(j))) then
             call mpp_write_meta(unit,files_write(i)%axes(num_axes),axisname,'none',axisname, &
                  data=axisdata(1:siz_y_axes(j)),domain=domain_y(y_axes(j)),cartesian='Y')
          else
!             if((mpp_pe()==mpp_root_pe().and.thread_w==MPP_SINGLE).or.(thread_w==MPP_MULTI)) &
                  call mpp_write_meta(unit,files_write(i)%axes(num_axes),axisname,'none',axisname, &
                  data=axisdata(1:siz_y_axes(j)),cartesian='Y')
          endif
          j=j+1
          num_y_axes=num_y_axes+1
          if (num_axes > max_axes) call mpp_error(FATAL,'# axes exceeded max_axes in fms_io')
       enddo

       num_z_axes=0
       j=1
       do while (z_axes(j) > 0) 
          if (j < 10) then
             write(axisname,'(a,i1)') 'zaxis_',j
          else
             write(axisname,'(a,i2)') 'zaxis_',j
          endif
          num_axes=num_axes+1
!          if((mpp_pe()==mpp_root_pe().and.thread_w==MPP_SINGLE).or.(thread_w==MPP_MULTI)) &
               call mpp_write_meta(unit,files_write(i)%axes(num_axes),axisname,'none',axisname, &
               data=axisdata(1:siz_z_axes(j)),cartesian='Z')
          j=j+1
          num_z_axes=num_z_axes+1
          if (num_axes > max_axes) call mpp_error(FATAL,'# axes exceeded max_axes in fms_io')
       enddo

       ! write time axis  (comment out if no time axis)
!       if((mpp_pe()==mpp_root_pe().and.thread_w==MPP_SINGLE).or.(thread_w==MPP_MULTI)) &
            call mpp_write_meta(unit,files_write(i)%axes(num_axes+1),&
            'Time','time level','Time',cartesian='T')
       t_axis_id = num_axes+1

       ! write metadata for fields
       do j = 1, files_write(i)%nvar
          size_field = files_write(i)%gsiz(j,:) 
          x_axis_id = lookup_axis(siz_x_axes,size_field(1),x_domains,domain_x(j))            
          y_axis_id = lookup_axis(siz_y_axes,size_field(2),y_domains,domain_y(j))+num_x_axes            
          z_axis_id = lookup_axis(siz_z_axes,size_field(3))+num_x_axes+num_y_axes            
!          x_axis_id = x_axes(j)
!          y_axis_id = y_axes(j)+num_x_axes
!          z_axis_id = z_axes(j)+num_x_axes+num_y_axes
!          if((mpp_pe()==mpp_root_pe().and.thread_w == MPP_SINGLE).or.(thread_w==MPP_MULTI)) &
               call mpp_write_meta(unit,files_write(i)%fields(j), &
               (/files_write(i)%axes(x_axis_id),&
               files_write(i)%axes(y_axis_id),&
               files_write(i)%axes(z_axis_id),&
               files_write(i)%axes(t_axis_id)/),files_write(i)%fieldname(j),&
               'none',files_write(i)%fieldname(j),pack=1)
       enddo

       ! write values for ndim of spatial axes
       do j = 1, num_axes
!          if((mpp_pe()==mpp_root_pe().and.thread_w == MPP_SINGLE).or.(thread_w==MPP_MULTI)) &
               call mpp_write(unit,files_write(i)%axes(j))
       enddo

       ! retrieve and write data of each field
       do k = 1, files_write(i)%max_ntime
          do j = 1, files_write(i)%nvar
             size_field = files_write(i)%siz(j,:) ! local size here
             global_size = files_write(i)%gsiz(j,:) ! global size here
             unit2 = files_write(i)%unit_tmpfile(j)
             if (k == 1) then
                temp_name = trim(files_write(i)%filename)//'_'//&
                     trim(files_write(i)%fieldname(j))//'_tmp'
                call mpp_close(unit2)
                call mpp_open(unit2,temp_name,form=MPP_NATIVE,nohdrs=.true.,threading=MPP_MULTI, &
                     fileset=MPP_MULTI, action=MPP_RDONLY) 
                if(thread_w.eq.MPP_SINGLE .and. mpp_pe() == mpp_root_pe()) then
                   allocate(files_write(i)%field_buffer(j)%buffer(global_size(1),global_size(2),&
                        global_size(3)))
                else
                   allocate(files_write(i)%field_buffer(j)%buffer(size_field(1),size_field(2),&
                        size_field(3)))
                endif
                files_write(i)%unit_tmpfile(j) = unit2
             endif
             if ( k <= size_field(4)) then
                if((mpp_pe() == mpp_root_pe().and.thread_w==MPP_SINGLE).or.(thread_w==MPP_MULTI)) &
                     read(unit2) files_write(i)%field_buffer(j)%buffer
             else
                files_write(i)%field_buffer(j)%buffer = 0.0
             endif

             tlev=k
             if(files_write(i)%domain_present(j)) then
!                if (num_x_axes > 1 .or. num_y_axes > 1) call mpp_error(FATAL,&
!                     'restart data need to be on same grid when domain flag present')
                domain=files_write(i)%domain(j)
                if (mpp_pe() == mpp_root_pe().and. thread_w==MPP_SINGLE) then
                   call mpp_write(unit,files_write(i)%fields(j),&
                        files_write(i)%field_buffer(j)%buffer,tlev)
                else if (thread_w == MPP_MULTI) then
                   call mpp_write(unit,files_write(i)%fields(j),domain,&
                        files_write(i)%field_buffer(j)%buffer,tlev)
                endif
             else
                if (thread_w == MPP_MULTI .or. ((mpp_pe() == mpp_root_pe()).and.thread_w == MPP_SINGLE)) then
                   call mpp_write(unit,files_write(i)%fields(j),&
                        files_write(i)%field_buffer(j)%buffer,tlev)
                endif
             endif
          enddo ! end j loop
       enddo ! end k loop
       call mpp_close(unit)
       do j = 1, files_write(i)%nvar
          call mpp_close(files_write(i)%unit_tmpfile(j), action = MPP_DELETE)
       enddo
       deallocate(files_write(i)%field_buffer)

    enddo ! end i loop
    call mpp_io_exit()

    module_is_initialized = .false.
    num_files_w = 0
    num_files_r = 0

  end subroutine fms_io_exit

  subroutine write_data_3d_new(filename, fieldname, data, domain)
    !<OVERVIEW>
    ! This subroutine performs writing "fieldname" to file "filename". All values of "fieldname" 
    ! will be written to a temporary file. The final NETCDF file will be created only at a later step
    ! when the user calls fms_io_exit. Therefore, make sure that fms_io_exit is called after all
    ! fields have been written by this subroutine.
    !</OVERVIEW>

    IMPLICIT NONE
    character(len=*), intent(in) :: filename, fieldname 
    real, dimension(:,:,:), intent(in) :: data
    type(domain2d), intent(in), optional :: domain
    real, dimension(:,:,:), pointer ::global_data

    character(len=128) :: name, temp_name, axisname, old_name,units,longname,field_name
    ! temp_name: name of the temporary file
    integer :: unit,siz_in(4), siz(4), i, j, k
    integer :: nfile  ! index of the currently open file in array files
    !      integer :: ndim  ! number of dimensions of the file
    integer :: nvar !index of current field (variable)
    integer :: var_dim ! number of dimension of one variable
    integer :: index_field ! position of the fieldname in the list of variables
    integer :: unit2 ! unit of temporary file
    integer :: iscomp, iecomp, jscomp, jecomp, cxsize, cysize,cxsize_max, cysize_max
    integer :: isdata, iedata, jsdata, jedata, dxsize, dysize,dxsize_max, dysize_max
    integer :: isglobal, ieglobal, jsglobal, jeglobal, gxsize, gysize,gxsize_max, &
         gysize_max
    logical :: data_is_global(max_fields),data_is_local,data_is_compute
    logical :: new_variable= .true.
    logical :: file_open = .false.

    type(fieldtype) :: newfield
    type(axistype), dimension(:), allocatable :: var_axes !axes of one variable

    ! Initialize files to default values
    if(.not.module_is_initialized) call mpp_error(FATAL,'need to initialize fms_io first')

    data_is_global=.true.;data_is_local=.false.;data_is_compute=.false.


    ! Check if filename has been open  or not

    if(num_files_w == 0) then
       file_open=.false.
    else
       file_open=.false.
       do i=1,num_files_w
          nfile = i
          if (trim(files_write(i)%filename) == trim(filename)) then
             file_open = .true. !file already open
             exit
          endif
       enddo
    endif

    ! defer opening netcdf file until fms_io_exit
    if (.not.file_open) then 
       ! record the file name in array files_write
       num_files_w=num_files_w + 1
       nfile = num_files_w           
       files_write(nfile)%filename = trim(filename)         


    endif

    ! check if the field is new or not      
    ! get position and dimension of the field
    index_field =  lookup_field_w(nfile,fieldname)
    new_variable = .false.
    if(index_field < 0) new_variable = .true.

    if(new_variable) then
       ! open temporary file for writing data only, each field is written to a separate file
       temp_name = trim(filename)//'_'//trim(fieldname)//'_tmp'
       call mpp_open(unit2,temp_name,form=MPP_NATIVE,nohdrs=.true.,threading=MPP_MULTI,&
            fileset=MPP_MULTI, action=MPP_OVERWR)          

       files_write(nfile)%nvar = files_write(nfile)%nvar +1
       index_field = files_write(nfile)%nvar
       files_write(nfile)%fieldname(index_field) = fieldname
       files_write(nfile)%siz(index_field,1) = size(data,1)
       files_write(nfile)%siz(index_field,2) = size(data,2)
       files_write(nfile)%siz(index_field,3) = size(data,3)
       files_write(nfile)%siz(index_field,4) = 1
       files_write(nfile)%unit_tmpfile(index_field) = unit2

       if(PRESENT(domain)) then
          files_write(nfile)%domain_present(index_field)=.true.
          files_write(nfile)%domain(index_field)=domain
          data_is_global(index_field) = .false.
       else if (ASSOCIATED(Current_domain)) then
          ! if domain flag not present, alternatively use domain
          ! defined by previous call to set_domain.  This is needed
          ! in the case of atmospheric physics modules which don't
          ! have access to the Domain2d information
          files_write(nfile)%domain_present(index_field)=.true.
          files_write(nfile)%domain(index_field)=Current_domain
          data_is_global(index_field) = .false.
       else
          files_write(nfile)%domain_present(index_field)=.false.
          data_is_global(index_field) = .true.
       endif

       if (files_write(nfile)%domain_present(index_field)) then
          call mpp_get_global_domain(files_write(nfile)%domain(index_field),xsize=gxsize,ysize=gysize)
          files_write(nfile)%gsiz(index_field,1) = gxsize
          files_write(nfile)%gsiz(index_field,2) = gysize
          files_write(nfile)%gsiz(index_field,3) = size(data,3)
       else
          files_write(nfile)%gsiz(index_field,1) = size(data,1)
          files_write(nfile)%gsiz(index_field,2) = size(data,2)
          files_write(nfile)%gsiz(index_field,3) = size(data,3)
       endif
    else     ! this field  already exists
       ! get previously stored unit2
       unit2=files_write(nfile)%unit_tmpfile(index_field)
       ! Increase time level of the fieldname
       files_write(nfile)%siz(index_field,4) =  files_write(nfile)%siz(index_field,4) +1
    endif


    if(thread_w.eq.MPP_SINGLE .and. files_write(nfile)%domain_present(index_field) ) then
       gxsize = files_write(nfile)%gsiz(index_field,1)
       gysize = files_write(nfile)%gsiz(index_field,2)
!       files_write(nfile)%siz(index_field,1) = gxsize
!       files_write(nfile)%siz(index_field,2) = gysize
       allocate (global_data(gxsize,gysize,size(data,3)))
       call mpp_global_field(files_write(nfile)%domain(index_field),data,global_data)
       if(mpp_pe() == mpp_root_pe()) write(unit2) global_data
       deallocate(global_data)
       ! write data to temporary storage without halos
    else 
       write(unit2) data ! write data to a temporary file
    endif

  end subroutine write_data_3d_new
  ! ...

  subroutine write_data_2d_new(filename, fieldname, data, domain)

    IMPLICIT NONE
    character(len=*), intent(in) :: filename, fieldname 
    real, dimension(:,:), intent(in) :: data
    real, dimension(size(data,1),size(data,2),1) :: data_3d
    type(domain2d), intent(in), optional :: domain

    ! Initialize files to default values
    if(.not.module_is_initialized) call mpp_error(FATAL,'need to initialize fms_io first')

    data_3d(:,:,1) = data(:,:)

    call write_data_3d_new(filename, fieldname, data_3d, domain)

  end subroutine write_data_2d_new

  ! ........................................................

  subroutine write_data_1d_new(filename, fieldname, data)

    IMPLICIT NONE
    character(len=*), intent(in) :: filename, fieldname 
    real, dimension(:), intent(in) :: data
    real, dimension(size(data),1,1) :: data_3d


    ! Initialize files to default values
    if(.not.module_is_initialized) call mpp_error(FATAL,'need to initialize fms_io first')

    data_3d(:,1,1) = data(:)

    call write_data_3d_new(filename, fieldname, data_3d)

  end subroutine write_data_1d_new
  ! ..............


  ! ...

  function lookup_field_w(nfile,fieldname)
    IMPLICIT NONE
    ! Given fieldname, this function returns the field position in the model's fields list

    integer, intent(in) :: nfile
    character(len=*), intent(in) :: fieldname
    integer :: lookup_field_w
    integer :: j
    character(len=128) :: name

    lookup_field_w=-1

    do j = 1, files_write(nfile)%nvar
       name = files_write(nfile)%fieldname(j)
       if (trim(name) == trim(fieldname)) then
          lookup_field_w = j
          exit 
       endif
    enddo
    return

  end function lookup_field_w

  function lookup_axis(axis_sizes,siz,domains,dom)
    ! Given axis size (global), this function returns the axis id  

    IMPLICIT NONE
    integer, intent(in) :: axis_sizes(:), siz
    type(domain1d), optional :: domains(:)
    type(domain1d), optional :: dom
    integer :: lookup_axis
    integer :: i,j
    character(len=128) :: name

    lookup_axis=-1
    do j=1,size(axis_sizes)
       if (siz == axis_sizes(j)) then
          if (PRESENT(domains)) then
             if (dom == domains(j)) then 
                lookup_axis = j
                exit
             endif
          else
             lookup_axis = j
             exit
          endif
       endif
    enddo
    if (lookup_axis == -1) call mpp_error(FATAL,'error in lookup_axis')               
  end function lookup_axis

  subroutine field_size(filename, fieldname, siz)
    ! Given filename and fieldname, this subroutine returns the size of field
    character(len=*), intent(in) :: filename, fieldname
    integer, intent(inout) :: siz(:)
    character(len=128) :: name
    integer :: i, nfile, unit, ndim, nvar, natt, ntime, siz_in(4)
    logical :: file_opened, found

    type(fieldtype) :: fields(max_fields)
    if (size(siz) < 4) call mpp_error(FATAL,'size array must be >=4 in get_field_size')

    ! Need to check if filename has been opened or not
    nfile = 0
    file_opened=.false.
    do i=1,num_files_r
       if (trim(files_read(i)%filename) == trim(filename))  then
          nfile = i
          file_opened = .true.
          exit ! file is already opened
       endif
    enddo

    !Need to open the file now     
    !Only works for single NetCDF files for now ...

    found= .false.
    siz=-1

    if (.not. file_opened) then
       call mpp_open(unit,trim(filename),form=MPP_NETCDF,action=MPP_RDONLY,threading=MPP_MULTI, &
            fileset=MPP_SINGLE)
       call mpp_get_info(unit,ndim,nvar,natt,ntime)
       call mpp_get_fields(unit,fields(1:nvar))
       do i=1, nvar
          call mpp_get_atts(fields(i),name=name,siz=siz_in)
          if (trim(name) == trim(fieldname)) then
             siz(1:4)=siz_in(1:4)
             found = .true.
             exit
          endif
       enddo
       call mpp_close(unit)
    else
       do i=1, files_read(nfile)%nvar
          if (trim(fieldname) == trim(files_read(nfile)%fieldname(i))) then
             found = .true.
             siz = files_read(nfile)%siz(i,:)
             exit
          endif
       enddo
       if (.not. found) then
          call mpp_get_info(files_read(nfile)%unit,ndim,nvar,natt,ntime)
          call mpp_get_fields(files_read(nfile)%unit,fields(1:nvar))
          do i=1, nvar
             call mpp_get_atts(fields(i),name=name,siz=siz_in)
             if (trim(name) == trim(fieldname)) then
                siz(1:4) = siz_in(1:4)
                found=.true.
                exit
             endif
          enddo
       endif
    endif
!    if (.not. found) call mpp_error(FATAL,'field not found in file')

    return
  end subroutine field_size


  subroutine read_data_3d_new(filename,fieldname,data,domain, timelevel)
    IMPLICIT NONE
    character(len=*), intent(in) :: filename, fieldname
    real, dimension(:,:,:), intent(out) :: data ! 3 dimensional data    
    type(domain2d), intent(in), optional :: domain
    integer, intent(in) , optional :: timelevel

    character(len=128) :: name
    integer :: unit, siz_in(4), siz(4), i, j, k
    integer :: nfile  ! index of the opened file in array files
    integer :: ndim, nvar, natt, ntime,var_dim, tlev=1
    integer :: index_field ! position of the fieldname in the list of variables
    integer :: iscomp, iecomp, jscomp, jecomp, cxsize, cysize,cxsize_max, cysize_max
    integer :: isdata, iedata, jsdata, jedata, dxsize, dysize,dxsize_max, dysize_max
    integer :: isglobal, ieglobal, jsglobal, jeglobal, gxsize, gysize,gxsize_max, &
         gysize_max
    logical :: data_is_global,data_is_local,data_is_compute
    logical ::  file_opened 


    ! Initialize files to default values
    if(.not.module_is_initialized) call mpp_error(FATAL,'need to initialize fms_io first')

    data_is_global=.false.;data_is_local=.false.;data_is_compute=.false.
    file_opened=.false.

    if (PRESENT(timelevel)) then
       tlev = timelevel
    else
       tlev = 1
    endif
    if (PRESENT(domain)) then
       call mpp_get_compute_domain(domain,iscomp,iecomp,jscomp,jecomp,cxsize, & 
            cxsize_max,cysize,cysize_max)
       call mpp_get_data_domain(domain,isdata,iedata,jsdata,jedata,dxsize,dxsize_max,&
            dysize,dysize_max)
       call mpp_get_global_domain(domain,isglobal,ieglobal,jsglobal,jeglobal,gxsize, &
            gxsize_max,gysize,gysize_max)
       if (gxsize == size(data,1) .and. gysize == size(data,2)) data_is_global = .true.
       if (dxsize == size(data,1) .and. dysize == size(data,2)) data_is_local = .true.
       if (cxsize == size(data,1) .and. cysize == size(data,2)) data_is_compute = .true.
    else  if (ASSOCIATED(Current_domain)) then
       gxsize=ieg-isg+1
       gysize=jeg-jsg+1
       dxsize=ied-isd+1
       dysize=jed-jsd+1
       cxsize=ie-is+1
       cysize=je-js+1
       if (gxsize == size(data,1) .and. gysize == size(data,2)) data_is_global = .true.
       if (dxsize == size(data,1) .and. dysize == size(data,2)) data_is_local = .true.
       if (cxsize == size(data,1) .and. cysize == size(data,2)) data_is_compute = .true.
    else 
       data_is_global =.true.
       gxsize = size(data,1)
       gysize = size(data,2)
       siz_in(3) = size(data,3)
    endif

    if (data_is_global .and. fset_r == MPP_MULTI) &
         call mpp_error(FATAL,'cant do global read on multi fileset ')

    ! Need to check if filename has been opened or not
    nfile = 0
    if(num_files_r ==0) then
       file_opened = .false.
    else
       file_opened=.false.
       do i=1,num_files_r
          if (files_read(i)%filename == trim(filename))  then
             nfile = i
             file_opened = .true.
             exit ! file is already opened
          endif
       enddo
    endif
    if (.not. file_opened) then !Need to open the file now     
       if (fset_r == MPP_MULTI .and. thread_r == MPP_SINGLE) then
          call mpp_error(FATAL,'single-threaded reads from multi fileset not allowed')
       endif
       call mpp_open(unit,trim(filename),form=form,action=MPP_RDONLY,threading=thread_r, &
            fileset=fset_r)
       ! Increase num_files_r and set file_type 
       num_files_r=num_files_r + 1
       files_read(num_files_r)%filename = trim(filename)
       nfile = num_files_r
       files_read(nfile)%unit = unit
    else
       unit = files_read(nfile)%unit
    endif

    ! Get info of this file and field
    

    if ((thread_r == MPP_MULTI).or.(mpp_pe()==mpp_root_pe())) then
       call mpp_get_info(unit,files_read(nfile)%ndim, &
            files_read(nfile)%nvar, files_read(nfile)%natt,files_read(nfile)%max_ntime)
       if(files_read(nfile)%max_ntime < 1)  files_read(nfile)%max_ntime = 1
       nvar = files_read(nfile)%nvar  
       if (nvar > max_fields) call mpp_error(FATAL,'FMS_IO: max_fields too small')
       call mpp_get_fields(unit,files_read(nfile)%fields(1:nvar))     
       siz_in = 1
       index_field = -1
       do i=1,files_read(nfile)%nvar
          call mpp_get_atts(files_read(nfile)%fields(i),name=name,ndim=var_dim,siz=siz_in)
          if(var_dim .lt.3) then
             do j=var_dim+1,3
                siz_in(j)=1
             enddo
          endif
          if (trim(name) == trim(fieldname)) then ! found the variable
             index_field = i
             files_read(nfile)%fieldname(i) = fieldname
             files_read(nfile)%siz(i,:)  = siz_in
             files_read(nfile)%gsiz(i,:) = siz_in
             if (fset_r == MPP_SINGLE) then
                if (siz_in(1) /= gxsize .or. siz_in(2) /= gysize .or. siz_in(3) /= &
                     size(data,3)) call mpp_error(FATAL,'MPP_single:external field size mismatch')
             else if (fset_r == MPP_MULTI) then
                if (siz_in(1) /= dxsize .or. siz_in(2) /= dysize & 
                     .or. siz_in(3) /= size(data,3)) &
                     call mpp_error(FATAL,'MPP_multi: external field size mismatch')
             endif
             exit  !jump out of i loop
          endif
       enddo
       if(index_field <1) call mpp_error(FATAL, 'variable not found in read_data')
       if ( tlev < 1 .or. files_read(nfile)%max_ntime < tlev)    &
            call mpp_error(FATAL,'invalid time indices')
       if (data_is_global) then
          call mpp_read(unit,files_read(nfile)%fields(index_field),data(:,:,:),tlev)
       else 
          if (PRESENT(domain)) then
             call mpp_read(unit,files_read(nfile)%fields(index_field),domain,data,tlev)
          else
             call mpp_read(unit, files_read(nfile)%fields(index_field),Current_domain,data,tlev)
          endif
       endif
    endif

    return
  end subroutine read_data_3d_new
  !.... 

  subroutine read_data_2d_new(filename,fieldname,data,domain,timelevel)
    IMPLICIT NONE
    character(len=*), intent(in) :: filename, fieldname
    real, dimension(:,:), intent(out) :: data     !2 dimensional data 
    real, dimension(size(data,1),size(data,2),1) :: data_3d
    type(domain2d), intent(in), optional :: domain
    integer, intent(in) , optional :: timelevel

    data_3d = 0.0

    call read_data_3d_new(filename,fieldname,data_3d,domain,timelevel) 

    data(:,:) = data_3d(:,:,1)

  end subroutine read_data_2d_new

  subroutine read_data_1d_new(filename,fieldname,data,domain,timelevel)
    IMPLICIT NONE
    character(len=*), intent(in) :: filename, fieldname
    real, dimension(:), intent(out) :: data     !1 dimensional data 
    real, dimension(size(data),1,1) :: data_3d
    type(domain2d), intent(in), optional :: domain
    integer, intent(in) , optional :: timelevel

    data_3d = 0.0

    call read_data_3d_new(filename,fieldname,data_3d,domain,timelevel)

    data(:) = data_3d(:,1,1)

  end subroutine read_data_1d_new

  subroutine read_data_scalar_new(filename,fieldname,data,domain,timelevel)
    IMPLICIT NONE
    ! this subroutine is for reading a single number
    character(len=*), intent(in) :: filename, fieldname
    real, intent(out) :: data     !zero dimension data 
    real, dimension(1,1,1) :: data_3d
    type(domain2d), intent(in), optional :: domain
    integer, intent(in) , optional :: timelevel

    data_3d = 0.0

    call read_data_3d_new(filename,fieldname,data_3d,domain,timelevel) 

    data = data_3d(1,1,1)

  end subroutine read_data_scalar_new

  function unique_axes(array,dom)
    IMPLICIT NONE
    integer, dimension(:), intent(in) :: array
    type(domain1d), dimension(:), intent(in), optional :: dom
    integer, dimension(max_axes) :: unique_axes
    integer :: i,n,j
    logical :: dup

    unique_axes=-1
    n=1

    do i=1,size(array)
       dup = .false.
       if (array(i) == 0) exit 
       do j=i-1,1,-1
          if (array(i) == array(j)) then
             if (PRESENT(dom)) then
                if (dom(i) == dom(j)) then
                   dup = .true.
                   exit
                endif
             else
                dup = .true.
                exit
             endif
          endif
       enddo
       if (.not. dup) then
!          unique_axes(n) = array(i)
          unique_axes(n) = i
          n=n+1
       endif
    enddo

  end function unique_axes

  !#######################################################################
  !#######################################################################
  !   --------- routines for reading distributed data ---------
  ! before calling these routines the domain decompostion must be set
  ! by calling "set_domain" with the appropriate domain2d data type
  !
  ! reading can be done either by all PEs (default) or by only the root PE
  ! this is controlled by namelist variable "read_all_pe".
  !#######################################################################

  subroutine read_data_2d ( unit, data, end )

    integer, intent(in)                        :: unit
    real,    intent(out), dimension(isd:,jsd:) :: data
    logical, intent(out), optional             :: end

    real, dimension(isg:ieg,jsg:jeg) :: gdata
    integer :: len

    !------------------------------
    include "read_data_2d.inc"
    !------------------------------

  end subroutine read_data_2d

  !#######################################################################

  subroutine read_ldata_2d ( unit, data, end )

    integer, intent(in)                        :: unit
    logical, intent(out), dimension(isd:,jsd:) :: data
    logical, intent(out), optional             :: end

    logical, dimension(isg:ieg,jsg:jeg) :: gdata
    integer :: len

    !------------------------------
    include "read_data_2d.inc"
    !------------------------------

  end subroutine read_ldata_2d

  !#######################################################################

  subroutine read_idata_2d ( unit, data, end )

    integer, intent(in)                        :: unit
    integer, intent(out), dimension(isd:,jsd:) :: data
    logical, intent(out), optional             :: end

    integer, dimension(isg:ieg,jsg:jeg) :: gdata
    integer :: len

    !------------------------------
    include "read_data_2d.inc"
    !------------------------------

  end subroutine read_idata_2d

  !#######################################################################

  subroutine read_cdata_2d ( unit, data, end )

    integer, intent(in)                        :: unit
    complex,    intent(out), dimension(isd:,jsd:) :: data
    logical, intent(out), optional             :: end

    complex, dimension(isg:ieg,jsg:jeg) :: gdata
    integer :: len

    !------------------------------
    include "read_data_2d.inc"
    !------------------------------

  end subroutine read_cdata_2d

  !#######################################################################

  subroutine read_data_3d ( unit, data, end )

    integer, intent(in)                          :: unit
    real,    intent(out), dimension(isd:,jsd:,:) :: data
    logical, intent(out), optional               :: end

    real, dimension(isg:ieg,jsg:jeg,size(data,3)) :: gdata
    integer :: len

    !------------------------------
    include "read_data_3d.inc"
    !------------------------------

  end subroutine read_data_3d

  !#######################################################################

  subroutine read_cdata_3d ( unit, data, end )

    integer, intent(in)                          :: unit
    complex, intent(out), dimension(isd:,jsd:,:) :: data
    logical, intent(out), optional               :: end

    complex, dimension(isg:ieg,jsg:jeg,size(data,3)) :: gdata
    integer :: len

    !------------------------------
    include "read_data_3d.inc"
    !------------------------------

  end subroutine read_cdata_3d

  !#######################################################################

  subroutine read_data_4d ( unit, data, end )

    integer, intent(in)                            :: unit
    real,    intent(out), dimension(isd:,jsd:,:,:) :: data
    logical, intent(out), optional                 :: end

    real, dimension(isg:ieg,jsg:jeg,size(data,3),size(data,4)) :: gdata
    integer :: len

    ! WARNING: memory usage with this routine could be costly
    !------------------------------
    include "read_data_4d.inc"
    !------------------------------

  end subroutine read_data_4d

  !#######################################################################

  subroutine read_cdata_4d ( unit, data, end )

    integer, intent(in)                            :: unit
    complex, intent(out), dimension(isd:,jsd:,:,:) :: data
    logical, intent(out), optional                 :: end

    complex, dimension(isg:ieg,jsg:jeg,size(data,3),size(data,4)) :: gdata
    integer :: len

    ! WARNING: memory usage with this routine could be costly
    !------------------------------
    include "read_data_4d.inc"
    !------------------------------

  end subroutine read_cdata_4d

  !#######################################################################
  !     -------- routines for writing distributed data --------
  ! before calling these routines the domain decompostion must be set
  ! by calling "set_domain" with the appropriate domain2d data type
  !#######################################################################

  subroutine write_data_2d ( unit, data )

    integer, intent(in)                       :: unit
    real,    intent(in), dimension(isd:,jsd:) :: data

    real, dimension(isg:ieg,jsg:jeg) :: gdata

    !--------------------------
    include "write_data.inc"
    !--------------------------

  end subroutine write_data_2d

  !#######################################################################

  subroutine write_ldata_2d ( unit, data )

    integer, intent(in)                       :: unit
    logical, intent(in), dimension(isd:,jsd:) :: data

    logical, dimension(isg:ieg,jsg:jeg) :: gdata

    !--------------------------
    include "write_data.inc"
    !--------------------------

  end subroutine write_ldata_2d

  !#######################################################################

  subroutine write_idata_2d ( unit, data )

    integer, intent(in)                       :: unit
    integer, intent(in), dimension(isd:,jsd:) :: data

    integer, dimension(isg:ieg,jsg:jeg) :: gdata

    !--------------------------
    include "write_data.inc"
    !--------------------------

  end subroutine write_idata_2d

  !#######################################################################

  subroutine write_cdata_2d ( unit, data )

    integer, intent(in)                       :: unit
    complex, intent(in), dimension(isd:,jsd:) :: data

    complex, dimension(isg:ieg,jsg:jeg) :: gdata

    !--------------------------
    include "write_data.inc"
    !--------------------------

  end subroutine write_cdata_2d

  !#######################################################################

  subroutine write_data_3d ( unit, data )

    integer, intent(in) :: unit
    real,    intent(in), dimension(isd:,jsd:,:) :: data

    real, dimension(isg:ieg,jsg:jeg,size(data,3)) :: gdata

    !--------------------------
    include "write_data.inc"
    !--------------------------

  end subroutine write_data_3d

  !#######################################################################

  subroutine write_cdata_3d ( unit, data )

    integer, intent(in) :: unit
    complex, intent(in), dimension(isd:,jsd:,:) :: data

    complex, dimension(isg:ieg,jsg:jeg,size(data,3)) :: gdata

    !--------------------------
    include "write_data.inc"
    !--------------------------

  end subroutine write_cdata_3d

  !#######################################################################

  subroutine write_data_4d ( unit, data )

    integer, intent(in) :: unit
    real,    intent(in), dimension(isd:,jsd:,:,:) :: data

    real, dimension(isg:ieg,jsg:jeg,size(data,3),size(data,4)) :: gdata
    integer :: n

    !--------------------------

    if (.not.associated(Current_domain))  &
         call mpp_error(FATAL,'set_domain not called')

    ! get the global data and write only on root pe
    ! do this one field at a time to save memory
    do n = 1, size(data,4)
       call mpp_global_field ( Current_domain, data(:,:,:,n), gdata(:,:,:,n) )
    enddo
    if ( mpp_pe() == mpp_root_pe() ) write (unit) gdata

    !--------------------------

  end subroutine write_data_4d

  !#######################################################################

  subroutine write_cdata_4d ( unit, data )

    integer, intent(in) :: unit
    complex,    intent(in), dimension(isd:,jsd:,:,:) :: data

    complex, dimension(isg:ieg,jsg:jeg,size(data,3),size(data,4)) :: gdata
    integer :: n

    !--------------------------

    if (.not.associated(Current_domain)) call mpp_error(FATAL, 'set_domain not called')

    ! get the global data and write only on root pe
    ! do this one field at a time to save memory
    do n = 1, size(data,4)
       call mpp_global_field ( Current_domain, data(:,:,:,n), gdata(:,:,:,n) )
    enddo
    if ( mpp_pe() == mpp_root_pe() ) write (unit) gdata

    !--------------------------

  end subroutine write_cdata_4d

  !#######################################################################
  ! private routines (read_eof,do_read)
  ! this routine is called when an EOF is found while
  ! reading a distributed data file using read_data

  subroutine read_eof (end_found)
    logical, intent(out), optional :: end_found

    if (present(end_found))then
       end_found = .true.
    else
       call mpp_error(FATAL,'unexpected EOF')
    endif

  end subroutine read_eof

!#######################################################################
! determines if current pe should read data
! checks namelist variable read_all_pe

  function do_read ( )
    logical :: do_read

    do_read = mpp_pe() == mpp_root_pe() .or. read_all_pe

  end function do_read


!#######################################################################
!#######################################################################
!
! routines for opening specific types of files:
!
!                       form        action 
! open_namelist_file  MPP_ASCII   MPP_RDONLY  
! open restart_file   MPP_NATIVE
! open_ieee32_file    MPP_IEEE32
!
! all have: access=MPP_SEQUENTIAL, threading=MPP_SINGLE, nohdrs=.true.
! use the close_file interface to close these files
!
! if other types of files need to be opened the mpp_open and
! mpp_close routines in the mpp_io_mod should be used
!
!#######################################################################
! opens single namelist file for reading only by all PEs
! the default file opened is called "input.nml"

  function open_namelist_file (file) result (unit)
    character(len=*), intent(in), optional :: file
    integer :: unit

    if (.not.module_is_initialized) call fms_io_init ( )

    if (present(file)) then
       call mpp_open ( unit, file, form=MPP_ASCII, action=MPP_RDONLY, &
            access=MPP_SEQUENTIAL, threading=MPP_SINGLE )
    else
       call mpp_open ( unit, 'input.nml', form=MPP_ASCII, action=MPP_RDONLY, &
            access=MPP_SEQUENTIAL, threading=MPP_SINGLE )
    endif

  end function open_namelist_file

  !#######################################################################
  ! opens single restart file for reading by all PEs or
  ! writing by root PE only
  ! the file has native format and no mpp header records

  function open_restart_file (file, action) result (unit)
    character(len=*), intent(in) :: file, action
    integer :: unit

    integer :: mpp_action

    if (.not.module_is_initialized) call fms_io_init ( )

    !   --- action (read,write) ---

    select case (lowercase(trim(action)))
    case ('read')
       mpp_action = MPP_RDONLY
    case ('write')
       mpp_action = MPP_OVERWR
    case default
       call mpp_error(FATAL,'invalid option for argument action')
    end select

    call mpp_open ( unit, file, form=MPP_NATIVE, action=mpp_action, &
         access=MPP_SEQUENTIAL, threading=MPP_SINGLE, nohdrs=.true. )

  end function open_restart_file

  !#######################################################################
  ! opens single 32-bit ieee file for reading by all PEs or 
  ! writing by root PE only (writing is not recommended)
  ! the file has no mpp header records

  function open_ieee32_file (file, action) result (unit)
    character(len=*), intent(in) :: file, action
    integer :: unit

    integer :: mpp_action

    if (.not.module_is_initialized) call fms_io_init ( )

    !   --- action (read,write) ---

    select case (lowercase(trim(action)))
    case ('read')
       mpp_action = MPP_RDONLY
    case ('write')
       mpp_action = MPP_OVERWR
    case default
       call mpp_error (FATAL,'invalid option for argument action')
    end select

    if (iospec_ieee32(1:1) == ' ') then
       call mpp_open ( unit, file, form=MPP_IEEE32, action=mpp_action, &
            access=MPP_SEQUENTIAL, threading=MPP_SINGLE,    &
            nohdrs=.true. )
    else
       call mpp_open ( unit, file, form=MPP_IEEE32, action=mpp_action, &
            access=MPP_SEQUENTIAL, threading=MPP_SINGLE,    &
            nohdrs=.true., iospec=iospec_ieee32 )
    endif

  end function open_ieee32_file

  !#######################################################################

  subroutine close_file (unit, status)
    integer,          intent(in)           :: unit
    character(len=*), intent(in), optional :: status

    if (.not.module_is_initialized) call fms_io_init ( )

    if (unit == stdlog()) return

    if (present(status)) then
       if (lowercase(trim(status)) == 'delete') then
          call mpp_close (unit, action=MPP_DELETE)
       else
          call mpp_error(FATAL,'invalid value for status')
       endif
    else
       call mpp_close (unit)
    endif

  end subroutine close_file

  !#######################################################################
  !#######################################################################
  ! set_domain is called to save the domain2d data type prior to
  ! calling the distributed data I/O routines, read_data and write_data

  subroutine set_domain (Domain2)

    type(domain2D), intent(in), target :: Domain2

    if (.NOT.module_is_initialized) call fms_io_init ( )

    !  --- set_domain must be called before a read_data or write_data ---

    if (associated(Current_domain)) nullify (Current_domain)
    Current_domain => Domain2

    !  --- module indexing to shorten read/write routines ---

    call mpp_get_compute_domain (Current_domain,is ,ie ,js ,je )
    call mpp_get_data_domain    (Current_domain,isd,ied,jsd,jed)
    call mpp_get_global_domain  (Current_domain,isg,ieg,jsg,jeg)

    !-----------------------------------------------------------------------

  end subroutine set_domain

  subroutine nullify_domain ()


    if (.NOT.module_is_initialized) call fms_io_init ( )

    !  --- set_domain must be called before a read_data or write_data ---

    if (associated(Current_domain)) nullify (Current_domain)

    is=0;ie=0;js=0;je=0
    isd=0;ied=0;jsd=0;jed=0
    isg=0;ieg=0;jsg=0;jeg=0

  end subroutine nullify_domain

  !#######################################################################
  ! this will be a private routine with the next release
  ! users should get the domain decomposition from the domain2d data type

  subroutine get_domain_decomp ( x, y )

    integer, intent(out), dimension(4) :: x, y

    if (mpp_pe() == mpp_root_pe())  call mpp_error(NOTE, &
         'subroutine get_domain_decomp will be &
         &removed with the next release')
    x = (/ isg, ieg, is, ie /)
    y = (/ jsg, jeg, js, je /)

  end subroutine get_domain_decomp

 function lowercase (cs) 
 character(len=*), intent(in) :: cs
 character(len=len(cs))       :: lowercase 
 character :: ca(len(cs)) 

 integer, parameter :: co=iachar('a')-iachar('A') ! case offset
    
    ca = transfer(cs,"x",len(cs)) 
    where (ca >= "A" .and. ca <= "Z") ca = achar(iachar(ca)+co) 
    lowercase = transfer(ca,cs) 
    
 end function lowercase 

end module fms_io_mod


