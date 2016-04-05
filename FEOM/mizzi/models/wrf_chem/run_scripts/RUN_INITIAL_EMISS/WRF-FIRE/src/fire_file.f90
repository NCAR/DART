
   module fire_file

   implicit none

!  private
   public :: write_fire_file
   public :: max_fire_size
   public :: dealloc_fire_emis_glb_atts

   real, parameter :: lt_fac(0:23) = (/ .43, .43, .43, .43, .43, .43, .43, .43, &
                                        .43, 3., 6., 10., 14., 17., 14., 12., 9., &
                                         6., 3., .43, .43, .43, .43, .43 /)
   real :: max_fire_size = 2.e6      ! 2 km^2

   include "netcdf.inc"

   contains

   subroutine write_fire_file( domain, date, n_wrf_spc, m1, m2, &
                               lon, lon_ndx, lat_ndx, n_fire_spc, fire_size, &
                               fire_emissions, dxsqi, proj, n_files, fire_directory, &
                               fire_filename, file, ngatts, attrs )
!---------------------------------------------------------------------
!  write the netcdf bio emission file
!---------------------------------------------------------------------

    use wrf_utils
    use utils, only : wrf2fire_type, diag_level
    use srf_types, only : plant_cover, ntypes
    use attr_types, only : glb_att

!---------------------------------------------------------------------
!  dummy arguments
!---------------------------------------------------------------------
      integer, intent(in) :: domain
      integer, intent(in) :: n_wrf_spc
      integer, intent(in) :: n_fire_spc
      integer, intent(in) :: n_files
      integer, intent(in) :: m1, m2
      integer, intent(in) :: ngatts
      integer, intent(in) :: file
      integer, intent(in) :: lon_ndx(m1:m2)
      integer, intent(in) :: lat_ndx(m1:m2)
      real, intent(in)    :: dxsqi                                ! m^2
      real, intent(in)    :: lon(m1:m2)                           ! degrees east
      real, intent(in)    :: fire_size(m1:m2)                     ! m^2
      real, intent(in)    :: fire_emissions(n_fire_spc,m1:m2)
      character(len=10), intent(in) :: date
      character(len=*),  intent(in) :: fire_directory
      character(len=*),  intent(in) :: fire_filename(n_files)
      type(proj_info), intent(in)   :: proj
      type(glb_att), pointer :: attrs(:)
!---------------------------------------------------------------------
!  local variables
!---------------------------------------------------------------------
      real, parameter :: gas_confac  = 1.e6/24     ! convert from mole/day -> mole/km^2/hr
      real, parameter :: aero_confac = 1.e5/8.64   ! convert from kg/day   -> ug/m^2/s

      integer :: i, j
      integer :: k, m, n
      integer :: n1, n2
      integer :: ncid
      integer :: varid
      integer :: lon_id
      integer :: lat_id
      integer :: time_id
      integer :: zdim_stag_id
      integer :: string_id
      integer :: hr_ndx
      integer :: nhour
      integer :: lt_ndx
      integer :: dims(4)
      integer :: start_ndx(4)
      integer :: length(4)
      real    :: loc_hr
      real    :: wrk_sum
      real    :: rhour
      real    :: areai
      real    :: wrk_lon, latl,latu
      real    :: wrk_lt_fac(0:23)
      real    :: wrk_emiss_vec(n_wrf_spc)
      real    :: wrk_emiss(proj%ide,proj%jde,n_wrf_spc)
      real    :: wrk_fire_size(proj%ide,proj%jde)
      real    :: wrk(proj%ide,proj%jde)
      character(len=ntypes) :: suffix(4) = (/ 'AGEF', 'AGTF', 'AGSV', 'AGGR' /)
      character(len=132) :: outpname
      character(len=132) :: message
      character(len=32)  :: spc_name
      character(len=19)  :: date_time
      character(len=3)   :: numa
      character(len=3)   :: hra

      wrk_sum = 1./sum( lt_fac(:) )
      wrk_lt_fac(:) = wrk_sum*lt_fac(:)*24.

      write(numa,'(i3)') 100+domain
hour_loop : &
      do nhour = 0,23
        write(hra,'(i3)') 100+nhour
        date_time = date // '_' // hra(2:3) // ':00:00'
        outpname = 'wrffirechemi_d' // numa(2:3) // '_' // date_time
!-----------------------------------------------------------------------
!  create netcdf bio emission file and enter define mode
!-----------------------------------------------------------------------
        message = 'write_fire_file: Failed to create ' // trim( outpname )
        call handle_ncerr( nf_create( trim( outpname ), nf_clobber, ncid ), message )

!-----------------------------------------------------------------------
!  define the dimensions
!-----------------------------------------------------------------------
        call handle_ncerr( nf_def_dim( ncid, 'west_east', proj%ide, lon_id ), &
                           'write_fire_file: Failed to define longitude dimension' )
        call handle_ncerr( nf_def_dim( ncid, 'south_north', proj%jde, lat_id ), &
                           'write_fire_file: Failed to define latitude dimension' )
        call handle_ncerr( nf_def_dim( ncid, 'DateStrLen', 19, string_id ), &
                           'write_fire_file: Failed to define DateStrLen dimension' )
        call handle_ncerr( nf_def_dim( ncid, 'emissions_zdim_stag', 1, zdim_stag_id ), &
                           'write_fire_file: Failed to create emissions_zdim_stag dimension' )
        call handle_ncerr( nf_def_dim( ncid, 'Time', 1, time_id ), &
                           'write_fire_file: Failed to create Time dimension' )

!-----------------------------------------------------------------------
!  define the variables, set attributes
!-----------------------------------------------------------------------
        dims(1:2) = (/ string_id, time_id /)
        call handle_ncerr( nf_def_var( ncid, 'Times', nf_char, 2, dims(1:2), varid ), &
                           'write_fire_file: Failed to define Times variable' )

        dims(:) = (/ lon_id, lat_id, zdim_stag_id, time_id /)
        do n = 1,n_wrf_spc
          spc_name = 'ebu_in_' // trim(wrf2fire_type(n)%wrf_name)
          message = 'write_fire_file: Failed to define ' // trim(spc_name) // ' variable'
          call handle_ncerr( nf_def_var( ncid, trim(spc_name), nf_float, 4, dims, varid ), trim(message) )
          call write_fire_file_var_attributes( ncid, varid, spc_name, wrf2fire_type(n)%is_aerosol )
        end do
        dims(:3) = (/ lon_id, lat_id, time_id /)
        do n = 1,ntypes
          spc_name = 'MEAN_FCT_' // trim(suffix(n))
          message = 'write_fire_file: Failed to define ' // trim(spc_name) // ' variable'
          call handle_ncerr( nf_def_var( ncid, trim(spc_name), nf_float, 3, dims, varid ), trim(message) )
          call write_fire_file_var_attributes( ncid, varid, spc_name, .false., n )
        end do
        do n = 1,ntypes
          spc_name = 'FIRESIZE_' // trim(suffix(n))
          message = 'write_fire_file: Failed to define ' // trim(spc_name) // ' variable'
          call handle_ncerr( nf_def_var( ncid, trim(spc_name), nf_float, 3, dims, varid ), trim(message) )
          call write_fire_file_var_attributes( ncid, varid, spc_name, .false., n )
        end do
!---------------------------------------------------------------------
!  write the global attributes
!---------------------------------------------------------------------
        call write_fire_file_global_attributes( ncid, proj, n_files, fire_directory, fire_filename, &
                                                ngatts, attrs )
!-----------------------------------------------------------------------
!  leave define mode
!-----------------------------------------------------------------------
        call handle_ncerr( nf_enddef( ncid ), 'write_fire_file: Failed to leave define mode' )

!-----------------------------------------------------------------------
!  set time factor and calc fire size
!-----------------------------------------------------------------------
        wrk_emiss(:,:,:)   = 0.
        wrk_fire_size(:,:) = 0.
        rhour = real(nhour)
        do m = m1,m2
          n1 = lon_ndx(m)
          if( n1 /= 0 ) then
            n2 = lat_ndx(m)
            loc_hr = rhour + mod(lon(m)+360.,360.)/15.
            lt_ndx = mod( int(loc_hr),24 )
            wrk_fire_size(n1,n2) = wrk_fire_size(n1,n2) + fire_size(m)
            wrk_emiss_vec(:) = 0.
            if( proj%code /= PROJ_LATLON ) then
              areai = dxsqi
            else
              call ijll( real(n1), real(n2)-.5, proj, latl, wrk_lon )
              call ijll( real(n1), real(n2)+.5, proj, latu, wrk_lon )
              areai = 1./(proj%dx*proj%dy*(sin(latu*rad_per_deg) - sin(latl*rad_per_deg)))
            endif
wrf_spc_loop : &
            do n = 1,n_wrf_spc
              wrk_sum = 0.
              do k = 1,wrf2fire_type(n)%fire_cnt
                wrk_sum = wrk_sum + fire_emissions(wrf2fire_type(n)%fire_ndx(k,file),m) &
                                    * wrf2fire_type(n)%fire_wght(k)
              end do
              wrk_sum = wrk_sum * wrk_lt_fac(lt_ndx) * areai
              if( .not. wrf2fire_type(n)%is_aerosol ) then
                wrk_sum = wrk_sum*gas_confac
              else
                wrk_sum = wrk_sum*aero_confac
              endif
              wrk_emiss_vec(n) = wrk_emiss_vec(n) + wrk_sum
            end do wrf_spc_loop
            wrk_emiss(n1,n2,:) = wrk_emiss(n1,n2,:) + wrk_emiss_vec(:)
          endif
        end do
!-----------------------------------------------------------------------
!  limit fire size
!-----------------------------------------------------------------------
        wrk_fire_size(:,:) = min( max_fire_size,wrk_fire_size(:,:) )

        if( diag_level >= 300 ) then
          write(*,*) ' '
          write(*,'(''write_fire_file: max fire size = '',1pg15.7,'' @ (i,j) = '',2i3)') maxval(wrk_fire_size(:,:)),maxloc(wrk_fire_size(:,:))
        endif

!-----------------------------------------------------------------------
!  write the fire emission variables
!-----------------------------------------------------------------------
        start_ndx(1:2) = (/ 1,1 /)
        length(1:2)    = (/ 19, 1 /)
        call handle_ncerr( nf_inq_varid( ncid, 'Times', varid ), &
                           'write_fire_file: Failed to get Times variable id' )
        call handle_ncerr( nf_put_vara_text( ncid, varid, start_ndx(:2), length(:2), date_time ), &
                           'write_fire_file: Failed to write Times variable' )
        start_ndx(:) = 1
        length(:)    = (/ proj%ide, proj%jde, 1, 1 /)
        do n = 1,n_wrf_spc
          spc_name = 'ebu_in_' // trim(wrf2fire_type(n)%wrf_name)
          message = 'write_fire_file: Failed to get ' // trim(spc_name) // ' variable id'
          call handle_ncerr( nf_inq_varid( ncid, trim(spc_name), varid ), trim(message) )
          message = 'write_fire_file: Failed to write ' // trim(spc_name) // ' variable'
          call handle_ncerr( nf_put_vara_real( ncid, varid, start_ndx, length, wrk_emiss(:,:,n) ), trim(message) )
        end do
!-----------------------------------------------------------------------
!  write the srf cover fraction variables
!-----------------------------------------------------------------------
        start_ndx(:3) = 1
        length(:3)    = (/ proj%ide, proj%jde, 1 /)
        do n = 1,ntypes
          spc_name = 'MEAN_FCT_' // trim( suffix(n) )
          message = 'write_fire_file: Failed to get ' // trim(spc_name) // ' variable id'
          call handle_ncerr( nf_inq_varid( ncid, trim(spc_name), varid ), trim(message) )
          message = 'write_fire_file: Failed to write ' // trim(spc_name) // ' variable'
          call handle_ncerr( nf_put_vara_real( ncid, varid, start_ndx(1:3), length(1:3), plant_cover(domain)%type_frac(:,:,n) ), trim(message) )
        end do
!-----------------------------------------------------------------------
!  write the fire size variables per cover type
!-----------------------------------------------------------------------
        do n = 1,ntypes
          spc_name = 'FIRESIZE_' // trim( suffix(n) )
          message = 'write_fire_file: Failed to get ' // trim(spc_name) // ' variable id'
          call handle_ncerr( nf_inq_varid( ncid, trim(spc_name), varid ), trim(message) )
          message = 'write_fire_file: Failed to write ' // trim(spc_name) // ' variable'
          wrk(:,:) = plant_cover(domain)%type_frac(:,:,n) * wrk_fire_size(:,:)
          call handle_ncerr( nf_put_vara_real( ncid, varid, start_ndx(1:3), length(1:3), wrk ), trim(message) )
        end do
!---------------------------------------------------------------------
!  close wrf file
!---------------------------------------------------------------------
        message = 'Failed to close ' // trim(outpname)
        call handle_ncerr( nf_close( ncid ), message )       
      end do hour_loop

   end subroutine write_fire_file

   subroutine write_fire_file_var_attributes( ncid, varid, varname, is_aerosol, ndx )
!---------------------------------------------------------------------
!   write common variable attributes
!---------------------------------------------------------------------
    use srf_types, only : ntypes

!---------------------------------------------------------------------
!   dummy arguments
!---------------------------------------------------------------------
      integer, intent(in) :: ncid
      integer, intent(in) :: varid
      integer, optional, intent(in) :: ndx
      character(len=*), intent(in) :: varname
      logical, intent(in) :: is_aerosol
!---------------------------------------------------------------------
!   local variables
!---------------------------------------------------------------------
      integer :: m
      logical :: is_3d
      character(len=21)  :: desc(ntypes) = (/ 'extra tropical forest', 'tropical forest      ', &
                                              'savanna              ', 'grassland            ' /)
      character(len=132) :: message
      character(len=64)  :: description

      is_3d = .not. present( ndx )
      message = 'write_fire_file_var_attributes: Failed to create ' // trim(varname) // ' attribute'
      if( is_3d ) then
        call handle_ncerr( nf_put_att_text( ncid, varid, 'MemoryOrder', 3, 'XYZ' ), message )
        call handle_ncerr( nf_put_att_text( ncid, varid, 'description', 9, 'EMISSIONS' ), message )
        if( .not. is_aerosol ) then
          call handle_ncerr( nf_put_att_text( ncid, varid, 'units', 14, 'mole km-2 hr-1'), message )
        else
          call handle_ncerr( nf_put_att_text( ncid, varid, 'units', 10, 'ug m-2 s-1'), message )
        endif
        call handle_ncerr( nf_put_att_text( ncid, varid, 'stagger', 1, 'Z' ), message )
      else
        call handle_ncerr( nf_put_att_text( ncid, varid, 'MemoryOrder', 2, 'XY' ), message )
        if( varname(1:4) == 'MEAN' ) then
          description = 'mean fraction of ' // trim( desc(ndx) )
          call handle_ncerr( nf_put_att_text( ncid, varid, 'description', len_trim(description), trim(description) ), message )
          call handle_ncerr( nf_put_att_text( ncid, varid, 'units', 8, 'fraction' ), message )
        else
          description = 'mean firesize of ' // trim( desc(ndx) )
          call handle_ncerr( nf_put_att_text( ncid, varid, 'description', len_trim(description), trim(description) ), message )
          call handle_ncerr( nf_put_att_text( ncid, varid, 'units', 2, 'm2' ), message )
        endif
        call handle_ncerr( nf_put_att_text( ncid, varid, 'stagger', 1, ' ' ), message )
      endif
      m = 104
      call handle_ncerr( nf_put_att_int( ncid, varid, 'FieldType', nf_int, 1, m ), message )

   end subroutine write_fire_file_var_attributes

   subroutine write_fire_file_global_attributes( ncid, proj, n_files, fire_directory, fire_filename, &
                                                 ngatts, attrs )
!---------------------------------------------------------------------
!   write file global attributes
!---------------------------------------------------------------------
    use wrf_utils, only : proj_info
    use attr_types, only : glb_att
!---------------------------------------------------------------------
!  dummy arguments
!---------------------------------------------------------------------
   integer, intent(in)          :: ncid
   integer, intent(in)          :: ngatts
   integer, intent(in)          :: n_files
   character(len=*), intent(in) :: fire_directory
   character(len=*), intent(in) :: fire_filename(n_files)
   type(proj_info), intent(in)  :: proj
   type(glb_att), pointer :: attrs(:)
!---------------------------------------------------------------------
!  local variables
!---------------------------------------------------------------------
   integer            :: n
   character(len=132) :: text, message
   character(len=10)  :: ctime
   character(len=8)   :: cdate
   character(len=10)  :: t_string(2)

   message = 'global_attributes: Failed to write title'
   text    = 'FINNv1 hourly fire emissions'
   call handle_ncerr( nf_put_att_text( ncid, nf_global, 'Title', len_trim(text), trim(text) ), message )
   message = 'global_attributes: Failed to write History'
   call date_and_time( cdate, ctime )
   t_string(1) = cdate(1:4) // '-' // cdate(5:6) // '-' // cdate(7:8)
   t_string(2) = ctime(1:2) // ':' // ctime(3:4)
   text    = 'Created on ' // trim(t_string(1)) // ' at ' // trim(t_string(2))
   call handle_ncerr( nf_put_att_text( ncid, nf_global, 'History', len_trim(text), trim(text) ), message )

   message = 'global_attributes: Failed to write Directory'
   text    = trim( fire_directory )
   call handle_ncerr( nf_put_att_text( ncid, nf_global, 'Directory', len_trim(text), trim(text) ), message )
   message = 'global_attributes: Failed to write Files'
   text    = trim(fire_filename(1))
   if( n_files > 1 ) then
     text = trim(fire_filename(1)) // ','
     do n = 2,n_files
       if( n /= n_files ) then
         text(len_trim(text)+1:) = trim(fire_filename(n)) // ','
       else
         text(len_trim(text)+1:) = trim(fire_filename(n))
       endif
     end do
   endif
   call handle_ncerr( nf_put_att_text( ncid, nf_global, 'Files', len_trim(text), trim(text) ), message )
   message = 'global_attributes: Author to write Files'
   text    = 'fire_emis'
   call handle_ncerr( nf_put_att_text( ncid, nf_global, 'Author', len_trim(text), trim(text) ), message )

   if( ngatts > 0 ) then
     call set_fire_emis_glb_atts( ncid, ngatts, attrs )
   endif

   end subroutine write_fire_file_global_attributes

   subroutine set_fire_emis_glb_atts( ncid, ngatts, attrs )
!---------------------------------------------------------------------
!   set the global attributes
!---------------------------------------------------------------------
   use attr_types, only : glb_att

   implicit none

!---------------------------------------------------------------------
!   dummy arguments
!---------------------------------------------------------------------
   integer, intent(in)    :: ncid
   integer, intent(in)    :: ngatts
   type(glb_att), pointer :: attrs(:)

!---------------------------------------------------------------------
!   local variables
!---------------------------------------------------------------------
   integer              :: m
   integer              :: attr_len
   integer              :: attr_xtype
   integer              :: slen
   character(len=132)   :: message
   character(len=132)   :: attr_name

!---------------------------------------------------------------------
!   loop over glb attributes
!---------------------------------------------------------------------
glb_attr_loop : &
   do m = 1,ngatts
     attr_name = trim(attrs(m)%name)
     slen      = len_trim(attr_name)
     write(message,*) 'set_glb_att: Failed to define glb att ',trim(attr_name)
     attr_len   = attrs(m)%len
     attr_xtype = attrs(m)%type
     select case( attrs(m)%type )
       case( nf_byte )
         call handle_ncerr( nf_put_att_int1( ncid, nf_global, attr_name(:slen), attr_xtype, attr_len, attrs(m)%attr_byte ), message )       
       case( nf_char )
         call handle_ncerr( nf_put_att_text( ncid, nf_global, attr_name(:slen), attr_len, attrs(m)%attr_char ), message )       
       case( nf_short )
         call handle_ncerr( nf_put_att_int2( ncid, nf_global, attr_name(:slen), attr_xtype, attr_len, attrs(m)%attr_short ), message )       
       case( nf_int )
         call handle_ncerr( nf_put_att_int( ncid, nf_global, attr_name(:slen), attr_xtype, attr_len, attrs(m)%attr_int ), message )       
       case( nf_float )
         call handle_ncerr( nf_put_att_real( ncid, nf_global, attr_name(:slen), attr_xtype, attr_len, attrs(m)%attr_real ), message )       
       case( nf_double )
         call handle_ncerr( nf_put_att_double( ncid, nf_global, attr_name(:slen), attr_xtype, attr_len, attrs(m)%attr_dbl ), message )       
     end select
   end do glb_attr_loop

   end subroutine set_fire_emis_glb_atts

   subroutine dealloc_fire_emis_glb_atts( ngatts, attrs )
!---------------------------------------------------------------------
!   deallocate variables
!---------------------------------------------------------------------
   use attr_types, only : glb_att

!---------------------------------------------------------------------
!   dummy arguments
!---------------------------------------------------------------------
   integer, intent(in)    :: ngatts
   type(glb_att), pointer :: attrs(:)

   integer :: m

   do m = 1,ngatts
     select case( attrs(m)%type )
       case( nf_byte )
         deallocate( attrs(m)%attr_byte )
       case( nf_short )
         deallocate( attrs(m)%attr_short )
       case( nf_int )
         deallocate( attrs(m)%attr_int )
       case( nf_float )
         deallocate( attrs(m)%attr_real )
       case( nf_double )
         deallocate( attrs(m)%attr_dbl )
     end select
   end do

   end subroutine dealloc_fire_emis_glb_atts

   subroutine handle_ncerr( ret, mes )
!---------------------------------------------------------------------
!	... netcdf error handling routine
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!	... dummy arguments
!---------------------------------------------------------------------
   integer, intent(in) :: ret
   character(len=*), intent(in) :: mes

   if( ret /= nf_noerr ) then
      write(*,*) nf_strerror( ret )
      stop 'netcdf error'
   endif

   end subroutine handle_ncerr

   end module fire_file
