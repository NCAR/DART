program update_wrf_bc

! <next four lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
! $Id$
!
! program to update BC file from 3dvar or filter output.
! current version reads only wrf-netcdf file format.
!

use               types_mod, only : r8
use           utilities_mod, only : file_exist, open_file, close_file, &
                                    initialize_utilities, finalize_utilities, register_module, &
                                    logfileunit
use module_netcdf_interface, only : get_dims_cdf, get_gl_att_real_cdf, put_gl_att_real_cdf, &
                                    get_var_3d_real_cdf, get_var_2d_real_cdf, put_var_3d_real_cdf, &
                                    put_var_2d_real_cdf, get_times_cdf, put_time_cdf
use       module_couple_uvw
use         module_timediff, only : time_diff, find_time_index

implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

!-----------------------------------------------------------------------
! Model namelist parameters with default values.
!-----------------------------------------------------------------------

integer :: mp_physics = 3

namelist /physics/ mp_physics

!-----------------------------------------------------------------------

integer, parameter :: max_3d_variables = 20, &
                      max_2d_variables = 20, &
                      max_times        = 100

character(len=80) :: wrf_3dvar_output_file, &
                     wrf_bdy_file

character(len=20) :: var_pref, var_name

character(len=20) :: var3d(max_3d_variables), &
                     var2d(max_2d_variables)

character(len=10), dimension(4) :: bdyname, tenname

character(len=19), dimension(max_times) :: udtime, bdytime, thisbdytime, nextbdytime

integer           :: ids, ide, jds, jde, kds, kde
integer           :: num3d, num2d, ndims, nmoist
integer           :: i,j,k,l,m,n
integer           :: ntimes_bdy, ntimes_ud, itime

integer, dimension(4) :: dims

integer, external :: iargc

real(r8), allocatable, dimension(:,:) :: tend2d, scnd2d, frst2d

real(r8), allocatable, dimension(:,:,:) :: tend3d, scnd3d, frst3d, full3d

real(r8), allocatable, dimension(:,:,:) :: u, v, w

real(r8), allocatable, dimension(:,  :) :: mu, mub, msfu, msfv, msfm

integer :: east_end, north_end

logical, parameter :: debug = .false.

real(r8) :: bdyfrq_old, bdyfrq

integer :: io, iunit

!----------------------------------------------------------------------

call initialize_utilities
call register_module(source, revision, revdate)
write(logfileunit,*)'STARTING update_wrf_bc ...'

! Reading the namelist input
if(file_exist('namelist.input')) then

   iunit = open_file('namelist.input', action = 'read')
   read(iunit, nml = physics, iostat = io )
   call close_file(iunit)

endif

write(logfileunit , nml=physics)
write(     *      , nml=physics)

select case(mp_physics)
case (0) ;
   nmoist = 1
case (1) ;
   nmoist = 3
case (2) ;
   nmoist = 6
case (3) ;
   nmoist = 3
case (4) ;
   nmoist = 5
case (5) ;
   nmoist = 2
case (6) ;
   nmoist = 6
case (98) ;
   nmoist = 3
case (99) ;
   nmoist = 5
case default ;
   print *, 'Microphysics package unknown. mp_physics = ', mp_physics
end select

wrf_3dvar_output_file='wrfinput_d01'
wrf_bdy_file  ='wrfbdy_d01'

!--boundary variables
bdyname(1)='_BXS'
bdyname(2)='_BXE'
bdyname(3)='_BYS'
bdyname(4)='_BYE'

!--boundary tendancy variables
tenname(1)='_BTXS'
tenname(2)='_BTXE'
tenname(3)='_BTYS'
tenname(4)='_BTYE'

!--3D need update
num3d = 5 + nmoist
var3d(1)='U'
var3d(2)='V'
var3d(3)='W'
var3d(4)='PH'
var3d(5)='T'
var3d(6)='QVAPOR'
var3d(7)='QCLOUD'
var3d(8)='QRAIN'
var3d(9)='QICE'
var3d(10)='QSNOW'
var3d(11)='QGRAUPEL'

!--2D need update
num2d=5
var2d(1)='MUB'
var2d(2)='MU'
var2d(3)='MAPFAC_U'
var2d(4)='MAPFAC_V'
var2d(5)='MAPFAC_M'

!---------------------------------------------------------------------

east_end=0
north_end=0

!---------------------------------------------------------------------
!-- Current time in file
call get_times_cdf( wrf_3dvar_output_file, 'Times', udtime, ntimes_ud, max_times, debug )

if(debug) print*, 'udtime = ',udtime(1)

!-- list of boundary times
call get_times_cdf( wrf_bdy_file, 'Times', bdytime, ntimes_bdy, max_times, debug )

!-- Time index of current tendency - to grab from wrfbdy
call find_time_index(bdytime, udtime(1), ntimes_bdy, itime)

if(debug) print*, 'itime = ',itime
if(debug) print*, 'bdytime = ',bdytime(itime)

!---------------------------------------------------------------------
!--First, the boundary frequency.
call get_times_cdf( wrf_bdy_file, 'md___thisbdytimee_x_t_d_o_m_a_i_n_m_e_t_a_data_', &
     thisbdytime, ntimes_bdy, max_times, debug )
call get_times_cdf( wrf_bdy_file, 'md___nextbdytimee_x_t_d_o_m_a_i_n_m_e_t_a_data_', &
     nextbdytime, ntimes_bdy, max_times, debug )
call time_diff( thisbdytime(itime), nextbdytime(itime), bdyfrq_old ) 

if(debug) then
   write(unit=*, fmt='(a, f12.2)') &
        'BDYFRQ=', bdyfrq_old
endif
!-- time diff between this time and most recent (to shorten bdyfrq)
call time_diff( udtime(1), nextbdytime(itime), bdyfrq ) 

if(debug) print*, 'New bdyfrq = ',bdyfrq

!-- put in the new thisbdytime
thisbdytime(itime) = udtime(1)
call put_time_cdf( wrf_bdy_file, 'md___thisbdytimee_x_t_d_o_m_a_i_n_m_e_t_a_data_', &
     thisbdytime(itime), itime, debug )

!--For 2D variables
!--Get mu, mub, msfu, msfv, and msfm

do n=1,num2d
   call get_dims_cdf( wrf_3dvar_output_file, trim(var2d(n)), dims, ndims, debug )

   select case(trim(var2d(n)))
   case ('MU') ;
      allocate(mu(dims(1), dims(2)))

      call get_var_2d_real_cdf( wrf_3dvar_output_file, trim(var2d(n)), mu, &
           dims(1), dims(2), 1, debug )

      east_end=dims(1)+1
      north_end=dims(2)+1
   case ('MUB') ;
      allocate(mub(dims(1), dims(2)))

      call get_var_2d_real_cdf( wrf_3dvar_output_file, trim(var2d(n)), mub, &
           dims(1), dims(2), 1, debug )
   case ('MAPFAC_U') ;
      allocate(msfu(dims(1), dims(2)))

      call get_var_2d_real_cdf( wrf_3dvar_output_file, trim(var2d(n)), msfu, &
           dims(1), dims(2), 1, debug )

      if(debug) then
         do j=1,dims(2)
            write(unit=*, fmt='(2(a,i5), a, f12.8)') &
                 'msfu(', dims(1), ',', j, ')=', msfu(dims(1),j)
         enddo
      endif

   case ('MAPFAC_V') ;
      allocate(msfv(dims(1), dims(2)))

      call get_var_2d_real_cdf( wrf_3dvar_output_file, trim(var2d(n)), msfv, &
           dims(1), dims(2), 1, debug )

      if(debug) then
         do i=1,dims(1)
            write(unit=*, fmt='(2(a,i5), a, f12.8)') &
                 'msfv(', i, ',', dims(2), ')=', msfv(i,dims(2))
         enddo
      endif

   case ('MAPFAC_M') ;
      allocate(msfm(dims(1), dims(2)))

      call get_var_2d_real_cdf( wrf_3dvar_output_file, trim(var2d(n)), msfm, &
           dims(1), dims(2), 1, debug )

      if(debug) then
         do i=1,dims(1)
            write(unit=*, fmt='(2(a,i5), a, f12.8)') &
                 'msfm(', i, ',', dims(2), ')=', msfm(i,dims(2))
         enddo
      endif

   case default ;
      print *, 'It is impossible here. var2d(n)=', trim(var2d(n))
   end select
enddo

if(debug) then
   write(unit=*, fmt='(2(2x,a,e20.12)/2x,a,i2,2x,a,4i6)') &
        'Sample mu =', mu(dims(1)/2,dims(2)/2), &
        'Sample mub=', mub(dims(1)/2,dims(2)/2), &
        'ndims=', ndims, 'dims=', (dims(i), i=1,ndims)
endif

if(east_end < 1 .or. north_end < 1) then
   write(unit=*, fmt='(a)') 'Wrong data for Boundary.'
   stop
endif

!---------------------------------------------------------------------

   do m=1,4
      var_name='MU' // trim(bdyname(m))

      call get_dims_cdf( wrf_bdy_file, trim(var_name), dims, ndims, debug )

      allocate(frst2d(dims(1), dims(2)))
      allocate(scnd2d(dims(1), dims(2)))
      allocate(tend2d(dims(1), dims(2)))

!-----Calculate variable at second time level
!--------Get variable tendancy at first time level
      var_name='MU' // trim(tenname(m))
      call get_var_2d_real_cdf( wrf_bdy_file, trim(var_name), tend2d, &
           dims(1), dims(2), itime, debug )

!--------Get variable at first time level
      var_name='MU' // trim(bdyname(m))
      call get_var_2d_real_cdf( wrf_bdy_file, trim(var_name), frst2d, &
           dims(1), dims(2), itime, debug )

      scnd2d = tend2d*bdyfrq_old + frst2d

!      call get_var_2d_real_cdf( wrf_bdy_file, trim(var_name), scnd2d, &
!                                dims(1), dims(2), 2, debug )



!-----calculate variable at first time level
      select case(m)
         case (1) ;		! West boundary
            do l=1,dims(2)
            do j=1,dims(1)
               frst2d(j,l)=mu(l,j)
            enddo
            enddo
         case (2) ;		! East boundary
            do l=1,dims(2)
            do j=1,dims(1)
               frst2d(j,l)=mu(east_end-l,j)
            enddo
            enddo
         case (3) ;		! South boundary
            do l=1,dims(2)
            do i=1,dims(1)
               frst2d(i,l)=mu(i,l)
            enddo
            enddo
         case (4) ;		! North boundary
            do l=1,dims(2)
            do i=1,dims(1)
               frst2d(i,l)=mu(i,north_end-l)
            enddo
            enddo
         case default ;
            print *, 'It is impossible here. mu, m=', m
      end select

!-----calculate new tendancy 

      tend2d = (scnd2d - frst2d)/bdyfrq

      if(debug) then
         write(unit=*, fmt='(a,i2,2x,2a/3(a,e20.12,4x)/a,i2,2x,a,4i6)') &
              'No.', n, 'Variable: ', trim(var_name), &
              'Sampe frst2d=', frst2d(dims(1)/2,dims(2)/2), &
                    'scnd2d=', scnd2d(dims(1)/2,dims(2)/2), &
                    'tend2d=', tend2d(dims(1)/2,dims(2)/2), &
              'ndims=', ndims, 'dims=', (dims(i), i=1,ndims)
      endif

!-----output new variable at first time level
      var_name='MU' // trim(bdyname(m))
      call put_var_2d_real_cdf( wrf_bdy_file, trim(var_name), frst2d, &
                                dims(1), dims(2), itime, debug )
!-----output new tendancy 
      var_name='MU' // trim(tenname(m))
      call put_var_2d_real_cdf( wrf_bdy_file, trim(var_name), tend2d, &
                                dims(1), dims(2), itime, debug )

      deallocate(frst2d)
      deallocate(scnd2d)
      deallocate(tend2d)
   enddo

!---------------------------------------------------------------------
!--For 3D variables

!--Get U
   call get_dims_cdf( wrf_3dvar_output_file, 'U', dims, ndims, debug )

   allocate(u(dims(1), dims(2), dims(3)))

   ids=1
   ide=dims(1)-1
   jds=1
   jde=dims(2)
   kds=1
   kde=dims(3)

   call get_var_3d_real_cdf( wrf_3dvar_output_file, 'U', u, &
                             dims(1), dims(2), dims(3), 1, debug )

   if(debug) then
      do j=1,dims(2)
         write(unit=*, fmt='(2(a,i5), a, f12.8)') &
              'u(', dims(1), ',', j, ',1)=', u(dims(1),j,1)
      enddo
   endif

!--Get V
   call get_dims_cdf( wrf_3dvar_output_file, 'V', dims, ndims, debug )

   allocate(v(dims(1), dims(2), dims(3)))

   call get_var_3d_real_cdf( wrf_3dvar_output_file, 'V', v, &
                             dims(1), dims(2), dims(3), 1, debug )

   if(debug) then
      do i=1,dims(1)
         write(unit=*, fmt='(2(a,i5), a, f12.8)') &
              'v(', i, ',', dims(2), ',1)=', v(i,dims(2),1)
      enddo
   endif

!--Get W
   call get_dims_cdf( wrf_3dvar_output_file, 'W', dims, ndims, debug )

   allocate(w(dims(1), dims(2), dims(3)))

   call get_var_3d_real_cdf( wrf_3dvar_output_file, 'W', w, &
                             dims(1), dims(2), dims(3), 1, debug )

   if(debug) then
      do i=1,dims(1)
         write(unit=*, fmt='(2(a,i5), a, f12.8)') &
              'w(', i, ',', dims(2), ',1)=', w(i,dims(2),1)
      enddo
   endif

   if(debug) then
      write(unit=*, fmt='(a,e20.12,4x)') &
           'Before couple Sampe u=', u(dims(1)/2,dims(2)/2,dims(3)/2), &
           'Before couple Sampe v=', v(dims(1)/2,dims(2)/2,dims(3)/2), &
           'Before couple Sampe w=', w(dims(1)/2,dims(2)/2,dims(3)/2)
   endif

!---------------------------------------------------------------------
!--Couple u, v, w.
   call couple_uvw ( u, v, w, mu, mub, msfu, msfv, msfm, ids, ide, jds, jde, kds, kde )

   if(debug) then
      write(unit=*, fmt='(a,e20.12,4x)') &
           'After  couple Sampe u=', u(dims(1)/2,dims(2)/2,dims(3)/2), &
           'After  couple Sampe v=', v(dims(1)/2,dims(2)/2,dims(3)/2)
   endif

!---------------------------------------------------------------------
!--For 3D variables

   do n=1,num3d
      call get_dims_cdf( wrf_3dvar_output_file, trim(var3d(n)), dims, ndims, debug )

      allocate(full3d(dims(1), dims(2), dims(3)))

      east_end=dims(1)+1
      north_end=dims(2)+1

      select case(trim(var3d(n)))
         case ('U') ;		! U
            var_pref='R' // trim(var3d(n))
            full3d(:,:,:)=u(:,:,:)
         case ('V') ;		! V
            var_pref='R' // trim(var3d(n))
            full3d(:,:,:)=v(:,:,:)
         case ('W') ;		! W
            var_pref='R' // trim(var3d(n))
            full3d(:,:,:)=w(:,:,:)
         case ('T', 'PH') ;
            var_pref=trim(var3d(n))

            call get_var_3d_real_cdf( wrf_3dvar_output_file, trim(var3d(n)), full3d, &
                                      dims(1), dims(2), dims(3), 1, debug )

            if(debug) then
               write(unit=*, fmt='(3a,e20.12,4x)') &
                    'Before couple Sampe ', trim(var3d(n)), &
                    '=', full3d(dims(1)/2,dims(2)/2,dims(3)/2)
            endif

            do j=1,dims(2)
            do k=1,dims(3)
            do i=1,dims(1)
               full3d(i,j,k)=full3d(i,j,k)*(mu(i,j)+mub(i,j))
            enddo
            enddo
            enddo

            if(debug) then
               write(unit=*, fmt='(3a,e20.12,4x)') &
                    'After  couple Sampe ', trim(var3d(n)), &
                    '=', full3d(dims(1)/2,dims(2)/2,dims(3)/2)
            endif
         case ('QVAPOR', 'QCLOUD', 'QRAIN', 'QICE', 'QSNOW', 'QGRAUPEL') ;
            var_pref='R' // var3d(n)(1:2)

            call get_var_3d_real_cdf( wrf_3dvar_output_file, trim(var3d(n)), full3d, &
                                      dims(1), dims(2), dims(3), 1, debug )

            if(debug) then
               write(unit=*, fmt='(3a,e20.12,4x)') &
                    'Before couple Sampe ', trim(var3d(n)), &
                    '=', full3d(dims(1)/2,dims(2)/2,dims(3)/2)
            endif

            do j=1,dims(2)
            do k=1,dims(3)
            do i=1,dims(1)
               full3d(i,j,k)=full3d(i,j,k)*(mu(i,j)+mub(i,j))
            enddo
            enddo
            enddo

            if(debug) then
               write(unit=*, fmt='(3a,e20.12,4x)') &
                    'After  couple Sampe ', trim(var3d(n)), &
                    '=', full3d(dims(1)/2,dims(2)/2,dims(3)/2)
            endif
         case default ;
            print *, 'It is impossible here. var3d(', n, ')=', trim(var3d(n))
      end select

      do m=1,4
         var_name=trim(var_pref) // trim(bdyname(m))
         call get_dims_cdf( wrf_bdy_file, trim(var_name), dims, ndims, debug )

         allocate(frst3d(dims(1), dims(2), dims(3)))
         allocate(scnd3d(dims(1), dims(2), dims(3)))
         allocate(tend3d(dims(1), dims(2), dims(3)))

!-----Calculate variable at second time level
!--------Get variable tendancy at first time level
        var_name=trim(var_pref) // trim(tenname(m))
        call get_var_3d_real_cdf( wrf_bdy_file, trim(var_name), tend3d, &
                                  dims(1), dims(2), dims(3), itime, debug )

!--------Get variable at first time level
        var_name=trim(var_pref) // trim(bdyname(m))
        call get_var_3d_real_cdf( wrf_bdy_file, trim(var_name), frst3d, &
                               dims(1), dims(2), dims(3), itime, debug )

       scnd3d = tend3d*bdyfrq_old + frst3d

!--------Get variable at second time level
!         call get_var_3d_real_cdf( wrf_bdy_file, trim(var_name), scnd3d, &
!                                   dims(1), dims(2), dims(3), 2, debug )

         select case(trim(bdyname(m)))
            case ('_BXS') ;		! West boundary
               do l=1,dims(3)
               do k=1,dims(2)
               do j=1,dims(1)
                  frst3d(j,k,l)=full3d(l,j,k)
               enddo
               enddo
               enddo
            case ('_BXE') ;		! East boundary
               do l=1,dims(3)
               do k=1,dims(2)
               do j=1,dims(1)
                  frst3d(j,k,l)=full3d(east_end-l,j,k)
               enddo
               enddo
               enddo
            case ('_BYS') ;		! South boundary
               do l=1,dims(3)
               do k=1,dims(2)
               do i=1,dims(1)
                  frst3d(i,k,l)=full3d(i,l,k)
               enddo
               enddo
               enddo
            case ('_BYE') ;		! North boundary
               do l=1,dims(3)
               do k=1,dims(2)
               do i=1,dims(1)
                  frst3d(i,k,l)=full3d(i,north_end-l,k)
               enddo
               enddo
               enddo
            case default ;
               print *, 'It is impossible here.'
               print *, 'bdyname(', m, ')=', trim(bdyname(m))
               stop
         end select

!--------calculate new tendancy

         tend3d = (scnd3d - frst3d)/bdyfrq

         if(debug) then
            write(unit=*, fmt='(a,i2,2x,2a/3(a,e20.12,4x)/a,i2,2x,a,4i6)') &
                 'No.', n, 'Variable: ', trim(var_name), &
                 'Sampe o frst3d=', frst3d(dims(1)/2,dims(2)/2,dims(3)/2), &
                         'scnd3d=', scnd3d(dims(1)/2,dims(2)/2,dims(3)/2), &
                         'tend3d=', tend3d(dims(1)/2,dims(2)/2,dims(3)/2), &
                 'ndims=', ndims, 'dims=', (dims(i), i=1,ndims)
         endif

         var_name=trim(var_pref) // trim(tenname(m))
         call put_var_3d_real_cdf( wrf_bdy_file, trim(var_name), tend3d, &
                                   dims(1), dims(2), dims(3), itime, debug )

         deallocate(frst3d)
         deallocate(scnd3d)
         deallocate(tend3d)
      enddo

      deallocate(full3d)
   enddo

   deallocate(mu)
   deallocate(mub)
   deallocate(msfu)
   deallocate(msfv)
   deallocate(msfm)
   deallocate(u)
   deallocate(v)
   deallocate(w)

   write(logfileunit,*)'FINISHED update_wrf_bc.'
   write(logfileunit,*)

   call finalize_utilities ! closes the log file.
 
   write(*,*) 'update_wrf_bc terminated normally.'

end program update_wrf_bc

