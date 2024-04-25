! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program update_wrf_bc

! program to update BC file from 3dvar or filter output.
! current version reads only wrf-netcdf file format.

! Input files: wrfinput_d01, wrfbdy_d01, and wrfinput_mean (optional).

use               types_mod, only : r8
use           utilities_mod, only : file_exist, open_file, close_file, &
                                    initialize_utilities, finalize_utilities, register_module, &
                                    error_handler, E_ERR, E_MSG, timestamp
use module_netcdf_interface, only : get_dims_cdf, get_gl_att_real_cdf, put_gl_att_real_cdf, &
                                    get_var_3d_real_cdf, get_var_2d_real_cdf, put_var_3d_real_cdf, &
                                    put_var_2d_real_cdf, get_times_cdf, put_time_cdf, variable_exist
use        module_couple_uv
use         module_timediff, only : time_diff, find_time_index

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!-----------------------------------------------------------------------
! Model namelist parameters with default values.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------

integer, parameter :: max_3d_variables = 20, &
                      max_2d_variables = 20, &
                      max_times        = 100

character(len=80) :: wrf_output_file, wrf_mean_output_file, &
                     wrf_bdy_file, time_name

character(len=20) :: var_name

character(len=20) :: var3d(max_3d_variables), &
                     var2d(max_2d_variables)

character(len=10), dimension(4) :: bdyname, tenname

character(len=19), dimension(max_times) :: udtime, bdytime, thisbdytime, nextbdytime

character(len=129) :: msgstring

integer           :: ids, ide, jds, jde, kds, kde
integer           :: num3d, num2d, ndims, nmoist
integer           :: i,j,k,l,m,n
integer           :: ntimes_bdy, ntimes_ud, itime

integer, dimension(4) :: dims

real(r8), allocatable, dimension(:,:) :: tend2d, scnd2d, frst2d

real(r8), allocatable, dimension(:,:,:) :: tend3d, scnd3d, frst3d, full3d, full3d_mean

real(r8), allocatable, dimension(:,:,:) :: u, v, w
real(r8), allocatable, dimension(:,:,:) :: u_mean, v_mean, w_mean

real(r8), allocatable, dimension(:,  :) :: mu, mu_mean, mub, msfu, msfv, msfm

integer :: east_end, north_end

logical, parameter :: debug = .false.

real(r8) :: bdyfrq_old, bdyfrq, infl

!----------------------------------------------------------------------

call initialize_utilities('update_wrf_bc')
call register_module(source, revision, revdate)

wrf_output_file='wrfinput_d01'
wrf_bdy_file  ='wrfbdy_d01'

if((.not. file_exist(wrf_output_file)) .or. (.not. file_exist(wrf_bdy_file))) then
   write(msgstring, *)'wrfinput_d01 or wrfbdy_d01 is absent.'
   call error_handler(E_ERR,'update_wrf_bc',msgstring,source,revision,revdate)
endif

! If the file 'wrfinput_mean' exists,
! input the ensemble mean to calculate deviations from the mean.
! These perturbations are then added to the fields at the end of the interval.

if(file_exist('wrfinput_mean')) then
   wrf_mean_output_file='wrfinput_mean'
   write(*,*) 'Input real coefficient multiplying (wrfinput_d01 - wrfinput_mean):'
   write(*,*) 'Result will be added to the fields at the end of the interval.'
   read(*,*) infl
else
   wrf_mean_output_file='wrfinput_d01'
   infl = 0.0_r8
endif

!--boundary variables
bdyname(1)='_BXS'
bdyname(2)='_BXE'
bdyname(3)='_BYS'
bdyname(4)='_BYE'

!--boundary tendency variables
tenname(1)='_BTXS'
tenname(2)='_BTXE'
tenname(3)='_BTYS'
tenname(4)='_BTYE'

!--3D need update
nmoist = 7 
num3d  = 5 + nmoist
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
var3d(11)='QGRAUP'
var3d(12)='QNICE'

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
time_name = 'Times'
call get_times_cdf( wrf_output_file, time_name, udtime, ntimes_ud, max_times, debug )

if(debug) print*, 'udtime = ',udtime(1)

!-- list of boundary times
call get_times_cdf( wrf_bdy_file, time_name, bdytime, ntimes_bdy, max_times, debug )

!-- Time index of current tendency - to grab from wrfbdy
call find_time_index(bdytime, udtime(1), ntimes_bdy, itime)

if(debug) print*, 'itime = ',itime
if(debug) print*, 'bdytime = ',bdytime(itime)

!---------------------------------------------------------------------
!--First, the boundary frequency.
time_name = 'md___thisbdytimee_x_t_d_o_m_a_i_n_m_e_t_a_data_'
call get_times_cdf( wrf_bdy_file, time_name, &
     thisbdytime, ntimes_bdy, max_times, debug )
time_name = 'md___nextbdytimee_x_t_d_o_m_a_i_n_m_e_t_a_data_'
call get_times_cdf( wrf_bdy_file, time_name, &
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
time_name = 'md___thisbdytimee_x_t_d_o_m_a_i_n_m_e_t_a_data_'
call put_time_cdf( wrf_bdy_file, time_name, &
     thisbdytime(itime), itime, debug )

!--For 2D variables
!--Get mu, mub, msfu, msfv, and msfm

do n=1,num2d
   call get_dims_cdf( wrf_output_file, trim(var2d(n)), dims, ndims, debug )

   select case(trim(var2d(n)))
   case ('MU') ;
      allocate(mu(dims(1), dims(2)))

      call get_var_2d_real_cdf( wrf_output_file, trim(var2d(n)), mu, &
           dims(1), dims(2), 1, debug )

      allocate(mu_mean(dims(1), dims(2)))

      call get_var_2d_real_cdf( wrf_mean_output_file, trim(var2d(n)), mu_mean, &
           dims(1), dims(2), 1, debug )

      east_end=dims(1)+1
      north_end=dims(2)+1
   case ('MUB') ;
      allocate(mub(dims(1), dims(2)))

      call get_var_2d_real_cdf( wrf_output_file, trim(var2d(n)), mub, &
           dims(1), dims(2), 1, debug )
   case ('MAPFAC_U') ;
      allocate(msfu(dims(1), dims(2)))

      call get_var_2d_real_cdf( wrf_output_file, trim(var2d(n)), msfu, &
           dims(1), dims(2), 1, debug )

      if(debug) then
         do j=1,dims(2)
            write(unit=*, fmt='(2(a,i5), a, f12.8)') &
                 'msfu(', dims(1), ',', j, ')=', msfu(dims(1),j)
         enddo
      endif

   case ('MAPFAC_V') ;
      allocate(msfv(dims(1), dims(2)))

      call get_var_2d_real_cdf( wrf_output_file, trim(var2d(n)), msfv, &
           dims(1), dims(2), 1, debug )

      if(debug) then
         do i=1,dims(1)
            write(unit=*, fmt='(2(a,i5), a, f12.8)') &
                 'msfv(', i, ',', dims(2), ')=', msfv(i,dims(2))
         enddo
      endif

   case ('MAPFAC_M') ;
      allocate(msfm(dims(1), dims(2)))

      call get_var_2d_real_cdf( wrf_output_file, trim(var2d(n)), msfm, &
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
!--------Get variable tendency
      var_name='MU' // trim(tenname(m))
      call get_var_2d_real_cdf( wrf_bdy_file, trim(var_name), tend2d, &
           dims(1), dims(2), itime, debug )

!--------Get variable at first time level
      var_name='MU' // trim(bdyname(m))
      call get_var_2d_real_cdf( wrf_bdy_file, trim(var_name), frst2d, &
           dims(1), dims(2), itime, debug )

      scnd2d = tend2d*bdyfrq_old + frst2d

!--------Get variable at second time level
!      call get_var_2d_real_cdf( wrf_bdy_file, trim(var_name), scnd2d, &
!                                dims(1), dims(2), 2, debug )



!-----Add BC perturbation at second time level
      select case(m)
         case (1) ;             ! West boundary
            do l=1,dims(2)
            do j=1,dims(1)
               frst2d(j,l)=mu(l,j)
               scnd2d(j,l) = scnd2d(j,l) + infl*(mu(l,j)-mu_mean(l,j))
            enddo
            enddo
         case (2) ;             ! East boundary
            do l=1,dims(2)
            do j=1,dims(1)
               frst2d(j,l)=mu(east_end-l,j)
               scnd2d(j,l) = scnd2d(j,l) + infl*(mu(east_end-l,j)-mu_mean(east_end-l,j))
            enddo
            enddo
         case (3) ;             ! South boundary
            do l=1,dims(2)
            do i=1,dims(1)
               frst2d(i,l)=mu(i,l)
               scnd2d(i,l) = scnd2d(i,l) + infl*(mu(i,l)-mu_mean(i,l))
            enddo
            enddo
         case (4) ;             ! North boundary
            do l=1,dims(2)
            do i=1,dims(1)
               frst2d(i,l)=mu(i,north_end-l)
               scnd2d(i,l) = scnd2d(i,l) + infl*(mu(i,north_end-l)-mu_mean(i,north_end-l))
            enddo
            enddo
         case default ;
            print *, 'It is impossible here. mu, m=', m
      end select

!-----calculate new tendency 

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
!-----output new tendency 
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
   call get_dims_cdf( wrf_output_file, 'U', dims, ndims, debug )

   allocate(u(dims(1), dims(2), dims(3)))

   ids=1
   ide=dims(1)-1
   jds=1
   jde=dims(2)
   kds=1
   kde=dims(3)

   call get_var_3d_real_cdf( wrf_output_file, 'U', u, &
                             dims(1), dims(2), dims(3), 1, debug )

   allocate(u_mean(dims(1), dims(2), dims(3)))

   call get_var_3d_real_cdf( wrf_mean_output_file, 'U', u_mean, &
                             dims(1), dims(2), dims(3), 1, debug )

   if(debug) then
      do j=1,dims(2)
         write(unit=*, fmt='(2(a,i5), a, f12.8)') &
              'u(', dims(1), ',', j, ',1)=', u(dims(1),j,1)
      enddo
   endif

!--Get V
   call get_dims_cdf( wrf_output_file, 'V', dims, ndims, debug )

   allocate(v(dims(1), dims(2), dims(3)))

   call get_var_3d_real_cdf( wrf_output_file, 'V', v, &
                             dims(1), dims(2), dims(3), 1, debug )

   allocate(v_mean(dims(1), dims(2), dims(3)))

   call get_var_3d_real_cdf( wrf_mean_output_file, 'V', v_mean, &
                             dims(1), dims(2), dims(3), 1, debug )

   if(debug) then
      do i=1,dims(1)
         write(unit=*, fmt='(2(a,i5), a, f12.8)') &
              'v(', i, ',', dims(2), ',1)=', v(i,dims(2),1)
      enddo
   endif

!--Get W
   call get_dims_cdf( wrf_output_file, 'W', dims, ndims, debug )

   allocate(w(dims(1), dims(2), dims(3)))

   call get_var_3d_real_cdf( wrf_output_file, 'W', w, &
                             dims(1), dims(2), dims(3), 1, debug )

   allocate(w_mean(dims(1), dims(2), dims(3)))

   call get_var_3d_real_cdf( wrf_mean_output_file, 'W', w_mean, &
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

   call couple_uvw ( u_mean, v_mean, w_mean, mu_mean, mub, msfu, msfv, msfm, ids, ide, jds, jde, kds, kde )

   if(debug) then
      write(unit=*, fmt='(a,e20.12,4x)') &
           'After  couple Sampe u=', u(dims(1)/2,dims(2)/2,dims(3)/2), &
           'After  couple Sampe v=', v(dims(1)/2,dims(2)/2,dims(3)/2)
   endif

!---------------------------------------------------------------------
!--For 3D variables

   loop3d : do n=1,num3d

      if ( .not. variable_exist(wrf_output_file, trim(var3d(n))) ) cycle loop3d

      call get_dims_cdf( wrf_output_file, trim(var3d(n)), dims, ndims, debug )

      allocate(full3d(dims(1), dims(2), dims(3)))

      allocate(full3d_mean(dims(1), dims(2), dims(3)))

      east_end=dims(1)+1
      north_end=dims(2)+1

      select case(trim(var3d(n)))
         case ('U') ;           ! U
            full3d(:,:,:)=u(:,:,:)
            full3d_mean(:,:,:)=u_mean(:,:,:)
         case ('V') ;           ! V
            full3d(:,:,:)=v(:,:,:)
            full3d_mean(:,:,:)=v_mean(:,:,:)
         case ('W') ;           ! W
            full3d(:,:,:)=w(:,:,:)
            full3d_mean(:,:,:)=w_mean(:,:,:)
         case ('T', 'PH', 'QVAPOR', 'QCLOUD', 'QRAIN', 'QICE', 'QSNOW', 'QGRAUP', 'QNICE') ;

            call get_var_3d_real_cdf( wrf_output_file, trim(var3d(n)), full3d, &
                                      dims(1), dims(2), dims(3), 1, debug )

            call get_var_3d_real_cdf( wrf_mean_output_file, trim(var3d(n)), full3d_mean, &
                                      dims(1), dims(2), dims(3), 1, debug )

            if(debug) then
               write(unit=*, fmt='(3a,e20.12,4x)') &
                    'Before couple Sampe ', trim(var3d(n)), &
                    '=', full3d(dims(1)/2,dims(2)/2,dims(3)/2)
            endif

            do k=1,dims(3)
            do j=1,dims(2)
            do i=1,dims(1)
               full3d(i,j,k)=full3d(i,j,k)*(mu(i,j)+mub(i,j))
               full3d_mean(i,j,k)=full3d_mean(i,j,k)*(mu_mean(i,j)+mub(i,j))
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
         var_name=trim(var3d(n)) // trim(bdyname(m))
         call get_dims_cdf( wrf_bdy_file, trim(var_name), dims, ndims, debug )

         allocate(frst3d(dims(1), dims(2), dims(3)))
         allocate(scnd3d(dims(1), dims(2), dims(3)))
         allocate(tend3d(dims(1), dims(2), dims(3)))

!-----Calculate variable at second time level
!--------Get variable tendency
        var_name=trim(var3d(n)) // trim(tenname(m))
        call get_var_3d_real_cdf( wrf_bdy_file, trim(var_name), tend3d, &
                                  dims(1), dims(2), dims(3), itime, debug )

!--------Get variable at first time level
        var_name=trim(var3d(n)) // trim(bdyname(m))
        call get_var_3d_real_cdf( wrf_bdy_file, trim(var_name), frst3d, &
                               dims(1), dims(2), dims(3), itime, debug )

       scnd3d = tend3d*bdyfrq_old + frst3d

!--------Get variable at second time level
!         call get_var_3d_real_cdf( wrf_bdy_file, trim(var_name), scnd3d, &
!                                   dims(1), dims(2), dims(3), 2, debug )

!-----Add BC perturbation at second time level
         select case(trim(bdyname(m)))
            case ('_BXS') ;             ! West boundary
               do l=1,dims(3)
               do k=1,dims(2)
               do j=1,dims(1)
                  frst3d(j,k,l)=full3d(l,j,k)
                  scnd3d(j,k,l) = scnd3d(j,k,l) + infl*(full3d(l,j,k)-full3d_mean(l,j,k))
               enddo
               enddo
               enddo
            case ('_BXE') ;             ! East boundary
               do l=1,dims(3)
               do k=1,dims(2)
               do j=1,dims(1)
                  frst3d(j,k,l)=full3d(east_end-l,j,k)
                  scnd3d(j,k,l) = scnd3d(j,k,l) + infl*(full3d(east_end-l,j,k)-full3d_mean(east_end-l,j,k))
               enddo
               enddo
               enddo
            case ('_BYS') ;             ! South boundary
               do l=1,dims(3)
               do k=1,dims(2)
               do i=1,dims(1)
                  frst3d(i,k,l)=full3d(i,l,k)
                  scnd3d(i,k,l) = scnd3d(i,k,l) + infl*(full3d(i,l,k) - full3d_mean(i,l,k))
               enddo
               enddo
               enddo
            case ('_BYE') ;             ! North boundary
               do l=1,dims(3)
               do k=1,dims(2)
               do i=1,dims(1)
                  frst3d(i,k,l)=full3d(i,north_end-l,k)
                  scnd3d(i,k,l) = scnd3d(i,k,l) + infl*(full3d(i,north_end-l,k)-full3d_mean(i,north_end-l,k))
               enddo
               enddo
               enddo
            case default ;
               print *, 'It is impossible here.'
               print *, 'bdyname(', m, ')=', trim(bdyname(m))
               stop
         end select

!--------Make sure that microphysics variables at the end of the interval are not negatives.

         if(n >= 6) then
            scnd3d(:,:,:) = max(0.0_r8,scnd3d(:,:,:))
         endif

!--------calculate new tendency

         tend3d = (scnd3d - frst3d)/bdyfrq

         if(debug) then
            write(unit=*, fmt='(a,i2,2x,2a/3(a,e20.12,4x)/a,i2,2x,a,4i6)') &
                 'No.', n, 'Variable: ', trim(var_name), &
                 'Sampe o frst3d=', frst3d(dims(1)/2,dims(2)/2,dims(3)/2), &
                         'scnd3d=', scnd3d(dims(1)/2,dims(2)/2,dims(3)/2), &
                         'tend3d=', tend3d(dims(1)/2,dims(2)/2,dims(3)/2), &
                 'ndims=', ndims, 'dims=', (dims(i), i=1,ndims)
         endif

!-----output new variable at first time level
         var_name=trim(var3d(n)) // trim(bdyname(m))
         call put_var_3d_real_cdf( wrf_bdy_file, trim(var_name), frst3d, &
                                   dims(1), dims(2), dims(3), itime, debug )

!-----output new tendency 
         var_name=trim(var3d(n)) // trim(tenname(m))
         call put_var_3d_real_cdf( wrf_bdy_file, trim(var_name), tend3d, &
                                   dims(1), dims(2), dims(3), itime, debug )

         deallocate(frst3d)
         deallocate(scnd3d)
         deallocate(tend3d)
      enddo

      deallocate(full3d)
      deallocate(full3d_mean)
   enddo loop3d

   deallocate(mu)
   deallocate(mu_mean)
   deallocate(mub)
   deallocate(msfu)
   deallocate(msfv)
   deallocate(msfm)
   deallocate(u)
   deallocate(v)
   deallocate(w)
   deallocate(u_mean)
   deallocate(v_mean)
   deallocate(w_mean)

   call error_handler(E_MSG, 'update_wrf_bc', 'update_wrf_bc terminated normally.')
   call error_handler(E_MSG, 'update_wrf_bc', 'FINISHED update_wrf_bc.')
   call error_handler(E_MSG, 'update_wrf_bc', 'Finished successfully.',source,revision,revdate)
   call finalize_utilities()
 
end program update_wrf_bc

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
