! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

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

  USE module_netcdf_interface
  USE module_couple_uv

   implicit none

   integer, parameter :: max_3d_variables = 20, &
                         max_2d_variables = 20
 
   character(len=80) :: wrf_3dvar_output_file, &
                        wrf_bdy_file
 
   character(len=20) :: var_pref, var_name

   character(len=20) :: var3d(max_3d_variables), &
                        var2d(max_2d_variables)

   character(len=10), dimension(4) :: bdyname, tenname

   integer           :: ids, ide, jds, jde, kds, kde
   integer           :: num3d, num2d, ndims
   integer           :: i,j,k,l,m,n

   integer, dimension(4) :: dims
 
   integer, external :: iargc
 
   real, allocatable, dimension(:,:,:) :: tend3d, scnd3d, frst3d, full3d

   real, allocatable, dimension(:,:,:) :: u, v

   real, allocatable, dimension(:,  :) :: mu, mub, msfu, msfv
 
   integer :: east_end, north_end

   logical :: debug

   real :: bdyfrq

   debug = .false. 

!---------------------------------------------------------------------
   wrf_3dvar_output_file='wrfinput_d01'
   wrf_bdy_file  ='wrfbdy_d01'
!---------------------------------------------------------------------
!  if (iargc() /= 2) then
!     write(6,*) ' usage: update_wrf_bc wrf_3dvar_output_file wrf_bdy_file '
!     stop
!  else
!     call getarg(1,wrf_3dvar_output_file)
!     call getarg(2,wrf_bdy_file)

!     if(debug) then
!        write(unit=*, fmt='(2a)') &
!             'input file is ', trim(wrf_3dvar_output_file), &
!             'bdy   file is ', trim(wrf_bdy_file)
!     endif
!  endif
!---------------------------------------------------------------------

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
!---------------------------------------------------------------------

!--3D need update
   num3d=5
   var3d(1)='U'
   var3d(2)='V'
   var3d(3)='T'
   var3d(4)='PH'
   var3d(5)='QVAPOR'

!--2D need update
   num2d=4
   var2d(1)='MUB'
   var2d(2)='MU'
   var2d(3)='MAPFAC_U'
   var2d(4)='MAPFAC_V'

!---------------------------------------------------------------------
   east_end=0
   north_end=0
!---------------------------------------------------------------------
  
!---------------------------------------------------------------------
!--First, the boundary frequency.
   call get_gl_att_real_cdf( wrf_bdy_file, 'BDYFRQ', bdyfrq, debug )

   if(debug) then
      write(unit=*, fmt='(a, f12.2)') &
             'BDYFRQ=', bdyfrq
   endif
!---------------------------------------------------------------------

!--For 2D variables
!--Get mu, mub, msfu, and msfv

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

            do j=1,dims(2)
               write(unit=*, fmt='(2(a,i5), a, f12.8)') &
                    'msfu(', dims(1), ',', j, ')=', msfu(dims(1),j)
            enddo

         case ('MAPFAC_V') ;
            allocate(msfv(dims(1), dims(2)))

            call get_var_2d_real_cdf( wrf_3dvar_output_file, trim(var2d(n)), msfv, &
                                      dims(1), dims(2), 1, debug )

            do i=1,dims(1)
               write(unit=*, fmt='(2(a,i5), a, f12.8)') &
                    'msfv(', i, ',', dims(2), ')=', msfv(i,dims(2))
            enddo

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

      allocate(frst3d(dims(1), dims(2), dims(3)))
      allocate(scnd3d(dims(1), dims(2), dims(3)))
      allocate(tend3d(dims(1), dims(2), dims(3)))

!-----Calculate variable at second time level
!--------Get variable tendancy at first time level
        var_name='MU' // trim(tenname(m))
        call get_var_3d_real_cdf( wrf_bdy_file, trim(var_name), tend3d, &
                                  dims(1), dims(2), dims(3), 1, debug )

!--------Get variable at first time level
        var_name='MU' // trim(bdyname(m))
        call get_var_3d_real_cdf( wrf_bdy_file, trim(var_name), frst3d, &
                               dims(1), dims(2), dims(3), 1, debug )

       scnd3d = tend3d*bdyfrq + frst3d

!      call get_var_3d_real_cdf( wrf_bdy_file, trim(var_name), scnd3d, &
!                                dims(1), dims(2), dims(3), 2, debug )



!-----calculate variable at first time level
      select case(m)
         case (1) ;		! West boundary
            do l=1,dims(3)
            do j=1,dims(1)
               frst3d(j,1,l)=mu(l,j)
            enddo
            enddo
         case (2) ;		! East boundary
            do l=1,dims(3)
            do j=1,dims(1)
               frst3d(j,1,l)=mu(east_end-l,j)
            enddo
            enddo
         case (3) ;		! South boundary
            do l=1,dims(3)
            do i=1,dims(1)
               frst3d(i,1,l)=mu(i,l)
            enddo
            enddo
         case (4) ;		! North boundary
            do l=1,dims(3)
            do i=1,dims(1)
               frst3d(i,1,l)=mu(i,north_end-l)
            enddo
            enddo
         case default ;
            print *, 'It is impossible here. mu, m=', m
      end select

!-----calculate new tendancy 
      do l=1,dims(3)
      do i=1,dims(1)
         tend3d(i,1,l)=(scnd3d(i,1,l)-frst3d(i,1,l))/bdyfrq
      enddo
      enddo

      if(debug) then
         write(unit=*, fmt='(a,i2,2x,2a/3(a,e20.12,4x)/a,i2,2x,a,4i6)') &
              'No.', n, 'Variable: ', trim(var_name), &
              'Sampe frst3d=', frst3d(dims(1)/2,1,dims(3)/2), &
                    'scnd3d=', scnd3d(dims(1)/2,1,dims(3)/2), &
                    'tend3d=', tend3d(dims(1)/2,1,dims(3)/2), &
              'ndims=', ndims, 'dims=', (dims(i), i=1,ndims)
      endif

!-----output new variable at first time level
      var_name='MU' // trim(bdyname(m))
      call put_var_3d_real_cdf( wrf_bdy_file, trim(var_name), frst3d, &
                                dims(1), dims(2), dims(3), 1, debug )
!-----output new tendancy 
      var_name='MU' // trim(tenname(m))
      call put_var_3d_real_cdf( wrf_bdy_file, trim(var_name), tend3d, &
                                dims(1), dims(2), dims(3), 1, debug )

      deallocate(frst3d)
      deallocate(scnd3d)
      deallocate(tend3d)
   enddo

!---------------------------------------------------------------------
!--For 3D variables

!--Get U
   call get_dims_cdf( wrf_3dvar_output_file, 'U', dims, ndims, debug )

!  call get_att_cdf( wrf_3dvar_output_file, 'U', debug )

   allocate(u(dims(1), dims(2), dims(3)))

   ids=1
   ide=dims(1)-1
   jds=1
   jde=dims(2)
   kds=1
   kde=dims(3)

   call get_var_3d_real_cdf( wrf_3dvar_output_file, 'U', u, &
                             dims(1), dims(2), dims(3), 1, debug )

   do j=1,dims(2)
      write(unit=*, fmt='(2(a,i5), a, f12.8)') &
           'u(', dims(1), ',', j, ',1)=', u(dims(1),j,1)
   enddo

!--Get V
   call get_dims_cdf( wrf_3dvar_output_file, 'V', dims, ndims, debug )

!  call get_att_cdf( wrf_3dvar_output_file, 'V', debug )

   allocate(v(dims(1), dims(2), dims(3)))

   call get_var_3d_real_cdf( wrf_3dvar_output_file, 'V', v, &
                             dims(1), dims(2), dims(3), 1, debug )

   do i=1,dims(1)
      write(unit=*, fmt='(2(a,i5), a, f12.8)') &
           'v(', i, ',', dims(2), ',1)=', v(i,dims(2),1)
   enddo

   if(debug) then
      write(unit=*, fmt='(a,e20.12,4x)') &
           'Before couple Sampe u=', u(dims(1)/2,dims(2)/2,dims(3)/2), &
           'Before couple Sampe v=', v(dims(1)/2,dims(2)/2,dims(3)/2)
   endif

!---------------------------------------------------------------------
!--Couple u, v.
   call couple_uv ( u, v, mu, mub, msfu, msfv, ids, ide, jds, jde, kds, kde )

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
         case ('QVAPOR') ;
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
                                  dims(1), dims(2), dims(3), 1, debug )

!--------Get variable at first time level
        var_name=trim(var_pref) // trim(bdyname(m))
        call get_var_3d_real_cdf( wrf_bdy_file, trim(var_name), frst3d, &
                               dims(1), dims(2), dims(3), 1, debug )

       scnd3d = tend3d*bdyfrq + frst3d

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
         do l=1,dims(3)
         do k=1,dims(2)
         do i=1,dims(1)
            tend3d(i,k,l)=(scnd3d(i,k,l)-frst3d(i,k,l))/bdyfrq
         enddo
         enddo
         enddo
   
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
                                   dims(1), dims(2), dims(3), 1, debug )

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
   deallocate(u)
   deallocate(v)

end program update_wrf_bc

