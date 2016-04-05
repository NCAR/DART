
      program main_bc_wrfchem

      use module_wrfchem_lib
      use module_mozart_lib
      use mo_calendar, only : diffdat
      use utils,       only : mapper
      use utils,       only : wrf2mz_map

      implicit none

!-----------------------------------------------------------------
!     parameters
!-----------------------------------------------------------------
      integer, parameter :: specmax = 500
      integer, parameter :: maxsize = 180
      integer, parameter :: kin = 15

!-----------------------------------------------------------------
!     control variables
!-----------------------------------------------------------------
      integer                :: nspec
      logical                :: do_bc = .false.            ! do boundary conditions
      logical                :: do_ic = .false.            ! do initial condition
      logical                :: def_missing_var = .false.  ! define missing variable in wrf data file
      character(len=maxsize) :: dir_moz
      character(len=maxsize) :: dir_wrf
      character(len=maxsize) :: init_cond_file_prefix = 'wrfinput'
      character(len=maxsize) :: bdy_cond_file_prefix  = 'wrfbdy'
      character(len=maxsize) :: fnb_wrf, fni_wrf
      character(len=maxsize) :: fn_moz
      character(len=3)       :: numa
      character(len=164)     :: spc_map(specmax)

      namelist /control/ dir_moz, dir_wrf, init_cond_file_prefix, bdy_cond_file_prefix, fn_moz, &
                         do_bc, do_ic, spc_map, moz_var_suffix, met_file_prefix, &
                         met_file_suffix, met_file_separator, surf_press_name, &
                         domain, def_missing_var

!-----------------------------------------------------------------
!     wrf variables 
!-----------------------------------------------------------------
      real :: p_top
      real, allocatable :: znu(:)
      real, allocatable :: xlon(:,:)
      real, allocatable :: xlat(:,:)

!-----------------------------------------------------------------
!     species vmr
!-----------------------------------------------------------------
      real, allocatable :: spc_ic(:,:,:)
      real, allocatable :: bcxs(:,:,:)
      real, allocatable :: bcxe(:,:,:)
      real, allocatable :: bcys(:,:,:)
      real, allocatable :: bcye(:,:,:)
!-----------------------------------------------------------------
!     species vmr old
!-----------------------------------------------------------------
      real, allocatable :: bcxso(:,:,:,:)
      real, allocatable :: bcxeo(:,:,:,:)
      real, allocatable :: bcyso(:,:,:,:)
      real, allocatable :: bcyeo(:,:,:,:)
!-----------------------------------------------------------------
!     species vmr tendencies
!-----------------------------------------------------------------
      real, allocatable :: btxs(:,:,:)
      real, allocatable :: btxe(:,:,:)
      real, allocatable :: btys(:,:,:)
      real, allocatable :: btye(:,:,:)
!-----------------------------------------------------------------
!     species vmr temp1
!-----------------------------------------------------------------
      real, allocatable :: txs(:,:,:)
      real, allocatable :: txe(:,:,:)
      real, allocatable :: tys(:,:,:)
      real, allocatable :: tye(:,:,:)

!-----------------------------------------------------------------
!     other working variables
!-----------------------------------------------------------------
      integer            :: ns, it
      integer            :: astat, istat
      real               :: tfac
      real               :: dt_wrfbdy
      real               :: rdiv
      real               :: f_int
      real               :: tinterp
      real, allocatable  :: ps_wrf(:,:)
      character(maxsize) :: fnb, fni
      logical            :: found
      logical            :: write_conc
      logical            :: write_tend

!-----------------------------------------------------------------
!     read control variables
!-----------------------------------------------------------------
      print *, 'APM: in mozbc'
      spc_map(:) = ' ' 
      read(*,nml=control,iostat=istat)
      if( istat /= 0 ) then
        write(*,*) 'main_bc_wrfchem: failed to read namelist; error = ',istat
        stop
      end if
!-----------------------------------------------------------------
!     check namelist inputs
!-----------------------------------------------------------------
      if( domain < 1 ) then
        write(*,*) 'main_bc_wrfchem: domain must be >= 1'
        stop
      end if
      write(numa,'(i3)') 100+domain
      domain_name = 'd' // numa(2:3)
      fni_wrf = trim( init_cond_file_prefix) // '_' // domain_name
      fnb_wrf = trim( bdy_cond_file_prefix) // '_' // domain_name
      if( .not. do_ic .and. .not. do_bc ) then
        write(*,*) 'main_bc_wrfchem: neither ic or bc requested'
        stop
      else
        if( do_bc ) then
          if( domain == 1 ) then
            write(*,*) 'main_bc_wrfchem: will set boundary conditions in ',trim(fnb_wrf)
          else
            write(*,*) 'main_bc_wrfchem: can only set boundary conditions for domain 1'
            if( .not. do_ic ) then
              write(*,*) 'main_bc_wrfchem: ic not requested and can not do bc'
              stop
            end if
            do_bc = .false.
          end if
        end if
        if( do_ic ) then
          write(*,*) 'main_bc_wrfchem: will set initial condition in ',trim(fni_wrf)
        end if
      end if
      call mapper( nspec, spc_map, specmax )
      write(*,*) 'main: nspec = ',nspec

!-----------------------------------------------------------------
!     intialize module_wrfchem_lib
!-----------------------------------------------------------------
      fnb  = trim(dir_wrf) // adjustl(fnb_wrf)
      fni  = trim(dir_wrf) // adjustl(fni_wrf)
      call init_wrfchem_lib( fnb, fni, dir_wrf, nspec, &
                             def_missing_var, do_bc, do_ic )

!-----------------------------------------------------------------
!     read eta values on half (mass) levels
!-----------------------------------------------------------------
      allocate( znu(nz),stat=astat )
      if( astat /= 0 ) then
        write(*,*) 'main_bc_wrfchem: failed to allocate znu; error = ',astat
        stop
      end if
      if( do_ic ) then
        call wrfchem_readscalar( 'P_TOP', p_top )
        call wrfchem_read2d( 'ZNU', znu )
      else
        call wrfchem_readscalar( 'P_TOP', p_top, trim(fni) )
        call wrfchem_read2d( 'ZNU', znu, trim(fni) )
      endif
      write(*,*) 'main_bc_wrfchem: read p_top'
      write(*,*) 'main_bc_wrfchem: read eta values on half (mass) levels'

!-----------------------------------------------------------------
!     read wrf longitudes and latitudes
!-----------------------------------------------------------------
      allocate( xlon(nx,ny), xlat(nx,ny), stat=astat )
      if( astat /= 0 ) then
        write(*,*) 'main_bc_wrfchem: failed to allocate xlon,xlat; error = ',astat
        stop
      end if
      if( do_ic ) then
        call wrfchem_read3d( 'XLONG', xlon )
        call wrfchem_read3d( 'XLAT', xlat )
      else
        call wrfchem_read3d( 'XLONG', xlon, trim(fni) )
        call wrfchem_read3d( 'XLAT', xlat, trim(fni) )
      end if
      write(*,*) 'main_bc_wrfchem: read wrf longitues and latitudes'

!-----------------------------------------------------------------
!     initialize module_mozart_lib
!-----------------------------------------------------------------
      fni = trim(dir_moz) // adjustl(fn_moz)
      call init_mozart_lib( dir_moz, fn_moz, xlon, xlat, wrf_date(1), &
                            wrf_datesec(1), nx, ny, nspec )

!-----------------------------------------------------------------
!     allocate arrays
!-----------------------------------------------------------------
      istat = 0
      if( do_ic ) then
        allocate( spc_ic(nx,ny,nz),stat=astat )
      end if
      istat = istat + astat
      allocate( bcxs(ny,nz,nw),stat=astat )
      istat = istat + astat
      allocate( bcxe(ny,nz,nw),stat=astat )
      istat = istat + astat
      allocate( bcys(nx,nz,nw),stat=astat )
      istat = istat + astat
      allocate( bcye(nx,nz,nw),stat=astat )
      istat = istat + astat
      allocate( bcxso(ny,nz,nw,nspec),stat=astat )
      istat = istat + astat
      allocate( bcxeo(ny,nz,nw,nspec),stat=astat )
      istat = istat + astat
      allocate( bcyso(nx,nz,nw,nspec),stat=astat )
      istat = istat + astat
      allocate( bcyeo(nx,nz,nw,nspec),stat=astat )
      istat = istat + astat
      allocate( btxs(ny,nz,nw),stat=astat )
      istat = istat + astat
      allocate( btxe(ny,nz,nw),stat=astat )
      istat = istat + astat
      allocate( btys(nx,nz,nw),stat=astat )
      istat = istat + astat
      allocate( btye(nx,nz,nw),stat=astat )
      istat = istat + astat
      allocate( txs(ny,nz,nw),stat=astat )
      istat = istat + astat
      allocate( txe(ny,nz,nw),stat=astat )
      istat = istat + astat
      allocate( tys(nx,nz,nw),stat=astat )
      istat = istat + astat
      allocate( tye(nx,nz,nw),stat=astat )
      istat = istat + astat
      allocate( ps_wrf(nx,ny),stat=astat )
      istat = istat + astat
      if( istat /= 0 ) then
        write(*,*) 'main_bc_wrfchem: failed to allocate bcxs ... ps_wrf; error = ',istat
        stop
      end if

      if( ntime > 1 ) then
        dt_wrfbdy = 86400. * diffdat( wrf_date(1), wrf_datesec(1), wrf_date(2), wrf_datesec(2) )
        write(*,*) 'dt_wrfbdy = ',dt_wrfbdy
      end if

!-----------------------------------------------------------------
!     read mozart data, interpolate and then write BC file
!-----------------------------------------------------------------
time_loop : &
      do it = 1,ntime
        write(*,*) 'main; doing time ',wrf_date(it),' : ',wrf_datesec(it)
        call get_moz_time_ndx( dir_moz, fn_moz, wrf_date(it), wrf_datesec(it), nspec )
        call wrf_ps( dir_wrf, wrf_date(it), wrf_datesec(it), ps_wrf )
species_loop : &
        do ns = 1, nspec
!-----------------------------------------------------------------
!      IC
!-----------------------------------------------------------------
          if( do_ic .and. it == 1 ) then
            write(*,*) 'main: setting initial condition for ',trim(wrf2mz_map(ns)%wrf_name)
            call ic_interpolate4d( ns, spc_ic, ps_wrf, znu, p_top, &
                                   nx, ny, nz )
            call wrfchem_ic_write4d( ns, spc_ic, it )
          end if
has_bndy_cnd : &
          if( do_bc ) then
            write(*,*) 'main: setting boundary condition for ',trim(wrf2mz_map(ns)%wrf_name)
            call bc_interpolate4d( ns, bcxs, bcxe, bcys, bcye, &
                                   ps_wrf, znu, p_top, nx, ny, &
                                   nz, nw )
!-----------------------------------------------------------------
!      tendency
!-----------------------------------------------------------------
            if( it >= 2 ) then
              rdiv = 1./dt_wrfbdy
              btxs(:,:,:) = (bcxs(:,:,:) - bcxso(:,:,:,ns))*rdiv
              btxe(:,:,:) = (bcxe(:,:,:) - bcxeo(:,:,:,ns))*rdiv
              btys(:,:,:) = (bcys(:,:,:) - bcyso(:,:,:,ns))*rdiv
              btye(:,:,:) = (bcye(:,:,:) - bcyeo(:,:,:,ns))*rdiv
            end if
            if( it < ntime ) then
              bcxso(:,:,:,ns) = bcxs(:,:,:)
              bcxeo(:,:,:,ns) = bcxe(:,:,:)
              bcyso(:,:,:,ns) = bcys(:,:,:)
              bcyeo(:,:,:,ns) = bcye(:,:,:)
            end if
              write_conc = it < ntime
              write_tend = it > 1
              call wrfchem_bc_write4d( ns, bcxs, bcxe, bcys, bcye, &
                                       btxs, btxe, btys, btye, it, write_conc, write_tend )
          end if has_bndy_cnd
        end do species_loop
      end do time_loop

!-----------------------------------------------------------------
!     exit from libs
!-----------------------------------------------------------------
      call exit_wrfchem_lib
      call exit_mozart_lib

!-----------------------------------------------------------------
!     deallocate memory space
!-----------------------------------------------------------------
      if( do_ic ) then
        deallocate( spc_ic )
      end if
      deallocate( znu, xlon, xlat, bcxs, bcxe, bcys, bcye,    &
                  bcxso, bcxeo, bcyso, bcyeo, &
                  btxs, btxe, btys, btye, &
                  txs, txe, tys, tye )

      write(*,*) 'bc_wrfchem completed successfully'

!-----------------------------------------------------------------
!     create success signal file
!-----------------------------------------------------------------
      open(unit=10,file='run_mozbc_done',iostat=istat)
      if( istat /= 0 ) then
        write(*,*) 'main_bc_wrfchem: failed to open signal file; error = ',istat
        stop 'Opn Err'
      end if
      close(10)

    end program main_bc_wrfchem
