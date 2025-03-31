!  SVN:$Id: ice_restart_driver.F90 607 2013-03-29 15:49:42Z eclare $
!=======================================================================

! Read and write ice model restart files
!
! authors Elizabeth C. Hunke, LANL
!         William H. Lipscomb LANL
!         David Bailey, NCAR
!
! 2004-05: Block structure added by William Lipscomb
!          Restart module separated from history module
! 2006 ECH: Accepted some CCSM code into mainstream CICE
!           Converted to free source form (F90) 
! 2008 ECH: Rearranged order in which internal stresses are written and read
! 2010 ECH: Changed eice, esno to qice, qsno
! 2012 ECH: Added routines for reading/writing extended grid
! 2013 DAB: Added generic interfaces for writing restart fields.

      module ice_restart_driver

      use ice_kinds_mod
      use ice_restart_shared, only: &
          restart, restart_ext, restart_dir, restart_file, pointer_file, &
          runid, runtype, use_restart_time, restart_format, lcdf64, lenstr
      use ice_restart

      implicit none
      private
      public :: dumpfile, restartfile, restartfile_v4
      save

!=======================================================================

      contains

!=======================================================================

!=======================================================================
!---! these subroutines write/read Fortran unformatted data files ..
!=======================================================================

! Dumps all values needed for a restart
! author Elizabeth C. Hunke, LANL

      subroutine dumpfile(filename_spec)

      use ice_blocks, only: nx_block, ny_block
      use ice_calendar, only: sec, month, mday, nyr, istep1, &
                              time, time_forc, year_init
      use ice_communicate, only: my_task, master_task
      use ice_constants, only: c0, c1
      use ice_domain, only: nblocks
      use ice_domain_size, only: nilyr, nslyr, ncat, max_blocks
      use ice_fileunits, only: nu_diag, nu_rst_pointer, nu_dump
      use ice_flux, only: scale_factor, swvdr, swvdf, swidr, swidf, &
          strocnxT, strocnyT, sst, frzmlt, iceumask, coszen, &
          stressp_1, stressp_2, stressp_3, stressp_4, &
          stressm_1, stressm_2, stressm_3, stressm_4, &
          stress12_1, stress12_2, stress12_3, stress12_4
      use ice_ocean, only: oceanmixed_ice
      use ice_read_write, only: ice_open, ice_write
      use ice_state, only: aicen, vicen, vsnon, trcrn, &
          nt_Tsfc, nt_sice, nt_qice, nt_qsno, uvel, vvel

      character(len=char_len_long), intent(in), optional :: filename_spec

      ! local variables

      integer (kind=int_kind) :: &
          i, j, k, n, iblk, &     ! counting indices
          iyear, imonth, iday     ! year, month, day

      character(len=char_len_long) :: filename

      logical (kind=log_kind) :: diag

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         work1

      character (len=3) :: nchar

      if (present(filename_spec)) then
         call init_restart_write(filename_spec)
      else
         call init_restart_write
      endif

      diag = .true.

      !-----------------------------------------------------------------
      ! state variables
      ! Tsfc is the only tracer written to binary files.  All other
      ! tracers are written to their own dump/restart binary files.
      !-----------------------------------------------------------------

      call write_restart_field(nu_dump,0,aicen(:,:,:,:),'ruf8','aicen',ncat,diag)
      call write_restart_field(nu_dump,0,vicen(:,:,:,:),'ruf8','vicen',ncat,diag)
      call write_restart_field(nu_dump,0,vsnon(:,:,:,:),'ruf8','vsnon',ncat,diag)
      call write_restart_field(nu_dump,0,trcrn(:,:,nt_Tsfc,:,:),'ruf8','Tsfcn',ncat,diag)

      do k=1,nilyr
         write(nchar,'(i3.3)') k
         call write_restart_field(nu_dump,0,trcrn(:,:,nt_sice+k-1,:,:),'ruf8', &
                                 'sice'//trim(nchar),ncat,diag)
      enddo

      do k=1,nilyr
         write(nchar,'(i3.3)') k
         call write_restart_field(nu_dump,0,trcrn(:,:,nt_qice+k-1,:,:),'ruf8', &
                                 'qice'//trim(nchar),ncat,diag)
      enddo

      do k=1,nslyr
         write(nchar,'(i3.3)') k
         call write_restart_field(nu_dump,0,trcrn(:,:,nt_qsno+k-1,:,:),'ruf8', &
                                 'qsno'//trim(nchar),ncat,diag)
      enddo

      !-----------------------------------------------------------------
      ! velocity
      !-----------------------------------------------------------------
      call write_restart_field(nu_dump,0,uvel,'ruf8','uvel',1,diag)
      call write_restart_field(nu_dump,0,vvel,'ruf8','vvel',1,diag)

      !-----------------------------------------------------------------
      ! radiation fields
      !-----------------------------------------------------------------
#ifdef CESMCOUPLED
      call write_restart_field(nu_dump,0,coszen,'ruf8','coszen',1,diag)
#endif
      call write_restart_field(nu_dump,0,scale_factor,'ruf8','scale_factor',1,diag)

      call write_restart_field(nu_dump,0,swvdr,'ruf8','swvdr',1,diag)
      call write_restart_field(nu_dump,0,swvdf,'ruf8','swvdf',1,diag)
      call write_restart_field(nu_dump,0,swidr,'ruf8','swidr',1,diag)
      call write_restart_field(nu_dump,0,swidf,'ruf8','swidf',1,diag)

      !-----------------------------------------------------------------
      ! ocean stress (for bottom heat flux in thermo)
      !-----------------------------------------------------------------
      call write_restart_field(nu_dump,0,strocnxT,'ruf8','strocnxT',1,diag)
      call write_restart_field(nu_dump,0,strocnyT,'ruf8','strocnyT',1,diag)

      !-----------------------------------------------------------------
      ! internal stress
      !-----------------------------------------------------------------
      call write_restart_field(nu_dump,0,stressp_1,'ruf8','stressp_1',1,diag)
      call write_restart_field(nu_dump,0,stressp_3,'ruf8','stressp_3',1,diag)
      call write_restart_field(nu_dump,0,stressp_2,'ruf8','stressp_2',1,diag)
      call write_restart_field(nu_dump,0,stressp_4,'ruf8','stressp_4',1,diag)

      call write_restart_field(nu_dump,0,stressm_1,'ruf8','stressm_1',1,diag)
      call write_restart_field(nu_dump,0,stressm_3,'ruf8','stressm_3',1,diag)
      call write_restart_field(nu_dump,0,stressm_2,'ruf8','stressm_2',1,diag)
      call write_restart_field(nu_dump,0,stressm_4,'ruf8','stressm_4',1,diag)

      call write_restart_field(nu_dump,0,stress12_1,'ruf8','stress12_1',1,diag)
      call write_restart_field(nu_dump,0,stress12_3,'ruf8','stress12_3',1,diag)
      call write_restart_field(nu_dump,0,stress12_2,'ruf8','stress12_2',1,diag)
      call write_restart_field(nu_dump,0,stress12_4,'ruf8','stress12_4',1,diag)

      !-----------------------------------------------------------------
      ! ice mask for dynamics
      !-----------------------------------------------------------------
      
      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            work1(i,j,iblk) = c0
            if (iceumask(i,j,iblk)) work1(i,j,iblk) = c1
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      call write_restart_field(nu_dump,0,work1,'ruf8','iceumask',1,diag)

      ! for mixed layer model
      if (oceanmixed_ice) then
         call write_restart_field(nu_dump,0,sst,'ruf8','sst',1,diag)
         call write_restart_field(nu_dump,0,frzmlt,'ruf8','frzmlt',1,diag)
      endif

      end subroutine dumpfile

!=======================================================================

! Restarts from a dump
! author Elizabeth C. Hunke, LANL

      subroutine restartfile (ice_ic)

      use ice_boundary, only: ice_HaloUpdate_stress
      use ice_broadcast, only: broadcast_scalar
      use ice_blocks, only: nghost, nx_block, ny_block, block, get_block
      use ice_calendar, only: istep0, istep1, time, time_forc, calendar, npt, dt
      use ice_communicate, only: my_task, master_task
      use ice_constants, only: c0, p5, &
          field_loc_center, field_loc_NEcorner, &
          field_type_scalar, field_type_vector
      use ice_domain, only: nblocks, blocks_ice, distrb_info, halo_info
      use ice_domain_size, only: nilyr, nslyr, ncat, nx_global, ny_global, &
          max_ntrcr, max_blocks
      use ice_fileunits, only: nu_diag, nu_rst_pointer, nu_restart
      use ice_flux, only: scale_factor, swvdr, swvdf, swidr, swidf, &
          strocnxT, strocnyT, sst, frzmlt, iceumask, coszen, &
          faero_ocn,  fiso_ocn,  &
          stressp_1, stressp_2, stressp_3, stressp_4, &
          stressm_1, stressm_2, stressm_3, stressm_4, &
          stress12_1, stress12_2, stress12_3, stress12_4
      use ice_gather_scatter, only: scatter_global_stress
      use ice_grid, only: tmask, grid_type
      use ice_itd, only: aggregate, cleanup_itd
      use ice_ocean, only: oceanmixed_ice
      use ice_read_write, only: ice_open, ice_read, ice_read_global
      use ice_state, only: trcr_depend, aice, vice, vsno, trcr, &
          aice0, aicen, vicen, vsnon, trcrn, aice_init, &
          nt_Tsfc, nt_sice, nt_qice, nt_qsno, uvel, vvel, ntrcr, &
          nbtrcr
      use ice_zbgc_shared, only: first_ice
      use ice_therm_shared, only: heat_capacity

      character (*), optional :: ice_ic

      ! local variables

      integer (kind=int_kind) :: &
         i, j, k, n, iblk, &     ! counting indices
         ilo, ihi, jlo, jhi, & ! beginning and end of physical domain  
         iignore                 ! dummy variable

      type (block) :: &
         this_block      ! block information for current block        
 
      logical (kind=log_kind) :: &
         l_stop          ! if true, abort model                                   
 
      integer (kind=int_kind) :: &
         istop, jstop    ! indices of grid cell where model aborts     

      real (kind=real_kind) :: &
         rignore                 ! dummy variable

      character(len=char_len_long) :: &
         filename, filename0

      logical (kind=log_kind) :: &
         diag

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         work1

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         work_g1, work_g2

      character (len=3) :: nchar

      call init_restart_read(ice_ic)

      diag = .true.

      !-----------------------------------------------------------------
      ! state variables
      ! Tsfc is the only tracer read in this file.  All other
      ! tracers are in their own dump/restart files.
      !-----------------------------------------------------------------
      if (my_task == master_task) &
           write(nu_diag,*) ' min/max area, vol ice, vol snow, e ice, e snow, Tsfc'

      call read_restart_field(nu_restart,0,aicen,'ruf8', &
              'aicen',ncat,diag,field_loc_center, field_type_scalar)
      call read_restart_field(nu_restart,0,vicen,'ruf8', &
              'vicen',ncat,diag,field_loc_center, field_type_scalar)
      call read_restart_field(nu_restart,0,vsnon,'ruf8', &
              'vsnon',ncat,diag,field_loc_center, field_type_scalar)
      call read_restart_field(nu_restart,0,trcrn(:,:,nt_Tsfc,:,:),'ruf8', &
              'Tsfcn',ncat,diag,field_loc_center, field_type_scalar)

      if (my_task == master_task) &
         write(nu_diag,*) 'min/max sice for each layer'
      do k=1,nilyr
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart,0,trcrn(:,:,nt_sice+k-1,:,:),'ruf8', &
              'sice'//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo

      if (my_task == master_task) &
         write(nu_diag,*) 'min/max qice for each layer'
      do k=1,nilyr
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart,0,trcrn(:,:,nt_qice+k-1,:,:),'ruf8', &
              'qice'//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo

      if (my_task == master_task) &
         write(nu_diag,*) 'min/max qsno for each layer'
      do k=1,nslyr
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart,0,trcrn(:,:,nt_qsno+k-1,:,:),'ruf8', &
              'qsno'//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo

      !-----------------------------------------------------------------
      ! velocity
      !-----------------------------------------------------------------
      if (my_task == master_task) &
           write(nu_diag,*) 'min/max velocity components'

      call read_restart_field(nu_restart,0,uvel,'ruf8', &
           'uvel',1,diag,field_loc_NEcorner, field_type_vector)
      call read_restart_field(nu_restart,0,vvel,'ruf8', &
           'vvel',1,diag,field_loc_NEcorner, field_type_vector)

      !-----------------------------------------------------------------
      ! radiation fields
      !-----------------------------------------------------------------

      if (my_task == master_task) &
         write(nu_diag,*) 'radiation fields'

#ifdef CESMCOUPLED
      call read_restart_field(nu_restart,0,coszen,'ruf8', &
           'coszen',1,diag, field_loc_center, field_type_scalar)
#endif
      call read_restart_field(nu_restart,0,scale_factor,'ruf8', &
           'scale_factor',1,diag, field_loc_center, field_type_scalar)
      call read_restart_field(nu_restart,0,swvdr,'ruf8', &
           'swvdr',1,diag,field_loc_center, field_type_scalar)
      call read_restart_field(nu_restart,0,swvdf,'ruf8', &
           'swvdf',1,diag,field_loc_center, field_type_scalar)
      call read_restart_field(nu_restart,0,swidr,'ruf8', &
           'swidr',1,diag,field_loc_center, field_type_scalar)
      call read_restart_field(nu_restart,0,swidf,'ruf8', &
           'swidf',1,diag,field_loc_center, field_type_scalar)

      !-----------------------------------------------------------------
      ! ocean stress
      !-----------------------------------------------------------------
      if (my_task == master_task) &
           write(nu_diag,*) 'min/max ocean stress components'

      call read_restart_field(nu_restart,0,strocnxT,'ruf8', &
           'strocnxT',1,diag,field_loc_center, field_type_vector)
      call read_restart_field(nu_restart,0,strocnyT,'ruf8', &
           'strocnyT',1,diag,field_loc_center, field_type_vector)

      !-----------------------------------------------------------------
      ! internal stress
      ! The stress tensor must be read and scattered in pairs in order
      ! to properly match corner values across a tripole grid cut.
      !-----------------------------------------------------------------
      if (my_task == master_task) write(nu_diag,*) &
           'internal stress components'
      
      call read_restart_field(nu_restart,0,stressp_1,'ruf8', &
           'stressp_1',1,diag,field_loc_center,field_type_scalar) ! stressp_1
      call read_restart_field(nu_restart,0,stressp_3,'ruf8', &
           'stressp_3',1,diag,field_loc_center,field_type_scalar) ! stressp_3
      call read_restart_field(nu_restart,0,stressp_2,'ruf8', &
           'stressp_2',1,diag,field_loc_center,field_type_scalar) ! stressp_2
      call read_restart_field(nu_restart,0,stressp_4,'ruf8', &
           'stressp_4',1,diag,field_loc_center,field_type_scalar) ! stressp_4

      call read_restart_field(nu_restart,0,stressm_1,'ruf8', &
           'stressm_1',1,diag,field_loc_center,field_type_scalar) ! stressm_1
      call read_restart_field(nu_restart,0,stressm_3,'ruf8', &
           'stressm_3',1,diag,field_loc_center,field_type_scalar) ! stressm_3
      call read_restart_field(nu_restart,0,stressm_2,'ruf8', &
           'stressm_2',1,diag,field_loc_center,field_type_scalar) ! stressm_2
      call read_restart_field(nu_restart,0,stressm_4,'ruf8', &
           'stressm_4',1,diag,field_loc_center,field_type_scalar) ! stressm_4

      call read_restart_field(nu_restart,0,stress12_1,'ruf8', &
           'stress12_1',1,diag,field_loc_center,field_type_scalar) ! stress12_1
      call read_restart_field(nu_restart,0,stress12_3,'ruf8', &
           'stress12_3',1,diag,field_loc_center,field_type_scalar) ! stress12_1

      call read_restart_field(nu_restart,0,stress12_2,'ruf8', &
           'stress12_2',1,diag,field_loc_center,field_type_scalar) ! stress12_2
      call read_restart_field(nu_restart,0,stress12_4,'ruf8', &
           'stress12_4',1,diag,field_loc_center,field_type_scalar) ! stress12_4

      if (trim(grid_type) == 'tripole') then
         call ice_HaloUpdate_stress(stressp_1, stressp_3, halo_info, &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_3, stressp_1, halo_info, &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_2, stressp_4, halo_info, &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_4, stressp_2, halo_info, &
                                    field_loc_center,  field_type_scalar)

         call ice_HaloUpdate_stress(stressm_1, stressm_3, halo_info, &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_3, stressm_1, halo_info, &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_2, stressm_4, halo_info, &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_4, stressm_2, halo_info, &
                                    field_loc_center,  field_type_scalar)

         call ice_HaloUpdate_stress(stress12_1, stress12_3, halo_info, &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stress12_3, stress12_1, halo_info, &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stress12_2, stress12_4, halo_info, &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stress12_4, stress12_2, halo_info, &
                                    field_loc_center,  field_type_scalar)
      endif

      !-----------------------------------------------------------------
      ! ice mask for dynamics
      !-----------------------------------------------------------------
      if (my_task == master_task) &
           write(nu_diag,*) 'ice mask for dynamics'

      call read_restart_field(nu_restart,0,work1,'ruf8', &
           'iceumask',1,diag,field_loc_center, field_type_scalar)

      iceumask(:,:,:) = .false.
      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            if (work1(i,j,iblk) > p5) iceumask(i,j,iblk) = .true.
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      ! for mixed layer model
      if (oceanmixed_ice) then

         if (my_task == master_task) &
              write(nu_diag,*) 'min/max sst, frzmlt'

         call read_restart_field(nu_restart,0,sst,'ruf8', &
              'sst',1,diag,field_loc_center, field_type_scalar)
         call read_restart_field(nu_restart,0,frzmlt,'ruf8', &
              'frzmlt',1,diag,field_loc_center, field_type_scalar)
      endif

      !-----------------------------------------------------------------
      ! Ensure unused stress values in west and south ghost cells are 0
      !-----------------------------------------------------------------
      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks
         do j = 1, nghost
         do i = 1, nx_block
            stressp_1 (i,j,iblk) = c0
            stressp_2 (i,j,iblk) = c0
            stressp_3 (i,j,iblk) = c0
            stressp_4 (i,j,iblk) = c0
            stressm_1 (i,j,iblk) = c0
            stressm_2 (i,j,iblk) = c0
            stressm_3 (i,j,iblk) = c0
            stressm_4 (i,j,iblk) = c0
            stress12_1(i,j,iblk) = c0
            stress12_2(i,j,iblk) = c0
            stress12_3(i,j,iblk) = c0
            stress12_4(i,j,iblk) = c0
         enddo
         enddo
         do j = 1, ny_block
         do i = 1, nghost
            stressp_1 (i,j,iblk) = c0
            stressp_2 (i,j,iblk) = c0
            stressp_3 (i,j,iblk) = c0
            stressp_4 (i,j,iblk) = c0
            stressm_1 (i,j,iblk) = c0
            stressm_2 (i,j,iblk) = c0
            stressm_3 (i,j,iblk) = c0
            stressm_4 (i,j,iblk) = c0
            stress12_1(i,j,iblk) = c0
            stress12_2(i,j,iblk) = c0
            stress12_3(i,j,iblk) = c0
            stress12_4(i,j,iblk) = c0
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      !-----------------------------------------------------------------
      ! Ensure ice is binned in correct categories
      ! (should not be necessary unless restarting from a run with
      !  different category boundaries).
      !
      ! If called, this subroutine does not give exact restart.
      !-----------------------------------------------------------------

      l_stop = .false.

      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks
        this_block = get_block(blocks_ice(iblk),iblk)
        ilo = this_block%ilo
        ihi = this_block%ihi
        jlo = this_block%jlo
        jhi = this_block%jhi

        call cleanup_itd (nx_block,    ny_block,   &
                        ilo, ihi,    jlo, jhi,   &
                        dt,          ntrcr,      &
                        aicen   (:,:,:,iblk), trcrn (:,:,1:ntrcr,:,iblk),   &
                        vicen   (:,:,:,iblk), vsnon (:,:,  :,iblk), &
                        aice0   (:,:,  iblk), aice      (:,:,iblk), &
                        trcr_depend(1:ntrcr), &
                        tr_aero=.false.,      &  
                        tr_iso=.false.,        &
                        tr_pond_topo=.false., &        
                        heat_capacity=heat_capacity,  &
                        nbtrcr=nbtrcr, &
                        first_ice=first_ice(:,:,:,iblk), &
                        l_stop=l_stop,                   &
                        istop=istop,         jstop=jstop)

      enddo
      !$OMP END PARALLEL DO

      !-----------------------------------------------------------------
      ! compute aggregate ice state and open water area
      !-----------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1, nblocks

         call aggregate (nx_block, ny_block, &
                         aicen(:,:,:,iblk),  &
                         trcrn(:,:,:,:,iblk),&
                         vicen(:,:,:,iblk),  &
                         vsnon(:,:,:,iblk),  &
                         aice (:,:,  iblk),  &
                         trcr (:,:,:,iblk),  &
                         vice (:,:,  iblk),  &
                         vsno (:,:,  iblk),  &
                         aice0(:,:,  iblk),  &
                         tmask(:,:,  iblk),  &
                         max_ntrcr,          &
                         trcr_depend)

         aice_init(:,:,iblk) = aice(:,:,iblk)

      enddo
      !$OMP END PARALLEL DO

      ! if runid is bering then need to correct npt for istep0
      if (trim(runid) == 'bering') then
         npt = npt - istep0
      endif

      end subroutine restartfile

!=======================================================================

! Restarts from a CICE v4.1 (binary) dump
! author Elizabeth C. Hunke, LANL

      subroutine restartfile_v4 (ice_ic)

      use ice_broadcast, only: broadcast_scalar
      use ice_blocks, only: nghost, nx_block, ny_block
      use ice_calendar, only: istep0, istep1, time, time_forc, calendar, npt
      use ice_communicate, only: my_task, master_task
      use ice_constants, only: c0, p5, &
          field_loc_center, field_loc_NEcorner, &
          field_type_scalar, field_type_vector
      use ice_domain, only: nblocks, distrb_info
      use ice_domain_size, only: nilyr, nslyr, ncat, nx_global, ny_global, &
          max_ntrcr, max_blocks
      use ice_fileunits, only: nu_diag, nu_rst_pointer, nu_restart
      use ice_flux, only: scale_factor, swvdr, swvdf, swidr, swidf, &
          strocnxT, strocnyT, sst, frzmlt, iceumask, &
          stressp_1, stressp_2, stressp_3, stressp_4, &
          stressm_1, stressm_2, stressm_3, stressm_4, &
          stress12_1, stress12_2, stress12_3, stress12_4
      use ice_gather_scatter, only: scatter_global_stress
      use ice_grid, only: tmask
      use ice_itd, only: aggregate
      use ice_ocean, only: oceanmixed_ice
      use ice_read_write, only: ice_open, ice_read, ice_read_global
      use ice_state, only: trcr_depend, aice, vice, vsno, trcr, &
          aice0, aicen, vicen, vsnon,trcrn, aice_init, &
          nt_Tsfc, nt_sice, nt_qice, nt_qsno, uvel, vvel


      character (*), optional :: ice_ic

      ! local variables

      integer (kind=int_kind) :: &
         i, j, k, n, iblk, &     ! counting indices
         iignore                 ! dummy variable

      real (kind=real_kind) :: &
         rignore                 ! dummy variable

      character(len=char_len_long) :: &
         filename, filename0

      logical (kind=log_kind) :: &
         diag

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         work1

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         work_g1, work_g2

      if (present(ice_ic)) then
         filename = ice_ic
      elseif (my_task == master_task) then
         open(nu_rst_pointer,file=pointer_file)
         read(nu_rst_pointer,'(a)') filename0
         filename = trim(filename0)
         close(nu_rst_pointer)
         write(nu_diag,*) 'Read ',pointer_file(1:lenstr(pointer_file))
      endif

      call ice_open(nu_restart,filename,0)

      if (my_task == master_task) &
         write(nu_diag,*) 'Using restart dump=', trim(filename)

      if (use_restart_time) then

         if (my_task == master_task) then
            read (nu_restart) istep0,time,time_forc
            write(nu_diag,*) 'Restart read at istep=',istep0,time,time_forc
         endif
         call broadcast_scalar(istep0,master_task)
         istep1 = istep0
         call broadcast_scalar(time,master_task)
         call broadcast_scalar(time_forc,master_task)
         call calendar(time)

      else

         if (my_task == master_task) &
            read (nu_restart) iignore,rignore,rignore

      endif

      diag = .true.     ! write min/max diagnostics for field

      !-----------------------------------------------------------------
      ! state variables
      ! Tsfc is the only tracer read in this file.  All other
      ! tracers are in their own dump/restart files.
      !-----------------------------------------------------------------
      do n=1,ncat
         if (my_task == master_task) &
              write(nu_diag,*) 'cat ',n, &
                               ' min/max area, vol ice, vol snow, Tsfc'

         call ice_read(nu_restart,0,aicen(:,:,n,:),'ruf8',diag, &
                       field_loc_center, field_type_scalar)
         call ice_read(nu_restart,0,vicen(:,:,n,:),'ruf8',diag, &
                       field_loc_center, field_type_scalar)
         call ice_read(nu_restart,0,vsnon(:,:,n,:),'ruf8',diag, &
                       field_loc_center, field_type_scalar)
         call ice_read(nu_restart,0,trcrn(:,:,nt_Tsfc,n,:),'ruf8',diag, &
                       field_loc_center, field_type_scalar)

         if (my_task == master_task) &
              write(nu_diag,*) 'cat ',n, 'min/max sice for each layer'
         do k=1,nilyr
            call ice_read(nu_restart,0,trcrn(:,:,nt_sice+k-1,n,:),'ruf8',diag, &
                          field_loc_center, field_type_scalar)
         enddo

         if (my_task == master_task) &
              write(nu_diag,*) 'cat ',n, 'min/max qice for each layer'
         do k=1,nilyr
            call ice_read(nu_restart,0,trcrn(:,:,nt_qice+k-1,n,:),'ruf8',diag, &
                          field_loc_center, field_type_scalar)
         enddo

         if (my_task == master_task) &
              write(nu_diag,*) 'cat ',n, 'min/max qsno for each layer'
         do k=1,nslyr
            call ice_read(nu_restart,0,trcrn(:,:,nt_qsno+k-1,n,:),'ruf8',diag, &
                          field_loc_center, field_type_scalar)
         enddo
      enddo ! ncat

      !-----------------------------------------------------------------
      ! velocity
      !-----------------------------------------------------------------
      if (my_task == master_task) &
           write(nu_diag,*) 'min/max velocity components'

      call ice_read(nu_restart,0,uvel,'ruf8',diag, &
                       field_loc_NEcorner, field_type_vector)
      call ice_read(nu_restart,0,vvel,'ruf8',diag, &
                       field_loc_NEcorner, field_type_vector)

      !-----------------------------------------------------------------
      ! radiation fields
      !-----------------------------------------------------------------

      if (my_task == master_task) &
         write(nu_diag,*) 'radiation fields'

      call ice_read(nu_restart,0,scale_factor,'ruf8',diag, &
                    field_loc_center, field_type_scalar)
      call ice_read(nu_restart,0,swvdr,'ruf8',diag, &
                    field_loc_center, field_type_scalar)
      call ice_read(nu_restart,0,swvdf,'ruf8',diag, &
                    field_loc_center, field_type_scalar)
      call ice_read(nu_restart,0,swidr,'ruf8',diag, &
                    field_loc_center, field_type_scalar)
      call ice_read(nu_restart,0,swidf,'ruf8',diag, &
                    field_loc_center, field_type_scalar)

      !-----------------------------------------------------------------
      ! ocean stress
      !-----------------------------------------------------------------
      if (my_task == master_task) &
           write(nu_diag,*) 'min/max ocean stress components'

      call ice_read(nu_restart,0,strocnxT,'ruf8',diag, &
                       field_loc_center, field_type_vector)
      call ice_read(nu_restart,0,strocnyT,'ruf8',diag, &
                       field_loc_center, field_type_vector)

      !-----------------------------------------------------------------
      ! internal stress
      ! The stress tensor must be read and scattered in pairs in order
      ! to properly match corner values across a tripole grid cut.
      !-----------------------------------------------------------------
      if (my_task == master_task) write(nu_diag,*) &
           'internal stress components'
      
      allocate (work_g1(nx_global,ny_global), &
                work_g2(nx_global,ny_global))

      call ice_read_global(nu_restart,0,work_g1,'ruf8',diag) ! stressp_1
      call ice_read_global(nu_restart,0,work_g2,'ruf8',diag) ! stressp_3
      call scatter_global_stress(stressp_1, work_g1, work_g2, &
                                 master_task, distrb_info)
      call scatter_global_stress(stressp_3, work_g2, work_g1, &
                                 master_task, distrb_info)

      call ice_read_global(nu_restart,0,work_g1,'ruf8',diag) ! stressp_2
      call ice_read_global(nu_restart,0,work_g2,'ruf8',diag) ! stressp_4
      call scatter_global_stress(stressp_2, work_g1, work_g2, &
                                 master_task, distrb_info)
      call scatter_global_stress(stressp_4, work_g2, work_g1, &
                                 master_task, distrb_info)

      call ice_read_global(nu_restart,0,work_g1,'ruf8',diag) ! stressm_1
      call ice_read_global(nu_restart,0,work_g2,'ruf8',diag) ! stressm_3
      call scatter_global_stress(stressm_1, work_g1, work_g2, &
                                 master_task, distrb_info)
      call scatter_global_stress(stressm_3, work_g2, work_g1, &
                                 master_task, distrb_info)

      call ice_read_global(nu_restart,0,work_g1,'ruf8',diag) ! stressm_2
      call ice_read_global(nu_restart,0,work_g2,'ruf8',diag) ! stressm_4
      call scatter_global_stress(stressm_2, work_g1, work_g2, &
                                 master_task, distrb_info)
      call scatter_global_stress(stressm_4, work_g2, work_g1, &
                                 master_task, distrb_info)

      call ice_read_global(nu_restart,0,work_g1,'ruf8',diag) ! stress12_1
      call ice_read_global(nu_restart,0,work_g2,'ruf8',diag) ! stress12_3
      call scatter_global_stress(stress12_1, work_g1, work_g2, &
                                 master_task, distrb_info)
      call scatter_global_stress(stress12_3, work_g2, work_g1, &
                                 master_task, distrb_info)

      call ice_read_global(nu_restart,0,work_g1,'ruf8',diag) ! stress12_2
      call ice_read_global(nu_restart,0,work_g2,'ruf8',diag) ! stress12_4
      call scatter_global_stress(stress12_2, work_g1, work_g2, &
                                 master_task, distrb_info)
      call scatter_global_stress(stress12_4, work_g2, work_g1, &
                                 master_task, distrb_info)

      deallocate (work_g1, work_g2)

      !-----------------------------------------------------------------
      ! ice mask for dynamics
      !-----------------------------------------------------------------
      if (my_task == master_task) &
           write(nu_diag,*) 'ice mask for dynamics'

      call ice_read(nu_restart,0,work1,'ruf8',diag, &
                    field_loc_center, field_type_scalar)

      iceumask(:,:,:) = .false.
      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            if (work1(i,j,iblk) > p5) iceumask(i,j,iblk) = .true.
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      ! for mixed layer model
      if (oceanmixed_ice) then

         if (my_task == master_task) &
              write(nu_diag,*) 'min/max sst, frzmlt'

         call ice_read(nu_restart,0,sst,'ruf8',diag, &
                       field_loc_center, field_type_scalar)
         call ice_read(nu_restart,0,frzmlt,'ruf8',diag, &
                       field_loc_center, field_type_scalar)
      endif

      if (my_task == master_task) close(nu_restart)

      !-----------------------------------------------------------------
      ! Ensure unused stress values in west and south ghost cells are 0
      !-----------------------------------------------------------------
      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks
         do j = 1, nghost
         do i = 1, nx_block
            stressp_1 (i,j,iblk) = c0
            stressp_2 (i,j,iblk) = c0
            stressp_3 (i,j,iblk) = c0
            stressp_4 (i,j,iblk) = c0
            stressm_1 (i,j,iblk) = c0
            stressm_2 (i,j,iblk) = c0
            stressm_3 (i,j,iblk) = c0
            stressm_4 (i,j,iblk) = c0
            stress12_1(i,j,iblk) = c0
            stress12_2(i,j,iblk) = c0
            stress12_3(i,j,iblk) = c0
            stress12_4(i,j,iblk) = c0
         enddo
         enddo
         do j = 1, ny_block
         do i = 1, nghost
            stressp_1 (i,j,iblk) = c0
            stressp_2 (i,j,iblk) = c0
            stressp_3 (i,j,iblk) = c0
            stressp_4 (i,j,iblk) = c0
            stressm_1 (i,j,iblk) = c0
            stressm_2 (i,j,iblk) = c0
            stressm_3 (i,j,iblk) = c0
            stressm_4 (i,j,iblk) = c0
            stress12_1(i,j,iblk) = c0
            stress12_2(i,j,iblk) = c0
            stress12_3(i,j,iblk) = c0
            stress12_4(i,j,iblk) = c0
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      !-----------------------------------------------------------------
      ! Ensure ice is binned in correct categories
      ! (should not be necessary unless restarting from a run with
      !  different category boundaries).
      !
      ! If called, this subroutine does not give exact restart.
      !-----------------------------------------------------------------

!!!      call  cleanup_itd

      !-----------------------------------------------------------------
      ! compute aggregate ice state and open water area
      !-----------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks

         call aggregate (nx_block, ny_block, &
                         aicen(:,:,:,iblk),  &
                         trcrn(:,:,:,:,iblk),&
                         vicen(:,:,:,iblk),  &
                         vsnon(:,:,:,iblk),  &
                         aice (:,:,  iblk),  &
                         trcr (:,:,:,iblk),  &
                         vice (:,:,  iblk),  &
                         vsno (:,:,  iblk),  &
                         aice0(:,:,  iblk),  &
                         tmask(:,:,  iblk),  &
                         max_ntrcr,          &
                         trcr_depend)

         aice_init(:,:,iblk) = aice(:,:,iblk)

      enddo
      !$OMP END PARALLEL DO

      ! creates netcdf if restart_format = 'nc'
      filename = trim(restart_dir) // '/iced.converted'
      call dumpfile(filename) 
      call final_restart
      ! stop

      ! if runid is bering then need to correct npt for istep0
      if (trim(runid) == 'bering') then
         npt = npt - istep0
      endif

      end subroutine restartfile_v4

!=======================================================================

      end module ice_restart_driver

!=======================================================================
