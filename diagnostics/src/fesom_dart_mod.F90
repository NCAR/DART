!
! Includes tools to process ensemble and dart outputs
! Tools are used to produce the analysis in:
! Aydoğdu, A., Hoar, T. J., Vukicevic, T., Anderson, J. L., Pinardi, N., Karspeck, A., Hendricks, J., Collins, N., Macchia, F., and
! Özsoy, E.: OSSE for a sustainable marine observing network in the Sea of Marmara, Nonlin. Processes Geophys., 25, 537-551,
! https://doi.org/10.5194/npg-25-537-2018, 2018
!
! Provided by ali.aydogdu@cmcc.it


module fesom_dart_mod

  use g_config,   only : iniday, endday, runyear, runid, day2ext, &
                         save_count, resultpath, level_number,    &
                         ensid, ensmem, dart_days, dart_secs
  use o_param,    only : rad
  use o_mesh,     only : coord_nod3d, coord_nod2d, nod3d_below_nod2d
  use o_array,    only : ssh, uf, tracer
  use g_PARFE,    only : mydim_nod2d, myDim_nod3d, mype, myDim_elem2D
  use utilities,  only : r8, i4

  implicit none

  public ::  read_increment,            & ! reads increments from the diff increment.nc of output_mean.nc and preassim_mean.nc
             read_ensemble_from_netcdf, & ! reads all ensemble members from fesom outputs to compute mean and variance
             calc_ensemble_mean,        & ! calculates ensemble mean
             calc_ensemble_variance,    & ! calculates ensemble variance
             read_section_from_NR_diff, & ! calculates diff between nature run and prior mean if a NR is provided
             read_section_from_inc        ! reads a horizontal section from increment.nc

  integer(i4)                 :: i, j, k
  character(len=*), parameter :: filename="inputfile.nc"

  contains

! read NR_diff file generated after postproc
    subroutine read_NR_diff

      #include "netcdf.inc"

      integer                   :: status, ncid, dimid_rec, nrec
      integer                   :: ssh_varid, tra_varid(2)
      integer                   :: istart(2), icount(2), n3
      character(100)            :: filename
      character(1)              :: trind
      real(r8), allocatable     :: aux3(:)

      n3   = myDim_nod3D
      nrec = day2ext
      allocate(aux3(myDim_nod3D))

      ! open files
      filename=trim(ResultPath)//'FILTER/NR_prior_diff.nc'
      print*, filename
      status = nf_open(filename, nf_nowrite, ncid)
      if (status .ne. nf_noerr) call handle_err(status)
      status = nf_inq_varid(ncid, 'temp', tra_varid(1))
      if (status .ne. nf_noerr) call handle_err(status)
      status = nf_inq_varid(ncid, 'salt', tra_varid(2))
      if (status .ne. nf_noerr) call handle_err(status)

      ! 3d fields
      istart=(/1,nrec/)
      icount=(/myDim_nod3D,1/)

      do j=1,2
         status=nf_get_vara_double(ncid, tra_varid(j), istart, icount, aux3)
         if (status .ne. nf_noerr) call handle_err(status)
         tracer(:,j)=aux3(1:myDim_nod3D)
      end do

      status=nf_close(ncid)
      if (status .ne. nf_noerr) call handle_err(status)

      ! the next record to be saved
      save_count=nrec+1

      deallocate(aux3)

    end subroutine read_NR_diff


! read increment file generated after postproc
    subroutine read_increment

    #include "netcdf.inc"

      integer                   :: status, ncid, dimid_rec, nrec
      integer                   :: ssh_varid, tra_varid(2)
      integer                   :: copy
      integer                   :: istart(3), icount(3), n3
      character(100)            :: filename
      character(1)              :: trind
      real(r8), allocatable ::  aux3(:,:)

      copy = 1
      n3   = myDim_nod3D
      nrec = day2ext
      allocate(aux3(myDim_nod3D,copy))

      ! open files
      filename=trim(ResultPath)//'FILTER/Increment.nc'
      print*, filename
      status = nf_open(filename, nf_nowrite, ncid)
      if (status .ne. nf_noerr) call handle_err(status)
      status = nf_inq_varid(ncid, 'temp', tra_varid(1))
      if (status .ne. nf_noerr) call handle_err(status)
      status = nf_inq_varid(ncid, 'salt', tra_varid(2))
      if (status .ne. nf_noerr) call handle_err(status)

      ! 3d fields
      istart=(/1,1,nrec/)
      icount=(/myDim_nod3D,copy,1/)

      do j=1,2
         status=nf_get_vara_double(ncid, tra_varid(j), istart, icount, aux3)
         if (status .ne. nf_noerr) call handle_err(status)
         tracer(:,j)=aux3(1:myDim_nod3D,1)
      end do

      status=nf_close(ncid)
      if (status .ne. nf_noerr) call handle_err(status)

      ! the next record to be saved
      save_count=nrec+1

      deallocate(aux3)

    end subroutine read_increment

! compute ensemble mean and spread from fesom outputs
    subroutine read_ensemble_from_netcdf

      character(6)       :: DAYNUM,LEVNUM
      character(7)       :: OUTDIR
      character(100)     :: OUTFILENAM,MaiNPaTH,ReSuPaTH
      integer            :: layer, realday

      OUTDIR='ENSMEAN'
      layer=level_number
      write(LEVNUM,'(a,i3.3)')'LEV',layer
      write(MaiNPaTH,'(a)')ResultPath

DAYLOOP:do day2ext=iniday,endday
        realday=(day2ext-1)/4+1
ENSLOOP:do j=1,ensmem
          write(runid,'(a,i2.2)')'ENS',j
          write(ReSuPaTH,'(a)')ResultPath
          write(ResultPath,'(2a,i2.2,a)')trim(ReSuPaTH),'ENS',j,'/'
          write(DAYNUM,'(a,i3.3)')'DAY',realday
          print*,DAYNUM
          call oce_input
          OUTFILENAM=trim(OUTDIR)//'/'//trim(OUTDIR)//'_'//runid//'_'//runyear//'_'//DAYNUM//'_'//LEVNUM//'.asc'
          open(unit=101,file=OUTFILENAM,status='replace',access='append',form='formatted')

NODLOOP:  do i=1,myDim_nod2D
            if ( nod3D_below_nod2D(layer,i).ge.1.and.nod3D_below_nod2D(layer,i).le.myDim_nod3D ) then
              write(101,'(7F15.9,1X,I7.7)') &
              coord_nod2D(1,i)/rad,coord_nod2D(2,i)/rad, &
              tracer((nod3D_below_nod2D(layer,i)),1), &
              tracer((nod3D_below_nod2D(layer,i)),2), &
              uf(nod3d_below_nod2d(layer,i)),    &
              uf(nod3d_below_nod2d(layer,i)+myDim_nod3D), ssh(i), &
              nod3D_below_nod2D(layer,i)
            else
              write(101,'(7F15.9,1X,a)') &
              coord_nod2D(1,i)/rad,coord_nod2D(2,i)/rad, &
              0,0,0,0,0,'-999'
            end if
          end do NODLOOP

          write(ResultPath,'(a)')MaiNPaTH
       end do ENSLOOP

      call calc_ensemble_mean(DAYNUM,LEVNUM)
      call calc_ensemble_variance(DAYNUM,LEVNUM)
      end do DAYLOOP
    end subroutine read_ensemble_from_netcdf

! calculate ensemble mean
    subroutine calc_ensemble_mean(DAYNUM,LEVNUM)

      character(6)         :: DAYNUM,LEVNUM
      character(7)         :: IDIR='ENSMEAN'
      character(7)         :: ODIR='ENSMEAN'
      character(100)       :: IFILENAM,OFILENAM
      real(r8),allocatable :: whole_ascii_data(:,:),ensem_ascii_data(:,:)

      allocate(ensem_ascii_data(8,myDim_nod2d))
      OFILENAM=trim(ODIR)//'/'//trim(ODIR)//'_'//ensid//'_EMEAN_'//runyear//'_'//DAYNUM//'_'//LEVNUM//'.asc'
      ensem_ascii_data=0
      do j=1,ensmem
      allocate(whole_ascii_data(8,myDim_nod2d))
        write(runid,'(a,i2.2)')'ENS',j
        IFILENAM=trim(IDIR)//'/'//trim(IDIR)//'_'//runid//'_'//runyear//'_'//DAYNUM//'_'//LEVNUM//'.asc'
        open(unit=102,file=IFILENAM,status='old')
        read(102,*) whole_ascii_data
        ensem_ascii_data=ensem_ascii_data+whole_ascii_data
        close(102)
      deallocate(whole_ascii_data)
      end do
      open(unit=103,file=OFILENAM,status='replace',access='append',form='formatted')
      write(103,'(8F20.9)') ensem_ascii_data/ensmem
      close(103)
      deallocate(ensem_ascii_data)

    end subroutine calc_ensemble_mean

! calculate variance
    subroutine calc_ensemble_variance(DAYNUM,LEVNUM)

      character(6)         :: DAYNUM,LEVNUM
      character(7)         :: IDIR='ENSMEAN'
      character(7)         :: MDIR='ENSMEAN'
      character(7)         :: ODIR='ENSMEAN'
      character(100)       :: IFILENAM,OFILENAM,MFILENAM
      real(r8),allocatable :: whole_ascii_data(:,:),ensem_ascii_data(:,:) , &
                            meanm_ascii_data(:,:)

      allocate(ensem_ascii_data(8,myDim_nod2d))
      allocate(meanm_ascii_data(8,myDim_nod2d))
      MFILENAM=trim(ODIR)//'/'//trim(ODIR)//'_'//ensid//'_EMEAN_'//runyear//'_'//DAYNUM//'_'//LEVNUM//'.asc'
      open(unit=101,file=MFILENAM,status='old')
        read(101,*) meanm_ascii_data
      OFILENAM=trim(ODIR)//'/'//trim(ODIR)//'_'//ensid//'_STDEV_'//runyear//'_'//DAYNUM//'_'//LEVNUM//'.asc'
      ensem_ascii_data=0
      do j=1,ensmem
      allocate(whole_ascii_data(8,myDim_nod2d))
        write(runid,'(a,i2.2)')'ENS',j
        IFILENAM=trim(IDIR)//'/'//trim(IDIR)//'_'//runid//'_'//runyear//'_'//DAYNUM//'_'//LEVNUM//'.asc'
        open(unit=102,file=IFILENAM,status='old')
        read(102,*) whole_ascii_data
        whole_ascii_data=(whole_ascii_data-meanm_ascii_data)**2
        ensem_ascii_data=ensem_ascii_data+whole_ascii_data
        close(102)
      deallocate(whole_ascii_data)
      end do
      ensem_ascii_data=ensem_ascii_data/ensmem
      open(unit=103,file=OFILENAM,status='replace',access='append',form='formatted')
      do i=1,myDim_nod2d
        write(103,'(8F20.9)') meanm_ascii_data(1:2,i),ensem_ascii_data(3:7,i)/ensmem, &
                             meanm_ascii_data(8,i)
      end do
      close(103)
      deallocate(meanm_ascii_data)
      deallocate(ensem_ascii_data)

    end subroutine calc_ensemble_variance


! extract a cross-section from NR_diff
    subroutine read_section_from_NR_diff

      character(20)     :: DAYNUM,LEVNUM
      character(8)      :: OUTDIR
      character(100)    :: OUTFILENAM
      integer           :: layer, nlayer
      real(r8)          :: tracer1(myDim_nod2D),tracer2(myDim_nod2D)

     OUTDIR='.'
     day2ext=1
       call read_NR_diff
       write(DAYNUM,'(i6.6,a,i5.5)')dart_days,'_',dart_secs
       print*,DAYNUM
       OUTFILENAM=trim(OUTDIR)//'/INO_'//runid//'_'//runyear//'_'//trim(DAYNUM)//'.asc'
       open(unit=101,file=trim(OUTFILENAM),status='replace',access='append',form='formatted')
       do i=1,myDim_nod2D
        tracer1 = 0
        tracer2 = 0
        nlayer  = 0
       do layer=1,10
         if ( nod3D_below_nod2D(layer,i).ge.1.and.nod3D_below_nod2D(layer,i).le.myDim_nod3D ) then
                tracer1(i)=tracer1(i)+(tracer((nod3D_below_nod2D(layer,i)),1)**2)
                tracer2(i)=tracer2(i)+(tracer((nod3D_below_nod2D(layer,i)),2)**2)
                nlayer=nlayer+1
                else
                nlayer=nlayer
         end if
       end do
       tracer1(i)=sqrt(tracer1(i)/nlayer)
       tracer2(i)=sqrt(tracer2(i)/nlayer)
       write(101,'(4F15.9,1X,I7.7)') &
           coord_nod2D(1,i)/rad,coord_nod2D(2,i)/rad, &
           tracer1(i), tracer2(i)
       end do
       close(101)
    end subroutine read_section_from_NR_diff


! extract a cross-section from increment
    subroutine read_section_from_inc

      character(8)       :: OUTDIR
      character(20)      :: DAYNUM,LEVNUM
      character(100)     :: OUTFILENAM
      integer            :: layer

     OUTDIR='.'
     do day2ext=iniday,endday
       call read_increment
       write(DAYNUM,'(i6.6,a,i5.5)')dart_days,'_',dart_secs
       print*,DAYNUM
       do j=1,3
       layer=level_number
       write(LEVNUM,'(a,i3.3)')'LEV',layer
       print*,LEVNUM
       OUTFILENAM=trim(OUTDIR)//'/INC_'//runid//'_'//runyear//'_'//trim(DAYNUM)//'_'//trim(LEVNUM)//'.asc'
       open(unit=101,file=trim(OUTFILENAM),status='replace',access='append',form='formatted')
       do i=1,myDim_nod2D
         if ( nod3D_below_nod2D(layer,i).ge.1.and.nod3D_below_nod2D(layer,i).le.myDim_nod3D ) then
           write(101,'(4F15.9,1X,I7.7)') &
           coord_nod2D(1,i)/rad,coord_nod2D(2,i)/rad, &
           tracer((nod3D_below_nod2D(layer,i)),1), &
           tracer((nod3D_below_nod2D(layer,i)),2), &
           nod3D_below_nod2D(layer,i)
         else
           write(101,'(7F15.9,1X,I7.7)') &
           coord_nod2D(1,i)/rad,coord_nod2D(2,i)/rad, &
           0,0,-999
         end if
       end do
       close(101)
       end do
     end do
    end subroutine read_section_from_inc


end module fesom_dart_mod
