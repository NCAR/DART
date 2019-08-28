module fesom_dart_diag

  use g_config,   only : iniday, endday, runyear, runid, day2ext, &
                         save_count, resultpath, level_number,    &
                         ensid, ensmem, dart_days, dart_secs
  use o_param,    only : rad
  use o_mesh,     only : coord_nod3d, coord_nod2d, nod3d_below_nod2d
  use o_array,    only : ssh, uf, tracer
  use g_PARFE,    only : mydim_nod2d, myDim_nod3d, mype, myDim_elem2D
  use utilities,  only : r8, i4

  implicit none

  public ::  read_increment,            &
             read_ensemble_from_netcdf, &
             calc_ensemble_mean,        &
             calc_ensemble_variance,    &
             read_section_from_ino,     &
             read_section_from_inc,     &
             dart_obs_seq_proc

  integer(i4)                 :: i, j, k
  character(len=*), parameter :: filename="inputfile.cdf"

  contains

! read innovation file generated after postproc
    subroutine read_innovation

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
      filename=trim(ResultPath)//'FILTER/Innovation.cdf'
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

    end subroutine read_innovation


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
      filename=trim(ResultPath)//'FILTER/Increment.cdf'
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


! extract a cross-section from innovation
    subroutine read_section_from_ino

      character(20)     :: DAYNUM,LEVNUM
      character(8)      :: OUTDIR
      character(100)    :: OUTFILENAM
      integer           :: layer, nlayer
      real(r8)          :: tracer1(myDim_nod2D),tracer2(myDim_nod2D)

     OUTDIR='.'
     day2ext=1
       call read_innovation
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
    end subroutine read_section_from_ino


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

    subroutine dart_obs_seq_proc

      use utilities, only : skip_read_line

      integer, parameter        ::   iread=101
      integer, parameter        ::   iwrite1=103
      integer, parameter        ::   iwrite2=104
      character(13)             ::   obsfilename="obs_seq.final"
      character(13)             ::   soutfilename="obs_seq.salt"
      character(13)             ::   toutfilename="obs_seq.temp"
      character(100)            ::   header_obs_seq
      character(100)            ::   dummy_text

      integer                   ::   tot_kind, num_copy, num_qc
      integer                   ::   num_obs, tot_obs, obs_num
      integer                   ::   first_obs, last_obs
      integer                   ::   prev_obs, next_obs
      integer                   ::   flag1, flag2
      integer                   ::   sec, day
      integer,allocatable       ::   num_kind(:)
      integer                   ::   obs_kind, iwrite
      character(20),allocatable ::   nam_kind(:)

      real(r8),allocatable      ::   copy(:)
      real(r8)                  ::   lat_rad, lon_rad, depth
      real(r8)                  ::   obs_err

      open(unit=iwrite1,file=soutfilename,status='replace',access='append',form='formatted')
      open(unit=iwrite2,file=toutfilename,status='replace',access='append',form='formatted')
      open(unit=iread, file=obsfilename)
      call skip_read_line(iread,2)
      read(iread,*) tot_kind
      allocate(num_kind(tot_kind), nam_kind(tot_kind))
      do i=1,2
        read(iread,*) num_kind(i), nam_kind(i)
      end do
      read(iread,*) dummy_text,num_copy,dummy_text,num_qc
      allocate(copy(num_copy+num_qc))
      read(iread,*) dummy_text,num_obs,dummy_text,tot_obs


      call skip_read_line(iread,num_copy+num_qc)
      read(iread,*) dummy_text,first_obs,dummy_text,last_obs
LOOP_OBS: do j=1,tot_obs
        read(iread,*) dummy_text, obs_num
LOOP_COPY:do i=1,num_copy+num_qc
          read(iread,*) copy(i)
        end do LOOP_COPY

        read(iread,*) prev_obs, next_obs, flag1
        call skip_read_line(iread,2)
        read(iread,*) lat_rad, lon_rad, depth, flag2
        call skip_read_line(iread,1)
        read(iread,*) obs_kind
        read(iread,*) sec, day
        read(iread,*) obs_err
        if (obs_kind.eq.1) iwrite=iwrite1
        if (obs_kind.eq.2) iwrite=iwrite2

        write(iwrite,'(I5,3F10.5,2I5,8F10.5)') &
                obs_kind, lat_rad, lon_rad, depth, sec, day, &
                copy(1), copy(2), copy(3), copy(4), copy(5), &
                copy(num_copy+1), copy(num_copy+2), obs_err
      end do LOOP_OBS

      deallocate(num_kind, nam_kind, copy)
      close(iread)
      close(iwrite2)
      close(iwrite1)

    end subroutine dart_obs_seq_proc

end module fesom_dart_diag
