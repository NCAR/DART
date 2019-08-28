!
! Includes tools to process forcing files
! Tools are used to produce the analysis in:
! Aydoğdu, A., Pinardi, N., Özsoy, E., Danabasoglu, G., Gürses, Ö., and Karspeck, A.: Circulation of the Turkish
! Straits System under interannual atmospheric forcing, Ocean Sci., 14, 999-1019, https://doi.org/10.5194/os-14-999-2018, 2018b.
!
! Provided by ali.aydogdu@cmcc.it
!
module fesom_ocean_mod

use g_config,   only : runyear, day2ext, resultpath, runid, save_count,      &
                       level_number, thalweg_directory, flux_section,        &
                       iniday, endday
use o_param,    only : km3yr2m3sec, rad, g, pi,                              &
                       marmax_lon, marmin_lon, marmax_lat, marmin_lat
use o_elements, only : voltriangle, voltetra, cluster_area_2d,               &
                       elem2d_nodes, elem3d_nodes, bafux_3d, bafuy_3d,       &
                       marm_area, bosp_area, dard_area, ocean_area
use o_mesh,     only : nod3D_below_nod2D, coord_nod2d, coord_nod3d,          &
                       max_num_layers, layerdepth, num_layers_below_nod2D
use o_array,    only : ssh, tracer, uf, density_insitu
!use g_clock
use g_parfe,    only : mydim_elem2d, mydim_elem3d, mydim_nod2d, mydim_nod3d, &
                       mype
use utilities,  only : r8, i4

implicit none

public   :: find_surface_area,          &   ! computes surface area in a polygon (should be provided)
            basin_mean_evolution,       &   ! computes the mean of variables in the whole domain
            marmara_mean_evolution,     &   ! computes the mean of variables in a region (marmax_lon, marmin_lon, marmax_lat, marmin_lat)
            read_section_from_netcdf,   &   ! reads and extracts a level from fesom ocean outputs
            read_thalweg_from_nc,       &   ! reads and extracts a transect from fesom ocean outputs (file to be provided)
            calc_section_annual_mean,   &   ! computes yearly average of a level
            calc_section_monthly_mean,  &   ! computes monthly average of a level
            calc_thalweg_annual_mean,   &   ! computes yearly average of a transect
            calc_thalweg_monthly_mean,  &   ! computes monthly average of a transect
            velocity_at_the_exit,       &   ! computes velocity at strait exits (a transect should be provided)
            dardanelles_for_MFS,        &   ! extracts a region in the exit of the Dardanelles (a region to be provided)
            bosphorus_for_blk_mfs,      &   ! extracts a region in the exit of the Bosphorus (a region to be provided)
            surface_kinetic_energy,     &   ! computes surface kinetic energy for 2D vis.
            total_kinetic_energy,       &   ! computes total kinetic energy for 2D vis
            compute_vorticity               ! computes vorticity for 2D vis

integer  :: i, j, k

contains

subroutine find_surface_area

  marm_area = surface_area('BOUND_MARM')
  print*, 'Marmara Sea Area: ', marm_area
!  bosp_area = surface_area('BOUND_BOSP')
!  print*, 'Bosphorus Strait Area: ', bosp_area
!  dard_area = surface_area('BOUND_DARD')
!  print*, 'Dardanelles Strait Area: ', dard_area

end subroutine find_surface_area

function surface_area(BONDFILE)

  use utilities, only : inside

  logical           :: check
  integer           :: elem, elnodes(3)
  integer           :: nline, inc
  real(r8)          :: x0, y0, vol, surface_area
  real(r8),parameter:: inv3=1.0/3.0_8
  character(100)    :: POLYFILE
  character(10)     :: BONDFILE
  real(r8), allocatable, dimension(:,:) :: polygon

!     allocate(cluster_area_2d(myDim_nod2D))
  cluster_area_2d=0.0
  ocean_area = 0.0

  POLYFILE=trim(thalweg_directory)//'/'//trim(BONDFILE)//'.lst'
  print*, 'polygon file: ',POLYFILE
  open(unit=111,file=POLYFILE)
  read(111,*) nline

  allocate(polygon(2,nline))
  do inc = 1,nline
  read(111,*) polygon(1,inc), polygon(2,inc)
  end do
  close(111)

  do elem=1,myDim_elem2d
     elnodes=elem2d_nodes(:,elem)
     vol=voltriangle(elem)*inv3
     cluster_area_2d(elnodes)=cluster_area_2d(elnodes)+vol
  end do

  vol=0.0
  do i=1,myDim_nod2d
    x0=(coord_nod2D(1,i)/rad)
    y0=(coord_nod2D(2,i)/rad)
    check = inside(x0, y0, polygon(1,:), polygon(2,:), nline)
    if ( check ) then
      vol=vol+cluster_area_2d(i)
    endif
  end do
  surface_area = vol
  return
end function surface_area

subroutine basin_mean_evolution

  integer           :: ii,n3
  real(r8)          :: inv3=1.0/3.0_8
  real(r8)          :: inv4=1.0/4.0_8
  real(r8)          :: inv12=1.0/12.0_8
  real(r8)          :: sshsrf, sshelem, totvol, toarea, nodarea
  real(r8)          :: sstsrf, temelem, temvol, sstelem
  real(r8)          :: ssssrf, salelem, salvol, ssselem
  real(r8)          :: wsrf, welem
  character(20)     :: outdir
  character(46)     :: filenam
  real(kind=8), allocatable :: aux3(:)


  n3=myDim_nod3D
  outdir='.'
  filenam=trim(outdir)//'/MEANEVOLUTION'//runid//runyear//'3D.asc'

  open(unit=101,file=filenam,status='replace',access='append',form='formatted')

DAYLOOP:  do j=iniday,endday
    allocate(aux3(myDim_nod3D))
    print*, 'day: ',j, runid
    day2ext=j
    call oce_input
    aux3(1:myDim_nod3D)=uf(1+2*n3:3*n3)
    wsrf=0
    sshsrf=0
    sstsrf=0
    ssssrf=0
    wsrf=0
    toarea=0
    nodarea=0
N2LOOP:      do i=1,myDim_elem2D
       toarea=toarea+voltriangle(i)
       sshelem=0
       sstelem=0
       ssselem=0
       welem=0
E2LOOP:          do k=1,3
               sshelem=sshelem+(ssh(elem2D_nodes(k,i)))*voltriangle(i)
               sstelem=sstelem+(tracer(elem2D_nodes(k,i),1)*voltriangle(i))
               ssselem=ssselem+(tracer(elem2D_nodes(k,i),2)*voltriangle(i))
               welem=welem+(aux3(elem2D_nodes(k,i))*voltriangle(i))
          end do E2LOOP
      sshsrf=sshsrf+sshelem*inv3
      sstsrf=sstsrf+sstelem*inv3
      ssssrf=ssssrf+ssselem*inv3
      wsrf=wsrf+welem*inv3
      end do N2LOOP
    salvol=0
    temvol=0
    totvol=0
       call update_mesh
N3LOOP:      do i=1,myDim_elem3D
       totvol=totvol+voltetra(i)
       temelem=0
       salelem=0
E3LOOP:          do k=1,4
          temelem=temelem+((tracer(elem3D_nodes(k,i),1)*voltetra(i)))
          salelem=salelem+((tracer(elem3D_nodes(k,i),2)*voltetra(i)))
          end do E3LOOP
          salvol=salvol+salelem*inv4
          temvol=temvol+temelem*inv4
      end do N3LOOP
     write(101,10) sshsrf/toarea, sstsrf/toarea, ssssrf/toarea, &
                   temvol/totvol, salvol/totvol, wsrf/toarea
    10   FORMAT(1X, F17.9 ,1X, F17.9, 1X, F17.9, 1X, F17.9, 1X, F17.9, 1X, F17.9)
    deallocate(aux3)
end do DAYLOOP

close(101)
end subroutine basin_mean_evolution

subroutine marmara_mean_evolution

  real(r8)            :: inv3=1.0/3.0_8
  real(r8)            :: inv4=1.0/4.0_8
  real(r8)            :: inv12=1.0/12.0_8
  integer             :: elnodes2(3),elnodes3(4)
  real(r8)            :: sshsrf, sshelem, totvol, toarea, nodarea, sshnode
  real(r8)            :: sstsrf, temelem, sstvol, sstelem
  real(r8)            :: ssssrf, salelem, sssvol, ssselem
  character(20)       :: outdir
  character(46)       :: filenam

  outdir='.'
  filenam=trim(outdir)//'/MARMEVOLUTION'//runid//runyear//'3D.asc'

  open(unit=101,file=filenam,status='replace',access='append',form='formatted')

  do j=iniday,endday
    print*, 'day: ',j, runid
    day2ext=j
    call oce_input
    sshsrf=0
    sstsrf=0
    ssssrf=0
    toarea=0
    nodarea=0
    sshnode=0
      do i=1,myDim_elem2D
       elnodes2=elem2D_nodes(:,i)
       do k=1,3
         if ( (coord_nod2D(2,elnodes2(k))/rad).gt.marmax_lat .or. (coord_nod2D(1,elnodes2(k))/rad).gt.marmax_lon &
         .or. (coord_nod2D(2,elnodes2(k))/rad).lt.marmin_lat .or. (coord_nod2D(1,elnodes2(k))/rad).lt.marmin_lon ) then
            go to 1232
         endif
       enddo

       toarea=toarea+voltriangle(i)
       sshelem=0
       sstelem=0
       ssselem=0
          do k=1,3
               sshelem=sshelem+(ssh(elem2D_nodes(k,i)))*voltriangle(i)
               sstelem=sstelem+(tracer(elem2D_nodes(k,i),1)*voltriangle(i))
               ssselem=ssselem+(tracer(elem2D_nodes(k,i),2)*voltriangle(i))
          end do
      sshsrf=sshsrf+sshelem*inv3
      sstsrf=sstsrf+sstelem*inv3
      ssssrf=ssssrf+ssselem*inv3
        1232 continue
      end do
    sssvol=0
    sstvol=0
    totvol=0
      do i=1,myDim_elem3D
       elnodes3=elem3D_nodes(:,i)
       do k=1,4
         if ( (coord_nod2D(2,elnodes3(k))/rad).gt.marmax_lat .or. (coord_nod2D(1,elnodes3(k))/rad).gt.marmax_lon &
         .or. (coord_nod2D(2,elnodes3(k))/rad).lt.marmin_lat .or. (coord_nod2D(1,elnodes3(k))/rad).lt.marmin_lon ) then
            go to 1233
         endif
       enddo
       totvol=totvol+voltetra(i)
       temelem=0
       salelem=0
          do k=1,4
          temelem=temelem+((tracer(elem3D_nodes(k,i),1)*voltetra(i)))
          salelem=salelem+((tracer(elem3D_nodes(k,i),2)*voltetra(i)))
          end do
          sssvol=sssvol+salelem*inv4
          sstvol=sstvol+temelem*inv4
        1233 continue
      end do
     write(101,10) sshsrf/toarea, sstsrf/toarea, ssssrf/toarea, &
                   sstvol/totvol, sssvol/totvol
    10   FORMAT(1X, F17.9 ,1X, F17.9, 1X, F17.9, 1X, F17.9, 1X, F17.9)
  end do

  close(101)
end subroutine marmara_mean_evolution

subroutine read_section_from_netcdf

  character(6)       :: DAYNUM,LEVNUM
  character(7)       :: OUTDIR
  character(100)     :: OUTFILENAM
  integer            :: layer
  OUTDIR='SECTION'
  layer=level_number
  write(LEVNUM,'(a,i3.3)')'LEV',layer
  do day2ext=iniday,endday
    write(DAYNUM,'(a,i3.3)')'DAY',day2ext
    print*,DAYNUM
    call oce_input
    OUTFILENAM=trim(OUTDIR)//'/'//trim(OUTDIR)//'_'//runid//'_'//runyear//'_'//DAYNUM//'_'//LEVNUM//'.asc'
    open(unit=101,file=OUTFILENAM,status='replace',access='append',form='formatted')
    do i=1,myDim_nod2D
      if ( nod3D_below_nod2D(layer,i).ge.1.and.nod3D_below_nod2D(layer,i).le.myDim_nod3D ) then
        write(101,'(7F15.9,1X,I7.7)')  &
        coord_nod2D(1,i)/rad,coord_nod2D(2,i)/rad, &
        tracer((nod3D_below_nod2D(layer,i)),1), &
        tracer((nod3D_below_nod2D(layer,i)),2), &
        uf(nod3d_below_nod2d(layer,i)),    &
        uf(nod3d_below_nod2d(layer,i)+myDim_nod3D), ssh(i), &
        nod3D_below_nod2D(layer,i)
      else
        write(101,'(7F15.9,1X,I7.7)')  &
        coord_nod2D(1,i)/rad,coord_nod2D(2,i)/rad, &
        0,0,0,0,0,-999
      end if
    end do
  end do
end subroutine read_section_from_netcdf

subroutine read_thalweg_from_nc

  integer            :: hnode,layer
  real(r8)           :: thalweglon,thalweglat,dist
  character(6)       :: DAYNUM
  character(7)       :: OUTDIR
  character(100)     :: OUTFILENAM,THALWEGFILE

  do day2ext=iniday,endday
  write(DAYNUM,'(a,i3.3)') 'DAY',day2ext
  print*,DAYNUM
  call oce_input
  OUTDIR='THALWEG'
  OUTFILENAM=trim(OUTDIR)//'/'//trim(OUTDIR)//'_'//runid//'_'//runyear//'_'//DAYNUM//'.asc'
  THALWEGFILE=thalweg_directory//'/tssthalweg.lst'
  open(unit=101,file=OUTFILENAM,status='replace',access='append',form='formatted')
  open(unit=102,file=THALWEGFILE)
  do i=1,636
     read(102,*) hnode,thalweglon,thalweglat,dist
     do layer=1,110
       if ( nod3D_below_nod2D(layer,hnode).ne.-999 ) then
         write(101,'(2i10,3F10.5,3i10,2F10.5)') i,hnode,dist,thalweglon,thalweglat,layer,nod3D_below_nod2D(layer,hnode), &
              layerdepth(layer),tracer((nod3D_below_nod2D(layer,hnode)),1:2)
       else
         write(101,'(2i10,3F10.5,3i10,A17)') i,hnode,dist,thalweglon,thalweglat,layer,nod3D_below_nod2D(layer,hnode), &
              layerdepth(layer)," -9999.9 -9999.9 "
       endif
     end do
  end do
  close(102)
  close(101)
  end do
end subroutine read_thalweg_from_nc


subroutine calc_section_annual_mean

  character(6)      :: DAYNAM,LEVNAM
  character(4)      :: YEAR
  character(9)      :: OUTDIR
  character(120)    :: OUTFILENAM
  integer           :: layer, totalday
  real(kind=8), allocatable, dimension(:)    :: uvel, vvel, elev, KE
  real(kind=8), allocatable, dimension(:)    :: temp, salt, dens
  allocate(temp(myDim_nod2D))
  allocate(salt(myDim_nod2D))
  allocate(dens(myDim_nod2D))
  allocate(uvel(myDim_nod2D))
  allocate(vvel(myDim_nod2D))
  allocate(elev(myDim_nod2D))
  allocate(KE(myDim_nod2D))

  OUTDIR='.'
  layer=level_number
  totalday=endday
  temp=0; salt=0; dens=0; uvel=0; vvel=0; elev=0

  write(LEVNAM,'(a,i3.3)')'LEV',layer
  do day2ext=1,totalday
    write(DAYNAM,'(a,i3.3)')'DAY',day2ext
    print*,'DAY: ', DAYNAM
    call oce_input
    call compute_density
    do i=1,myDim_nod2D
        temp(i)=temp(i)+tracer((nod3D_below_nod2D(layer,i)),1)
        salt(i)=salt(i)+tracer((nod3D_below_nod2D(layer,i)),2)
        dens(i)=dens(i)+density_insitu(nod3D_below_nod2D(layer,i))
        uvel(i)=uvel(i)+uf(nod3d_below_nod2d(layer,i))
        vvel(i)=vvel(i)+uf(nod3d_below_nod2d(layer,i)+myDim_nod3D)
        elev(i)=elev(i)+ssh(i)
        KE(i)=KE(i)+(uf(nod3d_below_nod2d(layer,i))**2 + uf(nod3d_below_nod2d(layer,i)+myDim_nod3D)**2)/2
    end do
  end do
    OUTFILENAM=trim(OUTDIR)//'/CALCSMEAN_'//runid//'_'//runyear//'_ANNUAL_'//LEVNAM//'.asc'
    open(unit=101,file=OUTFILENAM,status='replace',access='append',form='formatted')
    do i=1,myDim_nod2D
      if ( nod3D_below_nod2D(layer,i).ge.1.and.nod3D_below_nod2D(layer,i).le.myDim_nod3D ) then
        write(101,'(9F15.9,1X,I24)')  &
        coord_nod2D(1,i)/rad,coord_nod2D(2,i)/rad, &
        temp(i)/totalday, salt(i)/totalday, dens(i)/totalday, &
        uvel(i)/totalday, vvel(i)/totalday, elev(i)/totalday, &
        KE(i)/totalday, nod3D_below_nod2D(layer,i)
      else
        write(101,'(2F15.9,1X,8I7.5)')  &
        coord_nod2D(1,i)/rad,coord_nod2D(2,i)/rad, &
        0,0,1000,0,0,0,0,-999
      end if
    end do
  deallocate(temp)
  deallocate(salt)
  deallocate(dens)
  deallocate(uvel)
  deallocate(vvel)
  deallocate(elev)
  deallocate(KE)

end subroutine calc_section_annual_mean

subroutine calc_section_monthly_mean

  character(6)      :: DAYNAM,LEVNAM,MONNAM
  character(4)      :: YEAR
  character(9)      :: OUTDIR
  character(120)    :: OUTFILENAM
  integer           :: layer,totalday, FSTDAY, LSTDAY
  integer           :: MNT(13), YR
  real(kind=8), allocatable, dimension(:)    :: uvel, vvel, elev
  real(kind=8), allocatable, dimension(:)    :: temp, salt, dens
  allocate(temp(myDim_nod2D))
  allocate(salt(myDim_nod2D))
  allocate(dens(myDim_nod2D))
  allocate(uvel(myDim_nod2D))
  allocate(vvel(myDim_nod2D))
  allocate(elev(myDim_nod2D))
  write(YR,'(a4)') runyear
  if ( runyear == '2008' ) then
          MNT=(/ 1, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
  else
          MNT=(/ 1, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
  end if
  OUTDIR='.'
  layer=level_number
  do j=1,12
  temp=0
  salt=0
  dens=0
  uvel=0
  vvel=0
  elev=0
  FSTDAY=sum(MNT(1:j))
  LSTDAY=sum(MNT(2:j+1))
  totalday=LSTDAY-FSTDAY+1
  write(LEVNAM,'(a,i3.3)')'LEV',layer
  write(MONNAM,'(a,i2.2)')'MNTH',j
  print*, MONNAM
  do day2ext=FSTDAY,LSTDAY
    write(DAYNAM,'(a,i3.3)')'DAY',day2ext
    print*,'DAY: ', DAYNAM
    call oce_input
    call compute_density
    do i=1,myDim_nod2D
        temp(i)=temp(i)+tracer((nod3D_below_nod2D(layer,i)),1)
        salt(i)=salt(i)+tracer((nod3D_below_nod2D(layer,i)),2)
        dens(i)=dens(i)+density_insitu(nod3D_below_nod2D(layer,i))
        uvel(i)=uvel(i)+uf(nod3d_below_nod2d(layer,i))
        vvel(i)=vvel(i)+uf(nod3d_below_nod2d(layer,i)+myDim_nod3D)
        elev(i)=elev(i)+ssh(i)
    end do
  end do
    OUTFILENAM=trim(OUTDIR)//'/CALCSMEAN_'//runid//'_'//runyear//'_'//MONNAM//'_'//LEVNAM//'.asc'
    open(unit=101,file=OUTFILENAM,status='replace',access='append',form='formatted')
    do i=1,myDim_nod2D
      if ( nod3D_below_nod2D(layer,i).ge.1.and.nod3D_below_nod2D(layer,i).le.myDim_nod3D ) then
        write(101,'(8F15.9,1X,I24)')  &
        coord_nod2D(1,i)/rad,coord_nod2D(2,i)/rad, &
        temp(i)/totalday, salt(i)/totalday, dens(i)/totalday, &
        uvel(i)/totalday, vvel(i)/totalday, elev(i)/totalday, &
        nod3D_below_nod2D(layer,i)
      else
        write(101,'(2F15.9,1X,7I7.5)')  &
        coord_nod2D(1,i)/rad,coord_nod2D(2,i)/rad, &
        0,0,1000,0,0,0,-999
      end if
    end do
end do
  deallocate(temp)
  deallocate(salt)
  deallocate(dens)
  deallocate(uvel)
  deallocate(vvel)
  deallocate(elev)

end subroutine calc_section_monthly_mean

subroutine calc_thalweg_annual_mean

  integer           :: hnode,layer, totalday
  real(r8)          :: thalweglon,thalweglat,dist
  character(6)      :: DAYNAM,MONNAM
  character(9)      :: OUTDIR
  character(100)    :: OUTFILENAM,THALWEGFILE
  real(kind=8), allocatable, dimension(:,:)  :: temp,salt, dens
  integer, parameter :: npoint=1063
  integer, parameter :: nlevel=110

  allocate(temp(npoint,nlevel))
  allocate(salt(npoint,nlevel))
  allocate(dens(npoint,nlevel))
  OUTDIR='.'
  THALWEGFILE=thalweg_directory//'/tssthalweg.lst'
  temp = 0; salt = 0; dens = 0
  totalday=endday

  do day2ext=1,totalday
  write(DAYNAM,'(a,i3.3)')'DAY',day2ext
  print*,DAYNAM
  call oce_input
  call compute_density
  open(unit=102,file=THALWEGFILE)
  do i=1,npoint
     read(102,*) hnode,thalweglon,thalweglat
     do layer=1,nlevel
       if ( nod3D_below_nod2D(layer,hnode).ne.-999 ) then
       temp(i,layer)=temp(i,layer)+tracer((nod3D_below_nod2D(layer,hnode)),1)
       salt(i,layer)=salt(i,layer)+tracer((nod3D_below_nod2D(layer,hnode)),2)
       dens(i,layer)=dens(i,layer)+density_insitu(nod3D_below_nod2D(layer,hnode))
       else
               temp(i,layer)=temp(i,layer)-9999.9
               salt(i,layer)=salt(i,layer)-9999.9
               dens(i,layer)=dens(i,layer)-9999.9
       end if
     end do
  end do
  close(102)
  end do
  OUTFILENAM=trim(OUTDIR)//'/CALCTMEAN_'//runid//'_'//runyear//'_ANNUAL.asc'
  open(unit=102,file=THALWEGFILE)
  open(unit=101,file=OUTFILENAM,status='replace',access='append',form='formatted')
  do i=1,npoint
     read(102,*) hnode,thalweglon,thalweglat,dist
     do layer=1,nlevel
         write(101,'(2i10,3F10.5,3i10,3F16.5)') i,hnode, dist, &
              thalweglon,thalweglat,layer,nod3D_below_nod2D(layer,hnode), &
              layerdepth(layer),temp(i,layer)/totalday,salt(i,layer)/totalday,dens(i,layer)/totalday
     end do
  end do
  close(102)
  close(101)

  deallocate(temp)
  deallocate(salt)
  deallocate(dens)
end subroutine calc_thalweg_annual_mean

subroutine calc_thalweg_monthly_mean

  integer           :: hnode,layer
  real(r8)          :: thalweglon,thalweglat,dist
  character(6)      :: DAYNAM,MONNAM
  character(9)      :: OUTDIR
  character(100)    :: OUTFILENAM,THALWEGFILE
  integer           :: totalday, FSTDAY, LSTDAY
  integer           :: MNT(13), YR
  real(kind=8), allocatable, dimension(:,:)  :: temp,salt, dens
  integer, parameter :: npoint=1063
  integer, parameter :: nlevel=110

  allocate(temp(npoint,nlevel))
  allocate(salt(npoint,nlevel))
  allocate(dens(npoint,nlevel))
  write(YR,'(a4)') runyear
  if ( runyear == '2008'  ) then
          MNT=(/ 1, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
  else
          MNT=(/ 1, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
  end if
  OUTDIR='.'
  THALWEGFILE=thalweg_directory//'/tssthalweg.lst'
  do j=1,12
  temp = 0
  salt = 0
  dens = 0
  FSTDAY=sum(MNT(1:j))
  LSTDAY=sum(MNT(2:j+1))
  totalday=LSTDAY-FSTDAY+1
  write(MONNAM,'(a,i2.2)')'MNTH',j
  print*, MONNAM
  do day2ext=FSTDAY,LSTDAY
  write(DAYNAM,'(a,i3.3)')'DAY',day2ext
  print*,DAYNAM
  call oce_input
  call compute_density
  open(unit=102,file=THALWEGFILE)
  do i=1,npoint
     read(102,*) hnode,thalweglon,thalweglat
     do layer=1,nlevel
       if ( nod3D_below_nod2D(layer,hnode).ne.-999 ) then
       temp(i,layer)=temp(i,layer)+tracer((nod3D_below_nod2D(layer,hnode)),1)
       salt(i,layer)=salt(i,layer)+tracer((nod3D_below_nod2D(layer,hnode)),2)
       dens(i,layer)=dens(i,layer)+density_insitu(nod3D_below_nod2D(layer,hnode))
       else
               temp(i,layer)=temp(i,layer)-9999.9
               salt(i,layer)=salt(i,layer)-9999.9
               dens(i,layer)=dens(i,layer)-9999.9
       end if
     end do
  end do
  close(102)
  end do
  OUTFILENAM=trim(OUTDIR)//'/CALCTMEAN_'//runid//'_'//runyear//'_'//MONNAM//'.asc'
  open(unit=102,file=THALWEGFILE)
  open(unit=101,file=OUTFILENAM,status='replace',access='append',form='formatted')
  do i=1,npoint
     read(102,*) hnode,thalweglon,thalweglat,dist
     do layer=1,nlevel
         write(101,'(2i10,3F10.5,3i10,3F16.5)') i,hnode, dist, &
              thalweglon,thalweglat,layer,nod3D_below_nod2D(layer,hnode), &
              layerdepth(layer),temp(i,layer)/totalday,salt(i,layer)/totalday,dens(i,layer)/totalday
     end do
  end do
  close(102)
  close(101)

  end do
  deallocate(temp)
  deallocate(salt)
  deallocate(dens)
end subroutine calc_thalweg_monthly_mean

subroutine velocity_at_the_exit

  use  fesom_observation_mod, only   : find_close_2D, close_node

  integer           :: region,hnode(4),layer
  real(r8)          :: vlon(4),vlat(4), swap_uv
  real(r8)          :: uvel, vvel, tvel, theta, rtheta
  character(3)      :: DAYNUM
  character(7)      :: OUTDIR
  character(100)    :: OUTFILENAM,COORDINATES, LOCATION

  OUTDIR='.'
  OUTFILENAM=trim(OUTDIR)//'/VVEL_'//runid//'_'//runyear//'.lst'
  COORDINATES="/users/home/ans051/FEOM_POSTPROC/MESH_READ/COORD_VEL.lst"
  open(unit=101,file=OUTFILENAM,status='replace',access='append',form='formatted')
  do day2ext=iniday,endday
  write(DAYNUM,'(i3.3)') day2ext
  call oce_input
  open(unit=102,file=COORDINATES)
  do region=1,4
     if ( region.eq.1 ) LOCATION="Northern Bosphorus"
     if ( region.eq.2 ) LOCATION="Southern Bosphorus"
     if ( region.eq.3 ) LOCATION="Northern Dardanelles"
     if ( region.eq.4 ) LOCATION="Southern Dardanelles"
     read(102,*) vlon(region),vlat(region)
     call find_close_2D(vlon(region),vlat(region))
     hnode(region) = close_node
     do layer = 1,70
       uvel   = uf(nod3d_below_nod2d(layer,hnode(region)))
       vvel   = uf(nod3d_below_nod2d(layer,hnode(region))+myDim_nod3D)
!          tvel   = sqrt(uvel*uvel+vvel*vvel)
!          theta  = asin(vvel/tvel)
         if ( region .eq. 1 ) then
              rtheta = -pi/6
              uvel = uvel * cos(rtheta) - vvel * sin(rtheta)
              vvel = uvel * sin(rtheta) + vvel * cos(rtheta)
         else if ( region .eq. 2 ) then
              uvel = uvel
              vvel = vvel
         else if ( region .eq. 3 ) then
              rtheta = -pi/4
              uvel = uvel * cos(rtheta) - vvel * sin(rtheta)
              vvel = uvel * sin(rtheta) + vvel * cos(rtheta)
         else if ( region .eq. 4 ) then
              swap_uv = uvel
              uvel = vvel
              vvel = swap_uv
         end if


       if ( nod3D_below_nod2D(layer,hnode(region)).ne.-999 ) then
         write(101,'(3I5,A30,6F10.5)') day2ext,layer,region,trim(LOCATION), &
         uvel, vvel, theta, rtheta, &
         tracer(nod3D_below_nod2D(layer,hnode(region)),1), &
         tracer(nod3D_below_nod2D(layer,hnode(region)),2)
       else
         write(101,'(3I5,A30,A17,2F10.5,A17)') day2ext,layer,region,trim(LOCATION), &
         " NaN NaN ",theta, rtheta, " NaN NaN "
       endif
     end do
  end do
  close(102)
  end do
  close(101)
end subroutine velocity_at_the_exit

subroutine dardanelles_for_MFS

  integer           :: hnode(4),layer
  character(3)      :: DAYNUM
  character(11)     :: OUTDIR
  character(100)    :: OUTFILENAM,COORDINATES, LOCATION
  character(100)                 :: SALT_NODE_FILE
  integer                        :: salt_nodes, node, nodes
  integer,allocatable            :: snode(:)
  real(kind=8), allocatable      :: dummy1(:), dummy2(:)
  real(kind=8), allocatable      :: saltlon(:), saltlat(:)

  OUTDIR='DARDANELLES'
  SALT_NODE_FILE=trim(thalweg_directory)//'/SECTION_SD.lst'
  OUTFILENAM=trim(OUTDIR)//'/'//trim(OUTDIR)//'_'//runid//'_'//runyear//'.lst'

  open(unit=101,file=OUTFILENAM,status='replace',access='append',form='formatted')
  write(101,*)  'year: ', runyear
  write(101,*)  'day  lon  lat  depth  temp salt'
day_loop: do day2ext = iniday, endday
  write(DAYNUM,'(i3.3)') day2ext

  call oce_input
  open(unit=122,file=trim(SALT_NODE_FILE))
  read(122,*)  salt_nodes
  allocate(saltlon(salt_nodes), saltlat(salt_nodes))
  allocate(dummy1(salt_nodes), dummy2(salt_nodes))
  allocate(snode(salt_nodes))

  do nodes = 1,salt_nodes
   read(122,*) saltlon(nodes),saltlat(nodes), &
               dummy1(nodes),dummy2(nodes),snode(nodes)
  end do

node_loop: do nodes = 1,salt_nodes

  node = snode(nodes)

layer_loop: do layer = 1,110
       if ( nod3d_below_nod2d(layer,node).ne.-999 ) then
         write(101,'(I5,2F10.5,I7,2F10.5)')  day2ext, coord_nod2d(1,node)/rad, coord_nod2d(2,node)/rad, layerdepth(layer), tracer(nod3D_below_nod2D(layer,node),1), tracer(nod3D_below_nod2D(layer,node),2)
       endif
     end do layer_loop
  end do node_loop
  close(122)
  deallocate(saltlon, saltlat, dummy1, dummy2, snode)
  end do day_loop
  close(101)
end subroutine dardanelles_for_MFS

subroutine surface_kinetic_energy

  integer           :: hnode(4),layer, elnodes2(3)
  real(r8)          :: vlatlon(4)
  real(r8)          :: uvel, vvel, tke, tkn, totvol
  real(r8)          :: inv3=1.0/3.0_8
  real(r8)          :: inv4=1.0/4.0_8
  real(r8)          :: inv12=1.0/12.0_8
  character(3)      :: DAYNUM
  character(3)      :: OUTDIR
  character(100)    :: OUTFILENAM,COORDINATES, LOCATION

  OUTDIR='SKE'
  OUTFILENAM=trim(OUTDIR)//'/'//trim(OUTDIR)//'_'//runid//'_'//runyear//'.lst'

  open(unit=101,file=OUTFILENAM,status='replace',access='append',form='formatted')

  day_loop:   do day2ext = iniday, endday
  write(DAYNUM,'(i3.3)') day2ext
  totvol=0
  tke=0

  call oce_input

  elem_loop: do i=1,myDim_elem2D
       elnodes2=elem2D_nodes(:,i)
       do k=1,3
         if ( (coord_nod2D(2,elnodes2(k))/rad).gt.marmax_lat .or. (coord_nod2D(1,elnodes2(k))/rad).gt.marmax_lon &
         .or. (coord_nod2D(2,elnodes2(k))/rad).lt.marmin_lat .or. (coord_nod2D(1,elnodes2(k))/rad).lt.marmin_lon ) then
            go to 1234
         endif
       enddo
          uvel=0
          vvel=0
          do k=1,3
          if ( (coord_nod2D(2,elnodes2(k))/rad).gt.41 .and. (coord_nod2D(1,elnodes2(k))/rad).gt.28.7 ) then
                  go to 1234
          else
          uvel = uvel + uf(elnodes2(k))
          vvel = vvel + uf(elnodes2(k)+myDim_nod3D)
          endif
          enddo
          totvol = totvol + sum(cluster_area_2d(elnodes2(:))) * inv3
          uvel = uvel * inv3
          vvel = vvel * inv3
          tke  = tke + ( uvel**2 + vvel**2 ) / 2 * sum(cluster_area_2d(elnodes2(:))) * inv3
        1234 continue
  enddo elem_loop
     write(101,'(I5,F20.15)')  day2ext, tke / totvol
  enddo day_loop
  close(101)
end subroutine surface_kinetic_energy


subroutine total_kinetic_energy

  integer           :: hnode(4),layer, elnodes3(4)
  real(r8)          :: vlatlon(4)
  real(r8)          :: uvel, vvel, tke, tkn, totvol
  real(r8)          :: inv4=1.0/4.0_8
  real(r8)          :: inv12=1.0/12.0_8
  character(3)      :: DAYNUM
  character(3)      :: OUTDIR
  character(100)    :: OUTFILENAM,COORDINATES, LOCATION

  OUTDIR='TKE'
  OUTFILENAM=trim(OUTDIR)//'/'//trim(OUTDIR)//'_'//runid//'_'//runyear//'.lst'

  open(unit=101,file=OUTFILENAM,status='replace',access='append',form='formatted')

  day_loop:   do day2ext = iniday, endday
  write(DAYNUM,'(i3.3)') day2ext

  call oce_input

  tke=0
  totvol=0
  elem_loop: do i=1,myDim_elem3D
       elnodes3=elem3D_nodes(:,i)
       uvel = 0
       vvel = 0
       do k=1,4
         if ( (coord_nod3D(2,elnodes3(k))/rad).gt.marmax_lat .or. (coord_nod3D(1,elnodes3(k))/rad).gt.marmax_lon &
         .or. (coord_nod3D(2,elnodes3(k))/rad).lt.marmin_lat .or. (coord_nod3D(1,elnodes3(k))/rad).lt.marmin_lon ) then
            go to 1235
         endif
       enddo
          do k=1,4
          if ( (coord_nod3D(2,elnodes3(k))/rad).gt.41 .and. (coord_nod3D(1,elnodes3(k))/rad).gt.28.7 ) then
                  uvel = 0
                  vvel = 0
                  go to 1235
          else
                 uvel = uvel + uf(elnodes3(k))
                 vvel = vvel + uf(elnodes3(k)+myDim_nod3D)
          endif
          end do
          totvol = totvol + sum(voltetra(elnodes3(:))) * inv4
          uvel   = uvel * inv4
          vvel   = vvel * inv4
          tke    = tke + ( uvel**2 + vvel**2 ) /2  * sum(voltetra(elnodes3(:))) * inv4
        1235 continue
  enddo elem_loop
     write(101,'(I5,F20.15)')  day2ext, tke / totvol
  enddo day_loop
  close(101)
end subroutine total_kinetic_energy

subroutine compute_vorticity

  integer                        :: m, elem, el2, elnodes(4)
  integer                        :: row, row2, ind, n3, q
  real(kind=8)                   :: aux, vc(4), dx(4), dy(4), npx, npy
  real(kind=8)                   :: udx, udy, u_el(4), v_el(4)
  real(kind=8)                   :: vdx, vdy
  real(kind=8)                   :: inv4
  real(kind=8),allocatable       :: vorticity(:), dxvdyu(:), vort(:)
  real(kind=8)                   :: cori_p, dparam, beta, gamma

  character(6)      :: DAYNAM,LEVNAM
  character(4)      :: YEAR
  character(9)      :: OUTDIR
  character(120)    :: OUTFILENAM
  integer           :: layer, totalday

  totalday=endday
  OUTDIR='.'
  layer=level_number

  inv4=0.25_8


  allocate(vorticity(myDim_nod3d))
  allocate(vort(myDim_nod2d))
  allocate(dxvdyu(myDim_elem3d))

  vort=0
  vorticity=0

  write(LEVNAM,'(a,i3.3)')'LEV',layer
  do day2ext=1,totalday
    write(DAYNAM,'(a,i3.3)')'DAY',day2ext
    print*,'DAY: ', DAYNAM
    call oce_input
  do elem=1, myDim_elem3d

     elnodes=elem3D_nodes(:,elem)

     dx=bafux_3d(:, elem)
     dy=bafuy_3d(:, elem)
     u_el=uf(elnodes)
     v_el=uf(myDim_nod3D+elnodes)

!        udx=sum(dx*u_el)
     udy=sum(dy*u_el)
     vdx=sum(dx*v_el)
!        vdy=sum(dy*v_el)
!
!        dxvdyu(elem)=sum(u_el)
     do q=1,4
       row=elnodes(q)
       vorticity(row) = vdx - udy
     end do

   end do


    do i=1,myDim_nod2D
        vort(i)=vort(i)+vorticity(nod3D_below_nod2D(layer,i))
    end do
   end do
   OUTFILENAM=trim(OUTDIR)//'/VORTSMEAN_'//runid//'_'//runyear//'_ANNUAL_'//LEVNAM//'.asc'
    open(unit=101,file=OUTFILENAM,status='replace',access='append',form='formatted')
    do i=1,myDim_nod2D
      if ( nod3D_below_nod2D(layer,i).ge.1.and.nod3D_below_nod2D(layer,i).le.myDim_nod3D ) then
        write(101,'(3F20.12,1X,I24)')  &
        coord_nod2D(1,i)/rad,coord_nod2D(2,i)/rad, &
        vort(i)/totalday, nod3D_below_nod2D(layer,i)
      else
        write(101,'(2F20.15,1X,2I7.5)')  &
        coord_nod2D(1,i)/rad,coord_nod2D(2,i)/rad, &
        0,-999
      end if
    end do

   deallocate (vorticity, dxvdyu, vort)
   close(101)
end subroutine compute_vorticity

 subroutine compute_volume_transport

  use utilities, only   : earthdistance, to_radian

  character(6)                   :: DAYNAM,LEVNAM, DAYNUM
  character(100)                 :: FLUX_NODE_FILE
  character(4)                   :: YEAR
  character(9)                   :: OUTDIR
  character(122)                 :: OUTFILENAM
  integer                        :: nodes, dep
  integer                        :: layer
  integer                        :: minlev(3)
  integer                        :: flux_nodes
  integer,allocatable            :: snode(:)
  integer                        :: totalday
  real(kind=8)                   :: inv2, inv3, inv4
  real(kind=8)                   :: depth_of_layer, AREA
  real(kind=8)                   :: lon1,lon2,lat1,lat2,ndist,totdist
  real(kind=8)                   :: slope, alfa
  real(kind=8), allocatable      :: dummy1(:), dummy2(:)
  real(kind=8), allocatable      :: fluxlon(:), fluxlat(:)
  real(kind=8)                   :: locflux, netflux, upperflux, lowerflux
  real(kind=8)                   :: tnflux, tlflux, tuflux
  real(kind=8)                   :: salt(4), salt_interp, halocline
  real(kind=8)                   :: uv(2,4), uv_rot(2,4), uv_interp(2)

  inv4=0.25_8
  inv3=1/3
  inv2=0.5_8

  totalday=endday
  OUTDIR='.'
  layer=level_number

  tnflux=0
  tuflux=0
  tlflux=0

  if ( flux_section == 'NB' ) halocline = 19
  if ( flux_section == 'MB' ) halocline = 22
  if ( flux_section == 'SB' ) halocline = 25
  if ( flux_section == 'ND' ) halocline = 28
  if ( flux_section == 'SD' ) halocline = 36

OUTFILENAM=trim(OUTDIR)//'/FLX_'//runid//'_'//runyear//'_'//flux_section//'_VT11.lst'
  if(mype==0) write(*,*) "output file: ",OUTFILENAM
FLUX_NODE_FILE=trim(thalweg_directory)//'/SECTION_'//flux_section//'.lst'
open(unit=121,file=OUTFILENAM,status='replace',access='append',form='formatted')
  do day2ext=1,totalday

  netflux=0
  upperflux=0
  lowerflux=0

write(DAYNUM,'(a,i3.3)')'DAY',day2ext
print*,DAYNUM

call oce_input

open(unit=122,file=trim(FLUX_NODE_FILE))
read(122,*)  flux_nodes
allocate(fluxlon(flux_nodes), fluxlat(flux_nodes))
allocate(dummy1(flux_nodes), dummy2(flux_nodes))
allocate(snode(flux_nodes))

do nodes = 1,flux_nodes
 read(122,*) fluxlon(nodes),fluxlat(nodes), &
             dummy1(nodes),dummy2(nodes),snode(nodes)
end do

totdist = 0
do nodes = 1,flux_nodes-1
!   do nodes = 1,4
ndist = 0

lon1 = fluxlon(nodes)
lon2 = fluxlon(nodes+1)
lat1 = fluxlat(nodes)
lat2 = fluxlat(nodes+1)

slope = (to_radian(lat2)-to_radian(lat1))/(to_radian(lon2)-to_radian(lon1))
alfa  = atan( slope )
if ( day2ext == 2 .and. mype==0) write(*,*) "----------------------------"
if ( day2ext == 2 .and. mype==0) write(*,*) "slope: ",slope
if ( day2ext == 2 .and. mype==0) write(*,*) "alfa: ", alfa, alfa / rad
if (flux_section == "NB" ) alfa  = alfa
if (flux_section == "SB" ) alfa  = alfa
if (flux_section == "ND" ) alfa  = alfa
if (flux_section == "SD" .and. alfa < 0 ) alfa  = alfa
if (flux_section == "SD" .and. alfa > 0 ) alfa  = -(pi - alfa)
if ( day2ext == 2 .and. mype==0) write(*,*) "alfa: ", alfa, alfa / rad

if ( day2ext == 2 .and. mype==0) write(*,*) "lon1: ", lon1, "lat1: ", lat1
if ( day2ext == 2 .and. mype==0) write(*,*) "lon2: ", lon2, "lat2: ", lat2
if ( day2ext == 2 .and. mype==0) write(*,*) "node: ", snode(nodes),"-", snode(nodes+1)

ndist = earthdistance(lat1,lat2,lon1,lon2)
totdist = totdist + ndist
if ( day2ext == 2 .and. mype==0) write(*,*) "dist: ", ndist

minlev(1) = num_layers_below_nod2D(snode(nodes))
minlev(2) = num_layers_below_nod2D(snode(nodes+1))
minlev(3) = min(minlev(1), minlev(2))
if ( day2ext == 2 .and. mype==0) write(*,*) "minimum level: ", minlev(:)
do dep = 1,max_num_layers
 if ( dep .gt. minlev(3) ) then
    go to 1236
 endif
 if ( dep == 1 ) then
   depth_of_layer=layerdepth(dep+1)-layerdepth(dep)+(ssh(snode(nodes))+ssh(snode(nodes+1)))/2
   AREA = depth_of_layer * ndist
 elseif ( dep == minlev(3) ) then
   depth_of_layer=abs(layerdepth(minlev(1))-layerdepth(minlev(2)))
   AREA = depth_of_layer * ndist / 2
 else
   depth_of_layer=layerdepth(dep+1)-layerdepth(dep)
   AREA = depth_of_layer * ndist
 endif
 if ( day2ext == 2 .and. mype==0) write(*,*) "AREA: ", AREA
 if ( day2ext == 2 .and. mype==0) write(*,*) "layer depth: ", depth_of_layer, 'depth', dep

! rotate and interpolate the uv velocities
 salt(1) = tracer(nod3D_below_nod2D(dep,snode(nodes)),2)
 salt(2) = tracer(nod3D_below_nod2D(dep,snode(nodes+1)),2)
 salt(3) = tracer(nod3D_below_nod2D(dep+1,snode(nodes)),2)
 salt(4) = tracer(nod3D_below_nod2D(dep+1,snode(nodes+1)),2)
 uv(1,1) = uf(nod3D_below_nod2D(dep,snode(nodes)))
 uv(1,2) = uf(nod3D_below_nod2D(dep,snode(nodes+1)))
 uv(1,3) = uf(nod3D_below_nod2D(dep+1,snode(nodes)))
 uv(1,4) = uf(nod3D_below_nod2D(dep+1,snode(nodes+1)))
 uv(2,1) = uf(myDim_nod3d+nod3D_below_nod2D(dep,snode(nodes)))
 uv(2,2) = uf(myDim_nod3d+nod3D_below_nod2D(dep,snode(nodes+1)))
 uv(2,3) = uf(myDim_nod3d+nod3D_below_nod2D(dep+1,snode(nodes)))
 uv(2,4) = uf(myDim_nod3d+nod3D_below_nod2D(dep+1,snode(nodes+1)))
!
 if ( dep >= minlev(3) )  then
   uv      = 0.0
!      uv(1,1) = uf(nod3D_below_nod2D(minlev(1),snode(nodes)))
!      uv(1,2) = uf(nod3D_below_nod2D(minlev(2),snode(nodes+1)))
!      uv(1,3) = 0
!      uv(1,4) = 0
!      uv(2,1) = uf(myDim_nod3d+nod3D_below_nod2D(minlev(1),snode(nodes)))
!      uv(2,2) = uf(myDim_nod3d+nod3D_below_nod2D(minlev(2),snode(nodes+1)))
!      uv(2,3) = 0
!      uv(2,4) = 0
 endif
!
 if ( day2ext == 2 .and. mype==0) write(*,*) 'uv11',uv(1,1)
 if ( day2ext == 2 .and. mype==0) write(*,*) 'uv12',uv(1,2)
 if ( day2ext == 2 .and. mype==0) write(*,*) 'uv13',uv(1,3)
 if ( day2ext == 2 .and. mype==0) write(*,*) 'uv14',uv(1,4)
 if ( day2ext == 2 .and. mype==0) write(*,*) 'uv21',uv(2,1)
 if ( day2ext == 2 .and. mype==0) write(*,*) 'uv22',uv(2,2)
 if ( day2ext == 2 .and. mype==0) write(*,*) 'uv23',uv(2,3)
 if ( day2ext == 2 .and. mype==0) write(*,*) 'uv24',uv(2,4)
 do j=1,4
     uv_rot(1,j) = uv(1,j) * cos(alfa) + uv(2,j) * sin(alfa)
     uv_rot(2,j) = -uv(1,j) * sin(alfa) + uv(2,j) * cos(alfa)
 end do

 uv_interp(1) = sum(uv_rot(1,:))*inv4
 uv_interp(2) = sum(uv_rot(2,:))*inv4
 salt_interp  = sum(salt(:))*inv4

 if ( dep == minlev(3) )  then
      uv_interp(1) = sum(uv_rot(1,:))*inv2
      uv_interp(2) = sum(uv_rot(2,:))*inv2
      salt_interp  = sum(salt(1:2))*inv2
 endif

 if ( day2ext == 2 .and. mype==0) write(*,*) 'u_interp', uv_interp(1)
 if ( day2ext == 2 .and. mype==0) write(*,*) 'v_interp', uv_interp(2)
 if ( day2ext == 2 .and. mype==0) write(*,*) 'salt_interp', salt_interp

 if ( uv_interp(2) < 0 ) then
         locflux = AREA * uv_interp(2)
         upperflux = upperflux + locflux
 else if ( uv_interp(2) > 0 ) then
         locflux = AREA * uv_interp(2)
         lowerflux = lowerflux + locflux
 else
         locflux = 0
 endif
 locflux = AREA * uv_interp(2)
!
 netflux = netflux + locflux

 end do
 1236 continue
  end do

  write(121,'(3F20.9)') netflux/km3yr2m3sec,   &
                     upperflux/km3yr2m3sec, &
                     lowerflux/km3yr2m3sec

  tnflux = tnflux + netflux
  tuflux = tuflux + upperflux
  tlflux = tlflux + lowerflux


!    first determine a cross-section by listing the surface nodes
!    second save the u,v corresponding to those nodes and nodes3d corresponding
!    third find the distance between consecutive nodes
!    calculate the surface area by multiplying with the corresp. depth
!    find the angle of the surface to the horizontal
!    interpolate u,v to the center of the rotated surface
!    compute the flux by (u,v).AREA

  close(122)

  deallocate(fluxlon, fluxlat)
  deallocate(dummy1, dummy2)
  deallocate(snode)

enddo
 if ( day2ext == 2 .and. mype==0) write(*,*)  "total distance: ", totdist
 if ( day2ext == 2 .and. mype==0) write(*,*)  "-------------------------"


 write(121,'(3F20.9)') tnflux/km3yr2m3sec/totalday, &
                    tuflux/km3yr2m3sec/totalday, &
                    tlflux/km3yr2m3sec/totalday
 close(121)

 end subroutine compute_volume_transport

 subroutine compute_net_flux

  use utilities, only   : earthdistance, to_radian

  character(6)                   :: DAYNAM,LEVNAM, DAYNUM
  character(100)                 :: FLUX_NODE_FILE
  character(4)                   :: YEAR
  character(9)                   :: OUTDIR
  character(120)                 :: OUTFILENAM
  integer                        :: nodes, dep
  integer                        :: layer,minlev
  integer                        :: flux_nodes
  integer,allocatable            :: snode(:)
  integer                        :: totalday
  real(kind=8)                   :: inv2,inv4
  real(kind=8)                   :: depth_of_layer, AREA
  real(kind=8)                   :: lon1,lon2,lat1,lat2,ndist,totdist
  real(kind=8)                   :: slope, alfa
  real(kind=8), allocatable      :: dummy1(:), dummy2(:)
  real(kind=8), allocatable      :: fluxlon(:), fluxlat(:)
  real(kind=8)                   :: locflux, netflux, upperflux, lowerflux
  real(kind=8)                   :: tnflux, tlflux, tuflux
  real(kind=8)                   :: salt(4), salt_interp, halocline
  real(kind=8)                   :: uv(2,4), uv_rot(2,4), uv_interp(2)

  inv4=0.25_8
  inv2=0.5_8

  totalday=endday
  OUTDIR='.'
  layer=level_number

  tnflux=0
  tuflux=0
  tlflux=0

  if ( flux_section == 'NB' ) halocline = 19
  if ( flux_section == 'MB' ) halocline = 22
  if ( flux_section == 'SB' ) halocline = 25
  if ( flux_section == 'ND' ) halocline = 28
  if ( flux_section == 'SD' ) halocline = 36

OUTFILENAM=trim(OUTDIR)//'/FLX_'//runid//'_'//runyear//'_'//flux_section//'.lst'
  if(mype==0) write(*,*) "output file: ",OUTFILENAM
FLUX_NODE_FILE=trim(thalweg_directory)//'/SECTION_'//flux_section//'.lst'
open(unit=121,file=OUTFILENAM,status='replace',access='append',form='formatted')
  do day2ext=1,totalday

  netflux=0
  upperflux=0
  lowerflux=0

write(DAYNUM,'(a,i3.3)')'DAY',day2ext
print*,DAYNUM

call oce_input

open(unit=122,file=trim(FLUX_NODE_FILE))
read(122,*)  flux_nodes
allocate(fluxlon(flux_nodes), fluxlat(flux_nodes))
allocate(dummy1(flux_nodes), dummy2(flux_nodes))
allocate(snode(flux_nodes))

do nodes = 1,flux_nodes
 read(122,*) fluxlon(nodes),fluxlat(nodes), &
             dummy1(nodes),dummy2(nodes),snode(nodes)
end do

totdist = 0
do nodes = 1,flux_nodes-1
!   do nodes = 1,4
ndist = 0

lon1 = fluxlon(nodes)
lon2 = fluxlon(nodes+1)
lat1 = fluxlat(nodes)
lat2 = fluxlat(nodes+1)

slope = (to_radian(lat2)-to_radian(lat1))/(to_radian(lon2)-to_radian(lon1))
alfa  = atan( slope )
if (flux_section == "NB" ) alfa  = -alfa
if (flux_section == "MB" ) alfa  = -alfa
if (flux_section == "SB" ) alfa  = -alfa
if (flux_section == "ND" ) alfa  = -alfa
if (flux_section == "SD" .and. alfa < 0 ) alfa  = -alfa
if (flux_section == "SD" .and. alfa > 0 ) alfa  = pi - alfa

if ( day2ext == 2 .and. mype==0) write(*,*) "angle: ", alfa, alfa / rad
if ( day2ext == 2 .and. mype==0) write(*,*) "slope: ",slope
if ( day2ext == 2 .and. mype==0) write(*,*) "lon1: ", lon1, "lat1: ", lat1
if ( day2ext == 2 .and. mype==0) write(*,*) "lon2: ", lon2, "lat2: ", lat2

ndist = earthdistance(lat1,lat2,lon1,lon2)

if ( day2ext == 2 .and. mype==0) write(*,*) "dist: ", ndist, "between ",snode(nodes), snode(nodes+1)
if ( day2ext == 2 .and. mype==0) write(*,*) "----------------------------"

minlev = min(num_layers_below_nod2D(snode(nodes)),num_layers_below_nod2D(snode(nodes+1)))
totdist = totdist + ndist
    do dep = 1,max_num_layers
    depth_of_layer=layerdepth(dep+1)-layerdepth(dep)
    AREA = depth_of_layer * ndist
!       if ( dep == minlev ) then
!               AREA = depth_of_layer * ndist / 2
!       endif
if ( day2ext == 2 .and. mype==0) write(*,*) "AREA: ", AREA
if ( day2ext == 2 .and. mype==0) write(*,*) "depth: ", dep

! rotate and interpolate the uv velocities
    salt(1) = tracer(nod3D_below_nod2D(dep,snode(nodes)),2)
    salt(2) = tracer(nod3D_below_nod2D(dep,snode(nodes+1)),2)
    salt(3) = tracer(nod3D_below_nod2D(dep+1,snode(nodes)),2)
    salt(4) = tracer(nod3D_below_nod2D(dep+1,snode(nodes+1)),2)
    uv(1,1) = uf(nod3D_below_nod2D(dep,snode(nodes)))
    uv(1,2) = uf(nod3D_below_nod2D(dep,snode(nodes+1)))
    uv(1,3) = uf(nod3D_below_nod2D(dep+1,snode(nodes)))
    uv(1,4) = uf(nod3D_below_nod2D(dep+1,snode(nodes+1)))
    uv(2,1) = uf(myDim_nod3d+nod3D_below_nod2D(dep,snode(nodes)))
    uv(2,2) = uf(myDim_nod3d+nod3D_below_nod2D(dep,snode(nodes+1)))
    uv(2,3) = uf(myDim_nod3d+nod3D_below_nod2D(dep+1,snode(nodes)))
    uv(2,4) = uf(myDim_nod3d+nod3D_below_nod2D(dep+1,snode(nodes+1)))


    do j=1,4
     uv_rot(1,j) = uv(1,j) * cos(alfa) - uv(2,j) * sin(alfa)
     uv_rot(2,j) = uv(1,j) * sin(alfa) + uv(2,j) * cos(alfa)
    end do

    uv_interp(1) = sum(uv_rot(1,:))*inv4
    uv_interp(2) = sum(uv_rot(2,:))*inv4
    salt_interp  = sum(salt(:))*inv4

!       if ( dep == minlev )  then
!               uv_interp(1) = sum(uv_rot(1,1:2))*inv2
!               uv_interp(2) = sum(uv_rot(2,1:2))*inv2
!               salt_interp  = sum(salt(1:2))*inv2
!       endif

if ( day2ext == 2 .and. mype==0) write(*,*) uv_interp(1)
if ( day2ext == 2 .and. mype==0) write(*,*) uv_interp(2)
if ( day2ext == 2 .and. mype==0) write(*,*) salt_interp

    if ( uv_interp(2) < 0 .and. salt_interp < halocline ) then
            locflux = AREA * uv_interp(2)
            upperflux = upperflux + locflux
    else if ( uv_interp(2) > 0 .and. salt_interp > halocline ) then
            locflux = AREA * uv_interp(2)
            lowerflux = lowerflux + locflux
    else
            locflux = 0
    endif
!
    netflux = netflux + locflux

    end do
 end do

 write(121,'(3F20.9)') netflux/km3yr2m3sec,   &
                       upperflux/km3yr2m3sec, &
                       lowerflux/km3yr2m3sec

 tnflux = tnflux + netflux
 tuflux = tuflux + upperflux
 tlflux = tlflux + lowerflux


!    first determine a cross-section by listing the surface nodes
!    second save the u,v corresponding to those nodes and nodes3d corresponding
!    third find the distance between consecutive nodes
!    calculate the surface area by multiplying with the corresp. depth
!    find the angle of the surface to the horizontal
!    interpolate u,v to the center of the rotated surface
!    compute the flux by (u,v).AREA

close(122)

deallocate(fluxlon, fluxlat)
deallocate(dummy1, dummy2)
deallocate(snode)

enddo
if ( day2ext == 1 ) print*, "total distance: ", totdist

if ( day2ext == 1 ) print*, "----------------------------"


write(121,'(3F20.9)') tnflux/km3yr2m3sec/totalday, &
                      tuflux/km3yr2m3sec/totalday, &
                      tlflux/km3yr2m3sec/totalday
close(121)

 end subroutine compute_net_flux

 function marmara_volume()


integer            ::  elem, elnodes3(4)
real(r8)           ::  vol, lat_o, lon_o
real(r8)           ::  marmara_volume

vol = 0

do elem = 1,myDim_elem3d

elnodes3 = elem3D_nodes(:,elem)

 do k=1,4
   lat_o = coord_nod3D(2,elnodes3(k))/rad
   lon_o = coord_nod3D(1,elnodes3(k))/rad

   if ( lat_o.gt.marmax_lat.or.lon_o.gt.marmax_lon &
      .or.lat_o.lt.marmin_lat.or.lon_o.lt.marmin_lon) then
    go to 1234
   endif
 end do

 do k=1,4
    lat_o = coord_nod3D(2,elnodes3(k))/rad
    lon_o = coord_nod3D(1,elnodes3(k))/rad

    if ( lat_o.gt.41.and.lon_o.gt.28.7 ) go to 1234
 end do

  vol = vol + voltetra(elem)

   1234 continue
end do
  marmara_volume = vol



 end function marmara_volume

 subroutine bosphorus_for_blk_mfs

integer           :: hnode(4),layer
character(3)      :: DAYNUM
character(13)     :: OUTDIR
character(100)    :: OUTFILENAM,COORDINATES, LOCATION
integer           :: node
real(kind=8)      :: latmin, latmax, lonmin, lonmax, sssh

OUTDIR='BOSPHORUS_OBC'

day_loop: do day2ext = iniday, endday
write(DAYNUM,'(i3.3)') day2ext
OUTFILENAM=trim(OUTDIR)//'/'//trim(OUTDIR)//'_'//runid//'_'//runyear//'_'//DAYNUM//'.lst'
open(unit=101,file=OUTFILENAM,status='replace',access='append',form='formatted')
write(101,*)  'year: ', runyear
write(101,*)  'day  lon  lat  depth temp salt u v ssh'

call oce_input

latmin = 41.0
latmax = 41.4
lonmin = 28.8
lonmax = 29.6

node_loop: do node = 1,myDim_nod2D
if ( coord_nod2D(1,node)/rad .ge. lonmin .and. &
     coord_nod2D(1,node)/rad .le. lonmax ) then
  if ( coord_nod2D(2,node)/rad .ge. latmin .and. &
       coord_nod2D(2,node)/rad .le. latmax ) then

layer_loop: do layer = 1,110
     if ( layer == 1 ) then
             sssh = ssh(node)
     else
             sssh = -999.
     endif
     if ( nod3d_below_nod2d(layer,node).ne.-999 ) then
       write(101,'(I5,2F10.5,I7,4F8.3,F10.3)') day2ext, &
            coord_nod2d(1,node)/rad, coord_nod2d(2,node)/rad, layerdepth(layer), &
            tracer(nod3D_below_nod2D(layer,node),1), tracer(nod3D_below_nod2D(layer,node),2), &
            uf(nod3D_below_nod2D(layer,node)), uf(myDim_nod3d+nod3D_below_nod2D(layer,node)), &
            sssh
     endif
   end do layer_loop
  endif
endif
end do node_loop
close(101)
end do day_loop
 end subroutine bosphorus_for_blk_mfs

end module fesom_ocean_mod
