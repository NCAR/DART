!
! Includes tools to process forcing files
! Tools are used to produce the analysis in:
! Aydoğdu, A., Pinardi, N., Özsoy, E., Danabasoglu, G., Gürses, Ö., and Karspeck, A.: Circulation of the Turkish
! Straits System under interannual atmospheric forcing, Ocean Sci., 14, 999-1019, https://doi.org/10.5194/os-14-999-2018, 2018b.
!
! Provided by ali.aydogdu@cmcc.it
!
module fesom_observation_mod

use g_config,   only : runid, day2ext, runyear
use o_param,    only : rad
use o_mesh,     only : coord_nod2d, nod3d_below_nod2d, layerdepth, &
                       max_num_layers
use o_array,    only : tracer
use g_parfe,    only : mydim_nod2d

use utilities, only : r8, i4

implicit none

public    :: read_obs_nml,             & ! read observation namelist nml/namelist.obs
             read_ctd_data,            & ! compare model outputs with ctd profiles
             read_ship_track,          & ! compare model outputs with ship tracks
             calc_misfit,              & ! calculate misfits with CTD profiles
             calc_rms_latlon,          & ! calculate RMSE against CTD
             find_close_2D,            & ! find closest 2D location on the mesh to the obs
             find_close_layer,         & ! find closest level on the mesh to the obs
             profile_from_netcdf,      & ! extract profile from fesom ocean outputs
             synthetic_ferrybox_from_nr  ! generate ferrybox data from nature run (if exists, a ship track is needed)

character(len=100)   :: ctddir,wrkdir
character(len=100)   :: obslist, ctd_input_file
character(len=100)   :: fb_outdir, ferrylist
character(len=100)   :: synt_proutdir, synt_profilst
character(len=100)   :: write_tk_file, w_tk_format
character(len=100)   :: write_ms_file, w_ms_format
character(len=100)   :: write_vr_file, w_vr_format
character(len=100)   :: ship_dir
character(len=6)     :: surveyno
character(len=8)     :: hour
character(len=4)     :: year
character(len=3)     :: ship_mnth
character(len=4)     :: ship_year

integer              :: surveyyr
integer              :: debug,nobs
integer              :: w_tk_unit, w_ms_unit, w_vr_unit

integer              :: oyear,omnt,oday,mday
real(r8)               :: clon,clat,mlon,mlat,obslon,obslat

real(r8)               :: odep, idep, fdep
real(r8)               :: pres, temp, sigt, salt
real(r8)               :: mtem, msal
real(r8)               :: rms_tem, rms_sal
real(r8)               :: stmean, stdev_tem, stdev_sal

real(r8)               :: nodist, distno
integer              :: i,j,close_node, lev, upper, lower, c_lev
integer              :: ulayer, blayer
integer              :: var_index

integer,parameter    :: R=6373
integer,parameter    :: MISSING_DEPTH=-999
real(r8),parameter     :: MISSING_COORD=-999.0
integer, parameter   :: nvar=2
integer, parameter   :: nlev=50
real(r8)               :: bounds(4)
real(r8)               :: profile(nvar*3,nlev)
integer              :: depexist(nlev)
integer              :: fb_layer


contains

!----------------------------------------------------------
! subroutine to read the namelist.obs
!----------------------------------------------------------
  subroutine read_obs_nml

  character(len=100)   :: nmlfile
  namelist /debugging/    debug
  namelist /ctd_obs/      ctddir, wrkdir, surveyyr, surveyno, obslist
  namelist /fb_obs/       fb_layer, fb_outdir, ferrylist
  namelist /synt_profile/ synt_proutdir, synt_profilst, ulayer, blayer
  namelist /perturbation/ stmean, stdev_tem, stdev_sal
  namelist /region/       bounds
  namelist /ship_track/   ship_mnth, ship_year, ship_dir

  nmlfile ='namelist.obs'    ! name of general observation  namelist file
  open (20,file=nmlfile)
  read (20,NML=debugging)
  read (20,NML=ctd_obs)
  read (20,NML=fb_obs)
  read (20,NML=synt_profile)
  read (20,NML=perturbation)
  read (20,NML=region)
  read (20,NML=ship_track)
  close(20)


    if (debug > 2) then
            print*, "bounds: ",bounds(1), bounds
    endif

  end subroutine read_obs_nml

!----------------------------------------------------------
! subroutine to read the ctd data and find closest model node
!----------------------------------------------------------
  subroutine read_ctd_data

  integer              :: r_lst_unit,r_tk_unit,rcio
  integer              :: i

  call read_obs_nml
  rms_tem=0
  rms_sal=0
  profile=0
  depexist=0
  write(year,'(i4)') surveyyr
  r_lst_unit=13
  r_tk_unit=14
  w_tk_unit=15
  w_ms_unit=16
  w_vr_unit=17
  open(file=obslist, unit=r_lst_unit)

  write_ms_file=runid//'_ms_'//trim(surveyno)//'.asc'
  open(unit=w_ms_unit,file=write_ms_file,status='replace',access='append',form='formatted')
  obsFile_loop: do
    read(r_lst_unit, "(A)", iostat=rcio) ctd_input_file
    if (rcio /= 0) then
      if (debug > 2) print *, 'no other files, rcio = ', rcio
        exit obsFile_loop
    endif

    open(file=ctd_input_file,unit=r_tk_unit)
    read(r_tk_unit,*,iostat=rcio) clon,clat
      if ( clon == MISSING_COORD .OR. clat == MISSING_COORD ) CYCLE obsFile_loop
      if ( clon < bounds(1) .OR. clon > bounds(2) ) CYCLE obsFile_loop
      if ( clat < bounds(3) .OR. clat > bounds(4) ) CYCLE obsFile_loop
    read(r_tk_unit,*,iostat=rcio) oyear, omnt, oday, mday, hour
    read(r_tk_unit,*,iostat=rcio) idep, fdep

    if (debug > 4) then
      print*, ctd_input_file
      print*, 'lon: ',clon,'lat: ',clat
      print*, 'oday: ',oday,' mnt: ',omnt,' year: ',oyear,' hour: ',hour,' mday: ',mday
    endif

    day2ext=mday
    call oce_input
    call find_close_2D(clon,clat)

    write_tk_file=runid//'_'//trim(ctd_input_file)
    open(unit=w_tk_unit,file=write_tk_file,status='replace',access='append',form='formatted')

    nobs=0
    obs_loop: do i=INT(idep),INT(fdep)
      nobs=nobs+1
      read(r_tk_unit,*, iostat=rcio) odep, pres, temp, sigt, salt
       if (rcio /= 0) then
        if (debug > 2 ) print *, 'no other data, rcio = ', rcio
        exit obs_loop
       endif
       call find_close_layer(odep)
       if ( nod3D_below_nod2D(c_lev,close_node) == MISSING_DEPTH ) EXIT obs_loop
      if (debug > 5) then
        print*, 'idep: ',idep,' fdep: ',fdep
        print*, 'close_node: ',close_node,' c_lev: ',c_lev
        print*, 'VALUES:  ', odep, pres, temp, sigt, salt
      endif
        mtem=tracer((nod3D_below_nod2D(c_lev,close_node)),1)
        msal=tracer((nod3D_below_nod2D(c_lev,close_node)),2)
      w_tk_format="(3I5,1X,F12.5,1X,I10,1X,8F12.5)"
        write(w_tk_unit, w_tk_format) oyear, omnt, oday, odep, &
        nod3D_below_nod2D(c_lev,close_node), &
        clon, clat, coord_nod2d(1,close_node)/rad, coord_nod2d(2,close_node)/rad, &
        temp, mtem, salt, msal
        call calc_misfit

        if ( odep >= idep .AND. odep < fdep .AND. i <= nlev) then
                depexist(i)=depexist(i)+1
                profile(1,i) = profile(1,i) + (temp-mtem)**2
                profile(2,i) = profile(2,i) + (salt-msal)**2
                profile(3,i) = profile(3,i) + temp
                profile(4,i) = profile(4,i) + salt
                profile(5,i) = profile(5,i) + mtem
                profile(6,i) = profile(6,i) + msal
        endif
    end do obs_loop
    call calc_rms_latlon
  close(w_tk_unit)
  close (r_tk_unit)
  end do obsFile_loop
  close (w_ms_unit)
    print*, depexist
  write_vr_file=runid//'_vr_'//trim(surveyno)//'.asc'
  w_vr_format="(6F12.5)"
  open(unit=w_vr_unit,file=write_vr_file,status='replace',access='append',form='formatted')
  do i=1,nlev
    do var_index=1,6
      profile(var_index,i)=profile(var_index,i)/depexist(i)
    end do
        write(w_vr_unit, w_vr_format) sqrt(profile(1,i)), sqrt(profile(2,i)),  &
             profile(3,i), profile(4,i), profile(5,i), profile(6,i)
  end do
  close (w_vr_unit)
  close (r_lst_unit)
  end subroutine read_ctd_data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_ship_track

    character(len=100)   :: ship_filename, write_ship_file
    integer              :: r_tk_unit, w_tk_unit, rcio
    real(r8)               :: ltime1, ltime2, slat, slon, stime,sdep, stem
    integer              :: sday, smon, syr, mday, m_lev

    r_tk_unit = 1001
    w_tk_unit = 1002

    m_lev = 1

     call read_obs_nml

     ship_filename = trim(ship_dir)//'/marmara_'//ship_mnth//ship_year//'.txt'

     open(file=ship_filename,unit=r_tk_unit)

     w_tk_format="(3I5,1X,F12.5,1X,I10,1X,8F12.5)"
     write_ship_file='marmara_'//ship_mnth//ship_year//'.asc'
     open(unit=w_tk_unit,file=write_ship_file,status='replace',access='append',form='formatted')

     obs_loop: do
       read(r_tk_unit,*, iostat=rcio) ltime1, slat, slon, stime, syr, smon, sday, ltime2, sdep, stem, mday
       if (rcio /= 0) then
         if (debug > 2) print *, 'no other files, rcio = ', rcio
           exit obs_loop
       endif

      day2ext=mday
      call oce_input
      call find_close_2D(slon,slat)

       mtem=tracer((nod3D_below_nod2D(m_lev,close_node)),1)
       write(w_tk_unit, w_tk_format)  syr, smon, sday, sdep, nod3D_below_nod2D(m_lev,close_node),  &
                                      slon, slat, coord_nod2d(1,close_node)/rad, coord_nod2d(2,close_node)/rad, &
                                      stem, mtem
      end do obs_loop

      close(w_tk_unit)
      close(r_tk_unit)

  end subroutine read_ship_track
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine calc_misfit

     integer,parameter    :: ndep=2
     integer        :: k
     real(r8)         :: tem_ms, sal_ms
     real(r8)         :: dplv(ndep)

  data dplv /0, 50/

  tem_ms=0
  sal_ms=0

  do k=2,ndep
  if ( odep > dplv(k-1) .AND. odep <= dplv(k) ) then
          tem_ms=temp-mtem
          sal_ms=salt-msal
          rms_tem=rms_tem+(tem_ms)**2
          rms_sal=rms_sal+(sal_ms)**2
  endif
  end do
  end subroutine calc_misfit

  subroutine calc_rms_latlon

      rms_tem=sqrt(rms_tem/nobs)
      rms_sal=sqrt(rms_sal/nobs)
      w_ms_format="(3I5,1X,6F12.5)"
      write(w_ms_unit, w_ms_format) oyear, omnt, oday, &
      clon, clat, coord_nod2d(1,close_node)/rad, coord_nod2d(2,close_node)/rad, &
      rms_tem, rms_sal

  end subroutine calc_rms_latlon

  subroutine find_close_2D(obslon, obslat)

   real(r8),intent(in)      :: obslon, obslat
   real(r8)                 :: mindist, distance
    mindist=999999999.0
  node_loop:   do j=1,myDim_nod2D
       distance=((obslon*rad-coord_nod2d(1,j))**2 + (obslat*rad-coord_nod2d(2,j))**2)**0.5
       if ( distance < mindist ) then
        if (debug > 5) then
         print*, 'distance: ', j, distance, coord_nod2d(1,j)/rad, coord_nod2d(2,j)/rad
        endif
         mindist=distance
         close_node=INT(j)
       endif
    end do node_loop
  return
  end subroutine find_close_2D

  subroutine find_close_layer(odepth)

    real(r8)     ::  odepth

  do lev=1,max_num_layers-1
    if ( odepth > layerdepth(lev) .AND. odepth <= layerdepth(lev+1) ) then
     lower=lev+1
     upper=lev
     if ( (layerdepth(upper) - odepth) < (layerdepth(lower) - odepth) ) then
      c_lev=upper
     else
      c_lev=lower
     endif
    endif
  enddo
  return
  end subroutine find_close_layer
!----------------------------------------------------------
!-- CREATE SYNTHETIC PROFILE OBSERVATIONS FOR OSSE --------
!----------------------------------------------------------
  subroutine profile_from_netcdf
    use random_perturbation

    integer           :: check_EOF
    integer           :: fileID, k, hnode, fileno
    integer           :: min_begin, yr, mnt, day, hr, minute, second
    integer           :: prolayer
    real(r8)            :: lon, lat
    character*6       :: DAYNUM
    character*100     :: OUTFILENAM, profilename
    real(r8)            :: TEM, SAL

  call read_obs_nml
    OUTFILENAM = trim(fb_outdir)//'/Profile_'//runid//'_'//runyear//'.asc'

    open(unit=101,file=OUTFILENAM,status='replace',access='append',form='formatted')
      write(101,'(A)') 'ob_No  iprof  ilev   ob_kind  ob_value   ob_lat     ob_lon  ob_depth  iyy imm idd ihh imin iss   QC   ob_unit   ob_type   file'
    open(unit=103,file=synt_profilst)
    do
    print*, synt_profilst
       read(103,*,iostat=check_EOF) profilename
       print*, check_EOF
       if ( check_EOF /= 0 )  EXIT
       print*, trim(profilename)
       open(unit=102,file=trim(profilename))
       do
       read(102,*,iostat=check_EOF) lon, lat, min_begin, fileno, &
                                     yr, mnt, day, hr, minute, second
                                     print*, yr, second
       if ( check_EOF /=  0 )  EXIT
       call find_close_2D(lon,lat)
       hnode   = close_node
       day2ext = ( day - 1 ) * 24 + hr
        call oce_input
        do prolayer=ulayer,blayer
         if ( nod3D_below_nod2D(prolayer,hnode).eq.-999 ) EXIT
           TEM = tracer((nod3D_below_nod2D(prolayer,hnode)),1)
           write(101,'(A,3F10.5,2I5,5I4.2,A,1X,A)') "1  1  1  TT  ",TEM,                      &
                                                   lat,lon,layerdepth(prolayer),             &
                                                    yr, mnt, day, hr, minute, second,        &
                                                    " 1 degC Validation", trim(profilename)
           SAL = tracer((nod3D_below_nod2D(prolayer,hnode)),2)
           write(101,'(A,3F10.5,2I5,5I4.2,A,1X,A)') "1  1  1  SS  ",SAL,                      &
                                                   lat,lon,layerdepth(prolayer),             &
                                                    yr, mnt, day, hr, minute, second,        &
                                                    " 1  psu Validation", trim(profilename)
        end do
    end do
    close(102)
    end do
    close(103)
    close(101)
  end subroutine profile_from_netcdf


!----------------------------------------------------------
!-- CREATE SYNTHETIC FERRYBOX OBSERVATIONS FOR OSSE -------
!----------------------------------------------------------
  subroutine synthetic_ferrybox_from_nr
    use random_perturbation

    integer           :: check_EOF
    integer           :: fileID,k,hnode,fileno
    integer           :: min_begin, yr, mnt, day, hr, minute, second
    real(r8)          :: lon,lat
    character(len=6)  :: DAYNUM
    character(len=100):: OUTFILENAM,OBSCREAFILE,ferry_filename

    real(r8)            :: TEM, SAL

  call read_obs_nml
    OUTFILENAM = trim(fb_outdir)//'/Ferrybox_'//runid//'_'//runyear//'.asc'

    open(unit=101,file=OUTFILENAM,status='replace',access='append',form='formatted')
      write(101,'(A)') 'ob_No  iprof  ilev   ob_kind  ob_value   ob_lat     ob_lon ob_depth  iyy imm idd ihh imin iss   QC   ob_unit   ob_type   file'
    open(unit=103,file=ferrylist)
    do
    print*, ferrylist
       read(103,*,iostat=check_EOF) ferry_filename
       print*, check_EOF
       if ( check_EOF /= 0 )  EXIT
       print*, trim(ferry_filename)
       open(unit=102,file=trim(ferry_filename))
       do
       read(102,*,iostat=check_EOF) lon, lat, min_begin, fileno, &
                                     yr, mnt, day, hr, minute, second
                                     print*, yr, second
       if ( check_EOF /=  0 )  EXIT
       call find_close_2D(lon,lat)
       hnode=close_node
       day2ext = ( day - 1 ) * 24 + hr
        call oce_input
         if ( nod3D_below_nod2D(fb_layer,hnode).eq.-999 ) EXIT
           TEM = tracer((nod3D_below_nod2D(fb_layer,hnode)),1) + &
                 random_around_mean(stmean,stdev_tem)
           write(101,'(A,3F10.5,2I5,5I4.2,A,1X,A)') "1  1  1  FT  ",TEM,                      &
                                                   lat,lon,layerdepth(fb_layer),             &
                                                    yr, mnt, day, hr, minute, second,        &
                                                    " 1 degC Ferrybox", trim(ferry_filename)
           SAL = tracer((nod3D_below_nod2D(fb_layer,hnode)),2) +  &
                 random_around_mean(stmean,stdev_sal)
           write(101,'(A,3F10.5,2I5,5I4.2,A,1X,A)') "1  1  1  FS  ",SAL,                      &
                                                   lat,lon,layerdepth(fb_layer),             &
                                                    yr, mnt, day, hr, minute, second,        &
                                                    " 1  psu Ferrybox", trim(ferry_filename)
    end do
    close(102)
    end do
    close(103)
    close(101)
  end subroutine synthetic_ferrybox_from_nr

end module fesom_observation_mod
