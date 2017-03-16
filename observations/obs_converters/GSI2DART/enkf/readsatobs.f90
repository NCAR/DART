module readsatobs
!$$$  module documentation block
!
! module: readsatobs                   read data from satellite radiance
!                                      diag* files.
!
! prgmmr: whitaker         org: esrl/psd               date: 2009-02-23
!
! abstract: read data from satellite radiance diag* files written out
!  by GSI forward operator code.
!
! Public Subroutines:
!  get_num_satobs: determine the number of observations to read.
!  get_satobs_data: read the data.
!   
! Public Variables: 
!
! Modules Used: read_diag
!
! program history log:
!   2009-02-23  Initial version.
!
! attributes:
!   language: f95
!
!$$$
  
use kinds, only: r_kind,i_kind,r_single
use read_diag, only: diag_data_fix_list,diag_header_fix_list,diag_header_chan_list, &
    diag_data_chan_list,diag_data_extra_list,read_radiag_data,read_radiag_header, &
    diag_data_name_list
use params, only: nsats_rad, nsatmax_rad, dsis, sattypes_rad

implicit none

private
public :: get_satobs_data, get_num_satobs

! just amsu, hirs, airs and ssmi data.
! for now allow 200 satellite data files, to increase, adjust 200 in next 3 lines

contains


subroutine get_num_satobs(obspath,datestring,num_obs_tot,id)
    use radinfo, only: iuse_rad,nusis,jpch_rad,nuchan,npred
    character (len=500), intent(in) :: obspath
    character(len=500) obsfile
    character(len=10), intent(in) :: id, datestring
    character(len=20) ::  sat_type
    integer(i_kind), intent(out) :: num_obs_tot
    integer(i_kind) iunit, iflag, nsat, ios,n,nkeep, i, jpchstart,indxsat
    integer(i_kind) npred_radiag
    logical fexist,lretrieval,lverbose
    real(r_kind) :: errorlimit,errorlimit2

    type(diag_header_fix_list )         :: header_fix0
    type(diag_header_chan_list),allocatable :: header_chan0(:)
    type(diag_data_fix_list   )         :: data_fix0
    type(diag_data_chan_list  ),allocatable :: data_chan0(:)
    type(diag_data_extra_list) ,allocatable :: data_extra0(:,:)
    type(diag_data_name_list)           :: data_name0

!  make consistent with screenobs
    errorlimit=1._r_kind/sqrt(1.e9_r_kind)
    errorlimit2=1._r_kind/sqrt(1.e-6_r_kind)
    iunit = 7
    lretrieval=.false.
    npred_radiag=npred
    lverbose=.false.

    num_obs_tot = 0
    do nsat=1,nsats_rad
        jpchstart=0
        do i=1,jpch_rad
          write(sat_type,'(a20)') adjustl(dsis(nsat))
          if(sat_type == trim(nusis(i)) .and. iuse_rad(i) > 0) then
            jpchstart=i
            exit
          end if
        end do
        if(jpchstart == 0) cycle
        nkeep = 0
        obsfile = trim(adjustl(obspath))//"diag_"//trim(sattypes_rad(nsat))//"_ges."//datestring//'_'//trim(adjustl(id))
        inquire(file=obsfile,exist=fexist)
        if (.not. fexist .or. datestring .eq. '0000000000') then
        obsfile = trim(adjustl(obspath))//"diag_"//trim(sattypes_rad(nsat))//"_ges."//trim(adjustl(id))
        endif

        inquire(file=obsfile,exist=fexist)
        if (.not.fexist) goto 900

        open(iunit,form="unformatted",file=obsfile,iostat=ios)
        rewind(iunit)
        call read_radiag_header(iunit,npred_radiag,lretrieval,header_fix0,header_chan0,data_name0,iflag,lverbose)

        do
           call read_radiag_data(iunit,header_fix0,lretrieval,data_fix0,data_chan0,data_extra0,iflag )
           if( iflag /= 0 )exit
           chan: do n=1,header_fix0%nchan
             if(header_chan0(n)%iuse<1) cycle chan
             indxsat=header_chan0(n)%iochan
             if(data_chan0(n)%qcmark < 0. .or. data_chan0(n)%errinv < errorlimit &
                      .or. data_chan0(n)%errinv > errorlimit2 &
                      .or. indxsat == 0) cycle chan
             if(data_extra0(1,n)%extra <= 0.001_r_kind .or.  &
                data_extra0(1,n)%extra > 1200._r_kind  .or. &
                abs(data_chan0(n)%tbobs) > 1.e9_r_kind) cycle chan
             nkeep = nkeep + 1
           end do chan
        enddo
        num_obs_tot = num_obs_tot + nkeep
900     continue
        close(iunit)
        write(6,100) nsat,trim(sattypes_rad(nsat)),nkeep,num_obs_tot
100     format(2x,i3,2x,a20,2x,'nkeep= ',i9,2x,'num_obs_tot= ',i9)
    enddo
end subroutine get_num_satobs

subroutine get_satobs_data(obspath, datestring, nobs_max, h_x, h_xnobc, x_obs, x_err, &
           x_lon, x_lat, x_press, x_time, x_channum, x_errorig, x_type, x_biaspred, x_indx,id,id2)
  use radinfo, only: iuse_rad,nusis,jpch_rad,nuchan,npred,adp_anglebc,emiss_bc
  character*500, intent(in) :: obspath
  character*500 obsfile,obsfile2
  character(len=10), intent(in) :: id,id2

  real(r_single), dimension(nobs_max) :: h_x,h_xnobc,x_obs,x_err,x_lon,&
                               x_lat,x_press,x_time,x_errorig
  real(r_single), dimension(npred+1,nobs_max) :: x_biaspred
  integer(i_kind), dimension(nobs_max) ::  x_channum,x_indx
  character(len=20), dimension(nobs_max) ::  x_type
  character(len=20) ::  sat_type
  character(len=10), intent(in) ::  datestring

  integer(i_kind) nobs_max, iunit, iunit2,iflag, nobs, n, nsat, i,jpchstart,indxsat
  integer(i_kind) npred_radiag,iflag2
  logical twofiles,fexist1,fexist2,lretrieval,lverbose
  real(r_kind) :: errorlimit,errorlimit2

  type(diag_header_fix_list )         :: header_fix1,header_fix2
  type(diag_header_chan_list),allocatable :: header_chan1(:),header_chan2(:)
  type(diag_data_fix_list   )         :: data_fix1,data_fix2
  type(diag_data_chan_list  ),allocatable :: data_chan1(:),data_chan2(:)
  type(diag_data_extra_list) ,allocatable :: data_extra1(:,:),data_extra2(:,:)
  type(diag_data_name_list)           :: data_name1,data_name2

! make consistent with screenobs
  errorlimit=1._r_kind/sqrt(1.e9_r_kind)
  errorlimit2=1._r_kind/sqrt(1.e-6_r_kind)

  iunit = 7
  iunit2 = 17
  lretrieval=.false.
  npred_radiag=npred
  lverbose=.false.

  nobs = 0
  twofiles = id /= id2

  do nsat=1,nsats_rad
     jpchstart=0
     do i=1,jpch_rad
       write(sat_type,'(a20)') adjustl(dsis(nsat))
       if(sat_type == trim(nusis(i)) .and. iuse_rad(i) > 0) then
         jpchstart = i
         exit
       end if
     end do
     if(jpchstart == 0) cycle
     obsfile = trim(adjustl(obspath))//"diag_"//trim(sattypes_rad(nsat))//"_ges."//datestring//'_'//trim(adjustl(id))
     inquire(file=obsfile,exist=fexist1)
     if (.not. fexist1 .or. datestring .eq. '0000000000') then
     obsfile = trim(adjustl(obspath))//"diag_"//trim(sattypes_rad(nsat))//"_ges."//trim(adjustl(id))
     endif
     inquire(file=obsfile,exist=fexist1)
     if(.not.fexist1) goto 900

     open(iunit,form="unformatted",file=obsfile)
     rewind(iunit)
     call read_radiag_header(iunit,npred_radiag,lretrieval,header_fix1,header_chan1,data_name1,iflag,lverbose)

     if(twofiles)then
       obsfile2=trim(adjustl(obspath))//"diag_"//trim(sattypes_rad(nsat))//"_ges."//datestring//'_'//trim(adjustl(id2))
       inquire(file=obsfile2,exist=fexist2)
       if (.not. fexist2 .or. datestring .eq. '0000000000') then
       obsfile2 = trim(adjustl(obspath))//"diag_"//trim(sattypes_rad(nsat))//"_ges."//trim(adjustl(id2))
       endif
       inquire(file=obsfile2,exist=fexist2)
       if(.not.fexist2) goto 900

       open(iunit2,form="unformatted",file=obsfile2)
       rewind(iunit2)
       call read_radiag_header(iunit2,npred_radiag,lretrieval,header_fix2,header_chan2,data_name2,iflag2,lverbose)
     end if

     do
      call read_radiag_data(iunit,header_fix1,lretrieval,data_fix1,data_chan1,data_extra1,iflag )
      if( iflag /= 0 ) then
       exit
      end if
      if(twofiles)then
         call read_radiag_data(iunit2,header_fix2,lretrieval,data_fix2,data_chan2,data_extra2,iflag2 )
        if( header_fix1%nchan /= header_fix2%nchan .or. abs(data_fix1%lat-data_fix2%lat) .gt. 1.e-5 .or.  &
            abs(data_fix1%lon-data_fix2%lon) .gt. 1.e-5 .or. abs(data_fix1%obstime-data_fix2%obstime) .gt. 1.e-5) then
           write(6,*) 'inconsistent files',trim(obsfile2)
           write(6,*) 'nchan',header_fix1%nchan,header_fix2%nchan
           write(6,*) 'lat',data_fix1%lat,data_fix2%lat
           write(6,*) 'lon',data_fix1%lon,data_fix2%lon
           write(6,*) 'obstim',data_fix1%obstime,data_fix2%obstime
           call stop2(-99)
        end if
      end if
      chan:do n=1,header_fix1%nchan
         if(header_chan1(n)%iuse<1) cycle chan
         indxsat=header_chan1(n)%iochan
         if(data_chan1(n)%qcmark < 0. .or. data_chan1(n)%errinv < errorlimit &
                  .or. data_chan1(n)%errinv > errorlimit2 &
                  .or. indxsat == 0) cycle chan
         if(data_extra1(1,n)%extra <= 0.001_r_kind .or.  &
            data_extra1(1,n)%extra > 1200._r_kind  .or. &
            abs(data_chan1(n)%tbobs) > 1.e9_r_kind) cycle chan
         nobs = nobs + 1 
         if (nobs > nobs_max) then
             print *,'warning:  exceeding array bounds in readinfo_from_file',&
             nobs,nobs_max
         end if
         x_type(nobs)= sat_type
         x_channum(nobs) = n
         x_indx(nobs) = indxsat
         x_lon(nobs) = data_fix1%lon
         x_lat(nobs) = data_fix1%lat
         x_time(nobs) = data_fix1%obstime
         x_obs(nobs) = data_chan1(n)%tbobs 
         ! bias corrected Hx
         h_x(nobs) = x_obs(nobs) - data_chan1(n)%omgbc 
         ! un-bias corrected Hx
         if(twofiles)then
           h_xnobc(nobs) = x_obs(nobs) - data_chan2(n)%omgnbc
         else
           h_xnobc(nobs) = x_obs(nobs) - data_chan1(n)%omgnbc
         end if
         ! data_chan%errinv is inverse error variance.
         x_errorig(nobs) = header_chan1(n)%varch**2
         x_err(nobs) = (1._r_kind/data_chan1(n)%errinv)**2
         x_press(nobs) = data_extra1(1,n)%extra

!! DTK:  **NOTE**
!!       The bifix term will need to be expanded if/when the GSI/GDAS goes to using
!!       a higher polynomial version of the angle dependent bias correction (if
!!       and when it is moved into part of the varbc)
!!         x_biaspred(1,nobs) = data_chan1(n)%bifix! fixed angle dependent bias
         x_biaspred(1,nobs) = data_chan1(n)%bifix(1) ! fixed angle dependent bias
         x_biaspred(2,nobs) = data_chan1(n)%bicons ! constant bias correction
         x_biaspred(3,nobs) = data_chan1(n)%biang ! scan angle bias correction
         x_biaspred(4,nobs) = data_chan1(n)%biclw ! CLW bias correction
         x_biaspred(5,nobs) = data_chan1(n)%bilap2 ! square lapse rate bias corr
         x_biaspred(6,nobs) = data_chan1(n)%bilap ! lapse rate bias correction
         if (npred == 7) then
           x_biaspred(7,nobs) = data_chan1(n)%bicos ! node*cos(lat) bias correction for SSMIS
           x_biaspred(8,nobs) = data_chan1(n)%bisin ! sin(lat) bias correction for SSMIS                    
         endif
         if (emiss_bc) x_biaspred(9,nobs) = data_chan1(n)%biemis

         if (adp_anglebc) then
            x_biaspred( 1,nobs)  = data_chan1(n)%bifix(5) ! fixed angle dependent bias correction
            x_biaspred(npred-2,nobs)  = data_chan1(n)%bifix(1) ! 4th order scan angle (predictor)
            x_biaspred(npred-1,nobs)  = data_chan1(n)%bifix(2) ! 3rd order scan angle (predictor)
            x_biaspred(npred,nobs)  = data_chan1(n)%bifix(3) ! 2nd order scan angle (predictor)
            x_biaspred(npred+1,nobs)    = data_chan1(n)%bifix(4) ! 1st order scan angle (predictor)
         endif

      enddo chan
     enddo

     cycle

900  continue
     close(iunit)
     if(twofiles)close(iunit2)
 enddo

 end subroutine get_satobs_data

end module readsatobs
