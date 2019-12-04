! This code is not protected by the DART copyright agreement.
! DART $Id$

!     Nick Pedatella
!     Program to convert Aura Temperature data to DART Observation
!     sequence files. 

program convert_aura

use        types_mod, only : r8
use time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, set_time,&
                             increment_time, get_time, set_date, operator(-),  &
                             print_date, operator(+)
use    utilities_mod, only : initialize_utilities, find_namelist_in_file,    &
                             check_namelist_read, nmlfileunit, do_nml_file,    &
                             get_next_filename, error_handler, E_ERR, E_MSG, &
                             find_textfile_dims, finalize_utilities, &
                             timestamp,do_nml_term
use  netcdf_utilities_mod, only : nc_check
use     location_mod, only : VERTISPRESSURE, set_location ! pressure ccordinates for SABER
use obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq,       &
                             static_init_obs_sequence, init_obs, destroy_obs, &
                             write_obs_seq, init_obs_sequence, get_num_obs,   &
                             insert_obs_in_seq, destroy_obs_sequence,         &
                             set_copy_meta_data, set_qc_meta_data, set_qc,    & 
                             set_obs_values, set_obs_def, insert_obs_in_seq
use obs_def_mod,      only : obs_def_type, set_obs_def_time, set_obs_def_type_of_obs, &
                             set_obs_def_error_variance, set_obs_def_location, &
                             set_obs_def_key
use obs_utilities_mod, only : create_3d_obs,add_obs_to_seq
use     obs_kind_mod, only : AURAMLS_TEMPERATURE

use           netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

integer, parameter :: num_copies = 1,   &   ! number of copies in sequence
                      num_qc     = 1        ! number of QC entries

! variables
character (len=130) :: msgstring, next_infile
character (len=80) :: name

integer :: iyear,iday,imonth,idoy,ihr,imin,isec,nalt,nevent,ialt,ievent, &
           nfiles,ncid,varid,filenum,io,iunit,obs_num,num_new_obs,k,osec,oday
logical :: file_exist, first_obs, did_obs, from_list = .false.

real(r8) :: qc,oerr,pres_out

real(r8), allocatable :: lat(:), lon(:), pres(:),temp(:,:),qual(:),conv(:),time(:)
real(r8), allocatable :: localsolartime(:)

type(obs_def_type)      :: obs_def
type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: time_obs, time_anal, prev_time


! namelist parameters
character(len=128) :: aura_netcdf_file     = 'aura_input.nc'
character(len=128) :: aura_netcdf_filelist = 'aura_input_list'
character(len=128) :: aura_outfile        = 'obs_seq.aura'
integer :: aura_yr=2011,aura_doy = 1

namelist /convert_aura_nml/ aura_netcdf_file, aura_netcdf_filelist, &
                             aura_outfile,aura_yr,aura_doy

! initialize some values
obs_num = 1
qc = 1.0_r8
first_obs = .true.
call set_calendar_type(GREGORIAN)

!  read the necessary parameters from input.nml
call initialize_utilities()

call find_namelist_in_file("input.nml", "convert_aura_nml", iunit)
read(iunit, nml = convert_aura_nml, iostat = io)
call check_namelist_read(iunit, io, "convert_aura_nml")

! Record the namelist values used for the run 
if (do_nml_file()) write(nmlfileunit, nml=convert_aura_nml)
if (do_nml_term()) write(     *     , nml=convert_aura_nml)

! namelist checks for sanity

! cannot have both a single filename and a list; the namelist must
! shut one off.
if (aura_netcdf_file /= '' .and. aura_netcdf_filelist /= '') then
  call error_handler(E_ERR, 'convert_aura_cdf',                     &
                     'One of aura_netcdf_file or filelist must be NULL', &
                     source, revision, revdate)
endif
if (aura_netcdf_filelist /= '') from_list = .true.

! Define Number of observations:
if (from_list) then
  num_new_obs = (400 * 40000) * nfiles
else
  num_new_obs = (400 * 40000 )
endif

! output date
print *, 'OUTPUTING FOR YR=',aura_yr,'DOY=',aura_doy

!  either read existing obs_seq or create a new one
call static_init_obs_sequence()
call init_obs(obs, num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)
inquire(file=aura_outfile, exist=file_exist)
if ( file_exist ) then

   print *, "found existing obs_seq file, overwriting ", trim(aura_outfile)
  print *, "max entries = ", num_new_obs
  call init_obs_sequence(obs_seq, num_copies, num_qc, num_new_obs)

  do k = 1, num_copies
    call set_copy_meta_data(obs_seq, k, 'observation')
  end do
  do k = 1, num_qc
    call set_qc_meta_data(obs_seq, k, 'Data QC')
  end do

else

  print *, "no existing obs_seq file, creating ", trim(aura_outfile)
  print *, "max entries = ", num_new_obs
  call init_obs_sequence(obs_seq, num_copies, num_qc, num_new_obs)

  do k = 1, num_copies
    call set_copy_meta_data(obs_seq, k, 'observation')
  end do
  do k = 1, num_qc
    call set_qc_meta_data(obs_seq, k, 'Data QC')
  end do

end if


did_obs = .false.
! main loop that does either a single file or a list of files

filenum = 1

fileloop: do  ! until out of files

   ! get the single name, or the next name from a list
   if (from_list) then 
      next_infile = get_next_filename(aura_netcdf_filelist, filenum)
   else
      next_infile = aura_netcdf_file
      if (filenum > 1) next_infile = ''
   endif
   if (next_infile == '') exit fileloop

   ! open the file
   call nc_check( nf90_open(next_infile, nf90_nowrite, ncid), 'file open', next_infile)

   ! get the number of events and altitude:
   call nc_check( nf90_inq_dimid(ncid,"levels",varid), 'inq dimid altitude')
   call nc_check( nf90_inquire_dimension(ncid,varid,name,nalt), 'inq dim altitude')
   
   call nc_check( nf90_inq_dimid(ncid,"times",varid), 'inq dimid event')
   call nc_check( nf90_inquire_dimension(ncid,varid,name,nevent), 'inq dim event')
   print *, 'nalt=',nalt
   allocate( time(nevent) )
   allocate( lat(nevent) )
   allocate( lon(nevent) )
   allocate( pres(nalt) )
   allocate( temp(nalt,nevent) )
   allocate( qual(nevent) )
   allocate( conv(nevent) )
   allocate( localsolartime(nevent) )

   ! read in the latitude array
   call nc_check( nf90_inq_varid(ncid, "Latitude", varid) ,'inq varid Lat')
   call nc_check( nf90_get_var(ncid, varid, lat)     ,'get var   Lat')

   ! read the longitude array
   call nc_check( nf90_inq_varid(ncid, "Longitude", varid) ,'inq varid Lon')
   call nc_check( nf90_get_var(ncid, varid, lon)     ,'get var   Lon')

   ! read the altitude array
   call nc_check( nf90_inq_varid(ncid, "Pressure", varid) ,'inq varid pressure')
   call nc_check( nf90_get_var(ncid, varid, pres)        ,'get_var   pressure')

   ! read the temperature
   call nc_check( nf90_inq_varid(ncid, "Temperature", varid) ,'inq varid ktemp')
   call nc_check( nf90_get_var(ncid, varid, temp)        ,'get_var   ktemp')
   
   ! quality
   call nc_check( nf90_inq_varid(ncid, "Quality", varid) ,'inq varid quality')
   call nc_check( nf90_get_var(ncid, varid, qual)        ,'get_var   quality')
   ! convergence
   call nc_check( nf90_inq_varid(ncid, "Convergence", varid) ,'inq varid conv')
   call nc_check( nf90_get_var(ncid, varid, conv)        ,'get_var   conv')
   ! LST
   call nc_check( nf90_inq_varid(ncid, "LocalSolarTime", varid) ,'inq varid lst')
   call nc_check( nf90_get_var(ncid, varid, localsolartime)        ,'get_var   lst')

   ! read the date/time
   call nc_check( nf90_inq_varid(ncid, "Time", varid) ,'inq varid time')
   call nc_check( nf90_get_var(ncid, varid, time)        ,'get_var   time')
   
   ! finished reading the data 
   call nc_check( nf90_close(ncid), 'close file')

   ! now loop through each event / altitude and add to obs sequence file
   do ievent=1,nevent

     do ialt=1,nalt

       !calculate seconds from LST
       isec = int( (localsolartime(ievent)-(lon(ievent)/15.))/24.*86400. )
       if (isec.le.0) isec = isec + 86400       
       if (isec.ge.86400) isec = isec - 86400       

       iyear = aura_yr ! yr & day are inputs
       idoy = aura_doy 
       call convert_day(iyear,idoy,imonth,iday)

       if (temp(ialt,ievent).ne.-999. .and. & ! data is missing
           qual(ievent).ge.0.65 .and. & !specified quality check
           conv(ievent).le.1.2 .and. & ! specified convergence
           pres(ialt).ge.0.001 .and. pres(ialt).le.260 .and. & !good altitude range
           isec .lt.86400 ) then 
           

         ihr = int( isec / 3600 )
         imin = int( (isec-ihr*3600) / 60 )
         isec = int( (isec-ihr*3600)-imin*60 )
         time_obs = set_date(iyear,imonth,iday,ihr,imin,isec)
         call get_time(time_obs,osec,oday)

         ! need to set the error (for now just take SABER error):
         call saber_error(pres(ialt),oerr)

         ! convert pressure from mbar(hPa) -> Pa
         pres_out = pres(ialt)*100.

         if (lon(ievent).le.0) lon(ievent) = lon(ievent)+360.

         ! now creat the observation:
         call create_3d_obs(lat(ievent),lon(ievent),pres_out, &
                            VERTISPRESSURE,temp(ialt,ievent), &
                          AURAMLS_TEMPERATURE,oerr,oday,osec,qc,obs)
         call add_obs_to_seq(obs_seq,obs,time_obs,prev_obs,prev_time,first_obs)

         if(.not. did_obs) did_obs = .true.
       else

       end if
     end do ! end alt loop
   end do ! end event loop

   ! clean up and loop if there is another input file
   deallocate( lat,lon,pres,qual,conv,time,temp)
 
   filenum = filenum + 1

end do fileloop

! done with main loop. if we added any obs to the sequence, write it out:
if (did_obs) then
!print *, 'ready to write, nobs = ', get_num_obs(obs_seq)
   if (get_num_obs(obs_seq) > 0) &
      call write_obs_seq(obs_seq, aura_outfile)

   ! minor stab at cleanup, in the off chance this will someday get turned
   ! into a subroutine in a module.  probably not all that needs to be done,
   ! but a start.
   call destroy_obs(obs)
! if obs == prev_obs then you can't delete the same obs twice.
! buf if they differ, then it's a leak.  for now, don't delete prev
! since the program is exiting here anyway.
   !call destroy_obs(prev_obs)
   if (get_num_obs(obs_seq) > 0) call destroy_obs_sequence(obs_seq)
endif

! END OF MAIN ROUTINE

contains

subroutine convert_day(iyear,idoy,imonth,iday)
! NMP - Subroutine to convert from yr,doy to yr,mon,day


 integer, intent(in)   :: iyear
 integer, intent(in)   :: idoy
 integer, intent(out)  :: imonth
 integer, intent(out)  :: iday

 integer :: t

 t = 0
 if(modulo(iyear, 4) == 0) t = 1

 if(modulo(iyear, 400) /= 0 .and. modulo(iyear, 100) == 0) t = 0

 iday = idoy
 if(idoy > 59+t) iday = iday + 2 - t
 imonth = ((iday+91)*100)/3055
 iday = (iday+91) - (imonth*3055)/100
 imonth = imonth - 2

 if(imonth >= 1 .and. imonth <= 12) return
 write(unit=*,fmt="(a,i11,a)")  &
 "$$CALEND: DAY OF THE YEAR INPUT =",idoy," IS OUT OF RANGE."
 stop

end subroutine


subroutine saber_error(pres,err)

real(r8), intent(in)  :: pres
real(r8), intent(out) :: err

! NMP - simple function to calculate the error based on height
!  function loosely based on Remsberg et al (2008JGR)

err = 1.6718_r8 - 0.4281_r8*log(pres) + 0.0977_r8*log(pres)**2

end subroutine


end program

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
