
   program fire_emis
!-----------------------------------------------------------------
!     transform raw fire emissions to wrf,mozart,cam-chem domain
!-----------------------------------------------------------------

   use utils
   use wrf_utils
   use glb_utils, only : glb_file, write_glb_fire_file, glb_file_final
   use fire_file, only : write_fire_file, max_fire_size, dealloc_fire_emis_glb_atts
   use srf_types, only : fire_srf_init, fire_srf_types, fire_srf_final
   use srf_types, only : ntypes
   use attr_types

   implicit none

   integer, parameter :: map_size = 500
   integer, parameter :: maxsize  = 132
   integer, parameter :: max_files = 5

   integer :: nlines
   integer :: n_infiles
   integer :: astat, istat
   integer :: slen, comma, chr
   integer :: k, m, m1, m2, n, n1
   integer :: file, filep1
   integer :: n_fire_spc
   integer :: n_wrf_spc
   integer :: genveg
   integer :: domain
   integer :: last_time_ndx
   integer :: ndx_cnt
   integer :: day_cnt
   integer :: day_ndx
   integer :: match_day
   integer :: start_julday, end_julday
   integer :: start_file, end_file
   integer :: start_mnth, end_mnth
   integer :: start_yr, end_yr
   integer :: n_mnths
   integer :: n_total_days
   integer :: ngatts
   integer, allocatable :: lat_ndx(:)
   integer, allocatable :: lon_ndx(:)
   integer, allocatable :: gen_ndx(:)
   integer, allocatable :: beg_day_ndx(:)
   integer, allocatable :: day(:,:)
   integer :: julyr(2)
   integer :: julday(2)
   integer :: file_yr(max_files)
   integer :: file_lincnt(max_files)
   real    :: min_x, max_x
   real    :: min_y, max_y
   real    :: x, y
   real    :: timeod
   real    :: wrk_emiss
   real    :: dxsqi
   real    :: avrg_fire_size = 2.5e5    ! m^2, .25 km^2
   real, allocatable    :: lon(:)
   real, allocatable    :: lat(:)
   real, allocatable    :: fire_size(:)
   real, allocatable    :: fire_emissions(:,:)
   character(len=10)    :: wrk_date
   character(len=10)    :: output_start_date
   character(len=32)    :: dummy
   character(len=200)   :: filespec
   character(len=1024)  :: buffer
   character(len=16), allocatable   :: fire_species(:)
   logical              :: has_fire_size(max_files)
   logical              :: is_wrf, is_glb
   logical              :: fexist
   logical              :: found
   type(proj_info), allocatable :: proj(:)


   type(dom_glb_att), pointer :: dom_glb_attrs(:)

!-----------------------------------------------------------------
!     control variables
!-----------------------------------------------------------------
   integer                :: domains = 1
   character(len=maxsize) :: fire_directory = ' '
   character(len=maxsize) :: wrf_directory = ' '
   character(len=maxsize) :: fire_filename(max_files) = ' '
   character(len=maxsize) :: srf_type_filename(ntypes) = ' '
   character(len=164)     :: wrf2fire_map(map_size) = ' '
   character(len=164)     :: glb2fire_map(map_size) = ' '
   character(len=16)      :: resol = 'WRF'
   character(len=10)      :: start_date
   character(len=10)      :: end_date
   character(len=7)       :: output_timing = 'daily  '

   namelist /control/ domains, fire_directory, wrf_directory, fire_filename, &
                      srf_type_filename, wrf2fire_map, start_date, end_date, &
                      avrg_fire_size, max_fire_size, diag_level, resol,      &
                      output_timing, glb2fire_map

!-------------------------------------------------------------------------
!  read control variables
!-------------------------------------------------------------------------
   read(*,nml=control,iostat=istat)
   if( istat /= 0 ) then
     write(*,*) 'fire_emis: failed to read namelist; error = ',istat
     stop
   end if

   is_wrf = trim( resol ) == 'WRF'
   is_glb = .not. is_wrf
!-----------------------------------------------------------------
!     check namelist inputs
!-----------------------------------------------------------------
   if( domains < 1 ) then
     write(*,*) 'fire_emis: domains must be >= 1'
     stop 'Namelist error'
   else if( is_glb .and. domains > 1 ) then
     write(*,*) 'fire_emis: only one global domain allowed'
     stop 'Namelist error'
   end if

   if( diag_level >= 200 ) then
     write(*,*) 'fire_emis: start,end dates = ',start_date,' ',end_date
   endif

   if( diffdat( (/ start_date, end_date /) ) < 0. ) then
     write(*,*) 'fire_emis: start date > end date'
     stop 'Date error'
   endif

is_montly : &
   if( is_glb .and. output_timing == 'monthly' ) then
     wrk_date = start_date
!-----------------------------------------------------------------
!  check start date is beginning of month
!-----------------------------------------------------------------
     n = 0
     do
       read(wrk_date(9:10),*) m1
       if( m1 == 1 ) then
         start_date = wrk_date
         exit
       endif
       n = n + 24
       call geth_newdate( wrk_date, start_date, n )
       if( diffdat( (/ wrk_date, end_date /) ) <= 0. ) then
         write(*,*) 'fire_emis: monthly global output - but can not find'
         write(*,*) '           whole month between start and end dates'
         stop 'Namelist error'
       endif
     end do
!-----------------------------------------------------------------
!  check end date is end of month
!-----------------------------------------------------------------
     call geth_newdate( wrk_date, end_date, 24 )
     read(wrk_date(9:10),*) n1
     if( n1 /= 1 ) then
       read(end_date(6:7),*) m2
       n = 0
       do
         n = n - 24
         call geth_newdate( wrk_date, end_date, n )
         if( diffdat( (/ start_date, wrk_date /) ) <= 0. ) then
           write(*,*) 'fire_emis: monthly global output - but can not find'
           write(*,*) '           whole month between start and end dates'
           stop 'Namelist error'
         endif
         read(wrk_date(6:7),*) m1
         if( m2 /= m1 ) then
           end_date = wrk_date
           exit
         endif
       end do
     endif
   endif is_montly

   call get_julgmt( start_date, julyr(1), julday(1), timeod )
   if( diag_level >= 200 ) then
     write(*,*) 'fire_emis: start yr,day = ',julyr(1),julday(1)
   endif
   call get_julgmt( end_date, julyr(2), julday(2), timeod )
   if( diag_level >= 200 ) then
     write(*,*) 'fire_emis: end yr,day = ',julyr(2),julday(2)
   endif

!-------------------------------------------------------------------------
!  check start, end dates for proper order
!-------------------------------------------------------------------------
   n_total_days = int( diffdat( (/ start_date, end_date /) ) )
   if( n_total_days < 0 ) then
     write(*,*) 'fire_emis: start date > end date'
     stop 'Date error'
   endif
   n_total_days = n_total_days + 1

   if( is_glb .and. output_timing == 'monthly' ) then
     read(start_date(1:4),*) start_yr
     read(start_date(1:4),*) end_yr
     read(start_date(6:7),*) start_mnth
     read(end_date(6:7),*)   end_mnth
     n_mnths = (end_yr - start_yr)*12 + end_mnth - start_mnth + 1
   endif

!-------------------------------------------------------------------------
!  get total file count
!-------------------------------------------------------------------------
   do n_infiles = 1,max_files
     if( trim(fire_filename(n_infiles)) == ' ' ) then
       exit
     endif      
   end do
   n_infiles = n_infiles - 1
   if( n_infiles == 0 ) then
     write(*,*) 'fire_emis: no fire emission files specified in namelist'
     stop 'Namelist error'
   endif

!-------------------------------------------------------------------------
!  check for basic fire emission file existence
!-------------------------------------------------------------------------
   do file = 1,n_infiles
     filespec = trim(fire_directory) // trim(fire_filename(file))
     inquire( file=trim(filespec),exist=fexist )
     if( .not. fexist ) then
       write(*,*) ' '
       write(*,*) 'fire_emis: file ',trim(filespec),' does not exist'
       stop 'File error'
     endif
   end do

!-------------------------------------------------------------------------
!  check for end,start dates requiring 2 or more files
!-------------------------------------------------------------------------
   if( julyr(2) > julyr(1) .and. n_infiles == 1 ) then
     write(*,*) 'fire_emis: end year > start year but only one fire input file'
     stop 'Namelist error'
   endif

   if( is_glb ) then
     wrf2fire_map(:) = glb2fire_map(:)
   endif
!-----------------------------------------------------------------
!     form wrf to fire emission map
!-----------------------------------------------------------------
   call mapper( n_wrf_spc, wrf2fire_map, map_size, n_infiles )
   write(*,*) ' '
   write(*,'(''fire_emis: there are '',i3,'' output species'')') n_wrf_spc

   if( is_glb ) then
     allocate( proj(domains),stat=astat )
     if( astat /= 0 ) then
       write(*,*) 'fire_emis: failed to allocate proj; error = ',astat
       stop 'Alloc error'
     endif
     call glb_file( resol, proj(1), n_wrf_spc, n_total_days, start_date, &
                    n_infiles, fire_directory, fire_filename, output_timing, n_mnths )
   endif

file_loop : &
   do file = 1,n_infiles
!-------------------------------------------------------------------------
!  open basic fire emissions file
!-------------------------------------------------------------------------
     filespec = trim(fire_directory) // trim(fire_filename(file))
     open( unit=10,file=trim(filespec),iostat=istat )
     if( istat /= 0 ) then
       write(*,*) 'fire_emis: failed to open base fire emissions file ',trim(filespec),'; error = ',istat
       stop 'Open error'
     endif
!-------------------------------------------------------------------------
!  get length of fire emission file(s)
!-------------------------------------------------------------------------
     write(*,*) ' '
     write(*,*) 'fire_emis: getting fire record count for ',trim(filespec)
     write(*,*) 'fire_emis: this usually takes about 30 seconds'
     nlines = get_file_size( 10 ) - 1
     file_lincnt(file) = nlines
!-------------------------------------------------------------------------
!  read and process the header record
!-------------------------------------------------------------------------
     read(unit=10,fmt='(a)',iostat=istat) buffer
     if( istat /= 0 ) then
       write(*,*) 'fire_emis: failed to read ',trim(filespec),'; error = ',istat
       stop 'Read error'
     endif
!-------------------------------------------------------------------------
!  does fire emiss file have fire size?
!-------------------------------------------------------------------------
     has_fire_size(file) = index( trim(buffer),'AREA' ) /= 0
     if( is_wrf ) then
       write(*,*) ' '
       if( has_fire_size(file) ) then
         write(*,*) 'fire_emis: ',trim(filespec),' has fire size data'
       else
         write(*,*) 'fire_emis: ',trim(filespec),' does NOT have fire size data'
       endif
     endif

     comma = iachar( ',' )
     n = len_trim( buffer )
     allocate( gen_ndx(n),stat=astat )
     if( astat /= 0 ) then
       write(*,*) 'fire_emis: failed to allocate gen_ndx; error = ',astat
       stop 'Alloc error'
     endif
     do m = 1,n
       gen_ndx(m) = iachar( buffer(m:m) )
     enddo
     if( has_fire_size(file) ) then
       n_fire_spc = count(gen_ndx(:) == comma) - 5
     else
       n_fire_spc = count(gen_ndx(:) == comma) - 4
     endif
     deallocate( gen_ndx )
     allocate( fire_species(n_fire_spc),stat=astat )
     if( astat /= 0 ) then
       write(*,*) 'fire_emis: failed to allocate species; error = ',astat
       stop 'Alloc error'
     endif

     rewind(10)
     if( has_fire_size(file) ) then
       read(unit=10,fmt=*,iostat=istat) dummy,dummy,dummy,dummy,dummy,dummy,fire_species(1:n_fire_spc)
     else
       read(unit=10,fmt=*,iostat=istat) dummy,dummy,dummy,dummy,dummy,fire_species(1:n_fire_spc)
     endif
     if( istat /= 0 ) then
       write(*,*) 'fire_emis: failed to read fire.data; error = ',istat
       stop 'Read error'
     endif

     slen = len_trim(fire_species(n_fire_spc))
     chr  = iachar(fire_species(n_fire_spc)(slen:slen))
     if( .not. ((iachar('a') <= chr .and. chr <= iachar('z')) .or. &
                (iachar('A') <= chr .and. chr <= iachar('Z')) .or. &
                (iachar('0') <= chr .and. chr <= iachar('9'))) ) then
       fire_species(n_fire_spc) = fire_species(n_fire_spc)(1:len_trim(fire_species(n_fire_spc))-1)
     endif
     write(*,*) ' '
     write(*,*) 'fire_emis: fire species in file ',trim(fire_filename(file))
     write(*,*) ' '
     do m = 1,n_fire_spc,5
       write(*,'(5(1x,a16))') fire_species(m:min(m+4,n_fire_spc))
     enddo
!-----------------------------------------------------------------
!     check mapping
!-----------------------------------------------------------------
     call chk_map( n_wrf_spc, n_fire_spc, fire_species, file )
     deallocate( fire_species )
     write(*,*) ' '
     close( 10 )
   end do file_loop

!-------------------------------------------------------------------------
!  allocate day array
!-------------------------------------------------------------------------
   if( allocated( day ) ) then
     deallocate( day )
   endif
   allocate( day(maxval(file_lincnt(1:n_infiles)),n_infiles),stat=astat )
   if( astat /= 0 ) then
     write(*,*) 'fire_emis: failed to allocate day array; error = ',astat
     stop 'Alloc error'
   endif

file_loop_a : &
   do file = 1,n_infiles
     filespec = trim(fire_directory) // trim(fire_filename(file))
     open( unit=10,file=trim(filespec),iostat=istat )
     if( istat /= 0 ) then
       write(*,*) 'fire_emis: failed to open base fire emissions file ',trim(filespec),'; error = ',istat
       stop 'Open error'
     endif
     read(unit=10,fmt='(a)',iostat=istat) buffer
     if( istat /= 0 ) then
       write(*,*) 'fire_emis: failed to read ',trim(filespec),'; error = ',istat
       stop 'Read error'
     endif
     write(*,*) ' '
     write(*,*) 'fire_emis: getting julian days for ',trim(filespec)
     write(*,*) 'fire_emis: this usually takes about 30 seconds'
!-------------------------------------------------------------------------
!  read the days
!-------------------------------------------------------------------------
     do n = 1,file_lincnt(file)
       read(unit=10,fmt=*,iostat=istat) day(n,file)
       if( istat /= 0 ) then
         write(*,'(''fire_emis: failed to read '',a,'' at record = '',i10,''; error = '',i5)') trim(fire_filename(file)),n,istat
         close( 10 )
         stop 'Read error'
       endif
     end do
     close( 10 )
!-------------------------------------------------------------------------
!  check days
!-------------------------------------------------------------------------
     if( minval( day(1:file_lincnt(file),file) ) < 1 ) then
       write(*,*) 'fire_emis: julian days in ',trim(fire_filename(file)),' improperly ordered'
       stop 'Day numbering error'
     endif
     if( any( day(2:file_lincnt(file),file) < day(1:file_lincnt(file)-1,file) ) ) then
       write(*,*) 'fire_emis: julian days in ',trim(fire_filename(file)),' improperly ordered'
       stop 'Day numbering error'
     endif
   end do file_loop_a

!-------------------------------------------------------------------------
!  check output start, end dates
!-------------------------------------------------------------------------
   file_yr(1) = julyr(1)
file_loop_b : &
   do file = 1,n_infiles-1
     filep1 = file + 1
     if( .not. is_leap_year( file_yr(file) ) ) then
       if( day(file_lincnt(file),file) == 365 ) then
         if( day(1,filep1) /= 1 ) then
           write(*,'(''fire_emis: input files '',i1,'', '',i1,'' do NOT have consecutive julian day numbering'')') file,filep1
           stop 'Day numbering error'
         endif
         file_yr(filep1) = file_yr(file) + 1
       else
         if( day(1,filep1) /= day(file_lincnt(file),file)+1 ) then
           write(*,'(''fire_emis: input files '',i1,'', '',i1,'' do NOT have consecutive julian day numbering'')') file,filep1
           stop 'Day numbering error'
         endif
         file_yr(filep1) = file_yr(file)
       endif
     else
       if( day(file_lincnt(file),file) == 366 ) then
         if( day(1,filep1) /= 1 ) then
           write(*,'(''fire_emis: input files '',i1,'', '',i1,'' do NOT have consecutive julian day numbering'')') file,filep1
           stop 'Day numbering error'
         endif
         file_yr(filep1) = file_yr(file) + 1
       else
         if( day(1,filep1) /= day(file_lincnt(file),file)+1 ) then
           write(*,'(''fire_emis: input files '',i1,'', '',i1,'' do NOT have consecutive julian day numbering'')') file,filep1
           stop 'Day numbering error'
         endif
         file_yr(filep1) = file_yr(file)
       endif
     endif
   end do file_loop_b

   found = .false.
   do file = 1,n_infiles
     if( file_yr(file) == julyr(1) ) then
       if( julday(1) >= day(1,file) .and. julday(1) <= day(file_lincnt(file),file) ) then
         found = .true.
         exit
       endif
     endif
   end do
   if( .not. found ) then
     write(*,*) 'fire_emis: start output date not in any fire input file '
     stop 'Date error'
   endif
   start_file = file

   found = .false.
   do file = n_infiles,1,-1
     if( file_yr(file) == julyr(2) ) then
       if( julday(2) >= day(1,file) .and. julday(2) <= day(file_lincnt(file),file) ) then
         found = .true.
         exit
       endif
     endif
   end do
   if( .not. found ) then
     write(*,*) 'fire_emis: end output date not in any fire input file '
     stop 'Date error'
   endif
   end_file = file

   call geth_newdate( output_start_date, start_date, -24 )

   if( is_wrf ) then
     allocate( proj(domains),stat=astat )
     if( astat /= 0 ) then
       write(*,*) 'fire_emis: failed to allocate proj; error = ',astat
       stop 'Alloc error'
     endif
     call fire_srf_init( domains )
     allocate( dom_glb_attrs(domains),stat=astat )
     if( astat /= 0 ) then
       write(*,*) 'fire_emis: failed to allocate dom_glb_attrs; error = ',astat
       stop 'Alloc error'
     endif
   endif

file_loop_c : &
   do file = start_file,end_file
     if( file == start_file ) then
       start_julday = julday(1)
     else
       start_julday = day(1,file)
     endif
     if( file == end_file ) then
       end_julday = julday(2)
     else
       end_julday = day(file_lincnt(file),file)
     endif
     day_cnt = end_julday - start_julday + 2
     allocate( beg_day_ndx(day_cnt),stat=astat )
     if( astat /= 0 ) then
       write(*,*) 'fire_emis: failed to allocate day_ndx; error = ',astat
       stop 'Alloc error'
     endif
!-------------------------------------------------------------------------
!  find begin,end julian day indices in fire file
!-------------------------------------------------------------------------
     found = .false.
     do n = 1,file_lincnt(file)
       if( day(n,file) == start_julday ) then
         found = .true.
         beg_day_ndx(1) = n
         exit
       endif
     end do
     if( .not. found ) then
       write(*,*) 'fire_emis: requested start date ',start_date,' is not in dataset'
       call get_date( start_date, file_yr(file), day(1,file) )
       call get_date( end_date, file_yr(file), day(file_lincnt(file),file) )
       write(*,*) 'fire_emis: file start,end dates = ',start_date,' ',end_date
       stop 'Date error'
     endif
     found = .false.
     do last_time_ndx = file_lincnt(file),beg_day_ndx(1),-1
       if( day(last_time_ndx,file) == end_julday ) then
         found = .true.
         exit
       endif
     end do
     if( .not. found ) then
       write(*,*) 'fire_emis: requested end date ',end_date,' is not in dataset'
       call get_date( start_date, file_yr(file), day(1,file) )
       call get_date( end_date, file_yr(file), day(file_lincnt(file),file) )
       write(*,*) 'fire_emis: file start,end dates = ',start_date,' ',end_date
       stop 'Date error'
     endif

     if( diag_level >= 200 ) then
       write(*,*) 'fire_emis: beg_day_ndx,last_time_ndx = ',beg_day_ndx(1),last_time_ndx
     endif

     match_day = start_julday - 1
     do n = 2,day_cnt-1
       match_day = match_day + 1
       beg_day_ndx(n) = beg_day_ndx(n-1) + count( day(1:file_lincnt(file),file) == match_day )
     end do
     match_day = match_day + 1
     beg_day_ndx(day_cnt) = beg_day_ndx(day_cnt-1) + count( day(1:file_lincnt(file),file) == match_day )
     filespec = trim(fire_directory) // trim(fire_filename(file))
     open( unit=10,file=trim(filespec),iostat=istat )
     if( istat /= 0 ) then
       write(*,*) 'fire_emis: failed to open base fire emissions file ',trim(filespec),'; error = ',istat
       stop 'Open error'
     endif

domain_loop : &
     do domain = 1,domains
!-------------------------------------------------------------------------
!  read the wrf file parameters
!-------------------------------------------------------------------------
       if( is_wrf .and. file == start_file ) then
         call wrf_file( domain, wrf_directory, dxsqi, proj(domain), &
                        dom_glb_attrs(domain)%ngatts, dom_glb_attrs(domain)%attrs )
         dxsqi = 1./(dxsqi*dxsqi)
       endif
       min_x = .5
       max_x = real(proj(domain)%ide) + .5
       if( is_wrf ) then
         min_y = .5
         max_y = real(proj(domain)%jde) + .5
       else
         min_y = 1.
         max_y = real(proj(domain)%jde)
       endif
!-----------------------------------------------------------------------
!  form the srf type variables
!-----------------------------------------------------------------------
       if( is_wrf .and. file == start_file ) then
         call fire_srf_types( fire_directory, domain, proj(domain) )
       endif
!-------------------------------------------------------------------------
!  position the fire emission file before first output day
!-------------------------------------------------------------------------
       if( domain > 1 ) then
         rewind( unit = 10 )
       endif
       read(unit=10,fmt=*,iostat=istat) buffer
       if( istat /= 0 ) then
         write(*,*) 'fire_emis: failed to read ',trim(filespec),'; error = ',istat
         stop 'Read error'
       endif
       do n = 1,beg_day_ndx(1)-1
         read(unit=10,fmt=*,iostat=istat) timeod
         if( istat /= 0 ) then
           write(*,*) 'fire_emis: failed to read ',trim(filespec),'; error = ',istat
           stop 'Read error'
         endif
       end do
!-------------------------------------------------------------------------
!  read fire emissions for output day
!-------------------------------------------------------------------------
       wrk_date = output_start_date
day_loop : &
       do day_ndx = 1,day_cnt-1
         call geth_newdate( wrk_date, wrk_date, 24 )
         m1 = beg_day_ndx(day_ndx)
         m2 = beg_day_ndx(day_ndx+1)-1
         allocate( lon(m1:m2),lat(m1:m2),lon_ndx(m1:m2),lat_ndx(m1:m2), &
                   fire_size(m1:m2), fire_emissions(n_fire_spc,m1:m2),stat=astat )
         if( astat /= 0 ) then
           write(*,*) 'fire_emis: failed to allocate lon ... fire_emissions arrays; error = ',astat
           stop 'Alloc error'
         endif
         lon_ndx(:) = 0
         do m = m1,m2
           if( has_fire_size(file) ) then
             read(unit=10,fmt=*,iostat=istat) timeod,timeod,genveg,lat(m),lon(m),fire_size(m),fire_emissions(:,m)
           else
             read(unit=10,fmt=*,iostat=istat) timeod,timeod,genveg,lat(m),lon(m),fire_emissions(:,m)
             fire_size(m) = avrg_fire_size
           endif
           if( istat /= 0 ) then
             write(*,*) 'fire_emis: failed to read fire.data; error = ',istat
             stop 'Read error'
           endif
           if( is_glb .and. lon(m) < 0. ) then
             lon(m) = mod( lon(m)+360.,360. )
           endif
           call llij( lat(m), lon(m), proj(domain), x, y )
           if( is_glb .and. x > max_x ) then
             x = x - real(proj(domain)%ide)
           endif
           if( is_wrf ) then
             if( min_x <= x .and. x < max_x .and. &
                 min_y <= y .and. y < max_y ) then
               lon_ndx(m) = nint( x )
               lat_ndx(m) = nint( y )
             endif
           else
             lon_ndx(m) = mod( nint( x ) - 1,proj(domain)%ide ) + 1
             lat_ndx(m) = mod( nint( y ) - 1,proj(domain)%jde ) + 1
           endif
         end do
         if( is_wrf ) then
           call write_fire_file( domain, wrk_date, n_wrf_spc, m1, m2, &
                                 lon, lon_ndx, lat_ndx, n_fire_spc, fire_size, &
                                 fire_emissions, dxsqi, proj(domain), n_infiles, fire_directory, &
                                 fire_filename, file, dom_glb_attrs(domain)%ngatts, dom_glb_attrs(domain)%attrs )
         else
           call write_glb_fire_file( file, n_wrf_spc, m1, m2, lon_ndx, &
                                     lat_ndx, n_fire_spc, fire_emissions, proj(domain), output_timing, &
                                     wrk_date )
         endif
         deallocate( lon,lat,lon_ndx,lat_ndx,fire_size,fire_emissions,stat=astat )
         if( astat /= 0 ) then
           write(*,*) 'fire_emis: failed to deallocate lon ... fire_emissions arrays; error = ',astat
           stop 'Dealloc error'
         endif
       end do day_loop
     end do domain_loop
     close( unit=10 )
     deallocate( beg_day_ndx,stat=astat )
     if( astat /= 0 ) then
       write(*,*) 'fire_emis: failed to deallocate beg_day_ndx array; error = ',astat
       stop 'Dealloc error'
     endif
     output_start_date = wrk_date
   end do file_loop_c

   if( is_wrf ) then
     do domain = 1,domains
       if( dom_glb_attrs(domain)%ngatts > 0 ) then
         call dealloc_fire_emis_glb_atts( dom_glb_attrs(domain)%ngatts, dom_glb_attrs(domain)%attrs )
         deallocate( dom_glb_attrs(domain)%attrs )
       endif
     end do
     deallocate( dom_glb_attrs )
     call fire_srf_final( domains )
   else
     call glb_file_final
   endif

   write(*,*) ' '
   write(*,*) '================================='
   write(*,*) 'fire_emis: Completed successfully'
   write(*,*) '================================='

   contains

   integer function get_file_size( unitno )
!-----------------------------------------------------------------
!  get file size
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!  dummy arguments
!-----------------------------------------------------------------
   integer, intent(in) :: unitno
!-----------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------
   integer :: nl
   character :: c

   nl = 0
   do
     read(unitno,*,iostat=istat) c
     if( istat == 0 ) then
       nl = nl + 1
     else
       exit
     endif
   end do
   
   rewind( unitno )
   get_file_size = nl

   end function get_file_size

   end program fire_emis
