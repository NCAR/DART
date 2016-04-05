
       module utils

       implicit none

       private
       public :: wrf2mz_time
       public :: mz2wrf_time
       public :: mapper
       public :: species_map
       public :: wrf2mz_map

       type species_map
          integer           :: moz_cnt                         ! count of mozart species
          real              :: wrf_wght                        ! overall multiplier
          real              :: wrf_conc                        ! wrf species assigned concentration
          character(len=32) :: wrf_name                        ! wrf species name
          real, pointer     :: moz_wght(:)                     ! multiplier for each mozart species
          character(len=8), pointer :: moz_names(:)           ! mozart species names
          logical, pointer  :: moz_ext(:)                      ! mozart species uses preset suffix
       end type species_map

       type(species_map), allocatable  :: wrf2mz_map(:)

       contains

       subroutine wrf2mz_time( wrf_time, ncdate, ncsec )
!------------------------------------------------------------------
!     convert wrf time string to mozart date,sec format
!------------------------------------------------------------------

!------------------------------------------------------------------
!     dummy arguments
!------------------------------------------------------------------
       character(len=*), intent(in) :: wrf_time
       integer, intent(out)         :: ncdate
       integer, intent(out)         :: ncsec

!------------------------------------------------------------------
!     local variables
!------------------------------------------------------------------
      integer :: year, month, day
      integer :: hour, minute, second
      integer :: date_time_sep
      integer :: nstring
      integer :: istat

      nstring = len( wrf_time )
      date_time_sep = index( wrf_time, '_' )
      if( date_time_sep < 1 .or. date_time_sep >= nstring ) then
         write(*,*) 'wrf2mz_time: wrf date ',wrf_time(:nstring),' is invalid'
      end if
!------------------------------------------------------------------
!     read year, month and day
!------------------------------------------------------------------
       read(wrf_time(:date_time_sep-1),'(i4,1x,i2,1x,i2)',iostat=istat) year, month, day
       if( istat /= 0 ) then
         write(*,*) 'wrf2mz_time: failed to read year,month,day; iostat = ',istat
       end if
       read(wrf_time(date_time_sep+1:nstring),'(i2,1x,i2,1x,i2)',iostat=istat) hour, minute, second
       if( istat /= 0 ) then
         write(*,*) 'wrf2mz_time: failed to read hour,minute,second; iostat = ',istat
       end if
       ncdate = 100*(year*100 + month) + day
       ncsec  = 60*(60*hour + minute) + second

      end subroutine wrf2mz_time

      subroutine mz2wrf_time( wrf_time, ncdate, ncsec )
!------------------------------------------------------------------
!     convert wrf time string to mozart date,sec format
!------------------------------------------------------------------

!------------------------------------------------------------------
!     dummy arguments
!------------------------------------------------------------------
       character(len=*), intent(out) :: wrf_time
       integer, intent(in)           :: ncdate
       integer, intent(in)           :: ncsec

!------------------------------------------------------------------
!     local variables
!------------------------------------------------------------------
      integer :: year, month, day
      integer :: hour, minute, second
      integer :: nstring
      integer :: istat
      character(len=3) :: numa

      wrf_time = ' '
      year  = ncdate/10000
      month = mod( ncdate,10000 )/100
      day   = mod( ncdate,100)
      write(wrf_time,'(i4,''-'')',iostat=istat) year
      write(numa,'(i3)') 100+month
      wrf_time(6:8) = numa(2:3) // '-'
      write(numa,'(i3)') 100+day
      wrf_time(9:11) = numa(2:3) // '_'
      hour   = ncsec/3600
      minute = (ncsec - hour*3600)/60
      second = ncsec - 60*(hour*60 + minute)
      write(numa,'(i3)') 100+hour
      wrf_time(12:14) = numa(2:3) // ':'
      write(numa,'(i3)') 100+minute
      wrf_time(15:17) = numa(2:3) // ':'
      write(numa,'(i3)') 100+second
      wrf_time(18:19) = numa(2:3)

      end subroutine mz2wrf_time

      subroutine mapper( nspec, spc_map, specmax )
!------------------------------------------------------------------
!     decode wrf to mozart mappings
!------------------------------------------------------------------

!------------------------------------------------------------------
!     dummy arguments
!------------------------------------------------------------------
      integer, intent(in)          :: specmax
      integer, intent(out)         :: nspec
      character(len=*), intent(in) :: spc_map(specmax)

!------------------------------------------------------------------
!     local variables
!------------------------------------------------------------------
      integer :: i, m
      integer :: status
      integer :: sl, su, sr
      integer :: sbl, sbu
      integer :: send
      integer :: snum
      integer :: swght
      integer :: mzcnt
      integer :: slen
      character(len=164) :: mstring
      logical :: one2one

      nspec = 0
      do m = 1,specmax
         if( spc_map(m) == ' ' ) then
            nspec = m - 1
            exit
         else
            write(*,*) 'mapper: spc_map = ',trim(spc_map(m))
         end if
      end do

      if( nspec == 0 ) then
         write(*,*) 'mapper: No valid maps'
         stop 'mapper'
      end if
      allocate( wrf2mz_map(nspec),stat=status )
      if( status /= 0 ) then
         write(*,*) 'mapper: failed to allocate wrf2mz_map array; error = ',status
         stop 'mapper'
      end if

wrf_spec_loop : &
      do m = 1,nspec
         wrf2mz_map(m)%moz_cnt = 0
         sl = 0
         mstring = ' '
         wrf2mz_map(m)%wrf_name = ' '
!------------------------------------------------------------------
!	... remove blanks from string
!------------------------------------------------------------------
         do i = 1,len_trim( spc_map(m) )
            if( spc_map(m)(i:i) /= ' ' ) then
               sl = sl + 1
               mstring(sl:sl) = spc_map(m)(i:i)
            end if
         end do
         write(*,*) 'mapper: mstring = ',trim(mstring)
         sl = index( mstring, '->' )
         if( sl < 1 ) then
           sl = index( mstring, '=' )
           if( sl > 1 ) then
             wrf2mz_map(m)%wrf_name = mstring(:sl-1)
             read(mstring(sl+1:),*,iostat=status) wrf2mz_map(m)%wrf_conc
             if( status /= 0 ) then
                write(*,*) 'mapper; wrf concentration ',mstring(sl+1:),' is an invalid number'
                stop 'mapper'
             end if
             write(*,'('' mapper; wrf concentration '',a8,'' = '',g15.7)') mstring(sl+1:),wrf2mz_map(m)%wrf_conc
             cycle wrf_spec_loop
           else
             write(*,*) 'mapper: ',trim(spc_map(m)),' is invalid map'
             stop 'mapper'
           end if 
         else if( sl == 1 ) then
           write(*,*) 'mapper: ',trim(spc_map(m)),' is invalid map'
           stop 'mapper'
         end if 
         wrf2mz_map(m)%wrf_name = mstring(:sl-1)
         wrf2mz_map(m)%wrf_wght = 1.e6
         mzcnt = 1
         sr = sl + 2
         sl = sr
         do
            su = index( mstring(sl:), '+' )
            if( su > 0 ) then
               mzcnt = mzcnt + 1
               sl = sl + su
            else
               exit
            end if
         end do
         write(*,*) 'mapper: there are ',mzcnt,' mozart species for ',wrf2mz_map(m)%wrf_name
         allocate( wrf2mz_map(m)%moz_names(mzcnt), &
                   wrf2mz_map(m)%moz_ext(mzcnt), &
                   wrf2mz_map(m)%moz_wght(mzcnt),stat=status )
         if( status /= 0 ) then
            write(*,*) 'mapper: failed to allocate wrf2mz_map%moz_names; error = ',status
            stop 'mapper'
         end if
!------------------------------------------------------------------
!	... default values
!------------------------------------------------------------------
         wrf2mz_map(m)%moz_names(:mzcnt) = ' '
         wrf2mz_map(m)%moz_wght(:mzcnt)  = 1.
         wrf2mz_map(m)%moz_ext(:mzcnt)   = .true.
         wrf2mz_map(m)%moz_cnt = mzcnt
         mzcnt = 1
         sl    = sr
         swght = index( mstring, ';' )
         if( swght > 0 ) then
            send  = swght - 1
         else
            send  = len_trim( mstring )
         end if
         write(*,*) 'mapper: send = ',send
moz_spec_loop : &
         do
            su = index( mstring(sl:send), '+' )
            if( su > 0 ) then
               write(*,*) 'mapper: su,string = ',su,mstring(sl:sl+su-2)
               slen = su - 1
               one2one = .false.
            else
               slen = send - sl + 1
               su   = send + 2
               one2one = .true.
            end if
!------------------------------------------------------------------
!	... check for coefficient
!------------------------------------------------------------------
            snum = index( mstring(sl:sl+su-2), '*' ) - 1
            if( snum > 0 ) then
               read(mstring(sl:sl+snum-1),*,iostat=status) wrf2mz_map(m)%moz_wght(mzcnt)
               if( status /= 0 ) then
                  write(*,*) 'mapper; coefficient ',mstring(sl:sl+snum-1),' is an invalid number'
                  stop 'mapper'
               end if
               write(*,*) 'mapper; coefficient ',mstring(sl:sl+snum-1),' = ',wrf2mz_map(m)%moz_wght(mzcnt)
               sl   = sl + snum + 1
               slen = slen - (snum + 1)
               su   = su - (snum + 1)
            end if
!------------------------------------------------------------------
!	... check for species with no mozart suffix
!------------------------------------------------------------------
            sbl = index( mstring(sl:sl+su-2), '[' )
            if( sbl > 0 ) then
               sbu = index( mstring(sl:sl+su-2), ']' )
               if( sbu > 0 ) then
                  slen = slen - 2
                  write(*,*) 'mapper: string = ',mstring(sl:sl+su-2)
                  wrf2mz_map(m)%moz_names(mzcnt)(:slen) = mstring(sl+1:sl+su-3)
                  wrf2mz_map(m)%moz_ext(mzcnt)          = .false.
               end if
            else
               wrf2mz_map(m)%moz_names(mzcnt)(:slen) = mstring(sl:sl+su-2)
            end if
            if( one2one ) then
               exit
            else
               mzcnt = mzcnt + 1
               if( sbl > 0 ) then
                  sl = sl + slen + 3
               else
                  sl = sl + slen + 1
               end if
               write(*,*) 'mapper: sl,send = ',sl,send
               write(*,*) 'mapper: string  = ',mstring(sl:send)
            end if
         end do moz_spec_loop
!------------------------------------------------------------------
!	... check for wrf conversion factor
!------------------------------------------------------------------
         if( swght > 0 ) then
            read(mstring(swght+1:),*,iostat=status) wrf2mz_map(m)%wrf_wght
            if( status /= 0 ) then
               write(*,*) 'mapper; wrf_wght ',mstring(swght+1:),' is an invalid number'
               stop 'mapper'
            end if
            write(*,*) 'mapper; wrf_wght ',mstring(swght+1:),' = ',wrf2mz_map(m)%wrf_wght
         end if
      end do wrf_spec_loop

      do m = 1,nspec
        mzcnt = wrf2mz_map(m)%moz_cnt
        if( mzcnt > 0 ) then
          write(*,*)'mapper: ',m,trim(wrf2mz_map(m)%wrf_name),mzcnt,wrf2mz_map(m)%moz_names(:mzcnt), &
                                wrf2mz_map(m)%moz_wght(:mzcnt),wrf2mz_map(m)%moz_ext(:mzcnt), &
                                wrf2mz_map(m)%wrf_wght
        else
          write(*,*)'mapper: ',m,trim(wrf2mz_map(m)%wrf_name),mzcnt,wrf2mz_map(m)%wrf_conc
        end if
      end do

      end subroutine mapper

      end module utils
