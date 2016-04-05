
       module utils

       implicit none

       private
       public :: diffdat
       public :: heapsort
       public :: mapper
       public :: order
       public :: chk_map
       public :: get_julgmt
       public :: get_date
       public :: geth_newdate
       public :: is_leap_year
       public :: species_map
       public :: wrf2fire_type
       public :: diag_level

       integer :: diag_level = 100

       type species_map
          integer           :: fire_cnt                        ! count of fire species
          integer,pointer   :: fire_ndx(:,:)                   ! fire species indices
          real              :: wrf_conc                        ! wrf species assigned concentration
          real              :: m_wght                          ! molecular weight for aerosols (g/mole)
          character(len=32) :: wrf_name                        ! wrf species name
          real, pointer     :: fire_wght(:)                    ! multiplier for each fire species
          character(len=8), pointer :: fire_names(:)           ! fire species names
          logical           :: is_aerosol                      ! wrf species and fire species are aerosols
       end type species_map

       type(species_map), allocatable  :: wrf2fire_type(:)

       contains

      subroutine mapper( nspec, spc_map, specmax, n_files )
!------------------------------------------------------------------
!     decode wrf to fire mappings
!------------------------------------------------------------------

!------------------------------------------------------------------
!     dummy arguments
!------------------------------------------------------------------
      integer, intent(in)          :: specmax
      integer, intent(in)          :: n_files
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
      integer :: saerosol
      integer :: firecnt
      integer :: slen
      character(len=164) :: mstring
      logical :: one2one

      write(*,*) ' '
      nspec = 0
      do m = 1,specmax
         if( spc_map(m) == ' ' ) then
            nspec = m - 1
            exit
         else if( diag_level >= 200 ) then
            write(*,*) 'mapper: spc_map = ',trim(spc_map(m))
         end if
      end do

      if( nspec == 0 ) then
         write(*,*) 'mapper: No valid maps'
         stop 'Mapper error'
      end if
      allocate( wrf2fire_type(nspec),stat=status )
      if( status /= 0 ) then
         write(*,*) 'mapper: failed to allocate wrf2fire_type array; error = ',status
         stop 'Alloc error'
      end if

wrf_spec_loop : &
      do m = 1,nspec
         wrf2fire_type(m)%fire_cnt = 0
         sl = 0
         mstring = ' '
         wrf2fire_type(m)%wrf_name = ' '
!------------------------------------------------------------------
!	... remove blanks from string
!------------------------------------------------------------------
         do i = 1,len_trim( spc_map(m) )
            if( spc_map(m)(i:i) /= ' ' ) then
               sl = sl + 1
               mstring(sl:sl) = spc_map(m)(i:i)
            end if
         end do
         if( diag_level >= 200 ) then
           write(*,*) 'mapper: mstring = ',trim(mstring)
         endif
         sl = index( mstring, '->' )
         if( sl < 1 ) then
           sl = index( mstring, '=' )
           if( sl > 1 ) then
             wrf2fire_type(m)%wrf_name = mstring(:sl-1)
             read(mstring(sl+1:),*,iostat=status) wrf2fire_type(m)%wrf_conc
             if( status /= 0 ) then
                write(*,*) 'mapper; wrf concentration ',mstring(sl+1:),' is an invalid number'
                stop 'Read mapper'
             end if
             write(*,'('' mapper; wrf concentration '',a8,'' = '',g15.7)') mstring(sl+1:),wrf2fire_type(m)%wrf_conc
             cycle wrf_spec_loop
           else
             write(*,*) 'mapper: ',trim(spc_map(m)),' is invalid map'
             stop 'Mapper syntax error'
           end if 
         else if( sl == 1 ) then
           write(*,*) 'mapper: ',trim(spc_map(m)),' is invalid map'
           stop 'Mapper syntax error'
         end if 
         wrf2fire_type(m)%wrf_name = mstring(:sl-1)
         firecnt = 1
         sr = sl + 2
         sl = sr
         do
            su = index( mstring(sl:), '+' )
            if( su > 0 ) then
               firecnt = firecnt + 1
               sl = sl + su
            else
               exit
            end if
         end do
         write(*,'(''mapper: there are '',i2,'' fire emission species mapped to '',a)') firecnt,wrf2fire_type(m)%wrf_name
         allocate( wrf2fire_type(m)%fire_names(firecnt), &
                   wrf2fire_type(m)%fire_ndx(firecnt,n_files), &
                   wrf2fire_type(m)%fire_wght(firecnt),stat=status )
         if( status /= 0 ) then
            write(*,*) 'mapper: failed to allocate wrf2fire_type%fire_names; error = ',status
            stop 'Alloc mapper'
         end if
!------------------------------------------------------------------
!	... default values
!------------------------------------------------------------------
         wrf2fire_type(m)%fire_names(:firecnt) = ' '
         wrf2fire_type(m)%fire_ndx(:firecnt,:) = 0
         wrf2fire_type(m)%fire_wght(:firecnt)  = 1.
         wrf2fire_type(m)%fire_cnt = firecnt
         wrf2fire_type(m)%m_wght   = 12.
         wrf2fire_type(m)%is_aerosol = .false.
         firecnt = 1
         sl    = sr
         saerosol = index( mstring, ';' )
         if( saerosol > 0 ) then
            send  = saerosol - 1
         else
            send  = len_trim( mstring )
         end if
fire_spec_loop : &
         do
            su = index( mstring(sl:send), '+' )
            if( su > 0 ) then
               if( diag_level >= 400 ) then
                 write(*,*) 'mapper: su,string = ',su,mstring(sl:sl+su-2)
               endif
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
               read(mstring(sl:sl+snum-1),*,iostat=status) wrf2fire_type(m)%fire_wght(firecnt)
               if( status /= 0 ) then
                  write(*,*) 'mapper; coefficient ',mstring(sl:sl+snum-1),' is an invalid number'
                  stop 'Mapper syntax error'
               end if
               if( diag_level >= 300 ) then
                 write(*,*) 'mapper; coefficient ',mstring(sl:sl+snum-1),' = ',wrf2fire_type(m)%fire_wght(firecnt)
               end if
               sl   = sl + snum + 1
               slen = slen - (snum + 1)
               su   = su - (snum + 1)
            end if
!------------------------------------------------------------------
!	... check for species with no suffix
!------------------------------------------------------------------
            sbl = index( mstring(sl:sl+su-2), '[' )
            if( sbl > 0 ) then
               sbu = index( mstring(sl:sl+su-2), ']' )
               if( sbu > 0 ) then
                  slen = slen - 2
                  if( diag_level >= 400 ) then
                    write(*,*) 'mapper: string = ',mstring(sl:sl+su-2)
                  end if
                  wrf2fire_type(m)%fire_names(firecnt)(:slen) = mstring(sl+1:sl+su-3)
               end if
            else
               wrf2fire_type(m)%fire_names(firecnt)(:slen) = mstring(sl:sl+su-2)
            end if
            if( one2one ) then
               exit
            else
               firecnt = firecnt + 1
               if( sbl > 0 ) then
                  sl = sl + slen + 3
               else
                  sl = sl + slen + 1
               end if
               if( diag_level >= 400 ) then
                 write(*,*) 'mapper: sl,send = ',sl,send
                 write(*,*) 'mapper: string  = ',mstring(sl:send)
               end if
            end if
         end do fire_spec_loop
!------------------------------------------------------------------
!	... aerosol?
!------------------------------------------------------------------
         if( saerosol > 0 ) then
           saerosol = saerosol + 1
           if( mstring(saerosol:saerosol+6) == 'aerosol' ) then
             wrf2fire_type(m)%is_aerosol = .true.
             write(*,'(''mapper; '',a,'' is an aerosol'')') trim(wrf2fire_type(m)%wrf_name)
             sl = index( mstring(saerosol+7:), '(' )
             if( sl /= 0 ) then
               su = index( mstring(saerosol+7:), ')', back=.true. )
               sl = saerosol + 7 + sl
               su = saerosol + 5 + su
               if( su >= sl ) then
                 read(mstring(sl:su),*,iostat=status) wrf2fire_type(m)%m_wght
                 if( status /= 0 ) then
                   write(*,'(''mapper; molecular weight '',a,'' is an invalid number'')') mstring(sl:su)
                   stop 'Mapper syntax error'
                 endif
               endif
             endif
           end if
         end if
      end do wrf_spec_loop

      write(*,*) ' '
      do m = 1,nspec
        firecnt = wrf2fire_type(m)%fire_cnt
        if( firecnt > 0 ) then
          write(*,*)'mapper: ',m,trim(wrf2fire_type(m)%wrf_name),firecnt,wrf2fire_type(m)%fire_names(:firecnt), &
                               wrf2fire_type(m)%fire_wght(:firecnt),wrf2fire_type(m)%is_aerosol
        else
          write(*,*)'mapper: ',m,trim(wrf2fire_type(m)%wrf_name),firecnt,wrf2fire_type(m)%is_aerosol
        end if
      end do

      end subroutine mapper

      subroutine chk_map( n_wrf_spc, n_fire_spc, fire_species, file )
!------------------------------------------------------------------
!     check wrf to fire mappings
!------------------------------------------------------------------

!------------------------------------------------------------------
!     dummy arguments
!------------------------------------------------------------------
      integer, intent(in)           :: n_wrf_spc
      integer, intent(in)           :: n_fire_spc
      integer, intent(in)           :: file
      character(len=16), intent(in) :: fire_species(n_fire_spc)

!------------------------------------------------------------------
!     local variables
!------------------------------------------------------------------
      integer :: i, m, n
      character(len=16) :: spcnam
      logical :: found

      write(*,*) ' '
      do n = 1,n_wrf_spc
        write(*,'(''chk_map: checking output variable '',a)') trim(wrf2fire_type(n)%wrf_name)
        do i = 1,wrf2fire_type(n)%fire_cnt
          spcnam = trim(wrf2fire_type(n)%fire_names(i))
          found = .false.
          do m = 1,n_fire_spc
            if( trim(spcnam) == trim(fire_species(m)) ) then
              found = .true.
              wrf2fire_type(n)%fire_ndx(i,file) = m
              exit
            endif
          end do
          if( .not. found ) then
            write(*,*) 'chk_map: could not find ',spcnam,' in fire species'
            stop 'Mapper species error'
          end if
        end do
      end do

      end subroutine chk_map

   subroutine heapsort( n, ra, ndx )
 
   integer, intent(in) :: n
   integer, intent(inout) :: ndx(n)
   real, intent(inout)    :: ra(n)

   integer :: i, ir, j, l
   integer :: nh
   real    :: rra

   if( n >= 2 ) then
     l  = n/2 + 1
     ir = n
outer_loop : &
     do
       if( l > 1 ) then
         l = l - 1
         rra = ra(l)
         nh = ndx(l)
       else
         rra = ra(ir)
         nh = ndx(ir)
         ra(ir) = ra(1)
         ndx(ir) = ndx(1)
         ir = ir - 1
         if( ir == 1 ) then
           ra(1) = rra
           ndx(1) = nh
           exit outer_loop
         endif
       endif
       i = l
       j = l + l
inner_loop : &
       do
         if( j <= ir ) then
           if( j < ir .and. ra(j) < ra(j+1) ) then
             j = j + 1
           endif
           if( rra < ra(j) ) then
             ra(i) = ra(j)
             ndx(i) = ndx(j)
             i = j
             j = 2*j
           else
             j = ir + 1
           endif
         else
           exit inner_loop
         endif
         ra(i) = rra
         ndx(i) = nh
       end do inner_loop
     end do outer_loop
   endif

   end subroutine heapsort

   subroutine order( ns, src, nd, dst, ndx )

   integer, intent(in) :: ns, nd
   integer, intent(in) :: ndx(nd)
   real,    intent(in) :: src(ns)
   real,    intent(out) :: dst(nd)

   dst(:) = src(ndx(:))

   end subroutine order

   SUBROUTINE get_julgmt( date_str, julyr, julday, gmt )

! dummy Arguments
     INTEGER, INTENT(OUT  ) :: julyr
     INTEGER, INTENT(OUT  ) :: julday
     REAL   , INTENT(OUT  ) :: gmt
     CHARACTER(LEN=*) , INTENT(IN) :: date_str

! Local variables
     INTEGER :: ny , nm , nd , nh , ni , ns
     INTEGER :: my1, my2, my3
     INTEGER :: mmd(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
     CHARACTER(LEN=19) :: ldate_str

     ldate_str = date_str // '_00:00:00'
     CALL split_date_char ( ldate_str , ny , nm , nd , nh , ni , ns )

     GMT = nh + real(ni)/60. + real(ns)/3600.
     MY1 = MOD(ny,4)
     MY2 = MOD(ny,100)
     MY3 = MOD(ny,400)
     JULDAY = nd + sum( MMD(1:nm-1) )
     JULYR  = ny
     IF( nm > 2 .and. (MY1 == 0 .AND. MY2 /= 0 .OR. MY3 == 0) ) then
       julday = julday + 1
     endIF

   END SUBROUTINE get_julgmt

   SUBROUTINE get_date( date_str, julyr, julday )

! dummy Arguments
     INTEGER, INTENT(IN  ) :: julyr
     INTEGER, INTENT(IN  ) :: julday
     CHARACTER(LEN=*) , INTENT(OUT) :: date_str

! Local variables
     INTEGER :: m
     INTEGER :: my1, my2, my3
     INTEGER :: nd
     INTEGER :: mmd(12)
     CHARACTER(LEN=3) :: numa

     mmd(:) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
     if( is_leap_year( julyr ) ) then
       mmd(2) = 29
     endIF

     nd = julday
     do m = 1,11
       nd = nd - mmd(m)
       if( nd <= 0 ) then
         nd = nd + mmd(m)
         exit
       endif
     end do
     
     write(date_str(1:4),'(i4)') julyr
     write(numa,'(i3)') m+100
     date_str(5:7) = '-' // numa(2:3)
     write(numa,'(i3)') nd+100
     date_str(8:10) = '-' // numa(2:3)

   END SUBROUTINE get_date

   logical function is_leap_year( year )

   integer, intent(in) :: year

   integer :: my1, my2, my3

   my1 = mod(year,4)
   my2 = mod(year,100)
   my3 = mod(year,400)

   is_leap_year = (MY1 == 0 .AND. MY2 /= 0) .OR. MY3 == 0

   end function is_leap_year

   SUBROUTINE geth_newdate (ndate, odate, idt)
   
      !  From old date ('YYYY-MM-DD HH:MM:SS.ffff') and 
      !            [or ('YYYY-DDDDD HH:MM:SS.ffff')]
      !  delta-time, compute the new date.
   
      !  on entry     -  odate  -  the old hdate.
      !                  idt    -  the change in time
   
      !  on exit      -  ndate  -  the new hdate.
      
      INTEGER , INTENT(IN)           :: idt
      CHARACTER (LEN=*) , INTENT(OUT) :: ndate
      CHARACTER (LEN=*) , INTENT(IN)  :: odate
      
       
      !  Local Variables
       
      !  yrold    -  indicates the year associated with "odate"
      !  moold    -  indicates the month associated with "odate"
      !  dyold    -  indicates the day associated with "odate"
      !  hrold    -  indicates the hour associated with "odate"
      !  miold    -  indicates the minute associated with "odate"
      !  scold    -  indicates the second associated with "odate"
       
      !  yrnew    -  indicates the year associated with "ndate"
      !  monew    -  indicates the month associated with "ndate"
      !  dynew    -  indicates the day associated with "ndate"
      !  hrnew    -  indicates the hour associated with "ndate"
      !  minew    -  indicates the minute associated with "ndate"
      !  scnew    -  indicates the second associated with "ndate"
       
      !  mday     -  a list assigning the number of days in each month
      
      !  i        -  loop counter
      !  nday     -  the integer number of days represented by "idt"
      !  nhour    -  the integer number of hours in "idt" after taking out
      !              all the whole days
      !  nmin     -  the integer number of minutes in "idt" after taking out
      !              all the whole days and whole hours.
      !  nsec     -  the integer number of minutes in "idt" after taking out
      !              all the whole days, whole hours, and whole minutes.
       
      INTEGER :: nlen, olen
      INTEGER :: yrnew, monew, dynew, hrnew, minew, scnew, frnew
      INTEGER :: yrold, moold, dyold, hrold, miold, scold, frold
      INTEGER :: mday(12), nday, nhour, nmin, nsec, nfrac, i, ifrc
      LOGICAL :: opass
      CHARACTER (LEN=10) :: hfrc
      CHARACTER (LEN=1) :: sp
      ! INTEGER, EXTERNAL :: nfeb  ! in the same module now
      
      !  Assign the number of days in a months
      
      mday( 1) = 31
      mday( 2) = 28
      mday( 3) = 31
      mday( 4) = 30
      mday( 5) = 31
      mday( 6) = 30
      mday( 7) = 31
      mday( 8) = 31
      mday( 9) = 30
      mday(10) = 31
      mday(11) = 30
      mday(12) = 31
      
      !  Break down old hdate into parts
      
      hrold = 0
      miold = 0
      scold = 0
      frold = 0
      olen = LEN(odate)
      IF (olen.GE.11) THEN
         sp = odate(11:11)
      else
         sp = ' '
      END IF
      
      !  Use internal READ statements to convert the CHARACTER string
      !  date into INTEGER components.
   
      READ(odate(1:4),  '(I4)') yrold
      READ(odate(6:7),  '(I2)') moold
      READ(odate(9:10), '(I2)') dyold
      IF (olen.GE.13) THEN
         READ(odate(12:13),'(I2)') hrold
         IF (olen.GE.16) THEN
            READ(odate(15:16),'(I2)') miold
            IF (olen.GE.19) THEN
               READ(odate(18:19),'(I2)') scold
               IF (olen.GT.20) THEN
                  READ(odate(21:olen),'(I2)') frold
               END IF
            END IF
         END IF
      END IF
      
      !  Set the number of days in February for that year.
      
      mday(2) = nfeb(yrold)
      
      !  Check that ODATE makes sense.
      
      opass = .TRUE.
      
      !  Check that the month of ODATE makes sense.
      
      IF ((moold.GT.12).or.(moold.LT.1)) THEN
         WRITE(*,*) 'GETH_NEWDATE:  Month of ODATE = ', moold
         opass = .FALSE.
      END IF
      
      !  Check that the day of ODATE makes sense.
      
      IF ((dyold.GT.mday(moold)).or.(dyold.LT.1)) THEN
         WRITE(*,*) 'GETH_NEWDATE:  Day of ODATE = ', dyold
         opass = .FALSE.
      END IF
      !  Check that the hour of ODATE makes sense.
      
      IF ((hrold.GT.23).or.(hrold.LT.0)) THEN
         WRITE(*,*) 'GETH_NEWDATE:  Hour of ODATE = ', hrold
         opass = .FALSE.
      END IF
      
      !  Check that the minute of ODATE makes sense.
      
      IF ((miold.GT.59).or.(miold.LT.0)) THEN
         WRITE(*,*) 'GETH_NEWDATE:  Minute of ODATE = ', miold
         opass = .FALSE.
      END IF
      
      !  Check that the second of ODATE makes sense.
      
      IF ((scold.GT.59).or.(scold.LT.0)) THEN
         WRITE(*,*) 'GETH_NEWDATE:  Second of ODATE = ', scold
         opass = .FALSE.
      END IF
      
      !  Check that the fractional part  of ODATE makes sense.
      
      
      IF (.not.opass) THEN
         WRITE(*,*) 'module_date_time: GETH_NEWDATE: Bad ODATE: ', odate(1:olen), olen
         stop
      END IF
      
      !  Date Checks are completed.  Continue.
      
      
      !  Compute the number of days, hours, minutes, and seconds in idt
      
      IF (olen.GT.20) THEN !idt should be in fractions of seconds
         ifrc = olen-20
         ifrc = 10**ifrc
         nday   = ABS(idt)/(86400*ifrc)
         nhour  = MOD(ABS(idt),86400*ifrc)/(3600*ifrc)
         nmin   = MOD(ABS(idt),3600*ifrc)/(60*ifrc)
         nsec   = MOD(ABS(idt),60*ifrc)/(ifrc)
         nfrac = MOD(ABS(idt), ifrc)
      ELSE IF (olen.eq.19) THEN  !idt should be in seconds
         ifrc = 1
         nday   = ABS(idt)/86400 ! Integer number of days in delta-time
         nhour  = MOD(ABS(idt),86400)/3600
         nmin   = MOD(ABS(idt),3600)/60
         nsec   = MOD(ABS(idt),60)
         nfrac  = 0
      ELSE IF (olen.eq.16) THEN !idt should be in minutes
         ifrc = 1
         nday   = ABS(idt)/1440 ! Integer number of days in delta-time
         nhour  = MOD(ABS(idt),1440)/60
         nmin   = MOD(ABS(idt),60)
         nsec   = 0
         nfrac  = 0
      ELSE IF (olen.eq.13) THEN !idt should be in hours
         ifrc = 1
         nday   = ABS(idt)/24 ! Integer number of days in delta-time
         nhour  = MOD(ABS(idt),24)
         nmin   = 0
         nsec   = 0
         nfrac  = 0
      ELSE IF (olen.eq.10) THEN !idt should be in days
         ifrc = 1
         nday   = ABS(idt)/24 ! Integer number of days in delta-time
         nhour  = 0
         nmin   = 0
         nsec   = 0
         nfrac  = 0
      ELSE
         WRITE(*,*) 'module_date_time: GETH_NEWDATE: Strange length for ODATE: ',olen
         stop
      END IF
      
      IF (idt.GE.0) THEN
      
         frnew = frold + nfrac
         IF (frnew.GE.ifrc) THEN
            frnew = frnew - ifrc
            nsec = nsec + 1
         END IF
      
         scnew = scold + nsec
         IF (scnew .GE. 60) THEN
            scnew = scnew - 60
            nmin  = nmin + 1
         END IF
      
         minew = miold + nmin
         IF (minew .GE. 60) THEN
            minew = minew - 60
            nhour  = nhour + 1
         END IF
      
         hrnew = hrold + nhour
         IF (hrnew .GE. 24) THEN
            hrnew = hrnew - 24
            nday  = nday + 1
         END IF
      
         dynew = dyold
         monew = moold
         yrnew = yrold
         DO i = 1, nday
            dynew = dynew + 1
            IF (dynew.GT.mday(monew)) THEN
               dynew = dynew - mday(monew)
               monew = monew + 1
               IF (monew .GT. 12) THEN
                  monew = 1
                  yrnew = yrnew + 1
                  ! If the year changes, recompute the number of days in February
                  mday(2) = nfeb(yrnew)
               END IF
            END IF
         END DO
      
      ELSE IF (idt.LT.0) THEN
      
         frnew = frold - nfrac
         IF (frnew .LT. 0) THEN
            frnew = frnew + ifrc
            nsec = nsec - 1
         END IF
      
         scnew = scold - nsec
         IF (scnew .LT. 00) THEN
            scnew = scnew + 60
            nmin  = nmin + 1
         END IF
      
         minew = miold - nmin
         IF (minew .LT. 00) THEN
            minew = minew + 60
            nhour  = nhour + 1
         END IF
      
         hrnew = hrold - nhour
         IF (hrnew .LT. 00) THEN
            hrnew = hrnew + 24
            nday  = nday + 1
         END IF
      
         dynew = dyold
         monew = moold
         yrnew = yrold
         DO i = 1, nday
            dynew = dynew - 1
            IF (dynew.eq.0) THEN
               monew = monew - 1
               IF (monew.eq.0) THEN
                  monew = 12
                  yrnew = yrnew - 1
                  ! If the year changes, recompute the number of days in February
                  mday(2) = nfeb(yrnew)
               END IF
               dynew = mday(monew)
            END IF
         END DO
      END IF
      
      !  Now construct the new mdate
      
      nlen = LEN(ndate)
      
      IF (nlen.GT.20) THEN
         WRITE(ndate(1:19),19) yrnew, monew, dynew, hrnew, minew, scnew
         WRITE(hfrc,'(I10)') frnew+1000000000
         ndate = ndate(1:19)//'.'//hfrc(31-nlen:10)
      
      ELSE IF (nlen.eq.19.or.nlen.eq.20) THEN
         WRITE(ndate(1:19),19) yrnew, monew, dynew, hrnew, minew, scnew
      19   format(I4,'-',I2.2,'-',I2.2,'_',I2.2,':',I2.2,':',I2.2)
         IF (nlen.eq.20) ndate = ndate(1:19)//'.'
      
      ELSE IF (nlen.eq.16) THEN
         WRITE(ndate,16) yrnew, monew, dynew, hrnew, minew
      16   format(I4,'-',I2.2,'-',I2.2,'_',I2.2,':',I2.2)
      
      ELSE IF (nlen.eq.13) THEN
         WRITE(ndate,13) yrnew, monew, dynew, hrnew
      13   format(I4,'-',I2.2,'-',I2.2,'_',I2.2)
      
      ELSE IF (nlen.eq.10) THEN
         WRITE(ndate,10) yrnew, monew, dynew
      10   format(I4,'-',I2.2,'-',I2.2)
      
      END IF
      
      IF (olen.GE.11) ndate(11:11) = sp

   END SUBROUTINE geth_newdate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   FUNCTION nfeb ( year ) RESULT (num_days)
   
      ! Compute the number of days in February for the given year
   
      INTEGER :: year
      INTEGER :: num_days
   
      num_days = 28 ! By default, February has 28 days ...
      IF (MOD(year,4).eq.0) THEN  
         num_days = 29  ! But every four years, it has 29 days ...
         IF (MOD(year,100).eq.0) THEN
            num_days = 28  ! Except every 100 years, when it has 28 days ...
            IF (MOD(year,400).eq.0) THEN
               num_days = 29  ! Except every 400 years, when it has 29 days.
            END IF
         END IF
      END IF
   
   END FUNCTION nfeb

   SUBROUTINE split_date_char ( date , century_year , month , day , hour , minute , second )
     
!  dummy args
   
      INTEGER , INTENT(OUT) :: century_year , month , day , hour , minute , second
      CHARACTER(LEN=*) , INTENT(IN) :: date 

      INTEGER :: istat
      
      READ(date,'(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)',iostat=istat) century_year,month,day,hour,minute,second
!     if( istat /= 0 ) then
!       write(*,*) 'split_date_char: failed to read date; error = ',istat
!       stop 'date_error'
!     endif
   
   END SUBROUTINE split_date_char

      real function DIFFDAT( dates )
!-----------------------------------------------------------------------
! 	... Compute the difference: dates(2) - dates(1)  in days.
!           Return value:
!           dates(2) - dates(1)  in days.
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
! 	... input arguments
!-----------------------------------------------------------------------
      character(len=*), intent(in) :: dates(2)   ! dates in yyyy-mm-dd format
!-----------------------------------------------------------------------
! 	... local variables
!-----------------------------------------------------------------------
      integer :: yr, mnth, day
      integer :: n
      integer :: istat
      integer :: dash_pos(2)
      integer :: dats(2)   ! dates in yyyymmdd format

      do n = 1,2
        dash_pos(1) = index( dates(n), '-' )
        dash_pos(2) = index( dates(n), '-', back=.true. )
        if( any( dash_pos(:) <= 0 ) .or. dash_pos(1) >= dash_pos(2) ) then
          write(*,*) 'diffdat: improper date format in one or both arguments'
          write(*,*) 'diffdat: date1 = ',trim(dates(1))
          write(*,*) 'diffdat: date2 = ',trim(dates(2))
          stop 'DateError'
        endif
        read(dates(n)(1:dash_pos(1)-1),*,iostat=istat) yr
        if( istat /= 0 ) then
          write(*,*) 'diffdat: date, ',trim(dates(n)),' has improper year'
          stop 'DateError'
        endif
        read(dates(n)(dash_pos(1)+1:dash_pos(2)-1),*,iostat=istat) mnth
        if( istat /= 0 ) then
          write(*,*) 'diffdat: date, ',trim(dates(n)),' has improper month'
          stop 'DateError'
        endif
        read(dates(n)(dash_pos(2)+1:),*,iostat=istat) day
        if( istat /= 0 ) then
          write(*,*) 'diffdat: date, ',trim(dates(n)),' has improper day'
          stop 'DateError'
        endif
        dats(n) = 10000*yr + 100*mnth + day
      end do

      DIFFDAT = DIFFDATGRG( dats(1), 0, dats(2), 0 )

      end function DIFFDAT

      real function DIFFDATGRG( dat1, sec1, dat2, sec2 )
!-----------------------------------------------------------------------
! 	... Compute the difference: (dat2,sec2) - (dat1,sec1)  in days.
!           N.B. Assume Gregorian calendar.
!           Return value:
!           (dat2,sec2) - (dat1,sec1)  in days.
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
! 	... Dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) :: &
        dat1, &   ! date in yyyymmdd format
        sec1, &   ! seconds relative to dat1
        dat2, &   ! date in yyyymmdd format
        sec2      ! seconds relative to dat2


!-----------------------------------------------------------------------
! 	... Local variables
!-----------------------------------------------------------------------
      integer :: lastd, lasts, firstd, firsts, ndays
      real    :: sign, days

!-----------------------------------------------------------------------
!     	... Dates equal?
!-----------------------------------------------------------------------
      if( dat1 == dat2 .and. sec1 == sec2 ) then
         DIFFDATGRG = 0.
         return
      end if

!-----------------------------------------------------------------------
!     	... Which date is later?
!-----------------------------------------------------------------------
      if( dat2 > dat1 ) then
         sign   = 1.
         lastd  = dat2
         lasts  = sec2
         firstd = dat1
         firsts = sec1
      else if( dat2 < dat1 ) then
         sign   = -1.
         lastd  = dat1
         lasts  = sec1
         firstd = dat2
         firsts = sec2
      else
         if( sec2 > sec1 ) then
            sign   = 1.
            lastd  = dat2
            lasts  = sec2
            firstd = dat1
            firsts = sec1
         else
            sign   = -1.
            lastd  = dat1
            lasts  = sec1
            firstd = dat2
            firsts = sec2
         end if
      end if

!-----------------------------------------------------------------------
!     	... Compute number of days between lastd and firstd
!-----------------------------------------------------------------------
      ndays = GREG2JDAY( lastd ) - GREG2JDAY( firstd )

!-----------------------------------------------------------------------
!     	... Adjust for remaining seconds
!-----------------------------------------------------------------------
      days = REAL( ndays ) + REAL( lasts - firsts )/86400.

!-----------------------------------------------------------------------
!     	... Adjust sign
!-----------------------------------------------------------------------
      DIFFDATGRG = sign * days

      end function DIFFDATGRG
!-----------------------------------------------------------------------
! Currently the routines greg2jday and jday2greg are working in sense of
! j = greg2jday( jday2greg( j ) ) down to jday 1684595 which corresponds
! approx. to the date -1000228.
!-----------------------------------------------------------------------

      integer function GREG2JDAY( date )
!-----------------------------------------------------------------------
! 	... Return Julian day number given Gregorian date.
!
! Algorithm from Hatcher,D.A., Simple Formulae for Julian Day Numbers
! and Calendar Dates, Q.Jl.R.astr.Soc. (1984) v25, pp 53-55.
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!	... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) :: date

!-----------------------------------------------------------------------
!	... Local variables
!-----------------------------------------------------------------------
      integer :: yy, mm, dd
      integer :: ap, mp
      integer :: y, d, n, g

!-----------------------------------------------------------------------
!     	... Extract year, month, and day from date
!-----------------------------------------------------------------------
      yy = date / 10000
      mm = MOD( ABS(date),10000 ) / 100
      dd = MOD( ABS(date),100 )

!-----------------------------------------------------------------------
!     	... Modify year and month numbers
!-----------------------------------------------------------------------
      ap = yy - (12 - mm)/10
      mp = MOD( mm-3,12 )
      if( mp < 0 ) then
         mp = mp + 12
      end if

!-----------------------------------------------------------------------
!     	... Julian day
!-----------------------------------------------------------------------
      y = INT( 365.25*( ap + 4712 ) )
      d = INT( 30.6*mp + .5 )
      n = y + d + dd  + 59
      g = INT( .75*INT( ap/100 + 49 ) ) - 38
      GREG2JDAY = n - g

      end function GREG2JDAY

      integer function JDAY2GREG( day )
!-----------------------------------------------------------------------
! 	... Return Gregorian date given Julian day number.
!
! Algorithm from Hatcher,D.A., Simple Formulae for Julian Day Numbers
! and Calendar Dates, Q.Jl.R.astr.Soc. (1984) v25, pp 53-55.
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!	... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) :: day

!-----------------------------------------------------------------------
!	... Local variables
!-----------------------------------------------------------------------
      integer :: g, n
      integer :: yy, dp, mm, dd

      g = INT( .75*INT( (day - 4479.5)/36524.25 ) + .5 ) - 37
      n = day + g

      yy = INT( n/365.25 ) - 4712
      dp = INT( MOD( n-59.25, 365.25 ) )
      mm = MOD( INT( (dp+.5)/30.6 )+2, 12 ) + 1
      dd = INT( MOD( dp+.5, 30.6 ) ) + 1

      JDAY2GREG = SIGN( ABS(yy)*10000 + mm*100 +dd,yy )

      end function JDAY2GREG

   end module utils
