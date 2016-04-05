
      module utils

      use anthro_types, only : anthro_map_type, dates
      use anthro_types, only : linemax, linsize, maxsrc, namsize

      implicit none

      private
      public :: mapper
      public :: upcase
      public :: wrf2mz_time, mz2wrf_time
      public :: cleanup_grid
      public :: start_data, stop_data
      public :: anthro_map
      public :: molecw
      public :: diag_level
      public :: n_sub_cats
      public :: sub_cats
      public :: src_names, nsrc_names, src_active

      integer :: nsrc_names
      integer :: n_sub_cats = 0
      integer :: diag_level = 100

      character(len=namsize), allocatable :: sub_cats(:)

      character(len=4) :: sub_cats_defs(9) =        &
                             (/ 'agr ', 'awb ', 'dom ', &
                                'ene ', 'ind ', 'slv ', &
                                'tra ', 'wst ', 'ship' /)
      logical           :: src_active(maxsrc) = .false.
      character(len=namsize) :: src_names(maxsrc)
      real              :: molecw(maxsrc) = 0.               ! src species molecular weight (amu)

      type(anthro_map_type), allocatable  :: anthro_map(:)

      type(dates) :: start_data
      type(dates) :: stop_data

      contains

      subroutine mapper( nemis, emis_map, sub_categories )
!------------------------------------------------------------------
!     decode src to anthro emission mappings
!------------------------------------------------------------------

!------------------------------------------------------------------
!     dummy arguments
!------------------------------------------------------------------
      integer, intent(out)         :: nemis
      character(len=linsize), intent(in) :: emis_map(linemax)
      character(len=namsize), intent(in) :: sub_categories(maxsrc)

!------------------------------------------------------------------
!     local variables
!------------------------------------------------------------------
      integer, parameter :: maxtok = 50
      real, parameter    :: cwght = 1.

      integer :: i, j, k, m, n, nl, nu
      integer :: astat, ios
      integer :: sl, su, sr
      integer :: sbl, sbu
      integer :: send
      integer :: snum
      integer :: nsrc
      integer :: ncat
      integer :: saerosol
      integer :: nlines
      integer :: slen
      integer :: tokcnt
      integer :: cat_tokcnt
      integer :: ndx
      integer :: cat_ndx
      integer :: rpar_cnt, lpar_cnt, plus_cnt
      integer :: toklen(maxtok)
      integer :: rpar_ndx(maxsrc)
      integer :: lpar_ndx(maxsrc)
      integer :: plus_ndx(maxsrc)
      integer, allocatable :: spc_lineno(:)
      integer, allocatable :: cat_toklen(:)
      character(len=linsize) :: wstring, sub_string, cp_string
      character(len=34) :: tst_strng
      character(len=namsize) :: tokens(maxtok)
      character(len=namsize), allocatable :: cat_tokens(:)
      logical :: found

!------------------------------------------------------------------
!     species count
!------------------------------------------------------------------
      write(*,*) ' '
      nlines = 0
      nemis  = 0
      do m = 1,linemax
        if( emis_map(m) == ' ' ) then
          nlines = m - 1
          exit
        else if( diag_level >= 200 ) then
          write(*,*) 'mapper: emis_map = ',trim(emis_map(m))
        endif
        if( index( emis_map(m), '->' ) /= 0 ) then
          nemis = nemis + 1
        endif
      end do

      if( nemis == 0 ) then
         write(*,*) 'mapper: No valid maps'
         stop 'Mapper error'
      end if
      allocate( anthro_map(nemis),spc_lineno(nemis+1),stat=astat )
      if( astat /= 0 ) then
        write(*,*) 'mapper: failed to allocate anthro_map array; error = ',astat
        stop 'Alloc error'
      end if

!------------------------------------------------------------------
!     check src names
!------------------------------------------------------------------
      do nsrc_names = 1,maxsrc
        if( src_names(nsrc_names) == ' ' ) then
          exit
        endif
      end do
      nsrc_names = min( nsrc_names-1,maxsrc )
      if( nsrc_names < 1 ) then
        write(*,*) 'mapper: source species names not specified via src_names variable'
        stop 'Namelist error'
      endif
!------------------------------------------------------------------
!     check for src molecular weight
!------------------------------------------------------------------
      do n = 1,nsrc_names
        sbl = index( src_names(n), '(' )
        if( sbl > 1 ) then
          sbu = index( src_names(n), ')' )
          if( sbu > sbl ) then
             molecw(n) = get_number( src_names(n)(sbl+1:sbu-1) )
             src_names(n)(sbl:) = ' '
          else
            write(*,*) 'mapper: molecular weight improperly specified'
            write(*,*) '        input = ',trim(src_names(n) )
            stop 'Input error'
          endif
        endif
      end do
!------------------------------------------------------------------
!     begin, end line for each emission species map
!------------------------------------------------------------------
      n = 0
      do m = 1,nlines
        ndx = index( emis_map(m), '->' )
        if( ndx > 0 ) then
          if( ndx > 1 ) then
            n = n + 1
            spc_lineno(n) = m
          else
            write(*,*) 'mapper: emission species name invalid'
            stop 'Syntax error'
          endif
        endif
      end do
      spc_lineno(nemis+1) = nlines + 1
!------------------------------------------------------------------
!     sub category count
!------------------------------------------------------------------
      n_sub_cats  = 0
      do m = 1,maxsrc
        if( sub_categories(m) == ' ' ) then
          n_sub_cats = m - 1
          exit
        endif
      end do

      found = .true.
      if( n_sub_cats == 0 ) then
        n_sub_cats = 9
        found = .false.
      end if
      allocate( sub_cats(n_sub_cats),cat_tokens(n_sub_cats),cat_toklen(n_sub_cats),stat=astat )
      if( astat /= 0 ) then
        write(*,*) 'mapper: failed to allocate sub_cats .. cat_toklen array; error = ',astat
        stop 'Alloc error'
      end if
      if( found ) then
        sub_cats(:n_sub_cats) = sub_categories(:n_sub_cats)
      else
        sub_cats(:n_sub_cats) = sub_cats_defs(:n_sub_cats)
      endif
      do m = 1,nemis
        allocate( anthro_map(m)%cat_wght(n_sub_cats,nsrc_names),stat=astat )
        if( astat /= 0 ) then
          write(*,*) 'mapper: failed to allocate cat_wght array; error = ',astat
          stop 'Alloc error'
        end if
      end do

!------------------------------------------------------------------
!     determine number of source species in emission species mapping
!------------------------------------------------------------------
emis_spec_loop : &
      do m = 1,nemis
        nl = spc_lineno(m)
        nu = spc_lineno(m+1) - 1
        nsrc = 0
        anthro_map(m)%src_cnt = 0
        anthro_map(m)%is_gas  = .true.
line_loop : &
        do n = nl,nu
          cp_string = emis_map(n)
          call remove_blanks( cp_string )
          if( n == nl ) then
            ndx = index( cp_string, '->' )
            wstring = cp_string(ndx+2:)
            su = index( cp_string(:ndx-1), '(a)' )
            if( su == 0 ) then
              anthro_map(m)%emis_name = cp_string(:ndx-1)
            else
              anthro_map(m)%is_gas = .false.
              anthro_map(m)%emis_name = cp_string(:ndx-4)
            endif
          else
            wstring = cp_string
            if( wstring(1:1) /= '+' ) then
              write(*,*) 'mapper: input continuation line'
              write(*,*) trim(wstring)
              write(*,*) '        does not begin with "+"'
              stop 'Syntax err'
            endif
            wstring = cp_string(2:)
          endif
          call get_indices( '(', wstring, lpar_cnt, lpar_ndx )
          call get_indices( ')', wstring, rpar_cnt, rpar_ndx )
          call get_indices( ')+', wstring, plus_cnt, plus_ndx )
!------------------------------------------------------------------
!     basic syntax checks
!------------------------------------------------------------------
          if( lpar_cnt /= rpar_cnt ) then
            write(*,*) 'mapper: unbalanced () in emis_map'
            stop 'Syntax error'
          endif
          do i = 1,lpar_cnt
            if( lpar_ndx(i) >= rpar_ndx(i) ) then
              write(*,*) 'mapper: unbalanced () in emis_map'
              stop 'Syntax error'
            endif
          end do
          do i = 1,lpar_cnt-1
            k = rpar_ndx(i) + 1
            if( wstring(k:k) /= '+' ) then
              write(*,*) 'mapper: input line'
              write(*,*) trim(wstring)
              stop 'Syntax error'
            endif
          end do
!------------------------------------------------------------------
!     if present remove category specification
!------------------------------------------------------------------
          if( lpar_cnt > 0 ) then
            sub_string = ' '
            sl = 1
            su = 1
            do i = 1,lpar_cnt
              sub_string(su:) = wstring(sl:lpar_ndx(i)-1)
              su = len_trim(sub_string)+1 
              sl = rpar_ndx(i)+1
            end do
            sub_string(su:) = wstring(sl:)
          else
            sub_string = wstring(:)
          endif

!------------------------------------------------------------------
!     set the sources
!------------------------------------------------------------------
          call gettokens( sub_string, len_trim(sub_string), '+', 32, tokens, &
                          toklen, maxtok, tokcnt )

          ncat = 0
src_tok_loop : &
          do k = 1,tokcnt
            nsrc = nsrc + 1
            ndx  = index( tokens(k), '*' )
            if( ndx /= 0 ) then
              anthro_map(m)%src_wght(nsrc) = get_number( tokens(k)(:ndx-1) )
            else
              anthro_map(m)%src_wght(nsrc) = 1.
            endif
            ndx  = ndx + 1
            anthro_map(m)%src_var(nsrc) = trim( tokens(k)(ndx:) )
            found = .false.
            do j = 1,nsrc_names
              if( trim(anthro_map(m)%src_var(nsrc)) == trim(src_names(j)) ) then
                found = .true.
                exit
              endif
            end do
            if( .not. found ) then
              write(*,*) 'mapper: ',trim(anthro_map(m)%src_var(nsrc)),' not in source species list'
              stop 'Missing file err'
            endif
!------------------------------------------------------------------
!     set the source categories
!------------------------------------------------------------------
            tst_strng = trim( tokens(k) ) // '('
            sbl = index( wstring, trim(tst_strng) )
has_cat :   if( sbl > 0 ) then
              ncat = ncat + 1
              anthro_map(m)%cat_wght(:,nsrc) = 0.
!------------------------------------------------------------------
!     set the source categories
!------------------------------------------------------------------
              sbl = lpar_ndx(ncat)+1
              sbu = rpar_ndx(ncat)-1
              call gettokens( wstring(sbl:sbu), sbu-sbl+1, '+', namsize, cat_tokens, &
                              cat_toklen, n_sub_cats, cat_tokcnt )
              do j = 1,cat_tokcnt
                ndx  = index( cat_tokens(j), '*' )
                cat_ndx = get_cat_ndx( cat_tokens(j)(ndx+1:) )
                if( ndx /= 0 ) then
                  anthro_map(m)%cat_wght(cat_ndx,nsrc) = get_number( cat_tokens(j)(:ndx-1) )
                else
                  anthro_map(m)%cat_wght(cat_ndx,nsrc) = cwght
                endif
              end do
            else has_cat
              anthro_map(m)%cat_wght(:,nsrc) = cwght
            endif has_cat
            src_active(get_src_ndx(trim(anthro_map(m)%src_var(nsrc)))) = .true.
          end do src_tok_loop
          anthro_map(m)%src_cnt = anthro_map(m)%src_cnt + tokcnt
        end do line_loop
      end do emis_spec_loop

      deallocate( spc_lineno, cat_tokens, cat_toklen )

      end subroutine mapper

      subroutine gettokens( string, ls, delim, maxlen, tokens, &
                            toklen, maxtok, tokcnt )

      implicit none

!-----------------------------------------------------------------------     
!     Input arguments:
!        string - character string to crack into tokens
!        ls     - length of string
!        delim  - token delimiter character
!        maxlen - maximum length of any single token
!        maxtok - maximum number of tokens
!     Output arguments:
!        tokcnt - number of actual tokens
!                 < 0 => hit maxtok before end of string
!                 = 0 => error in input string
!        toklen - array containing length of each token
!        tokens - character array of tokens
!-----------------------------------------------------------------------     

      integer, intent(in)  ::  ls, maxlen, maxtok
      integer, intent(out) ::  tokcnt
      integer, intent(out) ::  toklen(*)
      
      character(len=*), intent(in)  :: string
      character(len=*), intent(out) :: tokens(*)
      character(len=1), intent(in)  :: delim
      
!-----------------------------------------------------------------------     
!	... Local variables
!-----------------------------------------------------------------------     
      integer  ::   marker, i, length, endpos

      tokcnt = 0
      marker = 1
      do i = 1,ls
         if( string(i:i) == delim .or. i == ls ) then
            if( i == ls ) then
               if( string(i:i) == delim ) then
                  tokcnt = 0
                  exit
               end if
               length = i - marker + 1
               endpos = i
            else
               length = i - marker
               endpos = i - 1
            end if
            if( length < 1 .or. length > maxlen ) then
               tokcnt = 0
               exit
            end if
            tokcnt = tokcnt + 1
            if( tokcnt > maxtok ) then
               tokcnt = -tokcnt
               exit
            end if
            tokens(tokcnt) = ' '
            tokens(tokcnt)(:length) = string(marker:endpos)
            toklen(tokcnt) = length
            marker = i + 1
         end if
      end do

      if( tokcnt < 1 ) then
        write(*,*) 'gettokens: input string'
        write(*,*) trim(string)
        write(*,*) '           has syntax error(s)'
        stop 'Syntax err'
      endif
      
      end subroutine gettokens


   subroutine get_indices( match_strng, string, match_cnt, match_ndx )
!------------------------------------------------------------------
!     get all indices of match char in string
!------------------------------------------------------------------

!------------------------------------------------------------------
!     dummy arguments
!------------------------------------------------------------------
   integer, intent(out) :: match_cnt
   integer, optional, intent(out) :: match_ndx(:)
   character(len=*), intent(in) :: match_strng
   character(len=*), intent(in) :: string

!------------------------------------------------------------------
!     local variables
!------------------------------------------------------------------
   integer :: pos, ndx
   logical :: set_ndx

   if( present( match_ndx ) ) then
     set_ndx = .true.
   else
     set_ndx = .false.
   endif
   ndx       = 0
   match_cnt = 0
   do
     pos = index( string(ndx+1:), match_strng )
     if( pos == 0 ) then
       exit
     endif
     ndx       = ndx + pos 
     match_cnt = match_cnt + 1 
     if( set_ndx ) then
       match_ndx(match_cnt) = ndx
     endif
   end do

   end subroutine get_indices

   real function get_number( string )
!------------------------------------------------------------------
!  encode incoming string into a number
!------------------------------------------------------------------

!------------------------------------------------------------------
!  dummy arguments
!------------------------------------------------------------------
   character(len=*), intent(in) :: string

!------------------------------------------------------------------
!  local variables
!------------------------------------------------------------------
   integer :: ios, slen
   real    :: coeff

   slen = len_trim( string )
   read(string(:slen),fmt=*,iostat=ios) coeff
   if( ios /= 0 ) then
     write(*,*) 'get_number: ',trim(string),' is not a valid number'
     stop 'Syntax err'
   endif

   get_number = coeff

   end function get_number

   integer function get_src_ndx( src )
!------------------------------------------------------------------
!  retrieve the sub-category index
!------------------------------------------------------------------

!------------------------------------------------------------------
!  dummy arguments
!------------------------------------------------------------------
   character(len=*), intent(in) :: src

!------------------------------------------------------------------
!  local variables
!------------------------------------------------------------------
   integer :: m

   do m = 1,nsrc_names
     if( trim(src_names(m)) == trim(src) ) then
       exit       
     endif
   end do

   if( m > nsrc_names ) then
     write(*,*) 'get_src_ndx: ',trim(src),' is not a valid source'
     stop 'Syntax err'
   endif

   get_src_ndx = m

   end function get_src_ndx

   integer function get_cat_ndx( cat )
!------------------------------------------------------------------
!  retrieve the sub-category index
!------------------------------------------------------------------

!------------------------------------------------------------------
!  dummy arguments
!------------------------------------------------------------------
   character(len=*), intent(in) :: cat

!------------------------------------------------------------------
!  local variables
!------------------------------------------------------------------
   integer :: m

   do m = 1,n_sub_cats
     if( trim(sub_cats(m)) == trim(cat) ) then
       exit       
     endif
   end do

   if( m > n_sub_cats ) then
     write(*,*) 'get_cat_ndx: ',trim(cat),' is not a valid sub-category'
     stop 'Syntax err'
   endif

   get_cat_ndx = m

   end function get_cat_ndx
 
   subroutine remove_blanks( string )

   character(len=*), intent(inout) :: string

   integer :: i, j
   character(len=linsize) :: tstrng

   j = 0
   tstrng = ' '
   do i = 1,len_trim(string)
     if( string(i:i) /= ' ' ) then
       j = j + 1     
       tstrng(j:j) = string(i:i)
     endif
   end do

   string = trim(tstrng)

   end subroutine remove_blanks

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

   subroutine upcase( in_string, out_string )
!------------------------------------------------------------------
!  convert string to upper case
!------------------------------------------------------------------

!------------------------------------------------------------------
!  dummy arguments
!------------------------------------------------------------------
   character(len=*), intent(in)  :: in_string
   character(len=*), intent(out) :: out_string

!------------------------------------------------------------------
!  local variables
!------------------------------------------------------------------
   integer :: incr
   integer :: m

   incr = ichar( 'A' ) - ichar( 'a' )
   out_string(:) = ' '
   do m = 1,min(len(in_string),len(out_string))
     if( lge( in_string(m:m),'a' ) .and. lle( in_string(m:m),'z' ) ) then
       out_string(m:m) = char( ichar( in_string(m:m) ) + incr )
     else
       out_string(m:m) = in_string(m:m)
     endif
   end do

   end subroutine upcase

   subroutine cleanup_grid( grid, ide, jde )
!------------------------------------------------------------------
!  cleanup allocated variables in grid structures
!------------------------------------------------------------------

   use mapper_types

!------------------------------------------------------------------
!  dummy arguments
!------------------------------------------------------------------
   integer, intent(in) :: ide
   integer, intent(in) :: jde
   type(grid_type), intent(inout) :: grid

!------------------------------------------------------------------
!  local variables
!------------------------------------------------------------------
   integer :: i, j

   if( allocated( grid%ix ) ) then
     deallocate( grid%ix, grid%jy, grid%ax, grid%by )
   endif
   if( allocated( grid%lon ) ) then
     deallocate( grid%lon, grid%lat, grid%xedge, grid%yedge )
   endif
   if( allocated( grid%model_area_type ) ) then
     do j = 1,jde
       do i = 1,ide
         if( allocated( grid%model_area_type(i,j)%dcell_lon_ndx ) ) then
           deallocate( grid%model_area_type(i,j)%dcell_lon_ndx )
           deallocate( grid%model_area_type(i,j)%dcell_lat_ndx )
         endif
         if( allocated( grid%model_area_type(i,j)%wght ) ) then
           deallocate( grid%model_area_type(i,j)%wght )
         endif
       end do
     end do
     deallocate( grid%model_area_type )
   endif

   end subroutine cleanup_grid

   end module utils
