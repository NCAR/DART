! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download


!>  This program takes as input one or more sets of N obs_seq filenames. 
!>  N is usually 2, but can be any number of files to compare together.
!>  It opens the obs_seq files N at a time and does an obs-by-obs comparison.  
!>  It is intended to pull out a common subset of obs which were successfully
!>  assimilated in all experiments so the subsequent diagnostics can give
!>  an apples-to-apples comparison of the relative error or fit to obs
!>  between different experiments.
!> 
!>  This tool does not search in the input obs_seq files to try to match obs
!>  in different orders or to skip over obs present in only one input file,
!>  so all experiments should be done with identical input obs_seq files.
!>  The obs to be assimilated or evaluated can be selected by namelist
!>  at run time if they are different, and different experiment settings
!>  can be selected by other namelist values.
!> 
!>  If the next input obs don't match exactly (identical types, times, locations,
!>  and DART QC) they are discarded.  If they match, and if the DART QC indicates 
!>  they have been successfully assimilated in all experiments, then each of the 
!>  N input obs are copied through to the N output obs_seq files.  At the end of 
!>  running this tool it creates sets of N output obs_seq files which contain a 
!>  consistent set of assimilated obs across the set.
!> 
!>  The outputs are suitable for input to the obs_diag or obs_seq_to_netcdf tools
!>  and there are matlab scripts for comparing the results of N experiments.
!>  Differences in the output diagnostics will now be because of the experiment 
!>  differences, not differences in the number of obs assimilated.
!> 
!>  Running this program creates sets of N new output files corresponding to each
!>  of N input files.  The output names are the input filenames with a suffix appended.  
!>  It can take a list of obs_seq files, either given directly in the namelist or
!>  a list of N filenames where each file contains one input filename per line.  
!>  If the lists contain more than 1 filename, this program will loop over each
!>  set of filenames in the lists until the end of the file.If more than N files
!> 
!>  There is an option to allow observations which have been either successfully 
!>  assimilated or evaluated to match each other.  The default requires them to
!>  have been treated identically in all experiments.

program obs_common_subset

use        types_mod, only : r8, metadatalength
use    utilities_mod, only : initialize_utilities,            &
                             find_namelist_in_file, check_namelist_read,       &
                             error_handler, E_ERR, E_MSG, nmlfileunit,         &
                             do_nml_file, do_nml_term, get_next_filename,      &
                             finalize_utilities, logfileunit
use     location_mod, only : location_type, operator(/=)
use      obs_def_mod, only : obs_def_type, get_obs_def_time, get_obs_def_type_of_obs,     &
                             get_obs_def_location
use     obs_kind_mod, only : max_defined_types_of_obs, get_name_for_type_of_obs
use time_manager_mod, only : time_type, print_date, print_time, set_time,      &
                             set_calendar_type, get_calendar_type,             &
                             operator(/=), operator(>), NO_CALENDAR
use obs_sequence_mod, only : obs_sequence_type, obs_type, write_obs_seq,       &
                             init_obs, assignment(=), get_obs_def, get_obs_key,&
                             init_obs_sequence, static_init_obs_sequence,      &
                             read_obs_seq_header, read_obs_seq, get_num_obs,   &
                             get_first_obs, get_next_obs, get_num_key_range,   &
                             insert_obs_in_seq, get_num_copies, get_num_qc,    &
                             get_copy_meta_data, get_qc_meta_data, get_qc,     &
                             set_copy_meta_data, set_qc_meta_data,             &
                             destroy_obs, destroy_obs_sequence

implicit none

character(len=*), parameter :: source = 'obs_common_subset.f90'

! max_num_input_files : maximum total number of input sequence files to be processed.
integer, parameter   :: max_num_input_files = 5000

! max number of files to compare against each other.  if you need to compare 
! more than 50 files together, let us know what you are doing because i'd be 
! very interested to hear.  but you can make the maxcomp larger in that case.

integer, parameter      :: maxcomp = 50  ! max number of files compared at once

type(obs_sequence_type) :: seq_in(maxcomp)
type(obs_sequence_type) :: seq_out(maxcomp)
type(obs_type)          :: obs_in(maxcomp), next_obs_in(maxcomp)
type(obs_type)          :: obs_out(maxcomp), prev_obs_out(maxcomp)
logical                 :: is_this_last(maxcomp)
logical                 :: wanted
integer                 :: size_seq_in(maxcomp), num_copies_in(maxcomp), num_qc_in(maxcomp)
integer                 :: size_seq_out, num_inserted, iunit, io, i, j, k, nextfile
integer                 :: max_num_obs(maxcomp), file_id, atonce, nsets
integer                 :: num_rejected_badqc, num_rejected_diffqc, num_rejected_other
integer                 :: num_mismatch_loc, num_mismatch_time, num_mismatch_type
character(len=129)      :: read_format
logical                 :: pre_I_format, cal
character(len=512)      :: msgstring, msgstring1, msgstring2, msgstring3
character(len=256)      :: filename_in(maxcomp)
character(len=384)      :: filename_out(maxcomp)     ! filename_in + . + suffix

character(len=metadatalength), parameter :: dart_qc_meta_data = 'DART quality control'
character(len=metadatalength) :: meta_data

integer                 :: qc_index
integer                 :: num_input_sets = 0

character(len=256) :: temp_filelist(max_num_input_files) = ''

!----------------------------------------------------------------
! Namelist input with default values

character(len=256) :: filename_seq(max_num_input_files) = ''
character(len=256) :: filename_seq_list(maxcomp)  = ''
character(len=32)  :: filename_out_suffix = '.common'

integer :: num_to_compare_at_once = 2
logical :: print_only = .false.

! not expected to be changed often, but its here if needed.
integer :: print_every = 10000
integer :: dart_qc_threshold = 3     ! ok if DART QC <= this
logical :: eval_and_assim_can_match = .false. 

character(len=32) :: calendar = 'Gregorian'

namelist /obs_common_subset_nml/ &
         num_to_compare_at_once,          &
         filename_seq, filename_seq_list, &
         filename_out_suffix, print_only, &
         print_every, dart_qc_threshold,  &
         calendar, eval_and_assim_can_match

!----------------------------------------------------------------
! Start of the program:
!
! Process each N input observation sequence file in turn,
! selecting observations to insert into the output sequence files.
!----------------------------------------------------------------

call setup()

! Read the namelist entry
call find_namelist_in_file("input.nml", "obs_common_subset_nml", iunit)
read(iunit, nml = obs_common_subset_nml, iostat = io)
call check_namelist_read(iunit, io, "obs_common_subset_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=obs_common_subset_nml)
if (do_nml_term()) write(     *     , nml=obs_common_subset_nml)

! the logic here is slightly different than the obs_sequence_tool.
! the user gives us a list of obs_seq files; either in the namelist
! or as the name of a file which contains the list, one per line.
! the files are compared N at a time, until the list is exhausted.
! as in the other tools, it is an error to specify both an explicit
! list and the name of file for input.

write(msgstring, '(A,I3,A)') "Asking to compare ", num_to_compare_at_once, " obs_seq files at a time"
call error_handler(E_MSG, "obs_common_subset", msgstring, source)

if ((num_to_compare_at_once < 2) .or. (num_to_compare_at_once > maxcomp)) then
   write(msgstring, *) "num_to_compare_at_once must be >= 2 and <= ", maxcomp
   call error_handler(E_ERR, "obs_common_subset", msgstring, source)
endif

! if there is only a single file from each experiment to compare, the order
! doesn't matter.  if each experiment has a list of input obs_seq.final files from
! a multi-step assimilation, list all the files from experiment 1 in list-file 1,
! all files from experiment 2 in list-file 2, etc.
call handle_filenames(num_to_compare_at_once, maxcomp, filename_seq, filename_seq_list, &
                      num_input_sets)

if (num_input_sets > 1) then
   write(msgstring, '(A,I3,A)') "Will loop ", num_input_sets, " times to process the entire input filelist"
   call error_handler(E_MSG, "obs_common_subset", msgstring, source)
endif

! because i am a lazy typist
atonce = num_to_compare_at_once
nsets = num_input_sets

! the default is a gregorian calendar.  if you are using a different type
! set it in the namelist.  this only controls how it prints out the first
! and last timestamps in the obs_seq files.
call set_calendar_type(calendar)

! set a logial to see if we have a calendar or not
cal = (get_calendar_type() /= NO_CALENDAR)

! end of namelist processing and setup


! process the files in sets, making sure the metadata matches. keep matching
! obs which share the same QC values, as long as the QCs are <= threshold.
! since this tool has very little utility without the DART QC values, fail
! if they aren't present.

! count of input sets was set in the handle_filenames() routine above.
nextfile = 1
NUMSETS: do j = 1, nsets

   ! if num_to_compare_at_once = 3 and num_input_sets = 100, then this
   ! assumes the input file list is organized in this order:
   !       a1,b1,c1, a2,b2,c2, ..., a100,b100,c100
   ! where a, b, c are the 3 separate experiments to compare, 
   ! and 1, 2, ... 100 are the obs_seq.final output files from 100 steps 
   ! of a multi-step assimilation.  handle_filenames() arranges this
   ! order if list files were used.

   ! fill in the names of the current set of obs_seq files
   do i = 1, atonce
      filename_in(i) = filename_seq(nextfile)
      nextfile = nextfile + 1
   enddo

   ! now filename_in(:) has all the names for this set.  read the headers in.
   do i = 1, atonce
  
      call read_obs_seq_header(filename_in(i), num_copies_in(i), num_qc_in(i), &
         size_seq_in(i), max_num_obs(i), file_id, read_format, pre_I_format, &
         close_the_file = .true.)

      if (size_seq_in(i) == 0) then
         write(msgstring1, *) 'This tool cannot process empty obs_seq files'
         write(msgstring,*) 'No obs found in input sequence file ', trim(filename_in(i))
         call error_handler(E_ERR,'obs_common_subset',msgstring, &
            source, text2=msgstring1)
      endif

   enddo

   do i = 2, atonce
      if (size_seq_in(1) /= size_seq_in(i)) then
         write(msgstring,*) 'Input obs_seq files cannot have different observation counts'
         write(msgstring2,*) 'Counts for file 1 and file ', i,' are ', size_seq_in(1), &
                             ' obs and ', size_seq_in(i), ' obs'
         call error_handler(E_ERR,'obs_common_subset',msgstring, &
            source, text2=msgstring2)
      endif
   enddo

   write(     *,      *) 'Starting to process input sequence files: '
   write(logfileunit, *) 'Starting to process input sequence files: '

   do i = 1, atonce
      write(     *     , *)  trim(filename_in(i))
      write(logfileunit, *)  trim(filename_in(i))

      call read_obs_seq(filename_in(i), 0, 0, 0, seq_in(i))
   enddo

   ! make sure the files have compatible metadata.  this errors out here if not.
   call compare_metadata(atonce, seq_in, filename_in)

   ! sanity check - ensure the linked list times are in increasing time order
   do i = 1, atonce
      call validate_obs_seq_time(seq_in(i), filename_in(i))
   enddo

   ! the output can have no more than the number of input obs, and we already
   ! tested that all sequences have the same obs count.
   size_seq_out = size_seq_in(1)

   ! find which DART QC copy has the DART quality control.  this version of
   ! the tool fails if one is not found.  there seems very little utility in
   ! this tool if it tries to process files without a DART QC.
   qc_index = -1
   QCLOOP: do i=1, get_num_qc(seq_in(1))
      if(index(get_qc_meta_data(seq_in(1),i), dart_qc_meta_data) > 0) then
         qc_index = i
         exit QCLOOP
      endif
   enddo QCLOOP

   ! if we do decide there is utility in running these on files without a
   ! DART QC, remove this error exit here.   nsc 24 oct 2013
   if (qc_index < 0) then
      write(msgstring,*)  'DART QC was not found in the input obs_seq files.'
      write(msgstring2,*) 'This tool only works on output files from an assimilation.'
      call error_handler(E_ERR,'obs_common_subset',msgstring, &
            source, text2=msgstring2)
   endif


   ! blank line, start of actually creating output file
   call error_handler(E_MSG,' ',' ')

   do i = 1, atonce
      ! Initialize individual observation variables
      call init_obs(      obs_in(i), num_copies_in(1), num_qc_in(1))
      call init_obs( next_obs_in(i), num_copies_in(1), num_qc_in(1))
      call init_obs(     obs_out(i), num_copies_in(1), num_qc_in(1))
      call init_obs(prev_obs_out(i), num_copies_in(1), num_qc_in(1))

      ! create the output sequences
      call init_obs_sequence(seq_out(i), num_copies_in(1), num_qc_in(1), size_seq_out)

      do k=1, num_copies_in(1)
         meta_data = get_copy_meta_data(seq_in(1), k)
         call set_copy_meta_data(seq_out(i), k, meta_data)
      enddo
      do k=1, num_qc_in(1)
         meta_data = get_qc_meta_data(seq_in(1), k)
         call set_qc_meta_data(seq_out(i), k, meta_data)
      enddo

   enddo

   ! Start to insert obs from sequence_in into sequence_out
   ! NOTE: insert_obs_in_seq CHANGES the obs passed in.
   !       Must pass a copy of incoming obs to insert_obs_in_seq.

   num_inserted = 0
   num_rejected_badqc = 0
   num_rejected_diffqc = 0
   num_rejected_other = 0
   num_mismatch_time = 0
   num_mismatch_type = 0
   num_mismatch_loc = 0

   ! get the first obs from each file.  since we have already checked
   ! for empty sequence files, this is not expected to fail at this point.
   do i = 1, atonce
      if ( .not. get_first_obs(seq_in(i), obs_in(i))) then
         write(msgstring, *)'no first observation in ',trim(filename_in(i))
         call error_handler(E_ERR,'obs_common_subset',msgstring,source)
      endif
      is_this_last(i) = .false.
      next_obs_in(i) = obs_in(i)
   enddo


   ! since we've checked to be sure that all the input obs_seq files have
   ! the same obs counts, i can test the result from the first seq only.
   ObsLoop : do while ( .not. is_this_last(1) )

      ! move the next obs into position to be worked on
      do i = 1, atonce
         obs_in(i) = next_obs_in(i)
      enddo

      ! see if we want to keep it or not
      wanted = all_good(obs_in, atonce, qc_index, dart_qc_threshold)

      if ( wanted ) then
         do i = 1, atonce
            obs_out(i) = obs_in(i)
         enddo

         ! Since the stride through the observation sequence file is always
         ! guaranteed to be in temporally-ascending order, we can use the
         ! 'previous' observation as the starting point to search for the
         ! correct insertion point.  This speeds up the insert code a lot.

         if (num_inserted == 0) then
            do i = 1, atonce
               call insert_obs_in_seq(seq_out(i), obs_out(i))
            enddo
         else
            do i = 1, atonce
               call insert_obs_in_seq(seq_out(i), obs_out(i), prev_obs_out(i))
            enddo
         endif

         ! update position in seq for next insert
         do i = 1, atonce
            prev_obs_out(i) = obs_out(i)
         enddo

         num_inserted = num_inserted + 1

         if (print_every > 0) then
            if (mod(num_inserted,print_every) == 0) then
               print*, 'inserted obs number ',num_inserted,' of ',size_seq_out
            endif
         endif

      endif

      do i = 1, atonce
         call get_next_obs(seq_in(i), obs_in(i), next_obs_in(i), is_this_last(i))
      enddo

   enddo ObsLoop


   ! finished creating output sequences.  print info and write out files

   write(msgstring, *) 'Number of obs in each output seq file :', get_num_key_range(seq_out(1))
   call error_handler(E_MSG,'obs_common_subset',msgstring)

   do i = 1, atonce
      filename_out(i) = trim(filename_in(i))//trim(filename_out_suffix)
      call print_obs_seq(seq_out(i), filename_out(i))
   enddo

   if (print_only)  then
      write(msgstring,*) 'Output sequence files not created; print_only in namelist is .true.'
      call error_handler(E_MSG,'', msgstring)
   else

      print*, '---------  Obs seq file set  #  :         ', j
      print*, 'Number of obs in each input     :         ', size_seq_in(1)
      print*, 'Number of obs rejected diff qc  :         ', num_rejected_diffqc
      print*, 'Number of obs rejected bad qc   :         ', num_rejected_badqc
      print*, 'Number of obs rejected other    :         ', num_rejected_other
      if (num_rejected_other > 0) then
      print*, ' Number of mismatches in time   :         ', num_mismatch_time
      print*, ' Number of mismatches in type   :         ', num_mismatch_type
      print*, ' Number of mismatches in loc    :         ', num_mismatch_loc
      endif
      print*, 'Number of obs copied to output  :         ', num_inserted
      print*, '---------------------------------------------------------'

      write(msgstring, *) 'Starting to write output sequence files'
      call error_handler(E_MSG,'obs_common_subset',msgstring)

      do i = 1, atonce
         call write_obs_seq(seq_out(i), filename_out(i))
      enddo
   endif

   ! clean up - prev_obs_out is a copy of obs_out, so don't call destroy on it
   ! because it's already been deleted.

   do i = 1, atonce
      call destroy_obs_sequence(seq_in(i))
      call destroy_obs_sequence(seq_out(i))
      call destroy_obs(obs_in(i))
      call destroy_obs(obs_out(i))
      call destroy_obs(next_obs_in(i))
   enddo

enddo NUMSETS

call shutdown()

!---------------------------------------------------------------------
! end of main program.
!---------------------------------------------------------------------


contains


!---------------------------------------------------------------------

subroutine setup()

! Initialize modules used that require it
call initialize_utilities('obs_common_subset')
call static_init_obs_sequence()

end subroutine setup

!---------------------------------------------------------------------

subroutine shutdown()

call finalize_utilities('obs_common_subset')

end subroutine shutdown

!---------------------------------------------------------------------

subroutine handle_filenames(setsize, maxsets, filename_seq, filename_seq_list, &
                            num_input_sets)

! sort out the input lists, set the length as a return in num_input_sets,
! fill in filename_seq if list was given, and make sure what's specified is consistent.

integer,            intent(in)    :: setsize
integer,            intent(in)    :: maxsets
character(len=*),   intent(inout) :: filename_seq(:)
character(len=*),   intent(in)    :: filename_seq_list(:)
integer,            intent(out)   :: num_input_sets

integer :: indx, num_input_files, i, ifiles, fcount, this_set_filecount
integer :: files_per_set, nlists, iset, src, dst
logical :: from_file
character(len=32) :: nname

! if the user specifies neither filename_seq nor filename_seq_list, error
! if the user specifies both, error.
! if list length isn't a whole multiple of num_to_compare_at_once, error.
! if user gives a filelist name, make sure the length is not more than 
! maxfiles and read it into the explicit list and continue.
! set num_input_sets to the count of sets in the list

if (filename_seq(1) == '' .and. filename_seq_list(1) == '') then
   msgstring2='One of filename_seq or filename_seq_list must be specified in namelist'
   call error_handler(E_ERR,'handle_filenames',            &
           'no filenames specified as input', source, text2=msgstring2)
endif

! make sure the namelist specifies one or the other but not both
if (filename_seq(1) /= '' .and. filename_seq_list(1) /= '') then
   call error_handler(E_ERR,'handle_filenames', &
       'cannot specify both filename_seq and filename_seq_list in namelist', source)
endif

! count of all filenames, either explicitly in the namelist
! or specified by a list of files which contain lists of names.
num_input_files = -1
nlists = -1

if (filename_seq(1) /= '') then

   ! if the files are given in the namelist, just need to make sure
   ! there are at least N given (where N == num_to_compare_at_once).

   nname = 'filename_seq'
   from_file = .false.

   FILELOOP: do indx = 1, max_num_input_files
      ! an empty name ends the list
      if (filename_seq(indx) == '') then
         if (indx < num_to_compare_at_once) then
            write(msgstring, '(A,I5,A)') trim(nname)//' must contain at least ', &
                num_to_compare_at_once, ' obs_seq filenames'
            call error_handler(E_ERR,'handle_filenames', msgstring, source)
         endif
         num_input_files = indx - 1 
         exit FILELOOP
      endif
   enddo FILELOOP

   ! make sure the list is an even multiple of setsize
   if (modulo(num_input_files, setsize) /= 0) then
      write(msgstring, *) 'number of input files must be an even multiple of ', setsize
      write(msgstring1, *) 'found ', num_input_files, ' input files'
      call error_handler(E_ERR,'handle_filenames', msgstring, source, text2=msgstring1)
   endif

   num_input_sets = num_input_files / setsize
else

   ! if they have specified a file which contains a list, read it into
   ! a temp array, set the counts, and do some initial error checks.

   nname = 'filename_seq_list'
   from_file = .true.

   ! this is the total running count of files across all lists
   fcount = 1   

   EXPLOOP: do iset = 1, maxsets
      
      ! an empty name ends the list of lists
      if (filename_seq_list(iset) == '') then
         num_input_files = fcount - 1 
         nlists = iset - 1
         if (nlists /= setsize) then
            write(msgstring1, '(A)') 'number of input file lists in "'//trim(nname)// &
                                     '" must equal "num_to_compare_at_once"'
            write(msgstring2, '(2(A,I5))') trim(nname)//' has ', nlists, &
                  ' items while num_to_compare_at_once is ', setsize
         
            call error_handler(E_ERR,'handle_filenames', msgstring1, &
                               source, text2=msgstring2)
         endif
         exit EXPLOOP
      endif
   
      LISTLOOP: do ifiles=1, max_num_input_files
         if (fcount >= max_num_input_files) then
            write(msgstring, *)  'cannot specify more than ',max_num_input_files,' total files'
            write(msgstring1, *) 'for more, change "max_num_input_files" in the source and recompile'
            call error_handler(E_ERR,'handle_filenames', msgstring, source, &
                               text2=msgstring1)
         endif

         temp_filelist(fcount) = get_next_filename(filename_seq_list(iset), ifiles)
     
         ! an empty name ends the list
         if (temp_filelist(fcount) == '') then
            if (ifiles == 1) then
               call error_handler(E_ERR,'handle_filenames', &
                   trim(filename_seq_list(iset))//' contains no input obs_seq filenames', &
                   source)
            endif
            this_set_filecount = ifiles-1
            exit LISTLOOP
         endif
         fcount = fcount + 1
      enddo LISTLOOP

      ! make sure each file in the list contains the same number
      ! of filenames to be compared.  set the target the first time
      ! through and then compare subsequent counts to be sure they match.
      if (iset == 1) then
         files_per_set = this_set_filecount
      else
         if (files_per_set /= this_set_filecount) then
            write(msgstring1, '(A)') trim(filename_seq_list(iset))//' does not contain ' // &
               'the same number of obs_seq filenames as previous lists'
            write(msgstring2, '(3(A,I5),A)') 'list 1 contains ', files_per_set, ' filenames while list ', iset, &
               ' contains ', this_set_filecount, ' filenames'
            call error_handler(E_ERR,'handle_filenames', msgstring1, &
                    source, text2=msgstring2)
         endif
      endif
   
   enddo EXPLOOP
   
   num_input_sets = files_per_set
endif

if (num_input_files >= max_num_input_files) then
   write(msgstring, *)  'cannot specify more than ',max_num_input_files,' files'
   write(msgstring1, *) 'for more, change "max_num_input_files" in the source and recompile'
   call error_handler(E_ERR,'handle_filenames', msgstring, source, text2=msgstring1)
endif

if (from_file) then
   ! now that we know we have the right number of filenames, copy them
   ! from the temp list into filename_seq in the right order for processing:  
   ! first all file 1s from each experiment, then file 2s, etc.
   ! the temp list is all files from exp1, all files from exp2, etc
   dst = 0
   do i=1, num_input_sets
      do j=1, setsize
         src = (j-1)*num_input_sets + i
         dst = dst + 1
         filename_seq(dst) = temp_filelist(src)
      enddo
   enddo
endif

end subroutine handle_filenames

!---------------------------------------------------------------------

subroutine compare_metadata(count, seq, fnames)

! This subroutine compares the metadata for 'count' different observation
! sequences and terminates the program if they are not conformable.
! In order to be processed all observation sequences must have the same
! number of qc values, the same number of copies, and the same metadata names.

integer,                 intent(in) :: count
type(obs_sequence_type), intent(in) :: seq(:)
character(len=*),        intent(in) :: fnames(:)

integer :: num_copies1, num_copiesN, num_qc1, num_qcN
integer :: i, j
character(len=metadatalength) :: str1, strN

num_copies1 = get_num_copies(seq(1))
num_qc1 = get_num_qc(seq(1))

do i = 2, count
   num_copiesN = get_num_copies(seq(i))
   num_qcN = get_num_qc(seq(i))

   ! get this ready in case we have to use it in some later error string.
   write(msgstring3,*)'Sequence files ', trim(fnames(1)), ' and ', trim(fnames(i)), &
                      ' are not compatible'

   ! compare the counts first - compare all against sequence 1 because they
   ! all have to be identical to succeed.
   if ( num_copies1 /= num_copiesN ) then
      write(msgstring2,*)'Different numbers of data copies found: ', &
                         num_copies1, ' vs ', num_copiesN
      num_copies1 = -1
   endif
   if ( num_qc1 /= num_qcN ) then
      write(msgstring2,*)'Different different numbers of QCs found: ', &
                         num_qc1, ' vs ', num_qcN
      num_qc1 = -1
   endif
   if ( num_copies1 < 0 .or. num_qc1 < 0 ) then
      call error_handler(E_ERR, 'compare_metadata', msgstring3, &
                                 source, text2=msgstring2)
   endif

   ! make sure the metadata matches exactly for both the copies and qcs
   ! in all files.
   CopyMetaData : do j=1, num_copies1
      str1 = get_copy_meta_data(seq(1),j)
      strN = get_copy_meta_data(seq(i),j)

      ! if they match, continue.
      ! if they don't match, it's a fatal error.
      if( str1 /= strN ) then
         write(msgstring1,*)'copy metadata value mismatch. ', trim(str1)
         write(msgstring2,*)'copy metadata value mismatch. ', trim(strN)
         call error_handler(E_ERR, 'compare_metadata', msgstring3, &
                    source, text2=msgstring1, text3=msgstring2)
      endif

   enddo CopyMetaData

   QCMetaData : do j=1, num_qc1
      str1 = get_qc_meta_data(seq(1),j)
      strN = get_qc_meta_data(seq(i),j)

      ! if they match, continue.
      ! if they don't match, it's a fatal error.
      if( str1 /= strN ) then
         write(msgstring1,*)'qc metadata value mismatch.  ', trim(str1)
         write(msgstring2,*)'qc metadata value mismatch.  ', trim(strN)
         call error_handler(E_ERR, 'compare_metadata', msgstring3, &
                    source, text2=msgstring1, text3=msgstring2)
      endif

   enddo QCMetaData

enddo

end subroutine compare_metadata

!---------------------------------------------------------------------

subroutine print_obs_seq(seq_in, filename)

! you can get more info by running the obs_diag program, but this
! prints out a quick table of obs types and counts, overall start and
! stop times, and metadata strings and counts.

type(obs_sequence_type), intent(in) :: seq_in
character(len=*),        intent(in) :: filename

type(obs_type)          :: obs, next_obs
type(obs_def_type)      :: this_obs_def
logical                 :: is_there_one, is_this_last
integer                 :: size_seq_in
integer                 :: i
integer                 :: this_obs_kind
! max_defined_types_of_obs is a public from obs_kind_mod.f90 and really is
! counting the max number of types, not kinds
integer                 :: type_count(max_defined_types_of_obs), identity_count


! Initialize input obs_types
do i = 1, max_defined_types_of_obs
   type_count(i) = 0
enddo
identity_count = 0

! make sure there are obs left to process before going on.
! num_obs should be ok since we just constructed this seq so it should
! have no unlinked obs.  if it might for some reason, use this instead:
! size_seq_in = get_num_key_range(seq_in)     !current size of seq_in

size_seq_in = get_num_obs(seq_in)
if (size_seq_in == 0) then
   msgstring = 'Obs_seq file '//trim(filename)//' is empty.'
   call error_handler(E_MSG,'print_obs_seq',msgstring)
   return
endif

! Initialize individual observation variables
call init_obs(     obs, get_num_copies(seq_in), get_num_qc(seq_in))
call init_obs(next_obs, get_num_copies(seq_in), get_num_qc(seq_in))

! blank line
call error_handler(E_MSG,'',' ')

write(msgstring,*) 'Processing sequence file ', trim(filename)
call error_handler(E_MSG,'',msgstring)

!  for an obs_seq final file which has lots of ensemble members this
!  is really long and obscures useful data.  disable it for now.
!  at some point we could make a verbose option which could re-enable it.
if (.false.) call print_metadata(seq_in, filename)


! Start to process obs from seq_in

is_there_one = get_first_obs(seq_in, obs)

if ( .not. is_there_one )  then
   write(msgstring,*)'no first observation in ',trim(filename)
   call error_handler(E_MSG,'print_obs_seq', msgstring)
endif

! process it here
is_this_last = .false.

call get_obs_def(obs, this_obs_def)
call print_time(get_obs_def_time(this_obs_def), ' First timestamp: ')
! does not work with NO_CALENDAR
if (cal) call print_date(get_obs_def_time(this_obs_def), '   calendar Date: ')

ObsLoop : do while ( .not. is_this_last)

   call get_obs_def(obs, this_obs_def)
   this_obs_kind = get_obs_def_type_of_obs(this_obs_def)
   if (this_obs_kind < 0) then
      identity_count = identity_count + 1
   else
      type_count(this_obs_kind) = type_count(this_obs_kind) + 1
   endif

   call get_next_obs(seq_in, obs, next_obs, is_this_last)
   if (.not. is_this_last) then
      obs = next_obs
   else
      call print_time(get_obs_def_time(this_obs_def), '  Last timestamp: ')
      if (cal) call print_date(get_obs_def_time(this_obs_def), '   calendar Date: ')
   endif

enddo ObsLoop


write(msgstring, *) 'Number of obs processed  :          ', size_seq_in
call error_handler(E_MSG, '', msgstring)
write(msgstring, *) '---------------------------------------------------------'
call error_handler(E_MSG, '', msgstring)
do i = 1, max_defined_types_of_obs
   if (type_count(i) > 0) then
      write(msgstring, '(a32,i8,a)') trim(get_name_for_type_of_obs(i)), &
                                     type_count(i), ' obs'
      call error_handler(E_MSG, '', msgstring)
   endif
enddo
if (identity_count > 0) then
   write(msgstring, '(a32,i8,a)') 'Identity observations', &
                                  identity_count, ' obs'
   call error_handler(E_MSG, '', msgstring)
endif

! another blank line
call error_handler(E_MSG, '', ' ')

! Time to clean up

call destroy_obs(     obs)
call destroy_obs(next_obs)

end subroutine print_obs_seq

!---------------------------------------------------------------------

subroutine validate_obs_seq_time(seq, filename)

! this eventually belongs in the obs_seq_mod code, but for now
! try it out here.  we just fixed a hole in the interactive create
! routine which would silently let you create out-of-time-order
! linked lists, which gave no errors but didn't assimilate the
! right obs at the right time when running filter.   this runs
! through the times in the entire sequence, ensuring they are
! monotonically increasing in time.  this should help catch any
! bad files which were created with older versions of code.

type(obs_sequence_type), intent(in) :: seq
character(len=*),        intent(in) :: filename

type(obs_type)          :: obs, next_obs
type(obs_def_type)      :: this_obs_def
logical                 :: is_there_one, is_this_last
integer                 :: size_seq, obs_count
integer                 :: key
type(time_type)         :: last_time, this_time


! make sure there are obs left to process before going on.
size_seq = get_num_obs(seq)
if (size_seq == 0) then
   msgstring = 'Obs_seq file '//trim(filename)//' is empty.'
   call error_handler(E_MSG,'validate_obs_seq_time',msgstring)
   return
endif

! Initialize individual observation variables
call init_obs(     obs, get_num_copies(seq), get_num_qc(seq))
call init_obs(next_obs, get_num_copies(seq), get_num_qc(seq))

obs_count = 0


! Start to process obs from seq

is_there_one = get_first_obs(seq, obs)

! we already tested for 0 obs above, so there should be a first obs here.
if ( .not. is_there_one )  then
   write(msgstring,*)'no first obs in sequence ' // trim(filename)
   call error_handler(E_ERR,'validate_obs_seq_time', &
                      msgstring, source)
   return
endif

is_this_last = .false.
last_time = set_time(0, 0)
ObsLoop : do while ( .not. is_this_last)

   call get_obs_def(obs, this_obs_def)
   this_time = get_obs_def_time(this_obs_def)

   if (last_time > this_time) then
      ! bad time order of observations in linked list
      call print_time(last_time, ' previous timestamp: ')
      if (cal) call print_date(last_time, '      calendar date: ')
      call print_time(this_time, '     next timestamp: ')
      if (cal) call print_date(this_time, '      calendar date: ')

      key = get_obs_key(obs)
      write(msgstring1,*)'obs number ', key, ' has earlier time than previous obs'
      write(msgstring2,*)'observations must be in increasing time order, file ' // trim(filename)
      call error_handler(E_ERR,'validate_obs_seq_time', msgstring2, &
                         source, text2=msgstring1)
   endif

   last_time = this_time
   obs_count = obs_count + 1

   call get_next_obs(seq, obs, next_obs, is_this_last)
   if (.not. is_this_last) obs = next_obs

enddo ObsLoop

! clean up
call destroy_obs(     obs)
call destroy_obs(next_obs)

! technically not a time validation, but easy to check.  obs_count should never
! be larger than size_seq - that's a fatal error.  obs_count < size_seq would
! suggest there are obs in the file that aren't part of the linked list.
! this does not necessarily indicate a fatal error but it's not a common
! situation and might indicate someone should check on the file.
if (obs_count /= size_seq) then
   write(msgstring,*) 'input sequence ', trim(filename)
   call error_handler(E_MSG,'validate_obs_seq_time', msgstring)

   write(msgstring,*) 'total obs in file: ', size_seq, '  obs in linked list: ', obs_count
   if (obs_count > size_seq) then
      ! this is a fatal error
      write(msgstring1,*) 'linked list obs_count > total size_seq, should not happen'
      call error_handler(E_ERR,'validate_obs_seq_time', msgstring, &
                         source, text2=msgstring1)
   else
      ! just warning msg
      write(msgstring1,*) 'only observations in linked list will be processed'
      call error_handler(E_MSG,'validate_obs_seq_time', msgstring, &
                         source, text2=msgstring1)
   endif
endif

end subroutine validate_obs_seq_time

!---------------------------------------------------------------------

subroutine print_metadata(seq, fname)

! print out the metadata strings, trimmed

type(obs_sequence_type), intent(in) :: seq
character(len=*),        intent(in) :: fname

integer :: num_copies , num_qc, i
character(len=metadatalength) :: str

num_copies = get_num_copies(seq)
num_qc     = get_num_qc(    seq)

if ( num_copies < 0 .or. num_qc < 0 ) then
   write(msgstring,*)' illegal copy or obs count in file '//trim(fname)
   call error_handler(E_ERR, 'print_metadata', msgstring, source)
endif

MetaDataLoop : do i=1, num_copies
   str = get_copy_meta_data(seq,i)

   write(msgstring,*)'Data Metadata: ',trim(str)
   call error_handler(E_MSG, '', msgstring)

enddo MetaDataLoop

QCMetaData : do i=1, num_qc
   str = get_qc_meta_data(seq,i)

   write(msgstring,*)'  QC Metadata: ', trim(str)
   call error_handler(E_MSG, '', msgstring)

enddo QCMetaData

end subroutine print_metadata

!---------------------------------------------------------------------

function all_good(obs_in, count, indx, threshold)

! compare location, time, type, QC, plus make sure QC <= threshold

type(obs_type),           intent(in) :: obs_in(:)
integer,                  intent(in) :: count, indx, threshold
logical                              :: all_good

type(obs_def_type)  :: test1_obs_def,  testN_obs_def
integer             :: test1_obs_type, testN_obs_type
type(time_type)     :: test1_obs_time, testN_obs_time
type(location_type) :: test1_obs_loc,  testN_obs_loc
integer             :: test1_qc,       testN_qc
real(r8)            :: temp(1)
integer             :: i, key

! assume failure so we can return as soon as we know they don't match.
all_good = .false.

! get the info we need from sequence 1
call get_obs_def(obs_in(1), test1_obs_def)
test1_obs_loc  = get_obs_def_location(test1_obs_def)
test1_obs_time = get_obs_def_time(test1_obs_def)
test1_obs_type = get_obs_def_type_of_obs(test1_obs_def)
call get_qc(obs_in(1), temp, indx)
test1_qc = nint(temp(1))


! compare the sets against the entries in sequence 1
do i = 2, count
   call get_obs_def(obs_in(i), testN_obs_def)
   testN_obs_loc  = get_obs_def_location(testN_obs_def)
   testN_obs_time = get_obs_def_time(testN_obs_def)
   testN_obs_type = get_obs_def_type_of_obs(testN_obs_def)
   call get_qc(obs_in(i), temp, indx)
   testN_qc = nint(temp(1))

   ! we expect the types, locations, and times to match.  test them
   ! first so we can report those errors if they don't, since that is
   ! not expected.

   if ((test1_obs_type /= testN_obs_type) .or.  &
       (test1_obs_time /= testN_obs_time) .or. &
       (test1_obs_loc  /= testN_obs_loc)) then
 
      ! count of rejected obs.
      num_rejected_other = num_rejected_other + 1

      key = get_obs_key(obs_in(1))
      write(msgstring,*) 'obs number ', key, ' does not match one or more of type/time/location'
      call error_handler(E_MSG, '', msgstring)

      ! more than one reason can be true for a single rejected ob
      if (test1_obs_type /= testN_obs_type) &
         num_mismatch_type = num_mismatch_type + 1
      if (test1_obs_time /= testN_obs_time) &
         num_mismatch_time = num_mismatch_time + 1
      if (test1_obs_loc  /= testN_obs_loc ) &
         num_mismatch_loc  = num_mismatch_loc  + 1

      return
   endif
   
   ! check first to see if either of the DART QCs indicate
   ! that this obs was rejected for some reason.  this code
   ! is testing the first QC multiple times if N > 2 but
   ! it makes the logic simpler so just repeat the test.
   if (test1_qc > threshold .or. testN_qc > threshold) then
      num_rejected_badqc = num_rejected_badqc + 1
      return
   endif

   ! this is the test we expect to be the reason most obs fail:
   ! if they've been assimilated in one experiment but not another.
   if (.not. match_qc(test1_qc, testN_qc)) then
      num_rejected_diffqc = num_rejected_diffqc + 1
      return
   endif

enddo

! all match - good return.
all_good = .true.

end function all_good

!---------------------------------------------------------------------

function match_qc(qc1, qc2)

! normal behavior is to return true if the values are equal,
! and false if not.  but if the "it's ok for evaluated and
! assimilated obs to match", then allow DART QC value 0 to 
! match 1, and 2 to match 3, but no other combinations are valid.

integer, intent(in) :: qc1
integer, intent(in) :: qc2
logical             :: match_qc

! if the combination of values is ok, return early
match_qc = .true.

if (qc1 == qc2) return

if (eval_and_assim_can_match) then
   if (qc1 == 0 .and. qc2 == 1) return
   if (qc1 == 1 .and. qc2 == 0) return

   if (qc1 == 2 .and. qc2 == 3) return
   if (qc1 == 3 .and. qc2 == 2) return
endif

! not one of the allowed matches.  these qcs really
! are not compatible and the match failed.

match_qc = .false.

end function match_qc


!---------------------------------------------------------------------
end program obs_common_subset

