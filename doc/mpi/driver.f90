! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program driver

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

use commtest_mod, only:  setup, transpose_state_compl, select_statevars, &
                         compute_priors, compute_qc, broadcast_obs, &
                         broadcast_increments, advance_model, &
                         transpose_ensemb_compl, takedown

!-----------------------------------------------------------------------------
! ---- namelist (saved in file commtest.nml)
! (will need more options eventually; high/low latency algorithm; which
! flavor of transpose, etc.)  this is enough to get started.
!
integer :: state_vector_size, ensemble_size, &
           observation_count, repeat_count, &
           state_complete_algorithm, ensemble_complete_algorithm, &
           select_algorithm, broadcast_algorithm


namelist /commtest_nml/ &
           state_vector_size, ensemble_size, &
           observation_count, repeat_count, &
           state_complete_algorithm, ensemble_complete_algorithm, &
           select_algorithm, broadcast_algorithm

!-----------------------------------------------------------------------------


   integer :: N, M, K, R
   integer :: i, j
   integer, pointer :: L(:)
   double precision, pointer :: P(:), QC(:)
 

   print *, "driver program start"

   nullify(L, P, QC)

   call read_namelist()
   print *, "namelist read in"

   ! shorthand which matches the paper and makes it easier to call the
   ! subroutines below.  remove them if it gets more confusing than helpful.
   N = ensemble_size
   M = state_vector_size
   R = repeat_count
   K = observation_count
 

   call setup(M, N, K)
   print *, "returned from setup"
  
   call broadcast_obs(K)

   print *, "start of main loop, R=", R
   MainLoop: do i=1, R
   
      ! only for timing
      call advance_model(M, N)

      ! timing plus allocate the P(:) array
      call compute_priors(M, K, N, P)

      ! timing plus allocate the QC(:) array
      call compute_qc(K, N, QC)

      ! decide which state vars get collected together across the ensemble
      ! the index list gets allocated here and returned in L(:)
      call select_statevars(M, K, L, select_algorithm)

      ! transpose to collect the selected state vars together.
      ! multiple possible algorithm choices here.
      call transpose_ensemb_compl(M, N, L, P, QC, ensemble_complete_algorithm)

      ! this loop might have to drop down into the broadcast routine.
      AssimLoop: do j=1, K

         ! each owner computes and sends in turn the increments
         call broadcast_increments(j, K, L, N, broadcast_algorithm)

      enddo AssimLoop

      ! transpose back.  again, multiple possible algorithm choices.
      call transpose_state_compl(M, N, L, P, QC, state_complete_algorithm)
    
   enddo MainLoop

   call takedown()
   print *, "returned from takedown"

   if (associated(L))  deallocate(L, stat=i)
   if (associated(P))  deallocate(P, stat=i)
   if (associated(QC)) deallocate(QC, stat=i)

   print *, "driver program end"

contains

subroutine read_namelist
   integer :: iunit, ioerr
   character(len=132) :: fname
  
   iunit = 19   ! something unlikely to be in use
   fname = "commtest.nml"
   ! on at least one machine, somehow the namelist is not found
   ! unless it is a full pathname. 
   !fname = "/home/nancy/dart/DART/doc/mpi/commtest.nml"
   open(iunit, status="old", file=fname, action="read", iostat=ioerr)
   if (ioerr /= 0) then
      print *, "namelist open returned error, code = ", ioerr
      print *, "unable to open namelist file <",trim(fname),">"
      print *, "stopping"
      stop
   endif

   read (iunit, nml = commtest_nml, iostat = ioerr)
   if (ioerr /= 0) then
      print *, "namelist read returned error, code = ", ioerr
      print *, "unable to read namelist file <",trim(fname),">"
      print *, "stopping"
      stop
   endif
   close (iunit)

   ! output the namelist to document the settings for this run 
   print *, "namelist read in, values are:"
   write(*, nml=commtest_nml)

end subroutine read_namelist

end program driver

