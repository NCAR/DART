! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program make_netcdf_table

! sampling_error_correction_table.nc
!
! netcdf sampling_error_correction_table {
! dimensions:
!         bins = 200 ;
!         ens_sizes = UNLIMITED ; // (40 currently)
! variables:
!         int count(ens_sizes, bins) ;
!                 count:description = "number of samples in each bin" ;
!         double true_corr_mean(ens_sizes, bins) ;
!         double alpha(ens_sizes, bins) ;
!                 alpha:description = "sampling error correction factors" ;
!         int ens_sizes(ens_sizes) ;
!                 ens_sizes:description = "ensemble size used for calculation" ;
! 
! // global attributes:
!                 :num_samples = 100000000 ;
!                 :title = "Sampling Error Corrections for fixed ensemble sizes." ;
!                 :reference = "Anderson, J., 2012: Localization and Sampling Error Correction in Ensemble Kalman Filter Data Assimilation. Mon. Wea. Rev., 140, 2359-2371, doi: 10.1175/MWR-D-11-00013.1." ;
!                 :version = "$Id$" ;
! 
!
! NF90_OPEN             ! open existing netCDF dataset
!    NF90_redef         ! put into define mode 
!    NF90_def_dim       ! define additional dimensions (if any)
!    NF90_def_var       ! define variables: from name, type, and dims
!    NF90_put_att       ! assign attribute values
! NF90_ENDDEF           ! end definitions: leave define mode
!    NF90_put_var       ! provide values for variable
! NF90_CLOSE            ! close: save updated netCDF dataset
!
!

use netcdf
use typesizes
use types_mod
use utilities_mod

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

integer, parameter :: bins = 200
integer, parameter :: nens = 40

integer :: ncid, binDimID, ensDimID
integer :: countVarID, corrVarID, alphaVarID, sizesVarID

integer :: iunit, io, ifile, iline, mysize, bincounter

integer :: ens_sizes(nens) = (/ 5,   6,   7,   8,   9,  10,  12,  14,  15,  16, &
                               18,  20,  22,  24,  28,  30,  32,  36,  40,  44, &
                               48,  49,  50,  52,  56,  60,  64,  70,  72,  80, &
                               84,  88,  90,  96, 100, 120, 140, 160, 180, 200 /)

integer  ::      counts(bins,nens)
real(r8) :: correlation(bins,nens)
real(r8) ::       alpha(bins,nens)

character(len=512) :: fname

io = nf90_create(path='sampling_error_correction_table.nc', cmode=NF90_SHARE, ncid=ncid)
io = nf90_put_att(ncid, NF90_GLOBAL, "num_samples", 100000000)
io = nf90_put_att(ncid, NF90_GLOBAL, "title", 'Sampling Error Corrections for fixed ensemble sizes.')
io = nf90_put_att(ncid, NF90_GLOBAL, "version", 'original values as from final_full.nn')
io = nf90_put_att(ncid, NF90_GLOBAL, "reference", "Anderson, J., 2012: Localization and Sampling Error Correction in Ensemble Kalman Filter Data Assimilation. Mon. Wea. Rev., 140, 2359-2371, doi: 10.1175/MWR-D-11-00013.1.")

io = nf90_def_dim(ncid, name="bins",      len=bins,           dimid = binDimID)
io = nf90_def_dim(ncid, name="ens_sizes", len=NF90_UNLIMITED, dimid = ensDimID)

io = nf90_def_var(ncid, name='count', xtype=NF90_INT, dimids=(/binDimID,ensDimID/), varid=countVarID )
io = nf90_put_att(ncid, countVarID, "description", "number of samples in each bin")

io = nf90_def_var(ncid, name='true_corr_mean', xtype=NF90_DOUBLE, dimids=(/binDimID,ensDimID/), varid=corrVarID )

io = nf90_def_var(ncid, name='alpha', xtype=NF90_DOUBLE, dimids=(/binDimID,ensDimID/), varid=alphaVarID )
io = nf90_put_att(ncid, alphaVarID, "description", "sampling error correction factors")

io = nf90_def_var(ncid, name='ens_sizes', xtype=NF90_INT, dimids=(/ensDimID/), varid=sizesVarID )
io = nf90_put_att(ncid, sizesVarID, "description", "ensemble size used for calculation")

io = nf90_enddef(ncid)

do ifile = 1,nens   ! 40 files to read

   mysize = ens_sizes(ifile)

   if (mysize < 10) then
      write(fname,'(''/Users/thoar/svn/DART/trunk/system_simulation/final_full_precomputed_tables/final_full.'',i1)')mysize
   elseif (mysize < 100) then
      write(fname,'(''/Users/thoar/svn/DART/trunk/system_simulation/final_full_precomputed_tables/final_full.'',i2)')mysize
   else
      write(fname,'(''/Users/thoar/svn/DART/trunk/system_simulation/final_full_precomputed_tables/final_full.'',i3)')mysize
   endif

   write(*,*)'Reading "'//trim(fname)//'"'

   iunit = open_file(fname,'formatted','read')

   ! 4 columns == bin    count     correlation   alpha
   READLOOP : do iline = 1,bins  ! 200 lines per file

      read(iunit, *, iostat=io) bincounter, counts(iline,ifile), correlation(iline,ifile), alpha(iline,ifile)
      if (io /= 0) then
         print *, 'Reading file ',ifile, 'ens size is ',mysize
         print *, 'got bad read code from input file, rcio = ', io, 'line ',iline
         stop
      endif

   enddo READLOOP

   call close_file(iunit)

enddo

io = nf90_put_var(ncid, sizesVarID,   ens_sizes)
io = nf90_put_var(ncid, alphaVarID,       alpha)
io = nf90_put_var(ncid,  corrVarID, correlation)
io = nf90_put_var(ncid, countVarID,      counts)

io = nf90_sync(ncid)
io = nf90_close(ncid)

end program make_netcdf_table

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
