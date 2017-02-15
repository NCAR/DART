! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id: make_single_netcdf_file.f90 10790 2016-12-15 20:55:16Z nancy@ucar.edu $

! To compile : assuming that you have your NETCDF environment variable set.
! gfortran -o make_single_file make_single_file.f90 -I $NETCDF/include \ 
!          -L $NETCDF/lib -lnetcdff -lnetcdf

program make_single_netcdf_file

use netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL: https://svn-dares-dart.cgd.ucar.edu/DART/branches/rma_single_file/models/lorenz_96make_single_netcdf_file $"
character(len=32 ), parameter :: revision = "$Revision: 10790 $"
character(len=128), parameter :: revdate  = "$Date: 2016-12-15 13:55:16 -0700 (Thu, 15 Dec 2016) $"


  ! This is the name of the data file we will read.
  character (len = *), parameter :: FILE_NAME = "work/filter_input.nc"

  ! When we create netCDF files, variables and dimensions, we get back
  ! an ID for each one.
  integer :: ncid, varid
  integer :: x_dimid, y_dimid
  
  integer :: ios, i, imem
  real(SELECTED_REAL_KIND(12)) , dimension(80,40) :: data_out
  real(SELECTED_REAL_KIND(12)) :: val
  character(len = 64) :: line

  ! Member Data
  open (unit = 27, file="work/filter_ics")
  
  i = 0
  imem = 0
  do
     read(27,'(A)',iostat = ios) line

     if (ios < 0) exit

     if (imem == 41) imem = 0

     if (imem == 0) then
        i = i + 1
        write(*,*) "new_member" , i
     else
        read(line,*) val 
        write(*,*) 'val = ', val
        data_out(i, imem)  = val
     endif

     imem = imem + 1

  enddo
  close(27) 

  call check( nf90_open(FILE_NAME, NF90_WRITE, ncid) )
  
  ! Get the varid of the data variable, based on its name.
  call check( nf90_inq_varid(ncid, "state", varid) )

  do i = 1,80
     ! Write the pretend data to the file. Although netCDF supports
     ! reading and writing subsets of data, in this case we write all the
     ! data in one operation.
     call check( nf90_put_var(ncid, varid, data_out(i,:), start=(/1,i,1/)) )
  enddo
  ! Close the file. This frees up any internal netCDF resources
  ! associated with the file, and flushes any buffers.
  call check( nf90_close(ncid) )

  print *, "*** SUCCESS writing example file ", FILE_NAME

contains
  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check 

end program
