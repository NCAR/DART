c This code is not protected by the DART copyright agreement.
c DART $Id$

      program read_qscat_info

c====================================================================
c     Filename: read_qscat_info.f
c
c     Usage: 
c
c       To run this program, use the following command:
c     
c       your_computer% read_qscat_info <filename>
c
c       where "<filename>" is the name of the QuikSCAT input file
c
c     Description:
c
c       This file contains 2 subroutines in order to read the 
c       Hierarchical Data Format (HDF) attributes in the 
c       QuikSCAT data.  The subroutines are as follows.
c
c       1. GET_TYPE:  a subroutine which returns a string 
c                     containing the data type in ascii given 
c                     the definition value (see HDF documentation 
c                     for conversions). 
c
c       2. READ_ATTRIBUTE:  a subroutine to read the name and 
c                           value(s) of a global attribute referenced 
c                           by its index.
c
c     NOTES:
c     1. Please refer all questions concerning this program and
c        QuikSCAT data obtained from the JPL PO.DAAC to
c        qscat@podaac.jpl.nasa.gov.
c
c     2. The HDF library must be installed before this program will 
c        work properly.  The HDF library and further information 
c        about HDF may be obtained from the National Center for 
c        Supercomputing Applications (NCSA) at http://hdf.ncsa.uiuc.edu.
c
c  7/1/1999 R.S. Dunbar, K.L. Perry
c  Copyright 1999, California Institute of Technology
c
c  DART $Id$
c
c====================================================================

c     Set Parameters
      integer DFACC_RDONLY
      parameter (DFACC_RDONLY = 1)

c     Define Variables
      character infile*128,name*256,string_data_type*10
      character attr_name*256,char_buffer*256
      character*256 values(500)

      integer sd_id,sds_id,retn
      integer n_datasets,n_file_attrs
      integer rank,data_type,dim_sizes(3),nattrs,num_type
      integer attr_index,adata_type,count,n_values
      integer i,j
      integer sfstart,sffinfo,sfselect,sdginfo,sfgcal
      integer sfgainfo,sfrattr,sfendacc,sfend

      double precision cal,cal_err,off,off_err

c     Read the input filename.
      call GETARG(1,infile)

c     Open the HDF input file and initiate the SD interface
      sd_id=sfstart(infile,DFACC_RDONLY)

c     Determine the contents of the input file
      retn=sffinfo(sd_id,n_datasets,n_file_attrs)

c     Print information
      write(*,*) 'HDF info for file: ',infile
      write(*,5) n_datasets
 5    format('Number of Datasets= ',i3)
      write(*,6) n_file_attrs
 6    format('Number of Global Attributes= ',i3)
      write(*,*) ' '
      write(*,*) 'Dataset/Name          ','Rank   ','Dimensions      ',
     &     'Scale        ','Offset     ','Data Type'
      write(*,*) '------------------------------------------------',
     &     '------------------------------'

c     Local Attributes
      if (n_datasets.gt.0) then

c     Access and print the names of every data set in the file
         do i=0,n_datasets-1

c     Initialize dim_sizes
            do j=1,3
               dim_sizes(j)=0
            enddo
            
            sds_id=sfselect(sd_id,i)
            retn=sfginfo(sds_id,name,rank,dim_sizes,data_type,
     &           nattrs)

            call get_type(data_type,string_data_type)

c     Get scaling and offset factors
      retn = sfgcal(sds_id,cal,cal_err,off,off_err,num_type)
            
c     Print local attributes to the screen
            write(*,10) i,name,rank,dim_sizes(1),dim_sizes(2),
     &           dim_sizes(3),cal,off,string_data_type
 10         format(i2,2x,a19,i3,2x,3(i5,x),2x,2(f11.9,2x),a10)

            retn=sfendacc(sds_id)
         enddo
      endif

      write(*,*) ' '
      write(*,*) 'Index/Attribute Name                   ',
     &     'Values'
      write(*,*) '------------------------------------------------',
     &     '------------------------------'

c     Global Attributes
      if (n_file_attrs.gt.0) then
         do i=0,n_file_attrs-1
            call read_attribute(sd_id,i,attr_name,num_type,
     &           n_values,values)
            write(*,110) i,attr_name,values(1)
            if (n_values.gt.1) then
               do j=2,n_values
                  write(*,111) values(j)
               enddo
            endif
 110        format(1x,i3,2x,a30,4x,a60)
 111        format(40x,a60)
         enddo
      endif

      retn=sfend(sd_id)
      end

c====================================================================
c    GET_TYPE:  a subroutine which returns a string containing the 
c               data type in ascii given the definition value 
c               (see HDF documentation for conversions).
c    
c    7/22/1999 K.L. Perry
c====================================================================
      subroutine get_type(data_type,string_data_type)

      integer data_type
      character string_data_type*10

      if (data_type.eq.3) then
         string_data_type='uchar'
      else if (data_type.eq.4) then
         string_data_type='char'
      else if (data_type.eq.5) then
         string_data_type='float32'
      else if (data_type.eq.20) then
         string_data_type='int8'
      else if (data_type.eq.21) then
         string_data_type='uint8'
      else if (data_type.eq.22) then
         string_data_type='int16'
      else if (data_type.eq.23) then
         string_data_type='uint16'
      else if (data_type.eq.24) then
         string_data_type='int32'
      else if (data_type.eq.25) then
         string_data_type='uint32'
      else if (data_type.eq.26) then
         string_data_type='int64'
      else if (data_type.eq.27) then
         string_data_type='uint64'
      else
         string_data_type='UNKNOWN'
         write(*,*) 'Abort Program:  Data Type Unknown'
         stop
      endif

      end

c====================================================================
c    READ_ATTRIBUTE:  a subroutine to read the name and value(s) 
c                     of a global attribute referenced by its 
c                     index.
c    
c    5/14/1998 R.S. Dunbar
c====================================================================

      subroutine read_attribute(sd_id,attr_index,attr_name,
     &     adata_type,n_values,fvalues)
      
      integer MAX_NC_NAME
      parameter (MAX_NC_NAME=256)

      integer sd_id,retn,adata_type,n_values
      integer attr_index,count,retn,n,oldn
      integer sfgainfo,sfrattr
      
      character*(*) fvalues(*)
      character attr_name*(MAX_NC_NAME),attr_data*512
      character*(MAX_NC_NAME) values(25)
      character cr
 
c     Get information about the  file attribute
      retn = sfgainfo(sd_id,attr_index,attr_name,adata_type,count)

c     Read the attribute data
      retn = sfrattr(sd_id,attr_index,attr_data)

      cr = char(10)
      ival = 0
      oldn = 1
 5    continue

c     QuikSCAT attributes have atleast three lines: 
c     metadata type, array size and metadata contents
c     Use "blank spaces" to identify the end of a line

      n = index(attr_data(oldn:(count-1)),cr)

c     Read all of the metadata lines
      if (n .eq. 0) then
         ival=ival+1
         values(ival) = attr_data(oldn:(count-1))
         goto 99
      else
         ival=ival+1
         values(ival) = attr_data(oldn:(oldn+n-2))
      endif
      oldn=n+oldn
      goto 5

 99   continue

      n_values = ival - 2
      do i=1,n_values
         fvalues(i) = values(i+2)
      enddo

      return
      end

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
