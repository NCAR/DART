c This code is not protected by the DART copyright agreement.
c DART $Id$

      program read_qscat2b

c====================================================================
c     Filename: read_qscat2b.f
c
c     Usage: 
c
c       To run this program, use the following command:
c     
c       your_computer% read_qscat2b <filename>
c
c       where "<filename>" is the name of the QuikSCAT Level 2B
c       input file
c
c     Description:
c
c       This file contains 3 subroutines in order to read the 
c       QuikSCAT Level 2B data in Hierarchical Data Format (HDF).  
c       The subroutines are as follows.
c
c       1. read_attrib_byname():  a subroutine to read the name 
c                                 and value(s) of a global attribute
c                                 referenced by its name.
c       
c       2. read_timetags():  a subroutine to read the timetag info
c                            contained in the HDF VDATA
c
c       3. extract_sds():  a subroutine to read the contents of an
c                          SDS from an HDF file
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
c     3. The L2B data are read in their entirety.  Examples of reading
c        the QuikSCAT data by slabs can be found in read_qscat1b.f
c        read_qscat2a.f.
c
c  7/1/1999 R.S. Dunbar, K.L. Perry
c
c  Modifications:
c
c  11/18/1999 Added capability to read new SDSs,
c             wind_speed_selection and wind_dir_selection.
c             K.L. Perry
c
c   1/12/2000 Added capability to read new SDSs, 
c             mp_rain_probability and nof_rain_index.  K.L. Perry
c
c   3/10/2000 Corrected handling of 8-bit unsigned integers
c             in extract_sds subroutine.  K.L. Perry
c
c   5/15/2000 Corrected handling of 32-bit unsigned integers
c             in extract_sds subroutine.  This change should
c             not affect the Level 2A or 2B software.  K.L. Perry
c
c   7/13/2006 sdrad_rain_rate and adapted code for 12.5 km product reading
C             Ted Lungu
c  Copyright 1999-2006, California Institute of Technology
c
c  DART $Id$
c
c====================================================================

c     Set Parameters

      integer MAX_ROWS,MAX_CELLS,MAX_SIG

      parameter (MAX_ROWS = 1624)
      parameter (MAX_WVC = 76)
c      parameter (MAX_ROWS = 3248)
c      parameter (MAX_WVC = 152)
      parameter (MAX_AMBIG = 4)

      integer DFACC_RDONLY
      parameter (DFACC_RDONLY = 1)

c     Define Variables

      character l2b_file*128,product*8
      character*21 TimeTags(MAX_ROWS)

      integer sd_id,retn,sfstart,sfend
      integer irec1,irec2,itmp,irow,iwvc,iamb

      real wvc_row(MAX_ROWS)
      real wvc_lat(MAX_WVC,MAX_ROWS),wvc_lon(MAX_WVC,MAX_ROWS)
      real wvc_index(MAX_WVC,MAX_ROWS)
      real num_in_fore(MAX_WVC,MAX_ROWS),num_in_aft(MAX_WVC,MAX_ROWS)
      real num_out_fore(MAX_WVC,MAX_ROWS),num_out_aft(MAX_WVC,MAX_ROWS)
      real wvc_quality_flag(MAX_WVC,MAX_ROWS)
      real atten_corr(MAX_WVC,MAX_ROWS),model_speed(MAX_WVC,MAX_ROWS)
      real model_dir(MAX_WVC,MAX_ROWS),num_ambigs(MAX_WVC,MAX_ROWS)
      real wind_speed(MAX_AMBIG,MAX_WVC,MAX_ROWS)
      real wind_dir(MAX_AMBIG,MAX_WVC,MAX_ROWS)
      real wind_speed_err(MAX_AMBIG,MAX_WVC,MAX_ROWS)
      real wind_dir_err(MAX_AMBIG,MAX_WVC,MAX_ROWS)
      real max_likelihood_est(MAX_AMBIG,MAX_WVC,MAX_ROWS)
      real wvc_selection(MAX_WVC,MAX_ROWS)
      real wind_speed_selection(MAX_WVC,MAX_ROWS)
      real wind_dir_selection(MAX_WVC,MAX_ROWS)
      real mp_rain_probability(MAX_WVC,MAX_ROWS)
      real nof_rain_index(MAX_WVC,MAX_ROWS)
      real srad_rain_rate(MAX_WVC,MAX_ROWS)

c     Read the input filename.
      call GETARG(1,l2b_file)
      if (l2b_file .eq. ' ') then
         print *,'Usage: read_qscat2b <Level 2B file>'
         stop
      endif
      write(*,*)
      write(*,*) 'FILENAME: ',l2b_file

c     Open the HDF input file and initiate the SD interface
      sd_id=sfstart(l2b_file,DFACC_RDONLY)

c     Make sure that the file is a QuikSCAT Level 2B file
      call read_attrib_byname(sd_id,'ShortName',ntype,nval,product)
      if (product.ne.'QSCATL2B') then
         print *,'The input file is not a QuikSCAT Level 2B file'
         print *,'*** Aborting program ***'
         stop
      endif

c     Read the timetag info contained in the HDF VDATA
      call read_timetags(l2b_file, TimeTags)

c     Read each SDS in its entirety.  For an example of reading
c     the QuikSCAT SDS data in slabs, please refer to read_qscat2a.f.

      irow=1
      call extract_sds(sd_id,'wvc_row',irow,MAX_ROWS,wvc_row)
      call extract_sds(sd_id,'wvc_lat',irow,MAX_ROWS,wvc_lat)
      call extract_sds(sd_id,'wvc_lon',irow,MAX_ROWS,wvc_lon)
      call extract_sds(sd_id,'wvc_index',irow,MAX_ROWS,wvc_index)
      call extract_sds(sd_id,'num_in_fore',irow,MAX_ROWS,num_in_fore)
      call extract_sds(sd_id,'num_in_aft',irow,MAX_ROWS,num_in_aft)
      call extract_sds(sd_id,'num_out_fore',irow,MAX_ROWS,num_out_fore)
      call extract_sds(sd_id,'num_out_aft',irow,MAX_ROWS,num_out_aft)
      call extract_sds(sd_id,'wvc_quality_flag',irow,MAX_ROWS,
     &     wvc_quality_flag)
      call extract_sds(sd_id,'atten_corr',irow,MAX_ROWS,atten_corr)
      call extract_sds(sd_id,'model_speed',irow,MAX_ROWS,model_speed)
      call extract_sds(sd_id,'model_dir',irow,MAX_ROWS,model_dir)
      call extract_sds(sd_id,'num_ambigs',irow,MAX_ROWS,num_ambigs)
      call extract_sds(sd_id,'wind_speed',irow,MAX_ROWS,wind_speed)
      call extract_sds(sd_id,'wind_dir',irow,MAX_ROWS,wind_dir)
      call extract_sds(sd_id,'wind_speed_err',irow,MAX_ROWS,
     &     wind_speed_err)
      call extract_sds(sd_id,'wind_dir_err',irow,MAX_ROWS,wind_dir_err)
      call extract_sds(sd_id,'max_likelihood_est',irow,MAX_ROWS,
     &     max_likelihood_est)
      call extract_sds(sd_id,'wvc_selection',irow,MAX_ROWS,
     &     wvc_selection)
      call extract_sds(sd_id,'wind_speed_selection',irow,MAX_ROWS,
     &     wind_speed_selection)
      call extract_sds(sd_id,'wind_dir_selection',irow,MAX_ROWS,
     &     wind_dir_selection)
      call extract_sds(sd_id,'mp_rain_probability',irow,MAX_ROWS,
     &     mp_rain_probability)
      call extract_sds(sd_id,'nof_rain_index',irow,MAX_ROWS,
     &     nof_rain_index)
      call extract_sds(sd_id,'srad_rain_rate',irow,MAX_ROWS,
     &     srad_rain_rate)

c     Select the wind vector cell rows to be read
      write(*,12)MAX_ROWS
 12   format('Enter the first and last record numbers [1-', i4,']: '$)
      read(*,*) irec1,irec2
      write(*,*) '  '

c     Make sure that "first comes before last"
      if (irec1.gt.irec2) then
         itmp=irec1
         irec1=irec2
         irec2=itmp
      endif

c     Check to see if selected rows are within limits
      if ((irec1.lt.1).or.(irec2.gt.MAX_ROWS)) then
         write(*,*) 'Number of rows must be between 1 and ', MAX_ROWS
         print *,'*** Aborting program ***'
         stop
      endif

c     Print results to screen
      do irow=irec1,irec2

         write(*,*) ' '
         write(*,*) 'TIME: ', TimeTags(irow)
         write(*,100) wvc_row(irow)
 100     format('WVC ROW: ',f5.0)

        write(*,105)
 105    format('WVC#','  WVC_Qual',
     &       '  WVC Latitude/Longitude','  Selected Wind Vector',
     &       ' NWP Wind Vector','  Num/Sel ambig','  DRE Wind Vector',
     &       ' MUDH Prob',' NOF Index',' SRR')

        do iwvc = 1,MAX_WVC
           if (num_ambigs(iwvc,irow).gt.0) then
              iamb=wvc_selection(iwvc,irow)
              write(*,110) wvc_index(iwvc,irow),
     &             int(wvc_quality_flag(iwvc,irow)),wvc_lat(iwvc,irow),
     &             wvc_lon(iwvc,irow),wind_speed(iamb,iwvc,irow),
     &             wind_dir(iamb,iwvc,irow),model_speed(iwvc,irow),
     &             model_dir(iwvc,irow),num_ambigs(iwvc,irow),
     &             wvc_selection(iwvc,irow),
     &             wind_speed_selection(iwvc,irow),
     &             wind_dir_selection(iwvc,irow),
     &             mp_rain_probability(iwvc,irow),
     &             nof_rain_index(iwvc,irow),
     &             srad_rain_rate(iwvc,irow)

 110          format(f3.0,4x,"0X",z4.4,8x,f6.2,3x,f6.2,6x,f6.2,3x,f6.2,
     &             4x,f6.2,3x,f6.2,7x,f2.0,3x,f2.0,3x,f6.2,3x,f6.2,
     &             3x,f6.2,3x,f6.0,2x, f6.2)
           endif
        enddo
      enddo

      retn=sfend(sd_id)
      end

c====================================================================
c    READ_ATTRIB_BYNAME:  a subroutine to read the name and
c                         value(s) of a global attribute
c                         referenced by its name.
c    
c    5/14/1998 R.S. Dunbar
c====================================================================

      subroutine read_attrib_byname(sd_id,in_attr_name,
     $     num_type,n_values,fvalues)
      
      integer MAX_NC_NAME
      parameter (MAX_NC_NAME=256)

      integer sd_id,num_type,n_values
      integer attr_index,count,retn,n,oldn
      integer sffattr,sfgainfo,sfrattr
      character*(*) in_attr_name
      character*(*) fvalues(*)
      character attr_name*(MAX_NC_NAME),attr_data*512
      character*(MAX_NC_NAME) values(20)
      character cr
 
c     Find the attribute assigned to in_attr_name
      attr_index = sffattr(sd_id,in_attr_name)

c     Get information about the  file attribute
      retn = sfgainfo(sd_id,attr_index,attr_name,num_type,count)

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

c====================================================================
c    READ_TIMETAGS:  a subroutine to read the timetag info
c                    contained in the HDF VDATA
c    
c    5/1998 R.S. Dunbar
c
c    Revisions:
c    7/1999 Code adapted to read timetags in their entirety.
c           Commenter were also added.  K.L. Perry
c====================================================================
      subroutine read_timetags(filename,timetags)

      character*80 filename
      character*21 timetags(*)
      character*60 fields
      character vdata_name*30
      integer file_id,vdata_ref,vdata_id
      integer n_records,interlace,vdata_size
      integer hopen,vsfgid,vsfatch,vsfinq,vsfread,vfsdtch,hclose

      integer DFACC_RDONLY,FULL_INTERLACE
      parameter(DFACC_RDONLY=1)
      parameter(FULL_INTERLACE=0)

c     Open the HDF file
      file_id = hopen(filename,DFACC_RDONLY,0)

c     Initialize the VS interface
      call vfstart(file_id)

c     Get the reference number for the first vdata in the file
      vdata_ref = -1
      vdata_ref = vsfgid(file_id,vdata_ref)

c     Attach to the vdata for reading if it is found, otherwise 
c     exit the program.
      if (vdata_ref.eq.0) then
         print *,'No Timetags were found in the HDF VDATA'
         print *,'*** Aborting program ***'
         stop
      endif

      vdata_id = vsfatch(file_id,vdata_ref,'r')

c     Get n_records
      retn=vsfinq(vdata_id,n_records,interlace,fields,
     &     vdata_size,vdata_name)

c     Read the timetags
      retn = vsfread(vdata_id,timetags,n_records,FULL_INTERLACE)

c     Terminate access to the vdata and to the VS interface, 
c     then close the HDF file.

      retn =  vsfdtch(vdata_id)
      call vfend(file_id)
      retn = hclose(file_id)

      return
      end

c====================================================================
c    EXTRACT_SDS:  a subroutine to read the contents of an
c                  SDS from an HDF file
c    
c    5/12/1998 R.S. Dunbar
c
c    Revisions:
c    7/1999   Code adapted to read input in bytes as well as ints 
c             and floats.  Comments were also added.  K.L. Perry
c
c    3/2000   Corrected code for 8-bit unsigned integers.  "buffer"
c             was used instead of "buffer2".  K.L. Perry
c
c    5/2000   Changed MAX_BUF_SIZE from 1000000 to 10000000.
c             Corrected code for 32-bit unsigned integers.  Created
c             "buffer3" array of int*4 to correctly read in uint32.
c             K.L. Perry
c
c====================================================================
      subroutine extract_sds(sd_id,in_var,irec,slab_size,out_var)

      integer MAX_BUF_SIZE
      parameter (MAX_BUF_SIZE=10000000)
      integer sd_id,sds_index,sds_id,retn
      integer rank,dim_sizes(3),data_type,nattrs,num_type
      integer edge(3),stride(3),start(3),irec,slab_size
      double precision cal,cal_err,off,off_err
      integer iprod,i,itmp

      character*(*) in_var
      character name*256
      integer sfn2index,sfselect,sfginfo,sfrdata,sfgcal,sfendacc

      integer*2 buffer(MAX_BUF_SIZE)
      byte buffer2(MAX_BUF_SIZE)
      integer*4 buffer3(MAX_BUF_SIZE)
      real out_var(MAX_BUF_SIZE)

c     Search for the index of "in_var"
      sds_index = sfn2index(sd_id, in_var)

c     Select data set corresponding to the returned index
      sds_id = sfselect(sd_id,sds_index)
      retn = sfginfo(sds_id,name,rank,dim_sizes,data_type,nattrs)

      do i=1,rank
         edge(i)=dim_sizes(i)
         start(i)=0
         stride(i)=1
      enddo
      edge(rank)=slab_size
      start(rank)=irec-1

      iprod=1
      do i=1,rank
         iprod=iprod*edge(i)
      enddo

c     Get the calibration and offset values of input
      retn = sfgcal(sds_id,cal,cal_err,off,off_err,num_type)

c     Read Arrays which are not float32 or int8 or uint8 or uint32
      if ((data_type.ne.5).and.(data_type.ne.20).and.
     &     (data_type.ne.21).and.(data_type.ne.25)) then

c     Read the data set into the "buffer" array
         retn=sfrdata(sds_id,start,stride,edge,buffer)

c     Calibrate the output
         do i=1,iprod
c     Correct for 16-bit unsigned integers
            if ((data_type.eq.23).and.(buffer(i).lt.0)) then
               out_var(i)=buffer(i)+65536.0

c     No correction needed for signed or positive unsigned integers
            else
               out_var(i)=buffer(i)
            endif

            out_var(i)=out_var(i)*cal
         enddo

c     Read uint32 arrays.
      else if (data_type.eq.25) then

c     Read the data set into the "buffer3" uint32 array
         retn=sfrdata(sds_id,start,stride,edge,buffer3)

c     Calibrate the output
         do i=1,iprod
c     Correct for 32-bit unsigned integers
            if ((data_type.eq.25).and.(buffer3(i).lt.0)) then
               out_var(i)=buffer3(i)+4294967296.0
            else
               out_var(i)=buffer3(i)
            endif
            out_var(i)=out_var(i)*cal
         enddo

c     Read int8 and uint8 arrays. 
      else if ((data_type.eq.20).or.(data_type.eq.21)) then

c     Read the data set into the "buffer2" byte-array
         retn=sfrdata(sds_id,start,stride,edge,buffer2)

c     Calibrate the output
         do i=1,iprod

c     Correct for 8-bit unsigned integers
            itmp=buffer2(i)
            if ((data_type.eq.21).and.(buffer2(i).lt.0)) then
               itmp=itmp+256
               out_var(i)=itmp
               
c     No correction needed for signed or positive unsigned integers
            else
               out_var(i)=itmp
            endif

            out_var(i)=out_var(i)*cal
         enddo

      else
c     Read float32 arrays directly into the "out_var" array
         retn=sfrdata(sds_id,start,stride,edge,out_var)

c     Calibrate the output
         do i=1,iprod
            out_var(i)=out_var(i)*cal
         enddo
      endif

c     Terminate access to the array data set.
      retn = sfendacc(sds_id)
      end

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
