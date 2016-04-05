!
! code to perturb the wrfchem emission files
!
! ifort -C perturb_chem_emiss.f90 -o perturb_chem_emiss.exe -lgfortran -lnetcdff -lnetcdf
!
          program main
             implicit none
             integer,parameter                        :: nx=100
             integer,parameter                        :: ny=40
             integer,parameter                        :: nz=33
             integer,parameter                        :: nz_chem=11
             integer,parameter                        :: nchem_spc=39
             integer,parameter                        :: nfire_spc=34
             integer,parameter                        :: nbio_spc=1
             integer                                  :: unit,isp,nseed,idate,itime,icnt
             integer,allocatable,dimension(:)         :: seed
             integer,dimension(8)                     :: values
             real                                     :: pi,ran1,ran2,tr_mean,tr_stdev,ens_member
             real                                     :: z_stdn
             real,allocatable,dimension(:,:)          :: chem_data2d
             real,allocatable,dimension(:,:,:)        :: chem_data3d
             real                                     :: z_ran,y_ran
             character(len=5)                         :: czone
             character(len=8)                         :: cdate
             character(len=10)                        :: ctime
             character(len=150)                       :: wrfchemi,wrffirechemi,wrfbiochemi
             character(len=150),dimension(nchem_spc)  :: ch_chem_spc 
             character(len=150),dimension(nfire_spc)  :: ch_fire_spc 
             character(len=150),dimension(nbio_spc)   :: ch_bio_spc 
             logical                                  :: pert_chem,pert_fire,pert_bio
             namelist /perturb_chem_emiss_nml/idate,ens_member,tr_mean,tr_stdev,wrfchemi, &
                       wrffirechemi,wrfbiochemi,pert_chem,pert_fire,pert_bio
!
! Assign constants
             pi=4.*atan(1.)
             z_stdn=2.58     ! 99%
             z_stdn=1.96     ! 95%
!
! Read namelist
             unit=20
             open(unit=unit,file='perturb_chem_emiss_nml.nl',form='formatted', &
             status='old',action='read')
             read(unit,perturb_chem_emiss_nml)
             close(unit)
!             print *, 'date         ',idate
!             print *, 'ens_member   ',ens_member
!             print *, 'tr_mean      ',tr_mean
!             print *, 'tr_stdev     ',tr_stdev
!             print *, 'wrfchemi     ',trim(wrfchemi)
!             print *, 'wrffirechemi ',trim(wrffirechemi)
!             print *, 'wrfbiochemi  ',trim(wrfbiochemi)
!             print *, 'pert_chem    ',pert_chem
!             print *, 'pert_fire    ',pert_fire
!             print *, 'pert_bio     ',pert_bio
!
! Assign emission species names
             ch_chem_spc(1)='E_CO'
             ch_chem_spc(2)='E_CO02'
             ch_chem_spc(3)='E_CO03'
             ch_chem_spc(4)='E_CO04'
             ch_chem_spc(5)='E_CO05'
             ch_chem_spc(6)='E_CO06'
             ch_chem_spc(7)='E_CO07'
             ch_chem_spc(8)='E_CO08'
             ch_chem_spc(9)='E_CO09'
             ch_chem_spc(10)='E_CO10'
             ch_chem_spc(11)='E_CO11'
             ch_chem_spc(12)='E_CO12'
             ch_chem_spc(13)='E_CO13'
             ch_chem_spc(14)='E_CO14'
             ch_chem_spc(15)='E_NO2'
             ch_chem_spc(16)='E_BIGALK'
             ch_chem_spc(17)='E_BIGENE'
             ch_chem_spc(18)='E_C2H4'
             ch_chem_spc(19)='E_C2H5OH'
             ch_chem_spc(20)='E_C2H6'
             ch_chem_spc(21)='E_C3H6'
             ch_chem_spc(22)='E_C3H8'
             ch_chem_spc(23)='E_CH2O'
             ch_chem_spc(24)='E_CH3CHO'
             ch_chem_spc(25)='E_CH3COCH3'
             ch_chem_spc(26)='E_CH3OH'
             ch_chem_spc(27)='E_MEK'
             ch_chem_spc(28)='E_SO2'
             ch_chem_spc(29)='E_TOLUENE'
             ch_chem_spc(30)='E_NH3'
             ch_chem_spc(31)='E_ISOP'
             ch_chem_spc(32)='E_C10H16'
             ch_chem_spc(33)='E_XNO'
             ch_chem_spc(34)='E_XNO2'
             ch_chem_spc(35)='E_sulf'
             ch_chem_spc(36)='E_PM_25'
             ch_chem_spc(37)='E_PM_10'
             ch_chem_spc(38)='E_OC'
             ch_chem_spc(39)='E_BC'
!
             ch_fire_spc(1)='ebu_in_co'
             ch_fire_spc(2)='ebu_in_no'
             ch_fire_spc(3)='ebu_in_so2'
             ch_fire_spc(4)='ebu_in_bigalk'
             ch_fire_spc(5)='ebu_in_bigene'
             ch_fire_spc(6)='ebu_in_c2h4'
             ch_fire_spc(7)='ebu_in_c2h5oh'
             ch_fire_spc(8)='ebu_in_c2h6'
             ch_fire_spc(9)='ebu_in_c3h8'
             ch_fire_spc(10)='ebu_in_c3h6'
             ch_fire_spc(11)='ebu_in_ch2o'
             ch_fire_spc(12)='ebu_in_ch3cho'
             ch_fire_spc(13)='ebu_in_ch3coch3'
             ch_fire_spc(14)='ebu_in_ch3oh'
             ch_fire_spc(15)='ebu_in_mek'
             ch_fire_spc(16)='ebu_in_toluene'
             ch_fire_spc(17)='ebu_in_nh3'
             ch_fire_spc(18)='ebu_in_no2'
             ch_fire_spc(19)='ebu_in_open'
             ch_fire_spc(20)='ebu_in_c10h16'
             ch_fire_spc(21)='ebu_in_ch3cooh'
             ch_fire_spc(22)='ebu_in_cres'
             ch_fire_spc(23)='ebu_in_glyald'
             ch_fire_spc(24)='ebu_in_mgly'
             ch_fire_spc(25)='ebu_in_gly'
             ch_fire_spc(26)='ebu_in_acetol'
             ch_fire_spc(27)='ebu_in_isop'
             ch_fire_spc(28)='ebu_in_macr'
             ch_fire_spc(29)='ebu_in_mvk'
             ch_fire_spc(30)='ebu_in_oc'
             ch_fire_spc(31)='ebu_in_bc'
             ch_fire_spc(32)='ebu_in_sulf'
             ch_fire_spc(33)='ebu_in_pm25'
             ch_fire_spc(34)='ebu_in_pm10'
!
             ch_bio_spc(1)='MSEBIO_ISOP'
!
! Allocate arrays
             allocate(chem_data3d(nx,ny,nz_chem),chem_data2d(nx,ny))
!
! Loop through chem species and perturb
             call random_seed(size=nseed)
             allocate(seed(nseed))
             do isp=1,nseed
                call date_and_time(cdate,ctime,czone,values)
                read(ctime,*) itime
                seed(isp)=idate+int(real(itime)/(pi**real(isp))/real(ens_member))
             enddo
             print *, idate,itime,ens_member,seed
             call random_seed(put=seed)
             deallocate(seed)
!
             if(pert_chem) then
                do isp=1,nchem_spc
!                   print *, 'processing species ',isp
                   icnt=0
10                 continue
                   call random_number(ran1)
                   call random_number(ran2)
                   y_ran=sqrt(-2.*log(ran1))*cos(2.*pi*ran2)   ! y_ran = N(0,1)
                   z_ran=tr_stdev*y_ran+tr_mean                ! z_ran = N(tr_mean,tr_sdev)
!
! bound the tails and require positive definite scaling
                   if((-z_stdn.gt.y_ran .or. z_stdn.lt.y_ran .or. z_ran.le.0.) .and. tr_stdev.gt.0.) then
                      icnt=icnt+1
                      if(icnt.ge.10) then
                         print *, 'Tail cutoff error loop 10'
                         stop
                      endif
                      go to 10
                   endif
                   call get_WRFCHEM_emiss_data(wrfchemi,ch_chem_spc(isp),chem_data3d,nx,ny,nz_chem)
                   chem_data3d(:,:,:)=z_ran*chem_data3d(:,:,:)
                   call put_WRFCHEM_emiss_data(wrfchemi,ch_chem_spc(isp),chem_data3d,nx,ny,nz_chem)
                enddo
             endif
             if(pert_fire) then
                do isp=1,nfire_spc
!                   print *, 'processing species ',isp
                   icnt=0
20                 continue
                   call random_number(ran1)
                   call random_number(ran2)
                   y_ran=sqrt(-2.*log(ran1))*cos(2.*pi*ran2)
                   z_ran=tr_stdev*y_ran+tr_mean
!
! bound the tails and require positive definite scaling
                   if((-z_stdn.gt.y_ran .or. z_stdn.lt.y_ran .or. z_ran.le.0.) .and. tr_stdev.gt.0) then
                      icnt=icnt+1
                      if(icnt.ge.10) then
                         print *, 'Tail cutoff error loop 20'
                         stop
                      endif
                      go to 20
                   endif
                   call get_WRFCHEM_emiss_data(wrffirechemi,ch_fire_spc(isp),chem_data2d,nx,ny,1)
                   chem_data2d(:,:)=z_ran*chem_data2d(:,:)
                   call put_WRFCHEM_emiss_data(wrffirechemi,ch_fire_spc(isp),chem_data2d,nx,ny,1)
                enddo
             endif
             if(pert_bio) then
                do isp=1,nbio_spc
!                   print *, 'processing species ',isp
                   icnt=0
30                 continue
                   call random_number(ran1)
                   call random_number(ran2)
                   y_ran=sqrt(-2.*log(ran1))*cos(2.*pi*ran2)
                   z_ran=tr_stdev*y_ran+tr_mean
!
! bound the tails and require positive definite scaling
                   if((-z_stdn.gt.y_ran .or. z_stdn.lt.y_ran .or. z_ran.le.0.) .and. tr_stdev.gt.0) then
                      icnt=icnt+1
                      if(icnt.ge.10) then
                         print *, 'Tail cutoff error loop 30'
                         stop
                      endif
                      go to 30
                   endif
                   call get_WRFCHEM_emiss_data(wrfbiochemi,ch_bio_spc(isp),chem_data2d,nx,ny,1)
                   chem_data2d(:,:)=z_ran*chem_data2d(:,:)
                   call put_WRFCHEM_emiss_data(wrfbiochemi,ch_bio_spc(isp),chem_data2d,nx,ny,1)
                enddo
             endif
!
! Deallocate arrays
             deallocate(chem_data3d,chem_data2d)
          end program main
!
          subroutine get_WRFCHEM_emiss_data(file,name,data,nx,ny,nz_chem)
             implicit none
             include 'netcdf.inc'
             integer, parameter                    :: maxdim=6
             integer                               :: nx,ny,nz_chem
             integer                               :: i,rc
             integer                               :: f_id
             integer                               :: v_id,v_ndim,typ,natts
             integer,dimension(maxdim)             :: one
             integer,dimension(maxdim)             :: v_dimid
             integer,dimension(maxdim)             :: v_dim
             real,dimension(nx,ny,nz_chem)         :: data
             character(len=150)                    :: v_nam
             character*(*)                         :: name
             character*(*)                         :: file
!
! open netcdf file
             rc = nf_open(trim(file),NF_NOWRITE,f_id)
!             print *, trim(file)
             if(rc.ne.0) then
                print *, 'nf_open error ',trim(file)
                call abort
             endif
!
! get variables identifiers
             rc = nf_inq_varid(f_id,trim(name),v_id)
!             print *, v_id
             if(rc.ne.0) then
                print *, 'nf_inq_varid error ', v_id
                call abort
             endif
!
! get dimension identifiers
             v_dimid=0
             rc = nf_inq_var(f_id,v_id,v_nam,typ,v_ndim,v_dimid,natts)
!             print *, v_dimid
             if(rc.ne.0) then
                print *, 'nf_inq_var error ', v_dimid
                call abort
             endif
!
! get dimensions
             v_dim(:)=1
             do i=1,v_ndim
                rc = nf_inq_dimlen(f_id,v_dimid(i),v_dim(i))
             enddo
!             print *, v_dim
             if(rc.ne.0) then
                print *, 'nf_inq_dimlen error ', v_dim
                call abort
             endif
!
! check dimensions
             if(nx.ne.v_dim(1)) then
                print *, 'ERROR: nx dimension conflict ',nx,v_dim(1)
                call abort
             else if(ny.ne.v_dim(2)) then
                print *, 'ERROR: ny dimension conflict ',ny,v_dim(2)
                call abort
             else if(nz_chem.ne.v_dim(3)) then             
                print *, 'ERROR: nz_chem dimension conflict ',nz_chem,v_dim(3)
                call abort
             else if(1.ne.v_dim(4)) then             
                print *, 'ERROR: time dimension conflict ',1,v_dim(4)
                call abort
             endif
!
! get data
             one(:)=1
             rc = nf_get_vara_real(f_id,v_id,one,v_dim,data)
             if(rc.ne.0) then
                print *, 'nf_get_vara_real ', data(1,1,1)
                call abort
             endif
             rc = nf_close(f_id)
             return
          end subroutine get_WRFCHEM_emiss_data   
!
          subroutine put_WRFCHEM_emiss_data(file,name,data,nx,ny,nz_chem)
             implicit none
             include 'netcdf.inc'
             integer, parameter                    :: maxdim=6
             integer                               :: nx,ny,nz_chem
             integer                               :: i,rc
             integer                               :: f_id
             integer                               :: v_id,v_ndim,typ,natts
             integer,dimension(maxdim)             :: one
             integer,dimension(maxdim)             :: v_dimid
             integer,dimension(maxdim)             :: v_dim
             real,dimension(nx,ny,nz_chem)         :: data
             character(len=150)                    :: v_nam
             character*(*)                         :: name
             character*(*)                         :: file
!
! open netcdf file
             rc = nf_open(trim(file),NF_WRITE,f_id)
             if(rc.ne.0) then
                print *, 'nf_open error ',trim(file)
                call abort
             endif
!
! get variables identifiers
             rc = nf_inq_varid(f_id,trim(name),v_id)
!             print *, v_id
             if(rc.ne.0) then
                print *, 'nf_inq_varid error ', v_id
                call abort
             endif
!
! get dimension identifiers
             v_dimid=0
             rc = nf_inq_var(f_id,v_id,v_nam,typ,v_ndim,v_dimid,natts)
!             print *, v_dimid
             if(rc.ne.0) then
                print *, 'nf_inq_var error ', v_dimid
                call abort
             endif
!
! get dimensions
             v_dim(:)=1
             do i=1,v_ndim
                rc = nf_inq_dimlen(f_id,v_dimid(i),v_dim(i))
             enddo
!             print *, v_dim
             if(rc.ne.0) then
                print *, 'nf_inq_dimlen error ', v_dim
                call abort
             endif
!
! check dimensions
             if(nx.ne.v_dim(1)) then
                print *, 'ERROR: nx dimension conflict ',nx,v_dim(1)
                call abort
             else if(ny.ne.v_dim(2)) then
                print *, 'ERROR: ny dimension conflict ',ny,v_dim(2)
                call abort
             else if(nz_chem.ne.v_dim(3)) then             
                print *, 'ERROR: nz_chem dimension conflict ',nz_chem,v_dim(3)
                call abort
             else if(1.ne.v_dim(4)) then             
                print *, 'ERROR: time dimension conflict ',1,v_dim(4)
                call abort
             endif
!
! put data
             one(:)=1
             rc = nf_put_vara_real(f_id,v_id,one,v_dim,data)
             if(rc.ne.0) then
                print *, 'nf_put_vara_real ', data(1,1,1)
                call abort
             endif
             rc = nf_close(f_id)
             return
          end subroutine put_WRFCHEM_emiss_data   




