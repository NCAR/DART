MODULE module_model_utils

! module for  wrf1d driver to be compatible with the DART code
! based on model_mod.f90 

! DART modules
  USE types_mod, ONLY : r8, missing_r8
  USE time_manager_mod, ONLY : time_type, set_time, get_time, &
       increment_time, print_time, &
       set_calendar_type, NO_CALENDAR, &
       OPERATOR(==), OPERATOR(<=)
  USE utilities_mod, ONLY : file_exist, open_file, &
       &close_file,register_module, error_handler, E_ERR, E_MSG,&
       &logfileunit
  
! WRF1d modules

  USE module_model_constants
  USE module_initialize
  USE module_namelist
  USE module_wrf

  IMPLICIT NONE

  CHARACTER(len=100) :: source='WRF1d', revision='1', revdate='22.12.2004'

  INTEGER :: wrf_rnd_seed
  TYPE(time_type) :: time_step
  TYPE(time_type) :: last_advance_time
  LOGICAL :: allocate_wrf = .TRUE.

CONTAINS
  
  SUBROUTINE static_init_model()
!------------------------------------------------------------------
! Initializes class data for this model. For now, simply outputs the
! identity info, sets the location of the state variables, and initializes
! the time type for the time stepping (is this general enough for time???)
    
    IMPLICIT NONE
    
    LOGICAL :: is_it_there = .FALSE.
    INTEGER :: unit_nml
    
    OPEN(logfileunit,file='wrf_log.out')

    ! Begin by reading the namelist input
    IF(file_exist('wrf1d_namelist.input')) THEN
       unit_nml = open_file(fname = 'wrf1d_namelist.input', &
            &action = 'read')
       CALL do_namelist_wrf1d(unit_nml,logfileunit)
       CLOSE(unit_nml)
    ENDIF

  
    CLOSE(logfileunit)

! initialize the random sequence
    wrf_rnd_seed = rnd_seed_val
  
! need some dimension information
    CALL static_init_wrf()

! initialize some timing stuff
    time_step = set_time(INT(dt), 0)
    
! initialize the last time
    last_advance_time = set_time(start_forecast*3600,0)
    
  END SUBROUTINE static_init_model
  
  SUBROUTINE init_conditions()
!------------------------------------------------------------------
! subroutine init_conditions
!
! Initial conditions for PBL 1D model is achived via module_wrf
    IMPLICIT NONE

    CALL init_wrf(allocate_wrf,wrf_rnd_seed)
    
    allocate_wrf = .FALSE.
    
! count the ensemble members
  END SUBROUTINE init_conditions
  
  SUBROUTINE adv_1step(dart_time)
!------------------------------------------------------------------
! subroutine adv_1step(time)
!
! Does single time step advance for PBL model - this is a call to wrf.F
! The Time argument is needed for compatibility with more complex models
! that need to know the time to compute their time tendency and is not
! used in L04. Is there a better way to do this in F90 than to just hang
! this argument out everywhere?
    IMPLICIT NONE
    
    TYPE(time_type), INTENT(inout) :: dart_time
   
    INTEGER                  :: dart_seconds,dart_days
   
    CALL get_time(dart_time,dart_seconds,dart_days)
    
    CALL wrf(dart_seconds,dart_days)
    
    last_advance_time = dart_time
    
  END SUBROUTINE adv_1step
  

  SUBROUTINE init_time(time)
!------------------------------------------------------------------
!
! Gets the initial time for a state from the model. Where should this info
! come from in the most general case?
    
    TYPE(time_type), INTENT(out) :: time
    
    INTEGER                      :: start_seconds
    
    start_seconds = start_forecast * 3600
! based on namelist
    time = set_time(start_seconds, 0)
    
  END SUBROUTINE init_time
  

  FUNCTION get_model_time_step()
!------------------------------------------------------------------
! function get_model_time_step()
!
! Returns the the time step of the model. In the long run should be repalced
! by a more general routine that returns details of a general time-stepping
! capability.
    
    TYPE(time_type) :: get_model_time_step
    
    get_model_time_step = time_step
    
  END FUNCTION get_model_time_step
  
  
  
  FUNCTION nc_write_model_atts( ncFileID ) RESULT (ierr)
!------------------------------------------------------------------
! Writes the model-specific attributes to a netCDF file
! TJH Jan 24 2003
!
! TJH 29 July 2003 -- for the moment, all errors are fatal, so the
! return code is always '0 == normal', since the fatal errors stop execution.
!
! JPH 22 Nov 04 -- for the PBL model, this is a whole bunch of info
! about specific model configs
!
! assim_model_mod:init_diag_output uses information from the location_mod
!     to define the location dimension and variable ID. All we need to do
!     is query, verify, and fill ...
!
! Typical sequence for adding new dimensions,variables,attributes:
! NF90_OPEN             ! open existing netCDF dataset
!    NF90_redef         ! put into define mode
!    NF90_def_dim       ! define additional dimensions (if any)
!    NF90_def_var       ! define variables: from name, type, and dims
!    NF90_put_att       ! assign attribute values
! NF90_ENDDEF           ! end definitions: leave define mode
!    NF90_put_var       ! provide values for variable
! NF90_CLOSE            ! close: save updated netCDF dataset

    USE typeSizes
    USE netcdf

    INTEGER, INTENT(in)  :: ncFileID      ! netCDF file identifier
    INTEGER              :: ierr          ! return value of function

!--------------------------------------------------------------------
! General netCDF variables
!--------------------------------------------------------------------
    
    INTEGER :: nDimensions, nVariables, nAttributes, unlimitedDimID
    
!--------------------------------------------------------------------
! netCDF variables for Location
!--------------------------------------------------------------------
    
    INTEGER :: LocationVarID
    INTEGER :: StateVarDimID, StateVarVarID
    INTEGER :: StateVarID, MemberDimID, TimeDimID
    
!------------------------------------------
! same for physical space
!------------------------------------------
    INTEGER :: ZVarVarID, SLVarVarID
    INTEGER :: ZVarDimID, SLVarDimID
    INTEGER :: ZVarID(6), SLVarID, ScrVarID(4)
! U, V, Z, T, QV, TKE
    
!--------------------------------------------------------------------
! local variables
!--------------------------------------------------------------------

    CHARACTER(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
    CHARACTER(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
    CHARACTER(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
    INTEGER, DIMENSION(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
    CHARACTER(len=NF90_MAX_NAME) :: str1,str2
    
    INTEGER             :: i
    ierr = 0                      ! assume normal termination
    
!--------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file
!--------------------------------------------------------------------
    
    CALL check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))
    CALL check(nf90_sync(ncFileID)) ! Ensure netCDF file is current
    CALL check(nf90_Redef(ncFileID))
    
!--------------------------------------------------------------------
! Determine ID's from stuff already in the netCDF file
!--------------------------------------------------------------------
    
! make sure time is unlimited dimid
    
    CALL check(nf90_inq_dimid(ncFileID,"copy",dimid=MemberDimID))
    CALL check(nf90_inq_dimid(ncFileID,"time",dimid=TimeDimID))
    
!--------------------------------------------------------------------
! Write Global Attributes
!--------------------------------------------------------------------
    CALL DATE_AND_TIME(crdate,crtime,crzone,values)
    WRITE(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
         values(1), values(2), values(3), values(5), values(6), values(7)
    
    CALL check(nf90_put_att(ncFileID, NF90_GLOBAL, "creation_date",str1))
    CALL check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_source", source ))
    CALL check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revision", revision ))
    CALL check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revdate", revdate ))
    CALL check(nf90_put_att(ncFileID, NF90_GLOBAL, "model", "PBL_1D"))
    CALL check(nf90_put_att(ncFileID, NF90_GLOBAL, "PBL_type", pbltype ))
    CALL check(nf90_put_att(ncFileID, NF90_GLOBAL, "sfc_type", sfctype ))
    CALL check(nf90_put_att(ncFileID, NF90_GLOBAL, "surface_physics", surface_physics ))
    CALL check(nf90_put_att(ncFileID, NF90_GLOBAL, "dt", dt ))
    CALL check(nf90_put_att(ncFileID, NF90_GLOBAL, "dx", dx ))
    CALL check(nf90_put_att(ncFileID, NF90_GLOBAL, "deep_soil_moisture", deep_soil_moisture ))
    CALL check(nf90_put_att(ncFileID, NF90_GLOBAL, "mavail_ref", mavail_ref ))
    CALL check(nf90_put_att(ncFileID, NF90_GLOBAL, "ivgtyp_ref", ivgtyp_ref ))
    CALL check(nf90_put_att(ncFileID, NF90_GLOBAL, "isltyp_ref", isltyp_ref ))
    CALL check(nf90_put_att(ncFileID, NF90_GLOBAL, "zo_ref", zo_ref ))
    CALL check(nf90_put_att(ncFileID, NF90_GLOBAL, "emiss_ref", emiss_ref ))
    CALL check(nf90_put_att(ncFileID, NF90_GLOBAL, "thc_ref", thc_ref ))
    CALL check(nf90_put_att(ncFileID, NF90_GLOBAL, "albedo_ref", albedo_ref ))
    CALL check(nf90_put_att(ncFileID, NF90_GLOBAL, "ts_ref", ts_ref ))
    CALL check(nf90_put_att(ncFileID, NF90_GLOBAL, "tmn_ref", tmn_ref ))
    CALL check(nf90_put_att(ncFileID, NF90_GLOBAL, "ps_ref", ps_ref ))
    CALL check(nf90_put_att(ncFileID, NF90_GLOBAL, "lat_ref", lat_ref ))
    CALL check(nf90_put_att(ncFileID, NF90_GLOBAL, "lon_ref", lon_ref ))
    CALL check(nf90_put_att(ncFileID, NF90_GLOBAL, "julday", julday ))
    
    
!--------------------------------------------------------------------
! Define the model size, state variable dimension ... whatever ...
!--------------------------------------------------------------------

    CALL check(nf90_def_dim(ncid=ncFileID, name="z_level", &
         len=nz, dimid = ZVarDimID))
    CALL check(nf90_def_dim(ncid=ncFileID, name="sl_level", &
            len=num_soil_layers, dimid = SLVarDimID))
    
!--------------------------------------------------------------------
! Define the Location Variable and add Attributes
! Some of the atts come from location_mod (via the USE: stmnt)
! CF standards for Locations:
! http://www.cgd.ucar.edu/cms/eaton/netcdf/CF-working.html#ctype
!--------------------------------------------------------------------
    
!--------------------------------------------------------------------
! Define either the "state vector" variables -OR- the "prognostic" variables.
!--------------------------------------------------------------------
    
    CALL check(NF90_def_var(ncFileID, name="z_level", xtype=nf90_int, &
         dimids = ZVarDimID, varid=ZVarVarID) )
    CALL check(NF90_def_var(ncFileID, name="sl_level", xtype=nf90_int, &
         dimids = SLVarDimID, varid=SLVarVarID) )
    CALL check(nf90_def_var(ncid=ncFileID, name="U", xtype=nf90_double, &
         dimids = (/ ZVarDimID, MemberDimID, TimeDimID /), varid=ZVarID(1)))
    CALL check(nf90_put_att(ncFileID, ZVarID(1), "long_name", "Zonal Wind"))
    CALL check(nf90_put_att(ncFileID, ZVarID(1), "units", "m/s"))
    CALL check(nf90_def_var(ncid=ncFileID, name="V", xtype=nf90_double, &
         dimids = (/ ZVarDimID, MemberDimID, TimeDimID /), varid=ZVarID(2)))
    CALL check(nf90_put_att(ncFileID, ZVarID(2), "long_name", "Meridional Wind"))
    CALL check(nf90_put_att(ncFileID, ZVarID(2), "units", "m/s"))
    CALL check(nf90_def_var(ncid=ncFileID, name="Z", xtype=nf90_double, &
         dimids = (/ ZVarDimID, MemberDimID, TimeDimID /), varid=ZVarID(3)))
    CALL check(nf90_put_att(ncFileID, ZVarID(3), "long_name", "Height AMSL"))
    CALL check(nf90_put_att(ncFileID, ZVarID(3), "units", "m"))
    CALL check(nf90_def_var(ncid=ncFileID, name="T", xtype=nf90_double, &
         dimids = (/ ZVarDimID, MemberDimID, TimeDimID /), varid=ZVarID(4)))
    CALL check(nf90_put_att(ncFileID, ZVarID(4), "long_name", "Temperature"))
    CALL check(nf90_put_att(ncFileID, ZVarID(4), "units", "K"))
    CALL check(nf90_def_var(ncid=ncFileID, name="QV", xtype=nf90_double, &
         dimids = (/ ZVarDimID, MemberDimID, TimeDimID /), varid=ZVarID(5)))
    CALL check(nf90_put_att(ncFileID, ZVarID(5), "long_name", "Mixing Ratio"))
    CALL check(nf90_put_att(ncFileID, ZVarID(5), "units", "kg/kg"))
    CALL check(nf90_def_var(ncid=ncFileID, name="TKE", xtype=nf90_double, &
         dimids = (/ ZVarDimID, MemberDimID, TimeDimID /), varid=ZVarID(6)))
    CALL check(nf90_put_att(ncFileID, ZVarID(6), "long_name", "Turb. Kinetic Energy"))
    CALL check(nf90_put_att(ncFileID, ZVarID(6), "units", "m2/s2"))
    
    CALL check(nf90_def_var(ncid=ncFileID, name="TSLB", xtype=nf90_double, &
         dimids = (/ SLVarDimID, MemberDimID, TimeDimID /), varid=SLVarID))
    CALL check(nf90_put_att(ncFileID, SLVarID, "long_name", "Soil Temperature"))
    CALL check(nf90_put_att(ncFileID, SLVarID, "units", "K"))
    
    CALL check(nf90_def_var(ncid=ncFileID, name="U10", xtype=nf90_double, &
         dimids = (/ MemberDimID, TimeDimID /), varid=ScrVarID(1)))
    CALL check(nf90_put_att(ncFileID, ScrVarID(1), "long_name", "10-m U-wind"))
    CALL check(nf90_put_att(ncFileID, ScrVarID(1), "units", "m/s"))
    CALL check(nf90_def_var(ncid=ncFileID, name="V10", xtype=nf90_double, &
         dimids = (/ MemberDimID, TimeDimID /), varid=ScrVarID(2)))
    CALL check(nf90_put_att(ncFileID, ScrVarID(2), "long_name", "10-m V-wind"))
    CALL check(nf90_put_att(ncFileID, ScrVarID(2), "units", "m/s"))
    CALL check(nf90_def_var(ncid=ncFileID, name="T2", xtype=nf90_double, &
         dimids = (/ MemberDimID, TimeDimID /), varid=ScrVarID(3)))
    CALL check(nf90_put_att(ncFileID, ScrVarID(3), "long_name", "2-m T"))
    CALL check(nf90_put_att(ncFileID, ScrVarID(3), "units", "m/s"))
    CALL check(nf90_def_var(ncid=ncFileID, name="Q2", xtype=nf90_double, &
         dimids = (/ MemberDimID, TimeDimID /), varid=ScrVarID(4)))
    CALL check(nf90_put_att(ncFileID, ScrVarID(4), "long_name", "2-m Mixing Ratio"))
    CALL check(nf90_put_att(ncFileID, ScrVarID(4), "units", "kg/kg"))
 
! Leave define mode so we can fill
    CALL check(nf90_enddef(ncfileID))
    
! Fill the state variable coordinate variable
    CALL check(nf90_put_var(ncFileID, ZVarVarID, (/ (i,i=1,nz) /) ))
    CALL check(nf90_put_var(ncFileID, SLVarVarID, (/ (i,i=1,num_soil_layers) /) ))
    
!--------------------------------------------------------------------
! Fill the location variable
!--------------------------------------------------------------------
    
    
!--------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!--------------------------------------------------------------------
    CALL check(nf90_sync(ncFileID))
    
    WRITE (*,*)'Model attributes written, netCDF file synched ...'
    
  CONTAINS
! Internal subroutine - checks error status after each netcdf, prints
!                       text message each time an error code is returned.
    SUBROUTINE check(istatus)
      INTEGER, INTENT ( in) :: istatus
      IF(istatus /= nf90_noerr) CALL error_handler(E_ERR,'nc_write_model_atts',&
           TRIM(nf90_strerror(istatus)), source, revision, revdate)
    END SUBROUTINE check
  END FUNCTION nc_write_model_atts
  
  FUNCTION nc_write_model_vars( ncFileID, statevec, copyindex, timeindex ) RESULT (ierr)         
!------------------------------------------------------------------
! Writes the model-specific attributes to a netCDF file
! TJH 23 May 2003
!
! TJH 29 July 2003 -- for the moment, all errors are fatal, so the
! return code is always '0 == normal', since the fatal errors stop execution.
!
! For the lorenz_04 model, each state variable is at a separate location.
! that's all the model-specific attributes I can think of ...
!
! assim_model_mod:init_diag_output uses information from the location_mod
!     to define the location dimension and variable ID. All we need to do
!     is query, verify, and fill ...
!
! Typical sequence for adding new dimensions,variables,attributes:
! NF90_OPEN             ! open existing netCDF dataset
!    NF90_redef         ! put into define mode
!    NF90_def_dim       ! define additional dimensions (if any)
!    NF90_def_var       ! define variables: from name, type, and dims
!    NF90_put_att       ! assign attribute values
! NF90_ENDDEF           ! end definitions: leave define mode
!    NF90_put_var       ! provide values for variable
! NF90_CLOSE            ! close: save updated netCDF dataset
    
    USE typeSizes
    USE netcdf
    
    INTEGER,                INTENT(in) :: ncFileID      ! netCDF file identifier
    REAL(r8), DIMENSION(:), INTENT(in) :: statevec
    INTEGER,                INTENT(in) :: copyindex
    INTEGER,                INTENT(in) :: timeindex
    INTEGER                            :: ierr          ! return value of function
    
!-------------------------------------------------------------------------------
! General netCDF variables
!-------------------------------------------------------------------------------
    
    INTEGER :: nDimensions, nVariables, nAttributes, unlimitedDimID
    INTEGER :: StateVarID, ZVarID(6), ScrVarID(4),SLVarID
    
!-------------------------------------------------------------------------------
! local variables
!-------------------------------------------------------------------------------

    ierr = 0                      ! assume normal termination
    
!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file
!-------------------------------------------------------------------------------
    
    CALL check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))
    
    
    CALL check(NF90_inq_varid(ncFileID, "U", ZVarID(1)) )
    CALL check(NF90_put_var(ncFileID, ZVarID(1), u_phy(1,:,1),  &
         start=(/ 1, copyindex, timeindex /)))
    CALL check(NF90_inq_varid(ncFileID, "V", ZVarID(2)) )
    CALL check(NF90_put_var(ncFileID, ZVarID(2), v_phy(1,:,1),  &
         start=(/ 1, copyindex, timeindex /)))
    CALL check(NF90_inq_varid(ncFileID, "Z", ZVarID(3)) )
    CALL check(NF90_put_var(ncFileID, ZVarID(3), z(1,:,1),  &
         start=(/ 1, copyindex, timeindex /)))
    CALL check(NF90_inq_varid(ncFileID, "T", ZVarID(4)) )
    CALL check(NF90_put_var(ncFileID, ZVarID(4), t_phy(1,:,1),  &
         start=(/ 1, copyindex, timeindex /)))
    CALL check(NF90_inq_varid(ncFileID, "QV", ZVarID(5)) )
    CALL check(NF90_put_var(ncFileID, ZVarID(5), moist(1,:,1,1),  &
         start=(/ 1, copyindex, timeindex /)))
    CALL check(NF90_inq_varid(ncFileID, "TKE", ZVarID(6)) )
    CALL check(NF90_put_var(ncFileID, ZVarID(6), tke_myj(1,:,1),  &
         start=(/ 1, copyindex, timeindex /)))
    
    CALL check(NF90_inq_varid(ncFileID, "TSLB", SLVarID) )
    CALL check(NF90_put_var(ncFileID, SLVarID, tslb(1,:,1),  &
         start=(/ 1, copyindex, timeindex /)))
    
    CALL check(NF90_inq_varid(ncFileID, "U10", ScrVarID(1)) )
    CALL check(NF90_put_var(ncFileID, ScrVarID(1), u10(1,1),  &
         start=(/ copyindex, timeindex /)))
    CALL check(NF90_inq_varid(ncFileID, "V10", ScrVarID(2)) )
    CALL check(NF90_put_var(ncFileID, ScrVarID(2), v10(1,1),  &
         start=(/ copyindex, timeindex /)))
    CALL check(NF90_inq_varid(ncFileID, "T2", ScrVarID(3)) )
    CALL check(NF90_put_var(ncFileID, ScrVarID(3), t2(1,1),  &
         start=(/ copyindex, timeindex /)))
    CALL check(NF90_inq_varid(ncFileID, "Q2", ScrVarID(4)) )
    CALL check(NF90_put_var(ncFileID, ScrVarID(4), q2(1,1),  &
         start=(/ copyindex, timeindex /)))

! write (*,*)'Finished filling variables ...'
    CALL check(nf90_sync(ncFileID))
! write (*,*)'netCDF file is synched ...'
    
  CONTAINS
! Internal subroutine - checks error status after each netcdf, prints
!                       text message each time an error code is returned.
    SUBROUTINE check(istatus)
      INTEGER, INTENT ( in) :: istatus
      IF(istatus /= nf90_noerr) CALL error_handler(E_ERR,'nc_write_model_vars',&
           TRIM(nf90_strerror(istatus)), source, revision, revdate)
    END SUBROUTINE check
  END FUNCTION nc_write_model_vars
  
END MODULE module_model_utils
  
