MODULE module_soil_pre

CONTAINS

   SUBROUTINE process_percent_cat ( xland , &
                                landuse_frac_input , soil_top_cat_input , soil_bot_cat_input , &
                                isltyp_input , ivgtyp_input , &
                                ids , ide , jds , jde , kds , kde , &
                                ims , ime , jms , jme , kms , kme , &
                                its , ite , jts , jte , kts , kte , &
                                iswater )

      IMPLICIT NONE


      INTEGER , INTENT(IN) :: ids , ide , jds , jde , kds , kde , &
                              ims , ime , jms , jme , kms , kme , &
                              its , ite , jts , jte , kts , kte , &
                              iswater
      REAL , DIMENSION(:,:,:) , INTENT(IN):: landuse_frac_input , soil_top_cat_input , soil_bot_cat_input
      INTEGER , DIMENSION(:,:), INTENT(OUT) :: isltyp_input , ivgtyp_input
      REAL , DIMENSION(ims:ime,jms:jme) , INTENT(OUT) :: xland

      INTEGER :: i , j , l , dominant_index , num_soil_cat , num_veg_cat
      REAL :: dominant_value

      INTEGER , PARAMETER :: iswater_soil = 14
      INTEGER :: iforce

      num_veg_cat = SIZE ( landuse_frac_input , DIM=3 )

      !  Compute the dominant VEGETATION INDEX.

      DO i = 1 , ite
         DO j = 1 , jte
            dominant_value = landuse_frac_input(i,j,1)
            dominant_index = 1
            DO l = 2 , num_veg_cat
               IF      ( ( l .EQ. iswater ) .AND. ( landuse_frac_input(i,j,l) .GT.            0.5 ) ) THEN
                  dominant_value = soil_top_cat_input(i,j,l)
                  dominant_index = l
               ELSE IF ( ( l .NE. iswater ) .AND. ( landuse_frac_input(i,j,l) .GT. dominant_value ) ) THEN
                  dominant_value = landuse_frac_input(i,j,l)
                  dominant_index = l
               END IF
            END DO
            IF      ( dominant_index .EQ. iswater ) THEN
               xland(i,j) = 2.
            ELSE IF ( dominant_index .NE. iswater ) THEN
               xland(i,j) = 1.
            END IF
            ivgtyp_input(i,j) = dominant_index
         END DO
      END DO

      num_soil_cat = SIZE ( soil_top_cat_input , DIM=3 )

      !  Compute the dominant SOIL TEXTURE INDEX, TOP.

      iforce = 0
      DO i = 1 , ite
         DO j = 1 , jte
            dominant_value = soil_top_cat_input(i,j,1)
            dominant_index = 1
            IF ( xland(i,j) .LT. 1.5 ) THEN
               DO l = 2 , num_soil_cat
                  IF ( ( l .NE. iswater_soil ) .AND. ( soil_top_cat_input(i,j,l) .GT. dominant_value ) ) THEN
                     dominant_value = soil_top_cat_input(i,j,l)
                     dominant_index = l
                  END IF
               END DO
               IF ( dominant_value .LT. 0.01 ) THEN
                  iforce = iforce + 1
!print *,'forcing a soil value over land'
!print *,iforce,NINT(soil_top_cat_input(i,j,:))
                  dominant_index = 8
               END IF
            ELSE
               dominant_index = iswater_soil
            END IF
            isltyp_input(i,j) = dominant_index
         END DO
      END DO
if(iforce.ne.0)then
print *,'forcing artificial silty clay loam at ',iforce,' points, out of ',(ite)*(jte)
endif

   END SUBROUTINE process_percent_cat

   SUBROUTINE process_soil_real ( tsk , tmn , xland , &
                                landmask_input , sst_input , &
                                st_input , sm_input , st_levels_input , sm_levels_input , &
                                zs , dzs , tslb , smois , &
                                flag_sst , &
                                st000010_input , st010040_input , st040100_input , st100200_input , &
                                st010200_input , &
                                sm000010_input , sm010040_input , sm040100_input , sm100200_input , &
                                sm010200_input , &
                                ids , ide , jds , jde , kds , kde , &
                                ims , ime , jms , jme , kms , kme , &
                                its , ite , jts , jte , kts , kte , &
                                bl_surface_physics , num_soil_layers , real_data_init_type , &
                                num_st_levels_input , num_sm_levels_input )

      IMPLICIT NONE

      INTEGER , INTENT(IN) :: ids , ide , jds , jde , kds , kde , &
                              ims , ime , jms , jme , kms , kme , &
                              its , ite , jts , jte , kts , kte , &
                              bl_surface_physics , num_soil_layers , real_data_init_type , &
                              num_st_levels_input , num_sm_levels_input

      LOGICAL , INTENT(IN) :: flag_sst

      REAL , DIMENSION(:,:) , INTENT(IN) :: landmask_input , sst_input

      REAL , DIMENSION(:,:,:) , INTENT(INOUT) :: st_input , sm_input
      INTEGER , DIMENSION(:) , INTENT(INOUT) :: st_levels_input , sm_levels_input

      REAL, DIMENSION(1:num_soil_layers), INTENT(OUT)  ::  zs,dzs
      REAL , DIMENSION(ims:ime,num_soil_layers,jms:jme) , INTENT(OUT) :: tslb , smois

      REAL , DIMENSION(ims:ime,jms:jme) , INTENT(IN) :: tsk , tmn , xland

      REAL , DIMENSION(:,:) , INTENT(IN) :: st000010_input , st010040_input , st040100_input , st100200_input , &
                                            st010200_input , &
                                            sm000010_input , sm010040_input , sm040100_input , sm100200_input , &
                                            sm010200_input 

      !  Initialize the soil depth, and the soil temperature and moisture.
   
      IF      ( ( bl_surface_physics .EQ. 1 ) .AND. ( num_soil_layers .GT. 1 ) ) THEN
         CALL init_soil_depth_1 ( zs , dzs , num_soil_layers )
         CALL init_soil_1_real ( tsk , tmn , tslb , zs , dzs , num_soil_layers , real_data_init_type , &
                                 landmask_input , sst_input , flag_sst , &
                                 ids , ide , jds , jde , kds , kde , &
                                 ims , ime , jms , jme , kms , kme , &
                                 its , ite , jts , jte , kts , kte )

      ELSE IF ( ( bl_surface_physics .EQ. 2 ) .AND. ( num_soil_layers .GT. 1 ) ) THEN
         CALL init_soil_depth_2 ( zs , dzs , num_soil_layers )
         CALL init_soil_2_real ( tsk , tmn , smois , tslb , &
                                 st_input , sm_input , landmask_input , sst_input , &
                                 zs , dzs , &
                                 st_levels_input , sm_levels_input , &
                                 num_soil_layers , num_st_levels_input , num_sm_levels_input , &
                                 flag_sst , &
                                 ids , ide , jds , jde , kds , kde , &
                                 ims , ime , jms , jme , kms , kme , &
                                 its , ite , jts , jte , kts , kte )
!        CALL init_soil_old ( tsk , tmn , &
!                             smois , tslb , zs , dzs , num_soil_layers , &
!                             st000010_input , st010040_input , st040100_input , st100200_input , &
!                             st010200_input , &
!                             sm000010_input , sm010040_input , sm040100_input , sm100200_input , &
!                             sm010200_input , &
!                             landmask_input , sst_input , &
!                             ids , ide , jds , jde , kds , kde , &
!                             ims , ime , jms , jme , kms , kme , &
!                             its , ite , jts , jte , kts , kte )
      ELSE IF ( ( bl_surface_physics .EQ. 3 ) .AND. ( num_soil_layers .GT. 1 ) ) THEN
         CALL init_soil_depth_3 ( zs , dzs , num_soil_layers )
         CALL init_soil_3_real ( tsk , tmn , smois , tslb , &
                                 st_input , sm_input , landmask_input , sst_input , &
                                 zs , dzs , &
                                 st_levels_input , sm_levels_input , &
                                 num_soil_layers , num_st_levels_input , num_sm_levels_input , &
                                 flag_sst , &
                                 ids , ide , jds , jde , kds , kde , &
                                 ims , ime , jms , jme , kms , kme , &
                                 its , ite , jts , jte , kts , kte )
      END IF

   END SUBROUTINE process_soil_real

   SUBROUTINE process_soil_ideal ( xland,xice,vegfra,snow,canwat,  &
                                   ivgtyp,isltyp,tslb,smois, &
                                   tsk,tmn,zs,dzs,           &
                                   num_soil_layers,          &
                                   bl_surface_physics ,      &
                                   ids,ide, jds,jde, kds,kde,&
                                   ims,ime, jms,jme, kms,kme,&
                                   its,ite, jts,jte, kts,kte )

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::ids,ide, jds,jde, kds,kde,  &
                            ims,ime, jms,jme, kms,kme,  &
                            its,ite, jts,jte, kts,kte
  
      INTEGER, INTENT(IN) :: num_soil_layers , bl_surface_physics

      REAL, DIMENSION( ims:ime, num_soil_layers, jms:jme ) , INTENT(INOUT) :: smois, tslb

      REAL, DIMENSION(num_soil_layers), INTENT(OUT) :: dzs,zs

      REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT) :: tsk, tmn
      REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(OUT) :: xland, snow, canwat, xice, vegfra
      INTEGER, DIMENSION( ims:ime, jms:jme ) , INTENT(OUT) :: ivgtyp, isltyp

      !  Local variables.

      INTEGER :: itf,jtf

      itf=MIN(ite,ide)
      jtf=MIN(jte,jde)

      IF      ( ( bl_surface_physics .EQ. 1 ) .AND. ( num_soil_layers .GT. 1 ) ) THEN
         CALL init_soil_depth_1 ( zs , dzs , num_soil_layers )
         CALL init_soil_1_ideal(tsk,tmn,tslb,xland,ivgtyp,zs,dzs,num_soil_layers,     &
                                ids,ide, jds,jde, kds,kde,               &
                                ims,ime, jms,jme, kms,kme,               &
                                its,ite, jts,jte, kts,kte                )
      ELSE IF ( ( bl_surface_physics .EQ. 2 ) .AND. ( num_soil_layers .GT. 1 ) ) THEN
         CALL init_soil_depth_2 ( zs , dzs , num_soil_layers )
         CALL init_soil_2_ideal ( xland,xice,vegfra,snow,canwat,         &
                                  ivgtyp,isltyp,tslb,smois,tmn,          &
                                  num_soil_layers,                       &
                                  ids,ide, jds,jde, kds,kde,             &
                                  ims,ime, jms,jme, kms,kme,             &
                                  its,ite, jts,jte, kts,kte              )
      END IF

   END SUBROUTINE process_soil_ideal

   SUBROUTINE adjust_soil_temp ( tmn , bl_surface_physics , &
                                 tsk , t_annual_avg_input , ter_input , toposoil_input , &
                                 st000010_input , st010040_input , st040100_input , st100200_input , st010200_input , &
                                 flag_st000010 , flag_st010040 , flag_st040100 , flag_st100200 , flag_st010200 , & 
                                 ids , ide , jds , jde , kds , kde , &
                                 ims , ime , jms , jme , kms , kme , &
                                 its , ite , jts , jte , kts , kte )

      IMPLICIT NONE
   
      INTEGER , INTENT(IN) :: ids , ide , jds , jde , kds , kde , &
                              ims , ime , jms , jme , kms , kme , &
                              its , ite , jts , jte , kts , kte 

      REAL , DIMENSION(:,:) , INTENT(IN) :: ter_input , toposoil_input
      REAL , DIMENSION(ims:ime,jms:jme) , INTENT(INOUT) :: tmn , tsk
      REAL , DIMENSION(:,:) , INTENT(INOUT) :: t_annual_avg_input , &
                              st000010_input , st010040_input , st040100_input , st100200_input , st010200_input
      LOGICAL , INTENT(IN) :: flag_st000010 , flag_st010040 , flag_st040100 , flag_st100200 , flag_st010200
      INTEGER , INTENT(IN) :: bl_surface_physics
 
      INTEGER :: i , j

      IF ( bl_surface_physics .EQ. 1 ) THEN
         DO j = 1 , jde
            DO i = 1 , ide
               tmn(i,j) = tmn(i,j) - 0.0065 * ( ter_input(i,j) - toposoil_input(i,j) )
            END DO
         END DO
      END IF

      DO j = 1 , jde
         DO i = 1 , ide
            tsk(i,j) = tsk(i,j) - 0.0065 * ( ter_input(i,j) - toposoil_input(i,j) )
!           t_annual_avg_input(i,j) = t_annual_avg_input(i,j) - 0.0065 * ter_input(i,j) ! handled by SI
            IF ( flag_st000010 ) THEN
               st000010_input(i,j) = st000010_input(i,j) - 0.0065 * ( ter_input(i,j) - toposoil_input(i,j) )
            END IF
            IF ( flag_st010040 ) THEN
               st010040_input(i,j) = st010040_input(i,j) - 0.0065 * ( ter_input(i,j) - toposoil_input(i,j) )
            END IF
            IF ( flag_st040100 ) THEN
               st040100_input(i,j) = st040100_input(i,j) - 0.0065 * ( ter_input(i,j) - toposoil_input(i,j) )
            END IF
            IF ( flag_st100200 ) THEN
               st100200_input(i,j) = st100200_input(i,j) - 0.0065 * ( ter_input(i,j) - toposoil_input(i,j) )
            END IF
            IF ( flag_st010200 ) THEN
               st010200_input(i,j) = st010200_input(i,j) - 0.0065 * ( ter_input(i,j) - toposoil_input(i,j) )
            END IF
         END DO
      END DO

   END SUBROUTINE adjust_soil_temp

   SUBROUTINE adjust_soil_temp_3 ( tmn , bl_surface_physics , &
                                 tsk , t_annual_avg_input , ter_input , toposoil_input , &
      soilt000_input , soilt005_input , soilt020_input , soilt040_input , soilt160_input , soilt300_input , &
      flag_soilt000 , flag_soilt005 , flag_soilt020 , flag_soilt040 , flag_soilt160 , flag_soilt300 , &
                                 ids , ide , jds , jde , kds , kde , &
                                 ims , ime , jms , jme , kms , kme , &
                                 its , ite , jts , jte , kts , kte )

      IMPLICIT NONE
   
      INTEGER , INTENT(IN) :: ids , ide , jds , jde , kds , kde , &
                              ims , ime , jms , jme , kms , kme , &
                              its , ite , jts , jte , kts , kte 

      REAL , DIMENSION(:,:) , INTENT(IN) :: ter_input , toposoil_input
      REAL , DIMENSION(ims:ime,jms:jme) , INTENT(INOUT) :: tmn , tsk
      REAL , DIMENSION(:,:) , INTENT(INOUT) :: t_annual_avg_input , &
      soilt000_input , soilt005_input , soilt020_input , soilt040_input , soilt160_input , soilt300_input
      LOGICAL , INTENT(IN) :: flag_soilt000 , flag_soilt005 , flag_soilt020 , flag_soilt040 , flag_soilt160 , flag_soilt300
      INTEGER , INTENT(IN) :: bl_surface_physics
 
      INTEGER :: i , j

      DO j = 1 , jde
         DO i = 1 , ide
            tsk(i,j) = tsk(i,j) - 0.0065 * ( ter_input(i,j) - toposoil_input(i,j) )
!           t_annual_avg_input(i,j) = t_annual_avg_input(i,j) - 0.0065 * ter_input(i,j) ! handled by SI
            IF ( flag_soilt000 ) THEN
               soilt000_input(i,j) = soilt000_input(i,j) - 0.0065 * ( ter_input(i,j) - toposoil_input(i,j) )
            END IF
            IF ( flag_soilt005 ) THEN
               soilt005_input(i,j) = soilt005_input(i,j) - 0.0065 * ( ter_input(i,j) - toposoil_input(i,j) )
            END IF
            IF ( flag_soilt020 ) THEN
               soilt020_input(i,j) = soilt020_input(i,j) - 0.0065 * ( ter_input(i,j) - toposoil_input(i,j) )
            END IF
            IF ( flag_soilt040 ) THEN
               soilt040_input(i,j) = soilt040_input(i,j) - 0.0065 * ( ter_input(i,j) - toposoil_input(i,j) )
            END IF
            IF ( flag_soilt160 ) THEN
               soilt160_input(i,j) = soilt160_input(i,j) - 0.0065 * ( ter_input(i,j) - toposoil_input(i,j) )
            END IF
            IF ( flag_soilt300 ) THEN
               soilt300_input(i,j) = soilt300_input(i,j) - 0.0065 * ( ter_input(i,j) - toposoil_input(i,j) )
            END IF
         END DO
      END DO

   END SUBROUTINE adjust_soil_temp_3

   SUBROUTINE init_soil_depth_0 ( zs, dzs, num_soil_layers, thc, cs )

! for force restore - depth based on depth of diurnal variation
! dependent on soil parameters
! d=2/.95 * cg/cs
! cg is thermal capacity (see Garrat p. 227-229)
! cs is volumetric heat capacity
! thc is thermal intertia from the namelist
! alpha is magic number= 3.288e6 converting thermal inertia to thermal capacity
! taken from routine slab
! cg= alpha*thc
     

      IMPLICIT NONE

      REAL, DIMENSION(:,:), INTENT(IN) :: thc, cs
   
      INTEGER, INTENT(IN) :: num_soil_layers
   
      REAL, DIMENSION(1:num_soil_layers), INTENT(OUT)  ::  zs, dzs
   
      REAL :: cg, alpha=3.288e6

      cg=alpha*thc(1,1)

      IF ( num_soil_layers .NE. 1 ) THEN
         PRINT '(A)','Usually, force restore diffusion uses 1 layer.  Change this in the namelist.'
         STOP 'force restore in module_soil_pre.F'
      END IF
   
      dzs(1)=2.11*cg/cs(1,1)
      zs(1)=.5*dzs(1)

   END SUBROUTINE init_soil_depth_0

   SUBROUTINE init_soil_depth_1 ( zs , dzs , num_soil_layers )

      IMPLICIT NONE
   
      INTEGER, INTENT(IN) :: num_soil_layers
   
      REAL, DIMENSION(1:num_soil_layers), INTENT(OUT)  ::  zs,dzs
   
      INTEGER                   ::      l

      !  Define layers (top layer = 0.01 m).  Double the thicknesses at each step (dzs values).
      !  The distance from the ground level to the midpoint of the layer is given by zs.

      !    -------   Ground Level   ----------      ||      ||   ||  || 
      !                                             ||      ||   ||  || zs(1) = 0.005 m
      !    --  --  --  --  --  --  --  --  --       ||      ||   ||  \/
      !                                             ||      ||   ||
      !    -----------------------------------      ||  ||  ||   \/   dzs(1) = 0.01 m
      !                                             ||  ||  || 
      !                                             ||  ||  || zs(2) = 0.02
      !    --  --  --  --  --  --  --  --  --       ||  ||  \/
      !                                             ||  ||
      !                                             ||  ||
      !    -----------------------------------  ||  ||  \/   dzs(2) = 0.02 m
      !                                         ||  || 
      !                                         ||  ||
      !                                         ||  || 
      !                                         ||  || zs(3) = 0.05
      !    --  --  --  --  --  --  --  --  --   ||  \/
      !                                         ||
      !                                         ||
      !                                         ||
      !                                         ||
      !    -----------------------------------  \/   dzs(3) = 0.04 m

      IF ( num_soil_layers .NE. 5 ) THEN
         PRINT '(A)','Usually, the 5-layer diffusion uses 5 layers.  Change this in the namelist.'
         STOP '5-layer_diffusion_uses_5_layers'
      END IF
   
      dzs(1)=.01
      zs(1)=.5*dzs(1)
   
      DO l=2,num_soil_layers
         dzs(l)=2*dzs(l-1)
         zs(l)=zs(l-1)+.5*dzs(l-1)+.5*dzs(l)
      ENDDO

   END SUBROUTINE init_soil_depth_1

   SUBROUTINE init_soil_depth_2 ( zs , dzs , num_soil_layers )

      IMPLICIT NONE
   
      INTEGER, INTENT(IN) :: num_soil_layers
   
      REAL, DIMENSION(1:num_soil_layers), INTENT(OUT)  ::  zs,dzs
   
      INTEGER                   ::      l

      dzs = (/ 0.1 , 0.3 , 0.6 , 1.0 /)

      IF ( num_soil_layers .NE. 4 ) THEN
         PRINT '(A)','Usually, the LSM uses 4 layers.  Change this in the namelist.'
         STOP 'LSM_uses_4_layers'
      END IF

      zs(1)=.5*dzs(1)
   
      DO l=2,num_soil_layers
         zs(l)=zs(l-1)+.5*dzs(l-1)+.5*dzs(l)
      ENDDO

   END SUBROUTINE init_soil_depth_2

   SUBROUTINE init_soil_depth_3 ( zs , dzs , num_soil_layers )

      IMPLICIT NONE
   
      INTEGER, INTENT(IN) :: num_soil_layers
   
      REAL, DIMENSION(1:num_soil_layers), INTENT(OUT)  ::  zs,dzs
   
      zs  = (/ 0.00 , 0.05 , 0.20 , 0.40 , 1.60 , 3.00 /)
      dzs = (/ 0.00 , 0.125, 0.175 , 0.70 , 1.30 , 1.40 /)

      IF ( num_soil_layers .NE. 6 ) THEN
         PRINT '(A)','Usually, the RUC LSM uses 6 layers.  Change this in the namelist.'
         STOP 'LSM_uses_6_layers'
      END IF

   END SUBROUTINE init_soil_depth_3

   SUBROUTINE init_soil_0_real ( tsk , tmn , tslb_init_f, zs , dzs , zs_f, &
        num_soil_layers, ns_f, real_data_init_type , &
        landmask_input , sst_input , flag_sst , &
        ids , ide , jds , jde , kds , kde , &
        ims , ime , jms , jme , kms , kme , &
        its , ite , jts , jte , kts , kte )
     
     IMPLICIT NONE
     
     INTEGER , INTENT(IN) :: num_soil_layers , ns_f,real_data_init_type , &
          ids , ide , jds , jde , kds , kde , &
          ims , ime , jms , jme , kms , kme , &
          its , ite , jts , jte , kts , kte 
     
     LOGICAL , INTENT(IN) :: flag_sst
     
     REAL , DIMENSION(:,:) , INTENT(IN) :: landmask_input , sst_input
     REAL , DIMENSION(ims:ime,jms:jme) :: tsk , tmn
     REAL , DIMENSION(ns_f) ::  tslb_init_f
     REAL , DIMENSION(num_soil_layers) :: zs , dzs
     REAL , DIMENSION(ns_f) :: zs_f
     
     INTEGER :: i , j , l
     
! deep soil linearly interpolated form input data possible
     
     DO j = 1 , jte
        DO i = 1 , ite
           IF ( ( real_data_init_type .EQ. 1 ) .AND. landmask_input(i,j) .GT. 0.5 ) THEN
              IF (zs(1) <= zs_f(1)) THEN
                 tmn(i,j)=tsk(i,j) + (tslb_init_f(1)-tsk(i,j))/zs_f(1)*zs(1)
              ELSEIF (zs(1) > zs_f(ns_f)) THEN
                 tmn(i,j)=tslb_init_f(ns_f)
              ELSE
                 DO l = 1 , ns_f
                    IF (zs_f(l) >= zs(1)) EXIT
                 ENDDO
                 IF (ABS(tslb_init_f(l)-tslb_init_f(l-1)) < .01) THEN 
                 tmn(i,j)=tslb_init_f(l-1)
              ELSE
                 tmn(i,j)=tslb_init_f(l-1)+(zs(1)-zs_f(l-1))/(zs_f(l)-zs_f(l-1))*&
                      &(tslb_init_f(l)-tslb_init_f(l-1))
              ENDIF
              
           ENDIF
           
        ELSE
           IF ( ( real_data_init_type .EQ. 1 ) .AND. ( flag_sst ) ) THEN
              tmn(i,j)= sst_input(i,j)
           END IF
        END IF
     END DO
  END DO
  
END SUBROUTINE init_soil_0_real


   SUBROUTINE init_soil_1_real ( tsk , tmn , tslb , zs , dzs , &
                                 num_soil_layers , real_data_init_type , &
                                 landmask_input , sst_input , flag_sst , &
                                 ids , ide , jds , jde , kds , kde , &
                                 ims , ime , jms , jme , kms , kme , &
                                 its , ite , jts , jte , kts , kte )

      IMPLICIT NONE

      INTEGER , INTENT(IN) :: num_soil_layers , real_data_init_type , &
                              ids , ide , jds , jde , kds , kde , &
                              ims , ime , jms , jme , kms , kme , &
                              its , ite , jts , jte , kts , kte 

      LOGICAL , INTENT(IN) :: flag_sst

      REAL , DIMENSION(:,:) , INTENT(IN) :: landmask_input , sst_input
      REAL , DIMENSION(ims:ime,jms:jme) , INTENT(IN) :: tsk , tmn
      REAL , DIMENSION(num_soil_layers) :: zs , dzs

      REAL , DIMENSION(ims:ime,num_soil_layers,jms:jme) , INTENT(OUT) :: tslb

      INTEGER :: i , j , l

      !  Soil temperature is linearly interpolated between the skin temperature (taken to be at a
      !  depth of 0.5 cm) and the deep soil, annual temperature (taken to be at a depth of 23 cm).
      !  The tslb(i,1,j) is the skin temperature, and the tslb(i,num_soil_layers,j) level is the 
      !  annual mean temperature.

      DO j = 1 , jte
         DO i = 1 , ite

            IF ( landmask_input(i,j) .GT. 0.5 ) THEN
               DO l = 1 , num_soil_layers
                  tslb(i,l,j)= ( tsk(i,j) * ( zs(num_soil_layers) - zs(l) )   + &
                                 tmn(i,j) * ( zs(              l) - zs(1) ) ) / &
                                            ( zs(num_soil_layers) - zs(1) )
               END DO
            ELSE
               IF ( ( real_data_init_type .EQ. 1 ) .AND. ( flag_sst ) ) THEN
                  DO l = 1 , num_soil_layers
                     tslb(i,l,j)= sst_input(i,j)
                  END DO
               ELSE
                  DO l = 1 , num_soil_layers
                     tslb(i,l,j)= tsk(i,j)
                  END DO
               END IF
            END IF
         END DO
      END DO

   END SUBROUTINE init_soil_1_real

   SUBROUTINE init_soil_1_ideal(tsk,tmn,tslb,xland,             &
                       ivgtyp,zs,dzs,num_soil_layers,           &
                       ids,ide, jds,jde, kds,kde,               &
                       ims,ime, jms,jme, kms,kme,               &
                       its,ite, jts,jte, kts,kte                )

      IMPLICIT NONE
   
      INTEGER, INTENT(IN   )    ::      ids,ide, jds,jde, kds,kde, &
                                        ims,ime, jms,jme, kms,kme, &
                                        its,ite, jts,jte, kts,kte
   
      INTEGER, INTENT(IN   )    ::      num_soil_layers
   
      REAL, DIMENSION( ims: , 1: , jms: ), INTENT(OUT) :: tslb
      REAL, DIMENSION( ims:ime , jms:jme ), INTENT(OUT) :: xland
      INTEGER, DIMENSION( ims:ime , jms:jme ), INTENT(OUT) :: ivgtyp
   
      REAL, DIMENSION(1:), INTENT(IN) :: dzs,zs
   
      REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(IN) :: tsk, tmn

      !  Lcal variables.
   
      INTEGER :: l,j,i,itf,jtf
   
      itf=MIN(ite,ide)
      jtf=MIN(jte,jde)
   
      IF (num_soil_layers.NE.1)THEN
         DO j=jts,jtf
            DO l=1,num_soil_layers
               DO i=its,itf
                 tslb(i,l,j)=( tsk(i,j)*(zs(num_soil_layers)-zs(l)) + tmn(i,j)*(zs(l)-zs(1)) ) / &
                             ( zs(num_soil_layers)-zs(1) )
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      DO j=jts,jtf
         DO i=its,itf
           xland(i,j)  = 2
           ivgtyp(i,j) = 7
         ENDDO
      ENDDO

   END SUBROUTINE init_soil_1_ideal


   SUBROUTINE init_soil_2_real ( tsk , tmn , smois , tslb , &
                                 st_input , sm_input , landmask_input , sst_input , &
                                 zs , dzs , &
                                 st_levels_input , sm_levels_input , &
                                 num_soil_layers , num_st_levels_input , num_sm_levels_input ,  &
                                 flag_sst , &
                                 ids , ide , jds , jde , kds , kde , &
                                 ims , ime , jms , jme , kms , kme , &
                                 its , ite , jts , jte , kts , kte )

      IMPLICIT NONE

      INTEGER , INTENT(IN) :: num_soil_layers , num_st_levels_input , num_sm_levels_input , &
                              ids , ide , jds , jde , kds , kde , &
                              ims , ime , jms , jme , kms , kme , &
                              its , ite , jts , jte , kts , kte 

      LOGICAL , INTENT(IN) :: flag_sst

      INTEGER , DIMENSION(:) , INTENT(INOUT) :: st_levels_input , sm_levels_input

      REAL , DIMENSION(:,:,:) , INTENT(INOUT) :: st_input , sm_input
      REAL , DIMENSION(:,:) , INTENT(IN) :: landmask_input , sst_input 

      REAL , DIMENSION(ims:ime,jms:jme) , INTENT(IN) :: tsk , tmn
      REAL , DIMENSION(num_soil_layers) :: zs , dzs

      REAL , DIMENSION(ims:ime,num_soil_layers,jms:jme) , INTENT(OUT) :: tslb , smois

      REAL , ALLOCATABLE , DIMENSION(:) :: zhave

      INTEGER :: i , j , l , lout , lin , lwant , lhave
      REAL :: temp

      !  Allocate the soil layer array used for interpolating.

      IF ( ( num_st_levels_input .LE. 0 ) .OR. & 
           ( num_sm_levels_input .LE. 0 ) ) THEN
         PRINT '(A)','No input soil level data (either temperature or moisture, or both are missing).  Required for LSM.'
         STOP 'no soil data'
      ELSE
         ALLOCATE ( zhave( MAX(num_st_levels_input,num_sm_levels_input) +2) )
      END IF

      !  Sort the levels for temperature.

      outert : DO lout = 1 , num_st_levels_input-1
         innert : DO lin = lout+1 , num_st_levels_input
            IF ( st_levels_input(lout) .GT. st_levels_input(lin) ) THEN
               temp = st_levels_input(lout) 
               print *,'!!!',temp,lout,lin,num_st_levels_input
               st_levels_input(lout) = st_levels_input(lin)
               st_levels_input(lin) = NINT(temp)
               DO j = 1 , jte
                  DO i = 1 , ite
                     temp = st_input(i,j,lout+1)
                     st_input(i,j,lout+1) = st_input(i,j,lin+1)
                     st_input(i,j,lin+1) = temp
                  END DO
               END DO
            END IF
         END DO innert
      END DO outert
      DO j = 1 , jte
         DO i = 1 , ite
            st_input(i,j,1) = tsk(i,j)
            st_input(i,j,num_st_levels_input+2) = tmn(i,j)
         END DO
      END DO
      print *,'st-input',st_input

      !  Sort the levels for moisture.

      outerm: DO lout = 1 , num_sm_levels_input-1
         innerm : DO lin = lout+1 , num_sm_levels_input
            IF ( sm_levels_input(lout) .GT. sm_levels_input(lin) ) THEN
               temp = sm_levels_input(lout) 
               sm_levels_input(lout) = sm_levels_input(lin)
               sm_levels_input(lin) = NINT(temp)
               DO j = 1 , jte
                  DO i = 1 , ite
                     temp = sm_input(i,j,lout+1)
                     sm_input(i,j,lout+1) = sm_input(i,j,lin+1)
                     sm_input(i,j,lin+1) = temp
                  END DO
               END DO
            END IF
         END DO innerm
      END DO outerm
      DO j = 1 , jte
         DO i = 1 , ite
            sm_input(i,j,1) = sm_input(i,j,2)
            sm_input(i,j,num_sm_levels_input+2) = sm_input(i,j,num_sm_levels_input+1)
         END DO
      END DO
      print *,'sm-input',sm_input

      !  Here are the levels that we have from the input for temperature.  The input levels plus
      !  two more: the skin temperature at 0 cm, and the annual mean temperature at 300 cm.

      zhave(1) = 0.
      DO l = 1 , num_st_levels_input
         zhave(l+1) = st_levels_input(l) / 100.
      END DO
      zhave(num_st_levels_input+2) = 300. / 100.

      !  Interpolate between the layers we have (zhave) and those that we want (zs).

      z_wantt : DO lwant = 1 , num_soil_layers
         z_havet : DO lhave = 1 , num_st_levels_input +2 -1
            IF ( ( zs(lwant) .GE. zhave(lhave  ) ) .AND. &
                 ( zs(lwant) .LE. zhave(lhave+1) ) ) THEN
               DO j = 1 , jte
                  DO i = 1 , ite
                     tslb(i,lwant,j)= ( st_input(i,j,lhave  ) * ( zhave(lhave+1) - zs   (lwant) ) + &
                                        st_input(i,j,lhave+1) * ( zs   (lwant  ) - zhave(lhave) ) ) / &
                                                                ( zhave(lhave+1) - zhave(lhave) )
                  END DO
               END DO
               EXIT z_havet
            END IF
         END DO z_havet
      END DO z_wantt

      !  Here are the levels that we have from the input for moisture.  The input levels plus
      !  two more: a value at 0 cm and one at 300 cm.  The 0 cm value is taken to be identical
      !  to the most shallow layer's value.  Similarly, the 300 cm value is taken to be the same
      !  as the most deep layer's value.

      zhave(1) = 0.
      DO l = 1 , num_sm_levels_input
         zhave(l+1) = sm_levels_input(l) / 100.
      END DO
      zhave(num_sm_levels_input+2) = 300. / 100.

      !  Interpolate between the layers we have (zhave) and those that we want (zs).

      z_wantm : DO lwant = 1 , num_soil_layers
         z_havem : DO lhave = 1 , num_sm_levels_input +2 -1
            IF ( ( zs(lwant) .GE. zhave(lhave  ) ) .AND. &
                 ( zs(lwant) .LE. zhave(lhave+1) ) ) THEN
               DO j = 1 , jte
                  DO i = 1 , ite
                     smois(i,lwant,j)= ( sm_input(i,j,lhave  ) * ( zhave(lhave+1) - zs   (lwant) ) + &
                                         sm_input(i,j,lhave+1) * ( zs   (lwant  ) - zhave(lhave) ) ) / &
                                                                 ( zhave(lhave+1) - zhave(lhave) )
                  END DO
               END DO
               EXIT z_havem
            END IF
         END DO z_havem
      END DO z_wantm

      !  Over water, put in reasonable values for soil temperature and moisture.  These won't be
      !  used, but they will make a more continuous plot.

      IF ( flag_sst ) THEN
         DO j = 1 , jte
            DO i = 1 , ite
               IF ( landmask_input(i,j) .LT. 0.5 ) THEN
                  DO l = 1 , num_soil_layers
                     tslb(i,l,j)= sst_input(i,j)
                     smois(i,l,j)= 1.0
                  END DO
               END IF
            END DO
         END DO
      ELSE
         DO j = 1 , jte
            DO i = 1 , ite
               IF ( landmask_input(i,j) .LT. 0.5 ) THEN
                  DO l = 1 , num_soil_layers
                     tslb(i,l,j)= tsk(i,j)
                     smois(i,l,j)= 1.0
                  END DO
               END IF
            END DO
         END DO
      END IF

      DEALLOCATE (zhave)

   END SUBROUTINE init_soil_2_real

   SUBROUTINE init_soil_2_ideal ( xland,xice,vegfra,snow,canwat,     &
                     ivgtyp,isltyp,tslb,smois,tmn,                  &
                     num_soil_layers,                               &
                     ids,ide, jds,jde, kds,kde,                     &
                     ims,ime, jms,jme, kms,kme,                     &
                     its,ite, jts,jte, kts,kte                      )

      IMPLICIT NONE 
   
      INTEGER, INTENT(IN) ::ids,ide, jds,jde, kds,kde,  &
                            ims,ime, jms,jme, kms,kme,  &
                            its,ite, jts,jte, kts,kte
   
      INTEGER, INTENT(IN) ::num_soil_layers
   
      REAL, DIMENSION( ims:ime, num_soil_layers, jms:jme ) , INTENT(OUT) :: smois, tslb 
   
      REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(OUT)  :: xland, snow, canwat, xice, vegfra, tmn
   
      INTEGER, DIMENSION( ims:ime, jms:jme ) , INTENT(OUT) :: ivgtyp, isltyp
   
      INTEGER :: icm,jcm,itf,jtf
      INTEGER ::  i,j,l
   
      itf=min0(ite,ide)
      jtf=min0(jte,jde)
   
      icm = ide/2
      jcm = jde/2
   
      DO j=jts,jtf
         DO l=1,num_soil_layers
            DO i=its,itf
   
               smois(i,1,j)=0.10
               smois(i,2,j)=0.10
               smois(i,3,j)=0.10
               smois(i,4,j)=0.10
      
               tslb(i,1,j)=295.
               tslb(i,2,j)=297.
               tslb(i,3,j)=293.
               tslb(i,4,j)=293. 

            ENDDO
         ENDDO
      ENDDO                                 

      DO j=jts,jtf
         DO i=its,itf
            xland(i,j)  =   2
            tmn(i,j)    = 294. 
            xice(i,j)   =   0.
            vegfra(i,j) =   0. 
            snow(i,j)   =   0.
            canwat(i,j) =   0.
            ivgtyp(i,j) =   7
            isltyp(i,j) =   8
         ENDDO
      ENDDO

   END SUBROUTINE init_soil_2_ideal

   SUBROUTINE init_soil_3_real ( tsk , tmn , smois , tslb , &
                                 st_input , sm_input , landmask_input , sst_input , &
                                 zs , dzs , &
                                 st_levels_input , sm_levels_input , &
                                 num_soil_layers , num_st_levels_input , num_sm_levels_input ,  &
                                 flag_sst , &
                                 ids , ide , jds , jde , kds , kde , &
                                 ims , ime , jms , jme , kms , kme , &
                                 its , ite , jts , jte , kts , kte )

      IMPLICIT NONE

      INTEGER , INTENT(IN) :: num_soil_layers , num_st_levels_input , num_sm_levels_input , &
                              ids , ide , jds , jde , kds , kde , &
                              ims , ime , jms , jme , kms , kme , &
                              its , ite , jts , jte , kts , kte 

      LOGICAL , INTENT(IN) :: flag_sst

      INTEGER , DIMENSION(:) , INTENT(INOUT) :: st_levels_input , sm_levels_input

      REAL , DIMENSION(:,:,:) , INTENT(INOUT) :: st_input , sm_input
      REAL , DIMENSION(:,:) , INTENT(IN) :: landmask_input , sst_input 

      REAL , DIMENSION(ims:ime,jms:jme) , INTENT(IN) :: tsk , tmn
      REAL , DIMENSION(num_soil_layers) :: zs , dzs

      REAL , DIMENSION(ims:ime,num_soil_layers,jms:jme) , INTENT(OUT) :: tslb , smois

      REAL , ALLOCATABLE , DIMENSION(:) :: zhave

      INTEGER :: i , j , l , lout , lin , lwant , lhave
      REAL :: temp

      !  Allocate the soil layer array used for interpolating.

      IF ( ( num_st_levels_input .LE. 0 ) .OR. & 
           ( num_sm_levels_input .LE. 0 ) ) THEN
         PRINT '(A)','No input soil level data (either temperature or moisture, or both are missing).  Required for RUC LSM.'
         STOP 'no soil data'
      ELSE
         ALLOCATE ( zhave( MAX(num_st_levels_input,num_sm_levels_input) ) )
      END IF

      !  Sort the levels for temperature.

      outert : DO lout = 1 , num_st_levels_input-1
         innert : DO lin = lout+1 , num_st_levels_input
            IF ( st_levels_input(lout) .GT. st_levels_input(lin) ) THEN
               temp = st_levels_input(lout) 
               st_levels_input(lout) = st_levels_input(lin)
               st_levels_input(lin) = NINT(temp)
               DO j = 1 , jte
                  DO i = 1 , ite
                     temp = st_input(i,j,lout)
                     st_input(i,j,lout) = st_input(i,j,lin)
                     st_input(i,j,lin) = temp
                  END DO
               END DO
            END IF
         END DO innert
      END DO outert

      !  Sort the levels for moisture.

      outerm: DO lout = 1 , num_sm_levels_input-1
         innerm : DO lin = lout+1 , num_sm_levels_input
            IF ( sm_levels_input(lout) .GT. sm_levels_input(lin) ) THEN
               temp = sm_levels_input(lout) 
               sm_levels_input(lout) = sm_levels_input(lin)
               sm_levels_input(lin) = NINT(temp)
               DO j = 1 , jte
                  DO i = 1 , ite
                     temp = sm_input(i,j,lout)
                     sm_input(i,j,lout) = sm_input(i,j,lin)
                     sm_input(i,j,lin) = temp
                  END DO
               END DO
            END IF
         END DO innerm
      END DO outerm

      !  Here are the levels that we have from the input for temperature.

      DO l = 1 , num_st_levels_input
         zhave(l) = st_levels_input(l) / 100.
      END DO

      !  Interpolate between the layers we have (zhave) and those that we want (zs).

      z_wantt : DO lwant = 1 , num_soil_layers
         z_havet : DO lhave = 1 , num_st_levels_input -1
            IF ( ( zs(lwant) .GE. zhave(lhave  ) ) .AND. &
                 ( zs(lwant) .LE. zhave(lhave+1) ) ) THEN
               DO j = 1 , jte
                  DO i = 1 , ite
                     tslb(i,lwant,j)= ( st_input(i,j,lhave  ) * ( zhave(lhave+1) - zs   (lwant) ) + &
                                        st_input(i,j,lhave+1) * ( zs   (lwant  ) - zhave(lhave) ) ) / &
                                                                ( zhave(lhave+1) - zhave(lhave) )
                  END DO
               END DO
               EXIT z_havet
            END IF
         END DO z_havet
      END DO z_wantt

      !  Here are the levels that we have from the input for moisture.

      DO l = 1 , num_sm_levels_input
         zhave(l) = sm_levels_input(l) / 100.
      END DO

      !  Interpolate between the layers we have (zhave) and those that we want (zs).

      z_wantm : DO lwant = 1 , num_soil_layers
         z_havem : DO lhave = 1 , num_sm_levels_input -1
            IF ( ( zs(lwant) .GE. zhave(lhave  ) ) .AND. &
                 ( zs(lwant) .LE. zhave(lhave+1) ) ) THEN
               DO j = 1 , jte
                  DO i = 1 , ite
                     smois(i,lwant,j)= ( sm_input(i,j,lhave  ) * ( zhave(lhave+1) - zs   (lwant) ) + &
                                         sm_input(i,j,lhave+1) * ( zs   (lwant  ) - zhave(lhave) ) ) / &
                                                                 ( zhave(lhave+1) - zhave(lhave) )
                  END DO
               END DO
               EXIT z_havem
            END IF
         END DO z_havem
      END DO z_wantm

      !  Over water, put in reasonable values for soil temperature and moisture.  These won't be
      !  used, but they will make a more continuous plot.

      IF ( flag_sst ) THEN
         DO j = 1 , jte
            DO i = 1 , ite
               IF ( landmask_input(i,j) .LT. 0.5 ) THEN
                  DO l = 1 , num_soil_layers
                     tslb(i,l,j)= sst_input(i,j)
                     smois(i,l,j)= 1.0
                  END DO
               END IF
            END DO
         END DO
      ELSE
         DO j = 1 , jte
            DO i = 1 , ite
               IF ( landmask_input(i,j) .LT. 0.5 ) THEN
                  DO l = 1 , num_soil_layers
                     tslb(i,l,j)= tsk(i,j)
                     smois(i,l,j)= 1.0
                  END DO
               END IF
            END DO
         END DO
      END IF

      DEALLOCATE (zhave)

   END SUBROUTINE init_soil_3_real

   SUBROUTINE init_soil_old_real ( tsk , tmn , &
                                 smois , tslb , zs , dzs , num_soil_layers , &
                                 st000010_input , st010040_input , st040100_input , st100200_input , &
                                 st010200_input , &
                                 sm000010_input , sm010040_input , sm040100_input , sm100200_input , &
                                 sm010200_input , &
                                 landmask_input , sst_input , &
                                 ids , ide , jds , jde , kds , kde , &
                                 ims , ime , jms , jme , kms , kme , &
                                 its , ite , jts , jte , kts , kte )

      !  This is the old version of init_soil_temp_2.  Here we directly assign the
      !  soil t and moisture levels to WRF levels - no interpolation.

      IMPLICIT NONE

      INTEGER , INTENT(IN) :: num_soil_layers , &
                              ids , ide , jds , jde , kds , kde , &
                              ims , ime , jms , jme , kms , kme , &
                              its , ite , jts , jte , kts , kte 

      REAL , DIMENSION(ims:ime,jms:jme) , INTENT(IN) :: tsk , tmn
      REAL , DIMENSION(num_soil_layers) :: zs , dzs

      REAL , DIMENSION(:,:) , INTENT(IN) :: st000010_input , st010040_input , st040100_input , st100200_input , &
                                            st010200_input , &
                                            sm000010_input , sm010040_input , sm040100_input , sm100200_input , &
                                            sm010200_input , &
                                            landmask_input , sst_input

      REAL , DIMENSION(ims:ime,num_soil_layers,jms:jme) , INTENT(OUT) :: tslb , smois

      INTEGER :: i , j , l

      !  Soil temperature is linearly interpolated between the skin temperature (taken to be at a
      !  depth of 0 cm) and the various input temperature levels.

      DO j = 1 , jte
         DO i = 1 , ite
            IF ( landmask_input(i,j) .EQ. 1 ) THEN
               tslb(i,1,j)= st000010_input(i,j)
               tslb(i,2,j)= st010040_input(i,j)
               tslb(i,3,j)= st040100_input(i,j)
               tslb(i,4,j)= st100200_input(i,j)
!tslb(i,4,j)= st010200_input(i,j)
               smois(i,1,j)= sm000010_input(i,j)
               smois(i,2,j)= sm010040_input(i,j)
               smois(i,3,j)= sm040100_input(i,j)
               smois(i,4,j)= sm100200_input(i,j)
!smois(i,4,j)= sm010200_input(i,j)
            ELSE
               DO l = 1 , num_soil_layers
                  tslb(i,l,j)= sst_input(i,j)
                  smois(i,l,j)= 1.0
               END DO
            END IF
         END DO
      END DO

   END SUBROUTINE init_soil_old_real

   FUNCTION char2int1( string3 ) RESULT ( int1 )
      CHARACTER (LEN=3) , INTENT(IN) :: string3
      INTEGER :: i1 , int1
      READ(string3,fmt='(I3)') i1
      int1 = i1
   END FUNCTION char2int1

   FUNCTION char2int2( string6 ) RESULT ( int1 )
      CHARACTER (LEN=6) , INTENT(IN) :: string6
      INTEGER :: i2 , i1 , int1
      READ(string6,fmt='(I3,I3)') i1,i2
      int1 = ( i2 + i1 ) / 2
   END FUNCTION char2int2

END MODULE module_soil_pre
