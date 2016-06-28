! DART software - Copyright 2004 - 2011 UCAR. This open source software
! is
! provided by UCAR, "as is", without charge, subject to all terms of use
! at
! http://www.image.ucar.edu/DAReS/DART/DART_download

! FIXME: check to see if obs are of volume or thickness - for now we
! will assume volume.

! FIXME: do we want to identify the satellite? (yes)
!  AMSRE is a passive microwave

! BEGIN DART PREPROCESS KIND LIST
!SAT_SEAICE_AGREG_CONCENTR,       KIND_SEAICE_AGREG_CONCENTR,     COMMON_CODE
!SAT_SEAICE_AGREG_VOLUME,         KIND_SEAICE_AGREG_VOLUME,       COMMON_CODE
!SAT_SEAICE_AGREG_SNOWVOLUME,     KIND_SEAICE_AGREG_SNOWVOLUME,   COMMON_CODE
!SAT_SEAICE_AGREG_THICKNESS,      KIND_SEAICE_AGREG_THICKNESS,    COMMON_CODE
!SAT_SEAICE_AGREG_SNOWDEPTH,      KIND_SEAICE_AGREG_SNOWDEPTH,    COMMON_CODE
!SAT_U_SEAICE_COMPONENT,          KIND_U_SEAICE_COMPONENT,        COMMON_CODE
!SAT_V_SEAICE_COMPONENT,          KIND_V_SEAICE_COMPONENT,        COMMON_CODE
!SAT_SEAICE_CONCENTR,             KIND_SEAICE_CONCENTR,           COMMON_CODE
!SAT_SEAICE_VOLUME,               KIND_SEAICE_VOLUME,             COMMON_CODE
!SAT_SEAICE_SNOWVOLUME,           KIND_SEAICE_SNOWVOLUME,         COMMON_CODE
! END DART PREPROCESS KIND LIST

