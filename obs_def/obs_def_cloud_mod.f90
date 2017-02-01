! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

! BEGIN DART PREPROCESS KIND LIST
!CLOUD_LIQUID_WATER, KIND_CLOUD_LIQUID_WATER,  COMMON_CODE
!CLOUD_ICE,          KIND_CLOUD_ICE,           COMMON_CODE
!CLOUD_FRACTION,     KIND_CLOUD_FRACTION,      COMMON_CODE
! END DART PREPROCESS KIND LIST

! eventually these should become the following to be consistent
! with other types of mixing ratios.  this is just a proposed
! rename of the KIND; any code using it would remain the same.
! but multiple models and forward operators could be using the old
! names; the change must be coordinated across the entire project.
!!KIND_CLOUD_LIQUID_WATER -> KIND_CLOUDWATER_MIXING_RATIO
!!KIND_CLOUD_ICE -> KIND_ICE_MIXING_RATIO

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
