! This code may (or may not) be part of the MPAS distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id$

module get_reconstruct_mod

  use types_mod, only : r8
  use get_coeff_mod

  implicit none

  public :: get_reconstruct_init, get_reconstruct

  contains

  subroutine get_reconstruct_init(nData, xData, yData, zData, &
              xReconstruct, yReconstruct, zReconstruct, normalDirectionData, &
              dataTangentPlane, coeffs_reconstruct)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! Purpose: pre-compute coefficients used by the reconstruct() routine
  !
  ! Input: nData, {x,y,z}Data, {x,y,z}Reconstruct, normalDirectionData,
  !        dataTangentPlane
  !
  ! Output: coeffs_reconstruct - coefficients used to reconstruct 
  !                              velocity vectors at reconstruct 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

    implicit none

    integer, intent(in) :: nData
    real(kind=r8), intent(in)    :: xReconstruct, yReconstruct, zReconstruct
    real(kind=r8), dimension(nData), intent(in)    :: xData, yData, zData
    real(kind=r8), dimension(3,nData), intent(in)  :: normalDirectionData
    real(kind=r8), dimension(3,2), intent(in)      :: dataTangentPlane
    real(kind=r8), dimension(3,nData), intent(out) :: coeffs_reconstruct

    ! local vars
    real(kind=r8) :: r, cellCenter(3), alpha, tangentPlane(2,3), Reconstruct(3)
    real(kind=r8), allocatable, dimension(:,:) :: dataLocations, &
      coeffs, normals
    integer :: iData, nData8


    ! init arrays
    coeffs_reconstruct = 0.0_r8

    allocate(dataLocations(nData,3))
    allocate(coeffs(nData,3))
    allocate(normals(nData,3))

    ! put reconstruct location into a vector
    Reconstruct(1) = xReconstruct
    Reconstruct(2) = yReconstruct
    Reconstruct(3) = zReconstruct

    do iData=1,nData
      dataLocations(iData, 1)  = xData(iData)
      dataLocations(iData, 2)  = yData(iData)
      dataLocations(iData, 3)  = zData(iData)
    end do

    alpha = 0.0_r8
    do iData=1,nData
      r = sqrt(sum((Reconstruct - dataLocations(iData,:))**2))
      alpha = alpha + r
    enddo
    alpha = alpha/nData

    normals(:,1) = normalDirectionData(1,:)
    normals(:,2) = normalDirectionData(2,:)
    normals(:,3) = normalDirectionData(3,:)

    tangentPlane(1,:) = dataTangentPlane(:,1)
    tangentPlane(2,:) = dataTangentPlane(:,2)

    ! the main call...
    nData8=nData
    call mpas_rbf_interp_func_3D_plane_vec_const_dir_comp_coeffs(nData8, &
        dataLocations, normals, &
        Reconstruct, alpha, tangentPlane, coeffs)
      
    do iData=1,nData
      coeffs_reconstruct(:,iData) = coeffs(iData,:)
    end do

    deallocate(dataLocations)
    deallocate(coeffs)
    deallocate(normals)

  end subroutine get_reconstruct_init



  subroutine get_reconstruct(nData, latReconstruct, lonReconstruct, &
                 coeffs_reconstruct, on_a_sphere, &
                 u, uReconstructX, uReconstructY, uReconstructZ, &
                 uReconstructZonal, uReconstructMeridional)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Purpose: reconstruct vector field at reconstruct point based on radial basis functions
  !
  ! Input: grid meta data and vector component data residing at reconstruct point
  !
  ! Output: reconstructed vector field (measured in X,Y,Z) located at the reconstruct point
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    integer,       intent(in)                     :: nData
    real(kind=r8), intent(in)                     :: latReconstruct, lonReconstruct
    real(kind=r8), dimension(3,nData), intent(in) :: coeffs_reconstruct
    logical, intent(in)                           :: on_a_sphere
    real(kind=r8), dimension(:), intent(in)       :: u
    real(kind=r8), intent(out)                    :: uReconstructX, uReconstructY, uReconstructZ
    real(kind=r8), intent(out)                    :: uReconstructZonal, uReconstructMeridional

    ! local vars
    integer       :: iData
    real(kind=r8) :: clat, slat, clon, slon

    ! init the intent(out)
    uReconstructX = 0.0_r8
    uReconstructY = 0.0_r8
    uReconstructZ = 0.0_r8

    ! a more efficient reconstruction where rbf_values*matrix_reconstruct has been precomputed
    ! in coeffs_reconstruct
    do iData=1,nData
      uReconstructX = uReconstructX &
        + coeffs_reconstruct(1,iData) * u(iData)
      uReconstructY = uReconstructY &
        + coeffs_reconstruct(2,iData) * u(iData)
      uReconstructZ = uReconstructZ &
        + coeffs_reconstruct(3,iData) * u(iData)
    enddo

    if(on_a_sphere) then
      clat = cos(latReconstruct)
      slat = sin(latReconstruct)
      clon = cos(lonReconstruct)
      slon = sin(lonReconstruct)
      uReconstructZonal = -uReconstructX*slon + uReconstructY*clon
      uReconstructMeridional = -(uReconstructX*clon &
        + uReconstructY*slon)*slat &
        + uReconstructZ*clat
   else
      uReconstructZonal = uReconstructX
      uReconstructMeridional = uReconstructY
    end if

  end subroutine get_reconstruct

end module get_reconstruct_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
