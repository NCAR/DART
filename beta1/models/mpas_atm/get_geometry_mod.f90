! This code may (or may not) be part of the MPAS distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id$

module get_geometry_mod

   use types_mod, only : r8

   implicit none
   private
   save

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Purpose: perform interpolation of scalar and vector functions in 2D
!   and 3D using Radial Basis Functions (RBFs).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

   ! Initialize the geometry that will be useful from interpolation
  public :: get_geometry

  contains

  subroutine get_geometry(nData, xData, yData, zData, &
             xReconstruct, yReconstruct, zReconstruct, normalDirectionData, &
             on_a_sphere, dataTangentPlane)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! Purpose: compute geometric fields that will be potentially useful for calling
  !          the interpolation routines
  !
  ! Input: the grid
  !
  ! Output: 
  !  normalDirectionData - the unit vector at of each data point tangent to the sphere
  !  dataTangentPlane - 2 orthogonal unit vectors in the tangent plane of reconstruct
  !                     The first unit vector is chosen to point toward the first
  !                     data point.
  !  localVerticalUnitVectors - the unit normal vector of the tangent plane at reconsruct 
  !       
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

    implicit none

    integer :: nData
    real(kind=r8), dimension(:) :: xData, yData, zData
    real(kind=r8) :: xReconstruct, yReconstruct, zReconstruct
    real(kind=r8), dimension(:,:) :: normalDirectionData
    logical :: on_a_sphere

    real(kind=r8), dimension(3,2) :: dataTangentPlane

    integer :: iData
    real(kind=r8), dimension(3) :: localVerticalUnitVectors
    real(kind=r8), dimension(3) :: xHatPlane, yHatPlane, rHat
    real(kind=r8) :: normalDotRHat

    ! init arrays
    localVerticalUnitVectors = 0

    if(on_a_sphere) then
       localVerticalUnitVectors(1) = xReconstruct
       localVerticalUnitVectors(2) = yReconstruct
       localVerticalUnitVectors(3) = zReconstruct
       call mpas_unit_vec_in_r3(localVerticalUnitVectors)
    else ! on a plane
       localVerticalUnitVectors(:) = (/ 0.0_r8, 0.0_r8, 1.0_r8 /)
    end if

    iData = 1
    ! xHat and yHat are a local basis in the plane of the horizontal cell
    ! we arbitrarily choose xHat to point toward the first data point
    rHat = localVerticalUnitVectors
    normalDotRHat = sum(normalDirectionData(:,iData)*rHat)
    xHatPlane = normalDirectionData(:,iData) - normalDotRHat*rHat
    call mpas_unit_vec_in_r3(xHatPlane)

    call mpas_cross_product_in_r3(rHat, xHatPlane, yHatPlane)
    call mpas_unit_vec_in_r3(yHatPlane) ! just to be sure...
    dataTangentPlane(:,1) = xHatPlane
    dataTangentPlane(:,2) = yHatPlane
    !write(6,*)
    !write(6,*) ' dataTangentPlane '
    !write(6,10) dataTangentPlane(:,1)
    !write(6,10) dataTangentPlane(:,2)
    !write(6,*)
    10 format(3e20.10)

  end subroutine get_geometry

  subroutine mpas_unit_vec_in_r3(xin)
    implicit none
    real (kind=r8), intent(inout) :: xin(3)
    real (kind=r8) :: mag
    mag = sqrt(xin(1)**2+xin(2)**2+xin(3)**2)
    xin(:) = xin(:) / mag
  end subroutine mpas_unit_vec_in_r3

  subroutine mpas_cross_product_in_r3(p_1,p_2,p_out)
    real (kind=r8), intent(in) :: p_1 (3), p_2 (3)
    real (kind=r8), intent(out) :: p_out (3)

    p_out(1) = p_1(2)*p_2(3)-p_1(3)*p_2(2)
    p_out(2) = p_1(3)*p_2(1)-p_1(1)*p_2(3)
    p_out(3) = p_1(1)*p_2(2)-p_1(2)*p_2(1)
  end subroutine mpas_cross_product_in_r3

end module get_geometry_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
