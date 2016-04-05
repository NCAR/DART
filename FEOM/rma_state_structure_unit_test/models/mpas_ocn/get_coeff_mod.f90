! This code may (or may not) be part of the MPAS distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id$

  module get_coeff_mod

  use types_mod, only : r8
  use utilities_mod, only : error_handler, E_ERR

  implicit none

  private
  save
  public :: mpas_rbf_interp_func_3D_plane_vec_const_dir_comp_coeffs

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

  character(len=256) :: string1

  contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! Purpose: Compute interpolation coefficients in 3D that can be used to
  !  interpolate a number of vector functions at a given locations. This is useful
  !  if the interpolation location does not change with time, or if several
  !  fields are to be interpolated at a given time step.  (If both the vector fields
  !  and the interpolation locations vary with time, there is no clear advantage in
  !  using either this method or the method for 2D interpoaltion above; for simplicity
  !  and because we foresee more uses for the method of this subroutine, we have not
  !  implemented a 3D version of the fixed field, variable interpolation location method
  !  as we have in 2D.) Coefficients are produced for handling Dirichlet boundary
  !  conditions (or no boundaries).  The interpolation is performed with basis functions
  !  that are RBFs plus a constant.
  ! Input:
  !  pointCount - the number of "source" points and functionValues supplied
  !  sourcePoints - the location of the "source" points in the 3D space where the values of
  !    the function are known.  The sourcePoints are projected into the plane given by
  !    planeBasisVectors
  !  unitVectors - the unit vectors associated with each of the sourcePoints.  Interpolation
  !    is performed by supplying the value of the vector function dotted into each of these unit
  !    vectors.  If multiple unit vectors are supplied at the same sourcePoint, they *must* be
  !    orthogonal for the interpolation to succeed.  The unitVectors are projected into the
  !    plane given by planeBasisVectors
  !  destinationPoint - the point where the interpolation will be performed.  The destinationPoint
  !    is projected into the plane given by planeBasisVectors
  !  alpha - a constant that give the characteristic length scale of the RBFs,
  !    should be on the order of the distance between points
  !  planeBasisVectors - the basis fectors for the plane where interpolation is performed.
  !    All points are projected into this plane. 
  ! Output:
  !  coefficients - the coefficients used to interpolate a function with Dirichlet
  !    boundary conditions to the specified destinationPoint
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  subroutine mpas_rbf_interp_func_3D_plane_vec_const_dir_comp_coeffs(pointCount, &
    sourcePoints, unitVectors, destinationPoint, &
    alpha, planeBasisVectors, coefficients)

    integer, intent(in) :: pointCount
    real(kind=r8), dimension(pointCount,3), intent(in) :: sourcePoints
    real(kind=r8), dimension(pointCount,3), intent(in) :: unitVectors
    real(kind=r8), dimension(3), intent(in) :: destinationPoint
    real(kind=r8), intent(in) :: alpha
    real(kind=r8), dimension(2,3) :: planeBasisVectors
    real(kind=r8), dimension(pointCount, 3), intent(out) :: coefficients

    integer :: i, j
    integer :: matrixSize

    real(kind=r8) :: rSquared, rbfValue, unitVectorDotProduct

    real(kind=r8), dimension(pointCount,2) :: planarSourcePoints
    real(kind=r8), dimension(pointCount,2) :: planarUnitVectors
    real(kind=r8), dimension(2) :: planarDestinationPoint

    real(kind=r8), dimension(:,:), pointer :: matrix, submatrix, matrixCopy
    real(kind=r8), dimension(:,:), pointer :: rhs, subrhs, coeffs
    integer, dimension(:), pointer :: pivotIndices

    matrixSize = pointCount+2 ! space for constant vector in plane

    allocate(matrix(matrixSize,matrixSize))  
    allocate(submatrix(pointCount,pointCount))
    allocate(matrixCopy(matrixSize,matrixSize))  
    allocate(rhs(matrixSize,2))
    allocate(subrhs(pointCount,2))
    allocate(coeffs(matrixSize,2))  
    allocate(pivotIndices(matrixSize)) 

    matrix = 0.0
    submatrix = 0.0
    rhs = 0.0
    subrhs = 0.0
    coeffs = 0.0

    do i = 1, pointCount
      planarSourcePoints(i,1) = sum(sourcePoints(i,:)*planeBasisVectors(1,:)) 
      planarSourcePoints(i,2) = sum(sourcePoints(i,:)*planeBasisVectors(2,:)) 
      planarUnitVectors(i,1) = sum(unitVectors(i,:)*planeBasisVectors(1,:)) 
      planarUnitVectors(i,2) = sum(unitVectors(i,:)*planeBasisVectors(2,:)) 
    end do
    planarDestinationPoint(1) = sum(destinationPoint*planeBasisVectors(1,:)) 
    planarDestinationPoint(2) = sum(destinationPoint*planeBasisVectors(2,:)) 

    call mpas_set_up_vector_dirichlet_rbf_matrix_and_rhs(pointCount, 2, &
      planarSourcePoints, planarUnitVectors, planarDestinationPoint, &
      alpha, submatrix, subrhs)

    matrix(1:pointCount,1:pointCount) = submatrix
    rhs(1:pointCount,:) = subrhs

    do i = 1, pointCount
      matrix(i,pointCount+1:pointCount+2) = planarUnitVectors(i,:) 
      matrix(pointCount+1:pointCount+2,i) = matrix(i,pointCount+1:pointCount+2)
    end do
    do i = 1,2 
      rhs(pointCount+i,i) = 1.0 ! the unit vector in the ith direction
    end do

    ! solve each linear system
    matrixCopy = matrix
    call mpas_legs(matrix, matrixSize, rhs(:,1), coeffs(:,1), pivotIndices)
    call mpas_legs(matrixCopy, matrixSize, rhs(:,2), coeffs(:,2), pivotIndices)

    do i = 1,3 
      coefficients(:,i) = planeBasisVectors(1,i)*coeffs(1:pointCount,1) &
        + planeBasisVectors(2,i)*coeffs(1:pointCount,2) 
    end do

    deallocate(matrix)
    deallocate(submatrix)
    deallocate(matrixCopy)  
    deallocate(rhs)
    deallocate(subrhs)
    deallocate(coeffs)  
    deallocate(pivotIndices)

  end subroutine mpas_rbf_interp_func_3D_plane_vec_const_dir_comp_coeffs 

 
  subroutine mpas_set_up_vector_dirichlet_rbf_matrix_and_rhs(pointCount, dimensions, &
    sourcePoints, unitVectors, destinationPoint, &
    alpha, matrix, rhs)

    integer, intent(in) :: pointCount, dimensions
    real(kind=r8), dimension(pointCount,dimensions), intent(in) :: sourcePoints
    real(kind=r8), dimension(pointCount,dimensions), intent(in) :: unitVectors
    real(kind=r8), dimension(dimensions), intent(in) :: destinationPoint
    real(kind=r8), intent(in) :: alpha
    real(kind=r8), dimension(pointCount,pointCount), intent(out) :: matrix
    real(kind=r8), dimension(pointCount,dimensions), intent(out) :: rhs

    integer :: i, j

    real(kind=r8) :: rSquared, rbfValue, unitVectorDotProduct

    do j = 1, pointCount
      do i = j, pointCount
        rSquared = sum((sourcePoints(i,:)-sourcePoints(j,:))**2)/alpha**2
        rbfValue = evaluate_rbf(rSquared)
        unitVectorDotProduct = sum(unitVectors(i,:)*unitVectors(j,:))
        matrix(i,j) = rbfValue*unitVectorDotProduct
        matrix(j,i) = matrix(i,j)
      end do
    end do

    do j = 1, pointCount
      rSquared = sum((destinationPoint-sourcePoints(j,:))**2)/alpha**2
      rhs(j,:) = evaluate_rbf(rSquared)*unitVectors(j,:)
    end do

  end subroutine mpas_set_up_vector_dirichlet_rbf_matrix_and_rhs


subroutine mpas_legs (A,N,B,X,INDX)
!
! subroutine to solve the equation A(N,N)*X(N) = B(N) with the
! partial-pivoting Gaussian elimination scheme.
! Copyright (c) Tao Pang 2001.
!
  IMPLICIT NONE
  integer, INTENT (IN) :: N
  integer :: I,J
  integer, INTENT (OUT), DIMENSION (N) :: INDX
  real(kind=r8), INTENT (INOUT), DIMENSION (N,N) :: A
  real(kind=r8), INTENT (INOUT), DIMENSION (N) :: B
  real(kind=r8), INTENT (OUT), DIMENSION (N) :: X
!
  CALL elgs (A,N,INDX)
!
  DO I = 1, N-1
    DO J = I+1, N
      B(INDX(J)) = B(INDX(J))-A(INDX(J),I)*B(INDX(I))
    END DO
  END DO
!
  X(N) = B(INDX(N))/A(INDX(N),N)
  DO I = N-1, 1, -1
    X(I) = B(INDX(I))
    DO J = I+1, N
      X(I) = X(I)-A(INDX(I),J)*X(J)
    END DO
    X(I) =  X(I)/A(INDX(I),I)
  END DO
!
END subroutine mpas_legs

function evaluate_rbf(rSquared) result(rbfValue)
  real(kind=r8), intent(in) :: rSquared
  real(kind=r8) :: rbfValue

  ! inverse multiquadratic
  rbfValue = 1/sqrt(1 + rSquared)

end function evaluate_rbf

subroutine elgs (A,N,INDX)
!
! subroutine to perform the partial-pivoting Gaussian elimination.
! A(N,N) is the original matrix in the input and transformed matrix
! plus the pivoting element ratios below the diagonal in the output.
! INDX(N) records the pivoting order.  Copyright (c) Tao Pang 2001.
!
  IMPLICIT NONE
  integer, INTENT (IN) :: N
  integer :: I,J,K,ITMP
  integer, INTENT (OUT), DIMENSION (N) :: INDX
  real(kind=r8) :: C1,PI,PI1,PJ
  real(kind=r8), INTENT (INOUT), DIMENSION (N,N) :: A
  real(kind=r8), DIMENSION (N) :: C
!
! Initialize the index
!
  DO I = 1, N
    INDX(I) = I
  END DO
!
! Find the rescaling factors, one from each row
!
  DO I = 1, N
    C1= 0.0
    DO J = 1, N
      !C1 = AMAX1(C1,ABS(A(I,J)))
      C1 = MAX(C1,ABS(A(I,J)))
    END DO
    C(I) = C1
  END DO
!
! Search the pivoting (largest) element from each column
!
  DO J = 1, N-1
    PI1 = 0.0
    K = 0
    DO I = J, N
      PI = ABS(A(INDX(I),J))/C(INDX(I))
      IF (PI.GT.PI1) THEN
        PI1 = PI
        K   = I
      ENDIF
    END DO
! This should never happen, but just in case:
    IF (K == 0) THEN
      write(string1,*)'K was never initialized!'
      call error_handler(E_ERR,'elgs',string1,source,revision,revdate)
    ENDIF

!
! Interchange the rows via INDX(N) to record pivoting order
!
    ITMP    = INDX(J)
    INDX(J) = INDX(K)
    INDX(K) = ITMP
    DO I = J+1, N
      PJ  = A(INDX(I),J)/A(INDX(J),J)
!
! Record pivoting ratios below the diagonal
!
      A(INDX(I),J) = PJ
!
! Modify other elements accordingly
!
      DO K = J+1, N
        A(INDX(I),K) = A(INDX(I),K)-PJ*A(INDX(J),K)
      END DO
    END DO
  END DO
!
END subroutine elgs

end module get_coeff_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
