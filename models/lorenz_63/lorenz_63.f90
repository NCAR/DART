module model_mod

private
public model_size, init_conditions, adv_1step, advance

integer, parameter :: model_size = 3

!  define model parameters
double precision, PARAMETER :: sigma = 10., r = 28., b = 8./3.
double precision, parameter :: deltat = 0.01

contains

!
!==================================================================
!
subroutine comp_dt(x, dt)

!  COMPUTES TIME TENDENCY OF THE LORENZ 1963 3-VARIABLE MODEL GIVEN 
!  CURRENT STATE

implicit none

double precision, intent(in) :: x(:)
double precision, intent(out) :: dt(:)

!  COMPUTE THE LORENZ MODEL DT FROM STANDARD EQUATIONS
dt(1) = sigma * (x(2) - x(1))
dt(2) = -1.0*x(1)*x(3) + r*x(1) - x(2)
dt(3) = x(1)*x(2) - b*x(3)

return
end subroutine comp_dt
!
!===================================================================
!

!  ADVANCES THE 3 VARIABLE LORENZ-63 MODEL BY A GIVEN NUMBER OF STEPS
!  CURRENT STATE IN X, NEW STATE IN XNEW, NUM TIME STEPS ADVANCED

subroutine advance(x, num, xnew)

   IMPLICIT NONE

   double precision, intent(in) :: x(:)
   double precision, intent(out) :: xnew(:)
   integer, intent(in) :: num

   integer i

!  COPY INITIAL CONDITIONS TO AVOID OVERWRITE
xnew = x
   
!  ADVANCE THE APPROPRIATE NUMBER OF STEPS
DO I = 1, NUM
   CALL ADV_1STEP(XNEW)
END DO

return
end subroutine advance
!
!===================================================================
!
   SUBROUTINE INIT_CONDITIONS(X)

!  OFF-ATTRACTOR INITIAL CONDITIONS FOR LORENZ 63

   IMPLICIT NONE

   DOUBLE PRECISION, intent(out) :: X(:)
   INTEGER I

   X = 0.10

   RETURN
   end subroutine init_conditions
!
!====================================================================
!
   SUBROUTINE LINEAR_DT(X, DX, DT)

!  OLD VERSION OF LINEARIZED LORENZ 63 MODEL TIME TENDENCY COMPUTATION
   
   IMPLICIT NONE
   
   DOUBLE PRECISION, intent(in) :: X(:), DX(:)
   double precision, intent(out) :: dt(:)
!  COMPUTE LINEAR MODEL LORENZ TIME TENDENCY
   DT(1) = -1.0*SIGMA*DX(1) + SIGMA*DX(2)
   DT(2) = (R - X(3))*DX(1) - DX(2) - X(1)*DX(3)
   DT(3) = X(2)*DX(1) + X(1)*DX(2) - B*DX(3)

   RETURN
   END subroutine linear_dt
!
!====================================================================
!
!  DOES SINGLE TIME STEP ADVANCE FOR LORENZ CONVECTIVE 3 VARIABLE MODEL
!  USING TWO STEP RK TIME STEP

   SUBROUTINE ADV_1STEP(X)

   IMPLICIT NONE

   DOUBLE PRECISION, intent(inout) :: X(:)
   double precision :: fract = 1.0

   FRACT = 1.0
   CALL ADV_SINGLE(X, FRACT)

   RETURN
   END  subroutine adv_1step

!
!====================================================================
!
!  DOES SINGLE TIME STEP ADVANCE FOR LORENZ CONVECTIVE 3 VARIABLE MODEL
!  USING TWO STEP RK TIME STEP

   SUBROUTINE ADV_SINGLE(X, FRACT)

   IMPLICIT NONE

DOUBLE PRECISION, intent(inout) :: X(:)
double precision, intent(in) :: fract
double precision :: x1(3), x2(3), dx(3)

INTEGER I

!  COMPUTE THE FIRST INTERMEDIATE STEP
   CALL COMP_DT(X, DX)
   X1 = X + FRACT * DELTAT * DX
!  COMPUTE THE SECOND INTERMEDIATE STEP
   CALL COMP_DT(X1, DX)
   X2 = X1 + FRACT * DELTAT * DX

!  NEW VALUE FOR X IS AVERAGE OF ORIGINAL VALUE AND SECOND INTERMEDIATE
   X = (X + X2) / 2.0

   RETURN
   END subroutine adv_single

!
!=====================================================================
!
!  DOES SINGLE TIME STEP ADVANCE FOR LORENZ CONVECTIVE 3 VARIABLE MODEL
!  USING FOUR STEP RK TIME STEP

   SUBROUTINE ADV_SINGLE_RK4(X, FRACT)

   IMPLICIT NONE

   DOUBLE PRECISION, intent(inout) :: X(:)
double precision, intent(in) :: FRACT
double precision :: X1(3), X2(3), X3(3), X4(3), DX(3), inter(3)
   INTEGER I

!  COMPUTE THE FIRST INTERMEDIATE STEP
   CALL COMP_DT(X, DX)
      X1 = FRACT * DELTAT * DX
      INTER = X + X1 / 2.0
!  COMPUTE THE SECOND INTERMEDIATE STEP
   CALL COMP_DT(INTER, DX)
      X2 = FRACT * DELTAT * DX
      INTER = X + X2 / 2.0
!  COMPUTE THE THIRD INTERMEDIATE STEP
   CALL COMP_DT(INTER, DX)
      X3 = FRACT * DELTAT * DX
      INTER = X + X3
!  COMPUTE FOURTH INTERMEDIATE STEP
   CALL COMP_DT(INTER, DX)
      X4 = FRACT * DELTAT * DX

!  COMPUTE NEW VALUE FOR X
      X = X + X1/6.0 + X2/3.0 + X3/3.0 + X4/6.0

   RETURN
   END subroutine adv_single_rk4

!
!=====================================================================
!
   SUBROUTINE INV_LINEAR_DT(X, DX, PX)

   IMPLICIT NONE

   DOUBLE PRECISION, intent(in) :: X(:), DX(:)
   double precision, intent(out) :: PX(3)
   DOUBLE PRECISION A(3, 3), FACT, TDX(3)
   INTEGER I

!  COMPUTE INV LINEAR MODEL LORENZ TIME TENDENCY (SEE NOTES 13MAR94)
!  FOR NOW ASSUMES STUPID LEAP FROG, WILL THIS BE SUFFICIENT?
   A(1, 1) = -SIGMA * DELTAT + 1.0
   A(1, 2) = SIGMA * DELTAT
   A(1, 3) = 0.0
   A(2, 1) = (R - X(3)) * DELTAT
   A(2, 2) = -1.0 * DELTAT + 1.0
   A(2, 3) = -X(1) * DELTAT
   A(3, 1) = X(2) * DELTAT
   A(3, 2) = X(1) * DELTAT
   A(3, 3) = -B * DELTAT + 1.0
   
!  INITIALIZE COPY OF DX
   write(*, *) 'this routine is not up to date'
   if(1 == 1) stop
!      TDX(I) = DX(I)
      TDX = DX

!  GET RID OF A(2, 3)
   FACT = A(2, 3) / A(3, 3)
      A(2, :) = A(2, :) - FACT * A(3, :)
   TDX(2) = TDX(2) - FACT * TDX(3)

!  GET RID OF A(1, 2)
   FACT = A(1, 2) / A(2, 2)
      A(1, :) = A(1, :) - FACT * A(2, :)
   TDX(1) = TDX(1) - FACT * TDX(2)

!  SOLVE FOR THE PREVIOUS STEP LINEAR PERTURBATION
   PX(1) = TDX(1) / A(1, 1)
   PX(2) = (TDX(2) - A(2, 1) * PX(1)) / A(2, 2)
   PX(3) = (TDX(3) - A(3, 1) * PX(1) - A(3, 2) * PX(2)) / A(3, 3)

   END subroutine inv_linear_dt

!========================================================================
!
        SUBROUTINE LINEARIZE(NL, L)
        IMPLICIT NONE

!  COMPUTE LINEAR OPERATOR AROUND STATE NL
        DOUBLE PRECISION NL(3), L(3, 3)

        L(1, 1) = -1.0 * SIGMA * DELTAT + 1.0
        L(1, 2) = SIGMA * DELTAT
        L(1, 3) = 0.0 * DELTAT
        L(2, 1) = (R - NL(3)) * DELTAT
        L(2, 2) = -1.0 * DELTAT + 1.0
        L(2, 3) = -1.0 * NL(1) * DELTAT
        L(3, 1) = NL(2) * DELTAT
        L(3, 2) = NL(1) * DELTAT
        L(3, 3) = -1.0 * B * DELTAT + 1.0
        RETURN
        END subroutine linearize
   
!=========================================================================

end module model_mod
