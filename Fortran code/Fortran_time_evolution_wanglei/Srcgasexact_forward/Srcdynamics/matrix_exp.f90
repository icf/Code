SUBROUTINE matrix_exponential_r(A,Exp_A,tau,n)
!
!Purpose: If matrix_exponential is called with argument types A=real, Exp_A=real, tau=Real
!         Then compute Exp_A=EXP(-tau*A)
!
!Based on routines by schneider, b. i.(nsf)
!
use prec
IMPLICIT NONE
INTEGER                                :: n, i, j, k, info
REAL(KIND=rKind), DIMENSION(:,:)                 :: A
REAL(KIND=rKind), DIMENSION(:,:)                 :: Exp_A
REAL(KIND=rKind), DIMENSION(:),   ALLOCATABLE    :: Eigen_Values
REAL(KIND=rKind), DIMENSION(:,:), ALLOCATABLE    :: Eigen_Vectors
REAL(KIND=rKind), DIMENSION(:,:), ALLOCATABLE    :: Temp
REAL(KIND=rKind), DIMENSION(:),   ALLOCATABLE    :: Scratch
REAL(KIND=rKind)                                 :: tau
REAL(KIND=rKind):: expeig
integer:: statint
! Allocate storage for diagonalization routine.
ALLOCATE ( Eigen_Values(n), Eigen_Vectors(n,n), Scratch(10*n), Temp(n,n) , STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate matrix_exp variables'
			END IF 
Eigen_Vectors = A
!Call LAPACK routine to diagonalize double precision real symmetric matrix
CALL DSYEV('v','l',n,Eigen_Vectors,n,Eigen_Values,           &
              Scratch,10*n,info)
			  
! Form the matrix with exponentials of the eigenvalues on the diagonal
! Then similarity transform back into the original basis
DO i=1,n
  	DO j=1,n
	Exp_A(i,j)=0.0_rKind
		DO k=1,n
  	  	expeig=exp(-tau*Eigen_Values(k))
        Exp_A(i,j) = Exp_A(i,j) + Eigen_Vectors(i,k)*expeig*Eigen_Vectors(j,k)
		END DO
    END DO
END DO

! Deallocate the unneeded storage
DEALLOCATE ( Eigen_Values, Eigen_Vectors, Scratch, Temp  , STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate matrix_exp variables'
			END IF   
END SUBROUTINE matrix_exponential_r


SUBROUTINE matrix_exponential_c(A,Exp_A,t,n)
!
!Purpose: If matrix_exponential is called with argument types A=complex, Exp_A=complex, t=complex
!         Then compute Exp_A=EXP(-i*t*A)
!
!Based on routines by schneider, b. i.(nsf)
!
use prec
IMPLICIT NONE
INTEGER                                    :: n, i, j, k, info
COMPLEX(KIND=rKind), DIMENSION(:,:)                 :: A
COMPLEX(KIND=rKind), DIMENSION(:,:)                 :: Exp_A
COMPLEX(KIND=rKind), DIMENSION(:,:), ALLOCATABLE    :: Eigen_Vectors
REAL(KIND=rKind),     DIMENSION(:),   ALLOCATABLE    :: Eigen_Values
COMPLEX(KIND=rKind), DIMENSION(:),   ALLOCATABLE    :: Workv
COMPLEX(KIND=rKind), DIMENSION(:,:), ALLOCATABLE    :: Temp
REAL(KIND=rKind),     DIMENSION(:),   ALLOCATABLE    :: Rworkv
COMPLEX(KIND=rKind)                                 :: t
COMPLEX(KIND=rKind)                                 :: eye=(0.0_rKind,1.0_rKind)
COMPLEX(KIND=rKind) :: expeig
integer:: statint
! Allocate some storage for diagonalization routine.
ALLOCATE ( Eigen_Values(n), Eigen_Vectors(n,n), Workv(10*n), Rworkv(10*n), Temp(n,n)  , STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate matrix_exp variables'
			END IF   
print *, 'test', A
Eigen_Vectors = A
!Call LAPACK routine to diagonalize double precision hermitian matrix
CALL ZHEEV('v','l',n,Eigen_Vectors,n,Eigen_Values,              &
              Workv,10*n,Rworkv,info)

! Form the matrix with exponentials of the eigenvalues on the diagonal
! Then similarity transform back into the original basis
DO i=1,n
	DO j=1,n
	Exp_A(i,j)=CMPLX(0.0,KIND=rKind)
		DO k=1,n
  	  	expeig=exp(-eye*t*Eigen_Values(k))
        Exp_A(i,j) = Exp_A(i,j) + Eigen_Vectors(i,k)*expeig*CONJG(Eigen_Vectors(j,k))
		END DO
    END DO
END DO
  ! Deallocate the unneeded storage
DEALLOCATE ( Eigen_Values, Eigen_Vectors, Workv, Rworkv, Temp  , STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate matrix_exp variables'
			END IF     
END SUBROUTINE matrix_exponential_c


