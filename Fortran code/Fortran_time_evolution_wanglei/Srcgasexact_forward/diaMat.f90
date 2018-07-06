

!====================================================================
!    PROGRAM: diaMat
!    TYPE   : 
!    PURPOSE: This subroutine to diagolize the Nambu Mat, keep the eigenstate and eigenvalue
!    I/O    :
!    VERSION: 24-Aug-07
!    COMMENT:  (1)  
!              (2)  
!              (3) 
!                  
!======================================================================
subroutine diaMat_realsymm(matsize, mat, energy, eigenstate)   ! To generate the kin term in the matrix
use prec         
implicit none 

integer,intent(IN):: matsize
real(kind=RKind),intent(IN):: mat(matsize, matsize)
real(kind=RKind),intent(OUT):: energy(matsize), eigenstate(matsize, matsize)


	!-----------for RS routine-
     integer:: ierr
	 real (kind=Rkind),dimension(matsize):: fv1, fv2 
	!------------


    !CALL dEVcsf(Matsize, mat, Matsize , energy,eigenstate, Matsize)       ! real symm mattrix,abs(eigenvalue) from large to small 


	call RS(Matsize,Matsize,mat,energy,1,eigenstate, fv1,fv2,ierr)   ! real symm matrix, eigenvaule from small to large

    ! if (ierr .ne. 0) then
     !  print *, 'There is problem in diagonalizing Hamilt &
     !           & using RS() in subroutine hamdiag().'
     !  stop
     !endif
  

end subroutine diaMat_realsymm
!====================================================================


!====================================================================
!    PROGRAM: diaMat
!    TYPE   : 
!    PURPOSE: This subroutine to diagolize the Nambu Mat, keep the eigenstate and eigenvalue
!    I/O    :
!    VERSION: 24-Aug-07
!    COMMENT:  (1)  
!              (2)  
!              (3) 
!                  
!======================================================================
subroutine diaMat_tri(matsize,a,b, energy, eigenstate)   ! To generate the kin term in the matrix
use prec         
implicit none 

integer,intent(IN):: matsize
real(kind=RKind),intent(IN):: a(matsize), b(matsize)
real(kind=RKind),intent(OUT):: energy(matsize), eigenstate(matsize, matsize)

     integer:: ierr, ii
    
	eigenstate=0.0_RKind
	do ii=1, matsize
	 eigenstate(ii, ii)=1.0_RKind
	enddo 
      
	call TQL2(Matsize,Matsize,a,b,eigenstate,ierr)   ! real tri matrix, eigenvaule from small to large

     if (ierr .ne. 0) then
       print *, 'There is problem in diagonalizing Hamilt &
                & using TQL2()'
       stop
     endif
  
  energy=a

end subroutine diaMat_tri
!====================================================================


!C========+=========+=========+=========+=========+=========+=========+=$
!C PROGRAM: a few subroutines from slatec
!C          which were minimally modified to compute 
!C          double precision eigensystems.
!C TYPE   : main
!C PURPOSE: 
!C I/O    :
!C VERSION: 
!C COMMENT: to learn about netlib, send an otherwise empty
!C          message to netlib@research.att.com
!C          containing 'send index' in the subject header)
!C          The WWW address of netlib is
!C          http://netlib.att.com/netlib/search.html
!Cnoprint=+=========+=========+=========+=========+=========+=========+=$
!*DECK RS
      SUBROUTINE RS (NM, N, A, W, MATZ, Z, FV1, FV2, IERR)
!C***BEGIN PROLOGUE  RS
!C***PURPOSE  Compute the eigenvalues and, optionally, the eigenvectors
!C            of a real symmetric matrix.
!C***LIBRARY   SLATEC (EISPACK)
!C***CATEGORY  D4A1
!C***TYPE      SINGLE PRECISION (RS-S, CH-C)
!C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!C***AUTHOR  Smith, B. T., et al.
!C***DESCRIPTION
!C
!C     This subroutine calls the recommended sequence of
!C     subroutines from the eigensystem subroutine package (EISPACK)
!C     to find the eigenvalues and eigenvectors (if desired)
!C     of a DOUBLE PRECISION SYMMETRIC matrix.
!C
!C     On Input
!C
!C        NM must be set to the row dimension of the two-dimensional
!C          array parameters, A and Z, as declared in the calling
!C          program dimension statement.  NM is an INTEGER variable.
!C
!C        N is the order of the matrix A.  N is an INTEGER variable.
!C          N must be less than or equal to NM.
!C
!C        A contains the real symmetric matrix.  A is a two-dimensional
!C          DOUBLE PRECISION array, dimensioned A(NM,N).
!C
!C        MATZ is an INTEGER variable set equal to zero if only
!C          eigenvalues are desired.  Otherwise, it is set to any
!C          non-zero integer for both eigenvalues and eigenvectors.
!C
!C     On Output
!C
!C        A is unaltered.
!C
!C        W contains the eigenvalues in ascending order.  W is a one-
!C          dimensional DOUBLE PRECISION array, dimensioned W(N).
!C
!C        Z contains the eigenvectors if MATZ is not zero.  The
!C          eigenvectors are orthonormal.  Z is a two-dimensional
!C          DOUBLE PRECISION array, dimensioned Z(NM,N).
!C
!C        IERR is an INTEGER flag set to
!C          Zero       for normal return,
!C          10*N       if N is greater than NM,
!C          J          if the J-th eigenvalue has not been
!C                     determined after 30 iterations.
!C                     The eigenvalues, and eigenvectors if requested,
!C                     should be correct for indices 1, 2, ..., IERR-1.
!C
!C  FV1 and FV2 are one-dimensional DOUBLE PRECISION arrays used for temporary
!C          storage, dimensioned FV1(N) and FV2(N).
!C
!C     Questions and comments should be directed to B. S. Garbow,
!C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!C     ------------------------------------------------------------------
!C
!C***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!C                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!C                 system Routines - EISPACK Guide, Springer-Verlag,
!C                 1976.
!C***ROUTINES CALLED  TQL2, TQLRAT,  TRED2
!C***REVISION HISTORY  (YYMMDD)
!C   760101  DATE WRITTEN
!C   890831  Modified array declarations.  (WRB)
!C   890831  REVISION DATE from Version 3.2
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   920501  Reformatted the REFERENCES section.  (WRB)
!C***END PROLOGUE  RS
!C
!===== test to change it into Fortran90 format.=========
!      INTEGER N,NM,IERR,MATZ
!      DOUBLE PRECISION A(NM,*),W(*),Z(NM,*),FV1(*),FV2(*)      
!===========
    use prec  
    implicit none
    integer,intent(IN):: N, NM, MATZ
    integer,intent(OUT):: IERR
    real(kind=rkind),dimension(NM,N),intent(IN):: A
    real(kind=rkind),dimension(NM,N),intent(OUT)::Z
    real(kind=rkind),dimension(N),intent(OUT)::W
    real(kind=rkind),dimension(N),intent(IN)::FV1,FV2
!========================================================
!   print *, (A(3,i),i=1,8)
!   pause
!C
!C***FIRST EXECUTABLE STATEMENT  RS
      IERR = 10 * N
!C     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
      CALL  TRED2(NM,N,A,W,FV1,Z)
      CALL  TQL2(NM,N,W,FV1,Z,IERR)
!      END
      end subroutine RS




!*DECK TQL2
!
      SUBROUTINE TQL2 (NM, N, D, E, Z, IERR)
!C***BEGIN PROLOGUE  TQL2
!C***PURPOSE  Compute the eigenvalues and eigenvectors of symmetric
!C            tridiagonal matrix.
!C***LIBRARY   SLATEC (EISPACK)
!C***CATEGORY  D4A5, D4C2A
!C***TYPE      SINGLE PRECISION (TQL2-S)
!C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!C***AUTHOR  Smith, B. T., et al.
!C***DESCRIPTION
!C
!C     This subroutine is a translation of the ALGOL procedure TQL2,
!C     NUM. MATH. 11, 293-306(1968) by Bowdler, Martin, Reinsch, and
!C     Wilkinson.
!C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
!C
!C     This subroutine finds the eigenvalues and eigenvectors
!C     of a SYMMETRIC TRIDIAGONAL matrix by the QL method.
!C     The eigenvectors of a FULL SYMMETRIC matrix can also
!C     be found if  TRED2  has been used to reduce this
!C     full matrix to tridiagonal form.
!C
!C     On Input
!C
!C        NM must be set to the row dimension of the two-dimensional
!C          array parameter, Z, as declared in the calling program
!C          dimension statement.  NM is an INTEGER variable.
!C
!C        N is the order of the matrix.  N is an INTEGER variable.
!C          N must be less than or equal to NM.
! C
!C        D contains the diagonal elements of the symmetric tridiagonal
!!  matrix.  D is a one-dimensional DOUBLE PRECISION array, dimensioned D(N).
!C
!C        E contains the subdiagonal elements of the symmetric
!C          tridiagonal matrix in its last N-1 positions.  E(1) is
!C          arbitrary.  E is a one-dimensional DOUBLE PRECISION array,
!C           dimensioned
!C          E(N).
!C
!C        Z contains the transformation matrix produced in the
!C          reduction by  TRED2, if performed.  If the eigenvectors
!C          of the tridiagonal matrix are desired, Z must contain
!C          the identity matrix.  Z is a two-dimensional DOUBLE PRECISION
!C          array dimensioned Z(NM,N).
!C
!C      On Output
!C
!C        D contains the eigenvalues in ascending order.  If an
!C          error exit is made, the eigenvalues are correct but
!C          unordered for indices 1, 2, ..., IERR-1.
!C
!C        E has been destroyed.
!C
!C        Z contains orthonormal eigenvectors of the symmetric
!C          tridiagonal (or full) matrix.  If an error exit is made,
!C          Z contains the eigenvectors associated with the stored
!C          eigenvalues.
!C
!C        IERR is an INTEGER flag set to
!C          Zero       for normal return,
!C          J          if the J-th eigenvalue has not been
!C                     determined after 30 iterations.
!C
!C     Calls PYTHAG(A,B) for sqrt(A**2 + B**2).
!C
!C     Questions and comments should be directed to B. S. Garbow,
!C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!C     ------------------------------------------------------------------
!C
!C***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!C                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!C                 system Routines - EISPACK Guide, Springer-Verlag,
!C                 1976.
!C***ROUTINES CALLED  PYTHAG
!C***REVISION HISTORY  (YYMMDD)
!C   760101  DATE WRITTEN
!C   890831  Modified array declarations.  (WRB)
!C   890831  REVISION DATE from Version 3.2
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   920501  Reformatted the REFERENCES section.  (WRB)
!C***END PROLOGUE  TQL2
!C
!====== test to change it to Fortran 90 format.==========
!      INTEGER I,J,K,L,M,N,II,L1,L2,NM,MML,IERR
!      DOUBLE PRECISION D(*),E(*),Z(NM,*)
!    DOUBLE PRECISION B,C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2
!      DOUBLE PRECISION PYTHAG
!========
      use prec
      implicit none
      integer,intent(IN)::NM,N
      integer,intent(INOUT)::IERR
      integer:: I,J,K,L,M,II,L1,L2,MML
      real(kind=rkind),dimension(N),intent(INOUT):: D,E
      real(kind=rkind),dimension(NM,N),intent(INOUT):: Z
      real(kind=rkind),external::PYTHAG
      real(kind=rkind):: B,C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2
!===========================================================
!C
!C***FIRST EXECUTABLE STATEMENT  TQL2
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
!C
      DO 100 I = 2, N
  100 E(I-1) = E(I)
!C
      F = 0.0_rkind
      B = 0.0_rkind
      E(N) = 0.0_rkind
!C
      DO 240 L = 1, N
         J = 0
         H = ABS(D(L)) + ABS(E(L))
         IF (B .LT. H) B = H
!C     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
         DO 110 M = L, N
            IF (B + ABS(E(M)) .EQ. B) GO TO 120
!C     .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
!C                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
!C
  120    IF (M .EQ. L) GO TO 220
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
!C     .......... FORM SHIFT ..........
         L1 = L + 1
         L2 = L1 + 1
         G = D(L)
         P = (D(L1) - G) / (2.0_rkind * E(L))
         R = PYTHAG(P,1.0_rkind)
         D(L) = E(L) / (P + SIGN(R,P))
         D(L1) = E(L) * (P + SIGN(R,P))
         DL1 = D(L1)
         H = G - D(L)
         IF (L2 .GT. N) GO TO 145
!C
         DO 140 I = L2, N
  140    D(I) = D(I) - H
!C
  145    F = F + H
!C     .......... QL TRANSFORMATION ..........
         P = D(M)
         C = 1.0_rkind
         C2 = C
         EL1 = E(L1)
         S = 0.0_rkind
         MML = M - L
!C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            C3 = C2
            C2 = C
            S2 = S
            I = M - II
            G = C * E(I)
            H = C * P
            IF (ABS(P) .LT. ABS(E(I))) GO TO 150
            C = E(I) / P
            R = SQRT(C*C+1.0_rkind)
            E(I+1) = S * P * R
            S = C / R
            C = 1.0_rkind/ R
            GO TO 160
  150       C = P / E(I)
            R = SQRT(C*C+1.0_rkind)
            E(I+1) = S * E(I) * R
            S = 1.0_rkind / R
            C = C * S
  160       P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
!C     .......... FORM VECTOR ..........
            DO 180 K = 1, N
               H = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * H
               Z(K,I) = C * Z(K,I) - S * H
  180       CONTINUE
!C
  200    CONTINUE
!C
         P = -S * S2 * C3 * EL1 * E(L) / DL1
         E(L) = S * P
         D(L) = C * P
         IF (B + ABS(E(L)) .GT. B) GO TO 130
  220    D(L) = D(L) + F
  240 CONTINUE
!C     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
!C
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
!C
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
!C
         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  280    CONTINUE
!C
  300 CONTINUE
!C
      GO TO 1001
!C     .......... SET ERROR -- NO CONVERGENCE TO AN
!C                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
!
! 1001 RETURN
!      END
 1001  end subroutine TQL2




!*DECK TRED2
      SUBROUTINE TRED2 (NM, N, A, D, E, Z)
!C***BEGIN PROLOGUE  TRED2
!C***PURPOSE  Reduce a real symmetric matrix to a symmetric tridiagonal
!C            matrix using and accumulating orthogonal transformations.
!C***LIBRARY   SLATEC (EISPACK)
!C***CATEGORY  D4C1B1
!C***TYPE      SINGLE PRECISION (TRED2-S)
!C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!C***AUTHOR  Smith, B. T., et al.
!C***DESCRIPTION
!C
!C     This subroutine is a translation of the ALGOL procedure TRED2,
!C     NUM. MATH. 11, 181-195(1968) by Martin, Reinsch, and Wilkinson.
!C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
!C
!C     This subroutine reduces a DOUBLE PRECISION SYMMETRIC matrix to a
!C     symmetric tridiagonal matrix using and accumulating
!C     orthogonal similarity transformations.
!C
!C     On Input
!C
!C        NM must be set to the row dimension of the two-dimensional
!C          array parameters, A and Z, as declared in the calling
!C          program dimension statement.  NM is an INTEGER variable.
!C
!C        N is the order of the matrix A.  N is an INTEGER variable.
!C          N must be less than or equal to NM.
!C
!C        A contains the real symmetric input matrix.  Only the lower
!C          triangle of the matrix need be supplied.  A is a two-
!C          dimensional DOUBLE PRECISION array, dimensioned A(NM,N).
!C
!C     On Output
!C
!C        D contains the diagonal elements of the symmetric tridiagonal
!C          matrix.  D is a one-dimensional DOUBLE PRECISION array,
!C          dimensioned D(N).
!C
!C        E contains the subdiagonal elements of the symmetric
!C          tridiagonal matrix in its last N-1 positions.  E(1) is set
!C          to zero.  E is a one-dimensional DOUBLE PRECISION array,
!C          dimensioned E(N).
!C
!C        Z contains the orthogonal transformation matrix produced in
!C          the reduction.  Z is a two-dimensional DOUBLE PRECISION array,
!C          dimensioned Z(NM,N).
!C
!C        A and Z may coincide.  If distinct, A is unaltered.
!C
!C     Questions and comments should be directed to B. S. Garbow,
!C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!C     ------------------------------------------------------------------
!C
!C***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!C                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!C                 system Routines - EISPACK Guide, Springer-Verlag,
!C                 1976.
!C***ROUTINES CALLED  (NONE)
!C***REVISION HISTORY  (YYMMDD)
!C   760101  DATE WRITTEN
!C   890831  Modified array declarations.  (WRB)
!C   890831  REVISION DATE from Version 3.2
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   920501  Reformatted the REFERENCES section.  (WRB)
!C***END PROLOGUE  TRED2
!C
!======= test to change it into Fortran 90 format.======
!      INTEGER I,J,K,L,N,II,NM,JP1
!      DOUBLE PRECISION A(NM,*),D(*),E(*),Z(NM,*)
!      DOUBLE PRECISION F,G,H,HH,SCALE
!============
      use prec
      implicit none
      integer,intent(IN)::NM,N
      integer:: I,J,K,L,II,JP1
      real(kind=rkind),dimension(NM,N),intent(IN):: A
      real(kind=rkind),dimension(NM,N),intent(OUT):: Z
      real(kind=rkind),dimension(N),intent(OUT):: D,E
      real(kind=rkind):: F,G,H,HH,SCALE
!========================================================

!C
!C***FIRST EXECUTABLE STATEMENT  TRED2
      DO 100 I = 1, N
!C
         DO 100 J = 1, I
            Z(I,J) = A(I,J)
  100 CONTINUE
!C
      IF (N .EQ. 1) GO TO 320
!C     .......... FOR I=N STEP -1 UNTIL 2 DO -- ..........
      DO 300 II = 2, N
         I = N + 2 - II
         L = I - 1
         H = 0.0_rkind
         SCALE = 0.0_rkind
         IF (L .LT. 2) GO TO 130
!C     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
         DO 120 K = 1, L
  120    SCALE = SCALE + ABS(Z(I,K))
!C
         IF (SCALE .NE. 0.0_rkind) GO TO 140
  130    E(I) = Z(I,L)
         GO TO 290
!C
  140    DO 150 K = 1, L
            Z(I,K) = Z(I,K) / SCALE
            H = H + Z(I,K) * Z(I,K)
  150    CONTINUE
!C
         F = Z(I,L)
         G = -SIGN(SQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         Z(I,L) = F - G
         F = 0.0_rkind
!C
         DO 240 J = 1, L
            Z(J,I) = Z(I,J) / H
            G = 0.0_rkind
!C     .......... FORM ELEMENT OF A*U ..........
            DO 180 K = 1, J
  180       G = G + Z(J,K) * Z(I,K)
!C
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
!C
            DO 200 K = JP1, L
  200       G = G + Z(K,J) * Z(I,K)
!C     .......... FORM ELEMENT OF P ..........
  220       E(J) = G / H
            F = F + E(J) * Z(I,J)
  240    CONTINUE
!C
         HH = F / (H + H)
!C     .......... FORM REDUCED A ..........
         DO 260 J = 1, L
            F = Z(I,J)
            G = E(J) - HH * F
            E(J) = G
!C
            DO 260 K = 1, J
               Z(J,K) = Z(J,K) - F * E(K) - G * Z(I,K)
  260    CONTINUE
!C
  290    D(I) = H
  300 CONTINUE
!C
  320 D(1) = 0.0_rkind
      E(1) = 0.0_rkind
!C     .......... ACCUMULATION OF TRANSFORMATION MATRICES ..........
      DO 500 I = 1, N
         L = I - 1
         IF (D(I) .EQ. 0.0_rkind) GO TO 380
!C
         DO 360 J = 1, L
            G = 0.0_rkind
!C
            DO 340 K = 1, L
  340       G = G + Z(I,K) * Z(K,J)
!C
            DO 360 K = 1, L
               Z(K,J) = Z(K,J) - G * Z(K,I)
  360    CONTINUE
!C
  380    D(I) = Z(I,I)
         Z(I,I) = 1.0_rkind
         IF (L .LT. 1) GO TO 500
!C
         DO 400 J = 1, L
            Z(I,J) = 0.0_rkind
            Z(J,I) = 0.0_rkind
  400    CONTINUE
!C
  500 CONTINUE
!C
!      RETURN
!      END
      end subroutine TRED2




!====== test to change it into Fortran90 format.====
!      double precision function pythag(a,b)
!      double precision a,b
!      double precision p,r,s,t,u
!==========
       function pythag(a,b)
       use prec
       implicit none
       real(kind=rkind),intent(IN):: a,b
       real(kind=rkind):: p,r,s,t,u,pythag
!====================================================
!c
!c     finds dsqrt(a**2+b**2) without overflow or destructive underflow
!c

      p = max(abs(a),abs(b))
      if (p .eq. 0.0_rkind) go to 20
      r = (min(abs(a),abs(b))/p)**2
   10 continue
         t = 4.0_rkind + r
         if (t .eq. 4.0_rkind) go to 20
         s = r/t
         u = 1.0_rkind + 2.0_rkind*s
         p = u*p
         r = (s/u)**2.0_rkind * r
      go to 10
   20 pythag = p
!      return
!      end
    end function PYTHAG


