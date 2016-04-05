       SUBROUTINE ALPBET2(A,B,C,D,K,N,IR,IC,ALPH,BET,AUX,FX,IPVT,Z)

C  COMPUTES ALPHAS AND BETAS FOR USE IN LINDZEN-KUO ALGORITHM
C  AND STORES ON UNIT 4

C  Version 2.0 created by J. Forbes 2/3/91
C  Modified on 3/14/91 to include NCAR/LINPACK routines
C  Modified on 3/30/04 to use either LINPACK or LAPACK routines (TJH)

        INTEGER IPVT(IC)
        COMPLEX Z(IC)
        COMPLEX A(IR,IC),B(IR,IC),C(IR,IC),D(IR)
        COMPLEX ALPH(IR,IC),BET(IR),AUX(IR,IC),FX(IR)

        integer info

C  UPPER BOUNDARY CONDITION; COMPUTE ALPHA ZERO AND BETA ZERO

        IF(K-1) 1,1,2

#ifdef LINPACK

        ! IR == lda
        ! IC ==   n
1       CALL CGECO(A,IR,IC,IPVT,RCOND,Z)
        CALL CGESL(A,IR,IC,IPVT,D,0)
        DO 20 J=1,IC
        CALL CGESL(A,IR,IC,IPVT,B(1,J),0)
20      CONTINUE

#elif LAPACK


        ! IR == lda == m
        ! IC == n

        !ALL CGETRF( m,  n, A, lda, ipvt, info)
1       CALL CGETRF(IR, IC, A,  IR, IPVT, info)
        if ( info .ne. 0 ) then
           write(*,*)'ERROR from LAPACK:cgetrf ... info is ',info 
           stop
        endif

        !    CGETRs('n',  n, nrhs, a, lda, ipiv, b, ldb, info)
        CALL CGETRS('N', IC,    1, A,  IR, IPVT, D,  IR, info)
        CALL CGETRS('N', IC,   IC, A,  IR, IPVT, B,  IR, info)

#endif


        DO 11 KI=1,IR
        BET(KI)=D(KI)
        DO 11 JI=1,IC
11      ALPH(KI,JI)=-B(KI,JI)
        WRITE(UNIT=4,REC=K) ((ALPH(I,J),I=1,IR),J=1,IC),(BET(I),I=1,IR)
        RETURN

C  INTERIOR POINTS

2       CONTINUE

C  ALPHA(N-1) AND BETA(N-1) WERE SAVED FROM PREVIOUS CALL AND NOW USED:

        CALL MPRDD(A,ALPH,AUX,IR,IC,IC)
        CALL MPRDDV(A,BET,FX,IR,IC)
        CALL GMADDM(AUX,B,AUX,IR,IC,1.) 
        CALL GMADDV(D,FX,BET,IR,-1.)

C       Factor and solve a linear system ... TJH

#ifdef LINPACK

c     sgeco factors a real matrix by gaussian elimination
c     and estimates the condition of the matrix.
c
c     if  rcond  is not needed, sgefa is slightly faster.
c     to solve  a*x = b , follow sgeco by sgesl.
c     to compute  inverse(a)*c , follow sgeco by sgesl.
c     to compute  determinant(a) , follow sgeco by sgedi.
c     to compute  inverse(a) , follow sgeco by sgedi.
c
c     on entry
c
c        a       real(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        rcond   real
c                an estimate of the reciprocal condition of  a .
c                for the system  a*x = b , relative perturbations
c                in  a  and  b  of size  epsilon  may cause
c                relative perturbations in  x  of size  epsilon/rcond .
c                if  rcond  is so small that the logical expression
c                           1.0 + rcond .eq. 1.0
c                is true, then  a  may be singular to working
c                precision.  in particular,  rcond  is zero  if
c                exact singularity is detected or the estimate
c                underflows.
c
c        z       real(n)
c                a work vector whose contents are usually unimportant.
c                if  a  is close to a singular matrix, then  z  is
c                an approximate null vector in the sense that
c                norm(a*z) = rcond*norm(a)*norm(z) 


c     sgesl solves the real system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by sgeco or sgefa.
c
c     on entry
c
c        a       real(lda, n)
c                the output from sgeco or sgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from sgeco or sgefa.
c
c        b       real(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b  where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if sgeco has set rcond .gt. 0.0
c        or sgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call sgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call sgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue

        CALL CGECO(AUX,IR,IC,IPVT,RCOND,Z)
        CALL CGESL(AUX,IR,IC,IPVT,BET,0)
        DO 10 J=1,IC
           CALL CGESL(AUX,IR,IC,IPVT,C(1,J),0)
10      CONTINUE

#elif  LAPACK

        CALL CGETRF(IR, IC, AUX, IR, IPVT, info)
        if ( info .ne. 0 ) then
           write(*,*)'ERROR from LAPACK:cgetrf ... info is ',info 
           stop
        endif

        CALL CGETRS('N', IC,    1, AUX,  IR, IPVT, BET,  IR, info)
        CALL CGETRS('N', IC,   IC, AUX,  IR, IPVT,   C,  IR, info)

#else

        call NEED_A_CHOLESKY_AND_SOLVE_ROUTINE()

#endif


        DO 12 KI=1,IR
        DO 12 JI=1,IC
12      ALPH(KI,JI)=-C(KI,JI)

        WRITE(UNIT=4,REC=K) ((ALPH(I,J),I=1,IR),J=1,IC),(BET(I),I=1,IR)

        RETURN
        END
