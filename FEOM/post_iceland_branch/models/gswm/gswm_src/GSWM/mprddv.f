	SUBROUTINE MPRDDV(A,B,R,N,M)

C  Version 2.0 created by J. Forbes 2/6/91

	COMPLEX A(N,M),B(M),R(N)
	DO 1 IR=1,N
	R(IR)=(0.0,0.0)
	DO 1 IC=1,M
    1   R(IR)=A(IR,IC)*B(IC)+R(IR)
	RETURN
	END
