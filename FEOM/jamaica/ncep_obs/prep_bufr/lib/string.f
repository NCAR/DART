C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE STRING(STR,LUN,I1,IO)
 
      PARAMETER (MXS=1000,JCONS=52)
 
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)
      COMMON /STCACH/ MSTR,NSTR,LSTR,LUX(MXS,2),USR(MXS),ICON(52,MXS)
      COMMON /USRSTR/ JCON(JCONS)
      COMMON /STORDS/ IORD(MXS),IORX(MXS)
 
      CHARACTER*(*) STR
      CHARACTER*80  USR,UST
 
C----------------------------------------------------------------------
C----------------------------------------------------------------------

      NXT = 0
      UST = STR
      IND = INODE(LUN)
      IF(LEN(STR).GT.80) GOTO 900
 
C  SEE IF STRING IS IN THE CACHE
C  -----------------------------

      DO N=1,NSTR
      IF(LUX(IORD(N),2).EQ.IND) THEN
         IORX(NXT+1) = IORD(N)
         NXT = NXT+1
      ENDIF
      ENDDO
      DO N=1,NXT 
      IF(UST.EQ.USR(IORX(N)))GOTO1
      ENDDO
      GOTO2
 
C  IF IT IS COPY PARAMETERS FROM THE CACHE
C  ---------------------------------------
 
1     DO J=1,JCONS
      JCON(J) = ICON(J,IORX(N))
      ENDDO
      GOTO 100
 
C  IF NOT PARSE IT AND PUT IT THERE
C  --------------------------------
 
2     CALL PARUSR(STR,LUN,I1,IO)
      LSTR = MAX(MOD(LSTR+1,MSTR+1),1)
      NSTR = MIN(NSTR+1,MSTR)
      LUX(LSTR,1) = LUN
      LUX(LSTR,2) = IND
      USR(LSTR) = STR
      DO J=1,JCONS
      ICON(J,LSTR) = JCON(J)
      ENDDO
 
C  REARRANGE THE CACHE ORDER AFTER AN UPDATE
C  -----------------------------------------
 
      DO N=NSTR,2,-1
      IORD(N) = IORD(N-1)
      ENDDO
      IORD(1) = LSTR
 
C  NORMAL AND ERROR EXITS
C  ----------------------
 
100   IF(JCON(1).GT.I1) GOTO 901
      RETURN
900   CALL BORT('STRING - USER STRING > 80 CHARS         :'//UST)
901   CALL BORT('STRING - MUST BE AT LEAST I1 STORE NODES:'//UST)
      END
