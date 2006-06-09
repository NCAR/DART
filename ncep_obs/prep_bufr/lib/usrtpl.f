C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE USRTPL(LUN,INVN,NBMP)                                  
                                                                        
      PARAMETER (MAXINV=15000)                                          
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
                                                                        
      CHARACTER*10 TAG                                                  
      CHARACTER*3  TYP                                                  
      DIMENSION    ITMP(MAXINV),VTMP(MAXINV)                            
      LOGICAL      DRP,DRS,DRB,DRX                                      
      REAL*8       VAL,VTMP                                             
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     PRINT*,'USRTPL:',LUN,':',INVN,':',NBMP,':',tag(inode(lun))        
                                                                        
      IF(NBMP.LE.0) RETURN                                              
                                                                        
      DRP = .FALSE.                                                     
      DRS = .FALSE.                                                     
      DRX = .FALSE.                                                     
                                                                        
C  SET UP A NODE EXPANSION                                              
C  -----------------------                                              
                                                                        
      IF(INVN.EQ.1) THEN                                                
         NODI = INODE(LUN)                                              
         INV(1,LUN) = NODI                                              
         NVAL(LUN)  = 1                                                 
         IF(NBMP.NE.1) GOTO 900                                         
      ELSEIF(INVN.GT.0 .AND. INVN.LE.NVAL(LUN)) THEN                    
         NODI = INV(INVN,LUN)                                           
         DRP  = TYP(NODI) .EQ. 'DRP'                                    
         DRS  = TYP(NODI) .EQ. 'DRS'                                    
         DRB  = TYP(NODI) .EQ. 'DRB'                                    
         DRX  = DRP .OR. DRS .OR. DRB                                   
         IVAL = VAL(INVN,LUN)                                           
         JVAL = 2**IBT(NODI)-1                                          
         VAL(INVN,LUN) = IVAL+NBMP                                      
         IF(DRB.AND.NBMP.NE.1) GOTO 900                                 
         IF(.NOT.DRX         ) GOTO 901                                 
         IF(IVAL.LT.0.       ) GOTO 902                                 
         IF(IVAL+NBMP.GT.JVAL) GOTO 903                                 
      ELSE                                                              
         GOTO 904                                                       
      ENDIF                                                             
                                                                        
C  RECALL A PRE-FAB NODE EXPANSION SEGMENT                              
C  ---------------------------------------                              
                                                                        
      NEWN = 0                                                          
      N1 = ISEQ(NODI,1)                                                 
      N2 = ISEQ(NODI,2)                                                 
                                                                        
      IF(N1.EQ.0          ) GOTO 905                                    
      IF(N2-N1+1.GT.MAXINV) GOTO 906                                    
                                                                        
      DO N=N1,N2                                                        
      NEWN = NEWN+1                                                     
      ITMP(NEWN) = JSEQ(N)                                              
      VTMP(NEWN) = VALI(JSEQ(N))                                        
      if(vtmp(newn).gt.10e9) vtmp(newn) = 10e10                         
      ENDDO                                                             
                                                                        
C  MOVE OLD NODES - STORE NEW ONES                                      
C  -------------------------------                                      
                                                                        
      IF(NVAL(LUN)+NEWN*NBMP.GT.MAXINV) print*,'@:',nval(lun)+newn*nbmp 
      IF(NVAL(LUN)+NEWN*NBMP.GT.MAXINV) GOTO 907                        
                                                                        
CDIR$ IVDEP                                                             
      DO J=NVAL(LUN),INVN+1,-1                                          
      INV(J+NEWN*NBMP,LUN) = INV(J,LUN)                                 
      VAL(J+NEWN*NBMP,LUN) = VAL(J,LUN)                                 
      ENDDO                                                             
                                                                        
      IF(DRP.OR.DRS) VTMP(1) = NEWN                                     
      KNVN = INVN                                                       
                                                                        
      DO I=1,NBMP                                                       
      DO J=1,NEWN                                                       
      KNVN = KNVN+1                                                     
      INV(KNVN,LUN) = ITMP(J)                                           
      VAL(KNVN,LUN) = VTMP(J)                                           
      ENDDO                                                             
      ENDDO                                                             
                                                                        
C  RESET POINTERS AND COUNTERS                                          
C  ---------------------------                                          
                                                                        
      NVAL(LUN) = NVAL(LUN) + NEWN*NBMP                                 
                                                                        
C     print*,tag(inv(invn,lun)),' ',newn,' ',nbmp,' ',nval(lun)         
C     DO I=1,NEWN                                                       
C     PRINT*,TAG(ITMP(I))                                               
C     ENDDO                                                             
                                                                        
                                                                        
      IF(DRX) THEN                                                      
         NODE = NODI                                                    
         INVR = INVN                                                    
4        NODE = JMPB(NODE)                                              
         IF(NODE.GT.0) THEN                                             
            IF(ITP(NODE).EQ.0) THEN                                     
               DO INVR=INVR-1,1,-1                                      
               IF(INV(INVR,LUN).EQ.NODE) THEN                           
                  VAL(INVR,LUN) = VAL(INVR,LUN)+NEWN*NBMP               
                  GOTO 4                                                
               ENDIF                                                    
               ENDDO                                                    
               GOTO 909                                                 
            ELSE                                                        
               GOTO 4                                                   
            ENDIF                                                       
         ENDIF                                                          
      ENDIF                                                             
                                                                        
C  NORMAL EXIT                                                          
C  -----------                                                          
                                                                        
      RETURN                                                            
                                                                        
C  ERROR EXITS                                                          
C  -----------                                                          
                                                                        
900   CALL BORT('USRTPL - NBMP <> 1 FOR        : '//TAG(NODI))         
901   CALL BORT('USRTPL - NODE NOT SUB,DRP,DRS : '//TAG(NODI))         
902   CALL BORT('USRTPL - NEGATIVE REP FACTOR  : '//TAG(NODI))         
903   CALL BORT('USRTPL - REP FACTOR OVERFLOW  : '//TAG(NODI))         
904   CALL BORT('USRTPL - INVENTORY INDEX OUT OF BOUNDS     ')         
905   CALL BORT('USRTPL - UNSET EXPANSION SEG  : '//TAG(NODI))         
906   CALL BORT('USRTPL - TEMP ARRAY OVERFLOW  : '//TAG(NODI))         
907   CALL BORT('USRTPL - INVENTORY OVERFLOW   : '//TAG(NODI))         
908   CALL BORT('USRTPL - TPL CACHE OVERFLOW   : '//TAG(NODI))         
909   CALL BORT('USRTPL - BAD BACKUP STRATEGY  : '//TAG(NODI))         
      END                                                               
