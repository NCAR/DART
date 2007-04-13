C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE DRSTPL(INOD,LUN,INV1,INV2,INVN)                        
                                                                        
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
                                                                        
      CHARACTER*10 TAG                                                  
      CHARACTER*3  TYP                                                  
      REAL*8       VAL                                                  
                                                                        
C-----------------------------------------------------------------------
CFPP$ EXPAND (INVWIN,USRTPL,NEWWIN)                                     
C-----------------------------------------------------------------------
                                                                        
1     NODE = INOD                                                       
2     NODE = JMPB(NODE)                                                 
      IF(NODE.EQ.0) RETURN                                              
      IF(TYP(NODE).EQ.'DRS' .OR. TYP(NODE).EQ.'DRB') THEN               
         INVN = INVWIN(NODE,LUN,INV1,INV2)                              
         IF(INVN.GT.0) THEN                                             
            CALL USRTPL(LUN,INVN,1)                                     
            CALL NEWWIN(LUN,INV1,INV2)                                  
            INVN = INVWIN(INOD,LUN,INVN,INV2)                           
            IF(INVN.GT.0) RETURN                                        
            GOTO 1                                                      
         ENDIF                                                          
      ENDIF                                                             
      GOTO 2                                                            
900   CALL BORT('DRSTPL - CANT FIND NODE:'//TAG(INOD))                 
C     print*,'drstpl:',tag(inod),':',tag(node),inv1,inv2                
C     print'(5a10)',(tag(inv(i,lun)),i=inv1,inv2)                       
      END                                                               
