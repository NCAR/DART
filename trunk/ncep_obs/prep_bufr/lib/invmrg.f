C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE INVMRG(LUBFI,LUBFJ)                                    
                                                                        
      COMMON /MRGCOM/ NRPL,NMRG,NAMB,NTOT                               
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
                                                                        
      CHARACTER*10  TAG                                                 
      CHARACTER*3   TYP                                                 
      LOGICAL       HEREI,HEREJ,MISSI,MISSJ,SAMEI                       
      REAL*8        VAL,BMISS                                           
                                                                        
      DATA BMISS /10E10/                                                
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      IS = 1                                                            
      JS = 1                                                            
                                                                        
C  GET THE UNIT POINTERS                                                
C  ---------------------                                                
                                                                        
      CALL STATUS(LUBFI,LUNI,IL,IM)                                     
      CALL STATUS(LUBFJ,LUNJ,JL,JM)                                     
                                                                        
C  STEP THROUGH THE BUFFERS COMPARING THE INVENTORY AND MERGING DATA    
C  -----------------------------------------------------------------    
                                                                        
      DO WHILE(IS.LE.NVAL(LUNI))                                        
                                                                        
C  CHECK TO SEE WE ARE AT THE SAME NODE IN EACH BUFFER                  
C  ---------------------------------------------------                  
                                                                        
      NODE = INV(IS,LUNI)                                               
      NODJ = INV(JS,LUNJ)                                               
      IF(NODE.NE.NODJ) CALL BORT('TABULAR MISMATCH')                   
                                                                        
      ITYP = ITP(NODE)                                                  
                                                                        
C  FOR TYPE 1 NODES DO AN ENTIRE SEQUENCE REPLACEMENT                   
C  --------------------------------------------------                   
                                                                        
      IF(ITYP.EQ.1) THEN                                                
         IF(TYP(NODE).EQ.'DRB') IOFF = 0                                
         IF(TYP(NODE).NE.'DRB') IOFF = 1                                
         IWRDS = NWORDS(IS,LUNI)+IOFF                                   
         JWRDS = NWORDS(JS,LUNJ)+IOFF                                   
         IF(IWRDS.GT.IOFF .AND. JWRDS.EQ.IOFF) THEN                     
CDIR$       IVDEP                                                       
            DO N=NVAL(LUNJ),JS+1,-1                                     
            INV(N+IWRDS-JWRDS,LUNJ) = INV(N,LUNJ)                       
            VAL(N+IWRDS-JWRDS,LUNJ) = VAL(N,LUNJ)                       
            ENDDO                                                       
            DO N=0,IWRDS                                                
            INV(JS+N,LUNJ) = INV(IS+N,LUNI)                             
            VAL(JS+N,LUNJ) = VAL(IS+N,LUNI)                             
            ENDDO                                                       
            NVAL(LUNJ) = NVAL(LUNJ)+IWRDS-JWRDS                         
            JWRDS = IWRDS                                               
            NRPL = NRPL+1                                               
         ENDIF                                                          
         IS = IS+IWRDS                                                  
         JS = JS+JWRDS                                                  
      ENDIF                                                             
                                                                        
C  FOR TYPES 2 AND 3 FILL MISSINGS                                      
C  -------------------------------                                      
                                                                        
      IF(ITYP.EQ.2) THEN                                                
         HEREI = VAL(IS,LUNI).LT.BMISS                                  
         HEREJ = VAL(JS,LUNJ).LT.BMISS                                  
         MISSI = VAL(IS,LUNI).GE.BMISS                                  
         MISSJ = VAL(JS,LUNJ).GE.BMISS                                  
         SAMEI = VAL(IS,LUNI).EQ.VAL(JS,LUNJ)                           
         IF(HEREI.AND.MISSJ) THEN                                       
            VAL(JS,LUNJ) = VAL(IS,LUNI)                                 
            NMRG = NMRG+1                                               
         ELSEIF(HEREI.AND.HEREJ.AND..NOT.SAMEI) THEN                    
            NAMB = NAMB+1                                               
         ENDIF                                                          
      ENDIF                                                             
                                                                        
      IF(ITYP.EQ.3) THEN                                                
         HEREI = VAL(IS,LUNI).NE.BMISS                                  
         HEREJ = VAL(JS,LUNJ).NE.BMISS                                  
         MISSI = VAL(IS,LUNI).EQ.BMISS                                  
         MISSJ = VAL(JS,LUNJ).EQ.BMISS                                  
         SAMEI = VAL(IS,LUNI).EQ.VAL(JS,LUNJ)                           
         IF(HEREI.AND.MISSJ) THEN                                       
            VAL(JS,LUNJ) = VAL(IS,LUNI)                                 
            NMRG = NMRG+1                                               
         ELSEIF(HEREI.AND.HEREJ.AND..NOT.SAMEI) THEN                    
            NAMB = NAMB+1                                               
         ENDIF                                                          
      ENDIF                                                             
                                                                        
C  BUMP THE COUNTERS AND GO CHECK THE NEXT PAIR                         
C  --------------------------------------------                         
                                                                        
      IS = IS + 1                                                       
      JS = JS + 1                                                       
      ENDDO                                                             
                                                                        
      NTOT = NTOT+1                                                     
                                                                        
      RETURN                                                            
                                                                        
C  ENTRY MRGINV PRINTS A SUMMARY OF MERGE ACTIVITY                      
C  -----------------------------------------------                      
                                                                        
      ENTRY MRGINV                                                      
                                                                        
      PRINT*                                                            
      PRINT*,'-----------------------------'                            
      PRINT*,'INVENTORY FROM INVMRG PROCESS'                            
      PRINT*,'-----------------------------'                            
      PRINT*,'DRB EXPANSION= ',NRPL                                     
      PRINT*,'MERGES       = ',NMRG                                     
      PRINT*,'AMBIGUOUS    = ',NAMB                                     
      PRINT*,'-----------------------------'                            
      PRINT*,'TOTAL VISITS = ',NTOT                                     
      PRINT*,'-----------------------------'                            
      PRINT*                                                            
                                                                        
      RETURN                                                            
      END                                                               
