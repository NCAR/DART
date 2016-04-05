C---------------------------------------------------------------------- 
C  UPDATE THE MESSAGE BUFFER WITH NEW SUBSET                            
C---------------------------------------------------------------------- 
      SUBROUTINE CPYUPD(LUNIT,LIN,LUN,IBYT)                             
                                                                        
      COMMON /MSGPTR/ NBY0,NBY1,NBY2,NBY3,NBY4,NBY5                     
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)     
                                                                        
      CHARACTER*8 CBAY                                                  
      EQUIVALENCE (CBAY,JBAY)                                           
                                                                        
C-----------------------------------------------------------------------
CFPP$ EXPAND(MVB,PKB)                                                   
C-----------------------------------------------------------------------
                                                                        
                                                                        
C  SEE IF THE NEW SUBSET FITS                                           
C  --------------------------                                           
                                                                        
      IF(MBYT(LUN)+IBYT.GT.MAXBYT) THEN                                 
         CALL MSGWRT(LUNIT,MBAY(1,LUN),MBYT(LUN))                       
         CALL MSGINI(LUN)                                               
      ENDIF                                                             
                                                                        
      IF(IBYT.GT.MAXBYT-MBYT(LUN)) GOTO 900                             
                                                                        
C  TRANSFER SUBSET FROM ONE MESSAGE TO THE OTHER                        
C  ---------------------------------------------                        
                                                                        
      CALL MVB(MBAY(1,LIN),MBYT(LIN)+1,MBAY(1,LUN),MBYT(LUN)-3,IBYT)    
                                                                        
C  UPDATE THE SUBSET AND BYTE COUNTERS                                  
C  --------------------------------------                               
                                                                        
      MBYT(LUN)   = MBYT(LUN)   + IBYT                                  
      NSUB(LUN)   = NSUB(LUN)   + 1                                     
                                                                        
      LBIT = (NBY0+NBY1+NBY2+4)*8                                       
      CALL PKB(NSUB(LUN),16,MBAY(1,LUN),LBIT)                           
                                                                        
      LBYT = NBY0+NBY1+NBY2+NBY3                                        
      NBYT = IUPB(MBAY(1,LUN),LBYT+1,24)                                
      LBIT = LBYT*8                                                     
      CALL PKB(NBYT+IBYT,24,MBAY(1,LUN),LBIT)                           
                                                                        
      RETURN                                                            
                                                                        
900   CALL BORT('CPYUPD - SUBSET LONGER THAN ANY POSSIBLE MESSAGE')    
      END                                                               
