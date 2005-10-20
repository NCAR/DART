C-----------------------------------------------------------------------
C  SUBROUTINE STANDARD WILL REWRITE A MESSAGE WRITTEN BY THE NCEP BUFR  
C  INTERFACE PROGRAM INTO A MORE STANDARD BUFR FORM. SECTION THREE IS   
C  REWRITTEN AS THE EXPANSION (1 LEVEL DEEP) OF THE NCEP SUBSET SEQUENCE
C  DESCRIPTOR. SECTION FOUR IS REWRITTEN TO CONFORM TO THE NEW SECTION  
C  THREE DESCRIPTOR LIST. ALL SUBSET BYTE COUNTERS AND BIT PADS ARE REMO
C  FROM THE DATA. IF THE 1 LEVEL EXPANSION OF THE SUBSET SEQUENCE DESCRI
C  CONTAINS ONLY STANDARD DESCRIPTORS, THEN THE NEW MESSAGE IS ENTIRLY  
C  AND STRICTLY STANDARD.                                               
C                                                                       
C  THE SUBROUTINE ARGUMENTS ARE:                                        
C                                                                       
C  INPUT:  LUNIT - UNIT OPENED WITH BUFR TABLES USING OPENBF            
C          MSGIN - ARRAY CONTAINING AN NCEP BUFR MESSAGE                
C                                                                       
C  OUTPUT: MSGOT - ARRAY CONTAINING STANDARDIZED FORM OF THE INPUT MESSA
C                                                                       
C  ----------------------------------------------                       
C  NOTE: MSGIN AND MSGOT MUST BE SEPARATE ARRAYS.                       
C  ----------------------------------------------                       
C                                                                       
C-----------------------------------------------------------------------
      SUBROUTINE STANDARD(LUNIT,MSGIN,MSGOT)                            
                                                                        
      DIMENSION MSGIN(*),MSGOT(*)                                       
                                                                        
      CHARACTER*8 SUBSET                                                
      CHARACTER*4 SEVN                                                  
      CHARACTER*1 TAB                                                   
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  LUNIT MUST POINT TO AN OPEN BUFR FILE                                
C  -------------------------------------                                
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) CALL BORT('STANDARD - NO OPEN BUFR FILE!!')          
                                                                        
C  IDENTIFY THE SECTION LENGTHS AND ADDRESSES IN MSGIN                  
C  ---------------------------------------------------                  
                                                                        
      IAD0 = 0                                                          
      LEN0 = 8                                                          
      LENN = IUPB(MSGIN,IAD0+5,24)                                      
                                                                        
      IAD1 = IAD0+LEN0                                                  
      LEN1 = IUPB(MSGIN,IAD1+1,24)                                      
      LEN2 = IUPB(MSGIN,IAD1+8,1)                                       
                                                                        
      IAD2 = IAD1+LEN1                                                  
      LEN2 = IUPB(MSGIN,IAD2+1,24)*LEN2                                 
                                                                        
      IAD3 = IAD2+LEN2                                                  
      LEN3 = IUPB(MSGIN,IAD3+1,24)                                      
                                                                        
      IAD4 = IAD3+LEN3                                                  
      LEN4 = IUPB(MSGIN,IAD4+1,24)                                      
                                                                        
      LENM = LEN0+LEN1+LEN2+LEN3+LEN4+4                                 
                                                                        
      IF(LENN.NE.LENM) CALL BORT('STANDARD - BAD INPUT BYTE COUNTS')   
                                                                        
      MBIT = (LENN-4)*8                                                 
      CALL UPC(SEVN,4,MSGIN,MBIT)                                       
      IF(SEVN.NE.'7777') CALL BORT('STANDARD - CANT FIND 7777')        
                                                                        
C  COPY SECTIONS 0 THROUGH PART OF SECTION 3 INTO MSGOT                 
C  ----------------------------------------------------                 
                                                                        
      CALL MVB(MSGIN,1,MSGOT,1,LEN0+LEN1+LEN2+7)                        
                                                                        
C  REWRITE NEW SECTION 3 IN A "STANDARD" FORM                           
C  ------------------------------------------                           
                                                                        
      NSUB = IUPB(MSGIN,IAD3+ 5,16)                                     
      ISUB = IUPB(MSGIN,IAD3+10,16)                                     
      IBIT = (IAD3+7)*8                                                 
                                                                        
C  LOOK UP THE SUBSET DESCRIPTOR AND ITS LENGTH IN DESCRIPTORS          
C  -----------------------------------------------------------          
                                                                        
      CALL NUMTAB(LUN,ISUB,SUBSET,TAB,ITAB)                             
      IF(ITAB.EQ.0) CALL BORT('STANDARD - UNKNOWN SUBSET DESCRIPTOR')  
      CALL UPTDD(ITAB,LUN,0,NSEQ)                                       
                                                                        
C  COPY EACH DESCRIPTOR IN THE SUBSET SEQUENCE INTO THE NEW SECTION 3   
C  ------------------------------------------------------------------   
                                                                        
      DO N=1,NSEQ                                                       
      CALL UPTDD(ITAB,LUN,N,IDSC)                                       
      CALL PKB(IDSC,16,MSGOT,IBIT)                                      
      IF(N.EQ.NSEQ) CALL PKB(0,8,MSGOT,IBIT)                            
      ENDDO                                                             
                                                                        
      IBIT = IAD3*8                                                     
      LEN3 = 8+NSEQ*2                                                   
      NAD4 = IAD3+LEN3                                                  
      CALL PKB(LEN3,24,MSGOT,IBIT)                                      
                                                                        
C  NOW THE TRICKY PART - NEW SECTION 4                                  
C  -----------------------------------                                  
                                                                        
      IBIT = (IAD4+4)*8                                                 
      JBIT = (NAD4+4)*8                                                 
                                                                        
C  COPY THE SUBSETS, MINUS THE BYTE COUNTER AND PAD, INTO THE NEW SECTIO
C  ---------------------------------------------------------------------
                                                                        
      DO 10 I=1,NSUB                                                    
      CALL UPB(LSUB,16,MSGIN,IBIT)                                      
                                                                        
      DO L=1,LSUB-2                                                     
      CALL UPB(NVAL,8,MSGIN,IBIT)                                       
      CALL PKB(NVAL,8,MSGOT,JBIT)                                       
      ENDDO                                                             
                                                                        
      DO K=1,8                                                          
      KBIT = IBIT-K-8                                                   
      CALL UPB(KVAL,8,MSGIN,KBIT)                                       
      IF(KVAL.EQ.K) THEN                                                
         JBIT = JBIT-K-8                                                
         GOTO 10                                                        
      ENDIF                                                             
      ENDDO                                                             
      CALL BORT('STANDARD - KBIT ERROR')                               
                                                                        
10    ENDDO                                                             
                                                                        
C  MAKE SURE NEW SECTION 4 HAS AN EVEN NUMBER OF BYTES AND ENTER THE COU
C  ---------------------------------------------------------------------
                                                                        
      DO WHILE(.NOT.(MOD(JBIT,8).EQ.0 .AND. MOD(JBIT/8,2).EQ.0))        
         CALL PKB(0,1,MSGOT,JBIT)                                       
      ENDDO                                                             
                                                                        
      IBIT = NAD4*8                                                     
      LEN4 = JBIT/8 - NAD4                                              
      CALL PKB(LEN4,24,MSGOT,IBIT)                                      
                                                                        
C  FINISH THE NEW MESSAGE WITH AN UPDATED SECTION-0 BYTE COUNT          
C  -----------------------------------------------------------          
                                                                        
      IBIT = 32                                                         
      LENM = LEN0+LEN1+LEN2+LEN3+LEN4+4                                 
      CALL PKB(LENM,24,MSGOT,IBIT)                                      
                                                                        
C  NORMAL EXIT                                                          
C  -----------                                                          
                                                                        
      CALL PKC('7777', 4,MSGOT,JBIT)                                    
                                                                        
      RETURN                                                            
      END                                                               
