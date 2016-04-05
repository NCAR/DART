C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      FUNCTION IRDERM(LUNIT,MBAY)                                       
                                                                        
      CHARACTER*4  SEVN                                            
      CHARACTER*1  BAY(5000*8)
      CHARACTER*36 SEC013,SECSAV
      DIMENSION    MBAY(5000),KBAY(5000)                                
      EQUIVALENCE  (BAY(1),KBAY(1),SEC013)                       
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      DO I=1,5000                                                       
      KBAY(I) = 0                                                       
      ENDDO                                                             
      IBSKIP = 0
                                                                        
C  FIND A BUFR MESSAGE                                                  
C  -------------------                                                  
                                                                        
1     READ(LUNIT,END=100,ERR=100) SEC013(1:32)
      LASTRD = 32
2     IBUFR = INDEX(SEC013(1:32),'BUFR')
      IF(IBUFR.EQ.0) THEN                                  
         IBSKIP = IBSKIP + LASTRD
         IF(INDEX(SEC013(1:32),'B').EQ.0) THEN
            READ(LUNIT,END=100,ERR=100) SEC013(1:32)                    
            LASTRD = 32
         ELSE IF(SEC013(32:32).EQ.'B') THEN
            SEC013(1:1) = 'B'
            READ(LUNIT,END=100,ERR=100) SEC013(2:32)
            LASTRD = 31
         ELSE IF(INDEX(SEC013(1:32),'BU').EQ.0) THEN
            READ(LUNIT,END=100,ERR=100) SEC013(1:32)
            LASTRD = 32
         ELSE IF(SEC013(31:32).EQ.'BU') THEN
            SEC013(1:2) = 'BU'
            READ(LUNIT,END=100,ERR=100) SEC013(3:32)
            LASTRD = 30
         ELSE IF(INDEX(SEC013(1:32),'BU').LT.30) THEN
            READ(LUNIT,END=100,ERR=100) SEC013(1:32)                    
            LASTRD = 32
         ELSE IF(SEC013(30:32).EQ.'BUF') THEN
            SEC013(1:3) = 'BUF'
            READ(LUNIT,END=100,ERR=100) SEC013(4:32)
            LASTRD = 29
         ELSE
            READ(LUNIT,END=100,ERR=100) SEC013(1:32)
            LASTRD = 32
         ENDIF
         GOTO 2                                                         
      ELSE IF(IBUFR.GT.1) THEN
         SECSAV = SEC013
         SEC013(1:33-IBUFR) = SECSAV(IBUFR:32)
         READ(LUNIT,END=100,ERR=100) SEC013(34-IBUFR:32)
         IBSKIP = IBSKIP + IBUFR - 1
      ENDIF                                                             
                                                                        
C  IF THIS IS BUFR RELEASE 0, THE FOLLOWING WILL ACCOUNT FOR SECTION-1
C  -------------------------------------------------------------------

      J = IUPM(BAY(5),24)
      I = J + 8

      IF(J.LE.32) THEN

C  DETERMINE WHETHER BYTE COUNT IS FOR SECTION-0 OR SECTION-1
C  ----------------------------------------------------------

         SEVN = SEC013(J-3:J)
         IF(SEVN.EQ.'7777') THEN
            IF(J.LT.32) THEN
               SECSAV = SEC013
               SEC013(1:32-J) = SECSAV(J+1:32)
            ENDIF
            READ(LUNIT,END=100,ERR=100) SEC013(33-J:32)
            PRINT*,'SHORT RECORD SKIPPED'                    
            GOTO 2
         ENDIF

C  IF THIS IS BUFR RELEASE 0, SHIFT BYTES 5 AND UP 4 BYTES TO THE RIGHT 
C  --------------------------------------------------------------------

         SECSAV = SEC013
         SEC013(9:36) = SECSAV(5:32)

C  IF SECTION-2 IS ABSENT, BEGIN WITH SECTION-3
C  --------------------------------------------

         IF(IUPM(BAY(16),1).GT.0) THEN         
            KCONT = 2
         ELSE
            KCONT = 3
         ENDIF

C  DETERMINE WHETHER SECTION-2 AND SECTION-3 ARE IN BYTES 9-36
C  -----------------------------------------------------------

         KREST = 0
         DO K=KCONT,3
         IF(I.LE.33) THEN
            I = I + IUPM(BAY(I+1),24)                                
            IF(I.GT.36) THEN
               READ(LUNIT,END=100,ERR=100) (BAY(J),J=37,I)
            ENDIF
            KREST = K + 1
         ELSE IF(I.LE.35) THEN
            READ(LUNIT,END=100,ERR=100) (BAY(J),J=37,I+3),
     .                  (BAY(J),J=I+4,I+IUPM(BAY(I+1),24))           
            I = I + IUPM(BAY(I+1),24)                                
            KREST = K + 1
         ENDIF                                                       
         ENDDO                                                          
         IF(KREST.NE.0) KCONT = KREST
      ELSE
         READ(LUNIT,END=100,ERR=100) (BAY(K),K=33,J)                

C  DETERMINE WHETHER BYTE COUNT IS FOR SECTION-0 OR SECTION-1
C  ----------------------------------------------------------

         SEVN = BAY(J-3)//BAY(J-2)//BAY(J-1)//BAY(J)
         IF(SEVN.EQ.'7777') GOTO 50

C  IF THIS IS BUFR RELEASE 0, SHIFT BYTES 5 AND UP 4 BYTES TO THE RIGHT 
C  --------------------------------------------------------------------

         READ(LUNIT,END=100,ERR=100) (BAY(K),K=J+5,J+8)                
         DO K=J,33,-1
         BAY(K+4) = BAY(K)
         ENDDO
         SECSAV = SEC013
         SEC013(9:36) = SECSAV(5:32)

C  IF SECTION-2 IS ABSENT, BEGIN WITH SECTION-3
C  --------------------------------------------

         IF(IUPM(BAY(16),1).GT.0) THEN         
            KCONT = 2
         ELSE
            KCONT = 3
         ENDIF

      ENDIF

C  FOR REMAINING SECTIONS (UP TO SECTION-4) READ BYTE COUNT AND BYTES
C  ------------------------------------------------------------------

      DO K=KCONT,4                                                      
      READ(LUNIT,END=100,ERR=100) (BAY(J),J=I+1,I+3),                
     .            (BAY(J),J=I+4,I+IUPM(BAY(I+1),24))                 
      I = I + IUPM(BAY(I+1),24)                                      
      ENDDO                                                             
                                                                        
C  CHECK ON SECTION 5 FOR BAD RECORD INDICATOR                          
C  -------------------------------------------                          
                                                                        
      READ(LUNIT,END=100,ERR=100) SEVN                                  
      IF(SEVN.NE.'7777') THEN
         PRINT*,'BAD RECORD SKIPPED'                    
         GOTO 1                                         
      ENDIF
      I = I+4                                                           
                                                                        
C  FILL IN THE ARRAY TO RETURN                                          
C  ---------------------------                                          
                                                                        
50    DO I=1,5000                                                       
      MBAY(I) = KBAY(I)                                                 
      ENDDO                                                             
                                                                        
      IRDERM =  0                                                       
      RETURN                                                            
                                                                        
100   IRDERM = -1                                                       
      RETURN                                                            
      END                                                               
