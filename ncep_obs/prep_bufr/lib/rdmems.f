C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE RDMEMS(ISUB,IRET)                                      
                                                                        
      PARAMETER (MAXMSG=50000,MAXMEM=16000000)                           
                                                                        
      COMMON /MSGMEM/ MUNIT,MLAST,MSGP(0:MAXMSG),MSGS(MAXMEM)           
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)     
      COMMON /UNPTYP/ MSGUNP(32)
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      CALL STATUS(MUNIT,LUN,IL,IM)                                      
      IF(IL       .GE.0) GOTO 900                                       
      IF(IM       .EQ.0) GOTO 901                                       
      IF(NSUB(LUN).NE.0) GOTO 902                                       
                                                                        
      IF(ISUB.GT.MSUB(LUN)) IRET = -1                                   
      IF(ISUB.GT.MSUB(LUN)) RETURN                                      
                                                                        
      MBYM = MBYT(LUN)                                                  
      NBYT = 0                                                          
                                                                        
C  POSITION TO READ SUBSET #ISUB FROM MEMORY MESSAGE #IMSG            
C  -------------------------------------------------------                          
      IF(MSGUNP(LUN).EQ.0) THEN
         NSUB(LUN) = ISUB-1
         DO I=1,ISUB-1
         MBYT(LUN) = MBYT(LUN) + IUPB(MBAY(1,LUN),MBYT(LUN)+1,16)
         ENDDO
      ELSEIF(MSGUNP(LUN).EQ.1) THEN
         DO I=1,ISUB-1
         CALL READSB(MUNIT,IRET)
         ENDDO
      ELSEIF(MSGUNP(LUN).EQ.2) THEN
         NSUB(LUN) = ISUB-1
      ENDIF

C  READ SUBSET #ISUB FROM MEMORY MESSAGE #IMSG            
C  -------------------------------------------                          

      CALL READSB(MUNIT,IRET)                                           
      IF(IRET.NE.0) GOTO 903                                            

C  RESET SUBSET POINTER BACK TO ZERO AND RETURN
C  --------------------------------------------

      MBYT(LUN) = MBYM                                                  
      NSUB(LUN) = 0                                                     
      RETURN                                                            

C  ERROR EXITS
C  -----------
                                                                        
900   CALL BORT('RDMEMS - MEMORY FILE NOT OPEN FOR INPUT')             
901   CALL BORT('RDMEMS - MEMORY MESSAGE NOT OPEN       ')             
902   CALL BORT('RDMEMS - MEMORY MESSAGE NOT AT BEGINING')             
903   CALL BORT('RDMEMS - ERROR CALLING READSB          ')             
      END                                                               
