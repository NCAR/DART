      SUBROUTINE CHEKSTAB(LUN)                                          

C************************************************************************
C* CHEKSTAB								*
C*									*
C* This subroutine checks that an internal BUFR tables representation	*
C* is self-consistent and fully defined.  If any errors are found, then	*
C* an appropriate call is made to subroutine BORT.			*
C*									*
C* CHEKSTAB  ( LUN )							*
C*									*
C* Input parameters:							*
C*	LUN		INTEGER		I/O stream index into internal	*
C*					arrays for tables to be checked	*
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		05/01	Added documentation			*
C************************************************************************
                                                                        
      COMMON /TABABD/ NTBA(0:32),NTBB(0:32),NTBD(0:32),MTAB(50,32),     
     .                IDNA(50,32,2),IDNB(250,32),IDND(250,32),          
     .                TABA(50,32),TABB(250,32),TABD(250,32)             
      CHARACTER*600 TABD                                                
      CHARACTER*128 TABB                                                
      CHARACTER*128 TABA                                                
      CHARACTER*24  UNIT                                                
      CHARACTER*8   NEMO,NEMS(250)                                      
      CHARACTER*1   TAB                                                 
      DIMENSION     IRPS(250),KNTS(250)                                 
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  THERE MUST BE ENTRIES IN TABLES A, B, AND D                          
C  -------------------------------------------                          
                                                                        
      IF(NTBA(LUN).EQ.0) GOTO 900                                       
      IF(NTBB(LUN).EQ.0) GOTO 901                                       
      IF(NTBD(LUN).EQ.0) GOTO 902                                       
                                                                        
C  MAKE SURE EACH TABLE A ENTRY DEFINED AS A SEQUENCE                   
C  --------------------------------------------------                   
                                                                        
      DO I=1,NTBA(LUN)                                                  
      NEMO = TABA(I,LUN)(4:11)                                          
      CALL NEMTAB(LUN,NEMO,IDN,TAB,IRET)                                
      IF(TAB.NE.'D') GOTO 903                                           
      ENDDO                                                             
                                                                        
C  CHECK TABLE B CONTENTS                                               
C  ----------------------                                               
                                                                        
      DO ITAB=1,NTBB(LUN)                                               
      CALL NEMTBB(LUN,ITAB,UNIT,ISCL,IREF,IBIT)                         
      ENDDO                                                             
                                                                        
C  CHECK TABLE D CONTNETS                                               
C  ----------------------                                               
                                                                        
      DO ITAB=1,NTBD(LUN)                                               
      CALL NEMTBD(LUN,ITAB,NSEQ,NEMS,IRPS,KNTS)                         
      ENDDO                                                             
                                                                        
      RETURN                                                            
900   CALL BORT('CHEKSTAB - EMPTY TABLE A')                            
901   CALL BORT('CHEKSTAB - EMPTY TABLE B')                            
902   CALL BORT('CHEKSTAB - EMPTY TABLE D')                            
903   CALL BORT('CHEKSTAB - NO SEQUENCE DEFINED FOR TABLE A: '//NEMO)  
      END                                                               
