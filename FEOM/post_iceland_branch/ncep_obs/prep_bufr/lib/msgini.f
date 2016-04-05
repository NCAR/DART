      SUBROUTINE MSGINI(LUN)                                            

C************************************************************************
C* MSGINI								*
C*									*
C* This subroutine initializes, within the internal arrays, a new BUFR	*
C* message for output.							*
C*									*
C* MSGINI  ( LUN )							*
C*									*
C* Input parameters:							*
C*	LUN		INTEGER		I/O stream index into internal	*
C*					arrays for new BUFR message	*
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		05/01	Added documentation			*
C************************************************************************
                                                                        
      COMMON /PADESC/ IBCT,IPD1,IPD2,IPD3,IPD4                          
      COMMON /MSGPTR/ NBY0,NBY1,NBY2,NBY3,NBY4,NBY5                     
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)     
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
                                                                        
      CHARACTER*10 TAG                                                  
      CHARACTER*8  SUBTAG                                               
      CHARACTER*4  BUFR,SEVN                                            
      CHARACTER*3  TYP                                                  
      CHARACTER*1  TAB                                                  
                                                                        
      DATA BUFR/'BUFR'/                                                 
      DATA SEVN/'7777'/                                                 
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  GET THE MESSAGE TAG AND TYPE, AND BREAK UP THE DATE                  
C  ---------------------------------------------------                  
                                                                        
      SUBTAG = TAG(INODE(LUN))                                          
      CALL NEMTBA(LUN,SUBTAG,MTYP,MSBT,INOD)                            
      IF(INODE(LUN).NE.INOD) GOTO 900                                   
      CALL NEMTAB(LUN,SUBTAG,ISUB,TAB,IRET)                             
      IF(IRET.EQ.0) GOTO 901                                            
                                                                        
C  DATE CAN BE YYMMDDHH OR YYYYMMDDHH
C  ----------------------------------                                   
                                                                        
      MCEN = MOD(IDATE(LUN)/10**8,100)+1
      MEAR = MOD(IDATE(LUN)/10**6,100)                             
      MMON = MOD(IDATE(LUN)/10**4,100)                             
      MDAY = MOD(IDATE(LUN)/10**2,100)                             
      MOUR = MOD(IDATE(LUN)      ,100)                             
      MMIN = 0                                                       

      IF(MCEN.EQ.1) CALL BORT('MSGINI - IDATE(LUN) = 0!')
      IF(MEAR.EQ.0) MCEN = MCEN-1
      IF(MEAR.EQ.0) MEAR = 100
                                                                        
C  INITIALIZE THE MESSAGE                                               
C  ----------------------                                               
                                                                        
      MBIT = 0                                                          
      NBY0 = 8                                                          
      NBY1 = 18                                                         
      NBY2 = 0                                                          
      NBY3 = 20                                                         
      NBY4 = 4                                                          
      NBY5 = 4                                                          
      NBYT = NBY0+NBY1+NBY2+NBY3+NBY4+NBY5                              
                                                                        
C  SECTION 0                                                            
C  ---------                                                            
                                                                        
      CALL PKC(BUFR ,  4 , MBAY(1,LUN),MBIT)                            
      CALL PKB(NBYT , 24 , MBAY(1,LUN),MBIT)                            
      CALL PKB(   3 ,  8 , MBAY(1,LUN),MBIT)                            
                                                                        
C  SECTION 1                                                            
C  ---------                                                            
                                                                        
      CALL PKB(NBY1 , 24 , MBAY(1,LUN),MBIT)                            
      CALL PKB(   0 ,  8 , MBAY(1,LUN),MBIT)                            
      CALL PKB(   3 ,  8 , MBAY(1,LUN),MBIT)                            
      CALL PKB(   7 ,  8 , MBAY(1,LUN),MBIT)                            
      CALL PKB(   0 ,  8 , MBAY(1,LUN),MBIT)                            
      CALL PKB(   0 ,  8 , MBAY(1,LUN),MBIT)                            
      CALL PKB(MTYP ,  8 , MBAY(1,LUN),MBIT)                            
      CALL PKB(MSBT ,  8 , MBAY(1,LUN),MBIT)                            
      CALL PKB(   4 ,  8 , MBAY(1,LUN),MBIT)                            
      CALL PKB(   0 ,  8 , MBAY(1,LUN),MBIT)                            
      CALL PKB(MEAR ,  8 , MBAY(1,LUN),MBIT)                            
      CALL PKB(MMON ,  8 , MBAY(1,LUN),MBIT)                            
      CALL PKB(MDAY ,  8 , MBAY(1,LUN),MBIT)                            
      CALL PKB(MOUR ,  8 , MBAY(1,LUN),MBIT)                            
      CALL PKB(MMIN ,  8 , MBAY(1,LUN),MBIT)                            
      CALL PKB(MCEN ,  8 , MBAY(1,LUN),MBIT)                            
                                                                        
C  SECTION 3                                                            
C  ---------                                                            
                                                                        
      CALL PKB(NBY3 , 24 , MBAY(1,LUN),MBIT)                            
      CALL PKB(   0 ,  8 , MBAY(1,LUN),MBIT)                            
      CALL PKB(   0 , 16 , MBAY(1,LUN),MBIT)                            
      CALL PKB(2**7 ,  8 , MBAY(1,LUN),MBIT)                            
      CALL PKB(IBCT , 16 , MBAY(1,LUN),MBIT)                            
      CALL PKB(ISUB , 16 , MBAY(1,LUN),MBIT)                            
      CALL PKB(IPD1 , 16 , MBAY(1,LUN),MBIT)                            
      CALL PKB(IPD2 , 16 , MBAY(1,LUN),MBIT)                            
      CALL PKB(IPD3 , 16 , MBAY(1,LUN),MBIT)                            
      CALL PKB(IPD4 , 16 , MBAY(1,LUN),MBIT)                            
      CALL PKB(   0 ,  8 , MBAY(1,LUN),MBIT)                            
                                                                        
C  SECTION 4                                                            
C  ---------                                                            
                                                                        
      CALL PKB(NBY4 , 24 , MBAY(1,LUN),MBIT)                            
      CALL PKB(   0 ,  8 , MBAY(1,LUN),MBIT)                            
                                                                        
C  SECTION 5                                                            
C  ---------                                                            
                                                                        
      CALL PKC(SEVN ,  4 , MBAY(1,LUN),MBIT)                            
                                                                        
C  DOUBLE CHECK INITIAL MESSAGE LENGTH                                  
C  -----------------------------------                                  
                                                                        
      IF(MOD(MBIT,8).NE.0) GOTO 902                                     
      IF(MBIT/8.NE.NBYT  ) GOTO 903                                     
                                                                        
      NMSG(LUN) = NMSG(LUN)+1                                           
      NSUB(LUN) = 0                                                     
      MBYT(LUN) = NBYT                                                  
                                                                        
      RETURN                                                            

C  ENTRY MINIMG WILL WRITE MINUTES INTO A MESSAGE SECTION ONE 
C  ----------------------------------------------------------

      ENTRY MINIMG(LUNIT,MINI)
C************************************************************************
C* MINIMG								*
C*									*
C* This entry point should only be called when logical unit LUNIT has	*
C* been opened for output operations.  It packs the value of MINI into	*
C* byte 17 of Section 1 of the BUFR message that is currently open	*
C* within memory for LUNIT, so that this value then becomes the minutes	*
C* component of the Section 1 date-time for the message.		*
C*									*
C* MINIMG  ( LUNIT, MINI )						*
C*									*
C* Input parameters:							*
C*	LUNIT		INTEGER		FORTRAN logical unit number	*
C*	MINI		INTEGER		Minutes value to be packed	*
C************************************************************************

      CALL STATUS(LUNIT,LUT,IL,IM)                                      
      IF(IL.LE.0) GOTO 904                                              
      IF(IM.EQ.0) GOTO 905                                              
      MBIT = 24*8
      CALL PKB(MINI,8,MBAY(1,LUT),MBIT)                            

      RETURN

C  ERROR EXITS
C  -----------

900   CALL BORT('MSGINI - INODE MISMATCH                      ')       
901   CALL BORT('MSGINI - CANT FIND SUBSET IDN: '//SUBTAG      )       
902   CALL BORT('MSGINI - INITIAL MESSAGE OUT OF BYTE BOUNDARY')       
903   CALL BORT('MSGINI - INITIAL MESSAGE COUNT FAILS CHECK   ')       
904   CALL BORT('MINIMG - FILE IS NOT OPEN FOR OUTPUT         ')       
905   CALL BORT('MINIMG - NO MESSGAGE OPEN                    ')       
      END                                                               
