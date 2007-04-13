C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE DXINIT(LUN,IOI)                                        
                                                                        
      COMMON /PADESC/ IBCT,IPD1,IPD2,IPD3,IPD4                          
      COMMON /REPTAB/ IDNR(5,2),TYPS(5,2),REPS(5,2),LENS(5)             
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /TABABD/ NTBA(0:32),NTBB(0:32),NTBD(0:32),MTAB(50,32),     
     .                IDNA(50,32,2),IDNB(250,32),IDND(250,32),          
     .                TABA(50,32),TABB(250,32),TABD(250,32)             
                                                                        
      CHARACTER*600 TABD                                                
      CHARACTER*128 TABB                                                
      CHARACTER*128 TABA                                                
      CHARACTER*8   INIB(6,5),INID(5)                                   
      CHARACTER*6   ADN30                                               
      CHARACTER*3   TYPS                                                
      CHARACTER*1   REPS                                                
                                                                        
      DATA INIB   /'------','BYTCNT  ','BYTES  ','+0','+0','16',        
     .             '------','BITPAD  ','NONE   ','+0','+0','1 ',        
     .             '031000','DRF1BIT ','NUMERIC','+0','+0','1 ',        
     .             '031001','DRF8BIT ','NUMERIC','+0','+0','8 ',        
     .             '031002','DRF16BIT','NUMERIC','+0','+0','16'/        
      DATA NINIB  /5/                                                   
                                                                        
      DATA INID   /'        ',                                          
     .             'DRP16BIT',                                          
     .             'DRP8BIT ',                                          
     .             'DRPSTAK ',                                          
     .             'DRP1BIT '/                                          
      DATA NINID  /5/                                                   
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  CLEAR OUT A MESSAGE CONTROL WORD PARTITION                           
C  ------------------------------------------                           
                                                                        
      NMSG(LUN)  = 0                                                    
      NSUB(LUN)  = 0                                                    
      MSUB(LUN)  = 0                                                    
      INODE(LUN) = 0                                                    
      IDATE(LUN) = 0                                                    
                                                                        
C  CLEAR OUT A TABLE PARTITION                                          
C  ---------------------------                                          
                                                                        
      NTBA(LUN) = 0                                                     
      DO I=1,NTBA(0)                                                    
      TABA(I,LUN) = ' '                                                 
      MTAB(I,LUN) = 0                                                   
      ENDDO                                                             
                                                                        
      NTBB(LUN) = 0                                                     
      DO I=1,NTBB(0)                                                    
      TABB(I,LUN) = ' '                                                 
      ENDDO                                                             
                                                                        
      NTBD(LUN) = 0                                                     
      DO I=1,NTBD(0)                                                    
      TABD(I,LUN) = ' '                                                 
      CALL PKTDD(I,LUN,0,IRET)                                          
      ENDDO                                                             
                                                                        
      IF(IOI.EQ.0) RETURN                                               
                                                                        
C  INITIALIZE TABLE WITH APRIORI TABLE B AND D ENTRIES                  
C  ---------------------------------------------------                  
                                                                        
      INIB(1,1) = ADN30(IBCT,6)                                         
      INIB(1,2) = ADN30(IPD4,6)                                         
                                                                        
      DO I=1,NINIB                                                      
      NTBB(LUN) = NTBB(LUN)+1                                           
      IDNB(I,LUN) = IFXY(INIB(1,I))                                     
      TABB(I,LUN)(  1:  6) = INIB(1,I)                                  
      TABB(I,LUN)(  7: 70) = INIB(2,I)                                  
      TABB(I,LUN)( 71: 94) = INIB(3,I)                                  
      TABB(I,LUN)( 95: 98) = INIB(4,I)                                  
      TABB(I,LUN)( 99:109) = INIB(5,I)                                  
      TABB(I,LUN)(110:112) = INIB(6,I)                                  
      ENDDO                                                             
                                                                        
      DO I=2,NINID                                                      
      N = NTBD(LUN)+1                                                   
      IDND(N,LUN) = IDNR(I,1)                                           
      TABD(N,LUN)(1: 6) = ADN30(IDNR(I,1),6)                            
      TABD(N,LUN)(7:70) = INID(I)                                       
      CALL PKTDD(N,LUN,IDNR(1,1),IRET)                                  
      CALL PKTDD(N,LUN,IDNR(I,2),IRET)                                  
      NTBD(LUN) = N                                                     
      ENDDO                                                             
                                                                        
      RETURN                                                            
      END                                                               
