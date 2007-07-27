C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE CPYMEM(LUNOT)                                          
                                                                        
      PARAMETER (MAXMSG=50000,MAXMEM=16000000)                           
                                                                        
      COMMON /MSGMEM/ MUNIT,MLAST,MSGP(0:MAXMSG),MSGS(MAXMEM)           
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)     
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
                                                                        
      CHARACTER*10 TAG                                                  
      CHARACTER*8  SUBSET                                               
      CHARACTER*3  TYP                                                  
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  CHECK THE FILE STATUSES                                              
C  -----------------------                                              
                                                                        
      CALL STATUS(MUNIT,LIN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IL.GT.0) GOTO 901                                              
      IF(IM.LE.0) GOTO 902                                              
                                                                        
      CALL STATUS(LUNOT,LOT,IL,IM)                                      
      IF(IL.EQ.0) GOTO 903                                              
      IF(IL.LT.0) GOTO 904                                              
      IF(IM.NE.0) GOTO 905                                              
                                                                        
C  MAKE SURE BOTH FILES HAVE THE SAME TABLES                            
C  -----------------------------------------                            
                                                                        
      SUBSET = TAG(INODE(LIN))                                          
      CALL NEMTBA(LOT,SUBSET,MTYP,MSBT,INOD)                            
      IF(INODE(LIN).NE.INOD) GOTO 906                                   
                                                                        
C  EVERYTHING OKAY, COPY A MESSAGE                                      
C  -------------------------------                                      
                                                                        
      MBYM = IUPB(MBAY(1,LIN),5,24)                                     
      CALL MSGWRT(LUNOT,MBAY(1,LIN),MBYM)                               
                                                                        
C  SET THE MESSAGE CONTROL WORLDS FOR LUNOT                             
C  ----------------------------------------                             
                                                                        
      NMSG (LOT) = NMSG(LOT) + 1                                        
      NSUB (LOT) = MSUB(LIN)                                            
      MSUB (LOT) = MSUB(LIN)                                            
      IDATE(LOT) = IDATE(LIN)                                           
      INODE(LOT) = INODE(LIN)                                           
                                                                        
C  NORMAL EXIT                                                          
C  -----------                                                          
                                                                        
      RETURN                                                            
                                                                        
C  ERROR EXITS                                                          
C  -----------                                                          
                                                                        
900   CALL BORT('CPYMEM - INPUT FILE IS CLOSED                    ')   
901   CALL BORT('CPYMEM - INPUT FILE IS OPEN FOR OUTPUT           ')   
902   CALL BORT('CPYMEM - NO INPUT FILE MESSAGE OPEN              ')   
903   CALL BORT('CPYMEM - OUTPUT FILE IS CLOSED                   ')   
904   CALL BORT('CPYMEM - OUTPUT FILE IS OPEN FOR OUTPUT          ')   
905   CALL BORT('CPYMEM - OUTPUT FILE MESSAGE OPEN                ')   
906   CALL BORT('CPYMEM - INPUT/OUTPUT FILES HAVE DIFFERENT TABLES')   
      END                                                               
