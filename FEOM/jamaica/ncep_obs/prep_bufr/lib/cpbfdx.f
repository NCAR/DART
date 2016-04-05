C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE CPBFDX(LUD,LUN)                                        
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /TABABD/ NTBA(0:32),NTBB(0:32),NTBD(0:32),MTAB(50,32),     
     .                IDNA(50,32,2),IDNB(250,32),IDND(250,32),          
     .                TABA(50,32),TABB(250,32),TABD(250,32)             
                                                                        
      CHARACTER*600 TABD                                                
      CHARACTER*128 TABB                                                
      CHARACTER*128 TABA                                                
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  INITIALIZE THE DX-TABLE PARTITION                                    
C  ---------------------------------                                    
                                                                        
      CALL DXINIT(LUN,0)                                                
                                                                        
C  COPY ONE TABLE PARTITION TO ANOTHER                                  
C  -----------------------------------                                  
                                                                        
      INODE(LUN) = INODE(LUD)                                           
                                                                        
      NTBA(LUN) = NTBA(LUD)                                             
      NTBB(LUN) = NTBB(LUD)                                             
      NTBD(LUN) = NTBD(LUD)                                             
                                                                        
      DO I=1,NTBA(LUD)                                                  
      IDNA(I,LUN,1) = IDNA(I,LUD,1)                                     
      IDNA(I,LUN,2) = IDNA(I,LUD,2)                                     
      TABA(I,LUN) = TABA(I,LUD)                                         
      MTAB(I,LUN) = MTAB(I,LUD)                                         
      ENDDO                                                             
                                                                        
      DO I=1,NTBB(LUD)                                                  
      IDNB(I,LUN) = IDNB(I,LUD)                                         
      TABB(I,LUN) = TABB(I,LUD)                                         
      ENDDO                                                             
                                                                        
      DO I=1,NTBD(LUD)                                                  
      IDND(I,LUN) = IDND(I,LUD)                                         
      TABD(I,LUN) = TABD(I,LUD)                                         
      ENDDO                                                             
                                                                        
      RETURN                                                            
      END                                                               
