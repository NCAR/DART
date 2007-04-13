C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE DXMINI(LUN,MBAY,MBYT,MB4,MBA,MBB,MBD)                  
                                                                        
      COMMON /TABABD/ NTBA(0:32),NTBB(0:32),NTBD(0:32),MTAB(50,32),     
     .                IDNA(50,32,2),IDNB(250,32),IDND(250,32),          
     .                TABA(50,32),TABB(250,32),TABD(250,32)             
      COMMON /DXTAB / MAXDX,IDXV,NXSTR(10),LDXA(10),LDXB(10),LDXD(10),  
     .                LD30(10),DXSTR(10)                                
                                                                        
      CHARACTER*600 TABD                                                
      CHARACTER*128 TABB                                                
      CHARACTER*128 TABA                                                
      CHARACTER*56  DXSTR                                               
      DIMENSION     MBAY(5000)                                          
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C  ENTRY POINTS FOR SEPARATING TABLES A B D
C  ----------------------------------------

      MSBT = IDXV                                                       
      GOTO 1
      ENTRY DXMINA(LUN,MBAY,MBYT,MB4,MBA,MBB,MBD)                  
      MSBT = 1                                                          
      GOTO 1
      ENTRY DXMINB(LUN,MBAY,MBYT,MB4,MBA,MBB,MBD)                  
      MSBT = 2                                                          
      GOTO 1
      ENTRY DXMIND(LUN,MBAY,MBYT,MB4,MBA,MBB,MBD)                  
      MSBT = 4                                                          
1     CONTINUE
      
                                                                        
C  INITIALIZE THE MESSAGE                                               
C  ----------------------                                               
                                                                        
      MBIT = 0                                                          
      DO I=1,5000                                                       
      MBAY(I) = 0                                                       
      ENDDO                                                             
                                                                        
      IH   = 0                                                          
      ID   = 0                                                          
      IM   = 0                                                          
      IY   = 0                                                          
      MTYP = 11                                                         
      NSUB = 1                                                          
      IDXS = IDXV+1                                                     
      LDXS = NXSTR(IDXS)                                                
                                                                        
      NBY0 = 8                                                          
      NBY1 = 18                                                         
      NBY2 = 0                                                          
      NBY3 = 7 + NXSTR(IDXS) + 1                                        
      NBY4 = 7                                                          
      NBY5 = 4                                                          
      MBYT = NBY0+NBY1+NBY2+NBY3+NBY4+NBY5                              
                                                                        
      IF(MOD(NBY3,2).NE.0) GOTO 900                                     
                                                                        
C  SECTION 0                                                            
C  ---------                                                            
                                                                        
      CALL PKC('BUFR' ,  4 , MBAY,MBIT)                                 
      CALL PKB(  MBYT , 24 , MBAY,MBIT)                                 
      CALL PKB(     3 ,  8 , MBAY,MBIT)                                 
                                                                        
C  SECTION 1                                                            
C  ---------                                                            
                                                                        
      CALL PKB(  NBY1 , 24 , MBAY,MBIT)                                 
      CALL PKB(     0 ,  8 , MBAY,MBIT)                                 
      CALL PKB(     3 ,  8 , MBAY,MBIT)                                 
      CALL PKB(     7 ,  8 , MBAY,MBIT)                                 
      CALL PKB(     0 ,  8 , MBAY,MBIT)                                 
      CALL PKB(     0 ,  8 , MBAY,MBIT)                                 
      CALL PKB(  MTYP ,  8 , MBAY,MBIT)                                 
      CALL PKB(  MSBT ,  8 , MBAY,MBIT)                                 
      CALL PKB(     4 ,  8 , MBAY,MBIT)                                 
      CALL PKB(  IDXV ,  8 , MBAY,MBIT)                                 
      CALL PKB(    IY ,  8 , MBAY,MBIT)                                 
      CALL PKB(    IM ,  8 , MBAY,MBIT)                                 
      CALL PKB(    ID ,  8 , MBAY,MBIT)                                 
      CALL PKB(    IH ,  8 , MBAY,MBIT)                                 
      CALL PKB(     0 ,  8 , MBAY,MBIT)                                 
      CALL PKB(     0 ,  8 , MBAY,MBIT)                                 
                                                                        
C  SECTION 3                                                            
C  ---------                                                            
                                                                        
      CALL PKB(       NBY3 ,   24 , MBAY,MBIT)                          
      CALL PKB(          0 ,    8 , MBAY,MBIT)                          
      CALL PKB(          1 ,   16 , MBAY,MBIT)                          
      CALL PKB(       2**7 ,    8 , MBAY,MBIT)                          
      DO I=1,LDXS                                                       
      CALL PKB(IUPM(DXSTR(IDXS)(I:I),8),8,MBAY,MBIT)                    
      ENDDO                                                             
      CALL PKB(          0 ,    8 , MBAY,MBIT)                          
                                                                        
C  SECTION 4                                                            
C  ---------                                                            
                                                                        
      MB4 = MBIT/8+1                                                    
      CALL PKB(NBY4 , 24 , MBAY,MBIT)                                   
      CALL PKB(   0 ,  8 , MBAY,MBIT)                                   
      MBA = MBIT/8+1                                                    
      CALL PKB(   0 ,  8 , MBAY,MBIT)                                   
      MBB = MBIT/8+1                                                    
      CALL PKB(   0 ,  8 , MBAY,MBIT)                                   
      MBD = MBIT/8+1                                                    
      CALL PKB(   0 ,  8 , MBAY,MBIT)                                   
                                                                        
      IF(MBIT/8+NBY5.NE.MBYT) GOTO 901                                  
                                                                        
      RETURN                                                            
900   CALL BORT('DXMINI - UNEVEN SECTION 3')                           
901   CALL BORT('DXMINI - BYTCNT IS OFF   ')                           
      END                                                               
