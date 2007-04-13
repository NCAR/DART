C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      SUBROUTINE UFBDMP(LUNIT,LUPRT)                                    
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
                                                                        
      CHARACTER*10 TAG,TG                                               
      CHARACTER*8  VC                                                   
      CHARACTER*3  TYP,TP                                               
      CHARACTER*1  YOU                                                  
      EQUIVALENCE  (VL,VC)                                              
      REAL*8       VAL,VL,BMISS
                                                                        
      DATA BMISS /10E10/                                                
                                                                        
C---------------------------------------------------------------------- 
CFFP$ EXPAND (STATUS)                                                   
C---------------------------------------------------------------------- 
                                                                        
      if(luprt.eq.0) luout = 6                                          
      if(luprt.ne.0) luout = luprt                                      
                                                                        
C  CHECK THE FILE STATUS AND I-NODE                                     
C  --------------------------------                                     
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IM.EQ.0) GOTO 901                                              
      IF(INODE(LUN).NE.INV(1,LUN)) GOTO 902                             
                                                                        
C  DUMP THE CONTENTS OF COMMON /USRINT/ FOR UNIT LUNIT                  
C  ---------------------------------------------------                  
                                                                        
      DO NV=1,NVAL(LUN)                                                 
      if(luprt.eq.0 .and. mod(nv,20).eq.0) then                         
         print*,'(MORE)'                                                
         read(5,'(a1)') you                                             
         if(you.eq.'q') return                                          
      endif                                                             
      ND = INV (NV,LUN)                                                 
      VL = VAL (NV,LUN)                                                 
      TG = TAG (ND)                                                     
      TP = TYP (ND)                                                     
      IT = ITP (ND)                                                     
      IB = IBT (ND)                                                     
      IS = ISC (ND)                                                     
      IR = IRF (ND)                                                     
      JP = JUMP(ND)                                                     
      LK = LINK(ND)                                                     
      JB = JMPB(ND)                                                     
      RJ = RJUST(TG)                                                    
      IF(TP.NE.'CHR' .AND. VL.LT.BMISS) THEN                      
         WRITE(LUOUT,1) NV,TP,IT,TG,VL,IB,IS,IR,ND,JP,LK,JB             
      ELSE                                                              
         IF(VL.EQ.BMISS) THEN 
            VC = 'MISSING'                                 
         ELSE
            VC = VC(1:IB/8)                                
         ENDIF
         RJ = RJUST(VC)                                                 
         WRITE(LUOUT,2) NV,TP,IT,TG,VC,IB,IS,IR,ND,JP,LK,JB             
      ENDIF                                                             
      ENDDO                                                             
                                                                        
1     FORMAT(I5,1X,A3,'-',I1,1X,A10,1X,F10.1,7(1X,I5))                  
2     FORMAT(I5,1X,A3,'-',I1,1X,A10,1X,A10  ,7(1X,I5))                  
                                                                        
                                                                        
C  EXITS                                                                
C  -----                                                                
                                                                        
      RETURN                                                            
900   CALL BORT('UFBDMP - FILE IS CLOSED                     ')        
901   CALL BORT('UFBDMP - NO MESSAGE OPEN                    ')        
902   CALL BORT('UFBDMP - I-NODE MISMATCH                    ')        
      END                                                               
