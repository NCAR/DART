      SUBROUTINE MAKESTAB                                               

C************************************************************************
C* MAKESTAB								*
C*									*
C* Using the information within the internal BUFR table arrays (within	*
C* COMMON /TABABD/) for *all* of the LUN (i.e. I/O stream index) values	*
C* that are currently defined to the BUFRLIB software, this subroutine	*
C* constructs an internal jump/link table within COMMON /TABLES/.	*
C* Note that the entire jump/link table will always be completely	*
C* reconstructed from scratch, even if some of the information within	*
C* the internal BUFR table arrays already existed there at the time of	*
C* the previous call to this subroutine, because there may have been	*
C* other events that have taken place since the previous call to this	*
C* subroutine that have *not* yet been reflected within the internal	*
C* jump/link table, such as, e.g. the unlinking of an LUN value from	*
C* the BUFRLIB software via a call to subroutine CLOSBF.		*
C*									*
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		05/01	Added documentation			*
C************************************************************************
                                                                        
      COMMON /QUIET/  IPRT                                              
      COMMON /STBFR/  IOLUN(32),IOMSG(32)                               
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /USRINT/ NVAL(32),INV(15000,32),VAL(15000,32)              
      COMMON /TABABD/ NTBA(0:32),NTBB(0:32),NTBD(0:32),MTAB(50,32),     
     .                IDNA(50,32,2),IDNB(250,32),IDND(250,32),          
     .                TABA(50,32),TABB(250,32),TABD(250,32)             
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
                                                                        
      CHARACTER*600 TABD                                                
      CHARACTER*128 TABB                                                
      CHARACTER*128 TABA                                                
      CHARACTER*10  TAG                                                 
      CHARACTER*8   NEMO                                                
      CHARACTER*3   TYP                                                 
      DIMENSION     LUS(32)                                             
      LOGICAL       EXPAND,PRTTAB                                       
      REAL*8        VAL,BMISS                                           
                                                                        
      DATA PRTTAB /.FALSE./                                            
      DATA BMISS  /10E10/                                               
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      PRTTAB = IPRT.GE.2                                                
                                                                        
C  RESET POINTER TABLE AND STRING CACHE                                 
C  ------------------------------------                                 
                                                                        
      NTAB = 0                                                          
      CALL STRCLN                                                       
                                                                        
C  FIGURE OUT WHICH UNITS SHARE TABLES                                  
C  -----------------------------------                                  
                                                                        
C     First, determine how many LUN values are currently being used and,
C     for each such one, whether it uses the same DX table information as
C     any other LUN values that we have examined so far.  If so, then set
C     LUS(LUN) to a nonzero value.

C     Note that, for each LUN value, the MTAB(*,LUN) array contains
C     pointer indices into the internal jump/link table for each of the
C     Table A mnemonics that is currently defined for that LUN value.
C     Thus, the code within the following DO loop is simply checking
C     whether the first Table A mnemonic is the same for two different
C     LUN values as the determination of whether those LUN values indeed
C     share the same DX tables.

      DO LUN=1,32                                                       
      LUS(LUN) = 0                                                      
      IF(IOLUN(LUN).NE.0) THEN                                          
         DO LUM=1,LUN-1                                                 
         IF(MTAB(1,LUN).EQ.MTAB(1,LUM)) LUS(LUN) = LUM                  
         ENDDO                                                          
      ENDIF                                                             
      ENDDO                                                             
                                                                        
C  INITIALIZE JUMP-LINK TABLES WITH SUBSETS/SEQUENCES/ELEMENTS          
C  ----------------------------------------------------------           
                                                                        
      DO LUN=1,32                                                       
                                                                        
      IF(IOLUN(LUN).NE.0) THEN                                          
                                                                        
C  RESET ANY EXISTING INVENTORY POINTERS                                
C  -------------------------------------                                
                                                                        
         IF(IOMSG(LUN).NE.0) THEN                                       
            IF(LUS(LUN).EQ.0) INC = (NTAB+1)-MTAB(1,LUN)                
            IF(LUS(LUN).NE.0) INC = MTAB(1,LUS(LUN))-MTAB(1,LUN)        
            DO N=1,NVAL(LUN)                                            
            INV(N,LUN) = INV(N,LUN)+INC                                 
            ENDDO                                                       
         ENDIF                                                          
                                                                        
C  CREATE NEW TABLE ENTRIES IF THIS UNIT DOESN'T SHARE EXISTING ONES    
C  -----------------------------------------------------------------    
                                                                        
         IF(LUS(LUN).EQ.0) THEN                                         

C	    The DX table information corresponding to this LUN has not
C	    yet been written into the internal jump/link table, so add
C	    it in now.

            CALL CHEKSTAB(LUN)                                          
            DO ITBA=1,NTBA(LUN)                                         
            INOD = NTAB+1                                               
            NEMO = TABA(ITBA,LUN)(4:11)                                 
            CALL TABSUB(LUN,NEMO)                                       
            MTAB(ITBA,LUN) = INOD                                       
            ISC(INOD)      = NTAB                                       

C**** note that the following lines are commented out****
C           DO N1=INOD,ISC(INOD)-1                                      
C           DO N2=N1+1,ISC(INOD)                                        
C           IF(TAG(N1).EQ.TAG(N2)) GOTO 900                             
C           ENDDO                                                       
C           ENDDO                                                       
C********************************************************

            ENDDO                                                       
         ENDIF                                                          
                                                                        
      ENDIF                                                             
      ENDDO                                                             
                                                                        
C  STORE TYPES AND INITIAL VALUES AND COUNTS                            
C  -----------------------------------------                            
                                                                        
      DO NODE=1,NTAB                                                    
      IF(TYP(NODE).EQ.'SUB') THEN                                       
         VALI(NODE) = 0                                                 
         KNTI(NODE) = 1                                                 
         ITP (NODE) = 0                                                 
      ELSEIF(TYP(NODE).EQ.'SEQ') THEN                                   
         VALI(NODE) = 0                                                 
         KNTI(NODE) = 1                                                 
         ITP (NODE) = 0                                                 
      ELSEIF(TYP(NODE).EQ.'RPC') THEN                                   
         VALI(NODE) = 0                                                 
         KNTI(NODE) = 0                                                 
         ITP (NODE) = 0                                                 
      ELSEIF(TYP(NODE).EQ.'RPS') THEN                                   
         VALI(NODE) = 0                                                 
         KNTI(NODE) = 0                                                 
         ITP (NODE) = 0                                                 
      ELSEIF(TYP(NODE).EQ.'REP') THEN                                   
         VALI(NODE) = BMISS                                             
         KNTI(NODE) = IRF(NODE)                                         
         ITP (NODE) = 0                                                 
      ELSEIF(TYP(NODE).EQ.'DRS') THEN                                   
         VALI(NODE) = 0                                                 
         KNTI(NODE) = 1                                                 
         ITP (NODE) = 1                                                 
      ELSEIF(TYP(NODE).EQ.'DRP') THEN                                   
         VALI(NODE) = 0                                                 
         KNTI(NODE) = 1                                                 
         ITP (NODE) = 1                                                 
      ELSEIF(TYP(NODE).EQ.'DRB') THEN                                   
         VALI(NODE) = 0                                                 
         KNTI(NODE) = 0                                                 
         ITP (NODE) = 1                                                 
      ELSEIF(TYP(NODE).EQ.'NUM') THEN                                   
         VALI(NODE) = BMISS                                             
         KNTI(NODE) = 1                                                 
         ITP (NODE) = 2                                                 
      ELSEIF(TYP(NODE).EQ.'CHR') THEN                                   
         VALI(NODE) = BMISS                                             
         KNTI(NODE) = 1                                                 
         ITP (NODE) = 3                                                 
      ELSE                                                              
         GOTO 901                                                       
      ENDIF                                                             
      ENDDO                                                             
                                                                        
C  SET UP EXPANSION SEGMENTS FOR TYPE 'SUB', 'DRP', AND 'DRS' NODES     
C  ----------------------------------------------------------------     
                                                                        
      NEWN = 0                                                          
                                                                        
      DO N=1,NTAB                                                       
      ISEQ(N,1) = 0                                                     
      ISEQ(N,2) = 0                                                     
      EXPAND = TYP(N).EQ.'SUB' .OR. TYP(N).EQ.'DRP' .OR. TYP(N).EQ.'DRS'
     .                         .OR. TYP(N).EQ.'REP' .OR. TYP(N).EQ.'DRB'
      IF(EXPAND) THEN                                                   
         ISEQ(N,1) = NEWN+1                                             
         NODA = N                                                       
         NODE = N+1                                                     
         DO K=1,15000                                                   
         KNT(K) = 0                                                     
         ENDDO                                                          
         IF(TYP(NODA).EQ.'REP') KNT(NODE) = KNTI(NODA)                  
         IF(TYP(NODA).NE.'REP') KNT(NODE) = 1                           
                                                                        
1        NEWN = NEWN+1                                                  
         IF(NEWN.GT.15000) GOTO 902                                     
         JSEQ(NEWN) = NODE                                              
         KNT(NODE) = MAX(KNTI(NODE),KNT(NODE))                          
2        IF(JUMP(NODE)*KNT(NODE).GT.0) THEN                             
            NODE = JUMP(NODE)                                           
            GOTO 1                                                      
         ELSE IF(LINK(NODE).GT.0) THEN                                  
            NODE = LINK(NODE)                                           
            GOTO 1                                                      
         ELSE                                                           
            NODE = JMPB(NODE)                                           
            IF(NODE.EQ.NODA) GOTO 3                                     
            IF(NODE.EQ.0   ) GOTO 903                                   
            KNT(NODE) = MAX(KNT(NODE)-1,0)                              
            GOTO 2                                                      
         ENDIF                                                          
3        ISEQ(N,2) = NEWN                                               
      ENDIF                                                             
      ENDDO                                                             
                                                                        
C  PRINT THE SEQUENCE TABLES                                            
C  ------------------------                                             
                                                                        
      IF(PRTTAB) THEN                                                   
         PRINT*                                                         
         DO I=1,NTAB                                                    
         PRINT99,I,                                                     
     .   TAG(I),TYP(I),JMPB(I),JUMP(I),LINK(I),IBT(I),IRF(I),ISC(I)     
         ENDDO                                                          
         PRINT*                                                         
99       FORMAT(I5,2X,A10,A5,6I8)                                       
      ENDIF                                                             
                                                                        
      RETURN                                                            
900   CALL BORT('MAKESTAB - DUP IN SUBSET: '//TAG(N1)//':'//NEMO)      
901   CALL BORT('MAKESTAB - UNKNOWN TYPE : '         //TYP(NODE))      
902   CALL BORT('MAKESTAB - JSEQ OVERFLOW       : '  //TAG(N   ))      
903   CALL BORT('MAKESTAB - FAILED TO CIRCULATE : '  //TAG(N   ))      
      END                                                               
