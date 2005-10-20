      SUBROUTINE SEQSDX(CARD,LUN)                                       

C************************************************************************
C* SEQSDX								*
C*									*
C* This subroutine decodes the Table D sequence information from a	*
C* mnemonic definition card that was previously read from a user DX	*
C* table by subroutine RDUSDX.  This decoded information is then added	*
C* to the already-existing entry for that mnemonic within the internal	*
C* BUFR table D array TABD(*,LUN).					*
C*									*
C* SEQSDX  ( CARD, LUN )						*
C*									*
C* Input parameters:							*
C*	CARD		CHARACTER*80	Mnemonic definition card that	*
C*					was read from a user DX table	*
C*	LUN		INTEGER		I/O stream index into internal	*
C*					arrays for this user DX table	*
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		05/01	Added documentation			*
C************************************************************************
                                                                        
      COMMON /REPTAB/ IDNR(5,2),TYPS(5,2),REPS(5,2),LENS(5)             
                                                                        
      CHARACTER*80  CARD,SEQS                                           
      CHARACTER*12  ATAG,TAGS(250)                                      
      CHARACTER*8   NEMO,NEMA,NEMB                                      
      CHARACTER*3   TYPS                                                
      CHARACTER*1   REPS,TAB                                            
                                                                        
      DATA MAXTGS /250/                                                 
      DATA MAXTAG /12/                                                  
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  FIND THE SEQUENCE TAG IN TABLE D AND PARSE THE SEQUENCE STRING       
C  --------------------------------------------------------------       
                                                                        
      NEMO = CARD( 3:10)                                                
      SEQS = CARD(14:78)                                                

C     Note that an entry for this mnemonic should already exist within
C     the internal BUFR table D array TABD(*,LUN); this entry should
C     have been created by subroutine RDUSDX when the mnemonic and its
C     associated FXY value and description were initially defined within
C     a card read from near the top of the user DX table.  Now, we need
C     to retrieve the positional index for that entry within TABD(*,LUN)
C     so that we can access the entry and then add the decoded sequence
C     information to it.
                                                                        
      CALL NEMTAB(LUN,NEMO,IDN,TAB,ISEQ)                                
      CALL PARSEQ(SEQS,TAGS,MAXTGS,NTAG)                                
      IF(TAB.NE.'D') GOTO 900                                           
      IF(NTAG.EQ.0 ) GOTO 900                                           
                                                                        
      DO N=1,NTAG                                                       
      ATAG = TAGS(N)                                                    
      IREP = 0                                                          
                                                                        
C  CHECK FOR REPLICATOR                                                 
C  --------------------                                                 
                                                                        
      DO I=1,5                                                          
      IF(ATAG(1:1).EQ.REPS(I,1)) THEN                                   

C	 Note that REPS(*,*), which contains all of the symbols used to
C	 denote all of the various replication schemes that are possible
C	 within a user DX table, was previously defined within
C	 subroutine BFRINI.

         DO J=2,MAXTAG                                                  
         IF(ATAG(J:J).EQ.REPS(I,2)) THEN                                
            IF(J.EQ.MAXTAG) GOTO 901                                    

C	    Note that subroutine STRNUM will return NUMR = 0 if the
C	    string passed to it contains all blanks (as *should* be the
C	    case whenever I = 2, 3, 4, or 5).

C	    However, when I = 1, then subroutine STRNUM will return
C	    NUMR = (the number of replications for the mnemonic using
C	    F=1 "regular" (i.e. non-delayed) replication).

            CALL STRNUM(ATAG(J+1:MAXTAG),NUMR)                          
            IF(I.EQ.1 .AND. NUMR.LE.0  ) GOTO 901                       
            IF(I.EQ.1 .AND. NUMR.GT.255) GOTO 901                       
            IF(I.NE.1 .AND. NUMR.NE.0  ) GOTO 901                       
            ATAG = ATAG(2:J-1)                                          
            IREP = I                                                    
            GOTO 1                                                      
         ENDIF                                                          
         ENDDO                                                          
         GOTO 901                                                       
      ENDIF                                                             
      ENDDO                                                             
                                                                        
C  CHECK FOR VALID TAG                                                  
C  -------------------                                                  
                                                                        
1     IF(NEMOCK(ATAG).NE.0) GOTO 901                                    
      CALL NEMTAB(LUN,ATAG,IDN,TAB,IRET)                                
      IF(IRET.GT.0) THEN                                                

C	 Note that the next code line checks that we are not trying to
C	 replicate a Table B mnemonic (which is currently not allowed!).
C	 The logic works because, for replicated mnemonics, IREP = I =
C	 (the index within REPS(*,*) of the symbol associated with the
C	 type of replication in question (e.g. "{, "<", etc.))

         IF(TAB.EQ.'B' .AND. IREP.NE.0) GOTO 902                        
         IF(ATAG(1:1).EQ.'.') THEN                                      

C	    This mnemonic is a "following value" mnemonic
C	    (i.e. it relates to the mnemonic that immediately
C	    follows it within the user DX table sequence), so
C	    confirm that it contains, as a substring, this
C	    mnemonic that immediately follows it.

            NEMB = TAGS(N+1)                                            
            CALL NUMTAB(LUN,IDN,NEMA,TAB,ITAB)                          
            CALL NEMTAB(LUN,NEMB,JDN,TAB,IRET)                          
            CALL RSVFVM(NEMA,NEMB)                                      
            IF(NEMA.NE.ATAG) GOTO 903                                   
            IF(N.GT.NTAG ) GOTO 905                                     
            IF(TAB.NE.'B') GOTO 906                                     
         ENDIF                                                          
      ELSE                                                              
         GOTO 903                                                       
      ENDIF                                                             
                                                                        
C  WRITE THE DESCRIPTOR STRING INTO TABD ARRAY                          
C  -------------------------------------------                          
                                                                        
10    IF(IREP.GT.0) CALL PKTDD(ISEQ,LUN,IDNR(IREP,1)+NUMR,IRET)         
      IF(IRET.LT.0) GOTO 904                                            
      CALL PKTDD(ISEQ,LUN,IDN,IRET)                                     
      IF(IRET.LT.0) GOTO 904                                            
                                                                        
      ENDDO                                                             
                                                                        
      RETURN                                                            
                                                                        
C  ERROR EXITS                                                          
C  -----------                                                          
                                                                        
900   CALL BORT('SEQSDX - UNDEFINED SEQUENCE: '             //   NEMO) 
901   CALL BORT('SEQSDX - BAD TAG IN SEQUENCE: '            //TAGS(N)) 
902   CALL BORT('SEQSDX - REPLICATED ELEMENTS NOT ALLOWED:' //TAGS(N)) 
903   CALL BORT('SEQSDX - UNDEFINED TAG: '                  //TAGS(N)) 
904   CALL BORT('SEQSDX - TOO MANY DESCRIPTORS IN STRING:'  //   NEMO) 
905   CALL BORT('SEQSDX - FOLLOWING-VALUE LAST IN STRING:'  //   NEMA) 
906   CALL BORT('SEQSDX - FOLLOWING VALUE NOT FROM TABLEB:' //   NEMB) 
      END                                                               
