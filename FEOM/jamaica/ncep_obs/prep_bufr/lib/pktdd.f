      SUBROUTINE PKTDD(ID,LUN,IDN,IRET)                                 

C************************************************************************
C* PKTDD								*
C*									*
C* Given a Table D sequence mnemonic (i.e. a "parent mnemonic") as	*
C* well as a mnemonic contained within that sequence (i.e. a "child	*
C* mnemonic", as determined within subroutine SEQSDX), this subroutine	*
C* stores information about the child mnemonic within the internal	*
C* BUFR Table D entry for the parent mnemonic.				*
C*									*
C* PKTDD  ( ID, LUN, IDN, IRET )					*
C*									*
C* Input parameters:							*
C*	ID		INTEGER		Positional index of parent	*
C*					mnemonic within internal BUFR	*
C*					Table D array TABD(*,*)		*
C*	LUN		INTEGER		I/O stream index of parent	*
C*					mnemonic within internal BUFR	*
C*					Table D array TABD(*,*)		*
C*	IDN		INTEGER		Bit-wise representation of FXY	*
C*					value corresponding to child	*
C*					mnemonic			*
C*					  0 = delete all information	*
C*					      about all child mnemonics	*
C*					      from within TABD(ID,LUN)	*
C* Output parameters:							*
C*	IRET		INTEGER		Total number of child mnemonics	*
C*					stored thus far (including IDN)	*
C*					for the parent mnemonic given	*
C*					by TABD(ID,LUN)			*
C*					  0 = information was cleared	*
C*					      from TABD(ID,LUN) because	*
C*					      input IDN value was 0	*
C*					 -1 = bad counter value *or*	*
C*					      maximum number of child	*
C*					      mnemonics already stored	*
C*					      for this parent mnemonic	*
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		05/01	Added documentation			*
C************************************************************************
                                                                        
      COMMON /TABABD/ NTBA(0:32),NTBB(0:32),NTBD(0:32),MTAB(50,32),     
     .                IDNA(50,32,2),IDNB(250,32),IDND(250,32),          
     .                TABA(50,32),TABB(250,32),TABD(250,32)             
      COMMON /DXTAB / MAXDX,IDXV,NXSTR(10),LDXA(10),LDXB(10),LDXD(10),  
     .                LD30(10),DXSTR(10)                                
                                                                        
      CHARACTER*600 TABD                                                
      CHARACTER*128 TABB                                                
      CHARACTER*128 TABA                                                
      CHARACTER*56  DXSTR                                               
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      LDD = LDXD(IDXV+1)+1                                              

C     LDD points to the byte within TABD(ID,LUN) which contains (in
C     packed integer format) a count of the number of child mnemonics
C     stored thus far for this parent mnemonic.
                                                                        
C  ZERO THE COUNTER IF IDN IS ZERO                                      
C  -------------------------------                                      
                                                                        
      IF(IDN.EQ.0) THEN                                                 
         CALL IPKM(TABD(ID,LUN)(LDD:LDD),1,0)                           
         IRET = 0                                                       
         RETURN                                                         
      ENDIF                                                             
                                                                        
C  UPDATE THE STORED DESCRIPTOR COUNT FOR THIS TABLE D ENTRY            
C  ---------------------------------------------------------            
                                                                        
      ND = IUPM(TABD(ID,LUN)(LDD:LDD),8)                                

C     ND is the (unpacked) count of the number of child mnemonics
C     stored thus far for this parent mnemonic.
                                                                        
      IF(ND.LT.0 .OR. ND.EQ.250) THEN                                   
         IRET = -1                                                      
         RETURN                                                         
      ELSE                                                              
         ND = ND+1                                                      
         CALL IPKM(TABD(ID,LUN)(LDD:LDD),1,ND)                          
         IRET = ND                                                      
      ENDIF                                                             
                                                                        
C  PACK AND STORE THE DESCRIPTOR                                        
C  -----------------------------                                        
                                                                        
      IDM = LDD+1 + (ND-1)*2                                            

C     IDM points to the starting byte within TABD(ID,LUN) at which
C     the IDN value for this child mnemonic will be stored (as a
C     packed integer of width = 2 bytes).

      CALL IPKM(TABD(ID,LUN)(IDM:IDM),2,IDN)                            
                                                                        
      RETURN                                                            
      END                                                               
