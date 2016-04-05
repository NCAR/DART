      SUBROUTINE INCTAB(ATAG,ATYP,NODE)                                 

C************************************************************************
C* INCTAB								*
C*									*
C* This subroutine returns the next available positional index for	*
C* writing into the internal jump/link table, and it also uses that	*
C* index to store ATAG and ATYP within, respectively, the internal	*
C* jump/link table arrays TAG(*) and TYP(*).  If there is no more room	*
C* for additional entries within the internal jump/link table, then an	*
C* appropriate call is made to subroutine BORT.				*
C*									*
C* INCTAB  ( ATAG, ATYP, NODE )						*
C*									*
C* Input parameters:							*
C*	ATAG		CHARACTER*(*)	Mnemonic name			*
C*	ATYP		CHARACTER*(*)	Mnemonic type			*
C*									*
C* Output parameters:							*
C*	NODE		INTEGER		Next available positional index	*
C*					for writing into the internal	*
C*					jump/link table			*
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		05/01	Added documentation			*
C************************************************************************
                                                                        
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
                                                                        
      CHARACTER*(*) ATAG,ATYP                                           
      CHARACTER*10  TAG                                                 
      CHARACTER*3   TYP                                                 
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      NTAB = NTAB+1                                                     
      IF(NTAB.GT.MAXTAB) CALL BORT('INCTAB - TOO MANY ENTRIES')        
                                                                        
      TAG(NTAB) = ATAG                                                  
      TYP(NTAB) = ATYP                                                  
      NODE = NTAB                                                       
                                                                        
      RETURN                                                            
      END                                                               
