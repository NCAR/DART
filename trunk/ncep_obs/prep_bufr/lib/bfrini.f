      SUBROUTINE BFRINI                                                 

C************************************************************************
C* BFRINI								*
C*									*
C* This subroutine is called only one time (during the first call to	*
C* subroutine OPENBF) in order to initialize some global variables	*
C* and arrays within several COMMON blocks.				*
C*									*
C**									*
C* Log:									*
C* J. Woollen/NCEP	??/??						*
C* J. Ator/NCEP		05/01	Added documentation			*
C************************************************************************

      PARAMETER (NFILES=32)
                                                                        
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)     
      COMMON /PADESC/ IBCT,IPD1,IPD2,IPD3,IPD4                          
      COMMON /REPTAB/ IDNR(5,2),TYPS(5,2),REPS(5,2),LENS(5)             
      COMMON /STBFR / IOLUN(32),IOMSG(32)                               
      COMMON /TABABD/ NTBA(0:32),NTBB(0:32),NTBD(0:32),MTAB(50,32),     
     .                IDNA(50,32,2),IDNB(250,32),IDND(250,32),          
     .                TABA(50,32),TABB(250,32),TABD(250,32)             
      COMMON /DXTAB / MAXDX,IDXV,NXSTR(10),LDXA(10),LDXB(10),LDXD(10),  
     .                LD30(10),DXSTR(10)                                
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
      COMMON /BUFRMG/ MSGLEN,MSGTXT(5000)                               
      COMMON /MRGCOM/ NRPL,NMRG,NAMB,NTOT                               
      COMMON /DATELN/ LENDAT
      COMMON /QUIET / IPRT                                              
      COMMON /ACMODE/ IAC
                                                                        
      CHARACTER*600 TABD                                                
      CHARACTER*128 TABB                                                
      CHARACTER*128 TABA                                                
      CHARACTER*56  DXSTR                                               
      CHARACTER*10  TAG                                                 
      CHARACTER*6   ADSN(5,2),DNDX(25,10)                               
      CHARACTER*3   TYPX(5,2),TYPS,TYP                                  
      CHARACTER*1   REPX(5,2),REPS                                      
      DIMENSION     NDNDX(10),NLDXA(10),NLDXB(10),NLDXD(10),NLD30(10)   
      DIMENSION     LENX(5)                                             
                                                                        
      DATA ADSN   / '101000','360001','360002','360003','360004' ,      
     .              '101255','031002','031001','031001','031000' /      
      DATA TYPX   /    'REP',   'DRP',   'DRP',   'DRS' ,  'DRB' ,      
     .                 'SEQ',   'RPC',   'RPC',   'RPS' ,  'SEQ' /      
      DATA REPX   /      '"',     '(',     '{',     '[' ,    '<' ,      
     .                   '"',     ')',     '}',     ']' ,    '>' /      
      DATA LENX   /       0 ,     16 ,      8 ,      8  ,     1  /      
                                                                        
      DATA (DNDX(I,1),I=1,25)/                                          
     .'102000','031001','000001','000002',                              
     .'110000','031001','000010','000011','000012','000013','000015',   
     .                  '000016','000017','000018','000019','000020',   
     .'107000','031001','000010','000011','000012','000013','101000',   
     .                  '031001','000030'/                              
                                                                        
      DATA (DNDX(I,2),I=1,15)/                                          
     .'103000','031001','000001','000002','000003',                     
     .'101000','031001','300004',                                       
     .'105000','031001','300003','205064','101000','031001','000030'/   
                                                                        
      DATA NDNDX /  25 ,  15 , 8*0 /                                    
      DATA NLDXA /  35 ,  67 , 8*0 /                                    
      DATA NLDXB /  80 , 112 , 8*0 /                                    
      DATA NLDXD /  38 ,  70 , 8*0 /                                    
      DATA NLD30 /   5 ,   6 , 8*0 /                                    
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  INITIALIZE /BITBUF/                                                  
C  -------------------                                                  
                                                                        
      MAXBYT = 9970                                                     
                                                                        
C  INITIALIZE /PADESC/                                                  
C  -------------------                                                  
                                                                        
      IBCT = IFXY('063000')                                             
      IPD1 = IFXY('102000')                                             
      IPD2 = IFXY('031001')                                             
      IPD3 = IFXY('206001')                                             
      IPD4 = IFXY('063255')                                             
                                                                        
C  INITIALIZE /STBFR/                                                   
C  ------------------                                                   
                                                                        
      DO I=1,NFILES                                                     
      IOLUN(I) = 0                                                      
      IOMSG(I) = 0                                                      
      ENDDO                                                             
                                                                        
C  INITIALIZE /REPTAB/                                                  
C  -------------------                                                  
                                                                        
      DO I=1,5                                                          
      LENS(I) = LENX(I)                                                 
      DO J=1,2                                                          
      IDNR(I,J) = IFXY(ADSN(I,J))                                       
      TYPS(I,J) = TYPX(I,J)                                             
      REPS(I,J) = REPX(I,J)                                             
      ENDDO                                                             
      ENDDO                                                             
                                                                        
C  INITIALIZE /TABABD/                                                  
C  -------------------                                                  
                                                                        
C     NTBA(0) is the maximum number of entries within internal BUFR table A.

      NTBA(0) = 50                                                      
                                                                        
                                                                        
C     NTBB(0) is the maximum number of entries within internal BUFR table B.

      NTBB(0) = 250                                                     
                                                                        
                                                                        
C     NTBD(0) is the maximum number of entries within internal BUFR table D.

      NTBD(0) = 250                                                     
                                                                        
C  INITIALIZE /DXTAB/                                                   
C  ------------------                                                   
                                                                        
      MAXDX = MAXBYT                                                    
      IDXV  = 1                                                         
                                                                        
      DO J=1,10                                                         
      LDXA(J)  = NLDXA(J)                                               
      LDXB(J)  = NLDXB(J)                                               
      LDXD(J)  = NLDXD(J)                                               
      LD30(J)  = NLD30(J)                                               
      DXSTR(J) = '      '                                               
      NXSTR(J) = NDNDX(J)*2                                             
      DO I=1,NDNDX(J)                                                   
      I1 = I*2-1                                                        
      CALL IPKM(DXSTR(J)(I1:I1),2,IFXY(DNDX(I,J)))                      
      ENDDO                                                             
      ENDDO                                                             
                                                                        
C  INITIALIZE /TABLES/                                                  
C  -------------------                                                  
                                                                        
      MAXTAB = 15000                                                    
                                                                        
C  INITIALIZE /BUFRMG/                                                  
C  -------------------                                                  
                                                                        
      MSGLEN = 0                                                        
                                                                        
C  INITIALIZE /MRGCOM/                                                  
C  -------------------                                                  
                                                                        
      NRPL = 0                                                          
      NMRG = 0                                                          
      NAMB = 0                                                          
      NTOT = 0                                                          
                                                                        
C  INITIALIZE /DATELN/                                                  
C  -------------------                                                   
                                                                        
      IF(LENDAT.NE.10) LENDAT = 8                                                        
                                                                        
C  INITIALIZE /QUIET/                                                   
C  ------------------                                                   
                                                                        
      IPRT = 0                                                          

C  INITIALIZE /ACMODE/                                                  
C  ------------------                                                   
                                                                        
      IAC = 0                                                          

                                                                        
      RETURN                                                            
      END                                                               
