C---------------------------------------------------------------------- 
C  ENTRY TO RESET STRING CACHE                                          
C---------------------------------------------------------------------- 
      SUBROUTINE STRCLN                                                 
      PARAMETER(MXS=1000)                                               
      COMMON /STCACH/ MSTR,NSTR,LSTR,LUNS(MXS,2),USRS(MXS),ICON(52,MXS)
      CHARACTER*80 USRS                                                 
                                                                        
      MSTR = MXS                                                        
      NSTR = 0                                                          
      LSTR = 0                                                          
      RETURN                                                            
      END                                                               
