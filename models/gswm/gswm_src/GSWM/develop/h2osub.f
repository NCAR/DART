C**********************************************************************
      SUBROUTINE QDH2O(xh,Z,NFREQ,NZONAL,NMODE,Q)
      COMPLEX Q     
      DIMENSION QR(21),QI(21),XN(21),S1(4),S2(4)  
      DIMENSION QRJAN(5,7,21),QRAPR(5,7,21),QRJUL(5,7,21),QROCT(5,7,21)         
      DIMENSION QIJAN(5,7,21),QIAPR(5,7,21),QIJUL(5,7,21),QIOCT(5,7,21)         
      COMMON/QSEZUN/IQSEAS    
C 
C     Coefficients are given in units of mW/kg; Conversion to units of
C	J/kg/s is done in the SR (multiply by .001)	(M.Hagan --11/4/93)
C	(N.B. mW==>milliWatt and Watt=Newton*meter/sec=Joule/sec)
C 
C     This sub-routine provides various seasonal (mode by mode)
C	contributions to the heating rate due to water vapor  
C       in the lower troposphere using the formulation of Groves.  
C     The latitude structure is defined by the Hough Function and
C	accounted for in the main program. 
C     N.B. Groves Notation (m,s,n): m= 1 (diurnal), 2 (semidiurnal);
C		s= zonal wave # (=m for migrating modes); n= modal index
C     Coefficients are specified for 0-deg longitude and reckoned to 0GMT
C	must invoke a sign change for consistency in our model where fields
C	are reckoned to noon. No long. differences for migrating modes
C 
C NON-MIGRATING DIURNAL COMPONENTS: COMPLEX COEFFICIENTS
C         
C MODE (1,-1,1)
C         
      DATA (QRJAN(1,1,J),J=1,21)/-.30,-.17,-.11,-.05,-.05,-.03,-.08,  
     1 -.06,.04,.06,-.05,-.06,.03,.04,.02,6*0.0/
      DATA (QRAPR(1,1,J),J=1,21)/.05,.15,.04,.03,-.08,-.05,-.06,-.02, 
     1 .01,.01,-.09,-.18,-.22,-.16,-.09,-.00,.03,.01,.00,.00,-.00/    
      DATA (QRJUL(1,1,J),J=1,21)/.21,.24,.11,.04,-.11,-.09,-.08,.00,  
     1 -.05,-.09,-.22,-.29,-.41,-.47,-.34,-.19,-.09,-.04,-.02,-.01,-.00/        
      DATA (QROCT(1,1,J),J=1,21)/-.00,.02,-.10,-.09,-.11,-.04,-.03,.07,         
     1 .15,.17,-.10,-.28,-.21,-.12,-.07,-.04,-.03,-.01,-.01,-.00,-.00/
      DATA (QIJAN(1,1,J),J=1,21)/.10,.16,.18,.19,.19,.14,.07,-.03,-.09,         
     1 -.11,-.06,-.10,-.21,-.21,-.13,-.07,-.03,-.01,-.01,-.00,-.00/   
      DATA (QIAPR(1,1,J),J=1,21)/.28,.27,.10,.13,.17,.16,.16,.18,.15, 
     1 .13,.09,.05,-.08,-.15,-.12,-.15,-.12,-.06,-.02,-.01,-.00/      
      DATA (QIJUL(1,1,J),J=1,21)/.36,.35,.11,.11,.09,.05,.04,.05,.12, 
     1 .12,.13,.23,.36,.35,.23,.13,.06,.03,.01,.00,.00/     
      DATA (QIOCT(1,1,J),J=1,21)/.19,.17,.08,.08,.09,.05,.01,-.03,.00,
     1 .01,.05,.13,.20,.19,.13,.07,.03,.01,.01,.00,-.00/    
C         
C MODE (1,3,1) this is actually mode (1,3,3) using our notation      
C         
      DATA (QRJAN(5,1,J),J=1,21)/-.33,-.26,-.24,-.11,.03,.07,.09,.23, 
     1 .33,.30,.14,.21,.28,.20,.10,.04,.02,.01,.00,.00,.00/ 
      DATA (QRAPR(5,1,J),J=1,21)/-.08,-.09,-.20,-.11,-.03,.01,.04,.19,
     1 .27,.25,.18,.25,.29,.26,.19,.17,.12,.05,.02,.01,.00/ 
      DATA (QRJUL(5,1,J),J=1,21)/-.23,-.26,-.38,-.27,-.16,-.09,-.01,  
     1 .16,.29,.36,.24,.25,.42,.37,.23,.12,.05,.02,.01,.01,.00/       
      DATA (QROCT(5,1,J),J=1,21)/-.22,-.28,-.38,-.29,-.15,-.05,.03,.30,         
     1 .48,.51,.31,.28,.42,.37,.22,.11,.05,.02,.01,.00,.00/ 
      DATA (QIJAN(5,1,J),J=1,21)/-.28,-.25,-.20,-.24,-.26,-.18,-.12,  
     1 -.04,.05,.10,.09,.14,.27,.26,.17,.09,.04,.02,-.01,.00,.00/     
      DATA (QIAPR(5,1,J),J=1,21)/-.12,-.06,-.05,-.16,-.25,-.26,-.30,  
     1 -.35,-.25,-.15,-.08,.07,.35,.45,.38,.33,.22,.10,.04,.01,.01/   
      DATA (QIJUL(5,1,J),J=1,21)/-.16,-.15,-.15,-.20,-.22,-.15,-.08,.02,        
     1 .05,.07,-.02,-.03,-.04,-.06,-.05,-.03,-.02,-.01,.00,.00,.00/   
      DATA (QIOCT(5,1,J),J=1,21)/-.11,-.12,-.11,-.14,-.16,-.11,-.06,  
     1 -.02,-.00,-.00,-.01,.01,-.01,-.04,-.04,-.02,-.01,.00,.00,.00,.00/        
C         
C MODE (1,0,1)
C         
      DATA (QRJAN(2,1,J),J=1,21)/.14,.17,.17,.19,.20,.08,-.11,-.32,-.26,        
     1 -.15,-.22,-.40,-.30,-.13,-.05,-.02,-.01,-.01,.00,.00,.00/      
      DATA (QRAPR(2,1,J),J=1,21)/.29,.24,.22,.26,.29,.23,.20,.16,-.01,
     1 -.09,-.33,-.75,-.95,-.80,-.53,-.35,-.21,-.10,-.04,-.02,-.01/   
      DATA (QRJUL(2,1,J),J=1,21)/.28,.10,-.02,-.08,-.08,-.04,.04,.30, 
     1 .30,.27,-.05,-.30,-.33,-.27,-.17,-.10,-.05,-.03,-.01,-.01,.00/ 
      DATA (QROCT(2,1,J),J=1,21)/.17,.12,.05,-.06,-.12,-.12,-.11,.01, 
     1 .10,.13,-.14,-.27,-.17,-.09,-.05,-.03,-.02,-.01,-.01,.00,.00/  
      DATA (QIJAN(2,1,J),J=1,21)/-.01,.12,.02,-.07,-.21,-.14,-.13,-.11,         
     1 .04,.17,.09,.09,.31,.33,.21,.11,.05,.02,.01,.00,.00/ 
      DATA (QIAPR(2,1,J),J=1,21)/.22,.30,.30,.23,-.02,-.04,-.17,-.37, 
     1 -.23,-.02,-.10,-.17,.04,.07,.00,-.11,-.11,-.05,-.02,.00,.00/   
      DATA (QIJUL(2,1,J),J=1,21)/-.06,.06,.15,.11,-.05,.00,.01,.07,.12,         
     1 .08,-.18,-.21,-.19,-.26,-.21,-.12,-.06,-.03,-.01,-.01,.00/     
      DATA (QIOCT(2,1,J),J=1,21)/-.10,.00,.04,-.01,-.13,-.07,-.03,.04,
     1 .10,.22,.11,-.05,.08,.15,.11,.05,.02,.01,.00,.00,.00/
C         
C MODE (1,1,2) this is actually mode (1,2,2) using our notation      
C         
      DATA (QRJAN(4,1,J),J=1,21)/-.01,.10,.10,.09,.08,.02,-.07,-.16,  
     1 -.09,.00,-.03,-.07,.09,.18,.14,.08,.04,.02,.01,.00,.00/        
      DATA (QRAPR(4,1,J),J=1,21)/.12,.18,.18,.24,.27,.23,.24,.22,.05, 
     1 -.05,-.16,-.44,-.66,-.60,-.41,-.29,-.18,-.08,-.03,-.01,-.01/   
      DATA (QRJUL(4,1,J),J=1,21)/.17,.08,-.12,-.12,-.10,-.06,.01,.22, 
     1 .29,.29,.11,.06,.11,.09,.05,.02,.01,.00,.00,.00,.00/ 
      DATA (QROCT(4,1,J),J=1,21)/.04,.01,-.08,-.17,-.22,-.18,-.11,.07,
     1 .21,.25,.08,.13,.20,.23,.14,.07,.03,.01,.01,.00,.00/ 
      DATA (QIJAN(4,1,J),J=1,21)/.09,.08,.17,.17,.19,.13,.12,.06,-.13,
     1 -.29,-.19,-.21,-.52,-.54,-.35,-.18,-.08,-.04,-.02,-.01,.00/    
      DATA (QIAPR(4,1,J),J=1,21)/-.24,-.19,-.14,-.08,.06,.07,.16,.27, 
     1 .07,-.17,-.04,-.01,-.31,-.33,-.18,.02,.08,.04,.01,.00,.00/     
      DATA (QIJUL(4,1,J),J=1,21)/.11,.09,.09,.13,.20,.14,.12,.03,-.20,
     1 -.32,-.12,-.18,-.45,-.42,-.25,-.13,-.06,-.03,-.01,-.01,.00/    
      DATA (QIOCT(4,1,J),J=1,21)/-.02,-.01,.04,.07,.13,.10,.08,.01,-.11,        
     1 -.27,-.13,-.02,-.24,-.31,-.21,-.11,-.05,-.02,-.01,.00,.00/     
C         
C MIGRATING DIURNAL COMPONENTS: ONLY REAL COEFFICIENTS
C         
C MODE (1,1,1)      
C         
      DATA (QRJAN(3,1,J),J=1,21)/-2.70,-2.60,-2.68,-3.15,-3.59,-3.83, 
     1 -4.36,-5.39,-5.86,-5.90,-5.86,-6.55,-6.49,-4.86,-2.95,-1.59,   
     1 -.83,-.44,-.24,-.14,-.09/        
      DATA (QRAPR(3,1,J),J=1,21)/-2.70,-2.65,-2.60,-2.99,-3.35,-3.71, 
     1 -4.44,-5.88,-6.36,-6.29,-6.30,-7.14,-6.98,-5.19,-3.16,-1.83,   
     1 -1.01,-.52,-.27,-.15,-.10/       
      DATA (QRJUL(3,1,J),J=1,21)/-1.92,-2.13,-2.33,-2.81,-3.22,-3.49, 
     1 -4.00,-4.90,-5.30,-5.33,-5.44,-6.21,-6.05,-4.37,-2.60,-1.40,   
     1 -.74,-.40,-.22,-.14,-.09/        
      DATA (QROCT(3,1,J),J=1,21)/-2.06,-2.04,-2.17,-2.81,-3.44,-3.73, 
     1 -4.26,-5.34,-6.09,-6.24,-6.32,-7.41,-7.73,-5.90,-3.58,-1.92,   
     1 -.99,-.51,-.27,-.16,-.10/        
C         
C MODE (1,1,2)      
C         
      DATA (QRJAN(3,2,J),J=1,21)/-.86,-.91,-.82,-.60,-.25,-.13,-.08,  
     1 .16,.56,.74,.40,.48,.99,.98,.63,.32,.14,.06,.02,.01,.00/       
      DATA (QRAPR(3,2,J),J=1,21)/-.78,-.70,-.43,-.29,-.04,.02,.04,.12,
     1 .34,.48,.26,.30,.64,.69,.49,.36,.23,.11,.04,.02,.01/ 
      DATA (QRJUL(3,2,J),J=1,21)/.46,.37,.47,.47,.46,.35,.28,-.05,-.47,         
     1 -.50,-.20,-.44,-.85,-.73,-.43,-.21,-.09,-.04,-.01,-.01,.00/    
      DATA (QROCT(3,2,J),J=1,21)/-.01,-.11,-.04,-.04,.01,.00,.05,.17, 
     1 .08,-.08,.05,.18,-.05,-.18,-.15,-.08,-.04,-.02,-.01,.00,.00/   
C         
C MODE (1,1,3)      
C         
      DATA (QRJAN(3,3,J),J=1,21)/.88,.75,.84,.92,.90,.94,1.05,1.22,1.31,        
     1 1.37,1.43,1.54,1.55,1.21,.76,.41,.22,.11,.06,.04,.02/
      DATA (QRAPR(3,3,J),J=1,21)/1.73,1.52,1.4,1.24,.97,.90,1.,1.1,   
     1 1.,1.01,1.28,1.37,1.21,.91,.6,.3,.14,.07,.05,.03,.02/
      DATA (QRJUL(3,3,J),J=1,21)/.54,.65,1.01,.95,.82,.87,.98,1.09,   
     1 1.1,1.21,1.38,1.41,1.28,.95,.59,.32,.17,.09,.05,.03,.02/       
      DATA (QROCT(3,3,J),J=1,21)/.79,.71,.79,1.,1.07,1.08,1.27,1.54,  
     1 1.46,1.24,1.4,1.71,1.42,.89,.5,.27,.15,.08,.05,.03,.02/        
C         
C MODE (1,1,-1)      
C         
      DATA (QRJAN(3,4,J),J=1,21)/1.10,2.06,2.77,2.97,3.02,3.34,3.80,
     1 4.39,4.45,4.73,5.23,5.45,4.73,3.28,1.93,1.05,.57,.31,.18,.11,
     2 .08/
      DATA (QRAPR(3,4,J),J=1,21)/-1.13,-.89,-1.16,-1.02,-1.33,-1.41,
     1 -1.56,-1.83,-1.71,-1.70,-1.30,-.65,-.13,.06,.04,.01,-.01,-.02,
     2 -.03,-.03,-.03/
      DATA (QRJUL(3,4,J),J=1,21)/-1.55,-1.43,-1.94,-2.10,-2.68,-3.16,
     1 -3.84,-4.72,-4.90,-5.46,-5.73,-5.61,-5.19,-3.97,-2.46,-1.34,
     2 -.69,-.36,-.20,-.12,-.08/
      DATA (QROCT(3,4,J),J=1,21)/.51,1.47,1.88,1.44,.89,.94,1.02,1.22,
     1 1.14,1.18,1.06,.51,-.06,-.23,-.17,-.08,-.02,.01,.02,.02,.02/
C         
C         
C MODE (1,1,-2)      
C         
      DATA (QRJAN(3,5,J),J=1,21)/-5.35,-6.97,-8.52,-9.21,-10.10,-10.69,
     1 -11.45,-12.59,-12.69,-13.46,-14.37,-14.19,-11.90,-8.25,-5.00,
     2 -2.84,-1.60,-.94,-.58,-.39,-.29/
      DATA (QRAPR(3,5,J),J=1,21)/-5.69,-7.01,-8.78,-9.51,-10.58,-11.05,
     1 -11.73,-12.87,-12.84,-13.45,-13.88,-13.17,-10.72,-7.36,-4.51,
     2 -2.64,-1.53,-.90,-.56,-.38,-.28/
      DATA (QRJUL(3,5,J),J=1,21)/-4.93,-6.12,-7.84,-8.52,-9.55,-10.24,
     1 -11.18,-12.59,-12.70,-13.48,-13.85,-13.35,-11.55,-8.36,-5.15,
     2 -2.90,-1.60,-.91,-.55,-.36,-.26/
      DATA (QROCT(3,5,J),J=1,21)/-5.18,-7.22,-9.14,-9.58,-10.24,-10.96,
     1 -11.96,-13.49,-13.50,-14.08,-14.30,-13.62,-11.26,-7.76,-4.69,
     2 -2.68,-1.53,-.91,-.57,-.39,-.29/
C         
C MODE (1,1,-3)
C         
      DATA (QRJAN(3,6,J),J=1,21)/-.76,-.78,-.71,-.63,-.55,-.42,-.19,.09,
     1 .07,-.12,-.04,.43,.72,.59,.34,.15,.06,.01,.00,-.01,-.01/
      DATA (QRAPR(3,6,J),J=1,21)/.12,.18,.37,.43,.46,.40,.40,.47,.48,
     1 .35,.24,.26,.31,.26,.17,.10,.06,.04,.03,.02,.02/
      DATA (QRJUL(3,6,J),J=1,21)/.24,.50,.55,.64,.75,.69,.48,.20,.15,
     1 .13,-.01,-.25,-.59,-.63,-.42,-.21,-.08,-.02,.00,.01,.01/
      DATA (QROCT(3,6,J),J=1,21)/-.61,-.52,-.47,-.33,-.28,-.22,-.26,
     1 -.40,-.40,-.21,.04,-.01,-.19,-.21,-.15,-.09,-.05,-.03,-.02,
     2 -.01,-.01/
C         
C MODE (1,1,-4)
C         
      DATA (QRJAN(3,7,J),J=1,21)/-1.68,-2.04,-2.58,-2.90,-3.36,-3.73,
     1 -4.20,-4.78,-4.73,-4.87,-5.34,-5.70,-5.06,-3.58,-2.16,-1.20,
     2 -.65,-.36,-.21,-.13,-.09/
      DATA (QRAPR(3,7,J),J=1,21)/-2.01,-2.36,-2.77,-3.07,-3.48,-3.92,
     1 -4.38,-4.82,-4.72,-5.08,-5.71,-5.88,-4.98,-3.45,-2.08,-1.17,
     2 -.64,-.35,-.20,-.12,-.08/
      DATA (QRJUL(3,7,J),J=1,21)/-1.75,-2.08,-2.68,-2.88,-3.10,-3.38,
     1 -3.79,-4.39,-4.41,-4.66,-4.97,-5.03,-4.57,-3.41,-2.12,-1.18,
     2 -.63,-.35,-.20,-.13,-.09/
      DATA (QROCT(3,7,J),J=1,21)/-1.43,-1.91,-2.48,-2.98,-3.53,-4.00,
     1 -4.53,-5.10,-5.04,-5.39,-5.92,-6.11,-5.26,-3.64,-2.16,-1.18,
     2 -.63,-.34,-.20,-.12,-.08/
C         
C MIGRATING SEMIDIURNAL COMPONENTS: ONLY REAL COEFFICIENTS
C         
C MODE (2,2,2)      
C         
      DATA (QRJAN(2,2,J),J=1,21)/2.69,2.96,3.43,3.74,4.13,4.39, 
     1 4.81,5.47,5.60,5.81,6.09,6.19,5.29,3.54,1.99,1.03,   
     1 .54,.29,.17,.11,.08/        
      DATA (QRAPR(2,2,J),J=1,21)/2.89,3.15,3.64,3.94,4.33,4.63, 
     1 5.11,5.89,5.95,6.12,6.37,6.37,5.35,3.57,2.04,1.11,   
     1 .60,.32,.18,.11,.08/       
      DATA (QRJUL(2,2,J),J=1,21)/2.30,2.64,3.24,3.52,3.87,4.15, 
     1 4.57,5.23,5.34,5.58,5.74,5.72,4.96,3.40,1.92,.99,   
     1 .51,.28,.16,.10,.07/        
      DATA (QROCT(2,2,J),J=1,21)/2.45,2.94,3.54,3.90,4.33,4.70, 
     1 5.25,6.09,6.23,6.44,6.62,6.72,5.84,3.97,2.25,1.18,   
     1 .62,.34,.20,.13,.09/        
C         
C MODE (2,2,3)     ****READY******
C         
     DATA (QRJAN(2,3,J),J=1,21)/0.37,0.31,0.23,0.19,0.20,0.13,0.01,-.21,
     1 -.31,-.28,-.25,-.45,-.58,-.41,-.19,-.06,-.01,0.01,0.01,0.01, 
     1 0.01/
      DATA (QRAPR(2,3,J),J=1,21)/0.29,0.07,-.03,-.16,-.17,-.17,
     1 -.16,-.15,-.21,-.25,-.39,-.61,-.72,-.59,-.37,-.21,-.12,-.07,
     1 -.04,-.02,-.02/
      DATA (QRJUL(2,3,J),J=1,21)/0.01,-.30,-.36,-.42,-.38,-.31,-.15,
     1 0.10,0.25,0.33,0.34,0.42,0.58,0.47,0.23,0.07,0.00,-.02,-.02
     1 -.02,-.02/
      DATA (QROCT(2,3,J),J=1,21)/0.13,-.08,-.16,-.03,0.15,0.14,0.14,
     1 0.14,0.20,0.20,0.14,0.31,0.50,0.44,0.28,0.15,0.07,0.04,0.02,
     1 0.01,0.01/
C         
C MODE (2,2,4)     ****READY******
C         
      DATA (QRJAN(2,4,J),J=1,21)/-.29,-.56,-.84,-.85,-.96,-1.01,
     1 -.99,-.90,-.77,-.90,-1.05,-.76,-.31,-.10,-.06,-.06,-.05,
     1 -.04,-.03,-.03,-.02/
      DATA (QRAPR(2,4,J),J=1,21)/-.49,-.82,-1.23,-1.28,-1.41,-1.43,
     1 -1.38,-1.18,-.99,-1.16,-1.32,-.91,-.36,-.12,-.07,-.04,-.03,
     1 -.03,-.03,-.03,-.02/
      DATA (QRJUL(2,4,J),J=1,21)/-.45,-.64,-.93,-.91,-.93,-.95,-.93,
     1 -.89,-.78,-.88,-.88,-.51,-.18,-.09,-.06,-.03,-.02,-.02,-.02,
     1 -.02,-.01/
      DATA (QROCT(2,4,J),J=1,21)/-.57,-.99,-1.39,-1.34,-1.36,-1.46,
     1 -1.56,-1.58,-1.31,-1.40,-1.48,-.99,-.31,-.03,-.02,-.04,-.05,
     1 -.05,-.04,-.03,-.03/
C         
C MODE (2,2,6)     ****READY******
C         
      DATA (QRJAN(2,6,J),J=1,21)/-.11,0.10,0.15,0.23,0.40,0.41,0.39,
     1 0.42,0.47,0.52,0.49,0.42,0.32,0.22,0.14,0.09,0.06,0.04,0.03,
     1 0.02,0.01/
      DATA (QRAPR(2,6,J),J=1,21)/-.19,0.14,0.34,0.49,0.72,0.74,0.68,
     1 0.70,0.81,0.92,0.83,0.70,0.55,0.38,0.23,0.16,0.11,0.07,0.04,
     1 0.03,0.02/
      DATA (QRJUL(2,6,J),J=1,21)/-.03,0.11,0.12,0.23,0.36,0.34,0.33,
     1 0.42,0.46,0.44,0.36,0.36,0.34,0.22,0.12,0.07,0.04,0.02,0.02,
     1 0.01,0.01/
      DATA (QROCT(2,6,J),J=1,21)/0.02,0.35,0.57,0.57,0.65,0.69,0.65,
     1 0.66,0.78,0.95,0.87,0.67,0.58,0.44,0.28,0.17,0.10,0.06,0.04,
     1 0.03,0.02/
C
C     Set the water vapour heating rate Q to 0, then change to appropriate 
C     value if heating rate is specified for that mode and height less than
C     approximately 2 scale heights.  
C
  101 format(1x,"NMODE=",i2,1x,"NM=",i2)
      Q=(0.0,0.0)   
      IF(xh.GT.2.) GO TO 200  
C     Evaluate if current mode has a specified water vapour heating rate.  
C     If yes, then go to 50 and determine heating rate, else go to 200.  
	NM=NMODE
	IF (NMODE .LT. 0.) THEN 
	NM=3.+(-1*NMODE)
	ENDIF
c     print 101, nmode, nm
      IF(NFREQ .EQ. 1) THEN
       NS = NZONAL + 2 
       IF(NZONAL .EQ. 1) THEN 
        IF((NM .ge. 1.) .and. (NM .le. 7.)) THEN 
         GO TO 50 
        ENDIF 
       ENDIF 
       IF(NZONAL .EQ. 2) THEN 
        IF(NMODE .EQ. 2) THEN
	 NS=2
	 NM=NMODE
         GO TO 50 
        ENDIF
       ENDIF 
       IF(NZONAL .EQ. 3) THEN 
        IF(NMODE .EQ. 3) THEN 
         GO TO 50
        ENDIF 
       ENDIF 
      ELSE 
       IF(NFREQ .EQ. -1) THEN 
        IF(NZONAL .EQ. 1) THEN 
         GO TO 50 
        ENDIF 
       ENDIF 
      ENDIF 
C     if you get to here, then current mode does not have H2O heating rate.  
      GO TO 200 
C
 50   CONTINUE 
      DO 100 I=1,21 
       XN(I)=FLOAT(I-1)*.1
       GO TO(1,2,3,4),IQSEAS   
  1    QR(I)=QRJAN(NS,NM,I)
       QI(I)=QIJAN(NS,NM,I)   
       GO TO 100     
  2    QR(I)=QRAPR(NS,NM,I)
       QI(I)=QIAPR(NS,NM,I)   
       GO TO 100     
  3    QR(I)=QRJUL(NS,NM,I)
       QI(I)=QIJUL(NS,NM,I)   
       GO TO 100     
  4    QR(I)=QROCT(NS,NM,I)
       QI(I)=QIOCT(NS,NM,I)   
100   CONTINUE      
      CALL ATSM(Xh,XN,QR,21,1,S1,S2,4)   
      CALL ALI(Xh,S1,S2,QQR,4,1.0E-03,IER)         
      CALL ATSM(Xh,XN,QI,21,1,S1,S2,4)   
      CALL ALI(Xh,S1,S2,QQI,4,1.0E-03,IER)         
      Q=CMPLX(QQR,QQI)       
C Introduce sign change due to noon/midnight convention discussed
C	in comments above:
      Q=Q*(-1.**NZONAL)
      Q=.001*Q
C
200   RETURN        
      END 
C
C*************************************************************************
      SUBROUTINE CALHGHB(SIGMR,IS,IN,HF,DF,VERTWAVE)   
      PARAMETER (NPOINT=46,IA=180,IP=2*IA,NMOD=2)
      DIMENSION H(NPOINT),HZ(NPOINT),HZZ(NPOINT)
      DIMENSION AA(IA,IA),B(IA),PL(IP),UNITY(IA,IA) 
      DIMENSION WR(IA),WI(IA),IPNT(IA),XA(IA),WROLD(IA)    
      DIMENSION HF(2*NPOINT-1,2*NMOD+1),VERTWAVE(2*NMOD,3)  
      DIMENSION DF(2*NPOINT-1,2*NMOD+1) 
C       OPEN(UNIT=2,FILE='hout.dat',STATUS='new')
C      OPEN(UNIT=3,FILE='hf.dat',STATUS='formated') 
C
C      ****** CALCULATE THE HOUGH FUNCTIONS ******
C
C           ***** USER INPUT PARAMETERS  *****
C
C     SIGMR=FREQUENCY
C       POSITIVE VALUE MEANS WESTWARD PROPAGATION
C       NEGATIVE VALUE MEANS EASTWARD PROPAGATION
C
C     IS=ZONAL WAVENUMBER (ALWAYS POSITIVE)
C
C     IN=SYMMETERIC PARAMETER
C       MODE IS SYMMETRIC IF IN=IS
C       MODE IS ANTI-SYMMETRIC IF IN=IS+1
C
C     NPOINT=NUMBER OF POINTS TO CALCULATE IN LATITUDE DIRECTION
C       THE FIRST VALUE IS THE POLE
C       THE LAST VALUE IS THE EQUATOR 
C
C     NMOD=NUMBER OF MODES CALCULATED FOR EACH SET OF INPUT CONDITIONS
C
C
C       INITIAL CONDITIONS      
C       TEST= 0.003  
C       SIGMR=7.2722E-5
C       IS=1
C       IN=IS+1 
C      NMOD=4   
C
C  BASICALLY: YOU RUN THE PROGRAM FOR A FREQUENCY, A ZONAL WAVENUMBER
C  AND A GIVEN WANTED SYMETRY
C  THE ABOVE VALUES ARE FOR SYMETRIC PROPAGATING DIURNAL MODES....
C  HOUGH MODES ARE CALCULATED FROM POLE TO EQUATOR AND THERE ARE NPOINT's
C  FROM POLE TO EQUATOR (INCLUSIVE).
C  HR ARE THE EQUIVALENT DEPTHS IN KM
C  NMOD IS  THE NUMBER OF MODES CALCULATED...
C  ALSO THE GROVES CONVENTION FOR HOUGH MODE SIGNS IS USED. 
C  H ARE THE HOUGH MODES. FIRST & SECOND DERIVATIVES ARE ALSO CALCULATED
C  (THEY ARE HZ AND HZZ).
C
C **********************************************************************
      OMG=7.2921E-5
      S=FLOAT(IS)
      IDE=1
      IF(IN.EQ.IS+1) IDE=IDE+1
      IF(IN .EQ. IS+1) ISYMFLAG = 1 
      IF(IN .EQ. IS)   ISYMFLAG = 2 
      F=SIGMR/(2.*OMG)
C
C  INITIALIZE MATRICES TO ZEROS      
C      
      JIF=175   
      DO 3 I=1,IA
       WROLD(I) = 0. 
3     CONTINUE
C
C WHILE JIF LESS THAN 175, DO... 
11    CONTINUE
C
C INITIALIZE MATRIX TO ZEROS
C
      DO 1 I=1,IA
        DO 1 J=1,IA
        AA(I,J)=0.
1     CONTINUE
C
C  GENERATE THE TRIDIAGNAL MATRIX BASED ON THE INITIAL CONDITIONS
C  AND THE INTERATION STEP JIF
C
      XN=FLOAT(IN)
      JIF=JIF+5 
      CALL CHRG(AA,XN,S,F,JIF,IA)
      N=JIF
      K=1
      L=N
C
C CONSTRUCT IDENITY MATRIX, EIGENVECTORS WILL REPLACE THIS MATRIX  
C
      DO 23 JJ = 1,IA 
       DO 23 KK = 1,IA
        UNITY(KK,JJ) = 0.
        IF(KK .EQ. JJ) UNITY(KK,JJ) = 1. 
23    CONTINUE 
      WRITE(6,332)
332   FORMAT(1X,'Solving tri-diagonal matrix.....') 
C      WRITE(6,333) JIF,N
333   FORMAT(1X,'JIF= ',I5,2X,'N= ',I5) 
C
C CALL EIGENVALUE/EIGENVECTOR SOLVER
C
       CALL HQR2(IA,N,1,N,AA,WR,WI,UNITY,IERR) 
C        WRITE(6,334) IERR  
334   FORMAT('   Made it out of HQR2, IERR=',I5) 
C
C PUT WR INTO INCREASING ORDER (WR IS REAL, WI IS ZERO) 
C KEEPING TRACK OF ORIGINAL POSITION TO RELATE TO EIGENVECTOR.  
C EIGEN-ROUTINE CALCULATES NEGATIVE VALUES wrt GROVES
C CONVENTION.  THIS ORDERING IS FROM LOWEST GRAVITATION MODE TO 
C HIGHTEST GRAVITATION MODE, THEN ZEROS, THEN HIGHEST ROTATIONAL
C MODES THROUGH LOWEST ROTATIONAL MODES. (LOWER THE MODE, GREATER
C THE EQUIVALENT DEPTH) 
C 
      DO 40 KK = 1,IA
       IPNT(KK) = KK
40    CONTINUE 
      DO 41 JJ = 1,IA 
       DO 41 KK = 1,IA-1 
        IF(WR(KK) .GT. WR(KK+1)) THEN 
         TEMP = WR(KK+1)
         WR(KK+1) = WR(KK) 
         WR(KK) = TEMP 
         TEMP = IPNT(KK+1) 
         IPNT(KK+1) = IPNT(KK) 
         IPNT(KK) = TEMP 
        ENDIF 
41    CONTINUE 
C
C  TEST THE CONVERGENCE OF THE SOLUTIONS
C
      Do 30 KK = 1,20 
       XA(KK)=WR(KK)-WROLD(KK)
       WROLD(KK) = WR(KK) 
C       WRITE(2,31) JIF, XA(KK) 
30    CONTINUE
31    FORMAT(1X,'JIF= ',I4,' ERROR XA(KK) =',F20.17) 
C
C  TEST: IF ENOUGH ITERATION, THEN GO TO 11 
      IF(JIF .LT. 175) GOTO 11
C 
C  CALCULATE THE HOUGH FUNCTIONS THAT SATISFY THE CONVERGENCE TEST
C
C  XB = CONSTANT = -(4*OMG^2*a^2)/g [EXPRESSED IN km]      
      XB=-1./.011349
C      WRITE(2,1001) 
      DO 22 I=1,JIF
       HR=XB*REAL(WR(I))
C       WRITE(2,1002) HR
1001   FORMAT(' EQUIVALENT DEPTHS: ')
1002   FORMAT(1X,'HR= ',F25.17)   
C       PRINT 1000,HR
22    CONTINUE
C      PRINT 1100
C      PRINT 1400
C
C
C PUT THE LATITUDE VALUES INTO FIRST COLUMN OF MATRIX HF 
      DO 90 KK = 1,2*NPOINT-1 
       HF(KK,1) = (FLOAT(KK-1)*180./((NPOINT-1.)*2.))
       DF(KK,1) = (FLOAT(KK-1)*180./((NPOINT-1.)*2.)) 
90    CONTINUE 
      NCOL = 2 
      MODENUM = 0 
C
C FIND THE HOUGH FUNCTIONS FOR THE FIRST NMOD MODES OF FIRST
C THE GRAVITATIONAL MODES THEN THE ROTATIONAL MODES.  IN BOTH
C CASES, THE LOWEST ORDER MODE IS CALCULATED FIRST. 
C
      HTEST=ABS(XB*REAL(WR(JIF)))
      DO 100 LT=1,2
C
C TEST ON THE EXISTENCE OF SIGNIFICANT NEGATIVES MODES...
C      IF(LT.EQ.2.AND.HTEST.LT.0.01) STOP
      IF(LT.EQ.2.AND.HTEST.LT.0.01) RETURN  
C
C IF POSITIVE MODE, LOOK AT FIRST THRU NMOD FUNCTIONS
      IF(LT.EQ.1) THEN
        IDEB=1
        IFIN=NMOD
        ISTEP=1
        MODECNT = 0 
      ELSE
C IF NEGATIVE MODE, LOOK AT LAST THRU (LAST-NMOD) FUNCTIONS
        IDEB=JIF
        IFIN=JIF-NMOD+1
        ISTEP=-1
        MODECNT = 0 
      ENDIF
C
C  INITIALIZE H, HZ, HZZ and EIGENVECTOR B 
      DO 100 J=IDEB,IFIN,ISTEP
       MODENUM = MODENUM + 1 
       MODECNT = MODECNT + ISTEP 
      HR=XB*REAL(WR(J))
      DO 2 I=1,NPOINT
        H(I)=0.
        HZ(I)=0.
        HZZ(I)=0.
2     CONTINUE
      DO 63 I=1,JIF
        B(I)=REAL(UNITY(I,IPNT(J)))
63    CONTINUE
C
C  CALL SUBROUTINE TO CALCULATE THE NORMALIZED COEFFICIENTS
      CALL NRMCOF(B,JIF,IS,IDE,IA,IP)
C
      DO 3333 JJ=1,NPOINT
C  CALL SUBROUTINE TO CALCULATE THE ASSOCIATED LEGENDRE FUNCTIONS
C  THEN CALL SUBROUTINE TO NORMALIZE LEGENDRE FUNCTIONS
C  THEN CALL SUBROUTINE TO CALCULATE H(LAT), dH/d(LAT), d^2H/d(LAT)^2
        CALL CPL(IS,PL,NPOINT,IP,JJ)
        CALL NORM(IS,IDE,PL,NPOINT,IP,JIF)
        CALL FHOUGH(PL,B,H,JIF,NPOINT,IA,JJ,IP,IDE,HZ,HZZ,IS)
3333  CONTINUE
C
C USE GROVES CONVENTION FOR MODE SIGN
C
      IF(H(NPOINT-1).LT.0.) THEN
        DO 80 I=1,NPOINT
          H(I)=-H(I)
          HZ(I)=-HZ(I)
          HZZ(I)=-HZZ(I)
80      CONTINUE
        DO 81 I=1,JIF
          B(I)=-B(I)
81      CONTINUE
      ENDIF
C
C OUTPUT HOUGH FUNCTIONS AND EACH'S DERIVATIVES... 
C
      VERTWAVE(MODENUM,1) = MODECNT  
      VERTWAVE(MODENUM,2) = HR 
      TEMP = SQRT(ABS((1.-4.*(2./7.)*7.6/HR)/4.)) 
      VERTWAVE(MODENUM,3) = 2.*ACOS(-1.)*7.6/TEMP       
C      IF(HR .LT. 0) THEN  
C       VERTWAVE(MODENUM,3) = (-1.)*VERTWAVE(MODENUM,3) 
C      ENDIF 
C      PRINT 1400
      PRINT 1500,HR
C      WRITE(2,1500) HR 
      LL=0
      PI=ACOS(-1.)
      DX=PI/(2.*(NPOINT-1))
      DO 64 LL=1,NPOINT
       HI1=(H(LL+1)-H(LL-1))/(2.*DX)
       HI2=(H(LL+1)+H(LL-1)-2.*H(LL))/(DX*DX)
C       PRINT 3000,LL,H(LL),HZ(LL),HI1,HZZ(LL),HI2
C       WRITE(2,3000) LL-1,H(LL),HZ(LL),HI1,HZZ(LL),HI2 
C       PRINT 3000,LL,H(LL)
C       WRITE(10) H(LL)
 64   CONTINUE
      DO 95 KK = 1,NPOINT 
       HF(KK,NCOL) = H(KK)
       DF(KK,NCOL) = HZ(KK) 
95    CONTINUE 
      IF(ISYMFLAG .EQ. 1) THEN
       NROW = NPOINT-1 
       DO 96 KK = NPOINT+1,2*NPOINT-1 
        HF(KK,NCOL) = (-1.)*H(NROW) 
        DF(KK,NCOL) = HZ(NROW) 
        NROW = NROW - 1 
 96    CONTINUE 
      ENDIF 
      IF(ISYMFLAG .EQ. 2) THEN 
       NROW = NPOINT-1 
       DO 97 KK= NPOINT+1,2*NPOINT-1 
        HF(KK,NCOL) = H(NROW)   
        DF(KK,NCOL) = (-1.)*HZ(NROW) 
        NROW = NROW - 1
 97    CONTINUE 
      ENDIF 
      NCOL = NCOL + 1 
100   CONTINUE
C
C SAVE HF TO DISK 
C
C      DO 200 KK = 1,2*NPOINT-1 
C 200   WRITE(3,210) (HF(KK,JJ), JJ = 1,NCOL-1) 
210   FORMAT(1X,9(1X,F8.4)) 
C
1000  FORMAT(10X,2(E13.5))
1100  FORMAT(//)
1200  FORMAT(2X,6(E13.5))
1300  FORMAT(1X,I3,//)
1400  FORMAT(1H1)
1500  FORMAT(1X,'EQUIVALENT DEPTH=',F11.5,' km')
2000  FORMAT(15X,2(E13.5))
3000  FORMAT(1X,I4,5(2X,E13.5,2X))
4000  FORMAT(E13.6)
C
C
C       STOP  
      END
      SUBROUTINE CHRG(AA,XN,S,F,JIF,IA)
      DIMENSION AA(IA,IA)
C
C   CALCULATION OF THE MATRIX TO BE INVERTED
C
C   if zonal wavenumber, S, == 0 (or less than 1.0) and 
C      the first meridional index, XN, == 0,  
C   then increment XN by 2. 
      IF((S .LT. 1.0) .AND. (XN .LT. 1.0)) THEN 
       XN = XN + 2 
      ENDIF 
      JII=1
      DO 10 I=JII,JIF
      ALFAN=(XN-S+1.)*(XN-S+2.)/((2.*XN+1.)*(2.*XN+3.))
      ALFAN=ALFAN/((XN+1.)*(XN+2.)-S/F)
      GAMAN=(XN+S+1.)*(XN+S+2.)/((2.*XN+3.)*(2.*XN+5.))
      GAMAN=GAMAN/((XN+1.)*(XN+2.)-S/F)
      BETAN=-F*F*(XN*(XN+1.)-S/F)/((XN*(XN+1.))**2)
C   Don't execute the next two lines if (XN-1.0) == 0 
      IF((XN-1.) .GT. 0.0) THEN 
       XA=(XN-1.)*(XN-1.)*(XN-S)*(XN+S)/(XN*XN*(2.*XN+1.)*(2.*XN-1.))
       BETAN=BETAN+XA/(XN*(XN-1.)-S/F)
      ENDIF 
      XA=((XN+2.)**2)*(XN-S+1)*(XN+S+1)/(((XN+1.)**2)*(2*XN+1)*(2*XN+3))
      BETAN=BETAN+XA/((XN+1.)*(XN+2.)-S/F)
      AA(I,I)=BETAN
      AA(I,I+1)=GAMAN
      AA(I+1,I)=ALFAN
      XN=XN+2.
10    CONTINUE
      RETURN
      END
C
C
      SUBROUTINE NRMCOF(B,JIF,IS,IDE,IA,NZ)
      DIMENSION B(IA)
C
C  CALCULATION OF THE NORMALISED COEFFICIENTS FOR HOUGH MODES COEF.
C
      F=0.
      JPL=0
      NA=2*JIF
      DO 1 J=IDE,NA,2
      JPL=JPL+1
      M=2*IS
      N=J-1
      XNORM=2.*FACT(M,N)/(2.*IS+2.*J-1)
      XNORM=SQRT(XNORM)
      B(JPL)=B(JPL)*XNORM
      IF(ABS(B(JPL)).LT.1E-10) B(JPL)=0.
      F=F+ABS(B(JPL))**2
1     CONTINUE
      F=SQRT(F)
      DO 2 J=1,JIF
      B(J)=B(J)/F
2     CONTINUE
      RETURN
      END
C
C
      SUBROUTINE CPL(IS,PL,NX,NZ,I)
C
C  CALCULATE THE ASSOCIATED LEGENDRE FUNCTIONS
C
      DIMENSION PL(NZ)
      PI=ACOS(-1.)
      XX=(I-1)*PI/(2.*NX-2.)
      X=COS(XX)
      S=FLOAT(IS)
C
      FS=1.
      F2S=1.
      IFLAG = 0 
      IF(IS .EQ. 0) THEN 
       IS = 1
       IFLAG = 1
      ENDIF 
      DO 10 J=1,IS
      FS=FS*FLOAT(J)
10    CONTINUE
      DO 11 J=1,2*IS
      F2S=F2S*FLOAT(J)
11    CONTINUE
      IF(IFLAG .EQ. 1) IS = 0 
C
      PL(1)=(-SQRT(1.-X*X))**IS
      PL(1)=F2S*PL(1)/(FS*2.**IS)
      PL(2)=(2.*S+1.)*X*PL(1)
C
      JJ=NZ-1
      DO 1 J=2,JJ
      PL(J+1)  =((2.*S+2.*J-1.)*PL(J)  *X-(2.*S+J-1.)*PL(J-1)  )/J
1     CONTINUE
      RETURN
      END
C      
C
      SUBROUTINE NORM(IS,IDE,PL,NX,NZ,JIF)
      DIMENSION PL(NZ)
C
C  NORMALIZE THE LEGENDRE COEFFICIENTS
C
      JPL=0
      NA=JIF*2
C     DO 1 J=IDE,NA,2
      DO 1 J=1,NA
      JPL=JPL+1
      M=2*IS
      N=J-1
      XNORM=2.*FACT(M,N)/(2.*IS+2.*J-1)
      XNORM=SQRT(XNORM)
      PL(JPL)=PL(J)/XNORM
1     CONTINUE
      RETURN
      END
C
C
      FUNCTION FACT(M,N)
C
C  RATIO BETWEEN 2 FACTORIEL FUNCTIONS
C
      FACT=1.
      IF(M.NE.0) GOTO 1
      IF(N.EQ.0) RETURN
1     CONTINUE
      DO 2 I=1,M
      FACT=FACT*(M+N+1-I)
2     CONTINUE
      RETURN
      END
C     
C
      SUBROUTINE FHOUGH(PLP,B,H,JIF,NX,IA,I,NZ,IDE,HZ,HZZ,IS)
      DIMENSION PLP(NZ)
      DIMENSION H(NX),B(IA),HZ(NX),HZZ(NX)
C
C CALCULATION OF HOUGH MODE AND ITS DERIVATIVE
C
      PI=ACOS(-1.)
      L=IDE
      DO 2 J=1,JIF
      H(I)=H(I)+B(J)*PLP(L)
      L=L+2
2     CONTINUE
C
C
      HZ(I)=0.
      HZZ(I)=0.
      IF(I.EQ.1) RETURN
C
      SX=SIN((I-1)*PI/(2.*(NX-1)))
      CX=COS((I-1)*PI/(2.*(NX-1)))
      XM=FLOAT(IS+IDE-1)
C
      IF(IDE.EQ.1) THEN
          PLX=XM*CX*PLP(1)/SX
          PLXX=-XM*PLP(1)+(XM-1)*CX*PLX/SX
      ELSE
          AXM=SQRT((2.*XM+1)*(XM*XM-IS*IS)/(2.*XM-1))
          PLX=(XM*CX*PLP(2)-AXM*PLP(1))/SX
          PLXX=-XM*PLP(2)+(XM-1)*CX*PLX/SX -AXM*(XM-1)*CX*PLP(1)/(SX*SX)
      ENDIF
C
      HZ(I)=B(1)*PLX
      HZZ(I)=B(1)*PLXX
      L=IDE
C
      DO 3 J=2,JIF
      L=L+2
      XM=FLOAT(IS+L-1)
      XM1=XM-1
      AXM=SQRT((2.*XM+1)*(XM*XM-IS*IS)/(2.*XM-1))
      AXM1=SQRT((2.*XM1+1)*(XM1*XM1-IS*IS)/(2.*XM1-1))
      PLX=(XM*CX*PLP(L)-AXM*PLP(L-1))/SX
      PLXX=-XM*PLP(L)+(XM-1)*CX*PLX/SX-AXM*((XM-1)*CX*PLP(L-1)
     1-AXM1*PLP(L-2))/(SX*SX)
      HZ(I)=HZ(I)+B(J)*PLX
      HZZ(I)=HZZ(I)+B(J)*PLXX
3     CONTINUE
C
      RETURN
      END

      SUBROUTINE HQR2(NM,N,LOW,IGH,H,WR,WI,Z,IERR)
C
      INTEGER I,J,K,L,M,N,EN,II,JJ,LL,MM,NA,NM,NN,
     X        IGH,ITN,ITS,LOW,MP2,ENM2,IERR
      REAL H(NM,N),WR(N),WI(N),Z(NM,N)
      REAL P,Q,R,S,T,W,X,Y,RA,SA,VI,VR,ZZ,NORM,TST1,TST2
      LOGICAL NOTLAS
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE HQR2,
C     NUM. MATH. 16, 181-204(1970) BY PETERS AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
C     OF A REAL UPPER HESSENBERG MATRIX BY THE QR METHOD.  THE
C     EIGENVECTORS OF A REAL GENERAL MATRIX CAN ALSO BE FOUND
C     IF  ELMHES  AND  ELTRAN  OR  ORTHES  AND  ORTRAN  HAVE
C     BEEN USED TO REDUCE THIS GENERAL MATRIX TO HESSENBERG FORM
C     AND TO ACCUMULATE THE SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C
C        H CONTAINS THE UPPER HESSENBERG MATRIX.
C
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED BY  ELTRAN
C          AFTER THE REDUCTION BY  ELMHES, OR BY  ORTRAN  AFTER THE
C          REDUCTION BY  ORTHES, IF PERFORMED.  IF THE EIGENVECTORS
C          OF THE HESSENBERG MATRIX ARE DESIRED, Z MUST CONTAIN THE
C          IDENTITY MATRIX.
C
C     ON OUTPUT
C
C        H HAS BEEN DESTROYED.
C
C        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVALUES.  THE EIGENVALUES
C          ARE UNORDERED EXCEPT THAT COMPLEX CONJUGATE PAIRS
C          OF VALUES APPEAR CONSECUTIVELY WITH THE EIGENVALUE
C          HAVING THE POSITIVE IMAGINARY PART FIRST.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES SHOULD BE CORRECT
C          FOR INDICES IERR+1,...,N.
C
C        Z CONTAINS THE REAL AND IMAGINARY PARTS OF THE EIGENVECTORS.
C          IF THE I-TH EIGENVALUE IS REAL, THE I-TH COLUMN OF Z
C          CONTAINS ITS EIGENVECTOR.  IF THE I-TH EIGENVALUE IS COMPLEX
C          WITH POSITIVE IMAGINARY PART, THE I-TH AND (I+1)-TH
C          COLUMNS OF Z CONTAIN THE REAL AND IMAGINARY PARTS OF ITS
C          EIGENVECTOR.  THE EIGENVECTORS ARE UNNORMALIZED.  IF AN
C          ERROR EXIT IS MADE, NONE OF THE EIGENVECTORS HAS BEEN FOUND.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE LIMIT OF 30*N ITERATIONS IS EXHAUSTED
C                     WHILE THE J-TH EIGENVALUE IS BEING SOUGHT.
C
C     CALLS CDIV FOR COMPLEX DIVISION.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      NORM = 0.0E0
      K = 1
C     .......... STORE ROOTS ISOLATED BY BALANC
C                AND COMPUTE MATRIX NORM ..........
      DO 50 I = 1, N
C
         DO 40 J = K, N
   40    NORM = NORM + ABS(H(I,J))
C
         K = I
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 50
         WR(I) = H(I,I)
         WI(I) = 0.0E0
   50 CONTINUE
C
      EN = IGH
      T = 0.0E0
      ITN = 30*N
C     .......... SEARCH FOR NEXT EIGENVALUES ..........
   60 IF (EN .LT. LOW) GO TO 340
      ITS = 0
      NA = EN - 1
      ENM2 = NA - 1
C     .......... LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
C                FOR L=EN STEP -1 UNTIL LOW DO -- ..........
   70 DO 80 LL = LOW, EN
         L = EN + LOW - LL
         IF (L .EQ. LOW) GO TO 100
         S = ABS(H(L-1,L-1)) + ABS(H(L,L))
         IF (S .EQ. 0.0E0) S = NORM
         TST1 = S
         TST2 = TST1 + ABS(H(L,L-1))
         IF (TST2 .EQ. TST1) GO TO 100
   80 CONTINUE
C     .......... FORM SHIFT ..........
  100 X = H(EN,EN)
      IF (L .EQ. EN) GO TO 270
      Y = H(NA,NA)
      W = H(EN,NA) * H(NA,EN)
      IF (L .EQ. NA) GO TO 280
      IF (ITN .EQ. 0) GO TO 1000   
      IF (ITS .NE. 10 .AND. ITS .NE. 20) GO TO 130
C     .......... FORM EXCEPTIONAL SHIFT ..........
      T = T + X
C
      DO 120 I = LOW, EN
  120 H(I,I) = H(I,I) - X
C
      S = ABS(H(EN,NA)) + ABS(H(NA,ENM2))
      X = 0.75E0 * S
      Y = X
      W = -0.4375E0 * S * S
  130 ITS = ITS + 1
      ITN = ITN - 1
C     .......... LOOK FOR TWO CONSECUTIVE SMALL
C                SUB-DIAGONAL ELEMENTS.
C                FOR M=EN-2 STEP -1 UNTIL L DO -- ..........
      DO 140 MM = L, ENM2
         M = ENM2 + L - MM
         ZZ = H(M,M)
         R = X - ZZ
         S = Y - ZZ
         P = (R * S - W) / H(M+1,M) + H(M,M+1)
         Q = H(M+1,M+1) - ZZ - R - S
         R = H(M+2,M+1)
         S = ABS(P) + ABS(Q) + ABS(R)
         P = P / S
         Q = Q / S
         R = R / S
         IF (M .EQ. L) GO TO 150
         TST1 = ABS(P)*(ABS(H(M-1,M-1)) + ABS(ZZ) + ABS(H(M+1,M+1)))
         TST2 = TST1 + ABS(H(M,M-1))*(ABS(Q) + ABS(R))
         IF (TST2 .EQ. TST1) GO TO 150
  140 CONTINUE
C
  150 MP2 = M + 2
C
      DO 160 I = MP2, EN
         H(I,I-2) = 0.0E0
         IF (I .EQ. MP2) GO TO 160
         H(I,I-3) = 0.0E0
  160 CONTINUE
C     .......... DOUBLE QR STEP INVOLVING ROWS L TO EN AND
C                COLUMNS M TO EN ..........
      DO 260 K = M, NA
         NOTLAS = K .NE. NA
         IF (K .EQ. M) GO TO 170
         P = H(K,K-1)
         Q = H(K+1,K-1)
         R = 0.0E0
         IF (NOTLAS) R = H(K+2,K-1)
         X = ABS(P) + ABS(Q) + ABS(R)
         IF (X .EQ. 0.0E0) GO TO 260
         P = P / X
         Q = Q / X
         R = R / X
  170    S = SIGN(SQRT(P*P+Q*Q+R*R),P)
         IF (K .EQ. M) GO TO 180
         H(K,K-1) = -S * X
         GO TO 190
  180    IF (L .NE. M) H(K,K-1) = -H(K,K-1)
  190    P = P + S
         X = P / S
         Y = Q / S
         ZZ = R / S
         Q = Q / P
         R = R / P
         IF (NOTLAS) GO TO 225
C     .......... ROW MODIFICATION ..........
         DO 200 J = K, N
            P = H(K,J) + Q * H(K+1,J)
            H(K,J) = H(K,J) - P * X
            H(K+1,J) = H(K+1,J) - P * Y
  200    CONTINUE
C
         J = MIN0(EN,K+3)
C     .......... COLUMN MODIFICATION ..........
         DO 210 I = 1, J
            P = X * H(I,K) + Y * H(I,K+1)
            H(I,K) = H(I,K) - P
            H(I,K+1) = H(I,K+1) - P * Q
  210    CONTINUE
C     .......... ACCUMULATE TRANSFORMATIONS ..........
         DO 220 I = LOW, IGH
            P = X * Z(I,K) + Y * Z(I,K+1)
            Z(I,K) = Z(I,K) - P
            Z(I,K+1) = Z(I,K+1) - P * Q
  220    CONTINUE
         GO TO 255
  225    CONTINUE
C     .......... ROW MODIFICATION ..........
         DO 230 J = K, N
            P = H(K,J) + Q * H(K+1,J) + R * H(K+2,J)
            H(K,J) = H(K,J) - P * X
            H(K+1,J) = H(K+1,J) - P * Y
            H(K+2,J) = H(K+2,J) - P * ZZ
  230    CONTINUE
C
         J = MIN0(EN,K+3)
C     .......... COLUMN MODIFICATION ..........
         DO 240 I = 1, J
            P = X * H(I,K) + Y * H(I,K+1) + ZZ * H(I,K+2)
            H(I,K) = H(I,K) - P
            H(I,K+1) = H(I,K+1) - P * Q
            H(I,K+2) = H(I,K+2) - P * R
  240    CONTINUE
C     .......... ACCUMULATE TRANSFORMATIONS ..........
         DO 250 I = LOW, IGH
            P = X * Z(I,K) + Y * Z(I,K+1) + ZZ * Z(I,K+2)
            Z(I,K) = Z(I,K) - P
            Z(I,K+1) = Z(I,K+1) - P * Q
            Z(I,K+2) = Z(I,K+2) - P * R
  250    CONTINUE
  255    CONTINUE
C
  260 CONTINUE
C
      GO TO 70
C     .......... ONE ROOT FOUND ..........
  270 H(EN,EN) = X + T
      WR(EN) = H(EN,EN)
      WI(EN) = 0.0E0
      EN = NA
      GO TO 60
C     .......... TWO ROOTS FOUND ..........
  280 P = (Y - X) / 2.0E0
      Q = P * P + W
      ZZ = SQRT(ABS(Q))
      H(EN,EN) = X + T
      X = H(EN,EN)
      H(NA,NA) = Y + T
      IF (Q .LT. 0.0E0) GO TO 320
C     .......... REAL PAIR ..........
      ZZ = P + SIGN(ZZ,P)
      WR(NA) = X + ZZ
      WR(EN) = WR(NA)
      IF (ZZ .NE. 0.0E0) WR(EN) = X - W / ZZ
      WI(NA) = 0.0E0
      WI(EN) = 0.0E0
      X = H(EN,NA)
      S = ABS(X) + ABS(ZZ)
      P = X / S
      Q = ZZ / S
      R = SQRT(P*P+Q*Q)
      P = P / R
      Q = Q / R
C     .......... ROW MODIFICATION ..........
      DO 290 J = NA, N
         ZZ = H(NA,J)
         H(NA,J) = Q * ZZ + P * H(EN,J)
         H(EN,J) = Q * H(EN,J) - P * ZZ
  290 CONTINUE
C     .......... COLUMN MODIFICATION ..........
      DO 300 I = 1, EN
         ZZ = H(I,NA)
         H(I,NA) = Q * ZZ + P * H(I,EN)
         H(I,EN) = Q * H(I,EN) - P * ZZ
  300 CONTINUE
C     .......... ACCUMULATE TRANSFORMATIONS ..........
      DO 310 I = LOW, IGH
         ZZ = Z(I,NA)
         Z(I,NA) = Q * ZZ + P * Z(I,EN)
         Z(I,EN) = Q * Z(I,EN) - P * ZZ
  310 CONTINUE
C
      GO TO 330
C     .......... COMPLEX PAIR ..........
  320 WR(NA) = X + P
      WR(EN) = X + P
      WI(NA) = ZZ
      WI(EN) = -ZZ
  330 EN = ENM2
      GO TO 60
C     .......... ALL ROOTS FOUND.  BACKSUBSTITUTE TO FIND
C                VECTORS OF UPPER TRIANGULAR FORM ..........
  340 IF (NORM .EQ. 0.0E0) GO TO 1001
C     .......... FOR EN=N STEP -1 UNTIL 1 DO -- ..........
      DO 800 NN = 1, N
         EN = N + 1 - NN
         P = WR(EN)
         Q = WI(EN)
         NA = EN - 1
         IF (Q) 710, 600, 800
C     .......... REAL VECTOR ..........
  600    M = EN
         H(EN,EN) = 1.0E0
         IF (NA .EQ. 0) GO TO 800
C     .......... FOR I=EN-1 STEP -1 UNTIL 1 DO -- ..........
         DO 700 II = 1, NA
            I = EN - II
            W = H(I,I) - P
            R = 0.0E0
C
            DO 610 J = M, EN
  610       R = R + H(I,J) * H(J,EN)
C
            IF (WI(I) .GE. 0.0E0) GO TO 630
            ZZ = W
            S = R
            GO TO 700
  630       M = I
            IF (WI(I) .NE. 0.0E0) GO TO 640
            T = W
            IF (T .NE. 0.0E0) GO TO 635
               TST1 = NORM
               T = TST1
  632          T = 0.01E0 * T
               TST2 = NORM + T
               IF (TST2 .GT. TST1) GO TO 632
  635       H(I,EN) = -R / T
            GO TO 680
C     .......... SOLVE REAL EQUATIONS ..........
  640       X = H(I,I+1)
            Y = H(I+1,I)
            Q = (WR(I) - P) * (WR(I) - P) + WI(I) * WI(I)
            T = (X * S - ZZ * R) / Q
            H(I,EN) = T
            IF (ABS(X) .LE. ABS(ZZ)) GO TO 650
            H(I+1,EN) = (-R - W * T) / X
            GO TO 680
  650       H(I+1,EN) = (-S - Y * T) / ZZ
C
C     .......... OVERFLOW CONTROL ..........
  680       T = ABS(H(I,EN))
            IF (T .EQ. 0.0E0) GO TO 700
            TST1 = T
            TST2 = TST1 + 1.0E0/TST1
            IF (TST2 .GT. TST1) GO TO 700
            DO 690 J = I, EN
               H(J,EN) = H(J,EN)/T
  690       CONTINUE
C
  700    CONTINUE
C     .......... END REAL VECTOR ..........
         GO TO 800
C     .......... COMPLEX VECTOR ..........
  710    M = NA
C     .......... LAST VECTOR COMPONENT CHOSEN IMAGINARY SO THAT
C                EIGENVECTOR MATRIX IS TRIANGULAR ..........
         IF (ABS(H(EN,NA)) .LE. ABS(H(NA,EN))) GO TO 720
         H(NA,NA) = Q / H(EN,NA)
         H(NA,EN) = -(H(EN,EN) - P) / H(EN,NA)
         GO TO 730
  720    CALL CDIV(0.0E0,-H(NA,EN),H(NA,NA)-P,Q,H(NA,NA),H(NA,EN))
  730    H(EN,NA) = 0.0E0
         H(EN,EN) = 1.0E0
         ENM2 = NA - 1
         IF (ENM2 .EQ. 0) GO TO 800
C     .......... FOR I=EN-2 STEP -1 UNTIL 1 DO -- ..........
         DO 795 II = 1, ENM2
            I = NA - II
            W = H(I,I) - P
            RA = 0.0E0
            SA = 0.0E0
C
            DO 760 J = M, EN
               RA = RA + H(I,J) * H(J,NA)
               SA = SA + H(I,J) * H(J,EN)
  760       CONTINUE
C
            IF (WI(I) .GE. 0.0E0) GO TO 770
            ZZ = W
            R = RA
            S = SA
            GO TO 795
  770       M = I
            IF (WI(I) .NE. 0.0E0) GO TO 780
            CALL CDIV(-RA,-SA,W,Q,H(I,NA),H(I,EN))
            GO TO 790
C     .......... SOLVE COMPLEX EQUATIONS ..........
  780       X = H(I,I+1)
            Y = H(I+1,I)
            VR = (WR(I) - P) * (WR(I) - P) + WI(I) * WI(I) - Q * Q
            VI = (WR(I) - P) * 2.0E0 * Q
            IF (VR .NE. 0.0E0 .OR. VI .NE. 0.0E0) GO TO 784
               TST1 = NORM * (ABS(W) + ABS(Q) + ABS(X)
     X                      + ABS(Y) + ABS(ZZ))
               VR = TST1
  783          VR = 0.01E0 * VR
               TST2 = TST1 + VR
               IF (TST2 .GT. TST1) GO TO 783
  784       CALL CDIV(X*R-ZZ*RA+Q*SA,X*S-ZZ*SA-Q*RA,VR,VI,
     X                H(I,NA),H(I,EN))
            IF (ABS(X) .LE. ABS(ZZ) + ABS(Q)) GO TO 785
            H(I+1,NA) = (-RA - W * H(I,NA) + Q * H(I,EN)) / X
            H(I+1,EN) = (-SA - W * H(I,EN) - Q * H(I,NA)) / X
            GO TO 790
  785       CALL CDIV(-R-Y*H(I,NA),-S-Y*H(I,EN),ZZ,Q,
     X                H(I+1,NA),H(I+1,EN))
C
C     .......... OVERFLOW CONTROL ..........
  790       T = AMAX1(ABS(H(I,NA)), ABS(H(I,EN)))
            IF (T .EQ. 0.0E0) GO TO 795
            TST1 = T
            TST2 = TST1 + 1.0E0/TST1
            IF (TST2 .GT. TST1) GO TO 795
            DO 792 J = I, EN
               H(J,NA) = H(J,NA)/T
               H(J,EN) = H(J,EN)/T
  792       CONTINUE
C
  795    CONTINUE
C     .......... END COMPLEX VECTOR ..........
  800 CONTINUE
C     .......... END BACK SUBSTITUTION.
C                VECTORS OF ISOLATED ROOTS ..........
      DO 840 I = 1, N
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 840
C
         DO 820 J = I, N
  820    Z(I,J) = H(I,J)
C
  840 CONTINUE
C     .......... MULTIPLY BY TRANSFORMATION MATRIX TO GIVE
C                VECTORS OF ORIGINAL FULL MATRIX.
C                FOR J=N STEP -1 UNTIL LOW DO -- ..........
      DO 880 JJ = LOW, N
         J = N + LOW - JJ
         M = MIN0(J,IGH)
C
         DO 880 I = LOW, IGH
            ZZ = 0.0E0
C
            DO 860 K = LOW, M
  860       ZZ = ZZ + Z(I,K) * H(K,J)
C
            Z(I,J) = ZZ
  880 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- ALL EIGENVALUES HAVE NOT
C                CONVERGED AFTER 30*N ITERATIONS ..........
 1000 IERR = EN
 1001 RETURN
      END
      SUBROUTINE CDIV(AR,AI,BR,BI,CR,CI)
      REAL AR,AI,BR,BI,CR,CI
C
C     COMPLEX DIVISION, (CR,CI) = (AR,AI)/(BR,BI)
C
      REAL S,ARS,AIS,BRS,BIS
      S = ABS(BR) + ABS(BI)
      ARS = AR/S
      AIS = AI/S
      BRS = BR/S
      BIS = BI/S
      S = BRS**2 + BIS**2
      CR = (ARS*BRS + AIS*BIS)/S
      CI = (AIS*BRS - ARS*BIS)/S
      RETURN
      END


