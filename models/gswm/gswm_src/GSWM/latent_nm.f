      SUBROUTINE LATENT_HEAT_KM(Z,HLJ,IJ)

c----- This subroutine calculates the latitude/height distribution of latent
c----- heating (UNITS = W/Kg or Joules/sec/Kg) for a given season and zonal
c----- wavenumber. Complex heating rates as a function of latitude are
c----- read in the SR SETGCI (called from MAIN_2)  and passed in common.
c----- These data are stored in ~hagan/GSWM/dat/gciforce (M. Hagan....7/97)

      complex HLJ(91),gciforce(91)

      COMMON/LATVECT/DTHETA,CLATT(91),XLAT(91),SNCLAT(91),CSCLAT(91),
     1    TNCLAT(91),CTNCLAT(91),SIN2I(91),GMLAT(91),DIP(91)
      COMMON /MODEi/NZONAL,MOIS,NSS
      COMMON /MODEr/PERIOD,FREQ,FLUX
      COMMON/gci/gciforce

      DO 2 I=1,IJ
      HLJ(I)=(0.0,0.0)
2     CONTINUE
      IF(Z.GT.20.) GO TO 99
c     go to 99

c------------------------------------------------------------------------------
c  The following expression for the vertical variation of latent heating due to
c  cumulus convection is from Hong and Wang (1980, Bull. Geophys., 19, 56-84),
c  which is based on the work of Reed and Recker (1971) and Nitta (1972). The
c  .00534 factor (W/Kg) translates to a rainfall rate of 1.0 mm/day.

      FH=.00534*(EXP(-(Z-6.5)*(Z-6.5)/29.05)-0.23*EXP(-Z/1.31))

      DO 1 I=1,IJ
      
         HLJ(I)=FH*gciforce(i)

1     CONTINUE

99    RETURN
      END

