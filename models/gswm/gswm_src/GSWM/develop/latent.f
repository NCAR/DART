      SUBROUTINE LATENT_HEAT_TIDES(Z,HJC,IJ)

c----- This subroutine calculates the latitude/height distribution of latent
c----- heating (UNITS = W/Kg or Joules/sec/Kg) corresponding to the diurnal
c----- and semidiurnal tides.  The vertical profile of heating is normalized
c----- to be consistent with a rainfall rate of 1.0 mm/day.  The latitudinal
c----- distribution of the amplitude and phase of the heating is calibrated
c----- according to the global rainfall analysis.

c----- The GSWM, due to historical reasons, defines the heating which maximizes
c----- at noon to be a real number.  That is, the forcing is given by

c			A cos w(t-tn)

c----- where w = frequency and t = hrs reckoned from local noon, and tn = time
c----- of maximum forcing reckoned from local noon.  In this notation, the
c----- complex heating must be specified as follows in the GSWM:

c		   A cos w(tn) - j A sin w(tn)

c----- (The conversion to phase defined as time of maximum reckoned from local
c----- midnight is made in the AMPHZ subroutine prior to printing the output.)

c----- The output array is defined to be complex in accord with the above 
c----- definitions:

C  Modifications to allow for 7-yr annual average migrating forcing
C	as well as 7-yr monthly mean migrating forcing
C					M. Hagan (5/2/96)

      integer heatmodel, dblhme11, diurno3, peako3, skipo3, skipir
	integer heatnm
	integer pwforce,pwheat,pwdelw
      integer call_late,o3conc
      real HJRa(2,91),HJIa(2,91)
      real hjrm,hjim
      COMPLEX XI
      COMPLEX HJC(91)

c  The following integer common block contains flags for the heat forcing
c       models

      common /heatforce/heatmodel,dblhme11,diurno3,peako3,skipo3,
     +     skipir,call_late,o3conc,heatnm,pwforce,pwheat,pwdelw

      COMMON/LATVECT/DTHETA,CLATT(91),XLAT(91),SNCLAT(91),CSCLAT(91),
     1    TNCLAT(91),CTNCLAT(91),SIN2I(91),GMLAT(91),DIP(91)

      COMMON/MODE/NZONAL,PERIOD,FREQ,MOIS,NSS,FLUX

c  Arrays containing 7-yr monthly mean latitude expansions are
C  	passed in common from main_2.f which calls SR getlate
      common/lateforce/HJRm(91),HJIm(91)

C Migrating Diurnal Forcing: 7-year annual average 
C	~/spectral/results/interplatent.pro output with
C	~jforbes/latent_heat/T24-S1-annavg.dat input

      data (hjra(1,i),i=1,91)/0.000572743,0.000773122,0.00104361,
     +	0.00140872,0.00190157,0.00256686,0.00346489,0.00467711,
     +	0.00631344,0.00852226,0.0115038,0.0155286,0.0209614,
     +	0.0282949,0.0381941,0.0515567,0.0695942,0.0939424,0.127287,
     +	0.171995,0.214653,0.181478,0.187683,0.181369,0.156932,
     +	0.330722,0.450905,0.521700,0.361210,0.276349,0.179461,
     +	0.311632,0.480853,0.385978,0.213010,0.108221,0.0630620,
     +	0.0952637,0.131077,0.135777,0.100586,0.0745162,0.0552029,
     +	0.0408953,0.0302960,0.0224438,0.0166268,0.0123174,.00912498,
     +	0.00675995,0.00500790,0.00370994,0.00274839,0.00203606,
     +	0.00150835,0.00111741,0.000827800,0.000613249,0.000454306,
     +	32*0.0/

      data (hjia(1,i),i=1,91)/-0.000737576,-0.000995623,-0.00134395,
     +	-0.00181414,-0.00244884,-0.00330559,-0.00446207,-0.00602317,
     +	-0.00813043,-0.0109749,-0.0148146,-0.0199976,-0.0269940,
     +	-0.0364380,-0.0491862,-0.0663945,-0.0896231,-0.120979,
     +	-0.163920,-0.221495,-0.182601,-0.198194,-0.313619,-0.369593,
     +	-0.431212,-0.533658,-0.585141,-0.461598,-0.489609,-0.731609,
     +	-0.810040,-0.715676,-0.874705,-0.732490,-0.706121,-0.592043,
     +	-0.475536,-0.323367,-0.228546,-0.146006,-0.108164,
     +	-0.0801296,-0.0593615,-0.0439761,-0.0325783,-0.0241346,
     +	-0.0178793,-0.0132453,-0.00981239,-0.00726920,-0.00538515,
     +	-0.00398942,-0.00295543,-0.00218944,-0.00162198,-0.00120159,
     +	-0.000890160,-0.000659446,-0.000488530,
     +	32*0.0/

C Migrating Semidiurnal Forcing: 7-year annual average 
C	~/spectral/results/interplatent.pro output with
C	~jforbes/latent_heat/T12-S2-annavg.dat input

      data (hjra(2,i),i=1,91)/-2.63896e-05,-3.56222e-05,-4.80850e-05,
     +	-6.49079e-05,-8.76165e-05,-0.000118270,-0.000159648,-0.000215502,
     +	-0.000290897,-0.000392670,-0.000530049,-0.000715492,-0.000965813,
     +	-0.00130371,-0.00175983,-0.00237552,-0.00320661,-0.00432847,
     +	-0.00586486,-0.00792482,-0.0123926,0.00401927,0.0422615,
     +	0.0175071,0.0130455,0.0244497,0.0315178,0.0488436,
     +	0.100110,0.129512,0.130922,0.0781160,0.0583962,
     +	0.0103068,0.0316945,0.0256340,0.0128404,-0.00354310,
     +	-0.0135051,-0.0303314,-0.0224700,-0.0166462,-0.0123318,
     +	-0.00913563,-0.00676784,-0.00501374,-0.00371427,-0.00275160,
     +	-0.00203843,-0.00151011,-0.00111872,-0.000828766,-0.000613965,
     +	-0.000454836,-0.000336951,-0.000249619,-0.000184923,-0.000136994,
     +	-0.000101488,
     +	32*0.0/

	data (hjia(2,i),i=1,91)/0.000474477,0.000640476,
     +	0.000864553,0.00116702,0.00157532,0.00212646,
     +	0.00287042,0.00387466,0.00523024,0.00706008,0.00953012,
     +	0.0128643,0.0173650,0.0234403,0.0316411,0.0427110,
     +	0.0576539,0.0778246,0.105448,0.142486,0.159429,
     +	0.166440,0.191228,0.232602,0.262046,0.375450,0.536147,
     +	0.577759,0.551971,0.574424,0.590298,0.581600,0.577796,
     +	0.440608,0.328338,0.223968,0.150446,0.124989,0.129116,
     +	0.119259,0.0883493,0.0654507,0.0484871,0.0359201,
     +	0.0266103,0.0197134,0.0146040,0.0108189,0.00801486,
     +	0.00593756,0.00439865,0.00325860,0.00241403,0.00178836,
     +	0.00132485,0.000981472,0.000727092,0.000538643,
     +	0.000399037,32*0.0/


      DO 2 I=1,IJ
      HJC(I)=(0.0,0.0)
2     CONTINUE
      IF(Z.GT.20.) GO TO 99
c     go to 99

      PI=ACOS(-1.)
      XI=(0.0,1.0)
      TPI=2.*PI

c------------------------------------------------------------------------------
c  The following expression for the vertical variation of latent heating due to
c  cumulus convection is from Hong and Wang (1980, Bull. Geophys., 19, 56-84),
c  which is based on the work of Reed and Recker (1971) and Nitta (1972). The
c  .00534 factor (W/Kg) translates to a rainfall rate of 1.0 mm/day.

      FH=.00534*(EXP(-(Z-6.5)*(Z-6.5)/29.05)-0.23*EXP(-Z/1.31))

c------------------------------------------------------------------------------
      DO 1 I=1,IJ
      X=CLATT(I)
      Y=XLAT(I)
      THETA=Y*180./PI
      SINY=SIN(Y)
      COSY=COS(Y)

C Turn off monthly mean (0) and turn on special case (1) with call_late flag
C For special case of 7-year annual average:

      if(call_late.eq.1)then  !special case annual tables above

         HJC(I)=FH*(HJRa(NZONAL,I)-XI*HJIa(NZONAL,I))

C For 7-year monthly mean     !normal case--call get_latent from main_2
      elseif(call_late.eq.0)then

         HJC(I)=FH*(HJRm(I)-XI*HJIm(I))

      endif

1     CONTINUE

99    RETURN
      END

