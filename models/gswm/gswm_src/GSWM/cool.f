	SUBROUTINE COOL(Z,AZ)

C  Newtonian cooling due to O3, from Salby (1981, JAS, 38, p.1809;
C  see also his reference to Hartmann, 1978 (peaks at .35 /day):


	Arg1=((Z-50.)/21.3)*((Z-50.)/21.3)
	AZ=4.051E-06*EXP(-Arg1)

C  Newtonian cooling due to CO2 (cf. Zhu and Strobel, JAS, 48, 184)
C  for vertical wavelengths between about 40 km and infinity (peaks
C  at .20/day):

	Arg2=((Z-80.)/20.)*((Z-80.)/20.)
	AZ=AZ+2.315E-06*EXP(-Arg2)

	RETURN
	END

