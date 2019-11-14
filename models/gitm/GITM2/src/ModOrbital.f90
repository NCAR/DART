module ModOrbital
  !Updated orbital elements for all GITM bodies.
  !Orbital Elements taken from Keplerian Elements for the Approximate Positions
  !of the Major Planets, E Standish, Solar System Dynamics, JPL, https://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf

  implicit none


!         Orbital Elements
!            a              e               I                L            long.peri.      long.node.
!         AU, AU/Cy     rad, rad/Cy     deg, deg/Cy      deg, deg/Cy      deg, deg/Cy     deg, deg/Cy
!-----------------------------------------------------------------------------------------------------------
!Mercury   0.38709927      0.20563593      7.00497902      252.25032350     77.45779628     48.33076593
!Venus     0.72333566      0.00677672      3.39467605      181.97909950    131.60246718     76.67984255
!EM Bary   1.00000261      0.01671123     -0.00001531      100.46457166    102.93768193      0.0
!Mars      1.52371034      0.09339410      1.84969142       -4.55343205    -23.94362959     49.55953891
!Jupiter   5.20288700      0.04838624      1.30439695       34.39644051     14.72847983    100.47390909
!Saturn    9.53667594      0.05386179      2.48599187       49.95424423     92.59887831    113.66242448
!Uranus   19.18916464      0.04725744      0.77263783      313.23810451    170.95427630     74.01692503
!Neptune  30.06992276      0.00859048      1.77004347      -55.12002969     44.96476227    131.78422574
!Pluto    39.48211675      0.24882730     17.14001206      238.92903833    224.06891629    110.30393684

!         Rates (per century)
!            a              e               I                L            long.peri.      long.node.
!         AU, AU/Cy     rad, rad/Cy     deg, deg/Cy      deg, deg/Cy      deg, deg/Cy     deg, deg/Cy
!-----------------------------------------------------------------------------------------------------------
!Mercury   0.00000037      0.00001906     -0.00594749   149472.67411175      0.16047689     -0.12534081
!Venus     0.00000390     -0.00004107     -0.00078890    58517.81538729      0.00268329     -0.27769418
!EM Bary   0.00000562     -0.00004392     -0.01294668    35999.37244981      0.32327364      0.0
!Mars      0.00001847      0.00007882     -0.00813131    19140.30268499      0.44441088     -0.29257343
!Jupiter   -0.00011607     -0.00013253     -0.00183714     3034.74612775      0.21252668      0.20469106
!Saturn    -0.00125060     -0.00050991      0.00193609     1222.49362201     -0.41897216     -0.28867794
!Uranus    -0.00196176     -0.00004397     -0.00242939      428.48202785      0.40805281      0.04240589
!Neptune    0.00026291      0.00005105      0.00035372      218.45945325     -0.32241464     -0.00508664
!Pluto     -0.00031596      0.00005170      0.00004818      145.20780515     -0.04062942     -0.01183482



  !Venus
  real, parameter :: semimajor_Venus = 0.72333566
  real, parameter :: eccentricity_Venus = 0.00677672
  real, parameter :: inclination_Venus = 3.39467605
  real, parameter :: longitudeNode_Venus = 76.67984255
  real, parameter :: longitudePerihelion_Venus = 131.60246718
  real, parameter :: meanLongitude_Venus = 181.97909950

  real, parameter :: semimajordot_Venus = 0.00000390
  real, parameter :: eccentricitydot_Venus = -0.00004107
  real, parameter :: inclinationdot_Venus = -0.00078890
  real, parameter :: longitudeNodedot_Venus = -0.27769418
  real, parameter :: longitudePeriheliondot_Venus =  0.00268329
  real, parameter :: meanLongitudedot_Venus = 58517.81538729


  !Earth

  real, parameter :: semimajor_Earth = 1.00000261
  real, parameter ::   eccentricity_Earth = 0.01671123
  real, parameter ::   inclination_Earth = -0.00001531
  real, parameter ::   longitudeNode_Earth =  0.0
  real, parameter ::   longitudePerihelion_Earth = 102.93768193
  real, parameter ::   meanLongitude_Earth = 100.46457166

  real, parameter ::   semimajordot_Earth =   0.00000562
  real, parameter ::   eccentricitydot_Earth = -0.00004392
  real, parameter ::   inclinationdot_Earth = -0.01294668
  real, parameter ::   longitudeNodedot_Earth = 0.0
  real, parameter ::   longitudePeriheliondot_Earth =  0.32327364
  real, parameter ::   meanLongitudedot_Earth = 35999.37244981

  !Mars
  real, parameter :: semimajor_Mars = 1.52371034
  real, parameter :: eccentricity_Mars = 0.09339410
  real, parameter :: inclination_Mars = 1.84969142
  real, parameter :: longitudeNode_Mars = 49.55953891
  real, parameter :: longitudePerihelion_Mars = -23.94362959
  real, parameter :: meanLongitude_Mars = -4.55343205

  real, parameter :: semimajordot_Mars = 0.00001847
  real, parameter :: eccentricitydot_Mars = 0.00007882
  real, parameter :: inclinationdot_Mars = -0.00813131
  real, parameter :: longitudeNodedot_Mars = -0.29257343
  real, parameter :: longitudePeriheliondot_Mars =  0.44441088
  real, parameter :: meanLongitudedot_Mars = 19140.30268499

  !Jupiter
  real, parameter :: semimajor_Jupiter = 5.20288700
  real, parameter :: eccentricity_Jupiter = 0.04838624
  real, parameter :: inclination_Jupiter = 1.30439695
  real, parameter :: longitudeNode_Jupiter = 100.47390909
  real, parameter :: longitudePerihelion_Jupiter = 14.72847983
  real, parameter :: meanLongitude_Jupiter = 34.39644051

  real, parameter :: semimajordot_Jupiter = -0.00011607
  real, parameter :: eccentricitydot_Jupiter =  -0.00013253
  real, parameter :: inclinationdot_Jupiter =  -0.00183714
  real, parameter :: longitudeNodedot_Jupiter = 0.20469106
  real, parameter :: longitudePeriheliondot_Jupiter =   0.21252668
  real, parameter :: meanLongitudedot_Jupiter = 3034.74612775

  !Saturn
  real, parameter :: semimajor_Saturn = 9.53667594
  real, parameter :: eccentricity_Saturn = 0.05386179
  real, parameter :: inclination_Saturn =  2.48599187
  real, parameter :: longitudeNode_Saturn = 113.66242448
  real, parameter :: longitudePerihelion_Saturn = 92.59887831
  real, parameter :: meanLongitude_Saturn = 49.95424423

  real, parameter :: semimajordot_Saturn = -0.00125060
  real, parameter :: eccentricitydot_Saturn = -0.00050991
  real, parameter :: inclinationdot_Saturn = 0.00193609
  real, parameter :: longitudeNodedot_Saturn = -0.28867794
  real, parameter :: longitudePeriheliondot_Saturn =  -0.41897216
  real, parameter :: meanLongitudedot_Saturn = 1222.49362201


end module ModOrbital
