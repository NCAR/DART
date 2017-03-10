BEGIN  { print "0" ;
         print "0" ;
       }
         
# extra metadata
/^COSMOS_NEUTRON_INTENSITY/  {  
                               ecount = 8;
                               evalues[1] = "10" ;
                               evalues[2] = "10" ;
                               evalues[3] = "10" ;
                               evalues[4] = "10" ;
                               evalues[5] = "10" ;
                               evalues[6] = "10" ;
                               evalues[7] = "10" ;
                               evalues[8] = "10" ;
                            }

/^DOPPLER_RADIAL_VELOCITY/    {
                               ecount = 7 ;
                               evalues[1] = "2" ;
                               evalues[2] = "400" ;
                               evalues[3] = "245" ;
                               evalues[4] = "42" ;
                               evalues[5] = "45" ;
                               evalues[6] = "45" ;
                               evalues[7] = "15" ;
                            }

/^GEO_CO_ASI|^GEO_CO_EUR|^GEO_CO_NAM|^IASI_CO_RETRIEVAL|^MOPITT_CO_RETRIEVAL/   {  
                   ecount = 12 ;
                   evalues[1]  = 100.0 ;
                   evalues[2]  = 1000.0 ;
                   evalues[3]  = 0.1 ;
                   evalues[4]  = 0.1 ;
                   evalues[5]  = 0.1 ;
                   evalues[6]  = 0.1 ;
                   evalues[7]  = 0.1 ;
                   evalues[8]  = 0.1 ;
                   evalues[9]  = 0.1 ;
                   evalues[10] = 0.1 ;
                   evalues[11] = 0.1 ;
                   evalues[12] = 0.1 ;
                }

/^GPSRO_REFRACTIVITY/  {
                   ecount = 1 ;
                   evalues[1]  = 1 ;
                   specialvert = 3;
                   specialval = 2000 ;
               }

/^BOB/  {
        }

# special vert restrictions?


# the common values for any obs

       { print "0" ;
         print $1  ;
         if ( ecount > 0 ) {
           for (i=1; i<=ecount; i++)
              print evalues[i];
         }
         if ( specialvert == 0) {
            print "2" ;
            print "500" ;
         } else {
            print specialvert ;
            print specialval ;
         }
         print "260" ;
         print "40" ;
         print "2015 10 30 15 30 0" ;
         print "1.0" ;
         ecount = 0 ;
         specialvert = 0 ;
       }

END    { print "obs_seq.in" }

