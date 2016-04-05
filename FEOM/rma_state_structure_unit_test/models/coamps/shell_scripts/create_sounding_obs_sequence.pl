#!/usr/bin/perl
#
# This code may (or may not) be part of the COAMPS distribution,
# So it is not protected by the DART copyright agreement.
#
# DART $Id$
#
# AUTHOR:   T. R. Whitcomb
#           Naval Research Laboratory
#
# Creates a file that can be used as input to create_obs_sequence that
# defines a sequence consisting of temperature/wind/moisture soundings
# at standard pressure levels at certain latitude/longitude points.
# 
# There are no command line arguments - the script reads in latitude
# and longitude data from a "latlon.dat" file.
#
# This script is based on the MITGCM coseq.pl script by Jim Hansen
######

# Open Lat/Lon input file
open(LLIN,'<latlon.dat');

# Open output file - may modify later to be user-specified
open(OUT,'>coseq.dat');

# Define standard sounding pressure levels
@plevs = (1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100,
      70, 50);

# Define variable types
# These may be changed depending on the forward operators available
# See DART documentation on supplying observations
@varkinds= ( "AIR_TEMPERATURE", 
         "U_WIND_COMPONENT", 
         "V_WIND_COMPONENT", 
         "LOG_SPECIFIC_HUMIDITY" );
print "There are $#varkinds variables\n";

# Define error covariances
@varcovs = (2.5, 5.0, 5.0, .005);

# Output a large number for memory allocation purposes.  If the
# subsequent running of ./create_obs_sequence doesn't end with the
# line: 
#  write_obs_seq  opening formatted file set_def.out
# then this number must be made larger.
print OUT "9000000\n";

# Output the number of copies of data.  Since only creating
# a definition here, output 0
print OUT "0\n";

# Output the number of quality control values per field.  There
# is not QC used here.
print OUT "0\n";

# Loop over latitude/longitude locations
foreach $latlon (<LLIN>)
{
    chomp($latlon);
    @llinfo = split(/,/,$latlon);
    $lat = $llinfo[0];
    $lon = $llinfo[1];

    print "Processing sounding for ($lat, $lon)\n";

    # Loop over pressure levels
    foreach $pressure (@plevs)
    {
    print "\tPressure $pressure mb\n";
    # Loop over variables
    for ($ii=0; $ii < $#varkinds; $ii++)
    {
        $kind = $varkinds[$ii];
        $cov  = $varcovs[$ii];

            # For each variable, enter the following information:
        # 0 (enter another observation)
        # variable type
        # 2 (pressure coordinates)
            # pressure level
            # longitude
            # latitude
            # 0 0 (time - create_fixed_network_seq will handle this)
            # error variance
        print OUT "0\n";
        print OUT "$kind\n";
        print OUT "2\n";
        print OUT "$pressure\n";
        print OUT "$lon\n";
        print OUT "$lat\n";
        print OUT "0 0\n";
        print OUT "$cov\n";
    }
    }
}

# Output that we are done entering obs info
print OUT "-1\n";

# Output the input file name
print OUT "set_def.out\n";

# Clean up
print "Done creating sounding sequence.\n";

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

