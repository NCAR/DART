#!/usr/bin/perl -w
#
# This code may (or may not) be part of the COAMPS distribution,
# So it is not protected by the DART copyright agreement.
#
# DART $Id$
#
########################################################################
#
# SCRIPT:   coamps_innovation_to_obs.pl
# AUTHOR:   T. R. Whitcomb
#           Naval Research Laboratory
#
# Read in a COAMPS innovation file, obtain the relevant data from the 
# observations, and translate this into a format that is written to a 
# file to be used as an input to DART's create_obs_sequence program
# to generate an observing network.
#
# The command line arguments specify the name of the innovation file
# and the name of the output file.
#
# N.B. This is just a bare-bones effort to get something in place that
#      works, but this doesn't handle very many observations yet.
######
use POSIX qw(floor);
print "Usage: innovation_to_obs.pl innovation_file output_file\n" and
    die unless $#ARGV>=1;

# Copy over input variables
$innovation_file_name = $ARGV[0];
$output_file_name     = $ARGV[1];

# Open the files
die "Couldn't open $innovation_file_name!\n" unless 
    open(INNOV,"<$innovation_file_name");
die "Couldn't open $output_file_name!\n" unless 
    open(OUTPUT,">$output_file_name");

print "Innovation file $innovation_file_name opened successfully.\n";
print "Output file $output_file_name opened successfully.\n";

# Do some pre-processing - create the hash tables that we will use to
# translate between the NAVDAS codes and the DART codes for
# observation types
&build_navdas_hash;
&build_dart_hash;

# Read  in the header information from the innovation file - for now,
# this is just discarded as I don't see how it gets used within the
# DART framework
print "\nProcessing header...\n";
do {$line = <INNOV>} until $line =~ m/number of obs/;
chomp($line);
$line =~ m/\D*(\d+)\D*(\d+)\D*(\d+)\D*/;
$n_obs = $1;
$cdtg = $2;
$tau = $3;
print "Number of observations in sequence : $n_obs\n";
print "Date-time group of forecast cycle  : $cdtg\n";
print "Lead time of current analysis      : $tau\n";

# Write out the boilerplate for create_obs_sequence
print OUTPUT "$n_obs\n";    # Upper bound on obs
print OUTPUT "1\n";         # Set 1 copy of the data
print OUTPUT "0\n";         # no QC fields
print OUTPUT "INOV FILE\n"; # metadata for sequence

# Skip the title line
<INNOV>;
print "\n";

# Start observation loop
while(<INNOV>)
{
    # Note that this split assumes that everything is just separated
    # by spaces - for some observations, the "platform ID" has a space
    # in it and sometimes it does not.  Currently, since the only
    # things after that are the database address of the original
    # value, the surface pressure change and the background humidity,
    # this should not pose a problem at this point but it is something
    # to keep an eye out for later in which case I'll have to replace
    # this with a more complicated regex. - this is evidently picking
    # up something null at the beginning since the 0th entry of the
    # array is blank
    @obs_fields = split(/\s+/);
    print "Processing observation $obs_fields[1]...\n";

    # Convert the innovation file's variable type to the DART type
    $dart_type = &navdas_to_obs_kind($obs_fields[11]);
    if ($dart_type == -999)
    {
    print "\tSkipping due to incompatible types...\n";
    next;
    }

    # Calculate the "time" of the observation - considering that the
    # DTG matches what we have, then the "time" is given by the lead
    # time minus the "time from analysis time" field
    $tau_secs = $tau * 3600;
    $obs_time = $tau_secs + $obs_fields[15];

    # Break up the time into days and seconds (86400 secs/day)
    $obs_days = floor($obs_time / 86400);
    $obs_secs = $obs_time - ($obs_days * 86400);
    print "\tObservation time is $obs_days days $obs_secs seconds\n";

    # Write out everything to the output file
    print OUTPUT "0\n";                   # Enter another observation
    print OUTPUT "$dart_type\n";          # Variable kind
    print OUTPUT "2\n";                   # Pressure coordinates
    print OUTPUT "$obs_fields[10]\n";     # Pressure level
    print OUTPUT "$obs_fields[9]\n";      # Longitude
    print OUTPUT "$obs_fields[8]\n";      # Latitude
    print OUTPUT "$obs_days $obs_secs\n"; # Observation time
    print OUTPUT "$obs_fields[6]\n";      # Error variance
    print OUTPUT "$obs_fields[2]\n";      # Actual observation
}

# There are no more observations - finish up
print "\nProcessing footer...\n";
print OUTPUT "-1\n"; 
print OUTPUT "set_def.out\n";

close INNOV;
close OUTPUT;

print "Done writing $output_file_name\n";

# navdas_to_obs_kind
# ------------------
# Provide translation between the NAVDAS observation keys and the DART
# obs_kind keys for the various types of observations.  Takes 1
# argument - the NAVDAS variable type.  Returns the DART number or
# -999 if it is "not found"
sub navdas_to_obs_kind
{
    $navdas_kind = $_[0];
    $obs_name = "";

    # Iterate through the hash to find the variable name we want
    while (($key,$val) = each(%navdas_hash))
    {
    if ($val == $navdas_kind)
    {
        $obs_name = $key;
        print "\tIdentified observation of type $obs_name\n";
    }
    }

    # Grab the DART observation number and return it
    $dart_num = $dart_hash{$obs_name}?$dart_hash{$obs_name}:-999;
    return $dart_num;
}

# build_navdas_hash
# -----------------
# Build a hash table with the keys for the NAVDAS variable numbers
sub build_navdas_hash
{
    $navdas_hash{'UWND'} = 3;   # Zonal wind (m/s)
    $navdas_hash{'VWND'} = 4;   # Meridional wind (m/s)
    $navdas_hash{'PSUR'} = 10;  # Surface pressure (mb)
    $navdas_hash{'TEMP'} = 2;   # Air temperature (K)
    $navdas_hash{'SHUM'} = 15;  # Specific humidity (kg/kg)
    $navdas_hash{'PRES'} = 11;  # Station pressure (mb)
    $navdas_hash{'HGHT'} = 1;   # Geopotential height (m)
    $navdas_hash{'TDPD'} = 12;  # Dewpoint depression (K)
    $navdas_hash{'POTT'} = 16;  # Potential temperature (K)
}

# build_dart_hash
# -----------------
# Build a hash table with the keys for the DART obs_kind numbers
sub build_dart_hash
{
    $dart_hash{'UWND'} = 1;    # Zonal wind (m/s)
    $dart_hash{'VWND'} = 2;    # Meridional wind (m/s)
    $dart_hash{'PSUR'} = 3;    # Surface pressure (mb)
    $dart_hash{'TEMP'} = 4;    # Air temperature (K)
    $dart_hash{'SHUM'} = 5;    # Specific humidity (kg/kg)
    $dart_hash{'PRES'} = 6;    # Station pressure (mb)
    $dart_hash{'HGHT'} = 206;  # Geopotential height (m)
    $dart_hash{'TDPD'} = 205;  # Dewpoint depression (K)
    $dart_hash{'POTT'} = 207;  # Potential temperature (K)
}

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

