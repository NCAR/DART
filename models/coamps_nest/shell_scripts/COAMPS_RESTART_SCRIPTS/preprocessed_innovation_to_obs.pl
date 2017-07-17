#!/usr/bin/perl
#
# This code may (or may not) be part of the COAMPS distribution,
# So it is not protected by the DART copyright agreement.
#
# DART $Id$
#
########################################################################
#
# SCRIPT:	preprocessed_innovation_to_obs.pl
# AUTHOR:	T. R. Whitcomb
#
# Script to handle the output of the innovation file preprocessor.
# This means the file that we deal with here has had all extraneous
# obs removed and contains only the things within the domain that we
# are going to assimilate.  In the future, we will just preprocess the
# obs files from both COAMPS and NOGAPS and run them through this
# script. 
#
# The values read in from this file are:
# 1. Observation Value
# 2. Observation Error
# 3. Observation Latitude
# 4. Observation Longitude
# 5. Observation Pressure
# 6. Observation Type
# 7. Observation QC flag
# 8. Observation Time offset
#
# The only command line argument for this file is the path to the
# preprocessed innovation file.
######

# This is a package to do the \|/-\ spinner to keep the user aware
# that we are working.
# From: http://home.san.rr.com/lanois/perl/20010201_000000.html

package Spinner;

sub new
{
    my ($class) = @_;

    bless {
        _position => 0,
        _picture  => [ '|', '/', '-', '\\', '|', '/', '-', '\\' ]
        }, $class;
}

sub spin()
{
    $_[0]->{_position} = ($_[0]->{_position} + 1) % 8;
    print $_[0]->{_picture}[$_[0]->{_position}];
    print "\r";
}

### Main program
package main;
my $spinner = Spinner->new;

print "Usage: preprocessed innovation_to_obs.pl preprocessed_file" . 
  " obs_file\n" and die unless $#ARGV>=1;

$innov_file_name  = $ARGV[0];
$obs_file_name = $ARGV[1];

die "Couldn't open $innovation_file_name!\n" unless 
    open(INNOV,"<$innov_file_name");
die "Couldn't open $output_file_name!\n" unless 
    open(OBS,">$obs_file_name");

@total_obs = split(/\s+/, `wc -l $innov_file_name`);

print "Dealing with a total of $total_obs[0] obs...\n";

$total_obs[0]++;
# Write out the boilerplate for create_obs_sequence
print OBS "$total_obs[0]\n";    # Upper bound on obs
print OBS "1\n";         # Set 1 copy of the data
print OBS "0\n";         # no QC fields
print OBS "INOV FILE\n"; # metadata for sequence

&build_type_array;
while (<INNOV>)
{
  chomp;
  @obs_data = split(/\s+/);
  $spinner->spin;

  # Write out everything to the output file
  print OBS "0\n";		      # Enter another observation
  print OBS "$types[$obs_data[6]]\n"; # Variable kind
  print OBS "2\n";		      # Pressure coordinates
  print OBS "$obs_data[5]\n";         # Pressure level
  print OBS "$obs_data[4]\n";         # Longitude
  print OBS "$obs_data[3]\n";         # Latitude
  print OBS "0 $obs_data[8]\n";       # Observation time
  print OBS "$obs_data[2]\n";         # Error variance
  print OBS "$obs_data[1]\n";         # Actual observation
}

print OBS "-1\n"; 
print OBS "navdas_obs.out\n";

close(INNOV);
close(OBS);

# build_type_array
# ---------------- 
# Build this as an array that translates the numbered
# obs type into a string that DART knows
sub build_type_array
{
  $types[1]  = 'GEOPOTENTIAL_HEIGHT';
  $types[2]  = 'AIR_TEMPERATURE';
  $types[3]  = 'U_WIND_COMPONENT';
  $types[4]  = 'V_WIND_COMPONENT';
  $types[5]  = 'LOG_SPECIFIC_HUMIDITY';
  
  # There are also more that are not used at present
  $types[6]  = 'OZONE_MIXING_RATIO';
  $types[7]  = 'WIND_DIRECTION';
  $types[8]  = 'WIND_SPEED';
  $types[9]  = 'THICKNESS';
  $types[10] = 'MEAN_SEA_LEVEL_PRESSURE';
  $types[11] = 'STATION_PRESSURE';
  $types[12] = 'DEWPOINT_DEPRESSION';
  $types[13] = 'BRIGHTNESS_TEMPERATURE';
  $types[14] = 'TOTAL_PRECIPITABLE_WATER';
  $types[15] = 'SPECIFIC_HUMIDITY';
  $types[16] = 'POTENTIAL_TEMPERATURE';
}

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

