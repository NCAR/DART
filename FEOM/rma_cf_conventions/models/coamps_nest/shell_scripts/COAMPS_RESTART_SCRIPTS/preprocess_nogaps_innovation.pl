#!/usr/bin/perl
#
# This code may (or may not) be part of the COAMPS distribution,
# So it is not protected by the DART copyright agreement.
#
# DART $Id$
#
########################################################################
#
# SCRIPT:   preprocess_nogaps_innovation.pl
# AUTHOR:   T. R. Whitcomb
#           Naval Research Laboratory
#
# Reads in a NOGAPS innovation file and preprocesses it to handle 
# conversion to the DART observation file set.  This means that
# we will:
#  1. Remove observations outside the domain of interest
#  2. Remove observations based on QC flags
#  3. Remove observations we can't assimilate (like brightness temp)
#
# The command line arguments specify the name of the innovation file
# and the name of the output file.
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


package main;

# Use DateTime library - may have to install from CPAN
# These paths are for NRL Hovermale cluster
use lib '/users/whitcomb/tools/modules';
use DateTime;
use DateTime::Duration;

print "Usage: preprocess_nogaps_innovation.pl innovation_file output_file\n"
  and die unless $#ARGV>=1;

$innov_file_name  = $ARGV[0];
$output_file_name = $ARGV[1];

# Useful constants
$MINUTES_TO_SECONDS = 60;
$HOURS_TO_SECONDS   = 60 * $MINUTES_TO_SECONDS;
$DAYS_TO_SECONDS    = 24 * $HOURS_TO_SECONDS;

# Set up the base DTG (for testing - this will be a parameter later)
$baseDtg = DateTime->new( year  => 2005,
              month => 06,
              day   => 24,
              hour  => 00
            );

# Set the location of the datahd File (for testing - this will be a
# paramter later) 
$datahdfile = 'datahd_sfc_000000_000000_1a2000x0001_2005080418_00000000_infofld';
$datahdout = 'domain.dat';
$check_in_grid_command = '/users/whitcomb/DART/models/coamps/work/check_in_grid';

# Open the files
die "Couldn't open $datahdfile!\n" unless 
    open(DATAHD,"<$datahdfile");
die "Couldn't open $innovation_file_name!\n" unless 
    open(INNOV,"<$innov_file_name");
die "Couldn't open $output_file_name!\n" unless 
    open(OUTPUT,">$output_file_name");

print "Innovation file $innov_file_name opened successfully.\n";
print "Output file $output_file_name opened successfully.\n";

# Read the grid information from the datahd file
@gridfull = <DATAHD>;
$gridlong = join(' ',@gridfull);
@gridfields = split(/\s+/,$gridlong);
@gridinfo = @gridfields[3..10,30..35];

# Read the header information - what we really want is the DTG
# information 
print "\nProcessing header...\n";
do {$line = <INNOV>} until $line =~ m/number of obs/;
chomp($line);
$line =~ m/\D*(\d+)\D*(\d+)\D*(\d+)\D*/;
$n_obs = $1;
$cdtg  = $2;
$tau   = $3;
print "Number of observations in sequence : $n_obs\n";
print "Date-time group of forecast cycle  : $cdtg\n";
print "Lead time of current analysis      : $tau\n";

# Convert the DTG to a Perl DateTime object
# Dtg is in the form YYYYMMDDHH
$cycleDtg = DateTime->new( year  => substr($cdtg,0,4),
               month => substr($cdtg,4,2),
               day   => substr($cdtg,6,2),
               hour  => substr($cdtg,8,2)
             );
$lead = DateTime::Duration->new (hours=>$tau);
$analDtg = $cycleDtg + $lead;
print "Analysis date-time-group           : " .
  $analDtg->strftime('%Y%m%d%H') . "\n";
print "Base date-time-group               : " .
  $baseDtg->strftime('%Y%m%d%H') . "\n";

# Skip the column labels
<INNOV>;

# The FORTRAN format string to write this out is
# i7,5(1x,f8.2),3(1x,f9.2),1x,f8.2,2(1x,i4),
# 1x,i5,1x,i4,1x,i7,1x,a17,1x,a11,1x,i4
$template  = 'A7 ';
$template .= 'A9 ' x 5;
$template .= 'A10 ' x 3;
$template .= 'A9 ';
$template .= 'A5 ' x 2;
$template .= 'A6 A5 A8 A18 A12 A5';

# Define the maximum and minimum lon/lat pairs that we'll consider
# checking - these not be exact for the grid since we'll do another
# test on the points within this rectangle to make sure that they are
# on the grid
$minlat = 9;
$maxlat = 61;
$minlon = 205;
$maxlon = 324;

# Initialize the counters that keep track of rejected obs
$sat_skipped = 0;
$lon_skipped = 0;
$lat_skipped = 0;
$out_skipped = 0;
$qc_skipped  = 0;

# Initialize the counters for the various types
@typecount = (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);

# Initialize the progress spinner
my $spinner = Spinner->new;

print "\nProcessing observations...\n";
while (<INNOV>)
{
  chomp;
  @innov_data = unpack($template,$_);

  $spinner->spin;

  # DART cares about:
  # 0. Observation Value
  # 1. Observation Error
  # 2. Observation Latitude
  # 3. Observation Longitude
  # 4. Observation Pressure
  # 5. Observation Type
  # 6. Observation QC flag
  # 7. Observation Time offset
  @dart_data  = @innov_data[1,5,7..10,13..14];

  # Now, we have a bunch of checks for what should get included and
  # what should get tossed out

  # Don't use if if this is a satellite obs
  if ($dart_data[5] == 13)
    {
      $sat_skipped++;
      next;
    }

  # Check if this is a ridiculous value of lat/lon for this run
  if ($dart_data[2] < $minlat or $dart_data[2] > $maxlat)
    {
      $lat_skipped++;
      next;
    }
  if ($dart_data[3] < $minlon or $dart_data[3] > $maxlon)
    {
      $lon_skipped++;
      next;
    }

  # Check QC flag
  if ($dart_data[6] != 0)
      {
    $qc_skipped++;
    next;
      }

  # Check if this value of lat/lon is within the domain - this is done
  # by calling a FORTRAN whose specific job is to return T or F based
  # on the grid values.  However, this should really be implemented in
  # this script to speed this computation up as it requires an
  # open/write/close EVERY TIME and is *very* slow.
  die "Couldn't open $datahdout!\n" unless 
    open(GRIDOUT,">$datahdout");
  print GRIDOUT join("\n",@gridinfo);
  print GRIDOUT $dart_data[2]."\n";
  print GRIDOUT $dart_data[3]."\n";
  close(GRIDOUT);
  $in_grid_flag = `$check_in_grid_command`;
  if ($in_grid_flag =~ /F/i)
    {
      $out_skipped++;
      next;
    }

  # If we got here, this is an ob we need to process

  ### TIME COMPUTATIONS
  # The actual observation time
  $dt    = DateTime::Duration->new(seconds => $dart_data[7]);
  $obDtg = $analDtg + $dt;

  # Figure the actual observation time in seconds relative to our base
  # date-time group (which is static)
  $baseObOffset = $obDtg - $baseDtg;
  @offsetData   = $baseObOffset->in_units('days','hours','minutes',
                      'seconds');
  $offset_seconds = ($offsetData[0] * $DAYS_TO_SECONDS) +
                    ($offsetData[1] * $HOURS_TO_SECONDS) +
                    ($offsetData[2] * $MINUTES_TO_SECONDS) +
                ($offsetData[3]);

  # Save it back to our array
  $dart_data[7] = $offset_seconds;

  # Increment the obs type counter
  $typecount[$dart_data[5]-1]++;

  print OUTPUT "@dart_data\n";
}

print "Satellite points rejected          : $sat_skipped\n";

$out_coarse_skipped = $lat_skipped + $lon_skipped;
print "Gross outside points rejected      : $out_coarse_skipped\n";
print "Quality control points rejected    : $qc_skipped\n";
print "Fine outside points rejected       : $out_skipped\n";

$obs_remain = $n_obs - $sat_skipped - $out_coarse_skipped;
print "\nObservations remaining             : $obs_remain\n";


$ii = 1;
print "\nObservation counts";
print "\n------------------\n";
foreach (@typecount)
  {

    print "Count for obs $ii \t: $_\n";
    $ii++;
  }
# Cleanup
close DATAHD;
close INNOV;
close OUTPUT;

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

