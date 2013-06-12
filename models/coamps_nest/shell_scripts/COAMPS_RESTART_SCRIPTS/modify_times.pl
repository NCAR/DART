#!/usr/bin/perl
#
# This code may (or may not) be part of the COAMPS distribution,
# So it is not protected by the DART copyright agreement.
#
# DART $Id$
#
########################################################################
#
# SCRIPT:   modify_times.pl
# AUTHOR:   T. R. Whitcomb
#           Naval Research Laboratory
#
# Accepts two file names - a namelist file to adjust the times in, and
# a file containing the new start and end times as a single comma-separated
# list
######

#
# Process the forecast times
#
open (TIMES, "<$ARGV[0]");
my $dart_times = <TIMES>;
chomp $dart_times;
my @tau_times = split(/\s*,\s*/, $dart_times);
my $tau_start = join(', ', @tau_times[0..2]);
my $tau_end   = join(', ', @tau_times[3..5]);
print "Advancing COAMPS Forecast:\n";
print "$tau_start ---> $tau_end\n";
close(TIMES);

#
# Modify the COAMPS namelist
#
open (NAMELIST, "<$ARGV[1]");
my @entire_namelist = <NAMELIST>;
close(NAMELIST);
my $namelist_data = join("", @entire_namelist);

# The start time is the current time
$namelist_data =~ s/(\s*ktaust\s*=)\s*(?:\d+\s*,[\s\n]*)+\n/\1 $tau_start,\n/;

# The "Destination time" is whatever DART says it is - needs one HMS triad for *each* nest
# Adding an extra spacing at the beginning of the line is unnecessary, but it makes it prettier
$namelist_data =~ m/nnest\s*=\s*(\d+)/;
my $num_nests  = $1;
my $nest_ktauf =  (" " x 10 . "$tau_end,\n") x ($num_nests-1);
$namelist_data =~ s/(\s*ktauf\s*=)\s*(?:\d+\s*,[\s\n]*)+\n/\1 $tau_end,\n$nest_ktauf/;

# Ensure that a restart file will be written by setting it to write out at the target time
$namelist_data =~ s/(\s*ksavea\s*=)\s*(?:\d+\s*,[\s\n]*)+\n/\1 $tau_end\n/;

# Do an in-place edit
open(NAMELIST, ">$ARGV[1]");
print NAMELIST $namelist_data;
close(NAMELIST);

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

