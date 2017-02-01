#!/usr/bin/perl
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# This perl script generates a file that can be redirected into
# create_obs_seequence.  It requires only two pieces of information 
# from STDIN.  As it currently stands, it assumes that the output 
# filename is set_def.out, and it assumes a single uniqe observation 
# set.
#
# The user inputs are the spacing between obs (the script figures
# out the resulting observation locations under the assumption
# that the first observation is at location 1) and the observational
# variance (the script assumes all observations have the same error
# variance).
#
# After running this script, the user runs the command:
#
#    ./create_obs_sequence < coseq.dat

# open output file
open(OUT,'>coseq.dat');

# Output a large number for memory allocation purposes
print OUT "10000000\n";

# Output the number of copies of data.  Since only creating
# a definition here, output 0
print OUT "0\n";

# Output the number of quality control values per field.  There
# is not QC used here.
print OUT "0\n";

# Obtain spacing between observations, assuming you start at point 1
print "Will assume that the first observation is at x(1).\n";
print "Input the subsequent spacing between obs: ";
$spaceobs=<STDIN>;

# Figure out how many observation locations this implies
$a=1;
$numobs=0;
while($a+$spaceobs < 960) {
 $a=$a+$spaceobs;
 $numobs++;
}
$numobs=$numobs;
print "Resulting number of obs is: ";
print "$numobs\n";

# Obtain error variances
print "Will assume that all observations have the same error variance.\n";
print "Input this error variance: ";
$var=<STDIN>;

# Put variance and observation location information
# into output file
$loc=1;
print OUT "0\n";		# enter another ob info
print OUT "-$loc\n";		# identity ob at $loc
print OUT "0, 0\n";		# time (meaningless here)
print OUT $var;			# error variance
for ($i==0; $i < $numobs; $i++) {
 print OUT "0\n";		# enter another ob info
 $loc = $loc + $spaceobs;
 print OUT "-$loc\n";		# identity ob at $loc
 print OUT "0, 0\n";		# time (meaningless here)
 print OUT $var;		# error variance
}

# Output that we are done entering obs info
print OUT "-1\n";

# Output the input file name
print OUT "set_def.out\n";

# Clean up
print "Done.\n";
print "\n";
print "Note, if you are entering huge numbers of observations, it might be\n";
print "necessary to enter coseq.pl and alter the first print OUT statement.\n";
close(OUT);

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

