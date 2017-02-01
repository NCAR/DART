#!/usr/bin/perl
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# This perl script generates a file that can be redirected into
# create_obs_seequence.  For the MITgcm annulus the script 
# requires 10 pieces of information from STDIN.  As the MITgcm
# configuration currently stands, there are 5 different variables
# that are assimilated, and the user is free to specify the 
# observational density and error covariance for each of the five
# variables independently.  The script assumes that the output
# filename is set_def.out, and it assumes a single uniqe observation 
# set.
#
# The user inputs are the spacing between observations for each
# of the five variable types (the script figures
# out the resulting observation locations under the assumption
# that the first observation is at location 1) and the observational
# variance for each of the five variables (the script assumes all 
# observations of that variable have the same error variance).
#
# After running this script, the user runs the command:
#
#    ./create_obs_sequence < coseq.dat

# specify the dimensions in the different directions
$naz=120;
$nrad=31;
$nzed=29;

# specify the number of elements per variable = naz*nrad*nzed
$statesize = $naz*$nrad*$nzed;

# construct a mask for where the fluid is!
$icount=0;
for ($k = 1; $k <= $nzed; $k++) {
 for ($j = 1; $j <= 8; $j++) {
  for ($i = 1; $i <= $naz;  $i++) {
   $icount++;
   $mask[$icount]=0;
  }
 }
 for ($j = 9; $j <= $nrad; $j++) {
  for ($i = 1; $i <= $naz;  $i++) {
   $icount++;
   $mask[$icount]=1;
  }
 }
}

# open output file
open(OUT,'>coseq.dat');

# Output a large number for memory allocation purposes.
# If the subsequent running of ./create_obs_sequence doesn't
# end with the line:
# write_obs_seq  opening formatted file set_def.out
# then this number must be made larger.
print OUT "100000\n";

# Output the number of copies of data.  Since only creating
# a definition here, output 0
print OUT "0\n";

# Output the number of quality control values per field.  There
# is not QC used here.
print OUT "0\n";

# Obtain spacing between observations, assuming you start at point 1
print "Enter observation information for the U variable.\n";
print "This script assumes that the first observation is at U(1).\n";
print "Input the subsequent spacing between obs: ";
$spaceobs=<STDIN>;

# Figure out how many U observation locations this implies
$a=1;
$numobs=0;
while($a + $spaceobs < $statesize) {
 $a = $a + $spaceobs;
 $numobs++;
}
$numobsout=$numobs+1;
print "Resulting number of obs is: ";
print "$numobsout.\n";

# Obtain error variances
print "This script assumes that all U observations have the same error variance.\n";
print "Input this observational error variance: ";
$var=<STDIN>;

# Put variance and observation location information
# into output file
$loc=1;
if ($mask[$loc] > 0) {
 print OUT "0\n";		# enter another ob info
 print OUT "-$loc\n";		# identity ob at $loc
 print OUT "0, 0\n";		# time (meaningless here)
 print OUT $var;		# error variance
}
for ($i=0; $i < $numobs; $i++) {
 $loc = $loc + $spaceobs;
 if ($mask[$loc] > 0) {
  print OUT "0\n";		# enter another ob info
  print OUT "-$loc\n";		# identity ob at $loc
  print OUT "0, 0\n";		# time (meaningless here)
  print OUT $var;		# error variance
 }
}

#-------

# Obtain spacing between observations, assuming you start at point 1
print "Enter observation information for the V variable.\n";
print "This script assumes that the first observation is at V(1).\n";
print "Input the subsequent spacing between obs: ";
$spaceobs=<STDIN>;

# Figure out how many V observation locations this implies
$a=1;
$numobs=0;
while($a + $spaceobs < $statesize) {
 $a = $a + $spaceobs;
 $numobs++;
}
$numobsout=$numobs+1;
print "Resulting number of obs is: ";
print "$numobsout.\n";

# Obtain error variances
print "This script assumes that all V observations have the same error variance.\n";
print "Input this observational error variance: ";
$var=<STDIN>;

# Put variance and observation location information
# into output file
$loc = $statesize + 1;
if ($mask[$loc] > 0) {
 print OUT "0\n";		# enter another ob info
 print OUT "-$loc\n";		# identity ob at $loc
 print OUT "0, 0\n";		# time (meaningless here)
 print OUT $var;		# error variance
} 
for ($i=0; $i < $numobs; $i++) {
 $loc = $loc + $spaceobs;
 if ($mask[$loc] > 0) {
  print OUT "0\n";		# enter another ob info
  print OUT "-$loc\n";		# identity ob at $loc
  print OUT "0, 0\n";		# time (meaningless here)
  print OUT $var;		# error variance
 }
}

#-------

# Obtain spacing between observations, assuming you start at point 1
print "Enter observation information for the W variable.\n";
print "This script assumes that the first observation is at W(1).\n";
print "Input the subsequent spacing between obs: ";
$spaceobs=<STDIN>;

# Figure out how many W observation locations this implies
$a=1;
$numobs=0;
while($a + $spaceobs < $statesize) {
 $a = $a + $spaceobs;
 $numobs++;
}
$numobsout=$numobs+1;
print "Resulting number of obs is: ";
print "$numobsout.\n";

# Obtain error variances
print "This script assumes that all W observations have the same error variance.\n";
print "Input this observational error variance: ";
$var=<STDIN>;

# Put variance and observation location information
# into output file
$loc = 2*$statesize + 1;
if ($mask[$loc] > 0) {
 print OUT "0\n";		# enter another ob info
 print OUT "-$loc\n";		# identity ob at $loc
 print OUT "0, 0\n";		# time (meaningless here)
 print OUT $var;		# error variance
}
for ($i=0; $i < $numobs; $i++) {
 $loc = $loc + $spaceobs;
 if ($mask[$loc] > 0) {
  print OUT "0\n";		# enter another ob info
  print OUT "-$loc\n";		# identity ob at $loc
  print OUT "0, 0\n";		# time (meaningless here)
  print OUT $var;		# error variance
 }
}

#-------

# Obtain spacing between observations, assuming you start at point 1
print "Enter observation information for the T variable.\n";
print "This script assumes that the first observation is at T(1).\n";
print "Input the subsequent spacing between obs: ";
$spaceobs=<STDIN>;

# Figure out how many T observation locations this implies
$a=1;
$numobs=0;
while($a + $spaceobs < $statesize) {
 $a = $a + $spaceobs;
 $numobs++;
}
$numobsout=$numobs+1;
print "Resulting number of obs is: ";
print "$numobsout.\n";

# Obtain error variances
print "This script assumes that all T observations have the same error variance.\n";
print "Input this observational error variance: ";
$var=<STDIN>;

# Put variance and observation location information
# into output file
$loc = 3*$statesize + 1;
if ($mask[$loc] > 0) {
 print OUT "0\n";		# enter another ob info
 print OUT "-$loc\n";		# identity ob at $loc
 print OUT "0, 0\n";		# time (meaningless here)
 print OUT $var;		# error variance
}
for ($i=0; $i < $numobs; $i++) {
 $loc = $loc + $spaceobs;
 if ($mask[$loc] > 0) {
  print OUT "0\n";		# enter another ob info
  print OUT "-$loc\n";		# identity ob at $loc
  print OUT "0, 0\n";		# time (meaningless here)
  print OUT $var;		# error variance
 }
}

#-------

# Obtain spacing between observations, assuming you start at point 1
print "Enter observation information for the P variable.\n";
print "This script assumes that the first observation is at P(1).\n";
print "Input the subsequent spacing between obs: ";
$spaceobs=<STDIN>;

# Figure out how many P observation locations this implies
$a=1;
$numobs=0;
while($a + $spaceobs < $statesize) {
 $a = $a + $spaceobs;
 $numobs++;
}
$numobsout=$numobs+1;
print "Resulting number of obs is: ";
print "$numobsout.\n";

# Obtain error variances
print "This script assumes that all P observations have the same error variance.\n";
print "Input this observational error variance: ";
$var=<STDIN>;

# Put variance and observation location information
# into output file
$loc = 4*$statesize + 1;
if ($mask[$loc] > 0) {
 print OUT "0\n";		# enter another ob info
 print OUT "-$loc\n";		# identity ob at $loc
 print OUT "0, 0\n";		# time (meaningless here)
 print OUT $var;		# error variance
}
for ($i=0; $i < $numobs; $i++) {
 $loc = $loc + $spaceobs;
 if ($mask[$loc] > 0) {
  print OUT "0\n";		# enter another ob info
  print OUT "-$loc\n";		# identity ob at $loc
  print OUT "0, 0\n";		# time (meaningless here)
  print OUT $var;		# error variance
 }
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

