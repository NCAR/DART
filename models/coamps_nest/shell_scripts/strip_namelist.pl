#!/usr/bin/perl
#
# This code may (or may not) be part of the COAMPS distribution,
# So it is not protected by the DART copyright agreement.
#
# DART $Id$
#
########################################################################
#
# AUTHOR:   T. R. Whitcomb
#           Naval Research Laboratory
#
# Takes an existing (presumably large) namelist and pulls out a
# user-defined subset of fields to create a new (presumably smaller)
# namelist.  The name of the output file is derived from the name of
# the stripped namelist.
#
# The command line inputs are the full namelist file name, the name
# of the file that contains the variables to pull out, and the name
# that the new namelist should be given.
######

print "Usage: strip_namelist full_namelist var_file namelist_out\n"
    and die unless $#ARGV >= 2;

# Copy over input variables
$full_listfile = $ARGV[0];
$variable_file = $ARGV[1];
$stripped_list = $ARGV[2];
$stripped_file = "$stripped_list.nml";

# Check and make sure that our file names exist
print "Using output file $stripped_file...\n";
print "Checking if $full_listfile & $variable_file exist...\n";
if ((! -e $full_listfile) || (! -e $variable_file))
{
    print "Argument specifies a non-existent file!\n";
    exit(1);
}

# Open the files
die "Couldn't open $full_listfile\n" unless
    open(FULL_LIST,"<$full_listfile"); 
die "Couldn't open $variable_file\n" unless
    open(VAR_LIST,"<$variable_file");
die "Couldn't open $stripped_file\n" unless
    open(STRIP_LIST,">$stripped_file");

# Read in the variables we want and strip off extra characters
@vars = <VAR_LIST>;
chomp(@vars);

# Read in the full namelist file
@namelist = <FULL_LIST>;

# Start writing out our stripped namelist file
print STRIP_LIST "&$stripped_list\n";

# Loop through our array of variables that we want to read in, and
# match it to the full namelist file (if possible)
# Elements in a namelist file are NAME = VALUE,\n
foreach $var (@vars)
{
    # Define our separator character - Namelists use ,
    $sep_char = ',';
    foreach (@namelist)
    {
      m|($var)\s*=\s*(.*?)[/,]\s*$|;
    if ($1)
    {
        # Write matched variable to the namelist file and move on
        print " Matched $1\n";
        print STRIP_LIST "  $1 = $2$sep_char\n";
        last;
    }
    }
}

# Last character of a namelist should be a /
print STRIP_LIST "/\n";

# Close everything up
close(FULL_LIST);
close(VAR_LIST);
close(STRIP_LIST);

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

