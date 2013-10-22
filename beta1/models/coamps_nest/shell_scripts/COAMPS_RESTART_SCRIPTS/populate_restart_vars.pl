#!/usr/bin/perl
#
# This code may (or may not) be part of the COAMPS distribution,
# So it is not protected by the DART copyright agreement.
#
# DART $Id$
#
########################################################################
#
# SCRIPT:   populate_restart_vars.pl
# AUTHOR:   T. R. Whitcomb
#           Naval Research Laboratory
#
# Reads in a user-specified data file containing the details of the
# variables that will go in the restart.vars file.  This file:  
#  - accepts comments, signified by a line starting with #
#  - has one variable per line
#  - has information separated by whitespace
#
# The command line argument specifies the name of the data file.
######

$DEBUG_FLAG = 0;

# Hold all data lines from the file
my @AllData;

# Check and make sure that the inputs are proper
print "Usage: populate_restart_vars restart_data_file\n" and die
    unless $#ARGV >=0;
$restart_data_file_name = $ARGV[0];
print "  Using data file $restart_data_file_name...\n";
if (! -e "$restart_data_file_name")
{
    print "  Argument specifies a non-existent file!\n";
    exit(1);
}
print "\n";

# Open up the files that we need
$restart_vars_file_name = "restart.vars";
die "Couldn't open file $restart_data_file_name\n" unless 
    open(DAT_FILE,"<$restart_data_file_name\n");
die "Couldn't open file $restart_vars_file_name\n" unless 
    open(VAR_FILE,">$restart_vars_file_name\n");

# Keep track of how many entries will go in the restart.vars file
$num_entries = 0;

print "-"x50 . "\n";
# Read in each line of the data file
while(<DAT_FILE>)
{
    # Skip comment lines and blank lines
    if (m/^\#/)
    {
    print "Skipping comment line...\n" if $DEBUG_FLAG;
    next;
    }
    if (m/^\s*$/)
    {
    print "Skipping blank line...\n" if $DEBUG_FLAG;
    next;
    }

    # Break up this line by whitespace
    @temp_dat = split;
    print "Read in $temp_dat[1] levels for $temp_dat[0]...\n";
    
    # Add the number of levels for this variable to the total count of
    # lines that will go into the restart.vars file
    $num_entries += $temp_dat[1];

    # Add this line to our array and inform the user what we read
    push @AllData,$_;

}
print "-"x50 . "\n";
print "Will process $num_entries entries into restart.vars\n";


### Output to restart.vars

# Put in the total number of variables that we will write
print VAR_FILE "$num_entries\n";

# Now, loop through the array that we built from the data file and put
# each level in the restart.vars file
foreach(@AllData)
{
    @cur_data = split;
    $levels = $cur_data[1];
    print "Writing $levels levels for $cur_data[0]...\n";
    for($ii = 1; $ii<=$levels; $ii++)
    {
    # Set the current level
    $cur_data[1] = sprintf('%03u',$ii);
    
    $data_string = join(' ',@cur_data);
    print VAR_FILE "$data_string\n";
    }

}
print "-"x50 . "\n";
close(DAT_FILE);
close(VAR_FILE)

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

