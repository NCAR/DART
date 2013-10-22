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
# Reads in a user-specified data file containing the details of the
# variables that will go in the state.vars file.  This file:  
#  - accepts comments, signified by a line starting with #
#  - has one variable per line
#  - has information separated by whitespace
#
# The command line argument specifies the name of the data file.
######

$DEBUG_FLAG = 0;
$flat_files = true;

# Hold all data lines from the file
my @AllData;

# Check and make sure that the inputs are proper
print "Usage: populate_state_vars state_data_file flat_files\n" and die
    unless $#ARGV >=0;
$state_data_file_name = $ARGV[0];
$flat_files = $ARGV[1];

$flat_file_io=$flat_files =~ m/true/i;
print "flat_files = $flat_files flat_file_io = $flat_file_io\n";
print "  Using data file $state_data_file_name...\n";
if (! -e "$state_data_file_name")
{
    print "  Argument specifies a non-existent file!\n";
    exit(1);
}
print "\n";


# Open up the files that we need
$state_vars_file_name = "state.vars";
die "Couldn't open file $state_data_file_name\n" unless 
    open(DAT_FILE,"<$state_data_file_name\n");
die "Couldn't open file $state_vars_file_name\n" unless 
    open(VAR_FILE,">$state_vars_file_name\n");

# Keep track of how many entries will go in the state.vars file
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
    if($flat_file_io)
	{
		if ($temp_dat[7] =~ m/true/i)
		{
		print "Skipping mean field...\n" if $DEBUG_FLAG;
		next;
		}
    	$num_entries += 1;
	}
    else
	{
    	# Add the number of levels for this variable to the total count of
    	# lines that will go into the state.vars file
    	$num_entries += $temp_dat[1];
    }
    print "Read in $temp_dat[1] levels for $temp_dat[0]...\n";

    # Add this line to our array and inform the user what we read
    push @AllData,$_;

}
print "-"x50 . "\n";
print "Will process $num_entries entries into state.vars\n";


### Output to state.vars

# Put in the total number of variables that we will write
if($flat_file_io)
{
  print VAR_FILE "FLAT_FILE\n";
}
else
{
  print VAR_FILE "RESTART_FILE\n";
}

# Now, loop through the array that we built from the data file and put
# each level in the state.vars file
foreach(@AllData)
{
    @cur_data = split;
    $levels = $cur_data[1];
    print "Writing $levels levels for $cur_data[0]...\n";
	if($flat_file_io)
	{
    $data_string = join(' ',@cur_data);
    print VAR_FILE "$data_string\n";
	}
    else
	{
    for($ii = 1; $ii<=$levels; $ii++)
    {
    # Set the current level
    $cur_data[1] = sprintf('%03u',$ii);
    
    $data_string = join(' ',@cur_data);
    print VAR_FILE "$data_string\n";
    }
	}

}
print "-"x50 . "\n";
close(DAT_FILE);
close(VAR_FILE)

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

