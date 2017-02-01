#!/usr/bin/perl
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

use File::Copy;

my $ensemble_size = 100;
my $cutoff = 400.0;              #cutoff radius for Schur product
my $cov_inflate = 1.1;           #covariance inflation factor
my $num_groups = 10;             #number of groups
my $first_obs = "0 43200";       #days, seconds
my $obs_period = "0 43200";      #days, seconds
my $num_periods = "24";          #number of obs times

#my $startobs = 1;                #first obs in set of X
#my $nobs = 36;                   #identity for X
#my $error_variance = 1.0;        #better for X
my $startobs = 37;               #first obs in set of Y
my $nobs = 360;                  #identity for Y
my $error_variance = 0.1;        #better for Y

my $max_obs = $nobs * $num_periods;

# create a set def
my $command = "mkmf_create_obs_sequence\n";
system "$command";
copy "input.nml.create_obs_sequence_default", "input.nml";

open INFILE, ">tmp";
print INFILE "$max_obs\n";                  #max obs
print INFILE "0\n";                         #0 copies for def
print INFILE "0\n";                         #no QC
for $iob ( 1 .. $nobs ) {
  print INFILE "$iob\n";                    #denotes identity
  $thisob = -($iob + $startobs - 1);        #identity location
  print INFILE "$thisob\n";                 #denotes identity
  print INFILE "0 0\n";                     #secs days
  print INFILE "$error_variance\n";         #error variance
}
print INFILE "-1\n";                       #done
print INFILE "obs_seq_def.out\n";          #def file name
close INFILE;

$command = "create_obs_sequence < tmp\n";
system "$command";
unlink "tmp";

# create an identity obs sequence
# propagate through times
my $command = "mkmf_create_fixed_network_seq\n";
system "$command";
copy "input.nml.create_fixed_network_seq_default", "input.nml";

open INFILE, ">tmp";
print INFILE "obs_seq_def.out\n";
print INFILE "1\n";                             #flag regular interval
print INFILE "$num_periods\n";                  
print INFILE "$first_obs\n";                    #time of initial ob 
print INFILE "$obs_period\n";                   #time between obs 
print INFILE "obs_seq.in\n";                    #output file
close INFILE;

$command = "create_fixed_network_seq < tmp \n";
system "$command\n";
unlink "tmp";

# create the obs
$command = "mkmf_perfect_model_obs";
system "$command\n";
my $template = "input.nml.perfect_model_obs_default";
open INFILE, $template;
my $nlfile = "input.nml";
open OUTFILE, ">$nlfile";
while (<INFILE>) {
  s/start_from_restart\s+=\s+.\w.+/start_from_restart = .true./;
  s/output_restart\s+=\s+.\w.+/output_restart = .true./;
  s/restart_in_file_name\s+=\s+"\w+"/restart_in_file_name = \"perfect_ics\"/;
  s/restart_out_file_name\s+=\s+"\w+"/restart_out_file_name = \"perfect_restart\"/;
  print OUTFILE;
}
close OUTFILE;
close INFILE;

my $command = "perfect_model_obs";
system "$command";

#run the filter
$command = "mkmf_filter";
system "$command\n";
$template = "input.nml.filter_default";
open INFILE, $template;
$nlfile = "input.nml";
open OUTFILE, ">$nlfile";
while (<INFILE>) {
  s/ens_size\s+=\s+\w+/ens_size = $ensemble_size/;
  s/cutoff\s+=\s+\w+.\w+/cutoff = $cutoff/;
  s/cov_inflate\s+=\s+\w+.\w+/cov_inflate = $cov_inflate/;
  s/start_from_restart\s+=\s+.\w.+/start_from_restart = .true./;
  s/output_restart\s+=\s+.\w.+/output_restart = .true./;
  s/restart_in_file_name\s+=\s+"\w+"/restart_in_file_name = \"filter_ics\"/;
  s/restart_out_file_name\s+=\s+"\w+"/restart_out_file_name = \"filter_restart\"/;
  s/num_output_state_members\s+=\s+\w+/num_output_state_members = $ensemble_size/;
  s/num_groups\s+=\s+\w+/num_groups = $num_groups/;
  print OUTFILE;
}
close OUTFILE;
close INFILE;

my $command = "filter";
system "$command";

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Date$
# $Revision$

