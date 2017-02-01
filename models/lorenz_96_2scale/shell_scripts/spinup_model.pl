#!/usr/bin/perl
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

use File::Copy;
use File::Path;

my $num_spinup_days = 1000;

my $ensemble_size = 100;
my $error_variance = 1.0;

my $command = "mkmf_create_obs_sequence\n";
system "$command";
copy "input.nml.create_obs_sequence_default", "input.nml";

# create a set def
open INFILE, ">tmp";
print INFILE "$num_spinup_days\n";          #max obs
print INFILE "0\n";                         #0 copies for def
print INFILE "0\n";                         #no QC
print INFILE "1\n";                         #only one obs
print INFILE "-1\n";                         #identity
print INFILE "0 0\n";                         #secs days
print INFILE "100000\n";                      #error variance
print INFILE "-1\n";                         #done
print INFILE "obs_seq_def.out\n";               #def file name
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
print INFILE "1\n";                             #regularly repeating times
print INFILE "$num_spinup_days\n";              #days
print INFILE "1,0\n";                           #time of initial ob (1 dy)
print INFILE "1,0\n";                           #time between obs (1/dy)
print INFILE "obs_seq.in\n";                    #output file
close INFILE;

my $command = "create_fixed_network_seq < tmp \n";
system "$command";
unlink "tmp";

# init the model onto the attrctor
$command = "mkmf_perfect_model_obs";
system "$command\n";
$template = "input.nml.perfect_model_obs_default";
open INFILE, $template;
my $nlfile = "input.nml";
open OUTFILE, ">$nlfile";
while (<INFILE>) {
  s/start_from_restart\s+=\s+.\w.+/start_from_restart = .false./;
  s/output_restart\s+=\s+.\w.+/output_restart = .true./;
  s/restart_in_file_name\s+=\s+"\w+"/restart_in_file_name = \"perfect_ics\"/;
  s/restart_out_file_name\s+=\s+"\w+"/restart_out_file_name = \"perfect_restart\"/;
  print OUTFILE;
}
close OUTFILE;
close INFILE;

my $command = "perfect_model_obs";
system "$command";

# generate a set of ICs
open INFILE, $template;
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

# save deterministic spun-up ICs 
copy "perfect_restart","perfect_ics";
# actually create the truth run
$command = "perfect_model_obs";
system "$command";

# Generate an ensemble
$command = "mkmf_filter";
system "$command\n";
$template = "input.nml.filter_default";
open INFILE, $template;
my $nlfile = "input.nml";
open OUTFILE, ">$nlfile";
while (<INFILE>) {
  s/ens_size\s+=\s+\w+/ens_size = $ensemble_size/;
  s/cutoff\s+=\s+\w+.\w+/cutoff = 0.0/;
  s/cov_inflate\s+=\s+\w+.\w+/cov_inflate = 1.0/;
  s/start_from_restart\s+=\s+.\w.+/start_from_restart = .false./;
  s/output_restart\s+=\s+.\w.+/output_restart = .true./;
  s/restart_in_file_name\s+=\s+"\w+"/restart_in_file_name = \"perfect_ics\"/;
  s/restart_out_file_name\s+=\s+"\w+"/restart_out_file_name = \"filter_restart\"/;
  s/num_output_state_members\s+=\s+\w+/num_output_state_members = $ensemble_size/;
  print OUTFILE;
}
close OUTFILE;
close INFILE;

# spin up
$command = "filter";
system "$command";
# this will be the ICs for the truth run
copy "perfect_restart","perfect_ics";
# these are the ICs to use in the ensemble
copy "filter_restart","filter_ics";

# save the diagnostic files so we can see the spinup
mkpath (["spinup"]);
copy "Prior_Diag.nc","spinup/Prior_Diag.nc";
copy "Posterior_Diag.nc","spinup/Posterior_Diag.nc";
copy "Truth.nc","spinup/Truth.nc";

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

