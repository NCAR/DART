#!/usr/bin/perl
#
# This code may (or may not) be part of the COAMPS distribution,
# So it is not protected by the DART copyright agreement.
#
# DART $Id$
#
########################################################################
#
# SCRIPT:   parse_field_list.pl
# AUTHOR:   T. R. Whitcomb
#           Naval Research Laboratory
#
# Takes a list of fields that have been pulled out of the COAMPS 
# restart file writing routines and calculates the information that
# DART needs to do the conversion to a state vector.  These are:
#  1.  The field's dimension (2d/3d on mass or w sigma levels)
#  2.  The field's position in listings of that type
#  3.  Whether or not an I/O processor is used.
#
# Command line arguments are the first fields in each of the three
# dimension types, the last field to process, the type of I/O (this
# should be either SINGLEIO or MULTIIO) and the file that contains
# the field listing (the easiest way to do this is generate it with
# the generate_restart_field_list.sh script).
######

$start_2D        = $ARGV[0];
$start_3D        = $ARGV[1];
$start_3DW       = $ARGV[2];
$last_field      = $ARGV[3];
$io_type         = $ARGV[4];
$field_list_file = $ARGV[5];

@dim_types    = ( 'DIM_TYPE_2D', 'DIM_TYPE_3D', 'DIM_TYPE_3DW' );
@num_vars     = ( 0            , 0            , 0              );
@first_fields = ( $start_2D    , $start_3D    , $start_3DW     );

$current_dim_type  = -1;

open(FIELD_LIST,"<$field_list_file");

while (<FIELD_LIST>)
{
        chomp;

        if (m/adom\(\w+\)%(\w+)\s*\(/)
        {
                $current_var_name = $1;
        
        if ($current_var_name eq $first_fields[$current_dim_type+1])
        {
            $current_dim_type++;
        }

        # Now actually process
        print "$current_var_name\t";
        print "$dim_types[$current_dim_type]\t";
        print ++$num_vars[$current_dim_type];
        print "\t$io_type";
        print "\n";

        if ($current_var_name eq $last_field)
        {
            last;
        }
        }
}

close(FIELD_LIST)

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

