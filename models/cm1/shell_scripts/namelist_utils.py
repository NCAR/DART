#!/usr/bin/env python
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# CREDIT: This script was donated to DART by Luke Madaus during his time
# at the University of Washington. Thanks Luke!
#-----------------------------------------------------------------------------

from __future__ import print_function, division
from collections import OrderedDict

def read_namelist(infname='namelist.input'):
    """
    This reads in the namelist file given by
    infname and returns a dictionary of
    the namelist contents organized by subsection
    and then variable name
    """
    nmld = OrderedDict()
    # Putting this in a "with" statement will close the
    # file when we are done
    with open(infname, 'r') as infile:
        for line in infile.readlines():
            cleanline = line.strip() # This is like Fortran 'trim'
            # Skip comment lines 
            if cleanline.startswith('!'):
                continue
            # Skip blank lines and section ends
            if cleanline in ['', '/']:
                continue
            # If we start with an ampersand, this is a section header
            if cleanline[0] == '&':
                sectionhead = cleanline[1:]
                # Make a section in the nmld
                nmld[sectionhead] = {}
                continue
            # If we're not blank, not a section end, and
            # we dont begin with '&', this is a variable
            # Split into to variables on equals sign
            try:
                varname, value = cleanline.split('=')
                varname = varname.strip()
                value = value.strip()
            except ValueError:
                # If the above split failed, we are still
                # in the previous variable
                # Hold over varname and append to current value
                if isinstance(nmld[sectionhead][varname], str):
                    oldvalue = nmld[sectionhead][varname]
                else:
                    oldvalue = ','.join([str(x) for x in nmld[sectionhead][varname]])
                value = ','.join([oldvalue,cleanline.strip()])
            # Remove comma at end, if it exists
            value = value.strip(',')
            # If there are other commas in this, it should
            # be stored as a list.  Trying to split
            # on commas will returna list
            valsplit = value.split(',')
            # Convert to correct types
            outvals = str_to_value(valsplit)
            # If this is length 1, don't store
            # the list
            if len(outvals) == 1:
                nmld[sectionhead][varname] = outvals[0]
            else:
                nmld[sectionhead][varname] = outvals
            # This will automatically continue to the next line here
    # Exit the with statement to close the file
    # And return the namelist dictionary
    return nmld
            
def str_to_value(instr):
    """
    Function to determine if the contents
    of a list or single string should stay a
    string, or be converted to a float or integer
    """
    if isinstance(instr, str):
        # Enclose this in a list for our usage
        instr = [strval]
    out = []
    for x in instr:
        try:
            # If this fails, this must be a string
            # and it will jump to the exception below
            float(x)
            # If there's a decimal point, store as float
            if '.' in x:
                out.append(float(x))
            else:
                out.append(int(x))

        except ValueError:
            # Remove any leading/trailing space and quotes
            strval = x.strip().strip("'").strip('"')
            # Check if it's a boolean
            if strval in ['.true.','T']:
                strval = True
            elif strval in ['.false.','T']:
                strval = False
            out.append(strval)
    return out
        

def format_single(v):
    """ 
    This function formats single values
    for the namelist writeout 
    """
    if isinstance(v, int):
        # Boolean is a subclass of int
        if isinstance(v, bool):
            if v:
                return ".true."
            else:
                return ".false."
        else:
            # Format as integer
            return "{:d}".format(v)
    elif isinstance(v, float):
        # Truncate at 4 decimals
        return "{:5.4f}".format(v)
    else:
        # Return the stringvalue with single quotes
        return "'{:s}'".format(v)

def var_format(invar):
    """ 
    Decide how to format this namelist value depending on what it is
    """
    # If it's a list, format each item individually and
    # join with commas
    if isinstance(invar, list):
        if len(invar) == 1:
            return format_single(invar[0])
        values = [format_single(v) for v in invar]
        """
        # Put on multiple lines if longer than a certain amount (5 here)
        if len(values) > 5:
            print(values)
            total_rows = int(len(values) / 5)
            allrows = []
            for r in xrange(total_rows):
                thisrow = values[r*5:(r+1)*5]
                rowjoin = ', '.join(thisrow)
                allrows.append(thisrow)
                print(allrows)
            return '\n                  '.join(allrows)
        else:
            return ', '.join(values)
        """ 
        return ', '.join(values)
    else:
        # Return the single format
        return format_single(invar)

def write_namelist(nmld, outfname='namelist.input'):
    """ 
    Formats a namelist from a dictionary (nmld)
    and writes it to file given by outfname
    """
    # Open the output file for writing
    # Set the newline character and a default indentation
    nl = '\n'
    indent = ' '
    with open(outfname, 'w') as outfile:
        print('', file=outfile)
        # Loop through all the sections of the namelist
        #section_names = nmld.keys()
        #section_names.sort()
        #section_names.reverse()
        for section_name, section in nmld.items():#section_names:
            #section = nmld[section_name]
            # Write the header
            #outfile.write(''.join((' &',section_name,nl)))
            print(''.join((' &',section_name)), file=outfile)
            # Now loop through all variables and write
            for varname, value in section.items():
                print(''.join((indent,varname.ljust(24), '= ',
                               var_format(value), ',')), file=outfile)
            # Close the section
            print(' /', file=outfile)
            print('', file=outfile)

if __name__ == '__main__':
    nmld = read_namelist('input.nml')
    write_namelist(nmld,'new_input.nml') 
    print(nmld['model_nml'])
