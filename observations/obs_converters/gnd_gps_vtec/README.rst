GND GPS VTEC
============

This is a modification of a standard "text" converter that comes with DART. 

``gnd_gps_vtec_text_to_obs.f90`` reads VTEC text files 
(from OpenMadrigal at http://madrigal.haystack.mit.edu/)
and outputs DART obs_seq.out files.

Please examine ``work/input.nml:&text_to_obs_nml`` as it specifies the name of the input and the output files

The provided file work/gps021201g.002.txt is only for example 
(only 2 datapoints are shown) and not for real estimation.

