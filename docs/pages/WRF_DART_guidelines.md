--- 
title: WRF and DART
layout: default
---

### The WRF Weather Model and DART

<span id="WRFDART"></span>

> How does DART interact with running WRF?

Most users with large WRF domains run a single cycle of filter to do
assimilation, and then advance each ensemble member of WRF from a
script, possibly submitting them in a batch to the job queues.  
  
For smaller WRF runs, if WRF can be compiled without MPI (the 'serial'
configuration) then filter can cycle inside the same program, advancing
multiple ensemble members in parallel. See the WRF documentation pages
for more
details.

> I have completed running *filter* and I have the *filter_restart.\#\#\#\#* files.
> Can you refer me to the utility to convert them back to a set of *wrfinput_d01* files?

<!-- TJH FIXME ... this section no longer appropriate for Manhattan ... -->

If you are using the *advance_model.csh* script that is distributed
with DART, it will take care of converting the filter output files back
to the WRF input files for the next model advance.  
  
If you are setting up a free run or doing something different than what
the basic script supports, read on to see what must be done.  
  
When you finish running DART it will have created a set of
*sssss.\#\#\#\#* restart files, where the *sssss* part of the filename
comes from the setting of *&filter_nml :: restart_out_file_name*
(and is frequently *filter_restart*). The *.\#\#\#\#* is a 4 digit
number appended by filter based on the ensemble number. These files
contain the WRF state vector data that was used in the assimilation,
which is usually a subset of all the fields in a *wrfinput_d01* file.  
  
*dart_to_wrf* is the standard utility to insert the DART state
information into a WRF input file, e.g. *wrfinput_d01*. For multiple
WRF domains, a single run of the converter program will update the
*_d02*, *_d03*, ..., files at the same time as the *_d01* file.  
  
In the *input.nml* file, set the following:

~~~
&dart_to_wrf_nml
   model_advance_file = .false.
   dart_restart_name  = 'filter_restart.####',
/
~~~

where '\#\#\#\#' is the ensemble member number. There is no option to
alter the input/output WRF filename. Run *dart_to_wrf*. Remember to
preserve each *wrfinput_d01* file or you will simply keep overwriting
the information in the same output file. Repeat for each ensemble member
and you will be ready to run WRF to make ensemble forecasts.  
  
If filter is advancing the WRF model, and you want to spawn forecasts
from intermediate assimilation steps:  
Use the *assim_model_state_ic.\#\#\#\#* files instead of the
*filter_restart.\#\#\#\#* files, and set the *model_advance_file*
namelist item to be *.true.* .  

\[[top](#)\]

-----
