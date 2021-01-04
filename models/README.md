
<span id="TOP" class="anchor"></span>

![DARTlogo](https://github.com/NCAR/DART/blob/Manhattan/docs/images/Dartboard7.png)

## This document is really quite out-of-date. TJH 22 Apr 2020

## Hints for porting a new model to DART:

copy this template directory into a DART/models/xxx
directory for your new model.

if the coordinate system for the model is 1d, you're ok as-is.
if model coordinates are 3d, edit the work/path_names_* files
and change location/oned/* to location/threed_sphere/*

if your model is closer to the simpler examples (e.g. lorenz),
the existing model_mod.f90 is a good place to start.
if your model is a full 3d geophysical one (e.g. like cam, pop, etc)
then rename full_model_mod.f90 to model_mod.f90 and start there.

edit all the work/path_names_* files and change models/template/xxx
to use the name of the directory for your model.

try ./quickbuild.csh and everything should compile at this point.

the required subroutines are these:
```
public :: get_model_size,         &
          adv_1step,              &
          get_state_meta_data,    &
          model_interpolate,      &
          get_model_time_step,    &
          end_model,              &
          static_init_model,      &
          init_time,              &
          init_conditions,        &
          nc_write_model_atts,    &
          nc_write_model_vars,    &
          pert_model_state,       &          
          get_close_maxdist_init, &
          get_close_obs_init,     &
          get_close_obs,          &
          ens_mean_for_model
```

in addition, model_mod can contain subroutines that are used
for other utility programs and we recommend at least the following
routines be added to `model_mod.f90`:
```
public :: model_file_to_dart_vector, &     ! converter
          dart_vector_to_model_file, &     ! converter
          get_gridsize,              &     ! called by everyone
          get_model_filename,        &     ! called by both (set_model_filename?)
          get_state_time,            &     ! model_to_sv, static_init_model
          set_state_time  !(?)             ! sv_to_model, trans_time
```

edit the model mod and fill in the routines in this order:

1. `static_init_model()` - make it read in the grid information
  and the number of variables that will be in the state vector
 (fill in the progvar derived type).   fill in the model_size
  variable.  as part of this work, fill in the `get_gridsize()`
  code.

  after number 1 is done, `get_model_size()` and 
  `get_model_time_step()` from the template should be ok as-is.

2. model_file_to_dart_vector() - given a model data file, read in
  the fields and put them into the 1D DART state vector.  make
  sure the order and values match the progvar array.  

3. dart_vector_to_model_file() - do the reverse of the previous step.

4. `get_state_meta_data()` - given an index number into the state vector
    return the location and type.  the code which loops over the
    progvar should already do the type, but code to compute what
    lon, lat, and vertical (for a 3d model) or x location (1d)
    corresponds to this item must be written.

5. `model_interpolate()` - given a location (lon/lat/vert in 3d, x in 1d)
   and a state QTY_xxx kind, return the interpolated value the field
   has at that location.   this is probably one of the routines that
   will take the most code to write.

6. `nc_write_model_atts(), nc_write_model_vars()` - when `filter` runs
   it calls these routines to output model data into a netcdf diagnostic
   file which is unrelated to the model data files.  it is possible to
   have the ensemble data just be dumped as a single 1D vector but
   that makes the files less useful.  generally it's most useful to
   put the grid info and dump each variable from the state vector
   into a separate netcdf variable.  the diagnostic files contain the
   ensemble mean, the ensemble stddev, the inflation information, and
   then optionally some or all of the individual ensemble members.

for now, ignore these routines:
```
   get_close_maxdist_init()
   get_close_obs_init()
   get_close_obs()
   ens_mean_for_model()
   end_model()
```
 
if you have data in a dart initial condition/restart file, then you
can ignore these routines:
```
   init_time()
   init_conditions()
```
otherwise, have them return an initial time and an initial default
ensemble state.

if your model is NOT subroutine callable, you can ignore this routine:
```
   adv_1step()
```
otherwise have it call the interface to your model and add the files
necessary to build your model to all the `work/path_names_*` files.
add the model source to a src/ directory.

if you want to let filter add gaussian noise to a single state vector
to generate an ensemble, you can ignore this routine
```
   pert_model_state()
```
otherwise fill in code that does whatever perturbation makes sense
to have an initial ensemble of states.  in some cases that means
adding a different range of values to each different field in the
state vector.

at this point you should have enough code to test and run simple
experiments.  the `model_mod_check` utility program can be used 
during this process to check your implementation.


the general flow is:

   1) `./model_to_dart` - read model data and convert it into a dart state vector file
   2) `./create_obs_sequence` - make a file with a single observation in it
   3) `./perfect_model_obs` - should interpolate a value for the obs
   4) `./dart_to_model` - convert the dart vector back into a model data file

   5) generate an ensemble of states, or set 'start_from_restart' to .false.
   6) run `./filter` with the single observation 
   7) look at the preassim.nc and analysis.nc files
        diff them with ncdiff:  ncdiff analysis.nc preassim.nc Innov.nc
        plot it, with ncview if possible:  ncview Innov.nc
        the difference between the two is the impact of that single observation
        see if it's at the right location and if the differences seem reasonable

more details on each of these 7 steps follows.

### 1) model_to_dart
this program needs to read the output file from the model,
whatever format that is (many of our supported models use netcdf).
it needs to create a 1d array of values in whatever order it chooses.
the model_mod code must be able to take any index into that array
(say array(25)) and be able to return what location and variable kind
that corresponds so, so the mapping from 2d or 3d array to this 1d array
has to be kept track of in the model_mod code so it can be inverted
on demand.

### 2) create_obs_sequence
you can make a synthetic observation (or a series of them) with this
interactive program and use them for testing.  before running, make sure
the observation types you want to use are in the input.nml file in the
&obs_kind_nml section, either in the assimilate or evaluate lists.

then run the program.  give the total number of obs you want to create
(start with 1).  answer 0 to number of data items and 0 to number of
quality control items.  answer 0 when it says enter -1 to quit.  you
will be prompted for an observation number to select what type of
observation you are going to test.  

give it a location that should be inside your domain, someplace where
you can compute (by hand) what the correct value should be.  when it
asks for time, give it a time that is the same as the time on your
model data.

when it asks for error variance, at this point it doesn't matter.
give it something like 10% of the expected data value.  later on
this is going to matter a lot, but for testing the interpolation of
a single synthetic obs, this will do.

for an output filename, it suggests 'seq_def.out' but in this case,
tell it 'obs_seq.in'.

### 3) perfect_model_obs
if you have run the model_to_dart and created a state vector, make sure
the name matches the input name in the input.nml file, the &perfect_model_obs_nml
namelist.  make sure the input obs_sequence is still set to 'obs_seq.in'.
run perfect_model_obs.  something bad will happen, most likely.  fix it.
eventually it will run and you will get an 'obs_seq.out' file.  for these
tests, make sure &obs_sequence_nml : write_binary_obs_sequence = .false.
in the input.nml file.  the sequence files will be short and in ascii.
you can check to see what the interpolated value is.  if it's right, yay.
if not, debug the interpolation code in the model_mod.f90 file.

### 4) dart_to_model
if you have run perfect_model_obs, you have not changed the dart
state vector in any way.  however, it's a good test to make a copy of
the model input file, then run 'model_to_dart' and then 'dart_to_model'
(you may need to set 'is_model_advance_file' to .false. for the
&dart_to_model_nml namelist) and you should get identical values back
in the model input file as you started with.

### 5) running from a single input state
in the &filter_nml namelist, set 'start_from_restart' to .false.
this tells filter that you have not generated N initial conditions,
that you are only going to supply one and it needs to perturb that
one to generate an initial ensemble.  make sure the 'restart_in_file_name' 
matches the name of the single state vector file you have.  you can
use the 'obs_seq.out' file from the perfect_model run because now
it has data for that observation.  later on you will need to decide
on how to generate a real set of initial states, and then you will
set 'start_from_restart' back to .true. and supply N files instead of one.
you may need to set the &ensemble_manager_nml : perturbation_amplitude
down to something smaller than 0.2 for these tests - 0.00001 is a good
first guess for adding small perturbations to a state.

### 6) filter
set the ens_size to something small for testing - between 4 and 10 is
usually a good range.  make sure your observation type is in the
'assimilate_these_obs_types' list and not in the evaluate list.
run filter.  find bugs and fix them until the output 'obs_seq.final' 
seems to have reasonable values.  running filter will generate two
netcdf diagnostic files: preassim.nc and analysis.nc

### 7) diagnostics
run 'ncdiff analysis.nc preassim.nc differences.nc' and use
your favorite netcdf plotting tool to see if there are any differences
between the 2 files.  for modules using a regular lat/lon grid 'ncview'
is a quick and dirty way to scan files.  for something on an irregular
grid a more complicated tool will have to be used.  if the files are
identical the assimilation didn't do anything.  check to see if there
is a non-zero DART quality control value in the obs_seq.final file.
check to see if there are errors in the dart_log.out file.  figure out
why there's no change.  if there is a difference, it should be at
the location of the observation and extend out from it for a short
distance.  if it isn't in the right location, look at your get_state_meta_data()
code.  if it doesn't have a reasonable value, look at your model_interpolate() code.

there's lots more to say about this, but this is a quick pointer to
how to get started.

