---
title: Documentation
layout: default
---

## Adding your efforts to DART.

Please let us know if you have suggestions or code to contribute to
DART. We're a small group, but we are willing to listen and will make
every effort to incorporate improvements to the code. Email us at
<dart@ucar.edu>.

# DART documentation

The DART distribution includes a full set of documentation. Once you
download DART, you may view the documentation offline by opening the
*index.html* file in the top-level DART directory. If you want to
explore the documentation page without downloading DART, you may 
[view the documentation for the Manhattan release](https://ncar.github.io/DART/api/v2.1.10/index.html).

## Links to major sections of this document:
- [Downloadable datasets for DART.](#datasets)
- [Creating initial conditions for DART](#creating_ics)
- ['Perfect Model' or 'OSSE' experiments](#perfect_osse)
- [Adding a model to DART](#adding_a_model)

----- 

## DART-supported models:

There are two broad classes of models supported by DART. Some are
'low-order' models, generally single-threaded, subroutine-callable, and
idealized: there are no **real** observations of these systems. The
other class of models are 'high-order' models. There **are** real
observations of these systems. Or at least, we like to think so ...

### Models that are ready to use with Manhattan:
   [lorenz_63](#lorenz_63)
   [lorenz_84](#lorenz_84)
   [9var](#NINEvar)
   [lorenz_96](#lorenz_96)
   [lorenz_96_2scale](#lorenz_96_2scale)
   [forced_lorenz_96](#forced_lorenz_96)
   [lorenz_04](#lorenz_04)
   [simple_advection](#simple_advection)
   [bgrid_solo](#bgrid_solo)
   [WRF](#wrf)
   [MPAS ATM](#mpas_atm)
   [ROMS](#ROMS)
   [CESM](#CESM)
   [CAM-FV](#cam-fv)
   [CAM-CHEM](#cam-chem)
   [WACCM](#WACCM) 
   [WACCM-X](#WACCM-X)
   [CICE](#CICE)
   [POP](#POP)
   [CM1](#CM1)
   [FESOM](#fesom)
   [NOAH-MP](#noah-mp)
   [WRF-Hydro](#wrf-hydro) 

### Models supported in Lanai:
   [GCCOM](#GCCOM)
   [LMDZ](#LMDZ) 
   [MITgcm_ocean](#MITgcm_ocean)
   [NAAPS](#NAAPS)
   [AM2](#AM2)
   [CAM-SE](#cam-se)
   [CLM](#CLM) 
   [COAMPS](#COAMPS)
   [COSMO](#COSMO)
   [dynamo](#dynamo) 
   [gitm](#gitm)
   [ikeda](#ikeda)
   [jules](#jules) 
   [mpas_ocean](#mpas_ocean)
   [null](#null) 
   [openggcm](#openggcm)
   [parflow](#parflow)
   [sqg](#sqg) 
   [tiegcm](#tiegcm)
   [wrf-chem](#wrf-chem)

### Models that were used a long time ago (these may not take that much work to revive):
   [ECHAM](#ECHAM)
   [PBL_1d](#PBL_1d) 
   [MITgcm_annulus](#MITgcm_annulus)
   [forced_barot](#forced_barot) 
   [pe2lyr](#pe2lyr)
   [ROSE](#rose) 
   [CABLE](#cable) 

<span id="models" class="anchor"></span>  

\[[top](#)\]

----- 

## The 'Manhattan-ready' models in DART.

<span id="lorenz_63" class="anchor"></span>

### lorenz_63

This is the 3-variable model as described in: Lorenz, E. N.
1963. Deterministic nonperiodic flow. *J. Atmos.
Sci.* **20**, 130-141.  
The system of equations is:

~~~ 
X' = -sigma*X + sigma*Y
Y' = -XZ + rX - Y
Z' =  XY -bZ
~~~

<span id="lorenz_84" class="anchor"></span>

### lorenz_84

This model is based on:   Lorenz E. N., 1984: Irregularity: A
fundamental property of the atmosphere. *Tellus*,  **36A**, 98-110.  
The system of equations is:

~~~ 
X' = -Y^2 - Z^2  - aX  + aF
Y' =  XY  - bXZ  - Y   + G
Z' = bXY  +  XZ  - Z
~~~

Where a, b, F, and G are the model parameters.

<span id="NINEvar" class="anchor"></span>

### 9var

This model provides interesting off-attractor transients that behave
something like gravity waves.

<span id="lorenz_96" class="anchor"></span>

### lorenz_96

This is the model we use to become familiar with new architectures,
i.e., it is the one we use 'first'. It can be called as a subroutine or
as a separate executable. We can test this model both single-threaded
and mpi-enabled.  
  
Quoting from the **Lorenz 1998** paper:

> ... the authors introduce a model consisting of 40 ordinary
> differential equations, with the dependent variables representing
> values of some atmospheric quantity at 40 sites spaced equally about a
> latitude circle. The equations contain quadratic, linear, and constant
> terms representing advection, dissipation, and external forcing.
> Numerical integration indicates that small errors (differences between
> solutions) tend to double in about 2 days. Localized errors tend to
> spread eastward as they grow, encircling the globe after about 14
> days.  
> ...  
> We have chosen a model with J variables, denoted by X<sub>1</sub>,
> ..., X<sub>J</sub>; in most of our experiments we have let J = 40. The
> governing equations are:  
>   
> 
> ~~~ 
> dXj/dt = (Xj+1 - Xj-2)Xj-1 - Xj + F         (1)
> ~~~
> 
>   
> for *j* = 1, ..., J. To make Eq. (1) meaningful for all values of *j*
> we define X<sub>-1</sub> = X<sub>J-1</sub>, X<sub>0</sub> =
> X<sub>J</sub>, and X<sub>J+1</sub> = X<sub>1</sub>, so that the
> variables form a cyclic chain, and may be looked at as values of some
> unspecified scalar meteorological quantity, perhaps vorticity or
> temperature, at J equally spaced sites extending around a latitude
> circle. Nothing will simulate the atmosphere's latitudinal or vertical
> extent.

<span id="lorenz_96_2scale" class="anchor"></span>

### lorenz_96_2scale

This is the Lorenz 96 2-scale model, documented in Lorenz (1995). It
also has the option of the variant on the model from Smith (2001), which
is invoked by setting *local_y = .true.* in the namelist. The time
step, coupling, forcing, number of X variables, and the number of Ys per
X are all specified in the namelist. Defaults are chosen depending on
whether the Lorenz or Smith option is specified in the namelist. Lorenz
is the default model. Interface written by **Josh Hacker**. Thanks
Josh\!

<span id="forced_lorenz_96" class="anchor"></span>

### forced_lorenz_96

The forced_lorenz_96 model implements the standard L96 equations
except that the forcing term, F, is added to the state vector and is
assigned an independent value at each gridpoint. The result is a model
that is twice as big as the standard L96 model. The forcing can be
allowed to vary in time or can be held fixed so that the model looks
like the standard L96 but with a state vector that includes the constant
forcing term. An option is also included to add random noise to the
forcing terms as part of the time tendency computation which can help in
assimilation performance. If the random noise option is turned off (see
namelist) the time tendency of the forcing terms is 0.
  
<span id="lorenz_04" class="anchor"></span>

### lorenz_04

The reference for these models is Lorenz, E.N., 2005: Designing chaotic
models. *J. Atmos. Sci.*, **62**, 1574-1587.  
Model II is a single-scale model, similar to Lorenz 96, but with spatial
continuity in the waves. Model III is a two-scale model. It is
fudamentally different from the Lorenz 96 two-scale model because of the
spatial continuity and the fact that both scales are projected onto a
single variable of integration. The scale separation is achived by a
spatial filter and is therefore not perfect (i.e. there is leakage). The
slow scale in model III is model II, and thus model II is a deficient
form of model III. The basic equations are documented in Lorenz (2005)
and also in the model_mod.f90 code. The user is free to choose model II
or III with a Namelist variable.

<span id="simple_advection" class="anchor"></span>

### simple_advection

This model is on a periodic one-dimensional domain. A wind field is
modeled using Burger's Equation with an upstream semi-lagrangian
differencing. This diffusive numerical scheme is stable and forcing is
provided by adding in random gaussian noise to each wind grid variable
independently at each timestep. An Eulerian option with
centered-in-space differencing is also provided. The Eulerian
differencing is both numerically unstable and subject to shock
formation. However, it can sometimes be made stable in assimilation mode
(see recent work by Majda and collaborators).

<span id="bgrid_solo" class="anchor"></span>

### bgrid_solo

This is a dynamical core for B-grid dynamics using the Held-Suarez
forcing. The resolution is configurable, and the entire model can be run
as a subroutine. *Status:* supported.

<span id="wrf" class="anchor"></span>

### WRF

The [Weather Research and Forecasting (WRF)
Model](http://www.wrf-model.org/) is a next-generation mesoscale
numerical weather prediction system designed to serve both operational
forecasting and atmospheric research needs. More people are using DART
with WRF than any other model. Note: The actual WRF code is not
distributed with DART. *Status:* supported.

<span id="mpas_atm" class="anchor"></span>

### MPAS ATM

[Model Prediction Across Scales -
atmosphere](https://mpas-dev.github.io/) *Status:* active
 
<span id="ROMS" class="anchor"></span>

### ROMS

[Regional Ocean Modelling System](https://www.myroms.org/) *Status:* active

<span id="CESM" class="anchor"></span>

### CESM

There are several 
[supported versions](http://www.cesm.ucar.edu/models/current.html) 
of the Community Earth System Model (CESM) and its ancestors
([CCSM4.0](http://www.cesm.ucar.edu/models/ccsm4.0)). Contact us for
support for unreleased, developmental versions of CESM. Not all are
supported because each requires modification of some subroutines and
setup scripts in order to work with DART. The supported versions depend
to some degree on the CESM component(s) which will be used as the
assimilating model(s). See CAM-FV, POP, and CICE, below, and
[CESM setup guidelines](CESM_setup_guidelines.md). In
general, later versions of CESM can build the latest component models
plus any earlier versions of the component models. For example, CESM2.0
can build CAM-FV version 6, 5, 4, ... while CESM1.2.1 can build CAM-FV
5, 4, ..., but not 6. Note: the source code for CESM component models is
not distributed with DART.

<span id="cam-fv" class="anchor"></span>
<span id="WACCM" class="anchor"></span>
<span id="WACCM-X" class="anchor"></span>
<span id="cam-chem" class="anchor"></span>

### CAM-FV CAM-Chem WACCM WACCM-X

See CESM, above. CAM-FV has been the "work horse" atmospheric climate
model for several generations of CESM releases. The Manhattan release of
DART provides interfaces to CESM1.5 and CESM2.0. The setup scripts for
those CESMs will currently build only the finite volume dynamical core
of CAM. This works for all of the variants of CAM-FV; CAM-Chem, WACCM,
WACCM-X. An interface between DART and the spectral element dy-core of
CAM is available in DART Classic and will be brought into the Manhattan
release when needed.
[(CAM5)](http://www.cesm.ucar.edu/models/cesm1.0/cam); Some SourceMods
and initial file ensembles for older and lower-resolution CAM-FVs are
available in 
[DART/CAM datasets](http://www.image.ucar.edu/pub/DART/CAM/) *Status:* available
for community use.

<span id="CICE" class="anchor"></span>

### CICE (pronounced 'sea ice')

See CESM, above. The sea-ice component of
[CESM](http://www.cesm.ucar.edu/models/current.html) The interface of
[CESM-CICE](http://www.cesm.ucar.edu/models/cesm1.2/cice/) to DART is
through CESM1.5. **Cecilia Bitz** and **Yongfei Zhang** created the
interfaces for DART.  *Status:* throroughly
beta-tested, full support awaiting the CESM2.0 release.

<span id="POP" class="anchor"></span>

### POP

See CESM, above. The Parallel Ocean Program
[(POP)](http://www.cesm.ucar.edu/models/cesm1.0/pop2/) was originally created by
the Los Alamos National Laboratory and has been modified to run in the [NCAR
Community Earth System Model
(CESM)](http://www.cesm.ucar.edu/models/current.html) framework. Additional
modifications are necessary for data assimilation and center around the need to
perform an adjustment upon restart to account for the fact that the input ocean
state has been modified by the assimilation. There are interfaces for CESM1.1.1
and CESM1.2.1. *Status:* available for community use.

<span id="CM1" class="anchor"></span>

### CM1

Cloud Model 1 (CM1) version 18 (CM1r18) is a non-hydrostatic numerical model in
Cartesian 3D coordinates designed for the study of micro- to mesoscale
atmospheric phenomena in idealized to semi-idealized simulations. The CM1 model
was developed and is maintained by George Bryan at the National Center for
Atmospheric Research (NCAR) Mesoscale and Microscale Meteorology Laboratory
(MMM). The model code is freely available from the CM1 website:
<http://www2.mmm.ucar.edu/people/bryan/cm1> and must be downloaded and compiled
outside of DART.
This model interface and scripting support were created by **Luke Madaus**.

<span id="fesom" class="anchor"></span>

### FESOM

<span id="noah-mp" class="anchor"></span>

[FESOM](https://fesom.de/models/fesom14) is an unstructured mesh global 
ocean model using finite element methods to solve the hydrostatic 
primitive equations with the Boussinesq approximation.
The [FESOM model interface](../../models/FESOM/Readme.md),
scripting support and some diagnostic routines were 
contributed by **Ali Aydoğdu**.
*Status:* available for community use.

### NOAH-MP

<span id="wrf-hydro" class="anchor"></span>

### WRF-HYDRO

The WRF-Hydro assimilation support has its own (private) GitHub repository
[NCAR/wrf_hydro_dart](https://github.com/NCAR/wrf_hydro_dart) that supports
the channel-only configuration of WRF-Hydro. Originally, this was almost 
entirely the work of **James McCreight** of NCAR's Research Applications 
Laboratory (RAL). The DAReS team has been working with RAL to incorporate new
features such as localization restricted to watersheds, new inflation algorithms
and variable transformations that provide much better results when assimilating
non-gaussian quantities such as streamflow. The private wrf_hydro_dart repository
is expected to be released in a public version very soon.


\[[top](#)\]

-----

<span id="models_in_progress" class="anchor"></span>  

## Models that are supported by DART Classic and could be ported to the Manhattan release if needed.

<span id="GCCOM" class="anchor"></span>

### General Curvilinear Coastal Ocean Model - GCCOM

GCCOM is a three-dimensional, nonhydrostatic Large Eddy Simulation (LES), 
rigid lid model that has the ability to run in a fully three-dimensional 
general curvilinear coordinate system. Much of the work of supporting GCCOM in
DART was by **Mariangel Garcia** while she was at San Diego State University.
One article is 
["Interfacing an ensemble Data Assimilation system with a 3D nonhydrostatic Coastal Ocean Model, an OSSE experiment"](https://ieeexplore.ieee.org/abstract/document/7760992)

<span id="LMDZ" class="anchor"></span>

### LMDZ

The DART interfaces were prototyped by **Tarkeshwar Singh** of the
Centre for Atmospheric Sciences, Indian Institute of Technology (IIT) Delhi.
From the LMDZ homepage:
> LMDZ is a general circulation model (or global climate model) developed 
> since the 70s at the "Laboratoire de Météorologie Dynamique", which 
> includes various variants for the Earth and other planets (Mars, 
> Titan, Venus, Exoplanets). The 'Z' in LMDZ stands for "zoom" 
> (and the 'LMD' is for  'Laboratoire de Météorologie Dynamique").


<span id="MITgcm_ocean" class="anchor"></span>

### MITgcm_ocean

The [MIT ocean GCM](http://mitgcm.org/) version 'checkpoint59a' is the
foundation of this implementation. It was modified by **Ibrahim Hoteit**
(then of Scripps) to accomodate the interfaces needed by DART. *Status:*
supported, and currently being ported to Manhattan.
 
<span id="NAAPS" class="anchor"></span>

### NAAPS
 
<span id="AM2" class="anchor"></span>

### AM2

The [FMS AM2](http://data1.gfdl.noaa.gov/~arl/pubrel/m/am2/doc/) model
is GFDL's atmosphere-only code using observed sea surface temperatures,
time-varying radiative forcings (including volcanos) and time-varying
land cover type. This version of AM2 (also called AM2.1) uses the
finite-volume dynamical core (Lin 2004). **Robert Pincus** (CIRES/NOAA ESRL
PSD1) and **Patrick Hoffman** (NOAA) wrote the DART interface and are
currently using the system for research. Note: the model code is not
distributed with DART. *Status:* supported

<span id="cable" class="anchor"></span>

### CABLE

The Community Atmosphere Biosphere Land Exchange (CABLE) model is a 
land surface model,used to calculate the fluxes of momentum, energy, 
water and carbon between the land surface and the atmosphere and to 
model the major biogeochemical cycles of the land ecosystem. The DART
interfaces for the standalone version of CABLE have preliminary support
and needs to be updated to be consistent with the Manhattan release.


<span id="cam-se" class="anchor"></span>

### CAM-SE

<span id="CLM" class="anchor"></span>

### CLM

Assimilation with the [Community Land Model](http://www.cesm.ucar.edu/models/clm/)
is well supported and the system has been used for many research interests, from 
biogeochemistry to snow, ice, soil moisture and more. DART/CLM has many research
branches and guidance for which branch is most appropriate is provided upon request. 
There is support for CLM under the Lanai release and several development branches
that are consistent with the Manhattan release. The version distributed with the 
Manhattan release is not as fully functional as the development branches.
Much of the original DART/CLM support was written by **Yongfei Zhang** while 
she was at the University of Texas at Austin.

<span id="COAMPS" class="anchor"></span>

### coamps_nest

The DART interface was originally written and supported by **Tim Whitcomb**.
The following model description is taken from the [COAMPS overview web
page:](http://www.nrlmry.navy.mil/coamps-web/web/view)

> The Coupled Ocean/Atmosphere Mesoscale Prediction System (COAMPS) has
> been developed by the Marine Meteorology Division (MMD) of the Naval
> Research Laboratory (NRL). The atmospheric components of COAMPS,
> described below, are used operationally by the U.S. Navy for
> short-term numerical weather prediction for various regions around the
> world.

Note: the model code is not distributed with DART. *Status:* supported

<span id="COSMO" class="anchor"></span>

### COSMO

<span id="dynamo" class="anchor"></span>

### dynamo

A Flux-Transport Dynamo model from **Mausumi Dikpati**.
The goal of this interface is to estimate the time variation of 
velocities to match given spatio-temporal observation of magnetic fields.

<span id="ikeda" class="anchor"></span>

### ikeda

The Ikeda model is a 2D chaotic map useful for visualization data
assimilation updating directly in state space. There are three
parameters: a, b, and mu. The state is 2D, x = \[X Y\]. The equations
are:

~~~ 
X(i+1) = 1 + mu * ( X(i) * cos( t ) - Y(i) * sin( t ) )
Y(i+1) =     mu * ( X(i) * sin( t ) + Y(i) * cos( t ) ),
~~~

where

~~~ 
t = a - b / ( X(i)**2 + Y(i)**2 + 1 )
~~~

Note the system is time-discrete already, meaning there is no delta_t.
The system stems from nonlinear optics 
(Ikeda 1979, Optics Communications).
Interface written by **Greg Lawson**. Thanks Greg\!

<span id="jules" class="anchor"></span>

### JULES

<span id="mpas_ocean" class="anchor"></span>

### MPAS Ocean

<span id="null" class="anchor"></span>

### null_model

This model provides very simple models for evaluating filtering
algorithms. It can provide simple linear growth around a fixed point, a
random draw from a Gaussian, or combinations of the two.

<span id="openggcm" class="anchor"></span>

### OpenGCCM

<span id="parflow" class="anchor"></span>

### PARFLOW

<span id="sqg" class="anchor"></span>

### SQG

<span id="tiegcm" class="anchor"></span>

### TIEGCM

The DART interfaces to the Thermosphere Ionosphere Electrodynamic General 
Circulation Model [TIEGCM](http://www.hao.ucar.edu/modeling/tgcm/tie.php)
are fully supported in the Lanai release.
TIEGCM is a community model developed at the NCAR High Altitude Observatory and
is widely used by the space physics and aeronomy community.
DART/TIEGCM has been used to assimilate neutral mass density retrieved from 
satellite-borne accelerometers and electon density obtained from ground-based 
and space-based GNSS signals. TIEGCM2 is not yet supported, and the existing
interfaces need to be updated to work under the Manhattan release.


\[[top](#)\]

-----

<span id="orphans" class="anchor"></span>  

## Models that were used a long time ago

<span id="ECHAM" class="anchor"></span>

### ECHAM

*Status:* orphaned.

<span id="PBL_1d" class="anchor"></span>

### PBL_1d

The PBL model is a single column version of the WRF model. The
functionality for this has been folded into the regular WRF model_mod
interface so this version is no longer needed. See the WRF model_mod
namelist documentation for how to use the single-column features.
*Status:* orphaned, obsolete.

<span id="MITgcm_annulus" class="anchor"></span>

### MITgcm_annulus

The MITgcm annulus model as configured for this application within DART
is a non-hydrostatic, rigid lid, C-grid, primitive equation model
utilizing a cylindrical coordinate system. For detailed information
about the MITgcm, see http://mitgcm.org *Status:* orphaned.

<span id="forced_barot" class="anchor"></span>

### FORCED_BAROT

*Status:* orphaned.

<span id="gitm" class="anchor"></span>

### GITM

*Status:* orphaned.

<span id="pe2lyr" class="anchor"></span>

### pe2lyr

This model is a 2-layer, isentropic, primitive equation model on a
sphere. *Status:* orphaned.

<span id="rose" class="anchor"></span>

### rose

The rose model is for the stratosphere-mesosphere and was used by Tomoko
Matsuo (now at CU-Boulder) for research in the assimilation of
observations of the Mesosphere Lower-Thermosphere (MLT). Note: the model
code is not distributed with DART. *Status:* orphaned.

<span id="datasets" class="anchor"></span>
<span id="observation_sets" class="anchor"></span>
<span id="verification" class="anchor"></span>  

\[[top](#)\]

-----

# Downloadable datasets for DART.

The code distribution was getting cluttered with datasets, boundary
conditions, intial conditions, ... large files that were not necessarily
interesting to all people who downloaded the DART code. This is
compounded by the fact subversion makes a local (hidden) copy of the
original repository contents, so the penalty for being large is doubled.
It just made sense to make all the large files available on as
'as-needed' basis.  
  
To keep the size of the DART distribution down we have a separate
www-site to provide some observation sequences, initial conditions, and
general datasets. It is our intent to populate this site with some
'verification' results, i.e. assimilations that were known to be 'good'
and that should be fairly reproducible - appropriate to test the DART
installation.  
  
Please be patient as I make time to populate this directory.
(yes, 'make', all my 'found' time is taken ...)  
  
Observation sequences can be found at 
[www.image.ucar.edu/pub/DART/Obs_sets](http://www.image.ucar.edu/pub/DART/Obs_sets)
  
<span id="initial_conditions" class="anchor"></span>
Useful bits for CAM can be found at 
[www.image.ucar.edu/pub/DART/CAM](http://www.image.ucar.edu/pub/DART/CAM).  
  
Useful bits for WRF can be found at
[www.image.ucar.edu/pub/DART/WRF](http://www.image.ucar.edu/pub/DART/WRF).  
  
Useful bits for MPAS_ocn can be found at
[www.image.ucar.edu/pub/DART/MPAS_OCN](http://www.image.ucar.edu/pub/DART/MPAS_OCN)
  
Useful bits for CICE can be found at
[www.image.ucar.edu/pub/DART/CICE](http://www.image.ucar.edu/pub/DART/CICE)
  
Verification experiments will be posted to
[www.image.ucar.edu/pub/DART/VerificationData](http://www.image.ucar.edu/pub/DART/VerificationData) 
as soon as I can
get to it. These experiments will consist of initial conditions files
for testing different high-order models like CAM, WRF, POP ... The
low-order models are still distributed with verification data in their
*work* directories.

<span id="creating_ics" class="anchor"></span> 

\[[top](#)\]

-----

# Creating initial conditions for DART

The idea is to generate an ensemble that has sufficient 'spread' to
cover the range of possible solutions. Insufficient spread can (and
usually will) lead to poor assimilations. Think 'filter divergence'.  
  
Generating an ensemble of initial conditions can be done in lots of
ways, only a couple of which will be discussed here. The first is to
generate a single initial condition and let DART perturb it with noise
of a nature you specify to generate as many ensemble members as you
like. The second is to take some existing collection of model states and
convert them to DART initial conditions files and then use the 
[NCO operators](http://nco.sourceforge.net/) to set the proper date in the
files. The hard part is then coming up with the original collection of
model state(s).

### Adding noise to a single model state

This method works well for some models, and fails miserably for others.
As it stands, DART supplies a routine that can add gaussian noise to
every element of a state vector. This can cause some models to be
numerically unstable. You can supply your own
*model_mod:pert_model_copies()* if you want a more sophisticated
perturbation scheme.

### Using a collection of model states.

Simply collect the filenames of all the model netCDF files - one per
ensemble member - and specify them through the *input.nml* 
```&filter_nml:input_state_file_list = "restarts_in.txt"```  
  
Frequently, the initial ensemble of restart files is some climatological
collection. For CAM experiments, we usually start with *N* different
'January 1' states ... from *N* different years. The timestamp in those
files can be ignored through namelist control. Experience has shown that
it takes less than a week of assimilating 4x/day to achieve a steady
ensemble spread. WRF has its own method of generating an initial
ensemble. For that, it is best to go to contact someone familiar with
WRF/DART.

<span id="low_order_spinup" class="anchor"></span>

### Initial conditions for the low-order models.

In general, there are 'restart files' for the low-order models that
already exist as ASCII sources for netCDF files. These files are usually
called ```work/filter_input.cdl``` and can be converted to netCDF files by
using the *ncgen -o* unix command.  
  
You can generate your own ensemble by adding noise to a single
```perfect_input.nc``` file and run *filter*. 
The way to specify the input state file is to use the
*input_state_file_list* mechanism. Simply put the name of the file
into the file referenced by *input_state_file_list*. In this example,
```filter_input_list.txt```  would contain exactly one line - the string
"perfect_input.nc". You will also need an observation sequence file,
and you may want to explicitly state the start/stop times.

~~~
&filter_nml
    ...
    input_state_files            = 'null',
    input_state_file_list        = 'filter_input_list.txt'
    perturb_from_single_instance = .true.
    ens_size                     = [whatever you want]
    init_time_days               = 0
    init_time_seconds            = 0
    first_obs_days               = 0
    first_obs_seconds            = 0
    last_obs_days                = 0
    last_obs_seconds             = 0
    output_state_files           = 'null',
    output_state_file_list       = 'filter_output_list.txt'
    ...
~~~

In this example, the ensemble will be created with whatever file name
you put in ```filter_output_list.txt```.

<span id="perfect_osse" class="anchor"></span>  
<span id="osse_simple" class="anchor"></span>

\[[top](#)\]

-----

# 'Perfect Model' or 'OSSE' experiments.

Once a model is compatible with the DART facility all of the
functionality of DART is available. This includes 'perfect model'
experiments (also called Observing System Simulation Experiments -
OSSEs). Essentially, the model is run forward from a known state and, at
predefined times, an observation forward operator is applied to the
model state to harvest synthetic observations. This model trajectory is
known as the 'true state'. The synthetic observations are then used in
an assimilation experiment. The assimilation performance can then be
evaluated precisely because the true state (of the model) is known.
Since the same forward operator is used to harvest the synthetic
observations as well as during the assimilation, the
'representativeness' error of the assimilation system is not an issue.  
  
The example described in this section uses low-order models, but the
logic and procedure is **exactly** the same for high-order models; the
complication is usually that researchers want more sophisticated
observation networks than those described here. All you have to do is
use the DART tools to create an observation sequence file (even a REAL
observation sequence file), and use that instead of creating one by hand
with *create_obs_sequence* and *create_fixed_network_sequence*.
*perfect_model_obs* will simply ignore the actual observation values
in this case and only use the observation metadata. Take care that the
observation error values in the file are appropriate for an OSSE - the
converters usually assume some sort of representativeness error in the
observation error specification.  
  
There are a set of MATLAB® functions to help explore the assimilation
performance in state-space as well as in observation-space.

<!-- TJH FIXME [Exploring the results of a Lorenz '96 OSSE](Research/Lorenz96/index.html). -->

## Perfect Model Experiment Overview

There are four fundamental steps to running an OSSE from within DART:

1.  [Create a blueprint](#obs_blueprint) of what, where, and when you
    want observations. Essentially, define the metadata of the
    observations without actually specifying the observation values. The
    default filename for the blueprint is *obs_seq.in*. For simple
    cases, this is just running
    [create_obs_sequence](https://ncar.github.io/DART/api/v2.1.10/program/create_obs_sequence.html)
    and
    [create_fixed_network_seq](https://ncar.github.io/DART/api/v2.1.10/program/create_fixed_network_seq.html).
    You can also use real observation sequences as long as you take care
    to specify observation error variances that do not incorporate
    representativeness error.  
2.  [Harvest the synthetic observations](#run_pmo) from the true model
    state by running
    [perfect_model_obs](https://ncar.github.io/DART/api/v2.1.10/program/perfect_model_obs.html)
    to advance the model from a known initial condition and apply the
    forward observation operator based on the observation 'blueprint'.
    The observation will have noise added to it based on a draw from a
    random normal distribution with the variance specified in the
    observation blueprint. The noise-free 'truth' and the noisy
    'observation' are recorded in the output observation sequence file.
    The entire time-history of the true state of the model is recorded
    in ```perfect_output.nc```. The default filename for the
    'observations' is ```obs_seq.out```.  
3.  [Assimilate the synthetic observations](#run_filter) with
    [filter](https://ncar.github.io/DART/api/v2.1.10/program/filter.html) in
    the usual way. The prior/forecast states are preserved in
    ```preassim.nc``` and the posterior/analysis states are preserved in
    ```filter_output.nc``` . The default filename for the file with the
    observations and (optionally) the ensemble estimates of the
    observations is ```obs_seq.final``` .  
4.  [Check to make sure the assimilation was effective\!](#run_diagnostics)
    Ensemble DA is not a black box\!
    YOU must check to make sure you are making effective use of the
    information in the observations\!

<span id="obs_blueprint" class="anchor"></span>

## 1. Defining the observation metadata - the 'blueprint'.

There are lots of ways to define an observation sequence that DART can
use as input for a perfect model experiment. If you have observations in
DART format already, you can simply use them. If you have observations
in one of the formats already supported by the DART converters
(check [DART/observations/obs_converters/observations.html](obs_converters_observations.html)),
convert it to a DART observation sequence. You may need to use the
[obs_sequence_tool](https://ncar.github.io/DART/api/v2.1.10/program/obs_sequence_tool.html)
to combine multiple observation sequence files into observation sequence
files for the perfect model experiment. Any existing observation values
and quality control information will be ignored by *perfect_model_obs*;
only the time and location information are used. In fact, any and all
existing observation and QC values will be removed.  
  
GENERAL COMMENT ABOUT THE INTERPLAY BETWEEN THE MODEL STOP/START
FREQUENCY AND THE IMPACT ON THE OBSERVATION FREQUENCY: There is usually
a very real difference between the dynamical timestep of the model and
when it is safe to stop and restart the model. The assimilation window
is (usually) required to be a multiple of the safe stop/start frequency.
For example, an atmospheric model may have a dynamical timestep of a few
seconds, but may be constrained such that it is only possible to
stop/restart every hour. In this case, the assimilation window is a
multiple of 3600 seconds. Trying to get observations at a finer
timescale is not possible, we only have access to the model state when
the model stops.

If you do not have an input observation sequence, it is simple to create one.

1.  Run
    [create_obs_sequence](https://ncar.github.io/DART/api/v2.1.10/program/create_obs_sequence.html)
    to generate the blueprint for the types of observations and
    observation error variances for whatever locations are desired.  
2.  Run
    [create_fixed_network_seq](https://ncar.github.io/DART/api/v2.1.10/program/create_fixed_network_seq.html)
    to define the temporal distribution of the desired observations.

Both *create_obs_sequence* and *create_fixed_network_seq*
interactively prompt you for the information they require. This can be
quite tedious if you want a spatially dense set of observations. People
have been known to actually write programs to generate the input to
*create_obs_sequence* and simply pipe or redirect the information into
the program. There are several examples of these in the
*models/bgrid_solo* directory: ```column_rand.f90,
id_set_def_stdin.f90, ps_id_stdin.f90,``` and ```ps_rand_local.f90``` .
Be advised that some observation types have different input
requirements, so a 'one size fits all' program is a waste of time.

NOTE: only the observation kinds in the ```input.nml 
&obs_kind_nml:assimilate_these_obs_types,evaluate_these_obs``` 
will be available to the *create_obs_sequence* program.

DEVELOPERS TIP: You can specify 'identity' observations as input to
*perfect_model_obs*. Identity observations are the model values AT the
exact gridcell location, there is no interpolation at all. Just a
straight table-lookup. This can be useful as you develop your model
interfaces; you can test many of the routines and scripts without having
a working *model_interpolate()*.

More information about creating observation sequence files for OSSE's is
available in the 
[Synthetic Observations section](Observations.md#obs_synthetic).

<span id="run_pmo" class="anchor"></span>

## 2. Generating the true state and harvesting the observation values - *perfect_model_obs*

[perfect_model_obs](https://ncar.github.io/DART/api/v2.1.10/program/perfect_model_obs.html)
reads the blueprint and an initial state and applies the appropriate
forward observation operator for each and every observation in the
current 'assimilation window'. If necessary, the model is advanced until
the next set of observations is desired. When it has run out of
observations or reached the stop time defined by the namelist control,
the program stops and writes out restarts, diagnostics, observation
sequences, and a log file. This is fundamentally a single deterministic
forecast for 'as long as it takes' to harvest all the observations.

<table>
<thead>
<tr class="header">
<th>default filename</th>
<th>format</th>
<th>contents</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><em>perfect_input.nc</em></td>
<td>netCDF</td>
<td>The DART model state to start from. If the variables have a <code>time</code> dimension, The last timestep will be used as the starting point.</td>
</tr>
<tr class="even">
<td><em>perfect_output.nc</em></td>
<td>netCDF</td>
<td>The DART model state at every assimilation timestep. This file has but one 'copy' - the truth. Dump the metadata and the time:<br />
<em>ncdump -v time,MemberMetadata perfect_output.nc</em></td>
</tr>
<tr class="odd">
<td><em>obs_seq.out</em></td>
<td>ASCII or binary  <br />
DART-specific linked list</td>
<td>This file has the observations - the result of the forward 
observation operator. This observation sequence file has two 'copies' 
of the observation: the noisy 'copy' and the noise-free 'copy'. 
The noisy copy is designated as the 'observation', the noise-free 
copy is the truth. The observation-space diagnostic program 
<em>obs_diag</em> has special options for using the true copy instead 
of the observation copy. See the 
<a href="https://ncar.github.io/DART/api/v2.1.10/program/obs_diag.html">obs_diag.html</a> for details.</td>
</tr>
<tr class="even">
<td><em>dart_log.out</em></td>
<td>ASCII</td>
<td>The run-time output of <em>perfect_model_obs</em> .</td>
</tr>
</tbody>
</table>

Each model may define the assimilation window differently, but
conceptually, all the observations plus or minus half the assimilation
window are considered to be simultaneous and a single model state
provides the basis for all those observations. For example: if the
blueprint requires temperature observations every 30 seconds, the
initial model time is noon (12:00) and the assimilation window is 1
hour; all the observations from 11:30 to 12:30 will use the same state
as input for the forward observation operator. The fact that you have a
blueprint for observations every 30 seconds means a lot of those
observations may have the same value (if they are in the same location).

*perfect_model_obs* uses the *input.nml* for its control. A subset of
the namelists and variables of particular interest for
*perfect_model_obs* are summarized here. Each namelist is fully
described by the corresponding module document.

~~~
&perfect_model_obs_nml
    ...
    read_input_state_from_file  = .true.              # some models can start from preset ICs
    single_file_in              = .true               # some models have nested domains ... 
    input_state_files           = 'perfect_input.nc'  # list of files ... for each domain
    write_output_state_to_file  = .true.
    single_file_out             = .true.
    output_state_files          = 'perfect_output.nc' # the time-evolution of the true state
    async                       = 0                   # totally depends on the model
    adv_ens_command             = './advance_ens.csh' #         depends on the model
    obs_seq_in_file_name        = 'obs_seq.in'
    obs_seq_out_file_name       = 'obs_seq.out'
    init_time_days              = -1                  # negative means use the time in ...
    init_time_seconds           = -1                  # the 'restart_in_file_name' file
    first_obs_days              = -1                  # negative means start at the first time in ...
    first_obs_seconds           = -1                  # the 'obs_seq_in_file_name' file.
    last_obs_days               = -1                  # negative means to stop with the last ...
    last_obs_seconds            = -1                  # observation in the file.
   /
    
&obs_sequence_nml
    write_binary_obs_sequence = .false.       #.false. will create ASCII - easy to check.
   /
    
&obs_kind_nml
    ...
    assimilate_these_obs_types = 'RADIOSONDE_TEMPERATURE',
    ...                                       # list all the synthetic observation
    ...                                       # types you want
   /
    
&model_nml
    ...
    time_step_days = 0,                       # some models call this 'assimilation_period_days'
    time_step_seconds = 3600                  # some models call this 'assimilation_period_seconds'
                                              # use what is appropriate for the model
   /
    
&utilities_nml
    ...
    termlevel   = 1                           # your choice
    logfilename = 'dart_log.out'              # your choice
   /
~~~

Since *perfect_model_obs* generally requires advancing the model, and
the model may use MPI or require special ancillary files or forcing
files or ..., it is not possible to provide a single example that will
cover all possibilities. The subroutine-callable models (i.e. the
low-order models) can run *perfect_model_obs* very simply:

> ./perfect_model_obs

<span id="run_filter" class="anchor"></span>

## 3. Performing the assimilation experiment - *filter*

This step is done with the program
[filter](https://ncar.github.io/DART/api/v2.1.10/program/filter.html), which
also uses ```input.nml``` for input and run-time control. A successful
assimilation will depend on many things: an approprite initial ensemble,
monitoring and perhaps correcting the ensemble spread, localization,
etc. It is simply not possible to design a one-size-fits-all system that
will work for all cases. It is **critically important** to analyze the
results of the assimilation and explore ways of making the assimilation
more effective. 
The [DART tutorial](dart_tutorial.md) and the
[DART_LAB](dart_lab.md) exercises
are an invaluable resource to learn and understand how to determine the
effectiveness of, and improve upon, an assimilation experiment. The
concepts learned with the low-order models are directly applicable to
the most complicated models.  
  
**It is important to remember that if *filter* 'terminates normally', it
does not necessarily mean the assimilation was effective\!**  
  
The Manhattan release of DART allows for a very high degree of
customization when it comes to output. To stay focused on the concepts,
I will restrict the examples to models that have
```single_file_in=.true.```, ```single_file_out=.true.```, and
```stages_to_write='preassim','output'```.  
  
*filter* generally produces at least two state-space output diagnostic
files (```preassim.nc``` and ```filter_output.nc```) which contains values of
the ensemble mean, ensemble spread, perhaps the inflation values, and
(optionally) ensemble members for the duration of the experiment.
*filter* also creates an observation sequence file that contains the
input observation information as well as the prior and posterior
ensemble mean estimates of that observation, the prior and posterior
ensemble spread for that observation, and (optionally), the actual prior
and posterior ensemble estimates of that observation. Rather than
replicate the observation metadata for each of these, the single
metadata is shared for all these 'copies' of the observation. See
[An overview of the observation sequence](Observations.md#obs_seq_overview)
for more detail.
*filter* also produces a run-time log file that can greatly aid in
determining what went wrong if the program terminates abnormally.  
  
A very short description of some of the most important namelist
variables is presented here. Basically, I am only discussing the
settings necessary to get *filter* to run. I can guarantee these
settings WILL NOT generate the BEST assimilation. Again, see the module
documentation for a full description of each namelist.

~~~
&filter_nml  <--- link to the full namelist description!
    async                        = 0
    ens_size                     = 40                 # something ≥ 20, please
    num_output_state_members     = 40                 # of FULL DART model states to put in state-space output files
    num_output_obs_members       = 40                 # of ensemble member estimates of observation to save
    obs_sequence_in_name         = 'obs_seq.out'      # output from perfect_model_obs
    obs_sequence_out_name        = 'obs_seq.final'
    init_time_days               = -1                 # the time in the restart file is correct
    init_time_seconds            = -1
    first_obs_days               = -1                 # same interpretation as with perfect_model_obs
    first_obs_seconds            = -1
    last_obs_days                = -1                 # same interpretation as with perfect_model_obs
    last_obs_seconds             = -1
    
    single_file_in               = .true.
    input_state_file_list        = 'filter_input_list.txt'   file containing the list of input files - 1 per domain
    stages_to_write              = 'preassim', 'output'
    single_file_out              = .true.
    output_state_file_list       = 'filter_output_list.txt'  file containing the list of (desired) output files - 1 per domain
    write_all_stages_at_end      = .false.
    
    inf_flavor               = 0,                       0    0 is 'do not inflate'
    ...
   /
    

&quality_control_nml
   input_qc_threshold       =  3.0,
    outlier_threshold       =  3.0               # Observation rejection criterion!
   /

&assim_tools_nml
    filter_kind             = 1             1 is EAKF, 2 is EnKF ...
    cutoff                  = 0.2           this is your localization - units depend on type of 'location_mod'
   /
    
&obs_kind_nml
    assimilate_these_obs_types = 'RAW_STATE_VARIABLE'    Again, use a list ... appropriate for your model
   /
    
&model_nml
    assimilation_perior_days    = 0                      the assimilation interval is up to you
    assimilation_perior_seconds = 3600
   /
~~~

Once the namelist is set, execute *filter* to integrate the ensemble
forward with the final ensemble state written to the files in
```filter_output_list.txt```. For the low-order models and *bgrid_solo*
(i.e. the models that can be run with *single_file_in = .true.* and
```single_file_out = .true.```) the default filenames will be
```preassim.nc``` and ```filter_output.nc``` and will contain values 
for 40 ensemble members once a day.

~~~
mpirun ./filter        -OR-
    
mpirun.lsf ./filter    -OR-
    
./filter               -OR-
    
however YOU run filter on your system!
~~~

<span id="run_diagnostics" class="anchor"></span>

## 4. ASSESS THE PERFORMANCE!

All the concepts of spread, rmse, rank histograms that were taught in
the DART tutorial and in DART_LAB should be applied now. Try the
techniques described in the 
[Did my experiment work?](Diagnostics.md#DidItWork) section.
The 'big three' state-space diagnostics are repeated here because 
they are so important.
The first two require the ```perfect_output.nc```.

| | |
| --- | --- |
| *plot_bins.m* | plots the rank histograms for a set of state variables. This requires you to have all or most of the ensemble members available in the ```preassim.nc``` or ```filter_output.nc``` files. |
| *plot_total_err.m* | plots the evolution of the error (un-normalized) and ensemble spread of all state variables. |
| *plot_ens_mean_time_series.m* | plots the evolution of a set of state variables - just the ensemble mean (and Truth, if available). *plot_ens_time_series.m* is actually a better choice if you can afford to write all/most of the ensemble members to the ```preassim.nc``` and ```filter_output.nc``` files. |

## DON'T FORGET ABOUT THE OBSERVATION-SPACE DIAGNOSTICS\!

<span id="adding_a_model" class="anchor"></span>  
<span id="addingAmodelOverview" class="anchor"></span>

\[[top](#)\]

-----

# Adding a model to DART - Overview

DART is designed to work with many models *without* modifications to the
DART routines **or** the model source code. DART can 'wrap around' your
model in two ways. One can be used if your model can be called as a
subroutine, the other is for models that are separate executables.
Either way, there are some steps that are common to both paths.  
  
Please be aware that several of the high-order models (CAM and WRF, in
particular) have been used for years and their scripts have incorporated
failsafe procedures and restart capabilities that have proven to be
useful but make the scripts complex - more complex than need be for the
initial attempts. Truly, some of the complexity is no longer required
for available platforms. Then again, we're not running one instance of a
highly complicated computer model, we're running N of them.  
  
*NEW* The DART Manhattan release provides native netCDF read/write
support. Consequently, there is no need for translation routines that we
have traditionally been calling *model_to_dart* or *dart_to_model*.
If, however, your model does not use netCDF for I/O, these programs must
be written. We have a lot of experience writing these converters - you
should not be afraid to ask for advice or for code to start from.  
  
*NEW* Manhattan provides a program to help test the required interfaces:
[assimilation_code/programs/model_mod_check/model_mod_check.f90](https://ncar.github.io/DART/api/v2.1.10/program/model_mod_check.html).
Many models start with this and modify it to suit their needs. Be aware
that some of the model-specific *model_mod_check.f90* programs use
deprecated features. Focus on the ones for Manhattan-compliant
models.

## The basic steps to include your model in DART - each of these topics has its own section farther down.

1.  Copy the *models/template* directory and files to your own DART
    model directory.  
2.  Modify the *model_mod.f90* file to return specifics about your
    model. This module MUST contain all the required interfaces (no
    surprise) but it can also contain more interfaces as is convenient.
    The required interfaces calling syntax (argument list) should not be
    modified in any way.  
3.  *If your model cannot be called as a subroutine:* Modify
    *shell_scripts/advance_model.csh* to collect all the input files
    needed to advance the model into a clean, temporary directory,
    convert the state vector file into input to your model, run your
    model, and convert your model output to the expected format for
    another assimilation by DART. DART will write out a control file
    that contains some information that must be passed to
    *advance_model.csh*: for example, the number of ensemble members,
    the input and output filenames for each ensemble member, etc.
    1.  Prepare a directory (or multiple directories) with the contents
        needed to advance your model.
    2.  Modify the input to your model communicating the run-time
        settings necessary to integrate your model from one time to
        another *arbitrary* time in the future.
    3.  Convert (if necessary) your input file to netCDF.
    4.  Run the model (you may need to watch the MPI syntax)
    5.  Convert (if necessary) the model output to a netCDF file DART
        can use for the next assimilation.
4.  *If a single instance of your model needs to advance using all the
    MPI tasks*, there is one more script that needs to work -
    *shell_scripts/run_filter.csh*. This script must do quite a lot.
    Find some examples in the *models/\*/shell_scripts* directories.  
5.  \[optional step\] Modify the MATLAB® routines to know about the
    specifics of the netCDF files produces by your model (sensible
    defaults, for the most part.)  
6.  *Test*. Generally, it is a good strategy to use DART to create a
    synthetic observation sequence with *ONE* observation location - and
    *ONE* observation type - for several assimilation periods. With
    that, it is possible to run *perfect_model_obs* and then *filter*
    without having to debug too much stuff at once. A separate document
    will address how to test your model with DART.  

## Programming style

\#1 Don't shoot the messenger. We have a lot of experience trying to
write portable/reproducible code and offer these suggestions. All of
these suggestions are for the standalone DART components. We are not
asking you to rewrite your model. If your model is a separate
executable, leaving it untouched is fine. **Writing portable code for
the DART components will allow us to include your model in the nightly
builds and reduces the risk of us making changes that adversely affect
the integration with your model.** There are some routines that have to
play with the core DART routines, these are the ones we are asking you
to write using these few simple guidelines.

  - Use explicit typing, do not use or rely on the 'autopromote' flag on
    your compiler.
  - Use the *intent()* attribute.
  - Use the *use, xxx_mod, only : bob, sally* statements for routines
    from other modules. This really helps us track down things and
    ensures you're using what you think you're using.
  - Use Fortran namelists for I/O if possible.
  - Check out the existing parameters/routines in
    *assimilation_code/modules/utilities/types_mod.f90*,
    *assimilation_code/modules/utilities/utilities_mod.f90*, and
    *assimilation_code/modules/utilities/time_manager_mod.f90*. You
    are free to use these and are encouraged to do so. No point
    reinventing the wheel and these routines have been tested
    extensively.

Hopefully, you have no idea how difficult it is to build each model with
'unique' compile options on N different platforms. Fortran90 provides a
nice mechanism to specify the type of variable, please do not use
vendor-specific extensions. (To globally autopromote 32bit reals to
64bit reals, for example. That is a horrible thing to do, since vendors
are not consistent about what happens to explicitly-typed variables.
Trust me. They lie. It also defeats the generic procedure interfaces
that are designed to use a single interface as a front-end to multiple
'type-specific' routines.) Compilers do abide by the standard, however,
so DART code looks like:

~~~ 
   character(len=8)      :: crdate
   integer, dimension(8) :: values
   ...
   real(r4) :: a,b
   real(r8) :: bob
   integer  :: istatus, itype
   ...
   real(r8),            intent(in)  :: x(:)
   type(location_type), intent(in)  :: location
   integer,             intent(in)  :: itype
   integer,             intent(out) :: istatus
   real(r8),            intent(out) :: obs_val  
~~~

depending on the use. The *r4* and *r8* types are explicitly defined in
*assimilation_code/modules/utilities/types_mod.f90* to accurately
represent what we have come to expect from 32bit and 64bit floating
point real variables, respectively. If you like, you can redefine *r8*
to be the same as *r4* to shrink your memory requirement. The people who
run with WRF frequently do this. Do not redefine the *digits12*
parameter, that one must provide 64bit precision, and is used in
precious few places.

-----

<span id="addingAmodelSpecific" class="anchor"></span>

## Adding a model to DART - Specifics

If your model is a separate executable, it would be wise to look at the
heavily commented template script
[models/template/shell_scripts/advance_model.csh](https://github.com/ncar/dart/models/template/shell_scripts/advance_model.csh)
and then a few higher-order models to see how they do it. Become
familiar with [DART's use of
MPI](dart_mpi.html), the [options for
parallelism](filter_async_modes.html), and
the *filter* namelist parameter
[*async*](https://ncar.github.io/DART/api/v2.1.10/program/filter.html).

<span id="Copying" class="anchor"></span>

### 1. Copying the template directory

A little explanation/motivation is warranted. If the model uses the
standard layout, it is much easier to include the model in the nightly
builds and testing. For this reason alone, please try to use the
recommended directory layout. Simply looking at the *DART/models*
directory should give you a pretty good idea of how things should be
laid out. Copy the *template* directory and its contents. 
The point of copying this directory is to get a ```model_mod.f90``` that works
as-is and you can modify/debug the routines one at a time.

The destination directory (your model directory) should be in the
*DART/models* directory to keep life simple. Moving them around will
cause problems for the ```work/mkmf_xxxxx``` configuration files. Each
model directory should have a *work* and *shell_scripts* directories,
and may have a *matlab* directory, a *src* directory, or anything else
you may find convenient.  
  
Now, you must change all the ```work/path_names_xxxxx``` file contents to
reflect the location of your ```model_mod.f90```.

<span id="model_mod" class="anchor"></span>

### 2. model_mod.f90

We have templates, examples, and a document describing the required
interfaces in the DART code tree -
[DART/models/template/model_mod.html](Manhattan/models/template/model_mod.html).
Every(?) user-visible DART program/module is intended to have a matching
piece of documentation that is distributed along with the code. The DART
code tree always has the most current documentation.  
  
Check out *time_manager_mod.f90* and *utilities_mod.f90* for
general-purpose routines ...  
  
Use Fortran namelists for I/O if possible.  
  
Modify the *model_mod.f90* file to return specifics about your model.
This module MUST contain all the required interfaces (no surprise) but
it can also contain many more interfaces as is convenient. This module
should be written with the understanding that print statements and error
terminations will be executed by multiple processors/tasks. To restrict
print statements to be written once (by the master task), it is
necessary to preface the print as in this example:  
~~~
if (do_output()) write(*,*)'model_mod:namelist cal_NML',startDate_1,startDate_2
~~~

<span id="requiredinterfaces" class="anchor"></span>

#### Required Interfaces in model_mod.f90

No matter the complexity of the model, the DART software requires a few
interface routines in a model-specific Fortran90 module *model_mod.f90*
file. The *models/template/model_mod.f90* file has extended comment
blocks at the heads of each of these routines that go into much more
detail for what is to be provided. **You cannot change the types or
number of required arguments to any of the required interface
routines.** You can add optional arguments, but you cannot go back
throught the DART tree to change the gazillion calls to the mandatory
routines. It is absolutely appropriate to look at existing models to get
ideas about how to implement the interfaces. Finding a model
implementation that is functionally close to yours always helps.  
  
The table of the mandatory interfaces and expected programming
degree-of-difficulty
is:

| subroutine callable | separate executable | routine | description |
| ------------------- | ------------------- | ------- | ----------- |
| easy | easy | [get_model_size](Manhattan/models/template/model_mod.html#get_model_size) | This function returns the size of all the model variables (prognostic or diagnosed or ...) that are packed into the 1D DART state vector. That is, it returns the length of the DART state vector as a single scalar integer. |
| depends | trivial | [adv_1step](Manhattan/models/template/model_mod.html#adv_1step) | For subroutine-callable models, this routine is the one to actually advance the model 1 timestep (see *models/bgrid_solo/model_mod.f90* for an example). For non-subroutine-callable models, this is a NULL interface. Easy. |
| depends | depends | [get_state_meta_data](Manhattan/models/template/model_mod.html#get_state_meta_data) | This routine takes as input an integer into the DART state vector and returns the associated location and (optionally) variable type from *obs_kind/obs_kind_mod.f90*. (See *models/\*/model_mod.f90* for examples.) Since DART uses netCDF and is responsible for the storage order, this is generally pretty easy. |
| easy | easy | [shortest_time_between_assimilations](Manhattan/models/template/model_mod.html#shortest_time_between_assimilations) | This routine returns the smallest increment in time (in seconds) that the model is capable of advancing the state in a given implementation. For example, the dynamical timestep of a model is 20 minutes, but there are reasons you don't want to (or cannot) restart at this interval and would like to restart AT MOST every 6 hours. For this case, *shortest_time_between_assimilations* should return *21600*, i.e. 6\*60\*60. This is also interpreted as the nominal **assimilation period**. This interface is required for all applications. |
| easy | easy | [end_model](Manhattan/models/template/model_mod.html#end_model) | Performs any shutdown and cleanup needed. Good form would dictate that you should deallocate any storage allocated when you instantiated the model (from *static_init_model*, for example). |
| depends | depends | [static_init_model](Manhattan/models/template/model_mod.html#static_init_model) | Called to do one-time initialization of the model. This generally includes setting the grid information, calendar, etc. |
| trivial | trivial | [init_time](Manhattan/models/template/model_mod.html#init_time) | Returns a time that is somehow appropriate for starting up a long integration of the model **IFF** the *&perfect_model_obs_nml* namelist parameter *read_input_state_from_file = .false.* If this option is not to be used in *perfect_model_obs*, this can be a routine that simply throws a fatal error. |
| easy | easy | [init_conditions](Manhattan/models/template/model_mod.html#init_conditions) | Companion interface to *init_time*. Returns a model state vector that is somehow appropriate for starting up a long integration of the model. Only needed **IFF** the *&perfect_model_obs_nml* namelist parameter *read_input_state_from_file = .false.* |
| trivial - difficult | trivial - difficult | [nc_write_model_atts](Manhattan/models/template/model_mod.html#nc_write_model_atts) | This routine is used to write the model-specific attributes to netCDF files created by DART. If you are simply updating existing (template) netCDF files, this routine is very easy. The subroutine in the *models/template/model_mod.f90* WILL WORK for new models but does not know anything about prognostic variables or geometry or ... Still, it is enough to get started without doing anything. More meaningful coordinate variables etc. are needed to supplant the default template. This can be as complicated as you like - see existing models for examples. |
| trivial - medium | trivial - medium | [nc_write_model_vars](Manhattan/models/template/model_mod.html#nc_write_model_vars) | This routine is currently unused but anticipated for future enhancements. The default routine in *default_model_mod* should be referenced. |
| trivial | trivial | [get_close_obs](Manhattan/models/template/model_mod.html#get_close_obs) | This is the routine that takes a single observation location and a list of other observation locations, returns the indices of all observation locations close to the single observation along with the number and the distances for the close ones. This is generally a 'pass-through' routine to a routine of the same name in the location module. |
| trivial | trivial | [get_close_state](Manhattan/models/template/model_mod.html#get_close_state) | This is the routine that takes a single observation location and a list of state locations, returns the indices of all the state locations close to the observation as well as the number and the distances for the close ones. This is generally a 'pass-through' routine to a routine of the same name in the location module. |
| depends | hard | [model_interpolate](Manhattan/models/template/model_mod.html#model_interpolate) | This is one of the more difficult routines. Given a DART state vector, a location, and a desired generic 'quantity' (like QTY_SURFACE_PRESSURE, QTY_TEMPERATURE, QTY_SPECIFIC_HUMIDITY, QTY_PRESSURE, ... ); return the desired scalar quantity and set the return status accordingly. This is what enables the model to use observation-specific 'forward operators' that are part of the common DART code. |
| depends | trivial - difficult | [convert_vertical_obs](Manhattan/models/template/model_mod.html#convert_vertical_obs) | If needed, the difficulty lies in the complexity of the model vertical coordinate system. |
| depends | trivial - difficult | [convert_vertical_state](Manhattan/models/template/model_mod.html#convert_vertical_state) | If needed, the difficulty lies in the complexity of the model vertical coordinate system. |
| depends | trivial - difficult | [pert_model_copies](Manhattan/models/template/model_mod.html#pert_model_copies) | This routine is used to generate initial ensembles. This may be a NULL interface if you can tolerate the default perturbation strategy of adding noise to every state element **or** if you generate your own ensembles outside the DART framework. There are other ways of generating ensembles ... climatological distributions, bred singular vectors, voodoo ... |
| easy | easy | [read_model_time](Manhattan/models/template/model_mod.html#read_model_time) | This routine simply stores a copy of the ensemble mean of the state vector within the *model_mod*. The ensemble mean may be needed for some calculations (like converting model sigma levels to the units of the observation - pressure levels, for example). |
| easy | easy | [write_model_time](Manhattan/models/template/model_mod.html#read_model_time) | This routine simply stores a copy of the ensemble mean of the state vector within the *model_mod*. The ensemble mean may be needed for some calculations (like converting model sigma levels to the units of the observation - pressure levels, for example). |

### If your model is subroutine-callable - you're done\!

<span id="add_a_simple_model" class="anchor"></span>

## Adding support for a "simple" model

### overview

each routine includes usual things it often has to do for subroutine-callable
models which can manufacture an initial condition state vector.
[model_mod_check.f90](https://ncar.github.io/DART/api/v2.1.10/program/model_mod_check.html)
can be used to test these routines individually before you run it with *filter*.
start with all defaults from other modules and add, in order the following
routines:

1. init_time()

2. init_conditions()  
for a "cold start" fill in an empty state vector with initial conditions and set the initial time. if the state vector is all 0s and the time is 0, you can use the default routines.

3. get_model_size()  
return number of items in the state vector 

**if you have only a single type of variable in your state vector, use the next two:**

4. static_init_model()  
often your model_size is set by namelist. allocate an array of that size and precompute all the locations for each state vector item. call add_domain() with the model size so dart knows how long the state vector is.

5. get_state_meta_data()  
return QTY_STATE_VARIABLE as the quantity if present, and return the location for that index by looking it up in a location array. 

**if you have more than a single type of variable in the state vector:**

4. static_init_model()  
read namelist to see how many fields are going to be read in for the state
vector. use add_domain() to indicate which netcdf vars should be read. read in
any auxiliary data needed by interpolation code (eg. topology). read template
file to set grid locations. use get_domain_size() to compute model_size.

5. get_state_meta_data()  
call get_model_variable_indices() and get_state_kind() to figure
out the i,j,k indices and which variable this offset is. use the i,j,k
to compute the grid location and return it along with the quantity.

**now continue**

6. end_model()  
deallocate any arrays allocated in static_init_model() 

**at this point you can assimilate identity obs at the model time**

7. adv_1step()  
if possible, embed the code that computes x(t+1) = F(x(t)). or call a
separate subroutine to advance the model state from one time to another.

8. shortest_time_between_assimilations()  
return a namelist or a fixed value for the minimum model advance time. 

**at this point you can assimilate a time series of identity obs**

9. model_interpolate()  
where the bulk of the work often is. this routine gets passed the location and quantity of the observation. find the indices which enclose that location and interpolate to get an array of expected values. 

**at this point you can assimilate obs at locations other than state vector points.**

10. nc_write_model_atts()  
add attributes to the output diagnostic files. 

**anything below here generally can use the default routines in other modules:**

11. read_model_time()

12. write_model_time()  
generally can use the system defaults, unless you have a model restart file that already stores time in a particular format.

13. pert_model_copies()  
the default is to add gaussian noise to the entire model state. if you
want to only perturb a single variable, or perturb it with different
noise ranges you can add code here. used to generate an ensemble from a
single model state for *filter*.

14. convert_vertical_obs()

15. convert_vertical_state()  
unused in models without vertical coordinate choices

16. get_close_obs()

17. get_close_state()  
often unused unless you want to modify the localization behavior

18. nc_write_model_vars()  
not currently called, leave it using the default routine. here for
possible future implementation.

<span id="add_a_complex_model" class="anchor"></span>

## Adding support for a "complex" model

### overview

Each routine includes usual things it often has to do for a large
geophysical model. this is different from the low order models.
[model_mod_check.f90](https://ncar.github.io/DART/api/v2.1.10/program/model_mod_check.html)
can be used to test these routines individually before you run it with
*filter*. start with all defaults from other modules and add, in order
the following routines:

1. static_init_model()  
read namelist to see how many fields are going to be read in for the
state vector. use add_domain() to indicate which netcdf vars should be
read. read in any auxiliary data needed by interpolation code (eg.
topology). read template file to set grid locations. use
get_domain_size() to compute model_size.

2. end_model()  
deallocate any arrays allocated in static_init_model()

3. get_model_size()  
return model_size computed in static_init_model()

4. shortest_time_between_assimilations()  
return a namelist or a fixed value for the minimum model advance time.

5. read_model_time()

6. write_model_time()  
if the time is stored in the netcdf files, supply routines that can read
and write it in the correct format. we have default routines that work
if it matches what those routines expect: a time variable with an
optional calendar variable. if none, it's fractional days. if the time
variable is an array, read/write the last one.

7. get_state_meta_data()  
call get_model_variable_indices() and get_state_kind() to figure
out the i,j,k indices and which variable this offset is. use the i,j,k
to compute the grid location and return it along with the quantity.

**now you can assim identity obs**

8. model_interpolate()  
where the bulk of the work will be. get the location and quantity of the
observation. find the i,j,k indices which enclose that location, or
search for the cell number. can compute i,j in a regular lat/lon grid,
have to search in a deformed grid. if multiple vertical options,
different ensemble members may result in more than a single level. use
get_state() to get the ensemble-sized array of values for each offset
into the state vector, and do interpolation to get an array of expected
values.

**other obs**

9. nc_write_model_atts()  
can leave for later. eventually add grid info to the diag files for
plotting.

10. convert_vertical_obs()

11. convert_vertical_state()  
if this model has a choice of multiple vertical coordinates (e.g.
pressure, height, etc) add code here to convert between the possible
verticals.

12. get_close_obs()

13. get_close_state()  
if you want to change the impact based on something other than the type
or kind, put code here. should test for vertical and do the conversion
on demand if it hasn't already been done.

14. pert_model_copies()  
the default is to add gaussian noise to the entire model state. if you
want to only perturb a single variable, or perturb it with different
noise ranges you can add code here. used to generate an ensemble from a
single model state for *filter*.

15. init_time()

16. init_conditions()

17. adv_1step()  
generally not used in large geophysical models, but if you can generate
a single model state without reading in a file, supply code in
init_conditions. if you can advance the model via a subroutine, add the
code to adv_1step.

18. nc_write_model_vars() 
not currently called, leave it using the default routine. here for
possible future implementation.

<span id="cycling_a_model" class="anchor"></span>

## Cycling Models

### Overview

For simple models which can be advanced by a subroutine call, the
*filter* program is driven by the input observation sequence. It
assimilates all observations in the current assimilation window and then
advances the model state until the window includes the next available
observation. When it runs out of observations, *filter* exits.  
  
For complex models which are themselves an MPI program or have
complicated scripting to run the model here are some simplified
considerations for scripting an experiment. A "cycling script" would
need to:

1.  \[if needed\] Run the ensemble of models forward to the time of the
    first observation.
2.  The input observation sequence file should be created (or trimmed)
    to only include observations in the current window.
3.  Run *filter* to assimilate all observations in the current window.
4.  Save a copy of the output files and diagnostic files. Often a
    timestamp is used as part of the filename or subdirectory name to
    make it unique.
5.  Run the ensemble of models forward in time.
6.  Run *filter* again.
7.  Repeat until all observations have been assimilated.

### In More Detail

The *filter* program requires an ensemble of model output files in
NetCDF file format as input. If the model does not use NetCDF a
translation step from the model native format to NetCDF is needed. The
files are often named with the ensemble number as part of the name and
also with a timestamp as part of the filename or part of a subdirectory
name which contains all the files for that timestep. Symbolic links can
be used to link a common simpler name to a file with a timestamp in the
filename or directory name.  
  
The *filter* program also requires an input observation sequence file.
Often these are named with a timestamp to indicate the central time of
the observations, e.g. *obs_seq.2010-10-04.00:00:00* and then a common
name (e.g. *obs_seq.out*) is used with a symbolic link to indicate the
right file for input.  
  
If adaptive inflation is being used the *filter* program also requires
inflation input files. Again, timestamps in the names with a common
symbolic link name are often used here.  
  
The *filter* program runs.  
  
The output of the *filter* program include updated model files using one
of three different workflows:

1.  The *filter* program directly overwrites the input files.
      - Advantages: uses the least amount of disk usage and minimizes
        file copying.
      - Disadvantages: if something crashes the files can be left in an
        indeterminate state making restarting more complicated.
2.  The script copies the input files to the output names, and the
    *filter* program updates the existing files.
      - Advantages: The *filter* program can easily be restarted in case
        of problems because the original input files are unchanged. The
        output files are immediately available to be used as input to
        the model.
      - Disadvantages: uses more disk space.
3.  The *filter* program creates new output files from scratch.
      - Advantages: The output files are smaller since they only contain
        the state vector and no other grid or auxiliary information. The
        *filter* program can easily be restarted in case of problems.
      - Disadvantages: generally requires a post-processing step to
        insert the updated state information into full model restart
        files.

The script should also save the *obs_seq.final* diagnostic file,
possibly with a timestamp in the filename or subdirectory name, and the
updated inflation files in the case where adaptive inflation is used.  
  
The script can run the ensemble of models forward in time in many ways.
A few of the ways we're aware of are:

1.  If a queuing system is available, the ensemble of models can be
    submitted either as independent jobs or using the batch system's job
    array syntax. They run as soon as resources are available. The
    disadvantage is it can be complicated to know when all the jobs have
    finished successfully.
2.  On smaller clusters the ensemble members can be advanced one after
    the other in a loop. There is no question about when the last member
    has been advanced and it requires no more resources than running a
    single copy of the model. The disadvantage is this is the slowest
    wall-clock way to advance the ensemble.

<span id="scriptingbackground" class="anchor"></span> <span id="advance_model" class="anchor"></span>

### 3. advance_model.csh

### The Big Picture for models advanced as separate executables.

The normal sequence of events is that DART reads in its own restart file
(do not worry about where this comes from right now) and eventually
determines it needs to advance the model. DART needs to be able to take
its internal representation of each model state vector, the valid time
of that state, and the amount of time to advance the state - and
communicate that to the model. When the model has advanced the state to
the requested time, the output must be ingested by DART and the cycle
begins again. DART is entirely responsible for reading the observations
and there are several programs for creating and manipulating the
observation sequence files.  
  
There are a couple of ways to exploit parallel architectures with DART,
and these have an immediate bearing on the design of the script(s) that
control how the model instances (each model *copy*) are advanced.
Perhaps the conceptually simplest method is when each model instance is
advanced by a single processor element. DART calls this *async = 2*. It
is generally efficient to relate the ensemble size to the number of
processors being used.  
  
The alternative is to advance every model instance one after another
using all available processors for each instance of the model. DART
calls this *async = 4*, and requires an additional script. For
portability reasons, DART uses the same processor set for both the
assimilation and the model advances. For example, if you advance the
model with 96 processors, all 96 processors will be employed to
assimilate. If your model requires 2000 processors, all 2000 will be
employed for the assimilation. Some people exploit the queueing systems
on their large machines to allow for the explicit customization of how
many tasks are used for each model advance and for an assimilation.  
  
*advance_model.csh* is invoked in one of two ways: 1) if *async = 2*
then *filter* uses a *system()* call, or 2) if *async = 4* then
*run_filter.csh* makes the call. Either way there are three arguments.

1.  the process number of the caller - could be the master task ID
    (zero) or (especially if *async = 2*) a process id that gets related
    to the copy. When multiple copies are being advanced simultaneously,
    each of the advances happens in its own run-time directory.
2.  the number of state copies belonging to that process
3.  the name of the (ASCII) filter_control_file for that process. The
    filter_control file contains the following information (one per
    line): the ensemble member, the name of the input file (containing
    the DART state vector), and the name of the output file from the
    model containing the new DART state vector. For example,  
      
    <div class="routine">
    1  
    assim_model_state_ic.0001  
    assim_model_state_ud.0001  
    2  
    assim_model_state_ic.0002  
    assim_model_state_ud.0002  
    ...
    </div>

#### async = 2 ... advancing many copies at the same time

Modify *shell_scripts/advance_model.csh* to:

1.  Collect all the input files needed to advance the model into a
    clean, temporary directory.
2.  Create a routine or set of routines to modify the input to your
    model communicating the run-time settings necessary to integrate
    your model from one time to another *arbitrary* time in the future.
    These routines are called in the *advance_model.csh* script. Every
    model is controlled differently, so writing detailed descriptions
    here is pointless.
3.  Determine how many tasks you have, and how many ensemble members you
    have. Determine how many 'batches' of ensemble members must be done
    to advance all of them. With 20 tasks and 80 ensemble members, you
    will need to loop 4 times, for example. clean, temporary directory,
    and
4.  Loop over the following steps - each loop advances one ensemble
    member:
    1.  If necessary, convert the DART (posterior) file into input for
        your model, After DART has assimilated the observations and
        created new (posterior) states, it may be necessary to
        post-process these states to impose model-specific limitations.
        There is nothing in the ensemble filter methodology that
        restricts posteriors to be physically meaningful (soil moistures
        could be slightly negative, for example), or that related
        quantities in the state have been conserved. Since each model
        handles these situations differently, it is up to the user to
        write any post-processing routines that may be necessary to
        check for viable input to the model. Frequently this is done by
        a program called *dart_to_model.f90*. If you need to do this,
        you will also need to create/modify a *mkmf_dart_to_model*
        and *path_names_dart_to_model* specific to your model.
    2.  run your model, and
    3.  if necessary, convert your model output (the prior) to netCDF or
        simply rename or link to the appropriate filename for another
        assimilation by DART.

During this initial phase, it may be useful to _leave_ the temporary
directory intact until you verify everything is as it should be.

#### async = 4 ... advancing each copy one at a time

In addition to modifying *shell_scripts/advance_model.csh* as
described above, you must also modify *shell_scripts/run_filter.csh*
in the following way: THIS PART NEEDS TO BE FILLED IN

<span id="testing_strategies" class="anchor"></span>

### 4. Testing Strategies - under construction

Generally testing when you add a new model to DART includes:  
Checking the converter programs.  
Checking the model advance control.  
Starting with one observation in a known location, with a known value
and error specification.  
Performing a 'perfect model' experiment for a few timesteps.  
Looking at what the assimilation did with:

> work % ncdiff filter_output.nc preassim.nc Innov.nc
> work % ncview Innov.nc

<span id="model_matlab_support" class="anchor"></span>  
  
### 5. Adding MATLAB® support for your own model - under construction.

Only needed for state-space diagnostics.  
Define a structure with required elements.  
Examples exist in the *diagnostics/matlab/private* directory.

<span id="examples" class="anchor"></span> 

-----

## Examples - under construction

1.  observation location/value plots
2.  [a brief explanation of 'localization'](https://ncar.github.io/DART/api/v2.1.10/module/assim_tools_mod.html)
3.  namelist settings for damped adaptive spatially-varying group filter

<span id="namelists" class="anchor"></span> 

-----

## Namelists

Many DART programs have namelists to specify run-time control. Some
programs use one or more modules - each module may have its own
namelist. As a consequence, we find it convenient to have one file
(called ```input.nml```) specifying all the namelists.  
  
For a complete list of namelists, see the following document:
[master list of namelists](Manhattan/documentation/index.html#Namelists)  
  
An example namelist for each program is automatically built when the
makefile is generated by *mkmf_xxxxx*. The example namelist is named
```input.nml.xxxxx_default``` where *xxxxx* is the name of the program. The
example namelists have default values, which may not be appropriate for
your use. The default ```input.nml``` in each *work* directory generally has
better values. As usual, the documentation for each module is the best
place to get information about the namelist settings.

\[[top](#)\]

-----
