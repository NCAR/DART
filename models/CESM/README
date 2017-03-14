# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

2 April 2013 

This is the start of an interface for assimilating observations
into the fully-coupled CESM whole-system model.  It makes use
of the existing POP (ocean), CLM (land), and CAM (atmosphere) 
model_mod codes.

Mon Nov  9 13:44:50 MST 2015
Modified/updated 
Mar 10 2017
> > > See ./doc/setup_guidelines.html for much more information. < < < 

We have adopted some terminology to help us keep things straight.
1) CESM already uses the term 'fully-coupled', so we use that in
   reference to CESM components only.  It means that CESM has
   an active atmosphere and ocean, ignoring other components.  
   In CESM an active atmosphere almost always implies 
   an active land, but that is not necessary for it 
   to be called 'fully coupled', and, by itself, is not 'fully coupled'.
2) We use the term 'single-component' to denote the
   situation in which the assimilations are performed for
   ONE active model component. Atmospheric obs directly impact 
   the atmosphere, OR land obs directly impact the land, etc..
   Any impact from the atmosphere to the land
   happens through interaction with the CESM coupler.
3) We use the term 'multi-component' to denote the
   situation in which the assimilations are performed separately for
   more than one active model component. Atmospheric obs directly impact 
   the atmosphere AND ocean obs directly impact the ocean, etc..
   Any impact from the atmosphere to the ocean
   happens through interaction with the CESM coupler.
4) 'cross-component' is used to specify the case
   when observations of one component can directly impact any/all of
   the other active components without going through the coupler.

Prior to 9 Nov 2015, models/CESM  had programs that were an attempt to
achieve the cross-component, fully-coupled data assimilation. Since
this is being implemented with the Remote Memory Access (RMA) strategy
that is not consistent with the current SVN trunk, the files that allow
that usage pattern are being removed from the SVN trunk.
However, there are scripts in ./shell_scripts which enable 
the 'multi-component' assimilation mode, which does not require
a models/CESM/model_mod.f90, since it uses a separate filter for
each component (cam-fv, pop, ...).

The models/CESM/work directory has nothing of use in it, since there 
are no programs to interact with a cross-component DART state vector
(a DART state that consists of atmosphere and/or ocean and/or land).

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
