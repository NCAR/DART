ROSE
====

Overview
--------

The rose model is an atmospheric model for the Mesosphere Lower-Thermosphere
(MLT). The DART interface was developed by Tomoko Matsuo (now at CU-Boulder).

The source code for rose is not distributed with DART, thus the
DART/models/rose/work/workshop_setup.csh script is SUPPOSED to fail without the
rose code.

The rose model is a research model that is still being developed. The DART
components here are simply to help the rose developers with the DART framework.

As of Mon Mar 22 17:23:20 MDT 2010 the rose project has been substantially
streamlined. There is no need for the trans_time and build_nml routines.
dart_to_model has assumed those responsibilities.
