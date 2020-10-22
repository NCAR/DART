###########
ROSE README
###########

Contents
========

#. `Overview`_
#. `Terms of Use`_

Overview
========

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

Terms of Use
============

|Copyright| University Corporation for Atmospheric Research

Licensed under the `Apache License, Version 2.0
<http://www.apache.org/licenses/LICENSE-2.0>`__. Unless required by applicable
law or agreed to in writing, software distributed under this license is
distributed on an "as is" basis, without warranties or conditions of any kind,
either express or implied.

.. |Copyright| unicode:: 0xA9 .. copyright sign
