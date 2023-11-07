Downloading DART
================

The DART source code is distributed on the GitHub repository
`NCAR/DART <https://github.com/NCAR/DART>`_ with the documentation
served through :ref:`readthedocs <Welcome page>`.

Go to https://github.com/NCAR/DART and clone the repository or get the
ZIP file according to your preference. See the `github help page on
cloning <https://help.github.com/en/github/creating-cloning-and-archiving-repositories/cloning-a-repository>`_
for more information on how to clone a repository. Take note of the
directory you installed into, which is referred to as ``DART`` throughout 
this documentation.

To checkout the latest release of DART:

.. code:: 

   git clone https://github.com/NCAR/DART.git

If you have forked the DART repository, replace ``NCAR`` with your
Github username.

.. note::

   If you are interested in contributing to DART, see the
   :doc:`contributors-guide` for more information. In short, you
   will need to be familiar with the
   `GitHub workflow <https://guides.github.com/introduction/flow/>`_.


Unzip or clone the distribution in your desired directory, which we refer to as
``DART`` in this document. Compiling the code in this tree (as is usually the
case) may require a large amount of additional disk space (up to the 1 Gb
required for DART), so be aware of any disk quota restrictions before
continuing.

Organization of the repository
------------------------------

The top level DART source code tree contains the following directories and
files:

+------------------------+------------------------------------------------------------+
| Directory              | Purpose                                                    |
+========================+============================================================+
| ``assimilation_code/`` | assimilation tools and programs                            |
+------------------------+------------------------------------------------------------+
| ``build_templates/``   | Configuration files for installation                       |
+------------------------+------------------------------------------------------------+
| ``developer_tests/``   | regression testing                                         |
+------------------------+------------------------------------------------------------+
| ``diagnostics/``       | routines to diagnose assimilation performance              |
+------------------------+------------------------------------------------------------+
| ``guide/``             | General documentation and DART_LAB tutorials               |
+------------------------+------------------------------------------------------------+
| ``models/``            | the interface routines for the models                      |
+------------------------+------------------------------------------------------------+
| ``observations/``      | routines for converting observations and forward operators |
+------------------------+------------------------------------------------------------+
| ``theory/``            | pedagogical material discussing data assimilation theory   |
+------------------------+------------------------------------------------------------+
| **Files**              | **Purpose**                                                |
+------------------------+------------------------------------------------------------+
| ``CHANGELOG.rst``      | Brief summary of recent changes                            |
+------------------------+------------------------------------------------------------+
| ``copyright.rst``      | terms of use and copyright information                     |
+------------------------+------------------------------------------------------------+
| ``README.rst``         | Basic Information about DART                               |
+------------------------+------------------------------------------------------------+
