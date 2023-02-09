How to test your model_mod routines
===================================

The program ``model_mod_check.f90`` can be used to test model_mod routines
individually before running them with filter. Add ``model_mod_check``
to the list of programs in ``DART/models/your_model/work/quickbuid.sh`` to
build ``model_mod_check`` with ``quickbuild.sh``.

For more information on the tests in ``model_mod_check`` see :doc:`model_mod_check 
<../assimilation_code/programs/model_mod_check/model_mod_check>` 

For more information on quickbuild.sh see :ref:`DART build system`.


