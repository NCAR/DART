How to test your model_mod routines
===================================

The program ``model_mod_check.f90`` can be used to test the routines
individually before running them with *filter*. Add a ``mkmf_model_mod_check``
and ``path_names_model_mod_check`` to your ``DART/models/your_model/work``
subdirectory. You might find it helpful to consult another model matching your
model type (simple or complex). See the documentation for ``model_mod_check`` in
``DART/assimilation_code/programs/model_mod_check`` for more information on the
tests available.
