Forward Operator
================

In Lanai the forward operator is performed by the first ens_size processors. This was because access to the whole state
vector is required for the forward operator, and only the first ens_size processors had the whole state vector. The
distributed state forward operator has a diffent loop structure to Lanai because all processors can do the foward
operator for their observations.

The forward operator is performed in ``get_obs_ens_distrb_state``. A limited call tree for ``get_obs_ens_distrb_state``
is shown below.

|image1|

The QC_LOOP is in ``get_obs_ens_distrb_state`` because the qc across the ensemble is known. This removes the need for a
transpose of the forward_op_ens_handle. Note this is different from Lanai. The window opening and closing in
``get_obs_ens_distrb_state`` is as follows:

#. State window created (processors can access other processor's memory)
#. Forward operator called
#. QC calculated
#. State window destroyed (processors can no longer access other processor's memory)

However, there may be occasions where having only the first ens_size processors perform the forward operator. For
example, if the forward operator is being read from a file, or the forward operator uses a large portion of the state.
Or when debugging it may be easier to have 1 task per ensemble member.

To transpose and do the forward operators like Lanai, you can use the filter_nml namelist option distribute_state =
.false. The process is the same as above except the window creation and destruction are transposing the state.

#. State window created (state ensemble is transposed var complete)
#. Forward operator called
#. QC calculated
#. State window destroyed (state ensemble is tranaposed to copy complete)

Note, that if you have fewer tasks than ensemble members some tasks will still be doing vectorized forward operators
(because they own more than one ensemble member).

State access
------------

Model_mod routines no longer get an array containing the state. The state is accessed through the function
``get_state``.

``x = get_state(i, state_handle)``

where x is the state at index i. ``state_handle`` is passed from above. During model_interpolate ``get_state`` returns
an array. Durring ``get_state`` returns a single value (the mean state).

.. |image1| image:: Graphs/forward_operator.gv.svg
