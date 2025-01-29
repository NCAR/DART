Default QCEFF files in lorenz_96/work
======================================

Three default QCEFF files are included in lorenz_96/work:

#. eakf_qceff_table.csv: This is the same as the default behavior but provides an example 
   of a completed QCEFF file. All variables are assumed to be normally distributed and unbounded.
#. bnrhf_qceff_table.csv: This file uses the bounded normal rank histogram distribution for all 
   assimilation distributions; the bounds in this case are at plus and minus infinity. 
#. enkf_qceff_table.csv: This file applies the traditional ensemble Kalman filter (reference), 
   a stochastic algorithm whose observation space behavior cannot be represented in the QCEFF 
   context. It does normal distributions for all transforms. 

For example, the following modification to the algorithm_info namelist entry results in applying 
the bounded normal rank histogram algorithms, give it a try:

.. code-block:: text

   &algorithm_info_nml
      qceff_table_filename = 'bnrhf_qceff_table.csv'
   /

Try changing back to 80 ensemble members with the BNRHF
