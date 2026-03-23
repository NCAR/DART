CrocoLake 
=========

This observation converter reads data from the CrocoLake, which is a database of oceanographic
observations developed and maintained in the framework of the NSF-sponsored project
`CROCODILE <https://github.com/CROCODILE-CESM>`__ (CESM Regional Ocean and Carbon 
cOnfigurator with Data assimilation and Embedding).

More details about CrocoLake can be found at the Woods Hole Oceanographic Institution
Biogeochemical Ocean Observing and Modeling Lab (`boom-lab <https://github.com/boom-lab>`__): `crocolake-python <https://github.com/boom-lab/crocolake-python>`__.

Required Python packages
------------------------

- dask[dataframe]
- gsw
- numpy
- pandas

To install the required packages (except for standard library modules):

.. code-block:: text

   pip install "dask[dataframe]" gsw numpy pandas


Example scripts
------------------

Two example Python scripts are provided to demonstrate how to select data from CrocoLake and 
convert it into DART observation sequence format.

.. Note:: 

   The example scripts assume that you have downloaded CrocoLake.
   The 'crocolake_path' in the scripts should be replaced with your own path to CrocoLake.


To run the example scripts from the command line:

.. code-block:: bash

   python3 crocolake_to_obsseq_example1.py
   python3 crocolake_to_obsseq_example2.py

The arguments for that can be passed to the CrocoLake class ObsSequence are:


.. code-block:: text


   crocolake_path (str):  path to desired CrocoLake database
   selected_vars (list):  list of variables to be extracted from the database
   db_filters (list):     list of db_filters to be applied to the database
   fill_na_qc (int):      replace value for NA in QC flags (default: None)
   fill_na_error (float): replace value for NA in error variables (default: None)
   obs_seq_out (str):     obs_seq file name
   loose (bool):          if True, store observation values also when
                          their QC and error are not present (default: False)
