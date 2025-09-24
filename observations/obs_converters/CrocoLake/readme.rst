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


Run the example scripts from the command line:

.. code-block:: bash

   python3 crocolake_to_obsseq_example1.py
   python3 crocolake_to_obsseq_example2.py