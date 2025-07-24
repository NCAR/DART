.. _pytools:

=========================
DART Python tools
=========================

DART provides a set of Python tools that can be used to work with DART data and models.
These tools are located in the DART/pytools directory.

We recommend creating a `virtual environment <https://docs.python.org/3/library/venv.html>`__ 
to install DART Python tools. This helps to keep the DART Python tools and their dependencies 
isolated from other Python packages you may have installed.

.. code-block:: text

    python3 -m venv py-dart
    source py-dart/bin/activate

You can then install the DART Python tools by navigating to the DART/pytools directory and running:

.. code-block:: text

    cd DART/pytools
    pip install -r pyDART.txt


When you want to use the DART Python tools in a new session you will need to activate the virtual environment:

.. code-block:: text

    source py-dart/bin/activate