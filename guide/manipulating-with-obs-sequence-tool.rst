Manipulating obs_seq files with the obs_sequence_tool
=====================================================

Please see the
:doc:`../assimilation_code/programs/obs_sequence_tool/obs_sequence_tool`
document for detailed information and examples.

The ``obs_sequence_tool`` is the primary means to manipulate observation
sequence files.

Observations sequence files are linked lists of observations organized by
time. The observations may appear in any order in the file, but traversing the
linked list will result in observations ordered by time.

The ``obs_sequence_tool`` can be used to combine observation sequences, convert
from ASCII to binary or vice-versa, extract a subset of observations, etc.

When you are testing your DA application, you should use the
``obs_sequence_tool`` to extract one or a small number of observations from an
existing observation sequence file for assimilation. Testing your application 
using a small number of observations will allow you to test and troubleshoot 
problems much faster than performing a full-scale assimilation.
