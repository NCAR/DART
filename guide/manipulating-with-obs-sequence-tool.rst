Manipulating obs_seq files with the obs_sequence_tool
=====================================================

First and foremost, check out the
`obs_sequence_tool.html <../../assimilation_code/programs/obs_sequence_tool/obs_sequence_tool.html>`__
document for detailed information and examples.

*obs_sequence_tool* is the primary tool for manipulating observation sequence
files. Observations sequence files are linked lists of observations organized by
time. That is to say, the observations may appear in any order in the file, but
traversing the linked list will result in observations ordered by time.
*obs_sequence_tool* can be used to combine observation sequences, convert from
ASCII to binary or vice-versa, extract a subset of observations, etc.

For testing, it is terribly useful to extract a small number of observations
(like ONE) from an existing observation sequence file.
