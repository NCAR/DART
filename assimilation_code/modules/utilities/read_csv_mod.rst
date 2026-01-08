.. index :: csv

MODULE read_csv_mod
===================

.. contents::
   :depth: 2
   :local:

Overview
--------

The ``read_csv_mod`` module provides lightweight utility routines for reading
simple ASCII, CSV-like tabular data products.  It is intended for use in DART
observation converters and preprocessing workflows that need column-based
access to non-NetCDF data without additional external dependencies.

The module supports files with exactly one header row followed by data rows,
where fields are separated by a single-character delimiter (typically comma or
semicolon).

Public types
------------

.. rubric:: ``type(csv_file_type)``

A handle that stores cached information about an open CSV file.  All
components are private and must be accessed via module procedures.

.. container:: routine

   ::

      type csv_file_type
         private

         character(len=256) :: filename 
         integer            :: nrows    
         integer            :: ncols    
         integer            :: iunit    
         character          :: delim    
         character(len=512) :: fields(MAX_NUM_FIELDS)
         logical            :: is_open  
      end type csv_file_type

Public interfaces and routines
------------------------------

.. rubric:: ``csv_open(fname, cf, forced_delim, context)``

Opens a delimited text file, reads and caches the header, determines the number
of data rows, and initializes the CSV handle.

.. rubric:: ``csv_close(cf)``

Closes the file unit (if open) and resets the handle.

.. rubric:: ``csv_get_nrows(cf)``

Returns the number of data rows in the file (excluding the header).

.. rubric:: ``csv_get_field(cf, varname, varvals, context)``

Generic interface to read a column by header name.  Dispatches to type-specific
implementations for character, integer, or real output arrays.

.. note::
   
   - The output array must be sized exactly to ``csv_get_nrows(cf)``.
   - Header name matching is case-insensitive.
   - The file is rewound and the header skipped internally for each call.

.. rubric:: ``csv_get_field_index(cf, varname)``

Returns the 1-based column index of ``varname`` or -1 if the field is not found.

.. rubric:: ``csv_field_exists(cf, varname)``

Returns ``.true.`` if ``varname`` exists in the cached header.

.. rubric:: ``csv_print_header(cf)``

Prints the cached header fields and their indices to standard output.

.. rubric:: ``get_csv_words_from_string(inline, delim, wordcount, words)``

Splits a single line into delimiter-separated fields using CSV-like parsing
rules.

Typical usage
-------------

A CSV file is opened using a CSV-type handle. Output arrays are then allocated
using the number of data rows in the file. Data are read one column at a time
as follows:

.. code-block:: fortran

   use read_csv_mod, only : csv_file_type, csv_open, csv_close, csv_get_nrows, csv_get_field
   use types_mod,    only : r8

   type(csv_file_type) :: cf
   real(r8), allocatable :: lat(:)

   call csv_open('input.csv', cf)

   allocate(lat(csv_get_nrows(cf)))
   call csv_get_field(cf, 'lat', lat)

   call csv_close(cf)

What this module does
---------------------

- Reads delimited ASCII tables with:

  - one header row containing field names
  - one data record per subsequent line
  - a single-character delimiter

- Opens the file once and caches metadata in a ``csv_file_type`` handle:

  - filename and Fortran unit number
  - detected or user-specified delimiter
  - header field names and number of columns
  - number of data rows (excluding the header)

- Provides column-based access by header name using the generic interface
  ``csv_get_field()``, returning:

  - character strings (raw field contents)
  - integers (via ``string_to_integer``)
  - reals (via ``string_to_real``)

- Allows repeated access to different columns without reopening the file.
  Internally, the file is rewound and the header skipped for each read.

- Handles numeric conversion failures non-fatally:

  - values that cannot be converted to integer are returned as ``MISSING_I``
  - values that cannot be converted to real are returned as ``MISSING_R8``

- Treats backslash (``\``) as an escape character, preventing interpretation of the following character during parsing.

What this module does not do
----------------------------

- It does not implement the full CSV specification (e.g., RFC 4180).  The parser
  is intentionally simple and designed for common, well-behaved tabular files.

- It does not support:

  - multiple header lines
  - comment lines or metadata blocks
  - embedded newlines inside quoted fields

- It does not infer or mix data types within a single read.  Each call to
  ``csv_get_field()`` reads a *single column* into a single output type
  (character, integer, or real).  If a column contains mixed content
  (e.g., numeric and non-numeric values), numeric conversion will produce
  missing values for the non-conforming entries rather than raising an error.

- It does not read multiple columns in a single pass.  Each call to
  ``csv_get_field()`` rewinds and scans the file again, which may be inefficient
  when many columns are required.

Delimiter behavior: forced vs detected
--------------------------------------

By default, ``csv_open()`` detects the delimiter from the header line using a
simple heuristic that distinguishes between comma (`,`) and semicolon (`;`)
characters.  If semicolons occur more frequently than commas in the header,
the delimiter is assumed to be a semicolon; otherwise a comma is used.

If the optional ``forced_delim`` argument is provided to ``csv_open()``,
delimiter detection is skipped and the specified delimiter is used instead.
