diag_table_tk
=============

Overview
--------

.. container::

   The script diag_table_tk is a GUI written in Perl/Tk for building diagnostics tables, which are used by the
   :doc:`./models/bgrid_solo/fms_src/shared/diag_manager/diag_manager` for run-time specification of diagnostics.

.. container::

   The diagnostics table allows users to specify sampling rates and the choice of fields at run time. The table consists
   of comma-separated ASCII values and may be hand-edited. The preferred method of building a table is to use the
   provided GUI interface diag_table_tk. A default diag table is provided with each runscript.

   The table is separated into three sections.

   #. **Global section:** The first two lines of the table contain the experiment title and base date. The base date is
      the reference time used for the time units. The base date must be greater than or equal to the model start date.
      The date consists of six space-separated integers: year, month, day, hour, minute, and second.

   #. **File section:** File lines contain 6 fields - file name, output frequency, output frequency units, file format
      (currently only support NetCDF), time units and long name for time axis. The format is as follows:

      ::

         "file_name", output_freq, "output_freq_units", format, "time_units", "time_long_name"


         output_freq:  
                  > 0  output frequency in "output_units"
                  = 0  output frequency every time step
                  =-1  output frequency at end of run

         output_freq_units = units used for output frequency
                 (years, months, days, minutes, hours, seconds)

         format:   1 NetCDF

         time_units   = units used to label the time axis
                  (days, minutes, hours, seconds)

   #. **Field section:** Field lines contain 8 fields - module name, field name, output field name, file name, time
      sampling (for averaging, currently only support all timesteps), time average, other operations (currently not
      implemented) and pack value (1,2,4 or 8). The format is as follows:

      ::

         "module_name", "field_name", "output_name", "file_name" "time_sampling", 
         time_avg, "other_opts", packing

         module_name :  e.g. "atmos_mod", "land_mod"

         time_avg = .true. or .false.

         packing  = 1  double precision
                  = 2  float
                  = 4  packed 16-bit integers
                  = 8  packed 1-byte (not tested?)

Installation
------------

.. container::

   diag_table_tk requires the following perl modules:
   ::

      use English;
      use Tk;
      use Cwd;
      require Tk::FileSelect;
      require Tk::Text;
      use Tk::widgets qw/Dialog ErrorDialog ROText/;
      use Tk::FileDialog;
      use Tk::Balloon;
      use File::Find;

   Most of these are built by default in perl 5.004 and above; however, you may need to install the perl Tk modules.

   | Obtain Tk and Tk-FileDialog from:
   | http://www.cpan.org

   | Obtain Tk and Tk-FileDialog in RPM (red hat package manager) format from:
   | http://rpmfind.net/linux/rpm2html/search.php?query=perl-Tk
   | http://rpmfind.net/linux/rpm2html/search.php?query=perl-Tk-FileDialog

Usage
-----

#. **Load and edit previously saved diag tables.**
   Choose "Load Table" from the "File" menu.
#. **Quick parsing of f90 source code for fields which may be registered. Fields are grouped by module name.**
   To obtain a list of available diagnostic fields, choose "Modify Output Field Entry" from the main menu. Enter the
   path to the directory containing your source code, and click "Search". After the search is complete, you can look in
   the "Field" menu for the list of available fields.
#. **Easy table editing, including ability to delete or edit selected lines.**
   To edit the text of an entry, click the "Show Table" button, Select the entry you wish to edit by clicking on the
   "Entry List" button, then click "Edit Entry". A new window will open in which you can make changes. Click "Save
   Changes" when you are finished.
#. **Error checks to help ensure that your diag table will work properly.**
   Ensures proper spacing and formatting.
#. **Online Help is available.**
   Choose "Help" from the menubar.

Bugs and future plans
---------------------

.. container::

   The "cancel" button doesn't seem to work.
   The Show Table window should be opened by default when you click "modify table a" or "modify table b". Visual
   feedback is good.

   | It should warn you if you make changes and quit without saving.
