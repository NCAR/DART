<!--

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

-->
<html>
<head>
<META http-equiv="Content-Type" content="text/html; charset=EUC-JP">
<title>Module field_manager_mod</title>
<link type="text/css" href="http://www.gfdl.noaa.gov/~fms/style/doc.css" rel="stylesheet">
<STYLE TYPE="text/css">
          .fixed {
            font-size:medium;
            font-family:monospace;
            border-style:none;
            border-width:0.1em;
            padding:0.1em;
            color:#663366;
          }
        </STYLE>
</head>
<body>
<a name="TOP"></a><font class="header" size="1"><a href="#PUBLIC INTERFACE">PUBLIC INTERFACE </a>~
          <a href="#PUBLIC DATA">PUBLIC DATA </a>~
          <a href="#PUBLIC ROUTINES">PUBLIC ROUTINES </a>~
          <a href="#NAMELIST">NAMELIST </a>~
          <a href="#DIAGNOSTIC FIELDS">DIAGNOSTIC FIELDS </a>~
          <a href="#ERROR MESSAGES">ERROR MESSAGES </a>~
          <a href="#REFERENCES">REFERENCES </a>~ 
          <a href="#NOTES">NOTES</a></font>
<hr>
<h2>module field_manager_mod</h2>
<a name="HEADER"></a>
<!-- BEGIN HEADER -->
<div>
<b>Contact:&nbsp;</b><a href="mailto:wfc@gfdl.noaa.gov">   William Cooke </a>,&nbsp;
    <a href="mailto:mh2@gfdl.noaa.gov">   Matthew Harrison </a>
<br>
<b>Reviewers:&nbsp;</b>
<br>
<b>Change History:&nbsp;</b><a href="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi">WebCVS Log home</a>
<br>
<br>
</div>
<!-- END HEADER -->
<a name="OVERVIEW"></a>
<hr>
<h4>OVERVIEW</h4>
<!-- BEGIN OVERVIEW -->
<p class="text">   The field manager reads entries from a field table and stores this information along with the type 
   of field it belongs to. This allows the component models to query the field manager to see if 
   non-default methods of operation are desired. In essence the field table is a powerful type of namelist.
   Default values can be provided for all the fields through a namelist, individual fields can be modified 
   through the field table however. </p>
<!-- END OVERVIEW -->
<a name="DESCRIPTION"></a>
<!-- BEGIN DESCRIPTION -->
<div>   The field table consists of entries in the following format.
   <br>
<br>
   The first line of an entry should consist of three quoted strings.
   The first quoted string will tell the field manager what type of field it is.
   At present the supported types of fields are 
   "tracer" for tracers, 
   "xland_mix" for cross-land mixing, and,
   "checkerboard" for checkerboard null mode.
   <br>
<br>
   The second quoted string will tell the field manager which model the field is 
   being applied to.
   The supported types at present are
   "atmos_mod" for the atmosphere model,
   "ocean_mod" for the ocean model,
   "land_mod" for the land model, and,
   "ice_mod" for the ice model.
   <br>
<br>
   The third quoted string should be a unique name that can be used as a query.
   <br>
<br>
   The second and following lines of each entry are called methods in this context.
   Methods can be developed within any module and these modules can query the field manager to
   find any methods that are supplied in the field table.
   <br>
<br>
   These lines can consist of two or three quoted strings. The first string will be an identifier 
   that the querying module will ask for. The second string will be a name that the querying module 
   can use to set up values for the module. The third string, if present, can supply parameters to 
   the calling module that can be parsed and used to further modify values.
   <br>
<br>
   An entry is ended with a backslash (/) as the final character in a row.
   <br>
<br>
   Comments can be inserted in the field table by having a # as the first character in the line.
   <br>
<br>
   An example of a field table entry could be <pre>"tracer","atmos_mod","sphum"/
"tracer","atmos_mod","sf6"
"longname","sulf_hex"
"advection_scheme_horiz","2nd_order"
"Profile_type","Fixed","surface_value = 0.0E+00"/</pre>   In this example we have two field entries. 
   <br>
<br>
   The first is a simple declaration of a tracer called "sphum". 
   <br>
<br>
   The second is for a tracer called "sf6". Methods that are being applied to this tracer include
   initiating the long name of the tracer to be "sulf_hex", changing the horizontal advection scheme 
   to be second order, and initiating the tracer with a profile with fixed values, in this example all zero.
   <br>
<br> 
</div>
<br>
<!-- END DESCRIPTION -->
<a name="OTHER MODULES USED"></a>
<hr>
<h4>OTHER MODULES USED</h4>
<!-- BEGIN OTHER MODULES USED -->
<div>
<pre>   mpp_mod<br>mpp_io_mod<br>   fms_mod</pre>
</div>
<!-- END OTHER MODULES USED -->
<!-- BEGIN PUBLIC INTERFACE -->
<a name="PUBLIC INTERFACE"></a>
<hr>
<h4>PUBLIC INTERFACE</h4>
<div>
<pre>use field_manager_mod [, only:  field_manager_init,<br>                                field_manager_end,<br>                                find_field_index,<br>                                get_field_info,<br>                                get_field_method,<br>                                get_field_methods,<br>                                parse ]</pre>
<dl>
<dt>
<a href="#field_manager_init">field_manager_init</a>:</dt>
<dd>   Routine to initialize the field manager. </dd>
<dt>
<a href="#field_manager_end">field_manager_end</a>:</dt>
<dd>   Destructor for field manager. </dd>
<dt>
<a href="#find_field_index">find_field_index</a>:</dt>
<dd>   Function to return the index of the field. </dd>
<dt>
<a href="#get_field_info">get_field_info</a>:</dt>
<dd>   This routine allows access to field information given an index. </dd>
<dt>
<a href="#get_field_method">get_field_method</a>:</dt>
<dd>   A routine to get a specified method. </dd>
<dt>
<a href="#get_field_methods">get_field_methods</a>:</dt>
<dd>   A routine to obtain all the methods associated with a field. </dd>
<dt>
<a href="#parse">parse</a>:</dt>
<dd>   A function to parse an integer or an array of integers, 
   a real or an array of reals, a string or an array of strings. </dd>
</dl>
</div>
<br>
<!-- END PUBLIC INTERFACE -->
<a name="PUBLIC DATA"></a>
<hr>
<h4>PUBLIC DATA</h4>
<!-- BEGIN PUBLIC DATA -->
<div>
<table align="center" cellspacing="2" CELLPADDING="2" BORDER="2">
<tr>
<th> Name  </th><th> Type  </th><th> Value  </th><th> Units  </th><th> Description  </th>
</tr>
<tr>
<td> NUM_MODELS  </td><td> integer, parameter  </td><td> 5  </td><td> ---  </td><td>    Number of models.   </td>
</tr>
<tr>
<td> module_is_initialized  </td><td> logical  </td><td> .false.  </td><td> ---  </td><td>    Field manager is initialized.   </td>
</tr>
<tr>
<td> MODEL_ATMOS  </td><td> integer, parameter  </td><td> 1  </td><td> ---  </td><td>    Atmospheric model.   </td>
</tr>
<tr>
<td> MODEL_OCEAN  </td><td> integer, parameter  </td><td> 2  </td><td> ---  </td><td>    Ocean model.   </td>
</tr>
<tr>
<td> MODEL_LAND  </td><td> integer, parameter  </td><td> 3  </td><td> ---  </td><td>    Land model.   </td>
</tr>
<tr>
<td> MODEL_ICE  </td><td> integer, parameter  </td><td> 4  </td><td> ---  </td><td>    Ice model.   </td>
</tr>
</table>
<br>
</div>
<!-- END PUBLIC DATA -->
<a name="PUBLIC ROUTINES"></a>
<hr>
<h4>PUBLIC ROUTINES</h4>
<!-- BEGIN PUBLIC ROUTINES -->
<ol type="a">
<li>
<a name="field_manager_init"></a>
<h4>field_manager_init</h4>
<pre>
<b>call field_manager_init </b>(nfields, table_name)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   This routine reads from a file containing formatted strings. 
   These formatted strings contain information on which schemes are needed within
   various modules. The field manager does not initialize any of those schemes 
   however. It simply holds the information and is queried by the appropriate 
   module. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>table_name&nbsp;&nbsp;&nbsp;</tt></td><td>   The name of the field table. The default name is field_table. <br>&nbsp;&nbsp;&nbsp;<span class="type">[character, optional, dimension(len=128)]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>nfields&nbsp;&nbsp;&nbsp;</tt></td><td>   The number of fields. <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="field_manager_end"></a>
<h4>field_manager_end</h4>
<pre>
<b>call field_manager_end </b>
</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   This subroutine writes to the logfile that the user is exiting field_manager and 
   changes the initialized flag to false. </dd>
<br>
<br>
</dl>
</li>
<li>
<a name="find_field_index"></a>
<h4>find_field_index</h4>
<pre>value= <b>find_field_index</b> ( model, field_name )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   This function when passed a model number and a field name will 
   return the index of the field within the field manager. This index 
   can be used to access other information from the field manager. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>model&nbsp;&nbsp;&nbsp;</tt></td><td>   The number indicating which model is used. <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="get_field_info"></a>
<h4>get_field_info</h4>
<pre>
<b>call get_field_info </b>( n,fld_type,fld_name,model,num_methods )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   When passed an index, this routine will return the type of field, 
   the name of the field, the model which the field is associated and 
   the number of methods associated with the field. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>n&nbsp;&nbsp;&nbsp;</tt></td><td>   The field index. <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>fld_type&nbsp;&nbsp;&nbsp;</tt></td><td>   The field type. <br>&nbsp;&nbsp;&nbsp;<span class="type">[character, dimension(*)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>fld_name&nbsp;&nbsp;&nbsp;</tt></td><td>   The name of the field. <br>&nbsp;&nbsp;&nbsp;<span class="type">[character, dimension(*)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>model&nbsp;&nbsp;&nbsp;</tt></td><td>   The number indicating which model is used. <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>num_methods&nbsp;&nbsp;&nbsp;</tt></td><td>   The number of methods. <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="get_field_method"></a>
<h4>get_field_method</h4>
<pre>
<b>call get_field_method </b>( n,m,method )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   This routine, when passed a field index and a method index will 
   return the method text associated with the field(n) method(m). </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>n&nbsp;&nbsp;&nbsp;</tt></td><td>   The field index. <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>m&nbsp;&nbsp;&nbsp;</tt></td><td>   The method index. <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="get_field_methods"></a>
<h4>get_field_methods</h4>
<pre>
<b>call get_field_methods </b>( n,methods )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   When passed a field index, this routine will return the text 
   associated with all the methods attached to the field. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>n&nbsp;&nbsp;&nbsp;</tt></td><td>   The field index. <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="parse"></a>
<h4>parse</h4>
<pre>number = <b>parse</b> (text, label, value)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Parse is an integer function that decodes values from a text string.
   The text string has the form: "label=list" where "label" is an
   arbitrary user defined label describing the values being decoded,
   and "list" is a list of one or more values separated by commas.
   The values may be integer, real, or character.
   Parse returns the number of values decoded. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>text&nbsp;&nbsp;&nbsp;</tt></td><td>   The text string from which the values will be parsed. <br>&nbsp;&nbsp;&nbsp;<span class="type">[character(len=*)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>label&nbsp;&nbsp;&nbsp;</tt></td><td>   A label which describes the values being decoded. <br>&nbsp;&nbsp;&nbsp;<span class="type">[character(len=*)]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>value&nbsp;&nbsp;&nbsp;</tt></td><td>   The value or values that have been decoded. <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, real, character(len=*)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>parse&nbsp;&nbsp;&nbsp;</tt></td><td>   The number of values that have been decoded. This allows 
   a user to define a large array and fill it partially with 
   values from a list. This should be the size of the value array. <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
</ol>
<!-- END PUBLIC ROUTINES -->
<a name="NAMELIST"></a>
<!-- BEGIN NAMELIST -->
<!-- END NAMELIST --><a name="DIAGNOSTIC FIELDS"></a>
<!-- BEGIN DIAGNOSTIC FIELDS -->
<!-- END DIAGNOSTIC FIELDS --><a name="DATA SETS"></a>
<!-- BEGIN DATA SETS -->
<hr>
<h4>DATA SETS</h4>
<div>None.<br>
<br>
</div>
<!-- END DATA SETS -->
<a name="ERROR MESSAGES"></a>
<!-- BEGIN ERROR MESSAGES -->
<hr>
<h4>ERROR MESSAGES</h4>
<div>
<dl>
<dt>
<b>NOTE in field_manager_init</b>
</dt>
<dd>
<span class="errmsg">No field table available, so no fields are being registered.</span>
</dd>
<dd>   The field table does not exist. </dd>
<dt>
<b>FATAL in field_manager_init</b>
</dt>
<dd>
<span class="errmsg">max fields exceeded</span>
</dd>
<dd>   Maximum number of fields for this module has been exceeded. </dd>
<dt>
<b>FATAL in field_manager_init</b>
</dt>
<dd>
<span class="errmsg">Too many fields in tracer entry.</span>
</dd>
<dd>   There are more that 3 fields in the tracer entry. This is probably due
   to separating the parameters entry into multiple strings. 
   The entry should look like <br>   "Type","Name","Control1=XXX,Control2=YYY" <br>   and not like<br>   "Type","Name","Control1=XXX","Control2=YYY" </dd>
<dt>
<b>FATAL in field_manager_init</b>
</dt>
<dd>
<span class="errmsg">Maximum number of methods for field exceeded</span>
</dd>
<dd>   Maximum number of methods allowed for entries in the field table has been exceeded. </dd>
<dt>
<b>NOTE in field_manager_init</b>
</dt>
<dd>
<span class="errmsg">field with identical name and model name duplicate found, skipping</span>
</dd>
<dd>   The name of the field and the model name are identical. Skipping that field. </dd>
<dt>
<b>FATAL in field_manager_init</b>
</dt>
<dd>
<span class="errmsg">error reading field table</span>
</dd>
<dd>   There is an error in reading the field table. </dd>
<dt>
<b>FATAL in get_field_info</b>
</dt>
<dd>
<span class="errmsg">invalid field index</span>
</dd>
<dd>   The field index is invalid because it is less than 1 or greater than the 
   number of fields. </dd>
<dt>
<b>FATAL in get_field_method</b>
</dt>
<dd>
<span class="errmsg">invalid field index</span>
</dd>
<dd>   The field index is invalid because it is less than 1 or greater than the 
   number of fields. </dd>
<dt>
<b>FATAL in get_field_method</b>
</dt>
<dd>
<span class="errmsg">invalid method index</span>
</dd>
<dd>   The method index is invalid because it is less than 1 or greater than 
   the number of methods. </dd>
<dt>
<b>FATAL in get_field_methods</b>
</dt>
<dd>
<span class="errmsg">invalid field index</span>
</dd>
<dd>   The field index is invalid because it is less than 1 or greater than the 
   number of fields. </dd>
<dt>
<b>FATAL in get_field_methods</b>
</dt>
<dd>
<span class="errmsg">method array too small</span>
</dd>
<dd>   The method array is smaller than the number of methods. </dd>
</dl>
<br>
</div>
<!-- END ERROR MESSAGES -->
<hr>
<div align="right">
<font size="-1"><a href="#TOP">top</a></font>
</div>
</body>
</html>
