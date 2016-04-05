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
<title>Module tracer_manager_mod</title>
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
<h2>module tracer_manager_mod</h2>
<a name="HEADER"></a>
<!-- BEGIN HEADER -->
<div>
<b>Contact:</b>&nbsp;<a href="mailto:wfc@gfdl.noaa.gov">
   William Cooke</a>
<br>
<b>Reviewers:</b>&nbsp;<a href="mailto:mjh@gfdl.noaa.gov">
   Matt Harrison</a>,&nbsp;
    <a href="mailto:bw@gfdl.noaa.gov">
   Bruce Wyman</a>
<br>
<b>Change History: </b>&nbsp;<a href="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/">WebCVS Log</a>
<br>
<b>Last Modified:</b>&nbsp;2002/06/11 18:44:40<br>
<br>
</div>
<!-- END HEADER -->
<a name="OVERVIEW"></a>
<hr>
<h4>OVERVIEW</h4>
<!-- BEGIN OVERVIEW -->
<p class="text">
   Code to manage the simple addition of tracers to the FMS code.
   This code keeps track of the numbers and names of tracers included
   in a tracer table.</p>
<!-- END OVERVIEW -->
<a name="DESCRIPTION"></a>
<!-- BEGIN DESCRIPTION -->
<div>
   This code is a grouping of calls which will allow the simple
   introduction of tracers into the FMS framework. It is designed to
   allow users of a variety of component models interact easily with
   the dynamical core of the model. 
   <br>
<br>
   In calling the tracer manager routines the user must provide a
   parameter identifying the model that the user is working with. This
   parameter is defined within field_manager as MODEL_X 
   where X is one of [ATMOS, OCEAN, LAND, ICE].
   <br>
<br>
   In many of these calls the argument list includes model and tracer_index. These 
   are the parameter corresponding to the component model and the tracer_index N is 
   the Nth tracer within the component model. Therefore a call with MODEL_ATMOS and 5 
   is different from a call with MODEL_OCEAN and 5.
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
<pre>          mpp_mod<br>       mpp_io_mod<br>          fms_mod<br>field_manager_mod</pre>
</div>
<!-- END OTHER MODULES USED -->
<!-- BEGIN PUBLIC INTERFACE -->
<a name="PUBLIC INTERFACE"></a>
<hr>
<h4>PUBLIC INTERFACE</h4>
<div>
<pre>use tracer_manager_mod [, only:  tracer_manager_init,<br>                                 register_tracers,<br>                                 get_number_tracers,<br>                                 get_tracer_indices,<br>                                 get_tracer_index,<br>                                 assign_tracer_field,<br>                                 tracer_manager_end,<br>                                 get_tracer_field,<br>                                 get_tracer_tlevels,<br>                                 get_tracer_tendency,<br>                                 get_tracer_names,<br>                                 get_family_name,<br>                                 check_if_prognostic,<br>                                 find_family_members,<br>                                 add_members_to_family,<br>                                 split_family_into_members,<br>                                 set_tracer_profile,<br>                                 query_method,<br>                                 query_combined,<br>                                 set_tracer_atts ]</pre>
<dl>
<dt>
<a href="#tracer_manager_init">tracer_manager_init</a>:</dt>
<dd>
   Routine to initialize the tracer manager</dd>
<dt>
<a href="#register_tracers">register_tracers</a>:</dt>
<dd>
   A routine to register the tracers included in a component model.</dd>
<dt>
<a href="#get_number_tracers">get_number_tracers</a>:</dt>
<dd>
   A routine to return the number of tracers included in a component model.</dd>
<dt>
<a href="#get_tracer_indices">get_tracer_indices</a>:</dt>
<dd>
   Routine to return the component model tracer indices as defined within
   the tracer manager.</dd>
<dt>
<a href="#get_tracer_index">get_tracer_index</a>:</dt>
<dd>
   Function which returns the number assigned to the tracer name.</dd>
<dt>
<a href="#assign_tracer_field">assign_tracer_field</a>:</dt>
<dd>
   Routine to point the appropriate field within the tracer_type to the 
   appropriate field within the component model.</dd>
<dt>
<a href="#tracer_manager_end">tracer_manager_end</a>:</dt>
<dd>
   Routine to write to the log file that the tracer manager is ending.</dd>
<dt>
<a href="#get_tracer_field">get_tracer_field</a>:</dt>
<dd>
   A function to retrieve the present timestep data.</dd>
<dt>
<a href="#get_tracer_tlevels">get_tracer_tlevels</a>:</dt>
<dd>
   A function to retrieve the three time levels  data.</dd>
<dt>
<a href="#get_tracer_tendency">get_tracer_tendency</a>:</dt>
<dd>
   A function to retrieve the tendency data.</dd>
<dt>
<a href="#get_tracer_names">get_tracer_names</a>:</dt>
<dd>
   Routine to find the names associated with a tracer number.</dd>
<dt>
<a href="#get_family_name">get_family_name</a>:</dt>
<dd>
   Routine to return the family name for tracer n.</dd>
<dt>
<a href="#check_if_prognostic">check_if_prognostic</a>:</dt>
<dd>
   Function to see if a tracer is prognostic or diagnostic.</dd>
<dt>
<a href="#find_family_members">find_family_members</a>:</dt>
<dd>
   Subroutine to find which tracers are members of family family_name.</dd>
<dt>
<a href="#add_members_to_family">add_members_to_family</a>:</dt>
<dd>
   Routine to sum up the members of a family of tracers so that they may
   be advected and diffused as one tracer.</dd>
<dt>
<a href="#split_family_into_members">split_family_into_members</a>:</dt>
<dd>
   Subroutine that sets the present value of the member of a tracer 
   family according to the fraction of the family that it was in the 
   previous step.</dd>
<dt>
<a href="#set_tracer_profile">set_tracer_profile</a>:</dt>
<dd>
   Subroutine to set the tracer field to the wanted profile.</dd>
<dt>
<a href="#query_method">query_method</a>:</dt>
<dd>
   A function to query the "methods" associated with each tracer.</dd>
<dt>
<a href="#query_combined">query_combined</a>:</dt>
<dd>
   A function to query whether families of tracers have been combined already.</dd>
<dt>
<a href="#set_tracer_atts">set_tracer_atts</a>:</dt>
<dd>
   A subroutine to allow the user set the tracer longname and units from the 
   tracer initialization routine.</dd>
</dl>
</div>
<br>
<!-- END PUBLIC INTERFACE -->
<a name="PUBLIC DATA"></a>
<hr>
<h4>PUBLIC DATA</h4>
<!-- BEGIN PUBLIC DATA -->
<div>None.<br>
<br>
</div>
<!-- END PUBLIC DATA -->
<a name="PUBLIC ROUTINES"></a>
<hr>
<h4>PUBLIC ROUTINES</h4>
<!-- BEGIN PUBLIC ROUTINES -->
<ol type="a">
<li>
<a name="tracer_manager_init"></a>
<h4>tracer_manager_init</h4>
<pre>
<b>call tracer_manager_init </b>
</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>
   This routine writes the version and tagname to the logfile and 
   sets the module initialization flag.</dd>
<br>
<br>
</dl>
</li>
<li>
<a name="register_tracers"></a>
<h4>register_tracers</h4>
<pre>
<b>call register_tracers </b>(model, num_tracers,num_prog,num_diag,num_family)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>
   This routine returns the total number of valid tracers, the number of
   prognostic and diagnostic tracers and the number of families of
   tracers.</dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>model&nbsp;&nbsp;&nbsp;</tt></td><td>A parameter to identify which model is being used.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
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
<td valign="top" align="left"><tt>num_tracers&nbsp;&nbsp;&nbsp;</tt></td><td>The total number of valid tracers within the component model.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>num_prog&nbsp;&nbsp;&nbsp;</tt></td><td>The number of prognostic tracers within the component model.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>num_diag&nbsp;&nbsp;&nbsp;</tt></td><td>The number of diagnostic tracers within the component model.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>num_family&nbsp;&nbsp;&nbsp;</tt></td><td>The number of family tracers within the component model.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="get_number_tracers"></a>
<h4>get_number_tracers</h4>
<pre>
<b>call get_number_tracers </b>(model, num_tracers,num_prog,num_diag,num_family)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>
   This routine returns the total number of valid tracers, the number of
   prognostic and diagnostic tracers and the number of families of
   tracers.</dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>model&nbsp;&nbsp;&nbsp;</tt></td><td>A parameter to identify which model is being used.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
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
<td valign="top" align="left"><tt>num_tracers&nbsp;&nbsp;&nbsp;</tt></td><td>The total number of valid tracers within the component model.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, optional]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>num_prog&nbsp;&nbsp;&nbsp;</tt></td><td>The number of prognostic tracers within the component model.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, optional]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>num_diag&nbsp;&nbsp;&nbsp;</tt></td><td>The number of diagnostic tracers within the component model.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, optional]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>num_family&nbsp;&nbsp;&nbsp;</tt></td><td>The number of family tracers within the component model.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, optional]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="get_tracer_indices"></a>
<h4>get_tracer_indices</h4>
<pre>
<b>call get_tracer_indices </b>(model, ind, prog_ind, diag_ind, fam_ind)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>
   If several models are being used or redundant
   tracers have been written to the tracer_table, then the indices in
   the component model and the tracer manager may not have a one to one
   correspondence. Therefore the component model needs to know what index
   to pass to calls to tracer_manager routines in order that the correct
   tracer information be accessed.</dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>model&nbsp;&nbsp;&nbsp;</tt></td><td>A parameter to identify which model is being used.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
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
<td valign="top" align="left"><tt>ind&nbsp;&nbsp;&nbsp;</tt></td><td>An array containing the tracer manager defined indices for
         all the tracers within the component model.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, optional, dimension(:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>prog_ind&nbsp;&nbsp;&nbsp;</tt></td><td>An array containing the tracer manager defined indices for
              the prognostic tracers within the component model.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, optional, dimension(:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>diag_ind&nbsp;&nbsp;&nbsp;</tt></td><td>An array containing the tracer manager defined indices for
              the diagnostic tracers within the component model.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, optional, dimension(:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>fam_ind&nbsp;&nbsp;&nbsp;</tt></td><td>An array containing the tracer manager defined indices for
             the family tracers within the component model.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, optional, dimension(:)]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="get_tracer_index"></a>
<h4>get_tracer_index</h4>
<pre>value= <b>get_tracer_index</b> (model, name, indices, verbose)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>
   This is a function which returns the index, as implied within the component model.</dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>model&nbsp;&nbsp;&nbsp;</tt></td><td>A parameter to identify which model is being used.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>name&nbsp;&nbsp;&nbsp;</tt></td><td>The name of the tracer (as assigned in the field table).<br>&nbsp;&nbsp;&nbsp;<span class="type">[character]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>indices&nbsp;&nbsp;&nbsp;</tt></td><td>An array of the component model indices. This array can be found by 
             calling get_tracer_indices.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, optional, dimension(:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>verbose&nbsp;&nbsp;&nbsp;</tt></td><td>A flag to allow the message saying that a tracer with this name has not 
             been found. This should only be used for debugging purposes.<br>&nbsp;&nbsp;&nbsp;<span class="type">[logical, optional]</span></td>
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
<td valign="top" align="left"><tt>get_tracer_index&nbsp;&nbsp;&nbsp;</tt></td><td>The index of the tracer named "name". 
                      If indices is passed then the result is the array index which 
                      corresponds to tracer named "name".<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="assign_tracer_field"></a>
<h4>assign_tracer_field</h4>
<pre>
<b>call assign_tracer_field </b>(model,index, data, data_tlevels, tendency)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>
   The generality provided here is that one can point the three
   dimensional tracer field at either a two time level scheme [data and
   tendency] or a three time level scheme [data_tlevels]. The tracer manager
   points the appropriate tracer_type field at the data supplied from the 
   component model.</dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>model&nbsp;&nbsp;&nbsp;</tt></td><td>A parameter representing the component model in use.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>index&nbsp;&nbsp;&nbsp;</tt></td><td>The tracer number that you wish to assign a tracer
           field for.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>data&nbsp;&nbsp;&nbsp;</tt></td><td>The 3D field that is associated with the present time 
          step in the component model.<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, target, optional, dimension(:,:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>tendency&nbsp;&nbsp;&nbsp;</tt></td><td>The 3D field that is associated with the tendency time
              step in the component model.<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, target, optional, dimension(:,:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>data_tlevels&nbsp;&nbsp;&nbsp;</tt></td><td>The 4D field that is associated with the tracer field 
                  in the component model.<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, target, optional, dimension(:,:,:,:)]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="tracer_manager_end"></a>
<h4>tracer_manager_end</h4>
<pre>
<b>call tracer_manager_end </b>
</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>
   Routine to write to the log file that the tracer manager is ending.</dd>
<br>
<br>
</dl>
</li>
<li>
<a name="get_tracer_field"></a>
<h4>get_tracer_field</h4>
<pre>array= <b>get_tracer_field</b> (model, tracer_index)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>
   Function to point to the 3D field associated with a tracer.</dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>model&nbsp;&nbsp;&nbsp;</tt></td><td>A parameter representing the component model in use.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>tracer_index&nbsp;&nbsp;&nbsp;</tt></td><td>The tracer number within the component model.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
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
<td valign="top" align="left"><tt>data&nbsp;&nbsp;&nbsp;</tt></td><td>The tracer field is returned in this array.<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, pointer, dimension(:,:,:)]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="get_tracer_tlevels"></a>
<h4>get_tracer_tlevels</h4>
<pre>array= <b>get_tracer_tlevels</b> (model, tracer_index)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>
   Function to point to the 4D field associated with a tracer.</dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>model&nbsp;&nbsp;&nbsp;</tt></td><td>A parameter representing the component model in use.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>tracer_index&nbsp;&nbsp;&nbsp;</tt></td><td>The tracer number within the component model.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
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
<td valign="top" align="left"><tt>data&nbsp;&nbsp;&nbsp;</tt></td><td>The tracer field is returned in this array.<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, pointer, dimension(:,:,:,:)]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="get_tracer_tendency"></a>
<h4>get_tracer_tendency</h4>
<pre>array= <b>get_tracer_tendency</b> (model, tracer_index)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>
   Function to point to the 3D field associated with a tracer.</dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>model&nbsp;&nbsp;&nbsp;</tt></td><td>A parameter representing the component model in use.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>tracer_index&nbsp;&nbsp;&nbsp;</tt></td><td>The tracer number within the component model.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
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
<td valign="top" align="left"><tt>data&nbsp;&nbsp;&nbsp;</tt></td><td>The tracer tendency field is returned in this array.<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, pointer, dimension(:,:,:)]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="get_tracer_names"></a>
<h4>get_tracer_names</h4>
<pre>
<b>call get_tracer_names </b>(model,n,name,longname, units)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>
   This routine can return the name, long name and units associated
   with a tracer.</dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>model&nbsp;&nbsp;&nbsp;</tt></td><td>A parameter representing the component model in use.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>n&nbsp;&nbsp;&nbsp;</tt></td><td>Tracer number.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
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
<td valign="top" align="left"><tt>name&nbsp;&nbsp;&nbsp;</tt></td><td>Field name associated with tracer number.<br>&nbsp;&nbsp;&nbsp;<span class="type">[character]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>longname&nbsp;&nbsp;&nbsp;</tt></td><td>The long name associated with tracer number.<br>&nbsp;&nbsp;&nbsp;<span class="type">[character, optional]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>units&nbsp;&nbsp;&nbsp;</tt></td><td>The units associated with tracer number.<br>&nbsp;&nbsp;&nbsp;<span class="type">[character, optional]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="get_family_name"></a>
<h4>get_family_name</h4>
<pre>
<b>call get_family_name </b>(model,n,name)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>
   You may wish to use this routine to retrieve the name of the family
   that a tracer belongs to.</dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>model&nbsp;&nbsp;&nbsp;</tt></td><td>A parameter representing the component model in use.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>n&nbsp;&nbsp;&nbsp;</tt></td><td>Tracer number that you want the family name for.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
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
<td valign="top" align="left"><tt>name&nbsp;&nbsp;&nbsp;</tt></td><td>The family name.<br>&nbsp;&nbsp;&nbsp;<span class="type">[character]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="check_if_prognostic"></a>
<h4>check_if_prognostic</h4>
<pre>logical = <b>check_if_prognostic</b> (model, n)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>
   All tracers are assumed to be prognostic when read in from the field_table
   However a tracer can be changed to a diagnostic tracer by adding the line
   "tracer_type","diagnostic"
   to the tracer description in field_table.</dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>model&nbsp;&nbsp;&nbsp;</tt></td><td>A parameter representing the component model in use.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>n&nbsp;&nbsp;&nbsp;</tt></td><td>Tracer number that you want the family name for.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
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
<td valign="top" align="left"><tt>check_if_prognostic&nbsp;&nbsp;&nbsp;</tt></td><td>A logical flag set TRUE if the tracer is 
                         prognostic.<br>&nbsp;&nbsp;&nbsp;<span class="type">[logical]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="find_family_members"></a>
<h4>find_family_members</h4>
<pre>
<b>call find_family_members </b>(model, family_name,is_family_member)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>
   Subroutine to find which tracers are members of family family_name.
   This will return a logical array where the array positions 
   corresponding to the tracer numbers for family members are set .TRUE.</dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>model&nbsp;&nbsp;&nbsp;</tt></td><td>A parameter representing the component model in use.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>family_name&nbsp;&nbsp;&nbsp;</tt></td><td>The family name of the members one is seeking.<br>&nbsp;&nbsp;&nbsp;<span class="type">[character]</span></td>
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
<td valign="top" align="left"><tt>is_family_member&nbsp;&nbsp;&nbsp;</tt></td><td>A logical array where the tracer number is used as 
                      the index to signify which tracer is part of the family.
                      i.e. If tracers 1, 3, and 7 are part of the same family
                      then is_family_member(1), is_family_member(3), and 
                      is_family_member(7) are set TRUE.<br>&nbsp;&nbsp;&nbsp;<span class="type">[logical, dimension(:)]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="add_members_to_family"></a>
<h4>add_members_to_family</h4>
<pre>
<b>call add_members_to_family </b>(model,family_name, cur, prev, next)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>
   Routine to sum up the members of a family of tracers so that they may
   be advected and diffused as one tracer. This should only be used in
   conjunction with split_family_into_members and should be placed before
   the advection scheme is called.</dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>model&nbsp;&nbsp;&nbsp;</tt></td><td>A parameter representing the component model in use.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>n&nbsp;&nbsp;&nbsp;</tt></td><td>Tracer number.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>cur&nbsp;&nbsp;&nbsp;</tt></td><td>Array index for the current time step. This is only of use 
         with a three timestep model.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, optional]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>prev&nbsp;&nbsp;&nbsp;</tt></td><td>Array index for the previous time step. This is only of use 
          with a three timestep model.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, optional]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>next&nbsp;&nbsp;&nbsp;</tt></td><td>Array index for the next time step. This is only of use 
          with a three timestep model.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, optional]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>NOTE</b>
</dt>
<dd>
   This should be used with extreme caution. 
   Unless the family member distributions are similar to each other spatially, 
   advection as one tracer and subsequent splitting will result in a different
   result to advecting each tracer separately. The user should understand the 
   possible repercussions of this before using it.</dd>
<br>
<br>
</dl>
</li>
<li>
<a name="split_family_into_members"></a>
<h4>split_family_into_members</h4>
<pre>
<b>call split_family_into_members </b>(model,family_name,cur,prev,next)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>
   Subroutine that sets the present value of the member of a tracer 
   family according to the fraction of the family that it was in the 
   previous step.
   <br>
<br>
   This splits the transported family into the constituent members. This
   should only be used in conjunction with &lt;I&gt;add_members_to_family&lt;/I&gt; and
   should be placed after the advection scheme is called.
   <br>
<br>
</dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>model&nbsp;&nbsp;&nbsp;</tt></td><td>A parameter representing the component model in use.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>family_name&nbsp;&nbsp;&nbsp;</tt></td><td>The name of the family of tracers that you would 
                 like to split up.<br>&nbsp;&nbsp;&nbsp;<span class="type">[character]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>cur&nbsp;&nbsp;&nbsp;</tt></td><td>Array index for the current time step. This is only of use 
         with a three timestep model.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, optional]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>prev&nbsp;&nbsp;&nbsp;</tt></td><td>Array index for the previous time step. This is only of use 
          with a three timestep model.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, optional]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>next&nbsp;&nbsp;&nbsp;</tt></td><td>Array index for the next time step. This is only of use 
          with a three timestep model.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, optional]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>NOTE</b>
</dt>
<dd>
   This should be used with extreme caution. 
   Unless the family member distributions are similar to each other spatially, 
   advection as one tracer and subsequent splitting will result in a different
   result to advecting each tracer separately. The user should understand the 
   possible repercussions of this before using it.</dd>
<br>
<br>
</dl>
</li>
<li>
<a name="set_tracer_profile"></a>
<h4>set_tracer_profile</h4>
<pre>
<b>call set_tracer_profile </b>(model, n, surf_value, multiplier)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>
   If the profile type is 'fixed' then the tracer field values are set 
   equal to the surface value.
   If the profile type is 'profile' then the top/bottom of model and
   surface values are read and an exponential profile is calculated,
   with the profile being dependent on the number of levels in the
   component model. This should be called from the part of the dynamical
   core where tracer restarts are called in the event that a tracer
   restart file does not exist.
   <br>
<br>
   This can be activated by adding a method to the field_table
   e.g.
   "profile_type","fixed","surface_value = 1e-12"
   would return values of surf_value = 1e-12 and a multiplier of 1.0
   One can use these to initialize the entire field with a value of 1e-12.
   <br>
<br>
   "profile_type","profile","surface_value = 1e-12, top_value = 1e-15"
   In a 15 layer model this would return values of surf_value = 1e-12 and 
   multiplier = 0.6309573 i.e 1e-15 = 1e-12*(0.6309573^15)
   In this case the model should be MODEL_ATMOS as you have a "top" value.
   <br>
<br>
   If you wish to initialize the ocean model, one can use bottom_value instead
   of top_value.</dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>model&nbsp;&nbsp;&nbsp;</tt></td><td>A parameter representing the component model in use.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>n&nbsp;&nbsp;&nbsp;</tt></td><td>Tracer number.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
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
<td valign="top" align="left"><tt>surf_value&nbsp;&nbsp;&nbsp;</tt></td><td>The surface value that will be initialized for the tracer<br>&nbsp;&nbsp;&nbsp;<span class="type">[real]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>multiplier&nbsp;&nbsp;&nbsp;</tt></td><td>The vertical multiplier for the tracer
                Level(k-1) = multiplier * Level(k)<br>&nbsp;&nbsp;&nbsp;<span class="type">[real]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="query_method"></a>
<h4>query_method</h4>
<pre>logical = <b>query_method</b> (method_type, model, n, name, control)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>
   A function to query the "methods" associated with each tracer. The
   "methods" are the parameters of the component model that can be
   adjusted by user by placing formatted strings, associated with a
   particular tracer, within the field table.
   These methods can control the advection, wet deposition, dry
   deposition or initial profile of the tracer in question. Any
   parametrization can use this function as long as a routine for parsing
   the name and control strings are provided by that routine.</dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>method_type&nbsp;&nbsp;&nbsp;</tt></td><td>The method that is being requested.<br>&nbsp;&nbsp;&nbsp;<span class="type">[character]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>model&nbsp;&nbsp;&nbsp;</tt></td><td>A parameter representing the component model in use.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>n&nbsp;&nbsp;&nbsp;</tt></td><td>Tracer number that you want the family name for.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
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
<td valign="top" align="left"><tt>name&nbsp;&nbsp;&nbsp;</tt></td><td>A string containing the modified name to be used with
          method_type. i.e. "2nd_order" might be the default for 
          advection. One could use "4th_order" here to modify 
          that behaviour.<br>&nbsp;&nbsp;&nbsp;<span class="type">[character]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>control&nbsp;&nbsp;&nbsp;</tt></td><td>A string containing the modified parameters that are 
             associated with the method_type and name.<br>&nbsp;&nbsp;&nbsp;<span class="type">[character, optional]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>query_method&nbsp;&nbsp;&nbsp;</tt></td><td>A flag to show whether method_type exists with regard to
                  tracer n. If method_type is not present then one must
                  have default values.<br>&nbsp;&nbsp;&nbsp;<span class="type">[logical]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>NOTE</b>
</dt>
<dd>
   At present the tracer manager module allows the initialization of a tracer
   profile if a restart does not exist for that tracer. 
   Options for this routine are as follows
   <br>
<br>
   Tracer profile setup
   ==================================================================
   |method_type  |method_name  |method_control                      |
   ==================================================================
   |profile_type |fixed        |surface_value = X                   |
   |profile_type |profile      |surface_value = X, top_value = Y    |(atmosphere)
   |profile_type |profile      |surface_value = X, bottom_value = Y |(ocean)
   ==================================================================
   <br>
<br>
</dd>
<br>
<br>
</dl>
</li>
<li>
<a name="query_combined"></a>
<h4>query_combined</h4>
<pre>logical = <b>query_combined</b> (model, index)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>
   A function to query whether families of tracers have been combined already.
   This function should only be used in conjunction with add_members_to_family 
   and split_family_into_members.</dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>model&nbsp;&nbsp;&nbsp;</tt></td><td>A parameter representing the component model in use.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>index&nbsp;&nbsp;&nbsp;</tt></td><td>Tracer number.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
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
<td valign="top" align="left"><tt>query_combined&nbsp;&nbsp;&nbsp;</tt></td><td>A flag to show whether the tracer family has been combined.<br>&nbsp;&nbsp;&nbsp;<span class="type">[logical]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="set_tracer_atts"></a>
<h4>set_tracer_atts</h4>
<pre>
<b>call set_tracer_atts </b>(model, name, longname, units)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>
   A function to allow the user set the tracer longname and units from the 
   tracer initialization routine. It seems sensible that the user who is 
   coding the tracer code will know what units they are working in and it 
   is probably safer to set the value in the tracer code rather than in 
   the field table.</dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>model&nbsp;&nbsp;&nbsp;</tt></td><td>A parameter representing the component model in use.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>name&nbsp;&nbsp;&nbsp;</tt></td><td>Tracer name.<br>&nbsp;&nbsp;&nbsp;<span class="type">[character]</span></td>
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
<td valign="top" align="left"><tt>longname&nbsp;&nbsp;&nbsp;</tt></td><td>A string describing the longname of the tracer for output to NetCDF files<br>&nbsp;&nbsp;&nbsp;<span class="type">[character, optional]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>units&nbsp;&nbsp;&nbsp;</tt></td><td>A string describing the units of the tracer for output to NetCDF files<br>&nbsp;&nbsp;&nbsp;<span class="type">[character, optional]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>set_tracer_atts&nbsp;&nbsp;&nbsp;</tt></td><td>A flag to show that<br>&nbsp;&nbsp;&nbsp;<span class="type">[character, optional]</span></td>
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
<b>FATAL in register_tracers</b>
</dt>
<dd>
<span class="errmsg">invalid model type</span>
</dd>
<dd>
   The index for the model type is invalid.</dd>
<dt>
<b>NOTE in register_tracers</b>
</dt>
<dd>
<span class="errmsg">No tracers are available to be registered.</span>
</dd>
<dd>
   No tracers are available to be registered. This means that the field
   table does not exist or is empty.</dd>
<dt>
<b>FATAL in register_tracers</b>
</dt>
<dd>
<span class="errmsg">MAX_TRACER_FIELDS exceeded</span>
</dd>
<dd>
   The maximum number of tracer fields has been exceeded.</dd>
<dt>
<b>NOTE in register_tracers</b>
</dt>
<dd>
<span class="errmsg">There is only 1 tracer for tracer family X. Making an orphan.</span>
</dd>
<dd>
   A tracer has been given a family name but that family has only this member. Therefore it should be an orphan.</dd>
<dt>
<b>FATL in register_tracers</b>
</dt>
<dd>
<span class="errmsg">MAX_TRACER_FIELDS needs to be increased</span>
</dd>
<dd>
   The number of tracer fields has exceeded the maximum allowed. 
   The parameter MAX_TRACER_FIELDS needs to be increased.</dd>
<dt>
<b>FATAL in get_number_tracers</b>
</dt>
<dd>
<span class="errmsg">Model number is invalid.</span>
</dd>
<dd>
   The index of the component model is invalid.</dd>
<dt>
<b>Fatal in get_tracer_indices</b>
</dt>
<dd>
<span class="errmsg">index array size too small in get_tracer_indices</span>
</dd>
<dd>
   The global index array is too small and cannot contain all the tracer numbers.</dd>
<dt>
<b>FATAL in get_tracer_indices</b>
</dt>
<dd>
<span class="errmsg">family array size too small in get_tracer_indices</span>
</dd>
<dd>
   The family index array is too small and cannot contain all the tracer numbers.</dd>
<dt>
<b>FATAL in get_tracer_indices</b>
</dt>
<dd>
<span class="errmsg">prognostic array size too small in get_tracer_indices</span>
</dd>
<dd>
   The prognostic index array is too small and cannot contain all the tracer numbers.</dd>
<dt>
<b>FATAL in get_tracer_indices</b>
</dt>
<dd>
<span class="errmsg">diagnostic array size too small in get_tracer_indices</span>
</dd>
<dd>
   The diagnostic index array is too small and cannot contain all the tracer numbers.</dd>
<dt>
<b>NOTE in get_tracer_index</b>
</dt>
<dd>
<span class="errmsg">tracer with this name not found: X</span>
</dd>
<dd>

</dd>
<dt>
<b>FATAL in assign_tracer_field</b>
</dt>
<dd>
<span class="errmsg">invalid index</span>
</dd>
<dd>
   The index that has been passed to this routine is invalid.</dd>
<dt>
<b>FATAL in assign_tracer_field</b>
</dt>
<dd>
<span class="errmsg">At least one of data, data_tlevels or tendency must be passed in here.</span>
</dd>
<dd>
   At least one of data, data_tlevels or tendency must be passed to assign_tracer_field
   Otherwise there is not much point in calling this routine.</dd>
<dt>
<b>FATAL in get_tracer_field</b>
</dt>
<dd>
<span class="errmsg">invalid index</span>
</dd>
<dd>
   The index that has been passed to this routine is invalid.
   Check the index that is being passed corresponds to a valid
   tracer name.</dd>
<dt>
<b>FATAL in get_tracer_field</b>
</dt>
<dd>
<span class="errmsg">invalid index</span>
</dd>
<dd>
   The index that has been passed to this routine is invalid.
   Check the index that is being passed corresponds to a valid
   tracer name.</dd>
<dt>
<b>FATAL in get_tracer_field</b>
</dt>
<dd>
<span class="errmsg">tracer field array not allocated</span>
</dd>
<dd>
   The tracer array has not been allocated. This means that a
   call to assign_tracer_field is absent in the code.</dd>
<dt>
<b>FATAL in get_tracer_tlevels</b>
</dt>
<dd>
<span class="errmsg">invalid index</span>
</dd>
<dd>
   The index that has been passed to this routine is invalid.
   Check the index that is being passed corresponds to a valid
   tracer name.</dd>
<dt>
<b>FATAL in get_tracer_tlevels</b>
</dt>
<dd>
<span class="errmsg">invalid index</span>
</dd>
<dd>
   The index that has been passed to this routine is invalid.
   Check the index that is being passed corresponds to a valid
   tracer name.</dd>
<dt>
<b>FATAL in get_tracer_tlevels</b>
</dt>
<dd>
<span class="errmsg">tracer field array not allocated</span>
</dd>
<dd>
   The tracer array has not been allocated. This means that a
   call to assign_tracer_field is absent in the code.</dd>
<dt>
<b>FATAL in get_tracer_tendency</b>
</dt>
<dd>
<span class="errmsg">invalid index</span>
</dd>
<dd>
   The index that has been passed to this routine is invalid.
   Check the index that is being passed corresponds to a valid
   tracer name.</dd>
<dt>
<b>FATAL in get_tracer_tendency</b>
</dt>
<dd>
<span class="errmsg">invalid index</span>
</dd>
<dd>
   The index that has been passed to this routine is invalid.
   Check the index that is being passed corresponds to a valid
   tracer name.</dd>
<dt>
<b>FATAL in get_tracer_tendency</b>
</dt>
<dd>
<span class="errmsg">tracer tendency field array not allocated</span>
</dd>
<dd>
   The tracer array has not been allocated. This means that a
   call to assign_tracer_field is absent in the code.</dd>
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
