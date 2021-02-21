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
<!DOCTYPE html>
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>module bgrid_diagnostics_mod</title>
</head>
<body bgcolor="#AABBCC" text="#332211">
<div align="center"><font size="-2"><a href="#INTERFACE">PUBLIC
INTERFACE</a> / <a href="#ROUTINES">ROUTINES</a> / <a href=
"#DIAGNOSTICS">DIAGNOSTICS</a> / <a href="#ERRORS">ERRORS</a> /
<a href="#NOTES">NOTES</a></font><br>
<br></div>
<hr>
<h1>module bgrid_diagnostics_mod</h1>
<a name="OVERVIEW" id="OVERVIEW">
<hr>
<h2>OVERVIEW</h2>
<!-- BEGIN OVERVIEW -->
<pre>

     Handles the archiving of B-grid prognostic variables and
     several other related diagnostics fields to NetCDF files.

</pre></a><!-- END OVERVIEW -->
<!--------------------------------------------------------------------->
 <a name="DESCRIPTION" id="DESCRIPTION"><!-- BEGIN DESCRIPTION -->
<pre>
     This version provides an interface to the general model diagnostics
     manager (diag_manager) which provides support for NetCDF and 
     distributed memory systems.  This also allows dynamical variables
     to be saved with variables from other parts of the (atmospheric) model.

     In this version the user control of file names, output interval,
     time averaging, packing, and many other features is handled through
     the diagnostics manager's input file called diag_table.

     One of the functions of this module is to initialize the atmospheric
     model's diagnostic axes. The axis identification returned by this
     module is then used throughout the atmospheric model wherever
     diagnostic fields are saved (e.g., the physics).

     A second function of this module is to initialize and save 
     the dynamical core's prognostic fields and time tendencies of the
     prognostic fields.

     The fields saved are:

        variables:   bk  vertical coordinate sigma/eta values
                     pk  vertical coordinate pressure values (pascals)
                  zsurf  height of the surface (m)
                    res  reciprocal of eta at the surface
                     ps  surface pressure (pascals)
                
                  ucomp  zonal wind component (m/sec)
                  vcomp  meridional wind component (m/sec)
                   temp  temperature (deg_k)
                  omega  omega vertical velocity (pascals/sec)
                tracers  multiple number of tracers fields

                   wspd  wind speed (m/s)
                    div  divergence (1/s)
                    vor  relative vorticity (1/s)

                udt_dyn  zonal wind tendency for dynamics (m/s2)
                vdt_dyn  meridional wind tendency for dynamics (m/s2)
                tdt_dyn  temperature tendency for dynamics (deg_k/sec)
        tracers"_dt_dyn" tracers tendencies for dynamics

     NOTES

     1)  The static variables bk, pk, zsurf, and res do not vary
         in time, they are written once and have no time axis.

</pre></a><!-- END DESCRIPTION -->
<!--------------------------------------------------------------------->
 <a name="MODULES_USED" id="MODULES_USED">
<hr>
<h2>OTHER MODULES USED</h2>
<!-- BEGIN MODULES_USED -->
<pre>

   bgrid_horiz_mod
   bgrid_vert_mod
   bgrid_masks_mod
   bgrid_prog_var_mod
   bgrid_change_grid_mod
   diag_manager_mod
   time_manager_mod
   fms_mod
   constants_mod
   field_manager_mod
   tracer_manager_mod

</pre></a><!-- END MODULES_USED -->
<!--------------------------------------------------------------------->
 <a name="INTERFACE" id="INTERFACE">
<hr>
<h2>PUBLIC INTERFACE</h2>
<!-- BEGIN INTERFACE -->
<pre>

<b>use bgrid_diagnostics_mod</b> [, only: bgrid_diagnostics_init,
                                 bgrid_diagnostics,
                                 bgrid_diagnostics_tend  ]

    <a href="#bgrid_diagnostics_init">bgrid_diagnostics_init</a>
         Must be called to initialize the module.
         Also, axis identifiers (for the diag_manager) are returned.
         The static fields, zsurf and res are saved to the output file.
 
    <a href="#bgrid_diagnostics">bgrid_diagnostics</a>
         Should be called at the end of every time step.
         Outputs the prognostics variables, omega, wind speed, divergence,
         and vorticity.  The time passed to bgrid_diagnostics will used to
         determine if the output fields should be written.
 
    <a href="#bgrid_diagnostics_tend">bgrid_diagnostics_tend</a>
         Should be called every time step immediately after dynamics..
         Outputs the tendencies of the prognostics variables due to the
         dynamics. 
 
</pre></a><!-- END INTERFACE -->
<!--------------------------------------------------------------------->
 <a name="ROUTINES" id="ROUTINES">
<hr>
<h2>PUBLIC ROUTINES</h2>
<!-- BEGIN ROUTINES -->
<pre>
<a name="bgrid_diagnostics_init" id="bgrid_diagnostics_init">

<b>call bgrid_diagnostics_init</b> ( Time, Hgrid, Vgrid, Var, fis, res, mass_axes, vel_axes ) INPUT Time The current time. [time_type] Hgrid Derived-type variable containing horizontal grid constants for the B-grid dynamical core. [horiz_grid_type] Vgrid Derived-type variable containing vertical grid constants for the B-grid dynamical core. [vert_grid_type] Var Derived-type variable containing the model's prognostics variables. [prog_var_type] fis Surface geopotential height, a static field that is passed to the diagnostics manager and saved to the NetCDF file as a static field (i.e., no time axis). [real, dimension(Hgrid%ilb:,Hgrid%jlb:)] res Reciprocal of eta at the surface, saved to the NetCDF file as a static field (i.e., no time axis). [real, dimension(Hgrid%ilb:,Hgrid%jlb:)] OUTPUT mass_axes Integer identifiers for the mass axes (x,y,pfull,phalf). [integer, dimension(4)] vel_axes Integer identifiers for the velocity axes (x,y,pfull,phalf). Note that the vertical axes (pfull,phalf) will be the same as for the mass axes. [integer, dimension(4)] NOTES 1) Arguments fis and res must have dimensions consistent with the dynamical core (in Hgrid), the storage can be allocated using module prog_var_mod. 2) Static fields fis and res are written to the output (NetCDF) file if requested by the diagnostics manager. ------------------------------------------------------------------ </a><a name="bgrid_diagnostics"
id="bgrid_diagnostics">

<b>call bgrid_diagnostics</b> ( Hgrid, Vgrid, Var, Masks, Time, omega ) INPUT Hgrid Derived-type variable containing horizontal grid constants for the B-grid dynamical core. [horiz_grid_type] Vgrid Derived-type variable containing vertical grid constants for the B-grid dynamical core. [vert_grid_type] Var Derived-type variable containing the model's prognostics variables. [prog_var_type] Masks Derived-type variable containing grid mask information. [grid_mask_type] Time The current time. [time_type] omega Omega diagnostic (in Pascals/sec). [real, dimension(Hgrid%ilb:,Hgrid%jlb:,:)] ------------------------------------------------------------------ </a><a name="bgrid_diagnostics_tend"
id="bgrid_diagnostics_tend">

<b>call bgrid_diagnostics_tend</b> ( Hgrid, Var_dt, Masks, Time ) INPUT Hgrid Derived-type variable containing horizontal grid constants for the B-grid dynamical core. [horiz_grid_type] Var_dt Derived-type variable containing the tendencies of the model's prognostics variables. [prog_var_type] Masks Derived-type variable containing grid mask information. [grid_mask_type] Time The current time. [time_type] </a></pre></a><!-- END ROUTINES -->
<!--------------------------------------------------------------------->
 <a name="DIAGNOSTICS" id="DIAGNOSTICS">
<hr>
<h2>DIAGNOSTIC FIELDS</h2>
<pre>
Diagnostic fields may be output to a netcdf file by specifying the
module name <b>dynamics</b> and the desired field names (given below)
in file <b>diag_table</b>. See the documentation for diag_manager.
</pre>
<!-- BEGIN DIAGNOSTICS -->
<pre>

Diagnostic fields for module name: <b>dynamics</b>

   field name   field description
   ----------   -----------------

    bk          vertical coordinate eta value
    pk          vertical coordinate pressure value (Pascals)
    zsurf       height of the surface (m)
    res         reciprocal of eta at the surface

    ps          surface pressure (Pascals)
    ucomp       zonal wind component (m/s)
    vcomp       meridional wind component (m/s)
    temp        temperature (deg_k)

    omega       omega vertical velocity (pascals/sec)
    wspd        wind speed (m/s)
    div         divergence (1/s)
    vor         relative vorticity (1/s)

    ucomp_sq    zonal wind component squared (m2/s2)
    vcomp_sq    meridional wind component squared (m2/s2)
    temp_sq     temperature squared (deg_K**2)
    omega_sq    omega vertical velocity squared (Pa**2/s**2)
    ucomp_vcomp zonal times meridional wind components (m2/s2)
    omega_temp  omega vertical velocity times temperature (Pascals*deg_K/sec)

    udt_dyn     zonal wind tendency for dynamics (m/s2)
    vdt_dyn     meridional wind tendency for dynamics (m/s2)
    tdt_dyn     temperature tendency for dynamics (deg_k/sec)

    trname         tracers concentrations
    trname_dt_dyn  tracers tendencies for dynamics

Notes:

   1) bk, pk, zsurf, res are static variables

   2) In this version, div, vor are recomputed for diagnostic

   3) The names for tracer tendencies are determined from the tracer names.
      The short names will be appended with "_dt_dyn", the long names are
      modified to indicated a tendency for dynamics, and the units are
      appended with "/s".

   4) All available diagnostic tracer fields will be printed to the logfile.

</pre></a><!-- END DIAGNOSTICS -->
<!--------------------------------------------------------------------->
 <a name="ERRORS" id="ERRORS">
<hr>
<h2>ERROR MESSAGES</h2>
<!-- BEGIN ERRORS -->
<pre>

   None.

</pre></a><!-- END ERRORS -->
<!--------------------------------------------------------------------->
 <a name="BUGS" id="BUGS">
<hr>
<h2>KNOWN BUGS</h2>
<!-- BEGIN BUGS -->
<pre>

     None.

</pre></a><!-- END BUGS -->
<!--------------------------------------------------------------------->
 <a name="NOTES" id="NOTES">
<hr>
<h2>NOTES</h2>
<!-- BEGIN NOTES -->
<pre>

   None.

</pre></a><!-- END NOTES -->
<!--------------------------------------------------------------------->
 <a name="PLANS" id="PLANS">
<hr>
<h2>FUTURE PLANS</h2>
<!-- BEGIN PLANS -->
<pre>

   None.

</pre></a><!-- END PLANS -->
<!--------------------------------------------------------------------->
<hr>
</body>
</html>
