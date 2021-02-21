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
<title>module bgrid_polar_filter_mod</title>
</head>
<body bgcolor="#AABBCC" text="#332211">
<div align="center"><font size="-2"><a href="#INTERFACE">PUBLIC
INTERFACE</a> / <a href="#DATA_TYPES">DATA</a> / <a href=
"#ROUTINES">ROUTINES</a> / <a href="#ERRORS">ERRORS</a> / <a href=
"#REFERENCES">REFERENCES</a> / <a href=
"#NOTES">NOTES</a></font><br>
<br></div>
<hr>
<h1>module bgrid_polar_filter_mod</h1>
<a name="OVERVIEW" id="OVERVIEW">
<hr>
<h2>OVERVIEW</h2>
<!-- BEGIN OVERVIEW -->
<pre>

     Provides polar filtering routines for the B-grid dynamical core.
     This is an MPP version that can be run with 1-D or 2-D decomposition.

</pre></a><!-- END OVERVIEW -->
<!--------------------------------------------------------------------->
 <a name="DESCRIPTION" id="DESCRIPTION"><!-- BEGIN DESCRIPTION -->
<pre>

     The polar filtering scheme is used at high latitudes to damp the shortest
     resolvable waves so that a longer time step can be taken.  Filtering is
     applied to the mass divergence, horizontal omega-alpha tendency,
     horizontal advective tendency of temperature and prognostic tracers,
     and the momentum components.  The momentum components are transformed to
     stereographic coordinates before they are filtered to minimize distortion
     near the poles.  The filtering scheme conserves mass and tracer mass,
     but does not conserve energy (Takacs and Balgovind, 1983).

     The fields are filtered by transforming a full latitude circle of data
     to Fourier components using a fast Fourier transform (FFT). The Fourier
     components are damped (i.e., multiplied) by a given function of wave
     number and latitude (see <a href=
"#NOTES">NOTES</a>), and then transformed back to grid point
     space using the inverse FFT.

     Load balancing for multiple processor runs is handled by gathering up full
     latitude rows on data and transmitting them equally to all processors.

</pre></a><!-- END DESCRIPTION -->
<!--------------------------------------------------------------------->
 <a name="MODULES_USED" id="MODULES_USED">
<hr>
<h2>OTHER MODULES USED</h2>
<!-- BEGIN MODULES_USED -->
<pre>

     bgrid_horiz_mod
     fft_mod
     fms_mod
     mpp_mod
     mpp_domains_mod

</pre></a><!-- END MODULES_USED -->
<!--------------------------------------------------------------------->
 <a name="INTERFACE" id="INTERFACE">
<hr>
<h2>PUBLIC INTERFACE</h2>
<!-- BEGIN INTERFACE -->
<pre>

  <b>use bgrid_polar_filter_mod</b> [,only: pfilt_control_type,
                                     polar_filter_init ,
                                     polar_filter      ,
                                     polar_filter_wind ]

  <a href="#DATA_TYPES">pfilt_control_type</a>
       Derived-type variable containing constants needed by the polar filtering
       interfaces. It is returned by the function polar_filter_init and is
       used by all other routine in this module.

  <a href="#polar_filter_init">polar_filter_init</a>
       Function that initializes a variable of derived-type pfilt_control_type.

  <a href="#polar_filter">polar_filter</a>
       Performs polar filtering of either one or two B-grid fields on either
       the mass or momentum grid.

  <a href="#polar_filter_wind">polar_filter_wind</a>
       Performs polar filtering of the velocity components. The components
       are transformed to stereographic components for the actually filtering.

</pre></a><!-- END INTERFACE -->
<!--------------------------------------------------------------------->
 <a name="DATA_TYPES" id="DATA_TYPES">
<hr>
<h2>PUBLIC DATA</h2>
<!-- BEGIN DATA_TYPES -->
<pre>

<b>type (pfilt_control_type)</b>

    The internal contents of this type are private and cannot be accessed
    by the user.

</pre></a><!-- END DATA_TYPES -->
<!--------------------------------------------------------------------->
 <a name="ROUTINES" id="ROUTINES">
<hr>
<h2>PUBLIC ROUTINES</h2>
<!-- BEGIN ROUTINES -->
<pre>
<a name="polar_filter_init" id="polar_filter_init">

Control = <b>polar_filter_init</b> ( Hgrid, reflat, weight, sigma, verbose ) INPUT Hgrid horiz_grid_type (see horiz_grid_mod) OPTIONAL INPUT reflat The reference latitude in degrees. This is the latitude where polar filtering begins. Equator-ward from this latitude there is no filtering. [real, default: reflat=60.] weight A weight to modify the overall strength of the filter without increasing the number of latitudes filtered. [integer, default: weight=1] sigma Logical flag the sigma or eta coordinate. Setting this flag equal true when running a sigma coordinate model will make the code more efficient. [logical, default: sigma=.false.] verbose verbose flag, a larger number increases the printed output [integer, default: verbose=0] RETURNS Control A variable of derived-type pfilt_control_type (see above). ----------------------------------------------------------------- </a><a name="polar_filter"
id="polar_filter">

<b>call polar_filter</b> (Control, data, grid, mask) OR <b>call polar_filter</b> (Control, u, v, grid, mask) INPUT Control Derived-type variable returned by a previous call to polar_filter_init. [type(pfilt_control_type)] grid Identifier for the grid that the data in on. Possible values are: grid = 1, mass grid = 2, velocity grid, u comp = 3, velocity grid, v comp [integer, scalar] INPUT/OUTPUT data, u, v The data arrays to be filtered. The returned value will be the filtered data. The horizontal dimensions (1 and 2) must be consistent with the model's <b>data domain</b>. The third and fourth dimensions are arbitrary. The halo regions of the returned fields have not been updated. The arguments u and v do not necessarily have to be on the momentum grid. [real, dimension(:,:) or (:,:,:) or (:,:,:,:)] OPTIONAL INPUT mask data mask, 0.0 (data not present) or 1.0 (data present) must have the same first 3 dimensions as data [real, dimension(:,:) or (:,:,:)] ----------------------------------------------------------------- </a><a name="polar_filter_wind"
id="polar_filter_wind">

<b>call polar_filter_wind</b> ( Control, u, v, mask ) INPUT Control Derived-type variable returned by a previous call to polar_filter_init. [type(pfilt_control_type)] INPUT/OUTPUT u, v velocity components to be filtered [real, dimension(:,:) or (:,:,:)] OPTIONAL INPUT mask data mask, 0.0 (data not present) or 1.0 (data present) must have the same first 3 dimensions as data [real, dimension(:,:) or (:,:,:)] </a></pre></a><!-- END ROUTINES -->
<!--------------------------------------------------------------------->
 <a name="ERRORS" id="ERRORS">
<hr>
<h2>ERROR MESSAGES</h2>
<!-- BEGIN ERRORS -->
<pre>

<b>FATAL error in polar_filter or polar_filter_wind</b>

    <b>input array has the wrong dimensions</b>
        The horizontal dimensions (1 and 2) have the wrong size.
        The dimensions must include any halo row points in the 
        model's data domain grid.

<b>FATAL error in polar_filter_......</b>

    <b>invalid grid arg</b>
        The argument "grid" must have the value 1,2, or 3.
        Look at the documentation of these routines for details.

<b>FATAL error in polar_filter_init</b>

    <b>reflat must lie between 0 and 90.</b>
        The reflat (reference latitude in degrees) was not in the
        required range.  Check the namelist value you supplied (for the
        B-grid model check namelist &amp;dynamics_driver).  Also consider
        whether round-off error may caused the problem, if so slightly
        adjust the value.

<b>FATAL error in load_balance_filter</b>

    <b>cannot get and put fft rows on the same pe</b>
       There may be a code error. Check with the delveloper.

</pre></a><!-- END ERRORS -->
<!--------------------------------------------------------------------->
 <a name="REFERENCES" id="REFERENCES">
<hr>
<h2>REFERENCES</h2>
<!-- BEGIN REFERENCES -->
<pre>

Arakawa, A. and R. Lamb, 1977: Computational design of the basic
    dynamical processes of the UCLA general circulation model. 
    Methods in Computational Physics, Vol. 17. Academic Press, 173-265.

Takacs, L. L. and R. C. Balgovind, 1983: High-latitude filtering
    in global grid-point models.  Mon. Wea. Rev., 111, 2005-2015.

Renner, V., 1981: Zonal filtering experiments with a barotropic
    model.  Contrib. Atmos. Phys., 54, 453-464.

</pre></a><!-- END REFERENCES -->
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

 Filter Response Function
 ------------------------

     Damping factors as a function of wave number (k) and latitude index (j)
     are defined as:                                        
                                                            
        SS(k,j) = { cos_PH(j) / [ cos_PH(ref) * sin_X ] } ** weight

            where  X = k * dx / 2.
                  dx = longitudinal grid spacing in radians

                  cos_PH(j)   = cosine of latitude at row j
                  cos_PH(ref) = cosine of reference latitude
                  weight      = optional parameter for increasing the overall
                                strength of the filter (default = 1)

     The reference latitude is the first non-filtering latitude.


 MPP Implementation
 ------------------

 Approach for one-dimensional decomposition in Y

    An algorithm was devised that transmits latitude rows of data from PEs
    that too many rows to PEs that have too few rows. This algorithm is applied
    in reverse after the filtering has been done.

 Approach for two-dimensional decomposition in X and Y

    An additional step is applied to the one-dimensional algorithm where
    the first PE in a latitude row gathers the data (using mpp_transmit)
    from the other PEs at that latitude.  At this point the problem is
    identical to the one-dimensional case.

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
