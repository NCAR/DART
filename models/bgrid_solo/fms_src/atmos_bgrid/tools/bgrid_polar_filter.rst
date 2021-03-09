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
<HTML>
<TITLE>module bgrid_polar_filter_mod</TITLE>
<BODY BGCOLOR="#AABBCC" TEXT="#332211" >

<DIV ALIGN="CENTER"> <FONT SIZE=-2>
<A HREF="#INTERFACE">PUBLIC INTERFACE</A> / 
<A HREF="#DATA_TYPES">DATA</A> / 
<A HREF="#ROUTINES">ROUTINES</A> / 
<A HREF="#ERRORS">ERRORS</A> / 
<A HREF="#REFERENCES">REFERENCES</A> / 
<A HREF="#NOTES">NOTES</A> 
</FONT>
<BR><BR></DIV><HR>


<H2>module bgrid_polar_filter_mod</H2>
<A NAME="HEADER">
<PRE>
     <B>Contact:</B>   B. Wyman
     <B>Reviewers:</B>
     <B>Change history:</B> <A HREF="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS
/atmos/bgrid/tools/bgrid_polar_filter.f90">WebCVS Log for bgrid_polar_filter.f90</A>

</PRE>
</A><!-- END HEADER -->
<!--------------------------------------------------------------------->
<A NAME="OVERVIEW">
<HR>
<H4>OVERVIEW</H4>
<!-- BEGIN OVERVIEW -->
<PRE>

     Provides polar filtering routines for the B-grid dynamical core.
     This is an MPP version that can be run with 1-D or 2-D decomposition.

</PRE>
</A><!-- END OVERVIEW -->
<!--------------------------------------------------------------------->
<A NAME="DESCRIPTION">
<!-- BEGIN DESCRIPTION -->
<PRE>

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
     number and latitude (see <a href="#NOTES">NOTES</a>), and then transformed back to grid point
     space using the inverse FFT.

     Load balancing for multiple processor runs is handled by gathering up full
     latitude rows on data and transmitting them equally to all processors.

</PRE>
</A><!-- END DESCRIPTION -->
<!--------------------------------------------------------------------->
<A NAME="MODULES_USED">
<HR>
<H4>OTHER MODULES USED</H4>
<!-- BEGIN MODULES_USED -->
<PRE>

     bgrid_horiz_mod
     fft_mod
     fms_mod
     mpp_mod
     mpp_domains_mod

</PRE>
</A><!-- END MODULES_USED -->
<!--------------------------------------------------------------------->
<A NAME="INTERFACE">
<HR>
<H4>PUBLIC INTERFACE</H4>
<!-- BEGIN INTERFACE -->
<PRE>

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

</PRE>
</A><!-- END INTERFACE -->
<!--------------------------------------------------------------------->
<A NAME="DATA_TYPES">
<HR>
<H4>PUBLIC DATA</H4>
<!-- BEGIN DATA_TYPES -->
<PRE>

<b>type (pfilt_control_type)</b>

    The internal contents of this type are private and cannot be accessed
    by the user.

</PRE>
</A><!-- END DATA_TYPES -->
<!--------------------------------------------------------------------->
<A NAME="ROUTINES">
<HR>
<H4>PUBLIC ROUTINES</H4>
<!-- BEGIN ROUTINES -->
<PRE>
<a name="polar_filter_init">

Control = <b>polar_filter_init</b> ( Hgrid, reflat, weight, sigma, verbose )

  INPUT

     Hgrid     horiz_grid_type (see horiz_grid_mod)

  OPTIONAL INPUT

     reflat    The reference latitude in degrees. This is the latitude
               where polar filtering begins. Equator-ward from this latitude
               there is no filtering.
                 [real, default: reflat=60.]

     weight    A weight to modify the overall strength of the filter without
               increasing the number of latitudes filtered.
                 [integer, default: weight=1]

     sigma     Logical flag the sigma or eta coordinate.
               Setting this flag equal true when running a sigma coordinate
               model will make the code more efficient.
                 [logical, default: sigma=.false.]

     verbose   verbose flag, a larger number increases the printed output
                 [integer, default: verbose=0]

  RETURNS

     Control   A variable of derived-type pfilt_control_type (see above).

-----------------------------------------------------------------
<a name="polar_filter">

<b>call polar_filter</b> (Control, data, grid, mask)
          OR
<b>call polar_filter</b> (Control, u, v, grid, mask)

  INPUT

    Control   Derived-type variable returned by a previous call
              to polar_filter_init. 
                 [type(pfilt_control_type)]

    grid      Identifier for the grid that the data in on.
              Possible values are:  grid = 1, mass grid
                                         = 2, velocity grid, u comp
                                         = 3, velocity grid, v comp
                 [integer, scalar]

  INPUT/OUTPUT

    data, u, v   The data arrays to be filtered.  The returned value
                 will be the filtered data.  The horizontal dimensions (1 and 2)
                 must be consistent with the model's <b>data domain</b>.
                 The third and fourth dimensions are arbitrary.
                 The halo regions of the returned fields have not been updated.
                 The arguments u and v do not necessarily have to be on the
                 momentum grid.
                   [real, dimension(:,:) or (:,:,:) or (:,:,:,:)]

  OPTIONAL INPUT

    mask      data mask, 0.0 (data not present) or 1.0 (data present)
              must have the same first 3 dimensions as data
                [real, dimension(:,:) or (:,:,:)]

-----------------------------------------------------------------
<a name="polar_filter_wind">

<b>call polar_filter_wind</b> ( Control, u, v, mask )

  INPUT

    Control   Derived-type variable returned by a previous call
              to polar_filter_init. 
                 [type(pfilt_control_type)]

  INPUT/OUTPUT

    u, v      velocity components to be filtered
                [real, dimension(:,:) or (:,:,:)]

  OPTIONAL INPUT

    mask      data mask, 0.0 (data not present) or 1.0 (data present)
              must have the same first 3 dimensions as data
                [real, dimension(:,:) or (:,:,:)]

</PRE>
</A><!-- END ROUTINES -->
<!--------------------------------------------------------------------->
<A NAME="ERRORS">
<HR>
<H4>ERROR MESSAGES</H4>
<!-- BEGIN ERRORS -->
<PRE>

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
        B-grid model check namelist &dynamics_driver).  Also consider
        whether round-off error may caused the problem, if so slightly
        adjust the value.

<b>FATAL error in load_balance_filter</b>

    <b>cannot get and put fft rows on the same pe</b>
       There may be a code error. Check with the delveloper.

</PRE>
</A><!-- END ERRORS -->
<!--------------------------------------------------------------------->
<A NAME="REFERENCES">
<HR>
<H4>REFERENCES</H4>
<!-- BEGIN REFERENCES -->
<PRE>

Arakawa, A. and R. Lamb, 1977: Computational design of the basic
    dynamical processes of the UCLA general circulation model. 
    Methods in Computational Physics, Vol. 17. Academic Press, 173-265.

Takacs, L. L. and R. C. Balgovind, 1983: High-latitude filtering
    in global grid-point models.  Mon. Wea. Rev., 111, 2005-2015.

Renner, V., 1981: Zonal filtering experiments with a barotropic
    model.  Contrib. Atmos. Phys., 54, 453-464.

</PRE>
</A><!-- END REFERENCES -->
<!--------------------------------------------------------------------->
<A NAME="BUGS">
<HR>
<H4>KNOWN BUGS</H4>
<!-- BEGIN BUGS -->
<PRE>

     None.

</PRE>
</A><!-- END BUGS -->
<!--------------------------------------------------------------------->
<A NAME="NOTES">
<HR>
<H4>NOTES</H4>
<!-- BEGIN NOTES -->
<PRE>

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

</PRE>
</A><!-- END NOTES -->
<!--------------------------------------------------------------------->
<A NAME="PLANS">
<HR>
<H4>FUTURE PLANS</H4>
<!-- BEGIN PLANS -->
<PRE>

     None.

</PRE>
</A><!-- END PLANS -->
<!--------------------------------------------------------------------->

<HR>
</BODY>
</HTML>
