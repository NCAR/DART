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
<TITLE>module atmosphere_mod</TITLE>
<BODY BGCOLOR="#AABBCC" TEXT="#332211" >

<DIV ALIGN="CENTER"> <FONT SIZE=1>
<A HREF="#INTERFACE">PUBLIC INTERFACE</A> / 
<A HREF="#ROUTINES">ROUTINES</A> / 
<A HREF="#NAMELIST">NAMELIST</A> / 
<A HREF="#ERRORS">ERRORS</A>
</FONT>
<BR><BR></DIV><HR>

<H2>module atmosphere_mod</H2>
<A NAME="HEADER">
<PRE>
<PRE>
     <B>Contact:</B>   Bruce Wyman
     <B>Reviewers:</B>
     <B>Change history:</B> <A HREF="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS
/atmos_bgrid/driver/solo/atmosphere.f90">WebCVS Log for solo/atmosphere.f90</A>

</PRE>
</A><!-- END HEADER -->
<!-- ------------------------------------------------------------------>
<A NAME="OVERVIEW">
<HR>
<H4>OVERVIEW</H4>
<!-- BEGIN OVERVIEW -->
<PRE>

     Atmospheric driver for the dry B-grid dynamical core and Held-Suarez
     benchmark (aka simple_phyhsics).
</PRE>
</A><!-- END OVERVIEW -->
<!-- ------------------------------------------------------------------>
<A NAME="DESCRIPTION">
<!-- BEGIN DESCRIPTION -->
<PRE>
     This module provides a standard interface to the B-grid dynamical
     core and B-grid interface to the simple physics.
     A similar interface for the spectral dynamical core exists that
     may be easily switched with this interface.

</PRE>
</A><!-- END DESCRIPTION -->
<!-- ------------------------------------------------------------------>
<A NAME="OTHER MODULES USED">
<HR>
<H4>OTHER MODULES USED</H4>
<!-- BEGIN OTHER MODULES USED -->
<PRE>

   bgrid_horiz_mod
   bgrid_vert_mod
   bgrid_prog_var_mod
   bgrid_halo_mod
   bgrid_grid_change_mod
   bgrid_core_driver_mod
   hs_forcing_mod
   time_manager_mod
   fms_mod

</PRE>
</A><!-- END MODULES_USED -->
<!-- ------------------------------------------------------------------>
<A NAME="PUBLIC INTERFACE">
<HR>
<H4>PUBLIC INTERFACE</H4>
<!-- BEGIN PUBLIC INTERFACE -->
<PRE>

  <b>use atmosphere_mod</b> [,only: atmosphere_init,       atmosphere_end,
                             atmosphere,
                             atmosphere_resolution, atmosphere_boundary,
                             get_atmosphere_axes  ]

  NOTES

     1)  Optional namelist interface <b>&atmosphere_nml</b> may be
         read from file <b>input.nml</b>.
                                
</PRE>
</A><!-- END PUBLIC INTERFACE -->
<!-- ------------------------------------------------------------------>
<A NAME="PUBLIC ROUTINES">
<HR>
<H4>PUBLIC ROUTINES</H4>
<!-- BEGIN PUBLIC ROUTINES -->
<PRE>

<b>call atmosphere_init</b> ( Time_init, Time, Time_step )

DESCRIPTION
   Initialization call for running the bgrid dynamical with the
   Held-Suarez GCM forcing.

INPUT
   Time_init   The initial (or base) time.  <FONT SIZE=-1 COLOR="#000099">[time_type]</FONT>

   Time        The current time.  <FONT SIZE=-1 COLOR="#000099">[time_type]</FONT>

   Time_step   The atmospheric model/physics time step.  <FONT SIZE=-1 COLOR="#000099">[time_type]</FONT>

<!------------------------>
<HR WIDTH="50%" ALIGN=LEFT>

<b>call atmosphere_end</b>

DESCRIPTION
   Termination call for the bgrid dynamical with Held-Suarez
   GCM forcing.  There are no arguments to this routine.

<!------------------------>
<HR WIDTH="50%" ALIGN=LEFT>

<b>call atmosphere</b> ( Time )

DESCRIPTION
   Advances the B-grid prognostic variables one time step forward.
   The dynamical core, Held-Suarez forcing, diagnostics, and time
   differencing are all called.  This routine should only be called
   once per time step.

INPUT
   Time    The current time.  <FONT SIZE=-1 COLOR="#000099">[time_type]</FONT>

NOTE
   The prognostic variables are stored in a private derived-type 
   variabledefined in the header section of this module.

---------------------------------------------------------------------

<b>call get_atmosphere_axes</b> ( axes )

OUTPUT

   axes        axis identifiers for the atmospheric grids
               The size of axes at least 1 but not greater than 4.
               The axes returned are ordered (/ x, y, p_full, p_half /).
                 [integer, dimension(:)]

<!------------------------>
<HR WIDTH="50%" ALIGN=LEFT>

<b>call atmosphere_resolution</b> ( nlon, nlat <FONT COLOR="#007700">[, global]</FONT> )

DESCRIPTION
   Returns the resolution of compute domain for either the
   current processor or the global domain.

OUTPUT
   nlon   The number of longitude points in the compute domain.
             <FONT SIZE=-1 COLOR="#000099">[integer]</FONT>

   nlat   The number of latitude points in the compute domain.
             <FONT SIZE=-1 COLOR="#000099">[integer]</FONT>

OPTIONAL INPUT

   global  Flag that specifies whether the returned compute domain size is
           for the global grid (TRUE) or for the current processor (FALSE).
              <FONT SIZE=-1 COLOR="#000099">[logical, default: FALSE]</FONT>
           
<!------------------------>
<HR WIDTH="50%" ALIGN=LEFT>

<b>call atmosphere_boundary</b> ( blon, blat <FONT COLOR="#007700">[, global]</FONT> )

DESCRIPTION
   Returns the grid box edges of compute domain for either the
   current processor or the global domain.

OUTPUT
   blon    The west-to-east longitude edges of grid boxes (in radians).
              <FONT SIZE=-1 COLOR="#000099">[real, dimension(nlon+1)]</FONT>

   blat    The south-to-north latitude edges of grid boxes (in radians).
              <FONT SIZE=-1 COLOR="#000099">[real, dimension(nlat+1)]</FONT>

OPTIONAL INPUT
   global  Flag that specifies whether the returned grid box edges are
           for the global grid (TRUE) or for the current processor (FALSE).
              <FONT SIZE=-1 COLOR="#000099">[logical, default: FALSE]</FONT>
           
NOTE
   The size of the output arguments, blon and blat, must be +1 more than the
   output arguments for call atmosphere_resolution, nlon+1 and nlat+1, respectively.

</PRE>
</A><!-- END PUBLIC ROUTINES -->
<!-- ------------------------------------------------------------------>
<A NAME="NAMELIST">
<HR>
<H4>NAMELIST</H4>
<!-- BEGIN NAMELIST -->
<PRE>

<b>&atmosphere_nml</b>

 physics_window  The number of "i" and "j" rows processed each time
                 the modular physics is called. To process the entire
                 domain use physics_window = 0,0.
                    <FONT SIZE=-1 COLOR="#000099">[integer, default: physics_window = 0,0]</FONT>

</PRE>
</A><!-- END NAMELIST -->
<!-- ------------------------------------------------------------------>
<A NAME="ERRORS">
<HR>
<H4>ERROR MESSAGES</H4>
<!-- BEGIN ERRORS -->
<PRE>

<b>FATAL errors from get_atmosphere_axes in atmosphere_mod</b>

    <b>size of argument is incorrect</b>
        The size of the argument to get_atmosphere_axes must be
        between 1 and 4.


</PRE>
</A><!-- END ERRORS -->
<!-- ------------------------------------------------------------------>

<HR>
</BODY>
</HTML>
