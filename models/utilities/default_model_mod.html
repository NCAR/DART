<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN"
          "http://www.w3.org/TR/html4/strict.dtd">
<HTML>
<HEAD>
<TITLE>module model_mod</TITLE>
<link rel="stylesheet" type="text/css" href="../../docs/html/doc.css" />
<link href="../../docs/images/dart.ico" rel="shortcut icon" />
</HEAD>
<BODY>
<A NAME="TOP"></A>

<H1>MODULE model_mod</H1>

<table border=0 summary="" cellpadding=5>
<tr>
    <td valign=middle>
    <img src="../../docs/images/Dartboard7.png" alt="DART project logo" height=70 />
    </td>
    <td>Jump to <a href="../../docs/index.html">DART Documentation Main Index</a></td>
</tr>
</table>

<A HREF="#Namelist">NAMELIST</A> /
<A HREF="#Interface">INTERFACES</A> /
<A HREF="#FilesUsed">FILES</A> /
<A HREF="#References">REFERENCES</A> /
<A HREF="#Errors">ERRORS</A> /
<A HREF="#FuturePlans">PLANS</A> /
<A HREF="#PrivateComponents">PRIVATE COMPONENTS</A> /
<A HREF="#Legalese">TERMS OF USE</A>

<H2>Overview</H2>

<P>
Every model that is DART compliant must provide an set of interfaces
that will be called by DART code.  For models which have no special
code for some of these routines, they can pass through the call to this
default module, which satisfies the call but does no work.
To use these routines in a <em class=file>model_mod.f90</em>,
add at the top:<br>
<pre>
use default_model_mod, only : xxx, yyy
</pre>
and then leave them in the public list.
</P>

<!--=====================================================================-->
<!--===================== DESCRIPTION OF A NAMELIST =====================-->
<!--=====================================================================-->

<A NAME="Namelist"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>NAMELIST</H2>
<P>
The default routines have no namelist.
</P>

<!--==================================================================-->
<P><!-- useless spacer to make the top behave correctly --></P>

<A NAME="Interface"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>OTHER MODULES USED</H2>
<PRE>
types_mod
time_manager_mod
location_mod
utilities_mod
netcdf_utilities_mod
ensemble_manager_mod
dart_time_io_mod
</PRE>

<!--==================================================================-->
<!-- Declare all public entities ...                                  -->
<!-- duplicate public routines template as many times as necessary    -->
<!-- make sure you replace all yyyroutine?? strings                   -->
<!--==================================================================-->
<!--Note to authors. The first row of the table is different.         -->
<!--==================================================================-->

<div class="top">[<a href="#">top</a>]</div><hr />
<H2>PUBLIC INTERFACES</H2>

<TABLE>
<TR><TD><em class=call>use model_mod, only : </em></TD>
                   <TD><A HREF="#get_model_size">get_model_size</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#adv_1step">adv_1step</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#get_state_meta_data">get_state_meta_data</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#model_interpolate">model_interpolate</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#shortest_time_between_assimilations">shortest_time_between_assimilations</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#static_init_model">static_init_model</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#init_time">init_time</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#fail_init_time">fail_init_time</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#init_conditions">init_conditions</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#fail_init_conditions">fail_init_conditions</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#nc_write_model_atts">nc_write_model_atts</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#nc_write_model_vars">nc_write_model_vars</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#pert_model_copies">pert_model_copies</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#get_close_obs">get_close_obs</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#get_close_state">get_close_state</A></TD></TR>

<TR><TD>&nbsp;</TD><TD><A HREF="#convert_vertical_obs">convert_vertical_obs</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#convert_vertical_state">convert_vertical_state</A></TD></TR>

<TR><TD>&nbsp;</TD><TD><A HREF="#read_model_time">read_model_time</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#write_model_time">write_model_time</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#end_model">end_model</A></TD></TR>
</TABLE>

<P>
   A note about documentation style.
   Optional arguments are enclosed in brackets
   <em class=optionalcode>[like this]</em>.
</P>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="get_model_size"></A>
<br>
<div class=routine>
<em class=call>model_size = get_model_size( )</em>
<pre>
integer(i8) :: <em class=code>get_model_size</em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
Returns the length of the model state vector as 1.
Probably not what you want.  The model_mod should
set this to the right size and not use this routine.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>model_size</em></TD>
    <TD>The length of the model state vector.</TD></TR>

</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="adv_1step"></A>
<br>
<div class=routine>
<em class=call>call adv_1step(x, time)</em>
<pre>
real(r8), dimension(:), intent(inout) :: <em class=code>x</em>
type(time_type),        intent(in)    :: <em class=code>time</em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
Throws a fatal error.  If the model_mod can advance the model it
should provide a real routine.  This default routine is intended
for use by models which cannot advance themselves from inside filter.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>x</em></TD>
    <TD>State vector of length model_size.</TD></TR>

<TR><TD valign=top><em class=code>time</em></TD>
    <TD>Current time of model state.</TD></TR>

</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="get_state_meta_data"></A>
<br>
<div class=routine>
<em class=call>call get_state_meta_data (index_in, location, 
                          <em class=optionalcode>[,&nbsp;var_type]</em> )</em>
<pre>
integer,             intent(in)  :: <em class=code>index_in</em>
type(location_type), intent(out) :: <em class=code>location</em>
integer, optional,   intent(out) :: <em class=optionalcode> var_type </em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
Sets the location to missing and the variable type to 0.
The model_mod should provide a routine that sets a real location
and a state vector type for the requested item in the state vector.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>index_in</em></TD>
    <TD>Index of state vector element about which information is requested.</TD></TR>

<TR><TD valign=top><em class=code>location</em></TD>
    <TD>The location of state variable element.</TD></TR>

<TR><TD valign=top><em class=optionalcode>var_type</em></TD>
    <TD>The generic quantity of the state variable element.</TD></TR>

</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="model_interpolate"></A>
<br>
<div class=routine>
<em class=call>call model_interpolate(state_handle, ens_size, location, obs_quantity, expected_obs, istatus)</em>
<pre>
type(ensemble_type),    intent(in)  :: <em class=code>state_handle</em>
integer,                intent(in)  :: <em class=code>ens_size</em>
type(location_type),    intent(in)  :: <em class=code>location</em>
integer,                intent(in)  :: <em class=code>obs_quantity</em>
real(r8),               intent(out) :: <em class=code>expected_obs(ens_size)</em>
integer,                intent(out) :: <em class=code>istatus(ens_size)</em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
Sets the expected obs to missing and returns an error code
for all obs.  This routine should be supplied by the model_mod.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>state_handle</em></TD>
    <TD>The handle to the state structure containing information about 
        the state vector about which information is requested.</TD></TR>

<TR><TD valign=top><em class=code>ens_size</em></TD>
    <TD>The ensemble size.</TD></TR>

<TR><TD valign=top><em class=code>location</em></TD>
    <TD>Location to which to interpolate.</TD></TR>

<TR><TD valign=top><em class=code>obs_quantity</em></TD>
    <TD>Quantity of state field to be interpolated.</TD></TR>

<TR><TD valign=top><em class=code>expected_obs</em></TD>
    <TD>The interpolated values from the model.</TD></TR>

<TR><TD valign=top><em class=code>istatus</em></TD>
    <TD>Integer values return 0 for success.
        Other positive values can be defined for various failures.</TD></TR>

</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="shortest_time_between_assimilations"></A>
<br>
<div class=routine>
<em class=call>var = shortest_time_between_assimilations()</em>
<pre>
type(time_type) :: <em class=code>shortest_time_between_assimilations</em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
Returns 1 day.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>var</em></TD>
    <TD>Smallest advance time of the model.</TD></TR>

</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="static_init_model"></A>
<br>
<div class=routine>
<em class=call>call static_init_model()</em>
</div>

<div class=indent1>
<!-- Description -->

<P>
Does nothing.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>
</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="init_time"></A>
<br>
<div class=routine>
<em class=call>call init_time(time)</em>
<pre>
type(time_type), intent(out) :: <em class=code>time</em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
Returns a time of 0.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>time</em></TD>
    <TD>Initial model time.</TD></TR>

</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="fail_init_time"></A>
<br>
<div class=routine>
<em class=call>call fail_init_time(time)</em>
<pre>
type(time_type), intent(out) :: <em class=code>time</em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
Throws a fatal error.
This is appropriate for models that cannot start from arbitrary initial conditions.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>time</em></TD>
    <TD>NOT SET. Initial model time.</TD></TR>

</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="init_conditions"></A>
<br>
<div class=routine>
<em class=call>call init_conditions(x)</em>
<pre>
real(r8), dimension(:), intent(out) :: <em class=code>x</em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
Returns x(:) = 0.0
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>x</em></TD>
    <TD>Initial conditions for state vector.</TD></TR>

</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="fail_init_conditions"></A>
<br>
<div class=routine>
<em class=call>call fail_init_conditions(x)</em>
<pre>
real(r8), dimension(:), intent(out) :: <em class=code>x</em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
Throws a fatal error.
This is appropriate for models that cannot start from arbitrary initial conditions.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>x</em></TD>
    <TD>NOT SET: Initial conditions for state vector.</TD></TR>

</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="nc_write_model_atts"></A>
<br>
<div class=routine>
<em class=call>call nc_write_model_atts(ncFileID, domain_id)</em>
<pre>
integer, intent(in) :: <em class=code>ncFileID</em>
integer, intent(in) :: <em class=code>domain_id</em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
Does nothing.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>ncFileID</em></TD>
    <TD>Integer file descriptor to previously-opened netCDF file.</TD></TR>

<TR><TD valign=top><em class=code>domain_id</em></TD>
        <TD>integer describing the domain (which can be a nesting level, a component 
        model ...) Models with nested grids are decomposed into 'domains' in DART.
        The concept is extended to refer to 'coupled' models where one model component
        may be the atmosphere, another component may be the ocean, or land, or 
        ionosphere ... these would be referenced as different domains. </TD></TR>

</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="nc_write_model_vars"></A>
<br>
<div class=routine>
<em class=call>call nc_write_model_vars(ncFileID, domain_id, state_ens_handle 
  <em class=optionalcode>[,&nbsp;memberindex]</em>
  <em class=optionalcode>[,&nbsp;timeindex]</em>)</em>
<pre>
integer,             intent(in) :: <em class=code>ncFileID</em>
integer,             intent(in) :: <em class=code>domain_id</em>
type(ensemble_type), intent(in) :: <em class=code>state_ens_handle</em>
integer, optional,   intent(in) :: <em class=optionalcode>memberindex</em>
integer, optional,   intent(in) :: <em class=optionalcode>timeindex</em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
Does nothing
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>ncFileID</em></TD>
    <TD>file descriptor to previously-opened netCDF file.</TD></TR>

<TR><TD valign=top><em class=code>domain_id</em></TD>
        <TD>integer describing the domain (which can be a nesting level, 
        a component model ...)</TD></TR>

<TR><TD valign=top><em class=code>state_ens_handle</em></TD>
    <TD>The handle to the state structure containing information about
        the state vector about which information is requested.</TD></TR>

<TR><TD valign=top><em class=code>memberindex</em></TD>
    <TD> Integer index of ensemble member to be written.</TD></TR>

<TR><TD valign=top><em class=code>timeindex</em></TD>
    <TD>The timestep counter for the given state.</TD></TR>

</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="pert_model_copies"></A>
<br>
<div class=routine>
<em class=call>call pert_model_copies(state_ens_handle, ens_size, pert_amp, interf_provided)</em>
<pre>
type(ensemble_type), intent(inout) :: <em class=code>state_ens_handle</em>
integer,             intent(in)    :: <em class=code>ens_size</em>
real(r8),            intent(in)    :: <em class=code>pert_amp</em>
logical,             intent(out)   :: <em class=code>interf_provided</em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
Returns 'interface provided' flag as false, so the default
perturb routine in DART will add small amounts of gaussian noise
to all parts of the state vector.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>state_ens_handle</em></TD>
    <TD>The handle containing an ensemble of state vectors to be perturbed.</TD></TR>

<TR><TD valign=top><em class=code>ens_size</em></TD>
    <TD>The number of ensemble members to perturb.</TD></TR>

<TR><TD valign=top><em class=code>pert_amp</em></TD>
    <TD>the amplitude of the perturbations. The interpretation is based
        on the model-specific implementation.
    </TD></TR>

<TR><TD valign=top><em class=code>interf_provided</em></TD>
    <TD>Returns false if model_mod cannot do this, else true.</TD></TR>

</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<!-- these definitions came from the threed_sphere/location_mod.f90 -->

<A NAME="get_close_obs"></A>
<br>
<div class=routine>
<em class=call>call get_close_obs(gc, base_loc, base_type,
  locs, loc_qtys, loc_types, num_close, close_ind <em class=optionalcode>[,&nbsp;dist] [,&nbsp;state_handle</em>) </em>
<pre>
type(get_close_type),          intent(in)  :: <em class=code>gc</em>
type(location_type),           intent(in)  :: <em class=code>base_loc</em>
integer,                       intent(in)  :: <em class=code>base_type</em>
type(location_type),           intent(in)  :: <em class=code>locs(:)</em>
integer,                       intent(in)  :: <em class=code>loc_qtys(:)</em>
integer,                       intent(in)  :: <em class=code>loc_types(:)</em>
integer,                       intent(out) :: <em class=code>num_close</em>
integer,                       intent(out) :: <em class=code>close_ind(:)</em>
real(r8),            optional, intent(out) :: <em class=optionalcode>dist(:)</em>
type(ensemble_type), optional, intent(in)  :: <em class=optionalcode>state_handle</em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
Passes the call through to the location module code.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>gc</em></TD>
    <TD>The get_close_type which stores precomputed information
        about the locations to speed up searching</TD></TR>

<TR><TD valign=top><em class=code>base_loc</em></TD>
    <TD>Reference location.  The distances will be computed
        between this location and every other location in the obs list</TD></TR>

<TR><TD valign=top><em class=code>base_type</em></TD>
    <TD>The DART quantity at the <em class=code>base_loc</em></TD></TR>

<TR><TD valign=top><em class=code>locs(:)</em></TD>
    <TD>Compute the distance between the <em class=code>base_loc</em> and each
        of the locations in this list</TD></TR>

<TR><TD valign=top><em class=code>loc_qtys(:)</em></TD>
    <TD>The corresponding quantity of each item in the <em class=code>locs</em> list</TD></TR>

<TR><TD valign=top><em class=code>loc_types(:)</em></TD>
    <TD>The corresponding type of each item in the <em class=code>locs</em> list.
        This is not available in the default implementation but may be used in 
        custom implementations.</TD></TR>

<TR><TD valign=top><em class=code>num_close</em></TD>
    <TD>The number of items from the <em class=code>locs</em> list
        which are within maxdist of the base location</TD></TR>

<TR><TD valign=top><em class=code>close_ind(:)</em></TD>
    <TD>The list of index numbers from the <em class=code>locs</em> list 
        which are within maxdist of the base location</TD></TR>

<TR><TD valign=top><em class=code>dist(:)</em></TD>
    <TD>If present, return the distance between each entry
        in the close_ind list and the base location.  If not
        present, all items in the obs list which are closer
        than maxdist will be added to the list but the overhead
        of computing the exact distances will be skipped.</TD></TR>

<TR><TD valign=top><em class=code>state_handle</em></TD>
    <TD>The handle to the state structure containing information about 
        the state vector about which information is requested.</TD></TR>

</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="get_close_state"></A>
<br>
<div class=routine>
<em class=call>call get_close_state(gc, base_loc, base_type,
  state_loc, state_qtys, state_indx, num_close, close_ind, dist, state_handle</em>) </em>
<pre>
type(get_close_type), intent(in)    :: <em class=code>gc</em>
type(location_type),  intent(inout) :: <em class=code>base_loc</em>
integer,              intent(in)    :: <em class=code>base_type</em>
type(location_type),  intent(inout) :: <em class=code>state_loc(:)</em>
integer,              intent(in)    :: <em class=code>state_qtys(:)</em>
integer(i8),          intent(in)    :: <em class=code>state_indx(:)</em>
integer,              intent(out)   :: <em class=code>num_close</em>
integer,              intent(out)   :: <em class=code>close_ind(:)</em>
real(r8),             intent(out)   :: <em class=code>dist(:)</em>
type(ensemble_type),  intent(in)    :: <em class=code>state_handle</em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
Passes the call through to the location module code.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>gc</em></TD>
    <TD>The get_close_type which stores precomputed information
        about the locations to speed up searching</TD></TR>

<TR><TD valign=top><em class=code>base_loc</em></TD>
    <TD>Reference location.  The distances will be computed
        between this location and every other location in the obs list</TD></TR>

<TR><TD valign=top><em class=code>base_type</em></TD>
    <TD>The DART quantity at the  <em class=code>base_loc</em></TD></TR>

<TR><TD valign=top><em class=code>state_loc(:)</em></TD>
    <TD>Compute the distance between the <em class=code>base_loc</em> and each
        of the locations in this list</TD></TR>

<TR><TD valign=top><em class=code>state_qtys(:)</em></TD>
    <TD>The corresponding quantity of each item in the <em class=code>state_loc</em> list</TD></TR>

<TR><TD valign=top><em class=code>state_indx(:)</em></TD>
    <TD>The corresponding DART index of each item in the <em class=code>state_loc</em> list.
        This is not available in the default implementation but may be used in
        custom implementations.</TD></TR>

<TR><TD valign=top><em class=code>num_close</em></TD>
    <TD>The number of items from the <em class=code>state_loc</em> list
        which are within maxdist of the base location</TD></TR>

<TR><TD valign=top><em class=code>close_ind(:)</em></TD>
    <TD>The list of index numbers from the <em class=code>state_loc</em> list 
        which are within maxdist of the base location</TD></TR>

<TR><TD valign=top><em class=code>dist(:)</em></TD>
    <TD>If present, return the distance between each entry
        in the <em class=code>close_ind</em> list and the base location.  If not
        present, all items in the <em class=code>state_loc</em> list which are closer
        than maxdist will be added to the list but the overhead
        of computing the exact distances will be skipped.</TD></TR>

<TR><TD valign=top><em class=code>state_handle</em></TD>
    <TD>The handle to the state structure containing information about 
        the state vector about which information is requested.</TD></TR>

</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="convert_vertical_obs"></A>
<br>
<div class=routine>
<em class=call>call convert_vertical_obs(state_handle, num, locs, loc_qtys, loc_types, which_vert, status)</em>
<pre>
type(ensemble_type), intent(in)  :: <em class=code>state_handle</em>
integer,             intent(in)  :: <em class=code>num</em>
type(location_type), intent(in)  :: <em class=code>locs(:)</em>
integer,             intent(in)  :: <em class=code>loc_qtys(:)</em>
integer,             intent(in)  :: <em class=code>loc_types(:)</em>
integer,             intent(in)  :: <em class=code>which_vert</em>
integer,             intent(out) :: <em class=code>status(:)</em>
</pre>
</div>

<div class=indent1><!-- Description -->

<P>
Passes the call through to the location module code.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>state_handle</em></TD>
    <TD>The handle to the state.</TD></TR>

<TR><TD valign=top><em class=code>num</em></TD>
    <TD>the number of observation locations</TD></TR>

<TR><TD valign=top><em class=code>locs</em></TD>
    <TD>the array of observation locations</TD></TR>

<TR><TD valign=top><em class=code>loc_qtys</em></TD>
    <TD>the array of observation quantities.</TD></TR>

<TR><TD valign=top><em class=code>loc_types</em></TD>
    <TD>the array of observation types.</TD></TR>

<TR><TD valign=top><em class=code>which_vert</em></TD>
    <TD>the desired vertical coordinate system. There is a table
        in the <em class=file>location_mod.f90</em> that relates
        integers to vertical coordinate systems.</TD></TR>

<TR><TD valign=top><em class=code>status</em></TD>
    <TD>Success or failure of the vertical conversion.
        If <em class=code>istatus&nbsp;=&nbsp;0</em>, the conversion was
        a success. Any other value is a failure.
    </TD></TR>

</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->


<A NAME="convert_vertical_state"></A>
<br>
<div class=routine>
<em class=call>call convert_vertical_state(state_handle, num, locs, loc_qtys, loc_types, which_vert, status)</em>
<pre>
type(ensemble_type), intent(in)  :: <em class=code>state_handle</em>
integer,             intent(in)  :: <em class=code>num</em>
type(location_type), intent(in)  :: <em class=code>locs(:)</em>
integer,             intent(in)  :: <em class=code>loc_qtys(:)</em>
integer,             intent(in)  :: <em class=code>loc_types(:)</em>
integer,             intent(in)  :: <em class=code>which_vert</em>
integer,             intent(out) :: <em class=code>status(:)</em>
</pre>
</div>

<div class=indent1><!-- Description -->

<P>
Passes the call through to the location module code.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>state_handle</em></TD>
    <TD>The handle to the state.</TD></TR>

<TR><TD valign=top><em class=code>num</em></TD>
    <TD>the number of state locations</TD></TR>

<TR><TD valign=top><em class=code>locs</em></TD>
    <TD>the array of state locations</TD></TR>

<TR><TD valign=top><em class=code>loc_qtys</em></TD>
    <TD>the array of state quantities.</TD></TR>

<TR><TD valign=top><em class=code>loc_types</em></TD>
    <TD>the array of state types.</TD></TR>

<TR><TD valign=top><em class=code>which_vert</em></TD>
    <TD>the desired vertical coordinate system. There is a table
        in the <em class=file>location_mod.f90</em> that relates
        integers to vertical coordinate systems.</TD></TR>

<TR><TD valign=top><em class=code>status</em></TD>
    <TD>Success or failure of the vertical conversion.
        If <em class=code>istatus&nbsp;=&nbsp;0</em>, the conversion was
        a success. Any other value is a failure.
    </TD></TR>

</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="read_model_time"></A>
<br>
<div class=routine>
<em class=call>model_time = read_model_time(filename)</em>
<pre>
character(len=*), intent(in) :: <em class=code>filename</em>
type(time_type)              :: <em class=code>model_time</em>
</pre>
</div>

<div class=indent1><!-- Description -->

<P>
Passes the call through to the dart_time_io module code.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>filename</em></TD>
    <TD>netCDF file name</TD></TR>

<TR><TD valign=top><em class=code>model_time</em></TD>
    <TD>The current time of the model state.</TD></TR>

</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="write_model_time"></A>
<br>
<div class=routine>
<em class=call>call write_model_time(ncid, dart_time)</em>
<pre>
integer,          intent(in) :: <em class=code>ncid</em>
type(time_type),  intent(in) :: <em class=code>dart_time</em>
</pre>
</div>

<div class=indent1><!-- Description -->

<P>
Passes the call through to the dart_time_io module code.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>ncid</em></TD>
    <TD>handle to an open netCDF file</TD></TR>

<TR><TD valign=top><em class=code>dart_time</em></TD>
    <TD>The current time of the model state.</TD></TR>

</TABLE>

</div>
<br>


<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="end_model"></A>
<br>
<div class=routine>
<em class=call>call end_model()</em>
</div>

<div class=indent1>
<!-- Description -->

<P>
Does nothing.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>
</TABLE>
</div>

<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
<P><!-- useless spacer to make the top behave correctly --></P>
<div class="top">[<a href="#">top</a>]</div><hr />

<A NAME="FilesUsed"></A>
<H2>FILES</H2>

<P>
none
</P>

<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->

<A NAME="References"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>REFERENCES</H2>
<ol>
<li> none </li>
</ol>

<!--==================================================================-->
<!-- Describe all the error conditions and codes.                     -->
<!--==================================================================-->

<A NAME="Errors"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>ERROR CODES and CONDITIONS</H2>

<P>
Standard errors.
</P>

<H2>KNOWN BUGS</H2>
<P>
none at this time
</P>

<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->

<A NAME="FuturePlans"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>FUTURE PLANS</H2>
<P>
none at this time.
</P>

<!--==================================================================-->
<!-- PrivateComponents                                                -->
<!--==================================================================-->

<A NAME="PrivateComponents"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>PRIVATE COMPONENTS</H2>
<P>
N/A
</P>

<!--==================================================================-->
<!-- Legalese & Metadata                                              -->
<!--==================================================================-->

<A NAME="Legalese"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>Terms of Use</H2>

<P>
DART software - Copyright UCAR. This open source software is provided
by UCAR, "as is", without charge, subject to all terms of use at
<a href="http://www.image.ucar.edu/DAReS/DART/DART_download">
http://www.image.ucar.edu/DAReS/DART/DART_download</a>
</P>

<!--==================================================================-->

</BODY>
</HTML>
