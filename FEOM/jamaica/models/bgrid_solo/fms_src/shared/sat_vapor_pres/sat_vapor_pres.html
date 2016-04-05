<html>
<head>
<META http-equiv="Content-Type" content="text/html; charset=EUC-JP">
<title>Module sat_vapor_pres_mod</title>
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
<h2>Module sat_vapor_pres_mod</h2>
<a name="HEADER"></a>
<!-- BEGIN HEADER -->
<div>
<b>Contact:&nbsp;</b><a href="mailto:bw@gfdl.noaa.gov">   Bruce Wyman </a>
<br>
<b>Reviewers:&nbsp;</b>
<br>
<b>Change History:&nbsp;</b><a href="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/shared/sat_vapor_pres">WebCVS Log</a>
<br>
<br>
</div>
<!-- END HEADER -->
<a name="OVERVIEW"></a>
<hr>
<h4>OVERVIEW</h4>
<!-- BEGIN OVERVIEW -->
<p class="text">   Routines for determining the saturation vapor pressure (<tt>ES</tt>) 
   and the derivative of <tt>ES</tt>   with respect to temperature. </p>
<!-- END OVERVIEW -->
<a name="DESCRIPTION"></a>
<!-- BEGIN DESCRIPTION -->
<div>   This module contains routines for determining the saturation vapor
   pressure (<tt>ES</tt>) from lookup tables constructed using equations given
   in the Smithsonian tables.  The <tt>ES</tt>   lookup tables are valid between
   -160C and +100C (approx 113K to 373K).
   The values of <tt>ES</tt>   are computed over ice from -160C to -20C,
   over water from 0C to 100C, and a blended value (over water and ice)
   from -20C to 0C.
   This version was written for non-vector machines.
   See the <a href="#NOTES">notes</a>   section for details on vectorization. </div>
<br>
<!-- END DESCRIPTION -->
<a name="OTHER MODULES USED"></a>
<hr>
<h4>OTHER MODULES USED</h4>
<!-- BEGIN OTHER MODULES USED -->
<div>
<pre>constants_mod<br>      fms_mod</pre>
</div>
<!-- END OTHER MODULES USED -->
<!-- BEGIN PUBLIC INTERFACE -->
<a name="PUBLIC INTERFACE"></a>
<hr>
<h4>PUBLIC INTERFACE</h4>
<div>   Description summarizing public interface. <dl>
<dt>
<a href="#lookup_es">lookup_es</a>:</dt>
<dd>   For the given temperatures, returns the saturation vapor pressures. </dd>
<dt>
<a href="#lookup_des">lookup_des</a>:</dt>
<dd>   For the given temperatures, returns the derivative of saturation vapor pressure
   with respect to temperature. </dd>
<dt>
<a href="#compute_es">compute_es</a>:</dt>
<dd>   For the given temperatures, computes the saturation vapor pressures. </dd>
<dt>
<a href="#sat_vapor_pres_init">sat_vapor_pres_init</a>:</dt>
<dd>   Initializes the lookup tables for saturation vapor pressure. </dd>
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
<a name="lookup_es"></a>
<h4>lookup_es</h4>
<pre>
<b>call lookup_es </b>( temp, esat )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   For the given temperatures these routines return the
   saturation vapor pressure (esat). The return values are derived from
   lookup tables (see notes below). </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>temp&nbsp;&nbsp;&nbsp;</tt></td><td>   Temperature in degrees Kelvin. <br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(scalar)]</span>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:)]</span>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:)]</span>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:,:)]</span></td>
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
<td valign="top" align="left"><tt>esat&nbsp;&nbsp;&nbsp;</tt></td><td>   Saturation vapor pressure in pascals.
   May be a scalar, 1d, 2d, or 3d array.
   Must have the same order and size as temp. <br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(scalar)]</span>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:)]</span>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:)]</span>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:,:)]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="lookup_des"></a>
<h4>lookup_des</h4>
<pre>
<b>call lookup_des </b>( temp, desat )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   For the given temperatures these routines return the derivative of esat w.r.t.
   temperature (desat). The return values are derived from
   lookup tables (see notes below). </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>temp&nbsp;&nbsp;&nbsp;</tt></td><td>   Temperature in degrees Kelvin. <br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(scalar)]</span>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:)]</span>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:)]</span>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:,:)]</span></td>
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
<td valign="top" align="left"><tt>desat&nbsp;&nbsp;&nbsp;</tt></td><td>   Derivative of saturation vapor pressure w.r.t. temperature
   in pascals/degree. May be a scalar, 1d, 2d, or 3d array.
   Must have the same order and size as temp. <br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(scalar)]</span>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:)]</span>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:)]</span>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:,:)]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="compute_es"></a>
<h4>compute_es</h4>
<pre>es = <b>compute_es</b> ( temp )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Computes saturation vapor pressure for the given temperature using
   the equations given in the Smithsonian Meteorological Tables.
   Between -20C and 0C a blended value over ice and water is returned. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>temp&nbsp;&nbsp;&nbsp;</tt></td><td>   Temperature in degrees Kelvin. <br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:)]</span>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(scalar)]</span>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:)]</span>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:,:)]</span></td>
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
<td valign="top" align="left"><tt>es&nbsp;&nbsp;&nbsp;</tt></td><td>   Saturation vapor pressure in pascals.
   May be a scalar, 1d, 2d, or 3d array.
   Must have the same order and size as temp. <br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:)]</span>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(scalar)]</span>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:)]</span>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:,:)]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="sat_vapor_pres_init"></a>
<h4>sat_vapor_pres_init</h4>
<pre>
<b>call sat_vapor_pres_init </b>
</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Initializes the lookup tables for saturation vapor pressure.
   This routine will be called automatically the first time <b>lookup_es</b>   or <b>lookup_des</b>   is called,
   the user does not need to call this routine.
   There are no arguments. </dd>
<br>
<br>
</dl>
</li>
</ol>
<!-- END PUBLIC ROUTINES -->
<a name="PUBLIC TYPES"></a>
<!-- BEGIN PUBLIC TYPES -->
<!-- END PUBLIC TYPES --><a name="NAMELIST"></a>
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
<a name="PUBLIC CODE"></a>
<!-- BEGIN PUBLIC CODE -->
<!-- END PUBLIC CODE --><a name="ERROR MESSAGES"></a>
<!-- BEGIN ERROR MESSAGES -->
<hr>
<h4>ERROR MESSAGES</h4>
<div>
<dl>
<dt>
<b>FATAL in lookup_es</b>
</dt>
<dd>
<span class="errmsg">table overflow, nbad=##</span>
</dd>
<dd>   Temperature(s) provided to the saturation vapor pressure lookup
   are outside the valid range of the lookup table (-160 to 100 deg C).
   This may be due to a numerical instability in the model.
   Information should have been printed to standard output to help
   determine where the instability may have occurred.
   If the lookup table needs a larger temperature range,
   then parameters in the module header must be modified. </dd>
</dl>
<br>
</div>
<!-- END ERROR MESSAGES -->
<a name="REFERENCES"></a>
<hr>
<h4>REFERENCES</h4>
<!-- BEGIN REFERENCES -->
<div>
<ol>
<li>   Smithsonian Meteorological Tables Page 350. </li>
</ol>
</div>
<br>
<!-- END REFERENCES -->
<a name="COMPILER SPECIFICS"></a>
<hr>
<h4>COMPILER SPECIFICS</h4>
<!-- BEGIN COMPILER SPECIFICS -->
<div>
        None.
      </div>
<br>
<!-- END COMPILER SPECIFICS -->
<a name="PRECOMPILER OPTIONS"></a>
<hr>
<h4>PRECOMPILER OPTIONS</h4>
<!-- BEGIN PRECOMPILER OPTIONS -->
<div>
        None.
      </div>
<br>
<!-- END PRECOMPILER OPTIONS -->
<a name="LOADER OPTIONS"></a>
<hr>
<h4>LOADER OPTIONS</h4>
<!-- BEGIN LOADER -->
<div>None.<br>
<br>
</div>
<!-- END LOADER OPTIONS -->
<a name="TEST PROGRAM"></a>
<hr>
<h4>TEST PROGRAM</h4>
<!-- BEGIN TEST PROGRAM -->
<div>
<dl>
<dt>test_sat_vapor_pres</dt>
<dd> 
<pre>use sat_vapor_pres_mod
implicit none

integer, parameter :: ipts=500, jpts=100, kpts=50, nloop=1
real, dimension(ipts,jpts,kpts) :: t,es,esn,des,desn
integer :: n

 generate temperatures between 120K and 340K
  call random_number (t)
  t = 130. + t * 200.

 initialize the tables (optional)
  call sat_vapor_pres_init

 compute actual es and "almost" actual des
   es = compute_es  (t)
  des = compute_des (t)

do n = 1, nloop
 es and des
  call lookup_es  (t, esn)
  call lookup_des (t,desn)
enddo

 terminate, print deviation from actual
  print *, 'size=',ipts,jpts,kpts,nloop
  print *, 'err es  = ', sum((esn-es)**2)
  print *, 'err des = ', sum((desn-des)**2)

contains

----------------------------------
 routine to estimate derivative

 function compute_des (tem) result (des)
 real, intent(in) :: tem(:,:,:)
 real, dimension(size(tem,1),size(tem,2),size(tem,3)) :: des,esp,esm
 real, parameter :: tdel = .01
    esp = compute_es (tem+tdel)
    esm = compute_es (tem-tdel)
    des = (esp-esm)/(2*tdel)
 end function compute_des
----------------------------------

end program test_sat_vapor_pres</pre> 
</dd>
</dl>
</div>
<br>
<!-- END TEST PROGRAM -->
<a name="KNOWN BUGS"></a>
<hr>
<h4>KNOWN BUGS</h4>
<!-- BEGIN KNOWN BUGS -->
<div>
<p>   No error checking is done to make sure that the size of the
   input and output fields match. </p>
</div>
<br>
<!-- END KNOWN BUGS -->
<a name="NOTES"></a>
<hr>
<h4>NOTES</h4>
<!-- BEGIN NOTES -->
<div>   1. <b>Vectorization</b>
<br>   To create a vector version the lookup routines need to be modified.
   The local variables: tmp, del, ind, should be changed to arrays
   with the same size and order as input array temp.
   <br>
<br>
   2. <b>Construction of the <tt>ES</tt>   tables</b>
<br>   The tables are constructed using the saturation vapor pressure (<tt>ES</tt>)
   equations in the Smithsonian tables. The tables are valid between
   -160C to +100C with increments at 1/10 degree. Between -160C and -20C
   values of <tt>ES</tt>   over ice are used, between 0C and 100C values of<tt>   ES</tt>   over water are used, between -20C and 0C blended values of <tt>ES</tt>   (over water and over ice) are used.
   <br>
<br>
   There are three tables constructed: <tt>ES</tt>, first derivative 
   (<tt>ES'</tt>), and
   second derivative (<tt>ES</tt>'').  The ES table is constructed directly from
   the equations in the Smithsonian tables. The <tt>ES</tt>' table is constructed
   by bracketing temperature values at +/- 0.01 degrees. The <tt>ES</tt>'' table
   is estimated by using centered differencing of the <tt>ES</tt>' table.
   <br>
<br>
   3. <b>Determination of <tt>es</tt>   and <tt>es'</tt>   from lookup tables</b>
<br>   Values of the saturation vapor pressure (<tt>es</tt>) and the 
   derivative (<tt>es'</tt>) are determined at temperature (T) from the lookup 
   tables (<tt>ES</tt>, <tt>ES'</tt>, <tt>ES''</tt>)
   using the following formula. <pre>    es (T) = ES(t) + ES'(t) * dt + 0.5 * ES''(t) * dt**2
    es'(T) = ES'(t) + ES''(t) * dt

    where     t = lookup table temperature closest to T
             dt = T - t</pre>   4. Internal (private) parameters<br>   These parameters can be modified to increase/decrease the size/range
   of the lookup tables. <pre>    tcmin   The minimum temperature (in deg C) in the lookup tables.
              [integer, default: tcmin = -160]

    tcmax   The maximum temperature (in deg C) in the lookup tables.
              [integer, default: tcmin = +100]</pre> 
</div>
<br>
<!-- END NOTES -->
<a name="FUTURE PLANS"></a>
<hr>
<h4>FUTURE PLANS</h4>
<!-- BEGIN FUTURE PLANS -->
<div>
        None.
      </div>
<br>
<!-- END FUTURE PLANS -->
<hr>
<div align="right">
<font size="-1"><a href="#TOP">top</a></font>
</div>
</body>
</html>
