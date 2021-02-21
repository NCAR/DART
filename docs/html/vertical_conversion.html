<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>Vertical Conversion of Observations</title>
<link rel="stylesheet" type="text/css" href="../html/doc.css">
<link href="../images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>Vertical Conversion of Observations</h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src="../images/Dartboard7.png" alt=
"DART project logo" height="70"></td>
<td>Jump to <a href="../index.html">DART Documentation Main
Index</a></td>
</tr>
</table>
In Lanai vertical conversion of observations occurs in
get_close_obs. The Lanai code in filter_assim is as follows:
<pre>
SEQUENTIAL_OBS do i = 1, obs_ens_handle%num_vars
   ...
   broadcast increments for observation i
   call get_close_obs()
   call convert_vertical_location(observation(i))
enddo</pre>
<p>If this algorithm was followed in RMA DART, all processors would
have to communicate to calculate the location of observation i.
This is a huge amount of contention since all processors are doing
exactly the same calculation so need to access exactly the same
state elements. This causes the code to run very slowly, for
example 1 minute for 1000 observations, versus 5 seconds for 1000
observations for Lanai.</p>
<p>However, there is no need to calculate the vertical conversion
inside the SEQUENTIAL_OBS do loop, since the mean state vector used
is not updated during the loop. (In Lanai it is the array passed to
the model_mod by ens_mean_for_model in filter_main). Also this
calculation does not scale, because all processors do the same
calculation.</p>
<p>In DART RMA the owner of an observation converts the vertical
location of an observation and broacasts it to all other processors
as part of the broadcast in the SEQUENTIAL_OBS do loop.</p>
<p>The DART RMA code calculates the vertical of all observations
before the loop. This potentially scales better because processors
only calculate their own observation conversions, but does require
model_mod interfaces for vertical conversion.</p>
<p>The DART RMA code in filter_assim is as follows:</p>
<pre>
do i =, obs_ens_handle%my_num_vars
   call convert_vertical_location(my_obs_loc(i))
end do
SEQUENTIAL_OBS do i = 1, obs_ens_handle%num_vars
   ...
   broadcast increments and vertical location for observation i
   ...
enddo
</pre>
<h3>Bitwise Problem</h3>
<p>Moving the <code>convert_vertical_location</code> changes the
number of <code>get/set location</code> calls. There is a bitwise
creep of the location when you do this. This is in the conversion
from degrees to radians and back again. If you want to do the exact
number of <code>get/set location</code> you can change the line
lanai_bitwise = .false. to lanai_bitwise = .true. in
assim_tools_mod.f90. Note this is not a namelist option because
production code should not be run with lanai_bitwise = .true. For
more detail on running bitwise with Lanai see <a href=
"bitwise_considerations.html">bitwise considerations</a>.</p>
<p><!-- make sure the 'top' is aligned correctly --></p>
<!--==================================================================-->
<!-- Legalese & Metadata                                              -->
<!--==================================================================-->
<a name="Legalese" id="Legalese"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>Terms of Use</h2>
<p>DART software - Copyright UCAR. This open source software is
provided by UCAR, "as is", without charge, subject to all terms of
use at <a href=
"http://www.image.ucar.edu/DAReS/DART/DART_download">http://www.image.ucar.edu/DAReS/DART/DART_download</a></p>
<!--==================================================================-->
</body>
</html>
