<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>module obs_def_rttov_mod</title>
<link rel="stylesheet" type="text/css" href=
"../../docs/html/doc.css">
<link href="../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>MODULE <em class="program">obs_def_rttov_mod</em></h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src="../../docs/images/Dartboard7.png"
alt="DART project logo" height="70"></td>
<td>Jump to <a href="../../docs/index.html">DART Documentation Main
Index</a></td>
</tr>
</table>
<a href="#Interface">INTERFACES</a> / <a href=
"#Namelist">NAMELIST</a> / <a href="#FilesUsed">FILES</a> /
<a href="#References">REFERENCES</a> / <a href="#Errors">ERRORS</a>
/ <a href="#FuturePlans">PLANS</a> / <a href=
"#PrivateComponents">PRIVATE COMPONENTS</a> / <a href=
"#Legalese">TERMS OF USE</a>
<h2>Overview</h2>
<p>DART RTTOV observation module, including the observation
operators for the two primary RTTOV-observation types --
visible/infrared radiances and microwave radiances/brightness
temperatures.<br>
<br>
This module acts as a pass-through for RTTOV version 12.3. For more
information, see <a href=
"https://www.nwpsaf.eu/site/software/rttov/documentation/">the
RTTOV site</a>.<br>
<br>
DART supports both RTTOV-direct for visible/infrared/microwave as
well as RTTOV-scatt for microwave computations. The code, in
principle, supports all features of version 12.3 as a pass-through
from the model to RTTOV, includes aerosols, trace gases, clouds,
and atmospheric variables. The code also includes directly
specifying scattering properties.<br>
<br>
However, a model may not have all of the variables necessary for
these functions depending on your model's setup.<br>
<br>
For example, DART can use any of the RTTOV clw or ice schemes, but
the WRF model is not directly compatible with the IR default cloud
classification of marine/continental stratus/cumulus clean/dirty.
We also offer a simple classification based on maximum vertical
velocity in the column and land type, but due to lack of aerosol
information, WRF/DART cannot differentiate between clean and dirty
cumulus. This may have some impact on the forward calculations -
but in experience the difference in cloud phase (ice versus water)
makes a much larger difference.<br>
<br>
Trace gases and aerosols may be important for actual observation
system experiments using visible/infrared; this may depend on the
precise frequencies you wish to use. Although a model may not have
the necessary inputs by itself, although the defaults in RTTOV
based on climatology can be used. The impact on the quality of the
results should be investigated.<br>
<br>
Known issues:</p>
<ul>
<li>DART does not yet provide any type of bias correction</li>
<li>Cross-channel error correlations are not yet supported. A
principal component approach has been discussed. For now, the best
bet is to use a subset of channels that are nearly independent of
one another.</li>
<li>Vertical localization will need to be tuned. Turning off
vertical localization may work well if you have a large number of
ensemble members. Using the maximum peak of the weighting function
or the cloud-top may be appropriate. There are also other potential
approaches being investigated.</li>
</ul>
<br>
<br>
Author and Contact information:
<ul>
<li>DART Code: Jeff Steward, jsteward at ucar.edu</li>
<li>Original DART/RTTOV work: Nancy Collins, Johnny Hendricks</li>
</ul>
<h3>Backward compatibility note:</h3>
<!--==================================================================-->
<a name="OtherModulesUsed" id="OtherModulesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>OTHER MODULES USED</h2>
<pre>
types_mod
utilities_mod
location_mod (threed_sphere)
assim_model_mod
obs_def_utilitie_mod
ensemble_manager_mod
utilities_mod
parkind1 (from RTTOV)
rttov_types (from RTTOV)
obs_kind_mod
</pre>
<!--==================================================================-->
<!-- Declare all public entities ...                                  -->
<!-- duplicate public routines template as many times as necessary    -->
<!-- make sure you replace all yyyroutine?? strings                   -->
<!--==================================================================-->
<!--Note to authors. The first row of the table is different.         -->
<!--==================================================================-->
<a name="Interface" id="Interface"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>PUBLIC INTERFACES</h2>
<table>
<tr>
<td><em class="call">use obs_def_rttov_mod, only :</em></td>
<td><a href="#set_visir_metadata">set_visir_metadata</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#set_mw_metadata">set_mw_metadata</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_expected_radiance">get_expected_radiance</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#get_rttov_option_logical">get_rttov_option_logical</a></td>
</tr>
</table>
<p>Namelist interface <a href="#Namelist"><em class=
"code">&amp;obs_def_rttov_mod_nml</em></a> is read from file
<em class="file">input.nml</em>.</p>
<p>A note about documentation style. Optional arguments are
enclosed in brackets <em class="optionalcode">[like this]</em>.</p>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="set_visir_metadata" id="set_visir_metadata"></a><br>
<div class="routine"><em class="call">call set_visir_metadata(key,
sat_az, sat_ze, sun_az, sun_ze, &amp; platform_id, sat_id,
sensor_id, channel, specularity)</em>
<pre>
integer,  intent(out) :: <em class="code">key</em>
real(r8), intent(in)  :: <em class="code">sat_az</em>
real(r8), intent(in)  :: <em class="code">sat_ze</em>
real(r8), intent(in)  :: <em class="code">sun_az</em>
real(r8), intent(in)  :: <em class="code">sun_ze</em>
integer,  intent(in)  :: <em class=
"code">platform_id, sat_id, sensor_id, channel</em>
real(r8), intent(in)  :: <em class="code">specularity</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Visible / infrared observations have several auxillary metadata
variables. Other than the key, which is standard DART fare, the
RTTOV satellite azimuth and satellite zenith angle must be
specified. See the RTTOV user guide for more information (in
particular, see figure 4). If the <em class="code">addsolar</em>
namelist value is set to true, then the solar azimuth and solar
zenith angles must be specified - again see the RTTOV user guide.
In addition to the platform/satellite/ sensor ID numbers, which are
the RTTOV unique identifiers, the channel specifies the chanenl
number in the RTTOV coefficient file. Finally, if <em class=
"code">do_lambertian</em> is true, specularity must be specified
here. Again, see the RTTOV user guide for more information.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">key  </em></td>
<td>The DART observation key.</td>
</tr>
<tr>
<td valign="top"><em class="code">sat_az</em></td>
<td>The satellite azimuth angle.</td>
</tr>
<tr>
<td valign="top"><em class="code">sat_ze</em></td>
<td>The satellite zenith angle.</td>
</tr>
<tr>
<td valign="top"><em class="code">sun_az</em></td>
<td>The solar azimuth angle. Only relevant if addsolar is
true.</td>
</tr>
<tr>
<td valign="top"><em class="code">sun_ze</em></td>
<td>The solar zenith angle. Only relevant if addsolar is true.</td>
</tr>
<tr>
<td valign="top"><em class="code">platform_id</em></td>
<td>The RTTOV platform ID.</td>
</tr>
<tr>
<td valign="top"><em class="code">sat_id</em></td>
<td>The RTTOV satellite ID.</td>
</tr>
<tr>
<td valign="top"><em class="code">sensor_id</em></td>
<td>The RTTOV sensor ID.</td>
</tr>
<tr>
<td valign="top"><em class="code">channel</em></td>
<td>The RTTOV channel number.</td>
</tr>
<tr>
<td valign="top"><em class="code">specularity</em></td>
<td>The surface specularity. Only relevant if do_lambertian is
true.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="set_mw_metadata" id="set_mw_metadata"></a><br>
<div class="routine"><em class="call">call set_mw_metadata(key,
sat_az, sat_ze, platform_id, sat_id, sensor_id, channel, mag_field,
cosbk, fastem_p1, fastem_p2, fastem_p3, fastem_p4, fastem_p5)</em>
<pre>
integer,  intent(out) :: <em class="code">key</em>
real(r8), intent(in)  :: <em class="code">sat_az</em>
real(r8), intent(in)  :: <em class="code">sat_ze</em>
integer,  intent(in)  :: <em class=
"code">platform_id, sat_id, sensor_id, channel</em>
real(r8), intent(in)  :: <em class="code">mag_field</em>
real(r8), intent(in)  :: <em class="code">cosbk</em>
real(r8), intent(in)  :: <em class="code">fastem_p[1-5]</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Microwave observations have several auxillary metadata
variables. Other than the key, which is standard DART fare, the
RTTOV satellite azimuth and satellite zenith angle must be
specified. See the RTTOV user guide for more information (in
particular, see figure 4). In addition to the platform/satellite/
sensor ID numbers, which are the RTTOV unique identifiers, the
channel specifies the chanenl number in the RTTOV coefficient file.
In addition, if <em class="code">use_zeeman</em> is true, the
magnetic field and cosine of the angle between the magnetic field
and angle of propagation must be specified. See the RTTOV user
guide for more information. Finally, the fastem parameters for land
must be specified here. This may be difficult for observations to
set, so default values (see table 21 in the RTTOV user guide) can
be used until a better solution is devised.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">key  </em></td>
<td>The DART observation key.</td>
</tr>
<tr>
<td valign="top"><em class="code">sat_az</em></td>
<td>The satellite azimuth angle.</td>
</tr>
<tr>
<td valign="top"><em class="code">sat_ze</em></td>
<td>The satellite zenith angle.</td>
</tr>
<tr>
<td valign="top"><em class="code">platform_id</em></td>
<td>The RTTOV platform ID.</td>
</tr>
<tr>
<td valign="top"><em class="code">sat_id</em></td>
<td>The RTTOV satellite ID.</td>
</tr>
<tr>
<td valign="top"><em class="code">sensor_id</em></td>
<td>The RTTOV sensor ID.</td>
</tr>
<tr>
<td valign="top"><em class="code">channel</em></td>
<td>The RTTOV channel number.</td>
</tr>
<tr>
<td valign="top"><em class="code">mag_field</em></td>
<td>The strength of the magnetic field. Only relevant if add_zeeman
is true.</td>
</tr>
<tr>
<td valign="top"><em class="code">cosbk</em></td>
<td>The cosine of the angle between the magnetic field and
direction of EM propagation. Only relevant if add_zeeman is
true.</td>
</tr>
<tr>
<td valign="top"><em class="code">fastem_p[1-5]</em></td>
<td>The five parameters used for fastem land/sea ice emissivities.
For ocean emissivities, an internal model is used based on the
value of fastem_version.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_expected_radiance" id=
"get_expected_radiance"></a><br>
<div class="routine"><em class="call">call
get_expected_radiance(obs_kind_ind, state_handle, ens_size,
location, key, val, istatus)</em>
<pre>
integer,             intent(in)  :: <em class=
"code">obs_kind_ind</em>
type(ensemble_type), intent(in)  :: <em class=
"code">state_handle</em>
integer,             intent(in)  :: <em class="code">ens_size</em>
type(location_type), intent(in)  :: <em class="code">location</em>
integer,             intent(in)  :: <em class="code">key</em>
real(r8),            intent(out) :: <em class=
"code">val(ens_size)</em>
integer,             intent(out) :: <em class=
"code">istatus(ens_size)</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Given a location and the state vector from one of the ensemble
members, compute the model-predicted satellite observation. This
can be either in units of radiance (mW/cm-1/sr/sq.m) or a
brightness temperature (in K), depending on if this is a
visible/infrared observation or a microwave observation.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"code">obs_kind_ind  </em></td>
<td>The index of the observation kind; since many observation kinds
are handled by this module, this can be used to determine precisely
which observation kind is being used.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">state_handle  </em></td>
<td>The ensemble of model states to be used for the observation
operator calculations.</td>
</tr>
<tr>
<td valign="top"><em class="code">location</em></td>
<td>Location of this observation</td>
</tr>
<tr>
<td valign="top"><em class="code">key</em></td>
<td>Unique identifier associated with this satellite
observation</td>
</tr>
<tr>
<td valign="top"><em class="code">val</em></td>
<td>The returned observation in units of either radiance or
brightness temperature.</td>
</tr>
<tr>
<td valign="top"><em class="code">istatus</em></td>
<td>Returned integer status code describing problems with applying
forward operator. 0 is a good value; any positive value indicates
an error; negative values are reserved for internal DART use
only.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_rttov_option_logical" id=
"get_rttov_option_logical"></a><br>
<div class="routine"><em class="call">p =
get_rttov_option_logical(field_name)</em>
<pre>
character(len=*),           intent(in)  :: <em class=
"code">field_name</em>
logical,                    result      :: <em class="code">p</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Return the logical value of the RTTOV parameter associated with
the field_name.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">field_name  </em></td>
<td>The name of the RTTOV parameter from the namelist.</td>
</tr>
<tr>
<td valign="top"><em class="code">p</em></td>
<td>The logical return value associated with the parameter.</td>
</tr>
</table>
</div>
<br>
<!--==================================================================-->
<!--=================== DESCRIPTION OF A NAMELIST  ===================-->
<!--==================================================================-->
 <a name="Namelist" id="Namelist"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>NAMELIST</h2>
<p>This namelist is read from the file <em class=
"file">input.nml</em>. Namelists start with an ampersand '&amp;'
and terminate with a slash '/'. Character strings that contain a
'/' must be enclosed in quotes to prevent them from prematurely
terminating the namelist.</p>
<div class="namelist">
<pre>
&amp;obs_def_rttov_nml
   rttov_sensor_db_file   = 'rttov_sensor_db.csv'
   first_lvl_is_sfc       = .true. 
   mw_clear_sky_only      = .false.
   interp_mode            = 1 
   do_checkinput          = .true.
   apply_reg_limits       = .true.
   verbose                = .true.
   fix_hgpl               = .false.
   do_lambertian          = .false.
   lambertian_fixed_angle = .true.
   rad_down_lin_tau       = .true.
   use_q2m                = .true.
   use_uv10m              = .true.
   use_wfetch             = .false.
   use_water_type         = .false.
   addrefrac              = .false.
   plane_parallel         = .false.
   use_salinity           = .false.
   apply_band_correction  = .true.
   cfrac_data             = .true.
   clw_data               = .true.
   rain_data              = .true.
   ciw_data               = .true.
   snow_data              = .true.
   graupel_data           = .true.
   hail_data              = .false.
   w_data                 = .true.
   clw_scheme             = 1
   clw_cloud_top          = 322.
   fastem_version         = 6
   supply_foam_fraction   = .false.
   use_totalice           = .true.
   use_zeeman             = .false.
   cc_threshold           = 0.05
   ozone_data             = .false.
   co2_data               = .false.
   n2o_data               = .false.
   co_data                = .false.
   ch4_data               = .false.
   so2_data               = .false.
   addsolar               = .false.
   rayleigh_single_scatt  = .true.
   do_nlte_correction     = .false.
   solar_sea_brdf_model   = 2
   ir_sea_emis_model      = 2
   use_sfc_snow_frac      = .false.
   add_aerosl             = .false.
   aerosl_type            = 1
   add_clouds             = .true.
   ice_scheme             = 1
   use_icede              = .false.
   idg_scheme             = 2
   user_aer_opt_param     = .false.
   user_cld_opt_param     = .false.
   grid_box_avg_cloud     = .true.
   cldstr_threshold       = -1.0
   cldstr_simple          = .false.
   cldstr_low_cloud_top   = 750.0
   ir_scatt_model         = 2
   vis_scatt_model        = 1
   dom_nstreams           = 8
   dom_accuracy           = 0.0
   dom_opdep_threshold    = 0.0
   addpc                  = .false.
   npcscores              = -1
   addradrec              = .false.
   ipcreg                 = 1
   use_htfrtc             = .false.
   htfrtc_n_pc            = -1
   htfrtc_simple_cloud    = .false.
   htfrtc_overcast        = .false.
/
</pre></div>
<br>
<br>
<div>1
<table border="0" cellpadding="10" width="100%" summary=
'namelist description'>
<thead align="left">
<tr>
<th>Item</th>
<th>Type</th>
<th>Description</th>
</tr>
</thead>
<tbody valign="top">
<tr>
<td>rttov_sensor_db_file</td>
<td>character(len=512)</td>
<td>The location of the RTTOV sensor database. The format for the
database is a comma-separated file. The columns of the database are
the DART observation-kind, the platform/satellite/sensor ID, the
observation type, the coefficient file, and a comma-separated list
of RTTOV channels to use for this observation type.</td>
</tr>
<tr>
<td>first_lvl_is_sfc</td>
<td>logical</td>
<td>Whether the first level of the model represents the surface
(true) or the top of the atmosphere (false).</td>
</tr>
<tr>
<td>mw_clear_sky_only</td>
<td>logical</td>
<td>If microwave calculations should be "clear-sky" only (although
cloud-liquid water absorption/emission is considered; see the RTTOV
user guide).</td>
</tr>
<tr>
<td>interp_mode</td>
<td>integer</td>
<td>The interpolation mode (see the RTTOV user guide).</td>
</tr>
<tr>
<td>do_checkinput</td>
<td>logical</td>
<td>Whether to check the input for reasonableness (see the RTTOV
user guide).</td>
</tr>
<tr>
<td>apply_reg_limits</td>
<td>logical</td>
<td>Whether to clamp the atmospheric values to the RTTOV bounds
(see the RTTOV user guide).</td>
</tr>
<tr>
<td>verbose</td>
<td>logical</td>
<td>Whether to output lots of additional output (see the RTTOV user
guide).</td>
</tr>
<tr>
<td>fix_hgpl</td>
<td>logical</td>
<td>Whether the surface pressure represents the surface or the 2
meter value (see the RTTOV user guide).</td>
</tr>
<tr>
<td>do_lambertian</td>
<td>logical</td>
<td>Whether to include the effects of surface specularity (see the
RTTOV user guide).</td>
</tr>
<tr>
<td>lambertian_fixed_angle</td>
<td>logical</td>
<td>Whether to include a fixed angle for the lambertian effect (see
the RTTOV user guide).</td>
</tr>
<tr>
<td>rad_down_lin_tau</td>
<td>logical</td>
<td>Whether to use the linear-in-tau approximation (see the RTTOV
user guide).</td>
</tr>
<tr>
<td>use_q2m</td>
<td>logical</td>
<td>Whether to use 2m humidity information (see the RTTOV user
guide). If true, the QTY_2M_SPECIFIC_HUMIDITY will be requested
from the model.</td>
</tr>
<tr>
<td>use_q2m</td>
<td>logical</td>
<td>Whether to use 2m humidity information (see the RTTOV user
guide). If true, the QTY_2M_SPECIFIC_HUMIDITY will be requested
from the model.</td>
</tr>
<tr>
<td>use_uv10m</td>
<td>logical</td>
<td>Whether to use 10m wind speed information (see the RTTOV user
guide). If true, the QTY_10M_U_WIND_COMPONENT and
QTY_10M_V_WIND_COMPONENTS will be requested from the model.</td>
</tr>
<tr>
<td>use_wfetch</td>
<td>logical</td>
<td>Whether to use wind fetch information (see the RTTOV user
guide). If true, the QTY_WIND_FETCH will be requested from the
model.</td>
</tr>
<tr>
<td>use_water_type</td>
<td>logical</td>
<td>Whether to use water-type information (0 = fresh, 1 = ocean;
see the RTTOV user guide). If true, the QTY_WATER_TYPE will be
requested from the model.</td>
</tr>
<tr>
<td>addrefrac</td>
<td>logical</td>
<td>Whether to enable atmospheric refraction (see the RTTOV user
guide).</td>
</tr>
<tr>
<td>plane_parallel</td>
<td>logical</td>
<td>Whether to treat the atmosphere as plane parallel (see the
RTTOV user guide).</td>
</tr>
<tr>
<td>use_salinity</td>
<td>logical</td>
<td>Whether to use salinity (see the RTTOV user guide). If true,
the QTY_SALINITY will be requested from the model.</td>
</tr>
<tr>
<td>apply_band_correction</td>
<td>logical</td>
<td>Whether to apply band correction from the coefficient field for
microwave data (see the RTTOV user guide).</td>
</tr>
<tr>
<td>cfrac_data</td>
<td>logical</td>
<td>Whether to use the cloud fraction from 0 to 1 (see the RTTOV
user guide). If true, the QTY_CLOUD_FRACTION will be requested from
the model.</td>
</tr>
<tr>
<td>clw_data</td>
<td>logical</td>
<td>Whether to use cloud-liquid water data (see the RTTOV user
guide). If true, the QTY_CLOUDWATER_MIXING_RATIO will be requested
from the model.</td>
</tr>
<tr>
<td>rain_data</td>
<td>logical</td>
<td>Whether to use precipitating water data (see the RTTOV user
guide). If true, the QTY_RAINWATER_MIXING_RATIO will be requested
from the model.</td>
</tr>
<tr>
<td>ciw_data</td>
<td>logical</td>
<td>Whether to use non-precipiting ice information (see the RTTOV
user guide). If true, the QTY_ICE_MIXING_RATIO will be requested
from the model.</td>
</tr>
<tr>
<td>snow_data</td>
<td>logical</td>
<td>Whether to use precipitating fluffy ice (see the RTTOV user
guide). If true, the QTY_SNOW_MIXING_RATIO will be requested from
the model.</td>
</tr>
<tr>
<td>graupel_data</td>
<td>logical</td>
<td>Whether to use precipting small, hard ice (see the RTTOV user
guide). If true, the QTY_GRAUPEL_MIXING_RATIO will be requested
from the model.</td>
</tr>
<tr>
<td>hail_data</td>
<td>logical</td>
<td>Whether to use precipitating large, hard ice (see the RTTOV
user guide). If true, the QTY_HAIL_MIXING_RATIO will be requested
from the model.</td>
</tr>
<tr>
<td>w_data</td>
<td>logical</td>
<td>Whether to use vertical velocity information. This will be used
to crudely classify if a cloud is cumulus or stratiform for the
purpose of visible/infrared calculations. If true, the
QTY_VERTICAL_VELOCITY will be requested from the model.</td>
</tr>
<tr>
<td>clw_scheme</td>
<td>integer</td>
<td>The clw_scheme to use (see the RTTOV user guide).</td>
</tr>
<tr>
<td>clw_cloud_top</td>
<td>real(r8)</td>
<td>Lower hPa limit for clw calculations (see the RTTOV user
guide).</td>
</tr>
<tr>
<td>fastem_version</td>
<td>integer</td>
<td>Which FASTEM version to use (see the RTTOV user guide).</td>
</tr>
<tr>
<td>supply_foam_fraction</td>
<td>logical</td>
<td>Whether to use sea-surface foam fraction (see the RTTOV user
guide). If true, the QTY_FOAM_FRAC will be requested from the
model.</td>
</tr>
<tr>
<td>use_totalice</td>
<td>logical</td>
<td>Whether to use totalice instead of precip/non-precip ice for
microwave (see the RTTOV user guide).</td>
</tr>
<tr>
<td>use_zeeman</td>
<td>logical</td>
<td>Whether to use the Zeeman effect (see the RTTOV user guide). If
true, the magnetic field and cosine of bk will be used from the
observation metadata.</td>
</tr>
<tr>
<td>cc_threshold</td>
<td>real(r8)</td>
<td>Cloud-fraction value to treat as clear-sky (see the RTTOV user
guide).</td>
</tr>
<tr>
<td>ozone_data</td>
<td>logical</td>
<td>Whether to use ozone (O3) profiles (see the RTTOV user guide).
If true, the QTY_O3 will be requested from the model.</td>
</tr>
<tr>
<td>co2_data</td>
<td>logical</td>
<td>Whether to use carbon dioxide (CO2) profiles (see the RTTOV
user guide). If true, the QTY_CO2 will be requested from the
model.</td>
</tr>
<tr>
<td>n2o_data</td>
<td>logical</td>
<td>Whether to use nitrous oxide (N2O) profiles (see the RTTOV user
guide). If true, the QTY_N2O will be requested from the model.</td>
</tr>
<tr>
<td>co_data</td>
<td>logical</td>
<td>Whether to use carbon monoxide (CO) profiles (see the RTTOV
user guide). If true, the QTY_CO will be requested from the
model.</td>
</tr>
<tr>
<td>ch4_data</td>
<td>logical</td>
<td>Whether to use methane (CH4) profiles (see the RTTOV user
guide). If true, the QTY_CH4 will be requested from the model.</td>
</tr>
<tr>
<td>so2_data</td>
<td>logical</td>
<td>Whether to use sulfur dioxide (SO2) (see the RTTOV user guide).
If true, the QTY_SO2 will be requested from the model.</td>
</tr>
<tr>
<td>addsolar</td>
<td>logical</td>
<td>Whether to use solar angles (see the RTTOV user guide). If
true, the sun_ze and sun_az from the observation metadata will be
used for visible/infrared.</td>
</tr>
<tr>
<td>rayleigh_single_scatt</td>
<td>logical</td>
<td>Whether to use only single scattering for Rayleigh scattering
for visible calculations (see the RTTOV user guide).</td>
</tr>
<tr>
<td>do_nlte_correction</td>
<td>logical</td>
<td>Whether to include non-LTE bias correction for HI-RES sounder
(see the RTTOV user guide).</td>
</tr>
<tr>
<td>solar_sea_brdf_model</td>
<td>integer</td>
<td>The solar sea BRDF model to use (see the RTTOV user
guide).</td>
</tr>
<tr>
<td>ir_sea_emis_model</td>
<td>logical</td>
<td>The infrared sea emissivity model to use (see the RTTOV user
guide).</td>
</tr>
<tr>
<td>use_sfc_snow_frac</td>
<td>logical</td>
<td>Whether to use the surface snow fraction (see the RTTOV user
guide). If true, the QTY_SNOWCOVER_FRAC will be requested from the
model.</td>
</tr>
<tr>
<td>add_aerosl</td>
<td>logical</td>
<td>Whether to use aerosols (see the RTTOV user guide).</td>
</tr>
<tr>
<td>aerosl_type</td>
<td>integer</td>
<td>Whether to use OPAC or CAMS aerosols (see the RTTOV user
guide).</td>
</tr>
<tr>
<td>add_clouds</td>
<td>logical</td>
<td>Whether to enable cloud scattering for visible/infrared (see
the RTTOV user guide).</td>
</tr>
<tr>
<td>ice_scheme</td>
<td>integer</td>
<td>The ice scheme to use (see the RTTOV user guide).</td>
</tr>
<tr>
<td>use_icede</td>
<td>logical</td>
<td>Whether to use the ice effective diameter for visible/infrared
(see the RTTOV user guide). If true, the QTY_CLOUD_ICE_DE will be
requested from the model.</td>
</tr>
<tr>
<td>idg_scheme</td>
<td>integer</td>
<td>The ice water effective diameter scheme to use (see the RTTOV
user guide).</td>
</tr>
<tr>
<td>user_aer_opt_param</td>
<td>logical</td>
<td>Whether to directly specify aerosol scattering properties (see
the RTTOV user guide). Not yet supported.</td>
</tr>
<tr>
<td>user_cld_opt_param</td>
<td>logical</td>
<td>Whether to directly specify cloud scattering properties (see
the RTTOV user guide). Not yet supported.</td>
</tr>
<tr>
<td>grid_box_avg_cloud</td>
<td>logical</td>
<td>Whether to cloud concentrations are grid box averages (see the
RTTOV user guide).</td>
</tr>
<tr>
<td>cldstr_threshold</td>
<td>real(r8)</td>
<td>Threshold for cloud stream weights for scattering (see the
RTTOV user guide).</td>
</tr>
<tr>
<td>cldstr_simple</td>
<td>logical</td>
<td>Whether to use one clear and one cloudy column (see the RTTOV
user guide).</td>
</tr>
<tr>
<td>cldstr_low_cloud_top</td>
<td>real(r8)</td>
<td>Cloud fraction maximum in layers from the top of the atmosphere
down to the specified hPa (see the RTTOV user guide).</td>
</tr>
<tr>
<td>ir_scatt_model</td>
<td>integer</td>
<td>Which infrared scattering method to use (see the RTTOV user
guide).</td>
</tr>
<tr>
<td>vis_scatt_model</td>
<td>integer</td>
<td>Which visible scattering method to use (see the RTTOV user
guide).</td>
</tr>
<tr>
<td>dom_nstreams</td>
<td>integer</td>
<td>The number of streams to use with DOM (see the RTTOV user
guide).</td>
</tr>
<tr>
<td>dom_accuracy</td>
<td>real(r8)</td>
<td>The convergence criteria for DOM (see the RTTOV user
guide).</td>
</tr>
<tr>
<td>dom_opdep_threshold</td>
<td>real(r8)</td>
<td>Ignore layers below this optical depth (see the RTTOV user
guide).</td>
</tr>
<tr>
<td>addpc</td>
<td>logical</td>
<td>Whether to do principal component calculations (see the RTTOV
user guide).</td>
</tr>
<tr>
<td>npcscores</td>
<td>integer</td>
<td>Number of principal components to use for addpc (see the RTTOV
user guide).</td>
</tr>
<tr>
<td>addradrec</td>
<td>logical</td>
<td>Reconstruct the radiances using addpc (see the RTTOV user
guide).</td>
</tr>
<tr>
<td>ipcreg</td>
<td>integer</td>
<td>Number of predictors to use with addpc (see the RTTOV user
guide).</td>
</tr>
<tr>
<td>use_htfrtc</td>
<td>logical</td>
<td>Whether to use HTFRTC (see the RTTOV user guide).</td>
</tr>
<tr>
<td>htfrtc_n_pc</td>
<td>integer</td>
<td>Number of PCs to use with HTFRTC (see the RTTOV user
guide).</td>
</tr>
<tr>
<td>htfrtc_simple_cloud</td>
<td>logical</td>
<td>Whether to use simple cloud scattering with htfrtc (see the
RTTOV user guide).</td>
</tr>
<tr>
<td>htfrtc_overcast</td>
<td>logical</td>
<td>Whether to calculate overcast radiances with HTFRTC (see the
RTTOV user guide).</td>
</tr>
</tbody>
</table>
</div>
<br>
<br>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
 <a name="FilesUsed" id="FilesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FILES</h2>
<ul>
<li>A DART observation sequence file containing Radar obs.</li>
</ul>
<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->
<a name="References" id="References"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>REFERENCES</h2>
<ul>
<li><a href=
"https://www.nwpsaf.eu/site/software/rttov/documentation/">RTTOV
user guide</a></li>
</ul>
<!--==================================================================-->
<!-- Describe all the error conditions and codes.                     -->
<!--==================================================================-->
<a name="Errors" id="Errors"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>ERROR CODES and CONDITIONS</h2>
<div class="errors">
<table border="1" cellspacing="1" cellpadding="10" width="100%">
<tr>
<th>Routine</th>
<th>Message</th>
<th>Comment</th>
</tr>
<tr><!-- routine -->
<td valign="top">initialize_module</td>
<!-- message -->
<td valign="top">initial allocation failed for satellite
observation data</td>
<!-- comment -->
<td valign="top">Need to increase MAXrttovkey</td>
</tr>
<tr><!-- routine -->
<td valign="top">initialize_rttov_sensor_runtime</td>
<!-- message -->
<td valign="top">Module or sensor is not initialized</td>
<!-- comment -->
<td valign="top">Both the module and the sensor must be initialized
before calling this routine.</td>
</tr>
<tr><!-- routine -->
<td valign="top">get_visir_metadata</td>
<!-- message -->
<td valign="top">The key exceeds the size of the metadata arrays,
or the key is not a VIS/IR type</td>
<!-- comment -->
<td valign="top">The number of satellite observations exceeds the
array size allocated in the module. Check the input and/or increase
MAXrttovkey.</td>
</tr>
<tr><!-- routine -->
<td valign="top">get_mw_metadata</td>
<!-- message -->
<td valign="top">The key exceeds the size of the metadata arrays,
or the key is not a MW type</td>
<!-- comment -->
<td valign="top">The number of satellite observations exceeds the
array size allocated in the module. Check the input and/or increase
MAXrttovkey.</td>
</tr>
<tr><!-- routine -->
<td valign="top">read_rttov_metadata</td>
<!-- message -->
<td valign="top">bad value for RTTOV fields</td>
<!-- comment -->
<td valign="top">The format of the input obs_seq file is not
consistent.</td>
</tr>
<tr><!-- routine -->
<td valign="top">get_expected_radiance</td>
<!-- message -->
<td valign="top">Could not find the platform/satellite/sensor id
combination in the RTTOV sensor database file.</td>
<!-- comment -->
<td valign="top">An unknown RTTOV instrument ID was encountered.
Check the database and/or the observation metadata.</td>
</tr>
</table>
</div>
<h2>KNOWN BUGS</h2>
<p>none at this time</p>
<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->
<a name="FuturePlans" id="FuturePlans"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FUTURE PLANS</h2>
<p>none at this time</p>
<!--==================================================================-->
<!-- PrivateComponents                                                -->
<!--==================================================================-->
<a name="PrivateComponents" id="PrivateComponents"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>PRIVATE COMPONENTS</h2>
<table>
<tr>
<td><em class="call">use obs_def_rttov_mod, only :</em></td>
<td><a href="#initialize_module">initialize_module</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#initialize_rttov_sensor_runtime">initialize_rttov_sensor_runtime</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#initialize_rttov_sensor_runtime">initialize_rttov_sensor_runtime</a></td>
</tr>
</table>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="initialize_module" id="initialize_module"></a><br>
<div class="routine"><em class="call">call
initialize_module()</em></div>
<div class="indent1"><!-- Description -->
<p>Reads the namelist, allocates space for the auxiliary data
associated wtih satellite observations, initializes the constants
used in subsequent computations (possibly altered by values in the
namelist), and prints out the list of constants and the values in
use.</p>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="initialize_rttov_sensor_runtime" id=
"initialize_rttov_sensor_runtime"></a><br>
<div class="routine"><em class="call">call
initialize_rttov_sensor_runtime(sensor,ens_size,nlevels)</em>
<pre>
type(rttov_sensor_type), pointer    :: <em class="code">sensor</em>
integer,                 intent(in) :: <em class=
"code">ens_size</em>
integer,                 intent(in) :: <em class=
"code">nlevels</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Initialize a RTTOV sensor runtime. A rttov_sensor_type instance
contains information such as options and coefficients that are
initialized in a "lazy" fashion only when it will be used for the
first time.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">sensor  </em></td>
<td>The sensor type to be initialized</td>
</tr>
<tr>
<td valign="top"><em class="code">ens_size</em></td>
<td>The size of the ensemble</td>
</tr>
<tr>
<td valign="top"><em class="code">nlevels</em></td>
<td>The number of vertical levels in the atmosphere</td>
</tr>
</table>
</div>
<br>
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
