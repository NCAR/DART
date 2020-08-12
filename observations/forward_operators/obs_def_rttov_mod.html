<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN"
          "http://www.w3.org/TR/html4/strict.dtd">
<HTML>
<HEAD>
<TITLE>module obs_def_rttov_mod</TITLE>
<link rel="stylesheet" type="text/css" href="../../docs/html/doc.css" />
<link href="../../docs/images/dart.ico" rel="shortcut icon" />
</HEAD>
<BODY>
<A NAME="TOP"></A>

<H1>MODULE <em class=program>obs_def_rttov_mod</em></H1>

<table border=0 summary="" cellpadding=5>
<tr>
    <td valign=middle>
    <img src="../../docs/images/Dartboard7.png" alt="DART project logo" height=70 />
    </td>
    <td>Jump to <a href="../../docs/index.html">DART Documentation Main Index</a></td>
</tr>
</table>

<A HREF="#Interface">INTERFACES</A> /
<A HREF="#Namelist">NAMELIST</A> /
<A HREF="#FilesUsed">FILES</A> /
<A HREF="#References">REFERENCES</A> /
<A HREF="#Errors">ERRORS</A> /
<A HREF="#FuturePlans">PLANS</A> /
<A HREF="#PrivateComponents">PRIVATE COMPONENTS</A> /
<A HREF="#Legalese">TERMS OF USE</A>

<H2>Overview</H2>

<P>
DART RTTOV observation module, including the observation operators for the
two primary RTTOV-observation types -- visible/infrared radiances and 
microwave radiances/brightness temperatures. 
<br>
<br>
This module acts as a pass-through for RTTOV version 12.3. For more information,
see <a href="https://www.nwpsaf.eu/site/software/rttov/documentation/">the RTTOV site</a>.
<br>
<br>
DART supports both RTTOV-direct for visible/infrared/microwave as well as
RTTOV-scatt for microwave computations. The code, in principle, supports all
features of version 12.3 as a pass-through from the model to RTTOV, includes
aerosols, trace gases, clouds, and atmospheric variables. The code also includes
directly specifying scattering properties.
<br>
<br>
However, a model may not have all of the variables necessary for these functions
depending on your model's setup.
<br>
<br>
For example, DART can use any of the RTTOV clw or ice schemes, but the WRF model
is not directly compatible with the IR default cloud classification of
marine/continental stratus/cumulus clean/dirty. We also offer a simple
classification based on maximum vertical velocity in the column and land type,
but due to lack of aerosol information, WRF/DART cannot differentiate between
clean and dirty cumulus. This may have some impact on the forward calculations -
but in experience the difference in cloud phase (ice versus water) makes a much
larger difference. 
<br>
<br>
Trace gases and aerosols may be important for actual observation system
experiments using visible/infrared; this may depend on the precise frequencies
you wish to use. Although a model may not have the necessary inputs by itself,
although the defaults in RTTOV based on climatology can be used. The impact on
the quality of the results should be investigated.
<br>
<br>
Known issues:
<ul>
<li>DART does not yet provide any type of bias correction
<li>Cross-channel error correlations are not yet supported. A principal
    component approach has been discussed. For now, the best bet is to use a
    subset of channels that are nearly independent of one another.
<li> Vertical localization will need to be tuned. Turning off vertical localization 
    may work well if you have a large number of ensemble members. Using the maximum peak 
    of the weighting function or the cloud-top may be appropriate. There are also other 
    potential approaches being investigated.
</ul>
<br>
<br>
Author and Contact information: 
</P>
<UL>
<LI>DART Code:  Jeff Steward, jsteward at ucar.edu</LI>
<LI>Original DART/RTTOV work:  Nancy Collins, Johnny Hendricks</LI>
</UL>

<H3>Backward compatibility note:</H3>

<!--==================================================================-->

<A NAME="OtherModulesUsed"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>OTHER MODULES USED</H2>
<PRE>
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
</PRE>

<!--==================================================================-->
<!-- Declare all public entities ...                                  -->
<!-- duplicate public routines template as many times as necessary    -->
<!-- make sure you replace all yyyroutine?? strings                   -->
<!--==================================================================-->
<!--Note to authors. The first row of the table is different.         -->
<!--==================================================================-->

<A NAME="Interface"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>PUBLIC INTERFACES</H2>
<TABLE>
<TR><TD><em class=call>use obs_def_rttov_mod, only : </em></TD>
                   <TD><A HREF="#set_visir_metadata">set_visir_metadata</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#set_mw_metadata">set_mw_metadata</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#get_expected_radiance">get_expected_radiance</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#get_rttov_option_logical">get_rttov_option_logical</A></TD></TR>
</TABLE>

<P>
Namelist interface
<A HREF="#Namelist"> <em class=code>&amp;obs_def_rttov_mod_nml</em> </A>
is read from file <em class=file>input.nml</em>.
</P>

<P>
   A note about documentation style.
   Optional arguments are enclosed in brackets
   <em class=optionalcode>[like this]</em>.
</P>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="set_visir_metadata"></A>
<br>
<div class=routine>
<em class=call>call set_visir_metadata(key, sat_az, sat_ze, sun_az, sun_ze, &
   platform_id, sat_id, sensor_id, channel, specularity)</em>
<pre>
integer,  intent(out) :: <em class=code>key</em>
real(r8), intent(in)  :: <em class=code>sat_az</em>
real(r8), intent(in)  :: <em class=code>sat_ze</em>
real(r8), intent(in)  :: <em class=code>sun_az</em>
real(r8), intent(in)  :: <em class=code>sun_ze</em>
integer,  intent(in)  :: <em class=code>platform_id, sat_id, sensor_id, channel</em>
real(r8), intent(in)  :: <em class=code>specularity</em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
Visible / infrared observations have several auxillary metadata variables.
Other than the key, which is standard DART fare, the RTTOV satellite azimuth 
and satellite zenith angle must be specified. See the RTTOV user guide for 
more information (in particular, see figure 4). If the <em class=code>addsolar</em> 
namelist value is set to true, then the solar azimuth and solar zenith angles must
be specified - again see the RTTOV user guide. In addition to the platform/satellite/
sensor ID numbers, which are the RTTOV unique identifiers, the channel specifies
the chanenl number in the RTTOV coefficient file. Finally, if 
<em class=code>do_lambertian</em> is true, specularity must be specified here. 
Again, see the RTTOV user guide for more information.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>
<TR><TD valign=top><em class=code>key&nbsp;&nbsp;</em></TD>
    <TD>The DART observation key.</TD></TR>
<TR><TD valign=top><em class=code>sat_az</em></TD>
    <TD>The satellite azimuth angle.</TD></TR>
<TR><TD valign=top><em class=code>sat_ze</em></TD>
    <TD>The satellite zenith angle.</TD></TR>
<TR><TD valign=top><em class=code>sun_az</em></TD>
    <TD>The solar azimuth angle. Only relevant if addsolar is true.</TD></TR>
<TR><TD valign=top><em class=code>sun_ze</em></TD>
    <TD>The solar zenith angle. Only relevant if addsolar is true.</TD></TR>
<TR><TD valign=top><em class=code>platform_id</em></TD>
    <TD>The RTTOV platform ID.</TD></TR>
<TR><TD valign=top><em class=code>sat_id</em></TD>
    <TD>The RTTOV satellite ID.</TD></TR>
<TR><TD valign=top><em class=code>sensor_id</em></TD>
    <TD>The RTTOV sensor ID.</TD></TR>
<TR><TD valign=top><em class=code>channel</em></TD>
    <TD>The RTTOV channel number.</TD></TR>
<TR><TD valign=top><em class=code>specularity</em></TD>
    <TD>The surface specularity. Only relevant if do_lambertian is true.</TD></TR>
</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="set_mw_metadata"></A>
<br>
<div class=routine>
<em class=call>call set_mw_metadata(key, sat_az, sat_ze, platform_id, sat_id, sensor_id, 
   channel, mag_field, cosbk, fastem_p1, fastem_p2, fastem_p3,         
   fastem_p4, fastem_p5)</em>
<pre>
integer,  intent(out) :: <em class=code>key</em>
real(r8), intent(in)  :: <em class=code>sat_az</em>
real(r8), intent(in)  :: <em class=code>sat_ze</em>
integer,  intent(in)  :: <em class=code>platform_id, sat_id, sensor_id, channel</em>
real(r8), intent(in)  :: <em class=code>mag_field</em>
real(r8), intent(in)  :: <em class=code>cosbk</em>
real(r8), intent(in)  :: <em class=code>fastem_p[1-5]</em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
Microwave observations have several auxillary metadata variables.
Other than the key, which is standard DART fare, the RTTOV satellite azimuth 
and satellite zenith angle must be specified. See the RTTOV user guide for 
more information (in particular, see figure 4). In addition to the platform/satellite/
sensor ID numbers, which are the RTTOV unique identifiers, the channel specifies
the chanenl number in the RTTOV coefficient file. In addition, if 
<em class=code>use_zeeman</em> is true, the magnetic field and cosine of the angle
between the magnetic field and angle of propagation must be specified. 
See the RTTOV user guide for more information. Finally, the fastem parameters for
land must be specified here. This may be difficult for observations to set, so 
default values (see table 21 in the RTTOV user guide) can be used until a better 
solution is devised. 
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>
<TR><TD valign=top><em class=code>key&nbsp;&nbsp;</em></TD>
    <TD>The DART observation key.</TD></TR>
<TR><TD valign=top><em class=code>sat_az</em></TD>
    <TD>The satellite azimuth angle.</TD></TR>
<TR><TD valign=top><em class=code>sat_ze</em></TD>
    <TD>The satellite zenith angle.</TD></TR>
<TR><TD valign=top><em class=code>platform_id</em></TD>
    <TD>The RTTOV platform ID.</TD></TR>
<TR><TD valign=top><em class=code>sat_id</em></TD>
    <TD>The RTTOV satellite ID.</TD></TR>
<TR><TD valign=top><em class=code>sensor_id</em></TD>
    <TD>The RTTOV sensor ID.</TD></TR>
<TR><TD valign=top><em class=code>channel</em></TD>
    <TD>The RTTOV channel number.</TD></TR>
<TR><TD valign=top><em class=code>mag_field</em></TD>
    <TD>The strength of the magnetic field. Only relevant if add_zeeman is true.</TD></TR>
<TR><TD valign=top><em class=code>cosbk</em></TD>
    <TD>The cosine of the angle between the magnetic field and direction of EM propagation. Only relevant if add_zeeman is true.</TD></TR>
<TR><TD valign=top><em class=code>fastem_p[1-5]</em></TD>
    <TD>The five parameters used for fastem land/sea ice emissivities. For ocean emissivities, an internal model is used based on the value of fastem_version.</TD></TR>
</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="get_expected_radiance"></A>
<br>
<div class=routine>
<em class=call> call get_expected_radiance(obs_kind_ind, state_handle, ens_size, location, 
         key, val, istatus) </em>
<pre>
integer,             intent(in)  :: <em class=code>obs_kind_ind</em>
type(ensemble_type), intent(in)  :: <em class=code>state_handle</em>
integer,             intent(in)  :: <em class=code>ens_size</em>
type(location_type), intent(in)  :: <em class=code>location</em>
integer,             intent(in)  :: <em class=code>key</em>
real(r8),            intent(out) :: <em class=code>val(ens_size)</em>
integer,             intent(out) :: <em class=code>istatus(ens_size)</em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
Given a location and the state vector from one of the ensemble members,
compute the model-predicted satellite observation. This can be either 
in units of radiance (mW/cm-1/sr/sq.m) or a brightness temperature 
(in K), depending on if this is a visible/infrared observation or a
microwave observation.
</P>
<TABLE width=100% border=0 summary="" cellpadding=3>
<TR><TD valign=top><em class=code>obs_kind_ind&nbsp;&nbsp;</em></TD>
    <TD>The index of the observation kind; since many observation kinds are
        handled by this module, this can be used to determine precisely which
        observation kind is being used.</TD></TR>
<TR><TD valign=top><em class=code>state_handle&nbsp;&nbsp;</em></TD>
    <TD>The ensemble of model states to be used for the observation operator
        calculations.</TD></TR>
<TR><TD valign=top><em class=code>location</em></TD>
    <TD>Location of this observation</TD></TR>
<TR><TD valign=top><em class=code>key</em></TD>
    <TD>Unique identifier associated with this satellite observation</TD></TR>
<TR><TD valign=top><em class=code>val</em></TD>
    <TD>The returned observation in units of either radiance or brightness temperature.</TD>
<TR><TD valign=top><em class=code>istatus</em></TD>
    <TD>Returned integer status code describing problems with applying
        forward operator.  0 is a good value; any positive value
        indicates an error; negative values are reserved for
        internal DART use only.</TD></TR>
</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="get_rttov_option_logical"></A>
<br>
<div class=routine>
<em class=call>p = get_rttov_option_logical(field_name)</em>
<pre>
character(len=*),           intent(in)  :: <em class=code>field_name</em>
logical,                    result      :: <em class=code>p</em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
Return the logical value of the RTTOV parameter associated with
the field_name.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>
<TR><TD valign=top><em class=code>field_name&nbsp;&nbsp;</em></TD>
    <TD>The name of the RTTOV parameter from the namelist.</TD></TR>
<TR><TD valign=top><em class=code>p</em></TD>
    <TD>The logical return value associated with the parameter.</TD></TR>
</TABLE>

</div>
<br>

<!--==================================================================-->
<!--=================== DESCRIPTION OF A NAMELIST  ===================-->
<!--==================================================================-->

<A NAME="Namelist"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>NAMELIST</H2>
<P>
This namelist is read from the file <em class=file>input.nml</em>.
Namelists start with an ampersand
'&amp;' and terminate with a slash '/'.
Character strings that contain a '/' must be
enclosed in quotes to prevent them from 
prematurely terminating the namelist.
</P>

<div class=namelist>
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
</pre>
</div>

<br />
<br />

<div>
<TABLE border=0 cellpadding=10 width=100% summary='namelist description'>
<THEAD align=left>
<TR><TH> Item </TH>
    <TH> Type </TH>
    <TH> Description </TH> </TR>
</THEAD>

<TBODY valign=top>

<TR><TD>rttov_sensor_db_file</TD>
    <TD>character(len=512)</TD>
    <TD> The location of the RTTOV sensor database. The format for the database
        is a comma-separated file. The columns of the database are the DART
        observation-kind, the platform/satellite/sensor ID, the observation
        type, the coefficient file, and a comma-separated list of RTTOV channels
        to use for this observation type.
</TD></TR>1

<TR><TD>first_lvl_is_sfc</TD>
    <TD>logical</TD>
    <TD>Whether the first level of the model represents the surface (true) or
        the top of the atmosphere (false).
</TD></TR>

<TR><TD>mw_clear_sky_only</TD>
    <TD>logical</TD>
    <TD>If microwave calculations should be "clear-sky" only (although
        cloud-liquid water absorption/emission is considered; see the RTTOV user
        guide).
</TD></TR>

<TR><TD>interp_mode</TD>
    <TD>integer</TD>
    <TD>The interpolation mode (see the RTTOV user guide).
</TD></TR>

<TR><TD>do_checkinput</TD>
    <TD>logical</TD>
    <TD>Whether to check the input for reasonableness (see the RTTOV user guide).
</TD></TR>

<TR><TD>apply_reg_limits</TD>
    <TD>logical</TD>
    <TD>Whether to clamp the atmospheric values to the RTTOV bounds (see the
        RTTOV user guide).
</TD></TR>

<TR><TD>verbose</TD>
    <TD>logical</TD>
    <TD>Whether to output lots of additional output (see the RTTOV user guide).
</TD></TR>

<TR><TD>fix_hgpl</TD>
    <TD>logical</TD>
    <TD>Whether the surface pressure represents the surface or the 2 meter value
        (see the RTTOV user guide).
</TD></TR>

<TR><TD>do_lambertian</TD>
    <TD>logical</TD>
    <TD>Whether to include the effects of surface specularity (see the RTTOV
        user guide).
</TD></TR>

<TR><TD>lambertian_fixed_angle</TD>
    <TD>logical</TD>
    <TD>Whether to include a fixed angle for the lambertian effect (see the
        RTTOV user guide).
</TD></TR>

<TR><TD>rad_down_lin_tau</TD>
    <TD>logical</TD>
    <TD>Whether to use the linear-in-tau approximation (see the RTTOV user guide).
</TD></TR>

<TR><TD>use_q2m</TD>
    <TD>logical</TD>
    <TD>Whether to use 2m humidity information (see the RTTOV user guide). If
        true, the QTY_2M_SPECIFIC_HUMIDITY will be requested from the model.
</TD></TR>

<TR><TD>use_q2m</TD>
    <TD>logical</TD>
    <TD>Whether to use 2m humidity information (see the RTTOV user guide). If
        true, the QTY_2M_SPECIFIC_HUMIDITY will be requested from the model.
</TD></TR>

<TR><TD>use_uv10m</TD>
    <TD>logical</TD>
    <TD>Whether to use 10m wind speed information (see the RTTOV user guide). If
        true, the QTY_10M_U_WIND_COMPONENT and QTY_10M_V_WIND_COMPONENTS will be
        requested from the model.
</TD></TR>

<TR><TD>use_wfetch</TD>
    <TD>logical</TD>
    <TD>Whether to use wind fetch information (see the RTTOV user guide). If
        true, the QTY_WIND_FETCH will be requested from the model.
</TD></TR>

<TR><TD>use_water_type</TD>
    <TD>logical</TD>
    <TD>Whether to use water-type information (0 = fresh, 1 = ocean; see the
        RTTOV user guide). If true, the QTY_WATER_TYPE will be requested from
        the model.
</TD></TR>

<TR><TD>addrefrac</TD>
    <TD>logical</TD>
    <TD>Whether to enable atmospheric refraction (see the RTTOV user guide). 
</TD></TR>

<TR><TD>plane_parallel</TD>
    <TD>logical</TD>
    <TD>Whether to treat the atmosphere as plane parallel (see the RTTOV user
        guide). 
</TD></TR>

<TR><TD>use_salinity</TD>
    <TD>logical</TD>
    <TD>Whether to use salinity (see the RTTOV user guide). If true, the
        QTY_SALINITY will be requested from the model.
</TD></TR>

<TR><TD>apply_band_correction</TD>
    <TD>logical</TD>
    <TD>Whether to apply band correction from the coefficient field for
        microwave data (see the RTTOV user guide). 
</TD></TR>

<TR><TD>cfrac_data</TD>
    <TD>logical</TD>
    <TD>Whether to use the cloud fraction from 0 to 1 (see the RTTOV user
        guide). If true, the QTY_CLOUD_FRACTION will be requested from the
        model.
</TD></TR>

<TR><TD>clw_data</TD>
    <TD>logical</TD>
    <TD>Whether to use cloud-liquid water data (see the RTTOV user guide). If
        true, the QTY_CLOUDWATER_MIXING_RATIO will be requested from the model.
</TD></TR>

<TR><TD>rain_data</TD>
    <TD>logical</TD>
    <TD>Whether to use precipitating water data (see the RTTOV user guide). If
        true, the QTY_RAINWATER_MIXING_RATIO will be requested from the model.
</TD></TR>

<TR><TD>ciw_data</TD>
    <TD>logical</TD>
    <TD>Whether to use non-precipiting ice information (see the RTTOV user
        guide). If true, the QTY_ICE_MIXING_RATIO will be requested from the
        model.
</TD></TR>

<TR><TD>snow_data</TD>
    <TD>logical</TD>
    <TD>Whether to use precipitating fluffy ice (see the RTTOV user guide). If
        true, the QTY_SNOW_MIXING_RATIO will be requested from the model.
</TD></TR>

<TR><TD>graupel_data</TD>
    <TD>logical</TD>
    <TD>Whether to use precipting small, hard ice (see the RTTOV user guide). If
        true, the QTY_GRAUPEL_MIXING_RATIO will be requested from the model.
</TD></TR>

<TR><TD>hail_data</TD>
    <TD>logical</TD>
    <TD>Whether to use precipitating large, hard ice (see the RTTOV user guide).
        If true, the QTY_HAIL_MIXING_RATIO will be requested from the model.
</TD></TR>

<TR><TD>w_data</TD>
    <TD>logical</TD>
    <TD>Whether to use vertical velocity information. This will be used to
        crudely classify if a cloud is cumulus or stratiform for the purpose of
        visible/infrared calculations. If true, the QTY_VERTICAL_VELOCITY will
        be requested from the model.
</TD></TR>

<TR><TD>clw_scheme</TD>
    <TD>integer</TD>
    <TD>The clw_scheme to use (see the RTTOV user guide).
</TD></TR>

<TR><TD>clw_cloud_top</TD>
    <TD>real(r8)</TD>
    <TD>Lower hPa limit for clw calculations (see the RTTOV user guide). 
</TD></TR>

<TR><TD>fastem_version</TD>
    <TD>integer</TD>
    <TD>Which FASTEM version to use (see the RTTOV user guide). 
</TD></TR>

<TR><TD>supply_foam_fraction</TD>
    <TD>logical</TD>
    <TD>Whether to use sea-surface foam fraction (see the RTTOV user guide). If
        true, the QTY_FOAM_FRAC will be requested from the model.
</TD></TR>

<TR><TD>use_totalice</TD>
    <TD>logical</TD>
    <TD>Whether to use totalice instead of precip/non-precip ice for microwave
        (see the RTTOV user guide). 
</TD></TR>

<TR><TD>use_zeeman</TD>
    <TD>logical</TD>
    <TD>Whether to use the Zeeman effect (see the RTTOV user guide). If true,
        the magnetic field and cosine of bk will be used from the observation
        metadata.
</TD></TR>

<TR><TD>cc_threshold</TD>
    <TD>real(r8)</TD>
    <TD>Cloud-fraction value to treat as clear-sky (see the RTTOV user guide). 
</TD></TR>

<TR><TD>ozone_data</TD>
    <TD>logical</TD>
    <TD>Whether to use ozone (O3) profiles (see the RTTOV user guide). If true,
        the QTY_O3 will be requested from the model.
</TD></TR>

<TR><TD>co2_data</TD>
    <TD>logical</TD>
    <TD>Whether to use carbon dioxide (CO2) profiles (see the RTTOV user guide).
        If true, the QTY_CO2 will be requested from the model.
</TD></TR>

<TR><TD>n2o_data</TD>
    <TD>logical</TD>
    <TD>Whether to use nitrous oxide (N2O) profiles (see the RTTOV user guide).
        If true, the QTY_N2O will be requested from the model.
</TD></TR>

<TR><TD>co_data</TD>
    <TD>logical</TD>
    <TD>Whether to use carbon monoxide (CO) profiles (see the RTTOV user guide).
        If true, the QTY_CO will be requested from the model.
</TD></TR>

<TR><TD>ch4_data</TD>
    <TD>logical</TD>
    <TD>Whether to use methane (CH4) profiles (see the RTTOV user guide). If
        true, the QTY_CH4 will be requested from the model.
</TD></TR>

<TR><TD>so2_data</TD>
    <TD>logical</TD>
    <TD>Whether to use sulfur dioxide (SO2) (see the RTTOV user guide). If true,
        the QTY_SO2 will be requested from the model.
</TD></TR>

<TR><TD>addsolar</TD>
    <TD>logical</TD>
    <TD>Whether to use solar angles (see the RTTOV user guide). If true, the
        sun_ze and sun_az from the observation metadata will be used for
        visible/infrared.
</TD></TR>

<TR><TD>rayleigh_single_scatt</TD>
    <TD>logical</TD>
    <TD>Whether to use only single scattering for Rayleigh scattering for
        visible calculations (see the RTTOV user guide). 
</TD></TR>

<TR><TD>do_nlte_correction</TD>
    <TD>logical</TD>
    <TD>Whether to include non-LTE bias correction for HI-RES sounder (see the
        RTTOV user guide). 
</TD></TR>

<TR><TD>solar_sea_brdf_model</TD>
    <TD>integer</TD>
    <TD>The solar sea BRDF model to use (see the RTTOV user guide). 
</TD></TR>

<TR><TD>ir_sea_emis_model</TD>
    <TD>logical</TD>
    <TD>The infrared sea emissivity model to use (see the RTTOV user guide). 
</TD></TR>

<TR><TD>use_sfc_snow_frac</TD>
    <TD>logical</TD>
    <TD>Whether to use the surface snow fraction (see the RTTOV user guide). If
        true, the QTY_SNOWCOVER_FRAC will be requested from the model.
</TD></TR>

<TR><TD>add_aerosl</TD>
    <TD>logical</TD>
    <TD>Whether to use aerosols (see the RTTOV user guide). 
</TD></TR>

<TR><TD>aerosl_type</TD>
    <TD>integer</TD>
    <TD>Whether to use OPAC or CAMS aerosols (see the RTTOV user guide). 
</TD></TR>

<TR><TD>add_clouds</TD>
    <TD>logical</TD>
    <TD>Whether to enable cloud scattering for visible/infrared (see the RTTOV user guide).
</TD></TR>

<TR><TD>ice_scheme</TD>
    <TD>integer</TD>
    <TD>The ice scheme to use (see the RTTOV user guide). 
</TD></TR>

<TR><TD>use_icede</TD>
    <TD>logical</TD>
    <TD>Whether to use the ice effective diameter for visible/infrared (see the
        RTTOV user guide). If true, the QTY_CLOUD_ICE_DE will be requested from
        the model.
</TD></TR>

<TR><TD>idg_scheme</TD>
    <TD>integer</TD>
    <TD>The ice water effective diameter scheme to use (see the RTTOV user guide). 
</TD></TR>

<TR><TD>user_aer_opt_param</TD>
    <TD>logical</TD>
    <TD>Whether to directly specify aerosol scattering properties (see the RTTOV
        user guide). Not yet supported.
</TD></TR>

<TR><TD>user_cld_opt_param</TD>
    <TD>logical</TD>
    <TD>Whether to directly specify cloud scattering properties (see the RTTOV
        user guide). Not yet supported.
</TD></TR>

<TR><TD>grid_box_avg_cloud</TD>
    <TD>logical</TD>
    <TD>Whether to cloud concentrations are grid box averages (see the RTTOV
        user guide). 
</TD></TR>

<TR><TD>cldstr_threshold</TD>
    <TD>real(r8)</TD>
    <TD>Threshold for cloud stream weights for scattering (see the RTTOV user
        guide).
</TD></TR>

<TR><TD>cldstr_simple</TD>
    <TD>logical</TD>
    <TD>Whether to use one clear and one cloudy column (see the RTTOV user
        guide). 
</TD></TR>

<TR><TD>cldstr_low_cloud_top</TD>
    <TD>real(r8)</TD>
    <TD>Cloud fraction maximum in layers from the top of the atmosphere down to
        the specified hPa (see the RTTOV user guide). 
</TD></TR>

<TR><TD>ir_scatt_model</TD>
    <TD>integer</TD>
    <TD>Which infrared scattering method to use (see the RTTOV user guide). 
</TD></TR>

<TR><TD>vis_scatt_model</TD>
    <TD>integer</TD>
    <TD>Which visible scattering method to use (see the RTTOV user guide).
</TD></TR>

<TR><TD>dom_nstreams</TD>
    <TD>integer</TD>
    <TD>The number of streams to use with DOM (see the RTTOV user guide). 
</TD></TR>

<TR><TD>dom_accuracy</TD>
    <TD>real(r8)</TD>
    <TD>The convergence criteria for DOM (see the RTTOV user guide). 
</TD></TR>

<TR><TD>dom_opdep_threshold</TD>
    <TD>real(r8)</TD>
    <TD>Ignore layers below this optical depth (see the RTTOV user guide). 
</TD></TR>

<TR><TD>addpc</TD>
    <TD>logical</TD>
    <TD>Whether to do principal component calculations (see the RTTOV user
        guide).
</TD></TR>

<TR><TD>npcscores</TD>
    <TD>integer</TD>
    <TD>Number of principal components to use for addpc (see the RTTOV user
        guide). 
</TD></TR>

<TR><TD>addradrec</TD>
    <TD>logical</TD>
    <TD>Reconstruct the radiances using addpc (see the RTTOV user guide). 
</TD></TR>

<TR><TD>ipcreg</TD>
    <TD>integer</TD>
    <TD>Number of predictors to use with addpc (see the RTTOV user guide).
</TD></TR>

<TR><TD>use_htfrtc</TD>
    <TD>logical</TD>
    <TD>Whether to use HTFRTC (see the RTTOV user guide). 
</TD></TR>

<TR><TD>htfrtc_n_pc</TD>
    <TD>integer</TD>
    <TD>Number of PCs to use with HTFRTC (see the RTTOV user guide). 
</TD></TR>

<TR><TD>htfrtc_simple_cloud</TD>
    <TD>logical</TD>
    <TD>Whether to use simple cloud scattering with htfrtc (see the RTTOV user
        guide). 
</TD></TR>

<TR><TD>htfrtc_overcast</TD>
    <TD>logical</TD>
    <TD>Whether to calculate overcast radiances with HTFRTC (see the RTTOV user
        guide).
</TD></TR>

</TBODY> 
</TABLE>
</div>

<br />
<br />

<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->

<A NAME="FilesUsed"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>FILES</H2>
<UL><LI>A DART observation sequence file containing Radar obs.
</UL>

<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->

<A NAME="References"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>REFERENCES</H2>
<UL>
<LI> <a href="https://www.nwpsaf.eu/site/software/rttov/documentation/">RTTOV user guide</a>
</UL>

<!--==================================================================-->
<!-- Describe all the error conditions and codes.                     -->
<!--==================================================================-->

<A NAME="Errors"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>ERROR CODES and CONDITIONS</H2>
<div class=errors>
<TABLE border=1 cellspacing=1 cellpadding=10 width=100%>
<TR><TH>Routine</TH><TH>Message</TH><TH>Comment</TH></TR>

<TR><!-- routine --><TD VALIGN=top>initialize_module</TD>
    <!-- message --><TD VALIGN=top>initial allocation failed for satellite
        observation data</TD>
    <!-- comment --><TD VALIGN=top>Need to increase MAXrttovkey</TD>
</TR>
<TR><!-- routine --><TD VALIGN=top>initialize_rttov_sensor_runtime</TD>
    <!-- message --><TD VALIGN=top>Module or sensor is not initialized</TD>
    <!-- comment --><TD VALIGN=top>Both the module and the sensor must be
        initialized before calling this routine.</TD>
</TR>
<TR><!-- routine --><TD VALIGN=top>get_visir_metadata</TD>
    <!-- message --><TD VALIGN=top>The key exceeds the size of the metadata
        arrays, or the key is not a VIS/IR type</TD>
    <!-- comment --><TD VALIGN=top>The number of satellite observations exceeds 
                                   the array size allocated in the module. 
                                   Check the input and/or increase MAXrttovkey.</TD>
</TR>
<TR><!-- routine --><TD VALIGN=top>get_mw_metadata</TD>
    <!-- message --><TD VALIGN=top>The key exceeds the size of the metadata
        arrays, or the key is not a MW type</TD>
    <!-- comment --><TD VALIGN=top>The number of satellite observations exceeds 
                                   the array size allocated in the module. 
                                   Check the input and/or increase MAXrttovkey.</TD>
</TR>
<TR><!-- routine --><TD VALIGN=top>read_rttov_metadata</TD>
    <!-- message --><TD VALIGN=top>bad value for RTTOV fields</TD>
    <!-- comment --><TD VALIGN=top>The format of the input obs_seq file is not consistent.</TD>
</TR>
<TR><!-- routine --><TD VALIGN=top>get_expected_radiance</TD>
    <!-- message --><TD VALIGN=top>Could not find the platform/satellite/sensor
        id combination in the RTTOV sensor database file.</TD>
    <!-- comment --><TD VALIGN=top>An unknown RTTOV instrument ID was
        encountered. Check the database and/or the observation metadata.</TD>
</TR>
</TABLE>
</div>

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
none at this time
</P>

<!--==================================================================-->
<!-- PrivateComponents                                                -->
<!--==================================================================-->

<A NAME="PrivateComponents"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>PRIVATE COMPONENTS</H2>

<TABLE>
<TR><TD><em class=call>use obs_def_rttov_mod, only : </em></TD>
                   <TD><A HREF="#initialize_module">initialize_module</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#initialize_rttov_sensor_runtime">initialize_rttov_sensor_runtime</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#initialize_rttov_sensor_runtime">initialize_rttov_sensor_runtime</A></TD></TR>
</TABLE>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="initialize_module"></A>
<br>
<div class=routine>
<em class=call>call initialize_module()</em>
</div>

<div class=indent1>
<!-- Description -->

<P>
Reads the namelist, allocates space for the auxiliary data associated
wtih satellite observations, initializes the constants used in
subsequent computations (possibly altered by values in the namelist),
and prints out the list of constants and the values in use.
</P>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="initialize_rttov_sensor_runtime"></A>
<br>
<div class=routine>
<em class=call>call initialize_rttov_sensor_runtime(sensor,ens_size,nlevels)</em>
<pre>
type(rttov_sensor_type), pointer    :: <em class=code>sensor</em>
integer,                 intent(in) :: <em class=code>ens_size</em>
integer,                 intent(in) :: <em class=code>nlevels</em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<p>
Initialize a RTTOV sensor runtime. A rttov_sensor_type instance contains 
information such as options and coefficients that are initialized in a "lazy"
fashion only when it will be used for the first time.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>
<TR><TD valign=top><em class=code>sensor&nbsp;&nbsp;</em></TD>
    <TD>The sensor type to be initialized</TD></TR>
<TR><TD valign=top><em class=code>ens_size</em></TD>
    <TD>The size of the ensemble</TD></TR>
<TR><TD valign=top><em class=code>nlevels</em></TD>
    <TD>The number of vertical levels in the atmosphere</TD></TR>
</TABLE>

</div>
<br>

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
