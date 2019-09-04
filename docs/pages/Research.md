---
title: "Research"
layout: default
menu: [main]
---

## **R**esearch

DART research is broadly categorized along three avenues: one is the
research toward data assimilation algorithmic and computational
efficiency, another is toward implementing and exploiting the
information in observations with one model or another, and another is to
learn more about the behavior of an individual model - perhaps by
looking at systematic features of the
increments.

### If you would like to add your research highlight using DART, please let us know\!

<span id="with_dart" class="anchor"></span>  

[[top](#)]

-----

### Research Projects involving DART

  - [Model Performance](#modelperformance)
  - [WACCM](#waccm)
  - [OpenGGCM](#space)
  - [Chemical Transport](#chemistry)
  - [\[novel\] observations like GPS RO](#gpsro)
  - [Sensitivity Analyses](#sensitivity)
  - [Carbon Monoxide](#cmaqco)
  - [Data Assimilation on Mars](#mars)
  - [Inertio-Gravity waves](#gravity)
  - [Boundary Layer projects](#pbl1d)
  - [Lorenz '96](#loworder)
  - \[\[Damped\] Adaptive\] Inflation Algorithms
  - Radar Reflectivity
  - Ocean Modeling
  - Tropical Cyclones

<span id="modelperformance" class="anchor"></span>

[[top](#)]

-----

### Model Performance

<table>
<tbody>
<tr class="odd">
<td><img src="../images/science_nuggets/DARTspaghettiSquare_full.png" alt="DART/CAM spaghetti plot" /></td>
<td>Ensemble Data assimilation can provide qualitative and quantitative uncertainty for quantities of interest to weather forecasters. <a href="Research/CAM_Raeder/index.html">[link to more information]</a><br />
<br />
Kevin Raeder, raeder@ucar.edu</td>
</tr>
</tbody>
</table>

<span id="waccm" class="anchor"></span> 

[[top](#)]

-----

### Whole Atmosphere Community Climate Model (WACCM)

<table>
<colgroup>
<col style="width: 50%" />
<col style="width: 50%" />
</colgroup>
<tbody>
<tr class="odd">
<td><img src="../images/science_nuggets/waccm_dart_ssw_t.png" alt="WACCM graphic" width="1500" height="1500"/>  </td>
<td><h5 id="this-is-an-excerpt-from-applying-the-science-and-technology-of-data-assimilation-by-brian-bevirt-07112017-as-part-of-a-cisl-news-series-describing-the-many-ways-cisl-improves-modeling-beyond-providing-supercomputing-systems-and-facilities.">This is an excerpt from <a href="https://www2.cisl.ucar.edu/news/applying-the-science-and-technology-of-data-assimilation">Applying the science and technology of data assimilation</a> by <em>Brian Bevirt</em> 07/11/2017 as part of a <em>CISL News</em> series describing the many ways CISL improves modeling beyond providing supercomputing systems and facilities.</h5>
These plots show measured and modeled zonal mean temperatures between 70N and 90N during the January 2009 sudden warming of the stratosphere. The bottom plot shows the observed temperatures (in degrees Kelvin, see legend at right), the center plot shows how this state of the atmosphere was simulated by the specified-dynamics version of the WACCM model, and the top plot shows WACCM's improved result after using DART to assimilate middle atmosphere observations.<br />
<br />
The key point in this figure is that WACCM+DART captures both the stratosphere warming and mesosphere cooling that are seen in the observations. Also seen in the specified-dynamics version of WACCM, the elevated stratopause that forms at high altitudes around day 30 descends too fast compared to the observations. The elevated stratopause is maintained at a high altitude in the WACCM+DART simulation. This has implications for the descent of species from the mesosphere into the stratosphere. Accurate representation of the mesosphere dynamics is important for the ionosphere variability during sudden stratosphere warming events. (Figure courtesy of Nick Pedatella, HAO)<br />
<br />
<a href="Research/WACCM_Pedatella/index.html">[link to more information]</a><br />
<br />
Nick Pedatella, nickp@ucar.edu Hanli Liu, liuh@ucar.edu Jing Liu, jingliu@ucar.edu</td>
</tr>
</tbody>
</table>

<span id="space" class="anchor"></span>

[[top](#)]

-----

### Open Geospace General Circulation Model (thermosphere/ionosphere/magnetosphere)

<table>
<tbody>
<tr class="odd">
<td><img src="../images/science_nuggets/OpenGGCM_graphic.png" alt="OpenGGCM graphic" width="1500" height="1500"/></td>
<td>The primary goal of this project is to combine the OpenGGCM (Open Geospace General Circulation Model) with the NCAR Data Assimilation Research Testbed (DART), which implements an Ensemble Kalman Filter (EnKF) and will enable the sequential assimilation of ionosphere, thermosphere, and magetosphere data.<br />
<br />
We will also optimize model parameters by including them into the state vector. This will improve model accuracy even when no data are assimilated. <a href="Research/OpenGGCM_Raeder/index.html">[link to more information]</a><br />
<br />
Jimmy Raeder, J.Raeder@unh.edu</td>
</tr>
</tbody>
</table>

<span id="chemistry" class="anchor"></span>

[[top](#)]

-----

### Chemical Data Assimilation

<table>
<tbody>
<tr class="odd">
<td><img src="../images/science_nuggets/ChemicalWeather_NHprojection_full.png" alt="CAM-Chem/DART CO Column" width="1500" height="1500"/></td>
<td>We are currently applying an ensemble-based chemical data assimilation system, consisting of regional to global chemical transport models (CAM-Chem, WRF-Chem) in conjunction with DART, for a joint assimilation of meteorological observations and satellite-derived CO measurements from MOPITT and aerosol optical depth (AOD) measurements from MODIS. The chemical data assimilation system has been recently used for near-real time chemical forecasting (see <a href="https://espo.nasa.gov/arctas/" class="uri">https://espo.nasa.gov/arctas/</a>) to support flight planning during the NASA Arctic Research of the Composition of the Troposphere from Aircraft and Satellites (ARCTAS). <a href="Research/Chem_Arellano/index.html">[link to more information]</a><br />
<br />
Ave Arellano, arellano@ucar.edu</td>
</tr>
</tbody>
</table>

<span id="gpsro" class="anchor"></span>

[[top](#)]

-----

### GPS RO Observations and Tropical Cyclone Forecasting

<table>
<tbody>
<tr class="odd">
<td><img src="../images/science_nuggets/GPS_LEO_full.png" alt="GPS RO geometry schematic" width="1500" height="1500"/></td>
<td>Profiles of atmospheric quantities deduced from GPS Radio Occultation data are available in otherwise data-sparse regions and provide information used to forecast the behavior of tropical cyclones. The <a href="http://www.cosmic.ucar.edu/">COSMIC/FORMOSAT-3</a> mission has been providing about 2000 data profiles per day since September 2007. <a href="Research/GPS_Liu/index.html">[link to more information]</a><br />
<br />
Hui Liu, hliu@ucar.edu</td>
</tr>
</tbody>
</table>

<span id="sensitivity" class="anchor"></span>

[[top](#)]

-----

### Sensitivity Analyses

<table>
<tbody>
<tr class="odd">
<td><img src="../images/science_nuggets/2005082500_f048_DLMU_sens_full.png" alt="DART/WRF sensitivity plot" width="1500" height="1500"/></td>
<td>Forecast sensitivity analysis provides an objective means of evaluating how initial condition errors affect a forecast and where to gather additional observations to reduce forecast errors. Most sensitivity studies use the adjoint of a linearized forecast model to determine the gradient of a forecast metric with respect to the initial conditions. Adjoints suffer from a number of difficulties including coding, linearity assumptions, and moist processes. Ensemble-based sensitivity analysis provides an attractive alternative to adjoint-based methods because it combines data assimilation and sensitivity analysis in a consistent manner. This image illustrates the effect of zonal winds aloft on the position of Hurricane Katrina. <a href="Research/Katrina_Torn/index.html">[link to more information]</a><br />
<br />
<a href="http://www.atmos.albany.edu/index.php?d=faculty.torn">Ryan Torn</a>, torn@atmos.albany.edu</td>
</tr>
</tbody>
</table>

<span id="cmaqco" class="anchor"></span>

[[top](#)]

-----

### Assimilation of CO

<table>
<tbody>
<tr class="odd">
<td><img src="../images/science_nuggets/CMAQ_teasergraphic_full.png" alt="CMAQ CO assimilation result" width="1500" height="1500" /></td>
<td>This project describes an integrated approach to modeling atmospheric chemistry with trace gas data assimilation. Specifically, we ran CMAQ from within DART to assimilate both synthetic and real observations of CO for the period of June 2001. <a href="Research/CMAQ_Zubrow/index.html">[link to more information]</a><br />
<br />
Alexis Zubrow, azubrow@unc.edu</td>
</tr>
</tbody>
</table>

<span id="mars" class="anchor"></span> 

[[top](#)]

-----

### Assimilation on MARS

<table>
<tbody>
<tr class="odd">
<td><img src="../images/science_nuggets/MarsWRF_teaser_full.png" alt="DART/MARS_Lawson graphic" width="1500" height="1500" /></td>
<td>The planetary atmospheres group at Caltech has produced a global and planetary version of NCAR's WRF (Weather Research and Forecasting) Model. We are using DART to attempt data assimilation within the Mars atmosphere using the Mars version of WRF, MarsWRF, as our GCM. <a href="Research/MARS_Lawson/index.html">[link to more information]</a><br />
<br />
<a href="http://www.gps.caltech.edu/people/wglawson/profile">Greg Lawson</a>, wglawson@gps.caltech.edu</td>
</tr>
</tbody>
</table>

<span id="gravity" class="anchor"></span>  

[[top](#)]

-----

### Inertio-Gravity waves

<table>
<tbody>
<tr class="odd">
<td><img src="../images/science_nuggets/Zagar_1_full.png" alt="DART/CAM normal modes" width="1500" height="1500"/></td>
<td>This project aims at quantifying the impact of various motion types in analysis and forecast fields by using normal modes. The DART/CAM is the main analysis system used in the project. The first question addressed is about how large part of the atmospheric energy is associated with the inertio-gravity motions, an important part of the global circulation primarily because of their role in the tropical system. <a href="Research/CAM_Zagar/index.html">[link to more information]</a><br />
<br />
<a href="http://www.fmf.uni-lj.si/~zagarn">Nedjeljka Zagar</a>, nedjeljka.zagar@fmf.uni-lj.si</td>
</tr>
</tbody>
</table>

<span id="pbl1d" class="anchor"></span> 

[[top](#)]

-----

### Planetary Boundary Layer

<table>
<tbody>
<tr class="odd">
<td><img src="../images/science_nuggets/JoshDorita_graphic_oneframe_full.png" alt="DART/CAM spaghetti plot" width="1500" height="1500"/></td>
<td>A long-term goal of this work is to find an efficient system for probabilistic planetary boundary layer (PBL) nowcasting that can be employed wherever surface observations are present. One approach showing promise is the use of a single column model (SCM) and ensemble filter data assimilation techniques. <a href="Research/PBL_JoshDorita/index.html">[link to more information]</a><br />
<br />
Dorita Rostkier-Edelstein, rostkier@ucar.edu<br />
Josh Hacker, hacker@ucar.edu</td>
</tr>
</tbody>
</table>

<span id="loworder" class="anchor"></span>  

[[top](#)]

-----

### The Lorenz '96 model

#### lorenz_96

Lorenz, E. N., and K. A. Emanuel, 1998:  
Optimal sites for supplementary weather observations: Simulations with a
small
model.  
*J. Atmos. Sci.*, **55**, 399-414.  
[10.1175/1520-0469(1998)055\<0399:OSFSWO\>2.0.CO;2](http://dx.doi.org/10.1175/1520-0469\(1998\)055%3C0399%3AOSFSWO%3E2.0.CO%3B2)

<table>
<colgroup>
<col style="width: 33%" />
<col style="width: 33%" />
<col style="width: 33%" />
</colgroup>
<tbody>
<tr class="odd">
<td><img src="../images/science_nuggets/L96BadRankHistogram_full.png" alt="bad rank histogram"/></td>
<td>The Lorenz '96 model is one of our favorite models. In our implementation, it is a 40-variable model that can be used to test inflation algorithms, the effects of localization schemes, the integrity of the DART installation itself, the state-space diagnostic routines; it is extensively used in the tutorial, <strong>and</strong> can even be run as a standalone executable to test the MPI support on a machine. <a href="Research/Lorenz96/index.html">[link to more information]</a><br />
<br />

Jeff Anderson, jla@ucar.edu, and<br />
Tim Hoar, thoar@ucar.edu</td>
<td><img src="../images/science_nuggets/L96GoodRankHistogram_full.png" alt="good rank histogram"/></td>
</tr>
</tbody>
</table>

<span id="suggestions" class="anchor"></span> 

[[top](#)]

-----

  

### Please [suggest ways for us to improve DART.](https://www2.cisl.ucar.edu/image/software/dart/suggestion-box)
