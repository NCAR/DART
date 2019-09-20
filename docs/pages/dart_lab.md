---
title: DART_LAB
layout: default
---

# DART_LAB Tutorial

## Overview
DART_LAB is a MATLAB®-based tutorial to demonstrate the principles of 
ensemble data assimilation.
DART_LAB consists of PDF tutorial materials and MATLAB® exercises.
See below for links to the PDF files and a list of the
corresponding MATLAB scripts.
The DART_LAB tutorial begins at a more introductory level than the materials in
the tutorial directory, and includes hands-on exercises at several points.
<!-- In a workshop setting, the full tutorial materials and exercises took
about 1.5 days to complete. -->

-----

<table>
<tbody>
<tr><td>Section 1</td><td><a href="../DART_LAB/presentation/DART_LAB_Section_01.pdf">The basics in 1D.</a></td></tr>
<tr><td>Section 2</td><td><a href="../DART_LAB/presentation/DART_LAB_Section_02.pdf">How should observations of a state variable impact an unobserved state variable? Multivariate assimilation.</a></td></tr>
<tr><td>Section 3</td><td><a href="../DART_LAB/presentation/DART_LAB_Section_03.pdf">Sampling error and localization.</a></td></tr>
<tr><td>Section 4</td><td><a href="../DART_LAB/presentation/DART_LAB_Section_04.pdf">The Ensemble Kalman Filter (Perturbed Observations).</a></td></tr>
<tr><td>Section 5</td><td><a href="../DART_LAB/presentation/DART_LAB_Section_05.pdf">Adaptive Inflation.</a></td></tr>
</tbody>
</table>

<span id="Matlab"></span>

-----

## MATLAB® Hands-On Exercises

In the `matlab` subdirectory are a set of MATLAB scripts and GUI
(graphical user interface) programs which are exercises that go with the
tutorial. Each is interactive with settings that can be changed and
rerun to explore various options. A valid
[MATLAB](http:/v0.0.2/www.mathworks.com/products/matlab/)
license is needed to run these scripts.

The exercises use the following functions:

| function            | description |
| ---                 | :---         |
| `gaussian_product`  | graphical representation of the product of two gaussians |
| `oned_ensemble`     | explore the details of ensemble data assimilation for a scalar |
| `oned_model`        | simple ensemble data assimilation example |
| `oned_model_inf`    | simple ensemble data assimilation example *with inflation* |
| `run_lorenz_63`     | ensemble DA with the 3-variable Lorenz '63 dynamical model - the "butterfly" model |
| `run_lorenz_96`     | ensemble DA with the 40-variable Lorenz '96 dynamical model |
| `run_lorenz_96_inf` | ensemble DA with the 40-variable Lorenz '96 dynamical model *with inflation* |
| `twod_ensemble`     | demonstrates the impact of observations on unobserved state variables |

To run these, cd into the `DART_LAB/matlab` directory, start matlab, and
type the names at the prompt.

-----
