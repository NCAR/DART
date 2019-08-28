## **Finite Element Sea-ice Ocean Model (FESOM)**

FESOM is an unstructured mesh global ocean model using finite element methods to solve the
hydro-static primitive equations with the Boussinesq approximation (Danilov et al., 2004; Wang et
al., 2008). FESOM v1.4 is interfaced with DART by Aydoğdu et al. (2018a) using a regional
implementation in Turkish Straits System (Gürses et al. 2016, Aydoğdu et al. 2018).

There is a recent version of the model called the Finite-volumE Sea ice–Ocean Model (FESOM2,
Danilov et al. 2017). A version for coastal applications FESOM-C v.2 (Androsov et al., 2019) is also
presented recently.


## **FESOM/DART interface**
### model_mod.f90

A module called *fesom_modules* is provided to pass the information from FESOM to DART.
fesom_modules.f90 includes fortran routines adopted from FESOM v1.4 to read the mesh, set the
variables and dimensions of the arrays. *fesom_modules* should have access to *nod2d.out*, *nod3d.out*,
*elem2d.out*, *elem3d.out*, *aux3d.out*, *depth.out* and *m3d.ini* mesh files.

Forward operators use an interpolation using closing model node in the horizontal, given that the
application in Aydoğdu et al. (2018a) uses a very high-resolution mesh. In the vertical, a linear
interpolation is performed between two enclosing model layers. Interpolation in model_interpolate
routine can be improved, if needed.

### Shell scripts

Shell scripts are written in bash for lsf queuing system. They should be modified to work with
others such as slurm. FESOM executables are called externally detached from DART therefore no need
for an advance model.


| Script                       | Queue  | Definition    |
| ---------------------------- | ------:|---------------|
| **environment.load**         | serial | Includes environment variables, relevant directories, experiment specifications. This file is sourced by every other script below. |
| **ensemble.launch**          | serial | Main script which modifies ens_members.template.lsf and calling ens_members.${EXPINFO}.lsf. Keeps also an experiment specific brief information which should be modified before launching the scripts to keep track. |
| **ens_members.template.lsf** | serial | Calls and submits initialize.template, forward_model.template check_ensemble.template subsequently. |
| **initialize.template**      | serial | Called only once at the beginning of the experiment. Sets the experiment directory, copies initial ensemble, namelists. |
| **forward_model.template**   |parallel| Submits a job array for all ensemble members. |
| **check_ensemble.template**  | serial | Checks if the forwarding for all members is finished. If so, first calls filter.template and then calls finalize.template to conclude current assimilation cycle. |
| **filter.template**          |parallel| Runs the filter after all members are forwarded to performs the analysis. |
| **finalize.template**        | serial | Checks if the whole experiment is finished. If so, stops. Otherwise, resubmits ens_members.${EXPINFO}.lsf for the next assimilation cycle. |

### Diagnostics

A toolbox for diagnostics is provided. Some are written for a specific regional application using
Ferrybox observations. However, it shouldn't be difficult to add new tools following the present
ones. A fortran toolbox post-processes the FESOM outputs and visualization is done using Generic
Mapping Tools (GMT). DART outputs are plotted using FERRET. Please see the description inside respective code.

| Directory   | code file                       |  description  |
| ------------|:--------------------------------|---------------|
| src/        |                                 ||
|             |fesom_post_main.F90              | main fortran routine calling each tool by a switch |
|             |fesom_ocean_mod.F90              | ocean diagnostic routines |
|             |fesom_dart_mod.F90               | dart outputs diagnostic routines |
|             |fesom_forcing_mod.F90            | forcing diagnostic routines |
|             |fesom_observation_mod.F90        | observation diagnostic routines |
|             |gen_input.F90                    | fesom routines for I/O |
|             |gen_modules_clock.F90            | fesom routines for timing |
|             |gen_modules_config.F90           | fesom routines for configuration |
|             |mesh_read.F90                    | fesom routines for reading the mesh |
|             |Makefile                         | Makefile reading mkmf_template |
|             |oce_dens_press.F90               | fesom routines to compute density and pressure|
|             |oce_mesh_setup.F90               | fesom routines for mesh setup |
|             |oce_modules.F90                  | fesom routines for ocean modules |
|             |random_perturbation.F90          | random perturbation to observation sampling |
|             |utilities.F90                    | various utilities |
| script/     |                                 ||
|             |compute_ensemble_mean            | computes ensemble mean and extracts a transect or level |
|             |compute_increment                | computes increment using preassim_mean.nc and postassim_mean.nc |
|             |compute_innovation               | computes innovation if there is a nature run |
|             |dart_obs_seq_diag                | dart statistics from obs_epoch.nc and obs_diag.nc|
|             |dart.postproc.env                | dart environment variables |
|             |fesom.postproc.env               | fesom environment variables|
|             |observe_nature_run               | sampling observations if there is a nature run |
|             |transect_daily_mean              | extracts and plots a transect of a individual ensemble member|
|             |zlevel_daily_mean                | extracts and plots a level of an individual ensemble member|
| gmt/        |                                 ||
|             |plot_ensemble_mean.gmt           | plots ensemble mean called by compute_ensemble_mean|
|             |plot_increment.gmt               | plots increment called by compute_increment|
|             |plot_innovation.gmt              | plots innovation called by compute_innovation|
|             |transect_daily_mean.gmt          | plots transects called by transect_daily_mean|
|             |zlevel_yearly_mean.gmt           | plots levels called by zlevel_daily_mean|
| ferret/     |                                 ||
|             |frt.obs_diag_TeMPLaTe.jnl        | plot dart diags called by dart_obs_seq_diag|
|             |frt.obs_epoch_TeMPLaTe.jnl       | plots dart diags called by dart_obs_seq_diag|

### References

- Androsov, A., Fofonova, V., Kuznetsov, I., Danilov, S., Rakowsky, N., Harig, S., Brix, H., and
Wiltshire, K. H.: FESOM-C v.2: coastal dynamics on hybrid unstructured meshes, Geosci. Model Dev.,
12, 1009-1028, https://doi.org/10.5194/gmd-12-1009-2019, 2019.

- Aydoğdu, A., Hoar, T. J., Vukicevic, T., Anderson, J. L., Pinardi, N., Karspeck, A., Hendricks,
J., Collins, N., Macchia, F., and Özsoy, E.: OSSE for a sustainable marine observing network in the
Sea of Marmara, Nonlin. Processes Geophys., 25, 537-551, https://doi.org/10.5194/npg-25-537-2018,
2018a.

- Aydoğdu, A., Pinardi, N., Özsoy, E., Danabasoglu, G., Gürses, Ö., and Karspeck, A.: Circulation of
the Turkish Straits System under interannual atmospheric forcing, Ocean Sci., 14, 999-1019,
https://doi.org/10.5194/os-14-999-2018, 2018b.

- Danilov, S., Kivman, G., and Schröter, J.: A finite-element ocean model: principles and
evaluation, Ocean Modell., 6, 125–150, 2004.

- Danilov, S., Sidorenko, D., Wang, Q., and Jung, T.: The Finite-volumE Sea ice–Ocean Model
(FESOM2), Geosci. Model Dev., 10, 765-789, https://doi.org/10.5194/gmd-10-765-2017, 2017.

- Gürses, Ö., Aydoğdu, A., Pinardi, N., and Özsoy, E.: A finite element modeling study of the
Turkish Straits System, in: The Sea of Marmara – Marine Biodiversity, Fisheries, Conservations and
Governance, edited by: Özsoy E., Çaǧatay, M. N., Balkis, N., and
Öztürk, B., TUDAV Publication, 169–184, 2016.

- Wang, Q., Danilov, S., and Schröter, J.: Finite element ocean circulation model based on
triangular prismatic elements, with application in studying the effect of topography representation,
J. Geophys. Res.-Oceans (1978–2012), 113, C05015, https://doi.org/10.1029/2007JC004482, 2008.
