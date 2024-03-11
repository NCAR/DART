SEIR
====

Overview
--------

The extended SEIR Model with Vaccination was first proposed by Ghostine et al. (2021) [1]_
to simulate the novel coronavirus disease (COVID-19) spread. The model considers 7
stages of infection:

  1. Susceptible (S),
  2. Exposed (E),
  3. Infected (I),
  4. Quarantined (Q),
  5. Recovered (R),
  6. Deaths (D),
  7. Vaccinated (V).

There are several parameters that can be changed to study different cases and regions:

  - :math:`\theta`: New births and new residents per unit of time,
  - :math:`\beta`: Transmission rate divided by the population size,
  - :math:`\alpha`: Vaccination rate,
  - :math:`\mu`: Natural death rate,
  - :math:`\gamma`: Average latent time,
  - :math:`\delta`: Average quarantine time,
  - :math:`\kappa`: Mortality rate, 
  - :math:`\lambda`: Average days until recovery, 
  - :math:`\rho`: Average days until death,
  - :math:`\sigma`: Vaccine in-efficacy (:math:`0 \leq \sigma \leq 1`).

Earth system models are often descritized in space. The state in these models represents
variables at different spatial locations. The variables of the SEIR model describe the 
stage/phase of the disease and they do not have a physical location. To this end, 
techniques such as spatial localization are not applicable in this model. DART assumes 
that all 7 variables belong to the same *virtual* point in space. Any assimilated 
observation will impact all 7 variables.

The SEIR model uses identity observations. Typical observations that can be assimilated
are:
 
  *Recovered*, *Death* and *Vaccinated*. 

Some agencies provide data for "*Confirmed*" cases. This can be used to compute and 
assimilate the number of active (which is equivelant to quarantined) cases as shown: 

*Active/Quarantined (Q) = Confirmed - Recovered (R) - Deaths (D)*

Initial versions of the model were tested using DART_LAB. This was conducted by  
**Shaniah Reece** as part of her SIParCS internship at NSF NCAR (2022).

Namelist
--------

The ``&model_nml`` namelist is read from the ``input.nml`` file. Namelists
start with an ampersand ``&`` and terminate with a slash ``/``. Character
strings that contain a ``/`` must be enclosed in quotes to prevent them from
prematurely terminating the namelist.

.. code-block:: fortran

  &model_nml
     model_size        = 40,
     delta_t           = 0.04167,
     time_step_days    = 0,
     time_step_seconds = 3600,
     num_pop           = 331996199,
     pert_size         = 0.5, 
     t_incub           = 5.6,
     t_infec           = 3.8,
     t_recov           = 14.0,
     t_death           = 7.0,
     alpha             = 0.000001,
     theta             = 12467,
     mu                = 0.000025,
     sigma             = 0.05,
     beta              = 0.00000000136,
     kappa             = 0.00308,
  /

Description of each namelist entry
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

+-------------------+----------+-------------------------------------------+
| Item              | Type     | Description                               |      
+===================+==========+===========================================+
| model_size        | integer  | Number of variables in model.             |   
+-------------------+----------+-------------------------------------------+
| delta_t           | real(r8) | Non-dimensional timestep. This is         |
|                   |          | mapped to the dimensional timestep        |
|                   |          | specified by time_step_days and           |
|                   |          | time_step_seconds.                        |
+-------------------+----------+-------------------------------------------+
| time_step_days    | integer  | Number of days for dimensional            |
|                   |          | timestep, mapped to delta_t.              |
+-------------------+----------+-------------------------------------------+
| time_step_seconds | integer  | Number of seconds for dimensional         |
|                   |          | timestep, mapped to delta_t.              |
+-------------------+----------+-------------------------------------------+
| num_pop           | integer  | Population size.                          |   
+-------------------+----------+-------------------------------------------+
| pert_size         | real(r8) | Size of perturbation used to create       |
|                   |          | an ensemble using a lognormal pdf.        |  
+-------------------+----------+-------------------------------------------+
| t_incub           | real(r8) | Incubation period                         |
|                   |          | :math:`\equiv 1/\gamma`.                  |  
+-------------------+----------+-------------------------------------------+
| t_infec           | real(r8) | Infection time                            |   
|                   |          | :math:`\equiv 1/\delta`.                  | 
+-------------------+----------+-------------------------------------------+  
| t_recov           | real(r8) | Recovery period                           |   
|                   |          | :math:`\equiv 1/\lambda`.                 | 
+-------------------+----------+-------------------------------------------+
| t_death           | real(r8) | Time until death                          |   
|                   |          | :math:`\equiv 1/\rho`.                    | 
+-------------------+----------+-------------------------------------------+  
| alpha             | real(r8) | Vaccination rate. If study period         |
|                   |          | starts before vaccination is              | 
|                   |          | available, this must be set to 0.         | 
+-------------------+----------+-------------------------------------------+  
| theta             | integer  | New birth and new residents.              |   
+-------------------+----------+-------------------------------------------+  
| mu                | real(r8) | Natural death rate.                       |   
+-------------------+----------+-------------------------------------------+ 
| sigma             | real(r8) | Vaccination inefficacy (e.g., if the      |
|                   |          | vaccine is 95% effective, then            |
|                   |          | :math:`\sigma = 1-0.95 = 0.05`).          |   
+-------------------+----------+-------------------------------------------+
| beta              | real(r8) | Transmission rate divided by population   |
|                   |          | size.                                     |
+-------------------+----------+-------------------------------------------+ 
| kappa             | real(r8) | Mortality rate.                           |   
+-------------------+----------+-------------------------------------------+ 
  
References
----------

.. [1] Ghostine, R.; Gharamti, M.; Hassrouny, S.; Hoteit, I. An Extended SEIR Model with Vaccination for Forecasting the COVID-19 Pandemic in Saudi Arabia Using an Ensemble Kalman Filter. Mathematics 2021, 9, 636. https://dx.doi.org/10.3390/math9060636.
