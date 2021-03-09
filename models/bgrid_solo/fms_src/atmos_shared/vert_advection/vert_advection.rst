module vert_advection_mod
=========================

Overview
--------

Computes the tendency due to vertical advection for an arbitrary quantity.

.. container::

   The advective tendency may be computed in *advective form* (for use in spectral models) or *flux form*. The spatial
   differencing may be *centered* (second or fourth order) or *finite volume* (van Leer) using a piecewise linear
   method.

Other modules used
------------------

.. container::

   ::

           fms_mod

Public interface
----------------

.. container::

   ::


      use vert_advection_mod [, only:  vert_advection,
                                       SECOND_CENTERED, FOURTH_CENTERED, VAN_LEER_LINEAR,
                                       FLUX_FORM, ADVECTIVE_FORM  ]

   `vert_advection <vert_advection.html#vert_advection>`__:
      Computes the vertical advective tendency for an arbitrary quantity. There is no initialization routine necessary.

| 

Public routines
---------------

a. .. rubric:: Vert_advection
      :name: vert_advection

   ::

      call vert_advection ( dt, w, dz, r, rdt [, mask, scheme, form] )

      DESCRIPTION
         This routine computes the vertical advective tendency for
         an arbitrary quantity. The tendency can be computed using
         one of several different choices for spatial differencing
         and numerical form of the equations. These choices are
         controlled through optional arguments. 
         There is no initialization routine necessary.

      INPUT
         dt   time step in seconds [real]

         w    advecting velocity at the vertical boundaries of the grid boxes
              does not assume velocities at top and bottom are zero
              units = [units of dz / second]
              [real, dimension(:,:,:)]
              [real, dimension(:,:)]
              [real, dimension(:)]

         dz   depth of model layers in arbitrary units (usually pressure)
              [real, dimension(:,:,:)]
              [real, dimension(:,:)]
              [real, dimension(:)]

         r    advected quantity in arbitrary units
              [real, dimension(:,:,:)]
              [real, dimension(:,:)]
              [real, dimension(:)]

      OUTPUT
         rdt  advective tendency for quantity "r" weighted by depth of layer
              units = [units of r * units of dz / second]
              [real, dimension(:,:,:)]
              [real, dimension(:,:)]
              [real, dimension(:)]

      OPTIONAL INPUT
         mask    mask for below ground layers,
                 where mask > 0 for layers above ground
                 [real, dimension(:,:,:)]
                 [real, dimension(:,:)]
                 [real, dimension(:)]

         scheme  spatial differencing scheme, use one of these values:
                    SECOND_CENTERED = second-order centered
                    FOURTH_CENTERED = fourth-order centered
                    VAN_LEER_LINEAR = piecewise linear, finite volume (van Leer)
                 [integer, default=VAN_LEER_LINEAR]

         form    form of equations, use one of these values:
                    FLUX_FORM      = solves for -d(w*r)/dt
                    ADVECTIVE_FORM = solves for -w*d(r)/dt
                 [integer, default=FLUX_FORM]

      NOTE
         The input/output arrays may be 1d, 2d, or 3d.
         The last dimension is always the vertical dimension.
         For the vertical dimension the following must be true:
         size(w,3) == size(dz,3)+1 == size(r,3)+1 == size(rdt,3)+1 == size(mask,3)+1
         All horizontal dimensions must have the same size (no check is done).

Error messages
--------------

.. container::

   **Errors in vert_advection_mod**
      vertical dimension of input arrays inconsistent
      The following was not true: size(w,3) = size(r,3)+1.
      invalid value for optional argument scheme
      The value of optional argument scheme must be one of the public parameters SECOND_CENTERED, FOURTH_CENTERED, or
      VAN_LEER_LINEAR.
      invalid value for optional argument form
      The value of optional argument form must be one of the public parameters FLUX_FORM or ADVECTIVE_FORM.

| 

References
----------

.. container::

   #. Lin, S.-J., W.C. Chao, Y.C. Sud, and G.K. Walker, 1994: A class of the van Leer-type transport schemes and its
      application to the moisture in a general circulation model. *Mon. Wea. Rev.*, **122**, 1575-1593.

| 

Notes
-----

.. container::

   None.

| 
