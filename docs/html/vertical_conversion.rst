Vertical Conversion of Observations
===================================

In Lanai vertical conversion of observations occurs in get_close_obs. The Lanai code in filter_assim is as follows:

If this algorithm was followed in RMA DART, all processors would have to communicate to calculate the location of
observation i. This is a huge amount of contention since all processors are doing exactly the same calculation so need
to access exactly the same state elements. This causes the code to run very slowly, for example 1 minute for 1000
observations, versus 5 seconds for 1000 observations for Lanai.

However, there is no need to calculate the vertical conversion inside the SEQUENTIAL_OBS do loop, since the mean state
vector used is not updated during the loop. (In Lanai it is the array passed to the model_mod by ens_mean_for_model in
filter_main). Also this calculation does not scale, because all processors do the same calculation.

In DART RMA the owner of an observation converts the vertical location of an observation and broacasts it to all other
processors as part of the broadcast in the SEQUENTIAL_OBS do loop.

The DART RMA code calculates the vertical of all observations before the loop. This potentially scales better because
processors only calculate their own observation conversions, but does require model_mod interfaces for vertical
conversion.

The DART RMA code in filter_assim is as follows:

::

   do i =, obs_ens_handle%my_num_vars
      call convert_vertical_location(my_obs_loc(i))
   end do
   SEQUENTIAL_OBS do i = 1, obs_ens_handle%num_vars
      ...
      broadcast increments and vertical location for observation i
      ...
   enddo

Bitwise problem
~~~~~~~~~~~~~~~

Moving the ``convert_vertical_location`` changes the number of ``get/set location`` calls. There is a bitwise creep of
the location when you do this. This is in the conversion from degrees to radians and back again. If you want to do the
exact number of ``get/set location`` you can change the line lanai_bitwise = .false. to lanai_bitwise = .true. in
assim_tools_mod.f90. Note this is not a namelist option because production code should not be run with lanai_bitwise =
.true. For more detail on running bitwise with Lanai see :doc:`./bitwise_considerations`.
