The Lorenz 63 model and its relevance to data assimilation
==========================================================

This section describes a consequential model in the development of humanity's 
understanding of the limits of predicting nature: the three-variable model of
Lorenz (1963). [1]_ This model captures the essence of chaotic systems and will
serve as an example to deepen your understanding of DART and data assimilation.

In 1963, Edward Lorenz developed a simplified three-variable model to
investigate atmospheric convection. By making several simplifications to the
Boussinesq approximation, the Lorenz model was derived for a single thin layer
of fluid uniformly heated from below and cooled from above. The original paper 
has been cited over 20,000 times. The relatively simple, yet nonlinear, system
of ordinary differential equations is:

.. math::

    \frac{dx}{dt} = \sigma(y-x)

    \frac{dy}{dt} = x(r-z)-y
    
    \frac{dz}{dt} = xy-bz

Here, :math:`x` is proportional to the rate of convection, :math:`y` is related
to the horizontal temperature variation, and :math:`z` is the vertical
temperature variation.

There are three constant parameters:

.. math::

    \sigma=10, r=28, b=8/3

- :math:`\sigma` relates to the Prandtl number
- :math:`r` relates to the Rayleigh number
- :math:`b` relates to the physical dimensions of the layer

Note that two of the equations have nonlinear terms: :math:`\frac{dy}{dt}` has 
the :math:`-xz` term and :math:`\frac{dz}{dt}` has the :math:`xy` term.

Lorenz 63 is a consequential model in the history of science because the
numerical investigation of the **chaos** arising from this system of ordinary
differential equations unexpectedly launched a revolution in humanity's
understanding of nature. These investigations lead to numerous mathematical and
scientific breakthroughs.

While the chaotic nature of certain systems such as the three-body problem had
been investigated previously, it was the electronic computer, which could
compute thousands of calculations per second, that allowed these ideas to be
formalized.

In particular, Lorenz's model made it clear for the first time how an
infinitesimally small change in the initial conditions of a system could end up
having a dramatic effect on the subsequent behavior of the system. Lorenz
discussed the strange behavior of this model in a popular science lecture, 
*The Essence of Chaos* [2]_:

    At one point I decided to repeat some of the computations in order to
    examine what was happening in greater detail. I stopped the computer,
    typed in a line of numbers that it had printed out a while earlier, and
    set it running again. I went down the hall for a cup of coffee and
    returned after about an hour, during which the computer had simulated
    about two months of weather. The numbers being printed out were nothing
    like the old ones. I immediately suspected a weak vacuum tube or some
    other computer trouble, which was not uncommon, but before calling for
    service I decided to see just where the mistake had occurred, knowing that
    this could speed up the servicing process. Instead of a sudden break, I
    found that the new values at first repeated the old ones, but soon
    afterward had differed by one and then several units in the last decimal
    place...The numbers I had typed in were not the exact original numbers,
    but were the rounded-off values that appeared in the original printout.
    The initial round-off errors were the culprits; they were steadily
    amplifying until they dominated the solution. In today's terminology,
    there was chaos.

Lorenz discovered that even in a model with just three variables, a very small
change in the initial conditions (in this case, the numbers he typed back into
the computer, which were very slightly different from the original numbers)
could cause the entire large-scale behavior to change. Lorenz's discovery has
many important practical implications:

1. If tiny changes can grow to dominate a system, it is no longer possible to
   find the one set of "perfect" initial conditions and hope to allow the
   system to run forever with perfect forecasts. Instead, forecasting chaotic
   systems must be approached **statistically**.
2. There is a practical limit of predictability inherent in chaotic systems. In
   other words, the nonlinear dynamics of a chaotic model are inherently
   difficult to predict. Multiple evaluations (an **ensemble**) can be run with
   different plausible initial conditions to quantify this error growth.
3. In order to forecast chaotic systems effectively, periodic observations of
   the state are required to effectively guide the forecast and narrow the
   uncertainty. Since in real-world applications observations are almost always
   sparse compared to the number of state variables, merging observations and
   forecasts (i.e. **data assimilation**) is required to effectively forecast
   chaotic systems.

While Lorenz 63 is a simple example of a chaotic system, there are many other
chaotic systems of real practical interest in areas such as weather prediction,
climate, oceanography, hydrology, ecology, biology, and many other disciplines.

In short, while the Lorenz model is a simple set of equations that can easily
be run on even the most basic of computers today, it is representative of the
same problem of predictability that can be found throughout science. DART
supports the investigation of forecasting chaotic systems in *any* field where
periodic observations can be used to constrain the uncertainty using an
ensemble.

References
----------

.. [1] Lorenz, Edward N., 1963: Deterministic Nonperiodic Flow. *Journal of the
       Atmospheric Sciences*, **20**, 130-141,
       `doi:0.1175/1520-0469(1963)020\<0130:DNF\>2.0.CO;2
       <https://doi.org/10.1175/1520-0469(1963)020\<0130:DNF\>2.0.CO;2>`__

.. [2] Lorenz, Edward N. *The Essence of Chaos*. University of Washington Press, 1995.