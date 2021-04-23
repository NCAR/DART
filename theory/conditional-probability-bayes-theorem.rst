Conditional probability and Bayes' theorem
==========================================

This section introduces two prerequisite concepts for understanding data
assimilation theory: conditional probability and Bayes' theorem. 

Conditional probability
-----------------------

Most real-world events involve uncertainty because the occurence of a specific
outcome isn't guaranteed. You can sense that in situations in which these are
possible outcomes:

- your flight departs on time
- you keep your New Year's resolution
- your car needs repairs in the next 6 months

there is a chance that the opposite outcome might occur. Describing such
situations accurately requires making probabilistic statements.

In mathematical notation, the probability of an event, :math:`A`, is denoted by
:math:`P(A)`. If the event :math:`A` means that your flight departs on time,
you can write:

.. math::

    P(A) = likely

since most flights do actually depart on time. 

Events usually occur in conjuction with other events, so it is useful to assign
conditional probabilities, or the probability that an outcome occurs if another
event also occurs.

If the event :math:`B` is that a blizzard approaches the airport an hour before
your scheduled departure you can write a conditional probability as 
:math:`P(A|B)`, or the probability that :math:`A` occurs, given that :math:`B`
also occurs. In this case, you can assign the probability that your flight
departs on time given that a blizzard approaches the airport an hour before
your scheduled departure as:

.. math::

    P(A|B) = unlikely

since it is unlikely that your flight departs on time in a blizzard. These
examples use informal, subjective probabilities. But the mathematical notation
can also be used to assign formal, quantitative probabilities as well.

Bayes' theorem
--------------

Imagine you are in a house and the carbon monoxide detector has set off its
alarm. Carbon monoxide is colorless and odorless, so you evacuate the house,
but you don't know whether there are actually significant concentrations of
carbon monoxide inside or if your detector is faulty.

In the United States, 100,000 carbon monoxide exposure events occur in houses
annually and the manufacturer of your detector claims that its detectors have a
0.1% error rate. Bayes' theorem allows you to calculate the quantitative
probability of whether or not there is a carbon monoxide exposure event in the
house, given that the carbon monoxide detector has set off its alarm.

Probability theory allows you to keep track of specific conditions and events.
The names of the relevant terms, and what they represent in this example are:

- the *prior*, :math:`P(A)` - the probability of a carbon monoxide
  exposure event in your house
- the *likelihood*, :math:`P(B|A)` - the probability your detector sets off its
  alarm given that there is a carbon monoxide exposure event in your house
- the *normalization*, :math:`P(B)` - the probablity your detector sets off its
  alarm
- the *posterior*, :math:`P(A|B)` - the probability of a carbon monoxide
  exposure event in your house given that your detector sets of its alarm

If this is your first experience with probability theory, you may be
unaccustomed to the terminology and level of nuance that the theory affords.
Take your time to think through each of the probabilities and conditions.
Notice, for example, the difference between :math:`P(B|A)` and :math:`P(A|B)`.

Bayes' theorem allows you to calculate the probability you want to know,
the posterior, :math:`P(A|B)`. The theorem is:

.. math::

    posterior = \frac{likelihood \times prior}{normalization}

or:

.. math::

    P(A|B) = \frac{P(B|A)P(A)}{P(B)}

To compute the right hand side of the equation you'll need to estimate the
prior, the likelihood, and the normalization.

Prior
~~~~~

You can estimate the probability of a carbon monoxide exposure event in your
house, :math:`P(A)`, by dividing the number of carbon monoxide exposure events
that occur annually in houses by the total number of houses in the United
States, which is 140 million houses:

.. math::

   P(A)=100,000 \div 140,000,000=7.1 \times 10^{-4}

Likelihood
~~~~~~~~~~

You can estimate the probability your detector sets off its alarm given that
there is a carbon monoxide exposure event in your house, :math:`P(B|A)`, since
you know the error rate of the detector, 0.1%:

.. math::

    P(B|A) = 1-0.001 = 0.999

Normalization
~~~~~~~~~~~~~

Estimating the probablity your detector sets off its alarm, :math:`P(B)`,
requires estimating two cases: the probability of a false alarm,
:math:`P(B^-)`, and the probability of a true alarm, :math:`P(B^+)`.

The probability of a false alarm is the portion of the population that does not
experience a carbon monoxide exposure event times the error rate of the
detector:

.. math::

    P(B^-) = \frac{(140,000,000-100,000)}{140,000,000} \times 0.001 = 9.9 \times 10^{-4}

The probability of a true alarm is the portion of the population that
experiences a carbon monoxide exposure event times the rate that the detector
will correctly set off its alarm:

.. math::

    P(B^+) = \frac{100,000}{140,000,000} \times (1-0.001) = 7.1 \times 10^{-4}

:math:`P(B)` is the sum of :math:`P(B^-)` and :math:`P(B^+)`:

.. math::

   P(B)= 9.9 \times 10^{-4} + 7.1 \times 10^{-4} = 1.7\times 10^{-3}

Posterior
~~~~~~~~~

You now have all of the necessary probabilities to estimate the probability of
a carbon monoxide exposure event in your house given that your detector sets
off its alarm, :math:`P(A|B)`:

.. math::

   P(A|B) = \frac{P(B|A)P(A)}{P(B)} = \frac{0.999 \times 7.1 \times 10^{-4}}{1.7\times 10^{-3}} = 0.42

Thus, the posterior probability is 0.42.

Bayesian inference
------------------

One of the primary benefits of Bayes' theorem is that it can be applied
multiple times to update a probability when new information is available. This
process is best illustrated by continuing the example.

While standing outside, you call the fire department. A fire engine arrives and
firefighters enter the house with a carbon monoxide meter. This meter is more
accurate than the one installed in the house. It has an error rate of 0.01%.

The meter detects dangerous levels carbon monoxide in the house. You know
intuitively that it is now highly probable that there are dangerous levels of
carbon monoxide in the house. Bayes' theorem provides a rigorous framework to
support your intuition.

You can apply Bayes' theorem again to update your estimate of the probability
of a carbon monoxide exposure event in the house. This updating process is
called Bayesian inference. 

When applying Bayes' theorem a second time, the process is the same but the 
probabilities involved are different.

Prior
~~~~~

In the first part of the example, you estimated the prior by dividing the
number of carbon monoxide exposure events that occur annually in houses by the
total number of houses in the United States. That was the correct approach at
first. But now your prior is the posterior from the first part:

.. math::

   P(A) = 0.42

since that is the probability of a carbon monoxide exposure event in your
house.

Likelihood
~~~~~~~~~~

Since the firefighters' carbon monoxide meter has a lower error rate than the
detector installed in the house, :math:`P(B|A)` is also different:

.. math::

    P(B|A) = 1-0.0001 = 0.9999

Normalization
~~~~~~~~~~~~~

The probablity that the meter detects carbon monoxide is still comprised of two
parts, the probability of a false detection, :math:`P(B^-)`, and the
probability of a true detection, :math:`P(B^+)`. But since the error rate of
the firefighters' meter is lower and your detector has also set off its alarm,
the normalization is different.

The probability of a false detection is the probability that there isn't a 
carbon monoxide exposure event in the house times the error rate of the meter:

.. math::

    P(B^-) = (1-0.42) \times 0.0001 = 5.8 \times 10^{-5}

The probability of a true detection is the probability that there is a carbon
monoxide exposure event in the house times the rate that the meter will
correctly detect it:

.. math::

    P(B^+) = 0.42 \times 0.9999 = 0.42

:math:`P(B)` is the sum of :math:`P(B^-)` and :math:`P(B^+)`:

.. math::
    
    P(B) = 5.8 \times 10^{-5} + 0.42 = 0.42

Posterior
~~~~~~~~~
    
You have all of the necessary probabilities to estimate the probability of
a carbon monoxide exposure event in your house given that both your detector
set off its alarm and the firefighters' meter also detected carbon monoxide,
:math:`P(A|B)`:

.. math::

    P(A|B) = \frac{P(B|A)P(A)}{P(B)} = \frac{0.9999 \times 0.42}{0.42} = 0.9999

Thus, the second posterior probability is 0.9999. This makes sense intuitively:
it is extremely likely that there is a carbon monoxide exposure event in the
house if both your alarm and the firefighters' meter detect carbon monoxide.

It also demonstrates the ability of Bayes' theorem to update the probability of
an event when new information becomes available.

With these concepts you can now begin the :doc:`readme`.
