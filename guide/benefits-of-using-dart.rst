The benefits of using DART
==========================

A common pitfall for graduate students and professionals alike is to look at the
simplicity of data assimilation, in particular ensemble data assimilation, and
decide they can easily write their own DA system. Indeed, this is true. After
learning of the core algorithms, a talented programmer using their favorite
language could write a functional DA system in a manner of weeks if not days.
However, he or she will soon find that while the core of DA systems are easy to
write, the more “real” the system needs to be, the more complex it will become.
Writing a parallel DA system that can efficiently utilize multiple cores with
MPI is not straight-forward, and adding covariance localization, observation
operators, multiple models, and auxiliary tools such as quality control and
pre-processing will quickly dwarf the amount of core DA code, not to mention the
headaches involved in supporting multiple computing environments, compilers,
etc.

DART employs a modular programming approach to apply an algorithm to move the
underlying models toward a state that is more consistent with information from a
set of observations. Models may be swapped in and out, as can different DA
algorithms. The method requires running multiple instances of a model to
generate an ensemble of states. A forward operator appropriate for the type of
observation being assimilated is applied to each of the states to generate the
model’s estimate of the observation.

DART remains the top choice for scientists, educators, and mathematicians
seeking mature and robust ensemble DA solutions without reinventing the wheel.
Here are some of the many benefits of using DART:

1. DART is **freely available, open source, and released under the**
   `Apache 2.0 License <https://www.apache.org/licenses/LICENSE-2.0>`__ **.** In short
   this means that you are granted a copyright license stating you are free to
   use, modify, and redistribute any derivative works derived from the DART
   system provided that you maintain the license and copyright information. Of
   course, we also ask that you credit DART in your publications, and kindly ask
   that you contribute your modifications so that other users may benefit. See
   `How should I cite DART? <#citeDart>`__ and `How can I contribute to
   DART? <#ContributeToDart>`__ for more information.
2. DART is **fully parallel and carefully engineered** to run on systems ranging
   from single-core research computers to the top performing multicore
   supercomputers in the world. Writing scalable parallel code is arguably the
   most difficult and time-consuming task in scientific computing today, but
   DART has already carefully implemented and tested this project, and the code
   is available for you to use out-of-the-box. For more information on how DART
   was written (and continues to be developed), see `DART’s design
   philosophy <#dartDesign>`__.
3. DART contains **numerous tools that accelerate getting started** on both
   research and “real-world” problems. Multiple rigorously tested inflation,
   localization, perturbation, and other auxiliary data assimilation algorithms
   are available for immediate use and testing. See `Important capabilities of
   DART <#dartCapabilities>`__ for more information.
4. DART **makes adding a new model straightforward**. A new model only needs to
   implement a list of (at most) 18 core functions or use the default behavior
   if applicable to take advantage of DART’s mature and robust DA algorithms. A
   basic data assimilation system for a large model can be built in
   person-weeks, and comprehensive systems have been built in a few months. See
   `How do I run DART with my model? <#RunWithMyModel>`__ for more information.
5. DART **makes it easy to add new observations** in order to test their
   potential beneficial impact. Incorporating new observation types only
   requires creating a forward operator that computes the expected value of an
   observation given a model’s state. See `How do I add my observations to
   DART? <#RunWithMyObs>`__ for more information.
6. DART **can be used to test new DA algorithms**. Many such algorithms have
   been successfully implemented, tested, and published using DART. This is not
   covered in this getting started guide as this is an “advanced user”
   functionality, so for this purpose it is best to first get in touch with the
   DART team at dart @ ucar.edu to make the process as smooth as possible.
7. Finally, and perhaps most importantly, DART **has world-class support**
   available from the DART team at NSF NCAR. A talented team of dedicated software
   engineers and data assimilation scientists work together to continually
   improve DART and support user needs. See the `About page <https://dart.ucar.edu/about/>`__ for
   more information about the DART team.
